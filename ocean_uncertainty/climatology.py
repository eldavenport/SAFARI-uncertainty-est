"""Monthly climatology of upper-ocean temperature error: products C1-C3 (docs/SERVER_PLAN.md).

One task per (year, month) downloads a sample of days, accumulates the reanalysis-ensemble variance and matches
that month's Argo against the ensemble mean; the build step pools the tasks by calendar month into:

  C1  sigma(month, depth, lat, lon) = sqrt(mean over sampled days of the 3-member variance), and the 0-50 m layer mean;
  C2  Argo minus the reanalysis ensemble mean, by calendar month (bias, RMS, the C1 sigma at the profiles, the
      Hollingsworth-Lonnberg split and correlation length, and RMS maps);
  C3  Argo minus the WOA23 climatology, by calendar month (the same statistics; an upper bound).

Written for a plain multi-core server (no Slurm): ``run`` fans the tasks out across a bounded process pool. Every
task and product skips outputs already on disk, so a run resumes. Run ``python -m ocean_uncertainty.climatology --help``.
"""

import argparse
import calendar
import concurrent.futures
import datetime
import pathlib

import numpy as np
import pandas as pd
import xarray as xr

from ocean_uncertainty import sources
from ocean_uncertainty.fetch import fetch_argo, fetch_copernicus, fetch_woa23
from ocean_uncertainty.observations import day_index, load_argo, sample, to_levels
from ocean_uncertainty.products import (innovation_statistics, innovation_table, reanalysis_members, write)
from ocean_uncertainty.regrid import to_grid


def sampled_days(year: int, month: int, step_days: int) -> list[datetime.date]:
    """Days of the month to sample: 1, 1 + step, ... up to day 28.

    Capped at 28 so every calendar month, February included, is sampled on the same days and adds the same degrees
    of freedom (SERVER_PLAN.md, "Open choices"): a 5-day step gives 1, 6, 11, 16, 21, 26; a 10-day step 1, 11, 21.
    """
    return [datetime.date(year, month, d) for d in range(1, 29, step_days)]


def _save(ds: xr.Dataset, path: pathlib.Path) -> pathlib.Path:
    """Write an intermediate file in full (double) precision, atomically. Products use products.write instead."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(".part")
    ds.to_netcdf(partial)
    partial.rename(path)
    return path


def _sample_daily(table: pd.DataFrame, fields: list[xr.DataArray]) -> np.ndarray:
    """Sample a list of daily fields at the table rows whose day each field covers (NaN elsewhere)."""
    out = np.full(len(table), np.nan)
    for field in fields:
        index = day_index(field["time"].values, table["time"])
        rows = np.flatnonzero(index >= 0)
        out[rows] = sample(field, table["lon"].values[rows], table["lat"].values[rows],
                           table["depth"].values[rows], index[rows])
    return out


def month_task(year: int, month: int, region, max_depth: float, max_pressure: float, step_days: int,
               data_dir: pathlib.Path, variable: str = "thetao") -> tuple[pathlib.Path, pathlib.Path]:
    """Build ``spread_YYYYMM.nc`` and ``innovations_YYYYMM.nc`` for one (year, month); keep that month's Argo file.

    The GREP ensemble is downloaded one sampled day at a time (geo-series). The 3-member variance (N - 1, double
    precision) is summed over the days with the day count, and the sampled days' Argo innovations against the
    ensemble mean are kept with that day's sigma. Existing outputs are skipped, so the task resumes.
    """
    data_dir = pathlib.Path(data_dir)
    tag = f"{year}{month:02d}"
    spread_path, innov_path = data_dir / f"spread_{tag}.nc", data_dir / f"innovations_{tag}.nc"

    member_names = [f"{variable}_{m}" for m in sources.GREP["members"]]
    day_paths = []
    for d in sampled_days(year, month, step_days):
        day_paths += fetch_copernicus(sources.GREP["members_dataset"], member_names, region, f"{d:%Y-%m-%d}",
                                      f"{d:%Y-%m-%d}", data_dir / "grep", max_depth, daily=True)

    last = calendar.monthrange(year, month)[1]
    fetch_argo(region, f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last:02d}", max_pressure, data_dir / "argo")

    if spread_path.exists() and innov_path.exists():   # resume: nothing left to compute for this month
        return spread_path, innov_path

    sum_variance, n_days, mean_fields, sigma_fields = None, 0, [], []
    for path in day_paths:
        members = reanalysis_members(path, variable)
        variance = members.var("member", ddof=1).astype(np.float64)   # (time=1, depth, lat, lon)
        day = variance.sum("time")
        sum_variance = day if sum_variance is None else sum_variance + day
        n_days += members.sizes["time"]
        mean_fields.append(members.mean("member"))
        sigma_fields.append(np.sqrt(variance))
    _save(xr.Dataset({"sum_variance": sum_variance, "n_days": n_days}), spread_path)

    with xr.open_dataset(day_paths[0]) as first:
        levels, lon0 = first["depth"].values, first["longitude"].values[0]
    argo_files = sorted((data_dir / "argo").glob(f"argo_*_{tag}??_*.nc"))
    table = to_levels(load_argo(argo_files), levels)
    table["lon"] = (table["lon"] - lon0) % 360 + lon0   # background grid's own longitude range (dateline-safe)
    table["background"] = _sample_daily(table, mean_fields)
    table["sigma"] = _sample_daily(table, sigma_fields)
    table["innovation"] = table["obs"] - table["background"]
    matched = table[np.isfinite(table["innovation"])]   # sampled days only
    _save(xr.Dataset.from_dataframe(matched.reset_index(drop=True)), innov_path)
    return spread_path, innov_path


def _calendar_month(path: pathlib.Path) -> int:
    """Calendar month (1-12) from a spread/innovations/argo filename ending in ...YYYYMM[DD]_... ."""
    for token in pathlib.Path(path).stem.split("_"):
        if token.isdigit() and len(token) >= 6:
            return int(token[4:6])
    raise ValueError(f"no YYYYMM in {path}")


def build_c1(spread_files, out: pathlib.Path, max_depth: float, variable: str = "thetao") -> pathlib.Path:
    """Pool the partial sums by calendar month into sigma = sqrt(Sigma variance / Sigma days), with the layer mean."""
    pooled: dict[int, list] = {}
    for path in spread_files:
        with xr.open_dataset(path) as ds:
            s, n = ds["sum_variance"].load(), int(ds["n_days"])
        acc = pooled.setdefault(_calendar_month(path), [None, 0])
        acc[0] = s if acc[0] is None else acc[0] + s
        acc[1] += n
    months = sorted(pooled)
    variance = xr.concat([pooled[m][0] / pooled[m][1] for m in months], pd.Index(months, name="month"))
    layer = to_grid(variance, depth_edges=[0.0, max_depth]).squeeze("depth", drop=True)
    ds = xr.Dataset({f"{variable}_std": np.sqrt(variance),
                     f"{variable}_std_0_{int(max_depth)}m": np.sqrt(layer),
                     "n_days": ("month", [pooled[m][1] for m in months])})
    return write(ds, out, "C1 reanalysis-spread climatology", [sources.GREP["citation"]],
                 "Per calendar month, the 3-member variance (N - 1) summed over the sampled days and divided by the "
                 "day count, square-rooted; the layer value is the 0-%d m mean variance, square-rooted." % max_depth,
                 "Three members sharing forcing, model family and observations: the spread probably underestimates "
                 "the error. Each sampled day adds 2 degrees of freedom (SERVER_PLAN.md).")


def _by_month_statistics(groups, out: pathlib.Path, product: str, citations, processing, limitations,
                         extra=None, **statistics) -> pathlib.Path:
    """Run innovation_statistics on each (month, table) in ``groups`` and concatenate along a ``month`` dimension."""
    per_month = []
    for month, table in groups:
        stats = innovation_statistics(table, **statistics)
        if extra is not None:
            stats = extra(stats, table)
        per_month.append(stats.expand_dims(month=[month]))
    # join="outer": each month's RMS maps cover different boxes, so take the union of map coordinates.
    return write(xr.concat(per_month, "month", join="outer"), out, product, citations, processing, limitations)


def build_c2(innovation_files, out: pathlib.Path, **statistics) -> pathlib.Path:
    """Argo minus the reanalysis ensemble mean by calendar month, with the C1 sigma at the profiles and the ratio."""
    frames = []
    for path in innovation_files:
        with xr.open_dataset(path) as ds:
            frames.append(ds.to_dataframe().assign(month=_calendar_month(path)))
    table = pd.concat(frames, ignore_index=True)
    # Profile numbers were per file; renumber globally (as observations.to_levels does) so pooling years of one
    # calendar month does not merge different profiles in the vertical-correlation pivot.
    table["profile"] = table.groupby(["platform", "cycle", "direction"], sort=False).ngroup()

    def add_sigma(stats, t):
        t = t[np.isfinite(t["innovation"]) & np.isfinite(t["sigma"])]
        rms = np.sqrt(t.assign(s2=t["sigma"] ** 2).groupby("depth")["s2"].mean())
        c1 = xr.DataArray(rms.values, dims="depth", coords={"depth": rms.index.values}).reindex(depth=stats["depth"])
        return stats.assign(c1_sigma=c1, c1_over_rms=c1 / stats["rms"])

    return _by_month_statistics(table.groupby("month"), out, "C2 Argo minus reanalysis by month",
                                [sources.ARGO["citation"], sources.GREP["citation"]],
                                "Argo (research mode) interpolated in depth to the reanalysis levels; innovation "
                                "against the ensemble mean on the profile's sampled day. Per calendar month: bias, "
                                "RMS, the C1 sigma at the profiles and their ratio, the Hollingsworth-Lonnberg "
                                "split and correlation length, and RMS maps.",
                                "The reanalysis assimilates Argo, so these differences underestimate its error. The "
                                "fitted length is the horizontal scale for a random perturbation field.",
                                extra=add_sigma, **statistics)


def build_c3(argo_files, woa_files, out: pathlib.Path, variable: str = "thetao", **statistics) -> pathlib.Path:
    """Argo minus the WOA23 climatology by calendar month (an upper bound), the same statistics as C2."""
    by_month: dict[int, list] = {}
    for path in argo_files:
        by_month.setdefault(_calendar_month(path), []).append(path)

    def tables():
        for month in sorted(by_month):
            table, _ = innovation_table(by_month[month], "woa23", list(woa_files), variable)
            yield month, table

    return _by_month_statistics(tables(), out, "C3 Argo minus climatology by month",
                                [sources.ARGO["citation"], sources.WOA23["citation"]],
                                "Argo interpolated in depth to the WOA23 levels; background interpolated linearly "
                                "between mid-months. Same per-month statistics as C2.",
                                "A climatology carries no day information: an upper bound on the analysis error.",
                                **statistics)


def _run_one(kw):
    month_task(**kw)
    return kw["year"], kw["month"]


def run(years, months, jobs, **task_kwargs) -> None:
    """Fan the (year, month) tasks out across ``jobs`` worker processes (bounded, to be a polite shared-box neighbour)."""
    tasks = [dict(year=y, month=m, **task_kwargs) for y in years for m in months]
    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as pool:
        for future in concurrent.futures.as_completed(pool.submit(_run_one, kw) for kw in tasks):
            y, m = future.result()   # re-raise any worker error
            print(f"done {y}-{m:02d}", flush=True)


STATS = dict(max_distance_km=1000.0, bin_km=25.0, max_days_apart=3.0, min_pairs=100, map_bin_deg=2.0)


def _add_stats_args(parser) -> None:
    for name, default in STATS.items():
        parser.add_argument("--" + name.replace("_", "-"), type=type(default), default=default)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    def task_args(parser):
        parser.add_argument("--region", type=float, nargs=4, default=[140, 240, 20, 60],
                            metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"))
        parser.add_argument("--max-depth", type=float, default=50.0)
        parser.add_argument("--max-pressure", type=float, default=100.0)
        parser.add_argument("--step-days", type=int, default=5)
        parser.add_argument("--data-dir", type=pathlib.Path, required=True)

    one = sub.add_parser("month", help="one (year, month) task")
    one.add_argument("--year", type=int, required=True)
    one.add_argument("--month", type=int, required=True)
    task_args(one)

    many = sub.add_parser("run", help="all tasks across a process pool")
    many.add_argument("--years", type=int, nargs="+", default=list(range(2011, 2023)))
    many.add_argument("--months", type=int, nargs="+", default=list(range(1, 13)))
    many.add_argument("--jobs", type=int, default=12, help="worker processes (also roughly concurrent downloads)")
    task_args(many)

    b = sub.add_parser("build", help="pool the tasks into C1-C3 (after they finish)")
    b.add_argument("--data-dir", type=pathlib.Path, required=True)
    b.add_argument("--out-dir", type=pathlib.Path, required=True)
    b.add_argument("--region", type=float, nargs=4, default=[140, 240, 20, 60])
    b.add_argument("--max-depth", type=float, default=50.0)
    b.add_argument("--woa-period", default="decav91C0")
    b.add_argument("--woa-resolution", default="0.25")
    _add_stats_args(b)

    a = p.parse_args()
    if a.command == "month":
        month_task(a.year, a.month, a.region, a.max_depth, a.max_pressure, a.step_days, a.data_dir)
    elif a.command == "run":
        run(a.years, a.months, a.jobs, region=a.region, max_depth=a.max_depth, max_pressure=a.max_pressure,
            step_days=a.step_days, data_dir=a.data_dir)
    else:
        stats = {k: getattr(a, k) for k in STATS}
        spread = sorted(a.data_dir.glob("spread_*.nc"))
        innov = sorted(a.data_dir.glob("innovations_*.nc"))
        argo = sorted((a.data_dir / "argo").glob("argo_*.nc"))
        woa = fetch_woa23(a.woa_period, a.woa_resolution, list(range(1, 13)), a.region, a.max_depth,
                          a.data_dir / "woa23")
        build_c1(spread, a.out_dir / "c1.nc", a.max_depth)
        build_c2(innov, a.out_dir / "c2.nc", **stats)
        build_c3(argo, woa, a.out_dir / "c3.nc", **stats)


if __name__ == "__main__":
    main()
