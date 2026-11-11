"""Monthly climatology of upper-ocean temperature error: products C1 and C2 (docs/SERVER_PLAN.md).

One task per (year, month) downloads a sample of days, accumulates the reanalysis-ensemble spread, and matches
that month's Argo against every ensemble member and the ensemble mean; the build step pools the tasks by calendar
month into:

  C1  reanalysis spread: sigma(month, depth, lat, lon) = sqrt(mean over sampled days of the 3-member variance),
      the 0-50 m layer mean, and per-member 0-50 m anomaly (member - mean) maps -- the perturbation realisations;
  C2  Argo minus each member and minus the ensemble mean, by calendar month (bias, RMS, std, counts, RMS maps, and
      the Hollingsworth-Lonnberg split as a secondary diagnostic), with the C1 spread sampled at the profiles.

The validation question is whether member - mean is a realistic perturbation: for a reliable 3-member ensemble,
RMS(Argo - member) ~ sqrt(2) * spread and RMS(Argo - mean) ~ sqrt(4/3) * spread.

Written for a plain multi-core server (no Slurm): ``run`` fans the tasks out across a bounded process pool. Every
task and product skips outputs already on disk, so a run resumes. Run ``python -m ocean_uncertainty.climatology --help``.
"""

import argparse
import calendar
import concurrent.futures
import datetime
import pathlib

import gsw
import numpy as np
import pandas as pd
import xarray as xr

from ocean_uncertainty import sources
from ocean_uncertainty.fetch import fetch_argo, fetch_copernicus, select_region
from ocean_uncertainty.observations import day_index, load_argo, sample, to_levels
from ocean_uncertainty.products import innovation_statistics, reanalysis_members, write
from ocean_uncertainty.regrid import to_grid

MEMBERS = list(sources.GREP["members"])          # ["glor", "cglo", "oras"]
BACKGROUNDS = ["mean"] + MEMBERS                  # what Argo is differenced against (+ "rg" once RG is added)


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

    The GREP ensemble is downloaded one sampled day at a time (geo-series). Accumulated over the days: the 3-member
    variance (N - 1), the day count, and each member's squared anomaly from the mean. The sampled days' Argo is
    differenced against every member and the ensemble mean, and kept with the C1 spread at each profile. Existing
    outputs are skipped, so the task resumes.
    """
    data_dir = pathlib.Path(data_dir)
    tag = f"{year}{month:02d}"
    spread_path, innov_path = data_dir / f"spread_{tag}.nc", data_dir / f"innovations_{tag}.nc"

    member_names = [f"{variable}_{m}" for m in MEMBERS]
    day_paths = []
    for d in sampled_days(year, month, step_days):
        day_paths += fetch_copernicus(sources.GREP["members_dataset"], member_names, region, f"{d:%Y-%m-%d}",
                                      f"{d:%Y-%m-%d}", data_dir / "grep", max_depth, daily=True)

    last = calendar.monthrange(year, month)[1]
    fetch_argo(region, f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last:02d}", max_pressure, data_dir / "argo")

    if spread_path.exists() and innov_path.exists():   # resume: nothing left to compute for this month
        return spread_path, innov_path

    sum_variance = sum_anom2 = None
    n_days = 0
    mean_fields, sigma_fields = [], []
    member_fields = {m: [] for m in MEMBERS}
    for path in day_paths:
        members = reanalysis_members(path, variable)             # (member, time=1, depth, lat, lon)
        mean = members.mean("member")
        variance = members.var("member", ddof=1).astype(np.float64)
        anom2 = ((members - mean) ** 2).astype(np.float64).sum("time")   # (member, depth, lat, lon)
        day_var = variance.sum("time")
        sum_variance = day_var if sum_variance is None else sum_variance + day_var
        sum_anom2 = anom2 if sum_anom2 is None else sum_anom2 + anom2
        n_days += members.sizes["time"]
        mean_fields.append(mean)
        sigma_fields.append(np.sqrt(variance))
        for m in MEMBERS:
            member_fields[m].append(members.sel(member=m))
    _save(xr.Dataset({"sum_variance": sum_variance, "sum_anom2": sum_anom2, "n_days": n_days}), spread_path)

    with xr.open_dataset(day_paths[0]) as first:
        levels, lon0 = first["depth"].values, first["longitude"].values[0]
    argo_files = sorted((data_dir / "argo").glob(f"argo_*_{tag}??_*.nc"))
    table = to_levels(load_argo(argo_files), levels)
    table["lon"] = (table["lon"] - lon0) % 360 + lon0   # background grid's own longitude range (dateline-safe)
    table["sigma"] = _sample_daily(table, sigma_fields)
    table["innov_mean"] = table["obs"] - _sample_daily(table, mean_fields)
    for m in MEMBERS:
        table[f"innov_{m}"] = table["obs"] - _sample_daily(table, member_fields[m])
    matched = table[np.isfinite(table["innov_mean"])]   # sampled days only
    _save(xr.Dataset.from_dataframe(matched.reset_index(drop=True)), innov_path)
    return spread_path, innov_path


def _calendar_month(path: pathlib.Path) -> int:
    """Calendar month (1-12) from a spread/innovations/argo filename ending in ...YYYYMM[DD]_... ."""
    for token in pathlib.Path(path).stem.split("_"):
        if token.isdigit() and len(token) >= 6:
            return int(token[4:6])
    raise ValueError(f"no YYYYMM in {path}")


def _pool_by_month(spread_files, name):
    """Sum a per-(year, month) field and the day count by calendar month."""
    field, days = {}, {}
    for path in spread_files:
        m = _calendar_month(path)
        with xr.open_dataset(path) as ds:
            if name not in ds:
                return None, None
            v, n = ds[name].load(), int(ds["n_days"])
        field[m] = v if m not in field else field[m] + v
        days[m] = n if m not in days else days[m] + n
    return field, days


def build_c1(spread_files, out: pathlib.Path, max_depth: float, variable: str = "thetao") -> pathlib.Path:
    """Pool the partial sums by calendar month into the spread sigma and per-member 0-50 m anomaly maps."""
    var_sum, days = _pool_by_month(spread_files, "sum_variance")
    months = sorted(var_sum)
    variance = xr.concat([var_sum[m] / days[m] for m in months], pd.Index(months, name="month"))
    layer = to_grid(variance, depth_edges=[0.0, max_depth]).squeeze("depth", drop=True)
    data = {f"{variable}_std": np.sqrt(variance),
            f"{variable}_std_0_{int(max_depth)}m": np.sqrt(layer),
            "n_days": ("month", [days[m] for m in months])}

    anom_sum, _ = _pool_by_month(spread_files, "sum_anom2")
    if anom_sum is not None:   # per-member (member - mean) 0-50 m RMS: the perturbation realisations
        anom = xr.concat([anom_sum[m] / days[m] for m in months], pd.Index(months, name="month"))
        anom_layer = to_grid(anom, depth_edges=[0.0, max_depth]).squeeze("depth", drop=True)
        data[f"{variable}_member_anomaly_0_{int(max_depth)}m"] = np.sqrt(anom_layer)

    return write(xr.Dataset(data), out, "C1 reanalysis-spread climatology", [sources.GREP["citation"]],
                 "Per calendar month: the 3-member variance (N - 1) summed over the sampled days / day count, "
                 "square-rooted; the 0-%d m mean variance square-rooted; and per member the 0-%d m RMS of "
                 "(member - mean), the error realisation used for perturbation." % (max_depth, max_depth),
                 "Three members sharing forcing, model family and observations: the spread probably underestimates "
                 "the error. Each sampled day adds 2 degrees of freedom (SERVER_PLAN.md).")


# ---------------------------------------------------------------------------------------------------------------
# Roemmich-Gilson Argo climatology: a model-free reference. Argo - RG isolates representativeness error.

RG_MEAN, RG_ANOM = "ARGO_TEMPERATURE_MEAN", "ARGO_TEMPERATURE_ANOMALY"


def rg_temperature(rg_dir: pathlib.Path, year: int, month: int, levels, region) -> xr.DataArray:
    """The Roemmich-Gilson mapped temperature (mean + that month's anomaly), cut to ``region`` and interpolated to
    ``levels`` (reanalysis depths). RG is 1 deg, monthly, Argo-only -- no ocean model. Months after 2018 come from
    the monthly extension files; the long-term mean always comes from the base file.
    """
    rg_dir = pathlib.Path(rg_dir)
    index = (year - 2004) * 12 + (month - 1)
    with xr.open_dataset(rg_dir / "RG_ArgoClim_Temperature_2019.nc", decode_times=False) as base:
        mean = base[RG_MEAN].load()
        if index < base.sizes["TIME"]:
            anomaly = base[RG_ANOM].isel(TIME=index).load()
        else:
            with xr.open_dataset(rg_dir / f"RG_ArgoClim_{year}{month:02d}_2019.nc", decode_times=False) as ext:
                anomaly = ext[RG_ANOM].isel(TIME=0).load()
    field = (mean + anomaly).rename(LATITUDE="latitude", LONGITUDE="longitude")
    # Pressure (dbar) to depth (m) at the region's mid-latitude; the two differ by under 1% in the upper ocean.
    depth = -gsw.z_from_p(field["PRESSURE"].values, 0.5 * (region[2] + region[3]))
    field = field.assign_coords(PRESSURE=depth).rename(PRESSURE="depth").sortby("depth")
    field = select_region(field, region, lon="longitude", lat="latitude")
    return field.interp(depth=levels).expand_dims(time=[0])


def add_rg_innovations(innovation_files, rg_dir: pathlib.Path, region) -> None:
    """Add an ``innov_rg`` column (Argo minus Roemmich-Gilson) to each innovations file, in place."""
    for path in innovation_files:
        tag = next(t for t in pathlib.Path(path).stem.split("_") if t.isdigit() and len(t) >= 6)
        year, month = int(tag[:4]), int(tag[4:6])
        with xr.open_dataset(path) as ds:
            table = ds.to_dataframe().reset_index(drop=True)
        field = rg_temperature(rg_dir, year, month, np.sort(table["depth"].unique()), region)
        rg = sample(field, table["lon"].values, table["lat"].values, table["depth"].values,
                    np.zeros(len(table), int))
        table["innov_rg"] = table["obs"].values - rg
        _save(xr.Dataset.from_dataframe(table), path)


def build_c2(innovation_files, out: pathlib.Path, **statistics) -> pathlib.Path:
    """Argo minus each member and minus the mean, by calendar month, with the C1 spread sampled at the profiles.

    For every background in {mean, members} and every calendar month, the innovation statistics (bias, RMS, std,
    counts, RMS maps, and the Hollingsworth-Lonnberg split) are computed and stacked on a ``background`` dimension.
    """
    frames = []
    for path in innovation_files:
        with xr.open_dataset(path) as ds:
            frames.append(ds.to_dataframe().assign(month=_calendar_month(path)))
    table = pd.concat(frames, ignore_index=True)
    # Profile numbers were per file; renumber globally (as observations.to_levels does) so pooling years of one
    # calendar month does not merge different profiles in the vertical-correlation pivot.
    table["profile"] = table.groupby(["platform", "cycle", "direction"], sort=False).ngroup()

    backgrounds = [b for b in BACKGROUNDS + ["rg"] if f"innov_{b}" in table.columns]
    per_background = []
    for background in backgrounds:
        per_month = []
        for month, group in table.groupby("month"):
            stats = innovation_statistics(group.assign(innovation=group[f"innov_{background}"]), **statistics)
            per_month.append(stats.expand_dims(month=[month]))
        per_background.append(xr.concat(per_month, "month", join="outer").expand_dims(background=[background]))
    ds = xr.concat(per_background, "background", join="outer")

    # C1 spread sampled at the profiles, per calendar month and depth (independent of the background).
    sigma = []
    for month, group in table.groupby("month"):
        t = group[np.isfinite(group["sigma"])]
        rms = np.sqrt(t.assign(s2=t["sigma"] ** 2).groupby("depth")["s2"].mean())
        sigma.append(xr.DataArray(rms.values, dims="depth", coords={"depth": rms.index.values},
                                  name="c1_sigma").expand_dims(month=[month]))
    ds["c1_sigma"] = xr.concat(sigma, "month").reindex(depth=ds["depth"])

    citations = [sources.ARGO["citation"], sources.GREP["citation"]]
    if "rg" in backgrounds:
        citations.append(sources.RG["citation"])
    return write(ds, out, "C2 Argo minus reanalysis members and mean, by month", citations,
                 "Argo (research mode) interpolated in depth to the reanalysis levels; innovation against each "
                 "member and the ensemble mean on the profile's sampled day. Per calendar month and background: "
                 "bias, RMS, std, counts, 2 deg RMS maps, and the Hollingsworth-Lonnberg split (a secondary "
                 "diagnostic). c1_sigma is the C1 spread at the same profiles.",
                 "The reanalysis assimilates Argo, so these differences underestimate its error. member - mean is a "
                 "realistic perturbation if RMS(Argo - member) ~ sqrt(2) * spread (see module docstring).")


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

    b = sub.add_parser("build", help="pool the tasks into C1 and C2 (after they finish)")
    b.add_argument("--data-dir", type=pathlib.Path, required=True)
    b.add_argument("--out-dir", type=pathlib.Path, required=True)
    b.add_argument("--max-depth", type=float, default=50.0)
    b.add_argument("--region", type=float, nargs=4, default=[140, 240, 20, 60])
    b.add_argument("--rg-dir", type=pathlib.Path, help="Roemmich-Gilson files; if given, add the Argo - RG background")
    _add_stats_args(b)

    a = p.parse_args()
    if a.command == "month":
        month_task(a.year, a.month, a.region, a.max_depth, a.max_pressure, a.step_days, a.data_dir)
    elif a.command == "run":
        run(a.years, a.months, a.jobs, region=a.region, max_depth=a.max_depth, max_pressure=a.max_pressure,
            step_days=a.step_days, data_dir=a.data_dir)
    else:
        stats = {k: getattr(a, k) for k in STATS}
        innovations = sorted(a.data_dir.glob("innovations_*.nc"))
        if a.rg_dir:
            add_rg_innovations(innovations, a.rg_dir, a.region)
        build_c1(sorted(a.data_dir.glob("spread_*.nc")), a.out_dir / "c1.nc", a.max_depth)
        build_c2(innovations, a.out_dir / "c2.nc", **stats)


if __name__ == "__main__":
    main()
