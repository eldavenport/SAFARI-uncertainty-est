"""Checks of the products against each other at the Argo profiles (docs/PRODUCT_CATALOGUE.md, "Checks").

1. Consistency (Desroziers et al. 2005): the mean squared innovation against the reanalysis-ensemble mean should
   equal the ensemble variance plus the observation-error variance, at each level.
2. Ordering: reanalysis spread (U1) <~ Argo minus reanalysis (U2) << Argo minus climatology (U3), with the SST
   analysis uncertainty (U4) alongside at every level.
2b. Patterns: the sigma maps of U1-U5 on U2's map boxes, compared by size and by rank correlation.
3. Hold-out: each profile's innovation predicted from the other profiles of its window by the U6 analysis, for
   each choice of error model. Tests both the increments and the uncertainty the analysis claims for them.

Results are reported, never used to adjust a product. Run ``python -m ocean_uncertainty.checks --help``.
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
import xarray as xr

from scipy.linalg import cho_factor, cho_solve

from ocean_uncertainty.analysis import BACKGROUND_MODELS, OBSERVATION_MODELS, Window, load_window
from ocean_uncertainty.observations import day_index, sample
from ocean_uncertainty.stats import distance_km


def obs_table(path: pathlib.Path) -> pd.DataFrame:
    """The matched observation table of a U2 or U3 file, with the "obs_" prefix removed."""
    with xr.open_dataset(path) as ds:
        names = [v for v in ds.data_vars if ds[v].dims == ("obs",)]
        table = ds[names].to_dataframe()
    table.columns = [c.removeprefix("obs_") for c in table.columns]
    return table[np.isfinite(table["innovation"])]


def at_obs(field_files, variable: str, table: pd.DataFrame) -> np.ndarray:
    """A daily field (one or more files) sampled at each observation on its day."""
    out = np.full(len(table), np.nan)
    for path in field_files:
        with xr.open_dataset(path) as ds:
            field = ds[variable].load()
        if "time" not in field.dims:
            field = field.expand_dims("time")
            index = np.zeros(len(table), int)
        else:
            index = day_index(field["time"].values, table["time"])
        rows = np.flatnonzero(index >= 0)
        if "depth" not in field.dims:
            field = field.expand_dims(depth=[0.0])
            depth = np.zeros(rows.size)
        else:
            depth = table["depth"].values[rows]
        out[rows] = sample(field, table["lon"].values[rows], table["lat"].values[rows], depth, index[rows])
    return out


def consistency(u2_file, u1_files, u5_file, variable: str) -> pd.DataFrame:
    """Per level: mean innovation squared (after removing the bias) against U1 variance + U5 observation variance."""
    table = obs_table(u2_file)
    table["u1_var"] = at_obs(u1_files, f"{variable}_std", table) ** 2
    with xr.open_dataset(u5_file) as u5:
        instrument = u5.attrs["instrument_variance"]
    table["u5_var"] = at_obs([u5_file], "representativeness_variance", table) + instrument
    table = table.dropna(subset=["u1_var", "u5_var"])
    rows = []
    for depth, t in table.groupby("depth"):
        d = t["innovation"] - t["innovation"].mean()
        expected = t["u1_var"].mean() + t["u5_var"].mean()
        rows.append(dict(depth=depth, n=len(t), innovation_var=np.mean(d**2), u1_var=t["u1_var"].mean(),
                         u5_var=t["u5_var"].mean(), ratio=np.mean(d**2) / expected))
    return pd.DataFrame(rows)


def ordering(u2_file, u3_file, u1_files, u4_file, variable: str) -> pd.DataFrame:
    """Per U2 level: U1 spread, U2 and U3 innovation standard deviations, U4, and whether U1 <= U2 < U3 holds.

    U3 sits on the climatology's levels; each U2 level is compared with the nearest U3 level. U4 is a surface value
    (OSTIA estimates the temperature below the layer of daily surface warming), repeated at each level so the
    depth where it stops matching U1 and U2 can be seen.
    """
    u2, u3 = obs_table(u2_file), obs_table(u3_file)
    u2["u1_var"] = at_obs(u1_files, f"{variable}_std", u2) ** 2
    u2["u4_var"] = at_obs([u4_file], "sst_std", u2) ** 2
    u3_std = u3.groupby("depth")["innovation"].std()
    rows = []
    for depth, t in u2.groupby("depth"):
        nearest = u3_std.index[np.argmin(np.abs(u3_std.index - depth))]
        s1, s2, s3 = np.sqrt(t["u1_var"].mean()), t["innovation"].std(), u3_std[nearest]
        rows.append(dict(depth=depth, n=len(t), u1_std=s1, u2_std=s2, u3_depth=nearest, u3_std=s3,
                         u4_std=np.sqrt(t["u4_var"].mean()),
                         holds=bool(s1 <= s2 < s3) if np.isfinite([s1, s2, s3]).all() else None))
    return pd.DataFrame(rows)


def box_mean(field: xr.DataArray, box_deg: float) -> pd.Series:
    """Mean of a gridded field over the boxes of U2's maps (centres at floor(x / box) * box + box / 2)."""
    centre = lambda x: np.floor(x / box_deg) * box_deg + box_deg / 2
    boxed = field.assign_coords(map_longitude=centre(field["longitude"]), map_latitude=centre(field["latitude"]))
    return boxed.groupby(["map_latitude", "map_longitude"]).mean(...).to_series()


def patterns(u1_files, u2_file, u3_file, u4_file, u5_file, depth: float, min_count: int, variable: str):
    """Sigma of each product at the level nearest ``depth``, on U2's map boxes holding >= ``min_count`` Argo values
    in both U2 and U3: box means, median ratio to U2, and Spearman rank correlations between the maps."""
    with xr.open_dataset(u2_file) as u2, xr.open_dataset(u3_file) as u3:
        box = u2.attrs["map_bin_deg"]
        k2 = int(np.argmin(np.abs(u2["depth"].values - depth)))
        k3 = int(np.argmin(np.abs(u3["depth"].values - depth)))
        level = float(u2["depth"][k2])
        maps = {"U2": (u2["rms_map"].isel(depth=k2) ** 2).where(u2["count_map"].isel(depth=k2) >= min_count),
                "U3": (u3["rms_map"].isel(depth=k3) ** 2).where(u3["count_map"].isel(depth=k3) >= min_count)}
        maps = {k: v.load().to_series() for k, v in maps.items()}
    total, days = 0.0, 0
    for path in u1_files:
        with xr.open_dataset(path) as ds:
            total = total + (ds[f"{variable}_std"].sel(depth=level, method="nearest") ** 2).sum("time").load()
            days += ds.sizes["time"]
    maps["U1"] = box_mean(total / days, box)
    with xr.open_dataset(u4_file) as u4:
        maps["U4"] = box_mean((u4["sst_std"] ** 2).mean("time").load(), box)
    with xr.open_dataset(u5_file) as u5:
        maps["U5"] = box_mean(u5["representativeness_variance"].sel(depth=level, method="nearest").load(), box)
    sigma = np.sqrt(pd.DataFrame(maps)[["U1", "U2", "U3", "U4", "U5"]].dropna())
    summary = pd.DataFrame({"box_mean_sigma_K": sigma.mean(), "median_ratio_to_U2": sigma.div(sigma["U2"], axis=0).median()})
    return level, len(sigma), summary, sigma.corr(method="spearman")


def hold_out(window: Window, n_folds: int) -> pd.DataFrame:
    """Predict each held-out profile's innovations from the profiles of the other folds (profile number mod folds).

    Returns one row per held-out observation: the innovation d, the prediction H dx from the other profiles, the
    variance the error model expects for d - H dx (H A H^T + R on the diagonal), the distance (km) to the nearest
    profile used, and the background standard deviation at the observation. Comparing (d - H dx)^2 with the expected variance tests whether the analysis is over- or
    under-confident.
    """
    obs, points, errors = window.obs, window.obs_points, window.errors
    profile = obs["profile"].values
    fold = profile % n_folds
    rows = []
    for f in range(n_folds):
        held, used = fold == f, fold != f
        if not held.any() or not used.any():
            continue
        S = errors.background(points[used], points[used]) + errors.observation(
            profile[used], points.level[used], window.obs_variance[used])
        factor = cho_factor(S, lower=True)
        cross = errors.background(points[held], points[used])
        gain = cho_solve(factor, cross.T).T
        d = obs["innovation"].values
        predicted = gain @ d[used]
        expected = (points[held].sigma ** 2 - np.sum(gain * cross, axis=1)
                    + np.diag(errors.observation(profile[held], points.level[held], window.obs_variance[held])))
        nearest = np.min(distance_km(points.lon[held], points.lat[held], points.lon[used], points.lat[used]),
                         axis=1)
        rows.append(pd.DataFrame(dict(profile=profile[held], depth=obs["depth"].values[held], innovation=d[held],
                                      predicted=predicted, expected_var=expected, nearest_km=nearest,
                                      background_std=points[held].sigma)))
    return pd.concat(rows, ignore_index=True)


def hold_out_summary(table: pd.DataFrame, distance_edges_km, n_sigma_classes: int) -> pd.DataFrame:
    """Hold-out scores by distance to the nearest profile used and by background standard deviation.

    - skill: 1 - sum (d - H dx)^2 / sum d^2 (0: no better than no update; dominated by the largest errors);
    - median_ratio: median |d - H dx| / median |d| (below 1: the typical innovation is predicted);
    - calibration: mean (d - H dx)^2 / mean expected variance (1: the stated uncertainty is right).
    Background standard deviation classes are quantiles of the table, so each holds the same number of values.
    """
    sigma_edges = table["background_std"].quantile(np.linspace(0, 1, n_sigma_classes + 1)).values
    table = table.assign(distance=pd.cut(table["nearest_km"], distance_edges_km),
                         sigma=pd.cut(table["background_std"], sigma_edges, include_lowest=True))
    groups = [(("all", "all"), table)] + list(table.groupby(["distance", "sigma"], observed=True))
    out = []
    for (distance, sigma), t in groups:
        residual = t["innovation"] - t["predicted"]
        out.append(dict(nearest_km=str(distance), background_std_K=str(sigma), n=len(t),
                        skill=1 - (residual**2).sum() / (t["innovation"] ** 2).sum(),
                        median_ratio=residual.abs().median() / t["innovation"].abs().median(),
                        calibration=(residual**2).mean() / t["expected_var"].mean()))
    return pd.DataFrame(out)


def u1_file_for(date: pd.Timestamp, u1_files) -> pathlib.Path:
    """The U1 file whose days include ``date``."""
    for path in u1_files:
        with xr.open_dataset(path) as ds:
            if date.floor("D") in pd.DatetimeIndex(ds["time"].values):
                return path
    raise ValueError(f"no U1 file covers {date.date()}")


def compare_error_models(u2_file, u1_files, u5_file, dates, window_days: float, n_folds: int, variable: str,
                         distance_edges_km, n_sigma_classes: int) -> pd.DataFrame:
    """Hold-out skill and calibration for every background model x observation model, pooled over ``dates``."""
    rows = []
    for background_model in BACKGROUND_MODELS:
        for obs_model in OBSERVATION_MODELS:
            tables = []
            for date in dates:
                window = load_window(u2_file, u1_file_for(date, u1_files), u5_file, str(date.date()), window_days,
                                     variable, background_model, obs_model)
                tables.append(hold_out(window, n_folds).assign(date=date))
            summary = hold_out_summary(pd.concat(tables, ignore_index=True), distance_edges_km, n_sigma_classes)
            rows.append(summary.assign(model=f"{background_model} / {obs_model}"))
    table = pd.concat(rows, ignore_index=True)
    return table.sort_values(["nearest_km", "background_std_K", "model"], kind="stable")[
        ["nearest_km", "background_std_K", "model", "n", "skill", "median_ratio", "calibration"]]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="check", required=True)
    rep = sub.add_parser("report", help="checks 1 and 2")
    rep.add_argument("--u3", type=pathlib.Path, required=True)
    rep.add_argument("--u4", type=pathlib.Path, required=True)
    rep.add_argument("--pattern-depths", type=float, nargs="+", required=True, help="depths (m) of the map comparison")
    rep.add_argument("--min-count", type=int, required=True, help="fewest Argo values in a map box to compare it")
    hold = sub.add_parser("holdout", help="check 3, for every error model")
    hold.add_argument("--start", required=True, help="centre of the first window")
    hold.add_argument("--end", required=True, help="no window is centred after this date")
    hold.add_argument("--window-days", type=float, required=True,
                      help="window length; windows are centred this many days apart, so none overlap")
    hold.add_argument("--folds", type=int, required=True, help="number of hold-out groups of profiles")
    hold.add_argument("--distance-edges-km", type=float, nargs="+", required=True,
                      help="edges of the distance-to-nearest-profile classes in the summary")
    hold.add_argument("--sigma-classes", type=int, required=True,
                      help="number of equally populated background standard deviation classes")
    for parser in (rep, hold):
        parser.add_argument("--u1", type=pathlib.Path, nargs="+", required=True)
        parser.add_argument("--u2", type=pathlib.Path, required=True, help="U2 against the GREP ensemble mean")
        parser.add_argument("--u5", type=pathlib.Path, required=True, help="U5 built on the GREP grid")
        parser.add_argument("--out", type=pathlib.Path, required=True, help="Markdown report")
        parser.add_argument("--variable", default="thetao")
    a = p.parse_args()

    if a.check == "report":
        sections = {
            "Consistency: innovation variance / (U1 variance + U5 observation variance); 1 if consistent, "
            "> 1 if U1 is under-dispersive": consistency(a.u2, a.u1, a.u5, a.variable),
            "Ordering: U1 <= U2 < U3 (U4 for comparison)": ordering(a.u2, a.u3, a.u1, a.u4, a.variable),
        }
        for depth in a.pattern_depths:
            level, n, summary, corr = patterns(a.u1, a.u2, a.u3, a.u4, a.u5, depth, a.min_count, a.variable)
            sections[f"Patterns at {level:.1f} m ({n} boxes): size"] = summary.reset_index(names="product")
            sections[f"Patterns at {level:.1f} m: rank correlation"] = corr.reset_index(names="product")
    else:
        dates = pd.date_range(a.start, a.end, freq=pd.Timedelta(days=a.window_days))
        sections = {
            f"Hold-out over {len(dates)} windows of {a.window_days:g} days ({a.start} to {dates[-1].date()}), "
            f"{a.folds} folds (scores: checks.hold_out_summary)":
            compare_error_models(a.u2, a.u1, a.u5, dates, a.window_days, a.folds, a.variable, a.distance_edges_km,
                                 a.sigma_classes),
        }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as f:
        for title, df in sections.items():
            f.write(f"## {title}\n\n{df.to_markdown(index=False, floatfmt='.3g')}\n\n")
    print(a.out.read_text())


if __name__ == "__main__":
    main()
