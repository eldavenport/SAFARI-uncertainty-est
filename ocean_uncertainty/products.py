"""Build the uncertainty products U1-U6 described in docs/PRODUCT_CATALOGUE.md from downloaded sources.

Run ``python -m ocean_uncertainty.products <product> --help``. Each product is a NetCDF file whose attributes give
its sources, processing and known limitations.
"""

import argparse
import datetime
import pathlib

import numpy as np
import pandas as pd
import xarray as xr

from ocean_uncertainty import sources
from ocean_uncertainty.observations import day_index, load_argo, month_weights, sample, to_levels
from ocean_uncertainty.regrid import cell_bounds, level_interfaces, subgrid_covariance, to_grid
from ocean_uncertainty.analysis import BACKGROUND_MODELS, OBSERVATION_MODELS, load_window
from ocean_uncertainty.stats import (analysis_update, ensemble_spread, fit_background_and_observation_error,
                                     pair_covariance)


def write(ds: xr.Dataset, path: pathlib.Path, product: str, citations: list[str], processing: str,
          limitations: str) -> pathlib.Path:
    """Save a product with its provenance in the global attributes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.attrs.update(product=product, sources=" | ".join(citations), processing=processing,
                    limitations=limitations, catalogue="docs/PRODUCT_CATALOGUE.md",
                    created=datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"))
    # Single precision is ample for storage (computations that need double precision convert first); compressed.
    encoding = {name: {"dtype": "float32", "zlib": True, "complevel": 1}
                for name, var in ds.data_vars.items() if np.issubdtype(var.dtype, np.floating)}
    partial = path.with_suffix(".part")
    ds.to_netcdf(partial, encoding=encoding)
    partial.rename(path)
    return path


# ---------------------------------------------------------------------------------------------------------------
# U1: spread of the reanalysis ensemble

def reanalysis_members(path: pathlib.Path, variable: str) -> xr.DataArray:
    """The GREP members of ``variable`` in one downloaded file, stacked along a ``member`` dimension."""
    with xr.open_dataset(path) as ds:
        names = list(sources.GREP["members"])
        return xr.concat([ds[f"{variable}_{m}"].load() for m in names], dim=pd.Index(names, name="member"))


def build_u1(grep_files: list[pathlib.Path], variable: str, out_dir: pathlib.Path) -> list[pathlib.Path]:
    """Ensemble mean, member anomalies and standard deviation (N - 1) of the GREP reanalyses, one file per input.

    Output files that already exist are kept, so an interrupted build resumes.
    """
    written = []
    for path in grep_files:
        target = out_dir / f"u1_{path.name}"
        if target.exists():   # resume: months already built are kept
            written.append(target)
            continue
        members = reanalysis_members(path, variable)
        anomalies, std = ensemble_spread(members)
        ds = xr.Dataset({f"{variable}_mean": members.mean("member"), f"{variable}_anomaly": anomalies,
                         f"{variable}_std": std})
        ds[f"{variable}_std"].attrs.update(long_name="standard deviation across reanalyses (divided by N - 1)",
                                           units=members.attrs.get("units", ""))
        ds["member"].attrs["description"] = "; ".join(f"{k}: {v}" for k, v in sources.GREP["members"].items())
        written.append(write(
            ds, target,
            "U1 reanalysis-ensemble spread", [sources.GREP["citation"]],
            "Members stacked; anomalies from the ensemble mean; standard deviation with N - 1 in the denominator.",
            "Three members sharing forcing, model family and observations: the spread probably underestimates the "
            "error, and a single day's standard deviation has only 2 degrees of freedom."))
    return written


# ---------------------------------------------------------------------------------------------------------------
# U2 and U3: Argo minus a reanalysis, and Argo minus a climatology

def innovations_daily(table: pd.DataFrame, fields) -> pd.DataFrame:
    """Add the background value from daily fields (an iterable of (time, depth, latitude, longitude) arrays)."""
    background = np.full(len(table), np.nan)
    for field in fields:
        index = day_index(field["time"].values, table["time"])
        rows = np.flatnonzero(index >= 0)
        background[rows] = sample(field, table["lon"].values[rows], table["lat"].values[rows],
                                  table["depth"].values[rows], index[rows])
    return table.assign(background=background, innovation=table["obs"] - background)


def innovations_monthly(table: pd.DataFrame, climatology: xr.DataArray) -> pd.DataFrame:
    """Add the background from a monthly climatology (time = month 1-12), interpolated between mid-months."""
    first, second, w = month_weights(table["time"])
    months = climatology["time"].values

    def index(month):
        # -1 marks a month that was not downloaded; sample() returns NaN for it.
        return np.where(np.isin(month, months), np.searchsorted(months, month), -1)

    args = (table["lon"].values, table["lat"].values, table["depth"].values)
    background = (1 - w) * sample(climatology, *args, index(first)) + w * sample(climatology, *args, index(second))
    return table.assign(background=background, innovation=table["obs"] - background)


def grep_mean_fields(grep_files, variable):
    for path in grep_files:
        yield reanalysis_members(path, variable).mean("member")


def glorys12_fields(glorys_files, variable):
    for path in glorys_files:
        with xr.open_dataset(path) as ds:
            yield ds[variable].load()


def woa23_climatology(woa_files) -> xr.DataArray:
    """WOA23 monthly objectively analysed means with time = month number and CF dimension names."""
    fields = []
    for path in sorted(woa_files):
        with xr.open_dataset(path, decode_times=False) as ds:
            month = int(path.stem.split("_t")[-1].split("_")[0])
            fields.append(ds["t_an"].load().assign_coords(time=[month]))
    return xr.concat(fields, "time").rename(lat="latitude", lon="longitude")


def innovation_table(argo_files, background: str, files, variable: str) -> tuple[pd.DataFrame, str]:
    """Argo values on the background's levels with the background and the innovation (observation minus background).

    ``background`` is "grep_mean" (U2), "glorys12" (U2) or "woa23" (U3).
    """
    points = load_argo(argo_files)
    if background == "woa23":
        climatology = woa23_climatology(files)
        levels, lon_axis = climatology["depth"].values, climatology["longitude"].values
    else:
        with xr.open_dataset(files[0]) as first:
            levels, lon_axis = first["depth"].values, first["longitude"].values
    table = to_levels(points, levels)
    # Longitudes in the background grid's own range, so maps of a region across the dateline stay in one piece.
    table["lon"] = (table["lon"] - lon_axis[0]) % 360 + lon_axis[0]
    if background == "woa23":
        return innovations_monthly(table, climatology), sources.WOA23["citation"]
    if background == "grep_mean":
        return innovations_daily(table, grep_mean_fields(files, variable)), sources.GREP["citation"]
    return innovations_daily(table, glorys12_fields(files, variable)), sources.GLORYS12["citation"]


def innovation_statistics(table: pd.DataFrame, *, max_distance_km: float, bin_km: float, max_days_apart: float,
                          min_pairs: int, map_bin_deg: float) -> xr.Dataset:
    """Per-level bias, spread and the background/observation split; maps of RMS; correlation between levels."""
    table = table[np.isfinite(table["innovation"])]
    levels = np.sort(table["depth"].unique())
    per_level, curves = [], []
    for depth in levels:
        t = table[table["depth"] == depth]
        binned = pair_covariance(t["lon"].values, t["lat"].values, t["time"].values, t["innovation"].values,
                                 max_distance_km=max_distance_km, bin_km=bin_km,
                                 max_time_apart=np.timedelta64(int(max_days_apart * 86400), "s"))
        try:
            fit = fit_background_and_observation_error(binned, min_pairs)
        except (ValueError, RuntimeError):
            fit = dict.fromkeys(["background_variance", "length_km", "observation_variance",
                                 "background_variance_se", "length_km_se"], np.nan)
        per_level.append(dict(bias=binned["bias"], rms=np.sqrt(np.mean(t["innovation"] ** 2)),
                              std=np.sqrt(binned["variance"]), count=binned["points"], **fit))
        curves.append((binned["covariance"], binned["pairs"]))

    ds = xr.Dataset({k: ("depth", [row[k] for row in per_level]) for k in per_level[0]},
                    coords={"depth": levels, "distance_km": binned["distance_km"]})
    ds["pair_covariance"] = (("depth", "distance_km"), np.array([c for c, _ in curves]))
    ds["pair_count"] = (("depth", "distance_km"), np.array([n for _, n in curves]))

    lon_bin = np.floor(table["lon"] / map_bin_deg) * map_bin_deg + map_bin_deg / 2
    lat_bin = np.floor(table["lat"] / map_bin_deg) * map_bin_deg + map_bin_deg / 2
    grouped = table.assign(lon_bin=lon_bin, lat_bin=lat_bin, sq=table["innovation"] ** 2).groupby(
        ["depth", "lat_bin", "lon_bin"])
    maps = xr.Dataset.from_dataframe(pd.DataFrame({"rms_map": np.sqrt(grouped["sq"].mean()),
                                                   "count_map": grouped["sq"].size()}))
    ds = xr.merge([ds, maps.rename(lat_bin="map_latitude", lon_bin="map_longitude")])

    wide = table.pivot_table(index="profile", columns="depth", values="innovation")
    ds["vertical_correlation"] = (("depth", "depth_2"), wide.corr(min_periods=min_pairs).reindex(
        index=levels, columns=levels).values)
    ds = ds.assign_coords(depth_2=levels)
    ds.attrs.update(max_distance_km=max_distance_km, bin_km=bin_km, max_days_apart=max_days_apart,
                    min_pairs=min_pairs, map_bin_deg=map_bin_deg)
    return ds


def build_innovations(argo_files, background: str, files, variable: str, out: pathlib.Path,
                      **statistics) -> pathlib.Path:
    """U2 (background grep_mean or glorys12) or U3 (background woa23): the matched table and its statistics."""
    table, citation = innovation_table(argo_files, background, files, variable)
    stats = innovation_statistics(table, **statistics)
    # Table columns get an "obs_" prefix in the file so they cannot clash with the statistics' coordinates.
    obs = xr.Dataset.from_dataframe(table.reset_index(drop=True).rename_axis("obs").add_prefix("obs_"))
    ds = xr.merge([stats, obs], compat="override")
    ds.attrs["background"] = background
    product = "U3 Argo minus climatology" if background == "woa23" else "U2 Argo minus reanalysis"
    limitation = ("A climatology carries no information about the day: this is an upper bound for analysis error."
                  if background == "woa23" else
                  "The reanalysis assimilates Argo, so these differences underestimate its error (a lower bound).")
    return write(ds, out, product, [sources.ARGO["citation"], citation],
                 "Argo (argopy research mode) interpolated linearly in depth to the background levels (pressure "
                 "to depth with TEOS-10); background bilinear in the horizontal on the profile's day (daily) or "
                 "linearly between mid-months (monthly). Pair covariances by distance at each level; Gaussian fit "
                 "to distances > 0 splits background and observation error (Hollingsworth & Lonnberg 1986).",
                 limitation + " The Gaussian covariance shape is an assumption; the binned curves are included.")


# ---------------------------------------------------------------------------------------------------------------
# U4: SST analysis uncertainty, and its documented extension below the surface

def build_u4(ostia_files, out: pathlib.Path) -> pathlib.Path:
    """OSTIA's analysis error standard deviation for the downloaded period, unchanged apart from renaming."""
    with xr.open_mfdataset(ostia_files) as ds:
        sigma = ds[sources.OSTIA["variable"]].load().rename("sst_std")
    return write(sigma.to_dataset(), out, "U4 SST analysis uncertainty", [sources.OSTIA["citation"]],
                 "OSTIA analysis_error copied as sst_std.",
                 "Self-assessed by the analysis; calibration against independent data is mixed (catalogue U4).")


def build_u4_profile(ostia_files, grep_files, out: pathlib.Path) -> pathlib.Path:
    """SST uncertainty extended down to the mixed-layer depth on the reanalysis grid (an assumption, kept separate).

    OSTIA variance is averaged onto the GREP grid (mean point-error variance in each cell), and held constant from
    the surface to the ensemble-mean mixed-layer depth, zero below.
    """
    pieces = []
    with xr.open_mfdataset(ostia_files) as sst:
        variance = sst[sources.OSTIA["variable"]].load() ** 2
    for path in grep_files:
        mld = reanalysis_members(path, "mlotst").mean("member")
        with xr.open_dataset(path) as ds:
            depth = ds["depth"].values
            lon, lat = ds["longitude"].values, ds["latitude"].values
        on_grid = to_grid(variance.sel(time=mld["time"]), cell_bounds(lon), cell_bounds(lat))
        inside = xr.DataArray(depth, dims="depth", coords={"depth": depth}) <= mld
        pieces.append(np.sqrt(on_grid.assign_coords(longitude=lon, latitude=lat)).where(inside, 0.0))
    sigma = xr.concat(pieces, "time").rename("temperature_std").transpose("time", "depth", "latitude", "longitude")
    return write(sigma.to_dataset(), out, "U4 SST uncertainty extended over the mixed layer",
                 [sources.OSTIA["citation"], sources.GREP["citation"]],
                 "OSTIA variance averaged onto the GREP grid; constant to the GREP ensemble-mean mixed-layer depth.",
                 "The vertical extension is an assumption to be tested against U2's vertical innovation profile.")


# ---------------------------------------------------------------------------------------------------------------
# U5: observation-error variance of an Argo temperature for a given target grid

def build_u5(glorys_files, grid_file: pathlib.Path, variable: str, out: pathlib.Path) -> pathlib.Path:
    """Representativeness error of an Argo value for the grid of ``grid_file`` (its longitude, latitude, depth).

    GLORYS12 (1/12 degree) is first averaged onto the target layers. Within each target cell, the variance of
    what remains at each level is the representativeness variance, and the covariance between levels gives its
    correlation between levels (pooled over the region and the days). Within a profile this error is strongly
    correlated in depth (an unresolved eddy shifts the whole mixed layer), so treating a profile's levels as
    independent would count one piece of information many times.
    """
    with xr.open_dataset(grid_file) as grid:
        target = {d: grid[d].values for d in ("depth", "latitude", "longitude")}
    lon_edges, lat_edges = cell_bounds(target["longitude"]), cell_bounds(target["latitude"])
    depth_edges = level_interfaces(target["depth"])
    n = target["depth"].size
    variance_maps, pooled = [], np.zeros((n, n))
    for path in glorys_files:
        with xr.open_dataset(path) as ds:
            # Double precision so the pooled covariance between levels stays positive semi-definite.
            layers = to_grid(ds[variable].load().astype(np.float64), depth_edges=depth_edges).assign_coords(
                depth=target["depth"])
        layers = layers.squeeze("time", drop=True) if "time" in layers.dims else layers
        variance_k = []
        for k in range(n):
            cov = subgrid_covariance(layers.isel(depth=k), layers, lon_edges, lat_edges)
            pooled[k] += cov.sum(("latitude", "longitude")).values
            variance_k.append(cov.isel(depth=k).clip(min=0.0))
        variance_maps.append(xr.concat(variance_k, "depth"))
    variance = xr.concat(variance_maps, "time").mean("time").assign_coords(target)
    scale = np.sqrt(np.diag(pooled))
    corr = pooled / np.outer(scale, scale)
    instrument = sources.ARGO["temperature_accuracy_K"] ** 2
    ds = xr.Dataset({"representativeness_variance": variance.transpose("depth", "latitude", "longitude"),
                     "representativeness_level_correlation": (("depth", "depth_2"), corr)},
                    coords={"depth_2": target["depth"]})
    ds.attrs["instrument_variance"] = instrument
    return write(ds, out, "U5 Argo temperature observation-error covariance",
                 [sources.ARGO["citation"], sources.GLORYS12["citation"]],
                 f"GLORYS12 averaged onto the target layers; within-cell variance and covariance between levels on "
                 f"the target grid, averaged over {len(variance_maps)} days. Instrument accuracy "
                 f"{sources.ARGO['temperature_accuracy_K']} K (Wong et al. 2020) squared, in instrument_variance. "
                 "R between two levels of one profile: sqrt(var_k var_l) * correlation(k, l) (+ instrument "
                 "variance on the diagonal); different profiles uncorrelated.",
                 "Variability finer than 1/12 degree, faster than daily, or within a target layer is missing: a "
                 "lower bound. The correlation between levels is one value for the whole region.")


# ---------------------------------------------------------------------------------------------------------------
# U6: increments that real Argo profiles would make to the reanalysis-ensemble mean

def build_u6(u2_file: pathlib.Path, u1_file: pathlib.Path, u5_file: pathlib.Path, date: str, window_days: float,
             out: pathlib.Path, variable: str, background_model: str, obs_model: str, length_km: float | None = None,
             profile: int | None = None, block_rows: int | None = None) -> pathlib.Path:
    """Optimal-interpolation increment to the GREP ensemble mean from the Argo profiles near ``date``.

    The error models are described in :mod:`ocean_uncertainty.analysis` and in the catalogue (U6).
    ``block_rows`` limits memory by processing that many state points at a time; results do not depend on it.
    """
    w = load_window(u2_file, u1_file, u5_file, date, window_days, variable, background_model, obs_model, length_km)
    obs, points, variance = w.obs, w.obs_points, w.obs_variance
    if profile is not None:
        keep = obs["profile"].values == profile
        if not keep.any():
            raise ValueError(f"profile {profile} is not in the window")
        obs, points, variance = obs[keep], points[keep], variance[keep]
    level = points.level
    R = w.errors.observation(obs["profile"].values, level, variance)
    n = w.state_points.lon.size
    rows = block_rows or n
    blocks = (w.errors.background(w.state_points[slice(i, i + rows)], points) for i in range(0, n, rows))
    increment, reduction = analysis_update(obs["innovation"].values, R, w.errors.background(points, points), blocks)
    background_var = w.state_points.sigma**2

    def to_field(flat):
        full = np.full(w.wet.size, np.nan)
        full[w.wet] = flat
        return xr.DataArray(full.reshape(w.state.shape[:3]), dims=("depth", "latitude", "longitude"),
                            coords={d: w.state[d] for d in ("depth", "latitude", "longitude")})

    ds = xr.Dataset({"increment": to_field(increment), "background_std": to_field(np.sqrt(background_var)),
                     "analysis_std": to_field(np.sqrt(np.maximum(background_var - reduction, 0.0)))})
    ds["obs_profile"] = ("obs", obs["profile"].values)
    ds["obs_innovation"] = ("obs", obs["innovation"].values)
    ds["level_correlation"] = (("depth", "depth_2"), w.errors.level_corr)
    ds = ds.assign_coords(depth_2=w.state["depth"].values)
    ds.attrs.update(date=str(pd.Timestamp(date).date()), window_days=window_days, background_model=background_model,
                    obs_model=obs_model, horizontal_length_km=w.errors.length_km, n_obs=len(obs),
                    n_profiles=int(obs["profile"].nunique()), single_profile=-1 if profile is None else profile)
    return write(ds, out, "U6 Argo-driven analysis increment", [sources.GREP["citation"], sources.ARGO["citation"]],
                 f"Optimal interpolation with the {background_model} background covariance and the {obs_model} "
                 "observation-error variance (correlated between the levels of a profile, U5); increments to the "
                 "GREP ensemble mean.",
                 "Ensemble: three members (rank 2 before localisation), so distant covariances are noisy. Separable: "
                 "the Gaussian horizontal shape and one length for all levels are assumptions. Both: one day's "
                 "ensemble used for the whole window; unbiased, uncorrelated errors assumed.")


# ---------------------------------------------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="product", required=True)

    def files(parser, name, help_text):
        parser.add_argument(name, type=pathlib.Path, nargs="+", required=True, help=help_text)

    u1 = sub.add_parser("u1", help="reanalysis-ensemble spread")
    files(u1, "--grep", "GREP files from fetch")
    u1.add_argument("--out-dir", type=pathlib.Path, required=True)

    inn = sub.add_parser("innovations", help="U2 (vs reanalysis) or U3 (vs climatology)")
    files(inn, "--argo", "Argo files from fetch")
    inn.add_argument("--background", choices=["grep_mean", "glorys12", "woa23"], required=True)
    files(inn, "--background-files", "files of the background source")
    inn.add_argument("--out", type=pathlib.Path, required=True)
    inn.add_argument("--max-distance-km", type=float, required=True, help="largest pair separation")
    inn.add_argument("--bin-km", type=float, required=True, help="width of the distance bins")
    inn.add_argument("--max-days-apart", type=float, required=True, help="largest time between paired profiles")
    inn.add_argument("--min-pairs", type=int, required=True, help="fewest pairs for a bin or correlation to count")
    inn.add_argument("--map-bin-deg", type=float, required=True, help="box size of the RMS maps")

    u4 = sub.add_parser("u4", help="SST analysis uncertainty")
    files(u4, "--ostia", "OSTIA files from fetch")
    u4.add_argument("--grep", type=pathlib.Path, nargs="+", help="GREP files with mlotst: adds the mixed-layer "
                                                                  "extension as a separate file")
    u4.add_argument("--out", type=pathlib.Path, required=True)

    u5 = sub.add_parser("u5", help="observation-error variance for a target grid")
    files(u5, "--glorys12", "GLORYS12 daily files from fetch")
    u5.add_argument("--grid", type=pathlib.Path, required=True, help="any NetCDF on the target grid")
    u5.add_argument("--out", type=pathlib.Path, required=True)

    u6 = sub.add_parser("u6", help="Argo-driven analysis increment")
    for name in ("--u2", "--u1", "--u5", "--out"):
        u6.add_argument(name, type=pathlib.Path, required=True)
    u6.add_argument("--date", required=True)
    u6.add_argument("--window-days", type=float, required=True,
                    help="observation window centred on the date; GREP members use 5-7 days (QUID 001-031 Table 1)")
    u6.add_argument("--background-model", choices=BACKGROUND_MODELS, required=True,
                    help="background-error covariance (ocean_uncertainty.analysis.ErrorModel)")
    u6.add_argument("--obs-model", choices=OBSERVATION_MODELS, required=True,
                    help="observation-error variance: U5 representativeness, or U2's innovation-based estimate")
    u6.add_argument("--length-km", type=float, help="horizontal correlation or localisation length; "
                                                    "default: median U2 background correlation length")
    u6.add_argument("--profile", type=int, help="use only this U2 profile number")
    u6.add_argument("--block-rows", type=int, help="state points per block, to limit memory (results unchanged)")

    for parser in sub.choices.values():
        parser.add_argument("--variable", default="thetao", help="model variable name")
    a = p.parse_args()

    if a.product == "u1":
        build_u1(a.grep, a.variable, a.out_dir)
    elif a.product == "innovations":
        build_innovations(a.argo, a.background, a.background_files, a.variable, a.out,
                          max_distance_km=a.max_distance_km, bin_km=a.bin_km, max_days_apart=a.max_days_apart,
                          min_pairs=a.min_pairs, map_bin_deg=a.map_bin_deg)
    elif a.product == "u4":
        build_u4(a.ostia, a.out)
        if a.grep:
            build_u4_profile(a.ostia, a.grep, a.out.with_name(a.out.stem + "_mixed_layer.nc"))
    elif a.product == "u5":
        build_u5(a.glorys12, a.grid, a.variable, a.out)
    else:
        build_u6(a.u2, a.u1, a.u5, a.date, a.window_days, a.out, a.variable, a.background_model, a.obs_model,
                 a.length_km, a.profile, a.block_rows)


if __name__ == "__main__":
    main()
