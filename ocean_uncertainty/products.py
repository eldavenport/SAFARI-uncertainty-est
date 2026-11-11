"""Shared helpers for building the error-climatology products (ocean_uncertainty.climatology).

``write`` saves a product with provenance; ``reanalysis_members`` loads the GREP members of a file;
``innovation_statistics`` reduces a matched Argo-minus-background table to per-level and mapped statistics.
"""

import datetime
import pathlib

import numpy as np
import pandas as pd
import xarray as xr

from ocean_uncertainty import sources
from ocean_uncertainty.stats import fit_background_and_observation_error, pair_covariance


def write(ds: xr.Dataset, path: pathlib.Path, product: str, citations: list[str], processing: str,
          limitations: str) -> pathlib.Path:
    """Save a product with its provenance in the global attributes."""
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.attrs.update(product=product, sources=" | ".join(citations), processing=processing,
                    limitations=limitations,
                    created=datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"))
    # Single precision is ample for storage (computations that need double precision convert first); compressed.
    encoding = {name: {"dtype": "float32", "zlib": True, "complevel": 1}
                for name, var in ds.data_vars.items() if np.issubdtype(var.dtype, np.floating)}
    partial = path.with_suffix(".part")
    ds.to_netcdf(partial, encoding=encoding)
    partial.rename(path)
    return path


def reanalysis_members(path: pathlib.Path, variable: str) -> xr.DataArray:
    """The GREP members of ``variable`` in one downloaded file, stacked along a ``member`` dimension."""
    with xr.open_dataset(path) as ds:
        names = list(sources.GREP["members"])
        return xr.concat([ds[f"{variable}_{m}"].load() for m in names], dim=pd.Index(names, name="member"))


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
