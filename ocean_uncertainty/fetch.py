"""Download the source datasets for a region, a period and a depth range.

Each source goes to its own directory, one file per month (or per day for single dates), and existing files are
skipped, so an interrupted download can be resumed. Run ``python -m ocean_uncertainty.fetch --help``.

A region is (lon_min, lon_max, lat_min, lat_max) in degrees, with lon_max > lon_min; longitudes may run past 180
(e.g. 140 to 240 for the North Pacific).
"""

import argparse
import pathlib

import argopy
import copernicusmarine
import numpy as np
import pandas as pd
import xarray as xr

from ocean_uncertainty import sources
from ocean_uncertainty.regrid import level_interfaces


def periods(start: str, end: str, daily: bool = False) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """Split [start, end] (dates, inclusive) into calendar months, or into single days if ``daily``."""
    start, end = pd.Timestamp(start), pd.Timestamp(end)
    if daily:
        return [(d, d) for d in pd.date_range(start, end, freq="D")]
    firsts = pd.date_range(start.to_period("M").to_timestamp(), end, freq="MS")
    return [(max(f, start), min(f + pd.offsets.MonthEnd(0), end)) for f in firsts]


def levels_reaching(depth: np.ndarray, max_depth: float) -> np.ndarray:
    """The levels whose layers start above ``max_depth``, so a layer mean down to ``max_depth`` is fully covered."""
    return depth[level_interfaces(depth)[:-1] < max_depth]


def fetch_copernicus(dataset_id: str, variables: list[str], region, start, end, out_dir: pathlib.Path,
                     max_depth: float | None = None, daily: bool = False) -> list[pathlib.Path]:
    """Save a Copernicus Marine dataset over ``region`` and [start, end], keeping the levels down to ``max_depth``."""
    out_dir.mkdir(parents=True, exist_ok=True)
    lon_min, lon_max, lat_min, lat_max = region
    written = []
    for first, last in periods(start, end, daily):
        path = out_dir / f"{dataset_id}_{first:%Y%m%d}_{last:%Y%m%d}.nc"
        if not path.exists():
            ds = copernicusmarine.open_dataset(
                dataset_id=dataset_id, variables=variables, minimum_longitude=lon_min, maximum_longitude=lon_max,
                minimum_latitude=lat_min, maximum_latitude=lat_max, start_datetime=f"{first:%Y-%m-%d}",
                end_datetime=f"{last:%Y-%m-%d}")
            if max_depth is not None and "depth" in ds.dims:
                ds = ds.sel(depth=levels_reaching(ds.depth.values, max_depth))
            partial = path.with_suffix(".part")
            ds.to_netcdf(partial)
            partial.rename(path)
        written.append(path)
    return written


def longitude_pieces(lon_min: float, lon_max: float) -> list[tuple[float, float]]:
    """Split a longitude range into pieces within -180 to 180, as argopy requires."""
    if lon_max - lon_min >= 360:
        return [(-180.0, 180.0)]
    lo = (lon_min + 180) % 360 - 180
    hi = lo + (lon_max - lon_min)
    return [(lo, hi)] if hi <= 180 else [(lo, 180.0), (-180.0, hi - 360)]


def fetch_argo(region, start, end, max_pressure: float, out_dir: pathlib.Path, mode: str = "research",
               src: str = "erddap") -> list[pathlib.Path]:
    """Save Argo temperature and salinity points (one file per month) from argopy.

    ``mode="research"`` keeps delayed-mode data with quality flag 1 and pressure error below 20 dbar (argopy
    user-mode documentation). ``max_pressure`` (dbar) should reach below the deepest level to be compared, so
    profiles can be interpolated to it.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    lon_min, lon_max, lat_min, lat_max = region
    argopy.set_options(mode=mode)
    written = []
    for first, last in periods(start, end):
        path = out_dir / f"argo_{mode}_{first:%Y%m%d}_{last:%Y%m%d}.nc"
        if not path.exists():
            parts = []
            for lo, hi in longitude_pieces(lon_min, lon_max):
                # argopy's upper time bound is exclusive (argopy data-selection documentation).
                box = [lo, hi, lat_min, lat_max, 0.0, max_pressure, f"{first:%Y-%m-%d}",
                       f"{last + pd.Timedelta(days=1):%Y-%m-%d}"]
                try:
                    parts.append(argopy.DataFetcher(src=src).region(box).load().data)
                except argopy.errors.DataNotFound:
                    continue
            ds = xr.concat(parts, "N_POINTS") if parts else xr.Dataset()
            # Keep provenance but not the local user name argopy records.
            ds.attrs = {k: v for k, v in (parts[0].attrs if parts else {}).items() if k != "Fetched_by"}
            ds.attrs.update(citation=sources.ARGO["citation"], acknowledgement=sources.ARGO["acknowledgement"],
                            argopy_mode=mode)
            partial = path.with_suffix(".part")
            ds.to_netcdf(partial)
            partial.rename(path)
        written.append(path)
    return written


def select_region(ds: xr.Dataset, region, lon: str = "lon", lat: str = "lat") -> xr.Dataset:
    """Cut ``ds`` to ``region``, rewriting longitudes into the region's own range so it may cross the dateline."""
    lon_min, lon_max, lat_min, lat_max = region
    shifted = (ds[lon].values - lon_min) % 360 + lon_min
    keep = np.flatnonzero(shifted <= lon_max)
    keep = keep[np.argsort(shifted[keep])]
    ds = ds.isel({lon: keep}).assign_coords({lon: shifted[keep]})
    return ds.sel({lat: slice(lat_min, lat_max)})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", choices=["grep", "argo"])
    parser.add_argument("--out", type=pathlib.Path, required=True, help="directory for this source's files")
    parser.add_argument("--region", type=float, nargs=4, metavar=("LON_MIN", "LON_MAX", "LAT_MIN", "LAT_MAX"))
    parser.add_argument("--start", help="first date, YYYY-MM-DD")
    parser.add_argument("--end", help="last date, YYYY-MM-DD (inclusive)")
    parser.add_argument("--daily", action="store_true", help="one file per day instead of per month")
    parser.add_argument("--max-depth", type=float, help="deepest depth (m) the products will use")
    parser.add_argument("--variables", nargs="+", default=["thetao"],
                        help="model variable names (GREP adds the member suffixes)")
    parser.add_argument("--max-pressure", type=float, help="Argo: deepest pressure to fetch (dbar)")
    parser.add_argument("--argo-mode", default="research", choices=["research", "standard", "expert"])
    a = parser.parse_args()

    if a.source == "grep":
        names = [f"{v}_{m}" for v in a.variables for m in sources.GREP["members"]]
        fetch_copernicus(sources.GREP["members_dataset"], names, a.region, a.start, a.end, a.out, a.max_depth,
                         a.daily)
    else:
        fetch_argo(a.region, a.start, a.end, a.max_pressure, a.out, a.argo_mode)


if __name__ == "__main__":
    main()
