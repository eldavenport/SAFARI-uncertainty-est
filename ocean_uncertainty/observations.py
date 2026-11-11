"""Argo profiles on a gridded product's levels, and gridded fields sampled at the profiles.

Comparing an observation with a gridded field needs both at the same place: profiles are interpolated in depth to
the field's levels (never extrapolated), and the field is interpolated bilinearly in the horizontal to each
profile's position, on the profile's day.
"""

import pathlib

import gsw
import numpy as np
import pandas as pd
import xarray as xr

PROFILE_KEY = ["PLATFORM_NUMBER", "CYCLE_NUMBER", "DIRECTION"]


def load_argo(files: list[pathlib.Path], variable: str = "TEMP") -> pd.DataFrame:
    """Argo points from :func:`ocean_uncertainty.fetch.fetch_argo` files, one row per measurement."""
    frames = []
    for path in files:
        with xr.open_dataset(path) as ds:
            if "N_POINTS" not in ds.dims:
                continue
            frames.append(ds[PROFILE_KEY + ["PRES", variable]].to_dataframe().reset_index(drop=True))
    points = pd.concat(frames, ignore_index=True)
    return points.rename(columns={"LONGITUDE": "lon", "LATITUDE": "lat", "TIME": "time", variable: "value"})


def to_levels(points: pd.DataFrame, levels: np.ndarray) -> pd.DataFrame:
    """Interpolate each profile linearly in depth to ``levels`` (m) lying between its shallowest and deepest point.

    Pressure is converted to depth with TEOS-10 (gsw.z_from_p, at the profile's latitude). The ``gap_m`` column
    gives the distance between the two measurements used for each value, so users can drop values interpolated
    across large gaps.
    """
    points = points.assign(depth=-gsw.z_from_p(points["PRES"].to_numpy(), points["lat"].to_numpy()))
    rows = []
    for key, p in points.groupby(PROFILE_KEY, sort=False):
        p = p.drop_duplicates("depth").sort_values("depth")
        depth, value = p["depth"].to_numpy(), p["value"].to_numpy()
        inside = levels[(levels >= depth[0]) & (levels <= depth[-1])]
        if inside.size == 0:
            continue
        upper = np.clip(np.searchsorted(depth, inside), 1, depth.size - 1)
        rows.append(pd.DataFrame({
            "platform": key[0], "cycle": key[1], "direction": key[2],
            "time": p["time"].iloc[0], "lon": p["lon"].iloc[0], "lat": p["lat"].iloc[0],
            "depth": inside, "obs": np.interp(inside, depth, value),
            "gap_m": depth[upper] - depth[upper - 1],
        }))
    table = pd.concat(rows, ignore_index=True)
    table["profile"] = table.groupby(["platform", "cycle", "direction"], sort=False).ngroup()
    return table


def sample(field: xr.DataArray, lon, lat, depth, time_index) -> np.ndarray:
    """Bilinear horizontal interpolation of ``field`` (time, depth, latitude, longitude) at each observation.

    ``depth`` values must be levels of the field and ``time_index`` the position of each observation's time on the
    field's time axis (-1 where the field does not cover it). Observations outside the grid, or next to a dry (NaN) grid point, get NaN: near coasts the
    gridded value is not representative of the open-water measurement.
    """
    lon_axis, lat_axis = field["longitude"].values, field["latitude"].values
    lon = (np.asarray(lon, float) - lon_axis[0]) % 360 + lon_axis[0]
    lat = np.asarray(lat, float)
    k = np.searchsorted(field["depth"].values, depth)
    k = np.clip(k, 0, field["depth"].size - 1)
    if not np.allclose(field["depth"].values[k], depth):
        raise ValueError("observation depths must be levels of the field")

    i = np.clip(np.searchsorted(lon_axis, lon) - 1, 0, lon_axis.size - 2)
    j = np.clip(np.searchsorted(lat_axis, lat) - 1, 0, lat_axis.size - 2)
    wx = (lon - lon_axis[i]) / (lon_axis[i + 1] - lon_axis[i])
    wy = (lat - lat_axis[j]) / (lat_axis[j + 1] - lat_axis[j])
    values = field.transpose("time", "depth", "latitude", "longitude").values
    t = np.asarray(time_index)
    out = ((1 - wx) * (1 - wy) * values[t, k, j, i] + wx * (1 - wy) * values[t, k, j, i + 1]
           + (1 - wx) * wy * values[t, k, j + 1, i] + wx * wy * values[t, k, j + 1, i + 1])
    outside = (wx < 0) | (wx > 1) | (wy < 0) | (wy > 1) | (t < 0)
    return np.where(outside, np.nan, out)


def day_index(field_time: np.ndarray, obs_time) -> np.ndarray:
    """Position of each observation's day on a daily-mean time axis (Copernicus daily means are stamped 00:00)."""
    days = pd.DatetimeIndex(obs_time).floor("D").to_numpy()
    index = np.searchsorted(field_time, days)
    found = (index < field_time.size) & (field_time[np.minimum(index, field_time.size - 1)] == days)
    return np.where(found, index, -1)


def month_weights(obs_time) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Linear interpolation in time between monthly climatological means placed at mid-month.

    Returns, for each observation, the two month numbers (1-12) and the weight of the second. A monthly mean is
    taken to represent the middle of its month; between two mid-months the climatology is assumed to change
    linearly.
    """
    t = pd.DatetimeIndex(obs_time)

    def middle(months: pd.PeriodIndex) -> pd.DatetimeIndex:
        return months.start_time + ((months + 1).start_time - months.start_time) / 2

    current = t.to_period("M")
    first = pd.PeriodIndex(np.where(t < middle(current), current - 1, current), freq="M")
    weight = (t - middle(first)) / (middle(first + 1) - middle(first))
    return np.asarray(first.month), np.asarray((first + 1).month), np.asarray(weight, float)
