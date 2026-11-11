import gsw
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from ocean_uncertainty.observations import day_index, month_weights, sample, to_levels


def profile_points(pressures, temperatures, lat=35.0):
    n = len(pressures)
    return pd.DataFrame({"PLATFORM_NUMBER": [1] * n, "CYCLE_NUMBER": [7] * n, "DIRECTION": ["A"] * n,
                         "time": [pd.Timestamp("2015-01-10 06:00")] * n, "lon": [-160.0] * n, "lat": [lat] * n,
                         "PRES": pressures, "value": temperatures})


def test_profile_is_interpolated_in_depth_without_extrapolation():
    pressures = np.array([4.0, 10.0, 20.0, 40.0])
    depth = -gsw.z_from_p(pressures, 35.0)
    # Temperature linear in depth, so linear interpolation is exact.
    points = profile_points(pressures, 18.0 - 0.05 * depth)
    levels = np.array([0.5, 5.0, 15.0, 30.0, 45.0])
    table = to_levels(points, levels)
    # 0.5 m lies above the shallowest measurement and 45 m below the deepest: both are left out.
    np.testing.assert_allclose(table["depth"], [5.0, 15.0, 30.0])
    np.testing.assert_allclose(table["obs"], 18.0 - 0.05 * np.array([5.0, 15.0, 30.0]))
    np.testing.assert_allclose(table["gap_m"], np.diff(depth))


def grid_field(values_fn, lon, lat, depth, times):
    lo, la = np.meshgrid(lon, lat)
    data = np.stack([np.stack([values_fn(lo, la) + 0.1 * d + t for d in range(depth.size)])
                     for t in range(times.size)])
    return xr.DataArray(data, dims=("time", "depth", "latitude", "longitude"),
                        coords=dict(time=times, depth=depth, latitude=lat, longitude=lon))


def test_bilinear_sampling_is_exact_for_a_bilinear_field():
    lon, lat = np.arange(140.0, 250.0, 0.25), np.arange(20.0, 40.0, 0.25)
    times = pd.date_range("2015-01-01", periods=3).to_numpy()
    field = grid_field(lambda x, y: 0.3 * x - 0.2 * y, lon, lat, np.array([0.5, 1.5]), times)
    # -160 E is 200 on this grid's longitudes.
    got = sample(field, lon=[-160.0, 145.1], lat=[30.13, 21.7], depth=[1.5, 0.5], time_index=[2, 0])
    np.testing.assert_allclose(got, [0.3 * 200 - 0.2 * 30.13 + 0.1 + 2, 0.3 * 145.1 - 0.2 * 21.7], rtol=1e-12)


def test_sampling_next_to_land_or_outside_is_missing():
    lon, lat = np.arange(0.0, 5.0), np.arange(0.0, 5.0)
    times = pd.date_range("2015-01-01", periods=1).to_numpy()
    field = grid_field(lambda x, y: x + y, lon, lat, np.array([1.0]), times)
    field[0, 0, 2, 2] = np.nan
    got = sample(field, lon=[2.5, 0.5, 10.0, 0.5], lat=[2.5, 0.5, 1.0, 0.5], depth=[1.0] * 4,
                 time_index=[0, 0, 0, -1])
    assert np.isnan(got[0]) and np.isfinite(got[1]) and np.isnan(got[2]) and np.isnan(got[3])


def test_day_index_matches_daily_means():
    field_time = pd.date_range("2015-01-01", periods=31).to_numpy()
    obs = pd.to_datetime(["2015-01-01 23:59", "2015-01-31 00:01", "2015-02-01 01:00"])
    np.testing.assert_array_equal(day_index(field_time, obs), [0, 30, -1])


def test_month_weights_between_mid_months():
    # Mid-December is 16 Dec 12:00 and mid-January 16 Jan 12:00 (31 days apart); 1 Jan 00:00 is 15.5 days on.
    first, second, w = month_weights(pd.to_datetime(["2015-01-01 00:00", "2015-01-16 12:00", "2015-06-30 00:00"]))
    np.testing.assert_array_equal(first, [12, 1, 6])
    np.testing.assert_array_equal(second, [1, 2, 7])
    assert w[0] == pytest.approx(15.5 / 31)
    assert w[1] == pytest.approx(0.0)
    # Mid-June is 16 Jun 00:00, mid-July 16 Jul 12:00: 14 days of 30.5.
    assert w[2] == pytest.approx(14 / 30.5)
