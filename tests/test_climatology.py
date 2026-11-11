import datetime

import numpy as np
import xarray as xr

from ocean_uncertainty.climatology import build_c1, sampled_days


def test_sampled_days_follow_the_stated_rule():
    # 5-day step: 1, 6, 11, 16, 21, 26; capped at 28 so February is sampled the same.
    assert [d.day for d in sampled_days(2015, 1, 5)] == [1, 6, 11, 16, 21, 26]
    assert [d.day for d in sampled_days(2016, 2, 5)] == [1, 6, 11, 16, 21, 26]   # leap February, unchanged
    assert [d.day for d in sampled_days(2015, 7, 10)] == [1, 11, 21]


def _spread_file(path, variance_days):
    """A partial-sum file like month_task writes: Sigma of daily variances and the day count."""
    coords = {"depth": [0.0, 10.0], "latitude": [30.0, 31.0], "longitude": [150.0, 151.0]}
    total = sum(variance_days)
    xr.Dataset({"sum_variance": (("depth", "latitude", "longitude"), total), "n_days": len(variance_days)},
               coords=coords).to_netcdf(path)


def test_pooled_sigma_equals_sigma_from_all_days_at_once(tmp_path):
    rng = np.random.default_rng(0)
    shape = (2, 2, 2)
    # Nine January days split across two partial files for the same calendar month.
    days = [rng.uniform(0.01, 0.5, shape) for _ in range(9)]
    _spread_file(tmp_path / "spread_201101.nc", days[:4])
    _spread_file(tmp_path / "spread_201201.nc", days[4:])

    build_c1(sorted(tmp_path.glob("spread_*.nc")), tmp_path / "c1.nc", max_depth=50.0)
    with xr.open_dataset(tmp_path / "c1.nc") as c1:
        pooled = c1["thetao_std"].sel(month=1).values

    expected = np.sqrt(np.mean(days, axis=0))   # sqrt(mean variance) over all nine days at once
    np.testing.assert_allclose(pooled, expected, rtol=1e-5)
