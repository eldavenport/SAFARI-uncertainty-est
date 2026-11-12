import numpy as np
import pandas as pd
import pytest
import xarray as xr

from ocean_uncertainty.increment import increment


def _tiny_products(sigma_b=0.4, obs_var=0.09, length_km=50.0):
    lon = np.arange(0.0, 7.0)          # ~111 km per degree: clusters 2 deg apart are independent at L = 50 km
    lat = np.arange(0.0, 3.0)
    sb = np.full((1, 1, lat.size, lon.size), sigma_b)
    c1 = xr.Dataset({"thetao_std": (("month", "depth", "latitude", "longitude"), sb)},
                    coords={"month": [1], "depth": [10.0], "latitude": lat, "longitude": lon})
    c2 = xr.Dataset({"length_km": (("background", "month", "depth"), [[[length_km]]]),
                     "observation_variance": (("background", "month", "depth"), [[[obs_var]]])},
                    coords={"background": ["mean"], "month": [1], "depth": [10.0]})
    return c1, c2


def _table(d=1.0):
    rows = []
    for lon in (1.0, 3.0, 5.0):        # three separated clusters, 4 identical profiles each
        for _ in range(4):
            rows.append(dict(depth=10.0, lon=lon, lat=1.0, innov_mean=d))
    return pd.DataFrame(rows)


def test_increment_is_the_textbook_gain_and_never_overshoots():
    sigma_b, obs_var, d = 0.4, 0.09, 1.0
    c1, c2 = _tiny_products(sigma_b, obs_var)
    ds = increment(c1, c2, _table(d), month=1)
    inc = ds["increment"].sel(depth=10.0)

    gain = sigma_b**2 / (sigma_b**2 + obs_var)            # single co-located obs: dx = gain * d
    assert inc.sel(latitude=1.0, longitude=1.0).item() == pytest.approx(gain * d, abs=0.03)
    # No overshoot anywhere, and the analysis is more certain than the background at the obs.
    assert np.nanmax(np.abs(inc.values)) <= d + 1e-9
    assert ds["analysis_std"].sel(depth=10.0, latitude=1.0, longitude=1.0).item() < sigma_b
    # Far from every float the increment vanishes.
    assert abs(inc.sel(latitude=0.0, longitude=0.0).item()) < 0.05
