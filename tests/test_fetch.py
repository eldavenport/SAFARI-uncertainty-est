import numpy as np
import pandas as pd

from ocean_uncertainty.fetch import levels_reaching, longitude_pieces, periods


def test_longitude_pieces():
    assert longitude_pieces(140, 240) == [(140, 180.0), (-180.0, -120)]
    assert longitude_pieces(-170, -160) == [(-170, -160)]
    assert longitude_pieces(170, 180) == [(170, 180)]
    assert longitude_pieces(200, 220) == [(-160, -140)]
    assert longitude_pieces(0, 360) == [(-180.0, 180.0)]


def test_periods_split_on_month_boundaries():
    got = periods("2011-11-15", "2012-01-03")
    assert got == [(pd.Timestamp("2011-11-15"), pd.Timestamp("2011-11-30")),
                   (pd.Timestamp("2011-12-01"), pd.Timestamp("2011-12-31")),
                   (pd.Timestamp("2012-01-01"), pd.Timestamp("2012-01-03"))]
    assert len(periods("2015-01-30", "2015-02-02", daily=True)) == 4


def test_levels_reaching_include_the_layer_across_the_bottom():
    # Levels 5, 15, 25, 35 m have layers starting at 0, 10, 20, 30 m. A 0-22 m mean needs the third layer too.
    np.testing.assert_array_equal(levels_reaching(np.array([5.0, 15.0, 25.0, 35.0]), 22.0), [5.0, 15.0, 25.0])


def test_select_region_across_the_dateline():
    import xarray as xr
    from ocean_uncertainty.fetch import select_region
    lon = np.arange(-179.5, 180, 1.0)
    ds = xr.Dataset({"t": ("lon", lon)}, coords={"lon": lon, "lat": [0.0]})
    out = select_region(ds, (170.0, 190.0, -1.0, 1.0))
    np.testing.assert_array_equal(out.lon.values, np.arange(170.5, 190, 1.0))
    # Values travel with their longitudes: the cell at 189.5 E is the source cell at -170.5.
    assert float(out.t.sel(lon=189.5)) == -170.5
