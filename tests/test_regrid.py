import numpy as np
import xarray as xr

from ocean_uncertainty.regrid import cell_bounds, level_interfaces, subgrid_variance, to_grid


def field_on(lon, lat, values, depth=None):
    dims = ("latitude", "longitude") if depth is None else ("depth", "latitude", "longitude")
    coords = dict(latitude=lat, longitude=lon) | ({} if depth is None else dict(depth=depth))
    return xr.DataArray(values, dims=dims, coords=coords)


def cell_areas(lon, lat):
    """Exact areas (unit sphere) of uniform lat-lon cells, from their edges."""
    dlon = np.radians(np.diff(cell_bounds(lon)))
    dsin = np.diff(np.sin(np.radians(cell_bounds(lat))))
    return dsin[:, None] * dlon[None, :]


def test_area_integral_is_conserved():
    rng = np.random.default_rng(0)
    lon, lat = np.arange(0.125, 10, 0.25), np.arange(20.125, 30, 0.25)
    src = field_on(lon, lat, rng.normal(size=(lat.size, lon.size)))
    lon_edges, lat_edges = np.arange(0, 10.1, 1.0), np.arange(20, 30.1, 1.0)
    dst = to_grid(src, lon_edges, lat_edges)
    src_total = np.sum(src.values * cell_areas(lon, lat))
    dst_total = np.sum(dst.values * cell_areas(dst.longitude.values, dst.latitude.values))
    np.testing.assert_allclose(dst_total, src_total, rtol=1e-12)


def test_coastal_cell_is_mean_of_its_wet_area():
    # One 1x1 degree target cell over a 2x2 source block. Wet: three cells, values 1 (dry neighbour), 2 and 4.
    # The cells have (nearly) equal area, so the answer is the area-weighted mean of the three wet values only.
    # Averaging one axis at a time would instead give the latitude strips equal weight.
    lon, lat = np.array([0.25, 0.75]), np.array([0.25, 0.75])
    values = np.array([[1.0, np.nan], [2.0, 4.0]])
    dst = to_grid(field_on(lon, lat, values), np.array([0.0, 1.0]), np.array([0.0, 1.0]))
    area = cell_areas(lon, lat)
    expected = np.nansum(values * area) / area[np.isfinite(values)].sum()
    np.testing.assert_allclose(dst.values[0, 0], expected, rtol=1e-12)


def test_target_longitudes_may_cross_the_dateline():
    # Source on -180..180; the target cell 179-181 E straddles the seam, so it takes source cells from both ends.
    lon = np.arange(-179.5, 180, 1.0)
    lat = np.array([0.5])
    src = field_on(lon, lat, np.where(lon > 0, 1.0, 3.0)[None, :])
    dst = to_grid(src, np.array([179.0, 181.0]), None)
    np.testing.assert_allclose(dst.values, [[2.0]])


def test_vertical_layer_mean_is_thickness_weighted():
    # Levels at 5, 15, 25 m have layers 0-10, 10-20, 20-30 m (midpoints; level_interfaces). A target layer of
    # 5-30 m covers half of the first layer and all of the other two: (1*5 + 2*10 + 4*10) / 25.
    depth = np.array([5.0, 15.0, 25.0])
    np.testing.assert_allclose(level_interfaces(depth), [0, 10, 20, 30])
    src = field_on(np.array([0.5]), np.array([0.5]), np.array([1.0, 2.0, 4.0])[:, None, None], depth)
    dst = to_grid(src, depth_edges=np.array([5.0, 30.0]))
    np.testing.assert_allclose(dst.values.ravel(), [(5 + 20 + 40) / 25])


def test_subgrid_variance_of_a_field_linear_in_longitude():
    # n equal-width source cells across a target cell, values a + s*x_i at spacing h: the variance of n equally
    # spaced, equally weighted values is s^2 h^2 (n^2 - 1) / 12.
    n, h, s = 4, 0.25, 0.8
    lon = np.arange(h / 2, 1.0, h)
    lat = np.array([10.0, 10.25])
    src = field_on(lon, lat, np.broadcast_to(3.0 + s * lon, (lat.size, lon.size)))
    var = subgrid_variance(src, np.array([0.0, 1.0]), None)
    np.testing.assert_allclose(var.values, s**2 * h**2 * (n**2 - 1) / 12, rtol=1e-12)


def test_subgrid_variance_is_zero_on_the_source_grid():
    rng = np.random.default_rng(1)
    lon, lat = np.arange(0.5, 5), np.arange(0.5, 4)
    src = field_on(lon, lat, rng.normal(size=(lat.size, lon.size)))
    var = subgrid_variance(src, cell_bounds(lon), cell_bounds(lat))
    np.testing.assert_allclose(var.values, 0.0, atol=1e-12)


def test_subgrid_covariance_between_two_fields():
    # b = 2a + c (c constant within the cell): cov(a, b) = 2 var(a). A point dry in b is left out of both.
    from ocean_uncertainty.regrid import subgrid_covariance
    lon, lat = np.array([0.25, 0.75]), np.array([0.25, 0.75])
    a_values = np.array([[1.0, 2.0], [4.0, 8.0]])
    a = field_on(lon, lat, a_values)
    b = field_on(lon, lat, 2 * a_values + 5.0)
    edges = np.array([0.0, 1.0])
    area = cell_areas(lon, lat)
    w = area / area.sum()
    var_a = np.sum(w * a_values**2) - np.sum(w * a_values) ** 2
    np.testing.assert_allclose(subgrid_covariance(a, b, edges, edges).values, [[2 * var_a]], rtol=1e-12)

    b_dry = b.copy()
    b_dry[1, 1] = np.nan
    wet = np.isfinite(b_dry.values)
    w = np.where(wet, area, 0.0) / area[wet].sum()
    mean_a = np.sum(w * a_values)
    var_wet = np.sum(w * a_values**2) - mean_a**2
    np.testing.assert_allclose(subgrid_covariance(a, b_dry, edges, edges).values, [[2 * var_wet]], rtol=1e-12)
