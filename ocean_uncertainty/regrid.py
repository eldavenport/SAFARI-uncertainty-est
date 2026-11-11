"""Conservative regridding of gridded ocean fields onto any target latitude-longitude grid and set of layers.

Products in this repository stay on their source grids. A user moves them to a model grid with :func:`to_grid`,
passing the target cell edges. Averages are taken over the wet (non-NaN) part of each target cell, which is right
for intensive quantities such as temperature: a coastal cell gets the mean of its water.

Note on uncertainty fields: averaging an error *realisation* (e.g. a reanalysis-ensemble anomaly) and then
computing its spread gives the error of the cell mean. Averaging a variance gives the mean point-error variance,
which is larger unless errors are fully correlated across the cell. Regrid realisations when they exist.
"""

import numpy as np
import xarray as xr

LON, LAT, DEPTH = "longitude", "latitude", "depth"


def cell_bounds(centres: np.ndarray) -> np.ndarray:
    """Cell edges for a 1-D axis of cell centres: midpoints between neighbours, and half a step beyond each end.

    Exact for uniform axes. Midpoints rather than a fixed step also absorb the rounding of coordinates stored in
    single precision (e.g. OSTIA latitudes), whose steps differ by ~1e-6 degree.
    """
    centres = np.asarray(centres, dtype=float)
    mid = 0.5 * (centres[1:] + centres[:-1])
    return np.concatenate([[centres[0] - (mid[0] - centres[0])], mid, [centres[-1] + (centres[-1] - mid[-1])]])


def level_interfaces(depth: np.ndarray) -> np.ndarray:
    """Layer edges (m) for values given at level depths: the surface, the midpoints, and half a step below the last.

    This assumes each value represents the layer around it, which is how z-level ocean models such as NEMO store
    temperature. Copernicus Marine files give level depths only, not layer edges.
    """
    depth = np.asarray(depth, dtype=float)
    if np.any(np.diff(depth) <= 0) or depth[0] < 0:
        raise ValueError("depths must be positive and increasing")
    mid = 0.5 * (depth[1:] + depth[:-1])
    bottom = depth[-1] + (depth[-1] - mid[-1]) if depth.size > 1 else 2 * depth[0]
    return np.concatenate([[0.0], mid, [bottom]])


def overlap_matrix(dst_edges: np.ndarray, src_edges: np.ndarray, *, period: float | None = None,
                   measure=lambda x: x) -> np.ndarray:
    """``W[i, j]``: the size of the overlap of target cell ``i`` with source cell ``j``.

    ``measure`` turns a coordinate into the quantity whose differences are cell sizes: the identity for longitude
    and depth, sin(latitude) for latitude (a latitude-longitude cell's area is R^2 dlon d(sin lat)). With
    ``period`` (360 for longitude) source cells are also compared shifted by one period either way, so target
    longitudes may use 0-360 while the source uses -180-180, or cross the dateline.
    """
    d0, d1 = np.asarray(dst_edges, float)[:-1, None], np.asarray(dst_edges, float)[1:, None]
    src_edges = np.asarray(src_edges, float)
    weights = np.zeros((d0.size, src_edges.size - 1))
    for shift in ((-period, 0.0, period) if period else (0.0,)):
        s0, s1 = src_edges[None, :-1] + shift, src_edges[None, 1:] + shift
        lo, hi = np.maximum(d0, s0), np.minimum(d1, s1)
        weights += np.where(hi > lo, measure(hi) - measure(lo), 0.0)
    return weights


def _sin_degrees(latitude):
    return np.sin(np.radians(np.clip(latitude, -90.0, 90.0)))


def to_grid(field: xr.DataArray, lon_edges=None, lat_edges=None, depth_edges=None) -> xr.DataArray:
    """Average ``field`` over the cells of a target grid.

    ``field`` has dimensions named ``longitude``, ``latitude`` and optionally ``depth`` (m, positive down), plus any
    others. Each ``*_edges`` argument gives the target cell edges along that axis (n + 1 values, increasing); an
    axis left as None is not regridded. Source cell edges are taken from :func:`cell_bounds` horizontally and
    :func:`level_interfaces` vertically. The average is over the wet volume of each target cell; cells with no
    water are NaN.
    """
    axes = []
    if depth_edges is not None:
        axes.append((DEPTH, overlap_matrix(depth_edges, level_interfaces(field[DEPTH].values)), depth_edges))
    if lat_edges is not None:
        axes.append((LAT, overlap_matrix(lat_edges, cell_bounds(field[LAT].values), measure=_sin_degrees),
                     lat_edges))
    if lon_edges is not None:
        axes.append((LON, overlap_matrix(lon_edges, cell_bounds(field[LON].values), period=360.0), lon_edges))

    wet = np.isfinite(field.values)
    total, size = np.where(wet, field.values, 0.0), wet.astype(float)
    # Sum the field and the wet volume separately and divide once, so each source cell keeps its own weight.
    for dim, weights, _ in axes:
        axis = field.get_axis_num(dim)
        total = np.moveaxis(np.moveaxis(total, axis, -1) @ weights.T, -1, axis)
        size = np.moveaxis(np.moveaxis(size, axis, -1) @ weights.T, -1, axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        values = np.where(size > 0, total / size, np.nan)

    regridded = {dim for dim, _, _ in axes}
    coords = {k: v for k, v in field.coords.items() if not regridded & set(v.dims)}
    for dim, _, edges in axes:
        edges = np.asarray(edges, float)
        coords[dim] = 0.5 * (edges[1:] + edges[:-1])
    return xr.DataArray(values, dims=field.dims, coords=coords, attrs=field.attrs, name=field.name)


def subgrid_covariance(a: xr.DataArray, b: xr.DataArray, lon_edges, lat_edges, depth_edges=None) -> xr.DataArray:
    """Covariance of two high-resolution fields within each target cell: mean of a*b minus product of the means.

    With b = a this is the within-cell variance: the representativeness error of a point observation compared with
    a cell average, i.e. the part of the point's variability a field on the target grid cannot hold. With a and b
    at two different depths it is the covariance of that error between the depths. Only points wet in both
    fields count. Variability finer than the source grid is missing, so this is a lower bound.
    """
    a, b = a.where(b.notnull()), b.where(a.notnull())
    mean_product = to_grid(a * b, lon_edges, lat_edges, depth_edges)
    return mean_product - to_grid(a, lon_edges, lat_edges, depth_edges) * to_grid(b, lon_edges, lat_edges,
                                                                                  depth_edges)


def subgrid_variance(field: xr.DataArray, lon_edges, lat_edges, depth_edges=None) -> xr.DataArray:
    """Within-cell variance of ``field`` (see :func:`subgrid_covariance`); rounding below zero is clipped."""
    return subgrid_covariance(field, field, lon_edges, lat_edges, depth_edges).clip(min=0.0)
