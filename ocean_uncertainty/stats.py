"""Statistics behind the products: ensemble spread, error covariances from observation pairs, and the analysis update.

Notation: d = y - H x_b is the innovation, B and R the background- and observation-error covariances.
"""

from collections.abc import Iterable

import numpy as np
import xarray as xr
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree

# Mean Earth radius R1 of the Geodetic Reference System 1980 (Moritz 2000, J. Geodesy 74, 128-133).
EARTH_RADIUS_KM = 6371.0088


def ensemble_spread(members: xr.DataArray, dim: str = "member") -> tuple[xr.DataArray, xr.DataArray]:
    """Anomalies from the ensemble mean, and the standard deviation with N - 1 in the denominator.

    N - 1 makes the variance unbiased for independent members; with the three GREP members this matters (a
    factor 1.5 in variance compared with dividing by N).
    """
    anomalies = members - members.mean(dim)
    return anomalies, members.std(dim, ddof=1)


def unit_vectors(lon, lat) -> np.ndarray:
    """Points on the unit sphere, shape (n, 3), from longitude and latitude in degrees."""
    lon, lat = np.radians(np.asarray(lon, float)), np.radians(np.asarray(lat, float))
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def distance_km(lon_a, lat_a, lon_b, lat_b) -> np.ndarray:
    """Great-circle distances (km) between every point of set a and every point of set b, shape (n_a, n_b)."""
    cos_angle = unit_vectors(lon_a, lat_a) @ unit_vectors(lon_b, lat_b).T
    return EARTH_RADIUS_KM * np.arccos(np.clip(cos_angle, -1.0, 1.0))


def gaspari_cohn(distance, half_width) -> np.ndarray:
    """The compactly supported correlation function of Gaspari & Cohn (1999, QJRMS 125, 723-757, eq. 4.10).

    It is 1 at zero distance and 0 beyond 2 * ``half_width``. Its curvature at zero equals that of a Gaussian
    exp(-r^2 / 2L^2) when half_width = L * sqrt(10/3) (see :func:`gaspari_cohn_half_width`).
    """
    z = np.abs(np.asarray(distance, float)) / half_width
    inner = ((-0.25 * z + 0.5) * z + 0.625) * z**3 - 5.0 / 3.0 * z**2 + 1.0
    with np.errstate(divide="ignore"):
        outer = ((((z / 12.0 - 0.5) * z + 0.625) * z + 5.0 / 3.0) * z - 5.0) * z + 4.0 - 2.0 / (3.0 * z)
    return np.where(z <= 1.0, inner, np.where(z < 2.0, outer, 0.0))


def gaspari_cohn_half_width(length):
    """Gaspari-Cohn half-width matching a Gaussian of length ``length`` near zero distance.

    Near zero the Gaussian is 1 - r^2 / (2 L^2) and Gaspari-Cohn is 1 - 5 r^2 / (3 c^2); equal when c^2 = 10 L^2 / 3.
    """
    return np.sqrt(10.0 / 3.0) * length


def pair_covariance(lon, lat, time, values, *, max_distance_km: float, bin_km: float,
                    max_time_apart: np.timedelta64) -> dict:
    """Covariance of ``values`` between pairs of points, binned by their distance apart.

    Used on innovations at one depth: if observation errors are uncorrelated between different profiles, the
    covariance of innovations at distance r > 0 is the background-error covariance, while the variance (r = 0)
    also contains the observation-error variance (Hollingsworth & Lonnberg 1986, Tellus 38A, 111-136).

    Pairs are profiles at most ``max_distance_km`` apart and ``max_time_apart`` apart in time. The sample mean
    (the bias) is removed first and returned. Returns bin centres (km), the mean product in each bin, the number of
    pairs per bin, the variance, the bias and the number of points.
    """
    values = np.asarray(values, float)
    ok = np.isfinite(values)
    lon, lat, time, values = (np.asarray(a)[ok] for a in (lon, lat, time, values))
    bias = values.mean()
    anomaly = values - bias

    # On the unit sphere, an arc of angle a is a chord of 2 sin(a / 2).
    chord = 2.0 * np.sin(max_distance_km / EARTH_RADIUS_KM / 2.0)
    pairs = cKDTree(unit_vectors(lon, lat)).query_pairs(chord, output_type="ndarray")
    i, j = pairs[:, 0], pairs[:, 1]
    keep = np.abs(time[i] - time[j]) <= max_time_apart
    i, j = i[keep], j[keep]

    edges = np.arange(0.0, max_distance_km + bin_km, bin_km)
    dist = EARTH_RADIUS_KM * np.arccos(np.clip(np.sum(unit_vectors(lon[i], lat[i]) * unit_vectors(lon[j], lat[j]),
                                                      axis=1), -1.0, 1.0))
    which = np.digitize(dist, edges) - 1
    inside = (which >= 0) & (which < edges.size - 1)
    counts = np.bincount(which[inside], minlength=edges.size - 1)
    sums = np.bincount(which[inside], weights=(anomaly[i] * anomaly[j])[inside], minlength=edges.size - 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        covariance = np.where(counts > 0, sums / counts, np.nan)
    return dict(distance_km=0.5 * (edges[1:] + edges[:-1]), covariance=covariance, pairs=counts,
                variance=float(np.mean(anomaly**2)), bias=float(bias), points=int(values.size))


def fit_background_and_observation_error(binned: dict, min_pairs: int) -> dict:
    """Split innovation variance into background and observation parts, following Hollingsworth & Lonnberg (1986).

    A Gaussian, s_b^2 exp(-r^2 / 2L^2), is fitted to the pair covariances at r > 0 (bins with at least ``min_pairs``
    pairs, weighted by their pair counts). Its value at r = 0 is the background-error variance s_b^2. The rest of
    the innovation variance is the observation-error variance (instrument plus representativeness):
    sigma_o^2 = variance - s_b^2. The Gaussian shape is an assumption; the binned values are returned with the fit
    so it can be checked.
    """
    use = (binned["pairs"] >= min_pairs) & np.isfinite(binned["covariance"])
    r, c, n = binned["distance_km"][use], binned["covariance"][use], binned["pairs"][use]
    if r.size < 3:
        raise ValueError(f"only {r.size} distance bins have at least {min_pairs} pairs; need 3 to fit")

    def gaussian(r, variance, length):
        return variance * np.exp(-0.5 * (r / length) ** 2)

    if c[0] <= 0:
        raise ValueError("the nearest distance bin has no positive covariance: no correlated background error to fit")
    # Starting guess: the first bin's covariance, and the distance where the covariance falls below half of it.
    below = np.flatnonzero(c < 0.5 * c[0])
    start = (c[0], r[below[0]] if below.size else r[-1])
    (variance_b, length), cov = curve_fit(gaussian, r, c, p0=start, sigma=1.0 / np.sqrt(n), absolute_sigma=False,
                                          bounds=([0.0, 1e-6], [np.inf, np.inf]))
    return dict(background_variance=float(variance_b), length_km=float(abs(length)),
                observation_variance=float(binned["variance"] - variance_b),
                background_variance_se=float(np.sqrt(cov[0, 0])), length_km_se=float(np.sqrt(cov[1, 1])))


def ensemble_covariance(anomalies_a: np.ndarray, anomalies_b: np.ndarray, localisation: np.ndarray) -> np.ndarray:
    """Ensemble covariance between two sets of points, (n_a, N) and (n_b, N) anomalies, times a localisation.

    With few members the raw covariance has spurious values far from the observation; multiplying by a
    correlation function that goes to zero with distance removes them (Houtekamer & Mitchell 2001, Mon. Wea. Rev.
    129, 123-137).
    """
    n_members = anomalies_a.shape[1]
    return (anomalies_a @ anomalies_b.T) / (n_members - 1) * localisation


def analysis_update(innovation: np.ndarray, obs_error_cov: np.ndarray, cov_obs_obs: np.ndarray,
                    cov_state_obs: Iterable[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """The optimal interpolation (OI) update, delta_x = B H^T (H B H^T + R)^-1 d, and the variance it removes.

    ``cov_obs_obs`` is H B H^T (n_obs, n_obs), ``obs_error_cov`` is R (n_obs, n_obs), and ``cov_state_obs``
    yields B H^T in blocks of state points (each (n_block, n_obs)), so a large state never has to be held at once.
    Returns delta_x and diag(B H^T (H B H^T + R)^-1 H B); the analysis-error variance is diag(B) minus the latter.
    Assumes unbiased, uncorrelated background and observation errors and a linear H (Kalnay 2003, ch. 5).
    """
    S = cov_obs_obs + obs_error_cov
    factor = cho_factor(S, lower=True)
    weights = cho_solve(factor, innovation)
    increments, reductions = [], []
    for block in cov_state_obs:
        increments.append(block @ weights)
        # diag(C S^-1 C^T) = column sums of (L^-1 C^T)^2, where S = L L^T.
        reductions.append(np.sum(solve_triangular(factor[0], block.T, lower=True) ** 2, axis=0))
    return np.concatenate(increments), np.concatenate(reductions)


def level_correlation(anomalies: np.ndarray) -> np.ndarray:
    """Correlation between levels of ensemble anomalies, (n_levels, n_samples), pooled over all samples.

    Anomalies are deviations from the ensemble mean, so they are not centred again. The sum is taken in double
    precision: in single precision (as the files store anomalies) rounding makes the matrix slightly indefinite.
    """
    anomalies = np.asarray(anomalies, dtype=np.float64)
    cov = anomalies @ anomalies.T
    scale = np.sqrt(np.diag(cov))
    return cov / np.outer(scale, scale)


def separable_covariance(sigma_a: np.ndarray, sigma_b: np.ndarray, distance: np.ndarray, length_km: float,
                         level_corr: np.ndarray) -> np.ndarray:
    """sigma_a sigma_b exp(-r^2 / 2L^2) times a correlation between the two points' levels.

    A full-rank covariance model: variances and the vertical structure from the ensemble, the horizontal
    correlation a Gaussian of length L (e.g. fitted to Argo innovations).
    """
    return np.outer(sigma_a, sigma_b) * np.exp(-0.5 * (distance / length_km) ** 2) * level_corr
