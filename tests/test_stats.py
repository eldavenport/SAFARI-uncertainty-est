import numpy as np
import pytest
import xarray as xr

from ocean_uncertainty.stats import (EARTH_RADIUS_KM, analysis_update, distance_km, ensemble_spread,
                                     fit_background_and_observation_error, gaspari_cohn, gaspari_cohn_half_width,
                                     pair_covariance)


def test_ensemble_spread_matches_sample_covariance():
    rng = np.random.default_rng(0)
    members = xr.DataArray(rng.normal(size=(3, 7)), dims=("member", "x"))
    anomalies, sigma = ensemble_spread(members)
    np.testing.assert_allclose(sigma.values**2, np.diag(np.cov(members.values, rowvar=False)), rtol=1e-12)
    np.testing.assert_allclose(anomalies.sum("member").values, 0.0, atol=1e-12)


def test_distance_quarter_meridian():
    # Equator to pole is a quarter of a great circle.
    np.testing.assert_allclose(distance_km([0.0], [0.0], [0.0], [90.0]), [[np.pi / 2 * EARTH_RADIUS_KM]])


def test_gaspari_cohn_shape():
    c = 100.0
    assert gaspari_cohn(0.0, c) == pytest.approx(1.0)
    assert gaspari_cohn(2 * c, c) == pytest.approx(0.0, abs=1e-12)
    assert gaspari_cohn(3 * c, c) == 0.0
    # Continuous where the two polynomial pieces meet (the published function is continuous there).
    assert gaspari_cohn(c * (1 - 1e-9), c) == pytest.approx(gaspari_cohn(c * (1 + 1e-9), c), abs=1e-7)
    # Its curvature at zero equals that of the Gaussian with the matching length.
    length = 70.0
    h = 1e-3
    gc = gaspari_cohn(np.array([-h, 0.0, h]), gaspari_cohn_half_width(length))
    gauss = np.exp(-0.5 * (np.array([-h, 0.0, h]) / length) ** 2)
    np.testing.assert_allclose(np.diff(gc, 2) / h**2, np.diff(gauss, 2) / h**2, rtol=1e-4)


def synthetic_innovations(seed, n_points, n_times, sigma_b, length_km, sigma_o):
    """Innovations = a Gaussian-correlated background error plus white observation error, at random points.

    Each time is a separate draw, 10 days apart, so pairs are only formed within one draw.
    """
    rng = np.random.default_rng(seed)
    out = {k: [] for k in ("lon", "lat", "time", "d")}
    for t in range(n_times):
        lon, lat = rng.uniform(160, 200, n_points), rng.uniform(30, 50, n_points)
        cov = sigma_b**2 * np.exp(-0.5 * (distance_km(lon, lat, lon, lat) / length_km) ** 2)
        background = np.linalg.cholesky(cov + 1e-10 * np.eye(n_points)) @ rng.normal(size=n_points)
        out["lon"].append(lon)
        out["lat"].append(lat)
        out["time"].append(np.full(n_points, np.datetime64("2015-01-01") + np.timedelta64(10 * t, "D")))
        out["d"].append(background + sigma_o * rng.normal(size=n_points))
    return {k: np.concatenate(v) for k, v in out.items()}


def test_background_and_observation_errors_are_recovered():
    truth = dict(sigma_b=0.4, length_km=300.0, sigma_o=0.2)
    x = synthetic_innovations(seed=0, n_points=500, n_times=30, **truth)
    binned = pair_covariance(x["lon"], x["lat"], x["time"], x["d"], max_distance_km=1200.0, bin_km=50.0,
                             max_time_apart=np.timedelta64(0, "D"))
    fit = fit_background_and_observation_error(binned, min_pairs=100)
    # Tolerances are 3 standard deviations of each estimate over 20 seeds of this set-up (seeds 0-19, run
    # 2026-11-10): 2.2 % for sigma_b, 5.6 % for the length, 5.2 % for sigma_o.
    assert np.sqrt(fit["background_variance"]) == pytest.approx(truth["sigma_b"], rel=0.07)
    assert fit["length_km"] == pytest.approx(truth["length_km"], rel=0.17)
    assert np.sqrt(fit["observation_variance"]) == pytest.approx(truth["sigma_o"], rel=0.16)


def test_scalar_update_is_the_textbook_weighting():
    # One state value observed directly: increment = b/(b+r) d, analysis variance = b r / (b + r).
    b, r, d = 0.16, 0.04, 0.5
    dx, reduction = analysis_update(np.array([d]), np.array([[r]]), np.array([[b]]), [np.array([[b]])])
    assert dx[0] == pytest.approx(b / (b + r) * d)
    assert b - reduction[0] == pytest.approx(b * r / (b + r))


def test_exact_and_useless_observations():
    rng = np.random.default_rng(2)
    n_state, n_obs = 12, 4
    root = rng.normal(size=(n_state, n_state))
    B = root @ root.T + n_state * np.eye(n_state)
    observed = np.array([1, 4, 7, 10])
    d = rng.normal(size=n_obs)
    # With R -> 0 the analysis matches the observations where they are.
    dx, _ = analysis_update(d, 1e-12 * np.eye(n_obs), B[np.ix_(observed, observed)], [B[:, observed]])
    np.testing.assert_allclose(dx[observed], d, atol=1e-8)
    # With R -> infinity the observations change nothing.
    dx, reduction = analysis_update(d, 1e12 * np.eye(n_obs), B[np.ix_(observed, observed)], [B[:, observed]])
    np.testing.assert_allclose(dx, 0.0, atol=1e-9)
    np.testing.assert_allclose(reduction, 0.0, atol=1e-9)


def test_blocked_update_matches_dense_formula():
    rng = np.random.default_rng(3)
    n_state, n_obs = 30, 6
    root = rng.normal(size=(n_state, n_state))
    B = root @ root.T
    H = rng.normal(size=(n_obs, n_state))
    root_r = rng.normal(size=(n_obs, n_obs))
    R = root_r @ root_r.T + np.eye(n_obs)   # correlated observation errors
    d = rng.normal(size=n_obs)
    gain = B @ H.T @ np.linalg.inv(H @ B @ H.T + R)
    A = (np.eye(n_state) - gain @ H) @ B
    BHt = B @ H.T
    dx, reduction = analysis_update(d, R, H @ B @ H.T, [BHt[:11], BHt[11:25], BHt[25:]])
    np.testing.assert_allclose(dx, gain @ d, rtol=1e-10)
    np.testing.assert_allclose(np.diag(B) - reduction, np.diag(A), rtol=1e-10, atol=1e-10)
    assert np.all(np.diag(A) <= np.diag(B) + 1e-12)


def test_separable_covariance_is_a_valid_covariance():
    # The product of a Gaussian correlation in distance and a correlation between levels is positive
    # semi-definite (Schur product theorem), so the analysis S = H B H^T + R can always be factorised.
    from ocean_uncertainty.stats import level_correlation, separable_covariance
    rng = np.random.default_rng(4)
    n, n_levels = 60, 5
    lon, lat = rng.uniform(150, 170, n), rng.uniform(30, 40, n)
    level = rng.integers(0, n_levels, n)
    corr = level_correlation(rng.normal(size=(n_levels, 40)))
    np.testing.assert_allclose(np.diag(corr), 1.0)
    sigma = rng.uniform(0.1, 0.5, n)
    B = separable_covariance(sigma, sigma, distance_km(lon, lat, lon, lat), 200.0, corr[np.ix_(level, level)])
    np.testing.assert_allclose(np.diag(B), sigma**2)
    assert np.linalg.eigvalsh(B).min() > -1e-10


def test_level_correlation_of_single_precision_anomalies_is_positive_semidefinite():
    # Nearly identical levels (as in a winter mixed layer) stored in float32: the correlation must still have no
    # negative eigenvalues beyond double-precision rounding.
    from ocean_uncertainty.stats import level_correlation
    rng = np.random.default_rng(5)
    common = rng.normal(size=20000)
    levels = np.stack([common + 1e-3 * rng.normal(size=common.size) for _ in range(18)]).astype(np.float32)
    corr = level_correlation(levels)
    assert np.linalg.eigvalsh(corr).min() > -1e-12
