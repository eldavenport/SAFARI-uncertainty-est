import numpy as np
import pandas as pd

from ocean_uncertainty.analysis import ErrorModel, Points, Window
from ocean_uncertainty.checks import hold_out, hold_out_summary


def synthetic_window(seed: int) -> Window:
    """Innovations drawn from exactly the B + R the analysis assumes (80 profiles of 4 levels)."""
    rng = np.random.default_rng(seed)
    n_profiles, n_levels, n_members = 80, 4, 50
    lon = np.repeat(rng.uniform(160, 175, n_profiles), n_levels)
    lat = np.repeat(rng.uniform(30, 40, n_profiles), n_levels)
    level = np.tile(np.arange(n_levels), n_profiles)
    corr = np.array([[1, .9, .8, .6], [.9, 1, .9, .8], [.8, .9, 1, .9], [.6, .8, .9, 1.]])
    points = Points(lon, lat, level, 0.3 * rng.normal(size=(lon.size, n_members)))
    errors = ErrorModel("separable", 150.0, corr, corr, 1e-6)
    profile = np.repeat(np.arange(n_profiles), n_levels)
    variance = np.full(lon.size, 0.04)
    cov = errors.background(points, points) + errors.observation(profile, level, variance)
    d = np.linalg.cholesky(cov) @ rng.normal(size=lon.size)
    obs = pd.DataFrame(dict(profile=profile, depth=level.astype(float), innovation=d))
    return Window(obs, points, variance, None, None, None, errors)


def test_hold_out_is_calibrated_when_the_error_model_is_right():
    # If d really has covariance B + R, the held-out residuals have exactly the variance the analysis expects.
    # Over 30 seeds (run 2026-11-10) one seed's calibration had mean 1.006 and standard deviation 0.17, so the mean
    # of 10 seeds has a standard deviation of about 0.054; the tolerance is 3 of those.
    calibration = [hold_out_summary(hold_out(synthetic_window(seed), 5), [0, 1e9], 1).loc[0, "calibration"]
                   for seed in range(10)]
    assert abs(np.mean(calibration) - 1.0) < 0.16


def test_skill_is_positive_for_correlated_profiles():
    # Background errors correlated over 150 km between profiles a few hundred km apart: neighbours carry
    # information, so predicting held-out innovations must beat predicting zero.
    skill = [hold_out_summary(hold_out(synthetic_window(seed), 5), [0, 1e9], 1).loc[0, "skill"] for seed in range(10)]
    assert np.mean(skill) > 0
