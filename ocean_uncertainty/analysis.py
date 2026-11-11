"""Error covariances for an analysis of Argo profiles against the reanalysis-ensemble mean (used by U6 and checks).

An analysis needs the background-error covariance between any two sets of points (B, from the reanalysis
ensemble U1 and the Argo statistics U2) and the observation-error covariance R of the profiles (from U5, or
scaled to U2's estimate). Both are built here so the U6 product and the hold-out check use the same models.
"""

import dataclasses
import pathlib

import numpy as np
import pandas as pd
import xarray as xr

from ocean_uncertainty.observations import sample
from ocean_uncertainty.stats import (distance_km, ensemble_covariance, gaspari_cohn, gaspari_cohn_half_width,
                                     level_correlation, separable_covariance)

BACKGROUND_MODELS = ("separable", "ensemble")
OBSERVATION_MODELS = ("u5", "innovations", "u5_scaled")


@dataclasses.dataclass(frozen=True)
class Points:
    """Positions, level numbers and ensemble anomalies (n, members) of a set of points."""
    lon: np.ndarray
    lat: np.ndarray
    level: np.ndarray
    anomalies: np.ndarray

    @property
    def sigma(self) -> np.ndarray:
        return np.sqrt(np.sum(self.anomalies**2, axis=1) / (self.anomalies.shape[1] - 1))

    def __getitem__(self, index) -> "Points":
        return Points(self.lon[index], self.lat[index], self.level[index], self.anomalies[index])


@dataclasses.dataclass(frozen=True)
class ErrorModel:
    """The background-error model (B) and the correlation structure of the observation error (R)."""
    background_model: str
    length_km: float
    level_corr: np.ndarray
    obs_level_corr: np.ndarray
    instrument_variance: float

    def background(self, a: Points, b: Points) -> np.ndarray:
        """B between the points of ``a`` and ``b``.

        - "separable": sigma_a sigma_b exp(-r^2 / 2L^2) C(level_a, level_b), with sigma and C from the ensemble;
        - "ensemble": the ensemble covariance times a Gaspari-Cohn localisation matched to the Gaussian of length L.
        """
        dist = distance_km(a.lon, a.lat, b.lon, b.lat)
        if self.background_model == "ensemble":
            return ensemble_covariance(a.anomalies, b.anomalies,
                                       gaspari_cohn(dist, gaspari_cohn_half_width(self.length_km)))
        return separable_covariance(a.sigma, b.sigma, dist, self.length_km,
                                    self.level_corr[np.ix_(a.level, b.level)])

    def observation(self, profile: np.ndarray, level: np.ndarray, variance: np.ndarray) -> np.ndarray:
        """R: correlated between the levels of one profile, independent between profiles, plus instrument error."""
        same_profile = profile[:, None] == profile[None, :]
        return (np.sqrt(np.outer(variance, variance)) * self.obs_level_corr[np.ix_(level, level)] * same_profile
                + self.instrument_variance * np.eye(profile.size))


@dataclasses.dataclass
class Window:
    """Everything an analysis on one date needs: the observations, the state and the error model."""
    obs: pd.DataFrame
    obs_points: Points
    obs_variance: np.ndarray
    state: xr.DataArray
    state_points: Points
    wet: np.ndarray
    errors: ErrorModel


def load_window(u2_file: pathlib.Path, u1_file: pathlib.Path, u5_file: pathlib.Path, date: str,
                window_days: float, variable: str, background_model: str, obs_model: str,
                length_km: float | None = None) -> Window:
    """Innovations within ``window_days`` centred on ``date``, and the error models for them.

    ``obs_model`` sets the observation-error variance:
    - "u5": U5's representativeness variance at each observation;
    - "innovations": U2's fitted observation-error variance at each level (Hollingsworth-Lonnberg split), the same
      everywhere;
    - "u5_scaled": U5's spatial pattern, scaled at each level so its mean over all U2 observations equals U2's
      fitted observation-error variance (pattern from the high-resolution reanalysis, size from the innovations).
    All use U5's correlation between levels. ``length_km`` defaults to the median of U2's fitted lengths.
    """
    if background_model not in BACKGROUND_MODELS or obs_model not in OBSERVATION_MODELS:
        raise ValueError(f"background_model in {BACKGROUND_MODELS}, obs_model in {OBSERVATION_MODELS}")
    when = pd.Timestamp(date)
    with xr.open_dataset(u2_file) as u2:
        if u2.attrs.get("background") != "grep_mean":
            raise ValueError("needs U2 innovations against the GREP ensemble mean")
        columns = ["profile", "time", "lon", "lat", "depth", "innovation"]
        table = u2[[f"obs_{c}" for c in columns]].to_dataframe()
        table.columns = columns
        fitted = u2[["length_km", "observation_variance"]].load()
    all_obs = table[np.isfinite(table["innovation"])]
    obs = table[((table["time"] - when).abs() <= pd.Timedelta(days=window_days / 2))
                & np.isfinite(table["innovation"])]
    if length_km is None:
        if not np.any(np.isfinite(fitted["length_km"].values)):
            raise ValueError("U2 has no fitted correlation length (too few pairs); pass length_km")
        length_km = float(np.nanmedian(fitted["length_km"].values))

    with xr.open_dataset(u1_file) as u1:
        anomaly = u1[f"{variable}_anomaly"].sel(time=when.floor("D")).load()
    with xr.open_dataset(u5_file) as u5:
        representativeness = u5["representativeness_variance"].load()
        obs_level_corr = u5["representativeness_level_correlation"].values
        instrument = u5.attrs["instrument_variance"]
    for d in ("depth", "latitude", "longitude"):
        if not np.allclose(representativeness[d].values, anomaly[d].values):
            raise ValueError(f"U5 must be built on the U1 grid ({d} differs)")

    depth = anomaly["depth"].values
    args = (obs["lon"].values, obs["lat"].values, obs["depth"].values, np.zeros(len(obs), int))
    anomalies = np.stack([sample(anomaly.sel(member=m).expand_dims(time=1), *args)
                          for m in anomaly["member"].values], axis=1)
    fitted_variance = fitted["observation_variance"].sel(depth=obs["depth"].values, method="nearest").values
    if obs_model == "innovations":
        variance = fitted_variance
    else:
        variance = sample(representativeness.expand_dims(time=1), *args)
    if obs_model == "u5_scaled":
        everywhere = sample(representativeness.expand_dims(time=1), all_obs["lon"].values, all_obs["lat"].values,
                            all_obs["depth"].values, np.zeros(len(all_obs), int))
        mean_u5 = pd.Series(everywhere).groupby(all_obs["depth"].values).mean()
        variance = variance * fitted_variance / mean_u5.reindex(obs["depth"].values).values
    usable = np.all(np.isfinite(anomalies), axis=1) & np.isfinite(variance) & (variance > 0)
    obs, anomalies, variance = obs[usable], anomalies[usable], variance[usable]
    if obs.empty:
        raise ValueError("no innovation in the window lies inside the U1 and U5 grids")
    obs_points = Points(obs["lon"].values, obs["lat"].values, np.searchsorted(depth, obs["depth"].values),
                        anomalies)

    state = anomaly.transpose("depth", "latitude", "longitude", "member")
    values = state.values.reshape(-1, state.sizes["member"])
    wet = np.all(np.isfinite(values), axis=1)
    shape = state.shape[:3]
    lon, lat = np.meshgrid(state["longitude"].values, state["latitude"].values)
    level = np.broadcast_to(np.arange(shape[0])[:, None, None], shape).reshape(-1)
    state_points = Points(np.broadcast_to(lon, shape).reshape(-1)[wet], np.broadcast_to(lat, shape).reshape(-1)[wet],
                          level[wet], values[wet])
    columns = state.values.reshape(shape[0], -1)
    level_corr = level_correlation(columns[:, np.all(np.isfinite(columns), axis=0)])

    errors = ErrorModel(background_model, length_km, level_corr, obs_level_corr, instrument)
    return Window(obs, obs_points, variance, state, state_points, wet, errors)
