"""C3 — how much a month of Argo profiles adjusts the upper ocean: the analysis increment.

    delta_x = B H^T (H B H^T + R)^-1 d

computed as an independent 2-D optimal interpolation at each depth level (vertical correlation neglected). The
pieces come from the other products:

  sigma_b  background-error std   -> C1 (thetao_std) at the calendar month, bilinear at each point;
  L        horizontal length      -> C2 (length_km), median over levels unless given;
  sigma_o  observation error       -> C2 (observation_variance), per level unless given;
  d        innovation Argo - bg    -> the matched-profile files (innov_mean, or innov_<member>).

``build_c3`` runs it for every (year, month) and pools by calendar month into the typical adjustment and the
uncertainty it removes. Run ``python -m ocean_uncertainty.increment --help``.
"""

import argparse
import pathlib

import numpy as np
import xarray as xr

from ocean_uncertainty import sources
from ocean_uncertainty.products import write
from ocean_uncertainty.stats import analysis_update, distance_km


def _gaussian_cov(sig_a, lon_a, lat_a, sig_b, lon_b, lat_b, length_km):
    """Separable background covariance between two point sets: sigma_a sigma_b exp(-r^2 / 2L^2)."""
    return np.outer(sig_a, sig_b) * np.exp(-0.5 * (distance_km(lon_a, lat_a, lon_b, lat_b) / length_km) ** 2)


def increment(c1: xr.Dataset, c2: xr.Dataset, table, month: int, background: str = "mean",
              length_km: float | None = None, obs_error_std: float | None = None, superob_deg: float = 0.5,
              gross_error_sigma: float = 3.0, block: int = 4000) -> xr.Dataset:
    """2-D optimal-interpolation increment to the ensemble-mean background from the profiles in ``table``.

    Returns ``increment``, ``background_std`` and ``analysis_std`` on the C1 grid for one (year, month). Profiles
    within ``superob_deg`` are averaged first: nearby Argo is not independent (correlated representativeness), and
    treating it as independent makes H B H^T + R ill-conditioned and the increment overshoot.
    """
    col = f"innov_{background}"
    sigma = c1["thetao_std"].sel(month=month)                 # (depth, latitude, longitude)
    depths = sigma["depth"].values
    lon, lat = sigma["longitude"].values, sigma["latitude"].values
    LON, LAT = np.meshgrid(lon, lat)

    if length_km is None:
        # C2 length_km (the HL fit) is too noisy per month to drive the OI (NaN to >1000 km); use the documented
        # ~50 km winter correlation length as a robust prescribed value. Override with --length-km.
        length_km = 50.0
    obs_var = c2["observation_variance"].sel(background=background, month=month).values

    def sigma_o(k):
        if obs_error_std is not None:
            return obs_error_std
        v = obs_var[k]
        if not np.isfinite(v) or v <= 0:
            v = np.nanmedian(obs_var[obs_var > 0]) if np.isfinite(obs_var).any() and (obs_var > 0).any() else 0.09
        return max(np.sqrt(v), 0.05)

    inc = np.full(sigma.shape, np.nan)
    analysis = np.full(sigma.shape, np.nan)
    n_used = 0
    for k, depth in enumerate(depths):
        s2d = sigma.isel(depth=k)
        rows = table[np.isclose(table["depth"], depth) & np.isfinite(table[col])]
        if len(rows) < 5:
            continue
        # Gross-error check: drop frontal/representativeness outliers the smooth background cannot carry.
        dd = rows[col].to_numpy()
        rows = rows[np.abs(dd - np.median(dd)) <= gross_error_sigma * np.std(dd)]
        # Superob: average profiles sharing a superob_deg cell (they are not independent observations).
        rows = rows.assign(_lb=np.round(rows["lon"].to_numpy() / superob_deg),
                           _tb=np.round(rows["lat"].to_numpy() / superob_deg))
        g = rows.groupby(["_lb", "_tb"])
        lon_o, lat_o, d = g["lon"].mean().to_numpy(), g["lat"].mean().to_numpy(), g[col].mean().to_numpy()
        if lon_o.size < 3:
            continue
        sig_o = s2d.interp(longitude=xr.DataArray(lon_o, dims="o"), latitude=xr.DataArray(lat_o, dims="o")).values
        ok = np.isfinite(sig_o)
        if ok.sum() < 3:
            continue
        lon_o, lat_o, sig_o, d = lon_o[ok], lat_o[ok], sig_o[ok], d[ok]
        n_used += d.size

        cov_obs = _gaussian_cov(sig_o, lon_o, lat_o, sig_o, lon_o, lat_o, length_km)
        R = np.diag(np.full(d.size, sigma_o(k) ** 2))
        wet = np.isfinite(s2d.values)
        sig_s, lon_s, lat_s = s2d.values[wet], LON[wet], LAT[wet]

        def blocks():
            for i in range(0, sig_s.size, block):
                sl = slice(i, i + block)
                yield _gaussian_cov(sig_s[sl], lon_s[sl], lat_s[sl], sig_o, lon_o, lat_o, length_km)

        dx, reduction = analysis_update(d, R, cov_obs, blocks())
        field = np.full(s2d.shape, np.nan); field[wet] = dx
        inc[k] = field
        a = np.full(s2d.shape, np.nan); a[wet] = np.sqrt(np.maximum(sig_s ** 2 - reduction, 0.0))
        analysis[k] = a

    ds = xr.Dataset(
        {"increment": (("depth", "latitude", "longitude"), inc),
         "background_std": (("depth", "latitude", "longitude"), sigma.values),
         "analysis_std": (("depth", "latitude", "longitude"), analysis)},
        coords={"depth": depths, "latitude": lat, "longitude": lon})
    ds.attrs.update(month=int(month), background=background, length_km=float(length_km), n_obs=int(n_used))
    return ds


def _month_of(path: pathlib.Path) -> int:
    tag = next(t for t in pathlib.Path(path).stem.split("_") if t.isdigit() and len(t) >= 6)
    return int(tag[4:6])


def build_c3(innovation_files, c1: xr.Dataset, c2: xr.Dataset, out: pathlib.Path, background: str = "mean",
             length_km: float | None = None, obs_error_std: float | None = None) -> pathlib.Path:
    """Pool per-(year, month) increments by calendar month: the typical Argo adjustment and the uncertainty removed."""
    by_month: dict[int, list] = {}
    for path in innovation_files:
        with xr.open_dataset(path) as ds:
            table = ds.to_dataframe()
        inc = increment(c1, c2, table, _month_of(path), background, length_km, obs_error_std)
        by_month.setdefault(inc.attrs["month"], []).append(inc)

    per_month = []
    for month in sorted(by_month):
        stack = xr.concat(by_month[month], "year")
        per_month.append(xr.Dataset({
            "increment_mean": stack["increment"].mean("year"),                 # systematic adjustment
            "increment_rms": np.sqrt((stack["increment"] ** 2).mean("year")),  # typical magnitude
            "background_std": stack["background_std"].isel(year=0, drop=True),
            "analysis_std": stack["analysis_std"].mean("year"),
        }).expand_dims(month=[month]))
    ds = xr.concat(per_month, "month")
    ds["uncertainty_reduction"] = 1.0 - ds["analysis_std"] / ds["background_std"]
    return write(ds, out, "C3 climatological Argo analysis increment",
                 [sources.GREP["citation"], sources.ARGO["citation"]],
                 "Per (year, month): delta_x = B H^T (H B H^T + R)^-1 d, a 2-D optimal interpolation at each depth "
                 "level, with sigma_b from C1, horizontal length and observation error from C2, and innovation "
                 f"Argo - {background}. Pooled by calendar month: increment_mean (systematic), increment_rms "
                 "(typical magnitude), and analysis_std vs background_std (the uncertainty Argo removes).",
                 "Vertical correlation neglected (levels analysed independently); a single horizontal length; "
                 "climatological sigma_b applied to each month's own profiles. A realistic magnitude/reach of the "
                 "Argo adjustment, not a specific operational analysis.")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", type=pathlib.Path, required=True, help="dir of innovations_YYYYMM.nc files")
    p.add_argument("--c1", type=pathlib.Path, required=True)
    p.add_argument("--c2", type=pathlib.Path, required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    p.add_argument("--background", default="mean", help="mean (default) or a member: glor/cglo/oras")
    p.add_argument("--length-km", type=float, help="override the C2 correlation length")
    p.add_argument("--obs-error-std", type=float, help="override the C2 observation error (K)")
    a = p.parse_args()
    with xr.open_dataset(a.c1) as c1, xr.open_dataset(a.c2) as c2:
        out = build_c3(sorted(a.data_dir.glob("innovations_*.nc")), c1.load(), c2.load(), a.out, a.background,
                       a.length_km, a.obs_error_std)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
