"""Validation figures for the monthly error climatology (C1-C3). Reads the products, writes PNGs.

    python scripts/plot_climatology.py --products DIR --out DIR

No cartopy: maps are pcolormesh in plain longitude/latitude. Land and empty cells are blank (NaN).
"""

import argparse
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

SEASONS = {"DJF": [12, 1, 2], "MAM": [3, 4, 5], "JJA": [6, 7, 8], "SON": [9, 10, 11]}
LAYER = "thetao_std_0_50m"
DEPTH_M = 9.8   # the level shown in single-depth comparisons (near the mixed layer)


def region_mean(field, lat_name="latitude", lon_name="longitude"):
    """Area-weighted (cos latitude) mean over the horizontal dimensions."""
    w = np.cos(np.deg2rad(field[lat_name]))
    return field.weighted(w).mean((lat_name, lon_name))


def nearest_depth(ds, target):
    return int(np.argmin(np.abs(ds["depth"].values - target)))


def season_mean(da, months):
    return da.sel(month=months).mean("month")


def _map(ax, field, title, vmin, vmax, cmap="viridis", x="longitude", y="latitude"):
    m = ax.pcolormesh(field[x], field[y], field, vmin=vmin, vmax=vmax, cmap=cmap, shading="auto")
    ax.set_title(title, fontsize=9)
    ax.set_aspect(1 / np.cos(np.deg2rad(40)))
    return m


def fig_seasonal_maps(c1, out):
    """C1 0-50 m layer sigma, mean of each season."""
    vmax = float(c1[LAYER].quantile(0.99))
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for ax, (name, months) in zip(axes.flat, SEASONS.items()):
        m = _map(ax, season_mean(c1[LAYER], months), f"{name}  0-50 m sigma", 0, vmax)
    fig.colorbar(m, ax=axes, label="K", shrink=0.8)
    fig.suptitle("C1 reanalysis-spread, 0-50 m layer temperature sigma by season")
    fig.savefig(out / "fig1_c1_seasonal_maps.png", dpi=130)


def fig_mean_min_max(c1, out):
    """Annual mean, and the month-wise minimum and maximum, of the 0-50 m layer sigma."""
    layer = c1[LAYER]
    panels = [("annual mean", layer.mean("month")), ("month-min", layer.min("month")),
              ("month-max", layer.max("month"))]
    vmax = float(panels[2][1].quantile(0.99))
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for ax, (name, field) in zip(axes, panels):
        m = _map(ax, field, f"0-50 m sigma: {name}", 0, vmax)
    fig.colorbar(m, ax=axes, label="K", shrink=0.8)
    fig.suptitle("C1 0-50 m sigma: mean and seasonal envelope (min/max over the 12 months)")
    fig.savefig(out / "fig2_c1_mean_min_max.png", dpi=130)


def fig_structure(c1, out):
    """Seasonal cycle, depth-month structure, and zonal/meridional variability of C1 sigma."""
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    # (a) seasonal cycle of the region-mean layer sigma
    rm = region_mean(c1[LAYER])
    ax[0, 0].plot(c1["month"], rm, "o-")
    ax[0, 0].set(xlabel="month", ylabel="K", title="(a) region-mean 0-50 m sigma by month", xticks=range(1, 13))
    ax[0, 0].grid(alpha=0.3)

    # (b) depth-month Hovmoller of the region-mean sigma
    hov = region_mean(c1["thetao_std"]).transpose("depth", "month")
    m = ax[0, 1].pcolormesh(c1["month"], c1["depth"], hov, shading="auto", cmap="viridis")
    ax[0, 1].invert_yaxis()
    ax[0, 1].set(xlabel="month", ylabel="depth (m)", title="(b) region-mean sigma(depth, month)")
    fig.colorbar(m, ax=ax[0, 1], label="K")

    # (c) sigma vs latitude, by season (longitude mean of the layer sigma)
    zonal = c1[LAYER].mean("longitude")
    for name, months in SEASONS.items():
        ax[1, 0].plot(c1["latitude"], zonal.sel(month=months).mean("month"), label=name)
    ax[1, 0].set(xlabel="latitude", ylabel="K", title="(c) 0-50 m sigma vs latitude (lon mean)")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)

    # (d) sigma vs longitude, by season (latitude mean of the layer sigma)
    merid = c1[LAYER].mean("latitude")
    for name, months in SEASONS.items():
        ax[1, 1].plot(c1["longitude"], merid.sel(month=months).mean("month"), label=name)
    ax[1, 1].set(xlabel="longitude", ylabel="K", title="(d) 0-50 m sigma vs longitude (lat mean)")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)

    fig.suptitle("C1 structure: season, depth, and horizontal variability")
    fig.savefig(out / "fig3_c1_structure.png", dpi=130)


def fig_ordering(c1, c2, c3, out):
    """C1 < C2 < C3 ordering, the correlation length, and the vertical profiles."""
    k2, k3 = nearest_depth(c2, DEPTH_M), nearest_depth(c3, DEPTH_M)
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    # (a) ordering by month at ~10 m
    ax[0, 0].plot(c2["month"], c2["c1_sigma"].isel(depth=k2), "o-", label="C1 sigma (background)")
    ax[0, 0].plot(c2["month"], c2["rms"].isel(depth=k2), "s-", label="C2 RMS (Argo - reanalysis)")
    ax[0, 0].plot(c3["month"], c3["rms"].isel(depth=k3), "^-", label="C3 RMS (Argo - WOA23)")
    ax[0, 0].set(xlabel="month", ylabel="K", title=f"(a) ordering at ~{DEPTH_M:.0f} m", xticks=range(1, 13))
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.3)

    # (b) C2 correlation length by month (perturbation scale)
    ax[0, 1].plot(c2["month"], c2["length_km"].isel(depth=k2), "o-")
    ax[0, 1].set(xlabel="month", ylabel="km", title=f"(b) C2 correlation length at ~{DEPTH_M:.0f} m",
                 xticks=range(1, 13))
    ax[0, 1].grid(alpha=0.3)

    # (c) C1 sigma / C2 RMS ratio by month
    ax[1, 0].plot(c2["month"], c2["c1_over_rms"].isel(depth=k2), "o-")
    ax[1, 0].axhline(1.0, color="k", lw=0.8, ls="--")
    ax[1, 0].set(xlabel="month", ylabel="ratio", title="(c) C1 sigma / C2 RMS by month", xticks=range(1, 13))
    ax[1, 0].grid(alpha=0.3)

    # (d) annual-mean vertical profiles
    ax[1, 1].plot(region_mean(c1["thetao_std"]).mean("month"), c1["depth"], "o-", label="C1 sigma")
    ax[1, 1].plot(c2["rms"].mean("month"), c2["depth"], "s-", label="C2 RMS")
    ax[1, 1].plot(c3["rms"].mean("month"), c3["depth"], "^-", label="C3 RMS")
    ax[1, 1].invert_yaxis()
    ax[1, 1].set(xlabel="K", ylabel="depth (m)", title="(d) annual-mean vertical profile")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)

    fig.suptitle("Validation: C1 (spread) < C2 (vs reanalysis) < C3 (vs climatology)")
    fig.savefig(out / "fig4_ordering_validation.png", dpi=130)


def fig_map_comparison(c1, c2, c3, out):
    """Winter (DJF) error maps at ~10 m: C1 on its grid, C2 and C3 on the 2 degree Argo boxes."""
    k1 = nearest_depth(c1, DEPTH_M)
    c1_map = season_mean(c1["thetao_std"].isel(depth=k1), SEASONS["DJF"])
    c2_map = season_mean(c2["rms_map"].isel(depth=nearest_depth(c2, DEPTH_M)), SEASONS["DJF"])
    c3_map = season_mean(c3["rms_map"].isel(depth=nearest_depth(c3, DEPTH_M)), SEASONS["DJF"])
    vmax = float(np.nanmax([c2_map.quantile(0.95), c1_map.quantile(0.95)]))
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
    _map(ax[0], c1_map, "C1 sigma (0.25 grid)", 0, vmax)
    _map(ax[1], c2_map, "C2 RMS (2 deg boxes)", 0, vmax, x="map_longitude", y="map_latitude")
    m = _map(ax[2], c3_map, "C3 RMS (2 deg boxes)", 0, vmax, x="map_longitude", y="map_latitude")
    fig.colorbar(m, ax=ax, label="K", shrink=0.8)
    fig.suptitle(f"Winter (DJF) temperature error at ~{DEPTH_M:.0f} m: C1 vs C2 vs C3")
    fig.savefig(out / "fig5_map_comparison_djf.png", dpi=130)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--products", type=pathlib.Path, required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    c1 = xr.open_dataset(a.products / "c1.nc")
    c2 = xr.open_dataset(a.products / "c2.nc")
    c3 = xr.open_dataset(a.products / "c3.nc")
    fig_seasonal_maps(c1, a.out)
    fig_mean_min_max(c1, a.out)
    fig_structure(c1, a.out)
    fig_ordering(c1, c2, c3, a.out)
    fig_map_comparison(c1, c2, c3, a.out)
    print("wrote figures to", a.out)


if __name__ == "__main__":
    main()
