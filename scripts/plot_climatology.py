"""Validation figures for the monthly error climatology (C1, C2). Reads the products, writes PNGs.

    python scripts/plot_climatology.py --products DIR --out DIR

No cartopy: maps are pcolormesh in plain longitude/latitude. Land and empty cells are blank (NaN). The C2 figures
test whether member - mean is a realistic perturbation: for a reliable 3-member ensemble RMS(Argo - member) should
sit near sqrt(2) * spread and RMS(Argo - mean) near sqrt(4/3) * spread.
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
ANOM = "thetao_member_anomaly_0_50m"
DEPTH_M = 9.8   # the level shown in single-depth comparisons (near the mixed layer)
SQRT2, MEAN_FACTOR = np.sqrt(2.0), np.sqrt(4.0 / 3.0)   # reliable-ensemble member/mean vs spread (N=3)


def region_mean(field, lat_name="latitude", lon_name="longitude"):
    w = np.cos(np.deg2rad(field[lat_name]))
    return field.weighted(w).mean((lat_name, lon_name))


def nearest_depth(ds, target):
    return int(np.argmin(np.abs(ds["depth"].values - target)))


def season_mean(da, months):
    return da.sel(month=months).mean("month")


def _map(ax, field, title, vmin, vmax, cmap="viridis", x="longitude", y="latitude"):
    m = ax.contourf(field[x], field[y], field, levels=np.linspace(vmin, vmax, 21), extend="max", cmap=cmap)
    ax.set_title(title, fontsize=9)
    ax.set_aspect(1 / np.cos(np.deg2rad(40)))
    return m


# --- C1: the spread product -------------------------------------------------------------------------------------

def fig_seasonal_maps(c1, out):
    vmax = float(c1[LAYER].quantile(0.99))
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for ax, (name, months) in zip(axes.flat, SEASONS.items()):
        m = _map(ax, season_mean(c1[LAYER], months), f"{name}  0-50 m sigma", 0, vmax)
    fig.colorbar(m, ax=axes, label="K", shrink=0.8)
    fig.suptitle("C1 reanalysis-spread, 0-50 m layer temperature sigma by season")
    fig.savefig(out / "fig1_c1_seasonal_maps.png", dpi=130)


def fig_mean_min_max(c1, out):
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
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    rm = region_mean(c1[LAYER])
    ax[0, 0].plot(c1["month"], rm, "o-")
    ax[0, 0].set(xlabel="month", ylabel="K", title="(a) region-mean 0-50 m sigma by month", xticks=range(1, 13))
    ax[0, 0].grid(alpha=0.3)
    hov = region_mean(c1["thetao_std"]).transpose("depth", "month")
    m = ax[0, 1].contourf(c1["month"], c1["depth"], hov, levels=21, cmap="viridis")
    ax[0, 1].invert_yaxis()
    ax[0, 1].set(xlabel="month", ylabel="depth (m)", title="(b) region-mean sigma(depth, month)")
    fig.colorbar(m, ax=ax[0, 1], label="K")
    zonal = c1[LAYER].mean("longitude")
    for name, months in SEASONS.items():
        ax[1, 0].plot(c1["latitude"], zonal.sel(month=months).mean("month"), label=name)
    ax[1, 0].set(xlabel="latitude", ylabel="K", title="(c) 0-50 m sigma vs latitude (lon mean)")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.3)
    merid = c1[LAYER].mean("latitude")
    for name, months in SEASONS.items():
        ax[1, 1].plot(c1["longitude"], merid.sel(month=months).mean("month"), label=name)
    ax[1, 1].set(xlabel="longitude", ylabel="K", title="(d) 0-50 m sigma vs longitude (lat mean)")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.3)
    fig.suptitle("C1 structure: season, depth, and horizontal variability")
    fig.savefig(out / "fig3_c1_structure.png", dpi=130)


# --- C2: is member - mean a realistic perturbation? -------------------------------------------------------------

def fig_triad(c2, members, out):
    """Spread vs RMS(Argo - member) vs RMS(Argo - mean) by month at ~10 m, against the reliable-ensemble factors."""
    k = nearest_depth(c2, DEPTH_M)
    s = c2["c1_sigma"].isel(depth=k)                                   # spread at profiles (month)
    rms_mean = c2["rms"].sel(background="mean").isel(depth=k)
    rms_mem = c2["rms"].sel(background=members).isel(depth=k)          # (background, month)
    mo = c2["month"]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), constrained_layout=True)

    ax[0].plot(mo, s, "o-", color="k", label="spread = std(member - mean)")
    ax[0].plot(mo, rms_mean, "s-", label="RMS(Argo - mean)")
    ax[0].plot(mo, rms_mem.mean("background"), "^-", label="RMS(Argo - member), mean of 3")
    ax[0].fill_between(mo, rms_mem.min("background"), rms_mem.max("background"), alpha=0.2, label="member range")
    if "rg" in c2["background"].values:   # obs-only reference (1 deg, monthly): an upper bound, not a tight nugget
        ax[0].plot(mo, c2["rms"].sel(background="rg").isel(depth=k), "D-", color="tab:red", ms=4,
                   label="RMS(Argo - RG), obs-only 1 deg")
    ax[0].plot(mo, SQRT2 * s, "--", color="gray", label=r"$\sqrt{2}$ * spread (reliable member)")
    ax[0].plot(mo, MEAN_FACTOR * s, ":", color="gray", label=f"$\\sqrt{{4/3}}$ * spread (reliable mean)")
    ax[0].set(xlabel="month", ylabel="K", title=f"(a) amplitudes at ~{DEPTH_M:.0f} m", xticks=range(1, 13))
    ax[0].legend(fontsize=7.5); ax[0].grid(alpha=0.3)

    ax[1].plot(mo, rms_mem.mean("background") / s, "^-", label="RMS(Argo - member) / spread")
    ax[1].plot(mo, rms_mean / s, "s-", label="RMS(Argo - mean) / spread")
    ax[1].axhline(SQRT2, ls="--", color="gray"); ax[1].axhline(MEAN_FACTOR, ls=":", color="gray")
    ax[1].text(12.2, SQRT2, r"$\sqrt{2}$", va="center", fontsize=8)
    ax[1].text(12.2, MEAN_FACTOR, "$\\sqrt{4/3}$", va="center", fontsize=8)
    ax[1].set(xlabel="month", ylabel="ratio", title="(b) ratio to spread (above line = under-dispersed + repr. error)",
              xticks=range(1, 13))
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.suptitle("Is (member - mean) a realistic perturbation? Argo-based spread-skill check")
    fig.savefig(out / "fig4_spread_skill_triad.png", dpi=130)


def fig_member_anomaly_maps(c1, out):
    """Winter (DJF) 0-50 m RMS of (member - mean): the perturbation realisation each member provides."""
    anom = season_mean(c1[ANOM], SEASONS["DJF"])            # (member, lat, lon)
    members = list(c1["member"].values)
    vmax = float(anom.quantile(0.99))
    fig, axes = plt.subplots(1, len(members), figsize=(5 * len(members), 4.2), constrained_layout=True)
    for ax, mem in zip(np.atleast_1d(axes), members):
        m = _map(ax, anom.sel(member=mem), f"{mem}: RMS(member - mean)", 0, vmax)
    fig.colorbar(m, ax=axes, label="K", shrink=0.8)
    fig.suptitle("C1 perturbation field: winter (DJF) 0-50 m amplitude of member - mean, per member")
    fig.savefig(out / "fig5_member_minus_mean_maps.png", dpi=130)


def fig_argo_maps(c2, members, out):
    """Winter (DJF) RMS maps at ~10 m of Argo - mean, Argo - each member, and (if present) Argo - RG (2 deg boxes)."""
    k = nearest_depth(c2, DEPTH_M)
    backgrounds = ["mean"] + members + (["rg"] if "rg" in c2["background"].values else [])
    maps = {b: season_mean(c2["rms_map"].sel(background=b).isel(depth=k), SEASONS["DJF"]) for b in backgrounds}
    vmax = float(np.nanmax([v.quantile(0.95) for v in maps.values()]))
    fig, axes = plt.subplots(1, len(backgrounds), figsize=(4.6 * len(backgrounds), 4.2), constrained_layout=True)
    for ax, b in zip(axes, backgrounds):
        m = _map(ax, maps[b], f"Argo - {b}", 0, vmax, x="map_longitude", y="map_latitude")
    fig.colorbar(m, ax=axes, label="K", shrink=0.8)
    fig.suptitle(f"Winter (DJF) RMS(Argo - background) at ~{DEPTH_M:.0f} m (2 deg boxes)")
    fig.savefig(out / "fig6_argo_minus_background_maps.png", dpi=130)


def fig_c3(c3, out):
    """C3: where and how much a month of Argo adjusts the upper ocean."""
    k = nearest_depth(c3, DEPTH_M)
    red = 100 * season_mean(c3["uncertainty_reduction"].isel(depth=k), SEASONS["DJF"])
    rms = season_mean(c3["increment_rms"].isel(depth=k), SEASONS["DJF"])
    fig, ax = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    m0 = _map(ax[0, 0], red, f"(a) DJF uncertainty reduction at ~{DEPTH_M:.0f} m", 0, float(red.quantile(0.99)))
    fig.colorbar(m0, ax=ax[0, 0], label="%")
    m1 = _map(ax[0, 1], rms, f"(b) DJF increment RMS at ~{DEPTH_M:.0f} m", 0, float(rms.quantile(0.99)))
    fig.colorbar(m1, ax=ax[0, 1], label="K")
    ax[1, 0].plot(c3["month"], 100 * region_mean(c3["uncertainty_reduction"].isel(depth=k)), "o-")
    ax[1, 0].set(xlabel="month", ylabel="%", title=f"(c) region-mean reduction at ~{DEPTH_M:.0f} m",
                 xticks=range(1, 13)); ax[1, 0].grid(alpha=0.3)
    ax[1, 1].plot(100 * region_mean(c3["uncertainty_reduction"]).mean("month"), c3["depth"], "o-")
    ax[1, 1].invert_yaxis()
    ax[1, 1].set(xlabel="%", ylabel="depth (m)", title="(d) annual region-mean reduction vs depth")
    ax[1, 1].grid(alpha=0.3)
    fig.suptitle("C3 Argo analysis increment: where and how much Argo adjusts the upper ocean")
    fig.savefig(out / "fig7_c3_increment.png", dpi=130)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--products", type=pathlib.Path, required=True)
    p.add_argument("--out", type=pathlib.Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    c1 = xr.open_dataset(a.products / "c1.nc")
    c2 = xr.open_dataset(a.products / "c2.nc")
    members = [b for b in c2["background"].values.tolist() if b not in ("mean", "rg")]
    fig_seasonal_maps(c1, a.out)
    fig_mean_min_max(c1, a.out)
    fig_structure(c1, a.out)
    fig_triad(c2, members, a.out)
    fig_member_anomaly_maps(c1, a.out)
    fig_argo_maps(c2, members, a.out)
    c3_path = a.products / "c3.nc"
    if c3_path.exists():
        fig_c3(xr.open_dataset(c3_path), a.out)
    print("wrote figures to", a.out)


if __name__ == "__main__":
    main()
