# Product catalogue: upper-ocean temperature uncertainty

This repository builds estimates of how well the upper-ocean temperature is known in a typical operational ocean
data assimilation (DA) system, and of how much real observations adjust it. It produces these products and
documents them. It does not apply them to any model.

**Status of facts.** Every value below cites its source. Some values are marked **UNVERIFIED** because no primary
source confirming them was found (research of 2026-11-09). Code must not depend on an unverified value until it is
resolved. Dataset IDs and variable names marked *(catalogue)* were read from the live Copernicus Marine catalogue
with `copernicusmarine.describe` (toolbox 2.3.0) on 2026-11-09.

---

## Notation

| Symbol | Meaning |
|---|---|
| $x_b$ | background state: the best estimate before new observations are used |
| $y$ | observations |
| $H$ | observation operator: maps a gridded state to the observation locations and depths |
| $d = y - Hx_b$ | innovation (observation minus background) |
| $B$, $R$ | background-error and observation-error covariances |
| $A$ | analysis-error covariance, $A = B - BH^\top(HBH^\top+R)^{-1}HB$ |
| $\delta x = BH^\top(HBH^\top+R)^{-1}d$ | analysis increment: the change the observations make to the background |

The update for $\delta x$ and $A$ is the best linear unbiased estimate, i.e. optimal interpolation (Kalnay 2003,
ch. 5; equation number **UNVERIFIED**).

---

## U1 — Spread of an ensemble of ocean reanalyses

**What it represents.** How much independent, state-of-the-art ocean reanalyses disagree on a given day. This
estimates background/analysis uncertainty in a typical DA system.

| Item | Value | Source |
|---|---|---|
| Product | Global Ocean Ensemble Physics Reanalysis (GREP), `GLOBAL_MULTIYEAR_PHY_ENS_001_031` | product page |
| DOI | 10.48670/moi-00024 | product page |
| Members (3) | GLORYS2V4 (Mercator; NEMO 3.1 with a SEEK filter, 7-day window); C-GLORS (CMCC; NEMO 3.4 with 3D-Var, 7-day window); ORAS5 (ECMWF; NEMO 3.4.1 with NEMOVAR 3D-Var, 5-day window) | QUID 001-031 Issue 1.3, Table 1 |
| FOAM | removed; the whole series was recomputed with three members | QUID 001-031 §V |
| Datasets | `cmems_mod_glo_phy-all_my_0.25deg_P1D-m` (members, daily); `cmems_mod_glo_phy-mnstd_my_0.25deg_P1D-m` (mean and std, daily); `_P1M-m` monthly versions | *(catalogue)* |
| Version | 202311 | *(catalogue)* |
| Period | 1993-01-01 to 2025-12-31, daily | *(catalogue)* |
| Grid | regular 0.25°, lon −180 to 179.75, lat −80 to 90; 75 levels | *(catalogue)*; PUM 001-031 |
| Levels in 0–100 m (m) | 0.50, 1.55, 2.66, 3.85, 5.14, 6.54, 8.09, 9.82, 11.77, 13.99, 16.52, 19.42, 22.75, 26.55, 30.87, 35.74, 41.18, 47.21, 53.85, 61.11, 69.02, 77.61, 86.92, 97.04 | PUM 001-031 p.9 |
| Variables | `thetao_glor`, `thetao_oras`, `thetao_cglo`; `thetao_mean`, `thetao_std`; the same pattern for `so`, `mlotst` | *(catalogue)* |
| Observations used | all members: SST, sea level anomaly, sea-ice concentration, in situ T/S profiles. Profile datasets: CORA (GLORYS2V4), EN4 (C-GLORS, ORAS5) | QUID 001-031 Table 1 |
| Atmospheric forcing | ERA-Interim, with ERA5 after 2019 | QUID 001-031 p.8 |

**Delivered form.** For each date: the member anomalies $e_k = x_k - \bar x$ ($k = 1, 2, 3$) and $\sigma$. With
$N = 3$, the ensemble covariance $\hat B = \frac{1}{N-1}\sum_k e_k e_k^\top$ has rank 2. A user can form
$g^\top \hat B g = \frac{1}{N-1}\sum_k \langle g, e_k\rangle^2$ without ever building $\hat B$.

**Known limitations.**
- **Few members.** For Gaussian errors, a 3-member variance has 2 degrees of freedom, so a single day's $\sigma$ is
  uncertain by a factor of order 2. Averaging over time or space, or using the monthly product, reduces this.
  Increments (U6) need B localised in space because of the low rank.
- **Shared inputs.** The members use the same forcing family, the same observation types and the same ocean model
  family (NEMO, ORCA025), and they share biases. The QUID notes "most temperature anomaly patterns are similar
  between individual members", including about +1 °C north of the Gulf Stream and Kuroshio (QUID 001-031 p.12–14).
  The spread therefore probably **underestimates** the true error. The QUID itself does not assess dispersion, so
  this is checked here (see Verification).
- **Reported size.** SST std is mostly 0.1–0.3 °C. It reaches 0.6 °C, and locally >1 °C, in western boundary
  currents (including the Kuroshio). The std peaks near 100 m and drops after 2004 with Argo
  (QUID 001-031 p.6, p.10).
- `thetao_std`'s normalisation (N or N−1) is **UNVERIFIED**. The members are therefore downloaded, and σ is
  recomputed with N−1.

**References.**
- Balmaseda, M. A., et al. (2015), J. Oper. Oceanogr., 8(sup1), S80–S97, doi:10.1080/1755876X.2015.1022329
  (the Ocean Reanalyses Intercomparison Project, ORA-IP).
- Palmer, M. D., et al. (2017), Clim. Dyn., 49, 909–930, doi:10.1007/s00382-015-2801-0 (heat-content spread across
  reanalyses).

---

## U2 — Argo profiles minus a reanalysis (innovation statistics)

**What it represents.** The difference between independent point observations and a gridded analysis. This is an
observation-based estimate of analysis error plus observation error.

**It is a lower bound.** Every reanalysis used here assimilates Argo, so the analysis is pulled toward these same
profiles.

| Item | Value | Source |
|---|---|---|
| Observations | Argo core profiles, temperature | Argo (2000), SEANOE, doi:10.17882/42182 (cite with the snapshot key) |
| Access | `argopy` `DataFetcher().region([lon_min, lon_max, lat_min, lat_max, p_min, p_max, t_start, t_end])`; the upper time bound is exclusive | argopy docs; Maze & Balem (2020), JOSS 5(53), 2425, doi:10.21105/joss.02425 |
| Quality selection | argopy "research" mode: delayed-mode only, QC = 1, pressure error < 20 dbar. "standard" mode (QC 1–2, including real-time) is the alternative | argopy user-mode docs |
| Field rule | data mode R uses `PARAM`; A or D uses `PARAM_ADJUSTED` | argo.ucsd.edu "how to use Argo files" |
| Instrument accuracy | 0.002 °C (T), 2.4 dbar (P), after delayed-mode adjustment | Wong et al. (2020), Front. Mar. Sci. 7:700, doi:10.3389/fmars.2020.00700 |
| Backgrounds | GLORYS12V1 `cmems_mod_glo_phy_my_0.083deg_P1D-m` (1/12°, daily, 1993-01-01 to 2026-09-22 *(catalogue)*, 7-day window with incremental analysis update, IAU; Lellouche et al. 2021, Front. Earth Sci. 9:698876, doi:10.3389/feart.2021.698876); and the GREP mean (U1) | |

**Delivered form.**
- A matched table: profile ID, time, location, depth, observed T, background T, difference.
- Binned σ(lon, lat, depth, season).
- Spatial and vertical correlation lengths fitted to the differences of profile pairs.

**Published comparison values.** The GLORYS12V1 QUID (001-030 Issue 1.6) shows Argo innovation statistics only as
figures (Figs 1.7–1.8). It gives no tables. The numbers stated in its text are:
- a surface warm bias of about 0.1 °C;
- an RMS difference below 1 °C in the upper 300 m in the Argo era.

Lellouche et al. (2021) give analysis-minus-observation RMS of about 0.45 °C (Argo era, global, full depth). These are
the only published values available for comparison. Region- and layer-specific published values are **UNVERIFIED**.

---

## U3 — Argo profiles minus a climatology

**What it represents.** The error of using a climatology as the ocean state, i.e. no information about the day.
This is an **upper bound** for DA-scale error.

| Climatology | Fit for purpose | Source |
|---|---|---|
| **WOA23 (recommended)** | CF netCDF; monthly, 0.25° (for A5B4 = 2005–2014, B5C2 = 2015–2022, the climate normals and "decav") or 1°; 57 levels over 0–1500 m, 5 m spacing in 0–100 m. Also provides `sd` (standard deviation about the statistical mean), `se`, `sdo` (objectively analysed standard deviation) and `sea` (standard error of the analysis) | Locarnini et al. (2024), NOAA Atlas NESDIS 89, doi:10.25923/54bh-1613; WOA23 Product Documentation Tables 1, 2, 4 |
| RG09 (Roemmich–Gilson) | Argo-only climatology, 2004–2018 baseline ("2019" version). Not CF-compliant (gzipped netCDF). Grid and levels **UNVERIFIED** | Roemmich & Gilson (2009), Prog. Oceanogr. 82, 81–100 (DOI **UNVERIFIED**); sio-argo.ucsd.edu/RG_Climatology.html |
| ISAS20 | 0.5°, 2002–2020; whether it has an error field is **UNVERIFIED** | Gaillard et al. (2016), J. Climate 29, 1305–1323, doi:10.1175/JCLI-D-15-0028.1; SEANOE doi:10.17882/52367 |

**Recommended choice:** WOA23 monthly 0.25°, climate normal 1991–2020 (`decav91C0`). It spans the period of
interest without being tuned to one decade, and its `sdo` field gives an independent, map-based σ_clim to compare
with the Argo-based one. The choice of normal is an option of the builder, not fixed in code.

---

## U4 — SST analysis uncertainty

**What it represents.** The error standard deviation that each SST analysis assigns to itself. It is a surface
value only.

| Product | Variable, units | Period, grid | How it is derived | Calibration evidence | Source |
|---|---|---|---|---|---|
| **OSTIA reprocessed** `METOFFICE-GLO-SST-L4-REP-OBS-SST` *(catalogue)* | `analysis_error`, K | 1981-10-01 to 2026-06-30 *(catalogue)*, 0.05° daily | From observation weights combined with prescribed background-error variances. 3 days of data per analysis, with the uncertainty of neighbouring days inflated | REP QUID v4.0 compared drifters and Argo, 2003–2019: drifter−OSTIA SD ≈ 0.4 K, Argo−OSTIA ≈ 0.35 K, vs analysis uncertainty ≈ 0.3 K, "reasonable" after allowing for in situ and representativeness error | PUM 010-011; QUID 010-011; Good et al. (2020), Remote Sens. 12, 720, doi:10.3390/rs12040720; Worsfold et al. (2024), Remote Sens. 16, 3358, doi:10.3390/rs16183358; DOI of the product 10.48670/moi-00168 |
| OISST v2.1 | `err`, °C | 1981-09 onward, 0.25° daily | Bias variance plus random/sampling variance (Reynolds et al. 2007, J. Climate 20, 5473–5496, doi:10.1175/2007JCLI1824.1) | none found (**UNVERIFIED**). It assimilates near-surface Argo, so Argo is not independent of it | Huang et al. (2021), J. Climate 34, 2923–2939, doi:10.1175/JCLI-D-20-0166.1; NCEI dataset doi:10.25921/RE9P-PT57 |
| ESA CCI SST v3 | `analysed_sst_uncertainty`, K | 1980–2021 (climate data record, CDR), 2022 onward (interim CDR, ICDR), 0.05° daily | Same method as OSTIA | Against drifters, the normalised spread is 20–30 % narrower than expected, i.e. it slightly overestimates the uncertainty (or overstates the in situ uncertainty) | Embury et al. (2024), Sci. Data 11, 326, doi:10.1038/s41597-024-03147-w; CEDA doi:10.5285/4a9654136a7148e39b7feb56f8bb02d2 |

**Recommended choice:** OSTIA reprocessed (one source covers the full period, a quality report is available, and it
is accessed with the same toolbox as U1). ESA CCI is the cross-check.

**Extension below the surface (an assumption, delivered separately).** σ is held equal to σ_SST from the surface to
the mixed-layer depth and set to zero below it. Mixed-layer depth comes from the reanalysis: GLORYS12V1 `mlotst`,
the depth where density exceeds its 10 m value by the equivalent of a 0.2 °C temperature drop (PUM 001-030), or
GREP `mlotst_*` (criterion per member **UNVERIFIED**). U2's vertical profile of innovation σ tests this assumption.
The extension is never merged into the U4 product itself.

**Known issues.**
- OSTIA near-real-time `analysis_error` was overestimated from 2013-01-17 to 2014-02-04 because of a bug (NRT PUM
  change log). The reprocessed product is not affected as far as is known (**UNVERIFIED**).
- The near-real-time QUID (v7.0, Table 6, 2014) found the error **underestimated** against near-surface Argo, e.g.
  SD 0.472 K vs a mean analysis error of 0.260 K. Evidence on calibration is mixed, which is why U4 is compared with
  U1 and U2 rather than trusted alone.

---

## U5 — Observation-error covariance for Argo temperature

**What it represents.** The error covariance a DA system would assign to the temperatures of one Argo profile,
i.e. how strongly the observations are believed. It is needed to build U6.

For two levels $k, l$ of one profile:

$$R_{kl} = \sigma_{repr,k}\,\sigma_{repr,l}\,C_{repr}(k,l) + \sigma_{instr}^2\,\delta_{kl},$$

and different profiles are taken as uncorrelated.

- $\sigma_{instr}$ = 0.002 °C (Wong et al. 2020). This is negligible next to $\sigma_{repr}$.
- $\sigma_{repr}^2$ (representativeness): the variability a point measurement sees that a grid cell cannot hold.
  It is computed **for a given target grid**. GLORYS12V1 (1/12°) is averaged onto the target layers, and the
  within-cell variance is taken at each level, averaged over the chosen days.
- $C_{repr}(k,l)$: the correlation of the same within-cell deviations between levels, pooled over the region and
  the days. **It matters.** In the smoke test (January 2015, 33–38 °N, 165–155 °W, GREP grid) it was 1.00 over
  0–47 m: an unresolved eddy shifts the whole winter mixed layer. With a diagonal R, the 14 levels of a profile
  counted as 14 independent observations of what B also treats as one quantity. Increments then reached
  2 background standard deviations; with the correlated R the largest was 0.6.
- Missing: variability finer than 1/12° or faster than daily, so this is a lower bound.
- The values that operational systems prescribe for Argo T (ORAS5: Zuo et al. 2019, Ocean Sci. 15, 779–808,
  doi:10.5194/os-15-779-2019, which inherits them from ECMWF Tech. Memo. 736, doi:10.21957/5awbusgo; GLORYS; FOAM)
  are **UNVERIFIED**. If found, they are reported next to the computed values for comparison.

---

## U6 — Increments driven by real Argo profiles

**What it represents.** What real Argo profiles near a date would change in a typical analysis:

$$\delta x = BH^\top(HBH^\top+R)^{-1}d,$$

with:
- $d$ from U2, using the GREP ensemble mean as background (the builder refuses any other background);
- $R$ from one of three models (a required choice), all correlated between the levels of a profile by U5's
  correlation:
  - **u5_scaled (recommended):** U5's spatial pattern of representativeness variance, scaled at each level so its
    mean over all U2 observations equals U2's fitted observation-error variance;
  - **u5:** U5 as it is. In 2014–15 it was a third of the fitted value, and the analysis was over-confident;
  - **innovations:** U2's fitted value, the same everywhere at each level.
- $B$ from one of two models (a required choice):
  - **separable (recommended):**
    $B_{ij} = \sigma_i\sigma_j \exp(-r_{ij}^2/2L^2)\,C_{ens}(k_i,k_j)$. Here $\sigma$ is the U1 standard
    deviation at each point, $C_{ens}$ the U1 correlation between levels pooled over the region, and $L$ U2's
    fitted background-error correlation length. It is full rank. The single Gaussian length for all levels is
    an assumption.
  - **ensemble:** the U1 covariance multiplied by a Gaspari–Cohn localisation (Gaspari & Cohn 1999, QJRMS 125,
    723–757), with half-width $\sqrt{10/3}\,L$ so that its curvature at zero matches the Gaussian. With three
    members it has rank 2 before localisation. In the smoke test it produced increments of opposite sign to every
    innovation far from the profiles (−0.42 K with all innovations > 0): this is sampling noise, not a feature of
    the ocean. It is kept for comparison with the separable model.

The hold-out evaluation of all six combinations for 2014–15 is in [EVALUATION_2014-15.md](EVALUATION_2014-15.md).
Separable / u5_scaled was chosen. In the frontal Kuroshio region no combination predicts the largest innovations,
and U6's stated analysis uncertainty there is too small by a factor of about 1.5 in variance.

The observation window is set equal to the members' assimilation window (5–7 days; QUID 001-031 Table 1), centred on
the date. It is an argument of the builder.

**Delivered form.** For each date:
- $\delta x$ from all profiles in the window, or from a single profile (`--profile`);
- background and analysis standard deviations, so a user can see the variance reduction;
- the correlation between levels used in B.

**Assumptions.**
- The errors are unbiased and Gaussian, and background and observation errors are uncorrelated.
- $H$ is linear: horizontal bilinear interpolation on the reanalysis levels.
- The ensemble on the centre date is used for the whole window.
- Biases shared by all members (U1) are not removed; they appear in $d$ and are reported in U2.

---

## Excluded, and why

| Candidate | Reason |
|---|---|
| A model forecast minus a reanalysis | It measures one free-running model's drift, not DA-scale error |
| A flux-forced ocean ensemble driven by ERA5 ensemble members | It needs an ocean model, and it estimates only the part of B caused by forcing error, which U1 already includes |
| A lagged-difference B ("NMC method"; Parrish & Derber 1992, Mon. Wea. Rev. 120, 1747–1763) | It measures how fast the state changes, not its error |

---

## Checks of the products

`python -m ocean_uncertainty.checks report` writes checks 1 and 2 (with a comparison of the σ maps).
`... checks holdout` scores the U6 error models on held-out profiles. Check 3 is read off U2 by hand.

1. **Consistency at Argo points** (Desroziers et al. 2005, QJRMS 131, 3385–3396, doi:10.1256/qj.05.108):
   $\overline{d d^\top} = HBH^\top + R$. At each level, the mean of $d^2$ after removing the bias (U2,
   GREP-mean background) is compared with $\sigma_{U1}^2 + \sigma_o^2$ (U5) at the same profiles. A ratio clearly
   above 1 flags U1 as under-dispersive. No inflation is applied silently.
2. **Ordering:** σ(U1) ≲ σ(U2) < σ(U3) at every level, with σ_SST (U4) alongside. Any violation is reported, not
   tuned away.
3. **Published values:** U2 against the GLORYS12V1 QUID text values (the surface bias of about 0.1 °C, RMS < 1 °C in
   the upper 300 m).

Smoke test (January 2015, 33–38 °N, 165–155 °W, 20 profiles; a test of the code, not a result):
- consistency ratio 0.79–0.93 over 5–47 m;
- ordering held at every level with data (U1 0.17, U2 0.23–0.25, U3 0.52–0.59, U4 0.24 K);
- U2 bias +0.13 to +0.20 K (Argo warmer than the GREP mean), in line with the QUID's GREP SST being about 0.2 °C
  colder than OSTIA (QUID 001-031 p.12–14).

---

## Open items

- Observation-error values for Argo T in ORAS5, GLORYS and FOAM (U5).
- Grid and levels of the RG09 files; whether ISAS has an error field (U3 alternatives only).
- The normalisation of GREP's `thetao_std` (U1; not blocking, since σ is recomputed from the members).
- Copernicus Marine licence terms (the acknowledgement is required; the full text was not read).
