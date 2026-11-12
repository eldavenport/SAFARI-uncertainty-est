# Methods — exactly what was done for each product

Authoritative method for the C1 and C2 products built by `ocean_uncertainty.climatology`. A copy ships with the
products (as `README.md` beside `c1.nc`/`c2.nc`); the figure guide is `figures/README.md`.

## Common settings (every case)

- **Region:** 140–240 °E, 20–60 °N. **Depth:** 0–50 m. **Variable:** potential temperature `thetao`.
- **Period / sampling:** 2011–2022. Within each month, days 1, 6, 11, 16, 21, 26 are used (every 5th day, capped at
  day 28 so February is sampled identically): 6 days/month × 12 years = **72 sampled days per calendar month**.
  Results are pooled **by calendar month** (all 12 Januaries together, etc.) into 12 monthly fields.
  - *Why this sampling.* With 3 members each day contributes 2 degrees of freedom, so 72 days give at most 144; the
    relative standard error of σ is ≥ 1/√(2·144) ≈ 6 % (larger in practice, since days 5 apart are not independent).
    Every 10th day would halve the cost at ≈ 8 %.
- **Reanalysis input:** GREP, the 3-member Copernicus ensemble (GLORYS2V4 = `glor`, C-GLORS = `cglo`,
  ORAS5 = `oras`), dataset `cmems_mod_glo_phy-all_my_0.25deg_P1D-m`, daily, 0.25°. One day per request.
- **Argo input:** argopy "research" mode (delayed-mode, QC flag 1), to 100 dbar. Each profile is converted
  pressure→depth (TEOS-10) and linearly interpolated in depth to the GREP levels; never extrapolated.
- **Cost (this server, cached inputs recompute in ~1 min):** the full download run is ~20 min across 12 workers
  (~15 MB per member-day); downloads dominate, compute is negligible.

## C1 — reanalysis spread (`c1.nc`)

For each sampled day, the variance across the 3 members is computed with N−1 = 2 in the denominator. Per calendar
month we accumulate the sum of those daily variances and the day count, then

    sigma(month, depth, lat, lon) = sqrt( (sum of daily 3-member variances) / (number of sampled days) )

i.e. the root-mean-square of the member anomalies `(member − mean)` across members, pooled over the sampled days.
This is a **standard deviation (spread)**, not a mean — the mean of `member − mean` is zero by construction.

**Variables.**
- `thetao_std (month, depth, lat, lon)` — spread on the 0.25° grid, 18 levels in 0–50 m.
- `thetao_std_0_50m (month, lat, lon)` — 0–50 m layer value: the 0–50 m mean variance, square-rooted
  (conservative vertical average, `regrid.to_grid`).
- `thetao_member_anomaly_0_50m (month, member, lat, lon)` — per member, the 0–50 m RMS of `(member − mean)` over the
  sampled days: the per-member perturbation realisation. Linked by `thetao_std² = 3/2 × mean_over_members(member_anomaly²)`.
- `n_days (month)` — sampled days pooled into each calendar month (72).

**What it is / limitation.** A model-based, Argo-independent error estimate — the realistic pattern and size to
perturb with. The 3 members share forcing, model family and assimilated Argo, so the spread captures only the part
of the error where they *disagree*; error **common** to all members cancels in `member − mean` and is invisible
here. C1 is therefore a lower bound (under-dispersive); C2 is the observational check on its size.

## C2 — Argo minus each background, by calendar month (`c2.nc`)

On each sampled day, every matched Argo profile is differenced against **each member**, the **ensemble mean**, and the
**Roemmich–Gilson** field (below). Only profiles on a sampled day are kept. Pooled by calendar month, then per
background and month, per level: bias, RMS, std, count; 2° RMS maps; and a Hollingsworth–Lönnberg split (a Gaussian
fitted to the innovation pair-covariance vs separation → background variance + correlation length + an
observation/representativeness nugget — a **secondary diagnostic**, noisy at single-month resolution).

Statistics parameters: max pair separation 1000 km, 25 km bins, pairs within 3 days, ≥100 pairs to fit, 2° maps.

The `background` dimension holds `mean`, `glor`, `cglo`, `oras`, `rg`. Key variables: `rms`, `std`, `bias`, `count`,
`rms_map`, `length_km`, `background_variance`, `observation_variance` (by `background, month, depth`), and
`c1_sigma (month, depth)` — the C1 spread sampled at the same profiles.

**What it tells us (validation).** Whether `member − mean` is a realistic perturbation. For a reliable 3-member
ensemble, `RMS(Argo − member) ≈ √2 × spread` and `RMS(Argo − mean) ≈ √(4/3) × spread`. Found: on these lines in
winter (spread is a realistic amplitude), higher in summer (traced mostly to representativeness — see RG).
**Irreducible limitation:** because every reanalysis assimilates Argo, error common to the reanalyses *and*
consistent with Argo is invisible to C1 **and** C2 — a floor set by verifying against assimilated observations.

### RG background — Roemmich–Gilson (model-free reference)

RG is a 1°, monthly, Argo-only objective analysis (no ocean model). The field for a year-month is
`ARGO_TEMPERATURE_MEAN + ARGO_TEMPERATURE_ANOMALY[month]` (base file 2004–2018; monthly extension files 2019–2022,
from <https://sio-argo.ucsd.edu/RG_Climatology.html>). Pressure→depth at the region mid-latitude, cut to region,
interpolated to the GREP levels; then `Argo − RG` is sampled at the stored profiles and added as the `rg` background.

Because RG is coarse (1°) and monthly, `Argo − RG` is *larger* than `Argo − member` — an inflated upper bound, not a
tight representativeness estimate. Its value is seasonal corroboration: `Argo − RG ≈ Argo − member` in winter but ≫
in summer, confirming the summer excess over the spread is fine-scale/sub-monthly representativeness (a point float
sees more than a grid box), not ensemble under-dispersion.

## C3 — Argo analysis increment (`c3.nc`)

How much a month of Argo adjusts the upper ocean, i.e. the DA increment `δx = B Hᵀ (H B Hᵀ + R)⁻¹ d`, with
`d = Argo − background`. For a point near one float this is `δx ≈ [σ_b²/(σ_b²+σ_o²)]·d` reaching out a length `L`,
so the spread σ_b both sets the gain and maps where Argo has leverage.

**What was done.** For each (year, month): a 2-D optimal interpolation at each depth level independently (vertical
correlation neglected). `σ_b` = C1 `thetao_std` for that calendar month, bilinear at each point; horizontal
background covariance `σ_b σ_b' exp(−r²/2L²)` with `L` = **prescribed 50 km** (C2's per-month HL `length_km` is too
noisy to use directly — NaN to >1000 km — so the documented ~50 km winter length is used; override with
`--length-km`); `R` diagonal from C2 `observation_variance` (or a chosen value); innovation `d` = `Argo − mean`. Profiles within 0.5° are
**superobbed** (nearby Argo is not independent), and innovations beyond 3σ per level are **gross-error-rejected**
(frontal/representativeness outliers the smooth grid cannot carry) — both prevent the increment overshooting. Pooled
by calendar month.

**Variables** (by `month, depth, lat, lon`): `increment_mean` (systematic adjustment), `increment_rms` (typical
magnitude), `background_std` (= C1 σ_b), `analysis_std`, and `uncertainty_reduction = 1 − analysis_std/background_std`.

**What it tells us.** The climatological (year-averaged) reduction reaches ~25–30 % at the best-sampled spots and
~1 % region-mean — sparse Argo at a ~50 km reach. In any *individual* month the reduction right at a float is much
larger (~70–90 %); it averages down here because floats sit in different places each year (use the per-(year, month)
`increment` for a specific case). **Caveat:** structural-spread σ_b and a single `L` make this a realistic *magnitude
and reach* of the adjustment, not a specific operational increment; the same assimilated-Argo floor as C2 applies.

## Deliberately not built

- **Argo − WOA23 climatology** — dropped; RG (monthly, obs-only) is the better model-free reference. (The "C3" name
  now denotes the increment above.)
- **U1–U6** (earlier per-winter daily products with DA increments) — retired with the per-winter plan; the 2014–15
  evaluation that chose the shared settings is kept as evidence in [EVALUATION_2014-15.md](EVALUATION_2014-15.md).
- **ORAS5 5-member ensemble** — considered and rejected: coarse/monthly, still Argo-suppressed, and not independent
  (ORAS5 is already one of the three GREP members).

## Citations

GREP: Copernicus Marine, DOI 10.48670/moi-00024. Argo: DOI 10.17882/42182. Roemmich & Gilson (2009),
Prog. Oceanogr. 82, 81–100, DOI 10.1016/j.pocean.2009.03.004. Each product file also carries a `sources` attribute.
