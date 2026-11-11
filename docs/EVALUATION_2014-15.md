# Evaluation of the product options: North Pacific, winter 2014–15

**Set-up**
- Region: 140–240 °E, 20–60 °N, 0–50 m, 1 November 2014 to 31 March 2015.
- Argo: 8062 delayed-mode profiles (QC 1).
- Products: U1–U5 built with `scripts/run_winter.sh` settings.
- Reproduce with `python -m ocean_uncertainty.checks report|holdout` (commands in that script).

The question was how different the options are, and which one to carry to all winters.

## 1. Error-size products (U1–U5): how different are they?

### Size at the Argo profiles

Ratios are relative to U2, the Argo-minus-reanalysis RMS, which is the only estimate made against independent
point observations. Medians are over 2° boxes with ≥ 10 Argo values. The first two columns are box means of σ.

| Product | σ at 9.8 m (K) | σ at 41.2 m (K) | Ratio to U2 at 9.8 m | Ratio to U2 at 41.2 m |
|---|---|---|---|---|
| U1 reanalysis spread | 0.35 | 0.41 | 0.95 | 1.02 |
| U2 Argo − reanalysis (RMS) | 0.38 | 0.42 | 1 | 1 |
| U3 Argo − climatology (RMS) | 1.01 | 1.02 | 2.66 | 2.35 |
| U4 OSTIA SST uncertainty | 0.29 | 0.29 | 0.89 | 0.80 |
| U5 representativeness (1/12° within 0.25°) | 0.18 | 0.18 | 0.49 | 0.46 |

### Spatial pattern

Spearman rank correlation of the box maps (291–301 boxes):

| | U1 | U2 | U3 | U4 | U5 |
|---|---|---|---|---|---|
| U1 | 1 | 0.76–0.79 | 0.33–0.43 | 0.75–0.78 | 0.83 |
| U2 | | 1 | 0.35–0.39 | 0.66–0.68 | 0.70 |
| U3 | | | 1 | 0.22–0.24 | 0.30–0.31 |

### Error budget at the profiles

Per level, from U2's Hollingsworth–Lönnberg split:
- innovation standard deviation 0.47–0.57 K;
- background part σ_b ≈ 0.32–0.37 K;
- observation part σ_o ≈ 0.35–0.40 K;
- correlation length L ≈ 50 ± 10 km.

These are robust to bin width (10, 25, 50 km) at a 3-day pairing window. L grows to about 65 km at 7 days,
because errors also decorrelate in time.

The consistency ratio of innovation variance to (U1 variance + U5 variance) is 1.28–1.34 below 5 m. U1's variance
matches σ_b², while U5's variance is a third of σ_o². **The reanalysis spread is the right size; the high-resolution
representativeness estimate is too small** (variability finer than 1/12° and GLORYS12's own errors are missing).

### Reading of section 1

- **U1, U2, U4 and U5 agree in pattern** (ρ 0.66–0.86): largest in the Kuroshio and its extension.
- **U1 and U2 also agree in size** (ratio 0.95–1.02).
- **U3 is a different kind of quantity.** It is 2.4–2.7 times larger, with a weakly related pattern (ρ 0.22–0.43):
  large-scale anomalies from the climatology rather than analysis error. It stays the upper bound.
- **U4 is a surface-only estimate.** It is 10–20 % smaller than U2, and its size does not change with depth.

## 2. Increment models (U6): how different are they?

Error models compared:
- **B:** `separable` (U1 σ, U1 correlation between levels, Gaussian with U2's L) or `ensemble` (U1 covariance,
  Gaspari–Cohn localisation);
- **R:** `u5` (U5 representativeness), `innovations` (U2's σ_o² per level, uniform in space) or `u5_scaled` (U5's
  pattern scaled per level to U2's σ_o²). All three use U5's correlation between levels.

**Method.**
- 21 non-overlapping 7-day windows over the winter.
- Each window's profiles are split into 5 folds. Each fold's innovations are predicted from the other four
  (`checks.hold_out`).
- Scores:
  - skill = 1 − Σ(d − Hδx)²/Σd²;
  - median ratio = median |d − Hδx| / median |d|;
  - calibration = mean (d − Hδx)² / mean expected variance (1 means the stated uncertainty is right).

**Profiles within 100 km of a profile used.** Background σ classes are terciles. In each cell: skill / median ratio / calibration.

| Model | σ_b ≤ 0.11 K | 0.11 < σ_b ≤ 0.22 K | σ_b > 0.22 K |
|---|---|---|---|
| separable / u5 | −0.03 / 0.96 / 2.82 | −0.03 / 0.84 / 2.86 | −0.08 / 0.81 / 2.31 |
| separable / innovations | 0.02 / 0.94 / 0.52 | 0.02 / 0.80 / 0.67 | −0.61 / 0.84 / 3.41 |
| **separable / u5_scaled** | 0.01 / 0.94 / 0.95 | 0.01 / 0.81 / 1.17 | −0.25 / 0.85 / 1.56 |
| ensemble / u5 | −0.01 / 0.94 / 2.75 | −0.15 / 0.85 / 3.15 | −0.24 / 0.88 / 2.60 |
| ensemble / innovations | 0.04 / 0.95 / 0.50 | −0.11 / 0.83 / 0.76 | −0.47 / 0.85 / 3.04 |
| ensemble / u5_scaled | 0.03 / 0.95 / 0.93 | −0.13 / 0.86 / 1.34 | −0.22 / 0.85 / 1.51 |

Profiles more than 100 km from any used profile are not changed by any model (skill about 0, median ratio about 1),
as expected with L ≈ 50 km. There, calibration tests R and σ_b alone:

| R model | Calibration range (> 100 km) |
|---|---|
| u5 | 1.25–2.26 |
| innovations | 0.39–1.06 |
| u5_scaled | 0.84–0.93 |

Full table: `holdout.md` of the 2014–15 build (regenerated on the server, docs/SERVER_PLAN.md step 2).

### Reading of section 2

- **The choice of B matters little.** Separable and ensemble differ by ≤ 0.04 in median ratio and by ≤ 0.15 in
  skill. The largest difference is in the middle σ class, where the ensemble loses (−0.13 against +0.01).
- **The choice of R matters a lot for the stated uncertainty.**
  - `u5` is over-confident everywhere (residuals 1.3–3.2 times their expected variance).
  - `innovations` is under-confident where errors are small and over-confident at fronts, because one value per
    level cannot follow their spatial variation.
  - `u5_scaled` is within 0.84–1.56 everywhere.
- **The typical innovation is predicted by every model** (median ratio 0.80–0.96 within 100 km).
- **No model predicts the largest errors.** These are in the frontal Kuroshio region (σ_b > 0.22 K), where squared
  skill is negative for all models. Profiles 2 km apart there differ in innovation by up to 4.8 K. Two reasons:
  - The errors there are front displacements, not smooth Gaussian fields.
  - Three-member σ (up to 4 K) scales the predictions.
- **The increments at fronts should be read with this in mind.** The analysis uncertainty U6 reports there is too
  small by a factor of about 1.5 in variance.

## 3. Choice for all winters

- **U1–U5:** built for every winter as they are. They are complementary, not alternatives:
  - U1 gives realisations and the right size;
  - U2 is the independent observation check;
  - U3 is the upper bound;
  - U4 is the surface analysis estimate;
  - U5 is the correlation structure of the observation error.
- **U6: `separable` B with `u5_scaled` R.** It has the best calibration of the six (1.18 overall; ensemble /
  u5_scaled 1.17). Its median ratio is within 0.003 of the best, it is full rank, and it avoids the ensemble's losses in
  the middle σ class.
- **Kept the same in every winter:** the statistics settings (25 km bins, ≤ 1000 km, ≤ 3 days apart, ≥ 100 pairs,
  2° maps) and the 7-day windows, so that winters can be compared.
- **Not done, recorded as possible improvements:**
  - a time correlation between profiles of one window;
  - smoothing the three-member σ;
  - flow-dependent (front-following) covariances.
