# Execution plan: monthly climatology of North Pacific upper-ocean temperature error (CURC Alpine)

**Revised 2026-11-10 (user):** the product is a monthly climatology, not a day-by-day product. All 12 calendar months
are covered, from 2011 to 2022 (the years of the AR cases). The aim is a realistic pattern and size of
upper-ocean (0–50 m) temperature error to perturb with. Region: 140–240 °E, 20–60 °N.

The earlier per-winter plan (U1–U6 daily, with weekly U6 increments) is retired. Its 2014–15 evaluation
([EVALUATION_2014-15.md](EVALUATION_2014-15.md)) stays as the evidence for the choices below.

## What is built

| Product | Content | Why |
|---|---|---|
| **C1** reanalysis spread climatology (main product) | σ(month, depth, lat, lon) = √(mean over sampled days of the 3-member variance), on the 0.25° GREP grid; also the 0–50 m layer mean | In 2014–15 U1 matched the independent Argo check in size (0.95–1.02) and pattern (ρ ≈ 0.77) |
| **C2** Argo minus reanalysis mean, by month | per month: bias, RMS, count; C1 σ at the same profiles and the ratio; Hollingsworth–Lönnberg split and correlation length; 2° RMS maps | Independent check of C1's size each month. The fitted length (≈ 50 km in winter) is the length to use for random perturbation fields |
| **C3** Argo minus WOA23 climatology, by month | as C2 | Upper bound (≈ 2.5 × C2 in 2014–15) |

**Dropped:** U4, U5, U6 and the hold-out comparison. They are needed only for date-specific increments. U4
(SST uncertainty) can be added later as a monthly climatology if wanted.

## Measured cost (laptop, 2026-11-10)

| Access | Time |
|---|---|
| Reanalysis ensemble (3 members), one day per request, `service="arco-geo-series"`, 0–60 m, whole region | 17.7 s per day (15 MB), including about 6 s to open |
| The same days taken from a month opened at once | 52–101 s per day, so not used |
| Argo (argopy "research", ERDDAP), one month, whole region | 50 s (47 000 points) |

**Sampling every 5th day:** 6 days per month × 12 years = 72 days per calendar month (876 days in all).
- Downloads in total: about 13 GB (reanalysis) + about 1 GB (Argo) + 6.7 GB (WOA23 global files, cut at once).
- Nothing large is kept: each day's variance is added to its month's total and the fields are discarded.
- Kept on disk:
  - about 1.3 GB of monthly partial sums (144 × 9 MB);
  - about 1 GB of Argo points;
  - the final products, under 0.2 GB.
- Run sequentially: 4.4 h (reanalysis) + 2 h (Argo), above the 3 h limit for the laptop.
- As a server array of 144 (year, month) tasks: about 3 min each.

## Laptop and server are different machines

Everything so far was developed and run on the laptop. The server is a different architecture, so nothing below may
be assumed to work there because it worked on the laptop.

| | Laptop (where the code and EVALUATION_2014-15.md were made) | Alpine `acpu` nodes |
|---|---|---|
| OS, CPU | macOS (Darwin), arm64, Apple M3 Max | Linux, x86_64 (confirm with `uname -m` and `lscpu`) |
| Python packages | pip wheels in `.venv` on Homebrew Python 3.12 | conda-forge Python 3.12 + pip (`scripts/setup_env.sh`); versions may resolve differently |
| Linear algebra (BLAS/LAPACK) | Apple Accelerate | whatever the Linux numpy/scipy build links (OpenBLAS or MKL; check `numpy.show_config()`) |
| Shell | bash 3.2; no `flock` | bash 4 or 5; `flock` present |
| Network | home/campus connection to Copernicus Marine and Ifremer ERDDAP | CURC network |

**What follows from this:**
- **New code is written and tested on the server, not on the laptop.** The laptop is not used to develop, smoke-test
  or run any of the steps below.
- **The tests must pass on the server.** They were only ever run on the laptop.
  - The random draws are reproducible across machines (NumPy's default generator).
  - Linear algebra may round differently, though.
  - If a test fails there, find out why rather than widening its tolerance.
- **The server environment is recorded.** Save `pip freeze`, `numpy.show_config()`, `uname -a` and `lscpu` from the
  env to `docs/ENVIRONMENT_alpine.txt` once it is built.
- **Laptop numbers carry no assumption to the server:**
  - Download times are re-measured in the pilot.
  - The 2014–15 numbers in EVALUATION_2014-15.md are compared with a tolerance, not bit for bit. Different BLAS
    builds change the last digits of the analysis and hold-out results.
  - The Argo archive may also have changed. Compare profile counts first.
- **The shell scripts have never run on Linux.** `setup_env.sh` and `paths.sh` were only syntax-checked with macOS
  bash 3.2, and no `sbatch` file has been submitted anywhere. The `flock` call in the retired `run_winter.sh` was
  never run, since macOS has no `flock`.

## Code to write on the server (before the pilot)

1. **`ocean_uncertainty/climatology.py`**, with three functions.
   - `month_task(year, month, step_days, ...)` runs one array task.
     - It downloads each sampled day, opening a one-day range on the geo-series service.
     - From the members it accumulates Σ variance (N − 1, double precision) and the number of days.
     - It fetches that month's Argo, interpolates the profiles to the reanalysis levels, and keeps the
       innovations against the ensemble mean on the sampled days, with that day's σ.
     - Outputs:
       - `spread_YYYYMM.nc`;
       - `argo_YYYYMM.nc`, all profiles, kept for C3;
       - `innovations_YYYYMM.nc`.
     - Existing files are skipped, so a task resumes.
   - `build_c1(spread_files)`: pools the partial sums by calendar month. It reuses
     `regrid.to_grid(depth_edges=[0, 50])` for the layer mean.
   - `build_c2(innovation_files)` and `build_c3(argo_files, woa_files)`: by calendar month, reusing
     `products.innovation_statistics`, `products.innovations_monthly` and `checks`-style comparisons.
2. **Tests:**
   - pooled monthly σ from partial sums equals σ computed from all the sampled days at once;
   - the sampled-day list for a month is the stated rule, e.g. days 1, 6, 11, 16, 21, 26 for a 5-day step.
3. **Smoke test on the server** (an interactive `acompile` or `sinteractive` session; not the login node, and not the
   laptop): one month, a 5° × 5° box. Check C1 against the scale recorded in EVALUATION_2014-15.md for 2014–15
   (0.35–0.41 K box mean). Those are winter values for the whole region, so this is a check of order of magnitude
   only.
4. **Job scripts:**
   - `scripts/climatology_month.sbatch`: array over 144 (year, month) tasks, 1 CPU.
     - Memory: one day of members, ~15 MB, plus one month of Argo.
     - Time: about 3 min each; request 30 min until measured.
   - `scripts/climatology_build.sbatch`: one task building C1–C3 after the array finishes (`afterok`).
     - WOA23 is fetched in this task, or once beforehand with `fetch woa23 --months 1 ... 12`.

`scripts/setup_env.sh` and `scripts/paths.sh` are kept. `scripts/run_winter.sh`, `fetch.sbatch` and
`build.sbatch` are retired with the per-winter plan and are removed when the climatology scripts replace them.

## Steps on Alpine

0. **Prerequisites:**
   - an Alpine association (`sacctmgr show user $USER withassoc` lists `alpine`);
   - the repository in `/projects/$USER/SAFARI-uncertainty-est`, copied with `rsync`, leaving out `.venv`;
   - Copernicus credentials (`copernicusmarine login`, or `~/.copernicusmarine_env`, `chmod 600`).
1. **Environment:** `acompile`, then `bash scripts/setup_env.sh`. It must end with all tests passing and "ENV OK".
   Record the environment (see "Laptop and server are different machines").
2. **Pilot:** a few months (`--array` over a handful of task indices).
   - Check that the reanalysis download takes about 18 s per day on Alpine.
   - Record `sacct` Elapsed and MaxRSS, and set `--time` and `--mem`.
   - Check that the pilot month's C1 and C2 agree with the server smoke test over the smoke-test box.
3. **Full array, then the build job:**
   - Check that C1's monthly σ and C2's RMS have a ratio near 1 every month. In 2014–15 it was 0.95–1.02. A
     month far from 1 is reported, not adjusted.
   - Check that C3 is above C2 in every month.

## Open choices

- **Day sampling.** Every 5th day is the plan. With three members, each day adds 2 degrees of freedom, so 72 days
  per calendar month gives at most 144. The relative standard error of σ is then ≥ 1/√(2·144) ≈ 6 %, larger in
  practice because errors on days 5 apart are not independent. Every 10th day halves the cost at about 8 %.
- **Download speed on Alpine.** It is unknown until the pilot. If it is much slower than the laptop's 18 s per day,
  the array's task time and the choice of day sampling are revisited.
