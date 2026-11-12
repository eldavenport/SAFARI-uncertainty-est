# ocean-uncertainty

Model-independent estimates of upper-ocean (0–50 m) temperature uncertainty for the North Pacific, built as a
**monthly climatology** for use as perturbations in adjoint / data-assimilation experiments. The repository builds
and documents these products; it does not apply them to any model.

Exact method for each product: [docs/METHODS.md](docs/METHODS.md). Region 140–240 °E, 20–60 °N; 2011–2022; pooled
by calendar month.

| Product | Content |
|---|---|
| **C1** `c1.nc` | Reanalysis spread: σ(month, depth, lat, lon) = std across the 3-member GREP ensemble (the perturbation amplitude/pattern), the 0–50 m layer mean, and per-member `member − mean` maps |
| **C2** `c2.nc` | Argo minus each member, the ensemble mean, and the Roemmich–Gilson Argo-only field, by calendar month: bias, RMS, std, counts, 2° RMS maps, and a Hollingsworth–Lönnberg split (secondary) |

C2 is the observational check on C1: for a reliable 3-member ensemble `RMS(Argo − member) ≈ √2 × spread`. See
[docs/METHODS.md](docs/METHODS.md) and the figure guide in the products' `figures/README.md`.

## Set-up

```bash
python3.12 -m venv .venv
.venv/bin/pip install -r requirements.txt -e .
.venv/bin/python -m pytest tests
```

Copernicus Marine downloads need an account (`copernicusmarine login`, once).

## Workflow

The pipeline is one module, `ocean_uncertainty.climatology`, run on a plain multi-core machine (no scheduler). Set a
data directory for the (large) downloads and intermediates, and an output directory for the small products.

```bash
D=/path/to/storage/ocean_uncertainty

# 1. Download + accumulate: one task per (year, month), fanned across a bounded process pool. Resumable.
python -m ocean_uncertainty.climatology run --data-dir $D/data --jobs 12

# 2. (optional) Roemmich–Gilson reference for the Argo − RG background: base file + 2019+ monthly extensions
#    from https://sio-argo.ucsd.edu/RG_Climatology.html, gunzipped into $D/data/rg/

# 3. Build the products (pass --rg-dir to add the Argo − RG background to C2).
python -m ocean_uncertainty.climatology build --data-dir $D/data --out-dir $D/products --rg-dir $D/data/rg

# 4. Validation figures.
python scripts/plot_climatology.py --products $D/products --out $D/products/figures
```

Defaults match the first use case (North Pacific, 0–50 m, every 5th day); override with `--region`, `--max-depth`,
`--step-days`, `--years`, `--months`. `python -m ocean_uncertainty.fetch grep|argo` downloads a single source
manually if needed.

## Using a product on a model grid

```python
from ocean_uncertainty.regrid import to_grid
on_model = to_grid(field, lon_edges, lat_edges, depth_edges)   # conservative, over the wet part of each cell
```

Perturb with the **error realisations** — the `member − mean` fields (C1's `thetao_member_anomaly_*`, or the member
anomalies themselves): regrid those and use them directly, since they carry realistic, balanced spatial structure.
Averaging a *variance* instead gives the mean point-error variance, which is larger than the error of a cell mean.

## History

The repository previously built per-winter daily products U1–U6 (reanalysis spread, Argo/climatology innovations,
SST uncertainty, observation-error covariance, and Argo-driven DA increments). That plan is retired; the 2014–15
evaluation that chose the shared settings is kept as historical evidence in
[docs/EVALUATION_2014-15.md](docs/EVALUATION_2014-15.md).
