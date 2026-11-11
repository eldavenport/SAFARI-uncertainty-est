# ocean-uncertainty

Model-independent estimates of upper-ocean temperature uncertainty in a typical ocean data assimilation (DA)
system, and of how much real Argo observations would adjust such an analysis. The repository builds and documents
these products. It does not apply them to any model.

What each product is, where its inputs come from, and its limitations: [docs/PRODUCT_CATALOGUE.md](docs/PRODUCT_CATALOGUE.md).

| Product | Content |
|---|---|
| U1 | Spread of the three-member Copernicus ocean reanalysis ensemble (GREP): member anomalies and standard deviation |
| U2 | Argo minus a reanalysis: matched table, bias and RMS by level, error split into background and observation parts, correlation lengths, maps |
| U3 | Argo minus the World Ocean Atlas 2023 climatology: the same statistics as U2 (an upper bound) |
| U4 | OSTIA SST analysis uncertainty, and separately its extension over the mixed layer |
| U5 | Observation-error covariance of an Argo profile for a chosen grid (representativeness + instrument) |
| U6 | Increment that the Argo profiles near a date would make to the reanalysis-ensemble mean, with analysis uncertainty |

## Set-up

```bash
python3.12 -m venv .venv
.venv/bin/pip install -r requirements.txt -e .
.venv/bin/python -m pytest tests
```

Copernicus Marine downloads need an account (`copernicusmarine login`, once).

## Workflow

`scripts/run_winter.sh fetch|build YEAR` (on a server: [docs/SERVER_PLAN.md](docs/SERVER_PLAN.md)) runs the first use case: the North Pacific (140–240 °E, 20–60 °N), 0–50 m,
November of YEAR to March of YEAR + 1. Its settings, and the evaluation that chose them, are in
[docs/EVALUATION_2014-15.md](docs/EVALUATION_2014-15.md). For another region, period or depth, copy the script and
change its settings. Every statistics argument is required, so each choice is recorded in the products.

The pipeline runs in four steps:
1. `python -m ocean_uncertainty.fetch <source>` downloads each source to its own directory.
2. `python -m ocean_uncertainty.products u1|innovations|u4|u5` builds the products.
3. `... checks report` compares them with each other at the Argo profiles. `... checks holdout` scores the U6
   error models on held-out profiles.
4. `... products u6 --background-model separable --obs-model u5_scaled` makes the increments for one date.

## Using a product on a model grid

```python
from ocean_uncertainty.regrid import to_grid, cell_bounds, level_interfaces
on_model = to_grid(field, lon_edges, lat_edges, depth_edges)   # conservative, over the wet part of each cell
```

Regrid error *realisations* (U1 anomalies, U6 increments) and recompute their spread on the model grid.
Averaging a variance instead gives the mean point-error variance, which is larger than the error of a cell mean.
