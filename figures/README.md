# What each figure shows

Plain-language guide to the figures in this folder. They come from two products:

- **C1 — the "spread" of the three ocean reanalyses.** We have three different ocean reconstructions (model estimates
  of past ocean temperature) for the same region. At every point we take how much the three disagree with each other
  (their standard deviation). That disagreement is our model-based estimate of how uncertain the ocean temperature is.
- **C2 — reanalysis vs. real floats.** We compare each reconstruction (and their average) against independent Argo
  float measurements. "Argo − member" means a real float measurement minus one reconstruction at the same place and
  time; "Argo − mean" uses the average of the three.
- **RG — a no-model reference (Roemmich–Gilson).** A smooth temperature map built only from Argo floats, no ocean
  model at all. It is coarse (1° grid) and monthly, so it captures only the large, slow features. "Argo − RG" is a
  real float minus that smooth float-only map.

Region is the North Pacific (140–240°E, 20–60°N), upper 0–50 m, every calendar month pooled over 2011–2022.

---

### fig1 — C1 spread, by season (four maps)
**Shows:** the typical size of the temperature uncertainty (in °C) across the map, for winter, spring, summer, autumn.
**Tells us:** where the ocean is most uncertain. It lights up along the Kuroshio Current and its extension off Japan —
the strong, eddy-rich current — and is largest in winter and spring. This is the pattern we would perturb with.

### fig2 — C1 spread: average, smallest month, largest month (three maps)
**Shows:** the year-round average uncertainty, and the quietest-month and busiest-month extremes at each location.
**Tells us:** how much the uncertainty swings through the year. A big gap between the min and max map means that place
has a strong seasonal cycle in how uncertain it is.

### fig3 — C1 structure (four panels)
**Shows:** (a) the region-average uncertainty month by month; (b) how it changes with depth through the year;
(c) how it varies with latitude; (d) how it varies with longitude.
**Tells us:** the uncertainty is largest near the surface and in late summer/early autumn, and is concentrated at
mid-latitudes (the current system) and in the west (near Japan).

### fig4 — Is "member − mean" a realistic perturbation? (two panels) — the key validation
**Shows:** month by month at ~10 m depth: the spread of the reconstructions (black), how far a real float sits from
one reconstruction (Argo − member), from their average (Argo − mean), and from the no-model float-only map
(Argo − RG, red). The grey lines are what we would expect if the spread were a perfect estimate of the error (a
reliable set of three gives Argo − member ≈ 1.41 × spread and Argo − mean ≈ 1.15 × spread).
**Tells us:** in winter the real misfits land almost exactly on the expected grey lines — so the spread is a realistic
perturbation size. In summer they run higher than expected. The red no-model line explains why: in summer it jumps
far above the reconstructions, meaning the summer upper ocean is full of small, fast features that neither a coarse
float-only map nor a model grid can capture. So the summer excess is mostly this fine-scale "a float sees more than a
grid box" effect — not the spread being too small. The reconstructions genuinely beat the float-only map (they add
real detail), especially in summer.

### fig5 — The perturbation field each reconstruction provides (three maps)
**Shows:** for each of the three reconstructions, how far it sits from the three-member average (member − mean),
as a winter 0–50 m map.
**Tells us:** this *is* the actual perturbation you would apply — a realistic, spatially-smooth temperature
disturbance concentrated in the Kuroshio region. Each member gives a slightly different but similar pattern.

### fig6 — Real floats minus each reference (five maps)
**Shows:** winter maps (on a coarse 2° grid, wherever floats sampled) of how far Argo floats sit from the average, from
each of the three reconstructions, and from the no-model float-only map (RG), at ~10 m.
**Tells us:** where each reference actually misses the real ocean. It is patchy because it only exists where floats
happened to be. The four model panels look alike and are smaller than the RG panel — the models track the floats more
closely than the smooth float-only map does.
