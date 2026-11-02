# Forecast plans, 2026-10-21 to 2027-03-31

162 consecutive nights scheduled on 2026-11-02, each from astronomical evening to
morning twilight at La Silla. History was bootstrapped from the top-level
`plans/*.csv` (as the nightly pipeline does) and carried forward night to night,
assuming every plan is executed. Weather losses are not accounted for.
Uses the phase-scaled Moon limit (10 deg floor, ~50 deg at full Moon) and the
-70 deg declination limit.

- `YYYYMMDD.csv`: plan for that local evening.
- `obsplans/`: telescope files, named for the following morning.
- `summary.csv`: nightly metrics.

This directory is outside the production history bootstrap, which reads only
top-level `plans/*.csv`.

Totals: 44,990 images in 22,495 complete pairs (1,170 SOLE); 1,107 of 1,933
fields visited; 21,421 deg² covered in g, i, and z.
