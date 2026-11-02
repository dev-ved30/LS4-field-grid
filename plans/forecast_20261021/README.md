# Forecast plans, 2026-10-21 to 2026-12-20

61 consecutive nights scheduled on 2026-11-02, each from astronomical evening to
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

Totals: 16,076 images in 8,038 complete pairs (6 SOLE); 678 of 1,933 fields
visited; 11,885 deg² covered in g, i, and z.
