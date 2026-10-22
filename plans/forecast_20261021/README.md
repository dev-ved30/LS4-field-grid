# Forecast plans, 2026-10-21 to 2026-12-20

61 consecutive nights scheduled on 2026-10-22, each from astronomical evening to
morning twilight at La Silla. History was bootstrapped from the top-level
`plans/*.csv` (as the nightly pipeline does) and carried forward night to night,
assuming every plan is executed. Weather losses are not accounted for.

- `YYYYMMDD.csv`: plan for that local evening.
- `obsplans/`: telescope files, named for the following morning.
- `summary.csv`: nightly metrics.

This directory is outside the production history bootstrap, which reads only
top-level `plans/*.csv`.

Totals: 16,056 images in 8,028 complete pairs (3 SOLE); 639 of 1,933 fields
visited; 11,440 deg² covered in g, i, and z.
