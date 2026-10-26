# LS4 field grid

Tools for building the LS4 field grid at La Silla Observatory, scheduling nightly
observing plans, and analyzing survey coverage and cadence.

Plans use a rectangular LS4 field of view, 60-second exposures, 40-second readout,
and airmass, night, and Moon-separation constraints. Every field is observed as an
initial pointing plus a dithered revisit at least 30 minutes later.

## Setup

```bash
conda activate ls4-fields
pip install -e .          # installs the `ls4` command and dependencies
```

Without installing, run the same commands from the repository root as `python -m ls4 ...`.

## Commands

| Command | What it does |
| --- | --- |
| `ls4 plan` | Schedule tonight (or `--date YYYY-MM-DD`) into `plans/YYYYMMDD.csv` and update the grid database |
| `ls4 export plans/X.csv` | Convert a plan to the telescope `.obsplan` format |
| `ls4 viz plans/X.csv` | Interactive plot; add `--gif`, `--visits-gif`, `--coverage`, or `--html` |
| `ls4 analyze [plans...]` | Per-field visits, cadence, nightly summary, and sky maps (default: all of `plans/`) |
| `ls4 bands [plans...]` | Sky coverage and cadence in g, i, and z (also written by `analyze` to `bands/`) |
| `ls4 simulate` | Schedule many nights on a fresh database (`--start`, `--nights`), or `--replay plans/` to re-plan saved nights and compare |
| `ls4 build-grid` | Rebuild the grid CSV and database (`--force` required if the database exists) |
| `ls4 site` | Build the observing-plan website in `outputs/site/` (`--releases` adds published nights, `--serve` previews it) |

Every command has `--help`.

### Nightly workflow

For the scheduled production pipeline, Slack notifications, nightly GIFs, and
La Silla weather reports, see [nightly operations](docs/nightly-operations.md).
The GitHub workflow runs daily at 16:17 UTC once activated on the default branch.
Cadence history currently assumes that generated plans were executed.

```bash
ls4 plan                                  # writes plans/<tonight>.csv
ls4 export plans/20260717.csv -d 20260718 # writes 20260718.obsplan
ls4 viz plans/20260717.csv --gif
```

`ls4 plan` uses the astronomical-twilight window that starts on the given local
evening date. It reads `assets/LS4_field_grid.db` and advances each field's
`last_scheduled_mjd` only when its complete initial/dither pair is in the plan.

### Simulations

```bash
ls4 simulate --start 2026-07-17 --nights 90 --name sim90 --gifs 5
ls4 simulate --replay plans --name replay
```

Each run writes to `outputs/<name>/`: one CSV per night, an isolated `grid.db`,
`summary.csv`, `REPORT.md`, `stats/`, and optional `gifs/`. The project
database and saved plans are never modified. A replay also records the saved plan's
metrics (`baseline_*` columns) next to the new ones.

### Observing-plan website

The website shows each night's plan and the survey's per-band coverage and cadence.
GitHub Actions rebuilds it on GitHub Pages after every nightly run (see
[nightly operations](docs/nightly-operations.md#website)). It has two views:

- **Night plan:** a dropdown of nights, defaulting to tonight. The sky map shows
  footprints colored by time or by g/i/z quadrant, with dithered revisits dashed. It also
  shows the airmass limit around the zenith, the Moon and its exclusion zone, and the
  sky observed on earlier nights. A time slider (with play) steps through the night. The
  page lists the night's metrics, the weather forecast and release for published nights,
  a table of pointings, and links to the plan CSV and `.obsplan`.
- **Coverage & cadence:** the band sky map from `ls4 bands` (nights observed, median
  gap, and days since last, per band) through the selected night. Click a pixel to see
  the nights it was observed in each band. The view also has charts of cumulative area,
  exposures per night, and the distribution of gaps between observations.

Plans are grouped into **Production** (top-level `plans/*.csv`, plus the newest
revision of each published nightly release) and one collection per `plans/forecast_*/`
directory. Links such as `#c=production&n=20261026&v=survey&band=all&metric=gap`
open a specific view.

```bash
ls4 site --releases --serve   # build outputs/site/ and preview it at http://127.0.0.1:8000/
```

## How the scheduler works

The night is split into pairs of 30-minute blocks: the first block takes initial
pointings and the second revisits the same fields, dithered, in the same order.
For each block the scheduler:

- picks a dense, overdue starting field and keeps the rest of the block within
  `max_cluster_radius` of it;
- checks both pointings against all constraints at their start and end times;
- requires each move to fit within readout plus any idle time (`slew_rate`, `read_out`);
- gives SOLE fields that are due priority until they fill 25% of the night's pair
  capacity, without leaving the cluster (see below).

Plans are validated (exactly one initial and one dither per field, at least 30 minutes
apart, no overlapping rows, all slews feasible) before they are saved.

## Dithers and bands

The camera's 2×2 quadrants carry i on one diagonal, z top right, and g bottom left. Each pair
is an initial pointing plus a dither half a grid step west in RA, so every point covered
by the pair gets i plus g or z (the bottom row of quadrants gives g and i, the top row gives i and z).
The pair alternates nightly between two dither phases (`dither_phase`): phase 1 shifts both
pointings half a grid step north in Dec, which swaps g and z. Over two consecutive nights
every part of a field gets i on both nights and g and z once each, matching SOLE's nominal
cadence. Each night uses a single phase so neighboring fields complete each other's coverage;
plans record it in `configuration`. `ls4 analyze` counts fields observed in both phases
(g, i, and z).

### Per-band coverage and cadence

`ls4 bands` (and `ls4 analyze`, into `<output>/bands/`) projects each exposure's four
quadrants onto HEALPix pixels (`--nside`, default 128, ~0.46° pixels) and writes:

- `band_sky_map.html`: interactive map with band (g, i, z, or all three), metric (nights
  observed, median gap, days since last), and a through-night slider. URL parameters
  pick the view, e.g. `band_sky_map.html?band=all&metric=gap&night=10`.
- `band_cadence_cdf.html`: gaps between nights a pixel was observed, per band.
- `band_area_growth.html`: cumulative area per band and in all three bands.
- `band_pixels.csv`, `band_nightly_area.csv`, `band_summary.txt`.

The figure's left/right orientation on the sky is assumed (`QUADRANT_BANDS` in
`ls4/config.py`); a camera flipped in RA mirrors the per-band maps.

## Survey programs

Each field's `program_id` comes from its Galactic coordinates (`ls4.grid.program_ids`,
parameters in `ls4/config.py`):

- **0, SOLE** (Stellar Oscillations, Lensing, and Eruptions): field centers in the plane,
  −72° ≤ l ≤ 36° and |b| ≤ 4°, or the bulge, |l| ≤ 20° and |b| ≤ 8°. That is 55 fields
  (~1,280 deg²), 47 of them above the −60° declination limit.
- **1, extragalactic**: everything else.

SOLE has a nominal 1-day cadence and about 25% of the survey. A SOLE field is due once half
a day has passed since its last pair. Until the night's SOLE budget is used, due SOLE fields are
preferred, and a new cluster starts on the due SOLE field that sets first, provided enough
fields nearby stay visible to fill the block. Programs are derived from `l` and `b` when the
database is read, so databases built with an older footprint need no rebuild.

## Plan CSV format

`target, start time (UTC), end time (UTC), duration (minutes), ra, dec, configuration,
start_time_mjd, block_number, ra_hr`. Dithered visits are named `<field>_dither`.
`configuration` holds `{'dither': bool, 'phase': 0 or 1}`; the `ra`/`dec` columns are the
actual pointings, which are offset from the grid center for dithers and phase 1.
Rows named `Unused Time` (or `TransitionBlock` in older plans) are gaps, not pointings.

## Survey analysis outputs

`ls4 analyze` writes `summary.txt`, `nightly_summary.csv`, `per_field_stats.csv`
(every grid field, including unvisited ones, with sampled visibility), Plotly HTML
histograms, and two self-contained sky maps: `visits_per_field_map.html` and
`visit_gap_percentile_map.html` (P0–P100 slider over gaps between initial visits).
With `--db`, it fails if any field's database history disagrees with the plans.

## Repository layout

```text
ls4/
  config.py     Telescope, field-of-view, constraint parameters, and paths
  grid.py       Grid construction and database access
  scheduler.py  Night window, visibility, and the pair scheduler
  metrics.py    Per-night metrics and multi-night survey statistics
  plots.py      Plotly figures
  animate.py    Sky-rotation GIFs
  bands.py      Per-band (g/i/z) coverage, cadence, and the band sky map
  sky_maps.py   Standalone HTML sky maps
  export.py     .obsplan writer
  website.py    Observing-plan website data (`ls4 site`); the page itself is in web/
  cli.py        The `ls4` command
assets/         Field-grid CSV (the .db is local and gitignored)
plans/          Real nightly plans; plans/archive/ holds older ones
outputs/        Generated runs, stats, and GIFs (gitignored)
tests/          Smoke test (`pytest` or `python tests/test_smoke.py`)
```

Changing values in `ls4/config.py` affects both the grid and the scheduler, so
regenerate dependent outputs afterwards.
