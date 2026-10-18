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
| `ls4 simulate` | Schedule many nights on a fresh database (`--start`, `--nights`), or `--replay plans/` to re-plan saved nights and compare |
| `ls4 build-grid` | Rebuild the grid CSV and database (`--force` required if the database exists) |

Every command has `--help`.

### Nightly workflow

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

## How the scheduler works

The night is split into pairs of 30-minute blocks: the first block takes initial
pointings and the second revisits the same fields, dithered, in the same order.
For each block the scheduler:

- picks a dense, overdue starting field and keeps the rest of the block within
  `max_cluster_radius` of it;
- checks both pointings against all constraints at their start and end times;
- requires each move to fit within readout plus any idle time (`slew_rate`, `read_out`);
- prefers galactic fields until they make up about 25% of pairs, without leaving the cluster.

Plans are validated (exactly one initial and one dither per field, at least 30 minutes
apart, no overlapping rows, all slews feasible) before they are saved.

## Plan CSV format

`target, start time (UTC), end time (UTC), duration (minutes), ra, dec, configuration,
start_time_mjd, block_number, ra_hr`. Dithered visits are named `<field>_dither`.
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
  sky_maps.py   Standalone HTML sky maps
  export.py     .obsplan writer
  cli.py        The `ls4` command
assets/         Field-grid CSV (the .db is local and gitignored)
plans/          Real nightly plans; plans/archive/ holds older ones
outputs/        Generated runs, stats, and GIFs (gitignored)
tests/          Smoke test (`pytest` or `python tests/test_smoke.py`)
```

Changing values in `ls4/config.py` affects both the grid and the scheduler, so
regenerate dependent outputs afterwards.
