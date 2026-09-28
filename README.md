# LS4 field grid

Tools for building a field grid and generating, visualizing, and evaluating observing plans for the LS4 instrument at La Silla Observatory.

The scheduling workflow uses a rectangular LS4 field of view, 60-second exposures, 40-second readout time, and visibility constraints based on airmass, astronomical night, and Moon separation. The default observing strategy schedules an initial pointing followed by a dithered revisit.

## Workflow

The main workflow is:

1. Build the field grid and its SQLite database.
2. Generate an observing plan for a night.
3. Visualize or animate the plan.
4. Aggregate multiple plans to inspect coverage and cadence.

## Requirements

The scripts require Python 3 and astronomy/scientific Python packages, including:

- `astropy` and `astroplan`
- `regions`
- `healpy`
- `m4opt`
- `numpy` and `pandas`
- `Pillow`, `plotly`, and `tqdm`

Install the packages in the Python environment used to run the scripts. This repository does not currently include a `requirements.txt` or packaging configuration.

## Build the field grid

The checked-in grid is available as `assets/LS4_field_grid.csv`. The local SQLite
database at `assets/LS4_field_grid.db` is generated separately. To rebuild the grid
and create the database from scratch:

```bash
python build_field_grid.py
```

The script assigns each field a grid name, equatorial and Galactic coordinates, and a program ID. It also produces a test field-grid visualization. The database table used by the scheduler is named `grid`.

## Generate an observing plan

Generate a plan using the current date:

```bash
python generate_obs_plan.py
```

Specify an observing date with an MJD and optionally choose the output path:

```bash
python generate_obs_plan.py --mjd 61238.5 --output plans/example.csv
```

If `--output` is omitted, the plan is written to `plans/YYYYMMDD.csv`. The scheduler reads the field grid from `assets/LS4_field_grid.db` and updates `last_scheduled_mjd` after writing the plan.

The scheduler chooses fields separately for each pair of observing blocks. It keeps each block's fields within the `max_cluster_radius` set in `constants.py` and favors nearby successive pointings. It checks both pointings at their planned start and end times, schedules the revisit at least 30 minutes after the first visit, and advances a field's database history only when the complete pair is present. Consecutive observations must allow for exposure time plus the larger of slew time and readout time; `slew_rate` and `read_out` are configurable in `constants.py`. Use `--db` to select another grid database and `--no-plots` for batch runs.

Plan CSVs contain target names, UTC start/end times, durations, RA/Dec, dither configuration, MJD timestamps, block numbers, and RA in hours. Rows with `Unused Time` or `TransitionBlock` describe schedule gaps or transitions rather than science pointings.

## Visualize a plan

Create the interactive observing-plan plot:

```bash
python vizualize_obs_plan.py plans/20260717.csv
```

The filename contains the historical spelling `vizualize` and should be used as written.

Create a GIF showing the plan as the sky rotates:

```bash
python make_observing_plan_gif.py plans/20260717.csv \
  --output observing_plan.gif
```

Create a GIF showing cumulative visits to each field:

```bash
python make_visit_count_gif.py plans/20260717.csv \
  --output visit_counts.gif
```

Both GIF scripts accept `--cadence-minutes`, `--duration`, `--width`, `--height`, and `--no-galaxy`. They use `assets/LS4_field_grid.csv` by default; an alternative grid can be supplied with `--field-grid`.

## Replay and compare the scheduler

Replay the nights in `plans/` from a fresh copy of `assets/LS4_field_grid.csv`:

```bash
python replay_plans.py
```

The script writes replacement plans, an isolated SQLite database, `summary.csv`, and
full-period field statistics in `stats/` to `replay_results/`. It leaves the original
plans and database untouched. Use `--limit 1` for a quick check or `--output-dir` to
save the comparison elsewhere. The summary compares image counts, complete pairs,
short revisit intervals, unused minutes, and slew angles for each night. The stats
step checks every field history in the replay database against the generated plans.

Add `--gifs-first 10` to make both sky-rotation and visit-count GIFs for the first ten replayed nights in `replay_results/gifs/`.

When the 30 nightly baseline plans from 2026-07-17 through 2026-08-15 are present
in `plans/`, run their simulation with:

```bash
python replay_plans.py --limit 30 --output-dir simulation_30_days --gifs-first 10
```

This writes 30 nightly plans through 2026-08-15, an isolated scheduling database,
nightly and field-level statistics, both 2D sky-grid maps, and GIFs for the first
ten nights under `simulation_30_days/`.

Compare the original and simulated plans over those same 30 nights:

```bash
python compare_survey_maps.py --days 30 --simulation-dir simulation_30_days
```

This writes each plan set's statistics and maps to `comparison_30_days/`, plus
`visits_side_by_side.html` and `visit_gaps_side_by_side.html`. Both panels share
the same border-color scale, and the gap slider controls both panels. The folder
also contains a field-by-field CSV comparison.

## Analyze a survey of plans

Analyze every CSV in `plans/`:

```bash
python analyze_plans.py
```

Analyze selected plans and write results elsewhere:

```bash
python analyze_plans.py plans/20260717.csv plans/20260718.csv \
  --output-dir survey_stats
```

To update the full 32-night clustered replay statistics and verify its field history:

```bash
python analyze_plans.py clustered_results/20*.csv \
  --db clustered_results/replay_grid.db --output-dir survey_stats
```

The analysis writes interactive HTML figures, `summary.txt`, `nightly_summary.csv`, and
`per_field_stats.csv`. The per-field table includes every grid field, even those with
zero visits, plus the number of nights when a complete pair was visible at sampled
10-minute intervals. With `--db`, analysis fails if any field's last scheduled MJD
disagrees with the plans. `visits_per_field_map.html` shows the full field grid on
a flat RA–Dec chart, with the Galactic plane overlaid. Only each field's border is
colored; unvisited fields have gray borders. The companion
`visit_gap_percentile_map.html` has a P0–P100 slider for the gap between consecutive
initial visits to each field; the 30-minute dither does not count as a separate
survey revisit. The maps use a lightweight, self-contained 2D canvas and open
without a network connection. Add `--show` to open the remaining Plotly figures
while generating them.

## Repository layout

```text
assets/              Field-grid CSV and SQLite database
plans/               Nightly observing-plan CSVs
survey_stats/        Generated survey statistics and HTML figures
build_field_grid.py  Generate the field grid and database
generate_obs_plan.py Generate a nightly observing plan
analyze_plans.py     Analyze multiple plans
visualizations.py    Shared Plotly and coverage visualizations
constants.py         Telescope, field-of-view, and scheduling parameters
```

Additional notebooks, FITS/region files, archived plans, and generated GIFs are included for exploration and validation.

## Configuration

Common scheduling parameters are defined in [`constants.py`](constants.py), including the field-of-view dimensions, exposure/readout times, declination limit, maximum airmass, and Moon-separation constraint. Changes to these values affect both field-grid generation and plan scheduling, so regenerate dependent outputs when appropriate.
