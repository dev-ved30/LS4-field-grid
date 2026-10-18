"""Generate consecutive LS4 observing nights from astronomical twilight.

Example:
    python simulate_nights.py --start-date 2026-07-17 --nights 90 \
        --output-dir simulation_90_days
"""

import argparse
from datetime import date, timedelta
from pathlib import Path
import sqlite3
import subprocess
import sys

import pandas as pd
from astropy.time import Time
from astropy.utils import iers

from constants import LS4, LS4_field_grid_path, read_out, slew_rate
from generate_obs_plan import get_obs_plan
from replay_plans import metrics

iers.conf.auto_download = False
iers.conf.auto_max_age = None


def night_window(day):
    """Return the astronomical-twilight window following local evening."""
    noon_utc = Time(f'{day.isoformat()}T12:00:00', scale='utc')
    start = LS4.twilight_evening_astronomical(noon_utc, which='next')
    end = LS4.twilight_morning_astronomical(start, which='next')
    if end <= start:
        raise ValueError(f'Invalid twilight window for {day}')
    return start, end


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start-date', type=date.fromisoformat, default=date(2026, 7, 17))
    parser.add_argument('--nights', type=int, default=90)
    parser.add_argument('--output-dir', type=Path, default=Path('simulation_90_days'))
    args = parser.parse_args()
    if args.nights <= 0:
        parser.error('--nights must be positive')

    args.output_dir.mkdir(parents=True, exist_ok=True)
    db_path = args.output_dir / 'simulation_grid.db'
    grid = pd.read_csv(LS4_field_grid_path)
    grid['last_scheduled_mjd'] = 0.0
    with sqlite3.connect(db_path) as db:
        grid.to_sql('grid', db, if_exists='replace', index=False)

    rows = []
    paths = []
    for offset in range(args.nights):
        day = args.start_date + timedelta(days=offset)
        start, end = night_window(day)
        path = args.output_dir / f'{day:%Y%m%d}.csv'
        print(f'\nNight {offset + 1}/{args.nights}: {day} '
              f'{start.iso} to {end.iso}', flush=True)
        plan = get_obs_plan(start, end, path, db_path=db_path, show_plots=False)
        result = metrics(plan)
        if result['unpaired_images'] or result['short_pairs']:
            raise AssertionError(f'Pair or revisit violation in {path}')
        rows.append({'night': path.stem, 'start_utc': start.iso, 'end_utc': end.iso,
                     'night_hours': (end - start).to_value('hour'), **result})
        paths.append(path)
        pd.DataFrame(rows).to_csv(args.output_dir / 'summary.csv', index=False)

    stats_dir = args.output_dir / 'stats'
    subprocess.run([sys.executable, 'analyze_plans.py', *map(str, paths),
                    '--db', str(db_path), '--output-dir', str(stats_dir)], check=True)
    summary = pd.DataFrame(rows)
    fields = pd.read_csv(stats_dir / 'per_field_stats.csv')
    report = f"""# {args.nights}-night LS4 simulation

Consecutive nights from {args.start_date} through
{args.start_date + timedelta(days=args.nights - 1)}. Each observing window runs
from astronomical evening twilight to astronomical morning twilight at La Silla.
The scheduler started with a fresh grid history and used a configurable
{slew_rate.to_value('deg/s'):.2f}°/second slew rate with {read_out.to_value('s'):.0f}-second readout.

| Metric | Result |
| --- | ---: |
| Science images | {int(summary['images'].sum()):,} |
| Complete initial/dither pairs | {int(summary['complete_pairs'].sum()):,} |
| Unpaired images | {int(summary['unpaired_images'].sum())} |
| Revisit pairs under 30 minutes | {int(summary['short_pairs'].sum())} |
| Unused hours | {summary['unused_minutes'].sum() / 60:.2f} |
| Unique fields visited | {int((fields['visits'] > 0).sum()):,} of {len(fields):,} |
| Replay DB histories matching plans | {int(fields['history_matches'].sum()):,} of {len(fields):,} |

The `stats/` directory contains field-level statistics and standalone border-only
visits and gap-percentile sky maps. `summary.csv` contains the nightly metrics.
The original plans and project grid database were not modified.
"""
    (args.output_dir / 'REPORT.md').write_text(report)
    print(f'Completed {args.nights} nights; {summary["images"].sum():,} images in '
          f'{summary["complete_pairs"].sum():,} complete pairs')


if __name__ == '__main__':
    main()
