"""Replay existing plan nights from a clean grid and compare schedule metrics.

The original plans and assets database are read only. New plans, a temporary
SQLite scheduling history, and a per-night comparison are written to --output-dir.
"""

import argparse
from pathlib import Path
import sqlite3
import subprocess
import sys

import numpy as np
import pandas as pd
from astropy.time import Time
from astropy.utils import iers

iers.conf.auto_download = False
iers.conf.auto_max_age = None

from constants import LS4_field_grid_path, read_out, slew_rate
from generate_obs_plan import get_obs_plan


def metrics(plan):
    science = plan[~plan['target'].isin(['Unused Time', 'TransitionBlock'])].copy()
    science['base_field'] = science['target'].str.removesuffix('_dither')
    science['is_dither'] = science['target'].str.endswith('_dither')
    pairs = science.groupby('base_field').agg(visits=('target', 'size'),
                                               dithers=('is_dither', 'sum'))
    complete = pairs[(pairs['visits'] == 2) & (pairs['dithers'] == 1)]
    short = 0
    for name in complete.index:
        visits = science[science['base_field'] == name].sort_values('start_time_mjd')
        if (visits.iloc[1]['start_time_mjd'] - visits.iloc[0]['start_time_mjd']) * 1440 < 30 - 1e-5:
            short += 1
    ra = np.deg2rad(science['ra'].to_numpy())
    dec = np.deg2rad(science['dec'].to_numpy())
    if len(science) > 1:
        cosine = (np.sin(dec[:-1]) * np.sin(dec[1:]) +
                  np.cos(dec[:-1]) * np.cos(dec[1:]) * np.cos(ra[:-1] - ra[1:]))
        steps = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    else:
        steps = []
    slew_speed = slew_rate.to_value('deg/s')
    readout_seconds = read_out.to_value('s')
    return {
        'images': len(science),
        'complete_pairs': len(complete),
        'unpaired_images': len(science) - 2 * len(complete),
        'short_pairs': short,
        'unused_minutes': plan.loc[plan['target'] == 'Unused Time', 'duration (minutes)'].sum(),
        'median_slew_deg': float(np.median(steps)) if len(steps) else 0,
        'p90_slew_deg': float(np.percentile(steps, 90)) if len(steps) else 0,
        'max_slew_deg': float(np.max(steps)) if len(steps) else 0,
        'slews_over_readout': int(np.sum(np.asarray(steps) / slew_speed > readout_seconds)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path, default=Path('plans'))
    parser.add_argument('--output-dir', type=Path, default=Path('replay_results'))
    parser.add_argument('--limit', type=int, help='Only replay the first N nights (for a quick check)')
    parser.add_argument('--gifs-first', type=int, default=0,
                        help='Create observing-plan and visit-count GIFs for the first N replayed nights')
    args = parser.parse_args()
    paths = sorted(args.baseline_dir.glob('*.csv'))
    if args.limit:
        paths = paths[:args.limit]
    if not paths:
        parser.error(f'No CSV plans in {args.baseline_dir}')

    args.output_dir.mkdir(parents=True, exist_ok=True)
    db_path = args.output_dir / 'replay_grid.db'
    grid = pd.read_csv(LS4_field_grid_path)
    grid['last_scheduled_mjd'] = 0.0
    with sqlite3.connect(db_path) as db:
        grid.to_sql('grid', db, if_exists='replace', index=False)

    results = []
    for path in paths:
        original = pd.read_csv(path)
        start = Time(original.iloc[0]['start time (UTC)'], scale='utc')
        end = Time(original.iloc[-1]['end time (UTC)'], scale='utc')
        new_path = args.output_dir / path.name
        print(f'\nReplaying {path.name}: {start.iso} to {end.iso}', flush=True)
        updated = get_obs_plan(start, end, new_path, db_path=db_path, show_plots=False)
        baseline = metrics(original)
        new = metrics(updated)
        if new['unpaired_images'] or new['short_pairs']:
            raise AssertionError(f'Pair or cadence violation in {new_path}')
        results.append({'night': path.stem, **{f'baseline_{key}': value for key, value in baseline.items()},
                        **{f'new_{key}': value for key, value in new.items()}})
        pd.DataFrame(results).to_csv(args.output_dir / 'summary.csv', index=False)

    summary = pd.DataFrame(results)
    print('\nTotals across replayed nights:')
    for key in ('images', 'complete_pairs', 'unpaired_images', 'short_pairs', 'unused_minutes'):
        print(f'  {key}: {summary[f"baseline_{key}"].sum():.1f} -> {summary[f"new_{key}"].sum():.1f}')
    print(f'Wrote {args.output_dir / "summary.csv"}')

    stats_dir = args.output_dir / 'stats'
    subprocess.run([sys.executable, 'analyze_plans.py',
                    *[str(args.output_dir / path.name) for path in paths],
                    '--db', str(db_path), '--output-dir', str(stats_dir)], check=True)

    if args.gifs_first:
        gif_dir = args.output_dir / 'gifs'
        gif_dir.mkdir(parents=True, exist_ok=True)
        for path in paths[:args.gifs_first]:
            plan = args.output_dir / path.name
            for script, suffix in (('make_observing_plan_gif.py', ''),
                                   ('make_visit_count_gif.py', '_visits')):
                output = gif_dir / f'{path.stem}{suffix}.gif'
                subprocess.run([sys.executable, script, str(plan), '--output', str(output)],
                               check=True)


if __name__ == '__main__':
    main()
