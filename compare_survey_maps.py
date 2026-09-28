"""Compare original and replayed plans over the same observing nights.

Example:
    python compare_survey_maps.py --days 30 --simulation-dir simulation_30_days
"""

import argparse
from pathlib import Path
import subprocess
import sys

import pandas as pd

from analyze_plans import load_plans
from sky_grid_maps import write_comparison_maps


def gaps_for(paths):
    plans = load_plans(paths)
    initials = plans[~plans['is_dither']]
    return {name: times.diff().dropna().dt.total_seconds().to_numpy() / 86400
            for name, times in initials.groupby('base_field')['utc'] if len(times) > 1}


def analyze(paths, output_dir, db=None):
    command = [sys.executable, 'analyze_plans.py', *map(str, paths),
               '--output-dir', str(output_dir)]
    if db is not None:
        command += ['--db', str(db)]
    subprocess.run(command, check=True)
    return pd.read_csv(output_dir / 'per_field_stats.csv')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-dir', type=Path, default=Path('plans'))
    parser.add_argument('--simulation-dir', type=Path, default=Path('simulation_30_days'))
    parser.add_argument('--days', type=int, default=30)
    parser.add_argument('--output-dir', type=Path, default=Path('comparison_30_days'))
    args = parser.parse_args()
    if args.days <= 0:
        parser.error('--days must be positive')

    baseline = {path.stem: path for path in args.baseline_dir.glob('20*.csv')}
    simulation = {path.stem: path for path in args.simulation_dir.glob('20*.csv')}
    nights = sorted(baseline.keys() & simulation.keys())[:args.days]
    if len(nights) != args.days:
        parser.error(f'Found {len(nights)} matching nights; requested {args.days}')
    old_paths = [baseline[night] for night in nights]
    new_paths = [simulation[night] for night in nights]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    old_stats = analyze(old_paths, args.output_dir / 'baseline')
    db = args.simulation_dir / 'replay_grid.db'
    new_stats = analyze(new_paths, args.output_dir / 'simulation', db if db.exists() else None)
    old_gaps, new_gaps = gaps_for(old_paths), gaps_for(new_paths)
    write_comparison_maps(old_stats, old_gaps, new_stats, new_gaps, args.output_dir)

    columns = ['base_field', 'visits', 'complete_pairs', 'observed_nights',
               'median_gap_days', 'first_obs', 'last_obs']
    per_field = old_stats[columns].merge(new_stats[columns], on='base_field',
                                         suffixes=('_original', '_simulation'),
                                         validate='one_to_one')
    per_field['visit_difference'] = (per_field['visits_simulation'] -
                                     per_field['visits_original'])
    per_field.to_csv(args.output_dir / 'field_comparison.csv', index=False)

    old_nightly = pd.read_csv(args.output_dir / 'baseline/nightly_summary.csv')
    new_nightly = pd.read_csv(args.output_dir / 'simulation/nightly_summary.csv')
    rows = [
        ('nights', len(nights), len(nights)),
        ('science_images', int(old_stats['visits'].sum()), int(new_stats['visits'].sum())),
        ('complete_pairs', int(old_stats['complete_pairs'].sum()),
         int(new_stats['complete_pairs'].sum())),
        ('unique_fields', int((old_stats['visits'] > 0).sum()),
         int((new_stats['visits'] > 0).sum())),
        ('unpaired_images', int(old_stats['visits'].sum() - 2 * old_stats['complete_pairs'].sum()),
         int(new_stats['visits'].sum() - 2 * new_stats['complete_pairs'].sum())),
        ('unused_hours', old_nightly['unused_minutes'].sum() / 60,
         new_nightly['unused_minutes'].sum() / 60),
    ]
    summary = pd.DataFrame(rows, columns=['metric', 'original', 'simulation'])
    summary.to_csv(args.output_dir / 'comparison_summary.csv', index=False)
    report = f"""# {len(nights)}-night plan comparison

Original plans and clustered replay, {nights[0]} through {nights[-1]}.

| Metric | Original | Clustered replay |
| --- | ---: | ---: |
| Science images | {int(old_stats['visits'].sum()):,} | {int(new_stats['visits'].sum()):,} |
| Complete pairs | {int(old_stats['complete_pairs'].sum()):,} | {int(new_stats['complete_pairs'].sum()):,} |
| Unique fields visited | {int((old_stats['visits'] > 0).sum())} | {int((new_stats['visits'] > 0).sum())} |
| Unpaired images | {int(old_stats['visits'].sum() - 2 * old_stats['complete_pairs'].sum())} | {int(new_stats['visits'].sum() - 2 * new_stats['complete_pairs'].sum())} |
| Unused hours | {old_nightly['unused_minutes'].sum() / 60:.2f} | {new_nightly['unused_minutes'].sum() / 60:.2f} |

Open `visits_side_by_side.html` or `visit_gaps_side_by_side.html` to compare the
full grids with one color scale. The gap slider controls both panels; gaps are
measured between consecutive initial visits, excluding dithers. Gray borders
mark fields with no visits or no measured inter-night gap. The two individual
sets of maps and statistics are in `baseline/` and `simulation/`.

`field_comparison.csv` gives per-field counts and cadence statistics for both
plans. Visibility diagnostics appear in both `summary.txt` files; the simulated
database-history audit appears in `simulation/summary.txt`.
"""
    (args.output_dir / 'REPORT.md').write_text(report)
    print(f'Compared {len(nights)} nights, {nights[0]} through {nights[-1]}')
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
