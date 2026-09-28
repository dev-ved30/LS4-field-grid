"""Survey statistics across a collection of LS4 observing-plan CSVs.

Aggregates one or more observing plans (see plans/*.csv for the schema) and
produces per-field visit counts, cadence distributions, and missing/stale
field diagnostics. Figures are written as HTML files (one per plot) and a
text summary is printed to stdout.

Usage:
    python analyze_plans.py                # all plans in plans/
    python analyze_plans.py plan1.csv plan2.csv
    python analyze_plans.py --show         # also open the figures interactively
"""

import argparse
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils import iers

from constants import (LS4_field_grid_path, default_block_duration, exp,
                       half_field_offset_ra, min_declination, read_out)
from generate_obs_plan import _visibility_at_slots
from sky_grid_maps import write_sky_grid_maps

iers.conf.auto_download = False
iers.conf.auto_max_age = None


def argument_parser():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('plans', nargs='*', type=Path,
                        help="Observing-plan CSV files. Defaults to all plans/*.csv.")
    parser.add_argument('--limit', type=int,
                        help='Analyze only the first N plans after sorting by filename.')
    parser.add_argument('--output-dir', type=Path, default=Path('survey_stats'),
                        help="Directory for output figures (default: survey_stats/).")
    parser.add_argument('--grid', type=Path, default=Path(LS4_field_grid_path),
                        help='Field-grid CSV used to include fields with zero visits.')
    parser.add_argument('--db', type=Path,
                        help='Optional replay database; verify every field history against the plans.')
    parser.add_argument('--visibility-sample-minutes', type=int, default=10,
                        help='Spacing for sampled pair visibility checks (default: 10).')
    parser.add_argument('--show', action='store_true',
                        help="Open figures in the browser in addition to saving them.")
    return parser.parse_args()


def load_plans(paths):
    """Load all plan CSVs into one dataframe with a night column."""

    frames = []
    for path in paths:
        df = pd.read_csv(path)
        df['night'] = path.stem
        frames.append(df)

    plans = pd.concat(frames, ignore_index=True)

    # drop unused/slew rows, keep only real pointings
    plans = plans[~plans['target'].isin(['Unused Time', 'TransitionBlock'])].copy()
    plans = plans.dropna(subset=['ra', 'dec'])

    plans['base_field'] = plans['target'].str.removesuffix('_dither')
    plans['is_dither'] = plans['target'].str.endswith('_dither')
    plans['utc'] = pd.to_datetime(plans['start time (UTC)'])

    return plans.sort_values('utc').reset_index(drop=True)


def load_grid(path):
    """Read the complete field universe, including unscheduled fields."""
    grid = pd.read_csv(path).rename(columns={
        'Field Name': 'base_field', 'ra_deg': 'ra', 'dec_deg': 'dec'})
    if grid['base_field'].isna().any() or grid['base_field'].duplicated().any():
        raise ValueError(f'Grid field names must be unique and nonempty: {path}')
    return grid


def audit_history(per_field, db_path):
    """Check that each replay DB timestamp equals its last planned pointing."""
    with sqlite3.connect(db_path) as db:
        history = pd.read_sql_query(
            'SELECT "Field Name" AS base_field, last_scheduled_mjd FROM grid', db)
    if history['base_field'].duplicated().any():
        raise ValueError(f'Duplicate field names in {db_path}')
    checked = per_field.merge(history, on='base_field', how='outer', indicator=True,
                              validate='one_to_one')
    expected = checked['last_plan_mjd'].fillna(0).to_numpy(dtype=float)
    actual = checked['last_scheduled_mjd'].fillna(0).to_numpy(dtype=float)
    checked['history_matches'] = (checked['_merge'] == 'both') & np.isclose(
        expected, actual, rtol=0, atol=1e-6)
    mismatches = checked.loc[~checked['history_matches'], 'base_field'].tolist()
    if mismatches:
        raise ValueError(f'{len(mismatches)} replay DB histories disagree with plans; '
                         f'first fields: {mismatches[:10]}')
    return checked.drop(columns='_merge')


def sampled_visibility(grid, plan_paths, sample_minutes):
    """Count nights when a field has at least one sampled initial/dither opportunity."""
    if sample_minutes <= 0:
        raise ValueError('Visibility sample interval must be positive')
    eligible = grid['dec'].to_numpy() > min_declination.to_value(u.deg)
    visible_nights = np.zeros(len(grid), dtype=int)
    positions = np.flatnonzero(eligible)
    ra = grid.iloc[positions]['ra'].to_numpy()
    dec = grid.iloc[positions]['dec'].to_numpy()
    initial = SkyCoord(ra * u.deg, dec * u.deg)
    dither = SkyCoord((ra - half_field_offset_ra.to_value(u.deg)) * u.deg,
                      dec * u.deg)
    duration = exp + read_out
    step = sample_minutes * u.minute
    for path in plan_paths:
        plan = pd.read_csv(path)
        start = Time(plan.iloc[0]['start time (UTC)'], scale='utc')
        end = Time(plan.iloc[-1]['end time (UTC)'], scale='utc')
        available = end - start - default_block_duration - duration
        count = int(np.floor((available / step).decompose().value)) + 1
        if count <= 0:
            continue
        starts = start + np.arange(count) * step
        possible = (_visibility_at_slots(initial, starts, duration) &
                    _visibility_at_slots(dither, starts + default_block_duration, duration))
        visible_nights[positions] += possible.any(axis=1)
    return eligible, visible_nights


def save(fig, name, output_dir, show):
    path = output_dir / f"{name}.html"
    fig.write_html(path, include_plotlyjs='cdn')
    print(f"  wrote {path}")
    if show:
        fig.show()


def plot_visit_histogram(per_field, output_dir, show):
    """Histogram of visit counts across fields."""

    counts = per_field['visits']
    fig = go.Figure(go.Histogram(
        x=counts,
        nbinsx=int(counts.max()) + 1,
        marker_color='#4169ab',
    ))
    n_never = int((counts == 0).sum())
    fig.update_layout(
        title=f'Distribution of visits per field ({n_never} fields never visited)',
        xaxis_title='Visits',
        yaxis_title='Number of fields',
        bargap=0.05,
    )
    save(fig, 'visits_per_field_hist', output_dir, show)


def plot_cadence(cadence, per_field_cadence, output_dir, show):
    """Distribution of gaps between initial pointings to the same field."""

    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=cadence['gap_days'],
        nbinsx=60,
        marker_color='#8a63c4',
        name='initial-pointing gaps',
    ))
    fig.update_layout(
        title='Cadence: time between initial pointings to the same field',
        xaxis_title='Gap (days)',
        yaxis_title='Count (log)',
        yaxis_type='log',
        bargap=0.02,
    )
    save(fig, 'cadence_gaps_hist', output_dir, show)

    fig = go.Figure(go.Histogram(
        x=per_field_cadence['median_gap_days'],
        nbinsx=40,
        marker_color='#c48a63',
    ))
    fig.update_layout(
        title='Per-field median gap between initial pointings',
        xaxis_title='Median gap (days)',
        yaxis_title='Number of fields',
    )
    save(fig, 'median_cadence_per_field_hist', output_dir, show)


def plot_cadence_cdf(cadence, output_dir, show):
    """CDF of gaps between consecutive visits to the same field."""

    gaps = np.sort(cadence['gap_days'].to_numpy())
    percent = 100.0 * np.arange(1, len(gaps) + 1) / len(gaps)

    fig = go.Figure(go.Scatter(
        x=gaps,
        y=percent,
        mode='lines',
        line=dict(color='#4169ab', width=2),
        hovertemplate='%{x:.2f} d: %{y:.1f}% of gaps<extra></extra>',
    ))

    median = np.median(gaps)
    p90 = np.percentile(gaps, 90)
    for value, label in ((median, 'median'), (p90, '90th pct')):
        fig.add_hline(y=50 if label == 'median' else 90, line=dict(width=1, dash='dot',
                      color='#8a63c4'), opacity=0.5)
        fig.add_vline(x=value, line=dict(width=1, dash='dot', color='#8a63c4'),
                      annotation_text=f'{label}: {value:.2f} d', annotation_position='top',
                      opacity=0.5)

    fig.update_layout(
        title=f'Cadence CDF: initial-pointing gaps shorter than x ({len(gaps)} gaps)',
        xaxis_title='Gap between consecutive visits (days)',
        yaxis_title='% of gaps ≤ x',
        xaxis_type='log',
        height=600,
    )
    save(fig, 'cadence_gaps_cdf', output_dir, show)


def plot_nightly_summary(nightly, output_dir, show):
    """Images taken and time efficiency for each night."""

    nightly = nightly.sort_values('night').copy()
    nightly['time_efficiency_pct'] = 100 * nightly['used_minutes'] / nightly['night_minutes']

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=nightly['night'],
        y=nightly['images'],
        name='images',
        marker_color='#44b8e8',
    ))
    fig.add_trace(go.Scatter(
        x=nightly['night'],
        y=nightly['time_efficiency_pct'],
        name='time efficiency (%)',
        yaxis='y2',
        mode='lines+markers',
        line=dict(color='#e652a0'),
    ))
    fig.update_layout(
        title='Nightly summary',
        xaxis_title='Plan',
        yaxis_title='Images',
        yaxis2=dict(title='Time efficiency (%)', overlaying='y', side='right',
                    range=[0, 105]),
        legend=dict(x=0.01, y=0.99),
    )
    save(fig, 'nightly_summary', output_dir, show)


def main():

    args = argument_parser()

    plan_paths = sorted(args.plans) if args.plans else sorted(Path('plans').glob('*.csv'))
    if args.limit is not None:
        if args.limit <= 0 or len(plan_paths) < args.limit:
            raise SystemExit(f'--limit requires 1 to {len(plan_paths)} plans')
        plan_paths = plan_paths[:args.limit]
    if not plan_paths:
        raise SystemExit('No plan CSVs found.')

    print(f'Analyzing {len(plan_paths)} plans: {", ".join(p.name for p in plan_paths)}\n')

    plans = load_plans(plan_paths)
    grid = load_grid(args.grid)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ---- per-field visit counts ---------------------------------------
    night_fields = plans.groupby(['night', 'base_field']).agg(
        images=('target', 'size'), dithers=('is_dither', 'sum')).reset_index()
    complete = night_fields[(night_fields['images'] == 2) &
                            (night_fields['dithers'] == 1)]
    pair_counts = complete.groupby('base_field').size().rename('complete_pairs')
    per_field = plans.groupby('base_field').agg(
        visits=('target', 'count'),
        dithers=('is_dither', 'sum'),
        observed_nights=('night', 'nunique'),
        first_obs=('utc', 'min'),
        last_obs=('utc', 'max'),
        last_plan_mjd=('start_time_mjd', 'max'),
    ).reset_index()
    unknown = per_field.loc[~per_field['base_field'].isin(grid['base_field']), 'base_field']
    if len(unknown):
        raise ValueError(f'Plans contain fields absent from grid: {unknown.tolist()[:10]}')
    per_field = grid[['base_field', 'program_id', 'ra', 'dec']].merge(
        per_field, on='base_field', how='left', validate='one_to_one')
    for column in ('visits', 'dithers', 'observed_nights'):
        per_field[column] = per_field[column].fillna(0).astype(int)
    per_field['complete_pairs'] = per_field['base_field'].map(pair_counts).fillna(0).astype(int)
    per_field['above_declination_limit'], per_field['sampled_visible_nights'] = (
        sampled_visibility(per_field, plan_paths, args.visibility_sample_minutes))
    if args.db:
        per_field = audit_history(per_field, args.db)

    observed = per_field[per_field['visits'] > 0]

    # ---- cadence between separate nights (the dither is not a revisit) ---
    gap_rows = []
    gap_distributions = {}
    initials = plans[~plans['is_dither']]
    for name, times in initials.sort_values('utc').groupby('base_field')['utc']:
        gaps = times.diff().dropna().dt.total_seconds() / 86400.0
        if len(gaps):
            gap_distributions[name] = gaps.to_numpy()
        for gap in gaps:
            gap_rows.append({'base_field': name, 'gap_days': gap})
    cadence = pd.DataFrame(gap_rows)

    if len(cadence):
        per_field_cadence = cadence.groupby('base_field')['gap_days'].agg(
            median_gap_days='median',
            mean_gap_days='mean',
            min_gap_days='min',
            max_gap_days='max',
        ).reset_index()
    else:
        per_field_cadence = pd.DataFrame(columns=[
            'base_field', 'median_gap_days', 'mean_gap_days', 'min_gap_days', 'max_gap_days'])

    # ---- nightly summary -------------------------------------------------
    all_rows = []
    for path in plan_paths:
        df = pd.read_csv(path)
        used = df[~df['target'].isin(['Unused Time', 'TransitionBlock'])]
        nightly = {
            'night': path.stem,
            'images': len(used.dropna(subset=['ra'])),
            'used_minutes': used['duration (minutes)'].clip(lower=0).sum(),
            'unused_minutes': df[df['target'] == 'Unused Time']['duration (minutes)'].sum(),
        }
        nightly['night_minutes'] = nightly['used_minutes'] + nightly['unused_minutes']
        all_rows.append(nightly)
    nightly = pd.DataFrame(all_rows)

    # ---- text summary ----------------------------------------------------
    span_start = plans['utc'].min()
    span_end = plans['utc'].max()
    span_days = (span_end - span_start).total_seconds() / 86400.0

    never = per_field[per_field['visits'] == 0]
    visible_never = never[never['sampled_visible_nights'] > 0]
    summary_lines = [
        '================ survey summary ================',
        f'nights               : {len(plan_paths)}',
        f'span                 : {span_start:%Y-%m-%d} to {span_end:%Y-%m-%d} ({span_days:.1f} days)',
        f'grid fields          : {len(per_field)}',
        f'total visits         : {len(plans)} ({plans["is_dither"].sum()} dithered)',
        f'complete pairs       : {len(complete)}',
        f'unpaired images      : {len(plans) - 2 * len(complete)}',
        f'unique fields visited: {len(observed)}',
        f'fields never visited : {len(never)} '
        f'(by program: {never["program_id"].value_counts().to_dict()})',
        f'  below dec limit    : {(~never["above_declination_limit"]).sum()}',
        f'  pair visible sampled, unvisited: {len(visible_never)}',
        f'  no sampled pair opportunity: '
        f'{len(never) - len(visible_never) - (~never["above_declination_limit"]).sum()}',
        f'mean visits/field    : {observed["visits"].mean():.2f}',
    ]
    if len(cadence):
        summary_lines.append(f'initial-visit cadence: median {cadence["gap_days"].median():.2f} d, '
                             f'mean {cadence["gap_days"].mean():.2f} d, '
                             f'max {cadence["gap_days"].max():.2f} d')
    summary_lines.append(f'total time used      : {nightly["used_minutes"].sum()/60:.1f} h of '
                         f'{nightly["night_minutes"].sum()/60:.1f} h '
                         f'({100*nightly["used_minutes"].sum()/nightly["night_minutes"].sum():.1f}%)')
    if args.db:
        summary_lines.append(f'replay DB histories  : {len(per_field)} of {len(per_field)} match plans')

    stale = observed.copy()
    stale['days_since_last'] = (span_end - stale['last_obs']).dt.total_seconds() / 86400.0

    summary_lines.append('longest time since last observation (top 10):')
    top = stale.nlargest(10, 'days_since_last')[['base_field', 'visits', 'last_obs', 'days_since_last']]
    summary_lines.append(top.to_string(index=False))
    summary_lines.append(f'Visibility sampled every {args.visibility_sample_minutes} minutes; '
                         'an opportunity requires both pointings to pass all sky constraints '
                         '30 minutes apart, without accounting for other scheduled fields.')
    summary = '\n'.join(summary_lines) + '\n'
    print(summary)
    (args.output_dir / 'summary.txt').write_text(summary)

    # ---- figures -----------------------------------------------------------
    print('\nWriting figures...')
    write_sky_grid_maps(per_field, gap_distributions, args.output_dir)
    plot_visit_histogram(per_field, args.output_dir, args.show)
    if len(cadence):
        plot_cadence(cadence, per_field_cadence, args.output_dir, args.show)
        plot_cadence_cdf(cadence, args.output_dir, args.show)
    plot_nightly_summary(nightly, args.output_dir, args.show)

    # per-field table for further digging
    table_path = args.output_dir / 'per_field_stats.csv'
    export = per_field.merge(per_field_cadence, on='base_field', how='left')
    export.to_csv(table_path, index=False)
    print(f'  wrote {table_path}')
    nightly.to_csv(args.output_dir / 'nightly_summary.csv', index=False)


if __name__ == '__main__':
    main()
