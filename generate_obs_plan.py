import time
import argparse
import sqlite3
from pathlib import Path

import pandas as pd
import healpy as hp
import numpy as np

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from constants import *
from visualizations import plot_coverage_map

import warnings
from astropy.coordinates import NonRotationTransformationWarning

warnings.filterwarnings(
    "ignore",
    category=NonRotationTransformationWarning
)

# Get the current time
current_time = Time.now() # Fix this to a specific time for testing, e.g. Time("2024-06-01 00:00:00")

def argument_parser():

    parser = argparse.ArgumentParser(description="Generate an observing plan for LS4 based on the field grid and current visibility.")
    parser.add_argument('--mjd', required=False, default=None, help="MJD for which to generate the observing plan. If not provided, the current time will be used.")
    parser.add_argument("--output", required=False, type=str, default=None, help="Path to save the generated observing plan CSV file. Default will just save to plans/yyyymmdd.csv")
    parser.add_argument("--db", type=Path, default=Path(LS4_field_grid_db_path),
                        help="Field-grid database (default: assets/LS4_field_grid.db)")
    parser.add_argument("--no-plots", action="store_true", help="Skip the coverage plot and area calculation")
    return parser.parse_args()

def compute_theoretical_max_images_per_night(night_start, night_end):

    night_duration = (night_end - night_start).to_value(u.second) * u.second
    total_time_per_field = exp + read_out
    theoretical_max_fields = (night_duration / total_time_per_field).to_value()

    return theoretical_max_fields

def compute_time_efficiency(obs_plan, night_start, night_end):

    total_night_duration = (night_end - night_start).to_value(u.minute) * u.minute

    # unused time
    unused_time = obs_plan[obs_plan['target'] == 'Unused Time']
    total_unused_time = unused_time['duration (minutes)'].sum() * u.minute

    efficiency = (total_night_duration - total_unused_time) / total_night_duration

    print(f"Total night duration: {total_night_duration.to_value(u.hour):.2f} hours")
    print(f"Total unused time: {total_unused_time.to_value(u.hour):.2f} hours")
    print(f"Time efficiency: {efficiency:.2%}")

    return efficiency

def compute_theoretical_max_area_per_night(night_start, night_end):

    theoretical_max_fields = compute_theoretical_max_images_per_night(night_start, night_end)
    area_per_field = FOV_width * FOV_length
    theoretical_max_area = theoretical_max_fields * area_per_field
    return theoretical_max_area.to_value(u.deg**2) / 2 # divide by 2 to account for the dithered pointings that cover some of the same area

def compute_union_area(obs_plan, nside=2048, night_start=None, night_end=None):
    """
    Compute the union area (deg^2) of rectangular FoVs on the sky.

    Parameters
    ----------
    pointings_ra : array-like
        RA of pointings (degrees)
    pointings_dec : array-like
        Dec of pointings (degrees)
    fov_width_deg : float
        Width of FoV (degrees, along RA direction in tangent plane)
    fov_height_deg : float
        Height of FoV (degrees, along Dec direction in tangent plane)
    nside : int
        HEALPix resolution (higher = more accurate, slower)

    Returns
    -------
    area_deg2 : float
        Total union area in square degrees
    """
    df = obs_plan[(obs_plan['target'] != 'Unused Time') & (obs_plan['target'] != 'TransitionBlock')]
    pointings_ra, pointings_dec = df['ra'].to_numpy(), df['dec'].to_numpy()

    npix = hp.nside2npix(nside)
    observed = np.zeros(npix, dtype=bool)
    visits = np.zeros(npix, dtype=int)

    half_w =  FOV_width.to(u.deg).value/ 2
    half_h = FOV_length.to(u.deg).value/ 2

    for ra, dec in zip(pointings_ra, pointings_dec):
        center = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')

        # Define rectangle corners in tangent plane
        offsets = [
            (-half_w, -half_h),
            ( half_w, -half_h),
            ( half_w,  half_h),
            (-half_w,  half_h),
        ]

        corners = []
        for dx, dy in offsets:
            corner = center.spherical_offsets_by(dx*u.deg, dy*u.deg)
            corners.append(corner)


        # Convert to HEALPix vectors
        vecs = np.array([
            hp.ang2vec(np.radians(90 - c.dec.deg),
                       np.radians(c.ra.deg))
            for c in corners
        ])

        # Fill polygon
        pixels = hp.query_polygon(nside, vecs)
        observed[pixels] = True
        visits[pixels] += 1

    # Compute area
    pixel_area_sr = hp.nside2pixarea(nside)
    total_area_sr = observed.sum() * pixel_area_sr
    total_area_deg2 = total_area_sr * (180/np.pi)**2

    # maximum are possible for the night
    print(f"Total observed area: {total_area_deg2:.2f} deg^2")
    if night_start is not None and night_end is not None:
        max_area_deg2 = compute_theoretical_max_area_per_night(night_start, night_end)
        print(f"Maximum area possible for the night based on exposure time and readout time: {max_area_deg2:.2f} deg^2")
        print(f"Area efficiency: {100 * total_area_deg2 / max_area_deg2:.2f}%")

    # plot the ares covered by different number of visits to visualize the dither pattern and coverage
    visit_counts = np.bincount(visits)
    for i in range(1, min(10, len(visit_counts))):
        area_sr = (visits == i).sum() * pixel_area_sr
        area_deg2 = area_sr * (180/np.pi)**2
        print(f"Area covered by {i} visits: {area_deg2:.2f} deg^2")

    return visits

def load_field_grid(db_path):
    """Load candidates without limiting the pool before individual blocks are planned."""
    with sqlite3.connect(db_path) as db:
        grid = pd.read_sql_query('SELECT * FROM grid WHERE dec_deg > ?', db,
                                 params=(min_declination.to_value(u.deg),))
    if grid.empty:
        raise ValueError(f"No schedulable fields in {db_path}")
    return grid


def _visibility_at_slots(coords, starts, exposure_duration):
    """Check every constraint at the start and end of each proposed observation."""
    endpoints = np.column_stack((starts.mjd, (starts + exposure_duration).mjd)).ravel()
    times = Time(endpoints, format='mjd', scale='utc')
    allowed = np.logical_and.reduce([
        constraint(LS4, coords, times=times, grid_times_targets=True)
        for constraint in global_constraints
    ])
    return allowed[:, 0::2] & allowed[:, 1::2]


def _row(target, start, end, ra=None, dec=None, dither=None, block_number=0):
    return {
        'target': target,
        'start time (UTC)': start.to_datetime().isoformat(sep=' '),
        'end time (UTC)': end.to_datetime().isoformat(sep=' '),
        'duration (minutes)': (end - start).to_value(u.minute),
        'ra': ra,
        'dec': dec,
        'configuration': str({'dither': dither}) if dither is not None else '',
        'start_time_mjd': start.mjd,
        'block_number': block_number,
        'ra_hr': ra / 15 if ra is not None else None,
    }


def _block_rows(observations, block_start, block_end, block_number):
    """Emit science pointings and the actual gaps between them."""
    rows = []
    cursor = block_start
    for target, start, end, ra, dec, dither in sorted(observations, key=lambda item: item[1].mjd):
        if (start - cursor).to_value(u.second) > 1e-3:
            rows.append(_row('Unused Time', cursor, start, block_number=block_number))
        rows.append(_row(target, start, end, ra, dec, dither, block_number))
        cursor = end
    if (block_end - cursor).to_value(u.second) > 1e-3:
        rows.append(_row('Unused Time', cursor, block_end, block_number=block_number))
    return rows


def update_last_scheduled_mjd(obs_plan, db_path=LS4_field_grid_db_path):
    """Advance history only for fields with both scheduled pointings."""
    science = obs_plan[~obs_plan['target'].isin(['Unused Time', 'TransitionBlock'])].copy()
    science['base_field'] = science['target'].str.removesuffix('_dither')
    science['is_dither'] = science['target'].str.endswith('_dither')
    pairs = science.groupby('base_field').agg(
        visits=('target', 'size'), dithers=('is_dither', 'sum'),
        last_mjd=('start_time_mjd', 'max'))
    pairs = pairs[(pairs['visits'] == 2) & (pairs['dithers'] == 1)]
    with sqlite3.connect(db_path) as db:
        db.executemany('UPDATE grid SET last_scheduled_mjd = ? WHERE "Field Name" = ?',
                       [(float(row.last_mjd), name) for name, row in pairs.iterrows()])
    return len(pairs)


def validate_pairs(obs_plan, min_separation=default_block_duration):
    """Reject missing, duplicate, reversed, or too-close visits before saving."""
    starts = pd.to_datetime(obs_plan['start time (UTC)'])
    ends = pd.to_datetime(obs_plan['end time (UTC)'])
    if (ends <= starts).any() or (starts.iloc[1:].to_numpy() < ends.iloc[:-1].to_numpy()).any():
        raise ValueError('Observation rows overlap or have nonpositive duration')
    science = obs_plan[~obs_plan['target'].isin(['Unused Time', 'TransitionBlock'])].copy()
    science['base_field'] = science['target'].str.removesuffix('_dither')
    for name, visits in science.groupby('base_field'):
        if len(visits) != 2 or set(visits['target']) != {name, f'{name}_dither'}:
            raise ValueError(f'{name} does not have exactly one initial and one dithered visit')
        initial = visits.loc[visits['target'] == name].iloc[0]
        dither = visits.loc[visits['target'] == f'{name}_dither'].iloc[0]
        separation = (dither['start_time_mjd'] - initial['start_time_mjd']) * u.day
        if separation < min_separation - 1 * u.millisecond:
            raise ValueError(f'{name} revisited after only {separation.to_value(u.minute):.2f} minutes')
    return len(science) // 2


def validate_slews(obs_plan):
    """Require each move to fit in readout plus any scheduled idle time."""
    science = obs_plan[~obs_plan['target'].isin(['Unused Time', 'TransitionBlock'])]
    previous = None
    for _, row in science.iterrows():
        start = Time(row['start time (UTC)'])
        end = Time(row['end time (UTC)'])
        coord = SkyCoord(row['ra'] * u.deg, row['dec'] * u.deg)
        if previous is not None:
            previous_coord, previous_end = previous
            required = coord.separation(previous_coord) / slew_rate
            available = read_out + (start - previous_end)
            if required > available + 1 * u.millisecond:
                raise ValueError(f"Slew to {row['target']} needs {required.to_value(u.second):.1f} s; "
                                 f'only {available.to_value(u.second):.1f} s available')
        previous = coord, end


def get_obs_plan(night_start, night_end, output_path, db_path=LS4_field_grid_db_path,
                 show_plots=True):
    """Plan complete initial/dither pairs from the available grid for each block."""
    grid = load_field_grid(db_path)
    names = grid['Field Name'].to_numpy()
    programs = grid['program_id'].to_numpy()
    last_scheduled = grid['last_scheduled_mjd'].fillna(0).to_numpy(dtype=float)
    coords = SkyCoord(grid['ra_deg'].to_numpy() * u.deg,
                      grid['dec_deg'].to_numpy() * u.deg)
    unit_vectors = coords.cartesian.xyz.value.T
    dither_coords = SkyCoord((grid['ra_deg'].to_numpy() - half_field_offset_ra.to_value(u.deg)) * u.deg,
                             grid['dec_deg'].to_numpy() * u.deg)
    used = np.zeros(len(grid), dtype=bool)
    rows = []
    current_time = night_start
    exposure_duration = exp + read_out
    min_revisit = default_block_duration
    total_pairs = 0
    galactic_pairs = 0
    previous_coord = None
    previous_end = None
    block_number = 0

    while (night_end - current_time) >= 2 * min_revisit:
        remaining = night_end - current_time
        block_duration = min_revisit if remaining >= 4 * min_revisit else remaining / 2
        slots = int((block_duration / exposure_duration).decompose().value + 1e-9)
        offsets = np.arange(slots) * exposure_duration
        first_starts = current_time + offsets
        second_starts = first_starts + block_duration
        first_ok = _visibility_at_slots(coords, first_starts, exposure_duration)
        second_ok = _visibility_at_slots(dither_coords, second_starts, exposure_duration)
        sidereal_degrees = first_starts.sidereal_time('mean', longitude=LS4.location.lon).deg
        first_observations = []
        second_observations = []
        last_choice = None
        last_second_end = None
        anchor = None
        pair_ok = first_ok & second_ok

        for slot in range(slots):
            available = np.flatnonzero(pair_ok[:, slot] & ~used)
            if not len(available):
                continue

            if anchor is not None:
                in_cluster = coords[available].separation(coords[anchor]) <= max_cluster_radius
                available = available[in_cluster]
                if not len(available):
                    continue

            if previous_coord is None:
                slew_degrees = np.zeros(len(available))
            else:
                slew_degrees = coords[available].separation(previous_coord).deg
                idle = max(0, (first_starts[slot] - previous_end).to_value(u.second)) * u.second
                reachable = (slew_degrees * u.deg / slew_rate) <= read_out + idle
                available = available[reachable]
                slew_degrees = slew_degrees[reachable]
                if not len(available):
                    continue

            if anchor is not None:
                # Program allocation is a preference within the current sky
                # cluster; it must not force a jump to a distant field.
                preferred = 0 if galactic_pairs < 0.25 * (total_pairs + 1) else 1
                in_program = programs[available] == preferred
                if np.any(in_program):
                    available = available[in_program]
                    slew_degrees = slew_degrees[in_program]

            age_days = np.where(last_scheduled[available] <= 0, 30,
                                np.clip(first_starts[slot].mjd - last_scheduled[available], 0, 30))
            hour_angle = (sidereal_degrees[slot] - grid['ra_deg'].to_numpy()[available] + 180) % 360 - 180
            if anchor is None:
                # Start in a dense region so the rest of the block can remain
                # nearby. Prefer an overdue region when densities are similar.
                future = unit_vectors[np.flatnonzero(np.any(pair_ok, axis=1) & ~used)]
                neighbors = (unit_vectors[available] @ future.T >=
                             np.cos(max_cluster_radius.to_value(u.rad))).sum(axis=1)
                score = 2 * neighbors + 0.2 * age_days + 0.5 * hour_angle / 90 - 0.4 * slew_degrees
            else:
                # Inside a cluster, favor adjacent fields over distant ones.
                score = -1.5 * slew_degrees + 0.05 * age_days + 0.2 * hour_angle / 90
            choice = available[np.argmax(score)]
            if anchor is None:
                anchor = choice
            used[choice] = True
            last_choice = choice
            total_pairs += 1
            galactic_pairs += int(programs[choice] == 0)

            first_start = first_starts[slot]
            second_start = second_starts[slot]
            first_observations.append((names[choice], first_start, first_start + exposure_duration,
                                       coords[choice].ra.deg, coords[choice].dec.deg, False))
            second_observations.append((f'{names[choice]}_dither', second_start,
                                        second_start + exposure_duration,
                                        dither_coords[choice].ra.deg, dither_coords[choice].dec.deg, True))
            previous_coord = coords[choice]
            previous_end = first_start + exposure_duration
            last_second_end = second_start + exposure_duration

        first_end = current_time + block_duration
        second_end = first_end + block_duration
        rows.extend(_block_rows(first_observations, current_time, first_end, block_number))
        rows.extend(_block_rows(second_observations, first_end, second_end, block_number + 1))
        if last_choice is not None:
            previous_coord = dither_coords[last_choice]
            previous_end = last_second_end
        print(f'Blocks {block_number}-{block_number + 1}: {len(first_observations)} complete pairs')
        current_time = second_end
        block_number += 2

    if (night_end - current_time).to_value(u.second) > 1e-3:
        rows.append(_row('Unused Time', current_time, night_end, block_number=block_number))

    obs_plan = pd.DataFrame(rows)
    completed = validate_pairs(obs_plan)
    validate_slews(obs_plan)
    assert completed == total_pairs
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    obs_plan.to_csv(output_path, index=False)
    recorded = update_last_scheduled_mjd(obs_plan, db_path)
    assert recorded == completed, 'Only complete pairs should be recorded'
    print(f'Wrote {output_path}: {2 * completed} images in {completed} complete pairs')
    compute_time_efficiency(obs_plan, night_start, night_end)
    if show_plots:
        visits = compute_union_area(obs_plan, nside=32, night_start=night_start,
                                    night_end=night_end)
        plot_coverage_map(visits)
    return obs_plan


if __name__ == "__main__":

    start_time = time.time()

    args = argument_parser()
    if args.mjd:
        mjd = Time(args.mjd, format='mjd')
    else:
        mjd = current_time  

    mjd = mjd - (0.5 * u.day) # turn this on at night


    night_start = LS4.twilight_evening_astronomical(Time(mjd, format='mjd'), which='next')
    night_end = LS4.twilight_morning_astronomical(Time(mjd, format='mjd'), which='next')
    
    if args.output is not None:
        output_path = args.output
    else:
        output_path = f"plans/{Time(mjd, format='mjd').to_datetime().strftime('%Y%m%d')}.csv"

    get_obs_plan(night_start, night_end, output_path, db_path=args.db,
                 show_plots=not args.no_plots)

    end_time = time.time()
    print(f"Scheduling took {end_time - start_time:.2f} seconds.")
