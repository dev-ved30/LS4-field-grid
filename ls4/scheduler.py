"""Nightly scheduler: complete initial/dither pairs in clustered 30-minute blocks."""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import NonRotationTransformationWarning, SkyCoord
from astropy.time import Time

from ls4.config import (GRID_DB, LS4, NON_SCIENCE_TARGETS, default_block_duration, exp,
                        global_constraints, half_field_offset_ra, max_cluster_radius, read_out,
                        slew_rate)
from ls4.grid import load_schedulable, record_pairs

warnings.filterwarnings("ignore", category=NonRotationTransformationWarning)


def night_window(day):
    """Return the astronomical-twilight window starting on the evening of `day` (local)."""
    noon_utc = Time(f"{day.isoformat()}T12:00:00", scale="utc")
    start = LS4.twilight_evening_astronomical(noon_utc, which="next")
    end = LS4.twilight_morning_astronomical(start, which="next")
    if end <= start:
        raise ValueError(f"Invalid twilight window for {day}")
    return start, end


def science_pointings(plan):
    """Science rows of a plan, with `base_field` and `is_dither` columns added."""
    science = plan[~plan["target"].isin(NON_SCIENCE_TARGETS)].copy()
    science["base_field"] = science["target"].str.removesuffix("_dither")
    science["is_dither"] = science["target"].str.endswith("_dither")
    return science


def visibility_at_slots(coords, starts, exposure_duration):
    """Check every constraint at the start and end of each proposed observation."""
    endpoints = np.column_stack((starts.mjd, (starts + exposure_duration).mjd)).ravel()
    times = Time(endpoints, format="mjd", scale="utc")
    allowed = np.logical_and.reduce([
        constraint(LS4, coords, times=times, grid_times_targets=True)
        for constraint in global_constraints
    ])
    return allowed[:, 0::2] & allowed[:, 1::2]


def dither_coords(ra_deg, dec_deg):
    return SkyCoord((np.asarray(ra_deg) - half_field_offset_ra.to_value(u.deg)) * u.deg,
                    np.asarray(dec_deg) * u.deg)


def _row(target, start, end, ra=None, dec=None, dither=None, block_number=0):
    return {
        "target": target,
        "start time (UTC)": start.to_datetime().isoformat(sep=" "),
        "end time (UTC)": end.to_datetime().isoformat(sep=" "),
        "duration (minutes)": (end - start).to_value(u.minute),
        "ra": ra,
        "dec": dec,
        "configuration": str({"dither": dither}) if dither is not None else "",
        "start_time_mjd": start.mjd,
        "block_number": block_number,
        "ra_hr": ra / 15 if ra is not None else None,
    }


def _block_rows(observations, block_start, block_end, block_number):
    """Emit science pointings and the actual gaps between them."""
    rows = []
    cursor = block_start
    for target, start, end, ra, dec, dither in sorted(observations, key=lambda item: item[1].mjd):
        if (start - cursor).to_value(u.second) > 1e-3:
            rows.append(_row("Unused Time", cursor, start, block_number=block_number))
        rows.append(_row(target, start, end, ra, dec, dither, block_number))
        cursor = end
    if (block_end - cursor).to_value(u.second) > 1e-3:
        rows.append(_row("Unused Time", cursor, block_end, block_number=block_number))
    return rows


def complete_pairs(science):
    """Per-field visit summary for fields with exactly one initial and one dither."""
    pairs = science.groupby("base_field").agg(
        visits=("target", "size"), dithers=("is_dither", "sum"),
        last_mjd=("start_time_mjd", "max"))
    return pairs[(pairs["visits"] == 2) & (pairs["dithers"] == 1)]


def validate_pairs(obs_plan, min_separation=default_block_duration):
    """Reject missing, duplicate, reversed, or too-close visits before saving."""
    starts = pd.to_datetime(obs_plan["start time (UTC)"])
    ends = pd.to_datetime(obs_plan["end time (UTC)"])
    if (ends <= starts).any() or (starts.iloc[1:].to_numpy() < ends.iloc[:-1].to_numpy()).any():
        raise ValueError("Observation rows overlap or have nonpositive duration")
    science = science_pointings(obs_plan)
    for name, visits in science.groupby("base_field"):
        if len(visits) != 2 or set(visits["target"]) != {name, f"{name}_dither"}:
            raise ValueError(f"{name} does not have exactly one initial and one dithered visit")
        initial = visits.loc[visits["target"] == name].iloc[0]
        dither = visits.loc[visits["target"] == f"{name}_dither"].iloc[0]
        separation = (dither["start_time_mjd"] - initial["start_time_mjd"]) * u.day
        if separation < min_separation - 1 * u.millisecond:
            raise ValueError(f"{name} revisited after only {separation.to_value(u.minute):.2f} minutes")
    return len(science) // 2


def validate_slews(obs_plan):
    """Require each move to fit in readout plus any scheduled idle time."""
    science = science_pointings(obs_plan)
    previous = None
    for _, row in science.iterrows():
        start = Time(row["start time (UTC)"])
        end = Time(row["end time (UTC)"])
        coord = SkyCoord(row["ra"] * u.deg, row["dec"] * u.deg)
        if previous is not None:
            previous_coord, previous_end = previous
            required = coord.separation(previous_coord) / slew_rate
            available = read_out + (start - previous_end)
            if required > available + 1 * u.millisecond:
                raise ValueError(f"Slew to {row['target']} needs {required.to_value(u.second):.1f} s; "
                                 f"only {available.to_value(u.second):.1f} s available")
        previous = coord, end


def plan_night(night_start, night_end, output_path, db_path=GRID_DB):
    """Plan complete initial/dither pairs for one night, save the CSV, and update history."""
    grid = load_schedulable(db_path)
    names = grid["Field Name"].to_numpy()
    programs = grid["program_id"].to_numpy()
    last_scheduled = grid["last_scheduled_mjd"].fillna(0).to_numpy(dtype=float)
    coords = SkyCoord(grid["ra_deg"].to_numpy() * u.deg,
                      grid["dec_deg"].to_numpy() * u.deg)
    unit_vectors = coords.cartesian.xyz.value.T
    dithers = dither_coords(grid["ra_deg"].to_numpy(), grid["dec_deg"].to_numpy())
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
        first_ok = visibility_at_slots(coords, first_starts, exposure_duration)
        second_ok = visibility_at_slots(dithers, second_starts, exposure_duration)
        sidereal_degrees = first_starts.sidereal_time("mean", longitude=LS4.location.lon).deg
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
            hour_angle = (sidereal_degrees[slot] - grid["ra_deg"].to_numpy()[available] + 180) % 360 - 180
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
            second_observations.append((f"{names[choice]}_dither", second_start,
                                        second_start + exposure_duration,
                                        dithers[choice].ra.deg, dithers[choice].dec.deg, True))
            previous_coord = coords[choice]
            previous_end = first_start + exposure_duration
            last_second_end = second_start + exposure_duration

        first_end = current_time + block_duration
        second_end = first_end + block_duration
        rows.extend(_block_rows(first_observations, current_time, first_end, block_number))
        rows.extend(_block_rows(second_observations, first_end, second_end, block_number + 1))
        if last_choice is not None:
            previous_coord = dithers[last_choice]
            previous_end = last_second_end
        print(f"Blocks {block_number}-{block_number + 1}: {len(first_observations)} complete pairs")
        current_time = second_end
        block_number += 2

    if (night_end - current_time).to_value(u.second) > 1e-3:
        rows.append(_row("Unused Time", current_time, night_end, block_number=block_number))

    obs_plan = pd.DataFrame(rows)
    completed = validate_pairs(obs_plan)
    validate_slews(obs_plan)
    assert completed == total_pairs
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    obs_plan.to_csv(output_path, index=False)

    # Advance history only for fields with both scheduled pointings.
    pairs = complete_pairs(science_pointings(obs_plan))
    record_pairs(pairs["last_mjd"], db_path)
    assert len(pairs) == completed, "Only complete pairs should be recorded"
    print(f"Wrote {output_path}: {2 * completed} images in {completed} complete pairs")
    return obs_plan
