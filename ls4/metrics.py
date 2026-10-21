"""Per-night schedule metrics and multi-night survey statistics."""

from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from ls4.config import (FOV_length, FOV_width, GRID_CSV, SOLE, default_block_duration, exp,
                        min_declination, read_out, slew_rate)
from ls4.grid import load_grid_csv, load_history, program_ids
from ls4.scheduler import complete_pairs, dither_coords, science_pointings, visibility_at_slots


# ---- single night ---------------------------------------------------------

def night_metrics(plan):
    """Image counts, pairing, unused time, and slew statistics for one plan."""
    science = science_pointings(plan)
    complete = complete_pairs(science)
    short = 0
    for name in complete.index:
        visits = science[science["base_field"] == name].sort_values("start_time_mjd")
        if (visits.iloc[1]["start_time_mjd"] - visits.iloc[0]["start_time_mjd"]) * 1440 < 30 - 1e-5:
            short += 1
    ra = np.deg2rad(science["ra"].to_numpy())
    dec = np.deg2rad(science["dec"].to_numpy())
    if len(science) > 1:
        cosine = (np.sin(dec[:-1]) * np.sin(dec[1:]) +
                  np.cos(dec[:-1]) * np.cos(dec[1:]) * np.cos(ra[:-1] - ra[1:]))
        steps = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    else:
        steps = []
    grid = load_grid_csv()
    sole_fields = set(grid.loc[program_ids(grid["l"], grid["b"]) == SOLE, "base_field"])
    slew_speed = slew_rate.to_value("deg/s")
    readout_seconds = read_out.to_value("s")
    return {
        "images": len(science),
        "complete_pairs": len(complete),
        "sole_pairs": int(complete.index.isin(sole_fields).sum()),
        "unpaired_images": len(science) - 2 * len(complete),
        "short_pairs": short,
        "unused_minutes": plan.loc[plan["target"] == "Unused Time", "duration (minutes)"].sum(),
        "median_slew_deg": float(np.median(steps)) if len(steps) else 0,
        "p90_slew_deg": float(np.percentile(steps, 90)) if len(steps) else 0,
        "max_slew_deg": float(np.max(steps)) if len(steps) else 0,
        "slews_over_readout": int(np.sum(np.asarray(steps) / slew_speed > readout_seconds)),
    }


def time_efficiency(obs_plan, night_start, night_end):
    total = (night_end - night_start).to_value(u.minute) * u.minute
    unused = obs_plan.loc[obs_plan["target"] == "Unused Time", "duration (minutes)"].sum() * u.minute
    efficiency = (total - unused) / total
    print(f"Total night duration: {total.to_value(u.hour):.2f} hours")
    print(f"Total unused time: {unused.to_value(u.hour):.2f} hours")
    print(f"Time efficiency: {efficiency:.2%}")
    return efficiency


def theoretical_max_area(night_start, night_end):
    """Area (deg^2) if every slot took a new field; halved for the dithered revisit."""
    night_duration = (night_end - night_start).to_value(u.second) * u.second
    max_fields = (night_duration / (exp + read_out)).to_value()
    return (max_fields * FOV_width * FOV_length).to_value(u.deg**2) / 2


def union_area(obs_plan, nside=2048, night_start=None, night_end=None):
    """Print the union area of all footprints and return a HEALPix visit-count map."""
    science = science_pointings(obs_plan)
    visits = np.zeros(hp.nside2npix(nside), dtype=int)
    half_w = FOV_width.to_value(u.deg) / 2
    half_h = FOV_length.to_value(u.deg) / 2
    offsets = [(-half_w, -half_h), (half_w, -half_h), (half_w, half_h), (-half_w, half_h)]

    for ra, dec in zip(science["ra"].to_numpy(), science["dec"].to_numpy()):
        center = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
        corners = [center.spherical_offsets_by(dx * u.deg, dy * u.deg) for dx, dy in offsets]
        vecs = np.array([hp.ang2vec(np.radians(90 - c.dec.deg), np.radians(c.ra.deg)) for c in corners])
        visits[hp.query_polygon(nside, vecs)] += 1

    pixel_area_deg2 = hp.nside2pixarea(nside, degrees=True)
    total_area = (visits > 0).sum() * pixel_area_deg2
    print(f"Total observed area: {total_area:.2f} deg^2")
    if night_start is not None and night_end is not None:
        max_area = theoretical_max_area(night_start, night_end)
        print(f"Maximum area possible for the night based on exposure time and readout time: {max_area:.2f} deg^2")
        print(f"Area efficiency: {100 * total_area / max_area:.2f}%")
    for i in range(1, min(10, visits.max() + 1)):
        print(f"Area covered by {i} visits: {(visits == i).sum() * pixel_area_deg2:.2f} deg^2")
    return visits


# ---- multi-night survey ---------------------------------------------------

def load_plans(paths):
    """Science pointings from every plan, with a `night` column."""
    frames = []
    for path in paths:
        df = pd.read_csv(path)
        df["night"] = Path(path).stem
        frames.append(df)
    plans = science_pointings(pd.concat(frames, ignore_index=True)).dropna(subset=["ra", "dec"])
    plans["utc"] = pd.to_datetime(plans["start time (UTC)"])
    return plans.sort_values("utc").reset_index(drop=True)


def audit_history(per_field, db_path):
    """Check that each database timestamp equals its last planned pointing."""
    history = load_history(db_path)
    if history["base_field"].duplicated().any():
        raise ValueError(f"Duplicate field names in {db_path}")
    checked = per_field.merge(history, on="base_field", how="outer", indicator=True,
                              validate="one_to_one")
    expected = checked["last_plan_mjd"].fillna(0).to_numpy(dtype=float)
    actual = checked["last_scheduled_mjd"].fillna(0).to_numpy(dtype=float)
    checked["history_matches"] = (checked["_merge"] == "both") & np.isclose(
        expected, actual, rtol=0, atol=1e-6)
    mismatches = checked.loc[~checked["history_matches"], "base_field"].tolist()
    if mismatches:
        raise ValueError(f"{len(mismatches)} DB histories disagree with plans; "
                         f"first fields: {mismatches[:10]}")
    return checked.drop(columns="_merge")


def sampled_visibility(grid, plan_paths, sample_minutes):
    """Count nights when a field has at least one sampled initial/dither opportunity."""
    if sample_minutes <= 0:
        raise ValueError("Visibility sample interval must be positive")
    eligible = grid["dec"].to_numpy() > min_declination.to_value(u.deg)
    visible_nights = np.zeros(len(grid), dtype=int)
    positions = np.flatnonzero(eligible)
    ra = grid.iloc[positions]["ra"].to_numpy()
    dec = grid.iloc[positions]["dec"].to_numpy()
    initial = SkyCoord(ra * u.deg, dec * u.deg)
    dither = dither_coords(ra, dec)
    duration = exp + read_out
    step = sample_minutes * u.minute
    for path in plan_paths:
        plan = pd.read_csv(path)
        start = Time(plan.iloc[0]["start time (UTC)"], scale="utc")
        end = Time(plan.iloc[-1]["end time (UTC)"], scale="utc")
        available = end - start - default_block_duration - duration
        count = int(np.floor((available / step).decompose().value)) + 1
        if count <= 0:
            continue
        starts = start + np.arange(count) * step
        possible = (visibility_at_slots(initial, starts, duration) &
                    visibility_at_slots(dither, starts + default_block_duration, duration))
        visible_nights[positions] += possible.any(axis=1)
    return eligible, visible_nights


def nightly_summary(plan_paths):
    rows = []
    for path in plan_paths:
        df = pd.read_csv(path)
        used = science_pointings(df)
        row = {
            "night": Path(path).stem,
            "images": len(used.dropna(subset=["ra"])),
            "used_minutes": used["duration (minutes)"].clip(lower=0).sum(),
            "unused_minutes": df.loc[df["target"] == "Unused Time", "duration (minutes)"].sum(),
        }
        row["night_minutes"] = row["used_minutes"] + row["unused_minutes"]
        rows.append(row)
    return pd.DataFrame(rows)


def analyze_survey(plan_paths, output_dir, grid_csv=GRID_CSV, db_path=None,
                   sample_minutes=10, show=False):
    """Write per-field, cadence, and nightly statistics plus figures for a set of plans."""
    from ls4 import plots
    from ls4.sky_maps import write_sky_grid_maps

    plan_paths = [Path(path) for path in plan_paths]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f'Analyzing {len(plan_paths)} plans: {", ".join(p.name for p in plan_paths)}\n')
    plans = load_plans(plan_paths)
    grid = load_grid_csv(grid_csv)

    # ---- per-field visit counts ---------------------------------------
    night_fields = plans.groupby(["night", "base_field"]).agg(
        images=("target", "size"), dithers=("is_dither", "sum")).reset_index()
    complete = night_fields[(night_fields["images"] == 2) & (night_fields["dithers"] == 1)]
    pair_counts = complete.groupby("base_field").size().rename("complete_pairs")
    per_field = plans.groupby("base_field").agg(
        visits=("target", "count"),
        dithers=("is_dither", "sum"),
        observed_nights=("night", "nunique"),
        first_obs=("utc", "min"),
        last_obs=("utc", "max"),
        last_plan_mjd=("start_time_mjd", "max"),
    ).reset_index()
    unknown = per_field.loc[~per_field["base_field"].isin(grid["base_field"]), "base_field"]
    if len(unknown):
        raise ValueError(f"Plans contain fields absent from grid: {unknown.tolist()[:10]}")
    per_field = grid[["base_field", "program_id", "ra", "dec"]].merge(
        per_field, on="base_field", how="left", validate="one_to_one")
    for column in ("visits", "dithers", "observed_nights"):
        per_field[column] = per_field[column].fillna(0).astype(int)
    per_field["complete_pairs"] = per_field["base_field"].map(pair_counts).fillna(0).astype(int)
    per_field["above_declination_limit"], per_field["sampled_visible_nights"] = (
        sampled_visibility(per_field, plan_paths, sample_minutes))
    if db_path:
        per_field = audit_history(per_field, db_path)
    observed = per_field[per_field["visits"] > 0]

    # ---- cadence between separate nights (the dither is not a revisit) ---
    gap_rows = []
    gap_distributions = {}
    for name, times in plans[~plans["is_dither"]].sort_values("utc").groupby("base_field")["utc"]:
        gaps = times.diff().dropna().dt.total_seconds() / 86400.0
        if len(gaps):
            gap_distributions[name] = gaps.to_numpy()
        gap_rows.extend({"base_field": name, "gap_days": gap} for gap in gaps)
    cadence = pd.DataFrame(gap_rows)
    if len(cadence):
        per_field_cadence = cadence.groupby("base_field")["gap_days"].agg(
            median_gap_days="median", mean_gap_days="mean",
            min_gap_days="min", max_gap_days="max").reset_index()
    else:
        per_field_cadence = pd.DataFrame(columns=[
            "base_field", "median_gap_days", "mean_gap_days", "min_gap_days", "max_gap_days"])

    nightly = nightly_summary(plan_paths)

    # ---- text summary ----------------------------------------------------
    span_start, span_end = plans["utc"].min(), plans["utc"].max()
    span_days = (span_end - span_start).total_seconds() / 86400.0
    never = per_field[per_field["visits"] == 0]
    visible_never = never[never["sampled_visible_nights"] > 0]
    lines = [
        "================ survey summary ================",
        f"nights               : {len(plan_paths)}",
        f"span                 : {span_start:%Y-%m-%d} to {span_end:%Y-%m-%d} ({span_days:.1f} days)",
        f"grid fields          : {len(per_field)}",
        f'total visits         : {len(plans)} ({plans["is_dither"].sum()} dithered)',
        f"complete pairs       : {len(complete)}",
        f"unpaired images      : {len(plans) - 2 * len(complete)}",
        f"unique fields visited: {len(observed)}",
        f"fields never visited : {len(never)} "
        f'(by program: {never["program_id"].value_counts().to_dict()})',
        f'  below dec limit    : {(~never["above_declination_limit"]).sum()}',
        f"  pair visible sampled, unvisited: {len(visible_never)}",
        f"  no sampled pair opportunity: "
        f'{len(never) - len(visible_never) - (~never["above_declination_limit"]).sum()}',
        f'mean visits/field    : {observed["visits"].mean():.2f}',
    ]
    if len(cadence):
        sole_gaps = cadence.loc[cadence["base_field"].isin(
            per_field.loc[per_field["program_id"] == SOLE, "base_field"]), "gap_days"]
        if len(sole_gaps):
            lines.append(f'SOLE cadence         : median {sole_gaps.median():.2f} d, '
                         f'{(sole_gaps < 1.5).mean():.0%} of gaps under 1.5 d')
        lines.append(f'initial-visit cadence: median {cadence["gap_days"].median():.2f} d, '
                     f'mean {cadence["gap_days"].mean():.2f} d, '
                     f'max {cadence["gap_days"].max():.2f} d')
    lines.append(f'total time used      : {nightly["used_minutes"].sum()/60:.1f} h of '
                 f'{nightly["night_minutes"].sum()/60:.1f} h '
                 f'({100*nightly["used_minutes"].sum()/nightly["night_minutes"].sum():.1f}%)')
    if db_path:
        lines.append(f"replay DB histories  : {len(per_field)} of {len(per_field)} match plans")
    stale = observed.copy()
    stale["days_since_last"] = (span_end - stale["last_obs"]).dt.total_seconds() / 86400.0
    lines.append("longest time since last observation (top 10):")
    top = stale.nlargest(10, "days_since_last")[["base_field", "visits", "last_obs", "days_since_last"]]
    lines.append(top.to_string(index=False))
    lines.append(f"Visibility sampled every {sample_minutes} minutes; "
                 "an opportunity requires both pointings to pass all sky constraints "
                 "30 minutes apart, without accounting for other scheduled fields.")
    summary = "\n".join(lines) + "\n"
    print(summary)
    (output_dir / "summary.txt").write_text(summary)

    # ---- figures and tables ------------------------------------------------
    print("\nWriting figures...")
    write_sky_grid_maps(per_field, gap_distributions, output_dir)
    plots.visit_histogram(per_field, output_dir, show)
    if len(cadence):
        plots.cadence_histograms(cadence, per_field_cadence, output_dir, show)
        plots.cadence_cdf(cadence, output_dir, show)
    plots.nightly_summary(nightly, output_dir, show)

    table_path = output_dir / "per_field_stats.csv"
    per_field.merge(per_field_cadence, on="base_field", how="left").to_csv(table_path, index=False)
    print(f"  wrote {table_path}")
    nightly.to_csv(output_dir / "nightly_summary.csv", index=False)
    return per_field
