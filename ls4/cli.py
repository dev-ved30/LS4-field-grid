"""Command-line entry point: `python -m ls4 <command>` (or `ls4 <command>` once installed).

Commands:
    build-grid  Tessellate the sky and write the grid CSV and scheduling database
    plan        Schedule one night and record it in the grid database
    simulate    Schedule consecutive nights (or replay saved plans) on a fresh database
    analyze     Survey statistics, cadence, and sky maps for a set of plans
    bands       Sky coverage and cadence in g, i, and z for a set of plans
    viz         Interactive plot, coverage map, or GIFs for one plan
    export      Convert a plan CSV to the telescope's .obsplan format
"""

import argparse
import time
from datetime import date, datetime, timedelta, timezone
from pathlib import Path

import pandas as pd
from astropy.time import Time

from ls4.config import GRID_CSV, GRID_DB, OUTPUTS_DIR, PLANS_DIR, read_out, slew_rate


def tonight():
    """Local evening date at La Silla: before ~noon UTC we are still in last night."""
    return (datetime.now(timezone.utc) - timedelta(hours=12)).date()


def plan_paths(paths, limit=None, default_dir=PLANS_DIR):
    paths = sorted(Path(p) for p in paths) if paths else sorted(Path(default_dir).glob("*.csv"))
    if limit is not None:
        if limit <= 0 or len(paths) < limit:
            raise SystemExit(f"--limit requires 1 to {len(paths)} plans")
        paths = paths[:limit]
    if not paths:
        raise SystemExit("No plan CSVs found.")
    return paths


# ---- commands ----------------------------------------------------------------

def cmd_build_grid(args):
    from ls4 import grid as grid_mod

    if Path(args.db).exists() and not args.force:
        raise SystemExit(f"{args.db} exists and holds the scheduling history; "
                         "pass --force to overwrite it.")
    grid = grid_mod.build_grid()
    Path(args.csv).parent.mkdir(parents=True, exist_ok=True)
    grid.to_csv(args.csv, index=False)
    grid_mod.write_db(grid, args.db, overwrite=True)
    print(f"Wrote {len(grid)} fields to {args.csv} and {args.db}")
    if args.show:
        from ls4.plots import field_grid
        field_grid(grid)


def cmd_plan(args):
    from ls4.metrics import time_efficiency, union_area
    from ls4.scheduler import night_window, plan_night

    if args.mjd is not None:
        day = Time(args.mjd, format="mjd").to_datetime().date()
    else:
        day = args.date or tonight()
    start, end = night_window(day)
    output = Path(args.output) if args.output else PLANS_DIR / f"{day:%Y%m%d}.csv"
    print(f"Planning night of {day}: {start.iso} to {end.iso} UTC")
    t0 = time.time()
    plan = plan_night(start, end, output, db_path=args.db)
    time_efficiency(plan, start, end)
    if args.coverage:
        from ls4.plots import coverage_map
        coverage_map(union_area(plan, nside=32, night_start=start, night_end=end))
    print(f"Scheduling took {time.time() - t0:.2f} seconds.")


def cmd_simulate(args):
    from ls4.grid import fresh_db
    from ls4.metrics import analyze_survey, night_metrics
    from ls4.scheduler import night_window, plan_night

    if args.replay:
        baselines = plan_paths([], args.limit, default_dir=args.replay)
        nights = []
        for path in baselines:
            original = pd.read_csv(path)
            nights.append((path.stem,
                           Time(original.iloc[0]["start time (UTC)"], scale="utc"),
                           Time(original.iloc[-1]["end time (UTC)"], scale="utc"),
                           original))
    else:
        if args.nights <= 0:
            raise SystemExit("--nights must be positive")
        nights = []
        for offset in range(args.nights):
            day = args.start + timedelta(days=offset)
            nights.append((f"{day:%Y%m%d}", *night_window(day), None))

    out = Path(args.output_dir) if args.output_dir else OUTPUTS_DIR / args.name
    out.mkdir(parents=True, exist_ok=True)
    db_path = fresh_db(out / "grid.db")

    rows, paths = [], []
    for i, (name, start, end, original) in enumerate(nights, start=1):
        path = out / f"{name}.csv"
        print(f"\nNight {i}/{len(nights)}: {name} {start.iso} to {end.iso}", flush=True)
        result = night_metrics(plan_night(start, end, path, db_path=db_path))
        if result["unpaired_images"] or result["short_pairs"]:
            raise AssertionError(f"Pair or revisit violation in {path}")
        row = {"night": name, "start_utc": start.iso, "end_utc": end.iso,
               "night_hours": (end - start).to_value("hour"), **result}
        if original is not None:
            row.update({f"baseline_{key}": value for key, value in night_metrics(original).items()})
        rows.append(row)
        paths.append(path)
        pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)

    summary = pd.DataFrame(rows)
    fields = analyze_survey(paths, out / "stats", db_path=db_path)
    write_report(out, summary, fields, replay=args.replay)

    if args.gifs:
        from ls4.animate import make_gif
        for path in paths[:args.gifs]:
            make_gif(path, out / "gifs" / f"{path.stem}.gif", mode="pointings")
            make_gif(path, out / "gifs" / f"{path.stem}_visits.gif", mode="visits")
    print(f"Completed {len(summary)} nights; {summary['images'].sum():,} images in "
          f"{summary['complete_pairs'].sum():,} complete pairs -> {out}")


def write_report(out, summary, fields, replay=None):
    source = (f"Replay of the saved plans in `{replay}`, using each plan's start and end time."
              if replay else
              f"Consecutive nights {summary['night'].iloc[0]} through {summary['night'].iloc[-1]}, "
              "each from astronomical evening twilight to astronomical morning twilight at La Silla.")
    lines = [
        f"# {len(summary)}-night LS4 simulation", "", source,
        "The scheduler started with a fresh grid history and used a "
        f"{slew_rate.to_value('deg/s'):.2f}°/second slew rate with {read_out.to_value('s'):.0f}-second readout.",
        "", "| Metric | Result |", "| --- | ---: |",
        f"| Science images | {int(summary['images'].sum()):,} |",
        f"| Complete initial/dither pairs | {int(summary['complete_pairs'].sum()):,} |",
        f"| SOLE (Galactic) pairs | {int(summary['sole_pairs'].sum()):,} "
        f"({summary['sole_pairs'].sum() / summary['complete_pairs'].sum():.1%}) |",
        f"| Unpaired images | {int(summary['unpaired_images'].sum())} |",
        f"| Revisit pairs under 30 minutes | {int(summary['short_pairs'].sum())} |",
        f"| Unused hours | {summary['unused_minutes'].sum() / 60:.2f} |",
        f"| Unique fields visited | {int((fields['visits'] > 0).sum()):,} of {len(fields):,} |",
        f"| Fields seen in g, i, and z | {int((fields['dither_phases'] == 2).sum()):,} |",
        f"| DB histories matching plans | {int(fields['history_matches'].sum()):,} of {len(fields):,} |",
    ]
    if replay:
        lines += ["", "| Total | Saved plans | Simulated |", "| --- | ---: | ---: |"]
        for key in ("images", "complete_pairs", "unpaired_images", "short_pairs", "unused_minutes"):
            lines.append(f"| {key} | {summary[f'baseline_{key}'].sum():.1f} | {summary[key].sum():.1f} |")
    lines += ["", "`summary.csv` has the nightly metrics; `stats/` has field-level statistics and sky maps."]
    (out / "REPORT.md").write_text("\n".join(lines) + "\n")


def cmd_analyze(args):
    from ls4.metrics import analyze_survey

    analyze_survey(plan_paths(args.plans, args.limit), args.output_dir, grid_csv=args.grid,
                   db_path=args.db, sample_minutes=args.visibility_sample_minutes, show=args.show,
                   band_nside=args.band_nside)


def cmd_bands(args):
    from ls4.bands import write_band_report
    from ls4.metrics import load_plans

    write_band_report(load_plans(plan_paths(args.plans, args.limit)), args.output_dir,
                      nside=args.nside, show=args.show)


def cmd_viz(args):
    plan_csv = Path(args.plan)
    out_dir = Path(args.output_dir) if args.output_dir else OUTPUTS_DIR / "viz"
    any_output = args.gif or args.visits_gif or args.coverage or args.html
    if args.html or not any_output:
        from ls4.plots import observing_plan
        observing_plan(pd.read_csv(plan_csv),
                       output=out_dir / f"{plan_csv.stem}.html" if args.html else None,
                       show=not args.html)
    if args.coverage:
        from ls4.metrics import union_area
        from ls4.plots import coverage_map
        coverage_map(union_area(pd.read_csv(plan_csv), nside=32), title=f"Coverage: {plan_csv.stem}")
    gif_options = dict(grid_path=args.field_grid, cadence_minutes=args.cadence_minutes,
                       duration=args.duration, width=args.width, height=args.height,
                       show_galaxy=not args.no_galaxy)
    if args.gif or args.visits_gif:
        from ls4.animate import make_gif
        if args.gif:
            make_gif(plan_csv, out_dir / f"{plan_csv.stem}.gif", mode="pointings", **gif_options)
        if args.visits_gif:
            make_gif(plan_csv, out_dir / f"{plan_csv.stem}_visits.gif", mode="visits", **gif_options)


def cmd_export(args):
    from ls4.export import write_obsplan

    stamp = args.date or (datetime.now(timezone.utc) + timedelta(days=1)).strftime("%Y%m%d")
    write_obsplan(args.plan, args.output or f"{stamp}.obsplan")


# ---- parser ------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(prog="ls4", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("build-grid", help="Build the field grid CSV and scheduling database")
    p.add_argument("--csv", type=Path, default=GRID_CSV)
    p.add_argument("--db", type=Path, default=GRID_DB)
    p.add_argument("--force", action="store_true", help="Overwrite an existing database (erases history)")
    p.add_argument("--show", action="store_true", help="Plot the grid with example footprints")
    p.set_defaults(func=cmd_build_grid)

    p = sub.add_parser("plan", help="Schedule one night")
    when = p.add_mutually_exclusive_group()
    when.add_argument("--date", type=date.fromisoformat,
                      help="Local evening date, YYYY-MM-DD (default: tonight at La Silla)")
    when.add_argument("--mjd", type=float, help="Any MJD on the evening's UTC date")
    p.add_argument("--output", type=Path, help="Plan CSV (default: plans/YYYYMMDD.csv)")
    p.add_argument("--db", type=Path, default=GRID_DB, help="Grid database whose history is updated")
    p.add_argument("--coverage", action="store_true", help="Compute sky area and show a coverage map")
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("simulate", help="Schedule many nights on a fresh database")
    source = p.add_mutually_exclusive_group()
    source.add_argument("--start", type=date.fromisoformat, default=date(2026, 7, 17),
                        help="First evening, YYYY-MM-DD (default: 2026-07-17)")
    source.add_argument("--replay", type=Path, metavar="DIR",
                        help="Re-plan each saved plan's time window and compare against it")
    p.add_argument("--nights", type=int, default=30, help="Nights to simulate from --start (default: 30)")
    p.add_argument("--limit", type=int, help="With --replay, only the first N plans")
    p.add_argument("--name", default="simulation", help="Run name under outputs/ (default: simulation)")
    p.add_argument("--output-dir", type=Path, help="Explicit output directory (overrides --name)")
    p.add_argument("--gifs", type=int, default=0, metavar="N", help="Make GIFs for the first N nights")
    p.set_defaults(func=cmd_simulate)

    p = sub.add_parser("analyze", help="Survey statistics for a set of plans")
    p.add_argument("plans", nargs="*", type=Path, help="Plan CSVs (default: plans/*.csv)")
    p.add_argument("--limit", type=int, help="Only the first N plans after sorting by name")
    p.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "survey_stats")
    p.add_argument("--grid", type=Path, default=GRID_CSV, help="Grid CSV (includes unvisited fields)")
    p.add_argument("--db", type=Path, help="Verify every field's history in this database against the plans")
    p.add_argument("--visibility-sample-minutes", type=int, default=10)
    p.add_argument("--band-nside", type=int, default=128,
                   help="HEALPix nside for the per-band maps in bands/ (0 skips them)")
    p.add_argument("--show", action="store_true", help="Also open the Plotly figures")
    p.set_defaults(func=cmd_analyze)

    p = sub.add_parser("bands", help="Sky coverage and cadence in g, i, and z")
    p.add_argument("plans", nargs="*", type=Path, help="Plan CSVs (default: plans/*.csv)")
    p.add_argument("--limit", type=int, help="Only the first N plans after sorting by name")
    p.add_argument("--output-dir", type=Path, default=OUTPUTS_DIR / "bands")
    p.add_argument("--nside", type=int, default=128, help="HEALPix resolution (128 is ~0.46° pixels)")
    p.add_argument("--show", action="store_true", help="Also open the Plotly figures")
    p.set_defaults(func=cmd_bands)

    p = sub.add_parser("viz", help="Visualize one plan (interactive plot by default)")
    p.add_argument("plan", type=Path)
    p.add_argument("--html", action="store_true", help="Save the interactive plot as HTML instead of opening it")
    p.add_argument("--coverage", action="store_true", help="Show a HEALPix coverage map")
    p.add_argument("--gif", action="store_true", help="GIF colored by first pointing / dither")
    p.add_argument("--visits-gif", action="store_true", help="GIF colored by cumulative visits")
    p.add_argument("--output-dir", type=Path, help="Where to write files (default: outputs/viz)")
    p.add_argument("--field-grid", type=Path, default=GRID_CSV)
    p.add_argument("--cadence-minutes", type=int, default=15, help="GIF time between frames")
    p.add_argument("--duration", type=int, default=250, help="GIF milliseconds per frame")
    p.add_argument("--width", type=int, default=1920)
    p.add_argument("--height", type=int, default=1200)
    p.add_argument("--no-galaxy", action="store_true", help="Omit the galactic-plane overlay")
    p.set_defaults(func=cmd_viz)

    p = sub.add_parser("export", help="Write a telescope .obsplan from a plan CSV")
    p.add_argument("plan", type=Path)
    p.add_argument("-d", "--date", help="YYYYMMDD for the file name (default: tomorrow, UTC)")
    p.add_argument("--output", type=Path, help="Output path (default: <date>.obsplan)")
    p.set_defaults(func=cmd_export)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
