"""Production nightly packages: python -m ls4.nightly --help.

Science imports are deliberately deferred so `notify` and `check` remain usable
when installation or scheduling fails. Generation never touches the local grid DB.
"""

import argparse
from datetime import date, datetime, timedelta, timezone
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import warnings
from zoneinfo import ZoneInfo

from ls4.automation import (GitHub, Slack, footer_links, header, mention, section, sha256, slack_blocks,
                            slack_message, write_json)

ROOT = Path(__file__).resolve().parent.parent
HISTORY_NOTE = "Generated plans, assuming execution; not confirmed observing results."


def evening_date(value=None):
    """Use the CLI's noon-UTC boundary; explicitly pass backfill dates."""
    return date.fromisoformat(value) if value else (datetime.now(timezone.utc) - timedelta(hours=12)).date()


def github_client():
    return GitHub(os.environ["GITHUB_REPOSITORY"], os.environ["GITHUB_TOKEN"])


def provenance():
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()

    return {"commit": os.environ.get("GITHUB_SHA") or git("rev-parse", "HEAD"),
            "working_tree_dirty": bool(git("status", "--porcelain", "--untracked-files=normal")),
            "run_id": os.environ.get("GITHUB_RUN_ID"),
            "run_attempt": os.environ.get("GITHUB_RUN_ATTEMPT"),
            "python": sys.version.split()[0],
            "dependency_lock_sha256": sha256(ROOT / "requirements-nightly.lock")}


def restore_history(day, db, grid_csv, *, use_github):
    """Use a previous immutable snapshot, or bootstrap from saved nightly CSVs.

    Today's revisions and future plans are excluded. A network/corruption error
    aborts rather than quietly resetting cadence to an empty database.
    """
    import pandas as pd
    from ls4.grid import load_history, record_pairs
    from ls4.scheduler import complete_pairs, science_pointings

    fields = set(load_history(db)["base_field"])
    release = github_client().latest_nightly(day) if use_github else None
    if release:
        client = github_client()
        manifest = json.loads(client.read_asset(release, "manifest.json"))
        raw = client.read_asset(release, "history.json")
        import hashlib

        if hashlib.sha256(raw).hexdigest() != manifest["files"]["history.json"]:
            raise ValueError("Published scheduling history checksum mismatch")
        snapshot = json.loads(raw)
        if snapshot["grid_sha256"] != sha256(grid_csv):
            raise ValueError("Field grid changed; migrate history explicitly before production use")
        if (snapshot["source"] != "generated-plans" or manifest["generation"] != "success"
                or snapshot["through_night"] != manifest["evening_date"]
                or date.fromisoformat(snapshot["through_night"]) >= day):
            raise ValueError("Invalid scheduling history provenance")
        values = snapshot["fields"]
        if set(values) != fields:
            raise ValueError("Published history does not match the field universe")
        source = release["html_url"]
    else:
        values = dict.fromkeys(fields, 0.0)
        source = "Repository plans/*.csv (bootstrap)"
        for path in sorted((ROOT / "plans").glob("*.csv")):
            try:
                plan_day = datetime.strptime(path.stem, "%Y%m%d").date()
            except ValueError:
                raise ValueError(f"History plan needs a YYYYMMDD filename: {path.name}") from None
            if plan_day >= day:
                continue
            pairs = complete_pairs(science_pointings(pd.read_csv(path)))
            for field, mjd in pairs["last_mjd"].items():
                if field not in values:
                    raise ValueError(f"Historical field {field} is absent from the grid")
                values[field] = max(values[field], float(mjd))
    # Zero means never scheduled; reject future, NaN, and infinite history values.
    from astropy.time import Time

    cutoff = Time(f"{day.isoformat()}T12:00:00", scale="utc").mjd
    if any(not isinstance(v, (int, float)) or not math.isfinite(v) or not 0 <= v < cutoff
           for v in values.values()):
        raise ValueError("History contains invalid or future timestamps")
    record_pairs(values, db)
    return source


def prepare_iers(start, end):
    """Try current Earth-orientation data; surface any degraded precision."""
    import numpy as np
    from astropy.time import Time
    from astropy.utils import iers

    notes = []
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        iers.conf.auto_download = True
        iers.conf.auto_max_age = 30
        iers.conf.remote_timeout = 15
        table = iers.IERS_Auto.open()
        _, status = table.ut1_utc(Time([start.mjd, end.mjd], format="mjd"), return_status=True)
        iers.earth_orientation_table.set(table)
        if np.any(np.asarray(status) < 0):
            notes.append("Earth-orientation tables do not cover the full night; precision is degraded.")
        if caught:
            notes.append("Earth-orientation refresh reported warnings; review the run logs.")
            for warning in caught:
                print(f"IERS: {warning.message}", file=sys.stderr)
    return notes


def validate_plan(plan, start, end):
    import numpy as np
    import pandas as pd
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astropy import units as u
    from ls4.config import LS4, global_constraints
    from ls4.scheduler import science_pointings, validate_pairs, validate_slews

    pairs = validate_pairs(plan)
    validate_slews(plan)
    science = science_pointings(plan)
    if not pairs or not np.isfinite(science[["ra", "dec", "start_time_mjd"]].to_numpy(dtype=float)).all():
        raise ValueError("Plan is empty or has nonfinite pointings")
    if (science["ra"].lt(0) | science["ra"].ge(360) | science["dec"].abs().gt(90)).any():
        raise ValueError("Pointing coordinates are outside their valid ranges")
    starts = Time(pd.to_datetime(science["start time (UTC)"]).to_numpy(), scale="utc")
    ends = Time(pd.to_datetime(science["end time (UTC)"]).to_numpy(), scale="utc")
    if np.any(starts.mjd < start.mjd - 1e-8) or np.any(ends.mjd > end.mjd + 1e-8):
        raise ValueError("Pointings fall outside astronomical twilight")
    coords = SkyCoord(science["ra"].to_numpy() * u.deg, science["dec"].to_numpy() * u.deg)
    for constraint in global_constraints:
        for times in (starts, ends):
            if not np.all(constraint(LS4, coords, times=times, grid_times_targets=False)):
                raise ValueError(f"Plan violates {type(constraint).__name__}")
    return {"images": len(science), "complete_pairs": pairs}


def render_report(manifest):
    from ls4.weather import forecast_summary

    weather = manifest["weather"]
    lines = [f"# LS4 observing night {manifest['evening_date']}", "",
        f"Generated UTC: {manifest['generated_at_utc']}", "",
        f"Observing window (UTC): {manifest['start_utc']} to {manifest['end_utc']}",
        f"Observing window (Chile): {manifest['start_local']} to {manifest['end_local']}", "",
        f"Validated exposures: {manifest['images']}; complete pairs: {manifest['complete_pairs']}.", "",
        f"Telescope file: `{manifest['obsplan']}` (following morning's date).",
        "", "## Cadence history", "", HISTORY_NOTE,
        f"History input: {manifest['history_input']}", "",
        "## La Silla weather", "", forecast_summary(weather), "",
        "Forecasts are advisory. Telescope weather limits and site operations remain authoritative.", ""]
    if weather["status"] == "available":
        location = weather["requested_location"]
        lines.extend([weather["attribution"], "",
            f"Retrieved UTC: {weather['retrieved_at_utc']}. Model issue time is not supplied by this API.",
            f"Site: {location['latitude']:.5f}°, {location['longitude']:.5f}°, {location['elevation_m']:.0f} m.",
            "Hourly UTC samples bracket astronomical twilight; maxima include these boundary hours.", "",
            "| UTC | Temperature °C | Cloud % | Rain probability % | Humidity % | Wind m/s | Gust m/s |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"])
        for row in weather["hourly"]:
            lines.append(f"| {row['utc']} | " + " | ".join(f"{row[key]:.1f}" for key in (
                "temperature_2m", "cloud_cover", "precipitation_probability",
                "relative_humidity_2m", "wind_speed_10m", "wind_gusts_10m")) + " |")
        lines.extend(["", f"[Forecast request]({weather['source_url']})", ""])
    if manifest["warnings"]:
        lines.extend(["## Items requiring attention", ""])
        lines.extend("- " + warning for warning in manifest["warnings"])
    lines.extend(["", "## Provenance", "", f"Code: `{manifest['provenance']['commit']}`.",
        f"Grid SHA256: `{manifest['grid_sha256']}`.",
        f"Dependency lock SHA256: `{manifest['provenance']['dependency_lock_sha256']}`.",
        "`manifest.json` records file checksums, versions, warnings, and run identity.",
        "`history.json` is the cumulative scheduling snapshot used by subsequent nights.", ""])
    return "\n".join(lines)


def generate(day, directory, *, use_github=False, refresh_iers=True):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=False)
    manifest = {"schema_version": 1, "evening_date": day.isoformat(),
                "generated_at_utc": datetime.now(timezone.utc).isoformat(),
                "generation": "failed", "stage": "setup", "warnings": [],
                "history_source": "generated-plans", "history_note": HISTORY_NOTE}
    try:
        from ls4.animate import make_gif
        from ls4.config import GRID_CSV, LS4
        from ls4.export import write_obsplan
        from ls4.grid import fresh_db, load_history
        from ls4.scheduler import night_window, plan_night
        from ls4.weather import fetch_forecast

        manifest["provenance"] = provenance()
        manifest["grid_sha256"] = sha256(GRID_CSV)
        start, end = night_window(day)
        if refresh_iers:
            manifest["warnings"].extend(prepare_iers(start, end))
            start, end = night_window(day)
        start_dt = start.to_datetime(timezone=timezone.utc)
        end_dt = end.to_datetime(timezone=timezone.utc)
        manifest.update(start_utc=start_dt.isoformat(), end_utc=end_dt.isoformat(),
                        start_local=start_dt.astimezone(ZoneInfo("America/Santiago")).isoformat(),
                        end_local=end_dt.astimezone(ZoneInfo("America/Santiago")).isoformat())
        csv_name = f"{day:%Y%m%d}.csv"
        obs_name = f"{day + timedelta(days=1):%Y%m%d}.obsplan"
        manifest.update(plan_csv=csv_name, obsplan=obs_name)
        with tempfile.TemporaryDirectory(prefix="ls4-nightly-") as temp:
            db = fresh_db(Path(temp) / "grid.db")
            manifest["stage"] = "history"
            manifest["history_input"] = restore_history(day, db, GRID_CSV, use_github=use_github)
            manifest["stage"] = "scheduling"
            plan = plan_night(start, end, directory / csv_name, db_path=db)
            manifest["stage"] = "validation"
            manifest.update(validate_plan(plan, start, end))
            write_obsplan(directory / csv_name, directory / obs_name)
            if len((directory / obs_name).read_text().splitlines()) != manifest["images"]:
                raise ValueError("Exported exposure count differs from the validated CSV")
            history = load_history(db)
            write_json(directory / "history.json", {
                "schema_version": 1, "source": "generated-plans", "note": HISTORY_NOTE,
                "through_night": day.isoformat(), "grid_sha256": manifest["grid_sha256"],
                "fields": dict(zip(history["base_field"], history["last_scheduled_mjd"].astype(float))),
            })
        manifest["stage"] = "weather"
        try:
            manifest["weather"] = fetch_forecast(start_dt, end_dt, {
                "latitude": float(LS4.location.lat.deg), "longitude": float(LS4.location.lon.deg),
                "elevation_m": float(LS4.location.height.to_value("m")),
            })
        except Exception as exc:
            manifest["weather"] = {"status": "unavailable", "error": str(exc)}
            manifest["warnings"].append("Weather forecast unavailable; inspect conditions before use.")
        write_json(directory / "weather.json", manifest["weather"])
        manifest["stage"] = "gif"
        try:
            gif = make_gif(directory / csv_name, directory / "night.gif", width=1280, height=800)
            from PIL import Image

            with Image.open(gif) as image:
                if image.format != "GIF" or image.n_frames < 2:
                    raise ValueError("Nightly GIF has no animation")
                image.verify()
            manifest["gif"] = "night.gif"
        except Exception as exc:
            (directory / "night.gif").unlink(missing_ok=True)
            manifest["warnings"].append(f"Nightly GIF unavailable ({type(exc).__name__}).")
        manifest.update(generation="success", stage="complete")
        (directory / "REPORT.md").write_text(render_report(manifest))
        manifest["files"] = {p.name: sha256(p) for p in sorted(directory.iterdir()) if p.is_file()}
    except Exception as exc:
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        write_json(directory / "manifest.json", manifest)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("generate", help="Generate an isolated nightly package")
    p.add_argument("--date", help="Local evening YYYY-MM-DD; default tonight")
    p.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "nightly")
    p.add_argument("--github-history", action="store_true", help="Restore history from published releases")
    p = sub.add_parser("publish", help="Publish a validated package as an immutable revision")
    p.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "nightly")
    p = sub.add_parser("notify", help="Report generation and publication status to Slack")
    p.add_argument("--output-dir", type=Path, default=ROOT / "outputs" / "nightly")
    p.add_argument("--job-status", default="failure")
    p.add_argument("--date")
    p = sub.add_parser("check", help="Alert Slack if today's published plan is missing or incomplete")
    p.add_argument("--date")
    args = parser.parse_args(argv)
    if args.command == "generate":
        generate(evening_date(args.date), args.output_dir, use_github=args.github_history)
    elif args.command == "publish":
        manifest = json.loads((args.output_dir / "manifest.json").read_text())
        day = date.fromisoformat(manifest["evening_date"])
        tag = f"ls4-{day:%Y%m%d}-{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
        result = github_client().publish(args.output_dir, tag, os.environ["GITHUB_SHA"])
        write_json(args.output_dir / "publication.json", result)
        if os.environ.get("GITHUB_STEP_SUMMARY"):
            with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as stream:
                stream.write(f"[Published nightly package]({result['url']})\n\n")
                stream.write((args.output_dir / "REPORT.md").read_text())
    elif args.command == "notify":
        manifest_path = args.output_dir / "manifest.json"
        manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {
            # Preserve an invalid manual date as diagnostic text instead of
            # making the failure notifier fail on the same date parse again.
            "evening_date": args.date or evening_date().isoformat(),
            "stage": "setup or dependency installation"}
        pub_path = args.output_dir / "publication.json"
        publication = json.loads(pub_path.read_text()) if pub_path.exists() else None
        options = dict(job_status=args.job_status, run_url=os.environ.get("LS4_RUN_URL", ""),
                       alert_user=os.environ.get("SLACK_ALERT_USER"))
        gif = args.output_dir / manifest["gif"] if manifest.get("gif") else None
        Slack(os.environ.get("SLACK_BOT_TOKEN", ""), os.environ["SLACK_CHANNEL_ID"]).send(
            slack_message(manifest, publication, **options), gif,
            blocks=slack_blocks(manifest, publication, **options))
    elif args.command == "check":
        slack = Slack(os.environ.get("SLACK_BOT_TOKEN", ""), os.environ["SLACK_CHANNEL_ID"])
        day = evening_date(args.date)
        release = None
        try:
            client = github_client()
            release = client.latest_nightly(day, before=False)
            if release:
                manifest = json.loads(client.read_asset(release, "manifest.json"))
                if manifest.get("generation") == "success" and not manifest.get("warnings"):
                    print(f"Published plan exists for {day}: {release['html_url']}")
                    return
            problem = ("No complete nightly package is available, or it requires attention. "
                       "Check the nightly workflow and generate/review the plan manually.")
        except Exception as exc:
            problem = f"Unable to verify the nightly package ({type(exc).__name__})."
        package = release["html_url"] if release else None
        links = [(package, "Available package"),
                 (os.environ["LS4_RUN_URL"], "Deadline-check logs")]
        problem += mention(os.environ.get("SLACK_ALERT_USER"))
        text = f":alarm_clock: *LS4 deadline check — {day}*\n{problem}" + "".join(
            f"\n<{url}|{label}>" for url, label in links if url)
        slack.send(text, blocks=[header(f":alarm_clock: LS4 deadline check — {day}"),
                                 section(":rotating_light: " + problem), *footer_links(links)])
        raise RuntimeError("Nightly deadline check needs attention")


if __name__ == "__main__":
    main()
