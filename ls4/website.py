"""Static website of nightly plans and survey coverage: `ls4 site`, deployed to GitHub Pages.

The site groups nights into collections: "production" (top-level `plans/*.csv` plus,
with `releases=True`, the newest published revision of every nightly release) and one
per `plans/forecast_*/` directory. For each collection it writes

- `data/<collection>/<YYYYMMDD>.json`: the night's pointings, metrics, Moon and sidereal time;
- `data/<collection>/survey.json`: per-band HEALPix coverage by night (as `ls4 bands`);
- `files/<collection>/`: plan CSVs and telescope `.obsplan` files for download;

plus `data/index.json` and the page itself (`ls4/web/`). Serve the directory over HTTP
(`ls4 site --serve`); browsers do not let a page opened from disk read its JSON files.
"""

import contextlib
from datetime import datetime, timedelta, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess

import healpy as hp
import numpy as np
import pandas as pd

from ls4.config import (BANDS, FOV_length, FOV_width, LS4, PLANS_DIR, QUADRANT_BANDS, ROOT,
                        SOLE, global_constraints, max_airmass)

WEB = Path(__file__).resolve().parent / "web"
REPOSITORY = "dev-ved30/LS4-field-grid"
HISTORY_NOTE = ("Cadence and coverage assume every plan was executed as scheduled; "
                "weather losses and deviations are not reconciled.")


def stamp_date(stamp):
    return datetime.strptime(stamp, "%Y%m%d").date()


def obsplan_name(stamp):
    """Telescope files are named for the following morning."""
    return f"{stamp_date(stamp) + timedelta(days=1):%Y%m%d}.obsplan"


def is_stamp(name):
    try:
        stamp_date(name)
    except ValueError:
        return False
    return len(name) == 8


def readme_summary(directory):
    """First paragraph after the title of a directory's README, if any."""
    readme = Path(directory) / "README.md"
    if not readme.exists():
        return ""
    paragraphs = [p.strip() for p in readme.read_text().split("\n\n") if p.strip()]
    body = [p for p in paragraphs if not p.startswith("#")]
    return " ".join(body[0].replace("`", "").split()) if body else ""


# ---- sources ---------------------------------------------------------------

def repository_collections(plans_dir=PLANS_DIR):
    """Production plans saved in the repository, then one collection per forecast directory."""
    plans_dir = Path(plans_dir)
    collections = [{
        "key": "production", "label": "Production",
        "description": "Plans scheduled for the telescope: published nightly releases and plans "
                       "saved in the repository.",
        "nights": {p.stem: {"csv": p} for p in sorted(plans_dir.glob("*.csv")) if is_stamp(p.stem)},
    }]
    for directory in sorted(plans_dir.glob("forecast_*")):
        nights = {}
        for path in sorted(directory.glob("*.csv")):
            if is_stamp(path.stem):
                obsplan = directory / "obsplans" / obsplan_name(path.stem)
                nights[path.stem] = {"csv": path, **({"obsplan": obsplan} if obsplan.exists() else {})}
        if nights:
            start = directory.name.removeprefix("forecast_")
            label = f"Forecast {stamp_date(start):%Y-%m-%d}" if is_stamp(start) else directory.name
            collections.append({"key": directory.name, "label": label,
                                "description": readme_summary(directory), "nights": nights})
    return collections


def add_releases(production, download_dir, repository=REPOSITORY, token=None):
    """Add the newest published revision of every nightly release (it replaces a saved plan)."""
    from ls4.automation import GitHub

    client = GitHub(repository, token)
    download_dir = Path(download_dir)
    download_dir.mkdir(parents=True, exist_ok=True)
    releases = client.nightly_releases()
    for stamp, release in sorted(releases.items()):
        manifest = json.loads(client.read_asset(release, "manifest.json"))
        if manifest.get("generation") != "success":
            continue
        night = {}
        for key, name in (("csv", manifest["plan_csv"]), ("obsplan", manifest["obsplan"])):
            data = client.read_asset(release, name)
            if hashlib.sha256(data).hexdigest() != manifest["files"][name]:
                raise ValueError(f"Checksum mismatch for {name} in {release['tag_name']}")
            path = download_dir / name
            path.write_bytes(data)
            night[key] = path
        assets = {a["name"]: a["browser_download_url"] for a in release["assets"]}
        weather = manifest.get("weather") or {}
        night["release"] = {
            "tag": release["tag_name"], "url": release["html_url"],
            "generated_utc": manifest.get("generated_at_utc"),
            "warnings": manifest.get("warnings", []),
            "gif": assets.get("night.gif"), "report": assets.get("REPORT.md"),
            "weather": weather.get("summary") if weather.get("status") == "available" else None,
        }
        production["nights"][stamp] = night
    print(f"  {len(releases)} published nights from {repository}")


# ---- per night -------------------------------------------------------------

def _plain(value):
    value = value.item() if hasattr(value, "item") else value
    return round(value, 3) if isinstance(value, float) else value


def sky_samples(start, end, minutes=15):
    """Local sidereal time and Moon position through the night, for the zenith and Moon overlays."""
    from astropy import units as u
    from astropy.coordinates import get_body

    count = int(np.ceil((end - start).to_value(u.minute) / minutes)) + 1
    times = start + np.linspace(0, (end - start).to_value(u.day), count) * u.day
    lst = LS4.local_sidereal_time(times).deg
    moon = get_body("moon", times, LS4.location)
    return [[round(float(t), 6), round(float(s), 3), round(float(r), 3), round(float(d), 3)]
            for t, s, r, d in zip(times.mjd, lst, moon.ra.deg, moon.dec.deg)]


def night_data(stamp, plan, night_info, files, sole_fields, area):
    from astropy.time import Time
    from ls4.metrics import night_metrics
    from ls4.scheduler import night_window, science_pointings

    start, end = night_window(stamp_date(stamp))
    science = science_pointings(plan).dropna(subset=["ra", "dec"]).sort_values("start_time_mjd")
    ends = Time(pd.to_datetime(science["end time (UTC)"]).to_numpy(), scale="utc").mjd
    phase = int(science["configuration"].astype(str).str.contains("'phase': 1").any())
    middle = start + (end - start) / 2
    data = {
        "night": stamp, "date": stamp_date(stamp).isoformat(),
        "window": [round(start.mjd, 6), round(end.mjd, 6)],
        "phase": phase if science["configuration"].astype(str).str.contains("phase").any() else None,
        "metrics": {key: _plain(value) for key, value in night_metrics(plan).items()},
        "area": area,
        "moonIllumination": round(float(LS4.moon_illumination(middle)), 3),
        "sky": sky_samples(start, end),
        "targets": science["base_field"].tolist(),
        # [start MJD, end MJD, RA, Dec, dither, block, SOLE]
        "pointings": [[round(float(s), 6), round(float(e), 6), round(float(r), 4), round(float(d), 4),
                       int(dither), int(block), int(field in sole_fields)]
                      for s, e, r, d, dither, block, field in zip(
                          science["start_time_mjd"], ends, science["ra"], science["dec"],
                          science["is_dither"], science["block_number"], science["base_field"])],
        "files": files,
    }
    if "release" in night_info:
        data["release"] = night_info["release"]
    return data


# ---- per collection --------------------------------------------------------

def build_collection(collection, out, nside):
    from ls4.bands import band_map_data, band_visits, nightly_area, pixel_stats
    from ls4.grid import load_grid_csv, program_ids
    from ls4.export import write_obsplan
    from ls4.metrics import load_plans

    key = collection["key"]
    files_dir = out / "files" / key
    data_dir = out / "data" / key
    files_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    order = sorted(collection["nights"])
    print(f"{collection['label']}: {len(order)} nights ({order[0]} to {order[-1]})")

    # Copy plans and telescope files under the night's name; `load_plans` names nights by file stem.
    csvs = {}
    for stamp in order:
        night = collection["nights"][stamp]
        csvs[stamp] = files_dir / f"{stamp}.csv"
        if Path(night["csv"]).resolve() != csvs[stamp].resolve():
            shutil.copyfile(night["csv"], csvs[stamp])
        obsplan = files_dir / obsplan_name(stamp)
        if night.get("obsplan"):
            if Path(night["obsplan"]).resolve() != obsplan.resolve():
                shutil.copyfile(night["obsplan"], obsplan)
        else:
            with contextlib.redirect_stdout(io.StringIO()):
                write_obsplan(csvs[stamp], obsplan)

    plans = load_plans([csvs[stamp] for stamp in order])
    stats, nightly = pixel_stats(band_visits(plans, nside), nside)
    pixel_area = hp.nside2pixarea(nside, degrees=True)
    any_band = nightly.groupby("night")["pixel"].nunique() * pixel_area
    areas = nightly_area(nightly, order, nside).set_index("night")

    grid = load_grid_csv()
    sole_fields = set(grid.loc[program_ids(grid["l"], grid["b"]) == SOLE, "base_field"])
    summaries = []
    for stamp in order:
        plan = pd.read_csv(csvs[stamp])
        area = {"any": round(float(any_band.get(stamp, 0)), 1),
                **{band: round(float(areas.loc[stamp, f"{band}_area_deg2"]), 1) for band in BANDS}}
        files = {"csv": f"files/{key}/{stamp}.csv", "obsplan": f"files/{key}/{obsplan_name(stamp)}"}
        night = night_data(stamp, plan, collection["nights"][stamp], files, sole_fields, area)
        write_compact(data_dir / f"{stamp}.json", night)
        summaries.append({"night": stamp, "date": night["date"], "phase": night["phase"],
                          "images": night["metrics"]["images"],
                          "pairs": night["metrics"]["complete_pairs"],
                          "solePairs": night["metrics"]["sole_pairs"], "area": area,
                          "published": "release" in night})

    survey = band_map_data(nightly, order, plans, nside)
    survey["cumulative"] = [{"night": stamp, **{band: round(float(areas.loc[stamp, f"{band}_cumulative_deg2"]), 1)
                                                for band in BANDS},
                             "all": round(float(areas.loc[stamp, "all_bands_cumulative_deg2"]), 1)}
                            for stamp in order]
    survey["nside"] = nside
    write_compact(data_dir / "survey.json", survey)
    return {"key": key, "label": collection["label"], "description": collection["description"],
            "nights": summaries}


# ---- site --------------------------------------------------------------------

def write_compact(path, data):
    Path(path).write_text(json.dumps(data, separators=(",", ":"), allow_nan=False))


def galactic_plane(step=2):
    from ls4.animate import galactic_to_equatorial

    return [[round(c, 3) for c in galactic_to_equatorial(l % 360, 0)] for l in range(0, 361, step)]


def code_commit():
    try:
        return os.environ.get("GITHUB_SHA") or subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def prepare_output(out):
    """Rebuild into a fresh directory, refusing to delete anything that is not a previous site."""
    out = Path(out)
    if out.exists():
        if any(out.iterdir()) and not (out / "data" / "index.json").exists():
            raise SystemExit(f"{out} exists and is not an LS4 site; choose another --output.")
        shutil.rmtree(out)
    out.mkdir(parents=True)
    return out


def build_site(output, *, plans_dir=PLANS_DIR, releases=False, repository=REPOSITORY,
               token=None, nside=64):
    from astroplan.constraints import MoonSeparationConstraint

    out = prepare_output(output)
    collections = repository_collections(plans_dir)
    if releases:
        add_releases(collections[0], out / "files" / "production", repository, token)
    collections = [c for c in collections if c["nights"]]
    if not collections:
        raise SystemExit("No plans found.")
    index = [build_collection(c, out, nside) for c in collections]

    decs = [p[3] for c in index for n in c["nights"]
            for p in json.loads((out / "data" / c["key"] / f"{n['night']}.json").read_text())["pointings"]]
    half = np.hypot(FOV_width.to_value("deg"), FOV_length.to_value("deg")) / 2 + 2
    moon_sep = next(c.min for c in global_constraints if isinstance(c, MoonSeparationConstraint))
    write_compact(out / "data" / "index.json", {
        "generatedUtc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "commit": code_commit(), "repository": repository,
        "historyNote": HISTORY_NOTE,
        "site": {"name": LS4.name, "lat": float(LS4.location.lat.deg), "lon": float(LS4.location.lon.deg)},
        "fov": {"ra": FOV_width.to_value("deg"), "dec": FOV_length.to_value("deg")},
        "quadrants": [[sx, sy, band] for (sx, sy), band in QUADRANT_BANDS.items()],
        "bands": list(BANDS),
        "airmassRadius": float(np.degrees(np.arccos(1 / max_airmass))),
        "maxAirmass": max_airmass,
        "moonSeparation": float(moon_sep.to_value("deg")),
        "decRange": [max(-90, 10 * np.floor((min(decs) - half) / 10)),
                     min(90, 10 * np.ceil((max(decs) + half) / 10))],
        "galactic": galactic_plane(),
        "collections": index,
    })
    for asset in WEB.iterdir():
        if asset.is_file():
            shutil.copyfile(asset, out / asset.name)
    (out / ".nojekyll").touch()
    print(f"Wrote site to {out}")
    return out


def serve(directory, port=8000):
    import functools
    import http.server

    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(directory))
    with http.server.ThreadingHTTPServer(("127.0.0.1", port), handler) as server:
        print(f"Serving {directory} at http://127.0.0.1:{port}/ (Ctrl-C to stop)")
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            pass
