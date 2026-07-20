#!/usr/bin/env python3
"""Create a Mollweide GIF showing one LS4 observing plan as Earth rotates.

Frames are sampled through the selected plan.  The map is centered on La Silla's local
sidereal time at the frame timestamp, so the celestial sphere moves overhead.
Every scheduled pointing is drawn as an LS4 field-of-view outline: first
pointings are cyan and dithered revisits are magenta. The Milky Way's galactic
plane is drawn in the background so you can see at a glance whether the survey
is tracking through/around it as expected.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


# LS4 rectangular field of view from constants.py, converted from arcseconds.
FOV_RA_DEG = 16924.67 / 3600
FOV_DEC_DEG = 17826.13 / 3600
LA_SILLA_LONGITUDE_DEG = -70.7367
LA_SILLA_LATITUDE_DEG = -29.2612
MARGIN = 70
MJD_UNIX_EPOCH = 40587.0

# J2000 galactic-pole/node parameters used to convert galactic (l, b) to
# equatorial (RA, Dec) without requiring an extra astronomy dependency.
GALACTIC_NGP_RA_DEG = 192.859508
GALACTIC_NGP_DEC_DEG = 27.128336
GALACTIC_NCP_LON_DEG = 122.932


@dataclass(frozen=True)
class Pointing:
    target: str
    start_mjd: float
    end_mjd: float
    ra: float
    dec: float
    dither: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path, help="observing-plan CSV to animate")
    parser.add_argument("--output", type=Path,
                        help="GIF to create (default: <plan filename>.gif in the current directory)")
    parser.add_argument("--field-grid", type=Path, default=Path("assets/LS4_field_grid.csv"))
    parser.add_argument("--cadence-minutes", type=int, default=15,
                        help="time between frames within each plan (default: 15)")
    parser.add_argument("--duration", type=int, default=250,
                        help="milliseconds per GIF frame (default: 250)")
    parser.add_argument("--width", type=int, default=1920)
    parser.add_argument("--height", type=int, default=1200)
    parser.add_argument("--no-galaxy", action="store_true",
                        help="omit the galactic-plane overlay")
    return parser.parse_args()


def read_grid(path: Path) -> list[tuple[float, float]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        return [(float(row["ra_deg"]), float(row["dec_deg"])) for row in reader]


def utc_to_mjd(timestamp: datetime) -> float:
    """Convert a naive UTC datetime to modified Julian date."""
    unix_epoch = datetime(1970, 1, 1)
    return MJD_UNIX_EPOCH + (timestamp - unix_epoch).total_seconds() / 86400.0


def mjd_to_utc(mjd: float) -> datetime:
    """Convert modified Julian date to a naive UTC datetime."""
    return datetime(1970, 1, 1) + timedelta(days=mjd - MJD_UNIX_EPOCH)


def read_plan(path: Path) -> list[Pointing]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"target", "start time (UTC)", "end time (UTC)", "ra", "dec", "configuration"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        result = []
        for row in reader:
            if row["target"] == "TransitionBlock":
                continue
            try:
                start_utc = datetime.fromisoformat(row["start time (UTC)"])
                end_utc = datetime.fromisoformat(row["end time (UTC)"])
                start_mjd = float(row["start_time_mjd"]) if row.get("start_time_mjd") else utc_to_mjd(start_utc)
                result.append(Pointing(
                    target=row["target"],
                    start_mjd=start_mjd,
                    end_mjd=start_mjd + (end_utc - start_utc).total_seconds() / 86400.0,
                    ra=float(row["ra"]), dec=float(row["dec"]),
                    dither="True" in row["configuration"],
                ))
            except (TypeError, ValueError):
                continue
    if not result:
        raise ValueError(f"{path} contains no scheduled pointings")
    return sorted(result, key=lambda pointing: pointing.start_mjd)


def local_sidereal_time(mjd: float) -> float:
    """Approximate local sidereal time in degrees at La Silla from MJD UTC."""
    julian_date = mjd + 2400000.5
    gmst = 280.46061837 + 360.98564736629 * (julian_date - 2451545.0)
    return (gmst + LA_SILLA_LONGITUDE_DEG) % 360.0


def galactic_to_equatorial(l_deg: float, b_deg: float) -> tuple[float, float]:
    """Convert J2000 galactic (l, b) in degrees to equatorial (RA, Dec) in degrees."""
    l = math.radians(l_deg)
    b = math.radians(b_deg)
    ra_ngp = math.radians(GALACTIC_NGP_RA_DEG)
    dec_ngp = math.radians(GALACTIC_NGP_DEC_DEG)
    l_ncp = math.radians(GALACTIC_NCP_LON_DEG)

    sin_dec = math.sin(b) * math.sin(dec_ngp) + math.cos(b) * math.cos(dec_ngp) * math.cos(l_ncp - l)
    dec = math.asin(max(-1.0, min(1.0, sin_dec)))

    y = math.cos(b) * math.sin(l_ncp - l)
    x = math.cos(dec_ngp) * math.sin(b) - math.sin(dec_ngp) * math.cos(b) * math.cos(l_ncp - l)
    ra = math.degrees(ra_ngp + math.atan2(y, x))
    return ra % 360.0, math.degrees(dec)


def mollweide(ra: float, dec: float, center_ra: float) -> tuple[float, float]:
    """Return normalized Mollweide coordinates in the range x [-2, 2], y [-1, 1]."""
    longitude = math.radians(((ra - center_ra + 180) % 360) - 180)
    latitude = math.radians(dec)
    theta = latitude
    for _ in range(10):
        numerator = 2 * theta + math.sin(2 * theta) - math.pi * math.sin(latitude)
        denominator = 2 + 2 * math.cos(2 * theta)
        if abs(denominator) < 1e-10:
            break
        theta -= numerator / denominator
    return 2 * math.sqrt(2) / math.pi * longitude * math.cos(theta), math.sqrt(2) * math.sin(theta)


def project(ra: float, dec: float, center_ra: float, width: int, height: int) -> tuple[int, int]:
    x, y = mollweide(ra, dec, center_ra)
    plot_width, plot_height = width - 2 * MARGIN, height - 2 * MARGIN
    return round(MARGIN + (x + 2 * math.sqrt(2)) / (4 * math.sqrt(2)) * plot_width), round(
        MARGIN + (math.sqrt(2) - y) / (2 * math.sqrt(2)) * plot_height
    )


def fov_corners(pointing: Pointing) -> list[tuple[float, float]]:
    """Approximate the rectangular sky footprint using small-angle offsets."""
    half_ra = FOV_RA_DEG / 2 / max(math.cos(math.radians(pointing.dec)), 0.05)
    half_dec = FOV_DEC_DEG / 2
    return [
        ((pointing.ra - half_ra) % 360, pointing.dec - half_dec),
        ((pointing.ra + half_ra) % 360, pointing.dec - half_dec),
        ((pointing.ra + half_ra) % 360, pointing.dec + half_dec),
        ((pointing.ra - half_ra) % 360, pointing.dec + half_dec),
        ((pointing.ra - half_ra) % 360, pointing.dec - half_dec),
    ]


def draw_fov(draw: ImageDraw.ImageDraw, pointing: Pointing, center_ra: float, width: int, height: int,
             color: str, line_width: int) -> None:
    corners = fov_corners(pointing)
    projected = [project(ra, dec, center_ra, width, height) for ra, dec in corners]
    # Avoid drawing a false line across the Mollweide map's seam.
    for start, end in zip(projected, projected[1:]):
        if abs(start[0] - end[0]) < width / 2:
            draw.line((start, end), fill=color, width=line_width)


def draw_galactic_plane(draw: ImageDraw.ImageDraw, center_ra: float, width: int, height: int,
                        label_font: ImageFont.ImageFont) -> None:
    """Draw a soft band around the galactic plane plus its centerline and center marker."""
    band_color = "#2e2140"
    line_color = "#8a63c4"

    # Filled +-10 deg latitude band, built from small quads so we can skip any
    # quad that would otherwise be stretched across the Mollweide seam.
    l_step = 4
    longitudes = list(range(0, 361, l_step))
    for l0, l1 in zip(longitudes, longitudes[1:]):
        quad = []
        for l, b in ((l0, -10), (l1, -10), (l1, 10), (l0, 10)):
            ra, dec = galactic_to_equatorial(l, b)
            quad.append(project(ra, dec, center_ra, width, height))
        xs = [point[0] for point in quad]
        if max(xs) - min(xs) < width / 2:
            draw.polygon(quad, fill=band_color)

    # Solid centerline at b = 0.
    line_points = []
    for l in range(0, 361, 2):
        ra, dec = galactic_to_equatorial(l, 0)
        line_points.append(project(ra, dec, center_ra, width, height))
    for start, end in zip(line_points, line_points[1:]):
        if abs(start[0] - end[0]) < width / 2:
            draw.line((start, end), fill=line_color, width=2)

    # Galactic center marker (l=0, b=0).
    gc_ra, gc_dec = galactic_to_equatorial(0, 0)
    x, y = project(gc_ra, gc_dec, center_ra, width, height)
    radius = 5
    draw.ellipse((x - radius, y - radius, x + radius, y + radius), fill=line_color)
    draw.text((x + 10, y - 10), "Galactic center", fill=line_color, font=label_font)


def load_font(size: int) -> ImageFont.ImageFont:
    """Use a scalable font when Pillow provides one, with a safe fallback."""
    for font_path in (
        "/System/Library/Fonts/Supplemental/Verdana.ttf",
        "DejaVuSans.ttf",
    ):
        try:
            return ImageFont.truetype(font_path, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def draw_graticule(draw: ImageDraw.ImageDraw, center_ra: float, width: int, height: int,
                   label_font: ImageFont.ImageFont) -> None:
    for dec in range(-60, 61, 30):
        coordinates = [project(ra, dec, center_ra, width, height) for ra in range(0, 361, 3)]
        draw.line(coordinates, fill="#314252", width=1)
    for ra in range(0, 360, 30):
        coordinates = [project(ra, dec, center_ra, width, height) for dec in range(-89, 90)]
        draw.line(coordinates, fill="#314252", width=1)
    for ra in range(0, 360, 30):
        x, y = project(ra, 0, center_ra, width, height)
        draw.text((x - 12, y + 7), f"{ra // 15:02d}h", fill="#a9b7c4", font=label_font)
    for dec in range(-60, 61, 30):
        _, y = project((center_ra + 179.5) % 360, dec, center_ra, width, height)
        draw.text((8, y - 10), f"{dec:+d}", fill="#a9b7c4", font=label_font)
    draw.text((MARGIN, height - MARGIN + 13), "RA", fill="#a9b7c4", font=label_font)
    draw.text((8, MARGIN), "Dec", fill="#a9b7c4", font=label_font)


def draw_la_silla_marker(draw: ImageDraw.ImageDraw, center_ra: float, width: int, height: int,
                          label_font: ImageFont.ImageFont) -> None:
    """Mark La Silla's zenith, the sky coordinate directly above the observatory."""
    x, y = project(center_ra, LA_SILLA_LATITUDE_DEG, center_ra, width, height)
    radius = 9
    draw.ellipse((x - radius, y - radius, x + radius, y + radius), outline="#ffd166", width=3)
    draw.line((x - radius - 5, y, x + radius + 5, y), fill="#ffd166", width=2)
    draw.line((x, y - radius - 5, x, y + radius + 5), fill="#ffd166", width=2)
    draw.text((x + 15, y - 12), "La Silla zenith", fill="#ffd166", font=label_font)


def make_frame(grid: list[tuple[float, float]], observed: list[Pointing], timestamp_mjd: float,
               frame_number: int, frame_count: int, width: int, height: int, show_galaxy: bool) -> Image.Image:
    image = Image.new("RGB", (width, height), "#101820")
    draw = ImageDraw.Draw(image)
    title_font = load_font(30)
    time_font = load_font(40)
    lst_font = load_font(30)
    label_font = load_font(21)
    timestamp = mjd_to_utc(timestamp_mjd)
    center_ra = local_sidereal_time(timestamp_mjd)

    if show_galaxy:
        draw_galactic_plane(draw, center_ra, width, height, label_font)

    draw_graticule(draw, center_ra, width, height, label_font)

    for ra, dec in grid:
        x, y = project(ra, dec, center_ra, width, height)
        draw.point((x, y), fill="#526473")

    for pointing in observed:
        draw_fov(draw, pointing, center_ra, width, height,
                 "#e652a0" if pointing.dither else "#44b8e8", 2)
    draw_la_silla_marker(draw, center_ra, width, height, label_font)

    draw.text((MARGIN, 14), f"LS4 observing plan   frame {frame_number}/{frame_count}",
              fill="#f2f6f9", font=title_font)
    draw.text((MARGIN, 49), f"{timestamp:%Y-%m-%d %H:%M UTC}   MJD {timestamp_mjd:.5f}",
              fill="#f2f6f9", font=time_font)
    caption = "cyan: first pointing   magenta: dithered revisit   yellow: La Silla zenith"
    if show_galaxy:
        caption += "   purple: galactic plane (+-10 deg)"
    draw.text((MARGIN, height - 28), caption, fill="#d4dce3", font=label_font)
    draw.text((width - 465, 14), f"La Silla LST {center_ra / 15:04.1f} h", fill="#d4dce3", font=lst_font)
    return image


def frame_times(plan: list[Pointing], cadence_minutes: int) -> list[float]:
    times = []
    # Start after the first exposure so every GIF frame includes a footprint.
    timestamp = plan[0].end_mjd
    final_time = plan[-1].end_mjd
    cadence = cadence_minutes / 1440.0
    while timestamp < final_time:
        times.append(timestamp)
        timestamp += cadence
    if not times or times[-1] != final_time:
        times.append(final_time)
    return times


def main() -> None:
    args = parse_args()
    if args.cadence_minutes <= 0 or args.duration <= 0:
        raise ValueError("cadence-minutes and duration must be positive")
    if not args.plan.is_file() or not args.field_grid.is_file():
        raise FileNotFoundError("plan CSV or field-grid CSV was not found")

    plan = read_plan(args.plan)
    grid = read_grid(args.field_grid)
    timeline = frame_times(plan, args.cadence_minutes)

    frames = []
    for frame_number, timestamp_mjd in enumerate(timeline, start=1):
        current = [pointing for pointing in plan if pointing.end_mjd <= timestamp_mjd]
        frames.append(make_frame(grid, current, timestamp_mjd, frame_number, len(timeline),
                                  args.width, args.height, not args.no_galaxy))
        print(f"Added frame {frame_number}/{len(timeline)}: {mjd_to_utc(timestamp_mjd):%Y-%m-%d %H:%M UTC}")

    output = args.output or Path(f"{args.plan.stem}.gif")
    output.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(output, save_all=True, append_images=frames[1:], duration=args.duration, loop=0,
                   optimize=False)
    print(f"Wrote {output} ({len(frames)} frames)")


if __name__ == "__main__":
    main()