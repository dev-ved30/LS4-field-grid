"""Mollweide GIFs of an observing plan as the sky rotates over La Silla.

Frames are sampled through the plan. The map is centered on La Silla's local
sidereal time at the frame timestamp, so the celestial sphere moves overhead.
Each scheduled pointing is drawn as an LS4 field-of-view outline, colored either
by pointing type (``mode="pointings"``: cyan first pointing, magenta dither) or
by cumulative visits to its base field (``mode="visits"``, Viridis). The
galactic plane is drawn in the background.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

from ls4.config import FOV_length, FOV_width, GRID_CSV, LS4

# LS4 rectangular field of view, in degrees.
FOV_RA_DEG = float(FOV_width.to_value("deg"))
FOV_DEC_DEG = float(FOV_length.to_value("deg"))
LA_SILLA_LONGITUDE_DEG = float(LS4.location.lon.deg)
LA_SILLA_LATITUDE_DEG = float(LS4.location.lat.deg)
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



VIRIDIS_STOPS = (
    (68, 1, 84),
    (59, 82, 139),
    (33, 145, 140),
    (94, 201, 98),
    (253, 231, 37),
)


def base_field_name(pointing: Pointing) -> str:
    return pointing.target.removesuffix("_dither")


def viridis_color(value: float, maximum: int) -> tuple[int, int, int]:
    """Interpolate a Viridis color from the integer visit count."""
    fraction = 0.0 if maximum <= 1 else (value - 1) / (maximum - 1)
    fraction = max(0.0, min(1.0, fraction))
    position = fraction * (len(VIRIDIS_STOPS) - 1)
    index = min(int(position), len(VIRIDIS_STOPS) - 2)
    remainder = position - index
    start, end = VIRIDIS_STOPS[index], VIRIDIS_STOPS[index + 1]
    return tuple(round(a + (b - a) * remainder) for a, b in zip(start, end))


def draw_colorbar(draw: ImageDraw.ImageDraw, maximum: int, width: int, height: int) -> None:
    font = load_font(20)
    x0, x1 = width - 430, width - MARGIN
    y0, y1 = 110, 134
    for x in range(x0, x1):
        fraction = (x - x0) / max(1, x1 - x0 - 1)
        color = viridis_color(1 + fraction * (maximum - 1), maximum)
        draw.line((x, y0, x, y1), fill=color)
    draw.rectangle((x0, y0, x1, y1), outline="#d4dce3", width=1)
    draw.text((x0, y0 - 27), "Cumulative visits", fill="#d4dce3", font=font)
    draw.text((x0, y1 + 4), "1", fill="#d4dce3", font=font)
    draw.text((x1 - 32, y1 + 4), str(maximum), fill="#d4dce3", font=font)




def make_frame(grid: list[tuple[float, float]], observed: list[Pointing], timestamp_mjd: float,
               frame_number: int, frame_count: int, width: int, height: int, show_galaxy: bool,
               mode: str = "pointings", max_visits: int = 1) -> Image.Image:
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

    if mode == "visits":
        visits = Counter(base_field_name(pointing) for pointing in observed)
        for pointing in observed:
            draw_fov(draw, pointing, center_ra, width, height,
                     viridis_color(visits[base_field_name(pointing)], max_visits), 2)
    else:
        for pointing in observed:
            draw_fov(draw, pointing, center_ra, width, height,
                     "#e652a0" if pointing.dither else "#44b8e8", 2)
    draw_la_silla_marker(draw, center_ra, width, height, label_font)

    title = "LS4 field visit count" if mode == "visits" else "LS4 observing plan"
    draw.text((MARGIN, 14), f"{title}   frame {frame_number}/{frame_count}",
              fill="#f2f6f9", font=title_font)
    draw.text((MARGIN, 49), f"{timestamp:%Y-%m-%d %H:%M UTC}   MJD {timestamp_mjd:.5f}",
              fill="#f2f6f9", font=time_font)
    if mode == "visits":
        draw.text((width - 465, 14), f"La Silla LST {center_ra / 15:04.1f} h", fill="#d4dce3", font=lst_font)
        draw_colorbar(draw, max_visits, width, height)
    else:
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



def make_gif(plan_path: Path, output: Path, mode: str = "pointings", grid_path: Path = GRID_CSV,
             cadence_minutes: int = 15, duration: int = 250, width: int = 1920, height: int = 1200,
             show_galaxy: bool = True) -> Path:
    """Render one plan as a GIF. `mode` is "pointings" or "visits"."""
    if mode not in ("pointings", "visits"):
        raise ValueError(f"Unknown GIF mode: {mode}")
    if cadence_minutes <= 0 or duration <= 0:
        raise ValueError("cadence-minutes and duration must be positive")
    plan = read_plan(Path(plan_path))
    grid = read_grid(Path(grid_path))
    timeline = frame_times(plan, cadence_minutes)
    max_visits = max(Counter(base_field_name(pointing) for pointing in plan).values())

    frames = []
    for frame_number, timestamp_mjd in enumerate(timeline, start=1):
        observed = [pointing for pointing in plan if pointing.end_mjd <= timestamp_mjd]
        frames.append(make_frame(grid, observed, timestamp_mjd, frame_number, len(timeline),
                                 width, height, show_galaxy, mode, max_visits))
        print(f"Added frame {frame_number}/{len(timeline)}: {mjd_to_utc(timestamp_mjd):%Y-%m-%d %H:%M UTC}")

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(output, save_all=True, append_images=frames[1:], duration=duration, loop=0,
                   optimize=False)
    print(f"Wrote {output} ({len(frames)} frames)")
    return output
