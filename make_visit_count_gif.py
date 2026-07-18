#!/usr/bin/env python3
"""Create a rotating Mollweide GIF of cumulative LS4 field visit counts.

The selected observing plan is sampled through time. Field-of-view outlines
are colored with the Viridis scale according to the number of visits to each
base field; a dithered visit counts as another visit to its original field.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from PIL import Image, ImageDraw

from make_observing_plan_gif import (
    MARGIN,
    Pointing,
    draw_fov,
    draw_graticule,
    draw_la_silla_marker,
    frame_times,
    load_font,
    local_sidereal_time,
    mjd_to_utc,
    project,
    read_grid,
    read_plan,
)


VIRIDIS_STOPS = (
    (68, 1, 84),
    (59, 82, 139),
    (33, 145, 140),
    (94, 201, 98),
    (253, 231, 37),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path, help="observing-plan CSV to animate")
    parser.add_argument("--output", type=Path,
                        help="GIF to create (default: <plan filename>_visits.gif)")
    parser.add_argument("--field-grid", type=Path, default=Path("assets/LS4_field_grid.csv"))
    parser.add_argument("--cadence-minutes", type=int, default=15)
    parser.add_argument("--duration", type=int, default=250)
    parser.add_argument("--width", type=int, default=1920)
    parser.add_argument("--height", type=int, default=1200)
    return parser.parse_args()


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
               frame_number: int, frame_count: int, max_visits: int, width: int, height: int) -> Image.Image:
    image = Image.new("RGB", (width, height), "#101820")
    draw = ImageDraw.Draw(image)
    title_font = load_font(30)
    time_font = load_font(40)
    lst_font = load_font(30)
    label_font = load_font(21)
    timestamp = mjd_to_utc(timestamp_mjd)
    center_ra = local_sidereal_time(timestamp_mjd)

    draw_graticule(draw, center_ra, width, height, label_font)
    for ra, dec in grid:
        x, y = project(ra, dec, center_ra, width, height)
        draw.point((x, y), fill="#526473")

    visits = Counter(base_field_name(pointing) for pointing in observed)
    for pointing in observed:
        draw_fov(draw, pointing, center_ra, width, height,
                 viridis_color(visits[base_field_name(pointing)], max_visits), 2)
    draw_la_silla_marker(draw, center_ra, width, height, label_font)

    draw.text((MARGIN, 14), f"LS4 field visit count   frame {frame_number}/{frame_count}",
              fill="#f2f6f9", font=title_font)
    draw.text((MARGIN, 49), f"{timestamp:%Y-%m-%d %H:%M UTC}   MJD {timestamp_mjd:.5f}",
              fill="#f2f6f9", font=time_font)
    draw.text((width - 465, 14), f"La Silla LST {center_ra / 15:04.1f} h", fill="#d4dce3", font=lst_font)
    draw_colorbar(draw, max_visits, width, height)
    return image


def main() -> None:
    args = parse_args()
    if args.cadence_minutes <= 0 or args.duration <= 0:
        raise ValueError("cadence-minutes and duration must be positive")
    if not args.plan.is_file() or not args.field_grid.is_file():
        raise FileNotFoundError("plan CSV or field-grid CSV was not found")

    plan = read_plan(args.plan)
    grid = read_grid(args.field_grid)
    timeline = frame_times(plan, args.cadence_minutes)
    max_visits = max(Counter(base_field_name(pointing) for pointing in plan).values())

    frames = []
    for frame_number, timestamp_mjd in enumerate(timeline, start=1):
        observed = [pointing for pointing in plan if pointing.end_mjd <= timestamp_mjd]
        frames.append(make_frame(grid, observed, timestamp_mjd, frame_number, len(timeline), max_visits,
                                 args.width, args.height))
        print(f"Added frame {frame_number}/{len(timeline)}: {mjd_to_utc(timestamp_mjd):%Y-%m-%d %H:%M UTC}")

    output = args.output or Path(f"{args.plan.stem}_visits.gif")
    output.parent.mkdir(parents=True, exist_ok=True)
    frames[0].save(output, save_all=True, append_images=frames[1:], duration=args.duration, loop=0,
                   optimize=False)
    print(f"Wrote {output} ({len(frames)} frames; maximum visits: {max_visits})")


if __name__ == "__main__":
    main()
