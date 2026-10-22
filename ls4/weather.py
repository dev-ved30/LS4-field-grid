"""Hourly La Silla forecasts for a nightly report (advisory, never a safety gate)."""

from datetime import datetime, timedelta, timezone
import math
from urllib.parse import urlencode

from ls4.automation import request_json

VARIABLES = {
    "temperature_2m": "°C",
    "relative_humidity_2m": "%",
    "precipitation_probability": "%",
    "cloud_cover": "%",
    "wind_speed_10m": "m/s",
    "wind_gusts_10m": "m/s",
}


def fetch_forecast(start, end, location):
    """Fetch UTC hours bracketing the complete twilight window, at site elevation.

    The API exposes retrieval time, not the underlying models' issue times. Do not
    confuse its `generationtime_ms` (request duration) with a forecast timestamp.
    """
    first = start.replace(minute=0, second=0, microsecond=0)
    last = end.replace(minute=0, second=0, microsecond=0)
    if last < end:
        last += timedelta(hours=1)
    query = urlencode({
        "latitude": location["latitude"], "longitude": location["longitude"],
        "elevation": location["elevation_m"], "hourly": ",".join(VARIABLES),
        "timezone": "GMT", "wind_speed_unit": "ms", "temperature_unit": "celsius",
        "start_date": first.date().isoformat(), "end_date": last.date().isoformat(),
    })
    url = "https://api.open-meteo.com/v1/forecast?" + query
    raw = request_json(url)
    if raw.get("utc_offset_seconds") != 0:
        raise ValueError("Forecast did not use UTC")
    hourly = raw["hourly"]
    times = hourly["time"]
    rows = {}
    for i, stamp in enumerate(times):
        instant = datetime.fromisoformat(stamp).replace(tzinfo=timezone.utc)
        if first <= instant <= last:
            row = {"utc": instant.isoformat()}
            for variable, unit in VARIABLES.items():
                if raw["hourly_units"].get(variable) != unit:
                    raise ValueError(f"Unexpected forecast units for {variable}")
                value = hourly[variable][i]
                if not isinstance(value, (int, float)) or not math.isfinite(value):
                    raise ValueError(f"Missing forecast value for {variable} at {stamp}")
                if unit == "%" and not 0 <= value <= 100:
                    raise ValueError(f"Invalid forecast percentage for {variable}")
                if unit == "m/s" and value < 0:
                    raise ValueError(f"Invalid forecast wind for {variable}")
                row[variable] = value
            rows[instant] = row
    expected = [first + timedelta(hours=i)
                for i in range(int((last - first).total_seconds() // 3600) + 1)]
    if set(rows) != set(expected):
        raise ValueError("Forecast does not cover the complete observing window")
    return {
        "status": "available", "provider": "Open-Meteo", "source_url": url,
        "attribution": "Weather data by Open-Meteo (https://open-meteo.com/), CC BY 4.0",
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "model_issue_time_utc": None,
        "requested_location": location,
        "forecast_grid_location": {key: raw.get(key) for key in ("latitude", "longitude", "elevation")},
        "units": VARIABLES, "hourly": [rows[t] for t in expected],
        "summary": {
            "temperature_min_c": min(r["temperature_2m"] for r in rows.values()),
            "temperature_max_c": max(r["temperature_2m"] for r in rows.values()),
            "cloud_cover_max_percent": max(r["cloud_cover"] for r in rows.values()),
            "precipitation_probability_max_percent": max(r["precipitation_probability"] for r in rows.values()),
            "humidity_max_percent": max(r["relative_humidity_2m"] for r in rows.values()),
            "wind_max_ms": max(r["wind_speed_10m"] for r in rows.values()),
            "gust_max_ms": max(r["wind_gusts_10m"] for r in rows.values()),
        },
    }


def forecast_summary(weather):
    if weather.get("status") != "available":
        return "Weather unavailable: " + weather.get("error", "no forecast retrieved")
    s = weather["summary"]
    return (f"Temperature {s['temperature_min_c']:.1f}–{s['temperature_max_c']:.1f} °C; "
            f"max cloud {s['cloud_cover_max_percent']:.0f}%; "
            f"max rain probability {s['precipitation_probability_max_percent']:.0f}%; "
            f"max humidity {s['humidity_max_percent']:.0f}%; "
            f"max wind/gust {s['wind_max_ms']:.1f}/{s['gust_max_ms']:.1f} m/s.")
