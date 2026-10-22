"""Per-band sky coverage and cadence, from the camera's g/i/z quadrant layout."""

import json
from pathlib import Path

import healpy as hp
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

from ls4.config import BANDS, FOV_length, FOV_width, QUADRANT_BANDS

HALF_RA = FOV_width.to_value(u.deg) / 2  # The field is FOV_width wide in RA (see `region`).
HALF_DEC = FOV_length.to_value(u.deg) / 2


def band_visits(plans, nside=128):
    """One row per (HEALPix pixel, band, exposure) whose quadrant covers the pixel center.

    `plans` holds science pointings with `night`, `start_time_mjd`, `ra`, and `dec`
    (as returned by `ls4.metrics.load_plans`).
    """
    centers = SkyCoord(plans["ra"].to_numpy() * u.deg, plans["dec"].to_numpy() * u.deg)
    # Unit vectors at the 3x3 grid of quadrant corners (offsets of -1, 0, +1 half-fields).
    corner = {}
    for sx in (-1, 0, 1):
        for sy in (-1, 0, 1):
            point = centers.spherical_offsets_by(sx * HALF_RA * u.deg, sy * HALF_DEC * u.deg)
            corner[sx, sy] = point.cartesian.xyz.value.T

    pixels, bands, rows = [], [], []
    for (sx, sy), band in QUADRANT_BANDS.items():
        band_index = BANDS.index(band)
        for row in range(len(plans)):
            polygon = np.array([corner[0, 0][row], corner[sx, 0][row],
                                corner[sx, sy][row], corner[0, sy][row]])
            covered = hp.query_polygon(nside, polygon)
            pixels.append(covered)
            bands.append(np.full(len(covered), band_index))
            rows.append(np.full(len(covered), row))
    rows = np.concatenate(rows)
    return pd.DataFrame({
        "pixel": np.concatenate(pixels),
        "band": np.array(BANDS)[np.concatenate(bands)],
        "night": plans["night"].to_numpy()[rows],
        "mjd": plans["start_time_mjd"].to_numpy()[rows],
    })


def pixel_stats(visits, nside):
    """Per (pixel, band): nights observed and gaps between observed nights (days)."""
    nightly = visits.groupby(["pixel", "band", "night"])["mjd"].min().reset_index()
    nightly = nightly.sort_values(["pixel", "band", "mjd"])
    nightly["gap_days"] = nightly.groupby(["pixel", "band"])["mjd"].diff()
    stats = nightly.groupby(["pixel", "band"]).agg(
        nights=("night", "size"), median_gap_days=("gap_days", "median"),
        max_gap_days=("gap_days", "max"), last_mjd=("mjd", "max")).reset_index()
    theta, phi = hp.pix2ang(nside, stats["pixel"].to_numpy())
    stats.insert(1, "ra", np.degrees(phi))
    stats.insert(2, "dec", 90 - np.degrees(theta))
    return stats, nightly


def nightly_area(nightly, night_order, nside):
    """Area observed each night per band, and cumulative area per band and in all three."""
    pixel_area = hp.nside2pixarea(nside, degrees=True)
    rows = []
    seen = {band: set() for band in BANDS}
    for night in night_order:
        tonight = nightly[nightly["night"] == night]
        row = {"night": night}
        for band in BANDS:
            pixels = set(tonight.loc[tonight["band"] == band, "pixel"])
            seen[band] |= pixels
            row[f"{band}_area_deg2"] = len(pixels) * pixel_area
            row[f"{band}_cumulative_deg2"] = len(seen[band]) * pixel_area
        row["all_bands_cumulative_deg2"] = len(set.intersection(*seen.values())) * pixel_area
        rows.append(row)
    return pd.DataFrame(rows)


def write_band_report(plans, output_dir, nside=128, show=False):
    """Per-band tables, an interactive sky map, cadence CDFs, and area growth."""
    from ls4 import plots

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    night_order = sorted(plans["night"].unique())
    legacy = sorted(plans.loc[~plans["configuration"].astype(str).str.contains("phase"), "night"].unique())
    if legacy:
        # Plans from before the alternating dither always used phase 0, so each half-row
        # of a field only ever gets g or z.
        print(f"WARNING: {len(legacy)} of {len(night_order)} plans predate the alternating dither "
              f"({legacy[0]} to {legacy[-1]}); their fields get g or z, never both.")
    print(f"Projecting {len(plans)} exposures onto g/i/z quadrants (nside={nside})...")
    stats, nightly = pixel_stats(band_visits(plans, nside), nside)
    area = nightly_area(nightly, night_order, nside)

    stats.to_csv(output_dir / "band_pixels.csv", index=False)
    area.to_csv(output_dir / "band_nightly_area.csv", index=False)
    print(f"  wrote {output_dir / 'band_pixels.csv'} and {output_dir / 'band_nightly_area.csv'}")
    write_band_sky_map(nightly, night_order, plans, nside, output_dir / "band_sky_map.html")
    plots.band_cadence_cdf(nightly, output_dir, show)
    plots.band_area_growth(area, output_dir, show)

    last = area.iloc[-1]
    lines = [f"{band}: {last[f'{band}_cumulative_deg2']:,.0f} deg^2 observed, median gap "
             f"{nightly.loc[nightly['band'] == band, 'gap_days'].median():.2f} d"
             for band in BANDS]
    lines.append(f"g, i and z: {last['all_bands_cumulative_deg2']:,.0f} deg^2")
    if legacy:
        lines.append(f"{len(legacy)} of {len(night_order)} plans predate the alternating dither "
                     "(fixed RA dither: g or z per half-row, never both)")
    (output_dir / "band_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return stats, area


def write_band_sky_map(nightly, night_order, plans, nside, output):
    """Self-contained HTML map with band, metric, and through-night controls."""
    from ls4.sky_maps import galactic_lines

    index = {night: i for i, night in enumerate(night_order)}
    night_mjd = plans.groupby("night")["start_time_mjd"].min().reindex(night_order)
    observed = {}
    for (pixel, band), group in nightly.groupby(["pixel", "band"]):
        observed.setdefault(pixel, {b: [] for b in BANDS})[band] = sorted(
            index[night] for night in group["night"])
    pixels = np.array(sorted(observed))
    theta, phi = hp.pix2ang(nside, pixels)
    data = {
        "bands": list(BANDS),
        "nights": [{"label": night, "mjd": round(float(mjd), 4)} for night, mjd in night_mjd.items()],
        "pixelDeg": float(np.degrees(hp.nside2resol(nside))),
        "pixels": [[round(float(np.degrees(p)), 3), round(float(90 - np.degrees(t)), 3),
                    *[observed[pixel][band] for band in BANDS]]
                   for pixel, t, p in zip(pixels, theta, phi)],
        "galactic": galactic_lines(),
    }
    output = Path(output)
    output.write_text(BAND_MAP_HTML.replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    print(f"  wrote {output}")


BAND_MAP_HTML = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LS4 sky coverage and cadence by band</title>
<style>
  .viz-root {
    color-scheme: light;
    --page: #f4f4f2; --surface-1: #fcfcfb; --border: #d9d8d4;
    --text-primary: #0b0b0b; --text-secondary: #52514e; --empty: #e7e6e2;
    --plane: #8a3ffc; --tooltip-bg: rgba(11,11,11,.92);
    --g: #1baf7a; --i: #2a78d6; --z: #eb6834;
    --ramp: #cde2fb, #86b6ef, #3987e5, #256abf, #104281;
  }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) .viz-root {
      color-scheme: dark;
      --page: #111110; --surface-1: #1a1a19; --border: #383835;
      --text-primary: #ffffff; --text-secondary: #c3c2b7; --empty: #2a2a28;
      --plane: #b39dff; --tooltip-bg: rgba(40,40,38,.96);
      --g: #199e70; --i: #3987e5; --z: #d95926;
      --ramp: #0d366b, #184f95, #2a78d6, #6da7ec, #cde2fb;
    }
  }
  body { margin: 0; }
  .viz-root { background: var(--page); color: var(--text-primary); min-height: 100vh;
              font-family: system-ui, -apple-system, sans-serif; }
  main { max-width: 1500px; margin: 0 auto; padding: 20px 24px 28px; }
  h1 { font-size: 1.3rem; margin: 0 0 4px; }
  p { color: var(--text-secondary); margin: 4px 0 12px; line-height: 1.4; font-size: .92rem; }
  .controls { display: flex; flex-wrap: wrap; align-items: center; gap: 10px 22px; margin: 10px 0 14px; }
  .group { display: flex; align-items: center; gap: 8px; }
  .group > span { font-weight: 600; font-size: .9rem; }
  .seg { display: inline-flex; border: 1px solid var(--border); border-radius: 6px; overflow: hidden; }
  .seg button { border: 0; background: var(--surface-1); color: var(--text-primary); padding: 5px 11px;
                font: inherit; font-size: .88rem; cursor: pointer; display: inline-flex; align-items: center; gap: 6px; }
  .seg button + button { border-left: 1px solid var(--border); }
  .seg button[aria-pressed="true"] { background: var(--text-primary); color: var(--surface-1); }
  .swatch { width: 10px; height: 10px; border-radius: 2px; display: inline-block; }
  .night { flex: 1; min-width: 320px; }
  .night input { flex: 1; accent-color: var(--i); }
  output { font-variant-numeric: tabular-nums; min-width: 9rem; font-size: .88rem; }
  .chart { position: relative; background: var(--surface-1); border: 1px solid var(--border);
           border-radius: 8px; overflow: hidden; }
  canvas { display: block; width: 100%; height: auto; aspect-ratio: 2 / 1; }
  .tooltip { display: none; position: absolute; z-index: 2; pointer-events: none; background: var(--tooltip-bg);
             color: #fff; border-radius: 6px; padding: 8px 10px; font-size: .82rem; line-height: 1.5; white-space: nowrap; }
  .tooltip .row { display: flex; align-items: center; gap: 6px; }
  .legend { display: flex; flex-wrap: wrap; align-items: center; gap: 10px; margin-top: 12px;
            font-size: .85rem; color: var(--text-secondary); }
  .ramp { width: 200px; height: 12px; border-radius: 2px; }
  .key { width: 14px; height: 14px; border-radius: 2px; background: var(--empty); }
  .plane { display: inline-block; width: 22px; border-top: 2px solid var(--plane); }
</style>
</head>
<body>
<div class="viz-root">
<main>
  <h1>LS4 sky coverage and cadence by band</h1>
  <p>Each cell is a HEALPix pixel colored by the selected band and metric, using the camera's
     g/i/z quadrants for every exposure through the selected night. Hover a cell for all three bands.
     With g + i + z, "nights observed" counts bands seen and the cadence metrics use the slowest band.</p>
  <div class="controls">
    <div class="group"><span>Band</span><div class="seg" id="band"></div></div>
    <div class="group"><span>Metric</span><div class="seg" id="metric"></div></div>
    <div class="group night"><span>Through</span>
      <input id="night" type="range" min="0" step="1"><output id="night-label"></output></div>
  </div>
  <div class="chart"><canvas id="sky" width="1400" height="700" role="img"
    aria-label="Sky map of LS4 coverage by band"></canvas><div class="tooltip" id="tooltip"></div></div>
  <div class="legend"><span id="minimum"></span><span class="ramp" id="ramp"></span><span id="maximum"></span>
    <span id="unit"></span><span class="key"></span><span id="empty-label"></span>
    <span class="plane"></span><span>Galactic plane</span></div>
  <p id="note"></p>
</main>
</div>
<script>
const DATA = __DATA__;
const BAND_NAMES = {g: 'g', i: 'i', z: 'z', all: 'g + i + z'};
const METRICS = {
  nights: {label: 'Nights observed', unit: 'nights', empty: 'only in other bands'},
  gap: {label: 'Median gap', unit: 'days between observed nights', empty: 'fewer than two nights'},
  since: {label: 'Days since last', unit: 'days since last observation', empty: 'only in other bands'},
};
// Optional URL parameters, e.g. band_sky_map.html?band=all&metric=gap&night=10
const params = new URLSearchParams(location.search);
const state = {band: params.get('band') || 'i', metric: params.get('metric') || 'nights',
               night: Math.min(DATA.nights.length, Number(params.get('night')) || DATA.nights.length) - 1};
const canvas = document.getElementById('sky'), ctx = canvas.getContext('2d');
const tooltip = document.getElementById('tooltip');
const css = name => getComputedStyle(document.querySelector('.viz-root')).getPropertyValue(name).trim();
let plot, width, height, values, colorMax, ramp;

function button(group, key, text, swatch) {
  const b = document.createElement('button');
  b.type = 'button'; b.dataset.key = key;
  if (swatch) { const s = document.createElement('span'); s.className = 'swatch';
    s.style.background = `var(--${swatch})`; b.append(s); }
  b.append(text);
  b.addEventListener('click', () => { state[group] = key; sync(); });
  document.getElementById(group).append(b);
}
DATA.bands.forEach(band => button('band', band, band, band));
button('band', 'all', 'g + i + z');
Object.entries(METRICS).forEach(([key, m]) => button('metric', key, m.label));
const slider = document.getElementById('night');
slider.max = DATA.nights.length - 1; slider.value = state.night;
slider.addEventListener('input', () => { state.night = Number(slider.value); sync(); });

function bandNights(pixel, band) {
  return pixel[2 + DATA.bands.indexOf(band)].filter(n => n <= state.night);
}
function metricFor(pixel, band) {
  const nights = bandNights(pixel, band);
  if (state.metric === 'nights') return nights.length || null;
  if (state.metric === 'since') return nights.length
    ? DATA.nights[state.night].mjd - DATA.nights[nights[nights.length - 1]].mjd : null;
  if (nights.length < 2) return null;
  const gaps = nights.slice(1).map((n, k) => DATA.nights[n].mjd - DATA.nights[nights[k]].mjd).sort((a, b) => a - b);
  const mid = (gaps.length - 1) / 2;
  return (gaps[Math.floor(mid)] + gaps[Math.ceil(mid)]) / 2;
}
function value(pixel) {
  if (state.band !== 'all') return metricFor(pixel, state.band);
  if (state.metric === 'nights') {
    // Nights through which every band has been seen at least once: count bands seen.
    const seen = DATA.bands.filter(b => bandNights(pixel, b).length).length;
    return seen || null;
  }
  // For all bands, the slowest band sets the cadence.
  const per = DATA.bands.map(b => metricFor(pixel, b));
  return per.some(v => v === null) ? null : Math.max(...per);
}
function shade(v) {
  const t = Math.max(0, Math.min(1, v / colorMax)) * (ramp.length - 1);
  const a = ramp[Math.floor(t)], b = ramp[Math.ceil(t)], w = t - Math.floor(t);
  const ch = i => Math.round(parseInt(a.slice(i, i + 2), 16) * (1 - w) + parseInt(b.slice(i, i + 2), 16) * w);
  return `rgb(${ch(1)},${ch(3)},${ch(5)})`;
}
function sync() {
  document.querySelectorAll('.seg button').forEach(b =>
    b.setAttribute('aria-pressed', String(state[b.parentElement.id] === b.dataset.key)));
  document.getElementById('night-label').value =
    `${DATA.nights[state.night].label} (night ${state.night + 1} of ${DATA.nights.length})`;
  ramp = css('--ramp').split(',').map(s => s.trim());
  values = DATA.pixels.map(value);
  const shown = values.filter(v => v !== null).sort((a, b) => a - b);
  const allNights = state.band === 'all' && state.metric === 'nights';
  colorMax = allNights ? 3 : state.metric === 'nights' ? state.night + 1
           : (shown.length ? shown[Math.floor(0.99 * (shown.length - 1))] : 1) || 1;
  const m = METRICS[state.metric];
  document.getElementById('minimum').textContent = allNights ? '1' : '0';
  document.getElementById('maximum').textContent = state.metric === 'nights' ? String(colorMax) : colorMax.toFixed(1);
  document.getElementById('unit').textContent = allNights ? 'of g, i, z seen' :
    state.band === 'all' ? `${m.unit} (slowest band)` : m.unit;
  document.getElementById('empty-label').textContent = m.empty;
  document.getElementById('ramp').style.background = `linear-gradient(to right, ${ramp.join(',')})`;
  document.getElementById('note').textContent = state.metric === 'gap'
    ? 'Gaps are measured between nights on which the pixel was observed in that band; the scale tops out at the 99th percentile.'
    : 'Pixels never observed in the selected band are left blank; gray pixels were observed only in other bands.';
  if (plot) draw();
}

const x = ra => plot.left + ((ra + 180) % 360 + 360) % 360 / 360 * plot.width;
const xLon = lon => plot.left + (lon + 180) / 360 * plot.width;
const y = dec => plot.top + (90 - dec) / 180 * plot.height;
function drawCurve(points, lineWidth, dash) {
  ctx.beginPath(); let open = false;
  for (const p of points) {
    if (!p) { open = false; continue; }
    if (!open) { ctx.moveTo(xLon(p[0]), y(p[1])); open = true; } else ctx.lineTo(xLon(p[0]), y(p[1]));
  }
  ctx.strokeStyle = css('--plane'); ctx.lineWidth = lineWidth; ctx.setLineDash(dash); ctx.stroke(); ctx.setLineDash([]);
}
function draw() {
  ctx.fillStyle = css('--surface-1'); ctx.fillRect(0, 0, width, height);
  const empty = css('--empty'), size = DATA.pixelDeg;
  DATA.pixels.forEach((pixel, k) => {
    const v = values[k];
    const any = DATA.bands.some(b => bandNights(pixel, b).length);
    if (v === null && !any) return;
    ctx.fillStyle = v === null ? empty : shade(v);
    const w = Math.min(360, size / Math.max(0.05, Math.cos(pixel[1] * Math.PI / 180))) / 360 * plot.width;
    const h = size / 180 * plot.height;
    ctx.fillRect(x(pixel[0]) - w / 2, y(pixel[1]) - h / 2, w + 0.6, h + 0.6);
  });
  drawCurve(DATA.galactic['0'], 1.5, []);
  ctx.fillStyle = css('--text-secondary'); ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'center';
  for (let lon = -180; lon <= 180; lon += 30) ctx.fillText(`${lon}°`, xLon(lon), plot.top + plot.height + 18);
  ctx.textAlign = 'right';
  for (let dec = -90; dec <= 90; dec += 30) ctx.fillText(`${dec}°`, plot.left - 7, y(dec) + 4);
  ctx.textAlign = 'center';
  ctx.fillText('Right ascension (wrapped degrees)', plot.left + plot.width / 2, height - 3);
}
function resize() {
  const rect = canvas.getBoundingClientRect(), ratio = window.devicePixelRatio || 1;
  width = rect.width; height = rect.height;
  canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  plot = {left: 48, top: 10, width: width - 60, height: height - 46};
  draw();
}
// Hover: bucket pixels by whole degree for a quick nearest-pixel lookup.
const buckets = new Map();
DATA.pixels.forEach((p, k) => {
  const key = `${Math.floor(p[0])},${Math.floor(p[1])}`;
  if (!buckets.has(key)) buckets.set(key, []); buckets.get(key).push(k);
});
canvas.addEventListener('mousemove', event => {
  const rect = canvas.getBoundingClientRect(), mx = event.clientX - rect.left, my = event.clientY - rect.top;
  const ra = (((mx - plot.left) / plot.width * 360 - 180) % 360 + 360) % 360, dec = 90 - (my - plot.top) / plot.height * 180;
  let best = -1, bestD = DATA.pixelDeg * 0.8;
  for (let da = -2; da <= 2; da++) for (let dd = -1; dd <= 1; dd++) {
    for (const k of buckets.get(`${(Math.floor(ra) + da + 360) % 360},${Math.floor(dec) + dd}`) || []) {
      const p = DATA.pixels[k];
      const d = Math.hypot(((p[0] - ra + 540) % 360 - 180) * Math.cos(dec * Math.PI / 180), p[1] - dec);
      if (d < bestD) { bestD = d; best = k; }
    }
  }
  if (best < 0) { tooltip.style.display = 'none'; return; }
  const p = DATA.pixels[best], saved = state.band;
  const fmt = (b, metric) => { state.band = b; const keep = state.metric; state.metric = metric;
    const v = metricFor(p, b); state.metric = keep; return v; };
  const rows = DATA.bands.map(b => {
    const n = bandNights(p, b).length, gap = fmt(b, 'gap');
    return `<div class="row"><span class="swatch" style="background:var(--${b})"></span>${b}: ` +
      `${n} night${n === 1 ? '' : 's'}${gap === null ? '' : `, median gap ${gap.toFixed(2)} d`}</div>`;
  }).join('');
  state.band = saved;
  tooltip.innerHTML = `<strong>RA ${p[0].toFixed(1)}°, Dec ${p[1].toFixed(1)}°</strong>${rows}`;
  tooltip.style.left = `${Math.min(mx + 14, width - 250)}px`;
  tooltip.style.top = `${Math.max(5, my - 70)}px`;
  tooltip.style.display = 'block';
});
canvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', sync);
new ResizeObserver(resize).observe(canvas);
sync();
</script>
</body>
</html>'''
