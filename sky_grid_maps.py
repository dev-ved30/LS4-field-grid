"""Fast, self-contained 2D sky-grid maps for monthly survey statistics."""

import json
from pathlib import Path

import numpy as np

from make_observing_plan_gif import galactic_to_equatorial


def wrap_lon(ra):
    return (float(ra) + 180) % 360 - 180


def grid_cells(per_field, gap_distributions):
    """Represent each field by its grid-ring cell in flat RA/Dec coordinates."""
    declinations = np.sort(per_field['dec'].unique())
    lower = np.r_[-90.0, (declinations[:-1] + declinations[1:]) / 2]
    upper = np.r_[(declinations[:-1] + declinations[1:]) / 2, 90.0]
    cells = []
    for dec, bottom, top in zip(declinations, lower, upper):
        ring = per_field[per_field['dec'] == dec]
        half_width = 180.0 / len(ring)
        for row in ring.itertuples():
            lon = wrap_lon(row.ra)
            left, right = lon - half_width, lon + half_width
            if left < -180 - 1e-8:
                parts = [[left + 360, 180], [-180, right]]
            elif right > 180 + 1e-8:
                parts = [[left, 180], [-180, right - 360]]
            else:
                parts = [[max(-180, left), min(180, right)]]
            parts = [part for part in parts if part[1] - part[0] > 1e-8]
            cells.append({
                'name': row.base_field, 'ra': round(float(row.ra), 3),
                'dec': round(float(row.dec), 3), 'bottom': float(bottom),
                'top': float(top), 'parts': parts, 'visits': int(row.visits),
                'nights': int(row.observed_nights),
                'gaps': sorted(round(float(gap), 6) for gap in
                               gap_distributions.get(row.base_field, [])),
            })
    return cells


def galactic_lines():
    lines = {}
    for latitude in (-10, 0, 10):
        points = []
        previous = None
        for longitude in range(0, 361, 2):
            ra, dec = galactic_to_equatorial(longitude % 360, latitude)
            lon = wrap_lon(ra)
            if previous is not None and abs(lon - previous) > 180:
                points.append(None)
            points.append([round(lon, 4), round(dec, 4)])
            previous = lon
        lines[str(latitude)] = points
    return lines


HTML = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root { color-scheme: light; font-family: system-ui, -apple-system, sans-serif; }
  body { margin: 0; background: #f1f5f9; color: #172554; }
  main { max-width: 1500px; margin: 0 auto; padding: 20px 24px 28px; }
  h1 { font-size: 1.35rem; margin: 0 0 4px; }
  p { color: #475569; margin: 4px 0 14px; line-height: 1.4; }
  .controls { display: flex; align-items: center; gap: 16px; margin: 10px 0 16px; }
  .controls[hidden] { display: none; }
  .controls label { white-space: nowrap; font-weight: 600; }
  input[type=range] { flex: 1; accent-color: #7e22ce; cursor: pointer; }
  output { min-width: 3.2rem; text-align: right; font-variant-numeric: tabular-nums; }
  .chart { position: relative; background: white; border: 1px solid #cbd5e1;
           border-radius: 8px; overflow: hidden; }
  canvas { display: block; width: 100%; height: auto; aspect-ratio: 2 / 1; }
  .tooltip { display: none; position: absolute; z-index: 2; pointer-events: none;
             background: rgba(15,23,42,.94); color: white; border-radius: 5px;
             padding: 8px 10px; font-size: .84rem; line-height: 1.45; white-space: nowrap; }
  .legend { display: flex; align-items: center; gap: 10px; margin-top: 12px;
            font-size: .85rem; color: #334155; }
  .ramp { width: 190px; height: 12px; border: 3px solid transparent; background: white; }
  .key { width: 15px; height: 15px; border: 2px solid #b7c3d5; background: white; }
  .plane { display: inline-block; width: 22px; border-top: 3px solid #d946ef; }
  .note { font-size: .82rem; margin-top: 12px; }
</style>
</head>
<body>
<main>
  <h1>__TITLE__</h1>
  <p>Flat equatorial RA–Dec field grid. Hover over a cell for its field statistics.</p>
  <div class="controls" id="controls">
    <label for="percentile">Gap percentile</label>
    <input id="percentile" type="range" min="0" max="100" step="1" value="50">
    <output id="percentile-value" for="percentile">P50</output>
  </div>
  <div class="chart"><canvas id="sky" width="1400" height="700"
    role="img" aria-label="__TITLE__"></canvas><div class="tooltip" id="tooltip"></div></div>
  <div class="legend"><span id="minimum">0</span><span class="ramp" id="ramp"></span>
    <span id="maximum"></span><span id="unit"></span>
    <span class="key"></span><span id="empty-label"></span>
    <span class="plane"></span><span>Galactic plane (dotted: ±10°)</span></div>
  <p class="note" id="note"></p>
</main>
<script>
const DATA = __DATA__;
const MODE = '__MODE__';
const canvas = document.getElementById('sky');
const chart = canvas.parentElement;
const ctx = canvas.getContext('2d');
const slider = document.getElementById('percentile');
const tooltip = document.getElementById('tooltip');
const palette = MODE === 'gaps'
  ? ['#0d0887', '#7e03a8', '#cc4678', '#f89540', '#f0f921']
  : ['#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#fde725'];
document.getElementById('controls').hidden = MODE !== 'gaps';
document.getElementById('unit').textContent = MODE === 'gaps' ? 'days' : 'visits';
document.getElementById('empty-label').textContent = MODE === 'gaps'
  ? 'fewer than two initial visits' : 'no visits';
document.getElementById('note').textContent = MODE === 'gaps'
  ? 'Only field borders are colored. The percentile uses each field’s gaps between consecutive initial visits; dithers are excluded. The scale spans 0 to the 99th percentile of displayed values.'
  : 'Only field borders are colored. Cells follow the catalog’s declination rings and RA spacing; they approximate tiles, not exact camera footprints.';
document.getElementById('ramp').style.borderImage = `linear-gradient(to right, ${palette.join(',')}) 1`;

let width, height, plot, values, colorMax;
const x = lon => plot.left + (lon + 180) / 360 * plot.width;
const y = dec => plot.top + (90 - dec) / 180 * plot.height;

function percentile(sorted, p) {
  if (!sorted.length) return null;
  const index = (sorted.length - 1) * p / 100;
  const lo = Math.floor(index), hi = Math.ceil(index);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (index - lo);
}
function shade(value) {
  const t = Math.max(0, Math.min(1, value / colorMax)) * (palette.length - 1);
  const a = palette[Math.floor(t)], b = palette[Math.min(palette.length - 1, Math.ceil(t))];
  const weight = t - Math.floor(t);
  const channel = i => Math.round(parseInt(a.slice(i, i + 2), 16) * (1 - weight) +
                                  parseInt(b.slice(i, i + 2), 16) * weight);
  return `rgb(${channel(1)},${channel(3)},${channel(5)})`;
}
function updateValues() {
  const p = Number(slider.value);
  document.getElementById('percentile-value').value = `P${p}`;
  values = DATA.cells.map(cell => MODE === 'gaps' ? percentile(cell.gaps, p) :
                           (cell.visits ? cell.visits : null));
  const nonempty = values.filter(value => value !== null).sort((a, b) => a - b);
  colorMax = MODE === 'gaps' ? (percentile(nonempty, 99) || 1) : DATA.maxVisits;
  document.getElementById('minimum').textContent = '0';
  document.getElementById('maximum').textContent = MODE === 'gaps'
    ? colorMax.toFixed(1) : String(colorMax);
}
function drawCurve(points, color, lineWidth, dash) {
  ctx.beginPath();
  let open = false;
  for (const point of points) {
    if (!point) { open = false; continue; }
    if (!open) { ctx.moveTo(x(point[0]), y(point[1])); open = true; }
    else ctx.lineTo(x(point[0]), y(point[1]));
  }
  ctx.strokeStyle = 'rgba(255,255,255,.92)'; ctx.lineWidth = lineWidth + 2;
  ctx.setLineDash(dash); ctx.stroke();
  ctx.strokeStyle = color; ctx.lineWidth = lineWidth; ctx.stroke();
  ctx.setLineDash([]);
}
function draw() {
  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, width, height);
  for (let i = 0; i < DATA.cells.length; i++) {
    const cell = DATA.cells[i];
    ctx.strokeStyle = values[i] === null ? '#b7c3d5' : shade(values[i]);
    ctx.lineWidth = 1.8;
    for (const part of cell.parts) {
      const left = x(part[0]), right = x(part[1]);
      const top = y(cell.top), bottom = y(cell.bottom);
      ctx.strokeRect(left + 1, top + 1, Math.max(0, right - left - 2),
                     Math.max(0, bottom - top - 2));
    }
  }
  drawCurve(DATA.galactic['-10'], '#d946ef', 1.3, [5, 4]);
  drawCurve(DATA.galactic['10'], '#d946ef', 1.3, [5, 4]);
  drawCurve(DATA.galactic['0'], '#a21caf', 2.4, []);
  ctx.fillStyle = '#334155'; ctx.font = '12px system-ui, sans-serif';
  ctx.textAlign = 'center';
  for (let lon = -180; lon <= 180; lon += 30)
    ctx.fillText(`${lon}°`, x(lon), plot.top + plot.height + 19);
  ctx.textAlign = 'right';
  for (let dec = -90; dec <= 90; dec += 30)
    ctx.fillText(`${dec}°`, plot.left - 7, y(dec) + 4);
  ctx.textAlign = 'center';
  ctx.fillText('Right ascension (wrapped degrees)', plot.left + plot.width / 2, height - 3);
}
function resize() {
  const rect = canvas.getBoundingClientRect();
  width = rect.width; height = rect.height;
  const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  plot = {left: 53, top: 12, width: width - 68, height: height - 51};
  draw();
}
slider.addEventListener('input', () => { updateValues(); draw(); });
canvas.addEventListener('mousemove', event => {
  const rect = canvas.getBoundingClientRect();
  const mx = event.clientX - rect.left, my = event.clientY - rect.top;
  if (mx < plot.left || mx > plot.left + plot.width ||
      my < plot.top || my > plot.top + plot.height) { tooltip.style.display = 'none'; return; }
  const lon = (mx - plot.left) / plot.width * 360 - 180;
  const dec = 90 - (my - plot.top) / plot.height * 180;
  const index = DATA.cells.findIndex(cell => dec >= cell.bottom && dec <= cell.top &&
                                  cell.parts.some(part => lon >= part[0] && lon <= part[1]));
  if (index < 0) { tooltip.style.display = 'none'; return; }
  const cell = DATA.cells[index], value = values[index];
  tooltip.innerHTML = `<strong>${cell.name}</strong><br>RA ${cell.ra}° · Dec ${cell.dec}°<br>` +
    (MODE === 'gaps' ? `${cell.gaps.length} inter-night gaps<br>` +
      (value === null ? 'No measured gap' : `P${slider.value} gap: ${value.toFixed(2)} days`)
      : `${cell.visits} visits across ${cell.nights} nights`);
  tooltip.style.left = `${Math.min(mx + 12, width - 205)}px`;
  tooltip.style.top = `${Math.max(5, my - 60)}px`;
  tooltip.style.display = 'block';
});
canvas.addEventListener('mouseleave', () => { tooltip.style.display = 'none'; });
updateValues();
new ResizeObserver(resize).observe(canvas);
</script>
</body>
</html>'''


def write_sky_grid_maps(per_field, gap_distributions, output_dir):
    """Write visit-count and percentile-selectable cadence maps without Plotly."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data = {'cells': grid_cells(per_field, gap_distributions),
            'galactic': galactic_lines(),
            'maxVisits': int(per_field['visits'].max())}
    payload = json.dumps(data, separators=(',', ':'))
    for mode, title, name in (
        ('visits', 'LS4 field grid · total visits', 'visits_per_field_map.html'),
        ('gaps', 'LS4 field grid · gap percentile', 'visit_gap_percentile_map.html'),
    ):
        html = HTML.replace('__TITLE__', title).replace('__DATA__', payload).replace('__MODE__', mode)
        path = output_dir / name
        path.write_text(html)
        print(f'  wrote {path}')


COMPARISON_HTML = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root { color-scheme: light; font-family: system-ui, -apple-system, sans-serif; }
  body { margin: 0; background: #f1f5f9; color: #172554; }
  main { max-width: 1900px; margin: auto; padding: 20px 22px 30px; }
  h1 { font-size: 1.4rem; margin: 0 0 6px; }
  p { color: #475569; margin: 4px 0 15px; }
  .controls { display: flex; gap: 14px; align-items: center; margin-bottom: 17px; }
  .controls[hidden] { display: none; }
  .controls label { white-space: nowrap; font-weight: 600; }
  input[type=range] { flex: 1; accent-color: #7e22ce; cursor: pointer; }
  output { min-width: 3rem; font-variant-numeric: tabular-nums; }
  .maps { display: grid; grid-template-columns: repeat(2, minmax(0,1fr)); gap: 16px; }
  section { min-width: 0; position: relative; background: white;
            border: 1px solid #cbd5e1; border-radius: 8px; padding: 10px; }
  h2 { font-size: 1rem; margin: 0 0 5px; }
  canvas { display: block; width: 100%; height: auto; aspect-ratio: 2 / 1; }
  .tooltip { display: none; position: absolute; z-index: 2; pointer-events: none;
             background: rgba(15,23,42,.94); color: white; border-radius: 5px;
             padding: 7px 9px; font-size: .8rem; line-height: 1.4; white-space: nowrap; }
  .legend { display: flex; align-items: center; gap: 10px; margin-top: 13px;
            font-size: .85rem; color: #334155; }
  .ramp { width: 190px; height: 12px; border: 3px solid transparent; background: white; }
  .key { width: 15px; height: 15px; border: 2px solid #b7c3d5; background: white; }
  .plane { display: inline-block; width: 22px; border-top: 3px solid #d946ef; }
  .note { font-size: .82rem; margin-top: 13px; }
  @media (max-width: 1000px) { .maps { grid-template-columns: 1fr; } }
</style>
</head>
<body>
<main>
  <h1>__TITLE__</h1>
  <p>Same 30 observing nights and field grid. Both panels use the same border-color scale.</p>
  <div class="controls" id="controls">
    <label for="percentile">Gap percentile</label>
    <input id="percentile" type="range" min="0" max="100" step="1" value="50">
    <output id="percentile-value" for="percentile">P50</output>
  </div>
  <div class="maps">
    <section><h2>Original plan</h2><canvas id="old" width="800" height="400"></canvas>
      <div class="tooltip" id="old-tooltip"></div></section>
    <section><h2>Clustered replay</h2><canvas id="new" width="800" height="400"></canvas>
      <div class="tooltip" id="new-tooltip"></div></section>
  </div>
  <div class="legend"><span>0</span><span class="ramp" id="ramp"></span>
    <span id="maximum"></span><span id="unit"></span>
    <span class="key"></span><span id="empty-label"></span>
    <span class="plane"></span><span>Galactic plane (dotted: ±10°)</span></div>
  <p class="note" id="note"></p>
</main>
<script>
const DATA = __DATA__;
const MODE = '__MODE__';
const slider = document.getElementById('percentile');
const palette = MODE === 'gaps'
  ? ['#0d0887', '#7e03a8', '#cc4678', '#f89540', '#f0f921']
  : ['#440154', '#414487', '#2a788e', '#22a884', '#7ad151', '#fde725'];
document.getElementById('controls').hidden = MODE !== 'gaps';
document.getElementById('unit').textContent = MODE === 'gaps' ? 'days' : 'visits';
document.getElementById('empty-label').textContent = MODE === 'gaps'
  ? 'fewer than two initial visits' : 'no visits';
document.getElementById('note').textContent = MODE === 'gaps'
  ? 'The shared slider selects each field’s percentile of gaps between consecutive initial visits. Dithers are excluded. The scale ends at the pooled 99th percentile of displayed values.'
  : 'Field interiors are white. The shared scale covers both plans; cells follow catalog grid spacing rather than exact camera footprints.';
document.getElementById('ramp').style.borderImage = `linear-gradient(to right, ${palette.join(',')}) 1`;

const canvases = [document.getElementById('old'), document.getElementById('new')];
const lists = [DATA.oldCells, DATA.newCells];
const geometry = [];
let values, colorMax;
function percentile(sorted, p) {
  if (!sorted.length) return null;
  const index = (sorted.length - 1) * p / 100;
  const lo = Math.floor(index), hi = Math.ceil(index);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (index - lo);
}
function shade(value) {
  const t = Math.max(0, Math.min(1, value / colorMax)) * (palette.length - 1);
  const a = palette[Math.floor(t)], b = palette[Math.min(palette.length - 1, Math.ceil(t))];
  const fraction = t - Math.floor(t);
  const channel = i => Math.round(parseInt(a.slice(i, i + 2), 16) * (1 - fraction) +
                                  parseInt(b.slice(i, i + 2), 16) * fraction);
  return `rgb(${channel(1)},${channel(3)},${channel(5)})`;
}
function updateValues() {
  const p = Number(slider.value);
  document.getElementById('percentile-value').value = `P${p}`;
  values = lists.map(cells => cells.map(cell => MODE === 'gaps'
    ? percentile(cell.gaps, p) : (cell.visits ? cell.visits : null)));
  const pooled = values.flat().filter(value => value !== null).sort((a, b) => a - b);
  colorMax = MODE === 'gaps' ? (percentile(pooled, 99) || 1) : DATA.maxVisits;
  document.getElementById('maximum').textContent = MODE === 'gaps'
    ? colorMax.toFixed(1) : String(colorMax);
}
function drawMap(index) {
  const canvas = canvases[index], ctx = canvas.getContext('2d');
  const rect = canvas.getBoundingClientRect(), width = rect.width, height = rect.height;
  const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  const plot = {left: 45, top: 9, width: width - 57, height: height - 40};
  geometry[index] = plot;
  const x = lon => plot.left + (lon + 180) / 360 * plot.width;
  const y = dec => plot.top + (90 - dec) / 180 * plot.height;
  ctx.fillStyle = '#ffffff'; ctx.fillRect(0, 0, width, height);
  for (let i = 0; i < lists[index].length; i++) {
    const cell = lists[index][i];
    ctx.strokeStyle = values[index][i] === null ? '#b7c3d5' : shade(values[index][i]);
    ctx.lineWidth = 1.5;
    for (const part of cell.parts) {
      const left = x(part[0]), right = x(part[1]);
      const top = y(cell.top), bottom = y(cell.bottom);
      ctx.strokeRect(left + .7, top + .7, Math.max(0, right - left - 1.4),
                     Math.max(0, bottom - top - 1.4));
    }
  }
  for (const latitude of ['-10', '10', '0']) {
    ctx.beginPath(); let open = false;
    for (const point of DATA.galactic[latitude]) {
      if (!point) { open = false; continue; }
      if (!open) { ctx.moveTo(x(point[0]), y(point[1])); open = true; }
      else ctx.lineTo(x(point[0]), y(point[1]));
    }
    ctx.setLineDash(latitude === '0' ? [] : [4, 3]);
    ctx.strokeStyle = 'rgba(255,255,255,.95)'; ctx.lineWidth = latitude === '0' ? 4 : 3;
    ctx.stroke();
    ctx.strokeStyle = latitude === '0' ? '#a21caf' : '#d946ef';
    ctx.lineWidth = latitude === '0' ? 2 : 1;
    ctx.stroke(); ctx.setLineDash([]);
  }
  ctx.fillStyle = '#334155'; ctx.font = '11px system-ui, sans-serif';
  ctx.textAlign = 'center';
  for (let lon = -180; lon <= 180; lon += 60)
    ctx.fillText(`${lon}°`, x(lon), plot.top + plot.height + 17);
  ctx.textAlign = 'right';
  for (let dec = -90; dec <= 90; dec += 30)
    ctx.fillText(`${dec}°`, plot.left - 6, y(dec) + 4);
}
function render() { updateValues(); drawMap(0); drawMap(1); }
slider.addEventListener('input', render);
window.addEventListener('resize', render);
canvases.forEach((canvas, index) => {
  const tip = document.getElementById(index ? 'new-tooltip' : 'old-tooltip');
  canvas.addEventListener('mousemove', event => {
    const rect = canvas.getBoundingClientRect();
    const mx = event.clientX - rect.left, my = event.clientY - rect.top;
    const plot = geometry[index];
    if (mx < plot.left || mx > plot.left + plot.width ||
        my < plot.top || my > plot.top + plot.height) { tip.style.display = 'none'; return; }
    const lon = (mx - plot.left) / plot.width * 360 - 180;
    const dec = 90 - (my - plot.top) / plot.height * 180;
    const idx = lists[index].findIndex(cell => dec >= cell.bottom && dec <= cell.top &&
      cell.parts.some(part => lon >= part[0] && lon <= part[1]));
    if (idx < 0) { tip.style.display = 'none'; return; }
    const cell = lists[index][idx], value = values[index][idx];
    tip.textContent = MODE === 'gaps'
      ? `${cell.name} · ${cell.gaps.length} gaps · ` +
        (value === null ? 'no measured gap' : `P${slider.value}: ${value.toFixed(2)} days`)
      : `${cell.name} · ${cell.visits} visits across ${cell.nights} nights`;
    tip.style.left = `${Math.min(mx + 12, rect.width - 210)}px`;
    tip.style.top = `${Math.max(5, my - 45)}px`;
    tip.style.display = 'block';
  });
  canvas.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
});
render();
</script>
</body>
</html>'''


def write_comparison_maps(old_fields, old_gaps, new_fields, new_gaps, output_dir):
    """Write paired maps with a shared scale and synchronized gap slider."""
    if set(old_fields['base_field']) != set(new_fields['base_field']):
        raise ValueError('Both plans must use the same field grid')
    new_fields = new_fields.set_index('base_field').loc[old_fields['base_field']].reset_index()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data = {'oldCells': grid_cells(old_fields, old_gaps),
            'newCells': grid_cells(new_fields, new_gaps),
            'galactic': galactic_lines(),
            'maxVisits': int(max(old_fields['visits'].max(), new_fields['visits'].max()))}
    payload = json.dumps(data, separators=(',', ':'))
    for mode, title, name in (
        ('visits', '30-night comparison · visits per field', 'visits_side_by_side.html'),
        ('gaps', '30-night comparison · visit-gap percentile', 'visit_gaps_side_by_side.html'),
    ):
        html = COMPARISON_HTML.replace('__TITLE__', title).replace('__DATA__', payload).replace('__MODE__', mode)
        path = output_dir / name
        path.write_text(html)
        print(f'  wrote {path}')
