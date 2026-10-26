// LS4 observing-plan site. Data are written by ls4/website.py (`ls4 site`).
'use strict';

const $ = id => document.getElementById(id);
const root = document.querySelector('.viz-root');
const css = name => getComputedStyle(root).getPropertyValue(name).trim();
const D2R = Math.PI / 180;
const CHILE = 'America/Santiago';
const wrap180 = a => ((a + 180) % 360 + 360) % 360 - 180;
const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
const mjdDate = mjd => new Date((mjd - 40587) * 864e5);
const hhmm = (mjd, timeZone = 'UTC') => mjdDate(mjd).toLocaleTimeString('en-GB',
  {timeZone, hour: '2-digit', minute: '2-digit', hourCycle: 'h23'});
const dayLabel = (iso, opts = {weekday: 'short', year: 'numeric'}) =>
  new Date(iso + 'T12:00:00Z').toLocaleDateString('en-US', {month: 'short', day: 'numeric', timeZone: 'UTC', ...opts});
const shortDay = mjd => mjdDate(mjd).toLocaleDateString('en-US', {month: 'short', day: 'numeric', timeZone: 'UTC'});
const stampIso = stamp => stamp.replace(/(\d{4})(\d\d)(\d\d)/, '$1-$2-$3');
const num = (v, digits = 0) => Number(v).toLocaleString('en-US',
  {maximumFractionDigits: digits, minimumFractionDigits: digits});
const raHours = ra => {
  const minutes = Math.round((((ra % 360) + 360) % 360) / 15 * 60) % 1440;
  return `${String(Math.floor(minutes / 60)).padStart(2, '0')}h${String(minutes % 60).padStart(2, '0')}m`;
};
const tonightStamp = () => new Date(Date.now() - 12 * 3600e3).toISOString().slice(0, 10).replaceAll('-', '');

function el(tag, attrs = {}, ...children) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') e.className = v; else if (k === 'style') e.style.cssText = v; else e.setAttribute(k, v);
  }
  for (const c of children) if (c !== null && c !== undefined) e.append(c);
  return e;
}
function svg(tag, attrs = {}, text) {
  const e = document.createElementNS('http://www.w3.org/2000/svg', tag);
  for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  return e;
}
const swatch = name => el('span', {class: 'swatch', style: `background:var(--${name})`});

const cache = new Map();
function load(path) {
  if (!cache.has(path)) {
    cache.set(path, fetch(path).then(r => {
      if (!r.ok) throw new Error(`${path}: HTTP ${r.status}`);
      return r.json();
    }));
  }
  return cache.get(path);
}

// ---- sky geometry -------------------------------------------------------------

const vec = (ra, dec) => [Math.cos(dec * D2R) * Math.cos(ra * D2R), Math.cos(dec * D2R) * Math.sin(ra * D2R),
                          Math.sin(dec * D2R)];
const radec = v => [((Math.atan2(v[1], v[0]) / D2R) + 360) % 360, Math.asin(clamp(v[2], -1, 1)) / D2R];
function basis(ra, dec) {
  const a = ra * D2R, d = dec * D2R;
  return {c: vec(ra, dec), east: [-Math.sin(a), Math.cos(a), 0],
          north: [-Math.sin(d) * Math.cos(a), -Math.sin(d) * Math.sin(a), Math.cos(d)]};
}
// Point at (dLon, dLat) degrees in the frame centered on (ra, dec): astropy's spherical_offsets_by,
// which ls4/bands.py uses for the quadrant corners.
function offset(ra, dec, dLon, dLat) {
  const {c, east, north} = basis(ra, dec), l = dLon * D2R, b = dLat * D2R;
  const x = Math.cos(b) * Math.cos(l), y = Math.cos(b) * Math.sin(l), z = Math.sin(b);
  return radec([0, 1, 2].map(k => x * c[k] + y * east[k] + z * north[k]));
}
function circle(ra, dec, radius, n = 120) {
  const {c, east, north} = basis(ra, dec), r = radius * D2R, points = [];
  for (let k = 0; k <= n; k++) {
    const t = 2 * Math.PI * k / n;
    points.push(radec([0, 1, 2].map(i => Math.cos(r) * c[i] + Math.sin(r) * (Math.cos(t) * north[i] + Math.sin(t) * east[i]))));
  }
  return points;
}

// Plate carrée map with RA increasing to the right, matching the camera diagram and QUADRANT_BANDS.
class SkyMap {
  constructor(canvas, draw) {
    this.canvas = canvas; this.ctx = canvas.getContext('2d'); this.draw = draw;
    this.center = 0; this.dec = INDEX.decRange; this.width = 0;
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
  }
  resize() {
    const w = this.canvas.clientWidth;
    if (!w || w === this.width) return;
    this.width = w;
    this.m = {left: 44, right: 12, top: 10, bottom: 34};
    this.scale = (w - this.m.left - this.m.right) / 360;
    this.plotH = (this.dec[1] - this.dec[0]) * this.scale;
    this.height = Math.round(this.m.top + this.m.bottom + this.plotH);
    this.canvas.style.height = `${this.height}px`;
    const ratio = window.devicePixelRatio || 1;
    this.canvas.width = Math.round(w * ratio); this.canvas.height = Math.round(this.height * ratio);
    this.ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    this.draw();
  }
  get ready() { return this.width > 0; }
  x(ra) { return this.m.left + (180 + wrap180(ra - this.center)) * this.scale; }
  y(dec) { return this.m.top + (this.dec[1] - dec) * this.scale; }
  toSky(px, py) {
    return [((this.center - 180 + (px - this.m.left) / this.scale) % 360 + 360) % 360,
            this.dec[1] - (py - this.m.top) / this.scale];
  }
  // Add a polygon or polyline to the current path, unwrapped in RA so it crosses the map edge cleanly.
  path(points, closed) {
    const ctx = this.ctx, d = [wrap180(points[0][0] - this.center)];
    for (let k = 1; k < points.length; k++) d.push(d[k - 1] + wrap180(points[k][0] - points[k - 1][0]));
    const lo = Math.min(...d), hi = Math.max(...d);
    for (const shift of [0, 360, -360]) {
      if (lo + shift >= 180 || hi + shift <= -180) continue;
      points.forEach((p, k) => {
        const X = this.m.left + (180 + d[k] + shift) * this.scale, Y = this.y(p[1]);
        if (k) ctx.lineTo(X, Y); else ctx.moveTo(X, Y);
      });
      if (closed) ctx.closePath();
    }
  }
  pixel(p, size) {
    const w = Math.min(360, size / Math.max(0.05, Math.cos(p[1] * D2R))) * this.scale, h = size * this.scale;
    const X = this.x(p[0]), Y = this.y(p[1]), right = this.m.left + 360 * this.scale;
    this.ctx.fillRect(X - w / 2, Y - h / 2, w + 0.6, h + 0.6);
    if (X - w / 2 < this.m.left) this.ctx.fillRect(X - w / 2 + 360 * this.scale, Y - h / 2, w + 0.6, h + 0.6);
    if (X + w / 2 > right) this.ctx.fillRect(X - w / 2 - 360 * this.scale, Y - h / 2, w + 0.6, h + 0.6);
  }
  begin() {
    const c = this.ctx, {left, top} = this.m, W = 360 * this.scale;
    c.setLineDash([]); c.globalAlpha = 1;
    c.fillStyle = css('--surface-1'); c.fillRect(0, 0, this.width, this.height);
    c.strokeStyle = css('--grid'); c.lineWidth = 1; c.beginPath();
    for (let h = 0; h < 24; h += 2) { const X = this.x(h * 15); c.moveTo(X, top); c.lineTo(X, top + this.plotH); }
    for (let dec = Math.ceil(this.dec[0] / 15) * 15; dec <= this.dec[1]; dec += 15) {
      c.moveTo(left, this.y(dec)); c.lineTo(left + W, this.y(dec));
    }
    c.stroke();
    c.save(); c.beginPath(); c.rect(left, top, W, this.plotH); c.clip();
  }
  galactic() {
    const c = this.ctx;
    c.beginPath(); this.path(INDEX.galactic, false);
    c.strokeStyle = css('--plane'); c.lineWidth = 1.5; c.setLineDash([]); c.stroke();
  }
  end() {
    const c = this.ctx, {left, top} = this.m, W = 360 * this.scale;
    c.restore();
    c.strokeStyle = css('--axis'); c.lineWidth = 1; c.setLineDash([]); c.strokeRect(left, top, W, this.plotH);
    c.fillStyle = css('--text-secondary'); c.font = '11px system-ui, sans-serif';
    c.textAlign = 'center';
    for (let h = 0; h < 24; h += 2) {
      const X = this.x(h * 15);
      if (X > left + 8 && X < left + W - 8) c.fillText(`${h}h`, X, top + this.plotH + 14);
    }
    c.fillText('Right ascension', left + W / 2, this.height - 4);
    c.textAlign = 'right';
    for (let dec = Math.ceil(this.dec[0] / 15) * 15; dec <= this.dec[1]; dec += 15) {
      c.fillText(`${dec}°`, left - 6, this.y(dec) + 4);
    }
  }
}

function rampColor(stops, t) {
  t = clamp(t, 0, 1) * (stops.length - 1);
  const a = stops[Math.floor(t)], b = stops[Math.ceil(t)], w = t - Math.floor(t);
  const ch = i => Math.round(parseInt(a.slice(i, i + 2), 16) * (1 - w) + parseInt(b.slice(i, i + 2), 16) * w);
  return `rgb(${ch(1)},${ch(3)},${ch(5)})`;
}
const stopsOf = name => css(name).split(',').map(s => s.trim());

// ---- state ---------------------------------------------------------------------

let INDEX;
const state = {collection: null, nightIdx: 0, view: 'night', color: 'time', underlay: true,
               band: 'i', metric: 'nights', night: null, survey: null, clock: 0, playing: false,
               pixel: null, loads: 0};
const nights = () => state.collection.nights;
const current = () => nights()[state.nightIdx];

function defaultNight(list) {
  const stamp = tonightStamp();
  let best = 0;
  list.forEach((n, i) => { if (n.night <= stamp) best = i; });
  return best;
}

function writeHash() {
  const p = new URLSearchParams({c: state.collection.key, n: current().night, v: state.view});
  if (state.view === 'survey') { p.set('band', state.band); p.set('metric', state.metric); }
  else if (state.color !== 'time') p.set('color', state.color);
  history.replaceState(null, '', `#${p}`);
}

function segmented(id, options, get, set) {
  const box = $(id);
  box.replaceChildren(...options.map(o => {
    const b = el('button', {type: 'button'}, o.swatch ? swatch(o.swatch) : null, o.label);
    b.dataset.key = o.key;
    b.addEventListener('click', () => set(o.key));
    return b;
  }));
  const sync = () => box.querySelectorAll('button').forEach(b =>
    b.setAttribute('aria-pressed', String(get() === b.dataset.key)));
  sync();
  return sync;
}

// ---- night view ------------------------------------------------------------------

function prepareNight(N) {
  const hx = INDEX.fov.ra / 2, hy = INDEX.fov.dec / 2;
  N.geo = N.pointings.map(p => {
    const g = {};
    for (const sx of [-1, 0, 1]) for (const sy of [-1, 0, 1]) g[`${sx},${sy}`] = offset(p[2], p[3], sx * hx, sy * hy);
    return {
      outline: [[-1, -1], [0, -1], [1, -1], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0]].map(([a, b]) => g[`${a},${b}`]),
      quads: INDEX.quadrants.map(([sx, sy, band]) => ({band, points: [g['0,0'], g[`${sx},0`], g[`${sx},${sy}`], g[`0,${sy}`]]})),
    };
  });
  const starts = N.pointings.map(p => p[0]), ends = N.pointings.map(p => p[1]);
  N.span = [Math.min(N.window[0], ...starts), Math.max(N.window[1], ...ends)];
  return N;
}

function skyAt(N, mjd) {
  const s = N.sky;
  let k = 0;
  while (k < s.length - 2 && s[k + 1][0] < mjd) k++;
  const a = s[k], b = s[k + 1] || a, w = b[0] > a[0] ? clamp((mjd - a[0]) / (b[0] - a[0]), 0, 1) : 0;
  const angle = i => (a[i] + w * wrap180(b[i] - a[i]) + 360) % 360;
  return {lst: angle(1), moonRa: angle(2), moonDec: a[3] + w * (b[3] - a[3])};
}

const clockMjd = () => state.night.span[0] + state.clock / 1440;
function visibleCount(mjd) {
  const p = state.night.pointings;
  let n = 0;
  while (n < p.length && p[n][0] <= mjd + 1e-9) n++;
  return n;
}

let priorCache = {key: null, pixels: []};
function priorPixels() {
  const S = state.survey, key = `${state.collection.key}/${state.nightIdx}`;
  if (!S) return [];
  if (priorCache.key !== key) {
    priorCache = {key, pixels: S.pixels.filter(p => [2, 3, 4].some(i => p[i].length && p[i][0] < state.nightIdx))};
  }
  return priorCache.pixels;
}

function drawNight() {
  const m = nightMap;
  if (!m.ready || state.view !== 'night') return;
  const c = m.ctx, N = state.night;
  m.begin();
  if (N) {
    if (state.underlay && state.survey) {
      c.fillStyle = css('--prior');
      for (const p of priorPixels()) m.pixel(p, state.survey.pixelDeg);
    }
    m.galactic();
    const t = clockMjd(), sky = skyAt(N, t);
    c.lineWidth = 1.2; c.setLineDash([5, 4]); c.strokeStyle = css('--text-secondary');
    c.beginPath(); m.path(circle(sky.lst, INDEX.site.lat, INDEX.airmassRadius), true); c.stroke();
    c.strokeStyle = css('--muted');
    c.beginPath(); m.path(circle(sky.moonRa, sky.moonDec, INDEX.moonSeparation), true); c.stroke();
    c.setLineDash([]);
    const zx = m.x(sky.lst), zy = m.y(INDEX.site.lat);
    c.strokeStyle = css('--text-secondary'); c.lineWidth = 1.5; c.beginPath();
    c.moveTo(zx - 6, zy); c.lineTo(zx + 6, zy); c.moveTo(zx, zy - 6); c.lineTo(zx, zy + 6); c.stroke();
    c.fillStyle = css('--muted'); c.beginPath(); c.arc(m.x(sky.moonRa), m.y(sky.moonDec), 6, 0, 2 * Math.PI); c.fill();

    const n = visibleCount(t), time = stopsOf('--time');
    if (state.color === 'bands') {
      c.globalAlpha = 0.42;
      const colors = Object.fromEntries(INDEX.bands.map(b => [b, css(`--${b}`)]));
      for (let k = 0; k < n; k++) for (const q of N.geo[k].quads) {
        c.fillStyle = colors[q.band]; c.beginPath(); m.path(q.points, true); c.fill();
      }
      c.globalAlpha = 0.45;
    }
    c.lineWidth = 1.3;
    for (let k = 0; k < n; k++) {
      const p = N.pointings[k];
      c.strokeStyle = state.color === 'bands' ? css('--text-secondary')
        : rampColor(time, (p[0] - N.span[0]) / (N.span[1] - N.span[0]));
      c.setLineDash(p[4] ? [4, 3] : []);
      c.beginPath(); m.path(N.geo[k].outline, true); c.stroke();
    }
    c.globalAlpha = 1; c.setLineDash([]);
    if (n > 0 && n < N.pointings.length) {
      c.lineWidth = 2.5; c.strokeStyle = css('--text-primary');
      c.beginPath(); m.path(N.geo[n - 1].outline, true); c.stroke();
    }
  }
  m.end();
}

function updateClock() {
  const N = state.night, t = clockMjd(), n = visibleCount(t);
  $('clock').value = state.clock;
  $('clock-label').value = `${hhmm(t)} UTC · ${hhmm(t, CHILE)} Chile · ${n} of ${N.pointings.length} exposures`;
  const rows = $('pointings').querySelectorAll('tbody tr');
  rows.forEach((r, k) => r.classList.toggle('current', k === n - 1));
  drawNight();
}

function tile(label, value, sub) {
  return el('div', {class: 'tile'}, el('div', {class: 'label'}, label), el('div', {class: 'value'}, value),
            sub ? el('div', {class: 'sub'}, sub) : null);
}

function renderNightInfo() {
  const N = state.night, M = N.metrics, [start, end] = N.window;
  const hours = (end - start) * 24;
  const tiles = [
    tile('Observing window (Chile)', `${hhmm(start, CHILE)}–${hhmm(end, CHILE)}`,
         `${hhmm(start)}–${hhmm(end)} UTC · ${num(hours, 1)} h between astronomical twilights`),
    tile('Exposures', num(M.images), `${num(M.complete_pairs)} initial + dither pairs · ${num(M.sole_pairs)} SOLE`),
    tile('Area tonight', `${num(N.area.any)} deg²`, `g ${num(N.area.g)} · i ${num(N.area.i)} · z ${num(N.area.z)} deg²`),
    tile('Dither phase', N.phase === null ? '—' : String(N.phase),
         N.phase === null ? 'Plan predates the alternating dither' : 'Alternates nightly, swapping g and z'),
    tile('Unused time', `${num(M.unused_minutes)} min`, `Median slew ${num(M.median_slew_deg, 1)}°, max ${num(M.max_slew_deg, 1)}°`),
    tile('Moon', `${num(N.moonIllumination * 100)}%`, `illuminated; pointings stay ≥ ${num(INDEX.moonSeparation)}° away`),
  ];
  const w = N.release && N.release.weather;
  if (w) {
    tiles.push(tile('Weather forecast', `${num(w.cloud_cover_max_percent)}% cloud`,
      `max; rain chance ${num(w.precipitation_probability_max_percent)}% · humidity ${num(w.humidity_max_percent)}% · ` +
      `${num(w.temperature_min_c, 1)} to ${num(w.temperature_max_c, 1)} °C`));
  }
  $('tiles').replaceChildren(...tiles);

  const links = $('links');
  links.replaceChildren();
  if (N.release) {
    const warnings = N.release.warnings || [];
    links.append(el('span', {class: 'status'},
      el('span', {class: 'dot', style: `background:var(--${warnings.length ? 'warn' : 'good'})`}),
      warnings.length ? 'Published with warnings' : 'Published'));
    links.append(el('a', {href: N.release.url}, `Release ${N.release.tag}`));
  } else if (state.collection.key === 'production') {
    links.append(el('span', {class: 'status'}, 'Saved plan from the repository (plans/)'));
  } else {
    links.append(el('span', {class: 'status'}, 'Forecast; not an operational plan'));
  }
  links.append(el('a', {href: N.files.csv, download: ''}, 'Plan CSV'));
  links.append(el('a', {href: N.files.obsplan, download: ''}, `Telescope file (${N.files.obsplan.split('/').pop()})`));
  if (N.release && N.release.gif) links.append(el('a', {href: N.release.gif}, 'Night GIF'));
  if (N.release && (N.release.warnings || []).length) {
    links.append(el('ul', {class: 'warnings'}, ...N.release.warnings.map(w => el('li', {}, w))));
  }

  const legend = $('night-legend');
  const item = (...children) => el('span', {class: 'item'}, ...children);
  const items = [];
  if (state.color === 'bands') {
    items.push(...INDEX.bands.map(b => item(swatch(b), `${b} quadrant`)));
  } else {
    items.push(item(`${hhmm(N.span[0])} UTC`),
      el('span', {class: 'ramp', style: `background:linear-gradient(to right, ${stopsOf('--time').join(',')})`}),
      item(`${hhmm(N.span[1])} UTC`));
  }
  items.push(item(el('span', {class: 'line-key'}), 'First pointing'),
    item(el('span', {class: 'line-key dashed'}), 'Dithered revisit'),
    item(el('span', {class: 'swatch', style: 'background:var(--prior)'}), 'Observed on earlier nights'),
    item(el('span', {class: 'line-key dashed'}), `Airmass ${INDEX.maxAirmass} limit (around the zenith +)`),
    item(el('span', {class: 'swatch', style: 'background:var(--muted);border-radius:50%'}), `Moon, ${num(INDEX.moonSeparation)}° exclusion`),
    item(el('span', {class: 'line-key', style: 'border-color:var(--plane)'}), 'Galactic plane'));
  legend.replaceChildren(...items);

  const summary = $('pointings').querySelector('summary');
  summary.textContent = `Pointings (${num(N.pointings.length)})`;
  $('pointings').querySelector('tbody').replaceChildren();
  if ($('pointings').open) fillTable();
}

function fillTable() {
  const N = state.night, body = $('pointings').querySelector('tbody');
  if (body.childElementCount) return;
  body.append(...N.pointings.map((p, k) => {
    const row = el('tr', {}, ...[k + 1, hhmm(p[0]), hhmm(p[0], CHILE), N.targets[k], p[4] ? 'Dither' : 'Initial',
      p[5], raHours(p[2]), `${num(p[3], 2)}°`, p[6] ? 'SOLE' : 'Extragalactic'].map(v => el('td', {}, String(v))));
    row.addEventListener('click', () => { stopPlay(); state.clock = (p[0] - N.span[0]) * 1440 + 1e-3; updateClock(); });
    return row;
  }));
  updateClock();
}

async function showNight() {
  const token = ++state.loads, info = current();
  const N = prepareNight(await load(`data/${state.collection.key}/${info.night}.json`));
  if (token !== state.loads) return;
  stopPlay();
  state.night = N;
  const minutes = Math.round((N.span[1] - N.span[0]) * 1440);
  $('clock').max = minutes;
  state.clock = minutes;
  nightMap.center = skyAt(N, (N.window[0] + N.window[1]) / 2).lst;
  renderNightInfo();
  updateClock();
}

let playFrame = null;
function stopPlay() {
  state.playing = false; $('play').textContent = 'Play';
  if (playFrame) cancelAnimationFrame(playFrame);
  playFrame = null;
}
function togglePlay() {
  if (state.playing) return stopPlay();
  const max = Number($('clock').max);
  if (state.clock >= max) state.clock = 0;
  state.playing = true; $('play').textContent = 'Pause';
  let last = performance.now();
  const step = now => {
    state.clock = Math.min(max, state.clock + (now - last) / 12000 * max);  // whole night in ~12 s
    last = now;
    updateClock();
    if (state.clock >= max) return stopPlay();
    playFrame = requestAnimationFrame(step);
  };
  playFrame = requestAnimationFrame(step);
}

function nightTooltip(event) {
  const N = state.night, tip = $('night-tooltip');
  if (!N) return;
  const rect = nightMap.canvas.getBoundingClientRect(), mx = event.clientX - rect.left, my = event.clientY - rect.top;
  const [ra, dec] = nightMap.toSky(mx, my), hx = INDEX.fov.ra / 2, hy = INDEX.fov.dec / 2;
  const hits = [];
  for (let k = visibleCount(clockMjd()) - 1; k >= 0; k--) {
    const p = N.pointings[k];
    if (Math.abs(dec - p[3]) < hy && Math.abs(wrap180(ra - p[2]) * Math.cos(p[3] * D2R)) < hx) hits.push(k);
  }
  if (!hits.length) { tip.style.display = 'none'; return; }
  tip.replaceChildren(el('strong', {}, hits.length === 1 ? N.targets[hits[0]] : `${hits.length} exposures here`));
  for (const k of hits.slice(0, 6)) {
    const p = N.pointings[k];
    tip.append(el('div', {class: 'row'}, `${hhmm(p[0])} UTC · ${N.targets[k]} · ${p[4] ? 'dithered revisit' : 'first pointing'}` +
      ` · block ${p[5]}${p[6] ? ' · SOLE' : ''}`));
  }
  if (hits.length > 6) tip.append(el('div', {class: 'row'}, `and ${hits.length - 6} earlier`));
  placeTooltip(tip, mx, my, nightMap.width);
}

function placeTooltip(tip, mx, my, width) {
  tip.style.display = 'block';
  const w = tip.offsetWidth;
  tip.style.left = `${mx + 16 + w > width ? Math.max(4, mx - 16 - w) : mx + 16}px`;
  tip.style.top = `${Math.max(4, my - 20)}px`;
}

// ---- survey view -------------------------------------------------------------------

const METRICS = {
  nights: {label: 'Nights observed', unit: 'nights', empty: 'observed only in other bands'},
  gap: {label: 'Median gap', unit: 'days between observed nights', empty: 'fewer than two nights'},
  since: {label: 'Days since last', unit: 'days since last observation', empty: 'observed only in other bands'},
};
const BAND_INDEX = {g: 2, i: 3, z: 4};

function bandNights(pixel, band, through = state.nightIdx) {
  const list = pixel[BAND_INDEX[band]];
  let n = 0;
  while (n < list.length && list[n] <= through) n++;
  return list.slice(0, n);
}
function gapsOf(nightList) {
  const mjd = state.survey.nights;
  return nightList.slice(1).map((n, k) => mjd[n].mjd - mjd[nightList[k]].mjd);
}
function median(values) {
  if (!values.length) return null;
  const s = [...values].sort((a, b) => a - b), mid = (s.length - 1) / 2;
  return (s[Math.floor(mid)] + s[Math.ceil(mid)]) / 2;
}
function metricFor(pixel, band, metric = state.metric) {
  const list = bandNights(pixel, band);
  if (metric === 'nights') return list.length || null;
  if (metric === 'since') return list.length
    ? state.survey.nights[state.nightIdx].mjd - state.survey.nights[list[list.length - 1]].mjd : null;
  return list.length < 2 ? null : median(gapsOf(list));
}
function pixelValue(pixel) {
  if (state.band !== 'all') return metricFor(pixel, state.band);
  if (state.metric === 'nights') return INDEX.bands.filter(b => bandNights(pixel, b).length).length || null;
  const per = INDEX.bands.map(b => metricFor(pixel, b));  // the slowest band sets the cadence
  return per.some(v => v === null) ? null : Math.max(...per);
}

function surveyCenter(S) {
  // Put the map's edge in the emptiest 30° of right ascension.
  const counts = new Array(360).fill(0);
  for (const p of S.pixels) counts[Math.floor(p[0]) % 360]++;
  let best = 0, bestCount = Infinity;
  for (let start = 0; start < 360; start++) {
    let total = 0;
    for (let k = 0; k < 30; k++) total += counts[(start + k) % 360];
    if (total < bestCount) { bestCount = total; best = start; }
  }
  return (best + 15 + 180) % 360;
}

let surveyScale = {max: 1};
function drawSurvey() {
  const m = surveyMap, S = state.survey;
  if (!m.ready || state.view !== 'survey') return;
  m.begin();
  if (S) {
    const c = m.ctx, ramp = stopsOf('--ramp'), empty = css('--empty');
    const values = S.pixels.map(pixelValue);
    const shown = values.filter(v => v !== null).sort((a, b) => a - b);
    const allNights = state.band === 'all' && state.metric === 'nights';
    const max = allNights ? 3 : state.metric === 'nights'
      ? Math.max(1, shown.length ? shown[shown.length - 1] : 1)
      : (shown.length ? shown[Math.floor(0.99 * (shown.length - 1))] : 1) || 1;
    surveyScale = {max, allNights};
    S.pixels.forEach((p, k) => {
      const v = values[k];
      if (v === null && !INDEX.bands.some(b => bandNights(p, b).length)) return;
      c.fillStyle = v === null ? empty : rampColor(ramp, allNights ? (v - 1) / 2 : v / max);
      m.pixel(p, S.pixelDeg);
    });
    m.galactic();
    if (state.pixel !== null) {
      const p = S.pixels[state.pixel], w = S.pixelDeg / Math.max(0.05, Math.cos(p[1] * D2R)) * m.scale;
      c.strokeStyle = css('--text-primary'); c.lineWidth = 2;
      c.strokeRect(m.x(p[0]) - w / 2 - 3, m.y(p[1]) - S.pixelDeg * m.scale / 2 - 3, w + 6, S.pixelDeg * m.scale + 6);
    }
  }
  m.end();
  renderSurveyLegend();
}

function renderSurveyLegend() {
  const {max, allNights} = surveyScale, M = METRICS[state.metric];
  const item = (...children) => el('span', {class: 'item'}, ...children);
  $('survey-legend').replaceChildren(
    item(allNights ? '1' : '0'),
    el('span', {class: 'ramp', style: `background:linear-gradient(to right, ${stopsOf('--ramp').join(',')})`}),
    item(state.metric === 'nights' ? String(max) : num(max, 1)),
    item(allNights ? 'of g, i, z seen' : state.band === 'all' ? `${M.unit} (slowest band)` : M.unit),
    item(el('span', {class: 'swatch', style: 'background:var(--empty)'}), M.empty),
    item(el('span', {class: 'line-key', style: 'border-color:var(--plane)'}), 'Galactic plane'),
    ...(state.metric === 'gap' ? [item('Scale tops out at the 99th percentile.')] : []));
}

function findPixel(event) {
  const S = state.survey;
  const rect = surveyMap.canvas.getBoundingClientRect(), mx = event.clientX - rect.left, my = event.clientY - rect.top;
  const [ra, dec] = surveyMap.toSky(mx, my);
  let best = null, bestD = S.pixelDeg * 0.8;
  for (let da = -3; da <= 3; da++) for (let dd = -2; dd <= 2; dd++) {
    for (const k of S.buckets.get(`${(Math.floor(ra) + da + 360) % 360},${Math.floor(dec) + dd}`) || []) {
      const p = S.pixels[k];
      const d = Math.hypot(wrap180(p[0] - ra) * Math.cos(dec * D2R), p[1] - dec);
      if (d < bestD) { bestD = d; best = k; }
    }
  }
  return {pixel: best, mx, my};
}

function bandSummary(p, band) {
  const list = bandNights(p, band), gap = list.length > 1 ? median(gapsOf(list)) : null;
  const since = list.length ? state.survey.nights[state.nightIdx].mjd - state.survey.nights[list[list.length - 1]].mjd : null;
  return `${list.length} night${list.length === 1 ? '' : 's'}` + (gap === null ? '' : ` · median gap ${num(gap, 1)} d`) +
    (since === null ? '' : ` · last ${num(since, 1)} d before`);
}

function surveyTooltip(event) {
  const tip = $('survey-tooltip');
  if (!state.survey) return;
  const {pixel, mx, my} = findPixel(event);
  if (pixel === null) { tip.style.display = 'none'; return; }
  const p = state.survey.pixels[pixel];
  tip.replaceChildren(el('strong', {}, `RA ${raHours(p[0])}, Dec ${num(p[1], 1)}°`),
    ...INDEX.bands.map(b => el('div', {class: 'row'}, el('span', {class: 'key', style: `border-color:var(--${b})`}),
                                `${b}: ${bandSummary(p, b)}`)));
  placeTooltip(tip, mx, my, surveyMap.width);
}

function renderPixelDetail() {
  const box = $('pixel-detail'), S = state.survey;
  if (!S || state.pixel === null) {
    box.replaceChildren('Click a part of the sky to see when it was observed in each band.');
    return;
  }
  const p = S.pixels[state.pixel], first = S.nights[0].mjd, last = S.nights[S.nights.length - 1].mjd;
  const through = S.nights[state.nightIdx].mjd;
  const rows = INDEX.bands.map(band => {
    const strip = svg('svg', {viewBox: '0 0 1000 18', preserveAspectRatio: 'none', role: 'img',
                              'aria-label': `${band}-band nights observed`});
    const X = mjd => 6 + (mjd - first) / Math.max(1, last - first) * 988;
    strip.append(svg('line', {x1: 0, x2: 1000, y1: 9, y2: 9, stroke: 'var(--grid)', 'stroke-width': 2}));
    for (const n of p[BAND_INDEX[band]]) {
      strip.append(svg('rect', {x: X(S.nights[n].mjd) - 2, y: 2, width: 4, height: 14, rx: 1,
        fill: `var(--${band})`, opacity: n <= state.nightIdx ? 1 : 0.3}));
    }
    strip.append(svg('line', {x1: X(through), x2: X(through), y1: 0, y2: 18, stroke: 'var(--text-primary)',
                              'stroke-width': 2, 'stroke-dasharray': '3 2'}));
    return el('div', {class: 'band-row'}, el('span', {class: 'name'}, swatch(band), band), strip,
              el('span', {}, bandSummary(p, band)));
  });
  const label = k => dayLabel(stampIso(S.nights[k].label), {});
  box.replaceChildren(el('h2', {}, `RA ${raHours(p[0])}, Dec ${num(p[1], 1)}° · nights observed, ` +
    `${label(0)} to ${label(S.nights.length - 1)} (dashed: selected night; faded: later nights)`), ...rows);
}

// ---- small charts --------------------------------------------------------------------

function niceTicks(lo, hi, count = 4) {
  const raw = (hi - lo) / count || 1, mag = 10 ** Math.floor(Math.log10(raw));
  const step = [1, 2, 5, 10].map(k => k * mag).find(s => s >= raw);
  const ticks = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) ticks.push(Number(v.toFixed(10)));
  return ticks;
}
const compact = v => Math.abs(v) >= 1000 ? `${num(v / 1000, v % 1000 ? 1 : 0)}k` : num(v, Math.abs(v) < 1 && v ? 2 : 0);

function frame(container, legend) {
  const width = container.clientWidth || 400, height = 240;
  const m = {left: 44, right: legend ? 48 : 12, top: 10, bottom: 26};
  const s = svg('svg', {viewBox: `0 0 ${width} ${height}`, role: 'img'});
  const tip = el('div', {class: 'tooltip'});
  container.replaceChildren(s, tip);
  return {s, tip, m, width, height};
}

function dateTicks(lo, hi, count) {
  const step = [1, 2, 7, 14, 28, 56, 112, 364].find(s => (hi - lo) / s <= count) || 364;
  const ticks = [];
  for (let v = Math.ceil(lo); v <= hi; v += step) ticks.push(v);  // UTC midnight, early in the Chilean evening
  return ticks;
}

function axes(f, X, Y, xDomain, yDomain, xFormat, dates = false, yFormat = compact) {
  const {s, m, width, height} = f;
  for (const v of niceTicks(yDomain[0], yDomain[1])) {
    s.append(svg('line', {class: 'gridline', x1: m.left, x2: width - m.right, y1: Y(v), y2: Y(v)}));
    s.append(svg('text', {x: m.left - 6, y: Y(v) + 4, 'text-anchor': 'end'}, yFormat(v)));
  }
  s.append(svg('line', {class: 'baseline', x1: m.left, x2: width - m.right, y1: height - m.bottom, y2: height - m.bottom}));
  const count = Math.max(2, Math.floor((width - m.left - m.right) / 70));
  const ticks = (dates ? dateTicks : niceTicks)(xDomain[0], xDomain[1], count);
  for (const v of ticks) s.append(svg('text', {x: X(v), y: height - m.bottom + 16, 'text-anchor': 'middle'}, xFormat(v)));
}

function hoverLayer(f, onMove) {
  const {s, m, width, height, tip} = f;
  const hair = svg('line', {class: 'crosshair', y1: m.top, y2: height - m.bottom, visibility: 'hidden'});
  const hit = svg('rect', {x: m.left, y: 0, width: width - m.left - m.right, height, fill: 'transparent'});
  s.append(hair, hit);
  hit.addEventListener('pointermove', event => {
    const box = s.getBoundingClientRect(), px = (event.clientX - box.left) * width / box.width;
    const x = onMove(px);
    if (x === null) return;
    hair.setAttribute('x1', x); hair.setAttribute('x2', x); hair.setAttribute('visibility', 'visible');
    placeTooltip(tip, x * box.width / width, (event.clientY - box.top), box.width);
  });
  hit.addEventListener('pointerleave', () => { hair.setAttribute('visibility', 'hidden'); tip.style.display = 'none'; });
  return hit;
}

function tipRows(tip, title, rows) {
  tip.replaceChildren(el('div', {}, title), ...rows.map(([color, value, label]) => el('div', {class: 'row'},
    el('span', {class: 'key', style: `border-color:${color}`}), el('strong', {}, value), label)));
}

function lineChart(container, {series, xDomain, yDomain, xFormat, yFormat, marker, snap, valueAt, title, dates}) {
  const f = frame(container, true), {s, m, width, height} = f;
  const X = v => m.left + (v - xDomain[0]) / (xDomain[1] - xDomain[0] || 1) * (width - m.left - m.right);
  const Y = v => height - m.bottom - (v - yDomain[0]) / (yDomain[1] - yDomain[0] || 1) * (height - m.top - m.bottom);
  axes(f, X, Y, xDomain, yDomain, xFormat, dates, yFormat);
  if (marker !== undefined) s.append(svg('line', {class: 'selected-line', x1: X(marker), x2: X(marker), y1: m.top, y2: height - m.bottom}));
  // Direct labels at line ends, nudged apart so they never overlap.
  const ends = [];
  for (const line of series) {
    if (!line.points.length) continue;
    const d = line.points.map(([x, y], k) => `${k ? 'L' : 'M'}${X(x).toFixed(1)},${Y(y).toFixed(1)}`).join('');
    s.append(svg('path', {d, fill: 'none', stroke: line.color, 'stroke-width': 2, 'stroke-linejoin': 'round'}));
    const lastPoint = line.points[line.points.length - 1];
    ends.push({y: Y(lastPoint[1]), x: X(lastPoint[0]), name: line.name});
  }
  ends.sort((a, b) => a.y - b.y).forEach((e, k, all) => { if (k && e.y - all[k - 1].y < 12) e.y = all[k - 1].y + 12; });
  for (const e of ends) s.append(svg('text', {x: e.x + 6, y: e.y + 4}, e.name));
  hoverLayer(f, px => {
    const x = snap(xDomain[0] + (px - m.left) / (width - m.left - m.right) * (xDomain[1] - xDomain[0]));
    if (x === null) return null;
    tipRows(f.tip, title(x), series.map(line => [line.color, valueAt(line, x), line.name]));
    return X(x);
  });
}

function renderCharts() {
  const S = state.survey;
  if (!S || state.view !== 'survey') return;
  const mjds = S.nights.map(n => n.mjd), first = mjds[0], last = mjds[mjds.length - 1];
  const through = mjds[state.nightIdx];
  const colors = Object.fromEntries([...INDEX.bands, 'all'].map(b => [b, css(`--${b}`)]));
  const nearestNight = x => mjds.reduce((best, v, k) => Math.abs(v - x) < Math.abs(mjds[best] - x) ? k : best, 0);
  const dateTick = v => shortDay(v - 0.5);  // label by local evening, as nights are named

  // Cumulative area per band and in all three.
  const areaSeries = [...INDEX.bands, 'all'].map(b => ({name: b === 'all' ? 'g+i+z' : b, key: b, color: colors[b],
    points: S.cumulative.map((row, k) => [mjds[k], row[b]])}));
  lineChart($('area-chart'), {
    series: areaSeries, xDomain: [first, Math.max(last, first + 1)],
    yDomain: [0, Math.max(...S.cumulative.map(r => Math.max(r.g, r.i, r.z))) * 1.05 || 1],
    xFormat: dateTick, marker: through, dates: true,
    snap: x => mjds[nearestNight(x)],
    valueAt: (line, x) => `${num(line.points[nearestNight(x)][1])} deg²`,
    title: x => dayLabel(stampIso(S.cumulative[nearestNight(x)].night)),
  });

  // Exposures per night; click to open that night.
  const list = nights(), f = frame($('exposure-chart'), false), {s, m, width, height} = f;
  const span = Math.max(last - first, 1);
  const X = v => m.left + 6 + (v - first) / span * (width - m.left - m.right - 12);
  const yMax = Math.max(...list.map(n => n.images)) * 1.05 || 1;
  const Y = v => height - m.bottom - v / yMax * (height - m.top - m.bottom);
  axes(f, X, Y, [first, last], [0, yMax], dateTick, true);
  const spacing = mjds.length > 1 ? Math.min(...mjds.slice(1).map((v, k) => v - mjds[k])) : span;
  const barW = clamp(spacing / span * (width - m.left - m.right) - 2, 2, 24);
  list.forEach((n, k) => {
    const bar = svg('rect', {x: X(mjds[k]) - barW / 2, y: Y(n.images), width: barW, height: Y(0) - Y(n.images),
      rx: Math.min(2, barW / 2), fill: k === state.nightIdx ? 'var(--text-primary)' : 'var(--bar)', cursor: 'pointer'});
    const hit = svg('rect', {x: X(mjds[k]) - Math.max(barW, 10) / 2, y: m.top, width: Math.max(barW, 10),
      height: height - m.top - m.bottom, fill: 'transparent', cursor: 'pointer', tabindex: 0});
    const show = event => {
      tipRows(f.tip, dayLabel(n.date), [['var(--bar)', num(n.images), ' exposures'],
        ['var(--bar)', num(n.pairs), ' pairs'], ['var(--bar)', `${num(n.area.any)} deg²`, ' area']]);
      const box = s.getBoundingClientRect();
      placeTooltip(f.tip, X(mjds[k]) * box.width / width, event.clientY ? event.clientY - box.top : 40, box.width);
    };
    hit.addEventListener('pointermove', show);
    hit.addEventListener('focus', show);
    hit.addEventListener('pointerleave', () => { f.tip.style.display = 'none'; });
    hit.addEventListener('click', () => { selectNight(k); setView('night'); });
    s.append(bar, hit);
  });

  // Distribution of gaps between observed nights, through the selected night.
  const gapSeries = INDEX.bands.map(b => {
    const gaps = [];
    for (const p of S.pixels) gaps.push(...gapsOf(bandNights(p, b)));
    gaps.sort((a, b) => a - b);
    return {name: b, color: colors[b], gaps};
  });
  const all = gapSeries.flatMap(g => g.gaps).sort((a, b) => a - b);
  const xMax = Math.max(4, all.length ? Math.ceil(all[Math.floor(0.95 * (all.length - 1))]) : 4);
  const cdfAt = (gaps, x) => {
    let lo = 0, hi = gaps.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (gaps[mid] <= x) lo = mid + 1; else hi = mid; }
    return gaps.length ? lo / gaps.length : 0;
  };
  const gapContainer = $('gap-chart');
  if (!all.length) {
    gapContainer.replaceChildren(el('p', {class: 'note'}, 'No pixel has been observed on two nights yet.'));
    return;
  }
  const grid = Array.from({length: 201}, (_, k) => k / 200 * xMax);
  lineChart(gapContainer, {
    series: gapSeries.map(g => ({...g, points: grid.map(x => [x, cdfAt(g.gaps, x) * 100])})),
    xDomain: [0, xMax], yDomain: [0, 100], yFormat: v => `${v}%`, xFormat: v => `${num(v, v < 10 && v % 1 ? 1 : 0)} d`,
    snap: x => clamp(Math.round(x * 2) / 2, 0, xMax),
    valueAt: (line, x) => `${num(cdfAt(line.gaps, x) * 100)}%`,
    title: x => `Gaps of ${num(x, 1)} days or less`,
  });
}

// ---- wiring ------------------------------------------------------------------------

let nightMap, surveyMap, syncColor, syncBand, syncMetric, syncView, syncCollection;

function setView(view) {
  state.view = view;
  $('night-view').hidden = view !== 'night';
  $('survey-view').hidden = view !== 'survey';
  syncView();
  writeHash();
  if (view === 'night') {
    nightMap.width = 0; nightMap.resize();
    if (!state.night || state.night.night !== current().night) showNight();
  } else {
    stopPlay();
    surveyMap.width = 0; surveyMap.resize();
    renderCharts(); renderPixelDetail();
  }
}

function selectNight(k) {
  state.nightIdx = clamp(k, 0, nights().length - 1);
  $('night').value = String(state.nightIdx);
  $('prev').disabled = state.nightIdx === 0;
  $('next').disabled = state.nightIdx === nights().length - 1;
  $('through').value = state.nightIdx;
  const info = current();
  $('through-label').value = `${dayLabel(info.date)} (night ${state.nightIdx + 1} of ${nights().length})`;
  writeHash();
  if (state.view === 'night') showNight();
  else { drawSurvey(); renderPixelDetail(); renderCharts(); }
}

async function selectCollection(key, nightStamp) {
  state.collection = INDEX.collections.find(c => c.key === key) || INDEX.collections[0];
  syncCollection();
  $('collection-note').textContent = state.collection.description;
  const stamp = tonightStamp();
  $('night').replaceChildren(...nights().map((n, k) => el('option', {value: String(k)},
    `${dayLabel(n.date)}${n.night === stamp ? ' (tonight)' : ''} · ${num(n.images)} exposures`)));
  $('through').max = nights().length - 1;
  state.survey = null; state.pixel = null; state.night = null;
  const k = nights().findIndex(n => n.night === nightStamp);
  selectNight(k >= 0 ? k : defaultNight(nights()));
  const collection = state.collection, S = await load(`data/${collection.key}/survey.json`);
  if (state.collection !== collection) return;  // another collection was picked meanwhile
  S.buckets = S.buckets || S.pixels.reduce((map, p, k2) => {
    const b = `${Math.floor(p[0])},${Math.floor(p[1])}`;
    if (!map.has(b)) map.set(b, []);
    map.get(b).push(k2);
    return map;
  }, new Map());
  surveyMap.center = surveyCenter(S);
  state.survey = S;
  drawNight(); drawSurvey(); renderPixelDetail(); renderCharts();
}

// URL hash, e.g. #c=production&n=20261026&v=survey&band=all&metric=gap
function readHash() {
  const hash = new URLSearchParams(location.hash.slice(1));
  state.view = hash.get('v') === 'survey' ? 'survey' : 'night';
  state.color = hash.get('color') === 'bands' ? 'bands' : 'time';
  if (['g', 'i', 'z', 'all'].includes(hash.get('band'))) state.band = hash.get('band');
  if (METRICS[hash.get('metric')]) state.metric = hash.get('metric');
  return hash;
}

async function applyHash() {
  const hash = readHash();
  syncColor(); syncBand(); syncMetric();
  $('night-view').hidden = state.view !== 'night';
  $('survey-view').hidden = state.view !== 'survey';
  syncView();
  const key = hash.get('c') || INDEX.collections[0].key;
  if (state.collection && state.collection.key === key) {
    const k = nights().findIndex(n => n.night === hash.get('n'));
    setView(state.view);
    if (k >= 0 && k !== state.nightIdx) selectNight(k);
    else if (state.night) { renderNightInfo(); drawNight(); }
  } else {
    await selectCollection(key, hash.get('n'));
  }
}

async function main() {
  INDEX = await load('data/index.json');
  readHash();

  nightMap = new SkyMap($('night-map'), drawNight);
  surveyMap = new SkyMap($('survey-map'), drawSurvey);

  syncCollection = segmented('collection', INDEX.collections.map(c => ({key: c.key, label: c.label})),
    () => state.collection && state.collection.key, key => selectCollection(key, current().night));
  syncView = segmented('view', [{key: 'night', label: 'Night plan'}, {key: 'survey', label: 'Coverage & cadence'}],
    () => state.view, setView);
  syncColor = segmented('color', [{key: 'time', label: 'Time'}, {key: 'bands', label: 'g/i/z quadrants'}],
    () => state.color, key => { state.color = key; syncColor(); writeHash(); renderNightInfo(); drawNight(); });
  syncBand = segmented('band', [...INDEX.bands.map(b => ({key: b, label: b, swatch: b})), {key: 'all', label: 'g + i + z'}],
    () => state.band, key => { state.band = key; syncBand(); writeHash(); drawSurvey(); });
  syncMetric = segmented('metric', Object.entries(METRICS).map(([key, m]) => ({key, label: m.label})),
    () => state.metric, key => { state.metric = key; syncMetric(); writeHash(); drawSurvey(); });

  $('night').addEventListener('change', e => selectNight(Number(e.target.value)));
  $('prev').addEventListener('click', () => selectNight(state.nightIdx - 1));
  $('next').addEventListener('click', () => selectNight(state.nightIdx + 1));
  $('tonight').addEventListener('click', () => selectNight(defaultNight(nights())));
  $('through').addEventListener('input', e => selectNight(Number(e.target.value)));
  $('underlay').addEventListener('change', e => { state.underlay = e.target.checked; drawNight(); });
  $('clock').addEventListener('input', e => { stopPlay(); state.clock = Number(e.target.value); updateClock(); });
  $('play').addEventListener('click', togglePlay);
  $('pointings').addEventListener('toggle', () => { if ($('pointings').open) fillTable(); });
  $('night-map').addEventListener('mousemove', nightTooltip);
  $('night-map').addEventListener('mouseleave', () => { $('night-tooltip').style.display = 'none'; });
  $('survey-map').addEventListener('mousemove', surveyTooltip);
  $('survey-map').addEventListener('mouseleave', () => { $('survey-tooltip').style.display = 'none'; });
  $('survey-map').addEventListener('click', event => {
    if (!state.survey) return;
    state.pixel = findPixel(event).pixel;
    drawSurvey(); renderPixelDetail();
  });
  document.addEventListener('keydown', event => {
    if (event.target.matches('input, select, textarea')) return;
    if (event.key === 'ArrowLeft') selectNight(state.nightIdx - 1);
    if (event.key === 'ArrowRight') selectNight(state.nightIdx + 1);
  });
  let chartWidth = 0;
  new ResizeObserver(entries => {
    const w = entries[0].contentRect.width;
    if (w !== chartWidth) { chartWidth = w; renderCharts(); }
  }).observe(document.querySelector('.charts'));
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
    if (state.night) renderNightInfo();
    drawNight(); drawSurvey(); renderCharts();
  });

  const commit = INDEX.commit ? INDEX.commit.slice(0, 7) : null;
  $('footer').replaceChildren(`${INDEX.historyNote} Generated ${INDEX.generatedUtc.replace('T', ' ').replace('+00:00', ' UTC')}`,
    commit ? ' from ' : '', commit ? el('a', {href: `https://github.com/${INDEX.repository}/commit/${INDEX.commit}`}, commit) : '',
    '. ', el('a', {href: `https://github.com/${INDEX.repository}`}, 'Source and documentation'), '.');

  window.addEventListener('hashchange', applyHash);
  await applyHash();
}

main().catch(error => {
  document.querySelector('main').append(el('p', {class: 'note'}, `Could not load the site data: ${error.message}`));
  throw error;
});
