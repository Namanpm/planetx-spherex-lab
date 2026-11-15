'use strict';
/* SPHEREx Epoch Lab front end. Plain DOM + SVG + canvas; every number shown comes from the API. */

const $ = (id) => document.getElementById(id);
const el = (tag, attrs = {}, ...children) => {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'class') node.className = v;
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) if (c !== null && c !== undefined && c !== false) node.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return node;
};
const SVG = 'http://www.w3.org/2000/svg';
const sv = (tag, attrs = {}, ...children) => {
  const node = document.createElementNS(SVG, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== undefined && v !== null) node.setAttribute(k, v);
  for (const c of children.flat()) if (c) node.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return node;
};
const fmt = (v, d = 3) => (v === null || v === undefined || Number.isNaN(v)) ? '—' : (typeof v === 'number' ? (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(2) : v.toFixed(d)) : String(v));
const short = (id) => (id || '').split('?').pop();
const mjdToDate = (mjd) => new Date((mjd - 40587) * 86400000);
const isoDay = (mjd) => mjdToDate(mjd).toISOString().slice(0, 10);

async function api(path, options = {}) {
  const response = await fetch(path, options);
  let body = null;
  try { body = await response.json(); } catch (_) { /* non-JSON */ }
  if (!response.ok) {
    const d = body && body.detail;
    let message = typeof d === 'string' ? d : (d && (d.message || d.code)) || `HTTP ${response.status}`;
    if (Array.isArray(d)) message = d.map((e) => e.msg).join('; ');
    const err = new Error(message); err.status = response.status; err.code = d && d.code; throw err;
  }
  return body;
}

/* ------------------------------------------------------------------ tooltip */
const tip = $('tooltip');
function showTip(evt, html) { tip.innerHTML = html; tip.hidden = false; moveTip(evt); }
function moveTip(evt) { const pad = 14; let x = evt.clientX + pad, y = evt.clientY + pad; const r = tip.getBoundingClientRect(); if (x + r.width > innerWidth - 8) x = evt.clientX - r.width - pad; if (y + r.height > innerHeight - 8) y = evt.clientY - r.height - pad; tip.style.left = x + 'px'; tip.style.top = y + 'px'; }
function hideTip() { tip.hidden = true; }
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

/* ------------------------------------------------------------------ navigation */
function showView(name) {
  for (const b of document.querySelectorAll('.tabs [data-view]')) b.setAttribute('aria-selected', String(b.dataset.view === name));
  for (const v of ['explore', 'runs', 'method']) $('view-' + v).hidden = v !== name;
  if (name === 'runs') loadRuns();
  if (name === 'method') renderMethod();
  history.replaceState(null, '', '#' + name);
}
document.querySelectorAll('.tabs [data-view]').forEach((b) => b.addEventListener('click', () => showView(b.dataset.view)));
document.addEventListener('click', (e) => { const g = e.target.closest('[data-goto]'); if (g) { e.preventDefault(); showView(g.dataset.goto); } });

/* ------------------------------------------------------------------ field selection */
const state = { field: null, screening: null, selected: new Set(), targets: [] };

function setField(field) {
  state.field = field;
  $('field-label').textContent = field ? `${field.label} · RA ${field.ra.toFixed(5)}°, Dec ${field.dec.toFixed(5)}°` : 'No field selected.';
  $('screen').disabled = !field;
}
async function loadTargets() {
  try {
    state.targets = await api('/api/lab/targets');
    for (const t of state.targets) $('preset').append(el('option', { value: t.id }, t.name));
  } catch (e) { $('field-label').textContent = 'Could not load presets: ' + e.message; }
}
$('preset').addEventListener('change', () => {
  const t = state.targets.find((x) => x.id === $('preset').value);
  if (!t) return setField(null);
  $('ra').value = ''; $('dec').value = ''; $('name').value = '';
  setField({ target_id: t.id, ra: t.ra_deg, dec: t.dec_deg, label: `${t.name} (preset, proper motion applied)` });
});
function coordsChanged() {
  const ra = parseFloat($('ra').value), dec = parseFloat($('dec').value);
  $('preset').value = '';
  if (Number.isFinite(ra) && Number.isFinite(dec) && ra >= 0 && ra < 360 && dec >= -90 && dec <= 90) setField({ ra, dec, label: 'Fixed coordinates' });
  else setField(null);
}
$('ra').addEventListener('input', coordsChanged); $('dec').addEventListener('input', coordsChanged);
async function resolveName() {
  const name = $('name').value.trim(); if (!name) return;
  $('resolve').disabled = true; $('field-label').textContent = `Resolving “${name}”…`;
  try {
    const r = await api('/api/lab/resolve?name=' + encodeURIComponent(name));
    if (r.target_id) { $('preset').value = r.target_id; $('preset').dispatchEvent(new Event('change')); }
    else { $('preset').value = ''; $('ra').value = r.ra_deg.toFixed(6); $('dec').value = r.dec_deg.toFixed(6); setField({ ra: r.ra_deg, dec: r.dec_deg, label: `${r.name} (${r.source})` }); }
  } catch (e) { $('field-label').textContent = e.message; } finally { $('resolve').disabled = false; }
}
$('resolve').addEventListener('click', resolveName);
$('name').addEventListener('keydown', (e) => { if (e.key === 'Enter') resolveName(); });
function fieldBody() { const f = state.field; return f.target_id ? { target_id: f.target_id } : { ra: f.ra, dec: f.dec }; }

/* ------------------------------------------------------------------ jobs */
const STAGES = {
  screen: [['archive_discovery', 'Query IRSA SPHEREx DR1 archive'], ['calibration_maps', 'Resolve detector wavelength maps (DataLink)'], ['screening', 'Estimate wavelength at target for every record'], ['completed', 'Rank same-wavelength epochs']],
  compare: [['archive_discovery', 'Query archive records'], ['catalog', 'Gaia DR3 field catalog'], ['download', 'Retrieve real FITS cutouts'], ['calibration', 'Own-calibration wavelength maps'], ['comparison', 'Background, PSF match, register, difference, validate'], ['products', 'Write FITS / PNG / report'], ['completed', 'Done']],
  series: [['archive_discovery', 'Query archive records'], ['catalog', 'Gaia DR3 field catalog'], ['download', 'Retrieve real FITS cutouts'], ['calibration', 'Own-calibration wavelength maps'], ['observer_ephemeris', 'JPL Horizons SPHEREx position'], ['measurement', 'ePSF photometry + astrometry'], ['completed', 'Done']],
};
let activeJob = null;
function renderProgress(kind, job) {
  $('job-panel').hidden = false;
  $('job-title').textContent = { screen: 'Screening archive', compare: 'Comparing two epochs', series: 'Building time series' }[kind];
  const stage = job.progress && job.progress.stage;
  const steps = STAGES[kind]; let reached = steps.findIndex((s) => s[0] === stage);
  if (job.status === 'COMPLETED') reached = steps.length;
  const list = $('job-steps'); list.replaceChildren();
  steps.forEach(([key, text], i) => {
    const cls = i < reached ? 'done' : (i === reached && job.status !== 'FAILED' ? 'active' : '');
    let label = text;
    if (i === reached && job.progress && job.progress.observation_id) label += ` · ${short(job.progress.observation_id)}`;
    list.append(el('li', { class: cls }, label));
  });
  const badge = { COMPLETED: ['good', 'Completed'], FAILED: ['critical', 'Failed'], INTERRUPTED: ['critical', 'Interrupted'], RUNNING: ['neutral', 'Running'], QUEUED: ['neutral', 'Queued'] }[job.status] || ['neutral', job.status];
  $('job-badge').replaceChildren(el('span', { class: 'badge ' + badge[0] }, badge[1]));
  $('job-error').replaceChildren(job.error ? el('div', { class: 'notice error' }, el('b', {}, job.error.code + ': '), job.error.message) : '');
}
async function runJob(kind, body) {
  if (activeJob) return;
  setBusy(true);
  try {
    const sub = await api('/api/lab/' + kind, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
    activeJob = sub.id;
    renderProgress(kind, { status: 'QUEUED', progress: { stage: 'queued' } });
    for (;;) {
      await new Promise((r) => setTimeout(r, 1200));
      let job;
      try { job = await api('/api/lab/jobs/' + sub.id); } catch (e) { continue; }
      renderProgress(kind, job);
      if (job.status === 'COMPLETED') return { job, base: `/api/lab/jobs/${sub.id}/files/` };
      if (['FAILED', 'INTERRUPTED'].includes(job.status)) return null;
    }
  } catch (e) {
    $('job-panel').hidden = false; $('job-steps').replaceChildren();
    $('job-error').replaceChildren(el('div', { class: 'notice error' }, e.message));
    return null;
  } finally { activeJob = null; setBusy(false); }
}
function setBusy(busy) { $('screen').disabled = busy || !state.field; updateSelection(busy); }

/* ------------------------------------------------------------------ screening */
$('screen').addEventListener('click', async () => {
  const out = await runJob('screen', fieldBody());
  if (out) { state.selected.clear(); renderScreening(out.job.result, $('screening'), true); }
});

function renderScreening(s, host, interactive) {
  $('welcome').hidden = interactive ? true : $('welcome').hidden;
  host.hidden = false; host.replaceChildren();
  state.screening = s;
  const screened = s.observations.filter((o) => o.status === 'SCREENED');
  const groups3 = s.matched_groups.filter((g) => g.epoch_count >= 3);
  const panel = el('section', { class: 'panel' },
    el('header', {}, el('div', {}, el('h2', {}, 'Archive coverage at this position'),
      el('p', {}, `SPHEREx DR1 records whose footprint contains RA ${fmt(s.query.ra_deg, 5)}°, Dec ${fmt(s.query.dec_deg, 5)}° (search radius ${s.query.radius_deg}°). Wavelengths are screening estimates (footprint model, ≤6 px); every comparison re-checks them on the downloaded cutouts.`))),
    el('div', { class: 'tiles' },
      tile(s.observation_count, 'archive records', `${s.rejected_records.length} outside footprint`),
      tile(s.screened_count, 'wavelength screened', 'detectors D1–D6'),
      tile(s.matched_pairs.length, 'same-wavelength pairs', `|Δλ| ≤ ${s.policy.tolerance_fraction} × bandwidth`),
      tile(groups3.length, 'groups with ≥3 epochs', 'light-curve candidates')));
  const chartHost = el('div', { class: 'chart', style: 'margin-top:16px' });
  panel.append(chartHost);
  panel.append(el('div', { class: 'legend' }, el('span', {}, el('i', { style: 'background:var(--series-1)' }), 'Observation (wavelength sampled at the target)'),
    el('span', {}, el('i', { style: 'background:var(--series-2)' }), 'Member of a same-wavelength pair'),
    el('span', { class: 'muted' }, 'Bands D1–D6 are the six detectors; click points to select')));
  host.append(panel);
  coverageChart(chartHost, screened, s.matched_pairs, interactive);

  const pairs = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Same-wavelength pairs'),
    el('p', {}, 'Same detector, distinct dates, wavelength at the target agreeing within the tolerance. Ranked by the predicted fraction of the 0.3° field that is spectrally matched.'))));
  if (!s.matched_pairs.length) pairs.append(el('div', { class: 'empty' }, 'No same-wavelength pair at this position in the searched records. Different-wavelength epochs are not compared.'));
  else {
    const rows = s.matched_pairs.slice(0, 60).map((p) => el('tr', {},
      el('td', { class: 'mono' }, short(p.observation_a)), el('td', { class: 'mono' }, short(p.observation_b)), el('td', { class: 'num' }, 'D' + p.detector),
      el('td', { class: 'num' }, fmt(p.wavelength_um, 4)), el('td', { class: 'num' }, fmt(p.delta_over_bandwidth, 3)), el('td', { class: 'num' }, fmt(p.separation_days, 1)),
      el('td', { class: 'num' }, p.predicted_matched_fraction === undefined ? '—' : (100 * p.predicted_matched_fraction).toFixed(0) + '%'),
      el('td', {}, interactive ? el('button', { class: 'btn small primary', onclick: () => startCompare(p.observation_a, p.observation_b) }, 'Compare') : '')));
    pairs.append(el('div', { class: 'table-wrap' }, el('table', {}, el('thead', {}, el('tr', {}, ...['Epoch A', 'Epoch B', 'Det.', 'λ (µm)', '|Δλ|/band', 'Δt (days)', 'Matched area', ''].map((h, i) => el('th', { class: i >= 2 && i < 7 ? 'num' : '' }, h)))), el('tbody', {}, rows))));
  }
  host.append(pairs);

  const gp = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Same-wavelength groups'),
    el('p', {}, 'Several epochs on one detector at one wavelength: a light curve and an astrometric trajectory become possible with three or more.'))));
  if (!s.matched_groups.length) gp.append(el('div', { class: 'empty' }, 'No multi-epoch group at a common wavelength.'));
  else gp.append(el('div', { class: 'table-wrap' }, el('table', {}, el('thead', {}, el('tr', {}, ...['Detector', 'λ (µm)', 'Epochs', 'Dates', ''].map((h) => el('th', {}, h)))),
    el('tbody', {}, s.matched_groups.slice(0, 40).map((g) => el('tr', {}, el('td', { class: 'num' }, 'D' + g.detector), el('td', { class: 'num' }, fmt(g.wavelength_um, 4)), el('td', { class: 'num' }, g.epoch_count),
      el('td', { class: 'mono small' }, g.mjd.map(isoDay).join(', ')),
      el('td', {}, interactive ? el('button', { class: 'btn small', onclick: () => startSeries(g.observation_ids.slice(0, 8)) }, 'Time series') : '')))))));
  host.append(gp);

  const obs = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'All observations'),
    el('p', {}, interactive ? 'Tick observations to compare two of them or build a time series from several on one detector.' : 'Screened archive records.'))));
  obs.append(el('div', { class: 'table-wrap' }, el('table', {}, el('thead', {}, el('tr', {}, ...['', 'Date', 'Observation', 'Det.', 'λ at target (µm)', 'Bandwidth (µm)', '± screening (µm)', 'Status'].map((h, i) => el('th', { class: i >= 3 && i <= 6 ? 'num' : '' }, h)))),
    el('tbody', { id: interactive ? 'obs-body' : null }, s.observations.map((o) => {
      const box = interactive && o.status === 'SCREENED' ? el('input', { type: 'checkbox', 'data-id': o.observation_id, 'aria-label': 'select ' + short(o.observation_id) }) : '';
      if (box) box.addEventListener('change', () => toggleSelect(o.observation_id, box.checked));
      return el('tr', { 'data-row': o.observation_id }, el('td', {}, box), el('td', { class: 'mono' }, o.mjd ? isoDay(o.mjd) : '—'), el('td', { class: 'mono' }, short(o.observation_id)),
        el('td', { class: 'num' }, o.detector ? 'D' + o.detector : '—'), el('td', { class: 'num' }, fmt(o.wavelength_um, 4)), el('td', { class: 'num' }, fmt(o.bandwidth_um, 4)),
        el('td', { class: 'num' }, fmt(o.wavelength_uncertainty_um, 4)), el('td', { class: 'small secondary' }, o.status === 'SCREENED' ? 'screened' : (o.reason || o.status)));
    })))));
  host.append(obs);
  if (interactive) { $('selection-panel').hidden = false; updateSelection(); }
}
function tile(value, label, sub) { return el('div', { class: 'tile' }, el('div', { class: 'value' }, value), el('div', { class: 'label' }, label), sub ? el('div', { class: 'sub' }, sub) : ''); }

function toggleSelect(id, on) {
  if (on) state.selected.add(id); else state.selected.delete(id);
  const box = document.querySelector(`#obs-body input[data-id="${CSS.escape(id)}"]`); if (box) box.checked = on;
  document.querySelectorAll('.cov-point').forEach((c) => c.setAttribute('stroke', state.selected.has(c.dataset.id) ? '#ffffff' : 'none'));
  updateSelection();
}
function updateSelection(busy = !!activeJob) {
  const ids = [...state.selected]; const dets = new Set(ids.map((i) => i.slice(-1)));
  $('selection-text').textContent = ids.length ? `${ids.length} selected${dets.size > 1 ? ' — mixed detectors (not comparable)' : ' on D' + [...dets][0]}` : 'Select observations in the table or chart.';
  $('compare-selected').disabled = busy || ids.length !== 2 || dets.size !== 1;
  $('series-selected').disabled = busy || ids.length < 2 || ids.length > 8 || dets.size !== 1;
}
$('compare-selected').addEventListener('click', () => { const [a, b] = [...state.selected]; startCompare(a, b); });
$('series-selected').addEventListener('click', () => startSeries([...state.selected]));

async function startCompare(a, b) {
  const s = state.screening; const ord = (id) => (s.observations.find((o) => o.observation_id === id) || {}).mjd || 0;
  const [first, second] = ord(a) <= ord(b) ? [a, b] : [b, a];
  const out = await runJob('compare', { ...fieldBody(), observation_a: first, observation_b: second });
  if (out) renderComparison(out.job.result, out.base, $('result'));
}
async function startSeries(ids) {
  const out = await runJob('series', { ...fieldBody(), observation_ids: ids });
  if (out) renderSeries(out.job.result, out.base, $('result'));
}

/* ------------------------------------------------------------------ chart helpers */
function chartWidth(host, fallback) { const w = host.getBoundingClientRect().width || (host.parentElement && host.parentElement.getBoundingClientRect().width); return Math.max(560, Math.round(w || fallback)); }
function monthTicks(x0, x1, maxTicks = 10) {
  const start = mjdToDate(x0); const d = new Date(Date.UTC(start.getUTCFullYear(), start.getUTCMonth() + 1, 1)); const months = [];
  while ((d.getTime() / 86400000 + 40587) < x1) { months.push(new Date(d)); d.setUTCMonth(d.getUTCMonth() + 1); }
  const step = Math.max(1, Math.ceil(months.length / maxTicks));
  return months.filter((_, i) => i % step === 0).map((md) => ({ mjd: md.getTime() / 86400000 + 40587, label: md.toISOString().slice(0, 7) }));
}

/* ------------------------------------------------------------------ coverage chart */
function coverageChart(host, obs, pairs, interactive) {
  if (!obs.length) { host.append(el('div', { class: 'empty' }, 'No screened observations.')); return; }
  const W = chartWidth(host, 960), H = 340, m = { l: 64, r: 16, t: 10, b: 38 };
  const mjds = obs.map((o) => o.mjd); let x0 = Math.min(...mjds) - 10, x1 = Math.max(...mjds) + 10;
  const y0 = 0.7, y1 = 5.1;
  const X = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r), Y = (v) => H - m.b - (v - y0) / (y1 - y0) * (H - m.t - m.b);
  const svg = sv('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'Wavelength sampled at the target versus observation date' });
  const bands = [[0.75, 1.09, 'D1'], [1.10, 1.62, 'D2'], [1.63, 2.41, 'D3'], [2.42, 3.82, 'D4'], [3.83, 4.41, 'D5'], [4.42, 5.0, 'D6']];
  bands.forEach(([lo, hi, label], i) => {
    if (i % 2 === 0) svg.append(sv('rect', { x: m.l, y: Y(hi), width: W - m.l - m.r, height: Y(lo) - Y(hi), fill: '#ffffff', 'fill-opacity': 0.03 }));
    svg.append(sv('text', { x: W - m.r - 4, y: Y((lo + hi) / 2) + 4, 'text-anchor': 'end', class: 'band-label' }, label));
  });
  const grid = sv('g', { class: 'grid' }), axis = sv('g', { class: 'axis' });
  for (let v = 1; v <= 5; v += 1) { grid.append(sv('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v) })); axis.append(sv('text', { x: m.l - 8, y: Y(v) + 4, 'text-anchor': 'end' }, v.toFixed(1))); }
  monthTicks(x0, x1, Math.floor(W / 95)).forEach((t) => { axis.append(sv('text', { x: X(t.mjd), y: H - m.b + 18, 'text-anchor': 'middle' }, t.label)); grid.append(sv('line', { x1: X(t.mjd), x2: X(t.mjd), y1: m.t, y2: H - m.b })); });
  axis.append(sv('text', { x: 14, y: (H - m.b) / 2, transform: `rotate(-90 14 ${(H - m.b) / 2})`, 'text-anchor': 'middle', class: 'axis-title' }, 'Wavelength at target (µm)'));
  svg.append(grid, axis);
  const byId = new Map(obs.map((o) => [o.observation_id, o]));
  const paired = new Set(); pairs.forEach((p) => { paired.add(p.observation_a); paired.add(p.observation_b); });
  const links = sv('g', {});
  pairs.slice(0, 200).forEach((p) => { const a = byId.get(p.observation_a), b = byId.get(p.observation_b); if (a && b) links.append(sv('line', { x1: X(a.mjd), y1: Y(a.wavelength_um), x2: X(b.mjd), y2: Y(b.wavelength_um), stroke: 'var(--series-2)', 'stroke-width': 1, 'stroke-opacity': 0.35 })); });
  svg.append(links);
  const pts = sv('g', {});
  obs.forEach((o) => {
    const c = sv('circle', { cx: X(o.mjd), cy: Y(o.wavelength_um), r: 4.5, fill: paired.has(o.observation_id) ? 'var(--series-2)' : 'var(--series-1)', stroke: state.selected.has(o.observation_id) ? '#ffffff' : 'none', 'stroke-width': 2, class: 'cov-point', 'data-id': o.observation_id, style: interactive ? 'cursor:pointer' : '' });
    c.addEventListener('mouseenter', (e) => showTip(e, `<b>${esc(short(o.observation_id))}</b><br>${isoDay(o.mjd)} · D${o.detector}<br>λ ${fmt(o.wavelength_um, 4)} µm ± ${fmt(o.wavelength_uncertainty_um, 4)} (screening)<br><span class="muted">bandwidth ${fmt(o.bandwidth_um, 4)} µm${paired.has(o.observation_id) ? ' · in a same-wavelength pair' : ''}</span>`));
    c.addEventListener('mousemove', moveTip); c.addEventListener('mouseleave', hideTip);
    if (interactive) c.addEventListener('click', () => toggleSelect(o.observation_id, !state.selected.has(o.observation_id)));
    pts.append(c);
  });
  svg.append(pts);
  host.append(svg);
}

/* ------------------------------------------------------------------ comparison result */
const STATUS_BADGE = { VALIDATED: 'good', PARTIAL: 'warning', REJECTED: 'critical', PASS: 'good', FAIL: 'critical', UNVERIFIED: 'warning' };
const GATE_TEXT = { spectral_overlap: 'Spectral overlap', photometric_consistency: 'Photometric scale', noise_model: 'Pixel noise model', source_noise_model: 'Source noise model', registration: 'Registration' };
const CLASS_TEXT = {
  CANDIDATE_VARIABLE: 'Candidate variable (catalog source)', UNCATALOGED_VARIABLE_CANDIDATE: 'Candidate variable (not in Gaia)', TRANSIENT_CANDIDATE: 'Transient candidate (one epoch only)',
  KNOWN_PROPER_MOTION_RESIDUAL: 'Known proper motion', TARGET_RESIDUAL: 'Residual at target', DIPOLE_PSF_OR_REGISTRATION: 'Dipole: PSF/registration', CONSISTENT_WITH_SYSTEMATIC_FLOOR: 'Within systematic floor',
  PERSISTENT_SOURCE_RESIDUAL: 'Persistent source, no significant change', UNASSOCIATED_RESIDUAL: 'Low-significance residual', CATALOG_SOURCE_UNMEASURED: 'Catalog source, unmeasurable',
  EDGE_ARTIFACT_LIKELY: 'Edge artifact', UNCLASSIFIED_NO_CATALOG: 'Unclassified (no catalog)'
};
const INTERESTING = new Set(['CANDIDATE_VARIABLE', 'UNCATALOGED_VARIABLE_CANDIDATE', 'TRANSIENT_CANDIDATE']);

async function renderComparison(report, base, host) {
  host.hidden = false; host.replaceChildren();
  const r = report; const ea = r.epochs.a, eb = r.epochs.b;
  const head = el('section', { class: 'panel' },
    el('header', {}, el('div', {}, el('h2', {}, `${short(r.observation_a)}  →  ${short(r.observation_b)}`),
      el('p', {}, `${ea.date.slice(0, 10)} → ${eb.date.slice(0, 10)} · ${fmt(r.separation_days, 1)} days · detector D${ea.detector} · matched λ ${r.spectral.wavelength_range_um ? fmt(r.spectral.wavelength_range_um[0], 4) + '–' + fmt(r.spectral.wavelength_range_um[1], 4) + ' µm' : 'none'}`)),
      el('span', { class: 'badge ' + (STATUS_BADGE[r.status] || 'neutral') }, r.status)),
    el('p', { class: 'secondary small', style: 'margin:0 0 12px' }, r.reason || ''),
    el('div', { style: 'display:flex;gap:8px;flex-wrap:wrap' }, Object.entries(r.gates || {}).map(([k, v]) => el('span', { class: 'badge ' + (STATUS_BADGE[v] || 'neutral'), title: k }, `${GATE_TEXT[k] || k}: ${v}`))));
  host.append(head);
  const viewerPanel = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Images'), el('p', {}, 'Background-subtracted, PSF-matched images on epoch A’s pixel grid. Pixels outside the spectrally matched region are dimmed: the two epochs sampled different wavelengths there.'))));
  const viewerHost = el('div', {}, el('div', { class: 'empty' }, 'Loading image data…'));
  viewerPanel.append(viewerHost); host.append(viewerPanel);
  host.append(validationPanel(r));
  host.append(targetPanel(r));
  const candHost = el('section', { class: 'panel' }); host.append(candHost);
  host.append(variabilityPanel(r));
  host.append(el('section', { class: 'panel' }, el('header', {}, el('h2', {}, 'Limitations')), el('ul', { class: 'plain' }, (r.limitations || []).map((l) => el('li', {}, l)))));
  host.append(downloadsPanel(base, [['comparison.fits', 'FITS products'], ['report.json', 'Full report (JSON)'], ['candidates.csv', 'Residual table (CSV)'], ['comparison.png', 'Summary figure (PNG)']], r));
  let viewer = null;
  try {
    const data = await api(base + 'viewer.json');
    viewerHost.replaceChildren();
    viewer = new ImageViewer(viewerHost, data, r);
  } catch (e) { viewerHost.replaceChildren(el('div', { class: 'notice error' }, 'Image data unavailable: ' + e.message)); }
  candidatesPanel(candHost, r, viewer);
  host.scrollIntoView({ behavior: 'smooth', block: 'start' });
}

function check(title, status, metric, text) {
  return el('div', { class: 'check' }, el('div', { class: 'head' }, el('h4', {}, title), status ? el('span', { class: 'badge ' + (STATUS_BADGE[status] || 'neutral') }, status) : ''), el('div', { class: 'metric' }, metric), el('p', {}, text));
}
function validationPanel(r) {
  const g = r.gates || {}; const sp = r.spectral || {}; const reg = r.registration || {}; const ph = r.photometry || {}; const nt = r.null_test || {}; const sf = r.systematic_floor || {}; const inj = r.injection_recovery || {};
  const cards = el('div', { class: 'check-grid' },
    check('Spectral match', g.spectral_overlap, `${sp.matched_pixels} px · ${(100 * (sp.matched_fraction_of_overlap || 0)).toFixed(1)}%`, `of overlapping pixels have |Δλ| ≤ ${r.policy.spectral_tolerance_fraction} × CBAND (median ${fmt(sp.median_delta_over_band_matched, 3)}). Elsewhere median ${fmt(sp.median_delta_over_band_overlap, 2)}.`),
    check('Registration', g.registration, reg.median_offset_arcsec ? `${fmt(reg.median_offset_arcsec[0], 3)}″, ${fmt(reg.median_offset_arcsec[1], 3)}″` : '—', reg.stars ? `Median calibration-source centroid offset after resampling (${reg.stars} stars; scatter ${fmt(reg.robust_scatter_px && reg.robust_scatter_px[0], 3)} px).` : (reg.reason || '')),
    check('Photometric scale', g.photometric_consistency, `B/A = ${fmt(ph.median_ratio_b_over_a, 4)} ± ${fmt(ph.median_ratio_standard_error, 4)}`, `${ph.scale_stars || 0} stars (${ph.scale_star_selection || '—'}). ${ph.decision || ''}`),
    check('Blank-sky null test', g.noise_model, nt.robust_std_sigma !== undefined ? `σ = ${fmt(nt.robust_std_sigma, 2)} · median ${fmt(nt.median_sigma, 2)}` : '—', nt.robust_std_sigma !== undefined ? `${nt.pixels} source-free pixels; uncertainties multiplied by ${fmt(nt.noise_scale_applied, 2)} (inflate only).` : (nt.reason || '')),
    check('Source noise model', g.source_noise_model, sf.alpha_fractional !== undefined && sf.alpha_fractional !== null ? `α ${(100 * sf.alpha_fractional).toFixed(1)}% · β ${fmt(sf.beta_aperture_mjy_sr_px, 3)}` : '—', `${sf.stars} calibration sources; calibrated z-scatter ${fmt(sf.calibration_check && sf.calibration_check.all && sf.calibration_check.all.robust_std_z, 2)} (1 = calibrated). α multiplicative, β additive (MJy/sr·px).`),
    check('Injection recovery', null, inj.completeness !== undefined && inj.completeness !== null ? `${inj.recovered}/${inj.injected} = ${(100 * inj.completeness).toFixed(0)}%` : '—', inj.status === 'COMPLETED' ? `Synthetic PSF sources at aperture S/N ${inj.aperture_snr} added to the real difference image; detection rerun with the ${r.policy.detection_sigma}σ production threshold.` : (inj.reason || '')));
  const bg = r.difference_background || {};
  return el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Validation'), el('p', {}, 'Each gate is computed from this pair’s own data. UNVERIFIED means the data could not test it, not that it passed.'))), cards,
    el('p', { class: 'small muted', style: 'margin:12px 0 0' }, `PSF: ${r.psf.method}; zones ${r.psf.zone_a}/${r.psf.zone_b}; interpolation smoothing ${r.psf.interpolation_smoothing_variance_px2 ? r.psf.interpolation_smoothing_variance_px2.map((v) => fmt(v, 3)).join('/') : '—'} px² matched on A. Difference pedestal ${bg.applied ? fmt(bg.median_mjy_sr, 4) + ' MJy/sr removed' : 'not fitted'}. Catalog: ${(r.catalog && r.catalog.mode) || 'none'}. Calibration sources: ${r.calibration_sources || '—'}.`));
}
function targetPanel(r) {
  const t = r.target; const p = el('section', { class: 'panel' }, el('header', {}, el('h2', {}, 'Target')));
  if (!t) { p.append(el('p', { class: 'secondary' }, 'Fixed-coordinate field: no target measurement requested.')); return p; }
  if (t.status === 'MEASURED') p.append(el('dl', { class: 'kv' }, el('dt', {}, 'Fractional flux change (B/A − 1)'), el('dd', {}, `${fmt(t.fractional_change, 4)} ± ${fmt(t.fractional_change_sigma, 4)}`), el('dt', {}, 'Significance'), el('dd', {}, `${fmt(t.change_z, 2)} σ`), el('dt', {}, 'S/N in A'), el('dd', {}, fmt(t.snr_a, 1))));
  else p.append(el('div', { class: 'notice' }, el('b', {}, t.status.replaceAll('_', ' ') + ': '), t.reason || ''));
  return p;
}
function candidatesPanel(host, r, viewer) {
  const cands = r.candidates || [];
  host.append(el('header', {}, el('div', {}, el('h2', {}, 'Residuals in the difference image'), el('p', {}, `Connected ≥${r.policy.detection_min_pixels}-pixel components at ≥${r.policy.detection_sigma}σ, classified against the propagated Gaia catalog and the calibrated noise model. Click a row to centre it.`))));
  const summary = Object.entries(r.candidate_summary || {}).sort((a, b) => b[1] - a[1]);
  host.append(el('div', { style: 'display:flex;gap:8px;flex-wrap:wrap;margin-bottom:12px' }, summary.length ? summary.map(([k, n]) => el('span', { class: 'badge ' + (INTERESTING.has(k) ? 'warning' : 'neutral') }, `${CLASS_TEXT[k] || k}: ${n}`)) : el('span', { class: 'muted' }, 'No residual components above threshold.')));
  if (!cands.length) return;
  const filter = el('select', { style: 'width:auto;margin-bottom:10px' }, el('option', { value: 'interesting' }, 'Candidates needing follow-up'), el('option', { value: 'all' }, 'All residuals'));
  const wrap = el('div', { class: 'table-wrap' });
  const draw = () => {
    const rows = cands.filter((c) => filter.value === 'all' || INTERESTING.has(c.classification));
    wrap.replaceChildren(rows.length ? el('table', {}, el('thead', {}, el('tr', {}, ...['#', 'Classification', 'Peak σ', 'RA', 'Dec', 'Gaia G', 'ΔF/F', 'Reason'].map((h, i) => el('th', { class: [2, 3, 4, 5, 6].includes(i) ? 'num' : '' }, h)))),
      el('tbody', {}, rows.map((c) => el('tr', { class: 'clickable', onclick: () => viewer && viewer.focus(c.x, c.y) },
        el('td', { class: 'num' }, c.id), el('td', {}, CLASS_TEXT[c.classification] || c.classification), el('td', { class: 'num' }, fmt(c.peak_sigma, 1)),
        el('td', { class: 'num' }, fmt(c.ra_deg, 5)), el('td', { class: 'num' }, fmt(c.dec_deg, 5)), el('td', { class: 'num' }, fmt(c.catalog_g_mag, 2)),
        el('td', { class: 'num' }, c.fractional_change === undefined ? '—' : `${fmt(c.fractional_change, 3)} ± ${fmt(c.fractional_change_sigma, 3)}`), el('td', { class: 'small secondary' }, c.reason || '')))))
      : el('div', { class: 'empty' }, 'No residual needs follow-up: every component is explained by known motion, PSF/registration dipoles, the calibrated systematic floor or a persistent unchanged source.'));
  };
  filter.addEventListener('change', draw); draw();
  host.append(filter, wrap);
}
function variabilityPanel(r) {
  const v = r.variability_candidates || []; const s = r.variability_summary || {};
  const p = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Source-level flux comparison'),
    el('p', {}, `${s.stars_tested || 0} Gaia stars measured in matched 5×5 apertures in both epochs; ${s.beyond_threshold || 0} exceed ${r.policy.detection_sigma}σ under the calibrated noise model; ${s.passed_vetting || 0} pass vetting (isolation from flagged pixels, distance from region edge, no catalog blend, PSF-shaped residual).`))));
  if (!v.length) { p.append(el('div', { class: 'empty' }, 'No star changed beyond the calibrated noise model.')); return p; }
  p.append(el('div', { class: 'table-wrap' }, el('table', {}, el('thead', {}, el('tr', {}, ...['Gaia DR3 source', 'G', 'ΔF/F', 'z', 'Verdict', 'Failed checks'].map((h, i) => el('th', { class: [1, 2, 3].includes(i) ? 'num' : '' }, h)))),
    el('tbody', {}, v.map((c) => el('tr', {}, el('td', { class: 'mono' }, c.source_id), el('td', { class: 'num' }, fmt(c.g_mag, 2)), el('td', { class: 'num' }, `${fmt(c.fractional_change, 3)} ± ${fmt(c.fractional_change_sigma, 3)}`),
      el('td', { class: 'num' }, fmt(c.change_z, 1)), el('td', {}, el('span', { class: 'badge ' + (c.classification === 'SOURCE_LEVEL_VARIABILITY_CANDIDATE' ? 'warning' : 'neutral') }, c.classification === 'SOURCE_LEVEL_VARIABILITY_CANDIDATE' ? 'Candidate' : 'Rejected')),
      el('td', { class: 'small secondary' }, (c.vetting && c.vetting.failed.length) ? c.vetting.failed.join(', ').replaceAll('_', ' ') : '—')))))));
  return p;
}
function downloadsPanel(base, files, report) {
  const p = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Downloads & provenance'))));
  p.append(el('div', { class: 'downloads' }, files.map(([f, label]) => el('a', { class: 'btn small', href: base + f, download: f }, label))));
  if (report && report.software_sha256) p.append(el('p', { class: 'small muted', style: 'margin:12px 0 0' }, `Software SHA-256 ${report.software_sha256.slice(0, 16)}… · science files: ${Object.entries(report.science_files || {}).map(([k, v]) => short(k) + ' (' + (v.sha256 || '').slice(0, 10) + '…)').join(', ') || 'see report'}`));
  return p;
}

/* ------------------------------------------------------------------ image viewer */
function decode(b64) { const bin = atob(b64); const bytes = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i); return new Float32Array(bytes.buffer); }
function percentile(values, q) { const v = values.filter(Number.isFinite).sort((a, b) => a - b); if (!v.length) return 0; return v[Math.min(v.length - 1, Math.max(0, Math.floor(q * (v.length - 1))))]; }
function lerp(a, b, t) { return a + (b - a) * t; }
const hex = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
const DIV = { neg: [hex('#cde2fb'), hex('#3987e5')], mid: hex('#383835'), pos: [hex('#e66767'), hex('#fbd5d5')] };
function divergingColor(t) { // t in [-1, 1]
  t = Math.max(-1, Math.min(1, t));
  if (t < 0) { const u = -t; return u < 0.6 ? DIV.mid.map((c, i) => lerp(c, DIV.neg[1][i], u / 0.6)) : DIV.neg[1].map((c, i) => lerp(c, DIV.neg[0][i], (u - 0.6) / 0.4)); }
  return t < 0.6 ? DIV.mid.map((c, i) => lerp(c, DIV.pos[0][i], t / 0.6)) : DIV.pos[0].map((c, i) => lerp(c, DIV.pos[1][i], (t - 0.6) / 0.4));
}

class ImageViewer {
  constructor(host, data, report) {
    this.report = report; this.h = data.shape[0]; this.w = data.shape[1];
    this.planes = Object.fromEntries(Object.entries(data.planes).map(([k, v]) => [k, decode(v)]));
    this.matched = decode(data.masks.matched_mask); this.ra = decode(data.ra); this.dec = decode(data.dec); this.wave = decode(data.wave_um);
    const both = [...this.planes.image_a, ...this.planes.image_b_aligned].filter(Number.isFinite);
    this.lo = percentile(both, 0.02); this.hi = percentile(both, 0.995);
    this.dlim = percentile([...this.planes.difference].filter(Number.isFinite).map(Math.abs), 0.99) || 1;
    this.mode = 'side'; this.alpha = 0.5; this.blinkOn = false; this.showMarks = true;
    this.zoom = 1; this.cx = this.w / 2; this.cy = this.h / 2;
    this.offscreen = {};
    for (const key of ['image_a', 'image_b_aligned', 'difference', 'significance', 'delta_over_band']) this.offscreen[key] = this.paint(key);
    this.build(host); this.render();
  }
  value(key, i) { return this.planes[key][i]; }
  colorFor(key, v) {
    if (!Number.isFinite(v)) return null;
    if (key === 'image_a' || key === 'image_b_aligned') { const s = Math.asinh((v - this.lo) / ((this.hi - this.lo) / 12 || 1)) / Math.asinh(12); const g = 255 * Math.max(0, Math.min(1, s)); return [g, g, g]; }
    if (key === 'difference') return divergingColor(v / this.dlim);
    if (key === 'significance') return divergingColor(v / 6);
    if (key === 'delta_over_band') { const t = Math.min(1, v / 1.0); const ok = v <= this.report.policy.spectral_tolerance_fraction; return ok ? [57, 135, 229] : [lerp(70, 217, t), lerp(70, 89, t), lerp(70, 38, t)]; }
    return [0, 0, 0];
  }
  paint(key) {
    const c = document.createElement('canvas'); c.width = this.w; c.height = this.h; const ctx = c.getContext('2d'); const img = ctx.createImageData(this.w, this.h);
    for (let y = 0; y < this.h; y++) for (let x = 0; x < this.w; x++) {
      const src = y * this.w + x, dst = ((this.h - 1 - y) * this.w + x) * 4; // FITS row 0 at the bottom
      const col = this.colorFor(key, this.value(key, src));
      if (!col) { img.data[dst + 3] = 255; img.data[dst] = 18; img.data[dst + 1] = 18; img.data[dst + 2] = 20; continue; }
      const dim = key !== 'delta_over_band' && !this.matched[src] ? 0.35 : 1;
      img.data[dst] = col[0] * dim; img.data[dst + 1] = col[1] * dim; img.data[dst + 2] = col[2] * dim; img.data[dst + 3] = 255;
    }
    ctx.putImageData(img, 0, 0); return c;
  }
  build(host) {
    const modes = [['side', 'Side by side'], ['overlay', 'Overlay'], ['blink', 'Blink'], ['difference', 'Difference'], ['significance', 'Significance'], ['spectral', 'Spectral match']];
    this.seg = el('div', { class: 'seg', role: 'group', 'aria-label': 'View mode' }, modes.map(([m, label]) => el('button', { 'aria-pressed': String(m === this.mode), onclick: () => this.setMode(m) }, label)));
    this.slider = el('input', { type: 'range', min: 0, max: 100, value: 50, 'aria-label': 'Overlay opacity of epoch B', style: 'width:140px' });
    this.slider.addEventListener('input', () => { this.alpha = this.slider.value / 100; this.render(); });
    this.sliderWrap = el('label', { class: 'small secondary', style: 'display:none;align-items:center;gap:6px' }, 'A', this.slider, 'B');
    const marks = el('label', { class: 'small secondary', style: 'display:flex;align-items:center;gap:6px;margin-left:auto' }, el('input', { type: 'checkbox', checked: true, onchange: (e) => { this.showMarks = e.target.checked; this.render(); } }), 'Residual markers');
    const reset = el('button', { class: 'btn small', onclick: () => { this.zoom = 1; this.cx = this.w / 2; this.cy = this.h / 2; this.render(); } }, 'Reset view');
    this.stage = el('div', { class: 'stage split' });
    this.canvas = el('canvas', {}); this.stage.append(this.canvas);
    this.capL = el('div', { class: 'caption' }); this.capR = el('div', { class: 'caption right' }); this.stage.append(this.capL, this.capR);
    this.colorbar = el('div', {}); this.readout = el('dl', { class: 'readout' });
    const side = el('div', { class: 'stack' }, el('h4', {}, 'Pixel readout'), this.readout, this.colorbar, el('p', { class: 'small muted' }, 'Scroll to zoom, drag to pan. Grey-scale images share one asinh stretch so brightness is comparable between epochs.'));
    host.append(el('div', { class: 'viewer-toolbar' }, this.seg, this.sliderWrap, marks, reset), el('div', { class: 'viewer' }, this.stage, side));
    this.readout.append(el('dt', {}, 'Hover'), el('dd', {}, 'over the image'));
    const ro = new ResizeObserver(() => this.render()); ro.observe(this.stage);
    this.stage.addEventListener('wheel', (e) => { e.preventDefault(); const f = e.deltaY < 0 ? 1.25 : 0.8; this.zoom = Math.max(1, Math.min(16, this.zoom * f)); this.render(); }, { passive: false });
    let drag = null;
    this.stage.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY, cx: this.cx, cy: this.cy }; this.stage.setPointerCapture(e.pointerId); });
    this.stage.addEventListener('pointerup', () => { drag = null; });
    this.stage.addEventListener('pointermove', (e) => {
      if (drag) { const s = this.scale(); this.cx = drag.cx - (e.clientX - drag.x) / s; this.cy = drag.cy + (e.clientY - drag.y) / s; this.render(); }
      this.hover(e);
    });
    this.stage.addEventListener('pointerleave', () => { this.render(); });
  }
  setMode(m) {
    this.mode = m; [...this.seg.children].forEach((b, i) => b.setAttribute('aria-pressed', String(['side', 'overlay', 'blink', 'difference', 'significance', 'spectral'][i] === m)));
    this.stage.classList.toggle('split', m === 'side'); this.sliderWrap.style.display = m === 'overlay' ? 'inline-flex' : 'none';
    clearInterval(this.timer); if (m === 'blink') this.timer = setInterval(() => { this.blinkOn = !this.blinkOn; this.render(); }, 700);
    this.render();
  }
  panels() { return this.mode === 'side' ? 2 : 1; }
  scale() { const r = this.stage.getBoundingClientRect(); return Math.min(r.width / this.panels() / this.w, r.height / this.h) * this.zoom; }
  layers() {
    const ea = this.report.epochs.a.date.slice(0, 10), eb = this.report.epochs.b.date.slice(0, 10);
    switch (this.mode) {
      case 'side': return [[['image_a', 1]], [['image_b_aligned', 1]], `A · ${ea}`, `B · ${eb}`];
      case 'overlay': return [[['image_a', 1], ['image_b_aligned', this.alpha]], null, `A ${(100 - this.alpha * 100).toFixed(0)}% / B ${(this.alpha * 100).toFixed(0)}%`, ''];
      case 'blink': return [[[this.blinkOn ? 'image_b_aligned' : 'image_a', 1]], null, this.blinkOn ? `B · ${eb}` : `A · ${ea}`, ''];
      case 'difference': return [[['difference', 1]], null, 'B − sA (MJy/sr)', ''];
      case 'significance': return [[['significance', 1]], null, 'Significance (σ, calibrated)', ''];
      default: return [[['delta_over_band', 1]], null, '|Δλ| / bandwidth', ''];
    }
  }
  render() {
    const rect = this.stage.getBoundingClientRect(); if (!rect.width) return;
    const dpr = window.devicePixelRatio || 1; this.canvas.width = Math.round(rect.width * dpr); this.canvas.height = Math.round(rect.height * dpr);
    const ctx = this.canvas.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.imageSmoothingEnabled = false; ctx.fillStyle = '#000'; ctx.fillRect(0, 0, rect.width, rect.height);
    const [left, right, capL, capR] = this.layers(); const pw = rect.width / this.panels(); const s = this.scale();
    const drawPanel = (layers, ox) => {
      ctx.save(); ctx.beginPath(); ctx.rect(ox, 0, pw, rect.height); ctx.clip();
      const x0 = ox + pw / 2 - this.cx * s, y0 = rect.height / 2 - (this.h - this.cy) * s;
      for (const [key, a] of layers) { ctx.globalAlpha = a; ctx.drawImage(this.offscreen[key], x0, y0, this.w * s, this.h * s); }
      ctx.globalAlpha = 1;
      const residualView = this.mode === 'difference' || this.mode === 'significance';
      if (this.showMarks) for (const c of this.report.candidates || []) {
        if (!residualView && !INTERESTING.has(c.classification)) continue;
        const px = x0 + (c.x + 0.5) * s, py = y0 + (this.h - c.y - 0.5) * s;
        ctx.beginPath(); ctx.arc(px, py, Math.max(6, 2.5 * s), 0, 2 * Math.PI); ctx.lineWidth = INTERESTING.has(c.classification) ? 2 : 1;
        ctx.strokeStyle = INTERESTING.has(c.classification) ? '#fab219' : 'rgba(195,194,183,0.55)'; ctx.stroke();
      }
      ctx.restore();
    };
    drawPanel(left, 0); if (right) { drawPanel(right, pw); ctx.fillStyle = '#34352f'; ctx.fillRect(pw - 1, 0, 2, rect.height); }
    this.capL.textContent = capL; this.capR.textContent = capR; this.capR.hidden = !capR;
    this.drawColorbar();
  }
  drawColorbar() {
    const key = this.layers()[0][0][0]; let grad, lo, hi;
    if (key === 'difference' || key === 'significance') { const stops = []; for (let i = 0; i <= 10; i++) { const c = divergingColor(-1 + i / 5); stops.push(`rgb(${c.map(Math.round).join(',')}) ${i * 10}%`); } grad = `linear-gradient(90deg, ${stops.join(',')})`; const L = key === 'difference' ? this.dlim : 6; lo = fmt(-L, 2); hi = fmt(L, 2); }
    else if (key === 'delta_over_band') { grad = `linear-gradient(90deg, rgb(57,135,229) 0%, rgb(57,135,229) 10%, rgb(70,70,70) 10%, rgb(217,89,38) 100%)`; lo = '0'; hi = '≥1 (different colour)'; }
    else { grad = 'linear-gradient(90deg, #000, #fff)'; lo = fmt(this.lo, 2); hi = fmt(this.hi, 2) + ' MJy/sr'; }
    this.colorbar.replaceChildren(el('div', { class: 'colorbar', style: `background:${grad}` }), el('div', { class: 'colorbar-labels' }, el('span', {}, lo), el('span', {}, hi)));
  }
  hover(e) {
    const rect = this.stage.getBoundingClientRect(); const pw = rect.width / this.panels(); const s = this.scale();
    const ox = (e.clientX - rect.left) >= pw && this.panels() === 2 ? pw : 0;
    const x = Math.floor(this.cx + (e.clientX - rect.left - ox - pw / 2) / s), y = Math.floor(this.cy - (e.clientY - rect.top - rect.height / 2) / s);
    if (x < 0 || y < 0 || x >= this.w || y >= this.h) return;
    const i = y * this.w + x; const rows = [['x, y', `${x}, ${y}`], ['RA', fmt(this.ra[i], 5) + '°'], ['Dec', fmt(this.dec[i], 5) + '°'], ['λ (A)', fmt(this.wave[i], 4) + ' µm'],
      ['|Δλ|/band', fmt(this.planes.delta_over_band[i], 3)], ['Matched', this.matched[i] ? 'yes' : 'no'], ['A', fmt(this.planes.image_a[i], 3)], ['B', fmt(this.planes.image_b_aligned[i], 3)], ['B − sA', fmt(this.planes.difference[i], 3)], ['σ', fmt(this.planes.significance[i], 2)]];
    this.readout.replaceChildren(...rows.flatMap(([k, v]) => [el('dt', {}, k), el('dd', {}, v)]));
  }
  focus(x, y) { this.cx = x + 0.5; this.cy = y + 0.5; this.zoom = 6; if (this.mode === 'side') this.setMode('significance'); else this.render(); this.stage.scrollIntoView({ behavior: 'smooth', block: 'center' }); }
}

/* ------------------------------------------------------------------ series result */
function renderSeries(r, base, host) {
  host.hidden = false; host.replaceChildren();
  const ph = r.photometry, as = r.astrometry; const pts = ph.points || [];
  const phBadge = { CONSISTENT_WITH_CONSTANT: 'good', VARIABILITY_CANDIDATE: 'warning' }[ph.status] || 'neutral';
  const asBadge = { MOTION_MEASURED: 'good', NO_SIGNIFICANT_MOTION: 'neutral' }[as.status] || 'neutral';
  host.append(el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, `Time series · ${r.epoch_count} epochs on D${r.detector}`),
    el('p', {}, 'Same-wavelength photometry and astrometry of the target from measured-ePSF fits on native pixels. Photometry and motion are analysed separately.')),
    el('div', { style: 'display:flex;gap:8px' }, el('span', { class: 'badge ' + phBadge }, 'Light curve: ' + ph.status.replaceAll('_', ' ')), el('span', { class: 'badge ' + asBadge }, 'Motion: ' + as.status.replaceAll('_', ' '))))));
  const lc = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Same-wavelength light curve'), el('p', {}, ph.interpretation || ph.reason || ''))));
  const used = pts.filter((p) => p.status === 'USED');
  host.append(lc);
  if (used.length) { const c = el('div', { class: 'chart' }); lc.append(c); lightCurve(c, pts, ph); }
  lc.append(el('div', { class: 'table-wrap', style: 'margin-top:12px' }, el('table', {}, el('thead', {}, el('tr', {}, ...['Date', 'Observation', 'Status', 'λ (µm)', 'Flux (mJy)', 'σ (mJy)', 'Zero point', 'Note'].map((h, i) => el('th', { class: [3, 4, 5, 6].includes(i) ? 'num' : '' }, h)))),
    el('tbody', {}, pts.map((p) => el('tr', {}, el('td', { class: 'mono' }, p.date.slice(0, 10)), el('td', { class: 'mono' }, short(p.observation_id)), el('td', {}, el('span', { class: 'badge ' + (p.status === 'USED' ? 'good' : 'neutral') }, p.status)),
      el('td', { class: 'num' }, fmt(p.wavelength_um, 4)), el('td', { class: 'num' }, p.calibrated_flux_jy ? fmt(p.calibrated_flux_jy * 1e3, 3) : '—'), el('td', { class: 'num' }, p.flux_sigma_jy ? fmt(p.flux_sigma_jy * 1e3, 3) : '—'),
      el('td', { class: 'num' }, fmt(p.zero_point, 4)), el('td', { class: 'small secondary' }, p.reason || '')))))));
  if (ph.constant_flux_test) lc.append(el('p', { class: 'small muted', style: 'margin:10px 0 0' }, `Constant-flux test: χ² = ${fmt(ph.constant_flux_test.chi2_constant, 2)} for ${ph.constant_flux_test.dof} dof (p = ${fmt(ph.constant_flux_test.pvalue, 3)}); weighted mean ${fmt(ph.constant_flux_test.constant_flux * 1e3, 3)} mJy. Field-star floor ${fmt(ph.fractional_floor, 3)}.`));
  const ap = el('section', { class: 'panel' }, el('header', {}, el('div', {}, el('h2', {}, 'Astrometric trajectory'), el('p', {}, as.reason || 'Target positions corrected by the per-epoch frame offset of image-detected field stars; error bars include the empirical field-star floor.'))));
  host.append(ap);
  if ((as.positions || []).length) { const c = el('div', { class: 'chart' }); ap.append(c); trajectory(c, as); }
  const sol = as.astrometric_solution;
  if (sol && sol.status === 'SOLVED') ap.append(el('div', { class: 'check-grid', style: 'margin-top:12px' },
    check('Proper motion (measured)', null, `${fmt(sol.pm_east_arcsec_per_year, 2)}, ${fmt(sol.pm_north_arcsec_per_year, 2)} ″/yr`, `± ${fmt(sol.pm_sigma_arcsec_per_year[0], 2)}, ${fmt(sol.pm_sigma_arcsec_per_year[1], 2)} (Birge ratio ${fmt(sol.birge_ratio_applied, 2)} applied); detected at ${fmt(sol.pm_significance, 1)}σ.`),
    check('Catalog proper motion', null, sol.catalog_pm_arcsec_per_year ? `${fmt(sol.catalog_pm_arcsec_per_year[0], 2)}, ${fmt(sol.catalog_pm_arcsec_per_year[1], 2)} ″/yr` : '—', `Difference ${fmt(sol.pm_difference_from_catalog_sigma, 2)}σ. SIMBAD values; a known-mover test, not a discovery.`),
    check('Parallax', null, `${fmt(sol.parallax_arcsec, 3)} ± ${fmt(sol.parallax_sigma_arcsec, 3)}″`, `Catalog ${fmt(sol.catalog_parallax_arcsec, 3)}″ (${fmt(sol.parallax_difference_from_catalog_sigma, 1)}σ apart). Observer: JPL Horizons SPHEREx vectors.`),
    check('Astrometric floor', null, `${fmt(as.astrometric_floor_arcsec, 3)}″`, `Epoch-to-epoch repeatability of ${as.field_stars} field stars beyond fit errors.`)));
  host.append(el('section', { class: 'panel' }, el('header', {}, el('h2', {}, 'Limitations')), el('ul', { class: 'plain' }, (r.limitations || []).map((l) => el('li', {}, l)))));
  host.append(downloadsPanel(base, [['series.json', 'Full report (JSON)'], ['light_curve.csv', 'Light curve (CSV)'], ['series.png', 'Summary figure (PNG)']], r));
  host.scrollIntoView({ behavior: 'smooth', block: 'start' });
}
function lightCurve(host, pts, ph) {
  const used = pts.filter((p) => p.status === 'USED'); const W = chartWidth(host, 900), H = 280, m = { l: 70, r: 20, t: 14, b: 40 };
  const xs = pts.map((p) => p.mjd); const x0 = Math.min(...xs) - 15, x1 = Math.max(...xs) + 15;
  const ys = used.flatMap((p) => [p.calibrated_flux_jy - p.flux_sigma_jy, p.calibrated_flux_jy + p.flux_sigma_jy].map((v) => v * 1e3)); let y0 = Math.min(...ys), y1 = Math.max(...ys); const pad = (y1 - y0) * 0.35 || y1 * 0.05; y0 -= pad; y1 += pad;
  const X = (v) => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r), Y = (v) => H - m.b - (v - y0) / (y1 - y0) * (H - m.t - m.b);
  const svg = sv('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'Same-wavelength light curve' });
  const grid = sv('g', { class: 'grid' }), axis = sv('g', { class: 'axis' });
  for (let i = 0; i <= 4; i++) { const v = y0 + (y1 - y0) * i / 4; grid.append(sv('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v) })); axis.append(sv('text', { x: m.l - 8, y: Y(v) + 4, 'text-anchor': 'end' }, v.toFixed(2))); }
  monthTicks(x0, x1, Math.floor(W / 95)).forEach((t) => { axis.append(sv('text', { x: X(t.mjd), y: H - m.b + 18, 'text-anchor': 'middle' }, t.label)); grid.append(sv('line', { x1: X(t.mjd), x2: X(t.mjd), y1: m.t, y2: H - m.b })); });
  axis.append(sv('text', { x: 16, y: (H - m.b) / 2, transform: `rotate(-90 16 ${(H - m.b) / 2})`, 'text-anchor': 'middle', class: 'axis-title' }, 'Flux density (mJy)'));
  svg.append(grid, axis);
  if (ph.constant_flux_test) { const c = ph.constant_flux_test.constant_flux * 1e3; svg.append(sv('line', { x1: m.l, x2: W - m.r, y1: Y(c), y2: Y(c), stroke: 'var(--text-muted)', 'stroke-dasharray': '4 4', 'stroke-width': 1.5 })); svg.append(sv('text', { x: m.l + 6, y: Y(c) - 6, fill: 'var(--text-muted)', 'font-size': 11 }, 'constant-flux fit')); }
  used.forEach((p) => {
    const f = p.calibrated_flux_jy * 1e3, s = p.flux_sigma_jy * 1e3; const g = sv('g', {});
    g.append(sv('line', { x1: X(p.mjd), x2: X(p.mjd), y1: Y(f - s), y2: Y(f + s), stroke: 'var(--series-1)', 'stroke-width': 2 }));
    g.append(sv('circle', { cx: X(p.mjd), cy: Y(f), r: 5, fill: 'var(--series-1)', stroke: 'var(--surface-1)', 'stroke-width': 2 }));
    g.append(sv('circle', { cx: X(p.mjd), cy: Y(f), r: 14, fill: 'transparent' }));
    g.addEventListener('mouseenter', (e) => showTip(e, `<b>${p.date.slice(0, 10)}</b> · ${esc(short(p.observation_id))}<br>${fmt(f, 3)} ± ${fmt(s, 3)} mJy at ${fmt(p.wavelength_um, 4)} µm<br><span class="muted">zero point ${fmt(p.zero_point, 4)}</span>`));
    g.addEventListener('mousemove', moveTip); g.addEventListener('mouseleave', hideTip); svg.append(g);
  });
  host.append(svg);
}
function trajectory(host, as) {
  const rows = as.positions; const W = chartWidth(host, 900), H = 360, m = { l: 70, r: 20, t: 14, b: 44 };
  const ex = rows.flatMap((r) => [r.east_arcsec, r.expected_east_arcsec].filter(Number.isFinite)), ny = rows.flatMap((r) => [r.north_arcsec, r.expected_north_arcsec].filter(Number.isFinite));
  const cx = (Math.max(...ex) + Math.min(...ex)) / 2, cy = (Math.max(...ny) + Math.min(...ny)) / 2;
  const span = Math.max(Math.max(...ex) - Math.min(...ex), Math.max(...ny) - Math.min(...ny)) * 0.65 + 0.6;
  const aspect = (W - m.l - m.r) / (H - m.t - m.b);
  const X = (v) => m.l + (W - m.l - m.r) / 2 - (v - cx) / (span * aspect) * (W - m.l - m.r) / 2; // east to the left
  const Y = (v) => m.t + (H - m.t - m.b) / 2 - (v - cy) / span * (H - m.t - m.b) / 2;
  const svg = sv('svg', { viewBox: `0 0 ${W} ${H}`, role: 'img', 'aria-label': 'Astrometric trajectory: measured versus catalog prediction' });
  const grid = sv('g', { class: 'grid' }), axis = sv('g', { class: 'axis' });
  for (let i = -2; i <= 2; i++) { const v = cy + span * i / 2; grid.append(sv('line', { x1: m.l, x2: W - m.r, y1: Y(v), y2: Y(v) })); axis.append(sv('text', { x: m.l - 8, y: Y(v) + 4, 'text-anchor': 'end' }, (v - cy).toFixed(1))); }
  for (let i = -4; i <= 4; i++) { const v = cx + span * aspect * i / 4; grid.append(sv('line', { x1: X(v), x2: X(v), y1: m.t, y2: H - m.b })); axis.append(sv('text', { x: X(v), y: H - m.b + 16, 'text-anchor': 'middle' }, (v - cx).toFixed(1))); }
  axis.append(sv('text', { x: (W + m.l) / 2, y: H - 8, 'text-anchor': 'middle', class: 'axis-title' }, 'East offset (arcsec, east to the left)'));
  axis.append(sv('text', { x: 16, y: (H - m.b) / 2, transform: `rotate(-90 16 ${(H - m.b) / 2})`, 'text-anchor': 'middle', class: 'axis-title' }, 'North offset (arcsec)'));
  svg.append(grid, axis);
  if (Number.isFinite(rows[0].expected_east_arcsec)) {
    svg.append(sv('polyline', { points: rows.map((r) => `${X(r.expected_east_arcsec)},${Y(r.expected_north_arcsec)}`).join(' '), fill: 'none', stroke: 'var(--series-2)', 'stroke-width': 2, 'stroke-dasharray': '5 4' }));
    rows.forEach((r) => svg.append(sv('rect', { x: X(r.expected_east_arcsec) - 5, y: Y(r.expected_north_arcsec) - 5, width: 10, height: 10, rx: 2, fill: 'var(--surface-1)', stroke: 'var(--series-2)', 'stroke-width': 2 })));
  }
  rows.forEach((r) => {
    const g = sv('g', {}); const x = X(r.east_arcsec), y = Y(r.north_arcsec);
    const sx = Math.abs(X(r.east_arcsec + r.sigma_east_arcsec) - x), sy = Math.abs(Y(r.north_arcsec + r.sigma_north_arcsec) - y);
    g.append(sv('line', { x1: x - sx, x2: x + sx, y1: y, y2: y, stroke: 'var(--series-1)', 'stroke-width': 2 }), sv('line', { x1: x, x2: x, y1: y - sy, y2: y + sy, stroke: 'var(--series-1)', 'stroke-width': 2 }));
    g.append(sv('circle', { cx: x, cy: y, r: 5, fill: 'var(--series-1)', stroke: 'var(--surface-1)', 'stroke-width': 2 }));
    g.append(sv('text', { x: x + 8, y: y - 8, fill: 'var(--text-secondary)', 'font-size': 11 }, r.date.slice(0, 10)));
    g.append(sv('circle', { cx: x, cy: y, r: 14, fill: 'transparent' }));
    g.addEventListener('mouseenter', (e) => showTip(e, `<b>${r.date.slice(0, 10)}</b><br>measured ${fmt(r.east_arcsec - cx, 3)}″ E, ${fmt(r.north_arcsec - cy, 3)}″ N ± ${fmt(r.sigma_east_arcsec, 3)}″${Number.isFinite(r.expected_east_arcsec) ? `<br><span class="muted">catalog ${fmt(r.expected_east_arcsec - cx, 3)}″ E, ${fmt(r.expected_north_arcsec - cy, 3)}″ N</span>` : ''}`));
    g.addEventListener('mousemove', moveTip); g.addEventListener('mouseleave', hideTip); svg.append(g);
  });
  host.append(svg, el('div', { class: 'legend' }, el('span', {}, el('i', { style: 'background:var(--series-1)' }), 'Measured (SPHEREx, frame-corrected)'), el('span', {}, el('i', { class: 'square', style: 'background:var(--series-2)' }), 'Catalog prediction (proper motion + parallax)')));
}

/* ------------------------------------------------------------------ recorded runs */
let runsLoaded = false;
async function loadRuns() {
  if (runsLoaded) return; const host = $('runs');
  try {
    const runs = await api('/api/lab/demos'); runsLoaded = true; host.replaceChildren();
    if (!runs.length) { host.append(el('div', { class: 'empty' }, 'No recorded runs yet. Generate them with: python make_demos.py')); return; }
    for (const run of runs) {
      const ok = run.status === 'COMPLETED'; const s = run.summary || {};
      const facts = run.kind === 'compare' && ok ? [['Status', s.status], ['Matched pixels', s.matched_pixels], ['Δt (days)', fmt(s.separation_days, 1)], ['Target', (s.target_status || '').replaceAll('_', ' ')]]
        : run.kind === 'series' && ok ? [['Light curve', s.photometry_status], ['Motion', s.astrometry_status], ['PM (″/yr)', s.proper_motion_arcsec_per_year ? s.proper_motion_arcsec_per_year.map((v) => fmt(v, 2)).join(', ') : '—'], ['Catalog PM', s.catalog_proper_motion ? s.catalog_proper_motion.map((v) => fmt(v, 2)).join(', ') : '—']]
        : run.kind === 'screen' && ok ? Object.entries(s).map(([t, v]) => [t, `${v.matched_pairs} pairs · ${v.groups_3plus} groups`]) : [['Error', run.error || '—']];
      host.append(el('article', { class: 'card' },
        el('div', { style: 'display:flex;justify-content:space-between;gap:8px;align-items:start' }, el('h3', {}, run.title), el('span', { class: 'badge ' + (ok ? 'good' : 'critical') }, ok ? run.kind : 'failed')),
        el('p', {}, run.why), el('dl', { class: 'kv' }, facts.flatMap(([k, v]) => [el('dt', {}, k), el('dd', {}, v === undefined ? '—' : v)])),
        el('p', { class: 'small muted' }, `Generated ${run.generated_utc} · ${run.elapsed_seconds}s · ${run.catalog_mode || ''}`),
        ok ? el('button', { class: 'btn primary', onclick: () => openRun(run) }, 'Open') : ''));
    }
  } catch (e) { host.replaceChildren(el('div', { class: 'notice error' }, 'Could not load recorded runs: ' + e.message)); }
}
async function openRun(run) {
  const base = `/api/lab/demos/${run.name}/files/`; const host = $('run-result'); host.replaceChildren(el('div', { class: 'empty' }, 'Loading…'));
  try {
    if (run.kind === 'compare') await renderComparison(await api(base + 'report.json'), base, host);
    else if (run.kind === 'series') renderSeries(await api(base + 'series.json'), base, host);
    else {
      host.replaceChildren(); const pick = el('div', { class: 'seg', style: 'margin-bottom:14px' }); const inner = el('div', { class: 'content' });
      const files = run.files.filter((f) => f.startsWith('screening_'));
      const show = async (f, btn) => { [...pick.children].forEach((b) => b.setAttribute('aria-pressed', String(b === btn))); inner.replaceChildren(el('div', { class: 'empty' }, 'Loading…')); const data = await api(base + f); inner.replaceChildren(); renderScreening(data, inner, false); };
      files.forEach((f, i) => { const b = el('button', { 'aria-pressed': 'false' }, f.replace('screening_', '').replace('.json', '')); b.addEventListener('click', () => show(f, b)); pick.append(b); if (i === 0) setTimeout(() => show(f, b)); });
      host.append(pick, inner);
    }
  } catch (e) { host.replaceChildren(el('div', { class: 'notice error' }, 'Could not open run: ' + e.message)); }
}

/* ------------------------------------------------------------------ method */
function renderMethod() {
  const m = $('method'); if (m.childElementCount) return;
  m.innerHTML = `
  <h2>The problem</h2>
  <p>NASA's SPHEREx maps the whole sky in 102 colours between 0.75 and 5.0 µm by sweeping linear variable filters across the sky. The same star is therefore observed many times, but almost always <em>through a different colour</em>. Comparing two SPHEREx images pixel by pixel compares different wavelengths unless that is explicitly ruled out, and a star's colour difference can then masquerade as variability.</p>
  <h2>What the Epoch Lab does</h2>
  <ol>
    <li><b>Screens the archive.</b> Every DR1 record covering the position is mapped onto its detector with the archive footprint polygon (median error 3.1 px, maximum 5.7 px against 26 real cutouts) and the detector's DataLink wavelength map. The result is the wavelength SPHEREx actually sampled at the target in each visit.</li>
    <li><b>Finds same-wavelength epochs.</b> Pairs on the same detector, on different dates, whose wavelengths agree within 10% of the local bandwidth. Revisits of the same pointing (weeks or a year apart) give whole-field matches; rotated revisits leave only a strip.</li>
    <li><b>Compares real images.</b> Cutouts and each observation's own calibration are downloaded. Each epoch gets a clipped background plane; the measured ePSFs are cross-convolved on common sky axes; epoch B is reprojected onto A with squared-weight variance and A receives the same interpolation smoothing. Only pixels whose own wavelengths match are compared.</li>
    <li><b>Validates every result.</b> Registration from hundreds of Gaia stars, a photometric scale from high-S/N stars, a blank-sky null test that can only inflate uncertainties, an empirical source noise model (multiplicative + additive), and injection–recovery of synthetic sources in the real difference image.</li>
    <li><b>Classifies residuals conservatively.</b> Known proper motion, PSF/registration dipoles, the systematic floor, persistent unchanged sources and edge effects are explained first. Only what remains is called a candidate, and candidates must pass vetting: isolated from flagged pixels, not blended, PSF-shaped.</li>
    <li><b>Builds time series.</b> For three or more same-wavelength epochs: a light curve with field-star zero points, and an astrometric trajectory compared with catalog proper motion and parallax using JPL Horizons positions of SPHEREx.</li>
  </ol>
  <h2>What it will not claim</h2>
  <p>A residual is not a discovery. Moving-object identification requires a measured multi-epoch trajectory; catalog-forced positions are not detections; different-wavelength epochs are rejected rather than corrected with an invented spectral model. Saturated or flagged targets are reported as unmeasurable. When an external service (IRSA, Gaia, JPL) fails, the job fails with that reason — nothing is substituted.</p>
  <h2>Policies</h2>
  <p>Spectral tolerance 0.10 × CBAND; detection 5σ with ≥2 connected pixels; registration tolerance 0.3 px; aperture 5×5 native pixels; null test needs ≥200 source-free pixels; source noise model needs ≥30 stars. These are analysis choices recorded in every report, not mission calibration constants.</p>
  <h2>Data</h2>
  <p>NASA/IPAC IRSA SPHEREx Data Release 1 (SIA search, cutout service, DataLink calibration products and ePSF tables), ESA Gaia DR3, SIMBAD/CDS Sesame and JPL Horizons. See the repository's <code>docs/</code> for full methodology and reproduction instructions.</p>`;
}

/* ------------------------------------------------------------------ boot */
loadTargets();
const routeFromHash = () => { const v = (location.hash || '#explore').slice(1); showView(['explore', 'runs', 'method'].includes(v) ? v : 'explore'); };
window.addEventListener('hashchange', routeFromHash);
routeFromHash();
