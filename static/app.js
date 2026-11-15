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

/* A static copy (scripts/build_static_demo.py) has no server: recorded runs are plain files and live analysis is unavailable. */
const STATIC = document.documentElement.getAttribute('data-static') === 'true';
function staticPath(path) {
  if (path === '/api/lab/targets') return 'data/targets.json';
  if (path === '/api/lab/demos') return 'data/demos.json';
  const m = /^\/api\/lab\/demos\/([\w-]+)\/files\/(.+)$/.exec(path);
  return m ? `data/demos/${m[1]}/${m[2]}` : null;
}
const demoBase = (name) => (STATIC ? `data/demos/${name}/` : `/api/lab/demos/${name}/files/`);
function apiError(message, code, status) { const err = new Error(message); err.code = code; err.status = status; return err; }

async function api(path, options = {}) {
  if (STATIC && path.startsWith('/api/')) {
    path = staticPath(path);
    if (!path) throw apiError('Live analysis needs the Epoch Lab server; this public copy shows recorded runs only.', 'STATIC_DEMO');
  }
  let response;
  try { response = await fetch(path, options); } catch (e) { throw apiError(`Could not reach ${STATIC ? 'the page files' : 'the Epoch Lab server'} (${e.message}).`, 'SERVER_UNREACHABLE'); }
  let body = null;
  try { body = await response.json(); } catch (_) { /* non-JSON */ }
  if (!response.ok) {
    const d = body && body.detail;
    let message = typeof d === 'string' ? d : (d && (d.message || d.code)) || `HTTP ${response.status}`;
    if (Array.isArray(d)) message = d.map((e) => e.msg).join('; ');
    throw apiError(message, (d && d.code) || (response.status === 422 ? 'INVALID_INPUT' : `HTTP_${response.status}`), response.status);
  }
  return body;
}

/* ------------------------------------------------------------------ tooltip */
const tip = $('tooltip');
function showTip(evt, html) { tip.innerHTML = html; tip.hidden = false; moveTip(evt); }
function moveTip(evt) { const pad = 14; let x = evt.clientX + pad, y = evt.clientY + pad; const r = tip.getBoundingClientRect(); if (x + r.width > innerWidth - 8) x = evt.clientX - r.width - pad; if (y + r.height > innerHeight - 8) y = evt.clientY - r.height - pad; tip.style.left = x + 'px'; tip.style.top = y + 'px'; }
function hideTip() { tip.hidden = true; }
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const human = (s) => { const t = String(s || '').replaceAll('_', ' ').toLowerCase(); return t.charAt(0).toUpperCase() + t.slice(1); };
const count = (n) => (typeof n === 'number' ? n.toLocaleString('en-US') : n);
const microns = (s) => String(s || '').replace(/(\d) um\b/g, '$1 µm');

/* ------------------------------------------------------------------ shared parts */
function tabset(defs) {
  const bar = el('div', { class: 'tabbar', role: 'tablist' }), body = el('div', {});
  const show = (id) => { for (const d of defs) { d.button.setAttribute('aria-selected', String(d.id === id)); d.node.hidden = d.id !== id; } };
  for (const d of defs) {
    d.button = el('button', { role: 'tab', type: 'button', onclick: () => show(d.id) }, d.label, d.count === undefined ? null : el('span', { class: 'count' + (d.flag ? ' flag' : '') }, count(d.count)));
    bar.append(d.button); body.append(d.node);
  }
  show(defs[0].id);
  return { bar, body, show };
}
function tableOf(heads, numeric, rows, bodyAttrs = {}) {
  return el('div', { class: 'table-wrap' }, el('table', {}, el('thead', {}, el('tr', {}, heads.map((h, i) => el('th', { class: numeric.includes(i) ? 'num' : null }, h)))), el('tbody', bodyAttrs, rows)));
}
function backBar(back) { return el('button', { class: 'back', type: 'button', onclick: back.go }, '← ' + back.label); }
function clearHost(host) { if (host.viewer) clearInterval(host.viewer.timer); host.viewer = null; host.token = null; host.replaceChildren(); }
/* Plain-language meaning of each error code; the server's own message is shown underneath. */
const ERROR_HELP = {
  ARCHIVE_UNAVAILABLE: 'NASA’s IRSA archive is not answering right now. The outage is on their side and usually clears within hours.',
  ARCHIVE_TIMEOUT: 'A NASA/IRSA request took too long and was stopped after a limited number of retries.',
  CALIBRATION_ARCHIVE_UNAVAILABLE: 'IRSA could not supply the wavelength calibration for an observation, so it cannot be compared safely.',
  GAIA_CATALOG_UNAVAILABLE: 'The ESA Gaia archive, used to check the comparison against catalog stars, is not answering.',
  QUEUE_FULL: 'Several analyses are already running. Wait for one to finish, then try again.',
  WORKER_RESTARTED: 'The server restarted before this analysis finished. Run it again.',
  ANALYSIS_REJECTED: 'The analysis stopped because the data did not meet one of its requirements.',
  INVALID_INPUT: 'The request was not valid.',
  SERVER_UNREACHABLE: 'The Epoch Lab server is not responding. Check that it is still running in its terminal window.',
  STATIC_DEMO: 'This public copy shows recorded runs only.',
};
const OUTAGES = new Set(['ARCHIVE_UNAVAILABLE', 'ARCHIVE_TIMEOUT', 'CALIBRATION_ARCHIVE_UNAVAILABLE', 'GAIA_CATALOG_UNAVAILABLE']);
function errorNotice(message, code) {
  const help = ERROR_HELP[code] || 'The analysis could not be completed.';
  const outage = OUTAGES.has(code);
  return el('div', { class: 'notice error' },
    el('p', {}, el('b', {}, help), outage ? ' Nothing was substituted for the missing data. Recorded runs work without the network.' : ''),
    message && message !== help ? el('p', { class: 'detail' }, message) : null,
    el('div', { class: 'notice-foot' }, el('code', {}, code || ''),
      state.retry ? el('button', { class: 'btn small', type: 'button', onclick: () => state.retry() }, 'Try again') : null,
      outage ? el('a', { class: 'btn small', href: '#runs' }, 'Open a recorded run') : null,
      el('button', { class: 'btn ghost small', type: 'button', onclick: () => { $('job-panel').hidden = true; } }, 'Dismiss')));
}

/* ------------------------------------------------------------------ navigation */
const VIEWS = ['explore', 'runs', 'method'];
function activate(name) {
  for (const b of document.querySelectorAll('.tabs [data-view]')) b.setAttribute('aria-selected', String(b.dataset.view === name));
  for (const v of VIEWS) $('view-' + v).hidden = v !== name;
  if (name === 'method') renderMethod();
}
function route() {
  const [view, arg] = (location.hash || '#explore').slice(1).split('/');
  const name = VIEWS.includes(view) ? view : 'explore';
  activate(name);
  if (name === 'runs') showRuns(arg); else window.scrollTo({ top: 0 });
}
function go(hash) { if (location.hash === '#' + hash) route(); else location.hash = hash; }
document.querySelectorAll('.tabs [data-view]').forEach((b) => b.addEventListener('click', () => go(b.dataset.view)));
document.addEventListener('click', (e) => { const g = e.target.closest('[data-goto]'); if (g) { e.preventDefault(); go(g.dataset.goto); } });

/* ------------------------------------------------------------------ field selection */
const state = { field: null, screening: null, selected: new Set(), targets: [], stage: 'start', retry: null, resolving: false };
const HINT = 'Pick a preset, type a name such as Vega, or enter coordinates in degrees such as 269.452 4.693.';
const PRESET_NOTE = 'preset, proper motion applied';
const COORDS = /^\s*([+-]?\d+(?:\.\d+)?)\s*°?\s*[,\s]\s*([+-]?\d+(?:\.\d+)?)\s*°?\s*$/;
let targetsReady = Promise.resolve();
const targetName = (id) => (state.targets.find((t) => t.id === id) || {}).name || id;

function fieldStatus(text, tone) { const p = $('field-label'); p.className = 'field-status' + (tone ? ' ' + tone : ''); p.textContent = text; }
function syncFinder(busy = !!activeJob) { $('screen').disabled = busy || !(state.field || $('query').value.trim()); }
function setField(field) {
  state.field = field;
  for (const c of document.querySelectorAll('#presets .chip')) c.setAttribute('aria-pressed', String(!!field && c.dataset.id === field.target_id));
  if (!field) fieldStatus(HINT);
  else { fieldStatus('', 'set'); $('field-label').append(el('b', {}, field.label), ` · RA ${field.ra.toFixed(5)}°, Dec ${field.dec.toFixed(5)}°${field.note ? ' · ' + field.note : ''}`); }
  syncFinder();
}
async function loadTargets() {
  try {
    state.targets = await api('/api/lab/targets');
    $('presets').replaceChildren(...state.targets.map((t) => el('button', { class: 'chip', type: 'button', 'data-id': t.id, 'aria-pressed': 'false', onclick: () => pickPreset(t) }, t.name)));
  } catch (e) { fieldStatus('Could not load presets: ' + e.message, 'error'); }
}
function pickPreset(t) {
  $('query').value = t.name;
  setField({ target_id: t.id, ra: t.ra_deg, dec: t.dec_deg, label: t.name, note: PRESET_NOTE });
  $('screen').focus();
}
function parseCoords(text) {
  const m = COORDS.exec(text.replace(/[−–]/g, '-')); if (!m) return null;
  const ra = parseFloat(m[1]), dec = parseFloat(m[2]);
  return { ra, dec, valid: ra >= 0 && ra < 360 && dec >= -90 && dec <= 90 };
}
$('query').addEventListener('input', () => {
  const text = $('query').value.trim(); const c = parseCoords(text);
  if (c && c.valid) return setField({ ra: c.ra, dec: c.dec, label: 'Fixed coordinates' });
  setField(null);
  if (c) fieldStatus('RA must be between 0 and 360°, Dec between −90 and 90°.', 'error');
  else if (text) fieldStatus(`Press Enter to look up “${text}”.`);
});
async function resolveQuery() {
  if (state.field) return true;
  const name = $('query').value.trim(); if (!name || parseCoords(name)) return false;
  fieldStatus(`Looking up “${name}”…`);
  try {
    const r = await api('/api/lab/resolve?name=' + encodeURIComponent(name));
    setField(r.target_id ? { target_id: r.target_id, ra: r.ra_deg, dec: r.dec_deg, label: r.name, note: PRESET_NOTE } : { ra: r.ra_deg, dec: r.dec_deg, label: r.name, note: r.source });
    return true;
  } catch (e) { fieldStatus(e.message, 'error'); return false; }
}
$('query').addEventListener('keydown', (e) => { if (e.key === 'Enter') findEpochs(); });
function fieldBody(f = state.field) { return f.target_id ? { target_id: f.target_id } : { ra: f.ra, dec: f.dec }; }

/* ------------------------------------------------------------------ jobs */
const STAGES = {
  screen: [['archive_discovery', 'Query IRSA SPHEREx DR1 archive'], ['calibration_maps', 'Resolve detector wavelength maps (DataLink)'], ['screening', 'Estimate wavelength at target for every record'], ['completed', 'Rank same-wavelength epochs']],
  compare: [['archive_discovery', 'Query archive records'], ['catalog', 'Gaia DR3 field catalog'], ['download', 'Retrieve real FITS cutouts'], ['calibration', 'Own-calibration wavelength maps'], ['comparison', 'Background, PSF match, register, difference, validate'], ['products', 'Write FITS / PNG / report'], ['completed', 'Done']],
  series: [['archive_discovery', 'Query archive records'], ['catalog', 'Gaia DR3 field catalog'], ['download', 'Retrieve real FITS cutouts'], ['calibration', 'Own-calibration wavelength maps'], ['observer_ephemeris', 'JPL Horizons SPHEREx position'], ['measurement', 'ePSF photometry + astrometry'], ['completed', 'Done']],
};
const JOB_TITLE = { screen: 'Screening the archive', compare: 'Comparing two epochs', series: 'Building the time series' };
const JOB_BADGE = { COMPLETED: ['good', 'Completed'], PARTIAL: ['warning', 'Partial result'], FAILED: ['critical', 'Failed'], INTERRUPTED: ['critical', 'Interrupted'], RUNNING: ['neutral', 'Running'], QUEUED: ['neutral', 'Queued'] };
let activeJob = null;
function renderProgress(kind, job) {
  const steps = STAGES[kind]; const failed = ['FAILED', 'INTERRUPTED'].includes(job.status); const done = ['COMPLETED', 'PARTIAL'].includes(job.status);
  let reached = steps.findIndex((s) => s[0] === (job.progress && job.progress.stage)); if (done) reached = steps.length;
  const panel = $('job-panel'); panel.hidden = false; panel.className = 'job' + (failed ? ' failed' : done ? '' : ' running');
  $('job-title').textContent = JOB_TITLE[kind];
  const badge = JOB_BADGE[job.status] || ['neutral', job.status];
  $('job-badge').replaceChildren(el('span', { class: 'badge ' + badge[0] }, badge[1]));
  $('job-fill').style.width = Math.round(100 * Math.max(reached, 0.25) / steps.length) + '%';
  const step = steps[Math.min(Math.max(reached, 0), steps.length - 1)];
  const obs = job.progress && job.progress.observation_id ? ` · ${short(job.progress.observation_id)}` : '';
  $('job-now').textContent = done ? 'Finished.' : reached < 0 ? (failed ? 'Stopped before the first step.' : 'Waiting for the worker…')
    : `${failed ? 'Stopped at step' : 'Step'} ${reached + 1} of ${steps.length} · ${step[1]}${obs}`;
  $('job-steps').replaceChildren(...steps.map(([, text], i) => el('li', { class: i < reached ? 'done' : (i === reached && !failed ? 'active' : null) }, text)));
  $('job-error').replaceChildren(...(job.error ? [errorNotice(job.error.message || String(job.error), job.error.code)] : []));
}
async function runJob(kind, body) {
  if (activeJob) return null;
  activeJob = 'submitting'; setBusy(true);
  try {
    const sub = await api('/api/lab/' + kind, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) });
    activeJob = sub.id;
    renderProgress(kind, { status: 'QUEUED', progress: { stage: 'queued' } });
    let last = { status: 'QUEUED', progress: { stage: 'queued' } };
    for (let misses = 0; ;) {
      await new Promise((r) => setTimeout(r, 1200));
      let job;
      try { job = await api('/api/lab/jobs/' + sub.id); misses = 0; } catch (e) {
        // Tolerate a few dropped polls; after ~15 s of silence report it instead of spinning forever.
        if (++misses < 12) continue;
        renderProgress(kind, { ...last, status: 'FAILED', error: { message: e.message, code: e.code } });
        return null;
      }
      last = job; renderProgress(kind, job);
      if (job.status === 'COMPLETED') { $('job-panel').hidden = true; return { job, base: `/api/lab/jobs/${sub.id}/files/` }; }
      if (['FAILED', 'INTERRUPTED'].includes(job.status)) return null;
    }
  } catch (e) {
    renderProgress(kind, { status: 'FAILED', progress: {}, error: { message: e.message, code: e.code } });
    $('job-badge').replaceChildren(el('span', { class: 'badge critical' }, 'Not started'));
    return null;
  } finally { activeJob = null; setBusy(false); }
}
function setBusy(busy) { document.body.classList.toggle('busy', busy); syncFinder(busy); updateSelection(busy); }

/* ------------------------------------------------------------------ screening */
function showStage(stage) {
  state.stage = stage;
  $('welcome').hidden = $('start').hidden = stage !== 'start';
  $('screening').hidden = stage !== 'epochs';
  $('result').hidden = stage !== 'result';
  updateSelection();
}
async function findEpochs() {
  if (activeJob || state.resolving) return;
  state.resolving = true; syncFinder(true);
  const ok = await resolveQuery();
  state.resolving = false; syncFinder();
  if (!ok) return;
  state.retry = findEpochs;
  const field = state.field;
  const out = await runJob('screen', fieldBody(field));
  if (!out) return;
  state.screenedField = field;
  state.selected.clear(); showStage('epochs');
  renderScreening(out.job.result, $('screening'), true);
  window.scrollTo({ top: 0 });
}
$('screen').addEventListener('click', findEpochs);

function renderScreening(s, host, interactive) {
  host.hidden = false; host.replaceChildren();
  if (interactive) state.screening = s;
  const screened = s.observations.filter((o) => o.status === 'SCREENED');
  const groups3 = s.matched_groups.filter((g) => g.epoch_count >= 3);
  const name = s.target_id ? targetName(s.target_id) : (interactive && state.screenedField ? state.screenedField.label : 'Fixed coordinates');
  host.append(el('header', { class: 'page-head' }, el('p', { class: 'eyebrow' }, 'Archive coverage'), el('h2', {}, name),
    el('p', { class: 'sub' }, `RA ${fmt(s.query.ra_deg, 5)}°, Dec ${fmt(s.query.dec_deg, 5)}° · search radius ${s.query.radius_deg}° · SPHEREx DR1`)),
  el('div', { class: 'stats' },
    stat(s.observation_count, 'archive records', `${s.rejected_records.length} outside the footprint`),
    stat(s.screened_count, 'wavelength screened', 'detectors D1 to D6'),
    stat(s.matched_pairs.length, 'same-wavelength pairs', `|Δλ| ≤ ${s.policy.tolerance_fraction} × bandwidth`),
    stat(groups3.length, 'groups of 3 or more epochs', 'light-curve candidates')));
  if (!s.observation_count) {
    host.append(el('div', { class: 'notice info' }, el('p', {}, el('b', {}, 'SPHEREx DR1 has no observation covering this position. '),
      `The archive search found no record whose footprint contains RA ${fmt(s.query.ra_deg, 4)}°, Dec ${fmt(s.query.dec_deg, 4)}°, so there is nothing to compare. `,
      'This describes the archive’s current coverage, not the sky: try another position or open a recorded run.')));
    return;
  }
  const chartHost = el('div', { class: 'chart' });
  host.append(el('section', { class: 'panel' },
    el('div', { class: 'panel-head' }, el('h3', {}, 'Wavelength sampled at the target, by date'),
      el('div', { class: 'legend' }, el('span', {}, el('i', { style: 'background:var(--series-1)' }), 'Observation'), el('span', {}, el('i', { style: 'background:var(--series-2)' }), 'In a same-wavelength pair'))),
    chartHost,
    el('p', { class: 'footnote' }, `Bands D1 to D6 are the six detectors.${interactive ? ' Click points to select observations.' : ''} Wavelengths are screening estimates from the footprint model (≤${s.policy.footprint_error_px} px); every comparison re-checks them on the downloaded cutouts.`)));
  coverageChart(chartHost, screened, s.matched_pairs, interactive);
  const tabs = tabset([
    { id: 'pairs', label: 'Pairs', count: s.matched_pairs.length, node: pairsPanel(s, interactive) },
    { id: 'groups', label: 'Groups', count: s.matched_groups.length, node: groupsPanel(s, interactive) },
    { id: 'all', label: 'All observations', count: s.observations.length, node: observationsPanel(s, interactive) },
  ]);
  host.append(tabs.bar, tabs.body);
  if (interactive) updateSelection();
}
function stat(value, label, sub) { return el('div', { class: 'stat' }, el('div', { class: 'value' }, count(value)), el('div', { class: 'label' }, label), sub ? el('div', { class: 'sub' }, sub) : null); }
function pairsPanel(s, interactive) {
  const shown = s.matched_pairs.slice(0, 60);
  const panel = el('div', { class: 'tabpanel' }, el('p', { class: 'intro' }, 'Same detector, different dates, wavelength at the target within the tolerance. Ranked by the predicted share of the 0.3° field that matches.'));
  if (!shown.length) { panel.append(el('div', { class: 'empty' }, 'No same-wavelength pair at this position in the searched records. Different-wavelength epochs are not compared.')); return panel; }
  panel.append(tableOf(['Epoch A', 'Epoch B', 'Det.', 'λ (µm)', '|Δλ|/band', 'Δt (days)', 'Matched area', ''], [2, 3, 4, 5, 6], shown.map((p) => {
    const pct = p.predicted_matched_fraction == null ? null : Math.round(100 * p.predicted_matched_fraction);
    return el('tr', {}, el('td', { class: 'mono' }, short(p.observation_a)), el('td', { class: 'mono' }, short(p.observation_b)), el('td', { class: 'num' }, 'D' + p.detector),
      el('td', { class: 'num' }, fmt(p.wavelength_um, 4)), el('td', { class: 'num' }, fmt(p.delta_over_bandwidth, 3)), el('td', { class: 'num' }, fmt(p.separation_days, 1)),
      el('td', { class: 'num' }, pct === null ? '—' : el('span', { class: 'meter' }, el('i', {}, el('b', { style: `width:${pct}%` })), pct + '%')),
      el('td', { class: 'num' }, interactive ? el('button', { class: 'btn small needs-idle', type: 'button', onclick: () => startCompare(p.observation_a, p.observation_b) }, 'Compare') : null));
  })));
  if (s.matched_pairs.length > shown.length) panel.append(el('p', { class: 'footnote' }, `Showing the ${shown.length} best of ${s.matched_pairs.length} pairs.`));
  return panel;
}
function groupsPanel(s, interactive) {
  const shown = s.matched_groups.slice(0, 40);
  const panel = el('div', { class: 'tabpanel' }, el('p', { class: 'intro' }, 'Several epochs on one detector at one wavelength. Three or more make a light curve and an astrometric trajectory possible.'));
  if (!shown.length) { panel.append(el('div', { class: 'empty' }, 'No multi-epoch group at a common wavelength.')); return panel; }
  panel.append(tableOf(['Det.', 'λ (µm)', 'Epochs', 'Dates', ''], [0, 1, 2], shown.map((g) => el('tr', {},
    el('td', { class: 'num' }, 'D' + g.detector), el('td', { class: 'num' }, fmt(g.wavelength_um, 4)), el('td', { class: 'num' }, g.epoch_count), el('td', { class: 'mono small' }, g.mjd.map(isoDay).join(', ')),
    el('td', { class: 'num' }, interactive ? el('button', { class: 'btn small needs-idle', type: 'button', onclick: () => startSeries(g.observation_ids.slice(0, 8)) }, 'Time series') : null)))));
  if (s.matched_groups.length > shown.length) panel.append(el('p', { class: 'footnote' }, `Showing ${shown.length} of ${s.matched_groups.length} groups.`));
  return panel;
}
function observationsPanel(s, interactive) {
  const rows = s.observations.map((o) => {
    const box = interactive && o.status === 'SCREENED' ? el('input', { type: 'checkbox', 'data-id': o.observation_id, 'aria-label': 'Select ' + short(o.observation_id) }) : null;
    if (box) box.addEventListener('change', () => toggleSelect(o.observation_id, box.checked));
    return el('tr', { 'data-row': o.observation_id }, el('td', {}, box), el('td', { class: 'mono' }, o.mjd ? isoDay(o.mjd) : '—'), el('td', { class: 'mono' }, short(o.observation_id)),
      el('td', { class: 'num' }, o.detector ? 'D' + o.detector : '—'), el('td', { class: 'num' }, fmt(o.wavelength_um, 4)), el('td', { class: 'num' }, fmt(o.bandwidth_um, 4)),
      el('td', { class: 'num' }, fmt(o.wavelength_uncertainty_um, 4)), el('td', { class: 'small secondary' }, o.status === 'SCREENED' ? 'screened' : (o.reason || o.status)));
  });
  return el('div', { class: 'tabpanel' }, el('p', { class: 'intro' }, interactive ? 'Tick two observations on one detector to compare them, or up to eight for a time series.' : 'Every screened archive record.'),
    tableOf(['', 'Date', 'Observation', 'Det.', 'λ at target (µm)', 'Bandwidth (µm)', '± screening (µm)', 'Status'], [3, 4, 5, 6], rows, { id: interactive ? 'obs-body' : null }));
}

function toggleSelect(id, on) {
  if (on) state.selected.add(id); else state.selected.delete(id);
  const box = document.querySelector(`#obs-body input[data-id="${CSS.escape(id)}"]`);
  if (box) { box.checked = on; box.closest('tr').classList.toggle('selected', on); }
  document.querySelectorAll('.cov-point').forEach((c) => c.setAttribute('stroke', state.selected.has(c.dataset.id) ? '#ffffff' : 'none'));
  updateSelection();
}
function updateSelection(busy = !!activeJob) {
  const ids = [...state.selected]; const dets = new Set(ids.map((i) => i.slice(-1)));
  $('selection-panel').hidden = !ids.length || state.stage !== 'epochs';
  $('selection-text').textContent = dets.size > 1 ? `${ids.length} selected on different detectors, which cannot be compared` : `${ids.length} selected on D${[...dets][0]}`;
  $('compare-selected').disabled = busy || ids.length !== 2 || dets.size !== 1;
  $('series-selected').disabled = busy || ids.length < 2 || ids.length > 8 || dets.size !== 1;
}
$('compare-selected').addEventListener('click', () => { const [a, b] = [...state.selected]; startCompare(a, b); });
$('series-selected').addEventListener('click', () => startSeries([...state.selected]));
$('clear-selected').addEventListener('click', () => { for (const id of [...state.selected]) toggleSelect(id, false); });

const backToEpochs = { label: 'Back to epochs', go: () => { showStage('epochs'); window.scrollTo({ top: 0 }); } };
async function startCompare(a, b) {
  const s = state.screening; const ord = (id) => (s.observations.find((o) => o.observation_id === id) || {}).mjd || 0;
  const [first, second] = ord(a) <= ord(b) ? [a, b] : [b, a];
  state.retry = () => startCompare(a, b);
  const out = await runJob('compare', { ...fieldBody(state.screenedField), observation_a: first, observation_b: second });
  if (out) { showStage('result'); renderComparison(out.job.result, out.base, $('result'), backToEpochs, liveSource(out.job)); }
}
async function startSeries(ids) {
  state.retry = () => startSeries(ids);
  const out = await runJob('series', { ...fieldBody(state.screenedField), observation_ids: ids });
  if (out) { showStage('result'); renderSeries(out.job.result, out.base, $('result'), backToEpochs, liveSource(out.job)); }
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

async function renderComparison(report, base, host, back, source) {
  clearHost(host); host.hidden = false;
  const token = host.token = {};
  const r = report; const ea = r.epochs.a, eb = r.epochs.b; const range = r.spectral.wavelength_range_um;
  const cands = r.candidates || [], vs = r.variability_summary || {};
  if (back) host.append(backBar(back));
  host.append(el('header', { class: 'verdict' },
    el('div', { class: 'verdict-top' }, el('span', { class: 'badge big ' + (STATUS_BADGE[r.status] || 'neutral') }, r.status), el('p', { class: 'verdict-reason' }, r.reason || ''), sourceChip(source)),
    el('h2', {}, el('span', {}, short(r.observation_a)), el('span', { class: 'arrow' }, '→'), el('span', {}, short(r.observation_b))),
    el('p', { class: 'sub' }, `${ea.date.slice(0, 10)} to ${eb.date.slice(0, 10)} · ${fmt(r.separation_days, 1)} days apart · detector D${ea.detector} · matched λ ${range ? fmt(range[0], 4) + '–' + fmt(range[1], 4) + ' µm' : 'none'}`),
    el('div', { class: 'gates' }, Object.entries(r.gates || {}).map(([k, v]) => el('span', { class: 'gate ' + (STATUS_BADGE[v] || 'neutral') }, GATE_TEXT[k] || k, el('em', {}, v))))));
  const target = targetStrip(r); if (target) host.append(target);
  const viewerHost = el('div', {}, el('div', { class: 'empty' }, 'Loading image data…'));
  const residuals = el('div', { class: 'tabpanel' });
  const tabs = tabset([
    { id: 'images', label: 'Images', node: el('div', { class: 'tabpanel' }, el('p', { class: 'intro' }, 'Background-subtracted, PSF-matched images on epoch A’s pixel grid. Dimmed pixels sampled different wavelengths in the two epochs and are not compared.'), viewerHost) },
    { id: 'validation', label: 'Validation', node: validationPanel(r) },
    { id: 'residuals', label: 'Residuals', count: cands.length, flag: cands.some((c) => INTERESTING.has(c.classification)), node: residuals },
    { id: 'stars', label: 'Star fluxes', count: vs.beyond_threshold || 0, flag: (vs.passed_vetting || 0) > 0, node: variabilityPanel(r) },
    { id: 'details', label: 'Provenance and files', node: el('div', { class: 'tabpanel' }, provenancePanel(r, source),
      el('h3', {}, 'Limitations'), el('ul', { class: 'plain' }, (r.limitations || []).map((l) => el('li', {}, l))),
      downloadsPanel(base, [['comparison.fits', 'FITS products'], ['report.json', 'Full report (JSON)'], ['candidates.csv', 'Residual table (CSV)'], ['comparison.png', 'Summary figure (PNG)']])) },
  ]);
  host.append(tabs.bar, tabs.body);
  host.scrollIntoView({ block: 'start' });
  let viewer = null;
  try {
    const data = await api(base + 'viewer.json');
    if (host.token !== token) return;
    viewerHost.replaceChildren();
    viewer = host.viewer = new ImageViewer(viewerHost, data, r);
  } catch (e) {
    if (host.token !== token) return;
    viewerHost.replaceChildren(el('div', { class: 'notice error' }, 'Image data unavailable: ' + e.message));
  }
  candidatesPanel(residuals, r, viewer && ((c) => { tabs.show('images'); viewer.focus(c.x, c.y); }));
}

function targetStrip(r) {
  const t = r.target; if (!t) return null;
  if (t.status !== 'MEASURED') return el('div', { class: 'notice' }, el('b', {}, human(t.status) + '. '), t.reason || '');
  return el('div', { class: 'target' }, el('h4', {}, 'Target'), el('dl', { class: 'facts' },
    el('div', {}, el('dt', {}, 'Flux change (B/A − 1)'), el('dd', {}, `${fmt(t.fractional_change, 4)} ± ${fmt(t.fractional_change_sigma, 4)}`)),
    el('div', {}, el('dt', {}, 'Significance'), el('dd', {}, `${fmt(t.change_z, 2)} σ`)),
    el('div', {}, el('dt', {}, 'S/N in A'), el('dd', {}, fmt(t.snr_a, 1)))));
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
  return el('div', { class: 'tabpanel' }, el('p', { class: 'intro' }, 'Each gate is computed from this pair’s own data. UNVERIFIED means the data could not test it, not that it passed.'), cards,
    el('p', { class: 'footnote' }, `PSF: ${r.psf.method}; zones ${r.psf.zone_a}/${r.psf.zone_b}; interpolation smoothing ${r.psf.interpolation_smoothing_variance_px2 ? r.psf.interpolation_smoothing_variance_px2.map((v) => fmt(v, 3)).join('/') : '—'} px² matched on A. Difference pedestal ${bg.applied ? fmt(bg.median_mjy_sr, 4) + ' MJy/sr removed' : 'not fitted'}. Catalog: ${(r.catalog && r.catalog.mode) || 'none'}. Calibration sources: ${r.calibration_sources || '—'}.`));
}
function candidatesPanel(host, r, focus) {
  const cands = r.candidates || []; const follow = cands.filter((c) => INTERESTING.has(c.classification));
  host.append(el('p', { class: 'intro' }, `Connected groups of at least ${r.policy.detection_min_pixels} pixels at ≥${r.policy.detection_sigma}σ in the difference image, classified against the propagated Gaia catalog and the calibrated noise model.${focus ? ' Click a row to centre it in the image.' : ''}`));
  const summary = Object.entries(r.candidate_summary || {}).sort((a, b) => b[1] - a[1]);
  if (!cands.length) { host.append(el('div', { class: 'empty' }, 'No residual components above threshold.')); return; }
  host.append(el('div', { class: 'badge-row' }, summary.map(([k, n]) => el('span', { class: 'badge ' + (INTERESTING.has(k) ? 'warning' : 'neutral') }, `${CLASS_TEXT[k] || k}: ${n}`))));
  const seg = el('div', { class: 'seg', role: 'group', 'aria-label': 'Residual filter' }), wrap = el('div', {});
  const draw = (mode) => {
    [...seg.children].forEach((b) => b.setAttribute('aria-pressed', String(b.dataset.mode === mode)));
    const rows = mode === 'all' ? cands : follow;
    wrap.replaceChildren(rows.length ? tableOf(['#', 'Classification', 'Peak σ', 'RA', 'Dec', 'Gaia G', 'ΔF/F', 'Reason'], [0, 2, 3, 4, 5, 6], rows.map((c) => el('tr', { class: focus ? 'clickable' : null, onclick: focus ? () => focus(c) : null },
      el('td', { class: 'num' }, c.id), el('td', {}, CLASS_TEXT[c.classification] || c.classification), el('td', { class: 'num' }, fmt(c.peak_sigma, 1)),
      el('td', { class: 'num' }, fmt(c.ra_deg, 5)), el('td', { class: 'num' }, fmt(c.dec_deg, 5)), el('td', { class: 'num' }, fmt(c.catalog_g_mag, 2)),
      el('td', { class: 'num' }, c.fractional_change === undefined ? '—' : `${fmt(c.fractional_change, 3)} ± ${fmt(c.fractional_change_sigma, 3)}`), el('td', { class: 'small secondary' }, c.reason || ''))))
      : el('div', { class: 'empty' }, 'No residual needs follow-up. Every component is explained by known motion, a PSF or registration dipole, the calibrated systematic floor, or a persistent unchanged source.'));
  };
  seg.append(el('button', { type: 'button', 'data-mode': 'follow', onclick: () => draw('follow') }, `Needs follow-up · ${follow.length}`), el('button', { type: 'button', 'data-mode': 'all', onclick: () => draw('all') }, `All · ${cands.length}`));
  host.append(seg, wrap); draw('follow');
}
function variabilityPanel(r) {
  const v = r.variability_candidates || []; const s = r.variability_summary || {};
  const p = el('div', { class: 'tabpanel' }, el('p', { class: 'intro' }, `${s.stars_tested || 0} Gaia stars measured in matched 5×5 apertures in both epochs. ${s.beyond_threshold || 0} exceed ${r.policy.detection_sigma}σ under the calibrated noise model and ${s.passed_vetting || 0} pass vetting (clear of flagged pixels, away from the region edge, no catalog blend, PSF-shaped residual).`));
  if (!v.length) { p.append(el('div', { class: 'empty' }, 'No star changed beyond the calibrated noise model.')); return p; }
  p.append(tableOf(['Gaia DR3 source', 'G', 'ΔF/F', 'z', 'Verdict', 'Failed checks'], [1, 2, 3], v.map((c) => {
    const pass = c.classification === 'SOURCE_LEVEL_VARIABILITY_CANDIDATE';
    return el('tr', {}, el('td', { class: 'mono' }, c.source_id), el('td', { class: 'num' }, fmt(c.g_mag, 2)), el('td', { class: 'num' }, `${fmt(c.fractional_change, 3)} ± ${fmt(c.fractional_change_sigma, 3)}`),
      el('td', { class: 'num' }, fmt(c.change_z, 1)), el('td', {}, el('span', { class: 'badge ' + (pass ? 'warning' : 'neutral') }, pass ? 'Candidate' : 'Rejected')),
      el('td', { class: 'small secondary' }, (c.vetting && c.vetting.failed.length) ? c.vetting.failed.join(', ').replaceAll('_', ' ') : '—'));
  })));
  return p;
}
/* Where a result came from: a saved run made by scripts/make_demos.py, or a job this server just ran. */
const recordedSource = (run) => ({ kind: 'recorded', label: 'Recorded run', detail: `generated ${(run.generated_utc || '').replace('T', ' ').slice(0, 16)} UTC by scripts/make_demos.py`, reproduce: run.reproduce });
const liveSource = (job) => ({ kind: 'live', label: 'Live run', detail: `job ${String(job.id || '').slice(0, 8)}${job.finished_utc ? ` · finished ${job.finished_utc.replace('T', ' ').slice(0, 16)} UTC` : ''}` });
function sourceChip(source) { return source ? el('span', { class: 'source ' + source.kind, title: source.detail }, source.label) : null; }
const hash = (v) => (v ? el('span', { class: 'mono', title: v }, v.slice(0, 12) + '…') : '—');
const fileName = (url) => (url ? String(url).split('?')[0].split('/').pop() : '—');

function provenancePanel(report, source) {
  const r = report; const cat = r.catalog || {};
  const epochs = r.epochs ? (Array.isArray(r.epochs) ? r.epochs : [r.epochs.a, r.epochs.b]) : [];
  const files = r.science_files || {};
  const rows = epochs.map((e) => {
    const cal = e.calibration || {};
    return el('tr', {}, el('td', { class: 'mono' }, String(e.date || '').slice(0, 10)), el('td', { class: 'mono' }, short(e.observation_id)), el('td', { class: 'num' }, 'D' + e.detector),
      el('td', { class: 'mono small' }, fileName(cal.calibration_url || e.calibration_url)), el('td', {}, hash(cal.calibration_sha256)), el('td', {}, hash((files[e.observation_id] || {}).sha256)));
  });
  for (const u of r.unavailable_observations || []) {
    rows.push(el('tr', {}, el('td', {}, '—'), el('td', { class: 'mono' }, short(u.observation_id)), el('td', { class: 'num' }, 'D' + String(u.observation_id).slice(-1)),
      el('td', { class: 'small secondary', colspan: '3' }, `Not used: ${u.reason || u.status}`)));
  }
  const facts = [
    ['Source', source ? `${source.label} · ${source.detail}` : 'Not recorded'],
    source && source.reproduce ? ['Reproduce', el('code', {}, source.reproduce)] : null,
    ['Processing code', r.software_sha256 ? el('span', { class: 'mono' }, `SHA-256 ${r.software_sha256.slice(0, 16)}…`) : '—'],
    ['Star catalog', cat.mode ? `${cat.catalog || 'Gaia DR3'} · ${cat.mode === 'UNAVAILABLE' ? 'unavailable, so image-detected sources were used' : cat.mode}${cat.retrieved_utc ? ` · retrieved ${cat.retrieved_utc.slice(0, 10)}` : ''}` : 'none'],
    r.cutout ? ['Cutout', `${fmt(r.cutout.size_deg, 2)}° centred on RA ${fmt(r.cutout.center_deg[0], 5)}°, Dec ${fmt(r.cutout.center_deg[1], 5)}°`] : null,
  ].filter(Boolean);
  return [el('h3', {}, 'Provenance'),
    el('dl', { class: 'kv wide' }, facts.flatMap(([k, v]) => [el('dt', {}, k), el('dd', {}, v)])),
    el('div', { style: 'margin-top:14px' }, tableOf(['Date', 'Observation', 'Det.', 'Calibration (own, via DataLink)', 'Calibration SHA-256', 'Science file SHA-256'], [2], rows)),
    el('p', { class: 'footnote' }, 'Each observation is calibrated with the wavelength map IRSA DataLink associates with it; no other epoch’s calibration is substituted. Direct equality of science and calibration UUIDs is not verified.')];
}
function downloadsPanel(base, files) {
  return [el('h3', {}, 'Downloads'), el('div', { class: 'downloads' }, files.map(([f, label]) => el('a', { class: 'btn small', href: base + f, download: f }, label)))];
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
    this.sliderWrap = el('label', { style: 'display:none' }, 'A', this.slider, 'B');
    const marks = el('label', { class: 'push' }, el('input', { type: 'checkbox', checked: true, onchange: (e) => { this.showMarks = e.target.checked; this.render(); } }), 'Residual markers');
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
function renderSeries(r, base, host, back, source) {
  clearHost(host); host.hidden = false;
  const ph = r.photometry, as = r.astrometry; const pts = ph.points || [];
  const phBadge = { CONSISTENT_WITH_CONSTANT: 'good', VARIABILITY_CANDIDATE: 'warning' }[ph.status] || 'neutral';
  const asBadge = { MOTION_MEASURED: 'good', NO_SIGNIFICANT_MOTION: 'neutral' }[as.status] || 'neutral';
  const wave = (pts.find((p) => Number.isFinite(p.wavelength_um)) || {}).wavelength_um;
  if (back) host.append(backBar(back));
  host.append(el('header', { class: 'verdict' },
    el('div', { class: 'verdict-top' }, el('span', { class: 'badge big ' + phBadge }, 'Light curve: ' + human(ph.status)), el('span', { class: 'badge big ' + asBadge }, 'Motion: ' + human(as.status)), sourceChip(source)),
    el('h2', {}, `Time series of ${r.epoch_count} epochs`),
    el('p', { class: 'sub' }, `Detector D${r.detector}${wave ? ` · ${fmt(wave, 3)} µm` : ''} · measured-ePSF fits on native pixels · brightness and motion are analysed separately`)));
  if (r.analysis_status === 'PARTIAL') {
    host.append(el('div', { class: 'notice info' }, el('p', {}, el('b', {}, `Partial result: ${r.epoch_count} of ${r.requested_epoch_count} observations calibrated. `), 'Failed epochs are listed here and were not assigned a substitute calibration.'),
      el('ul', { class: 'plain', style: 'margin-top:10px' }, (r.unavailable_observations || []).map((o) => el('li', {}, el('span', { class: 'mono' }, short(o.observation_id)), ` · ${o.stage}: ${o.reason}`)))));
  }
  const measured = (r.epochs || []).filter((e) => e.target && e.target.status === 'MEASURED').length;
  if (!pts.length && !(as.positions || []).length) host.append(el('div', { class: 'notice' }, el('p', {}, el('b', {}, 'Nothing to plot for this target. '),
    `The target was measured in ${measured} of ${(r.epochs || []).length} epochs, so there is no light curve or trajectory to draw. The reason for each epoch is listed in the light-curve section.`)));
  host.append(seriesAssistant(r, base));

  const lc = el('section', { class: 'panel' }, el('div', { class: 'panel-head' }, el('div', {}, el('h3', {}, 'Same-wavelength light curve'), el('p', {}, ph.interpretation || ph.reason || ''))));
  host.append(lc);
  if (pts.some((p) => p.status === 'USED')) { const c = el('div', { class: 'chart' }); lc.append(c); lightCurve(c, pts, ph); }
  if (!pts.length) lc.append(el('div', { class: 'empty' }, 'No light curve to draw: no epoch has a usable measurement of the target.'), epochTable(r));
  else lc.append(el('div', { style: 'margin-top:14px' }, tableOf(['Date', 'Observation', 'Status', 'λ (µm)', 'Flux (mJy)', 'σ (mJy)', 'Zero point', 'Note'], [3, 4, 5, 6], pts.map((p) => el('tr', {},
    el('td', { class: 'mono' }, p.date.slice(0, 10)), el('td', { class: 'mono' }, short(p.observation_id)), el('td', {}, el('span', { class: 'badge ' + (p.status === 'USED' ? 'good' : 'neutral') }, human(p.status))),
    el('td', { class: 'num' }, fmt(p.wavelength_um, 4)), el('td', { class: 'num' }, p.calibrated_flux_jy ? fmt(p.calibrated_flux_jy * 1e3, 3) : '—'), el('td', { class: 'num' }, p.flux_sigma_jy ? fmt(p.flux_sigma_jy * 1e3, 3) : '—'),
    el('td', { class: 'num' }, fmt(p.zero_point, 4)), el('td', { class: 'small secondary' }, p.reason || ''))))));
  const test = ph.constant_flux_test;
  if (test) lc.append(el('p', { class: 'footnote' }, `Constant-flux test: χ² = ${fmt(test.chi2_constant, 2)} for ${test.dof} dof (p = ${fmt(test.pvalue, 3)}); weighted mean ${fmt(test.constant_flux * 1e3, 3)} mJy. Field-star floor ${fmt(ph.fractional_floor, 3)}.`));

  const ap = el('section', { class: 'panel' }, el('div', { class: 'panel-head' }, el('div', {}, el('h3', {}, 'Astrometric trajectory'), el('p', {}, as.reason || 'Target positions corrected by the per-epoch frame offset of image-detected field stars. Error bars include the empirical field-star floor.'))));
  host.append(ap);
  if ((as.positions || []).length) { const c = el('div', { class: 'chart' }); ap.append(c); trajectory(c, as); }
  else ap.append(el('div', { class: 'empty' }, 'No trajectory to draw: the target position was not measured in enough epochs.'));
  const sol = as.astrometric_solution;
  if (sol && sol.status === 'SOLVED') ap.append(el('div', { class: 'check-grid' },
    check('Proper motion (measured)', null, `${fmt(sol.pm_east_arcsec_per_year, 2)}, ${fmt(sol.pm_north_arcsec_per_year, 2)} ″/yr`, `± ${fmt(sol.pm_sigma_arcsec_per_year[0], 2)}, ${fmt(sol.pm_sigma_arcsec_per_year[1], 2)} (Birge ratio ${fmt(sol.birge_ratio_applied, 2)} applied); detected at ${fmt(sol.pm_significance, 1)}σ.`),
    check('Catalog proper motion', null, sol.catalog_pm_arcsec_per_year ? `${fmt(sol.catalog_pm_arcsec_per_year[0], 2)}, ${fmt(sol.catalog_pm_arcsec_per_year[1], 2)} ″/yr` : '—', `Difference ${fmt(sol.pm_difference_from_catalog_sigma, 2)}σ. SIMBAD values; a known-mover test, not a discovery.`),
    check('Parallax', null, `${fmt(sol.parallax_arcsec, 3)} ± ${fmt(sol.parallax_sigma_arcsec, 3)}″`, `Catalog ${fmt(sol.catalog_parallax_arcsec, 3)}″ (${fmt(sol.parallax_difference_from_catalog_sigma, 1)}σ apart). Observer: JPL Horizons SPHEREx vectors.`),
    check('Astrometric floor', null, `${fmt(as.astrometric_floor_arcsec, 3)}″`, `Epoch-to-epoch repeatability of ${as.field_stars} field stars beyond fit errors.`)));
  host.append(el('section', { class: 'panel tabpanel' }, provenancePanel(r, source),
    el('h3', {}, 'Limitations'), el('ul', { class: 'plain' }, (r.limitations || []).map((l) => el('li', {}, l))),
    downloadsPanel(base, [['series.json', 'Full report (JSON)'], ['light_curve.csv', 'Light curve (CSV)'], ['series.png', 'Summary figure (PNG)']])));
  host.scrollIntoView({ block: 'start' });
}

function epochTable(r) {
  return tableOf(['Date', 'Observation', 'Target', 'Reason'], [], (r.epochs || []).map((e) => {
    const t = e.target || {};
    return el('tr', {}, el('td', { class: 'mono' }, (e.date || '').slice(0, 10)), el('td', { class: 'mono' }, short(e.observation_id)),
      el('td', {}, el('span', { class: 'badge ' + (t.status === 'MEASURED' ? 'good' : 'neutral') }, human(t.status || 'unknown'))), el('td', { class: 'small secondary' }, t.reason || ''));
  }));
}

function seriesAssistant(report, base) {
  const input = el('input', { type: 'text', maxlength: '240', placeholder: 'Ask about flux, motion, uncertainty or excluded epochs', 'aria-label': 'Ask about this time series' });
  const ask = el('button', { class: 'btn primary', type: 'button' }, 'Ask');
  const response = el('div', { class: 'notice info answer', 'aria-live': 'polite', hidden: true });
  const jobId = (base.match(/\/jobs\/([^/]+)\/files\//) || [])[1];
  const submit = async (question = input.value) => {
    response.hidden = false;
    if (!question.trim()) { response.textContent = 'Type a question first.'; return; }
    input.value = question;
    ask.disabled = true; response.textContent = 'Asking the local model…';
    try {
      if (!jobId) throw new Error('This result has no live job record.');
      const result = await api(`/api/lab/jobs/${jobId}/assistant`, {
        method: 'POST', headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ question }),
      });
      response.textContent = `${result.answer}\n\nSource: ${result.model} · ${result.source}`;
    } catch (error) {
      response.textContent = `Local LLM unavailable; report-based fallback: ${explainSeries(report, question)} (${error.message})`;
    } finally { ask.disabled = false; }
  };
  const quick = el('div', { class: 'assistant-quick' },
    ...['Summarize this run', 'Did the brightness change?', 'What was excluded?', 'What does the motion show?'].map((q) => el('button', { class: 'btn small', type: 'button', onclick: () => submit(q) }, q)));
  const submitInput = () => submit(input.value);
  ask.addEventListener('click', submitInput);
  input.addEventListener('keydown', (event) => { if (event.key === 'Enter') submitInput(); });
  return el('section', { class: 'panel ask' },
    el('div', { class: 'panel-head' }, el('div', {}, el('h3', {}, 'Ask about this run'), el('p', {}, 'A small local model explains this run’s measured report. It cannot fill missing epochs or make new detections.'))),
    el('div', { class: 'row' }, input, ask), quick, response);
}

function explainSeries(report, question) {
  const q = String(question || '').toLowerCase();
  const phot = report.photometry || {}, ast = report.astrometry || {};
  const used = (phot.points || []).filter((p) => p.status === 'USED');
  const excluded = (phot.points || []).filter((p) => p.status !== 'USED');
  const missing = report.unavailable_observations || [];
  const test = phot.constant_flux_test;
  if (!q.trim()) return 'Type a question about the results first.';
  if (/excluded|missing|fail|unavailable|skip/.test(q)) {
    const rows = [...missing.map((x) => `${short(x.observation_id)} unavailable during ${x.stage}: ${x.reason}`),
      ...excluded.map((x) => `${short(x.observation_id)} excluded: ${x.reason || x.status}`)];
    return rows.length ? rows.join(' ') : 'No unavailable or excluded epochs were recorded in this report.';
  }
  if (/motion|moving|position|astrom/.test(q)) {
    const fit = ast.astrometric_solution;
    if (fit && fit.status === 'SOLVED') return `The measured motion is ${fmt(fit.pm_east_arcsec_per_year, 2)} east and ${fmt(fit.pm_north_arcsec_per_year, 2)} north arcsec/year. The fit significance is ${fmt(fit.pm_significance, 1)}σ. Catalog motion differs by ${fmt(fit.pm_difference_from_catalog_sigma, 2)}σ. This is a known-target consistency check, not a new-object discovery.`;
    return `Astrometry status: ${ast.status || 'not available'}. ${ast.reason || 'The report has no solved motion fit.'}`;
  }
  if (/uncertain|error bar|precision|sigma|σ/.test(q)) {
    const values = used.map((p) => `${short(p.observation_id)}: ${fmt(p.flux_sigma_jy * 1e3, 3)} mJy (total reported uncertainty)`);
    if (values.length) return `Reported flux uncertainties for used epochs: ${values.join('; ')}. They include the statistical term and available zero-point/floor terms in the pipeline report.`;
    return 'No used photometric points have uncertainty estimates in this report.';
  }
  if (/bright|flux|variab|constant|change/.test(q)) {
    if (test) return `${phot.interpretation || phot.status}. ${used.length} epochs were used; weighted mean flux is ${fmt(test.constant_flux * 1e3, 3)} mJy. Constant-flux test: χ²=${fmt(test.chi2_constant, 2)} for ${test.dof} degrees of freedom, p=${fmt(test.pvalue, 3)}. This test applies only to the wavelength-matched, calibrated points.`;
    return `Photometry status: ${phot.status || 'not available'}. ${phot.reason || `Only ${used.length} usable same-wavelength epochs are available; no constant-flux test was produced.`}`;
  }
  if (/summar|result|overall|what happened|explain/.test(q)) {
    const coverage = `${report.epoch_count} of ${report.requested_epoch_count ?? report.epoch_count} requested epochs were calibrated`;
    const variability = test ? ` The ${used.length} usable points are ${phot.status === 'CONSISTENT_WITH_CONSTANT' ? 'consistent with constant flux' : phot.status === 'VARIABILITY_CANDIDATE' ? 'flagged as a variability candidate' : phot.status}; p=${fmt(test.pvalue, 3)}.` : ` Photometry: ${phot.status || 'not available'}.`;
    const motion = ast.astrometric_solution?.status === 'SOLVED' ? ` Motion fit: ${fmt(ast.astrometric_solution.pm_significance, 1)}σ.` : ` Astrometry: ${ast.status || 'not available'}.`;
    return `${report.analysis_status === 'PARTIAL' ? 'Partial run: ' : ''}${coverage}.${variability}${motion} A variability candidate is not a discovery; review the limitations and unavailable-epoch list.`;
  }
  return 'I can answer from this report about the light curve/flux, astrometric motion, uncertainty, and excluded or unavailable epochs. Ask one of those directly; I do not query new data or infer missing measurements.';
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
const KIND = { compare: 'Comparison', series: 'Time series', screen: 'Screening' };
const perYear = (v) => (v ? v.map((x) => fmt(x, 2)).join(', ') + ' ″/yr' : '—');
let runsRequest = null;
function fetchRuns() { if (!runsRequest) runsRequest = api('/api/lab/demos').catch((e) => { runsRequest = null; throw e; }); return runsRequest; }
function runCard(run) {
  const ok = run.status === 'COMPLETED'; const s = run.summary || {};
  let badge = ['critical', 'Failed'], facts = [['Error', run.error]];
  if (ok && run.kind === 'compare') {
    badge = [STATUS_BADGE[s.status] || 'neutral', s.status];
    facts = [['Matched pixels', count(s.matched_pixels)], ['Baseline', fmt(s.separation_days, 0) + ' days'], ['Matched λ', s.wavelength_range_um ? `${fmt(s.wavelength_range_um[0], 2)}–${fmt(s.wavelength_range_um[1], 2)} µm` : '—'], ['Target', human(s.target_status)]];
  } else if (ok && run.kind === 'series') {
    badge = ['neutral', `${(run.observations || []).length} epochs`];
    facts = [['Light curve', human(s.photometry_status)], ['Motion', human(s.astrometry_status)], ['Measured motion', perYear(s.proper_motion_arcsec_per_year)], ['Catalog motion', perYear(s.catalog_proper_motion)]];
  } else if (ok && run.kind === 'screen') {
    badge = ['neutral', `${Object.keys(s).length} targets`];
    facts = Object.entries(s).map(([t, v]) => [targetName(t), `${v.matched_pairs} pairs · ${v.groups_3plus} groups`]);
  } else if (ok) { badge = ['good', 'Completed']; facts = []; }
  const catalog = run.catalog_mode === 'UNAVAILABLE' ? 'Gaia catalog unavailable' : run.catalog_mode;
  return el(ok ? 'a' : 'article', { class: 'card', href: ok ? '#runs/' + run.name : null },
    el('div', { class: 'card-top' }, el('p', { class: 'eyebrow' }, KIND[run.kind] || 'Run'), el('span', { class: 'badge ' + badge[0] }, badge[1])),
    el('h3', {}, microns(run.title || run.name)), run.why ? el('p', {}, run.why) : null,
    el('dl', { class: 'kv' }, facts.flatMap(([k, v]) => [el('dt', {}, k), el('dd', {}, v || '—')])),
    el('div', { class: 'card-foot' }, el('span', {}, [(run.generated_utc || '').slice(0, 10), catalog].filter(Boolean).join(' · ')), ok ? el('b', {}, 'Open →') : null));
}
async function showRuns(open) {
  const index = $('runs-index'), result = $('run-result'), cards = $('run-cards');
  let runs;
  try { await targetsReady; runs = await fetchRuns(); } catch (e) {
    index.hidden = false; result.hidden = true;
    cards.replaceChildren(el('div', { class: 'notice error' }, 'Could not load recorded runs: ' + e.message)); return;
  }
  if (!cards.querySelector('.card')) cards.replaceChildren(...(runs.length ? runs.map(runCard) : [el('div', { class: 'empty' }, 'No recorded runs yet. Generate them with: python -m scripts.make_demos')]));
  const run = open && runs.find((r) => r.name === open && r.status === 'COMPLETED');
  index.hidden = !!run; result.hidden = !run;
  window.scrollTo({ top: 0 });
  if (run) await openRun(run, result); else clearHost(result);
}
async function openRun(run, host) {
  const base = demoBase(run.name), source = recordedSource(run); const back = { label: 'All recorded runs', go: () => go('runs') };
  clearHost(host); host.append(el('div', { class: 'empty' }, 'Loading…'));
  try {
    if (run.kind === 'compare') await renderComparison(await api(base + 'report.json'), base, host, back, source);
    else if (run.kind === 'series') renderSeries(await api(base + 'series.json'), base, host, back, source);
    else {
      const pick = el('div', { class: 'seg', role: 'group', 'aria-label': 'Target' }), inner = el('div', { class: 'content' });
      const show = async (file, button) => {
        [...pick.children].forEach((b) => b.setAttribute('aria-pressed', String(b === button)));
        inner.replaceChildren(el('div', { class: 'empty' }, 'Loading…'));
        try { renderScreening(await api(base + file), inner, false); } catch (e) { inner.replaceChildren(el('div', { class: 'notice error' }, 'Could not load this screening: ' + e.message)); }
      };
      for (const file of (run.files || []).filter((f) => f.startsWith('screening_'))) pick.append(el('button', { type: 'button', 'aria-pressed': 'false', onclick: (e) => show(file, e.currentTarget) }, targetName(file.replace('screening_', '').replace('.json', ''))));
      clearHost(host); host.append(backBar(back), el('p', { class: 'sub' }, microns(run.title)), pick, inner);
      if (pick.firstChild) pick.firstChild.click();
    }
  } catch (e) { clearHost(host); host.append(backBar(back), el('div', { class: 'notice error' }, 'Could not open run: ' + e.message)); }
}
async function loadQuickRuns() {
  try {
    const runs = (await fetchRuns()).filter((r) => r.status === 'COMPLETED' && r.title);
    $('quick-runs').replaceChildren(...runs.map((run) => el('a', { class: 'quick', href: '#runs/' + run.name }, el('span', { class: 'kind' }, KIND[run.kind] || 'Run'), el('span', {}, microns(run.title)), el('span', { class: 'go' }, '→'))));
    $('shortcuts').hidden = !runs.length;
  } catch (_) { $('shortcuts').hidden = true; }
}

/* ------------------------------------------------------------------ method */
function renderMethod() {
  const m = $('method-body'); if (m.childElementCount) return;
  m.innerHTML = `
  <p class="eyebrow">Method</p>
  <h1>Compare like with like</h1>
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
  <p>A residual is not a discovery. Moving-object identification requires a measured multi-epoch trajectory; catalog-forced positions are not detections; different-wavelength epochs are rejected rather than corrected with an invented spectral model. Saturated or flagged targets are reported as unmeasurable. When an external service (IRSA, Gaia, JPL) fails, the job fails with that reason. Nothing is substituted.</p>
  <h2>Policies</h2>
  <p>Spectral tolerance 0.10 × CBAND; detection 5σ with ≥2 connected pixels; registration tolerance 0.3 px; aperture 5×5 native pixels; null test needs ≥200 source-free pixels; source noise model needs ≥30 stars. These are analysis choices recorded in every report, not mission calibration constants.</p>
  <h2>Data</h2>
  <p>NASA/IPAC IRSA SPHEREx Data Release 1 (SIA search, cutout service, DataLink calibration products and ePSF tables), ESA Gaia DR3, SIMBAD/CDS Sesame and JPL Horizons. See the repository's <code>docs/</code> for full methodology and reproduction instructions.</p>`;
}

/* ------------------------------------------------------------------ boot */
if (STATIC) {
  // No server behind a static copy: show the note instead of the search, and drop the API docs link.
  $('static-note').hidden = false;
  document.querySelector('.finder').hidden = true;
  document.querySelector('.toplink').hidden = true;
}
targetsReady = loadTargets();
loadQuickRuns();
setField(null);
window.addEventListener('hashchange', route);
route();
