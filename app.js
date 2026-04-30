// === CONSTANTS & CONFIG ===

const POLAR_URL =
  'https://raw.githubusercontent.com/XCSoar/XCSoar/master/src/Polar/PolarStore.cpp';

const MARGIN = { top: 30, right: 24, bottom: 56, left: 70 };

// Conversion factors from internal units (km/h, m/s) to display units
const SPEED_UNITS = {
  kts: { label: 'kts',  factor: 1 / 1.852 },
  kmh: { label: 'km/h', factor: 1 },
  mph: { label: 'mph',  factor: 1 / 1.60934 },
};

const SINK_UNITS = {
  kts: { label: 'kts',  factor: 1 / 0.51444 },
  fts: { label: 'ft/s', factor: 1 / 0.3048 },
  ms:  { label: 'm/s',  factor: 1 },
};

// Resolved once from CSS variables — canvas can't use CSS vars directly
const C = {
  polar:  '#1976d2',
  mc:     '#e53935',
  hover:  '#2e7d32',
  grid:   '#e0e0e0',
  axis:   '#424242',
  text:   '#212121',
  minsink: '#9c27b0',
};

// === APPLICATION STATE ===

const state = {
  polars: [],
  selectedIndex: 0,
  coeffs: null,       // { a, b, c } — quadratic fit in km/h and m/s
  airmass_ms: 0,      // m/s, positive = lift
  mc_ms: 0,           // m/s, MacCready setting (positive)
  speedUnit: 'kts',
  sinkUnit: 'kts',
  hoverPoint: null,   // { v_kmh, w_ms } | null
  ranges: null,       // cached axis ranges
  canvasW: 0,
  canvasH: 0,
};

// === UNIT CONVERSION ===

function convertSpeed(v_kmh, unit) {
  return v_kmh * SPEED_UNITS[unit].factor;
}

// w_ms is negative (sink); returns positive display value
function convertSink(w_ms, unit) {
  return Math.abs(w_ms) * SINK_UNITS[unit].factor;
}

function speedLabel() { return SPEED_UNITS[state.speedUnit].label; }
function sinkLabel()  { return SINK_UNITS[state.sinkUnit].label; }

function ktsToMs(kts) { return kts * 0.51444; }

// === DATA FETCHING & PARSING ===

async function fetchPolars() {
  const resp = await fetch(POLAR_URL);
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  const text = await resp.text();
  return parsePolars(text);
}

function parsePolars(cpp) {
  const results = [];

  // Match each { "Name" [optional /* comment */], num, num, ...(12 numbers) }
  const ENTRY_RE =
    /\{\s*"([^"]+)"\s*(?:\/\*[^*]*\*+(?:[^/*][^*]*\*+)*\/)?\s*(?:,\s*(-?[\d.]+)){12}\s*\}/g;

  // The above won't capture each number separately in all JS engines due to repeated group capture.
  // Use a two-step approach: first find the whole entry, then extract numbers.
  const BLOCK_RE =
    /\{\s*"([^"]+)"\s*(?:\/\*[^*]*\*+(?:[^/*][^*]*\*+)*\/)?\s*,([^}]+)\}/g;

  let m;
  while ((m = BLOCK_RE.exec(cpp)) !== null) {
    const name = m[1].trim();
    const nums = m[2].match(/-?[\d]+\.?[\d]*/g);
    if (!nums || nums.length < 12) continue;

    const [empty_mass, max_ballast, v1, w1, v2, w2, v3, w3, wing_area, ref_mass, v_no, max_speed]
      = nums.map(Number);

    // Sanity: speeds must be positive and in realistic range, sink must be negative
    if (v1 <= 0 || v2 <= 0 || v3 <= 0) continue;
    if (w1 >= 0 || w2 >= 0 || w3 >= 0) continue;

    results.push({ name, empty_mass, max_ballast, v1, w1, v2, w2, v3, w3,
                   wing_area, ref_mass, v_no, max_speed });
  }

  results.sort((a, b) => a.name.localeCompare(b.name));

  // Deduplicate by name (keep first occurrence after sort)
  const seen = new Set();
  return results.filter(p => {
    if (seen.has(p.name)) return false;
    seen.add(p.name);
    return true;
  });
}

// === POLAR MATH ===

// Fit quadratic w = a*v^2 + b*v + c through 3 points using Cramer's rule
function fitPolar(entry) {
  const { v1, w1, v2, w2, v3, w3 } = entry;

  function det3(M) {
    return (
      M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
      M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
      M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    );
  }

  function replaceCol(rows, col, vec) {
    return rows.map((row, i) => row.map((v, j) => j === col ? vec[i] : v));
  }

  const A = [
    [v1 * v1, v1, 1],
    [v2 * v2, v2, 1],
    [v3 * v3, v3, 1],
  ];
  const rhs = [w1, w2, w3];
  const detA = det3(A);
  if (Math.abs(detA) < 1e-12) return null;

  const a = det3(replaceCol(A, 0, rhs)) / detA;
  const b = det3(replaceCol(A, 1, rhs)) / detA;
  const c = det3(replaceCol(A, 2, rhs)) / detA;

  if (a >= 0) return null; // degenerate — real polars have a < 0 (parabola opens downward)
  return { a, b, c };
}

function polarSink(coeffs, v_kmh) {
  return coeffs.a * v_kmh * v_kmh + coeffs.b * v_kmh + coeffs.c;
}

// MacCready optimal speed to fly.
// Tangent from (v=0, w=mc_ms) to the airmass-shifted polar w_s = w + airmass_ms.
// v_opt = sqrt((c + airmass_ms - mc_ms) / a)  [a < 0, numerator also negative for valid result]
function computeMcOptimal(coeffs, mc_ms, airmass_ms) {
  const discriminant = (coeffs.c + airmass_ms - mc_ms) / coeffs.a;
  if (discriminant <= 0) return null;

  const v_opt = Math.sqrt(discriminant);
  const w_opt = polarSink(coeffs, v_opt) + airmass_ms;

  return { v_opt, w_opt };
}

// Min-sink speed: dw/dv = 0 → v = -b/(2a)
function minSinkSpeed(coeffs) {
  return -coeffs.b / (2 * coeffs.a);
}

// === CHART GEOMETRY ===

function chartArea() {
  return {
    left:   MARGIN.left,
    right:  state.canvasW - MARGIN.right,
    top:    MARGIN.top,
    bottom: state.canvasH - MARGIN.bottom,
  };
}

function computeRanges(entry, coeffs) {
  const v_min_kmh = entry.v1 * 0.85;
  const v_max_kmh = Math.min(entry.max_speed, entry.v_no * 1.05);

  let maxSink_ms = 0;
  for (let v = v_min_kmh; v <= v_max_kmh; v += 2) {
    const w = polarSink(coeffs, v) + state.airmass_ms;
    if (w < maxSink_ms) maxSink_ms = w; // w negative, more negative = more sink
  }

  const sinkExtent = Math.abs(maxSink_ms) * 1.2;

  return {
    v_min_kmh,
    v_max_kmh,
    sink_max_ms: sinkExtent, // positive — the maximum magnitude of sink shown
    v_min_disp:  convertSpeed(v_min_kmh, state.speedUnit),
    v_max_disp:  convertSpeed(v_max_kmh, state.speedUnit),
    sink_min_disp: 0,
    sink_max_disp: sinkExtent * SINK_UNITS[state.sinkUnit].factor,
  };
}

function toCanvasX(sinkDisp, ranges) {
  const { left, right } = chartArea();
  return left + (sinkDisp - ranges.sink_min_disp) /
    (ranges.sink_max_disp - ranges.sink_min_disp) * (right - left);
}

function toCanvasY(speedDisp, ranges) {
  const { top, bottom } = chartArea();
  // speed increases upward → high speed = low pixel Y
  return bottom - (speedDisp - ranges.v_min_disp) /
    (ranges.v_max_disp - ranges.v_min_disp) * (bottom - top);
}

function fromCanvasY(py, ranges) {
  const { top, bottom } = chartArea();
  return ranges.v_min_disp + (bottom - py) /
    (bottom - top) * (ranges.v_max_disp - ranges.v_min_disp);
}

// === CANVAS SETUP ===

const canvas = document.getElementById('polar-canvas');
const ctx = canvas.getContext('2d');

function resizeCanvas() {
  const dpr = window.devicePixelRatio || 1;
  const rect = canvas.parentElement.getBoundingClientRect();
  state.canvasW = rect.width;
  state.canvasH = rect.height;
  canvas.width  = rect.width  * dpr;
  canvas.height = rect.height * dpr;
  canvas.style.width  = rect.width  + 'px';
  canvas.style.height = rect.height + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
}

function resolveColours() {
  const s = getComputedStyle(document.documentElement);
  const get = v => s.getPropertyValue(v).trim();
  const apply = (key, cssVar) => { const v = get(cssVar); if (v) C[key] = v; };
  apply('polar',  '--color-polar');
  apply('mc',     '--color-mc');
  apply('hover',  '--color-hover');
  apply('grid',   '--color-grid');
  apply('axis',   '--color-axis');
  apply('text',   '--color-text');
}

// === DRAWING ===

function niceTickInterval(range, targetTicks) {
  const rough = range / targetTicks;
  const mag = Math.pow(10, Math.floor(Math.log10(rough)));
  const normalised = rough / mag;
  let nice;
  if (normalised < 1.5) nice = 1;
  else if (normalised < 3.5) nice = 2;
  else if (normalised < 7.5) nice = 5;
  else nice = 10;
  return nice * mag;
}

function drawGrid(ranges) {
  const { left, right, top, bottom } = chartArea();
  ctx.save();
  ctx.strokeStyle = C.grid;
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 4]);
  ctx.font = '11px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.fillStyle = C.axis;
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';

  // Horizontal lines — speed axis (Y)
  const speedRange = ranges.v_max_disp - ranges.v_min_disp;
  const speedTick = niceTickInterval(speedRange, 6);
  const speedStart = Math.ceil(ranges.v_min_disp / speedTick) * speedTick;
  for (let v = speedStart; v <= ranges.v_max_disp + speedTick * 0.01; v += speedTick) {
    const py = toCanvasY(v, ranges);
    if (py < top - 1 || py > bottom + 1) continue;
    ctx.beginPath();
    ctx.moveTo(left, py);
    ctx.lineTo(right, py);
    ctx.stroke();
    ctx.fillText(v.toFixed(speedTick < 1 ? 1 : 0), left - 6, py);
  }

  // Vertical lines — sink axis (X)
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  const sinkRange = ranges.sink_max_disp - ranges.sink_min_disp;
  const sinkTick = niceTickInterval(sinkRange, 6);
  const sinkStart = Math.ceil(ranges.sink_min_disp / sinkTick) * sinkTick;
  for (let s = sinkStart; s <= ranges.sink_max_disp + sinkTick * 0.01; s += sinkTick) {
    const px = toCanvasX(s, ranges);
    if (px < left - 1 || px > right + 1) continue;
    ctx.beginPath();
    ctx.moveTo(px, top);
    ctx.lineTo(px, bottom);
    ctx.stroke();
    ctx.fillText(s.toFixed(sinkTick < 0.1 ? 2 : sinkTick < 1 ? 1 : 0), px, bottom + 4);
  }

  ctx.restore();
}

function drawAxes(ranges) {
  const { left, right, top, bottom } = chartArea();
  ctx.save();
  ctx.strokeStyle = C.axis;
  ctx.lineWidth = 1.5;
  ctx.setLineDash([]);

  // X axis (bottom)
  ctx.beginPath();
  ctx.moveTo(left, bottom);
  ctx.lineTo(right, bottom);
  ctx.stroke();

  // Y axis (left)
  ctx.beginPath();
  ctx.moveTo(left, top);
  ctx.lineTo(left, bottom);
  ctx.stroke();

  ctx.fillStyle = C.axis;

  // X axis label
  ctx.font = '12px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText(`Sink Rate (${sinkLabel()})`, left + (right - left) / 2, state.canvasH - 4);

  // Y axis label (rotated)
  ctx.save();
  ctx.translate(12, top + (bottom - top) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText(`Speed (${speedLabel()})`, 0, 0);
  ctx.restore();

  ctx.restore();
}

function drawPolarCurve(ranges) {
  const { coeffs } = state;
  const { left, right, top, bottom } = chartArea();

  ctx.save();
  // Clip to chart area to prevent overdraw outside axes
  ctx.beginPath();
  ctx.rect(left, top, right - left, bottom - top);
  ctx.clip();

  ctx.beginPath();
  ctx.strokeStyle = C.polar;
  ctx.lineWidth = 2.5;
  ctx.setLineDash([]);
  ctx.lineJoin = 'round';

  const steps = 300;
  const dv = (ranges.v_max_kmh - ranges.v_min_kmh) / steps;
  let first = true;

  for (let i = 0; i <= steps; i++) {
    const v_kmh = ranges.v_min_kmh + i * dv;
    const w_ms  = polarSink(coeffs, v_kmh) + state.airmass_ms;

    const sinkDisp  = convertSink(w_ms, state.sinkUnit);
    const speedDisp = convertSpeed(v_kmh, state.speedUnit);

    const px = toCanvasX(sinkDisp, ranges);
    const py = toCanvasY(speedDisp, ranges);

    if (first) { ctx.moveTo(px, py); first = false; }
    else ctx.lineTo(px, py);
  }
  ctx.stroke();

  ctx.restore();
}

function drawMcLine(ranges) {
  const mc = computeMcOptimal(state.coeffs, state.mc_ms, state.airmass_ms);
  if (!mc) return;

  const { v_opt, w_opt } = mc;
  const { left, right, top, bottom } = chartArea();

  // Optimal point in display units
  const x_opt = convertSink(w_opt, state.sinkUnit);  // positive value
  const y_opt = convertSpeed(v_opt, state.speedUnit);

  // MC line anchor: in (v, w) space the line starts at (0, mc_ms) — a positive w value
  // representing the thermal climb rate. In display space (X = |w|, Y = speed), this maps
  // to (X = -mc_ms * factor, Y = 0), i.e. off the left edge of the chart.
  // Airmass shifts the polar curve but does NOT move the anchor point.
  const anchor_sink_disp = -state.mc_ms * SINK_UNITS[state.sinkUnit].factor;

  // Slope in display space
  if (Math.abs(x_opt - anchor_sink_disp) < 1e-9) return;
  const slope = (y_opt - 0) / (x_opt - anchor_sink_disp);

  // Compute y at left and right chart edges
  const sink_at_left  = ranges.sink_min_disp;
  const sink_at_right = ranges.sink_max_disp;
  const spd_at_left   = slope * (sink_at_left  - anchor_sink_disp);
  const spd_at_right  = slope * (sink_at_right - anchor_sink_disp);

  ctx.save();
  ctx.beginPath();
  ctx.rect(left, top, right - left, bottom - top);
  ctx.clip();

  ctx.strokeStyle = C.mc;
  ctx.lineWidth = 1.5;
  ctx.setLineDash([6, 5]);
  ctx.beginPath();
  ctx.moveTo(toCanvasX(sink_at_left,  ranges), toCanvasY(spd_at_left,  ranges));
  ctx.lineTo(toCanvasX(sink_at_right, ranges), toCanvasY(spd_at_right, ranges));
  ctx.stroke();
  ctx.setLineDash([]);

  // Tangent point marker
  const px = toCanvasX(x_opt, ranges);
  const py = toCanvasY(y_opt, ranges);
  ctx.beginPath();
  ctx.arc(px, py, 5, 0, Math.PI * 2);
  ctx.fillStyle = C.mc;
  ctx.fill();

  ctx.restore();

  // Annotation label beside the tangent point
  const spd = convertSpeed(v_opt, state.speedUnit).toFixed(1);
  const snk = convertSink(w_opt, state.sinkUnit).toFixed(2);
  ctx.save();
  ctx.font = 'bold 11px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.fillStyle = C.mc;
  ctx.textBaseline = 'middle';
  // Position label to the right of the point, or left if near right edge
  const labelX = px + 8;
  const labelY = py - 12;
  ctx.textAlign = 'left';
  ctx.fillText(`STF: ${spd} ${speedLabel()}`, labelX, labelY);
  ctx.fillText(`Sink: ${snk} ${sinkLabel()}`, labelX, labelY + 14);
  ctx.restore();
}

function drawMinSinkMarker(ranges) {
  const { coeffs } = state;
  const v_ms_kmh = minSinkSpeed(coeffs);
  if (v_ms_kmh < ranges.v_min_kmh || v_ms_kmh > ranges.v_max_kmh) return;

  const w_ms = polarSink(coeffs, v_ms_kmh) + state.airmass_ms;
  const px = toCanvasX(convertSink(w_ms, state.sinkUnit), ranges);
  const py = toCanvasY(convertSpeed(v_ms_kmh, state.speedUnit), ranges);

  const { left, right, top, bottom } = chartArea();
  if (px < left || px > right || py < top || py > bottom) return;

  ctx.save();
  ctx.beginPath();
  ctx.arc(px, py, 4, 0, Math.PI * 2);
  ctx.strokeStyle = C.minsink;
  ctx.lineWidth = 2;
  ctx.stroke();
  ctx.restore();
}

function drawHoverMarker(ranges) {
  const { v_kmh, w_ms } = state.hoverPoint;
  const px = toCanvasX(convertSink(w_ms, state.sinkUnit), ranges);
  const py = toCanvasY(convertSpeed(v_kmh, state.speedUnit), ranges);

  const { left, right, top, bottom } = chartArea();
  if (px < left || px > right || py < top || py > bottom) return;

  ctx.save();
  // Crosshair lines
  ctx.strokeStyle = C.hover;
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 3]);
  ctx.globalAlpha = 0.6;
  ctx.beginPath();
  ctx.moveTo(left, py);
  ctx.lineTo(px, py);
  ctx.moveTo(px, py);
  ctx.lineTo(px, bottom);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.globalAlpha = 1;

  // Point circle
  ctx.beginPath();
  ctx.arc(px, py, 5, 0, Math.PI * 2);
  ctx.fillStyle = C.hover;
  ctx.fill();
  ctx.restore();
}

// === MAIN REDRAW ===

function redraw() {
  if (!state.coeffs || state.canvasW === 0) return;

  ctx.clearRect(0, 0, state.canvasW, state.canvasH);

  const entry = state.polars[state.selectedIndex];
  const ranges = computeRanges(entry, state.coeffs);
  state.ranges = ranges;

  drawGrid(ranges);
  drawAxes(ranges);
  drawPolarCurve(ranges);
  drawMcLine(ranges);
  drawMinSinkMarker(ranges);
  if (state.hoverPoint) drawHoverMarker(ranges);
}

// === TOOLTIP ===

const tooltip = document.getElementById('tooltip');
const tipSpeed = document.getElementById('tip-speed');
const tipSink  = document.getElementById('tip-sink');

function showTooltip(px, py, v_kmh, w_ms) {
  tipSpeed.textContent = `${convertSpeed(v_kmh, state.speedUnit).toFixed(1)} ${speedLabel()}`;
  tipSink.textContent  = `Sink: ${convertSink(w_ms, state.sinkUnit).toFixed(2)} ${sinkLabel()}`;

  // Offset above and to the right of pointer; flip if near right/top edge
  const container = document.getElementById('chart-container');
  const cw = container.clientWidth;
  let lx = px + 14;
  let ly = py - 64;
  if (lx + 120 > cw) lx = px - 130;
  if (ly < 4) ly = py + 10;

  tooltip.style.left = `${lx}px`;
  tooltip.style.top  = `${ly}px`;
  tooltip.removeAttribute('hidden');
}

function hideTooltip() {
  tooltip.setAttribute('hidden', '');
}

// === INTERACTION ===

function handlePointer(e) {
  if (!state.coeffs || !state.ranges) return;

  const rect = canvas.getBoundingClientRect();
  const px = e.clientX - rect.left;
  const py = e.clientY - rect.top;
  const { left, right, top, bottom } = chartArea();

  if (px < left || px > right || py < top || py > bottom) {
    if (state.hoverPoint) { state.hoverPoint = null; redraw(); hideTooltip(); }
    return;
  }

  // Convert pointer Y to speed, clamp to range, then compute w at that speed
  const speedDisp  = fromCanvasY(py, state.ranges);
  const v_kmh = Math.max(
    state.ranges.v_min_kmh,
    Math.min(state.ranges.v_max_kmh,
             speedDisp / SPEED_UNITS[state.speedUnit].factor)
  );
  const w_ms = polarSink(state.coeffs, v_kmh) + state.airmass_ms;

  state.hoverPoint = { v_kmh, w_ms };
  redraw();
  showTooltip(px, py, v_kmh, w_ms);
}

canvas.addEventListener('pointermove', handlePointer);
canvas.addEventListener('pointerdown', handlePointer);
canvas.addEventListener('pointerleave', () => {
  state.hoverPoint = null;
  redraw();
  hideTooltip();
});

// === UI CONTROLS ===

function updateSliderLabels() {
  const airmassKts = parseFloat(document.getElementById('airmass-slider').value);
  const mcKts = parseFloat(document.getElementById('mc-slider').value);

  // Show slider values in the selected sink unit
  const airmassMs = ktsToMs(airmassKts);
  const mcMs = ktsToMs(mcKts);
  const factor = SINK_UNITS[state.sinkUnit].factor;

  document.getElementById('airmass-value').textContent =
    (airmassMs >= 0 ? '+' : '') + (airmassMs * factor).toFixed(2);
  document.getElementById('mc-value').textContent =
    (mcMs * factor).toFixed(2);

  document.querySelectorAll('.airmass-unit, .mc-unit').forEach(el => {
    el.textContent = sinkLabel();
  });
}

function initControls() {
  const select = document.getElementById('aircraft-select');

  select.innerHTML = state.polars
    .map((p, i) => `<option value="${i}">${p.name}</option>`)
    .join('');
  select.removeAttribute('disabled');

  // Restore default selection
  select.value = String(state.selectedIndex);

  select.addEventListener('change', () => {
    state.selectedIndex = parseInt(select.value, 10);
    const entry = state.polars[state.selectedIndex];
    state.coeffs = fitPolar(entry);
    state.hoverPoint = null;
    hideTooltip();
    redraw();
  });

  document.querySelectorAll('input[name="speed-unit"]').forEach(r => {
    r.addEventListener('change', () => {
      if (!r.checked) return;
      state.speedUnit = r.value;
      redraw();
    });
  });

  document.querySelectorAll('input[name="sink-unit"]').forEach(r => {
    r.addEventListener('change', () => {
      if (!r.checked) return;
      state.sinkUnit = r.value;
      updateSliderLabels();
      redraw();
    });
  });

  const airmassSlider = document.getElementById('airmass-slider');
  airmassSlider.addEventListener('input', () => {
    state.airmass_ms = ktsToMs(parseFloat(airmassSlider.value));
    updateSliderLabels();
    redraw();
  });

  const mcSlider = document.getElementById('mc-slider');
  mcSlider.addEventListener('input', () => {
    state.mc_ms = ktsToMs(parseFloat(mcSlider.value));
    updateSliderLabels();
    redraw();
  });

  updateSliderLabels();
}

// === INITIALISATION ===

function debounce(fn, ms) {
  let timer;
  return (...args) => { clearTimeout(timer); timer = setTimeout(() => fn(...args), ms); };
}

async function init() {
  resolveColours();
  resizeCanvas();

  window.addEventListener('resize', debounce(() => {
    resizeCanvas();
    redraw();
  }, 150));

  // Also re-resolve colours if system theme changes
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
    resolveColours();
    redraw();
  });

  try {
    state.polars = await fetchPolars();
    if (state.polars.length === 0) throw new Error('No polars parsed');

    // Default to Discus or first entry
    const defaultIdx = state.polars.findIndex(p => p.name.startsWith('Discus'));
    state.selectedIndex = defaultIdx >= 0 ? defaultIdx : 0;
    state.coeffs = fitPolar(state.polars[state.selectedIndex]);

    initControls();
    redraw();
  } catch (err) {
    console.error('Failed to load polars:', err);
    const select = document.getElementById('aircraft-select');
    select.innerHTML = '<option>Error loading polars</option>';
    select.disabled = true;
  }
}

init();
