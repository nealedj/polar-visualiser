// === CONSTANTS & CONFIG ===

const POLAR_URL =
  'https://raw.githubusercontent.com/XCSoar/XCSoar/master/src/Polar/PolarStore.cpp';

const MARGIN = { top: 24, right: 30, bottom: 52, left: 72 };

// Conversion factors from internal units (km/h, m/s) to display units
const SPEED_UNITS = {
  kts: { label: 'kts',  factor: 1 / 1.852 },
  kmh: { label: 'km/h', factor: 1 },
  mph: { label: 'mph',  factor: 1 / 1.60934 },
};

const SINK_UNITS = {
  kts: { label: 'kts',   factor: 1 / 0.51444 },
  fts: { label: 'ft/min', factor: 60 / 0.3048 },
  ms:  { label: 'm/s',   factor: 1 },
};

// Slider config per sink unit: step, airmass range, MacCready max
const SINK_SLIDER = {
  kts: { step: 0.5,  airmassMin: -10,   airmassMax: 10,   mcMax: 10,   decimals: 1 },
  fts: { step: 100,  airmassMin: -2000,  airmassMax: 2000, mcMax: 2000, decimals: 0 },
  ms:  { step: 0.2,  airmassMin: -5.5,   airmassMax: 5.5,  mcMax: 5.5,  decimals: 1 },
};

// Resolved once from CSS variables — canvas can't use CSS vars directly
const C = {
  polar:   '#1976d2',
  mc:      '#e53935',
  hover:   '#2e7d32',
  grid:    '#e0e0e0',
  axis:    '#424242',
  text:    '#212121',
  zeroline:'#888888',
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

// Signed vertical rate: negative = sink, positive = lift/climb
function convertRate(w_ms, unit) {
  return w_ms * SINK_UNITS[unit].factor;
}

function speedLabel() { return SPEED_UNITS[state.speedUnit].label; }
function sinkLabel()  { return SINK_UNITS[state.sinkUnit].label; }

function ktsToMs(kts) { return kts * 0.51444; }

function msToSinkDisp(ms) { return ms * SINK_UNITS[state.sinkUnit].factor; }
function sinkDispToMs(val) { return val / SINK_UNITS[state.sinkUnit].factor; }

// === DATA FETCHING & PARSING ===

async function fetchPolars() {
  const resp = await fetch(POLAR_URL);
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  const text = await resp.text();
  return parsePolars(text);
}

function parsePolars(cpp) {
  const results = [];

  // Match each { "Name" [optional /* comment */], num x12 }
  const BLOCK_RE =
    /\{\s*"([^"]+)"\s*(?:\/\*[^*]*\*+(?:[^/*][^*]*\*+)*\/)?\s*,([^}]+)\}/g;

  let m;
  while ((m = BLOCK_RE.exec(cpp)) !== null) {
    const name = m[1].trim();
    const nums = m[2].match(/-?[\d]+\.?[\d]*/g);
    if (!nums || nums.length < 12) continue;

    const [empty_mass, max_ballast, v1, w1, v2, w2, v3, w3, wing_area, ref_mass, v_no, max_speed]
      = nums.map(Number);

    // Sanity: speeds positive, sink negative
    if (v1 <= 0 || v2 <= 0 || v3 <= 0) continue;
    if (w1 >= 0 || w2 >= 0 || w3 >= 0) continue;

    results.push({ name, empty_mass, max_ballast, v1, w1, v2, w2, v3, w3,
                   wing_area, ref_mass, v_no, max_speed });
  }

  results.sort((a, b) => a.name.localeCompare(b.name));

  const seen = new Set();
  return results.filter(p => {
    if (seen.has(p.name)) return false;
    seen.add(p.name);
    return true;
  });
}

// === POLAR MATH ===

// Fit quadratic w = a*v^2 + b*v + c through 3 points using Cramer's rule.
// For real polars a < 0 (parabola opens downward — sink worsens at both speed extremes).
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

  if (a >= 0) return null;
  return { a, b, c };
}

function polarSink(coeffs, v_kmh) {
  return coeffs.a * v_kmh * v_kmh + coeffs.b * v_kmh + coeffs.c;
}

// MacCready optimal speed: tangent from (v=0, w=mc_ms) to the shifted polar.
// v_opt = sqrt((c + airmass_ms - mc_ms) / a)
function computeMcOptimal(coeffs, mc_ms, airmass_ms) {
  const discriminant = (coeffs.c + airmass_ms - mc_ms) / coeffs.a;
  if (discriminant <= 0) return null;

  const v_opt = Math.sqrt(discriminant);
  const w_opt = polarSink(coeffs, v_opt) + airmass_ms;

  return { v_opt, w_opt };
}

// Min-sink speed: vertex of parabola at dw/dv = 0 → v = -b/(2a)
function minSinkSpeed(coeffs) {
  return -coeffs.b / (2 * coeffs.a);
}

// === CHART GEOMETRY ===
// X axis = speed (left → right, increasing)
// Y axis = vertical rate (bottom → top, negative sink below 0, positive lift above 0)

function chartArea() {
  return {
    left:   MARGIN.left,
    right:  state.canvasW - MARGIN.right,
    top:    MARGIN.top,
    bottom: state.canvasH - MARGIN.bottom,
  };
}

function computeRanges(entry, coeffs) {
  const v_minsink = minSinkSpeed(coeffs);
  const v_min_kmh = Math.max(40, Math.min(entry.v1 * 0.75, v_minsink * 0.75));
  const v_max_kmh = entry.max_speed;

  // Find worst (most negative) sink over the displayed speed range
  let w_min_ms = 0;
  for (let v = v_min_kmh; v <= v_max_kmh; v += 2) {
    const w = polarSink(coeffs, v) + state.airmass_ms;
    if (w < w_min_ms) w_min_ms = w;
  }

  // Y axis bottom: most negative value + 15% padding
  const w_min_disp = convertRate(w_min_ms, state.sinkUnit) * 1.15;
  // Y axis top: positive region — enough to show 0 clearly and the MC anchor
  const totalNeg = Math.abs(w_min_disp);
  const mc_disp  = state.mc_ms * SINK_UNITS[state.sinkUnit].factor;
  const w_max_disp = Math.max(totalNeg * 0.3, mc_disp + totalNeg * 0.06);

  return {
    v_min_kmh,
    v_max_kmh,
    v_min_disp: convertSpeed(v_min_kmh, state.speedUnit),
    v_max_disp: convertSpeed(v_max_kmh, state.speedUnit),
    w_min_disp,   // Y bottom (negative)
    w_max_disp,   // Y top (positive)
  };
}

// Speed → canvas X
function toCanvasX(speedDisp, ranges) {
  const { left, right } = chartArea();
  return left + (speedDisp - ranges.v_min_disp) /
    (ranges.v_max_disp - ranges.v_min_disp) * (right - left);
}

// Signed vertical rate → canvas Y (positive rate = up = lower pixel Y)
function toCanvasY(rateDisp, ranges) {
  const { top, bottom } = chartArea();
  return bottom - (rateDisp - ranges.w_min_disp) /
    (ranges.w_max_disp - ranges.w_min_disp) * (bottom - top);
}

// Canvas X → speed in display units
function fromCanvasX(px, ranges) {
  const { left, right } = chartArea();
  return ranges.v_min_disp + (px - left) /
    (right - left) * (ranges.v_max_disp - ranges.v_min_disp);
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
  const mag = Math.pow(10, Math.floor(Math.log10(Math.abs(rough))));
  const normalised = rough / mag;
  let nice;
  if (normalised < 1.5)      nice = 1;
  else if (normalised < 3.5) nice = 2;
  else if (normalised < 7.5) nice = 5;
  else                       nice = 10;
  return nice * mag;
}

function drawGrid(ranges) {
  const { left, right, top, bottom } = chartArea();
  ctx.save();
  ctx.font = '11px -apple-system, BlinkMacSystemFont, sans-serif';

  // --- Vertical grid lines: speed ticks ---
  const speedTick = niceTickInterval(ranges.v_max_disp - ranges.v_min_disp, 6);
  const speedStart = Math.ceil(ranges.v_min_disp / speedTick) * speedTick;

  ctx.strokeStyle = C.grid;
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 4]);
  ctx.fillStyle = C.axis;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';

  for (let v = speedStart; v <= ranges.v_max_disp + speedTick * 0.01; v += speedTick) {
    const px = toCanvasX(v, ranges);
    if (px < left - 1 || px > right + 1) continue;
    ctx.beginPath();
    ctx.moveTo(px, top);
    ctx.lineTo(px, bottom);
    ctx.stroke();
    ctx.fillText(v.toFixed(speedTick < 1 ? 1 : 0), px, bottom + 4);
  }

  // --- Horizontal grid lines: rate ticks ---
  const rateRange = ranges.w_max_disp - ranges.w_min_disp;
  const rateTick  = niceTickInterval(rateRange, 6);
  const rateStart = Math.ceil(ranges.w_min_disp / rateTick) * rateTick;

  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';

  for (let w = rateStart; w <= ranges.w_max_disp + rateTick * 0.01; w += rateTick) {
    if (Math.abs(w) < rateTick * 0.01) continue; // skip 0 — drawn separately
    const py = toCanvasY(w, ranges);
    if (py < top - 1 || py > bottom + 1) continue;
    ctx.strokeStyle = C.grid;
    ctx.setLineDash([3, 4]);
    ctx.beginPath();
    ctx.moveTo(left, py);
    ctx.lineTo(right, py);
    ctx.stroke();
    ctx.fillText(w.toFixed(rateTick < 0.1 ? 2 : rateTick < 1 ? 1 : 0), left - 6, py);
  }

  ctx.restore();
}

function drawAxes(ranges) {
  const { left, right, top, bottom } = chartArea();
  ctx.save();
  ctx.setLineDash([]);
  ctx.fillStyle = C.axis;
  ctx.font = '12px -apple-system, BlinkMacSystemFont, sans-serif';

  // Border axes
  ctx.strokeStyle = C.axis;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(left, top);
  ctx.lineTo(left, bottom);
  ctx.lineTo(right, bottom);
  ctx.stroke();

  // Zero-rate line — prominent horizontal line at w=0
  const py0 = toCanvasY(0, ranges);
  if (py0 >= top && py0 <= bottom) {
    ctx.strokeStyle = C.zeroline;
    ctx.lineWidth = 1.5;
    ctx.setLineDash([]);
    ctx.beginPath();
    ctx.moveTo(left, py0);
    ctx.lineTo(right, py0);
    ctx.stroke();
    // "0" label
    ctx.fillStyle = C.zeroline;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText('0', left - 6, py0);
  }

  // X axis label: Speed
  ctx.fillStyle = C.axis;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'bottom';
  ctx.fillText(`Speed (${speedLabel()})`, left + (right - left) / 2, state.canvasH - 4);

  // Y axis label: Climb / Sink Rate (rotated)
  ctx.save();
  ctx.translate(13, top + (bottom - top) / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText(`Climb / Sink Rate (${sinkLabel()})`, 0, 0);
  ctx.restore();

  ctx.restore();
}

function drawPolarCurve(ranges) {
  const { coeffs } = state;
  const { left, right, top, bottom } = chartArea();

  ctx.save();
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
  let started = false;

  for (let i = 0; i <= steps; i++) {
    const v_kmh = ranges.v_min_kmh + i * dv;
    const w_ms  = polarSink(coeffs, v_kmh) + state.airmass_ms;

    const px = toCanvasX(convertSpeed(v_kmh, state.speedUnit), ranges);
    const py = toCanvasY(convertRate(w_ms,   state.sinkUnit),  ranges);

    if (!started) { ctx.moveTo(px, py); started = true; }
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

  // In (speed, rate) display space the MC line passes through:
  //   anchor: (0, mc_ms_disp) — at zero speed, the thermal climb rate (off left of chart)
  //   tangent point: (v_opt_disp, w_opt_disp)
  // Line equation: rate = mc_disp + slope * speed
  const mc_disp    = state.mc_ms * SINK_UNITS[state.sinkUnit].factor;
  const v_opt_disp = convertSpeed(v_opt, state.speedUnit);
  const w_opt_disp = convertRate(w_opt, state.sinkUnit);

  if (Math.abs(v_opt_disp) < 1e-9) return;
  const slope = (w_opt_disp - mc_disp) / v_opt_disp;

  const rate_at_left  = mc_disp + slope * ranges.v_min_disp;
  const rate_at_right = mc_disp + slope * ranges.v_max_disp;

  // Anchor dot on Y axis at mc_disp height (at v=0, outside chart area)
  const anchorY = toCanvasY(mc_disp, ranges);
  if (anchorY >= top && anchorY <= bottom) {
    ctx.save();
    ctx.beginPath();
    ctx.arc(left, anchorY, 4, 0, Math.PI * 2);
    ctx.fillStyle = C.mc;
    ctx.fill();
    ctx.restore();
  }

  ctx.save();
  ctx.beginPath();
  ctx.rect(left, top, right - left, bottom - top);
  ctx.clip();

  ctx.strokeStyle = C.mc;
  ctx.lineWidth = 1.5;
  ctx.setLineDash([6, 5]);
  ctx.beginPath();
  ctx.moveTo(toCanvasX(ranges.v_min_disp, ranges), toCanvasY(rate_at_left,  ranges));
  ctx.lineTo(toCanvasX(ranges.v_max_disp, ranges), toCanvasY(rate_at_right, ranges));
  ctx.stroke();
  ctx.setLineDash([]);

  // Tangent point dot
  const px = toCanvasX(v_opt_disp, ranges);
  const py = toCanvasY(w_opt_disp, ranges);
  ctx.beginPath();
  ctx.arc(px, py, 5, 0, Math.PI * 2);
  ctx.fillStyle = C.mc;
  ctx.fill();
  ctx.restore();

  // Annotation
  const spd = v_opt_disp.toFixed(1);
  const snk = Math.abs(w_opt_disp).toFixed(2);
  ctx.save();
  ctx.font = 'bold 11px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.fillStyle = C.mc;
  const labelX = Math.min(px + 8, state.canvasW - MARGIN.right - 80);
  const labelY = py - 18;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`STF: ${spd} ${speedLabel()}`, labelX, labelY);
  ctx.fillText(`Sink: ${snk} ${sinkLabel()}`, labelX, labelY + 14);
  ctx.restore();
}

function drawMinSinkMarker(ranges) {
  const { coeffs } = state;
  const v_ms_kmh = minSinkSpeed(coeffs);
  if (v_ms_kmh < ranges.v_min_kmh || v_ms_kmh > ranges.v_max_kmh) return;

  const w_ms = polarSink(coeffs, v_ms_kmh) + state.airmass_ms;
  const px = toCanvasX(convertSpeed(v_ms_kmh, state.speedUnit), ranges);
  const py = toCanvasY(convertRate(w_ms, state.sinkUnit), ranges);

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
  const px = toCanvasX(convertSpeed(v_kmh, state.speedUnit), ranges);
  const py = toCanvasY(convertRate(w_ms, state.sinkUnit), ranges);

  const { left, right, top, bottom } = chartArea();
  if (px < left || px > right || py < top || py > bottom) return;

  ctx.save();
  ctx.strokeStyle = C.hover;
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 3]);
  ctx.globalAlpha = 0.6;
  ctx.beginPath();
  // Vertical line at hover speed
  ctx.moveTo(px, top);
  ctx.lineTo(px, bottom);
  // Horizontal line at hover rate
  ctx.moveTo(left, py);
  ctx.lineTo(right, py);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.globalAlpha = 1;

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

  const entry  = state.polars[state.selectedIndex];
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

const tooltip  = document.getElementById('tooltip');
const tipSpeed = document.getElementById('tip-speed');
const tipSink  = document.getElementById('tip-sink');

function showTooltip(px, py, v_kmh, w_ms) {
  const rateDisp = convertRate(w_ms, state.sinkUnit);
  const sign = rateDisp >= 0 ? '+' : '';
  tipSpeed.textContent = `${convertSpeed(v_kmh, state.speedUnit).toFixed(1)} ${speedLabel()}`;
  tipSink.textContent  = `${sign}${rateDisp.toFixed(2)} ${sinkLabel()}`;

  const container = document.getElementById('chart-container');
  const cw = container.clientWidth;
  let lx = px + 14;
  let ly = py - 64;
  if (lx + 130 > cw) lx = px - 140;
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

  // Convert pointer X → speed → find sink on polar at that speed
  const speedDisp = fromCanvasX(px, state.ranges);
  const v_kmh = Math.max(
    state.ranges.v_min_kmh,
    Math.min(state.ranges.v_max_kmh, speedDisp / SPEED_UNITS[state.speedUnit].factor)
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
  const airmassDisp = parseFloat(document.getElementById('airmass-slider').value);
  const mcDisp      = parseFloat(document.getElementById('mc-slider').value);
  const dec = SINK_SLIDER[state.sinkUnit].decimals;

  document.getElementById('airmass-value').textContent =
    (airmassDisp >= 0 ? '+' : '') + airmassDisp.toFixed(dec);
  document.getElementById('mc-value').textContent = mcDisp.toFixed(dec);

  document.querySelectorAll('.airmass-unit, .mc-unit').forEach(el => {
    el.textContent = sinkLabel();
  });
}

function updateSliderConfigs() {
  const cfg = SINK_SLIDER[state.sinkUnit];
  const airmassSlider = document.getElementById('airmass-slider');
  const mcSlider = document.getElementById('mc-slider');

  airmassSlider.min   = cfg.airmassMin;
  airmassSlider.max   = cfg.airmassMax;
  airmassSlider.step  = cfg.step;
  airmassSlider.value = msToSinkDisp(state.airmass_ms);

  mcSlider.min   = 0;
  mcSlider.max   = cfg.mcMax;
  mcSlider.step  = cfg.step;
  mcSlider.value = msToSinkDisp(state.mc_ms);
}

function initControls() {
  const select = document.getElementById('aircraft-select');

  select.innerHTML = state.polars
    .map((p, i) => `<option value="${i}">${p.name}</option>`)
    .join('');
  select.removeAttribute('disabled');
  select.value = String(state.selectedIndex);

  select.addEventListener('change', () => {
    state.selectedIndex = parseInt(select.value, 10);
    state.coeffs = fitPolar(state.polars[state.selectedIndex]);
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
      updateSliderConfigs();
      updateSliderLabels();
      redraw();
    });
  });

  document.getElementById('airmass-slider').addEventListener('input', function () {
    state.airmass_ms = sinkDispToMs(parseFloat(this.value));
    updateSliderLabels();
    redraw();
  });

  document.getElementById('mc-slider').addEventListener('input', function () {
    state.mc_ms = sinkDispToMs(parseFloat(this.value));
    updateSliderLabels();
    redraw();
  });

  updateSliderConfigs();
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

  window.addEventListener('resize', debounce(() => { resizeCanvas(); redraw(); }, 150));
  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
    resolveColours(); redraw();
  });

  try {
    state.polars = await fetchPolars();
    if (state.polars.length === 0) throw new Error('No polars parsed');

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
