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
  stall:   '#e65100',
  compare: '#00897b',
};

// === APPLICATION STATE ===

const state = {
  polars: [],
  selectedIndex: 0,
  coeffs: null,        // { a, b, c } — quadratic fit at ref_mass (unballasted)
  activeCoeffs: null,  // coeffs scaled for current ballast
  ballast_kg: 0,
  airmass_ms: 0,       // m/s, positive = lift
  mc_ms: 0,            // m/s, MacCready setting (positive)
  speedUnit: 'kts',
  sinkUnit: 'kts',
  hoverPoint: null,    // { v_kmh, w_ms } | null
  ranges: null,        // cached axis ranges
  canvasW: 0,
  canvasH: 0,
  // comparison mode
  compareMode: false,
  compareIndex: 0,
  compareCoeffs: null,
  compareActiveCoeffs: null,
  compareBallast_kg: 0,
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

    // Struct order: reference_mass, max_ballast, v1,w1,v2,w2,v3,w3,
    //               wing_area, v_no_ms (m/s, 0=unknown), contest_handicap, empty_mass
    const [reference_mass, max_ballast, v1, w1, v2, w2, v3, w3,
           wing_area, v_no_ms, contest_handicap, empty_mass]
      = nums.map(Number);

    // Sanity: speeds positive, sink negative
    if (v1 <= 0 || v2 <= 0 || v3 <= 0) continue;
    if (w1 >= 0 || w2 >= 0 || w3 >= 0) continue;

    // v_no is stored in m/s; convert to km/h for the chart upper bound.
    // If unknown (0), extend 50% beyond the last measured point.
    const v_max_kmh = v_no_ms > 0 ? v_no_ms * 3.6 : v3 * 1.5;

    results.push({ name, reference_mass, max_ballast, v1, w1, v2, w2, v3, w3,
                   wing_area, v_max_kmh });
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

// Scale polar coefficients for added ballast (wing-loading scaling).
// Both axes scale by k = sqrt((ref_mass + ballast) / ref_mass):
//   a_new = a/k,  b_new = b,  c_new = c*k
function applyBallast(coeffs, ref_mass, ballast_kg) {
  if (ballast_kg <= 0 || !ref_mass) return coeffs;
  const k = Math.sqrt((ref_mass + ballast_kg) / ref_mass);
  return { a: coeffs.a / k, b: coeffs.b, c: coeffs.c * k };
}

function updateActiveCoeffs() {
  const entry = state.polars[state.selectedIndex];
  state.activeCoeffs = applyBallast(state.coeffs, entry.reference_mass, state.ballast_kg);
}

function updateCompareCoeffs() {
  if (!state.compareMode || state.compareCoeffs === null) return;
  const entry = state.polars[state.compareIndex];
  state.compareActiveCoeffs = applyBallast(
    state.compareCoeffs, entry.reference_mass, state.compareBallast_kg
  );
}

// Estimate stall speed in km/h using wing loading formula, or fallback to v1 × 0.55.
// CLmax ≈ 2.0 matches published stall speeds for typical competition gliders.
function computeStallSpeed(entry, ballast_kg) {
  const g = 9.81, rho = 1.225, CLmax = 2.0;
  const mass = entry.reference_mass + ballast_kg;
  if (entry.wing_area > 0) {
    return Math.sqrt(2 * mass * g / (rho * entry.wing_area * CLmax)) * 3.6;
  }
  // No wing area recorded: ~48% of first measured speed, scaled for mass
  const k = Math.sqrt(mass / entry.reference_mass);
  return entry.v1 * 0.48 * k;
}

function polarSink(coeffs, v_kmh) {
  return coeffs.a * v_kmh * v_kmh + coeffs.b * v_kmh + coeffs.c;
}

// MacCready optimal speed: tangent from (v=0, w=mc_ms) to the shifted polar.
// v_opt = sqrt((c + airmass_ms - mc_ms) / a)
function computeMcOptimal(coeffs, mc_ms, airmass_ms) {
  const discriminant = (coeffs.c + airmass_ms - mc_ms) / coeffs.a;
  if (discriminant <= 0) return null;

  // Never draw STF below min-sink speed — clamp to the vertex of the parabola
  const v_opt = Math.max(Math.sqrt(discriminant), minSinkSpeed(coeffs));
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
  // X axis starts at v=0 so the Y axis IS the zero-speed axis.
  // This makes the best-glide line (MC=0) originate from (0,0) at the left axis.
  const v_min_kmh = 0;
  let v_max_kmh = entry.v_max_kmh;
  const stall_kmh = computeStallSpeed(entry, state.ballast_kg);

  // Find worst (most negative) and best (most positive) rate over the visible range,
  // starting from stall speed so sub-stall extrapolations don't distort the axis.
  const scan_start = Math.max(stall_kmh, 5);
  let w_min_ms = polarSink(coeffs, v_max_kmh) + state.airmass_ms;
  let w_max_ms = polarSink(coeffs, scan_start) + state.airmass_ms;
  for (let v = scan_start; v <= v_max_kmh; v += 2) {
    const w = polarSink(coeffs, v) + state.airmass_ms;
    if (w < w_min_ms) w_min_ms = w;
    if (w > w_max_ms) w_max_ms = w;
  }

  // Envelope comparison glider's curve into the axis ranges
  if (state.compareMode && state.compareActiveCoeffs) {
    const cEntry = state.polars[state.compareIndex];
    const cVmax  = cEntry.v_max_kmh;
    const cStall = computeStallSpeed(cEntry, state.compareBallast_kg);
    const cStart = Math.max(cStall, 5);
    if (cVmax > v_max_kmh) v_max_kmh = cVmax;
    for (let v = cStart; v <= cVmax; v += 2) {
      const w = polarSink(state.compareActiveCoeffs, v) + state.airmass_ms;
      if (w < w_min_ms) w_min_ms = w;
      if (w > w_max_ms) w_max_ms = w;
    }
  }

  // Y axis bottom: worst visible sink + 15% padding
  const w_min_disp = convertRate(w_min_ms, state.sinkUnit) * 1.15;
  // Y axis top: highest of — curve peak, MC anchor, or 15% of range for breathing room
  const totalNeg = Math.abs(w_min_disp);
  const mc_disp        = state.mc_ms * SINK_UNITS[state.sinkUnit].factor;
  const curve_top_disp = convertRate(w_max_ms, state.sinkUnit);
  const w_max_disp = Math.max(totalNeg * 0.15, mc_disp + totalNeg * 0.06, curve_top_disp * 1.15);

  return {
    v_min_kmh,
    v_max_kmh,
    stall_kmh,
    v_min_disp: 0,   // 0 speed is 0 in every unit
    v_max_disp: convertSpeed(v_max_kmh, state.speedUnit),
    w_min_disp,
    w_max_disp,
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
  apply('polar',   '--color-polar');
  apply('mc',      '--color-mc');
  apply('hover',   '--color-hover');
  apply('grid',    '--color-grid');
  apply('axis',    '--color-axis');
  apply('text',    '--color-text');
  apply('compare', '--color-compare');
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

function drawPolarCurve(ranges, coeffs, color, stallKmh, curveVmax, label) {
  const { left, right, top, bottom } = chartArea();

  ctx.save();
  ctx.beginPath();
  ctx.rect(left, top, right - left, bottom - top);
  ctx.clip();

  ctx.beginPath();
  ctx.strokeStyle = color;
  ctx.lineWidth = 2.5;
  ctx.setLineDash([]);
  ctx.lineJoin = 'round';

  // Start curve at estimated stall speed — leave a gap to the left axis
  const v_start = Math.max(stallKmh, ranges.v_min_kmh);
  const steps = 300;
  const dv = (curveVmax - v_start) / steps;
  let started = false;
  let lastPx, lastPy;

  for (let i = 0; i <= steps; i++) {
    const v_kmh = v_start + i * dv;
    const w_ms  = polarSink(coeffs, v_kmh) + state.airmass_ms;

    const px = toCanvasX(convertSpeed(v_kmh, state.speedUnit), ranges);
    const py = toCanvasY(convertRate(w_ms,   state.sinkUnit),  ranges);

    if (!started) { ctx.moveTo(px, py); started = true; }
    else ctx.lineTo(px, py);
    lastPx = px; lastPy = py;
  }
  ctx.stroke();

  if (label && lastPx !== undefined) {
    ctx.font = 'bold 10px -apple-system, BlinkMacSystemFont, sans-serif';
    ctx.fillStyle = color;
    ctx.textAlign = 'right';
    ctx.textBaseline = 'bottom';
    ctx.fillText(label, lastPx - 4, lastPy - 4);
  }

  ctx.restore();
}

function drawStallMarker(ranges) {
  const stall_disp = convertSpeed(ranges.stall_kmh, state.speedUnit);
  const px = toCanvasX(stall_disp, ranges);
  const { left, right, top, bottom } = chartArea();
  if (px < left || px > right) return;

  ctx.save();

  // Faint vertical dashed line at estimated stall speed
  ctx.strokeStyle = C.stall;
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 4]);
  ctx.globalAlpha = 0.45;
  ctx.beginPath();
  ctx.moveTo(px, top);
  ctx.lineTo(px, bottom);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.globalAlpha = 1;

  // Small label at top of the line
  ctx.font = '10px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.fillStyle = C.stall;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Est. stall', px, top + 3);
  ctx.fillText(`~${stall_disp.toFixed(0)} ${speedLabel()}`, px, top + 15);

  ctx.restore();
}

function drawMcLine(ranges, coeffs, color, labelBelow = false) {
  const mc = computeMcOptimal(coeffs, state.mc_ms, state.airmass_ms);
  if (!mc) return;

  const { v_opt, w_opt } = mc;
  const { left, right, top, bottom } = chartArea();

  // MC line passes through (v=0, w=mc_disp) and the tangent point.
  // Line equation in display space: rate(v) = mc_disp + slope * v
  const mc_disp    = state.mc_ms * SINK_UNITS[state.sinkUnit].factor;
  const v_opt_disp = convertSpeed(v_opt, state.speedUnit);
  const w_opt_disp = convertRate(w_opt, state.sinkUnit);

  if (Math.abs(v_opt_disp) < 1e-9) return;
  const slope = (w_opt_disp - mc_disp) / v_opt_disp;

  // The true line originates at (v=0, w=mc_disp). v=0 is left of the chart,
  // so we compute the line's canvas coords for both the left chart boundary
  // and the right chart boundary, then draw it clipped to the chart area.
  const rate_at_left  = mc_disp + slope * ranges.v_min_disp;
  const rate_at_right = mc_disp + slope * ranges.v_max_disp;

  ctx.save();
  ctx.beginPath();
  ctx.rect(left, top, right - left, bottom - top);
  ctx.clip();

  ctx.strokeStyle = color;
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
  ctx.fillStyle = color;
  ctx.fill();
  ctx.restore();

  // Annotation: "Best Glide" at MC=0, "STF" otherwise
  const label  = state.mc_ms < 1e-9 ? 'Best Glide' : 'STF';
  const spd = v_opt_disp.toFixed(1);
  ctx.save();
  ctx.font = 'bold 11px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.fillStyle = color;
  const labelX = Math.min(px + 8, state.canvasW - MARGIN.right - 90);
  const labelY = labelBelow ? py + 8 : py - 10;
  ctx.textAlign = 'left';
  ctx.textBaseline = 'top';
  ctx.fillText(`${label}: ${spd} ${speedLabel()}`, labelX, labelY);
  ctx.restore();

  // Cross-country speed: where the MC line crosses y=0 (zero-rate axis).
  // From rate(v) = mc_disp + slope*v = 0  →  v_cc = -mc_disp / slope
  // At MC=0 this is 0 (origin); at MC>0 it is a positive speed.
  const v_cc_disp = Math.abs(slope) > 1e-12 ? -mc_disp / slope : 0;
  const x_cc = toCanvasX(v_cc_disp, ranges);
  const y0   = toCanvasY(0, ranges);

  if (x_cc >= left && x_cc <= right && y0 >= top && y0 <= bottom) {
    ctx.save();
    // Diamond marker on the zero line
    ctx.beginPath();
    ctx.moveTo(x_cc,     y0 - 6);
    ctx.lineTo(x_cc + 5, y0);
    ctx.lineTo(x_cc,     y0 + 6);
    ctx.lineTo(x_cc - 5, y0);
    ctx.closePath();
    ctx.fillStyle = color;
    ctx.fill();
    // Speed label just above the zero line
    ctx.font = '10px -apple-system, BlinkMacSystemFont, sans-serif';
    ctx.fillStyle = color;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    const v_cc_kmh = v_cc_disp / SPEED_UNITS[state.speedUnit].factor;
    const line1 = state.speedUnit === 'kmh'
      ? `${v_cc_kmh.toFixed(0)} km/h`
      : `${v_cc_disp.toFixed(1)} ${speedLabel()} (${v_cc_kmh.toFixed(0)} km/h)`;
    ctx.fillText(line1, x_cc, y0 - 7);
    ctx.restore();
  }
}

function drawMinSinkMarker(ranges) {
  const coeffs = state.activeCoeffs;
  const v_ms_kmh = minSinkSpeed(coeffs);
  if (v_ms_kmh < ranges.v_min_kmh || v_ms_kmh > ranges.v_max_kmh) return;

  const w_ms = polarSink(coeffs, v_ms_kmh) + state.airmass_ms;
  const v_ms_disp = convertSpeed(v_ms_kmh, state.speedUnit);
  const px = toCanvasX(v_ms_disp, ranges);
  const py = toCanvasY(convertRate(w_ms, state.sinkUnit), ranges);

  const { left, right, top, bottom } = chartArea();
  if (px < left || px > right || py < top || py > bottom) return;

  ctx.save();

  // Faint vertical tick to X axis
  ctx.strokeStyle = C.minsink;
  ctx.lineWidth = 1;
  ctx.setLineDash([3, 3]);
  ctx.globalAlpha = 0.5;
  ctx.beginPath();
  ctx.moveTo(px, py);
  ctx.lineTo(px, bottom);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.globalAlpha = 1;

  // Circle on the polar curve
  ctx.beginPath();
  ctx.arc(px, py, 4, 0, Math.PI * 2);
  ctx.strokeStyle = C.minsink;
  ctx.lineWidth = 2;
  ctx.stroke();

  // Speed label below the circle
  ctx.font = '10px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.fillStyle = C.minsink;
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  ctx.fillText('Min sink', px, py + 7);
  ctx.fillText(`${v_ms_disp.toFixed(1)} ${speedLabel()}`, px, py + 19);

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
  if (!state.activeCoeffs || state.canvasW === 0) return;

  ctx.clearRect(0, 0, state.canvasW, state.canvasH);

  const entry  = state.polars[state.selectedIndex];
  const ranges = computeRanges(entry, state.activeCoeffs);
  state.ranges = ranges;

  drawGrid(ranges);
  drawAxes(ranges);

  const primaryStall = computeStallSpeed(entry, state.ballast_kg);

  if (!state.compareMode) drawStallMarker(ranges);

  drawPolarCurve(ranges, state.activeCoeffs, C.polar,
                 primaryStall, entry.v_max_kmh,
                 state.compareMode ? entry.name : null);

  if (state.compareMode && state.compareActiveCoeffs) {
    const cEntry = state.polars[state.compareIndex];
    const cStall = computeStallSpeed(cEntry, state.compareBallast_kg);
    drawPolarCurve(ranges, state.compareActiveCoeffs, C.compare,
                   cStall, cEntry.v_max_kmh, cEntry.name);
  }

  drawMcLine(ranges, state.activeCoeffs, C.mc);
  if (state.compareMode && state.compareActiveCoeffs) {
    drawMcLine(ranges, state.compareActiveCoeffs, C.compare, true);
  }
  if (!state.compareMode) drawMinSinkMarker(ranges);
  if (state.hoverPoint) drawHoverMarker(ranges);
}

// === TOOLTIP ===

const tooltip       = document.getElementById('tooltip');
const tipSpeed      = document.getElementById('tip-speed');
const tipSink       = document.getElementById('tip-sink');
const tipCompareSink = document.getElementById('tip-compare-sink');
const tipLd         = document.getElementById('tip-ld');

function showTooltip(px, py, v_kmh, w_ms) {
  const rateDisp = convertRate(w_ms, state.sinkUnit);
  const sign = rateDisp >= 0 ? '+' : '';
  tipSpeed.textContent = `${convertSpeed(v_kmh, state.speedUnit).toFixed(1)} ${speedLabel()}`;
  tipSink.textContent  = `${sign}${rateDisp.toFixed(2)} ${sinkLabel()}`;

  if (state.compareMode && state.compareActiveCoeffs) {
    const cw_ms = polarSink(state.compareActiveCoeffs, v_kmh) + state.airmass_ms;
    const cRateDisp = convertRate(cw_ms, state.sinkUnit);
    const cSign = cRateDisp >= 0 ? '+' : '';
    tipCompareSink.textContent = `${cSign}${cRateDisp.toFixed(2)} ${sinkLabel()}`;
    tipCompareSink.removeAttribute('hidden');
    tooltip.classList.add('comparing');
  } else {
    tipCompareSink.setAttribute('hidden', '');
    tooltip.classList.remove('comparing');
  }

  // L/D: forward speed (m/s) divided by sink speed (m/s), only meaningful when sinking
  const v_ms = v_kmh / 3.6;
  const ld = w_ms < -0.01 ? (v_ms / Math.abs(w_ms)).toFixed(1) : '—';
  tipLd.textContent = `L/D: ${ld}`;

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
  if (!state.activeCoeffs || !state.ranges) return;

  const rect = canvas.getBoundingClientRect();
  const px = e.clientX - rect.left;
  const py = e.clientY - rect.top;
  const { left, right, top, bottom } = chartArea();

  if (px < left || px > right || py < top || py > bottom) {
    if (state.hoverPoint) { state.hoverPoint = null; redraw(); hideTooltip(); }
    return;
  }

  // Convert pointer X → speed → find sink on polar at that speed.
  // Clamp to [stall_kmh, v_max_kmh] so hover follows only the drawn curve.
  const speedDisp = fromCanvasX(px, state.ranges);
  const v_kmh = Math.max(
    state.ranges.stall_kmh,
    Math.min(state.ranges.v_max_kmh, speedDisp / SPEED_UNITS[state.speedUnit].factor)
  );
  const w_ms = polarSink(state.activeCoeffs, v_kmh) + state.airmass_ms;

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
    state.ballast_kg = 0;
    updateActiveCoeffs();
    updateBallastSlider();
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

  document.getElementById('ballast-slider').addEventListener('input', function () {
    state.ballast_kg = parseInt(this.value, 10);
    document.getElementById('ballast-value').textContent = state.ballast_kg;
    updateActiveCoeffs();
    redraw();
  });

  updateBallastSlider();
  updateSliderConfigs();
  updateSliderLabels();
  initCompareControls();
}

function initCompareControls() {
  const toggle      = document.getElementById('compare-toggle');
  const group       = document.getElementById('compare-aircraft-group');
  const selectEl    = document.getElementById('compare-aircraft-select');
  const ballastGrp  = document.getElementById('compare-ballast-group');
  const ballastSldr = document.getElementById('compare-ballast-slider');
  const ballastVal  = document.getElementById('compare-ballast-value');
  const ballastMax  = document.getElementById('compare-ballast-max');

  selectEl.innerHTML = state.polars
    .map((p, i) => `<option value="${i}">${p.name}</option>`)
    .join('');
  state.compareIndex = state.polars.length > 1 ? 1 : 0;
  selectEl.value = String(state.compareIndex);
  selectEl.removeAttribute('disabled');

  function loadCompareGlider() {
    const cEntry = state.polars[state.compareIndex];
    state.compareBallast_kg = 0;
    state.compareCoeffs = fitPolar(cEntry);
    updateCompareCoeffs();
    const max = cEntry.max_ballast || 0;
    ballastSldr.max      = max;
    ballastSldr.value    = 0;
    ballastSldr.disabled = max === 0;
    ballastVal.textContent = '0';
    ballastMax.textContent = max > 0 ? `/ ${max} kg` : '(none)';
  }

  toggle.addEventListener('click', () => {
    state.compareMode = !state.compareMode;
    toggle.setAttribute('aria-pressed', String(state.compareMode));
    group.hidden      = !state.compareMode;
    ballastGrp.hidden = !state.compareMode;
    if (state.compareMode) {
      loadCompareGlider();
    } else {
      state.compareCoeffs       = null;
      state.compareActiveCoeffs = null;
    }
    state.hoverPoint = null;
    hideTooltip();
    redraw();
  });

  selectEl.addEventListener('change', () => {
    state.compareIndex = parseInt(selectEl.value, 10);
    loadCompareGlider();
    redraw();
  });

  ballastSldr.addEventListener('input', function () {
    state.compareBallast_kg = parseInt(this.value, 10);
    ballastVal.textContent = state.compareBallast_kg;
    updateCompareCoeffs();
    redraw();
  });
}

function updateBallastSlider() {
  const entry  = state.polars[state.selectedIndex];
  const max    = entry.max_ballast || 0;
  const slider = document.getElementById('ballast-slider');
  slider.max   = max;
  slider.value = 0;
  slider.disabled = max === 0;
  document.getElementById('ballast-value').textContent = '0';
  document.getElementById('ballast-max').textContent   = max > 0 ? `/ ${max} kg` : '(none)';
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
    updateActiveCoeffs();

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
