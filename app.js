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
  // Ranges must be multiples of step so 0 sits on the slider's step grid
  ms:  { step: 0.2,  airmassMin: -5.6,   airmassMax: 5.6,  mcMax: 5.6,  decimals: 1 },
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
  surface: '#ffffff',
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

const POLAR_CACHE_KEY = 'polarStore.v1';
const POLAR_CACHE_TTL = 7 * 24 * 3600 * 1000; // 1 week

function readPolarCache(maxAge) {
  try {
    const raw = localStorage.getItem(POLAR_CACHE_KEY);
    if (!raw) return null;
    const { ts, polars } = JSON.parse(raw);
    if (!Array.isArray(polars) || polars.length === 0) return null;
    if (maxAge !== Infinity && Date.now() - ts > maxAge) return null;
    return polars;
  } catch {
    return null;
  }
}

function writePolarCache(polars) {
  try {
    localStorage.setItem(POLAR_CACHE_KEY, JSON.stringify({ ts: Date.now(), polars }));
  } catch { /* storage full or unavailable — cache is best-effort */ }
}

async function fetchPolars() {
  const resp = await fetch(POLAR_URL);
  if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
  const text = await resp.text();
  return parsePolars(text);
}

// Fresh cache → instant load. Otherwise fetch and cache; if the fetch
// fails, fall back to a stale cache rather than showing an error.
async function loadPolars() {
  const fresh = readPolarCache(POLAR_CACHE_TTL);
  if (fresh) return fresh;
  try {
    const polars = await fetchPolars();
    if (polars.length > 0) writePolarCache(polars);
    return polars;
  } catch (err) {
    const stale = readPolarCache(Infinity);
    if (stale) {
      console.warn('Using stale polar cache:', err);
      return stale;
    }
    throw err;
  }
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
  // Y axis top: curve peak or 15% of range for breathing room. The MC anchor
  // point at (0, mc) is deliberately NOT enveloped — at high MC settings it
  // would push the polar into the bottom half of the chart; the dashed MC
  // line simply enters from the top edge instead.
  const totalNeg = Math.abs(w_min_disp);
  const curve_top_disp = convertRate(w_max_ms, state.sinkUnit);
  const w_max_disp = Math.max(totalNeg * 0.15, curve_top_disp * 1.15);

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
  apply('surface', '--color-surface');
  apply('zeroline','--color-zeroline');
  apply('minsink', '--color-minsink');
  apply('stall',   '--color-stall');
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

function drawPolarCurve(ranges, coeffs, color, stallKmh, curveVmax) {
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

  for (let i = 0; i <= steps; i++) {
    const v_kmh = v_start + i * dv;
    const w_ms  = polarSink(coeffs, v_kmh) + state.airmass_ms;

    const px = toCanvasX(convertSpeed(v_kmh, state.speedUnit), ranges);
    const py = toCanvasY(convertRate(w_ms,   state.sinkUnit),  ranges);

    if (!started) { ctx.moveTo(px, py); started = true; }
    else ctx.lineTo(px, py);
  }
  ctx.stroke();

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
  // Primary: label to the right of tangent point, above.
  // Compare: label to the left of tangent point, below — avoids overlap.
  const label = state.mc_ms < 1e-9 ? 'Best Glide' : 'STF';
  const spd = v_opt_disp.toFixed(1);
  ctx.save();
  ctx.font = 'bold 11px -apple-system, BlinkMacSystemFont, sans-serif';
  ctx.fillStyle = color;
  ctx.textBaseline = 'top';
  if (labelBelow) {
    ctx.textAlign = 'right';
    const labelX = Math.max(px - 8, MARGIN.left + 4);
    ctx.fillText(`${label}: ${spd} ${speedLabel()}`, labelX, py + 8);
  } else {
    ctx.textAlign = 'left';
    const labelX = Math.min(px + 8, state.canvasW - MARGIN.right - 4);
    ctx.fillText(`${label}: ${spd} ${speedLabel()}`, labelX, py - 18);
  }
  ctx.restore();

  // Cross-country speed: where the MC line crosses y=0 (zero-rate axis).
  // From rate(v) = mc_disp + slope*v = 0  →  v_cc = -mc_disp / slope
  // At MC=0 the crossing is the origin and carries no information — skip it.
  if (state.mc_ms < 1e-9) return;
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
    // Speed label: primary anchors left (text extends right), above zero line.
    // Compare anchors right (text extends left), below zero line.
    ctx.font = '10px -apple-system, BlinkMacSystemFont, sans-serif';
    ctx.fillStyle = color;
    const v_cc_kmh = v_cc_disp / SPEED_UNITS[state.speedUnit].factor;
    const ccLabel = state.speedUnit === 'kmh'
      ? `${v_cc_kmh.toFixed(0)} km/h`
      : `${v_cc_disp.toFixed(1)} ${speedLabel()} (${v_cc_kmh.toFixed(0)} km/h)`;
    if (labelBelow) {
      ctx.textAlign = 'right';
      ctx.textBaseline = 'top';
      ctx.fillText(ccLabel, x_cc - 8, y0 + 8);
    } else {
      ctx.textAlign = 'left';
      ctx.textBaseline = 'bottom';
      ctx.fillText(ccLabel, x_cc + 8, y0 - 8);
    }
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

// Best glide ratio in still air: tangent from origin → v = sqrt(c/a)
function bestGlideRatio(coeffs) {
  const disc = coeffs.c / coeffs.a;
  if (disc <= 0) return null;
  const v_kmh = Math.sqrt(disc);
  const w_ms = polarSink(coeffs, v_kmh);
  if (w_ms >= -1e-9) return null;
  return (v_kmh / 3.6) / Math.abs(w_ms);
}

// Fixed legend box (top-right) for compare mode — replaces fragile labels at
// the curve tails, which sat against the chart edge and could be clipped.
function drawLegend() {
  const items = [
    { color: C.polar,   entry: state.polars[state.selectedIndex], coeffs: state.activeCoeffs },
    { color: C.compare, entry: state.polars[state.compareIndex],  coeffs: state.compareActiveCoeffs },
  ];

  const { right, top } = chartArea();
  const pad = 10, rowH = 18, swatchW = 18, gap = 7;

  ctx.save();
  ctx.font = '12px -apple-system, BlinkMacSystemFont, sans-serif';

  const texts = items.map(it => {
    const ld = it.coeffs ? bestGlideRatio(it.coeffs) : null;
    return ld ? `${it.entry.name} — L/D ${Math.round(ld)}` : it.entry.name;
  });
  const textW = Math.max(...texts.map(t => ctx.measureText(t).width));
  const boxW = pad + swatchW + gap + textW + pad;
  const boxH = pad * 2 + rowH * items.length - 4;
  const x = right - boxW - 8;
  const y = top + 8;

  ctx.globalAlpha = 0.92;
  ctx.fillStyle = C.surface;
  ctx.strokeStyle = C.grid;
  ctx.lineWidth = 1;
  ctx.beginPath();
  if (ctx.roundRect) ctx.roundRect(x, y, boxW, boxH, 6);
  else ctx.rect(x, y, boxW, boxH);
  ctx.fill();
  ctx.stroke();
  ctx.globalAlpha = 1;

  items.forEach((it, i) => {
    const cy = y + pad + rowH * i + 5;
    ctx.strokeStyle = it.color;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    ctx.moveTo(x + pad, cy);
    ctx.lineTo(x + pad + swatchW, cy);
    ctx.stroke();

    ctx.fillStyle = C.text;
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    ctx.fillText(texts[i], x + pad + swatchW + gap, cy);
  });

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
                 primaryStall, entry.v_max_kmh);

  if (state.compareMode && state.compareActiveCoeffs) {
    const cEntry = state.polars[state.compareIndex];
    const cStall = computeStallSpeed(cEntry, state.compareBallast_kg);
    drawPolarCurve(ranges, state.compareActiveCoeffs, C.compare,
                   cStall, cEntry.v_max_kmh);
  }

  drawMcLine(ranges, state.activeCoeffs, state.compareMode ? C.polar : C.mc);
  if (state.compareMode && state.compareActiveCoeffs) {
    drawMcLine(ranges, state.compareActiveCoeffs, C.compare, true);
  }
  if (!state.compareMode) drawMinSinkMarker(ranges);
  if (state.compareMode) drawLegend();
  if (state.hoverPoint) drawHoverMarker(ranges);
}

// === TOOLTIP ===

const tooltip       = document.getElementById('tooltip');
const tipSpeed      = document.getElementById('tip-speed');
const tipSink       = document.getElementById('tip-sink');
const tipCompareSink = document.getElementById('tip-compare-sink');
const tipLd          = document.getElementById('tip-ld');
const tipCompareLd   = document.getElementById('tip-compare-ld');

function showTooltip(px, py, v_kmh, w_ms) {
  const rateDisp = convertRate(w_ms, state.sinkUnit);
  const sign = rateDisp >= 0 ? '+' : '';
  tipSpeed.textContent = `${convertSpeed(v_kmh, state.speedUnit).toFixed(1)} ${speedLabel()}`;
  tipSink.textContent  = `${sign}${rateDisp.toFixed(2)} ${sinkLabel()}`;

  const v_ms = v_kmh / 3.6;

  if (state.compareMode && state.compareActiveCoeffs) {
    const cw_ms = polarSink(state.compareActiveCoeffs, v_kmh) + state.airmass_ms;
    const cRateDisp = convertRate(cw_ms, state.sinkUnit);
    const cSign = cRateDisp >= 0 ? '+' : '';
    tipCompareSink.textContent = `${cSign}${cRateDisp.toFixed(2)} ${sinkLabel()}`;
    tipCompareSink.removeAttribute('hidden');
    const cLd = cw_ms < -0.01 ? (v_ms / Math.abs(cw_ms)).toFixed(1) : '—';
    tipCompareLd.textContent = `L/D: ${cLd}`;
    tipCompareLd.removeAttribute('hidden');
    tooltip.classList.add('comparing');
  } else {
    tipCompareSink.setAttribute('hidden', '');
    tipCompareLd.setAttribute('hidden', '');
    tooltip.classList.remove('comparing');
  }

  // L/D: forward speed (m/s) divided by sink speed (m/s), only meaningful when sinking
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

// Touch: tap-to-show/hide, drag-to-track. Mouse: hover as normal.
const DRAG_THRESHOLD = 8;
let touchStartPos = null;
let touchDragging  = false;

canvas.addEventListener('pointerdown', (e) => {
  if (e.pointerType === 'mouse') { handlePointer(e); return; }
  touchStartPos = { x: e.clientX, y: e.clientY };
  touchDragging  = false;
  canvas.setPointerCapture(e.pointerId);
});

canvas.addEventListener('pointermove', (e) => {
  if (e.pointerType === 'mouse') { handlePointer(e); return; }
  if (!touchStartPos) return;
  const dx = e.clientX - touchStartPos.x;
  const dy = e.clientY - touchStartPos.y;
  if (!touchDragging && Math.hypot(dx, dy) > DRAG_THRESHOLD) touchDragging = true;
  if (touchDragging) handlePointer(e);
});

canvas.addEventListener('pointerup', (e) => {
  if (e.pointerType !== 'touch') return;
  if (!touchDragging) {
    // Tap: toggle visibility
    if (state.hoverPoint) {
      state.hoverPoint = null; redraw(); hideTooltip();
    } else {
      handlePointer(e);
    }
  }
  touchStartPos = null;
  touchDragging  = false;
});

canvas.addEventListener('pointercancel', () => {
  touchStartPos = null; touchDragging = false;
});

canvas.addEventListener('pointerleave', (e) => {
  if (e.pointerType === 'touch') return; // touch manages its own lifecycle
  state.hoverPoint = null; redraw(); hideTooltip();
});

// === UI CONTROLS ===

function updateSliderLabels() {
  const airmassSlider = document.getElementById('airmass-slider');
  const mcSlider      = document.getElementById('mc-slider');
  const airmassDisp = parseFloat(airmassSlider.value);
  const mcDisp      = parseFloat(mcSlider.value);
  const dec = SINK_SLIDER[state.sinkUnit].decimals;

  const airmassText = (airmassDisp >= 0 ? '+' : '') + airmassDisp.toFixed(dec);
  const mcText      = mcDisp.toFixed(dec);
  document.getElementById('airmass-value').textContent = airmassText;
  document.getElementById('mc-value').textContent = mcText;
  airmassSlider.setAttribute('aria-valuetext', `${airmassText} ${sinkLabel()}`);
  mcSlider.setAttribute('aria-valuetext', `${mcText} ${sinkLabel()}`);

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

  // The browser snaps slider values to the step grid — read them back so
  // the chart always matches what the labels show.
  state.airmass_ms = sinkDispToMs(parseFloat(airmassSlider.value));
  state.mc_ms      = sinkDispToMs(parseFloat(mcSlider.value));
}

// Glider pickers are text inputs backed by a shared <datalist>, so the long
// list (hundreds of gliders) can be filtered by typing.
function gliderIndexByName(text) {
  const needle = text.trim().toLowerCase();
  return state.polars.findIndex(p => p.name.toLowerCase() === needle);
}

function populateGliderList() {
  const dl = document.getElementById('glider-list');
  dl.replaceChildren(...state.polars.map(p => {
    const o = document.createElement('option');
    o.value = p.name;
    return o;
  }));
}

// Wire a datalist-backed picker: commit on a valid name, revert otherwise.
function initGliderPicker(input, getIndex, onPick) {
  input.removeAttribute('disabled');
  input.value = state.polars[getIndex()].name;
  input.addEventListener('focus', () => input.select());
  const revert = () => { input.value = state.polars[getIndex()].name; };
  input.addEventListener('change', () => {
    const idx = gliderIndexByName(input.value);
    if (idx >= 0 && idx !== getIndex()) onPick(idx);
    revert();
  });
  input.addEventListener('blur', revert);
}

function initControls() {
  populateGliderList();

  initGliderPicker(
    document.getElementById('aircraft-select'),
    () => state.selectedIndex,
    (idx) => {
      state.selectedIndex = idx;
      state.coeffs = fitPolar(state.polars[idx]);
      state.ballast_kg = 0;
      updateActiveCoeffs();
      updateBallastSlider();
      state.hoverPoint = null;
      hideTooltip();
      redraw();
      scheduleHashUpdate();
    }
  );

  document.querySelectorAll('input[name="speed-unit"]').forEach(r => {
    r.checked = r.value === state.speedUnit;
    r.addEventListener('change', () => {
      if (!r.checked) return;
      state.speedUnit = r.value;
      redraw();
      scheduleHashUpdate();
    });
  });

  document.querySelectorAll('input[name="sink-unit"]').forEach(r => {
    r.checked = r.value === state.sinkUnit;
    r.addEventListener('change', () => {
      if (!r.checked) return;
      state.sinkUnit = r.value;
      updateSliderConfigs();
      updateSliderLabels();
      redraw();
      scheduleHashUpdate();
    });
  });

  document.getElementById('airmass-slider').addEventListener('input', function () {
    state.airmass_ms = sinkDispToMs(parseFloat(this.value));
    updateSliderLabels();
    redraw();
    scheduleHashUpdate();
  });

  document.getElementById('mc-slider').addEventListener('input', function () {
    state.mc_ms = sinkDispToMs(parseFloat(this.value));
    updateSliderLabels();
    redraw();
    scheduleHashUpdate();
  });

  document.getElementById('ballast-slider').addEventListener('input', function () {
    state.ballast_kg = parseInt(this.value, 10);
    document.getElementById('ballast-value').textContent = state.ballast_kg;
    this.setAttribute('aria-valuetext', `${state.ballast_kg} kg of ${this.max} kg`);
    updateActiveCoeffs();
    redraw();
    scheduleHashUpdate();
  });

  updateBallastSlider();
  updateSliderConfigs();
  updateSliderLabels();
  initCompareControls();
}

function initCompareControls() {
  const toggle      = document.getElementById('compare-toggle');
  const group       = document.getElementById('compare-aircraft-group');
  const inputEl     = document.getElementById('compare-aircraft-select');
  const ballastGrp  = document.getElementById('compare-ballast-group');
  const ballastSldr = document.getElementById('compare-ballast-slider');
  const ballastVal  = document.getElementById('compare-ballast-value');
  const ballastMax  = document.getElementById('compare-ballast-max');

  // Default compare glider, unless one was already restored from the URL
  if (!state.compareMode) {
    state.compareIndex = state.polars.length > 1 ? 1 : 0;
  }

  // Preserves state.compareBallast_kg (clamped to the new glider's max)
  function loadCompareGlider() {
    const cEntry = state.polars[state.compareIndex];
    const max = cEntry.max_ballast || 0;
    state.compareBallast_kg = Math.min(state.compareBallast_kg, max);
    state.compareCoeffs = fitPolar(cEntry);
    updateCompareCoeffs();
    ballastSldr.max      = max;
    ballastSldr.value    = state.compareBallast_kg;
    ballastSldr.disabled = max === 0;
    ballastSldr.setAttribute('aria-valuetext', `${state.compareBallast_kg} kg of ${max} kg`);
    ballastVal.textContent = String(state.compareBallast_kg);
    ballastMax.textContent = max > 0 ? `/ ${max} kg` : '(none)';
  }

  initGliderPicker(
    inputEl,
    () => state.compareIndex,
    (idx) => {
      state.compareIndex = idx;
      state.compareBallast_kg = 0;
      loadCompareGlider();
      redraw();
      scheduleHashUpdate();
    }
  );

  function setCompareVisible(on) {
    toggle.setAttribute('aria-pressed', String(on));
    group.hidden      = !on;
    ballastGrp.hidden = !on;
  }

  toggle.addEventListener('click', () => {
    state.compareMode = !state.compareMode;
    setCompareVisible(state.compareMode);
    if (state.compareMode) {
      loadCompareGlider();
    } else {
      state.compareCoeffs       = null;
      state.compareActiveCoeffs = null;
    }
    state.hoverPoint = null;
    hideTooltip();
    redraw();
    scheduleHashUpdate();
  });

  ballastSldr.addEventListener('input', function () {
    state.compareBallast_kg = parseInt(this.value, 10);
    ballastVal.textContent = state.compareBallast_kg;
    this.setAttribute('aria-valuetext', `${state.compareBallast_kg} kg of ${this.max} kg`);
    updateCompareCoeffs();
    redraw();
    scheduleHashUpdate();
  });

  // Restore compare mode from the URL hash
  if (state.compareMode) {
    setCompareVisible(true);
    loadCompareGlider();
  }
}

// Preserves state.ballast_kg (clamped to the glider's max) so restored URL
// state survives; callers that want a reset zero state.ballast_kg first.
function updateBallastSlider() {
  const entry  = state.polars[state.selectedIndex];
  const max    = entry.max_ballast || 0;
  const slider = document.getElementById('ballast-slider');
  state.ballast_kg = Math.min(state.ballast_kg, max);
  slider.max   = max;
  slider.value = state.ballast_kg;
  slider.disabled = max === 0;
  slider.setAttribute('aria-valuetext', `${state.ballast_kg} kg of ${max} kg`);
  document.getElementById('ballast-value').textContent = String(state.ballast_kg);
  document.getElementById('ballast-max').textContent   = max > 0 ? `/ ${max} kg` : '(none)';
}

// === URL STATE ===
// The full configuration lives in the hash so any view can be shared as a
// link, e.g. #g=Discus&su=kts&zu=ms&mc=1.5&b=50&cmp=1&cg=LS-4

function applyHashState() {
  const params = new URLSearchParams(location.hash.slice(1));

  const su = params.get('su');
  if (su && SPEED_UNITS[su]) state.speedUnit = su;
  const zu = params.get('zu');
  if (zu && SINK_UNITS[zu]) state.sinkUnit = zu;

  const g = params.get('g');
  if (g) {
    const idx = gliderIndexByName(g);
    if (idx >= 0) state.selectedIndex = idx;
  }

  const cfg = SINK_SLIDER[state.sinkUnit];
  const am = parseFloat(params.get('am'));
  if (Number.isFinite(am)) {
    state.airmass_ms = Math.max(sinkDispToMs(cfg.airmassMin),
                                Math.min(sinkDispToMs(cfg.airmassMax), am));
  }
  const mc = parseFloat(params.get('mc'));
  if (Number.isFinite(mc)) {
    state.mc_ms = Math.max(0, Math.min(sinkDispToMs(cfg.mcMax), mc));
  }

  const b = parseInt(params.get('b'), 10);
  if (Number.isFinite(b)) {
    const max = state.polars[state.selectedIndex].max_ballast || 0;
    state.ballast_kg = Math.max(0, Math.min(max, b));
  }

  if (params.get('cmp') === '1') {
    state.compareMode = true;
    const cg = params.get('cg');
    if (cg) {
      const idx = gliderIndexByName(cg);
      if (idx >= 0) state.compareIndex = idx;
    }
    const cb = parseInt(params.get('cb'), 10);
    if (Number.isFinite(cb)) {
      const max = state.polars[state.compareIndex].max_ballast || 0;
      state.compareBallast_kg = Math.max(0, Math.min(max, cb));
    }
  }
}

function updateHash() {
  const p = new URLSearchParams();
  p.set('g', state.polars[state.selectedIndex].name);
  p.set('su', state.speedUnit);
  p.set('zu', state.sinkUnit);
  if (state.airmass_ms !== 0) p.set('am', state.airmass_ms.toFixed(3));
  if (state.mc_ms > 0)        p.set('mc', state.mc_ms.toFixed(3));
  if (state.ballast_kg > 0)   p.set('b', String(state.ballast_kg));
  if (state.compareMode) {
    p.set('cmp', '1');
    p.set('cg', state.polars[state.compareIndex].name);
    if (state.compareBallast_kg > 0) p.set('cb', String(state.compareBallast_kg));
  }
  history.replaceState(null, '', '#' + p.toString());
}

// === INITIALISATION ===

function debounce(fn, ms) {
  let timer;
  return (...args) => { clearTimeout(timer); timer = setTimeout(() => fn(...args), ms); };
}

const scheduleHashUpdate = debounce(updateHash, 250);

const chartStatus = document.getElementById('chart-status');

function setStatus(text) {
  chartStatus.textContent = text;
  chartStatus.removeAttribute('hidden');
}

function clearStatus() {
  chartStatus.setAttribute('hidden', '');
}

async function init() {
  resolveColours();
  resizeCanvas();

  // Redraw whenever the chart container changes size — not just on window
  // resize: toggling compare mode grows the control bar and shrinks the chart.
  const ro = new ResizeObserver(debounce(() => { resizeCanvas(); redraw(); }, 100));
  ro.observe(canvas.parentElement);

  window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
    resolveColours(); redraw();
  });

  setStatus('Loading polar data…');
  try {
    state.polars = await loadPolars();
    if (state.polars.length === 0) throw new Error('No polars parsed');

    const defaultIdx = state.polars.findIndex(p => p.name.startsWith('Discus'));
    state.selectedIndex = defaultIdx >= 0 ? defaultIdx : 0;

    applyHashState();

    state.coeffs = fitPolar(state.polars[state.selectedIndex]);
    updateActiveCoeffs();

    clearStatus();
    initControls();
    redraw();
  } catch (err) {
    console.error('Failed to load polars:', err);
    setStatus('Failed to load polar data — check your connection and reload the page.');
    const input = document.getElementById('aircraft-select');
    input.placeholder = 'Error loading polars';
    input.disabled = true;
  }
}

init();
