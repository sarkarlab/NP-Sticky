// npsticky.js — client-side port of NP-Sticky CGI back-end logic
// All thermodynamic models are faithful ports of the original Python CGI scripts.
// Licensed under CC BY-NC 4.0 (same as original).

'use strict';

// ─────────────────────────────────────────────
// DNA sequence helpers
// ─────────────────────────────────────────────

function basesToNum(base) {
  const lbase = base.toLowerCase();
  const bUpNum = [], bLowNum = [];
  let flag = true;
  for (const ch of lbase) {
    if      (ch === 'a') { bUpNum.push(0); bLowNum.push(1); }
    else if (ch === 't') { bUpNum.push(1); bLowNum.push(0); }
    else if (ch === 'c') { bUpNum.push(2); bLowNum.push(3); }
    else if (ch === 'g') { bUpNum.push(3); bLowNum.push(2); }
    else { flag = false; break; }
  }
  bLowNum.reverse();
  return { bUpNum, bLowNum, flag };
}

function numToBase(nums) {
  return nums.map(n => ['a','t','c','g'][n]).join('');
}

// ─────────────────────────────────────────────
// Thermodynamic scoring matrices
// ─────────────────────────────────────────────

// HSCORE[row][col] where row = a*4+b, col = c*4+d
const HSCORE = [
  [ 2.28,  0.69,  2.64,  1.89,  0.61, -1.00,  0.88,  0.14,  2.64,  1.33,  3.01,  2.26,  1.89,  0.74,  2.26,  1.50],
  [ 0.61,  2.23,  2.52,  1.63, -0.88,  0.69,  0.73,  0.07,  0.77,  2.59,  2.88,  2.00,  0.02,  1.84,  2.13,  1.25],
  [ 2.64,  2.52,  2.77,  0.17,  0.77,  0.64,  1.33, -1.44,  3.01,  2.88,  3.14,  0.47,  2.26,  2.13,  2.38, -0.52],
  [ 1.89,  1.74,  0.43,  1.58,  0.02,  0.71, -1.28, -0.13,  2.26,  2.10,  0.79,  1.94,  1.51,  1.35,  0.11,  1.19],
  [ 0.69, -0.58,  0.92,  0.42,  2.23,  0.68,  2.59,  1.84,  2.52,  0.97,  2.88,  2.13,  1.74,  0.43,  2.10,  1.35],
  [-1.00,  0.68,  0.75,  0.34,  0.69,  2.18,  2.47,  1.58,  0.64,  2.47,  2.76,  1.87,  0.71,  1.69,  1.98,  0.74],
  [ 1.33,  0.97,  1.05, -1.30,  2.59,  2.47,  2.72,  0.45,  2.88,  2.75,  3.01,  0.62,  2.10,  1.97,  2.23,  0.08],
  [ 0.74,  0.43, -1.45,  0.44,  1.84,  1.69, -0.12,  1.53,  2.13,  1.97,  0.62,  1.81,  1.35,  0.52, -0.47,  1.03],
  [ 2.64,  0.92,  3.01,  2.26,  2.52,  0.75,  2.88,  2.13,  2.77,  1.05,  3.14,  2.38,  0.43, -1.45,  0.75,  0.03],
  [ 0.88,  2.59,  2.88,  2.00,  0.73,  2.47,  2.76,  1.87,  1.33,  2.72,  3.01,  2.13, -1.28, -0.12,  0.40, -0.32],
  [ 3.01,  2.88,  3.14,  0.81,  2.88,  2.76,  3.01,  0.98,  3.14,  3.01,  3.26,  0.79,  0.79,  0.62,  0.70, -1.84],
  [ 2.26,  2.10,  0.75,  1.94,  2.13,  1.98,  0.40,  1.82,  2.38,  2.23,  0.70,  2.07,  0.11, -0.47, -2.17, -0.11],
  [ 1.89,  0.42,  2.26,  1.50,  1.63,  0.34,  2.00,  1.24,  0.17, -1.30,  0.81, -0.25,  1.58,  0.44,  1.94,  1.19],
  [ 0.14,  1.84,  2.13,  1.24,  0.07,  1.58,  1.87,  1.15, -1.44,  0.45,  0.98, -0.59, -0.13,  1.53,  1.82,  0.93],
  [ 2.26,  2.13,  2.38, -0.25,  2.00,  1.87,  2.13, -0.59,  0.47,  0.62,  0.79, -2.24,  1.94,  1.81,  2.07, -1.11],
  [ 1.50,  1.35,  0.03,  1.19,  1.25,  0.74, -0.32,  0.93, -0.52,  0.08, -1.84, -1.11,  1.19,  1.03, -0.11,  0.87],
];

const STACK_FORIAN = [
  [-0.8, -0.5, -0.3, -0.9],
  [-0.5, -0.4, -0.2, -0.5],
  [-0.3, -0.2, -0.6, -0.6],
  [-0.9, -0.5, -0.6, -1.1],
];

const STACK_PROTOZANOVA = [
  [-1.11, -1.34, -1.81, -1.06],
  [-0.19, -1.11, -1.43, -0.55],
  [-0.55, -1.06, -1.44, -0.91],
  [-1.43, -1.81, -2.17, -1.44],
];

// ─────────────────────────────────────────────
// Core thermodynamic functions
// ─────────────────────────────────────────────

function hBond(bUpNum, bLowNum) {
  let delGHbond = 0;
  const n = Math.min(bUpNum.length, bLowNum.length);
  for (let i = 0; i < n - 1; i++) {
    const a = bUpNum[i];
    const b = bUpNum[i + 1];
    const c = bLowNum[n - i - 1];
    const d = bLowNum[n - i - 2];
    delGHbond += HSCORE[a * 4 + b][c * 4 + d];
  }
  return delGHbond;
}

// stackOpt: 'forian' | 'protozanova' | 'off'
// For ligation calculators that don't expose stacking option, always pass 'forian'
function stackEnergy(stackOpt, left1, left2, right1, right2) {
  let S;
  if      (stackOpt === 'forian')      S = STACK_FORIAN;
  else if (stackOpt === 'protozanova') S = STACK_PROTOZANOVA;
  else return 0; // 'off'
  return S[left1][left2] + S[right1][right2];
}

function delGCal(seUp, seLow, seLeftStr, seRightStr, seType, stackOpt, delGLig) {
  const { bUpNum: seLeftUp, bLowNum: seLeftLow } = basesToNum(seLeftStr);
  const { bUpNum: seRightUp, bLowNum: seRightLow } = basesToNum(seRightStr);
  let delGse;
  if (seType === '5') {
    delGse = hBond(seUp, seLow) +
             stackEnergy(stackOpt, seLeftUp[0], seUp[0], seRightLow[0], seLow[0]) +
             delGLig;
  } else { // '3'
    const n = seLow.length - 1;
    delGse = hBond(seUp, seLow) +
             stackEnergy(stackOpt, seLow[n], seLeftLow[0], seUp[n], seRightUp[0]) +
             delGLig;
  }
  return delGse;
}

function delGToK(delG) {
  const T = 25 + 273.15;
  const R = 1.9858775;
  return Math.exp((delG * 1e3) / (R * T)) * 1e9; // nM
}

function massToMolar(massConc, bp, vol, totalVol) {
  const MW = bp * 607.4 + 157.9;
  return ((massConc * vol) / MW) / (totalVol * 1e-6);
}

function massToMolarStock(massConc, bp) {
  const MW = bp * 607.4 + 157.9;
  return (massConc * 0.001) / MW * 1e9;
}

function KcirCal(K, bp) {
  const MW = bp * 607.4 + 157.9;
  const j = 3.6 * Math.pow(30.8e6 / MW, 1.5);
  const Aeff = j / 6.022;
  return K / Aeff;
}

// ─────────────────────────────────────────────
// Linear algebra
// ─────────────────────────────────────────────

function gaussianElim(A_orig, b_orig) {
  const n = b_orig.length;
  // Build augmented matrix
  const M = A_orig.map((row, i) => row.map(v => v).concat([b_orig[i]]));
  for (let col = 0; col < n; col++) {
    // Partial pivoting
    let maxRow = col, maxVal = Math.abs(M[col][col]);
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(M[row][col]) > maxVal) { maxVal = Math.abs(M[row][col]); maxRow = row; }
    }
    [M[col], M[maxRow]] = [M[maxRow], M[col]];
    if (Math.abs(M[col][col]) < 1e-14) return null;
    for (let row = col + 1; row < n; row++) {
      const f = M[row][col] / M[col][col];
      for (let k = col; k <= n; k++) M[row][k] -= f * M[col][k];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}

function newtonSolve(F, x0, maxIter = 300, tol = 1e-8) {
  const n = x0.length;
  const h = 1e-7;
  let x = x0.slice();
  for (let iter = 0; iter < maxIter; iter++) {
    const fx = F(x);
    const rnorm = Math.sqrt(fx.reduce((s, v) => s + v * v, 0));
    if (!isFinite(rnorm)) return { x, ok: false };
    if (rnorm < tol) break;
    // Numerical Jacobian: J[i][j] = dF_i/dx_j
    const J = Array.from({ length: n }, () => new Array(n));
    for (let j = 0; j < n; j++) {
      const xh = x.slice(); xh[j] += h;
      const fxh = F(xh);
      for (let i = 0; i < n; i++) J[i][j] = (fxh[i] - fx[i]) / h;
    }
    const dx = gaussianElim(J, fx.map(v => -v));
    if (!dx || dx.some(v => !isFinite(v))) return { x, ok: false };
    // Backtracking line search
    let alpha = 1.0;
    for (let ls = 0; ls < 8; ls++) {
      const xnew = x.map((xi, i) => xi + alpha * dx[i]);
      const fxnew = F(xnew);
      const rnew = Math.sqrt(fxnew.reduce((s, v) => s + v * v, 0));
      if (rnew < rnorm) { x = xnew; break; }
      alpha *= 0.5;
    }
  }
  const fx = F(x);
  const ok = x.every(v => v >= -1e-6) && fx.every(v => Math.abs(v) < 1e-4);
  return { x, ok };
}

// Solve with multiple random restarts
function solveWithRestarts(F, n, scale, numRestarts = 50) {
  for (let trial = 0; trial < numRestarts; trial++) {
    const x0 = Array.from({ length: n }, () => Math.random() * scale);
    const { x, ok } = newtonSolve(F, x0);
    if (ok) return x;
  }
  return null;
}

// ─────────────────────────────────────────────
// 2-piece linear ligation — analytical solution
// ─────────────────────────────────────────────
// Returns {A, B, AB, incorp, eff} or null

function solveLinear2(K, c1, c2, libIdx /* 1 or 2 */) {
  const s = c1 + c2 + K;
  const disc = s * s - 4 * c1 * c2;
  if (disc < 0) return null;
  const AB = (s - Math.sqrt(disc)) / 2;
  if (AB < 0) return null;
  const A = c1 - AB;
  const B = c2 - AB;
  const incorp = (AB / (libIdx === 1 ? c1 : c2)) * 100;
  const eff    = (AB / (c1 + c2)) * 2 * 100;
  return { conc: [A, B, AB], incorp, eff };
}

// ─────────────────────────────────────────────
// 3-piece linear ligation — Newton's method
// ─────────────────────────────────────────────
// x = [A, B, C, AB, BC, ABC1, ABC2]
// F(x) = [A*B - K12*AB,  B*C - K23*BC,  AB*C - K23*ABC2,
//          BC*A - K12*ABC1,
//          A+AB+ABC1+ABC2-c1, B+AB+BC+ABC1+ABC2-c2, C+BC+ABC1+ABC2-c3]

function solveLinear3(K12, K23, c1, c2, c3, libIdx) {
  const conc = [c1, c2, c3];
  const F = (x) => [
    x[0] * x[1] - K12 * x[3],
    x[1] * x[2] - K23 * x[4],
    x[3] * x[2] - K23 * x[6],
    x[4] * x[0] - K12 * x[5],
    x[0] + x[3] + x[5] + x[6] - c1,
    x[1] + x[3] + x[4] + x[5] + x[6] - c2,
    x[2] + x[4] + x[5] + x[6] - c3,
  ];
  const scale = Math.max(c1, c2, c3);
  const x = solveWithRestarts(F, 7, scale);
  if (!x) return null;
  const incorp = ((x[5] + x[6]) / conc[libIdx - 1]) * 100;
  const eff    = ((x[5] + x[6]) / (c1 + c2 + c3)) * 3 * 100;
  return { conc: x, incorp, eff };
}

// ─────────────────────────────────────────────
// Circular ligation (vector-insert) — Newton's method
// ─────────────────────────────────────────────
// x[0..15] = [V, I, VI1, VI2, VIc1, VIc2, IVI1, IVI2, VIV1, VIV2,
//             VIVI1a, VIVI2a, VIVI1b, VIVI2b, VIVIc1, VIVIc2]
// conc = [V_total, I_total]

function solveCircular(Ki, KcirVI, KcirVIVI, concV, concI) {
  const F = (x) => [
    x[0]*x[1]  - Ki[0]*x[2],
    x[0]*x[1]  - Ki[1]*x[3],
    x[2]       - KcirVI[1]*x[5],
    x[3]       - KcirVI[0]*x[4],
    x[2]*x[1]  - Ki[1]*x[7],
    x[3]*x[1]  - Ki[0]*x[6],
    x[2]*x[0]  - Ki[1]*x[9],
    x[3]*x[0]  - Ki[0]*x[8],
    (x[6]+x[7])*x[0] - Ki[0]*x[10],
    (x[6]+x[7])*x[0] - Ki[1]*x[11],
    (x[8]+x[9])*x[1] - Ki[0]*x[12],
    (x[8]+x[9])*x[1] - Ki[1]*x[13],
    (x[10]+x[12]) - KcirVIVI[1]*x[15],
    (x[11]+x[13]) - KcirVIVI[0]*x[14],
    x[0]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+2*(x[8]+x[9]+x[10]+x[11]+x[12]+x[13]+x[14]+x[15]) - concV,
    x[1]+x[2]+x[3]+x[4]+x[5]+x[8]+x[9]+2*(x[6]+x[7]+x[10]+x[11]+x[12]+x[13]+x[14]+x[15]) - concI,
  ];
  const scale = Math.max(concV, concI);
  const x = solveWithRestarts(F, 16, scale);
  if (!x) return null;
  const incorp = (x[4] + x[5]) / concI * 100;
  const eff    = (x[4] + x[5]) / (concV + concI) * 2 * 100;
  const circTotal = x[4] + x[5] + x[14] + x[15];
  const accur = circTotal > 0 ? (x[4] + x[5]) / circTotal * 100 : 0;
  return { conc: x, incorp, eff, accur };
}

// ─────────────────────────────────────────────
// Recommendation grid computation (30×30)
// ─────────────────────────────────────────────

function logspace(start, stop, num, base) {
  const arr = [];
  const step = (stop - start) / (num - 1);
  for (let i = 0; i < num; i++) arr.push(Math.pow(base, start + i * step));
  return arr;
}

function linspace(start, stop, num) {
  const arr = [];
  const step = (stop - start) / (num - 1);
  for (let i = 0; i < num; i++) arr.push(start + i * step);
  return arr;
}

// 2-piece linear recommendation  (analytical → fast)
function computeLinRecom2(K, libConc, libOpt /* 1 or 2 */) {
  const n = 30;
  const ratio    = logspace(-3, 3, n, 2);
  const libRange = linspace(libConc / 10, libConc, n);
  const incorp = Array.from({ length: n }, () => new Array(n).fill(0));
  const eff    = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      const concArr = libOpt === 1
        ? [libRange[j], ratio[i] * libRange[j]]
        : [ratio[i] * libRange[j], libRange[j]];
      const res = solveLinear2(K, concArr[0], concArr[1], libOpt);
      if (res) { incorp[j][i] = res.incorp; eff[j][i] = res.eff; }
    }
  }
  return { ratio, libRange, incorp, eff };
}

// 3-piece linear recommendation  (Newton's)
function computeLinRecom3(K12, K23, libConc, libOpt /* 1,2,3 */) {
  const n = 30;
  const ratio    = logspace(-3, 3, n, 2);
  const libRange = linspace(libConc / 10, libConc, n);
  const incorp = Array.from({ length: n }, () => new Array(n).fill(0));
  const eff    = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      let concs;
      if      (libOpt === 1) concs = [libRange[j], ratio[i]*libRange[j], ratio[i]*libRange[j]];
      else if (libOpt === 2) concs = [ratio[i]*libRange[j], libRange[j], ratio[i]*libRange[j]];
      else                   concs = [ratio[i]*libRange[j], ratio[i]*libRange[j], libRange[j]];
      const res = solveLinear3(K12, K23, concs[0], concs[1], concs[2], libOpt);
      if (res) { incorp[j][i] = res.incorp; eff[j][i] = res.eff; }
    }
  }
  return { ratio, libRange, incorp, eff };
}

// Circular recommendation  (Newton's)
function computeCirRecom(Ki, KcirVI, KcirVIVI, libConc) {
  const n = 30;
  const ratio    = logspace(-3, 3, n, 2);
  const libRange = linspace(libConc / 10, libConc, n);
  const incorp = Array.from({ length: n }, () => new Array(n).fill(0));
  const eff    = Array.from({ length: n }, () => new Array(n).fill(0));
  const accur  = Array.from({ length: n }, () => new Array(n).fill(0));
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < n; i++) {
      const concV = ratio[i] * libRange[j];
      const concI = libRange[j];
      const res = solveCircular(Ki, KcirVI, KcirVIVI, concV, concI);
      if (res) { incorp[j][i] = res.incorp; eff[j][i] = res.eff; accur[j][i] = res.accur; }
    }
  }
  return { ratio, libRange, incorp, eff, accur };
}

// ─────────────────────────────────────────────
// SE evaluation (mismatch calculations)
// ─────────────────────────────────────────────

function delGMmCal(seUp, seLow, seLeftStr, seRightStr, seType, stackOpt, delGLig) {
  const { bUpNum: seLeftUp, bLowNum: seLeftLow } = basesToNum(seLeftStr);
  const { bUpNum: seRightUp, bLowNum: seRightLow } = basesToNum(seRightStr);
  let uu, ll;
  const n = seUp.length - 1;
  if (seType === '5') {
    uu = hBond(seUp, seUp) + stackEnergy(stackOpt, seRightLow[0], seUp[0], seRightLow[0], seUp[0]) + delGLig;
    ll = hBond(seLow, seLow) + stackEnergy(stackOpt, seLeftUp[0], seLow[0], seLeftUp[0], seLow[0]) + delGLig;
  } else {
    const nl = seLow.length - 1;
    uu = hBond(seUp, seUp) + stackEnergy(stackOpt, seUp[n], seLeftLow[0], seUp[n], seLeftLow[0]) + delGLig;
    ll = hBond(seLow, seLow) + stackEnergy(stackOpt, seLow[nl], seRightUp[0], seLow[nl], seRightUp[0]) + delGLig;
  }
  const uu_gap = hBond(seUp.slice(0, -1), seUp.slice(0, -1)) + delGLig;
  const ll_gap = hBond(seLow.slice(0, -1), seLow.slice(0, -1)) + delGLig;
  const ul_gap = hBond(seUp.slice(0, -1), seLow.slice(0, -1)) + delGLig;
  return [uu, ll, uu_gap, ll_gap, ul_gap];
}

function delG2SeMmCal(se1Up, se1Low, se2Up, se2Low,
                      se1LeftStr, se1RightStr, se2LeftStr, se2RightStr,
                      seType, stackOpt, delGLig) {
  const { bUpNum: se1Lu, bLowNum: se1Ll } = basesToNum(se1LeftStr);
  const { bUpNum: se1Ru, bLowNum: se1Rl } = basesToNum(se1RightStr);
  const { bUpNum: se2Lu, bLowNum: se2Ll } = basesToNum(se2LeftStr);
  const { bUpNum: se2Ru, bLowNum: se2Rl } = basesToNum(se2RightStr);
  let u1u2, l1l2, u1l2, l1u2;
  const n1 = se1Up.length - 1, n2 = se2Up.length - 1;
  if (seType === '5') {
    u1u2 = hBond(se1Up, se2Up) + stackEnergy(stackOpt, se2Rl[0], se1Up[0], se1Rl[0], se2Up[0]) + delGLig;
    l1l2 = hBond(se1Low, se2Low) + stackEnergy(stackOpt, se2Lu[0], se1Low[0], se1Lu[0], se2Low[0]) + delGLig;
    u1l2 = hBond(se1Up, se2Low) + stackEnergy(stackOpt, se2Lu[0], se1Up[0], se1Rl[0], se2Low[0]) + delGLig;
    l1u2 = hBond(se1Low, se2Up) + stackEnergy(stackOpt, se1Lu[0], se2Up[0], se2Rl[0], se1Low[0]) + delGLig;
  } else {
    u1u2 = hBond(se1Up, se2Up) + stackEnergy(stackOpt, se1Up[n1], se2Ll[0], se2Up[n2], se1Ll[0]) + delGLig;
    l1l2 = hBond(se1Low, se2Low) + stackEnergy(stackOpt, se1Low[n1], se2Ru[0], se2Low[n2], se1Ru[0]) + delGLig;
    u1l2 = hBond(se1Up, se2Low) + stackEnergy(stackOpt, se1Up[n1], se2Ru[0], se2Low[n2], se1Ll[0]) + delGLig;
    l1u2 = hBond(se1Low, se2Up) + stackEnergy(stackOpt, se2Up[n2], se1Ru[0], se1Low[n1], se2Ll[0]) + delGLig;
  }
  const u1u2_gap = hBond(se1Up.slice(0,-1), se2Up.slice(0,-1)) + delGLig;
  const l1l2_gap = hBond(se1Low.slice(0,-1), se2Low.slice(0,-1)) + delGLig;
  const u1l2_gap = hBond(se1Up.slice(0,-1), se2Low.slice(0,-1)) + delGLig;
  const l1u2_gap = hBond(se1Low.slice(0,-1), se2Up.slice(0,-1)) + delGLig;
  return [u1u2, l1l2, u1l2, l1u2, u1u2_gap, l1l2_gap, u1l2_gap, l1u2_gap];
}

// ─────────────────────────────────────────────
// SE Finder scoring
// ─────────────────────────────────────────────

function mean(arr) { return arr.reduce((s,v) => s+v, 0) / arr.length; }
function std(arr)  {
  const m = mean(arr);
  return Math.sqrt(arr.reduce((s,v) => s + (v-m)**2, 0) / arr.length);
}

function scoreCal(delG, delGMm) {
  const mm_mean = mean(delGMm);
  const mm_std  = std(delGMm);
  if (mm_std === 0) return 0;
  return ((-1)*delG - (-1)*mm_mean) / mm_std;
}

function scoreCal2Se(d1G, d1Mm, d2G, d2Mm, d12Mm) {
  const all1 = [...d1Mm, ...d12Mm];
  const all2 = [...d2Mm, ...d12Mm];
  const m1 = mean(all1), s1 = std(all1);
  const m2 = mean(all2), s2 = std(all2);
  const sc1 = s1 > 0 ? ((-1)*d1G - (-1)*m1) / s1 : 0;
  const sc2 = s2 > 0 ? ((-1)*d2G - (-1)*m2) / s2 : 0;
  return sc1 + sc2;
}

// Find all sticky ends in a DNA sequence and return sorted result
// Returns array of {start, end, seqStr, score}
function findStickyEnds(dnaStr, overhang, seType, stackOpt, delGLig) {
  const { bUpNum, bLowNum, flag } = basesToNum(dnaStr);
  if (!flag || bUpNum.length < overhang + 2) return null;
  const results = [];
  const len = bUpNum.length;
  for (let i = 0; i <= len - overhang - 2; i++) {
    const se_l   = dnaStr[i];
    const seUp   = bUpNum.slice(i + 1, i + overhang + 1);
    const seLowr = bLowNum.slice(len - overhang - i - 1, len - i - 1);
    const se_r   = dnaStr[i + overhang + 1];
    const delG   = delGCal(seUp, seLowr, se_l, se_r, seType, stackOpt, delGLig);
    const delGMm = delGMmCal(seUp, seLowr, se_l, se_r, seType, stackOpt, delGLig);
    const score  = scoreCal(delG, delGMm);
    results.push({ start: i+2, end: i+overhang+1, seqStr: numToBase(seUp),
                   seUp: seUp.slice(), seLow: seLowr.slice(),
                   seLeftStr: se_l, seRightStr: se_r, seType,
                   delG, delGMm, score });
  }
  // Stable sort by descending score; ties broken by ascending start (mimics numpy argsort)
  results.sort((a, b) => b.score - a.score || a.start - b.start);
  return results;
}

// ─────────────────────────────────────────────
// SE fitting (grid search)
// ─────────────────────────────────────────────

function fitDelG(icArr, omrArr, bp, gLow, gHigh) {
  // ic: array of [c1,c2], omr: array of [mr1,mr2,mr3] normalized
  const T = 25 + 273.15, R = 1.9858775;
  const gSpace = [];
  for (let g = gLow; g <= gHigh + 0.005; g += 0.01) gSpace.push(+g.toFixed(5));
  const bpArr = [bp[0], bp[1], bp[0]+bp[1]];
  const bpRatio = bpArr.map(v => v / (bp[0]+bp[1]));
  const residue = new Array(gSpace.length).fill(0);
  for (let j = 0; j < gSpace.length; j++) {
    const K = Math.exp(gSpace[j]*1e3/(R*T))*1e9;
    let sumSq = 0;
    for (let i = 0; i < icArr.length; i++) {
      const [c1,c2] = icArr[i];
      const s = c1 + c2 + K;
      const AB = (s - Math.sqrt(s*s - 4*c1*c2)) / 2;
      if (!isFinite(AB) || AB < 0) continue;
      const A = c1 - AB, B = c2 - AB;
      const molArr = [A/AB, B/AB, 1.0];
      const massArr = molArr.map((m,k) => m * bpRatio[k]);
      for (let k = 0; k < 3; k++) sumSq += (massArr[k] - omrArr[i][k])**2;
    }
    residue[j] = sumSq;
  }
  let bestIdx = 0;
  for (let j = 1; j < residue.length; j++) if (residue[j] < residue[bestIdx]) bestIdx = j;
  // Window ±100 points around best fit (matches original Python behaviour)
  const wA = Math.max(0, bestIdx - 100);
  const wB = Math.min(gSpace.length, bestIdx + 100);
  return { bestDelG: gSpace[bestIdx],
           gSpace: gSpace.slice(wA, wB),
           residue: residue.slice(wA, wB) };
}

function defaultGRange(overhang) {
  const ranges = { 2: [-12.81,-6.88], 3: [-14.65,-7.76], 4: [-17.22,-8.34], 5: [-19.06,-9.22] };
  return ranges[overhang] || [-17.22,-8.34];
}

// ─────────────────────────────────────────────
// Canvas heatmap renderer
// ─────────────────────────────────────────────

function heatColor(t) {
  // blue (0) → cyan → green (0.5) → yellow → red (1)
  t = Math.max(0, Math.min(1, t));
  let r, g, b;
  if (t < 0.25) {
    r = 0; g = Math.round(t * 4 * 255); b = 255;
  } else if (t < 0.5) {
    r = 0; g = 255; b = Math.round((0.5 - t) * 4 * 255);
  } else if (t < 0.75) {
    r = Math.round((t - 0.5) * 4 * 255); g = 255; b = 0;
  } else {
    r = 255; g = Math.round((1 - t) * 4 * 255); b = 0;
  }
  return `rgb(${r},${g},${b})`;
}

/**
 * Draw a heatmap on a canvas.
 * @param {HTMLCanvasElement} canvas
 * @param {string} title
 * @param {string} xLabel
 * @param {string} yLabel
 * @param {number[]} xVals  - x axis values (ratio, may use log2 scale)
 * @param {number[]} yVals  - y axis values (linear)
 * @param {number[][]} data - data[j][i], j=y index, i=x index
 * @param {boolean} logX    - whether x-axis uses log2 scale
 * @param {string} unit     - '%' etc.
 */
function drawHeatmap(canvas, title, xLabel, yLabel, xVals, yVals, data, logX, unit) {
  const ctx  = canvas.getContext('2d');
  const W    = canvas.width;
  const H    = canvas.height;
  const ml = 72, mr = 70, mt = 36, mb = 56;
  const pw = W - ml - mr;
  const ph = H - mt - mb;

  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#fff';
  ctx.fillRect(0, 0, W, H);

  const nx = xVals.length, ny = yVals.length;

  // -- Find data range --
  let dMin = Infinity, dMax = -Infinity;
  for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
    if (data[j][i] < dMin) dMin = data[j][i];
    if (data[j][i] > dMax) dMax = data[j][i];
  }
  const dRange = dMax - dMin || 1;

  const mapX = logX
    ? (r) => ml + (Math.log2(r) - Math.log2(xVals[0])) / (Math.log2(xVals[nx-1]) - Math.log2(xVals[0])) * pw
    : (r) => ml + (r - xVals[0]) / (xVals[nx-1] - xVals[0]) * pw;
  const mapY = (v) => mt + ph - (v - yVals[0]) / (yVals[ny-1] - yVals[0]) * ph;

  // -- Draw cells --
  for (let j = 0; j < ny - 1; j++) {
    for (let i = 0; i < nx - 1; i++) {
      const t = (data[j][i] - dMin) / dRange;
      ctx.fillStyle = heatColor(t);
      const x0 = mapX(xVals[i]);
      const x1 = mapX(xVals[i + 1]);
      const y0 = mapY(yVals[j]);
      const y1 = mapY(yVals[j + 1]);
      ctx.fillRect(x0, y1, x1 - x0 + 1, y0 - y1 + 1);
    }
  }

  // -- Border --
  ctx.strokeStyle = '#333';
  ctx.lineWidth = 1;
  ctx.strokeRect(ml, mt, pw, ph);

  // -- X axis ticks --
  ctx.fillStyle = '#000';
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'center';
  const xTicks = logX
    ? [[1/8,'1/8'],[1/4,'1/4'],[1/2,'1/2'],[1,'1'],[2,'2'],[4,'4'],[8,'8']]
    : xVals.filter((_, i) => i % Math.floor(nx/6) === 0).map(v => [v, v.toFixed(2)]);
  for (const [xv, xl] of xTicks) {
    const px = mapX(xv);
    ctx.beginPath(); ctx.moveTo(px, mt+ph); ctx.lineTo(px, mt+ph+5); ctx.stroke();
    ctx.fillText(xl, px, mt + ph + 17);
  }

  // -- Y axis ticks --
  ctx.textAlign = 'right';
  const yTickCount = 5;
  const yStep = (yVals[ny-1] - yVals[0]) / yTickCount;
  for (let k = 0; k <= yTickCount; k++) {
    const yv = yVals[0] + k * yStep;
    const py = mapY(yv);
    ctx.beginPath(); ctx.moveTo(ml, py); ctx.lineTo(ml-5, py); ctx.stroke();
    ctx.fillText(yv.toFixed(1), ml - 8, py + 4);
  }

  // -- Axis labels --
  ctx.textAlign = 'center';
  ctx.font = '12px sans-serif';
  ctx.fillText(xLabel, ml + pw/2, H - 4);
  ctx.save();
  ctx.translate(14, mt + ph/2);
  ctx.rotate(-Math.PI/2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();

  // -- Title --
  ctx.font = 'bold 13px sans-serif';
  ctx.fillText(title, ml + pw/2, mt - 10);

  // -- Colorbar --
  const cbX = W - mr + 10, cbW = 18, cbH = ph;
  for (let k = 0; k < cbH; k++) {
    const t = 1 - k / cbH;
    ctx.fillStyle = heatColor(t);
    ctx.fillRect(cbX, mt + k, cbW, 1);
  }
  ctx.strokeStyle = '#333';
  ctx.strokeRect(cbX, mt, cbW, cbH);
  ctx.fillStyle = '#000';
  ctx.textAlign = 'left';
  ctx.font = '10px sans-serif';
  for (let k = 0; k <= 4; k++) {
    const val = dMin + (dMax - dMin) * (1 - k/4);
    const py  = mt + (k/4) * cbH;
    ctx.fillText(val.toFixed(1) + (unit||''), cbX + cbW + 3, py + 4);
  }
}

// -- Line chart for fitting results --
function drawLineChart(canvas, title, xLabel, yLabel, xVals, yVals, bestX) {
  const ctx = canvas.getContext('2d');
  const W = canvas.width, H = canvas.height;
  const ml = 70, mr = 20, mt = 36, mb = 56;
  const pw = W - ml - mr, ph = H - mt - mb;
  ctx.clearRect(0, 0, W, H);
  ctx.fillStyle = '#fff'; ctx.fillRect(0,0,W,H);

  const xMin = xVals[0], xMax = xVals[xVals.length-1];
  let yMin = Infinity, yMax = -Infinity;
  for (const v of yVals) { if (v < yMin) yMin = v; if (v > yMax) yMax = v; }
  const yRange = yMax - yMin || 1;
  const mapX = x => ml + (x - xMin)/(xMax - xMin)*pw;
  const mapY = y => mt + ph - (y - yMin)/yRange*ph;

  // Draw line
  ctx.beginPath(); ctx.strokeStyle = '#2255cc'; ctx.lineWidth = 2;
  for (let i = 0; i < xVals.length; i++) {
    if (i===0) ctx.moveTo(mapX(xVals[i]), mapY(yVals[i]));
    else       ctx.lineTo(mapX(xVals[i]), mapY(yVals[i]));
  }
  ctx.stroke();

  // Best point
  if (bestX !== undefined) {
    const idx = xVals.findIndex(v => Math.abs(v - bestX) < 0.015);
    if (idx >= 0) {
      ctx.beginPath(); ctx.fillStyle = 'red';
      ctx.arc(mapX(xVals[idx]), mapY(yVals[idx]), 5, 0, 2*Math.PI); ctx.fill();
    }
  }

  // Axes
  ctx.strokeStyle = '#333'; ctx.lineWidth = 1;
  ctx.strokeRect(ml, mt, pw, ph);
  ctx.fillStyle = '#000'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
  const xStep = (xMax - xMin) / 6;
  for (let k = 0; k <= 6; k++) {
    const xv = xMin + k*xStep, px = mapX(xv);
    ctx.beginPath(); ctx.moveTo(px, mt+ph); ctx.lineTo(px, mt+ph+5); ctx.stroke();
    ctx.fillText(xv.toFixed(1), px, mt+ph+17);
  }
  ctx.textAlign = 'right';
  const yStep = yRange / 4;
  for (let k = 0; k <= 4; k++) {
    const yv = yMin + k*yStep, py = mapY(yv);
    ctx.beginPath(); ctx.moveTo(ml, py); ctx.lineTo(ml-5, py); ctx.stroke();
    ctx.fillText(yv.toFixed(2), ml-8, py+4);
  }

  ctx.textAlign = 'center'; ctx.font = '12px sans-serif';
  ctx.fillText(xLabel, ml + pw/2, H - 4);
  ctx.save(); ctx.translate(14, mt+ph/2); ctx.rotate(-Math.PI/2);
  ctx.fillText(yLabel, 0, 0); ctx.restore();
  ctx.font = 'bold 13px sans-serif'; ctx.fillText(title, ml + pw/2, mt - 10);
}
