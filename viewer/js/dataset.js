// Dataset loading + camera-frustum geometry for the standalone viewer.
//
// The heavy lifting (parsing COLMAP / Nerfstudio / Metashape) runs in the WASM
// bridge (dataset_bridge.cpp) over MEMFS; this module orchestrates the mount ->
// enumerate -> parse -> readback flow and builds the camera-frustum wireframe on
// the CPU from the per-camera intrinsics + distortion, porting
// fill_frustum_segments_kernel / projection_utils.cuh (generate_ray).
//
// Frustum vertices are emitted as (apex, unit-offset) pairs so the frustum-size
// slider is a cheap GPU uniform (see LINE_VS): apex = camera center (native
// frame), offset = the size-1 border/anchor point rotated into the native frame.

import {
  dsMount, dsEnumerate, dsParse, dsReadCameras, dsReadPoints,
  dsSummary, dsFitSphere, dsLastError,
} from './wasm.js';

// Engine CameraModelType ints (ParsedDataset.camera_models / Common.cuh:130).
const M_PINHOLE = 0, M_FISHEYE = 1, M_EQUISOLID = 2, M_EQUIRECT = 3;
const NSEG = 16;                       // segments per image edge (Visualizer.cu:93)
const FRUSTUM_COLOR = [110, 175, 255]; // base wireframe color (light blue)

// ---- distortion (OpenCV + thin-prism + b1/b2 mix), dist = [k1 k2 k3 k4 p1 p2 s1 s2 b1 b2]
function hasDistortion(d) { for (let i = 0; i < 10; i++) if (d[i]) return true; return false; }

function distort(u, v, d) {
  const r2 = u*u + v*v;
  const radial = 1 + r2*(d[0] + r2*(d[1] + r2*(d[2] + r2*d[3])));
  const xd = u*radial + 2*d[4]*u*v + d[5]*(r2 + 2*u*u) + d[6]*r2;
  const yd = v*radial + 2*d[5]*u*v + d[4]*(r2 + 2*v*v) + d[7]*r2;
  return [xd + d[8]*xd + d[9]*yd, yd];   // b1/b2 mix into x (projection_utils.cuh:1166)
}

// Solve distort(q) = (u,v) for q (undistorted), Newton with numerical
// Jacobian. Mirrors undistort_point_0 (projection_utils.cuh:1170-1263)
// INCLUDING its failure conditions: after the iterations the solution must
// have a well-posed forward Jacobian (min(det, J00, J11) > 0 — the
// is_valid_distortion criterion) and re-distort to within 0.01 of the input.
// Outside the calibrated domain (e.g. image corners of a >180° fisheye whose
// polynomial has gone non-monotonic) this returns null and the caller bisects
// toward the image center, exactly like fill_frustum_segments_kernel.
function jacobian(qx, qy, d) {
  const e = 1e-5;
  const f = distort(qx, qy, d);
  const fx = distort(qx + e, qy, d), fy = distort(qx, qy + e, d);
  return [ (fx[0]-f[0])/e, (fx[1]-f[1])/e,     // j00 j10
           (fy[0]-f[0])/e, (fy[1]-f[1])/e ];   // j01 j11
}
function undistort(u, v, d) {
  let qx = u, qy = v;
  for (let it = 0; it < 8; it++) {
    const f = distort(qx, qy, d);
    const rx = f[0] - u, ry = f[1] - v;
    const [j00, j10, j01, j11] = jacobian(qx, qy, d);
    const det = j00*j11 - j01*j10;
    if (Math.abs(det) < 1e-12) break;
    const inv = 1/det;
    qx -= ( j11*rx - j01*ry)*inv;
    qy -= (-j10*rx + j00*ry)*inv;
  }
  if (!isFinite(qx) || !isFinite(qy)) return null;
  const [J00, J10, J01, J11] = jacobian(qx, qy, d);
  const det = J00*J11 - J01*J10;
  if (Math.min(det, J00, J11) <= 0) return null;         // ill-posed / folded
  const f = distort(qx, qy, d);
  if (Math.hypot(f[0]-u, f[1]-v) >= 0.01) return null;   // did not converge
  return [qx, qy];
}

// Unproject a normalized image point uv=((px-cx)/fx,(py-cy)/fy) to a unit ray in
// camera space (CV convention: +Z forward, +Y down). Returns null when outside
// the valid domain. Mirrors generate_ray (projection_utils.cuh:1368).
function generateRay(u, v, model, d) {
  if (model === M_EQUIRECT) {
    if (Math.abs(u) > Math.PI || Math.abs(v) > Math.PI/2) return null;
    const cl = Math.cos(v);
    return [cl*Math.sin(u), Math.sin(v), cl*Math.cos(u)];
  }
  let uu = u, vv = v;
  if (hasDistortion(d)) {
    const q = undistort(u, v, d);
    if (!q) return null;               // outside the distortion's valid domain
    uu = q[0]; vv = q[1];
  }
  const r = Math.hypot(uu, vv);
  if (model === M_FISHEYE) {                       // equidistant: theta = r
    if (r >= Math.PI) return null;
    const s = r < 1e-3 ? (1 - r*r/6) : (Math.sin(r)/r);
    return norm3(uu*s, vv*s, Math.cos(r));
  }
  if (model === M_EQUISOLID) {                     // r = 2 sin(theta/2)
    if (r >= 2.0) return null;
    const k = Math.sqrt(Math.max(0, 1 - 0.25*r*r));
    return norm3(uu*k, vv*k, 1 - 0.5*r*r);
  }
  return norm3(uu, vv, 1.0);                        // pinhole
}

function norm3(x, y, z) { const l = Math.hypot(x, y, z) || 1; return [x/l, y/l, z/l]; }

// Camera-space size-1 frustum wireframe (CV convention).
// key: model + resolution + intrinsics + distortion.
//
// Pinhole keeps the classic look: image-border rectangle + 4 apex->corner
// anchors. Wide models (fisheye / equisolid / equirectangular) sit on a sphere
// shell where a lone border loop reads as a scribble (a >180-degree fisheye's
// border is a ring *behind* the camera), so they additionally get image-aligned
// interior gridlines — the wireframe of the textured uv-grid surface the
// training viewer rasterizes (fill_frustum_segments_kernel faces) — making the
// dome / globe shape and the view orientation legible. For equirectangular the
// corner-interpolated border is degenerate (poles collapse, left/right edges
// coincide), so the border is dropped entirely in favor of the pixel-grid
// meridians + parallels.
//
// Returns { lines: [{pts, closed, dim}], anchors: [point], maxR }.
function frustumTemplate(cam) {
  const { fx, fy, model, dist } = cam;
  const w = cam.w || Math.round(2*cam.cx) || 1;
  const h = cam.h || Math.round(2*cam.cy) || 1;
  const cx = cam.cx, cy = cam.cy;
  // Depth-placement scale factors (Visualizer.cu:123-127) at size = 1.
  const a = Math.sqrt(fx*fy/(w*h));
  const rs = 1/Math.sqrt(a);
  const sxy = a*rs;              // = sqrt(a)
  const sz  = 1/rs;             // = sqrt(a)  (equal; kept general)
  const wide = (model === M_FISHEYE || model === M_EQUISOLID || model === M_EQUIRECT);
  const shell = Math.sqrt((2*sxy*sxy + sz*sz)/3);
  const place = (dir) => {
    if (wide) return [dir[0]*shell, dir[1]*shell, dir[2]*shell];   // dir is unit
    if (Math.abs(dir[2]) < 1e-6) return [0,0,0];
    return [dir[0]*sxy/dir[2], dir[1]*sxy/dir[2], sz];             // z = sz plane
  };
  // Unproject one normalized image point; outside the valid domain, shrink uv
  // toward the principal point until it re-enters (Visualizer.cu:146-157).
  const ray = (u, v) => {
    let dir = generateRay(u, v, model, dist);
    if (!dir) {
      let t0 = 0, t1 = 1, best = null;
      for (let k = 0; k < 12; k++) {
        const s = 0.5*(t0+t1);
        const rr = generateRay(u*s, v*s, model, dist);
        if (rr) { t0 = s; best = rr; } else t1 = s;
      }
      dir = best || [0, 0, 1];
    }
    return dir;
  };
  // Polyline between two pixel coords, n segments, placed on the frustum.
  // Wide-model gridlines span up to a full circle, so they get 2x the border's
  // per-edge sampling to stay smooth close up.
  const sample = (x0, y0, x1, y1, n = NSEG) => {
    const pts = [];
    for (let i = 0; i <= n; i++) {
      const t = i/n;
      const u = (x0 + (x1-x0)*t - cx)/fx, v = (y0 + (y1-y0)*t - cy)/fy;
      pts.push(place(ray(u, v)));
    }
    return pts;
  };
  const degenerate = (pts) => {  // e.g. an equirect border row collapsed to a pole
    const p0 = pts[0];
    return pts.every(p => Math.hypot(p[0]-p0[0], p[1]-p0[1], p[2]-p0[2]) < 1e-5*shell + 1e-12);
  };

  const lines = [], anchors = [];
  if (model === M_EQUIRECT) {
    // Lat/long wire globe over the pixel grid (covers partial panoramas too:
    // the i=0/4, j=0/4 lines *are* the image border). Center meridian + center
    // parallel bright as the view-direction cue, the rest dim.
    for (let i = 0; i <= 4; i++) {
      const mer = sample(w*i/4, 0, w*i/4, h, 2*NSEG);
      if (!degenerate(mer)) lines.push({ pts: mer, closed: false, dim: i !== 2 });
      const par = sample(0, h*i/4, w, h*i/4, 2*NSEG);
      if (!degenerate(par)) lines.push({ pts: par, closed: false, dim: i !== 2 });
    }
    anchors.push(place(ray((w/2-cx)/fx, (h/2-cy)/fy)));      // view direction
  } else {
    // image-border loop (4 edges, corners at indices 0, NSEG, 2*NSEG, 3*NSEG)
    const corners = [[0,0], [w,0], [w,h], [0,h]];
    const border = [];
    for (let e = 0; e < 4; e++) {
      const [xa, ya] = corners[e], [xb, yb] = corners[(e+1)%4];
      border.push(...sample(xa, ya, xb, yb).slice(0, NSEG));
    }
    lines.push({ pts: border, closed: true, dim: false });
    if (wide) {
      // interior gridlines -> wire dome
      for (let i = 1; i <= 3; i++) {
        lines.push({ pts: sample(w*i/4, 0, w*i/4, h, 2*NSEG), closed: false, dim: true });
        lines.push({ pts: sample(0, h*i/4, w, h*i/4, 2*NSEG), closed: false, dim: true });
      }
    }
    // Anchors: apex->corner when the corners are in the front hemisphere
    // (classic pyramid); for wider cameras they cross the dome interior and
    // scribble, so use a single apex->view-direction anchor instead.
    const cornerPts = [border[0], border[NSEG], border[2*NSEG], border[3*NSEG]];
    if (!wide || cornerPts.every(p => p[2] > 0))
      anchors.push(...cornerPts);
    else
      anchors.push(place(ray((w/2-cx)/fx, (h/2-cy)/fy)));
  }

  let maxR = 0;
  for (const l of lines) for (const p of l.pts) maxR = Math.max(maxR, Math.hypot(p[0], p[1], p[2]));
  for (const p of anchors) maxR = Math.max(maxR, Math.hypot(p[0], p[1], p[2]));
  return { lines, anchors, maxR };
}

// GL line-list vertex count of a template (used to preallocate buffers).
function templateVertCount(tmpl) {
  const ASEG = 4;
  let n = 0;
  for (const l of tmpl.lines) n += 2*(l.pts.length - 1 + (l.closed ? 1 : 0));
  return n + tmpl.anchors.length*ASEG*2;
}

// Build the frustum line geometry for all cameras.
// Returns { apex:Float32Array, offset:Float32Array, colors:Uint8Array, count,
//           ranges, pickR:Float32Array } — pickR[i] is the size-1 bounding-
//           sphere radius around camera i's center (multiply by the frustum-
//           size slider for ray picking, see pickCamera in main.js).
export function buildFrustums(cameras) {
  const N = cameras.length;
  if (!N) return { apex: new Float32Array(0), offset: new Float32Array(0), colors: new Uint8Array(0), count: 0, ranges: [], pickR: new Float32Array(0) };
  const ASEG = 4;                     // subdivisions per apex anchor line

  const templates = new Map();
  const keyOf = (c) => c.model+'|'+c.w+'|'+c.h+'|'+c.fx.toFixed(3)+'|'+c.fy.toFixed(3)+
                       '|'+c.cx.toFixed(3)+'|'+c.cy.toFixed(3)+'|'+c.dist.join(',');
  const camTmpl = new Array(N);
  let total = 0;
  for (let i = 0; i < N; i++) {
    const key = keyOf(cameras[i]);
    let tmpl = templates.get(key);
    if (!tmpl) { tmpl = frustumTemplate(cameras[i]); templates.set(key, tmpl); }
    camTmpl[i] = tmpl;
    total += templateVertCount(tmpl);
  }

  const apex = new Float32Array(total*3);
  const offset = new Float32Array(total*3);
  const colors = new Uint8Array(total*3);
  const ranges = new Array(N);
  const pickR = new Float32Array(N);
  const DIM = FRUSTUM_COLOR.map(c => Math.round(c*0.5));   // interior gridlines

  let vp = 0;   // vertex pointer
  for (let i = 0; i < N; i++) {
    const cam = cameras[i];
    const tmpl = camTmpl[i];
    pickR[i] = tmpl.maxR;

    // c2w [3,4] row-major, OpenGL/nerfstudio. CV cam->world rotation = negate
    // the y,z columns; translation is the camera center.
    const m = cam.c2w;
    const R = [
      [m[0], -m[1], -m[2]],
      [m[4], -m[5], -m[6]],
      [m[8], -m[9], -m[10]],
    ];
    const tx = m[3], ty = m[7], tz = m[11];
    // rotate a CV cam-space point into the native-frame offset (apex-relative)
    const off = (p) => [
      R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2],
      R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2],
      R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2],
    ];

    const start = vp;
    const push = (o, col) => {
      apex[vp*3] = tx; apex[vp*3+1] = ty; apex[vp*3+2] = tz;
      offset[vp*3] = o[0]; offset[vp*3+1] = o[1]; offset[vp*3+2] = o[2];
      colors[vp*3] = col[0]; colors[vp*3+1] = col[1]; colors[vp*3+2] = col[2];
      vp++;
    };
    for (const line of tmpl.lines) {
      const pts = line.pts.map(off);
      const col = line.dim ? DIM : FRUSTUM_COLOR;
      const n = pts.length;
      for (let j = 0; j + 1 < n; j++) { push(pts[j], col); push(pts[j+1], col); }
      if (line.closed) { push(pts[n-1], col); push(pts[0], col); }
    }
    // anchor lines: apex -> corner / view direction, subdivided so they curve
    // correctly (and split at seams) under nonlinear display cameras
    for (const p of tmpl.anchors) {
      const o = off(p);
      for (let j = 0; j < ASEG; j++) {
        push([o[0]*j/ASEG, o[1]*j/ASEG, o[2]*j/ASEG], FRUSTUM_COLOR);
        push([o[0]*(j+1)/ASEG, o[1]*(j+1)/ASEG, o[2]*(j+1)/ASEG], FRUSTUM_COLOR);
      }
    }
    ranges[i] = { start, count: vp - start };
  }
  return { apex, offset, colors, count: total, ranges, pickR };
}

// ---------------------------------------------------------------------------
// Orchestration: mount files, enumerate components, parse one, read back.
// entries: [{file, path}]. token: optional component token (default first).
// Returns everything main.js needs; caller uploads points/frustum to the GPU
// then frees the wasm-side dataset.
// ---------------------------------------------------------------------------
export async function loadDatasetFiles(entries) {
  await dsMount(entries);
  const components = dsEnumerate();
  if (!components.length) throw new Error('no COLMAP / Nerfstudio / Metashape / point-cloud data found');
  const first = components[0].token;
  const res = parseDatasetComponent(first);
  return { components, selectedToken: first, ...res };
}

// Parse (or re-parse) a specific component; the files are already mounted.
export function parseDatasetComponent(token) {
  const ok = dsParse(token);
  const error = dsLastError();          // "" = clean; non-empty + ok = partial
  if (!ok) throw new Error(error || 'failed to parse dataset component');
  const cameras = dsReadCameras();      // copied
  const summary = dsSummary();          // copied
  const fit = dsFitSphere();            // copied
  const frustum = buildFrustums(cameras);
  const points = dsReadPoints();        // live heap views — upload before dsFree
  return { cameras, summary, fit, frustum, points, error };
}
