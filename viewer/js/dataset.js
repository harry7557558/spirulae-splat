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

// Camera-space border points of a size-1 frustum (CV convention).
// key: model + resolution + intrinsics + distortion.
function frustumTemplate(cam) {
  const { fx, fy, model, dist } = cam;
  const w = cam.w || Math.round(2*cam.cx) || 1;
  const h = cam.h || Math.round(2*cam.cy) || 1;
  const cx = cam.cx, cy = cam.cy;
  const corners = [
    [-cx/fx,      -cy/fy],
    [(w-cx)/fx,   -cy/fy],
    [(w-cx)/fx,   (h-cy)/fy],
    [-cx/fx,      (h-cy)/fy],
  ];
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
  const pts = [];
  for (let i = 0; i < 4*NSEG; i++) {
    const ci = Math.floor(i/NSEG), t = (i % NSEG)/NSEG;
    const u = corners[ci][0] + (corners[(ci+1)%4][0] - corners[ci][0])*t;
    const v = corners[ci][1] + (corners[(ci+1)%4][1] - corners[ci][1])*t;
    let dir = generateRay(u, v, model, dist);
    if (!dir) {
      // shrink uv toward the image center until it re-enters the valid domain
      let t0 = 0, t1 = 1, best = null;
      for (let k = 0; k < 12; k++) {
        const s = 0.5*(t0+t1);
        const rr = generateRay(u*s, v*s, model, dist);
        if (rr) { t0 = s; best = rr; } else t1 = s;
      }
      dir = best || [0, 0, 1];
    }
    pts.push(place(dir));
  }
  return pts;   // 4*NSEG points; corners at indices 0, NSEG, 2*NSEG, 3*NSEG
}

// Build the frustum line geometry for all cameras.
// Returns { apex:Float32Array, offset:Float32Array, colors:Uint8Array, count, ranges }.
export function buildFrustums(cameras) {
  const N = cameras.length;
  if (!N) return { apex: new Float32Array(0), offset: new Float32Array(0), colors: new Uint8Array(0), count: 0, ranges: [] };
  const ASEG = 4;                     // subdivisions per apex->corner anchor line
  const VPC = 4*NSEG*2 + 4*ASEG*2;    // vertices per camera (border loop + 4 anchors)
  const total = N*VPC;
  const apex = new Float32Array(total*3);
  const offset = new Float32Array(total*3);
  const colors = new Uint8Array(total*3);
  const ranges = new Array(N);

  const templates = new Map();
  const keyOf = (c) => c.model+'|'+c.w+'|'+c.h+'|'+c.fx.toFixed(3)+'|'+c.fy.toFixed(3)+
                       '|'+c.cx.toFixed(3)+'|'+c.cy.toFixed(3)+'|'+c.dist.join(',');

  let vp = 0;   // vertex pointer
  for (let i = 0; i < N; i++) {
    const cam = cameras[i];
    const key = keyOf(cam);
    let tmpl = templates.get(key);
    if (!tmpl) { tmpl = frustumTemplate(cam); templates.set(key, tmpl); }

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
    const worldOff = tmpl.map(off);

    const start = vp;
    const push = (o) => {
      apex[vp*3] = tx; apex[vp*3+1] = ty; apex[vp*3+2] = tz;
      offset[vp*3] = o[0]; offset[vp*3+1] = o[1]; offset[vp*3+2] = o[2];
      colors[vp*3] = FRUSTUM_COLOR[0]; colors[vp*3+1] = FRUSTUM_COLOR[1]; colors[vp*3+2] = FRUSTUM_COLOR[2];
      vp++;
    };
    // border loop
    const n = 4*NSEG;
    for (let j = 0; j < n; j++) { push(worldOff[j]); push(worldOff[(j+1)%n]); }
    // 4 anchor lines: apex -> each image corner, subdivided so they curve
    // correctly (and split at seams) under nonlinear display cameras
    for (let k = 0; k < 4; k++) {
      const o = worldOff[k*NSEG];
      for (let j = 0; j < ASEG; j++) {
        push([o[0]*j/ASEG, o[1]*j/ASEG, o[2]*j/ASEG]);
        push([o[0]*(j+1)/ASEG, o[1]*(j+1)/ASEG, o[2]*(j+1)/ASEG]);
      }
    }
    ranges[i] = { start, count: vp - start };
  }
  return { apex, offset, colors, count: total, ranges };
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
