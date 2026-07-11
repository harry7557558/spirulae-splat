// Application entry: wires the WASM loader, WebGL renderer, camera controller,
// and the UI panel together.

import { initWasm, loadModel, sortSplats, freeSplatSh, splatHistogram, meshHistogram, meshEdgeCount, fitSphere, raycastMesh } from './wasm.js';
import { Renderer } from './renderer.js';
import { Camera, Nav } from './camera.js';
import { v3, mat3 } from './linalg.js';
import { GAMUTS, toColMajor, hexToRgb } from './colors.js';
import { drawHistogram } from './histogram.js';

const $ = (id) => document.getElementById(id);
const canvas = $('glcanvas');

let renderer, camera, nav;
let dirty = true;
let model = null;          // { type:'splat'|'mesh', ... meta }
let lastSortDir = [0,0,0];
const opts = {
  primitive: 0, gamut: toColMajor(GAMUTS['Rec.709']), isLinear: false,
  shDegree: 0, exposure: 1.0, opacityScale: 1.0,
  cameraModel: 'perspective', upAxis: 'z', showGrid: true,
  background: hexToRgb('0a0b0e'), shade: true, flatShade: false, meshColor: true,
};

function setStatus(text, cls='') { $('status-text').textContent = text; $('status-dot').className = cls; }

// ---------------------------------------------------------------------------
// Boot
// ---------------------------------------------------------------------------
async function boot() {
  setStatus('Loading engine…');
  try {
    await initWasm();
    renderer = new Renderer(canvas);
  } catch (e) {
    setStatus('Init failed: ' + e.message, 'err');
    console.error(e);
    return;
  }
  camera = new Camera();
  nav = new Nav(camera);
  nav.onChange = () => { dirty = true; };
  wireInput();
  wireControls();
  syncControls();       // apply any browser-restored control values on refresh
  resize();
  new ResizeObserver(resize).observe($('viewport'));
  setStatus('Ready — drop a model', 'ok');
  requestAnimationFrame(loop);

  // Optional auto-load: index.html?model=<url> (handy for hosting a default
  // model or for automated testing).
  const url = new URLSearchParams(location.search).get('model');
  if (url) {
    try { await load(await fetchModelFiles(url)); }
    catch (e) { setStatus('Auto-load failed: ' + e.message, 'err'); }
  }
}

// Fetch a model URL plus any sibling files it references (so hosting a
// textured .gltf/.obj with external .bin/.mtl/image works from a link).
async function fetchModelFiles(url) {
  const dir = url.slice(0, url.lastIndexOf('/')+1);
  const name = url.split('/').pop().split('?')[0];
  const fetchFile = async (u, n) => new File([await (await fetch(u)).blob()], n);
  const files = [await fetchFile(url, name)];
  const lname = name.toLowerCase();
  const add = async (n) => { try { files.push(await fetchFile(dir+n, n)); } catch {} };
  if (lname.endsWith('.gltf')) {
    const g = JSON.parse(await files[0].text());
    for (const b of g.buffers||[]) if (b.uri && !b.uri.startsWith('data:')) await add(b.uri.split('/').pop());
    for (const im of g.images||[]) if (im.uri && !im.uri.startsWith('data:')) await add(im.uri.split('/').pop());
  } else if (lname.endsWith('.obj')) {
    const t = await files[0].text();
    const mtl = t.match(/^mtllib\s+(.+)$/m);
    if (mtl) { const mn = mtl[1].trim(); await add(mn);
      const mf = files.find(f=>f.name===mn);
      if (mf) { const mm = (await mf.text()).match(/^\s*map_Kd\s+(.+)$/m); if (mm) await add(mm[1].trim().split(/[\\/]/).pop()); } }
  }
  return files;
}

// expose for automated testing / embedding
window.__viewer = {
  load: (files) => load(files),
  get model() { return model; },
  get opts() { return opts; },
  setOpt: (k, v) => { opts[k] = v; dirty = true; },
  get nav() { return nav; },
  get camera() { return camera; },
  snapshot: () => { if (model && model.type==='splat') maybeSort(true); renderer.render(camera, opts); return renderer.snapshot(); },
};

// ---------------------------------------------------------------------------
// Render loop
// ---------------------------------------------------------------------------
let lastT = performance.now(), frames = 0, fpsT = performance.now();
function loop(now) {
  const dt = Math.min((now - lastT)/1000, 0.1); lastT = now;
  if (nav.tick(dt)) dirty = true;
  if (nav.gamepadTick(dt)) dirty = true;

  if (dirty) {
    if (model && model.type === 'splat') maybeSort(false);
    renderer.render(camera, opts);
    dirty = false;
  }
  // fps
  frames++;
  if (now - fpsT > 500) { $('fps').textContent = (frames*1000/(now-fpsT)|0) + ' fps'; frames=0; fpsT=now; }
  requestAnimationFrame(loop);
}

function upTransform() {
  return opts.upAxis === 'y' ? [1,0,0, 0,0,1, 0,-1,0] : [1,0,0, 0,1,0, 0,0,1];
}
// Fit the camera to the current model using a robust (median-based) center
// and radius — trained splats and generated meshes often contain far-away
// outliers that make a bounding-box fit useless.
function fitModel() {
  if (!model) return;
  const fs = fitSphere();                       // [cx,cy,cz, medianDist] (native frame)
  const c = mat3.mulVec(upTransform(), [fs[0], fs[1], fs[2]]);
  const r = 2.0 * fs[3] || 1;
  nav.fitSphere(c, r, [0,0,1]);
}
// Switch the up axis: the model (and its native-frame grid) visibly rotates
// into the new canonical orientation; the camera stays put, so there is no
// viewport reset and the toggle takes effect immediately.
function setUpAxis(axis) {
  if (axis === opts.upAxis) return;
  opts.upAxis = axis;
  lastSortDir = [0,0,0];
  dirty = true;
}

function forwardNative() {
  // camera forward in the model's native frame: U^T * forward_canonical
  return v3.norm(mat3.mulVec(mat3.transpose(upTransform()), camera.forward()));
}
function camPosNative() {
  return mat3.mulVec(mat3.transpose(upTransform()), camera.pos);
}
// Sorting-depth mode (must match ssv_sort): planar depth for perspective /
// orthographic; camera distance for equirectangular; a smooth |z|/radial
// blend for fisheye so >180° views with splats behind the camera sort right.
function sortMode() {
  return camera.model === 'equirectangular' ? 1
       : (camera.model === 'fisheye' || camera.model === 'equisolid') ? 2 : 0;
}
let lastSortPos = [0,0,0], lastSortMode = -1;
function maybeSort(force) {
  if (!model || model.type !== 'splat') return;
  const d = forwardNative();
  const mode = sortMode();
  const c = camPosNative();
  if (!force && mode === lastSortMode &&
      v3.dot(d, lastSortDir) > Math.cos(0.06) &&
      // nonplanar depths depend on the camera position too
      !(mode > 0 && v3.len(v3.sub(c, lastSortPos)) > 0.01 * nav._sceneScale)) return;
  const order = sortSplats(mode, c[0], c[1], c[2], d[0], d[1], d[2]);
  if (order) renderer.updateSplatOrder(order);
  lastSortDir = d; lastSortPos = c; lastSortMode = mode;
}

// ---------------------------------------------------------------------------
// Resize (device-pixel aware, capped for performance)
// ---------------------------------------------------------------------------
function resize() {
  const vp = $('viewport');
  const dpr = window.devicePixelRatio || 1;
  let w = Math.round(vp.clientWidth * dpr), h = Math.round(vp.clientHeight * dpr);
  const maxPix = 2500000;
  const scale = Math.min(1, Math.sqrt(maxPix / Math.max(w*h, 1)));
  w = Math.max(1, Math.round(w*scale)); h = Math.max(1, Math.round(h*scale));
  if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; dirty = true; }
}

// ---------------------------------------------------------------------------
// Model loading
// ---------------------------------------------------------------------------
async function load(files) {
  setStatus('Loading model…', '');
  try {
    const res = await loadModel(files);
    const firstModel = !model;
    if (res.kind === 'splat') {
      renderer.setSplat(res.data);
      freeSplatSh();   // SH lives on the GPU now; drop the WASM copy
      model = { type:'splat', count: res.data.count, shDegree: res.data.shDegree };
      showSplatUI(res.data);
    } else {
      renderer.setMesh(res.data, res.image || null);
      model = { type:'mesh', nv: res.data.nv, nt: res.data.nt };
      showMeshUI(res.data);
    }
    // Only fit the camera for the first model of the session: replacing the
    // model (e.g. dropping the mesh of the same object after a splat) keeps
    // the current view. The scene scale (move/zoom speed) is still refreshed.
    if (firstModel) fitModel();
    else { const fs = fitSphere(); nav._sceneScale = 2.0 * fs[3] || 1; }
    lastSortDir = [0,0,0];
    histCache.clear();
    updateHistParams();
    $('drop-hint').style.display = 'none';
    const name = Array.from(files).find(f=>/\.(ply|obj|gltf|glb)$/i.test(f.name))?.name || '';
    $('st-file').textContent = name.length>22 ? '…'+name.slice(-21) : name;
    setStatus(model.type==='splat' ? 'Splat model loaded' : 'Mesh loaded', 'ok');
    dirty = true;
  } catch (e) {
    console.error(e);
    setStatus('Load error: ' + e.message, 'err');
  }
}

function showSplatUI(data) {
  $('sec-splat').classList.remove('hidden');
  $('sec-mesh').classList.add('hidden');
  $('st-type').textContent = '3D Gaussian Splats';
  $('st-splat-stats').style.display = '';
  $('st-mesh-stats').style.display = 'none';
  $('st-count').textContent = data.count.toLocaleString();
  $('st-sh').textContent = data.shDegree;
  // slider range 0..model degree, defaulting to the model's actual degree
  const sh = $('shdeg'); sh.max = data.shDegree; sh.value = data.shDegree;
  opts.shDegree = data.shDegree; $('v-shdeg').textContent = sh.value;
}
function showMeshUI(data) {
  $('sec-mesh').classList.remove('hidden');
  $('sec-splat').classList.add('hidden');
  $('st-type').textContent = data.uv ? 'Textured mesh' : (data.color ? 'Vertex-color mesh' : 'Mesh');
  $('st-splat-stats').style.display = 'none';
  $('st-mesh-stats').style.display = '';
  $('st-verts').textContent = data.nv.toLocaleString();
  $('st-faces').textContent = data.nt.toLocaleString();
  $('st-edges').textContent = '…';
  // edge count can be a touch slow for huge meshes; compute async-ish
  setTimeout(() => { try { $('st-edges').textContent = meshEdgeCount().toLocaleString(); } catch {} }, 0);
}

// ---------------------------------------------------------------------------
// Histograms
// ---------------------------------------------------------------------------
const SPLAT_PARAMS = [
  ['Opacity',0],['Scale (log₁₀)',1],['Effective rank',2],
  ['Color R',3],['Color G',4],['Color B',5],
  ['Scale X (log₁₀)',6],['Scale Y (log₁₀)',7],['Scale Z (log₁₀)',8],['Anisotropy (log₁₀)',9],
];
const MESH_PARAMS = [
  ['Edge length (log₁₀)',0],['Triangle area (log₁₀)',1],['X',2],['Y',3],['Z',4],
];
const histCache = new Map();
function updateHistParams() {
  const sel = $('hist-param'); sel.innerHTML = '';
  const params = model && model.type==='mesh' ? MESH_PARAMS : SPLAT_PARAMS;
  for (const [name, id] of params) { const o=document.createElement('option'); o.value=id; o.textContent=name; sel.appendChild(o); }
  const ctx = $('hist-canvas').getContext('2d'); ctx.clearRect(0,0,9999,9999);
}
function computeHistogram() {
  if (!model) return;
  const param = +$('hist-param').value;
  const key = model.type + ':' + param;
  let hist = histCache.get(key);
  if (!hist) {
    hist = model.type==='mesh' ? meshHistogram(param, 64) : splatHistogram(param, 64);
    histCache.set(key, hist);
  }
  const label = ($('hist-param').selectedOptions[0]||{}).textContent || '';
  drawHistogram($('hist-canvas'), hist, label);
}

// ---------------------------------------------------------------------------
// Input handling (pointer / wheel / touch / keyboard)
// ---------------------------------------------------------------------------
function wireInput() {
  const vp = $('viewport');
  const pointers = new Map();
  let prevPinch = null;

  vp.addEventListener('pointerdown', (e) => {
    vp.setPointerCapture(e.pointerId);
    pointers.set(e.pointerId, { x:e.clientX, y:e.clientY, button:e.button });
    if (pointers.size === 2) prevPinch = pinchState(pointers);
  });
  vp.addEventListener('pointermove', (e) => {
    const p = pointers.get(e.pointerId); if (!p) return;
    const dx = e.clientX - p.x, dy = e.clientY - p.y;
    p.x = e.clientX; p.y = e.clientY;
    if (pointers.size === 1) {
      const panning = e.shiftKey || p.button === 1 || p.button === 2;
      nav.drag(dx, dy, panning);
      dirty = true;
    } else if (pointers.size === 2) {
      const ps = pinchState(pointers);
      if (prevPinch) {
        nav.dolly((prevPinch.dist - ps.dist) * 4);
        nav.pan(ps.cx - prevPinch.cx, ps.cy - prevPinch.cy);
        let dang = ps.ang - prevPinch.ang;
        dang = (dang + Math.PI) % (2*Math.PI) - Math.PI;
        nav.roll(dang);
        dirty = true;
      }
      prevPinch = ps;
    }
  });
  const up = (e) => { pointers.delete(e.pointerId); if (pointers.size < 2) prevPinch = null; };
  vp.addEventListener('pointerup', up);
  vp.addEventListener('pointercancel', up);
  vp.addEventListener('wheel', (e) => { e.preventDefault(); nav.dolly(e.deltaY); dirty = true; }, { passive:false });
  vp.addEventListener('dblclick', (e) => { e.preventDefault(); pickRecenter(e); });
  vp.addEventListener('contextmenu', (e) => e.preventDefault());

  window.addEventListener('keydown', (e) => {
    if (e.target.tagName === 'SELECT' || e.target.tagName === 'INPUT') return;
    const k = e.key.toLowerCase();
    if (['w','a','s','d','q','e','arrowup','arrowdown','arrowleft','arrowright'].includes(k)) { e.preventDefault(); nav.keyDown(k); }
  });
  window.addEventListener('keyup', (e) => nav.keyUp(e.key.toLowerCase()));
  window.addEventListener('blur', () => nav._keys.clear());

  // drag & drop / file picker
  vp.addEventListener('dragover', (e) => { e.preventDefault(); vp.classList.add('dragover'); });
  vp.addEventListener('dragleave', () => vp.classList.remove('dragover'));
  vp.addEventListener('drop', (e) => { e.preventDefault(); vp.classList.remove('dragover'); if (e.dataTransfer.files.length) load(e.dataTransfer.files); });
  $('file-input').addEventListener('change', (e) => { if (e.target.files.length) load(e.target.files); });
  $('drop-hint').addEventListener('pointerup', () => $('file-input').click());
}
// ---------------------------------------------------------------------------
// Double-click view centering (MeshLab-style): find the 3D point under the
// cursor — GPU depth accumulation for splats, WASM brute-force raycast for
// meshes — then pan so it sits on the optical axis and make it the orbit
// pivot. Works with every camera model.
// ---------------------------------------------------------------------------
function canvasPixel(e) {
  const rect = canvas.getBoundingClientRect();
  return [ (e.clientX - rect.left) * canvas.width / rect.width,
           canvas.height - (e.clientY - rect.top) * canvas.height / rect.height ];
}
// camera-space ray direction for pixel-plane uv (matches unprojectCam in GLSL)
function unprojectDir(uv) {
  const m = camera.model;
  if (m === 'perspective' || m === 'orthographic') return v3.norm([uv[0], uv[1], -1]);
  if (m === 'equirectangular') {
    const cp = Math.cos(uv[1]);
    return [cp*Math.sin(uv[0]), Math.sin(uv[1]), -cp*Math.cos(uv[0])];
  }
  const r = Math.hypot(uv[0], uv[1]);
  const d = r > 1e-9 ? [uv[0]/r, uv[1]/r] : [0, 0];
  const th = m === 'equisolid' ? 2*Math.asin(Math.max(-1, Math.min(1, 0.5*r))) : r;
  return [d[0]*Math.sin(th), d[1]*Math.sin(th), -Math.cos(th)];
}
function pickRecenter(e) {
  if (!model) return;
  const [px, py] = canvasPixel(e);
  const intr = renderer._intrinsics(camera, canvas.width, canvas.height);
  const uv = [(px - intr.cx)/intr.fx, (py - intr.cy)/intr.fy];
  const R = camera.rotMat3();            // cam -> world (canonical)
  let p = null;
  if (model.type === 'splat') {
    const depth = renderer.pickSplatDepth(camera, opts, px, py);
    if (depth != null) {
      let pc;
      if (camera.model === 'orthographic') pc = [uv[0], uv[1], -depth];
      else if (camera.model === 'perspective') pc = [uv[0]*depth, uv[1]*depth, -depth]; // axial depth
      else pc = v3.scale(unprojectDir(uv), depth);                                     // ray length
      p = v3.add(camera.pos, mat3.mulVec(R, pc));
    }
  } else {
    let ro, rd;
    if (camera.model === 'orthographic') {
      ro = v3.add(camera.pos, mat3.mulVec(R, [uv[0], uv[1], 0]));
      rd = mat3.mulVec(R, [0, 0, -1]);
    } else {
      ro = camera.pos;
      rd = mat3.mulVec(R, unprojectDir(uv));
    }
    const Ut = mat3.transpose(upTransform());
    const on = mat3.mulVec(Ut, ro), dn = mat3.mulVec(Ut, rd);
    const t = raycastMesh(on[0], on[1], on[2], dn[0], dn[1], dn[2]);
    if (t > 0) p = v3.add(ro, v3.scale(rd, t));
  }
  if (p) recenterAt(p);
}
// Pan laterally so p lies on the optical axis (view orientation preserved) and
// set it as the orbit pivot; if p is behind the camera plane (possible with
// >180° cameras) rotate toward it instead.
function recenterAt(p) {
  const fwd = camera.forward();
  const rel = v3.sub(p, camera.pos);
  const d = v3.dot(rel, fwd);
  if (d > 0) camera.pos = v3.add(camera.pos, v3.sub(rel, v3.scale(fwd, d)));
  else camera.lookAt(camera.pos, p, camera.up());
  camera.target = p.slice();
  dirty = true;
}

function pinchState(pointers) {
  const [a,b] = [...pointers.values()];
  // angle convention matches the training viewer (atan2(dx, dy)) so two-finger
  // rotation rolls in the same direction
  return { dist: Math.hypot(a.x-b.x, a.y-b.y), cx:(a.x+b.x)/2, cy:(a.y+b.y)/2, ang: Math.atan2(a.x-b.x, a.y-b.y) };
}
// ---------------------------------------------------------------------------
// UI control wiring
// ---------------------------------------------------------------------------
function wireControls() {
  const bind = (id, ev, fn) => $(id).addEventListener(ev, fn);
  bind('nav-mode','change', e => nav.mode = e.target.value);
  bind('move-speed','input', e => { nav.speedExp = +e.target.value; $('v-speed').textContent = Math.pow(10,+e.target.value).toFixed(2)+'×'; });
  bind('btn-reset','click', () => { if (model) fitModel(); else nav.reset(); });

  bind('primitive','change', e => { opts.primitive = +e.target.value; dirty=true; });
  bind('gamut','change', e => { opts.gamut = toColMajor(GAMUTS[e.target.value]); dirty=true; });
  bind('linear','change', e => { opts.isLinear = e.target.checked; dirty=true; });
  bind('shdeg','input', e => { opts.shDegree = +e.target.value; $('v-shdeg').textContent = e.target.value; dirty=true; });
  bind('exposure','input', e => { opts.exposure = Math.pow(10,+e.target.value); $('v-exposure').textContent = opts.exposure.toFixed(2); dirty=true; });
  bind('shade','change', e => { opts.shade = e.target.checked; dirty=true; });

  bind('camera-model','change', e => {
    opts.cameraModel = e.target.value; camera.model = e.target.value;
    // FOV is meaningless for equirectangular (full 360×180)
    $('fov').disabled = (e.target.value === 'equirectangular');
    $('fov').parentElement.style.opacity = $('fov').disabled ? 0.4 : 1;
    dirty=true;
  });
  bind('fov','input', e => { camera.fov = +e.target.value*Math.PI/180; $('v-fov').textContent = e.target.value+'°'; dirty=true; });
  bind('up-axis','change', e => setUpAxis(e.target.checked ? 'y' : 'z'));
  bind('flat-shade','change', e => { opts.flatShade = e.target.checked; dirty=true; });
  bind('mesh-color','change', e => { opts.meshColor = e.target.checked; dirty=true; });
  bind('grid','change', e => { opts.showGrid = e.target.checked; dirty=true; });
  bind('bg-color','input', e => { opts.background = hexToRgb(e.target.value); dirty=true; });

  bind('hist-param','change', () => computeHistogram());
  bind('btn-hist','click', () => computeHistogram());
  bind('panel-toggle','click', () => $('panel').classList.toggle('collapsed'));

  updateNavLegend();
  $('nav-mode').addEventListener('change', updateNavLegend);
}

// Read every control's current DOM value and apply it. Browsers (Firefox in
// particular) restore form values on refresh without firing events, which
// would leave opts/camera out of sync; dispatching the events fixes that.
function syncControls() {
  const panel = $('panel');
  panel.querySelectorAll('select, input[type=color]').forEach(el => el.dispatchEvent(new Event(el.tagName==='SELECT'?'change':'input')));
  panel.querySelectorAll('input[type=range]').forEach(el => el.dispatchEvent(new Event('input')));
  panel.querySelectorAll('input[type=checkbox]').forEach(el => el.dispatchEvent(new Event('change')));
}
function updateNavLegend() {
  const m = $('nav-mode').value;
  const L = {
    turntable: 'LMB Orbit · MMB/Shift Pan · Wheel Zoom · DblClick Center',
    trackball: 'LMB Rotate · MMB/Shift Pan · Wheel Zoom · DblClick Center',
    fps: 'WASD Move · E/Q Up/Down · LMB Look · DblClick Center',
    fly: 'WASD Move · E/Q Roll · LMB Look · DblClick Center',
  };
  $('nav-legend').textContent = L[m];
}

boot();
