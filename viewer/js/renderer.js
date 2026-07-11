// WebGL2 renderer: 3D Gaussian splats (instanced, HDR + tonemap) and meshes
// (depth-tested), plus a screen-space grid/axes overlay. Only one model is
// resident at a time; loading another frees the previous GPU resources.

import { SPLAT_VS, SPLAT_FS, SPLAT_FS_3D, PICK_FS, PICK_FS_3D, FULLSCREEN_VS, GRID_FS, TONEMAP_FS, MESH_VS, MESH_FS } from './shaders.js';
import { mat3, mat4, quat } from './linalg.js';

function compile(gl, type, src) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS))
    throw new Error('shader compile: ' + gl.getShaderInfoLog(s) + '\n' + src.split('\n').slice(0,40).join('\n'));
  return s;
}
function link(gl, vs, fs) {
  const p = gl.createProgram();
  gl.attachShader(p, compile(gl, gl.VERTEX_SHADER, vs));
  gl.attachShader(p, compile(gl, gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS))
    throw new Error('program link: ' + gl.getProgramInfoLog(p));
  return p;
}
function uniforms(gl, prog) {
  const u = {};
  const n = gl.getProgramParameter(prog, gl.ACTIVE_UNIFORMS);
  for (let i=0;i<n;i++) { const info = gl.getActiveUniform(prog, i); u[info.name] = gl.getUniformLocation(prog, info.name); }
  return u;
}
// inject "#define K V" lines right after the #version directive
function withDefines(src, defs) {
  const lines = Object.entries(defs).map(([k, v]) => `#define ${k} ${v}`).join('\n');
  return src.replace(/^(#version[^\n]*\n)/, `$1${lines}\n`);
}

// choose a (W,H,L) layout for `count` texels. Layers are capped at 2048x2048
// (4M texels) so uploads never need a padded full-size copy, growing toward
// maxTex only if the hardware layer limit would be exceeded.
function layout(count, maxTex, maxLayers = 256) {
  count = Math.max(count, 1);
  const capW = Math.min(2048, maxTex);
  if (count <= capW * capW) {          // small model: single compact layer
    let W = 1; while (W*W < count) W <<= 1;
    W = Math.min(W, maxTex);
    return { W, H: Math.max(Math.ceil(count / W), 1), L: 1 };
  }
  let W = capW, H = capW;
  while (Math.ceil(count / (W*H)) > maxLayers && (W < maxTex || H < maxTex)) {
    if (H <= W) H = Math.min(H*2, maxTex); else W = Math.min(W*2, maxTex);
  }
  return { W, H, L: Math.ceil(count / (W*H)) };
}

// Allocate a (W,H,L) array texture and upload `count` texels straight from a
// (typically WASM-heap) view — full layers, then full rows, then the row tail.
// No padded intermediate copy: for multi-GB models that copy would double
// peak memory. RGB16F/RGB formats upload 3-component data as-is.
function uploadLayered(gl, view, count, comps, internalFormat, format, type, W, H, L) {
  const tex = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D_ARRAY, tex);
  gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
  gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D_ARRAY, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
  gl.texImage3D(gl.TEXTURE_2D_ARRAY, 0, internalFormat, W, H, L, 0, format, type, null);
  const perLayer = W*H;
  let done = 0, l = 0;
  while (done < count) {
    const inLayer = Math.min(perLayer, count - done);
    const rows = Math.floor(inLayer / W);
    if (rows > 0)
      gl.texSubImage3D(gl.TEXTURE_2D_ARRAY, 0, 0, 0, l, W, rows, 1, format, type,
                       view.subarray(done*comps, (done + rows*W)*comps));
    const tail = inLayer - rows*W;
    if (tail > 0)
      gl.texSubImage3D(gl.TEXTURE_2D_ARRAY, 0, 0, rows, l, tail, 1, 1, format, type,
                       view.subarray((done + rows*W)*comps, (done + inLayer)*comps));
    done += inLayer; l++;
  }
  return tex;
}

export class Renderer {
  constructor(canvas) {
    const gl = canvas.getContext('webgl2', {
      antialias: false, alpha: false, premultipliedAlpha: false,
      preserveDrawingBuffer: false, powerPreference: 'high-performance',
    });
    if (!gl) throw new Error('WebGL2 not available');
    this.gl = gl; this.canvas = canvas;
    if (!gl.getExtension('EXT_color_buffer_float'))
      console.warn('EXT_color_buffer_float unavailable — HDR splat rendering may fail');
    gl.getExtension('OES_texture_float_linear');
    this.maxTex = Math.min(gl.getParameter(gl.MAX_TEXTURE_SIZE), 16384);

    // splat/grid/mesh programs are compiled lazily per (primitive, camera
    // model) combination — see _prog(); only the tonemap pass is fixed
    this._progs = {};
    this.pTone = link(gl, FULLSCREEN_VS, TONEMAP_FS); this.uTone = uniforms(gl, this.pTone);

    this.emptyVao = gl.createVertexArray();
    this.hdr = null; this.hdrTex = null; this.hdrW = 0; this.hdrH = 0;
    this.model = null;   // { type, ... }
  }

  // compile-time specialized program cache
  _prog(name, vsSrc, fsSrc, defs) {
    const key = name + ':' + Object.entries(defs).map(([k,v])=>k+'='+v).join(',');
    let p = this._progs[key];
    if (!p) {
      const prog = link(this.gl, withDefines(vsSrc, defs), withDefines(fsSrc, defs));
      p = this._progs[key] = { prog, u: uniforms(this.gl, prog) };
    }
    return p;
  }

  // ---- resource management ----
  _freeModel() {
    const gl = this.gl, m = this.model;
    if (!m) return;
    if (m.type === 'splat') {
      for (const t of [m.posT,m.rotT,m.scaleT,m.colT,...(m.shT||[])]) if (t) gl.deleteTexture(t);
      gl.deleteBuffer(m.indexBuf); gl.deleteVertexArray(m.vao);
    } else if (m.type === 'mesh') {
      for (const b of [m.posB,m.nrmB,m.colB,m.uvB]) if (b) gl.deleteBuffer(b);
      for (const p of (m.idxParts||[])) gl.deleteBuffer(p.buf);
      if (m.tex) gl.deleteTexture(m.tex);
      gl.deleteVertexArray(m.vao);
    }
    this.model = null;
  }

  // ---- splat model ----
  // data.posop is f32; quat/scale/fdc/sh are half floats (Uint16Array views)
  // and upload directly as HALF_FLOAT with no conversion or padding.
  setSplat(data) {
    this._freeModel();
    const gl = this.gl;
    const N = data.count;
    const maxLayers = gl.getParameter(gl.MAX_ARRAY_TEXTURE_LAYERS);
    const lay = layout(N, this.maxTex, maxLayers);
    const posT   = uploadLayered(gl, data.posop, N, 4, gl.RGBA32F, gl.RGBA, gl.FLOAT,      lay.W, lay.H, lay.L);
    const rotT   = uploadLayered(gl, data.quat,  N, 4, gl.RGBA16F, gl.RGBA, gl.HALF_FLOAT, lay.W, lay.H, lay.L);
    const scaleT = uploadLayered(gl, data.scale, N, 3, gl.RGB16F,  gl.RGB,  gl.HALF_FLOAT, lay.W, lay.H, lay.L);
    const colT   = uploadLayered(gl, data.fdc,   N, 3, gl.RGB16F,  gl.RGB,  gl.HALF_FLOAT, lay.W, lay.H, lay.L);

    // SH is by far the largest attribute; split it across up to 4 textures by
    // layer ranges so no single GL resource exceeds ~512MB (SwiftShader's
    // Vulkan and D3D11 both cap per-resource allocations).
    let shT = null, shLay = { W:1, H:1, L:1 }, shChunkLayers = 1;
    if (data.shRest > 0 && data.sh && data.sh.length) {
      const K = data.shRest, M = N*K;
      shLay = layout(M, this.maxTex, maxLayers);
      const layerBytes = shLay.W * shLay.H * 8;          // RGB16F is stored padded
      const nChunks = Math.min(4, Math.max(1, Math.ceil(shLay.L * layerBytes / (512 << 20))));
      shChunkLayers = Math.ceil(shLay.L / nChunks);
      shT = [];
      for (let l0 = 0; l0 < shLay.L; l0 += shChunkLayers) {
        const Lc = Math.min(shChunkLayers, shLay.L - l0);
        const t0 = l0 * shLay.W * shLay.H;
        const cnt = Math.min(M - t0, shLay.W * shLay.H * Lc);
        shT.push(uploadLayered(gl, data.sh.subarray(t0*3), cnt, 3,
                               gl.RGB16F, gl.RGB, gl.HALF_FLOAT, shLay.W, shLay.H, Lc));
      }
    }

    // sorted-index instanced buffer (initialized 0..N-1)
    const idx = new Uint32Array(N); for (let i=0;i<N;i++) idx[i]=i;
    const indexBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, indexBuf);
    gl.bufferData(gl.ARRAY_BUFFER, idx, gl.DYNAMIC_DRAW);
    const vao = gl.createVertexArray();
    gl.bindVertexArray(vao);
    gl.bindBuffer(gl.ARRAY_BUFFER, indexBuf);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribIPointer(0, 1, gl.UNSIGNED_INT, 0, 0);
    gl.vertexAttribDivisor(0, 1);
    gl.bindVertexArray(null);

    this.model = {
      type:'splat', count:N, shDegree:data.shDegree, shRest:data.shRest,
      posT, rotT, scaleT, colT, shT, shChunkLayers, lay, shLay, indexBuf, vao,
    };
  }

  updateSplatOrder(orderU32) {
    if (!this.model || this.model.type !== 'splat') return;
    const gl = this.gl;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.model.indexBuf);
    gl.bufferSubData(gl.ARRAY_BUFFER, 0, orderU32);
  }

  // ---- mesh model ----
  setMesh(data, texImage) {
    this._freeModel();
    const gl = this.gl;
    const vao = gl.createVertexArray();
    gl.bindVertexArray(vao);
    const mkBuf = (arr, loc, comps) => {
      const b = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, b);
      gl.bufferData(gl.ARRAY_BUFFER, arr, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(loc);
      gl.vertexAttribPointer(loc, comps, gl.FLOAT, false, 0, 0);
      return b;
    };
    const posB = mkBuf(data.pos, 0, 3);
    let nrmB=null, colB=null, uvB=null;
    if (data.normal) nrmB = mkBuf(data.normal, 1, 3);
    else { gl.disableVertexAttribArray(1); gl.vertexAttrib3f(1, 0, 0, 0); }  // zero => FS uses face normals
    if (data.color) {
      // normalized ubyte attribute — no float conversion copy (large meshes)
      colB = gl.createBuffer();
      gl.bindBuffer(gl.ARRAY_BUFFER, colB);
      gl.bufferData(gl.ARRAY_BUFFER, data.color, gl.STATIC_DRAW);
      gl.enableVertexAttribArray(2);
      gl.vertexAttribPointer(2, 3, gl.UNSIGNED_BYTE, true, 0, 0);
    } else { gl.disableVertexAttribArray(2); gl.vertexAttrib3f(2, 0.78, 0.78, 0.78); }
    if (data.uv) uvB = mkBuf(data.uv, 3, 2);
    else { gl.disableVertexAttribArray(3); gl.vertexAttrib2f(3, 0, 0); }

    // split the index buffer into ~256MB chunks (per-resource size limits);
    // the chunk size MUST be a multiple of 3 to not split a triangle
    const idxParts = [];
    const CHUNK = 63 << 20;   // indices per part (252MB), divisible by 3
    for (let i0 = 0; i0 < data.idx.length; i0 += CHUNK) {
      const part = data.idx.subarray(i0, Math.min(i0 + CHUNK, data.idx.length));
      const buf = gl.createBuffer();
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, buf);
      gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, part, gl.STATIC_DRAW);
      idxParts.push({ buf, count: part.length });
    }
    gl.bindVertexArray(null);

    let tex = null, mode = 0;
    if (texImage && data.uv) {
      tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texImage);
      gl.generateMipmap(gl.TEXTURE_2D);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      mode = 2;
    } else if (data.color) mode = 1;

    this.model = { type:'mesh', nt:data.nt, vao, posB, nrmB, colB, uvB, idxParts, tex, mode };
  }

  // ---- HDR target ----
  _ensureHdr(w, h) {
    const gl = this.gl;
    if (this.hdrW === w && this.hdrH === h && this.hdr) return;
    if (this.hdr) { gl.deleteFramebuffer(this.hdr); gl.deleteTexture(this.hdrTex); }
    this.hdrTex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, this.hdrTex);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA16F, w, h, 0, gl.RGBA, gl.HALF_FLOAT, null);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    this.hdr = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.hdr);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.hdrTex, 0);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    this.hdrW = w; this.hdrH = h;
  }

  // camera params -> intrinsics
  _intrinsics(cam, w, h) {
    let fx, fy;
    const half = Math.min(w, h) / 2;
    if (cam.model === 'orthographic') {
      const halfH = cam.dist() * Math.tan(cam.fov/2);
      fy = h / (2*halfH); fx = fy;
    } else if (cam.model === 'fisheye') {           // equidistant: r = f*theta
      fx = fy = half / Math.max(cam.fov/2, 1e-3);
    } else if (cam.model === 'equisolid') {         // r = 2 f sin(theta/2)
      fx = fy = half / (2*Math.sin(cam.fov/4) + 1e-6);
    } else if (cam.model === 'equirectangular') {   // full 360x180
      fx = w / (2*Math.PI); fy = h / Math.PI;
    } else {                                         // perspective (pinhole)
      fy = (h/2) / Math.tan(cam.fov/2); fx = fy;
    }
    return { fx, fy, cx: w/2, cy: h/2 };
  }

  render(cam, opts) {
    const gl = this.gl;
    const w = this.canvas.width, h = this.canvas.height;
    if (w === 0 || h === 0) return;
    const cameraModel = CAM_MODEL[cam.model] ?? 0;
    const intr = this._intrinsics(cam, w, h);
    const viewR = cam.viewR ? cam.viewR() : mat3.transpose(quat.toMat3(cam.rot));
    const camPos = cam.pos;
    // up-axis world transform U (native -> canonical z-up)
    const U = opts.upAxis === 'y' ? [1,0,0, 0,0,1, 0,-1,0] : [1,0,0, 0,1,0, 0,0,1];

    const bg = opts.background;
    gl.viewport(0, 0, w, h);

    if (this.model && this.model.type === 'splat') {
      this._ensureHdr(w, h);
      gl.bindFramebuffer(gl.FRAMEBUFFER, this.hdr);
      gl.clearColor(0, 0, 0, 1);   // alpha = transmittance (starts fully transparent bg)
      gl.clear(gl.COLOR_BUFFER_BIT);
      gl.disable(gl.DEPTH_TEST);

      if (opts.showGrid) this._drawGrid(viewR, camPos, U, intr, cameraModel, true);

      // splat pass (over blending; alpha accumulates coverage)
      gl.enable(gl.BLEND);
      gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA, gl.ZERO, gl.ONE_MINUS_SRC_ALPHA);
      this._drawSplats(viewR, camPos, U, intr, cameraModel, w, h, opts);
      gl.disable(gl.BLEND);

      // tonemap HDR -> canvas
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.viewport(0, 0, w, h);
      gl.useProgram(this.pTone);
      gl.activeTexture(gl.TEXTURE0);
      gl.bindTexture(gl.TEXTURE_2D, this.hdrTex);
      gl.uniform1i(this.uTone.uHdr, 0);
      gl.uniformMatrix3fv(this.uTone.uGamut, false, opts.gamut);
      gl.uniform1i(this.uTone.uIsLinear, opts.isLinear ? 1 : 0);
      gl.uniform3fv(this.uTone.uBackground, bg);
      gl.uniform1f(this.uTone.uExposure, opts.exposure ?? 1.0);
      gl.bindVertexArray(this.emptyVao);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      gl.bindVertexArray(null);
    } else if (this.model && this.model.type === 'mesh') {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.clearColor(bg[0], bg[1], bg[2], 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      if (opts.showGrid) { gl.disable(gl.DEPTH_TEST); this._drawGrid(viewR, camPos, U, intr, cameraModel, false); }
      this._drawMesh(cam, viewR, camPos, U, w, h, opts);
    } else {
      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
      gl.clearColor(bg[0], bg[1], bg[2], 1);
      gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
      if (opts.showGrid) { gl.disable(gl.DEPTH_TEST); this._drawGrid(viewR, camPos, U, intr, cameraModel, false); }
    }
  }

  _drawGrid(viewR, camPos, U, intr, cameraModel, intoHdr) {
    const gl = this.gl;
    const { prog, u } = this._prog('grid', FULLSCREEN_VS, GRID_FS, { CAM_MODEL: cameraModel });
    gl.useProgram(prog);
    gl.uniformMatrix3fv(u.uViewR, false, viewR);
    gl.uniform3fv(u.uCamPos, camPos);
    // grid/axes are evaluated in the model's native frame (Y-up/Z-up aware)
    gl.uniformMatrix3fv(u.uUpT, false, mat3.transpose(U));
    gl.uniform2f(u.uFocal, intr.fx, intr.fy);
    gl.uniform2f(u.uCenter, intr.cx, intr.cy);
    gl.enable(gl.BLEND);
    // into HDR: alpha is transmittance -> multiply it down (ZERO, 1-src_a).
    gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA,
                         intoHdr?gl.ZERO:gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.bindVertexArray(this.emptyVao);
    gl.drawArrays(gl.TRIANGLES, 0, 3);
    gl.bindVertexArray(null);
    gl.disable(gl.BLEND);
  }

  // set the shared splat-program uniforms + texture binds (render and pick)
  _splatUniforms(u, viewR, camPos, U, intr, w, h, opts) {
    const gl = this.gl, m = this.model;
    // uViewR_native = viewR * U ; uCamPos_native = U^T * camPos
    const viewRnat = mat3.mul(viewR, U);
    const camPosNat = mat3.mulVec(mat3.transpose(U), camPos);
    gl.uniformMatrix3fv(u.uViewR, false, viewRnat);
    gl.uniform3fv(u.uCamPos, camPosNat);
    gl.uniform2f(u.uFocal, intr.fx, intr.fy);
    gl.uniform2f(u.uCenter, intr.cx, intr.cy);
    gl.uniform2f(u.uViewport, w, h);
    gl.uniform1f(u.uNear, 0.01);
    gl.uniform1f(u.uOpacityScale, opts.opacityScale ?? 1.0);
    gl.uniform1i(u.uShDegree, opts.shDegree != null ? Math.min(opts.shDegree, m.shDegree) : m.shDegree);
    gl.uniform1i(u.uShRest, m.shRest);
    gl.uniform2i(u.uDim, m.lay.W, m.lay.H);
    gl.uniform2i(u.uShDim, m.shLay.W, m.shLay.H);
    gl.uniform1i(u.uShChunkLayers, Math.max(m.shChunkLayers || 1, 1));
    const bind = (loc, tex, unit) => { gl.activeTexture(gl.TEXTURE0+unit); gl.bindTexture(gl.TEXTURE_2D_ARRAY, tex); gl.uniform1i(loc, unit); };
    bind(u.uPos, m.posT, 0); bind(u.uRot, m.rotT, 1); bind(u.uScale, m.scaleT, 2); bind(u.uCol, m.colT, 3);
    const sh = m.shT || [];
    for (let c = 0; c < 4; c++) bind(u['uSh'+c], sh[c] || sh[0] || m.colT, 4+c);
  }

  _drawSplats(viewR, camPos, U, intr, cameraModel, w, h, opts) {
    const gl = this.gl, m = this.model;
    const prim = opts.primitive | 0;
    const { prog, u } = this._prog('splat', SPLAT_VS, prim === 2 ? SPLAT_FS_3D : SPLAT_FS,
                                   { CAM_MODEL: cameraModel, PRIM: prim });
    gl.useProgram(prog);
    this._splatUniforms(u, viewR, camPos, U, intr, w, h, opts);
    gl.bindVertexArray(m.vao);
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, m.count);
    gl.bindVertexArray(null);
  }

  // Depth under a canvas pixel for the current splat model (double-click view
  // centering): render the splats once into a 1x1 float target with the
  // viewport shifted so pixel (px,py) lands on it, accumulating premultiplied
  // depth over the sorted order. Returns projCam-convention depth, or null on
  // a background click. O(1) memory, one instanced draw.
  pickSplatDepth(cam, opts, px, py) {
    const gl = this.gl, m = this.model;
    if (!m || m.type !== 'splat') return null;
    const w = this.canvas.width, h = this.canvas.height;
    const cameraModel = CAM_MODEL[cam.model] ?? 0;
    const prim = opts.primitive | 0;
    const intr = this._intrinsics(cam, w, h);
    const viewR = cam.viewR();
    const U = opts.upAxis === 'y' ? [1,0,0, 0,0,1, 0,-1,0] : [1,0,0, 0,1,0, 0,0,1];

    if (!this.pickFbo) {
      this.pickTex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, this.pickTex);
      gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA32F, 1, 1, 0, gl.RGBA, gl.FLOAT, null);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      this.pickFbo = gl.createFramebuffer();
      gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickFbo);
      gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, this.pickTex, 0);
    }
    const ox = Math.floor(px), oy = Math.floor(py);
    const { prog, u } = this._prog('pick', SPLAT_VS, prim === 2 ? PICK_FS_3D : PICK_FS,
                                   { CAM_MODEL: cameraModel, PRIM: prim });
    gl.bindFramebuffer(gl.FRAMEBUFFER, this.pickFbo);
    gl.viewport(-ox, -oy, w, h);        // canvas pixel (ox,oy) -> FBO pixel (0,0)
    gl.clearColor(0, 0, 0, 0);
    gl.clear(gl.COLOR_BUFFER_BIT);
    gl.disable(gl.DEPTH_TEST);
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);   // premultiplied over, back-to-front
    gl.useProgram(prog);
    this._splatUniforms(u, viewR, cam.pos, U, intr, w, h, opts);
    gl.uniform2f(u.uPixelOffset, ox, oy);
    gl.bindVertexArray(m.vao);
    gl.drawArraysInstanced(gl.TRIANGLE_STRIP, 0, 4, m.count);
    gl.bindVertexArray(null);
    gl.disable(gl.BLEND);
    const buf = new Float32Array(4);
    gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.FLOAT, buf);
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    const cov = buf[3];
    if (!(cov > 0.1)) return null;      // background / negligible coverage
    return buf[0] / cov;                // alpha-weighted expected depth
  }

  _drawMesh(cam, viewR, camPos, U, w, h, opts) {
    const gl = this.gl, m = this.model;
    const cameraModel = CAM_MODEL[cam.model] ?? 0;
    const { prog, u } = this._prog('mesh', MESH_VS, MESH_FS, { CAM_MODEL: cameraModel });
    gl.enable(gl.DEPTH_TEST);
    gl.depthFunc(gl.LEQUAL);
    const near = Math.max(cam.dist()*0.001, 1e-3), far = cam.dist()*100 + 100;
    const view = mat4.view(mat3.transpose(viewR), camPos);
    const proj = cam.model === 'orthographic'
      ? (() => { const hh = cam.dist()*Math.tan(cam.fov/2), hw = hh*w/h; return mat4.ortho(-hw,hw,-hh,hh,near,far); })()
      : mat4.perspective(cam.fov, w/h, near, far);
    // model matrix = U (mat3 -> mat4)
    const model4 = new Float32Array([U[0],U[1],U[2],0, U[3],U[4],U[5],0, U[6],U[7],U[8],0, 0,0,0,1]);
    const mvp = mat4.mul(mat4.mul(proj, view), model4);
    const intr = this._intrinsics(cam, w, h);
    gl.useProgram(prog);
    gl.uniformMatrix4fv(u.uMVP, false, mvp);
    gl.uniformMatrix3fv(u.uNormalMat, false, U);
    gl.uniform1i(u.uMode, m.mode);
    gl.uniform1i(u.uShade, opts.shade ? 1 : 0);
    gl.uniform1i(u.uFlat, opts.flatShade ? 1 : 0);
    gl.uniform1i(u.uColorOn, opts.meshColor !== false ? 1 : 0);
    // uniforms for the nonlinear-camera path (fisheye / equirect)
    gl.uniformMatrix3fv(u.uModelR, false, U);
    gl.uniformMatrix3fv(u.uViewR, false, viewR);
    gl.uniform3fv(u.uCamPos, camPos);
    gl.uniform2f(u.uFocal, intr.fx, intr.fy);
    gl.uniform2f(u.uCenter, intr.cx, intr.cy);
    gl.uniform2f(u.uViewport, w, h);
    gl.uniform2f(u.uNearFar, near, far);
    if (m.mode === 2) { gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, m.tex); gl.uniform1i(u.uTex, 0); }
    gl.bindVertexArray(m.vao);
    for (const p of m.idxParts) {
      gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, p.buf);
      gl.drawElements(gl.TRIANGLES, p.count, gl.UNSIGNED_INT, 0);
    }
    gl.bindVertexArray(null);
    gl.disable(gl.DEPTH_TEST);
  }
}

// camera model string -> shader int (must match PROJ_GLSL in shaders.js)
const CAM_MODEL = { perspective:0, orthographic:1, fisheye:2, equisolid:3, equirectangular:4 };

Renderer.prototype.snapshot = function() {
  const gl = this.gl, w = this.canvas.width, h = this.canvas.height;
  const buf = new Uint8Array(w*h*4);
  gl.bindFramebuffer(gl.FRAMEBUFFER, null);
  gl.readPixels(0, 0, w, h, gl.RGBA, gl.UNSIGNED_BYTE, buf);
  let nonBlack = 0; const sum = [0,0,0];
  for (let i=0;i<buf.length;i+=4) { sum[0]+=buf[i]; sum[1]+=buf[i+1]; sum[2]+=buf[i+2];
    if (buf[i]>2||buf[i+1]>2||buf[i+2]>2) nonBlack++; }
  const n = w*h;
  return { w, h, nonBlack, frac: +(nonBlack/n).toFixed(3), avg: sum.map(s=>Math.round(s/n)) };
};
