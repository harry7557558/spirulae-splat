// Camera + navigation controller. Ported/adapted from the training
// viewer (viewer.html): quaternion pose, orbit/pan/dolly/roll, and four modes
// (turntable / trackball / fps / fly). Renderer-agnostic — it just maintains a
// world-space pose that the renderer reads.

import { v3, quat, mat3 } from './linalg.js';

export class Camera {
  constructor() {
    this.pos = [0, 0, 3];
    this.rot = quat.identity();          // camera-to-world; camera looks down -Z
    this.target = [0, 0, 0];
    this.fov = 90 * Math.PI / 180;       // vertical FOV (radians)
    this.model = 'perspective';          // 'perspective' | 'orthographic'
  }
  right()   { return quat.rotVec(this.rot, [1, 0, 0]); }
  up()      { return quat.rotVec(this.rot, [0, 1, 0]); }
  forward() { return quat.rotVec(this.rot, [0, 0, -1]); }
  // cam-to-world rotation as mat3 (col-major)
  rotMat3() { return quat.toMat3(this.rot); }
  // world-to-camera rotation (transpose of cam-to-world)
  viewR()   { return mat3.transpose(quat.toMat3(this.rot)); }
  dist()    { return v3.len(v3.sub(this.pos, this.target)); }
  // orient the camera at `pos` looking toward `target`, with world up `up`
  lookAt(pos, target, up) {
    this.pos = pos.slice(); this.target = target.slice();
    const fwd = v3.norm(v3.sub(target, pos));           // camera -Z
    let right = v3.cross(fwd, up);
    if (v3.len(right) < 1e-6) right = v3.cross(fwd, [0,1,0]);
    right = v3.norm(right);
    const camUp = v3.cross(right, fwd);
    // cam-to-world matrix columns: right, up, -forward (col-major)
    const m = [ right[0],right[1],right[2], camUp[0],camUp[1],camUp[2], -fwd[0],-fwd[1],-fwd[2] ];
    this.rot = quat.fromMat3(m);
  }
}

export class Nav {
  constructor(camera) {
    this.cam = camera;
    this.mode = 'turntable';
    this._speedExp = 0;                  // 10^exp scale
    this._keys = new Set();
    this._pointers = new Map();
    this.onChange = () => {};
  }
  // Speed multiplier, identical to the training viewer (pure 10^exp).
  // Scene size is applied exactly once, in linear translations only, via
  // `_sceneScale` — never here (it would leak into the exponential orbit zoom).
  get speed() { return Math.pow(10, this._speedExp); }
  set speedExp(v) { this._speedExp = v; }
  _sceneScale = 1;

  // Fit the camera to a bounding sphere (robust center + radius from the
  // model, e.g. median-based so splat/mesh outliers don't blow up framing).
  fitSphere(c, r, up = [0,0,1]) {
    r = Math.max(r, 1e-6);
    this._sceneScale = r;
    const dist = r / Math.tan(this.cam.fov/2) * 1.25;
    // azimuth 35°, elevation 22° above the horizontal plane (perpendicular to up)
    const az = 35*Math.PI/180, el = 22*Math.PI/180;
    // build a basis where `up` is vertical
    let f = v3.norm(up);
    let a = Math.abs(f[2]) < 0.9 ? [0,0,1] : [1,0,0];
    let e0 = v3.norm(v3.cross(a, f));
    let e1 = v3.cross(f, e0);
    const dir = v3.add(v3.add(v3.scale(e0, Math.cos(el)*Math.cos(az)),
                               v3.scale(e1, Math.cos(el)*Math.sin(az))),
                       v3.scale(f, Math.sin(el)));
    this.cam.lookAt(v3.add(c, v3.scale(dir, dist)), c, up);
    this.onChange();
  }
  reset() { this.fitSphere([0,0,0], 1); }

  // --- navigation, ported verbatim from the training viewer (viewer.html);
  // the only additions are `_sceneScale` folded into LINEAR translations so
  // the feel is scale-invariant for arbitrarily sized models (the training
  // viewer works on pre-normalized unit-scale scenes). World up is +Z in this
  // canonical space; the Y-up/Z-up toggle rotates the model into it, so the
  // turntable axis follows the chosen up direction. ---
  orbit(dx, dy) {
    const sensitivity = 0.005;
    const right = this.cam.right();
    const up = this.mode === 'trackball' ? this.cam.up() : [0,0,1];
    const qa = quat.fromAxisAngle(up, -dx * sensitivity);
    const qb = quat.fromAxisAngle(right, -dy * sensitivity);
    this.cam.rot = quat.norm(quat.mul(quat.mul(qa, qb), this.cam.rot));
    const dist = v3.len(v3.sub(this.cam.pos, this.cam.target));
    const newFwd = quat.rotVec(this.cam.rot, [0,0,1]);
    this.cam.pos = v3.add(this.cam.target, v3.scale(newFwd, dist));
    this.onChange();
  }
  look(dx, dy) {
    const s = 0.003;
    const qa = quat.fromAxisAngle(this.cam.up(), -dx * s);
    const qb = quat.fromAxisAngle(this.cam.right(), -dy * s);
    if (this.mode === 'fps') {
      this.cam.rot = quat.norm(quat.mul(this.cam.rot, quat.fromAxisAngle([0,1,0], -dx*s)));
      this.cam.rot = quat.norm(quat.mul(qb, this.cam.rot));
    } else {
      this.cam.rot = quat.norm(quat.mul(quat.mul(qa, qb), this.cam.rot));
    }
    this.onChange();
  }
  pan(dx, dy) {
    const s = this.speed * 0.002 * this._sceneScale;
    const delta = v3.add(v3.scale(this.cam.right(), -dx*s), v3.scale(this.cam.up(), dy*s));
    this.cam.pos = v3.add(this.cam.pos, delta);
    this.cam.target = v3.add(this.cam.target, delta);
    this.onChange();
  }
  dolly(delta) {
    const s = this.speed * 0.2;
    if (this.mode === 'trackball' || this.mode === 'turntable') {
      this.cam.pos = v3.add(v3.scale(v3.sub(this.cam.pos, this.cam.target), Math.exp(delta * 0.004 * s)), this.cam.target);
    } else {
      const fwd = this.cam.forward();
      const d = -delta * 0.004 * s * this._sceneScale;
      this.cam.pos = v3.add(this.cam.pos, v3.scale(fwd, d));
      this.cam.target = v3.add(this.cam.target, v3.scale(fwd, d));
    }
    this.onChange();
  }
  roll(delta) {
    const q = quat.fromAxisAngle(this.cam.forward(), delta);
    this.cam.rot = quat.norm(quat.mul(q, this.cam.rot));
    this.onChange();
  }
  drag(dx, dy, panning) {
    if (panning) this.pan(dx, dy);
    else if (this.mode === 'turntable' || this.mode === 'trackball') this.orbit(dx, dy);
    else this.look(dx, dy);
  }

  // ---- keyboard / gamepad ticks (ported from the training viewer) ----
  keyDown(k) { this._keys.add(k); }
  keyUp(k) { this._keys.delete(k); }
  tick(dt) {
    const s = this.speed * dt * this._sceneScale;
    const sr = dt * 1.0;
    const K = this._keys;
    const fwd = this.cam.forward(), right = this.cam.right();
    const up = (this.mode === 'turntable' || this.mode === 'fps') ? [0,0,1] : this.cam.up();
    let move = [0,0,0], rotate = 0.0;
    if (K.has('w')||K.has('arrowup')) move = v3.add(move, v3.scale(fwd, s));
    if (K.has('s')||K.has('arrowdown')) move = v3.add(move, v3.scale(fwd, -s));
    if (K.has('a')||K.has('arrowleft')) move = v3.add(move, v3.scale(right, -s));
    if (K.has('d')||K.has('arrowright')) move = v3.add(move, v3.scale(right, s));
    if (this.mode === 'fly') {
      if (K.has('e')) rotate -= sr;
      if (K.has('q')) rotate += sr;
    } else {
      if (K.has('e')) move = v3.add(move, v3.scale(up, s));
      if (K.has('q')) move = v3.add(move, v3.scale(up, -s));
    }
    let moved = false;
    if (v3.len(move) > 1e-12) { this.cam.pos = v3.add(this.cam.pos, move); this.cam.target = v3.add(this.cam.target, move); moved = true; }
    if (rotate !== 0.0) { this.roll(rotate); moved = true; }
    if (moved) this.onChange();
    return moved;
  }
  // Ported verbatim from the training viewer: left stick moves, right stick
  // looks, the two front buttons (triggers, buttons 6/7) are always roll, and
  // while moving they double as down/up. `_sceneScale` scales translations.
  gamepadTick(dt) {
    const gamepads = navigator.getGamepads ? navigator.getGamepads() : [];
    let moved = false;
    for (const gp of gamepads) {
      if (!gp || !gp.connected) continue;
      const lx = gp.axes[0] || 0, ly = gp.axes[1] || 0; // left stick
      const rx = gp.axes[2] || 0, ry = gp.axes[3] || 0; // right stick
      const rollInput = (gp.buttons[6] ? gp.buttons[6].value : 0) - (gp.buttons[7] ? gp.buttons[7].value : 0);
      const deadzone = 0.1;
      // Movement
      if (Math.abs(lx) > deadzone || Math.abs(ly) > deadzone) {
        const s = this.speed * dt * 1.0 * this._sceneScale;
        const fwd = this.cam.forward();
        const right = this.cam.right();
        const up = [0,0,1];
        let move = v3.add(v3.scale(right, lx * s), v3.scale(fwd, -ly * s));
        // Triggers for up/down
        const lt = gp.buttons[6] ? gp.buttons[6].value : 0;
        const rt = gp.buttons[7] ? gp.buttons[7].value : 0;
        move = v3.add(move, v3.scale(up, (rt - lt) * s));
        this.cam.pos = v3.add(this.cam.pos, move);
        this.cam.target = v3.add(this.cam.target, move);
        moved = true;
      }
      // Look
      if (Math.abs(rx) > deadzone || Math.abs(ry) > deadzone) {
        this.look(rx * dt * 400, ry * dt * 400);
        moved = true;
      }
      // Roll (always available on the front buttons)
      if (Math.abs(rollInput) > deadzone) {
        const sr = dt * 1.0;
        const q = quat.fromAxisAngle(this.cam.forward(), rollInput * sr);
        this.cam.rot = quat.norm(quat.mul(q, this.cam.rot));
        moved = true;
      }
    }
    if (moved) this.onChange();
    return moved;
  }
}
