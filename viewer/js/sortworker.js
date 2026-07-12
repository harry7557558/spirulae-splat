// Depth-sort worker. Runs its own instance of the ssv_wasm module holding a
// copy of splat positions, so the native 16-bit counting sort (ssw_sort in
// viewer.cpp — same code path as ssv_sort) runs off the main thread and
// look-around never blocks rendering on multi-million-splat models. If the
// WASM module fails to load, falls back to an equivalent pure-JS sort.
//
// Protocol:
//   {type:'init', gen, count, positions: ArrayBuffer}   (x,y,z,opacity)*N f32
//   {type:'sort', gen, mode, c:[3], d:[3]}
//     -> {type:'sorted', gen, order: ArrayBuffer}       (empty if stale)
//   {type:'recycle', buffer: ArrayBuffer}                return for reuse
//
// Messages are handled strictly in order (serialized behind the module load),
// so an 'init' is always applied before the 'sort's that follow it.

let Module = null;       // wasm instance, or null -> JS fallback
let wposPtr = 0;         // heap pointer to the positions copy (wasm path)
let pos = null;          // Float32Array positions (JS fallback path)
let N = 0, gen = -1;
let key = null, cnt = null;  // JS fallback scratch
const pool = [];

const boot = (async () => {
  try {
    const factory = (await import('./ssv_wasm.js')).default;
    Module = await factory({
      locateFile: (p) => new URL(p, import.meta.url).href,
    });
  } catch (e) {
    Module = null;  // JS fallback
  }
})();

let queue = boot;
onmessage = (e) => { queue = queue.then(() => handle(e.data)); };

function handle(m) {
  if (m.type === 'init') {
    N = m.count; gen = m.gen;
    pool.length = 0;
    const src = new Float32Array(m.positions);
    if (Module) {
      wposPtr = Module.ccall('ssw_init', 'number', ['number'], [N]) >>> 0;
      if (wposPtr) {
        new Float32Array(Module.HEAPF32.buffer, wposPtr, N * 4).set(src);
        pos = null;  // heap copy owns the data now
        return;
      }
      // allocation failed (shouldn't happen — positions are N*16 bytes)
    }
    pos = src;
    key = new Uint16Array(N);
    cnt = new Uint32Array(65537);
    return;
  }
  if (m.type === 'recycle') {
    if (N && m.buffer.byteLength === N * 4) pool.push(m.buffer);
    return;
  }
  if (m.type !== 'sort') return;
  if (m.gen !== gen || N === 0 || (!pos && !wposPtr)) {
    // always reply so the main thread's in-flight flag clears
    postMessage({ type: 'sorted', gen: m.gen, order: new ArrayBuffer(0) }, [new ArrayBuffer(0)]);
    return;
  }

  const order = new Uint32Array(pool.pop() || new ArrayBuffer(N * 4));
  if (Module && wposPtr) {
    const p = Module.ccall('ssw_sort', 'number',
      ['number', 'number', 'number', 'number', 'number', 'number', 'number'],
      [m.mode | 0, m.c[0], m.c[1], m.c[2], m.d[0], m.d[1], m.d[2]]) >>> 0;
    order.set(new Uint32Array(Module.HEAPU32.buffer, p, N));
  } else {
    jsSort(m.mode | 0, m.c, m.d, order);
  }
  postMessage({ type: 'sorted', gen, order: order.buffer, native: !!(Module && wposPtr) }, [order.buffer]);
}

// Pure-JS replica of ssw_sort, used only if the WASM module failed to load.
function jsSort(mode, c, d, order) {
  const cx = c[0], cy = c[1], cz = c[2];
  const dx = d[0], dy = d[1], dz = d[2];
  const P = pos;
  // depth per ssv_sort: 0 planar, 1 distance^2 (equirect),
  // 2 sqrt(sqrt(z^4 + (x^2+y^2)^2/512)) (fisheye/equisolid)
  const depth = (i) => {
    const x = P[i * 4], y = P[i * 4 + 1], z = P[i * 4 + 2];
    if (mode === 0) return x * dx + y * dy + z * dz;
    const rx = x - cx, ry = y - cy, rz = z - cz;
    const r2 = rx * rx + ry * ry + rz * rz;
    if (mode === 1) return r2;
    const zc = rx * dx + ry * dy + rz * dz;
    let p2 = r2 - zc * zc; if (p2 < 0) p2 = 0;
    const z2 = zc * zc;
    return Math.sqrt(Math.sqrt(z2 * z2 + p2 * p2 * (1 / 512)));
  };

  let mn = Infinity, mx = -Infinity;
  for (let i = 0; i < N; i++) {
    const dep = depth(i);
    if (dep < mn) mn = dep;
    if (dep > mx) mx = dep;
  }
  let range = mx - mn; if (range < 1e-12) range = 1e-12;
  const B = 65536, scale = (B - 1) / range;
  cnt.fill(0);
  for (let i = 0; i < N; i++) {
    let q = ((depth(i) - mn) * scale) | 0;
    if (q < 0) q = 0; else if (q >= B) q = B - 1;
    const r = (B - 1) - q;          // reverse: farthest first
    key[i] = r;
    cnt[r + 1]++;
  }
  for (let i = 0; i < B; i++) cnt[i + 1] += cnt[i];
  for (let i = 0; i < N; i++) order[cnt[key[i]]++] = i;
}
