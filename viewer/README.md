# Spirulae Standalone Viewer

A dependency-free, client-side **WebGL2** viewer for 3D Gaussian Splats and
meshes. It runs entirely in the browser (no server) and can be hosted as a
static site (e.g. GitHub Pages). Performance-critical parsing, depth sorting,
and statistics run in **C++ compiled to WebAssembly** (built with CMake +
Emscripten); rendering is WebGL2.

It is intentionally **independent** from the rest of `spirulae_splat` — nothing
here imports from the training code, so it will not break when that code
changes.

## Features

- **3D Gaussian Splatting** (`.ply`, INRIA/spirulae layout, including very large
  files — the file is streamed straight into the WASM heap, so no 2 GB JS
  `ArrayBuffer` limit).
  - Primitive select: **3DGS**, **Mip** (antialiased), **3DGUT** (unscented
    transform projection; fragments evaluate the 3D Gaussian along per-pixel
    rays, "eval3d"). 3DGS/Mip use analytic projection Jacobians with the
    conventional out-of-bounds Jacobian clipping.
  - **Color gamut** (Rec.709 / DCI-P3 / Rec.2020 / AdobeRGB / ACEScg /
    ACES2065-1) and a **linear-color** toggle, matching the training pipeline's
    `rgb_to_srgb` conventions.
  - Spherical harmonics up to degree 4, exposure.
  - GPU-friendly: attributes live in `TEXTURE_2D_ARRAY`s laid out so arbitrarily
    large models render **even under a small `MAX_TEXTURE_SIZE`** (tiles into
    array layers).
- **Meshes** (`.ply`, `.obj`, `.gltf`, `.glb`) as produced by the meshing code:
  vertex colors or a base-color texture atlas, shading toggles (shaded /
  unshaded, flat / interpolated normals, color on/off) with a view-following
  headlight so the surface reads from every angle. GLTF/GLB are parsed in JS;
  drop the model together with its external `.bin` / `.mtl` / image files for
  textured `.gltf` / `.obj`.
- Orbit / trackball / first-person / free-fly navigation (mouse, touch, keyboard,
  gamepad — matched to the training viewer), Y-up ↔ Z-up toggle (switching keeps
  the current view), and camera models **perspective / orthographic / fisheye
  (equidistant) / fisheye (equisolid) / equirectangular 360°** — all with
  analytic Jacobians / sigma-point projection; mesh triangles crossing a
  projection discontinuity (the equirect seam, the fisheye backward point) are
  discarded, and fragments outside a fisheye image circle are clipped.
- A distance-adaptive **axes + grid** overlay (stays legible at any zoom) and a
  configurable background.
- **Statistics** (splat count / vertices / edges / faces) and on-demand
  **parameter histograms** (opacity, scale, effective rank, RGB, anisotropy for
  splats; edge length, triangle area, coordinates for meshes), computed in WASM
  and cached.
- **Double-click** the viewport to center the view on the point under the
  cursor (MeshLab-style: the point becomes the orbit pivot and slides onto the
  optical axis). Splats use a one-pixel GPU depth pass; meshes use a WASM
  raycast. Works with every camera model.
- One model at a time — dropping another replaces it and frees the previous GPU
  buffers. Replacing keeps the current viewpoint (the camera is only fitted for
  the first model; refresh the page to start over).

## Build

Requires the Emscripten SDK on `PATH` (`emcc`, `emcmake`) and CMake ≥ 3.16.

```bash
source ~/emsdk/emsdk_env.sh      # activate emsdk
./build.sh                       # emcmake cmake + cmake --build
```

This produces `js/ssv_wasm.js` and `js/ssv_wasm.wasm` (committed so the site is
directly hostable without a build step).

## Run

Serve the directory over HTTP (ES modules + WASM require it — `file://` will not
work):

```bash
python3 -m http.server -d .      # then open http://localhost:8000/
```

Drag a model onto the canvas, or click to browse. You can also deep-link a
hosted model: `index.html?model=<url>` (external `.bin`/`.mtl`/image siblings are
fetched automatically).

## Layout

```
index.html          UI + panel
css/style.css        theme (adapted from the training viewer)
js/
  main.js            app wiring, input, render loop, sorting
  renderer.js        WebGL2 renderer (splat HDR pass, mesh, grid, tonemap)
  shaders.js         GLSL (splat / grid / tonemap / mesh)
  camera.js          camera + navigation modes
  wasm.js            WASM bridge + streaming loader + GLTF/GLB (JS)
  colors.js          gamut matrices
  linalg.js          vec/quat/mat helpers
  histogram.js       canvas bar chart
  ssv_wasm.{js,wasm} built WASM module
src/viewer.cpp       parsers (PLY/OBJ), depth sort, histograms
CMakeLists.txt       Emscripten build
```
