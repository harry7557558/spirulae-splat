# Spirulae Standalone Viewer

A dependency-free, client-side **WebGL2** viewer for 3D Gaussian Splats,
meshes, and **MVS datasets** (COLMAP / Nerfstudio / Metashape). It runs
entirely in the browser (no server) and can be hosted as a static site (e.g.
GitHub Pages). Performance-critical parsing, depth sorting, and statistics run
in **C++ compiled to WebAssembly** (built with CMake + Emscripten); rendering
is WebGL2.

The viewer's own runtime code is **independent** from the rest of
`spirulae_splat` — nothing here imports from the training code at runtime. The
one deliberate exception is at build time: the WASM module compiles the
trainer's dataset parsers (`spirulae_splat/splat/cuda/csrc/app/*Parser.cpp`)
**in place** — referenced by relative path, not copied — so COLMAP/Nerfstudio/
Metashape parsing has exactly one implementation in the repo. Those files are
plain C++17 with no CUDA dependency (see `csrc/CameraModel.h`).

## Features

- **3D Gaussian Splatting** (`.ply`, INRIA/spirulae layout, including very
  large files — binary PLY is parsed **streaming** through a small chunk
  buffer, non-position attributes are stored as half floats, rest SH
  coefficients are **8-bit quantized** (Gaussian-wise scale, `RGB8_SNORM`
  textures — half the VRAM and heap of f16, visually lossless), and the SH
  texture / mesh index buffers are split into chunks below per-resource GPU
  limits; tested with a 6.4 GB, 39M-splat SH2 scan (~3 GB VRAM). Splats are
  **Morton-reordered** at load so the depth-sorted draw order touches
  attribute textures cache-coherently (several-fold frame-rate gain on
  multi-10M-splat models); the async sort worker keeps only an xyz copy.
  Deep-linked `?model=` URLs stream straight into the parser (no Blob
  buffering, so multi-GB hosted models load). If the GPU runs out of memory
  on the SH textures, the viewer drops one SH degree at a time and retries
  instead of failing.
  - Primitive select: **3DGS**, **Mip** (antialiased), **3DGUT** (unscented
    transform projection; fragments evaluate the 3D Gaussian along per-pixel
    rays, "eval3d"). 3DGS/Mip use analytic projection Jacobians with the
    conventional out-of-bounds Jacobian clipping.
  - **Color gamut** (Rec.709 / DCI-P3 / Rec.2020 / AdobeRGB / ACEScg /
    ACES2065-1) and a **linear-color** toggle, matching the training pipeline's
    `rgb_to_srgb` conventions.
  - Spherical harmonics up to degree 4, exposure.
  - Depth sorting matches the training code's `get_sorting_depth`: planar for
    perspective/orthographic, distance for equirectangular, and a smooth
    |z|/radial blend for fisheye — content behind the camera composites
    correctly in >180° views. Sorting runs **asynchronously in a Web Worker**
    that hosts its own instance of the WASM module (native counting sort over
    a positions copy), so looking around never blocks the render loop; results
    apply latest-wins when ready.
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
- A rasterized, depth-tested **axes + grid** overlay in the model's native
  frame (power-of-10 cells that adapt to the zoom level; the line patch
  follows the orbit target while staying on the global lattice) and a
  configurable background.
- **Statistics** (splat count / vertices / edges / faces) and on-demand
  **parameter histograms** (opacity, scale, effective rank, RGB, anisotropy for
  splats; edge length, triangle area, coordinates for meshes), computed in WASM
  and cached.
- **MVS datasets** — drag & drop a **COLMAP** reconstruction
  (`cameras`/`images`/`points3D`, binary `.bin` or text `.txt`), a
  **Nerfstudio** `transforms.json` (+ point cloud `.ply`), or a **Metashape**
  camera export (`.xml` + `.ply`, optional `.psx`) — as a folder or as
  individual files:
  - Renders the **seed point cloud** and **camera frustums** with the true
    per-camera projection: the image border is discretized and unprojected
    through the camera model (pinhole / fisheye / equisolid / equirectangular)
    **including OpenCV distortion** (Newton undistort, `k1–k4 p1 p2 s1 s2 b1
    b2`), so a fisheye camera's frustum visibly bulges. Wide cameras get an
    image-aligned wire dome (fisheye) or a lat/long wire globe
    (equirectangular) instead of a lone border ring. Point size and frustum
    size are adjustable; either layer can be hidden.
  - A compact info panel: source format, image / point counts, camera groups,
    per-model image counts.
  - **Component picker** when the dataset has several reconstructions
    (`sparse/0`, `sparse/1`, multiple Metashape `<component>` chunks, …).
  - **Hover** anywhere on a camera frustum (ray-picked, not just its apex) to
    see its intrinsics (model, fx/fy/cx/cy, distortion); **double-click** a
    camera to view the scene from it (average focal, cx/cy/distortion
    omitted); double-click the point cloud to recenter. Picking is
    depth-ordered: whatever is visually in front at the cursor — frustum or
    point cloud — wins.
  - Degrades gracefully: dropping only `sparse/`, only a `transforms.json`,
    only a Metashape `.xml`, or only a point-cloud `.ply` shows whatever is
    available (images are never required, or read at all — only the metadata
    files are parsed). Load failures and partial loads pop a visible error
    toast (no need to open the console).
- **Double-click** the viewport to center the view on the point under the
  cursor (MeshLab-style: the point becomes the orbit pivot and slides onto the
  optical axis). Splats use a one-pixel GPU depth pass; meshes use a WASM
  raycast; datasets pick the nearest point along the view ray. Works with every
  camera model.
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

Drag a model or a dataset folder onto the canvas, or click to browse (the
picker takes multiple files; folder drops walk the directory tree). You can
also deep-link a hosted model: `index.html?model=<url>` (external
`.bin`/`.mtl`/image siblings are fetched automatically).

## Test

`test/run.sh` generates synthetic COLMAP / Nerfstudio / Metashape datasets and
drives the viewer end-to-end in headless Chrome (needs `google-chrome` and
`node` ≥ 20):

```bash
test/run.sh            # all cases
test/run.sh metashape  # one case; see test/test_ds.html for the list
```

## Layout

```
index.html          UI + panel
css/style.css        theme (adapted from the training viewer)
js/
  main.js            app wiring, input, render loop, sorting
  renderer.js        WebGL2 renderer (splat HDR pass, mesh, dataset, grid lines, tonemap)
  shaders.js         GLSL (splat / tonemap / mesh / points / frustum+grid lines)
  camera.js          camera + navigation modes
  dataset.js         dataset load orchestration + frustum geometry (undistort)
  sortworker.js      async depth-sort worker (own WASM instance)
  wasm.js            WASM bridge + streaming loader + GLTF/GLB (JS) + MEMFS mount
  colors.js          gamut matrices
  linalg.js          vec/quat/mat helpers
  histogram.js       canvas bar chart
  ssv_wasm.{js,wasm} built WASM module
src/viewer.cpp       parsers (PLY/OBJ), depth sort, histograms
src/dataset_bridge.cpp  dataset C ABI over MEMFS (drives the trainer's parsers)
test/                end-to-end harness (headless Chrome) + data generators
CMakeLists.txt       Emscripten build (also compiles ../spirulae_splat/.../app parsers)
```

Frustum lines use the same seam handling as meshes: fragments whose
interpolated camera-space position re-projects far from the rasterized
position (a segment wrapped across the equirect ±180° seam or the fisheye
backward point) are discarded, and fisheye display models clip to the image
circle. Double-clicking a dataset camera switches the display projection to
that camera's model (pinhole → perspective, fisheye → equidistant, equisolid,
equirectangular) with the matching field of view.
