# Vulkan backend (step 5 of the backend plan — see ../README.md)

Design decisions here are derived from THIS codebase (primitives as types,
camera model / SH degree / distortion as value parameters, eval3d vs
non-eval3d sharing one kernel source, single-stream engine, Slang device-math
layer), not from VkSplat. VkSplat (Apache-2.0) is used as a reference for
Vulkan mechanics only; its research-prototype patterns (dedicated allocation
per buffer, busy-wait single fence, per-variant `-D` shader builds,
descriptor-set data plane, compiled-out error handling) are deliberately not
reproduced.

## Verified premises (probed 2026-07 with slang-2026.2.1, SDK 1.3.296;
## Vulkan backend now pins slang-2026.12.0.1 — see "Slang compiler notes")

- The existing `slang/*.slang` device functions compile to SPIR-V unchanged:
  a compute entry point that `#include`s `primitive_3dgs.slang` and calls
  `projection_3dgs<cam, ...>` validates for Vulkan 1.2 requiring only
  `PhysicalStorageBufferAddresses` + `Int64` (+ `Shader`).
- Raw `T*` parameters in an entry point's `uniform Params` struct become a
  push-constant block of 64-bit device addresses at natural C offsets — the
  byte layout of the corresponding C++ struct. This is what lets the Slang
  math (written against raw pointers for the CUDA target) be reused verbatim.
- `[SpecializationConstant]` + `switch` over instantiations of
  `projection_3dgs<let camera_model, ...>` folds at pipeline creation.

## Device baseline

Vulkan 1.2 core with features: `bufferDeviceAddress`, `shaderInt64`,
`timelineSemaphore` (all core-1.2 features; MoltenVK exposes all three on
Apple silicon). Optional, probed per device and reflected as pipeline
variants where needed:

- `VK_EXT_shader_atomic_float` (`shaderBufferFloat32AtomicAdd`) — training
  backward passes (~1,069 float atomicAdd sites). Absent on MoltenVK and
  llvmpipe → CAS-loop emulation variant selected AT RUNTIME (unlike
  VkSplat, where the emulated variants were never built or dispatched).
  Forward/render path does not need it (only float atomicMax on radii,
  which is integer-monotonic on non-negative floats → emulate with
  u32 atomicMax over the float bit pattern; exact, not a CAS loop).
- Subgroup size: never assumed 32. Kernels use `WaveActive*` /
  `WavePrefixSum` and size-agnostic reductions; AMD wave64 and Intel
  variable-width are first-class.

## Memory model

`backend::device_malloc` = one `VkBuffer` (+ dedicated `VkDeviceMemory`,
`STORAGE | TRANSFER_SRC | TRANSFER_DST | SHADER_DEVICE_ADDRESS`) per
allocation; returns the buffer device address as `void*`. Rationale:
`DevicePool` (Tensor.h) already pools with high-water-mark semantics above
this seam — allocation count is bounded by pool slots (≪
`maxMemoryAllocationCount` ≈ 4096), so a suballocator below the pool would
duplicate what the pool does. Revisit only if allocation count observably
grows.

- A global address-range map (base address → {VkBuffer, size, memory})
  backs `is_device_pointer` (range check), memcpy routing for
  `MemcpyKind::Auto`, and address→(buffer, offset) resolution for
  `vkCmdCopyBuffer`/`vkCmdFillBuffer`. Pointer arithmetic on device
  pointers stays valid: BDA addresses are linear within a buffer.
- `host_malloc_pinned` = dedicated host-visible|coherent buffer,
  persistently mapped; returns the MAPPED pointer (host-dereferenceable,
  like `cudaMallocHost`). Registered in a host-address map so async D2H/H2D
  copies to/from pinned memory record a direct `vkCmdCopyBuffer` — no
  staging hop, mirroring CUDA's pinned fast path.
- Pageable-host copies go through a persistently-mapped staging ring with
  timeline-based reclamation; sync copies wait, async-to-pageable degrades
  to sync (same as CUDA).
- `memset` = `vkCmdFillBuffer` with the byte replicated to a 32-bit word.
  Allocation sizes are rounded up to 4 so whole-tensor fills (the only
  pattern in the engine) never need a byte tail; non-4-aligned fills are
  rejected loudly rather than silently mis-set.
- `device_free` honors the documented contract (all in-flight work
  completes first) by flushing every stream and waiting the queue timeline
  before destroying.

## Streams, submission, synchronization

One compute queue. `backend::Stream` = pointer to a per-stream command-batch
context (`kDefaultStream`/nullptr = lazily created default context, matching
the engine's stream-0-only launch macros; the viewer's second stream maps to
a second context — its high-priority preemption nuance is intentionally
dropped for now and noted as a limitation).

- Each context: `VkCommandPool` + a small ring of command buffers.
  Dispatches/copies record into the open command buffer; a flush ends it
  and submits, signaling one shared queue-level timeline semaphore with a
  monotonically increasing value. Submission order on the single queue
  serializes streams — a strictly stronger ordering than CUDA streams,
  never weaker.
- Flush points: `stream_synchronize` / `device_synchronize` (wait timeline
  — `vkWaitSemaphores`, never fence busy-wait), sync memcpy/memset,
  `event_record`, staging-ring reclamation.
- Barriers: a conservative global `COMPUTE|TRANSFER → COMPUTE|TRANSFER`
  memory barrier after every dispatch/copy inside a batch reproduces CUDA
  stream ordering exactly. Per-resource narrowing is a later optimization,
  applied with profiling, not by default.
- `Event` = {timeline value at record, optional timestamp-query slot}.
  `event_elapsed_ms` reads the two timestamps × `timestampPeriod`.
- Errors: sticky `VkResult` + message; `last_error()` returns-and-clears
  (cudaGetLastError semantics). Nothing is compiled out: failures always
  set the sticky error, and allocation failures return nullptr.

## Kernel ABI

Per-kernel `struct XxxParams` defined once in the entry-point `.slang` file
and mirrored in the C++ launcher (static_assert'd layout). Members are
device addresses (the by-value pointer structs of the CUDA ABI —
`WorldBuffer`/`ScreenBuffer` `TensorArray<N>`, `DeviceTensor2D`,
`CameraDistortionCoeffsBuffer` — flatten to address + stride/size fields
with identical semantics) plus scalars.

- sizeof(Params) ≤ 128 → pushed directly (spec-guaranteed floor;
  the three local test devices all report 256).
- Larger → the struct is written to a mapped params ring and a single
  device address is pushed. Explicit per kernel; no silent switching.

Variant axes follow the CUDA instantiation structure
(`generate_kernel_instantiation.py`):

- **Primitive** (Vanilla3DGS / MipSplatting / Vanilla3DGUT) is a type
  parameter → separate entry point (buffer layouts differ).
- **Camera model, SH degree, DistortionType, output_median, VALUE_BITS**
  are value parameters → specialization constants on ONE SPIR-V module,
  folded at pipeline creation. Pipelines are created lazily and cached by
  (entry point, spec-constant tuple, device capability set).

## Shaders: build and shipping

Entry points live in `slang/vulkan/*.slang` and `#include` the existing
`slang/*.slang` module files (same trick as the CUDA side's
`#define CudaDeviceExport ForceInline` include). SPIR-V is NEVER committed
(unlike VkSplat): CMake compiles it at build time via
`slang/build_spirv.py` (parallel slangc `-target spirv`, per-entry checksum
cache in the build tree keyed over the transitive #include closure) and
embeds the blobs into one generated TU (byte arrays + name registry) via
`spirulae_splat/embed_spirv.py`, so binaries and the Python module are
self-contained. CMake locates slangc in PATH / `-DSSPLAT_SLANGC=`, checks
`slangc -v` against the pinned `SSPLAT_SLANG_VERSION`, and on miss or
mismatch downloads + extracts the pinned GitHub release into the build tree
(platform-detected archive; verified again after extraction).

## Slang compiler notes (2026-07: pinned 2026.12.0.1 for SPIR-V)

The Vulkan backend pins `SSPLAT_SLANG_VERSION = 2026.12.0.1` (root
CMakeLists); the CUDA emission flow stays on 2026.2.1. The upgrade was
needed because 2026.2.1 could not compile the projection-backward autodiff
(bwd_diff of the full projection chain). Slang autodiff has a history of
bugs (see shader-slang/slang #8777, #11160, both filed from this project) —
three were worked around in the CANONICAL slang sources (all verified to
still compile + emit identical CUDA exports with 2026.2.1):

- `primitive_3dgs.slang`: `_projection_3dgs_vjp` took the projection
  interface as an existential VALUE parameter and called
  `bwd_diff(iface.projection)` — 2026.12 segfaults IR-gen on this shape.
  Rewritten as a generic constraint (`<ProjT : _DiffProjection3DGS>` +
  `bwd_diff(ProjT::projection)`).
- `projection_utils.slang`: the `*_proj_jac` Jacobians used
  `fwd_diff(*_proj<true>)` where the projection writes through an `out`
  parameter — the exact bwd-of-fwd x out-param combination of issue #11160.
  Each camera projection now has a value-returning `*_proj_uv` core (used
  by both the bool wrapper and the Jacobians; identical primal ops in every
  path).
- `pixel_wise.slang`: 2026.12 SPIR-V codegen constant-folds the SCALAR's
  derivative of `scalar * vector` to zero under bwd_diff — SILENT wrong
  code (v_transmittance in blend_background came back all-zero; caught by
  pwtrain_parity). Worked around with an explicit
  `float3(t, t, t) * vector` splat. RULE: when a differentiable function
  needs the gradient of a scalar that multiplies a vector, write the splat
  explicitly, and treat parity coverage of every scalar gradient as
  mandatory.

## Sort / scan (backend/common/SortScan.h implementation)

Implemented in Slang (not the vendored GLSL): classic 3-kernel LSD radix
sort — per-partition histogram (upsweep), spine scan, rank-and-scatter
(downsweep) — 8 bits/pass, `begin_bit`/`end_bit` skipping unused passes,
subgroup-size-agnostic. **Keys are 64-bit or 32-bit** (two compiled variants;
the 64-bit requirement is DECIDED — tile keys are `(tile_id << 32) |
ordered_float_depth` and 32-bit packing overflows on large real-world
scenes). Values int32. Host side ping-pongs the caller's `DoubleBuffer`
exactly like CUB. `inclusive_sum`/`exclusive_sum` = two-level
workgroup-scan (int32/int64, plus a float inclusive_sum for the densify
MCMC sampling cumsum — deterministic per device, but its summation
association differs from CUB's so cross-backend float cumsums agree only
to rounding); `select_flagged` = scan + scatter + count readback.
Parity-tested against CPU references on every available device. Float
keys are sorted by emitting the order-preserving uint32 flip transform
from the producing kernel and sorting raw-unsigned via the int32 overload
(matches CUB's float-key order bit-exactly for identical key bits).

## Phase order (from ../README.md step 5) and status

| Phase | Contents | Status |
|---|---|---|
| 0 | Runtime (this doc's memory/stream model) + pipeline layer + smoke test | DONE 2026-07 |
| 1 | SortScan (64-bit radix, scans, select) + parity tests | DONE 2026-07 |
| 2 | Projection fwd (3 primitives × spec-const cams/SH) + parity vs CUDA | DONE 2026-07 incl. 3DGUT + packed + SH value-quant q8/q16 |
| 3 | Background SH fwd, tile intersect (key gen + offsets over SortScan) | DONE 2026-07 |
| 4 | Rasterization fwd (2D + eval3D shared source), visualizer blit | DONE 2026-07 (+ render-path PixelWise, engine-level end-to-end parity) |
| 5 | Training: losses, backwards (float-atomic variants), optimizer, densify | IN PROGRESS: optimizer/utility kernels, pixel-wise backwards, background-SH backward, rasterization backward, projection backward family (plain + quantgrad + fused-optimizer), densify/MCMC + remaining optimizer variants (adamtr color, fused 3DGS geometry), bilagrid/PPISP color pipeline (all five sampler families, TV/channel-mean, fused TV-Adam/AdaGrad, PPISP image transform + reg losses), multi-scale per-pixel loss stack (per-pixel fused losses, fused SSIM, image pyramid, edge-aware maps, quantile radix-select) DONE 2026-07 (optim_parity, pwtrain_parity, raster_bwd_parity, proj_bwd_parity, projqgrad_parity, fpbo_parity, optimgeo_parity, densify_parity, bilagrid_parity, ppisp_parity, msloss_parity); remaining stubs: dataset warp/conversion (12) |

Phase-0/1 hard-won portability rules (validated on NVIDIA proprietary,
Mesa ANV, llvmpipe — all tests + validation layers clean on all three):

- **Pin the subgroup size.** Compute pipelines are created with
  VK_EXT_subgroup_size_control's required-size + full-subgroups whenever the
  device supports it. Intel/ANV's default is a VARYING SIMD width, which
  silently breaks `tid / WaveGetLaneCount()` subgroup indexing (the radix
  sort produced garbage on Intel until pinned). Kernel workgroup sizes must
  stay multiples of every plausible subgroup size (8..128). Specifically the
  **X dimension** must be the multiple (VUID-02757 under
  REQUIRE_FULL_SUBGROUPS) — a 16x16 workgroup fails even though its total is
  256, so 2D kernels use flat 64/128/256-wide X layouts.
- **No 8-bit types in shaders.** `uint8_t` loads / `uint8_t*` push-constant
  members require shaderInt8 / 8-bit-storage features outside the baseline.
  Byte arrays are read as packed u32 words (`load_flag_byte` pattern);
  device allocations are 16-byte rounded so word reads past `n` stay in
  bounds. Revisit only if SH value-quant kernels measurably need native
  8-bit access (then feature-gate it).
- **Subgroup scratch sizing.** Ballot-ranking scratch is bounded by
  MAX_NSG = workgroup / 8 (llvmpipe's subgroup is 8); at workgroup 128 and
  RADIX 256 this is 17KB shared — inside every device's 32KB floor.

Phase-2 state (projection forward): `kernels/ProjectionFwd.cpp` +
`kernels/ProjectionPackedFwd.cpp` implement the full `ProjectionFwd.cuh` /
`ProjectionPackedFwd.cuh` APIs (fp32 SH) for all three primitives via
`slang/vulkan/projection_fwd.slang` — the same `projection_3dgs<cam,...>` /
`sh{N}_to_color` slang device functions the CUDA kernels use. Spec-constant
axes: camera model, SH degree, antialiased (MIP), eval3d (3DGUT; its
"conic" output lands in the screen scale slot, used by the eval3d
rasterizer for culling). Params exceeding the 128-byte push floor use the
params-ring pattern (`vk::params_alloc` + 8-byte address push); the packed
mask kernel fits and is pushed directly. Float `atomicMax` on radii is
`InterlockedMax` on the u32 bit pattern (exact for non-negative floats).
Packed projection substitutes an int32 0/1 mask (bool would need 8-bit
stores) scanned by `backend::inclusive_sum<int32>` with a 4-byte nnz
readback. Parity: `backend/tests/projection_parity.cpp` builds under BOTH
backends (CUDA `-DSSPLAT_BUILD_BACKEND_TESTS=ON` dumps, Vulkan compares);
30 fused + 6 packed fp32 configs plus 12 fused + 2 packed value-quant
configs, ~7.2M floats, zero tolerance violations on all three local devices
(packed nnz and compaction order match exactly).

**SH value-quant (q8/q16)** is implemented WITHOUT any 8/16-bit device
features: `slang/vulkan/sh_quant.slang` mirrors harmonics.slang's
`_sh_load_q{8,16}` + `sh_coeffs_to_color_q{8,16}` math but reads the packed
byte/short buffer as u32 words (the no-8-bit-types rule above), so quantized
checkpoints render on the plain baseline. `kShValueBits` (0/8/16) is a new
spec-constant axis on the projection entries (declared before `kEval3d`,
which is now spec id 4); bounds / base / stride plumb through the params
structs with the CUDA kernels' exact semantics (stride 0 -> the FPBO
per-splat-block default). Gotcha: the bounds-block division must run in
u32 — Intel/ANV effectively hangs (30+ min CPU-bound in the driver compiler)
lowering the slang-emitted 64-bit integer division; u32 is identical while
the packed buffer stays under 2^32 cells.

Phase-3 + raster-fwd state (full forward render): the whole
projection -> tile-intersect -> rasterize -> background chain now runs on
Vulkan.

- `kernels/BackgroundShFwd.cpp` + `slang/vulkan/background_sh.slang`:
  per-pixel skybox via the same `generate_ray` / `transform_ray_d` /
  `sh{N}_to_color_dir` exports the CUDA kernel wraps. SH degree is a spec
  constant; camera model is a runtime int (as in CUDA). The backward throws
  until the training phase (`TODO` in the file).
- `kernels/IntersectTile.cpp` + `slang/vulkan/intersect_tile.slang`:
  `do_intersect_tile_generic` = count -> `backend::inclusive_sum<int64>` ->
  8-byte n_isects readback -> key write -> `backend::sort_pairs<int64,int32>`
  (begin 0, end 32 + tile bits, exactly the CUB call) -> offset kernel.
  Ellipse-vs-AABB mode is a runtime null-check on proj_conic (the CUDA
  template bool only folds that check). `do_intersect_tile_post` has no
  callers and is not ported.
- `kernels/RasterizeFwd.cpp` + `slang/vulkan/rasterize_fwd.slang`: closely
  follows RasterizationEval3DFwd_kernel.cuh (the CUDA source of both the 2D
  and eval3D kernels). Two entries (2D shared by 3DGS/MIP; 3DGUT eval3d
  reusing `evaluate_alpha_3dgs` / `evaluate_color_3dgs` /
  `compute_3dgut_iscl_rot`), spec axes kCameraModel (eval3d) / kDistType
  (None, D, RGB_D) / kOutputMedian. 64-thread micro-tile workgroups;
  CUDA's `__syncthreads_count(done)` early-out is a monotonic groupshared
  counter (each thread adds exactly once, on its done transition; the
  loop-top barrier publishes it). Params via the ring.
- Parity: `backend/tests/render_parity.cpp` (background images + 8
  full-pipeline configs across primitives / packed / distortion / median,
  ~5.7M floats). NVIDIA: zero violations. Intel / llvmpipe: <= 0.0008%
  violations — near-tie depth keys sort differently when the projected
  depths differ in the last ulp across compilers, reordering the blend for
  a handful of pixels; the violation-fraction cap absorbs this.
- build_spirv.py's entry scanner now tolerates `//` comments between the
  [shader]/[numthreads] attributes and `void`.

Phase-4 completion + engine-level end-to-end (2026-07): the REAL engine
render path (`forward_3dgs` -> background blend -> color space ->
depth-normal -> `engine_blit_view`) runs on Vulkan, verified against CUDA at
the engine level.

- `kernels/PixelWiseRender.cpp` + `slang/vulkan/pixel_wise_render.slang`:
  the render-path PixelWise subset — `blend_background_forward`,
  `blend_background_noise_forward` (hash_uint3 ported bit-exact),
  `rgb_to_srgb_forward`, `depth_to_normal_forward{,_tv}` (16x16 shared-mem
  apron tile as a flat 256-thread workgroup per the X-dim subgroup rule).
- `kernels/Visualizer.cpp` + `slang/vulkan/visualizer.slang`: full
  visualizer — frustum geometry, LBVH build (morton keys over
  `backend::sort_pairs<int64,int32>` bits 0..63, Karras internal nodes,
  bottom-up AABBs), depth min/max, the annotate/colormap/overlay blit, the
  thumbnail capture, and both entry points (`engine_blit_view` cached path +
  legacy `blit_train_cameras_tensor`). Portability substitutions: float
  atomicMin/Max as InterlockedMin/Max over an ORDER-PRESERVING u32 mapping
  (exact; tree AABBs unmapped by a small pass, min_max decoded by the
  consumer), the 8x4 warp-ballot gray fit as a groupshared serial reduction
  (CUDA assumes warp==32), all uint8 I/O as u32 words with a separate RGB8
  pack kernel writing the [H,W,3] byte image in whole words (tail lands in
  the allocation's 4-byte rounding). Viewer work runs on the default stream
  (the high-priority viewer stream nuance is dropped; single-queue
  submission order already serializes). `SSPLAT_VK_DEBUG_SYNC=1` (legacy alias
  `SSPLAT_VIS_DEBUG_SYNC`) syncs + reports after EVERY kernel dispatch —
  the standard way to bisect device-lost/segfault failures to a kernel.
- **UB-loop portability rule**: Visualizer.cu's Karras split search
  (`for (tf = 2; (t = (l+tf-1)/tf) >= 1; tf <<= 1)`) terminates only through
  signed-overflow / divide-by-zero UB. NVIDIA and llvmpipe happen to exit;
  Intel spins forever -> VK_ERROR_DEVICE_LOST. The port uses the bounded
  ceil-halving form (identical probe sequence and converged split). Grep any
  future kernel port for loops whose exit depends on overflow.
- `EngineBackground.cu` became portable `EngineBackground.cpp` (its one raw
  kernel replaced by the existing `float_add_into` launcher, cudaMemcpy* ->
  backend::). It now compiles into BOTH backends unchanged.
- **Speculated-load portability rule (llvmpipe segfault)**: llvmpipe's LLVM
  JIT can SPECULATE a load that sits behind a runtime branch — a
  constant-offset load guarded by `if (p.ptr != nullptr)` executed with a
  null device address segfaulted the process (densify_update_weight's
  accum_weight_scalar; discrete GPUs shrug off null BDA reads, the CPU JIT
  does not). Rules, now applied across every kernel TU: a params pointer
  field a shader may load through is NEVER null — launchers substitute
  `vkk::or_fallback()` (a 4 KB zeroed allocation in KernelCommon.h);
  "null means zeros" inputs (dist_coeffs and friends) just read the
  fallback; "null selects a mode" moved to specialization constants
  (intersect_tile kEllipse/kHasXy/kPacked, rasterize_fwd kRasterPacked —
  spec-folding removes the dead load entirely) or explicit flag fields
  (blit has_lss/has_tri, thumbnails has_alpha, optimizer has_steps/
  grad_quant, densify use_opacs/use_w2). Guarded STORES are fine (never
  speculated).
- **Shared launcher helpers** live in `kernels/KernelCommon.h` (namespace
  `vkk`): fold_1d 65535-grid folding, dispatch/dispatch_flat/dispatch_ring
  (throwing, debug-sync hook), or_fallback, resolve_sh_quant. New kernel TUs
  should not hand-roll these.
- **Optimizer port (phase 5, first slice)**: `kernels/Optimizer.cpp` +
  `slang/vulkan/optimizer.slang` implement fused_adam_step (fp32),
  fused_adam_step_quantized (4/8-bit state), fused_adam_step_quantized_value
  (state x value x optional gradq int8 grad), fused_adagrad_step,
  float_add_into, increment_int32_inplace; `kernels/Densify.cpp` +
  `densify.slang` implement densify_update_weight. The quantized-state
  codecs live in `slang/vulkan/optim_quant.slang` (QuantizedAdamState /
  QuantizedTensor / gradq::Codec<8> mirrors): packed cells are read as u32
  words and written by cooperative word assembly in groupshared (a
  256-cell quant block is always word-aligned; CUDA leaves final-word tail
  bytes untouched, the port zeroes them — nothing reads them). CUDA's
  warp-shuffle min/max block reductions became a deterministic groupshared
  serial reduction; log1pf/expm1f are emulated with the compensated forms
  (parity: codes within +-1 step). Kernels with barriers mask work with
  `inside`/`blk_ok` instead of early returns so every thread reaches every
  barrier. `optim_parity` (dump/compare like the others) covers all of it:
  158k floats + 148k codes, 0 violations on NVIDIA/Intel/llvmpipe.
  densify_update_weight's Median mode draws contraction-sensitive hash14
  randomness — parity compares only its deterministic count channel.
- **Pixel-wise training + background-SH backward (phase 5, second slice)**:
  `kernels/PixelWiseTrain.cpp` + `pixel_wise_train.slang` implement the
  blend/noise/srgb backwards, overexposure_grad_add, depth_to_normal_backward
  (16x16 tile apron + neighbor scatter), linear_depth_to_ray_depth_inplace,
  and color_shift_reg_step (ColorShiftReg.cu; shared-mem float atomics ->
  groupshared serial reduce + one portable atomic add per channel/block).
  `atomic_float.slang` provides the portable float atomic add: a
  compare-exchange loop over the u32 bit image (no VK_EXT_shader_atomic_float
  dependency; same nondeterministic ordering as CUDA atomicAdd).
  render_background_sh_backward replaces CUDA's warp==32 block-atomic SH
  reduction (`sh_block_atomic_add_f3_512`) with an atomics-free design: a
  fixed grid of 16384 grid-striding threads each accumulates a private
  float3 scratch slice through the plain `sh{N}_to_color_dir_vjp_inplace`
  exports, then a single-workgroup deterministic reduce folds slices into
  v_sh_coeffs. `pwtrain_parity` covers all of it (incl. bg SH fwd+bwd at
  degrees 1/3): 221k floats, <=0.003% violations x 3 devices (the outliers
  are fisheye-boundary pixels of the undistort iteration).
- **Rasterization backward (phase 5, third slice)**: `kernels/RasterizeBwd.cpp`
  + `rasterize_bwd.slang` (both primitives; tile constants + the ellipse/box
  cull now shared with the forward via `raster_common.slang`). The CUDA
  kernel runs ONE WARP per 8x8 sub-tile with a back-to-front systolic sweep
  (thread k undoes pixel (t-k) for survivor splat k, __syncwarp per step)
  and __ballot_sync survivor compaction; the port keeps the exact schedule
  with a 32-thread workgroup but replaces every warp-width assumption:
  __syncwarp -> full workgroup barriers (the sweep body is an inlined
  function so its early-outs are returns, never skipping the caller's
  barrier), ballot compaction -> WavePrefixCountBits + per-wave totals in
  groupshared, float atomicAdd/atomicMax -> atomic_add_f32 CAS loop /
  atomic_max_f32 signed-int-max trick (atomic_float.slang) with CUDA's
  isfinite filter. Six spec-constant axes (camera, dist, median,
  accum-weight, viewmat-grad, packed). The 2D fragment vjp math is ported
  inline; 3DGUT goes through the evaluate_{alpha,color}_3dgs_vjp +
  compute_3dgut_iscl_rot_vjp slang exports. `raster_bwd_parity` builds real
  forward state through the ported projection/intersect/raster-fwd chain,
  then compares 285k gradient floats over 5 configs: NVIDIA 1 violation
  (0.0004%); Intel/llvmpipe ~0.05% — verified to be inherited upstream
  near-tie blend-order flips, not a port defect, by comparing those devices
  against an NVIDIA-Vulkan dump of the SAME shader (identical diffs).
  Note: Primitive.cuh zeros_pool had a leftover raw cudaMemset (host code)
  -> backend::memset_sync (both backends re-verified).
- **Projection backward family (phase 5, fourth slice)**: three kernels over
  shared machinery in `proj_bwd_common.slang` (splat/camera loaders, the
  per-intersection projection VJP wrappers over the primitive_3dgs.slang
  autodiff exports, SH color VJP dispatch, CAS scatter helpers) +
  `sh_vjp.slang` (register-array mirror of harmonics.slang's
  `_sh_coeffs_to_color_vjp` float3 form — SPIR-V physical pointers cannot
  address locals and the CUDA `_vjp_atomic` variant needs float
  InterlockedAdd, so the VJP accumulates into an `inout float3[24]` and
  callers either keep it in registers or scatter with atomic_add_f32; source
  coefficients decode once via fp32 or the word-based q8/q16 codec).
  - `kernels/ProjectionBwd.cpp` + `projection_bwd.slang`
    (projection_{3dgs,mip,3dgut}_backward): one thread per intersection,
    spec axes camera/degree/antialiased/value-bits/packed/viewmat-grad;
    world grads + optional v_viewmats accumulate with CAS adds (CUDA's
    labeled-partition warp reduce is a perf shortcut with the same
    order-nondeterministic sums). `proj_bwd_parity`: 1.24M gradient floats
    over 7 configs (3 primitives, 4 cameras, packed, q8/q16, viewmat grad).
  - `kernels/ProjectionBwdQuantGrad.cpp` + `projection_qgrad.slang`
    (projection_*_backward_quantgrad): splat-parallel camera loop in
    registers, then decode -> add -> block-reduce bounds -> encode onto the
    signed-symmetric gradq codec (16-bit non-SH / 8-bit SH, added to
    optim_quant.slang together with `vkq_store_u16/u8` — And+Or masked
    sub-word atomic stores for per-splat cell runs that need not be
    word-aligned; disjoint masks commute under any interleaving). Packed
    preprocessing (identity-perm radix sort by gaussian id via
    backend::sort_pairs + [N+1] camera ranges) lives in
    `vkk::build_packed_camera_ranges`, shared with FPBO. Active attributes
    are a spec bitmask (kGqMask) so inactive loads are dead at pipeline
    creation. `projqgrad_parity`: two sub-batch calls (the second decodes
    real state), 62k floats + 845k codes, +-1 quantum.
  - `kernels/FusedProjectionBwdOptim.cpp` + `fpbo.slang`
    (fused_projection_bwd_optimizer_*): the 1080-line CUDA kernel —
    camera-loop gradients + per-splat regularizers (per_splat_losses.slang)
    + in-place Adam on all six attributes with optional 16-bit qadam non-SH
    state, 8-bit qadam SH state and the 16-bit SH value codec (LEVEL axis).
    11 spec axes; block reduces via the deterministic oq_block_reduce
    helpers. `fpbo_parity`: two steps from the zero quantized state, 1.81M
    floats + 915k codes over 5 configs. Two test-design caveats mirrored
    from CUDA semantics: N must be a multiple of the 256 block (the CUDA
    kernel reads tail-thread aabb/screen entries it later discards — fine
    against pooled slack, not against exact-size test buffers), and
    scale-agnostic means x non-SH quant is excluded from comparison (mixed
    g1/g2 units make u = g1/sqrt(g2) unbounded for radii==0 splats, so +-1
    quantum decode differences amplify chaotically — inherent to the lossy
    codec, on CUDA too).
  - Porting this slice surfaced a HOST bug worth remembering:
    `backend::vk::SpecList` had a fixed 8-entry array and the 11-spec fpbo
    dispatch overflowed it in the initializer-list constructor (clobbering
    `count` and beyond) — NVIDIA segfaulted inside pipeline creation,
    llvmpipe silently ran with garbage spec values. kMax is now 16 with an
    assert; if a kernel ever needs more axes, raise kMax explicitly.
- **Densify/MCMC + remaining optimizer variants (phase 5, fifth slice)**:
  `kernels/Densify.cpp` implements the full densification API
  (densify_clip_scale, MCMC/revised noise, long-axis-split relocate/add,
  MCMC relocate/add) over `slang/vulkan/densify.slang`, which includes the
  canonical `slang/densify.slang` device code (long_axis_split_3dgs, the
  randn3 noise kernels) and mirrors the Densify.cu-only mcmc_relocation
  math; `kernels/Optimizer.cpp` gains the adamtr trust-region color
  variants (`slang/vulkan/optim_color.slang`) and the fused 3DGS geometry
  step (`slang/vulkan/optim_geometry.slang`, sharing the FPBO helpers
  refactored into optim_quant.slang: oq_state_read/oq_accum_us/
  oq_state_encode/oq_sqrt3/oq_sqrt4). Design notes:
  - **Weighted sampling without replacement**: the efraimidis-spirakis
    float key `w / log(1 - u)` is emitted through the order-preserving
    uint32 flip (`b ^= (b >> 31) ? 0xffffffff : 0x80000000`), so the
    raw-unsigned radix sort reproduces CUB's float-key ascending order for
    identical key bits. The `u` draw hashes (seed, blockIdx, threadIdx)
    with an integer hash (bit-exact both backends); the kernels keep
    CUDA's 256-thread linearized geometry so hashes match.
  - **Float cumsum**: `backend::inclusive_sum<float>` (new
    scan_blocks/spine/add_f32 entries in sort_scan.slang) feeds the MCMC
    binary-search draws. Its summation association differs from CUB's, so
    cross-backend cumsums agree only to rounding — rare boundary flips in
    the sampled indices are expected and the parity tool budgets for them
    (loose channel / code violation caps; each flip rewrites one splat's
    row).
  - **Relocation compaction** keeps CUDA's atomic-counter form, so the
    src→dst PAIRING is scheduling-dependent (true on CUDA alone);
    densify_parity's relocate configs are built with exactly one
    relocating splat, while the add path (same kernel, deterministic
    dst = cur + i, 300-600 concurrent dsts) covers the multi-splat
    behavior including neighbor-word And+Or masked stores in the 8/4-bit
    SH-state zeroing and 16/8-bit value-codec copies.
  - **Spec IDs follow declaration order per module** (kShOptimBits = 0,
    kShValueBits = 1 in densify.slang): the copy-kernel dispatch initially
    passed `SpecList{bits}`, which set kShOptimBits instead and left
    kShValueBits at its default 16 — in the 8-bit value config the copy
    then ran 16-bit halfword stores at 2x the byte offset, corrupting
    neighboring pool buffers. When a dispatch only cares about a later
    spec axis it must still supply the earlier IDs.
  - The adamtr kernels mirror two upstream CUDA quirks bit-for-bit
    (flagged, not fixed): they never write exp_avg/exp_avg_sq back, and
    the SH variant launches num_gs threads over an element-indexed
    [N, K, 3] tensor so only the first num_gs scalars update per call.
  - Parity: `optimgeo_parity` (1.04M floats + 86K codes; 4 geometry
    configs covering scale-agnostic means, per-splat steps, densify score,
    zero-init non-SH 16-bit qadam two-step protocol, 16-bit grad-quant
    decode; 4 adamtr variants) and `densify_parity` (24K tight + 3.24M
    loose floats + 2.8M codes; 10 configs over sh_optim_bits 32/8/4 x
    sh_value_bits 32/16/8 x bounds layouts x non-SH quant x degree
    warmup). All three devices pass; validation clean.
- **Bilagrid/PPISP color pipeline (phase 5, sixth slice)**:
  `kernels/Bilagrid.cpp` implements the full csrc/BilagridBindings.h raw
  launcher family (affine / PPISP / log-linear / depth / normal samplers:
  uniform + patched + per-pixel-coords + packed, forward + backward v1/v2)
  plus TV loss / channel mean and the EngineInternal.h fused
  TV-Adam/AdaGrad, q16 encode, identity init and scatter;
  `kernels/PpispImage.cpp` implements the PixelWise.cuh PPISP image
  transform + regularization APIs and the PpispInit.cu helpers. Device
  work: `slang/vulkan/bilagrid_{common,affine,ppisp,loglinear,depth,
  normal,tv,optim}.slang` + `ppisp_image.slang` (37 new entries). Design
  notes:
  - **BilagridReader** maps to three never-null pointer fields
    (`vkk::or_fallback`) + a `kValueQuant` spec constant; the
    `grid_indices == nullptr -> identity` semantic becomes a
    has_grid_indices flag field (llvmpipe speculation rule). Uniform vs
    patched layouts share one entry via a `kPatched` spec axis.
  - **PPISP math is the canonical slang** (`slang/ppisp.slang`) on both
    backends: the image transform + reg losses call the same
    apply_ppisp* / compute_*_loss* exports the CUDA build gets through
    generated/ppisp.cuh, so ppisp_parity's per-pixel values AND gradients
    are tight-channel. The bilagrid-PPISP samplers reuse
    apply_exposure/apply_color_correction_ppisp with Slang autodiff for
    the backward — analytically equal to the hand-written CUDA chain in
    BilagridPpispMathBwd.cuh but rounded differently, so those gradients
    are loose-channel. The normal family's axis-angle backward is a
    hand-written Rodrigues-form gradient upstream (NOT autodiff of its
    forward), so it is hand-mirrored term by term and stays tight.
  - **backward_v1 grid-grad kernels** keep CUDA's logical 8x8 block
    geometry / mult_x-mult_y tile decomposition exactly (the per-thread
    pixel-range partition determines the writeback structure), but run as
    flat 64-wide workgroups with an in-shader 8x8 decode — the pinned
    subgroup size requires workgroup X to be a multiple of 32 (VUID-02757;
    same reason the TV/channel-mean 4x4x4 blocks and the single-thread
    reg-loss finals are 64-/32-wide with logical decode / a tid-0 guard).
    The warp-shuffle block reduce becomes a groupshared tree
    (`bg_block_add_64`); CUDA's fast/slow writeback three-path split is
    preserved.
  - **Fused TV-Adam/AdaGrad**: one 256-cell codec block per workgroup;
    Adam moments through QuantizedAdamState<4|8>, the AdaGrad accumulator
    through a new QuantizedTensorLog<4|8> port in optim_quant.slang
    (qtlog_*; the CUDA 4-bit odd/even two-phase byte serialization is
    subsumed by the cooperative whole-word store). The cross-block TV
    neighbor reads race with value writeback on BOTH backends (documented
    partial-update semantic); bilagrid_parity therefore runs tv_weight > 0
    only in the fp32-state configs — through the QUANTIZED states'
    log-domain block bounds (log1p(sqrt(g2)/eps) of near-zero cells) the
    tiny race shift rescales entire blocks of codes, and llvmpipe's
    sequential block execution makes the race systematic.
  - A residual cross-backend sensitivity worth knowing: a single 4-bit
    state code on a rounding boundary (FMA contraction differences) can
    shift one cell's next-step value by ~lr * (state range / 15); if that
    cell is a value-quant block's min/max, the whole block's 16-bit codes
    shift by a fixed offset (~300 codes = ~0.007 in value space on
    NVIDIA). Same +-1-flip amplification class as fpbo_parity; the code
    violation caps budget for it.
  - Parity: `bilagrid_parity` (372K tight + 207K loose floats + 20K
    codes; all 30 sampler launchers incl. q16 readers, grid_indices,
    mult>1 block-reduce path, mi batching, packed int64 indices; TV /
    channel mean; 4 fused-optimizer configs x Adam/AdaGrad; encode/init/
    scatter utils) and `ppisp_parity` (63K tight + 0.9K loose; 3 param
    types x cam_indices on/off + reg fwd/bwd with synthetic deterministic
    bwd inputs). All three devices pass; validation clean.
- **Multi-scale per-pixel loss stack (phase 5, seventh slice)**:
  `kernels/PerPixelLoss.cpp` mirrors `compute_multi_scale_per_pixel_losses`
  (PerPixelLoss.cu) end to end — pyramid construction, per-scale grad
  aliasing/upsample accumulation, and the one-iteration-behind AsyncReadout
  loss/SSIM display path (DevicePool + AsyncReadout are backend-portable).
  `kernels/Quantile.cpp` implements `batch_quantile_masked_radix_select`
  for real (was a manual stub; the `<true>` instantiation is used by
  BilagridBindings). Device work in 4 new modules:
  - `multi_scale_loss.slang` includes the CANONICAL
    `slang/per_pixel_losses.slang` (CudaDeviceExport -> ForceInline, like
    ppisp), so the per-pixel loss math and its autodiff backward are the
    same functions CUDA emits -> per-pixel grads land tight. Bilinear
    GT-resolution sampling/scatter re-derived from Interpolation.cuh.
    Optional buffers become in_flags/out_flags bit fields over never-null
    (or_fallback) pointers; bool masks are read via u32 words and the bool
    2x-downsample writes one whole u32 word per thread (And/Or atomics for
    the tail word) since 8-bit storage isn't baseline.
  - `fused_ssim.slang` ports the memory-efficient inplace SSIM
    forward+backward. CUDA's 24x24 tile needs ~68 KB of shared overlays,
    over llvmpipe's 32 KB floor -> restructured to 16x16 logical tiles in
    flat 256-wide workgroups (~30 KB, same buffer1/buffer2 overlay reuse).
    Every per-pixel output is a pure function of the zero-padded input, so
    tile size does not change the map or dL/dimg1 (tight); only the
    display-only scalar's accumulation grouping differs (loose).
  - `canny.slang` (canny rgb/scalar + BT.601 residual luma + Tukey
    biweight) with globally-clamped loads -> per-pixel deterministic, so
    even the RobustEdgeAware densification map compares tight.
  - `quantile.slang`: 4-pass radix select over float bit images; all
    integer histograms/atomics -> bit-exact across backends and devices.
  - Two CUDA-reference quirks surfaced (both display-only, documented in
    msloss_parity.cpp): the SSIM scalar sums over TILE-GRID positions, so
    non-multiple-of-tile images pick up zero-padded out-of-image
    contributions that differ 24-vs-16; and Common.cuh's
    `blockAtomicAdd<576>` reduces its 18 warp partials with power-of-two
    shuffle strides, silently dropping warps 8 and 17 (~4% low on the
    CUDA side; the Vulkan value matches a numpy brute force of the same
    center set exactly). Gradients and loss values are unaffected.
  - Parity: `msloss_parity` (432K tight + 22K loose; 6 configs covering
    all 7 loss-map modes, equal-shape + scaled GT modalities, masks at
    equal/different resolutions, camera_indices, 1-3 pyramid scales, plus
    direct quantile checks with NaN/Inf/zero poisoning). All three
    devices pass (<= 0.018% tight — isolated branch-flip classes:
    normal_loss's 0.5/sqrt(cos_sim_loss) gradient near aligned normals
    under bilinear-resampled GT, and the Pearson-depth chain fed by
    atomic-order-dependent sums); validation clean.
- **Training stubs**: `kernels/TrainingStubs.gen.cpp` (generated by
  `spirulae_splat/generate_vulkan_stubs.py` from a link probe of
  csrc_portable vs the backend lib) + `TrainingStubsManual.cpp`
  (meshing::OccupancyEvaluator methods) provide throwing
  definitions for every still-CUDA-only launch function, so the FULL
  portable engine links against ssplat_backend_vulkan today. Regenerate
  after porting a module (a ported symbol must lose its stub). The
  kernel-level parity tools never pull these objects (static-lib
  granularity), so they stay engine-free.
- Parity: `backend/tests/engine/engine_render_parity.cpp` links the whole
  portable engine (CMake `backend/tests/engine/` group; CUDA branch builds
  it against csrc under SSPLAT_BUILD_BACKEND_TESTS). It drives
  set_data_3dgs / set_camera_params / engine_init_background_{sh,noise} /
  engine_init_color_space / forward_3dgs (4 primitive-camera configs +
  noise-mode config) / depth_to_normal / engine_viewer_init + set_grid +
  engine_blit_view. ~4.1M floats: NVIDIA zero violations; Intel/llvmpipe
  <= 0.0005% (the usual near-tie sort flips). 259K blit bytes: <= 0.014%
  of bytes differ by more than 2 (frustum-line edge pixels at rounding /
  gray-threshold boundaries). Validation clean on all three devices.

Each phase lands with tests that run on all local Vulkan devices (RTX 5070 =
NVIDIA proprietary, Intel iGPU = Mesa ANV, llvmpipe = CPU/no-float-atomics
stand-in for portability gaps).
