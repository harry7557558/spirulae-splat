# Vulkan backend (step 5 of the backend plan — see ../README.md)

Design decisions here are derived from THIS codebase (primitives as types,
camera model / SH degree / distortion as value parameters, eval3d vs
non-eval3d sharing one kernel source, single-stream engine, Slang device-math
layer), not from VkSplat. VkSplat (Apache-2.0) is used as a reference for
Vulkan mechanics only; its research-prototype patterns (dedicated allocation
per buffer, busy-wait single fence, per-variant `-D` shader builds,
descriptor-set data plane, compiled-out error handling) are deliberately not
reproduced.

## Verified premises (probed 2026-07 with slang-2026.2.1, SDK 1.3.296)

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

## Sort / scan (backend/common/SortScan.h implementation)

Implemented in Slang (not the vendored GLSL): classic 3-kernel LSD radix
sort — per-partition histogram (upsweep), spine scan, rank-and-scatter
(downsweep) — 8 bits/pass, `begin_bit`/`end_bit` skipping unused passes,
subgroup-size-agnostic. **Keys are 64-bit or 32-bit** (two compiled variants;
the 64-bit requirement is DECIDED — tile keys are `(tile_id << 32) |
ordered_float_depth` and 32-bit packing overflows on large real-world
scenes). Values int32. Host side ping-pongs the caller's `DoubleBuffer`
exactly like CUB. `inclusive_sum`/`exclusive_sum` = two-level
workgroup-scan; `select_flagged` = scan + scatter + count readback.
Parity-tested against CPU references on every available device.

## Phase order (from ../README.md step 5) and status

| Phase | Contents | Status |
|---|---|---|
| 0 | Runtime (this doc's memory/stream model) + pipeline layer + smoke test | DONE 2026-07 |
| 1 | SortScan (64-bit radix, scans, select) + parity tests | DONE 2026-07 |
| 2 | Projection fwd (3 primitives × spec-const cams/SH) + parity vs CUDA | DONE 2026-07 incl. 3DGUT + packed (SH value-quant pending) |
| 3 | Background SH fwd, tile intersect (key gen + offsets over SortScan) | — |
| 4 | Rasterization fwd (2D + eval3D shared source), visualizer blit | — |
| 5 | Training: losses, backwards (float-atomic variants), optimizer, densify | — |

Phase-0/1 hard-won portability rules (validated on NVIDIA proprietary,
Mesa ANV, llvmpipe — all tests + validation layers clean on all three):

- **Pin the subgroup size.** Compute pipelines are created with
  VK_EXT_subgroup_size_control's required-size + full-subgroups whenever the
  device supports it. Intel/ANV's default is a VARYING SIMD width, which
  silently breaks `tid / WaveGetLaneCount()` subgroup indexing (the radix
  sort produced garbage on Intel until pinned). Kernel workgroup sizes must
  stay multiples of every plausible subgroup size (8..128).
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
30 fused + 6 packed configs, ~5.2M floats, zero tolerance violations on all
three local devices (packed nnz and compaction order match exactly).
Pending in this phase: SH value-quant — the q8/q16 harmonics exports need
per-entry blobs gated on shaderInt8/16-bit-storage features (per-entry
compilation keeps the capabilities out of the baseline blobs), plus parity
inputs in the engine's real codec layout.

Each phase lands with tests that run on all local Vulkan devices (RTX 5070 =
NVIDIA proprietary, Intel iGPU = Mesa ANV, llvmpipe = CPU/no-float-atomics
stand-in for portability gaps).
