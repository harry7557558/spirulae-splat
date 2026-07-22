# Backends

Two interchangeable compute backends behind one seam. **Both must work.**

The authoritative documents live next to the code and are more detailed than
this page — this is the orientation, they are the reference:

- **`src/backend/README.md`** — the seam's design:
  what `api/` is, why the per-kernel `.cuh` headers are authoritative and the
  module headers are forwarders, and how `BackendRuntime.h` abstracts
  alloc/memcpy/memset/events.
- **`src/backend/vulkan/README.md`** — the Vulkan
  backend in depth: verified premises, device baseline and per-capability
  pipeline variants, memory model, atomics strategy, sort/scan, Slang compiler
  notes, and **the kernel coverage table**. Read it before touching anything
  under `backend/vulkan/`.

## The seam

```
Engine*.cpp  ──►  backend/api/*.h              (launch declarations)
                  backend/api/BackendRuntime.h (alloc / memcpy / memset / events)
                        │
        ┌───────────────┴───────────────┐
  CUDA: src/kernels/**/*.cu       Vulkan: src/backend/vulkan/
        + src/instantiations/           runtime + kernels/ + shaders/
        │                                │
        └────► backend/common/SortScan.h ◄────┘
               CUB (CUDA) / onesweep radix + prefix sum (Vulkan)
```

The engine layer is already torch-free *and* CUDA-free: it calls kernels only
through launch functions exported via `/*[AutoHeaderGeneratorExport]*/`
markers. The `api/` directory makes that implicit boundary explicit.

`backend/api/BackendTypes.h` provides portable CUDA-layout vector types and
empty `__device__`/`__host__` macros, so the whole declaration surface parses
under `-DSSPLAT_BACKEND_VULKAN` with no CUDA toolkit installed.

## Shared device math

`src/shaders/*.slang` is compiled **twice**: to CUDA headers (`src/generated/*.cuh`)
and to SPIR-V. The projection, SH, primitive and PPISP math therefore has one
source. `src/backend/vulkan/shaders/*.slang` adds the Vulkan-only compute entry
points and the compatibility layers (`int64_compat`, `int8_compat`, `atomic_float`).

Don't fork device math per backend. If a backend needs different code, it
belongs in a capability variant (below), not a copy.

## Vulkan capability variants

Rather than requiring optional features, each entry point is compiled into
several SPIR-V blobs and the pipeline layer picks per device at module load —
no in-shader branching:

| suffix | feature | fallback when absent |
|---|---|---|
| `.atomicadd` | `VK_EXT_shader_atomic_float` | CAS-loop emulation (MoltenVK, llvmpipe) |
| `.noint64` | `shaderInt64` | `uint2` word-pair emulation for keys, morton codes, indexing |
| `.int8` | `shaderInt8` + `storageBuffer8BitAccess` | packed u32 word access |

Baseline is Vulkan 1.2 core with `bufferDeviceAddress` + `timelineSemaphore`.
Subgroup size is **never** assumed to be 32 — AMD wave64 and Intel
variable-width are first-class. `SSPLAT_VK_NATIVE_ATOMICS=0`,
`SSPLAT_VK_NATIVE_INT64=0`, `SSPLAT_VK_NATIVE_INT8=0` force the fallback blobs
for A/B testing.

The Vulkan README documents hard-won specifics that are easy to regress —
e.g. the emulated i64 scan accumulator must stay a `uint2` rather than a
two-field struct, because one Intel Windows driver segfaults in
`vkCreateComputePipelines` on the struct form.

## Adding or changing a kernel

Three deliverables, always:

1. **CUDA**: `src/kernels/<family>/<Kernel>.cu` (+ `_kernel.cuh` if the device body is shared
   between fwd/bwd or eval3d/non-eval3d), with
   `/*[AutoHeaderGeneratorExport]*/` on the launcher. Rerun
   `generate_headers.py`, and `generate_kernel_instantiation.py` if it is
   templated over primitive / camera model / SH degree.
2. **Vulkan**: the Slang entry point in `backend/vulkan/shaders/*.slang` and its
   launcher in `backend/vulkan/kernels/*.cpp`.
3. **Parity test**: a tool in `backend/tests/` that builds under both backends
   and does `dump` / `compare`. See [testing.md](testing.md).

If step 2 isn't done yet, `generate_vulkan_stubs.py` emits a **throwing** stub
so the portable engine still links. That is a TODO marker, not a completed
port — and rerunning the generator after a port lands is required, or the
symbol ends up defined twice.

If a new subsystem appears, add its module header to `generate_backend_api.py`
and add the kernel family to `generate_headers.py`'s `header_names` list.

## Performance notes

Backend performance is not uniform, and the differences are large enough to
change design decisions:

- Write-coalescing, not atomic throughput, has repeatedly been the dominant
  cost in the fused projection-backward/optimizer path on Vulkan. Staged
  coalesced writeback plus a tree reduction fixed it.
- Vulkan driver choice matters (`amdvlk` vs RADV differ by an order of
  magnitude on some kernels). When comparing, state the driver.

Use `SSPLAT_PROFILE=1` for the per-stage breakdown before optimizing anything.
