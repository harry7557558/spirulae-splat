# Testing

Three suites, in descending order of how much you should trust them.

## 1. Native cross-backend parity tests (the important ones)

`csrc/backend/tests/*.cpp` — currently 17 tools covering projection (fwd, bwd,
quant-grad), rasterization bwd, tile intersect, warp, FPBO, optimizer (general
+ geometry), densify, per-pixel train, PPISP, bilagrid, multi-scale loss, plus
`backend/tests/engine/` which drives the *real* engine end to end
(render parity, train parity). `backend/vulkan/tests/` adds 3 Vulkan-only
smoke tests (runtime, pipeline, sort/scan).

**The same source builds under both backends.** Each file only touches
`backend::`, the generated launch declarations, and `Tensor.h`. The workflow
is dump-then-compare:

```bash
# on the CUDA machine
./build_cuda/projection_parity dump ref.bin
# on the target machine / device
./build/projection_parity compare ref.bin
```

Inputs are deterministic, and comparison is tolerance-based — fast-math
exp/sqrt chains legitimately differ across compilers, and borderline-cull
flips change whole rows, so a small allowance for those is built in.

Several tools also take a `*_DUMP_GOT` environment variable
(`FPBO_DUMP_GOT`, `PPISP_DUMP_GOT`, `MSLOSS_DUMP_GOT`, `PWTRAIN_DUMP_GOT`,
`BILAGRID_DUMP_GOT`, `DENSIFY_DUMP_GOT`) to write the *actual* values
alongside the reference, which is how you diff a mismatch numerically instead
of guessing.

### Building them

```bash
# CUDA branch: opt-in
bash build_develop.bash -B build_cuda -DSSPLAT_BACKEND=cuda -DSSPLAT_BUILD_BACKEND_TESTS=ON
# Vulkan branch: built unconditionally
bash build_develop.bash -DSSPLAT_BACKEND=vulkan
```

Each `.cpp` becomes an executable of the same base name in the build dir.

### Cross-machine / cross-vendor runs

The comparison target is often a different machine (e.g. an AMD GPU box), and
often offline. The pattern that works:

1. Transfer a matching `slangc` to the target and point `-DSSPLAT_SLANGC=` at
   it — SPIR-V is compiled at build time and never committed, so the target
   needs a compiler, and the version is pinned.
2. Dump references on the CUDA host.
3. Copy the `.bin` files over and run `compare` on the target.

Keep reference dumps out of git (`parity_refs/` is gitignored).

## 2. GUI / viewer checks

The web viewer can be driven headlessly over the Chrome DevTools Protocol.
Headless defaults to SwiftShader; to exercise a real GPU, run against a real
display (`DISPLAY=:0`).

If you script a Python training run that serves the viewer, pass
**`--no-keep-viewer-alive`** or the process hangs at exit waiting on it.

## 3. Python tests (`tests/`)

```bash
pytest tests/
```

`test_per_pixel.py`, `test_ppisp.py`, `test_rasterization.py`,
`test_render_background_sh.py`, `test_ssim.py`, plus `utils.py`.

**These are unmaintained and may fail.** They date from the PyTorch-first era,
compare CUDA kernels against `_torch_impl.py` reference implementations, and
have not been run recently. Treat a failure here as "needs investigation",
not as a release blocker — but do investigate, since they are the only
remaining check on the `_torch_impl.py` reference math.

There is also `csrc/tests/test_delaunay3d.py` (a Python test living in the C++
tree) and `csrc/tests/delaunay3d_bench.cpp`.

## What to run before calling a change done

| change | gate |
|---|---|
| any kernel | CUDA build + Vulkan build + the relevant parity test on both |
| engine logic | both builds + `engine_render_parity` + `engine_train_step`-level check |
| config field | rerun `generate_cli_config.py`; check `ssplat-train --help` |
| build system | all four modes in [build.md](build.md) |
| Python-facing | a short `spirulae-train` run with `--no-keep-viewer-alive` |
| anything | one short training run per backend on a public scene |

## Profiling

`SSPLAT_PROFILE=1` enables the env-gated per-stage timing breakdown
(H2D / D2H / D2D / memset / device / host). Header-only, works on both
backends — the right first tool when a backend is unexpectedly slow rather
than wrong.
