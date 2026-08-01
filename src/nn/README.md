# GPU inference (`src/nn/`)

A general-purpose inference layer on Vulkan + Slang: its own Vulkan runtime
(`nn/vk/`), a device tensor, an op set, and host image I/O. **It knows nothing
about any model** — `src/sam/` is the first thing built on it, and if a feature
detector, a matcher or a depth model is added later it builds on this
unchanged. Do not put model-specific constants here.

Vulkan-only, like `src/sfm/`: built by default for `SSPLAT_BACKEND=vulkan`
(`SSPLAT_BUILD_SAM`), absent from a CUDA build, and invisible to `setup.py` and
the pybind module. It carries its own device; converging the repository's three
Vulkan contexts onto one is `docs/notes/sfm-port-plan.md` phase 6.

Runtime knobs: `SSPLAT_NN_LOG=0..3`, `SSPLAT_VK_DEVICE`, `SSPLAT_PROFILE=1`,
`SSPLAT_VK_VALIDATION=1`, `SSPLAT_NN_DEBUG_SYNC=1`.

## The op layer

An op runs when you call it, writing into a caller-provided output — usually an
`Arena` tensor inside an `ArenaScope`. There is no computation graph, and that
is deliberate: a module's forward pass reads as straight-line C++ next to its
PyTorch reference, and the arena's scoping is the whole memory strategy.

### Layout conventions

Shapes are **PyTorch order** (row-major, slowest dim first). ggml's reversed
`ne[]` is undone once, in the weight loader.

| | |
|---|---|
| feature maps | `[H, W, C]` — channel-last. LayerNorm2d becomes LayerNorm; a 1×1 conv becomes a Linear. |
| token sequences | `[N, n_heads * head_dim]` — a fused qkv projection lands exactly where the attention kernel indexes, with no permute. |
| weights | `[out_features, in_features]`, the checkpoint's bytes unchanged. A `[Cout, Cin, kh, kw]` conv kernel is the same memory as the `[Cout, Cin*kh*kw]` matrix, and `Tensor::asMatrix()` is how ops see it. |

## What is fused, and why

Fusion here is not micro-optimization: at these shapes the extra pass costs
more than the arithmetic.

- **GEMM epilogue** — bias, activation and a residual add. SAM 3's forward pass
  is ~200 Linears and most are followed by all three. The residual may alias
  the output (one thread per element), so `x = W·h + b + x` is in place.
- **LayerNorm + residual** — the post-norm blocks are always
  `x = LN(x + sublayer(x))`.
- **Attention** — flash-style, online softmax. Not optional: the ViT's global
  blocks are 5184 × 5184 × 16 heads, 1.7 GiB of scores.
- **Resize + threshold + byte pack** — mask export, so the host never sees a
  float mask.

## Kernel selection notes

- **GEMM** has three kernels, picked in `OpGemm.cpp`:

  | when | kernel | tile |
  |---|---|---|
  | `M ≤ 4` | `gemm_nt_thin` | one workgroup per output element, reduce over K in registers |
  | `M ≥ 128`, `N ≥ 128`, fp16 weight | `gemm_nt_big` | 128×128, `BK`=16, 8×8 accumulators per thread |
  | otherwise | `gemm_nt` | 64×64, `BK`=16, 4×4 accumulators |

  Small `M` is everywhere in this model — MLP heads run a single token, the
  hypernetwork mask heads one each — and the big tile carries everything that
  matters for frame time. It keeps its weight tile in shared in the fp16 the
  checkpoint already holds, which is why it needs an fp16 weight and why an
  fp32 second operand (a matmul against another activation) falls back.

- **Attention** is one 64-query × 32-key tile whatever the head dim: Q and K are
  staged in 16-dim chunks, so shared memory no longer scales with `head_dim`
  and the per-thread work is a 4×2 register tile of scores plus
  `4 × ceil(head_dim/16)` output accumulators. Any head dim up to 256 works; one
  that does not divide 16 (Hiera runs 96 and 72) just leaves lanes idle in the
  last chunk, and the bounds tests fold away for the aligned case because
  `kHeadDim` is a specialization constant.

  When `AttnOpts::arena` is set and the problem produces too few workgroups to
  cover the device, the launcher splits the key range across extra workgroups
  and merges the partial softmaxes in a second pass — see below for why that is
  the lever this kernel responds to.

## What the kernels are actually bound by

Measured on an RTX 3070 Laptop with `test_ops --bench`, which runs both kernels
at the shapes SAM 3 uses. Worth reading before tuning either of them, because
the two are bound by different things and the obvious lever only works on one.

An Ampere SM retires 128 fp32 FMAs and reads 128 bytes of shared per clock, so
a kernel gets **one shared byte per FMA** and no more. `gemm_nt_big` is exactly
where that predicts: 4×4 accumulators spend 2 bytes/FMA and ran at 3.5 TFLOP/s;
8×8 with an fp32 weight tile spends 1.0 and ran at 4.7; 8×8 with the weight
tile kept fp16 spends 0.75 and runs at 5.3–6.4. Each step tracked the ratio.

`flash_attn` does not behave that way. Cutting its shared traffic per FMA by
1.7× (a 64-key tile with 4×4 score registers and fp16 probabilities) measured
*identical* — 2.99 vs 2.93 TFLOP/s on the ViT window shape. It is bound by
latency and occupancy, and the lever that moves it is more workgroups: the same
kernel over 5184 queries runs at 3.2 TFLOP/s with 16 heads to spread across and
2.0 with one head (81 workgroups on 46 SMs).

That prediction held, and **flash-decoding** is now implemented on the strength
of it. One head over a few thousand queries is the memory attention's shape and
nobody else's, so the launcher splits the key range only when the workgroup
count falls below `kTargetBlocks`:

| shape | single-pass | split-K |
|---|---|---|
| 5184 q × 5184 k, 1 head, d256 | 1.64 TFLOP/s | 2.36 |
| 4096 q × 28 672 k, 1 head, d256 (SAM 2 memory cross) | 2.02 | 2.60 |

which is what a tail-efficiency model predicts: 64 workgroups over 46 SMs is
1.4 waves at ~72% occupancy, and splitting eight ways takes it to 11 waves at
96%. Raising the target past 512 measured flat, so the win is the tail and
nothing else — 2.6 TFLOP/s is this kernel's steady-state rate. Going further
means fp16 operands and cooperative matrix, which is outside the Vulkan 1.2 core
baseline.

The combine pass costs one extra kernel (~0.4% of frame time) and a partial
buffer of `n_splits × nq × dim` floats, taken from the caller's arena, which is
why `AttnOpts` carries an `arena` at all.

Register prefetching and an L2-friendlier grid order were both tried on the
GEMM and both measured flat; they are still in `gemm_nt_big` because they cost
nothing, but do not expect them to buy anything either.
- **Convolution** is chunked im2col + GEMM, not a bespoke kernel. A full im2col
  of the seg head's 288×288×256 3×3 conv would be 764 MiB, so the output
  positions are chunked to keep the column buffer near 32 MiB. Depthwise is
  direct; transposed 2×2/stride-2 is a GEMM plus a scatter (the four kernel taps
  are four output-channel groups). One well-tuned matmul therefore carries ~95%
  of the model's FLOPs.

## Portability rules for new kernels

These are load-bearing; see `src/vk/README.md` for the device baseline.

1. Workgroup **X** must be a multiple of every plausible subgroup size. 2-D
   kernels use flat 64/128/256-wide X and decode `(x, y)` in-shader.
2. Shared memory budget is **16 KiB**, the Vulkan floor.
3. Never leave a params pointer field null if a shader may load through it —
   CPU Vulkan implementations speculate loads across branches. Use
   `vk::or_fallback()`.
4. No 8-bit or 64-bit integer types, and no `f16tof32`.
5. Grids cap at 65535 per dimension. Long flat ranges fold across `y`
   (`Stream::fold1D` + `fold_index`); a kernel that keeps out-of-range threads
   alive for a groupshared reduction must gate its block-indexed writes,
   because the last folded row is padding.
6. Reductions use a deterministic groupshared tree, not wave intrinsics, so
   wave32 and wave64 devices agree bit for bit.

## Testing

`tests/test_ops.cpp` checks every kernel against an independent scalar CPU
implementation. A new op lands with a case there — the model layer treats this
API as trustworthy, and that is only true because of that file.
