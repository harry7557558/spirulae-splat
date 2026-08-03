# ALIKED + LightGlue: the learned frontend

The plan for `src/sfm/README.md`'s unstarted item 15, and the record of what
was decided and why. Companion to [sfm-port-plan.md](sfm-port-plan.md) and
[segmentation-port.md](segmentation-port.md).

Goal: `--features aliked-n16rot|aliked-n32` and `--matcher lightglue` as
alternatives to GPU SIFT and the brute-force matcher, reaching the same
`FeatureSet` / `IFeatureMatcher` contracts, with **no new dependency** — no
onnxruntime, no PyTorch, no self-hosted weights.

## What is actually being ported

Nothing from COLMAP's `src/colmap/feature/aliked.cc`. That file is 341 lines of
ONNX plumbing around a black box, and `bruteforce-matcher.onnx` is a **5 KB**
graph doing cosine + ratio + cross-check, which `sfm/feature/Matcher.h` already
does. What gets ported is:

1. **Two networks, reimplemented on `src/nn/`** — the layer that exists for
   exactly this ("A learned feature detector for SfM … goes on top of it
   unchanged"). This is the work.
2. **COLMAP's pre/post-processing conventions, literally.** Image
   normalization, the pad-to-/32 replicate padder, and the keypoint coordinate
   un-normalization. These are five lines each and every one of them is a place
   to lose half a pixel silently.
3. **COLMAP's option semantics**, which are not the same shape as ours:
   `min_cossim = 0.85`, `max_ratio = 1.0` — i.e. the ratio test is *off* and an
   absolute similarity threshold is the only filter. Our `reduce()` has no such
   threshold today.

## Weights: the ONNX file is the checkpoint

**We do not host or convert anything.** `AlikedWeights` fetches the same
artifact COLMAP does, from the same URL, verifies the same SHA-256, and parses
it in process:

```
https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n16rot.onnx
  39c423d0a6f03d39ec89d3d1d61853765c2fb6a8b8381376c703e5758778a547   3.0 MB
https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n32.onnx
  a077728a02d2de1a775c66df6de8cfeb7c6b51ca57572c64c680131c988c8b3c   4.2 MB
https://github.com/colmap/colmap/releases/download/3.13.0/aliked-lightglue.onnx
  b9a5de7204648b18a8cf5dcac819f9d30de1a5961ef03756803c8b86c2dceb8d
```

(The URI triples live in COLMAP's `src/colmap/feature/resources.h`; ours mirror
them.) Consequences, all good:

- Bit-identical weights to what COLMAP runs, so a parity gate against
  `colmap feature_extractor --FeatureExtraction.type ALIKED_N16ROT` means
  something.
- No converter to keep in sync, nothing to publish, nothing to commit. The
  cache directory and the download path are `src/app/gui/ModelCache.cpp`'s
  business, unchanged.
- Licensing is the easy case: ALIKED is BSD-3-Clause and LightGlue is
  Apache-2.0, so this is a `needs_tick = false` family, unlike SAM 3.

### Parsing ONNX without a dependency

An ONNX file is a protobuf. We need **initializers only** — the graph structure
is hard-coded in our forward pass — so the reader is a ~200-line varint walk of
three nested messages and no schema:

```
ModelProto  field 7  -> GraphProto
GraphProto  field 5  -> repeated TensorProto     (the initializers)
TensorProto field 1 dims (packed or repeated), 2 data_type,
                  8 name, 9 raw_data, 4 float_data
```

Everything else in the file is skipped by wire type. Verified against both
checkpoints: 2519 nodes, **69 initializers, all f32**, 2.71 MB of a 3.00 MB
file. `src/aliked/model/Onnx.cpp` is that walk and nothing more.

### What the checkpoints say

The two variants have **identical graphs** — 2519 nodes, 69 initializers each —
and differ only in `M`, the number of SDDH sample positions, which is
`desc_head.agg_weights.shape[0]`. So: one code path, `M` read from the file,
`n16rot` vs `n32` is a URL. `n16rot` is architecturally `n16`, trained
rotation-robust.

```
block1   ConvBlock(3->16)     conv1[16,3,3,3]  bn1  conv2[16,16,3,3]  bn2
block2   ResBlock(16->32)     + downsample[32,16,1,1]+bias
block3   ResBlock(32->64)     conv{1,2} are DeformableConv2d:
                                offset_conv[18,Cin,3,3]+bias, regular_conv[Co,Cin,3,3]
block4   ResBlock(64->128)    same shape
conv1..4 [32,{16,32,64,128},1,1]     the dim/4 projections, no bias
score_head.{0,2,4,6}                 [8,128,1,1] [4,8,3,3] [4,4,3,3] [1,4,3,3], no bias
desc_head.offset_conv.0              [2M,128,3,3]+bias      (K=3 patch, VALID)
desc_head.offset_conv.2              [2M,2M,1,1]+bias
desc_head.sf_conv                    [128,128,1,1], no bias
desc_head.agg_weights                [M,128,128]
hw_grid                              [25,2]  the 5x5 soft-argmax grid
```

`offset_conv` out-channels are `2*k*k = 18`, not `3*k*k` — **no modulation
mask**, so the deformable convolution is the plain torchvision
`deform_conv2d(x, offset, w)` with `mask = None`. BatchNorm is *not* folded in
the export (8 `BatchNormalization` nodes, `running_*` present as initializers);
we fold it into the preceding conv at load, which is why the loader wants the
conv weight and its BN together.

## Module layout

```
src/aliked/                      # a sibling of src/sam/, on top of nn/
├── README.md
├── Aliked.h                     # the public surface: AlikedExtractor, LightGlueMatcher
├── model/Onnx.{h,cpp}           # the protobuf walk
├── model/Weights.{h,cpp}        # host tensors -> VRAM, BN folding, hparams from shapes
├── model/AlikedModel.cpp        # backbone, score head, DKD, SDDH
├── model/LightGlue.{h,cpp}      # 9 self/cross layers + assignment
├── shaders/aliked.slang         # SDDH, NMS, soft-argmax, the assignment tail
└── tests/aliked_test.cpp
src/sfm/feature/Extractor.{h,cpp}        # IFeatureExtractor + factory (SIFT is one impl)
src/sfm/feature/LearnedMatcher.{h,cpp}   # LightGlue behind IFeatureMatcher + factory
```

(As built. Neither sfm header includes an `aliked/` one, so `src/sfm/` still
compiles without the inference layer; the factories say so at run time.)

`ssplat_aliked` goes in `cmake/SsplatNn.cmake` next to `ssplat_sam`, with its
own `ssplat_nn_shaders(aliked …)` edge and its own declare/ensure pair — *a
static library whose only content is a static initializer is not linked*, so a
new shader directory that skips this comes back "no shader module". `SsplatSfm`
links it when `SSPLAT_BUILD_SAM AND SSPLAT_BUILD_SFM`.

Model-specific constants stay here, never in `nn/`.

## The nn/ op gap

Confirmed against the exported graph's op mix (`Conv, Selu, Sigmoid,
AveragePool, Resize, MaxPool, GridSample, TopK, Einsum, BatchNormalization`).

| op | status |
|---|---|
| `Act::Selu` | new — ALIKED's gate is `nn.SELU` everywhere. One enum row, one `apply_act` case. |
| `avgpool` (2×2 s2, 4×4 s4) | new; `maxpool2x2` exists, this is the same shape |
| `resize_bilinear` **align_corners=True** | ALIKED uses `nn.Upsample(align_corners=True)`; ours is hard-wired to torch's `align_corners=False` mapping. Add a flag; **do not** change the default, mask upsampling depends on it. |
| **deformable conv 3×3** | new, and the only substantial kernel. offset conv (existing `conv2d`, Cout=2·9) → clamp to ±max(H,W)/4 → `deform_im2col` feeding the existing chunked im2col+GEMM. Same chunking as `OpConv.cpp`, bilinear fetch instead of a direct one, **zero outside the image** (torchvision's rule, not clamp-to-edge). |
| **`grid_sample_points`** | new. Bilinear sample of a `[H,W,C]` map at N arbitrary positions → `[N,C]`. SDDH needs it; `roi_align` is close but not a substitute. |
| NMS + compaction | new, but it is `sfm/feature/Sift.h`'s pattern verbatim: 5×5 max-pool equality (`simple_nms(radius=2)`, twice), append through an atomic counter, read the count back, exact top-K on the host. A device sort is not worth it at a few thousand candidates. |
| soft-argmax | new, tiny — 5×5 patch, T=0.1, the `hw_grid` initializer is its coordinate table |
| **SDDH** | new. One workgroup per keypoint: gather the 3×3×C patch → two small matmuls → clamp → M bilinear samples → `sf_conv` → `einsum('ncp,pcd->nd')` over `agg_weights[M,C,C]` → L2 normalize. ~1 GFLOP at N=4k, M=16, C=128. |
| **LightGlue: nothing** | `linear` (bias/residual/act), `layer_norm`, `attention` (4 heads × 64), `rope` — our freqs layout `[n, hd/2, 2]` *is* LightGlue's rotary form. Only the tail needs a kernel: log-double-softmax over N₀×N₁ + mutual-NN + threshold. |

Every one lands with an `nn_ops_test` case against an independent scalar CPU
reference. That file is why the model layer can treat this API as trustworthy.

## sfm/ plumbing

1. **`IFeatureExtractor`.** `extractDirectory` news `SiftExtractor` directly
   today; a factory on the feature type also keeps sfm's Vulkan device
   uncreated when ALIKED is selected — which matters, see (7).
2. **`FeatureSet`**: `dtype = F32, dim = 128`. Already allowed (D1). ~1 MB per
   image at 2048 keypoints, the same order as SIFT at 8192.
3. **Persist the detection score.** ALIKED has no scale and no orientation.
   `PairSelection::topScaleSubset` ranks by `scale`, so with an all-zero scale
   it silently picks an arbitrary 512 keypoints — a wrong answer that looks
   like a working feature. `features.bin` goes to **v5** with a per-keypoint
   score (`Keypoint::response`, which exists and is unwritten), and the subset
   ranks by scale-or-score.
4. **Working resolution.** COLMAP's `EffMaxImageSize()` is **1600 for ALIKED**
   against 3200 for SIFT, and it means it: the aggregated map is 128 channels
   at *full* resolution, so 1600×1200 is ~1 GB fp32 for one tensor. `--quality`
   gets a per-frontend default.
5. **RGB in, `/255`, no mean/std.** `GrayImage::rgb` already exists behind
   `want_color = true`, which `extract` already sets for point colors.
6. **The coordinate chain**, copied from `aliked.cc` exactly: pad to /32 by
   replicating the right and bottom edges; keypoints come back normalized to
   [-1,1] over the *padded* extent; `px = (nx + 1) * 0.5 * (padded_w - 1) +
   0.5`; drop anything outside the original bounds. ALIKED's origin is the
   top-left pixel *center*, ours is its corner — that is what the `+0.5` is.
7. **Two Vulkan devices.** `nn::Device` and `sfm::VkContext` are separate
   contexts. Under `match` with LightGlue, pair selection (sfm) and LightGlue
   (nn) would both be live; sequence them rather than interleaving. Converging
   is sfm-port-plan phase 6 and does not block this.
8. **Matching.** Brute force stays the existing kernel: quantize the
   L2-normalized descriptors to int8 *on upload* (128 B/descriptor, exactly
   what the packed-dot path and the residency budget are sized for), keeping
   F32 on disk for LightGlue. Add the absolute-similarity threshold that
   COLMAP's defaults require. LightGlue is ~100+ GFLOP for a 2048×2048 pair —
   tens of milliseconds — so it runs **only behind pair selection**, never over
   a raw exhaustive list.
9. **Config**: `--features`, `--matcher`, `--aliked-max-features`,
   `--aliked-min-score`, `--lightglue-min-score`, `--aliked-model`. One
   `SFM_CONFIG_FIELDS` row each; the GUI's editor gets them for free.

## Phases, each with a gate

| | work | gate | |
|---|---|---|---|
| P0 | ONNX reader + weight store | every expected tensor present with the expected shape; `aliked_test` over both checkpoints | **done** |
| P1 | the new nn ops | `nn_ops_test` cases vs scalar CPU references | **done** |
| P2 | backbone + score head | score-map cosine vs an ORT reference dump on one image. 80% of the bugs die here. | **done** |
| P3 | DKD + SDDH | vs COLMAP's own extractor on ~20 images: keypoint repeatability, descriptor cosine ≥ 0.99 | **done** |
| P4 | sfm plumbing, brute force end to end | `tools/sfm/eval_poses.py` AUC against the SIFT baseline | **done** |
| P5 | LightGlue | same gate; the win should show on wide-baseline and low-texture captures | **done** |

`aliked-n32` was free after P3, as predicted: one constant read from the file,
no code change, and it hit parity on the first run. `aliked-t16` and the ONNX
brute-force matcher are deliberately not ported.

### What P0-P3 measured

`tools/aliked/compare_colmap.py` against `colmap feature_extractor
--FeatureExtraction.type ALIKED_N16ROT`, four images, both variants:
98.6-99.8% of keypoints within 1 px, mean offset below 0.003 px in both axes,
descriptor cosine 0.99995 mean. The residual is at the score threshold — peaks
within float noise of `min_score` survive on one side and not the other.
56-66 ms at ~1.1 MP against ~700 ms for COLMAP's onnxruntime CPU path. The
table lives in `src/aliked/README.md`.

### What P4-P5 measured

`ssplat sfm auto` on a 20-image wide-baseline capture, exhaustive matching, one
camera group. Registration and point count, not just match volume:

| frontend + matcher | features/img | pairs kept | inliers | registered | points |
|---|---|---|---|---|---|
| SIFT + brute force | 8192 | 68/190 | 7 703 | 55% | — |
| ALIKED + brute force | 4096 | 190/190 | 24 597 | 100% | 6 795 |
| ALIKED + LightGlue | 4096 | 75/75 (shortlist) | 52 981 | 100% | 17 343 |

LightGlue is ~53 ms for a 2029 x 2044 pair against ~3.5 ms for brute force, so
its preset switches `--pairs` to `prefilter`; exhaustive matching with it is
hours where pair selection is minutes. It matched ORT exactly on a real pair --
702 matches both sides, 100% identical partners, max score difference 0.024.

**COLMAP's matching defaults were measured and rejected.** Its ALIKED
brute-force settings are `min_cossim 0.85` with the ratio test off; on this
data an absolute 0.85 cosine rejects ~70% of mutual-nearest matches (their
median cosine is 0.726) and left 24/190 pairs. What the data wanted instead was
the ratio test kept and *loosened* to 0.92 -- a learned descriptor's second-best
distance sits much closer to its best (median ratio 0.826, just the wrong side
of SIFT's 0.8), which is the same property that makes LightGlue worth having.
`--min-similarity` still exists for anyone who wants COLMAP's shape of test.

Two things in this plan were wrong, both found by bisecting against the graph
rather than by reading the reference implementation, and both worth recording
because they are the shape of bug this port produces:

- **`ResBlock`'s second gate is outside the residual add**, not inside it.
  The plan asserted the opposite on the strength of a Selu *count* (17), which
  is consistent with both placements. Gating inside gave a detector that still
  found 16% of the same keypoints and descriptors that were uncorrelated --
  i.e. it looked like a working feature extractor. The exported graph's
  dataflow (`/block2/gate_1/Selu` consumes `/block2/Add_output_0`) settles it,
  and a node's *inputs* are worth more than its count.
- **A GEMM cannot be run in place.** `linear(feat, feat, w)` for the descriptor
  head's `sf_conv` races: every output column re-reads the whole input row.
  Elementwise ops in `nn/` do alias safely and `LinearOpts::residual` is
  explicitly allowed to, which is what made it look permissible.

## Risks, named up front

- **The deformable conv's sampling rule.** torchvision zero-fills outside the
  image and orders offsets `(dy, dx)` per tap; getting either wrong produces a
  score map that looks plausible and matches badly.
- **`grid_sample` conventions.** `align_corners` and the ±1 normalization
  appear three times (upsample, deform conv, SDDH) and are not the same
  convention in all three.
- **The coordinate chain** in (6): four frames — normalized, padded pixels,
  extraction pixels, source pixels — and `scaleKeypoints` is a fifth hop.
- **fp16.** SAM's policy is f16 for matmul/conv weights. ALIKED is 3 MB total;
  keep everything f32 until P3 passes, then measure. The activations are what
  cost memory here, not the weights.
