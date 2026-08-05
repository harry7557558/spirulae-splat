# ALIKED (`src/aliked/`)

The learned SfM frontend: ALIKED keypoints and 128-D float descriptors, on the
inference layer (`src/nn/`), with no onnxruntime, no PyTorch and no converter.

Status: **done and wired in.** The extractor and LightGlue both match COLMAP's
ONNX implementations, and `spirula sfm --features aliked-n16rot --matcher
lightglue` runs the whole pipeline.

## Rules

- **Vulkan only**, like `src/sfm/` and `src/sam/`. Built with
  `SS_BUILD_SAM` (which is what builds `ss_nn`).
- **Model-specific things live here, never in `nn/`.** The general ops this
  needed — deformable convolution, point-wise grid sample, average pooling, row
  L2 normalization, SELU — went into `nn/` and are tested there. What is left
  in `shaders/aliked.slang` is only what could not be general: the detector's
  suppression rule, its soft-argmax, and two coordinate conversions.
- **The checkpoint is COLMAP's, fetched from COLMAP's release URL and parsed in
  process.** We host nothing and convert nothing.

## The checkpoint is an ONNX file

`model/Onnx.cpp` is a ~300-line varint walk of protobuf that reads
*initializers only* — the graph structure is hard-coded in `AlikedModel.cpp`,
so nothing has to understand a node, an operator or a shape rule.

The point of this is parity. Both sides run the same bytes, so a difference
against `colmap feature_extractor --FeatureExtraction.type ALIKED_N16ROT` is
ours and can be bisected, which is exactly how the two bugs below were found.
It also means there is no export to keep in step with an upstream re-export.

`aliked-n16rot` and `aliked-n32` have **byte-identical graph structure** — 2519
nodes, 69 initializers each — and differ only in `M`, the number of SDDH sample
positions, which is `desc_head.agg_weights.shape[0]`. So there is one code
path, `M` comes from the file, and the variant is a URL. `n16rot` is
architecturally `n16`, trained rotation-robust.

BatchNorm is *not* folded in the export (8 `BatchNormalization` nodes,
`running_*` present as initializers); `model/Weights.cpp` folds each into the
conv before it, with the epsilon read from the node rather than assumed.

## Layout

```
Aliked.h              the public surface: Extractor, ExtractOptions, Features
Common.h              the names borrowed from nn/, as sam/Common.h does
model/Onnx.{h,cpp}    the protobuf walk
model/Fetch.{h,cpp}   URL + SHA-256 table, cache, download through curl
model/Weights.{h,cpp} initializers -> device tensors; BN folding; hparams
model/AlikedModel.cpp the forward pass, read against nets/aliked.py
model/LightGlue.{h,cpp}  the matcher: 9 self/cross layers + one assignment
shaders/aliked.slang  NMS, soft-argmax, the SDDH coordinate conversions,
                        and LightGlue's assignment reduction
tests/aliked_test.cpp checkpoint shapes, SHA-256, extraction and matching
```

On the SfM side, `sfm/feature/Extractor.h` and `sfm/feature/LearnedMatcher.h`
are the seams this plugs into; neither includes an `aliked/` header, so
`src/sfm/` still builds without the inference layer and says so at run time.

## Five conventions that are not guessable

Each of these was wrong at some point, and each failure looked like a working
network rather than a broken one.

1. **`ResBlock`'s second gate is outside the residual add** —
   `gate(bn2(conv2(gate(bn1(conv1(x))))) + downsample(x))`. Reading the module
   source suggests the gate sits inside the branch; the exported graph settles
   it (`/block2/gate_1/Selu` takes `/block2/Add_output_0`). Gating inside gave
   16% keypoint agreement instead of 99.6%, with uncorrelated descriptors.
2. **Two offset layouts in one network.** The deformable convolutions use
   torchvision's `(dy, dx)` interleaved per tap. The descriptor head's offset
   conv emits `2M` channels that ALIKED reads as `view(N, 2, M)` — the first
   `M` are all the x components, the second `M` all the y. Swapping either
   produces plausible-looking descriptors that match badly.
3. **`align_corners=True`** on every upsample and every normalized coordinate
   downstream of one — which is *not* the convention `nn::resize_bilinear`
   defaults to, because mask upsampling needs the other one.
4. **LightGlue's fused qkv is `[head][dim][3]`** — q, k and v interleaved per
   element, which no stride can express. The projection's output rows are
   permuted at load so it becomes `[3][head][dim]`, the layout
   `nn::attention`'s q/k/v strides already address.
5. **LightGlue's assignment score contains `sim` twice**, once from the row
   log-softmax and once from the column one. The arg-max is over
   `2*sim[i][j] - lse_row[i] - lse_col[j] + ls0[i] + ls1[j]`; dropping the
   factor changes which column wins, not just the value, and pushes every
   score far enough below zero that nothing matches at all.

## Memory

The aggregated feature map is `dim` (128) channels at **full resolution**, so
it dominates everything: ~1 GB at COLMAP's default working size for this
extractor (1600 px, against 3200 for SIFT). `plan_arena_bytes` sizes the whole
forward pass up front because the arena will not grow while anything is live.

Keeping that map in f16 would halve it and is the obvious next move; it is not
done because `l2_normalize_rows` would need an f16 output path, i.e. a second
numeric path under a kernel that currently has one.

## Testing

```bash
./build/aliked_test                    # cached checkpoints, or SKIP
./build/aliked_test --fetch            # download from COLMAP's releases
./build/aliked_test --image IMG.jpg --out /tmp/ours.bin
./build/aliked_test --match /tmp/a.bin /tmp/b.bin   # LightGlue on two dumps
```

The checkpoint gate is strict about *shapes* — we do not own these weights and
cannot embed a golden copy, so what it can check is that every tensor the
forward pass will ask for exists at the width the rest of the model assumes.

The gate that matters is parity against COLMAP, which needs a COLMAP built with
ONNX support:

```bash
colmap feature_extractor --database_path /tmp/db.db --image_path IMAGES \
    --FeatureExtraction.type ALIKED_N16ROT
./build/aliked_test --image IMAGES/x.jpg --out /tmp/ours.bin
python3 tools/aliked/compare_colmap.py /tmp/db.db /tmp/ours.bin
```

Measured on four images, `n16rot` and `n32`, on an RTX 5070 Laptop:

| | keypoints (ours / COLMAP) | matched ≤1 px | mean offset | descriptor cosine |
|---|---|---|---|---|
| n16rot | 2029 / 2029 | 99.6% | (+0.0002, −0.0001) px | 0.999983 |
| n16rot | 2044 / 2044 | 99.8% | (−0.0004, −0.0002) px | 0.999982 |
| n16rot | 1071 / 1079 | 99.2% | (+0.0022, −0.0010) px | 0.999953 |
| n16rot | 1541 / 1529 | 98.6% | (+0.0023, +0.0005) px | 0.999971 |
| n32    | 2048 / 2048 | 99.8% | (−0.0002, −0.0002) px | 0.999978 |

The residual disagreement is at the score threshold: peaks within float noise
of `min_score` survive on one side and not the other. Extraction is 56–66 ms at
~1.1 MP, against ~700 ms for COLMAP's onnxruntime CPU path on the same image.

`nn_ops_test`'s "Learned frontend" section covers the ops this needed, each
against an independent scalar CPU reference.

## LightGlue

Nine transformer layers, each self-attention within an image then
cross-attention between the two, and one assignment head. It needed **no new
general ops**: `linear`, `layer_norm`, `attention` and `rope` cover the forward
pass, and our RoPE frequency layout already *is* LightGlue's rotary encoding.

The assignment tail is the only kernel, and it exists to avoid materializing a
second N0 x N1 matrix (16 MB at 2048 keypoints each). Every term of the score
separates into a row constant and a column constant, so the two log-sum-exps
are vectors and the arg-max applies them on the fly -- three passes over `sim`
instead of five.

The export runs the assignment on the last layer only: LightGlue's early exit
and token pruning are not in the graph, so they are not implemented. Both are
speedups, not behaviour.

Measured against ORT on a real pair: **702 matches on both sides, 100%
identical partners**, max score difference 0.024 (fp32 accumulation order over
nine layers). 53 ms for a 2029 x 2044 pair.

## Not done yet

1. **f16 for the aggregated map**, per "Memory" above.
2. LightGlue's early exit and width pruning, worth roughly 2x and absent from
   the exported graph.
3. The candidate list is collected and top-K'd on the host. That is the same
   shape as GPU SIFT's stage boundary and has not been measured as a cost.
