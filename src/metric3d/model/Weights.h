#pragma once
// The Metric3D checkpoint: an .onnx file in, named device tensors out.
//
// Four things happen between those, and all four are properties of the *file*
// rather than of the forward pass:
//
//   * Linear weights are anonymous. Torch lowers nn.Linear to MatMul + Add, so
//     only the bias keeps its qualified name and the matrix arrives as
//     "onnx::MatMul_2537". The MatMul node's own name carries the module path,
//     which is the route that also works for Readout.project_learn -- it has
//     bias=False, so there is no Add to key off.
//   * ...and they are [in, out], the operand order `x @ W` wants. nn::linear
//     takes PyTorch's [out, in], so every one is transposed here.
//   * Big matrices go to the device as f16 and everything else as f32. The
//     tensor-core GEMM needs an f16 weight, and the published exports are
//     already f16, so this is a re-encode rather than a loss. Vectors --
//     LayerNorm affines, LayerScale gammas, biases, the depth anchors -- stay
//     f32: they are read once per element and rounding them buys nothing.
//   * The positional embedding stays on the HOST. It is stored for a 37x37
//     patch grid and has to be bicubically resampled to whichever grid the
//     input produces, which is a per-resolution operation, not a per-run one.
//
// Names are normalized on the way in: the `meta_arch.depth_model.` prefix goes,
// and the encoder's single BlockChunk is flattened, so the forward pass asks
// for `encoder.blocks.23.norm1.weight` and reads like the PyTorch source.

#include "metric3d/model/Hparams.h"
#include "nn/Tensor.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace metric3d {

class Weights {
public:
    Weights() = default;
    ~Weights();
    Weights(const Weights&) = delete;
    Weights& operator=(const Weights&) = delete;

    // Parses, validates, transforms and uploads. Throws nn::Error naming the
    // tensor on any missing or unexpectedly shaped weight -- a checkpoint that
    // is not Metric3D must fail here with a sentence, not later with a fault.
    void load(const std::string& onnx_path);

    bool loaded() const { return loaded_; }
    const Hparams& hparams() const { return hp_; }
    const std::string& path() const { return path_; }
    uint64_t deviceBytes() const { return device_bytes_; }

    // Throws when absent, listing near misses: a typo in a weight name is
    // otherwise a null tensor that faults somewhere else entirely.
    nn::Tensor get(const std::string& name) const;
    // printf-style, for the per-block names ("encoder.blocks.%d.norm1.weight").
    nn::Tensor getf(const char* fmt, ...) const;
    bool has(const std::string& name) const { return tensors_.count(name) != 0; }

    // The stored positional embedding, on the host. `patch_pos` is
    // [grid*grid, embed_dim] and `cls_pos` is [embed_dim]; the register tokens
    // get no positional embedding at all, which is why they are concatenated
    // after this is added and not before.
    int posGrid() const { return pos_grid_; }
    const std::vector<float>& patchPos() const { return patch_pos_; }
    const std::vector<float>& clsPos() const { return cls_pos_; }
    // The depth regressor's anchors, log-spaced over [min_depth, max_depth].
    const std::vector<float>& depthBins() const { return bins_; }

private:
    std::unordered_map<std::string, nn::Tensor> tensors_;
    Hparams                    hp_;
    std::string                path_;
    std::vector<nn::DevicePtr> chunks_;
    uint64_t                   device_bytes_ = 0;
    int                        pos_grid_ = 0;
    std::vector<float>         patch_pos_, cls_pos_, bins_;
    bool                       loaded_ = false;
};

}  // namespace metric3d
