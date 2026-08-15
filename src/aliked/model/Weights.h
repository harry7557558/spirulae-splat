#pragma once
// The ALIKED checkpoint: an .onnx file in, named device tensors out.
//
// Two things happen between those, and both are here rather than in the
// forward pass because both are properties of the *file*:
//
//   * BatchNorm folding. The export did not fold it (8 BatchNormalization
//     nodes, running_mean / running_var present as initializers), so every
//     conv-then-BN pair becomes one conv with a bias. The epsilon comes from
//     the node attribute, not from PyTorch's default -- nn/io/Onnx.cpp
//     reads it.
//   * Hyperparameters. `aliked-n16rot` and `aliked-n32` have byte-identical
//     graph structure and differ only in M, the number of SDDH sample
//     positions. So M is read from desc_head.agg_weights, not from the id, and
//     one code path serves both. Same for every channel count.
//
// Everything stays f32 on the device. The whole checkpoint is 2.7 MB; what
// costs memory in this model is the full-resolution activations, and halving
// 2.7 MB against that is not worth a second numeric path.

#include "aliked/Common.h"
#include "nn/Tensor.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace aliked {

// Read off the checkpoint's tensor shapes; see the class comment.
struct AlikedHparams {
    int c1 = 0, c2 = 0, c3 = 0, c4 = 0;  // encoder block widths
    int dim = 0;                          // aggregated feature width (4 * dim4)
    int dim4 = 0;                         // per-block projection width, dim / 4
    int K = 0;                            // SDDH patch size (3)
    int M = 0;                            // SDDH sample positions (16 or 32)
    int desc_dim = 0;                     // descriptor width (128)
};

class AlikedWeights {
public:
    AlikedWeights() = default;
    ~AlikedWeights();
    AlikedWeights(const AlikedWeights&) = delete;
    AlikedWeights& operator=(const AlikedWeights&) = delete;

    // Parses, validates, folds and uploads. Throws nn::Error naming the tensor
    // on any missing or unexpectedly shaped weight -- a checkpoint that is not
    // ALIKED must fail here with a sentence, not later with a fault.
    void load(const std::string& onnx_path);

    bool loaded() const { return loaded_; }
    const AlikedHparams& hparams() const { return hp_; }
    const std::string&   path() const { return path_; }
    uint64_t             deviceBytes() const { return device_bytes_; }

    // Throws when absent, listing near misses: a typo in a weight name is
    // otherwise a null tensor that faults somewhere else entirely.
    nn::Tensor get(const std::string& name) const;
    // printf-style, for the per-block names ("block%d.conv1.weight").
    nn::Tensor getf(const char* fmt, ...) const;
    bool       has(const std::string& name) const { return tensors_.count(name) != 0; }

private:
    std::unordered_map<std::string, nn::Tensor> tensors_;
    AlikedHparams hp_;
    std::string   path_;
    nn::DevicePtr blob_ = 0;
    uint64_t      device_bytes_ = 0;
    bool          loaded_ = false;
};

}  // namespace aliked
