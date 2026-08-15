#pragma once
// The loaded model and the two halves of its forward pass.
//
// Split across Encoder.cpp and Decoder.cpp because they are two different
// networks bolted together -- a DINOv2 ViT and a RAFT-style iterative decoder --
// and each is meant to be read against its own PyTorch source
// (ViT_DINO_reg.py and RAFTDepthNormalDPTDecoder5.py in YvanYin/Metric3D).
//
// There is no computation graph: an op runs when it is called and writes into
// an arena tensor. The arena is a bump allocator that refuses to grow while
// anything is live, so plan_arena_bytes() has to bound the whole pass up front.

#include "metric3d/model/Weights.h"
#include "nn/Tensor.h"
#include "nn/vk/Memory.h"

#include <string>
#include <vector>

namespace metric3d {

// What the encoder hands the decoder: four token tensors, each
// [1 + registers + gh*gw, embed_dim]. With a final LayerNorm all four are the
// same tensor; without one they are four intermediate blocks.
struct Features {
    nn::Tensor tokens[4];
    int64_t    gh = 0, gw = 0;
};

struct Model {
    Weights       weights;
    nn::vk::Arena arena{"metric3d"};

    // The positional embedding resampled for one patch grid, and the grid it
    // was built for. Bicubic resampling is a host operation done once per input
    // size, not per image, and a dataset is usually one size.
    nn::DevicePtr pos_blob = 0;
    nn::Tensor    pos_patch;   // [gh*gw, D]
    nn::Tensor    pos_cls;     // [D]
    int64_t       pos_gh = 0, pos_gw = 0;

    ~Model();

    const Hparams& hp() const { return weights.hparams(); }

    // `image` is [H, W, 3] f32 on the device, already ImageNet-normalized.
    Features encode(const nn::Tensor& image);

    // Writes `depth` [Hf*Wf], `normal` [Hf*Wf*3] and `kappa` [Hf*Wf] arena
    // tensors, at the 1/4-grid resolution times the convex upsample factor --
    // which is the input resolution. A null out-tensor is not computed.
    void decode(const Features& f, nn::Tensor* depth, nn::Tensor* normal,
                nn::Tensor* kappa);

    void ensurePosEmbed(int64_t gh, int64_t gw);

    uint64_t planArenaBytes(int64_t H, int64_t W) const;
};

// SS_METRIC3D_DUMP=<dir> writes one .npy per stage for
// tools/metric3d/compare_ort.py. Off unless the variable is set; the check is
// a getenv on the first call and a null pointer test after that.
bool  dump_enabled();
void  dump_tensor(const char* name, const nn::Tensor& t,
                  const std::vector<int64_t>& shape);

}  // namespace metric3d
