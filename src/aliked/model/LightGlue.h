#pragma once
// LightGlue: a learned matcher for ALIKED keypoints.
//
// Two sets of keypoints and descriptors in, a list of correspondences out. The
// network is nine transformer layers, each doing self-attention within each
// image and cross-attention between them, followed by one assignment head.
//
// It needed NO new general ops -- linear, layer_norm, attention and rope cover
// the whole forward pass, and our RoPE frequency layout is already the one
// LightGlue's rotary encoding uses. Only the assignment tail is here, in
// shaders/aliked.slang, and only because it is a reduction over a matrix that
// should never be materialized.
//
// Cost is the thing to know before using it: roughly 100 GFLOP for a
// 2048 x 2048 pair, i.e. tens of milliseconds. It is a matcher for a
// SHORTLIST, never for a raw exhaustive pair list -- run pair selection first.

#include <cstdint>
#include <string>
#include <vector>

namespace aliked {

struct MatchOptions {
    // Assignment confidence a pair must reach. COLMAP's
    // LightGlueONNXMatchingOptions default.
    float min_score = 0.1f;
};

struct Match {
    uint32_t i = 0, j = 0;   // indices into set 0 and set 1
    float    score = 0;
};

// One image's input. Keypoints are in pixels with the top-left pixel's CORNER
// at (0, 0) -- this repository's convention, which is half a pixel off
// LightGlue's own; the conversion happens inside.
struct MatchInput {
    const float* keypoints = nullptr;    // [n, 2] xy
    const float* descriptors = nullptr;  // [n, dim], L2-normalized
    uint32_t     n = 0;
    int          width = 0, height = 0;  // the image the keypoints refer to
};

class Matcher {
public:
    Matcher();
    ~Matcher();
    Matcher(const Matcher&) = delete;
    Matcher& operator=(const Matcher&) = delete;

    // "aliked-lightglue" (fetched and cached) or a path to an .onnx file.
    void load(const std::string& model);
    bool loaded() const;
    int  descriptorDim() const;

    std::vector<Match> match(const MatchInput& a, const MatchInput& b,
                             const MatchOptions& opts = {});

private:
    struct Impl;
    Impl* impl_ = nullptr;
};

}  // namespace aliked
