#pragma once
// Metric3D v2 monocular geometry: one image in, a depth map and a surface
// normal map out.
//
// The public surface of src/metric3d/. Everything below it -- the ONNX reader,
// the weight store, the ViT and the RAFT decoder -- is an implementation
// detail, and nothing here exposes Vulkan, nn::Tensor or the checkpoint's
// layout.
//
// Two conventions are worth stating because getting either wrong is silent:
//
//   * The input must be a PINHOLE image. The network estimates z-depth for a
//     perspective camera and was trained on nothing else; a fisheye frame fed
//     straight in produces a plausible-looking map that is wrong everywhere off
//     axis. metric3d/Warp.h is how a distorted camera reaches this.
//   * `depth` is CANONICAL, not metric: the network's camera has a focal length
//     of 1000 px, so metres = depth * fx / 1000 with `fx` in the pixels of the
//     image that was passed in. Anything that only cares about relative depth
//     can ignore that; anything that reports a distance cannot.

#include <cstdint>
#include <string>
#include <vector>

namespace metric3d {

// Which artifact to run. All three are the same architecture at three widths;
// `large` is the default everywhere because it is what the reference pipeline
// used and giant2 costs 3.3x the VRAM for a modest gain.
struct PredictOptions {
    bool want_depth = true;
    bool want_normal = true;
};

struct Prediction {
    int width = 0, height = 0;
    // [H*W] canonical depth, clamped to the network's own 0.1..200 range.
    std::vector<float> depth;
    // [H*W*3] unit normals in the camera frame, and [H*W] of the von
    // Mises-Fisher concentration the same head emits. `confidence` is filled
    // only when it was asked for; nothing in this repository consumes it yet.
    std::vector<float> normal;
    std::vector<float> confidence;
};

// Owns a checkpoint and the working buffers for one image at a time. Not
// thread-safe: one Predictor per thread, and they share the process-wide
// inference device.
class Predictor {
public:
    Predictor();
    ~Predictor();
    Predictor(const Predictor&) = delete;
    Predictor& operator=(const Predictor&) = delete;

    // `model` is a known id ("metric3d-vit-large"), which is fetched and
    // cached, or a path to an .onnx file. Throws nn::Error on anything that is
    // not a Metric3D v2 checkpoint.
    void load(const std::string& model);
    bool loaded() const;

    // The ViT's patch stride (14).
    static int patchSize();

    // The side length a prediction comes back at UNCHANGED, which is not the
    // patch stride: the decoder works on a floor(3.5 * side/14) grid and
    // upsamples it by 4, so a side of 14*(2k+1) comes back two pixels short.
    // 28 is what round-trips, and it is what a caller sizing its own input
    // should land on.
    static int sizeGranularity();

    // `rgb` is width*height*3 floats in 0..1, row-major, interleaved. Sizes
    // that are not multiples of patchSize() are cropped to one, and the
    // prediction comes back at whatever the decoder produced -- `out.width` /
    // `out.height` say what that was, which is the input size exactly when it
    // was a multiple of sizeGranularity().
    Prediction predict(const float* rgb, int width, int height,
                       const PredictOptions& opts = {});

    // Bytes this Predictor is holding on the device.
    uint64_t deviceBytes() const;

private:
    struct Impl;
    Impl* impl_ = nullptr;
};

}  // namespace metric3d
