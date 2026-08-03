#pragma once
// ALIKED keypoint extraction: an RGB image in, keypoints and 128-D float
// descriptors out.
//
// The public surface of src/aliked/. Everything below it -- the ONNX reader,
// the weight store, the forward pass -- is an implementation detail, and
// nothing here exposes Vulkan, nn::Tensor or the checkpoint's layout, so
// src/sfm/ can hold one of these without depending on the inference layer's
// headers.
//
// Conventions, all inherited from COLMAP's extractor so that a reconstruction
// built on our features and one built on COLMAP's are comparable
// (docs/notes/aliked-port-plan.md):
//
//   * Input is 8-bit interleaved RGB, normalized by 1/255 and nothing else.
//   * The image is padded to a multiple of 32 by replicating its right and
//     bottom edges; the network needs it, and keypoints found in the padding
//     are dropped.
//   * Output coordinates are in the ORIGINAL image's pixels with the top-left
//     pixel's *corner* at (0, 0) -- COLMAP's convention, half a pixel off
//     ALIKED's own, which puts the top-left pixel's centre there.
//   * Descriptors are L2-normalized floats, so a dot product is a cosine.

#include <cstdint>
#include <string>
#include <vector>

namespace aliked {

struct ExtractOptions {
    // Top-scoring keypoints kept. COLMAP's AlikedExtractionOptions default.
    int max_num_features = 2048;
    // Minimum detection score, applied to the sub-pixel-refined score exactly
    // where COLMAP applies it.
    float min_score = 0.2f;
    // DKD's NMS radius and soft-argmax temperature. Not exposed as flags
    // anywhere; they are the trained detector's constants, here so the
    // forward pass has one place to read them from.
    int   nms_radius = 2;
    float temperature = 0.1f;
    // Candidate list capacity. Every local maximum of the score map is
    // collected before the top-K cut, and a 1600x1200 image has tens of
    // thousands; saturating warns rather than failing, as GPU SIFT's lists do.
    uint32_t max_candidates = 262144;
};

struct Keypoint {
    float x = 0, y = 0;   // original-image pixels, corner origin
    float score = 0;
};

struct Features {
    int width = 0, height = 0;      // the image the coordinates refer to
    int desc_dim = 0;
    std::vector<Keypoint> keypoints;
    std::vector<float>    descriptors;   // keypoints.size() * desc_dim
};

// Owns a checkpoint and the working buffers for one image at a time.
// Not thread-safe: one Extractor per thread, and they share the process-wide
// inference device.
class Extractor {
public:
    Extractor();
    ~Extractor();
    Extractor(const Extractor&) = delete;
    Extractor& operator=(const Extractor&) = delete;

    // `model` is a known id ("aliked-n16rot", "aliked-n32"), which is fetched
    // and cached, or a path to an .onnx file. Throws nn::Error on anything
    // that is not an ALIKED checkpoint.
    void load(const std::string& model);
    bool loaded() const;

    // Descriptor width of the loaded checkpoint (128 for every released one).
    int descriptorDim() const;

    // `rgb` is width*height*3 interleaved bytes.
    Features extract(const uint8_t* rgb, int width, int height,
                     const ExtractOptions& opts = {});

private:
    struct Impl;
    Impl* impl_ = nullptr;
};

}  // namespace aliked
