#pragma once

// Camera-distortion inputs shared by the backend parity tools: the compiled
// (camera model, tier) set and one coefficient row per (tier, camera).

#include <core/Common.cuh>

#include <cstdint>
#include <vector>

namespace dist_fixture {

inline const char* const kTierNames[4] = {"NONE", "OPENCV", "THIN_PRISM",
                                          "RATIONAL"};

// Tiers compiled for each camera model, in CameraModelType order; see
// camera_distortion_is_compiled in core/CameraModel.h. Rows are padded with
// tier 0 so a fixed-width table indexes cleanly -- read only kNumTiers[m].
inline constexpr int kNumTiers[4] = {4, 3, 3, 1};
inline constexpr int kTiers[4][4] = {
    {0, 1, 2, 3},   // PINHOLE
    {0, 1, 2, 0},   // FISHEYE
    {0, 1, 2, 0},   // EQUISOLID
    {0, 0, 0, 0},   // EQUIRECTANGULAR
};

// [tier][camera][kCameraDistortionParams], each row written in ITS OWN tier's
// coefficient order (core/CameraModel.h) -- a row is only meaningful under the
// tier it was built for. Every tier fills the terms the cheaper tiers cannot
// express (tangential p1/p2, thin-prism sx1/sy1, the rational denominator
// k4..k6), so sweeping the tiers sweeps the distortion math and not just the
// dispatch. Magnitudes stay mild: a strong lens flips borderline culls, which
// rewrite whole rows of the comparison.
inline std::vector<float> distortion_rows(int64_t C) {
    static const float row[4][kCameraDistortionParams] = {
        {},
        {0.05f, -0.01f, 0.0012f, -0.0020f},
        {0.05f, -0.01f, 0.0020f, -0.0005f,
         0.0012f, -0.0020f, 0.0015f, -0.0008f},
        {0.05f, -0.01f, 0.0020f, 0.0300f, -0.0040f, 0.0010f,
         0.0012f, -0.0020f},
    };
    std::vector<float> out((size_t)4 * C * kCameraDistortionParams, 0.0f);
    for (int t = 0; t < 4; t++)
        for (int64_t c = 0; c < C; c++) {
            const float s = 1.0f - 0.3f * (float)(c % 3);
            for (int k = 0; k < kCameraDistortionParams; k++)
                out[((size_t)t * C + c) * kCameraDistortionParams + k] =
                    s * row[t][k];
        }
    return out;
}

// Offset in floats of camera `c`'s row under `tier`, into distortion_rows(C).
inline int64_t row_offset(int tier, int64_t C, int64_t c = 0) {
    return ((int64_t)tier * C + c) * kCameraDistortionParams;
}

}  // namespace dist_fixture
