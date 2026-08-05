#pragma once

// Image metrics for the eval pass. All host-side, all on [0,1] float RGB in
// HWC order, matching what the Python eval computed with torchmetrics.
//
// LPIPS is deliberately absent: it needs AlexNet + VGG16, and lives in
// reference/python/eval_lpips.py, which reads the PNGs `--save-eval-images`
// writes. Everything here is arithmetic with no model behind it.

#include <cstdint>
#include <vector>

namespace spirula {

// Mean absolute error.
float image_l1(const float* a, const float* b, int64_t n);

// Peak signal-to-noise ratio, data_range = 1 (torchmetrics PeakSignalNoiseRatio).
float image_psnr(const float* a, const float* b, int64_t n);

// Structural similarity, matching torchmetrics' StructuralSimilarityIndexMeasure
// defaults: 11x11 gaussian window, sigma 1.5, data_range 1, per-channel.
// Inputs are reflect-padded by 5 and the score is averaged over the FULL HxW
// map -- torchmetrics only crops the border back off on its
// return_contrast_sensitivity path, not for the plain score.
float image_ssim(const float* a, const float* b, int H, int W, int C);

// Warp `img` to match `ref`'s colours: least squares over a quadratic
// expansion of the input pixel, re-solved `num_iters` times while updating
// which pixels count as unsaturated. This is the standard Mip-NeRF 360 /
// bilagrid colour correction that the eval's cc_* metrics are computed on.
// Both inputs are [H*W*C] in [0,1]; the result is clipped to [0,1].
std::vector<float> color_correct(const float* img, const float* ref,
                                 int64_t num_pixels, int C,
                                 int num_iters = 5, float eps = 0.5f / 255.0f);

// Same, writing into a caller-owned buffer. Scoring a 4K view otherwise spends
// more time faulting in fresh scratch than doing the arithmetic, so the eval
// pass reuses one buffer per worker across views.
void color_correct_into(const float* img, const float* ref,
                        int64_t num_pixels, int C, std::vector<float>& out,
                        int num_iters = 5, float eps = 0.5f / 255.0f);

}  // namespace spirula
