#pragma once

// Fit an arbitrary source camera onto a compiled (CameraModelType,
// CameraDistortionType) pair on the CPU.
//
// Used when a dataset carries a lens model no tier represents exactly (COLMAP
// FOV, Metashape affinity/skew, EUCM/UCM, ...). The fit is in PIXELS against
// the source's own pixel map over its whole valid domain, so the fitted camera
// keeps the source image size and the same solid angle -- no inscribed-rectangle
// crop. `fit_camera` is pure; cache the result per camera group.

#include "core/CameraModel.h"

#include <functional>

namespace dsfit {

// A fit closer than this reproduces the source camera to better than bilinear
// resampling can resolve, so a parser that reaches it can keep the images as
// they are instead of re-distorting them.
constexpr double kExactFitPx = 0.1;

// Forward map of the SOURCE camera: a view-space direction (need not be
// normalized) -> pixel. false where the source has no image: behind the camera,
// or outside its valid domain.
using SourceProject = std::function<bool(double x, double y, double z,
                                         double* u, double* v)>;

struct FitTarget {
    CameraModelType      model = CameraModelType::PINHOLE;
    CameraDistortionType distortion = CameraDistortionType::ThinPrism;
    float fx = 0.0f, fy = 0.0f, cx = 0.0f, cy = 0.0f;
    float coeffs[kCameraDistortionParams] = {};
};

struct FitOptions {
    int    grid = 64;                // Chebyshev nodes per axis
    bool   refine_intrinsics = true; // false: fx/fy/cx/cy from the undistorted fit only
    double max_rms_px = 0.5;         // decides FitResult::ok, nothing else
    int    lawson_iters = 8;         // 0 = plain Chebyshev-weighted least squares
    // Which of the 8 ThinPrism slots the fit may use, low bit = k1, in the
    // order k1 k2 k3 k4 p1 p2 sx1 sy1. Retrying with the high-order terms
    // masked off is how a non-invertible fit is recovered.
    unsigned coeff_mask = 0xFFu;
};

struct FitResult {
    FitTarget target;
    // Measured on an independent uniform grid twice as dense as the fit grid,
    // not on the nodes the fit optimized.
    double rms_px = 0.0;
    double max_px = 0.0;   // the near-minimax objective
    int    samples = 0;
    double fov_deg = 0.0;  // full angle the source images, 2x the largest theta
    bool   invertible = true;  // fitted distortion passes is_valid_distortion everywhere sampled
    bool   ok = false;         // solved, invertible, and rms_px <= max_rms_px
};

// Fit `target_model` + ThinPrism coefficients to `src` over a `width` x `height`
// image. Returns ok == false for EQUIRECTANGULAR (carries no distortion) and
// when the source has no image around the axis.
FitResult fit_camera(const SourceProject& src, int width, int height,
                     CameraModelType target_model,
                     const FitOptions& opt = FitOptions());

// Same, choosing the camera model from the source's measured field of view and
// degrading the coefficient set until the distortion is invertible.
//
// Always returns a camera the engine can train with -- a dataset that took
// hours to reconstruct must not fail to load over a lens model. Read `max_px`
// for how faithful the choice is; `ok` is false only when the source projects
// nothing at all, and even then the target is a usable (if arbitrary) camera.
FitResult fit_camera_auto(const SourceProject& src, int width, int height,
                          const FitOptions& opt = FitOptions());

}  // namespace dsfit
