#pragma once
// Turning MoGe's affine point map into a depth map.
//
// The network's `points` is correct up to one unknown z shift: the camera-space
// point is (X, Y, Z + shift), and `shift` is what makes the map project through
// a pinhole of the given focal length. Recovering it is a one-dimensional least
// squares -- MoGe's own `solve_optimal_shift` -- and when the focal is unknown
// too it comes out of the same search, because the optimal focal at a given
// shift has a closed form.
//
// Host math with no Vulkan in it, so `moge_test --recover` can check it against
// an analytic plane.

namespace moge {

struct Recovered {
    float shift = 0.0f;
    // Relative to half the image diagonal, which is the unit MoGe's uv grid is
    // written in: focal = 2 * fx_pixels / hypot(width, height).
    float focal = 1.0f;
    bool  solved = false;
};

// `points` is [n][3], `uv` is [n][2]; a positive `focal` fixes it, zero recovers
// both. Minimized by a scan and golden section rather than the
// Levenberg-Marquardt scipy runs, which walks into the barrier at -min(z).
Recovered recover_shift(const float* points, const float* uv, int n, float focal);

}  // namespace moge
