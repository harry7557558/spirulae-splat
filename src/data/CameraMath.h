#pragma once

// The camera models on the CPU: lens distortion both ways, a pixel to a ray,
// and a ray back to a pixel.
//
// Host mirror of shaders/projection_utils.slang, exactly as
// data/SourceCamera.h is the host mirror of shaders/camera_source.slang. The
// two cannot be shared -- one runs on the GPU -- so this is the ONE host copy
// and everything that needs it on the CPU calls in here: the GUI's frustum
// wireframe and `spirula geometry`'s camera resampling.
//
// Coordinates are NORMALIZED image coordinates ((px - cx) / fx), not pixels,
// and rays are CV camera space (+Z forward, +Y down). The distortion tier
// travels with the coefficients because the slot order is per-tier; see
// core/CameraModel.h.

namespace camhost {

// Any coefficient set to something other than zero. A tier that is nominally
// present but all-zero is not worth a Newton solve.
bool has_distortion(int tier, const float* d);

// Lens distortion alone, in the tier's own slot order. `distort_lens` is
// closed form; `undistort_lens` is Newton, and returns false when the forward
// map folds over itself or the solve does not converge.
void distort_lens(double u, double v, int tier, const float* d, double out[2]);
bool undistort_lens(double u, double v, int tier, const float* d, double out[2]);

// The full model: perspective-normalized coordinates -> the source camera's
// normalized coordinates, i.e. the fisheye or equisolid angular mapping and
// then the lens. EQUIRECTANGULAR passes through -- its coordinates are already
// (lon, lat) and it carries no lens.
void distort_point(double u, double v, int model, int tier, const float* d, double out[2]);

// The inverse: the source camera's normalized coordinates -> perspective
// normalized coordinates. False where the model has no image for that point.
bool undistort_point(double u, double v, int model, int tier, const float* d,
                     double out[2]);

// The source camera's normalized coordinates -> a UNIT ray. False outside the
// model's valid domain, including past the angle where it folds back.
bool generate_ray(double u, double v, int model, int tier, const float* d, double out[3]);

// And back: a ray (need not be normalized) -> the source camera's normalized
// coordinates. False for a direction the camera cannot see. This is the
// forward projection the wide-to-pinhole resampling gathers through.
bool project_ray(const double ray[3], int model, int tier, const float* d, double out[2]);

// What fraction of this frame ONE pinhole undistortion at the same intrinsics
// keeps. Blind to the distortion coefficients and to cx / cy on purpose: a
// threshold that flips under -ffast-math trains differently on two machines.
double pinhole_coverage(int model, int width, int height, double fx, double fy);

// Whether `spirula geometry --split auto` resamples this lens into pinhole
// faces, and so whether the depth map it wrote is RAY depth. The trainer's
// --input-depth-is-ray-depth default asks this rather than re-spelling it.
bool splits_to_pinhole_faces(int model, int width, int height, double fx, double fy);

}  // namespace camhost
