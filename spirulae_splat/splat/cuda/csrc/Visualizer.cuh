#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>

#include <gsplat/Common.h>

#include "types.cuh"


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor blit_train_cameras_tensor(
    at::Tensor render_rgbs,  // [H, W, 3]
    at::Tensor render_depths,  // [H, W, 1]
    at::Tensor render_alphas,  // [H, W, 1]
    const bool view_is_fisheye,
    const at::Tensor view_intrins,  // [4]
    const at::Tensor view_viewmat,  // [4, 4]
    const CameraDistortionCoeffsTensor view_dist_coeffs,
    const at::Tensor intrins,  // [N, 4]
    const at::Tensor widths,  // [N]
    const at::Tensor heights,  // [N]
    const at::Tensor is_fisheye,  // [N]
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor camera_to_worlds,  // [N, 3, 4]
    const float camera_size
);
