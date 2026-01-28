#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <ATen/Tensor.h>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



at::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    std::string camera_model,
    const at::Tensor &intrins,  // fx, fy, cx, cy
    const at::Tensor &rotation,  // row major 3x3
    const unsigned sh_degree,
    const at::Tensor &sh_coeffs
);


std::tuple<
    at::Tensor,  // v_rotation
    at::Tensor  // v_sh_coeffs
> render_background_sh_backward_tensor(
    const unsigned w,
    const unsigned h,
    const std::string camera_model,
    const at::Tensor &intrins,  // fx, fy, cx, cy
    const at::Tensor &rotation,  // row major 3x3
    const unsigned sh_degree,
    const at::Tensor &sh_coeffs,
    const at::Tensor &out_color,
    const at::Tensor &v_out_color
);
