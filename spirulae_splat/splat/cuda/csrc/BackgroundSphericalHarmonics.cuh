#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



torch::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    std::string camera_model,
    const torch::Tensor &Ks,  // row major 3x3
    const torch::Tensor &rotation,  // row major 3x3
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs
);


std::tuple<
    torch::Tensor,  // v_rotation
    torch::Tensor  // v_sh_coeffs
> render_background_sh_backward_tensor(
    const unsigned w,
    const unsigned h,
    const std::string camera_model,
    const torch::Tensor &Ks,  // row major 3x3
    const torch::Tensor &rotation,  // row major 3x3
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs,
    const torch::Tensor &out_color,
    const torch::Tensor &v_out_color
);
