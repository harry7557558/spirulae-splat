#include <cuda.h>
#include <cuda_runtime.h>
#include <cstdint>

#include <torch/types.h>


// refactored function arguments

#define _ARGS_render_background_sh_forward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const float* rotation,  /* row major 3x3 */ \
    const unsigned sh_degree, \
    const float3* __restrict__ sh_coeffs_float3, \
    float3* __restrict__ out_img

#define _ARGS_render_background_sh_backward_kernel \
    const dim3 tile_bounds, \
    const dim3 img_size, \
    const float4 intrins, \
    const float2* __restrict__ undistortion_map, \
    const float* rotation,  /* row major 3x3 */ \
    const unsigned sh_degree, \
    const float3* __restrict__ sh_coeffs_float3, \
    const float3* __restrict__ out_color, \
    const float3* __restrict__ v_out_color, \
    float3* __restrict__ v_rotation, \
    float3* __restrict__ v_sh_coeffs


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



torch::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &rotation,
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
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &rotation,
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs,
    const torch::Tensor &out_color,
    const torch::Tensor &v_out_color
);
