#include "core/Tensor.h"

#include "backend/api/BackendTypes.h"
#include <cstdint>
#include <string>


/* == AUTO HEADER GENERATOR - DO NOT EDIT THIS LINE OR ANYTHING BELOW THIS LINE == */



void render_background_sh_forward(
    int w,
    int h,
    std::string camera_model,
    std::string distortion,
    int sh_degree,                       // actual SH degree (0..4)
    TorchTensorView viewmats,            // [B, 4, 4] row-major world->camera (per-batch)
    TorchTensorView intrins,             // [B, 4]
    TorchTensorView dist_coeffs,         // [B, 8]; null/empty -> zeros
    TorchTensorView sh_coeffs,           // [(sh_degree+1)^2, 3]
    TorchTensorView out_color            // [B, H, W, 3]  pre-allocated
);


void render_background_sh_backward(
    int w,
    int h,
    std::string camera_model,
    std::string distortion,
    int sh_degree,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs,
    TorchTensorView sh_coeffs,
    TorchTensorView out_color,
    TorchTensorView v_out_color,
    TorchTensorView v_sh_coeffs           // [(sh_degree+1)^2, 3] zeroed
);
