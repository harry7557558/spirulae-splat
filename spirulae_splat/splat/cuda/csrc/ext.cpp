#define SLANG_PRELUDE_EXPORT
#include "common.cuh"

#include "SphericalHarmonics.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PixelWise.cuh"
#include "SplatTileIntersector.cuh"

#define TORCH_INDUCTOR_CPP_WRAPPER
#include <torch/extension.h>

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {

    // constants
    m.attr("BLOCK_WIDTH") = py::int_(BLOCK_WIDTH);

    // background
    m.def("render_background_sh_forward", &render_background_sh_forward_tensor);
    m.def("render_background_sh_backward", &render_background_sh_backward_tensor);
    // m.def("render_undistortion_map", &render_undistortion_map_tensor);

    // projection
    m.def("compute_sh_forward", &compute_sh_forward_tensor);
    m.def("compute_sh_backward", &compute_sh_backward_tensor);

    // misc
    m.def("compute_per_splat_losses_forward", &compute_per_splat_losses_forward_tensor);
    m.def("compute_per_splat_losses_backward", &compute_per_splat_losses_backward_tensor);
    m.def("blend_background_forward", &blend_background_forward_tensor);
    m.def("blend_background_backward", &blend_background_backward_tensor);

    m.def("intersect_splat_tile", &SplatTileIntersector::intersect_splat_tile);
    
}
