#define SLANG_PRELUDE_EXPORT
#include "common.cuh"

#include "SphericalHarmonics.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PixelWise.cuh"
#include "SplatTileIntersector.cuh"
#include "ProjectionEWA3DGSHetero.cuh"

#define TORCH_INDUCTOR_CPP_WRAPPER
#include <torch/extension.h>

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {

    // py::enum_<gsplat::CameraModelType>(m, "SSplatCameraModelType")
    //     .value("PINHOLE", gsplat::CameraModelType::PINHOLE)
    //     .value("ORTHO", gsplat::CameraModelType::ORTHO)
    //     .value("FISHEYE", gsplat::CameraModelType::FISHEYE)
    //     .value("FTHETA", gsplat::CameraModelType::FTHETA)
    //     .export_values();

    m.attr("TILE_SIZE") = py::int_(TILE_SIZE);

    // SphericalHarmonics.cuh
    m.def("compute_sh_forward", &compute_sh_forward_tensor);
    m.def("compute_sh_backward", &compute_sh_backward_tensor);

    // BackgroundSphericalHarmonics.cuh
    m.def("render_background_sh_forward", &render_background_sh_forward_tensor);
    m.def("render_background_sh_backward", &render_background_sh_backward_tensor);

    // PerSplatLoss.cuh
    m.def("compute_per_splat_losses_forward", &compute_per_splat_losses_forward_tensor);
    m.def("compute_per_splat_losses_backward", &compute_per_splat_losses_backward_tensor);

    // PixelWise.cuh
    m.def("blend_background_forward", &blend_background_forward_tensor);
    m.def("blend_background_backward", &blend_background_backward_tensor);

    // SplatTileIntersector.cuh
    m.def("intersect_splat_tile", &SplatTileIntersector::intersect_splat_tile);
    
    // ProjectionEWA3DGSHetero.cuh
    m.def("projection_ewa_3dgs_hetero_forward", &projection_ewa_3dgs_hetero_forward_tensor);

}
