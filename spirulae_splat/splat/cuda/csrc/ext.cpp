#define SLANG_PRELUDE_EXPORT
#define SSPLAT_HOST_ONLY
#include "common.cuh"

#include "SphericalHarmonics.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PerPixelLoss.cuh"
#include "PixelWise.cuh"
#include "SplatTileIntersector.cuh"
#include "Projection.cuh"
#include "ProjectionEval3D.cuh"
#include "ProjectionEWA3DGSHetero.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "RasterizationSortedEval3DFwd.cuh"
#include "RasterizationSortedEval3DBwd.cuh"

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

    // PerPixelLoss.cuh
    m.def("compute_per_pixel_losses_forward", &compute_per_pixel_losses_forward_tensor);
    m.def("compute_per_pixel_losses_backward", &compute_per_pixel_losses_backward_tensor);

    // PixelWise.cuh
    m.def("blend_background_forward", &blend_background_forward_tensor);
    m.def("blend_background_backward", &blend_background_backward_tensor);
    m.def("depth_to_normal_forward", &depth_to_normal_forward_tensor);
    m.def("depth_to_normal_backward", &depth_to_normal_backward_tensor);
    m.def("ray_depth_to_linear_depth_forward", &ray_depth_to_linear_depth_forward_tensor);
    m.def("ray_depth_to_linear_depth_backward", &ray_depth_to_linear_depth_backward_tensor);
    m.def("distort_image", &distort_image_tensor);
    m.def("undistort_image", &undistort_image_tensor);

    // SplatTileIntersector.cuh
    m.def("intersect_splat_tile_3dgs", &intersect_splat_tile_3dgs);
    m.def("intersect_splat_tile_opaque_triangle", &intersect_splat_tile_opaque_triangle);
    
    // Projection.cuh
    m.def("projection_ewa_3dgs_forward", &projection_ewa_3dgs_forward_tensor);
    m.def("projection_ewa_3dgs_backward", &projection_ewa_3dgs_backward_tensor);
    m.def("projection_opaque_triangle_forward", &projection_opaque_triangle_forward_tensor);
    m.def("projection_opaque_triangle_backward", &projection_opaque_triangle_backward_tensor);

    // ProjectionEval3D.cuh
    m.def("projection_ewa_3dgs_eval3d_forward", &projection_ewa_3dgs_eval3d_forward_tensor);
    m.def("projection_ewa_3dgs_eval3d_backward", &projection_ewa_3dgs_eval3d_backward_tensor);
    m.def("projection_opaque_triangle_eval3d_forward", &projection_opaque_triangle_eval3d_forward_tensor);
    m.def("projection_opaque_triangle_eval3d_backward", &projection_opaque_triangle_eval3d_backward_tensor);

    // ProjectionEWA3DGSHetero.cuh
    m.def("projection_ewa_3dgs_hetero_forward", &projection_ewa_3dgs_hetero_forward_tensor);
    m.def("projection_ewa_3dgs_hetero_backward", &projection_ewa_3dgs_hetero_backward_tensor);
    m.def("projection_opaque_triangle_hetero_forward", &projection_opaque_triangle_hetero_forward_tensor);
    m.def("projection_opaque_triangle_hetero_backward", &projection_opaque_triangle_hetero_backward_tensor);

    // RasterizationFwd.cuh and RasterizationBwd.cuh
    m.def("rasterization_3dgs_forward", &rasterize_to_pixels_3dgs_fwd);
    m.def("rasterization_3dgs_backward", &rasterize_to_pixels_3dgs_bwd);
    m.def("rasterization_opaque_triangle_forward", &rasterize_to_pixels_opaque_triangle_fwd);
    m.def("rasterization_opaque_triangle_backward", &rasterize_to_pixels_opaque_triangle_bwd);

    // RasterizationEval3DFwd.cuh and RasterizationEval3DBwd.cuh
    m.def("rasterization_3dgs_eval3d_forward", &rasterize_to_pixels_3dgs_eval3d_fwd);
    m.def("rasterization_3dgs_eval3d_backward", &rasterize_to_pixels_3dgs_eval3d_bwd);
    // m.def("rasterization_opaque_triangle_eval3d_forward", &rasterize_to_pixels_opaque_triangle_eval3d_fwd);
    // m.def("rasterization_opaque_triangle_eval3d_backward", &rasterize_to_pixels_opaque_triangle_eval3d_bwd);
    m.def("rasterization_opaque_triangle_eval3d_forward", &rasterize_to_pixels_opaque_triangle_sorted_eval3d_fwd);
    m.def("rasterization_opaque_triangle_eval3d_backward", &rasterize_to_pixels_opaque_triangle_sorted_eval3d_bwd);

}
