#define SLANG_PRELUDE_EXPORT
#define SSPLAT_HOST_ONLY
#include "common.cuh"

#include "IntersectTile.cuh"
#include "SphericalHarmonics.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PerPixelLoss.cuh"
#include "PixelWise.cuh"
#include "FusedSSIM.cuh"
#include "SplatTileIntersector.cuh"
#include "SVHash.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionHeteroFwd.cuh"
#include "ProjectionHeteroBwd.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "RasterizationSortedEval3DFwd.cuh"
#include "RasterizationSortedEval3DBwd.cuh"
#include "Optimizer.cuh"

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

    // IntersectTile.cuh
    m.def("intersect_tile_3dgs", &intersect_tile_3dgs_tensor);
    m.def("intersect_tile_mip", &intersect_tile_mip_tensor);
    m.def("intersect_tile_3dgut", &intersect_tile_3dgut_tensor);
    m.def("intersect_tile_3dgut_sv", &intersect_tile_3dgut_sv_tensor);
    m.def("intersect_tile_opaque_triangle", &intersect_tile_opaque_triangle_tensor);
    m.def("intersect_tile_voxel", &intersect_tile_voxel_tensor);

    // SphericalHarmonics.cuh
    m.def("compute_sh_forward", &compute_sh_forward_tensor);
    m.def("compute_sh_backward", &compute_sh_backward_tensor);

    // BackgroundSphericalHarmonics.cuh
    m.def("render_background_sh_forward", &render_background_sh_forward_tensor);
    m.def("render_background_sh_backward", &render_background_sh_backward_tensor);

    // PerSplatLoss.cuh
    m.def("compute_per_splat_losses_forward", &compute_per_splat_losses_forward_tensor);
    m.def("compute_per_splat_losses_backward", &compute_per_splat_losses_backward_tensor);
    m.def("mcmc_add_noise_3dgs", &mcmc_add_noise_3dgs_tensor);

    // PerPixelLoss.cuh
    m.def("compute_per_pixel_losses_forward", &compute_per_pixel_losses_forward_tensor);
    m.def("compute_per_pixel_losses_backward", &compute_per_pixel_losses_backward_tensor);

    // PixelWise.cuh
    m.def("blend_background_forward", &blend_background_forward_tensor);
    m.def("blend_background_backward", &blend_background_backward_tensor);
    m.def("linear_rgb_to_srgb_forward", &linear_rgb_to_srgb_forward_tensor);
    m.def("linear_rgb_to_srgb_backward", &linear_rgb_to_srgb_backward_tensor);
    m.def("depth_to_normal_forward", &depth_to_normal_forward_tensor);
    m.def("depth_to_normal_backward", &depth_to_normal_backward_tensor);
    m.def("ray_depth_to_linear_depth_forward", &ray_depth_to_linear_depth_forward_tensor);
    m.def("ray_depth_to_linear_depth_backward", &ray_depth_to_linear_depth_backward_tensor);
    m.def("distort_image", &distort_image_tensor);
    m.def("undistort_image", &undistort_image_tensor);
    m.def("ppisp_forward", &ppisp_forward_tensor);
    m.def("ppisp_backward", &ppisp_backward_tensor);
    m.def("compute_ppsip_regularization_forward", &compute_ppsip_regularization_forward_tensor);
    m.def("compute_ppsip_regularization_backward", &compute_ppsip_regularization_backward_tensor);

    // FusedSSIM.cuh
    m.def("fused_ssim_forward", &fused_ssim_forward);
    m.def("fused_ssim_backward", &fused_ssim_backward);

    // SplatTileIntersector.cuh
    m.def("intersect_splat_tile_3dgs", &intersect_splat_tile_3dgs);
    m.def("intersect_splat_tile_opaque_triangle", &intersect_splat_tile_opaque_triangle);
    m.def("intersect_splat_tile_voxel", &intersect_splat_tile_voxel);

    // SVHash.cuh
    m.def("svhash_create_initial_volume", &svhashCreateInitialVolume);
    m.def("svhash_get_voxels", &svhashGetVoxels);
    m.def("svhash_split_voxels", &svhashSplitVoxels);
    
    // ProjectionFwd.cuh
    m.def("projection_3dgs_forward", &projection_3dgs_forward_tensor);
    m.def("projection_mip_forward", &projection_mip_forward_tensor);
    m.def("projection_3dgut_forward", &projection_3dgut_forward_tensor);
    m.def("projection_3dgut_sv_forward", &projection_3dgut_forward_tensor);
    m.def("projection_opaque_triangle_forward", &projection_opaque_triangle_forward_tensor);
    m.def("projection_voxel_forward", &projection_voxel_forward_tensor);

    // ProjectionBwd.cuh
    m.def("projection_3dgs_backward", &projection_3dgs_backward_tensor);
    m.def("projection_mip_backward", &projection_mip_backward_tensor);
    m.def("projection_3dgut_backward", &projection_3dgut_backward_tensor);
    m.def("projection_3dgut_backward_with_position_hessian_diagonal", &projection_3dgut_backward_with_position_hessian_diagonal_tensor);
    m.def("projection_3dgut_sv_backward", &projection_3dgut_backward_tensor);
    m.def("projection_opaque_triangle_backward", &projection_opaque_triangle_backward_tensor);
    m.def("projection_voxel_backward", &projection_voxel_backward_tensor);

    // ProjectionHeteroFwd.cuh
    m.def("projection_3dgs_hetero_forward", &projection_3dgs_hetero_forward_tensor);
    m.def("projection_mip_hetero_forward", &projection_mip_hetero_forward_tensor);
    m.def("projection_3dgut_hetero_forward", &projection_3dgut_hetero_forward_tensor);
    m.def("projection_opaque_triangle_hetero_forward", &projection_opaque_triangle_hetero_forward_tensor);

    // ProjectionHeteroBwd.cuh
    m.def("projection_3dgs_hetero_backward", &projection_3dgs_hetero_backward_tensor);
    m.def("projection_mip_hetero_backward", &projection_mip_hetero_backward_tensor);
    m.def("projection_3dgut_hetero_backward", &projection_3dgut_hetero_backward_tensor);
    m.def("projection_opaque_triangle_hetero_backward", &projection_opaque_triangle_hetero_backward_tensor);

    // RasterizationFwd.cuh and RasterizationBwd.cuh
    m.def("rasterization_3dgs_forward", &rasterize_to_pixels_3dgs_fwd);
    m.def("rasterization_3dgs_backward", &rasterize_to_pixels_3dgs_bwd);
    m.def("rasterization_mip_forward", &rasterize_to_pixels_mip_fwd);
    m.def("rasterization_mip_backward", &rasterize_to_pixels_mip_bwd);

    // RasterizationEval3DFwd.cuh and RasterizationEval3DBwd.cuh
    m.def("rasterization_3dgut_forward", &rasterize_to_pixels_3dgut_fwd);
    m.def("rasterization_3dgut_backward", &rasterize_to_pixels_3dgut_bwd);
    m.def("rasterization_3dgut_backward_with_hessian_diagonal", &rasterize_to_pixels_3dgut_bwd_with_hessian_diagonal);
    m.def("rasterization_3dgut_sv_forward", &rasterize_to_pixels_3dgut_fwd);
    m.def("rasterization_3dgut_sv_backward", &rasterize_to_pixels_3dgut_bwd);
    m.def("rasterization_opaque_triangle_forward", &rasterize_to_pixels_opaque_triangle_sorted_fwd);
    m.def("rasterization_opaque_triangle_backward", &rasterize_to_pixels_opaque_triangle_sorted_bwd);
    m.def("rasterization_voxel_forward", &rasterize_to_pixels_voxel_eval3d_fwd);
    m.def("rasterization_voxel_backward", &rasterize_to_pixels_voxel_eval3d_bwd);

    // Optimizer.cuh
    m.def("fused_adam", &fused_adam);
    m.def("fused_adam_multi", &fused_adam_multi);
    m.def("fused_newton", &fused_newton);
    m.def("fused_newton_multi", &fused_newton_multi);
    m.def("fused_3dgs2tr_mean_optim", &fused_3dgs2tr_mean_optim);
}
