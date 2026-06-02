#define SLANG_PRELUDE_EXPORT
#include "common.cuh"
#include "Engine.h"

#include "IntersectTile.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PerPixelLoss.cuh"
#include "PixelWise.cuh"
#include "FusedSSIM.cuh"
#include "SplatTileIntersector.cuh"
#include "SVHash.cuh"
// #include "Projection.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "RasterizationSortedEval3DFwd.cuh"
#include "RasterizationSortedEval3DBwd.cuh"
#include "Optimizer.cuh"
#include "FusedProjectionBwdOptim.cuh"
#include "Densify.cuh"
#include "BilagridUtils.cuh"
#include "Visualizer.cuh"

#define TORCH_INDUCTOR_CPP_WRAPPER
#include <torch/extension.h>

// at::Tensor convenience wrappers — delegate to TorchTensorView-based functions
// needed by Primitive*.cuh code inside #ifndef NO_TORCH blocks
static inline TorchTensorView _to_tv(const at::Tensor& t) {
    return TorchTensorView(
        (uint64_t)t.data_ptr(),
        (uint32_t)t.element_size(),
        t.sizes().vec()
    );
}

inline void set_zero_tensor(at::Tensor& x) {
    set_zero_tensor(_to_tv(x));
}

inline at::Tensor zeros_like_tensor(const at::Tensor& x) {
    at::Tensor y = at::empty_like(x);
    set_zero_tensor(y);
    return y;
}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {

    // py::enum_<ssplat::CameraModelType>(m, "SSplatCameraModelType")
    //     .value("PINHOLE", ssplat::CameraModelType::PINHOLE)
    //     .value("ORTHO", ssplat::CameraModelType::ORTHO)
    //     .value("FISHEYE", ssplat::CameraModelType::FISHEYE)
    //     .value("FTHETA", ssplat::CameraModelType::FTHETA)
    //     .export_values();

    m.attr("TILE_SIZE") = py::int_(TILE_SIZE);

    // IntersectTile.cuh
    m.def("intersect_tile", &do_intersect_tile_generic);
    // m.def("intersect_tile_3dgs", &intersect_tile_3dgs_tensor);
    // m.def("intersect_tile_mip", &intersect_tile_3dgs_tensor);
    // m.def("intersect_tile_3dgut", &intersect_tile_3dgut_tensor);
    // m.def("intersect_tile_3dgut_sv", &intersect_tile_3dgut_tensor);
    // m.def("intersect_tile_opaque_triangle", &intersect_tile_opaque_triangle_tensor);
    // m.def("intersect_tile_voxel", &intersect_tile_voxel_tensor);

    // BackgroundSphericalHarmonics.cuh — TorchTensorView API. Outputs are
    // caller-allocated; the C++ side fills them in place.
    m.def("render_background_sh_forward",  &render_background_sh_forward);
    m.def("render_background_sh_backward", &render_background_sh_backward);

    #if 0
    // PerSplatLoss.cuh
    m.def("compute_per_splat_losses_forward", &compute_per_splat_losses_forward);
    m.def("compute_per_splat_losses_backward", &compute_per_splat_losses_backward);
    m.def("compute_per_splat_losses_backward_with_hessian_diagonal", &compute_per_splat_losses_backward_with_hessian_diagonal);
    #endif

    // PerPixelLoss.cuh
    m.def("compute_multi_scale_per_pixel_losses", &compute_multi_scale_per_pixel_losses);

    // PixelWise.cuh
    m.def("uint8_image_to_float", &uint8_image_to_float_tensor);
    m.def("uint16_image_to_float", &uint16_image_to_float_tensor);
    m.def("rendered_depth_to_expected_depth_forward", &rendered_depth_to_expected_depth_forward);
    m.def("rendered_depth_to_expected_depth_backward", &rendered_depth_to_expected_depth_backward);
    m.def("blend_background_noise_forward", &blend_background_noise_forward);
    m.def("blend_background_noise_backward", &blend_background_noise_backward);
    m.def("blend_background_forward", &blend_background_forward);
    m.def("blend_background_backward", &blend_background_backward);
    m.def("rgb_to_srgb_forward", &rgb_to_srgb_forward);
    m.def("rgb_to_srgb_backward", &rgb_to_srgb_backward);
    m.def("depth_to_points_forward", &depth_to_points_forward);
    m.def("depth_to_points_backward", &depth_to_points_backward);
    m.def("depth_to_normal_forward", &depth_to_normal_forward_tv);
    m.def("depth_to_normal_backward", &depth_to_normal_backward_tv);
    m.def("depth_normal_loss_forward", &depth_normal_loss_forward);
    m.def("depth_normal_loss_backward", &depth_normal_loss_backward);
    m.def("ray_depth_to_linear_depth_forward", &ray_depth_to_linear_depth_forward);
    m.def("ray_depth_to_linear_depth_backward", &ray_depth_to_linear_depth_backward);
    m.def("distort_image", &distort_image_tensor);
    m.def("undistort_image", &undistort_image_tensor);
    m.def("warp_image_wide_to_pinhole", &warp_image_wide_to_pinhole_tensor);
    m.def("warp_image_equirectangular_to_pinhole", &warp_image_equirectangular_to_pinhole_tensor);
    m.def("warp_image_pinhole_to_wide", &warp_image_pinhole_to_wide_tensor);
    m.def("warp_linear_depth_pinhole_to_wide", &warp_linear_depth_pinhole_to_wide_tensor);
    m.def("warp_ray_depth_pinhole_to_wide", &warp_ray_depth_pinhole_to_wide_tensor);
    m.def("warp_points_pinhole_to_wide", &warp_points_pinhole_to_wide_tensor);
    m.def("warp_depth_pinhole_to_wide_scale_matrix", &warp_depth_pinhole_to_wide_scale_matrix_tensor);
    m.def("ppisp_forward", &ppisp_forward);
    m.def("ppisp_backward", &ppisp_backward);
    m.def("compute_ppsip_regularization_forward", &compute_ppsip_regularization_forward);
    m.def("compute_ppsip_regularization_backward", &compute_ppsip_regularization_backward);

    // FusedSSIM.cuh
    m.def("fused_ssim_forward", &fused_ssim_forward);
    m.def("fused_ssim_backward", &fused_ssim_backward);
    m.def("fused_ssim_inplace", &fused_ssim_inplace);

    #if 0
    // SplatTileIntersector.cuh
    m.def("intersect_splat_tile_3dgs", &intersect_splat_tile_3dgs);
    // m.def("intersect_splat_tile_opaque_triangle", &intersect_splat_tile_opaque_triangle);
    // m.def("intersect_splat_tile_voxel", &intersect_splat_tile_voxel);

    // SVHash.cuh
    m.def("svhash_create_initial_volume", &svhashCreateInitialVolume);
    m.def("svhash_get_voxels", &svhashGetVoxels);
    m.def("svhash_split_voxels", &svhashSplitVoxels);

    #endif

#if 0
    // PixelWise.cuh (partial)
    m.def("rendered_depth_to_expected_depth_forward", &rendered_depth_to_expected_depth_forward);
    m.def("rendered_depth_to_expected_depth_backward", &rendered_depth_to_expected_depth_backward);
    m.def("ray_depth_to_linear_depth_forward", &ray_depth_to_linear_depth_forward);
    m.def("ray_depth_to_linear_depth_backward", &ray_depth_to_linear_depth_backward);

    // ProjectionFwd.cuh
    m.def("projection_3dgs_forward", &projection_3dgs_forward);
    m.def("projection_mip_forward", &projection_mip_forward);
    m.def("projection_3dgut_forward", &projection_3dgut_forward);
    m.def("projection_3dgut_sv_forward", &projection_3dgut_forward);
    // m.def("projection_opaque_triangle_forward", &projection_opaque_triangle_forward);
    // m.def("projection_voxel_forward", &projection_voxel_forward);

    // ProjectionBwd.cuh
    m.def("projection_3dgs_backward", &projection_3dgs_backward);
    m.def("projection_mip_backward", &projection_mip_backward);
    m.def("projection_3dgut_backward", &projection_3dgut_backward);
    m.def("projection_3dgut_sv_backward", &projection_3dgut_backward);
    // m.def("projection_opaque_triangle_backward", &projection_opaque_triangle_backward);
    // m.def("projection_voxel_backward", &projection_voxel_backward);

    // ProjectionBwd.cuh, second order
    m.def("projection_3dgs_backward_with_position_hessian_diagonal", &projection_3dgs_backward_with_position_hessian_diagonal);
    m.def("projection_mip_backward_with_position_hessian_diagonal", &projection_mip_backward_with_position_hessian_diagonal);
    m.def("projection_3dgut_backward_with_position_hessian_diagonal", &projection_3dgut_backward_with_position_hessian_diagonal);
    m.def("projection_3dgs_backward_with_hessian_diagonal", &projection_3dgs_backward_with_hessian_diagonal);
    m.def("projection_mip_backward_with_hessian_diagonal", &projection_mip_backward_with_hessian_diagonal);
    m.def("projection_3dgut_backward_with_hessian_diagonal", &projection_3dgut_backward_with_hessian_diagonal);

    // ProjectionPackedFwd.cuh
    m.def("projection_3dgs_packed_forward", &projection_3dgs_packed_forward);
    m.def("projection_mip_packed_forward", &projection_mip_packed_forward);
    m.def("projection_3dgut_packed_forward", &projection_3dgut_packed_forward);
    m.def("projection_3dgut_sv_packed_forward", &projection_3dgut_packed_forward);
    // m.def("projection_opaque_triangle_packed_forward", &projection_opaque_triangle_packed_forward);
    // m.def("projection_voxel_packed_forward", &projection_voxel_packed_forward);
#endif

#if 0
    // ProjectionHeteroFwd.cuh
    m.def("projection_3dgs_hetero_forward", &projection_3dgs_hetero_forward);
    m.def("projection_mip_hetero_forward", &projection_mip_hetero_forward);
    m.def("projection_3dgut_hetero_forward", &projection_3dgut_hetero_forward);
    // m.def("projection_opaque_triangle_hetero_forward", &projection_opaque_triangle_hetero_forward);

    // ProjectionHeteroBwd.cuh
    m.def("projection_3dgs_hetero_backward", &projection_3dgs_hetero_backward);
    m.def("projection_mip_hetero_backward", &projection_mip_hetero_backward);
    m.def("projection_3dgut_hetero_backward", &projection_3dgut_hetero_backward);
    // m.def("projection_opaque_triangle_hetero_backward", &projection_opaque_triangle_hetero_backward);
#endif

#if 0
    // RasterizationFwd.cuh and RasterizationBwd.cuh
    m.def("rasterization_3dgs_forward", &rasterize_to_pixels_3dgs_fwd);
    m.def("rasterization_3dgs_backward", &rasterize_to_pixels_3dgs_bwd);
    m.def("rasterization_3dgs_backward_with_hessian_diagonal", &rasterize_to_pixels_3dgs_bwd_with_hessian_diagonal);
    m.def("rasterization_mip_forward", &rasterize_to_pixels_mip_fwd);
    m.def("rasterization_mip_backward", &rasterize_to_pixels_3dgs_bwd);
    m.def("rasterization_mip_backward_with_hessian_diagonal", &rasterize_to_pixels_3dgs_bwd_with_hessian_diagonal);

    // RasterizationEval3DFwd.cuh and RasterizationEval3DBwd.cuh
    m.def("rasterization_3dgut_forward", &rasterize_to_pixels_3dgut_fwd);
    m.def("rasterization_3dgut_backward", &rasterize_to_pixels_3dgut_bwd);
    m.def("rasterization_3dgut_backward_with_hessian_diagonal", &rasterize_to_pixels_3dgut_bwd_with_hessian_diagonal);
    m.def("rasterization_3dgut_sv_forward", &rasterize_to_pixels_3dgut_fwd);
    m.def("rasterization_3dgut_sv_backward", &rasterize_to_pixels_3dgut_bwd);
    // m.def("rasterization_opaque_triangle_forward", &rasterize_to_pixels_opaque_triangle_sorted_fwd);
    // m.def("rasterization_opaque_triangle_backward", &rasterize_to_pixels_opaque_triangle_sorted_bwd);
    // m.def("rasterization_voxel_forward", &rasterize_to_pixels_voxel_eval3d_fwd);
    // m.def("rasterization_voxel_backward", &rasterize_to_pixels_voxel_eval3d_bwd);
#endif

#if 0
    // Optimizer.cuh
    m.def("set_zero", &set_zero_tensor);
    m.def("zeros_like", &zeros_like_tensor);
    m.def("fused_adam", &fused_adam);
    m.def("offloaded_adam", &offloaded_adam);
    m.def("fused_adam_multi", &fused_adam_multi);
    m.def("fused_adam_riemannian_quat", &fused_adam_riemannian_quat);
    m.def("fused_newton", &fused_newton);
    m.def("fused_newton_multi", &fused_newton_multi);
    m.def("fused_adam_scale_agnostic_mean", &fused_adam_scale_agnostic_mean);
    m.def("fused_optim_3dgs_geometry", &fused_optim_3dgs_geometry);
    m.def("fused_adam_with_steps", &fused_adam_with_steps_tensor);
    m.def("fused_adam_with_steps_8bit", &fused_adam_with_steps_8bit_tensor);
    m.def("fused_3dgs2tr_mean_optim", &fused_3dgs2tr_mean_optim);
    m.def("fused_3dgs2tr_scale_optim", &fused_3dgs2tr_scale_optim);
    m.def("fused_3dgs2tr_color_optim", &fused_3dgs2tr_color_optim);
    m.def("fused_3dgs2tr_opacity_optim", &fused_3dgs2tr_opacity_optim);
    m.def("fused_3dgs2tr_quat_optim", &fused_3dgs2tr_quat_optim);
    m.def("fused_adam_linear_rgb_optim", &fused_adam_linear_rgb_optim);
    m.def("fused_adamtr_rgb_optim", &fused_adamtr_rgb_optim);
    m.def("fused_adamtr_linear_rgb_optim", &fused_adamtr_linear_rgb_optim);
    m.def("fused_adamtr_rgb_sh_optim", &fused_adamtr_rgb_sh_optim);
    m.def("fused_adamtr_linear_rgb_sh_optim", &fused_adamtr_linear_rgb_sh_optim);
#endif

    // FusedProjectionBwdOptim.cuh
    m.def("fused_projection_bwd_optimizer_3dgs",  &fused_projection_bwd_optimizer_3dgs);
    m.def("fused_projection_bwd_optimizer_mip",   &fused_projection_bwd_optimizer_mip);
    m.def("fused_projection_bwd_optimizer_3dgut", &fused_projection_bwd_optimizer_3dgut);

    // Visualizer.cuh
    m.def("blit_train_cameras", &blit_train_cameras_tensor);

    // Engine.h - config structs (built on the Python side, passed to engine_*_step).
    py::class_<LossConfig>(m, "LossConfig")
        .def(py::init<>())
        .def_readwrite("weights",          &LossConfig::weights)
        .def_readwrite("w_ssim",           &LossConfig::w_ssim)
        .def_readwrite("num_loss_scales",  &LossConfig::num_loss_scales)
        .def_readwrite("compute_loss_map", &LossConfig::compute_loss_map)
        .def_readwrite("structure_only_loss_map", &LossConfig::structure_only_loss_map)
        .def_readwrite("overexposure_reg_weight", &LossConfig::overexposure_reg_weight);

    py::class_<OptimConfig>(m, "OptimConfig")
        .def(py::init<>())
        .def_readwrite("lr_means",        &OptimConfig::lr_means)
        .def_readwrite("lr_quats",        &OptimConfig::lr_quats)
        .def_readwrite("lr_scales",       &OptimConfig::lr_scales)
        .def_readwrite("lr_opacities",    &OptimConfig::lr_opacities)
        .def_readwrite("lr_features_dc",  &OptimConfig::lr_features_dc)
        .def_readwrite("lr_features_sh",  &OptimConfig::lr_features_sh)
        .def_readwrite("max_gauss_ratio",                &OptimConfig::max_gauss_ratio)
        .def_readwrite("scale_regularization_weight",    &OptimConfig::scale_regularization_weight)
        .def_readwrite("mcmc_opacity_reg_weight",        &OptimConfig::mcmc_opacity_reg_weight)
        .def_readwrite("mcmc_scale_reg_weight",          &OptimConfig::mcmc_scale_reg_weight)
        .def_readwrite("erank_reg_weight",               &OptimConfig::erank_reg_weight)
        .def_readwrite("erank_reg_weight_s3",            &OptimConfig::erank_reg_weight_s3)
        .def_readwrite("quat_norm_reg_weight",           &OptimConfig::quat_norm_reg_weight)
        .def_readwrite("sh_reg_weight",                  &OptimConfig::sh_reg_weight)
        .def_readwrite("use_scale_agnostic_mean",        &OptimConfig::use_scale_agnostic_mean)
        .def_readwrite("quantize_sh",                    &OptimConfig::quantize_sh)
        .def_readwrite("use_per_splat_bias_correction",  &OptimConfig::use_per_splat_bias_correction)
        .def_readwrite("use_fused_proj_bwd_optim",       &OptimConfig::use_fused_proj_bwd_optim)
        .def_readwrite("use_color_trust_region",         &OptimConfig::use_color_trust_region)
        .def_readwrite("color_is_linear",                &OptimConfig::color_is_linear)
        .def_readwrite("eps_tr",                         &OptimConfig::eps_tr);

    py::class_<DensifyConfig>(m, "DensifyConfig")
        .def(py::init<>())
        .def_readwrite("refine_start_iter",              &DensifyConfig::refine_start_iter)
        .def_readwrite("refine_stop_num_iter",           &DensifyConfig::refine_stop_num_iter)
        .def_readwrite("refine_every",                   &DensifyConfig::refine_every)
        .def_readwrite("growth_factor",                  &DensifyConfig::growth_factor)
        .def_readwrite("min_opacity",                    &DensifyConfig::min_opacity)
        .def_readwrite("max_screen_size",                &DensifyConfig::max_screen_size)
        .def_readwrite("max_screen_size_clip_hardness",  &DensifyConfig::max_screen_size_clip_hardness)
        .def_readwrite("max_world_size",                 &DensifyConfig::max_world_size)
        .def_readwrite("noise_lr",                       &DensifyConfig::noise_lr)
        .def_readwrite("noise_lr_final",                 &DensifyConfig::noise_lr_final)
        .def_readwrite("use_revised_densification",      &DensifyConfig::use_revised_densification)
        .def_readwrite("score_mode",                     &DensifyConfig::score_mode);

    py::class_<BilagridStepConfig>(m, "BilagridStepConfig")
        .def(py::init<>())
        .def_readwrite("lr_rgb",         &BilagridStepConfig::lr_rgb)
        .def_readwrite("lr_depth",       &BilagridStepConfig::lr_depth)
        .def_readwrite("lr_normal",      &BilagridStepConfig::lr_normal)
        .def_readwrite("tv_weight_rgb",  &BilagridStepConfig::tv_weight_rgb)
        .def_readwrite("tv_weight_depth",  &BilagridStepConfig::tv_weight_depth)
        .def_readwrite("tv_weight_normal", &BilagridStepConfig::tv_weight_normal);

    py::class_<PpispStepConfig>(m, "PpispStepConfig")
        .def(py::init<>())
        .def_readwrite("lr",          &PpispStepConfig::lr)
        .def_readwrite("reg_weights", &PpispStepConfig::reg_weights);

    py::class_<BackgroundStepConfig>(m, "BackgroundStepConfig")
        .def(py::init<>())
        .def_readwrite("lr_dc",            &BackgroundStepConfig::lr_dc)
        .def_readwrite("lr_sh",            &BackgroundStepConfig::lr_sh)
        .def_readwrite("randomize_weight", &BackgroundStepConfig::randomize_weight)
        .def_readwrite("seed",             &BackgroundStepConfig::seed);

    py::class_<EngineStepConfig>(m, "EngineStepConfig")
        .def(py::init<>())
        .def_readwrite("loss",       &EngineStepConfig::loss)
        .def_readwrite("optim",      &EngineStepConfig::optim)
        .def_readwrite("densify",    &EngineStepConfig::densify)
        .def_readwrite("bilagrid",   &EngineStepConfig::bilagrid)
        .def_readwrite("ppisp",      &EngineStepConfig::ppisp)
        .def_readwrite("background", &EngineStepConfig::background);

    // Engine.h - unified forward/backward/optimize
    m.def("engine_reset", &engine_reset);
    m.def("set_data_3dgs", &set_data_3dgs);
    m.def("set_camera_params", &set_camera_params);
    m.def("set_training_data", &set_training_data);
    m.def("engine_forward_3dgs", &forward_3dgs);
    m.def("engine_compute_loss_backward", &engine_compute_loss_backward);
    m.def("engine_optim_step", &engine_optim_step);
    m.def("engine_densify_step", &engine_densify_step);
    m.def("engine_train_step", &engine_train_step);
    m.def("engine_debug_forward", &engine_debug_forward);
    m.def("engine_copy_accum_buffer", &engine_copy_accum_buffer);
    m.def("engine_get_cur_num_splats", &engine_get_cur_num_splats);
    m.def("engine_copy_render_to_host", &engine_copy_render_to_host);
    m.def("engine_copy_splats_to_host", &engine_copy_splats_to_host);
    m.def("engine_get_pool_breakdown", &engine_get_pool_breakdown);
    m.def("engine_get_scratch_bytes", &engine_get_scratch_bytes);
    m.def("engine_init_bilagrid_rgb", &engine_init_bilagrid_rgb);
    m.def("engine_init_bilagrid_depth", &engine_init_bilagrid_depth);
    m.def("engine_init_bilagrid_normal", &engine_init_bilagrid_normal);
    m.def("engine_bilagrid_forward", &engine_bilagrid_forward);
    m.def("engine_bilagrid_optim_step", &engine_bilagrid_optim_step);
    m.def("engine_init_ppisp", &engine_init_ppisp);
    m.def("engine_ppisp_forward", &engine_ppisp_forward);
    m.def("engine_ppisp_optim_step", &engine_ppisp_optim_step);
    m.def("engine_init_background_noise", &engine_init_background_noise);
    m.def("engine_init_background_sh",    &engine_init_background_sh);
    m.def("engine_set_background_step_params", &engine_set_background_step_params);
    m.def("engine_background_optim_step", &engine_background_optim_step);
    m.def("engine_copy_background_to_host", &engine_copy_background_to_host);
    m.def("engine_init_color_space", &engine_init_color_space);
    m.def("engine_save_checkpoint", &engine_save_checkpoint);

#if 0
    // Densify.cuh
    m.def("quantile_of_abs_of_finite_elements", &quantile_of_abs_of_finite_elements_tensor);
    m.def("normalize_by_median_inplace", &normalize_by_median_inplace_tensor);
    m.def("inplace_index", &inplace_index_tensor);
    m.def("inplace_scatter_add", &inplace_scatter_add_tensor);
    m.def("inplace_scatter_max", &inplace_scatter_max_tensor);
    m.def("weighted_sample_without_replacement", &weighted_sample_without_replacement_tensor);
    m.def("densify_clip_scale", &densify_clip_scale_tensor);
    m.def("densify_update_weight", &densify_update_weight);
    m.def("relocate_splats_with_long_axis_split", &relocate_splats_with_long_axis_split_tensor);
    m.def("add_splats_with_long_axis_split", &add_splats_with_long_axis_split_tensor);
    m.def("relocate_splats_mcmc", &relocate_splats_mcmc_tensor);
    m.def("add_splats_mcmc", &add_splats_mcmc_tensor);
    m.def("cov_scale_init", &cov_scale_init_tensor);
    m.def("mcmc_add_noise", &mcmc_add_noise_tensor);
    m.def("revised_add_noise", &revised_add_noise_tensor);
    // m.def("compute_relocation", &compute_relocation_tensor);
    m.def("long_axis_split", &long_axis_split_tensor);
    m.def("laplacian_edge_filter", &laplacian_edge_filter_tensor);
    m.def("smoothed_laplacian_edge_filter", &smoothed_laplacian_edge_filter_tensor);
    m.def("canny_edge_filter", &canny_edge_filter_tensor);

    // BilagridUtils.cuh
    m.def("dct3d_type1_ortho", &dct3d_type1_ortho_tensor);
#endif

    // Visualizer.cuh
    m.def("blit_train_cameras", &blit_train_cameras_tensor);
}
