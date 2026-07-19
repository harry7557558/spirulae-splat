#define SLANG_PRELUDE_EXPORT
#include "Camera.h"
#include "Common.cuh"
#include "Engine.h"

#include "IntersectTile.cuh"
#include "BackgroundSphericalHarmonics.cuh"
#include "PerSplatLoss.cuh"
#include "PerPixelLoss.cuh"
#include "PixelWise.cuh"
#include "FusedSSIM.cuh"
#include "SplatTileIntersector.cuh"
// #include "Projection.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionBwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"
#include "Optimizer.cuh"
#include "FusedProjectionBwdOptim.cuh"
#include "Densify.cuh"
#include "BilagridUtils.cuh"
#include "Visualizer.cuh"

#define TORCH_INDUCTOR_CPP_WRAPPER
#include <torch/extension.h>

#include "Delaunay3D.h"
#include "Meshing.h"
#include <cstring>

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

// Delaunay3D.h — self-contained multithreaded 3D Delaunay (CPU).
// points: (N,3) tensor (any float/double dtype, moved to CPU double internally).
// Returns {cell_vertices (M,4) int32, cell_adjacents (M,4) int32 (-1 on hull)}.
static std::vector<at::Tensor> delaunay3d_tensor(
    at::Tensor points, int64_t num_threads
) {
    TORCH_CHECK(points.dim() == 2 && points.size(1) == 3,
                "delaunay3d: points must have shape (N,3)");
    at::Tensor pts = points.to(at::kDouble).contiguous().cpu();
    int n = (int)pts.size(0);
    auto res = delaunay3d::compute_delaunay_3d(
        pts.data_ptr<double>(), n, (int)num_threads, /*adjacency=*/true
    );
    auto opts = at::TensorOptions().dtype(at::kInt);
    at::Tensor cells = at::empty({(int64_t)res.nb_cells, 4}, opts);
    at::Tensor adj   = at::empty({(int64_t)res.nb_cells, 4}, opts);
    if(res.nb_cells > 0) {
        std::memcpy(cells.data_ptr<int>(), res.cell_vertices.data(),
                    sizeof(int) * res.cell_vertices.size());
        std::memcpy(adj.data_ptr<int>(), res.cell_adjacents.data(),
                    sizeof(int) * res.cell_adjacents.size());
    }
    return {cells, adj};
}

// Meshing.h -- surface mesh extraction from a trained 3DGS model.
// Splat params are CPU float32 tensors in standard inria-PLY convention
// (raw quats, log scales, logit opacities). camera_positions may be empty
// (no dataset -> static density occupancy). Writes the mesh in one or more of
// PLY/OBJ/GLTF/GLB (comma-separated `formats`), colored per `color_mode`
// (none / vertex / texture; incompatible format+mode pairs raise).
static bool generate_mesh_tensor(
    at::Tensor means, at::Tensor quats, at::Tensor scales, at::Tensor opacities,
    at::Tensor features_dc, at::Tensor camera_positions,
    const std::string& output_path,
    double iso, double merge_factor, int64_t bisection_iters,
    int64_t max_cameras, int64_t max_grid_res, double grid_cell_factor,
    int64_t num_threads, bool verbose,
    at::Tensor viewmats, at::Tensor intrins, at::Tensor dist_coeffs,
    at::Tensor cam_widths, at::Tensor cam_heights, const std::string& camera_model,
    int64_t carve_k, bool cull_unseen,
    double merge_max_flip_deg, int64_t floater_min_faces, double fill_hole_ratio,
    int64_t fill_hole_max_edges, double degenerate_angle_deg,
    const std::string& color_mode, const std::string& formats,
    int64_t texture_size, int64_t tex_gutter_px, double chart_angle_deg,
    int64_t quality_iters
) {
    auto f32 = [](at::Tensor t) { return t.to(at::kFloat).contiguous().cpu(); };
    at::Tensor m = f32(means), q = f32(quats), s = f32(scales), o = f32(opacities);
    at::Tensor fdc = f32(features_dc);
    int N = (int)m.size(0);
    TORCH_CHECK(m.dim() == 2 && m.size(1) == 3, "means must be (N,3)");
    TORCH_CHECK(q.size(0) == N && q.size(1) == 4, "quats must be (N,4)");
    TORCH_CHECK(s.size(0) == N && s.size(1) == 3, "scales must be (N,3)");
    TORCH_CHECK(o.numel() == N, "opacities must be (N,)");
    TORCH_CHECK(fdc.dim() == 2 && fdc.size(0) == N && fdc.size(1) == 3,
                "features_dc must be (N,3)");

    int C = 0;
    const float* cam_ptr = nullptr;
    at::Tensor cam;
    if (camera_positions.numel() > 0) {
        cam = f32(camera_positions);
        TORCH_CHECK(cam.dim() == 2 && cam.size(1) == 3, "camera_positions must be (C,3)");
        C = (int)cam.size(0);
        cam_ptr = cam.data_ptr<float>();
    }

    // Optional camera intrinsics for the rasterize-and-sample path. Kept alive
    // in this scope (vm/in/dc) so the pointers passed in CameraParams stay valid
    // for the whole generate_mesh call.
    meshing::CameraParams cams;
    at::Tensor vm, in, dc, cw, ch;
    if (C > 0 && viewmats.numel() > 0 && intrins.numel() > 0) {
        vm = f32(viewmats);
        in = f32(intrins);
        TORCH_CHECK(vm.dim() == 3 && vm.size(0) == C && vm.size(1) == 4 && vm.size(2) == 4,
                    "viewmats must be (C,4,4) matching camera_positions");
        TORCH_CHECK(in.dim() == 2 && in.size(0) == C && in.size(1) == 4,
                    "intrins must be (C,4): fx, fy, cx, cy");
        cams.viewmats = vm.data_ptr<float>();
        cams.intrins  = in.data_ptr<float>();
        if (dist_coeffs.numel() > 0) {
            dc = f32(dist_coeffs);
            TORCH_CHECK(dc.dim() == 2 && dc.size(0) == C && dc.size(1) == 10,
                        "dist_coeffs must be (C,10) in engine layout");
            cams.dist_coeffs = dc.data_ptr<float>();
        }
        // Per-camera image size. Accept a scalar (broadcast to all cameras) or a
        // (C,) tensor so datasets with mixed resolutions are supported.
        auto i32 = [&](at::Tensor t) {
            t = t.to(at::kInt).contiguous().cpu();
            if (t.numel() == 1) t = t.expand({C}).contiguous();
            TORCH_CHECK(t.numel() == C, "cam_widths/cam_heights must be scalar or (C,)");
            return t;
        };
        cw = i32(cam_widths);
        ch = i32(cam_heights);
        cams.widths  = cw.data_ptr<int>();
        cams.heights = ch.data_ptr<int>();
        cams.camera_model = camera_model;
    }

    meshing::MeshingConfig cfg;
    cfg.iso = (float)iso;
    cfg.merge_factor = (float)merge_factor;
    cfg.bisection_iters = (int)bisection_iters;
    cfg.max_cameras = (int)max_cameras;
    cfg.max_grid_res = (int)max_grid_res;
    cfg.grid_cell_factor = (float)grid_cell_factor;
    cfg.num_threads = (int)num_threads;
    cfg.verbose = verbose;
    cfg.carve_k = (int)carve_k;
    cfg.cull_unseen = cull_unseen;
    cfg.merge_max_flip_deg = (float)merge_max_flip_deg;
    cfg.floater_min_faces = (int)floater_min_faces;
    cfg.fill_hole_ratio = (float)fill_hole_ratio;
    cfg.fill_hole_max_edges = (int)fill_hole_max_edges;
    cfg.degenerate_angle_deg = (float)degenerate_angle_deg;

    if (color_mode == "none")         cfg.color_mode = meshing::MeshColorMode::None;
    else if (color_mode == "vertex")  cfg.color_mode = meshing::MeshColorMode::Vertex;
    else if (color_mode == "texture") cfg.color_mode = meshing::MeshColorMode::Texture;
    else TORCH_CHECK(false, "color_mode must be none/vertex/texture, got ", color_mode);
    // throws on unknown formats / bad "+jpeg<q>" encodings / duplicates
    std::vector<meshing::MeshFormatSpec> specs = meshing::parse_mesh_formats(formats);
    TORCH_CHECK(!specs.empty(), "formats must name at least one of ply/obj/gltf/glb");
    cfg.formats.clear();
    for (const auto& spec : specs) {
        std::string err = meshing::check_export_support(spec, cfg.color_mode);
        TORCH_CHECK(err.empty(), err);
        cfg.formats.push_back(spec.token());
    }
    cfg.texture_size = (int)texture_size;
    cfg.tex_gutter_px = (int)tex_gutter_px;
    cfg.chart_angle_deg = (float)chart_angle_deg;
    cfg.quality_iters = (int)quality_iters;

    return meshing::generate_mesh(
        m.data_ptr<float>(), q.data_ptr<float>(), s.data_ptr<float>(), o.data_ptr<float>(),
        fdc.data_ptr<float>(), N, cam_ptr, C, cams, cfg, output_path);
}

PYBIND11_MODULE(TORCH_EXTENSION_NAME, m) {

    // Meshing.h
    m.def("generate_mesh", &generate_mesh_tensor,
          "Extract a surface mesh from a trained 3DGS model and save to PLY",
          pybind11::arg("means"), pybind11::arg("quats"),
          pybind11::arg("scales"), pybind11::arg("opacities"),
          pybind11::arg("features_dc"), pybind11::arg("camera_positions"),
          pybind11::arg("output_path"),
          pybind11::arg("iso") = 0.5,
          pybind11::arg("merge_factor") = 1.0,
          pybind11::arg("bisection_iters") = 3,
          pybind11::arg("max_cameras") = 64,
          pybind11::arg("max_grid_res") = 512,
          pybind11::arg("grid_cell_factor") = 2.0,
          pybind11::arg("num_threads") = 0,
          pybind11::arg("verbose") = true,
          pybind11::arg("viewmats") = at::Tensor(),
          pybind11::arg("intrins") = at::Tensor(),
          pybind11::arg("dist_coeffs") = at::Tensor(),
          pybind11::arg("cam_widths") = at::Tensor(),
          pybind11::arg("cam_heights") = at::Tensor(),
          pybind11::arg("camera_model") = std::string("PINHOLE"),
          pybind11::arg("carve_k") = 1,
          pybind11::arg("cull_unseen") = true,
          pybind11::arg("merge_max_flip_deg") = 60.0,
          pybind11::arg("floater_min_faces") = 100,
          pybind11::arg("fill_hole_ratio") = 0.05,
          pybind11::arg("fill_hole_max_edges") = 12,
          pybind11::arg("degenerate_angle_deg") = 2.0,
          pybind11::arg("color_mode") = std::string("vertex"),
          pybind11::arg("formats") = std::string("ply"),
          pybind11::arg("texture_size") = 0,  // 0 = auto
          pybind11::arg("tex_gutter_px") = 4,
          pybind11::arg("chart_angle_deg") = 60.0,
          pybind11::arg("quality_iters") = 3);

    // Delaunay3D.h
    m.def("delaunay3d", &delaunay3d_tensor,
          "3D Delaunay triangulation (CPU, multithreaded)",
          pybind11::arg("points"), pybind11::arg("num_threads") = 0);

    // IntersectTile.cuh
    m.def("intersect_tile", &do_intersect_tile_generic);
    // m.def("intersect_tile_3dgs", &intersect_tile_3dgs_tensor);
    // m.def("intersect_tile_mip", &intersect_tile_3dgs_tensor);
    // m.def("intersect_tile_3dgut", &intersect_tile_3dgut_tensor);

    // BackgroundSphericalHarmonics.cuh — TorchTensorView API. Outputs are
    // caller-allocated; the C++ side fills them in place.
    m.def("render_background_sh_forward",  &render_background_sh_forward);
    m.def("render_background_sh_backward", &render_background_sh_backward);

    #if 0
    // PerSplatLoss.cuh
    m.def("compute_per_splat_losses_forward", &compute_per_splat_losses_forward);
    m.def("compute_per_splat_losses_backward", &compute_per_splat_losses_backward);
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
    m.def("engine_viewer_init", &engine_viewer_init);
    m.def("engine_viewer_set_grid", &engine_viewer_set_grid);
    m.def("engine_blit_view",   &engine_blit_view);

    // Engine.h - config structs (built on the Python side, passed to engine_*_step).
    py::class_<LossConfig>(m, "LossConfig")
        .def(py::init<>())
        .def_readwrite("weights",          &LossConfig::weights)
        .def_readwrite("w_ssim",           &LossConfig::w_ssim)
        .def_readwrite("num_loss_scales",  &LossConfig::num_loss_scales)
        .def_readwrite("loss_scale_min_pixels", &LossConfig::loss_scale_min_pixels)
        .def_readwrite("compute_loss_map", &LossConfig::compute_loss_map)
        .def_readwrite("loss_map_mode", &LossConfig::loss_map_mode)
        .def_readwrite("robust_edge_aware_quantile", &LossConfig::robust_edge_aware_quantile)
        .def_readwrite("overexposure_reg_weight", &LossConfig::overexposure_reg_weight)
        .def_readwrite("color_shift_reg_weight",  &LossConfig::color_shift_reg_weight)
        .def_readwrite("color_shift_reg_beta",    &LossConfig::color_shift_reg_beta)
        .def_readwrite("input_depth_is_ray_depth", &LossConfig::input_depth_is_ray_depth);

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
        .def_readwrite("sh_optim_bits",                  &OptimConfig::sh_optim_bits)
        .def_readwrite("sh_value_bits",                  &OptimConfig::sh_value_bits)
        .def_readwrite("non_sh_optim_bits",              &OptimConfig::non_sh_optim_bits)
        .def_readwrite("quantization_level",          &OptimConfig::quantization_level)
        .def_readwrite("use_per_splat_bias_correction",  &OptimConfig::use_per_splat_bias_correction)
        .def_readwrite("use_fused_proj_bwd_optim",       &OptimConfig::use_fused_proj_bwd_optim)
        .def_readwrite("write_densify_world_grad_score", &OptimConfig::write_densify_world_grad_score)
        .def_readwrite("split_batch",     &OptimConfig::split_batch)
        .def_readwrite("use_color_trust_region",         &OptimConfig::use_color_trust_region)
        .def_readwrite("color_is_linear",                &OptimConfig::color_is_linear)
        .def_readwrite("eps_tr",                         &OptimConfig::eps_tr);

    py::class_<DensifyConfig>(m, "DensifyConfig")
        .def(py::init<>())
        .def_readwrite("refine_start_iter",              &DensifyConfig::refine_start_iter)
        .def_readwrite("refine_stop_num_iter",           &DensifyConfig::refine_stop_num_iter)
        .def_readwrite("refine_stop_iter",               &DensifyConfig::refine_stop_iter)
        .def_readwrite("refine_every",                   &DensifyConfig::refine_every)
        .def_readwrite("growth_factor",                  &DensifyConfig::growth_factor)
        .def_readwrite("min_opacity",                    &DensifyConfig::min_opacity)
        .def_readwrite("max_screen_size",                &DensifyConfig::max_screen_size)
        .def_readwrite("max_screen_size_clip_hardness",  &DensifyConfig::max_screen_size_clip_hardness)
        .def_readwrite("max_world_size",                 &DensifyConfig::max_world_size)
        .def_readwrite("noise_lr",                       &DensifyConfig::noise_lr)
        .def_readwrite("noise_lr_final",                 &DensifyConfig::noise_lr_final)
        .def_readwrite("use_revised_densification",      &DensifyConfig::use_revised_densification)
        .def_readwrite("score_mode",                     &DensifyConfig::score_mode)
        .def_readwrite("score_blend_world_grad",         &DensifyConfig::score_blend_world_grad)
        .def_readwrite("las_split_opacity_k_init",       &DensifyConfig::las_split_opacity_k_init)
        .def_readwrite("las_split_opacity_k_final",      &DensifyConfig::las_split_opacity_k_final)
        .def_readwrite("las_split_opacity_k_warmup",     &DensifyConfig::las_split_opacity_k_warmup);

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
        .def_readwrite("lr",                  &PpispStepConfig::lr)
        .def_readwrite("reg_weights",         &PpispStepConfig::reg_weights)
        .def_readwrite("run_before_bilagrid", &PpispStepConfig::run_before_bilagrid);

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
    m.def("set_training_data", &set_training_data,
          pybind11::arg("gt_rgb"), pybind11::arg("gt_depth"),
          pybind11::arg("gt_normal"), pybind11::arg("gt_alpha"),
          pybind11::arg("input_depth_is_ray_depth") = true);
    m.def("engine_forward_3dgs", &forward_3dgs);
    m.def("engine_compute_loss_backward", &engine_compute_loss_backward);
    m.def("engine_backward_from_render_grad", &engine_backward_from_render_grad);
    m.def("engine_optim_step", &engine_optim_step);
    m.def("engine_densify_step", &engine_densify_step);
    m.def("engine_train_step", &engine_train_step);

    // ---- DataManager + managed train step ---------------------------------
    py::enum_<CacheMode>(m, "DataManagerCacheMode")
        .value("CPU",  CacheMode::CPU)
        .value("DISK", CacheMode::DISK)
        .export_values();

    py::class_<DataManagerConfig>(m, "DataManagerConfig")
        .def(py::init<>())
        .def_readwrite("cache_mode",       &DataManagerConfig::cache_mode)
        .def_readwrite("load_masks",       &DataManagerConfig::load_masks)
        .def_readwrite("load_depths",      &DataManagerConfig::load_depths)
        .def_readwrite("load_normals",     &DataManagerConfig::load_normals)
        .def_readwrite("train_batch_size", &DataManagerConfig::train_batch_size)
        .def_readwrite("val_batch_size",   &DataManagerConfig::val_batch_size)
        .def_readwrite("seed",             &DataManagerConfig::seed)
        .def_readwrite("workers_rgb",      &DataManagerConfig::workers_rgb)
        .def_readwrite("workers_depth",    &DataManagerConfig::workers_depth)
        .def_readwrite("workers_normal",   &DataManagerConfig::workers_normal)
        .def_readwrite("prefetch_batches", &DataManagerConfig::prefetch_batches)
        .def_readwrite("mask_boundary_offset",
                                           &DataManagerConfig::mask_boundary_offset)
        .def_readwrite("warp_to_pinhole",  &DataManagerConfig::warp_to_pinhole);

    py::enum_<CameraModelType>(m, "CameraModelType")
        .value("PINHOLE",         CameraModelType::PINHOLE)
        .value("FISHEYE",         CameraModelType::FISHEYE)
        .value("EQUISOLID",       CameraModelType::EQUISOLID)
        .value("EQUIRECTANGULAR", CameraModelType::EQUIRECTANGULAR)
        .export_values();

    m.def("engine_setup_data_manager",   &engine_setup_data_manager);
    // Release the GIL while engine_train_step_managed runs. This is the
    // single longest-blocking C++ call in the hot path (~200 ms/step incl.
    // the synchronous reads needed to populate loss_dict). Holding the GIL
    // here starves the HTTP/render threads, turning a 250 ms render request
    // into a multi-second wait. None of the C++ body touches Python objects
    // (cfg is converted before entry; the returned std::map is converted
    // after the GIL is reacquired), so releasing it is safe.
    m.def("engine_train_step_managed",   &engine_train_step_managed,
          pybind11::call_guard<pybind11::gil_scoped_release>());

    // COLMAP / NerfStudio camera-model string -> enum, mirroring
    // dataparser.py's lookup table. Exposed so Trainer can build the
    // enum without re-implementing the table on the Python side.
    m.def("camera_model_from_name", &camera_model_from_name);
    m.def("engine_debug_forward", &engine_debug_forward);
    m.def("engine_copy_accum_buffer", &engine_copy_accum_buffer);
    m.def("engine_get_cur_num_splats", &engine_get_cur_num_splats);
    m.def("engine_copy_render_to_host", &engine_copy_render_to_host);
    m.def("engine_copy_distortion_to_host", &engine_copy_distortion_to_host);
    m.def("engine_get_gt_rgb_shape",     &engine_get_gt_rgb_shape);
    m.def("engine_get_gt_alpha_shape",   &engine_get_gt_alpha_shape);
    m.def("engine_get_render_rgb_shape", &engine_get_render_rgb_shape);
    m.def("engine_copy_gt_rgb_to_host",  &engine_copy_gt_rgb_to_host);
    m.def("engine_copy_gt_alpha_to_host",&engine_copy_gt_alpha_to_host);
    m.def("engine_copy_splats_to_host", &engine_copy_splats_to_host);
    m.def("engine_copy_grads_to_host", &engine_copy_grads_to_host);
    m.def("engine_get_pool_breakdown", &engine_get_pool_breakdown);
    m.def("engine_get_pool_breakdown_categorized", &engine_get_pool_breakdown_categorized);
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
    m.def("engine_load_checkpoint", &engine_load_checkpoint);

#if 0
    // Densify.cuh
    m.def("quantile_of_abs_of_finite_elements", &quantile_of_abs_of_finite_elements_tensor);
    m.def("normalize_by_median_inplace", &normalize_by_median_inplace_tensor);
    m.def("inplace_index", &inplace_index_tensor);
    m.def("inplace_scatter_add", &inplace_scatter_add_tensor);
    m.def("inplace_scatter_max", &inplace_scatter_max_tensor);
    m.def("weighted_sample_without_replacement", &weighted_sample_without_replacement_tensor);
    m.def(PoolSlot::DensifyClipScale, &densify_clip_scale_tensor);
    m.def(PoolSlot::DensifyUpdateWeight, &densify_update_weight);
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
    m.def("engine_viewer_init", &engine_viewer_init);
    m.def("engine_viewer_set_grid", &engine_viewer_set_grid);
    m.def("engine_blit_view",   &engine_blit_view);
}
