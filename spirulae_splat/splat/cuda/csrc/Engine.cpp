#include "Engine.h"
#include "Tensor.h"
#include "Primitive.cuh"
#include "ProjectionFwd.cuh"
#include "ProjectionPackedFwd.cuh"
#include "ProjectionBwd.cuh"
#include "IntersectTile.cuh"
#include "RasterizationFwd.cuh"
#include "RasterizationBwd.cuh"
#include "RasterizationEval3DFwd.cuh"
#include "RasterizationEval3DBwd.cuh"


namespace Buffers {

std::string primitive;
int64_t num_splats;
int sh_degree;

TorchTensorView world_means;
TorchTensorView world_quats;
TorchTensorView world_scales;
TorchTensorView world_opacities;
TorchTensorView world_features_dc;
TorchTensorView world_features_sh;

int32_t num_cameras;
int32_t image_width;
int32_t image_height;
ssplat::CameraModelType camera_model;
std::string camera_model_str;
TorchTensorView camera_viewmats;
TorchTensorView camera_intrins;
TorchTensorView camera_dist_coeffs;

bool packed;

// Forward intermediate results (kept for backward pass)
std::optional<DeviceVector<int32_t>> fwd_camera_ids;
std::optional<DeviceVector<int32_t>> fwd_gaussian_ids;
DeviceTensor2D<float4> fwd_aabb;   // [nnz,1] packed or [C,N] non-packed
std::vector<DeviceTensorFloatND> fwd_splats_w;
std::vector<DeviceTensorFloatND> fwd_splats_s;
DeviceTensor3D<int32_t> fwd_tile_offsets;
DeviceVector<int32_t> fwd_flatten_ids;
DeviceTensor3D<float> fwd_render_Ts;
DeviceTensor3D<int32_t> fwd_last_ids;
RenderOutput::TensorTuple fwd_renders;  // for 3dgut backward

} // namespace Buffers


// Create DeviceTensorFloatND from TorchTensorView by appending trailing 1 (for float-typed access)
static DeviceTensorFloatND tv_to_fnd(const TorchTensorView& tv) {
    auto shape = std::get<2>(tv);
    shape.push_back(1LL);
    return DeviceTensorFloatND(TorchTensorView{std::get<0>(tv), std::get<1>(tv), shape});
}

// Create DeviceTensor2D<float4> view from DeviceVector<float4> (for packed projection backward)
static DeviceTensor2D<float4> vec_to_2d_float4(const DeviceVector<float4>& vec) {
    TorchTensorView tv{(uint64_t)vec.data_ptr(), (uint32_t)sizeof(float), {vec.size(), 1LL, 4LL}};
    return DeviceTensor2D<float4>(tv);
}

// Copy pool-backed device tensor to a Python-allocated output tensor
template<typename T>
static void copy_result(const DeviceTensor3D<T>& src, const TorchTensorView& dst) {
    if (src.data_ptr() == nullptr || std::get<0>(dst) == 0) return;
    cudaMemcpy((void*)std::get<0>(dst), src.data_ptr(), src.numel() * sizeof(T), cudaMemcpyDeviceToDevice);
}


void set_data_3dgs(
    int64_t num_splats,
    TorchTensorView means,
    TorchTensorView quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView features_dc,
    TorchTensorView features_sh
) {
    Buffers::num_splats = num_splats;
    Buffers::world_means = means;
    Buffers::world_quats = quats;
    Buffers::world_scales = scales;
    Buffers::world_opacities = opacities;
    Buffers::world_features_dc = features_dc;
    Buffers::world_features_sh = features_sh;

    int64_t max_num_splats = std::get<2>(means)[0];
    if (std::get<2>(quats)[0] != max_num_splats ||
        std::get<2>(scales)[0] != max_num_splats ||
        std::get<2>(opacities)[0] != max_num_splats ||
        std::get<2>(features_dc)[0] != max_num_splats ||
        std::get<2>(features_sh)[0] != max_num_splats)
        throw std::runtime_error("setData3DGS: max_num_splats mismatch across splat tensors");
    if (max_num_splats < num_splats)
        throw std::runtime_error("setData3DGS: tensor size smaller than num_splats");
}


void set_camera_params(
    int width,
    int height,
    std::string camera_model,
    TorchTensorView viewmats,
    TorchTensorView intrins,
    TorchTensorView dist_coeffs
) {
    Buffers::image_width = width;
    Buffers::image_height = height;
    Buffers::camera_model = cmt(camera_model);
    Buffers::camera_model_str = camera_model;
    Buffers::camera_viewmats = viewmats;
    Buffers::camera_intrins = intrins;
    Buffers::camera_dist_coeffs = dist_coeffs;
    Buffers::num_cameras = std::get<2>(viewmats)[0];

    if (std::get<2>(intrins)[0] != Buffers::num_cameras ||
        std::get<2>(dist_coeffs)[0] != Buffers::num_cameras)
        throw std::runtime_error("setCameraParams: num_cameras mismatch");
}


void forward_3dgs(
    std::string primitive,   // "3dgs", "mip", "3dgut"
    int sh_degree,
    bool packed,
    TorchTensorView out_rgb,    // [C, H, W, 3]  float32, pre-allocated by Python
    TorchTensorView out_depth,  // [C, H, W, 1]  float32, pre-allocated by Python
    TorchTensorView out_Ts      // [C, H, W, 1]  float32, pre-allocated by Python
) {
    Buffers::primitive = primitive;
    Buffers::sh_degree = sh_degree;
    Buffers::packed = packed;

    std::vector<DeviceTensorFloatND> in_splats = {
        tv_to_fnd(Buffers::world_means),
        tv_to_fnd(Buffers::world_quats),
        tv_to_fnd(Buffers::world_scales),
        tv_to_fnd(Buffers::world_opacities),
        tv_to_fnd(Buffers::world_features_dc),
        tv_to_fnd(Buffers::world_features_sh),
    };
    Buffers::fwd_splats_w = in_splats;

    DeviceTensorFloatND aabb_nd, depths_nd;

    // --- Projection ---
    if (packed) {
        DeviceVector<int32_t> cam_ids, gauss_ids;
        DeviceVector<float4> aabb_vec;
        DeviceVector<float> depths_vec;

        if (primitive == "3dgs") {
            auto [a, b, c, d, e, f] = projection_3dgs_packed_forward(
                Buffers::num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = f;
        } else if (primitive == "mip") {
            auto [a, b, c, d, e, f] = projection_mip_packed_forward(
                Buffers::num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = f;
        } else if (primitive == "3dgut") {
            auto [a, b, c, d, e, f] = projection_3dgut_packed_forward(
                Buffers::num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = f;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        Buffers::fwd_camera_ids = std::make_optional(cam_ids);
        Buffers::fwd_gaussian_ids = std::make_optional(gauss_ids);
        Buffers::fwd_aabb = vec_to_2d_float4(aabb_vec);  // [nnz, 1] view for backward
        aabb_nd = DeviceTensorFloatND(aabb_vec);          // [nnz, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_vec);      // [nnz]   for intersect
    } else {
        DeviceTensor2D<float4> aabb_2d;
        DeviceTensor2D<float> depths_2d;  // [C, N] — sorting depths per camera

        if (primitive == "3dgs") {
            auto [a, b, c, d] = projection_3dgs_forward(
                Buffers::num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs);
            aabb_2d = a; depths_2d = b;  // b = sorting_depths[C,N], c = radii[N] (unused)
            Buffers::fwd_splats_s = d;
        } else if (primitive == "mip") {
            auto [a, b, c, d] = projection_mip_forward(
                Buffers::num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = d;
        } else if (primitive == "3dgut") {
            auto [a, b, c, d] = projection_3dgut_forward(
                Buffers::num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = d;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        Buffers::fwd_camera_ids = std::nullopt;
        Buffers::fwd_gaussian_ids = std::nullopt;
        Buffers::fwd_aabb = aabb_2d;                           // [C, N] for backward
        aabb_nd = DeviceTensorFloatND(aabb_2d);                // [C, N, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_2d);            // [C, N] → numel=C*N for intersect
    }

    // --- Tile intersection (AABB mode) ---
    DeviceVector<int32_t>* image_ids_ptr = nullptr;
    if (packed && Buffers::fwd_camera_ids.has_value())
        image_ids_ptr = &Buffers::fwd_camera_ids.value();

    auto [isect_ids, flatten_ids, tile_offsets] = do_intersect_tile_generic(
        aabb_nd, depths_nd,
        nullptr, nullptr, nullptr,
        (uint32_t)Buffers::num_cameras,
        Buffers::camera_intrins,
        (uint32_t)Buffers::image_width,
        (uint32_t)Buffers::image_height,
        image_ids_ptr
    );

    Buffers::fwd_tile_offsets = tile_offsets;
    Buffers::fwd_flatten_ids = flatten_ids;

    // --- Rasterization forward ---
    RenderOutput::TensorTuple renders;
    DeviceTensor3D<float> render_Ts;
    DeviceTensor3D<int32_t> last_ids;

    if (primitive == "3dgs") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_3dgs_fwd(
            in_splats, Buffers::fwd_splats_s, Buffers::fwd_gaussian_ids,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    } else if (primitive == "mip") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_mip_fwd(
            in_splats, Buffers::fwd_splats_s, Buffers::fwd_gaussian_ids,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    } else if (primitive == "3dgut") {
        auto [r, rTs, lids, r2, dist] = rasterize_to_pixels_3dgut_fwd(
            in_splats, Buffers::fwd_splats_s, Buffers::fwd_gaussian_ids,
            Buffers::camera_viewmats, Buffers::camera_intrins,
            Buffers::camera_model_str, Buffers::camera_dist_coeffs,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            tile_offsets, flatten_ids, false);
        renders = r; render_Ts = rTs; last_ids = lids;
    }

    Buffers::fwd_renders = renders;
    Buffers::fwd_render_Ts = render_Ts;
    Buffers::fwd_last_ids = last_ids;

    // Copy results to Python output tensors
    copy_result(std::get<0>(renders), out_rgb);    // rgb  [C, H, W] float3
    copy_result(std::get<1>(renders), out_depth);  // depth [C, H, W] float
    copy_result(render_Ts, out_Ts);                // Ts   [C, H, W] float
}


void backward_3dgs(
    // Gradient of rendered outputs (from Python, pointing into grad tensors)
    TorchTensorView v_rgb,    // [C, H, W, 3]  float32
    TorchTensorView v_depth,  // [C, H, W, 1]  float32
    TorchTensorView v_Ts,     // [C, H, W, 1]  float32
    // Pre-allocated gradient output tensors for world-space splats (zero-initialized by Python)
    TorchTensorView v_means,       // [N, 3]
    TorchTensorView v_quats,       // [N, 4]
    TorchTensorView v_scales,      // [N, 3]
    TorchTensorView v_opacities,   // [N, 1]
    TorchTensorView v_features_dc, // [N, 3]
    TorchTensorView v_features_sh  // [N, K, 3]
) {
    // v_rgb [C, H, W, 3]: DeviceTensor3D<float3> needs shape[3]*elem_size = 3*4 = sizeof(float3)
    // v_depth [C, H, W, 1]: DeviceTensor3D<float> needs shape[3]*elem_size = 1*4 = sizeof(float)
    RenderOutput::TensorTuple v_render_outputs = std::make_tuple(
        DeviceTensor3D<float3>(v_rgb),
        DeviceTensor3D<float>(v_depth),
        DeviceTensor3D<float3>()  // no normal gradient
    );
    DeviceTensor3D<float> v_render_Ts(v_Ts);

    // Pre-allocated world-space gradient outputs (zero-initialized by Python via zero_grad)
    std::vector<DeviceTensorFloatND> v_splats_w = {
        tv_to_fnd(v_means),
        tv_to_fnd(v_quats),
        tv_to_fnd(v_scales),
        tv_to_fnd(v_opacities),
        tv_to_fnd(v_features_dc),
        tv_to_fnd(v_features_sh),
    };

    std::vector<DeviceTensorFloatND> v_splats_w_out, v_splats_s_out;

    // --- Rasterization backward ---
    if (Buffers::primitive == "3dgs" || Buffers::primitive == "mip") {
        // mip uses the same rasterization backward kernel as 3dgs
        auto [vw, vs, accum_weight] = rasterize_to_pixels_3dgs_bwd(
            Buffers::num_splats,
            Buffers::fwd_splats_w,
            Buffers::fwd_splats_s,
            Buffers::fwd_gaussian_ids,
            (uint32_t)Buffers::image_width,
            (uint32_t)Buffers::image_height,
            Buffers::fwd_tile_offsets,
            Buffers::fwd_flatten_ids,
            Buffers::fwd_render_Ts,
            Buffers::fwd_last_ids,
            Buffers::fwd_renders,
            std::nullopt,  // no accum_weight_map
            v_render_outputs,
            v_render_Ts,
            std::make_optional(v_splats_w),  // pre-allocated: rasterizer writes SH grads here
            std::nullopt   // let rasterizer allocate v_splats_s
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
    } else if (Buffers::primitive == "3dgut") {
        auto [vw, vs, vviewmats, accum_weight] = rasterize_to_pixels_3dgut_bwd(
            Buffers::num_splats,
            Buffers::fwd_splats_w,
            Buffers::fwd_splats_s,
            Buffers::fwd_gaussian_ids,
            Buffers::camera_viewmats,
            Buffers::camera_intrins,
            Buffers::camera_model_str,
            Buffers::camera_dist_coeffs,
            (uint32_t)Buffers::image_width,
            (uint32_t)Buffers::image_height,
            Buffers::fwd_tile_offsets,
            Buffers::fwd_flatten_ids,
            Buffers::fwd_render_Ts,
            Buffers::fwd_last_ids,
            Buffers::fwd_renders,  // render_outputs needed for 3dgut
            std::nullopt,  // render2_outputs
            std::nullopt,  // loss_map
            std::nullopt,  // accum_weight_map
            v_render_outputs,
            v_render_Ts,
            std::nullopt,  // v_distortion_outputs
            std::make_optional(v_splats_w),  // pre-allocated
            std::nullopt,  // let rasterizer allocate v_splats_s
            false  // need_viewmat_grad
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
    } else {
        throw std::runtime_error("engine_backward: unknown primitive: " + Buffers::primitive);
    }

    // --- Projection backward: accumulates geometric gradients into v_splats_w_out ---
    if (Buffers::primitive == "3dgs") {
        projection_3dgs_backward(
            Buffers::num_splats,
            Buffers::sh_degree,
            Buffers::fwd_splats_w,
            Buffers::camera_viewmats,
            Buffers::camera_intrins,
            (uint32_t)Buffers::image_width,
            (uint32_t)Buffers::image_height,
            Buffers::camera_model_str,
            Buffers::camera_dist_coeffs,
            Buffers::fwd_camera_ids,
            Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb,
            v_splats_s_out,
            v_splats_w_out,
            nullptr  // v_viewmats: not needed
        );
    } else if (Buffers::primitive == "mip") {
        projection_mip_backward(
            Buffers::num_splats,
            Buffers::sh_degree,
            Buffers::fwd_splats_w,
            Buffers::camera_viewmats,
            Buffers::camera_intrins,
            (uint32_t)Buffers::image_width,
            (uint32_t)Buffers::image_height,
            Buffers::camera_model_str,
            Buffers::camera_dist_coeffs,
            Buffers::fwd_camera_ids,
            Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb,
            v_splats_s_out,
            v_splats_w_out,
            nullptr
        );
    } else if (Buffers::primitive == "3dgut") {
        projection_3dgut_backward(
            Buffers::num_splats,
            Buffers::sh_degree,
            Buffers::fwd_splats_w,
            Buffers::camera_viewmats,
            Buffers::camera_intrins,
            (uint32_t)Buffers::image_width,
            (uint32_t)Buffers::image_height,
            Buffers::camera_model_str,
            Buffers::camera_dist_coeffs,
            Buffers::fwd_camera_ids,
            Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb,
            v_splats_s_out,
            v_splats_w_out,
            nullptr
        );
    }
}

