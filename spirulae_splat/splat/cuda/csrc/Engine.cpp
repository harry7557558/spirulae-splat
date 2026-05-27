#include "Engine.h"
#include "Tensor.h"
#include "Primitive.cuh"
#include <map>
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
int64_t cur_num_splats;
int64_t max_num_splats;
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
DeviceVector<int32_t> fwd_camera_ids;
DeviceVector<int32_t> fwd_gaussian_ids;
DeviceTensor2D<float4> fwd_aabb;   // [nnz,1] packed or [C,N] non-packed
std::vector<DeviceTensorFloatND> fwd_splats_w;
std::vector<DeviceTensorFloatND> fwd_splats_s;
DeviceTensor3D<int32_t> fwd_tile_offsets;
DeviceVector<int32_t> fwd_flatten_ids;
DeviceTensor3D<float> fwd_render_Ts;
DeviceTensor3D<int32_t> fwd_last_ids;
RenderOutput::TensorTuple fwd_renders;  // for 3dgut backward
DeviceVector<float> fwd_accum_weight;  // [max_num_splats] per-splat score from raster bwd

// Training ground truth (set per batch via set_training_data)
TorchTensorView gt_rgb;
TorchTensorView gt_depth;
TorchTensorView gt_normal;
TorchTensorView gt_alpha;
TorchTensorView gt_rgb_mask;
TorchTensorView gt_depth_mask;
TorchTensorView gt_normal_mask;
TorchTensorView gt_alpha_mask;
bool has_gt = false;

// Pool-backed training state
int num_sh = 0;
bool optim_initialized = false;
bool train_quantize_sh = false;

// Gradients (pool-backed, zeroed each backward)
TorchTensorView grad_means, grad_quats, grad_scales, grad_opacities;
TorchTensorView grad_features_dc, grad_features_sh;

// Optimizer state g1/g2 (pool-backed, persistent across steps)
TorchTensorView g1_means, g1_quats, g1_scales, g1_opacities;
TorchTensorView g1_features_dc, g1_features_sh;
TorchTensorView g2_means, g2_quats, g2_scales, g2_opacities;
TorchTensorView g2_features_dc, g2_features_sh;

// Training aux (pool-backed)
DeviceVector<float> radii;         // [max_N] float
TorchTensorView accum_buffer;      // [max_N, 2] float
TorchTensorView quant_bounds_sh;   // [n_blocks, 4] float, or null
DeviceVector<int32_t> bias_correction_steps;  // [max_N] int32, or null
bool use_per_splat_bias_correction = false;

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
    Buffers::cur_num_splats = num_splats;
    Buffers::max_num_splats = std::get<2>(means)[0];
    Buffers::world_means = means;
    Buffers::world_quats = quats;
    Buffers::world_scales = scales;
    Buffers::world_opacities = opacities;
    Buffers::world_features_dc = features_dc;
    Buffers::world_features_sh = features_sh;

    auto sh_shape = std::get<2>(features_sh);
    Buffers::num_sh = (sh_shape.size() >= 2) ? (int)sh_shape[1] : 0;

    if (std::get<2>(quats)[0] != Buffers::max_num_splats ||
        std::get<2>(scales)[0] != Buffers::max_num_splats ||
        std::get<2>(opacities)[0] != Buffers::max_num_splats ||
        std::get<2>(features_dc)[0] != Buffers::max_num_splats ||
        std::get<2>(features_sh)[0] != Buffers::max_num_splats)
        throw std::runtime_error("setData3DGS: max_num_splats mismatch across splat tensors");
    if (Buffers::max_num_splats < num_splats)
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

    // Allocate and zero radii buffer before projection (kernel uses atomicMax)
    Buffers::radii.resize("eng.radii", Buffers::max_num_splats);
    Buffers::radii.zero();

    // --- Projection ---
    if (packed) {
        DeviceVector<int32_t> cam_ids, gauss_ids;
        DeviceVector<float4> aabb_vec;
        DeviceVector<float> depths_vec;

        if (primitive == "3dgs") {
            auto [a, b, c, d, e] = projection_3dgs_packed_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs,
                Buffers::radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = e;
        } else if (primitive == "mip") {
            auto [a, b, c, d, e] = projection_mip_packed_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs,
                Buffers::radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = e;
        } else if (primitive == "3dgut") {
            auto [a, b, c, d, e] = projection_3dgut_packed_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs,
                Buffers::radii);
            cam_ids = a; gauss_ids = b; aabb_vec = c; depths_vec = d;
            Buffers::fwd_splats_s = e;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        Buffers::fwd_camera_ids = cam_ids;
        Buffers::fwd_gaussian_ids = gauss_ids;
        Buffers::fwd_aabb = vec_to_2d_float4(aabb_vec);  // [nnz, 1] view for backward
        aabb_nd = DeviceTensorFloatND(aabb_vec);          // [nnz, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_vec);      // [nnz]   for intersect
    } else {
        DeviceTensor2D<float4> aabb_2d;
        DeviceTensor2D<float> depths_2d;  // [C, N] — sorting depths per camera

        if (primitive == "3dgs") {
            auto [a, b, c] = projection_3dgs_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs,
                Buffers::radii);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = c;
        } else if (primitive == "mip") {
            auto [a, b, c] = projection_mip_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs,
                Buffers::radii);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = c;
        } else if (primitive == "3dgut") {
            auto [a, b, c] = projection_3dgut_forward(
                Buffers::cur_num_splats, sh_degree, in_splats,
                Buffers::camera_viewmats, Buffers::camera_intrins,
                (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
                Buffers::camera_model_str, Buffers::camera_dist_coeffs,
                Buffers::radii);
            aabb_2d = a; depths_2d = b;
            Buffers::fwd_splats_s = c;
        } else {
            throw std::runtime_error("engine_forward: unknown primitive: " + primitive);
        }

        Buffers::fwd_camera_ids = DeviceVector<int32_t>();
        Buffers::fwd_gaussian_ids = DeviceVector<int32_t>();
        Buffers::fwd_aabb = aabb_2d;                           // [C, N] for backward
        aabb_nd = DeviceTensorFloatND(aabb_2d);                // [C, N, 4] for intersect
        depths_nd = DeviceTensorFloatND(depths_2d);            // [C, N] → numel=C*N for intersect
    }

    // --- Tile intersection (AABB mode) ---
    DeviceVector<int32_t>* image_ids_ptr = nullptr;
    if (packed && Buffers::fwd_camera_ids.data_ptr() != nullptr)
        image_ids_ptr = &Buffers::fwd_camera_ids;

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


void engine_debug_forward(
    TorchTensorView override_features_dc,
    int override_sh_degree,
    TorchTensorView out_rgb
) {
    if (std::get<0>(Buffers::fwd_renders).data_ptr() == nullptr)
        throw std::runtime_error("engine_debug_forward: forward_3dgs must be called first");

    // Swap in overrides
    TorchTensorView saved_dc = Buffers::world_features_dc;
    int saved_sh = Buffers::sh_degree;
    if (std::get<0>(override_features_dc) != 0)
        Buffers::world_features_dc = override_features_dc;
    if (override_sh_degree >= 0)
        Buffers::sh_degree = override_sh_degree;

    // Dummy depth/Ts outputs (we only care about rgb)
    TorchTensorView null_tv(0, 4, {0});
    forward_3dgs(Buffers::primitive, Buffers::sh_degree, Buffers::packed,
                 out_rgb, null_tv, null_tv);

    // Restore
    Buffers::world_features_dc = saved_dc;
    Buffers::sh_degree = saved_sh;
}


void engine_copy_accum_buffer(TorchTensorView dst) {
    if (std::get<0>(Buffers::accum_buffer) == 0 || std::get<0>(dst) == 0)
        return;
    auto& dst_shape = std::get<2>(dst);
    int64_t dst_n = dst_shape[0];
    int64_t src_n = std::get<2>(Buffers::accum_buffer)[0];
    int64_t n = std::min(dst_n, src_n);
    cudaMemcpy((void*)std::get<0>(dst),
               (void*)std::get<0>(Buffers::accum_buffer),
               n * 2 * sizeof(float), cudaMemcpyDeviceToDevice);
}


int64_t engine_get_cur_num_splats() {
    return Buffers::cur_num_splats;
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
            Buffers::cur_num_splats,
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
            DeviceTensor3D<float>(),  // no accum_weight_map
            v_render_outputs,
            v_render_Ts,
            std::make_optional(v_splats_w),
            std::nullopt
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
    } else if (Buffers::primitive == "3dgut") {
        auto [vw, vs, vviewmats, accum_weight] = rasterize_to_pixels_3dgut_bwd(
            Buffers::cur_num_splats,
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
            Buffers::fwd_renders,
            std::nullopt,  // render2_outputs
            DeviceTensor3D<float>(),  // loss_map
            DeviceTensor3D<float>(),  // accum_weight_map
            v_render_outputs,
            v_render_Ts,
            std::nullopt,  // v_distortion_outputs
            std::make_optional(v_splats_w),
            std::nullopt,
            false
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
    } else {
        throw std::runtime_error("engine_backward: unknown primitive: " + Buffers::primitive);
    }

    // --- Projection backward: accumulates geometric gradients into v_splats_w_out ---
    if (Buffers::primitive == "3dgs") {
        projection_3dgs_backward(
            Buffers::cur_num_splats,
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
            Buffers::cur_num_splats,
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
            Buffers::cur_num_splats,
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


void set_training_data(
    TorchTensorView gt_rgb,
    TorchTensorView gt_depth,
    TorchTensorView gt_normal,
    TorchTensorView gt_alpha,
    TorchTensorView gt_rgb_mask,
    TorchTensorView gt_depth_mask,
    TorchTensorView gt_normal_mask,
    TorchTensorView gt_alpha_mask
) {
    Buffers::gt_rgb = gt_rgb;
    Buffers::gt_depth = gt_depth;
    Buffers::gt_normal = gt_normal;
    Buffers::gt_alpha = gt_alpha;
    Buffers::gt_rgb_mask = gt_rgb_mask;
    Buffers::gt_depth_mask = gt_depth_mask;
    Buffers::gt_normal_mask = gt_normal_mask;
    Buffers::gt_alpha_mask = gt_alpha_mask;
    Buffers::has_gt = (std::get<0>(gt_rgb) != 0);
}


static inline TorchTensorView _tv_null() {
    return TorchTensorView(0, 4, {0});
}

static inline bool _tv_valid(const TorchTensorView& tv) {
    return std::get<0>(tv) != 0;
}

static inline TorchTensorView _pool_tv(const std::string& key, int64_t B, int64_t H, int64_t W, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(B * H * W * C));
    return TorchTensorView((uint64_t)p, 4, {B, H, W, C});
}

static inline TorchTensorView _pool_tv_zero(const std::string& key, int64_t B, int64_t H, int64_t W, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(B * H * W * C));
    cudaMemset(p, 0, B * H * W * C * sizeof(float));
    return TorchTensorView((uint64_t)p, 4, {B, H, W, C});
}

static inline TorchTensorView _pool_tv_1d(const std::string& key, int64_t N) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)N);
    return TorchTensorView((uint64_t)p, 4, {N});
}

static inline TorchTensorView _pool_tv_1d_zero(const std::string& key, int64_t N) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)N);
    cudaMemset(p, 0, N * sizeof(float));
    return TorchTensorView((uint64_t)p, 4, {N});
}

static inline TorchTensorView _pool_tv_2d(const std::string& key, int64_t N, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(N * C));
    return TorchTensorView((uint64_t)p, 4, {N, C});
}

static inline TorchTensorView _pool_tv_3d_f(const std::string& key, int64_t N, int64_t K, int64_t C) {
    float* p = DevicePool::global().acquire<float>(key, (size_t)(N * K * C));
    return TorchTensorView((uint64_t)p, 4, {N, K, C});
}

static inline TorchTensorView _pool_tv_3d_u8(const std::string& key, int64_t N, int64_t K, int64_t C) {
    uint8_t* p = DevicePool::global().acquire<uint8_t>(key, (size_t)(N * K * C));
    return TorchTensorView((uint64_t)p, 1, {N, K, C});
}

static void _zero_tv(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return;
    int64_t bytes = std::get<1>(tv);
    for (auto s : std::get<2>(tv)) bytes *= s;
    cudaMemset((void*)std::get<0>(tv), 0, bytes);
}

static void _alloc_grad_buffers() {
    int64_t N = Buffers::max_num_splats;
    int64_t K = Buffers::num_sh;
    Buffers::grad_means       = _pool_tv_2d("eng.v_means", N, 3);
    Buffers::grad_quats       = _pool_tv_2d("eng.v_quats", N, 4);
    Buffers::grad_scales      = _pool_tv_2d("eng.v_scales", N, 3);
    Buffers::grad_opacities   = _pool_tv_2d("eng.v_opacities", N, 1);
    Buffers::grad_features_dc = _pool_tv_2d("eng.v_features_dc", N, 3);
    Buffers::grad_features_sh = _pool_tv_3d_f("eng.v_features_sh", N, K, 3);
    _zero_tv(Buffers::grad_means);
    _zero_tv(Buffers::grad_quats);
    _zero_tv(Buffers::grad_scales);
    _zero_tv(Buffers::grad_opacities);
    _zero_tv(Buffers::grad_features_dc);
    _zero_tv(Buffers::grad_features_sh);
}

static void _ensure_optim_state(bool quantize_sh, bool use_per_splat_bias_correction = false) {
    if (Buffers::optim_initialized && Buffers::train_quantize_sh == quantize_sh
        && Buffers::use_per_splat_bias_correction == use_per_splat_bias_correction)
        return;

    int64_t N = Buffers::max_num_splats;
    int64_t K = Buffers::num_sh;
    Buffers::train_quantize_sh = quantize_sh;

    // g1 (exp_avg)
    Buffers::g1_means       = _pool_tv_2d("eng.g1_means", N, 3);
    Buffers::g1_quats       = _pool_tv_2d("eng.g1_quats", N, 4);
    Buffers::g1_scales      = _pool_tv_2d("eng.g1_scales", N, 3);
    Buffers::g1_opacities   = _pool_tv_2d("eng.g1_opacities", N, 1);
    Buffers::g1_features_dc = _pool_tv_2d("eng.g1_features_dc", N, 3);
    Buffers::g1_features_sh = quantize_sh
        ? _pool_tv_3d_u8("eng.g1_features_sh", N, K, 3)
        : _pool_tv_3d_f("eng.g1_features_sh", N, K, 3);

    // g2 (exp_avg_sq)
    Buffers::g2_means       = _pool_tv_2d("eng.g2_means", N, 3);
    Buffers::g2_quats       = _pool_tv_2d("eng.g2_quats", N, 4);
    Buffers::g2_scales      = _pool_tv_2d("eng.g2_scales", N, 3);
    Buffers::g2_opacities   = _pool_tv_2d("eng.g2_opacities", N, 1);
    Buffers::g2_features_dc = _pool_tv_2d("eng.g2_features_dc", N, 3);
    Buffers::g2_features_sh = quantize_sh
        ? _pool_tv_3d_u8("eng.g2_features_sh", N, K, 3)
        : _pool_tv_3d_f("eng.g2_features_sh", N, K, 3);

    // radii [max_N]
    Buffers::radii.resize("eng.radii", N);

    // accum_buffer [max_N, 2]
    Buffers::accum_buffer = _pool_tv_2d("eng.accum_buffer", N, 2);

    // quant_bounds_sh
    if (quantize_sh) {
        int64_t sh_numel = N * K * 3;
        int64_t n_blocks = (sh_numel + 255) / 256;
        Buffers::quant_bounds_sh = _pool_tv_2d("eng.quant_bounds_sh", n_blocks, 4);
    } else {
        Buffers::quant_bounds_sh = _tv_null();
    }

    // bias_correction_steps
    Buffers::use_per_splat_bias_correction = use_per_splat_bias_correction;
    if (use_per_splat_bias_correction) {
        Buffers::bias_correction_steps.resize("eng.bias_correction_steps", N);
    } else {
        Buffers::bias_correction_steps = DeviceVector<int32_t>();
    }

    // Zero everything on first init
    _zero_tv(Buffers::g1_means);       _zero_tv(Buffers::g2_means);
    _zero_tv(Buffers::g1_quats);       _zero_tv(Buffers::g2_quats);
    _zero_tv(Buffers::g1_scales);      _zero_tv(Buffers::g2_scales);
    _zero_tv(Buffers::g1_opacities);   _zero_tv(Buffers::g2_opacities);
    _zero_tv(Buffers::g1_features_dc); _zero_tv(Buffers::g2_features_dc);
    _zero_tv(Buffers::g1_features_sh); _zero_tv(Buffers::g2_features_sh);
    _zero_tv(Buffers::accum_buffer);
    if (_tv_valid(Buffers::quant_bounds_sh))
        _zero_tv(Buffers::quant_bounds_sh);
    Buffers::bias_correction_steps.zero();

    Buffers::optim_initialized = true;
}


std::map<std::string, float> engine_compute_loss_backward(
    int step,
    std::array<float, (int)LossWeightIndex::length> loss_weights,
    float w_ssim,
    int num_loss_scales,
    bool compute_loss_map
) {
    // Validate that forward was run
    if (std::get<0>(Buffers::fwd_renders) .data_ptr() == nullptr)
        throw std::runtime_error("engine_compute_loss_backward: forward_3dgs must be called first");
    if (!Buffers::has_gt)
        throw std::runtime_error("engine_compute_loss_backward: set_training_data must be called first");

    // Allocate and zero gradient buffers from pool
    _alloc_grad_buffers();

    int64_t C = Buffers::num_cameras;
    int64_t H = Buffers::image_height;
    int64_t W = Buffers::image_width;

    // Pool-allocate intermediates for loss computation
    TorchTensorView loss_map_buf = compute_loss_map ?
        _pool_tv_zero("eng.loss_map", C, H, W, 1) : _tv_null();
    TorchTensorView total_losses_buf = _pool_tv_1d_zero("eng.total_losses", (int)LossIndex::length);

    // v_losses: constant vector [1, 0, 1, 1, ...] (1 for all, 0 for psnr slot)
    // Initialized once; pool never shrinks so pointer is stable
    static bool v_losses_initialized = false;
    TorchTensorView v_losses_buf = _pool_tv_1d("eng.v_losses", (int)LossIndex::length);
    if (!v_losses_initialized) {
        float h_v[(int)LossIndex::length];
        for (int i = 0; i < (int)LossIndex::length; i++) h_v[i] = 1.0f;
        h_v[(int)LossIndex::RgbPSNR] = 0.0f;
        cudaMemcpy((void*)std::get<0>(v_losses_buf), h_v, sizeof(h_v), cudaMemcpyHostToDevice);
        v_losses_initialized = true;
    }

    // Render outputs from forward pass (pool-backed, already populated)
    TorchTensorView render_rgb = TorchTensorView(
        (uint64_t)std::get<0>(Buffers::fwd_renders).data_ptr(), 4, {C, H, W, 3});
    TorchTensorView render_depth = TorchTensorView(
        (uint64_t)std::get<1>(Buffers::fwd_renders).data_ptr(), 4, {C, H, W, 1});
    TorchTensorView render_Ts = TorchTensorView(
        (uint64_t)Buffers::fwd_render_Ts.data_ptr(), 4, {C, H, W, 1});

    // Render normal: not yet available from forward (needs distortion output)
    // TODO: pass distortion outputs from forward when enabled
    TorchTensorView render_normal = _tv_null();
    TorchTensorView depth_normal = _tv_null();
    TorchTensorView rgb_dist = _tv_null();
    TorchTensorView depth_dist = _tv_null();
    TorchTensorView normal_dist = _tv_null();

    std::vector<bool> needs_input_grad = {
        true,   // pred_rgb
        false,  // gt_rgb
        true,   // pred_depth
        false,  // gt_depth  (TODO: true when bilateral grid for geometry)
        true,   // pred_normal
        true,   // pred_depth_normal
        false,  // gt_normal (TODO: true when bilateral grid for geometry)
        true,   // pred_transmittance
        true, true, true,  // distortion
    };

    PerPixelGrads pixel_grads = {};

    // Allocate per-pixel gradient outputs
    pixel_grads.v_render_rgb  = _pool_tv("eng.v_rgb",   C, H, W, 3);
    pixel_grads.v_render_depth = _pool_tv("eng.v_depth", C, H, W, 1);
    pixel_grads.v_render_Ts   = _pool_tv("eng.v_Ts",    C, H, W, 1);
    // TODO: allocate normal/distortion grads when those features are enabled

    // --- Compute per-pixel losses + SSIM, get gradients ---
    auto [psnr, ssim] = compute_multi_scale_per_pixel_losses(
        num_loss_scales,
        render_rgb,
        Buffers::gt_rgb,
        render_depth,
        Buffers::gt_depth,
        render_normal,
        depth_normal,
        Buffers::gt_normal,
        render_Ts,
        rgb_dist,
        depth_dist,
        normal_dist,
        Buffers::gt_alpha,
        Buffers::gt_rgb_mask,
        Buffers::gt_depth_mask,
        Buffers::gt_normal_mask,
        Buffers::gt_alpha_mask,
        loss_weights,
        w_ssim,
        v_losses_buf,
        needs_input_grad,
        -1,  // num_train_images: -1 means use batch size
        _tv_null(),  // camera_indices: null means identity
        loss_map_buf,
        total_losses_buf,
        pixel_grads
    );

    // TODO: bilateral grid forward/backward
    // TODO: PPISP forward/backward
    // TODO: color space conversion (rgb_to_srgb) forward/backward
    // TODO: background blending forward/backward
    // TODO: depth_to_normal backward when depth_normal gradient is non-null

    // --- Rasterization + projection backward ---
    // Build v_render_outputs from pixel_grads
    RenderOutput::TensorTuple v_render_outputs = std::make_tuple(
        DeviceTensor3D<float3>(pixel_grads.v_render_rgb),
        DeviceTensor3D<float>(pixel_grads.v_render_depth),
        DeviceTensor3D<float3>()  // no normal gradient yet
    );
    DeviceTensor3D<float> v_render_Ts(pixel_grads.v_render_Ts);

    std::vector<DeviceTensorFloatND> v_splats_w = {
        tv_to_fnd(Buffers::grad_means),
        tv_to_fnd(Buffers::grad_quats),
        tv_to_fnd(Buffers::grad_scales),
        tv_to_fnd(Buffers::grad_opacities),
        tv_to_fnd(Buffers::grad_features_dc),
        tv_to_fnd(Buffers::grad_features_sh),
    };

    // Build accum_weight_map from loss_map (pixel-space -> per-splat mapping in raster bwd)
    DeviceTensor3D<float> accum_weight_map;
    if (compute_loss_map && _tv_valid(loss_map_buf)) {
        accum_weight_map = DeviceTensor3D<float>(loss_map_buf);
    }

    std::vector<DeviceTensorFloatND> v_splats_w_out, v_splats_s_out;

    if (Buffers::primitive == "3dgs" || Buffers::primitive == "mip") {
        auto [vw, vs, accum_weight] = rasterize_to_pixels_3dgs_bwd(
            Buffers::cur_num_splats,
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
            accum_weight_map,
            v_render_outputs,
            v_render_Ts,
            std::make_optional(v_splats_w),
            std::nullopt
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
        if (accum_weight.data_ptr()) Buffers::fwd_accum_weight = accum_weight;
    } else if (Buffers::primitive == "3dgut") {
        auto [vw, vs, vviewmats, accum_weight] = rasterize_to_pixels_3dgut_bwd(
            Buffers::cur_num_splats,
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
            Buffers::fwd_renders,
            std::nullopt,  // render2_outputs
            DeviceTensor3D<float>(),  // loss_map
            accum_weight_map,
            v_render_outputs,
            v_render_Ts,
            std::nullopt,  // v_distortion_outputs
            std::make_optional(v_splats_w),
            std::nullopt,
            false
        );
        v_splats_w_out = vw;
        v_splats_s_out = vs;
        if (accum_weight.data_ptr()) Buffers::fwd_accum_weight = accum_weight;
    } else {
        throw std::runtime_error("engine_compute_loss_backward: unknown primitive: " + Buffers::primitive);
    }

    // --- Projection backward ---
    if (Buffers::primitive == "3dgs") {
        projection_3dgs_backward(
            Buffers::cur_num_splats, Buffers::sh_degree, Buffers::fwd_splats_w,
            Buffers::camera_viewmats, Buffers::camera_intrins,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            Buffers::camera_model_str, Buffers::camera_dist_coeffs,
            Buffers::fwd_camera_ids, Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb, v_splats_s_out, v_splats_w_out, nullptr);
    } else if (Buffers::primitive == "mip") {
        projection_mip_backward(
            Buffers::cur_num_splats, Buffers::sh_degree, Buffers::fwd_splats_w,
            Buffers::camera_viewmats, Buffers::camera_intrins,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            Buffers::camera_model_str, Buffers::camera_dist_coeffs,
            Buffers::fwd_camera_ids, Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb, v_splats_s_out, v_splats_w_out, nullptr);
    } else if (Buffers::primitive == "3dgut") {
        projection_3dgut_backward(
            Buffers::cur_num_splats, Buffers::sh_degree, Buffers::fwd_splats_w,
            Buffers::camera_viewmats, Buffers::camera_intrins,
            (uint32_t)Buffers::image_width, (uint32_t)Buffers::image_height,
            Buffers::camera_model_str, Buffers::camera_dist_coeffs,
            Buffers::fwd_camera_ids, Buffers::fwd_gaussian_ids,
            Buffers::fwd_aabb, v_splats_s_out, v_splats_w_out, nullptr);
    }

    // --- Build loss dict for display ---
    float h_losses[(int)LossIndex::length];
    cudaMemcpy(h_losses, (void*)std::get<0>(total_losses_buf),
               sizeof(h_losses), cudaMemcpyDeviceToHost);

    auto sdiv = [](float x, float y) -> float { return y != 0.0f ? x / y : 0.0f; };

    float rgb_loss = h_losses[(int)LossIndex::RgbLoss];

    std::map<std::string, float> loss_dict;
    float psnr_from_losses = h_losses[(int)LossIndex::RgbPSNR];
    loss_dict["rgb_loss"] = rgb_loss + w_ssim * (1.0f - ssim);
    loss_dict["psnr"] = psnr_from_losses;
    loss_dict["ssim"] = ssim;
    loss_dict["depth_loss"] = sdiv(h_losses[(int)LossIndex::DepthSup],
        loss_weights[(int)LossWeightIndex::DepthSup]);
    loss_dict["normal_loss"] = sdiv(h_losses[(int)LossIndex::NormalSup],
        loss_weights[(int)LossWeightIndex::NormalSup]);
    loss_dict["alpha_loss"] = sdiv(h_losses[(int)LossIndex::AlphaSup],
        loss_weights[(int)LossWeightIndex::AlphaSup] + loss_weights[(int)LossWeightIndex::AlphaSupUnder]);
    loss_dict["normal_reg"] = sdiv(h_losses[(int)LossIndex::NormalReg],
        loss_weights[(int)LossWeightIndex::NormalReg]);
    loss_dict["alpha_reg"] = sdiv(h_losses[(int)LossIndex::AlphaReg],
        loss_weights[(int)LossWeightIndex::AlphaReg]);
    loss_dict["rgb_dist_reg"] = sdiv(h_losses[(int)LossIndex::RgbDistReg],
        loss_weights[(int)LossWeightIndex::RgbDistReg]);
    loss_dict["depth_dist_reg"] = sdiv(h_losses[(int)LossIndex::DepthDistReg],
        loss_weights[(int)LossWeightIndex::DepthDistReg]);
    loss_dict["normal_dist_reg"] = sdiv(h_losses[(int)LossIndex::NormalDistReg],
        loss_weights[(int)LossWeightIndex::NormalDistReg]);

    return loss_dict;
}


void engine_optim_step(
    int step,
    float lr_means, float lr_quats, float lr_scales, float lr_opacities,
    float lr_features_dc, float lr_features_sh,
    float max_gauss_ratio, float scale_regularization_weight,
    float mcmc_opacity_reg_weight, float mcmc_scale_reg_weight,
    float erank_reg_weight, float erank_reg_weight_s3,
    float quat_norm_reg_weight, float sh_reg_weight,
    bool use_scale_agnostic_mean, bool quantize_sh,
    bool use_per_splat_bias_correction
) {
    _ensure_optim_state(quantize_sh, use_per_splat_bias_correction);

    int64_t N = Buffers::cur_num_splats;

    // Per-splat bias correction: increment all steps by 1 before use
    DeviceVector<int32_t> per_splat_steps;
    if (use_per_splat_bias_correction) {
        increment_int32_inplace(Buffers::bias_correction_steps, Buffers::max_num_splats);
        per_splat_steps = Buffers::bias_correction_steps;
    }

    fused_optim_3dgs_geometry(
        N,
        DeviceVector<float3>(Buffers::world_means),     DeviceVector<float3>(Buffers::grad_means),
        DeviceVector<float3>(Buffers::g1_means),        DeviceVector<float3>(Buffers::g2_means),
        DeviceVector<float4>(Buffers::world_quats),     DeviceVector<float4>(Buffers::grad_quats),
        DeviceVector<float4>(Buffers::g1_quats),        DeviceVector<float4>(Buffers::g2_quats),
        DeviceVector<float3>(Buffers::world_scales),    DeviceVector<float3>(Buffers::grad_scales),
        DeviceVector<float3>(Buffers::g1_scales),       DeviceVector<float3>(Buffers::g2_scales),
        DeviceVector<float>(Buffers::world_opacities),  DeviceVector<float>(Buffers::grad_opacities),
        DeviceVector<float>(Buffers::g1_opacities),     DeviceVector<float>(Buffers::g2_opacities),
        Buffers::radii,
        lr_means, lr_quats, lr_scales, lr_opacities,
        max_gauss_ratio, scale_regularization_weight,
        mcmc_opacity_reg_weight, mcmc_scale_reg_weight,
        erank_reg_weight, erank_reg_weight_s3, quat_norm_reg_weight,
        use_scale_agnostic_mean,
        step + 1, per_splat_steps
    );

    fused_adam_step(N,
        tv_to_fnd(Buffers::world_features_dc), tv_to_fnd(Buffers::grad_features_dc),
        tv_to_fnd(Buffers::g1_features_dc), tv_to_fnd(Buffers::g2_features_dc),
        lr_features_dc, step + 1, per_splat_steps,
        sh_reg_weight, 0.5f / 0.28209479177387814f);

    if (quantize_sh && _tv_valid(Buffers::quant_bounds_sh)) {
        fused_adam_step_8bit(N,
            tv_to_fnd(Buffers::world_features_sh), tv_to_fnd(Buffers::grad_features_sh),
            (uint8_t*)std::get<0>(Buffers::g1_features_sh),
            (uint8_t*)std::get<0>(Buffers::g2_features_sh),
            (float4*)std::get<0>(Buffers::quant_bounds_sh),
            lr_features_sh, step + 1, per_splat_steps,
            sh_reg_weight, 0.0f);
    } else {
        fused_adam_step(N,
            tv_to_fnd(Buffers::world_features_sh), tv_to_fnd(Buffers::grad_features_sh),
            tv_to_fnd(Buffers::g1_features_sh), tv_to_fnd(Buffers::g2_features_sh),
            lr_features_sh, step + 1, per_splat_steps,
            sh_reg_weight, 0.0f);
    }
}


// Helper: build DeviceVector<T> from TorchTensorView (non-owning view)
template<typename T>
static inline DeviceVector<T> _dv(const TorchTensorView& tv) {
    return DeviceVector<T>(tv);
}
// For null-able TorchTensorViews, return default (null) DeviceVector
template<typename T>
static inline DeviceVector<T> _dv_or_null(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return DeviceVector<T>();
    return DeviceVector<T>(tv);
}
// For [N, K, C] tensors like features_sh: flatten to [N*K, C] then build DeviceVector
template<typename T>
static inline DeviceVector<T> _dv_flat(const TorchTensorView& tv) {
    if (std::get<0>(tv) == 0) return DeviceVector<T>();
    const auto& shape = std::get<2>(tv);
    int64_t flat_count = 1;
    for (size_t i = 0; i + 1 < shape.size(); i++) flat_count *= shape[i];
    TorchTensorView flat_tv(std::get<0>(tv), std::get<1>(tv),
        {flat_count, (int64_t)(sizeof(T) / std::get<1>(tv))});
    return DeviceVector<T>(flat_tv);
}


int engine_densify_step(
    int step,
    int max_steps,
    int refine_start_iter,
    int refine_stop_num_iter,
    int refine_every,
    float growth_factor,
    float min_opacity,
    float max_screen_size,
    float max_screen_size_clip_hardness,
    float max_world_size,
    float noise_lr,
    float noise_lr_final,
    float relocate_heuristic_weight
) {
    // _ensure_optim_state(Buffers::train_quantize_sh);

    int64_t cur_num_splats = Buffers::cur_num_splats;
    int64_t max_num_splats = Buffers::max_num_splats;
    bool quantize_sh = Buffers::train_quantize_sh;
    int num_sh = Buffers::num_sh;

    bool use_revised = (relocate_heuristic_weight >= 1.0f);
    bool densify_ongoing = (step < max_steps - refine_stop_num_iter);
    bool do_densify = densify_ongoing && (step > refine_start_iter && step % refine_every == 0);
    float progress = ((float)step + 0.5f) / (float)max_steps;

    // Build DeviceVectors from pool-backed and world TorchTensorViews
    auto dv_means = _dv<float3>(Buffers::world_means);
    auto dv_quats = _dv<float4>(Buffers::world_quats);
    auto dv_scales = _dv<float3>(Buffers::world_scales);
    auto dv_opacs = _dv<float>(Buffers::world_opacities);
    auto dv_features_dc = _dv<float3>(Buffers::world_features_dc);
    auto dv_features_sh = _dv_flat<float3>(Buffers::world_features_sh);
    auto dv_g1_means = _dv<float3>(Buffers::g1_means);
    auto dv_g1_quats = _dv<float4>(Buffers::g1_quats);
    auto dv_g1_scales = _dv<float3>(Buffers::g1_scales);
    auto dv_g1_opacs = _dv<float>(Buffers::g1_opacities);
    auto dv_g1_features_dc = _dv<float3>(Buffers::g1_features_dc);
    auto dv_g1_features_sh = _dv_flat<float3>(Buffers::g1_features_sh);
    auto dv_g2_means = _dv<float3>(Buffers::g2_means);
    auto dv_g2_quats = _dv<float4>(Buffers::g2_quats);
    auto dv_g2_scales = _dv<float3>(Buffers::g2_scales);
    auto dv_g2_opacs = _dv<float>(Buffers::g2_opacities);
    auto dv_g2_features_dc = _dv<float3>(Buffers::g2_features_dc);
    auto dv_g2_features_sh = _dv_flat<float3>(Buffers::g2_features_sh);
    auto& dv_radii = Buffers::radii;
    auto dv_accum_buf = _dv_or_null<float2>(Buffers::accum_buffer);
    auto dv_bias_steps = Buffers::bias_correction_steps;

    // Clip large splats
    if (std::isfinite(max_screen_size) || std::isfinite(max_world_size)) {
        densify_clip_scale_tensor(
            cur_num_splats,
            dv_radii, dv_scales,
            nullptr,  // don't clip opacities
            max_screen_size,
            max_screen_size_clip_hardness,
            max_world_size
        );
    }

    // Update densification score using accum_weight from rasterization backward
    if (densify_ongoing && use_revised && dv_accum_buf.data_ptr() != nullptr) {
        DeviceVector<float> score;
        if (Buffers::fwd_accum_weight.data_ptr() != nullptr) {
            score = Buffers::fwd_accum_weight;
        } else {
            auto dv_v_opacs = _dv_or_null<float>(Buffers::grad_opacities);
            score = dv_v_opacs.data_ptr() ? dv_v_opacs : dv_opacs;
        }

        densify_update_weight_tensor(
            cur_num_splats,
            dv_radii,
            nullptr,
            (float*)dv_opacs.data_ptr(),
            score,
            dv_accum_buf,
            false
        );
    }

    int num_added = 0;

    if (do_densify && use_revised) {
        // Revised relocation (long axis split)
        relocate_splats_with_long_axis_split_tensor(
            cur_num_splats, min_opacity,
            dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
            dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
            dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
            dv_accum_buf, dv_bias_steps,
            quantize_sh, num_sh, 2 * step + 0
        );

        // Add more splats
        int64_t n_target = std::min(max_num_splats, (int64_t)(growth_factor * cur_num_splats));
        num_added = (int)std::max((int64_t)0, n_target - cur_num_splats);
        if (num_added > 0) {
            add_splats_with_long_axis_split_tensor(
                cur_num_splats, num_added,
                dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
                dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
                dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
                dv_accum_buf, dv_bias_steps,
                quantize_sh, num_sh, 2 * step + 1
            );
        }
    } else if (do_densify) {
        // MCMC relocation
        relocate_splats_mcmc_tensor(
            cur_num_splats, min_opacity,
            dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
            dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
            dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
            dv_bias_steps,
            quantize_sh, num_sh, 2 * step + 0
        );

        // MCMC sample add
        int64_t n_target = std::min(max_num_splats, (int64_t)(growth_factor * cur_num_splats));
        num_added = (int)std::max((int64_t)0, n_target - cur_num_splats);
        if (num_added > 0) {
            add_splats_mcmc_tensor(
                cur_num_splats, num_added, min_opacity,
                dv_means, dv_quats, dv_scales, dv_opacs, dv_features_dc, dv_features_sh,
                dv_g1_means, dv_g1_quats, dv_g1_scales, dv_g1_opacs, dv_g1_features_dc, dv_g1_features_sh,
                dv_g2_means, dv_g2_quats, dv_g2_scales, dv_g2_opacs, dv_g2_features_dc, dv_g2_features_sh,
                dv_bias_steps,
                quantize_sh, num_sh, 2 * step + 1
            );
        }
    }

    // Reset accum buffer after densification step
    if (do_densify && dv_accum_buf.data_ptr() != nullptr) {
        cudaMemset(dv_accum_buf.data_ptr(), 0, dv_accum_buf.size() * sizeof(float2));
    }

    // Add MCMC noise
    if (noise_lr > 0.0f && noise_lr_final > 0.0f) {
        float noise_scalar = noise_lr * powf(noise_lr_final / noise_lr, progress);
        if (use_revised) {
            revised_add_noise_tensor(
                cur_num_splats + num_added, noise_scalar,
                dv_radii, dv_means, dv_scales, dv_quats, dv_opacs);
        } else {
            mcmc_add_noise_tensor(
                cur_num_splats + num_added, noise_scalar,
                dv_means, dv_scales, dv_quats, dv_opacs);
        }
    }

    Buffers::cur_num_splats = cur_num_splats + num_added;
    return num_added;
}

