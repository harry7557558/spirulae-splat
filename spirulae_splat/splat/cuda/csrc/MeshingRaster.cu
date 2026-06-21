/*
 * MeshingRaster.cu  -- see MeshingRaster.cuh.
 *
 * Per-camera render driver: marshals the raw splat params into the projection's
 * in_splats layout and runs
 *     projection_3dgut_forward -> do_intersect_tile_generic
 *     -> rasterize_moments_3dgut_fwd
 * for one camera, reusing the engine's optimized 3DGUT projection + tiling.
 */

#include "MeshingRaster.cuh"

#include <cuda_runtime.h>
#include <algorithm>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <vector>

#include <Tensor.h>
#include <Common.cuh>

#include "ProjectionFwd.cuh"           // projection_3dgut_forward
#include "IntersectTile.cuh"           // do_intersect_tile_generic
#include "RasterizationMomentsFwd.cuh" // rasterize_moments_3dgut_fwd

namespace meshing {

namespace {
template<typename T> T* dmalloc_copy(const T* host, size_t n) {
    T* p = nullptr;
    CHECK_DEVICE_ERROR(cudaMalloc(&p, n * sizeof(T)));
    if (host) CHECK_DEVICE_ERROR(cudaMemcpy(p, host, n * sizeof(T), cudaMemcpyHostToDevice));
    else CHECK_DEVICE_ERROR(cudaMemset(p, 0, n * sizeof(T)));
    return p;
}
inline TorchTensorView tv(const float* p, std::initializer_list<int64_t> shape) {
    return TorchTensorView((uint64_t)p, 4, shape);
}
} // namespace


struct RenderContext {
    int N = 0;
    int C = 0;
    std::vector<int> Ws, Hs;     // per-camera image size
    int Wmax = 0, Hmax = 0;      // max over cameras (scratch-buffer sizing)
    std::string model;

    // raw (un-activated) splat params on device
    float* d_means = nullptr;
    float* d_quats = nullptr;
    float* d_logsc = nullptr;
    float* d_logit = nullptr;   // [N] -> viewed as [N,1]
    float* d_fdc   = nullptr;

    // camera intrinsics on device (all cameras)
    float* d_viewmats = nullptr; // [C*16]
    float* d_intrins  = nullptr; // [C*4]
    float* d_dist     = nullptr; // [C*10]

    DeviceVector<float> radii;   // pool-backed scratch, size N

    int carve_k = 1;

    ~RenderContext() {
        for (void* p : {(void*)d_means,(void*)d_quats,(void*)d_logsc,(void*)d_logit,
                        (void*)d_fdc,(void*)d_viewmats,(void*)d_intrins,(void*)d_dist})
            if (p) cudaFree(p);
    }

    // in_splats in the projection's WorldBuffer layout (6 entries). f_sh is a
    // dummy [N,3] view (never read: we project with sh_degree=0 -> DC only).
    // DeviceTensorFloatND uses a channel-last view: the trailing 1 is the
    // per-scalar channel, so [N,K] is passed as the view shape {N,K,1}.
    std::vector<DeviceTensorFloatND> in_splats() const {
        auto nd = [](float* p, int64_t rows, int64_t cols) {
            return DeviceTensorFloatND(TorchTensorView((uint64_t)p, 4, {rows, cols, 1}));
        };
        return {
            nd(d_means, N, 3),
            nd(d_quats, N, 4),
            nd(d_logsc, N, 3),
            nd(d_logit, N, 1),
            nd(d_fdc,   N, 3),
            nd(d_fdc,   N, 3),  // dummy f_sh (unused at sh_degree=0)
        };
    }
};


RenderContext* render_context_create(
    const float* means, const float* quats, const float* log_scales,
    const float* logit_opac, const float* features_dc, int num_splats,
    const float* viewmats, const float* intrins, const float* dist,
    int num_cameras, const int* widths, const int* heights,
    const std::string& camera_model, int carve_k
) {
    RenderContext* ctx = new RenderContext();
    ctx->N = num_splats;
    ctx->C = num_cameras;
    ctx->Ws.assign(widths, widths + num_cameras);
    ctx->Hs.assign(heights, heights + num_cameras);
    for (int c = 0; c < num_cameras; ++c) {
        ctx->Wmax = std::max(ctx->Wmax, ctx->Ws[c]);
        ctx->Hmax = std::max(ctx->Hmax, ctx->Hs[c]);
    }
    ctx->model = camera_model;
    ctx->carve_k = (carve_k < 1) ? 1 : carve_k;

    const size_t N = (size_t)num_splats;
    ctx->d_means = dmalloc_copy(means, N * 3);
    ctx->d_quats = dmalloc_copy(quats, N * 4);
    ctx->d_logsc = dmalloc_copy(log_scales, N * 3);
    ctx->d_logit = dmalloc_copy(logit_opac, N);
    ctx->d_fdc   = dmalloc_copy(features_dc, N * 3);

    ctx->d_viewmats = dmalloc_copy(viewmats, (size_t)num_cameras * 16);
    ctx->d_intrins  = dmalloc_copy(intrins, (size_t)num_cameras * 4);
    ctx->d_dist     = dmalloc_copy(dist, (size_t)num_cameras * 10);  // null -> zeros

    ctx->radii.resize("meshing.render.radii", N);
    return ctx;
}

// ---------------------------------------------------------------------------
// Shared device helpers
// ---------------------------------------------------------------------------
namespace {

// Bilinear-sample the [H,W] float3 image at continuous pixel (u,v), where pixel
// (j,i) center is at (j+0.5, i+0.5) (matching the rasterizer).
__device__ __forceinline__ float3 bilinear3(
    const float3* img, int W, int H, float u, float v
) {
    float sx = u - 0.5f, sy = v - 0.5f;
    int j0 = (int)floorf(sx), i0 = (int)floorf(sy);
    float fx = sx - j0, fy = sy - i0;
    int j1 = j0 + 1, i1 = i0 + 1;
    j0 = min(max(j0, 0), W - 1); j1 = min(max(j1, 0), W - 1);
    i0 = min(max(i0, 0), H - 1); i1 = min(max(i1, 0), H - 1);
    float3 a = img[i0 * W + j0], b = img[i0 * W + j1];
    float3 c = img[i1 * W + j0], d = img[i1 * W + j1];
    float w00 = (1 - fx) * (1 - fy), w10 = fx * (1 - fy);
    float w01 = (1 - fx) * fy,       w11 = fx * fy;
    return make_float3(
        w00*a.x + w10*b.x + w01*c.x + w11*d.x,
        w00*a.y + w10*b.y + w01*c.y + w11*d.y,
        w00*a.z + w10*b.z + w01*c.z + w11*d.z);
    // return make_float3(
    //     fmaxf(fmaxf(a.x, b.x), fmaxf(c.x, d.x)),
    //     fmaxf(fmaxf(a.y, b.y), fmaxf(c.y, d.y)),
    //     fmaxf(fmaxf(a.z, b.z), fmaxf(c.z, d.z))
    // );
}

// Project p (world) into the camera. Returns false if behind the camera or out
// of frame; else fills pixel (u,v) and the along-ray depth z = |p_cam| (matches
// the rasterizer's evaluate_color depth metric). PINHOLE only for now.
__device__ __forceinline__ bool project_point(
    const float* viewmat, const float* intrin, int W, int H,
    float px, float py, float pz, float& u, float& v, float& z
) {
    float cx_ = viewmat[0]*px + viewmat[1]*py + viewmat[2]*pz + viewmat[3];
    float cy_ = viewmat[4]*px + viewmat[5]*py + viewmat[6]*pz + viewmat[7];
    float cz_ = viewmat[8]*px + viewmat[9]*py + viewmat[10]*pz + viewmat[11];
    if (cz_ <= 1e-6f) return false;
    u = intrin[0] * (cx_ / cz_) + intrin[2];
    v = intrin[1] * (cy_ / cz_) + intrin[3];
    if (u < 0.0f || u >= (float)W || v < 0.0f || v >= (float)H) return false;
    z = sqrtf(cx_*cx_ + cy_*cy_ + cz_*cz_);
    return true;
}

// Reconstruct occupancy 1-T(z) from a sampled moment
// Returns false (ABSTAIN) when the pixel's ray never reached tau_hi -- i.e. no
// confident surface (sky / grazing / too thin) -- so it does not carve.
__device__ __forceinline__ bool occ_from_moment(
    float3 m, float z, float& occ
) {
    float m0 = m.x, c0 = m.y, c1 = m.z;
    if (!(m0 > 0.0f))
        return false;
    occ = fminf(m0, c0 + c1 * z);
    return true;
}

} // namespace


// Render one camera: projection -> intersect -> moment (+ optional rgb) raster.
static void render_one(RenderContext* ctx, int cam_idx,
                       float3* d_moments, float3* d_rgb) {
    if (cam_idx < 0 || cam_idx >= ctx->C)
        throw std::runtime_error("render_one: cam_idx out of range");

    const uint32_t W = (uint32_t)ctx->Ws[cam_idx], H = (uint32_t)ctx->Hs[cam_idx];
    std::vector<DeviceTensorFloatND> in_splats = ctx->in_splats();

    // per-camera views (single image, I=1)
    TorchTensorView viewmats = tv(ctx->d_viewmats + (size_t)cam_idx * 16, {1, 4, 4});
    TorchTensorView intrins  = tv(ctx->d_intrins  + (size_t)cam_idx * 4,  {1, 4});
    TorchTensorView dist     = tv(ctx->d_dist     + (size_t)cam_idx * 10, {1, 10});

    ctx->radii.zero();  // projection accumulates via atomicMax

    // --- projection (3DGUT, sh_degree = 0 -> DC color only) ---
    auto [aabb_2d, depths_2d, splats_s] = projection_3dgut_forward(
        (int64_t)ctx->N, /*max_sh_degree=*/0, in_splats,
        viewmats, intrins, W, H, ctx->model, dist,
        ctx->radii,
        std::nullopt, std::nullopt, /*num_sh_buffer=*/0, /*sh_value_bits=*/32,
        /*sh_bounds_stride=*/0);

    // --- tile intersection (ellipse mode: conic = splats_s[0], opac = [1]) ---
    DeviceTensorFloatND aabb_nd(aabb_2d);
    DeviceTensorFloatND depths_nd(depths_2d);
    DeviceTensorFloatND proj_conic = splats_s[0];
    DeviceTensorFloatND proj_opac  = splats_s[1];
    auto [isect_ids, flatten_ids, tile_offsets] = do_intersect_tile_generic(
        aabb_nd, depths_nd, nullptr, &proj_conic, &proj_opac,
        /*I=*/1, intrins, W, H, nullptr);

    // --- moment (+ rgb) rasterization ---
    rasterize_moments_3dgut_fwd(
        (int64_t)ctx->N, in_splats, splats_s, DeviceVector<int32_t>(),
        viewmats, intrins, ctx->model, dist,
        aabb_2d, W, H, tile_offsets, flatten_ids,
        d_moments, d_rgb);
}

void render_camera_moments(RenderContext* ctx, int cam_idx, void* d_moments) {
    render_one(ctx, cam_idx, reinterpret_cast<float3*>(d_moments), nullptr);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
}

void render_context_destroy(RenderContext* ctx) { delete ctx; }
int render_context_width(const RenderContext* ctx, int cam_idx) { return ctx->Ws[cam_idx]; }
int render_context_height(const RenderContext* ctx, int cam_idx) { return ctx->Hs[cam_idx]; }
int render_context_num_cameras(const RenderContext* ctx) { return ctx->C; }


// ---------------------------------------------------------------------------
// Occupancy sampling (k-th-smallest over the cameras that see the point)
// ---------------------------------------------------------------------------
namespace {

__global__ void sample_occ_kernel(
    const float* __restrict__ xyz, int n,
    const float* __restrict__ viewmat, const float* __restrict__ intrin,
    const float3* __restrict__ moments, int W, int H,
    int k,
    float* __restrict__ occ_kmin, int* __restrict__ cnt
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float u, v, z;
    if (!project_point(viewmat, intrin, W, H, xyz[3*i], xyz[3*i+1], xyz[3*i+2], u, v, z))
        return;
    float occ;
    if (!occ_from_moment(bilinear3(moments, W, H, u, v), z, occ))
        return;  // abstain
    // insert occ into the k smallest seen so far (ascending, +inf padded)
    float* arr = occ_kmin + (size_t)i * k;
    if (occ < arr[k-1]) {
        int p = k - 1;
        while (p > 0 && arr[p-1] > occ) { arr[p] = arr[p-1]; --p; }
        arr[p] = occ;
    }
    cnt[i] += 1;
}

__global__ void finalize_occ_kernel(
    int n, const float* __restrict__ occ_kmin, const int* __restrict__ cnt,
    int k, float* __restrict__ occ
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int c = cnt[i];
    if (c <= 0) { occ[i] = 0.0f; return; }      // unseen -> free space
    int idx = min(k, c) - 1;                     // k-th smallest (or all if fewer)
    occ[i] = occ_kmin[(size_t)i * k + idx];
}

} // namespace

void render_evaluate_occupancy(
    RenderContext* ctx, const int* cam_indices, int num_cams,
    const float* d_xyz, int n, float* d_occ
) {
    if (n <= 0) return;
    const int k = ctx->carve_k;
    const size_t npix = (size_t)ctx->Wmax * ctx->Hmax;
    float3* d_moments = nullptr;
    float* d_occ_kmin = nullptr;
    int* d_cnt = nullptr;
    CHECK_DEVICE_ERROR(cudaMalloc(&d_moments, npix * sizeof(float3)));
    CHECK_DEVICE_ERROR(cudaMalloc(&d_occ_kmin, (size_t)n * k * sizeof(float)));
    CHECK_DEVICE_ERROR(cudaMalloc(&d_cnt, (size_t)n * sizeof(int)));
    {
        std::vector<float> big((size_t)n * k, 1e30f);
        CHECK_DEVICE_ERROR(cudaMemcpy(d_occ_kmin, big.data(),
            (size_t)n * k * sizeof(float), cudaMemcpyHostToDevice));
        CHECK_DEVICE_ERROR(cudaMemset(d_cnt, 0, (size_t)n * sizeof(int)));
    }

    const int TPB = 256;
    int blocks = (n + TPB - 1) / TPB;
    for (int ci = 0; ci < num_cams; ++ci) {
        int cam = cam_indices[ci];
        render_one(ctx, cam, d_moments, nullptr);
        sample_occ_kernel<<<blocks, TPB>>>(
            d_xyz, n,
            ctx->d_viewmats + (size_t)cam * 16, ctx->d_intrins + (size_t)cam * 4,
            d_moments, ctx->Ws[cam], ctx->Hs[cam], k,
            d_occ_kmin, d_cnt);
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    finalize_occ_kernel<<<blocks, TPB>>>(n, d_occ_kmin, d_cnt, k, d_occ);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());

    cudaFree(d_moments); cudaFree(d_occ_kmin); cudaFree(d_cnt);
}


// ---------------------------------------------------------------------------
// Color sampling (rendered DC color weighted by transmittance until the point)
// ---------------------------------------------------------------------------
namespace {

__global__ void sample_color_kernel(
    const float* __restrict__ xyz, int n,
    const float* __restrict__ viewmat, const float* __restrict__ intrin,
    const float3* __restrict__ moments, const float3* __restrict__ rgb,
    int W, int H,
    float3* __restrict__ num, float* __restrict__ den
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float u, v, z;
    if (!project_point(viewmat, intrin, W, H, xyz[3*i], xyz[3*i+1], xyz[3*i+2], u, v, z))
        return;
    float occ;
    if (!occ_from_moment(bilinear3(moments, W, H, u, v), z, occ))
        return;  // abstain (no confident surface from this view)
    float w = 1.0f - occ;            // transmittance until the point
    if (w <= 1e-4f) return;
    float3 c = bilinear3(rgb, W, H, u, v);
    num[i].x += w * c.x; num[i].y += w * c.y; num[i].z += w * c.z;
    den[i] += w;
}

__global__ void finalize_color_kernel(
    int n, const float3* __restrict__ num, const float* __restrict__ den,
    float* __restrict__ rgb
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float d = den[i];
    if (d > 1e-6f) {
        rgb[3*i+0] = num[i].x / d;
        rgb[3*i+1] = num[i].y / d;
        rgb[3*i+2] = num[i].z / d;
    } else {
        rgb[3*i+0] = -1.0f; rgb[3*i+1] = -1.0f; rgb[3*i+2] = -1.0f;  // fallback
    }
}

} // namespace

void render_evaluate_color(
    RenderContext* ctx, const int* cam_indices, int num_cams,
    const float* d_xyz, int n, float* d_rgb
) {
    if (n <= 0) return;
    const size_t npix = (size_t)ctx->Wmax * ctx->Hmax;
    float3* d_moments = nullptr; float3* d_rgbimg = nullptr;
    float3* d_num = nullptr; float* d_den = nullptr;
    CHECK_DEVICE_ERROR(cudaMalloc(&d_moments, npix * sizeof(float3)));
    CHECK_DEVICE_ERROR(cudaMalloc(&d_rgbimg, npix * sizeof(float3)));
    CHECK_DEVICE_ERROR(cudaMalloc(&d_num, (size_t)n * sizeof(float3)));
    CHECK_DEVICE_ERROR(cudaMalloc(&d_den, (size_t)n * sizeof(float)));
    CHECK_DEVICE_ERROR(cudaMemset(d_num, 0, (size_t)n * sizeof(float3)));
    CHECK_DEVICE_ERROR(cudaMemset(d_den, 0, (size_t)n * sizeof(float)));

    const int TPB = 256;
    int blocks = (n + TPB - 1) / TPB;
    for (int ci = 0; ci < num_cams; ++ci) {
        int cam = cam_indices[ci];
        render_one(ctx, cam, d_moments, d_rgbimg);
        sample_color_kernel<<<blocks, TPB>>>(
            d_xyz, n,
            ctx->d_viewmats + (size_t)cam * 16, ctx->d_intrins + (size_t)cam * 4,
            d_moments, d_rgbimg, ctx->Ws[cam], ctx->Hs[cam],
            d_num, d_den);
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    finalize_color_kernel<<<blocks, TPB>>>(n, d_num, d_den, d_rgb);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());

    cudaFree(d_moments); cudaFree(d_rgbimg); cudaFree(d_num); cudaFree(d_den);
}

} // namespace meshing
