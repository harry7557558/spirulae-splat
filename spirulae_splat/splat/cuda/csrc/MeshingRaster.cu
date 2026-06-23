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
#include <cub/cub.cuh>
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

#ifdef __CUDACC__
// {persp,fisheye,equisolid}_proj_nav -- match the (possibly distorted) camera
// model the occupancy/color images were rendered with. slang.cuh + FixedArray /
// Matrix come in via Common.cuh above.
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
#endif

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
// the rasterizer's evaluate_color depth metric). Uses the same projection +
// distortion model the occupancy/color images were rendered with, so the sample
// lands on the pixel the splats were actually rasterized to.
__device__ __forceinline__ bool project_point(
    const float* viewmat, const float* intrin, const float* dist,
    CameraModelType model, int W, int H,
    float px, float py, float pz, float& u, float& v, float& z
) {
    float cx_ = viewmat[0]*px + viewmat[1]*py + viewmat[2]*pz + viewmat[3];
    float cy_ = viewmat[4]*px + viewmat[5]*py + viewmat[6]*pz + viewmat[7];
    float cz_ = viewmat[8]*px + viewmat[9]*py + viewmat[10]*pz + viewmat[11];

    float3 p_cam = make_float3(cx_, cy_, cz_);
    float4 intr = make_float4(intrin[0], intrin[1], intrin[2], intrin[3]);
    CameraDistortionCoeffs dist_coeffs;
    #pragma unroll
    for (int t = 0; t < 10; ++t) dist_coeffs[t] = dist ? dist[t] : 0.0f;

    // proj_nav handles the behind-camera / invalid-distortion cases and returns
    // pixel-space uv (already scaled by fx,fy and offset by cx,cy).
    float2 uv;
    bool valid =
        (model == CameraModelType::FISHEYE)
            ? SlangProjectionUtils::fisheye_proj_nav(p_cam, intr, dist_coeffs, &uv) :
        (model == CameraModelType::EQUISOLID)
            ? SlangProjectionUtils::equisolid_proj_nav(p_cam, intr, dist_coeffs, &uv) :
            SlangProjectionUtils::persp_proj_nav(p_cam, intr, dist_coeffs, &uv);
    if (!valid) return false;
    u = uv.x; v = uv.y;
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
    #if 1
    // occ = fminf(m0, c0 + c1 * z);
    occ = fminf(m0, 1.0f - __expf(-powf(fmaxf(c0 + c1 * z, 0.0f), 2.0f)));
    #else
    if (c0 == 0.0f && c1 == 0.0f)
        return false;
    // occ = fminf(m0, 1.0f - expf(c0 + c1 * z));
    occ = fmaxf(fminf(m0, c0 + c1 * z), 0.0f);
    #endif
    return true;
}

// Bilinear-sample occupancy at continuous pixel (u,v): evaluate the moment ->
// occupancy at each of the 4 corner pixels first, then bilinearly blend the
// resulting occupancies. Corners that abstain (occ_from_moment false) are
// dropped and the weights renormalized over the valid corners; abstains only if
// every corner abstains. This avoids blending moments across a surface edge
// (sky m0=0 mixed into an object) before the nonlinear evaluation.
__device__ __forceinline__ bool occ_bilinear(
    const float3* img, int W, int H, float u, float v, float z, float& occ
) {
    float sx = u - 0.5f, sy = v - 0.5f;
    int j0 = (int)floorf(sx), i0 = (int)floorf(sy);
    float fx = sx - j0, fy = sy - i0;
    int j1 = j0 + 1, i1 = i0 + 1;
    j0 = min(max(j0, 0), W - 1); j1 = min(max(j1, 0), W - 1);
    i0 = min(max(i0, 0), H - 1); i1 = min(max(i1, 0), H - 1);
    const float3 m[4] = { img[i0*W+j0], img[i0*W+j1], img[i1*W+j0], img[i1*W+j1] };
    const float w[4] = { (1-fx)*(1-fy), fx*(1-fy), (1-fx)*fy, fx*fy };
    float acc = 0.0f, wsum = 0.0f;
    for (int t = 0; t < 4; ++t) {
        float o;
        if (occ_from_moment(m[t], z, o)) { acc += w[t] * o; wsum += w[t]; }
    }
    if (wsum <= 0.0f) return false;
    occ = (acc + 1e-6f) / (wsum + 1e-6f);
    occ = fminf(fmaxf(occ, 0.0f), 1.0f);
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
    const float* __restrict__ dist, CameraModelType model,
    const float3* __restrict__ moments, int W, int H,
    int k,
    float* __restrict__ occ_kmin, int* __restrict__ cnt
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float u, v, z;
    if (!project_point(viewmat, intrin, dist, model, W, H,
                       xyz[3*i], xyz[3*i+1], xyz[3*i+2], u, v, z))
        return;
    float occ;
    if (!occ_bilinear(moments, W, H, u, v, z, occ))
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
            ctx->d_dist + (size_t)cam * 10, cmt(ctx->model),
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
    const float* __restrict__ dist, CameraModelType model,
    const float3* __restrict__ moments, const float3* __restrict__ rgb,
    int W, int H,
    float3* __restrict__ num, float* __restrict__ den
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float u, v, z;
    if (!project_point(viewmat, intrin, dist, model, W, H,
                       xyz[3*i], xyz[3*i+1], xyz[3*i+2], u, v, z))
        return;
    float occ;
    if (!occ_bilinear(moments, W, H, u, v, z, occ))
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
#if 0
    if (d > 1e-6f) {
        rgb[3*i+0] = num[i].x / d;
        rgb[3*i+1] = num[i].y / d;
        rgb[3*i+2] = num[i].z / d;
    } else {
        // picked up by colorize_fallback_kernel
        rgb[3*i+0] = -1.0f; rgb[3*i+1] = -1.0f; rgb[3*i+2] = -1.0f;  // fallback
    }
#else
    constexpr float e = 1e-8f;
    d = fmaxf(d, 0.0f) + e;
    rgb[3*i+0] = (num[i].x + e) / d;
    rgb[3*i+1] = (num[i].y + e) / d;
    rgb[3*i+2] = (num[i].z + e) / d;
#endif
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
            ctx->d_dist + (size_t)cam * 10, cmt(ctx->model),
            d_moments, d_rgbimg, ctx->Ws[cam], ctx->Hs[cam],
            d_num, d_den);
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
    finalize_color_kernel<<<blocks, TPB>>>(n, d_num, d_den, d_rgb);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());

    cudaFree(d_moments); cudaFree(d_rgbimg); cudaFree(d_num); cudaFree(d_den);
}


// ---------------------------------------------------------------------------
// Visibility cull: drop mesh vertices seen by no camera.
//
// A vertex is "seen" by a camera when it projects inside that camera's frame
// (project_point) AND the segment from the vertex to that camera's center is not
// blocked by a mesh triangle that does NOT contain the vertex. Occlusion is
// tested against an LBVH built over the mesh triangles (Karras radix tree,
// mirroring the Gaussian LBVH in Meshing.cu). A vertex is kept iff >= 1 camera
// sees it. Camera centers are recovered from each view matrix (C = -R^T t).
// ---------------------------------------------------------------------------
namespace {

__device__ __forceinline__ float3 fmin3(float3 a, float3 b) {
    return make_float3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}
__device__ __forceinline__ float3 fmax3(float3 a, float3 b) {
    return make_float3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}

// 21-bit-per-axis Morton interleave (same as Meshing.cu's LBVH build).
__device__ __forceinline__ uint64_t expand_bits_21(uint64_t x) {
    x &= 0x1fffffULL;
    x = (x | (x << 32)) & 0x1f00000000ffffULL;
    x = (x | (x << 16)) & 0x1f0000ff0000ffULL;
    x = (x | (x << 8))  & 0x100f00f00f00f00fULL;
    x = (x | (x << 4))  & 0x10c30c30c30c30c3ULL;
    x = (x | (x << 2))  & 0x1249249249249249ULL;
    return x;
}
__device__ __forceinline__ uint64_t morton3D(float3 c /* in [0,1] */) {
    const float s = (float)((1u << 21) - 1);
    uint64_t xi = (uint64_t)fminf(fmaxf(c.x * s, 0.0f), s);
    uint64_t yi = (uint64_t)fminf(fmaxf(c.y * s, 0.0f), s);
    uint64_t zi = (uint64_t)fminf(fmaxf(c.z * s, 0.0f), s);
    return (expand_bits_21(xi) << 2) | (expand_bits_21(yi) << 1) | expand_bits_21(zi);
}

__device__ __forceinline__ float atomicMinF(float* addr, float val) {
    int* a = (int*)addr; int old = *a, assumed;
    do { assumed = old;
         if (__int_as_float(assumed) <= val) break;
         old = atomicCAS(a, assumed, __float_as_int(val));
    } while (assumed != old);
    return __int_as_float(old);
}
__device__ __forceinline__ float atomicMaxF(float* addr, float val) {
    int* a = (int*)addr; int old = *a, assumed;
    do { assumed = old;
         if (__int_as_float(assumed) >= val) break;
         old = atomicCAS(a, assumed, __float_as_int(val));
    } while (assumed != old);
    return __int_as_float(old);
}

// Slab test: does segment o + t*dir, t in [0,1], touch the AABB?
__device__ __forceinline__ bool seg_aabb(float3 o, float3 dir, float3 bmin, float3 bmax) {
    float t0 = 0.0f, t1 = 1.0f;
    #pragma unroll
    for (int a = 0; a < 3; ++a) {
        float oa = a==0?o.x:(a==1?o.y:o.z);
        float da = a==0?dir.x:(a==1?dir.y:dir.z);
        float lo = a==0?bmin.x:(a==1?bmin.y:bmin.z);
        float hi = a==0?bmax.x:(a==1?bmax.y:bmax.z);
        if (fabsf(da) < 1e-20f) { if (oa < lo || oa > hi) return false; }
        else {
            float inv = 1.0f / da;
            float ta = (lo - oa) * inv, tb = (hi - oa) * inv;
            if (ta > tb) { float tmp = ta; ta = tb; tb = tmp; }
            t0 = fmaxf(t0, ta); t1 = fminf(t1, tb);
            if (t0 > t1) return false;
        }
    }
    return true;
}

// Moeller-Trumbore: segment o + t*dir intersects triangle (a,b,c)? Reports t.
__device__ __forceinline__ bool ray_tri(
    float3 o, float3 dir, float3 a, float3 b, float3 c, float& tout
) {
    float3 e1 = b - a, e2 = c - a;
    float3 pv = cross(dir, e2);
    float det = dot(e1, pv);
    if (fabsf(det) < 1e-20f) return false;          // ray parallel to triangle
    float inv = 1.0f / det;
    float3 tv = o - a;
    float u = dot(tv, pv) * inv;
    if (u < -1e-6f || u > 1.0f + 1e-6f) return false;
    float3 qv = cross(tv, e1);
    float v = dot(dir, qv) * inv;
    if (v < -1e-6f || u + v > 1.0f + 1e-6f) return false;
    tout = dot(e2, qv) * inv;
    return true;
}

// Does triangle `tid` (skipped if it contains vertex `vid`) block the open
// segment o->cam? Parametric margins exclude intersections hugging the
// endpoints: the near margin suppresses self-occlusion acne from triangles
// adjacent to (but not containing) the vertex, biasing the borderline toward
// "visible" so the surface is not over-culled.
__device__ __forceinline__ bool tri_hits(
    float3 o, float3 dir, int tid, int vid,
    const float* __restrict__ verts, const int* __restrict__ faces
) {
    int ia = faces[3*tid], ib = faces[3*tid+1], ic = faces[3*tid+2];
    if (ia == vid || ib == vid || ic == vid) return false;   // connected -> skip
    float3 a = make_float3(verts[3*ia], verts[3*ia+1], verts[3*ia+2]);
    float3 b = make_float3(verts[3*ib], verts[3*ib+1], verts[3*ib+2]);
    float3 c = make_float3(verts[3*ic], verts[3*ic+1], verts[3*ic+2]);
    float t;
    return ray_tri(o, dir, a, b, c, t) && t > 1e-3f && t < 1.0f - 1e-4f;
}

inline constexpr int kCullStack = 64;

__device__ bool seg_blocked(
    float3 o, float3 cam, int vid,
    const float* __restrict__ verts, const int* __restrict__ faces, int nf,
    const float3* __restrict__ leafMin, const float3* __restrict__ leafMax,
    const int2* __restrict__ internal, const float3* __restrict__ nodeAABB
) {
    float3 dir = cam - o;
    if (nf < 2) {                       // no internal nodes: scan the (<=1) tris
        for (int t = 0; t < nf; ++t)
            if (tri_hits(o, dir, t, vid, verts, faces)) return true;
        return false;
    }
    int stack[kCullStack]; int sp = 0;
    if (!seg_aabb(o, dir, nodeAABB[0], nodeAABB[1])) return false;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        int2 ch = internal[ni];
        #pragma unroll
        for (int cc = 0; cc < 2; ++cc) {
            int child = cc == 0 ? ch.x : ch.y;
            if (child < 0) {
                int tid = ~child;
                if (seg_aabb(o, dir, leafMin[tid], leafMax[tid]) &&
                    tri_hits(o, dir, tid, vid, verts, faces))
                    return true;
            } else if (seg_aabb(o, dir, nodeAABB[2*child], nodeAABB[2*child+1])) {
                if (sp < kCullStack) stack[sp++] = child;
            }
        }
    }
    return false;
}

// Per-triangle leaf AABB + centroid Morton code (and an iota id).
__global__ void tri_prep_kernel(
    int nf, const int* __restrict__ faces, const float* __restrict__ verts,
    float3 bmin, float3 inv_ext,
    float3* leafMin, float3* leafMax, uint64_t* morton, int* iota
) {
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= nf) return;
    int ia = faces[3*t], ib = faces[3*t+1], ic = faces[3*t+2];
    float3 a = make_float3(verts[3*ia], verts[3*ia+1], verts[3*ia+2]);
    float3 b = make_float3(verts[3*ib], verts[3*ib+1], verts[3*ib+2]);
    float3 c = make_float3(verts[3*ic], verts[3*ic+1], verts[3*ic+2]);
    leafMin[t] = fmin3(fmin3(a, b), c);
    leafMax[t] = fmax3(fmax3(a, b), c);
    float3 ctr = (a + b + c) * (1.0f / 3.0f);
    float3 rc = make_float3((ctr.x - bmin.x) * inv_ext.x,
                            (ctr.y - bmin.y) * inv_ext.y,
                            (ctr.z - bmin.z) * inv_ext.z);
    morton[t] = morton3D(rc);
    iota[t] = t;
}

// Karras 2012 single-level radix tree (identical to Meshing.cu's bvh build).
__global__ void tri_internal_kernel(
    int n, const uint64_t* __restrict__ morton, const int* __restrict__ argsort,
    int2* __restrict__ internal, int* __restrict__ parent
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n - 1) return;

    #define delta(a, b) \
        (((b) < 0 || (b) >= n) ? -1 : \
         (morton[a] == morton[b] ? 64 + __clz((a) ^ (b)) \
                                 : __clzll(morton[a] ^ morton[b])))

    int d = delta(i, i+1) - delta(i, i-1);
    d = d > 0 ? 1 : (d < 0 ? -1 : 0);
    int delta_min = delta(i, i-d);
    int lmax = 2;
    while (delta(i, i + lmax*d) > delta_min) lmax <<= 1;
    int l = 0;
    for (int t = lmax >> 1; t >= 1; t >>= 1)
        if (delta(i, i + (l+t)*d) > delta_min) l += t;
    int j = i + l*d;
    int delta_node = delta(i, j);
    int sp_ = 0;
    for (int tf = 2, t; (t = (l + tf - 1) / tf) >= 1; tf <<= 1)
        if (delta(i, i + (sp_ + t)*d) > delta_node) sp_ += t;
    int gamma = i + sp_*d + min(d, 0);

    int left  = (min(i,j) == gamma)     ? ~argsort[gamma]     : gamma;
    int right = (max(i,j) == gamma + 1) ? ~argsort[gamma+1]   : gamma + 1;
    internal[i] = make_int2(left, right);
    if (left  >= 0) atomicMax(&parent[left],  i);
    if (right >= 0) atomicMax(&parent[right], i);
    #undef delta
}

__global__ void tri_initaabb_kernel(int n_internal, float3* nodeAABB) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_internal) return;
    nodeAABB[2*i]   = make_float3(1e30f, 1e30f, 1e30f);
    nodeAABB[2*i+1] = make_float3(-1e30f, -1e30f, -1e30f);
}

__global__ void tri_aabb_kernel(
    int n,
    const int2* __restrict__ internal, const int* __restrict__ parent,
    const float3* __restrict__ leafMin, const float3* __restrict__ leafMax,
    float3* __restrict__ nodeAABB
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n - 1) return;
    int2 ch = internal[i];
    if (ch.x >= 0 && ch.y >= 0) return;  // no leaf child -> filled from below

    float3 bmin = make_float3(1e30f, 1e30f, 1e30f);
    float3 bmax = make_float3(-1e30f, -1e30f, -1e30f);
    if (ch.x < 0) { bmin = fmin3(bmin, leafMin[~ch.x]); bmax = fmax3(bmax, leafMax[~ch.x]); }
    if (ch.y < 0) { bmin = fmin3(bmin, leafMin[~ch.y]); bmax = fmax3(bmax, leafMax[~ch.y]); }

    int node = i;
    do {
        bool covered =
            atomicMinF(&nodeAABB[2*node].x,   bmin.x) <= bmin.x &
            atomicMinF(&nodeAABB[2*node].y,   bmin.y) <= bmin.y &
            atomicMinF(&nodeAABB[2*node].z,   bmin.z) <= bmin.z &
            atomicMaxF(&nodeAABB[2*node+1].x, bmax.x) >= bmax.x &
            atomicMaxF(&nodeAABB[2*node+1].y, bmax.y) >= bmax.y &
            atomicMaxF(&nodeAABB[2*node+1].z, bmax.z) >= bmax.z;
        if (covered) break;
        node = parent[node];
    } while (node >= 0);
}

__global__ void cull_kernel(
    const float* __restrict__ verts, int nv,
    const int* __restrict__ faces, int nf,
    const float* __restrict__ viewmats, const float* __restrict__ intrins,
    const float* __restrict__ dist,
    const int* __restrict__ Ws, const int* __restrict__ Hs,
    CameraModelType model, int C,
    const float3* __restrict__ leafMin, const float3* __restrict__ leafMax,
    const int2* __restrict__ internal, const float3* __restrict__ nodeAABB,
    unsigned char* __restrict__ visible
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nv) return;
    float3 p = make_float3(verts[3*i], verts[3*i+1], verts[3*i+2]);
    for (int c = 0; c < C; ++c) {
        const float* vm = viewmats + (size_t)c * 16;
        float u, v, z;
        if (!project_point(vm, intrins + (size_t)c*4, dist + (size_t)c*10,
                           model, Ws[c], Hs[c], p.x, p.y, p.z, u, v, z))
            continue;                              // out of frame / behind camera
        // camera center in world: C = -R^T t (viewmat row-major world->cam 4x4)
        float3 cam = make_float3(
            -(vm[0]*vm[3] + vm[4]*vm[7] + vm[8]*vm[11]),
            -(vm[1]*vm[3] + vm[5]*vm[7] + vm[9]*vm[11]),
            -(vm[2]*vm[3] + vm[6]*vm[7] + vm[10]*vm[11]));
        if (!seg_blocked(p, cam, i, verts, faces, nf,
                         leafMin, leafMax, internal, nodeAABB)) {
            visible[i] = 1;                        // seen by this camera -> keep
            return;
        }
    }
    visible[i] = 0;                                // seen by no camera
}

} // namespace

void render_cull_unseen_vertices(
    RenderContext* ctx, const float* verts, int nv,
    const int* faces, int nf, unsigned char* visible
) {
    if (nv <= 0) return;

    float* d_verts = dmalloc_copy(verts, (size_t)nv * 3);
    int*   d_faces = (nf > 0) ? dmalloc_copy(faces, (size_t)nf * 3) : nullptr;
    int*   d_W = dmalloc_copy(ctx->Ws.data(), (size_t)ctx->C);
    int*   d_H = dmalloc_copy(ctx->Hs.data(), (size_t)ctx->C);
    unsigned char* d_vis = nullptr;
    CHECK_DEVICE_ERROR(cudaMalloc(&d_vis, (size_t)nv));

    // ---- triangle LBVH over the mesh faces ----
    float3 *d_leafMin = nullptr, *d_leafMax = nullptr, *d_nodeAABB = nullptr;
    int2   *d_internal = nullptr;
    if (nf >= 1) {
        // scene bbox over the mesh vertices (host) for Morton normalization. The
        // mesh is a bounded surface, so a plain linear normalization suffices
        // (no distant outliers like the Gaussian LBVH's coordinate remap).
        float bb[6] = {1e30f,1e30f,1e30f,-1e30f,-1e30f,-1e30f};
        for (int i = 0; i < nv; ++i)
            for (int a = 0; a < 3; ++a) {
                float x = verts[3*i+a];
                bb[a] = std::min(bb[a], x); bb[3+a] = std::max(bb[3+a], x);
            }
        float3 bmin = make_float3(bb[0], bb[1], bb[2]);
        float3 inv_ext = make_float3(1.0f/std::max(bb[3]-bb[0], 1e-12f),
                                     1.0f/std::max(bb[4]-bb[1], 1e-12f),
                                     1.0f/std::max(bb[5]-bb[2], 1e-12f));
        CHECK_DEVICE_ERROR(cudaMalloc(&d_leafMin, (size_t)nf * sizeof(float3)));
        CHECK_DEVICE_ERROR(cudaMalloc(&d_leafMax, (size_t)nf * sizeof(float3)));
        uint64_t* d_morton = nullptr; int* d_iota = nullptr;
        CHECK_DEVICE_ERROR(cudaMalloc(&d_morton, (size_t)nf * sizeof(uint64_t)));
        CHECK_DEVICE_ERROR(cudaMalloc(&d_iota,   (size_t)nf * sizeof(int)));
        tri_prep_kernel<<<_LAUNCH_ARGS_1D(nf, 256)>>>(
            nf, d_faces, d_verts, bmin, inv_ext, d_leafMin, d_leafMax, d_morton, d_iota);
        CHECK_DEVICE_ERROR(cudaGetLastError());

        if (nf >= 2) {
            uint64_t* d_morton_s = nullptr; int* d_argsort = nullptr;
            CHECK_DEVICE_ERROR(cudaMalloc(&d_morton_s, (size_t)nf * sizeof(uint64_t)));
            CHECK_DEVICE_ERROR(cudaMalloc(&d_argsort,  (size_t)nf * sizeof(int)));
            CUB_WRAPPER(cub::DeviceRadixSort::SortPairs,
                d_morton, d_morton_s, d_iota, d_argsort, nf);
            CHECK_DEVICE_ERROR(cudaGetLastError());

            CHECK_DEVICE_ERROR(cudaMalloc(&d_internal, (size_t)(nf - 1) * sizeof(int2)));
            int* d_parent = nullptr;
            CHECK_DEVICE_ERROR(cudaMalloc(&d_parent, (size_t)(nf - 1) * sizeof(int)));
            cudaMemset(d_parent, 0xff, (size_t)(nf - 1) * sizeof(int));
            tri_internal_kernel<<<_LAUNCH_ARGS_1D(nf - 1, 256)>>>(
                nf, d_morton_s, d_argsort, d_internal, d_parent);
            CHECK_DEVICE_ERROR(cudaGetLastError());

            CHECK_DEVICE_ERROR(cudaMalloc(&d_nodeAABB, (size_t)2 * (nf - 1) * sizeof(float3)));
            tri_initaabb_kernel<<<_LAUNCH_ARGS_1D(nf - 1, 256)>>>(nf - 1, d_nodeAABB);
            tri_aabb_kernel<<<_LAUNCH_ARGS_1D(nf - 1, 256)>>>(
                nf, d_internal, d_parent, d_leafMin, d_leafMax, d_nodeAABB);
            CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
            cudaFree(d_morton_s); cudaFree(d_argsort); cudaFree(d_parent);
        }
        cudaFree(d_morton); cudaFree(d_iota);
    }

    const int TPB = 256;
    int blocks = (nv + TPB - 1) / TPB;
    cull_kernel<<<blocks, TPB>>>(
        d_verts, nv, d_faces, nf,
        ctx->d_viewmats, ctx->d_intrins, ctx->d_dist, d_W, d_H, cmt(ctx->model), ctx->C,
        d_leafMin, d_leafMax, d_internal, d_nodeAABB, d_vis);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    CHECK_DEVICE_ERROR(cudaMemcpy(visible, d_vis, (size_t)nv, cudaMemcpyDeviceToHost));

    cudaFree(d_verts); if (d_faces) cudaFree(d_faces);
    cudaFree(d_W); cudaFree(d_H); cudaFree(d_vis);
    if (d_leafMin)  cudaFree(d_leafMin);
    if (d_leafMax)  cudaFree(d_leafMax);
    if (d_internal) cudaFree(d_internal);
    if (d_nodeAABB) cudaFree(d_nodeAABB);
}

} // namespace meshing
