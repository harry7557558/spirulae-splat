/*
 * Meshing.cu
 *
 * GPU side of the 3DGS surface meshing pipeline (see Meshing.h):
 *   - Gaussian activation (quat/scale/opacity) and principal-axis extraction.
 *   - 7-points-per-Gaussian point cloud sampling.
 *   - A uniform spatial grid over the Gaussians' support (count / scan / sort /
 *     range pattern, mirroring IntersectTile.cu).
 *   - The occupancy field: a static (density aggregation) variant and a dataset
 *     (per-camera closed-form max density along the point->camera segment,
 *     minimized over cameras) variant.
 *   - Per-edge crossing-point bisection.
 */

#include "Meshing.h"

#include "Common.cuh"

#include <cub/cub.cuh>

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

namespace meshing {

// ---------------------------------------------------------------------------
// Small device vector helpers (Common.cuh already provides +,-,*,/ on floatN).
// ---------------------------------------------------------------------------
__host__ __device__ __forceinline__ float dot3(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
__device__ __forceinline__ int imin(int a, int b) { return a < b ? a : b; }
__device__ __forceinline__ int imax(int a, int b) { return a > b ? a : b; }
__device__ __forceinline__ int iclamp(int x, int a, int b) {
    return x < a ? a : (x > b ? b : x);
}

// per-axis cap on how many grid cells a single Gaussian may stamp into
inline constexpr int kStampAxisCap = 31;

// ---------------------------------------------------------------------------
// Device-side view of the whole scene (Gaussians + grid + cameras). POD, passed
// by value into kernels.
// ---------------------------------------------------------------------------
struct GpuScene {
    // Gaussian SoA, indexed by ORIGINAL splat id [0, N)
    const float3* __restrict__ mean;
    const float3* __restrict__ ax0;   // unit principal axes (columns of R)
    const float3* __restrict__ ax1;
    const float3* __restrict__ ax2;
    const float3* __restrict__ invs2; // 1 / sigma^2 per axis
    const float*  __restrict__ opac;  // sigmoid opacity
    const float*  __restrict__ k2;    // 2 ln(opac / ALPHA): Mahalanobis cutoff^2

    // uniform grid
    float3 origin;
    float  cell;
    float  inv_cell;
    int3   res;
    const int* __restrict__ cellStart; // [numCells]
    const int* __restrict__ cellEnd;   // [numCells]
    const int* __restrict__ sortedGauss; // [numPairs] -> original splat id

    // cameras
    const float3* __restrict__ campos;
    int num_cameras;

    float iso;
};

__device__ __forceinline__ int3 cell_of(const GpuScene& s, float3 p) {
    return make_int3(
        (int)floorf((p.x - s.origin.x) * s.inv_cell),
        (int)floorf((p.y - s.origin.y) * s.inv_cell),
        (int)floorf((p.z - s.origin.z) * s.inv_cell));
}
__device__ __forceinline__ bool cell_in(const GpuScene& s, int3 c) {
    return c.x >= 0 && c.x < s.res.x &&
           c.y >= 0 && c.y < s.res.y &&
           c.z >= 0 && c.z < s.res.z;
}
__device__ __forceinline__ int cell_lin(const GpuScene& s, int3 c) {
    return (c.z * s.res.y + c.y) * s.res.x + c.x;
}

// Mahalanobis quadratic form delta^T Sigma^-1 delta for Gaussian g.
__device__ __forceinline__ float quad_form(const GpuScene& s, int g, float3 d) {
    float u0 = dot3(d, s.ax0[g]);
    float u1 = dot3(d, s.ax1[g]);
    float u2 = dot3(d, s.ax2[g]);
    float3 is = s.invs2[g];
    return u0 * u0 * is.x + u1 * u1 * is.y + u2 * u2 * is.z;
}

// Point density of Gaussian g at p (0 if below the ALPHA boundary).
__device__ __forceinline__ float density_at(const GpuScene& s, int g, float3 p) {
    float q = quad_form(s, g, p - s.mean[g]);
    if (q >= s.k2[g]) return 0.0f;
    return s.opac[g] * __expf(-0.5f * q);
}

// ---- static occupancy: aggregate densities of Gaussians in p's cell --------
__device__ float occ_static(const GpuScene& s, float3 p) {
    int3 c = cell_of(s, p);
    if (!cell_in(s, c)) return 0.0f;
    int lin = cell_lin(s, c);
    int a = s.cellStart[lin], b = s.cellEnd[lin];
    float S = 0.0f;  // sum of log(1 - density)
    for (int j = a; j < b; ++j) {
        int g = s.sortedGauss[j];
        float d = density_at(s, g, p);
        if (d > ALPHA_THRESHOLD)
            S += __logf(1.0f - fminf(d, 0.999f));
    }
    return 1.0f - __expf(S);
}

// ---- dataset occupancy -----------------------------------------------------
// Closed-form max density of Gaussian g along the segment o + t*dir, t in [0,1],
// returning the density and (via tstar) the closest-approach parameter.
__device__ __forceinline__ float seg_max_density(
    const GpuScene& s, int g, float3 o, float3 dir, float& tstar
) {
    float3 e = o - s.mean[g];
    float e0 = dot3(e, s.ax0[g]), e1 = dot3(e, s.ax1[g]), e2 = dot3(e, s.ax2[g]);
    float d0 = dot3(dir, s.ax0[g]), d1 = dot3(dir, s.ax1[g]), d2 = dot3(dir, s.ax2[g]);
    float3 is = s.invs2[g];
    float A = d0 * d0 * is.x + d1 * d1 * is.y + d2 * d2 * is.z;       // dir^T M dir
    float B = e0 * d0 * is.x + e1 * d1 * is.y + e2 * d2 * is.z;       // e^T M dir
    float C = e0 * e0 * is.x + e1 * e1 * is.y + e2 * e2 * is.z;       // e^T M e
    float t = (A > 1e-20f) ? (-B / A) : 0.0f;
    t = fminf(fmaxf(t, 0.0f), 1.0f);
    tstar = t;
    float q = A * t * t + 2.0f * B * t + C;
    if (q >= s.k2[g]) return 0.0f;
    return s.opac[g] * __expf(-0.5f * q);
}

// Aggregate (1 - prod(1-d)) of Gaussians whose closest-approach point to the
// segment p->cam lies in the cells the segment traverses. Each Gaussian is
// counted exactly once (in the cell containing its closest-approach point),
// which avoids double counting from multi-cell stamping. Returns sum log(1-d).
__device__ float seg_log_transmittance(const GpuScene& s, float3 p, float3 cam) {
    float3 dir = cam - p;
    float S = 0.0f;

    int3 c = cell_of(s, p);
    // step / tMax / tDelta for a 3D-DDA over a cubic grid, t in [0,1]
    int   stepx = dir.x > 0 ? 1 : -1, stepy = dir.y > 0 ? 1 : -1, stepz = dir.z > 0 ? 1 : -1;
    float invdx = fabsf(dir.x) > 1e-20f ? 1.0f / fabsf(dir.x) : 1e20f;
    float invdy = fabsf(dir.y) > 1e-20f ? 1.0f / fabsf(dir.y) : 1e20f;
    float invdz = fabsf(dir.z) > 1e-20f ? 1.0f / fabsf(dir.z) : 1e20f;
    float tDeltax = s.cell * invdx, tDeltay = s.cell * invdy, tDeltaz = s.cell * invdz;
    // distance (in t) from p to the first cell boundary in each axis
    float bx = s.origin.x + (float)(c.x + (stepx > 0 ? 1 : 0)) * s.cell;
    float by = s.origin.y + (float)(c.y + (stepy > 0 ? 1 : 0)) * s.cell;
    float bz = s.origin.z + (float)(c.z + (stepz > 0 ? 1 : 0)) * s.cell;
    float tMaxx = (bx - p.x) * (stepx > 0 ? invdx : -invdx);
    float tMaxy = (by - p.y) * (stepy > 0 ? invdy : -invdy);
    float tMaxz = (bz - p.z) * (stepz > 0 ? invdz : -invdz);

    const int max_steps = s.res.x + s.res.y + s.res.z + 3;
    float t = 0.0f;
    for (int it = 0; it < max_steps; ++it) {
        if (cell_in(s, c)) {
            int lin = cell_lin(s, c);
            int a = s.cellStart[lin], b = s.cellEnd[lin];
            for (int j = a; j < b; ++j) {
                int g = s.sortedGauss[j];
                float tstar;
                float d = seg_max_density(s, g, p, dir, tstar);
                if (d <= ALPHA_THRESHOLD) continue;
                float3 pstar = p + dir * tstar;
                int3 cs = cell_of(s, pstar);
                if (cs.x == c.x && cs.y == c.y && cs.z == c.z)
                    S += __logf(1.0f - fminf(d, 0.999f));
            }
        }
        // advance to the next cell
        if (tMaxx < tMaxy && tMaxx < tMaxz) {
            t = tMaxx; c.x += stepx; tMaxx += tDeltax;
        } else if (tMaxy < tMaxz) {
            t = tMaxy; c.y += stepy; tMaxy += tDeltay;
        } else {
            t = tMaxz; c.z += stepz; tMaxz += tDeltaz;
        }
        if (t > 1.0f) break;
    }
    return S;
}

__device__ float occ_dynamic(const GpuScene& s, float3 p) {
    float occ_min = 1.0f;
    for (int ci = 0; ci < s.num_cameras; ++ci) {
        float S = seg_log_transmittance(s, p, s.campos[ci]);
        float occ = 1.0f - __expf(S);
        occ_min = fminf(occ_min, occ);
    }
    return occ_min;
}

__device__ __forceinline__ float occ_eval(const GpuScene& s, float3 p, bool dynamic) {
    return dynamic ? occ_dynamic(s, p) : occ_static(s, p);
}

// ---------------------------------------------------------------------------
// Kernels
// ---------------------------------------------------------------------------
__global__ void activate_kernel(
    int N,
    const float* __restrict__ means,  const float* __restrict__ quats,
    const float* __restrict__ logsc,  const float* __restrict__ logit,
    float3* mean, float3* ax0, float3* ax1, float3* ax2,
    float3* invs2, float* opac, float* radius, float* k2, int* valid
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    mean[i] = make_float3(means[3 * i], means[3 * i + 1], means[3 * i + 2]);

    float w = quats[4 * i], x = quats[4 * i + 1], y = quats[4 * i + 2], z = quats[4 * i + 3];
    float n = rsqrtf(fmaxf(w * w + x * x + y * y + z * z, 1e-20f));
    w *= n; x *= n; y *= n; z *= n;
    // columns of the rotation matrix = principal axes
    ax0[i] = make_float3(1 - 2 * (y * y + z * z), 2 * (x * y + w * z),     2 * (x * z - w * y));
    ax1[i] = make_float3(2 * (x * y - w * z),     1 - 2 * (x * x + z * z), 2 * (y * z + w * x));
    ax2[i] = make_float3(2 * (x * z + w * y),     2 * (y * z - w * x),     1 - 2 * (x * x + y * y));

    float sx = __expf(logsc[3 * i]), sy = __expf(logsc[3 * i + 1]), sz = __expf(logsc[3 * i + 2]);
    invs2[i] = make_float3(1.0f / (sx * sx), 1.0f / (sy * sy), 1.0f / (sz * sz));

    float op = 1.0f / (1.0f + __expf(-logit[i]));
    opac[i] = op;
    if (op > ALPHA_THRESHOLD) {
        float kk2 = 2.0f * logf(op / ALPHA_THRESHOLD);
        k2[i] = kk2;
        float k = sqrtf(kk2);
        radius[i] = k * fmaxf(sx, fmaxf(sy, sz));
        valid[i] = 1;
    } else {
        k2[i] = 0.0f;
        radius[i] = 0.0f;
        valid[i] = 0;
    }
}

// 7 points per kept Gaussian: center + (+/- k*sigma) along each principal axis.
__global__ void pointcloud_kernel(
    int num_kept, const int* __restrict__ kept,
    const float3* mean, const float3* ax0, const float3* ax1, const float3* ax2,
    const float3* invs2, const float* k2,
    double* out  // [num_kept*7*3]
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_kept) return;
    int g = kept[i];
    float k = sqrtf(k2[g]);
    float3 m = mean[g];
    float3 is = invs2[g];
    float s0 = rsqrtf(is.x), s1 = rsqrtf(is.y), s2 = rsqrtf(is.z);
    float3 o0 = ax0[g] * (k * s0);
    float3 o1 = ax1[g] * (k * s1);
    float3 o2 = ax2[g] * (k * s2);
    float3 pts[7] = { m, m + o0, m - o0, m + o1, m - o1, m + o2, m - o2 };
    long base = (long)i * 7 * 3;
    for (int j = 0; j < 7; ++j) {
        out[base + 3 * j + 0] = (double)pts[j].x;
        out[base + 3 * j + 1] = (double)pts[j].y;
        out[base + 3 * j + 2] = (double)pts[j].z;
    }
}

// grid: per kept Gaussian count the cells its bounding-sphere AABB stamps into
__global__ void grid_count_kernel(
    int num_kept, const int* __restrict__ kept,
    const float3* mean, const float* radius,
    float3 origin, float inv_cell, int3 res,
    int* counts
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_kept) return;
    int g = kept[i];
    float3 m = mean[g];
    float r = radius[g];
    int lx = iclamp((int)floorf((m.x - r - origin.x) * inv_cell), 0, res.x - 1);
    int hx = iclamp((int)floorf((m.x + r - origin.x) * inv_cell), 0, res.x - 1);
    int ly = iclamp((int)floorf((m.y - r - origin.y) * inv_cell), 0, res.y - 1);
    int hy = iclamp((int)floorf((m.y + r - origin.y) * inv_cell), 0, res.y - 1);
    int lz = iclamp((int)floorf((m.z - r - origin.z) * inv_cell), 0, res.z - 1);
    int hz = iclamp((int)floorf((m.z + r - origin.z) * inv_cell), 0, res.z - 1);
    hx = imin(hx, lx + kStampAxisCap - 1);
    hy = imin(hy, ly + kStampAxisCap - 1);
    hz = imin(hz, lz + kStampAxisCap - 1);
    counts[i] = (hx - lx + 1) * (hy - ly + 1) * (hz - lz + 1);
}

__global__ void grid_fill_kernel(
    int num_kept, const int* __restrict__ kept,
    const float3* mean, const float* radius,
    float3 origin, float inv_cell, int3 res,
    const int64_t* __restrict__ cum,  // inclusive scan of counts
    int* keys, int* vals
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_kept) return;
    int g = kept[i];
    float3 m = mean[g];
    float r = radius[g];
    int lx = iclamp((int)floorf((m.x - r - origin.x) * inv_cell), 0, res.x - 1);
    int hx = iclamp((int)floorf((m.x + r - origin.x) * inv_cell), 0, res.x - 1);
    int ly = iclamp((int)floorf((m.y - r - origin.y) * inv_cell), 0, res.y - 1);
    int hy = iclamp((int)floorf((m.y + r - origin.y) * inv_cell), 0, res.y - 1);
    int lz = iclamp((int)floorf((m.z - r - origin.z) * inv_cell), 0, res.z - 1);
    int hz = iclamp((int)floorf((m.z + r - origin.z) * inv_cell), 0, res.z - 1);
    hx = imin(hx, lx + kStampAxisCap - 1);
    hy = imin(hy, ly + kStampAxisCap - 1);
    hz = imin(hz, lz + kStampAxisCap - 1);
    long pos = (long)(cum[i] - (int64_t)((hx - lx + 1) * (hy - ly + 1) * (hz - lz + 1)));
    for (int cz = lz; cz <= hz; ++cz)
        for (int cy = ly; cy <= hy; ++cy)
            for (int cx = lx; cx <= hx; ++cx) {
                keys[pos] = (cz * res.y + cy) * res.x + cx;
                vals[pos] = g;
                ++pos;
            }
}

__global__ void grid_range_kernel(
    long n_pairs, const int* __restrict__ keys, int* cellStart, int* cellEnd
) {
    long j = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n_pairs) return;
    int c = keys[j];
    if (j == 0 || keys[j - 1] != c) cellStart[c] = (int)j;
    if (j == n_pairs - 1 || keys[j + 1] != c) cellEnd[c] = (int)(j + 1);
}

__global__ void occ_kernel(
    GpuScene s, const double* __restrict__ pts, int n, float* occ, int dynamic
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float3 p = make_float3((float)pts[3 * i], (float)pts[3 * i + 1], (float)pts[3 * i + 2]);
    occ[i] = occ_eval(s, p, dynamic != 0);
}

__global__ void bisect_kernel(
    GpuScene s, const double* __restrict__ cloud,
    const int* __restrict__ ea, const int* __restrict__ eb,
    const float* __restrict__ oa, const float* __restrict__ ob,
    int n, int iters, int dynamic, double* out
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int a = ea[i], b = eb[i];
    float3 pa = make_float3((float)cloud[3 * a], (float)cloud[3 * a + 1], (float)cloud[3 * a + 2]);
    float3 pb = make_float3((float)cloud[3 * b], (float)cloud[3 * b + 1], (float)cloud[3 * b + 2]);
    bool a_in = oa[i] >= s.iso;  // is endpoint a on the "occupied" side
    for (int it = 0; it < iters; ++it) {
        float3 mid = (pa + pb) * 0.5f;
        float om = occ_eval(s, mid, dynamic != 0);
        if ((om >= s.iso) == a_in) pa = mid;
        else                       pb = mid;
    }
    float3 mid = (pa + pb) * 0.5f;
    out[3 * i + 0] = (double)mid.x;
    out[3 * i + 1] = (double)mid.y;
    out[3 * i + 2] = (double)mid.z;
}

// ---------------------------------------------------------------------------
// Impl
// ---------------------------------------------------------------------------
template<typename T> static T* dmalloc(size_t n) {
    T* p = nullptr;
    CHECK_DEVICE_ERROR(cudaMalloc(&p, n * sizeof(T)));
    return p;
}

struct OccupancyEvaluator::Impl {
    MeshingConfig cfg;
    int N = 0;

    // Gaussian SoA (device)
    float3 *mean = nullptr, *ax0 = nullptr, *ax1 = nullptr, *ax2 = nullptr, *invs2 = nullptr;
    float  *opac = nullptr, *radius = nullptr, *k2 = nullptr;
    int    *valid = nullptr;

    int *kept = nullptr;     // [num_kept] -> original id (device)
    int num_kept = 0;

    // grid
    GpuScene scene{};
    int *cellStart = nullptr, *cellEnd = nullptr, *sortedGauss = nullptr;
    int *keysBuf = nullptr;  // owned sort buffers
    long n_pairs = 0;

    float3 *campos = nullptr;
    int num_cameras = 0;

    ~Impl() {
        for (void* p : {(void*)mean,(void*)ax0,(void*)ax1,(void*)ax2,(void*)invs2,
                        (void*)opac,(void*)radius,(void*)k2,(void*)valid,(void*)kept,
                        (void*)cellStart,(void*)cellEnd,(void*)sortedGauss,(void*)campos})
            if (p) cudaFree(p);
    }

    GpuScene make_scene() const {
        GpuScene s = scene;
        s.mean = mean; s.ax0 = ax0; s.ax1 = ax1; s.ax2 = ax2; s.invs2 = invs2;
        s.opac = opac; s.k2 = k2;
        s.cellStart = cellStart; s.cellEnd = cellEnd; s.sortedGauss = sortedGauss;
        s.campos = campos; s.num_cameras = num_cameras;
        s.iso = cfg.iso;
        return s;
    }
};

OccupancyEvaluator::OccupancyEvaluator(
    const float* means, const float* quats,
    const float* log_scales, const float* logit_opac, int num_splats,
    const float* cam_pos, int num_cameras,
    const MeshingConfig& cfg
) {
    impl_ = new Impl();
    impl_->cfg = cfg;
    impl_->N = num_splats;
    const int N = num_splats;

    // ---- upload raw params + activate ----
    float *d_means = dmalloc<float>((size_t)N * 3);
    float *d_quats = dmalloc<float>((size_t)N * 4);
    float *d_logsc = dmalloc<float>((size_t)N * 3);
    float *d_logit = dmalloc<float>((size_t)N);
    cudaMemcpy(d_means, means, sizeof(float) * N * 3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_quats, quats, sizeof(float) * N * 4, cudaMemcpyHostToDevice);
    cudaMemcpy(d_logsc, log_scales, sizeof(float) * N * 3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_logit, logit_opac, sizeof(float) * N, cudaMemcpyHostToDevice);

    impl_->mean = dmalloc<float3>(N); impl_->ax0 = dmalloc<float3>(N);
    impl_->ax1 = dmalloc<float3>(N);  impl_->ax2 = dmalloc<float3>(N);
    impl_->invs2 = dmalloc<float3>(N); impl_->opac = dmalloc<float>(N);
    impl_->radius = dmalloc<float>(N); impl_->k2 = dmalloc<float>(N);
    impl_->valid = dmalloc<int>(N);

    activate_kernel<<<_LAUNCH_ARGS_1D(N, 256)>>>(
        N, d_means, d_quats, d_logsc, d_logit,
        impl_->mean, impl_->ax0, impl_->ax1, impl_->ax2, impl_->invs2,
        impl_->opac, impl_->radius, impl_->k2, impl_->valid);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaFree(d_means); cudaFree(d_quats); cudaFree(d_logsc); cudaFree(d_logit);

    // ---- download mean/radius/valid; build kept list + scene bbox on host ----
    std::vector<float> h_mean(N * 3), h_radius(N);
    std::vector<int> h_valid(N);
    cudaMemcpy(h_mean.data(), impl_->mean, sizeof(float) * N * 3, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_radius.data(), impl_->radius, sizeof(float) * N, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_valid.data(), impl_->valid, sizeof(int) * N, cudaMemcpyDeviceToHost);

    std::vector<int> kept;
    kept.reserve(N);
    float bmin[3] = {1e30f, 1e30f, 1e30f}, bmax[3] = {-1e30f, -1e30f, -1e30f};
    double rsum = 0.0;
    float rmax = 0.0f;
    for (int i = 0; i < N; ++i) {
        if (!h_valid[i]) continue;
        kept.push_back(i);
        float r = h_radius[i];
        rsum += r; rmax = std::max(rmax, r);
        for (int a = 0; a < 3; ++a) {
            float m = h_mean[3 * i + a];
            bmin[a] = std::min(bmin[a], m - r);
            bmax[a] = std::max(bmax[a], m + r);
        }
    }
    impl_->num_kept = (int)kept.size();
    num_kept_ = impl_->num_kept;
    num_points_ = impl_->num_kept * 7;

    if (impl_->num_kept == 0)
        throw std::runtime_error("meshing: no Gaussians above the opacity threshold");

    impl_->kept = dmalloc<int>(impl_->num_kept);
    cudaMemcpy(impl_->kept, kept.data(), sizeof(int) * impl_->num_kept, cudaMemcpyHostToDevice);

    // ---- grid params ----
    float mean_r = (float)(rsum / impl_->num_kept);
    float cell = std::max(cfg.grid_cell_factor * mean_r, 1e-8f);
    // pad bbox by one max-radius so every stamp lands inside the grid
    for (int a = 0; a < 3; ++a) { bmin[a] -= rmax * 0.5f; bmax[a] += rmax * 0.5f; }
    int3 res;
    int* resp = &res.x;
    float ext[3];
    for (int a = 0; a < 3; ++a) {
        ext[a] = std::max(bmax[a] - bmin[a], cell);
        int r = (int)std::ceil(ext[a] / cell);
        resp[a] = std::min(std::max(r, 1), cfg.max_grid_res);
    }
    // if any axis was capped, grow the cell so the grid still covers the bbox
    float cell_need = cell;
    for (int a = 0; a < 3; ++a) cell_need = std::max(cell_need, ext[a] / resp[a]);
    cell = cell_need;
    for (int a = 0; a < 3; ++a)
        resp[a] = std::min(std::max((int)std::ceil(ext[a] / cell), 1), cfg.max_grid_res);

    impl_->scene.origin = make_float3(bmin[0], bmin[1], bmin[2]);
    impl_->scene.cell = cell;
    impl_->scene.inv_cell = 1.0f / cell;
    impl_->scene.res = res;
    long numCells = (long)res.x * res.y * res.z;

    if (cfg.verbose)
        printf("[meshing] %d/%d Gaussians kept, grid %dx%dx%d (cell=%.4g), %ld points\n",
               impl_->num_kept, N, res.x, res.y, res.z, cell, (long)num_points_);

    // ---- build grid: count -> scan -> fill -> sort -> ranges ----
    int* d_counts = dmalloc<int>(impl_->num_kept);
    grid_count_kernel<<<_LAUNCH_ARGS_1D(impl_->num_kept, 256)>>>(
        impl_->num_kept, impl_->kept, impl_->mean, impl_->radius,
        impl_->scene.origin, impl_->scene.inv_cell, res, d_counts);
    CHECK_DEVICE_ERROR(cudaGetLastError());

    int64_t* d_cum = dmalloc<int64_t>(impl_->num_kept);
    // counts are int; widen during scan via a transform-less InclusiveSum into int64
    // (CUB handles mixed types by output type)
    CUB_WRAPPER(cub::DeviceScan::InclusiveSum, d_counts, d_cum, impl_->num_kept);
    CHECK_DEVICE_ERROR(cudaGetLastError());
    int64_t n_pairs = 0;
    cudaMemcpy(&n_pairs, d_cum + (impl_->num_kept - 1), sizeof(int64_t), cudaMemcpyDeviceToHost);
    impl_->n_pairs = n_pairs;

    int *keys_a = dmalloc<int>(n_pairs), *keys_b = dmalloc<int>(n_pairs);
    int *vals_a = dmalloc<int>(n_pairs), *vals_b = dmalloc<int>(n_pairs);
    grid_fill_kernel<<<_LAUNCH_ARGS_1D(impl_->num_kept, 256)>>>(
        impl_->num_kept, impl_->kept, impl_->mean, impl_->radius,
        impl_->scene.origin, impl_->scene.inv_cell, res, d_cum, keys_a, vals_a);
    CHECK_DEVICE_ERROR(cudaGetLastError());
    cudaFree(d_counts); cudaFree(d_cum);

    int nbits = 0; while ((1L << nbits) <= numCells) ++nbits;
    cub::DoubleBuffer<int> dk(keys_a, keys_b), dv(vals_a, vals_b);
    CUB_WRAPPER(cub::DeviceRadixSort::SortPairs, dk, dv, (int)n_pairs, 0, nbits);
    CHECK_DEVICE_ERROR(cudaGetLastError());
    int* sortedKeys = dk.Current();
    impl_->sortedGauss = dv.Current();
    // free the unused (non-current) buffers
    if (dk.Current() == keys_a) cudaFree(keys_b); else cudaFree(keys_a);
    if (dv.Current() == vals_a) cudaFree(vals_b); else cudaFree(vals_a);
    impl_->keysBuf = sortedKeys;  // keep alive (range kernel reads it); freed in dtor? no -> free now after ranges

    impl_->cellStart = dmalloc<int>(numCells);
    impl_->cellEnd   = dmalloc<int>(numCells);
    cudaMemset(impl_->cellStart, 0, sizeof(int) * numCells);
    cudaMemset(impl_->cellEnd,   0, sizeof(int) * numCells);
    grid_range_kernel<<<_LAUNCH_ARGS_1D(n_pairs, 256)>>>(
        n_pairs, sortedKeys, impl_->cellStart, impl_->cellEnd);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaFree(sortedKeys);
    impl_->keysBuf = nullptr;

    // ---- cameras (subsample to max_cameras) ----
    if (cam_pos != nullptr && num_cameras > 0) {
        int stride = std::max(1, (num_cameras + cfg.max_cameras - 1) / cfg.max_cameras);
        std::vector<float> cams;
        for (int c = 0; c < num_cameras; c += stride) {
            cams.push_back(cam_pos[3 * c]);
            cams.push_back(cam_pos[3 * c + 1]);
            cams.push_back(cam_pos[3 * c + 2]);
        }
        impl_->num_cameras = (int)(cams.size() / 3);
        impl_->campos = dmalloc<float3>(impl_->num_cameras);
        cudaMemcpy(impl_->campos, cams.data(), sizeof(float) * cams.size(), cudaMemcpyHostToDevice);
        if (cfg.verbose)
            printf("[meshing] using %d/%d cameras (dataset occupancy)\n",
                   impl_->num_cameras, num_cameras);
    }
}

OccupancyEvaluator::~OccupancyEvaluator() { delete impl_; }

void OccupancyEvaluator::generate_point_cloud(std::vector<double>& xyz_out) {
    long n = (long)impl_->num_kept * 7;
    xyz_out.resize((size_t)n * 3);
    double* d_out = dmalloc<double>((size_t)n * 3);
    pointcloud_kernel<<<_LAUNCH_ARGS_1D(impl_->num_kept, 256)>>>(
        impl_->num_kept, impl_->kept,
        impl_->mean, impl_->ax0, impl_->ax1, impl_->ax2, impl_->invs2, impl_->k2,
        d_out);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaMemcpy(xyz_out.data(), d_out, sizeof(double) * n * 3, cudaMemcpyDeviceToHost);
    cudaFree(d_out);
}

void OccupancyEvaluator::evaluate(const double* xyz, int n, float* occ_out) {
    if (n <= 0) return;
    double* d_xyz = dmalloc<double>((size_t)n * 3);
    float*  d_occ = dmalloc<float>(n);
    cudaMemcpy(d_xyz, xyz, sizeof(double) * n * 3, cudaMemcpyHostToDevice);
    GpuScene s = impl_->make_scene();
    int dynamic = impl_->num_cameras > 0 ? 1 : 0;
    occ_kernel<<<_LAUNCH_ARGS_1D(n, 128)>>>(s, d_xyz, n, d_occ, dynamic);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaMemcpy(occ_out, d_occ, sizeof(float) * n, cudaMemcpyDeviceToHost);
    cudaFree(d_xyz); cudaFree(d_occ);
}

void OccupancyEvaluator::bisect_edges(
    const double* cloud_xyz,
    const int32_t* edge_a, const int32_t* edge_b,
    const float* occ_a, const float* occ_b,
    int n_edges, double* xyz_out
) {
    if (n_edges <= 0) return;
    long ncloud = (long)num_points_;
    double* d_cloud = dmalloc<double>((size_t)ncloud * 3);
    cudaMemcpy(d_cloud, cloud_xyz, sizeof(double) * ncloud * 3, cudaMemcpyHostToDevice);
    int* d_ea = dmalloc<int>(n_edges); int* d_eb = dmalloc<int>(n_edges);
    float* d_oa = dmalloc<float>(n_edges); float* d_ob = dmalloc<float>(n_edges);
    double* d_out = dmalloc<double>((size_t)n_edges * 3);
    cudaMemcpy(d_ea, edge_a, sizeof(int) * n_edges, cudaMemcpyHostToDevice);
    cudaMemcpy(d_eb, edge_b, sizeof(int) * n_edges, cudaMemcpyHostToDevice);
    cudaMemcpy(d_oa, occ_a, sizeof(float) * n_edges, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ob, occ_b, sizeof(float) * n_edges, cudaMemcpyHostToDevice);
    GpuScene s = impl_->make_scene();
    int dynamic = impl_->num_cameras > 0 ? 1 : 0;
    bisect_kernel<<<_LAUNCH_ARGS_1D(n_edges, 128)>>>(
        s, d_cloud, d_ea, d_eb, d_oa, d_ob, n_edges, impl_->cfg.bisection_iters, dynamic, d_out);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaMemcpy(xyz_out, d_out, sizeof(double) * n_edges * 3, cudaMemcpyDeviceToHost);
    cudaFree(d_cloud); cudaFree(d_ea); cudaFree(d_eb);
    cudaFree(d_oa); cudaFree(d_ob); cudaFree(d_out);
}

} // namespace meshing
