/*
 * Meshing.cu
 *
 * GPU side of the 3DGS surface meshing pipeline (see Meshing.h):
 *   - Gaussian activation (quat/scale/opacity) and principal-axis extraction.
 *   - 7-points-per-Gaussian point cloud sampling.
 *   - An LBVH over the Gaussians' support (Karras radix tree + coordinate
 *     remap, mirroring SplatTileIntersector.cu), so unbounded scenes with
 *     distant splats stay fast -- unlike a uniform grid, whose cell size is
 *     dictated by the farthest splat.
 *   - The occupancy field: a static (density aggregation) variant and a dataset
 *     (per-camera closed-form max density along the point->camera segment,
 *     minimized over cameras) variant. Each Gaussian is a single BVH leaf, so
 *     it is visited exactly once per query -- no multi-cell double counting.
 *   - Per-edge crossing-point bisection.
 */

#include "Meshing.h"
#include "MeshingRaster.cuh"

#include "Common.cuh"

#include <cub/cub.cuh>

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

namespace meshing {

// ---------------------------------------------------------------------------
// Small helpers (Common.cuh already provides +,-,*,/ on floatN).
// ---------------------------------------------------------------------------
__host__ __device__ __forceinline__ float dot3(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
__device__ __forceinline__ float3 fmin3(float3 a, float3 b) {
    return make_float3(fminf(a.x,b.x), fminf(a.y,b.y), fminf(a.z,b.z));
}
__device__ __forceinline__ float3 fmax3(float3 a, float3 b) {
    return make_float3(fmaxf(a.x,b.x), fmaxf(a.y,b.y), fmaxf(a.z,b.z));
}

// Coordinate remap: identity near the origin, ~ x^(1/k) for large |x|. Applied
// only to Morton ordering so a few distant splats don't starve the near-origin
// region of spatial resolution. (Copied from SplatTileIntersector.cu.)
__host__ __device__ __forceinline__ float remap_coord(float x, float rel_scale) {
    constexpr float k = 2.5f;
    return k * sinhf((1.0f / k) * asinhf(x / rel_scale)) * rel_scale;
}

// 21-bit-per-axis Morton interleave.
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

// ---------------------------------------------------------------------------
// Host: pick k representative cameras out of C.
//
// For small k, uniform-interval (file-order) sampling clusters poorly and is
// order-dependent. Instead run k-means (farthest-point seeding + Lloyd) on the
// camera centers and return the MEDOID of each cluster -- the real camera
// nearest the centroid. Medoids (not centroids) matter: averaging cameras that
// surround an object lands near the object center, a useless ray origin.
// Returns the selected positions interleaved (xyz), up to k of them.
// ---------------------------------------------------------------------------
static std::vector<float> select_cameras_kmeans(
    const float* cam, int C, int k, int iters = 25,
    std::vector<int>* out_idx = nullptr
) {
    auto d2 = [&](int a, const float* c) {
        float dx = cam[3*a]-c[0], dy = cam[3*a+1]-c[1], dz = cam[3*a+2]-c[2];
        return dx*dx + dy*dy + dz*dz;
    };
    std::vector<float> cen(3 * k);
    // farthest-point (k-center greedy) seeding -- deterministic, good spread
    std::vector<float> mind(C, 1e30f);
    int pick = 0;
    for (int j = 0; j < k; ++j) {
        cen[3*j] = cam[3*pick]; cen[3*j+1] = cam[3*pick+1]; cen[3*j+2] = cam[3*pick+2];
        float best = -1.0f; int bi = 0;
        for (int c = 0; c < C; ++c) {
            float dd = d2(c, &cen[3*j]);
            if (dd < mind[c]) mind[c] = dd;
            if (mind[c] > best) { best = mind[c]; bi = c; }
        }
        pick = bi;
    }
    // Lloyd iterations
    std::vector<int> assign(C, 0);
    for (int it = 0; it < iters; ++it) {
        bool changed = false;
        for (int c = 0; c < C; ++c) {
            int best = 0; float bd = 1e30f;
            for (int j = 0; j < k; ++j) {
                float dd = d2(c, &cen[3*j]);
                if (dd < bd) { bd = dd; best = j; }
            }
            if (assign[c] != best) { assign[c] = best; changed = true; }
        }
        std::vector<double> sum(3 * k, 0.0);
        std::vector<int> cnt(k, 0);
        for (int c = 0; c < C; ++c) {
            int j = assign[c];
            sum[3*j] += cam[3*c]; sum[3*j+1] += cam[3*c+1]; sum[3*j+2] += cam[3*c+2];
            ++cnt[j];
        }
        for (int j = 0; j < k; ++j)
            if (cnt[j] > 0)
                for (int a = 0; a < 3; ++a) cen[3*j+a] = (float)(sum[3*j+a] / cnt[j]);
        if (!changed) break;
    }
    // medoid per non-empty cluster
    std::vector<int> medoid(k, -1);
    std::vector<float> mbest(k, 1e30f);
    for (int c = 0; c < C; ++c) {
        int j = assign[c];
        float dd = d2(c, &cen[3*j]);
        if (dd < mbest[j]) { mbest[j] = dd; medoid[j] = c; }
    }
    std::vector<float> out;
    for (int j = 0; j < k; ++j)
        if (medoid[j] >= 0) {
            int c = medoid[j];
            out.push_back(cam[3*c]); out.push_back(cam[3*c+1]); out.push_back(cam[3*c+2]);
            if (out_idx) out_idx->push_back(c);
        }
    return out;
}

// ---------------------------------------------------------------------------
// Device-side view of the whole scene (Gaussians + LBVH + cameras). POD.
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
    const float3* __restrict__ gcol;  // base RGB from SH DC, in [0,1]

    // LBVH over kept Gaussians (leaf index = "kept position" kp)
    const int*    __restrict__ kept;        // [num_kept] kp -> original id
    const float3* __restrict__ leafMin;     // [num_kept]
    const float3* __restrict__ leafMax;     // [num_kept]
    const int2*   __restrict__ internal;    // [num_kept-1] child links
    const float3* __restrict__ nodeAABB;    // [2*(num_kept-1)] (min,max) pairs
    int num_kept;

    // cameras
    const float3* __restrict__ campos;
    int num_cameras;

    float iso;
};

// ---- Gaussian density math -------------------------------------------------
__device__ __forceinline__ float quad_form(const GpuScene& s, int g, float3 d) {
    float u0 = dot3(d, s.ax0[g]);
    float u1 = dot3(d, s.ax1[g]);
    float u2 = dot3(d, s.ax2[g]);
    float3 is = s.invs2[g];
    return u0 * u0 * is.x + u1 * u1 * is.y + u2 * u2 * is.z;
}
__device__ __forceinline__ float density_at(const GpuScene& s, int g, float3 p) {
    float q = quad_form(s, g, p - s.mean[g]);
    if (q >= s.k2[g]) return 0.0f;
    return s.opac[g] * __expf(-0.5f * q);
}
// Closed-form max density of Gaussian g along the segment o + t*dir, t in [0,1].
__device__ __forceinline__ float seg_max_density(const GpuScene& s, int g, float3 o, float3 dir) {
    float3 e = o - s.mean[g];
    float e0 = dot3(e, s.ax0[g]), e1 = dot3(e, s.ax1[g]), e2 = dot3(e, s.ax2[g]);
    float d0 = dot3(dir, s.ax0[g]), d1 = dot3(dir, s.ax1[g]), d2 = dot3(dir, s.ax2[g]);
    float3 is = s.invs2[g];
    float A = d0*d0*is.x + d1*d1*is.y + d2*d2*is.z;
    float B = e0*d0*is.x + e1*d1*is.y + e2*d2*is.z;
    float C = e0*e0*is.x + e1*e1*is.y + e2*e2*is.z;
    float t = (A > 1e-20f) ? (-B / A) : 0.0f;
    t = fminf(fmaxf(t, 0.0f), 1.0f);
    float q = A*t*t + 2.0f*B*t + C;
    if (q >= s.k2[g]) return 0.0f;
    return s.opac[g] * __expf(-0.5f * q);
}

// ---- AABB tests ------------------------------------------------------------
__device__ __forceinline__ bool point_in_aabb(float3 p, float3 bmin, float3 bmax) {
    return p.x >= bmin.x && p.x <= bmax.x &&
           p.y >= bmin.y && p.y <= bmax.y &&
           p.z >= bmin.z && p.z <= bmax.z;
}
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

inline constexpr int kStackSize = 64;

// ---- static occupancy: density aggregation of Gaussians containing p -------
__device__ float occ_static(const GpuScene& s, float3 p) {
    float S = 0.0f;  // sum log(1 - density)
    if (s.num_kept == 1) {
        float d = density_at(s, s.kept[0], p);
        return d > ALPHA_THRESHOLD ? d : 0.0f;
    }
    int stack[kStackSize]; int sp = 0;
    if (!point_in_aabb(p, s.nodeAABB[0], s.nodeAABB[1])) return 0.0f;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        int2 ch = s.internal[ni];
        #pragma unroll
        for (int c = 0; c < 2; ++c) {
            int child = c == 0 ? ch.x : ch.y;
            if (child < 0) {
                int kp = ~child;
                if (point_in_aabb(p, s.leafMin[kp], s.leafMax[kp])) {
                    float d = density_at(s, s.kept[kp], p);
                    if (d > ALPHA_THRESHOLD) S += __logf(1.0f - fminf(d, 0.999f));
                }
            } else if (point_in_aabb(p, s.nodeAABB[2*child], s.nodeAABB[2*child+1])) {
                if (sp < kStackSize) stack[sp++] = child;
            }
        }
    }
    return 1.0f - __expf(S);
}

// ---- one camera: sum log(1-d) over Gaussians the segment p->cam crosses -----
__device__ float seg_log_transmittance(const GpuScene& s, float3 p, float3 cam) {
    float3 dir = cam - p;
    float S = 0.0f;
    if (s.num_kept == 1) {
        float d = seg_max_density(s, s.kept[0], p, dir);
        return d > ALPHA_THRESHOLD ? __logf(1.0f - fminf(d, 0.999f)) : 0.0f;
    }
    int stack[kStackSize]; int sp = 0;
    if (!seg_aabb(p, dir, s.nodeAABB[0], s.nodeAABB[1])) return 0.0f;
    stack[sp++] = 0;
    while (sp > 0) {
        int ni = stack[--sp];
        int2 ch = s.internal[ni];
        #pragma unroll
        for (int c = 0; c < 2; ++c) {
            int child = c == 0 ? ch.x : ch.y;
            if (child < 0) {
                int kp = ~child;
                if (seg_aabb(p, dir, s.leafMin[kp], s.leafMax[kp])) {
                    float d = seg_max_density(s, s.kept[kp], p, dir);
                    if (d > ALPHA_THRESHOLD) S += __logf(1.0f - fminf(d, 0.999f));
                }
            } else if (seg_aabb(p, dir, s.nodeAABB[2*child], s.nodeAABB[2*child+1])) {
                if (sp < kStackSize) stack[sp++] = child;
            }
        }
    }
    return S;
}

__device__ float occ_dynamic(const GpuScene& s, float3 p) {
    float occ_min = 1.0f;
    for (int ci = 0; ci < s.num_cameras; ++ci) {
        float S = seg_log_transmittance(s, p, s.campos[ci]);
        occ_min = fminf(occ_min, 1.0f - __expf(S));
    }
    return occ_min;
}

__device__ __forceinline__ float occ_eval(const GpuScene& s, float3 p, bool dynamic) {
    return dynamic ? occ_dynamic(s, p) : occ_static(s, p);
}

// ---------------------------------------------------------------------------
// Vertex coloring from the splats' DC color (gcol).
// ---------------------------------------------------------------------------
// Same as seg_max_density but also reports the (clamped) closest-approach
// parameter, needed to order splats front-to-back along a camera ray.
__device__ __forceinline__ float seg_color_sample(
    const GpuScene& s, int g, float3 o, float3 dir, float& tstar
) {
    float3 e = o - s.mean[g];
    float e0 = dot3(e, s.ax0[g]), e1 = dot3(e, s.ax1[g]), e2 = dot3(e, s.ax2[g]);
    float d0 = dot3(dir, s.ax0[g]), d1 = dot3(dir, s.ax1[g]), d2 = dot3(dir, s.ax2[g]);
    float3 is = s.invs2[g];
    float A = d0*d0*is.x + d1*d1*is.y + d2*d2*is.z;
    float B = e0*d0*is.x + e1*d1*is.y + e2*d2*is.z;
    float C = e0*e0*is.x + e1*e1*is.y + e2*e2*is.z;
    float t = (A > 1e-20f) ? (-B / A) : 0.0f;
    t = fminf(fmaxf(t, 0.0f), 1.0f);
    tstar = t;
    float q = A*t*t + 2.0f*B*t + C;
    if (q >= s.k2[g]) return 0.0f;
    return s.opac[g] * __expf(-0.5f * q);
}

// Density-weighted average color of the splats that contain p. Returns false
// if p is inside no splat's support (e.g. a vertex nudged into a gap by merge).
__device__ bool color_static(const GpuScene& s, float3 p, float3& out) {
    float wsum = 0.0f; float3 csum = make_float3(0.f, 0.f, 0.f);
    if (s.num_kept == 1) {
        float d = density_at(s, s.kept[0], p);
        if (d > 0.0f) { out = s.gcol[s.kept[0]]; return true; }
        return false;
    }
    int stack[kStackSize]; int sp = 0;
    if (!point_in_aabb(p, s.nodeAABB[0], s.nodeAABB[1])) return false;
    stack[sp++] = 0;
    while (sp > 0) {
        int2 ch = s.internal[stack[--sp]];
        #pragma unroll
        for (int c = 0; c < 2; ++c) {
            int child = c == 0 ? ch.x : ch.y;
            if (child < 0) {
                int kp = ~child;
                if (point_in_aabb(p, s.leafMin[kp], s.leafMax[kp])) {
                    int g = s.kept[kp];
                    float d = density_at(s, g, p);
                    if (d > 0.0f) { wsum += d; csum += s.gcol[g] * d; }
                }
            } else if (point_in_aabb(p, s.nodeAABB[2*child], s.nodeAABB[2*child+1])) {
                if (sp < kStackSize) stack[sp++] = child;
            }
        }
    }
    if (wsum > 0.0f) { out = csum * (1.0f / wsum); return true; }
    return false;
}

// Color of the splat whose center is nearest p (last-resort fallback).
__device__ float3 color_nearest(const GpuScene& s, float3 p) {
    if (s.num_kept <= 1) return s.gcol[s.kept[0]];
    float best = 1e30f; int bg = s.kept[0];
    int stack[kStackSize]; int sp = 0; stack[sp++] = 0;
    while (sp > 0) {
        int2 ch = s.internal[stack[--sp]];
        #pragma unroll
        for (int c = 0; c < 2; ++c) {
            int child = c == 0 ? ch.x : ch.y;
            if (child < 0) {
                int g = s.kept[~child];
                float3 d = p - s.mean[g];
                float d2 = dot3(d, d);
                if (d2 < best) { best = d2; bg = g; }
            } else {
                float3 lo = s.nodeAABB[2*child], hi = s.nodeAABB[2*child+1];
                float dx = fmaxf(fmaxf(lo.x - p.x, p.x - hi.x), 0.0f);
                float dy = fmaxf(fmaxf(lo.y - p.y, p.y - hi.y), 0.0f);
                float dz = fmaxf(fmaxf(lo.z - p.z, p.z - hi.z), 0.0f);
                if (dx*dx + dy*dy + dz*dz < best && sp < kStackSize) stack[sp++] = child;
            }
        }
    }
    return s.gcol[bg];
}

inline constexpr int kColorCap = 64;  // max splats composited per camera ray

// Camera-aware color: for each camera, alpha-composite the splats front-to-back
// along the segment camera->p (up to p), then average the per-camera colors
// weighted by accumulated opacity. Falls back to the static color, then to the
// nearest splat, when no camera ray deposits weight.
__device__ float3 color_camera(const GpuScene& s, float3 p) {
    float3 num = make_float3(0.f, 0.f, 0.f);
    float den = 0.0f;

    for (int ci = 0; ci < s.num_cameras; ++ci) {
        float3 cam = s.campos[ci];
        float3 dir = p - cam;             // segment cam->p, t in [0,1]

        // gather hits (keep the CAP closest-to-camera ones)
        float ts[kColorCap], dd[kColorCap]; float3 cc[kColorCap];
        int cnt = 0; float tmax = -1.0f; int tmax_i = 0;
        int stack[kStackSize]; int sp = 0;
        if (s.num_kept >= 2 && seg_aabb(cam, dir, s.nodeAABB[0], s.nodeAABB[1]))
            stack[sp++] = 0;
        else if (s.num_kept == 1) { /* handled below via single leaf */ }
        while (sp > 0) {
            int2 ch = s.internal[stack[--sp]];
            #pragma unroll
            for (int c = 0; c < 2; ++c) {
                int child = c == 0 ? ch.x : ch.y;
                if (child < 0) {
                    int kp = ~child;
                    if (!seg_aabb(cam, dir, s.leafMin[kp], s.leafMax[kp])) continue;
                    int g = s.kept[kp]; float tstar;
                    float d = seg_color_sample(s, g, cam, dir, tstar);
                    if (d <= ALPHA_THRESHOLD) continue;
                    if (cnt < kColorCap) {
                        ts[cnt] = tstar; dd[cnt] = d; cc[cnt] = s.gcol[g];
                        if (tstar > tmax) { tmax = tstar; tmax_i = cnt; }
                        ++cnt;
                    } else if (tstar < tmax) {
                        ts[tmax_i] = tstar; dd[tmax_i] = d; cc[tmax_i] = s.gcol[g];
                        tmax = -1.0f;
                        for (int j = 0; j < kColorCap; ++j)
                            if (ts[j] > tmax) { tmax = ts[j]; tmax_i = j; }
                    }
                } else if (seg_aabb(cam, dir, s.nodeAABB[2*child], s.nodeAABB[2*child+1])) {
                    if (sp < kStackSize) stack[sp++] = child;
                }
            }
        }
        if (cnt == 0) continue;

        // insertion sort by t (front to back)
        for (int a = 1; a < cnt; ++a) {
            float kt = ts[a], kd = dd[a]; float3 kc = cc[a]; int b = a - 1;
            while (b >= 0 && ts[b] > kt) { ts[b+1]=ts[b]; dd[b+1]=dd[b]; cc[b+1]=cc[b]; --b; }
            ts[b+1]=kt; dd[b+1]=kd; cc[b+1]=kc;
        }
        // composite front-to-back (alpha-blend, transmittance-weighted)
        float T = 1.0f; float3 Cc = make_float3(0.f, 0.f, 0.f);
        for (int a = 0; a < cnt; ++a) {
            float al = fminf(dd[a], 0.999f);
            Cc += cc[a] * (T * al);
            T *= (1.0f - al);
            if (T < 1e-3f) break;
        }
        num += Cc; den += (1.0f - T);
    }

    if (den > 1e-6f) return num * (1.0f / den);
    float3 cs;
    if (color_static(s, p, cs)) return cs;
    return color_nearest(s, p);
}

__global__ void colorize_kernel(
    GpuScene s, const float* __restrict__ verts, int n, float* rgb, int dynamic
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float3 p = make_float3(verts[3*i], verts[3*i+1], verts[3*i+2]);
    float3 col;
    if (dynamic && false) col = color_camera(s, p);
    else if (!color_static(s, p, col)) col = color_nearest(s, p);
    rgb[3*i+0] = col.x; rgb[3*i+1] = col.y; rgb[3*i+2] = col.z;
}

// Combine the camera occlusion term with the view-independent static density,
// mirroring the BVH field occ = 1 - (1-occ_static)(1-occ_front): occ_static
// protects real surfaces (high density on any surface, regardless of view) from
// being carved by an underestimating camera, while occ_front (the render term)
// fills occluded interior and carves visible free space.
__global__ void occ_combine_kernel(
    int n, float* __restrict__ occ, const float* __restrict__ occ_static
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float a = occ_static[i], b = occ[i];
    occ[i] = 1.0f - (1.0f - a) * (1.0f - b);
}

// Fallback for the render color path: verts the cameras couldn't color (rgb<0,
// e.g. seen by no camera or fully occluded) get the static density-weighted color.
__global__ void colorize_fallback_kernel(
    GpuScene s, const float* __restrict__ verts, int n, float* rgb
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (rgb[3*i+0] >= 0.0f) return;
    float3 p = make_float3(verts[3*i], verts[3*i+1], verts[3*i+2]);
    float3 col;
    if (!color_static(s, p, col)) col = color_nearest(s, p);
    rgb[3*i+0] = col.x; rgb[3*i+1] = col.y; rgb[3*i+2] = col.z;
}

// ---------------------------------------------------------------------------
// Kernels: activation / point cloud
// ---------------------------------------------------------------------------
__global__ void activate_kernel(
    int N,
    const float* __restrict__ means,  const float* __restrict__ quats,
    const float* __restrict__ logsc,  const float* __restrict__ logit,
    const float* __restrict__ fdc,
    float3* mean, float3* ax0, float3* ax1, float3* ax2,
    float3* invs2, float* opac, float* radius, float* k2, int* valid, float3* gcol
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= N) return;

    mean[i] = make_float3(means[3*i], means[3*i+1], means[3*i+2]);

    float w = quats[4*i], x = quats[4*i+1], y = quats[4*i+2], z = quats[4*i+3];
    float n = rsqrtf(fmaxf(w*w + x*x + y*y + z*z, 1e-20f));
    w *= n; x *= n; y *= n; z *= n;
    ax0[i] = make_float3(1 - 2*(y*y + z*z), 2*(x*y + w*z),     2*(x*z - w*y));
    ax1[i] = make_float3(2*(x*y - w*z),     1 - 2*(x*x + z*z), 2*(y*z + w*x));
    ax2[i] = make_float3(2*(x*z + w*y),     2*(y*z - w*x),     1 - 2*(x*x + y*y));

    float sx = __expf(logsc[3*i]), sy = __expf(logsc[3*i+1]), sz = __expf(logsc[3*i+2]);
    invs2[i] = make_float3(1.0f/(sx*sx), 1.0f/(sy*sy), 1.0f/(sz*sz));

    // base color from the SH band-0 (DC) coefficient
    const float C0 = 0.28209479177387814f;
    gcol[i] = make_float3(
        fminf(fmaxf(0.5f + C0 * fdc[3*i+0], 0.0f), 1.0f),
        fminf(fmaxf(0.5f + C0 * fdc[3*i+1], 0.0f), 1.0f),
        fminf(fmaxf(0.5f + C0 * fdc[3*i+2], 0.0f), 1.0f));

    float op = 1.0f / (1.0f + __expf(-logit[i]));
    opac[i] = op;
    if (op > ALPHA_THRESHOLD) {
        float kk2 = 2.0f * logf(op / ALPHA_THRESHOLD);
        k2[i] = kk2;
        radius[i] = sqrtf(kk2) * fmaxf(sx, fmaxf(sy, sz));
        valid[i] = 1;
    } else {
        k2[i] = 0.0f; radius[i] = 0.0f; valid[i] = 0;
    }
}

__global__ void pointcloud_kernel(
    int num_kept, const int* __restrict__ kept,
    const float3* mean, const float3* ax0, const float3* ax1, const float3* ax2,
    const float3* invs2, const float* k2, float* out
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_kept) return;
    int g = kept[i];
    float k = sqrtf(k2[g]);
    if (true) k *= 2.5f;  // address case when splats overlap
    float3 m = mean[g];
    float3 is = invs2[g];
    float3 o0 = ax0[g] * (k * rsqrtf(is.x));
    float3 o1 = ax1[g] * (k * rsqrtf(is.y));
    float3 o2 = ax2[g] * (k * rsqrtf(is.z));
    float3 pts[7] = { m, m+o0, m-o0, m+o1, m-o1, m+o2, m-o2 };
    long base = (long)i * 21;
    for (int j = 0; j < 7; ++j) {
        out[base + 3*j + 0] = pts[j].x;
        out[base + 3*j + 1] = pts[j].y;
        out[base + 3*j + 2] = pts[j].z;
    }
}

// ---------------------------------------------------------------------------
// Kernels: LBVH build
// ---------------------------------------------------------------------------
// Per kept Gaussian: its (real-space) leaf AABB, its Morton code (remapped),
// and an iota value.
__global__ void bvh_prep_kernel(
    int num_kept, const int* __restrict__ kept,
    const float3* mean, const float3* ax0, const float3* ax1, const float3* ax2,
    const float3* invs2, const float* k2,
    float3 remap_min, float3 remap_inv_ext, float rel_scale,
    float3* leafMin, float3* leafMax, uint64_t* morton, int* iota
) {
    int kp = blockIdx.x * blockDim.x + threadIdx.x;
    if (kp >= num_kept) return;
    int g = kept[kp];
    float3 m = mean[g];
    float3 is = invs2[g];
    float s0 = 1.0f/is.x, s1 = 1.0f/is.y, s2 = 1.0f/is.z;  // sigma^2 per axis
    float3 a0 = ax0[g], a1 = ax1[g], a2 = ax2[g];
    // diagonal of Sigma = R diag(sigma^2) R^T
    float cxx = a0.x*a0.x*s0 + a1.x*a1.x*s1 + a2.x*a2.x*s2;
    float cyy = a0.y*a0.y*s0 + a1.y*a1.y*s1 + a2.y*a2.y*s2;
    float czz = a0.z*a0.z*s0 + a1.z*a1.z*s1 + a2.z*a2.z*s2;
    float k = sqrtf(k2[g]);
    float3 bound = make_float3(k*sqrtf(cxx), k*sqrtf(cyy), k*sqrtf(czz));
    leafMin[kp] = m - bound;
    leafMax[kp] = m + bound;

    float3 rc = make_float3(
        (remap_coord(m.x, rel_scale) - remap_min.x) * remap_inv_ext.x,
        (remap_coord(m.y, rel_scale) - remap_min.y) * remap_inv_ext.y,
        (remap_coord(m.z, rel_scale) - remap_min.z) * remap_inv_ext.z);
    morton[kp] = morton3D(rc);
    iota[kp] = kp;
}

// Karras 2012 single-level radix tree over the sorted Morton array.
__global__ void bvh_internal_kernel(
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

__global__ void bvh_initaabb_kernel(int n_internal, float3* nodeAABB) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n_internal) return;
    nodeAABB[2*i]   = make_float3(1e30f, 1e30f, 1e30f);
    nodeAABB[2*i+1] = make_float3(-1e30f, -1e30f, -1e30f);
}

// Bottom-up node AABBs: seed from internal nodes that have a leaf child, then
// walk to the root via parent pointers (same merge pattern as the codebase).
__global__ void bvh_aabb_kernel(
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

// ---------------------------------------------------------------------------
// Kernels: occupancy / bisection
// ---------------------------------------------------------------------------
__global__ void occ_kernel(
    GpuScene s, const float* __restrict__ pts, int n, float* occ, int dynamic
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    float3 p = make_float3(pts[3*i], pts[3*i+1], pts[3*i+2]);
    occ[i] = occ_eval(s, p, dynamic != 0);
}

__global__ void bisect_kernel(
    GpuScene s, const float* __restrict__ cloud,
    const int* __restrict__ ea, const int* __restrict__ eb,
    const float* __restrict__ oa, const float* __restrict__ ob,
    int n, int iters, int dynamic, float* out
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int a = ea[i], b = eb[i];
    float3 pa = make_float3(cloud[3*a], cloud[3*a+1], cloud[3*a+2]);
    float3 pb = make_float3(cloud[3*b], cloud[3*b+1], cloud[3*b+2]);
    // bisection narrows the bracket [pa,pb] straddling iso; we also carry the
    // occupancy at each end so the final step is a linear interpolation rather
    // than a midpoint -- much more accurate for the same iteration count.
    float occa = oa[i], occb = ob[i];
    bool a_in = occa >= s.iso;
    for (int it = 0; it < iters; ++it) {
        float3 mid = (pa + pb) * 0.5f;
        float om = occ_eval(s, mid, dynamic != 0);
        if ((om >= s.iso) == a_in) { pa = mid; occa = om; }
        else                       { pb = mid; occb = om; }
    }
    float denom = occb - occa;
    float t = (fabsf(denom) > 1e-12f) ? (s.iso - occa) / denom : 0.5f;
    t = fminf(fmaxf(t, 0.0f), 1.0f);
    float3 cross = pa + (pb - pa) * t;
    out[3*i+0] = cross.x; out[3*i+1] = cross.y; out[3*i+2] = cross.z;
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

    float3 *mean=nullptr, *ax0=nullptr, *ax1=nullptr, *ax2=nullptr, *invs2=nullptr;
    float  *opac=nullptr, *radius=nullptr, *k2=nullptr;
    float3 *gcol=nullptr;
    int    *valid=nullptr;

    int *kept = nullptr;
    int num_kept = 0;

    // LBVH
    float3 *leafMin=nullptr, *leafMax=nullptr, *nodeAABB=nullptr;
    int2   *internal=nullptr;

    float3 *campos = nullptr;
    int num_cameras = 0;

    // Camera intrinsics for the rasterize-and-sample path. Host-side copies kept
    // for the lifetime of the evaluator (the caller's tensors may be freed after
    // construction); uploaded to device lazily by the render path (Phase 2+).
    // Empty when the static (LBVH) occupancy path is in use.
    std::vector<float> cam_viewmats;  // [C*16] row-major world->cam
    std::vector<float> cam_intrins;   // [C*4]  fx, fy, cx, cy
    std::vector<float> cam_dist;      // [C*10] engine distortion layout
    std::vector<int> cam_widths;      // [C] per-camera image width
    std::vector<int> cam_heights;     // [C] per-camera image height
    std::string cam_model;
    bool use_render = false;          // cams.valid() && cameras present
    RenderContext* rctx = nullptr;    // built when use_render (MeshingRaster.cu)
    std::vector<int> render_cam_indices;  // camera subset rendered per batch

    ~Impl() {
        if (rctx) render_context_destroy(rctx);
        for (void* p : {(void*)mean,(void*)ax0,(void*)ax1,(void*)ax2,(void*)invs2,
                        (void*)opac,(void*)radius,(void*)k2,(void*)gcol,(void*)valid,(void*)kept,
                        (void*)leafMin,(void*)leafMax,(void*)nodeAABB,(void*)internal,
                        (void*)campos})
            if (p) cudaFree(p);
    }

    GpuScene make_scene() const {
        GpuScene s{};
        s.mean=mean; s.ax0=ax0; s.ax1=ax1; s.ax2=ax2; s.invs2=invs2;
        s.opac=opac; s.k2=k2; s.gcol=gcol;
        s.kept=kept; s.leafMin=leafMin; s.leafMax=leafMax;
        s.internal=internal; s.nodeAABB=nodeAABB; s.num_kept=num_kept;
        s.campos=campos; s.num_cameras=num_cameras;
        s.iso = cfg.iso;
        return s;
    }
};

OccupancyEvaluator::OccupancyEvaluator(
    const float* means, const float* quats,
    const float* log_scales, const float* logit_opac, const float* features_dc,
    int num_splats,
    const float* cam_pos, int num_cameras,
    const CameraParams& cams,
    const MeshingConfig& cfg
) {
    impl_ = new Impl();
    impl_->cfg = cfg;
    impl_->N = num_splats;
    const int N = num_splats;

    // ---- upload + activate ----
    float *d_means = dmalloc<float>((size_t)N*3), *d_quats = dmalloc<float>((size_t)N*4);
    float *d_logsc = dmalloc<float>((size_t)N*3), *d_logit = dmalloc<float>((size_t)N);
    float *d_fdc = dmalloc<float>((size_t)N*3);
    cudaMemcpy(d_means, means, sizeof(float)*N*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_quats, quats, sizeof(float)*N*4, cudaMemcpyHostToDevice);
    cudaMemcpy(d_logsc, log_scales, sizeof(float)*N*3, cudaMemcpyHostToDevice);
    cudaMemcpy(d_logit, logit_opac, sizeof(float)*N, cudaMemcpyHostToDevice);
    cudaMemcpy(d_fdc, features_dc, sizeof(float)*N*3, cudaMemcpyHostToDevice);

    impl_->mean=dmalloc<float3>(N); impl_->ax0=dmalloc<float3>(N);
    impl_->ax1=dmalloc<float3>(N);  impl_->ax2=dmalloc<float3>(N);
    impl_->invs2=dmalloc<float3>(N); impl_->opac=dmalloc<float>(N);
    impl_->radius=dmalloc<float>(N); impl_->k2=dmalloc<float>(N);
    impl_->gcol=dmalloc<float3>(N); impl_->valid=dmalloc<int>(N);

    activate_kernel<<<_LAUNCH_ARGS_1D(N, 256)>>>(
        N, d_means, d_quats, d_logsc, d_logit, d_fdc,
        impl_->mean, impl_->ax0, impl_->ax1, impl_->ax2, impl_->invs2,
        impl_->opac, impl_->radius, impl_->k2, impl_->valid, impl_->gcol);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaFree(d_means); cudaFree(d_quats); cudaFree(d_logsc); cudaFree(d_logit); cudaFree(d_fdc);

    // ---- kept list + scene bbox (host) ----
    std::vector<float> h_mean(N*3);
    std::vector<int> h_valid(N);
    cudaMemcpy(h_mean.data(), impl_->mean, sizeof(float)*N*3, cudaMemcpyDeviceToHost);
    cudaMemcpy(h_valid.data(), impl_->valid, sizeof(int)*N, cudaMemcpyDeviceToHost);

    std::vector<int> kept; kept.reserve(N);
    float bmin[3]={1e30f,1e30f,1e30f}, bmax[3]={-1e30f,-1e30f,-1e30f};
    double csum[3]={0,0,0};
    for (int i = 0; i < N; ++i) {
        if (!h_valid[i]) continue;
        kept.push_back(i);
        for (int a = 0; a < 3; ++a) {
            float m = h_mean[3*i+a];
            bmin[a] = std::min(bmin[a], m); bmax[a] = std::max(bmax[a], m);
            csum[a] += m;
        }
    }
    impl_->num_kept = (int)kept.size();
    num_kept_ = impl_->num_kept;
    num_points_ = impl_->num_kept * 7;
    if (impl_->num_kept == 0)
        throw std::runtime_error("meshing: no Gaussians above the opacity threshold");

    // rel_scale ~ RMS spread of centroids (core scene scale for the remap)
    double cmean[3] = {csum[0]/impl_->num_kept, csum[1]/impl_->num_kept, csum[2]/impl_->num_kept};
    double var = 0.0;
    for (int kp = 0; kp < impl_->num_kept; ++kp) {
        int i = kept[kp];
        for (int a = 0; a < 3; ++a) { double d = h_mean[3*i+a] - cmean[a]; var += d*d; }
    }
    float rel_scale = (float)std::max(1e-6, std::sqrt(var / (3.0 * impl_->num_kept)));

    impl_->kept = dmalloc<int>(impl_->num_kept);
    cudaMemcpy(impl_->kept, kept.data(), sizeof(int)*impl_->num_kept, cudaMemcpyHostToDevice);

    if (cfg.verbose)
        printf("[meshing] %d/%d Gaussians kept, %d points, rel_scale=%.4g\n",
               impl_->num_kept, N, num_points_, rel_scale);

    // ---- LBVH build ----
    const int n = impl_->num_kept;
    impl_->leafMin = dmalloc<float3>(n);
    impl_->leafMax = dmalloc<float3>(n);
    uint64_t* d_morton = dmalloc<uint64_t>(n);
    int* d_iota = dmalloc<int>(n);

    // remapped root bounds (remap is monotone per-axis, so remap of the real
    // bbox bounds the remapped centroids)
    float3 remap_min = make_float3(0,0,0), remap_inv_ext = make_float3(1,1,1);
    {
        float rmin[3], rext[3];
        for (int a = 0; a < 3; ++a) {
            float lo = remap_coord(bmin[a], rel_scale);
            float hi = remap_coord(bmax[a], rel_scale);
            rmin[a] = lo;
            rext[a] = std::max(hi - lo, 1e-12f);
        }
        remap_min = make_float3(rmin[0], rmin[1], rmin[2]);
        remap_inv_ext = make_float3(1.0f/rext[0], 1.0f/rext[1], 1.0f/rext[2]);
    }

    bvh_prep_kernel<<<_LAUNCH_ARGS_1D(n, 256)>>>(
        n, impl_->kept, impl_->mean, impl_->ax0, impl_->ax1, impl_->ax2,
        impl_->invs2, impl_->k2, remap_min, remap_inv_ext, rel_scale,
        impl_->leafMin, impl_->leafMax, d_morton, d_iota);
    CHECK_DEVICE_ERROR(cudaGetLastError());

    if (n >= 2) {
        // sort (morton, iota) -> argsort
        uint64_t* d_morton_s = dmalloc<uint64_t>(n);
        int* d_argsort = dmalloc<int>(n);
        CUB_WRAPPER(cub::DeviceRadixSort::SortPairs,
            d_morton, d_morton_s, d_iota, d_argsort, n);
        CHECK_DEVICE_ERROR(cudaGetLastError());

        impl_->internal = dmalloc<int2>(n - 1);
        int* d_parent = dmalloc<int>(n - 1);
        cudaMemset(d_parent, 0xff, sizeof(int)*(n - 1));
        bvh_internal_kernel<<<_LAUNCH_ARGS_1D(n - 1, 256)>>>(
            n, d_morton_s, d_argsort, impl_->internal, d_parent);
        CHECK_DEVICE_ERROR(cudaGetLastError());

        impl_->nodeAABB = dmalloc<float3>(2*(n - 1));
        bvh_initaabb_kernel<<<_LAUNCH_ARGS_1D(n - 1, 256)>>>(n - 1, impl_->nodeAABB);
        bvh_aabb_kernel<<<_LAUNCH_ARGS_1D(n - 1, 256)>>>(
            n, impl_->internal, d_parent, impl_->leafMin, impl_->leafMax, impl_->nodeAABB);
        CHECK_DEVICE_ERROR(cudaDeviceSynchronize());

        cudaFree(d_morton_s); cudaFree(d_argsort); cudaFree(d_parent);
    }
    cudaFree(d_morton); cudaFree(d_iota);

    // ---- cameras: keep all, else select max_cameras representatives ----
    if (cam_pos != nullptr && num_cameras > 0) {
        std::vector<float> cams;
        std::vector<int>& sel = impl_->render_cam_indices;
        auto max_cameras = cfg.max_cameras;
        if (max_cameras <= 0)
            max_cameras = num_cameras;
        if (num_cameras <= max_cameras) {
            cams.assign(cam_pos, cam_pos + 3 * num_cameras);
            sel.resize(num_cameras);
            for (int c = 0; c < num_cameras; ++c) sel[c] = c;
        } else {
            cams = select_cameras_kmeans(cam_pos, num_cameras, max_cameras, 25, &sel);
        }
        impl_->num_cameras = (int)(cams.size()/3);
        impl_->campos = dmalloc<float3>(impl_->num_cameras);
        cudaMemcpy(impl_->campos, cams.data(), sizeof(float)*cams.size(), cudaMemcpyHostToDevice);
        if (cfg.verbose)
            printf("[meshing] using %d/%d cameras (k-means medoids, dataset occupancy)\n",
                   impl_->num_cameras, num_cameras);
    }

    // ---- camera intrinsics (rasterize-and-sample path) ----
    // Stored host-side now; the render path uploads + uses them in Phase 2+.
    // For now this only records availability so generate_mesh can branch.
    if (cams.valid() && num_cameras > 0) {
        impl_->cam_viewmats.assign(cams.viewmats, cams.viewmats + (size_t)num_cameras * 16);
        impl_->cam_intrins.assign(cams.intrins, cams.intrins + (size_t)num_cameras * 4);
        impl_->cam_dist.assign((size_t)num_cameras * 10, 0.0f);
        if (cams.dist_coeffs)
            impl_->cam_dist.assign(cams.dist_coeffs, cams.dist_coeffs + (size_t)num_cameras * 10);
        impl_->cam_widths.assign(cams.widths, cams.widths + num_cameras);
        impl_->cam_heights.assign(cams.heights, cams.heights + num_cameras);
        impl_->cam_model  = cams.camera_model;
        impl_->use_render = true;
        impl_->rctx = render_context_create(
            means, quats, log_scales, logit_opac, features_dc, N,
            impl_->cam_viewmats.data(), impl_->cam_intrins.data(),
            impl_->cam_dist.empty() ? nullptr : impl_->cam_dist.data(),
            num_cameras, impl_->cam_widths.data(), impl_->cam_heights.data(),
            cams.camera_model, cfg.carve_k);
        if (cfg.verbose) {
            int w0 = impl_->cam_widths[0], h0 = impl_->cam_heights[0];
            bool uniform = true;
            for (int c = 1; c < num_cameras; ++c)
                if (impl_->cam_widths[c] != w0 || impl_->cam_heights[c] != h0) { uniform = false; break; }
            if (uniform)
                printf("[meshing] camera intrinsics: %dx%d, model=%s (render path ready)\n",
                       w0, h0, cams.camera_model.c_str());
            else
                printf("[meshing] camera intrinsics: per-camera resolution, model=%s (render path ready)\n",
                       cams.camera_model.c_str());
        }
    }
}

OccupancyEvaluator::~OccupancyEvaluator() { delete impl_; }

bool OccupancyEvaluator::debug_render_moments(
    int cam_idx, std::vector<float>& out, int& width, int& height
) {
    if (!impl_->rctx) return false;
    width = render_context_width(impl_->rctx, cam_idx);
    height = render_context_height(impl_->rctx, cam_idx);
    size_t npix = (size_t)width * height;
    float3* d_mom = dmalloc<float3>(npix);
    render_camera_moments(impl_->rctx, cam_idx, d_mom);
    out.resize(npix * 3);
    cudaMemcpy(out.data(), d_mom, sizeof(float3) * npix, cudaMemcpyDeviceToHost);
    cudaFree(d_mom);
    return true;
}

void OccupancyEvaluator::generate_point_cloud(std::vector<float>& xyz_out) {
    long n = (long)impl_->num_kept * 7;
    xyz_out.resize((size_t)n * 3);
    float* d_out = dmalloc<float>((size_t)n * 3);
    pointcloud_kernel<<<_LAUNCH_ARGS_1D(impl_->num_kept, 256)>>>(
        impl_->num_kept, impl_->kept,
        impl_->mean, impl_->ax0, impl_->ax1, impl_->ax2, impl_->invs2, impl_->k2, d_out);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaMemcpy(xyz_out.data(), d_out, sizeof(float)*n*3, cudaMemcpyDeviceToHost);
    cudaFree(d_out);
}

void OccupancyEvaluator::evaluate(const float* xyz, int n, float* occ_out) {
    if (n <= 0) return;
    float* d_xyz = dmalloc<float>((size_t)n*3);
    float* d_occ = dmalloc<float>(n);
    cudaMemcpy(d_xyz, xyz, sizeof(float)*n*3, cudaMemcpyHostToDevice);
    if (impl_->use_render && !impl_->render_cam_indices.empty()) {
        render_evaluate_occupancy(impl_->rctx, impl_->render_cam_indices.data(),
            (int)impl_->render_cam_indices.size(), d_xyz, n, d_occ);
        // protect surfaces with the static density term (BVH-equivalent field)
        float* d_s = dmalloc<float>(n);
        GpuScene s = impl_->make_scene();
        occ_kernel<<<_LAUNCH_ARGS_1D(n, 128)>>>(s, d_xyz, n, d_s, 0);
        occ_combine_kernel<<<_LAUNCH_ARGS_1D(n, 128)>>>(n, d_occ, d_s);
        CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
        cudaFree(d_s);
    } else {
        GpuScene s = impl_->make_scene();
        int dynamic = impl_->num_cameras > 0 ? 1 : 0;
        occ_kernel<<<_LAUNCH_ARGS_1D(n, 128)>>>(s, d_xyz, n, d_occ, dynamic);
        CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    }
    cudaMemcpy(occ_out, d_occ, sizeof(float)*n, cudaMemcpyDeviceToHost);
    cudaFree(d_xyz); cudaFree(d_occ);
}

void OccupancyEvaluator::bisect_edges(
    const float* cloud_xyz,
    const int32_t* edge_a, const int32_t* edge_b,
    const float* occ_a, const float* occ_b,
    int n_edges, float* xyz_out
) {
    if (n_edges <= 0) return;

    // Render path: host-driven bisection (each iteration renders all cameras and
    // samples the midpoint batch), so the occupancy field is consistent with
    // evaluate(). Finishes with a linear interpolation between the final bracket.
    if (impl_->use_render && !impl_->render_cam_indices.empty()) {
        const float iso = impl_->cfg.iso;
        const int* idx = impl_->render_cam_indices.data();
        const int ncam = (int)impl_->render_cam_indices.size();
        std::vector<float> lo((size_t)n_edges*3), hi((size_t)n_edges*3);
        for (int e = 0; e < n_edges; ++e) {
            int a = edge_a[e], b = edge_b[e];
            for (int t = 0; t < 3; ++t) {
                lo[3*e+t] = cloud_xyz[3*a+t];
                hi[3*e+t] = cloud_xyz[3*b+t];
            }
        }
        std::vector<float> occ_lo(occ_a, occ_a + n_edges);
        std::vector<float> occ_hi(occ_b, occ_b + n_edges);
        std::vector<float> mid((size_t)n_edges*3), occ_mid(n_edges);
        float* d_mid = dmalloc<float>((size_t)n_edges*3);
        float* d_occ = dmalloc<float>(n_edges);
        float* d_occ_s = dmalloc<float>(n_edges);
        GpuScene s = impl_->make_scene();
        for (int it = 0; it < impl_->cfg.bisection_iters; ++it) {
            for (size_t k = 0; k < (size_t)n_edges*3; ++k) mid[k] = 0.5f*(lo[k]+hi[k]);
            cudaMemcpy(d_mid, mid.data(), sizeof(float)*(size_t)n_edges*3, cudaMemcpyHostToDevice);
            render_evaluate_occupancy(impl_->rctx, idx, ncam, d_mid, n_edges, d_occ);
            occ_kernel<<<_LAUNCH_ARGS_1D(n_edges, 128)>>>(s, d_mid, n_edges, d_occ_s, 0);
            occ_combine_kernel<<<_LAUNCH_ARGS_1D(n_edges, 128)>>>(n_edges, d_occ, d_occ_s);
            CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
            cudaMemcpy(occ_mid.data(), d_occ, sizeof(float)*n_edges, cudaMemcpyDeviceToHost);
            for (int e = 0; e < n_edges; ++e) {
                bool mid_neg = (occ_mid[e] - iso) < 0.0f;
                bool lo_neg  = (occ_lo[e]  - iso) < 0.0f;
                if (mid_neg == lo_neg) {
                    for (int t = 0; t < 3; ++t) lo[3*e+t] = mid[3*e+t];
                    occ_lo[e] = occ_mid[e];
                } else {
                    for (int t = 0; t < 3; ++t) hi[3*e+t] = mid[3*e+t];
                    occ_hi[e] = occ_mid[e];
                }
            }
        }
        for (int e = 0; e < n_edges; ++e) {
            float denom = occ_hi[e] - occ_lo[e];
            float t = (fabsf(denom) > 1e-12f) ? (iso - occ_lo[e]) / denom : 0.5f;
            t = std::min(std::max(t, 0.0f), 1.0f);
            for (int a = 0; a < 3; ++a)
                xyz_out[3*e+a] = lo[3*e+a] + t * (hi[3*e+a] - lo[3*e+a]);
        }
        cudaFree(d_mid); cudaFree(d_occ); cudaFree(d_occ_s);
        return;
    }

    long ncloud = (long)num_points_;
    float* d_cloud = dmalloc<float>((size_t)ncloud*3);
    cudaMemcpy(d_cloud, cloud_xyz, sizeof(float)*ncloud*3, cudaMemcpyHostToDevice);
    int* d_ea = dmalloc<int>(n_edges); int* d_eb = dmalloc<int>(n_edges);
    float* d_oa = dmalloc<float>(n_edges); float* d_ob = dmalloc<float>(n_edges);
    float* d_out = dmalloc<float>((size_t)n_edges*3);
    cudaMemcpy(d_ea, edge_a, sizeof(int)*n_edges, cudaMemcpyHostToDevice);
    cudaMemcpy(d_eb, edge_b, sizeof(int)*n_edges, cudaMemcpyHostToDevice);
    cudaMemcpy(d_oa, occ_a, sizeof(float)*n_edges, cudaMemcpyHostToDevice);
    cudaMemcpy(d_ob, occ_b, sizeof(float)*n_edges, cudaMemcpyHostToDevice);
    GpuScene s = impl_->make_scene();
    int dynamic = impl_->num_cameras > 0 ? 1 : 0;
    bisect_kernel<<<_LAUNCH_ARGS_1D(n_edges, 128)>>>(
        s, d_cloud, d_ea, d_eb, d_oa, d_ob, n_edges, impl_->cfg.bisection_iters, dynamic, d_out);
    CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    cudaMemcpy(xyz_out, d_out, sizeof(float)*n_edges*3, cudaMemcpyDeviceToHost);
    cudaFree(d_cloud); cudaFree(d_ea); cudaFree(d_eb);
    cudaFree(d_oa); cudaFree(d_ob); cudaFree(d_out);
}

void OccupancyEvaluator::colorize(const float* verts, int n, float* rgb_out) {
    if (n <= 0) return;
    float* d_v = dmalloc<float>((size_t)n*3);
    float* d_c = dmalloc<float>((size_t)n*3);
    cudaMemcpy(d_v, verts, sizeof(float)*n*3, cudaMemcpyHostToDevice);
    GpuScene s = impl_->make_scene();
    if (impl_->use_render && !impl_->render_cam_indices.empty()) {
        render_evaluate_color(impl_->rctx, impl_->render_cam_indices.data(),
            (int)impl_->render_cam_indices.size(), d_v, n, d_c);
        colorize_fallback_kernel<<<_LAUNCH_ARGS_1D(n, 128)>>>(s, d_v, n, d_c);
        CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    } else {
        int dynamic = impl_->num_cameras > 0 ? 1 : 0;
        colorize_kernel<<<_LAUNCH_ARGS_1D(n, 128)>>>(s, d_v, n, d_c, dynamic);
        CHECK_DEVICE_ERROR(cudaDeviceSynchronize());
    }
    cudaMemcpy(rgb_out, d_c, sizeof(float)*n*3, cudaMemcpyDeviceToHost);
    cudaFree(d_v); cudaFree(d_c);
}

} // namespace meshing
