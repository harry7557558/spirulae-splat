#include <cuda_runtime.h>
#include <cstdint>

#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}


#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"



inline constexpr int kNumFrustumSegments = 16;
inline constexpr int kNumFrustumFaces = 16;

__global__ void fill_frustum_segments_kernel(
    const float4* __restrict__ intrins, // [N, 4]
    const int32_t* __restrict__ widths, // [N]
    const int32_t* __restrict__ heights, // [N]
    const bool* __restrict__ is_fisheyes, // [N]
    CameraDistortionCoeffsBuffer dist_coeffs_buffer, // [N, ...]
    const float* __restrict__ camera_to_worlds,  // [N, 3, 4]
    float size,
    float4* __restrict__ lss_buffer,  // [N, 8*kNumFrustumSegments, 2*4]
    float4* __restrict__ tri_buffer  // [N, 4*kNumFrustumFaces*kNumFrustumFaces, 4*4]
) {
    uint32_t bid = blockIdx.x;
    uint32_t tid = threadIdx.x;

    // intrins
    float4 intrin = intrins[bid];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    float width = (float)widths[bid];
    float height = (float)heights[bid];
    bool is_fisheye = is_fisheyes[bid];
    CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(bid);
    float2 corners[4] = {
        {-cx / fx, -cy / fy},
        {(width-cx) / fx, -cy / fy},
        {(width-cx) / fx, (height-cy) / fy},
        {-cx / fx, (height-cy) / fy},
    };
    float sxy = sqrtf(fx * fy / (width * height));
    float sz = rsqrtf(sxy);
    sxy *= sz * size;
    const float r0 = 0.015f*size*sz, r1 = 0.03f*size*sz;
    sz = size / sz;

    // extrins
    camera_to_worlds += 12 * bid;
    float3x3 R = {  // row major
        camera_to_worlds[0], camera_to_worlds[1], camera_to_worlds[2],  // 1st row
        camera_to_worlds[4], camera_to_worlds[5], camera_to_worlds[6],  // 2nd row
        camera_to_worlds[8], camera_to_worlds[9], camera_to_worlds[10],  // 3rd row
    };
    float3 t = { camera_to_worlds[3], camera_to_worlds[7], camera_to_worlds[11] };

    // polylines
    __shared__ float3 polylines[4*kNumFrustumSegments];
    for (int i = tid; i < 4*kNumFrustumSegments; i += blockDim.x) {
        int corner_idx = i / kNumFrustumSegments;
        float2 uv = corners[corner_idx] + (corners[(corner_idx+1)%4] - corners[corner_idx])
             * ((float)(i % kNumFrustumSegments) / kNumFrustumSegments);
        float3 raydir = float3{NAN, NAN, NAN};
        bool valid = SlangProjectionUtils::generate_ray(uv, is_fisheye, dist_coeffs, &raydir);
        if (!valid) {
            // binary search for valid
            float t0 = 0.0f, t1 = 1.0f;
            for (int iter = 0; iter < 12; ++iter) {
                float t = 0.5f*(t0+t1);
                float3 temp;
                if (SlangProjectionUtils::generate_ray(uv*t, is_fisheye, dist_coeffs, &temp))
                    t0 = t, raydir = temp;
                else
                    t1 = t;
            }
        }
        float3 p = raydir;
        if (is_fisheye)
            p = normalize(p) * sqrtf((2.0f*sxy*sxy + sz*sz) / 3.0f);
        else
            p = p * float3{sxy, sxy, sz} / p.z;
        p = float3{dot(R[0], p), dot(R[1], p), dot(R[2], p)} + t;
        polylines[i] = p;
    }
    __syncthreads();

    // frame
    lss_buffer += bid * (8*kNumFrustumSegments*2);
    for (int i = tid; i < 4*kNumFrustumSegments; i += blockDim.x) {
        float3 p0 = polylines[i];
        float3 p1 = polylines[(i+1)%(4*kNumFrustumSegments)];
        if (!isfinite(dot(p0, p1)))
            p0 = p1 = t;
        lss_buffer[2*i] = {p0.x, p0.y, p0.z, r1};
        lss_buffer[2*i+1] = {p1.x, p1.y, p1.z, r1};
    }

    // anchor
    lss_buffer += 4*kNumFrustumSegments*2;
    for (int i = tid; i < 4*kNumFrustumSegments; i += blockDim.x) {
        float3 p = polylines[(i/kNumFrustumSegments)*kNumFrustumSegments];
        float s0 = (float)(i%kNumFrustumSegments) / (float)kNumFrustumSegments;
        float s1 = s0 + 1.0f / (float)kNumFrustumSegments;
        float3 p0 = t + (p - t) * s0;
        float3 p1 = t + (p - t) * s1;
        if (!isfinite(dot(p0, p1)))
            p0 = p1 = t;
        lss_buffer[2*i] = {p0.x, p0.y, p0.z, r0 + (r1 - r0) * s0};
        lss_buffer[2*i+1] = {p1.x, p1.y, p1.z, r0 + (r1 - r0) * s1};
    }

    // faces
    tri_buffer += bid * (4*kNumFrustumFaces*kNumFrustumFaces * 4);
    for (int i = tid; i < kNumFrustumFaces*kNumFrustumFaces; i += blockDim.x) {
        float u0 = (float)(i / kNumFrustumFaces) / kNumFrustumFaces;
        float v0 = (float)(i % kNumFrustumFaces) / kNumFrustumFaces;
        float u1 = u0 + 1.0f / kNumFrustumFaces;
        float v1 = v0 + 1.0f / kNumFrustumFaces;
        float uc = 0.5f * (u0 + u1);
        float vc = 0.5f * (v0 + v1);
        float2 uvs[5] = {
            {uc, vc}, {u0, v0}, {u0, v1}, {u1, v1}, {u1, v0},
        };
        float3 verts[5];
        #pragma unroll
        for (int pi = 0; pi < 5; ++pi) {
            float u = uvs[pi].x, v = uvs[pi].y;
            float2 uv =
                corners[0] * (1.0f-u)*(1.0f-v) +
                corners[1] * u*(1.0f-v) +
                corners[2] * u*v +
                corners[3] * (1.0f-u)*v;
            float3 raydir = float3{NAN, NAN, NAN};
            bool valid = SlangProjectionUtils::generate_ray(uv, is_fisheye, dist_coeffs, &raydir);
            if (!valid) {
                // binary search for valid
                float t0 = 0.0f, t1 = 1.0f;
                for (int iter = 0; iter < 12; ++iter) {
                    float t = 0.5f*(t0+t1);
                    float3 temp;
                    if (SlangProjectionUtils::generate_ray(uv*t, is_fisheye, dist_coeffs, &temp))
                        t0 = t, raydir = temp;
                    else
                        t1 = t;
                }
            }
            float3 p = raydir;
            if (is_fisheye)
                p = normalize(p) * sqrtf((2.0f*sxy*sxy + sz*sz) / 3.0f);
            else
                p = p * float3{sxy, sxy, sz} / p.z;
            p = float3{dot(R[0], p), dot(R[1], p), dot(R[2], p)} + t;
            verts[pi] = p;
        }
        #pragma unroll
        for (int pi = 0; pi < 4; ++pi) {
            float3 p0 = verts[0];
            float3 p1 = verts[pi+1];
            float3 p2 = verts[(pi+1)%4+1];
            float2 uv0 = uvs[0];
            float2 uv1 = uvs[pi+1];
            float2 uv2 = uvs[(pi+1)%4+1];
            tri_buffer[4*(4*i+pi)+0] = {p0.x, p0.y, p0.z, p1.x};
            tri_buffer[4*(4*i+pi)+1] = {p1.y, p1.z, p2.x, p2.y};
            tri_buffer[4*(4*i+pi)+2] = {p2.z, uv0.x, uv0.y, uv1.x};
            tri_buffer[4*(4*i+pi)+3] = {uv1.y, uv2.x, uv2.y, __uint_as_float(bid)};
        }
    }
}



inline __device__ void linear_swept_sphere_aabb(
    const float4* __restrict__ lss_buffer,
    uint32_t idx,
    float3 &bmin, float3 &bmax
) {
    float4 e0 = lss_buffer[2*idx];  // xyzr
    float4 e1 = lss_buffer[2*idx+1];  // xyzr
    bmin = {
        fminf(e0.x-e0.w, e1.x-e1.w),
        fminf(e0.y-e0.w, e1.y-e1.w),
        fminf(e0.z-e0.w, e1.z-e1.w)
    };
    bmax = {
        fmaxf(e0.x+e0.w, e1.x+e1.w),
        fmaxf(e0.y+e0.w, e1.y+e1.w),
        fmaxf(e0.z+e0.w, e1.z+e1.w)
    };
}

inline __device__ void unpack_triangle(
    const float4* __restrict__ buffer,
    const uint32_t idx,
    float3 verts[3],
    float2 uvs[3],
    uint32_t &cam_idx
) {
    float4 e0 = buffer[4*idx+0];
    float4 e1 = buffer[4*idx+1];
    float4 e2 = buffer[4*idx+2];
    float4 e3 = buffer[4*idx+3];
    verts[0] = {e0.x, e0.y, e0.z};
    verts[1] = {e0.w, e1.x, e1.y};
    verts[2] = {e1.z, e1.w, e2.x};
    uvs[0] = {e2.y, e2.z};
    uvs[1] = {e2.w, e3.x};
    uvs[2] = {e3.y, e3.z};
    cam_idx = __float_as_uint(e3.w);
}

inline __device__ void triangle_aabb(
    const float4* __restrict__ tri_buffer,
    uint32_t idx,
    float3 &bmin, float3 &bmax
) {
    float3 verts[3];
    float2 uvs[3];
    uint32_t cam_idx;
    unpack_triangle(tri_buffer, idx, verts, uvs, cam_idx);
    bmin = fmin(fmin(verts[0], verts[1]), verts[2]);
    bmax = fmax(fmax(verts[0], verts[1]), verts[2]);
}

inline constexpr uint kMortonBitsPerDim = 21;

__device__ __forceinline__ uint64_t insert_2_zeros_between_bits(uint64_t x) {
    x = (x | (x << 32)) & (uint64_t)0xFFFF00000000FFFFULL;
    x = (x | (x << 16)) & (uint64_t)0x00FF0000FF0000FFULL;
    x = (x | (x << 8))  & (uint64_t)0x100F00F00F00F00FULL;
    x = (x | (x << 4))  & (uint64_t)0x10C30C30C30C30C3ULL;
    x = (x | (x << 2))  & (uint64_t)0x1249249249249249ULL;
    return x;
}

__device__ __forceinline__ uint64_t get_morton_key(float3 pos) {
    constexpr uint64_t mask_comp = ((uint64_t)1 << (uint64_t)kMortonBitsPerDim) - 1;
    uint64_t x = uint64_t(pos.x * exp2f(kMortonBitsPerDim) + 0.5f) & mask_comp;
    uint64_t y = uint64_t(pos.y * exp2f(kMortonBitsPerDim) + 0.5f) & mask_comp;
    uint64_t z = uint64_t(pos.z * exp2f(kMortonBitsPerDim) + 0.5f) & mask_comp;
    constexpr uint64_t mask_full = (((uint64_t)1 << (uint64_t)(3*kMortonBitsPerDim)) - 1);
    x = insert_2_zeros_between_bits(x) & mask_full;
    y = insert_2_zeros_between_bits(y) & mask_full;
    z = insert_2_zeros_between_bits(z) & mask_full;
    return (x * 2 + y) * 2 + z;
}

enum class VisPrimitive {
    LinearSweptSphere,
    Triangle,
};

template<VisPrimitive primitive>
__global__ void compute_aabb_kernel(
    long num_elem,
    const float4* __restrict__ buffer,
    float3* __restrict__ aabb_reduced
) {
    unsigned gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= num_elem)
        return;

    float3 aabb_min, aabb_max;
    if (primitive == VisPrimitive::LinearSweptSphere)
        linear_swept_sphere_aabb(buffer, gid, aabb_min, aabb_max);
    else if (primitive == VisPrimitive::Triangle)
        triangle_aabb(buffer, gid, aabb_min, aabb_max);

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);

    aabb_min.x = cg::reduce(warp, aabb_min.x, cg::less<float>());
    if (warp.thread_rank() == 0) atomicMin((float*)aabb_reduced+0, aabb_min.x);
    aabb_min.y = cg::reduce(warp, aabb_min.y, cg::less<float>());
    if (warp.thread_rank() == 0) atomicMin((float*)aabb_reduced+1, aabb_min.y);
    aabb_min.z = cg::reduce(warp, aabb_min.z, cg::less<float>());
    if (warp.thread_rank() == 0) atomicMin((float*)aabb_reduced+2, aabb_min.z);

    aabb_max.x = cg::reduce(warp, aabb_max.x, cg::greater<float>());
    if (warp.thread_rank() == 0) atomicMax((float*)aabb_reduced+3, aabb_max.x);
    aabb_max.y = cg::reduce(warp, aabb_max.y, cg::greater<float>());
    if (warp.thread_rank() == 0) atomicMax((float*)aabb_reduced+4, aabb_max.y);
    aabb_max.z = cg::reduce(warp, aabb_max.z, cg::greater<float>());
    if (warp.thread_rank() == 0) atomicMax((float*)aabb_reduced+5, aabb_max.z);
}

template<VisPrimitive primitive>
__global__ void fill_sorting_keys_kernel(
    const uint32_t num_elem,
    const float4* __restrict__ buffer,
    float3 root_min, float3 root_max,
    uint64_t* __restrict__ splat_keys
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_elem) return;

    float3 aabb_min, aabb_max;
    if (primitive == VisPrimitive::LinearSweptSphere)
        linear_swept_sphere_aabb(buffer, tid, aabb_min, aabb_max);
    else if (primitive == VisPrimitive::Triangle)
        triangle_aabb(buffer, tid, aabb_min, aabb_max);
    float3 aabb_center = (aabb_min + aabb_max) * 0.5f;

    uint64_t key = get_morton_key((aabb_center - root_min) / (root_max - root_min));
    splat_keys[tid] = key;
}

__global__ void fill_lbvh_internal_nodes_kernel(
    const uint32_t num_elements,
    const uint64_t* __restrict__ morton,
    const int32_t* __restrict__ element_indices,
    int2* __restrict__ internal_nodes,
    int32_t* __restrict__ parent_nodes
) {
    // https://developer.nvidia.com/blog/parallelforall/wp-content/uploads/2012/11/karras2012hpg_paper.pdf
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= num_elements-1)
        return;

    #define delta(i, j) \
        (((j)<0 || (j)>=num_elements) ? -1 : morton[i] == morton[j] ? 64 + __clz(i ^ j) : __clzll(morton[i] ^ morton[j]))

    // Determine direction of the range (+1 or -1)
    int d = delta(i,i+1) - delta(i,i-1);
    d = d > 0 ? 1 : d < 0 ? -1 : 0;

    // Compute upper bound for the length of the range
    int delta_min = delta(i, i-d);
    int lmax = 2;
    while (delta(i, i+lmax*d) > delta_min)
        lmax <<= 1;

    // Find the other end using binary search
    int l = 0;
    for (int t = lmax>>1; t >= 1; t >>= 1)
        if (delta(i, i+(l+t)*d) > delta_min)
            l += t;
    int j = i + l * d;

    // Find the split position using binary search
    int delta_node = delta(i, j);
    int s = 0;
    for (int tf = 2, t; (t = (l+tf-1)/tf) >= 1; tf <<= 1)
        if (delta(i, i+(s+t)*d) > delta_node)
            s += t;
    int gamma = i + s*d + min(d, 0);

    // Output child pointers
    // regular for internal node, bit flip for leaf
    int left = min(i,j) == gamma ? ~element_indices[gamma] : gamma;
    int right = max(i,j) == gamma+1 ? ~element_indices[gamma+1] : gamma+1;
    internal_nodes[i] = make_int2(left, right);

    // Output parent pointers
    if (left >= 0)
        atomicMax(&parent_nodes[left], i);
    if (right >= 0)
        atomicMax(&parent_nodes[right], i);

    #undef delta
}

__global__ void fill_tree_init_aabb_kernel(
    unsigned num_cells,
    float3* __restrict__ treeAABB
) {
    unsigned tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_cells) return;
    treeAABB[2*tid+0] = make_float3(1e10f);
    treeAABB[2*tid+1] = -make_float3(1e10f);
}

template<VisPrimitive primitive>
__global__ void compute_lbvh_aabb_kernel(
    const uint32_t num_elements,
    const float4* __restrict__ buffer,
    const int2* __restrict__ internal_nodes,
    const int32_t* __restrict__ parent_nodes,
    float3* __restrict__ treeAABB
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_elements-1) return;

    int2 children = internal_nodes[i];
    if (children.x >= 0 && children.y >= 0)
        return;

    // find splat AABB
    float3 aabb_min, aabb_max;
    if (children.x < 0) {
        if (primitive == VisPrimitive::LinearSweptSphere)
            linear_swept_sphere_aabb(buffer, ~children.x, aabb_min, aabb_max);
        else if (primitive == VisPrimitive::Triangle)
            triangle_aabb(buffer, ~children.x, aabb_min, aabb_max);
    }
    if (children.y < 0) {
        float3 aabb_min1, aabb_max1;
        if (primitive == VisPrimitive::LinearSweptSphere)
            linear_swept_sphere_aabb(buffer, ~children.y, aabb_min1, aabb_max1);
        else if (primitive == VisPrimitive::Triangle)
            triangle_aabb(buffer, ~children.y, aabb_min1, aabb_max1);
        if (children.x < 0)
            aabb_min = fmin(aabb_min, aabb_min1),
            aabb_max = fmax(aabb_max, aabb_max1);
        else
            aabb_min = aabb_min1, aabb_max = aabb_max1;
    }

    // fill parent AABB
    do {
        if (atomicMin(&treeAABB[2*i].x, aabb_min.x) < aabb_min.x &
            atomicMin(&treeAABB[2*i].y, aabb_min.y) < aabb_min.y &
            atomicMin(&treeAABB[2*i].z, aabb_min.z) < aabb_min.z &
            atomicMax(&treeAABB[2*i+1].x, aabb_max.x) > aabb_max.x &
            atomicMax(&treeAABB[2*i+1].y, aabb_max.y) > aabb_max.y &
            atomicMax(&treeAABB[2*i+1].z, aabb_max.z) > aabb_max.z
        ) break;
    } while ((i = parent_nodes[i]) >= 0);

}



inline __device__ float ray_aabb_intersection(
    float3 bc, float3 br,
    float3 ro, float3 rd,
    float t0, float t1
) {
    float3 n = (bc - ro) / rd;
    float3 k = br / float3{fabsf(rd.x), fabsf(rd.y), fabsf(rd.z)};
    float tmin = fmaxf(n.x-k.x, fmaxf(n.y-k.y, n.z-k.z));
    float tmax = fminf(n.x+k.x, fminf(n.y+k.y, n.z+k.z));
    bool intersect = tmin < tmax && tmax > t0 && tmin < t1;
    // if (intersect) t0 = fmaxf(t0, tmin), t1 = fminf(t1, tmax);
    return intersect ? (tmin > t0 ? tmin : tmax) : -1.0f;
}

inline __device__ float ray_sphere_intersection(
    float3 c, float r,
    float3 ro, float3 rd,
    float t0, float t1
) {
    float3 q = ro - c;
    float b = dot(q, rd);
    float h = b*b + r*r - dot(q, q);
    float d = sqrtf(h);
    float tmin = -b - d, tmax = -b + d;
    return h > 0.0f ? (tmin > t0 ? tmin : tmax) : -1.0f;
}

inline __device__ float ray_linear_swept_sphere_intersection(
    float3 pa, float3 pb,
    float ra, float rb,
    float3 ro, float3 rd,
    float t0, float t1
) {
    // https://www.shadertoy.com/view/MlKfzm by Inigo Quilez, MIT license
    float3 ba = pb - pa, oa = ro - pa, ob = ro - pb;
    float rr = ra - rb;
    float m0 = dot(ba,ba), m1 = dot(ba,oa), m2 = dot(ba,rd), m3 = dot(rd,oa),
        m5 = dot(oa,oa), m6 = dot(ob,rd), m7 = dot(ob,ob);
    float d2 = m0 - rr*rr;
	float k2 = d2 - m2*m2;
    float k1 = d2*m3 - m1*m2 + m2*rr*ra;
    float k0 = d2*m5 - m1*m1 + m1*rr*ra*2.0f - m0*ra*ra;
	float h = k1*k1 - k0*k2;
	if (h < 0.0f)
        return -1.0f;
    float t = (-sqrtf(h)-k1)/k2;
    float y = m1 - ra*rr + t*m2;
    if (y > 0.0f && y < d2)
        return (t-t0)*(t-t1) < 0.0f ? t : -1.0f;
    float h1 = m3*m3 - m5 + ra*ra;
    float h2 = m6*m6 - m7 + rb*rb;
    if (fmaxf(h1, h2) < 0.0f)
        return -1.0f;
    float r = 1e20f;
    if (h1 > 0.0f)
    	r = -m3 - sqrtf(h1);
	if (h2 > 0.0f) {
    	float t = -m6 - sqrtf(h2);
        if ((t-t0)*(t-t1) < 0.0f && t < r) r = t;
    }
    return (r-t0)*(r-t1) < 0.0f ? r : -1.0f;
}

inline __device__ bool ray_triangle_intersection(
    const float4* __restrict__ tri_buffer,
    uint32_t idx,
    float3 ro, float3 rd,
    float t0, float t1,
    float &out_t, uint32_t &out_cam_idx, float2 &out_uv
) {
    float3 verts[3];
    float2 uvs[3];
    unpack_triangle(tri_buffer, idx, verts, uvs, out_cam_idx);
    float3 po = verts[0] - ro;
    float3 a = verts[1] - verts[0];
    float3 b = verts[2] - verts[0];
    float invdet = 1.0f / dot(rd, cross(a, b));
    float t = dot(po, cross(a, b)) * invdet;
    float u = dot(rd, cross(b, po)) * invdet;
    float v = dot(rd, cross(po, a)) * invdet;
    float w = 1.0f - u - v;
    bool intersect = (t > t0 && t < t1 && u > 0.0f && v > 0.0f && w > 0.0f);
    if (intersect) {
        out_t = t;
        out_uv = uvs[0] * w + uvs[1] * u + uvs[2] * v;
    }
    return intersect;
}


__global__ void blit_aabb_bvh_kernel(
    TensorView<float, 3> render_rgbs,  // [H, W, 3]
    const TensorView<float, 3> render_depths,  // [H, W, 1]
    const TensorView<float, 3> render_alphas,  // [H, W, 1]
    const bool view_is_fisheye,
    const int image_width,
    const int image_height,
    const float4* __restrict__ view_intrins,  // [4]
    const float* __restrict__ view_viewmat,  // [4, 4]
    const CameraDistortionCoeffsBuffer view_dist_coeffs,
    const int num_aabb,
    const float3* __restrict__ aabb_buffer  // [N, 2, 3]
) {
    const int pix_x = blockIdx.x * blockDim.x + threadIdx.x;
    const int pix_y = blockIdx.y * blockDim.y + threadIdx.y;
    if (pix_x >= image_width || pix_y >= image_height)
        return;

    constexpr int MSAA = 2;
    float3 ray_o[MSAA*MSAA], ray_d[MSAA*MSAA];
    bool inside = false;
    {
        float4 intrin = view_intrins[0];
        float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
        float3x3 R = {  // row major
            view_viewmat[0], view_viewmat[1], view_viewmat[2],  // 1st row
            view_viewmat[4], view_viewmat[5], view_viewmat[6],  // 2nd row
            view_viewmat[8], view_viewmat[9], view_viewmat[10],  // 3rd row
        };
        float3 t = { view_viewmat[3], view_viewmat[7], view_viewmat[11] };
        CameraDistortionCoeffs dist_coeffs = view_dist_coeffs.load(0);

        #pragma unroll
        for (int i = 0; i < MSAA*MSAA; ++i) {
            const float px = (float)pix_x + (i/MSAA + 0.5f) / (float)MSAA;
            const float py = (float)pix_y + (i%MSAA + 0.5f) / (float)MSAA;
            ray_o[i] = SlangProjectionUtils::transform_ray_o(R, t);
            inside |= SlangProjectionUtils::generate_ray(
                {(px-cx)/fx, (py-cy)/fy},
                view_is_fisheye, dist_coeffs,
                &ray_d[i]
            );
            if (!view_is_fisheye)
                ray_d[i] = ray_d[i] * (1.0f / ray_d[i].z);
            ray_d[i] = SlangProjectionUtils::transform_ray_d(R, ray_d[i]);
        }
    }

    if (!inside)
        return;

    float depth = render_depths.load1(pix_y, pix_x);
    float opac = render_alphas.load1(pix_y, pix_x);
    depth /= fminf(fmaxf(opac, 1e-6f) / 0.5f, 1.0f);

    float counts[MSAA*MSAA] = {0.0f};

    for (int li = 0; li < num_aabb; ++li) {
        float3 p0 = aabb_buffer[2*li];
        float3 p1 = aabb_buffer[2*li+1];

        for (int i = 0; i < MSAA*MSAA; ++i) {
            // if (counts[i] == 1.0f)
            //     continue;
            float ti = ray_aabb_intersection(
                0.5f*(p0+p1), 0.5f*(p1-p0),
                ray_o[i], ray_d[i], 0.0f, depth
            );
            if (ti > 0.0f)
                counts[i] += 1.0f;
        }
    }

    float alpha_final = 0.0f;
    #pragma unroll
    for (int i = 0; i < MSAA*MSAA; ++i) {
        float a = counts[i] / (float)(MSAA*MSAA);
        alpha_final += a;
    }
    alpha_final /= (log2f(num_aabb) + 1);
    alpha_final = alpha_final / (1.0f + alpha_final);

    if (alpha_final > 0.0f) {
        float3 rgb = render_rgbs.load3(pix_y, pix_x);
        rgb = rgb * (1.0f - alpha_final) + float3{1, 0, 0} * alpha_final;
        render_rgbs.store3(pix_y, pix_x, rgb);
    }
}


__constant__ float COEFFS_turbo[3][8] = {
    { 0.14637796, 2.94711014, -10.15061040, -90.46877071, 550.44382083, -1061.31232675, 874.85901369, -266.03287948 },
    { 0.08594198, 2.06520532, 9.64581379, -55.09180383, 149.40813531, -245.26823636, 202.98596195, -63.82163439 },
    { 0.23431523, 6.33792818, 9.26380342, -203.76964931, 604.07329733, -766.82678386, 447.99287080, -97.24397473 },
};

inline __device__ float3 colormap_turbo(float t) {
    float3 c{0.0f, 0.0f, 0.0f};
    #pragma unroll
    for (int k = 7; k >= 0; --k)
        c = c * t + float3{COEFFS_turbo[0][k], COEFFS_turbo[1][k], COEFFS_turbo[2][k]};
    return c;
}

__constant__ float COEFFS_viridis[3][4] = {
    { 0.29387218, 0.04824199, -2.14339482, 2.89151588 },
    { 0.01752027, 1.23905227, -0.16754949, -0.17213704 },
    { 0.34360395, 1.23133794, -1.71784815, 0.17750210 },
};

inline __device__ float3 colormap_viridis(float t) {
    float3 c{0.0f, 0.0f, 0.0f};
    #pragma unroll
    for (int k = 3; k >= 0; --k)
        c = c * t + float3{COEFFS_viridis[0][k], COEFFS_viridis[1][k], COEFFS_viridis[2][k]};
    return c;
}

__constant__ float COEFFS_magma[3][5] = {
    { -0.00157294, 0.45897387, 4.30604404, -5.49016037, 1.69082060 },
    { -0.01042786, 0.82685306, -2.93736985, 5.48987916, -2.36381252 },
    { -0.05490989, 3.73227094, -7.52578597, 3.89534203, 0.75796302 },
};

inline __device__ float3 colormap_magma(float t) {
    float3 c{0.0f, 0.0f, 0.0f};
    #pragma unroll
    for (int k = 4; k >= 0; --k)
        c = c * t + float3{COEFFS_magma[0][k], COEFFS_magma[1][k], COEFFS_magma[2][k]};
    return c;
}

__constant__ float COEFFS_inferno[3][5] = {
    { -0.01087222, 0.74909880, 3.64171238, -5.11462265, 1.66959821 },
    { -0.00208239, 0.51010061, -1.86975118, 4.46548484, -2.07978236 },
    { 0.00011545, 2.59796038, -2.30104808, -6.74156149, 7.14515755 },
};

inline __device__ float3 colormap_inferno(float t) {
    float3 c{0.0f, 0.0f, 0.0f};
    #pragma unroll
    for (int k = 4; k >= 0; --k)
        c = c * t + float3{COEFFS_inferno[0][k], COEFFS_inferno[1][k], COEFFS_inferno[2][k]};
    return c;
}

__constant__ float COEFFS_cividis[3][4] = {
    { -0.05784352, 1.43486697, -0.95976947, 0.59905315 },
    { 0.13417334, 0.69925032, -0.07321042, 0.14530728 },
    { 0.38674306, 0.05514949, 0.66155846, -0.89703048 },
};

inline __device__ float3 colormap_cividis(float t) {
    float3 c{0.0f, 0.0f, 0.0f};
    #pragma unroll
    for (int k = 3; k >= 0; --k)
        c = c * t + float3{COEFFS_cividis[0][k], COEFFS_cividis[1][k], COEFFS_cividis[2][k]};
    return c;
}


inline __device__ float3 get_thumbnail_bilinear(
    const TensorView<uint8_t, 4> image,  // [B, H, W, 4]
    uint32_t bid,
    float2 uv
) {
    const int W = image.shape[2],
        H = image.shape[1];
    float x = fminf(fmaxf(uv.x * (float)W + 0.5f, 0.0f), (float)(W-1));
    float y = fminf(fmaxf(uv.y * (float)H + 0.5f, 0.0f), (float)(H-1));
    int x0 = (long)floorf(x);
    int x1 = x0 + 1;
    int y0 = (long)floorf(y);
    int y1 = y0 + 1;
    int wx1 = (int)(256 * (x - x0));
    int wx0 = 256 - wx1;
    int wy1 = (int)(256 * (y - y0));
    int wy0 = 256 - wy1;

    auto c00 = image.load4(bid, y0, x0);
    auto c10 = image.load4(bid, y0, x1);
    auto c01 = image.load4(bid, y1, x0);
    auto c11 = image.load4(bid, y1, x1);

    float4 rgba = float4{
        (float)(c00.x * (wx0 * wy0) + c10.x * (wx1 * wy0) + c01.x * (wx0 * wy1) + c11.x * (wx1 * wy1)),
        (float)(c00.y * (wx0 * wy0) + c10.y * (wx1 * wy0) + c01.y * (wx0 * wy1) + c11.y * (wx1 * wy1)),
        (float)(c00.z * (wx0 * wy0) + c10.z * (wx1 * wy0) + c01.z * (wx0 * wy1) + c11.z * (wx1 * wy1)),
        (float)(c00.w * (wx0 * wy0) + c10.w * (wx1 * wy0) + c01.w * (wx0 * wy1) + c11.w * (wx1 * wy1)),
    } * (1.0f / (255.0f*256.0f*256.0f));
    return float3{
        rgba.x * rgba.w,
        rgba.y * rgba.w,
        rgba.z * rgba.w + (1.0f - rgba.w)
    };
}

__global__ void blit_with_bvh_kernel(
    const TensorView<float, 3> render_rgbs,  // [H, W, 3]
    const TensorView<float, 3> render_depths,  // [H, W, 1]
    const TensorView<float, 3> render_alphas,  // [H, W, 1]
    const bool view_is_fisheye,
    const int image_width,
    const int image_height,
    const float4* __restrict__ view_intrins,  // [4]
    const float* __restrict__ view_viewmat,  // [4, 4]
    const CameraDistortionCoeffsBuffer view_dist_coeffs,
    const float4* __restrict__ lss_buffer,
    const int2* __restrict__ lss_nodes,
    const float3* __restrict__ lss_aabb,
    const float4* __restrict__ tri_buffer,
    const int2* __restrict__ tri_nodes,
    const float3* __restrict__ tri_aabb,
    const TensorView<uint8_t, 4> thumbnails,
    const float* __restrict__ min_max,
    TensorView<uint8_t, 3> out_rgbs
) {
    const int pix_x = blockIdx.x * blockDim.x + threadIdx.x;
    const int pix_y = blockIdx.y * blockDim.y + threadIdx.y;

    float3 render_rgb;
    if (render_rgbs.shape[2] == 1) {
        float amin = min_max[0], amax = min_max[1];
        float z = render_rgbs.load1(min(pix_y, image_height-1), min(pix_x, image_width-1));
        z = (z - amin) / (amax - amin);
        render_rgb = colormap_turbo(z);
    }
    else render_rgb = render_rgbs.load3(min(pix_y, image_height-1), min(pix_x, image_width-1));

    // decide whether to render camera border black or white based on background color, with noise reduction
    auto warp = cg::tiled_partition<WARP_SIZE>(cg::this_thread_block());
    static constexpr float BLOCK_DIM_X = 8;
    static constexpr float BLOCK_DIM_Y = 4;
    static constexpr float3 rgb2gray = float3{0.299f, 0.587f, 0.114f};
    // static constexpr float gray_threshold = 0.5f - 0.5f*dot(rgb2gray, normalize(rgb2gray));
    static constexpr float gray_threshold = 0.1657f;
    float gray_a = 0.0f, gray_b = 0.0f, gray_c = dot(render_rgb, rgb2gray);
    int num_gray = __popc(__ballot_sync(~0u, gray_c >= gray_threshold));
    // static constexpr int pad_gray = 2;
    static constexpr int pad_gray = 0;
    if (num_gray > pad_gray && num_gray < WARP_SIZE-pad_gray) {
        float x = (float)threadIdx.x * (1.0f / (float)BLOCK_DIM_X);
        float y = (float)threadIdx.y * (1.0f / (float)BLOCK_DIM_Y);
        #if 0
        // re-calibration needed when block size changes
        // inv([sxx sxy sx; sxy syy sy; sx sy s])
        float sxx = cg::reduce(warp, x*x, cg::plus<float>());  // TODO: closed form?
        float sxy = cg::reduce(warp, x*y, cg::plus<float>());  // TODO: closed form?
        float syy = cg::reduce(warp, y*y, cg::plus<float>());  // TODO: closed form?
        float sx = cg::reduce(warp, x, cg::plus<float>());  // TODO: closed form?
        float sy = cg::reduce(warp, y, cg::plus<float>());  // TODO: closed form?
        float s = WARP_SIZE;
        if (pix_x == 0 && pix_y == 0) printf("%f %f %f %f %f %f\n", sxx, sxy, syy, sx, sy, s);
        #endif
        float sxg = cg::reduce(warp, x*gray_c, cg::plus<float>());
        float syg = cg::reduce(warp, y*gray_c, cg::plus<float>());
        float sg = cg::reduce(warp, gray_c, cg::plus<float>());
        gray_a = 3.809523809523809e-01f * sxg - 1.666666666666667e-01f * sg;
        gray_b = 3.999999999999999e-01f * syg - 1.500000000000000e-01f * sg;
        gray_c = -1.666666666666667e-01f * sxg - 1.500000000000000e-01f * syg + 1.604166666666667e-01f * sg;
    }

    if (pix_x >= image_width || pix_y >= image_height)
        return;

    float depth = render_depths.load1(pix_y, pix_x);
    float opac = render_alphas.load1(pix_y, pix_x);
    depth /= fminf(fmaxf(opac, 1e-6f) / 0.5f, 1.0f);

    constexpr int MSAA = 2;
    float4 intrin = view_intrins[0];
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    float3x3 R = {  // row major
        view_viewmat[0], view_viewmat[1], view_viewmat[2],  // 1st row
        view_viewmat[4], view_viewmat[5], view_viewmat[6],  // 2nd row
        view_viewmat[8], view_viewmat[9], view_viewmat[10],  // 3rd row
    };
    float3 t = { view_viewmat[3], view_viewmat[7], view_viewmat[11] };
    CameraDistortionCoeffs dist_coeffs = view_dist_coeffs.load(0);

    float alpha_final = 0.0f;
    float3 rgb_final = {0.0f, 0.0f, 0.0f};

    for (int i = 0; i < MSAA*MSAA; ++i) {
        const float px = (float)pix_x + (i/MSAA + 0.5f) / (float)MSAA;
        const float py = (float)pix_y + (i%MSAA + 0.5f) / (float)MSAA;
        float gray = gray_a * (px * (1.0f / (float)BLOCK_DIM_X) - blockIdx.x) +
            gray_b * (py * (1.0f / (float)BLOCK_DIM_Y) - blockIdx.y) + gray_c;
        float3 ray_o = SlangProjectionUtils::transform_ray_o(R, t);
        float3 ray_d;
        bool inside = SlangProjectionUtils::generate_ray(
            {(px-cx)/fx, (py-cy)/fy},
            view_is_fisheye, dist_coeffs,
            &ray_d
        );
        if (!view_is_fisheye)
            // ray_d = ray_d * (1.0f / ray_d.z);
            depth *= ray_d.z;
        ray_d = SlangProjectionUtils::transform_ray_d(R, ray_d);
        if (!inside)
            continue;

        float3 rgb = {0.0f, 0.0f, 0.0f};
        float alpha = 0.0f;
        float tmax = depth;

        constexpr uint MAX_STACK_SIZE = 8*sizeof(int32_t)+1;
        __shared__ uint32_t stack[WARP_SIZE][MAX_STACK_SIZE];
        const int thread_rank = threadIdx.y * blockDim.x + threadIdx.x;

        // BVH traversal for linear swept spheres
        if (lss_buffer != nullptr) {
            uint stackSize = 0;
            float3 bmin = lss_aabb[0], bmax = lss_aabb[1];
            if (ray_aabb_intersection(
                0.5f*(bmin+bmax), 0.5f*(bmax-bmin), ray_o, ray_d, 0.0f, tmax
            ) > 0.0f) {
                stack[thread_rank][stackSize] = (uint32_t)0;
                stackSize++;
            }
            for (uint _num_steps = 0; _num_steps < 65536; _num_steps++) {
                if (stackSize == 0)
                    break;

                --stackSize;
                uint32_t elem = stack[thread_rank][stackSize];
                int2 node = lss_nodes[elem];

                for (uint ci = 0; ci < 2; ci++) {
                    int childIdx = ci == 0 ? node.x : node.y;
                    // leaf
                    if (childIdx < 0) {
                        int idx = ~childIdx;
                        float4 e0 = lss_buffer[2*idx], e1 = lss_buffer[2*idx+1];
                        float t_temp = ray_linear_swept_sphere_intersection(
                            {e0.x, e0.y, e0.z}, {e1.x, e1.y, e1.z}, e0.w, e1.w,
                            ray_o, ray_d, 0.0f, tmax
                        );
                        if (t_temp > 0.0f) {
                            tmax = t_temp;
                            rgb = make_float3(gray >= gray_threshold ? 0.1f : 0.85f);
                            alpha = 1.0f;
                        }
                    }
                    // node
                    else {
                        float3 bmin = lss_aabb[2*childIdx], bmax = lss_aabb[2*childIdx+1];
                        if (ray_aabb_intersection(
                            0.5f*(bmin+bmax), 0.5f*(bmax-bmin), ray_o, ray_d, 0.0f, tmax
                        ) > 0.0f && stackSize < MAX_STACK_SIZE) {
                            stack[thread_rank][stackSize] = (uint32_t)childIdx;
                            stackSize += 1;
                        }
                    }
                }
            }
        }

        // BVH traversal for triangle
        if (tri_buffer != nullptr) {
            uint stackSize = 0;
            float3 bmin = tri_aabb[0], bmax = tri_aabb[1];
            if (ray_aabb_intersection(
                0.5f*(bmin+bmax), 0.5f*(bmax-bmin), ray_o, ray_d, 0.0f, tmax
            ) > 0.0f) {
                stack[thread_rank][stackSize] = (uint32_t)0;
                stackSize++;
            }
            for (uint _num_steps = 0; _num_steps < 65536; _num_steps++) {
                if (stackSize == 0)
                    break;

                --stackSize;
                uint32_t elem = stack[thread_rank][stackSize];
                int2 node = tri_nodes[elem];

                for (uint ci = 0; ci < 2; ci++) {
                    int childIdx = ci == 0 ? node.x : node.y;
                    // leaf
                    if (childIdx < 0) {
                        float t_temp;
                        uint32_t cam_idx;
                        float2 uv;
                        if (ray_triangle_intersection(
                            tri_buffer, ~childIdx,
                            ray_o, ray_d, 0.0f, tmax,
                            t_temp, cam_idx, uv
                        )) {
                            tmax = t_temp;
                            rgb = get_thumbnail_bilinear(thumbnails, cam_idx, uv);
                            alpha = 1.0f;
                        }
                    }
                    // node
                    else {
                        float3 bmin = tri_aabb[2*childIdx], bmax = tri_aabb[2*childIdx+1];
                        if (ray_aabb_intersection(
                            0.5f*(bmin+bmax), 0.5f*(bmax-bmin), ray_o, ray_d, 0.0f, tmax
                        ) > 0.0f && stackSize < MAX_STACK_SIZE) {
                            stack[thread_rank][stackSize] = (uint32_t)childIdx;
                            stackSize += 1;
                        }
                    }
                }
            }
        }

        float a = alpha / (float)(MSAA*MSAA);
        alpha_final += a;
        rgb_final += rgb * a;
    }

    if (alpha_final > 0.0f)
        render_rgb = render_rgb * (1.0f - alpha_final) + rgb_final;
    render_rgb = 255.0f * fmin(fmax(render_rgb, make_float3(0.0f)), make_float3(1.0f)) + 0.5f;
    out_rgbs.store3(pix_y, pix_x, {(uint8_t)render_rgb.x, (uint8_t)render_rgb.y, (uint8_t)render_rgb.z});
}

template<VisPrimitive primitive>
inline std::tuple<at::Tensor, at::Tensor> build_bvh(
    uint32_t num_elem,
    const at::Tensor& buffer,
    float3 rootAABBMin,
    float3 rootAABBMax,
    c10::TensorOptions options
) {
    // Morton sorting
    at::Tensor morton = at::empty({num_elem}, options.dtype(at::kLong));
    fill_sorting_keys_kernel<primitive>
    <<<_LAUNCH_ARGS_1D(num_elem, 256)>>>(
        num_elem,
        (float4*)buffer.data_ptr<float>(),
        rootAABBMin, rootAABBMax,
        (uint64_t*)morton.data_ptr<int64_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    auto [sorted_morton, argsort] = at::sort(morton);
    argsort = argsort.to(at::kInt);

    // Build tree
    at::Tensor internal_nodes = at::empty({num_elem-1, 2}, options.dtype(at::kInt));
    at::Tensor parent_nodes = at::empty({num_elem-1}, options.dtype(at::kInt));
    cudaMemset(parent_nodes.data_ptr<int32_t>(), 0xff, (num_elem-1)*sizeof(int32_t));
    CHECK_DEVICE_ERROR(cudaGetLastError());
    fill_lbvh_internal_nodes_kernel<<<_LAUNCH_ARGS_1D(num_elem-1, 256)>>>(
        num_elem,
        (uint64_t*)sorted_morton.data_ptr<int64_t>(),
        argsort.data_ptr<int32_t>(),
        (int2*)internal_nodes.data_ptr<int32_t>(),
        parent_nodes.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // Compute AABB
    at::Tensor treeAABB = at::empty({num_elem, 2, 3}, options);
    fill_tree_init_aabb_kernel<<<_LAUNCH_ARGS_1D(num_elem-1, 256)>>>(
        num_elem-1,
        (float3*)treeAABB.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    compute_lbvh_aabb_kernel<primitive>
    <<<_LAUNCH_ARGS_1D(num_elem-1, 256)>>>(
        num_elem,
        (float4*)buffer.data_ptr<float>(),
        (int2*)internal_nodes.data_ptr<int32_t>(),
        parent_nodes.data_ptr<int32_t>(),
        (float3*)treeAABB.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(internal_nodes, treeAABB);
}

__global__ void compute_min_max_kernel(
    TensorView<float, 3> depths,
    float* __restrict__ min_max  // [2]
) {
    uint32_t xid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t yid = blockIdx.y * blockDim.y + threadIdx.y;

    float amin = 1e38f, amax = -1e38f;
    if (xid < depths.shape[1] && yid < depths.shape[0]) {
        float z = depths.load1(yid, xid);
        if (isfinite(z))
            amin = fminf(amin, z), amax = fmaxf(amax, z);
    }

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);

    amin = cg::reduce(warp, amin, cg::less<float>());
    if (warp.thread_rank() == 0) atomicMin(&min_max[0], amin);

    amax = cg::reduce(warp, amax, cg::greater<float>());
    if (warp.thread_rank() == 0) atomicMax(&min_max[1], amax);
}


/*[AutoHeaderGeneratorExport]*/
at::Tensor blit_train_cameras_tensor(
    at::Tensor render_rgbs,  // [H, W, 3]
    at::Tensor render_depths,  // [H, W, 1]
    at::Tensor render_alphas,  // [H, W, 1]
    const bool view_is_fisheye,
    const at::Tensor view_intrins,  // [4]
    const at::Tensor view_viewmat,  // [4, 4]
    const CameraDistortionCoeffsTensor view_dist_coeffs,
    const at::Tensor intrins,  // [N, 4]
    const at::Tensor widths,  // [N]
    const at::Tensor heights,  // [N]
    const at::Tensor is_fisheye,  // [N]
    const CameraDistortionCoeffsTensor dist_coeffs,
    const at::Tensor camera_to_worlds,  // [N, 3, 4]
    at::Tensor thumbnails,
    const float camera_size,
    const bool show_training_cameras
) {
    DEVICE_GUARD(render_rgbs);
    CHECK_CUDA(render_rgbs);
    CHECK_CUDA(render_depths);
    CHECK_CUDA(render_alphas);
    CHECK_INPUT(view_intrins);
    CHECK_INPUT(view_viewmat);
    CHECK_INPUT(intrins);
    CHECK_INPUT(widths);
    CHECK_INPUT(heights);
    CHECK_INPUT(is_fisheye);
    CHECK_INPUT(camera_to_worlds);
    CHECK_INPUT(thumbnails);

    auto h = render_rgbs.size(-3);
    auto w = render_rgbs.size(-2);
    auto c = render_rgbs.size(-1);
    auto n = intrins.size(0);

    constexpr int kFloatPInfByte = 0x7f;  // 0x7f7f7f7f -> 3.39615e+38
    constexpr int kFloatNInfByte = 0xfe;  // 0xfefefefe -> -1.69474e+38

    at::Tensor min_max_tensor;
    if (c == 1) {
        min_max_tensor = at::empty({2}, render_rgbs.options());
        cudaMemset(min_max_tensor.data_ptr<float>()+0, kFloatPInfByte, sizeof(float));
        cudaMemset(min_max_tensor.data_ptr<float>()+1, kFloatNInfByte, sizeof(float));
        compute_min_max_kernel<<<_LAUNCH_ARGS_2D(w, h, 16, 16)>>>(
            tensor2view<float, 3>(render_rgbs),
            min_max_tensor.data_ptr<float>()
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    if (!show_training_cameras) {
        at::Tensor out_rgb = at::empty({h, w, 3}, render_rgbs.options().dtype(at::kByte));
        blit_with_bvh_kernel<<<_LAUNCH_ARGS_2D(w, h, 8, 4)>>>(
            tensor2view<float, 3>(render_rgbs),
            tensor2view<float, 3>(render_depths),
            tensor2view<float, 3>(render_alphas),
            view_is_fisheye,
            w, h,
            (float4*)view_intrins.data_ptr<float>(),
            view_viewmat.data_ptr<float>(),
            view_dist_coeffs,
            nullptr, nullptr, nullptr,
            nullptr, nullptr, nullptr,
            tensor2view<uint8_t, 4>(thumbnails),
            c == 1 ? min_max_tensor.data_ptr<float>() : nullptr,
            tensor2view<uint8_t, 3>(out_rgb)
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
        return out_rgb;
    }

    uint32_t num_lss = n*8*kNumFrustumSegments;
    uint32_t num_tri = n*4*kNumFrustumFaces*kNumFrustumFaces;
    at::Tensor lss_buffer = at::empty({num_lss, 8}, render_rgbs.options());
    at::Tensor tri_buffer = at::empty({num_tri, 16}, render_rgbs.options());
    fill_frustum_segments_kernel
    <<<_LAUNCH_ARGS_1D(n*4*kNumFrustumSegments, 4*kNumFrustumSegments)>>>(
        (float4*)intrins.data_ptr<float>(),
        widths.data_ptr<int32_t>(),
        heights.data_ptr<int32_t>(),
        is_fisheye.data_ptr<bool>(),
        dist_coeffs,
        camera_to_worlds.data_ptr<float>(),
        camera_size,
        (float4*)lss_buffer.data_ptr<float>(),
        (float4*)tri_buffer.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // compute global AABB
    at::Tensor root_aabb_tensor = at::empty({2, 3}, render_rgbs.options());
    cudaMemset(root_aabb_tensor.data_ptr<float>()+0, kFloatPInfByte, 3*sizeof(float));
    cudaMemset(root_aabb_tensor.data_ptr<float>()+3, kFloatNInfByte, 3*sizeof(float));
    compute_aabb_kernel<VisPrimitive::LinearSweptSphere>
    <<<_LAUNCH_ARGS_1D(num_lss, 256)>>>(
        num_lss,
        (float4*)lss_buffer.data_ptr<float>(),
        (float3*)root_aabb_tensor.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    compute_aabb_kernel<VisPrimitive::Triangle>
    <<<_LAUNCH_ARGS_1D(num_tri, 256)>>>(
        num_tri,
        (float4*)tri_buffer.data_ptr<float>(),
        (float3*)root_aabb_tensor.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    glm::vec3 rootAABBMin, rootAABBMax;
    {
        root_aabb_tensor = root_aabb_tensor.cpu();
        glm::vec3* root_aabb = (glm::vec3*)root_aabb_tensor.data_ptr<float>();
        rootAABBMin = root_aabb[0];
        rootAABBMax = root_aabb[1];
        glm::vec3 center = 0.5f * (rootAABBMax + rootAABBMin);
        glm::vec3 extend = 0.5f * (rootAABBMax - rootAABBMin);
        float max_size = 1.01f * fmax(extend.x, fmax(extend.y, extend.z));
        rootAABBMin = center - glm::vec3(max_size);
        rootAABBMax = center + glm::vec3(max_size);
    }
    // printf("AABB: %f %f %f  %f %f %f\n", rootAABBMin.x, rootAABBMin.y, rootAABBMin.z, rootAABBMax.x, rootAABBMax.y, rootAABBMax.z);

    // build tree for LSS
    auto [lss_nodes, lss_aabb] = build_bvh<VisPrimitive::LinearSweptSphere>(
        num_lss, lss_buffer,
        *(float3*)&rootAABBMin, *(float3*)&rootAABBMax,
        render_rgbs.options()
    );
    auto [tri_nodes, tri_aabb] = build_bvh<VisPrimitive::Triangle>(
        num_tri, tri_buffer,
        *(float3*)&rootAABBMin, *(float3*)&rootAABBMax,
        render_rgbs.options()
    );

    #if 0
    blit_aabb_bvh_kernel<<<_LAUNCH_ARGS_2D(w, h, 8, 4)>>>(
        tensor2view<float, 3>(render_rgbs),
        tensor2view<float, 3>(render_depths),
        tensor2view<float, 3>(render_alphas),
        view_is_fisheye,
        w, h,
        (float4*)view_intrins.data_ptr<float>(),
        view_viewmat.data_ptr<float>(),
        view_dist_coeffs,
        num_lss-2,
        (float3*)lss_aabb.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    #endif

    at::Tensor out_rgb = at::empty({h, w, 3}, render_rgbs.options().dtype(at::kByte));

    blit_with_bvh_kernel<<<_LAUNCH_ARGS_2D(w, h, 8, 4)>>>(
        tensor2view<float, 3>(render_rgbs),
        tensor2view<float, 3>(render_depths),
        tensor2view<float, 3>(render_alphas),
        view_is_fisheye,
        w, h,
        (float4*)view_intrins.data_ptr<float>(),
        view_viewmat.data_ptr<float>(),
        view_dist_coeffs,
        (float4*)lss_buffer.data_ptr<float>(),
        (int2*)lss_nodes.data_ptr<int32_t>(),
        (float3*)lss_aabb.data_ptr<float>(),
        (float4*)tri_buffer.data_ptr<float>(),
        (int2*)tri_nodes.data_ptr<int32_t>(),
        (float3*)tri_aabb.data_ptr<float>(),
        tensor2view<uint8_t, 4>(thumbnails),
        c == 1 ? min_max_tensor.data_ptr<float>() : nullptr,
        tensor2view<uint8_t, 3>(out_rgb)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_rgb;
}
