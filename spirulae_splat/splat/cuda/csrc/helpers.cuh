#pragma once

#include "common.cuh"
#include <cuda_runtime.h>
#include <cuda_fp16.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace cg = cooperative_groups;

template <typename T> using vec2 = glm::vec<2, T>;
template <typename T> using vec3 = glm::vec<3, T>;
template <typename T> using vec4 = glm::vec<4, T>;
template <typename T> using mat2 = glm::mat<2, 2, T>;
template <typename T> using mat3 = glm::mat<3, 3, T>;
template <typename T> using mat4 = glm::mat<4, 4, T>;
template <typename T> using mat3x2 = glm::mat<3, 2, T>;
template <typename T> using mat2x3 = glm::mat<2, 3, T>;

struct halfc {
    __half value;
    __host__ __device__ halfc() {}
    __host__ __device__ halfc(float f) : value(__float2half(f)) {}
    __host__ __device__ halfc(__half h) : value(h) {}

    __host__ __device__ operator float() const { return __half2float(value); }
    __host__ __device__ operator __half() const { return value; }
};

inline __device__ halfc atomicAdd(halfc* address, halfc val) {
    __half* address_v = (half*)address;
    __half val_v = val.value;
    __half v = atomicAdd(address_v, val_v);
    return halfc(v);
}


inline __device__ float safe_denom(float x, float e) {
    return abs(x)>e ? x : x>0.f ? e : -e;
}

inline __device__ float nan_to_num(float x, float v=0.0f) {
    return isfinite(x) ? x : v;
}
inline __device__ float2 nan_to_num(float2 x, float v=0.0f) {
    return { nan_to_num(x.x), nan_to_num(x.y) };
}
inline __device__ float3 nan_to_num(float3 x, float v=0.0f) {
    return { nan_to_num(x.x), nan_to_num(x.y), nan_to_num(x.z) };
}
inline __device__ float4 nan_to_num(float4 x, float v=0.0f) {
    return { nan_to_num(x.x), nan_to_num(x.y), nan_to_num(x.z), nan_to_num(x.w) };
}


#if SPLAT_KERNEL == 1

// "Gaussian" kernel
inline __device__ float visibility_kernel(const float r2) {
    // return r2 > 1.0f ? 0.0f : (1.0f-r2) * (1.0f-r2);
    return r2 > 1.0f ? 0.0f : 1.0f-r2;
}

// gradient of "Gaussian" kernel
inline __device__ float visibility_kernel_grad(const float r2) {
    // return r2 > 1.0f ? 0.0f : 2.0f * (r2-1.0f);
    return r2 > 1.0f ? 0.0f : -1.0f;
}

// radius of "Gaussian" kernel
inline __device__ float visibility_kernel_radius() {
    return 1.0f;
}

#elif SPLAT_KERNEL == 0

// Gaussian kernel
inline __device__ float visibility_kernel(const float r2) {
    return expf(-r2);
}

// gradient of Gaussian kernel
inline __device__ float visibility_kernel_grad(const float r2) {
    return -expf(-r2);
}

// radius of Gaussian kernel
inline __device__ float visibility_kernel_radius() {
    return 3.0f;
}

#endif


inline __device__ float depth_map(float z) {
    return z>0.0f ? logf(z+1.0f) : z;
}

inline __device__ float depth_map_vjp(float z, float v_z) {
    return z>0.0f ? v_z/(z+1.0f) : v_z;
}

inline __device__ float depth_inv_map(float z) {
    return z>0.0f ? expf(z)-1.0f : z;
}

inline __device__ float depth_inv_map_vjp(float z, float v_z) {
    return z>0.0f ? v_z*expf(z) : v_z;
}


inline __device__ bool get_alpha(
    const glm::vec2 uv,
    const float opac,
    // const glm::vec2 aniso,
    float &alpha
) {
    const float r2 = glm::dot(uv, uv);
    const float vis = visibility_kernel(r2);
    const float t = 0.0f; //glm::dot(uv, aniso);
    const float m = t<0.f ? 1.f : t>1.f ? 0.f :
        t*t*(2.0f*t-3.0f) + 1.0f;
    alpha = opac * vis * m;
    return r2 >= 0.f && alpha > 1e-3f;
}

inline __device__ void get_alpha_vjp(
    const glm::vec2 uv,
    const float opac,
    // const glm::vec2 aniso,
    const float v_alpha,
    glm::vec2 &v_uv,
    float &v_opac
    // glm::vec2 &v_aniso
) {
    const float r2 = glm::dot(uv, uv);
    const float vis = visibility_kernel(r2);
    const float t = 0.0f; //glm::dot(uv, aniso);
    const float m = t<0.f ? 1.f : t>1.f ? 0.f :
        t*t*(2.0f*t-3.0f) + 1.0f;
    // const float alpha = opac * vis * m;
    // const float v_m = opac * vis * v_alpha;
    // const float v_t = t<0.f||t>1.f ? 0.f : 6.0f*t*(t-1.0f) * v_m;
    const float v_vis = opac * m * v_alpha;
    const float v_r2 = v_vis * visibility_kernel_grad(r2);
    v_uv = 2.0f * uv * v_r2;// + v_t * aniso;
    v_opac = vis * m * v_alpha;
    // v_aniso = uv * v_t;
}

inline __device__ bool get_intersection(
    const glm::vec3 position,
    const glm::mat2x3 axis_uv,
    const glm::vec2 pos_2d,
    glm::vec3 &poi, glm::vec2 &uv
) {
    const float radius = visibility_kernel_radius();
    glm::mat3 A = glm::mat3(
        axis_uv[0], axis_uv[1],
        glm::vec3(pos_2d, 1.0f)
    );
    #if 0
    if (glm::determinant(A) == 0.0f) {
        uv = {-radius, -radius};
        return false;
    }
    #endif
    glm::vec3 uvt = -glm::inverse(A) * position;
    uv = {uvt.x, uvt.y};
    float r2 = glm::dot(uv, uv);
    float t = -uvt.z;
    if (r2 > radius*radius || isnan(r2) || !isfinite(t)) {
        uv = {-radius, -radius};
        return false;
    }
    poi = glm::vec3(pos_2d*t, t);
    return true;
}

inline __device__ void get_intersection_vjp(
    const glm::vec3 position,
    const glm::mat2x3 axis_uv,
    const glm::vec2 pos_2d,
    const glm::vec3 v_poi,
    const glm::vec2 v_uv,
    glm::vec3 &v_position,
    glm::mat2x3 &v_axis_uv
) {
    const float radius = visibility_kernel_radius();
    // forward
    glm::mat3 A = glm::mat3(
        axis_uv[0], axis_uv[1],
        glm::vec3(pos_2d, 1.0f)
    );
    glm::mat3 invA = glm::inverse(A);
    glm::vec3 uvt = -invA * position;
    glm::vec2 uv = {uvt.x, uvt.y};
    float t = -uvt.z;
    glm::vec3 poi = glm::vec3(pos_2d*t, t);
    // backward
    float v_t = glm::dot(v_poi, glm::vec3(pos_2d, 1.0f));
    glm::vec3 v_uvt = {v_uv.x, v_uv.y, -v_t};
    v_position = -glm::transpose(invA) * v_uvt;
    glm::mat3 v_invA = -glm::outerProduct(v_uvt, position);
    glm::mat3 v_A = -glm::transpose(invA) * v_invA * glm::transpose(invA);
    v_axis_uv += glm::mat2x3(v_A[0], v_A[1]);
}


inline __device__ glm::vec3 get_normal_from_axisuv(
    const glm::mat2x3 &axis_uv,
    const glm::vec3 &vdir
) {
    glm::vec3 normal = glm::cross(axis_uv[0], axis_uv[1]) / 
        sqrtf(fmaxf(glm::dot(axis_uv[0], axis_uv[0]) * glm::dot(axis_uv[1], axis_uv[1]), 1e-12f));
    float normal_mult = glm::dot(vdir, normal) > 0.0f ? -1.0f : 1.0f;
    return normal * normal_mult;
}

inline __device__ void get_normal_from_axisuv_vjp(
    const glm::mat2x3 &axis_uv,
    const glm::vec3 &vdir,
    const glm::vec3 &v_normal,
    glm::vec3 &normal,
    glm::mat2x3 &v_axis_uv
) {
    float norm_u = fmaxf(glm::length(axis_uv[0]), 1e-6f);
    float norm_v = fmaxf(glm::length(axis_uv[1]), 1e-6f);
    glm::vec3 u_norm = axis_uv[0] / norm_u;
    glm::vec3 v_norm = axis_uv[1] / norm_v;
    glm::vec3 normal0 = glm::cross(u_norm, v_norm);
    float normal_mult = glm::dot(vdir, normal0) > 0.0f ? -1.0f : 1.0f;
    normal = normal0 * normal_mult;
    glm::vec3 v_normal0 = v_normal * normal_mult;
    glm::vec3 v_u_norm = glm::cross(v_norm, v_normal0);
    glm::vec3 v_v_norm = glm::cross(v_normal0, u_norm);
    glm::vec3 v_u_nr = (v_u_norm - u_norm * glm::dot(u_norm, v_u_norm)) / norm_u;
    glm::vec3 v_v_nr = (v_v_norm - v_norm * glm::dot(v_norm, v_v_norm)) / norm_v;
    v_axis_uv = glm::mat2x3(v_u_nr, v_v_nr);
}



inline __device__ void get_bbox(
    const float2 center,
    const float2 dims,
    const dim3 img_size,
    int2 &bb_min,
    int2 &bb_max
) {
    // get bounding box with center and dims, within bounds
    // bounding box coords returned in tile coords, inclusive min, exclusive max
    // clamp between 0 and tile bounds
    bb_min.x = min(max(0, (int)(center.x - dims.x)), img_size.x);
    bb_max.x = min(max(0, (int)(center.x + dims.x + 1)), img_size.x);
    bb_min.y = min(max(0, (int)(center.y - dims.y)), img_size.y);
    bb_max.y = min(max(0, (int)(center.y + dims.y + 1)), img_size.y);
}

inline __device__ void get_tile_bbox(
    const float2 pix_center,
    const float2 pix_radius,
    const dim3 tile_bounds,
    int2 &tile_min, int2 &tile_max
) {
    // gets gaussian dimensions in tile space, i.e. the span of a gaussian in
    // tile_grid (image divided into tiles)
    float2 tile_center = {
        pix_center.x / (float)TILE_SIZE, pix_center.y / (float)TILE_SIZE
    };
    float2 tile_radius = {
        pix_radius.x / (float)TILE_SIZE, pix_radius.y / (float)TILE_SIZE
    };
    get_bbox(tile_center, tile_radius, tile_bounds, tile_min, tile_max);
}

// compute vjp from df/d_conic to df/c_cov2d
inline __device__ void cov2d_to_conic_vjp(
    const float3 &conic, const float3 &v_conic, float3 &v_cov2d
) {
    // conic = inverse cov2d
    // df/d_cov2d = -conic * df/d_conic * conic
    glm::mat2 X = glm::mat2(conic.x, conic.y, conic.y, conic.z);
    glm::mat2 G = glm::mat2(v_conic.x, v_conic.y / 2.f, v_conic.y / 2.f, v_conic.z);
    glm::mat2 v_Sigma = -X * G * X;
    v_cov2d.x = v_Sigma[0][0];
    v_cov2d.y = v_Sigma[1][0] + v_Sigma[0][1];
    v_cov2d.z = v_Sigma[1][1];
}

inline __device__ void cov2d_to_compensation_vjp(
    const float compensation, const float3 &conic, const float v_compensation, float3 &v_cov2d
) {
    // comp = sqrt(det(cov2d - 0.3 I) / det(cov2d))
    // conic = inverse(cov2d)
    // df / d_cov2d = df / d comp * 0.5 / comp * [ d comp^2 / d cov2d ]
    // d comp^2 / d cov2d = (1 - comp^2) * conic - 0.3 I * det(conic)
    float inv_det = conic.x * conic.z - conic.y * conic.y;
    float one_minus_sqr_comp = 1 - compensation * compensation;
    float v_sqr_comp = v_compensation * 0.5 / (compensation + 1e-6);
    v_cov2d.x += v_sqr_comp * (one_minus_sqr_comp * conic.x - 0.3 * inv_det);
    v_cov2d.y += 2 * v_sqr_comp * (one_minus_sqr_comp * conic.y);
    v_cov2d.z += v_sqr_comp * (one_minus_sqr_comp * conic.z - 0.3 * inv_det);
}

// helper for applying R^T * p for a ROW MAJOR 4x3 matrix [R, t], ignoring t
inline __device__ float3 transform_4x3_rot_only_transposed(const float *mat, const float3 p) {
    float3 out = {
        mat[0] * p.x + mat[4] * p.y + mat[8] * p.z,
        mat[1] * p.x + mat[5] * p.y + mat[9] * p.z,
        mat[2] * p.x + mat[6] * p.y + mat[10] * p.z,
    };
    return out;
}

inline __device__ float3 transform_4x3_rot_only(const float *mat, const float3 p) {
    float3 out = {
        mat[0] * p.x + mat[1] * p.y + mat[2] * p.z,
        mat[4] * p.x + mat[5] * p.y + mat[6] * p.z,
        mat[8] * p.x + mat[9] * p.y + mat[10] * p.z,
    };
    return out;
}

// helper for applying R * p + T, expect mat to be ROW MAJOR
inline __device__ float3 transform_4x3(const float *mat, const float3 p) {
    float3 out = {
        mat[0] * p.x + mat[1] * p.y + mat[2] * p.z + mat[3],
        mat[4] * p.x + mat[5] * p.y + mat[6] * p.z + mat[7],
        mat[8] * p.x + mat[9] * p.y + mat[10] * p.z + mat[11],
    };
    return out;
}

// helper to apply 4x4 transform to 3d vector, return homo coords
// expects mat to be ROW MAJOR
inline __device__ float4 transform_4x4(const float *mat, const float3 p) {
    float4 out = {
        mat[0] * p.x + mat[1] * p.y + mat[2] * p.z + mat[3],
        mat[4] * p.x + mat[5] * p.y + mat[6] * p.z + mat[7],
        mat[8] * p.x + mat[9] * p.y + mat[10] * p.z + mat[11],
        mat[12] * p.x + mat[13] * p.y + mat[14] * p.z + mat[15],
    };
    return out;
}

inline __device__ float2 project_pix(
    const float2 f, const float3 p_view, const float2 c
) {
    float rw = 1.f / (p_view.z + 1e-6f);
    float2 p_proj = { p_view.x * rw, p_view.y * rw };
    float2 p_pix = { p_proj.x * f.x + c.x, p_proj.y * f.y + c.y };
    return p_pix;
}

// given v_xy_pix, get v_xyz
inline __device__ float3 project_pix_vjp(
    const float2 fxfy, const float3 p_view, const float2 v_xy
) {
    float rw = 1.f / (p_view.z + 1e-6f);
    float2 v_proj = { fxfy.x * v_xy.x, fxfy.y * v_xy.y };
    float3 v_view = {
        v_proj.x * rw, v_proj.y * rw, -(v_proj.x * p_view.x + v_proj.y * p_view.y) * rw * rw
    };
    return v_view;
}

inline __device__ glm::mat3 quat_to_rotmat(const float4 quat) {
    // quat to rotation matrix
    float w = quat.x;
    float x = quat.y;
    float y = quat.z;
    float z = quat.w;

    // glm matrices are column-major
    return glm::mat3(
        1.f - 2.f * (y * y + z * z),  // [0][0]
        2.f * (x * y + w * z),  // [0][1]
        2.f * (x * z - w * y),  // [0][2]
        2.f * (x * y - w * z),  // [1][0]
        1.f - 2.f * (x * x + z * z),  // [1][1]
        2.f * (y * z + w * x),  // [1][2]
        2.f * (x * z + w * y),  // [2][0]
        2.f * (y * z - w * x),  // [2][1]
        1.f - 2.f * (x * x + y * y)  // [2][2]
    );
}

inline __device__ float4
quat_to_rotmat_vjp(const float4 quat, const glm::mat3 v_R) {
    float w = quat.x;
    float x = quat.y;
    float y = quat.z;
    float z = quat.w;

    float4 v_quat;
    // v_R is COLUMN MAJOR
    // w element stored in x field
    v_quat.x =
        2.f * (
                  // v_quat.w = 2.f * (
                  x * (v_R[1][2] - v_R[2][1]) + y * (v_R[2][0] - v_R[0][2]) +
                  z * (v_R[0][1] - v_R[1][0])
              );
    // x element in y field
    v_quat.y =
        2.f *
        (
            // v_quat.x = 2.f * (
            -2.f * x * (v_R[1][1] + v_R[2][2]) + y * (v_R[0][1] + v_R[1][0]) +
            z * (v_R[0][2] + v_R[2][0]) + w * (v_R[1][2] - v_R[2][1])
        );
    // y element in z field
    v_quat.z =
        2.f *
        (
            // v_quat.y = 2.f * (
            x * (v_R[0][1] + v_R[1][0]) - 2.f * y * (v_R[0][0] + v_R[2][2]) +
            z * (v_R[1][2] + v_R[2][1]) + w * (v_R[2][0] - v_R[0][2])
        );
    // z element in w field
    v_quat.w =
        2.f *
        (
            // v_quat.z = 2.f * (
            x * (v_R[0][2] + v_R[2][0]) + y * (v_R[1][2] + v_R[2][1]) -
            2.f * z * (v_R[0][0] + v_R[1][1]) + w * (v_R[0][1] - v_R[1][0])
        );
    return v_quat;
}

inline __device__ glm::mat3
scale_to_mat(const float3 scale) {
    glm::mat3 S = glm::mat3(1.f);
    S[0][0] = scale.x;
    S[1][1] = scale.y;
    S[2][2] = scale.z;
    return S;
}

// device helper for culling near points
inline __device__ bool clip_near_plane(
    const float3 p, const float *viewmat, float3 &p_view, float thresh
) {
    p_view = transform_4x3(viewmat, p);
    if (p_view.z <= thresh) {
        return true;
    }
    return false;
}



template<typename vec3>
inline __device__ void warpSum3(vec3& val, cg::thread_block_tile<32>& tile){
    val.x = cg::reduce(tile, val.x, cg::plus<float>());
    val.y = cg::reduce(tile, val.y, cg::plus<float>());
    val.z = cg::reduce(tile, val.z, cg::plus<float>());
}

template<typename vec2>
inline __device__ void warpSum2(vec2& val, cg::thread_block_tile<32>& tile){
    val.x = cg::reduce(tile, val.x, cg::plus<float>());
    val.y = cg::reduce(tile, val.y, cg::plus<float>());
}

inline __device__ void warpSum(float& val, cg::thread_block_tile<32>& tile){
    val = cg::reduce(tile, val, cg::plus<float>());
}
