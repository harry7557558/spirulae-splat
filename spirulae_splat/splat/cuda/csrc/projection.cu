#include "helpers.cuh"
#include <algorithm>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;




// device helper to get screen space depth gradient
inline __device__ float2 projected_depth_grad(
    const glm::vec3 p, const glm::mat3 R,
    const float fx, const float fy
) {
    glm::vec3 n1 = R[2];
    glm::mat3 invJ = glm::mat3(
        p.z/fx, 0.0f, 0.0f,
        0.0f, p.z/fy, 0.0f,
        p.x/p.z, p.y/p.z, 1.0f
    );
    glm::vec3 n = glm::transpose(invJ) * n1;
    // n.z = safe_denom(n.z, 1e-3f);
    // return { -n.x/n.z, -n.y/n.z };
    return {-n.x*n.z, -n.y*n.z};
}

inline __device__ void projected_depth_grad_vjp(
    const glm::vec3 p, const glm::mat3 R,
    const float fx, const float fy,
    const float2 v_depth_grad,
    glm::vec3 &v_p_view, glm::mat3 &v_R
) {
    // forward
    glm::vec3 n1 = R[2];
    glm::mat3 invJ = glm::mat3(
        p.z/fx, 0.0f, 0.0f,
        0.0f, p.z/fy, 0.0f,
        p.x/p.z, p.y/p.z, 1.0f
    );
    glm::vec3 n = glm::transpose(invJ) * n1;
    // n.z = safe_denom(n.z, 1e-2f);
    // glm::vec2 depth_grad = glm::vec2(-n.x/n.z, -n.y/n.z);
    glm::vec2 depth_grad = glm::vec2(-n.x*n.z, -n.y*n.z);

    // backward
    glm::vec3 v_n = glm::vec3(
        // -1.0f/n.z * v_depth_grad.x,
        // -1.0f/n.z * v_depth_grad.y,
        // (n.x*v_depth_grad.x + n.y*v_depth_grad.y) / safe_denom(n.z*n.z,1e-2f)
        -n.z * v_depth_grad.x,
        -n.z * v_depth_grad.y,
        -(n.x * v_depth_grad.x + n.y*v_depth_grad.y)
    );
    // rotation
    glm::vec3 v_n1 = invJ * v_n;
    v_R = glm::mat3(0.0);
    v_R[2] = v_n1;
    // view
    glm::mat3 v_invJ = glm::outerProduct(v_n, n1);
    glm::vec3 v_p = glm::vec3(
        v_invJ[0][2] / p.z,
        v_invJ[1][2] / p.z,
        v_invJ[0][0]/fx + v_invJ[1][1]/fy -
        (p.x*v_invJ[0][2]+p.y*v_invJ[1][2]) / safe_denom(p.z*p.z, 1e-2f)
    );
    v_p_view = {v_p.x, v_p.y, v_p.z};
}



// kernel function for projecting each gaussian on device
// each thread processes one gaussian
__global__ void project_gaussians_forward_kernel(
    const int num_points,
    const float3* __restrict__ means3d,
    const float2* __restrict__ scales,
    const float4* __restrict__ quats,
    const float* __restrict__ viewmat,
    const float4 intrins,
    const dim3 tile_bounds,
    const unsigned block_width,
    const float clip_thresh,
    int4* __restrict__ bounds,
    int32_t* __restrict__ num_tiles_hit,
    float3* __restrict__ positions,
    float3* __restrict__ axes_u,
    float3* __restrict__ axes_v
    // float2* __restrict__ depth_grads
) {
    unsigned idx = cg::this_grid().thread_rank(); // idx of thread within grid
    if (idx >= num_points) {
        return;
    }
    bounds[idx] = {0, 0, 0, 0};
    num_tiles_hit[idx] = 0;

    glm::mat3 R0 = glm::mat3(
        viewmat[0], viewmat[4], viewmat[8],
        viewmat[1], viewmat[5], viewmat[9],
        viewmat[2], viewmat[6], viewmat[10]
    );
    glm::vec3 T0 = { viewmat[3], viewmat[7], viewmat[11] };

    // world to view
    glm::vec3 p_world = *(glm::vec3*)&means3d[idx];
    glm::vec3 p_view = R0 * p_world + T0;
    if (!(p_view.z >= clip_thresh))
        return;

    // patch orientation
    float2 scale = scales[idx];
    float4 quat = quats[idx];
    glm::mat3 Rq = quat_to_rotmat(quat);
    glm::mat3 R = R0 * Rq;
    glm::vec3 V0 = scale.x * R[0];
    glm::vec3 V1 = scale.y * R[1];

    // project to 2d
    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    float2 center;
    float3 bound;
    const float kr = visibility_kernel_radius();
    project_ellipse_bound(p_view, kr*V0, kr*V1, fx, fy, cx, cy, center, bound);

    // compute the projected area
    int2 tile_min, tile_max;
    get_tile_bbox(center, {bound.x,bound.y}, tile_bounds, tile_min, tile_max, block_width);
    int32_t tile_area = (tile_max.x - tile_min.x) * (tile_max.y - tile_min.y);
    if (tile_area <= 0)
        return;

    // compute the depth gradient
    // float2 depth_grad = projected_depth_grad(p_view, R, fx, fy);

    // output
    bounds[idx] = {tile_min.x, tile_min.y, tile_max.x, tile_max.y};
    num_tiles_hit[idx] = tile_area;
    positions[idx] = {p_view.x, p_view.y, p_view.z};
    axes_u[idx] = {V0.x, V0.y, V0.z};
    axes_v[idx] = {V1.x, V1.y, V1.z};
    // depth_grads[idx] = {depth_grad.x, depth_grad.y};
}


__global__ void project_gaussians_backward_kernel(
    const int num_points,
    const float3* __restrict__ means3d,
    const float2* __restrict__ scales,
    const float4* __restrict__ quats,
    const float* __restrict__ viewmat,  // 3x4 row major
    const float4 intrins,
    const int* __restrict__ num_tiles_hit,
    const float3* __restrict__ v_positions,
    const float3* __restrict__ v_axes_u,
    const float3* __restrict__ v_axes_v,
    // const float2* __restrict__ v_depth_grads,
    float3* __restrict__ v_means3d,
    float2* __restrict__ v_scales,
    float4* __restrict__ v_quats,
    float* __restrict__ v_viewmat
) {
    unsigned idx = cg::this_grid().thread_rank(); // idx of thread within grid
    if (idx >= num_points || num_tiles_hit <= 0) {
        return;
    }

    glm::mat3 R0 = glm::mat3(
        viewmat[0], viewmat[4], viewmat[8],
        viewmat[1], viewmat[5], viewmat[9],
        viewmat[2], viewmat[6], viewmat[10]
    );
    glm::vec3 T0 = { viewmat[3], viewmat[7], viewmat[11] };

    // position
    glm::vec3 p_world = *(glm::vec3*)&means3d[idx];
    float fx = intrins.x;
    float fy = intrins.y;
    glm::vec3 p_view = R0 * p_world + T0;
    glm::vec3 v_p_view = *(glm::vec3*)&v_positions[idx];

    // forward
    float2 scale = scales[idx];
    float4 quat = quats[idx];
    glm::mat3 Rq = quat_to_rotmat(quat);
    glm::mat3 R = R0 * Rq;
    glm::vec3 V0 = scale.x * R[0];
    glm::vec3 V1 = scale.y * R[1];

    // scale
    glm::vec3 v_V0 = *(glm::vec3*)&v_axes_u[idx];
    glm::vec3 v_V1 = *(glm::vec3*)&v_axes_v[idx];
    float2 v_scale = {glm::dot(R[0], v_V0), glm::dot(R[1], v_V1)};

    // orientation
    glm::mat3 v_R = glm::mat3(0.0f);
    v_R[0] = scale.x * v_V0;
    v_R[1] = scale.y * v_V1;

    // depth_grad
    #if 0
    glm::mat3 v_R_dg;
    glm::vec3 v_p_view_dg = {0.f, 0.f, 0.f};
    projected_depth_grad_vjp(
        p_view, R, fx, fy, v_depth_grads[idx],
        v_p_view_dg, v_R_dg);
    v_R += v_R_dg;
    v_p_view += v_p_view_dg;
    #endif

    glm::vec3 v_p_world = glm::transpose(R0) * v_p_view;
    v_means3d[idx] = {v_p_world.x, v_p_world.y, v_p_world.z};
    v_scales[idx] = v_scale;
    glm::mat3 v_Rq = glm::transpose(R0) * v_R;
    float4 v_quat = quat_to_rotmat_vjp(quat, v_Rq);
    v_quats[idx] = v_quat;
    glm::mat3 v_R0 = v_R * glm::transpose(Rq) +
        glm::outerProduct(v_p_view, p_world);
    glm::vec3 v_T0 = v_p_view;

    atomicAdd(&v_viewmat[0], v_R0[0][0]);
    atomicAdd(&v_viewmat[1], v_R0[1][0]);
    atomicAdd(&v_viewmat[2], v_R0[2][0]);
    atomicAdd(&v_viewmat[3], v_T0[0]);
    atomicAdd(&v_viewmat[4], v_R0[0][1]);
    atomicAdd(&v_viewmat[5], v_R0[1][1]);
    atomicAdd(&v_viewmat[6], v_R0[2][1]);
    atomicAdd(&v_viewmat[7], v_T0[1]);
    atomicAdd(&v_viewmat[8], v_R0[0][2]);
    atomicAdd(&v_viewmat[9], v_R0[1][2]);
    atomicAdd(&v_viewmat[10], v_R0[2][2]);
    atomicAdd(&v_viewmat[11], v_T0[2]);
}



// kernel to map each intersection from tile ID and depth to a gaussian
// writes output to isect_ids and gaussian_ids
__global__ void map_gaussian_to_intersects(
    const int num_points,
    const float3* __restrict__ positions,
    int4* __restrict__ bounds,
    const int32_t* __restrict__ cum_tiles_hit,
    const dim3 tile_bounds,
    const unsigned block_width,
    int64_t* __restrict__ isect_ids,
    int32_t* __restrict__ gaussian_ids
) {
    unsigned idx = cg::this_grid().thread_rank();
    if (idx >= num_points)
        return;
    int4 bound = bounds[idx];
    if (min(bound.z-bound.x, bound.w-bound.y) <= 0)
        return;

    // update the intersection info for all tiles this gaussian hits
    int32_t cur_idx = (idx == 0) ? 0 : cum_tiles_hit[idx - 1];
    // printf("point %d starting at %d\n", idx, cur_idx);
    int64_t depth_id = (int64_t) * (int32_t *)&(positions[idx].z);
    for (int i = bound.y; i < bound.w; ++i) {
        for (int j = bound.x; j < bound.z; ++j) {
            // isect_id is tile ID and depth as int32
            int64_t tile_id = i * tile_bounds.x + j; // tile within image
            isect_ids[cur_idx] = (tile_id << 32) | depth_id; // tile | depth id
            gaussian_ids[cur_idx] = idx;                     // 3D gaussian id
            ++cur_idx; // handles gaussians that hit more than one tile
        }
    }
    // printf("point %d ending at %d\n", idx, cur_idx);
}


// kernel to map sorted intersection IDs to tile bins
// expect that intersection IDs are sorted by increasing tile ID
// i.e. intersections of a tile are in contiguous chunks
__global__ void get_tile_bin_edges(
    const int num_intersects, const int64_t* __restrict__ isect_ids_sorted, int2* __restrict__ tile_bins
) {
    unsigned idx = cg::this_grid().thread_rank();
    if (idx >= num_intersects)
        return;
    // save the indices where the tile_id changes
    int32_t cur_tile_idx = (int32_t)(isect_ids_sorted[idx] >> 32);
    if (idx == 0 || idx == num_intersects - 1) {
        if (idx == 0)
            tile_bins[cur_tile_idx].x = 0;
        if (idx == num_intersects - 1)
            tile_bins[cur_tile_idx].y = num_intersects;
    }
    if (idx == 0)
        return;
    int32_t prev_tile_idx = (int32_t)(isect_ids_sorted[idx - 1] >> 32);
    if (prev_tile_idx != cur_tile_idx) {
        tile_bins[prev_tile_idx].y = idx;
        tile_bins[cur_tile_idx].x = idx;
        return;
    }
}


// Compute relocation for MCMC update
__global__ void compute_relocation_kernel(
    const int num_points,
    const int num_scales,
    const float *opacities,
    const float *scales,
    const int *ratios,
    float *binoms,
    int n_max,
    float *new_opacities,
    float *new_scales
) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= num_points)
        return;

    int n_idx = ratios[idx];
    float denom_sum = 0.0f;

    float opac = opacities[idx];

#if 0
    // for Gaussian, minimize residual
    float new_opac = 1.0f - powf(1.0f - opac, 1.0f / n_idx);
    for (int i = 1; i <= n_idx; ++i) {
        for (int k = 0; k <= (i - 1); ++k) {
            float bin_coeff = binoms[(i - 1) * n_max + k];
            float term =
                (pow(-1, k) / sqrt(k + 1)) * pow(new_opac, k + 1);
            denom_sum += (bin_coeff * term);
        }
    }
    float sc = (opacities[idx] / denom_sum);
#else
    // for 1-r^2 splats, strictly less than

    // keep opacity, change scale
    //float new_opac = 1.0f - powf(1.0f - opac, 1.0f / n_idx);
    //float sc = sqrtf(n_idx/opac * new_opac * powf(1.0f-new_opac, n_idx-1.0f));

    // keep scale, change opacity
    float new_opac = opac / n_idx;
    float sc = 1.0f;

#endif

    new_opacities[idx] = new_opac;
    for (int i = 0; i < num_scales; ++i)
        new_scales[idx*num_scales+i] = sc * scales[idx*num_scales+i];
}


__global__ void compute_relocation_split_kernel(
    const int num_points,
    const float3 *positions,
    const float4 *quats,
    const float *opacities,
    const float2 *scales,
    float3 *new_position_offsets,
    float *new_opacities,
    float2 *new_scales
) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= num_points)
        return;

    // split into 2 splats, compute position offset (plus/minus)
    // practically the effect of drawing more than 2 splats should be minimal

    // https://www.desmos.com/3d/0ov5xkvbqx

    const float a = 0.25f;

    float opac = opacities[idx];

    const float cx0 = 1.12234f;
    const float cx1 = 0.394553f;
    const float cx2 = -0.120572f;
    float sx1 = cx0 + opac * (cx1 + cx2 * opac);

    const float cy0 = 0.986738f;
    const float cy1 = 0.0304732f;
    const float cy2 = 0.0145201f;
    float sy1 = cy0 + opac * (cy1 + cy2 * opac);

    const float co1 = 0.554027f;
    const float co2 = 0.295756f;
    new_opacities[idx] = opac * (co1 + co2 * opac);

    float4 quat = quats[idx];
    glm::mat3 R = quat_to_rotmat(quat);

    float2 scale = scales[idx];
    if (scale.x > scale.y) {
        glm::vec3 offset = a * scale.x * R[0];
        new_position_offsets[idx] = { offset.x, offset.y, offset.z };
        new_scales[idx] = { scale.x / sx1, scale.y / sy1 };
    }
    else {
        glm::vec3 offset = a * scale.y * R[1];
        new_position_offsets[idx] = { offset.x, offset.y, offset.z };
        new_scales[idx] = { scale.x / sy1, scale.y / sx1 };
    }

}
