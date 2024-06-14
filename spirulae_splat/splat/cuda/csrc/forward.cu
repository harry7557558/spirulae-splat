#include "forward.cuh"
#include "helpers.cuh"
#include "ch.cuh"
#include <algorithm>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <iostream>
#include <cuda_fp16.h>

namespace cg = cooperative_groups;

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
    float3* __restrict__ axes_v,
    float2* __restrict__ depth_grads
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
    float2 depth_grad = projected_depth_grad(p_view, R, fx, fy);

    // output
    bounds[idx] = {tile_min.x, tile_min.y, tile_max.x, tile_max.y};
    num_tiles_hit[idx] = tile_area;
    positions[idx] = {p_view.x, p_view.y, p_view.z};
    axes_u[idx] = {V0.x, V0.y, V0.z};
    axes_v[idx] = {V1.x, V1.y, V1.z};
    depth_grads[idx] = {depth_grad.x, depth_grad.y};
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


__global__ void rasterize_simple_forward(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float2* __restrict__ anisotropies,
    const float3& __restrict__ background,
    int* __restrict__ final_index,
    float3* __restrict__ out_img,
    float* __restrict__ out_alpha
) {
    // each thread draws one pixel, but also timeshares caching gaussians in a
    // shared tile

    auto block = cg::this_thread_block();
    int32_t tile_id =
        block.group_index().y * tile_bounds.x + block.group_index().x;
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    glm::vec2 pos_screen = { (float)j + 0.5f, (float)i + 0.5f };
    glm::vec2 pos_2d = { (pos_screen.x-cx)/fx, (pos_screen.y-cy)/fy };
    int32_t pix_id = i * img_size.x + j;

    // return if out of bounds
    // keep not rasterizing threads around for reading data
    bool inside = (i < img_size.y && j < img_size.x);
    bool done = !inside;

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 color_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 opacity_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    float T = 1.f;
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    int tr = block.thread_rank();
    float3 pix_out = {0.f, 0.f, 0.f};
    for (int b = 0; b < num_batches; ++b) {
        // resync all threads before beginning next batch
        // end early if entire tile is done
        if (__syncthreads_count(done) >= block_size) {
            break;
        }

        // each thread fetch 1 gaussian from front to back
        // index of gaussian to load
        int batch_start = range.x + block_size * b;
        int idx = batch_start + tr;
        if (idx < range.y) {
            int32_t g_id = gaussian_ids_sorted[idx];
            const float3 pos = positions[g_id];
            const float opac = opacities[g_id];
            const float2 aniso = anisotropies[g_id];
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = {aniso.x, aniso.y, opac};
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        int batch_size = min(block_size, range.y - batch_start);
        for (int t = 0; (t < batch_size) && !done; ++t) {
            glm::vec3 pos = position_batch[t];
            glm::vec2 aniso = {opacity_batch[t].x, opacity_batch[t].y};
            float opac = opacity_batch[t].z;
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            // if (!get_intersection(pos, axis_uv, pos_2d, poi, uv));
            //     continue;
            get_intersection(pos, axis_uv, pos_2d, poi, uv);
            if (glm::length(uv) > visibility_kernel_radius())
                continue;
            float alpha;
            if (!get_alpha(uv, opac, aniso, alpha))
                continue;

            const float next_T = T * (1.f - alpha);

            glm::vec3 color = color_batch[t];

            const float vis = alpha * T;
            pix_out.x = pix_out.x + color.x * vis;
            pix_out.y = pix_out.y + color.y * vis;
            pix_out.z = pix_out.z + color.z * vis;
            T = next_T;
            cur_idx = batch_start + t;
            if (T <= 1e-3f) {
                done = true;
                break;
            }
        }
    }

    if (inside) {
        final_index[pix_id] = cur_idx;
        float3 final_color;
        final_color.x = pix_out.x + T * background.x;
        final_color.y = pix_out.y + T * background.y;
        final_color.z = pix_out.z + T * background.z;
        out_img[pix_id] = final_color;
        out_alpha[pix_id] = 1.0f - T;
    }
}

__global__ void rasterize_depth_forward(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    const float2* __restrict__ anisotropies,
    int* __restrict__ final_index,
    float* __restrict__ out_depth,
    float2* __restrict__ out_visibility
) {
    // each thread draws one pixel, but also timeshares caching gaussians in a
    // shared tile

    auto block = cg::this_thread_block();
    int32_t tile_id =
        block.group_index().y * tile_bounds.x + block.group_index().x;
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    glm::vec2 pos_screen = { (float)j + 0.5f, (float)i + 0.5f };
    glm::vec2 pos_2d = { (pos_screen.x-cx)/fx, (pos_screen.y-cy)/fy };
    int32_t pix_id = i * img_size.x + j;

    // return if out of bounds
    // keep not rasterizing threads around for reading data
    bool inside = (i < img_size.y && j < img_size.x);
    bool done = !inside;

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 opacity_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    float T = 1.f;
    float interp = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    int tr = block.thread_rank();
    float median_depth = 0.0f;
    for (int b = 0; b < num_batches; ++b) {
        // resync all threads before beginning next batch
        // end early if entire tile is done
        if (__syncthreads_count(done) >= block_size) {
            break;
        }

        // each thread fetch 1 gaussian from front to back
        // index of gaussian to load
        int batch_start = range.x + block_size * b;
        int idx = batch_start + tr;
        if (idx < range.y) {
            int32_t g_id = gaussian_ids_sorted[idx];
            const float3 pos = positions[g_id];
            const float opac = opacities[g_id];
            const float2 aniso = anisotropies[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            opacity_batch[tr] = {aniso.x, aniso.y, opac};
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        int batch_size = min(block_size, range.y - batch_start);
        for (int t = 0; (t < batch_size) && !done; ++t) {
            glm::vec3 pos = position_batch[t];
            glm::vec2 aniso = {opacity_batch[t].x, opacity_batch[t].y};
            float opac = opacity_batch[t].z;
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            // if (!get_intersection(pos, axis_uv, pos_2d, poi, uv));
            //     continue;
            get_intersection(pos, axis_uv, pos_2d, poi, uv);
            if (glm::length(uv) > visibility_kernel_radius())
                continue;
            float alpha;
            if (!get_alpha(uv, opac, aniso, alpha))
                continue;

            const float next_T = T * (1.f - alpha);
            // const float next_depth = pos.z;
            const float next_depth = poi.z;
            if (next_T < 0.5f) {
                if (T < 0.99999f) {
                    // https://www.desmos.com/3d/4kuwygxuio
                    interp = (1.0f-alpha)/alpha * (2.0f*T-1.0f);
                    interp = glm::clamp(interp, 0.0f, 1.0f);
                    median_depth = median_depth + (next_depth-median_depth)*interp;
                }
                else {
                    median_depth = next_depth;
                }
                T = next_T;
                cur_idx = batch_start + t;
                done = true;
                break;
            }
            median_depth = next_depth;
            T = next_T;
            cur_idx = batch_start + t;
        }
    }

    if (inside) {
        final_index[pix_id] = cur_idx;
        out_depth[pix_id] = median_depth;
        out_visibility[pix_id] = {T, interp};
    }
}

__global__ void rasterize_forward(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const unsigned ch_degree_r,
    const unsigned ch_degree_phi,
    const float* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    const float2* __restrict__ anisotropies,
    const float3& __restrict__ background,
    const float2* __restrict__ depth_grads,
    const float3* __restrict__ depth_ref_im,
    int* __restrict__ final_index,
    float* __restrict__ out_alpha,
    float3* __restrict__ out_img,
    float4* __restrict__ out_depth_grad,
    float* __restrict__ out_reg_depth,
    float* __restrict__ out_reg_normal
) {
    // each thread draws one pixel, but also timeshares caching gaussians in a
    // shared tile

    auto block = cg::this_thread_block();
    int32_t tile_id =
        block.group_index().y * tile_bounds.x + block.group_index().x;
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    glm::vec2 pos_screen = { (float)j + 0.5f, (float)i + 0.5f };
    glm::vec2 pos_2d = { (pos_screen.x-cx)/fx, (pos_screen.y-cy)/fy };
    int32_t pix_id = i * img_size.x + j;

    // return if out of bounds
    // keep not rasterizing threads around for reading data
    bool inside = (i < img_size.y && j < img_size.x);
    bool done = !inside;

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ int32_t id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 color_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 opacity_batch[MAX_BLOCK_SIZE];
    // __shared__ glm::vec2 depth_grad_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    int tr = block.thread_rank();
    float T = 1.f;  // current/total visibility
    float2 g_sum = {0.f, 0.f};  // sum of "normals"
    float3 pix_out = {0.f, 0.f, 0.f};  // output radiance
    float vis_sum = 0.f;  // output alpha
    float depth_sum = 0.f;  // output depth
    float depth_squared_sum = 0.f;  // for L2 depth regularizer
    const float3 depth_ref_raw = inside ?
        depth_ref_im[pix_id] : make_float3(0.f, 0.f, 0.f);
    const float2 depth_normal_ref = {depth_ref_raw.x, depth_ref_raw.y};
    const float depth_ref = depth_ref_raw.z;
    float reg_depth = 0.f;  // output depth regularizer
    float reg_normal = 0.f;  // output normal regularizer
    for (int b = 0; b < num_batches; ++b) {
        // resync all threads before beginning next batch
        // end early if entire tile is done
        if (__syncthreads_count(done) >= block_size) {
            break;
        }
        // each thread fetch 1 gaussian from front to back
        // index of gaussian to load
        int batch_start = range.x + block_size * b;
        int idx = batch_start + tr;
        if (idx < range.y) {
            int32_t g_id = gaussian_ids_sorted[idx];
            id_batch[tr] = g_id;
            const float3 pos = positions[g_id];
            const float opac = opacities[g_id];
            const float2 aniso = anisotropies[g_id];
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            // const float2 depth_grad = depth_grads[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = {aniso.x, aniso.y, opac};
            // depth_grad_batch[tr] = {depth_grad.x, depth_grad.y};
        }
        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        int batch_size = min(block_size, range.y - batch_start);
        for (int t = 0; (t < batch_size) && !done; ++t) {
            glm::vec3 pos = position_batch[t];
            glm::vec2 aniso = {opacity_batch[t].x, opacity_batch[t].y};
            float opac = opacity_batch[t].z;
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            // if (!get_intersection(pos, axis_uv, pos_2d, poi, uv));
            //     continue;
            get_intersection(pos, axis_uv, pos_2d, poi, uv);
            if (glm::length(uv) > visibility_kernel_radius())
                continue;
            float alpha;
            if (!get_alpha(uv, opac, aniso, alpha))
                continue;

            const float next_T = T * (1.f - alpha);

            glm::vec3 color_0 = color_batch[t];
            glm::vec3 color;
            if (dim_ch > 0) {
                int32_t g_id = id_batch[t];
                const float* coeffs = &ch_coeffs[3*dim_ch*g_id];
                glm::vec3 ch_color;
                ch_coeffs_to_color(
                    ch_degree_r, ch_degree_phi,
                    coeffs, {uv.x, uv.y}, &ch_color.x
                );
                color = color_0 / (1.0f+glm::exp(-ch_color));
            }
            else color = color_0;

            const float vis = alpha * T;
            #if DEPTH_REG_L == 01
            const float depth = pos.z;
            #else
            const float depth = poi.z;
            #endif
            const glm::vec2 g_i = *(glm::vec2*)&depth_grads[id_batch[t]];
            const float g_i_norm = glm::length(g_i) + 1e-6f;
            const glm::vec2 n_i = g_i / g_i_norm;

            pix_out.x = pix_out.x + color.x * vis;
            pix_out.y = pix_out.y + color.y * vis;
            pix_out.z = pix_out.z + color.z * vis;
            #if DEPTH_REG_L == 01
            reg_depth += vis*depth * vis_sum - vis * depth_sum;
            #elif DEPTH_REG_L == 02
            reg_depth += vis * (vis_sum*depth*depth + depth_squared_sum - 2.0f*depth*depth_sum);
            #elif DEPTH_REG_L == 11
            reg_depth += vis * abs(depth - depth_ref);
            #elif DEPTH_REG_L == 12
            reg_depth += vis * (depth-depth_ref) * (depth-depth_ref);
            #endif
            reg_normal += vis * (1.0f - (n_i.x*depth_normal_ref.x+n_i.y*depth_normal_ref.y));
            vis_sum += vis;
            depth_sum += vis*depth;
            depth_squared_sum += vis*depth*depth;
            g_sum.x = g_sum.x + vis * g_i.x;
            g_sum.y = g_sum.y + vis * g_i.y;

            T = next_T;
            cur_idx = batch_start + t;
            if (T <= 1e-3f) {
                done = true;
                break;
            }
        }
    }

    if (inside) {
        final_index[pix_id] = cur_idx;
        out_alpha[pix_id] = 1.0f - T;
        float3 final_color;
        final_color.x = pix_out.x + T * background.x;
        final_color.y = pix_out.y + T * background.y;
        final_color.z = pix_out.z + T * background.z;
        out_img[pix_id] = final_color;
        out_depth_grad[pix_id] = {g_sum.x, g_sum.y, depth_sum, depth_squared_sum};
        out_reg_normal[pix_id] = reg_normal;
        out_reg_depth[pix_id] = reg_depth;
    }
}

// device helper to get screen space depth gradient
__device__ float2 projected_depth_grad(
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
