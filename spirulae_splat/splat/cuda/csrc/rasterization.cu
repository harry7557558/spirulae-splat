#include "rasterization.cuh"

#include "helpers.cuh"
#include "camera.cuh"
#include "ch.cuh"
#include <algorithm>



template<CameraType CAMERA_TYPE>
__global__ void rasterize_simple_forward_kernel(
    _ARGS_rasterize_simple_forward_kernel
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

    // camera distortion
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y))
            done = true;
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 color_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

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
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = opac;
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        int batch_size = min(block_size, range.y - batch_start);
        for (int t = 0; (t < batch_size) && !done; ++t) {
            glm::vec3 pos = position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            // if (!get_intersection(pos, axis_uv, pos_2d, poi, uv));
            //     continue;
            get_intersection(pos, axis_uv, pos_2d, poi, uv);
            if (glm::length(uv) > visibility_kernel_radius())
                continue;
            float alpha;
            if (!get_alpha(uv, opac, alpha))
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


__global__ void rasterize_simple_backward_kernel(
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
    const float3 __restrict__ background,
    const int* __restrict__ final_index,
    const float* __restrict__ output_alpha,
    const float3* __restrict__ v_output,
    const float* __restrict__ v_output_alpha,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float* __restrict__ v_opacities
) {
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
    // clamp this value to the last pixel
    const int32_t pix_id = min(i * img_size.x + j, img_size.x * img_size.y - 1);

    // keep not rasterizing threads around for reading data
    const bool inside = (i < img_size.y && j < img_size.x);

    // this is the T AFTER the last gaussian in this pixel
    float T_final = 1.0f - output_alpha[pix_id];
    float T = T_final;
    // the contribution from gaussians behind the current one
    float3 buffer = {0.f, 0.f, 0.f};
    // index of last gaussian to contribute to this pixel
    const int bin_final = inside? final_index[pix_id] : 0;

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    const int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    const int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ int32_t id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 color_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // df/d_out for this pixel
    const float3 v_out = nan_to_num(v_output[pix_id]);
    const float v_out_alpha = nan_to_num(v_output_alpha[pix_id]);

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing
    const int tr = block.thread_rank();
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(block);
    const int warp_bin_final = cg::reduce(warp, bin_final, cg::greater<int>());
    for (int b = 0; b < num_batches; ++b) {
        // resync all threads before writing next batch of shared mem
        block.sync();

        // each thread fetch 1 gaussian from back to front
        // 0 index will be furthest back in batch
        // index of gaussian to load
        // batch end is the index of the last gaussian in the batch
        const int batch_end = range.y - 1 - block_size * b;
        int batch_size = min(block_size, batch_end + 1 - range.x);
        const int idx = batch_end - tr;
        if (idx >= range.x) {
            int32_t g_id = gaussian_ids_sorted[idx];
            id_batch[tr] = g_id;
            const float3 pos = positions[g_id];
            const float opac = opacities[g_id];
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = opac;
        }
        // wait for other threads to collect the gaussians in batch
        block.sync();
        // process gaussians in the current batch for this pixel
        // 0 index is the furthest back gaussian in the batch
        for (int t = max(0,batch_end - warp_bin_final); t < batch_size; ++t) {
            int valid = inside;
            if (batch_end - t > bin_final) {
                valid = 0;
            }

            glm::vec3 pos = position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            if(valid){
                get_intersection(pos, axis_uv, pos_2d, poi, uv);
                if (glm::length(uv) > visibility_kernel_radius())
                    valid = 0;
            }
            float alpha;
            if (valid){
                if (!get_alpha(uv, opac, alpha))
                    valid = 0;
            }
            if(!warp.any(valid)){
                continue;
            }

            glm::vec3 v_position_local = {0.f, 0.f, 0.f};
            glm::vec2 v_position_xy_abs_local = {0.f, 0.f};
            glm::vec3 v_axis_u_local = {0.f, 0.f, 0.f};
            glm::vec3 v_axis_v_local = {0.f, 0.f, 0.f};
            glm::vec3 v_color_local = {0.f, 0.f, 0.f};
            float v_opacity_local = 0.f;
            //initialize everything to 0, only set if the lane is valid
            if(valid){
                // compute the current T for this gaussian
                const float ra = 1.f / (1.f - alpha);
                const float next_T = T * ra;
                const float vis = alpha * next_T;

                // update v_rgb for this gaussian
                float v_alpha = 0.f;
                v_color_local = {vis * v_out.x, vis * v_out.y, vis * v_out.z};

                const glm::vec3 color = color_batch[t];
                const float opacity = opacity_batch[t];
                // contribution from this pixel
                v_alpha += (color.x * T - buffer.x) * ra * v_out.x;
                v_alpha += (color.y * T - buffer.y) * ra * v_out.y;
                v_alpha += (color.z * T - buffer.z) * ra * v_out.z;

                v_alpha += T_final * ra * v_out_alpha;
                // contribution from background pixel
                v_alpha += -T_final * ra * background.x * v_out.x;
                v_alpha += -T_final * ra * background.y * v_out.y;
                v_alpha += -T_final * ra * background.z * v_out.z;
                // update the running sum
                buffer.x += color.x * vis;
                buffer.y += color.y * vis;
                buffer.z += color.z * vis;

                glm::vec2 v_uv;
                get_alpha_vjp(
                    uv, opacity,
                    v_alpha, v_uv, v_opacity_local
                );
                glm::mat2x3 v_axis_uv = glm::mat2x3(0.0f);
                get_intersection_vjp(
                    pos, axis_uv, pos_2d,
                    glm::vec3(0), v_uv,
                    v_position_local, v_axis_uv
                );
                v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local)) / glm::vec2(fx, fy);
                // v_position_xy_abs_local /= pos.z;
                v_axis_u_local = v_axis_uv[0];
                v_axis_v_local = v_axis_uv[1];

                T = next_T;
            }
            warpSum3(v_position_local, warp);
            warpSum2(v_position_xy_abs_local, warp);
            warpSum3(v_axis_u_local, warp);
            warpSum3(v_axis_v_local, warp);
            warpSum3(v_color_local, warp);
            warpSum(v_opacity_local, warp);
            if (warp.thread_rank() == 0) {
                int32_t g = id_batch[t];

                float* v_position_ptr = (float*)(v_positions);
                atomicAdd(v_position_ptr + 3*g + 0, v_position_local.x);
                atomicAdd(v_position_ptr + 3*g + 1, v_position_local.y);
                atomicAdd(v_position_ptr + 3*g + 2, v_position_local.z);
                float* v_positions_xy_abs_ptr = (float*)(v_positions_xy_abs);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 0, v_position_xy_abs_local.x);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 1, v_position_xy_abs_local.y);

                float* v_axis_u_ptr = (float*)(v_axes_u);
                atomicAdd(v_axis_u_ptr + 3*g + 0, v_axis_u_local.x);
                atomicAdd(v_axis_u_ptr + 3*g + 1, v_axis_u_local.y);
                atomicAdd(v_axis_u_ptr + 3*g + 2, v_axis_u_local.z);
                float* v_axis_v_ptr = (float*)(v_axes_v);
                atomicAdd(v_axis_v_ptr + 3*g + 0, v_axis_v_local.x);
                atomicAdd(v_axis_v_ptr + 3*g + 1, v_axis_v_local.y);
                atomicAdd(v_axis_v_ptr + 3*g + 2, v_axis_v_local.z);
                
                float* v_color_ptr = (float*)(v_colors);
                atomicAdd(v_color_ptr + 3*g + 0, v_color_local.x);
                atomicAdd(v_color_ptr + 3*g + 1, v_color_local.y);
                atomicAdd(v_color_ptr + 3*g + 2, v_color_local.z);
                
                atomicAdd(v_opacities + g, v_opacity_local);
            }
        }
    }
}





template <DepthMode DEPTH_MODE, CameraType CAMERA_TYPE>
__global__ void rasterize_depth_forward_kernel(
    _ARGS_rasterize_depth_forward_kernel
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

    // camera distortion
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y))
            done = true;
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    float T = 1.f;
    float interp = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    int tr = block.thread_rank();
    float output_depth = 0.0f;
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
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            opacity_batch[tr] = opac;
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        int batch_size = min(block_size, range.y - batch_start);
        for (int t = 0; (t < batch_size) && !done; ++t) {
            glm::vec3 pos = (glm::vec3)position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = (glm::mat2x3)axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            // if (!get_intersection(pos, axis_uv, pos_2d, poi, uv));
            //     continue;
            get_intersection(pos, axis_uv, pos_2d, poi, uv);
            if (glm::length(uv) > visibility_kernel_radius())
                continue;
            float alpha;
            if (!get_alpha(uv, opac, alpha))
                continue;

            const float next_T = T * (1.f - alpha);

            // mean depth
            if (DEPTH_MODE == DepthMode::Mean) {

                // const float depth_raw = pos.z;
                const float depth_raw = poi.z;
                const float depth = depth_map(depth_raw);
                float vis = alpha * T;
                output_depth += vis * depth;

            }  // DEPTH_MODE

            // median depth
            else if (DEPTH_MODE == DepthMode::Median) {

                const float next_depth_raw = poi.z;
                const float next_depth = depth_map(next_depth_raw);
                if (next_T < DEPTH_REG_MEDIAN_TH) {
                    if (T < 0.99999f) {
                        // https://www.desmos.com/3d/fttajoozww
                        interp = (1.0f-alpha)/alpha * (T-DEPTH_REG_MEDIAN_TH)/DEPTH_REG_MEDIAN_TH;
                        interp = glm::clamp(interp, 0.0f, 1.0f);
                        output_depth = output_depth + (next_depth-output_depth)*interp;
                    }
                    else {
                        output_depth = next_depth;
                    }
                    T = next_T;
                    cur_idx = batch_start + t;
                    done = true;
                    break;
                }
                output_depth = next_depth;

            }  // DEPTH_MODE

            T = next_T;
            cur_idx = batch_start + t;
        }
    }

    if (inside) {
        final_index[pix_id] = cur_idx;
        if (DEPTH_MODE == DepthMode::Mean) {
            float depth = T == 1.0f ? output_depth : output_depth / (1.0f-T);
            // out_depth[pix_id] = depth_inv_map(depth);
            out_depth[pix_id] = depth;
            out_visibility[pix_id] = {T, 1.0f-T};
        }
        else if (DEPTH_MODE == DepthMode::Median) {
            // out_depth[pix_id] = depth_inv_map(output_depth);
            out_depth[pix_id] = output_depth;
            out_visibility[pix_id] = {T, interp};
        }
    }
}


template <DepthMode DEPTH_MODE>
__global__ void rasterize_depth_backward_kernel(
    _ARGS_rasterize_depth_backward_kernel
) {
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
    // clamp this value to the last pixel
    const int32_t pix_id = min(i * img_size.x + j, img_size.x * img_size.y - 1);

    // keep not rasterizing threads around for reading data
    const bool inside = (i < img_size.y && j < img_size.x);

    // this is the T AFTER the last gaussian in this pixel
    glm::vec2 meta_out = *(glm::vec2*)&out_visibility[pix_id];
    float T_final = meta_out.x;
    float T = T_final;
    float v_T = 0.0f;
    const float interp = meta_out.y;
    // index of last gaussian to contribute to this pixel
    const int bin_final = inside? final_index[pix_id] : 0;

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    const int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    const int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ int32_t id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // df/d_out for this pixel
    float output_depth = out_depth[pix_id];
    float v_output_depth = nan_to_num(v_out_depth[pix_id]);
    float v_out_alpha = 0.0f;
    float v_depth = 0.f;
    float v_depth_next = 0.f;
    float v_alpha = 0.f;
    float v_interp = 0.f;
    if (DEPTH_MODE == DepthMode::Mean) {
        if (T != 1.0f) {
            float alpha = 1.0f-T;
            v_out_alpha = -output_depth / fmax(alpha, 1e-4f) * v_output_depth;
            output_depth *= alpha;
            v_output_depth /= fmax(alpha, 1e-4);
        }
    }

    float depth_buffer = 0.0f;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing
    const int tr = block.thread_rank();
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(block);
    const int warp_bin_final = cg::reduce(warp, bin_final, cg::greater<int>());
    for (int b = 0; b < num_batches; ++b) {
        // resync all threads before writing next batch of shared mem
        block.sync();

        // each thread fetch 1 gaussian from back to front
        // 0 index will be furthest back in batch
        // index of gaussian to load
        // batch end is the index of the last gaussian in the batch
        const int batch_end = range.y - 1 - block_size * b;
        int batch_size = min(block_size, batch_end + 1 - range.x);
        const int idx = batch_end - tr;
        if (idx >= range.x) {
            int32_t g_id = gaussian_ids_sorted[idx];
            id_batch[tr] = g_id;
            const float3 pos = positions[g_id];
            const float opac = opacities[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            opacity_batch[tr] = opac;
        }
        // wait for other threads to collect the gaussians in batch
        block.sync();
        // process gaussians in the current batch for this pixel
        // 0 index is the furthest back gaussian in the batch
        for (int t = max(0,batch_end - warp_bin_final); t < batch_size; ++t) {
            int valid = inside;
            if (batch_end - t > bin_final) {
                valid = 0;
            }

            glm::vec3 pos = (glm::vec3)position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = (glm::mat2x3)axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            if(valid){
                get_intersection(pos, axis_uv, pos_2d, poi, uv);
                if (glm::length(uv) > visibility_kernel_radius())
                    valid = 0;
            }
            float alpha;
            if (valid){
                if (!get_alpha(uv, opac, alpha))
                    valid = 0;
            }
            if(!warp.any(valid)){
                continue;
            }

            glm::vec3 v_position_local = {0.f, 0.f, 0.f};
            glm::vec2 v_position_xy_abs_local = {0.f, 0.f};
            glm::vec3 v_axis_u_local = {0.f, 0.f, 0.f};
            glm::vec3 v_axis_v_local = {0.f, 0.f, 0.f};
            float v_opacity_local = 0.f;
            //initialize everything to 0, only set if the lane is valid
            if(valid) {
                // compute the current T for this gaussian
                const float ra = 1.f / (1.f - alpha);
                const float next_T = T * ra;
                const float vis = alpha * next_T;

                float depth_raw = poi.z;
                float depth = depth_map(depth_raw);

                // mean depth
                if (DEPTH_MODE == DepthMode::Mean) {

                    v_depth = vis * v_output_depth;
                    v_alpha = (depth * T - depth_buffer) * ra * v_output_depth;
                    v_alpha += T_final * ra * v_out_alpha;
                    depth_buffer += depth * vis;

                }  // DEPTH_MODE

                // median depth
                else if (DEPTH_MODE == DepthMode::Median) {

                    // depth gradient
                    if (T == T_final) {
                        v_depth = v_output_depth * interp;
                        v_depth_next = v_output_depth * (1.0f-interp);
                    }
                    else {
                        v_depth = v_depth_next;
                        v_depth_next = 0.0f;
                    }

                    // alpha gradient
                    if (T == T_final && interp < 1.0f && interp > 0.0f) {
                        float depth_0 = (output_depth-depth*interp) / (1.0f-interp);
                        v_interp = (depth-depth_0) * v_output_depth;
                        v_alpha = (next_T-DEPTH_REG_MEDIAN_TH)/DEPTH_REG_MEDIAN_TH * \
                            v_interp / safe_denom(-alpha*alpha, 1e-3);
                        v_T = (1.0f-alpha)/alpha * v_interp / DEPTH_REG_MEDIAN_TH;
                    }
                    else {
                        v_alpha = v_T * (-next_T);
                        v_T = v_T * (1.0f-alpha);
                    }

                }  // DEPTH_MODE

                T = next_T;

                // backward
                const float opacity = opacity_batch[t];
                glm::vec2 v_uv;
                get_alpha_vjp(
                    uv, opacity,
                    v_alpha, v_uv, v_opacity_local
                );
                glm::mat2x3 v_axis_uv = glm::mat2x3(0.0f);
                float v_depth_raw = depth_map_vjp(depth_raw, v_depth);
                get_intersection_vjp(
                    pos, axis_uv, pos_2d,
                    {0.f, 0.f, v_depth_raw}, v_uv,
                    v_position_local, v_axis_uv
                );
                v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local)) / glm::vec2(fx, fy);
                v_axis_u_local = v_axis_uv[0];
                v_axis_v_local = v_axis_uv[1];
            }
            warpSum3(v_position_local, warp);
            warpSum2(v_position_xy_abs_local, warp);
            warpSum3(v_axis_u_local, warp);
            warpSum3(v_axis_v_local, warp);
            warpSum(v_opacity_local, warp);
            if (warp.thread_rank() == 0) {
                int32_t g = id_batch[t];

                float* v_position_ptr = (float*)(v_positions);
                atomicAdd(v_position_ptr + 3*g + 0, v_position_local.x);
                atomicAdd(v_position_ptr + 3*g + 1, v_position_local.y);
                atomicAdd(v_position_ptr + 3*g + 2, v_position_local.z);
                float* v_positions_xy_abs_ptr = (float*)(v_positions_xy_abs);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 0, v_position_xy_abs_local.x);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 1, v_position_xy_abs_local.y);

                float* v_axis_u_ptr = (float*)(v_axes_u);
                atomicAdd(v_axis_u_ptr + 3*g + 0, v_axis_u_local.x);
                atomicAdd(v_axis_u_ptr + 3*g + 1, v_axis_u_local.y);
                atomicAdd(v_axis_u_ptr + 3*g + 2, v_axis_u_local.z);
                float* v_axis_v_ptr = (float*)(v_axes_v);
                atomicAdd(v_axis_v_ptr + 3*g + 0, v_axis_v_local.x);
                atomicAdd(v_axis_v_ptr + 3*g + 1, v_axis_v_local.y);
                atomicAdd(v_axis_v_ptr + 3*g + 2, v_axis_v_local.z);
                
                float v_opacity_local_ = (float)v_opacity_local;
                atomicAdd(v_opacities + g, v_opacity_local_);
            }
        }
    }
}






__global__ void rasterize_forward_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const float depth_reg_pairwise_factor,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const float3* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    // const float3 __restrict__ background,
    const float* __restrict__ depth_ref_im,
    int* __restrict__ final_index,
    float* __restrict__ out_alpha,
    float3* __restrict__ out_img,
    float2* __restrict__ out_depth,
    float3* __restrict__ out_normal,
    float* __restrict__ out_reg_depth
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
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    int tr = block.thread_rank();
    float T = 1.f;  // current/total visibility
    float3 normal_out = {0.f, 0.f, 0.f};  // sum of normals
    float3 pix_out = {0.f, 0.f, 0.f};  // output radiance
    float vis_sum = 0.f;  // output alpha
    float depth_sum = 0.f;  // output depth
    float depth_squared_sum = 0.f;  // for L2 depth regularizer
    const float depth_ref = inside ? depth_ref_im[pix_id] : 0.f;
    float reg_depth_p = 0.f, reg_depth_i = 0.f;  // output depth regularizer
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
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            // const float2 depth_grad = depth_grads[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = opac;
            // depth_grad_batch[tr] = {depth_grad.x, depth_grad.y};
        }
        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        int batch_size = min(block_size, range.y - batch_start);
        for (int t = 0; (t < batch_size) && !done; ++t) {
            glm::vec3 pos = position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            // if (!get_intersection(pos, axis_uv, pos_2d, poi, uv));
            //     continue;
            get_intersection(pos, axis_uv, pos_2d, poi, uv);
            if (glm::length(uv) > visibility_kernel_radius())
                continue;
            float alpha;
            if (!get_alpha(uv, opac, alpha))
                continue;

            const float next_T = T * (1.f - alpha);

            glm::vec3 color_0 = color_batch[t];
            glm::vec3 color;
            if (dim_ch > 0) {
                int32_t g_id = id_batch[t];
                const glm::vec3* coeffs = (glm::vec3*)&ch_coeffs[dim_ch*g_id];
                glm::vec3 ch_color = ch_coeffs_to_color(
                    ch_degree_r, ch_degree_r_to_use,
                    ch_degree_phi, ch_degree_phi_to_use,
                    coeffs, {uv.x, uv.y}
                );
                color = color_0 / (1.0f+glm::exp(-ch_color));
            }
            else color = color_0;

            const float vis = alpha * T;
            #if DEPTH_REG_L == 01 && false
            const float depth_raw = pos.z;
            #else
            const float depth_raw = poi.z;
            #endif
            const float depth = depth_map(depth_raw);

            pix_out.x = pix_out.x + color.x * vis;
            pix_out.y = pix_out.y + color.y * vis;
            pix_out.z = pix_out.z + color.z * vis;

            // depth regularization
            {
                // float pairwise_l1 = vis*depth * vis_sum - vis * depth_sum;  // requires pos.z for depth
                float pairwise_l2 = vis * (vis_sum*depth*depth + depth_squared_sum - 2.0f*depth*depth_sum);
                float intersect_l1 = vis * abs(depth - depth_ref);
                // float intersect_l2 = vis * (depth-depth_ref) * (depth-depth_ref);
                reg_depth_p += pairwise_l2;
                reg_depth_i += intersect_l1;
            }
            vis_sum += vis;
            depth_sum += vis*depth;
            depth_squared_sum += vis*depth*depth;

            // normal regularization
            glm::vec3 normal = get_normal_from_axisuv(axis_uv, poi);
            normal_out.x = normal_out.x + normal.x * vis;
            normal_out.y = normal_out.y + normal.y * vis;
            normal_out.z = normal_out.z + normal.z * vis;

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
        // final_color.x = pix_out.x + T * background.x;
        // final_color.y = pix_out.y + T * background.y;
        // final_color.z = pix_out.z + T * background.z;
        final_color.x = pix_out.x;
        final_color.y = pix_out.y;
        final_color.z = pix_out.z;
        out_img[pix_id] = final_color;
        out_depth[pix_id] = {depth_sum, depth_squared_sum};
        out_normal[pix_id] = normal_out;
        out_reg_depth[pix_id] = reg_depth_i + (reg_depth_p-reg_depth_i) * depth_reg_pairwise_factor;
    }
}


__global__ void rasterize_backward_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const float depth_reg_pairwise_factor,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float3* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    // const float3 __restrict__ background,
    const float* __restrict__ depth_ref_im,
    const int* __restrict__ final_index,
    const float* __restrict__ output_alpha,
    const float2* __restrict__ output_depth,
    const float* __restrict__ v_output_alpha,
    const float3* __restrict__ v_output,
    const float2* __restrict__ v_output_depth,
    const float3* __restrict__ v_output_normal,
    const float* __restrict__ v_output_reg_depth,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float3* __restrict__ v_ch_coeffs,
    // float* __restrict__ v_ch_coeffs_abs,
    float* __restrict__ v_opacities,
    // float3* __restrict__ v_background,
    float* __restrict__ v_depth_ref_im
) {
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
    // clamp this value to the last pixel
    const int32_t pix_id = min(i * img_size.x + j, img_size.x * img_size.y - 1);

    // keep not rasterizing threads around for reading data
    const bool inside = (i < img_size.y && j < img_size.x);

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    const int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    const int num_batches = (range.y - range.x + block_size - 1) / block_size;

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);
    assert(dim_ch <= MAX_CH_FLOAT3);

    __shared__ int32_t id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 color_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];
    // __shared__ glm::vec2 depth_grad_batch[MAX_BLOCK_SIZE];

    // df/d_out for this pixel
    const float2 out_depth = output_depth[pix_id];
    const float3 v_out = nan_to_num(v_output[pix_id]);
    const float2 v_out_depth = nan_to_num(v_output_depth[pix_id]);
    const float3 v_out_normal = nan_to_num(v_output_normal[pix_id]);
    const float v_out_alpha = nan_to_num(v_output_alpha[pix_id]);
    const float v_out_reg_depth = nan_to_num(v_output_reg_depth[pix_id]);
    const float v_reg_depth_p = v_out_reg_depth * depth_reg_pairwise_factor;
    const float v_reg_depth_i = v_out_reg_depth * (1.0f-depth_reg_pairwise_factor);
    const float v_depth_sum = v_out_depth.x;
    const float v_depth_squared_sum = v_out_depth.y;

    // this is the T AFTER the last gaussian in this pixel
    float T_final = 1.0f - output_alpha[pix_id];
    float T = T_final;
    // index of last gaussian to contribute to this pixel
    const int bin_final = inside? final_index[pix_id] : 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing
    const int tr = block.thread_rank();

    // regularization
    const float depth_ref = inside ? depth_ref_im[pix_id] : 0.f;
    float v_depth_ref = 0.f;

    const float vis_sum_final = 1.0f - T_final;
    const float depth_sum_final = out_depth.x;
    const float depth_squared_sum_final = out_depth.y;
    float vis_sum = vis_sum_final;
    float depth_sum = depth_sum_final;
    float depth_squared_sum = depth_squared_sum_final;

    float3 buffer = {0.f, 0.f, 0.f};
    float2 buffer_depth = {0.f, 0.f};
    float3 buffer_normal = {0.f, 0.f, 0.f};
    float buffer_depth_reg = 0.f;
    
    // gradient
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(block);
    const int warp_bin_final = cg::reduce(warp, bin_final, cg::greater<int>());
    T = T_final;
    for (int b = 0; b < num_batches; ++b) {
        // resync all threads before writing next batch of shared mem
        block.sync();

        // each thread fetch 1 gaussian from back to front
        // 0 index will be furthest back in batch
        // index of gaussian to load
        // batch end is the index of the last gaussian in the batch
        const int batch_end = range.y - 1 - block_size * b;
        int batch_size = min(block_size, batch_end + 1 - range.x);
        const int idx = batch_end - tr;
        if (idx >= range.x) {
            int32_t g_id = gaussian_ids_sorted[idx];
            id_batch[tr] = g_id;
            const float3 pos = positions[g_id];
            const float opac = opacities[g_id];
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = opac;
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        // 0 index is the furthest back gaussian in the batch
        for (int t = max(0,batch_end - warp_bin_final); t < batch_size; ++t) {
            int valid = inside;
            if (batch_end - t > bin_final) {
                valid = 0;
            }

            glm::vec3 pos = position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            if(valid){
                get_intersection(pos, axis_uv, pos_2d, poi, uv);
                if (glm::length(uv) > visibility_kernel_radius())
                    valid = 0;
            }
            float alpha;
            if (valid){
                if (!get_alpha(uv, opac, alpha))
                    valid = 0;
            }
            if(!warp.any(valid)){
                continue;
            }

            glm::vec3 v_position_local = {0.f, 0.f, 0.f};
            glm::vec2 v_position_xy_abs_local = {0.f, 0.f};
            glm::vec3 v_axis_u_local = {0.f, 0.f, 0.f};
            glm::vec3 v_axis_v_local = {0.f, 0.f, 0.f};
            glm::vec3 v_color_local = {0.f, 0.f, 0.f};
            float v_opacity_local = 0.f;
            glm::vec3 v_ch_coeff_local[MAX_CH_FLOAT3];
            for (int i = 0; i < dim_ch; i++)
                v_ch_coeff_local[i] = {0.f, 0.f, 0.f};
            float v_ch_coeff_abs_local = 0.f;
            //initialize everything to 0, only set if the lane is valid
            if(valid){
                // compute the current T for this gaussian
                const float ra = 1.f / (1.f - alpha);
                const float next_T = T * ra;
                const float vis = alpha * next_T;

                // update accumulation
                float v_depth = 0.0f;
                #if DEPTH_REG_L == 01 && false
                const float depth_raw = pos.z;
                const float depth = depth_map(depth_raw);
                v_depth += vis * v_depth_sum;
                v_depth += vis * 2.0f*depth * v_depth_squared_sum;
                #else
                const float depth_raw = poi.z;
                const float depth = depth_map(depth_raw);
                v_depth += vis * v_depth_sum;
                v_depth += vis * 2.0f*depth * v_depth_squared_sum;
                #endif

                // update depth regularizer
                float vis_sum_next = vis_sum - vis;
                float depth_sum_next = depth_sum - vis*depth;
                float depth_squared_sum_next = depth_squared_sum - vis*depth*depth;
                #if 0  // pairwise L1, requires pos.z for depth
                v_depth += v_reg_depth_p * vis * (vis_sum_next - (vis_sum_final-vis_sum));
                float reg_depth_i_p = (
                    depth * vis_sum_next - depth_sum_next +
                    (depth_sum_final-depth_sum) - depth * (vis_sum_final-vis_sum)
                );
                v_position_local.z = depth_map_vjp(depth_raw, v_depth);
                #else  // pairwise L2
                v_depth += v_reg_depth_p * vis * 2.0f * (
                    vis_sum_final * depth - depth_sum_final);
                float reg_depth_i_p =
                    vis_sum_final*depth*depth + depth_squared_sum_final
                    - 2.0f*depth*depth_sum_final;
                #endif
                #if 1  // L1 with intersected depth
                float v_z = v_reg_depth_i * vis * glm::sign(depth-depth_ref);
                v_depth += v_z;
                v_depth_ref += (-v_z);
                float reg_depth_i_i = abs(depth-depth_ref);
                #else  // L2 with intersected depth
                float v_z = v_reg_depth_i * vis * 2.0f*(depth-depth_ref);
                v_depth += v_z;
                v_depth_ref += (-v_z);
                float reg_depth_i_i = (depth-depth_ref) * (depth-depth_ref);
                #endif
                float reg_depth_i = reg_depth_i_i + (reg_depth_i_p-reg_depth_i_i) * depth_reg_pairwise_factor;

                float v_depth_raw = depth_map_vjp(depth_raw, v_depth);
                glm::vec3 v_poi = {0.f, 0.f, v_depth_raw};

                // normal regularization
                glm::vec3 v_normal = {vis * v_out_normal.x, vis * v_out_normal.y, vis * v_out_normal.z};
                glm::mat2x3 v_axis_uv; glm::vec3 normal;
                get_normal_from_axisuv_vjp(axis_uv, poi, v_normal, normal, v_axis_uv);

                // update color
                glm::vec3 v_color_1 = {vis * v_out.x, vis * v_out.y, vis * v_out.z};
                const float opacity = opacity_batch[t];
                const glm::vec3 color_0 = color_batch[t];
                glm::vec3 color_1;
                glm::vec2 v_uv_ch = {0.f, 0.f};
                if (dim_ch > 0) {
                    glm::vec3 v_ch_color_sigmoid = v_color_1 * color_0;
                    #if 0
                    int32_t g_id = id_batch[t];
                    glm::vec3 ch_color = ch_coeffs_to_color(
                        ch_degree_r, ch_degree_r_to_use,
                        ch_degree_phi, ch_degree_phi_to_use,
                        (glm::vec3*)&ch_coeffs[dim_ch*g_id], {uv.x, uv.y}
                    );
                    glm::vec3 ch_color_sigmoid = 1.0f / (1.0f+glm::exp(-ch_color));
                    glm::vec3 v_ch_color = v_ch_color_sigmoid * ch_color_sigmoid*(1.0f-ch_color_sigmoid);
                    ch_coeffs_to_color_vjp(
                        ch_degree_r, ch_degree_r_to_use,
                        ch_degree_phi, ch_degree_phi_to_use,
                        (glm::vec3*)&ch_coeffs[dim_ch*g_id],
                        {uv.x, uv.y},
                        v_ch_color,
                        ch_color,
                        v_ch_coeff_local, v_ch_coeff_abs_local,
                        v_uv_ch
                    );
                    #else
                    // makes overall training 0.1x faster
                    int32_t g_id = id_batch[t];
                    glm::vec3 ch_color_sigmoid;
                    ch_coeffs_to_color_sigmoid_vjp(
                        ch_degree_r, ch_degree_r_to_use,
                        ch_degree_phi, ch_degree_phi_to_use,
                        (glm::vec3*)&ch_coeffs[dim_ch*g_id],
                        {uv.x, uv.y},
                        v_ch_color_sigmoid,
                        ch_color_sigmoid,
                        v_ch_coeff_local, v_ch_coeff_abs_local,
                        v_uv_ch
                    );
                    #endif
                    color_1 = color_0 * ch_color_sigmoid;
                    v_color_local = v_color_1 * ch_color_sigmoid;
                }
                else {
                    color_1 = color_0;
                    v_color_local = v_color_1;
                }

                float v_alpha = 0.0f;
                // contribution from this pixel
                v_alpha += (color_1.x * T - buffer.x) * ra * v_out.x;
                v_alpha += (color_1.y * T - buffer.y) * ra * v_out.y;
                v_alpha += (color_1.z * T - buffer.z) * ra * v_out.z;
                v_alpha += T_final * ra * v_out_alpha;
                // v_alpha += -T_final * ra * background.x * v_out.x;
                // v_alpha += -T_final * ra * background.y * v_out.y;
                // v_alpha += -T_final * ra * background.z * v_out.z;
                v_alpha += (depth * T - buffer_depth.x) * ra * v_depth_sum;
                v_alpha += (depth*depth * T - buffer_depth.y) * ra * v_depth_squared_sum;
                v_alpha += (reg_depth_i * T - buffer_depth_reg) * ra * v_out_reg_depth;
                v_alpha += (normal.x * T - buffer_normal.x) * ra * v_out_normal.x;
                v_alpha += (normal.y * T - buffer_normal.y) * ra * v_out_normal.y;
                v_alpha += (normal.z * T - buffer_normal.z) * ra * v_out_normal.z;

                // update the running sum
                buffer.x += color_1.x * vis;
                buffer.y += color_1.y * vis;
                buffer.z += color_1.z * vis;
                buffer_depth.x += depth * vis;
                buffer_depth.y += depth*depth * vis;
                buffer_depth_reg += reg_depth_i * vis;
                buffer_normal.x += normal.x * vis;
                buffer_normal.y += normal.y * vis;
                buffer_normal.z += normal.z * vis;

                // grad
                glm::vec2 v_uv;
                get_alpha_vjp(
                    uv, opacity,
                    v_alpha, v_uv, v_opacity_local
                );
                v_uv += v_uv_ch;
                glm::vec3 v_position_local_temp;
                get_intersection_vjp(
                    pos, axis_uv, pos_2d,
                    v_poi, v_uv,
                    v_position_local_temp, v_axis_uv
                );
                v_position_local += v_position_local_temp;
                v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local)) / glm::vec2(fx, fy);
                v_axis_u_local = v_axis_uv[0];
                v_axis_v_local = v_axis_uv[1];

                // next loop
                T = next_T;
                vis_sum = vis_sum_next;
                depth_sum = depth_sum_next;
                depth_squared_sum = depth_squared_sum_next;
            }
            warpSum3(v_position_local, warp);
            warpSum2(v_position_xy_abs_local, warp);
            warpSum3(v_axis_u_local, warp);
            warpSum3(v_axis_v_local, warp);
            warpSum3(v_color_local, warp);
            for (int i = 0; i < dim_ch; i++)
                warpSum3(v_ch_coeff_local[i], warp);
            warpSum(v_ch_coeff_abs_local, warp);
            warpSum(v_opacity_local, warp);
            if (warp.thread_rank() == 0) {
                int32_t g = id_batch[t];

                float* v_position_ptr = (float*)(v_positions);
                atomicAdd(v_position_ptr + 3*g + 0, v_position_local.x);
                atomicAdd(v_position_ptr + 3*g + 1, v_position_local.y);
                atomicAdd(v_position_ptr + 3*g + 2, v_position_local.z);
                float* v_positions_xy_abs_ptr = (float*)(v_positions_xy_abs);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 0, v_position_xy_abs_local.x);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 1, v_position_xy_abs_local.y);

                float* v_axis_u_ptr = (float*)(v_axes_u);
                atomicAdd(v_axis_u_ptr + 3*g + 0, v_axis_u_local.x);
                atomicAdd(v_axis_u_ptr + 3*g + 1, v_axis_u_local.y);
                atomicAdd(v_axis_u_ptr + 3*g + 2, v_axis_u_local.z);
                float* v_axis_v_ptr = (float*)(v_axes_v);
                atomicAdd(v_axis_v_ptr + 3*g + 0, v_axis_v_local.x);
                atomicAdd(v_axis_v_ptr + 3*g + 1, v_axis_v_local.y);
                atomicAdd(v_axis_v_ptr + 3*g + 2, v_axis_v_local.z);
                
                float* v_color_ptr = (float*)(v_colors);
                atomicAdd(v_color_ptr + 3*g + 0, v_color_local.x);
                atomicAdd(v_color_ptr + 3*g + 1, v_color_local.y);
                atomicAdd(v_color_ptr + 3*g + 2, v_color_local.z);
                float* v_ch_coeffs_ptr = (float*)(v_ch_coeffs);
                for (int i = 0; i < dim_ch; i++) {
                    atomicAdd(v_ch_coeffs_ptr + 3*dim_ch*g + 3*i + 0, v_ch_coeff_local[i].x);
                    atomicAdd(v_ch_coeffs_ptr + 3*dim_ch*g + 3*i + 1, v_ch_coeff_local[i].y);
                    atomicAdd(v_ch_coeffs_ptr + 3*dim_ch*g + 3*i + 2, v_ch_coeff_local[i].z);
                }
                // atomicAdd(v_ch_coeffs_abs + g, v_ch_coeff_abs_local);

                atomicAdd(v_opacities + g, v_opacity_local);
            }
        }
    }

    if (inside) {
        v_depth_ref_im[pix_id] = v_depth_ref;

        // background gradient
        #if 0
        float3 v_bkg = {
            v_out.x * T_final,
            v_out.y * T_final,
            v_out.z * T_final
        };
        atomicAdd((float*)v_background+0, v_bkg.x);
        atomicAdd((float*)v_background+1, v_bkg.y);
        atomicAdd((float*)v_background+2, v_bkg.z);
        #endif
    }

}



template<CameraType CAMERA_TYPE>
__global__ void rasterize_simplified_forward_kernel(
    _ARGS_rasterize_simplified_forward_kernel
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

    // camera distortion
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y))
            done = true;
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 color_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    int tr = block.thread_rank();
    float T = 1.f;  // current/total visibility
    float3 normal_out = {0.f, 0.f, 0.f};  // sum of normals
    float3 pix_out = {0.f, 0.f, 0.f};  // output radiance
    float vis_sum = 0.f;  // output alpha
    float depth_sum = 0.f;  // output depth
    float depth_squared_sum = 0.f;  // for L2 depth regularizer
    float reg_depth_p = 0.f;  // output depth regularizer
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
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = opac;
        }
        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        int batch_size = min(block_size, range.y - batch_start);
        for (int t = 0; (t < batch_size) && !done; ++t) {
            glm::vec3 pos = position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            // if (!get_intersection(pos, axis_uv, pos_2d, poi, uv));
            //     continue;
            get_intersection(pos, axis_uv, pos_2d, poi, uv);
            if (glm::length(uv) > visibility_kernel_radius())
                continue;
            float alpha;
            if (!get_alpha(uv, opac, alpha))
                continue;

            const float next_T = T * (1.f - alpha);
            const float vis = alpha * T;

            // color
            glm::vec3 color = color_batch[t];
            pix_out.x = pix_out.x + color.x * vis;
            pix_out.y = pix_out.y + color.y * vis;
            pix_out.z = pix_out.z + color.z * vis;

            // depth regularization
            const float depth_raw = poi.z;
            const float depth = depth_map(depth_raw);
            {
                // float pairwise_l1 = vis*depth * vis_sum - vis * depth_sum;  // requires pos.z for depth
                float pairwise_l2 = vis * (vis_sum*depth*depth + depth_squared_sum - 2.0f*depth*depth_sum);
                reg_depth_p += pairwise_l2;
            }
            vis_sum += vis;
            depth_sum += vis*depth;
            depth_squared_sum += vis*depth*depth;

            // normal regularization
            glm::vec3 normal = get_normal_from_axisuv(axis_uv, poi);
            normal_out.x = normal_out.x + normal.x * vis;
            normal_out.y = normal_out.y + normal.y * vis;
            normal_out.z = normal_out.z + normal.z * vis;

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
        out_img[pix_id] = pix_out;
        out_depth[pix_id] = { depth_sum, depth_squared_sum };
        out_normal[pix_id] = normal_out;
        out_depth_reg[pix_id] = reg_depth_p;
    }
}


template<CameraType CAMERA_TYPE>
__global__ void rasterize_simplified_backward_kernel(
    _ARGS_rasterize_simplified_backward_kernel
) {
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
    // clamp this value to the last pixel
    const int32_t pix_id = min(i * img_size.x + j, img_size.x * img_size.y - 1);

    // keep not rasterizing threads around for reading data
    bool inside = (i < img_size.y && j < img_size.x);

    // camera distortion
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y))
            inside = false;
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    // have all threads in tile process the same gaussians in batches
    // first collect gaussians between range.x and range.y in batches
    // which gaussians to look through in this tile
    const int2 range = tile_bins[tile_id];
    const int block_size = block.size();
    const int num_batches = (range.y - range.x + block_size - 1) / block_size;

    __shared__ int32_t id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 color_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // df/d_out for this pixel
    const float3 v_out = nan_to_num(v_output_img[pix_id]);
    const float2 v_out_depth = nan_to_num(v_output_depth[pix_id]);
    const float3 v_out_normal = nan_to_num(v_output_normal[pix_id]);
    const float v_out_alpha = nan_to_num(v_output_alpha[pix_id]);
    const float v_reg_depth_p = nan_to_num(v_output_depth_reg[pix_id]);
    const float v_depth_sum = v_out_depth.x;
    const float v_depth_squared_sum = v_out_depth.y;

    // this is the T AFTER the last gaussian in this pixel
    float T_final = 1.0f - output_alpha[pix_id];
    float T = T_final;
    // index of last gaussian to contribute to this pixel
    const int bin_final = inside? final_index[pix_id] : 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing
    const int tr = block.thread_rank();

    const float vis_sum_final = 1.0f - T_final;
    const float depth_sum_final = output_depth[pix_id].x;
    const float depth_squared_sum_final = output_depth[pix_id].y;
    float vis_sum = vis_sum_final;

    float3 buffer = {0.f, 0.f, 0.f};
    float2 buffer_depth = {0.f, 0.f};  // depth, depth^2
    float3 buffer_normal = {0.f, 0.f, 0.f};
    float buffer_depth_reg = 0.f;

    // gradient
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(block);
    const int warp_bin_final = cg::reduce(warp, bin_final, cg::greater<int>());
    T = T_final;
    for (int b = 0; b < num_batches; ++b) {
        // resync all threads before writing next batch of shared mem
        block.sync();

        // each thread fetch 1 gaussian from back to front
        // 0 index will be furthest back in batch
        // index of gaussian to load
        // batch end is the index of the last gaussian in the batch
        const int batch_end = range.y - 1 - block_size * b;
        int batch_size = min(block_size, batch_end + 1 - range.x);
        const int idx = batch_end - tr;
        if (idx >= range.x) {
            int32_t g_id = gaussian_ids_sorted[idx];
            id_batch[tr] = g_id;
            const float3 pos = positions[g_id];
            const float opac = opacities[g_id];
            const float3 color = colors[g_id];
            const float3 v0 = axes_u[g_id];
            const float3 v1 = axes_v[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_batch[tr] = {color.x, color.y, color.z};
            opacity_batch[tr] = opac;
        }

        // wait for other threads to collect the gaussians in batch
        block.sync();

        // process gaussians in the current batch for this pixel
        // 0 index is the furthest back gaussian in the batch
        for (int t = max(0,batch_end - warp_bin_final); t < batch_size; ++t) {
            int valid = inside;
            if (batch_end - t > bin_final) {
                valid = 0;
            }

            glm::vec3 pos = position_batch[t];
            float opac = opacity_batch[t];
            glm::mat2x3 axis_uv = axes_uv_batch[t];

            glm::vec3 poi;
            glm::vec2 uv;
            if(valid){
                get_intersection(pos, axis_uv, pos_2d, poi, uv);
                if (glm::length(uv) > visibility_kernel_radius())
                    valid = 0;
            }
            float alpha;
            if (valid){
                if (!get_alpha(uv, opac, alpha))
                    valid = 0;
            }
            if(!warp.any(valid)){
                continue;
            }

            glm::vec3 v_position_local = {0.f, 0.f, 0.f};
            glm::vec2 v_position_xy_abs_local = {0.f, 0.f};
            glm::vec3 v_axis_u_local = {0.f, 0.f, 0.f};
            glm::vec3 v_axis_v_local = {0.f, 0.f, 0.f};
            glm::vec3 v_color_local = {0.f, 0.f, 0.f};
            float v_opacity_local = 0.f;
            //initialize everything to 0, only set if the lane is valid
            if(valid){
                // compute the current T for this gaussian
                const float ra = 1.f / (1.f - alpha);
                const float next_T = T * ra;
                const float vis = alpha * next_T;

                // update accumulation
                const float depth_raw = poi.z;
                const float depth = depth_map(depth_raw);
                float v_depth = 0.0f;
                v_depth += vis * v_depth_sum;
                v_depth += vis * 2.0f*depth * v_depth_squared_sum;

                // update depth regularizer
                float vis_sum_next = vis_sum - vis;
                // pairwise L2
                v_depth += v_reg_depth_p * vis * 2.0f * (
                    vis_sum_final * depth - depth_sum_final);
                float reg_depth_i =
                    vis_sum_final*depth*depth + depth_squared_sum_final
                    - 2.0f*depth*depth_sum_final;

                float v_depth_raw = depth_map_vjp(depth_raw, v_depth);
                glm::vec3 v_poi = {0.f, 0.f, v_depth_raw};

                // update color
                const float opacity = opacity_batch[t];
                glm::vec3 color = color_batch[t];
                v_color_local = {vis * v_out.x, vis * v_out.y, vis * v_out.z};

                // normal regularization
                glm::vec3 v_normal = {vis * v_out_normal.x, vis * v_out_normal.y, vis * v_out_normal.z};
                glm::mat2x3 v_axis_uv; glm::vec3 normal;
                get_normal_from_axisuv_vjp(axis_uv, poi, v_normal, normal, v_axis_uv);

                float v_alpha = 0.0f;
                // contribution from this pixel
                v_alpha += (color.x * T - buffer.x) * ra * v_out.x;
                v_alpha += (color.y * T - buffer.y) * ra * v_out.y;
                v_alpha += (color.z * T - buffer.z) * ra * v_out.z;
                v_alpha += T_final * ra * v_out_alpha;
                // v_alpha += -T_final * ra * background.x * v_out.x;
                // v_alpha += -T_final * ra * background.y * v_out.y;
                // v_alpha += -T_final * ra * background.z * v_out.z;
                // float v_alpha_color_only = v_alpha;
                v_alpha += (depth * T - buffer_depth.x) * ra * v_depth_sum;
                v_alpha += (depth*depth * T - buffer_depth.y) * ra * v_depth_squared_sum;
                v_alpha += (reg_depth_i * T - buffer_depth_reg) * ra * v_reg_depth_p;
                v_alpha += (normal.x * T - buffer_normal.x) * ra * v_out_normal.x;
                v_alpha += (normal.y * T - buffer_normal.y) * ra * v_out_normal.y;
                v_alpha += (normal.z * T - buffer_normal.z) * ra * v_out_normal.z;

                // update the running sum
                buffer.x += color.x * vis;
                buffer.y += color.y * vis;
                buffer.z += color.z * vis;
                buffer_depth.x += depth * vis;
                buffer_depth.y += depth*depth * vis;
                buffer_depth_reg += reg_depth_i * vis;
                buffer_normal.x += normal.x * vis;
                buffer_normal.y += normal.y * vis;
                buffer_normal.z += normal.z * vis;

                // grad
                glm::vec2 v_uv;
                get_alpha_vjp(
                    uv, opacity,
                    v_alpha, v_uv, v_opacity_local
                );
                glm::vec3 v_position_local_temp;
                get_intersection_vjp(
                    pos, axis_uv, pos_2d,
                    v_poi, v_uv,
                    v_position_local_temp, v_axis_uv
                );
                v_position_local += v_position_local_temp;
                v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local)) / glm::vec2(fx, fy);
                v_axis_u_local = v_axis_uv[0];
                v_axis_v_local = v_axis_uv[1];

                // next loop
                T = next_T;
                vis_sum = vis_sum_next;
            }
            warpSum3(v_position_local, warp);
            warpSum2(v_position_xy_abs_local, warp);
            warpSum3(v_axis_u_local, warp);
            warpSum3(v_axis_v_local, warp);
            warpSum3(v_color_local, warp);
            warpSum(v_opacity_local, warp);
            if (warp.thread_rank() == 0) {
                int32_t g = id_batch[t];

                float* v_position_ptr = (float*)(v_positions);
                atomicAdd(v_position_ptr + 3*g + 0, v_position_local.x);
                atomicAdd(v_position_ptr + 3*g + 1, v_position_local.y);
                atomicAdd(v_position_ptr + 3*g + 2, v_position_local.z);
                float* v_positions_xy_abs_ptr = (float*)(v_positions_xy_abs);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 0, v_position_xy_abs_local.x);
                atomicAdd(v_positions_xy_abs_ptr + 2*g + 1, v_position_xy_abs_local.y);

                float* v_axis_u_ptr = (float*)(v_axes_u);
                atomicAdd(v_axis_u_ptr + 3*g + 0, v_axis_u_local.x);
                atomicAdd(v_axis_u_ptr + 3*g + 1, v_axis_u_local.y);
                atomicAdd(v_axis_u_ptr + 3*g + 2, v_axis_u_local.z);
                float* v_axis_v_ptr = (float*)(v_axes_v);
                atomicAdd(v_axis_v_ptr + 3*g + 0, v_axis_v_local.x);
                atomicAdd(v_axis_v_ptr + 3*g + 1, v_axis_v_local.y);
                atomicAdd(v_axis_v_ptr + 3*g + 2, v_axis_v_local.z);
                
                float* v_color_ptr = (float*)(v_colors);
                atomicAdd(v_color_ptr + 3*g + 0, v_color_local.x);
                atomicAdd(v_color_ptr + 3*g + 1, v_color_local.y);
                atomicAdd(v_color_ptr + 3*g + 2, v_color_local.z);

                atomicAdd(v_opacities + g, v_opacity_local);
            }
        }
    }

}




template<CameraType CAMERA_TYPE>
__global__ void render_background_sh_forward_kernel(
    _ARGS_render_background_sh_forward_kernel
) {
    unsigned i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    int32_t pix_id = i * img_size.x + j;

    if (i >= img_size.y || j >= img_size.x) return;

    float fx = intrins.x, fy = intrins.y;
    float cx = intrins.z, cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y)) {
            out_img[pix_id] = {0.0f, 0.0f, 0.0f};
            return;
        }
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    float xi = pos_2d.x;
    float yi = -pos_2d.y;
    float zi = -1.0f;
    float xr = rotation[0] * xi + rotation[1] * yi + rotation[2] * zi;
    float yr = rotation[3] * xi + rotation[4] * yi + rotation[5] * zi;
    float zr = rotation[6] * xi + rotation[7] * yi + rotation[8] * zi;
    float norm = sqrtf(fmaxf(xr * xr + yr * yr + zr * zr, 1e-12f));
    float x = isfinite(xr) ? xr / norm : 0.0f;
    float y = isfinite(yr) ? yr / norm : 0.0f;
    float z = isfinite(zr) ? zr / norm : 0.0f;

    float xx = x*x, yy = y*y, zz = z*z;

    glm::vec3 color = glm::vec3(0.0f);
    glm::vec3 *sh_coeffs = (glm::vec3*)sh_coeffs_float3;

    // l0
    color += 0.28209479177387814f * sh_coeffs[0];

    // l1
    if (sh_degree > 1) {
        color += 0.4886025119029199f * y * sh_coeffs[1];
        color += 0.4886025119029199f * z * sh_coeffs[2];
        color += 0.4886025119029199f * x * sh_coeffs[3];
    }

    // l2
    if (sh_degree > 2) {
        color += 1.0925484305920792f * x * y * sh_coeffs[4];
        color += 1.0925484305920792f * y * z * sh_coeffs[5];
        color += (0.9461746957575601f * zz - 0.31539156525251999f) * sh_coeffs[6];
        color += 1.0925484305920792f * x * z * sh_coeffs[7];
        color += 0.5462742152960396f * (xx - yy) * sh_coeffs[8];
    }

    // l3
    if (sh_degree > 3) {
        color += 0.5900435899266435f * y * (3.0f * xx - yy) * sh_coeffs[9];
        color += 2.890611442640554f * x * y * z * sh_coeffs[10];
        color += 0.4570457994644658f * y * (5.0f * zz - 1.0f) * sh_coeffs[11];
        color += 0.3731763325901154f * z * (5.0f * zz - 3.0f) * sh_coeffs[12];
        color += 0.4570457994644658f * x * (5.0f * zz - 1.0f) * sh_coeffs[13];
        color += 1.445305721320277f * z * (xx - yy) * sh_coeffs[14];
        color += 0.5900435899266435f * x * (xx - 3.0f * yy) * sh_coeffs[15];
    }

    // l4
    if (sh_degree > 4) {
        color += 2.5033429417967046f * x * y * (xx - yy) * sh_coeffs[16];
        color += 1.7701307697799304f * y * z * (3.0f * xx - yy) * sh_coeffs[17];
        color += 0.9461746957575601f * x * y * (7.0f * zz - 1.0f) * sh_coeffs[18];
        color += 0.6690465435572892f * y * z * (7.0f * zz - 3.0f) * sh_coeffs[19];
        color += 0.10578554691520431f * (35.0f * zz * zz - 30.0f * zz + 3.0f) * sh_coeffs[20];
        color += 0.6690465435572892f * x * z * (7.0f * zz - 3.0f) * sh_coeffs[21];
        color += 0.47308734787878004f * (xx - yy) * (7.0f * zz - 1.0f) * sh_coeffs[22];
        color += 1.7701307697799304f * x * z * (xx - 3.0f * yy) * sh_coeffs[23];
        color += 0.6258357354491761f * (xx * (xx - 3.0f * yy) - yy * (3.0f * xx - yy)) * sh_coeffs[24];
    }

    color.x = fmaxf(color.x + 0.5f, 0.0f);
    color.y = fmaxf(color.y + 0.5f, 0.0f);
    color.z = fmaxf(color.z + 0.5f, 0.0f);

    out_img[pix_id] = *(float3*)&color;
}


template<CameraType CAMERA_TYPE>
__global__ void render_background_sh_backward_kernel(
    _ARGS_render_background_sh_backward_kernel
) {
    unsigned i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    int32_t pix_id = i * img_size.x + j;

    bool inside = (i < img_size.y && j < img_size.x);

    unsigned idx = i * img_size.x + j;
    glm::vec3 v_color = glm::vec3(0.0);
    if (inside) {
        glm::vec3 color = ((glm::vec3*)out_color)[idx];
        v_color = ((glm::vec3*)v_out_color)[idx];
        if (color.x == 0.0f || !isfinite(v_color.x)) v_color.x = 0.0f;
        if (color.y == 0.0f || !isfinite(v_color.y)) v_color.y = 0.0f;
        if (color.z == 0.0f || !isfinite(v_color.z)) v_color.z = 0.0f;
        // v_color = glm::clamp(v_color, -glm::vec3(1e4f), glm::vec3(1e4f));
    }

    float fx = intrins.x, fy = intrins.y;
    float cx = intrins.z, cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted && inside) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y))
            inside = false;
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }
    if (__syncthreads_count(inside) == 0)
        return;

    auto block = cg::this_thread_block();
    cg::thread_block_tile<32> warp = cg::tiled_partition<32>(block);

    float xi = pos_2d.x;
    float yi = -pos_2d.y;
    float zi = -1.0f;
    float xr = rotation[0] * xi + rotation[1] * yi + rotation[2] * zi;
    float yr = rotation[3] * xi + rotation[4] * yi + rotation[5] * zi;
    float zr = rotation[6] * xi + rotation[7] * yi + rotation[8] * zi;
    float norm2 = xr * xr + yr * yr + zr * zr;
    float norm = sqrtf(fmaxf(norm2, 1e-12f));
    float x = inside && isfinite(xr) ? xr / norm : 0.0f;
    float y = inside && isfinite(yr) ? yr / norm : 0.0f;
    float z = inside && isfinite(zr) ? zr / norm : 0.0f;

    float xx = x*x, yy = y*y, zz = z*z;

    float v_x = 0.0f, v_y = 0.0f, v_z = 0.0f;
    float v_xx = 0.0f, v_yy = 0.0f, v_zz = 0.0f;

    glm::vec3 *sh_coeffs = (glm::vec3*)sh_coeffs_float3;

    glm::vec3 v_sh;
    #define _ATOMIC_ADD_SH_COEFFS(idx) \
        warpSum3(v_sh, warp); \
        if (warp.thread_rank() == idx) { \
            atomicAdd(&v_sh_coeffs[idx].x, v_sh.x); \
            atomicAdd(&v_sh_coeffs[idx].y, v_sh.y); \
            atomicAdd(&v_sh_coeffs[idx].z, v_sh.z); \
        }
        // __syncthreads();

    // l0
    float v_color_dot_sh_coeff = 0.0f;
    v_sh = 0.28209479177387814f * v_color;
    _ATOMIC_ADD_SH_COEFFS(0);

    // l1 - manually calculated
    if (sh_degree > 1) {

        // color += 0.4886025119029199f * y * sh_coeffs[1];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[1]);
        v_y += 0.4886025119029199f * v_color_dot_sh_coeff;
        v_sh = 0.4886025119029199f * y * v_color;
        _ATOMIC_ADD_SH_COEFFS(1);

        // color += 0.4886025119029199f * z * sh_coeffs[2];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[2]);
        v_z += 0.4886025119029199f * v_color_dot_sh_coeff;
        v_sh = 0.4886025119029199f * z * v_color;
        _ATOMIC_ADD_SH_COEFFS(2);

        // color += 0.4886025119029199f * x * sh_coeffs[3];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[3]);
        v_x += 0.4886025119029199f * v_color_dot_sh_coeff;
        v_sh = 0.4886025119029199f * x * v_color;
        _ATOMIC_ADD_SH_COEFFS(3);
    }

    // l2 - manually calculated
    if (sh_degree > 2) {

        // color += 1.0925484305920792f * x * y * sh_coeffs[4];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[4]);
        v_x += 1.0925484305920792f * y * v_color_dot_sh_coeff;
        v_y += 1.0925484305920792f * x * v_color_dot_sh_coeff;
        v_sh = 1.0925484305920792f * x * y * v_color;
        _ATOMIC_ADD_SH_COEFFS(4);

        // color += 1.0925484305920792f * y * z * sh_coeffs[5];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[5]);
        v_z += 1.0925484305920792f * y * v_color_dot_sh_coeff;
        v_y += 1.0925484305920792f * z * v_color_dot_sh_coeff;
        v_sh = 1.0925484305920792f * y * z * v_color;
        _ATOMIC_ADD_SH_COEFFS(5);

        // color += (0.9461746957575601f * zz - 0.31539156525251999f) * sh_coeffs[6];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[6]);
        v_zz += 0.9461746957575601f * v_color_dot_sh_coeff;
        v_sh = (0.9461746957575601f * zz - 0.31539156525251999f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(6);

        // color += 1.0925484305920792f * x * z * sh_coeffs[7];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[7]);
        v_x += 1.0925484305920792f * z * v_color_dot_sh_coeff;
        v_z += 1.0925484305920792f * x * v_color_dot_sh_coeff;
        v_sh = 1.0925484305920792f * x * z * v_color;
        _ATOMIC_ADD_SH_COEFFS(7);

        // color += 0.5462742152960396f * (xx - yy) * sh_coeffs[8];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[8]);
        v_xx += 0.5462742152960396f * v_color_dot_sh_coeff;
        v_yy -= 0.5462742152960396f * v_color_dot_sh_coeff;
        v_sh = 0.5462742152960396f * (xx - yy) * v_color;
        _ATOMIC_ADD_SH_COEFFS(8);
    }

    // l3 - AI generated, one incorrect line commented
    if (sh_degree > 3) {
        // color += 0.5900435899266435f * y * (3.0f * xx - yy) * sh_coeffs[9];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[9]);
        v_xx += 1.7701307697799305f * y * v_color_dot_sh_coeff;
        v_yy -= 0.5900435899266435f * y * v_color_dot_sh_coeff;
        v_y += 0.5900435899266435f * (3.0f * xx - yy) * v_color_dot_sh_coeff;
        v_sh = 0.5900435899266435f * y * (3.0f * xx - yy) * v_color;
        _ATOMIC_ADD_SH_COEFFS(9);

        // color += 2.890611442640554f * x * y * z * sh_coeffs[10];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[10]);
        v_x += 2.890611442640554f * y * z * v_color_dot_sh_coeff;
        v_y += 2.890611442640554f * x * z * v_color_dot_sh_coeff;
        v_z += 2.890611442640554f * x * y * v_color_dot_sh_coeff;
        v_sh = 2.890611442640554f * x * y * z * v_color;
        _ATOMIC_ADD_SH_COEFFS(10);

        // color += 0.4570457994644658f * y * (5.0f * zz - 1.0f) * sh_coeffs[11];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[11]);
        v_zz += 2.285228997322329f * y * v_color_dot_sh_coeff;
        v_y += 0.4570457994644658f * (5.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_sh = 0.4570457994644658f * y * (5.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(11);

        // color += 0.3731763325901154f * z * (5.0f * zz - 3.0f) * sh_coeffs[12];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[12]);
        v_z += 0.3731763325901154f * (5.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_zz += 1.865881662950577f * z * v_color_dot_sh_coeff;
        v_sh = 0.3731763325901154f * z * (5.0f * zz - 3.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(12);

        // color += 0.4570457994644658f * x * (5.0f * zz - 1.0f) * sh_coeffs[13];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[13]);
        v_x += 0.4570457994644658f * (5.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_zz += 2.285228997322329f * x * v_color_dot_sh_coeff;
        v_sh = 0.4570457994644658f * x * (5.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(13);

        // color += 1.445305721320277f * z * (xx - yy) * sh_coeffs[14];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[14]);
        v_xx += 1.445305721320277f * z * v_color_dot_sh_coeff;
        v_yy -= 1.445305721320277f * z * v_color_dot_sh_coeff;
        v_z += 1.445305721320277f * (xx - yy) * v_color_dot_sh_coeff;
        v_sh = 1.445305721320277f * z * (xx - yy) * v_color;
        _ATOMIC_ADD_SH_COEFFS(14);

        // color += 0.5900435899266435f * x * (xx - 3.0f * yy) * sh_coeffs[15];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[15]);
        // v_xx += 1.1800871798532870f * x * v_color_dot_sh_coeff;
        v_xx += 0.5900435899266435f * x * v_color_dot_sh_coeff;
        v_yy -= 1.7701307697799305f * x * v_color_dot_sh_coeff;
        v_x += 0.5900435899266435f * (xx - 3.0f * yy) * v_color_dot_sh_coeff;
        v_sh = 0.5900435899266435f * x * (xx - 3.0f * yy) * v_color;
        _ATOMIC_ADD_SH_COEFFS(15);
    }

    // l4 - AI generated, two incorrect lines commented
    if (sh_degree > 4) {
        // color += 2.5033429417967046f * x * y * (xx - yy) * sh_coeffs[16];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[16]);
        v_x += 2.5033429417967046f * y * (xx - yy) * v_color_dot_sh_coeff;
        v_y += 2.5033429417967046f * x * (xx - yy) * v_color_dot_sh_coeff;
        v_xx += 2.5033429417967046f * x * y * v_color_dot_sh_coeff;
        v_yy -= 2.5033429417967046f * x * y * v_color_dot_sh_coeff;
        v_sh = 2.5033429417967046f * x * y * (xx - yy) * v_color;
        _ATOMIC_ADD_SH_COEFFS(16);

        // color += 1.7701307697799304f * y * z * (3.0f * xx - yy) * sh_coeffs[17];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[17]);
        v_xx += 5.3103923093397912f * y * z * v_color_dot_sh_coeff;
        v_yy -= 1.7701307697799304f * y * z * v_color_dot_sh_coeff;
        v_y += 1.7701307697799304f * z * (3.0f * xx - yy) * v_color_dot_sh_coeff;
        v_z += 1.7701307697799304f * y * (3.0f * xx - yy) * v_color_dot_sh_coeff;
        v_sh = 1.7701307697799304f * y * z * (3.0f * xx - yy) * v_color;
        _ATOMIC_ADD_SH_COEFFS(17);

        // color += 0.9461746957575601f * x * y * (7.0f * zz - 1.0f) * sh_coeffs[18];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[18]);
        v_x += 0.9461746957575601f * y * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_y += 0.9461746957575601f * x * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_zz += 6.6232228703029207f * x * y * v_color_dot_sh_coeff;
        v_sh = 0.9461746957575601f * x * y * (7.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(18);

        // color += 0.6690465435572892f * y * z * (7.0f * zz - 3.0f) * sh_coeffs[19];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[19]);
        v_y += 0.6690465435572892f * z * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_z += 0.6690465435572892f * y * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_zz += 4.6833258049010244f * y * z * v_color_dot_sh_coeff;
        v_sh = 0.6690465435572892f * y * z * (7.0f * zz - 3.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(19);

        // color += 0.10578554691520431f * (35.0f * zz * zz - 30.0f * zz + 3.0f) * sh_coeffs[20];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[20]);
        v_zz += 0.10578554691520431f * (70.0f * zz - 30.0f) * v_color_dot_sh_coeff;
        v_sh = 0.10578554691520431f * (35.0f * zz * zz - 30.0f * zz + 3.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(20);

        // color += 0.6690465435572892f * x * z * (7.0f * zz - 3.0f) * sh_coeffs[21];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[21]);
        v_x += 0.6690465435572892f * z * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_z += 0.6690465435572892f * x * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_zz += 4.6833258049010244f * x * z * v_color_dot_sh_coeff;
        v_sh = 0.6690465435572892f * x * z * (7.0f * zz - 3.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(21);

        // color += 0.47308734787878004f * (xx - yy) * (7.0f * zz - 1.0f) * sh_coeffs[22];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[22]);
        v_xx += 0.47308734787878004f * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_yy -= 0.47308734787878004f * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_zz += 3.3116114351514603f * (xx - yy) * v_color_dot_sh_coeff;
        v_sh = 0.47308734787878004f * (xx - yy) * (7.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD_SH_COEFFS(22);

        // color += 1.7701307697799304f * x * z * (xx - 3.0f * yy) * sh_coeffs[23];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[23]);
        v_x += 1.7701307697799304f * z * (xx - 3.0f * yy) * v_color_dot_sh_coeff;
        v_z += 1.7701307697799304f * x * (xx - 3.0f * yy) * v_color_dot_sh_coeff;
        v_xx += 1.7701307697799304f * x * z * v_color_dot_sh_coeff;
        v_yy -= 5.3103923093397912f * x * z * v_color_dot_sh_coeff;
        v_sh = 1.7701307697799304f * x * z * (xx - 3.0f * yy) * v_color;
        _ATOMIC_ADD_SH_COEFFS(23);

        // color += 0.6258357354491761f * (xx * (xx - 3.0f * yy) - yy * (3.0f * xx - yy)) * sh_coeffs[24];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[24]);
        // v_xx += 0.6258357354491761f * (4.0f * xx - 6.0f * yy) * v_color_dot_sh_coeff;
        // v_yy += 0.6258357354491761f * (6.0f * yy - 12.0f * xx) * v_color_dot_sh_coeff;
        v_xx += 0.6258357354491761f * (2.0f * xx - 6.0f * yy) * v_color_dot_sh_coeff;
        v_yy += 0.6258357354491761f * (2.0f * yy - 6.0f * xx) * v_color_dot_sh_coeff;
        v_sh = 0.6258357354491761f * (xx * (xx - 3.0f * yy) - yy * (3.0f * xx - yy)) * v_color;
        _ATOMIC_ADD_SH_COEFFS(24);
    }

    #undef _ATOMIC_ADD_SH_COEFFS

    v_x += v_xx * 2.0f*x;
    v_y += v_yy * 2.0f*y;
    v_z += v_zz * 2.0f*z;

    glm::vec3 xyz = glm::vec3(x, y, z);
    glm::mat3 dp_dpr = (glm::mat3(1.0f) - glm::outerProduct(xyz, xyz)) / norm;
    glm::vec3 v_p = dp_dpr * glm::vec3(v_x, v_y, v_z);
    v_p *= (inside ? 1.0f : 0.0f);

    float tmp[9] = {
        v_p.x * xi, v_p.x * yi, v_p.x * zi,
        v_p.y * xi, v_p.y * yi, v_p.y * zi,
        v_p.z * xi, v_p.z * yi, v_p.z * zi
    };
    #pragma unroll
    for (int i = 0; i < 9; i++) {
        warpSum(tmp[i], warp);
        if (warp.thread_rank() == i)
            atomicAdd(&v_rotation[i], tmp[i]);
    }
}



template __global__ void rasterize_simple_forward_kernel<CameraType::Undistorted>(
    _ARGS_rasterize_simple_forward_kernel
);
template __global__ void rasterize_simple_forward_kernel<CameraType::GenericDistorted>(
    _ARGS_rasterize_simple_forward_kernel
);


template __global__ void rasterize_depth_forward_kernel
<DepthMode::Mean, CameraType::Undistorted>(
    _ARGS_rasterize_depth_forward_kernel
);
template __global__ void rasterize_depth_forward_kernel
<DepthMode::Mean, CameraType::GenericDistorted>(
    _ARGS_rasterize_depth_forward_kernel
);
template __global__ void rasterize_depth_forward_kernel
<DepthMode::Median, CameraType::Undistorted>(
    _ARGS_rasterize_depth_forward_kernel
);
template __global__ void rasterize_depth_forward_kernel
<DepthMode::Median, CameraType::GenericDistorted>(
    _ARGS_rasterize_depth_forward_kernel
);


template __global__ void rasterize_depth_backward_kernel<DepthMode::Mean>(
    _ARGS_rasterize_depth_backward_kernel
);
template __global__ void rasterize_depth_backward_kernel<DepthMode::Median>(
    _ARGS_rasterize_depth_backward_kernel
);


template __global__ void rasterize_simplified_forward_kernel<CameraType::Undistorted>(
    _ARGS_rasterize_simplified_forward_kernel
);
template __global__ void rasterize_simplified_forward_kernel<CameraType::GenericDistorted>(
    _ARGS_rasterize_simplified_forward_kernel
);


template __global__ void rasterize_simplified_backward_kernel<CameraType::Undistorted>(
    _ARGS_rasterize_simplified_backward_kernel
);
template __global__ void rasterize_simplified_backward_kernel<CameraType::GenericDistorted>(
    _ARGS_rasterize_simplified_backward_kernel
);


template __global__ void render_background_sh_forward_kernel<CameraType::Undistorted>(
    _ARGS_render_background_sh_forward_kernel
);
template __global__ void render_background_sh_forward_kernel<CameraType::GenericDistorted>(
    _ARGS_render_background_sh_forward_kernel
);
template __global__ void render_background_sh_backward_kernel<CameraType::Undistorted>(
    _ARGS_render_background_sh_backward_kernel
);
template __global__ void render_background_sh_backward_kernel<CameraType::GenericDistorted>(
    _ARGS_render_background_sh_backward_kernel
);
