#include "backward.cuh"
#include "helpers.cuh"
#include <cuda_fp16.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
namespace cg = cooperative_groups;

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

__global__ void rasterize_backward_kernel(
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
    const float3& __restrict__ background,
    const float2* __restrict__ depth_grads,
    const float2* __restrict__ depth_normal_ref_im,
    const int* __restrict__ final_index,
    const float* __restrict__ output_alpha,
    const float3* __restrict__ output_depth_grad,
    const float* __restrict__ v_output_alpha,
    const float3* __restrict__ v_output,
    const float3* __restrict__ v_output_depth_grad,
    const float* __restrict__ v_output_reg_depth,
    const float* __restrict__ v_output_reg_normal,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float3* __restrict__ v_colors,
    float* __restrict__ v_opacities,
    float2* __restrict__ v_depth_grad,
    float2* __restrict__ v_depth_normal_ref
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

    __shared__ int32_t id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec4 color_opacity_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec2 depth_grad_batch[MAX_BLOCK_SIZE];

    // df/d_out for this pixel
    const float3 out_depth_grad = output_depth_grad[pix_id];
    const float3 v_out = v_output[pix_id];
    const float3 v_out_depth_grad = v_output_depth_grad[pix_id];
    const float v_out_alpha = v_output_alpha[pix_id];
    const float v_out_reg_depth = v_output_reg_depth[pix_id];
    const float v_out_reg_normal = v_output_reg_normal[pix_id];
    const glm::vec2 v_g_sum = {v_out_depth_grad.x, v_out_depth_grad.y};
    const float v_depth_sum = v_out_depth_grad.z;

    // this is the T AFTER the last gaussian in this pixel
    float T_final = 1.0f - output_alpha[pix_id];
    float T = T_final;
    // index of last gaussian to contribute to this pixel
    const int bin_final = inside? final_index[pix_id] : 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing
    const int tr = block.thread_rank();

    // regularization
    const float2 depth_normal_ref = inside ?
        depth_normal_ref_im[pix_id] : make_float2(0.f, 0.f);
    glm::vec2 n_bar = {depth_normal_ref.x, depth_normal_ref.y};
    glm::vec2 v_n_bar = {0.f, 0.f};

    const float vis_sum_final = 1.0f - T_final;
    const float depth_sum_final = out_depth_grad.z;
    float vis_sum = vis_sum_final;
    float depth_sum = depth_sum_final;
    glm::vec2 g_sum = {out_depth_grad.x, out_depth_grad.y};

    float3 buffer = {0.f, 0.f, 0.f};
    float3 buffer_depth = {0.f, 0.f, 0.f};
    float buffer_depth_reg = 0.f;
    float buffer_normal_reg = 0.f;
    
    float v_sum_vis = v_out_alpha;

    glm::vec2 v_g_bar = {v_out_depth_grad.x, v_out_depth_grad.y};
    float v_depth_out = v_out_depth_grad.z;

    // second run through, full gradient calculation
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
            const float2 depth_grad = depth_grads[g_id];
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
            color_opacity_batch[tr] = {color.x, color.y, color.z, opac};
            depth_grad_batch[tr] = {depth_grad.x, depth_grad.y};
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
            float opac = color_opacity_batch[t].w;
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
            glm::vec2 v_depth_grad_local = {0.f, 0.f};
            float v_opacity_local = 0.f;
            //initialize everything to 0, only set if the lane is valid
            if(valid){
                // compute the current T for this gaussian
                const float ra = 1.f / (1.f - alpha);
                const float next_T = T * ra;
                const float vis = alpha * next_T;

                // update depth regularizer
                const glm::vec2 depth_grad = depth_grad_batch[t];
                const float depth = pos.z;
                float vis_sum_next = vis_sum - vis;
                float depth_sum_next = depth_sum - vis*depth;
                v_position_local.z += v_out_reg_depth * vis * (vis_sum_next - (vis_sum_final-vis_sum));
                float reg_depth_i = (
                    depth * vis_sum_next - depth_sum_next +
                    (depth_sum_final-depth_sum) - depth * (vis_sum_final-vis_sum)
                );

                // update normal regularizer
                glm::vec2 g_i = {depth_grad.x, depth_grad.y};
                float g_i_norm = glm::length(g_i) + 1e-6f;
                glm::vec2 n_i = g_i / g_i_norm;
                glm::mat2 J_i = (glm::mat2(1.0f) - glm::outerProduct(n_i, n_i)) / g_i_norm;
                float reg_normal_i = 1.0f - dot(n_i, n_bar);
                glm::vec2 v_normal_glm = v_out_reg_normal * (-vis) * J_i * n_bar;
                v_depth_grad_local.x += v_normal_glm.x;
                v_depth_grad_local.y += v_normal_glm.y;
                v_n_bar += vis * (-n_i) * v_out_reg_normal;

                // update v_rgb for this gaussian
                v_color_local = {vis * v_out.x, vis * v_out.y, vis * v_out.z};
                v_depth_grad_local.x += vis * v_g_sum.x;
                v_depth_grad_local.y += vis * v_g_sum.y;
                v_position_local.z += vis * v_depth_sum;

                float v_alpha = 0.0f;
                const glm::vec4 rgba = color_opacity_batch[t];
                // contribution from this pixel
                v_alpha += (rgba.x * T - buffer.x) * ra * v_out.x;
                v_alpha += (rgba.y * T - buffer.y) * ra * v_out.y;
                v_alpha += (rgba.z * T - buffer.z) * ra * v_out.z;
                v_alpha += (depth_grad.x * T - buffer_depth.x) * ra * v_g_sum.x;
                v_alpha += (depth_grad.y * T - buffer_depth.y) * ra * v_g_sum.y;
                v_alpha += (pos.z * T - buffer_depth.z) * ra * v_depth_sum;
                v_alpha += (reg_depth_i * T - buffer_depth_reg) * ra * v_out_reg_depth;
                v_alpha += (reg_normal_i * T - buffer_normal_reg) * ra * v_out_reg_normal;

                v_alpha += T_final * ra * v_out_alpha;
                // contribution from background pixel
                v_alpha += -T_final * ra * background.x * v_out.x;
                v_alpha += -T_final * ra * background.y * v_out.y;
                v_alpha += -T_final * ra * background.z * v_out.z;
                // update the running sum
                buffer.x += rgba.x * vis;
                buffer.y += rgba.y * vis;
                buffer.z += rgba.z * vis;
                buffer_depth.x += depth_grad.x * vis;
                buffer_depth.y += depth_grad.y * vis;
                buffer_depth.z += pos.z * vis;
                buffer_depth_reg += reg_depth_i * vis;
                buffer_normal_reg += reg_normal_i * vis;

                glm::vec2 v_uv;
                get_alpha_vjp(
                    uv, rgba.w, v_alpha,
                    v_uv, v_opacity_local
                );
                glm::mat2x3 v_axis_uv;
                glm::vec3 v_position_local_temp;
                get_intersection_vjp(
                    pos, axis_uv, pos_2d,
                    glm::vec3(0), v_uv,
                    v_position_local_temp, v_axis_uv
                );
                v_position_local += v_position_local_temp;
                v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local));
                // v_position_xy_abs_local /= pos.z;
                v_axis_u_local = v_axis_uv[0];
                v_axis_v_local = v_axis_uv[1];

                T = next_T;
                vis_sum = vis_sum_next;
                depth_sum = depth_sum_next;
            }
            warpSum3(v_position_local, warp);
            warpSum2(v_position_xy_abs_local, warp);
            warpSum3(v_axis_u_local, warp);
            warpSum3(v_axis_v_local, warp);
            warpSum3(v_color_local, warp);
            warpSum(v_opacity_local, warp);
            warpSum2(v_depth_grad_local, warp);
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

                float* v_depth_grad_ptr = (float*)(v_depth_grad);
                atomicAdd(v_depth_grad_ptr + 2*g + 0, v_depth_grad_local.x);
                atomicAdd(v_depth_grad_ptr + 2*g + 1, v_depth_grad_local.y);
            }
        }
    }

    if (inside) {
        v_depth_normal_ref[pix_id] = {v_n_bar.x, v_n_bar.y};
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
    const float3& __restrict__ background,
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
    __shared__ glm::vec4 color_opacity_batch[MAX_BLOCK_SIZE];

    // df/d_out for this pixel
    const float3 v_out = v_output[pix_id];
    const float v_out_alpha = v_output_alpha[pix_id];

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
            color_opacity_batch[tr] = {color.x, color.y, color.z, opac};
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
            float opac = color_opacity_batch[t].w;
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

                const glm::vec4 rgba = color_opacity_batch[t];
                // contribution from this pixel
                v_alpha += (rgba.x * T - buffer.x) * ra * v_out.x;
                v_alpha += (rgba.y * T - buffer.y) * ra * v_out.y;
                v_alpha += (rgba.z * T - buffer.z) * ra * v_out.z;

                v_alpha += T_final * ra * v_out_alpha;
                // contribution from background pixel
                v_alpha += -T_final * ra * background.x * v_out.x;
                v_alpha += -T_final * ra * background.y * v_out.y;
                v_alpha += -T_final * ra * background.z * v_out.z;
                // update the running sum
                buffer.x += rgba.x * vis;
                buffer.y += rgba.y * vis;
                buffer.z += rgba.z * vis;

                glm::vec2 v_uv;
                get_alpha_vjp(
                    uv, rgba.w, v_alpha,
                    v_uv, v_opacity_local
                );
                glm::mat2x3 v_axis_uv;
                get_intersection_vjp(
                    pos, axis_uv, pos_2d,
                    glm::vec3(0), v_uv,
                    v_position_local, v_axis_uv
                );
                v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local));
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

__global__ void project_gaussians_backward_kernel(
    const int num_points,
    const float3* __restrict__ means3d,
    const float2* __restrict__ scales,
    const float4* __restrict__ quats,
    const float* __restrict__ viewmat,
    const float4 intrins,
    const int* __restrict__ num_tiles_hit,
    const float3* __restrict__ v_positions,
    const float3* __restrict__ v_axes_u,
    const float3* __restrict__ v_axes_v,
    const float2* __restrict__ v_depth_grads,
    float3* __restrict__ v_means3d,
    float2* __restrict__ v_scales,
    float4* __restrict__ v_quats
) {
    unsigned idx = cg::this_grid().thread_rank(); // idx of thread within grid
    if (idx >= num_points || num_tiles_hit <= 0) {
        return;
    }

    // position
    float3 p_world = means3d[idx];
    float fx = intrins.x;
    float fy = intrins.y;
    float3 p_view = transform_4x3(viewmat, p_world);
    float3 v_p_view = v_positions[idx];

    // forward
    float2 scale = scales[idx];
    float4 quat = quats[idx];
    glm::mat3 R1 = glm::transpose(glm::mat3(
        viewmat[0], viewmat[1], viewmat[2],
        viewmat[4], viewmat[5], viewmat[6],
        viewmat[8], viewmat[9], viewmat[10]
    ));
    glm::mat3 R2 = quat_to_rotmat(quat);
    glm::mat3 R = R1 * R2;
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
    glm::mat3 v_R_dg;
    float3 v_p_view_dg = {0.f, 0.f, 0.f};
    projected_depth_grad_vjp(
        p_view, R, fx, fy, v_depth_grads[idx],
        v_p_view_dg, v_R_dg);
    v_R += v_R_dg;
    v_p_view.x += v_p_view_dg.x;
    v_p_view.y += v_p_view_dg.y;
    v_p_view.z += v_p_view_dg.z;

    float3 v_p_world = transform_4x3_rot_only_transposed(viewmat, v_p_view);
    v_means3d[idx] = v_p_world;
    v_scales[idx] = v_scale;
    glm::mat3 v_R2 = glm::transpose(R1) * v_R;
    float4 v_quat = quat_to_rotmat_vjp(quat, v_R2);
    v_quats[idx] = v_quat;
}

// output space: 2D covariance, input space: cov3d
__device__ void project_cov3d_ewa_vjp(
    const float3& __restrict__ mean3d,
    const float* __restrict__ cov3d,
    const float* __restrict__ viewmat,
    const float fx,
    const float fy,
    const float3& __restrict__ v_cov2d,
    float3& __restrict__ v_mean3d,
    float* __restrict__ v_cov3d
) {
    // viewmat is row major, glm is column major
    // upper 3x3 submatrix
    // clang-format off
    glm::mat3 W = glm::mat3(
        viewmat[0], viewmat[4], viewmat[8],
        viewmat[1], viewmat[5], viewmat[9],
        viewmat[2], viewmat[6], viewmat[10]
    );
    // clang-format on
    glm::vec3 p = glm::vec3(viewmat[3], viewmat[7], viewmat[11]);
    glm::vec3 t = W * glm::vec3(mean3d.x, mean3d.y, mean3d.z) + p;
    float rz = 1.f / t.z;
    float rz2 = rz * rz;

    // column major
    // we only care about the top 2x2 submatrix
    // clang-format off
    glm::mat3 J = glm::mat3(
        fx * rz,         0.f,             0.f,
        0.f,             fy * rz,         0.f,
        -fx * t.x * rz2, -fy * t.y * rz2, 0.f
    );
    glm::mat3 V = glm::mat3(
        cov3d[0], cov3d[1], cov3d[2],
        cov3d[1], cov3d[3], cov3d[4],
        cov3d[2], cov3d[4], cov3d[5]
    );
    // cov = T * V * Tt; G = df/dcov = v_cov
    // -> d/dV = Tt * G * T
    // -> df/dT = G * T * Vt + Gt * T * V
    glm::mat3 v_cov = glm::mat3(
        v_cov2d.x,        0.5f * v_cov2d.y, 0.f,
        0.5f * v_cov2d.y, v_cov2d.z,        0.f,
        0.f,              0.f,              0.f
    );
    // clang-format on

    glm::mat3 T = J * W;
    glm::mat3 Tt = glm::transpose(T);
    glm::mat3 Vt = glm::transpose(V);
    glm::mat3 v_V = Tt * v_cov * T;
    glm::mat3 v_T = v_cov * T * Vt + glm::transpose(v_cov) * T * V;

    // vjp of cov3d parameters
    // v_cov3d_i = v_V : dV/d_cov3d_i
    // where : is frobenius inner product
    v_cov3d[0] = v_V[0][0];
    v_cov3d[1] = v_V[0][1] + v_V[1][0];
    v_cov3d[2] = v_V[0][2] + v_V[2][0];
    v_cov3d[3] = v_V[1][1];
    v_cov3d[4] = v_V[1][2] + v_V[2][1];
    v_cov3d[5] = v_V[2][2];

    // compute df/d_mean3d
    // T = J * W
    glm::mat3 v_J = v_T * glm::transpose(W);
    float rz3 = rz2 * rz;
    glm::vec3 v_t = glm::vec3(
        -fx * rz2 * v_J[2][0],
        -fy * rz2 * v_J[2][1],
        -fx * rz2 * v_J[0][0] + 2.f * fx * t.x * rz3 * v_J[2][0] -
            fy * rz2 * v_J[1][1] + 2.f * fy * t.y * rz3 * v_J[2][1]
    );
    // printf("v_t %.2f %.2f %.2f\n", v_t[0], v_t[1], v_t[2]);
    // printf("W %.2f %.2f %.2f\n", W[0][0], W[0][1], W[0][2]);
    v_mean3d.x += (float)glm::dot(v_t, W[0]);
    v_mean3d.y += (float)glm::dot(v_t, W[1]);
    v_mean3d.z += (float)glm::dot(v_t, W[2]);
}

// given cotangent v in output space (e.g. d_L/d_cov3d) in R(6)
// compute vJp for scale and rotation
__device__ void scale_rot_to_cov3d_vjp(
    const float2 scale,
    const float4 quat,
    const float* __restrict__ v_cov3d,
    float2& __restrict__ v_scale,
    float4& __restrict__ v_quat
) {
    // cov3d is upper triangular elements of matrix
    // off-diagonal elements count grads from both ij and ji elements,
    // must halve when expanding back into symmetric matrix
    glm::mat3 v_V = glm::mat3(
        v_cov3d[0],
        0.5 * v_cov3d[1],
        0.5 * v_cov3d[2],
        0.5 * v_cov3d[1],
        v_cov3d[3],
        0.5 * v_cov3d[4],
        0.5 * v_cov3d[2],
        0.5 * v_cov3d[4],
        v_cov3d[5]
    );
    glm::mat3 R = quat_to_rotmat(quat);
    glm::mat3 S = scale_to_mat({ scale.x, scale.y, 0.0f });
    glm::mat3 M = R * S;
    // https://math.stackexchange.com/a/3850121
    // for D = W * X, G = df/dD
    // df/dW = G * XT, df/dX = WT * G
    glm::mat3 v_M = 2.f * v_V * M;
    // glm::mat3 v_S = glm::transpose(R) * v_M;
    v_scale.x = (float)glm::dot(R[0], v_M[0]);
    v_scale.y = (float)glm::dot(R[1], v_M[1]);
    // v_scale.z = (float)glm::dot(R[2], v_M[2]);

    glm::mat3 v_R = v_M * S;
    v_quat = quat_to_rotmat_vjp(quat, v_R);
}


__device__ void projected_depth_grad_vjp(
    const float3 p, const glm::mat3 R,
    const float fx, const float fy,
    const float2 v_depth_grad,
    float3 &v_p_view, glm::mat3 &v_R
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
