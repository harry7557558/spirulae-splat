#include "helpers.cuh"
#include "rasterization_sorted.cuh"
#include <algorithm>

#include "stdio.h"


__global__ void rasterize_indices_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ gaussian_ids_sorted,
    const int2* __restrict__ tile_bins,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    int* __restrict__ num_intersects,
    int32_t* __restrict__ sorted_indices_,
    float* __restrict__ sorted_depths_
) {
    // each thread draws one pixel, but also timeshares caching splats in a shared tile

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

    __shared__ int id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    float T = 1.f;
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    // outputs
    int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];
    float* sorted_depths = &sorted_depths_[pix_id*MAX_SORTED_SPLATS];
    int intersect_count = 0;

    // collect and process batches of gaussians
    // each thread loads one gaussian at a time before rasterizing its
    // designated pixel
    int tr = block.thread_rank();
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
            id_batch[tr] = g_id;
            position_batch[tr] = {pos.x, pos.y, pos.z};
            axes_uv_batch[tr] = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};
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
            if (!std::isfinite(poi.z))
                continue;
            float alpha;
            if (!get_alpha(uv, opac, alpha))
                continue;

            sorted_indices[intersect_count] = (int32_t)id_batch[t];
            sorted_depths[intersect_count] = poi.z;
            if (++intersect_count >= MAX_SORTED_SPLATS){
                done = true;
                break;
            }

            const float next_T = T * (1.f - alpha);

            T = next_T;
            cur_idx = batch_start + t;
            if (T <= 1e-3f) {
                done = true;
                break;
            }
        }
    }

    if (inside) {
        num_intersects[pix_id] = intersect_count;

        // mark this so no need for num_intersects in rendering
        if (intersect_count < MAX_SORTED_SPLATS)
            sorted_indices[intersect_count] = SORTED_INDEX_INF;

        #if 0
            while (intersect_count < MAX_SORTED_SPLATS)
                sorted_depths[intersect_count++] = 1e6f;
        #endif
    }

}




inline __device__ void sort_per_pixel_insertion(int n, int32_t* indices, float* depths) {
    for (int i = 1; i < n; ++i) {
        int gid = indices[i];
        float z = depths[i];
        int j = i-1;

        while (j >= 0 && depths[j] > z) {
            depths[j+1] = depths[j];
            indices[j+1] = indices[j];
            j--;
        }
        depths[j+1] = z;
        indices[j+1] = gid;
    }
}


inline __device__ void sort_per_pixel_quick(int32_t* indices, float* depths, int i0, int i1) {
    if (i1-i0+1 < 8) {
        return sort_per_pixel_insertion(i1-i0+1, indices+i0, depths+i0);
    }

    // split into two arrays by value at center
    int ic = (i0+i1) / 2;
    float z = depths[ic];
    int pi = i0-1;
    for (int j = i0; j <= i1; j++) {
        if (depths[j] < z) {
            if (++pi == j) continue;
            // now there's always pi < j and j != ic
            // swap pi and j
            float z_ = depths[pi]; depths[pi] = depths[j]; depths[j] = z_;
            int32_t i_ = indices[pi]; indices[pi] = indices[j]; indices[j] = i_;
            if (ic == pi) ic = j;
        }
    }
    ++pi;
    // move ic to pi
    if (pi != ic) {
        // float z_ = depths[pi]; depths[pi] = depths[ic]; depths[ic] = z_;
        depths[ic] = depths[pi]; depths[pi] = z;
        int32_t i_ = indices[pi]; indices[pi] = indices[ic]; indices[ic] = i_;
    }

    // debug
    #if 0
        if (depths[pi] != z) printf("@");
        if (pi > ic) printf("#");
        for (int i = i0; i <= pi-1; i++)
            if (depths[i] > depths[pi]) printf("l");
        for (int i = pi+1; i <= i1; i++)
            if (depths[i] < depths[pi]) printf("r");
    #endif

    // warning: may cause stack overflow on some devices in a worst case
    // i.e. "CUDA error: an illegal memory access was encountered"
    if ((pi-1) - i0 > 0)
        sort_per_pixel_quick(indices, depths, i0, pi-1);
    if (i1 - (pi+1) > 0)
        sort_per_pixel_quick(indices, depths, pi+1, i1);
}


inline __device__ void _sort_per_pixel_heap_heapify(int n, int32_t* indices, float* depths, int i) {
    float z_i = depths[i];

    while (true) {
        int largest = i;
        float z_largest = z_i;
        int l = 2*i+1;
        int r = 2*i+2;

        float z_l = depths[l];
        if (l < n && z_l > z_largest) {
            largest = l; z_largest = z_l;
        }

        float z_r = depths[r];
        if (r < n && z_r > z_largest) {
            largest = r; z_largest = z_r;
        }

        if (largest == i) break;

        // swap i and largest
        depths[i] = z_largest; depths[largest] = z_i;
        float i_ = indices[i]; indices[i] = indices[largest]; indices[largest] = i_;

        i = largest;
    }
}

inline __device__ void sort_per_pixel_heap(int n, int32_t* indices, float* depths) {
    for (int i = n/2-1; i >= 0; --i) {
        _sort_per_pixel_heap_heapify(n, indices, depths, i);
    }

    for (int i = n-1; i > 0; --i) {
        // swap 0 and i
        float z_ = depths[0]; depths[0] = depths[i]; depths[i] = z_;
        float i_ = indices[0]; indices[0] = indices[i]; indices[i] = i_;

        _sort_per_pixel_heap_heapify(i, indices, depths, 0);
    }
}


inline __device__ void sort_per_pixel_randomized_quick(int32_t* indices, float* depths, int i0, int i1) {
    if (i1-i0+1 < 8) {
        return sort_per_pixel_insertion(i1-i0+1, indices+i0, depths+i0);
    }

    // random number (FNV-1a hash)
    uint32_t hash = 2166136261u;
    hash ^= (uint32_t)i0;
    hash *= 16777619u;
    hash ^= (uint32_t)i1;
    hash *= 16777619u;
    uint32_t m = i1-i0+1;
    int ic = i0 + (int)(hash%m);

    // split into two arrays
    float z = depths[ic];
    int pi = i0-1;
    for (int j = i0; j <= i1; j++) {
        if (depths[j] < z) {
            if (++pi == j) continue;
            // now there's always pi < j and j != ic
            // swap pi and j
            float z_ = depths[pi]; depths[pi] = depths[j]; depths[j] = z_;
            int32_t i_ = indices[pi]; indices[pi] = indices[j]; indices[j] = i_;
            if (ic == pi) ic = j;
        }
    }
    ++pi;
    // move ic to pi
    if (pi != ic) {
        // float z_ = depths[pi]; depths[pi] = depths[ic]; depths[ic] = z_;
        depths[ic] = depths[pi]; depths[pi] = z;
        int32_t i_ = indices[pi]; indices[pi] = indices[ic]; indices[ic] = i_;
    }

    // warning: may cause stack overflow on some devices in a worst case
    // i.e. "CUDA error: an illegal memory access was encountered"
    if ((pi-1) - i0 > 0)
        sort_per_pixel_randomized_quick(indices, depths, i0, pi-1);
    if (i1 - (pi+1) > 0)
        sort_per_pixel_randomized_quick(indices, depths, pi+1, i1);
}



template <PerPixelSortType SORT_TYPE>
__global__ void sort_per_pixel_kernel(
    const dim3 tile_bounds,
    const dim3 img_size,
    const int* __restrict__ num_intersects,
    int32_t* __restrict__ indices_,
    float* __restrict__ depths_
) {
    auto block = cg::this_thread_block();
    int32_t tile_id =
        block.group_index().y * tile_bounds.x + block.group_index().x;
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    bool inside = (i < img_size.y && j < img_size.x);
    if (!inside) return;

    int32_t pix_id = i * img_size.x + j;

    int n = num_intersects[pix_id];
    if (n <= 1)
        return;
    int32_t* indices = &indices_[pix_id*MAX_SORTED_SPLATS];
    float* depths = &depths_[pix_id*MAX_SORTED_SPLATS];
    int intersect_count = 0;

    switch (SORT_TYPE)
    {
    case PerPixelSortType::InsertionSort:
        sort_per_pixel_insertion(n, indices, depths);
        break;
    case PerPixelSortType::QuickSort:
        sort_per_pixel_quick(indices, depths, 0, n-1);
        break;
    case PerPixelSortType::HeapSort:
        sort_per_pixel_heap(n, indices, depths);
        break;
    case PerPixelSortType::RandomizedQuickSort:
        sort_per_pixel_randomized_quick(indices, depths, 0, n-1);
        break;
    }
}



__global__ void rasterize_simple_sorted_forward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float3& __restrict__ background,
    float3* __restrict__ out_img,
    float* __restrict__ out_alpha
) {
    auto block = cg::this_thread_block();
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    bool inside = (i < img_size.y && j < img_size.x);
    if (!inside) return;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    glm::vec2 pos_screen = { (float)j + 0.5f, (float)i + 0.5f };
    glm::vec2 pos_2d = { (pos_screen.x-cx)/fx, (pos_screen.y-cy)/fy };
    int32_t pix_id = i * img_size.x + j;

    // list of indices
    const int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];

    // current visibility left to render
    float T = 1.f;

    // rasterize
    float3 pix_out = {0.f, 0.f, 0.f};
    for (int cur_idx = 0; cur_idx < MAX_SORTED_SPLATS; cur_idx++) {
        int g_id = sorted_indices[cur_idx];
        if (g_id == SORTED_INDEX_INF)
            break;

        const glm::vec3 pos = *(glm::vec3*)&positions[g_id];
        const float opac = opacities[g_id];
        const glm::vec3 color = *(glm::vec3*)&colors[g_id];
        const float3 v0 = axes_u[g_id];
        const float3 v1 = axes_v[g_id];
        glm::mat2x3 axis_uv = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};

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
        pix_out.x = pix_out.x + color.x * vis;
        pix_out.y = pix_out.y + color.y * vis;
        pix_out.z = pix_out.z + color.z * vis;
        T = next_T;
        if (T <= 1e-3f) break;
    }

    if (inside) {
        float3 final_color;
        final_color.x = pix_out.x + T * background.x;
        final_color.y = pix_out.y + T * background.y;
        final_color.z = pix_out.z + T * background.z;
        out_img[pix_id] = final_color;
        out_alpha[pix_id] = 1.0f - T;
    }
}


__global__ void rasterize_simple_sorted_backward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int* __restrict__ num_intersects,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float* __restrict__ opacities,
    const float3& __restrict__ background,
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
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    bool inside = (i < img_size.y && j < img_size.x);
    if (!inside) return;
    
    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    glm::vec2 pos_screen = { (float)j + 0.5f, (float)i + 0.5f };
    glm::vec2 pos_2d = { (pos_screen.x-cx)/fx, (pos_screen.y-cy)/fy };
    int32_t pix_id = i * img_size.x + j;
    
    int n = num_intersects[pix_id];
    if (n == 0) return;

    // this is the T AFTER the last gaussian in this pixel
    float T_final = 1.0f - output_alpha[pix_id];
    float T = T_final;
    // the contribution from gaussians behind the current one
    float3 buffer = {0.f, 0.f, 0.f};
    // index of last gaussian to contribute to this pixel

    // df/d_out for this pixel
    const float3 v_out = nan_to_num(v_output[pix_id]);
    const float v_out_alpha = nan_to_num(v_output_alpha[pix_id]);

    // rasterize
    const int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];
    for (int cur_idx = n-1; cur_idx >= 0; cur_idx--) {
        int g_id = sorted_indices[cur_idx];

        const glm::vec3 pos = *(glm::vec3*)&positions[g_id];
        const float opac = opacities[g_id];
        const glm::vec3 color = *(glm::vec3*)&colors[g_id];
        const float3 v0 = axes_u[g_id];
        const float3 v1 = axes_v[g_id];
        glm::mat2x3 axis_uv = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};

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

        glm::vec3 v_position_local = {0.f, 0.f, 0.f};
        glm::vec2 v_position_xy_abs_local = {0.f, 0.f};
        glm::vec3 v_axis_u_local = {0.f, 0.f, 0.f};
        glm::vec3 v_axis_v_local = {0.f, 0.f, 0.f};
        glm::vec3 v_color_local = {0.f, 0.f, 0.f};
        float v_opacity_local = 0.f;

        // compute the current T for this gaussian
        const float ra = 1.f / (1.f - alpha);
        const float next_T = T * ra;
        const float vis = alpha * next_T;

        // update v_rgb for this gaussian
        float v_alpha = 0.f;
        v_color_local = {vis * v_out.x, vis * v_out.y, vis * v_out.z};

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
            uv, opac,
            v_alpha, v_uv, v_opacity_local
        );
        glm::mat2x3 v_axis_uv = glm::mat2x3(0.0f);
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

        float* v_position_ptr = (float*)(v_positions);
        atomicAdd(v_position_ptr + 3*g_id + 0, v_position_local.x);
        atomicAdd(v_position_ptr + 3*g_id + 1, v_position_local.y);
        atomicAdd(v_position_ptr + 3*g_id + 2, v_position_local.z);
        float* v_positions_xy_abs_ptr = (float*)(v_positions_xy_abs);
        atomicAdd(v_positions_xy_abs_ptr + 2*g_id + 0, v_position_xy_abs_local.x);
        atomicAdd(v_positions_xy_abs_ptr + 2*g_id + 1, v_position_xy_abs_local.y);

        float* v_axis_u_ptr = (float*)(v_axes_u);
        atomicAdd(v_axis_u_ptr + 3*g_id + 0, v_axis_u_local.x);
        atomicAdd(v_axis_u_ptr + 3*g_id + 1, v_axis_u_local.y);
        atomicAdd(v_axis_u_ptr + 3*g_id + 2, v_axis_u_local.z);
        float* v_axis_v_ptr = (float*)(v_axes_v);
        atomicAdd(v_axis_v_ptr + 3*g_id + 0, v_axis_v_local.x);
        atomicAdd(v_axis_v_ptr + 3*g_id + 1, v_axis_v_local.y);
        atomicAdd(v_axis_v_ptr + 3*g_id + 2, v_axis_v_local.z);

        float* v_color_ptr = (float*)(v_colors);
        atomicAdd(v_color_ptr + 3*g_id + 0, v_color_local.x);
        atomicAdd(v_color_ptr + 3*g_id + 1, v_color_local.y);
        atomicAdd(v_color_ptr + 3*g_id + 2, v_color_local.z);
        
        atomicAdd(v_opacities + g_id, v_opacity_local);
    }

}



template <DepthMode DEPTH_MODE>
__global__ void rasterize_depth_sorted_forward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    int* __restrict__ final_index,
    float* __restrict__ out_depth,
    float2* __restrict__ out_visibility
) {
    auto block = cg::this_thread_block();
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    bool inside = (i < img_size.y && j < img_size.x);
    if (!inside) return;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    glm::vec2 pos_screen = { (float)j + 0.5f, (float)i + 0.5f };
    glm::vec2 pos_2d = { (pos_screen.x-cx)/fx, (pos_screen.y-cy)/fy };
    int32_t pix_id = i * img_size.x + j;

    // list of indices
    const int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];

    // current visibility left to render
    float T = 1.f;
    float interp = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    // rasterize
    float output_depth = 0.0f;
    float output_visibility = 0.0f;

    for (; cur_idx < MAX_SORTED_SPLATS; cur_idx++) {
        int g_id = sorted_indices[cur_idx];
        if (g_id == SORTED_INDEX_INF)
            break;

        const glm::vec3 pos = *(glm::vec3*)&positions[g_id];
        const float opac = opacities[g_id];
        const float3 v0 = axes_u[g_id];
        const float3 v1 = axes_v[g_id];
        glm::mat2x3 axis_uv = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};

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
                cur_idx++;
                break;
            }
            output_depth = next_depth;

        }  // DEPTH_MODE

        T = next_T;
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
__global__ void rasterize_depth_sorted_backward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const int* __restrict__ final_index,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    const float* __restrict__ out_depth,
    const float2* __restrict__ out_visibility,
    const float* __restrict__ v_out_depth,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float* __restrict__ v_opacities
) {
    auto block = cg::this_thread_block();
    unsigned i =
        block.group_index().y * block.group_dim().y + block.thread_index().y;
    unsigned j =
        block.group_index().x * block.group_dim().x + block.thread_index().x;

    bool inside = (i < img_size.y && j < img_size.x);
    if (!inside) return;
    
    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;
    glm::vec2 pos_screen = { (float)j + 0.5f, (float)i + 0.5f };
    glm::vec2 pos_2d = { (pos_screen.x-cx)/fx, (pos_screen.y-cy)/fy };
    int32_t pix_id = i * img_size.x + j;
    
    int n = final_index[pix_id];
    if (n == 0) return;

    // this is the T AFTER the last gaussian in this pixel
    glm::vec2 meta_out = *(glm::vec2*)&out_visibility[pix_id];
    float T_final = meta_out.x;
    float T = T_final;
    float v_T = 0.0f;
    const float interp = meta_out.y;

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

    // rasterize
    const int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];
    for (int cur_idx = n-1; cur_idx >= 0; cur_idx--) {
        int g_id = sorted_indices[cur_idx];

        const glm::vec3 pos = *(glm::vec3*)&positions[g_id];
        const float opac = opacities[g_id];
        const float3 v0 = axes_u[g_id];
        const float3 v1 = axes_v[g_id];
        glm::mat2x3 axis_uv = {v0.x, v0.y, v0.z, v1.x, v1.y, v1.z};

        glm::vec3 poi;
        glm::vec2 uv;
        get_intersection(pos, axis_uv, pos_2d, poi, uv);
        if (glm::length(uv) > visibility_kernel_radius())
            continue;
        float alpha;
        if (!get_alpha(uv, opac, alpha))
            continue;

        glm::vec3 v_position_local = {0.f, 0.f, 0.f};
        glm::vec2 v_position_xy_abs_local = {0.f, 0.f};
        glm::vec3 v_axis_u_local = {0.f, 0.f, 0.f};
        glm::vec3 v_axis_v_local = {0.f, 0.f, 0.f};
        float v_opacity_local = 0.f;

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
        glm::vec2 v_uv;
        get_alpha_vjp(
            uv, opac,
            v_alpha, v_uv, v_opacity_local
        );
        glm::mat2x3 v_axis_uv = glm::mat2x3(0.0f);
        float v_depth_raw = depth_map_vjp(depth_raw, v_depth);
        get_intersection_vjp(
            pos, axis_uv, pos_2d,
            {0.f, 0.f, v_depth_raw}, v_uv,
            v_position_local, v_axis_uv
        );
        v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local));
        v_axis_u_local = v_axis_uv[0];
        v_axis_v_local = v_axis_uv[1];
        
        float* v_position_ptr = (float*)(v_positions);
        atomicAdd(v_position_ptr + 3*g_id + 0, v_position_local.x);
        atomicAdd(v_position_ptr + 3*g_id + 1, v_position_local.y);
        atomicAdd(v_position_ptr + 3*g_id + 2, v_position_local.z);
        float* v_positions_xy_abs_ptr = (float*)(v_positions_xy_abs);
        atomicAdd(v_positions_xy_abs_ptr + 2*g_id + 0, v_position_xy_abs_local.x);
        atomicAdd(v_positions_xy_abs_ptr + 2*g_id + 1, v_position_xy_abs_local.y);

        float* v_axis_u_ptr = (float*)(v_axes_u);
        atomicAdd(v_axis_u_ptr + 3*g_id + 0, v_axis_u_local.x);
        atomicAdd(v_axis_u_ptr + 3*g_id + 1, v_axis_u_local.y);
        atomicAdd(v_axis_u_ptr + 3*g_id + 2, v_axis_u_local.z);
        float* v_axis_v_ptr = (float*)(v_axes_v);
        atomicAdd(v_axis_v_ptr + 3*g_id + 0, v_axis_v_local.x);
        atomicAdd(v_axis_v_ptr + 3*g_id + 1, v_axis_v_local.y);
        atomicAdd(v_axis_v_ptr + 3*g_id + 2, v_axis_v_local.z);
        
        float v_opacity_local_ = (float)v_opacity_local;
        atomicAdd(v_opacities + g_id, v_opacity_local_);
    }
}




template __global__ void sort_per_pixel_kernel<PerPixelSortType::InsertionSort>(
    const dim3 tile_bounds,
    const dim3 img_size,
    const int* __restrict__ num_intersects,
    int32_t* __restrict__ indices_,
    float* __restrict__ depths_
);

template __global__ void sort_per_pixel_kernel<PerPixelSortType::QuickSort>(
    const dim3 tile_bounds,
    const dim3 img_size,
    const int* __restrict__ num_intersects,
    int32_t* __restrict__ indices_,
    float* __restrict__ depths_
);

template __global__ void sort_per_pixel_kernel<PerPixelSortType::HeapSort>(
    const dim3 tile_bounds,
    const dim3 img_size,
    const int* __restrict__ num_intersects,
    int32_t* __restrict__ indices_,
    float* __restrict__ depths_
);

template __global__ void sort_per_pixel_kernel<PerPixelSortType::RandomizedQuickSort>(
    const dim3 tile_bounds,
    const dim3 img_size,
    const int* __restrict__ num_intersects,
    int32_t* __restrict__ indices_,
    float* __restrict__ depths_
);

template __global__ void rasterize_depth_sorted_forward_kernel<DepthMode::Mean>(
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    int* __restrict__ final_index,
    float* __restrict__ out_depth,
    float2* __restrict__ out_visibility
);

template __global__ void rasterize_depth_sorted_forward_kernel<DepthMode::Median>(
    const dim3 img_size,
    const float4 intrins,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    int* __restrict__ final_index,
    float* __restrict__ out_depth,
    float2* __restrict__ out_visibility
);

template __global__ void rasterize_depth_sorted_backward_kernel<DepthMode::Mean>(
    const dim3 img_size,
    const float4 intrins,
    const int* __restrict__ num_intersects,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    const float* __restrict__ out_depth,
    const float2* __restrict__ out_visibility,
    const float* __restrict__ v_out_depth,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float* __restrict__ v_opacities
);

template __global__ void rasterize_depth_sorted_backward_kernel<DepthMode::Median>(
    const dim3 img_size,
    const float4 intrins,
    const int* __restrict__ num_intersects,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float* __restrict__ opacities,
    const float* __restrict__ out_depth,
    const float2* __restrict__ out_visibility,
    const float* __restrict__ v_out_depth,
    float3* __restrict__ v_positions,
    float2* __restrict__ v_positions_xy_abs,
    float3* __restrict__ v_axes_u,
    float3* __restrict__ v_axes_v,
    float* __restrict__ v_opacities
);
