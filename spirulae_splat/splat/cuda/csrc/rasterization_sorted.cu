#include "helpers.cuh"
#include "ch.cuh"
#include "rasterization_sorted.cuh"
#include <algorithm>

#include "stdio.h"


template<CameraType CAMERA_TYPE>
__global__ void rasterize_indices_kernel(
    _ARGS_rasterize_indices_kernel
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

    __shared__ int id_batch[MAX_BLOCK_SIZE];
    __shared__ glm::vec3 position_batch[MAX_BLOCK_SIZE];
    __shared__ glm::mat2x3 axes_uv_batch[MAX_BLOCK_SIZE];
    __shared__ float opacity_batch[MAX_BLOCK_SIZE];

    // current visibility left to render
    float T = 1.f;

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


inline __device__ void _sort_per_pixel_network_compare_swap(int32_t *indices, float* depths, int i, int j) {
    float u = depths[i], v = depths[j];
    if (u > v) {
        depths[j] = u, depths[i] = v;
        int32_t k = indices[i]; indices[i] = indices[j], indices[j] = k;
    }
}

inline __device__ void sort_per_pixel_network(int n, int32_t* indices, float* depths) {
    // https://bertdobbelaere.github.io/sorting_networks.html
    // s = """<paste sorting network>"""
    // for i, j in sum([eval(x) for x in s.strip().split('\n')], []):
    //     print(f"_({i},{j})", end='')
    // print(' ')
    #define _(i,j) _sort_per_pixel_network_compare_swap(indices,depths,i,j);
    if (n <= 4) {
        if (n == 2) { _(0,1) }
        else if (n == 3) { _(0,2)_(0,1)_(1,2) }
        else { _(0,2)_(1,3)_(0,1)_(2,3)_(1,2) }
    }
    else if (n <= 8) {
        if (n == 5) { _(0,3)_(1,4)_(0,2)_(1,3)_(0,1)_(2,4)_(1,2)_(3,4)_(2,3) }
        else if (n == 6) { _(0,5)_(1,3)_(2,4)_(1,2)_(3,4)_(0,3)_(2,5)_(0,1)_(2,3)_(4,5)_(1,2)_(3,4) }
        else if (n == 7) { _(0,6)_(2,3)_(4,5)_(0,2)_(1,4)_(3,6)_(0,1)_(2,5)_(3,4)_(1,2)_(4,6)_(2,3)_(4,5)_(1,2)_(3,4)_(5,6) }
        else { _(0,2)_(1,3)_(4,6)_(5,7)_(0,4)_(1,5)_(2,6)_(3,7)_(0,1)_(2,3)_(4,5)_(6,7)_(2,4)_(3,5)_(1,4)_(3,6)_(1,2)_(3,4)_(5,6) }
    }
#if MAX_SORTED_SPLATS > 8
    else if (n <= 12) {
        if (n == 9) { _(0,3)_(1,7)_(2,5)_(4,8)_(0,7)_(2,4)_(3,8)_(5,6)_(0,2)_(1,3)_(4,5)_(7,8)_(1,4)_(3,6)_(5,7)_(0,1)_(2,4)_(3,5)_(6,8)_(2,3)_(4,5)_(6,7)_(1,2)_(3,4)_(5,6) }
        else if (n == 10) { _(0,8)_(1,9)_(2,7)_(3,5)_(4,6)_(0,2)_(1,4)_(5,8)_(7,9)_(0,3)_(2,4)_(5,7)_(6,9)_(0,1)_(3,6)_(8,9)_(1,5)_(2,3)_(4,8)_(6,7)_(1,2)_(3,5)_(4,6)_(7,8)_(2,3)_(4,5)_(6,7)_(3,4)_(5,6) }
        else if (n == 11) { _(0,9)_(1,6)_(2,4)_(3,7)_(5,8)_(0,1)_(3,5)_(4,10)_(6,9)_(7,8)_(1,3)_(2,5)_(4,7)_(8,10)_(0,4)_(1,2)_(3,7)_(5,9)_(6,8)_(0,1)_(2,6)_(4,5)_(7,8)_(9,10)_(2,4)_(3,6)_(5,7)_(8,9)_(1,2)_(3,4)_(5,6)_(7,8)_(2,3)_(4,5)_(6,7) }
        else { _(0,8)_(1,7)_(2,6)_(3,11)_(4,10)_(5,9)_(0,1)_(2,5)_(3,4)_(6,9)_(7,8)_(10,11)_(0,2)_(1,6)_(5,10)_(9,11)_(0,3)_(1,2)_(4,6)_(5,7)_(8,11)_(9,10)_(1,4)_(3,5)_(6,8)_(7,10)_(1,3)_(2,5)_(6,9)_(8,10)_(2,3)_(4,5)_(6,7)_(8,9)_(4,6)_(5,7)_(3,4)_(5,6)_(7,8) }
    }
#endif
#if MAX_SORTED_SPLATS > 12
    else if (n <= 16) {
        if (n == 13) { _(0,11)_(1,7)_(2,4)_(3,5)_(8,9)_(10,12)_(0,2)_(3,6)_(4,12)_(5,7)_(8,10)_(0,8)_(1,3)_(2,5)_(4,9)_(6,11)_(7,12)_(0,1)_(2,10)_(3,8)_(4,6)_(9,11)_(1,3)_(2,4)_(5,10)_(6,8)_(7,9)_(11,12)_(1,2)_(3,4)_(5,8)_(6,9)_(7,10)_(2,3)_(4,7)_(5,6)_(8,11)_(9,10)_(4,5)_(6,7)_(8,9)_(10,11)_(3,4)_(5,6)_(7,8)_(9,10) }
        else if (n == 14) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(0,2)_(1,3)_(4,8)_(5,9)_(10,12)_(11,13)_(0,4)_(1,2)_(3,7)_(5,8)_(6,10)_(9,13)_(11,12)_(0,6)_(1,5)_(3,9)_(4,10)_(7,13)_(8,12)_(2,10)_(3,11)_(4,6)_(7,9)_(1,3)_(2,8)_(5,11)_(6,7)_(10,12)_(1,4)_(2,6)_(3,5)_(7,11)_(8,10)_(9,12)_(2,4)_(3,6)_(5,8)_(7,10)_(9,11)_(3,4)_(5,6)_(7,8)_(9,10)_(6,7) }
        else if (n == 15) { _(1,2)_(3,10)_(4,14)_(5,8)_(6,13)_(7,12)_(9,11)_(0,14)_(1,5)_(2,8)_(3,7)_(6,9)_(10,12)_(11,13)_(0,7)_(1,6)_(2,9)_(4,10)_(5,11)_(8,13)_(12,14)_(0,6)_(2,4)_(3,5)_(7,11)_(8,10)_(9,12)_(13,14)_(0,3)_(1,2)_(4,7)_(5,9)_(6,8)_(10,11)_(12,13)_(0,1)_(2,3)_(4,6)_(7,9)_(10,12)_(11,13)_(1,2)_(3,5)_(8,10)_(11,12)_(3,4)_(5,6)_(7,8)_(9,10)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(5,6)_(7,8) }
        else { _(0,13)_(1,12)_(2,15)_(3,14)_(4,8)_(5,6)_(7,11)_(9,10)_(0,5)_(1,7)_(2,9)_(3,4)_(6,13)_(8,14)_(10,15)_(11,12)_(0,1)_(2,3)_(4,5)_(6,8)_(7,9)_(10,11)_(12,13)_(14,15)_(0,2)_(1,3)_(4,10)_(5,11)_(6,7)_(8,9)_(12,14)_(13,15)_(1,2)_(3,12)_(4,6)_(5,7)_(8,10)_(9,11)_(13,14)_(1,4)_(2,6)_(5,8)_(7,10)_(9,13)_(11,14)_(2,4)_(3,6)_(9,12)_(11,13)_(3,5)_(6,8)_(7,9)_(10,12)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(6,7)_(8,9) }
    }
#endif
#if MAX_SORTED_SPLATS > 16
    else if (n <= 20) {
        if (n == 17) { _(0,11)_(1,15)_(2,10)_(3,5)_(4,6)_(8,12)_(9,16)_(13,14)_(0,6)_(1,13)_(2,8)_(4,14)_(5,15)_(7,11)_(0,8)_(3,7)_(4,9)_(6,16)_(10,11)_(12,14)_(0,2)_(1,4)_(5,6)_(7,13)_(8,9)_(10,12)_(11,14)_(15,16)_(0,3)_(2,5)_(6,11)_(7,10)_(9,13)_(12,15)_(14,16)_(0,1)_(3,4)_(5,10)_(6,9)_(7,8)_(11,15)_(13,14)_(1,2)_(3,7)_(4,8)_(6,12)_(11,13)_(14,15)_(1,3)_(2,7)_(4,5)_(9,11)_(10,12)_(13,14)_(2,3)_(4,6)_(5,7)_(8,10)_(3,4)_(6,8)_(7,9)_(10,12)_(5,6)_(7,8)_(9,10)_(11,12)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13) }
        else if (n == 18) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(0,2)_(1,3)_(4,12)_(5,13)_(6,8)_(9,11)_(14,16)_(15,17)_(0,14)_(1,16)_(2,15)_(3,17)_(0,6)_(1,10)_(2,9)_(7,16)_(8,15)_(11,17)_(1,4)_(3,9)_(5,7)_(8,14)_(10,12)_(13,16)_(0,1)_(2,5)_(3,13)_(4,14)_(7,9)_(8,10)_(12,15)_(16,17)_(1,2)_(3,5)_(4,6)_(11,13)_(12,14)_(15,16)_(4,8)_(5,12)_(6,10)_(7,11)_(9,13)_(1,4)_(2,8)_(3,6)_(5,7)_(9,15)_(10,12)_(11,14)_(13,16)_(2,4)_(5,8)_(6,10)_(7,11)_(9,12)_(13,15)_(3,5)_(6,8)_(7,10)_(9,11)_(12,14)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14) }
        else if (n == 19) { _(0,12)_(1,4)_(2,8)_(3,5)_(6,17)_(7,11)_(9,14)_(10,13)_(15,16)_(0,2)_(1,7)_(3,6)_(4,11)_(5,17)_(8,12)_(10,15)_(13,16)_(14,18)_(3,10)_(4,14)_(5,15)_(6,13)_(7,9)_(11,17)_(16,18)_(0,7)_(1,10)_(4,6)_(9,15)_(11,16)_(12,17)_(13,14)_(0,3)_(2,6)_(5,7)_(8,11)_(12,16)_(1,8)_(2,9)_(3,4)_(6,15)_(7,13)_(10,11)_(12,18)_(1,3)_(2,5)_(6,9)_(7,12)_(8,10)_(11,14)_(17,18)_(0,1)_(2,3)_(4,8)_(6,10)_(9,12)_(14,15)_(16,17)_(1,2)_(5,8)_(6,7)_(9,11)_(10,13)_(14,16)_(15,17)_(3,6)_(4,5)_(7,9)_(8,10)_(11,12)_(13,14)_(15,16)_(3,4)_(5,6)_(7,8)_(9,10)_(11,13)_(12,14)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15) }
        else { _(0,3)_(1,7)_(2,5)_(4,8)_(6,9)_(10,13)_(11,15)_(12,18)_(14,17)_(16,19)_(0,14)_(1,11)_(2,16)_(3,17)_(4,12)_(5,19)_(6,10)_(7,15)_(8,18)_(9,13)_(0,4)_(1,2)_(3,8)_(5,7)_(11,16)_(12,14)_(15,19)_(17,18)_(1,6)_(2,12)_(3,5)_(4,11)_(7,17)_(8,15)_(13,18)_(14,16)_(0,1)_(2,6)_(7,10)_(9,12)_(13,17)_(18,19)_(1,6)_(5,9)_(7,11)_(8,12)_(10,14)_(13,18)_(3,5)_(4,7)_(8,10)_(9,11)_(12,15)_(14,16)_(1,3)_(2,4)_(5,7)_(6,10)_(9,13)_(12,14)_(15,17)_(16,18)_(1,2)_(3,4)_(6,7)_(8,9)_(10,11)_(12,13)_(15,16)_(17,18)_(2,3)_(4,6)_(5,8)_(7,9)_(10,12)_(11,14)_(13,15)_(16,17)_(4,5)_(6,8)_(7,10)_(9,12)_(11,13)_(14,15)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16) }
    }
#endif
#if MAX_SORTED_SPLATS > 20
    else if (n <= 24) {
        if (n == 21) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,14)_(13,15)_(16,18)_(17,19)_(0,8)_(1,9)_(2,10)_(3,11)_(4,12)_(5,13)_(6,14)_(7,15)_(0,4)_(1,5)_(3,7)_(6,20)_(8,12)_(9,13)_(10,14)_(15,19)_(2,6)_(3,18)_(7,20)_(2,16)_(3,6)_(5,18)_(7,17)_(11,20)_(0,2)_(3,8)_(6,12)_(7,10)_(9,16)_(11,15)_(13,17)_(14,18)_(19,20)_(1,7)_(2,3)_(4,9)_(10,11)_(13,16)_(15,18)_(17,19)_(1,4)_(5,10)_(6,13)_(7,8)_(11,14)_(12,16)_(15,17)_(18,19)_(1,2)_(3,4)_(5,6)_(10,12)_(11,13)_(14,16)_(17,18)_(2,3)_(4,5)_(6,9)_(10,11)_(12,13)_(14,15)_(16,17)_(6,7)_(8,9)_(15,16)_(4,6)_(7,8)_(9,12)_(13,15)_(3,4)_(5,7)_(8,10)_(9,11)_(12,14)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14) }
        else if (n == 22) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(0,2)_(1,3)_(4,6)_(5,7)_(8,12)_(9,13)_(14,16)_(15,17)_(18,20)_(19,21)_(0,4)_(1,5)_(2,6)_(3,7)_(8,10)_(9,12)_(11,13)_(14,18)_(15,19)_(16,20)_(17,21)_(0,14)_(1,15)_(2,18)_(3,19)_(4,16)_(5,17)_(6,20)_(7,21)_(9,11)_(10,12)_(2,8)_(3,11)_(6,9)_(10,18)_(12,15)_(13,19)_(0,2)_(1,10)_(3,16)_(5,18)_(6,14)_(7,15)_(8,12)_(9,13)_(11,20)_(19,21)_(2,6)_(3,10)_(4,8)_(5,12)_(9,16)_(11,18)_(13,17)_(15,19)_(1,4)_(7,13)_(8,14)_(9,12)_(17,20)_(1,2)_(3,8)_(4,6)_(7,11)_(10,14)_(13,18)_(15,17)_(19,20)_(2,4)_(5,10)_(7,9)_(11,16)_(12,14)_(17,19)_(5,6)_(7,8)_(9,11)_(10,12)_(13,14)_(15,16)_(3,5)_(6,7)_(8,10)_(9,12)_(11,13)_(14,15)_(16,18)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18) }
        else if (n == 23) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,14)_(13,15)_(16,18)_(17,19)_(21,22)_(0,4)_(1,5)_(2,6)_(3,7)_(8,12)_(9,13)_(10,14)_(11,15)_(17,21)_(18,20)_(19,22)_(0,8)_(1,9)_(2,10)_(3,11)_(4,12)_(5,13)_(6,14)_(7,15)_(1,2)_(5,18)_(7,19)_(9,16)_(10,21)_(12,20)_(15,22)_(5,9)_(6,7)_(10,18)_(11,21)_(12,17)_(13,20)_(14,15)_(3,17)_(6,16)_(7,14)_(8,12)_(15,19)_(20,21)_(3,4)_(5,8)_(6,10)_(9,12)_(13,16)_(14,15)_(17,18)_(19,21)_(0,5)_(1,8)_(2,12)_(3,9)_(4,10)_(7,13)_(11,17)_(14,16)_(18,20)_(2,6)_(3,5)_(4,8)_(7,11)_(10,12)_(13,18)_(14,17)_(15,20)_(1,3)_(2,5)_(6,9)_(7,10)_(11,13)_(12,14)_(15,18)_(16,17)_(19,20)_(2,3)_(4,6)_(8,9)_(11,12)_(13,14)_(15,16)_(17,19)_(3,4)_(5,6)_(7,8)_(9,10)_(12,13)_(14,15)_(17,18)_(4,5)_(6,7)_(8,9)_(10,11)_(16,17) }
        else { _(0,20)_(1,12)_(2,16)_(3,23)_(4,6)_(5,10)_(7,21)_(8,14)_(9,15)_(11,22)_(13,18)_(17,19)_(0,3)_(1,11)_(2,7)_(4,17)_(5,13)_(6,19)_(8,9)_(10,18)_(12,22)_(14,15)_(16,21)_(20,23)_(0,1)_(2,4)_(3,12)_(5,8)_(6,9)_(7,10)_(11,20)_(13,16)_(14,17)_(15,18)_(19,21)_(22,23)_(2,5)_(4,8)_(6,11)_(7,14)_(9,16)_(12,17)_(15,19)_(18,21)_(1,8)_(3,14)_(4,7)_(9,20)_(10,12)_(11,13)_(15,22)_(16,19)_(0,7)_(1,5)_(3,4)_(6,11)_(8,15)_(9,14)_(10,13)_(12,17)_(16,23)_(18,22)_(19,20)_(0,2)_(1,6)_(4,7)_(5,9)_(8,10)_(13,15)_(14,18)_(16,19)_(17,22)_(21,23)_(2,3)_(4,5)_(6,8)_(7,9)_(10,11)_(12,13)_(14,16)_(15,17)_(18,19)_(20,21)_(1,2)_(3,6)_(4,10)_(7,8)_(9,11)_(12,14)_(13,19)_(15,16)_(17,20)_(21,22)_(2,3)_(5,10)_(6,7)_(8,9)_(13,18)_(14,15)_(16,17)_(20,21)_(3,4)_(5,7)_(10,12)_(11,13)_(16,18)_(19,20)_(4,6)_(8,10)_(9,12)_(11,14)_(13,15)_(17,19)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18) }
    }
#endif
#if MAX_SORTED_SPLATS > 24
    else if (n <= 28) {
        if (n == 25) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(22,23)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,14)_(13,15)_(16,18)_(17,19)_(21,22)_(23,24)_(0,4)_(1,5)_(2,6)_(3,7)_(8,12)_(9,13)_(10,14)_(11,15)_(18,21)_(20,23)_(22,24)_(0,8)_(1,9)_(2,10)_(3,11)_(4,12)_(5,13)_(6,14)_(7,15)_(16,20)_(17,22)_(19,24)_(21,23)_(1,18)_(3,21)_(5,23)_(6,19)_(11,14)_(15,24)_(1,16)_(3,17)_(6,9)_(7,11)_(13,19)_(14,23)_(0,1)_(2,16)_(3,8)_(7,20)_(10,13)_(11,22)_(15,23)_(1,2)_(5,10)_(7,18)_(11,21)_(15,20)_(19,22)_(4,7)_(5,6)_(9,18)_(10,17)_(11,12)_(13,21)_(14,15)_(19,20)_(22,23)_(3,4)_(7,8)_(9,10)_(11,16)_(12,17)_(13,18)_(19,21)_(20,22)_(1,3)_(2,4)_(5,11)_(6,16)_(7,9)_(8,10)_(12,13)_(14,19)_(15,18)_(2,3)_(5,7)_(6,9)_(8,11)_(10,16)_(12,14)_(15,17)_(3,5)_(4,6)_(7,8)_(9,11)_(10,12)_(13,14)_(15,16)_(17,18)_(4,7)_(6,8)_(9,10)_(11,12)_(13,15)_(14,16)_(17,19)_(18,21)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21) }
        else if (n == 26) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(22,23)_(24,25)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(14,16)_(15,17)_(18,20)_(19,21)_(22,24)_(23,25)_(0,4)_(1,6)_(2,5)_(3,7)_(8,14)_(9,16)_(10,15)_(11,17)_(18,22)_(19,24)_(20,23)_(21,25)_(0,18)_(1,19)_(2,20)_(3,21)_(4,22)_(5,23)_(6,24)_(7,25)_(9,12)_(13,16)_(3,11)_(8,9)_(10,13)_(12,15)_(14,22)_(16,17)_(0,8)_(1,9)_(2,14)_(6,12)_(7,15)_(10,18)_(11,23)_(13,19)_(16,24)_(17,25)_(1,2)_(3,18)_(4,8)_(7,22)_(17,21)_(23,24)_(3,14)_(4,10)_(5,18)_(7,20)_(8,13)_(11,22)_(12,17)_(15,21)_(1,4)_(5,6)_(7,9)_(8,10)_(15,17)_(16,18)_(19,20)_(21,24)_(2,5)_(3,10)_(6,14)_(9,13)_(11,19)_(12,16)_(15,22)_(20,23)_(2,8)_(5,7)_(6,9)_(11,12)_(13,14)_(16,19)_(17,23)_(18,20)_(2,4)_(3,5)_(6,11)_(7,10)_(9,16)_(12,13)_(14,19)_(15,18)_(20,22)_(21,23)_(3,4)_(5,8)_(6,7)_(9,11)_(10,12)_(13,15)_(14,16)_(17,20)_(18,19)_(21,22)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18)_(19,20)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21) }
        else if (n == 27) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,14)_(15,16)_(17,18)_(19,20)_(21,22)_(23,24)_(25,26)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,13)_(15,17)_(16,18)_(19,21)_(20,22)_(23,25)_(24,26)_(0,23)_(1,24)_(2,25)_(3,26)_(4,8)_(5,9)_(6,10)_(7,11)_(13,14)_(15,19)_(16,20)_(17,21)_(18,22)_(0,4)_(1,6)_(2,19)_(3,20)_(5,13)_(9,21)_(11,14)_(12,16)_(17,23)_(18,24)_(22,26)_(5,17)_(6,16)_(7,22)_(9,25)_(10,24)_(12,15)_(13,20)_(14,26)_(1,12)_(4,15)_(7,23)_(10,19)_(11,16)_(13,18)_(20,24)_(22,25)_(0,1)_(6,12)_(8,11)_(9,15)_(10,17)_(14,24)_(16,21)_(18,19)_(1,4)_(2,8)_(3,11)_(12,15)_(14,20)_(16,22)_(21,25)_(2,5)_(3,17)_(8,13)_(11,23)_(21,22)_(24,25)_(1,2)_(3,10)_(5,6)_(7,13)_(11,15)_(14,21)_(18,23)_(20,22)_(4,5)_(6,9)_(7,8)_(13,17)_(14,16)_(19,23)_(22,24)_(2,4)_(3,6)_(5,7)_(8,12)_(9,10)_(11,13)_(14,18)_(15,17)_(16,19)_(21,23)_(3,5)_(6,8)_(7,9)_(10,12)_(11,14)_(13,16)_(15,18)_(17,19)_(20,21)_(22,23)_(5,6)_(8,11)_(9,10)_(12,14)_(13,15)_(17,18)_(19,21)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,20)_(21,22)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18)_(19,20) }
        else { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(22,23)_(24,25)_(26,27)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,14)_(13,15)_(16,18)_(17,19)_(20,22)_(21,23)_(24,26)_(25,27)_(0,4)_(1,5)_(2,6)_(3,7)_(8,12)_(9,13)_(14,18)_(15,19)_(20,24)_(21,25)_(22,26)_(23,27)_(0,20)_(1,21)_(2,22)_(3,23)_(4,24)_(5,25)_(6,26)_(7,27)_(9,17)_(10,18)_(11,15)_(12,16)_(1,2)_(4,20)_(5,6)_(7,23)_(8,12)_(9,16)_(10,14)_(11,18)_(13,17)_(15,19)_(21,22)_(25,26)_(0,8)_(1,9)_(2,12)_(3,5)_(4,10)_(6,16)_(7,13)_(11,21)_(14,20)_(15,25)_(17,23)_(18,26)_(19,27)_(22,24)_(2,4)_(3,7)_(5,17)_(8,14)_(9,11)_(10,22)_(13,19)_(16,18)_(20,24)_(23,25)_(1,8)_(3,9)_(5,11)_(6,10)_(7,15)_(12,20)_(16,22)_(17,21)_(18,24)_(19,26)_(1,2)_(4,6)_(5,9)_(10,16)_(11,17)_(12,14)_(13,15)_(18,22)_(21,23)_(25,26)_(4,8)_(6,12)_(7,11)_(10,14)_(13,17)_(15,21)_(16,20)_(19,23)_(2,4)_(6,8)_(7,16)_(9,14)_(10,12)_(11,20)_(13,18)_(15,17)_(19,21)_(23,25)_(3,10)_(5,12)_(7,9)_(11,13)_(14,16)_(15,22)_(17,24)_(18,20)_(3,6)_(5,8)_(7,10)_(9,12)_(11,14)_(13,16)_(15,18)_(17,20)_(19,22)_(21,24)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18)_(19,20)_(21,22)_(23,24) }
    }
    else {
        if (n == 29) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,28)_(17,26)_(18,25)_(19,23)_(21,27)_(22,24)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,14)_(13,15)_(17,22)_(18,19)_(20,27)_(23,25)_(24,26)_(0,4)_(1,5)_(2,6)_(3,7)_(8,12)_(9,13)_(10,14)_(11,15)_(16,20)_(17,18)_(19,22)_(23,24)_(25,26)_(27,28)_(0,8)_(1,9)_(2,10)_(3,11)_(4,12)_(5,13)_(6,14)_(7,15)_(16,18)_(20,22)_(21,25)_(24,27)_(26,28)_(1,8)_(2,24)_(4,19)_(5,20)_(6,21)_(7,27)_(9,18)_(10,23)_(11,26)_(13,22)_(14,25)_(15,28)_(16,17)_(0,6)_(2,4)_(3,24)_(5,10)_(12,19)_(13,18)_(14,21)_(15,25)_(20,23)_(26,27)_(0,16)_(1,6)_(3,12)_(4,8)_(5,17)_(7,24)_(14,20)_(15,26)_(18,21)_(19,23)_(25,27)_(1,5)_(2,16)_(3,10)_(6,9)_(7,18)_(8,17)_(11,19)_(13,14)_(15,22)_(21,23)_(25,26)_(1,2)_(3,5)_(4,8)_(6,16)_(7,11)_(9,17)_(10,12)_(14,20)_(15,18)_(19,24)_(22,27)_(4,6)_(9,16)_(10,13)_(11,19)_(12,14)_(20,21)_(22,26)_(23,24)_(2,4)_(3,6)_(7,16)_(8,9)_(11,17)_(15,19)_(18,23)_(24,25)_(3,4)_(5,9)_(7,10)_(11,13)_(12,16)_(14,17)_(15,20)_(19,21)_(22,24)_(5,8)_(6,7)_(9,12)_(10,11)_(13,14)_(15,16)_(17,20)_(18,19)_(21,23)_(24,25)_(5,6)_(7,8)_(9,10)_(11,12)_(13,15)_(14,16)_(17,18)_(19,20)_(21,22)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(22,23) }
        else if (n == 30) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(22,23)_(24,25)_(26,27)_(28,29)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(13,15)_(14,16)_(18,20)_(19,21)_(22,24)_(23,25)_(26,28)_(27,29)_(0,4)_(1,5)_(2,6)_(3,7)_(8,14)_(9,17)_(10,16)_(12,20)_(13,19)_(15,21)_(22,26)_(23,27)_(24,28)_(25,29)_(0,8)_(1,9)_(2,14)_(3,17)_(4,10)_(5,11)_(6,16)_(12,26)_(13,23)_(15,27)_(18,24)_(19,25)_(20,28)_(21,29)_(1,13)_(2,12)_(3,15)_(4,18)_(5,19)_(6,20)_(7,21)_(8,22)_(9,23)_(10,24)_(11,25)_(14,26)_(16,28)_(17,27)_(0,4)_(2,8)_(3,13)_(5,9)_(6,22)_(7,23)_(10,12)_(11,15)_(14,18)_(16,26)_(17,19)_(20,24)_(21,27)_(25,29)_(0,2)_(1,14)_(3,5)_(4,8)_(9,13)_(11,17)_(12,18)_(15,28)_(16,20)_(21,25)_(24,26)_(27,29)_(2,4)_(5,9)_(6,14)_(7,13)_(8,10)_(15,23)_(16,22)_(19,21)_(20,24)_(25,27)_(6,8)_(7,11)_(10,14)_(12,16)_(13,17)_(15,19)_(18,22)_(21,23)_(4,6)_(7,9)_(8,10)_(11,13)_(12,14)_(15,17)_(16,18)_(19,21)_(20,22)_(23,25)_(1,8)_(3,18)_(5,20)_(7,22)_(9,24)_(10,12)_(11,26)_(13,15)_(14,16)_(17,19)_(21,28)_(1,2)_(3,10)_(5,12)_(7,14)_(9,16)_(11,18)_(13,20)_(15,22)_(17,24)_(19,26)_(27,28)_(2,4)_(3,6)_(5,8)_(7,10)_(9,12)_(11,14)_(13,16)_(15,18)_(17,20)_(19,22)_(21,24)_(23,26)_(25,27)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18)_(19,20)_(21,22)_(23,24)_(25,26) }
        else if (n == 31) { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(22,23)_(24,25)_(26,27)_(28,29)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,14)_(13,15)_(16,18)_(17,19)_(20,22)_(21,23)_(24,26)_(25,27)_(28,30)_(0,4)_(1,5)_(2,6)_(3,7)_(8,12)_(9,13)_(10,14)_(11,15)_(16,20)_(17,21)_(18,22)_(19,23)_(24,28)_(25,29)_(26,30)_(0,8)_(1,9)_(2,10)_(3,11)_(4,12)_(5,13)_(6,14)_(7,15)_(16,24)_(17,25)_(18,26)_(19,27)_(20,28)_(21,29)_(22,30)_(0,16)_(1,8)_(2,4)_(3,12)_(5,10)_(6,9)_(7,14)_(11,13)_(17,24)_(18,20)_(19,28)_(21,26)_(22,25)_(23,30)_(27,29)_(1,2)_(3,5)_(4,8)_(6,22)_(7,11)_(9,25)_(10,12)_(13,14)_(17,18)_(19,21)_(20,24)_(23,27)_(26,28)_(29,30)_(1,17)_(2,18)_(3,19)_(4,20)_(5,10)_(7,23)_(8,24)_(11,27)_(12,28)_(13,29)_(14,30)_(21,26)_(3,17)_(4,16)_(5,21)_(6,18)_(7,9)_(8,20)_(10,26)_(11,23)_(13,25)_(14,28)_(15,27)_(22,24)_(1,4)_(3,8)_(5,16)_(7,17)_(9,21)_(10,22)_(11,19)_(12,20)_(14,24)_(15,26)_(23,28)_(27,30)_(2,5)_(7,8)_(9,18)_(11,17)_(12,16)_(13,22)_(14,20)_(15,19)_(23,24)_(26,29)_(2,4)_(6,12)_(9,16)_(10,11)_(13,17)_(14,18)_(15,22)_(19,25)_(20,21)_(27,29)_(5,6)_(8,12)_(9,10)_(11,13)_(14,16)_(15,17)_(18,20)_(19,23)_(21,22)_(25,26)_(3,5)_(6,7)_(8,9)_(10,12)_(11,14)_(13,16)_(15,18)_(17,20)_(19,21)_(22,23)_(24,25)_(26,28)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18)_(19,20)_(21,22)_(23,24)_(25,26)_(27,28) }
        else { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(14,15)_(16,17)_(18,19)_(20,21)_(22,23)_(24,25)_(26,27)_(28,29)_(30,31)_(0,2)_(1,3)_(4,6)_(5,7)_(8,10)_(9,11)_(12,14)_(13,15)_(16,18)_(17,19)_(20,22)_(21,23)_(24,26)_(25,27)_(28,30)_(29,31)_(0,4)_(1,5)_(2,6)_(3,7)_(8,12)_(9,13)_(10,14)_(11,15)_(16,20)_(17,21)_(18,22)_(19,23)_(24,28)_(25,29)_(26,30)_(27,31)_(0,8)_(1,9)_(2,10)_(3,11)_(4,12)_(5,13)_(6,14)_(7,15)_(16,24)_(17,25)_(18,26)_(19,27)_(20,28)_(21,29)_(22,30)_(23,31)_(0,16)_(1,8)_(2,4)_(3,12)_(5,10)_(6,9)_(7,14)_(11,13)_(15,31)_(17,24)_(18,20)_(19,28)_(21,26)_(22,25)_(23,30)_(27,29)_(1,2)_(3,5)_(4,8)_(6,22)_(7,11)_(9,25)_(10,12)_(13,14)_(17,18)_(19,21)_(20,24)_(23,27)_(26,28)_(29,30)_(1,17)_(2,18)_(3,19)_(4,20)_(5,10)_(7,23)_(8,24)_(11,27)_(12,28)_(13,29)_(14,30)_(21,26)_(3,17)_(4,16)_(5,21)_(6,18)_(7,9)_(8,20)_(10,26)_(11,23)_(13,25)_(14,28)_(15,27)_(22,24)_(1,4)_(3,8)_(5,16)_(7,17)_(9,21)_(10,22)_(11,19)_(12,20)_(14,24)_(15,26)_(23,28)_(27,30)_(2,5)_(7,8)_(9,18)_(11,17)_(12,16)_(13,22)_(14,20)_(15,19)_(23,24)_(26,29)_(2,4)_(6,12)_(9,16)_(10,11)_(13,17)_(14,18)_(15,22)_(19,25)_(20,21)_(27,29)_(5,6)_(8,12)_(9,10)_(11,13)_(14,16)_(15,17)_(18,20)_(19,23)_(21,22)_(25,26)_(3,5)_(6,7)_(8,9)_(10,12)_(11,14)_(13,16)_(15,18)_(17,20)_(19,21)_(22,23)_(24,25)_(26,28)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(13,14)_(15,16)_(17,18)_(19,20)_(21,22)_(23,24)_(25,26)_(27,28) }
    }
#endif
    #undef _
}


template<typename val4>
inline __device__ void _pps_memcpy(int n, val4* src, val4* dst) {
#if 0
    for (int i = 0; i < n; i++)
        dst[i] = src[i];
#elif 0
    int m = n / 4;
    float4* src4 = (float4*)src;
    float4* dst4 = (float4*)dst;
    for (int i = 0; i < m; i++) {
        dst4[i] = src4[i];
    }
    for (int i = 4*m; i < n; i++) {
        dst[i] = src[i];
    }
#else
    int m = (n + 3) / 4;
    float4* src4 = (float4*)src;
    float4* dst4 = (float4*)dst;
    for (int i = 0; i < m; i++) {
        dst4[i] = src4[i];
    }
#endif
}

template <PerPixelSortType SORT_TYPE>
__global__ void sort_per_pixel_kernel(
    _ARGS_sort_per_pixel_kernel
) {
    unsigned i = threadIdx.y + blockIdx.y * blockDim.y;
    unsigned j = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= img_height || j >= img_width)
        return;
    unsigned idx = i * img_width + j;
    unsigned tidx = threadIdx.y * blockDim.x + threadIdx.x;

    int n = num_intersects[idx];
    if (n <= 1)
        return;
    int32_t* indices_g = &indices_[idx*MAX_SORTED_SPLATS];
    float* depths_g = &depths_[idx*MAX_SORTED_SPLATS];

    __shared__ int32_t indices_s[MAX_SORTED_SPLATS*N_THREADS_PPS];
    __shared__ float depths_s[MAX_SORTED_SPLATS*N_THREADS_PPS];

    int32_t* indices = &indices_s[tidx*MAX_SORTED_SPLATS];
    float* depths = &depths_s[tidx*MAX_SORTED_SPLATS];
    _pps_memcpy<int32_t>(n, indices_g, indices);
    _pps_memcpy<float>(n, depths_g, depths);

    // int32_t* indices = indices_g;
    // float* depths = depths_g;

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
    case PerPixelSortType::NetworkSort:
        sort_per_pixel_network(n, indices, depths);
        break;
    }

    _pps_memcpy<int32_t>(n, indices, indices_g);
    _pps_memcpy<float>(n, depths, depths_g);
}



template<CameraType CAMERA_TYPE>
__global__ void rasterize_simple_sorted_forward_kernel(
    _ARGS_rasterize_simple_sorted_forward_kernel
) {
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= img_height || j >= img_width)
        return;
    int pix_id = i * img_width + j;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y)) {
            out_img[pix_id] = {0.0f, 0.0f, 0.0f};
            out_alpha[pix_id] = 0.0f;
            return;
        }
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

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
    }

    float3 final_color;
    final_color.x = pix_out.x + T * background.x;
    final_color.y = pix_out.y + T * background.y;
    final_color.z = pix_out.z + T * background.z;
    out_img[pix_id] = final_color;
    out_alpha[pix_id] = 1.0f - T;
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
    const float3 __restrict__ background,
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
        v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local)) / glm::vec2(fx, fy);
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



template <DepthMode DEPTH_MODE, CameraType CAMERA_TYPE>
__global__ void rasterize_depth_sorted_forward_kernel(
    _ARGS_rasterize_depth_sorted_forward_kernel
) {
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= img_height || j >= img_width)
        return;
    int pix_id = i * img_width + j;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y)) {
            out_depth[pix_id] = 0.0f;
            out_visibility[pix_id] = { 1.0f, 0.0f };
            return;
        }
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    // list of indices
    const int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];

    // current visibility left to render
    float T = 1.f;
    float interp = 1.0f;
    // index of most recent gaussian to write to this thread's pixel
    int cur_idx = 0;

    // rasterize
    float output_depth = 0.0f;

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


template <DepthMode DEPTH_MODE>
__global__ void rasterize_depth_sorted_backward_kernel(
    _ARGS_rasterize_depth_sorted_backward_kernel
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
        v_position_xy_abs_local = glm::abs(glm::vec2(v_position_local)) / glm::vec2(fx, fy);
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




__global__ void rasterize_sorted_forward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const float depth_reg_pairwise_factor,
    const int32_t* __restrict__ sorted_indices_,
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
    float* __restrict__ out_alpha,
    float3* __restrict__ out_img,
    float2* __restrict__ out_depth,
    float3* __restrict__ out_normal,
    float* __restrict__ out_reg_depth
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

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);

    float T = 1.f;  // current/total visibility
    float3 normal_out = {0.f, 0.f, 0.f};  // sum of normals
    float3 pix_out = {0.f, 0.f, 0.f};  // output radiance
    float vis_sum = 0.f;  // output alpha
    float depth_sum = 0.f;  // output depth
    float depth_squared_sum = 0.f;  // for L2 depth regularizer
    const float depth_ref = inside ? depth_ref_im[pix_id] : 0.f;
    float reg_depth_p = 0.f, reg_depth_i = 0.f;  // output depth regularizer

    // rasterize
    for (int cur_idx = 0; cur_idx < MAX_SORTED_SPLATS; cur_idx++) {
        int g_id = sorted_indices[cur_idx];
        if (g_id == SORTED_INDEX_INF)
            break;

        const glm::vec3 pos = *(glm::vec3*)&positions[g_id];
        const float opac = opacities[g_id];
        const glm::vec3 color_0 = *(glm::vec3*)&colors[g_id];
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

        glm::vec3 color;
        if (dim_ch > 0) {
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
    }

    if (inside) {
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


__global__ void rasterize_sorted_backward_kernel(
    const dim3 img_size,
    const float4 intrins,
    const unsigned ch_degree_r,
    const unsigned ch_degree_r_to_use,
    const unsigned ch_degree_phi,
    const unsigned ch_degree_phi_to_use,
    const float depth_reg_pairwise_factor,
    const int* __restrict__ num_intersects,
    const int32_t* __restrict__ sorted_indices_,
    const float3* __restrict__ positions,
    const float3* __restrict__ axes_u,
    const float3* __restrict__ axes_v,
    const float3* __restrict__ colors,
    const float3* __restrict__ ch_coeffs,
    const float* __restrict__ opacities,
    // const float3 __restrict__ background,
    const float* __restrict__ depth_ref_im,
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

    const int dim_ch = ch_degree_r * (2*ch_degree_phi+1);
    assert(dim_ch <= MAX_CH_FLOAT3);

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
    
    // rasterize
    const int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];
    for (int cur_idx = n-1; cur_idx >= 0; cur_idx--) {
        int g_id = sorted_indices[cur_idx];

        const glm::vec3 pos = *(glm::vec3*)&positions[g_id];
        const float opac = opacities[g_id];
        const glm::vec3 color_0 = *(glm::vec3*)&colors[g_id];
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
        glm::vec3 v_color_local = {0.f, 0.f, 0.f};
        float v_opacity_local = 0.f;
        glm::vec3 v_ch_coeff_local[MAX_CH_FLOAT3];
        for (int i = 0; i < dim_ch; i++)
            v_ch_coeff_local[i] = {0.f, 0.f, 0.f};
        float v_ch_coeff_abs_local = 0.f;
        //initialize everything to 0, only set if the lane is valid

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
            uv, opac,
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
        float* v_ch_coeffs_ptr = (float*)(v_ch_coeffs);
        for (int i = 0; i < dim_ch; i++) {
            atomicAdd(v_ch_coeffs_ptr + 3*dim_ch*g_id + 3*i + 0, v_ch_coeff_local[i].x);
            atomicAdd(v_ch_coeffs_ptr + 3*dim_ch*g_id + 3*i + 1, v_ch_coeff_local[i].y);
            atomicAdd(v_ch_coeffs_ptr + 3*dim_ch*g_id + 3*i + 2, v_ch_coeff_local[i].z);
        }
        // atomicAdd(v_ch_coeffs_abs + g, v_ch_coeff_abs_local);

        atomicAdd(v_opacities + g_id, v_opacity_local);
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
__global__ void rasterize_simplified_sorted_forward_kernel(
    _ARGS_rasterize_simplified_sorted_forward_kernel
) {
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= img_height || j >= img_width)
        return;
    int pix_id = i * img_width + j;

    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y)) {
            out_alpha[pix_id] = 0.0f;
            out_img[pix_id] = {0.0f, 0.0f, 0.0f};
            out_depth[pix_id] = {0.0f, 0.0f};
            out_normal[pix_id] = {0.0f, 0.0f, 0.0f};
            out_depth_reg[pix_id] = 0.0f;
            return;
        }
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    // list of indices
    const int32_t* sorted_indices = &sorted_indices_[pix_id*MAX_SORTED_SPLATS];

    float T = 1.f;  // current/total visibility
    float3 normal_out = {0.f, 0.f, 0.f};  // sum of normals
    float3 pix_out = {0.f, 0.f, 0.f};  // output radiance
    float vis_sum = 0.f;  // output alpha
    float depth_sum = 0.f;  // output depth
    float depth_squared_sum = 0.f;  // for L2 depth regularizer
    float reg_depth_p = 0.f;  // output depth regularizer

    // rasterize
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

        // color
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
    }

    out_alpha[pix_id] = 1.0f - T;
    out_img[pix_id] = pix_out;
    out_depth[pix_id] = { depth_sum, depth_squared_sum };
    out_normal[pix_id] = normal_out;
    out_depth_reg[pix_id] = reg_depth_p;
}


template<CameraType CAMERA_TYPE>
__global__ void rasterize_simplified_sorted_backward_kernel(
    _ARGS_rasterize_simplified_sorted_backward_kernel
) {
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= img_height || j >= img_width)
        return;
    int pix_id = i * img_width + j;

    int n = num_intersects[pix_id];
    if (n == 0) return;
    
    const float fx = intrins.x;
    const float fy = intrins.y;
    const float cx = intrins.z;
    const float cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y))
            return;
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

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

    const float vis_sum_final = 1.0f - T_final;
    const float depth_sum_final = output_depth[pix_id].x;
    const float depth_squared_sum_final = output_depth[pix_id].y;
    float vis_sum = vis_sum_final;

    float3 buffer = {0.f, 0.f, 0.f};
    float2 buffer_depth = {0.f, 0.f};  // depth, depth^2
    float3 buffer_normal = {0.f, 0.f, 0.f};
    float buffer_depth_reg = 0.f;

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
        //initialize everything to 0, only set if the lane is valid

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
            uv, opac,
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




template __global__ void rasterize_indices_kernel<CameraType::Undistorted>(
    _ARGS_rasterize_indices_kernel
);
template __global__ void rasterize_indices_kernel<CameraType::GenericDistorted>(
    _ARGS_rasterize_indices_kernel
);

template __global__ void sort_per_pixel_kernel<PerPixelSortType::InsertionSort>(
    _ARGS_sort_per_pixel_kernel
);
template __global__ void sort_per_pixel_kernel<PerPixelSortType::QuickSort>(
    _ARGS_sort_per_pixel_kernel
);
template __global__ void sort_per_pixel_kernel<PerPixelSortType::HeapSort>(
    _ARGS_sort_per_pixel_kernel
);
template __global__ void sort_per_pixel_kernel<PerPixelSortType::RandomizedQuickSort>(
    _ARGS_sort_per_pixel_kernel
);
template __global__ void sort_per_pixel_kernel<PerPixelSortType::NetworkSort>(
    _ARGS_sort_per_pixel_kernel
);

template __global__ void rasterize_simple_sorted_forward_kernel<CameraType::Undistorted>(
    _ARGS_rasterize_simple_sorted_forward_kernel
);
template __global__ void rasterize_simple_sorted_forward_kernel<CameraType::GenericDistorted>(
    _ARGS_rasterize_simple_sorted_forward_kernel
);

template __global__ void rasterize_depth_sorted_forward_kernel
<DepthMode::Mean, CameraType::Undistorted>(
    _ARGS_rasterize_depth_sorted_forward_kernel
);
template __global__ void rasterize_depth_sorted_forward_kernel
<DepthMode::Mean, CameraType::GenericDistorted>(
    _ARGS_rasterize_depth_sorted_forward_kernel
);
template __global__ void rasterize_depth_sorted_forward_kernel
<DepthMode::Median, CameraType::Undistorted>(
    _ARGS_rasterize_depth_sorted_forward_kernel
);
template __global__ void rasterize_depth_sorted_forward_kernel
<DepthMode::Median, CameraType::GenericDistorted>(
    _ARGS_rasterize_depth_sorted_forward_kernel
);

template __global__ void rasterize_depth_sorted_backward_kernel<DepthMode::Mean>(
    _ARGS_rasterize_depth_sorted_backward_kernel
);
template __global__ void rasterize_depth_sorted_backward_kernel<DepthMode::Median>(
    _ARGS_rasterize_depth_sorted_backward_kernel
);

template __global__ void rasterize_simplified_sorted_forward_kernel<CameraType::Undistorted>(
    _ARGS_rasterize_simplified_sorted_forward_kernel
);

template __global__ void rasterize_simplified_sorted_forward_kernel<CameraType::GenericDistorted>(
    _ARGS_rasterize_simplified_sorted_forward_kernel
);

template __global__ void rasterize_simplified_sorted_backward_kernel<CameraType::Undistorted>(
    _ARGS_rasterize_simplified_sorted_backward_kernel
);

template __global__ void rasterize_simplified_sorted_backward_kernel<CameraType::GenericDistorted>(
    _ARGS_rasterize_simplified_sorted_backward_kernel
);
