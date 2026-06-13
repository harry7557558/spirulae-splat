#include <cuda_runtime.h>
#include <cstdint>

#include <Common.cuh>

#ifndef NO_TORCH
#define NO_TORCH
#endif


#include <cooperative_groups.h>
namespace cg = cooperative_groups;



template<typename SplatPrimitive, CameraModelType camera_model, int VALUE_BITS = 32>
__global__ void projection_packed_mask_kernel(
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,  // [N, ...]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    bool *__restrict__ intersection_mask,  // [C, N]
    // SH VALUE-quant. The mask kernel discards splat_screen.rgb but still
    // calls project() which evaluates SH internally; passing the codec args
    // keeps the SH eval well-defined when fp32 features_sh is stale.
    const uint8_t* __restrict__ sh_value_packed = nullptr,
    const float2* __restrict__ sh_value_bounds  = nullptr,
    const uint32_t num_sh_buffer = 0,
    const int64_t  sh_bounds_stride = 0
) {
    // parallelize over C * N.
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= C * N) {
        return;
    }
    const uint32_t cid = (idx / N) % C; // camera id
    const uint32_t gid = idx % N; // gaussian id

    // Load camera
    viewmats += cid * 16;
    float4 intrin = intrins[cid];
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    ProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height,
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(cid);

    // Load splat
    // TODO: verify that SH is not loaded after compiler optimization
    typename SplatPrimitive::World splat_world;
    splat_world.load(splats_world, gid);

    // Projection
    float sorting_depth;
    float4 aabb;
    float radius;
    typename SplatPrimitive::Screen splat_screen;
    if constexpr (VALUE_BITS == 32) {
        splat_world.template project<camera_model, 32>(
            cam, splat_screen, aabb, sorting_depth, radius);
    } else {
        const int64_t sh_base = (int64_t)3 * (int64_t)num_sh_buffer * gid;
        const int64_t stride = (sh_bounds_stride > 0)
            ? sh_bounds_stride
            : (int64_t)256 * 3 * (int64_t)num_sh_buffer;
        splat_world.template project<camera_model, VALUE_BITS>(
            cam, splat_screen, aabb, sorting_depth, radius,
            const_cast<uint8_t*>(sh_value_packed),
            const_cast<float2*>(sh_value_bounds),
            sh_base, stride);
    }

    // Save results
    aabb.x = fminf(fmaxf(aabb.x, 0.0f), image_width-1.0f);
    aabb.y = fminf(fmaxf(aabb.y, 0.0f), image_height-1.0f);
    aabb.z = fminf(fmaxf(aabb.z, 0.0f), image_width-1.0f);
    aabb.w = fminf(fmaxf(aabb.w, 0.0f), image_height-1.0f);
    intersection_mask[idx] = (aabb.z - aabb.x > 1e-3f && aabb.w - aabb.y > 1e-3f);
}


template<typename SplatPrimitive, CameraModelType camera_model, int VALUE_BITS = 32>
__global__ void projection_packed_fwd_kernel(
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,  // [N, ...]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const int64_t* __restrict__ intersection_mask_scan,  // [C, N], inclusive scan
    // outputs
    int32_t *__restrict__ camera_ids,    // [nnz]
    int32_t *__restrict__ gaussian_ids,  // [nnz]
    float4 *__restrict__ aabbs,         // [nnz, 4]
    float *__restrict__ sorting_depths,  // [nnz]
    float *__restrict__ radii,  // [N]
    typename SplatPrimitive::ScreenBuffer splats_screen,  // [nnz, ...]
    const uint8_t* __restrict__ sh_value_packed = nullptr,
    const float2* __restrict__ sh_value_bounds  = nullptr,
    const uint32_t num_sh_buffer = 0,
    const int64_t  sh_bounds_stride = 0
) {
    // parallelize over C * N.
    uint32_t idx = cg::this_grid().thread_rank();
    if (idx >= C * N) {
        return;
    }

    int64_t out_idx = (idx == 0 ? 0 : intersection_mask_scan[idx-1]);
    int64_t out_idx_1 = intersection_mask_scan[idx];
    if (out_idx_1 == out_idx)
        return;

    const uint32_t cid = (idx / N) % C; // camera id
    const uint32_t gid = idx % N; // gaussian id

    // Load camera
    viewmats += cid * 16;
    float4 intrin = intrins[cid];
    float3x3 R = {
        viewmats[0], viewmats[1], viewmats[2],  // 1st row
        viewmats[4], viewmats[5], viewmats[6],  // 2nd row
        viewmats[8], viewmats[9], viewmats[10],  // 3rd row
    };
    float3 t = { viewmats[3], viewmats[7], viewmats[11] };
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    ProjCamera cam = {
        R, t, fx, fy, cx, cy,
        image_width, image_height,
       
    };
    cam.dist_coeffs = dist_coeffs_buffer.load(cid);

    // Load splat
    typename SplatPrimitive::World splat_world;
    splat_world.load(splats_world, gid);

    // Projection
    float sorting_depth;
    float4 aabb;
    float radius = 0.0f;
    typename SplatPrimitive::Screen splat_screen;
    if constexpr (VALUE_BITS == 32) {
        splat_world.template project<camera_model, 32>(
            cam, splat_screen, aabb, sorting_depth, radius);
    } else {
        const int64_t sh_base = (int64_t)3 * (int64_t)num_sh_buffer * gid;
        const int64_t stride = (sh_bounds_stride > 0)
            ? sh_bounds_stride
            : (int64_t)256 * 3 * (int64_t)num_sh_buffer;
        splat_world.template project<camera_model, VALUE_BITS>(
            cam, splat_screen, aabb, sorting_depth, radius,
            const_cast<uint8_t*>(sh_value_packed),
            const_cast<float2*>(sh_value_bounds),
            sh_base, stride);
    }

    // Save results
    camera_ids[out_idx] = (int32_t)cid;
    gaussian_ids[out_idx] = (int32_t)gid;
    aabb.x = fminf(fmaxf(aabb.x, 0.0f), image_width-1.0f);
    aabb.y = fminf(fmaxf(aabb.y, 0.0f), image_height-1.0f);
    aabb.z = fminf(fmaxf(aabb.z, 0.0f), image_width-1.0f);
    aabb.w = fminf(fmaxf(aabb.w, 0.0f), image_height-1.0f);
    aabbs[out_idx] = aabb;
    sorting_depths[out_idx] = sorting_depth;
    splat_screen.store(splats_screen, out_idx);
    atomicMax(&radii[gid], radius);
}


template<typename SplatPrimitive, CameraModelType camera_model>
void projection_packed_mask_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,  // [N, ...]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    // outputs
    bool *__restrict__ intersection_mask,  // [C, N]
    const uint8_t* __restrict__ sh_value_packed,
    const float2* __restrict__ sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_bounds_stride
) {
    constexpr uint block = 128;
    #define _LAUNCH(VB) \
        projection_packed_mask_kernel<SplatPrimitive, camera_model, VB> \
        <<<_CEIL_DIV(C*N, block), block, 0, stream>>>( \
            C, N, \
            splats_world, viewmats, intrins, dist_coeffs_buffer, \
            image_width, image_height, \
            intersection_mask, \
            sh_value_packed, sh_value_bounds, num_sh_buffer, sh_bounds_stride)
    if      (sh_value_bits == 8)  { _LAUNCH(8); }
    else if (sh_value_bits == 16) { _LAUNCH(16); }
    else                          { _LAUNCH(32); }
    #undef _LAUNCH
}

template<typename SplatPrimitive, CameraModelType camera_model>
void projection_packed_fwd_kernel_wrapper(
    cudaStream_t stream,
    const uint32_t C,
    const uint32_t N,
    typename SplatPrimitive::WorldBuffer splats_world,  // [N, ...]
    const float *__restrict__ viewmats, // [C, 4, 4]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer,
    const uint32_t image_width,
    const uint32_t image_height,
    const int64_t* __restrict__ intersection_mask_scan,  // [C, N], inclusive scan
    // outputs
    int32_t *__restrict__ camera_ids,    // [nnz]
    int32_t *__restrict__ gaussian_ids,  // [nnz]
    float4 *__restrict__ aabbs,         // [nnz, 4]
    float *__restrict__ sorting_depths,         // [nnz]
    float *__restrict__ radii,  // [N]
    typename SplatPrimitive::ScreenBuffer splats_screen,  // [nnz, ...]
    const uint8_t* __restrict__ sh_value_packed,
    const float2* __restrict__ sh_value_bounds,
    const uint32_t num_sh_buffer,
    const int sh_value_bits,
    const int64_t sh_bounds_stride
) {
    constexpr uint block = 128;
    #define _LAUNCH(VB) \
        projection_packed_fwd_kernel<SplatPrimitive, camera_model, VB> \
        <<<_CEIL_DIV(C*N, block), block, 0, stream>>>( \
            C, N, \
            splats_world, viewmats, intrins, dist_coeffs_buffer, \
            image_width, image_height, \
            intersection_mask_scan, \
            camera_ids, gaussian_ids, aabbs, sorting_depths, radii, splats_screen, \
            sh_value_packed, sh_value_bounds, num_sh_buffer, sh_bounds_stride)
    if      (sh_value_bits == 8)  { _LAUNCH(8); }
    else if (sh_value_bits == 16) { _LAUNCH(16); }
    else                          { _LAUNCH(32); }
    #undef _LAUNCH
}
