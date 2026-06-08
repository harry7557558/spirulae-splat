#include "Densify.cuh"

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangDensify {
#include "generated/set_namespace.cuh"
#include "generated/densify.cuh"
}
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
// #include "generated/projection_utils.cuh"
#include "generated/primitive_3dgs.cuh"
}
#endif

#include "Common.cuh"

#include <cub/cub.cuh>


// ================
// Helper: DeviceTensor3D<T> -> TensorView<U, 4>
// ================

template<typename U, typename T>
static TensorView<U, 4> _dt3d_to_tv4(const DeviceTensor3D<T>& dt) {
    static_assert(sizeof(T) % sizeof(U) == 0, "T must be a multiple of U");
    constexpr int C = sizeof(T) / sizeof(U);
    TensorView<U, 4> v;
    v.data = (U*)dt.data_ptr();
    v.shape[0] = dt.template size<0>();
    v.shape[1] = dt.template size<1>();
    v.shape[2] = dt.template size<2>();
    v.shape[3] = C;
    long HW = v.shape[1] * v.shape[2];
    v.strides[0] = HW * C;
    v.strides[1] = v.shape[2] * C;
    v.strides[2] = C;
    v.strides[3] = 1;
    return v;
}

// Build a TensorView<float,4> for [B,H,W,1] from a DeviceTensor3D<float>
static TensorView<float, 4> _dt3d_float_to_tv4_1ch(const DeviceTensor3D<float>& dt) {
    TensorView<float, 4> v;
    v.data = (float*)dt.data_ptr();
    v.shape[0] = dt.template size<0>();
    v.shape[1] = dt.template size<1>();
    v.shape[2] = dt.template size<2>();
    v.shape[3] = 1;
    long HW = v.shape[1] * v.shape[2];
    v.strides[0] = HW;
    v.strides[1] = v.shape[2];
    v.strides[2] = 1;
    v.strides[3] = 1;
    return v;
}


// ================
// Quantile / Median Normalization
// ================

template<bool invert_quantile>
cudaError_t batch_quantile_masked_radix_select(
    const float* d_x,
    int B,
    int N,
    float q,
    float* d_out,
    uint32_t* temp,
    cudaStream_t stream
);

// Internal helper: pool-allocated quantile computation
static void quantile_of_abs_of_finite_elements_internal(
    const float* inputs_ptr,
    int B,
    int N,
    float q,
    bool return_reciprocal,
    float* outputs_ptr
) {
    if (B == 0)
        return;
    float* temp_ptr = DevicePool::global().acquire<float>(
        "densify_quantile_temp", 1024 * B);

    (return_reciprocal ? batch_quantile_masked_radix_select<true> :
        batch_quantile_masked_radix_select<false>
    )(
        inputs_ptr,
        B, N, q,
        outputs_ptr,
        (uint32_t*)temp_ptr,
        (cudaStream_t)0
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void quantile_of_abs_of_finite_elements_tensor(
    DeviceVector<float> inputs,
    float q,
    bool return_reciprocal,
    DeviceVector<float> outputs
) {
    int64_t total = inputs.size();
    int64_t B = outputs.size();
    if (B == 0) return;
    int N = (int)(total / B);

    quantile_of_abs_of_finite_elements_internal(
        inputs.data_ptr(), B, N, q, return_reciprocal,
        outputs.data_ptr()
    );
}

__global__ void multiply_by_inverse_median_kernel(
    int B,
    int N,
    float* __restrict__ data,
    const float* __restrict__ quantiles
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= B*N)
        return;
    data[idx] *= quantiles[idx / N];
}

/*[AutoHeaderGeneratorExport]*/
void normalize_by_median_inplace_tensor(
    DeviceVector<float> data
) {
    int64_t total = data.size();
    // treat as 1-batch
    int B = 1;
    int N = (int)total;

    // pool-allocate inv_median [B]
    float* inv_median = DevicePool::global().acquire<float>(
        "densify_inv_median", B);

    quantile_of_abs_of_finite_elements_internal(
        data.data_ptr(), B, N, 0.5f, true, inv_median);

    multiply_by_inverse_median_kernel<<<_LAUNCH_ARGS_1D(B*N, 256)>>>(
        B, N,
        data.data_ptr(),
        inv_median
    );
}



// ================
// Indexing
// ================

__global__ void index_kernel(
    size_t numel,
    size_t stride,
    int32_t* __restrict__ indices,
    float* __restrict__ src,
    float* __restrict__ dst
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numel)
        return;

    size_t bidx = indices[idx / stride] * stride + idx % stride;
    dst[idx] = src[bidx];
}

/*[AutoHeaderGeneratorExport]*/
void inplace_index_tensor(
    DeviceVector<int32_t> indices,
    DeviceVector<float> src,
    DeviceVector<float> dst
) {
    if (indices.size() == 0)
        return;
    int64_t dst_numel = dst.size();
    int64_t stride = dst_numel / indices.size();
    index_kernel<<<_LAUNCH_ARGS_1D(dst_numel, 256)>>>(
        dst_numel,
        stride,
        indices.data_ptr(),
        src.data_ptr(),
        dst.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Scatter Reduce
// ================

__global__ void scatter_add_kernel(
    size_t numel,
    size_t stride,
    int32_t* __restrict__ indices,
    float* __restrict__ src,
    float* __restrict__ dst
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numel)
        return;

    size_t bidx = indices[idx / stride] * stride + idx % stride;
    float x = src[idx];
    if (x != 0.0f && isfinite(x))
        atomicAdd(&dst[bidx], src[idx]);
}

__global__ void scatter_max_kernel(
    int numel,
    int stride,
    int32_t* __restrict__ indices,
    float* __restrict__ src,
    float* __restrict__ dst
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numel)
        return;

    size_t bidx = indices[idx / stride] * stride + idx % stride;
    float x = src[idx];
    if (isfinite(x))
        atomicMax(&dst[bidx], src[idx]);
}

/*[AutoHeaderGeneratorExport]*/
void inplace_scatter_add_tensor(
    DeviceVector<int32_t> indices,
    DeviceVector<float> src,
    DeviceVector<float> dst
) {
    if (indices.size() == 0)
        return;
    int64_t src_numel = src.size();
    int64_t stride = src_numel / indices.size();
    scatter_add_kernel<<<_LAUNCH_ARGS_1D(src_numel, 256)>>>(
        src_numel,
        stride,
        indices.data_ptr(),
        src.data_ptr(),
        dst.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void inplace_scatter_max_tensor(
    DeviceVector<int32_t> indices,
    DeviceVector<float> src,
    DeviceVector<float> dst
) {
    if (indices.size() == 0)
        return;
    int64_t src_numel = src.size();
    int64_t stride = src_numel / indices.size();
    scatter_max_kernel<<<_LAUNCH_ARGS_1D(src_numel, 256)>>>(
        src_numel,
        stride,
        indices.data_ptr(),
        src.data_ptr(),
        dst.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Multinomial sample
// ================

__global__ void compute_efraimidis_spirakis_weight_kernel(
    int64_t numel,
    int stride,
    uint32_t seed,
    const float* weights,
    const bool* mask,
    float* out_weights
) {
    int64_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numel)
        return;

    float u = (float)hash_uint3(seed, blockIdx.x, threadIdx.x) * exp2f(-32.0f);
    float w = weights[idx * stride];
    w = w / __logf(fmaxf(1.0f - u, 1e-30f));  // larger w -> smaller (more negative) value
    if (mask != nullptr && !mask[idx])
        w = 0.0f;
    out_weights[idx] = w;
}

__global__ void iota_kernel(
    int64_t numel,
    int32_t* buffer
) {
    int64_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numel)
        return;
    buffer[idx] = idx;
}

// Internal helper for weighted sampling, returns pool-allocated indices
static int32_t* weighted_sample_without_replacement_internal(
    int64_t numel,
    float* weights_ptr,
    int64_t weights_numel,
    bool* masks_ptr,
    uint32_t num_sample,
    uint32_t seed
) {
    // int stride = (int)(weights_numel / numel);
    int stride = 2;  // above is incorrect during warmup

    float* sorting_values = DevicePool::global().acquire<float>(
        "densify_wswr_sorting_values", numel);

    compute_efraimidis_spirakis_weight_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        numel,
        stride,
        seed,
        weights_ptr,
        masks_ptr,
        sorting_values
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    int32_t* out_idx = DevicePool::global().acquire<int32_t>(
        "densify_wswr_out_idx", num_sample);

    float* d_keys_in = sorting_values;  // reuse
    float* d_keys_out = DevicePool::global().acquire<float>(
        "densify_wswr_keys_out", numel);

    int32_t* d_indices_in = DevicePool::global().acquire<int32_t>(
        "densify_wswr_indices_in", numel);
    int32_t* d_indices_out = DevicePool::global().acquire<int32_t>(
        "densify_wswr_indices_out", numel);

    iota_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        numel,
        d_indices_in
    );

    void* d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;

    cub::DeviceRadixSort::SortPairs(
        d_temp_storage,
        temp_storage_bytes,
        d_keys_in,
        d_keys_out,
        d_indices_in,
        d_indices_out,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    d_temp_storage = DeviceScratch::global().acquire(temp_storage_bytes);

    cub::DeviceRadixSort::SortPairs(
        d_temp_storage,
        temp_storage_bytes,
        d_keys_in,
        d_keys_out,
        d_indices_in,
        d_indices_out,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    cudaMemcpy(
        out_idx,
        d_indices_out,
        sizeof(int) * num_sample,
        cudaMemcpyDeviceToDevice
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_idx;
}

/*[AutoHeaderGeneratorExport]*/
void weighted_sample_without_replacement_tensor(
    int64_t numel,
    DeviceVector<float> weights,
    bool* masks_ptr,
    uint32_t num_sample,
    uint32_t seed,
    DeviceVector<int32_t> out_idx
) {
    if (numel == -1)
        numel = weights.size();

    int32_t* result = weighted_sample_without_replacement_internal(
        numel,
        weights.data_ptr(),
        weights.size(),
        masks_ptr,
        num_sample,
        seed
    );

    cudaMemcpy(
        out_idx.data_ptr(),
        result,
        sizeof(int32_t) * num_sample,
        cudaMemcpyDeviceToDevice
    );
}


// ================
// Covariance-Based Scale Initialization
// ================

__global__ void cov_scale_init_kernel(
    int64_t num_points,
    int32_t num_cameras,
    const float3* __restrict__ points,  // [N, 3]
    const bool* __restrict__ is_fisheye,  // [C]; TODO: equisolid
    const int2* __restrict__ sizes,  // [C, 2]
    const float4 *__restrict__ intrins,  // [C, 4], fx, fy, cx, cy
    const float4 *__restrict__ viewmats,  // [C, 4, 4]
    const CameraDistortionCoeffsBuffer dist_coeffs_buffer, // [C]
    float* __restrict__ log_scales  // [N]
) {
    int64_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_points)
        return;

    float3 p_world = points[idx];
    float max_log_scale = __logf(1e-30f);

    for (int32_t i = 0; i < num_cameras; ++i) {
        float4 intrin = intrins[i];
        int width = sizes[i].x, height = sizes[i].y;
        CameraDistortionCoeffs dist_coeffs = dist_coeffs_buffer.load(i);

        float4 p_wh = {p_world.x, p_world.y, p_world.z, 1.0f};
        float3 p_view = {
            dot(viewmats[4*i], p_wh),
            dot(viewmats[4*i+1], p_wh),
            dot(viewmats[4*i+2], p_wh),
        };

        bool valid = false;
        constexpr float eps = 1e-6f;
        float3x3 cov3d = {eps, 0, 0, 0, eps, 0, 0, 0, eps};
        float2x2 cov2d;
        float2 mean2d;
        if (is_fisheye[i]) {
            valid = SlangProjectionUtils::fisheye_proj_3dgs_nav(
                p_view, cov3d, intrin, dist_coeffs, &cov2d, &mean2d
            );
        }
        else {
            valid = SlangProjectionUtils::persp_proj_3dgs_nav(
                p_view, cov3d, intrin, dist_coeffs, width, height, &cov2d, &mean2d
            );
        }

        #pragma unroll
        for (int i = 0; i < 2; ++i) {
            cov2d[i].x = __fmul_rn(cov2d[i].x, 1.0f/eps);
            cov2d[i].y = __fmul_rn(cov2d[i].y, 1.0f/eps);
        }

        float det = cov2d[0].x * cov2d[1].y - cov2d[0].y * cov2d[1].x;
        if (valid && det > 0.0f) {
            float sc = 0.5f * __logf((float)(width * height) / det);
            max_log_scale = fmaxf(max_log_scale, sc);
        }
    }

    log_scales[idx] = max_log_scale;
}

/*[AutoHeaderGeneratorExport]*/
void cov_scale_init_tensor(
    DeviceVector<float3> points,  // [N, 3]
    DeviceVector<bool> is_fisheye,  // [C], bool
    DeviceVector<int2> sizes,  // [C, 2], int32
    DeviceVector<float4> intrins,  // [C, 4]
    DeviceVector<float4> viewmats,  // [C, 4, 4] as 4*C float4 elements
    TorchTensorView dist_coeffs, // [C]
    DeviceVector<float> log_scales  // [N, 1] output
) {
    int64_t N = points.size();
    int64_t C = intrins.size();

    cov_scale_init_kernel<<<_LAUNCH_ARGS_1D(N, 256)>>>(
        N, C,
        points.data_ptr(),
        is_fisheye.data_ptr(),
        sizes.data_ptr(),
        intrins.data_ptr(),
        viewmats.data_ptr(),
        dist_coeffs,
        log_scales.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Densification Parameter Update
// ================

__global__ void densify_clip_scale_kernel(
    long num_splats,
    float max_scale2d,
    float clip_hardness,
    float max_scale3d,
    const float* __restrict__ radii,
    float3* __restrict__ log_scales,
    float* __restrict__ logit_opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    if (isfinite(max_scale2d)) {
        float oversize_factor = radii[idx] / max_scale2d;
        if (isfinite(clip_hardness))
            oversize_factor = fminf(oversize_factor, clip_hardness);
        if (oversize_factor > 1.0f) {
            oversize_factor = __logf(oversize_factor);
            log_scales[idx] = log_scales[idx] - make_float3(oversize_factor);
            // this encourages being relocated to, may cause unintended effects in non-MCMC
            if (logit_opacs != nullptr)
                logit_opacs[idx] = fminf(logit_opacs[idx] + 3.0f * oversize_factor, fmaxf(logit_opacs[idx], 5.0f));
        }
    }

    if (isfinite(max_scale3d)) {
        max_scale3d = __logf(max_scale3d);
        log_scales[idx].x = fminf(log_scales[idx].x, max_scale3d);
        log_scales[idx].y = fminf(log_scales[idx].y, max_scale3d);
        log_scales[idx].z = fminf(log_scales[idx].z, max_scale3d);
        // don't touch opacity
    }
}

/*[AutoHeaderGeneratorExport]*/
void densify_clip_scale_tensor(
    int64_t num_splats,
    DeviceVector<float> radii,
    DeviceVector<float3> log_scales,
    float* logit_opacs_ptr,
    float max_scale2d,
    float clip_hardness,
    float max_scale3d
) {
    densify_clip_scale_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, max_scale2d, clip_hardness, max_scale3d,
        radii.data_ptr(),
        log_scales.data_ptr(),
        logit_opacs_ptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// https://www.shadertoy.com/view/4djSRW
__device__ __forceinline__ float hash14(float4 p4) {
    p4 = p4 * float4{.1031, .1030, .0973, .1099};
    p4.x = fmodf(p4.x, 1.0f); p4.y = fmodf(p4.y, 1.0f); p4.z = fmodf(p4.z, 1.0f); p4.w = fmodf(p4.w, 1.0f);
    p4 = p4 + dot(p4, float4{p4.w, p4.z, p4.x, p4.y} + 33.33f);
    return fmodf((p4.x + p4.y) * (p4.z + p4.w), 1.0f);
}

__global__ void densify_update_weight_kernel(
    long num_splats,
    int score_mode,
    const float* __restrict__ radii,  // [N]
    const float3* __restrict__ scales,  // [N, 3], optional
    const float* __restrict__ opacs,  // [N], optional
    const float* __restrict__ accum_weight_scalar,  // [1]
    const float* __restrict__ accum_weight,  // [N]
    float2* __restrict__ accum_buffer  // [N, 2]
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    if (!(radii[idx] > 0))
        return;

    float weight = fabsf(accum_weight[idx]);
    if (opacs)
        weight *= sigmoid(opacs[idx]);
    if (accum_weight_scalar != nullptr)
        weight *= accum_weight_scalar[0];
    if (weight == 0.0f)
        return;

    float2 accum = accum_buffer[idx];
    if (score_mode == (int)DensifyScoreMode::Max) {
        accum.x = fmaxf(accum.x, weight);
        accum.y = fmaxf(accum.y, 1.0f);
    } else if (score_mode == (int)DensifyScoreMode::Mean) {
        accum.x = (accum.x * accum.y + weight) / (accum.y + 1.0f);
        accum.y += 1.0f;
    } else if (score_mode == (int)DensifyScoreMode::Median) {
        float rand = hash14(1e2f * float4{weight, accum.x, accum.y, accum.y + 1.5f});
        if (accum.y == 0.0f) {
            accum.x = weight * exp2f(rand - 0.5f);
            accum.y = 1.0f;
        } else if (weight != 0.0f) {
            float s = weight > accum.x ? 1.0f : -1.0f;
            // accum.x *= exp2f(s / sqrtf(accum.y + 1.0f));
            accum.x *= exp2f(s * rand);
            accum.y += 1.0f;
        };
    } else if (score_mode == (int)DensifyScoreMode::Geom) {
        accum.x = __expf((__logf(fmaxf(accum.x, 1e-30f)) * accum.y +  __logf(fmaxf(weight, 1e-30f))) / (accum.y + 1.0f));
        accum.y += 1.0f;
    }
    accum_buffer[idx] = accum;
}


/*[AutoHeaderGeneratorExport]*/
void densify_update_weight(
    int64_t num_splats,
    DeviceVector<float> radii,
    float3* scales_ptr,
    float* opacs_ptr,
    DeviceVector<float> accum_weight,
    DeviceVector<float2> accum_buffer,
    int score_mode
) {
    densify_update_weight_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, score_mode,
        radii.data_ptr(),
        scales_ptr,
        opacs_ptr,
        nullptr,
        accum_weight.data_ptr(),
        accum_buffer.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Relocation
// ================

// Resetting a splat's quantized SH Adam state to "zero" means writing the
// (u_q, log_s_q) byte pair that decodes (with the CURRENT bounds slot for
// that cell) as closely as possible to (g1 = 0, g2 = 0). Naively writing 0
// bytes only matches that intent when the bounds happen to satisfy
// mm.x <= 0 <= mm.y AND mm.z == 0 -- otherwise decode(0,0,mm) lands at
// (mm.x, eps*expm1(mm.z)), which contaminates the new splat's first few
// optim steps. Saturation is acceptable: the next optim pass recomputes
// bounds and re-encodes.
//
// Bounds layout is one of:
//   bounds_per_splat=true  (FPBO): one float4 per 256 splats, covering all
//                                  3*num_sh cells per splat in the block.
//   bounds_per_splat=false (regular Optimizer.cu): one float4 per 256 cells,
//                                  so cells within a single splat may span
//                                  multiple bounds slots.
// Encode a single u or log_s value as the (BITS-bit) quantized nibble or byte
// that decodes back to `zero_val` within the [lo, hi] block bound. Used by
// densify to initialize dst splats' optim-state bytes to (g1 = 0, g2 = 0).
template<int BITS>
__device__ __forceinline__ uint8_t _quant_encode_zero_byte(float zero_val, float lo, float hi) {
    constexpr float qmax = (BITS == 8) ? 255.0f : 15.0f;  // BITS == 4 -> 15
    float range = fmaxf(hi - lo, 1e-30f);
    float qf = roundf(qmax * (zero_val - lo) / range);
    return (uint8_t)fminf(fmaxf(qf, 0.0f), qmax);
}

// SH VALUE-quant codec copy: decode src splat's SH cells against src bounds
// and encode into dst splat's cells against dst's CURRENT bounds. Cells
// outside dst's current quant range are clipped at byte 0 or kQMax.
//
// For the "add" case the dst slot's bound starts at (0, 0) -- everything
// clips to zero on the first densify step. The FPBO writeback's block-wide
// (min, max) reduction expands the bound on the next training step, so the
// child splat starts at SH=0 and adapts via Adam. Cleaner inheritance
// (atomic bound expansion + neighbor re-encoding) is left for later.
template<int BITS, bool BOUNDS_PER_SPLAT>
__device__ __forceinline__ void _copy_quant_sh_value_for_splat(
    uint8_t* __restrict__ packed,
    float2*  __restrict__ bounds,
    int64_t  src_splat,
    int64_t  dst_splat,
    int      num_sh,        // runtime SH count (degree-capped)
    int      num_sh_buffer  // buffer stride (model max)
) {
    using Codec = QuantizedTensor<BITS, 256>;
    int64_t src_base = src_splat * 3 * (int64_t)num_sh_buffer;
    int64_t dst_base = dst_splat * 3 * (int64_t)num_sh_buffer;

    float2 src_mm{0.f, 0.f}, dst_mm{0.f, 0.f};
    if constexpr (BOUNDS_PER_SPLAT) {
        src_mm = bounds[src_splat / 256];
        dst_mm = bounds[dst_splat / 256];
    }

    int64_t cells = (int64_t)3 * num_sh;
    #pragma unroll 1
    for (int64_t c = 0; c < cells; ++c) {
        int64_t src_cell = src_base + c;
        int64_t dst_cell = dst_base + c;
        if constexpr (!BOUNDS_PER_SPLAT) {
            src_mm = bounds[src_cell / 256];
            dst_mm = bounds[dst_cell / 256];
        }
        float v = Codec::decode_v(packed, src_cell, src_mm);
        Codec::encode_v(packed, dst_cell, v, dst_mm);
    }
}

template<int BITS, bool BOUNDS_PER_SPLAT>
__global__ void densify_copy_quant_sh_value_kernel(
    int64_t cur_num_splats,        // for default dst when dst_indices is null
    int64_t num_pairs,
    const int32_t* __restrict__ src_indices,
    const int32_t* __restrict__ dst_indices,  // nullptr -> dst = cur_num_splats + idx
    uint8_t* __restrict__ packed,
    float2*  __restrict__ bounds,
    int num_sh,
    int num_sh_buffer
) {
    int64_t pair = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (pair >= num_pairs) return;
    int32_t src = src_indices[pair];
    int32_t dst = (dst_indices == nullptr)
        ? (int32_t)(cur_num_splats + pair)
        : dst_indices[pair];
    _copy_quant_sh_value_for_splat<BITS, BOUNDS_PER_SPLAT>(
        packed, bounds, src, dst, num_sh, num_sh_buffer);
}

// Variant: src given by index_map[dst]; copy only when index_map[dst] != dst.
// Used by the MCMC relocate path (mcmc_update_relocation_kernel mirror).
template<int BITS, bool BOUNDS_PER_SPLAT>
__global__ void densify_copy_quant_sh_value_index_map_kernel(
    int64_t num_splats,
    const int32_t* __restrict__ index_map,   // [num_splats]
    uint8_t* __restrict__ packed,
    float2*  __restrict__ bounds,
    int num_sh,
    int num_sh_buffer
) {
    int64_t dst = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (dst >= num_splats) return;
    int32_t src = index_map[dst];
    if ((int64_t)src == dst) return;
    _copy_quant_sh_value_for_splat<BITS, BOUNDS_PER_SPLAT>(
        packed, bounds, src, dst, num_sh, num_sh_buffer);
}

static inline void _launch_densify_copy_quant_sh_value_index_map(
    int64_t num_splats,
    const int32_t* index_map,
    uint8_t* packed,
    float2*  bounds,
    int      num_sh,
    int      num_sh_buffer,
    int      bits,
    bool     bounds_per_splat
) {
    if (bits == 32 || num_splats == 0 || packed == nullptr || bounds == nullptr)
        return;
    constexpr int BLOCK = 256;
    int grid = (int)((num_splats + BLOCK - 1) / BLOCK);
    #define _LAUNCH(BB, BPS) \
        densify_copy_quant_sh_value_index_map_kernel<BB, BPS><<<grid, BLOCK>>>( \
            num_splats, index_map, packed, bounds, num_sh, num_sh_buffer)
    if      (bits == 8  &&  bounds_per_splat) { _LAUNCH(8,  true);  }
    else if (bits == 8  && !bounds_per_splat) { _LAUNCH(8,  false); }
    else if (bits == 16 &&  bounds_per_splat) { _LAUNCH(16, true);  }
    else if (bits == 16 && !bounds_per_splat) { _LAUNCH(16, false); }
    #undef _LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

// Host-side launcher: dispatches the codec copy based on bits + bounds layout.
// Safe no-op when bits == 32 (no value-quant) or num_pairs == 0.
static inline void _launch_densify_copy_quant_sh_value(
    int64_t cur_num_splats,
    int64_t num_pairs,
    const int32_t* src_indices,
    const int32_t* dst_indices,
    uint8_t* packed,
    float2*  bounds,
    int      num_sh,
    int      num_sh_buffer,
    int      bits,                  // 32 / 8 / 16
    bool     bounds_per_splat
) {
    if (bits == 32 || num_pairs == 0 || packed == nullptr || bounds == nullptr)
        return;
    constexpr int BLOCK = 256;
    int grid = (int)((num_pairs + BLOCK - 1) / BLOCK);
    #define _LAUNCH(BB, BPS) \
        densify_copy_quant_sh_value_kernel<BB, BPS><<<grid, BLOCK>>>( \
            cur_num_splats, num_pairs, src_indices, dst_indices, \
            packed, bounds, num_sh, num_sh_buffer)
    if      (bits == 8  &&  bounds_per_splat) { _LAUNCH(8,  true);  }
    else if (bits == 8  && !bounds_per_splat) { _LAUNCH(8,  false); }
    else if (bits == 16 &&  bounds_per_splat) { _LAUNCH(16, true);  }
    else if (bits == 16 && !bounds_per_splat) { _LAUNCH(16, false); }
    #undef _LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// Per-cell stores for the joint QuantizedAdamState<BITS, 256> codec.
//   BITS == 8: AoS layout, 2 bytes per cell -- byte[2k]=u_q, byte[2k+1]=log_s_q
//   BITS == 4: joint nibbles, 1 byte per cell -- low nibble = u_q, high = log_s_q
// Within densify both halves of every cell get zeroed for the dst splat, so
// these stores don't race with neighbors at the cell boundary (each thread
// owns its splat's full cell range).
template<int BITS>
__device__ __forceinline__ void _zero_quant_sh_store_cell(
    uint8_t* __restrict__ packed, int64_t cell, uint8_t u_q, uint8_t log_s_q
) {
    if constexpr (BITS == 8) {
        packed[cell * 2 + 0] = u_q;
        packed[cell * 2 + 1] = log_s_q;
    } else {  // BITS == 4
        packed[cell] = (uint8_t)((u_q & 0x0Fu) | ((log_s_q & 0x0Fu) << 4));
    }
}

// Non-SH 16-bit Adam-state codec mirrors QuantizedAdamState<16, 256> in
// Tensor.h: each cell is 4 bytes -- (uint16_t u_q)(uint16_t log_s_q). One
// per-splat-block float4 bound covers all PRIMS cells per splat for the
// attribute. Used by densify to initialize a relocated splat's optim-state
// bytes so they decode to (g1 = 0, g2 = 0) against the block's CURRENT bound;
// the next FPBO step's reduce expands the bound to cover the new value.
// 0 may fall outside the current bound (e.g. mm.x > 0); the encode clamps to
// the nearest endpoint then -- accepted small transient, identical to the
// SH-value densify-copy strategy.
__device__ __forceinline__ uint16_t _quant_encode_zero_word16(
    float zero_val, float lo, float hi
) {
    constexpr float qmax = 65535.0f;
    float range = fmaxf(hi - lo, 1e-30f);
    float qf = roundf(qmax * (zero_val - lo) / range);
    return (uint16_t)fminf(fmaxf(qf, 0.0f), qmax);
}

template<int PRIMS>
__device__ __forceinline__ void _zero_quant_non_sh_for_splat(
    uint8_t* __restrict__ packed,
    const float4* __restrict__ bounds,
    int64_t splat_idx
) {
    if (packed == nullptr || bounds == nullptr) return;
    float4 mm = bounds[splat_idx / 256];
    uint16_t u_q     = _quant_encode_zero_word16(0.0f, mm.x, mm.y);
    uint16_t log_s_q = _quant_encode_zero_word16(0.0f, mm.z, mm.w);
    int64_t base_cell = (int64_t)PRIMS * splat_idx;
    uint16_t* p = reinterpret_cast<uint16_t*>(packed);
    #pragma unroll
    for (int i = 0; i < PRIMS; ++i) {
        p[(base_cell + i) * 2 + 0] = u_q;
        p[(base_cell + i) * 2 + 1] = log_s_q;
    }
}

__device__ __forceinline__ void _zero_quant_non_sh_all(
    const NonShQuantState& non_sh, int64_t splat_idx
) {
    if (!non_sh.enabled) return;
    _zero_quant_non_sh_for_splat<3>(non_sh.means_packed,       non_sh.means_bounds,       splat_idx);
    _zero_quant_non_sh_for_splat<4>(non_sh.quats_packed,       non_sh.quats_bounds,       splat_idx);
    _zero_quant_non_sh_for_splat<3>(non_sh.scales_packed,      non_sh.scales_bounds,      splat_idx);
    _zero_quant_non_sh_for_splat<1>(non_sh.opacities_packed,   non_sh.opacities_bounds,   splat_idx);
    _zero_quant_non_sh_for_splat<3>(non_sh.features_dc_packed, non_sh.features_dc_bounds, splat_idx);
}


template<int BITS>
__device__ __forceinline__ void _zero_quant_sh_for_splat(
    uint8_t* __restrict__ packed,
    int64_t splat_idx,
    int num_sh,
    const float4* __restrict__ bounds,
    bool bounds_per_splat
) {
    int64_t cells_per_splat = (int64_t)num_sh * 3;
    int64_t base_cell = splat_idx * cells_per_splat;
    if (bounds_per_splat) {
        float4 mm = bounds[splat_idx / 256];
        uint8_t u_q     = _quant_encode_zero_byte<BITS>(0.0f, mm.x, mm.y);
        uint8_t log_s_q = _quant_encode_zero_byte<BITS>(0.0f, mm.z, mm.w);
        #pragma unroll 1
        for (int64_t c = 0; c < cells_per_splat; ++c) {
            _zero_quant_sh_store_cell<BITS>(packed, base_cell + c, u_q, log_s_q);
        }
    } else {
        #pragma unroll 1
        for (int64_t c = 0; c < cells_per_splat; ++c) {
            int64_t cell = base_cell + c;
            float4 mm = bounds[cell / 256];
            uint8_t u_q     = _quant_encode_zero_byte<BITS>(0.0f, mm.x, mm.y);
            uint8_t log_s_q = _quant_encode_zero_byte<BITS>(0.0f, mm.z, mm.w);
            _zero_quant_sh_store_cell<BITS>(packed, cell, u_q, log_s_q);
        }
    }
}

template<typename g_features_sh_t3>
__global__ void relocate_with_long_axis_split_kernel(
    int64_t cur_num_splats,
    int64_t num_new_splats,
    int32_t* __restrict__ src_indices,
    int32_t* __restrict__ dst_indices,
    float3*__restrict__ means, float3*__restrict__ g1_means, float3*__restrict__ g2_means,
    float4*__restrict__ quats, float4*__restrict__ g1_quats, float4*__restrict__ g2_quats,
    float3*__restrict__ scales, float3*__restrict__ g1_scales, float3*__restrict__ g2_scales,
    float*__restrict__ opacs, float*__restrict__ g1_opacs, float*__restrict__ g2_opacs,
    float3*__restrict__ features_dc, float3*__restrict__ g1_features_dc, float3*__restrict__ g2_features_dc,
    int num_sh, float3*__restrict__ features_sh, g_features_sh_t3*__restrict__ g1_features_sh, g_features_sh_t3*__restrict__ g2_features_sh,
    const float4* __restrict__ sh_quant_bounds,  // nullptr in non-quant mode
    bool sh_bounds_per_splat,
    NonShQuantState non_sh,
    float2*__restrict__ densify_accum_buffer,
    int32_t* __restrict__ bias_correction_steps
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_new_splats)
        return;

    int64_t idx_src = src_indices[idx];
    int64_t idx_dst = dst_indices == nullptr ? cur_num_splats + idx : dst_indices[idx];

    // geometry - long axis split
    float3 mean = means[idx_src], mean_delta;
    float3 scale = scales[idx_src], new_scale;
    float4 quat = quats[idx_src];
    float opac = opacs[idx_src], new_opac;
    SlangDensify::long_axis_split_3dgs(
        scale, opac, quat,
        &new_scale, &new_opac, &mean_delta
    );
    means[idx_src] = mean - mean_delta;
    means[idx_dst] = mean + mean_delta;
    scales[idx_src] = new_scale;
    scales[idx_dst] = new_scale;
    opacs[idx_src] = new_opac;
    opacs[idx_dst] = new_opac;
    quats[idx_dst] = quat;

    // appearance - copy
    features_dc[idx_dst] = features_dc[idx_src];
    // features_sh may be nullptr when SH-value quantization is on: the canonical
    // store lives in the packed buffer (copied below via _copy_quant_sh_value),
    // and the fp32 features_sh allocation is freed. Skip the fp32 copy then.
    if (features_sh) {
        for (int i = 0; i < num_sh; ++i)  // TODO: slow; more cache friendly way to do so?
            features_sh[num_sh*idx_dst+i] = features_sh[num_sh*idx_src+i];
    }

    // optimizer state - zero
#if 1
    // Non-SH Adam state: fp32 path zeros the g1_/g2_ buffers directly;
    // 16-bit quant path encodes (g1=0, g2=0) into the packed bytes against
    // the current per-splat-block bound (which may not include 0 -- the
    // encode clamps to the nearest endpoint and the next FPBO step's reduce
    // expands the bound).
    if (g1_means) {
        g1_means[idx_dst] = make_float3(0.0f);
        g2_means[idx_dst] = make_float3(0.0f);
        g1_quats[idx_dst] = make_float4(0.0f);
        g2_quats[idx_dst] = make_float4(0.0f);
        g1_scales[idx_dst] = make_float3(0.0f);
        g2_scales[idx_dst] = make_float3(0.0f);
        g1_opacs[idx_dst] = 0.0f;
        g2_opacs[idx_dst] = 0.0f;
        g1_features_dc[idx_dst] = make_float3(0.0f);
        g2_features_dc[idx_dst] = make_float3(0.0f);
    }
    _zero_quant_non_sh_all(non_sh, idx_dst);
    if constexpr (sizeof(g_features_sh_t3) == sizeof(short3)) {
        // 8-bit quant path: write encoded (u=0, log_s=0) bytes; g1 and g2
        // alias the same packed buffer, so one pass covers both.
        _zero_quant_sh_for_splat<8>(
            (uint8_t*)g1_features_sh, idx_dst, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else if constexpr (sizeof(g_features_sh_t3) == sizeof(uchar3)) {
        // 4-bit quant path: joint nibbles, 1 byte/cell. Same aliasing as 8-bit.
        _zero_quant_sh_for_splat<4>(
            (uint8_t*)g1_features_sh, idx_dst, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else {
        for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
            g1_features_sh[num_sh*idx_dst+i] = g_features_sh_t3{0,0,0};
            g2_features_sh[num_sh*idx_dst+i] = g_features_sh_t3{0,0,0};
        }
    }
    if (bias_correction_steps)
        bias_correction_steps[idx_dst] = 0;
    #if 0
    g1_means[idx_src] = make_float3(0.0f);
    g2_means[idx_src] = make_float3(0.0f);
    g1_quats[idx_src] = make_float4(0.0f);
    g2_quats[idx_src] = make_float4(0.0f);
    g1_scales[idx_src] = make_float3(0.0f);
    g2_scales[idx_src] = make_float3(0.0f);
    g1_opacs[idx_src] = 0.0f;
    g2_opacs[idx_src] = 0.0f;
    g1_features_dc[idx_src] = make_float3(0.0f);
    g2_features_dc[idx_src] = make_float3(0.0f);
    if constexpr (sizeof(g_features_sh_t3) == sizeof(short3)) {
        _zero_quant_sh_for_splat<8>(
            (uint8_t*)g1_features_sh, idx_src, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else if constexpr (sizeof(g_features_sh_t3) == sizeof(uchar3)) {
        _zero_quant_sh_for_splat<4>(
            (uint8_t*)g1_features_sh, idx_src, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else {
        for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
            g1_features_sh[num_sh*idx_src+i] = g_features_sh_t3{0,0,0};
            g2_features_sh[num_sh*idx_src+i] = g_features_sh_t3{0,0,0};
        }
    }
    if (bias_correction_steps)
        bias_correction_steps[idx_src] = 0;
    #endif
#else
    // to avoid messing up Adam bias correction
    constexpr float k = 1.0f;
    g1_means[idx_dst] = g1_means[idx_src] = k*g1_means[idx_src];
    g2_means[idx_dst] = g2_means[idx_src] = k*g2_means[idx_src];
    g1_quats[idx_dst] = g1_quats[idx_src] = k*g1_quats[idx_src];
    g2_quats[idx_dst] = g2_quats[idx_src] = k*g2_quats[idx_src];
    g1_scales[idx_dst] = g1_scales[idx_src] = k*g1_scales[idx_src];
    g2_scales[idx_dst] = g2_scales[idx_src] = k*g2_scales[idx_src];
    g1_opacs[idx_dst] = g1_opacs[idx_src] = k*g1_opacs[idx_src];
    g2_opacs[idx_dst] = g2_opacs[idx_src] = k*g2_opacs[idx_src];
    g1_features_dc[idx_dst] = g1_features_dc[idx_src] = k*g1_features_dc[idx_src];
    g2_features_dc[idx_dst] = g2_features_dc[idx_src] = k*g2_features_dc[idx_src];
    for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
        g1_features_sh[num_sh*idx_dst+i] = g1_features_sh[num_sh*idx_src+i] = k*g1_features_sh[num_sh*idx_src+i];
        g2_features_sh[num_sh*idx_dst+i] = g2_features_sh[num_sh*idx_src+i] = k*g2_features_sh[num_sh*idx_src+i];
    }
    if (bias_correction_steps)
        bias_correction_steps[idx_dst] = bias_correction_steps[idx_src];
#endif

    // will be cleared after densification anyway but doesn't hurt to do so
    if (densify_accum_buffer)
        densify_accum_buffer[idx_dst] = densify_accum_buffer[idx_src];
}


__global__ void compute_relocation_mask_kernel(
    int64_t num_splats,
    float min_opacity,
    const float3* __restrict__ means,
    const float4* __restrict__ quats,
    const float3* __restrict__ scales,
    const float* __restrict__ opacities,
    const float3* __restrict__ features_dc,
    bool* __restrict__ masks,
    int32_t* __restrict__ num_relocate,
    int32_t* __restrict__ relocate_indices
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    int num_relocate_delta = 0;

    if (idx < num_splats) {

        float3 mean = means[idx];
        float4 quat = quats[idx];
        float3 scale = scales[idx];
        float opac = opacities[idx];
        float3 feature_dc = features_dc[idx];

        bool is_low_opac = sigmoid(opac) <= min_opacity;

        bool is_finite = isfinite(
            dot(mean, mean) / dot(quat, quat) + dot(scale, feature_dc) * opac
        );
        bool relocate = is_low_opac || !is_finite;

        masks[idx] = !relocate;

        num_relocate_delta = (int)relocate;
    }

    // warpAtomicAdd(num_relocate, num_relocate_delta);

    // TODO: faster way to do so
    if (num_relocate_delta)
        relocate_indices[atomicAdd(num_relocate, num_relocate_delta)] = (int32_t)idx;
}


/*[AutoHeaderGeneratorExport]*/
void relocate_splats_with_long_axis_split_tensor(
    int64_t cur_num_splats,
    float min_opacity,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<float2> densify_accum_buffer,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    // SH-quant bounds buffer + layout flag used to encode (g1=0, g2=0) into
    // the dst splats' packed bytes. Null when sh_optim_bits == 32 (no quant).
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    // SH VALUE-quant codec copy params. When sh_value_bits != 32 we also do
    // a codec-aware src->dst copy of the SH coefficient bytes (decode src,
    // encode dst against current dst bounds; clipping is acceptable -- see
    // _copy_quant_sh_value_for_splat for the rationale).
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,             // 32 / 8 / 16
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,             // buffer stride for cell indexing
    // Non-SH Adam-state quant: when enabled, each relocated dst splat gets
    // its packed bytes set to codec-encoded zero against the current bound.
    NonShQuantState non_sh,
    uint32_t seed
) {
    int32_t* bias_correction_steps_ptr = bias_correction_steps.data_ptr();
    bool* mask = DevicePool::global().acquire<bool>(
        "densify_reloc_mask", cur_num_splats);
    int32_t* num_relocate_ptr = DevicePool::global().acquire<int32_t>(
        "densify_reloc_count", 1);
    cudaMemset(num_relocate_ptr, 0, sizeof(int32_t));

    int32_t* dst_indices = DevicePool::global().acquire<int32_t>(
        "densify_reloc_dst_indices", cur_num_splats);

    compute_relocation_mask_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        cur_num_splats,
        min_opacity,
        means.data_ptr(),
        quats.data_ptr(),
        scales.data_ptr(),
        opacs.data_ptr(),
        features_dc.data_ptr(),
        mask,
        num_relocate_ptr,
        dst_indices
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    int32_t num_relocate_host = 0;
    cudaMemcpy(&num_relocate_host, num_relocate_ptr, sizeof(int32_t), cudaMemcpyDeviceToHost);
    int64_t num_relocate = (int64_t)num_relocate_host;
    if (num_relocate == 0)
        return;

    int32_t* src_indices = weighted_sample_without_replacement_internal(
        cur_num_splats, (float*)densify_accum_buffer.data_ptr(),
        densify_accum_buffer.size() * 2, mask, num_relocate, seed);

    #define _DENSIFY_LAS_LAUNCH(T, bnd_ptr, bps) \
            relocate_with_long_axis_split_kernel<T><<<_LAUNCH_ARGS_1D(num_relocate, 256)>>>( \
                cur_num_splats, \
                num_relocate, \
                src_indices, \
                dst_indices, \
                means.data_ptr(), g1_means.data_ptr(), g2_means.data_ptr(), \
                quats.data_ptr(), g1_quats.data_ptr(), g2_quats.data_ptr(), \
                scales.data_ptr(), g1_scales.data_ptr(), g2_scales.data_ptr(), \
                opacs.data_ptr(), g1_opacs.data_ptr(), g2_opacs.data_ptr(), \
                features_dc.data_ptr(), g1_features_dc.data_ptr(), g2_features_dc.data_ptr(), \
                num_sh, features_sh.data_ptr(), (T*)g1_features_sh.data_ptr(), (T*)g2_features_sh.data_ptr(), \
                bnd_ptr, bps, \
                non_sh, \
                densify_accum_buffer.data_ptr(), \
                bias_correction_steps_ptr)
        if      (sh_optim_bits == 32) _DENSIFY_LAS_LAUNCH(float3, nullptr, false);
        else if (sh_optim_bits == 8)  _DENSIFY_LAS_LAUNCH(short3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
        else if (sh_optim_bits == 4)  _DENSIFY_LAS_LAUNCH(uchar3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
        #undef _DENSIFY_LAS_LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());

    _launch_densify_copy_quant_sh_value(
        cur_num_splats, num_relocate, src_indices, dst_indices,
        sh_value_packed.data_ptr(), sh_value_bounds.data_ptr(),
        num_sh, num_sh_buffer, sh_value_bits, sh_value_bounds_per_splat);
}

/*[AutoHeaderGeneratorExport]*/
void add_splats_with_long_axis_split_tensor(
    int64_t cur_num_splats,
    int64_t num_new_splats,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<float2> densify_accum_buffer,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,
    NonShQuantState non_sh,
    uint32_t seed
) {
    if (num_new_splats == 0)
        return;
    int32_t* bias_correction_steps_ptr = bias_correction_steps.data_ptr();

    int32_t* split_indices = weighted_sample_without_replacement_internal(
        cur_num_splats, (float*)densify_accum_buffer.data_ptr(),
        densify_accum_buffer.size() * 2, nullptr, num_new_splats, seed);

    #define _DENSIFY_LAS_ADD_LAUNCH(T, bnd_ptr, bps) \
        relocate_with_long_axis_split_kernel<T><<<_LAUNCH_ARGS_1D(num_new_splats, 256)>>>( \
            cur_num_splats, \
            num_new_splats, \
            split_indices, \
            nullptr, \
            means.data_ptr(), g1_means.data_ptr(), g2_means.data_ptr(), \
            quats.data_ptr(), g1_quats.data_ptr(), g2_quats.data_ptr(), \
            scales.data_ptr(), g1_scales.data_ptr(), g2_scales.data_ptr(), \
            opacs.data_ptr(), g1_opacs.data_ptr(), g2_opacs.data_ptr(), \
            features_dc.data_ptr(), g1_features_dc.data_ptr(), g2_features_dc.data_ptr(), \
            num_sh, features_sh.data_ptr(), (T*)g1_features_sh.data_ptr(), (T*)g2_features_sh.data_ptr(), \
            bnd_ptr, bps, \
            non_sh, \
            densify_accum_buffer.data_ptr(), \
            bias_correction_steps_ptr)
    if      (sh_optim_bits == 32) _DENSIFY_LAS_ADD_LAUNCH(float3, nullptr, false);
    else if (sh_optim_bits == 8)  _DENSIFY_LAS_ADD_LAUNCH(short3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
    else if (sh_optim_bits == 4)  _DENSIFY_LAS_ADD_LAUNCH(uchar3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
    #undef _DENSIFY_LAS_ADD_LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());

    _launch_densify_copy_quant_sh_value(
        cur_num_splats, num_new_splats, split_indices, /*dst_indices=*/nullptr,
        sh_value_packed.data_ptr(), sh_value_bounds.data_ptr(),
        num_sh, num_sh_buffer, sh_value_bits, sh_value_bounds_per_splat);
}


// ================
// MCMC Relocation
// ================

// largely from https://github.com/harry7557558/vksplat/blob/main/vksplat/slang/mcmc.slang

inline __device__ float log_factorial(float x) {
    // least squares for x = [0, 1... 50]
    // absolute error: max 4.6e-5, L1 9.3e-6
    // notebook: https://www.desmos.com/calculator/mik6cz7h5v
    static const float kA = 0.996763591291f;
    static const float kH = 2.73507778369f;
    static const float kB = 0.978531458767f;
    static const float kC2 = -0.0400131099548f;
    static const float kC1 = -0.596969060666f;
    static const float kD3 = 0.00251628691783f;
    static const float kD2 = 0.0959712063178f;
    static const float kD1 = 0.803207449544f;
    return (kA * __logf(x + kH) - kB +
        (kC2 * x + kC1) / (((kD3 * x + kD2) * x + kD1) * x + 1.0f)
    ) * x;
}

inline __device__ float binom(float n, float k) {
    return __expf(log_factorial(n) - log_factorial(k) - log_factorial(n-k));
}

inline __device__ void mcmc_relocation(float& opacity, float3& scale, int n_idx) {
    n_idx = min(n_idx, 50);  // log_factorial only fits to 50

    opacity = sigmoid(opacity);
    scale = {__expf(scale.x), __expf(scale.y), __expf(scale.z)};

    float new_opacity = 1.0f - powf(1.0f-opacity, 1.0f / n_idx);

    float denom_sum = 0.0f;
    for (int i = 1; i <= n_idx; ++i) {
        for (int k = 0; k <= (i - 1); ++k) {
            denom_sum += binom(i-1, k) *
                (cosf((float)M_PI*k) / sqrtf(k+1)) *  // (-1)^k / sqrt(k+1)
                powf(new_opacity, k+1);
        }
    }
    float coeff = (opacity / denom_sum);

    opacity = new_opacity;
    scale = coeff * scale;

    opacity = logit(opacity);
    scale = {__logf(scale.x), __logf(scale.y), __logf(scale.z)};
}


__global__ void mcmc_compute_relocation_probabilities_kernel(
    uint32_t num_splats,
    float min_opacity,
    const float* __restrict__ opacs,
    float* __restrict__ probs
) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_splats)
        return;

    float opac = sigmoid(opacs[tid]);

    if (opac <= min_opacity || !isfinite(opac))
        opac = 0.0f;

    probs[tid] = opac;
}


__global__ void mcmc_compute_relocation_index_map_kernel(
    float* __restrict__ sample_probs,  // [N]
    float* __restrict__ sample_probs_cumsum,  // [N]
    int32_t* __restrict__ index_map,  // [N]
    int32_t* __restrict__ n_idx_buffer,  // [N]
    uint32_t numel,
    uint32_t seed
) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= numel)
        return;

    float prob = sample_probs[tid];

    // not relocated, put same ID
    if (prob != 0.0f) {
        index_map[tid] = tid;
        return;
    }

    // binary search for random index
    // find the min idx such that sample_probs_cumsum[idx] >= rand
    float rand_max = sample_probs_cumsum[numel-1];
    float rand = rand_max * hash_uint3(seed, blockIdx.x, threadIdx.x) * exp2f(-32.0f);
    uint32_t lo = 0;
    uint32_t hi = numel - 1;
    while (lo < hi) {
        uint32_t mid = lo + ((hi - lo) >> 1);
        if (sample_probs_cumsum[mid] < rand)
            lo = mid + 1;
        else
            hi = mid;
    }
    index_map[tid] = lo;

    atomicAdd(&n_idx_buffer[lo], 1);
}


template<typename g_features_sh_t3>
__global__ void mcmc_compute_relocation_kernel(
    uint32_t num_splats,
    float min_opacity,
    int32_t* __restrict__ n_idx_buffer,  // [N]
    float3*__restrict__ means, float3*__restrict__ g1_means, float3*__restrict__ g2_means,
    float4*__restrict__ quats, float4*__restrict__ g1_quats, float4*__restrict__ g2_quats,
    float3*__restrict__ scales, float3*__restrict__ g1_scales, float3*__restrict__ g2_scales,
    float*__restrict__ opacs, float*__restrict__ g1_opacs, float*__restrict__ g2_opacs,
    float3*__restrict__ features_dc, float3*__restrict__ g1_features_dc, float3*__restrict__ g2_features_dc,
    int num_sh, float3*__restrict__ features_sh, g_features_sh_t3*__restrict__ g1_features_sh, g_features_sh_t3*__restrict__ g2_features_sh,
    const float4* __restrict__ sh_quant_bounds,
    bool sh_bounds_per_splat,
    NonShQuantState non_sh,
    int32_t* __restrict__ bias_correction_steps
) {
    uint32_t cur_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cur_idx >= num_splats)
        return;

    // only compute on indices being relocated to
    uint n_idx = n_idx_buffer[cur_idx] + 1;
    if (n_idx == 1)
        return;

    // compute relocation
    float3 scale = scales[cur_idx];
    float opac = opacs[cur_idx];
    mcmc_relocation(opac, scale, int(n_idx));
    scales[cur_idx] = scale;
    opacs[cur_idx] = opac;

    // set grad to zero (skip non-SH fp32 writes when those buffers are
    // freed by FPBO non-SH Adam-state quantization; see relocate kernel).
    if (g1_means) {
        g1_means[cur_idx] = make_float3(0.0f);
        g2_means[cur_idx] = make_float3(0.0f);
        g1_quats[cur_idx] = make_float4(0.0f);
        g2_quats[cur_idx] = make_float4(0.0f);
        g1_scales[cur_idx] = make_float3(0.0f);
        g2_scales[cur_idx] = make_float3(0.0f);
        g1_opacs[cur_idx] = 0.0f;
        g2_opacs[cur_idx] = 0.0f;
        g1_features_dc[cur_idx] = make_float3(0.0f);
        g2_features_dc[cur_idx] = make_float3(0.0f);
    }
    if constexpr (sizeof(g_features_sh_t3) == sizeof(short3)) {
        _zero_quant_sh_for_splat<8>(
            (uint8_t*)g1_features_sh, cur_idx, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else if constexpr (sizeof(g_features_sh_t3) == sizeof(uchar3)) {
        _zero_quant_sh_for_splat<4>(
            (uint8_t*)g1_features_sh, cur_idx, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else {
        for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
            g1_features_sh[num_sh*cur_idx+i] = g_features_sh_t3{0,0,0};
            g2_features_sh[num_sh*cur_idx+i] = g_features_sh_t3{0,0,0};
        }
    }
    _zero_quant_non_sh_all(non_sh, cur_idx);
    if (bias_correction_steps)
        bias_correction_steps[cur_idx] = 0;
}


__global__ void mcmc_update_relocation_kernel(
    uint32_t num_splats,
    int32_t* __restrict__ index_map,  // [N]
    float3*__restrict__ means,
    float4*__restrict__ quats,
    float3*__restrict__ scales,
    float*__restrict__ opacs,
    float3*__restrict__ features_dc,
    int num_sh, float3*__restrict__ features_sh
) {
    uint32_t id_dst = blockIdx.x * blockDim.x + threadIdx.x;
    if (id_dst >= num_splats)
        return;

    // only update on those relocated from
    uint32_t id_src = index_map[id_dst];
    if (id_src == id_dst)
        return;

    // copy components
    means[id_dst] = means[id_src];
    quats[id_dst] = quats[id_src];
    scales[id_dst] = scales[id_src];
    opacs[id_dst] = opacs[id_src];
    features_dc[id_dst] = features_dc[id_src];
    // features_sh may be null under SH-value quantization (canonical store is
    // the packed buffer; fp32 features_sh allocation is freed).
    if (features_sh) {
        for (int i = 0; i < num_sh; ++i)  // TODO: slow; more cache friendly way to do so?
            features_sh[num_sh*id_dst+i] = features_sh[num_sh*id_src+i];
    }
}


/*[AutoHeaderGeneratorExport]*/
void relocate_splats_mcmc_tensor(
    int64_t cur_num_splats,
    float min_opacity,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,
    NonShQuantState non_sh,
    uint32_t seed
) {
    int32_t* bias_correction_steps_ptr = bias_correction_steps.data_ptr();
    float* sample_probs = DevicePool::global().acquire<float>(
        "densify_mcmc_sample_probs", cur_num_splats);
    mcmc_compute_relocation_probabilities_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        cur_num_splats,
        min_opacity,
        opacs.data_ptr(),
        sample_probs
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // CUB inclusive sum for cumsum
    float* sample_probs_cumsum = DevicePool::global().acquire<float>(
        "densify_mcmc_sample_probs_cumsum", cur_num_splats);
    {
        size_t temp_bytes = 0;
        cub::DeviceScan::InclusiveSum(
            nullptr, temp_bytes,
            sample_probs, sample_probs_cumsum, (int)cur_num_splats);
        void* temp = DeviceScratch::global().acquire(temp_bytes);
        cub::DeviceScan::InclusiveSum(
            temp, temp_bytes,
            sample_probs, sample_probs_cumsum, (int)cur_num_splats);
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    int32_t* index_map = DevicePool::global().acquire<int32_t>(
        "densify_mcmc_index_map", cur_num_splats);
    int32_t* n_idx_buffer = DevicePool::global().acquire<int32_t>(
        "densify_mcmc_n_idx_buffer", cur_num_splats);
    cudaMemset(n_idx_buffer, 0, cur_num_splats * sizeof(int32_t));

    mcmc_compute_relocation_index_map_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        sample_probs,
        sample_probs_cumsum,
        index_map,
        n_idx_buffer,
        cur_num_splats,
        seed
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #define _DENSIFY_MCMC_RELOC_LAUNCH(T, bnd_ptr, bps) \
        mcmc_compute_relocation_kernel<T><<<_LAUNCH_ARGS_1D(cur_num_splats, 64)>>>( \
            cur_num_splats, \
            min_opacity, \
            n_idx_buffer, \
            means.data_ptr(), g1_means.data_ptr(), g2_means.data_ptr(), \
            quats.data_ptr(), g1_quats.data_ptr(), g2_quats.data_ptr(), \
            scales.data_ptr(), g1_scales.data_ptr(), g2_scales.data_ptr(), \
            opacs.data_ptr(), g1_opacs.data_ptr(), g2_opacs.data_ptr(), \
            features_dc.data_ptr(), g1_features_dc.data_ptr(), g2_features_dc.data_ptr(), \
            num_sh, features_sh.data_ptr(), (T*)g1_features_sh.data_ptr(), (T*)g2_features_sh.data_ptr(), \
            bnd_ptr, bps, \
            non_sh, \
            bias_correction_steps_ptr)
    if      (sh_optim_bits == 32) _DENSIFY_MCMC_RELOC_LAUNCH(float3, nullptr, false);
    else if (sh_optim_bits == 8)  _DENSIFY_MCMC_RELOC_LAUNCH(short3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
    else if (sh_optim_bits == 4)  _DENSIFY_MCMC_RELOC_LAUNCH(uchar3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
    #undef _DENSIFY_MCMC_RELOC_LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());

    mcmc_update_relocation_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 64)>>>(
        cur_num_splats,
        index_map,
        means.data_ptr(),
        quats.data_ptr(),
        scales.data_ptr(),
        opacs.data_ptr(),
        features_dc.data_ptr(),
        num_sh, features_sh.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    _launch_densify_copy_quant_sh_value_index_map(
        cur_num_splats, index_map,
        sh_value_packed.data_ptr(), sh_value_bounds.data_ptr(),
        num_sh, num_sh_buffer, sh_value_bits, sh_value_bounds_per_splat);
}



__global__ void mcmc_compute_add_index_map_kernel(
    float* __restrict__ sample_probs_cumsum,  // [N]
    int32_t* __restrict__ index_map,  // [dN]
    int32_t* __restrict__ n_idx_buffer,  // [N]
    uint32_t num_splats,
    uint32_t num_add,
    uint32_t seed
) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_add)
        return;

    // binary search for random index
    // find the min idx such that sample_probs_cumsum[idx] >= rand
    float rand_max = sample_probs_cumsum[num_splats-1];
    float rand = rand_max * hash_uint3(seed, blockIdx.x, threadIdx.x) * exp2f(-32.0f);
    uint32_t lo = 0;
    uint32_t hi = num_splats - 1;
    while (lo < hi) {
        uint32_t mid = lo + ((hi - lo) >> 1);
        if (sample_probs_cumsum[mid] <= rand)
            lo = mid + 1;
        else
            hi = mid;
    }
    index_map[tid] = lo;

    atomicAdd(&n_idx_buffer[lo], 1);
}


__global__ void mcmc_compute_add_kernel(
    uint32_t num_splats,
    float min_opacity,
    int32_t* __restrict__ n_idx_buffer,
    float3*__restrict__ scales,
    float*__restrict__ opacs
) {
    uint32_t cur_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (cur_idx >= num_splats)
        return;

    // only compute on indices being relocated to
    uint n_idx = n_idx_buffer[cur_idx] + 1;
    if (n_idx == 1)
        return;

    // compute relocation
    float3 scale = scales[cur_idx];
    float opac = opacs[cur_idx];
    mcmc_relocation(opac, scale, int(n_idx));
    scales[cur_idx] = scale;
    opacs[cur_idx] = opac;

    // don't change grad
}


template<typename g_features_sh_t3>
__global__ void mcmc_update_add_kernel(
    uint32_t num_splats,
    uint32_t num_add,
    int32_t* __restrict__ index_map,  // [dN]
    float3*__restrict__ means, float3*__restrict__ g1_means, float3*__restrict__ g2_means,
    float4*__restrict__ quats, float4*__restrict__ g1_quats, float4*__restrict__ g2_quats,
    float3*__restrict__ scales, float3*__restrict__ g1_scales, float3*__restrict__ g2_scales,
    float*__restrict__ opacs, float*__restrict__ g1_opacs, float*__restrict__ g2_opacs,
    float3*__restrict__ features_dc, float3*__restrict__ g1_features_dc, float3*__restrict__ g2_features_dc,
    int num_sh, float3*__restrict__ features_sh, g_features_sh_t3*__restrict__ g1_features_sh, g_features_sh_t3*__restrict__ g2_features_sh,
    const float4* __restrict__ sh_quant_bounds,
    bool sh_bounds_per_splat,
    NonShQuantState non_sh,
    int32_t* __restrict__ bias_correction_steps
) {
    uint32_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid >= num_add)
        return;

    uint32_t id_dst = tid + num_splats;
    uint32_t id_src = index_map[tid];

    // copy components
    means[id_dst] = means[id_src];
    quats[id_dst] = quats[id_src];
    scales[id_dst] = scales[id_src];
    opacs[id_dst] = opacs[id_src];
    features_dc[id_dst] = features_dc[id_src];
    // features_sh may be null under SH-value quantization (canonical store is
    // the packed buffer; fp32 features_sh allocation is freed).
    if (features_sh) {
        for (int i = 0; i < num_sh; ++i)  // TODO: slow; more cache friendly way to do so?
            features_sh[num_sh*id_dst+i] = features_sh[num_sh*id_src+i];
    }

    // set grad to zero (skip non-SH fp32 writes when freed by FPBO non-SH
    // Adam-state quantization).
    if (g1_means) {
        g1_means[id_dst] = make_float3(0.0f);
        g2_means[id_dst] = make_float3(0.0f);
        g1_quats[id_dst] = make_float4(0.0f);
        g2_quats[id_dst] = make_float4(0.0f);
        g1_scales[id_dst] = make_float3(0.0f);
        g2_scales[id_dst] = make_float3(0.0f);
        g1_opacs[id_dst] = 0.0f;
        g2_opacs[id_dst] = 0.0f;
        g1_features_dc[id_dst] = make_float3(0.0f);
        g2_features_dc[id_dst] = make_float3(0.0f);
    }
    if constexpr (sizeof(g_features_sh_t3) == sizeof(short3)) {
        _zero_quant_sh_for_splat<8>(
            (uint8_t*)g1_features_sh, id_dst, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else if constexpr (sizeof(g_features_sh_t3) == sizeof(uchar3)) {
        _zero_quant_sh_for_splat<4>(
            (uint8_t*)g1_features_sh, id_dst, num_sh,
            sh_quant_bounds, sh_bounds_per_splat);
    } else {
        for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
            g1_features_sh[num_sh*id_dst+i] = g_features_sh_t3{0,0,0};
            g2_features_sh[num_sh*id_dst+i] = g_features_sh_t3{0,0,0};
        }
    }
    _zero_quant_non_sh_all(non_sh, id_dst);
    if (bias_correction_steps)
        bias_correction_steps[id_dst] = 0;
}


/*[AutoHeaderGeneratorExport]*/
void add_splats_mcmc_tensor(
    int64_t cur_num_splats,
    int64_t num_add,
    float min_opacity,
    DeviceVector<float3> means, DeviceVector<float4> quats, DeviceVector<float3> scales, DeviceVector<float> opacs, DeviceVector<float3> features_dc, DeviceVector<float3> features_sh,
    DeviceVector<float3> g1_means, DeviceVector<float4> g1_quats, DeviceVector<float3> g1_scales, DeviceVector<float> g1_opacs, DeviceVector<float3> g1_features_dc, DeviceVector<float3> g1_features_sh,
    DeviceVector<float3> g2_means, DeviceVector<float4> g2_quats, DeviceVector<float3> g2_scales, DeviceVector<float> g2_opacs, DeviceVector<float3> g2_features_dc, DeviceVector<float3> g2_features_sh,
    DeviceVector<int32_t> bias_correction_steps,
    int sh_optim_bits,
    int num_sh,
    DeviceVector<float4> sh_quant_bounds,
    bool sh_bounds_per_splat,
    DeviceVector<uint8_t> sh_value_packed,
    DeviceVector<float2>  sh_value_bounds,
    int  sh_value_bits,
    bool sh_value_bounds_per_splat,
    int  num_sh_buffer,
    NonShQuantState non_sh,
    uint32_t seed
) {
    if (num_add == 0)
        return;
    int32_t* bias_correction_steps_ptr = bias_correction_steps.data_ptr();

    float* sample_probs = DevicePool::global().acquire<float>(
        "densify_mcmc_add_sample_probs", cur_num_splats);
    mcmc_compute_relocation_probabilities_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        cur_num_splats,
        min_opacity,
        opacs.data_ptr(),
        sample_probs
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // CUB inclusive sum for cumsum
    float* sample_probs_cumsum = DevicePool::global().acquire<float>(
        "densify_mcmc_add_sample_probs_cumsum", cur_num_splats);
    {
        size_t temp_bytes = 0;
        cub::DeviceScan::InclusiveSum(
            nullptr, temp_bytes,
            sample_probs, sample_probs_cumsum, (int)cur_num_splats);
        void* temp = DeviceScratch::global().acquire(temp_bytes);
        cub::DeviceScan::InclusiveSum(
            temp, temp_bytes,
            sample_probs, sample_probs_cumsum, (int)cur_num_splats);
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    int32_t* index_map = DevicePool::global().acquire<int32_t>(
        "densify_mcmc_add_index_map", num_add);
    int32_t* n_idx_buffer = DevicePool::global().acquire<int32_t>(
        "densify_mcmc_add_n_idx_buffer", cur_num_splats);
    cudaMemset(n_idx_buffer, 0, cur_num_splats * sizeof(int32_t));

    mcmc_compute_add_index_map_kernel<<<_LAUNCH_ARGS_1D(num_add, 256)>>>(
        sample_probs_cumsum,
        index_map,
        n_idx_buffer,
        cur_num_splats,
        num_add,
        seed
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    mcmc_compute_add_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 64)>>>(
        cur_num_splats,
        min_opacity,
        n_idx_buffer,
        scales.data_ptr(),
        opacs.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #define _DENSIFY_MCMC_ADD_LAUNCH(T, bnd_ptr, bps) \
        mcmc_update_add_kernel<T><<<_LAUNCH_ARGS_1D(num_add, 64)>>>( \
            cur_num_splats, \
            num_add, \
            index_map, \
            means.data_ptr(), g1_means.data_ptr(), g2_means.data_ptr(), \
            quats.data_ptr(), g1_quats.data_ptr(), g2_quats.data_ptr(), \
            scales.data_ptr(), g1_scales.data_ptr(), g2_scales.data_ptr(), \
            opacs.data_ptr(), g1_opacs.data_ptr(), g2_opacs.data_ptr(), \
            features_dc.data_ptr(), g1_features_dc.data_ptr(), g2_features_dc.data_ptr(), \
            num_sh, features_sh.data_ptr(), (T*)g1_features_sh.data_ptr(), (T*)g2_features_sh.data_ptr(), \
            bnd_ptr, bps, \
            non_sh, \
            bias_correction_steps_ptr)
    if      (sh_optim_bits == 32) _DENSIFY_MCMC_ADD_LAUNCH(float3, nullptr, false);
    else if (sh_optim_bits == 8)  _DENSIFY_MCMC_ADD_LAUNCH(short3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
    else if (sh_optim_bits == 4)  _DENSIFY_MCMC_ADD_LAUNCH(uchar3, sh_quant_bounds.data_ptr(), sh_bounds_per_splat);
    #undef _DENSIFY_MCMC_ADD_LAUNCH
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // mcmc_update_add: dst[i] = cur_num_splats + i, src[i] = index_map[i]
    _launch_densify_copy_quant_sh_value(
        cur_num_splats, num_add, index_map, /*dst_indices=*/nullptr,
        sh_value_packed.data_ptr(), sh_value_bounds.data_ptr(),
        num_sh, num_sh_buffer, sh_value_bits, sh_value_bounds_per_splat);
}


// ================
// MCMC Add Noise
// ================


__global__ void mcmc_add_noise_3dgs_kernel(
    long num_splats,
    float scaler,
    float3* __restrict__ means,
    const float3* __restrict__ log_scales,
    const float4* __restrict__ quats,
    const float* __restrict__ logit_opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::mcmc_add_noise_3dgs(
        scaler,
        &means[idx], log_scales[idx], quats[idx], logit_opacs[idx]
    );
}

__global__ void revised_add_noise_3dgs_kernel(
    long num_splats,
    float scaler,
    const float* __restrict__ radii,
    float3* __restrict__ means,
    const float3* __restrict__ log_scales,
    const float4* __restrict__ quats,
    const float* __restrict__ logit_opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::revised_add_noise_3dgs(
        scaler, radii[idx],
        &means[idx], log_scales[idx], quats[idx], logit_opacs[idx]
    );
}

// __global__ void mcmc_add_noise_triangle_kernel(
//     long num_splats,
//     float scaler, float min_opacity,
//     float3* __restrict__ means,
//     const float3* __restrict__ log_scales,
//     const float4* __restrict__ quats,
//     const float* __restrict__ opacs
// ) {
//     size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
//     if (idx >= num_splats)
//         return;

//     SlangDensify::mcmc_add_noise_triangle(
//         scaler, min_opacity,
//         &means[idx], log_scales[idx], quats[idx], opacs[idx]
//     );
// }


/*[AutoHeaderGeneratorExport]*/
void mcmc_add_noise_tensor(
    int64_t num_splats,
    float scaler,
    DeviceVector<float3> means,
    DeviceVector<float3> log_scales,
    DeviceVector<float4> quats,
    DeviceVector<float> opacs
) {
    mcmc_add_noise_3dgs_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, scaler,
        means.data_ptr(),
        log_scales.data_ptr(),
        quats.data_ptr(),
        opacs.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


/*[AutoHeaderGeneratorExport]*/
void revised_add_noise_tensor(
    int64_t num_splats,
    float scaler,
    DeviceVector<float> radii,
    DeviceVector<float3> means,
    DeviceVector<float3> log_scales,
    DeviceVector<float4> quats,
    DeviceVector<float> opacs
) {
    revised_add_noise_3dgs_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, scaler,
        radii.data_ptr(),
        means.data_ptr(),
        log_scales.data_ptr(),
        quats.data_ptr(),
        opacs.data_ptr()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Long axis split (https://arxiv.org/abs/2508.12313)
// ================

__global__ void long_axis_split_3dgs_kernel(
    long num_splats,
    const float3* __restrict__ log_scales,
    const float* __restrict__ logit_opacities,
    const float4* __restrict__ quats,
    float3* __restrict__ new_log_scales,
    float* __restrict__ new_logit_opacities,
    float3* __restrict__ mean_deltas
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::long_axis_split_3dgs(
        log_scales[idx], logit_opacities[idx], quats[idx],
        &new_log_scales[idx], &new_logit_opacities[idx], &mean_deltas[idx]
    );
}

/*[AutoHeaderGeneratorExport]*/
void long_axis_split_tensor(
    std::string primitive,
    DeviceVector<float3> log_scales,
    DeviceVector<float> logit_opacities,
    DeviceVector<float4> quats,
    DeviceVector<float3> new_log_scales,
    DeviceVector<float> new_logit_opacities,
    DeviceVector<float3> mean_deltas
) {
    const size_t num_splats = quats.size();

    if (num_splats == 0)
        return;

    if (primitive == "3dgs" || primitive == "mip" || primitive == "3dgut")
        long_axis_split_3dgs_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
            num_splats,
            log_scales.data_ptr(),
            logit_opacities.data_ptr(),
            quats.data_ptr(),
            new_log_scales.data_ptr(),
            new_logit_opacities.data_ptr(),
            mean_deltas.data_ptr()
        );
    else throw std::runtime_error("Unsupported primitive: " + primitive);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Image filter (https://arxiv.org/abs/2508.12313, https://arxiv.org/abs/2603.08661)
// ================

__constant__ float blur_filter_5x5[5][5] = {
    { 2.f/159.f,  4.f/159.f,  5.f/159.f,  4.f/159.f, 2.f/159.f },
    { 4.f/159.f,  9.f/159.f, 12.f/159.f,  9.f/159.f, 4.f/159.f },
    { 5.f/159.f, 12.f/159.f, 15.f/159.f, 12.f/159.f, 5.f/159.f },
    { 4.f/159.f,  9.f/159.f, 12.f/159.f,  9.f/159.f, 4.f/159.f },
    { 2.f/159.f,  4.f/159.f,  5.f/159.f,  4.f/159.f, 2.f/159.f },
};

__constant__ float laplacian_filter_3x3[3][3] = {
    { 1.f/6.f,  4.f/6.f, 1.f/6.f },
    { 4.f/6.f,-20.f/6.f, 4.f/6.f },
    { 1.f/6.f,  4.f/6.f, 1.f/6.f },
};

__constant__ float canny_filter_3x3[3][3] = {
    { -1.0f, 0.0f, 1.0f },
    { -2.0f, 0.0f, 2.0f },
    { -1.0f, 0.0f, 1.0f },
};

__global__ void laplacian_edge_filter_kernel(
    const TensorView<float, 4> img_in,
    TensorView<float, 4> img_out
) {
    constexpr int BLOCK = 32;
    constexpr int HALO = 1;

    const int32_t xid = blockIdx.x * BLOCK + threadIdx.x;
    const int32_t yid = blockIdx.y * BLOCK + threadIdx.y;
    const int32_t bid = blockIdx.z;
    const uint32_t H = img_in.shape[1], W = img_in.shape[2];

    // load pixels
    __shared__ float shared_pixels[BLOCK+2*HALO][BLOCK+2*HALO];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO)*(BLOCK+2*HALO); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO);
        int x = tid % (BLOCK+2*HALO);
        if (y < BLOCK+2*HALO) {
            int yi = min(max((int)(blockIdx.y * BLOCK) + y - HALO, 0), H-1);
            int xi = min(max((int)(blockIdx.x * BLOCK) + x - HALO, 0), W-1);
            shared_pixels[y][x] = dot(img_in.load3(bid, yi, xi), float3{0.299f, 0.587f, 0.114f});
        }
    }
    __syncthreads();

    // 3x3 Laplacian
    float total = 0.0f;
    for (int cy = -HALO; cy <= HALO; ++cy)
        for (int cx = -HALO; cx <= HALO; ++cx) {
            float conv_weight = laplacian_filter_3x3[cy+HALO][cx+HALO];
            int yi = threadIdx.y + cy;
            int xi = threadIdx.x + cx;
            total += conv_weight * shared_pixels[yi+HALO][xi+HALO];
        }
    total = fabsf(total);
    if (yid < H && xid < W)
        img_out.store1(bid, yid, xid, total);
}

/*[AutoHeaderGeneratorExport]*/
void laplacian_edge_filter_tensor(
    DeviceTensor3D<float3> img_in,
    DeviceTensor3D<float> img_out
) {
    int B  = img_in.template size<0>();
    int H  = img_in.template size<1>();
    int W  = img_in.template size<2>();

    laplacian_edge_filter_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        _dt3d_to_tv4<float>(img_in),
        _dt3d_float_to_tv4_1ch(img_out)
    );
}

__global__ void smoothed_laplacian_edge_filter_kernel(
    const TensorView<float, 4> img_in,
    TensorView<float, 4> img_out
) {
    constexpr int BLOCK = 32;
    constexpr int HALO = 3;
    constexpr int HALO1 = 1;

    const int32_t xid = blockIdx.x * BLOCK + threadIdx.x;
    const int32_t yid = blockIdx.y * BLOCK + threadIdx.y;
    const int32_t bid = blockIdx.z;
    const uint32_t H = img_in.shape[1], W = img_in.shape[2];

    // load pixels
    __shared__ float shared_pixels[BLOCK+2*HALO][BLOCK+2*HALO];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO)*(BLOCK+2*HALO); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO);
        int x = tid % (BLOCK+2*HALO);
        if (y < BLOCK+2*HALO) {
            int yi = min(max((int)(blockIdx.y * BLOCK) + y - HALO, 0), H-1);
            int xi = min(max((int)(blockIdx.x * BLOCK) + x - HALO, 0), W-1);
            shared_pixels[y][x] = dot(img_in.load3(bid, yi, xi), float3{0.299f, 0.587f, 0.114f});
        }
    }
    __syncthreads();

    // TODO: optimize into 7x7 horizontal + vertical convolution?

    // 5x5 blur
    __shared__ float shared_blurred[BLOCK+2*HALO1][BLOCK+2*HALO1];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO1)*(BLOCK+2*HALO1); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO1);
        int x = tid % (BLOCK+2*HALO1);
        if (y >= BLOCK+2*HALO1)
            continue;
        float total = 0.0f;
        for (int cy = -2; cy <= 2; ++cy)
            for (int cx = -2; cx <= 2; ++cx) {
                float conv_weight = blur_filter_5x5[cy+2][cx+2];
                int yi = y - HALO1 + cy;
                int xi = x - HALO1 + cx;
                total += conv_weight * shared_pixels[yi+HALO][xi+HALO];
            }
        shared_blurred[y][x] = total;
    }
    __syncthreads();

    // 3x3 Laplacian
    float total = 0.0f;
    for (int cy = -HALO1; cy <= HALO1; ++cy)
        for (int cx = -HALO1; cx <= HALO1; ++cx) {
            float conv_weight = laplacian_filter_3x3[cy+HALO1][cx+HALO1];
            int yi = threadIdx.y + cy;
            int xi = threadIdx.x + cx;
            total += conv_weight * shared_blurred[yi+HALO1][xi+HALO1];
        }
    total = fabsf(total);
    if (yid < H && xid < W)
        img_out.store1(bid, yid, xid, total);
}

/*[AutoHeaderGeneratorExport]*/
void smoothed_laplacian_edge_filter_tensor(
    DeviceTensor3D<float3> img_in,
    DeviceTensor3D<float> img_out
) {
    int B  = img_in.template size<0>();
    int H  = img_in.template size<1>();
    int W  = img_in.template size<2>();

    smoothed_laplacian_edge_filter_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        _dt3d_to_tv4<float>(img_in),
        _dt3d_float_to_tv4_1ch(img_out)
    );
}

__global__ void canny_edge_filter_kernel(
    TensorView<float, 4> img_in,
    const bool* __restrict__ mask_in,
    TensorView<float, 4> img_out
) {
    constexpr int BLOCK = 32;
    constexpr int HALO = 4;
    constexpr int HALO1 = 2;
    constexpr int HALO2 = 1;

    const int32_t xid = blockIdx.x * BLOCK + threadIdx.x;
    const int32_t yid = blockIdx.y * BLOCK + threadIdx.y;
    const int32_t bid = blockIdx.z;
    const uint32_t H = img_in.shape[1], W = img_in.shape[2];

    // load pixels
    __shared__ float shared_pixels[BLOCK+2*HALO][BLOCK+2*HALO];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO)*(BLOCK+2*HALO); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO);
        int x = tid % (BLOCK+2*HALO);
        if (y < BLOCK+2*HALO) {
            int yi = min(max((int)(blockIdx.y * BLOCK) + y - HALO, 0), H-1);
            int xi = min(max((int)(blockIdx.x * BLOCK) + x - HALO, 0), W-1);
            shared_pixels[y][x] = dot(img_in.load3(bid, yi, xi), float3{0.299f, 0.587f, 0.114f});
        }
    }
    __syncthreads();

    // 5x5 blur
    __shared__ float shared_blurred[BLOCK+2*HALO1][BLOCK+2*HALO1];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO1)*(BLOCK+2*HALO1); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO1);
        int x = tid % (BLOCK+2*HALO1);
        if (y >= BLOCK+2*HALO1)
            continue;
        float total = 0.0f;
        for (int cy = -2; cy <= 2; ++cy)
            for (int cx = -2; cx <= 2; ++cx) {
                float conv_weight = blur_filter_5x5[cy+2][cx+2];
                int yi = y - HALO1 + cy;
                int xi = x - HALO1 + cx;
                total += conv_weight * shared_pixels[yi+HALO][xi+HALO];
            }
        shared_blurred[y][x] = total;
    }
    __syncthreads();

    // 3x3 canny
    __shared__ float2 shared_filtered[BLOCK+2*HALO2][BLOCK+2*HALO2];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO2)*(BLOCK+2*HALO2); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO2);
        int x = tid % (BLOCK+2*HALO2);
        if (y >= BLOCK+2*HALO2)
            continue;
        float total1 = 0.0f, total2 = 0.0f;
        for (int cy = -1; cy <= 1; ++cy)
            for (int cx = -1; cx <= 1; ++cx) {
                float conv_weight_1 = canny_filter_3x3[cy+1][cx+1];
                float conv_weight_2 = canny_filter_3x3[cx+1][cy+1];
                int yi = y - HALO2 + cy;
                int xi = x - HALO2 + cx;
                float value = shared_blurred[yi+HALO1][xi+HALO1];
                total1 += conv_weight_1 * value;
                total2 += conv_weight_2 * value;
            }
        shared_filtered[y][x] = {total1, total2};
    }
    __syncthreads();

    // non-maximum suppression
    float2 total = shared_filtered[threadIdx.y+HALO2][threadIdx.x+HALO2];
    float mag = length(total);
    if (mag > 0.0f) {
        int dx = min(max((int)roundf(total.x / mag), -HALO2), HALO2);
        int dy = min(max((int)roundf(total.y / mag), -HALO2), HALO2);
        float total1 = length(shared_filtered[(int)threadIdx.y+dy+HALO2][(int)threadIdx.x+dx+HALO2]);
        float total2 = length(shared_filtered[(int)threadIdx.y-dy+HALO2][(int)threadIdx.x-dx+HALO2]);
        if (mag < total1 || mag < total2)
            mag = 0.0f;
    }
    if (yid < H && xid < W) {
        if (mask_in && !mask_in[(&img_in.at(bid, yid, xid, 0) - img_in.data)/3])
            mag = 0.0f;
        img_out.store1(bid, yid, xid, mag);
    }
}

/*[AutoHeaderGeneratorExport]*/
void canny_edge_filter_tensor(
    DeviceTensor3D<float3> img_in,
    bool* mask_in_ptr,
    DeviceTensor3D<float> img_out
) {
    int B  = img_in.template size<0>();
    int H  = img_in.template size<1>();
    int W  = img_in.template size<2>();

    canny_edge_filter_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        _dt3d_to_tv4<float>(img_in),
        mask_in_ptr,
        _dt3d_float_to_tv4_1ch(img_out)
    );
}


// ============================================================================
// Robust residual canny: RobustNeRF-style Tukey biweight on |render - GT| in
// BT.601 luma, followed by canny edge detection. Densification guidance that
// is simultaneously (a) zero where the render already matches GT, (b) tolerant
// to global luminance drift (DC residuals contribute no spatial gradient to
// canny), (c) suppressing distractor pixels (residuals above the per-image
// q-quantile get downweighted to zero).
//
// VRAM footprint per call: one [B,H,W] float scratch (pool-keyed
// "densify_robust_resid") + one [B] float Tukey threshold (pool-keyed
// "densify_tukey_c"). Both stable across calls/scales -- the pool's
// high-water-mark sizing handles the largest scale's allocation and reuses it
// for the smaller scales.
// ============================================================================

// Stage 1: per-pixel |render - GT| projected onto BT.601 luma. Luma matches
// canny's internal RGB->luma weighting, so the downstream edge detector sees
// the same scalar field a color-blind canny would have computed on the
// post-difference image.
__global__ void _robust_residual_luma_kernel(
    int B, int H, int W,
    const float3* __restrict__ render,
    const float3* __restrict__ ref,
    float* __restrict__ out  // [B, H, W]
) {
    int xid = blockIdx.x * blockDim.x + threadIdx.x;
    int yid = blockIdx.y * blockDim.y + threadIdx.y;
    int bid = blockIdx.z * blockDim.z + threadIdx.z;
    if (xid >= W || yid >= H || bid >= B) return;
    size_t idx = ((size_t)bid * H + yid) * W + xid;
    float3 d;
    d.x = render[idx].x - ref[idx].x;
    d.y = render[idx].y - ref[idx].y;
    d.z = render[idx].z - ref[idx].z;
    out[idx] = fabsf(0.299f * d.x + 0.587f * d.y + 0.114f * d.z);
}

// Stage 3: in-place Tukey biweight. Input `data[i]` holds |r_i| (luma resid);
// overwritten with |r_i| * w(r_i) where w(r) = (1 - (r/c)^2)^2 for |r| <= c
// else 0. The product r * w (the psi-style influence magnitude) is what the
// densification map should peak at: small for well-reconstructed regions
// (|r| ~ 0), large in the middle range, and back to zero past the cutoff
// (likely-distractor pixels). c is the per-image q-quantile of |r|.
__global__ void _robust_tukey_inplace_kernel(
    int B, int N,
    float* __restrict__ data,        // [B, N]
    const float* __restrict__ c_buf  // [B]
) {
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)B * (size_t)N;
    if (idx >= total) return;
    int b = (int)(idx / (size_t)N);
    float c = c_buf[b];
    float r = data[idx];
    if (!(c > 0.0f) || r >= c) {
        data[idx] = 0.0f;
        return;
    }
    float u = r / c;
    float w = 1.0f - u * u;
    data[idx] = r * w * w;
}

// Scalar-input canny: identical body to canny_edge_filter_kernel, only the
// input load differs (load1 vs RGB->luma dot product). Kept as a separate
// kernel rather than templatized because the body uses shared memory and
// the indirection through a __device__ functor would inhibit some loop /
// unroll heuristics.
__global__ void canny_edge_filter_kernel_scalar(
    TensorView<float, 4> img_in,    // [B, H, W, 1]
    const bool* __restrict__ mask_in,
    TensorView<float, 4> img_out    // [B, H, W, 1]
) {
    constexpr int BLOCK = 32;
    constexpr int HALO = 4;
    constexpr int HALO1 = 2;
    constexpr int HALO2 = 1;

    const int32_t xid = blockIdx.x * BLOCK + threadIdx.x;
    const int32_t yid = blockIdx.y * BLOCK + threadIdx.y;
    const int32_t bid = blockIdx.z;
    const uint32_t H = img_in.shape[1], W = img_in.shape[2];

    __shared__ float shared_pixels[BLOCK+2*HALO][BLOCK+2*HALO];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO)*(BLOCK+2*HALO); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO);
        int x = tid % (BLOCK+2*HALO);
        if (y < BLOCK+2*HALO) {
            int yi = min(max((int)(blockIdx.y * BLOCK) + y - HALO, 0), H-1);
            int xi = min(max((int)(blockIdx.x * BLOCK) + x - HALO, 0), W-1);
            shared_pixels[y][x] = img_in.load1(bid, yi, xi);
        }
    }
    __syncthreads();

    __shared__ float shared_blurred[BLOCK+2*HALO1][BLOCK+2*HALO1];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO1)*(BLOCK+2*HALO1); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO1);
        int x = tid % (BLOCK+2*HALO1);
        if (y >= BLOCK+2*HALO1)
            continue;
        float total = 0.0f;
        for (int cy = -2; cy <= 2; ++cy)
            for (int cx = -2; cx <= 2; ++cx) {
                float conv_weight = blur_filter_5x5[cy+2][cx+2];
                int yi = y - HALO1 + cy;
                int xi = x - HALO1 + cx;
                total += conv_weight * shared_pixels[yi+HALO][xi+HALO];
            }
        shared_blurred[y][x] = total;
    }
    __syncthreads();

    __shared__ float2 shared_filtered[BLOCK+2*HALO2][BLOCK+2*HALO2];
    #pragma unroll
    for (int batch = 0; batch < (BLOCK+2*HALO2)*(BLOCK+2*HALO2); batch += BLOCK*BLOCK) {
        int tid = batch + threadIdx.y * BLOCK + threadIdx.x;
        int y = tid / (BLOCK+2*HALO2);
        int x = tid % (BLOCK+2*HALO2);
        if (y >= BLOCK+2*HALO2)
            continue;
        float total1 = 0.0f, total2 = 0.0f;
        for (int cy = -1; cy <= 1; ++cy)
            for (int cx = -1; cx <= 1; ++cx) {
                float conv_weight_1 = canny_filter_3x3[cy+1][cx+1];
                float conv_weight_2 = canny_filter_3x3[cx+1][cy+1];
                int yi = y - HALO2 + cy;
                int xi = x - HALO2 + cx;
                float value = shared_blurred[yi+HALO1][xi+HALO1];
                total1 += conv_weight_1 * value;
                total2 += conv_weight_2 * value;
            }
        shared_filtered[y][x] = {total1, total2};
    }
    __syncthreads();

    float2 total = shared_filtered[threadIdx.y+HALO2][threadIdx.x+HALO2];
    float mag = length(total);
    if (mag > 0.0f) {
        int dx = min(max((int)roundf(total.x / mag), -HALO2), HALO2);
        int dy = min(max((int)roundf(total.y / mag), -HALO2), HALO2);
        float total1 = length(shared_filtered[(int)threadIdx.y+dy+HALO2][(int)threadIdx.x+dx+HALO2]);
        float total2 = length(shared_filtered[(int)threadIdx.y-dy+HALO2][(int)threadIdx.x-dx+HALO2]);
        if (mag < total1 || mag < total2)
            mag = 0.0f;
    }
    if (yid < H && xid < W) {
        if (mask_in && !mask_in[(&img_in.at(bid, yid, xid, 0) - img_in.data)])
            mag = 0.0f;
        img_out.store1(bid, yid, xid, mag);
    }
}

/*[AutoHeaderGeneratorExport]*/
void robust_canny_residual_tensor(
    DeviceTensor3D<float3> render,   // [B, H, W, 3]
    DeviceTensor3D<float3> ref,      // [B, H, W, 3]
    bool* mask_in_ptr,               // optional [B*H*W] mask; nullptr for none
    float quantile,                  // Tukey cutoff = per-image q-quantile of |r|
    DeviceTensor3D<float> img_out    // [B, H, W, 1] -- written (not added)
) {
    int B = render.template size<0>();
    int H = render.template size<1>();
    int W = render.template size<2>();
    if (B <= 0 || H <= 0 || W <= 0) return;
    size_t N = (size_t)H * (size_t)W;

    // Pool scratch buffers. Same keys across scales/calls -- pool's high-water
    // mark holds the largest, smaller scales reuse the allocation.
    float* resid = DevicePool::global().acquire<float>(
        "densify_robust_resid", (size_t)B * N);
    float* c_buf = DevicePool::global().acquire<float>(
        "densify_tukey_c", (size_t)B);

    // 1) per-pixel BT.601-luma |render - ref|.
    _robust_residual_luma_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        B, H, W,
        (const float3*)render.data_ptr(),
        (const float3*)ref.data_ptr(),
        resid
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // 2) per-image q-quantile of |r| -> Tukey cutoff c.
    quantile_of_abs_of_finite_elements_internal(
        resid, B, (int)N, quantile, /*return_reciprocal=*/false, c_buf);

    // 3) in-place Tukey biweight: resid <- |r| * (1 - (r/c)^2)^2, zero past c.
    _robust_tukey_inplace_kernel<<<_LAUNCH_ARGS_1D((size_t)B * N, 256)>>>(
        B, (int)N, resid, c_buf
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // 4) canny on the robustified scalar -> img_out (overwrites).
    TensorView<float, 4> resid_tv;
    resid_tv.data = resid;
    resid_tv.shape[0] = B;
    resid_tv.shape[1] = H;
    resid_tv.shape[2] = W;
    resid_tv.shape[3] = 1;
    resid_tv.strides[0] = (long)N;
    resid_tv.strides[1] = W;
    resid_tv.strides[2] = 1;
    resid_tv.strides[3] = 1;
    canny_edge_filter_kernel_scalar<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        resid_tv,
        mask_in_ptr,
        _dt3d_float_to_tv4_1ch(img_out)
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
