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

#include "common.cuh"

#include <cub/cub.cuh>


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

/*[AutoHeaderGeneratorExport]*/
at::Tensor quantile_of_abs_of_finite_elements_tensor(
    at::Tensor inputs,
    float q,
    bool return_reciprocal
) {
    DEVICE_GUARD(inputs);
    CHECK_INPUT(inputs);

    int B = inputs.ndimension() > 1 ? inputs.size(0) : 1;
    int N = inputs.numel() / B;
    at::Tensor outputs = at::empty({B}, inputs.options());
    if (B == 0)
        return outputs;
    at::Tensor temp = at::empty({1024*B}, inputs.options());

    (return_reciprocal ? batch_quantile_masked_radix_select<true> :
        batch_quantile_masked_radix_select<false>
    )(
        inputs.data_ptr<float>(),
        B, N, q,
        outputs.data_ptr<float>(),
        (uint32_t*)temp.data_ptr<float>(),
        at::cuda::getCurrentCUDAStream()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return outputs;
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
    at::Tensor data
) {
    at::Tensor inv_median = quantile_of_abs_of_finite_elements_tensor(data, 0.5, true);

    int B = data.ndimension() > 1 ? data.size(0) : 1;
    int N = data.numel() / B;

    multiply_by_inverse_median_kernel<<<_LAUNCH_ARGS_1D(B*N, 256)>>>(
        B, N,
        data.data_ptr<float>(),
        inv_median.data_ptr<float>()
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
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
) {
    DEVICE_GUARD(indices);
    CHECK_INPUT(indices);
    CHECK_INPUT(src);
    CHECK_INPUT(dst);

    if (indices.numel() == 0)
        return;
    index_kernel<<<_LAUNCH_ARGS_1D(dst.numel(), 256)>>>(
        dst.numel(),
        dst.numel() / indices.numel(),
        indices.data_ptr<int32_t>(),
        src.data_ptr<float>(),
        dst.data_ptr<float>()
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
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
) {
    DEVICE_GUARD(indices);
    CHECK_INPUT(indices);
    CHECK_INPUT(src);
    CHECK_INPUT(dst);

    if (indices.numel() == 0)
        return;
    scatter_add_kernel<<<_LAUNCH_ARGS_1D(src.numel(), 256)>>>(
        src.numel(),
        src.numel() / indices.numel(),
        indices.data_ptr<int32_t>(),
        src.data_ptr<float>(),
        dst.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void inplace_scatter_max_tensor(
    at::Tensor indices,
    at::Tensor src,
    at::Tensor dst
) {
    DEVICE_GUARD(indices);
    CHECK_INPUT(indices);
    CHECK_INPUT(src);
    CHECK_INPUT(dst);

    if (indices.numel() == 0)
        return;
    scatter_max_kernel<<<_LAUNCH_ARGS_1D(src.numel(), 256)>>>(
        src.numel(),
        src.numel() / indices.numel(),
        indices.data_ptr<int32_t>(),
        src.data_ptr<float>(),
        dst.data_ptr<float>()
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

/*[AutoHeaderGeneratorExport]*/
at::Tensor weighted_sample_without_replacement_tensor(
    int64_t numel,
    at::Tensor weights,
    std::optional<at::Tensor> masks,
    uint32_t num_sample,
    uint32_t seed
) {
    DEVICE_GUARD(weights);
    CHECK_INPUT(weights);
    if (masks.has_value())
        CHECK_INPUT(masks.value());

    if (numel == -1)
        numel = weights.ndimension() == 1 ?
            weights.numel() : weights.numel() / weights.size(-1);
    at::Tensor sorting_values = at::empty({numel}, weights.options());

    compute_efraimidis_spirakis_weight_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        numel,
        weights.numel() / numel,
        seed,
        weights.data_ptr<float>(),
        masks.has_value() ? masks.value().data_ptr<bool>() : nullptr,
        sorting_values.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    at::Tensor out_idx = at::empty({(long)num_sample}, weights.options().dtype(at::kInt));

    at::Tensor d_keys_in = sorting_values;  // reuse
    at::Tensor d_keys_out = at::empty_like(d_keys_in);

    at::Tensor d_indices_in  =
        at::empty({numel}, weights.options().dtype(at::kInt));
    at::Tensor d_indices_out =
        at::empty({numel}, weights.options().dtype(at::kInt));

    iota_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        numel,
        d_indices_in.data_ptr<int>()
    );

    void* d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;

    cub::DeviceRadixSort::SortPairs(
        d_temp_storage,
        temp_storage_bytes,
        d_keys_in.data_ptr<float>(),
        d_keys_out.data_ptr<float>(),
        d_indices_in.data_ptr<int>(),
        d_indices_out.data_ptr<int>(),
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    at::Tensor workspace =
        at::empty({(long)(temp_storage_bytes)},
                  weights.options().dtype(at::kByte));
    d_temp_storage = workspace.data_ptr<uint8_t>();

    cub::DeviceRadixSort::SortPairs(
        d_temp_storage,
        temp_storage_bytes,
        d_keys_in.data_ptr<float>(),
        d_keys_out.data_ptr<float>(),
        d_indices_in.data_ptr<int>(),
        d_indices_out.data_ptr<int>(),
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    cudaMemcpy(
        out_idx.data_ptr<int>(),
        d_indices_out.data_ptr<int>(),
        sizeof(int) * num_sample,
        cudaMemcpyDeviceToDevice
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_idx;
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
at::Tensor cov_scale_init_tensor(
    at::Tensor points,  // [N, 3]
    at::Tensor is_fisheye,  // [C], bool
    at::Tensor sizes,  // [C, 2], int32
    at::Tensor intrins,  // [C, 4]
    at::Tensor viewmats,  // [C, 4, 4]
    CameraDistortionCoeffsTensor dist_coeffs // [C]
) {
    DEVICE_GUARD(points);
    CHECK_INPUT(points);
    CHECK_INPUT(is_fisheye);
    CHECK_INPUT(sizes);
    CHECK_INPUT(intrins);
    CHECK_INPUT(viewmats);

    int64_t N = points.size(0);
    int64_t C = intrins.size(0);
    if (points.numel() != 3*N || points.size(-1) != 3)
        AT_ERROR("points shape must be (N, 3)");
    if (intrins.numel() != 4*C || intrins.size(-1) != 4)
        AT_ERROR("intrins shape must be (C, 4)");
    if (is_fisheye.numel() != C)
        AT_ERROR("is_fisheye shape must be (C,)");
    if (sizes.numel() != 2*C || sizes.size(-1) != 2)
        AT_ERROR("sizes shape must be (C, 2)");
    if (viewmats.numel() != C*4*4)
        AT_ERROR("viewmats shape must be (C, 4, 4)");

    at::Tensor log_scales = at::empty({N, 1}, points.options());

    cov_scale_init_kernel<<<_LAUNCH_ARGS_1D(N, 256)>>>(
        N, C,
        (float3*)points.data_ptr<float>(),
        is_fisheye.data_ptr<bool>(),
        (int2*)sizes.data_ptr<int32_t>(),
        (float4*)intrins.data_ptr<float>(),
        (float4*)viewmats.data_ptr<float>(),
        dist_coeffs,
        log_scales.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return log_scales;
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
    at::Tensor radii,
    at::Tensor log_scales,
    std::optional<at::Tensor> logit_opacs,
    float max_scale2d,
    float clip_hardness,
    float max_scale3d
) {
    DEVICE_GUARD(radii);
    CHECK_INPUT(radii);
    CHECK_INPUT(log_scales);
    if (logit_opacs.has_value())
        CHECK_INPUT(logit_opacs.value());

    densify_clip_scale_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, max_scale2d, clip_hardness, max_scale3d,
        radii.data_ptr<float>(),
        (float3*)log_scales.data_ptr<float>(),
        logit_opacs.has_value() ? logit_opacs.value().data_ptr<float>() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


inline __device__ float ellpise_surface_area(float3 scale) {
    return powf(
        powf(scale.x*scale.y, 1.6f) +
        powf(scale.x*scale.z, 1.6f) +
        powf(scale.y*scale.z, 1.6f), 1.0f / 1.6f);
}


__global__ void densify_update_weight_kernel(
    long num_splats,
    bool is_max_mode,
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
    if (true) // TODO
        // weight *= rsqrtf(radii[idx]);
        weight /= fmaxf(radii[idx], 1e-10f);
    if (scales) {
        float3 scale = scales[idx];
        scale = {__expf(scale.x), __expf(scale.y), __expf(scale.z)};
        float sqrt_visible_area_est =
            radii[idx] / fmaxf(fmaxf(scale.x, scale.y), scale.z) * sqrtf(ellpise_surface_area(scale));
        weight /= fmaxf(sqrt_visible_area_est, 1e-10f);
    }
    if (accum_weight_scalar != nullptr)
        weight *= accum_weight_scalar[0];
    float2 accum = accum_buffer[idx];
    accum.x *= accum.y;
    if (is_max_mode) {
        accum.x = fmaxf(accum.x, weight);
        accum.y = fmaxf(accum.y, 1.0f);
    } else {
        accum.x += weight;
        accum.y += 1.0f;
    }
    if (accum.y > 0.0f)
        accum.x /= accum.y;
    accum_buffer[idx] = accum;
}


/*[AutoHeaderGeneratorExport]*/
void densify_update_weight_tensor(
    int64_t num_splats,
    at::Tensor radii,
    std::optional<at::Tensor> scales,
    std::optional<at::Tensor> opacs,
    at::Tensor accum_weight,
    at::Tensor accum_buffer,
    bool is_max_mode
) {
    DEVICE_GUARD(radii);
    CHECK_INPUT(radii);
    CHECK_INPUT(accum_weight);
    CHECK_INPUT(accum_buffer);
    if (scales.has_value())
        CHECK_INPUT(scales.value());
    if (opacs.has_value())
        CHECK_INPUT(opacs.value());

    // at::Tensor inv_median = quantile_of_abs_of_finite_elements_tensor(
    //     accum_weight, 0.5f, true
    // );

    densify_update_weight_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        num_splats, is_max_mode,
        radii.data_ptr<float>(),
        scales.has_value() ? (float3*)scales.value().data_ptr<float>() : nullptr,
        opacs.has_value() ? opacs.value().data_ptr<float>() : nullptr,
        // inv_median.data_ptr<float>(),
        nullptr,
        accum_weight.data_ptr<float>(),
        (float2*)accum_buffer.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Relocation
// ================

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
    int num_sh, float3*__restrict__ features_sh, float3*__restrict__ g1_features_sh, float3*__restrict__ g2_features_sh,
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
    for (int i = 0; i < num_sh; ++i)  // TODO: slow; more cache friendly way to do so?
        features_sh[num_sh*idx_dst+i] = features_sh[num_sh*idx_src+i];

    // optimizer state - zero
#if 1
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
    for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
        g1_features_sh[num_sh*idx_dst+i] = make_float3(0.0f);
        g2_features_sh[num_sh*idx_dst+i] = make_float3(0.0f);
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
    for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
        g1_features_sh[num_sh*idx_src+i] = make_float3(0.0f);
        g2_features_sh[num_sh*idx_src+i] = make_float3(0.0f);
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
    at::Tensor means, at::Tensor quats, at::Tensor scales, at::Tensor opacs, at::Tensor features_dc, at::Tensor features_sh,
    at::Tensor g1_means, at::Tensor g1_quats, at::Tensor g1_scales, at::Tensor g1_opacs, at::Tensor g1_features_dc, at::Tensor g1_features_sh,
    at::Tensor g2_means, at::Tensor g2_quats, at::Tensor g2_scales, at::Tensor g2_opacs, at::Tensor g2_features_dc, at::Tensor g2_features_sh,
    at::Tensor densify_accum_buffer,
    std::optional<at::Tensor> bias_correction_steps,
    uint32_t seed
) {
    int num_sh = features_sh.size(-2);

    at::Tensor mask = at::empty({cur_num_splats}, densify_accum_buffer.options().dtype(at::kBool));
    at::Tensor num_relocate_tensor = at::empty({1}, densify_accum_buffer.options().dtype(at::kInt));
    set_zero<int32_t>(num_relocate_tensor);

    at::Tensor dst_indices = at::empty({cur_num_splats}, densify_accum_buffer.options().dtype(at::kInt));

    compute_relocation_mask_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        cur_num_splats,
        min_opacity,
        (float3*)means.data_ptr<float>(),
        (float4*)quats.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(),
        (float*)opacs.data_ptr<float>(),
        (float3*)features_dc.data_ptr<float>(),
        mask.data_ptr<bool>(),
        num_relocate_tensor.data_ptr<int32_t>(),
        dst_indices.data_ptr<int32_t>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    int64_t num_relocate = (int64_t)num_relocate_tensor.item<int32_t>();
    if (num_relocate == 0)
        return;

    at::Tensor src_indices = weighted_sample_without_replacement_tensor(
        cur_num_splats, densify_accum_buffer, mask, num_relocate, seed);

    relocate_with_long_axis_split_kernel<<<_LAUNCH_ARGS_1D(num_relocate, 256)>>>(
        cur_num_splats,
        num_relocate,
        src_indices.data_ptr<int32_t>(),
        dst_indices.data_ptr<int32_t>(),
        (float3*)means.data_ptr<float>(), (float3*)g1_means.data_ptr<float>(), (float3*)g2_means.data_ptr<float>(),
        (float4*)quats.data_ptr<float>(), (float4*)g1_quats.data_ptr<float>(), (float4*)g2_quats.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(), (float3*)g1_scales.data_ptr<float>(), (float3*)g2_scales.data_ptr<float>(),
        (float*)opacs.data_ptr<float>(), (float*)g1_opacs.data_ptr<float>(), (float*)g2_opacs.data_ptr<float>(),
        (float3*)features_dc.data_ptr<float>(), (float3*)g1_features_dc.data_ptr<float>(), (float3*)g2_features_dc.data_ptr<float>(),
        num_sh, (float3*)features_sh.data_ptr<float>(), (float3*)g1_features_sh.data_ptr<float>(), (float3*)g2_features_sh.data_ptr<float>(),
        (float2*)densify_accum_buffer.data_ptr<float>(),
        bias_correction_steps.has_value() ? bias_correction_steps.value().data_ptr<int32_t>() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void add_splats_with_long_axis_split_tensor(
    int64_t cur_num_splats,
    int64_t num_new_splats,
    at::Tensor means, at::Tensor quats, at::Tensor scales, at::Tensor opacs, at::Tensor features_dc, at::Tensor features_sh,
    at::Tensor g1_means, at::Tensor g1_quats, at::Tensor g1_scales, at::Tensor g1_opacs, at::Tensor g1_features_dc, at::Tensor g1_features_sh,
    at::Tensor g2_means, at::Tensor g2_quats, at::Tensor g2_scales, at::Tensor g2_opacs, at::Tensor g2_features_dc, at::Tensor g2_features_sh,
    at::Tensor densify_accum_buffer,
    std::optional<at::Tensor> bias_correction_steps,
    uint32_t seed
) {
    if (num_new_splats == 0)
        return;

    int num_sh = features_sh.size(-2);

    at::Tensor split_indices = weighted_sample_without_replacement_tensor(
        cur_num_splats, densify_accum_buffer, std::nullopt, num_new_splats, seed);

    relocate_with_long_axis_split_kernel<<<_LAUNCH_ARGS_1D(num_new_splats, 256)>>>(
        cur_num_splats,
        num_new_splats,
        split_indices.data_ptr<int32_t>(),
        nullptr,
        (float3*)means.data_ptr<float>(), (float3*)g1_means.data_ptr<float>(), (float3*)g2_means.data_ptr<float>(),
        (float4*)quats.data_ptr<float>(), (float4*)g1_quats.data_ptr<float>(), (float4*)g2_quats.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(), (float3*)g1_scales.data_ptr<float>(), (float3*)g2_scales.data_ptr<float>(),
        (float*)opacs.data_ptr<float>(), (float*)g1_opacs.data_ptr<float>(), (float*)g2_opacs.data_ptr<float>(),
        (float3*)features_dc.data_ptr<float>(), (float3*)g1_features_dc.data_ptr<float>(), (float3*)g2_features_dc.data_ptr<float>(),
        num_sh, (float3*)features_sh.data_ptr<float>(), (float3*)g1_features_sh.data_ptr<float>(), (float3*)g2_features_sh.data_ptr<float>(),
        (float2*)densify_accum_buffer.data_ptr<float>(),
        bias_correction_steps.has_value() ? bias_correction_steps.value().data_ptr<int32_t>() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
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


__global__ void mcmc_compute_relocation_kernel(
    uint32_t num_splats,
    float min_opacity,
    int32_t* __restrict__ n_idx_buffer,  // [N]
    float3*__restrict__ means, float3*__restrict__ g1_means, float3*__restrict__ g2_means,
    float4*__restrict__ quats, float4*__restrict__ g1_quats, float4*__restrict__ g2_quats,
    float3*__restrict__ scales, float3*__restrict__ g1_scales, float3*__restrict__ g2_scales,
    float*__restrict__ opacs, float*__restrict__ g1_opacs, float*__restrict__ g2_opacs,
    float3*__restrict__ features_dc, float3*__restrict__ g1_features_dc, float3*__restrict__ g2_features_dc,
    int num_sh, float3*__restrict__ features_sh, float3*__restrict__ g1_features_sh, float3*__restrict__ g2_features_sh,
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

    // set grad to zero
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
    for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
        g1_features_sh[num_sh*cur_idx+i] = make_float3(0.0f);
        g2_features_sh[num_sh*cur_idx+i] = make_float3(0.0f);
    }
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
    for (int i = 0; i < num_sh; ++i)  // TODO: slow; more cache friendly way to do so?
        features_sh[num_sh*id_dst+i] = features_sh[num_sh*id_src+i];
}


/*[AutoHeaderGeneratorExport]*/
void relocate_splats_mcmc_tensor(
    int64_t cur_num_splats,
    float min_opacity,
    at::Tensor means, at::Tensor quats, at::Tensor scales, at::Tensor opacs, at::Tensor features_dc, at::Tensor features_sh,
    at::Tensor g1_means, at::Tensor g1_quats, at::Tensor g1_scales, at::Tensor g1_opacs, at::Tensor g1_features_dc, at::Tensor g1_features_sh,
    at::Tensor g2_means, at::Tensor g2_quats, at::Tensor g2_scales, at::Tensor g2_opacs, at::Tensor g2_features_dc, at::Tensor g2_features_sh,
    std::optional<at::Tensor> bias_correction_steps,
    uint32_t seed
) {
    int num_sh = features_sh.size(-2);

    at::Tensor sample_probs = at::empty_like(opacs);
    mcmc_compute_relocation_probabilities_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        cur_num_splats,
        min_opacity,
        opacs.data_ptr<float>(),
        sample_probs.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    at::Tensor sample_probs_cumsum = at::cumsum(sample_probs, 0);

    at::Tensor index_map = at::empty({cur_num_splats}, opacs.options().dtype(at::kInt));
    at::Tensor n_idx_buffer = at::empty_like(index_map);
    set_zero<int32_t>(n_idx_buffer);

    mcmc_compute_relocation_index_map_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        sample_probs.data_ptr<float>(),
        sample_probs_cumsum.data_ptr<float>(),
        index_map.data_ptr<int32_t>(),
        n_idx_buffer.data_ptr<int32_t>(),
        cur_num_splats,
        seed
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // printf("%d relocate\n", n_idx_buffer.sum().item<int32_t>());

    mcmc_compute_relocation_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 64)>>>(
        cur_num_splats,
        min_opacity,
        n_idx_buffer.data_ptr<int32_t>(),
        (float3*)means.data_ptr<float>(), (float3*)g1_means.data_ptr<float>(), (float3*)g2_means.data_ptr<float>(),
        (float4*)quats.data_ptr<float>(), (float4*)g1_quats.data_ptr<float>(), (float4*)g2_quats.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(), (float3*)g1_scales.data_ptr<float>(), (float3*)g2_scales.data_ptr<float>(),
        (float*)opacs.data_ptr<float>(), (float*)g1_opacs.data_ptr<float>(), (float*)g2_opacs.data_ptr<float>(),
        (float3*)features_dc.data_ptr<float>(), (float3*)g1_features_dc.data_ptr<float>(), (float3*)g2_features_dc.data_ptr<float>(),
        num_sh, (float3*)features_sh.data_ptr<float>(), (float3*)g1_features_sh.data_ptr<float>(), (float3*)g2_features_sh.data_ptr<float>(),
        bias_correction_steps.has_value() ? bias_correction_steps.value().data_ptr<int32_t>() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    mcmc_update_relocation_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 64)>>>(
        cur_num_splats,
        index_map.data_ptr<int32_t>(),
        (float3*)means.data_ptr<float>(),
        (float4*)quats.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(),
        (float*)opacs.data_ptr<float>(),
        (float3*)features_dc.data_ptr<float>(),
        num_sh, (float3*)features_sh.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
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


__global__ void mcmc_update_add_kernel(
    uint32_t num_splats,
    uint32_t num_add,
    int32_t* __restrict__ index_map,  // [dN]
    float3*__restrict__ means, float3*__restrict__ g1_means, float3*__restrict__ g2_means,
    float4*__restrict__ quats, float4*__restrict__ g1_quats, float4*__restrict__ g2_quats,
    float3*__restrict__ scales, float3*__restrict__ g1_scales, float3*__restrict__ g2_scales,
    float*__restrict__ opacs, float*__restrict__ g1_opacs, float*__restrict__ g2_opacs,
    float3*__restrict__ features_dc, float3*__restrict__ g1_features_dc, float3*__restrict__ g2_features_dc,
    int num_sh, float3*__restrict__ features_sh, float3*__restrict__ g1_features_sh, float3*__restrict__ g2_features_sh,
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
    for (int i = 0; i < num_sh; ++i)  // TODO: slow; more cache friendly way to do so?
        features_sh[num_sh*id_dst+i] = features_sh[num_sh*id_src+i];

    // set grad to zero
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
    for (int i = 0; i < num_sh; ++i) {  // TODO: slow; more cache friendly way to do so?
        g1_features_sh[num_sh*id_dst+i] = make_float3(0.0f);
        g2_features_sh[num_sh*id_dst+i] = make_float3(0.0f);
    }
    if (bias_correction_steps)
        bias_correction_steps[id_dst] = 0;
}


/*[AutoHeaderGeneratorExport]*/
void add_splats_mcmc_tensor(
    int64_t cur_num_splats,
    int64_t num_add,
    float min_opacity,
    at::Tensor means, at::Tensor quats, at::Tensor scales, at::Tensor opacs, at::Tensor features_dc, at::Tensor features_sh,
    at::Tensor g1_means, at::Tensor g1_quats, at::Tensor g1_scales, at::Tensor g1_opacs, at::Tensor g1_features_dc, at::Tensor g1_features_sh,
    at::Tensor g2_means, at::Tensor g2_quats, at::Tensor g2_scales, at::Tensor g2_opacs, at::Tensor g2_features_dc, at::Tensor g2_features_sh,
    std::optional<at::Tensor> bias_correction_steps,
    uint32_t seed
) {
    if (num_add == 0)
        return;

    int num_sh = features_sh.size(-2);

    at::Tensor sample_probs = at::empty_like(opacs);
    mcmc_compute_relocation_probabilities_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 256)>>>(
        cur_num_splats,
        min_opacity,
        opacs.data_ptr<float>(),
        sample_probs.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    at::Tensor sample_probs_cumsum = at::cumsum(sample_probs, 0);

    at::Tensor index_map = at::empty({num_add}, opacs.options().dtype(at::kInt));
    at::Tensor n_idx_buffer = at::empty({cur_num_splats}, opacs.options().dtype(at::kInt));
    set_zero<int32_t>(n_idx_buffer);

    mcmc_compute_add_index_map_kernel<<<_LAUNCH_ARGS_1D(num_add, 256)>>>(
        sample_probs_cumsum.data_ptr<float>(),
        index_map.data_ptr<int32_t>(),
        n_idx_buffer.data_ptr<int32_t>(),
        cur_num_splats,
        num_add,
        seed
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    // printf("%d/%d add\n", n_idx_buffer.sum().item<int32_t>(), (int)num_add);

    mcmc_compute_add_kernel<<<_LAUNCH_ARGS_1D(cur_num_splats, 64)>>>(
        cur_num_splats,
        min_opacity,
        n_idx_buffer.data_ptr<int32_t>(),
        (float3*)scales.data_ptr<float>(),
        (float*)opacs.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    mcmc_update_add_kernel<<<_LAUNCH_ARGS_1D(num_add, 64)>>>(
        cur_num_splats,
        num_add,
        index_map.data_ptr<int32_t>(),
        (float3*)means.data_ptr<float>(), (float3*)g1_means.data_ptr<float>(), (float3*)g2_means.data_ptr<float>(),
        (float4*)quats.data_ptr<float>(), (float4*)g1_quats.data_ptr<float>(), (float4*)g2_quats.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(), (float3*)g1_scales.data_ptr<float>(), (float3*)g2_scales.data_ptr<float>(),
        (float*)opacs.data_ptr<float>(), (float*)g1_opacs.data_ptr<float>(), (float*)g2_opacs.data_ptr<float>(),
        (float3*)features_dc.data_ptr<float>(), (float3*)g1_features_dc.data_ptr<float>(), (float3*)g2_features_dc.data_ptr<float>(),
        num_sh, (float3*)features_sh.data_ptr<float>(), (float3*)g1_features_sh.data_ptr<float>(), (float3*)g2_features_sh.data_ptr<float>(),
        bias_correction_steps.has_value() ? bias_correction_steps.value().data_ptr<int32_t>() : nullptr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// MCMC Add Noise
// ================


__global__ void mcmc_add_noise_3dgs_kernel(
    long num_splats,
    float scaler, float min_opacity,
    float3* __restrict__ means,
    const float3* __restrict__ log_scales,
    const float4* __restrict__ quats,
    const float* __restrict__ opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::mcmc_add_noise_3dgs(
        scaler, min_opacity,
        &means[idx], log_scales[idx], quats[idx], opacs[idx]
    );
}

__global__ void mcmc_add_noise_triangle_kernel(
    long num_splats,
    float scaler, float min_opacity,
    float3* __restrict__ means,
    const float3* __restrict__ log_scales,
    const float4* __restrict__ quats,
    const float* __restrict__ opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::mcmc_add_noise_triangle(
        scaler, min_opacity,
        &means[idx], log_scales[idx], quats[idx], opacs[idx]
    );
}

/*[AutoHeaderGeneratorExport]*/
void mcmc_add_noise_tensor(
    std::string primitive,
    float scaler, float min_opacity,
    at::Tensor &means,
    at::Tensor &log_scales,
    at::Tensor &quats,
    at::Tensor &opacs
) {
    DEVICE_GUARD(means);
    CHECK_INPUT(means);
    CHECK_INPUT(log_scales);
    CHECK_INPUT(quats);
    CHECK_INPUT(opacs);

    const size_t num_splats = opacs.numel();

    if (primitive == "3dgs" || primitive == "mip" || primitive == "3dgut")
        mcmc_add_noise_3dgs_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
            num_splats, scaler, min_opacity,
            (float3*)means.data_ptr<float>(),
            (float3*)log_scales.data_ptr<float>(),
            (float4*)quats.data_ptr<float>(),
            opacs.data_ptr<float>()
        );
    else if (primitive == "opaque_triangle")
        mcmc_add_noise_triangle_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
            num_splats, scaler, min_opacity,
            (float3*)means.data_ptr<float>(),
            (float3*)log_scales.data_ptr<float>(),
            (float4*)quats.data_ptr<float>(),
            opacs.data_ptr<float>()
        );
    else throw std::runtime_error("Unsupported primitive: " + primitive);
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
std::tuple<at::Tensor, at::Tensor, at::Tensor>
long_axis_split_tensor(
    std::string primitive,
    at::Tensor &log_scales,
    at::Tensor &logit_opacities,
    at::Tensor &quats
) {
    DEVICE_GUARD(log_scales);
    CHECK_INPUT(log_scales);
    CHECK_INPUT(logit_opacities);
    CHECK_INPUT(quats);

    const size_t num_splats = quats.numel() / 4;

    at::Tensor new_log_scales = at::empty_like(log_scales);
    at::Tensor new_logit_opacities= at::empty_like(logit_opacities);
    at::Tensor mean_deltas = at::empty_like(log_scales);

    if (num_splats == 0)
        return std::make_tuple(new_log_scales, new_logit_opacities, mean_deltas);

    if (primitive == "3dgs" || primitive == "mip" || primitive == "3dgut")
        long_axis_split_3dgs_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
            num_splats,
            (float3*)log_scales.data_ptr<float>(),
            logit_opacities.data_ptr<float>(),
            (float4*)quats.data_ptr<float>(),
            (float3*)new_log_scales.data_ptr<float>(),
            new_logit_opacities.data_ptr<float>(),
            (float3*)mean_deltas.data_ptr<float>()
        );
    else throw std::runtime_error("Unsupported primitive: " + primitive);
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(new_log_scales, new_logit_opacities, mean_deltas);
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
at::Tensor laplacian_edge_filter_tensor(
    at::Tensor &img_in
) {
    DEVICE_GUARD(img_in);
    CHECK_INPUT(img_in);

    if (img_in.ndimension() != 4 || img_in.size(-1) != 3)
        AT_ERROR("img must be [B, H, W, 3]");
    int B  = img_in.size(0);
    int H  = img_in.size(1);
    int W  = img_in.size(2);

    auto img_out = at::empty({B, H, W, 1}, img_in.options());

    laplacian_edge_filter_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        tensor2view<float, 4>(img_in),
        tensor2view<float, 4>(img_out)
    );

    return img_out;
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
at::Tensor smoothed_laplacian_edge_filter_tensor(
    at::Tensor &img_in
) {
    DEVICE_GUARD(img_in);
    CHECK_INPUT(img_in);

    if (img_in.ndimension() != 4 || img_in.size(-1) != 3)
        AT_ERROR("img must be [B, H, W, 3]");
    int B  = img_in.size(0);
    int H  = img_in.size(1);
    int W  = img_in.size(2);

    auto img_out = at::empty({B, H, W, 1}, img_in.options());

    smoothed_laplacian_edge_filter_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        tensor2view<float, 4>(img_in),
        tensor2view<float, 4>(img_out)
    );

    return img_out;
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
    if (mask_in && !mask_in[(&img_in.at(bid, yid, xid, 0) - img_in.data)/3])
        mag = 0.0f;
    if (yid < H && xid < W)
        img_out.store1(bid, yid, xid, mag);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor canny_edge_filter_tensor(
    at::Tensor &img_in,
    at::optional<at::Tensor> &mask_in
) {
    DEVICE_GUARD(img_in);
    CHECK_INPUT(img_in);
    if (mask_in.has_value())
        CHECK_INPUT(mask_in.value());

    if (img_in.ndimension() != 4 || img_in.size(-1) != 3)
        AT_ERROR("img must be [B, H, W, 3]");
    int B  = img_in.size(0);
    int H  = img_in.size(1);
    int W  = img_in.size(2);

    auto img_out = at::empty({B, H, W, 1}, img_in.options());

    canny_edge_filter_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        tensor2view<float, 4>(img_in),
        mask_in.has_value() ? mask_in.value().data_ptr<bool>() : nullptr,
        tensor2view<float, 4>(img_out)
    );

    return img_out;
}
