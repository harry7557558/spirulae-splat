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
// Covariance-Based Scale Initialization
// ================

__global__ void cov_scale_init_kernel(
    int64_t num_points,
    int32_t num_cameras,
    const float3* __restrict__ points,  // [N, 3]
    const bool* __restrict__ is_fisheye,  // [C]
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
// MCMC Relocation (from GSplat)
// ================

__global__ void compute_relocation_kernel(
    int N,
    float* __restrict__ opacities,
    float3* __restrict__ scales,
    int *ratios,
    float *binoms,
    int n_max,
    float *new_opacities,
    float3 *new_scales
) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx >= N)
        return;

    int n_idx = ratios[idx];
    float denom_sum = 0.0f;

    // compute new opacity
    float old_opacity = opacities[idx];
    float new_opacity = 1.0f - powf(1.0f - old_opacity, 1.0f / n_idx);
    new_opacity = fmaxf(new_opacity, 0.0f);
    new_opacities[idx] = new_opacity;

    // compute new scale
    for (int i = 1; i <= n_idx; ++i) {
        for (int k = 0; k <= (i - 1); ++k) {
            float bin_coeff = binoms[(i - 1) * n_max + k];
            float term = (((k & 1) ? -1.0f : 1.0f) / sqrtf((float)(k + 1))) *
                         powf(new_opacity, (float)(k + 1));
            denom_sum += (bin_coeff * term);
        }
    }
    float coeff = (old_opacity / denom_sum);
    new_scales[idx] = coeff * scales[idx];
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor>
compute_relocation_tensor(
    // inputs
    at::Tensor opacities, // [N]
    at::Tensor scales,    // [N, 3]
    at::Tensor ratios,    // [N]
    at::Tensor binoms,    // [n_max, n_max]
    const int n_max
) {
    DEVICE_GUARD(opacities);
    CHECK_INPUT(opacities);
    CHECK_INPUT(scales);
    CHECK_INPUT(ratios);
    CHECK_INPUT(binoms);

    uint32_t N = opacities.size(0);

    int64_t n_elements = N;

    at::Tensor new_opacities = at::empty_like(opacities);
    at::Tensor new_scales = at::empty_like(scales);

    if (n_elements == 0)
        return std::make_tuple(new_opacities, new_scales);

    compute_relocation_kernel<<<_LAUNCH_ARGS_1D(N, 256)>>>(
        N,
        opacities.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(),
        ratios.data_ptr<int>(),
        binoms.data_ptr<float>(),
        n_max,
        new_opacities.data_ptr<float>(),
        (float3*)new_scales.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(new_opacities, new_scales);
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
    const TensorView<float, 4> img_in,
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
    if (yid < H && xid < W)
        img_out.store1(bid, yid, xid, mag);
}

/*[AutoHeaderGeneratorExport]*/
at::Tensor canny_edge_filter_tensor(
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

    canny_edge_filter_kernel<<<_LAUNCH_ARGS_3D(W, H, B, 32, 32, 1)>>>(
        tensor2view<float, 4>(img_in),
        tensor2view<float, 4>(img_out)
    );

    return img_out;
}
