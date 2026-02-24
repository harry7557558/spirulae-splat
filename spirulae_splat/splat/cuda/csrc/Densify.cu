#include "Densify.cuh"

#ifdef __CUDACC__
#include "generated/slang.cuh"
namespace SlangDensify {
#include "generated/set_namespace.cuh"
#include "generated/densify.cuh"
}
#endif

#include "common.cuh"


// ================
// MCMC Add Noise
// ================


__global__ void mcmc_add_noise_3dgs_kernel(
    long num_splats,
    float scaler, float min_opacity,
    float3* __restrict__ means,
    const float3* __restrict__ scales,
    const float4* __restrict__ quats,
    const float* __restrict__ opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::mcmc_add_noise_3dgs(
        scaler, min_opacity,
        &means[idx], scales[idx], quats[idx], opacs[idx]
    );
}

__global__ void mcmc_add_noise_triangle_kernel(
    long num_splats,
    float scaler, float min_opacity,
    float3* __restrict__ means,
    const float3* __restrict__ scales,
    const float4* __restrict__ quats,
    const float* __restrict__ opacs
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::mcmc_add_noise_triangle(
        scaler, min_opacity,
        &means[idx], scales[idx], quats[idx], opacs[idx]
    );
}

/*[AutoHeaderGeneratorExport]*/
void mcmc_add_noise_tensor(
    std::string primitive,
    float scaler, float min_opacity,
    at::Tensor &means,
    at::Tensor &scales,
    at::Tensor &quats,
    at::Tensor &opacs
) {
    DEVICE_GUARD(means);
    CHECK_INPUT(means);
    CHECK_INPUT(scales);
    CHECK_INPUT(quats);
    CHECK_INPUT(opacs);

    const size_t num_splats = opacs.numel();

    if (primitive == "3dgs" || primitive == "mip" || primitive == "3dgut")
        mcmc_add_noise_3dgs_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
            num_splats, scaler, min_opacity,
            (float3*)means.data_ptr<float>(),
            (float3*)scales.data_ptr<float>(),
            (float4*)quats.data_ptr<float>(),
            opacs.data_ptr<float>()
        );
    else if (primitive == "opaque_triangle")
        mcmc_add_noise_triangle_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
            num_splats, scaler, min_opacity,
            (float3*)means.data_ptr<float>(),
            (float3*)scales.data_ptr<float>(),
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
    new_opacities[idx] = 1.0f - powf(1.0f - opacities[idx], 1.0f / n_idx);

    // compute new scale
    for (int i = 1; i <= n_idx; ++i) {
        for (int k = 0; k <= (i - 1); ++k) {
            float bin_coeff = binoms[(i - 1) * n_max + k];
            float term = (pow(-1.0f, k) / sqrt(static_cast<float>(k + 1))) *
                         pow(new_opacities[idx], k + 1);
            denom_sum += (bin_coeff * term);
        }
    }
    float coeff = (opacities[idx] / denom_sum);
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
    const float3* __restrict__ scales,
    const float4* __restrict__ quats,
    float3* __restrict__ new_scales,
    float3* __restrict__ mean_deltas
) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_splats)
        return;

    SlangDensify::long_axis_split_3dgs(
        scales[idx], quats[idx], &new_scales[idx], &mean_deltas[idx]
    );
}

/*[AutoHeaderGeneratorExport]*/
std::tuple<at::Tensor, at::Tensor>
long_axis_split_tensor(
    std::string primitive,
    at::Tensor &scales,
    at::Tensor &quats
) {
    DEVICE_GUARD(scales);
    CHECK_INPUT(scales);
    CHECK_INPUT(quats);

    const size_t num_splats = quats.numel() / 4;

    at::Tensor new_scales = at::empty_like(scales);
    at::Tensor mean_deltas = at::empty_like(scales);

    if (num_splats == 0)
        return std::make_tuple(new_scales, mean_deltas);

    if (primitive == "3dgs" || primitive == "mip" || primitive == "3dgut")
        long_axis_split_3dgs_kernel<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
            num_splats,
            (float3*)scales.data_ptr<float>(),
            (float4*)quats.data_ptr<float>(),
            (float3*)new_scales.data_ptr<float>(),
            (float3*)mean_deltas.data_ptr<float>()
        );
    else throw std::runtime_error("Unsupported primitive: " + primitive);
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return std::make_tuple(new_scales, mean_deltas);
}

