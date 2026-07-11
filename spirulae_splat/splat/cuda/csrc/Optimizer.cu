#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}
namespace SlangPixelWise {
#include "generated/set_namespace.cuh"
#include "generated/pixel_wise.cuh"
}
namespace SlangPerSplatLosses {
#include "generated/set_namespace.cuh"
#include "generated/per_splat_losses.cuh"
}

#include <Common.cuh>
#include <Tensor.h>
#include "NonShQuantState.h"
#include <stdexcept>
#include <string>

#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#if defined(__INTEL_COMPILER) || defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER)
#if defined(__AVX2__) || defined(__SSE2__)
#include <immintrin.h>
#endif
#endif




inline constexpr float kSh0 = 0.28209479177387814f;

// Helper: compute total number of elements from a TorchTensorView
static inline int64_t _tv_numel(const TorchTensorView& tv) {
    int64_t n = 1;
    for (auto d : std::get<2>(tv)) n *= d;
    return n;
}
// Helper: get float* from TorchTensorView
static inline float* _tv_f(const TorchTensorView& tv) { return (float*)std::get<0>(tv); }

__forceinline__ __device__ float3 fmaxf(float3 v, float k) {
    return {
        fmaxf(v.x, k),
        fmaxf(v.y, k),
        fmaxf(v.z, k)
    };
}
__forceinline__ __device__ float3 sqrtf(float3 v) {
    return {
        sqrtf(fmaxf(v.x, 0.0f)),
        sqrtf(fmaxf(v.y, 0.0f)),
        sqrtf(fmaxf(v.z, 0.0f))
    };
}

__forceinline__ __device__ float4 sqrtf(float4 v) {
    return {
        sqrtf(fmaxf(v.x, 0.0f)),
        sqrtf(fmaxf(v.y, 0.0f)),
        sqrtf(fmaxf(v.z, 0.0f)),
        sqrtf(fmaxf(v.w, 0.0f))
    };
}

__global__ void set_zero_kernel(size_t numel, float* data) {
    const size_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    if (idx < numel)
        data[idx] = 0.0f;
}

__global__ void set_zero_kernel(size_t numel, float4* data) {
    const size_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    if (idx < numel)
        data[idx] = float4{0.0f, 0.0f, 0.0f, 0.0f};
}

/*[AutoHeaderGeneratorExport]*/
void set_zero_tensor(TorchTensorView x) {
    int64_t numel = _tv_numel(x);
    float* ptr = _tv_f(x);
    if (numel % 4 == 0)
        set_zero_kernel<<<_LAUNCH_ARGS_1D(numel/4, 512)>>>
            (numel/4, (float4*)ptr);
    else
        set_zero_kernel<<<_LAUNCH_ARGS_1D(numel, 512)>>>
            (numel, ptr);
}


// ================
// "Standard" Adam
// ================

template <typename scalar_t>
__global__ void fused_adam_kernel(
    scalar_t* __restrict__ param,
    const scalar_t* __restrict__ grad,
    scalar_t* __restrict__ exp_avg,
    scalar_t* __restrict__ exp_avg_sq,
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const int64_t numel
) {
    const int64_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    
    if (idx < numel) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;
        
        // Load values
        float grad_val = static_cast<float>(grad[idx]);
        if (!isfinite(grad_val))  // TODO: gradient clipping
            grad_val = 0.0f;
        float m_val = static_cast<float>(exp_avg[idx]);
        float v_val = static_cast<float>(exp_avg_sq[idx]);
        
        // Update biased first moment estimate: m_t = beta1 * m_{t-1} + (1 - beta1) * g_t
        m_val = beta1 * m_val + (1.0f - beta1) * grad_val;
        
        // Update biased second raw moment estimate: v_t = beta2 * v_{t-1} + (1 - beta2) * g_t^2
        v_val = beta2 * v_val + (1.0f - beta2) * grad_val * grad_val;
        
        // Compute update: theta_t = theta_{t-1} - step_size * m_t / (sqrt(v_t / bias_correction2) + eps)
        const float denom = sqrtf(v_val / bias_correction2) + eps;
        const float param_val = static_cast<float>(param[idx]) - step_size * (m_val / denom);
        
        // Write back
        param[idx] = static_cast<scalar_t>(param_val);
        exp_avg[idx] = static_cast<scalar_t>(m_val);
        exp_avg_sq[idx] = static_cast<scalar_t>(v_val);
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
) {
    const int64_t numel = _tv_numel(param);
    if (numel == 0)
        return;

    fused_adam_kernel<float><<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        _tv_f(param),
        _tv_f(grad),
        _tv_f(exp_avg),
        _tv_f(exp_avg_sq),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

void fused_adam_core(float* p, const float* g, float* m, float* v, 
                      int size, float step_size, float beta1, float beta2, 
                      float eps, float bias_corr2_sqrt) {
    int i = 0;
#if defined(__AVX2__)
    __m256 v_beta1 = _mm256_set1_ps(beta1);
    __m256 v_beta2 = _mm256_set1_ps(beta2);
    __m256 v_one_minus_beta1 = _mm256_set1_ps(1.0f - beta1);
    __m256 v_one_minus_beta2 = _mm256_set1_ps(1.0f - beta2);
    __m256 v_step_size = _mm256_set1_ps(step_size);
    __m256 v_eps = _mm256_set1_ps(eps);
    __m256 v_bc2_sqrt = _mm256_set1_ps(bias_corr2_sqrt);

    for (; i <= size - 8; i += 8) {
        __m256 grad = _mm256_loadu_ps(g + i);
        __m256 m_old = _mm256_loadu_ps(m + i);
        __m256 v_old = _mm256_loadu_ps(v + i);
        __m256 param = _mm256_loadu_ps(p + i);

        __m256 m_new = _mm256_add_ps(_mm256_mul_ps(v_beta1, m_old), _mm256_mul_ps(v_one_minus_beta1, grad));
        __m256 v_new = _mm256_add_ps(_mm256_mul_ps(v_beta2, v_old), _mm256_mul_ps(v_one_minus_beta2, _mm256_mul_ps(grad, grad)));
        
        _mm256_storeu_ps(m + i, m_new);
        _mm256_storeu_ps(v + i, v_new);

        __m256 denom = _mm256_add_ps(_mm256_sqrt_ps(_mm256_div_ps(v_new, v_bc2_sqrt)), v_eps);
        __m256 update = _mm256_mul_ps(v_step_size, _mm256_div_ps(m_new, denom));
        
        _mm256_storeu_ps(p + i, _mm256_sub_ps(param, update));
    }
#elif defined(__SSE2__)
    __m128 v_beta1 = _mm_set1_ps(beta1);
    __m128 v_beta2 = _mm_set1_ps(beta2);
    __m128 v_one_m_beta1 = _mm_set1_ps(1.0f - beta1);
    __m128 v_one_m_beta2 = _mm_set1_ps(1.0f - beta2);
    __m128 v_step_size = _mm_set1_ps(step_size);
    __m128 v_eps = _mm_set1_ps(eps);
    __m128 v_bc2 = _mm_set1_ps(bias_corr2_sqrt);

    for (; i <= size - 4; i += 4) {
        __m128 grad = _mm_loadu_ps(g + i);
        __m128 m_val = _mm_add_ps(_mm_mul_ps(v_beta1, _mm_loadu_ps(m + i)), _mm_mul_ps(v_one_m_beta1, grad));
        __m128 v_val = _mm_add_ps(_mm_mul_ps(v_beta2, _mm_loadu_ps(v + i)), _mm_mul_ps(v_one_m_beta2, _mm_mul_ps(grad, grad)));
        
        _mm_storeu_ps(m + i, m_val);
        _mm_storeu_ps(v + i, v_val);

        __m128 denom = _mm_add_ps(_mm_sqrt_ps(_mm_div_ps(v_val, v_bc2)), v_eps);
        __m128 p_val = _mm_sub_ps(_mm_loadu_ps(p + i), _mm_mul_ps(v_step_size, _mm_div_ps(m_val, denom)));
        _mm_storeu_ps(p + i, p_val);
    }
#endif
    for (; i < size; ++i) {
        m[i] = beta1 * m[i] + (1.0f - beta1) * g[i];
        v[i] = beta2 * v[i] + (1.0f - beta2) * g[i] * g[i];
        float denom = std::sqrt(v[i] / bias_corr2_sqrt) + eps;
        p[i] -= step_size * (m[i] / denom);
    }
}

/*[AutoHeaderGeneratorExport]*/
void offloaded_adam(
    TorchTensorView param,      // Device
    TorchTensorView grad,       // Device
    TorchTensorView exp_avg,    // Host
    TorchTensorView exp_avg_sq, // Host
    float lr, float beta1, float beta2, float eps, int step
) {
    const int64_t numel = _tv_numel(param);
    if (numel == 0) return;

    const float bias_correction1 = 1.0f - powf(beta1, step);
    const float bias_correction2 = 1.0f - powf(beta2, step);
    const float step_size = lr / bias_correction1;
    const float bc2_sqrt = bias_correction2;

    float* p_ptr = _tv_f(param);
    float* g_ptr = _tv_f(grad);
    float* m_ptr = _tv_f(exp_avg);
    float* v_ptr = _tv_f(exp_avg_sq);

    const int chunk_size = 1024 * 256; // ~1MB per buffer to stay in L3 cache
    cudaStream_t stream = (cudaStream_t)0;

    float *h_p, *h_g;
    cudaMallocHost(&h_p, chunk_size * sizeof(float) * 2); // double buffered
    cudaMallocHost(&h_g, chunk_size * sizeof(float) * 2);

    auto process_chunk = [&](int64_t offset, int size, int buf_idx) {
        float* curr_h_p = h_p + (buf_idx * chunk_size);
        float* curr_h_g = h_g + (buf_idx * chunk_size);

        cudaMemcpyAsync(curr_h_p, p_ptr + offset, size * sizeof(float), cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(curr_h_g, g_ptr + offset, size * sizeof(float), cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);

        fused_adam_core(curr_h_p, curr_h_g, m_ptr + offset, v_ptr + offset,
                         size, step_size, beta1, beta2, eps, bc2_sqrt);

        cudaMemcpyAsync(p_ptr + offset, curr_h_p, size * sizeof(float), cudaMemcpyHostToDevice, stream);
    };

    for (int64_t offset = 0; offset < numel; offset += chunk_size) {
        int current_chunk = std::min(chunk_size, (int)(numel - offset));
        process_chunk(offset, current_chunk, (int)offset % 2);
    }
    cudaStreamSynchronize(stream);

    cudaFreeHost(h_p);
    cudaFreeHost(h_g);
}

/*[AutoHeaderGeneratorExport]*/
void semi_offloaded_adam(
    TorchTensorView param,      // Device
    TorchTensorView grad,       // Device
    TorchTensorView exp_avg,    // Host
    TorchTensorView exp_avg_sq, // Host
    float lr, float beta1, float beta2, float eps, int step
) {
    const int64_t numel = _tv_numel(param);
    if (numel == 0) return;

    // 64MB VRAM; TODO: reduce this for smaller inputs
    constexpr int64_t chunk_size = 64 * 0x100000 / (2 * sizeof(float));

    float* p_ptr = _tv_f(param);
    float* g_ptr = _tv_f(grad);
    float* m_ptr_host = _tv_f(exp_avg);
    float* v_ptr_host = _tv_f(exp_avg_sq);

    float* buffer = DevicePool::global().acquire<float>(PoolSlot::SemiOffloadedAdamBuf, 2 * chunk_size);
    float* m_ptr_device = buffer;
    float* v_ptr_device = buffer + chunk_size;

    cudaStream_t stream = (cudaStream_t)0;

    const float bias_correction1 = 1.0f - powf(beta1, step);
    const float bias_correction2 = 1.0f - powf(beta2, step);

    auto process_chunk = [&](int64_t offset, size_t size) {

        cudaMemcpyAsync(m_ptr_device, m_ptr_host+offset, size * sizeof(float), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(v_ptr_device, v_ptr_host+offset, size * sizeof(float), cudaMemcpyHostToDevice, stream);
        cudaStreamSynchronize(stream);

        fused_adam_kernel<float><<<_LAUNCH_ARGS_1D(size, 256)>>>(
            p_ptr + offset,
            g_ptr + offset,
            m_ptr_device,
            v_ptr_device,
            lr,
            beta1,
            beta2,
            bias_correction1,
            bias_correction2,
            eps,
            size
        );

        cudaMemcpyAsync(m_ptr_host+offset, m_ptr_device, size * sizeof(float), cudaMemcpyDeviceToHost, stream);
        cudaMemcpyAsync(v_ptr_host+offset, v_ptr_device, size * sizeof(float), cudaMemcpyDeviceToHost, stream);
    };

    for (int64_t offset = 0; offset < numel; offset += chunk_size) {
        int64_t current_chunk = std::min(chunk_size, (numel - offset));
        process_chunk(offset, current_chunk);
    }
    cudaStreamSynchronize(stream);
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_multi(
    std::vector<TorchTensorView> params,
    std::vector<TorchTensorView> grads,
    std::vector<TorchTensorView> exp_avgs,
    std::vector<TorchTensorView> exp_avg_sqs,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
) {
    const int num_tensors = params.size();

    for (int i = 0; i < num_tensors; i++)
        fused_adam(
            params[i],
            grads[i],
            exp_avgs[i],
            exp_avg_sqs[i],
            lr,
            beta1,
            beta2,
            eps,
            step
        );
}



// ================
// Adam for quaternions
// ================

__global__ void fused_adam_riemannian_quat_kernel(
    float4* __restrict__ param,
    const float4* __restrict__ grad,
    float4* __restrict__ exp_avg,
    float4* __restrict__ exp_avg_sq,
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const int64_t numel
) {
    const int64_t qid = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    if (qid >= numel) return;

    float4 q = param[qid];
    float4 g = grad[qid];
    float4 m = exp_avg[qid];
    float4 v = exp_avg_sq[qid];

    float dot_qg = dot(q, g);
    g -= q * dot_qg;

    m = beta1 * m + (1.f - beta1) * g;
    v = beta2 * v + (1.f - beta2) * g*g;

    const float step_size = lr / bias_correction1;
    float4 d = -step_size * m / (sqrtf(v / bias_correction2) + eps);

    param[qid] = normalize(q + d);
    exp_avg[qid] = m;
    exp_avg_sq[qid] = v;
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_riemannian_quat(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
) {
    const int64_t numel = _tv_numel(param) / 4;
    if (numel == 0)
        return;

    fused_adam_riemannian_quat_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        (float4*)_tv_f(param),
        (float4*)_tv_f(grad),
        (float4*)_tv_f(exp_avg),
        (float4*)_tv_f(exp_avg_sq),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Newton
// ================

template <typename scalar_t>
__global__ void fused_newton_kernel(
    scalar_t* __restrict__ param,
    const scalar_t* __restrict__ grad,
    const scalar_t* __restrict__ hess_diag,
    scalar_t* __restrict__ exp_avg,
    scalar_t* __restrict__ exp_avg_sq,
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const int numel
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numel) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;
        
        // Load values
        const float grad_val = static_cast<float>(grad[idx]);
        const float hess_diag_val = static_cast<float>(hess_diag[idx]);
        float m_val = static_cast<float>(exp_avg[idx]);
        float v_val = static_cast<float>(exp_avg_sq[idx]);
        
        // Update biased first moment estimate: m_t = beta1 * m_{t-1} + (1 - beta1) * g_t
        m_val = beta1 * m_val + (1.0f - beta1) * grad_val;
        
        // Update biased second raw moment estimate: v_t = beta2 * v_{t-1} + (1 - beta2) * h_t
        v_val = beta2 * v_val + (1.0f - beta2) * hess_diag_val;
        
        // Compute update: theta_t = theta_{t-1} - step_size * m_t / (v_t / bias_correction2 + eps)
        const float denom = v_val / bias_correction2 + eps;
        const float param_val = static_cast<float>(param[idx]) - step_size * (m_val / denom);
        
        // Write back
        param[idx] = static_cast<scalar_t>(param_val);
        exp_avg[idx] = static_cast<scalar_t>(m_val);
        exp_avg_sq[idx] = static_cast<scalar_t>(v_val);
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_newton(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView hess_diag,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step1,
    int step2
) {
    const int numel = (int)_tv_numel(param);
    if (numel == 0)
        return;

    fused_newton_kernel<float><<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        _tv_f(param),
        _tv_f(grad),
        _tv_f(hess_diag),
        _tv_f(exp_avg),
        _tv_f(exp_avg_sq),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step1),
        1.0f - powf(beta2, step2),
        eps,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void fused_newton_multi(
    std::vector<TorchTensorView> params,
    std::vector<TorchTensorView> grads,
    std::vector<TorchTensorView> hess_diags,
    std::vector<TorchTensorView> exp_avgs,
    std::vector<TorchTensorView> exp_avg_sqs,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step1,
    int step2
) {
    const int num_tensors = params.size();

    for (int i = 0; i < num_tensors; i++)
        fused_newton(
            params[i],
            grads[i],
            hess_diags[i],
            exp_avgs[i],
            exp_avg_sqs[i],
            lr,
            beta1,
            beta2,
            eps,
            step1,
            step2
        );
}



// ================
// Adam Mean, scale agnostic
// ================

__global__ void fused_adam_scale_agnostic_mean_kernel(
    float3* __restrict__ param,
    const float3* __restrict__ grad,
    float3* __restrict__ exp_avg,
    float3* __restrict__ exp_avg_sq,
    const float3* __restrict__ scales,  // [N, 3], log space
    const float4* __restrict__ quats,  // [N, 4], unnormalized
    const float* __restrict__ opacities,  // [N, 1], logit space
    const float* __restrict__ radii,  // [N]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr,
    const int64_t numel
) {
    const int64_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    if (idx >= numel) return;

    float3 mean = param[idx];
    float3 v_mean = grad[idx];
    float3 g1_mean = exp_avg[idx];
    float3 g2_mean = exp_avg_sq[idx];

    float3 log_scale = scales[idx];
    float4 quat = quats[idx];
    float opac = opacities[idx];
    float3 sqrt_scale = float3{expf(0.5f*log_scale.x), expf(0.5f*log_scale.y), expf(0.5f*log_scale.z)};
    quat = normalize(quat);
    opac = sigmoid(opac);
    Matrix<float, 3, 3> sqrt_covar;
    SlangProjectionUtils::quat_scale_to_covar(quat, sqrt_scale, &sqrt_covar);
    Matrix<float, 3, 3> covar;
    SlangProjectionUtils::quat_scale_to_covar(quat, sqrt_scale*sqrt_scale, &covar);

    float3 v_mean_scaled_num = float3{dot(sqrt_covar[0], v_mean), dot(sqrt_covar[1], v_mean), dot(sqrt_covar[2], v_mean)};
    float v_mean_scaled_den = radii[idx] * 0.25f;
    #if 0
    v_mean_scaled_num = v_mean_scaled_num / (v_mean_scaled_den + (eps * (float)numel));
    v_mean_scaled_den = 1.0f;
    #endif
    g1_mean = beta1 * g1_mean + (1.f - beta1) * v_mean_scaled_num;
    g2_mean = beta2 * g2_mean + (1.f - beta2) * v_mean*v_mean * v_mean_scaled_den*v_mean_scaled_den;

    const float step_size = lr / bias_correction1;
    float3 delta = -step_size * g1_mean / (sqrtf(g2_mean / bias_correction2) + eps);

    // trust region clip
    #if 0
    float k = -8.0f * __logf(1.0f - eps_tr / fmaxf(opac, 1e-12f));
    float3 clip = {
        sqrtf(fmaxf(k * covar[0].x, 0.0f)),
        sqrtf(fmaxf(k * covar[1].y, 0.0f)),
        sqrtf(fmaxf(k * covar[2].z, 0.0f))
    };
    delta.x = fminf(fmaxf(delta.x, -clip.x), clip.x);
    delta.y = fminf(fmaxf(delta.y, -clip.y), clip.y);
    delta.z = fminf(fmaxf(delta.z, -clip.z), clip.z);
    delta.x = isfinite(delta.x) ? delta.x : 0.0f;
    delta.y = isfinite(delta.y) ? delta.y : 0.0f;
    delta.z = isfinite(delta.z) ? delta.z : 0.0f;
    #endif

    param[idx] = mean + delta;
    exp_avg[idx] = g1_mean;
    exp_avg_sq[idx] = g2_mean;
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_scale_agnostic_mean(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView scales,
    TorchTensorView quats,
    TorchTensorView opacities,
    TorchTensorView radii,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step
) {
    const int64_t numel = _tv_numel(param) / 3;
    if (numel == 0)
        return;

    fused_adam_scale_agnostic_mean_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        (float3*)_tv_f(param),
        (float3*)_tv_f(grad),
        (float3*)_tv_f(exp_avg),
        (float3*)_tv_f(exp_avg_sq),
        (float3*)_tv_f(scales),
        (float4*)_tv_f(quats),
        _tv_f(opacities),
        _tv_f(radii),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps,
        eps_tr,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Fused geometry optimizer
// ================

// Per-attribute non-SH Adam-state quantization helper. Identical math to the
// FPBO kernel's _NonShQ helper (kept inline here to avoid pulling FPBO
// templates into Optimizer.cu's TU). PRIMS is the per-splat cell count
// (means/scales/dc=3, quats=4, opacities=1).
template<int PRIMS>
struct _OptimNonShQ {
    using NQB = QuantizedAdamState<16, 256>;
    __device__ static inline void decode(
        const uint8_t* __restrict__ packed, const float4* __restrict__ bounds,
        int64_t gid, float* __restrict__ g1_out, float* __restrict__ g2_out
    ) {
        float4 mm_old = bounds[blockIdx.x];
        int64_t base = (int64_t)PRIMS * gid;
        #pragma unroll
        for (int i = 0; i < PRIMS; ++i) {
            float2 ab = NQB::decode_g1g2(packed, base + i, mm_old);
            g1_out[i] = ab.x;
            g2_out[i] = ab.y;
        }
    }
    __device__ static inline void accumulate(
        const float* __restrict__ g1, const float* __restrict__ g2, float4& mm
    ) {
        #pragma unroll
        for (int i = 0; i < PRIMS; ++i) {
            float2 us = NQB::g1g2_to_us(g1[i], g2[i]);
            mm.x = fminf(mm.x, us.x);
            mm.y = fmaxf(mm.y, us.x);
            mm.z = fminf(mm.z, us.y);
            mm.w = fmaxf(mm.w, us.y);
        }
    }
    __device__ static inline void encode(
        uint8_t* __restrict__ packed, int64_t gid,
        const float* __restrict__ g1, const float* __restrict__ g2, float4 mm
    ) {
        int64_t base = (int64_t)PRIMS * gid;
        #pragma unroll
        for (int i = 0; i < PRIMS; ++i) {
            NQB::encode_g1g2(packed, base + i, g1[i], g2[i], mm);
        }
    }
};

// Block reduction over a float4 (paired min/max for u and sqrt_g2).
template<int BLOCK_SIZE>
__device__ inline float4 _optim_block_reduce_minmax_f4(float4 mm) {
    cg::thread_block       block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    mm.x = cg::reduce(warp, mm.x, cg::less<float>());
    mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
    mm.z = cg::reduce(warp, mm.z, cg::less<float>());
    mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
    __shared__ float4 shared_reduce[BLOCK_SIZE / WARP_SIZE];
    if (threadIdx.x % WARP_SIZE == 0)
        shared_reduce[threadIdx.x / WARP_SIZE] = mm;
    __syncthreads();
    mm = (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
        ? shared_reduce[threadIdx.x]
        : float4{1e30f, -1e30f, 1e30f, -1e30f};
    mm.x = cg::reduce(warp, mm.x, cg::less<float>());
    mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
    mm.z = cg::reduce(warp, mm.z, cg::less<float>());
    mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
    __syncthreads();
    if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
        shared_reduce[threadIdx.x] = mm;
    __syncthreads();
    return shared_reduce[threadIdx.x / WARP_SIZE];
}


template<bool use_scale_agnostic_mean, bool zero_grad, bool non_sh_quant>
__global__ void fused_optim_3dgs_geometry_kernel(
    float3* __restrict__ means,
    float3* __restrict__ v_means,
    float3* __restrict__ g1_means,
    float3* __restrict__ g2_means,
    float4* __restrict__ quats,
    float4* __restrict__ v_quats,
    float4* __restrict__ g1_quats,
    float4* __restrict__ g2_quats,
    float3* __restrict__ scales,
    float3* __restrict__ v_scales,
    float3* __restrict__ g1_scales,
    float3* __restrict__ g2_scales,
    float* __restrict__ opacities,
    float* __restrict__ v_opacities,
    float* __restrict__ g1_opacities,
    float* __restrict__ g2_opacities,
    // features_dc -- only consumed when non_sh_quant is on; the non-quant path
    // continues to drive features_dc via the separate fused_adam_step call.
    float3* __restrict__ features_dc,
    float3* __restrict__ v_features_dc,
    const float* __restrict__ radii,
    // [N] optional densification score output: ||dL/dmean_world|| * max
    // post-exp world scale, written per splat when non-null (see
    // DensifyConfig::score_blend_world_grad).
    float* __restrict__ densify_score,
    float lr_means,
    float lr_quats,
    float lr_scales,
    float lr_opacs,
    float lr_features_dc,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const float sh_reg_weight,
    const float grad_scale,
    // Non-SH Adam-state quantization bundle. Only read when non_sh_quant is on.
    const NonShQuantState non_sh,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const int64_t numel
) {
    static constexpr float eps = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;
    [[maybe_unused]] static constexpr int   BLOCK_SIZE = 256;
    [[maybe_unused]] static constexpr float kSh0 = 0.28209479177387814f;  // 1 / (2*sqrt(pi))

    const int64_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    bool inside = (idx < numel);

    // Per-thread non-SH-quant state. Hoisted to function scope so the
    // post-update block reduction can see all attrs together; only used when
    // non_sh_quant is on.
    [[maybe_unused]] float3 nq_g1_scale = make_float3(0.f), nq_g2_scale = make_float3(0.f);
    [[maybe_unused]] float4 nq_g1_quat  = make_float4(0.f), nq_g2_quat  = make_float4(0.f);
    [[maybe_unused]] float  nq_g1_opac = 0.f, nq_g2_opac = 0.f;
    [[maybe_unused]] float3 nq_g1_mean = make_float3(0.f), nq_g2_mean = make_float3(0.f);
    [[maybe_unused]] float3 nq_g1_dc   = make_float3(0.f), nq_g2_dc   = make_float3(0.f);
    [[maybe_unused]] float4 nq_mm_scale = make_float4(1e30f, -1e30f, 1e30f, -1e30f);
    [[maybe_unused]] float4 nq_mm_quat  = make_float4(1e30f, -1e30f, 1e30f, -1e30f);
    [[maybe_unused]] float4 nq_mm_opac  = make_float4(1e30f, -1e30f, 1e30f, -1e30f);
    [[maybe_unused]] float4 nq_mm_mean  = make_float4(1e30f, -1e30f, 1e30f, -1e30f);
    [[maybe_unused]] float4 nq_mm_dc    = make_float4(1e30f, -1e30f, 1e30f, -1e30f);

    if (inside) {
        float step = (float)(steps != nullptr ? steps[idx] : scalar_step);
        float inv_bias_correction1 = 1.0f / (1.0f - powf(beta1, step));
        float inv_bias_correction2 = 1.0f / (1.0f - powf(beta2, step));
        lr_means *= inv_bias_correction1;
        lr_quats *= inv_bias_correction1;
        lr_scales *= inv_bias_correction1;
        lr_opacs *= inv_bias_correction1;
        lr_features_dc *= inv_bias_correction1;

        // load params
        float3 scale = scales[idx];
        float4 quat = normalize(quats[idx]);
        float opac = opacities[idx];

        // load gradient, with regularization
        static constexpr int kNumPerSplatLosses = 5;
        FixedArray<float, kNumPerSplatLosses> v_losses = {1.0f, 1.0f, 1.0f, 1.0f, 1.0f};
        float3 v_scale = make_float3(0.0f);
        float4 v_quat = make_float4(0.0f);
        float v_opac = 0.0f;
        SlangPerSplatLosses::per_splat_losses_bwd(
            scale, opac, quat, v_losses,
            &v_scale, &v_opac, &v_quat,
            mcmc_opacity_reg_weight,
            mcmc_scale_reg_weight,
            max_gauss_ratio,
            scale_regularization_weight,
            erank_reg_weight,
            erank_reg_weight_s3,
            quat_norm_reg_weight
        );
        v_scale += grad_scale * v_scales[idx];
        v_quat += grad_scale * v_quats[idx];
        v_opac += grad_scale * v_opacities[idx];
        if constexpr (zero_grad) {
            v_scales[idx] = make_float3(0.0f);
            v_quats[idx] = make_float4(0.0f);
            v_opacities[idx] = 0.0f;
        }

        // update scales
        float3 g1_scale, g2_scale;
        if constexpr (non_sh_quant) {
            _OptimNonShQ<3>::decode(non_sh.scales_packed, non_sh.scales_bounds, idx,
                                    (float*)&g1_scale, (float*)&g2_scale);
        } else {
            g1_scale = g1_scales[idx];
            g2_scale = g2_scales[idx];
        }
        g1_scale = beta1 * g1_scale + (1.f - beta1) * v_scale;
        g2_scale = beta2 * g2_scale + (1.f - beta2) * v_scale*v_scale;
        float3 updated_scale = scale - lr_scales * g1_scale / (sqrtf(g2_scale * inv_bias_correction2) + eps);
        scales[idx] = updated_scale;
        if constexpr (non_sh_quant) {
            nq_g1_scale = g1_scale;
            nq_g2_scale = g2_scale;
            _OptimNonShQ<3>::accumulate((float*)&g1_scale, (float*)&g2_scale, nq_mm_scale);
        } else {
            g1_scales[idx] = g1_scale;
            g2_scales[idx] = g2_scale;
        }

        // update quats (Riemannian)
        v_quat -= dot(quat, v_quat) * quat;
        float4 g1_quat, g2_quat;
        if constexpr (non_sh_quant) {
            _OptimNonShQ<4>::decode(non_sh.quats_packed, non_sh.quats_bounds, idx,
                                    (float*)&g1_quat, (float*)&g2_quat);
        } else {
            g1_quat = g1_quats[idx];
            g2_quat = g2_quats[idx];
        }
        g1_quat = beta1 * g1_quat + (1.f - beta1) * v_quat;
        g2_quat = beta2 * g2_quat + (1.f - beta2) * v_quat*v_quat;
        quats[idx] = normalize(quat - lr_quats * g1_quat / (sqrtf(g2_quat * inv_bias_correction2) + eps));
        if constexpr (non_sh_quant) {
            nq_g1_quat = g1_quat;
            nq_g2_quat = g2_quat;
            _OptimNonShQ<4>::accumulate((float*)&g1_quat, (float*)&g2_quat, nq_mm_quat);
        } else {
            g1_quats[idx] = g1_quat;
            g2_quats[idx] = g2_quat;
        }

        // update opacs
        float g1_opac, g2_opac;
        if constexpr (non_sh_quant) {
            _OptimNonShQ<1>::decode(non_sh.opacities_packed, non_sh.opacities_bounds, idx,
                                    &g1_opac, &g2_opac);
        } else {
            g1_opac = g1_opacities[idx];
            g2_opac = g2_opacities[idx];
        }
        g1_opac = beta1 * g1_opac + (1.f - beta1) * v_opac;
        g2_opac = beta2 * g2_opac + (1.f - beta2) * v_opac*v_opac;
        float updated_opac = opac - lr_opacs * g1_opac / (sqrtf(g2_opac * inv_bias_correction2) + eps);
        opacities[idx] = updated_opac;
        if constexpr (non_sh_quant) {
            nq_g1_opac = g1_opac;
            nq_g2_opac = g2_opac;
            _OptimNonShQ<1>::accumulate(&g1_opac, &g2_opac, nq_mm_opac);
        } else {
            g1_opacities[idx] = g1_opac;
            g2_opacities[idx] = g2_opac;
        }

        // update means (scale agnostic)
        float3 mean = means[idx];
        float3 v_mean = grad_scale * v_means[idx];
        if constexpr (zero_grad)
            v_means[idx] = make_float3(0.0f);
        // densification score: ||dL/dmean_world|| * max post-exp world scale.
        // v_mean is the (1/B-scaled) batch-accumulated data gradient; `scale`
        // holds the pre-update log scales loaded above. The grad_scale factor
        // is a global constant and doesn't affect relative ranking.
        if (densify_score != nullptr) {
            float s_max = fmaxf(fmaxf(scale.x, scale.y), scale.z);
            densify_score[idx] = expf(s_max) * sqrtf(dot(v_mean, v_mean));
        }
        float3 g1_mean, g2_mean;
        if constexpr (non_sh_quant) {
            _OptimNonShQ<3>::decode(non_sh.means_packed, non_sh.means_bounds, idx,
                                    (float*)&g1_mean, (float*)&g2_mean);
        } else {
            g1_mean = g1_means[idx];
            g2_mean = g2_means[idx];
        }

        float opac_post_sigmoid = sigmoid(opac);

        if constexpr (use_scale_agnostic_mean) {
            float3 v_mean_scaled_num = SlangProjectionUtils::apply_covar_to_vec(
                quat,
                {expf(0.5f*scale.x), expf(0.5f*scale.y), expf(0.5f*scale.z)},
                v_mean
            ) * sqrtf(2.0f * __logf(fmaxf(255.0f * opac_post_sigmoid, 1.00001f)));
            float v_mean_scaled_den = radii[idx] * 0.6f;
            g1_mean = beta1 * g1_mean + (1.f - beta1) * v_mean_scaled_num;
            g2_mean = beta2 * g2_mean + (1.f - beta2) * v_mean*v_mean * v_mean_scaled_den*v_mean_scaled_den;
        } else {
            g1_mean = beta1 * g1_mean + (1.f - beta1) * v_mean;
            g2_mean = beta2 * g2_mean + (1.f - beta2) * v_mean*v_mean;
        }

        means[idx] = mean - lr_means * g1_mean / (sqrtf(g2_mean * inv_bias_correction2) + eps);
        if constexpr (non_sh_quant) {
            nq_g1_mean = g1_mean;
            nq_g2_mean = g2_mean;
            _OptimNonShQ<3>::accumulate((float*)&g1_mean, (float*)&g2_mean, nq_mm_mean);
        } else {
            g1_means[idx] = g1_mean;
            g2_means[idx] = g2_mean;
        }

        // update features_dc (only when non_sh_quant is on; the fp32 path
        // keeps the separate fused_adam_step launch for features_dc).
        if constexpr (non_sh_quant) {
            float3 fdc  = features_dc[idx];
            float3 v_dc = grad_scale * v_features_dc[idx];
            if constexpr (zero_grad)
                v_features_dc[idx] = make_float3(0.0f);
            // L1-shrinkage regularization for SH/DC color: matches the
            // existing fused_adam_step(features_dc) launch's
            // l2_reg + l2_reg_offset = 0.5/kSh0 (mirror it here so non_sh_quant
            // doesn't silently drop the DC color reg). reg pushes toward 0
            // within [-0.5/kSh0, +0.5/kSh0] using a clamped grad.
            const float dc_off = 0.5f / kSh0;
            v_dc.x += sh_reg_weight * (fmaxf(fdc.x - dc_off, 0.f) + fminf(fdc.x + dc_off, 0.f));
            v_dc.y += sh_reg_weight * (fmaxf(fdc.y - dc_off, 0.f) + fminf(fdc.y + dc_off, 0.f));
            v_dc.z += sh_reg_weight * (fmaxf(fdc.z - dc_off, 0.f) + fminf(fdc.z + dc_off, 0.f));

            float3 g1_dc, g2_dc;
            _OptimNonShQ<3>::decode(non_sh.features_dc_packed, non_sh.features_dc_bounds, idx,
                                    (float*)&g1_dc, (float*)&g2_dc);
            g1_dc = beta1 * g1_dc + (1.f - beta1) * v_dc;
            g2_dc = beta2 * g2_dc + (1.f - beta2) * v_dc*v_dc;
            features_dc[idx] = fdc - lr_features_dc * g1_dc / (sqrtf(g2_dc * inv_bias_correction2) + eps);
            nq_g1_dc = g1_dc;
            nq_g2_dc = g2_dc;
            _OptimNonShQ<3>::accumulate((float*)&g1_dc, (float*)&g2_dc, nq_mm_dc);
        }
    }  // inside

    // Block-reduce per-attribute bounds (in u, sqrt_g2 basis) + encode.
    // All threads must participate in the reduces; outside-threads pass
    // sentinel min/max that don't perturb the result.
    if constexpr (non_sh_quant) {
        nq_mm_scale = _optim_block_reduce_minmax_f4<BLOCK_SIZE>(nq_mm_scale);
        if (threadIdx.x == 0) non_sh.scales_bounds[blockIdx.x] = nq_mm_scale;
        nq_mm_quat = _optim_block_reduce_minmax_f4<BLOCK_SIZE>(nq_mm_quat);
        if (threadIdx.x == 0) non_sh.quats_bounds[blockIdx.x] = nq_mm_quat;
        nq_mm_opac = _optim_block_reduce_minmax_f4<BLOCK_SIZE>(nq_mm_opac);
        if (threadIdx.x == 0) non_sh.opacities_bounds[blockIdx.x] = nq_mm_opac;
        nq_mm_mean = _optim_block_reduce_minmax_f4<BLOCK_SIZE>(nq_mm_mean);
        if (threadIdx.x == 0) non_sh.means_bounds[blockIdx.x] = nq_mm_mean;
        nq_mm_dc = _optim_block_reduce_minmax_f4<BLOCK_SIZE>(nq_mm_dc);
        if (threadIdx.x == 0) non_sh.features_dc_bounds[blockIdx.x] = nq_mm_dc;
        if (inside) {
            _OptimNonShQ<3>::encode(non_sh.scales_packed,      idx, (float*)&nq_g1_scale, (float*)&nq_g2_scale, nq_mm_scale);
            _OptimNonShQ<4>::encode(non_sh.quats_packed,       idx, (float*)&nq_g1_quat,  (float*)&nq_g2_quat,  nq_mm_quat);
            _OptimNonShQ<1>::encode(non_sh.opacities_packed,   idx, &nq_g1_opac,          &nq_g2_opac,          nq_mm_opac);
            _OptimNonShQ<3>::encode(non_sh.means_packed,       idx, (float*)&nq_g1_mean,  (float*)&nq_g2_mean,  nq_mm_mean);
            _OptimNonShQ<3>::encode(non_sh.features_dc_packed, idx, (float*)&nq_g1_dc,    (float*)&nq_g2_dc,    nq_mm_dc);
        }
    }
}


/*[AutoHeaderGeneratorExport]*/
void fused_optim_3dgs_geometry(
    int64_t num_splats,
    DeviceVector<float3> means, DeviceVector<float3> v_means, DeviceVector<float3> g1_means, DeviceVector<float3> g2_means,
    DeviceVector<float4> quats, DeviceVector<float4> v_quats, DeviceVector<float4> g1_quats, DeviceVector<float4> g2_quats,
    DeviceVector<float3> scales, DeviceVector<float3> v_scales, DeviceVector<float3> g1_scales, DeviceVector<float3> g2_scales,
    DeviceVector<float> opacities, DeviceVector<float> v_opacities, DeviceVector<float> g1_opacities, DeviceVector<float> g2_opacities,
    // features_dc + grad: consumed only when non_sh has its quant buffers
    // populated (via non_sh.enabled flagged in `non_sh`). The fp32 path
    // ignores them and runs the separate fused_adam_step features_dc call.
    DeviceVector<float3> features_dc, DeviceVector<float3> v_features_dc,
    DeviceVector<float> radii,
    // optional [N] densification world-grad score output; empty disables.
    DeviceVector<float> densify_score,
    const float lr_means, const float lr_quats, const float lr_scales, const float lr_opacs,
    const float lr_features_dc,
    const float max_gauss_ratio, const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight, const float mcmc_scale_reg_weight,
    const float erank_reg_weight, const float erank_reg_weight_s3, const float quat_norm_reg_weight,
    const float sh_reg_weight,
    bool use_scale_agnostic_mean,
    NonShQuantState non_sh,
    int32_t step, DeviceVector<int32_t> per_splat_steps,
    float grad_scale, bool zero_grad
) {
    if (num_splats == 0)
        return;

    // Dispatch over (use_scale_agnostic_mean, zero_grad, non_sh_quant) -> 8 instantiations.
    using KFn = void(*)(
        float3*, float3*, float3*, float3*,
        float4*, float4*, float4*, float4*,
        float3*, float3*, float3*, float3*,
        float*,  float*,  float*,  float*,
        float3*, float3*,
        const float*,
        float*,
        float, float, float, float, float,
        const float, const float, const float, const float,
        const float, const float, const float, const float,
        const float,
        const NonShQuantState,
        const int32_t, const int32_t*, const int64_t);
    KFn kfn = nullptr;
    const bool nq = non_sh.enabled;
    if (use_scale_agnostic_mean) {
        if (zero_grad) {
            kfn = nq ? fused_optim_3dgs_geometry_kernel<true,  true,  true>
                     : fused_optim_3dgs_geometry_kernel<true,  true,  false>;
        } else {
            kfn = nq ? fused_optim_3dgs_geometry_kernel<true,  false, true>
                     : fused_optim_3dgs_geometry_kernel<true,  false, false>;
        }
    } else {
        if (zero_grad) {
            kfn = nq ? fused_optim_3dgs_geometry_kernel<false, true,  true>
                     : fused_optim_3dgs_geometry_kernel<false, true,  false>;
        } else {
            kfn = nq ? fused_optim_3dgs_geometry_kernel<false, false, true>
                     : fused_optim_3dgs_geometry_kernel<false, false, false>;
        }
    }
    kfn<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        means.data_ptr(), v_means.data_ptr(), g1_means.data_ptr(), g2_means.data_ptr(),
        quats.data_ptr(), v_quats.data_ptr(), g1_quats.data_ptr(), g2_quats.data_ptr(),
        scales.data_ptr(), v_scales.data_ptr(), g1_scales.data_ptr(), g2_scales.data_ptr(),
        opacities.data_ptr(), v_opacities.data_ptr(), g1_opacities.data_ptr(), g2_opacities.data_ptr(),
        features_dc.data_ptr(), v_features_dc.data_ptr(),
        (const float*)radii.data_ptr(),
        densify_score.data_ptr(),
        lr_means, lr_quats, lr_scales, lr_opacs, lr_features_dc,
        max_gauss_ratio,
        scale_regularization_weight / (float)num_splats,
        mcmc_opacity_reg_weight / (float)num_splats,
        mcmc_scale_reg_weight / (float)num_splats,
        erank_reg_weight / (float)num_splats,
        erank_reg_weight_s3 / (float)num_splats,
        quat_norm_reg_weight / (float)num_splats,
        sh_reg_weight,
        grad_scale,
        non_sh,
        step, per_splat_steps.data_ptr(), num_splats
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Fused appearance optimizer
// ================

template<bool zero_grad>
__global__ void fused_adam_with_steps_kernel(
    float* __restrict__ param,
    float* __restrict__ grad,
    float* __restrict__ exp_avg,
    float* __restrict__ exp_avg_sq,
    const float lr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const float decay,
    const float decay_offset,
    const float grad_scale,
    const int64_t numel,
    const int stride
) {
    static constexpr float eps = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;

    const int64_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    if (idx >= numel) return;

    float step = (float)(steps ? steps[idx / stride] : scalar_step);
    float inv_bias_correction1 = 1.0f / (1.0f - powf(beta1, step));
    float inv_bias_correction2 = 1.0f / (1.0f - powf(beta2, step));

    float x = param[idx];
    float v = grad[idx];
    if (!isfinite(v))
        v = 0.0f;
    if constexpr (zero_grad)
        grad[idx] = 0.0f;
    v *= grad_scale;
    v += decay * (fmaxf(x - decay_offset, 0.0f) + fminf(x + decay_offset, 0.0f));
    float g1 = exp_avg[idx];
    float g2 = exp_avg_sq[idx];

    g1 = beta1 * g1 + (1.0f - beta1) * v;
    g2 = beta2 * g2 + (1.0f - beta2) * v * v;

    x -= lr * inv_bias_correction1 * g1 / (sqrtf(g2 * inv_bias_correction2) + eps);

    param[idx] = x;
    exp_avg[idx] = g1;
    exp_avg_sq[idx] = g2;
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_step(
    int64_t num_splats,
    DeviceTensorFloatND param,
    DeviceTensorFloatND grad,
    DeviceTensorFloatND exp_avg,
    DeviceTensorFloatND exp_avg_sq,
    float lr,
    int32_t step, DeviceVector<int32_t> per_splat_steps,
    float l2_reg,
    float l2_reg_offset,
    float grad_scale, bool zero_grad
) {
    int64_t param_numel = param.numel();
    if (param_numel == 0 || num_splats == 0)
        return;
    int stride = (int)(param_numel / num_splats);

    auto kfn = zero_grad ? fused_adam_with_steps_kernel<true>
                         : fused_adam_with_steps_kernel<false>;
    kfn<<<_LAUNCH_ARGS_1D(num_splats*stride, 256)>>>(
        param.data_ptr(),
        grad.data_ptr(),
        exp_avg.data_ptr(),
        exp_avg_sq.data_ptr(),
        lr,
        step, per_splat_steps.data_ptr(),
        2.0f*l2_reg/(float)(num_splats*stride),
        l2_reg_offset,
        grad_scale,
        num_splats*stride,
        stride
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// Unscheduled fp32 AdaGrad:
//     accum += g*g
//     x    -= lr * g / (sqrt(accum) + eps)
// Matches torch.optim.Adagrad(lr_decay=0, weight_decay=0,
// initial_accumulator_value=0, eps=1e-15).
__global__ void fused_adagrad_kernel(
    float* __restrict__ param,
    const float* __restrict__ grad,
    float* __restrict__ accum,
    const float lr,
    const int64_t numel
) {
    static constexpr float eps = 1e-15f;
    const int64_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    if (idx >= numel) return;
    float v = grad[idx];
    if (!isfinite(v)) v = 0.0f;
    float a = accum[idx] + v * v;
    param[idx] -= lr * v / (sqrtf(a) + eps);
    accum[idx] = a;
}

/*[AutoHeaderGeneratorExport]*/
void fused_adagrad_step(
    DeviceTensorFloatND param,
    DeviceTensorFloatND grad,
    DeviceTensorFloatND accum,
    float lr
) {
    int64_t numel = param.numel();
    if (numel == 0) return;
    fused_adagrad_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        param.data_ptr(),
        grad.data_ptr(),
        accum.data_ptr(),
        lr,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template<int BLOCK_SIZE, int QUANT_BITS, bool zero_grad>
__global__ void fused_adam_with_steps_8bit_kernel(
    float* __restrict__ param,
    float* __restrict__ grad,
    uint8_t* __restrict__ packed,       // AoS (u, sqrt_g2) packed cells
    float4* __restrict__ quant_bounds,  // (u_min, u_max, sqrt_g2_min, sqrt_g2_max)
    const float lr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const float decay,
    const float decay_offset,
    const float grad_scale,
    const int64_t numel,
    const int stride
) {
    using QState = QuantizedAdamState<QUANT_BITS, BLOCK_SIZE>;
    static constexpr float eps = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;

    const int64_t idx = (int64_t)blockIdx.x * BLOCK_SIZE + (int64_t)threadIdx.x;
    bool inside = (idx < numel);

    float step = (float)(steps ? steps[idx / stride] : scalar_step);
    float inv_bias_correction1 = 1.0f / (1.0f - powf(beta1, step));
    float inv_bias_correction2 = 1.0f / (1.0f - powf(beta2, step));

    float x = inside ? param[idx] : 0.0f;
    float v = inside ? grad[idx] : 0.0f;
    if (!isfinite(v))
        v = 0.0f;
    if constexpr (zero_grad) {
        if (inside) grad[idx] = 0.0f;
    }
    v *= grad_scale;
    v += decay * (fmaxf(x - decay_offset, 0.0f) + fminf(x + decay_offset, 0.0f));

    // Decode joint (u, sqrt(g2)) via QuantizedAdamState codec.
    float4 mm = quant_bounds[blockIdx.x];
    float g1, g2;
    if (inside) {
        float2 g1g2 = QState::decode_g1g2(packed, idx, mm);
        g1 = g1g2.x;
        g2 = g1g2.y;
    } else {
        g1 = 0.0f;
        g2 = 0.0f;
    }

    g1 = beta1 * g1 + (1.0f - beta1) * v;
    g2 = beta2 * g2 + (1.0f - beta2) * v*v;

    x -= lr * inv_bias_correction1 * g1 / (sqrtf(g2 * inv_bias_correction2) + eps);
    param[idx] = x;

    // Re-encode the new Adam state in the (u, sqrt(g2)) basis.
    float2 us_new = QState::g1g2_to_us(g1, g2);
    float u_new       = us_new.x;
    float sqrt_g2_new = us_new.y;

    cg::thread_block block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    mm = inside ? float4{u_new, u_new, sqrt_g2_new, sqrt_g2_new}
                : float4{1e30f, -1e30f, 1e30f, -1e30f};
    mm.x = cg::reduce(warp, mm.x, cg::less<float>());
    mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
    mm.z = cg::reduce(warp, mm.z, cg::less<float>());
    mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
    __shared__ float4 shared_reduce[BLOCK_SIZE / WARP_SIZE];
    if (threadIdx.x % WARP_SIZE == 0)
        shared_reduce[threadIdx.x / WARP_SIZE] = mm;
    __syncthreads();
    mm = (threadIdx.x < BLOCK_SIZE / WARP_SIZE) ?
        shared_reduce[threadIdx.x] : float4{1e30f, -1e30f, 1e30f, -1e30f};
    mm.x = cg::reduce(warp, mm.x, cg::less<float>());
    mm.y = cg::reduce(warp, mm.y, cg::greater<float>());
    mm.z = cg::reduce(warp, mm.z, cg::less<float>());
    mm.w = cg::reduce(warp, mm.w, cg::greater<float>());
    __syncthreads();
    if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
        shared_reduce[threadIdx.x] = mm;
    __syncthreads();
    mm = shared_reduce[threadIdx.x / WARP_SIZE];

    if (inside) {
        QState::encode_us(packed, idx, u_new, sqrt_g2_new, mm);
    }

    if (threadIdx.x == 0)
        quant_bounds[blockIdx.x] = mm;
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_step_quantized(
    int64_t num_splats,
    DeviceTensorFloatND param,
    DeviceTensorFloatND grad,
    uint8_t* packed,                    // AoS (u, sqrt_g2) packed
    float4* quant_bounds,
    float lr,
    int32_t step, DeviceVector<int32_t> per_splat_steps,
    float l2_reg,
    float l2_reg_offset,
    int bits,                           // 4 or 8 -- selects QuantizedAdamState<BITS, 256>
    float grad_scale, bool zero_grad
) {
    int64_t param_numel = param.numel();
    if (param_numel == 0 || num_splats == 0)
        return;
    int stride = (int)(param_numel / num_splats);
    constexpr int BLOCK_SIZE = 256;

    #define _ARGS_TAIL \
        param.data_ptr(), grad.data_ptr(), packed, quant_bounds, \
        lr, step, per_splat_steps.data_ptr(), \
        2.0f*l2_reg/(float)(num_splats*stride), l2_reg_offset, \
        grad_scale, \
        num_splats*stride, stride

    if (bits == 4) {
        auto kfn = zero_grad ? fused_adam_with_steps_8bit_kernel<BLOCK_SIZE, 4, true>
                             : fused_adam_with_steps_8bit_kernel<BLOCK_SIZE, 4, false>;
        kfn<<<_LAUNCH_ARGS_1D(num_splats*stride, BLOCK_SIZE)>>>(_ARGS_TAIL);
    } else if (bits == 8) {
        auto kfn = zero_grad ? fused_adam_with_steps_8bit_kernel<BLOCK_SIZE, 8, true>
                             : fused_adam_with_steps_8bit_kernel<BLOCK_SIZE, 8, false>;
        kfn<<<_LAUNCH_ARGS_1D(num_splats*stride, BLOCK_SIZE)>>>(_ARGS_TAIL);
    } else {
        throw std::runtime_error(
            "fused_adam_step_quantized: bits must be 4 or 8, got " +
            std::to_string(bits));
    }
    #undef _ARGS_TAIL
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Doubly-quantized fused Adam for SH features (non-FPBO).
// Quantizes BOTH the optimizer state (g1, g2) AND the parameter value (sh
// coeffs) against per-CELL-block (256 cells per block) min/max bounds:
//
//   * Optimizer state -- QuantizedAdamState<OPTIM_BITS, 256> in (u, sqrt_g2)
//     basis, one float4 bound per block (u_min, u_max, sqrt_g2_min, sqrt_g2_max).
//   * Parameter value -- QuantizedTensor<VALUE_BITS, 256>, one float2 bound
//     per block (v_min, v_max).
//
// FPBO already handles this combo per-splat-block; this kernel mirrors the
// math but indexes per-cell-block, matching the existing non-FPBO Adam
// kernel's launch geometry. Used by engine_optim_step when sh_value_bits != 32
// and use_fused_proj_bwd_optim is off (no fp32 features_sh exists then;
// the packed buffer IS the canonical storage).
// ================

template<int BLOCK_SIZE, int OPTIM_BITS, int VALUE_BITS, bool zero_grad>
__global__ void fused_adam_with_steps_qq_kernel(
    float* __restrict__ grad,           // fp32 grad (still dense)
    uint8_t* __restrict__ optim_packed, // optim-state packed (u, sqrt_g2)
    float4*  __restrict__ optim_bounds, // per-cell-block float4 (u_mm, sqrt_g2_mm)
    uint8_t* __restrict__ value_packed, // parameter-value packed
    float2*  __restrict__ value_bounds, // per-cell-block float2 (v_min, v_max)
    const float lr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const float decay,
    const float decay_offset,
    const float grad_scale,
    const int64_t numel,
    const int stride
) {
    using QState     = QuantizedAdamState<OPTIM_BITS, BLOCK_SIZE>;
    using QValue     = QuantizedTensor<VALUE_BITS, BLOCK_SIZE>;
    static constexpr float eps   = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;

    const int64_t idx = (int64_t)blockIdx.x * BLOCK_SIZE + (int64_t)threadIdx.x;
    bool inside = (idx < numel);

    float step = (float)(steps ? steps[idx / stride] : scalar_step);
    float inv_bias_correction1 = 1.0f / (1.0f - powf(beta1, step));
    float inv_bias_correction2 = 1.0f / (1.0f - powf(beta2, step));

    // Decode old value from packed bytes against the block's old bound. Skip
    // the load for out-of-range threads (they pass sentinels through the
    // reduction so the bound only reflects inside threads' values).
    float2 value_mm_old = inside ? value_bounds[blockIdx.x] : float2{0.f, 0.f};
    float x = inside ? QValue::decode_v(value_packed, idx, value_mm_old) : 0.0f;

    float v = inside ? grad[idx] : 0.0f;
    if (!isfinite(v))
        v = 0.0f;
    if constexpr (zero_grad) {
        if (inside) grad[idx] = 0.0f;
    }
    v *= grad_scale;
    v += decay * (fmaxf(x - decay_offset, 0.0f) + fminf(x + decay_offset, 0.0f));

    // Decode joint (u, sqrt_g2) for Adam state.
    float4 optim_mm = optim_bounds[blockIdx.x];
    float g1 = 0.0f, g2 = 0.0f;
    if (inside) {
        float2 g1g2 = QState::decode_g1g2(optim_packed, idx, optim_mm);
        g1 = g1g2.x;
        g2 = g1g2.y;
    }

    g1 = beta1 * g1 + (1.0f - beta1) * v;
    g2 = beta2 * g2 + (1.0f - beta2) * v * v;

    x -= lr * inv_bias_correction1 * g1 / (sqrtf(g2 * inv_bias_correction2) + eps);

    // ---- value bounds block-reduce (float2 min/max) ----
    float2 value_mm_new = inside ? float2{x, x} : float2{1e30f, -1e30f};

    cg::thread_block       block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    value_mm_new.x = cg::reduce(warp, value_mm_new.x, cg::less<float>());
    value_mm_new.y = cg::reduce(warp, value_mm_new.y, cg::greater<float>());
    __shared__ float2 shared_reduce_v[BLOCK_SIZE / WARP_SIZE];
    if (threadIdx.x % WARP_SIZE == 0)
        shared_reduce_v[threadIdx.x / WARP_SIZE] = value_mm_new;
    __syncthreads();
    value_mm_new = (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
        ? shared_reduce_v[threadIdx.x] : float2{1e30f, -1e30f};
    value_mm_new.x = cg::reduce(warp, value_mm_new.x, cg::less<float>());
    value_mm_new.y = cg::reduce(warp, value_mm_new.y, cg::greater<float>());
    __syncthreads();
    if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
        shared_reduce_v[threadIdx.x] = value_mm_new;
    __syncthreads();
    value_mm_new = shared_reduce_v[threadIdx.x / WARP_SIZE];

    // ---- optim state bounds block-reduce ((u, sqrt_g2) -> float4) ----
    float2 us_new = QState::g1g2_to_us(g1, g2);
    float u_new = us_new.x, sqrt_g2_new = us_new.y;
    float4 optim_mm_new = inside
        ? float4{u_new, u_new, sqrt_g2_new, sqrt_g2_new}
        : float4{1e30f, -1e30f, 1e30f, -1e30f};
    optim_mm_new.x = cg::reduce(warp, optim_mm_new.x, cg::less<float>());
    optim_mm_new.y = cg::reduce(warp, optim_mm_new.y, cg::greater<float>());
    optim_mm_new.z = cg::reduce(warp, optim_mm_new.z, cg::less<float>());
    optim_mm_new.w = cg::reduce(warp, optim_mm_new.w, cg::greater<float>());
    __shared__ float4 shared_reduce_o[BLOCK_SIZE / WARP_SIZE];
    if (threadIdx.x % WARP_SIZE == 0)
        shared_reduce_o[threadIdx.x / WARP_SIZE] = optim_mm_new;
    __syncthreads();
    optim_mm_new = (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
        ? shared_reduce_o[threadIdx.x]
        : float4{1e30f, -1e30f, 1e30f, -1e30f};
    optim_mm_new.x = cg::reduce(warp, optim_mm_new.x, cg::less<float>());
    optim_mm_new.y = cg::reduce(warp, optim_mm_new.y, cg::greater<float>());
    optim_mm_new.z = cg::reduce(warp, optim_mm_new.z, cg::less<float>());
    optim_mm_new.w = cg::reduce(warp, optim_mm_new.w, cg::greater<float>());
    __syncthreads();
    if (threadIdx.x < BLOCK_SIZE / WARP_SIZE)
        shared_reduce_o[threadIdx.x] = optim_mm_new;
    __syncthreads();
    optim_mm_new = shared_reduce_o[threadIdx.x / WARP_SIZE];

    // Encode + write back. Each inside thread encodes its cell; thread 0 writes
    // the per-block bounds.
    if (inside) {
        QValue::encode_v(value_packed, idx, x, value_mm_new);
        QState::encode_us(optim_packed, idx, u_new, sqrt_g2_new, optim_mm_new);
    }
    if (threadIdx.x == 0) {
        value_bounds[blockIdx.x] = value_mm_new;
        optim_bounds[blockIdx.x] = optim_mm_new;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_step_quantized_value(
    int64_t num_splats,
    int64_t param_numel,                // = num_splats * stride (e.g. 3 * num_sh)
    DeviceTensorFloatND grad,           // fp32 dense [num_splats, stride]
    uint8_t* optim_packed,
    float4*  optim_bounds,
    uint8_t* value_packed,
    float2*  value_bounds,
    float lr,
    int32_t step, DeviceVector<int32_t> per_splat_steps,
    float l2_reg,
    float l2_reg_offset,
    int optim_bits,                     // 4 or 8
    int value_bits,                     // 8 or 16
    float grad_scale, bool zero_grad
) {
    if (param_numel == 0 || num_splats == 0)
        return;
    int stride = (int)(param_numel / num_splats);
    constexpr int BLOCK_SIZE = 256;

    #define _ARGS_TAIL_QQ \
        grad.data_ptr(), optim_packed, optim_bounds, value_packed, value_bounds, \
        lr, step, per_splat_steps.data_ptr(), \
        2.0f*l2_reg/(float)(num_splats*stride), l2_reg_offset, \
        grad_scale, \
        num_splats*stride, stride

    // 4 instantiations: (optim_bits in {4, 8}) x (value_bits in {8, 16}) x
    // (zero_grad in {false, true}). Throw on any other combination.
    if (optim_bits == 4 && value_bits == 8) {
        auto kfn = zero_grad ? fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 4, 8,  true>
                             : fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 4, 8,  false>;
        kfn<<<_LAUNCH_ARGS_1D(num_splats*stride, BLOCK_SIZE)>>>(_ARGS_TAIL_QQ);
    } else if (optim_bits == 4 && value_bits == 16) {
        auto kfn = zero_grad ? fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 4, 16, true>
                             : fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 4, 16, false>;
        kfn<<<_LAUNCH_ARGS_1D(num_splats*stride, BLOCK_SIZE)>>>(_ARGS_TAIL_QQ);
    } else if (optim_bits == 8 && value_bits == 8) {
        auto kfn = zero_grad ? fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 8, 8,  true>
                             : fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 8, 8,  false>;
        kfn<<<_LAUNCH_ARGS_1D(num_splats*stride, BLOCK_SIZE)>>>(_ARGS_TAIL_QQ);
    } else if (optim_bits == 8 && value_bits == 16) {
        auto kfn = zero_grad ? fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 8, 16, true>
                             : fused_adam_with_steps_qq_kernel<BLOCK_SIZE, 8, 16, false>;
        kfn<<<_LAUNCH_ARGS_1D(num_splats*stride, BLOCK_SIZE)>>>(_ARGS_TAIL_QQ);
    } else {
        throw std::runtime_error(
            "fused_adam_step_quantized_value: optim_bits in {4, 8} and "
            "value_bits in {8, 16}; got optim_bits=" + std::to_string(optim_bits) +
            ", value_bits=" + std::to_string(value_bits));
    }
    #undef _ARGS_TAIL_QQ
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// 3DGS^2-TR Mean - https://arxiv.org/abs/2602.00395
// ================

__global__ void fused_3dgs2tr_mean_optim_kernel(
    const int num_gs,
    float3* __restrict__ means,  // [N, 3]
    const float3* __restrict__ vr_means,  // [N, 3]
    const float3* __restrict__ h_means,  // [N, 3]
    const float3* __restrict__ scales,  // [N, 3], log space
    const float4* __restrict__ quats,  // [N, 4], unnormalized
    const float* __restrict__ opacities,  // [N, 1], logit space
    float3* __restrict__ exp_avg_means,  // [N, 3]
    float3* __restrict__ exp_avg_sq_means,  // [N, 3]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_gs) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;

        // Update momentum
        float3 vr_mean = vr_means[idx];
        float3 h_mean = h_means[idx];
        float3 exp_avg_mean = exp_avg_means[idx];
        float3 exp_avg_sq_mean = exp_avg_sq_means[idx];
        exp_avg_mean = beta1 * exp_avg_mean + (1.0f - beta1) * vr_mean;
        exp_avg_sq_mean = beta2 * exp_avg_sq_mean + (1.0f - beta2) * h_mean;
        exp_avg_means[idx] = exp_avg_mean;
        exp_avg_sq_means[idx] = exp_avg_sq_mean;

        // Compute delta
        float3 delta = -step_size * exp_avg_mean / (exp_avg_sq_mean / bias_correction2 + eps);

        // Compute trust region
        float3 scale = scales[idx];
        float4 quat = quats[idx];
        float opac = opacities[idx];
        scale = {expf(scale.x), expf(scale.y), expf(scale.z)};
        quat = normalize(quat);
        opac = sigmoid(opac);
        Matrix<float, 3, 3> covar;
        SlangProjectionUtils::quat_scale_to_covar(quat, scale, &covar);
        float k = -8.0f * __logf(1.0f - eps_tr / fmaxf(opac, 1e-12f));
        float3 clip = {
            sqrtf(fmaxf(k * covar[0].x, 0.0f)),
            sqrtf(fmaxf(k * covar[1].y, 0.0f)),
            sqrtf(fmaxf(k * covar[2].z, 0.0f))
        };

        // clip and update
        delta.x = fminf(fmaxf(delta.x, -clip.x), clip.x);
        delta.y = fminf(fmaxf(delta.y, -clip.y), clip.y);
        delta.z = fminf(fmaxf(delta.z, -clip.z), clip.z);
        delta.x = isfinite(delta.x) ? delta.x : 0.0f;
        delta.y = isfinite(delta.y) ? delta.y : 0.0f;
        delta.z = isfinite(delta.z) ? delta.z : 0.0f;
        means[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_3dgs2tr_mean_optim(
    TorchTensorView means,
    TorchTensorView vr_means,
    TorchTensorView h_means,
    TorchTensorView scales,
    TorchTensorView quats,
    TorchTensorView opacities,
    TorchTensorView exp_avg_means,
    TorchTensorView exp_avg_sq_means,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    const int num_gs = (int)(_tv_numel(means) / 3);
    if (num_gs == 0)
        return;

    fused_3dgs2tr_mean_optim_kernel<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        (float3*)_tv_f(means),
        (float3*)_tv_f(vr_means),
        (float3*)_tv_f(h_means),
        (float3*)_tv_f(scales),
        (float4*)_tv_f(quats),
        _tv_f(opacities),
        (float3*)_tv_f(exp_avg_means),
        (float3*)_tv_f(exp_avg_sq_means),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step1),
        1.0f - powf(beta2, step2),
        eps,
        eps_tr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// 3DGS^2-TR Scale - https://arxiv.org/abs/2602.00395
// ================

__global__ void fused_3dgs2tr_scale_optim_kernel(
    const int num_gs,
    float3* __restrict__ scales,  // [N, 3]
    const float3* __restrict__ vr_scales,  // [N, 3]
    const float3* __restrict__ h_scales,  // [N, 3]
    const float* __restrict__ opacities,  // [N, 1], logit space
    float3* __restrict__ exp_avg_scales,  // [N, 3]
    float3* __restrict__ exp_avg_sq_scales,  // [N, 3]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_gs) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;

        // Update momentum
        float3 vr_scale = vr_scales[idx];
        float3 h_scale = h_scales[idx];
        float3 exp_avg_scale = exp_avg_scales[idx];
        float3 exp_avg_sq_scale = exp_avg_sq_scales[idx];
        exp_avg_scale = beta1 * exp_avg_scale + (1.0f - beta1) * vr_scale;
        exp_avg_sq_scale = beta2 * exp_avg_sq_scale + (1.0f - beta2) * h_scale;
        exp_avg_scales[idx] = exp_avg_scale;
        exp_avg_sq_scales[idx] = exp_avg_sq_scale;

        // Compute delta
        float3 delta = -step_size * exp_avg_scale / (exp_avg_sq_scale / bias_correction2 + eps);

        // Compute trust region
        float opac = opacities[idx];
        opac = sigmoid(opac);
        float clip = sqrtf(2.0f * eps_tr / fmaxf(opac, 1e-12f));

        // clip and update
        delta.x = fminf(fmaxf(delta.x, -clip), clip);
        delta.y = fminf(fmaxf(delta.y, -clip), clip);
        delta.z = fminf(fmaxf(delta.z, -clip), clip);
        delta.x = isfinite(delta.x) ? delta.x : 0.0f;
        delta.y = isfinite(delta.y) ? delta.y : 0.0f;
        delta.z = isfinite(delta.z) ? delta.z : 0.0f;
        scales[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_3dgs2tr_scale_optim(
    TorchTensorView scales,
    TorchTensorView vr_scales,
    TorchTensorView h_scales,
    TorchTensorView opacities,
    TorchTensorView exp_avg_scales,
    TorchTensorView exp_avg_sq_scales,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    const int num_gs = (int)(_tv_numel(scales) / 3);
    if (num_gs == 0)
        return;

    fused_3dgs2tr_scale_optim_kernel<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        (float3*)_tv_f(scales),
        (float3*)_tv_f(vr_scales),
        (float3*)_tv_f(h_scales),
        _tv_f(opacities),
        (float3*)_tv_f(exp_avg_scales),
        (float3*)_tv_f(exp_avg_sq_scales),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step1),
        1.0f - powf(beta2, step2),
        eps,
        eps_tr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// 3DGS^2-TR Color - https://arxiv.org/abs/2602.00395
// ================

__global__ void fused_3dgs2tr_color_optim_kernel(
    const int num_gs,
    float3* __restrict__ colors,  // [N, 3]
    const float3* __restrict__ vr_colors,  // [N, 3]
    const float3* __restrict__ h_colors,  // [N, 3]
    const float* __restrict__ opacities,  // [N, 1], logit space
    float3* __restrict__ exp_avg_colors,  // [N, 3]
    float3* __restrict__ exp_avg_sq_colors,  // [N, 3]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_gs) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;

        // Update momentum
        float3 vr_color = vr_colors[idx];
        float3 h_color = h_colors[idx];
        float3 exp_avg_color = exp_avg_colors[idx];
        float3 exp_avg_sq_color = exp_avg_sq_colors[idx];
        exp_avg_color = beta1 * exp_avg_color + (1.0f - beta1) * vr_color;
        exp_avg_sq_color = beta2 * exp_avg_sq_color + (1.0f - beta2) * h_color;
        exp_avg_colors[idx] = exp_avg_color;
        exp_avg_sq_colors[idx] = exp_avg_sq_color;

        // Compute delta
        float3 delta = -step_size * exp_avg_color / (exp_avg_sq_color / bias_correction2 + eps);

        // Compute trust region
        float opac = opacities[idx];
        opac = sigmoid(opac);
        // float3 color = fmaxf(kSh0 * colors[idx] + 0.5f, (0.5f/255.0f)*(0.5f/255.0f));
        // float3 color = fmaxf(kSh0 * colors[idx] + 0.5f, 0.5f/255.0f);
        float3 color = fmaxf(kSh0 * colors[idx] + 0.5f, (1.0f/255.0f)*(1.0f/255.0f));
        float3 clip = kSh0 * sqrtf(4.0f * eps_tr * color / fmaxf(opac, 1e-12f));

        // clip and update
        delta.x = fminf(fmaxf(delta.x, -clip.x), clip.x);
        delta.y = fminf(fmaxf(delta.y, -clip.y), clip.y);
        delta.z = fminf(fmaxf(delta.z, -clip.z), clip.z);
        delta.x = isfinite(delta.x) ? delta.x : 0.0f;
        delta.y = isfinite(delta.y) ? delta.y : 0.0f;
        delta.z = isfinite(delta.z) ? delta.z : 0.0f;
        colors[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_3dgs2tr_color_optim(
    TorchTensorView colors,
    TorchTensorView vr_colors,
    TorchTensorView h_colors,
    TorchTensorView opacities,
    TorchTensorView exp_avg_colors,
    TorchTensorView exp_avg_sq_colors,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    const int num_gs = (int)(_tv_numel(colors) / 3);
    if (num_gs == 0)
        return;

    fused_3dgs2tr_color_optim_kernel<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        (float3*)_tv_f(colors),
        (float3*)_tv_f(vr_colors),
        (float3*)_tv_f(h_colors),
        _tv_f(opacities),
        (float3*)_tv_f(exp_avg_colors),
        (float3*)_tv_f(exp_avg_sq_colors),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step1),
        1.0f - powf(beta2, step2),
        eps,
        eps_tr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// 3DGS^2-TR Opacity - https://arxiv.org/abs/2602.00395
// ================

__global__ void fused_3dgs2tr_opacity_optim_kernel(
    const int num_gs,
    float* __restrict__ opacities,  // [N, 1], logit space
    const float* __restrict__ vr_opacities,  // [N, 1]
    const float* __restrict__ h_opacities,  // [N, 1]
    float* __restrict__ exp_avg_opacities,  // [N, 1]
    float* __restrict__ exp_avg_sq_opacities,  // [N, 1]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_gs) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;

        // Update momentum
        float vr_opacity = vr_opacities[idx];
        float h_opacity = h_opacities[idx];
        float exp_avg_opacity = exp_avg_opacities[idx];
        float exp_avg_sq_opacity = exp_avg_sq_opacities[idx];
        exp_avg_opacity = beta1 * exp_avg_opacity + (1.0f - beta1) * vr_opacity;
        exp_avg_sq_opacity = beta2 * exp_avg_sq_opacity + (1.0f - beta2) * h_opacity;
        exp_avg_opacities[idx] = exp_avg_opacity;
        exp_avg_sq_opacities[idx] = exp_avg_sq_opacity;

        // Compute delta
        float delta = -step_size * exp_avg_opacity / (exp_avg_sq_opacity / bias_correction2 + eps);

        // Compute trust region
        float opac = opacities[idx];
        opac = sigmoid(opac);
        float clip = sqrtf(4.0f * eps_tr * opac) / fmaxf(opac * (1.0f - opac), 1e-12f);

        // clip and update
        delta = fminf(fmaxf(delta, -clip), clip);
        delta = isfinite(delta) ? delta : 0.0f;
        opacities[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_3dgs2tr_opacity_optim(
    TorchTensorView opacities,
    TorchTensorView vr_opacities,
    TorchTensorView h_opacities,
    TorchTensorView exp_avg_opacities,
    TorchTensorView exp_avg_sq_opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    const int num_gs = (int)_tv_numel(opacities);
    if (num_gs == 0)
        return;

    fused_3dgs2tr_opacity_optim_kernel<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        _tv_f(opacities),
        _tv_f(vr_opacities),
        _tv_f(h_opacities),
        _tv_f(exp_avg_opacities),
        _tv_f(exp_avg_sq_opacities),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step1),
        1.0f - powf(beta2, step2),
        eps,
        eps_tr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// 3DGS^2-TR Quaternion - https://arxiv.org/abs/2602.00395
// ================

__device__ __forceinline__ float3x3 operator*(float3x3 A, float k) {
    return {A[0]*k, A[1]*k, A[2]*k};
}

__device__ float4 compute_delta_q_bound(
    float4 q, 
    float3 S_diag, 
    float alpha, 
    float epsilon
) {
    float w = q.x, x = q.y, y = q.z, z = q.w;
    float r2 = w*w + x*x + y*y + z*z;
    float r = sqrtf(r2);
    float r4 = r2 * r2;
    float r6 = r4 * r2;

    float3x3 R = SlangProjectionUtils::normalized_quat_to_rotmat(q / r);
    float3x3 R_T = float3x3{
        float3{R[0].x, R[1].x, R[2].x},
        float3{R[0].y, R[1].y, R[2].y},
        float3{R[0].z, R[1].z, R[2].z},
    };
    float3x3 R_tilde = R * r2;
    
    float3 S_inv_diag = make_float3(1.0f/S_diag.x, 1.0f/S_diag.y, 1.0f/S_diag.z);

    float4 bounds;

    #pragma unroll
    for (int i = 0; i < 4; ++i) {
        float qc = (i == 0) ? w : (i == 1) ? x : (i == 2) ? y : z;

        float3x3 dR_tilde_dc;
        float3x3 d2R_tilde_dc2;

        if (i == 0) { // w
            dR_tilde_dc = float3x3{
                float3{2*w, -2*z, 2*y},
                float3{2*z, 2*w, -2*x},
                float3{-2*y, 2*x, 2*w} };
            d2R_tilde_dc2 = float3x3{float3{2, 0, 0}, float3{0, 2, 0}, float3{0, 0, 2}};
        } else if (i == 1) { // x
            dR_tilde_dc = float3x3{
                float3{2*x, 2*y, 2*z},
                float3{2*y, -2*x, -2*w},
                float3{2*z, 2*w, -2*x} };
            d2R_tilde_dc2 = float3x3{float3{2, 0, 0}, float3{0, -2, 0}, float3{0, 0, -2}};
        } else if (i == 2) { // y
            dR_tilde_dc = float3x3{
                float3{-2*y, 2*x, 2*w},
                float3{2*x, 2*y, 2*z},
                float3{-2*w, 2*z, -2*y} };
            d2R_tilde_dc2 = float3x3{float3{-2, 0, 0}, float3{0, 2, 0}, float3{0, 0, -2}};
        } else { // z
            dR_tilde_dc = float3x3{
                float3{-2*z, -2*w, 2*x},
                float3{2*w, -2*z, 2*y},
                float3{2*x, 2*y, 2*z} };
            d2R_tilde_dc2 = float3x3{float3{-2, 0, 0}, float3{0, -2, 0}, float3{0, 0, 2}};
        }

        // TODO: make sure following LLM generated code is correct;
        // (already manually checked above)
        float3x3 dE = R_T * (dR_tilde_dc * (1.0f/r2) - R_tilde * (2.0f * qc / r4));

        float3x3 d2E = R_T * (dR_tilde_dc * (-4.0f * qc / r4) + 
                             d2R_tilde_dc2 * (1.0f / r2) + 
                             R_tilde * (8.0f * qc * qc / r6 - 2.0f / r4));

        float frob_sq = 0;
        frob_sq += powf(S_diag.x * dE[0].x * S_inv_diag.x, 2);
        frob_sq += powf(S_diag.x * dE[0].y * S_inv_diag.y, 2);
        frob_sq += powf(S_diag.x * dE[0].z * S_inv_diag.z, 2);
        frob_sq += powf(S_diag.y * dE[1].x * S_inv_diag.x, 2);
        frob_sq += powf(S_diag.y * dE[1].y * S_inv_diag.y, 2);
        frob_sq += powf(S_diag.y * dE[1].z * S_inv_diag.z, 2);
        frob_sq += powf(S_diag.z * dE[2].x * S_inv_diag.x, 2);
        frob_sq += powf(S_diag.z * dE[2].y * S_inv_diag.y, 2);
        frob_sq += powf(S_diag.z * dE[2].z * S_inv_diag.z, 2);

        float trace_d2E = d2E[0].x + d2E[1].y + d2E[2].z;
        float beta_c = 2.0f * frob_sq + 2.0f * trace_d2E;

        float log_term = logf(1.0f - (epsilon / fmaxf(alpha, 1e-12f)));
        float val = sqrtf(fmaxf(0.0f, -8.0f / fmaxf(beta_c, 1e-12f) * log_term));

        if (i == 0) bounds.x = val;  // w
        else if (i == 1) bounds.y = val;  // x
        else if (i == 2) bounds.z = val;  // y
        else bounds.w = val;  // z
    }

    return bounds;
}

__global__ void fused_3dgs2tr_quat_optim_kernel(
    const int num_gs,
    float4* __restrict__ quats,  // [N, 4], unnormalized
    const float4* __restrict__ vr_quats,  // [N, 4]
    const float4* __restrict__ h_quats,  // [N, 4]
    const float3* __restrict__ scales,  // [N, 3], log space
    const float* __restrict__ opacities,  // [N, 1], logit space
    float4* __restrict__ exp_avg_quats,  // [N, 4]
    float4* __restrict__ exp_avg_sq_quats,  // [N, 4]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_gs) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;

        // Update momentum
        float4 vr_quat = vr_quats[idx];
        float4 h_quat = h_quats[idx];
        float4 exp_avg_quat = exp_avg_quats[idx];
        float4 exp_avg_sq_quat = exp_avg_sq_quats[idx];
        exp_avg_quat = beta1 * exp_avg_quat + (1.0f - beta1) * vr_quat;
        exp_avg_sq_quat = beta2 * exp_avg_sq_quat + (1.0f - beta2) * h_quat;
        exp_avg_quats[idx] = exp_avg_quat;
        exp_avg_sq_quats[idx] = exp_avg_sq_quat;

        // Compute delta
        float4 delta = -step_size * exp_avg_quat / (exp_avg_sq_quat / bias_correction2 + eps);

        // Compute trust region
        float4 quat = quats[idx];
        float3 scale = scales[idx];
        float opac = opacities[idx];
        scale = {expf(scale.x), expf(scale.y), expf(scale.z)};
        opac = sigmoid(opac);
        float4 clip = compute_delta_q_bound(quat, scale, opac, eps_tr);

        // clip and update
        delta.x = fminf(fmaxf(delta.x, -clip.x), clip.x);
        delta.y = fminf(fmaxf(delta.y, -clip.y), clip.y);
        delta.z = fminf(fmaxf(delta.z, -clip.z), clip.z);
        delta.w = fminf(fmaxf(delta.w, -clip.w), clip.w);
        delta.x = isfinite(delta.x) ? delta.x : 0.0f;
        delta.y = isfinite(delta.y) ? delta.y : 0.0f;
        delta.z = isfinite(delta.z) ? delta.z : 0.0f;
        delta.w = isfinite(delta.w) ? delta.w : 0.0f;
        // quats[idx] += delta;
        quats[idx] = normalize(quats[idx] + delta);
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_3dgs2tr_quat_optim(
    TorchTensorView quats,
    TorchTensorView vr_quats,
    TorchTensorView h_quats,
    TorchTensorView scales,
    TorchTensorView opacities,
    TorchTensorView exp_avg_quats,
    TorchTensorView exp_avg_sq_quats,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    const int num_gs = (int)(_tv_numel(quats) / 4);
    if (num_gs == 0)
        return;

    fused_3dgs2tr_quat_optim_kernel<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        (float4*)_tv_f(quats),
        (float4*)_tv_f(vr_quats),
        (float4*)_tv_f(h_quats),
        (float3*)_tv_f(scales),
        _tv_f(opacities),
        (float4*)_tv_f(exp_avg_quats),
        (float4*)_tv_f(exp_avg_sq_quats),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step1),
        1.0f - powf(beta2, step2),
        eps,
        eps_tr
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Adam for linear RGB
// ================

__global__ void fused_adam_linear_rgb_optim_kernel(
    const int numel,
    float* __restrict__ rgbs,  // [N]
    const float* __restrict__ grad,  // [N]
    float* __restrict__ exp_avg,  // [N]
    float* __restrict__ exp_avg_sq,  // [N]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numel) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;
        
        // Load values
        float x = rgbs[idx];
        float v = grad[idx];
        float m_val = exp_avg[idx];
        float v_val = exp_avg_sq[idx];

        // Convert gradient to linear color space
        v /= SlangPixelWise::linear_rgb_to_srgb_grad(kSh0 * x + 0.5f);
        
        // Update biased first moment estimate: m_t = beta1 * m_{t-1} + (1 - beta1) * g_t
        m_val = beta1 * m_val + (1.0f - beta1) * v;
        
        // Update biased second raw moment estimate: v_t = beta2 * v_{t-1} + (1 - beta2) * g_t^2
        v_val = beta2 * v_val + (1.0f - beta2) * v * v;
        
        // Compute update: theta_t = theta_{t-1} - step_size * m_t / (sqrt(v_t / bias_correction2) + eps)
        float denom = sqrtf(v_val / bias_correction2) + eps;
        x = x - step_size * (m_val / denom);
        
        // Write back
        rgbs[idx] = x;
        exp_avg[idx] = m_val;
        exp_avg_sq[idx] = v_val;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_linear_rgb_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
) {
    const int numel = (int)_tv_numel(param);
    if (numel == 0)
        return;

    fused_adam_linear_rgb_optim_kernel<<<_LAUNCH_ARGS_1D(numel, 256)>>>(
        numel,
        _tv_f(param),
        _tv_f(grad),
        _tv_f(exp_avg),
        _tv_f(exp_avg_sq),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Trust Region Adam for RGB
// ================

template<bool is_linear, bool zero_grad>
__global__ void fused_adamtr_rgb_optim_kernel(
    const int num_gs,
    float3* __restrict__ rgbs,  // [N, 3]
    float3* __restrict__ grad,  // [N, 3]
    float3* __restrict__ exp_avg,  // [N, 3]
    float3* __restrict__ exp_avg_sq,  // [N, 3]
    float* __restrict__ opacities,  // [N]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr,
    const float grad_scale,
    const int step
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < num_gs) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;

        // Load values
        float3 x = rgbs[idx];
        float3 v = grad[idx];
        if constexpr (zero_grad) grad[idx] = make_float3(0.0f);
        v = grad_scale * v;
        float3 m_val = exp_avg[idx];
        float3 v_val = exp_avg_sq[idx];

        // Convert gradient to linear color space
        if (is_linear) {
            v.x /= SlangPixelWise::linear_rgb_to_srgb_grad(kSh0 * x.x + 0.5f);
            v.y /= SlangPixelWise::linear_rgb_to_srgb_grad(kSh0 * x.y + 0.5f);
            v.z /= SlangPixelWise::linear_rgb_to_srgb_grad(kSh0 * x.z + 0.5f);
        }
        
        // Update momentum
        m_val = beta1 * m_val + (1.0f - beta1) * v;
        v_val = beta2 * v_val + (1.0f - beta2) * v * v;

        // Compute delta
        float3 denom = sqrtf(v_val / bias_correction2) + eps;
        float3 delta = -step_size * (m_val / denom);

        // Compute trust region
        float opac = opacities[idx];
        opac = sigmoid(opac);
        // float3 rgb = fmaxf(kSh0 * x + 0.5f, (0.5f/255.0f)*(0.5f/255.0f));
        // float3 rgb = fmaxf(kSh0 * x + 0.5f, 0.5f/255.0f);
        float3 rgb = fmaxf(kSh0 * x + 0.5f, (1.0f/255.0f)*(1.0f/255.0f));
        float3 clip = kSh0 * sqrtf(4.0f * eps_tr * rgb / fmaxf(opac, 1e-12f));

        // clip and update
        delta.x = fminf(fmaxf(delta.x, -clip.x), clip.x);
        delta.y = fminf(fmaxf(delta.y, -clip.y), clip.y);
        delta.z = fminf(fmaxf(delta.z, -clip.z), clip.z);
        delta.x = isfinite(delta.x) ? delta.x : 0.0f;
        delta.y = isfinite(delta.y) ? delta.y : 0.0f;
        delta.z = isfinite(delta.z) ? delta.z : 0.0f;
        rgbs[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_adamtr_linear_rgb_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step,
    float grad_scale, bool zero_grad
) {
    const int num_gs = (int)(_tv_numel(param) / 3);
    if (num_gs == 0)
        return;

    auto kfn = zero_grad ? fused_adamtr_rgb_optim_kernel<true,  true>
                         : fused_adamtr_rgb_optim_kernel<true,  false>;
    kfn<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        (float3*)_tv_f(param),
        (float3*)_tv_f(grad),
        (float3*)_tv_f(exp_avg),
        (float3*)_tv_f(exp_avg_sq),
        _tv_f(opacities),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps,
        eps_tr,
        grad_scale,
        step
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void fused_adamtr_rgb_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step,
    float grad_scale, bool zero_grad
) {
    const int num_gs = (int)(_tv_numel(param) / 3);
    if (num_gs == 0)
        return;

    auto kfn = zero_grad ? fused_adamtr_rgb_optim_kernel<false, true>
                         : fused_adamtr_rgb_optim_kernel<false, false>;
    kfn<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        (float3*)_tv_f(param),
        (float3*)_tv_f(grad),
        (float3*)_tv_f(exp_avg),
        (float3*)_tv_f(exp_avg_sq),
        _tv_f(opacities),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps,
        eps_tr,
        grad_scale,
        step
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Trust Region Adam for linear RGB, SH coefficients
// ================

template<bool is_linear, bool zero_grad>
__global__ void fused_adamtr_rgb_sh_optim_kernel(
    const int num_params,
    const int num_sh,
    float* __restrict__ param,  // [N, K, 3]
    float* __restrict__ grad,  // [N, K, 3]
    float* __restrict__ exp_avg,  // [N, K, 3]
    float* __restrict__ exp_avg_sq,  // [N, K, 3]
    float* __restrict__ rgbs,  // [N, 3]
    float* __restrict__ opacities,  // [N]
    const float lr,
    const float beta1,
    const float beta2,
    const float bias_correction1,
    const float bias_correction2,
    const float eps,
    const float eps_tr,
    const float grad_scale,
    const int step
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < num_params) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float step_size = lr / bias_correction1;

        // Load values
        // float x = param[idx];
        float v = grad[idx];
        if constexpr (zero_grad) grad[idx] = 0.0f;
        v = grad_scale * v;
        float m_val = exp_avg[idx];
        float v_val = exp_avg_sq[idx];

        // Convert gradient to linear color space
        float c = rgbs[(idx / (3*num_sh)) * 3 + (idx % 3)] * kSh0 + 0.5f;
        if (is_linear)
            v /= SlangPixelWise::linear_rgb_to_srgb_grad(c);
        
        // Update momentum
        m_val = beta1 * m_val + (1.0f - beta1) * v;
        v_val = beta2 * v_val + (1.0f - beta2) * v * v;

        // Compute delta
        float denom = sqrtf(v_val / bias_correction2) + eps;
        float delta = -step_size * (m_val / denom);

        // Compute trust region
        float opac = opacities[idx / (3*num_sh)];
        opac = sigmoid(opac);
        // c = fmaxf(c, (0.5f/255.0f)*(0.5f/255.0f));
        // c = fmaxf(c, 0.5f/255.0f);
        c = fmaxf(c, (1.0f/255.0f)*(1.0f/255.0f));
        float clip = kSh0 * sqrtf(4.0f * eps_tr * c / fmaxf(opac, 1e-12f));

        // clip and update
        delta = fminf(fmaxf(delta, -clip), clip);
        param[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_adamtr_linear_rgb_sh_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView colors,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step,
    float grad_scale, bool zero_grad
) {
    int64_t colors_numel = _tv_numel(colors);
    const int num_gs = (int)(colors_numel / 3);
    if (num_gs == 0)
        return;
    const int num_sh = (int)(_tv_numel(param) / colors_numel);
    if (num_sh == 0)
        return;

    auto kfn = zero_grad ? fused_adamtr_rgb_sh_optim_kernel<true,  true>
                         : fused_adamtr_rgb_sh_optim_kernel<true,  false>;
    kfn<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        num_sh,
        _tv_f(param),
        _tv_f(grad),
        _tv_f(exp_avg),
        _tv_f(exp_avg_sq),
        _tv_f(colors),
        _tv_f(opacities),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps,
        eps_tr,
        grad_scale,
        step
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void fused_adamtr_rgb_sh_optim(
    TorchTensorView param,
    TorchTensorView grad,
    TorchTensorView exp_avg,
    TorchTensorView exp_avg_sq,
    TorchTensorView colors,
    TorchTensorView opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step,
    float grad_scale, bool zero_grad
) {
    int64_t colors_numel = _tv_numel(colors);
    const int num_gs = (int)(colors_numel / 3);
    if (num_gs == 0)
        return;
    const int num_sh = (int)(_tv_numel(param) / colors_numel);
    if (num_sh == 0)
        return;

    auto kfn = zero_grad ? fused_adamtr_rgb_sh_optim_kernel<false, true>
                         : fused_adamtr_rgb_sh_optim_kernel<false, false>;
    kfn<<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
        num_gs,
        num_sh,
        _tv_f(param),
        _tv_f(grad),
        _tv_f(exp_avg),
        _tv_f(exp_avg_sq),
        _tv_f(colors),
        _tv_f(opacities),
        lr,
        beta1,
        beta2,
        1.0f - powf(beta1, step),
        1.0f - powf(beta2, step),
        eps,
        eps_tr,
        grad_scale,
        step
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



__global__ void increment_int32_kernel(int32_t* __restrict__ data, int64_t n) {
    int64_t idx = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) data[idx] += 1;
}

/*[AutoHeaderGeneratorExport]*/
void increment_int32_inplace(DeviceVector<int32_t> data, int64_t n) {
    if (n == 0 || data.data_ptr() == nullptr) return;
    increment_int32_kernel<<<(n + 255) / 256, 256>>>(data.data_ptr(), n);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// dst[i] += src[i] for i in [0, n). Used by the sub-batched training
// dispatcher to accumulate per-sub-batch raster-bwd accum_weight buffers
// into a persistent buffer that densify consumes at the end of the step.
__global__ void float_add_into_kernel(float* __restrict__ dst,
                                      const float* __restrict__ src,
                                      int64_t n) {
    int64_t idx = (int64_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) dst[idx] += src[idx];
}

/*[AutoHeaderGeneratorExport]*/
void float_add_into(DeviceVector<float> dst, DeviceVector<float> src, int64_t n) {
    if (n == 0 || dst.data_ptr() == nullptr || src.data_ptr() == nullptr) return;
    float_add_into_kernel<<<(n + 255) / 256, 256>>>(dst.data_ptr(), src.data_ptr(), n);
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
