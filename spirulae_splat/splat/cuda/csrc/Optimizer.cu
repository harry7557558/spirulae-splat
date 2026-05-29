#include "types.cuh"

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

#include "common.cuh"

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

    float* buffer = DevicePool::global().acquire<float>("semi_offloaded_adam_buf", 2 * chunk_size);
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

template<bool use_scale_agnostic_mean>
__global__ void fused_optim_3dgs_geometry_kernel(
    float3* __restrict__ means,
    const float3* __restrict__ v_means,
    float3* __restrict__ g1_means,
    float3* __restrict__ g2_means,
    float4* __restrict__ quats,
    const float4* __restrict__ v_quats,
    float4* __restrict__ g1_quats,
    float4* __restrict__ g2_quats,
    float3* __restrict__ scales,
    const float3* __restrict__ v_scales,
    float3* __restrict__ g1_scales,
    float3* __restrict__ g2_scales,
    float* __restrict__ opacities,
    const float* __restrict__ v_opacities,
    float* __restrict__ g1_opacities,
    float* __restrict__ g2_opacities,
    const float* __restrict__ radii,
    float lr_means,
    float lr_quats,
    float lr_scales,
    float lr_opacs,
    const float max_gauss_ratio,
    const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight,
    const float mcmc_scale_reg_weight,
    const float erank_reg_weight,
    const float erank_reg_weight_s3,
    const float quat_norm_reg_weight,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const int64_t numel
) {
    static constexpr float eps = 1e-15f;
    static constexpr float beta1 = 0.9f;
    static constexpr float beta2 = 0.999f;

    const int64_t idx = (int64_t)blockIdx.x * (int64_t)blockDim.x + (int64_t)threadIdx.x;
    if (idx >= numel) return;

    float step = (float)(steps != nullptr ? steps[idx] : scalar_step);
    float inv_bias_correction1 = 1.0f / (1.0f - powf(beta1, step));
    float inv_bias_correction2 = 1.0f / (1.0f - powf(beta2, step));
    lr_means *= inv_bias_correction1;
    lr_quats *= inv_bias_correction1;
    lr_scales *= inv_bias_correction1;
    lr_opacs *= inv_bias_correction1;

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
    v_scale += v_scales[idx];
    v_quat += v_quats[idx];
    v_opac += v_opacities[idx];

    // update scales
    float3 g1_scale = beta1 * g1_scales[idx] + (1.f - beta1) * v_scale;
    float3 g2_scale = beta2 * g2_scales[idx] + (1.f - beta2) * v_scale*v_scale;
    float3 updated_scale = scale - lr_scales * g1_scale / (sqrtf(g2_scale * inv_bias_correction2) + eps);
    scales[idx] = updated_scale;
    g1_scales[idx] = g1_scale;
    g2_scales[idx] = g2_scale;

    // update quats (Riemannian)
    v_quat -= dot(quat, v_quat) * quat;
    float4 g1_quat = beta1 * g1_quats[idx] + (1.f - beta1) * v_quat;
    float4 g2_quat = beta2 * g2_quats[idx] + (1.f - beta2) * v_quat*v_quat;
    quats[idx] = normalize(quat - lr_quats * g1_quat / (sqrtf(g2_quat * inv_bias_correction2) + eps));
    g1_quats[idx] = g1_quat;
    g2_quats[idx] = g2_quat;

    // update opacs
    float g1_opac = beta1 * g1_opacities[idx] + (1.f - beta1) * v_opac;
    float g2_opac = beta2 * g2_opacities[idx] + (1.f - beta2) * v_opac*v_opac;
    float updated_opac = opac - lr_opacs * g1_opac / (sqrtf(g2_opac * inv_bias_correction2) + eps);
    opacities[idx] = updated_opac;
    g1_opacities[idx] = g1_opac;
    g2_opacities[idx] = g2_opac;

    // update means (scale agnostic)
    float3 mean = means[idx];
    float3 v_mean = v_means[idx];
    float3 g1_mean = g1_means[idx];
    float3 g2_mean = g2_means[idx];

    float opac_post_sigmoid = sigmoid(opac);

    if constexpr (use_scale_agnostic_mean) {
        float3 v_mean_scaled_num = SlangProjectionUtils::apply_covar_to_vec(
            quat,
            {expf(0.5f*scale.x), expf(0.5f*scale.y), expf(0.5f*scale.z)},  // unit: L^0.5
            v_mean  // unit: L^-1
        ) * sqrtf(2.0f * __logf(fmaxf(255.0f * opac_post_sigmoid, 1.00001f)));  // unit: dimensionless
    #if 1
        float v_mean_scaled_den = radii[idx] * 0.6f;  // unit: dimensionless
    #else
        float3 scale_post_exp = float3{__expf(scale.x), __expf(scale.y), __expf(scale.z)};
        v_mean_scaled_num = v_mean_scaled_num * (
            fmaxf(scale_post_exp.x, fmaxf(scale_post_exp.y, scale_post_exp.z)) *
                powf(powf(scale_post_exp.x*scale_post_exp.y, 1.6f) +
                    powf(scale_post_exp.x*scale_post_exp.z, 1.6f) +
                    powf(scale_post_exp.y*scale_post_exp.z, 1.6f), -0.5f / 1.6f)  // proportional to longest axis to visible radius ratio
        );
        float v_mean_scaled_den = radii[idx] * 0.9f;
    #endif
    #if 0
        v_mean_scaled_num = v_mean_scaled_num / (v_mean_scaled_den + (eps * (float)numel));
        v_mean_scaled_den = 1.0f;
    #endif
        g1_mean = beta1 * g1_mean + (1.f - beta1) * v_mean_scaled_num;  // unit: dimensionless
        g2_mean = beta2 * g2_mean + (1.f - beta2) * v_mean*v_mean * v_mean_scaled_den*v_mean_scaled_den;  // unit: L^-2
    } else {
        g1_mean = beta1 * g1_mean + (1.f - beta1) * v_mean;  // unit: L^-1
        g2_mean = beta2 * g2_mean + (1.f - beta2) * v_mean*v_mean;  // unit: L^-2
    }

    means[idx] = mean - lr_means * g1_mean / (sqrtf(g2_mean * inv_bias_correction2) + eps); // unit: L or dimensionless for dimensionless lr_means
    g1_means[idx] = g1_mean;
    g2_means[idx] = g2_mean;
}


/*[AutoHeaderGeneratorExport]*/
void fused_optim_3dgs_geometry(
    int64_t num_splats,
    DeviceVector<float3> means, DeviceVector<float3> v_means, DeviceVector<float3> g1_means, DeviceVector<float3> g2_means,
    DeviceVector<float4> quats, DeviceVector<float4> v_quats, DeviceVector<float4> g1_quats, DeviceVector<float4> g2_quats,
    DeviceVector<float3> scales, DeviceVector<float3> v_scales, DeviceVector<float3> g1_scales, DeviceVector<float3> g2_scales,
    DeviceVector<float> opacities, DeviceVector<float> v_opacities, DeviceVector<float> g1_opacities, DeviceVector<float> g2_opacities,
    DeviceVector<float> radii,
    const float lr_means, const float lr_quats, const float lr_scales, const float lr_opacs,
    const float max_gauss_ratio, const float scale_regularization_weight,
    const float mcmc_opacity_reg_weight, const float mcmc_scale_reg_weight,
    const float erank_reg_weight, const float erank_reg_weight_s3, const float quat_norm_reg_weight,
    bool use_scale_agnostic_mean,
    int32_t step, DeviceVector<int32_t> per_splat_steps
) {
    if (num_splats == 0)
        return;

    (use_scale_agnostic_mean ?
        fused_optim_3dgs_geometry_kernel<true> :
        fused_optim_3dgs_geometry_kernel<false>
    )<<<_LAUNCH_ARGS_1D(num_splats, 256)>>>(
        means.data_ptr(), v_means.data_ptr(), g1_means.data_ptr(), g2_means.data_ptr(),
        quats.data_ptr(), v_quats.data_ptr(), g1_quats.data_ptr(), g2_quats.data_ptr(),
        scales.data_ptr(), v_scales.data_ptr(), g1_scales.data_ptr(), g2_scales.data_ptr(),
        opacities.data_ptr(), v_opacities.data_ptr(), g1_opacities.data_ptr(), g2_opacities.data_ptr(),
        (const float*)radii.data_ptr(),
        lr_means, lr_quats, lr_scales, lr_opacs,
        max_gauss_ratio,
        scale_regularization_weight / (float)num_splats,
        mcmc_opacity_reg_weight / (float)num_splats,
        mcmc_scale_reg_weight / (float)num_splats,
        erank_reg_weight / (float)num_splats,
        erank_reg_weight_s3 / (float)num_splats,
        quat_norm_reg_weight / (float)num_splats,
        step, per_splat_steps.data_ptr(), num_splats
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}



// ================
// Fused appearance optimizer
// ================

__global__ void fused_adam_with_steps_kernel(
    float* __restrict__ param,
    const float* __restrict__ grad,
    float* __restrict__ exp_avg,
    float* __restrict__ exp_avg_sq,
    const float lr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const float decay,
    const float decay_offset,
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
    float l2_reg_offset
) {
    int64_t param_numel = param.numel();
    if (param_numel == 0 || num_splats == 0)
        return;
    int stride = (int)(param_numel / num_splats);

    fused_adam_with_steps_kernel<<<_LAUNCH_ARGS_1D(num_splats*stride, 256)>>>(
        param.data_ptr(),
        grad.data_ptr(),
        exp_avg.data_ptr(),
        exp_avg_sq.data_ptr(),
        lr,
        step, per_splat_steps.data_ptr(),
        2.0f*l2_reg/(float)(num_splats*stride),
        l2_reg_offset,
        num_splats*stride,
        stride
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


template<int BLOCK_SIZE>
__global__ void fused_adam_with_steps_8bit_kernel(
    float* __restrict__ param,
    const float* __restrict__ grad,
    uint8_t* __restrict__ exp_avg,
    uint8_t* __restrict__ exp_avg_sq,
    float4* __restrict__ quant_bounds,  // g1 min, g1 max, g2 min, g2 max
    const float lr,
    const int32_t scalar_step,
    const int32_t* __restrict__ steps,
    const float decay,
    const float decay_offset,
    const int64_t numel,
    const int stride
) {
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
    v += decay * (fmaxf(x - decay_offset, 0.0f) + fminf(x + decay_offset, 0.0f));

    // Joint (u, sqrt(g2)) quantization:
    //   exp_avg     stores u = g1 / (sqrt(g2) + eps)
    //   exp_avg_sq  stores sqrt(g2)
    //   quant_bounds.{xy} = (u_min, u_max)
    //   quant_bounds.{zw} = (sqrt_g2_min, sqrt_g2_max)
    // This keeps the on-device storage budget at 2 bytes/cell (+ shared
    // float4 per 256-cell block) but preserves the Adam update direction
    // g1/sqrt(g2) much better than per-component encoding of g1 and sqrt(g2).
    // Endpoint-exact decode: q=0 -> lo, q=255 -> hi exactly. Eliminates
    // the +range/512 phantom-value bias that midpoint decoding (q+0.5)/256
    // induces for non-uniform-within-block distributions.
    float u_norm  = inside ? (float)exp_avg[idx]    * (1.0f / 255.0f) : 0.0f;
    float s_norm  = inside ? (float)exp_avg_sq[idx] * (1.0f / 255.0f) : 0.0f;
    float4 mm = quant_bounds[blockIdx.x];
    float u_val   = mm.x + (mm.y - mm.x) * u_norm;
    float sqrt_g2 = mm.z + (mm.w - mm.z) * s_norm;
    float g2 = sqrt_g2 * sqrt_g2;
    float g1 = u_val * (sqrt_g2 + eps);

    g1 = beta1 * g1 + (1.0f - beta1) * v;
    g2 = beta2 * g2 + (1.0f - beta2) * v*v;

    x -= lr * inv_bias_correction1 * g1 / (sqrtf(g2 * inv_bias_correction2) + eps);
    param[idx] = x;

    // Re-encode the new Adam state in the (u, sqrt(g2)) basis.
    float sqrt_g2_new = sqrtf(g2);
    float u_new       = g1 / (sqrt_g2_new + eps);

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

    float u_range = fmaxf(mm.y - mm.x, eps);
    float s_range = fmaxf(mm.w - mm.z, eps);
    // Endpoint-exact encode: q = round(255 * (x - lo) / range), clamped.
    exp_avg[idx]    = (uint8_t)fminf(fmaxf(roundf(255.0f * (u_new       - mm.x) / u_range), 0.0f), 255.0f);
    exp_avg_sq[idx] = (uint8_t)fminf(fmaxf(roundf(255.0f * (sqrt_g2_new - mm.z) / s_range), 0.0f), 255.0f);

    if (threadIdx.x == 0)
        quant_bounds[blockIdx.x] = mm;
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_step_8bit(
    int64_t num_splats,
    DeviceTensorFloatND param,
    DeviceTensorFloatND grad,
    uint8_t* exp_avg,
    uint8_t* exp_avg_sq,
    float4* quant_bounds,
    float lr,
    int32_t step, DeviceVector<int32_t> per_splat_steps,
    float l2_reg,
    float l2_reg_offset
) {
    int64_t param_numel = param.numel();
    if (param_numel == 0 || num_splats == 0)
        return;
    int stride = (int)(param_numel / num_splats);
    constexpr int BLOCK_SIZE = 256;

    fused_adam_with_steps_8bit_kernel<BLOCK_SIZE><<<_LAUNCH_ARGS_1D(num_splats*stride, BLOCK_SIZE)>>>(
        param.data_ptr(),
        grad.data_ptr(),
        exp_avg,
        exp_avg_sq,
        quant_bounds,
        lr,
        step, per_splat_steps.data_ptr(),
        2.0f*l2_reg/(float)(num_splats*stride),
        l2_reg_offset,
        num_splats*stride,
        stride
    );
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

template<bool is_linear>
__global__ void fused_adamtr_rgb_optim_kernel(
    const int num_gs,
    float3* __restrict__ rgbs,  // [N, 3]
    const float3* __restrict__ grad,  // [N, 3]
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
    int step
) {
    const int num_gs = (int)(_tv_numel(param) / 3);
    if (num_gs == 0)
        return;

    fused_adamtr_rgb_optim_kernel<true><<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
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
    int step
) {
    const int num_gs = (int)(_tv_numel(param) / 3);
    if (num_gs == 0)
        return;

    fused_adamtr_rgb_optim_kernel<false><<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
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
        step
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


// ================
// Trust Region Adam for linear RGB, SH coefficients
// ================

template<bool is_linear>
__global__ void fused_adamtr_rgb_sh_optim_kernel(
    const int num_params,
    const int num_sh,
    float* __restrict__ param,  // [N, K, 3]
    const float* __restrict__ grad,  // [N, K, 3]
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
    int step
) {
    int64_t colors_numel = _tv_numel(colors);
    const int num_gs = (int)(colors_numel / 3);
    if (num_gs == 0)
        return;
    const int num_sh = (int)(_tv_numel(param) / colors_numel);
    if (num_sh == 0)
        return;

    fused_adamtr_rgb_sh_optim_kernel<true><<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
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
    int step
) {
    int64_t colors_numel = _tv_numel(colors);
    const int num_gs = (int)(colors_numel / 3);
    if (num_gs == 0)
        return;
    const int num_sh = (int)(_tv_numel(param) / colors_numel);
    if (num_sh == 0)
        return;

    fused_adamtr_rgb_sh_optim_kernel<false><<<_LAUNCH_ARGS_1D(num_gs, 256)>>>(
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
