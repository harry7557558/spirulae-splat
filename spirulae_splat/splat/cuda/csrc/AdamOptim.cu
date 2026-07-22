// AdamOptim.cu -- Adam variants: standard, multi-tensor, stepped, 8-bit, and the Riemannian quaternion form.
//
// Part of the Optimizer family -- see OptimizerCommon.cuh.

#include "OptimizerCommon.cuh"

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


