// NewtonOptim.cu -- Newton optimizer (single- and multi-tensor).
//
// Part of the Optimizer family -- see OptimizerCommon.cuh.

#include "OptimizerCommon.cuh"

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



