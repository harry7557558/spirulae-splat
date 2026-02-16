#include "types.cuh"

#include "generated/slang.cuh"
namespace SlangAll {
#include "generated/set_namespace.cuh"
#include "generated/slang_all.cuh"
}

#include "common.cuh"


/* "Standard" Adam */

template <typename scalar_t>
__global__ void fused_adam_kernel(
    scalar_t* __restrict__ param,
    const scalar_t* __restrict__ grad,
    scalar_t* __restrict__ exp_avg,
    scalar_t* __restrict__ exp_avg_sq,
    const float lr,
    const float beta1,
    const float beta2,
    const float eps,
    const int step,
    const int numel
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numel) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float bias_correction1 = 1.0f - powf(beta1, step);
        const float bias_correction2 = 1.0f - powf(beta2, step);
        const float step_size = lr / bias_correction1;
        
        // Load values
        const float grad_val = static_cast<float>(grad[idx]);
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
    at::Tensor param,
    at::Tensor grad,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step
) {
    DEVICE_GUARD(param);
    CHECK_INPUT(param);
    CHECK_INPUT(grad);
    CHECK_INPUT(exp_avg);
    CHECK_INPUT(exp_avg_sq);

    const int numel = param.numel();
    const int threads = 256;
    const int blocks = (numel + threads - 1) / threads;
    
    fused_adam_kernel<float><<<blocks, threads>>>(
        param.data_ptr<float>(),
        grad.data_ptr<float>(),
        exp_avg.data_ptr<float>(),
        exp_avg_sq.data_ptr<float>(),
        lr,
        beta1,
        beta2,
        eps,
        step,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void fused_adam_multi(
    std::vector<at::Tensor> params,
    std::vector<at::Tensor> grads,
    std::vector<at::Tensor> exp_avgs,
    std::vector<at::Tensor> exp_avg_sqs,
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


/* Newton */

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
    const float eps,
    const int step1,
    const int step2,
    const int numel
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < numel) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float bias_correction1 = 1.0f - powf(beta1, step1);
        const float bias_correction2 = 1.0f - powf(beta2, step2);
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
    at::Tensor param,
    at::Tensor grad,
    at::Tensor hess_diag,
    at::Tensor exp_avg,
    at::Tensor exp_avg_sq,
    float lr,
    float beta1,
    float beta2,
    float eps,
    int step1,
    int step2
) {
    DEVICE_GUARD(param);
    CHECK_INPUT(param);
    CHECK_INPUT(grad);
    CHECK_INPUT(hess_diag);
    CHECK_INPUT(exp_avg);
    CHECK_INPUT(exp_avg_sq);

    const int numel = param.numel();
    const int threads = 256;
    const int blocks = (numel + threads - 1) / threads;
    
    fused_newton_kernel<float><<<blocks, threads>>>(
        param.data_ptr<float>(),
        grad.data_ptr<float>(),
        hess_diag.data_ptr<float>(),
        exp_avg.data_ptr<float>(),
        exp_avg_sq.data_ptr<float>(),
        lr,
        beta1,
        beta2,
        eps,
        step1,
        step2,
        numel
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

/*[AutoHeaderGeneratorExport]*/
void fused_newton_multi(
    std::vector<at::Tensor> params,
    std::vector<at::Tensor> grads,
    std::vector<at::Tensor> hess_diags,
    std::vector<at::Tensor> exp_avgs,
    std::vector<at::Tensor> exp_avg_sqs,
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



/* 3DGS Mean - https://arxiv.org/abs/2602.00395 */

__global__ void fused_3dgs2tr_mean_optim_kernel(
    const int num_gs,
    float3* __restrict__ means,  // [N, 3]
    const float3* __restrict__ v_means,  // [N, 3]
    const float3* __restrict__ h_means,  // [N, 3]
    const float3* __restrict__ scales,  // [N, 3], log space
    const float4* __restrict__ quats,  // [N, 4], unnormalized
    const float* __restrict__ opacities,  // [N, 1], logit space
    float3* __restrict__ exp_avg_means,  // [N, 3]
    float3* __restrict__ exp_avg_sq_means,  // [N, 3]
    const float lr,
    const float beta1,
    const float beta2,
    const float eps,
    const float eps_tr,
    const int step1,
    const int step2
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (idx < num_gs) {
        // Bias correction terms
        // TODO: proper bias correction after densification
        const float bias_correction1 = 1.0f - powf(beta1, step1);
        const float bias_correction2 = 1.0f - powf(beta2, step2);
        const float step_size = lr / bias_correction1;

        // Update momentum
        float3 v_mean = v_means[idx];
        float3 h_mean = h_means[idx];
        float3 exp_avg_mean = exp_avg_means[idx];
        float3 exp_avg_sq_mean = exp_avg_sq_means[idx];
        exp_avg_mean = beta1 * exp_avg_mean + (1.0f - beta1) * v_mean;
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
        opac = 1.0f / (1.0f + __expf(-opac));
        Matrix<float, 3, 3> covar;
        SlangAll::quat_scale_to_covar(quat, scale, &covar);
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
    at::Tensor means,
    at::Tensor v_means,
    at::Tensor h_means,
    at::Tensor scales,
    at::Tensor quats,
    at::Tensor opacities,
    at::Tensor exp_avg_means,
    at::Tensor exp_avg_sq_means,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    DEVICE_GUARD(means);
    CHECK_INPUT(means);
    CHECK_INPUT(v_means);
    CHECK_INPUT(h_means);
    CHECK_INPUT(scales);
    CHECK_INPUT(quats);
    CHECK_INPUT(opacities);
    CHECK_INPUT(exp_avg_means);
    CHECK_INPUT(exp_avg_sq_means);

    const int num_gs = means.numel() / 3;
    const int threads = 256;
    const int blocks = (num_gs + threads - 1) / threads;
    
    fused_3dgs2tr_mean_optim_kernel<<<blocks, threads>>>(
        num_gs,
        (float3*)means.data_ptr<float>(),
        (float3*)v_means.data_ptr<float>(),
        (float3*)h_means.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(),
        (float4*)quats.data_ptr<float>(),
        opacities.data_ptr<float>(),
        (float3*)exp_avg_means.data_ptr<float>(),
        (float3*)exp_avg_sq_means.data_ptr<float>(),
        lr,
        beta1,
        beta2,
        eps,
        eps_tr,
        step1,
        step2
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
