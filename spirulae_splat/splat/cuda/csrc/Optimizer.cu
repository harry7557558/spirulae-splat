#include "types.cuh"

#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
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



/* 3DGS^2-TR Mean - https://arxiv.org/abs/2602.00395 */

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
        opac = 1.0f / (1.0f + __expf(-opac));
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
    at::Tensor means,
    at::Tensor vr_means,
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
    CHECK_INPUT(vr_means);
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
        (float3*)vr_means.data_ptr<float>(),
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


/* 3DGS^2-TR Scale - https://arxiv.org/abs/2602.00395 */
/* This can also be used for colors (equivalent to with c=0.5, otherwise adjust eps_tr accordingly) */

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
        opac = 1.0f / (1.0f + __expf(-opac));
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
    at::Tensor scales,
    at::Tensor vr_scales,
    at::Tensor h_scales,
    at::Tensor opacities,
    at::Tensor exp_avg_scales,
    at::Tensor exp_avg_sq_scales,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    DEVICE_GUARD(scales);
    CHECK_INPUT(scales);
    CHECK_INPUT(vr_scales);
    CHECK_INPUT(h_scales);
    CHECK_INPUT(opacities);
    CHECK_INPUT(exp_avg_scales);
    CHECK_INPUT(exp_avg_sq_scales);

    const int num_gs = scales.numel() / 3;
    const int threads = 256;
    const int blocks = (num_gs + threads - 1) / threads;
    
    fused_3dgs2tr_scale_optim_kernel<<<blocks, threads>>>(
        num_gs,
        (float3*)scales.data_ptr<float>(),
        (float3*)vr_scales.data_ptr<float>(),
        (float3*)h_scales.data_ptr<float>(),
        opacities.data_ptr<float>(),
        (float3*)exp_avg_scales.data_ptr<float>(),
        (float3*)exp_avg_sq_scales.data_ptr<float>(),
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


/* 3DGS^2-TR Opacity - https://arxiv.org/abs/2602.00395 */

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
        opac = 1.0f / (1.0f + __expf(-opac));
        float clip = sqrtf(4.0f * eps_tr * opac) / fmaxf(opac * (1.0f - opac), 1e-12f);

        // clip and update
        delta = fminf(fmaxf(delta, -clip), clip);
        delta = isfinite(delta) ? delta : 0.0f;
        opacities[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_3dgs2tr_opacity_optim(
    at::Tensor opacities,
    at::Tensor vr_opacities,
    at::Tensor h_opacities,
    at::Tensor exp_avg_opacities,
    at::Tensor exp_avg_sq_opacities,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    DEVICE_GUARD(opacities);
    CHECK_INPUT(opacities);
    CHECK_INPUT(vr_opacities);
    CHECK_INPUT(h_opacities);
    CHECK_INPUT(exp_avg_opacities);
    CHECK_INPUT(exp_avg_sq_opacities);

    const int num_gs = opacities.numel();
    const int threads = 256;
    const int blocks = (num_gs + threads - 1) / threads;
    
    fused_3dgs2tr_opacity_optim_kernel<<<blocks, threads>>>(
        num_gs,
        opacities.data_ptr<float>(),
        vr_opacities.data_ptr<float>(),
        h_opacities.data_ptr<float>(),
        exp_avg_opacities.data_ptr<float>(),
        exp_avg_sq_opacities.data_ptr<float>(),
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


/* 3DGS^2-TR Quaternion - https://arxiv.org/abs/2602.00395 */

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
        opac = 1.0f / (1.0f + __expf(-opac));
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
        quats[idx] += delta;
    }
}

/*[AutoHeaderGeneratorExport]*/
void fused_3dgs2tr_quat_optim(
    at::Tensor quats,
    at::Tensor vr_quats,
    at::Tensor h_quats,
    at::Tensor scales,
    at::Tensor opacities,
    at::Tensor exp_avg_quats,
    at::Tensor exp_avg_sq_quats,
    float lr,
    float beta1,
    float beta2,
    float eps,
    float eps_tr,
    int step1,
    int step2
) {
    DEVICE_GUARD(quats);
    CHECK_INPUT(quats);
    CHECK_INPUT(vr_quats);
    CHECK_INPUT(h_quats);
    CHECK_INPUT(scales);
    CHECK_INPUT(opacities);
    CHECK_INPUT(exp_avg_quats);
    CHECK_INPUT(exp_avg_sq_quats);

    const int num_gs = quats.numel() / 4;
    const int threads = 256;
    const int blocks = (num_gs + threads - 1) / threads;
    
    fused_3dgs2tr_quat_optim_kernel<<<blocks, threads>>>(
        num_gs,
        (float4*)quats.data_ptr<float>(),
        (float4*)vr_quats.data_ptr<float>(),
        (float4*)h_quats.data_ptr<float>(),
        (float3*)scales.data_ptr<float>(),
        opacities.data_ptr<float>(),
        (float4*)exp_avg_quats.data_ptr<float>(),
        (float4*)exp_avg_sq_quats.data_ptr<float>(),
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

