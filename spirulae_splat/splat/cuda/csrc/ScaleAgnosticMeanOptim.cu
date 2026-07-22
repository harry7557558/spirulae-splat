// ScaleAgnosticMeanOptim.cu -- scale-agnostic mean Adam.
//
// Part of the Optimizer family -- see OptimizerCommon.cuh.

#include "OptimizerCommon.cuh"

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



