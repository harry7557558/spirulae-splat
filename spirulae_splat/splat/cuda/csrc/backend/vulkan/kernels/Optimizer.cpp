// Vulkan implementation of the optimizer launch APIs the portable engine
// references (csrc/Optimizer.cuh): fused_adam_step (fp32 / quantized-state /
// doubly-quantized), fused_adagrad_step, float_add_into,
// increment_int32_inplace. Device work: slang/vulkan/optimizer.slang. The
// remaining Optimizer.cu entry points (adamtr color variants, the fused
// 3DGS geometry step) land with the densify/MCMC port.

#include <Optimizer.cuh>

#include "KernelCommon.h"

namespace {

// Mirrors FusedAdamParams in slang/vulkan/optimizer.slang.
struct FusedAdamParams {
    uint64_t param, grad, exp_avg, exp_avg_sq, steps;
    float lr, decay, decay_offset, grad_scale;
    int32_t scalar_step;
    uint32_t has_steps, numel, stride, zero_grad, wgs_per_row;
};
static_assert(sizeof(FusedAdamParams) == 5 * 8 + 10 * 4,
              "params layout must match the slang struct");

// Mirrors AdagradParams.
struct AdagradParams {
    uint64_t param, grad, accum;
    float lr;
    uint32_t numel, wgs_per_row;
};
static_assert(sizeof(AdagradParams) == 3 * 8 + 3 * 4 + 4 /*pad*/,
              "params layout must match the slang struct");

// Mirrors QAdamParams.
struct QAdamParams {
    uint64_t param, grad, packed, quant_bounds, steps;
    float lr, decay, decay_offset, grad_scale;
    int32_t scalar_step;
    uint32_t has_steps, numel, stride, zero_grad, num_blocks, wgs_per_row;
    uint32_t _pad0;
};
static_assert(sizeof(QAdamParams) == 5 * 8 + 12 * 4,
              "params layout must match the slang struct");

// Mirrors QQAdamParams.
struct QQAdamParams {
    uint64_t grad, grad_q_packed, grad_q_bounds, optim_packed, optim_bounds,
        value_packed, value_bounds, steps;
    float lr, decay, decay_offset, grad_scale;
    int32_t scalar_step;
    uint32_t has_steps, grad_quant, numel, stride, zero_grad, num_blocks,
        wgs_per_row;
};
static_assert(sizeof(QQAdamParams) == 8 * 8 + 12 * 4,
              "params layout must match the slang struct");

// Mirrors FloatAddParams.
struct FloatAddParams {
    uint64_t dst, src;
    uint32_t n, wgs_per_row;
};
static_assert(sizeof(FloatAddParams) == 2 * 8 + 2 * 4, "layout");

// Mirrors IncI32Params.
struct IncI32Params {
    uint64_t data;
    uint32_t n, wgs_per_row;
};
static_assert(sizeof(IncI32Params) == 8 + 2 * 4, "layout");

uint32_t checked_u32_numel(int64_t numel, const char* what) {
    if (numel < 0 || numel > (int64_t)UINT32_MAX)
        throw std::runtime_error(std::string(what) +
                                 ": numel exceeds the u32 cell-index range");
    return (uint32_t)numel;
}

}  // namespace

/* API definitions matching csrc/Optimizer.cuh (engine-referenced subset) */

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
    int64_t stride = param_numel / num_splats;
    int64_t numel = num_splats * stride;

    FusedAdamParams p{};
    p.param = (uint64_t)param.data_ptr();
    p.grad = (uint64_t)grad.data_ptr();
    p.exp_avg = (uint64_t)exp_avg.data_ptr();
    p.exp_avg_sq = (uint64_t)exp_avg_sq.data_ptr();
    p.steps = vkk::or_fallback(per_splat_steps.data_ptr());
    p.has_steps = per_splat_steps.data_ptr() ? 1u : 0u;
    p.lr = lr;
    p.decay = 2.0f * l2_reg / (float)numel;
    p.decay_offset = l2_reg_offset;
    p.grad_scale = grad_scale;
    p.scalar_step = step;
    p.numel = checked_u32_numel(numel, "fused_adam_step");
    p.stride = (uint32_t)stride;
    p.zero_grad = zero_grad ? 1u : 0u;
    vkk::dispatch_flat("optimizer.fused_adam_fwd", backend::vk::SpecList{},
                       numel, 256, &p, sizeof(p), &p.wgs_per_row);
}

void fused_adagrad_step(
    DeviceTensorFloatND param,
    DeviceTensorFloatND grad,
    DeviceTensorFloatND accum,
    float lr
) {
    int64_t numel = param.numel();
    if (numel == 0) return;
    AdagradParams p{};
    p.param = (uint64_t)param.data_ptr();
    p.grad = (uint64_t)grad.data_ptr();
    p.accum = (uint64_t)accum.data_ptr();
    p.lr = lr;
    p.numel = checked_u32_numel(numel, "fused_adagrad_step");
    vkk::dispatch_flat("optimizer.fused_adagrad_fwd", backend::vk::SpecList{},
                       numel, 256, &p, sizeof(p), &p.wgs_per_row);
}

void fused_adam_step_quantized(
    int64_t num_splats,
    DeviceTensorFloatND param,
    DeviceTensorFloatND grad,
    uint8_t* packed,
    float4* quant_bounds,
    float lr,
    int32_t step, DeviceVector<int32_t> per_splat_steps,
    float l2_reg,
    float l2_reg_offset,
    int bits,
    float grad_scale, bool zero_grad
) {
    int64_t param_numel = param.numel();
    if (param_numel == 0 || num_splats == 0)
        return;
    if (bits != 4 && bits != 8)
        throw std::runtime_error(
            "fused_adam_step_quantized: bits must be 4 or 8, got " +
            std::to_string(bits));
    int64_t stride = param_numel / num_splats;
    int64_t numel = num_splats * stride;

    QAdamParams p{};
    p.param = (uint64_t)param.data_ptr();
    p.grad = (uint64_t)grad.data_ptr();
    p.packed = (uint64_t)packed;
    p.quant_bounds = (uint64_t)quant_bounds;
    p.steps = vkk::or_fallback(per_splat_steps.data_ptr());
    p.has_steps = per_splat_steps.data_ptr() ? 1u : 0u;
    p.lr = lr;
    p.decay = 2.0f * l2_reg / (float)numel;
    p.decay_offset = l2_reg_offset;
    p.grad_scale = grad_scale;
    p.scalar_step = step;
    p.numel = checked_u32_numel(numel, "fused_adam_step_quantized");
    p.stride = (uint32_t)stride;
    p.zero_grad = zero_grad ? 1u : 0u;
    p.num_blocks = (uint32_t)((numel + 255) / 256);
    // Spec IDs: 0 = kOptimBits (kValueBits unused by this entry).
    vkk::dispatch_flat("optimizer.fused_adam_q",
                       backend::vk::SpecList{(uint32_t)bits}, numel, 256, &p,
                       sizeof(p), &p.wgs_per_row);
}

void fused_adam_step_quantized_value(
    int64_t num_splats,
    int64_t param_numel,
    DeviceTensorFloatND grad,
    const uint8_t* grad_q_packed,
    const float2* grad_q_bounds,
    uint8_t* optim_packed,
    float4* optim_bounds,
    uint8_t* value_packed,
    float2* value_bounds,
    float lr,
    int32_t step, DeviceVector<int32_t> per_splat_steps,
    float l2_reg,
    float l2_reg_offset,
    int optim_bits,
    int value_bits,
    float grad_scale, bool zero_grad
) {
    if (param_numel == 0 || num_splats == 0)
        return;
    if ((optim_bits != 4 && optim_bits != 8) ||
        (value_bits != 8 && value_bits != 16))
        throw std::runtime_error(
            "fused_adam_step_quantized_value: optim_bits in {4, 8} and "
            "value_bits in {8, 16}; got optim_bits=" +
            std::to_string(optim_bits) +
            ", value_bits=" + std::to_string(value_bits));
    int64_t stride = param_numel / num_splats;
    int64_t numel = num_splats * stride;

    QQAdamParams p{};
    p.grad = vkk::or_fallback(grad.data_ptr());
    p.grad_q_packed = vkk::or_fallback(grad_q_packed);
    p.grad_q_bounds = vkk::or_fallback(grad_q_bounds);
    p.grad_quant = grad_q_packed ? 1u : 0u;
    p.optim_packed = (uint64_t)optim_packed;
    p.optim_bounds = (uint64_t)optim_bounds;
    p.value_packed = (uint64_t)value_packed;
    p.value_bounds = (uint64_t)value_bounds;
    p.steps = vkk::or_fallback(per_splat_steps.data_ptr());
    p.has_steps = per_splat_steps.data_ptr() ? 1u : 0u;
    p.lr = lr;
    p.decay = 2.0f * l2_reg / (float)numel;
    p.decay_offset = l2_reg_offset;
    p.grad_scale = grad_scale;
    p.scalar_step = step;
    p.numel = checked_u32_numel(numel, "fused_adam_step_quantized_value");
    p.stride = (uint32_t)stride;
    p.zero_grad = zero_grad ? 1u : 0u;
    p.num_blocks = (uint32_t)((numel + 255) / 256);
    // Spec IDs: 0 = kOptimBits, 1 = kValueBits.
    vkk::dispatch_flat(
        "optimizer.fused_adam_qq",
        backend::vk::SpecList{(uint32_t)optim_bits, (uint32_t)value_bits},
        numel, 256, &p, sizeof(p), &p.wgs_per_row);
}

void float_add_into(DeviceVector<float> dst, DeviceVector<float> src,
                    int64_t n) {
    if (n == 0 || dst.data_ptr() == nullptr || src.data_ptr() == nullptr)
        return;
    FloatAddParams p{};
    p.dst = (uint64_t)dst.data_ptr();
    p.src = (uint64_t)src.data_ptr();
    p.n = checked_u32_numel(n, "float_add_into");
    vkk::dispatch_flat("optimizer.float_add", backend::vk::SpecList{}, n, 256,
                       &p, sizeof(p), &p.wgs_per_row);
}

void increment_int32_inplace(DeviceVector<int32_t> data, int64_t n) {
    if (n == 0 || data.data_ptr() == nullptr) return;
    IncI32Params p{};
    p.data = (uint64_t)data.data_ptr();
    p.n = checked_u32_numel(n, "increment_int32_inplace");
    vkk::dispatch_flat("optimizer.increment_i32", backend::vk::SpecList{}, n,
                       256, &p, sizeof(p), &p.wgs_per_row);
}
