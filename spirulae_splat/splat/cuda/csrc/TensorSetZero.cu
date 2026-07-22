// TensorSetZero.cu -- set_zero_tensor: zero an arbitrary float tensor.
//
// Part of the Optimizer family -- see OptimizerCommon.cuh.

#include "OptimizerCommon.cuh"

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
