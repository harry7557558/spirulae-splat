#pragma once

#ifndef NO_TORCH
#include <torch/types.h>
#endif

#include <gsplat/Common.h>


#ifdef __CUDACC__
#include "generated/slang.cuh"
#endif

#ifndef __CUDACC__

template<typename T, size_t SIZE>
struct FixedArray
{
    T m_data[SIZE];
};

#endif



#ifndef NO_TORCH
typedef std::optional<at::Tensor> CameraDistortionCoeffsTensor;
#endif

#ifdef __CUDACC__
// k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2
typedef FixedArray<float, 10> CameraDistortionCoeffs;
#endif

struct CameraDistortionCoeffsBuffer {
    float* __restrict__ coeffs;

    #ifndef NO_TORCH
    CameraDistortionCoeffsBuffer(const CameraDistortionCoeffsTensor &tensors);
    #endif

    #ifdef __CUDACC__
    __device__ CameraDistortionCoeffs load(long idx) const {
        CameraDistortionCoeffs res;
        if (coeffs == nullptr) {
            #pragma unroll
            for (int i = 0; i < 10; i++)
                res[i] = 0.0f;
        } else {
            float* c = coeffs + 10 * idx;
            #pragma unroll
            for (int i = 0; i < 10; i++)
                res[i] = c[i];
        }
        return res;
    }
    #endif
};

enum class HessianDiagonalOutputMode {
    None,
    Position,
    AllReasonable
};
