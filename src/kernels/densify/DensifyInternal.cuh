#pragma once

// Densify-internal primitives implemented in DensifySampling.cu and reused by
// the other densify translation units. These are not part of the public
// Densify.cuh surface (no /*[AutoHeaderGeneratorExport]*/) -- they exist only
// so the parts of the densify family can share one implementation instead of
// each carrying a copy.

#include "kernels/densify/DensifyCommon.cuh"

// Quantile over |x| across the finite elements of each of B rows of N.
// Writes B values; `return_reciprocal` emits 1/q instead of q.
// Uses PoolSlot::DensifyQuantileTemp scratch.
void quantile_of_abs_of_finite_elements_internal(
    const float* inputs_ptr,
    int B,
    int N,
    float q,
    bool return_reciprocal,
    float* outputs_ptr
);

// Weighted sampling without replacement (Efraimidis-Spirakis). Returns a
// pool-allocated int32 index buffer of length `num_sample`; the caller does
// not own it (PoolSlot-backed).
int32_t* weighted_sample_without_replacement_internal(
    int64_t numel,
    float* weights_ptr,
    int64_t weights_numel,
    bool* masks_ptr,
    uint32_t num_sample,
    uint32_t seed
);
