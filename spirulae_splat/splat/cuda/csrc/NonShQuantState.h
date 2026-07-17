#pragma once

#include "backend/api/BackendTypes.h"
#include <cstdint>


// Bundle of non-SH Adam-state quantization buffers (means, quats, scales,
// opacities, features_dc). One QuantizedAdamState<16, 256> per attribute, in
// FPBO per-splat-block layout (one float4 bound per kFpboBlock=256 splats).
// When `enabled == false`, the kernels take the legacy fp32 g1/g2 path and
// all pointers are unused.
//
// Threaded through FPBO (decode+Adam+encode every step) and Densify (encode
// (g1=0, g2=0) bytes when a splat is relocated, against the current block
// bound).
struct NonShQuantState {
    bool     enabled            = false;
    uint8_t* means_packed       = nullptr;
    uint8_t* quats_packed       = nullptr;
    uint8_t* scales_packed      = nullptr;
    uint8_t* opacities_packed   = nullptr;
    uint8_t* features_dc_packed = nullptr;
    float4*  means_bounds       = nullptr;
    float4*  quats_bounds       = nullptr;
    float4*  scales_bounds      = nullptr;
    float4*  opacities_bounds   = nullptr;
    float4*  features_dc_bounds = nullptr;
};
