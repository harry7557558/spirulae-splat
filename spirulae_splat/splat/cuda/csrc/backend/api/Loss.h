#pragma once
// Backend API module — GENERATED forwarder by
// spirulae_splat/generate_backend_api.py. DO NOT EDIT.
//
// Fused SSIM, multi-scale per-pixel losses, per-splat losses.
//
// Authoritative declarations live in the per-kernel headers included below
// (declaration sections generated from /*[AutoHeaderGeneratorExport]*/
// markers by generate_headers.py; CUDA-include-free via BackendTypes.h).
// Implementations: backend/cuda (the existing .cu files), backend/vulkan
// (SPIR-V pipeline dispatch). See csrc/backend/README.md.

#include "../../FusedSSIM.cuh"
#include "../../PerPixelLoss.cuh"
#include "../../PerSplatLoss.cuh"
