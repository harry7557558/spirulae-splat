#pragma once
// Backend API module — GENERATED forwarder by
// spirulae_splat/generate_backend_api.py. DO NOT EDIT.
//
// Splat projection forward/backward, packed variants, quantized-gradient backward, and the fused projection-backward + optimizer path.
//
// Authoritative declarations live in the per-kernel headers included below
// (declaration sections generated from /*[AutoHeaderGeneratorExport]*/
// markers by generate_headers.py; CUDA-include-free via BackendTypes.h).
// Implementations: backend/cuda (the existing .cu files), backend/vulkan
// (SPIR-V pipeline dispatch). See csrc/backend/README.md.

#include "../../ProjectionFwd.cuh"
#include "../../ProjectionBwd.cuh"
#include "../../ProjectionPackedFwd.cuh"
#include "../../ProjectionBwdQuantGrad.cuh"
#include "../../FusedProjectionBwdOptim.cuh"
