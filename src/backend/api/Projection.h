#pragma once
// Backend API module — GENERATED forwarder by
// tools/codegen/generate_backend_api.py. DO NOT EDIT.
//
// Splat projection forward/backward, packed variants, quantized-gradient backward, and the fused projection-backward + optimizer path.
//
// Authoritative declarations live in the per-kernel headers included below
// (declaration sections generated from /*[AutoHeaderGeneratorExport]*/
// markers by generate_headers.py; CUDA-include-free via BackendTypes.h).
// Implementations: the .cu files next to each header (CUDA) and
// backend/vulkan/kernels/ (SPIR-V pipeline dispatch). See backend/README.md.

#include "kernels/projection/ProjectionFwd.cuh"
#include "kernels/projection/ProjectionBwd.cuh"
#include "kernels/projection/ProjectionPackedFwd.cuh"
#include "kernels/projection/ProjectionBwdQuantGrad.cuh"
#include "kernels/optim/FusedProjectionBwdOptim.cuh"
