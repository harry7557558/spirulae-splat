#pragma once
// Backend API module — GENERATED forwarder by
// spirulae_splat/generate_backend_api.py. DO NOT EDIT.
//
// Tile-based rasterization forward/backward (2D and eval-3D/3DGUT).
//
// Authoritative declarations live in the per-kernel headers included below
// (declaration sections generated from /*[AutoHeaderGeneratorExport]*/
// markers by generate_headers.py; CUDA-include-free via BackendTypes.h).
// Implementations: backend/cuda (the existing .cu files), backend/vulkan
// (SPIR-V pipeline dispatch). See csrc/backend/README.md.

#include "../../RasterizationFwd.cuh"
#include "../../RasterizationBwd.cuh"
#include "../../RasterizationEval3DFwd.cuh"
#include "../../RasterizationEval3DBwd.cuh"
#include "../../RasterizationMomentsFwd.cuh"
