#pragma once
// Render output types used in backend API signatures: RenderOutputType,
// DistortionType (+ dist_* helpers), and RenderOutput / RenderOutput::
// TensorTuple.
//
// These live in Primitive.cuh, which already compiles as host C++ (its device
// code is #ifdef __CUDACC__-guarded) — it stays the single authoritative
// definition and this header just forwards to it. What makes that portable is
// the pass described in backend/README.md: Primitive.cuh's transitive
// cuda_runtime.h dependency (via Common.cuh / Tensor.h) is replaced by
// BackendTypes.h / BackendRuntime.h, after which this include chain is valid
// in non-CUDA builds unchanged.

#include "backend/api/BackendTypes.h"

#include "primitives/Primitive.cuh"
