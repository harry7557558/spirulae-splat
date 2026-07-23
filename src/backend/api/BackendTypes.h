#pragma once
// Portable replacements for the CUDA vector types and function-scope macros
// that appear in backend API signatures and header-level host code.
//
// CUDA is the default backend: unless SSPLAT_BACKEND_VULKAN is defined, the
// CUDA definitions are authoritative and this header defers to them. The
// stand-alone definitions below MUST match the CUDA layouts and alignments
// exactly — they are ABI, shared with Slang-generated SPIR-V structured
// buffers.

#ifndef SSPLAT_BACKEND_VULKAN

#include <cuda_runtime.h>

#else  // non-CUDA backend

// cuda_runtime.h defines these as empty macros in host translation units;
// header-level code (Tensor.h, Primitive.cuh) relies on that.
#define __device__
#define __host__
#define __forceinline__ inline
// __restrict__ is understood natively by GCC/Clang, but MSVC spells it
// `__restrict`, so map it there.
#ifdef _MSC_VER
#define __restrict__ __restrict
#endif

// Alignments per CUDA vector_types.h: 2- and 4-component types align to
// their full size (uchar4 -> 4, float2 -> 8, float4/int4/uint4 -> 16);
// 3-component types have no extra alignment.
struct alignas(8)  float2  { float x, y; };
struct             float3  { float x, y, z; };
struct alignas(16) float4  { float x, y, z, w; };
struct alignas(8)  int2    { int x, y; };
struct             int3    { int x, y, z; };
struct alignas(16) int4    { int x, y, z, w; };
struct alignas(8)  uint2   { unsigned x, y; };
struct             uint3   { unsigned x, y, z; };
struct alignas(16) uint4   { unsigned x, y, z, w; };
struct alignas(2)  uchar2  { unsigned char x, y; };
struct             uchar3  { unsigned char x, y, z; };
struct alignas(4)  uchar4  { unsigned char x, y, z, w; };
struct alignas(4)  ushort2 { unsigned short x, y; };
struct             ushort3 { unsigned short x, y, z; };
struct alignas(8)  ushort4 { unsigned short x, y, z, w; };
struct alignas(16) double2 { double x, y; };

inline float2 make_float2(float x, float y) { return {x, y}; }
inline float3 make_float3(float x, float y, float z) { return {x, y, z}; }
inline float4 make_float4(float x, float y, float z, float w) { return {x, y, z, w}; }
inline int2   make_int2(int x, int y) { return {x, y}; }
inline int3   make_int3(int x, int y, int z) { return {x, y, z}; }
inline int4   make_int4(int x, int y, int z, int w) { return {x, y, z, w}; }
inline uint2  make_uint2(unsigned x, unsigned y) { return {x, y}; }
inline uint3  make_uint3(unsigned x, unsigned y, unsigned z) { return {x, y, z}; }
inline uint4  make_uint4(unsigned x, unsigned y, unsigned z, unsigned w) { return {x, y, z, w}; }
inline uchar4 make_uchar4(unsigned char x, unsigned char y, unsigned char z, unsigned char w) { return {x, y, z, w}; }

#endif  // backend select
