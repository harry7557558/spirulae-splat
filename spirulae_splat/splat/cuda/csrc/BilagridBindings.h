#pragma once

// Forward declarations of all bilagrid CUDA launchers ported from
// fused-bilagrid/fused_bilagrid/bindings.h. The original "tensor" wrappers used
// torch::Tensor; this header replaces them with TorchTensorView and pre-allocated
// output buffers (the caller is responsible for allocation).

#include <cuda_runtime.h>
#include <cstdint>

#include "Tensor.h"
#include "BilagridReader.cuh"


// ============================================================================
// Raw launchers (operate on float* / int* pointers).
// ============================================================================

void bilagrid_sample_forward(
    BilagridReader bilagrid,
    const float* coords,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
);

void bilagrid_sample_backward(
    BilagridReader bilagrid,
    const float* coords,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_coords,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
);

void bilagrid_uniform_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int h, int w,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_patched_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
);

void bilagrid_uniform_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int h, int w,
    const int target_tile_size,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_patched_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    const int target_tile_size,
    const int mi_batch_size,
    cudaStream_t stream
);

void bilagrid_uniform_sample_backward_v2(
    BilagridReader bilagrid,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int h, int w,
    cudaStream_t stream
);

void bilagrid_patched_sample_backward_v2(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
);

void bilagrid_ppisp_sample_forward(
    BilagridReader bilagrid,
    const float* coords,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
);

void bilagrid_ppisp_sample_backward(
    BilagridReader bilagrid,
    const float* coords,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_coords,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
);

void bilagrid_ppisp_packed_sample_forward(
    BilagridReader bilagrid,
    const int64_t* image_indices,
    const float* coords,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int nnz,
    cudaStream_t stream
);

void bilagrid_ppisp_packed_sample_backward(
    BilagridReader bilagrid,
    const int64_t* image_indices,
    const float* coords,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_coords,
    float* v_rgb,
    int N, int L, int H, int W,
    int nnz,
    cudaStream_t stream
);

void bilagrid_ppisp_uniform_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int h, int w,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_ppisp_patched_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
);

void bilagrid_ppisp_uniform_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int h, int w,
    const int target_tile_size,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_ppisp_patched_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    const int target_tile_size,
    const int mi_batch_size,
    cudaStream_t stream
);

void bilagrid_loglinear_uniform_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int h, int w,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_loglinear_patched_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
);

void bilagrid_loglinear_uniform_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int h, int w,
    const int target_tile_size,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_loglinear_patched_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    const int target_tile_size,
    const int mi_batch_size,
    cudaStream_t stream
);

void bilagrid_depth_uniform_sample_forward(
    BilagridReader bilagrid,
    const float* depth,
    const float* scalars,
    float* output,
    int N, int L, int H, int W,
    int h, int w,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_depth_patched_sample_forward(
    BilagridReader bilagrid,
    const float* depth,
    const float* scalars,
    const int* offsets,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
);

void bilagrid_depth_uniform_sample_backward_v1(
    BilagridReader bilagrid,
    const float* depth,
    const float* scalars,
    const float* v_output,
    float* v_bilagrid,
    float* v_depth,
    int N, int L, int H, int W,
    int h, int w,
    const int target_tile_size,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_depth_patched_sample_backward_v1(
    BilagridReader bilagrid,
    const float* depth,
    const float* scalars,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_depth,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    const int target_tile_size,
    const int mi_batch_size,
    cudaStream_t stream
);

void bilagrid_normal_uniform_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    float* output,
    int N, int L, int H, int W,
    int h, int w,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_normal_patched_sample_forward(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
);

void bilagrid_normal_uniform_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int h, int w,
    const int target_tile_size,
    cudaStream_t stream,
    const int* grid_indices = nullptr
);

void bilagrid_normal_patched_sample_backward_v1(
    BilagridReader bilagrid,
    const float* rgb,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_rgb,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    const int target_tile_size,
    const int mi_batch_size,
    cudaStream_t stream
);

void tv_loss_forward(
    BilagridReader bilagrid,
    float* tv_loss,
    int N, int C, int L, int H, int W,
    cudaStream_t stream
);

void tv_loss_backward(
    BilagridReader bilagrid,
    const float v_tv_loss,
    float* v_bilagrid,
    int N, int C, int L, int H, int W, bool inplace,
    cudaStream_t stream
);

void channel_mean_forward(
    BilagridReader bilagrid,
    float* channel_mean,
    int N, int C, int L, int H, int W,
    cudaStream_t stream
);

void channel_mean_backward(
    BilagridReader bilagrid,
    const float* v_channel_mean,
    float* v_bilagrid,
    int N, int C, int L, int H, int W, bool inplace,
    cudaStream_t stream
);


// ============================================================================
// TorchTensorView wrappers (ports of bindings.h's torch::Tensor wrappers).
// Outputs are caller-allocated; tensor shape conventions follow the
// corresponding original wrappers.
// ============================================================================

void bilagrid_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,12,L,H,W]
    TorchTensorView coords,    // [N,m,h,w,2]
    TorchTensorView rgb,       // [N,m,h,w,3]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,12,L,H,W]
    TorchTensorView coords,      // [N,m,h,w,2]
    TorchTensorView rgb,         // [N,m,h,w,3]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,12,L,H,W] (caller pre-zeroed)
    TorchTensorView v_coords,    // [N,m,h,w,2] (caller pre-zeroed; null to skip)
    TorchTensorView v_rgb        // [N,m,h,w,3]
);

void bilagrid_ppisp_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,9,L,H,W]
    TorchTensorView coords,    // [N,m,h,w,2]
    TorchTensorView rgb,       // [N,m,h,w,3]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_ppisp_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,9,L,H,W]
    TorchTensorView coords,      // [N,m,h,w,2]
    TorchTensorView rgb,         // [N,m,h,w,3]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,9,L,H,W]
    TorchTensorView v_coords,    // [N,m,h,w,2] (null to skip)
    TorchTensorView v_rgb        // [N,m,h,w,3]
);

void bilagrid_ppisp_packed_sample_forward_tensor(
    TorchTensorView bilagrid,        // [N,9,L,H,W]
    TorchTensorView image_indices,   // [nnz]   int64
    TorchTensorView coords,          // [nnz,2]
    TorchTensorView rgb,             // [nnz,3]
    TorchTensorView output           // [nnz,3]
);

void bilagrid_ppisp_packed_sample_backward_tensor(
    TorchTensorView bilagrid,        // [N,9,L,H,W]
    TorchTensorView image_indices,   // [nnz]
    TorchTensorView coords,          // [nnz,2]
    TorchTensorView rgb,             // [nnz,3]
    TorchTensorView v_output,        // [nnz,3]
    TorchTensorView v_bilagrid,      // [N,9,L,H,W]
    TorchTensorView v_coords,        // [nnz,2] (null to skip)
    TorchTensorView v_rgb            // [nnz,3]
);

void bilagrid_uniform_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,12,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_uniform_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,12,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,12,L,H,W]
    TorchTensorView v_rgb,       // [N,m,h,w,3]
    int version,
    int target_tile_size
);

void bilagrid_ppisp_uniform_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,9,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_ppisp_uniform_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,9,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,9,L,H,W]
    TorchTensorView v_rgb,       // [N,m,h,w,3]
    int version,
    int target_tile_size
);


void bilagrid_loglinear_uniform_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,9,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_loglinear_uniform_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,9,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,9,L,H,W]
    TorchTensorView v_rgb,       // [N,m,h,w,3]
    int version,
    int target_tile_size
);

void compute_depth_scalars_tensor(
    TorchTensorView depth,    // [N,m,h,w,1]
    bool patched,
    TorchTensorView output    // [B] = [1] when patched else [N]
);

void bilagrid_depth_uniform_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,2,L,H,W]
    TorchTensorView depth,     // [N,m,h,w,1]
    TorchTensorView scalars,   // [N]
    TorchTensorView output     // [N,m,h,w,1]
);

void bilagrid_depth_uniform_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,2,L,H,W]
    TorchTensorView depth,       // [N,m,h,w,1]
    TorchTensorView scalars,     // [N]
    TorchTensorView v_output,    // [N,m,h,w,1]
    TorchTensorView v_bilagrid,  // [N,2,L,H,W]
    TorchTensorView v_depth,     // [N,m,h,w,1]
    int version,
    int target_tile_size
);

void bilagrid_normal_uniform_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,3,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_normal_uniform_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,3,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,3,L,H,W]
    TorchTensorView v_rgb,       // [N,m,h,w,3]
    int version,
    int target_tile_size
);

void bilagrid_patched_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,12,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,   // [N,m,2] int32
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_patched_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,12,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,     // [N,m,2] int32
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,12,L,H,W]
    TorchTensorView v_rgb        // [N,m,h,w,3]
);

void bilagrid_ppisp_patched_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,9,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,   // [N,m,2]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_ppisp_patched_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,9,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,     // [N,m,2]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,9,L,H,W]
    TorchTensorView v_rgb        // [N,m,h,w,3]
);

void bilagrid_loglinear_patched_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,9,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,   // [N,m,2]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_loglinear_patched_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,9,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,     // [N,m,2]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,9,L,H,W]
    TorchTensorView v_rgb        // [N,m,h,w,3]
);

void bilagrid_depth_patched_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,2,L,H,W]
    TorchTensorView depth,     // [N,m,h,w,1]
    TorchTensorView scalars,   // [N]
    int h0, int w0,
    TorchTensorView offsets,   // [N,m,2]
    TorchTensorView output     // [N,m,h,w,1]
);

void bilagrid_depth_patched_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,2,L,H,W]
    TorchTensorView depth,       // [N,m,h,w,1]
    TorchTensorView scalars,     // [N]
    int h0, int w0,
    TorchTensorView offsets,     // [N,m,2]
    TorchTensorView v_output,    // [N,m,h,w,1]
    TorchTensorView v_bilagrid,  // [N,2,L,H,W]
    TorchTensorView v_depth      // [N,m,h,w,1]
);

void bilagrid_normal_patched_sample_forward_tensor(
    TorchTensorView bilagrid,  // [N,3,L,H,W]
    TorchTensorView rgb,       // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,   // [N,m,2]
    TorchTensorView output     // [N,m,h,w,3]
);

void bilagrid_normal_patched_sample_backward_tensor(
    TorchTensorView bilagrid,    // [N,3,L,H,W]
    TorchTensorView rgb,         // [N,m,h,w,3]
    int h0, int w0,
    TorchTensorView offsets,     // [N,m,2]
    TorchTensorView v_output,    // [N,m,h,w,3]
    TorchTensorView v_bilagrid,  // [N,3,L,H,W]
    TorchTensorView v_rgb        // [N,m,h,w,3]
);

void tv_loss_forward_tensor(
    TorchTensorView bilagrid,  // [N,C,L,H,W]
    TorchTensorView tv_loss    // scalar [] (caller pre-zeroed)
);

void tv_loss_backward_tensor(
    TorchTensorView bilagrid,    // [N,C,L,H,W]
    float v_tv_loss,             // scalar
    TorchTensorView v_bilagrid,  // [N,C,L,H,W]
    bool inplace                 // if true, accumulates into v_bilagrid
);

void channel_mean_forward_tensor(
    TorchTensorView bilagrid,      // [N,C,L,H,W]
    TorchTensorView channel_mean   // [C] (caller pre-zeroed)
);

void channel_mean_backward_tensor(
    TorchTensorView bilagrid,        // [N,C,L,H,W]
    TorchTensorView v_channel_mean,  // [C]
    TorchTensorView v_bilagrid,      // [N,C,L,H,W]
    bool inplace
);
