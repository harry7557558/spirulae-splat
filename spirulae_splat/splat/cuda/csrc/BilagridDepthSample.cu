#include "BilagridDepthSampleFwd_kernel.cuh"
#include "BilagridDepthSampleBwdV1_kernel.cuh"
// #include "uniform_sample_depth_backward_v2.cu"

#define PATCHED
#include "BilagridDepthSampleFwd_kernel.cuh"
#include "BilagridDepthSampleBwdV1_kernel.cuh"
// #include "uniform_sample_depth_backward_v2.cu"


void bilagrid_depth_uniform_sample_forward(
    const float* bilagrid,
    const float* depth,
    const float* scalars,
    float* output,
    int N, int L, int H, int W,
    int h, int w,
    cudaStream_t stream,
    const int* grid_indices
) {
    int total = N * h * w;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    bilagrid_depth_uniform_sample_forward_kernel<<<blocks, threads, 0, stream>>>(
        bilagrid, depth, scalars, output,
        N, L, H, W, h, w,
        grid_indices
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    // cudaDeviceSynchronize();
}


void bilagrid_depth_patched_sample_forward(
    const float* bilagrid,
    const float* depth,
    const float* scalars,
    const int* offsets,
    float* output,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
) {
    int total = N * m * h * w;
    int threads = 256;
    int blocks = (total + threads - 1) / threads;
    bilagrid_depth_patched_sample_forward_kernel<<<blocks, threads, 0, stream>>>(
        bilagrid, depth, scalars, output,
        N, L, H, W, m, h, w, h0, w0, offsets
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
    // cudaDeviceSynchronize();
}


void bilagrid_depth_uniform_sample_backward_v1(
    const float* bilagrid,
    const float* depth,
    const float* scalars,
    const float* v_output,
    float* v_bilagrid,
    float* v_depth,            // null: skip the v_pre kernel (GT isn't trained)
    int N, int L, int H, int W,
    int h, int w,
    const unsigned block_x, const unsigned block_y,
    const int target_tile_size,
    cudaStream_t stream,
    const int* grid_indices
) {
    // v_bilagrid (always needed: trains the bilagrid depth grid)
    {
        dim3 block = { block_x, block_y, 1 };

        int mult_x = (2*w+W)/(block.x*W*target_tile_size);
        int mult_y = (2*h+H)/(block.y*H*target_tile_size);
        if (mult_x * mult_y < 4)
            mult_x = mult_y = 1;
        else {
            mult_x = max(mult_x, 1) * block.x;
            mult_y = max(mult_y, 1) * block.y;
        }

        dim3 bounds = {
            (W*mult_x +block.x-1)/block.x,
            (H*mult_y +block.y-1)/block.y,
            (N*L +block.z-1)/block.z
        };
        bilagrid_depth_uniform_sample_backward_v1_kernel_bilagrid<<<bounds, block, 0, stream>>>(
            bilagrid, depth, scalars, v_output, v_bilagrid,
            N, L, H, W, h, w, mult_x, mult_y,
            grid_indices
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    // v_depth: gradient w.r.t. pre-bilagrid input depth (i.e., the raw GT).
    // Skipped when the caller passes null — GT isn't a trainable parameter.
    if (v_depth != nullptr) {
        int total = N * h * w;
        int threads = 256;
        int blocks = (total + threads - 1) / threads;
        bilagrid_depth_uniform_sample_backward_v1_kernel_depth<<<blocks, threads, 0, stream>>>(
            bilagrid, depth, scalars, v_output,
            v_depth,
            N, L, H, W, h, w,
            grid_indices
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
}


void bilagrid_depth_patched_sample_backward_v1(
    const float* bilagrid,
    const float* depth,
    const float* scalars,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_depth,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    const unsigned block_x, const unsigned block_y,
    const int target_tile_size,
    const int mi_batch_size,
    cudaStream_t stream
) {
    // v_bilagrid
    {
        dim3 block = { block_x, block_y, 1 };
    
        int mult_x = (2*w0+W)/(block.x*W*target_tile_size);
        int mult_y = (2*h0+H)/(block.y*H*target_tile_size);
        if (mult_x * mult_y < 4)
            mult_x = mult_y = 1;
        else {
            mult_x = max(mult_x, 1) * block.x;
            mult_y = max(mult_y, 1) * block.y;
        }
        // printf("mult_x: %d, mult_y: %d\n", mult_x, mult_y);

        // int W1 = min(W, (W*w + w0-1)/w0 + 3);
        // int H1 = min(H, (H*h + h0-1)/h0 + 3);
        int W1 = W, H1 = H;

        int num_m_batches = (m+mi_batch_size-1) / mi_batch_size;

        dim3 bounds = {
            (W1*mult_x +block.x-1)/block.x,
            (H1*mult_y +block.y-1)/block.y,
            (N*num_m_batches*L +block.z-1)/block.z
        };
        // printf("bounds: %u %u %u\n", bounds.x, bounds.y, bounds.z);
        bilagrid_depth_patched_sample_backward_v1_kernel_bilagrid<<<bounds, block, 0, stream>>>(
            bilagrid, depth, scalars, v_output, v_bilagrid,
            N, L, H, W, m, h, w, h0, w0, offsets, mult_x, mult_y, num_m_batches
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }

    // v_depth
    {
        int total = N * m * h * w;
        int threads = 256;
        int blocks = (total + threads - 1) / threads;
        bilagrid_depth_patched_sample_backward_v1_kernel_depth<<<blocks, threads, 0, stream>>>(
            bilagrid, depth, scalars, v_output, v_depth,
            N, L, H, W, m, h, w, h0, w0, offsets
        );
        CHECK_DEVICE_ERROR(cudaGetLastError());
    }
}


#if 0

void bilagrid_depth_uniform_sample_backward_v2(
    const float* bilagrid,
    const float* depth,
    const float* v_output,
    float* v_bilagrid,
    float* v_depth,
    int N, int L, int H, int W,
    int m, int h, int w,
    cudaStream_t stream
) {
    dim3 block = { 16, 16, 1 };
    dim3 bounds = {
        (w +block.x-1)/block.x,
        (h +block.y-1)/block.y,
        (N*m +block.z-1)/block.z
    };
    bilagrid_depth_uniform_sample_backward_v2_kernel<<<bounds, block, 0, stream>>>(
        bilagrid, depth, v_output,
        v_bilagrid, v_depth,
        N, L, H, W, m, h, w
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


void bilagrid_patched_sample_backward_v2(
    const float* bilagrid,
    const float* depth,
    const int* offsets,
    const float* v_output,
    float* v_bilagrid,
    float* v_depth,
    int N, int L, int H, int W,
    int m, int h, int w, int h0, int w0,
    cudaStream_t stream
) {
    // dim3 block = { 16, 16, 1 };
    // dim3 bounds = {
    //     (w +block.x-1)/block.x,
    //     (h +block.y-1)/block.y,
    //     (N*m +block.z-1)/block.z
    // };
    unsigned block = 256;
    unsigned bounds = (w*h*N*m +block-1)/block;
    bilagrid_depth_patched_sample_backward_v2_kernel<<<bounds, block, 0, stream>>>(
        bilagrid, depth, v_output,
        v_bilagrid, v_depth,
        N, L, H, W, m, h, w, h0, w0, offsets
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}

#endif
