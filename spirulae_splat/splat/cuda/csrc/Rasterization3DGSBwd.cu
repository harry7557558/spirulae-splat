// Modified from https://github.com/nerfstudio-project/gsplat/blob/main/gsplat/cuda/csrc/RasterizeToPixels3DGSBwd.cu

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include <cub/cub.cuh>

#include "Primitive.cuh"


constexpr uint BLOCK_SIZE = TILE_SIZE * TILE_SIZE;
constexpr uint SPLAT_BATCH_SIZE = 128;


template <uint32_t CDIM>
__global__ void rasterize_to_pixels_3dgs_bwd_kernel(
    const uint32_t I,
    const uint32_t N,
    const uint32_t n_isects,
    const bool packed,
    // fwd inputs
    Vanilla3DGS::Screen::Buffer splat_buffer,
    const float *__restrict__ colors,      // [..., N, CDIM] or [nnz, CDIM]
    const float *__restrict__ backgrounds, // [..., CDIM] or [nnz, CDIM]
    const bool *__restrict__ masks,           // [..., tile_height, tile_width]
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_width,
    const uint32_t tile_height,
    const int32_t *__restrict__ tile_offsets, // [..., tile_height, tile_width]
    const int32_t *__restrict__ flatten_ids,  // [n_isects]
    // fwd outputs
    const float
        *__restrict__ render_Ts,      // [..., image_height, image_width, 1]
    const int32_t *__restrict__ last_ids, // [..., image_height, image_width]
    // grad outputs
    const float *__restrict__ v_render_colors, // [..., image_height,
                                                  // image_width, CDIM]
    const float
        *__restrict__ v_render_alphas, // [..., image_height, image_width, 1]
    // grad inputs
    Vanilla3DGS::Screen::Buffer v_splat_buffer,
    float *__restrict__ v_colors   // [..., N, CDIM] or [nnz, CDIM]
) {
    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);
    uint32_t image_id = block.group_index().x;
    uint32_t tile_id = block.group_index().y * tile_width + block.group_index().z;
    uint32_t thread_id = block.thread_rank();

    tile_offsets += image_id * tile_height * tile_width;
    render_Ts += image_id * image_height * image_width;
    last_ids += image_id * image_height * image_width;
    v_render_colors += image_id * image_height * image_width * CDIM;
    v_render_alphas += image_id * image_height * image_width;
    if (backgrounds != nullptr) {
        backgrounds += image_id * CDIM;
    }
    if (masks != nullptr) {
        masks += image_id * tile_height * tile_width;
    }

    // when the mask is provided, do nothing and return if
    // this tile is labeled as False
    if (masks != nullptr && !masks[tile_id]) {
        return;
    }

    // load pixels
    __shared__ int32_t pix_bin_final[BLOCK_SIZE];
    __shared__ float2 pix_Ts_with_grad[BLOCK_SIZE];
    __shared__ float v_pix_colors[BLOCK_SIZE*CDIM];
    #pragma unroll
    for (uint pix_id0 = 0; pix_id0 < BLOCK_SIZE; pix_id0 += SPLAT_BATCH_SIZE) {
        static_assert(BLOCK_SIZE % SPLAT_BATCH_SIZE == 0);
        uint pix_id_local = pix_id0 + thread_id;
        int pix_x = block.group_index().z * TILE_SIZE + pix_id_local % TILE_SIZE;
        int pix_y = block.group_index().y * TILE_SIZE + pix_id_local / TILE_SIZE;
        uint pix_id_global = pix_y * image_width + pix_x;
        bool inside = (pix_x < image_width && pix_y < image_height);
        
        int32_t bin_final = (inside ? last_ids[pix_id_global] : 0);
        pix_bin_final[pix_id_local] = bin_final;
        pix_Ts_with_grad[pix_id_local] = {
            (inside ? render_Ts[pix_id_global] : 0.0f),
            (inside ? -v_render_alphas[pix_id_global] : 0.0f)
        };
        #pragma unroll
        for (uint k = 0; k < CDIM; k++) {
            v_pix_colors[pix_id_local * CDIM + k] =
                (inside ? v_render_colors[pix_id_global * CDIM + k] : 0.0f);
        }
    }
    static_assert(CDIM <= SPLAT_BATCH_SIZE);
    __shared__ float pix_background[CDIM];
    if (thread_id < CDIM && backgrounds != nullptr)
        pix_background[thread_id] = backgrounds[thread_id];
    block.sync();

    // threads fist load splats, then swept through pixels
    // do this in batches

    int32_t range_start = tile_offsets[tile_id];
    int32_t range_end =
        (image_id == I - 1) && (tile_id == tile_width * tile_height - 1)
            ? n_isects
            : tile_offsets[tile_id + 1];
    const uint32_t num_splat_batches =
        _CEIL_DIV(range_end - range_start, SPLAT_BATCH_SIZE);

    // if (warp.thread_rank() == 0)
    //     printf("range_start=%d range_end=%d num_splat_batches=%u\n", range_start, range_end, num_splat_batches);
    for (uint32_t splat_b = 0; splat_b < num_splat_batches; ++splat_b) {
        const int32_t splat_batch_end = range_end - 1 - SPLAT_BATCH_SIZE * splat_b;
        const int32_t splat_batch_size = min(SPLAT_BATCH_SIZE, splat_batch_end + 1 - range_start);
        const int32_t splat_idx = splat_batch_end - thread_id;

        // load splats
        Vanilla3DGS::Screen splat;
        uint32_t splat_gid;
        float splat_color[CDIM];
        if (splat_idx >= range_start) {
            splat_gid = flatten_ids[splat_idx]; // flatten index in [I * N] or [nnz]
            splat = Vanilla3DGS::Screen::load(splat_buffer, splat_gid);
            #pragma unroll
            for (uint32_t k = 0; k < CDIM; ++k)
                splat_color[k] = colors[splat_gid * CDIM + k];
        }

        // accumulate gradient
        float v_rgb_local[CDIM] = {0.f};
        #pragma unroll
        for (uint32_t k = 0; k < CDIM; ++k)
            v_rgb_local[k] = 0.f;
        Vanilla3DGS::Screen v_splat = Vanilla3DGS::Screen::zero();

        // thread 0 takes last splat, 1 takes second last, etc.
        // at t=0, thread 0 (splat -1) undo pixel 0
        // at t=1, thread 0 (splat -1) undo pixel 1, thread 1 (splat -2) undo pixel 0
        // ......

        // process gaussians in the current batch for this pixel
        // 0 index is the furthest back gaussian in the batch
        for (int t = 0; t < splat_batch_size + BLOCK_SIZE - 1; ++t) {
            int pix_id = t - thread_id;

            int pix_local_x = pix_id % TILE_SIZE;
            int pix_local_y = pix_id / TILE_SIZE;
            int pix_global_x = block.group_index().z * TILE_SIZE + pix_local_x;
            int pix_global_y = block.group_index().y * TILE_SIZE + pix_local_y;
            const float px = (float)pix_global_x + 0.5f;
            const float py = (float)pix_global_y + 0.5f;
            if ((pix_id >= 0 && pix_id < BLOCK_SIZE &&
                pix_global_x < image_width && pix_global_y < image_height &&
                splat_idx >= range_start) &&
                splat_idx <= pix_bin_final[pix_id]
            ) {

            // evaluate alpha and early skip
            float alpha = splat.evaluate_alpha(px, py);
            if (alpha >= ALPHA_THRESHOLD) {

            // printf("t=%d, thread %u, splat %d (%u), pix_id %d, pix %d %d\n", t, thread_id, splat_idx-range_start, splat_gid, pix_id, pix_global_x, pix_global_y);

            // forward:
            // \left(c_{1},T_{1}\right)=\left(c_{0}+\alpha_{i}T_{0}c_{i},\ T_{0}\left(1-\alpha_{i}\right)\right)
            float T1 = pix_Ts_with_grad[pix_id].x;
            float v_T1 = pix_Ts_with_grad[pix_id].y;

            // undo pixel:
            // T_{0}=\frac{T_{1}}{1-\alpha_{i}}
            float ra = 1.0f / (1.0f - alpha);
            float T0 = T1 * ra;

            // gradient to alpha:
            // \frac{dL}{d\alpha_{i}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{d\alpha_{i}}+\frac{dL}{dT_{1}}\frac{dT_{1}}{d\alpha_{i}}
            // = T_{0}\frac{dL}{dc_{1}}c_{i}-\frac{dL}{dT_{1}}T_{0}

            // gradient to color:
            // \frac{dL}{dc_{i}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{dc_{i}}
            // = \alpha_{i}T_{0}\frac{dL}{dc_{1}}

            // update pixel gradient:
            // \frac{dL}{dT_{0}}
            // = \frac{dL}{dc_{1}}\frac{dc_{1}}{dT_{0}}+\frac{dL}{dT_{1}}\frac{dT_{1}}{dT_{0}}
            // = \alpha_{i}\frac{dL}{dc_{1}}c_{i}+\frac{dL}{dT_{1}}\left(1-\alpha_{i}\right)

            float v_alpha = -v_T1 * T0;  // gradient to alpha
            float v_T0 = v_T1 * (1.0f - alpha);  // update pixel gradient
            #pragma unroll
            for (uint32_t k = 0; k < CDIM; ++k) {
                float c = splat_color[k];
                float v_c = v_pix_colors[pix_id * CDIM + k];
                v_alpha += c * v_c * T0;  // gradient to alpha
                v_rgb_local[k] += alpha * T0 * v_c;  // gradient to color
                v_T0 += c * v_c * alpha; // update pixel gradient
            }

            // backward diff splat
            v_splat += splat.evaluate_alpha_vjp(px, py, v_alpha);

            // update pixel states
            pix_Ts_with_grad[pix_id] = { T0, v_T0 };
            // v_pix_colors remains the same

            }}
            block.sync();
        }

        // accumulate gradient
        {
            if (splat_idx >= range_start) {

                float *v_rgb_ptr = (float *)(v_colors) + CDIM * splat_gid;
                #pragma unroll
                for (uint32_t k = 0; k < CDIM; ++k) {
                    if (v_rgb_local[k] != 0.0f)
                        atomicAdd(v_rgb_ptr + k, v_rgb_local[k]);
                }

                v_splat.atomicAddBuffer(v_splat_buffer, splat_gid);
            }
        }
    }
}


template <uint32_t CDIM>
void launch_rasterize_to_pixels_3dgs_bwd_kernel(
    // Gaussian parameters
    Vanilla3DGS::Screen::Tensor splats,
    const at::Tensor colors,                    // [..., N, 3] or [nnz, 3]
    const std::optional<at::Tensor> backgrounds, // [..., 3]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    // gradients of outputs
    const at::Tensor v_render_colors, // [..., image_height, image_width, 3]
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    // outputs
    Vanilla3DGS::Screen::Tensor v_splats,
    at::Tensor v_colors                    // [..., N, 3] or [nnz, 3]
) {
    bool packed = splats.isPacked();
    uint32_t N = packed ? 0 : splats.size(); // number of gaussians
    uint32_t I = render_Ts.numel() / (image_height * image_width); // number of images
    uint32_t tile_height = tile_offsets.size(-2);
    uint32_t tile_width = tile_offsets.size(-1);
    uint32_t n_isects = flatten_ids.size(0);

    // Each block covers a tile on the image. In total there are
    // I * tile_height * tile_width blocks.
    dim3 threads = {SPLAT_BATCH_SIZE, 1, 1};
    dim3 grid = {I, tile_height, tile_width};

    // int64_t shmem_size =
    //     TILE_SIZE * TILE_SIZE *
    //     (sizeof(int32_t) + sizeof(glm::vec3) + sizeof(glm::vec3) + sizeof(float) * CDIM);

    if (n_isects == 0) {
        // skip the kernel launch if there are no elements
        return;
    }

    // TODO: an optimization can be done by passing the actual number of
    // channels into the kernel functions and avoid necessary global memory
    // writes. This requires moving the channel padding from python to C side.
    // if (cudaFuncSetAttribute(
    //         rasterize_to_pixels_3dgs_bwd_kernel<CDIM>,
    //         cudaFuncAttributeMaxDynamicSharedMemorySize,
    //         shmem_size
    //     ) != cudaSuccess) {
    //     AT_ERROR(
    //         "Failed to set maximum shared memory size (requested ",
    //         shmem_size,
    //         " bytes), try lowering tile_size."
    //     );
    // }

    rasterize_to_pixels_3dgs_bwd_kernel<CDIM>
        // <<<grid, threads, shmem_size, at::cuda::getCurrentCUDAStream()>>>(
        <<<grid, threads>>>(
            I,
            N,
            n_isects,
            packed,
            splats,
            colors.data_ptr<float>(),
            backgrounds.has_value() ? backgrounds.value().data_ptr<float>()
                                    : nullptr,
            masks.has_value() ? masks.value().data_ptr<bool>() : nullptr,
            image_width,
            image_height,
            tile_width,
            tile_height,
            tile_offsets.data_ptr<int32_t>(),
            flatten_ids.data_ptr<int32_t>(),
            render_Ts.data_ptr<float>(),
            last_ids.data_ptr<int32_t>(),
            v_render_colors.data_ptr<float>(),
            v_render_alphas.data_ptr<float>(),
            v_splats,
            v_colors.data_ptr<float>()
        );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}


std::tuple<
    Vanilla3DGS::Screen::TensorTuple,
    at::Tensor,  // v_colors
    std::optional<at::Tensor>  // absgrad
> rasterize_to_pixels_3dgs_bwd(
    // Gaussian parameters
    Vanilla3DGS::Screen::TensorTuple splats_tuple,
    const at::Tensor colors,                    // [..., N, channels] or [nnz, channels]
    const std::optional<at::Tensor> backgrounds, // [..., channels]
    const std::optional<at::Tensor> masks,       // [..., tile_height, tile_width]
    // image size
    const uint32_t image_width,
    const uint32_t image_height,
    const uint32_t tile_size,
    // intersections
    const at::Tensor tile_offsets, // [..., tile_height, tile_width]
    const at::Tensor flatten_ids,  // [n_isects]
    // forward outputs
    const at::Tensor render_Ts, // [..., image_height, image_width, 1]
    const at::Tensor last_ids,      // [..., image_height, image_width]
    // gradients of outputs
    const at::Tensor v_render_colors, // [..., image_height, image_width, channels]
    const at::Tensor v_render_alphas, // [..., image_height, image_width, 1]
    // options
    bool absgrad
) {
    DEVICE_GUARD(colors);
    CHECK_INPUT(colors);
    CHECK_INPUT(tile_offsets);
    CHECK_INPUT(flatten_ids);
    CHECK_INPUT(render_Ts);
    CHECK_INPUT(last_ids);
    CHECK_INPUT(v_render_colors);
    CHECK_INPUT(v_render_alphas);
    if (backgrounds.has_value())
        CHECK_INPUT(backgrounds.value());
    if (masks.has_value())
        CHECK_INPUT(masks.value());

    if (tile_size != TILE_SIZE)
        AT_ERROR("Unsupported tile size");

    uint32_t channels = colors.size(-1);

    Vanilla3DGS::Screen::Tensor splats(splats_tuple);
    Vanilla3DGS::Screen::Tensor v_splats = splats.zeros_like(absgrad);

    at::Tensor v_colors = torch::zeros_like(colors, splats.options());

#define __LAUNCH_KERNEL__(N)                                                   \
    case N:                                                                    \
        launch_rasterize_to_pixels_3dgs_bwd_kernel<N>(                         \
            splats,                                                            \
            colors,                                                            \
            backgrounds,                                                       \
            masks,                                                             \
            image_width,                                                       \
            image_height,                                                      \
            tile_offsets,                                                      \
            flatten_ids,                                                       \
            render_Ts,                                                         \
            last_ids,                                                          \
            v_render_colors,                                                   \
            v_render_alphas,                                                   \
            v_splats,                                                          \
            v_colors                                                           \
        );                                                                     \
        break;

    // TODO: an optimization can be done by passing the actual number of
    // channels into the kernel functions and avoid necessary global memory
    // writes. This requires moving the channel padding from python to C side.
    switch (channels) {
        __LAUNCH_KERNEL__(1)
        __LAUNCH_KERNEL__(2)
        __LAUNCH_KERNEL__(3)
        __LAUNCH_KERNEL__(4)
        __LAUNCH_KERNEL__(5)
        // __LAUNCH_KERNEL__(8)
        // __LAUNCH_KERNEL__(9)
        // __LAUNCH_KERNEL__(16)
        // __LAUNCH_KERNEL__(17)
        // __LAUNCH_KERNEL__(32)
        // __LAUNCH_KERNEL__(33)
        // __LAUNCH_KERNEL__(64)
        // __LAUNCH_KERNEL__(65)
        // __LAUNCH_KERNEL__(128)
        // __LAUNCH_KERNEL__(129)
        // __LAUNCH_KERNEL__(256)
        // __LAUNCH_KERNEL__(257)
        // __LAUNCH_KERNEL__(512)
        // __LAUNCH_KERNEL__(513)
    default:
        AT_ERROR("Unsupported number of channels: ", channels);
    }
#undef __LAUNCH_KERNEL__

    return std::make_tuple(
        v_splats.tuple(), v_colors, v_splats.absgrad
    );
}


// Explicit Instantiation: this should match how it is being called in .cpp
// file.
// TODO: this is slow to compile, can we do something about it?
#define __INS__(CDIM)                                                          \
    template void launch_rasterize_to_pixels_3dgs_bwd_kernel<CDIM>(            \
        Vanilla3DGS::Screen::Tensor splats,                                                \
        const at::Tensor colors,                                               \
        const std::optional<at::Tensor> backgrounds,                            \
        const std::optional<at::Tensor> masks,                                  \
        uint32_t image_width,                                                  \
        uint32_t image_height,                                                 \
        const at::Tensor tile_offsets,                                         \
        const at::Tensor flatten_ids,                                          \
        const at::Tensor render_Ts,                                            \
        const at::Tensor last_ids,                                             \
        const at::Tensor v_render_colors,                                      \
        const at::Tensor v_render_alphas,                                      \
        Vanilla3DGS::Screen::Tensor v_splats,                                              \
        at::Tensor v_opacities                                                 \
    );

__INS__(1)
__INS__(2)
__INS__(3)
__INS__(4)
__INS__(5)
// __INS__(8)
// __INS__(9)
// __INS__(16)
// __INS__(17)
// __INS__(32)
// __INS__(33)
// __INS__(64)
// __INS__(65)
// __INS__(128)
// __INS__(129)
// __INS__(256)
// __INS__(257)
// __INS__(512)
// __INS__(513)
#undef __INS__
