// Engine ↔ DataManager glue.
//
// engine_setup_data_manager:
//   Build / replace the engine-owned DataManager from a parsed-out dataset.
//   Idempotent; an existing DataManager is destroyed first (and its disk
//   prefetch threads cleanly joined) before swap-in.
//
// engine_train_step_managed:
//   Pull the next training batch from engine().dm and dispatch the same
//   fused step as engine_train_step. The DataManager owns the per-step
//   host buffers; we hand pointers into those buffers to the engine as
//   TorchTensorViews. They stay alive until the *next* call to
//   next_train_batch().

#include "Camera.h"         // camera_model_to_string
#include "DataManager.h"
#include "Engine.h"
#include "EngineState.h"

#include <stdexcept>


void engine_setup_data_manager(
    DataManagerConfig         cfg,
    std::vector<int32_t>      camera_models,
    std::vector<std::string>  image_filenames,
    std::vector<std::string>  mask_filenames,
    std::vector<std::string>  depth_filenames,
    std::vector<std::string>  normal_filenames,
    std::vector<int32_t>      widths,
    std::vector<int32_t>      heights,
    std::vector<int32_t>      K_per_camera,
    std::vector<int32_t>      post_offsets,
    std::vector<float>        viewmats,
    std::vector<float>        intrins,
    std::vector<float>        dist_coeffs,
    std::vector<float>        input_intrins,
    std::vector<float>        input_dist_coeffs,
    std::vector<int32_t>      train_indices,
    std::vector<int32_t>      val_indices)
{
    // Drop the existing manager first so its scheduler / worker threads are
    // joined before we begin spinning up new ones.
    engine().dm.reset();

    engine().dm = std::make_unique<DataManager>(
        std::move(cfg), std::move(camera_models),
        std::move(image_filenames),  std::move(mask_filenames),
        std::move(depth_filenames),  std::move(normal_filenames),
        std::move(widths),           std::move(heights),
        std::move(K_per_camera),     std::move(post_offsets),
        std::move(viewmats),         std::move(intrins),
        std::move(dist_coeffs),
        std::move(input_intrins),    std::move(input_dist_coeffs),
        std::move(train_indices),    std::move(val_indices));
}


std::map<std::string, float> engine_train_step_managed(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    const EngineStepConfig& cfg)
{
    if (!engine().dm) {
        throw std::runtime_error(
            "engine_train_step_managed: DataManager not configured — "
            "call engine_setup_data_manager(...) first.");
    }

    const DecodedBatch& b = engine().dm->next_train_batch();

    // Build the per-step bilagrid cam-indices buffer at POST-split size.
    // For non-warp batches (K == 1) this is just `b.indices`. For warp
    // batches it expands to {post_offsets[j] + k for j in 0..B_in for k}.
    // Static buffer so we don't allocate per step.
    static std::vector<int32_t> _bg_idx_buf;
    _bg_idx_buf.resize((size_t)b.num);
    // bilagrid_cam_indices must be the POST-split camera id, not the input
    // dataset id. For mixed datasets (e.g. K=5 fisheye + K=1 pinhole) those
    // diverge after the first K>1 input, so reading b.indices[j] for K==1
    // batches points the thumbnail kernel + bilagrid + PPISP at fisheye
    // slots and leaves the real pinhole slots blank.
    // `b.post_offsets[j] = _post_offsets[b.indices[j]]` is the starting
    // post-split slot for input image j, already populated by
    // DataManager::allocate_batch.
    {
        int K = b.K;
        for (int j = 0; j < b.input_num; ++j) {
            int32_t off = b.post_offsets[j];
            for (int k = 0; k < K; ++k) _bg_idx_buf[(size_t)j * K + k] = off + k;
        }
    }
    TorchTensorView bilagrid_cam_indices(
        (uint64_t)_bg_idx_buf.data(), 4,
        {(int64_t)b.num, 1LL});

    if (b.K <= 1) {
        // Standard path: GT is already at engine shape, no warp needed.
        return engine_train_step(
            step, max_steps,
            std::move(primitive), sh_degree, packed,
            (int)b.width, (int)b.height,
            camera_model_to_string(b.model),
            b.viewmats_view, b.intrins_view, b.dist_coeffs_view,
            b.rgb_view, b.depth_view, b.normal_view, b.mask_view,
            bilagrid_cam_indices,
            cfg);
    }

    // ---- Warp path: fisheye / equisolid + warp_to_pinhole OR equirectangular
    if (std::get<0>(b.depth_view) != 0 || std::get<0>(b.normal_view) != 0) {
        throw std::runtime_error(
            "engine_train_step_managed: depth/normal GT not supported when "
            "warp_to_pinhole splitting is active.");
    }
    return engine_train_step_warped(
        step, max_steps,
        std::move(primitive), sh_degree, packed,
        (int)b.width, (int)b.height,
        b.viewmats_view, b.intrins_view, b.dist_coeffs_view,
        camera_model_to_string(b.input_model),
        b.input_num, b.input_height, b.input_width, b.K,
        b.input_intrins_view, b.input_dist_coeffs_view,
        b.rgb_view, b.mask_view,
        b.mask_height, b.mask_width,
        (uint64_t)b.axes_dev,
        bilagrid_cam_indices,
        cfg);
}
