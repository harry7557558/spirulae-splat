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
    CameraModelType           model,
    std::vector<std::string>  image_filenames,
    std::vector<std::string>  mask_filenames,
    std::vector<std::string>  depth_filenames,
    std::vector<std::string>  normal_filenames,
    std::vector<int32_t>      widths,
    std::vector<int32_t>      heights,
    std::vector<float>        viewmats,
    std::vector<float>        intrins,
    std::vector<float>        dist_coeffs,
    std::vector<int32_t>      train_indices,
    std::vector<int32_t>      val_indices)
{
    // Drop the existing manager first so its scheduler / worker threads are
    // joined before we begin spinning up new ones.
    engine().dm.reset();

    engine().dm = std::make_unique<DataManager>(
        std::move(cfg), model,
        std::move(image_filenames),  std::move(mask_filenames),
        std::move(depth_filenames),  std::move(normal_filenames),
        std::move(widths),           std::move(heights),
        std::move(viewmats),         std::move(intrins),
        std::move(dist_coeffs),
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

    // Bilagrid / PPISP cam selector — reuse the batch's dataset indices.
    // The DataManager guarantees indices is non-null and length == b.num.
    TorchTensorView bilagrid_cam_indices(
        (uint64_t)b.indices.data(), 4,
        {(int64_t)b.num, 1LL});

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
