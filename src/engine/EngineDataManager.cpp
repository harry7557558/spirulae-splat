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

#include "core/Camera.h"         // camera_model_to_string
#include "data/DataManager.h"
#include "engine/Engine.h"
#include "engine/EngineState.h"

#include <cstdio>
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


// Resolve the `split_batch` + `use_fused_proj_bwd_optim` conflict. Both
// turned on at once is an inconsistency: split_batch loops over single
// cameras and accumulates atomicAdd into world-grad buffers, whereas FPBO
// folds projection-bwd into the optim kernel and never materializes those
// buffers. If we know (from the DataManager) that the post-split batch
// size will always be 1, split_batch is a no-op and FPBO is the correct
// choice. Otherwise the user wants the memory win of split_batch and we
// must turn FPBO off. Prints a one-shot warning describing the choice.
static EngineStepConfig _resolve_split_vs_fpbo(const EngineStepConfig& cfg,
                                               int64_t max_batch_known) {
    if (!cfg.optim.split_batch || !cfg.optim.use_fused_proj_bwd_optim)
        return cfg;
    EngineStepConfig out = cfg;
    static bool warned = false;
    // max_batch_known == -1 -> caller couldn't determine ahead of time;
    // treat as "may exceed 1" so the user gets the memory-safe path.
    bool batch_le_one = (max_batch_known > 0 && max_batch_known <= 1);
    if (batch_le_one) {
        out.optim.split_batch = false;
        if (!warned) {
            fprintf(stderr,
                "[spirulae_splat] WARNING: both `split_batch` and "
                "`use_fused_proj_bwd_optim` are enabled, but the dataset's max "
                "post-split batch size is 1; split_batch would be a no-op. "
                "Disabling split_batch and keeping FPBO.\n");
            warned = true;
        }
    } else {
        out.optim.use_fused_proj_bwd_optim = false;
        if (!warned) {
            fprintf(stderr,
                "[spirulae_splat] WARNING: both `split_batch` and "
                "`use_fused_proj_bwd_optim` are enabled, but the post-split "
                "batch size can exceed 1 (max=%lld); FPBO is incompatible with "
                "sub-batched gradient accumulation. Disabling FPBO and keeping "
                "split_batch.\n",
                (long long)max_batch_known);
            warned = true;
        }
    }
    return out;
}


std::map<std::string, float> engine_train_step_managed(
    int step, int max_steps,
    std::string primitive,
    int sh_degree,
    bool packed,
    const EngineStepConfig& cfg_in)
{
    if (!engine().dm) {
        throw std::runtime_error(
            "engine_train_step_managed: DataManager not configured — "
            "call engine_setup_data_manager(...) first.");
    }

    // Resolve split_batch vs FPBO using the DataManager's view of the
    // dataset (max train batch size × max K). Done once per session via
    // the static warned flag inside the helper.
    const EngineStepConfig cfg = _resolve_split_vs_fpbo(
        cfg_in, engine().dm->max_input_batch_size());

    const TrainStep& stp = engine().dm->next_train_step();
    if (stp.subs.empty())
        throw std::runtime_error("engine_train_step_managed: empty training step");

    // Build a POST-split bilagrid cam-index buffer for one sub-batch.
    // bilagrid_cam_indices must be the POST-split camera id, not the input
    // dataset id: for mixed datasets (e.g. K=5 fisheye + K=1 pinhole) those
    // diverge after the first K>1 input, so reading b.indices[j] for K==1
    // batches would point the thumbnail kernel + bilagrid + PPISP at fisheye
    // slots. `b.post_offsets[j]` is the starting post-split slot for input
    // image j (populated by DataManager::allocate_batch).
    auto build_bg_idx = [](const DecodedBatch& b, std::vector<int32_t>& buf) {
        buf.resize((size_t)b.num);
        int K = b.K;
        for (int j = 0; j < b.input_num; ++j) {
            int32_t off = b.post_offsets[j];
            for (int k = 0; k < K; ++k) buf[(size_t)j * K + k] = off + k;
        }
    };

    // ---- Fast path: single homogeneous sub-batch (unchanged behavior) -----
    if (stp.subs.size() == 1) {
        const DecodedBatch& b = *stp.subs[0];
        static std::vector<int32_t> _bg_idx_buf;   // reused across steps
        build_bg_idx(b, _bg_idx_buf);
        TorchTensorView bilagrid_cam_indices(
            (uint64_t)_bg_idx_buf.data(), 4, {(int64_t)b.num, 1LL});

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

        // ---- Warp path: fisheye / equisolid + warp_to_pinhole OR equirect --
        // Depth is warped to per-face ray depth, normal is rotated into each
        // face's camera frame (see set_training_data_warped).
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
            b.depth_view, b.depth_height, b.depth_width,
            b.normal_view, b.normal_height, b.normal_width,
            (uint64_t)b.axes_dev,
            bilagrid_cam_indices,
            cfg);
    }

    // ---- Heterogeneous step: multiple non-warp sub-batches, accumulated ----
    // The DataManager only packs K == 1 (non-warp) groups across resolutions,
    // so every sub-batch here is a plain pinhole/model batch that goes through
    // the standard (non-warp) fwd/bwd. Per-sub bilagrid index buffers must
    // stay alive across the whole engine_train_step_hetero call.
    std::vector<std::vector<int32_t>> bg_bufs(stp.subs.size());
    std::vector<HeteroSubBatch> hsubs;
    hsubs.reserve(stp.subs.size());
    for (size_t i = 0; i < stp.subs.size(); ++i) {
        const DecodedBatch& b = *stp.subs[i];
        if (b.K > 1) {
            throw std::runtime_error(
                "engine_train_step_managed: warp (K>1) sub-batch inside a "
                "multi-sub-batch step — the scheduler must never pack these.");
        }
        build_bg_idx(b, bg_bufs[i]);

        HeteroSubBatch hs;
        hs.width        = (int)b.width;
        hs.height       = (int)b.height;
        hs.num          = (int)b.num;
        hs.camera_model = camera_model_to_string(b.model);
        hs.viewmats     = b.viewmats_view;
        hs.intrins      = b.intrins_view;
        hs.dist_coeffs  = b.dist_coeffs_view;
        hs.gt_rgb       = b.rgb_view;
        hs.gt_depth     = b.depth_view;
        hs.gt_normal    = b.normal_view;
        hs.gt_alpha     = b.mask_view;
        hs.bilagrid_cam_indices = TorchTensorView(
            (uint64_t)bg_bufs[i].data(), 4, {(int64_t)b.num, 1LL});
        hsubs.push_back(std::move(hs));
    }

    return engine_train_step_hetero(
        step, max_steps,
        std::move(primitive), sh_degree, packed,
        hsubs, cfg);
}


// ---------------------------------------------------------------------------
// Eval: pull a batch, install it, forward only. See Engine.h.
// ---------------------------------------------------------------------------
static TorchTensorView _tv_null() { return {0, 0, {}}; }

int engine_eval_forward(std::string primitive, int sh_degree, bool packed) {
    if (!engine().dm)
        throw std::runtime_error(
            "engine_eval_forward: DataManager not configured — call "
            "engine_setup_data_manager(...) with the eval split first.");

    const TrainStep& stp = engine().dm->next_train_step();
    if (stp.subs.empty()) return 0;
    if (stp.subs.size() != 1)
        throw std::runtime_error(
            "engine_eval_forward: expected one sub-batch per eval step "
            "(set train_batch_size = 1)");

    const DecodedBatch& b = *stp.subs[0];

    if (b.K <= 1) {
        set_camera_params((int)b.width, (int)b.height,
                          camera_model_to_string(b.model),
                          b.viewmats_view, b.intrins_view, b.dist_coeffs_view);
        // No depth/normal: eval compares RGB only, and skipping them avoids
        // the linear->ray depth conversion pass.
        set_training_data(b.rgb_view, _tv_null(), _tv_null(),
                          b.mask_view, true);
    } else {
        set_camera_params((int)b.width, (int)b.height, "PINHOLE",
                          b.viewmats_view, b.intrins_view, b.dist_coeffs_view);
        set_training_data_warped(
            camera_model_to_string(b.input_model),
            b.input_num, (int)b.input_height, (int)b.input_width,
            b.K, (int)b.height, (int)b.width,
            b.rgb_view, b.mask_view, (int)b.mask_height, (int)b.mask_width,
            _tv_null(), 0, 0,
            _tv_null(), 0, 0,
            true,
            b.input_intrins_view, b.input_dist_coeffs_view,
            (uint64_t)b.axes_dev);
    }

    forward_3dgs(std::move(primitive), sh_degree, packed);
    return (int)b.num;
}
