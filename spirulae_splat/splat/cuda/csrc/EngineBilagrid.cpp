// Engine bilagrid integration (RGB / depth / normal).
//
// Pool keys are scoped under "eng.bg.*" so the pool reports them as engine state.

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineInternal.h"
#include "EngineState.h"

#include "BilagridBindings.h"

#include <stdexcept>
#include <string>


static constexpr cudaStream_t kBilagridStream = (cudaStream_t)0;


void engine_init_bilagrid_rgb(int n_grids, std::string type, int L, int H, int W,
                              bool quantize_optim) {
    int C;
    if (type == "affine") C = 12;
    else if (type == "ppisp" || type == "loglinear") C = 9;
    else throw std::runtime_error("engine_init_bilagrid_rgb: unknown type " + type);
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_rgb: n_grids must be > 0");

    engine().bilagrid_rgb.type = type;
    engine().bilagrid_rgb.C = C;
    engine().bilagrid_rgb.quantize_optim = quantize_optim;
    engine().bilagrid_rgb.grids.resize("eng.bg.rgb.grids", n_grids, L, H, W, C);
    if (type == "affine") {
        engine().bilagrid_rgb.grids.zero();
        bilagrid_affine_identity_init(
            engine().bilagrid_rgb.grids.data_ptr(), n_grids, L, H, W, kBilagridStream);
    } else {
        engine().bilagrid_rgb.grids.zero();
    }
    engine().bilagrid_rgb.enabled = true;
    engine().bilagrid_rgb.optim_initialized = false;
}

void engine_init_bilagrid_depth(int n_grids, int L, int H, int W,
                                bool quantize_optim) {
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_depth: n_grids must be > 0");
    engine().bilagrid_depth.quantize_optim = quantize_optim;
    engine().bilagrid_depth.grids.resize("eng.bg.depth.grids", n_grids, L, H, W, 2);
    engine().bilagrid_depth.grids.zero();
    engine().bilagrid_depth.scalars.resize("eng.bg.depth.scalars", n_grids);
    engine().bilagrid_depth.scalars.zero();
    engine().bilagrid_depth.enabled = true;
    engine().bilagrid_depth.optim_initialized = false;
}

void engine_init_bilagrid_normal(int n_grids, int L, int H, int W,
                                 bool quantize_optim) {
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_normal: n_grids must be > 0");
    engine().bilagrid_normal.quantize_optim = quantize_optim;
    engine().bilagrid_normal.grids.resize("eng.bg.normal.grids", n_grids, L, H, W, 3);
    engine().bilagrid_normal.grids.zero();
    engine().bilagrid_normal.enabled = true;
    engine().bilagrid_normal.optim_initialized = false;
}

// Compute TV losses for each enabled bilagrid into a pool-backed 3-float
// device buffer at [rgb_tv, depth_tv, normal_tv]; unset entries stay 0.
void _engine_bilagrid_tv_into(float* tv_buf3_device) {
    cudaMemsetAsync(tv_buf3_device, 0, 3 * sizeof(float), kBilagridStream);
    auto run = [&](DeviceTensor5D<float>& grids, int slot) {
        if (grids.data_ptr() == nullptr) return;
        int N = (int)grids.size<0>();
        int L = (int)grids.size<1>(), H = (int)grids.size<2>(), W = (int)grids.size<3>();
        int C = (int)grids.size<4>();
        tv_loss_forward(grids.data_ptr(), tv_buf3_device + slot,
                        N, C, L, H, W, kBilagridStream);
    };
    if (engine().bilagrid_rgb.enabled) run(engine().bilagrid_rgb.grids, 0);
    if (engine().bilagrid_depth.enabled) run(engine().bilagrid_depth.grids, 1);
    if (engine().bilagrid_normal.enabled) run(engine().bilagrid_normal.grids, 2);
}

// Populate engine().bilagrid_cur_cam_indices from a host or device int32
// tensor. Empty/null tensor -> leave it null (kernels fall back to identity).
// Called by engine_bilagrid_forward / engine_ppisp_forward, and directly by
// engine_train_step before the forwards.
void _set_cur_cam_indices(TorchTensorView tv) {
    if (std::get<0>(tv) == 0 || std::get<2>(tv).empty()) {
        engine().bilagrid_cur_cam_indices = DeviceVector<int32_t>();
    } else {
        engine().bilagrid_cur_cam_indices = _hv_to_dv<int32_t>(
            "eng.bg.cam_indices", tv);
    }
}

// Apply bilagrid forward in-place for each enabled type for the current batch.
// Saves a pre-bilagrid copy used by the backward pass.
//
// The kernel uses indirect grid indexing: each image in the batch (index `ni`
// inside the kernel) reads from the grid slot `cam_indices[ni]`. This allows
// a multi-image batch with different per-image cam_idx values to be processed
// in a single kernel launch without any gather of grid params. When
// `cam_indices` is empty/null, kernels fall back to identity (ni -> ni).
void engine_bilagrid_forward(TorchTensorView cam_indices) {
    _set_cur_cam_indices(cam_indices);
    int H = engine().camera.height;
    int W = engine().camera.width;
    int C_batch = (int)engine().camera.num;
    if (C_batch <= 0) return;

    const int* cam_idx_dev = engine().bilagrid_cur_cam_indices.data_ptr();

    // --- RGB: applies to engine().fwd.renders.rgb ---
    if (engine().bilagrid_rgb.enabled) {
        if (std::get<0>(engine().fwd.renders).data_ptr() == nullptr)
            throw std::runtime_error("engine_bilagrid_forward: forward_3dgs must run first");

        // Save pre-bilagrid render and apply in place.
        engine().bilagrid_rgb.fwd_pre.resize("eng.bg.rgb.pre", C_batch, H, W);
        float3* render_rgb = std::get<0>(engine().fwd.renders).data_ptr();
        cudaMemcpyAsync(
            engine().bilagrid_rgb.fwd_pre.data_ptr(),
            render_rgb,
            (size_t)C_batch * H * W * sizeof(float3),
            cudaMemcpyDeviceToDevice, kBilagridStream);

        int L = (int)engine().bilagrid_rgb.grids.size<1>();
        int gH = (int)engine().bilagrid_rgb.grids.size<2>();
        int gW = (int)engine().bilagrid_rgb.grids.size<3>();
        // Pass the FULL grid table; kernel does the per-image indirect lookup.
        float* grid_ptr = engine().bilagrid_rgb.grids.data_ptr();

        // In-place: rgb input and output are the same buffer (the kernel reads
        // and writes index-by-index per pixel, so this is safe).
        if (engine().bilagrid_rgb.type == "affine") {
            bilagrid_uniform_sample_forward(
                grid_ptr, (const float*)render_rgb, (float*)render_rgb,
                /*N=*/C_batch, L, gH, gW, H, W,
                kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "ppisp") {
            bilagrid_ppisp_uniform_sample_forward(
                grid_ptr, (const float*)render_rgb, (float*)render_rgb,
                /*N=*/C_batch, L, gH, gW, H, W,
                kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "loglinear") {
            bilagrid_loglinear_uniform_sample_forward(
                grid_ptr, (const float*)render_rgb, (float*)render_rgb,
                /*N=*/C_batch, L, gH, gW, H, W,
                kBilagridStream, cam_idx_dev);
        }
    }

    // --- Depth (gt side) ---
    if (engine().bilagrid_depth.enabled && engine().gt.depth.data_ptr() != nullptr) {
        engine().bilagrid_depth.fwd_pre.resize("eng.bg.depth.pre", C_batch, H, W);
        float* gt_depth_ptr = engine().gt.depth.data_ptr();
        cudaMemcpyAsync(
            engine().bilagrid_depth.fwd_pre.data_ptr(),
            gt_depth_ptr,
            (size_t)C_batch * H * W * sizeof(float),
            cudaMemcpyDeviceToDevice, kBilagridStream);

        // Compute per-image scalar (median quantile) over pre-bilagrid depth
        // and scatter into bilagrid_depth_scalars at the cam_indices slots.
        // compute_depth_scalars_tensor writes scalars[batch_i] = quantile_i —
        // we want to land them at scalars[cam_indices[batch_i]] in the full
        // table. Simplest: compute into a temp [C_batch] buffer then scatter.
        float* scalar_full = engine().bilagrid_depth.scalars.data_ptr();
        float* tmp_scalars = DevicePool::global().acquire<float>(
            "eng.bg.depth.tmp_scalars", (size_t)C_batch);
        TorchTensorView depth_tv((uint64_t)engine().bilagrid_depth.fwd_pre.data_ptr(),
            4, {C_batch, H, W, 1, 1});
        TorchTensorView scalar_tv((uint64_t)tmp_scalars, 4, {C_batch});
        compute_depth_scalars_tensor(depth_tv, /*patched=*/true, scalar_tv);
        // Scatter tmp_scalars -> scalar_full at cam_indices.
        if (cam_idx_dev != nullptr) {
            bilagrid_scatter_floats(tmp_scalars, cam_idx_dev, C_batch,
                                    scalar_full, kBilagridStream);
        } else {
            // Identity: tmp_scalars[i] -> scalar_full[i].
            cudaMemcpyAsync(scalar_full, tmp_scalars,
                C_batch * sizeof(float), cudaMemcpyDeviceToDevice, kBilagridStream);
        }

        int L = (int)engine().bilagrid_depth.grids.size<1>();
        int gH = (int)engine().bilagrid_depth.grids.size<2>();
        int gW = (int)engine().bilagrid_depth.grids.size<3>();
        float* grid_ptr = engine().bilagrid_depth.grids.data_ptr();
        bilagrid_depth_uniform_sample_forward(
            grid_ptr,
            (const float*)engine().bilagrid_depth.fwd_pre.data_ptr(),
            scalar_full,  // kernel reads scalars[cam_indices[ni]]
            gt_depth_ptr,
            /*N=*/C_batch, L, gH, gW, H, W,
            kBilagridStream, cam_idx_dev);
    }

    // --- Normal (gt side) ---
    if (engine().bilagrid_normal.enabled && engine().gt.normal.data_ptr() != nullptr) {
        engine().bilagrid_normal.fwd_pre.resize("eng.bg.normal.pre", C_batch, H, W);
        float3* gt_normal_ptr = engine().gt.normal.data_ptr();
        cudaMemcpyAsync(
            engine().bilagrid_normal.fwd_pre.data_ptr(),
            gt_normal_ptr,
            (size_t)C_batch * H * W * sizeof(float3),
            cudaMemcpyDeviceToDevice, kBilagridStream);

        int L = (int)engine().bilagrid_normal.grids.size<1>();
        int gH = (int)engine().bilagrid_normal.grids.size<2>();
        int gW = (int)engine().bilagrid_normal.grids.size<3>();
        float* grid_ptr = engine().bilagrid_normal.grids.data_ptr();
        bilagrid_normal_uniform_sample_forward(
            grid_ptr,
            (const float*)engine().bilagrid_normal.fwd_pre.data_ptr(),
            (float*)gt_normal_ptr,
            /*N=*/C_batch, L, gH, gW, H, W,
            kBilagridStream, cam_idx_dev);
    }
}

// Bilagrid backward hook (Phase B): image-loss gradient now writes a
// SPARSE per-batch buffer (`eng.bg.*.image_grad` sized [C_batch, C, L, H, W])
// instead of the full N_grids-wide dense table. Reads of the bilagrid grid
// parameters still go through cam_indices for the true camera lookup. The
// per-batch buffer is zeroed up front so the kernel's atomicAdds start clean.
// The TV-loss contribution is no longer materialized here — it folds into the
// fused optimizer step (BilagridFusedAdam.cu).
//
// For RGB the kernel additionally overwrites v_render_rgb (post-bilagrid ->
// pre-bilagrid) so it can flow into rasterization backward.
static void _ensure_bilagrid_batch_grad(
    DeviceTensor5D<float>& grad,
    const DeviceTensor5D<float>& grids,
    int C_batch,
    const std::string& key
) {
    if (grids.data_ptr() == nullptr) return;
    int L = (int)grids.size<1>();
    int H = (int)grids.size<2>();
    int W = (int)grids.size<3>();
    int C = (int)grids.size<4>();
    grad.resize(key, C_batch, L, H, W, C);
    grad.zero();
}

void _engine_bilagrid_backward_hook(
    TorchTensorView v_render_rgb,   // [C, H, W, 3] — overwritten when RGB enabled
    TorchTensorView v_ref_depth,    // [C, H, W, 1] — input only
    TorchTensorView v_ref_normal    // [C, H, W, 3] — input only
) {
    int H = engine().camera.height;
    int W = engine().camera.width;
    int C_batch = (int)engine().camera.num;
    const int* cam_idx_dev = engine().bilagrid_cur_cam_indices.data_ptr();

    // --- RGB backward ---
    if (engine().bilagrid_rgb.enabled && std::get<0>(v_render_rgb) != 0 &&
        engine().bilagrid_rgb.fwd_pre.data_ptr() != nullptr)
    {
        int L = (int)engine().bilagrid_rgb.grids.size<1>();
        int gH = (int)engine().bilagrid_rgb.grids.size<2>();
        int gW = (int)engine().bilagrid_rgb.grids.size<3>();
        float* grid_ptr = engine().bilagrid_rgb.grids.data_ptr();
        _ensure_bilagrid_batch_grad(engine().bilagrid_rgb.image_grad,
                                    engine().bilagrid_rgb.grids, C_batch,
                                    "eng.bg.rgb.image_grad");
        float* grad_grid_ptr = engine().bilagrid_rgb.image_grad.data_ptr();

        // v_render_rgb is the gradient w.r.t. POST-bilagrid rgb (what entered loss).
        // The backward overwrites it with the gradient w.r.t. PRE-bilagrid rgb.
        float* v_rgb_ptr = (float*)std::get<0>(v_render_rgb);
        const float* rgb_pre_ptr = (const float*)engine().bilagrid_rgb.fwd_pre.data_ptr();

        if (engine().bilagrid_rgb.type == "affine") {
            bilagrid_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                /*block_x=*/8, /*block_y=*/8, /*target_tile_size=*/5,
                kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "ppisp") {
            bilagrid_ppisp_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                8, 8, 5, kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "loglinear") {
            bilagrid_loglinear_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                8, 8, 5, kBilagridStream, cam_idx_dev);
        }
    }

    // --- Depth backward (accumulate into bilagrid grad; discard input grad) ---
    // We pass v_depth=nullptr so the bilagrid backward skips the kernel that
    // would compute the ∂L/∂pre-bilagrid (= ∂L/∂gt_depth) buffer. That kernel's
    // output is never consumed downstream (GT isn't a parameter), and skipping
    // it eliminates a full-resolution C*H*W*4 byte scratch + one kernel launch.
    if (engine().bilagrid_depth.enabled && std::get<0>(v_ref_depth) != 0 &&
        engine().bilagrid_depth.fwd_pre.data_ptr() != nullptr)
    {
        int L = (int)engine().bilagrid_depth.grids.size<1>();
        int gH = (int)engine().bilagrid_depth.grids.size<2>();
        int gW = (int)engine().bilagrid_depth.grids.size<3>();
        float* grid_ptr = engine().bilagrid_depth.grids.data_ptr();
        _ensure_bilagrid_batch_grad(engine().bilagrid_depth.image_grad,
                                    engine().bilagrid_depth.grids, C_batch,
                                    "eng.bg.depth.image_grad");
        float* grad_grid_ptr = engine().bilagrid_depth.image_grad.data_ptr();
        float* scalars = engine().bilagrid_depth.scalars.data_ptr();
        bilagrid_depth_uniform_sample_backward_v1(
            grid_ptr,
            (const float*)engine().bilagrid_depth.fwd_pre.data_ptr(),
            scalars,
            (const float*)std::get<0>(v_ref_depth),
            grad_grid_ptr, /*v_depth=*/nullptr,
            /*N=*/C_batch, L, gH, gW, H, W,
            8, 8, 5, kBilagridStream, cam_idx_dev);
    }

    // --- Normal backward (accumulate; discard input grad) ---
    // Same reasoning as depth: v_rgb=nullptr skips the unused gt-side grad
    // kernel, saving a C*H*W*12 byte scratch + one kernel launch.
    if (engine().bilagrid_normal.enabled && std::get<0>(v_ref_normal) != 0 &&
        engine().bilagrid_normal.fwd_pre.data_ptr() != nullptr)
    {
        int L = (int)engine().bilagrid_normal.grids.size<1>();
        int gH = (int)engine().bilagrid_normal.grids.size<2>();
        int gW = (int)engine().bilagrid_normal.grids.size<3>();
        float* grid_ptr = engine().bilagrid_normal.grids.data_ptr();
        _ensure_bilagrid_batch_grad(engine().bilagrid_normal.image_grad,
                                    engine().bilagrid_normal.grids, C_batch,
                                    "eng.bg.normal.image_grad");
        float* grad_grid_ptr = engine().bilagrid_normal.image_grad.data_ptr();
        bilagrid_normal_uniform_sample_backward_v1(
            grid_ptr,
            (const float*)engine().bilagrid_normal.fwd_pre.data_ptr(),
            (const float*)std::get<0>(v_ref_normal),
            grad_grid_ptr, /*v_rgb=*/nullptr,
            /*N=*/C_batch, L, gH, gW, H, W,
            8, 8, 5, kBilagridStream, cam_idx_dev);
    }
}


// Ensure bilagrid Adam moment buffers are allocated and zeroed for each
// enabled bilagrid type. When quantize_optim is true, g1/g2 storage is uint8
// (1 byte/cell) plus a float4 per 256-cell block holding the {g1 min, g1 max,
// g2 min, g2 max} ranges — same scheme as the SH optimizer state.
//
// Note: the per-batch image_grad buffer is sized + zeroed per-iter in
// _ensure_bilagrid_batch_grad (called from _engine_bilagrid_backward_hook),
// since its size depends on the current batch size.
void _ensure_bilagrid_optim_state() {
    auto setup = [](DeviceTensor5D<float>& grids,
                    DeviceTensor5D<float>& g1,                 // used when !quantize
                    DeviceTensor5D<float>& g2,                 // used when !quantize
                    QuantizedAdamState<8, 256>& quant_state,   // used when quantize
                    bool quantize_optim,
                    const std::string& key_prefix,
                    bool& done) {
        if (done || grids.data_ptr() == nullptr) return;
        int64_t N = grids.size<0>();
        int64_t L = grids.size<1>(), H = grids.size<2>(), W = grids.size<3>();
        int64_t C = grids.size<4>();
        int64_t numel = N * L * H * W * C;
        if (quantize_optim) {
            constexpr int64_t BLOCK_SIZE = QuantizedAdamState<8, 256>::kBlockSize;
            int64_t n_blocks = (numel + BLOCK_SIZE - 1) / BLOCK_SIZE;
            quant_state.resize(key_prefix + ".bg_quant", numel, n_blocks);
            quant_state.zero();
        } else {
            g1.resize(key_prefix + ".g1", N, L, H, W, C);
            g1.zero();
            g2.resize(key_prefix + ".g2", N, L, H, W, C);
            g2.zero();
        }
        done = true;
    };
    if (engine().bilagrid_rgb.enabled) {
        setup(engine().bilagrid_rgb.grids,
              engine().bilagrid_rgb.g1, engine().bilagrid_rgb.g2,
              engine().bilagrid_rgb.quant_state,
              engine().bilagrid_rgb.quantize_optim,
              "eng.bg.rgb", engine().bilagrid_rgb.optim_initialized);
    }
    if (engine().bilagrid_depth.enabled) {
        setup(engine().bilagrid_depth.grids,
              engine().bilagrid_depth.g1, engine().bilagrid_depth.g2,
              engine().bilagrid_depth.quant_state,
              engine().bilagrid_depth.quantize_optim,
              "eng.bg.depth", engine().bilagrid_depth.optim_initialized);
    }
    if (engine().bilagrid_normal.enabled) {
        setup(engine().bilagrid_normal.grids,
              engine().bilagrid_normal.g1, engine().bilagrid_normal.g2,
              engine().bilagrid_normal.quant_state,
              engine().bilagrid_normal.quantize_optim,
              "eng.bg.normal", engine().bilagrid_normal.optim_initialized);
    }
}

void engine_bilagrid_optim_step(int step, const BilagridStepConfig& cfg) {
    _ensure_bilagrid_optim_state();
    int32_t adam_step = step + 1;
    int C_batch = (int)engine().camera.num;
    const int* cam_idx_dev = engine().bilagrid_cur_cam_indices.data_ptr();

    // Single fused pass per bilagrid type: reads the sparse per-batch
    // image_grad, computes the local TV gradient from `grids` neighbors, and
    // applies Adam (float or uint8-quantized) in one kernel. The previous
    // separate tv_loss_backward + zeroing of a dense grad table is gone —
    // TV is inlined, image_grad is the only persistent gradient state and
    // it gets re-zeroed at the start of the *next* iter's backward hook.
    auto run = [&](DeviceTensor5D<float>& grids,
                   DeviceTensor5D<float>& image_grad,
                   DeviceTensor5D<float>& g1,
                   DeviceTensor5D<float>& g2,
                   QuantizedAdamState<8, 256>& quant_state,
                   bool quantize_optim,
                   float lr, float tv_weight) {
        if (grids.data_ptr() == nullptr || lr <= 0.0f) return;
        if (image_grad.data_ptr() == nullptr) return;
        int N = (int)grids.size<0>();
        int L = (int)grids.size<1>(), H = (int)grids.size<2>(), W = (int)grids.size<3>();
        int C = (int)grids.size<4>();
        fused_bilagrid_tv_adam(
            grids.data_ptr(),
            quantize_optim ? nullptr : g1.data_ptr(),
            quantize_optim ? nullptr : g2.data_ptr(),
            quantize_optim ? quant_state.packed_ptr() : nullptr,
            quantize_optim ? quant_state.bounds_ptr() : nullptr,
            image_grad.data_ptr(),
            cam_idx_dev,
            N, C_batch, C, L, H, W,
            lr, tv_weight, adam_step,
            quantize_optim, kBilagridStream);
    };

    run(engine().bilagrid_rgb.grids, engine().bilagrid_rgb.image_grad,
        engine().bilagrid_rgb.g1, engine().bilagrid_rgb.g2,
        engine().bilagrid_rgb.quant_state, engine().bilagrid_rgb.quantize_optim,
        cfg.lr_rgb, cfg.tv_weight_rgb);
    run(engine().bilagrid_depth.grids, engine().bilagrid_depth.image_grad,
        engine().bilagrid_depth.g1, engine().bilagrid_depth.g2,
        engine().bilagrid_depth.quant_state, engine().bilagrid_depth.quantize_optim,
        cfg.lr_depth, cfg.tv_weight_depth);
    run(engine().bilagrid_normal.grids, engine().bilagrid_normal.image_grad,
        engine().bilagrid_normal.g1, engine().bilagrid_normal.g2,
        engine().bilagrid_normal.quant_state, engine().bilagrid_normal.quantize_optim,
        cfg.lr_normal, cfg.tv_weight_normal);
}
