// Engine bilagrid integration (RGB / depth / normal).
//
// Pool keys are scoped under "eng.bg.*" so the pool reports them as engine state.

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineInternal.h"
#include "EngineState.h"

#include "BilagridBindings.h"
#include "BilagridReader.cuh"

#include <stdexcept>
#include <string>


static constexpr cudaStream_t kBilagridStream = (cudaStream_t)0;


// Validate optim_bits: 32 = off, 4 or 8 = enabled.
static void _validate_optim_bits(const char* fn, int bits) {
    if (bits != 32 && bits != 8 && bits != 4) {
        throw std::runtime_error(
            std::string(fn) + ": optim_bits must be 32 (off), 4, or 8; got " +
            std::to_string(bits));
    }
}
// Validate value_bits: 32 = off, 16 = enabled.
static void _validate_value_bits(const char* fn, int bits) {
    if (bits != 32 && bits != 16) {
        throw std::runtime_error(
            std::string(fn) + ": value_bits must be 32 (off) or 16; got " +
            std::to_string(bits));
    }
}

// Allocate the 16-bit packed grid + per-256-cell float2 bounds for a bilagrid
// when value_bits == 16. The bound buffer is zeroed (mm = (0, 0) per block);
// the first Adam step's value-bound block-reduce installs real bounds.
template<typename BG>
static void _alloc_grids_quant(BG& bg, const char* key_prefix, int64_t cells) {
    constexpr int64_t BLOCK_SIZE = QuantizedTensor<16, 256>::kBlockSize;
    int64_t n_bounds = (cells + BLOCK_SIZE - 1) / BLOCK_SIZE;
    bg.grids_quant.resize(std::string(key_prefix) + ".q", cells, n_bounds);
    bg.grids_quant.zero();
}

// One-shot encode-then-free: encode the fp32 grids into the packed buffer +
// bounds, then drop the fp32 pool slot. Subsequent sampler / Adam calls read
// and write through `grids_quant` only. The DT5D's shape descriptor is kept
// (set_shape_no_alloc) so downstream size<>() queries still work.
template<typename BG>
static void _encode_and_free_fp32(BG& bg, const char* fp32_pool_key) {
    int64_t N = bg.grids.template size<0>();
    int64_t L = bg.grids.template size<1>();
    int64_t H = bg.grids.template size<2>();
    int64_t W = bg.grids.template size<3>();
    int64_t C = bg.grids.template size<4>();
    int64_t total_cells = N * L * H * W * C;
    bilagrid_encode_q16_launch(
        bg.grids.data_ptr(),
        (uint16_t*)bg.grids_quant.packed_ptr(),
        bg.grids_quant.bounds_ptr(),
        total_cells,
        (cudaStream_t)0);
    DevicePool::global().free(fp32_pool_key);
    bg.grids.set_shape_no_alloc(N, L, H, W, C);
}

// Construct a BilagridReader that dispatches on the bilagrid's value_bits.
template<typename BG>
static inline BilagridReader _bg_reader(const BG& bg) {
    if (bg.value_bits == 32) return BilagridReader(bg.grids.data_ptr());
    return BilagridReader::from_q16(
        (const uint16_t*)bg.grids_quant.packed_ptr(),
        bg.grids_quant.bounds_ptr());
}

void engine_init_bilagrid_rgb(int n_grids, std::string type, int L, int H, int W,
                              int optim_bits, int value_bits, bool use_adagrad) {
    _validate_optim_bits("engine_init_bilagrid_rgb", optim_bits);
    _validate_value_bits("engine_init_bilagrid_rgb", value_bits);
    int C;
    if (type == "affine") C = 12;
    else if (type == "ppisp" || type == "loglinear") C = 9;
    else throw std::runtime_error("engine_init_bilagrid_rgb: unknown type " + type);
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_rgb: n_grids must be > 0");

    engine().bilagrid_rgb.type = type;
    engine().bilagrid_rgb.C = C;
    engine().bilagrid_rgb.optim_bits = optim_bits;
    engine().bilagrid_rgb.value_bits = value_bits;
    engine().bilagrid_rgb.use_adagrad = use_adagrad;
    // fp32 grids is the working buffer regardless of value_bits during init
    // (identity-init and zero-init both need fp32 writes); when quantized we
    // encode it to packed below and then free the fp32 slot at first
    // optimizer step (handled in EngineLoss / EngineBilagrid bilagrid step).
    engine().bilagrid_rgb.grids.resize("eng.bg.rgb.grids", n_grids, L, H, W, C);
    if (type == "affine") {
        engine().bilagrid_rgb.grids.zero();
        bilagrid_affine_identity_init(
            engine().bilagrid_rgb.grids.data_ptr(), n_grids, L, H, W, kBilagridStream);
    } else {
        engine().bilagrid_rgb.grids.zero();
    }
    if (value_bits == 16) {
        int64_t cells = (int64_t)n_grids * L * H * W * C;
        _alloc_grids_quant(engine().bilagrid_rgb, "eng.bg.rgb.grids_q", cells);
        _encode_and_free_fp32(engine().bilagrid_rgb, "eng.bg.rgb.grids");
    } else {
        engine().bilagrid_rgb.grids_quant = QuantizedTensor<16, 256>();
    }
    engine().bilagrid_rgb.enabled = true;
    engine().bilagrid_rgb.optim_initialized = false;
}

void engine_init_bilagrid_depth(int n_grids, int L, int H, int W,
                                int optim_bits, int value_bits, bool use_adagrad) {
    _validate_optim_bits("engine_init_bilagrid_depth", optim_bits);
    _validate_value_bits("engine_init_bilagrid_depth", value_bits);
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_depth: n_grids must be > 0");
    engine().bilagrid_depth.optim_bits = optim_bits;
    engine().bilagrid_depth.value_bits = value_bits;
    engine().bilagrid_depth.use_adagrad = use_adagrad;
    engine().bilagrid_depth.grids.resize("eng.bg.depth.grids", n_grids, L, H, W, 2);
    engine().bilagrid_depth.grids.zero();
    engine().bilagrid_depth.scalars.resize("eng.bg.depth.scalars", n_grids);
    engine().bilagrid_depth.scalars.zero();
    if (value_bits == 16) {
        int64_t cells = (int64_t)n_grids * L * H * W * 2;
        _alloc_grids_quant(engine().bilagrid_depth, "eng.bg.depth.grids_q", cells);
        _encode_and_free_fp32(engine().bilagrid_depth, "eng.bg.depth.grids");
    } else {
        engine().bilagrid_depth.grids_quant = QuantizedTensor<16, 256>();
    }
    engine().bilagrid_depth.enabled = true;
    engine().bilagrid_depth.optim_initialized = false;
}

void engine_init_bilagrid_normal(int n_grids, int L, int H, int W,
                                 int optim_bits, int value_bits, bool use_adagrad) {
    _validate_optim_bits("engine_init_bilagrid_normal", optim_bits);
    _validate_value_bits("engine_init_bilagrid_normal", value_bits);
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_bilagrid_normal: n_grids must be > 0");
    engine().bilagrid_normal.optim_bits = optim_bits;
    engine().bilagrid_normal.value_bits = value_bits;
    engine().bilagrid_normal.use_adagrad = use_adagrad;
    engine().bilagrid_normal.grids.resize("eng.bg.normal.grids", n_grids, L, H, W, 3);
    engine().bilagrid_normal.grids.zero();
    if (value_bits == 16) {
        int64_t cells = (int64_t)n_grids * L * H * W * 3;
        _alloc_grids_quant(engine().bilagrid_normal, "eng.bg.normal.grids_q", cells);
        _encode_and_free_fp32(engine().bilagrid_normal, "eng.bg.normal.grids");
    } else {
        engine().bilagrid_normal.grids_quant = QuantizedTensor<16, 256>();
    }
    engine().bilagrid_normal.enabled = true;
    engine().bilagrid_normal.optim_initialized = false;
}

// Compute TV losses for each enabled bilagrid into a pool-backed 3-float
// device buffer at [rgb_tv, depth_tv, normal_tv]; unset entries stay 0.
void _engine_bilagrid_tv_into(float* tv_buf3_device) {
    cudaMemsetAsync(tv_buf3_device, 0, 3 * sizeof(float), kBilagridStream);
    auto run = [&](auto& bg, int slot) {
        if (bg.grids.template size<0>() == 0) return;
        int N = (int)bg.grids.template size<0>();
        int L = (int)bg.grids.template size<1>();
        int H = (int)bg.grids.template size<2>();
        int W = (int)bg.grids.template size<3>();
        int C = (int)bg.grids.template size<4>();
        BilagridReader br = _bg_reader(bg);
        tv_loss_forward(br, tv_buf3_device + slot,
                        N, C, L, H, W, kBilagridStream);
    };
    if (engine().bilagrid_rgb.enabled) run(engine().bilagrid_rgb, 0);
    if (engine().bilagrid_depth.enabled) run(engine().bilagrid_depth, 1);
    if (engine().bilagrid_normal.enabled) run(engine().bilagrid_normal, 2);
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
        auto& fwd_rgb_tensor = std::get<0>(engine().fwd.renders);
        if (fwd_rgb_tensor.data_ptr() == nullptr)
            throw std::runtime_error("engine_bilagrid_forward: forward_3dgs must run first");

        // pre = current renders.rgb (pointer alias, no D2D copy). The bilagrid
        // kernels already accept distinct in/out pointers, so we read pre and
        // write into a fresh "post" buffer, then re-point engine state.
        engine().bilagrid_rgb.fwd_pre = fwd_rgb_tensor;
        DeviceTensor3D<float3> post_rgb;
        post_rgb.resize("eng.bg.rgb.post", C_batch, H, W);
        const float* in_ptr  = (const float*)engine().bilagrid_rgb.fwd_pre.data_ptr();
        float*       out_ptr = (float*)post_rgb.data_ptr();

        int L = (int)engine().bilagrid_rgb.grids.size<1>();
        int gH = (int)engine().bilagrid_rgb.grids.size<2>();
        int gW = (int)engine().bilagrid_rgb.grids.size<3>();
        // Pass the FULL grid table; kernel does the per-image indirect lookup.
        BilagridReader grid_ptr = _bg_reader(engine().bilagrid_rgb);

        if (engine().bilagrid_rgb.type == "affine") {
            bilagrid_uniform_sample_forward(
                grid_ptr, in_ptr, out_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "ppisp") {
            bilagrid_ppisp_uniform_sample_forward(
                grid_ptr, in_ptr, out_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "loglinear") {
            bilagrid_loglinear_uniform_sample_forward(
                grid_ptr, in_ptr, out_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                kBilagridStream, cam_idx_dev);
        }
        fwd_rgb_tensor = post_rgb;
    }

    // --- Depth (gt side) ---
    if (engine().bilagrid_depth.enabled && engine().gt.depth.data_ptr() != nullptr) {
        // pre = current gt.depth (pointer alias). The bilagrid sample kernel
        // reads `pre` and writes a fresh "post" buffer; we re-point gt.depth
        // at the post so the loss reads the transformed depth and the bwd
        // hook still has the original gt via fwd_pre. The bilagrid acts on
        // the GT depth at GT's own (H, W), which may differ from render shape.
        engine().bilagrid_depth.fwd_pre = engine().gt.depth;
        int H_d = (int)engine().gt.depth.template size<1>();
        int W_d = (int)engine().gt.depth.template size<2>();
        DeviceTensor3D<float> post_depth;
        post_depth.resize("eng.bg.depth.post", C_batch, H_d, W_d);

        // Compute per-image scalar (median quantile) over pre-bilagrid depth
        // (the un-transformed gt.depth) and scatter into the per-camera
        // scalars table at the cam_indices slots.
        float* scalar_full = engine().bilagrid_depth.scalars.data_ptr();
        float* tmp_scalars = DevicePool::global().acquire<float>(
            "eng.bg.depth.tmp_scalars", (size_t)C_batch);
        TorchTensorView depth_tv((uint64_t)engine().bilagrid_depth.fwd_pre.data_ptr(),
            4, {C_batch, H_d, W_d, 1, 1});
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
        BilagridReader grid_ptr = _bg_reader(engine().bilagrid_depth);
        bilagrid_depth_uniform_sample_forward(
            grid_ptr,
            (const float*)engine().bilagrid_depth.fwd_pre.data_ptr(),
            scalar_full,  // kernel reads scalars[cam_indices[ni]]
            post_depth.data_ptr(),
            /*N=*/C_batch, L, gH, gW, H_d, W_d,
            kBilagridStream, cam_idx_dev);
        engine().gt.depth = post_depth;
    }

    // --- Normal (gt side) ---
    if (engine().bilagrid_normal.enabled && engine().gt.normal.data_ptr() != nullptr) {
        // Same pointer-swap trick as depth: alias pre, write post, re-point gt.
        engine().bilagrid_normal.fwd_pre = engine().gt.normal;
        int H_n = (int)engine().gt.normal.template size<1>();
        int W_n = (int)engine().gt.normal.template size<2>();
        DeviceTensor3D<float3> post_normal;
        post_normal.resize("eng.bg.normal.post", C_batch, H_n, W_n);

        int L = (int)engine().bilagrid_normal.grids.size<1>();
        int gH = (int)engine().bilagrid_normal.grids.size<2>();
        int gW = (int)engine().bilagrid_normal.grids.size<3>();
        BilagridReader grid_ptr = _bg_reader(engine().bilagrid_normal);
        bilagrid_normal_uniform_sample_forward(
            grid_ptr,
            (const float*)engine().bilagrid_normal.fwd_pre.data_ptr(),
            (float*)post_normal.data_ptr(),
            /*N=*/C_batch, L, gH, gW, H_n, W_n,
            kBilagridStream, cam_idx_dev);
        engine().gt.normal = post_normal;
    }
}

// Bilagrid backward hook (Phase B): image-loss gradient now writes a
// SPARSE per-batch buffer (`eng.bg.*.image_grad` sized [C_batch, C, L, H, W])
// instead of the full N_grids-wide dense table. Reads of the bilagrid grid
// parameters still go through cam_indices for the true camera lookup. The
// per-batch buffer is zeroed up front so the kernel's atomicAdds start clean.
// The TV-loss contribution is no longer materialized here -- it folds into the
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
    // grids may have data_ptr == nullptr when value_bits == 16 (canonical
    // store is the packed buffer); shape descriptor is still alive. Gate on
    // the SHAPE rather than data_ptr so the gradient buffer is still sized.
    if (grids.size<0>() == 0) return;
    int L = (int)grids.size<1>();
    int H = (int)grids.size<2>();
    int W = (int)grids.size<3>();
    int C = (int)grids.size<4>();
    grad.resize(key, C_batch, L, H, W, C);
    grad.zero();
}

void _engine_bilagrid_backward_hook(
    TorchTensorView v_render_rgb,   // [C, H, W, 3] -- overwritten when RGB enabled
    TorchTensorView v_ref_depth,    // [C, H, W, 1] -- input only
    TorchTensorView v_ref_normal    // [C, H, W, 3] -- input only
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
        BilagridReader grid_ptr = _bg_reader(engine().bilagrid_rgb);
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
                /*target_tile_size=*/5,
                kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "ppisp") {
            bilagrid_ppisp_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                5, kBilagridStream, cam_idx_dev);
        } else if (engine().bilagrid_rgb.type == "loglinear") {
            bilagrid_loglinear_uniform_sample_backward_v1(
                grid_ptr, rgb_pre_ptr, v_rgb_ptr,
                grad_grid_ptr, v_rgb_ptr,
                /*N=*/C_batch, L, gH, gW, H, W,
                5, kBilagridStream, cam_idx_dev);
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
        // GT depth may have its own (H, W) distinct from render shape. The
        // bilagrid kernel samples at the GT buffer's resolution; v_ref_depth
        // has already been scatter-summed by the per-pixel loss kernel at
        // exactly that resolution, so the two shapes match by construction.
        int H_d = (int)engine().gt.depth.template size<1>();
        int W_d = (int)engine().gt.depth.template size<2>();
        BilagridReader grid_ptr = _bg_reader(engine().bilagrid_depth);
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
            /*N=*/C_batch, L, gH, gW, H_d, W_d,
            5, kBilagridStream, cam_idx_dev);
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
        int H_n = (int)engine().gt.normal.template size<1>();
        int W_n = (int)engine().gt.normal.template size<2>();
        BilagridReader grid_ptr = _bg_reader(engine().bilagrid_normal);
        _ensure_bilagrid_batch_grad(engine().bilagrid_normal.image_grad,
                                    engine().bilagrid_normal.grids, C_batch,
                                    "eng.bg.normal.image_grad");
        float* grad_grid_ptr = engine().bilagrid_normal.image_grad.data_ptr();
        bilagrid_normal_uniform_sample_backward_v1(
            grid_ptr,
            (const float*)engine().bilagrid_normal.fwd_pre.data_ptr(),
            (const float*)std::get<0>(v_ref_normal),
            grad_grid_ptr, /*v_rgb=*/nullptr,
            /*N=*/C_batch, L, gH, gW, H_n, W_n,
            5, kBilagridStream, cam_idx_dev);
    }
}


// Ensure bilagrid Adam moment buffers are allocated and zeroed for each
// enabled bilagrid type. When quantize_optim is true, g1/g2 storage is uint8
// (1 byte/cell) plus a float4 per 256-cell block holding the {g1 min, g1 max,
// g2 min, g2 max} ranges -- same scheme as the SH optimizer state.
//
// Note: the per-batch image_grad buffer is sized + zeroed per-iter in
// _ensure_bilagrid_batch_grad (called from _engine_bilagrid_backward_hook),
// since its size depends on the current batch size.
void _ensure_bilagrid_optim_state() {
    auto setup = [](DeviceTensor5D<float>& grids,
                    DeviceTensor5D<float>& g1,                 // Adam, !quantize
                    DeviceTensor5D<float>& g2,                 // Adam, !quantize
                    QuantizedAdamState<8, 256>& quant_state,   // Adam, quantize
                    DeviceTensor5D<float>& accum_f,            // AdaGrad, !quantize
                    QuantizedTensorLog<8, 256>& adagrad_quant, // AdaGrad, quantize
                    int  optim_bits,                            // 32 = off, 4 or 8
                    bool use_adagrad,
                    const std::string& key_prefix,
                    bool& done) {
        bool quantize_optim = (optim_bits != 32);
        // Gate on the SHAPE descriptor, not data_ptr -- when value_bits=16 the
        // fp32 grids buffer is freed (canonical store is grids_quant) but the
        // optim state is sized off the same dims, so the allocation must
        // still run.
        if (done || grids.size<0>() == 0) return;
        int64_t N = grids.size<0>();
        int64_t L = grids.size<1>(), H = grids.size<2>(), W = grids.size<3>();
        int64_t C = grids.size<4>();
        int64_t numel = N * L * H * W * C;
        if (use_adagrad) {
            if (quantize_optim) {
                constexpr int64_t BLOCK_SIZE =
                    QuantizedTensorLog<8, 256>::kBlockSize;
                int64_t n_blocks = (numel + BLOCK_SIZE - 1) / BLOCK_SIZE;
                adagrad_quant.resize(key_prefix + ".bg_ag", numel, n_blocks);
                adagrad_quant.zero();
            } else {
                accum_f.resize(key_prefix + ".accum", N, L, H, W, C);
                accum_f.zero();
            }
        } else {
            if (quantize_optim) {
                constexpr int64_t BLOCK_SIZE =
                    QuantizedAdamState<8, 256>::kBlockSize;
                int64_t n_blocks = (numel + BLOCK_SIZE - 1) / BLOCK_SIZE;
                quant_state.resize(key_prefix + ".bg_quant", numel, n_blocks);
                quant_state.zero();
            } else {
                g1.resize(key_prefix + ".g1", N, L, H, W, C);
                g1.zero();
                g2.resize(key_prefix + ".g2", N, L, H, W, C);
                g2.zero();
            }
        }
        done = true;
    };
    if (engine().bilagrid_rgb.enabled) {
        setup(engine().bilagrid_rgb.grids,
              engine().bilagrid_rgb.g1, engine().bilagrid_rgb.g2,
              engine().bilagrid_rgb.quant_state,
              engine().bilagrid_rgb.accum_f,
              engine().bilagrid_rgb.adagrad_quant,
              engine().bilagrid_rgb.optim_bits,
              engine().bilagrid_rgb.use_adagrad,
              "eng.bg.rgb", engine().bilagrid_rgb.optim_initialized);
    }
    if (engine().bilagrid_depth.enabled) {
        setup(engine().bilagrid_depth.grids,
              engine().bilagrid_depth.g1, engine().bilagrid_depth.g2,
              engine().bilagrid_depth.quant_state,
              engine().bilagrid_depth.accum_f,
              engine().bilagrid_depth.adagrad_quant,
              engine().bilagrid_depth.optim_bits,
              engine().bilagrid_depth.use_adagrad,
              "eng.bg.depth", engine().bilagrid_depth.optim_initialized);
    }
    if (engine().bilagrid_normal.enabled) {
        setup(engine().bilagrid_normal.grids,
              engine().bilagrid_normal.g1, engine().bilagrid_normal.g2,
              engine().bilagrid_normal.quant_state,
              engine().bilagrid_normal.accum_f,
              engine().bilagrid_normal.adagrad_quant,
              engine().bilagrid_normal.optim_bits,
              engine().bilagrid_normal.use_adagrad,
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
    // applies Adam OR AdaGrad (float or uint8-quantized) in one kernel.
    auto run = [&](DeviceTensor5D<float>& grids,
                   DeviceTensor5D<float>& image_grad,
                   DeviceTensor5D<float>& g1,
                   DeviceTensor5D<float>& g2,
                   QuantizedAdamState<8, 256>& quant_state,
                   DeviceTensor5D<float>& accum_f,
                   QuantizedTensorLog<8, 256>& adagrad_quant,
                   QuantizedTensor<16, 256>& grids_quant,
                   int  optim_bits,
                   int  value_bits,
                   bool use_adagrad,
                   int N, int L, int H, int W, int C,
                   float lr, float tv_weight) {
        if (lr <= 0.0f) return;
        if (image_grad.data_ptr() == nullptr) return;
        const bool quantize_optim = (optim_bits != 32);
        const bool value_quantize = (value_bits != 32);
        // When value_quantize, grids.data_ptr() is null (freed after init);
        // the kernel reads/writes via grids_quant. Otherwise grids is fp32.
        float*    grids_fp32 = value_quantize ? nullptr            : grids.data_ptr();
        uint16_t* grids_q16  = value_quantize ? (uint16_t*)grids_quant.packed_ptr() : nullptr;
        float2*   value_bnd  = value_quantize ? grids_quant.bounds_ptr() : nullptr;
        if (!value_quantize && grids_fp32 == nullptr) return;
        if (use_adagrad) {
            fused_bilagrid_tv_adagrad(
                grids_fp32, grids_q16, value_bnd,
                quantize_optim ? nullptr : accum_f.data_ptr(),
                quantize_optim ? adagrad_quant.packed_ptr() : nullptr,
                quantize_optim ? adagrad_quant.bounds_ptr() : nullptr,
                image_grad.data_ptr(),
                cam_idx_dev,
                N, C_batch, C, L, H, W,
                lr, tv_weight,
                quantize_optim, optim_bits, value_quantize, kBilagridStream);
        } else {
            fused_bilagrid_tv_adam(
                grids_fp32, grids_q16, value_bnd,
                quantize_optim ? nullptr : g1.data_ptr(),
                quantize_optim ? nullptr : g2.data_ptr(),
                quantize_optim ? quant_state.packed_ptr() : nullptr,
                quantize_optim ? quant_state.bounds_ptr() : nullptr,
                image_grad.data_ptr(),
                cam_idx_dev,
                N, C_batch, C, L, H, W,
                lr, tv_weight, adam_step,
                quantize_optim, optim_bits, value_quantize, kBilagridStream);
        }
    };

    auto& brgb = engine().bilagrid_rgb;
    auto& bdep = engine().bilagrid_depth;
    auto& bnor = engine().bilagrid_normal;
    // Use the value-bit-aware shape source for dims (when value-quant, fp32
    // grids is empty so we need to fall back to grids_quant cell count).
    auto get_dims = [](auto& bg, int channels) {
        // Returns (N, L, H, W, C). When grids is alive use its size<>(); else
        // derive from grids_quant.n_cells and channels (grids cells = N*L*H*W*C).
        // We stash L/H/W in bg.value_LHW... actually simpler: keep grids
        // alive as a shape-only descriptor (data_ptr can be null).
        return std::make_tuple(
            (int)bg.grids.template size<0>(),
            (int)bg.grids.template size<1>(),
            (int)bg.grids.template size<2>(),
            (int)bg.grids.template size<3>(),
            (int)bg.grids.template size<4>());
    };
    {
        auto [N, L, H, W, C] = get_dims(brgb, brgb.C);
        run(brgb.grids, brgb.image_grad, brgb.g1, brgb.g2,
            brgb.quant_state, brgb.accum_f, brgb.adagrad_quant,
            brgb.grids_quant,
            brgb.optim_bits, brgb.value_bits, brgb.use_adagrad,
            N, L, H, W, C,
            cfg.lr_rgb, cfg.tv_weight_rgb);
    }
    {
        auto [N, L, H, W, C] = get_dims(bdep, 2);
        run(bdep.grids, bdep.image_grad, bdep.g1, bdep.g2,
            bdep.quant_state, bdep.accum_f, bdep.adagrad_quant,
            bdep.grids_quant,
            bdep.optim_bits, bdep.value_bits, bdep.use_adagrad,
            N, L, H, W, C,
            cfg.lr_depth, cfg.tv_weight_depth);
    }
    {
        auto [N, L, H, W, C] = get_dims(bnor, 3);
        run(bnor.grids, bnor.image_grad, bnor.g1, bnor.g2,
            bnor.quant_state, bnor.accum_f, bnor.adagrad_quant,
            bnor.grids_quant,
            bnor.optim_bits, bnor.value_bits, bnor.use_adagrad,
            N, L, H, W, C,
            cfg.lr_normal, cfg.tv_weight_normal);
    }
}
