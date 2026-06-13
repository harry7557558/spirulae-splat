// Engine PPISP integration (per-camera RGB-only photometric correction).

#include "Engine.h"
#include "EngineCommon.h"
#include "EngineInternal.h"
#include "EngineState.h"

#include "BilagridBindings.h"

#include <array>
#include <stdexcept>
#include <string>
#include <vector>


static constexpr cudaStream_t kPpispStream = (cudaStream_t)0;

// Build a TorchTensorView for the full [N_cam, P] PPISP parameter table.
static TorchTensorView _ppisp_params_full_tv(DeviceTensor2D<float>& table) {
    int64_t N = table.size<0>(), P = table.size<1>();
    return TorchTensorView(
        (uint64_t)table.data_ptr(), 4, {N, P});
}

// Build a TorchTensorView for the current batch's camera_intrins ([C_batch, 4]).
static TorchTensorView _ppisp_intrins_batch_tv() {
    int C_batch = (int)engine().camera.num;
    return TorchTensorView(
        (uint64_t)engine().camera.intrins.data_ptr(), 4, {(int64_t)C_batch, 4LL});
}

// Build a TorchTensorView wrapper around the current batch's cam_indices
// DeviceVector, or null when no cam_indices were supplied (identity mode).
static TorchTensorView _ppisp_cam_indices_tv() {
    const int* ptr = engine().bilagrid_cur_cam_indices.data_ptr();
    if (ptr == nullptr) return TorchTensorView(0, 4, {});
    int64_t n = engine().bilagrid_cur_cam_indices.size();
    return TorchTensorView((uint64_t)ptr, 4, {n});
}

void engine_init_ppisp(int n_grids, std::string param_type, bool use_adagrad) {
    if (n_grids <= 0)
        throw std::runtime_error("engine_init_ppisp: n_grids must be > 0");
    int P;
    if (param_type == "original" || param_type == "") P = 36;
    else if (param_type == "rqs") P = 39;
    else if (param_type == "no_crf") P = 24;
    else throw std::runtime_error(
        "engine_init_ppisp: unknown param_type \"" + param_type +
        "\", must be \"original\", \"rqs\", or \"no_crf\"");

    engine().ppisp.param_type = (param_type == "" ? std::string("original") : param_type);
    engine().ppisp.num_params = P;
    engine().ppisp.use_adagrad = use_adagrad;
    engine().ppisp.params.resize("eng.ppisp.params", n_grids, P);
    if (engine().ppisp.param_type == "original") {
        ppisp_original_default_init(
            engine().ppisp.params.data_ptr(), n_grids, kPpispStream);
    } else {
        engine().ppisp.params.zero();
    }
    engine().ppisp.enabled = true;
    engine().ppisp.optim_initialized = false;
}

// Apply PPISP forward in place on the current rendered RGB. Reads the
// per-image PPISP parameter slot via engine().bilagrid_cur_cam_indices (shared
// with bilagrid). When that index buffer is empty, fall back to identity.
void engine_ppisp_forward(TorchTensorView cam_indices) {
    if (!engine().ppisp.enabled) return;
    // Repopulating with the same tensor is cheap (DeviceVector<int32_t> takes
    // a view of device memory when the source is on device).
    _set_cur_cam_indices(cam_indices);

    int H = engine().camera.height;
    int W = engine().camera.width;
    int C_batch = (int)engine().camera.num;
    if (C_batch <= 0) return;
    if (std::get<0>(engine().fwd.renders).data_ptr() == nullptr)
        throw std::runtime_error("engine_ppisp_forward: forward_3dgs must run first");

    // pre = current renders.rgb (pointer alias, no D2D copy). ppisp_forward
    // already accepts distinct in/out tensors, so we read pre and write a
    // fresh "post" buffer, then re-point engine state.
    auto& fwd_rgb_tensor = std::get<0>(engine().fwd.renders);
    engine().ppisp.fwd_pre = fwd_rgb_tensor;

    DeviceTensor3D<float3> post_rgb;
    post_rgb.resize("eng.ppisp.rgb_post", C_batch, H, W);

    ppisp_forward(
        engine().ppisp.fwd_pre,
        _ppisp_params_full_tv(engine().ppisp.params),
        _ppisp_intrins_batch_tv(),
        (float)W, (float)H,
        engine().ppisp.param_type,
        _ppisp_cam_indices_tv(),
        post_rgb
    );
    fwd_rgb_tensor = post_rgb;
}

// Ensure PPISP optimizer state buffers are allocated + zeroed. Adam uses
// (g1, g2); AdaGrad uses (accum_f). grads is always allocated since both
// paths read the per-camera parameter gradient.
void _ensure_ppisp_optim_state() {
    if (engine().ppisp.optim_initialized) return;
    if (engine().ppisp.params.data_ptr() == nullptr) return;
    int64_t N = engine().ppisp.params.size<0>();
    int64_t P = engine().ppisp.params.size<1>();
    engine().ppisp.grads.resize("eng.ppisp.grads", N, P);
    engine().ppisp.grads.zero();
    if (engine().ppisp.use_adagrad) {
        engine().ppisp.accum_f.resize("eng.ppisp.accum", N, P);
        engine().ppisp.accum_f.zero();
    } else {
        engine().ppisp.g1.resize("eng.ppisp.g1", N, P);
        engine().ppisp.g1.zero();
        engine().ppisp.g2.resize("eng.ppisp.g2", N, P);
        engine().ppisp.g2.zero();
    }
    engine().ppisp.optim_initialized = true;
}

// Backward hook: rewrites v_render_rgb (post-PPISP -> pre-PPISP) in place and
// accumulates the per-camera PPISP parameter gradient into ppisp_grads. The
// kernel writes via atomicAdd, so duplicate cam_indices entries are accumulated
// correctly into the same per-camera slot.
void _engine_ppisp_backward_hook(TorchTensorView v_render_rgb) {
    int H = engine().camera.height;
    int W = engine().camera.width;
    int C_batch = (int)engine().camera.num;
    if (std::get<0>(v_render_rgb) == 0) return;
    if (engine().ppisp.fwd_pre.data_ptr() == nullptr) return;

    TorchTensorView pre_tv(
        (uint64_t)engine().ppisp.fwd_pre.data_ptr(), 4,
        {(int64_t)C_batch, H, W, 3});
    DeviceTensor3D<float3> in_img(pre_tv);
    DeviceTensor3D<float3> v_out_img(v_render_rgb);
    DeviceTensor3D<float3> v_in_img(v_render_rgb);  // overwrite in place

    ppisp_backward(
        in_img,
        _ppisp_params_full_tv(engine().ppisp.params),
        _ppisp_intrins_batch_tv(),
        (float)W, (float)H,
        v_out_img,
        engine().ppisp.param_type,
        _ppisp_cam_indices_tv(),
        v_in_img,
        _ppisp_params_full_tv(engine().ppisp.grads)
    );
}

// Compute the 6 PPISP regularization losses, scaled by their loss weights.
// Returns the [6]-float device buffer (pool-backed). When `compute_grad=true`,
// also accumulates the per-camera parameter gradient into engine().ppisp.grads
// using v_losses = ones (matching training_losses.py).
float* _engine_ppisp_reg_loss_into(
    const std::array<float, (int)PPISPRegLossIndex::length>& loss_weights,
    bool compute_grad
) {
    int N = (int)engine().ppisp.params.size<0>();
    int kRaw =
        (engine().ppisp.param_type == "rqs")    ? (int)RawPPISPRegLossIndexRQS::length :
        (engine().ppisp.param_type == "no_crf") ? (int)RawPPISPRegLossIndexNoCRF::length :
                                                  (int)RawPPISPRegLossIndex::length;
    int kLoss = (int)PPISPRegLossIndex::length;

    // Output losses (zeroed each call so the in-kernel write is a clean store).
    float* losses_buf = DevicePool::global().acquire<float>("eng.ppisp.reg_losses", kLoss);
    cudaMemsetAsync(losses_buf, 0, kLoss * sizeof(float), kPpispStream);

    // Raw losses scratch ([N+1, kRaw], pre-zeroed).
    float* raw_losses_buf = DevicePool::global().acquire<float>(
        "eng.ppisp.reg_raw_losses", (size_t)(N + 1) * kRaw);
    cudaMemsetAsync(raw_losses_buf, 0, (size_t)(N + 1) * kRaw * sizeof(float),
                    kPpispStream);

    TorchTensorView params_tv(
        (uint64_t)engine().ppisp.params.data_ptr(), 4,
        {(int64_t)N, (int64_t)engine().ppisp.num_params});
    TorchTensorView losses_tv((uint64_t)losses_buf, 4, {(int64_t)kLoss});
    TorchTensorView raw_tv((uint64_t)raw_losses_buf, 4,
        {(int64_t)(N + 1), (int64_t)kRaw});

    compute_ppsip_regularization_forward(
        params_tv, loss_weights, engine().ppisp.param_type,
        losses_tv, raw_tv);

    if (compute_grad) {
        // v_losses = ones[kLoss]: gradient flows back through reg-loss sum.
        float* v_losses = DevicePool::global().acquire<float>(
            "eng.ppisp.v_reg_losses", kLoss);
        std::vector<float> h_ones(kLoss, 1.0f);
        cudaMemcpyAsync(v_losses, h_ones.data(), kLoss * sizeof(float),
                        cudaMemcpyHostToDevice, kPpispStream);
        TorchTensorView v_losses_tv((uint64_t)v_losses, 4, {(int64_t)kLoss});

        // Accumulate into ppisp_grads (the regularization backward writes a
        // fresh tensor, so use a scratch buffer and add in afterwards).
        float* v_params_scratch = DevicePool::global().acquire<float>(
            "eng.ppisp.v_reg_params",
            (size_t)N * engine().ppisp.num_params);
        cudaMemsetAsync(v_params_scratch, 0,
            (size_t)N * engine().ppisp.num_params * sizeof(float), kPpispStream);
        TorchTensorView v_params_tv(
            (uint64_t)v_params_scratch, 4,
            {(int64_t)N, (int64_t)engine().ppisp.num_params});

        compute_ppsip_regularization_backward(
            params_tv, loss_weights, raw_tv, v_losses_tv,
            engine().ppisp.param_type, v_params_tv);

        // ppisp_grads += v_params_scratch (over all N * P floats).
        size_t total = (size_t)N * engine().ppisp.num_params;
        ppisp_add_into_grad(v_params_scratch,
            engine().ppisp.grads.data_ptr(), total, kPpispStream);
    }

    return losses_buf;
}

void engine_ppisp_optim_step(int step, const PpispStepConfig& cfg) {
    if (!engine().ppisp.enabled) return;
    _ensure_ppisp_optim_state();
    if (cfg.lr <= 0.0f) {
        // Still zero accumulated grads for next iteration to match other paths.
        engine().ppisp.grads.zero();
        return;
    }

    // Fold the 6 regularization losses' gradients into ppisp_grads.
    bool any_reg = false;
    for (float w : cfg.reg_weights) if (w > 0.0f) { any_reg = true; break; }
    if (any_reg) {
        (void)_engine_ppisp_reg_loss_into(cfg.reg_weights, /*compute_grad=*/true);
    }

    int64_t N = engine().ppisp.params.size<0>();
    int64_t P = engine().ppisp.params.size<1>();
    int64_t numel = N * P;
    auto flat_view = [numel](float* p) {
        TorchTensorView tv((uint64_t)p, 4, {numel, 1LL});
        return DeviceTensorFloatND(tv);
    };
    if (engine().ppisp.use_adagrad) {
        fused_adagrad_step(
            flat_view(engine().ppisp.params.data_ptr()),
            flat_view(engine().ppisp.grads.data_ptr()),
            flat_view(engine().ppisp.accum_f.data_ptr()),
            cfg.lr);
    } else {
        DeviceVector<int32_t> no_per_splat_steps;
        fused_adam_step(
            numel,
            flat_view(engine().ppisp.params.data_ptr()),
            flat_view(engine().ppisp.grads.data_ptr()),
            flat_view(engine().ppisp.g1.data_ptr()),
            flat_view(engine().ppisp.g2.data_ptr()),
            cfg.lr, step + 1, no_per_splat_steps,
            /*l2_reg=*/0.0f, /*l2_reg_offset=*/0.0f,
            /*grad_scale=*/1.0f, /*zero_grad=*/false);
    }
    engine().ppisp.grads.zero();
}
