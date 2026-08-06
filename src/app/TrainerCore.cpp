// TrainerCore.cpp -- see TrainerCore.h. Function-level comments still cite the
// Python trainer this began as (model.py, trainer.py, ...); those files are
// gone, and the citations are kept only as provenance for the numerics.

#include "app/TrainerCore.h"
#include "app/EvalMetrics.h"
#include "checkpoint/Adapt.h"
#include "checkpoint/Resume.h"
#include "i18n/catalog/Log.h"
#include "data/Knn.h"

#ifndef _WIN32
#include <ftw.h>
#endif

#include "external/stb_image_write.h"

#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <ctime>
#include <deque>
#include <numeric>
#include <random>
#include <stdexcept>
#include <thread>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef SS_BACKEND_VULKAN
#include <cuda_runtime.h>  // check_cuda_runtime() driver/runtime preflight
#endif

namespace fs = std::filesystem;
namespace lmsg = spirula::i18n::msg::log;

// Progress lines carrying a path or a count. See i18n/Message.h: whole
// sentences with {0} placeholders, never concatenated fragments.
static std::string lfmt(const spirula::i18n::Msg& m,
                        std::initializer_list<spirula::i18n::Arg> a) {
    return spirula::i18n::format(m, a);
}

namespace spirula {

// Recursive delete. Was nftw rather than std::filesystem::remove_all because
// libtorch interposed its own std::filesystem symbols; torch is gone, so this
// can become remove_all whenever someone wants to.
static void remove_tree(const std::filesystem::path& p) {
#ifndef _WIN32
    nftw(p.string().c_str(),
         [](const char* f, const struct stat*, int, struct FTW*) {
             return ::remove(f);
         }, 16, FTW_DEPTH | FTW_PHYS);
#else
    std::filesystem::remove_all(p);
#endif
}

// ===========================================================================
// Color-space handling
// ===========================================================================

Mat3f gamut_to_rec709(const std::string& name) {
    if (name.empty()) return {1,0,0, 0,1,0, 0,0,1};
    if (name == "ACES2065-1") return {
        2.5247180476f, -1.1325619434f, -0.3921561044f,
        -0.2776344819f, 1.3709123773f, -0.0932778953f,
        -0.0165202369f, -0.1479259606f, 1.1644461975f};
    if (name == "ACEScg") return {
        1.7072552160f, -0.6200352595f, -0.0872199564f,
        -0.1311566587f, 1.1391010566f, -0.0079443978f,
        -0.0245499075f, -0.1248045805f, 1.1493544880f};
    if (name == "Rec.2020") return {
        1.6604910021f, -0.5876411388f, -0.0728498633f,
        -0.1245504745f, 1.1328998971f, -0.0083494226f,
        -0.0181507634f, -0.1005788980f, 1.1187296614f};
    if (name == "AdobeRGB") return {
        1.3983671735f, -0.3983451225f, 0.0000054016f,
        -0.0000103176f, 0.9999916496f, -0.0000039459f,
        -0.0000003709f, -0.0429269510f, 1.0429319656f};
    if (name == "DCI-P3") return {
        1.1548337042f, -0.1451763523f, -0.0096573518f,
        -0.0393300117f, 1.0378282998f, 0.0015017119f,
        -0.0184786235f, -0.0689101110f, 1.0873887345f};
    throw std::runtime_error("unsupported color gamut: " + name);
}

Mat3f invert3x3(const Mat3f& m) {
    float a = m[0], b = m[1], c = m[2],
          d = m[3], e = m[4], g = m[5],
          h = m[6], i = m[7], j = m[8];
    float det = a*(e*j - g*i) - b*(d*j - g*h) + c*(d*i - e*h);
    if (std::abs(det) < 1e-20f) throw std::runtime_error("singular color matrix");
    float inv = 1.0f / det;
    return {(e*j - g*i)*inv, (c*i - b*j)*inv, (b*g - c*e)*inv,
            (g*h - d*j)*inv, (a*j - c*h)*inv, (c*d - a*g)*inv,
            (d*i - e*h)*inv, (b*h - a*i)*inv, (a*e - b*d)*inv};
}

// model.py populate_modules:612-623, applied to the config in place there;
// pure-function port here.
ColorResolution resolve_color(const TrainConfig& c) {
    ColorResolution r;
    r.image_gamut  = c.image_color_gamut;
    r.image_linear = c.image_color_is_linear;
    std::optional<bool> convert = c.convert_initial_point_cloud_color;

    r.splat_gamut = c.splat_color_gamut;
    if (r.splat_gamut.empty()) r.splat_gamut = r.image_gamut;
    else if (!convert.has_value()) convert = true;
    if (r.splat_gamut == "Rec.709") r.splat_gamut = "";

    if (c.splat_color_is_linear.has_value()) {
        r.splat_linear = *c.splat_color_is_linear;
        if (!convert.has_value()) convert = true;
    } else {
        r.splat_linear = r.image_linear;
    }
    r.convert_seed = convert.value_or(false);
    return r;
}


// ===========================================================================
// LR schedule -- port of OptimizerConfig.get_scheduled_lr (optimizer.py:54).
// `lr_final` mirrors getattr(name+"_lr_final", lr): pass the field when the
// Python class has one, std::nullopt when the attribute doesn't exist (which
// getattr defaults to lr -> same constant result as nullopt here).
// ===========================================================================

float scheduled_lr(int step, int max_steps, float lr,
                   std::optional<float> lr_final,
                   std::optional<int> warmup) {
    float s = lr;
    if (lr_final.has_value() && lr != 0.0f && *lr_final != 0.0f)
        s = lr * std::pow(*lr_final / lr,
                          std::min((float)step / (float)std::max(max_steps, 1), 1.0f));
    if (warmup.has_value())
        s = std::min(s, lr * std::min((float)step / (float)std::max(*warmup, 1), 1.0f));
    return s;
}


// ===========================================================================
// Splat seeding -- port of model.py populate_modules:546-651 (3dgs branch)
// ===========================================================================

SeedSplats seed_splats(const ColmapPoints3D& pts, const TrainConfig& cfg,
                       const ColorResolution& color) {
    std::mt19937 rng(42);
    std::normal_distribution<float> gauss(0.f, 1.f);
    std::uniform_real_distribution<float> uni(0.f, 1.f);

    float scale_init   = cfg.scale_init.value_or(0.5f);     // model.py:569-573
    float opacity_init = cfg.opacity_init.value_or(0.1f);

    // Resolve seed count into [min_init, cap_max] (model.py:547-559).
    int64_t n_src = pts.num();
    if (n_src == 0) throw std::runtime_error("seed_splats: empty point cloud");
    int64_t min_init = std::max<int64_t>(
        (int64_t)(std::min(cfg.min_init_fraction, 1.0f) * cfg.cap_max), 1);
    min_init = std::min<int64_t>(min_init, cfg.cap_max);

    std::vector<int64_t> pick;
    if (n_src > cfg.cap_max) {
        pick.resize(n_src);
        std::iota(pick.begin(), pick.end(), 0);
        std::shuffle(pick.begin(), pick.end(), rng);
        pick.resize(cfg.cap_max);
    } else {
        // Repeat modulo when under min_init. TODO: the Python path jitters
        // repeats toward a nearest neighbor (model.py:550-556) instead of
        // duplicating exactly.
        int64_t n = std::max(n_src, min_init);
        pick.resize(n);
        for (int64_t i = 0; i < n; i++) pick[i] = i % n_src;
    }
    const int64_t num = (int64_t)pick.size();
    const int64_t cap = cfg.preallocate_splat_tensors && cfg.use_mcmc
        ? std::max<int64_t>(num, cfg.cap_max) : num;
    const int64_t dim_sh = (int64_t)(cfg.sh_degree + 1) * (cfg.sh_degree + 1);

    SeedSplats s;
    s.num = num;
    s.means.assign(cap * 3, 0.f);
    s.quats.assign(cap * 4, 0.f);
    s.scales.assign(cap * 3, 0.f);
    s.opacities.assign(cap * 1, 0.f);
    s.features_dc.assign(cap * 3, 0.f);
    s.features_sh.assign(cap * (dim_sh - 1) * 3, 0.f);

    // means (*= relative_scale, model.py:561-562)
    float rescale = cfg.relative_scale.value_or(1.0f);
    for (int64_t i = 0; i < num; i++)
        for (int d = 0; d < 3; d++)
            s.means[i*3 + d] = pts.xyz[pick[i]*3 + d] * rescale;

    // log(scale_init * sqrt(mean d^2 of 4-NN)) over xyz (model.py:578-583).
    // TODO: suppress_initial_scales (model.py:584-585).
    std::vector<float> nn = knn::mean_knn_dist(s.means, num, 4);
    for (int64_t i = 0; i < num; i++) {
        float v = std::log(scale_init * nn[i] + 1e-8f);
        s.scales[i*3+0] = s.scales[i*3+1] = s.scales[i*3+2] = v;
    }

    for (int64_t i = 0; i < num; i++) {
        float q[4] = {gauss(rng), gauss(rng), gauss(rng), gauss(rng)};
        float qn = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
        for (int d = 0; d < 4; d++) s.quats[i*4 + d] = q[d] / std::max(qn, 1e-12f);
        s.opacities[i] = std::log(opacity_init / (1.f - opacity_init));
    }

    // Seed colors (model.py:596-635). Uniform-color clouds are randomized.
    bool all_same = true;
    for (int64_t i = 0; i < num && all_same; i++)
        for (int d = 0; d < 3; d++)
            if (pts.rgb[pick[i]*3 + d] != pts.rgb[pick[0]*3]) { all_same = false; break; }
    Mat3f inv_image = invert3x3(gamut_to_rec709(color.image_gamut));
    for (int64_t i = 0; i < num; i++) {
        float col[3];
        for (int d = 0; d < 3; d++)
            col[d] = all_same ? uni(rng) : pts.rgb[pick[i]*3 + d] / 255.f;
        if (color.convert_seed) {
            // sRGB -> linear (thresholds as in model.py:626/631).
            if (color.splat_linear || !color.splat_gamut.empty())
                for (int d = 0; d < 3; d++)
                    col[d] = col[d] < 0.055f ? col[d] / 12.92f
                        : std::pow(std::max((col[d] + 0.055f) / 1.055f, 0.f), 2.4f);
            if (!color.splat_gamut.empty()) {
                float t[3] = {col[0], col[1], col[2]};
                for (int r = 0; r < 3; r++)
                    col[r] = inv_image[r*3+0]*t[0] + inv_image[r*3+1]*t[1] + inv_image[r*3+2]*t[2];
                if (!color.splat_linear)
                    for (int d = 0; d < 3; d++)
                        col[d] = col[d] < 0.0031308f ? 12.92f * col[d]
                            : 1.055f * std::pow(std::max(col[d], 0.f), 1.f/2.4f) - 0.055f;
            }
        }
        for (int d = 0; d < 3; d++)
            s.features_dc[i*3 + d] = (col[d] - 0.5f) / 0.28209479177387814f;
    }
    return s;
}


// ===========================================================================
// Per-step EngineStepConfig
// ===========================================================================

namespace {

int densify_loss_map_mode_int(const std::string& mode) {
    // Mirrors _DENSIFY_LOSS_MAP_MODE_TO_INT (model.py:437).
    if (mode == "none")              return 0;
    if (mode == "loss_full")         return 1;
    if (mode == "ssim_full")         return 2;
    if (mode == "ssim_cs")           return 3;
    if (mode == "ssim_structure")    return 4;
    if (mode == "edge_aware")        return 5;
    if (mode == "robust_edge_aware") return 6;
    throw std::runtime_error("unknown densify_loss_map_mode: " + mode);
}

}  // namespace

// model.py:1483 _build_loss_weights.
std::array<float, (int)LossWeightIndex::length>
build_loss_weights(const TrainConfig& c, int step) {
    float dist_factor = std::min((float)step / std::max(c.distortion_reg_warmup, 1), 1.0f);
    float reg_active  = step >= c.reg_warmup_length ? 1.0f : 0.0f;
    float sup_active  = step > c.supervision_warmup ? 1.0f : 0.0f;
    float median_factor = std::min((float)step / std::max(c.median_warmup, 1), 1.0f);
    float alpha_reg_factor = c.alpha_reg_weight *
        std::min((float)step / std::max(c.alpha_reg_warmup, 1), 1.0f);
    float mask = c.apply_loss_for_mask ? 1.0f : 0.0f;

    float w_rgb_l1 = std::max(0.0f, c.l1_weight);
    float w_rgb_l2 = std::max(0.0f, c.l2_weight);
    float w_y_l1   = std::max(0.0f, c.l1_weight_y);
    float w_y_l2   = std::max(0.0f, c.l2_weight_y);
    float w_u_l2   = std::max(0.0f, c.l2_weight_u);
    float w_v_l2   = std::max(0.0f, c.l2_weight_v);
    float total = w_rgb_l1 + w_rgb_l2 + w_y_l1 + w_y_l2 + w_u_l2 + w_v_l2;
    float scale = total > 0.0f ? (1.0f - c.ssim_lambda) / total : 0.0f;

    std::array<float, (int)LossWeightIndex::length> w{};
    w[(int)LossWeightIndex::RgbSupL1]      = w_rgb_l1 * scale;
    w[(int)LossWeightIndex::RgbSupL2]      = w_rgb_l2 * scale;
    w[(int)LossWeightIndex::YSupL1]        = w_y_l1 * scale;
    w[(int)LossWeightIndex::YSupL2]        = w_y_l2 * scale;
    w[(int)LossWeightIndex::USupL2]        = w_u_l2 * scale;
    w[(int)LossWeightIndex::VSupL2]        = w_v_l2 * scale;
    w[(int)LossWeightIndex::DepthSup]      = sup_active * c.depth_supervision_weight;
    w[(int)LossWeightIndex::NormalSup]     = sup_active * c.normal_supervision_weight;
    w[(int)LossWeightIndex::AlphaSup]      = mask * c.alpha_loss_weight;
    w[(int)LossWeightIndex::AlphaSupUnder] = mask * c.alpha_loss_weight_under;
    w[(int)LossWeightIndex::NormalReg]     = reg_active * c.normal_reg_weight * dist_factor;
    w[(int)LossWeightIndex::AlphaReg]      = reg_active * alpha_reg_factor;
    w[(int)LossWeightIndex::RgbDistReg]    = reg_active * c.rgb_distortion_reg * dist_factor;
    w[(int)LossWeightIndex::DepthDistReg]  = reg_active * c.depth_distortion_reg * dist_factor;
    w[(int)LossWeightIndex::NormalDistReg] = reg_active * c.normal_distortion_reg * dist_factor;
    w[(int)LossWeightIndex::MeanMedianDepthSup]    = median_factor * c.mean_median_depth_weight;
    w[(int)LossWeightIndex::MedianDepthNormalReg]  = median_factor * c.median_depth_normal_reg_weight;
    w[(int)LossWeightIndex::MedianNormalSup]       = median_factor * c.median_normal_supervision_weight;
    w[(int)LossWeightIndex::MedianRenderNormalReg] = median_factor * c.median_render_normal_reg_weight;
    return w;
}

EngineStepConfig build_step_config(const TrainConfig& c, const RunState& st, int step) {
    int max_steps_lr = c.max_steps.value_or(c.num_iterations);
    float alpha = st.train_frame_scale;
    EngineStepConfig cfg;

    // ---- loss (model.py engine_train_step_managed:1903-1915) --------------
    cfg.loss.weights = build_loss_weights(c, step);
    cfg.loss.w_ssim = c.ssim_lambda;
    cfg.loss.num_loss_scales = c.num_loss_scales + 1;
    cfg.loss.loss_scale_min_pixels = c.loss_scale_min_pixels;
    int loss_map_mode = densify_loss_map_mode_int(c.densify_loss_map_mode);
    // blend >= 1: world-grad score only, loss map has no consumer.
    if (c.densify_score_blend_world_grad >= 1.0f) loss_map_mode = 0;
    cfg.loss.loss_map_mode = loss_map_mode;
    cfg.loss.compute_loss_map = (loss_map_mode != 0);
    cfg.loss.robust_edge_aware_quantile = c.densify_robust_edge_aware_quantile;
    cfg.loss.overexposure_reg_weight = c.overexposure_reg;
    if (st.bilagrid_rgb_init || st.ppisp_init) {
        cfg.loss.color_shift_reg_weight = c.color_shift_reg_weight;
        cfg.loss.color_shift_reg_beta =
            std::max(0.0f, 1.0f - 1.0f / std::max(c.color_shift_reg_ema_period, 1));
    }
    cfg.loss.input_depth_is_ray_depth = c.input_depth_is_ray_depth;

    // ---- optim (core.py _build_optim_config:280) ---------------------------
    float means_lr = scheduled_lr(step, max_steps_lr, c.means_lr, c.means_lr_final);
    if (!c.use_scale_agnostic_mean) means_lr *= alpha;
    cfg.optim.lr_means       = means_lr;
    cfg.optim.lr_quats       = scheduled_lr(step, max_steps_lr, c.quats_lr);
    cfg.optim.lr_scales      = scheduled_lr(step, max_steps_lr, c.scales_lr, c.scales_lr_final);
    cfg.optim.lr_opacities   = scheduled_lr(step, max_steps_lr, c.opacities_lr);
    cfg.optim.lr_features_dc = scheduled_lr(step, max_steps_lr, c.features_dc_lr);
    cfg.optim.lr_features_sh = scheduled_lr(step, max_steps_lr, c.features_sh_lr);
    cfg.optim.max_gauss_ratio             = c.max_gauss_ratio;
    cfg.optim.scale_regularization_weight = c.scale_regularization_weight;
    cfg.optim.mcmc_opacity_reg_weight     = c.opacity_reg;
    cfg.optim.mcmc_scale_reg_weight       = c.scale_reg / alpha;
    cfg.optim.erank_reg_weight            = c.erank_reg;
    cfg.optim.erank_reg_weight_s3         = c.erank_reg_s3;
    cfg.optim.quat_norm_reg_weight        = c.quat_norm_reg;
    cfg.optim.sh_reg_weight               = c.sh_reg;
    cfg.optim.use_scale_agnostic_mean     = c.use_scale_agnostic_mean;
    // Renderer.__init__ level -> bits mapping (core.py:95-102).
    cfg.optim.quantization_level = c.quantization_level;
    cfg.optim.sh_optim_bits      = c.quantization_level == 0 ? 32 : 8;
    cfg.optim.sh_value_bits      = c.quantization_level == 0 ? 32 : 16;
    cfg.optim.non_sh_optim_bits  = c.quantization_level == 0 ? 32 : 16;
    cfg.optim.use_per_splat_bias_correction = c.use_per_splat_bias_correction;
    cfg.optim.use_fused_proj_bwd_optim      = c.use_fused_proj_bwd_optim;
    cfg.optim.write_densify_world_grad_score =
        c.densify_score_blend_world_grad > 0.0f && c.use_revised_densification;
    cfg.optim.split_batch     = c.split_batch;
    cfg.optim.color_is_linear = st.splat_linear;
    cfg.optim.use_color_trust_region = st.splat_linear;
    cfg.optim.eps_tr = 1e-6f * std::pow(0.01f, (float)step / std::max(max_steps_lr, 1));

    // ---- densify (core.py _build_densify_config:321) -----------------------
    float noise_lr_scalar = c.use_revised_densification ? 1.0f : alpha;
    cfg.densify.refine_start_iter             = c.refine_start_iter;
    cfg.densify.refine_stop_num_iter          = c.refine_stop_num_iter;
    cfg.densify.refine_stop_iter              = c.refine_stop_iter;
    cfg.densify.refine_every                  = c.refine_every;
    cfg.densify.growth_factor                 = c.growth_factor;
    cfg.densify.min_opacity                   = c.min_opacity;
    cfg.densify.max_screen_size               = c.max_screen_size;
    cfg.densify.max_screen_size_clip_hardness = c.max_screen_size_clip_hardness;
    cfg.densify.max_world_size                = c.max_world_size * alpha;
    cfg.densify.noise_lr                      = c.noise_lr * noise_lr_scalar;
    cfg.densify.noise_lr_final                = c.noise_lr_final * noise_lr_scalar;
    cfg.densify.use_revised_densification     = c.use_revised_densification;
    cfg.densify.score_mode = c.densify_score_mode == "mean" ? 0
        : c.densify_score_mode == "max" ? 1
        : c.densify_score_mode == "median" ? 2 : 3;
    cfg.densify.score_blend_world_grad = c.densify_score_blend_world_grad;
    cfg.densify.las_split_opacity_k_init   = c.long_axis_split_opacity_k[0];
    cfg.densify.las_split_opacity_k_final  = c.long_axis_split_opacity_k[1];
    cfg.densify.las_split_opacity_k_warmup = (int)c.long_axis_split_opacity_k[2];

    // ---- bilagrid LRs + TV (model.py:1931-1961) ----------------------------
    if (st.bilagrid_rgb_init) {
        cfg.bilagrid.lr_rgb = c.use_adagrad_bilagrid_optim
            ? c.bilagrid_adagrad_lr
            : scheduled_lr(step, max_steps_lr, c.bilagrid_lr, c.bilagrid_lr_final,
                           c.bilagrid_lr_warmup);
        cfg.bilagrid.tv_weight_rgb = c.bilagrid_tv_loss_weight;
    }
    if (st.bilagrid_depth_init) {
        cfg.bilagrid.lr_depth = c.use_adagrad_bilagrid_optim
            ? c.bilagrid_adagrad_depth_lr
            : scheduled_lr(step, max_steps_lr, c.bilagrid_depth_lr,
                           c.bilagrid_depth_lr_final, c.bilagrid_depth_lr_warmup);
        cfg.bilagrid.tv_weight_depth = c.bilagrid_tv_loss_weight_geometry;
    }
    if (st.bilagrid_normal_init) {
        cfg.bilagrid.lr_normal = c.use_adagrad_bilagrid_optim
            ? c.bilagrid_adagrad_normal_lr
            : scheduled_lr(step, max_steps_lr, c.bilagrid_normal_lr,
                           c.bilagrid_normal_lr_final, c.bilagrid_normal_lr_warmup);
        cfg.bilagrid.tv_weight_normal = c.bilagrid_tv_loss_weight_geometry;
    }

    // ---- PPISP (model.py:1980-1997) ----------------------------------------
    if (st.ppisp_init) {
        cfg.ppisp.lr = c.use_adagrad_ppisp_optim
            ? c.ppisp_adagrad_lr
            : scheduled_lr(step, max_steps_lr, c.ppisp_lr, c.ppisp_lr_final,
                           c.ppisp_lr_warmup);
        // PPISPRegLossIndex order (core.py:449-456).
        cfg.ppisp.reg_weights = {
            c.ppisp_reg_exposure_mean, c.ppisp_reg_vig_center,
            c.ppisp_reg_vig_non_pos,   c.ppisp_reg_vig_channel_var,
            c.ppisp_reg_color_mean,    c.ppisp_reg_crf_channel_var,
        };
    }
    // Outside the guard: it is an ordering flag, not a rate, so it reflects
    // the config whether or not PPISP is live (matching the Python path).
    cfg.ppisp.run_before_bilagrid = c.apply_ppisp_before_bilagrid;

    // ---- background (model.py:1963-1977) -----------------------------------
    if (c.background_mode == "noise") {
        float rw = std::min((float)step / std::max(c.background_noise_warmup, 1), 1.0f);
        cfg.background.randomize_weight =
            1.0f - (1.0f - c.background_noise_pre_warmup) * (1.0f - rw);
    } else if (c.background_mode == "sh") {
        cfg.background.lr_dc = scheduled_lr(step, max_steps_lr, c.background_dc_lr);
        cfg.background.lr_sh = scheduled_lr(step, max_steps_lr, c.background_sh_lr);
    }
    cfg.background.seed = (uint32_t)(step & 0x7FFFFFFF);

    return cfg;
}


// ===========================================================================
// config.json dump
// ===========================================================================

namespace {

// JSON value emitters (Python-json compatible: float('inf') serializes as
// Infinity, None as null).
std::string json_str(bool v)               { return v ? "true" : "false"; }
std::string json_str(int v)                { return std::to_string(v); }
std::string json_str(float v) {
    if (std::isinf(v)) return v > 0 ? "Infinity" : "-Infinity";
    char buf[32]; std::snprintf(buf, sizeof buf, "%g", v); return buf;
}
std::string json_str(const std::string& v) {
    if (v.empty()) return "null";
    std::string s = "\"";
    for (char ch : v) { if (ch == '"' || ch == '\\') s += '\\'; s += ch; }
    return s + "\"";
}
template <typename T> std::string json_str(const std::optional<T>& v) {
    return v ? json_str(*v) : "null";
}
template <typename T, size_t N> std::string json_str(const std::array<T, N>& v) {
    std::string s = "[";
    for (size_t i = 0; i < N; i++) s += (i ? ", " : "") + json_str(v[i]);
    return s + "]";
}

}  // namespace

// Nested-by-group dump. Keys are train_json_key(flag) -- see the note on it
// in config/TrainConfig.h; this is a read-back format (`spirula mesh` and
// --resume parse it), so the key set is not free to change.
void save_config_json(const TrainConfig& c, const fs::path& out_dir,
                      const std::string& preset) {
    FILE* f = std::fopen((out_dir / "config.json").string().c_str(), "w");
    if (!f) throw std::runtime_error("cannot write config.json");
    std::fprintf(f, "{\n    \"preset\": \"%s\"", preset.c_str());
    const char* open_group = "";
#define SS_DUMP(type, member, default_, group, choices, help)                  \
    {                                                                          \
    const char* key = train_json_key(#member);                                \
    if (std::strcmp(group, "trainer") == 0) {                                  \
        std::fprintf(f, ",\n    \"%s\": %s", key, json_str(c.member).c_str()); \
    } else {                                                                   \
        if (std::strcmp(open_group, group) != 0) {                             \
            if (*open_group) std::fprintf(f, "\n    }");                       \
            std::fprintf(f, ",\n    \"%s\": {", group);                        \
            open_group = group;                                                \
            std::fprintf(f, "\n        \"%s\": %s", key, json_str(c.member).c_str()); \
        } else {                                                               \
            std::fprintf(f, ",\n        \"%s\": %s", key, json_str(c.member).c_str()); \
        }                                                                      \
    }                                                                          \
    }
    SS_CONFIG_FIELDS(SS_DUMP)
#undef SS_DUMP
    if (*open_group) std::fprintf(f, "\n    }");
    std::fprintf(f, "\n}\n");
    std::fclose(f);
}


// ===========================================================================
// TrainerSession
// ===========================================================================

void TrainerSession::log(const std::string& msg) {
    if (log_fn) { log_fn(msg); return; }
    std::printf("%s\n", msg.c_str());
    std::fflush(stdout);
}

// Unported-feature guards: fail early rather than ignore a flag.
void TrainerSession::check_config() {
    auto not_impl = [](const std::string& what) {
        throw std::runtime_error(what +
            " is not supported yet");
    };
    if (cfg.use_bvh)                    not_impl("--use-bvh");
    if (cfg.use_camera_optimizer)       not_impl("--use-camera-optimizer");
    if (cfg.deblur_training_images)     not_impl("--deblur-training-images");
    if (!cfg.optimizer_offload.empty()) not_impl("--optimizer-offload");
    if (cfg.cache_images == "gpu")      not_impl("--cache-images gpu");
    if (cfg.rescale_camera_to_fit < 0)  not_impl("--rescale-camera-to-fit auto-detect");
    if (cfg.num_downscales > 0)
        log("warning: --num-downscales is a Python-data-path "
            "feature; ignored by the managed engine path");
    if (cfg.validation_fraction > 0)
        log("warning: validation images are held out but "
            "early stopping / eval is not ported yet");
    if (cfg.orientation_method != "up" || cfg.center_method != "poses")
        log("warning: orientation/center method '" + cfg.orientation_method +
            "'/'" + cfg.center_method + "' approximated as 'up'/'poses' "
            "(affects only train_frame_scale; see docs/notes/pose-normalization.md "
            "for the unported reference implementation)");
    if (cfg.train_frame != "points")
        not_impl("--train-frame " + cfg.train_frame);
    if (cfg.primitive != "3dgs" && cfg.primitive != "mip" && cfg.primitive != "3dgut")
        not_impl("--primitive " + cfg.primitive);
    if (cfg.quantization_level != 0 && cfg.quantization_level != 1)
        throw std::runtime_error("quantization_level must be 0 or 1");
}

void TrainerSession::load_dataset() {
    DatasetParserConfig pcfg;
    pcfg.recon_dir            = cfg.colmap_recon_dir;
    pcfg.image_dir            = cfg.image_dir;
    pcfg.mask_dir             = cfg.mask_dir;
    pcfg.depth_dir            = cfg.depth_dir;
    pcfg.normal_dir           = cfg.normal_dir;
    pcfg.validation_fraction  = cfg.validation_fraction;
    pcfg.eval_mode            = cfg.eval_mode;
    pcfg.eval_interval        = cfg.eval_interval;
    pcfg.train_split_fraction = cfg.train_split_fraction;
    pcfg.outlier_threshold    = cfg.outlier_threshold;
    pcfg.rescale_camera_to_fit   = cfg.rescale_camera_to_fit;
    pcfg.downscale_rounding_mode = cfg.downscale_rounding_mode;
    pcfg.metashape_xml           = cfg.metashape_xml;
    pcfg.metashape_ply           = cfg.metashape_ply;
    pcfg.metashape_psx           = cfg.metashape_psx;
    ds = parse_dataset(cfg.data, pcfg, cfg.data_format);

    // relative_scale: scale the world like model.py:561 (means) +
    // trainer.py:407 (viewmat T; c2w t scaled here, pre-bake).
    // auto_scale_poses=False forces the normalized-frame scale to 1
    // (dataparser.py:466-468).
    if (cfg.relative_scale.has_value()) {
        float rs = *cfg.relative_scale;
        for (auto& v : ds.points.xyz) v *= rs;
        for (int64_t i = 0; i < ds.num_cameras; i++)
            for (int r = 0; r < 3; r++)
                ds.c2w[i*12 + r*4 + 3] *= rs;
    }
    if (!cfg.auto_scale_poses) ds.train_frame_scale = 1.0f;

    // POST-split camera bake (identity when no warp flag applies).
    post = bake_post_split(
        ds, cfg.warp_to_pinhole, cfg.warp_spherical_to_pinhole);

    // Warp-path feature guards (trainer.py:294-303).
    has_depth  = !ds.depth_filenames.empty()  && cfg.load_depths;
    has_normal = !ds.normal_filenames.empty() && cfg.load_normals;
    if (post.direct_equirect && (has_depth || has_normal))
        throw std::runtime_error(
            "Direct equirectangular training (warp_spherical_to_pinhole=0) "
            "does not support depth/normal supervision yet.");

    char scale[32];
    std::snprintf(scale, sizeof scale, "%.4g", ds.train_frame_scale);
    log(lfmt(lmsg::parsed_dataset,
             {(long long)ds.num_cameras, (long long)post.n_post,
              (long long)ds.points.num(), scale}));
}

// Pre-flight GPU check. A binary compiled by a newer CUDA toolkit than the
// installed driver supports links and loads fine, but every kernel launch then
// fails -- and cudaGetErrorString can't name the resulting error, so it shows
// up as the cryptic "CUDA Error ...: (null)". Detect the mismatch up front and
// report it with the numbers and the fix, instead of dying deep in a kernel.
// CUDA-backend-specific diagnostics by design; the Vulkan backend replaces
// this with instance/physical-device enumeration at the same call site.
#ifndef SS_BACKEND_VULKAN
static void check_cuda_runtime() {
    auto fmt = [](int v) {
        return std::to_string(v / 1000) + "." + std::to_string((v % 1000) / 10);
    };
    // cudaRuntimeGetVersion is answered by the statically-linked runtime and
    // does not touch the driver -- it tells us which CUDA toolkit built this
    // binary even when the driver is unusable.
    int runtime_ver = 0;
    cudaRuntimeGetVersion(&runtime_ver);
    std::string built_with = runtime_ver ? " (this binary was built with CUDA "
                                           + fmt(runtime_ver) + ")" : "";

    // The first driver-backed call triggers lazy CUDA init; if the driver is
    // too old for the runtime it fails here with cudaErrorInsufficientDriver
    // instead of much later inside a kernel launch (where cudaGetErrorString
    // returns null -> the cryptic "CUDA Error ...: (null)").
    int dev_count = 0;
    cudaError_t cerr = cudaGetDeviceCount(&dev_count);
    if (cerr == cudaErrorInsufficientDriver) {
        int driver_ver = 0;
        cudaDriverGetVersion(&driver_ver);  // 0 if the driver is far too old
        std::string drv = driver_ver
            ? "The installed NVIDIA driver supports only up to CUDA "
              + fmt(driver_ver) + "."
            : "The installed NVIDIA driver is too old.";
        throw std::runtime_error(
            "NVIDIA driver too old for this build. " + drv + built_with +
            " Update the GPU driver, or rebuild against an older CUDA toolkit that matches the driver.");
    }
    if (cerr != cudaSuccess) {
        const char* s = cudaGetErrorString(cerr);
        throw std::runtime_error(
            std::string("CUDA initialization failed: ") +
            (s ? s : "unknown error") + built_with +
            ". Check the NVIDIA driver / GPU installation.");
    }
    if (dev_count == 0)
        throw std::runtime_error("No CUDA-capable GPU detected.");
}
#endif  // SS_BACKEND_VULKAN

void TrainerSession::setup_engine() {
#ifndef SS_BACKEND_VULKAN
    check_cuda_runtime();
#endif

    // ---- Output dir (trainer.py _setup_output_dir:524) ----------------------
    if (!out_dir_override.empty()) {
        out_dir = fs::path(out_dir_override);
    } else if (!cfg.output_dir_name.empty()) {
        out_dir = fs::path(cfg.output_dir_prefix) / cfg.output_dir_name;
    } else {
        std::time_t t = std::time(nullptr);
        char stamp[32];
        std::strftime(stamp, sizeof stamp, "%Y%m%d-%H%M%S", std::localtime(&t));
        out_dir = fs::path(cfg.output_dir_prefix) /
                  (fs::path(cfg.data).stem().string() + "_" + stamp);
    }
    fs::create_directories(out_dir);
    if (write_config_json)
        save_config_json(cfg, out_dir, preset);
    log(lfmt(lmsg::output_directory, {fs::absolute(out_dir).string()}));

    // ---- Engine setup -------------------------------------------------
    engine_reset();

    ColorResolution color = resolve_color(cfg);
    SeedSplats seed = seed_splats(ds.points, cfg, color);
    int64_t cap = (int64_t)seed.opacities.size();
    int64_t dim_sh = (int64_t)(cfg.sh_degree + 1) * (cfg.sh_degree + 1);
    auto tv = [](std::vector<float>& v, std::vector<int64_t> shape) -> TorchTensorView {
        return {(uint64_t)(uintptr_t)v.data(), (uint32_t)sizeof(float), std::move(shape)};
    };
    set_data_3dgs(seed.num,
                  tv(seed.means,       {cap, 3}),
                  tv(seed.quats,       {cap, 4}),
                  tv(seed.scales,      {cap, 3}),
                  tv(seed.opacities,   {cap, 1}),
                  tv(seed.features_dc, {cap, 3}),
                  tv(seed.features_sh, {cap, dim_sh - 1, 3}));

    // Background blending (model.py:516-526).
    if (cfg.background_mode == "noise")
        engine_init_background_noise(color.splat_linear);
    else if (cfg.background_mode == "sh")
        engine_init_background_sh(cfg.background_sh_degree, color.splat_linear);

    // Linear / wide-gamut color space (model.py:532-543).
    bool splat_cs_on = color.splat_linear || !color.splat_gamut.empty();
    bool image_cs_on = color.image_linear || !color.image_gamut.empty();
    {
        auto vec = [](const Mat3f& m) { return std::vector<float>(m.begin(), m.end()); };
        engine_init_color_space(
            splat_cs_on, color.splat_linear,
            splat_cs_on ? vec(gamut_to_rec709(color.splat_gamut)) : std::vector<float>{},
            image_cs_on, color.image_linear,
            image_cs_on ? vec(gamut_to_rec709(color.image_gamut)) : std::vector<float>{});
    }

    // ---- DataManager (trainer.py _setup_cpp_data_manager) -------------------
    const int64_t N = ds.num_cameras;
    int64_t num_val = (int64_t)ds.val_indices.size();
    int64_t num_train = N - num_val;
    // Batch-size policy (datamanager.py:91-105, stochastic=False).
    double n_batch = std::max((double)num_train / std::max(cfg.max_batch_per_epoch, 1), 1.0);
    int train_bs = std::max(1, (int)(n_batch + 0.5));
    int val_bs = 1;
    if (num_val > 0)
        val_bs = std::max(1, (int)std::ceil(n_batch * (double)num_val / (double)num_train));

    DataManagerConfig dm;
    dm.cache_mode  = (cfg.cache_images == "disk") ? CacheMode::DISK : CacheMode::CPU;
    // Enable masks also when the fisheye warp is active with no on-disk
    // masks -- the DataManager synthesizes a 1x1 white placeholder per
    // such image so the wide-warp kernel produces the post-split FOV
    // mask (1 inside the lens circle, 0 outside). Without it the unseen
    // face regions train as black (trainer.py:452-461).
    dm.load_masks  = !ds.mask_filenames.empty() || post.any_fisheye_warp;
    dm.load_depths      = has_depth;
    dm.load_normals     = has_normal;
    dm.train_batch_size = train_bs;
    dm.val_batch_size   = val_bs;
    dm.mask_boundary_offset = cfg.mask_boundary_offset;
    dm.warp_to_pinhole  = cfg.warp_to_pinhole;
    engine_setup_data_manager(
        dm, ds.camera_models,
        ds.image_filenames, ds.mask_filenames,
        has_depth ? ds.depth_filenames : std::vector<std::string>{},
        has_normal ? ds.normal_filenames : std::vector<std::string>{},
        ds.widths, ds.heights,
        post.any_warp ? post.K_per_camera : std::vector<int32_t>{},
        post.any_warp ? post.post_offsets : std::vector<int32_t>{},
        post.viewmats, post.intrins, post.dist_coeffs,
        post.input_intrins, post.input_dist_coeffs,
        ds.train_indices, ds.val_indices);

    // ---- Bilagrid / PPISP init (model.py _maybe_init_bilagrid:823). --------
    // Enablement conditions are static here (dataset modalities known up
    // front), so the lazy per-step init collapses to setup time.
    st = RunState{};
    st.train_frame_scale = ds.train_frame_scale;
    st.splat_linear      = color.splat_linear;
    int optim_bits = cfg.quantization_level == 0 ? 32 : 8;
    int value_bits = cfg.quantization_level == 0 ? 32 : 16;
    // num_train_data resolves to the POST-split camera count -- the
    // bilagrid / PPISP tables have one slot per post camera and the
    // TV-loss normalization depends on it (trainer.py:486-494).
    int n_grids = (int)post.n_post;

    if (cfg.use_bilateral_grid &&
        (cfg.use_adagrad_bilagrid_optim ? cfg.bilagrid_adagrad_lr
                                        : cfg.bilagrid_lr) > 0.0f) {
        // bilagrid_shape is (X, Y, W) -> engine (L=W, H=Y, W=X) (model.py:847-851).
        engine_init_bilagrid_rgb(n_grids, cfg.bilagrid_type,
                                 cfg.bilagrid_shape[2], cfg.bilagrid_shape[1],
                                 cfg.bilagrid_shape[0],
                                 optim_bits, value_bits,
                                 cfg.use_adagrad_bilagrid_optim);
        st.bilagrid_rgb_init = true;
    }
    if (cfg.use_bilateral_grid_for_geometry && has_depth &&
        cfg.depth_supervision_weight > 0.0f &&
        (cfg.use_adagrad_bilagrid_optim ? cfg.bilagrid_adagrad_depth_lr
                                        : cfg.bilagrid_depth_lr) > 0.0f) {
        engine_init_bilagrid_depth(n_grids,
                                   cfg.bilagrid_shape_geometry[2],
                                   cfg.bilagrid_shape_geometry[1],
                                   cfg.bilagrid_shape_geometry[0],
                                   optim_bits, value_bits,
                                   cfg.use_adagrad_bilagrid_optim);
        st.bilagrid_depth_init = true;
    }
    if (cfg.use_bilateral_grid_for_geometry && has_normal &&
        cfg.normal_supervision_weight > 0.0f &&
        (cfg.use_adagrad_bilagrid_optim ? cfg.bilagrid_adagrad_normal_lr
                                        : cfg.bilagrid_normal_lr) > 0.0f) {
        engine_init_bilagrid_normal(n_grids,
                                    cfg.bilagrid_shape_geometry[2],
                                    cfg.bilagrid_shape_geometry[1],
                                    cfg.bilagrid_shape_geometry[0],
                                    optim_bits, value_bits,
                                    cfg.use_adagrad_bilagrid_optim);
        st.bilagrid_normal_init = true;
    }
    if (cfg.use_ppisp &&
        (cfg.use_adagrad_ppisp_optim ? cfg.ppisp_adagrad_lr : cfg.ppisp_lr) > 0.0f) {
        engine_init_ppisp(n_grids, cfg.ppisp_param_type, cfg.use_adagrad_ppisp_optim);
        st.ppisp_init = true;
    }

    // ---- Resume (resume.py + trainer.py _resume_from_checkpoint:510) --------
    // Last, because engine_load_checkpoint() overwrites the skeleton just
    // built: the world must already be allocated at max_num_splats and every
    // appearance channel the checkpoint carries must already exist as a
    // restore target, which is what everything above establishes.
    if (!cfg.resume.empty()) restore_checkpoint();
}

// Restore engine state from cfg.resume, adapting the checkpoint's buffers on
// the host first when its layout differs from the one just built (fewer
// splats, different SH degree, bilagrid/PPISP added or dropped).
void TrainerSession::restore_checkpoint() {
    ckpt::ResolvedCheckpoint r = ckpt::resolve_checkpoint(cfg.resume);
    ckpt::check_resumable(r.ckpt_dir);

    // The target is what the engine ACTUALLY holds, not what the config asks
    // for: a depth/normal grid also needs the dataset to carry those maps,
    // which setup_engine() resolved into `st`.
    ckpt::TargetLayout target;
    target.max_num_splats = engine_get_max_num_splats();
    target.num_sh         = (cfg.sh_degree + 1) * (cfg.sh_degree + 1);
    target.num_images     = (int)post.n_post;
    auto lhw = [](const std::array<int, 3>& xyw) {
        return std::array<int, 3>{xyw[2], xyw[1], xyw[0]};   // (X,Y,W)->(L,H,W)
    };
    if (st.bilagrid_rgb_init)    target.bilagrid_rgb    = lhw(cfg.bilagrid_shape);
    if (st.bilagrid_depth_init)  target.bilagrid_depth  = lhw(cfg.bilagrid_shape_geometry);
    if (st.bilagrid_normal_init) target.bilagrid_normal = lhw(cfg.bilagrid_shape_geometry);
    target.ppisp = st.ppisp_init;

    fs::path load_from = r.ckpt_dir;
    fs::path tmp;
    bool adapted = false;
    {
        JsonValue state = ckpt::read_state_json(r.ckpt_dir);
        if (ckpt::needs_adapt(state, target)) {
            tmp = out_dir / ".resume_adapt";
            log(lmsg::ckpt_adapting.get());
            adapted = ckpt::adapt_checkpoint(r.ckpt_dir, target, tmp);
            if (adapted) load_from = tmp;
        }
    }

    try {
        start_step = engine_load_checkpoint(load_from.string());
    } catch (const std::exception& e) {
        if (adapted) remove_tree(tmp);
        throw std::runtime_error(
            std::string("cannot resume from ") + r.ckpt_dir.string() + ": " +
            e.what());
    }
    if (adapted) remove_tree(tmp);
    log(lfmt(lmsg::resumed_from, {r.ckpt_dir.string(), (long long)start_step}));
}

// Checkpoint save (trainer.py save_checkpoint:939).
void TrainerSession::save_checkpoint(int step) {
    char name[32];
    std::snprintf(name, sizeof name, "step-%09d.ckpt", step);
    fs::path ckpt = out_dir / name;
    fs::create_directories(ckpt);
    engine_save_checkpoint(ckpt.string(), cfg.save_full_checkpoint, step);
    if (cfg.save_only_latest_checkpoint) {
        std::vector<fs::path> stale;
        for (const auto& e : fs::directory_iterator(out_dir)) {
            std::string b = e.path().filename().string();
            if (b.rfind("step-", 0) == 0 &&
                b.find(".ckpt") != std::string::npos && e.path() != ckpt)
                stale.push_back(e.path());
        }
        for (const auto& p : stale) remove_tree(p);
    }
}

// One step. Split out of train() so a front-end that keeps its own loop
// (the Python trainer, for resume / eval / profiling) shares the ported
// per-step config rather than rebuilding it -- see bindings/bind_trainer.cpp.
std::map<std::string, float> TrainerSession::train_step(int step) {
    int sh_degree_to_use = step / std::max(cfg.sh_degree_warmup_every, 1);
    EngineStepConfig sc = build_step_config(cfg, st, step);
    return engine_train_step_managed(
        step, cfg.num_iterations, cfg.primitive, sh_degree_to_use,
        cfg.packed || cfg.use_bvh, sc);
}

// The training loop (trainer.py train:741).
void TrainerSession::train(const TrainerCallbacks& cb) {
    _start_time = std::chrono::steady_clock::now();

    int step = start_step;
    for (; step < cfg.num_iterations; step++) {
        // Pause gate + render-fairness yield: give viewer render workers an
        // uncontended window to take the engine mutex.
        while (paused.load() && !stop_requested.load())
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        if (stop_requested.load()) break;
        while (render_pending.load())
            std::this_thread::sleep_for(std::chrono::microseconds(500));

        auto step_start = std::chrono::steady_clock::now();
        std::map<std::string, float> losses;
        {
            std::lock_guard<std::mutex> lk(engine_mutex);
            if (step > 0 && cfg.steps_per_save > 0 && step % cfg.steps_per_save == 0)
                save_checkpoint(step);
            losses = train_step(step);
        }
        cur_step = step + 1;
        double latency = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - step_start).count();
        {
            std::lock_guard<std::mutex> lk(_progress_mutex);
            _step_latencies.push_back(latency);
            if (_step_latencies.size() > 100) _step_latencies.pop_front();
        }

        if (cb.on_step) {
            TrainerProgress p;
            p.step = step;
            p.total_steps = cfg.num_iterations;
            p.step_latency = latency;
            p.num_splats = engine_get_cur_num_splats();
            p.losses = std::move(losses);
            cb.on_step(p);
        }
    }

    training_time_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - _start_time).count();
    // Pool capacities are a monotonic high-water mark, so reading them after
    // the loop gives the training-time peak.
    {
        size_t cap = 0;
        for (const auto& e : engine_get_pool_breakdown()) cap += std::get<2>(e);
        engine_vram_mb = (double)(cap + engine_get_scratch_bytes()) / (1024.0 * 1024.0);
    }

    if (cfg.steps_per_save != 0) {
        std::lock_guard<std::mutex> lk(engine_mutex);
        save_checkpoint(step);
        log(lfmt(lmsg::checkpoint_saved, {fs::absolute(out_dir).string()}));
    }
}

std::string TrainerSession::progress_json() {
    int step = cur_step.load();
    double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - _start_time).count();
    double avg = 0.0;
    size_t nlat = 0;
    {
        std::lock_guard<std::mutex> lk(_progress_mutex);
        for (double v : _step_latencies) avg += v;
        nlat = _step_latencies.size();
    }
    if (nlat) avg /= (double)nlat;
    char buf[256];
    if (nlat && step > 0) {
        double eta = (cfg.num_iterations - step) * avg;
        std::snprintf(buf, sizeof buf,
            "{\"step\": %d, \"total_steps\": %d, \"elapsed_time\": %.3f, "
            "\"eta\": %.3f, \"latency_ms\": %.3f, \"paused\": %s}",
            step, cfg.num_iterations, elapsed, eta, avg * 1000.0,
            paused.load() ? "true" : "false");
    } else {
        std::snprintf(buf, sizeof buf,
            "{\"step\": %d, \"total_steps\": %d, \"elapsed_time\": %.3f, "
            "\"eta\": null, \"latency_ms\": null, \"paused\": %s}",
            step, cfg.num_iterations, elapsed,
            paused.load() ? "true" : "false");
    }
    return buf;
}

ViewerRenderConfig TrainerSession::make_viewer_config() const {
    ViewerRenderConfig vc;
    vc.primitive = cfg.primitive;
    vc.packed = cfg.packed || cfg.use_bvh;
    vc.sh_degree_warmup_every = cfg.sh_degree_warmup_every;
    vc.relative_scale = cfg.relative_scale;
    vc.output_median = cfg.mean_median_depth_weight > 0.0f ||
                       cfg.median_depth_normal_reg_weight > 0.0f ||
                       cfg.median_normal_supervision_weight > 0.0f ||
                       cfg.median_render_normal_reg_weight > 0.0f;
    vc.distortion_reg_on = cfg.rgb_distortion_reg != 0.0f ||
                           cfg.depth_distortion_reg != 0.0f ||
                           cfg.normal_distortion_reg != 0.0f;
    vc.train_frame_scale = ds.train_frame_scale;
    vc.train_to_normalized = ds.train_to_normalized;
    vc.base_camera_size = viewer_base_camera_size;
    return vc;
}

ViewerHooks TrainerSession::make_viewer_hooks() {
    ViewerHooks hooks;
    hooks.engine_mutex = &engine_mutex;
    hooks.current_step = [this] { return cur_step.load(); };
    hooks.set_render_pending = [this](bool v) { render_pending = v; };
    hooks.pause_toggle = [this] {
        bool now = !paused.load();
        paused = now;
        return now;
    };
    hooks.progress_json = [this] { return progress_json(); };
    return hooks;
}



// ===========================================================================
// Eval (trainer.py eval:465 + model.py get_image_metrics_and_images:1597)
// ===========================================================================

namespace {

// Clip to [0,1] and quantize to 8-bit, which is what the saved PNGs hold and
// therefore what the Python LPIPS tool sees.
std::vector<uint8_t> to_png_bytes(const std::vector<float>& rgb) {
    std::vector<uint8_t> out(rgb.size());
    for (size_t i = 0; i < rgb.size(); i++) {
        float v = std::min(std::max(rgb[i], 0.0f), 1.0f);
        out[i] = (uint8_t)std::lround(v * 255.0f);
    }
    return out;
}

}  // namespace

void TrainerSession::eval() {
    // eval_mode "all" trains on every frame, so nothing is held out.
    if (cfg.eval_mode == "all") return;

    // Re-parse for the eval side of the split. The parser computes the split
    // over all frames, so this is the exact complement of what training saw.
    DatasetParserConfig pcfg;
    pcfg.recon_dir            = cfg.colmap_recon_dir;
    pcfg.image_dir            = cfg.image_dir;
    pcfg.mask_dir             = cfg.mask_dir;
    pcfg.depth_dir            = cfg.depth_dir;
    pcfg.normal_dir           = cfg.normal_dir;
    pcfg.validation_fraction  = 0.0f;      // no early-stop holdout inside eval
    pcfg.eval_mode            = cfg.eval_mode;
    pcfg.eval_interval        = cfg.eval_interval;
    pcfg.train_split_fraction = cfg.train_split_fraction;
    pcfg.outlier_threshold    = cfg.outlier_threshold;
    pcfg.rescale_camera_to_fit   = cfg.rescale_camera_to_fit;
    pcfg.downscale_rounding_mode = cfg.downscale_rounding_mode;
    pcfg.metashape_xml           = cfg.metashape_xml;
    pcfg.metashape_ply           = cfg.metashape_ply;
    pcfg.metashape_psx           = cfg.metashape_psx;
    pcfg.split                   = "eval";

    ParsedDataset eds = parse_dataset(cfg.data, pcfg, cfg.data_format);
    if (eds.num_cameras == 0) {
        log(lmsg::eval_split_empty.get());
        return;
    }
    // Same world scaling training applied, so the cameras line up with the
    // trained splats.
    if (cfg.relative_scale.has_value()) {
        float rs = *cfg.relative_scale;
        for (int64_t i = 0; i < eds.num_cameras; i++)
            for (int r = 0; r < 3; r++) eds.c2w[i*12 + r*4 + 3] *= rs;
    }
    PostSplitCameras epost = bake_post_split(
        eds, cfg.warp_to_pinhole, cfg.warp_spherical_to_pinhole);

    // One image per step: metrics are per-image, and the batch scheduler would
    // otherwise pack several resolutions into one step.
    DataManagerConfig dm;
    dm.cache_mode  = (cfg.cache_images == "disk") ? CacheMode::DISK : CacheMode::CPU;
    dm.load_masks  = !eds.mask_filenames.empty() || epost.any_fisheye_warp;
    dm.load_depths = false;
    dm.load_normals = false;
    dm.train_batch_size = 1;
    dm.val_batch_size   = 1;
    dm.mask_boundary_offset = cfg.mask_boundary_offset;
    dm.warp_to_pinhole  = cfg.warp_to_pinhole;
    std::vector<int32_t> all_idx((size_t)eds.num_cameras);
    std::iota(all_idx.begin(), all_idx.end(), 0);

    {
        std::lock_guard<std::mutex> lk(engine_mutex);
        engine_setup_data_manager(
            dm, eds.camera_models,
            eds.image_filenames, eds.mask_filenames, {}, {},
            eds.widths, eds.heights,
            epost.any_warp ? epost.K_per_camera : std::vector<int32_t>{},
            epost.any_warp ? epost.post_offsets : std::vector<int32_t>{},
            epost.viewmats, epost.intrins, epost.dist_coeffs,
            epost.input_intrins, epost.input_dist_coeffs,
            all_idx, {});
    }

    log(lfmt(lmsg::eval_views, {(long long)epost.n_post}));

    // Rendering is serial (one process-global engine), but scoring a view --
    // colour correction, the metrics, and the PNG encode -- is pure host work
    // that depends on nothing else, so it runs on a pool while the GPU gets on
    // with the next view. Results are written into a slot indexed by view, so
    // the per-image lists in metrics.json do not depend on who finished first.
    struct ViewJob {
        int64_t index = 0;
        int H = 0, W = 0, C = 0;
        std::vector<float> gt, pred;
    };
    struct ViewScore {
        bool  filled = false;
        float l1 = 0, psnr = 0, ssim = 0, cc_l1 = 0, cc_psnr = 0, cc_ssim = 0;
    };

    // Per worker: a GT + a render + the colour-corrected copy (3 x H*W*C
    // floats) plus SSIM's eleven H*W double buffers. ~725 MB at 4K, so the
    // pool is capped by a memory budget as well as by cores -- an eval run
    // must not be the thing that OOMs the machine.
    const int64_t view_pixels = eds.widths.empty()
        ? (int64_t)1 << 20
        : (int64_t)eds.widths[0] * eds.heights[0];
    const int64_t bytes_per_view = view_pixels * (3 * 3 * 4 + 11 * 8);
    const int64_t kBudget = (int64_t)4 << 30;
    int n_workers = (int)std::max(1u, std::thread::hardware_concurrency() / 4u);
    n_workers = (int)std::min<int64_t>(n_workers,
                                       std::max<int64_t>(1, kBudget / std::max<int64_t>(bytes_per_view, 1)));
    n_workers = std::min<int>(n_workers, 8);

    std::deque<ViewJob> queue;
    std::vector<ViewScore> scores;
    std::mutex qmu, smu;
    std::condition_variable qcv, spacecv;
    bool producing = true;
    std::exception_ptr worker_error;

    auto score_view = [&](ViewJob& j) {
        const int64_t view_px = (int64_t)j.H * j.W * j.C;
        const float* g = j.gt.data();
        thread_local std::vector<float> cc;   // reused across this worker's views
        color_correct_into(j.pred.data(), g, (int64_t)j.H * j.W, j.C, cc);
        ViewScore s;
        s.filled  = true;
        s.l1      = image_l1(g, j.pred.data(), view_px);
        s.psnr    = image_psnr(g, j.pred.data(), view_px);
        s.ssim    = image_ssim(g, j.pred.data(), j.H, j.W, j.C);
        s.cc_l1   = image_l1(g, cc.data(), view_px);
        s.cc_psnr = image_psnr(g, cc.data(), view_px);
        s.cc_ssim = image_ssim(g, cc.data(), j.H, j.W, j.C);
        {
            std::lock_guard<std::mutex> lk(smu);
            if ((size_t)j.index >= scores.size()) scores.resize((size_t)j.index + 1);
            scores[(size_t)j.index] = s;
        }
        if (cfg.save_eval_images) {
            char nm[64];
            auto png = [&](const char* kind, const std::vector<float>& img) {
                std::snprintf(nm, sizeof nm, "eval-%s-%05d.png", kind,
                              (int)j.index);
                std::vector<uint8_t> bytes = to_png_bytes(img);
                stbi_write_png((out_dir / nm).string().c_str(), j.W, j.H, j.C,
                               bytes.data(), j.W * j.C);
            };
            // Both sides, because the LPIPS tool needs the pair -- and the
            // 8-bit clip here is what it will score, so its numbers are
            // reproducible from the files alone.
            png("gt", j.gt);
            png("render", j.pred);
        }
    };

    std::vector<std::thread> workers;
    for (int t = 0; t < n_workers; t++) {
        workers.emplace_back([&] {
            // Threads inside the metrics would oversubscribe against the pool;
            // the per-view work is already the coarser and cheaper split.
#ifdef _OPENMP
            omp_set_num_threads(std::max(1, (int)std::thread::hardware_concurrency() / n_workers));
#endif
            for (;;) {
                ViewJob job;
                {
                    std::unique_lock<std::mutex> lk(qmu);
                    qcv.wait(lk, [&] { return !queue.empty() || !producing; });
                    if (queue.empty()) return;
                    job = std::move(queue.front());
                    queue.pop_front();
                }
                spacecv.notify_one();
                try {
                    score_view(job);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(smu);
                    if (!worker_error) worker_error = std::current_exception();
                    return;
                }
            }
        });
    }

    const int sh_deg = cfg.sh_degree;
    int64_t next_slot = 0;

    for (int64_t i = 0; i < eds.num_cameras; i++) {
        int64_t H = 0, W = 0, B = 0, Cc = 0;
        std::vector<float> gt, render;
        {
            std::lock_guard<std::mutex> lk(engine_mutex);
            int n_view = engine_eval_forward(cfg.primitive, sh_deg, cfg.packed);
            if (n_view == 0) break;
            auto shape = engine_get_render_rgb_shape();
            B = std::get<0>(shape); H = std::get<1>(shape);
            W = std::get<2>(shape); Cc = std::get<3>(shape);
            const int64_t npx = B * H * W * Cc;
            render.resize((size_t)npx);
            gt.resize((size_t)npx);
            engine_copy_render_to_host(
                TorchTensorView{(uint64_t)(uintptr_t)render.data(), 4, {B, H, W, Cc}},
                TorchTensorView{0, 0, {}}, TorchTensorView{0, 0, {}},
                TorchTensorView{0, 0, {}}, TorchTensorView{0, 0, {}});
            engine_copy_gt_rgb_to_host(
                TorchTensorView{(uint64_t)(uintptr_t)gt.data(), 4, {B, H, W, Cc}});
        }

        // Per POST-split view inside the batch (K faces of one input image).
        const int64_t view_px = H * W * Cc;
        for (int64_t v = 0; v < B; v++) {
            ViewJob job;
            job.index = next_slot++;
            job.H = (int)H; job.W = (int)W; job.C = (int)Cc;
            job.gt.assign(gt.begin() + (size_t)(v * view_px),
                          gt.begin() + (size_t)((v + 1) * view_px));
            job.pred.assign(render.begin() + (size_t)(v * view_px),
                            render.begin() + (size_t)((v + 1) * view_px));
            for (float& x : job.pred) x = std::min(std::max(x, 0.0f), 1.0f);
            {
                std::unique_lock<std::mutex> lk(qmu);
                spacecv.wait(lk, [&] { return (int)queue.size() < n_workers; });
                queue.push_back(std::move(job));
            }
            qcv.notify_one();
        }
    }

    {
        std::lock_guard<std::mutex> lk(qmu);
        producing = false;
    }
    qcv.notify_all();
    for (auto& t : workers) t.join();
    if (worker_error) std::rethrow_exception(worker_error);

    std::map<std::string, std::vector<float>> per_image;
    for (const ViewScore& s : scores) {
        if (!s.filled) continue;
        per_image["l1"].push_back(s.l1);
        per_image["psnr"].push_back(s.psnr);
        per_image["ssim"].push_back(s.ssim);
        per_image["cc_l1"].push_back(s.cc_l1);
        per_image["cc_psnr"].push_back(s.cc_psnr);
        per_image["cc_ssim"].push_back(s.cc_ssim);
    }

    if (per_image.empty()) {
        log(lmsg::eval_no_views.get());
        return;
    }

    std::map<std::string, float> avg;
    for (const auto& [k, v] : per_image) {
        double s = 0.0;
        for (float x : v) s += x;
        avg[k] = (float)(s / (double)v.size());
        log("  " + k + ": " + std::to_string(avg[k]));
    }

    // metrics.json: per-image lists plus avg_* scalars, the shape
    // ss_benchmark.py reads back.
    std::ofstream mf((out_dir / "metrics.json").string());
    if (!mf) throw std::runtime_error("cannot write metrics.json");
    mf << "{\n";
    bool first = true;
    for (const auto& [k, v] : per_image) {
        mf << (first ? "" : ",\n") << "    \"" << k << "\": [";
        for (size_t i = 0; i < v.size(); i++)
            mf << (i ? ", " : "") << v[i];
        mf << "]";
        first = false;
    }
    for (const auto& [k, v] : avg)
        mf << ",\n    \"avg_" << k << "\": " << v;
    mf << ",\n    \"num_eval_images\": " << per_image.begin()->second.size();
    // Per-run training stats: a benchmark launches each scene as its own
    // process, so metrics.json is the only place it can read these.
    mf << ",\n    \"training_time\": " << training_time_s;
    mf << ",\n    \"engine_vram\": " << engine_vram_mb;
    mf << "\n}\n";
    log(lfmt(lmsg::eval_metrics_written, {(out_dir / "metrics.json").string()}));
}

}  // namespace spirula
