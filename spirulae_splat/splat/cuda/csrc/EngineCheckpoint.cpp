// Engine checkpoint save (PLY + npy + meta.txt). Optional full/ dump of every
// persistent buffer for offline inspection.

#include "Engine.h"
#include "EngineState.h"

#include "npy.hpp"
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>


namespace {

namespace fs = std::filesystem;

static void _ckpt_mkdir(const fs::path& p) {
    std::error_code ec;
    fs::create_directories(p, ec);
    if (ec) throw std::runtime_error("create_directories failed: " + p.string() + " (" + ec.message() + ")");
}

// D->H copy into a freshly sized host buffer.
template<typename T>
static std::vector<T> _ckpt_d2h(const T* d_ptr, size_t n) {
    std::vector<T> host(n);
    if (n > 0 && d_ptr != nullptr) {
        cudaMemcpy(host.data(), d_ptr, n * sizeof(T), cudaMemcpyDeviceToHost);
    }
    return host;
}

// Save a device-side buffer of `numel` elements as a .npy with the given shape.
template<typename Scalar>
static void _ckpt_save_npy(const fs::path& path, const Scalar* d_ptr, size_t numel,
                           const npy::shape_t& shape) {
    if (d_ptr == nullptr || numel == 0) return;
    auto host = _ckpt_d2h<Scalar>(d_ptr, numel);
    npy::npy_data_ptr<Scalar> data{host.data(), shape, false};
    npy::write_npy(path.string(), data);
}

template<typename T>
static npy::shape_t _ckpt_shape_5d(const DeviceTensor5D<T>& t) {
    return {
        (unsigned long)t.template size<0>(),
        (unsigned long)t.template size<1>(),
        (unsigned long)t.template size<2>(),
        (unsigned long)t.template size<3>(),
        (unsigned long)t.template size<4>(),
    };
}

} // anon namespace


void engine_save_checkpoint(
    std::string output_dir,
    bool full_dump,
    int step
) {
    EngineState& s = engine();

    fs::path out_root(output_dir);
    _ckpt_mkdir(out_root);

    cudaDeviceSynchronize();

    const int64_t N = s.cur_num_splats;
    const int K = s.num_sh;

    // --- D->H copy world splat tensors at s.cur_num_splats ---
    auto h_means        = _ckpt_d2h<float3>(s.world.means.data_ptr(),        (size_t)N);
    auto h_quats        = _ckpt_d2h<float4>(s.world.quats.data_ptr(),        (size_t)N);
    auto h_scales       = _ckpt_d2h<float3>(s.world.scales.data_ptr(),       (size_t)N);
    auto h_opacities    = _ckpt_d2h<float >(s.world.opacities.data_ptr(),    (size_t)N);
    auto h_features_dc  = _ckpt_d2h<float3>(s.world.features_dc.data_ptr(),  (size_t)N);
    // features_sh is stored [max_N, K] contiguous; copying first N rows = N*K float3.
    std::vector<float3> h_features_sh;
    if (K > 0) {
        h_features_sh = _ckpt_d2h<float3>(s.world.features_sh.data_ptr(), (size_t)(N * K));
    }

    // --- Filter mask: NaN/Inf + low-opacity (logit < logit(1/255) ~= -5.5373) ---
    const float OPA_MIN = -5.5373f;
    auto fin1 = [](float v) { return std::isfinite(v); };
    auto fin3 = [&](const float3& v) { return fin1(v.x) && fin1(v.y) && fin1(v.z); };
    auto fin4 = [&](const float4& v) { return fin1(v.x) && fin1(v.y) && fin1(v.z) && fin1(v.w); };

    std::vector<char> keep((size_t)N, 1);
    int64_t kept = 0, nan_dropped = 0, lowopa_dropped = 0;
    for (int64_t i = 0; i < N; i++) {
        bool ok = fin3(h_means[i]) && fin4(h_quats[i]) && fin3(h_scales[i])
                  && fin1(h_opacities[i]) && fin3(h_features_dc[i]);
        if (ok && K > 0) {
            for (int j = 0; j < K; j++) {
                if (!fin3(h_features_sh[i*K + j])) { ok = false; break; }
            }
        }
        if (!ok) { keep[i] = 0; nan_dropped++; continue; }
        if (h_opacities[i] < OPA_MIN) { keep[i] = 0; lowopa_dropped++; continue; }
        kept++;
    }

    // --- PLY (binary little-endian) ---
    const int vertex_floats = 3 + 3 + 3 + 3 * K + 1 + 3 + 4;  // pos + normal + dc + sh + opa + scale + rot
    fs::path ply_path = out_root / "splat.ply";
    {
        std::ofstream ply(ply_path, std::ios::binary);
        if (!ply) throw std::runtime_error("Failed to open PLY for write: " + ply_path.string());
        ply << "ply\n";
        ply << "format binary_little_endian 1.0\n";
        ply << "element vertex " << kept << "\n";
        ply << "property float x\nproperty float y\nproperty float z\n";
        ply << "property float nx\nproperty float ny\nproperty float nz\n";
        ply << "property float f_dc_0\nproperty float f_dc_1\nproperty float f_dc_2\n";
        for (int j = 0; j < 3 * K; j++) ply << "property float f_rest_" << j << "\n";
        ply << "property float opacity\n";
        ply << "property float scale_0\nproperty float scale_1\nproperty float scale_2\n";
        ply << "property float rot_0\nproperty float rot_1\nproperty float rot_2\nproperty float rot_3\n";
        ply << "end_header\n";

        // Buffer multiple rows per write for throughput.
        constexpr int ROWS_PER_FLUSH = 1024;
        std::vector<float> buf((size_t)vertex_floats * ROWS_PER_FLUSH);
        int rows_in_buf = 0;
        auto flush = [&]() {
            if (rows_in_buf == 0) return;
            ply.write(reinterpret_cast<const char*>(buf.data()),
                      (std::streamsize)rows_in_buf * vertex_floats * sizeof(float));
            rows_in_buf = 0;
        };
        for (int64_t i = 0; i < N; i++) {
            if (!keep[i]) continue;
            float* row = buf.data() + (size_t)rows_in_buf * vertex_floats;
            int p = 0;
            row[p++] = h_means[i].x; row[p++] = h_means[i].y; row[p++] = h_means[i].z;
            row[p++] = 0.f; row[p++] = 0.f; row[p++] = 0.f;  // nx, ny, nz
            row[p++] = h_features_dc[i].x; row[p++] = h_features_dc[i].y; row[p++] = h_features_dc[i].z;
            // SH rest in (3, K) order: r0..r_{K-1}, g0..g_{K-1}, b0..b_{K-1}
            for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].x;
            for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].y;
            for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].z;
            row[p++] = h_opacities[i];
            row[p++] = h_scales[i].x; row[p++] = h_scales[i].y; row[p++] = h_scales[i].z;
            row[p++] = h_quats[i].x; row[p++] = h_quats[i].y; row[p++] = h_quats[i].z; row[p++] = h_quats[i].w;
            rows_in_buf++;
            if (rows_in_buf == ROWS_PER_FLUSH) flush();
        }
        flush();
    }

    // --- Top-level bilagrid / PPISP grids (param tables only, no optimizer state) ---
    if (s.bilagrid_rgb.enabled) {
        _ckpt_save_npy<float>(out_root / "bilagrid_rgb.npy",
                              s.bilagrid_rgb.grids.data_ptr(),
                              (size_t)s.bilagrid_rgb.grids.numel(),
                              _ckpt_shape_5d(s.bilagrid_rgb.grids));
    }
    if (s.bilagrid_depth.enabled) {
        _ckpt_save_npy<float>(out_root / "bilagrid_depth.npy",
                              s.bilagrid_depth.grids.data_ptr(),
                              (size_t)s.bilagrid_depth.grids.numel(),
                              _ckpt_shape_5d(s.bilagrid_depth.grids));
        if (s.bilagrid_depth.scalars.data_ptr() && s.bilagrid_depth.scalars.size() > 0) {
            _ckpt_save_npy<float>(out_root / "bilagrid_depth_scalars.npy",
                                  s.bilagrid_depth.scalars.data_ptr(),
                                  (size_t)s.bilagrid_depth.scalars.size(),
                                  {(unsigned long)s.bilagrid_depth.scalars.size()});
        }
    }
    if (s.bilagrid_normal.enabled) {
        _ckpt_save_npy<float>(out_root / "bilagrid_normal.npy",
                              s.bilagrid_normal.grids.data_ptr(),
                              (size_t)s.bilagrid_normal.grids.numel(),
                              _ckpt_shape_5d(s.bilagrid_normal.grids));
    }
    if (s.ppisp.enabled) {
        _ckpt_save_npy<float>(out_root / "ppisp.npy",
                              s.ppisp.params.data_ptr(),
                              (size_t)(s.ppisp.params.size<0>() * s.ppisp.params.size<1>()),
                              {(unsigned long)s.ppisp.params.size<0>(),
                               (unsigned long)s.ppisp.params.size<1>()});
    }

    // --- meta.txt ---
    {
        std::ofstream meta((out_root / "meta.txt").string());
        meta << "step=" << step << "\n";
        meta << "primitive=" << s.primitive << "\n";
        meta << "cur_num_splats=" << s.cur_num_splats << "\n";
        meta << "max_num_splats=" << s.max_num_splats << "\n";
        meta << "num_sh=" << s.num_sh << "\n";
        meta << "sh_degree=" << s.sh_degree << "\n";
        meta << "num_cameras=" << s.camera.num << "\n";
        meta << "image_width=" << s.camera.width << "\n";
        meta << "image_height=" << s.camera.height << "\n";
        meta << "camera_model=" << s.camera.model_str << "\n";
        meta << "train_quantize_sh=" << (s.optim.quantize_sh ? 1 : 0) << "\n";
        // QuantizedAdamState codec signature — tells offline tools which
        // byte-to-(g1, g2) inverse to use when reading the *_packed.npy
        // streams. Currently every quantized buffer (SH + bilagrid) uses the
        // same scheme: 8-bit AoS, 256-cell blocks, joint (u, log_s) primitives
        // with u linear and log_s = log1p(sqrt_g2 / eps), eps = 1e-15.
        meta << "quant_codec=joint_u_log_sqrt_g2_v1\n";
        meta << "quant_bits=8\n";
        meta << "quant_block_size=256\n";
        meta << "quant_eps=1e-15\n";
        meta << "use_per_splat_bias_correction=" << (s.optim.use_per_splat_bias_correction ? 1 : 0) << "\n";
        meta << "ply_kept=" << kept << "\n";
        meta << "ply_nan_dropped=" << nan_dropped << "\n";
        meta << "ply_low_opacity_dropped=" << lowopa_dropped << "\n";
        if (s.bilagrid_rgb.enabled) {
            meta << "bilagrid_rgb_enabled=1\n";
            meta << "bilagrid_rgb_type=" << s.bilagrid_rgb.type << "\n";
            meta << "bilagrid_rgb_C=" << s.bilagrid_rgb.C << "\n";
            meta << "bilagrid_rgb_L=" << s.bilagrid_rgb.grids.size<1>() << "\n";
            meta << "bilagrid_rgb_H=" << s.bilagrid_rgb.grids.size<2>() << "\n";
            meta << "bilagrid_rgb_W=" << s.bilagrid_rgb.grids.size<3>() << "\n";
            meta << "bilagrid_rgb_quantize_optim=" << (s.bilagrid_rgb.quantize_optim ? 1 : 0) << "\n";
        }
        if (s.bilagrid_depth.enabled) {
            meta << "bilagrid_depth_enabled=1\n";
            meta << "bilagrid_depth_L=" << s.bilagrid_depth.grids.size<1>() << "\n";
            meta << "bilagrid_depth_H=" << s.bilagrid_depth.grids.size<2>() << "\n";
            meta << "bilagrid_depth_W=" << s.bilagrid_depth.grids.size<3>() << "\n";
            meta << "bilagrid_depth_quantize_optim=" << (s.bilagrid_depth.quantize_optim ? 1 : 0) << "\n";
        }
        if (s.bilagrid_normal.enabled) {
            meta << "bilagrid_normal_enabled=1\n";
            meta << "bilagrid_normal_L=" << s.bilagrid_normal.grids.size<1>() << "\n";
            meta << "bilagrid_normal_H=" << s.bilagrid_normal.grids.size<2>() << "\n";
            meta << "bilagrid_normal_W=" << s.bilagrid_normal.grids.size<3>() << "\n";
            meta << "bilagrid_normal_quantize_optim=" << (s.bilagrid_normal.quantize_optim ? 1 : 0) << "\n";
        }
        if (s.ppisp.enabled) {
            meta << "ppisp_enabled=1\n";
            meta << "ppisp_param_type=" << s.ppisp.param_type << "\n";
            meta << "ppisp_num_params=" << s.ppisp.num_params << "\n";
        }
    }

    if (!full_dump) return;

    // --- full/ subfolder: every persistent buffer, one .npy per buffer ---
    fs::path full_root = out_root / "full";
    _ckpt_mkdir(full_root);

    // World splat parameters at s.max_num_splats (full table).
    auto save_dv_f3 = [&](const char* name, const DeviceVector<float3>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<float>(full_root / name,
                              reinterpret_cast<const float*>(v.data_ptr()),
                              (size_t)v.size() * 3,
                              {(unsigned long)v.size(), 3ul});
    };
    auto save_dv_f4 = [&](const char* name, const DeviceVector<float4>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<float>(full_root / name,
                              reinterpret_cast<const float*>(v.data_ptr()),
                              (size_t)v.size() * 4,
                              {(unsigned long)v.size(), 4ul});
    };
    auto save_dv_f1 = [&](const char* name, const DeviceVector<float>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<float>(full_root / name, v.data_ptr(),
                              (size_t)v.size(), {(unsigned long)v.size()});
    };
    auto save_dt2d_f3 = [&](const char* name, const DeviceTensor2D<float3>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<float>(full_root / name,
                              reinterpret_cast<const float*>(t.data_ptr()),
                              (size_t)t.size<0>() * t.size<1>() * 3,
                              {(unsigned long)t.size<0>(), (unsigned long)t.size<1>(), 3ul});
    };
    auto save_dt2d_uc3 = [&](const char* name, const DeviceTensor2D<uchar3>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<uint8_t>(full_root / name,
                                reinterpret_cast<const uint8_t*>(t.data_ptr()),
                                (size_t)t.size<0>() * t.size<1>() * 3,
                                {(unsigned long)t.size<0>(), (unsigned long)t.size<1>(), 3ul});
    };
    auto save_dt2d_f = [&](const char* name, const DeviceTensor2D<float>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<float>(full_root / name, t.data_ptr(),
                              (size_t)t.size<0>() * t.size<1>(),
                              {(unsigned long)t.size<0>(), (unsigned long)t.size<1>()});
    };
    auto save_5d = [&](const char* name, const DeviceTensor5D<float>& t) {
        if (!t.data_ptr()) return;
        _ckpt_save_npy<float>(full_root / name, t.data_ptr(),
                              (size_t)t.numel(), _ckpt_shape_5d(t));
    };
    auto save_dv_u8 = [&](const char* name, const DeviceVector<uint8_t>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<uint8_t>(full_root / name, v.data_ptr(),
                                (size_t)v.size(), {(unsigned long)v.size()});
    };

    // World params (s.max_num_splats sized).
    save_dv_f3("world_means.npy",        s.world.means);
    save_dv_f4("world_quats.npy",        s.world.quats);
    save_dv_f3("world_scales.npy",       s.world.scales);
    save_dv_f1("world_opacities.npy",    s.world.opacities);
    save_dv_f3("world_features_dc.npy",  s.world.features_dc);
    save_dt2d_f3("world_features_sh.npy", s.world.features_sh);

    // World Adam optimizer state.
    save_dv_f3("g1_means.npy",        s.optim.g1_means);
    save_dv_f3("g2_means.npy",        s.optim.g2_means);
    save_dv_f4("g1_quats.npy",        s.optim.g1_quats);
    save_dv_f4("g2_quats.npy",        s.optim.g2_quats);
    save_dv_f3("g1_scales.npy",       s.optim.g1_scales);
    save_dv_f3("g2_scales.npy",       s.optim.g2_scales);
    save_dv_f1("g1_opacities.npy",    s.optim.g1_opacities);
    save_dv_f1("g2_opacities.npy",    s.optim.g2_opacities);
    save_dv_f3("g1_features_dc.npy",  s.optim.g1_features_dc);
    save_dv_f3("g2_features_dc.npy",  s.optim.g2_features_dc);
    if (s.optim.quantize_sh) {
        // Joint (u, sqrt_g2) AoS: single packed byte stream + per-block bounds.
        if (s.optim.sh_quant_state.packed.data_ptr())
            save_dv_u8("sh_quant_packed.npy", s.optim.sh_quant_state.packed);
        if (s.optim.sh_quant_state.bounds.data_ptr())
            save_dv_f4("sh_quant_bounds.npy", s.optim.sh_quant_state.bounds);
    } else {
        save_dt2d_f3("g1_features_sh.npy", s.optim.g1_features_sh);
        save_dt2d_f3("g2_features_sh.npy", s.optim.g2_features_sh);
    }

    // Per-splat aux.
    if (s.optim.use_per_splat_bias_correction
        && s.optim.bias_correction_steps.data_ptr() && s.optim.bias_correction_steps.size() > 0) {
        _ckpt_save_npy<int32_t>(full_root / "bias_correction_steps.npy",
                                s.optim.bias_correction_steps.data_ptr(),
                                (size_t)s.optim.bias_correction_steps.size(),
                                {(unsigned long)s.optim.bias_correction_steps.size()});
    }
    if (s.optim.accum_buffer.data_ptr() && s.optim.accum_buffer.size() > 0) {
        _ckpt_save_npy<float>(full_root / "accum_buffer.npy",
                              reinterpret_cast<const float*>(s.optim.accum_buffer.data_ptr()),
                              (size_t)s.optim.accum_buffer.size() * 2,
                              {(unsigned long)s.optim.accum_buffer.size(), 2ul});
    }

    // Bilagrid grids + optimizer state, per type.
    auto save_bilagrid = [&](const char* prefix, bool enabled, bool quantize,
                             const DeviceTensor5D<float>& grids,
                             const DeviceTensor5D<float>& g1_f,
                             const DeviceTensor5D<float>& g2_f,
                             const QuantizedAdamState<8, 256>& quant_state) {
        if (!enabled) return;
        save_5d((std::string(prefix) + ".npy").c_str(), grids);
        if (quantize) {
            // AoS packed (u, sqrt_g2) byte stream + per-block bounds.
            save_dv_u8((std::string(prefix) + "_packed.npy").c_str(), quant_state.packed);
            save_dv_f4((std::string(prefix) + "_quant_bounds.npy").c_str(), quant_state.bounds);
        } else {
            save_5d((std::string(prefix) + "_g1.npy").c_str(), g1_f);
            save_5d((std::string(prefix) + "_g2.npy").c_str(), g2_f);
        }
    };
    save_bilagrid("bilagrid_rgb",   s.bilagrid_rgb.enabled,   s.bilagrid_rgb.quantize_optim,
                  s.bilagrid_rgb.grids, s.bilagrid_rgb.g1, s.bilagrid_rgb.g2,
                  s.bilagrid_rgb.quant_state);
    save_bilagrid("bilagrid_depth", s.bilagrid_depth.enabled, s.bilagrid_depth.quantize_optim,
                  s.bilagrid_depth.grids, s.bilagrid_depth.g1, s.bilagrid_depth.g2,
                  s.bilagrid_depth.quant_state);
    if (s.bilagrid_depth.enabled) {
        save_dv_f1("bilagrid_depth_scalars.npy", s.bilagrid_depth.scalars);
    }
    save_bilagrid("bilagrid_normal", s.bilagrid_normal.enabled, s.bilagrid_normal.quantize_optim,
                  s.bilagrid_normal.grids, s.bilagrid_normal.g1, s.bilagrid_normal.g2,
                  s.bilagrid_normal.quant_state);

    // PPISP table + Adam moments.
    if (s.ppisp.enabled) {
        save_dt2d_f("ppisp.npy",    s.ppisp.params);
        save_dt2d_f("ppisp_g1.npy", s.ppisp.g1);
        save_dt2d_f("ppisp_g2.npy", s.ppisp.g2);
    }
}
