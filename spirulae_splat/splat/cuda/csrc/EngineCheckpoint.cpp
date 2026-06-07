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
    // When sh_value_bits != 32 the canonical store is the packed buffer
    // (`features_sh_quant{8,16}{,_fpbo}`); decode it on the host into the
    // PLY row directly so we don't pay the ~N*K*12B fp32 staging allocation.
    // The packed copy is ~N*3*K bytes (8-bit) or 2x (16-bit), much smaller
    // than the fp32 staging would be.
    const bool quant_sh_value = (K > 0) && s.world.sh_value_quantize_enabled();
    const int  sh_value_bits  = s.world.sh_value_bits;

    std::vector<float3> h_features_sh;            // fp32 path only
    std::vector<uint8_t> h_value_packed;          // quant path only
    std::vector<float2>  h_value_bounds;          // quant path only
    int64_t  bounds_stride = 256;                 // quant path: cell-block vs fpbo
    if (K > 0 && !quant_sh_value) {
        h_features_sh = _ckpt_d2h<float3>(s.world.features_sh.data_ptr(), (size_t)(N * K));
    } else if (quant_sh_value) {
        // Pick the populated layout. FPBO (per-splat-block) takes precedence
        // since that's what `_ensure_optim_state` allocates when FPBO is on.
        bool fpbo_layout = false;
        const uint8_t* dev_packed = nullptr;
        const float2*  dev_bounds = nullptr;
        int64_t        dev_n_bounds = 0;
        if (sh_value_bits == 8) {
            if (s.world.features_sh_quant8_fpbo.initialized()) {
                dev_packed = s.world.features_sh_quant8_fpbo.packed.data_ptr();
                dev_bounds = s.world.features_sh_quant8_fpbo.bounds.data_ptr();
                dev_n_bounds = s.world.features_sh_quant8_fpbo.n_bounds;
                fpbo_layout = true;
            } else if (s.world.features_sh_quant8.initialized()) {
                dev_packed = s.world.features_sh_quant8.packed.data_ptr();
                dev_bounds = s.world.features_sh_quant8.bounds.data_ptr();
                dev_n_bounds = s.world.features_sh_quant8.n_bounds;
            }
        } else /* sh_value_bits == 16 */ {
            if (s.world.features_sh_quant16_fpbo.initialized()) {
                dev_packed = s.world.features_sh_quant16_fpbo.packed.data_ptr();
                dev_bounds = s.world.features_sh_quant16_fpbo.bounds.data_ptr();
                dev_n_bounds = s.world.features_sh_quant16_fpbo.n_bounds;
                fpbo_layout = true;
            } else if (s.world.features_sh_quant16.initialized()) {
                dev_packed = s.world.features_sh_quant16.packed.data_ptr();
                dev_bounds = s.world.features_sh_quant16.bounds.data_ptr();
                dev_n_bounds = s.world.features_sh_quant16.n_bounds;
            }
        }
        if (dev_packed != nullptr && dev_bounds != nullptr) {
            // Copy only the bytes needed for the first N splats (skip the
            // tail of the over-allocated max_N buffer).
            const int64_t cells_needed = (int64_t)N * 3 * K;
            const int bytes_per_cell = (sh_value_bits == 8) ? 1 : 2;
            h_value_packed.resize((size_t)(cells_needed * bytes_per_cell));
            cudaMemcpy(h_value_packed.data(), dev_packed,
                       (size_t)(cells_needed * bytes_per_cell),
                       cudaMemcpyDeviceToHost);
            // Bounds layout: cell-block bounds_stride = 256;
            //                fpbo (per-splat-block) bounds_stride = 256 * 3 * K.
            int64_t n_bounds_needed;
            if (fpbo_layout) {
                bounds_stride = (int64_t)256 * 3 * (int64_t)K;
                n_bounds_needed = (N + 255) / 256;
            } else {
                bounds_stride = 256;
                n_bounds_needed = (cells_needed + 255) / 256;
            }
            // Don't run off the allocated bounds array.
            n_bounds_needed = std::min(n_bounds_needed, dev_n_bounds);
            h_value_bounds.resize((size_t)n_bounds_needed);
            cudaMemcpy(h_value_bounds.data(), dev_bounds,
                       (size_t)n_bounds_needed * sizeof(float2),
                       cudaMemcpyDeviceToHost);
        }
    }

    // Decode one SH cell value on the host. Inlined into the PLY write loop
    // below; safe no-op (returns 0) if quant buffers weren't populated.
    auto decode_sh = [&](int64_t splat_i, int j, int ch) -> float {
        if (h_value_packed.empty() || h_value_bounds.empty()) return 0.0f;
        const int64_t cell = splat_i * (int64_t)3 * (int64_t)K + (int64_t)3 * j + (int64_t)ch;
        const int64_t b_idx = cell / bounds_stride;
        if (b_idx >= (int64_t)h_value_bounds.size()) return 0.0f;
        const float2 mm = h_value_bounds[(size_t)b_idx];
        if (sh_value_bits == 8) {
            const uint8_t q = h_value_packed[(size_t)cell];
            return mm.x + (mm.y - mm.x) * ((float)q * (1.0f / 255.0f));
        } else {
            const uint16_t* p16 = reinterpret_cast<const uint16_t*>(h_value_packed.data());
            const uint16_t q = p16[(size_t)cell];
            return mm.x + (mm.y - mm.x) * ((float)q * (1.0f / 65535.0f));
        }
    };

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
        // Quant-decoded SH values are finite by construction (byte to fp32
        // linear interp from finite bounds), so skip the NaN check on that
        // path. The fp32 path still validates each coef.
        if (ok && K > 0 && !quant_sh_value) {
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
            // SH rest in (3, K) order: r0..r_{K-1}, g0..g_{K-1}, b0..b_{K-1}.
            // Branch hoisted out of the inner loops so the codec decode only
            // runs when value-quant is active.
            if (quant_sh_value) {
                for (int j = 0; j < K; j++) row[p++] = decode_sh(i, j, 0);
                for (int j = 0; j < K; j++) row[p++] = decode_sh(i, j, 1);
                for (int j = 0; j < K; j++) row[p++] = decode_sh(i, j, 2);
            } else {
                for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].x;
                for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].y;
                for (int j = 0; j < K; j++) row[p++] = h_features_sh[i*K + j].z;
            }
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
        // Bit depth (32 = no quant, 4 or 8 = packed). Replaces the old
        // train_quantize_sh bool. Offline tools should treat absence of this
        // key as "8 if train_quantize_sh==1 else 32" for back-compat.
        meta << "sh_optim_bits=" << s.optim.sh_optim_bits << "\n";
        // SH VALUE-quant bit depth (32 = off, 8 or 16 = packed). When != 32,
        // world_features_sh.npy is OMITTED and replaced by world_features_sh_packed.npy
        // (uint8 raw, kBytesPerCell = 1 for 8-bit, 2 for 16-bit) plus
        // world_features_sh_bounds.npy (float2 per block). The bound layout
        // (cell-block vs per-splat-block) is signaled by world_features_sh_layout.
        meta << "sh_value_bits=" << s.world.sh_value_bits << "\n";
        // QuantizedAdamState codec signature -- tells offline tools which
        // byte-to-(g1, g2) inverse to use when reading the *_packed.npy
        // streams. Codec scheme: joint (u, log_s) AoS with u linear and
        // log_s = log1p(sqrt_g2 / eps), eps = 1e-15. The PER-BUFFER bit
        // depth is given by sh_optim_bits / bilagrid_*_optim_bits below.
        meta << "quant_codec=joint_u_log_sqrt_g2_v1\n";
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
            meta << "bilagrid_rgb_optim_bits=" << s.bilagrid_rgb.optim_bits << "\n";
            meta << "bilagrid_rgb_use_adagrad=" << (s.bilagrid_rgb.use_adagrad ? 1 : 0) << "\n";
        }
        if (s.bilagrid_depth.enabled) {
            meta << "bilagrid_depth_enabled=1\n";
            meta << "bilagrid_depth_L=" << s.bilagrid_depth.grids.size<1>() << "\n";
            meta << "bilagrid_depth_H=" << s.bilagrid_depth.grids.size<2>() << "\n";
            meta << "bilagrid_depth_W=" << s.bilagrid_depth.grids.size<3>() << "\n";
            meta << "bilagrid_depth_optim_bits=" << s.bilagrid_depth.optim_bits << "\n";
            meta << "bilagrid_depth_use_adagrad=" << (s.bilagrid_depth.use_adagrad ? 1 : 0) << "\n";
        }
        if (s.bilagrid_normal.enabled) {
            meta << "bilagrid_normal_enabled=1\n";
            meta << "bilagrid_normal_L=" << s.bilagrid_normal.grids.size<1>() << "\n";
            meta << "bilagrid_normal_H=" << s.bilagrid_normal.grids.size<2>() << "\n";
            meta << "bilagrid_normal_W=" << s.bilagrid_normal.grids.size<3>() << "\n";
            meta << "bilagrid_normal_optim_bits=" << s.bilagrid_normal.optim_bits << "\n";
            meta << "bilagrid_normal_use_adagrad=" << (s.bilagrid_normal.use_adagrad ? 1 : 0) << "\n";
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
    if (!s.world.sh_value_quantize_enabled()) {
        save_dt2d_f3("world_features_sh.npy", s.world.features_sh);
    } else {
        // SH VALUE-quant active: canonical storage is packed bytes + per-block
        // bounds. Only the populated layout (cell-block vs FPBO) gets written;
        // the unused half is an empty (default-constructed) tensor and the
        // saver skips zero-size buffers.
        auto save_quant_pair = [&](const char* layout_tag,
                                   const auto& q8, const auto& q16) {
            // Helper to dump uint8_t* packed.
            auto save_u8 = [&](const char* name, const uint8_t* data, size_t n) {
                if (!data || n == 0) return;
                _ckpt_save_npy<uint8_t>(full_root / name, data, n,
                                        {(unsigned long)n});
            };
            auto save_f2 = [&](const char* name, const float2* data, size_t n) {
                if (!data || n == 0) return;
                _ckpt_save_npy<float>(full_root / name,
                                      reinterpret_cast<const float*>(data),
                                      n * 2,
                                      {(unsigned long)n, 2ul});
            };
            if (s.world.sh_value_bits == 8 && q8.initialized()) {
                save_u8("world_features_sh_packed.npy",
                        q8.packed.data_ptr(), (size_t)q8.packed_bytes());
                save_f2("world_features_sh_bounds.npy",
                        q8.bounds.data_ptr(), (size_t)q8.n_bounds);
            } else if (s.world.sh_value_bits == 16 && q16.initialized()) {
                save_u8("world_features_sh_packed.npy",
                        q16.packed.data_ptr(), (size_t)q16.packed_bytes());
                save_f2("world_features_sh_bounds.npy",
                        q16.bounds.data_ptr(), (size_t)q16.n_bounds);
            }
            // Layout marker so offline tools know how to index bounds:
            //   cell  -> bounds[cell_idx / 256]
            //   fpbo  -> bounds[splat_idx / 256]
            std::ofstream layout(full_root / "world_features_sh_layout.txt");
            layout << layout_tag << "\n";
        };
        // Use whichever layout actually has bytes; fall back to FPBO if both
        // empty (shouldn't happen, but harmless).
        bool cell_block_active = s.world.features_sh_quant8.initialized()
                              || s.world.features_sh_quant16.initialized();
        if (cell_block_active) {
            save_quant_pair("cell", s.world.features_sh_quant8, s.world.features_sh_quant16);
        } else {
            save_quant_pair("fpbo", s.world.features_sh_quant8_fpbo, s.world.features_sh_quant16_fpbo);
        }
    }

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
    if (s.optim.sh_quantize_enabled()) {
        // Joint (u, sqrt_g2) AoS: single packed byte stream + per-block bounds.
        // FPBO mode uses a different per-bound granularity, so it lives in a
        // separate state buffer.
        if (s.optim.sh_quant_state.packed.data_ptr())
            save_dv_u8("sh_quant_packed.npy", s.optim.sh_quant_state.packed);
        if (s.optim.sh_quant_state.bounds.data_ptr())
            save_dv_f4("sh_quant_bounds.npy", s.optim.sh_quant_state.bounds);
        if (s.optim.sh_quant_state_fpbo.packed.data_ptr())
            save_dv_u8("sh_quant_packed_fpbo.npy", s.optim.sh_quant_state_fpbo.packed);
        if (s.optim.sh_quant_state_fpbo.bounds.data_ptr())
            save_dv_f4("sh_quant_bounds_fpbo.npy", s.optim.sh_quant_state_fpbo.bounds);
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

    // Bilagrid grids + optimizer state, per type. Three optimizer layouts:
    //   Adam float        : *_g1.npy, *_g2.npy
    //   Adam quantized    : *_packed.npy (AoS u,sqrt_g2) + *_quant_bounds.npy
    //   AdaGrad float     : *_accum.npy
    //   AdaGrad quantized : *_accum_packed.npy + *_accum_bounds.npy (float2)
    auto save_dv_f2 = [&](const char* name, const DeviceVector<float2>& v) {
        if (!v.data_ptr() || v.size() == 0) return;
        _ckpt_save_npy<float>(full_root / name,
                              reinterpret_cast<const float*>(v.data_ptr()),
                              (size_t)v.size() * 2,
                              {(unsigned long)v.size(), 2ul});
    };
    auto save_bilagrid = [&](const char* prefix, bool enabled, bool quantize,
                             bool use_adagrad,
                             const DeviceTensor5D<float>& grids,
                             const DeviceTensor5D<float>& g1_f,
                             const DeviceTensor5D<float>& g2_f,
                             const QuantizedAdamState<8, 256>& quant_state,
                             const DeviceTensor5D<float>& accum_f,
                             const QuantizedTensorLog<8, 256>& adagrad_quant) {
        if (!enabled) return;
        save_5d((std::string(prefix) + ".npy").c_str(), grids);
        if (use_adagrad) {
            if (quantize) {
                save_dv_u8((std::string(prefix) + "_accum_packed.npy").c_str(),
                           adagrad_quant.packed);
                save_dv_f2((std::string(prefix) + "_accum_bounds.npy").c_str(),
                           adagrad_quant.bounds);
            } else {
                save_5d((std::string(prefix) + "_accum.npy").c_str(), accum_f);
            }
        } else {
            if (quantize) {
                save_dv_u8((std::string(prefix) + "_packed.npy").c_str(),
                           quant_state.packed);
                save_dv_f4((std::string(prefix) + "_quant_bounds.npy").c_str(),
                           quant_state.bounds);
            } else {
                save_5d((std::string(prefix) + "_g1.npy").c_str(), g1_f);
                save_5d((std::string(prefix) + "_g2.npy").c_str(), g2_f);
            }
        }
    };
    save_bilagrid("bilagrid_rgb",
                  s.bilagrid_rgb.enabled, s.bilagrid_rgb.quantize_optim(),
                  s.bilagrid_rgb.use_adagrad,
                  s.bilagrid_rgb.grids, s.bilagrid_rgb.g1, s.bilagrid_rgb.g2,
                  s.bilagrid_rgb.quant_state,
                  s.bilagrid_rgb.accum_f, s.bilagrid_rgb.adagrad_quant);
    save_bilagrid("bilagrid_depth",
                  s.bilagrid_depth.enabled, s.bilagrid_depth.quantize_optim(),
                  s.bilagrid_depth.use_adagrad,
                  s.bilagrid_depth.grids, s.bilagrid_depth.g1, s.bilagrid_depth.g2,
                  s.bilagrid_depth.quant_state,
                  s.bilagrid_depth.accum_f, s.bilagrid_depth.adagrad_quant);
    if (s.bilagrid_depth.enabled) {
        save_dv_f1("bilagrid_depth_scalars.npy", s.bilagrid_depth.scalars);
    }
    save_bilagrid("bilagrid_normal",
                  s.bilagrid_normal.enabled, s.bilagrid_normal.quantize_optim(),
                  s.bilagrid_normal.use_adagrad,
                  s.bilagrid_normal.grids, s.bilagrid_normal.g1, s.bilagrid_normal.g2,
                  s.bilagrid_normal.quant_state,
                  s.bilagrid_normal.accum_f, s.bilagrid_normal.adagrad_quant);

    // PPISP table + Adam moments.
    if (s.ppisp.enabled) {
        save_dt2d_f("ppisp.npy",    s.ppisp.params);
        save_dt2d_f("ppisp_g1.npy", s.ppisp.g1);
        save_dt2d_f("ppisp_g2.npy", s.ppisp.g2);
    }
}
