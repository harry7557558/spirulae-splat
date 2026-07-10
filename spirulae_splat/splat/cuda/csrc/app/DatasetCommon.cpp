// DatasetCommon.cpp -- shared dataset-bake helpers (see DatasetParser.h)
// plus the POST-split camera expansion (bake_post_split), used by both the
// COLMAP and nerfstudio parsers and by main.cpp.

#include "DatasetParser.h"

#include "../Camera.h"   // camera_model_from_name

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace fs = std::filesystem;


namespace dsparse {

// ---------------------------------------------------------------------------
// Normalized-frame scale (port of the scale-relevant subset of
// camera_utils.auto_orient_and_center_poses + auto_scale, dataparser.py:461):
//   up      = normalize(mean of c2w Y columns)
//   R_align = rotation taking `up` to +Z (Rodrigues)
//   center  = mean camera position
//   scale_factor = 1 / max |R_align @ (pos - center)|
// TODO: "pca" / "vertical" / "gsplat" / "focus" methods.
// ---------------------------------------------------------------------------
double compute_normalized_transform(const std::vector<float>& c2w, int64_t n,
                                     double T_out[16]) {
    std::fill(T_out, T_out + 16, 0.0);
    T_out[0] = T_out[5] = T_out[10] = T_out[15] = 1.0;
    if (n <= 0) return 1.0;
    double up[3] = {0, 0, 0}, center[3] = {0, 0, 0};
    for (int64_t i = 0; i < n; i++) {
        for (int r = 0; r < 3; r++) {
            up[r]     += c2w[i*12 + r*4 + 1];
            center[r] += c2w[i*12 + r*4 + 3];
        }
    }
    double un = std::sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
    for (auto& u : up) u /= std::max(un, 1e-12);
    for (auto& c : center) c /= (double)n;

    double axis[3] = {up[1], -up[0], 0.0};              // up x z
    double s = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1]);
    double c = up[2];
    double R[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
    if (s > 1e-12) {
        for (auto& a : axis) a /= s;
        double C = 1.0 - c;
        double Rr[3][3] = {
            {c + axis[0]*axis[0]*C,         axis[0]*axis[1]*C - axis[2]*s, axis[1]*s},
            {axis[0]*axis[1]*C + axis[2]*s, c + axis[1]*axis[1]*C,        -axis[0]*s},
            {-axis[1]*s,                    axis[0]*s,                     c},
        };
        std::copy(&Rr[0][0], &Rr[0][0] + 9, &R[0][0]);
    } else if (c < 0.0) {
        R[1][1] = -1.0; R[2][2] = -1.0;                 // up == -z: flip
    }

    double max_abs = 0.0;
    for (int64_t i = 0; i < n; i++) {
        double d[3];
        for (int r = 0; r < 3; r++) d[r] = c2w[i*12 + r*4 + 3] - center[r];
        for (int r = 0; r < 3; r++)
            max_abs = std::max(max_abs,
                std::abs(R[r][0]*d[0] + R[r][1]*d[1] + R[r][2]*d[2]));
    }
    double scale_factor = 1.0 / std::max(max_abs, 1e-12);

    // T_n_from_camera = scale * [R_align | -R_align @ center]
    for (int r = 0; r < 3; r++) {
        double t = 0.0;
        for (int col = 0; col < 3; col++) {
            T_out[r*4 + col] = R[r][col] * scale_factor;
            t -= R[r][col] * center[col];
        }
        T_out[r*4 + 3] = t * scale_factor;
    }
    return scale_factor;
}

double compute_normalized_scale_factor(const std::vector<float>& c2w, int64_t n) {
    double T[16];
    return compute_normalized_transform(c2w, n, T);
}

void invert_affine4x4(const double in[16], double out[16]) {
    double a = in[0], b = in[1], c = in[2],
           d = in[4], e = in[5], g = in[6],
           h = in[8], i = in[9], j = in[10];
    double det = a*(e*j - g*i) - b*(d*j - g*h) + c*(d*i - e*h);
    if (std::abs(det) < 1e-30)
        throw std::runtime_error("invert_affine4x4: singular matrix");
    double inv = 1.0 / det;
    double Ai[3][3] = {
        {(e*j - g*i)*inv, (c*i - b*j)*inv, (b*g - c*e)*inv},
        {(g*h - d*j)*inv, (a*j - c*h)*inv, (c*d - a*g)*inv},
        {(d*i - e*h)*inv, (b*h - a*i)*inv, (a*e - b*d)*inv},
    };
    for (int r = 0; r < 3; r++) {
        for (int col = 0; col < 3; col++) out[r*4 + col] = Ai[r][col];
        out[r*4 + 3] = -(Ai[r][0]*in[3] + Ai[r][1]*in[7] + Ai[r][2]*in[11]);
    }
    out[12] = out[13] = out[14] = 0.0;
    out[15] = 1.0;
}


// ---------------------------------------------------------------------------
// eval_mode train subset (get_train_eval_split_*, dataparser.py:51-129)
// ---------------------------------------------------------------------------
std::vector<int64_t> train_subset(int64_t n, const std::vector<std::string>& names,
                                  const DatasetParserConfig& cfg) {
    std::vector<int64_t> keep;
    if (cfg.eval_mode == "all") {
        keep.resize(n);
        for (int64_t i = 0; i < n; i++) keep[i] = i;
    } else if (cfg.eval_mode == "fraction") {
        int64_t num_train = (int64_t)std::ceil((double)n * cfg.train_split_fraction);
        std::vector<char> flag(n, 0);
        for (int64_t k = 0; k < num_train; k++)
            flag[(num_train == 1) ? 0 : (int64_t)std::llround(
                (double)k * (double)(n - 1) / (double)(num_train - 1))] = 1;
        for (int64_t i = 0; i < n; i++) if (flag[i]) keep.push_back(i);
    } else if (cfg.eval_mode == "interval") {
        for (int64_t i = 0; i < n; i++)
            if (cfg.eval_interval <= 0 || i % cfg.eval_interval != 0)
                keep.push_back(i);
    } else if (cfg.eval_mode == "filename") {
        for (int64_t i = 0; i < n; i++) {
            std::string base = fs::path(names[i]).filename().string();
            if (base.find("train") != std::string::npos) keep.push_back(i);
            else if (base.find("eval") == std::string::npos)
                throw std::runtime_error(
                    "eval_mode=filename requires 'train'/'eval' in every image "
                    "name; got " + names[i]);
        }
    } else {
        throw std::runtime_error("unknown eval_mode " + cfg.eval_mode);
    }
    if (keep.empty())
        throw std::runtime_error("eval_mode split left no training images");
    return keep;
}


void assign_val_split(ParsedDataset& ds, float validation_fraction) {
    const int64_t n = ds.num_cameras;
    std::vector<char> is_val(n, 0);
    if (validation_fraction > 0.0f && n > 1) {
        int64_t num_train = (int64_t)std::ceil((double)n * (1.0 - validation_fraction));
        std::vector<char> is_train(n, 0);
        for (int64_t k = 0; k < num_train; k++) {
            int64_t idx = (num_train == 1) ? 0
                : (int64_t)std::llround((double)k * (double)(n - 1) / (double)(num_train - 1));
            is_train[idx] = 1;
        }
        for (int64_t i = 0; i < n; i++) is_val[i] = !is_train[i];
    }
    ds.train_indices.clear();
    ds.val_indices.clear();
    for (int64_t i = 0; i < n; i++)
        (is_val[i] ? ds.val_indices : ds.train_indices).push_back((int32_t)i);
}


// ---------------------------------------------------------------------------
// Auxiliary buffer discovery (dataparser.py _add_auxiliary_buffers, trimmed)
// ---------------------------------------------------------------------------
std::string find_aux_file(const std::string& aux_dir_s, const std::string& rel_name,
                          const char* suffix_tag) {
    fs::path aux_dir(aux_dir_s);
    if (!fs::is_directory(aux_dir)) return "";
    fs::path rel(rel_name);
    std::string stem_rel = (rel.parent_path() / rel.stem()).string();
    const std::string exts[] = {".png", ".PNG", ".jpg", ".JPG", ".jpeg", ".JPEG"};
    std::vector<std::string> candidates;
    for (const auto& e : exts) {
        candidates.push_back(rel_name + e);   // image.jpg.png
        candidates.push_back(stem_rel + e);   // image.png
    }
    candidates.push_back(stem_rel + "_" + suffix_tag + ".png");   // image_mask.png
    for (const auto& cand : candidates) {
        fs::path p = aux_dir / cand;
        if (fs::exists(p)) return p.string();
    }
    return "";
}


// ---------------------------------------------------------------------------
// Outlier rejection (geometric_median, dataparser.py:138-173 + 319-325)
// ---------------------------------------------------------------------------
namespace {
double median_of(std::vector<double> v) {
    if (v.empty()) return 0.0;
    size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double hi = v[mid];
    if (v.size() % 2 == 1) return hi;
    double lo = *std::max_element(v.begin(), v.begin() + mid);
    return 0.5 * (lo + hi);
}
}  // namespace

std::vector<char> outlier_keep_mask(const std::vector<double>& pos,
                                    int64_t n, float threshold) {
    std::vector<char> keep(n, 1);
    if (!(threshold < std::numeric_limits<float>::infinity()) || n == 0)
        return keep;

    // Geometric median via Weiszfeld with the zero-distance correction
    // (dataparser.py geometric_median, eps=0, maxiter=10).
    double y[3];
    for (int d = 0; d < 3; d++) {
        std::vector<double> col(n);
        for (int64_t i = 0; i < n; i++) col[i] = pos[i*3 + d];
        y[d] = median_of(col);
    }
    for (int it = 0; it < 10; it++) {
        double T[3] = {0, 0, 0}, Dinvs = 0.0;
        int64_t num_zeros = 0;
        for (int64_t i = 0; i < n; i++) {
            double dx = pos[i*3] - y[0], dy = pos[i*3+1] - y[1], dz = pos[i*3+2] - y[2];
            double D = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (D == 0.0) { num_zeros++; continue; }
            double w = 1.0 / D;
            Dinvs += w;
            for (int d = 0; d < 3; d++) T[d] += w * pos[i*3 + d];
        }
        double y1[3];
        if (num_zeros == n) break;
        if (Dinvs > 0) for (int d = 0; d < 3; d++) T[d] /= Dinvs;
        if (num_zeros == 0) {
            for (int d = 0; d < 3; d++) y1[d] = T[d];
        } else {
            double R[3], r = 0.0;
            for (int d = 0; d < 3; d++) { R[d] = (T[d] - y[d]) * Dinvs; r += R[d]*R[d]; }
            r = std::sqrt(r);
            double rinv = (r == 0.0) ? 0.0 : (double)num_zeros / r;
            double a = std::max(0.0, 1.0 - rinv), b = std::min(1.0, rinv);
            for (int d = 0; d < 3; d++) y1[d] = a * T[d] + b * y[d];
        }
        double diff = 0.0;
        for (int d = 0; d < 3; d++) {
            diff += (y[d] - y1[d]) * (y[d] - y1[d]);
            y[d] = y1[d];
        }
        if (diff == 0.0) break;
    }

    std::vector<double> dist(n);
    for (int64_t i = 0; i < n; i++) {
        double dx = pos[i*3] - y[0], dy = pos[i*3+1] - y[1], dz = pos[i*3+2] - y[2];
        dist[i] = std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    double mad = median_of(dist);
    for (int64_t i = 0; i < n; i++)
        keep[i] = dist[i] <= (double)threshold * mad;
    return keep;
}

}  // namespace dsparse


// ===========================================================================
// bake_post_split (trainer.py _setup_cpp_data_manager:257-414)
// ===========================================================================

namespace {

// 5- / 6-face cubemap axis tables. MUST match DataManager.cpp's
// kAxesFisheye5 / kAxesEquirect6 (and trainer.py:332-346) -- the engine's
// GPU warp kernel reads its own copies; these drive the camera poses.
const double kAxes5[5][3][3] = {
    {{1,0,0},{0,1,0},{0,0,1}},
    {{0,1,0},{0,0,1},{1,0,0}},
    {{-1,0,0},{0,0,1},{0,1,0}},
    {{0,-1,0},{0,0,1},{-1,0,0}},
    {{1,0,0},{0,0,1},{0,-1,0}},
};
const double kAxes6[6][3][3] = {
    {{1,0,0},{0,1,0},{0,0,1}},
    {{0,0,-1},{0,1,0},{1,0,0}},
    {{-1,0,0},{0,0,1},{0,1,0}},
    {{0,-1,0},{0,0,1},{-1,0,0}},
    {{1,0,0},{0,0,1},{0,-1,0}},
    {{1,0,0},{0,-1,0},{0,0,-1}},
};

}  // namespace

PostSplitCameras bake_post_split(const ParsedDataset& ds,
                                 bool warp_to_pinhole,
                                 bool warp_spherical_to_pinhole) {
    const int64_t N = ds.num_cameras;
    const int PINHOLE_V   = (int)camera_model_from_name("PINHOLE");
    const int FISHEYE_V   = (int)camera_model_from_name("FISHEYE");
    const int EQUISOLID_V = (int)camera_model_from_name("EQUISOLID");
    const int EQUIRECT_V  = (int)camera_model_from_name("EQUIRECTANGULAR");

    PostSplitCameras out;
    out.K_per_camera.resize(N);
    out.post_offsets.resize(N);
    int64_t n_post = 0;
    for (int64_t i = 0; i < N; i++) {
        int cm = ds.camera_models[i], K = 1;
        if (cm == EQUIRECT_V)
            K = warp_spherical_to_pinhole ? 6 : 1;
        else if (warp_to_pinhole && (cm == FISHEYE_V || cm == EQUISOLID_V))
            K = 5;
        out.K_per_camera[i] = K;
        out.post_offsets[i] = (int32_t)n_post;
        n_post += K;
        if (K > 1) {
            out.any_warp = true;
            if (cm != EQUIRECT_V) out.any_fisheye_warp = true;
        } else if (cm == EQUIRECT_V) {
            out.direct_equirect = true;
        }
    }
    out.n_post = n_post;

    out.viewmats.assign(n_post * 16, 0.f);
    out.intrins.assign(n_post * 4, 0.f);
    out.dist_coeffs.assign(n_post * 10, 0.f);
    out.c2w_flip.assign(n_post * 12, 0.f);
    out.post_widths.assign(n_post, 0);
    out.post_heights.assign(n_post, 0);
    out.post_models.assign(n_post, PINHOLE_V);

    // POST-split c2w staging (double; algebra matches trainer.py:355-401).
    std::vector<double> post_c2w(n_post * 12);
    static const double D[3] = {1.0, -1.0, -1.0};

    for (int64_t i = 0; i < N; i++) {
        const int K   = out.K_per_camera[i];
        const int off = out.post_offsets[i];
        double ci[3][4];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 4; c++)
                ci[r][c] = ds.c2w[i*12 + r*4 + c];

        if (K == 1) {
            std::copy(&ci[0][0], &ci[0][0] + 12, &post_c2w[off*12]);
            out.post_widths[off]  = ds.widths[i];
            out.post_heights[off] = ds.heights[i];
            out.post_models[off]  = ds.camera_models[i];
            if (ds.camera_models[i] == EQUIRECT_V) {
                // Direct equirectangular: canonical panorama intrinsics
                // matching the equirect projection kernel (trainer.py:361-369).
                double f = ds.widths[i] / (2.0 * M_PI);
                out.intrins[off*4 + 0] = (float)f;
                out.intrins[off*4 + 1] = (float)f;
                out.intrins[off*4 + 2] = ds.widths[i]  / 2.0f;
                out.intrins[off*4 + 3] = ds.heights[i] / 2.0f;
            } else {
                std::copy(&ds.intrins[i*4], &ds.intrins[i*4] + 4,
                          &out.intrins[off*4]);
            }
            std::copy(&ds.dist_coeffs[i*10], &ds.dist_coeffs[i*10] + 10,
                      &out.dist_coeffs[off*10]);
            continue;
        }

        // K > 1: cubemap-face expansion.
        //   R' = R_c2w * diag(1,-1,-1)   (columns; OpenGL -> OpenCV)
        //   R_post[k] = (R' @ axes[k]^T) * diag(1,-1,-1)
        const auto& axes = (K == 6) ? kAxes6 : kAxes5;
        int out_shape = (int)std::ceil(
            std::sqrt((double)ds.heights[i] * (double)ds.widths[i] / (double)K));
        double R[3][3];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                R[r][c] = ci[r][c] * D[c];
        for (int k = 0; k < K; k++) {
            double* pc = &post_c2w[(off + k) * 12];
            for (int r = 0; r < 3; r++) {
                for (int c = 0; c < 3; c++) {
                    // (R @ axes[k]^T)[r][c] = sum_m R[r][m] * axes[k][c][m]
                    double v = R[r][0]*axes[k][c][0] + R[r][1]*axes[k][c][1]
                             + R[r][2]*axes[k][c][2];
                    pc[r*4 + c] = v * D[c];
                }
                pc[r*4 + 3] = ci[r][3];
            }
            float s = out_shape / 2.0f;
            for (int c = 0; c < 4; c++) out.intrins[(off + k)*4 + c] = s;
            // dist_coeffs stay zero; faces render as canonical PINHOLE.
            out.post_widths[off + k]  = out_shape;
            out.post_heights[off + k] = out_shape;
        }
    }

    // c2w -> engine viewmat over the POST arrays (trainer.py:403-414):
    // R_v = R_post * diag(1,-1,-1) (columns); viewmat = [R_v^T | -R_v^T t].
    // [R_v | t] is also the y/z-flipped c2w form the viewer blit expects
    // (annotation.py:103-106).
    for (int64_t j = 0; j < n_post; j++) {
        const double* pc = &post_c2w[j*12];
        double Rv[3][3], t[3];
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) Rv[r][c] = pc[r*4 + c] * D[c];
            t[r] = pc[r*4 + 3];
        }
        float* vm = &out.viewmats[j*16];
        float* cf = &out.c2w_flip[j*12];
        for (int r = 0; r < 3; r++) {
            double ti = 0.0;
            for (int c = 0; c < 3; c++) {
                vm[r*4 + c] = (float)Rv[c][r];
                cf[r*4 + c] = (float)Rv[r][c];
                ti -= Rv[c][r] * t[c];
            }
            vm[r*4 + 3] = (float)ti;
            cf[r*4 + 3] = (float)t[r];
        }
        vm[15] = 1.f;
    }

    if (out.any_warp) {
        out.input_intrins     = ds.intrins;
        out.input_dist_coeffs = ds.dist_coeffs;
    }
    return out;
}
