// main.cpp -- standalone CLI trainer (no Python, no tyro).
//
// Command shape mirrors spirulae-train:
//     ssplat-train [<preset>] --data <dir> [--flag value ...]
// where <preset> is one of the tyro subcommands (3dgs = default,
// 360-camera, in-the-wild, linear-color, synthetic, meshing,
// academic-baseline). Flags are the FLATTENED training config fields
// (--sh-degree, not --model.sh-degree); '-' and '_' are interchangeable;
// booleans take a value (--warp-to-pinhole 1). The config struct, flag
// table, and preset appliers are code-generated from the Python config
// dataclasses -- see generated/cli_config.h and generate_cli_config.py.
//
// The per-step engine plumbing below is a direct port of the Python
// managed-path trainer:
//   model.py  _build_loss_weights / engine_train_step_managed
//   core.py   _build_optim_config / _build_densify_config
//   optimizer.py get_scheduled_lr
//   model.py  populate_modules (seeding) / _maybe_init_bilagrid /
//             background + color-space init
//   trainer.py train() save cadence / _setup_cpp_data_manager (non-warp)
//
// Not ported yet (each throws or warns at its check in main):
//   warp_to_pinhole camera splitting, resume, eval/early-stop, viewer,
//   camera optimizer, BVH primitive, nerfstudio/metashape data formats.

#include "../Engine.h"
#include "DatasetParser.h"
#include "Viewer.h"
#include "generated/cli_config.h"

#ifndef _WIN32
#include <ftw.h>
#endif

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <filesystem>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;


// ===========================================================================
// Generic flag parsing over the generated SSPLAT_CONFIG_FIELDS table
// ===========================================================================

namespace {

std::string normalize_key(std::string s) {
    while (!s.empty() && s[0] == '-') s.erase(s.begin());
    for (auto& c : s) if (c == '-') c = '_';
    return s;
}

[[noreturn]] void bad_value(const std::string& key, const std::string& val,
                            const std::string& why) {
    throw std::runtime_error("--" + key + ": invalid value '" + val + "' (" + why + ")");
}

// Single-token setters; overload resolution picks the right one per field.
void parse_token(bool& out, const std::string& key, const std::string& v) {
    if (v == "1" || v == "true" || v == "True")  { out = true;  return; }
    if (v == "0" || v == "false" || v == "False") { out = false; return; }
    bad_value(key, v, "expected 0/1/true/false");
}
void parse_token(int& out, const std::string& key, const std::string& v) {
    try { size_t p; out = std::stoi(v, &p); if (p == v.size()) return; } catch (...) {}
    bad_value(key, v, "expected an integer");
}
void parse_token(float& out, const std::string& key, const std::string& v) {
    if (v == "inf" || v == "Infinity") { out = std::numeric_limits<float>::infinity(); return; }
    try { size_t p; out = std::stof(v, &p); if (p == v.size()) return; } catch (...) {}
    bad_value(key, v, "expected a number");
}
void parse_token(std::string& out, const std::string&, const std::string& v) {
    out = (v == "none" || v == "None" || v == "null") ? "" : v;
}
template <typename T>
void parse_token(std::optional<T>& out, const std::string& key, const std::string& v) {
    if (v == "none" || v == "None" || v == "null") { out = std::nullopt; return; }
    T t{}; parse_token(t, key, v); out = t;
}

// Arity per field type (tuples consume N CLI tokens).
template <typename T> struct ArgArity { static constexpr int value = 1; };
template <typename T, size_t N> struct ArgArity<std::array<T, N>> {
    static constexpr int value = (int)N;
};

template <typename T>
void consume(T& out, const std::string& key, int argc, char** argv, int& i) {
    if (i + 1 >= argc) throw std::runtime_error("--" + key + ": missing value");
    parse_token(out, key, argv[++i]);
}
template <typename T, size_t N>
void consume(std::array<T, N>& out, const std::string& key, int argc, char** argv, int& i) {
    for (size_t k = 0; k < N; k++) {
        if (i + 1 >= argc)
            throw std::runtime_error("--" + key + ": expected " + std::to_string(N) + " values");
        parse_token(out[k], key, argv[++i]);
    }
}

// Choices validation for string fields. `choices` is "a|b|c" ('' free-form);
// a lone "none" marks an optional free-form string.
void check_choices(const std::string& value, const std::string& key, const char* choices) {
    std::string ch = choices;
    if (ch.empty() || ch == "none") return;
    std::string want = value.empty() ? "none" : value;
    size_t pos = 0;
    while (pos <= ch.size()) {
        size_t bar = ch.find('|', pos);
        std::string tok = ch.substr(pos, bar == std::string::npos ? std::string::npos : bar - pos);
        if (tok == want) return;
        if (bar == std::string::npos) break;
        pos = bar + 1;
    }
    bad_value(key, value, "expected one of: " + ch);
}
void check_choices(...) {}  // non-string fields

// Returns false if the key is unknown.
bool set_config_field(SsplatConfig& c, const std::string& key,
                      int argc, char** argv, int& i) {
#define SSPLAT_TRY_SET(member, cli_key, pyname, group, choices, help)          \
    if (key == cli_key) {                                                      \
        consume(c.member, key, argc, argv, i);                                 \
        check_choices(c.member, key, choices);                                 \
        return true;                                                           \
    }
    SSPLAT_CONFIG_FIELDS(SSPLAT_TRY_SET)
#undef SSPLAT_TRY_SET
    return false;
}

// ---- Help / config dump ---------------------------------------------------

std::string value_str(bool v)               { return v ? "1" : "0"; }
std::string value_str(int v)                { return std::to_string(v); }
std::string value_str(float v) {
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[32]; std::snprintf(buf, sizeof buf, "%g", v); return buf;
}
std::string value_str(const std::string& v) { return v.empty() ? "none" : v; }
template <typename T> std::string value_str(const std::optional<T>& v) {
    return v ? value_str(*v) : "none";
}
template <typename T, size_t N> std::string value_str(const std::array<T, N>& v) {
    std::string s;
    for (size_t i = 0; i < N; i++) s += (i ? " " : "") + value_str(v[i]);
    return s;
}

void print_help(const char* argv0, const SsplatConfig& c) {
    std::printf("usage: %s [<preset>] --data <dataset_dir> [--flag value ...]\n\n", argv0);
    std::printf("presets (tyro subcommands; default: 3dgs):\n");
    for (const auto& p : kSsplatPresets)
        std::printf("  %-18s %s\n", p.name, p.help);
    std::printf("\nflags ('-' and '_' interchangeable; bools take 0/1; 'none' clears "
                "optional values;\n defaults shown for the selected preset):\n");
    const char* cur_group = "";
#define SSPLAT_PRINT_HELP(member, cli_key, pyname, group, choices, help)       \
    if (std::strcmp(cur_group, group) != 0) {                                  \
        cur_group = group;                                                     \
        std::printf("\n  [%s]\n", group);                                      \
    }                                                                          \
    {                                                                          \
        std::string h = help;                                                  \
        size_t dot = h.find(". ");                                             \
        if (dot != std::string::npos) h = h.substr(0, dot + 1);                \
        if (h.size() > 110) h = h.substr(0, 107) + "...";                      \
        std::string ch = choices;                                              \
        std::string key_disp = cli_key;                                        \
        for (auto& ck : key_disp) if (ck == '_') ck = '-';                     \
        std::printf("  --%-38s [%s]%s%s\n      %s\n", key_disp.c_str(),        \
                    value_str(c.member).c_str(),                               \
                    (ch.empty() || ch == "none") ? "" : (" {" + ch + "}").c_str(), \
                    "", h.c_str());                                            \
    }
    SSPLAT_CONFIG_FIELDS(SSPLAT_PRINT_HELP)
#undef SSPLAT_PRINT_HELP
}

// JSON value emitters for the config.json dump (Python-json compatible:
// float('inf') serializes as Infinity, None as null).
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

// Nested-by-group dump, using the original Python field names. Close to the
// trainer's config.json (enough for humans / tooling; the Python --resume
// path needs a tyro round-trip and is not supported by the CLI yet).
void save_config_json(const SsplatConfig& c, const fs::path& out_dir,
                      const std::string& preset) {
    FILE* f = std::fopen((out_dir / "config.json").string().c_str(), "w");
    if (!f) throw std::runtime_error("cannot write config.json");
    std::fprintf(f, "{\n    \"preset\": \"%s\"", preset.c_str());
    const char* open_group = "";
#define SSPLAT_DUMP(member, cli_key, pyname, group, choices, help)             \
    if (std::strcmp(group, "trainer") == 0) {                                  \
        std::fprintf(f, ",\n    \"%s\": %s", pyname, json_str(c.member).c_str()); \
    } else {                                                                   \
        if (std::strcmp(open_group, group) != 0) {                             \
            if (*open_group) std::fprintf(f, "\n    }");                       \
            std::fprintf(f, ",\n    \"%s\": {", group);                        \
            open_group = group;                                                \
            std::fprintf(f, "\n        \"%s\": %s", pyname, json_str(c.member).c_str()); \
        } else {                                                               \
            std::fprintf(f, ",\n        \"%s\": %s", pyname, json_str(c.member).c_str()); \
        }                                                                      \
    }
    SSPLAT_CONFIG_FIELDS(SSPLAT_DUMP)
#undef SSPLAT_DUMP
    if (*open_group) std::fprintf(f, "\n    }");
    std::fprintf(f, "\n}\n");
    std::fclose(f);
}


// Debug dump for numeric verification against the Python dataparser /
// trainer camera algebra (SSPLAT_DUMP_CAMERAS=<path> env). Full precision.
void dump_cameras_json(const char* path, const ParsedDataset& ds,
                       const PostSplitCameras& post) {
    FILE* f = std::fopen(path, "w");
    if (!f) throw std::runtime_error(std::string("cannot write ") + path);
    auto arr_f = [&](const char* k, const std::vector<float>& v) {
        std::fprintf(f, "\"%s\": [", k);
        for (size_t i = 0; i < v.size(); i++)
            std::fprintf(f, "%s%.9g", i ? "," : "", v[i]);
        std::fprintf(f, "]");
    };
    auto arr_i = [&](const char* k, const std::vector<int32_t>& v) {
        std::fprintf(f, "\"%s\": [", k);
        for (size_t i = 0; i < v.size(); i++)
            std::fprintf(f, "%s%d", i ? "," : "", v[i]);
        std::fprintf(f, "]");
    };
    std::fprintf(f, "{\n\"num_cameras\": %lld,\n\"n_post\": %lld,\n"
                 "\"train_frame_scale\": %.9g,\n\"num_points\": %lld,\n",
                 (long long)ds.num_cameras, (long long)post.n_post,
                 ds.train_frame_scale, (long long)ds.points.num());
    std::fprintf(f, "\"image_filenames\": [");
    for (size_t i = 0; i < ds.image_filenames.size(); i++)
        std::fprintf(f, "%s\"%s\"", i ? "," : "", ds.image_filenames[i].c_str());
    std::fprintf(f, "],\n");
    arr_i("camera_models", ds.camera_models);   std::fprintf(f, ",\n");
    arr_i("widths", ds.widths);                 std::fprintf(f, ",\n");
    arr_i("heights", ds.heights);               std::fprintf(f, ",\n");
    arr_f("c2w", ds.c2w);                       std::fprintf(f, ",\n");
    arr_i("K_per_camera", post.K_per_camera);   std::fprintf(f, ",\n");
    arr_i("post_offsets", post.post_offsets);   std::fprintf(f, ",\n");
    arr_f("viewmats", post.viewmats);           std::fprintf(f, ",\n");
    arr_f("intrins", post.intrins);             std::fprintf(f, ",\n");
    arr_f("dist_coeffs", post.dist_coeffs);     std::fprintf(f, ",\n");
    arr_f("input_intrins", post.input_intrins); std::fprintf(f, ",\n");
    arr_f("input_dist_coeffs", post.input_dist_coeffs); std::fprintf(f, ",\n");
    std::vector<float> t2n(ds.train_to_normalized.begin(), ds.train_to_normalized.end());
    arr_f("train_to_normalized", t2n);          std::fprintf(f, ",\n");
    std::vector<float> pts_head(ds.points.xyz.begin(),
        ds.points.xyz.begin() + std::min<size_t>(ds.points.xyz.size(), 30));
    arr_f("points_head", pts_head);
    std::fprintf(f, "\n}\n");
    std::fclose(f);
}


// ===========================================================================
// Color-space handling (port of _wrapper_per_pixel.get_color_transform_matrix
// + the resolution logic in model.py populate_modules:612-631 / __init__:532)
// ===========================================================================

using Mat3f = std::array<float, 9>;

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

struct ColorResolution {
    std::string splat_gamut;   // "" = Rec.709 / none
    std::string image_gamut;
    bool splat_linear  = false;
    bool image_linear  = false;
    bool convert_seed  = false;  // convert_initial_point_cloud_color resolved
};

// model.py populate_modules:612-623, applied to the config in place there;
// pure-function port here.
ColorResolution resolve_color(const SsplatConfig& c) {
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
                   std::optional<float> lr_final = std::nullopt,
                   std::optional<int> warmup = std::nullopt) {
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

struct SeedSplats {
    int64_t num = 0;                 // live splats (<= cap rows)
    std::vector<float> means;        // [cap, 3]
    std::vector<float> quats;        // [cap, 4]
    std::vector<float> scales;       // [cap, 3]  (log)
    std::vector<float> opacities;    // [cap, 1]  (logit)
    std::vector<float> features_dc;  // [cap, 3]
    std::vector<float> features_sh;  // [cap, dim_sh-1, 3]  (zeros)
};

// Mean distance to the k nearest neighbors via a uniform spatial hash
// (27-cell neighborhood). Approximate stand-in for the Python
// k_nearest_neighbor; falls back to the cell size when a neighborhood is
// empty.
std::vector<float> knn_mean_dist(const std::vector<float>& xyz, int64_t n, int k) {
    float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
    for (int64_t i = 0; i < n; i++)
        for (int d = 0; d < 3; d++) {
            lo[d] = std::min(lo[d], xyz[i*3+d]);
            hi[d] = std::max(hi[d], xyz[i*3+d]);
        }
    double vol = std::max((double)(hi[0]-lo[0]) * (hi[1]-lo[1]) * (hi[2]-lo[2]), 1e-30);
    float cell = (float)std::cbrt(vol / std::max<int64_t>(n, 1) * 2.0);
    cell = std::max(cell, 1e-12f);

    auto cell_key = [&](float x, float y, float z) -> uint64_t {
        int64_t cx = (int64_t)std::floor((x - lo[0]) / cell);
        int64_t cy = (int64_t)std::floor((y - lo[1]) / cell);
        int64_t cz = (int64_t)std::floor((z - lo[2]) / cell);
        return ((uint64_t)(cx & 0x1FFFFF) << 42) |
               ((uint64_t)(cy & 0x1FFFFF) << 21) |
               ((uint64_t)(cz & 0x1FFFFF));
    };
    std::unordered_map<uint64_t, std::vector<int32_t>> grid;
    grid.reserve(n);
    for (int64_t i = 0; i < n; i++)
        grid[cell_key(xyz[i*3], xyz[i*3+1], xyz[i*3+2])].push_back((int32_t)i);

    std::vector<float> out(n);
    std::vector<float> d2s;
    for (int64_t i = 0; i < n; i++) {
        d2s.clear();
        float x = xyz[i*3], y = xyz[i*3+1], z = xyz[i*3+2];
        for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
        for (int dz = -1; dz <= 1; dz++) {
            auto it = grid.find(cell_key(x + dx*cell, y + dy*cell, z + dz*cell));
            if (it == grid.end()) continue;
            for (int32_t j : it->second) {
                if (j == (int32_t)i) continue;
                float ddx = xyz[j*3]-x, ddy = xyz[j*3+1]-y, ddz = xyz[j*3+2]-z;
                d2s.push_back(ddx*ddx + ddy*ddy + ddz*ddz);
            }
        }
        if (d2s.empty()) { out[i] = cell; continue; }
        int kk = std::min<int>(k, (int)d2s.size());
        std::partial_sort(d2s.begin(), d2s.begin() + kk, d2s.end());
        double acc = 0.0;
        for (int j = 0; j < kk; j++) acc += d2s[j];
        out[i] = (float)std::sqrt(acc / kk);   // sqrt(mean(d^2)), model.py:582
    }
    return out;
}

SeedSplats seed_splats(const ColmapPoints3D& pts, const SsplatConfig& cfg,
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
    std::vector<float> nn = knn_mean_dist(s.means, num, 4);
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
// Per-step EngineStepConfig -- ports of model.py _build_loss_weights /
// engine_train_step_managed and core.py _build_optim_config /
// _build_densify_config.
// ===========================================================================

struct RunState {
    float train_frame_scale = 1.0f;
    bool  splat_linear = false;
    bool  bilagrid_rgb_init    = false;
    bool  bilagrid_depth_init  = false;
    bool  bilagrid_normal_init = false;
    bool  ppisp_init           = false;
};

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

// model.py:1483 _build_loss_weights.
std::array<float, (int)LossWeightIndex::length>
build_loss_weights(const SsplatConfig& c, int step) {
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

EngineStepConfig build_step_config(const SsplatConfig& c, const RunState& st, int step) {
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
        cfg.ppisp.run_before_bilagrid = c.apply_ppisp_before_bilagrid;
    }

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

}  // namespace


// ===========================================================================
// main
// ===========================================================================

int main(int argc, char** argv) {
    try {
        // ---- Preset + flags ------------------------------------------------
        std::string preset = "3dgs";
        int argi = 1;
        if (argi < argc && argv[argi][0] != '-') preset = argv[argi++];

        SsplatConfig cfg;
        if (!ssplat_apply_preset(cfg, preset)) {
            std::string names;
            for (const auto& p : kSsplatPresets) names += std::string(" ") + p.name;
            throw std::runtime_error("unknown preset '" + preset + "'; expected one of:" + names);
        }

        for (int i = argi; i < argc; i++) {
            std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") { print_help(argv[0], cfg); return 0; }
            if (arg.rfind("--", 0) != 0)
                throw std::runtime_error("unexpected argument: " + arg + " (flags are --key value)");
            // --key=value form: re-parse via a 2-token mini-argv. Tuple
            // fields need N separate tokens, so = only supports arity 1.
            std::string key = arg;
            size_t eq = key.find('=');
            if (eq != std::string::npos) {
                static std::vector<std::string> keep;   // token storage
                keep.push_back(key.substr(eq + 1));
                key = key.substr(0, eq);
                char* mini[2] = {const_cast<char*>(""), keep.back().data()};
                int mi = 0;
                if (!set_config_field(cfg, normalize_key(key), 2, mini, mi))
                    throw std::runtime_error("unknown flag: " + key);
                continue;
            }
            if (!set_config_field(cfg, normalize_key(key), argc, argv, i))
                throw std::runtime_error("unknown flag: " + key +
                                         " (see --help for the full list)");
        }

#define SSPLAT_CHECK_REQUIRED(member)                                          \
        if (cfg.member.empty())                                                \
            throw std::runtime_error("--" #member " is required");
        SSPLAT_CONFIG_REQUIRED_FIELDS(SSPLAT_CHECK_REQUIRED)
#undef SSPLAT_CHECK_REQUIRED

        // ---- Unported-feature guards (fail early, like the Python trainer) -
        auto not_impl = [](const char* what) {
            throw std::runtime_error(std::string(what) +
                " is not supported by the standalone CLI yet; use spirulae-train");
        };
        if (!cfg.resume.empty())            not_impl("--resume");
        if (cfg.use_bvh)                    not_impl("--use-bvh");
        if (cfg.use_camera_optimizer)       not_impl("--use-camera-optimizer");
        if (cfg.deblur_training_images)     not_impl("--deblur-training-images");
        if (!cfg.optimizer_offload.empty()) not_impl("--optimizer-offload");
        if (cfg.save_eval_images)           not_impl("--save-eval-images");
        if (cfg.data_format == "metashape")
            not_impl("--data-format metashape");
        if (cfg.cache_images == "gpu")      not_impl("--cache-images gpu");
        if (cfg.rescale_camera_to_fit < 0)  not_impl("--rescale-camera-to-fit auto-detect");
        if (cfg.num_downscales > 0)
            std::fprintf(stderr, "warning: --num-downscales is a Python-data-path "
                         "feature; ignored by the managed engine path\n");
        if (cfg.validation_fraction > 0)
            std::fprintf(stderr, "warning: validation images are held out but "
                         "early stopping / eval is not ported yet\n");
        if (cfg.orientation_method != "up" || cfg.center_method != "poses")
            std::fprintf(stderr, "warning: orientation/center method '%s'/'%s' "
                         "approximated as 'up'/'poses' (affects only "
                         "train_frame_scale)\n",
                         cfg.orientation_method.c_str(), cfg.center_method.c_str());
        if (cfg.train_frame != "points")
            not_impl(("--train-frame " + cfg.train_frame).c_str());
        if (cfg.primitive != "3dgs" && cfg.primitive != "mip" && cfg.primitive != "3dgut")
            not_impl(("--primitive " + cfg.primitive).c_str());

        // ---- Parse dataset -------------------------------------------------
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
        ParsedDataset ds = parse_dataset(cfg.data, pcfg, cfg.data_format);

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
        PostSplitCameras post = bake_post_split(
            ds, cfg.warp_to_pinhole, cfg.warp_spherical_to_pinhole);

        // Warp-path feature guards (trainer.py:294-303).
        bool has_depth  = !ds.depth_filenames.empty()  && cfg.load_depths;
        bool has_normal = !ds.normal_filenames.empty() && cfg.load_normals;
        if (post.direct_equirect && (has_depth || has_normal))
            throw std::runtime_error(
                "Direct equirectangular training (warp_spherical_to_pinhole=0) "
                "does not support depth/normal supervision yet.");

        std::printf("Parsed %lld cameras (%lld post-split), %lld seed points "
                    "(train_frame_scale=%.4g)\n",
                    (long long)ds.num_cameras, (long long)post.n_post,
                    (long long)ds.points.num(), ds.train_frame_scale);

        // Hidden debug flag (set via env to avoid polluting the config):
        // dump parsed + post-split arrays as JSON and exit, for numeric
        // verification against the Python dataparser/trainer algebra.
        if (const char* dump = std::getenv("SSPLAT_DUMP_CAMERAS")) {
            dump_cameras_json(dump, ds, post);
            std::printf("Dumped cameras to %s\n", dump);
            return 0;
        }

        // ---- Output dir (trainer.py _setup_output_dir:524) ------------------
        fs::path out_dir;
        if (!cfg.output_dir_name.empty()) {
            out_dir = fs::path(cfg.output_dir_prefix) / cfg.output_dir_name;
        } else {
            std::time_t t = std::time(nullptr);
            char stamp[32];
            std::strftime(stamp, sizeof stamp, "%Y%m%d-%H%M%S", std::localtime(&t));
            out_dir = fs::path(cfg.output_dir_prefix) /
                      (fs::path(cfg.data).stem().string() + "_" + stamp);
        }
        fs::create_directories(out_dir);
        save_config_json(cfg, out_dir, preset);
        std::printf("Output directory: %s\n", fs::absolute(out_dir).string().c_str());

        // ---- Engine setup (Trainer.__init__ + SpirulaeSplatModel.__init__) --
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

        // ---- DataManager (trainer.py _setup_cpp_data_manager) ---------------
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

        // ---- Bilagrid / PPISP init (model.py _maybe_init_bilagrid:823). ----
        // Enablement conditions are static for the CLI (dataset modalities
        // known up front), so the lazy per-step init collapses to here.
        RunState st;
        st.train_frame_scale = ds.train_frame_scale;
        st.splat_linear      = color.splat_linear;
        if (cfg.quantization_level != 0 && cfg.quantization_level != 1)
            throw std::runtime_error("quantization_level must be 0 or 1");
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

        // ---- Checkpoint helper (trainer.py save_checkpoint:939) -------------
        // Recursive delete via POSIX nftw, NOT std::filesystem::remove_all:
        // libtorch.so (linked via csrc on this branch) exports its own
        // std::filesystem symbols and the dynamic linker binds remove_all to
        // that ABI-incompatible copy, which crashes. Goes away with no-torch.
        auto remove_tree = [](const fs::path& p) {
#ifndef _WIN32
            nftw(p.string().c_str(),
                 [](const char* f, const struct stat*, int, struct FTW*) {
                     return ::remove(f);
                 }, 16, FTW_DEPTH | FTW_PHYS);
#else
            std::filesystem::remove_all(p);
#endif
        };
        auto save_checkpoint = [&](int step) {
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
        };

        // ---- Viewer + train-loop coordination state -------------------------
        // Mirrors Trainer's lock + _render_pending fairness signal
        // (trainer.py:146-153) and the pause event (trainer.py:161-179).
        std::mutex engine_mutex;
        std::atomic<bool> render_pending{false};
        std::atomic<bool> paused{false};
        std::atomic<int>  cur_step{0};
        std::mutex progress_mutex;               // guards the latency window
        std::deque<double> step_latencies;       // last 100, seconds
        auto train_start = std::chrono::steady_clock::now();

        ViewerServer viewer;
        bool viewer_on = !cfg.disable_viewer;
        if (viewer_on) {
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

            ViewerHooks hooks;
            hooks.engine_mutex = &engine_mutex;
            hooks.current_step = [&] { return cur_step.load(); };
            hooks.set_render_pending = [&](bool v) { render_pending = v; };
            hooks.pause_toggle = [&] {
                bool now = !paused.load();
                paused = now;
                return now;
            };
            // Trainer.get_progress port (trainer.py:181-205).
            hooks.progress_json = [&]() -> std::string {
                int step = cur_step.load();
                double elapsed = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - train_start).count();
                double avg = 0.0;
                size_t nlat = 0;
                {
                    std::lock_guard<std::mutex> lk(progress_mutex);
                    for (double v : step_latencies) avg += v;
                    nlat = step_latencies.size();
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
            };

            viewer.start("0.0.0.0", cfg.viewer_port, vc, hooks, post);
            std::printf("Viewer at http://0.0.0.0:%d/ (forward the port for "
                        "remote boxes: ssh -L %d:localhost:%d <host>)\n",
                        cfg.viewer_port, cfg.viewer_port, cfg.viewer_port);
        }

        // ---- Train loop (trainer.py train:741) ------------------------------
        std::string primitive = cfg.primitive;
        bool packed = cfg.packed || cfg.use_bvh;
        auto t0 = std::chrono::steady_clock::now();
        for (int step = 0; step < cfg.num_iterations; step++) {
            // Pause gate + render-fairness yield: give the viewer's render
            // worker an uncontended window to take the engine mutex.
            while (paused.load())
                std::this_thread::sleep_for(std::chrono::milliseconds(50));
            while (render_pending.load())
                std::this_thread::sleep_for(std::chrono::microseconds(500));

            auto step_start = std::chrono::steady_clock::now();
            std::map<std::string, float> losses;
            {
                std::lock_guard<std::mutex> lk(engine_mutex);
                if (step > 0 && cfg.steps_per_save > 0 && step % cfg.steps_per_save == 0)
                    save_checkpoint(step);
                int sh_degree_to_use = step / std::max(cfg.sh_degree_warmup_every, 1);
                EngineStepConfig sc = build_step_config(cfg, st, step);
                losses = engine_train_step_managed(
                    step, cfg.num_iterations, primitive, sh_degree_to_use, packed, sc);
            }
            cur_step = step + 1;
            {
                std::lock_guard<std::mutex> lk(progress_mutex);
                step_latencies.push_back(std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - step_start).count());
                if (step_latencies.size() > 100) step_latencies.pop_front();
            }

            if (step % 100 == 0 || step == cfg.num_iterations - 1) {
                double dt = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - t0).count();
                std::printf("step %6d/%d  splats %lld  [%.1fs]",
                            step, cfg.num_iterations,
                            (long long)engine_get_cur_num_splats(), dt);
                for (const char* k : {"rgb_loss", "ssim", "psnr"}) {
                    auto it = losses.find(k);
                    if (it != losses.end()) std::printf("  %s=%.4g", k, it->second);
                }
                std::printf("\n");
                std::fflush(stdout);
            }
        }

        if (cfg.steps_per_save != 0) {
            std::lock_guard<std::mutex> lk(engine_mutex);
            save_checkpoint(cfg.num_iterations);
        }
        std::printf("Checkpoint saved to: %s\n", fs::absolute(out_dir).string().c_str());
        // TODO: eval pass over val_indices (PSNR/SSIM), early stopping.

        if (viewer_on && cfg.keep_viewer_alive) {
            std::printf("Training complete. Viewer still running -- press "
                        "Ctrl-C to exit.\n");
            std::fflush(stdout);
            for (;;) std::this_thread::sleep_for(std::chrono::seconds(3600));
        }

    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
