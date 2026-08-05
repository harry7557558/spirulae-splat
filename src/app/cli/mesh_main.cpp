// mesh_main.cpp -- standalone mesh-extraction CLI.
//
//     spirula-mesh <checkpoint> [--data <dir>] [--format ply,glb]
//                 [--color none|vertex|texture] [--flag value ...]
//
// <checkpoint> is a run directory (containing config.json and a step-*.ckpt/
// folder), a *.ckpt directory, or a splat.ply file directly. The dataset
// (for camera-based occupancy + colors) defaults to config.json's `data`;
// pass --no-data to mesh from Gaussian densities only.
//
// It loads the raw
// (un-activated) Gaussians from the checkpoint PLY, parses the dataset with
// the same C++ parsers the CLI trainer uses, and calls meshing::generate_mesh
// (Meshing.h) for the heavy lifting.

#include "app/Tools.h"

#include "mesh/Meshing.h"
#include "core/Camera.h"
#include "data/DatasetParser.h"
#include "data/Json.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ===========================================================================
// splat.ply loader (raw, un-activated Gaussians; float32 binary LE as written
// by EngineCheckpoint.cpp)
// ===========================================================================

struct SplatPly {
    int64_t n = 0;
    std::vector<float> means;        // [n*3]
    std::vector<float> quats;        // [n*4] (w,x,y,z as stored in rot_0..3)
    std::vector<float> scales;       // [n*3] log
    std::vector<float> opacities;    // [n]   logit
    std::vector<float> features_dc;  // [n*3]
};

SplatPly load_splat_ply(const fs::path& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path.string());

    std::string line;
    if (!std::getline(f, line) || line.rfind("ply", 0) != 0)
        throw std::runtime_error(path.string() + ": not a PLY file");

    int64_t n_vertex = -1;
    std::vector<std::string> props;   // property names, in file order
    bool in_vertex = false, binary_le = false;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream ss(line);
        std::string tok;
        ss >> tok;
        if (tok == "format") {
            ss >> tok;
            binary_le = tok == "binary_little_endian";
        } else if (tok == "element") {
            std::string name;
            int64_t cnt;
            ss >> name >> cnt;
            in_vertex = name == "vertex";
            if (in_vertex) n_vertex = cnt;
        } else if (tok == "property" && in_vertex) {
            std::string type, name;
            ss >> type >> name;
            if (type != "float" && type != "float32")
                throw std::runtime_error(path.string() + ": vertex property '" +
                                         name + "' is " + type + ", expected float");
            props.push_back(name);
        } else if (tok == "end_header") break;
    }
    if (!binary_le)
        throw std::runtime_error(path.string() + ": expected binary_little_endian");
    if (n_vertex < 0 || props.empty())
        throw std::runtime_error(path.string() + ": no vertex element");

    auto col = [&](const std::string& name) -> int {
        for (size_t i = 0; i < props.size(); ++i)
            if (props[i] == name) return (int)i;
        throw std::runtime_error(path.string() + ": missing property " + name);
    };
    const int cx = col("x"), cy = col("y"), cz = col("z");
    const int cr[4] = {col("rot_0"), col("rot_1"), col("rot_2"), col("rot_3")};
    const int cs[3] = {col("scale_0"), col("scale_1"), col("scale_2")};
    const int co = col("opacity");
    const int cd[3] = {col("f_dc_0"), col("f_dc_1"), col("f_dc_2")};

    SplatPly out;
    out.n = n_vertex;
    out.means.resize(n_vertex * 3);
    out.quats.resize(n_vertex * 4);
    out.scales.resize(n_vertex * 3);
    out.opacities.resize(n_vertex);
    out.features_dc.resize(n_vertex * 3);

    const size_t stride = props.size();
    const int64_t kRows = 1 << 16;
    std::vector<float> buf(stride * kRows);
    for (int64_t row = 0; row < n_vertex; row += kRows) {
        int64_t nr = std::min(kRows, n_vertex - row);
        f.read(reinterpret_cast<char*>(buf.data()),
               (std::streamsize)(nr * stride * sizeof(float)));
        if (!f) throw std::runtime_error(path.string() + ": truncated vertex data");
        for (int64_t i = 0; i < nr; ++i) {
            const float* r = buf.data() + i * stride;
            int64_t o = row + i;
            out.means[o*3+0] = r[cx]; out.means[o*3+1] = r[cy]; out.means[o*3+2] = r[cz];
            for (int k = 0; k < 4; ++k) out.quats[o*4+k] = r[cr[k]];
            for (int k = 0; k < 3; ++k) out.scales[o*3+k] = r[cs[k]];
            out.opacities[o] = r[co];
            for (int k = 0; k < 3; ++k) out.features_dc[o*3+k] = r[cd[k]];
        }
    }
    return out;
}

// ===========================================================================
// Checkpoint discovery (mirror of ss_meshing._find_ckpt + raw .ply support)
// ===========================================================================

// Returns (splat_ply_path, run_dir). run_dir is where config.json is searched.
std::pair<fs::path, fs::path> find_ckpt(const fs::path& checkpoint) {
    if (fs::is_regular_file(checkpoint) && checkpoint.extension() == ".ply")
        return {checkpoint, checkpoint.parent_path().parent_path()};
    if (fs::is_regular_file(checkpoint / "splat.ply"))
        return {checkpoint / "splat.ply", checkpoint.parent_path()};
    if (fs::is_directory(checkpoint)) {
        std::vector<fs::path> ckpts;
        for (const auto& e : fs::directory_iterator(checkpoint)) {
            std::string b = e.path().filename().string();
            if (e.is_directory() &&
                (b.rfind("step-", 0) == 0 || b.find(".ckpt") != std::string::npos))
                ckpts.push_back(e.path());
        }
        std::sort(ckpts.begin(), ckpts.end());
        for (auto it = ckpts.rbegin(); it != ckpts.rend(); ++it)
            if (fs::is_regular_file(*it / "splat.ply"))
                return {*it / "splat.ply", checkpoint};
    }
    throw std::runtime_error("no splat.ply found under " + checkpoint.string());
}

// ===========================================================================
// Camera loading via the CLI trainer's dataset parsers
// ===========================================================================

struct MeshCameras {
    std::vector<float>   positions;    // [C*3] camera centers, splat frame
    std::vector<float>   viewmats;     // [C*16] world->cam, engine convention
    std::vector<float>   intrins;      // [C*4]
    std::vector<float>   dist_coeffs;  // [C*10]
    std::vector<int32_t> widths, heights;
    std::string          model;
    int64_t num() const { return (int64_t)widths.size(); }
};

MeshCameras load_cameras(const JsonValue& run_cfg, const std::string& data_dir,
                         const std::string& data_format_override) {
    DatasetParserConfig pcfg;
    std::string data_format = data_format_override;
    const JsonValue* dp = run_cfg.find("dataparser");
    auto dp_str = [&](const char* key, const std::string& def) {
        if (!dp) return def;
        const JsonValue* v = dp->find(key);
        return (v && !v->is_null()) ? v->as_string() : def;
    };
    pcfg.recon_dir  = dp_str("colmap_recon_dir", "");
    pcfg.image_dir  = dp_str("image_dir", "images");
    pcfg.mask_dir   = dp_str("mask_dir", "masks");
    pcfg.metashape_xml = dp_str("metashape_xml", "");
    pcfg.metashape_ply = dp_str("metashape_ply", "");
    pcfg.metashape_psx = dp_str("metashape_psx", "");
    pcfg.downscale_rounding_mode = dp_str("downscale_rounding_mode", "floor");
    if (dp) {
        const JsonValue* v = dp->find("rescale_camera_to_fit");
        // bool auto-detect is unported; a number divides intrinsics
        if (v && v->type == JsonValue::Type::Number)
            pcfg.rescale_camera_to_fit = (float)v->as_double(0.0);
        const JsonValue* fmt = dp->find("data_format");
        if (data_format.empty() && fmt && !fmt->is_null()) data_format = fmt->as_string();
    }
    // Use every parsed frame (train + would-be eval): more cameras only make
    // the visual-hull occupancy and the bake better.
    pcfg.eval_mode = "all";
    pcfg.validation_fraction = 0.0f;

    ParsedDataset ds = parse_dataset(data_dir, pcfg, data_format);
    const int64_t C = ds.num_cameras;
    if (C == 0) throw std::runtime_error("dataset has no cameras: " + data_dir);

    // relative_scale: splats live in a frame scaled by this (model.py:561);
    // scale the c2w translations to match before inverting.
    float rel = 1.0f;
    if (const JsonValue* model = run_cfg.find("model")) {
        const JsonValue* rs = model->find("relative_scale");
        if (rs && rs->type == JsonValue::Type::Number) rel = (float)rs->as_double(1.0);
    }
    if (rel != 1.0f)
        for (int64_t i = 0; i < C; ++i)
            for (int r = 0; r < 3; ++r)
                ds.c2w[i*12 + r*4 + 3] *= rel;

    // Engine-convention world->cam (Y/Z flip + inverse) via the identity
    // (no-warp) post-split bake -- the same algebra the trainer uses.
    PostSplitCameras post = bake_post_split(ds, false, false);
    if (post.n_post != C)
        throw std::runtime_error("unexpected camera split in mesh CLI");

    MeshCameras out;
    out.viewmats = std::move(post.viewmats);
    out.intrins = std::move(post.intrins);
    out.dist_coeffs = std::move(post.dist_coeffs);
    out.widths = ds.widths;
    out.heights = ds.heights;
    out.positions.resize(C * 3);
    for (int64_t i = 0; i < C; ++i)
        for (int r = 0; r < 3; ++r)
            out.positions[i*3 + r] = ds.c2w[i*12 + r*4 + 3];

    int32_t m0 = ds.camera_models.empty() ? 0 : ds.camera_models[0];
    for (int32_t m : ds.camera_models)
        if (m != m0) {
            std::fprintf(stderr, "warning: mixed camera models in dataset; "
                         "using the first one for all cameras\n");
            break;
        }
    out.model = camera_model_to_string((CameraModelType)m0);
    return out;
}

// ===========================================================================
// Flags
// ===========================================================================

struct Options {
    std::string checkpoint;
    std::string data;
    bool no_data = false;
    std::string output;                 // base path (extension optional)
    std::string format = "ply";
    std::string color = "vertex";       // none | vertex | texture
    std::string data_format;            // "" = config/auto
    std::optional<float> iso;           // default depends on cameras
    meshing::MeshingConfig m;           // remaining knobs live here
};

void print_help(const char* argv0) {
    Options d;
    std::printf(
        "usage: %s <checkpoint> [--flag value ...]\n"
        "\n"
        "Extract a surface mesh from a trained 3DGS model.\n"
        "<checkpoint> is a run dir (with config.json + step-*.ckpt/), a *.ckpt\n"
        "dir, or a splat.ply file.\n"
        "\n"
        "  --data <dir>            dataset for camera-based occupancy/colors\n"
        "                          (default: config.json's data)\n"
        "  --no-data               mesh from Gaussian densities only\n"
        "  --data-format <fmt>     colmap|nerfstudio|metashape (default: auto)\n"
        "  --output <path>         output base path (default: <ckpt>/mesh);\n"
        "                          a known extension is stripped\n"
        "  --format <list>         comma-separated: ply,obj,gltf,glb [%s].\n"
        "                          With --color texture, a format may carry a\n"
        "                          texture encoding: glb+png (default),\n"
        "                          glb+jpg (JPEG q95), glb+jpeg75 (JPEG q75)\n"
        "  --color <mode>          none|vertex|texture [%s]\n"
        "                          (PLY+texture and OBJ+vertex are rejected)\n"
        "  --texture-size <n>      texture atlas resolution; 0 = auto from the\n"
        "                          observed-detail texel budget [%d]\n"
        "  --tex-gutter-px <n>     atlas spacing between UV charts [%d]\n"
        "  --chart-angle-deg <a>   max normal deviation within a UV chart [%g]\n"
        "\n"
        "  --iso <v>               isosurface level [0.5 with cameras, 0.2 without]\n"
        "  --merge-factor <v>      short-edge merge threshold multiplier [%g]\n"
        "  --bisection-iters <n>   bisection steps per cut edge [%d]\n"
        "  --max-cameras <n>       cap on cameras used, -1 = all [%d]\n"
        "  --max-grid-res <n>      acceleration grid cap [%d]\n"
        "  --grid-cell-factor <v>  grid cell size factor [%g]\n"
        "  --num-threads <n>       0 = all hardware threads [%d]\n"
        "  --carve-k <n>           k-th smallest occupancy over cameras [%d]\n"
        "  --cull-unseen <0|1>     drop mesh verts seen by no camera [%d]\n"
        "  --merge-max-flip-deg <a> fold guard for merge collapses [%g]\n"
        "  --floater-min-faces <n> drop components smaller than this [%d]\n"
        "  --fill-hole-ratio <v>   fill loops smaller than this x component [%g]\n"
        "  --fill-hole-max-edges <n> always fill loops with <= n edges [%d]\n"
        "  --degenerate-angle-deg <a> repair triangles with angles below [%g]\n"
        "  --quality-iters <n>     valence-flip + tangential-relax iterations [%d]\n"
        "  --quiet                 less progress output\n",
        argv0, d.format.c_str(), d.color.c_str(),
        d.m.texture_size, d.m.tex_gutter_px, (double)d.m.chart_angle_deg,
        (double)d.m.merge_factor, d.m.bisection_iters, d.m.max_cameras,
        d.m.max_grid_res, (double)d.m.grid_cell_factor, d.m.num_threads,
        d.m.carve_k, (int)d.m.cull_unseen, (double)d.m.merge_max_flip_deg,
        d.m.floater_min_faces, (double)d.m.fill_hole_ratio,
        d.m.fill_hole_max_edges, (double)d.m.degenerate_angle_deg,
        d.m.quality_iters);
}

float parse_f(const std::string& key, const char* v) {
    try { size_t p; float x = std::stof(v, &p); if (p == std::strlen(v)) return x; }
    catch (...) {}
    throw std::runtime_error("--" + key + ": expected a number, got '" + v + "'");
}
int parse_i(const std::string& key, const char* v) {
    try { size_t p; int x = std::stoi(v, &p); if (p == std::strlen(v)) return x; }
    catch (...) {}
    throw std::runtime_error("--" + key + ": expected an integer, got '" + v + "'");
}

Options parse_args(int argc, char** argv) {
    Options o;
    // meshing preset defaults match ss_meshing.py's argparse defaults
    o.m.floater_min_faces = 10;
    o.m.fill_hole_max_edges = 20;
    o.m.degenerate_angle_deg = 15.0f;
    o.m.max_cameras = -1;

    int i = 1;
    if (i < argc && argv[i][0] != '-') o.checkpoint = argv[i++];
    for (; i < argc; ++i) {
        std::string key = argv[i];
        if (key == "--help" || key == "-h") { print_help(argv[0]); std::exit(0); }
        while (!key.empty() && key[0] == '-') key.erase(key.begin());
        for (auto& c : key) if (c == '-') c = '_';
        auto next = [&]() -> const char* {
            if (i + 1 >= argc)
                throw std::runtime_error("--" + key + ": missing value");
            return argv[++i];
        };
        if      (key == "data")            o.data = next();
        else if (key == "no_data")         o.no_data = true;
        else if (key == "data_format")     o.data_format = next();
        else if (key == "output")          o.output = next();
        else if (key == "format")          o.format = next();
        else if (key == "color")           o.color = next();
        else if (key == "texture_size")    o.m.texture_size = parse_i(key, next());
        else if (key == "tex_gutter_px")   o.m.tex_gutter_px = parse_i(key, next());
        else if (key == "chart_angle_deg") o.m.chart_angle_deg = parse_f(key, next());
        else if (key == "iso")             o.iso = parse_f(key, next());
        else if (key == "merge_factor")    o.m.merge_factor = parse_f(key, next());
        else if (key == "bisection_iters") o.m.bisection_iters = parse_i(key, next());
        else if (key == "max_cameras")     o.m.max_cameras = parse_i(key, next());
        else if (key == "max_grid_res")    o.m.max_grid_res = parse_i(key, next());
        else if (key == "grid_cell_factor") o.m.grid_cell_factor = parse_f(key, next());
        else if (key == "num_threads")     o.m.num_threads = parse_i(key, next());
        else if (key == "carve_k")         o.m.carve_k = parse_i(key, next());
        else if (key == "cull_unseen")     o.m.cull_unseen = parse_i(key, next()) != 0;
        else if (key == "no_cull_unseen")  o.m.cull_unseen = false;
        else if (key == "merge_max_flip_deg") o.m.merge_max_flip_deg = parse_f(key, next());
        else if (key == "floater_min_faces")  o.m.floater_min_faces = parse_i(key, next());
        else if (key == "fill_hole_ratio")    o.m.fill_hole_ratio = parse_f(key, next());
        else if (key == "fill_hole_max_edges") o.m.fill_hole_max_edges = parse_i(key, next());
        else if (key == "degenerate_angle_deg") o.m.degenerate_angle_deg = parse_f(key, next());
        else if (key == "quality_iters")   o.m.quality_iters = parse_i(key, next());
        else if (key == "quiet")           o.m.verbose = false;
        else throw std::runtime_error("unknown flag: --" + key + " (see --help)");
    }
    if (o.checkpoint.empty())
        throw std::runtime_error("missing <checkpoint> argument (see --help)");
    return o;
}

}  // namespace


int spirula_mesh_main(int argc, char** argv) {
    try {
        Options o = parse_args(argc, argv);

        // ---- resolve color mode + formats, and validate BEFORE any work ----
        meshing::MeshColorMode mode;
        if      (o.color == "none")    mode = meshing::MeshColorMode::None;
        else if (o.color == "vertex")  mode = meshing::MeshColorMode::Vertex;
        else if (o.color == "texture") mode = meshing::MeshColorMode::Texture;
        else throw std::runtime_error("--color: expected none/vertex/texture, got '"
                                      + o.color + "'");
        o.m.color_mode = mode;
        std::vector<meshing::MeshFormatSpec> specs =
            meshing::parse_mesh_formats(o.format);
        if (specs.empty())
            throw std::runtime_error("--format: no valid format given");
        o.m.formats.clear();
        for (const auto& spec : specs) {
            std::string err = meshing::check_export_support(spec, mode);
            if (!err.empty()) throw std::runtime_error(err);
            o.m.formats.push_back(spec.token());
        }

        // ---- checkpoint ----
        auto [splat_ply, run_dir] = find_ckpt(o.checkpoint);
        std::printf("[meshing] loading %s\n", splat_ply.string().c_str());
        SplatPly splats = load_splat_ply(splat_ply);
        std::printf("[meshing] %lld Gaussians\n", (long long)splats.n);

        JsonValue run_cfg;
        fs::path cfg_path = run_dir / "config.json";
        if (fs::is_regular_file(cfg_path))
            run_cfg = json_parse_file(cfg_path.string());

        // ---- cameras ----
        MeshCameras cams;
        if (!o.no_data) {
            std::string data_dir = o.data;
            if (data_dir.empty()) {
                const JsonValue* d = run_cfg.find("data");
                if (d && !d->is_null()) {
                    fs::path cand = d->as_string();
                    if (cand.is_relative()) cand = run_dir / cand;
                    if (fs::exists(cand)) data_dir = cand.string();
                }
            }
            if (!data_dir.empty() && fs::exists(data_dir)) {
                std::printf("[meshing] loading cameras from %s\n", data_dir.c_str());
                cams = load_cameras(run_cfg, data_dir, o.data_format);
                std::printf("[meshing] %lld cameras (%dx%d, %s)\n",
                            (long long)cams.num(), cams.widths[0], cams.heights[0],
                            cams.model.c_str());
            }
        }
        const bool using_cameras = cams.num() > 0;
        if (!using_cameras)
            std::printf("\033[1;33m[meshing] No camera dataset in use: meshing "
                        "from Gaussian densities only.\n          Pass --data "
                        "<dataset_dir> for significantly better surfaces.\033[0m\n");

        o.m.iso = o.iso.value_or(using_cameras ? 0.5f : 0.2f);

        // ---- output base ----
        std::string out_base = o.output;
        if (out_base.empty())
            out_base = (splat_ply.parent_path() / "mesh").string();

        meshing::CameraParams cp;
        if (using_cameras) {
            cp.viewmats = cams.viewmats.data();
            cp.intrins = cams.intrins.data();
            cp.dist_coeffs = cams.dist_coeffs.data();
            cp.widths = cams.widths.data();
            cp.heights = cams.heights.data();
            cp.camera_model = cams.model;
        }
        bool ok = meshing::generate_mesh(
            splats.means.data(), splats.quats.data(), splats.scales.data(),
            splats.opacities.data(), splats.features_dc.data(), (int)splats.n,
            using_cameras ? cams.positions.data() : nullptr, (int)cams.num(),
            cp, o.m, out_base);
        if (!ok) return 1;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
