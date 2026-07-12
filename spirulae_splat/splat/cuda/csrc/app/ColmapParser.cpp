// ColmapParser.cpp -- COLMAP-format reader for DatasetParser.h. Port of
// modules/colmap_utils.py (binary readers) + the COLMAP branch of
// modules/dataparser.py (frame assembly, pose conventions). Shared bake
// helpers live in DatasetCommon.cpp.

#include "DatasetParser.h"
#include "FastFloat.h"

#include "../CameraModel.h"   // camera_model_from_name (CUDA-free)

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <stdexcept>

namespace fs = std::filesystem;


// ===========================================================================
// Little-endian binary readers (colmap_utils.py read_next_bytes)
// ===========================================================================

namespace {

struct BinReader {
    FILE* f = nullptr;
    std::string path;

    explicit BinReader(const std::string& p) : path(p) {
        f = std::fopen(p.c_str(), "rb");
        if (!f) throw std::runtime_error("ColmapParser: cannot open " + p);
    }
    ~BinReader() { if (f) std::fclose(f); }

    template <typename T>
    T read() {
        // COLMAP files are little-endian; so is every platform we target.
        T v;
        if (std::fread(&v, sizeof(T), 1, f) != 1)
            throw std::runtime_error("ColmapParser: truncated file " + path);
        return v;
    }
    void skip(size_t bytes) {
        if (std::fseek(f, (long)bytes, SEEK_CUR) != 0)
            throw std::runtime_error("ColmapParser: truncated file " + path);
    }
    std::string read_cstring() {
        std::string s;
        for (;;) {
            int c = std::fgetc(f);
            if (c == EOF) throw std::runtime_error("ColmapParser: truncated file " + path);
            if (c == '\0') break;
            s.push_back((char)c);
        }
        return s;
    }
};

// COLMAP model id -> (name, num_params). Matches colmap/src/base/camera_models.h.
struct ColmapModelInfo { const char* name; int num_params; };
const std::map<int, ColmapModelInfo>& colmap_model_table() {
    static const std::map<int, ColmapModelInfo> t = {
        {0,  {"SIMPLE_PINHOLE", 3}},
        {1,  {"PINHOLE", 4}},
        {2,  {"SIMPLE_RADIAL", 4}},
        {3,  {"RADIAL", 5}},
        {4,  {"OPENCV", 8}},
        {5,  {"OPENCV_FISHEYE", 8}},
        {6,  {"FULL_OPENCV", 12}},
        {7,  {"FOV", 5}},
        {8,  {"SIMPLE_RADIAL_FISHEYE", 4}},
        {9,  {"RADIAL_FISHEYE", 5}},
        {10, {"THIN_PRISM_FISHEYE", 12}},
    };
    return t;
}

}  // namespace


std::map<int32_t, ColmapCamera> read_cameras_binary(const std::string& recon_dir) {
    BinReader r(recon_dir + "/cameras.bin");
    std::map<int32_t, ColmapCamera> cameras;
    uint64_t n = r.read<uint64_t>();
    for (uint64_t i = 0; i < n; i++) {
        ColmapCamera cam;
        cam.camera_id = r.read<int32_t>();
        int model_id  = r.read<int32_t>();
        cam.width     = r.read<uint64_t>();
        cam.height    = r.read<uint64_t>();
        auto it = colmap_model_table().find(model_id);
        if (it == colmap_model_table().end())
            throw std::runtime_error("ColmapParser: unknown camera model id " +
                                     std::to_string(model_id));
        cam.model = it->second.name;
        cam.params.resize(it->second.num_params);
        for (auto& p : cam.params) p = r.read<double>();
        cameras[cam.camera_id] = std::move(cam);
    }
    return cameras;
}

std::map<int32_t, ColmapImage> read_images_binary(const std::string& recon_dir) {
    BinReader r(recon_dir + "/images.bin");
    std::map<int32_t, ColmapImage> images;
    uint64_t n = r.read<uint64_t>();
    for (uint64_t i = 0; i < n; i++) {
        ColmapImage im;
        im.image_id = r.read<int32_t>();
        for (auto& q : im.qvec) q = r.read<double>();
        for (auto& t : im.tvec) t = r.read<double>();
        im.camera_id = r.read<int32_t>();
        im.name = r.read_cstring();
        // Skip the 2D feature track (x, y, point3D_id) -- not needed.
        uint64_t num_points2D = r.read<uint64_t>();
        r.skip(num_points2D * (2 * sizeof(double) + sizeof(uint64_t)));
        images[im.image_id] = std::move(im);
    }
    return images;
}

ColmapPoints3D read_points3D_binary(const std::string& recon_dir) {
    BinReader r(recon_dir + "/points3D.bin");
    ColmapPoints3D pts;
    uint64_t n = r.read<uint64_t>();
    pts.xyz.reserve(n * 3);
    pts.rgb.reserve(n * 3);
    for (uint64_t i = 0; i < n; i++) {
        r.skip(sizeof(uint64_t));                       // point3D_id
        for (int k = 0; k < 3; k++) pts.xyz.push_back((float)r.read<double>());
        for (int k = 0; k < 3; k++) pts.rgb.push_back(r.read<uint8_t>());
        r.skip(sizeof(double));                         // reprojection error
        uint64_t track_len = r.read<uint64_t>();
        r.skip(track_len * 2 * sizeof(int32_t));        // (image_id, point2D_idx)
    }
    return pts;
}


// ===========================================================================
// Text-format readers (cameras.txt / images.txt / points3D.txt), port of
// colmap_utils.py read_*_text. '#' lines are comments; fields are
// whitespace-separated (image names therefore contain no spaces, matching
// COLMAP's own text reader); each image record is followed by one raw
// POINTS2D line (possibly empty) which is skipped.
// ===========================================================================

namespace {

struct TextReader {
    std::string data;
    size_t pos = 0;

    explicit TextReader(const std::string& p) {
        FILE* f = std::fopen(p.c_str(), "rb");
        if (!f) throw std::runtime_error("ColmapParser: cannot open " + p);
        std::fseek(f, 0, SEEK_END);
        long n = std::ftell(f);
        std::fseek(f, 0, SEEK_SET);
        data.resize(n > 0 ? (size_t)n : 0);
        size_t got = data.empty() ? 0 : std::fread(&data[0], 1, data.size(), f);
        std::fclose(f);
        if (got != data.size())
            throw std::runtime_error("ColmapParser: cannot read " + p);
    }
    // next non-comment, non-empty line; false at EOF
    bool next_line(const char** s, const char** e) {
        while (pos < data.size()) {
            size_t eol = data.find('\n', pos);
            if (eol == std::string::npos) eol = data.size();
            size_t b = pos, t = eol;
            pos = eol + 1;
            while (b < t && (data[b]==' '||data[b]=='\t'||data[b]=='\r')) b++;
            while (t > b && (data[t-1]=='\r'||data[t-1]==' '||data[t-1]=='\t')) t--;
            if (b >= t || data[b] == '#') continue;
            *s = data.c_str() + b;
            *e = data.c_str() + t;
            return true;
        }
        return false;
    }
    // consume exactly one raw line (an empty POINTS2D line still counts)
    void skip_raw_line() {
        size_t eol = data.find('\n', pos);
        pos = (eol == std::string::npos) ? data.size() : eol + 1;
    }
};

}  // namespace

std::map<int32_t, ColmapCamera> read_cameras_text(const std::string& recon_dir) {
    TextReader r(recon_dir + "/cameras.txt");
    std::map<int32_t, ColmapCamera> cameras;
    const char *s, *e;
    while (r.next_line(&s, &e)) {
        ColmapCamera cam;
        char* p;
        cam.camera_id = (int32_t)std::strtol(s, &p, 10);
        while (p < e && std::isspace((unsigned char)*p)) p++;
        const char* m0 = p;
        while (p < e && !std::isspace((unsigned char)*p)) p++;
        cam.model.assign(m0, (size_t)(p - m0));
        cam.width  = std::strtoull(p, &p, 10);
        cam.height = std::strtoull(p, &p, 10);
        for (;;) {
            char* q;
            double v = fast_strtod(p, &q);
            if (q == p || q > e) break;      // no more params on this line
            cam.params.push_back(v);
            p = q;
        }
        if (cam.model.empty())
            throw std::runtime_error("ColmapParser: malformed cameras.txt line");
        cameras[cam.camera_id] = std::move(cam);
    }
    return cameras;
}

std::map<int32_t, ColmapImage> read_images_text(const std::string& recon_dir) {
    TextReader r(recon_dir + "/images.txt");
    std::map<int32_t, ColmapImage> images;
    const char *s, *e;
    while (r.next_line(&s, &e)) {
        ColmapImage im;
        char* p;
        im.image_id = (int32_t)fast_strtol(s, &p);
        for (auto& q : im.qvec) q = fast_strtod(p, &p);
        for (auto& t : im.tvec) t = fast_strtod(p, &p);
        im.camera_id = (int32_t)fast_strtol(p, &p);
        while (p < e && std::isspace((unsigned char)*p)) p++;
        const char* n0 = p;
        while (p < e && !std::isspace((unsigned char)*p)) p++;
        im.name.assign(n0, (size_t)(p - n0));
        if (im.name.empty())
            throw std::runtime_error("ColmapParser: malformed images.txt line");
        r.skip_raw_line();                   // POINTS2D[] line (may be empty)
        images[im.image_id] = std::move(im);
    }
    return images;
}

ColmapPoints3D read_points3D_text(const std::string& recon_dir) {
    TextReader r(recon_dir + "/points3D.txt");
    ColmapPoints3D pts;
    const char *s, *e;
    while (r.next_line(&s, &e)) {
        char* p;
        fast_strtol(s, &p);                  // point3D_id
        for (int k = 0; k < 3; k++) pts.xyz.push_back((float)fast_strtod(p, &p));
        for (int k = 0; k < 3; k++) pts.rgb.push_back((uint8_t)fast_strtol(p, &p));
        // reprojection error + track: rest of line, skipped
    }
    return pts;
}


// ===========================================================================
// Camera-param normalization (port of parse_colmap_camera_params,
// colmap_utils.py:431) + engine dist_coeffs layout (dataparser.py
// DISTORTION_KEYS = k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2).
// ===========================================================================

namespace {

struct BakedIntrins {
    float fx, fy, cx, cy;
    std::array<float, 10> dist{};   // k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2
    std::string model_name;         // nerfstudio-normalized: OPENCV / OPENCV_FISHEYE / ...
};

BakedIntrins bake_colmap_intrins(const ColmapCamera& cam) {
    const auto& p = cam.params;
    auto P = [&](size_t i) { return (float)p[i]; };
    BakedIntrins o;

    if (cam.model == "SIMPLE_PINHOLE" || cam.model == "SIMPLE_RADIAL" ||
        cam.model == "RADIAL") {
        o.fx = o.fy = P(0); o.cx = P(1); o.cy = P(2);
        if (cam.model != "SIMPLE_PINHOLE") o.dist[0] = P(3);          // k1
        if (cam.model == "RADIAL")         o.dist[1] = P(4);          // k2
        o.model_name = "OPENCV";
    } else if (cam.model == "PINHOLE") {
        o.fx = P(0); o.fy = P(1); o.cx = P(2); o.cy = P(3);
        o.model_name = "OPENCV";
    } else if (cam.model == "OPENCV") {
        o.fx = P(0); o.fy = P(1); o.cx = P(2); o.cy = P(3);
        o.dist[0] = P(4); o.dist[1] = P(5);                           // k1 k2
        o.dist[4] = P(6); o.dist[5] = P(7);                           // p1 p2
        o.model_name = "OPENCV";
    } else if (cam.model == "FULL_OPENCV") {
        o.fx = P(0); o.fy = P(1); o.cx = P(2); o.cy = P(3);
        o.dist[0] = P(4); o.dist[1] = P(5);                           // k1 k2
        o.dist[4] = P(6); o.dist[5] = P(7);                           // p1 p2
        o.dist[2] = P(8); o.dist[3] = P(9);                           // k3 k4
        // k5/k6 have no slot in the engine's 10-coeff layout; the Python
        // dataparser drops them too (DISTORTION_KEYS has no k5/k6).
        o.model_name = "FULL_OPENCV";
    } else if (cam.model == "OPENCV_FISHEYE" || cam.model == "THIN_PRISM_FISHEYE") {
        o.fx = P(0); o.fy = P(1); o.cx = P(2); o.cy = P(3);
        o.dist[0] = P(4); o.dist[1] = P(5);                           // k1 k2
        if (cam.model == "OPENCV_FISHEYE") {
            o.dist[2] = P(6); o.dist[3] = P(7);                       // k3 k4
        } else {
            o.dist[4] = P(6); o.dist[5] = P(7);                       // p1 p2
            o.dist[2] = P(8); o.dist[3] = P(9);                       // k3 k4
            o.dist[6] = P(10); o.dist[7] = P(11);                     // sx1 sy1
        }
        o.model_name = "OPENCV_FISHEYE";
    } else if (cam.model == "SIMPLE_RADIAL_FISHEYE" || cam.model == "RADIAL_FISHEYE") {
        o.fx = o.fy = P(0); o.cx = P(1); o.cy = P(2);
        o.dist[0] = P(3);                                             // k1
        if (cam.model == "RADIAL_FISHEYE") o.dist[1] = P(4);          // k2
        o.model_name = "OPENCV_FISHEYE";
    } else {
        throw std::runtime_error("ColmapParser: unsupported camera model " + cam.model);
    }
    return o;
}

// colmap_utils.py qvec2rotmat (world->camera rotation from (w,x,y,z)).
void qvec2rotmat(const std::array<double, 4>& q, double R[3][3]) {
    double w = q[0], x = q[1], y = q[2], z = q[3];
    R[0][0] = 1 - 2*y*y - 2*z*z; R[0][1] = 2*x*y - 2*w*z;     R[0][2] = 2*x*z + 2*w*y;
    R[1][0] = 2*x*y + 2*w*z;     R[1][1] = 1 - 2*x*x - 2*z*z; R[1][2] = 2*y*z - 2*w*x;
    R[2][0] = 2*x*z - 2*w*y;     R[2][1] = 2*y*z + 2*w*x;     R[2][2] = 1 - 2*x*x - 2*y*y;
}

// dataparser.py:644-649 -- COLMAP w2c -> nerfstudio/OpenGL c2w:
//   c2w[:3,:3] = R^T with columns 1, 2 negated (OpenCV -> OpenGL axis flip)
//   c2w[:3,3]  = -R^T @ t
void colmap_to_c2w(const ColmapImage& im, float* out12) {
    double R[3][3];
    qvec2rotmat(im.qvec, R);
    static const double flip[3] = {1.0, -1.0, -1.0};
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++)
            out12[r*4 + c] = (float)(R[c][r] * flip[c]);
        double t = 0.0;
        for (int c = 0; c < 3; c++) t -= R[c][r] * im.tvec[c];
        out12[r*4 + 3] = (float)t;
    }
}

}  // namespace


// ===========================================================================
// parse_colmap_dataset
// ===========================================================================

ParsedDataset parse_colmap_dataset(const std::string& dataset_dir,
                                   const DatasetParserConfig& cfg) {
    // ---- Locate the reconstruction (dataparser.py:609-635) ---------------
    std::vector<std::string> probe;
    if (!cfg.recon_dir.empty()) probe = {cfg.recon_dir};
    else probe = {"sparse/0", "colmap/sparse/0", "sparse", "colmap", ""};

    // points3D is required only in strict (trainer) mode; the lenient viewer
    // accepts a cameras-only reconstruction (poses / frustums). Binary and
    // text (cameras.txt / images.txt / points3D.txt) formats both count.
    auto has_recon = [&](const fs::path& d, const char* ext) {
        return fs::exists(d / (std::string("cameras.") + ext)) &&
               fs::exists(d / (std::string("images.") + ext)) &&
               (fs::exists(d / (std::string("points3D.") + ext)) ||
                !cfg.require_image_files);
    };
    std::string recon_dir;
    bool text_format = false;
    for (const auto& rel : probe) {
        fs::path d = fs::path(dataset_dir) / rel;
        if (has_recon(d, "bin")) { recon_dir = d.string(); break; }
        if (has_recon(d, "txt")) { recon_dir = d.string(); text_format = true; break; }
    }
    if (recon_dir.empty())
        throw std::runtime_error(
            "ColmapParser: no COLMAP reconstruction (cameras/images/points3D"
            " .bin or .txt) found under " + dataset_dir);

    auto cameras = text_format ? read_cameras_text(recon_dir)
                               : read_cameras_binary(recon_dir);
    auto images  = text_format ? read_images_text(recon_dir)
                               : read_images_binary(recon_dir);

    // ---- Assemble frames sorted by image filename (dataparser.py:300-316) -
    std::vector<const ColmapImage*> frames;
    frames.reserve(images.size());
    for (const auto& [id, im] : images) frames.push_back(&im);
    std::sort(frames.begin(), frames.end(),
              [](const ColmapImage* a, const ColmapImage* b) { return a->name < b->name; });

    // ---- All-frame c2w (needed for outlier filter + train_frame_scale) ----
    int64_t n_all = (int64_t)frames.size();
    std::vector<float> c2w_all(n_all * 12);
    std::vector<double> positions(n_all * 3);
    for (int64_t i = 0; i < n_all; i++) {
        colmap_to_c2w(*frames[i], &c2w_all[i*12]);
        for (int r = 0; r < 3; r++) positions[i*3 + r] = c2w_all[i*12 + r*4 + 3];
    }

    // ---- Outlier rejection (dataparser.py:319-325) -------------------------
    {
        std::vector<char> keep = dsparse::outlier_keep_mask(
            positions, n_all, cfg.outlier_threshold);
        std::vector<const ColmapImage*> kept;
        std::vector<float> kept_c2w;
        for (int64_t i = 0; i < n_all; i++) {
            if (!keep[i]) continue;
            kept.push_back(frames[i]);
            kept_c2w.insert(kept_c2w.end(), &c2w_all[i*12], &c2w_all[i*12] + 12);
        }
        frames = std::move(kept);
        c2w_all = std::move(kept_c2w);
        n_all = (int64_t)frames.size();
    }

    // ---- train_frame_scale + viewer remap transform over ALL post-outlier
    // frames (train + eval, matching the Python dataparser, which splits
    // after normalization). No applied_transform on the COLMAP path, so
    // train_to_normalized = inv(T_n_from_camera). -----------------------------
    double T_n[16], T_inv[16];
    double scale_factor = dsparse::compute_normalized_transform(c2w_all, n_all, T_n);
    dsparse::invert_affine4x4(T_n, T_inv);
    float train_frame_scale = (float)(scale_factor != 0.0 ? 1.0 / scale_factor : 1.0);

    // ---- eval_mode train subset --------------------------------------------
    std::vector<std::string> names(n_all);
    for (int64_t i = 0; i < n_all; i++) names[i] = frames[i]->name;
    std::vector<int64_t> subset = dsparse::train_subset(n_all, names, cfg);

    fs::path image_dir = fs::path(dataset_dir) / cfg.image_dir;

    ParsedDataset ds;
    const int64_t N = (int64_t)subset.size();
    ds.num_cameras = N;
    ds.train_frame_scale = train_frame_scale;
    for (int k = 0; k < 16; k++) ds.train_to_normalized[k] = (float)T_inv[k];
    ds.camera_models.reserve(N);
    ds.image_filenames.reserve(N);
    ds.widths.reserve(N);
    ds.heights.reserve(N);
    ds.c2w.resize(N * 12);
    ds.intrins.resize(N * 4);
    ds.dist_coeffs.resize(N * 10);

    std::vector<std::string> mask_files(N), depth_files(N), normal_files(N);
    bool any_mask = false, any_depth = false, any_normal = false;

    for (int64_t j = 0; j < N; j++) {
        const int64_t i = subset[j];
        const ColmapImage& im = *frames[i];
        auto cam_it = cameras.find(im.camera_id);
        if (cam_it == cameras.end())
            throw std::runtime_error("ColmapParser: image " + im.name +
                                     " references missing camera id " +
                                     std::to_string(im.camera_id));
        const ColmapCamera& cam = cam_it->second;

        fs::path img_path = image_dir / im.name;
        if (cfg.require_image_files && !fs::exists(img_path))
            throw std::runtime_error("ColmapParser: " + img_path.string() +
                                     " does not exist (set --image-dir if needed)");
        ds.image_filenames.push_back(img_path.string());

        BakedIntrins bi = bake_colmap_intrins(cam);
        float W = (float)cam.width, H = (float)cam.height;
        if (cfg.rescale_camera_to_fit > 0.0f) {
            float s = cfg.rescale_camera_to_fit;
            bi.fx /= s; bi.fy /= s; bi.cx /= s; bi.cy /= s;
            auto round_dim = [&](float v) {
                if (cfg.downscale_rounding_mode == "ceil")  return std::ceil(v / s);
                if (cfg.downscale_rounding_mode == "round") return std::round(v / s);
                return std::floor(v / s);
            };
            W = round_dim(W); H = round_dim(H);
        }
        ds.widths.push_back((int32_t)W);
        ds.heights.push_back((int32_t)H);
        ds.intrins[j*4 + 0] = bi.fx;
        ds.intrins[j*4 + 1] = bi.fy;
        ds.intrins[j*4 + 2] = bi.cx;
        ds.intrins[j*4 + 3] = bi.cy;
        std::copy(bi.dist.begin(), bi.dist.end(), ds.dist_coeffs.begin() + j*10);

        CameraModelType model = camera_model_from_name(bi.model_name);
        if ((int)model < 0)
            throw std::runtime_error("ColmapParser: unmapped camera model " + bi.model_name);
        ds.camera_models.push_back((int32_t)model);

        std::copy(&c2w_all[i*12], &c2w_all[i*12] + 12, &ds.c2w[j*12]);

        // Auxiliary supervision buffers, discovered by filename convention.
        mask_files[j]   = dsparse::find_aux_file(
            (fs::path(dataset_dir) / cfg.mask_dir).string(),   im.name, "mask");
        depth_files[j]  = dsparse::find_aux_file(
            (fs::path(dataset_dir) / cfg.depth_dir).string(),  im.name, "depth");
        normal_files[j] = dsparse::find_aux_file(
            (fs::path(dataset_dir) / cfg.normal_dir).string(), im.name, "normal");
        any_mask   |= !mask_files[j].empty();
        any_depth  |= !depth_files[j].empty();
        any_normal |= !normal_files[j].empty();
    }
    if (any_mask)   ds.mask_filenames   = std::move(mask_files);
    if (any_depth)  ds.depth_filenames  = std::move(depth_files);
    if (any_normal) ds.normal_filenames = std::move(normal_files);

    // In lenient (viewer) mode, tolerate a missing points3D file so a
    // cameras-only reconstruction still yields camera poses / frustums.
    if (fs::exists(fs::path(recon_dir) / "points3D.bin"))
        ds.points = read_points3D_binary(recon_dir);
    else if (fs::exists(fs::path(recon_dir) / "points3D.txt"))
        ds.points = read_points3D_text(recon_dir);
    else if (cfg.require_image_files)
        ds.points = read_points3D_binary(recon_dir);   // throws "cannot open"
    dsparse::assign_val_split(ds, cfg.validation_fraction);
    return ds;
}


// ===========================================================================
// Format dispatch / auto-detect (dataparser.py parse():246-269, same probe
// order: nerfstudio first, then COLMAP, then Metashape).
// ===========================================================================

ParsedDataset parse_dataset(const std::string& dataset_dir,
                            const DatasetParserConfig& cfg,
                            const std::string& format) {
    if (format == "colmap")     return parse_colmap_dataset(dataset_dir, cfg);
    if (format == "nerfstudio") return parse_nerfstudio_dataset(dataset_dir, cfg);
    if (format == "metashape")  return parse_metashape_dataset(dataset_dir, cfg);
    if (!format.empty())
        throw std::runtime_error("unsupported data format: " + format);
    if (fs::exists(fs::path(dataset_dir) / "transforms.json"))
        return parse_nerfstudio_dataset(dataset_dir, cfg);
    try {
        return parse_colmap_dataset(dataset_dir, cfg);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "Failed to parse COLMAP data: %s\n", e.what());
        std::fprintf(stderr, "Attempting to parse Metashape data...\n");
    }
    return parse_metashape_dataset(dataset_dir, cfg);
}
