// ColmapParser.cpp -- see ColmapParser.h. Port of modules/colmap_utils.py
// (binary readers) + the COLMAP branch of modules/dataparser.py (frame
// assembly, aux-buffer discovery, pose conventions, normalization scale).

#include "ColmapParser.h"

#include "../Camera.h"   // camera_model_from_name

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <numeric>
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


// ===========================================================================
// Pose math
// ===========================================================================

using Mat3 = std::array<std::array<double, 3>, 3>;
using Vec3 = std::array<double, 3>;

// colmap_utils.py qvec2rotmat (world->camera rotation from (w,x,y,z)).
Mat3 qvec2rotmat(const std::array<double, 4>& q) {
    double w = q[0], x = q[1], y = q[2], z = q[3];
    return {{
        {1 - 2*y*y - 2*z*z, 2*x*y - 2*w*z,     2*x*z + 2*w*y},
        {2*x*y + 2*w*z,     1 - 2*x*x - 2*z*z, 2*y*z - 2*w*x},
        {2*x*z - 2*w*y,     2*y*z + 2*w*x,     1 - 2*x*x - 2*y*y},
    }};
}

// dataparser.py:644-649 -- COLMAP w2c -> nerfstudio/OpenGL c2w:
//   c2w[:3,:3] = R^T with columns 1, 2 negated (OpenCV -> OpenGL axis flip)
//   c2w[:3,3]  = -R^T @ t
struct C2W { Mat3 R; Vec3 t; };
C2W colmap_to_c2w(const ColmapImage& im) {
    Mat3 R = qvec2rotmat(im.qvec);
    C2W o;
    static const double flip[3] = {1.0, -1.0, -1.0};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            o.R[i][j] = R[j][i] * flip[j];
    for (int i = 0; i < 3; i++) {
        o.t[i] = 0.0;
        for (int j = 0; j < 3; j++) o.t[i] -= R[j][i] * im.tvec[j];
    }
    return o;
}

// Normalized-frame scale factor: auto_orient_and_center_poses(method="up",
// center_method="poses") followed by auto_scale_poses
// (dataparser.py:461-468). Only the resulting scalar matters for
// train_frame="points" -- poses/points themselves stay untouched -- so this
// implements just enough of camera_utils.auto_orient_and_center_poses to
// reproduce scale_factor:
//   up      = normalize(mean of c2w Y columns)
//   R_align = rotation taking `up` to +Z (Rodrigues)
//   center  = mean camera position
//   scale_factor = 1 / max |R_align @ (pos - center)|
// TODO: "pca" / "vertical" / "gsplat" orientation methods (camera_utils.py).
double compute_normalized_scale_factor(const std::vector<C2W>& poses) {
    Vec3 up{0, 0, 0}, center{0, 0, 0};
    for (const auto& p : poses) {
        for (int i = 0; i < 3; i++) up[i]     += p.R[i][1];
        for (int i = 0; i < 3; i++) center[i] += p.t[i];
    }
    double un = std::sqrt(up[0]*up[0] + up[1]*up[1] + up[2]*up[2]);
    for (auto& u : up) u /= std::max(un, 1e-12);
    for (auto& c : center) c /= (double)poses.size();

    // Rodrigues rotation aligning `up` with (0, 0, 1).
    Vec3 axis = {up[1], -up[0], 0.0};                  // up x z
    double s = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1]);
    double c = up[2];
    Mat3 R_align = {{{1,0,0},{0,1,0},{0,0,1}}};
    if (s > 1e-12) {
        for (auto& a : axis) a /= s;
        double C = 1.0 - c;
        R_align = {{
            {c + axis[0]*axis[0]*C,       axis[0]*axis[1]*C - axis[2]*s, axis[1]*s},
            {axis[0]*axis[1]*C + axis[2]*s, c + axis[1]*axis[1]*C,      -axis[0]*s},
            {-axis[1]*s,                   axis[0]*s,                    c},
        }};
    } else if (c < 0.0) {
        R_align = {{{1,0,0},{0,-1,0},{0,0,-1}}};       // up == -z: flip
    }

    double max_abs = 0.0;
    for (const auto& p : poses) {
        Vec3 d = {p.t[0] - center[0], p.t[1] - center[1], p.t[2] - center[2]};
        for (int i = 0; i < 3; i++) {
            double v = R_align[i][0]*d[0] + R_align[i][1]*d[1] + R_align[i][2]*d[2];
            max_abs = std::max(max_abs, std::abs(v));
        }
    }
    return 1.0 / std::max(max_abs, 1e-12);
}

// Auxiliary buffer discovery (dataparser.py _add_auxiliary_buffers, trimmed
// candidate list). Returns "" when no file matches.
std::string find_aux_file(const fs::path& aux_dir, const std::string& rel_name,
                          const char* suffix_tag) {
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

}  // namespace


// ===========================================================================
// parse_colmap_dataset
// ===========================================================================

ParsedDataset parse_colmap_dataset(const std::string& dataset_dir,
                                   const ColmapParserConfig& cfg) {
    // ---- Locate the reconstruction (dataparser.py:609-635) ---------------
    std::vector<std::string> probe;
    if (!cfg.recon_dir.empty()) probe = {cfg.recon_dir};
    else probe = {"sparse/0", "colmap/sparse/0", "sparse", "colmap", ""};

    std::string recon_dir;
    for (const auto& rel : probe) {
        fs::path d = fs::path(dataset_dir) / rel;
        if (fs::exists(d / "cameras.bin") && fs::exists(d / "images.bin") &&
            fs::exists(d / "points3D.bin")) {
            recon_dir = d.string();
            break;
        }
        // TODO: text-format fallback (cameras.txt / images.txt / points3D.txt),
        // port of colmap_utils.py read_*_text.
    }
    if (recon_dir.empty())
        throw std::runtime_error(
            "ColmapParser: no COLMAP reconstruction (cameras.bin/images.bin/"
            "points3D.bin) found under " + dataset_dir);

    auto cameras = read_cameras_binary(recon_dir);
    auto images  = read_images_binary(recon_dir);

    // ---- Assemble frames sorted by image filename (dataparser.py:300-316) -
    std::vector<const ColmapImage*> frames;
    frames.reserve(images.size());
    for (const auto& [id, im] : images) frames.push_back(&im);
    std::sort(frames.begin(), frames.end(),
              [](const ColmapImage* a, const ColmapImage* b) { return a->name < b->name; });

    // ---- Train/eval split: keep only the train subset (dataparser.py
    // eval_mode; get_train_eval_split_* at dataparser.py:51-129). ----------
    if (cfg.eval_mode != "all") {
        const int64_t n_all = (int64_t)frames.size();
        std::vector<const ColmapImage*> train;
        if (cfg.eval_mode == "fraction") {
            int64_t num_train = (int64_t)std::ceil((double)n_all * cfg.train_split_fraction);
            std::vector<char> keep(n_all, 0);
            for (int64_t k = 0; k < num_train; k++)
                keep[(num_train == 1) ? 0 : (int64_t)std::llround(
                    (double)k * (double)(n_all - 1) / (double)(num_train - 1))] = 1;
            for (int64_t i = 0; i < n_all; i++) if (keep[i]) train.push_back(frames[i]);
        } else if (cfg.eval_mode == "interval") {
            for (int64_t i = 0; i < n_all; i++)
                if (cfg.eval_interval <= 0 || i % cfg.eval_interval != 0)
                    train.push_back(frames[i]);
        } else if (cfg.eval_mode == "filename") {
            for (const auto* f : frames) {
                std::string base = fs::path(f->name).filename().string();
                if (base.find("train") != std::string::npos) train.push_back(f);
                else if (base.find("eval") == std::string::npos)
                    throw std::runtime_error(
                        "ColmapParser: eval_mode=filename requires 'train'/'eval' "
                        "in every image name; got " + f->name);
            }
        } else {
            throw std::runtime_error("ColmapParser: unknown eval_mode " + cfg.eval_mode);
        }
        if (train.empty())
            throw std::runtime_error("ColmapParser: eval_mode split left no training images");
        frames = std::move(train);
    }

    fs::path image_dir = fs::path(dataset_dir) / cfg.image_dir;

    ParsedDataset ds;
    const int64_t N = (int64_t)frames.size();
    ds.num_cameras = N;
    ds.camera_models.reserve(N);
    ds.image_filenames.reserve(N);
    ds.widths.reserve(N);
    ds.heights.reserve(N);
    ds.c2w.resize(N * 12);
    ds.viewmats.resize(N * 16);
    ds.intrins.resize(N * 4);
    ds.dist_coeffs.resize(N * 10);

    std::vector<std::string> mask_files(N), depth_files(N), normal_files(N);
    bool any_mask = false, any_depth = false, any_normal = false;

    std::vector<C2W> poses(N);

    for (int64_t i = 0; i < N; i++) {
        const ColmapImage& im = *frames[i];
        auto cam_it = cameras.find(im.camera_id);
        if (cam_it == cameras.end())
            throw std::runtime_error("ColmapParser: image " + im.name +
                                     " references missing camera id " +
                                     std::to_string(im.camera_id));
        const ColmapCamera& cam = cam_it->second;

        fs::path img_path = image_dir / im.name;
        if (!fs::exists(img_path))
            throw std::runtime_error("ColmapParser: " + img_path.string() +
                                     " does not exist (set --image-dir if needed)");
        ds.image_filenames.push_back(img_path.string());

        BakedIntrins bi = bake_colmap_intrins(cam);
        float W = (float)cam.width, H = (float)cam.height;
        if (cfg.rescale_camera_to_fit > 0.0f) {
            float s = cfg.rescale_camera_to_fit;
            bi.fx /= s; bi.fy /= s; bi.cx /= s; bi.cy /= s;
            W = std::floor(W / s); H = std::floor(H / s);   // TODO: rounding modes
        }
        ds.widths.push_back((int32_t)W);
        ds.heights.push_back((int32_t)H);
        ds.intrins[i*4 + 0] = bi.fx;
        ds.intrins[i*4 + 1] = bi.fy;
        ds.intrins[i*4 + 2] = bi.cx;
        ds.intrins[i*4 + 3] = bi.cy;
        std::copy(bi.dist.begin(), bi.dist.end(), ds.dist_coeffs.begin() + i*10);

        CameraModelType model = camera_model_from_name(bi.model_name);
        if ((int)model < 0)
            throw std::runtime_error("ColmapParser: unmapped camera model " + bi.model_name);
        ds.camera_models.push_back((int32_t)model);

        poses[i] = colmap_to_c2w(im);
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++)
                ds.c2w[i*12 + r*4 + c] = (float)poses[i].R[r][c];
            ds.c2w[i*12 + r*4 + 3] = (float)poses[i].t[r];
        }

        // Auxiliary supervision buffers, discovered by filename convention.
        mask_files[i]   = find_aux_file(fs::path(dataset_dir) / cfg.mask_dir,   im.name, "mask");
        depth_files[i]  = find_aux_file(fs::path(dataset_dir) / cfg.depth_dir,  im.name, "depth");
        normal_files[i] = find_aux_file(fs::path(dataset_dir) / cfg.normal_dir, im.name, "normal");
        any_mask   |= !mask_files[i].empty();
        any_depth  |= !depth_files[i].empty();
        any_normal |= !normal_files[i].empty();
    }
    if (any_mask)   ds.mask_filenames   = std::move(mask_files);
    if (any_depth)  ds.depth_filenames  = std::move(depth_files);
    if (any_normal) ds.normal_filenames = std::move(normal_files);

    // ---- Bake c2w -> engine viewmats (trainer.py:403-414) ----------------
    // R_v = c2w R with columns 1, 2 negated (OpenGL -> OpenCV); then
    // viewmat = [R_v^T | -R_v^T t; 0 0 0 1]. TODO: relative_scale on T.
    for (int64_t i = 0; i < N; i++) {
        static const double flip[3] = {1.0, -1.0, -1.0};
        double Rv[3][3];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                Rv[r][c] = poses[i].R[r][c] * flip[c];
        float* vm = ds.viewmats.data() + i*16;
        for (int r = 0; r < 3; r++) {
            double ti = 0.0;
            for (int c = 0; c < 3; c++) {
                vm[r*4 + c] = (float)Rv[c][r];          // R_v^T
                ti -= Rv[c][r] * poses[i].t[c];
            }
            vm[r*4 + 3] = (float)ti;
        }
        vm[12] = 0.f; vm[13] = 0.f; vm[14] = 0.f; vm[15] = 1.f;
    }

    // ---- train_frame_scale (train_frame="points": poses stay raw) --------
    double scale_factor = compute_normalized_scale_factor(poses);
    ds.train_frame_scale = (float)(scale_factor != 0.0 ? 1.0 / scale_factor : 1.0);

    // ---- Seed points ------------------------------------------------------
    ds.points = read_points3D_binary(recon_dir);

    // ---- Train / val split (eval_mode="all" + validation_fraction) -------
    // Port of get_train_eval_split_fraction with fraction =
    // 1 - validation_fraction: val indices are linspace-spread, train = all
    // (matching trainer.py:428-431, where train excludes the val set).
    std::vector<char> is_val(N, 0);
    if (cfg.validation_fraction > 0.0f && N > 1) {
        int64_t num_train = (int64_t)std::ceil((double)N * (1.0 - cfg.validation_fraction));
        std::vector<char> is_train(N, 0);
        for (int64_t k = 0; k < num_train; k++) {
            int64_t idx = (num_train == 1) ? 0
                : (int64_t)std::llround((double)k * (double)(N - 1) / (double)(num_train - 1));
            is_train[idx] = 1;
        }
        for (int64_t i = 0; i < N; i++) is_val[i] = !is_train[i];
    }
    for (int64_t i = 0; i < N; i++)
        (is_val[i] ? ds.val_indices : ds.train_indices).push_back((int32_t)i);

    return ds;
}
