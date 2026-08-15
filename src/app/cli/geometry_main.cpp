// `spirula geometry` -- depth and surface normals for a dataset. One
// monocular network per image; all the work around it is cameras, and that
// part is app/GeometryWarp.h.
//
// The output goes to `normals/` and `depths/` beside the images, which both
// dataset readers already probe by name. transforms.json is NOT rewritten: a
// tool that edits a user's dataset description can corrupt it.

#include "app/Tools.h"

#include "app/DepthPng.h"
#include "app/GeometryWarp.h"
#include "app/WriterPool.h"
#include "data/CameraMath.h"
#include "data/DatasetParser.h"
#include "i18n/Locale.h"
#include "i18n/catalog/Geometry.h"
#include "metric3d/Metric3D.h"
#include "metric3d/model/Fetch.h"
#include "nn/core/Error.h"
#include "nn/core/Log.h"
#include "nn/core/Parallel.h"
#include "nn/io/Image.h"
#include "nn/vk/Context.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace G = spirula::i18n::msg::geometry;

using spirula::i18n::format;

namespace {

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

enum class Tri { Auto, Yes, No };

struct Options {
    std::string dataset;
    std::string model = "metric3d-vit-large";
    // Metric3D's own inference size, 2.5x quicker than 1600 and no worse
    // (metric3d/README.md). A ceiling on one face, not on the frame.
    int  max_size = 1064;
    bool want_depth = false;
    bool want_normal = true;
    std::string normal_format = "png";
    int  jpeg_quality = 95;
    bool depth_mm = false;
    Tri  ray_depth = Tri::Auto;
    Tri  split = Tri::Auto;
    bool overwrite = false;
};

void help_row(const char* flags, const spirula::i18n::Msg& m, int col = 26) {
    std::string left = std::string("    ") + flags;
    if (spirula::i18n::display_width(left) >= col + 4) {
        std::fprintf(stderr, "%s\n", left.c_str());
        left.clear();
    }
    left = spirula::i18n::pad_to(left, col + 4);
    for (const std::string& line : spirula::i18n::wrap(m.get(), 86 - col - 4)) {
        std::fprintf(stderr, "%s%s\n", left.c_str(), line.c_str());
        left.assign((size_t)col + 4, ' ');
    }
}

void usage() {
    const std::string prog = app::program_name();
    std::fprintf(stderr, "%s -- %s\n\n", prog.c_str(), G::tagline.get());
    std::fprintf(stderr, "    %s <dataset> [options]\n\n", prog.c_str());
    for (const std::string& l : spirula::i18n::wrap(G::usage_target.get(), 80))
        std::fprintf(stderr, "    %s\n", l.c_str());
    std::fprintf(stderr, "\n%s\n", G::head_options.get());
    help_row("--model <id|file>", G::opt_model);
    help_row("--max-size <n>", G::opt_max_size);
    help_row("--depth", G::opt_depth);
    help_row("--normal-format png|jpg", G::opt_normal_format);
    help_row("--jpeg-quality <n>", G::opt_jpeg_quality);
    help_row("--depth-units relative|mm", G::opt_depth_units);
    help_row("--ray-depth auto|yes|no", G::opt_ray_depth);
    help_row("--split auto|yes|no", G::opt_split);
    help_row("--overwrite", G::opt_overwrite);
    // English, like the other deep diagnostics in this repository: what it
    // prints is a table of numerical errors, read by whoever changed the warp.
    std::fprintf(stderr, "    --check                   "
                         "run the camera round-trip self-test and exit\n");
    std::fprintf(stderr, "\n%s --device <index|name>  --lang <code>\n",
                 G::label_common.get());
    std::fprintf(stderr, "%s SS_NN_LOG=0..3  SS_VK_DEVICE  SS_PROFILE=1\n",
                 G::label_environment.get());
    std::fprintf(stderr, "\n    %s\n", metric3d::model_id_list().c_str());
}

bool tri(const char* v, Tri& out) {
    if (!std::strcmp(v, "auto")) out = Tri::Auto;
    else if (!std::strcmp(v, "yes") || !std::strcmp(v, "1")) out = Tri::Yes;
    else if (!std::strcmp(v, "no") || !std::strcmp(v, "0")) out = Tri::No;
    else return false;
    return true;
}

// ---------------------------------------------------------------------------
// Host image handling
// ---------------------------------------------------------------------------

// Area average, not the bilinear one next door in DataManager: at a 2x or 3x
// downscale bilinear reads two of every three source pixels, and aliases the
// fine structure the normals are made of.
std::vector<float> resize_area(const nn::Image& img, int dw, int dh) {
    std::vector<float> out((size_t)dw * dh * 3, 0.0f);
    nn::parallel_for(dh, [&](int64_t y0, int64_t y1) {
        for (int64_t y = y0; y < y1; ++y) {
            const int sy0 = (int)((double)y * img.height / dh);
            const int sy1 = std::max(sy0 + 1, (int)((double)(y + 1) * img.height / dh));
            for (int x = 0; x < dw; ++x) {
                const int sx0 = (int)((double)x * img.width / dw);
                const int sx1 = std::max(sx0 + 1, (int)((double)(x + 1) * img.width / dw));
                double acc[3] = {0, 0, 0};
                for (int sy = sy0; sy < sy1; ++sy)
                    for (int sx = sx0; sx < sx1; ++sx) {
                        const uint8_t* p = &img.data[((size_t)sy * img.width + sx) * 3];
                        for (int c = 0; c < 3; ++c) acc[c] += p[c];
                    }
                const double n = (double)(sy1 - sy0) * (sx1 - sx0) * 255.0;
                for (int c = 0; c < 3; ++c)
                    out[((size_t)y * dw + x) * 3 + c] = (float)(acc[c] / n);
            }
        }
    });
    return out;
}

// The scale a depth map is stored against. Not the maximum: the network
// clamps its own output at 200, and a few sky pixels there would leave a whole
// room inside 1% of the 16-bit range.
float depth_scale(const std::vector<float>& depth) {
    std::vector<float> v;
    v.reserve(depth.size());
    for (float d : depth)
        if (d > 0.0f) v.push_back(d);
    if (v.empty()) return 1.0f;
    const size_t at = (size_t)((double)(v.size() - 1) * 0.999);
    std::nth_element(v.begin(), v.begin() + (long)at, v.end());
    return std::fmax(v[at], 1e-6f);
}

std::string human_time(double ms) {
    char buf[64];
    if (ms < 60000) std::snprintf(buf, sizeof buf, "%.0fs", ms / 1000.0);
    else if (ms < 3600000) std::snprintf(buf, sizeof buf, "%.0fm", ms / 60000.0);
    else std::snprintf(buf, sizeof buf, "%.1fh", ms / 3600000.0);
    return buf;
}

// ---------------------------------------------------------------------------
// --check: the camera round trip, with no network in it
//
// English, like the other deep diagnostics in this repository. It is the only
// thing that tests the face rotation and the depth conversion: put an analytic
// PLANE in front of the camera, hand each face the depth and normal that plane
// would produce in ITS frame, and check what comes back in the source camera's
// frame is the plane again. A wrong rotation, a wrong z-vs-ray scale or a
// transposed map all fail here, and none of them is visible in a normal map
// that merely looks plausible.
// ---------------------------------------------------------------------------

int self_check() {
    struct Case {
        const char* name;
        int model, tier;
        int w, h;
        float fx, fy, cx, cy;
        float dist[kCameraDistortionParams];
        bool split;
    };
    const Case cases[] = {
        {"pinhole + opencv", (int)CameraModelType::PINHOLE,
         (int)CameraDistortionType::OpenCV, 1600, 1200, 1500, 1500, 800, 600,
         {-0.03f, 0.01f, 0.001f, -0.002f}, false},
        {"fisheye + prism, one face", (int)CameraModelType::FISHEYE,
         (int)CameraDistortionType::ThinPrism, 1000, 1500, 530, 530, 500, 750,
         {0.009f, 0.002f, -0.0006f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, false},
        {"fisheye + prism, split", (int)CameraModelType::FISHEYE,
         (int)CameraDistortionType::ThinPrism, 1000, 1500, 530, 530, 500, 750,
         {0.009f, 0.002f, -0.0006f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}, true},
        {"equisolid, split", (int)CameraModelType::EQUISOLID,
         (int)CameraDistortionType::None, 1400, 1400, 420, 420, 700, 700, {}, true},
        {"equirectangular, cube map", (int)CameraModelType::EQUIRECTANGULAR,
         (int)CameraDistortionType::None, 1024, 512,
         1024.0f / 6.28318530718f, 1024.0f / 6.28318530718f, 512, 256, {}, true},
    };
    // Three planes, none of them axis-aligned, so a swapped axis cannot pass.
    const double planes[3][3] = {{0.0, 0.0, 1.0}, {0.35, -0.25, 0.90},
                                 {-0.5, 0.4, 0.77}};
    const double dist_to_plane = 4.0;

    int failures = 0;
    for (const Case& c : cases) {
        app::GeometryCamera cam;
        cam.model = c.model;
        cam.distortion = c.tier;
        cam.width = c.w;
        cam.height = c.h;
        cam.fx = c.fx; cam.fy = c.fy; cam.cx = c.cx; cam.cy = c.cy;
        std::copy(std::begin(c.dist), std::end(c.dist), std::begin(cam.dist));

        app::GeometryWarp warp;
        const int patch = metric3d::Predictor::sizeGranularity();
        // Uncapped faces on purpose: what is left after a correct rotation is
        // bilinear interpolation of a curve, and that error falls with the
        // square of the face's pixels.
        warp.plan(cam, c.w / patch * patch, c.h / patch * patch, c.split, patch, 0);

        const double* axes = warp.faceAxes();
        // The plane as each face would see it: `fd` its own z-depth, `fn` its
        // own frame's normal.
        auto face_gt = [&](const double nn3[3], std::vector<std::vector<float>>& fd,
                           std::vector<std::vector<float>>& fn) {
            fd.assign((size_t)warp.faces(), {});
            fn.assign((size_t)warp.faces(), {});
            for (int k = 0; k < warp.faces(); ++k) {
                fd[k].assign((size_t)warp.faceWidth() * warp.faceHeight(), 0.0f);
                fn[k].assign(fd[k].size() * 3, 0.0f);
                for (int y = 0; y < warp.faceHeight(); ++y)
                    for (int x = 0; x < warp.faceWidth(); ++x) {
                        const double tx = -1.0 + 2.0 * (x + 0.5) / warp.faceWidth();
                        const double ty = -1.0 + 2.0 * (y + 0.5) / warp.faceHeight();
                        double r[3], ax[3], ay[3], az[3];
                        if (axes) {
                            const double* a = axes + (size_t)k * 9;
                            for (int m = 0; m < 3; ++m) {
                                ax[m] = a[m]; ay[m] = a[3+m]; az[m] = a[6+m];
                                r[m] = az[m] + tx * ax[m] + ty * ay[m];
                            }
                            // The extent lives in the axis lengths; the frame
                            // a normal is written in does not.
                            const double lx = std::sqrt(ax[0]*ax[0]+ax[1]*ax[1]+ax[2]*ax[2]);
                            const double ly = std::sqrt(ay[0]*ay[0]+ay[1]*ay[1]+ay[2]*ay[2]);
                            for (int m = 0; m < 3; ++m) { ax[m] /= lx; ay[m] /= ly; }
                        } else {
                            // One face: its frame is the camera's, and tx, ty
                            // are not what indexes it -- the face pixel is the
                            // camera's own normalized coordinate.
                            const double fx = cam.fx * (double)warp.outWidth() / c.w;
                            const double fy = cam.fy * (double)warp.outHeight() / c.h;
                            const double cx = cam.cx * (double)warp.outWidth() / c.w;
                            const double cy = cam.cy * (double)warp.outHeight() / c.h;
                            r[0] = (x + 0.5 - cx) / fx;
                            r[1] = (y + 0.5 - cy) / fy;
                            r[2] = 1.0;
                            ax[0]=1; ax[1]=0; ax[2]=0;
                            ay[0]=0; ay[1]=1; ay[2]=0;
                            az[0]=0; az[1]=0; az[2]=1;
                        }
                        const double dot = nn3[0]*r[0] + nn3[1]*r[1] + nn3[2]*r[2];
                        const size_t i = (size_t)y * warp.faceWidth() + x;
                        // The plane at face z-depth t, where the face's own z
                        // is r . az and r was built with r . az == 1.
                        fd[k][i] = (float)(dist_to_plane / std::fmax(dot, 0.05));
                        fn[k][i * 3 + 0] =
                            (float)(nn3[0]*ax[0] + nn3[1]*ax[1] + nn3[2]*ax[2]);
                        fn[k][i * 3 + 1] =
                            (float)(nn3[0]*ay[0] + nn3[1]*ay[1] + nn3[2]*ay[2]);
                        fn[k][i * 3 + 2] =
                            (float)(nn3[0]*az[0] + nn3[1]*az[1] + nn3[2]*az[2]);
                    }
            }
        };

        double worst_n = 0, worst_z = 0, worst_r = 0, worst_align = 0;
        for (const double(&n)[3] : planes) {
            const double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
            const double nn3[3] = {n[0]/len, n[1]/len, n[2]/len};

            std::vector<std::vector<float>> fd, fn;
            face_gt(nn3, fd, fn);

            std::vector<float> dz, dr, nrm;
            warp.gather(fd, fn, false, &dz, &nrm);
            warp.gather(fd, fn, true, &dr, nullptr);

            const double sx = (double)warp.outWidth() / c.w;
            const double sy = (double)warp.outHeight() / c.h;
            const double fx = cam.fx * sx, fy = cam.fy * sy;
            const double cx = cam.cx * sx, cy = cam.cy * sy;
            for (int y = 0; y < warp.outHeight(); ++y)
                for (int x = 0; x < warp.outWidth(); ++x) {
                    const size_t i = (size_t)y * warp.outWidth() + x;
                    if (dz[i] == 0.0f) continue;
                    double r[3];
                    if (!camhost::generate_ray((x + 0.5 - cx) / fx, (y + 0.5 - cy) / fy,
                                               cam.model, cam.distortion, cam.dist, r))
                        continue;
                    const double dot = nn3[0]*r[0] + nn3[1]*r[1] + nn3[2]*r[2];
                    if (dot < 0.25) continue;   // grazing: the plane is edge-on
                    const double want_ray = dist_to_plane / dot;
                    const double want_z = want_ray * r[2];
                    worst_r = std::fmax(worst_r, std::fabs(dr[i] - want_ray) / want_ray);
                    worst_z = std::fmax(worst_z, std::fabs(dz[i] - want_z) /
                                                     std::fmax(std::fabs(want_z), 1e-6));
                    const double c_ang = nrm[i*3+0]*nn3[0] + nrm[i*3+1]*nn3[1] +
                                         nrm[i*3+2]*nn3[2];
                    worst_n = std::fmax(worst_n,
                                        std::acos(std::fmin(1.0, std::fmax(-1.0, c_ang))));
                }

            // Hand the faces that same plane at deliberately different
            // scales and the blend should put them back on one. The gauge is
            // free, so this measures the ratio's SPREAD, not the ratio.
            if (!warp.split()) continue;
            for (int k = 0; k < warp.faces(); ++k) {
                const float s = (float)std::exp(0.25 * (k - warp.faces() * 0.5));
                for (float& v : fd[(size_t)k]) v *= s;
            }
            std::vector<float> misfit;
            warp.gather(fd, {}, true, &dr, nullptr);
            for (int y = 0; y < warp.outHeight(); ++y)
                for (int x = 0; x < warp.outWidth(); ++x) {
                    const size_t i = (size_t)y * warp.outWidth() + x;
                    if (dr[i] == 0.0f) continue;
                    double r[3];
                    if (!camhost::generate_ray((x + 0.5 - cx) / fx, (y + 0.5 - cy) / fy,
                                               cam.model, cam.distortion, cam.dist, r))
                        continue;
                    const double dot = nn3[0]*r[0] + nn3[1]*r[1] + nn3[2]*r[2];
                    if (dot < 0.25) continue;
                    misfit.push_back((float)std::log(dr[i] * dot / dist_to_plane));
                }
            if (misfit.size() > 100) {
                std::sort(misfit.begin(), misfit.end());
                worst_align = std::fmax(
                    worst_align,
                    misfit[misfit.size() * 99 / 100] - misfit[misfit.size() / 100]);
            }
        }
        // Bilinear resampling of a smooth but curved field is what is left; a
        // rotation or a scale error is orders of magnitude above these.
        const bool ok = worst_n < 0.02 && worst_z < 0.02 && worst_r < 0.02 &&
                        worst_align < 0.05;
        if (!ok) ++failures;
        std::printf("  %s %-30s faces %d  normal %.4f deg  z %.2e  ray %.2e"
                    "  realign %.3f\n",
                    ok ? "ok  " : "FAIL", c.name, warp.faces(),
                    worst_n * 57.2957795, worst_z, worst_r, worst_align);
    }

    // The closed form that decided how many faces there were, against
    // counting the pixels it integrates -- of the ideal lens, which is what it
    // claims to measure (data/CameraMath.h).
    std::printf("\n  --split auto / --input-depth-is-ray-depth threshold at 0.75\n");
    for (const Case& c : cases) {
        const float zero[kCameraDistortionParams] = {};
        const int tier = (int)CameraDistortionType::None;
        const double got = camhost::pinhole_coverage(c.model, c.w, c.h, c.fx, c.fy);
        // Only the two angular models integrate anything; a perspective lens
        // is 1 by definition and a panorama is 0 by fiat.
        if (c.model != (int)CameraModelType::FISHEYE &&
            c.model != (int)CameraModelType::EQUISOLID) {
            std::printf("  ok   %-30s coverage %.4f  (by definition)   %s\n", c.name,
                        got, got <= 0.75 ? "split" : "undistort");
            continue;
        }
        int64_t inside = 0, total = 0;
        for (int y = 0; y < 512; ++y)
            for (int x = 0; x < 512; ++x) {
                const double u = ((x + 0.5) / 512 * c.w - c.w * 0.5) / c.fx;
                const double v = ((y + 0.5) / 512 * c.h - c.h * 0.5) / c.fy;
                double r[3];
                ++total;
                if (!camhost::generate_ray(u, v, c.model, tier, zero, r)) continue;
                if (r[2] <= 0) continue;
                if (std::fabs(r[0] / r[2]) * c.fx <= c.w * 0.5 &&
                    std::fabs(r[1] / r[2]) * c.fy <= c.h * 0.5)
                    ++inside;
            }
        const double want = (double)inside / (double)total;
        const bool ok = std::fabs(got - want) < 0.01;
        if (!ok) ++failures;
        std::printf("  %s %-30s coverage %.4f  (counted %.4f)  %s\n",
                    ok ? "ok  " : "FAIL", c.name, got, want,
                    got <= 0.75 ? "split" : "undistort");
    }

    std::printf("\n%s\n", failures ? "FAILURES" : "camera round trip is exact");
    return failures ? 1 : 0;
}

// Every image whose camera is byte-identical shares one warp plan; a dataset
// usually has one camera and always has few.
std::string camera_key(const ParsedDataset& ds, int64_t i) {
    std::string k;
    k.reserve(96);
    auto num = [&](double v) { k += std::to_string(v) + ","; };
    num(ds.camera_models[(size_t)i]);
    num(ds.camera_distortions[(size_t)i]);
    num(ds.widths[(size_t)i]);
    num(ds.heights[(size_t)i]);
    for (int j = 0; j < 4; ++j) num(ds.intrins[(size_t)i * 4 + j]);
    for (int j = 0; j < kCameraDistortionParams; ++j)
        num(ds.dist_coeffs[(size_t)i * kCameraDistortionParams + j]);
    if (!ds.redistort.empty()) {
        num(ds.redistort[(size_t)i].source_model);
        for (int j = 0; j < 16; ++j) num(ds.redistort[(size_t)i].params[j]);
    }
    return k;
}

}  // namespace

int spirula_geometry_main(int argc, char** argv) {
    app::set_program_name(argc > 0 ? argv[0] : nullptr, "spirula geometry");
    Options o;
    std::string device;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "%s needs a value\n", a.c_str());
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--help" || a == "-h") { usage(); return 0; }
        else if (a == "--check") return self_check();
        else if (a == "--model") o.model = next();
        else if (a == "--max-size") o.max_size = std::atoi(next());
        else if (a == "--depth") o.want_depth = true;
        else if (a == "--no-normal") o.want_normal = false;
        else if (a == "--normal-format") o.normal_format = next();
        else if (a == "--jpeg-quality") o.jpeg_quality = std::atoi(next());
        else if (a == "--depth-units") o.depth_mm = std::string(next()) == "mm";
        else if (a == "--ray-depth") { if (!tri(next(), o.ray_depth)) { usage(); return 2; } }
        else if (a == "--split") { if (!tri(next(), o.split)) { usage(); return 2; } }
        else if (a == "--overwrite") o.overwrite = true;
        else if (a == "--device") device = next();
        else if (!a.empty() && a[0] == '-') {
            std::fprintf(stderr, "unknown option '%s'\n\n", a.c_str());
            usage();
            return 2;
        } else if (o.dataset.empty()) o.dataset = a;
        else { usage(); return 2; }
    }
    if (o.dataset.empty()) {
        std::fprintf(stderr, "%s\n",
                     format(G::err_no_dataset, {app::program_name()}).c_str());
        return 2;
    }
    if (!o.want_depth && !o.want_normal) return 0;
    if (!device.empty()) {
        // Creating the context here fixes the device for the process; the
        // model's own Context::get() then returns this one.
        nn::vk::ContextOptions vo;
        char* end = nullptr;
        const long idx = std::strtol(device.c_str(), &end, 10);
        if (end && *end == '\0') vo.device_index = (int)idx;
        else vo.device_match = device;
        nn::vk::Context::get(vo);
    }

    try {
        // ---- the dataset --------------------------------------------------
        DatasetParserConfig cfg;
        cfg.require_image_files = false;   // the point cloud is not our business
        const ParsedDataset ds = parse_dataset(o.dataset, cfg, "");
        const int64_t N = ds.num_cameras;
        NN_CHECK(N > 0, "'%s' holds no images", o.dataset.c_str());

        std::vector<std::string> keys((size_t)N);
        std::vector<int> group((size_t)N, 0);
        std::vector<int64_t> group_first;
        {
            std::vector<std::string> seen;
            for (int64_t i = 0; i < N; ++i) {
                keys[(size_t)i] = camera_key(ds, i);
                auto it = std::find(seen.begin(), seen.end(), keys[(size_t)i]);
                if (it == seen.end()) {
                    group[(size_t)i] = (int)seen.size();
                    seen.push_back(keys[(size_t)i]);
                    group_first.push_back(i);
                } else {
                    group[(size_t)i] = (int)(it - seen.begin());
                }
            }
        }
        std::printf("%s\n", format(G::log_dataset,
                                   {fs::path(o.dataset).filename().string(),
                                    (long long)N, (long long)group_first.size()})
                                .c_str());

        // ---- one warp plan per camera --------------------------------------
        const int patch = metric3d::Predictor::sizeGranularity();
        std::vector<app::GeometryWarp> warps(group_first.size());
        for (size_t g = 0; g < group_first.size(); ++g) {
            const int64_t i = group_first[g];
            app::GeometryCamera cam;
            cam.model = ds.camera_models[(size_t)i];
            cam.distortion = ds.camera_distortions[(size_t)i];
            cam.width = ds.widths[(size_t)i];
            cam.height = ds.heights[(size_t)i];
            cam.fx = ds.intrins[(size_t)i * 4 + 0];
            cam.fy = ds.intrins[(size_t)i * 4 + 1];
            cam.cx = ds.intrins[(size_t)i * 4 + 2];
            cam.cy = ds.intrins[(size_t)i * 4 + 3];
            for (int j = 0; j < kCameraDistortionParams; ++j)
                cam.dist[j] = ds.dist_coeffs[(size_t)i * kCameraDistortionParams + j];
            if (!ds.redistort.empty() && ds.redistort[(size_t)i].source_model >= 0) {
                cam.source_model = ds.redistort[(size_t)i].source_model;
                std::copy(std::begin(ds.redistort[(size_t)i].params),
                          std::end(ds.redistort[(size_t)i].params),
                          std::begin(cam.source_params));
                NN_LOG_WARN("%s\n",
                            format(G::warn_fitted_camera, {(long long)g}).c_str());
            }

            // The output resolution: the source frame under the longest-side
            // cap, on the patch grid so the undistort path is not cropped.
            double s = 1.0;
            if (o.max_size > 0 && std::max(cam.width, cam.height) > o.max_size)
                s = (double)o.max_size / std::max(cam.width, cam.height);
            const int ow = std::max(patch, (int)(cam.width * s) / patch * patch);
            const int oh = std::max(patch, (int)(cam.height * s) / patch * patch);

            // The rule the trainer's --input-depth-is-ray-depth resolves with
            // too: split when one undistortion would throw away more than a
            // quarter of the frame (data/CameraMath.h).
            const bool split =
                o.split == Tri::Auto
                    ? camhost::pinhole_coverage(cam.model, cam.width, cam.height,
                                                cam.fx, cam.fy) <= 0.75
                    : o.split == Tri::Yes;
            warps[g].plan(cam, ow, oh, split, patch, o.max_size);

            const char* model_name = camera_model_to_string(
                (CameraModelType)cam.model);
            if (warps[g].split())
                std::printf("  %s\n", format(G::log_camera_split,
                                             {(long long)g, model_name,
                                              (long long)warps[g].faces(),
                                              (long long)warps[g].faceWidth(),
                                              (long long)warps[g].faceHeight()})
                                          .c_str());
            else
                std::printf("  %s\n", format(G::log_camera_single,
                                             {(long long)g, model_name,
                                              (long long)warps[g].faceWidth(),
                                              (long long)warps[g].faceHeight()})
                                          .c_str());
        }

        // ---- where the output goes -----------------------------------------
        const fs::path root(o.dataset);
        const fs::path image_root = root / cfg.image_dir;
        const fs::path normal_dir = root / cfg.normal_dir;
        const fs::path depth_dir = root / cfg.depth_dir;
        std::error_code ec;
        if (o.want_normal) fs::create_directories(normal_dir, ec);
        if (o.want_depth) fs::create_directories(depth_dir, ec);

        // The output mirrors whatever shape `images/` has, subdirectories and
        // all -- a multi-camera capture is `images/cam0/...`, and a writer
        // that only made the top folder failed every single frame of one.
        auto out_path = [&](const fs::path& dir, int64_t i, const char* ext) {
            fs::path rel = fs::relative(fs::path(ds.image_filenames[(size_t)i]),
                                        image_root, ec);
            if (ec || rel.empty() || rel.native()[0] == '.')
                rel = fs::path(ds.image_filenames[(size_t)i]).filename();
            rel.replace_extension(ext);
            const fs::path out = dir / rel;
            if (out.has_parent_path()) fs::create_directories(out.parent_path(), ec);
            return out;
        };

        // ---- the model -------------------------------------------------------
        metric3d::Predictor pred;
        pred.load(o.model);

        app::WriterPool writers;
        const double t_start = nn::now_ms();
        int64_t written = 0, skipped = 0;
        double model_ms = 0;

        for (int64_t i = 0; i < N; ++i) {
            const app::GeometryWarp& warp = warps[(size_t)group[(size_t)i]];
            const fs::path np = out_path(normal_dir, i,
                                         o.normal_format == "jpg" ? ".jpg" : ".png");
            const fs::path dp = out_path(depth_dir, i, ".png");
            const bool need_normal = o.want_normal && (o.overwrite || !fs::exists(np, ec));
            const bool need_depth = o.want_depth && (o.overwrite || !fs::exists(dp, ec));
            if (!need_normal && !need_depth) { ++skipped; continue; }

            const nn::Image img = nn::load_image(ds.image_filenames[(size_t)i]);
            if (img.empty()) {
                NN_LOG_WARN("skipping %s\n", ds.image_filenames[(size_t)i].c_str());
                continue;
            }
            const std::vector<float> src =
                resize_area(img, warp.sampleWidth(), warp.sampleHeight());

            const double t0 = nn::now_ms();
            std::vector<std::vector<float>> face_depth, face_normal;
            std::vector<float> face_rgb;
            for (int k = 0; k < warp.faces(); ++k) {
                warp.sampleFace(k, src.data(), face_rgb);
                metric3d::PredictOptions po;
                po.want_depth = need_depth;
                po.want_normal = need_normal;
                metric3d::Prediction p =
                    pred.predict(face_rgb.data(), warp.faceWidth(), warp.faceHeight(), po);
                NN_CHECK(p.width == warp.faceWidth() && p.height == warp.faceHeight(),
                         "the network returned %dx%d for a %dx%d face", p.width,
                         p.height, warp.faceWidth(), warp.faceHeight());
                face_depth.push_back(std::move(p.depth));
                face_normal.push_back(std::move(p.normal));
            }
            model_ms += nn::now_ms() - t0;

            const bool ray = o.ray_depth == Tri::Auto ? warp.defaultRayDepth()
                                                      : o.ray_depth == Tri::Yes;
            std::vector<float> depth, normal;
            warp.gather(face_depth, face_normal, ray, need_depth ? &depth : nullptr,
                        need_normal ? &normal : nullptr);

            if (need_normal) {
                app::WriteJob job;
                job.path = np.string();
                // A quality INSIDE 0..100 is what selects JPEG; the extension
                // is not consulted (nn/io/Image.h). 100 wrote a JPEG called
                // .png, which is lossy right where the mask edge is.
                job.quality = o.normal_format == "jpg" ? o.jpeg_quality : -1;
                job.image.width = warp.outWidth();
                job.image.height = warp.outHeight();
                job.image.channels = 3;
                job.image.data.resize(normal.size());
                // BLACK is the trainer's "no normal here": it decodes
                // byte/127.5 - 1 and masks on sum <= -2.366
                // (shaders/per_pixel_losses.slang), which mid-grey passes as a
                // valid normal of zero. A unit normal sums to at least -1.732.
                for (size_t p = 0; p * 3 < normal.size(); ++p) {
                    const float* v = &normal[p * 3];
                    uint8_t* out = &job.image.data[p * 3];
                    if (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] < 0.25f) {
                        out[0] = out[1] = out[2] = 0;
                        continue;
                    }
                    for (int c = 0; c < 3; ++c)
                        out[c] = (uint8_t)std::lround(
                            std::fmin(255.0f, std::fmax(0.0f, 127.5f + 127.5f * v[c])));
                }
                writers.submit(std::move(job));
            }
            if (need_depth) {
                // 0 is the trainer's "no ground truth here", so a covered pixel
                // never rounds into it.
                const float inv = o.depth_mm ? (float)warp.faceFocal()
                                             : 1.0f / depth_scale(depth);
                app::WriteJob job;
                job.path = dp.string();
                job.depth_w = warp.outWidth();
                job.depth_h = warp.outHeight();
                job.depth.resize(depth.size());
                for (size_t j = 0; j < depth.size(); ++j) {
                    if (depth[j] <= 0.0f) { job.depth[j] = 0; continue; }
                    const float v = depth[j] * inv * (o.depth_mm ? 1.0f : 65535.0f);
                    job.depth[j] = (uint16_t)std::lround(
                        std::fmin(65535.0f, std::fmax(1.0f, v)));
                }
                writers.submit(std::move(job));
            }

            ++written;
            if (written % 10 == 0 || i + 1 == N) {
                const double each = model_ms / (double)written;
                std::printf("\r%s   ",
                            format(G::log_progress,
                                   {(long long)(written + skipped), (long long)N,
                                    (long long)std::lround(each),
                                    human_time(each * (double)(N - i - 1))})
                                .c_str());
                std::fflush(stdout);
            }
        }
        writers.finish();
        std::printf("\r%s\n",
                    format(G::log_done, {(long long)written, (long long)skipped,
                                         human_time(nn::now_ms() - t_start)})
                        .c_str());
        return writers.failures() ? 1 : 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "\nerror: %s\n", e.what());
        return 1;
    }
}
