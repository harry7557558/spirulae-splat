// ssplat-sfm: the SfM pipeline CLI. Subcommands are the stage graph
// (src/sfm/README.md), each reading and writing files on disk so any one of
// them can be replaced by COLMAP's equivalent to bisect a failure:
//
//   ssplat-sfm auto    <image_dir> -o <workspace>      all of the below
//   ssplat-sfm extract <image|dir> -o <features>
//   ssplat-sfm match   <features>  -o <matches.bin>
//   ssplat-sfm map     <matches.bin> <features> -o <sparse/>
//   ssplat-sfm merge   <sparse/>   -o <merged/>
//   ssplat-sfm ba      <bal_problem.txt>               solver benchmark
//
// `--help` on any subcommand prints its flags. The self-checks are separate
// binaries (src/sfm/tests/, one per area).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <functional>
#include <set>
#include <string>
#include <vector>

#include <random>
#include <atomic>
#include <thread>

#include "sfm/core/CameraSetup.h"
#include "sfm/core/Features.h"
#include "sfm/core/Image.h"
#include "sfm/core/ImageLoader.h"
#include "sfm/core/Mask.h"
#include "sfm/core/Matches.h"
#include "sfm/feature/Matcher.h"
#include "sfm/feature/PairSelection.h"
#include "sfm/feature/Pairing.h"
#include "sfm/feature/Sift.h"
#include "sfm/feature/Verification.h"
#include "sfm/geometry/TwoView.h"
#include "sfm/map/Manager.h"
#include "sfm/map/Mapper.h"
#include "sfm/map/Merge.h"

// `ssplat-sfm ba`, in sfm_ba.cpp.
int cmdBa(int argc, char** argv);

namespace fs = std::filesystem;
using namespace sfm;

static bool isImageExt(const std::string& e) {
    std::string s;
    for (char c : e) s += (char)std::tolower((unsigned char)c);
    return s == ".jpg" || s == ".jpeg" || s == ".png" || s == ".bmp" || s == ".tga" ||
           s == ".ppm" || s == ".pgm";
}

// Path of `p` relative to the directory it was enumerated from.
//
// NOT fs::relative(): that canonicalizes both operands, which resolves
// symlinks. A directory of symlinked images (the natural way to stage a
// dataset without copying it) then relativizes to the *link targets*, yielding
// "../../.." paths that escape the output directory -- `extract` would write
// features next to the source images instead of into -o. Everything we pass
// here came out of a directory_iterator over `root`, so it is literally
// prefixed by `root` and a lexical comparison is both correct and cheap.
static fs::path relativeTo(const fs::path& p, const fs::path& root) {
    fs::path r = root;
    // "dir/" iterates to "dir/x.jpg" but has a trailing empty element, which
    // lexically_relative would mismatch into "../x.jpg". Drop it.
    if (!r.empty() && r.filename().empty()) r = r.parent_path();
    return p.lexically_relative(r);
}

// Wall-clock seconds since an arbitrary epoch, for the stage timings `auto`
// reports.
static double now() {
    using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count();
}

// Mean/median reprojection error over a reconstruction's observations. The one
// number that says whether a model is sane; reported by both `map` and `auto`.
static void reprojStats(const Reconstruction& rec, const std::vector<FeatureSet>& feats,
                        double& mean, double& median, size_t& nobs) {
    std::vector<double> e;
    for (const auto& kv : rec.points3D)
        for (const TrackElement& t : kv.second.track) {
            auto img = rec.images.find(t.image_id);
            if (img == rec.images.end() || !img->second.registered) continue;
            auto cam = rec.cameras.find(img->second.camera_id);
            if (cam == rec.cameras.end()) continue;
            Vec3 pc = mul(img->second.pose.R, kv.second.xyz) + img->second.pose.t;
            if (pc.z <= 0) continue;
            Vec2 px = cam->second.project(pc);
            const Keypoint& k = feats[t.image_id].keypoints[t.point2D_idx];
            e.push_back(std::hypot(px.x - k.x, px.y - k.y));
        }
    nobs = e.size();
    mean = median = 0;
    if (e.empty()) return;
    double s = 0;
    for (double v : e) s += v;
    mean = s / e.size();
    std::nth_element(e.begin(), e.begin() + e.size() / 2, e.end());
    median = e[e.size() / 2];
}

// Sample an RGB color at every keypoint from the (already decoded) source image
// and store it in the feature set, so the point cloud can be colored without a
// second decode pass. No-op if the image was decoded without color.
static void sampleFeatureColors(FeatureSet& fs, const GrayImage& img) {
    if (!img.hasColor()) return;
    fs.colors.resize((size_t)fs.count() * 3);
    for (uint32_t i = 0; i < fs.count(); i++)
        sampleColor(img, fs.keypoints[i].x, fs.keypoints[i].y, &fs.colors[(size_t)i * 3]);
}

// Write every reconstruction the mapper produced as <dir>/0, <dir>/1, ...
// (D41). This is COLMAP's layout for a dataset that does not form one connected
// view graph -- `sparse/0` is the model with the most 3D points, and the rest
// are the components a later merge step can align onto it. A single-model
// dataset therefore still writes exactly `sparse/0`, as before.
// Put a freshly extracted set into the source image's frame and attach what
// EXIF said about the camera. Called once per image, after anything that
// indexes the *decoded* image (colors, masks) is done.
static void finishFeatures(FeatureSet& fs, const GrayImage& img) {
    scaleKeypoints(fs, img.orig_width, img.orig_height);
    fs.exif_focal = exifFocalPx(img.exif, fs.width, fs.height);
    fs.exif_camera = exifCameraKey(img.exif, fs.width, fs.height);
}

static void writeModels(const std::vector<Reconstruction>& models, const fs::path& dir,
                        bool verbose) {
    for (size_t i = 0; i < models.size(); i++) {
        fs::path p = dir / std::to_string(i);
        fs::create_directories(p);
        models[i].writeBinary(p.string());
        if (verbose)
            fprintf(stderr, "[map] wrote model %zu (%u images, %zu points) to %s\n", i,
                    models[i].numRegistered(), models[i].points3D.size(), p.string().c_str());
    }
    // A re-run that produces fewer models than the last one left numbered
    // directories behind, and `--resume` would read them back as if they were
    // part of this reconstruction. Remove them.
    if (!fs::is_directory(dir)) return;
    for (const auto& e : fs::directory_iterator(dir)) {
        if (!e.is_directory()) continue;
        const std::string name = e.path().filename().string();
        if (name.find_first_not_of("0123456789") != std::string::npos) continue;
        size_t idx = std::stoull(name);
        if (idx < models.size()) continue;
        if (!fs::exists(e.path() / "images.bin")) continue;
        std::error_code ec;
        fs::remove_all(e.path(), ec);
        if (verbose && !ec)
            fprintf(stderr, "[map] removed stale model directory %s from an earlier run\n",
                    e.path().string().c_str());
    }
}

// Distinct images covered by any model. Not the sum of the model sizes: models
// deliberately overlap by up to max_model_overlap images (D41), which is what a
// merge step aligns on, so summing double-counts the joins.
static size_t distinctRegistered(const std::vector<Reconstruction>& models) {
    std::set<uint32_t> ids;
    for (const Reconstruction& m : models)
        for (const auto& kv : m.images)
            if (kv.second.registered) ids.insert(kv.first);
    return ids.size();
}

// One line per model beyond the first, so a fragmented capture is visible in
// the summary rather than only in the directory listing.
static void printExtraModels(const std::vector<Reconstruction>& models,
                             const std::vector<FeatureSet>& feats, const char* indent) {
    for (size_t i = 1; i < models.size(); i++) {
        double mn = 0, md = 0;
        size_t nobs = 0;
        reprojStats(models[i], feats, mn, md, nobs);
        printf("%smodel %zu: %u images, %zu points, %.3f px mean\n", indent, i,
               models[i].numRegistered(), models[i].points3D.size(), mn);
    }
}

// Mean/median reprojection error straight from a model, with no FeatureSets in
// hand: Image::points2D already carries the keypoint coordinates. `merge`
// starts from models on disk, where that is all there is. Observations behind
// the camera are excluded rather than counted as infinite, as reprojStats does.
static void modelReprojStats(const Reconstruction& rec, double& mean, double& median,
                             size_t& nobs) {
    std::vector<double> e;
    for (const auto& kv : rec.points3D)
        for (const TrackElement& t : kv.second.track) {
            auto img = rec.images.find(t.image_id);
            if (img == rec.images.end() || !img->second.registered) continue;
            auto cam = rec.cameras.find(img->second.camera_id);
            if (cam == rec.cameras.end() || t.point2D_idx >= img->second.points2D.size()) continue;
            double r = reprojErrorAt(cam->second, img->second.pose,
                                     img->second.points2D[t.point2D_idx], kv.second.xyz);
            if (r < 1e29) e.push_back(r);
        }
    nobs = e.size();
    mean = median = 0;
    if (e.empty()) return;
    double s = 0;
    for (double v : e) s += v;
    mean = s / e.size();
    std::nth_element(e.begin(), e.begin() + e.size() / 2, e.end());
    median = e[e.size() / 2];
}

// ---- merging (D43) ----
// Everything `merge` and `auto` share: run the session to convergence, then
// bundle-adjust the models that absorbed something. A merged model is two
// independently optimized halves glued along a seam that has never been
// optimized as one, so the BA is not cosmetic -- it is where the shared images'
// two answers become one. Models that merged nothing are moved through
// untouched, so a dataset with nothing to merge is bit-identical to before.
struct MergeSummary {
    size_t before = 0, after = 0, merges = 0, refused = 0;
    double seconds = 0, ba_seconds = 0;
};

static std::vector<Reconstruction> mergeModels(std::vector<Reconstruction> models,
                                               const MergeOptions& mo, bool refine, int device,
                                               MergeSummary& sum) {
    sum.before = sum.after = models.size();
    if (models.size() < 2) return models;
    const double t0 = now();
    MergeSession session(std::move(models), mo);
    sum.merges = session.runAuto();
    for (const MergeAttempt& a : session.log())
        if (!a.merged) sum.refused++;
    if (refine && sum.merges) {
        const double tb = now();
        VkContext ctx;  // one device + pipeline set for every model, as the mapper does
        for (size_t i : session.changed()) {
            Reconstruction& m = session.modelMut(i);
            BundleOptions bo;
            bo.device = device;
            bo.verbose = false;
            bo.shared_ctx = &ctx;
            double cost = runGlobalBA(m, bo);
            size_t robs = 0, rpts = 0;
            filterModel(m, mo.filter_reproj_error, mo.min_tri_angle_deg, robs, rpts);
            if (mo.verbose)
                fprintf(stderr,
                        "[merge] model %zu re-bundled (cost %.3e), filtered %zu obs / %zu points, "
                        "%zu remain\n",
                        i, cost, robs, rpts, m.points3D.size());
        }
        sum.ba_seconds = now() - tb;
    }
    std::vector<Reconstruction> out = session.take();
    sum.after = out.size();
    sum.seconds = now() - t0;
    return out;
}

// ---- reading models back off disk ----
// An input is either one model directory or a `sparse/` holding 0, 1, 2, ...
// Numbered sub-directories are read in numeric order, so index 0 stays the
// model with the most structure and therefore the natural anchor.
static bool isModelDir(const fs::path& p) {
    return fs::exists(p / "cameras.bin") && fs::exists(p / "images.bin");
}

static bool collectModelDirs(const std::string& input, std::vector<fs::path>& out) {
    fs::path p(input);
    if (!fs::is_directory(p)) {
        fprintf(stderr, "%s is not a directory\n", input.c_str());
        return false;
    }
    if (isModelDir(p)) {
        out.push_back(p);
        return true;
    }
    std::vector<std::pair<long, fs::path>> numbered;
    for (const auto& e : fs::directory_iterator(p)) {
        if (!e.is_directory() || !isModelDir(e.path())) continue;
        const std::string n = e.path().filename().string();
        if (n.empty() || n.find_first_not_of("0123456789") != std::string::npos) continue;
        numbered.emplace_back(std::stol(n), e.path());
    }
    if (numbered.empty()) {
        fprintf(stderr, "%s holds no model (no cameras.bin, no numbered sub-model)\n",
                input.c_str());
        return false;
    }
    std::sort(numbered.begin(), numbered.end());
    for (auto& kv : numbered) out.push_back(kv.second);
    return true;
}

static bool readModels(const std::string& dir, std::vector<Reconstruction>& models, bool verbose) {
    std::vector<fs::path> dirs;
    if (!collectModelDirs(dir, dirs)) return false;
    for (const fs::path& d : dirs) {
        try {
            models.push_back(Reconstruction::readBinary(d.string()));
        } catch (const std::exception& e) {
            fprintf(stderr, "cannot read %s: %s\n", d.string().c_str(), e.what());
            return false;
        }
        if (verbose)
            printf("   read %s: %u images, %zu points\n", d.string().c_str(),
                   models.back().numRegistered(), models.back().points3D.size());
    }
    return true;
}

// Whether the final principal-point pass runs unless the CLI says otherwise.
// Both switches exist either way (`--final-principal-point` /
// `--no-final-principal-point`) so the default can move without the flags
// changing meaning.
static constexpr bool kFinalPrincipalPointDefault = true;

// One last global bundle adjustment per model with the principal point
// released (D51). COLMAP's documentation: hold it during reconstruction, where
// it is ill-posed, then "try to refine the principal point in global bundle
// adjustment" once every image is in and the problem is constrained --
// "especially when sharing intrinsic parameters between multiple images",
// which is why a camera group needs `pp_min_images` behind it to qualify.
//
// Runs after the manager, on models nothing else will touch, so a group whose
// principal point turns out to be badly conditioned can only spoil its own
// intrinsics -- there is no growth pass left to build on the result. Mapper::
// polish skips models with more than one camera group; measured, the gain
// tracks how far the *reference's* own principal point sits from the image
// centre (D51).
static std::vector<Reconstruction> polishModels(Mapper& mapper,
                                                std::vector<Reconstruction> models,
                                                bool verbose, double& secs) {
    const double t0 = now();
    for (Reconstruction& m : models)
        if (m.numRegistered() >= 2) m = mapper.polish(m);
    secs = now() - t0;
    if (verbose)
        fprintf(stderr, "[map] final principal-point refinement over %zu model(s), %.1f s\n",
                models.size(), secs);
    return models;
}

// The post-mapping loop (D44): merge, grow, prune, reseed until a round changes
// nothing. Shared by `auto` and `map`; needs the mapper because growing and
// refining are the mapper's own operations over the same feature database.
static std::vector<Reconstruction> manageModels(Mapper& mapper,
                                                std::vector<Reconstruction> models,
                                                const ManagerOptions& mgopt, ManagerStats& st,
                                                double& secs) {
    const double t0 = now();
    ModelManager mgr(mapper, mgopt);
    std::vector<Reconstruction> out = mgr.run(std::move(models));
    st = mgr.stats();
    secs = now() - t0;
    return out;
}

// Feature stems carry no extension; put the real filename (relative to the
// image dir, so sub-folders survive) back into every model for COLMAP tooling.
// The directory is walked once, not once per model -- a fragmented capture can
// produce dozens (D41).
static void resolveImageNames(std::vector<Reconstruction>& models, const std::string& imagedir) {
    if (imagedir.empty()) return;
    std::map<std::string, std::string> stem2name;
    for (const auto& e : fs::recursive_directory_iterator(imagedir))
        if (e.is_regular_file()) {
            fs::path rel = relativeTo(e.path(), imagedir);
            fs::path stem = rel;
            stem.replace_extension();
            stem2name[stem.generic_string()] = rel.generic_string();
        }
    for (Reconstruction& rec : models)
        for (auto& kv : rec.images) {
            auto it = stem2name.find(kv.second.name);
            if (it != stem2name.end()) kv.second.name = it->second;
        }
}

// How images are grouped into cameras (D17, D40). COLMAP's ImageReader keys its
// default on the EXIF make/model, which we cannot read yet (D5).
//
// Every mode splits on image resolution first and is only then free to group:
// two images of different sizes cannot share a principal point, so a camera
// spanning both fits neither (D40). That collapses what used to be two modes --
// `single` ("one camera for everything") and `auto` ("one per resolution") --
// into one, since honouring the resolution split is exactly what `auto` did;
// `auto` remains accepted as a spelling of `single`.
// Camera grouping (--camera-mode), the distortion model and the starting focal
// all live in sfm/core/CameraSetup.h, which the CLI drives; --camera-model and
// --focal accept either a dataset-wide value or PREFIX=VALUE for one group
// (D46). The accepted model names come from sfm/core/Camera.h's kCamModelInfo
// table, the single source of truth.
static bool parseCamModel(const std::string& v, CamModel& out) {
    return parseCamModelName(v, out);
}

// Report the grouping decision, one line per camera. Worth printing in full:
// a wrong --camera-mode is otherwise invisible until the intrinsics come out
// strange, and a mixed-model capture is exactly where it goes wrong.
static void printCameraSetup(const char* tag, const CameraSetup& cs,
                             const CameraSetupOptions& sopt, size_t nimages) {
    std::map<uint32_t, size_t> counts;
    for (uint32_t id : cs.ids) counts[id]++;
    fprintf(stderr, "[%s] camera mode %s -> %u camera(s) over %zu images\n", tag,
            cameraModeName(cs.mode_used), cs.count(), nimages);
    if (cs.mode_switched)
        fprintf(stderr, "[%s]   (%zu distinct frame sizes over %zu images: this is a photo "
                "collection, not one camera's capture -- --camera-mode folder overrides)\n",
                tag, cs.dim_buckets, nimages);
    if (cs.exif_focal_images)
        fprintf(stderr, "[%s] EXIF: %zu/%zu images carry a focal length%s\n", tag,
                cs.exif_focal_images, nimages,
                sopt.exif_focal ? "" : " (ignored, --no-exif-focal)");
    // Capped: --camera-mode image on an internet collection makes one camera
    // per image, and a thousand lines of stderr helps nobody.
    const size_t kMaxLines = 8;
    size_t shown = 0;
    for (const auto& kv : cs.cameras) {
        if (shown++ >= kMaxLines) {
            fprintf(stderr, "[%s]   ... and %zu more camera(s)\n", tag,
                    cs.cameras.size() - kMaxLines);
            break;
        }
        const Camera& c = kv.second;
        fprintf(stderr, "[%s]   camera %u: %zu image(s), %dx%d, %s, f=%.1f%s\n", tag, kv.first,
                counts[kv.first], c.width, c.height, camInfo(c.model).cli_name, c.focal(),
                cs.focal_known.count(kv.first)   ? " (prior)"
                : cs.focal_given.count(kv.first) ? " (given)"
                                                 : " (guess)");
    }
}

// -----------------------------------------------------------------------
// extract
// -----------------------------------------------------------------------

// Batch-extract every image in `imagedir` to `outdir`/<stem>.bin. Shared by
// `ssplat-sfm extract DIR` and `ssplat-sfm auto`.
struct ExtractStats {
    size_t images = 0, failed = 0, unreadable = 0;
    uint64_t features = 0;
    // Masking (D39), all zero when no mask directory was given.
    size_t masked_images = 0;     // images that found a mask file
    size_t unmasked_images = 0;   // images that did not
    size_t mask_unreadable = 0;   // found a mask file but could not decode it
    uint64_t masked_out = 0;      // keypoints dropped by masks
    std::string first_unmasked;   // an example, for the warning
    bool warned_empty = false;    // "this mask masked out everything", warned once
};

// Warn once per distinct (mask size, image size) pair when the two disagree in
// *aspect*. A size difference is expected and handled -- masks are sampled in
// uv (D39) -- but a different aspect ratio means the mask was made for a
// differently-cropped image, and uv sampling will happily stretch it over the
// wrong content. That is a dataset bug we can see and the user cannot.
static void checkMaskShape(const std::string& mask_path, const Mask& m,
                           const std::pair<int, int>& img_dims) {
    static std::set<std::array<int, 4>> warned;
    if (m.empty() || img_dims.first <= 0 || img_dims.second <= 0) return;
    const double am = (double)m.width / m.height;
    const double ai = (double)img_dims.first / img_dims.second;
    if (std::fabs(am - ai) <= 0.01 * ai) return;
    if (!warned.insert({m.width, m.height, img_dims.first, img_dims.second}).second) return;
    fprintf(stderr,
            "[extract] WARNING: mask %s is %dx%d (aspect %.3f) but the image is %dx%d "
            "(aspect %.3f); it will be stretched to fit. Suppressing further warnings for "
            "this size pair.\n",
            mask_path.c_str(), m.width, m.height, am, img_dims.first, img_dims.second, ai);
}

// A mask that drops nearly every keypoint in the dataset is far more often
// inverted than intended. Both conventions exist in the wild -- COLMAP's (and
// ours) is white = keep, but some tools export the *excluded* region as white,
// and one of the benchmark captures does exactly that: its masks are 6% white,
// covering only the ground below the horizon, so honouring them literally
// reconstructs the floor under a moving camera and nothing else.
//
// This is a warning and not an automatic flip, because a legitimately small
// mask exists too: an object-centric capture masks away everything but the
// object. The user is the one who knows which they have.
static void warnIfMasksLookInverted(const ExtractStats& st) {
    const uint64_t before = st.features + st.masked_out;
    if (!st.masked_images || before == 0) return;
    const double dropped = (double)st.masked_out / (double)before;
    if (dropped < 0.7) return;
    fprintf(stderr,
            "[extract] WARNING: masks dropped %.0f%% of all keypoints. If this capture is not "
            "object-centric, the masks are probably inverted -- this pipeline keeps what is "
            "*white*, as COLMAP does. Re-run with --no-masks to check.\n",
            100.0 * dropped);
}

static int extractDirectory(const std::string& imagedir, const fs::path& outdir,
                            const SiftOptions& opt, int max_image_size, int decode_threads,
                            int decode_budget_mb, const std::string& maskdir,
                            ExtractStats& stats) {
    // Recursive: datasets that carry per-folder intrinsics (ppisp) keep images
    // in images/<camera>/, and the folder is the grouping key (D17). Feature
    // files mirror the tree so stems from different folders cannot collide.
    //
    // A mask directory nested inside the image directory is skipped: masks are
    // themselves PNGs, and extracting features from them would silently double
    // the image count with garbage views (`ssplat-sfm auto DATASET` with the images
    // one level up is exactly the layout that trips this).
    std::error_code skip_ec;
    const bool skip_masks = !maskdir.empty() && fs::is_directory(maskdir, skip_ec);
    std::vector<fs::path> found;
    for (auto it = fs::recursive_directory_iterator(imagedir);
         it != fs::recursive_directory_iterator(); ++it) {
        if (skip_masks && it->is_directory() && fs::equivalent(it->path(), maskdir, skip_ec)) {
            it.disable_recursion_pending();
            continue;
        }
        if (it->is_regular_file() && isImageExt(it->path().extension().string()))
            found.push_back(it->path());
    }
    if (found.empty()) { fprintf(stderr, "no images in %s\n", imagedir.c_str()); return 1; }
    std::sort(found.begin(), found.end());

    // Probe every header once (the old comparator re-probed O(n log n) times),
    // both to order the batch and to size the decoder's memory.
    std::vector<fs::path> imgs;
    std::vector<std::pair<int, int>> dims;
    for (const fs::path& p : found) {
        int w = 0, h = 0;
        if (!imageSize(p.string(), w, h) || w <= 0 || h <= 0) {
            fprintf(stderr, "[extract] skipping %s: not a decodable image\n",
                    p.filename().string().c_str());
            stats.unreadable++;
            continue;
        }
        imgs.push_back(p);
        dims.emplace_back(w, h);
    }
    if (imgs.empty()) { fprintf(stderr, "no decodable images in %s\n", imagedir.c_str()); return 1; }

    // Largest-first so SiftExtractor allocates device buffers exactly once.
    std::vector<size_t> order(imgs.size());
    for (size_t i = 0; i < order.size(); i++) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return (int64_t)dims[a].first * dims[a].second > (int64_t)dims[b].first * dims[b].second;
    });
    std::vector<std::string> paths(order.size());
    std::vector<std::pair<int, int>> sorted_dims(order.size());
    for (size_t k = 0; k < order.size(); k++) {
        paths[k] = imgs[order[k]].string();
        sorted_dims[k] = dims[order[k]];
    }

    fs::create_directories(outdir);

    ImageLoadOptions lopt;
    lopt.max_image_size = max_image_size;
    lopt.num_threads = decode_threads;
    lopt.want_color = true;  // sample per-keypoint colors while the image is hot
    if (decode_budget_mb > 0) lopt.memory_budget_bytes = (size_t)decode_budget_mb << 20;

    // Resolve each image's mask (D39) up front, in the same order as `paths`,
    // so the decode pool can pick them up. Doing it here rather than inside the
    // consumer also means a mask directory that matches nothing is reported
    // before the GPU stage burns an hour on unmasked features.
    MaskIndex masks(maskdir);
    if (!maskdir.empty() && !masks.valid())
        fprintf(stderr, "[extract] WARNING: mask directory %s does not exist; "
                        "extracting without masks\n", maskdir.c_str());
    if (masks.valid()) {
        lopt.mask_paths.assign(paths.size(), std::string());
        for (size_t k = 0; k < paths.size(); k++) {
            std::string rel = relativeTo(paths[k], imagedir).generic_string();
            lopt.mask_paths[k] = masks.find(rel);
            if (lopt.mask_paths[k].empty()) {
                stats.unmasked_images++;
                if (stats.first_unmasked.empty()) stats.first_unmasked = rel;
            } else {
                stats.masked_images++;
            }
        }
        printf("   masks: %zu/%zu images matched in %s\n", stats.masked_images, paths.size(),
               maskdir.c_str());
        if (stats.masked_images == 0) {
            fprintf(stderr,
                    "[extract] ERROR: no mask in %s matches any image (tried e.g. \"%s.png\" "
                    "and \"<stem>.png\"). Fix --masks or drop it; continuing unmasked would "
                    "quietly produce a different reconstruction.\n",
                    maskdir.c_str(), stats.first_unmasked.c_str());
            return 1;
        }
        if (stats.unmasked_images)
            fprintf(stderr,
                    "[extract] WARNING: %zu image(s) have no mask (e.g. %s); their features "
                    "are kept in full\n",
                    stats.unmasked_images, stats.first_unmasked.c_str());
    }
    ImageLoadPlan plan = planImageLoad(sorted_dims, lopt);
    if (opt.verbose)
        fprintf(stderr,
                "[extract] %zu images, decoding on %d thread(s), window %d "
                "(~%zu MB peak in-flight)\n",
                paths.size(), plan.num_threads, plan.window,
                ((size_t)plan.num_threads * plan.decode_peak_bytes +
                 (size_t)plan.window * plan.held_bytes) >> 20);

    SiftExtractor ext(opt);
    loadImagesInOrder(
        paths, plan, lopt,
        [&](size_t k, GrayImage& img) {
            FeatureSet f = ext.extract(img);
            sampleFeatureColors(f, img);
            uint32_t dropped = 0;
            if (!lopt.mask_paths.empty() && !lopt.mask_paths[k].empty()) {
                if (img.mask.empty()) {
                    stats.mask_unreadable++;
                    fprintf(stderr, "[extract] WARNING: cannot decode mask %s; %s kept unmasked\n",
                            lopt.mask_paths[k].c_str(),
                            fs::path(paths[k]).filename().string().c_str());
                } else {
                    checkMaskShape(lopt.mask_paths[k], img.mask, sorted_dims[k]);
                    const uint32_t before = f.count();
                    dropped = applyMask(f, img.mask);
                    stats.masked_out += dropped;
                    // Masking away an entire image is what an inverted mask
                    // (or one whose keep value is 0) looks like, and it is
                    // otherwise invisible until the reconstruction is short of
                    // images. Warn once; the run still continues.
                    if (before && dropped == before && !stats.warned_empty) {
                        stats.warned_empty = true;
                        fprintf(stderr,
                                "[extract] WARNING: mask %s left no keypoints at all in %s. "
                                "Masks keep nonzero pixels and ignore zero ones -- an inverted "
                                "mask masks out the whole image.\n",
                                lopt.mask_paths[k].c_str(),
                                fs::path(paths[k]).filename().string().c_str());
                    }
                }
            }
            // Back to the source file's coordinates, so cameras.bin describes
            // the images the user has rather than the extractor's working copy
            // (D46). Everything that reads a keypoint against the decoded image
            // -- the color sampling above, the mask (sampled in uv) -- has
            // already run.
            finishFeatures(f, img);
            fs::path rel = relativeTo(paths[k], imagedir);
            fs::path out = outdir / rel;
            out.replace_extension(".bin");
            fs::create_directories(out.parent_path());
            writeFeatures(out.string(), f);
            stats.features += f.count();
            stats.images++;
            if (opt.verbose) {
                printf("[%zu/%zu] %s: %u features", stats.images, paths.size(),
                       fs::path(paths[k]).filename().string().c_str(), f.count());
                if (dropped) printf(" (%u masked out)", dropped);
                printf(" -> %s\n", out.string().c_str());
            }
        },
        [&](size_t k, const std::string& err) {
            fprintf(stderr, "[extract] FAILED %s: %s\n",
                    fs::path(paths[k]).filename().string().c_str(), err.c_str());
            stats.failed++;
        });
    return 0;
}

static int cmdExtract(int argc, char** argv) {
    std::string image, output, maskdir;
    SiftOptions opt;
    int max_image_size = 3200;
    int decode_threads = 0;    // batch decode pool; 0 = hardware_concurrency
    int decode_budget_mb = 0;  // 0 = ImageLoadOptions default
    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--output" || a == "-o") output = next();
        else if (a == "--masks" || a == "--mask-dir") maskdir = next();
        else if (a == "--max-features") opt.max_num_features = std::stoi(next());
        else if (a == "--octaves") opt.num_octaves = std::stoi(next());
        else if (a == "--peak-threshold") opt.peak_threshold = std::stod(next());
        else if (a == "--edge-threshold") opt.edge_threshold = std::stod(next());
        else if (a == "--max-orientations") opt.max_num_orientations = std::stoi(next());
        else if (a == "--max-image-size") max_image_size = std::stoi(next());
        else if (a == "--decode-threads") decode_threads = std::stoi(next());
        else if (a == "--decode-budget") decode_budget_mb = std::stoi(next());
        else if (a == "--device") opt.device = std::stoi(next());
        else if (a == "--profile") opt.profile = true;
        else if (a == "--quiet") opt.verbose = false;
        else if (a == "--spv-path") opt.spv_path = next();
        else if (a[0] != '-') image = a;
        else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    if (image.empty()) {
        fprintf(stderr,
                "usage: ssplat-sfm extract <image|dir> [--output FILE|DIR]\n"
                "  [--masks DIR] [--max-features N] [--octaves N] [--peak-threshold X]\n"
                "  [--edge-threshold X] [--max-orientations N] [--max-image-size N]\n"
                "  [--decode-threads N] [--decode-budget MB]\n"
                "  [--device I] [--profile] [--quiet] [--spv-path FILE]\n"
                "A directory extracts every image in it to <output>/<name>.bin (default dir\n"
                "'features'), reusing one GPU context and processing largest-first.\n"
                "Decoding runs on a pool sized to fit --decode-budget (default 1024 MB);\n"
                "--decode-threads 1 decodes inline. Peak memory does not grow with the\n"
                "number of images.\n"
                "--masks DIR drops keypoints on zero mask pixels; masks are matched by\n"
                "name (<image>.png, <stem>.png, <stem>_mask.png, ...) and sampled in uv,\n"
                "so a mask needs no particular resolution.\n"
                "Keypoints are written in the *source* image's coordinates even when\n"
                "--max-image-size downscaled it, and each file's EXIF focal length and\n"
                "camera identity are recorded alongside them, so `match` and `map` can\n"
                "build intrinsics for the images on disk without seeing them (D46).\n");
        return 1;
    }

    // ---- directory (batch) ----
    if (fs::is_directory(image)) {
        fs::path outdir = output.empty() ? fs::path("features") : fs::path(output);
        ExtractStats st;
        int rc = extractDirectory(image, outdir, opt, max_image_size, decode_threads,
                                  decode_budget_mb, maskdir, st);
        if (rc) return rc;
        printf("done: %zu images, %llu features total", st.images,
               (unsigned long long)st.features);
        if (st.masked_out)
            printf(" (%llu masked out over %zu masked images)",
                   (unsigned long long)st.masked_out, st.masked_images);
        if (st.failed || st.unreadable)
            printf(" (%zu failed to decode, %zu unreadable headers)", st.failed, st.unreadable);
        printf("\n");
        warnIfMasksLookInverted(st);
        return (st.failed || st.unreadable) ? 1 : 0;
    }

    // ---- single image ----
    // One image has no relative path to key on, so the mask is looked up by
    // filename -- MaskIndex's basename fallback resolves it either way.
    std::string maskpath;
    if (!maskdir.empty()) {
        maskpath = MaskIndex(maskdir).find(fs::path(image).filename().generic_string());
        if (maskpath.empty())
            fprintf(stderr, "[sfm] WARNING: no mask for %s in %s\n", image.c_str(),
                    maskdir.c_str());
    }
    GrayImage img = loadGrayImage(image, max_image_size, /*want_color=*/true, maskpath);
    if (opt.verbose) fprintf(stderr, "[sfm] %s -> %dx%d gray\n", image.c_str(), img.width, img.height);
    SiftExtractor ext(opt);
    FeatureSet fs = ext.extract(img);
    sampleFeatureColors(fs, img);
    uint32_t masked_out = applyMask(fs, img.mask);
    finishFeatures(fs, img);
    printf("features: %u  dim: %u  image: %dx%d", fs.count(), fs.dim, fs.width, fs.height);
    if (fs.exif_focal > 0) printf("  exif focal: %.1f px", fs.exif_focal);
    if (!img.mask.empty())
        printf("  mask: %dx%d, %u masked out", img.mask.width, img.mask.height, masked_out);
    printf("\n");
    if (!output.empty()) {
        writeFeatures(output, fs);
        if (opt.verbose) fprintf(stderr, "[sfm] wrote %s\n", output.c_str());
    }
    return 0;
}

// -----------------------------------------------------------------------
// match
// -----------------------------------------------------------------------

// Match + verify a feature directory into a MatchesDatabase. Shared by
// `ssplat-sfm match` and `ssplat-sfm auto`.
struct MatchStats {
    size_t images = 0, pairs = 0, kept = 0, scored = 0;
    uint64_t inliers = 0, putative = 0;
    double select_seconds = 0;
};

// Calibrated verification (D45). A fisheye pair verified on raw pixels asks a
// pinhole fundamental matrix to explain rays 100 deg off axis, which it cannot
// represent at all; the correspondences that carry the wide field of view are
// thrown away as outliers. With a camera model in hand we verify on unit
// bearings instead, where the epipolar constraint is exact at any FOV.
//
// The model needs a focal length before it can produce bearings, and the
// geometric default (diag/pi) is off by ~1.7x on a 200 deg lens, which is
// enough to warp the bearings and lose most of the benefit -- so unless one is
// given we search for it on a sample of pairs first (`bootstrapFocal`).
struct VerifyCalibration {
    CameraSetupOptions setup;
    size_t sample_pairs = 150;  // pairs per group used by the focal search
    // outputs
    CameraSetup cameras;        // what the grouping decided, for the caller to reuse
    bool used_bearings = false;
};

static int matchFeatureDir(const std::string& featdir, std::vector<FeatureSet>& feats,
                           MatchesDatabase& db, const MatchOptions& opt, PairMode mode,
                           int overlap, const PairSelectionOptions& popt, bool verify,
                           const TwoViewOptions& tvopt, int threads, bool verbose,
                           MatchStats& stats, VerifyCalibration* calib = nullptr) {
    // Load every features.bin under the directory (recursively -- the tree
    // mirrors the image tree), sorted by name for stable indices. The image
    // name is the relative path without ".bin", which is also COLMAP's
    // convention for `images.bin` names.
    std::vector<fs::path> files;
    for (const auto& e : fs::recursive_directory_iterator(featdir))
        if (e.is_regular_file() && e.path().extension() == ".bin") files.push_back(e.path());
    std::sort(files.begin(), files.end());
    if (files.size() < 2) {
        fprintf(stderr, "need >=2 feature files in %s\n", featdir.c_str());
        return 1;
    }

    feats.assign(files.size(), FeatureSet());
    db.images.resize(files.size());
    for (size_t i = 0; i < files.size(); i++) {
        feats[i] = readFeatures(files[i].string());
        fs::path rel = relativeTo(files[i], featdir);
        rel.replace_extension();
        db.images[i] = {rel.generic_string(), feats[i].count()};
    }
    stats.images = files.size();

    std::vector<std::pair<uint32_t, uint32_t>> pairs;
    if (mode == PairMode::Prefilter) {
        stats.scored = files.size() * (files.size() - 1) / 2;
        std::function<void(size_t, size_t)> sp;
        if (verbose)
            sp = [&](size_t done, size_t total) {
                fprintf(stderr, "\r[match] %zu/%zu pairs scored", done, total);
            };
        double t0 = now();
        pairs = prefilterPairs(feats, popt, sp);
        stats.select_seconds = now() - t0;
        if (verbose)
            fprintf(stderr,
                    "\n[match] pair selection kept %zu/%zu pairs "
                    "(top-%u features, %u neighbors, %.1f s)\n",
                    pairs.size(), stats.scored, popt.num_features, popt.num_neighbors,
                    stats.select_seconds);
    } else {
        pairs = generatePairs((uint32_t)files.size(), mode, overlap);
    }
    stats.pairs = pairs.size();
    if (verbose)
        fprintf(stderr, "[match] %zu images, %zu pairs (%s)\n", files.size(), pairs.size(),
                mode == PairMode::Exhaustive  ? "exhaustive"
                : mode == PairMode::Sequential ? "sequential"
                                               : "prefilter");
    if (verbose && mode == PairMode::Prefilter)
        fprintf(stderr, "[match] pair selection: top-%u features, %u neighbors\n",
                popt.num_features, popt.num_neighbors);

    BruteForceMatcher matcher(opt);
    auto matchFn = [&](size_t b, size_t e, std::vector<std::vector<FeatureMatch>>& mout) {
        matcher.matchBatch(feats, pairs, b, e, mout);
    };
    std::function<void(size_t, size_t)> progress;
    if (verbose)
        progress = [&](size_t done, size_t total_pairs) {
            if (done % 200 == 0 || done == total_pairs)
                fprintf(stderr, "\r[match] %zu/%zu pairs matched", done, total_pairs);
        };

    if (verify) {
        // Verification runs on a worker pool fed by the (serial, GPU-bound)
        // matcher -- it is the pipeline's dominant cost, see sfm/feature/Verification.h.
        VerificationOptions vopt;
        vopt.two_view = tvopt;
        vopt.num_threads = threads;
        vopt.match_batch_pairs = opt.batch_pairs;

        // Calibrated verification (D45), for fisheye captures only: see
        // VerifyCalibration. Everything else keeps the pixel path, where the
        // pinhole assumption is exact and results are long settled.
        BearingCache bc;
        if (calib) {
            CameraSetup& cs = calib->cameras;
            cs = buildCameras(db.images, feats, calib->setup);
            if (verbose) printCameraSetup("match", cs, calib->setup, feats.size());
            if (cs.anyWide()) {
                // Sample pairs evenly across the list -- taking a prefix would
                // sample one part of the capture, and pair lists are ordered.
                // The sample is per group, so it is scaled by the group count.
                const size_t want = calib->sample_pairs * std::max<size_t>(1, cs.count());
                std::vector<std::pair<uint32_t, uint32_t>> sample;
                size_t stride = std::max<size_t>(1, pairs.size() / std::max<size_t>(1, want));
                for (size_t p = 0; p < pairs.size() && sample.size() < want; p += stride)
                    sample.push_back(pairs[p]);
                std::vector<std::vector<FeatureMatch>> sm, chunk;
                for (size_t b = 0; b < sample.size(); b += 16) {
                    size_t e = std::min(b + 16, sample.size());
                    matcher.matchBatch(feats, sample, b, e, chunk);
                    for (size_t k = b; k < e; k++) sm.push_back(std::move(chunk[k - b]));
                }
                bootstrapGroupFocals(feats, cs.ids, sample, sm, cs.cameras, cs.focal_given, tvopt,
                                     calib->sample_pairs, threads, verbose);
                std::vector<Camera> percam(feats.size());
                for (size_t i = 0; i < feats.size(); i++) percam[i] = cs.cameras.at(cs.ids[i]);
                double t_b = now();
                // A mixed capture calibrates *everything*: its cross pairs have
                // a fisheye on one side, and the pixel path cannot represent
                // those at all. A wholly rectilinear capture keeps the pixel
                // path, where the pinhole assumption is exact and results are
                // long settled.
                bc = precomputeBearings(feats, percam, /*wide_only=*/!cs.mixed());
                vopt.bearings = &bc;
                calib->used_bearings = true;
                if (verbose) {
                    fprintf(stderr, "[match] calibrated verification on bearings (%.1f s, %.0f MB)",
                            now() - t_b, bc.bytes() / 1048576.0);
                    for (const auto& kv : cs.cameras)
                        fprintf(stderr, " | cam %u f=%.1f", kv.first, kv.second.focal());
                    fprintf(stderr, "\n");
                }
            }
        }
        // Hand the setup on through the database (D47), so `ssplat-sfm map` inherits
        // the grouping and the focals this stage measured instead of
        // re-deriving them from the inliers it is about to produce.
        if (calib) storeCameraSetup(db, calib->cameras);
        if (verbose)
            fprintf(stderr, "[match] verifying on %d thread(s)\n", verificationThreadCount(vopt));
        db.pairs = verifyPairs(feats, pairs, matchFn, vopt, &stats.putative, progress);
        for (const TwoViewMatches& tvm : db.pairs) stats.inliers += tvm.matches.size();
    } else {
        const size_t batch = std::max(1, opt.batch_pairs);
        std::vector<std::vector<FeatureMatch>> mout;
        for (size_t b = 0; b < pairs.size(); b += batch) {
            size_t e = std::min(b + batch, pairs.size());
            matchFn(b, e, mout);
            for (size_t p = b; p < e; p++) {
                uint32_t i = pairs[p].first, j = pairs[p].second;
                std::vector<FeatureMatch>& m = mout[p - b];
                stats.putative += m.size();
                if (!m.empty()) {
                    stats.inliers += m.size();
                    db.pairs.push_back({i, j, 0, std::move(m)});
                }
                if (progress) progress(p + 1, pairs.size());
            }
        }
    }
    stats.kept = db.pairs.size();
    if (verbose) fprintf(stderr, "\n");
    return 0;
}

static int cmdMatch(int argc, char** argv) {
    std::string featdir, output;
    MatchOptions opt;
    PairMode mode = PairMode::Exhaustive;
    PairSelectionOptions popt;
    int overlap = 10;
    bool verbose = true, verify = true;
    int threads = 0;  // 0 = hardware_concurrency
    TwoViewOptions tvopt;
    VerifyCalibration calib;
    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--output" || a == "-o") output = next();
        else if (a == "--ratio") opt.max_ratio = std::stof(next());
        else if (a == "--no-cross-check") opt.cross_check = false;
        else if (a == "--max-matches") opt.max_num_matches = (uint32_t)std::stoul(next());
        else if (a == "--pairs") {
            std::string m = next();
            if (m == "sequential") mode = PairMode::Sequential;
            else if (m == "prefilter") mode = PairMode::Prefilter;
            else if (m == "exhaustive") mode = PairMode::Exhaustive;
            else { fprintf(stderr, "unknown --pairs %s\n", m.c_str()); return 1; }
        }
        else if (a == "--overlap") overlap = std::stoi(next());
        else if (a == "--prefilter-features") popt.num_features = (uint32_t)std::stoul(next());
        else if (a == "--prefilter-train") popt.train_features = (uint32_t)std::stoul(next());
        else if (a == "--prefilter-neighbors") popt.num_neighbors = (uint32_t)std::stoul(next());
        else if (a == "--prefilter-min-score") popt.min_score = (uint32_t)std::stoul(next());
        else if (a == "--prefilter-ratio") popt.ratio = std::stof(next());
        else if (a == "--no-verify") verify = false;
        else if (a == "--camera-model") {
            std::string v = next();
            if (!parseCameraOverride(v, /*is_focal=*/false, calib.setup.overrides)) {
                fprintf(stderr, "bad --camera-model %s (MODEL or PREFIX=MODEL)\n", v.c_str());
                return 1;
            }
            if (v.find('=') == std::string::npos) parseCamModel(v, calib.setup.model);
        }
        else if (a == "--camera-mode") {
            if (!parseCameraMode(next(), calib.setup.mode)) { fprintf(stderr, "bad --camera-mode\n"); return 1; }
            calib.setup.mode_explicit = true;
        }
        else if (a == "--focal") {
            std::string v = next();
            if (!parseCameraOverride(v, /*is_focal=*/true, calib.setup.overrides)) {
                fprintf(stderr, "bad --focal %s (F or PREFIX=F)\n", v.c_str());
                return 1;
            }
            if (v.find('=') == std::string::npos) calib.setup.focal = std::stod(v);
        }
        else if (a == "--no-exif-focal") calib.setup.exif_focal = false;
        else if (a == "--exif-groups") calib.setup.exif_groups = true;
        else if (a == "--no-exif-groups") calib.setup.exif_groups = false;
        else if (a == "--max-error") tvopt.ransac.max_error = std::stod(next());
        else if (a == "--min-inliers") tvopt.min_num_inliers = std::stoi(next());
        else if (a == "--threads") threads = std::stoi(next());
        else if (a == "--device") opt.device = std::stoi(next());
        else if (a == "--quiet") verbose = false;
        else if (a[0] != '-') featdir = a;
        else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    if (featdir.empty()) {
        fprintf(stderr, "usage: ssplat-sfm match <features_dir> --output matches.bin\n"
                        "  [--pairs exhaustive|sequential|prefilter] [--overlap N] [--ratio X]\n"
                        "  [--prefilter-features K] [--prefilter-neighbors N] [--prefilter-min-score S]\n"
                        "  [--no-cross-check] [--max-matches N] [--no-verify]\n"
                        "  [--max-error PX]  inlier radius in pixels of the image SIFT\n"
                        "     ran on, not of the source file (default 3, D47)\n"
                        "  [--min-inliers N] [--threads N]\n"
                        "  [--camera-model M|PREFIX=M] [--camera-mode single|folder|image]\n"
                        "  [--focal F|PREFIX=F] [--no-exif-focal] [--no-exif-groups]\n"
                        "  [--device I] [--quiet]\n"
                        "Geometric verification (F/H RANSAC) runs by default; --no-verify\n"
                        "keeps raw putative matches. Verification is parallel over pairs\n"
                        "(--threads 0 = all cores, 1 = serial); results do not depend on it.\n"
                        "A fisheye --camera-model switches verification to unit bearings,\n"
                        "where the epipolar constraint holds at any field of view; the focal\n"
                        "is searched per camera group on a sample of pairs unless --focal or\n"
                        "EXIF gives one (D45/D46). --camera-model and --focal take either a\n"
                        "dataset-wide value or PREFIX=VALUE naming one group by path, so a\n"
                        "rig can mix models: --camera-model opencv --camera-model\n"
                        "cam0=thin-prism-fisheye --focal cam0=520\n"
                        "--pairs prefilter scores every pair by mini-matching the K (256)\n"
                        "largest-scale features and fully matches only each image's N (32)\n"
                        "best-scoring partners -- the escape from exhaustive O(N^2) matching.\n");
        return 1;
    }
    popt.device = opt.device;
    popt.batch_pairs = std::max(popt.batch_pairs, opt.batch_pairs);

    std::vector<FeatureSet> feats;
    MatchesDatabase db;
    MatchStats stats;
    if (int rc = matchFeatureDir(featdir, feats, db, opt, mode, overlap, popt, verify, tvopt,
                                 threads, verbose, stats, &calib))
        return rc;
    if (verify)
        printf("matched %zu/%zu pairs (>=1 inlier), %llu inlier / %llu putative matches\n",
               stats.kept, stats.pairs, (unsigned long long)stats.inliers,
               (unsigned long long)stats.putative);
    else
        printf("matched %zu/%zu pairs with >=1 match, %llu total matches\n", stats.kept,
               stats.pairs, (unsigned long long)stats.inliers);
    if (!output.empty()) {
        writeMatches(output, db);
        if (verbose) fprintf(stderr, "[match] wrote %s\n", output.c_str());
    }
    return 0;
}

// -----------------------------------------------------------------------
// selftest
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
static int cmdMap(int argc, char** argv) {
    std::string matchesPath, featdir, output, imagedir, resume;
    bool audit_first = false, check_only = false;
    MapperOptions opt;
    ManagerOptions mgopt;
    // CLI defaults, which are deliberately not MapperOptions' library defaults:
    // per-folder grouping and OpenCV distortion are what a real capture wants
    // (D40), while the library keeps the minimal RADIAL model so the synthetic
    // self-tests stay on the model they were written against.
    CameraSetupOptions sopt;
    sopt.model = opt.camera_model = CamModel::OpenCV;
    // Did the user say anything about cameras? If not, the setup recorded by
    // verification wins over anything this command would derive (D47).
    bool cam_args = false;
    bool final_pp = kFinalPrincipalPointDefault;   // one PP-free global BA at the end (D51)
    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--output" || a == "-o") output = next();
        else if (a == "--camera-mode") {
            if (!parseCameraMode(next(), sopt.mode)) { fprintf(stderr, "bad --camera-mode\n"); return 1; }
            sopt.mode_explicit = true;
            cam_args = true;
        }
        else if (a == "--camera-model") {
            std::string v = next();
            if (!parseCameraOverride(v, /*is_focal=*/false, sopt.overrides)) {
                fprintf(stderr, "bad --camera-model %s (MODEL or PREFIX=MODEL)\n", v.c_str());
                return 1;
            }
            if (v.find('=') == std::string::npos) parseCamModel(v, sopt.model);
            cam_args = true;
        }
        else if (a == "--features") featdir = next();
        else if (a == "--images") imagedir = next();
        else if (a == "--focal") {
            std::string v = next();
            if (!parseCameraOverride(v, /*is_focal=*/true, sopt.overrides)) {
                fprintf(stderr, "bad --focal %s (F or PREFIX=F)\n", v.c_str());
                return 1;
            }
            if (v.find('=') == std::string::npos) sopt.focal = std::stod(v);
            cam_args = true;
        }
        else if (a == "--no-exif-focal") { sopt.exif_focal = false; cam_args = true; }
        else if (a == "--exif-groups") { sopt.exif_groups = true; cam_args = true; }
        else if (a == "--no-exif-groups") { sopt.exif_groups = false; cam_args = true; }
        else if (a == "--max-error") opt.max_reproj_error = std::stod(next());
        else if (a == "--focal-trials") opt.focal_trials = std::stoi(next());
        else if (a == "--refine-principal-point") opt.refine_principal_point = true;
        else if (a == "--no-final-principal-point") final_pp = false;
        else if (a == "--final-principal-point") final_pp = true;
        else if (a == "--pp-min-images") opt.pp_min_images = (size_t)std::stoul(next());
        else if (a == "--min-tri-angle") opt.min_tri_angle_deg = std::stod(next());
        else if (a == "--min-pnp-inliers") opt.min_num_pnp_inliers = std::stoi(next());
        else if (a == "--min-pnp-ratio") opt.min_pnp_inlier_ratio = std::stod(next());
        else if (a == "--min-image-points") opt.min_image_points = std::stoi(next());
        else if (a == "--ba-loss") opt.ba_loss = next();
        else if (a == "--ba-loss-param") opt.ba_loss_param = std::stod(next());
        else if (a == "--retri-scale") opt.retri_scale = std::stod(next());
        else if (a == "--ba-growth") opt.ba_growth_ratio = std::stod(next());
        else if (a == "--ba-growth-rtol") opt.ba_growth_rtol = std::stod(next());
        else if (a == "--ba-growth-patience") opt.ba_growth_patience = std::stoi(next());
        else if (a == "--max-models") opt.max_num_models = std::stoi(next());
        else if (a == "--model-overlap") opt.max_model_overlap = std::stoi(next());
        else if (a == "--min-model-size") opt.min_model_size = std::stoi(next());
        else if (a == "--resume") resume = next();
        else if (a == "--audit") audit_first = true;
        else if (a == "--audit-evidence") opt.audit_min_evidence = std::stoi(next());
        else if (a == "--rounds") mgopt.max_rounds = std::stoi(next());
        else if (a == "--check") check_only = true;
        else if (a == "--no-manage") {
            mgopt.do_merge = mgopt.do_grow = mgopt.do_reseed = mgopt.do_audit = false;
        }
        else if (a == "--no-merge") mgopt.do_merge = false;
        else if (a == "--no-grow") mgopt.do_grow = false;
        else if (a == "--no-reseed") mgopt.do_reseed = false;
        else if (a == "--no-audit") mgopt.do_audit = false;
        else if (a == "--no-split") mgopt.do_split = false;
        else if (a == "--fold-split") mgopt.do_duplicate_split = true;
        else if (a == "--no-fold-split") mgopt.do_duplicate_split = false;
        else if (a == "--fold-max-cut") mgopt.duplicate.max_cut_fraction = std::stod(next());
        else if (a == "--no-joint-ba") mgopt.do_joint_ba = false;
        else if (a == "--device") opt.device = std::stoi(next());
        else if (a == "--quiet") opt.verbose = false;
        else if (a[0] != '-' && matchesPath.empty()) matchesPath = a;
        else if (a[0] != '-' && featdir.empty()) featdir = a;
        else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    if (matchesPath.empty() || featdir.empty()) {
        fprintf(stderr, "usage: ssplat-sfm map <matches.bin> <features_dir> --output sparse/\n"
                        "  [--camera-mode single|folder|image]  (default folder)\n"
                        "  [--camera-model M|PREFIX=M]  (default opencv; one of\n"
                        "     simple-pinhole|pinhole|radial|opencv|full-opencv|opencv-fisheye|\n"
                        "     thin-prism-fisheye|equirectangular) [--focal F|PREFIX=F]\n"
                        "     PREFIX names one camera group by image path, so a rig can mix\n"
                        "  [--refine-principal-point]  (refine cx,cy throughout; off, as in\n"
                        "     COLMAP: it is nearly the same parameter as a camera rotation, D50.\n"
                        "  [--no-final-principal-point] [--pp-min-images N]  one PP-free global\n"
                        "     BA on the finished model, when it has ONE camera group of >= N (20)\n"
                        "     images. On: worth 8-24 AUC where the lens really is decentred,\n"
                        "     ~1 where it is not, and skipped on rigs (D51)\n"
                        "     models and known with unknown focals (D46)\n"
                        "  [--no-exif-focal]  ignore the focal length EXIF recorded\n"
                        "  [--no-exif-groups]  do not split camera groups by what EXIF says\n"
                        "     the camera and its focal setting were\n"
                        "  [--max-error PX]  inlier radius in pixels of the image SIFT\n"
                        "     ran on, not of the source file (default 3, D47)\n"
                        "  [--focal-trials N]  trial reconstructions used to pick the focal\n"
                        "     when the motion cannot determine it (D48); 0 = off, default 5\n"
                        "  [--min-tri-angle DEG]\n"
                        "  [--ba-growth RATIO] [--ba-growth-rtol TOL (0 = solver default)]\n"
                        "  [--max-models N] (default 50, 1 = single model)\n"
                        "  [--model-overlap N] [--min-model-size N]\n"
                        "  [--resume DIR]  adopt the models in DIR instead of mapping from\n"
                        "     scratch, and run the merge/grow loop on them (D44)\n"
                        "  [--rounds N] [--no-manage] [--no-merge] [--no-grow] [--no-reseed]\n"
                        "  [--audit] [--no-audit] [--audit-evidence N]\n"
                        "  [--no-split] [--no-joint-ba]\n"
                        "  [--no-fold-split]  keep a model that has written two places on\n"
                        "     top of each other; the split is on by default and is taken only\n"
                        "     when it severs almost no co-visibility (D46), which is what\n"
                        "     separates a fold from a dense capture. `--check` reports both\n"
                        "  [--check]  with --resume: report how well each model agrees with\n"
                        "     the verified two-view geometries it was built from, and what\n"
                        "     splitting it would do; changes nothing and writes nothing\n"
                        "  [--images DIR] [--device I] [--quiet]\n"
                        "A capture that is not one connected view graph reconstructs as\n"
                        "several models, written to <output>/0, <output>/1, ... largest\n"
                        "first (COLMAP's layout). They are then merged where they share\n"
                        "images and grown into whatever else they can register.\n");
        return 1;
    }
    MatchesDatabase db = readMatches(matchesPath);
    std::vector<FeatureSet> feats(db.images.size());
    for (size_t i = 0; i < db.images.size(); i++)
        feats[i] = readFeatures(featdir + "/" + db.images[i].name + ".bin");

    // The camera setup, in order of authority: what the command line asked for,
    // else what verification recorded in matches.bin (D47), else derived here.
    CameraSetup cs;
    const bool from_db = !cam_args && loadCameraSetup(db, cs);
    if (!from_db) cs = buildCameras(db.images, feats, sopt);
    if (opt.verbose) {
        printCameraSetup("map", cs, sopt, db.images.size());
        if (from_db)
            fprintf(stderr, "[map] camera setup taken from %s (as verified)\n",
                    matchesPath.c_str());
        else if (cam_args && db.hasCameras())
            fprintf(stderr, "[map] camera setup rebuilt from the command line, ignoring the "
                            "one in %s\n", matchesPath.c_str());
    }

    // A fisheye group with no focal prior and none recorded (D45/D46/D47).
    // `ssplat-sfm auto` and `ssplat-sfm match` measure this before verifying, on raw putative
    // matches, and it now travels in the match database; this is the fallback
    // for a matches.bin that predates that or arrived from elsewhere. The
    // sample here is verified *inliers*, so if the verification was itself
    // calibrated the answer is pulled towards the focal it used -- which is
    // exactly why carrying the measured one is better than re-deriving it.
    if (!from_db && db.pairs.size() >= 8) {
        std::vector<std::pair<uint32_t, uint32_t>> sample;
        std::vector<std::vector<FeatureMatch>> sm;
        const size_t want = 200 * std::max<size_t>(1, cs.count());
        size_t stride = std::max<size_t>(1, db.pairs.size() / want);
        for (size_t p = 0; p < db.pairs.size() && sample.size() < want; p += stride) {
            sample.push_back({db.pairs[p].image1, db.pairs[p].image2});
            sm.push_back(db.pairs[p].matches);
        }
        bootstrapGroupFocals(feats, cs.ids, sample, sm, cs.cameras, cs.focal_given,
                             TwoViewOptions{}, 200, 0, opt.verbose);
    }
    opt.initial_cameras = cs.cameras;
    opt.known_focal_cameras = cs.focal_known;
    opt.given_focal_cameras = cs.focal_given;

    Mapper mapper(db, feats, opt, cs.ids);
    std::vector<Reconstruction> models;
    if (resume.empty()) {
        models = mapper.run();
    } else {
        // Adopt what a previous run wrote and work on it instead. The models
        // must come from this database (image ids are positions in it); adopt()
        // checks the names and says so if they do not.
        if (!readModels(resume, models, opt.verbose)) return 1;
        printf("resumed %zu model(s) covering %zu/%zu images\n", models.size(),
               distinctRegistered(models), db.images.size());
    }
    if (models.empty()) {
        fprintf(stderr, "map: no model to work with\n");
        return 1;
    }

    // Diagnostic (D45): does each model agree with the two-view geometries it
    // was built from? A model welded together at a wrong relative pose is
    // internally perfect and fails no other test, but the verified pairs that
    // span the weld do not hold, and the images fall into disconnected groups.
    if (check_only) {
        for (size_t i = 0; i < models.size(); i++) {
            Mapper::SplitStats ss;
            std::vector<Reconstruction> parts = mapper.splitInconsistent(
                models[i], mgopt.merge.max_reproj_error, mgopt.seam_min_pair_fraction,
                mgopt.split_min_matches, mgopt.split_min_group, &ss);
            printf("model %zu: %u images, %zu points | %zu/%zu inner pairs hold (%.0f%%)\n"
                   "    pair agreement p05/p10/p25/p50 %.2f/%.2f/%.2f/%.2f | "
                   "%zu group(s), largest %zu",
                   i, models[i].numRegistered(), models[i].points3D.size(), ss.pairs_agree,
                   ss.pairs_tested, ss.pairs_tested ? 100.0 * ss.pairs_agree / ss.pairs_tested : 0.0,
                   ss.percentile(0.05), ss.percentile(0.10), ss.percentile(0.25),
                   ss.percentile(0.50), ss.groups, ss.largest);
            if (parts.size() > 1) {
                printf(" -> would split into");
                for (const Reconstruction& p : parts) printf(" %u", p.numRegistered());
                if (ss.dropped_images) printf(" (%zu images dropped)", ss.dropped_images);
            }
            printf("\n");
            DuplicateReport dr =
                findDuplicateStructure(models[i], mgopt.duplicate, mapper.matchedPredicate());
            printf("    duplicate structure: %zu of %zu co-located pairs share no points "
                   "and never matched (%.0f%%; %zu more share none but did match)\n",
                   dr.conflicts, dr.colocated, 100.0 * dr.ratio(), dr.unmatched_but_seen);
            // The cut is reported whenever there are conflicts at all, not only
            // when the rate passes: the whole point is that the rate is not the
            // criterion (D46), so both numbers have to be visible side by side.
            if (!dr.pairs.empty()) {
                size_t dropped = 0;
                DuplicateCut cut;
                std::vector<Reconstruction> dp = splitDuplicateStructure(
                    models[i], dr, mgopt.split_min_group, &dropped, &cut);
                printf("    cut: %zu group(s), severs %.2f%% of co-visibility (%llu of %llu)%s",
                       cut.groups, 100.0 * cut.fraction(), (unsigned long long)cut.severed,
                       (unsigned long long)(cut.severed + cut.kept),
                       foldSplitAccepted(dr, cut, mgopt.duplicate) ? " -- folded" : " -- kept whole");
                if (dp.size() > 1) {
                    printf(" -> would split into");
                    for (const Reconstruction& p : dp) printf(" %u", p.numRegistered());
                    if (dropped) printf(" (%zu images dropped)", dropped);
                }
                printf("\n");
            }
        }
        return 0;
    }

    // Audit models that came from elsewhere before doing anything with them:
    // every image's pose is checked against the correspondence graph, and the
    // ones the model cannot support are re-registered (D44).
    if (audit_first) {
        for (size_t i = 0; i < models.size(); i++) {
            Mapper::AuditStats as;
            const uint32_t before = models[i].numRegistered();
            models[i] = mapper.audit(models[i], &as);
            printf("audit model %zu: %u images checked, %u unsupported, %u re-registered, "
                   "%u dropped -> %u images\n",
                   i, as.checked, as.unsupported, as.reregistered, as.deregistered,
                   models[i].numRegistered());
            (void)before;
        }
    }

    ManagerStats mst;
    double t_manage = 0;
    mgopt.verbose = opt.verbose;
    mgopt.merge.verbose = opt.verbose;
    mgopt.merge.filter_reproj_error = opt.max_reproj_error;
    mgopt.merge.min_tri_angle_deg = opt.min_tri_angle_deg;
    mgopt.merge.min_image_points = opt.min_image_points;
    if (mgopt.do_merge || mgopt.do_grow || mgopt.do_reseed || mgopt.do_split ||
        mgopt.do_duplicate_split)
        models = manageModels(mapper, std::move(models), mgopt, mst, t_manage);
    if (final_pp) {
        double t_pp = 0;
        models = polishModels(mapper, std::move(models), opt.verbose, t_pp);
    }
    if (mst.rounds)
        printf("manage: %zu round(s), %zu merge(s), %zu image(s) registered by growth, "
               "%zu repaired / %zu dropped by the audit, %zu model(s) covering %zu -> %zu "
               "images, %.1f s\n",
               mst.rounds, mst.merges, mst.grown_images, mst.audited_repaired, mst.audited_out,
               models.size(), mst.covered_before, mst.covered_after, t_manage);

    resolveImageNames(models, imagedir);

    const Reconstruction& rec = models.front();
    double mean = 0, median = 0;
    size_t nobs = 0;
    reprojStats(rec, feats, mean, median, nobs);
    printf("reconstruction: %u/%zu images registered, %zu 3D points, "
           "%.3f px mean reprojection (%.3f median, %zu obs)\n",
           rec.numRegistered(), db.images.size(), rec.points3D.size(), mean, median, nobs);
    if (models.size() > 1) {
        printf("  %zu models covering %zu/%zu distinct images:\n", models.size(),
               distinctRegistered(models), db.images.size());
        printExtraModels(models, feats, "  ");
    }
    if (!output.empty()) writeModels(models, output, opt.verbose);
    return 0;
}

// -----------------------------------------------------------------------
// merge: fold the models of a fragmented capture back together (D43)
// -----------------------------------------------------------------------

static int cmdMerge(int argc, char** argv) {
    std::vector<std::string> inputs;
    std::string output;
    MergeOptions mo;
    bool refine = true, in_place = false;
    int device = -1;
    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--output" || a == "-o") output = next();
        else if (a == "--max-error") mo.max_reproj_error = std::stod(next());
        else if (a == "--min-common") mo.min_common_images = std::stoi(next());
        else if (a == "--min-inlier-ratio") mo.min_inlier_ratio = std::stod(next());
        else if (a == "--filter-error") mo.filter_reproj_error = std::stod(next());
        else if (a == "--min-tri-angle") mo.min_tri_angle_deg = std::stod(next());
        else if (a == "--no-ba") refine = false;
        else if (a == "--in-place") in_place = true;
        else if (a == "--device") device = std::stoi(next());
        else if (a == "--quiet") mo.verbose = false;
        else if (a[0] != '-') inputs.push_back(a);
        else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    if (inputs.empty() || (output.empty() && !in_place)) {
        fprintf(stderr,
                "usage: ssplat-sfm merge <sparse_dir | model_dir> [more...] --output DIR\n"
                "  a directory holding 0/, 1/, ... is expanded to all of them\n"
                "  [--max-error PX]       alignment inlier threshold (default 8)\n"
                "  [--min-common N]       shared images an alignment needs (default 3)\n"
                "  [--min-inlier-ratio R] of the shared images (default 0.3)\n"
                "  [--filter-error PX] [--min-tri-angle DEG]  post-merge filtering\n"
                "  [--no-ba]              skip the bundle adjustment across the seams\n"
                "  [--in-place]           write back over the input directory\n"
                "  [--device I] [--quiet]\n"
                "Merges the models of a fragmented capture (D41) on the images they\n"
                "share: those give the similarity transform between the two gauges.\n"
                "Models built from different features cannot be merged.\n");
        return 1;
    }

    std::vector<fs::path> dirs;
    for (const std::string& in : inputs)
        if (!collectModelDirs(in, dirs)) return 1;
    if (in_place) {
        if (inputs.size() != 1) {
            fprintf(stderr, "merge: --in-place takes exactly one input directory\n");
            return 1;
        }
        if (output.empty()) output = inputs[0];
    }
    // Writing merged models over their own inputs is destructive (model 0 is
    // replaced and the absorbed ones are removed), so it has to be asked for.
    if (!in_place)
        for (const fs::path& d : dirs)
            if (fs::exists(output) && fs::equivalent(fs::path(output), d.parent_path())) {
                fprintf(stderr,
                        "merge: --output %s is where the input models live; pass --in-place if "
                        "that is what you want\n",
                        output.c_str());
                return 1;
            }

    std::vector<Reconstruction> models;
    for (const fs::path& d : dirs) {
        try {
            models.push_back(Reconstruction::readBinary(d.string()));
        } catch (const std::exception& e) {
            fprintf(stderr, "merge: cannot read %s: %s\n", d.string().c_str(), e.what());
            return 1;
        }
        printf("   read %s: %u images, %zu points\n", d.string().c_str(),
               models.back().numRegistered(), models.back().points3D.size());
    }
    if (models.size() < 2) {
        fprintf(stderr, "merge: need at least two models, found %zu\n", models.size());
        return 1;
    }
    const size_t covered_before = distinctRegistered(models);

    MergeSummary sum;
    models = mergeModels(std::move(models), mo, refine, device, sum);

    printf("merged %zu -> %zu model(s) in %.2f s (%zu merge(s), %zu refused",
           sum.before, sum.after, sum.seconds, sum.merges, sum.refused);
    if (sum.ba_seconds > 0) printf(", %.2f s bundle adjustment", sum.ba_seconds);
    printf(")\n");
    for (size_t i = 0; i < models.size(); i++) {
        double mean = 0, median = 0;
        size_t nobs = 0;
        modelReprojStats(models[i], mean, median, nobs);
        printf("  model %zu: %u images, %zu points, %.3f px mean (%.3f median, %zu obs)\n", i,
               models[i].numRegistered(), models[i].points3D.size(), mean, median, nobs);
    }
    const size_t covered_after = distinctRegistered(models);
    if (covered_after != covered_before)
        printf("  NOTE: %zu/%zu distinct images survived the merge\n", covered_after,
               covered_before);

    writeModels(models, fs::path(output), mo.verbose);
    // In place, the models that were absorbed must not stay behind as stale
    // directories claiming to be reconstructions.
    if (in_place)
        for (const fs::path& d : dirs) {
            const std::string name = d.filename().string();
            if (d.parent_path() != fs::path(output) || name.empty() ||
                name.find_first_not_of("0123456789") != std::string::npos)
                continue;
            long idx = std::stol(name);
            if (idx >= (long)models.size()) {
                fs::remove_all(d);
                printf("  removed %s (absorbed)\n", d.string().c_str());
            }
        }
    return 0;
}

// -----------------------------------------------------------------------
// auto: one command from an image directory to a COLMAP sparse model
// -----------------------------------------------------------------------
//
// The equivalent of COLMAP's `automatic_reconstructor`: pick sane settings from
// two knobs (quality, data type), run extract -> match -> map in one process,
// and lay the workspace out the way COLMAP does so existing tooling reads it.
//
// Presets mirror COLMAP's OptionManager::ModifyFor*Quality() as fractions of
// the SIFT extractor's 3200 px default, and ModifyForVideoData()'s halved
// initial triangulation angle. Above COLMAP's 200-image exhaustive cutoff --
// where COLMAP switches to vocabulary-tree retrieval -- we switch to GPU pair
// selection (sfm/feature/PairSelection.h, D35) unless --pairs was given explicitly.
enum class Quality { Low, Medium, High, Extreme };
enum class DataType { Individual, Video, Internet };

static int cmdAuto(int argc, char** argv) {
    std::string imagedir, maskdir, workspace;
    bool maskdir_explicit = false;
    Quality quality = Quality::High;  // COLMAP's default
    DataType data_type = DataType::Individual;
    CameraSetupOptions sopt_cam;
    bool cmode_explicit = false;
    int max_num_models = MapperOptions().max_num_models;
    // Merging is on by default (D43): a fragmented capture's models are only
    // separate because the mapper had no way to relate their gauges, and the
    // images they share are exactly that relation. It is a no-op on a capture
    // that reconstructs as one model.
    bool merge = true;
    MergeOptions mopt_merge;
    ManagerOptions mgopt;
    int overlap = 10, threads = 0, decode_threads = 0, decode_budget_mb = 0, device = -1;
    // The one geometric tolerance, in extraction pixels (D47): the two-view
    // verifier's inlier radius and the mapper's reprojection cap are the same
    // quantity measured in the same frame, and moving them apart has never
    // been what a caller wanted.
    double max_error = 0;   // 0 = each stage's own default
    int focal_trials = -1;  // -1 = the mapper's default
    bool refine_pp = false;
    bool final_pp = kFinalPrincipalPointDefault;   // one PP-free global BA at the end (D51)
    size_t pp_min_images = MapperOptions().pp_min_images;
    bool verbose = true;
    PairMode pairs_override = PairMode::Exhaustive;
    bool pairs_explicit = false;
    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&]() { return std::string(argv[++i]); };
        if (a == "--output" || a == "-o") workspace = next();
        else if (a == "--pairs") {
            std::string m = next();
            if (m == "sequential") pairs_override = PairMode::Sequential;
            else if (m == "prefilter") pairs_override = PairMode::Prefilter;
            else if (m == "exhaustive") pairs_override = PairMode::Exhaustive;
            else { fprintf(stderr, "unknown --pairs %s\n", m.c_str()); return 1; }
            pairs_explicit = true;
        }
        else if (a == "--quality") {
            std::string q = next();
            if (q == "low") quality = Quality::Low;
            else if (q == "medium") quality = Quality::Medium;
            else if (q == "high") quality = Quality::High;
            else if (q == "extreme") quality = Quality::Extreme;
            else { fprintf(stderr, "unknown --quality %s\n", q.c_str()); return 1; }
        }
        else if (a == "--data-type") {
            std::string d = next();
            if (d == "individual") data_type = DataType::Individual;
            else if (d == "video") data_type = DataType::Video;
            else if (d == "internet") data_type = DataType::Internet;
            else { fprintf(stderr, "unknown --data-type %s\n", d.c_str()); return 1; }
        }
        else if (a == "--camera-mode") {
            if (!parseCameraMode(next(), sopt_cam.mode)) { fprintf(stderr, "bad --camera-mode\n"); return 1; }
            sopt_cam.mode_explicit = true;
            cmode_explicit = true;
        }
        else if (a == "--camera-model") {
            std::string v = next();
            if (!parseCameraOverride(v, /*is_focal=*/false, sopt_cam.overrides)) {
                fprintf(stderr, "bad --camera-model %s (MODEL or PREFIX=MODEL)\n", v.c_str());
                return 1;
            }
            if (v.find('=') == std::string::npos) parseCamModel(v, sopt_cam.model);
        }
        else if (a == "--focal") {
            std::string v = next();
            if (!parseCameraOverride(v, /*is_focal=*/true, sopt_cam.overrides)) {
                fprintf(stderr, "bad --focal %s (F or PREFIX=F)\n", v.c_str());
                return 1;
            }
            if (v.find('=') == std::string::npos) sopt_cam.focal = std::stod(v);
        }
        else if (a == "--no-exif-focal") sopt_cam.exif_focal = false;
        else if (a == "--exif-groups") sopt_cam.exif_groups = true;
        else if (a == "--no-exif-groups") sopt_cam.exif_groups = false;
        else if (a == "--max-error") max_error = std::stod(next());
        else if (a == "--focal-trials") focal_trials = std::stoi(next());
        else if (a == "--refine-principal-point") refine_pp = true;
        else if (a == "--no-final-principal-point") final_pp = false;
        else if (a == "--final-principal-point") final_pp = true;
        else if (a == "--pp-min-images") pp_min_images = (size_t)std::stoul(next());
        else if (a == "--masks" || a == "--mask-dir") { maskdir = next(); maskdir_explicit = true; }
        else if (a == "--no-masks") { maskdir.clear(); maskdir_explicit = true; }
        else if (a == "--max-models") max_num_models = std::stoi(next());
        else if (a == "--no-merge") merge = false;
        else if (a == "--merge-max-error") mopt_merge.max_reproj_error = std::stod(next());
        else if (a == "--merge-min-common") mopt_merge.min_common_images = std::stoi(next());
        else if (a == "--rounds") mgopt.max_rounds = std::stoi(next());
        else if (a == "--no-grow") mgopt.do_grow = false;
        else if (a == "--no-reseed") mgopt.do_reseed = false;
        else if (a == "--no-audit") mgopt.do_audit = false;
        else if (a == "--no-split") mgopt.do_split = false;
        else if (a == "--fold-split") mgopt.do_duplicate_split = true;
        else if (a == "--no-fold-split") mgopt.do_duplicate_split = false;
        else if (a == "--fold-max-cut") mgopt.duplicate.max_cut_fraction = std::stod(next());
        else if (a == "--no-joint-ba") mgopt.do_joint_ba = false;
        else if (a == "--overlap") overlap = std::stoi(next());
        else if (a == "--threads") threads = std::stoi(next());
        else if (a == "--decode-threads") decode_threads = std::stoi(next());
        else if (a == "--decode-budget") decode_budget_mb = std::stoi(next());
        else if (a == "--device") device = std::stoi(next());
        else if (a == "--quiet") verbose = false;
        else if (a[0] != '-' && imagedir.empty()) imagedir = a;
        else { fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    if (workspace.empty()) {
        fprintf(stderr,
                "usage: ssplat-sfm auto [image_dir|dataset_dir] --output WORKSPACE\n"
                "  the positional defaults to `images`; if it names a dataset directory\n"
                "  containing `images/`, that sub-directory is used\n"
                "  [--masks DIR]   (default: `masks` beside the image directory)\n"
                "  [--no-masks]    ignore a masks directory even if one is there\n"
                "  [--quality low|medium|high|extreme]   (default high)\n"
                "  [--data-type individual|video|internet]  (default individual)\n"
                "     internet defaults --camera-mode to `image`\n"
                "  [--camera-mode single|folder|image]  (default folder)\n"
                "     single = one camera per distinct image resolution\n"
                "     folder = one camera per image sub-folder (and per resolution)\n"
                "     image  = one camera per image\n"
                "  [--camera-model M|PREFIX=M]  (default opencv; one of\n"
                "     simple-pinhole|pinhole|radial|opencv|full-opencv|opencv-fisheye|\n"
                "     thin-prism-fisheye|equirectangular) [--focal F|PREFIX=F]\n"
                "     PREFIX names one camera group by image path, so a rig can mix\n"
                "  [--refine-principal-point]  (refine cx,cy throughout; off, as in\n"
                "     COLMAP: it is nearly the same parameter as a camera rotation, D50.\n"
                "  [--no-final-principal-point] [--pp-min-images N]  one PP-free global\n"
                "     BA on the finished model, when it has ONE camera group of >= N (20)\n"
                "     images. On: worth 8-24 AUC where the lens really is decentred,\n"
                "     ~1 where it is not, and skipped on rigs (D51)\n"
                "     models and known with unknown focals (D46)\n"
                "  [--no-exif-focal]  ignore the focal length EXIF recorded\n"
                "  [--no-exif-groups]  do not split camera groups by what EXIF says\n"
                "     the camera and its focal setting were\n"
                "  [--max-error PX]  inlier radius for verification and mapping, in\n"
                "     pixels of the image SIFT ran on (D47); default 3\n"
                "  [--focal-trials N]  trial reconstructions used to pick the focal when\n"
                "     the capture's motion cannot determine it (D48); 0 = off, default 5\n"
                "  [--pairs exhaustive|sequential|prefilter]\n"
                "     default: video -> sequential; otherwise exhaustive below 100 images\n"
                "     and prefilter (GPU pair selection) at or above\n"
                "  [--max-models N]  (default 50, 1 = one model only)\n"
                "  [--no-merge]      keep the components separate instead of merging\n"
                "     the ones that share images (see `ssplat-sfm merge`)\n"
                "  [--merge-max-error PX] [--merge-min-common N]\n"
                "  [--rounds N] [--no-grow] [--no-reseed] [--no-audit]\n"
                "  [--no-split]     keep models the correspondence graph contradicts\n"
                "  [--no-fold-split]  keep a model that has written two places on top\n"
                "     of each other; cutting it is on by default (D46) and is taken only\n"
                "     when it severs almost no co-visibility -- see `ssplat-sfm map --check`\n"
                "  [--no-joint-ba]  refine each component separately instead of\n"
                "     sharing one set of intrinsics per camera group (D45)\n"
                "  [--overlap N] [--threads N] [--decode-threads N] [--decode-budget MB]\n"
                "  [--device I] [--quiet]\n"
                "Runs extract -> match -> map with COLMAP's automatic_reconstructor\n"
                "presets and writes WORKSPACE/{features,matches.bin,sparse/0}. A capture\n"
                "that is not one connected view graph also writes sparse/1, sparse/2, ...\n"
                "-- the remaining components, largest first, for a later merge.\n");
        return 1;
    }

    // ---- where the images and masks are (D39/D40) ----
    // The default layout is a dataset directory holding `images/` and `masks/`,
    // which is spirulae-splat's and nerfstudio's. Pointing straight at an image
    // directory still works: `masks` is then looked for as its sibling, so
    // `ssplat-sfm auto DATASET/images` and `ssplat-sfm auto DATASET` behave the same.
    if (imagedir.empty()) imagedir = "images";
    if (!fs::is_directory(imagedir)) {
        fprintf(stderr, "%s is not a directory\n", imagedir.c_str());
        return 1;
    }
    {
        fs::path nested = fs::path(imagedir) / "images";
        if (fs::is_directory(nested)) {
            printf("   %s contains images/, using %s as the image directory\n", imagedir.c_str(),
                   nested.string().c_str());
            imagedir = nested.string();
        }
    }
    if (!maskdir_explicit) {
        fs::path p = fs::path(imagedir);
        if (!p.empty() && p.filename().empty()) p = p.parent_path();  // drop a trailing '/'
        fs::path sibling = p.parent_path() / "masks";
        if (fs::is_directory(sibling)) maskdir = sibling.string();
    }

    // ---- presets ----
    SiftOptions sopt;
    sopt.device = device;
    sopt.verbose = verbose;
    int max_image_size = 3200;  // the extractor default the fractions apply to
    // Pair-selection breadth also follows the quality knob (D42): `k` is how
    // many partners each image keeps, and it is the one pair-selection
    // parameter that trades match time against how much of the view graph is
    // even offered to verification. `high` keeps sfm/feature/PairSelection.h's 32, so
    // the default path -- and every benchmark, which runs at high -- is
    // unchanged; `low`/`medium` buy back match time, `extreme` spends it.
    uint32_t prefilter_neighbors = PairSelectionOptions().num_neighbors;  // 32
    switch (quality) {
        case Quality::Low:    max_image_size = 1000; sopt.max_num_features = 2048;
                              prefilter_neighbors = 16; break;
        case Quality::Medium: max_image_size = 1600; sopt.max_num_features = 4096;
                              prefilter_neighbors = 24; break;
        case Quality::High:   max_image_size = 2400; sopt.max_num_features = 8192; break;
        case Quality::Extreme: max_image_size = 3200; sopt.max_num_features = 16384;
                              prefilter_neighbors = 48; break;
    }
    MatchOptions mopt;
    mopt.device = device;
    TwoViewOptions tvopt;
    MapperOptions mapopt;
    if (max_error > 0) {
        tvopt.ransac.max_error = max_error;
        mapopt.max_reproj_error = max_error;
    }
    if (focal_trials >= 0) mapopt.focal_trials = focal_trials;
    mapopt.refine_principal_point = refine_pp;
    mapopt.pp_min_images = pp_min_images;
    mapopt.device = device;
    mapopt.verbose = verbose;
    mapopt.camera_model = sopt_cam.model;
    mapopt.max_num_models = max_num_models;
    PairMode mode = PairMode::Exhaustive;
    if (data_type == DataType::Video) {
        mode = PairMode::Sequential;
        mapopt.init_min_tri_angle_deg /= 2;  // COLMAP ModifyForVideoData
    }
    if (pairs_explicit) mode = pairs_override;
    PairSelectionOptions popt;
    popt.device = device;
    popt.num_neighbors = prefilter_neighbors;
    // Camera grouping follows the data type unless asked for explicitly.
    // Internet photos are the case that matters: two Flickr uploads that happen
    // to be 1024x768 are two different cameras that were each downscaled, so
    // grouping them on resolution fuses unrelated intrinsics. Give every image
    // its own (D20).
    if (!cmode_explicit && data_type == DataType::Internet) {
        sopt_cam.mode = CameraMode::Image;
        sopt_cam.mode_explicit = true;   // asked for, just not by --camera-mode
    }

    fs::path ws(workspace);
    fs::create_directories(ws);
    const fs::path featdir = ws / "features";
    const fs::path matchpath = ws / "matches.bin";
    // COLMAP's layout: sparse/<i> per reconstruction, sparse/0 the largest.
    const fs::path sparsedir = ws / "sparse";

    printf("== ssplat-sfm auto: %s -> %s ==\n", imagedir.c_str(), workspace.c_str());
    printf("   quality %s (max image %d px, max %d features), data type %s\n",
           quality == Quality::Low ? "low" : quality == Quality::Medium ? "medium"
               : quality == Quality::High ? "high" : "extreme",
           max_image_size, sopt.max_num_features,
           data_type == DataType::Video      ? "video"
           : data_type == DataType::Internet ? "internet"
                                             : "individual");
    printf("   camera model %s\n", camInfo(sopt_cam.model).cli_name);

    // ---- 1. extract ----
    double t0 = now();
    ExtractStats est;
    if (int rc = extractDirectory(imagedir, featdir, sopt, max_image_size, decode_threads,
                                  decode_budget_mb, maskdir, est))
        return rc;
    double t_extract = now() - t0;
    if (est.images < 2) {
        fprintf(stderr, "auto: only %zu image(s) extracted; need at least 2\n", est.images);
        return 1;
    }
    warnIfMasksLookInverted(est);
    if (mode == PairMode::Exhaustive && est.images >= 100) {
        if (pairs_explicit)
            fprintf(stderr,
                    "[auto] WARNING: %zu images with exhaustive pairing = %zu pairs; this is "
                    "quadratic. Drop --pairs exhaustive to get pair selection.\n",
                    est.images, est.images * (est.images - 1) / 2);
        else {
            mode = PairMode::Prefilter;
            printf("   %zu images >= 100: switching to pair selection "
                   "(exhaustive would be %zu pairs; --pairs exhaustive forces it)\n",
                   est.images, est.images * (est.images - 1) / 2);
        }
    }

    // ---- 2. match + geometric verification ----
    t0 = now();
    std::vector<FeatureSet> feats;
    MatchesDatabase db;
    MatchStats mstats;
    VerifyCalibration calib;
    calib.setup = sopt_cam;
    if (int rc = matchFeatureDir(featdir.string(), feats, db, mopt, mode, overlap, popt, true,
                                 tvopt, threads, verbose, mstats, &calib))
        return rc;
    double t_match = now() - t0;
    writeMatches(matchpath.string(), db);

    // ---- 3. incremental mapping ----
    // The grouping and the focals the two-view stage settled on carry straight
    // into mapping: a focal it measured beats the geometric default the mapper
    // would otherwise start the group from, and every group starting from a
    // measurement is what stops small components inventing their own
    // intrinsics (D45/D46).
    const CameraSetup& cs = calib.cameras;
    printf("   camera mode %s -> %u camera(s)\n", cameraModeName(cs.mode_used), cs.count());
    mapopt.initial_cameras = cs.cameras;
    mapopt.known_focal_cameras = cs.focal_known;
    mapopt.given_focal_cameras = cs.focal_given;

    t0 = now();
    Mapper mapper(db, feats, mapopt, cs.ids);
    std::vector<Reconstruction> models = mapper.run();
    double t_map = now() - t0;

    // ---- 3b. merge, grow, prune, reseed until nothing changes (D43, D44) ----
    ManagerStats mst;
    double t_manage = 0;
    if (merge) {
        mgopt.verbose = verbose;
        mgopt.merge = mopt_merge;
        mgopt.merge.verbose = verbose;
        mgopt.merge.filter_reproj_error = mapopt.max_reproj_error;
        mgopt.merge.min_tri_angle_deg = mapopt.min_tri_angle_deg;
        mgopt.merge.min_image_points = mapopt.min_image_points;
        if (models.size() > 1 || distinctRegistered(models) < est.images)
            models = manageModels(mapper, std::move(models), mgopt, mst, t_manage);
    }
    if (final_pp) {
        double t_pp = 0;
        models = polishModels(mapper, std::move(models), verbose, t_pp);
        t_map += t_pp;
    }

    resolveImageNames(models, imagedir);
    writeModels(models, sparsedir, verbose);

    // ---- report ----
    const Reconstruction& rec = models.front();
    double mean = 0, median = 0;
    size_t nobs = 0;
    reprojStats(rec, feats, mean, median, nobs);
    const uint32_t reg = rec.numRegistered();
    printf("\n== summary ==\n");
    printf("  extract : %6.2f s  %zu images, %llu features\n", t_extract, est.images,
           (unsigned long long)est.features);
    if (est.masked_images) {
        const uint64_t before = est.features + est.masked_out;
        printf("  masks   :          %zu/%zu images masked, %llu keypoints dropped (%.1f%%)",
               est.masked_images, est.images, (unsigned long long)est.masked_out,
               before ? 100.0 * est.masked_out / before : 0.0);
        if (est.unmasked_images || est.mask_unreadable)
            printf(", %zu without a mask, %zu undecodable", est.unmasked_images,
                   est.mask_unreadable);
        printf("\n");
    }
    printf("  match   : %6.2f s  %zu/%zu pairs kept, %llu inlier / %llu putative\n", t_match,
           mstats.kept, mstats.pairs, (unsigned long long)mstats.inliers,
           (unsigned long long)mstats.putative);
    printf("  map     : %6.2f s  %u/%zu images registered, %zu points, %u camera(s)\n", t_map,
           reg, est.images, rec.points3D.size(), (uint32_t)rec.cameras.size());
    if (mst.rounds)
        printf("  manage  : %6.2f s  %zu -> %zu models, %zu merged, %zu refused, %zu grown, "
               "%zu repaired, %zu dropped, %zu -> %zu covered\n",
               t_manage, mst.models_before, mst.models_after, mst.merges, mst.merges_refused,
               mst.grown_images, mst.audited_repaired, mst.audited_out, mst.covered_before,
               mst.covered_after);
    printf("  total   : %6.2f s\n", t_extract + t_match + t_map + t_manage);
    printf("  model   : %.3f px mean reprojection (%.3f median) over %zu observations\n", mean,
           median, nobs);
    if (models.size() > 1) {
        // A fragmented capture: sparse/0 is the largest component, the rest are
        // separate reconstructions with no known transform between them (D41).
        printf("  models  : %zu components covering %zu/%zu distinct images\n", models.size(),
               distinctRegistered(models), est.images);
        printExtraModels(models, feats, "          ");
    }
    if (models.size() > 1)
        printf("  written : %s/{0..%zu}\n", sparsedir.string().c_str(), models.size() - 1);
    else
        printf("  written : %s\n", (sparsedir / "0").string().c_str());

    // Verdict, so a batch run can be scanned without reading every number.
    // Thresholds are deliberately loose -- this flags "obviously broken", not
    // "not as good as COLMAP".
    const double frac = est.images ? (double)reg / est.images : 0.0;
    if (reg < 2 || rec.points3D.empty()) {
        printf("  RESULT  : FAILED (no reconstruction)\n");
        return 2;
    }
    if (frac < 0.5 || mean > 2.0) {
        printf("  RESULT  : PARTIAL (%.0f%% of images registered, %.2f px)\n", 100 * frac, mean);
        return 3;
    }
    printf("  RESULT  : OK (%.0f%% of images registered, %.2f px)\n", 100 * frac, mean);
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: ssplat-sfm <command> [args]\n"
                        "  auto     image directory -> COLMAP sparse model, one command\n"
                        "  extract  detect SIFT features from an image or directory\n"
                        "  match    brute-force match a directory of feature files\n"
                        "  map      incremental reconstruction -> COLMAP sparse model\n"
                        "  merge    fuse the models of a fragmented capture into fewer\n"
                        "  ba       bundle-adjust a BAL problem (solver benchmark)\n"
                        "\n");
        return 1;
    }
    // Accept `--flag=value` as well as `--flag value`, everywhere. Only tokens
    // that start with `--` are split, and only at their first `=`, so a value
    // that itself contains one (`--camera-model cam0=opencv-fisheye`) is
    // untouched. Every subcommand parses positionally, so this is the one place
    // that has to know.
    std::vector<std::string> store;
    std::vector<char*> args;
    store.reserve((size_t)argc * 2);
    args.reserve((size_t)argc * 2);
    for (int i = 0; i < argc; i++) {
        std::string a = argv[i];
        size_t eq = a.find('=');
        if (i >= 2 && a.size() > 2 && a[0] == '-' && a[1] == '-' && eq != std::string::npos &&
            eq > 2) {
            store.push_back(a.substr(0, eq));
            store.push_back(a.substr(eq + 1));
        } else {
            store.push_back(std::move(a));
        }
    }
    for (std::string& a : store) args.push_back(&a[0]);
    argc = (int)args.size();
    argv = args.data();

    std::string cmd = argv[1];
    if (cmd == "auto") return cmdAuto(argc - 2, argv + 2);
    if (cmd == "extract") return cmdExtract(argc - 2, argv + 2);
    if (cmd == "match") return cmdMatch(argc - 2, argv + 2);
    if (cmd == "map") return cmdMap(argc - 2, argv + 2);
    if (cmd == "merge") return cmdMerge(argc - 2, argv + 2);
    if (cmd == "ba") return cmdBa(argc - 2, argv + 2);
    fprintf(stderr, "unknown command %s\n", cmd.c_str());
    return 1;
}
