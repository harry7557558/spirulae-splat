// Incremental Structure-from-Motion mapper (src/sfm/README.md).
//
// MVP control flow, ported in spirit from COLMAP's IncrementalMapper:
//   seed pair -> initialize (relative pose + triangulate)
//   loop: register-next (2D-3D via the correspondence graph, PnP RANSAC)
//         -> continue tracks + triangulate new points
//         -> periodic global BA + reprojection filtering
//
// Robustness (D36) follows COLMAP: strict seed acceptance with stepwise
// relaxation, PnP inlier-count and inlier-ratio gates with nonlinear pose
// refinement, iterated global BA + filtering (reprojection and triangulation
// angle) at 10% model growth, and de-registration of images the filtering
// hollows out -- the mapper's way of undoing a registration that stopped
// agreeing with the model.
//
// Remaining simplifications (src/sfm/README.md / D10): transitive tracks
// approximated by one-hop correspondences, global BA only (no local BA), no
// retriangulation/merge pass. Points are colored by averaging the
// per-keypoint colors sampled at extraction; each has a clear upgrade path.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "sfm/core/Features.h"
#include "sfm/core/Model.h"
#include "sfm/geometry/AbsolutePose.h"
#include "sfm/geometry/Triangulation.h"
#include "sfm/geometry/TwoView.h"
#include "sfm/map/Bundle.h"
#include "sfm/map/CorrespondenceGraph.h"
#include "sfm/map/Profile.h"

namespace sfm {

struct MapperOptions {
    double focal = 0;                  // 0 = COLMAP default (1.2 * max dim)
    // Every pixel threshold below is in *extraction* pixels -- the frame the
    // keypoints were measured in, not the frame they are stored in. They are
    // the same thing unless the extractor downscaled; Camera::pixel_scale is
    // the conversion, and it keeps `--max-error 4` meaning one thing across
    // quality presets and across a mixed-resolution capture (D47).
    double max_reproj_error = 3.0;     // extraction px (see above; D47)
    double min_tri_angle_deg = 1.5;
    // Seed acceptance (D36): thresholds are COLMAP's IncrementalMapper
    // defaults. If no candidate pair passes, initialize() relaxes them
    // stepwise rather than failing outright.
    double init_min_tri_angle_deg = 16.0;   // median angle over the seed points
    int init_min_inliers = 100;
    double init_max_forward_motion = 0.95;  // |baseline . viewing dir| cap
    // Registration gates: the ratio is COLMAP's abs_pose_min_inlier_ratio; it
    // rejects an image whose hundreds of 2D-3D correspondences agree only by
    // accident. The inlier count stays at 15 (COLMAP uses 30): sparse-match
    // sets live on 15-25 inlier registrations, and the transactional
    // refinement now catches the ones that turn out toxic (D36).
    int min_num_pnp_inliers = 15;
    double min_pnp_inlier_ratio = 0.25;
    int max_reg_trials = 3;            // per image, COLMAP's default
    // Focal search when a camera group registers its first image and its focal
    // is still the 1.2*max_dim guess (D18). COLMAP solves this with a P4Pf
    // minimal solver; a log-spaced sweep of P3P hypotheses is the cheaper
    // equivalent, and the ratio bounds are COLMAP's.
    int min_model_size = 10;           // COLMAP's default; smaller -> retry the seed
    // Stop retrying seeds once a model is both >= min_model_size and covers
    // this fraction of the images. Below that the model is "a" reconstruction
    // but not "the" reconstruction, and another seed is usually worth the time.
    double min_model_fraction = 0.5;
    int max_init_trials = 8;           // seed attempts before keeping the best
    // Multiple models (D41). A capture that does not form one connected view
    // graph -- two ends of a building with nothing joining them, a sequence
    // broken by a featureless corridor -- yields several reconstructions, and
    // there is no correct way to fuse them without knowing the transform
    // between them. COLMAP's answer is to emit them all as sparse/0, sparse/1,
    // ... for a later merge step, and this matches it:
    //   * the primary model is still the best of max_init_trials seed attempts
    //     (D19) -- unchanged, so a scene that reconstructs as one model behaves
    //     exactly as before;
    //   * afterwards, seeds are searched among images no kept model registered
    //     (COLMAP FindFirstInitialImage's num_registrations == 0 rule) and each
    //     resulting model is kept if it reaches min_model_size;
    //   * a sub-model may re-register images an earlier model already holds --
    //     that overlap is what a later merge aligns on -- but stops growing
    //     once it has taken max_model_overlap of them (COLMAP's rule).
    // max_num_models 1 restores single-model output.
    int max_num_models = 50;           // COLMAP's default
    int max_model_overlap = 20;        // COLMAP's default
    // Total seed attempts the further-models search may spend, counting the
    // ones discarded for being under min_model_size (COLMAP's init_num_trials,
    // same role and same default). Without it a dataset whose leftovers are
    // dust keeps seeding two-image models until the candidate pair list runs
    // out, which on a 2000-image scene is tens of thousands of attempts.
    int max_model_trials = 20;  // 200
    int focal_search_samples = 15;     // 0 disables the search
    double min_focal_ratio = 0.1, max_focal_ratio = 10.0;
    // Focal bootstrap (D48), for the captures where bundle adjustment cannot
    // recover the focal from the 1.2*max_dim guess: rotation-degenerate motion
    // (a dashcam, a dolly, a rail scan) and genuinely wide lenses. See
    // bootstrapFocalLength. 0 = off; the descent usually stops after 1-2.
    int focal_trials = 5;
    // Relative gain in mutually consistent observations that justifies moving
    // the focal. The guess is at least a *consistent* place to start and a young
    // model's observation count is a noisy statistic, so a hypothesis a few
    // percent ahead of it has shown nothing. Measured separations are far
    // larger on a few real-world datasets.
    double focal_min_gain = 0.15;
    size_t focal_model_size = 20;      // images per trial model
    // Reported, not a gate: the orientation spread below which the focal is a
    // free parameter of the reconstruction, so a bad one leaves no trace in the
    // residual. A straight KITTI drive spans 8.5 deg over 114 frames; the
    // tightest ordinary capture measured (mip-NeRF 360 garden's first 20 images,
    // one corner of an object orbit) spans 32.
    double focal_max_rot_spread_deg = 15.0;
    // Global-refinement cadence and image filtering (D36; D38). COLMAP's
    // trigger ratio, kept: D37's 1.25 measured fine on easy sets but cost
    // registration + pose accuracy on a real-world dataset -- fewer refine rounds
    // means fewer retriangulation passes, exactly the machinery D36 added for
    // such sets. Speed comes from convergence-adaptive iterations + the
    // persistent solver instead (D38).
    double ba_growth_ratio = 1.1;      // COLMAP ba_global_images_ratio
    int ba_max_refinements = 5;        // final pass; growth passes use 2
    double ba_refine_change = 0.0005;  // stop when changed-obs fraction is below
    // Growth-phase BAs stop when relative cost improvement stays below
    // ba_growth_rtol for ba_growth_patience accepted steps (D38): iteration
    // count adapts to convergence -- a shaky young model keeps iterating, a
    // converged large one stops after a few -- replacing D37's fixed
    // iteration cap, which starved exactly the models that still needed work.
    // Final passes keep the solver's tight defaults. 0 = solver default.
    double ba_growth_rtol = 1e-4;
    int ba_growth_patience = 5;
    // Reprojection acceptance for retriangulation/track completion, as a
    // fraction of max_reproj_error (the churn hysteresis); 0 disables the
    // retriangulation pass entirely.
    double retri_scale = 0.75;
    // Auditing an assembled model (D44). An image is put back only when the
    // structure it did *not* bring supports a competing pose: one that clears
    // the registration gates, explains `audit_alternative_factor` times as
    // many of those correspondences as the current pose does, and points
    // somewhere else. The evidence floor keeps the test off images the model
    // barely sees, where a lucky RANSAC on 20 stray correspondences would
    // otherwise unseat a perfectly good pose.
    int audit_min_evidence = 40;
    int audit_min_alternative = 25;
    double audit_alternative_factor = 3.0;
    double audit_min_rotation_deg = 5.0;
    // ... or a camera-center shift of this fraction of the model's own scale
    // (the RMS spread of its camera centers), which is what catches an image
    // that kept its orientation and moved.
    double audit_min_shift_frac = 0.01;
    int audit_ransac_trials = 1000;
    int min_image_points = 5;          // de-register images that fall below this
    double max_extra_param = 1.0;      // |distortion param| beyond this = bogus
    CamModel camera_model = CamModel::Radial;  // distortion model for new cameras (D29)
    // Starting intrinsics per camera id, built by sfm/core/CameraSetup.h from
    // --camera-model / --focal / EXIF (D46). A camera with no entry here falls
    // back to Camera::defaultFor(focal, camera_model), the old behaviour, so
    // the self-tests and any library caller are unaffected.
    std::map<uint32_t, Camera> initial_cameras;
    // Cameras whose focal came from a *per-group* prior (EXIF, or a --focal
    // that named the group) rather than a guess. The registration focal sweep
    // leaves these alone; a dataset-wide --focal is not in here, deliberately,
    // because it says nothing about which group it describes (D45/D46).
    std::set<uint32_t> known_focal_cameras;
    // Cameras whose focal was *supplied* -- a per-group prior, EXIF, or a
    // dataset-wide --focal (CameraSetup::focal_given). Superset of the above.
    // The focal bootstrap only searches cameras outside it: a supplied focal may
    // be worth refining, but it is never worth replacing with a guess (D48).
    std::set<uint32_t> given_focal_cameras;
    std::string ba_loss = "huber";     // robust loss for mapping-time BA (D36)
    double ba_loss_param = 2.0;        // Huber delta / Cauchy c, in extraction px
    // Let BA move each camera's principal point off the image centre. Off, as
    // in COLMAP: it is nearly the same parameter as a camera rotation, so on a
    // rig its per-group drift lands in the relative orientation of the lenses
    // (D50). See BundleOptions::refine_principal_point.
    bool refine_principal_point = false;
    // Release it for one global bundle adjustment at the very end, when the
    // model is complete -- COLMAP's documented advice, and a different question
    // from refining it *during* reconstruction (D51). `Mapper::polish` is what
    // runs it; the CLI drives that, not the mapper's own loop.
    size_t pp_min_images = 20;   // ... for groups with at least this many images
    RealCfg ba_real = RealCfg::F64;
    int device = -1;
    bool verbose = true;
};

class Mapper {
public:
    // `camera_ids[i]` is the (1-based) camera image i belongs to; images sharing
    // an id share intrinsics and one BA parameter group. Empty = one shared
    // camera for everything, the old behaviour. The CLI builds this from
    // --camera-mode (D17).
    Mapper(const MatchesDatabase& db, const std::vector<FeatureSet>& feats, MapperOptions opt,
           std::vector<uint32_t> camera_ids = {})
        : db_(db), feats_(feats), opt_(opt), cam_ids_(std::move(camera_ids)) {}

    // All reconstructions the dataset supports, largest first (by 3D point
    // count -- COLMAP's ReconstructionManager::Write ordering, so models[0] is
    // what lands in sparse/0). Never empty in practice: a dataset that seeds at
    // all yields at least one model, and one that does not yields a single
    // empty reconstruction so callers have something to report on.
    std::vector<Reconstruction> run() {
        auto prof_start = std::chrono::steady_clock::now();
        ensureSetup();
        std::vector<Reconstruction> models;

        // ---- phase 1: the primary model, best of max_init_trials seeds ----
        // A reconstruction is only as good as its seed: the densest verified
        // pair can sit in a small corner of the view graph, and then the model
        // stops after a handful of images with no way to notice. COLMAP's
        // answer is to discard a model smaller than `min_model_size` and start
        // again from the next candidate seed; we keep the largest attempt (D19)
        // -- and, since D41, *all* the attempts, because on a fragmented
        // capture the rejected ones are exactly the other components. They are
        // admitted below rather than rebuilt from scratch.
        // Before anything is built on it: is the focal something this capture
        // can actually determine? (D48; a no-op unless it cannot.)
        {
            ProfTimer pt(g_map_prof.init_seed);
            bootstrapFocalLength();
        }

        std::vector<Reconstruction> attempts;
        size_t seed_from = 0;
        for (int attempt = 0; attempt < std::max(1, opt_.max_init_trials); attempt++) {
            resetModel();
            bool seeded;
            {
                ProfTimer pt(g_map_prof.init_seed);
                seeded = initialize(seed_from);
            }
            if (!seeded) {
                if (opt_.verbose && attempt == 0) reportInitFailure();
                break;
            }
            globalRefine(false);  // COLMAP's two-view BA before any growth
            grow();
            uint32_t reg = rec_.numRegistered();
            if (reg) attempts.push_back(snapshotModel());
            uint32_t enough = (uint32_t)std::max(
                (double)opt_.min_model_size, opt_.min_model_fraction * db_.images.size());
            if (reg >= enough) break;
            if (opt_.verbose)
                fprintf(stderr, "[map] model too small (%u < %u), retrying from another seed\n",
                        reg, enough);
        }
        if (attempts.empty()) {
            // Nothing seeded at all. Hand back one empty reconstruction so the
            // caller has a model to report on, as before D41.
            resetModel();
            models.push_back(snapshotModel());
            return finishRun(models, prof_start);
        }

        // Largest attempt first, so the primary model is the one D19 would have
        // returned; each further attempt is admitted only if it is genuinely a
        // *different* component. Seeds land in the same component all the time,
        // and writing three views of one component as three models would be worse
        // than useless -- max_model_overlap is the same bound COLMAP uses to
        // decide when a sub-model has drifted into an existing one.
        std::stable_sort(attempts.begin(), attempts.end(),
                         [](const Reconstruction& a, const Reconstruction& b) {
                             return a.numRegistered() > b.numRegistered();
                         });
        for (size_t i = 0; i < attempts.size(); i++) {
            if ((int)models.size() >= std::max(1, opt_.max_num_models)) break;
            std::string why;
            if (i > 0 && !admitModel(attempts[i], why)) {
                if (opt_.verbose)
                    fprintf(stderr, "[map] seed attempt discarded: %s\n", why.c_str());
                continue;
            }
            claimImages(attempts[i]);
            recordCameras(attempts[i]);
            models.push_back(std::move(attempts[i]));
        }

        // ---- phase 2: further models from whatever is still unclaimed ----
        // Seeds come only from images no kept model registered, so this does
        // nothing at all when the primary model covers the dataset -- a
        // single-component capture takes the pre-D41 path exactly.
        seedFurtherModels(models);
        return finishRun(models, prof_start);
    }

    // ---- the engine, driven from outside (D44) ----------------------------
    //
    // run() is one policy over these; ModelManager (sfm/map/Manager.h) is
    // another, and a bottom-up hierarchical mapper would be a third. They all
    // need the same three operations on an *existing* model, which is why they
    // are public: keep growing it, refine it, or look for more models beside
    // it. Every one of them starts by adopting the model into `rec_`, so the
    // mapper never has to have built it itself -- it may equally have come off
    // disk or out of a merge.

    struct GrowStats {
        uint32_t before = 0, after = 0;
        uint32_t registered = 0;   // images this pass brought in
        bool refined = false;      // whether a final refinement ran
    };

    // Adopt `m` and keep registering into it until nothing else fits. Claims
    // are cleared for the pass: an image another model holds is a legitimate
    // target here, and the overlap it creates is what lets the two models
    // merge afterwards (D43).
    //
    // A pass that registers nothing returns `m` untouched -- not a re-refined
    // copy of it. That is what makes the manager's grow round free on a
    // dataset it cannot help: no BA runs, no numbers move.
    Reconstruction continueFrom(const Reconstruction& m, GrowStats* out = nullptr) {
        ensureSetup();
        resetModel();
        adopt(m);
        model_count_.clear();
        rebuildScores();
        GrowStats st;
        st.before = rec_.numRegistered();
        st.registered = growLoop();
        st.after = rec_.numRegistered();
        if (out) *out = st;
        if (!st.registered) return m;
        st.refined = true;
        checkedRefine(true);
        if (out) {
            out->refined = true;
            out->after = rec_.numRegistered();
        }
        return snapshotModel();
    }

    // Adopt and bundle-adjust, with the mapper's own filtering and
    // de-registration rules. This is what a merged model needs: two halves
    // glued along a seam that has never been optimized as one.
    Reconstruction refine(const Reconstruction& m) {
        ensureSetup();
        resetModel();
        adopt(m);
        rebuildScores();
        globalRefine(true);
        return snapshotModel();
    }

    // The same, with the principal point released -- for use once, on a model
    // that is finished (D51). Held during reconstruction the principal point is
    // an ill-posed parameter that trades against camera rotation (D50); on a
    // complete model with many images sharing one set of intrinsics the trade
    // is pinned down, which is exactly the distinction COLMAP's documentation
    // draws. Groups with fewer than `pp_min_images` images keep theirs: a group
    // of one image has nothing to share with, so moving its principal point is
    // just rotating that camera.
    Reconstruction polish(const Reconstruction& m) {
        // Only a single-camera-group model. With two or more groups each one's
        // principal point drifts its own way and the difference is a real error
        // in their relative orientation -- on a dual-fisheye rig, whose
        // inter-lens rotation is a physical constant, that cost 21 points of
        // AUC (D51). One group has no such difference to get wrong.
        std::set<uint32_t> groups;
        for (const auto& kv : m.images)
            if (kv.second.registered) groups.insert(kv.second.camera_id);
        if (groups.size() != 1) {
            if (opt_.verbose)
                fprintf(stderr,
                        "[map] %zu camera groups: skipping the final principal-point pass "
                        "(it would move them apart, D51)\n", groups.size());
            return m;
        }
        ensureSetup();
        resetModel();
        adopt(m);
        rebuildScores();
        std::map<uint32_t, Vec2> before;
        for (const auto& kv : rec_.cameras) before[kv.first] = {kv.second.cx, kv.second.cy};
        pp_free_ = true;
        globalRefine(true);
        pp_free_ = false;
        if (opt_.verbose)
            for (const auto& kv : rec_.cameras) {
                const Vec2& b = before[kv.first];
                const double d = std::hypot(kv.second.cx - b.x, kv.second.cy - b.y);
                if (d > 1e-9)
                    fprintf(stderr,
                            "[map] camera %u principal point %.1f,%.1f -> %.1f,%.1f (%.1f px)\n",
                            kv.first, b.x, b.y, kv.second.cx, kv.second.cy, d);
            }
        return snapshotModel();
    }

    struct AuditStats {
        uint32_t checked = 0, deregistered = 0, unsupported = 0, reregistered = 0;
    };

    // Ask the correspondence graph whether every image in `m` belongs where
    // `m` says it does, and move the ones that do not (D44).
    //
    // A model assembled from parts -- merged (D43), or stacked up by a
    // hierarchical mapper -- can place an image somewhere its own observations
    // still support, because those observations came along with it. What it
    // cannot fake is the rest of the model: if the structure the image did
    // *not* bring supports a different pose decisively better (see
    // poseContradicted), the image is in the wrong place, and no amount of
    // bundle adjustment will walk it back.
    //
    // A contradicted image is moved to the pose that won and re-triangulated
    // there, rather than deleted. Deleting was tried first and cost 60 images
    // on one a merge on a dataset: re-registration applies the full
    // registration gates, including an inlier *ratio* that a misplaced image's
    // huge junk-dominated correspondence pool cannot meet, so the ones the
    // audit unseated could not come back. The evidence that selected the new
    // pose is weaker than registration demands but decisively stronger than
    // what the old pose had, and the refinement that follows -- BA, filtering,
    // de-registration of anything hollow -- judges the result under the
    // ordinary rules.
    //
    // Measured on the same dataset: without this, merging models that have
    // each drifted along a long walk left ~7% of the rig frames (which share a
    // pose by construction, tools/rig_check.py) tens of degrees apart.
    Reconstruction audit(const Reconstruction& m, AuditStats* out = nullptr) {
        ensureSetup();
        resetModel();
        adopt(m);
        // As in continueFrom: this is not a sub-model being built beside the
        // others, so the claim bookkeeping (and the overlap break it drives)
        // must not stop the re-registration loop below.
        model_count_.clear();
        rebuildScores();
        AuditStats st;
        std::vector<std::pair<uint32_t, Pose>> repairs;
        for (const auto& kv : rec_.images) {
            if (!kv.second.registered) continue;
            st.checked++;
            Pose alt;
            if (poseContradicted(kv.first, alt)) repairs.emplace_back(kv.first, alt);
        }
        st.unsupported = (uint32_t)repairs.size();
        if (!repairs.empty()) {
            if (opt_.verbose)
                fprintf(stderr, "[map] audit: %u/%u image(s) sit where the rest of the model "
                        "contradicts them; moving\n", st.unsupported, st.checked);
            // Detach first, all of them: their old observations are evidence
            // for the old pose and must not survive it.
            for (const auto& r : repairs) deregisterImage(r.first);
            for (const auto& r : repairs) {
                Image& im = rec_.images[r.first];
                im.pose = r.second;
                im.registered = true;
            }
            rebuildScores();
            // Attach to what the new pose can see, then triangulate what
            // nothing sees yet -- registration's two steps, for a pose that
            // came from outside instead of from PnP. Without the first, a
            // repaired image owns no observations at all (everything it looks
            // at is already triangulated, so there is nothing to *create*) and
            // the next filtering pass de-registers it as hollow.
            for (const auto& r : repairs) attachExisting(r.first);
            for (const auto& r : repairs) triangulateForImage(r.first);
            rebuildScores();
            // Images that could not register before may be able to now, both
            // because the model changed and because their trial budget is
            // reset; that is a bonus, not a side effect to design around.
            reg_trials_.assign(db_.images.size(), 0);
            growLoop();
        }
        checkedRefine(true);
        for (const auto& r : repairs)
            if (rec_.images.at(r.first).registered) st.reregistered++;
        st.deregistered = st.unsupported - st.reregistered;
        if (out) *out = st;
        return snapshotModel();
    }

    // Seed and grow further models among images that `models` does not cover,
    // appending each admitted one. `restart_relaxation` re-arms the seed
    // threshold ladder, which a later pass needs: the first pass leaves it
    // exhausted, but by then the claimed set has changed and pairs that were
    // ineligible are not any more.
    void seedFurtherModels(std::vector<Reconstruction>& models, bool restart_relaxation = false) {
        ensureSetup();
        if (restart_relaxation) init_relax_ = 0;
        int trials = 0;
        while ((int)models.size() < std::max(1, opt_.max_num_models) && unclaimedImages() > 0) {
            if (++trials > std::max(1, opt_.max_model_trials)) {
                if (opt_.verbose)
                    fprintf(stderr, "[map] sub-model search hit its %d-attempt budget with %zu "
                            "image(s) unaccounted for\n", opt_.max_model_trials,
                            unclaimedImages());
                break;
            }
            resetModel();
            size_t from = 0;
            bool seeded;
            {
                ProfTimer pt(g_map_prof.init_seed);
                seeded = initialize(from);
            }
            if (!seeded) break;  // nothing left that can seed a model
            globalRefine(false);
            grow();
            Reconstruction sub = snapshotModel();
            std::string why;
            if (!admitModel(sub, why)) {
                // Not worth a directory of its own, and claiming its images
                // would only starve later passes. Leave them for the next seed;
                // `used_seeds_` stops this pair being tried again.
                if (opt_.verbose)
                    fprintf(stderr, "[map] sub-model discarded: %s\n", why.c_str());
                continue;
            }
            const uint32_t reg = sub.numRegistered();
            models.push_back(std::move(sub));
            claimImages(models.back());
            recordCameras(models.back());
            if (opt_.verbose)
                fprintf(stderr, "[map] sub-model %zu: %u images, %zu points (%zu images left)\n",
                        models.size() - 1, reg, models.back().points3D.size(), unclaimedImages());
        }
    }

    // Record `models` as the claimed set, replacing whatever was there.
    void claimAll(const std::vector<Reconstruction>& models) {
        model_count_.assign(db_.images.size(), 0);
        for (const Reconstruction& m : models) claimImages(m);
    }

    // ---- does a model agree with the two-view geometries it was built from? --
    //
    // Every acceptance test the merger can run by itself is computed from the
    // same evidence the alignment used: the images the two models share, and
    // the points they both triangulated. A capture that walks past the same
    // facade twice can satisfy all of it and still be glued together wrongly,
    // because the repeated structure is genuinely consistent -- locally.
    //
    // The correspondence graph knows something the models do not. Two images
    // were matched and verified long before any model existed, and their
    // two-view geometry does not care where a merge, or a registration, later
    // put them. So: compute the relative pose the model implies for a verified
    // pair and ask how much of that pair's own evidence it still explains. A
    // correct model reproduces those geometries; a wrong one cannot, and no
    // amount of internal self-consistency will save it.
    //
    // Bearings throughout, so this is meaningful for a fisheye (D45).

    // Fraction of a verified pair's matches the model's relative pose explains,
    // or -1 when the pair cannot be judged (an image missing, no camera).
    double pairAgreement(const Reconstruction& m, const TwoViewMatches& p, double max_error_px,
                         double model_scale) const {
        auto ia = m.images.find(p.image1);
        auto ib = m.images.find(p.image2);
        if (ia == m.images.end() || ib == m.images.end()) return -1;
        if (!ia->second.registered || !ib->second.registered || p.matches.empty()) return -1;
        auto ca = m.cameras.find(ia->second.camera_id);
        auto cb = m.cameras.find(ib->second.camera_id);
        if (ca == m.cameras.end() || cb == m.cameras.end()) return -1;
        Mat3 R = mul(ib->second.pose.R, transpose(ia->second.pose.R));
        Vec3 t = ib->second.pose.t - mul(R, ia->second.pose.t);
        // px -> rad, in the frame the keypoints were measured in (D47).
        const double thr =
            0.5 * (ca->second.errRad(max_error_px) + cb->second.errRad(max_error_px));
        const bool wide_baseline = model_scale > 0 && t.norm() > 1e-3 * model_scale;
        Mat3 E = mul(crossMatrix(t.normalized()), R);
        size_t ok = 0;
        for (const FeatureMatch& fm : p.matches) {
            Vec3 b1 = ca->second.bearing(kp(p.image1, fm.idx1));
            Vec3 b2 = cb->second.bearing(kp(p.image2, fm.idx2));
            double err2;
            if (wide_baseline) {
                err2 = sampsonSqBearing(E, b1, b2);
            } else {
                // No baseline to speak of: the epipolar constraint is vacuous,
                // so compare the rays directly through R.
                Vec3 pred = mul(R, b1);
                Vec3 c = pred.cross(b2);
                double ang = std::atan2(c.norm(), pred.dot(b2));
                err2 = ang * ang;
            }
            if (err2 < thr * thr) ok++;
        }
        return (double)ok / (double)p.matches.size();
    }

    struct SeamCheck {
        size_t cross_pairs = 0;    // verified pairs found across the seam
        size_t tested = 0;         // ... of which were actually evaluated
        size_t agree = 0;
        double median_frac = 0;    // median per-pair fraction of matches explained
    };

    SeamCheck checkSeam(const Reconstruction& m, const std::set<uint32_t>& src_side,
                        double max_error_px = 8.0, double min_pair_frac = 0.5,
                        size_t max_pairs = 600) const {
        SeamCheck sc;
        std::vector<const TwoViewMatches*> cross;
        for (const TwoViewMatches& p : db_.pairs) {
            auto ia = m.images.find(p.image1);
            auto ib = m.images.find(p.image2);
            if (ia == m.images.end() || ib == m.images.end()) continue;
            if (!ia->second.registered || !ib->second.registered) continue;
            if (src_side.count(p.image1) == src_side.count(p.image2)) continue;  // same side
            cross.push_back(&p);
        }
        sc.cross_pairs = cross.size();
        if (cross.empty()) return sc;
        // Prefer the pairs with the most evidence, then spread the sample over
        // them: a seam is only as good as its strongest links, and a few
        // hundred are plenty to tell a good merge from a bad one.
        std::stable_sort(cross.begin(), cross.end(), [](const TwoViewMatches* a,
                                                        const TwoViewMatches* b) {
            return a->matches.size() > b->matches.size();
        });
        if (cross.size() > max_pairs) cross.resize(max_pairs);

        const double scale = modelScaleOf(m);
        std::vector<double> fracs;
        for (const TwoViewMatches* p : cross) {
            double frac = pairAgreement(m, *p, max_error_px, scale);
            if (frac < 0) continue;
            fracs.push_back(frac);
            sc.tested++;
            if (frac >= min_pair_frac) sc.agree++;
        }
        if (!fracs.empty()) {
            std::sort(fracs.begin(), fracs.end());
            sc.median_frac = fracs[fracs.size() / 2];
        }
        return sc;
    }

    // ---- split a model along the geometries it violates (D45) --------------
    //
    // The same measurement, applied to every verified pair inside one model,
    // says more than pass/fail: it says *where* the model stops being true.
    // Keep only the pairs the model reproduces, and the images fall into
    // connected groups. One group means the model is coherent. Two large ones
    // joined by nothing means two pieces of the capture were welded at the
    // wrong relative pose -- by a merge, or by a chain of registrations through
    // repeated structure -- and no bundle adjustment will ever pull them apart,
    // because each piece is internally perfect.
    //
    // Splitting is the honest response: the pieces are real reconstructions,
    // and once separated they can be re-merged (this time against the seam
    // test) or grown independently. Images in groups too small to keep are
    // simply de-registered; growth will offer them a place again.
    struct SplitStats {
        size_t pairs_tested = 0, pairs_agree = 0;
        size_t groups = 0;          // connected groups of agreeing images
        size_t largest = 0;
        size_t dropped_images = 0;  // in groups too small to keep
        // Per-pair fraction explained, sorted. A model with a real problem has
        // a bimodal distribution -- most pairs near 1, a tail near 0 -- while a
        // threshold that is merely too tight shows a smooth spread.
        std::vector<double> fractions;
        double percentile(double q) const {
            if (fractions.empty()) return 0;
            size_t i = (size_t)(q * (double)(fractions.size() - 1));
            return fractions[i];
        }
    };

    std::vector<Reconstruction> splitInconsistent(const Reconstruction& m, double max_error_px,
                                                  double min_pair_frac, int min_matches,
                                                  size_t min_group, SplitStats* out = nullptr) const {
        SplitStats st;
        std::vector<uint32_t> ids;
        std::map<uint32_t, size_t> pos;
        for (const auto& kv : m.images)
            if (kv.second.registered) { pos[kv.first] = ids.size(); ids.push_back(kv.first); }
        if (ids.size() < 2) return {m};

        std::vector<size_t> parent(ids.size());
        for (size_t i = 0; i < parent.size(); i++) parent[i] = i;
        std::function<size_t(size_t)> find = [&](size_t x) {
            while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
            return x;
        };
        const double scale = modelScaleOf(m);
        for (const TwoViewMatches& p : db_.pairs) {
            if ((int)p.matches.size() < min_matches) continue;
            auto a = pos.find(p.image1);
            auto b = pos.find(p.image2);
            if (a == pos.end() || b == pos.end()) continue;
            double frac = pairAgreement(m, p, max_error_px, scale);
            if (frac < 0) continue;
            st.pairs_tested++;
            st.fractions.push_back(frac);
            if (frac < min_pair_frac) continue;
            st.pairs_agree++;
            size_t ra = find(a->second), rb = find(b->second);
            if (ra != rb) parent[ra] = rb;
        }

        std::sort(st.fractions.begin(), st.fractions.end());
        std::map<size_t, std::vector<uint32_t>> groups;
        for (size_t i = 0; i < ids.size(); i++) groups[find(i)].push_back(ids[i]);
        std::vector<std::vector<uint32_t>> gs;
        for (auto& kv : groups) gs.push_back(std::move(kv.second));
        std::sort(gs.begin(), gs.end(),
                  [](const std::vector<uint32_t>& a, const std::vector<uint32_t>& b) {
                      return a.size() > b.size();
                  });
        st.groups = gs.size();
        st.largest = gs.empty() ? 0 : gs[0].size();
        if (out) *out = st;
        if (gs.size() <= 1) return {m};

        std::vector<Reconstruction> parts;
        for (const std::vector<uint32_t>& g : gs) {
            if (g.size() < min_group) { st.dropped_images += g.size(); continue; }
            std::set<uint32_t> keep(g.begin(), g.end());
            parts.push_back(subsetModel(m, keep));
        }
        if (out) *out = st;
        if (parts.empty()) return {m};
        return parts;
    }

    // "Were these two images matched to each other with real support?" -- the
    // question findDuplicateStructure needs the correspondence graph for
    // (D45). Built once per call over the verified pair list.
    MatchedFn matchedPredicate(int min_matches = 15) const {
        auto index = std::make_shared<std::set<std::pair<uint32_t, uint32_t>>>();
        for (const TwoViewMatches& p : db_.pairs) {
            if ((int)p.matches.size() < min_matches) continue;
            uint32_t a = p.image1, b = p.image2;
            if (a > b) std::swap(a, b);
            index->insert({a, b});
        }
        return [index](uint32_t a, uint32_t b) {
            if (a > b) std::swap(a, b);
            return index->count({a, b}) > 0;
        };
    }

    // Bundle-adjust every component at once with one shared set of intrinsics
    // per camera group (D45; runJointBA). Cheaper than it sounds -- it is one
    // solve instead of N -- and it is the only place where a small component's
    // intrinsics are constrained by the big component's evidence. Each model is
    // then filtered and re-refined through the ordinary path, so the same
    // observation and image gates apply as after any other BA.
    void jointRefine(std::vector<Reconstruction>& models) {
        ensureSetup();
        size_t live = 0;
        for (const Reconstruction& m : models)
            if (m.numRegistered() >= 2) live++;
        if (live < 2) return;
        BundleOptions bo;
        bo.real = opt_.ba_real;
        bo.device = opt_.device;
        bo.verbose = false;
        bo.loss = opt_.ba_loss;
        bo.loss_param = (float)(opt_.ba_loss_param * medianPixelScale());
        bo.refine_principal_point = opt_.refine_principal_point || pp_free_;
        bo.pp_min_images = opt_.pp_min_images;
        bo.shared_ctx = &ba_ctx_;
        runJointBA(models, bo);
        // Deliberately no per-model refine afterwards: that would re-fit each
        // component's intrinsics to its own observations and undo the sharing
        // this pass exists for. Observations the shared solution no longer
        // explains are dropped by the next audit or growth pass, which refine
        // through the ordinary gates.
        clearCameraConsensus();
        for (const Reconstruction& m : models) recordCameras(m);
    }

    // ---- camera consensus across sub-models (D45) -------------------------
    //
    // Physically, one camera group is one lens: the same intrinsics whichever
    // component an image ended up in. The mapper used to forget that between
    // models -- every new sub-model started its cameras from the geometric
    // default and ran its own focal search, on a handful of images with almost
    // no parallax.
    //
    // So: a model that is admitted publishes its intrinsics, weighted by how
    // many images constrain them, and every later model starts from the best
    // set published so far instead of from the default. Weight, not recency,
    // decides -- the primary model is built first and is nearly always the
    // best-constrained, and a 12-image component never overwrites it.
    void recordCameras(const Reconstruction& m) {
        std::map<uint32_t, double> weight;
        for (const auto& kv : m.images)
            if (kv.second.registered) weight[kv.second.camera_id] += 1.0;
        for (const auto& kv : weight) {
            auto it = m.cameras.find(kv.first);
            if (it == m.cameras.end()) continue;
            auto cur = cam_consensus_.find(kv.first);
            if (cur == cam_consensus_.end() || kv.second > cur->second.second)
                cam_consensus_[kv.first] = {it->second, kv.second};
        }
    }

    // Drop what earlier models published (the manager re-publishes from the
    // current set, so a model that has since been split or repaired does not
    // keep voting).
    void clearCameraConsensus() { cam_consensus_.clear(); }

    // The intrinsics a fresh model would start from, for reporting.
    const std::map<uint32_t, std::pair<Camera, double>>& cameraConsensus() const {
        return cam_consensus_;
    }

    // Does the rest of the model support a *different* pose for this image
    // than the one it has?
    //
    // Asking whether the current pose is "supported" does not work, and the
    // measurement says why: only the features that carry no 3D point of their
    // own are evidence (an image that was moved wrongly kept its own tracks,
    // which reproject perfectly wherever it went), and that pool is mostly
    // junk -- one-hop graph approximations and matches an earlier filter
    // already rejected.
    //
    // Posing it as a competition does separate. Run the ordinary PnP RANSAC on
    // exactly that pool: if the outside structure has a pose for this image
    // that clears the registration gates and is somewhere else entirely, the
    // image is in the wrong place, and no amount of bundle adjustment will
    // walk it back. If the pool is noise, RANSAC finds nothing and the image
    // is left alone -- which is the common case and costs one failed RANSAC.
    //
    // On a hit, `alternative` is the pose that won, and the caller moves the
    // image there rather than throwing it away: this evidence is by
    // construction weaker than registration demands, but it is decisively
    // better than the pose in place, and the refinement that follows filters
    // the result honestly.
    bool poseContradicted(uint32_t img, Pose& alternative) const {
        std::vector<Vec3> X, br;
        const Image& im = rec_.images.at(img);
        for (uint32_t f = 0; f < feats_[img].count(); f++) {
            if (im.point3D_ids[f] != kInvalidPoint3D) continue;  // evidence it brought itself
            for (const Correspondence& c : graph_.at(img, f)) {
                if (c.image_id == img) continue;
                const Image& oi = rec_.images.at(c.image_id);
                if (!oi.registered) continue;
                uint64_t pid = oi.point3D_ids[c.feature_idx];
                if (pid == kInvalidPoint3D) continue;
                auto pt = rec_.points3D.find(pid);
                if (pt == rec_.points3D.end()) continue;
                X.push_back(pt->second.xyz);
                br.push_back(bearing(img, f));
                break;
            }
        }
        const int n = (int)X.size();
        if (n < opt_.audit_min_evidence) return false;  // nothing to contradict it with

        // How well the current pose explains that pool, for comparison.
        const double thr = camOf(img).errRad(opt_.max_reproj_error);
        const double thr2 = thr * thr;
        int cur = 0;
        for (int k = 0; k < n; k++) cur += pnpResidualSq(im.pose, X[k], br[k]) < thr2 ? 1 : 0;

        PnPResult r = ransacPnP(X, br, camOf(img).focal(), errPx(img), 0,
                                opt_.audit_ransac_trials);
        bool contradicted = false;
        double rot_deg = 0, shift = 0;
        // Note what is *not* required: a fraction of the pool. The pool of an
        // image that was placed wrongly is enormous precisely because none of
        // it was spliced, so a ratio gate hides the very case this exists for.
        // Absolute support, and dominance over the pose in place, are the honest tests.
        if (r.success && r.num_inliers >= opt_.audit_min_alternative &&
            r.num_inliers > (int)(opt_.audit_alternative_factor * cur)) {
            // A pose is only "different" if it is different: the alternative
            // usually *is* the current pose, recovered from the same geometry,
            // and finding it again is a confirmation rather than a problem.
            //
            // Both halves have to be tested. An image put down in the wrong
            // *place* -- a repeated facade, the same corridor one floor up --
            // keeps its orientation and only moves, which is precisely the
            // failure that survived a rotation-only test on a dataset.
            Mat3 D = mul(r.pose.R, transpose(im.pose.R));
            double tr = std::max(-1.0, std::min(1.0, (D[0] + D[4] + D[8] - 1) * 0.5));
            rot_deg = std::acos(tr) * 180.0 / M_PI;
            shift = (cameraCenter(r.pose) - cameraCenter(im.pose)).norm() / modelScale();
            contradicted =
                rot_deg > opt_.audit_min_rotation_deg || shift > opt_.audit_min_shift_frac;
            if (contradicted) alternative = r.pose;
        }
        if (audit_dump_)
            fprintf(stderr,
                    "[audit] %s: pool %d, current %d, alternative %d, rot %.1f deg, shift %.4f "
                    "-> %s\n", db_.images[img].name.c_str(), n, cur,
                    r.success ? r.num_inliers : 0, rot_deg, shift,
                    contradicted ? "CONTRADICTED" : "ok");
        return contradicted;
    }

    // A length to measure pose differences against, since a reconstruction has
    // no units: the RMS distance of the registered camera centers from their
    // centroid. Cached per adopted model -- the audit asks for it once per
    // image and the model does not move underneath it.
    double modelScale() const {
        if (scale_cache_ > 0) return scale_cache_;
        Vec3 c{0, 0, 0};
        size_t n = 0;
        for (const auto& kv : rec_.images)
            if (kv.second.registered) { c = c + cameraCenter(kv.second.pose); n++; }
        if (!n) return scale_cache_ = 1.0;
        c = c * (1.0 / (double)n);
        double s = 0;
        for (const auto& kv : rec_.images)
            if (kv.second.registered) {
                Vec3 d = cameraCenter(kv.second.pose) - c;
                s += d.dot(d);
            }
        return scale_cache_ = std::max(1e-12, std::sqrt(s / (double)n));
    }

    // The same scale for a model the mapper has not adopted (checkSeam works
    // on a candidate merge, which is nobody's `rec_` yet).
    static double modelScaleOf(const Reconstruction& m) {
        Vec3 c{0, 0, 0};
        size_t n = 0;
        for (const auto& kv : m.images)
            if (kv.second.registered) { c = c + cameraCenter(kv.second.pose); n++; }
        if (!n) return 1.0;
        c = c * (1.0 / (double)n);
        double s = 0;
        for (const auto& kv : m.images)
            if (kv.second.registered) {
                Vec3 d = cameraCenter(kv.second.pose) - c;
                s += d.dot(d);
            }
        return std::max(1e-12, std::sqrt(s / (double)n));
    }

    size_t unclaimed() const { return unclaimedImages(); }
    size_t numImages() const { return db_.images.size(); }
    const MapperOptions& options() const { return opt_; }
    MapperOptions& options() { return opt_; }

private:
    // Load an existing model into `rec_`: its poses, its cameras (whose focals
    // are then facts, not guesses), and its points re-added as fresh tracks.
    // Images the model does not hold keep the cleared state resetModel() left.
    void adopt(const Reconstruction& m) {
        size_t cam_mismatch = 0, missing = 0, name_mismatch = 0;
        for (const auto& kv : m.cameras) {
            rec_.cameras[kv.first] = kv.second;
            // cameras.bin cannot carry pixel_scale (it says nothing about the
            // lens), so an adopted model arrives with 1.0. Restore what the
            // features say, or every threshold silently reverts to source
            // pixels on --resume (D47).
            auto d = default_cams_.find(kv.first);
            if (d != default_cams_.end()) rec_.cameras[kv.first].pixel_scale = d->second.pixel_scale;
            focal_known_.insert(kv.first);
        }
        for (const auto& kv : m.images) {
            if (!kv.second.registered) continue;
            auto it = rec_.images.find(kv.first);
            if (it == rec_.images.end()) { missing++; continue; }
            if (it->second.camera_id != kv.second.camera_id) cam_mismatch++;
            // Image ids are positions in this database. A model from a
            // *different* database would adopt cleanly and silently reconstruct
            // nonsense, so the names are checked (the model's carry an
            // extension, the database's do not).
            if (!kv.second.name.empty()) {
                const std::string& want = it->second.name;
                if (kv.second.name.compare(0, want.size(), want) != 0 ||
                    (kv.second.name.size() > want.size() && kv.second.name[want.size()] != '.'))
                    name_mismatch++;
            }
            it->second.pose = kv.second.pose;
            it->second.registered = true;
        }
        rec_.points3D.clear();
        rec_.next_point3D_id = 1;
        for (const auto& kv : m.points3D) {
            std::vector<TrackElement> tr;
            for (const TrackElement& e : kv.second.track) {
                auto it = rec_.images.find(e.image_id);
                if (it == rec_.images.end() || !it->second.registered) continue;
                if (e.point2D_idx >= it->second.point3D_ids.size()) continue;
                if (it->second.point3D_ids[e.point2D_idx] != kInvalidPoint3D) continue;
                tr.push_back(e);
            }
            if (tr.size() >= 2) {
                uint64_t id = rec_.addPoint3D(kv.second.xyz, tr);
                rec_.points3D[id].rgb[0] = kv.second.rgb[0];
                rec_.points3D[id].rgb[1] = kv.second.rgb[1];
                rec_.points3D[id].rgb[2] = kv.second.rgb[2];
            }
        }
        if (name_mismatch)
            fprintf(stderr,
                    "[map] WARNING: %zu adopted image(s) have a different name than the database "
                    "entry with the same id -- the model was probably built from other matches, "
                    "and adopting it will produce nonsense\n", name_mismatch);
        if ((missing || cam_mismatch) && opt_.verbose)
            fprintf(stderr,
                    "[map] adopted model: %zu image(s) not in this database, %zu with a "
                    "different camera group\n", missing, cam_mismatch);
    }

    // Sort, report and return. Split out only because run() has two exits.
    std::vector<Reconstruction> finishRun(std::vector<Reconstruction>& models,
                                          std::chrono::steady_clock::time_point prof_start) {
        // COLMAP orders the written models by 3D point count, descending
        // (ReconstructionManager::Write); sparse/0 is therefore the model with
        // the most structure, not the first one found.
        std::stable_sort(models.begin(), models.end(),
                         [](const Reconstruction& a, const Reconstruction& b) {
                             return a.points3D.size() > b.points3D.size();
                         });
        if (opt_.verbose) {
            // Distinct, not the sum: models overlap by design (up to
            // max_model_overlap), so summing would double-count the joins.
            std::set<uint32_t> covered;
            for (const Reconstruction& m : models)
                for (const auto& kv : m.images)
                    if (kv.second.registered) covered.insert(kv.first);
            fprintf(stderr, "[map] done: %zu model(s) covering %zu/%zu distinct images\n",
                    models.size(), covered.size(), db_.images.size());
            for (size_t i = 0; i < models.size(); i++)
                fprintf(stderr, "[map]   model %zu: %u images, %zu points\n", i,
                        models[i].numRegistered(), models[i].points3D.size());
        }
        g_map_prof.report(std::chrono::duration<double>(
                              std::chrono::steady_clock::now() - prof_start).count());
        return models;
    }

private:
    // ---- multi-model bookkeeping (D41) ----
    // `rec_` as a standalone reconstruction: its registered images only, the
    // cameras they use, points colored. `rec_` itself keeps every image --
    // resetModel() reuses those records for the next attempt, so the pruning
    // happens on the copy.
    //
    // Pruning is what makes keeping several models affordable: `rec_.images`
    // carries a points2D / point3D_ids pair for *every* image in the dataset
    // (~200 kB each at 8192 features), so an unpruned 2000-image model costs
    // ~400 MB whether it registered 2000 images or 20.
    Reconstruction snapshotModel() const {
        Reconstruction m = rec_;
        for (auto it = m.images.begin(); it != m.images.end();)
            it = it->second.registered ? std::next(it) : m.images.erase(it);
        // Cameras that ended with no registered image (an unused resolution
        // bucket, or a group whose images were all rejected) would be written
        // to cameras.bin with their default-guess intrinsics -- drop them.
        std::set<uint32_t> used;
        for (const auto& kv : m.images) used.insert(kv.second.camera_id);
        for (auto it = m.cameras.begin(); it != m.cameras.end();)
            it = used.count(it->first) ? std::next(it) : m.cameras.erase(it);
        assignColors(m);
        return m;
    }

    // Record the images a kept model registered. They stop being seed
    // candidates; a later model may still re-register them, up to
    // max_model_overlap of them, which is the overlap a merge step aligns on.
    void claimImages(const Reconstruction& m) {
        if (model_count_.size() != db_.images.size()) model_count_.assign(db_.images.size(), 0);
        for (const auto& kv : m.images)
            if (kv.second.registered) model_count_[kv.first]++;
    }

    size_t unclaimedImages() const {
        size_t n = 0;
        for (uint32_t c : model_count_) n += (c == 0);
        return n;
    }

    bool claimed(uint32_t img) const {
        return img < model_count_.size() && model_count_[img] > 0;
    }

    // How many of `m`'s images some already-kept model also holds.
    size_t overlapWithKept(const Reconstruction& m) const {
        if (model_count_.empty()) return 0;
        size_t n = 0;
        for (const auto& kv : m.images)
            if (kv.second.registered && model_count_[kv.first] > 0) n++;
        return n;
    }

    // Is `m` worth a directory of its own, given what is already kept? Two
    // ways it can fail to be, and both were observed on real-world datasets.
    // So: bounded overlap *and* at least a model's worth of genuinely new
    // images. `why` gets a reason for the log when the answer is no.
    bool admitModel(const Reconstruction& m, std::string& why) const {
        const uint32_t reg = m.numRegistered();
        const size_t ov = overlapWithKept(m);
        const size_t fresh = reg - ov;
        char buf[160];
        if (ov > (size_t)std::max(0, opt_.max_model_overlap)) {
            snprintf(buf, sizeof buf, "%u images, %zu of them already in a kept model",
                     reg, ov);
            why = buf;
            return false;
        }
        if (fresh < (size_t)opt_.min_model_size) {
            snprintf(buf, sizeof buf, "%u images but only %zu not already covered", reg, fresh);
            why = buf;
            return false;
        }
        return true;
    }

    // Same question for the model under construction. Free during the primary
    // model (nothing is claimed yet); a linear scan afterwards, which is noise
    // next to the PnP it gates.
    size_t sharedRegistered() const {
        if (model_count_.empty()) return 0;
        size_t n = 0;
        for (const auto& kv : rec_.images)
            if (kv.second.registered && model_count_[kv.first] > 0) n++;
        return n;
    }

    // Color each point by the mean of its observations' keypoint colors (COLMAP's
    // ExtractColorsForAllImages, but the samples were taken once at extraction so
    // no image is decoded again here). Points whose features carry no color keep
    // the neutral-gray default.
    void assignColors(Reconstruction& rec) const {
        for (auto& kv : rec.points3D) {
            Point3D& p = kv.second;
            uint32_t acc[3] = {0, 0, 0}, n = 0;
            for (const TrackElement& e : p.track) {
                const FeatureSet& fs = feats_[e.image_id];
                if (!fs.hasColors()) continue;
                const uint8_t* c = &fs.colors[(size_t)e.point2D_idx * 3];
                acc[0] += c[0]; acc[1] += c[1]; acc[2] += c[2]; n++;
            }
            if (n) {
                p.rgb[0] = (uint8_t)((acc[0] + n / 2) / n);
                p.rgb[1] = (uint8_t)((acc[1] + n / 2) / n);
                p.rgb[2] = (uint8_t)((acc[2] + n / 2) / n);
            }
        }
    }

    // Grow the current model until nothing else registers, or until it has
    // absorbed max_model_overlap images that a previously-kept model already
    // holds (D41; COLMAP's max_model_overlap break). During the primary model
    // nothing is claimed yet, so the overlap test is inert and this is the
    // pre-D41 loop exactly.
    void grow() {
        growLoop();
        checkedRefine(true);
    }

    // The registration loop alone, returning how many images it brought in.
    // Split from grow() so a continuation pass can skip the final refinement
    // when it registered nothing -- refining a model that did not change is
    // both wasted time and a silent perturbation of a finished result.
    // `max_reg` caps the model's size (0 = grow until nothing registers); the
    // focal bootstrap uses it to build a model just big enough to score.
    uint32_t growLoop(uint32_t max_reg = 0) {
        uint32_t registered_here = 0;
        double next_ba = 3;
        recent_regs_.clear();
        rebuildScores();
        while (true) {
            if (max_reg && rec_.numRegistered() >= max_reg) break;
            if (sharedRegistered() >= (size_t)std::max(1, opt_.max_model_overlap)) {
                if (opt_.verbose)
                    fprintf(stderr, "[map] sub-model reached %d images shared with an earlier "
                            "model; stopping its growth\n", opt_.max_model_overlap);
                break;
            }
            // COLMAP's shape: rank all candidates, try them in order until one
            // registers, and only then recompute the ranking. A single failure
            // must not retire an image -- it usually just means not enough of
            // the scene is triangulated *yet* (see D15).
            std::vector<uint32_t> cands;
            {
                ProfTimer pt(g_map_prof.choose);
                g_map_prof.n_choose++;
                cands = chooseNextImages();
            }
            if (cands.empty()) break;
            bool registered = false;
            for (uint32_t img : cands) {
                reg_trials_[img]++;
                g_map_prof.n_reg_try++;
                bool ok;
                {
                    ProfTimer pt(g_map_prof.reg);
                    ok = registerImage(img);
                }
                if (ok) {
                    g_map_prof.n_reg_ok++;
                    registered = true;
                    registered_here++;
                    ProfTimer pt(g_map_prof.tri);
                    triangulateForImage(img);
                    recent_regs_.push_back(img);
                    break;
                }
            }
            if (!registered) break;  // nothing in the ranking can be registered
            if (rec_.numRegistered() >= next_ba) {
                checkedRefine(false);
                // Refinement mutates observations wholesale (filtering,
                // retriangulation, de-registration, possibly a snapshot
                // restore), so the incremental score cache starts over.
                rebuildScores();
                // De-registration may have shrunk the model; the next trigger
                // is always relative to what actually survived.
                next_ba = std::ceil(rec_.numRegistered() * opt_.ba_growth_ratio);
            }
        }
        return registered_here;
    }

    // Transactional global refinement (D36): one toxic registration reaching a
    // trivial-loss global BA can bend a small model so far that the filters
    // shred it, seed included. Snapshot first; if refinement collapses the
    // model, restore the snapshot and de-register the images added since the
    // last refinement -- the suspects -- instead of keeping the wreckage. This
    // is the mapper "going back": the registrations are undone, the images
    // keep their remaining trials, and growth continues from known-good state.
    void checkedRefine(bool final_pass) {
        uint32_t reg_before = rec_.numRegistered();
        Reconstruction snap;
        std::set<uint32_t> fk;
        {
            ProfTimer pt(g_map_prof.snapshot);
            snap = rec_;
            fk = focal_known_;
        }
        globalRefine(final_pass);
        // Collapse = losing half the registered images. Observation loss alone
        // is NOT a collapse: shredding most of a junk-heavy image's points
        // while every pose survives is the filters working as intended, and
        // vegetation-grade bootstraps routinely shed 2/3 of their observations.
        bool collapsed = 2 * rec_.numRegistered() < reg_before;
        if (collapsed && !recent_regs_.empty()) {
            if (opt_.verbose)
                fprintf(stderr,
                        "[map] refinement collapsed the model (%u -> %u images); undoing %zu "
                        "recent registration(s)\n", reg_before, rec_.numRegistered(),
                        recent_regs_.size());
            rec_ = std::move(snap);
            focal_known_ = std::move(fk);
            for (uint32_t img : recent_regs_) deregisterImage(img);
            resetOrphanCameras();
        }
        recent_regs_.clear();
    }

    // ---- helpers ----
    Vec2 kp(uint32_t img, uint32_t f) const {
        const Keypoint& k = feats_[img].keypoints[f];
        return {k.x, k.y};
    }
    // Images may have different intrinsics *and* different resolutions, so
    // every projection has to go through the image's own camera.
    const Camera& camOf(uint32_t img) const {
        return rec_.cameras.at(rec_.images.at(img).camera_id);
    }
    // A pixel threshold in image `img`'s own pixels. Thresholds are given in
    // extraction pixels, where feature noise lives; this is the only place that
    // knows the difference (D47, Camera::pixel_scale).
    double errPx(uint32_t img) const { return camOf(img).errPx(opt_.max_reproj_error); }
    // The same conversion for a whole-problem scalar the GPU solver can only
    // take one of: the median over the registered cameras. Exact on the common
    // case (one extraction scale for the capture) and a compromise only when a
    // dataset mixes resolutions *and* crosses --max-image-size in one run.
    double medianPixelScale() const {
        std::vector<double> v;
        for (const auto& kv : rec_.images)
            if (kv.second.registered) v.push_back(camOf(kv.first).pixel_scale);
        if (v.empty()) {
            for (const auto& kv : rec_.cameras) v.push_back(kv.second.pixel_scale);
            if (v.empty()) return 1.0;
        }
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        return v[v.size() / 2];
    }
    // Unit viewing ray for feature f of image img (the geometry core's
    // interchange type, D31; fisheye-ready, unlike z=1 normalized coords).
    Vec3 bearing(uint32_t img, uint32_t f) const { return camOf(img).bearing(kp(img, f)); }
    Mat34 Pmat(uint32_t img) const {
        const Pose& p = rec_.images.at(img).pose;
        return {p.R[0], p.R[1], p.R[2], p.t.x, p.R[3], p.R[4], p.R[5], p.t.y,
                p.R[6], p.R[7], p.R[8], p.t.z};
    }
    // A point is "in front" of a camera if it lies along the viewing ray. For
    // the pinhole family that is p.z > 0 (kept bit-identical); a wide-FOV camera
    // sees past 90 deg, where p.z < 0 is still valid, so the test becomes the
    // sign of p . bearing (D33). `b` is the observation's unit bearing (already
    // computed by the caller for the geometry, so this adds no unprojection).
    static bool inFront(const Camera& cam, const Vec3& pc, const Vec3& b, double zmin) {
        return cam.wideFov() ? pc.dot(b) > 0 : pc.z >= zmin;
    }
    double reprojErr(uint32_t img, uint32_t f, const Vec3& X) const {
        const Pose& p = rec_.images.at(img).pose;
        const Camera& cam = camOf(img);
        Vec3 pc = mul(p.R, X) + p.t;
        // Cheirality: pinhole path stays exactly `pc.z < 1e-8` (no bearing cost);
        // a wide-FOV camera sees past 90 deg, so it tests the sign along the ray
        // (D33). For a spherical camera every direction is in view, and the dot
        // test costs nothing beyond the bearing it already needs.
        if (cam.wideFov()) {
            if (pc.dot(cam.bearing(kp(img, f))) <= 0) return 1e30;
        } else if (pc.z < 1e-8) {
            return 1e30;
        }
        Vec2 px = cam.project(pc);
        Vec2 o = kp(img, f);
        return std::hypot(px.x - o.x, px.y - o.y);
    }

    // Triangulate the correspondence (a,fa)<->(b,fb); accept on cheirality,
    // angle and reprojection. Returns true and fills X on success.
    // `err_scale` tightens the reprojection acceptance (< 1 during
    // retriangulation, so re-created points must beat the filter's kill
    // threshold with margin -- the churn hysteresis, see D36).
    bool triangulatePair(uint32_t a, uint32_t fa, uint32_t b, uint32_t fb, Vec3& X,
                         double err_scale = 1.0) const {
        Vec3 ba = bearing(a, fa), bb = bearing(b, fb);
        X = triangulateDLT(Pmat(a), Pmat(b), ba, bb);
        Vec3 pa = mul(rec_.images.at(a).pose.R, X) + rec_.images.at(a).pose.t;
        Vec3 pb = mul(rec_.images.at(b).pose.R, X) + rec_.images.at(b).pose.t;
        if (!inFront(camOf(a), pa, ba, 1e-6) || !inFront(camOf(b), pb, bb, 1e-6)) return false;
        double ang = triangulationAngle(X, cameraCenter(rec_.images.at(a).pose),
                                        cameraCenter(rec_.images.at(b).pose));
        if (ang * 180.0 / M_PI < opt_.min_tri_angle_deg) return false;
        if (reprojErr(a, fa, X) > err_scale * errPx(a)) return false;
        if (reprojErr(b, fb, X) > err_scale * errPx(b)) return false;
        return true;
    }

    // Idempotent: every public entry point calls it, and only the first does
    // anything. The graph and the per-image records are properties of the
    // database, not of a particular model.
    void ensureSetup() {
        if (!setup_done_) {
            setup();
            setup_done_ = true;
        }
    }

    void setup() {
        if (cam_ids_.size() != db_.images.size()) cam_ids_.assign(db_.images.size(), 1);
        // One Camera per distinct id, sized from the first image that uses it.
        // The pristine default is kept so a camera whose images all get
        // de-registered can start over instead of retrying from bad intrinsics.
        for (uint32_t i = 0; i < db_.images.size(); i++) {
            uint32_t cid = cam_ids_[i];
            if (rec_.cameras.count(cid)) continue;
            auto proto = opt_.initial_cameras.find(cid);
            rec_.cameras[cid] = proto != opt_.initial_cameras.end()
                                    ? proto->second
                                    : Camera::defaultFor(cid, feats_[i].width, feats_[i].height,
                                                         opt_.focal, opt_.camera_model);
            rec_.cameras[cid].id = cid;
            default_cams_[cid] = rec_.cameras[cid];
        }
        // Priors survive every reset; see MapperOptions::known_focal_cameras.
        for (uint32_t cid : opt_.known_focal_cameras)
            if (rec_.cameras.count(cid)) focal_known_.insert(cid);
        for (uint32_t i = 0; i < db_.images.size(); i++) {
            Image im;
            im.id = i;
            im.camera_id = cam_ids_[i];
            im.name = db_.images[i].name;  // feature stem; CLI resolves the real filename
            im.points2D.resize(feats_[i].count());
            im.point3D_ids.assign(feats_[i].count(), kInvalidPoint3D);
            for (uint32_t f = 0; f < feats_[i].count(); f++) im.points2D[f] = kp(i, f);
            rec_.images[i] = im;
        }
        graph_.build(db_, [&] {
            std::vector<uint32_t> nf(db_.images.size());
            for (size_t i = 0; i < db_.images.size(); i++) nf[i] = feats_[i].count();
            return nf;
        }());
        reg_trials_.assign(db_.images.size(), 0);
        // Size the score bookkeeping now: attachObservation can run before the
        // first grow() (the post-seed globalRefine retriangulates), and must
        // never index unallocated rows. Values there are throwaway -- grow()
        // rebuilds before the first ranking.
        rebuildScores();
    }

    // Wipe everything a previous seed attempt built, so the next attempt starts
    // from the same state setup() left behind.
    void resetModel() {
        scale_cache_ = 0;
        rec_.points3D.clear();
        rec_.cameras.clear();
        focal_known_ = opt_.known_focal_cameras;
        for (uint32_t i = 0; i < db_.images.size(); i++) {
            uint32_t cid = cam_ids_[i];
            if (!rec_.cameras.count(cid)) {
                // An earlier model's refined intrinsics if there are any, and
                // then the focal search stays off: those parameters were fitted
                // on hundreds of images, and a new model's first registration
                // has no business overwriting them (D45).
                auto cons = cam_consensus_.find(cid);
                if (cons != cam_consensus_.end()) {
                    rec_.cameras[cid] = cons->second.first;
                    focal_known_.insert(cid);
                } else {
                    rec_.cameras[cid] = default_cams_.at(cid);
                }
            }
            Image& im = rec_.images[i];
            im.registered = false;
            im.pose = {mat3Identity(), {0, 0, 0}};
            im.point3D_ids.assign(feats_[i].count(), kInvalidPoint3D);
        }
        reg_trials_.assign(db_.images.size(), 0);
    }

    // ---- initialization ----
    // Why every candidate seed was turned down. "initialization failed" on its
    // own tells a user nothing they can act on, and the reasons call for
    // opposite responses: `few_inliers` means the matcher or the pair selection
    // came up short, `low_angle` means the capture has no wide baseline to
    // start from, `forward` means it is a dolly/driving shot. Printed once when
    // initialize() gives up.
    struct InitTally {
        size_t candidates = 0, no_pose = 0, config = 0, few_inliers = 0, forward = 0,
               few_points = 0, low_angle = 0;
        double best_angle = 0;      // best median angle any candidate reached
        double best_forward = 2.0;  // most sideways motion any candidate had
    };

    void reportInitFailure() const {
        const InitTally& t = init_tally_;
        fprintf(stderr, "[map] initialization failed: %zu candidate pair(s) tried; "
                "rejected %zu no pose, %zu wrong config, %zu too few inliers, "
                "%zu forward motion, %zu too few points, %zu low angle\n",
                t.candidates, t.no_pose, t.config, t.few_inliers, t.forward, t.few_points,
                t.low_angle);
        if (t.candidates) {
            fprintf(stderr, "[map]   best median triangulation angle %.2f deg, "
                    "most sideways baseline %.3f (cap %.2f)\n",
                    t.best_angle, t.best_forward, opt_.init_max_forward_motion);
            if (t.best_angle < opt_.init_min_tri_angle_deg / 8)
                fprintf(stderr, "[map]   no pair has enough parallax to triangulate on: "
                        "the capture may be a pure rotation, or too small a sweep\n");
        }
    }

    // ---- focal bootstrap for rotation-degenerate captures (D48) ----
    // The largest angle between any two registered cameras' orientations. This
    // is the quantity that decides whether a self-calibrating bundle adjustment
    // can see the focal length at all: for a set of cameras that all point the
    // same way, stretching the scene along the viewing axis and scaling the
    // focal to match reproduces every image exactly, so the focal is a free
    // parameter of the reconstruction (the classical critical motion sequence
    // -- pure translation -- and it does not matter whether the translation is
    // forwards or sideways). Only the distortion terms, whose radial
    // polynomial is not scale-covariant, break the tie, and they break it
    // weakly.
    //   Measured: a straight KITTI drive spans 8.5 deg, while all 18 benchmark
    // scenes span 180 -- there is no threshold-tuning problem here.
    double rotationSpreadDeg() const {
        std::vector<const Mat3*> Rs;
        for (const auto& kv : rec_.images)
            if (kv.second.registered) Rs.push_back(&kv.second.pose.R);
        // Sampled: this is O(n^2) in a pass that runs per focal hypothesis, and
        // the *maximum* of a 60-camera sample is within a degree of the true
        // maximum on every capture we have -- it is a yes/no question.
        const size_t kMax = 60;
        const size_t stride = std::max<size_t>(1, Rs.size() / kMax);
        double worst = 0;
        for (size_t i = 0; i < Rs.size(); i += stride)
            for (size_t j = i + stride; j < Rs.size(); j += stride) {
                const Mat3& A = *Rs[i];
                const Mat3& B = *Rs[j];
                // trace(A^T B), i.e. sum of the columnwise dot products.
                double tr = 0;
                for (int r = 0; r < 3; r++)
                    for (int c = 0; c < 3; c++) tr += A[r * 3 + c] * B[r * 3 + c];
                double cs = std::min(1.0, std::max(-1.0, 0.5 * (tr - 1.0)));
                worst = std::max(worst, std::acos(cs) * 180.0 / M_PI);
            }
        return worst;
    }

    struct FocalTrial {
        bool ok = false;
        double focal = 0;
        uint32_t registered = 0;
        size_t observations = 0;
        double rot_spread = 0;
    };

    // Build a small model from `pm` with every guessed camera group set to
    // `focal`, and report how much mutually consistent structure it supports.
    // The score is the observation count, not the reprojection residual: a
    // wrong focal in a degenerate capture fits its own reconstruction just as
    // tightly (0.29 px either way on KITTI) but cannot make as many
    // correspondences agree with one another, so tracks break and the filters
    // take them out. Measured on 25 KITTI frames, 73.5k observations at the
    // right focal against 52k at 1.9x it, with the residual flat throughout.
    FocalTrial focalTrial(const TwoViewMatches& pm, const std::vector<uint32_t>& cams,
                          double focal, uint32_t cap) {
        for (uint32_t cid : cams) {
            default_cams_[cid].setFocal(focal);
            default_cams_[cid].k1 = default_cams_[cid].k2 = 0;
            default_cams_[cid].p1 = default_cams_[cid].p2 = 0;
        }
        resetModel();
        FocalTrial t;
        t.focal = focal;
        if (!trySeedPair(pm, 0.0, 0, true, false)) return t;
        globalRefine(false);
        growLoop(cap);
        globalRefine(false);
        t.ok = rec_.numRegistered() >= 3;
        t.registered = rec_.numRegistered();
        t.observations = countObservations();
        t.rot_spread = rotationSpreadDeg();
        // The trial refined the focal itself; report where it landed, since
        // that is the value worth carrying forward.
        if (!cams.empty() && rec_.cameras.count(cams[0])) t.focal = rec_.cameras[cams[0]].focal();
        return t;
    }

    // Cameras whose focal is still the geometric guess -- nothing measured it,
    // so nothing is lost by searching it. A camera an earlier model already
    // published intrinsics for (cam_consensus_) counts as measured: those were
    // fitted over that model's images, and resetModel starts every later model
    // from them, so a trial could not depart from them anyway (D45).
    std::vector<uint32_t> guessedFocalCameras() const {
        std::vector<uint32_t> out;
        for (const auto& kv : rec_.cameras)
            if (!opt_.given_focal_cameras.count(kv.first) && !cam_consensus_.count(kv.first) &&
                !kv.second.isSpherical())  // no focal to search for (D49)
                out.push_back(kv.first);
        return out;
    }

    // Pick the starting focal by trial reconstruction whenever nothing measured
    // it. Runs before the primary model and leaves nothing behind but a focal.
    //
    // The search is a descent rather than a grid, because that is what makes it
    // cheap enough to always run: halve the focal, and keep halving only while
    // each step makes materially more of the correspondences agree. A capture
    // whose guess is already in the right basin pays one extra trial model and
    // stops. Everything is grown from the *same* seed pair, so the trials'
    // scores differ only by the focal.
    void bootstrapFocalLength() {
        if (opt_.focal_trials <= 0) return;
        std::vector<uint32_t> cams = guessedFocalCameras();
        if (cams.empty()) return;
        // One camera group only. With several, a single scalar hypothesis would
        // be applied to lenses that need different answers, and the score
        // cannot say which one was wrong -- exactly the reason bootstrapFocal
        // works per group (D46). Two guessed groups means a rig or a mixed
        // collection, and there the per-group EXIF prior or the calibration is
        // the answer.
        if (cams.size() > 1) return;
        const double f0 = rec_.cameras.at(cams[0]).focal();

        // Probe: the ordinary seed, grown just far enough to have an opinion.
        const uint32_t cap = (uint32_t)std::max<size_t>(
            8, std::min<size_t>(opt_.focal_model_size, db_.images.size()));
        size_t from = 0;
        resetModel();
        if (!initialize(from) || !seed_pair_) { restoreAfterBootstrap(f0, cams, f0); return; }
        const TwoViewMatches& pm = *seed_pair_;
        globalRefine(false);
        growLoop(cap);
        globalRefine(false);
        FocalTrial base;
        base.ok = true;
        base.focal = rec_.cameras.at(cams[0]).focal();
        base.registered = rec_.numRegistered();
        base.observations = countObservations();
        base.rot_spread = rotationSpreadDeg();
        if (opt_.verbose)
            fprintf(stderr, "[map] focal probe at %.0f: %u image(s), %zu observation(s), "
                    "orientations span %.1f deg%s\n", f0, base.registered, base.observations,
                    base.rot_spread,
                    base.rot_spread <= opt_.focal_max_rot_spread_deg
                        ? " (too little to determine the focal)" : "");

        FocalTrial best = base;
        bool moved = false;
        int used = 0;
        // Descend. The direction is not symmetric: a too-short focal splays the
        // bearings, which leaves a model BA can walk back up (every KITTI start
        // from 0.36x to 1.0x of the truth converged on it), while a too-long one
        // flattens the model into a self-consistent pancake it cannot climb out
        // of. So the search only ever needs to look down, and overshooting down
        // costs little.
        double f = f0;
        while (used < opt_.focal_trials) {
            f *= 0.5;
            used++;
            FocalTrial t = focalTrial(pm, cams, f, cap);
            const bool better = t.ok && (double)t.observations >
                                        (1.0 + opt_.focal_min_gain) * (double)best.observations;
            if (opt_.verbose)
                fprintf(stderr, "[map]   focal %.0f -> %.0f: %u image(s), %zu observation(s)%s\n",
                        f, t.focal, t.registered, t.observations,
                        !t.ok ? " (failed)" : better ? " (better)" : " (no better, stopping)");
            if (!better) break;
            best = t;
            moved = true;
        }
        // Nothing below the guess helped. Look once above it, in case the guess
        // is the short one -- a telephoto capture, where the guess is 2-4x low.
        if (!moved && used < opt_.focal_trials) {
            FocalTrial t = focalTrial(pm, cams, f0 * 1.6, cap);
            const bool better = t.ok && (double)t.observations >
                                        (1.0 + opt_.focal_min_gain) * (double)base.observations;
            if (opt_.verbose)
                fprintf(stderr, "[map]   focal %.0f -> %.0f: %u image(s), %zu observation(s)%s\n",
                        f0 * 1.6, t.focal, t.registered, t.observations,
                        !t.ok ? " (failed)" : better ? " (better)" : " (no better)");
            if (better) { best = t; moved = true; }
        }
        // Either way the probe's *refined* focal is what carries forward, not the
        // raw guess: the probe is a real reconstruction of `cap` images and its
        // bundle adjustment has already had a say. That matters where the guess
        // is far out but the descent finds no decisive gain.
        if (opt_.verbose)
            fprintf(stderr, "[map] focal %.0f -> %.0f (%s; %zu observations vs %zu at the "
                    "guess)\n", f0, best.focal,
                    moved ? "trial reconstruction" : "probe refinement only",
                    best.observations, base.observations);
        restoreAfterBootstrap(best.focal, cams, f0);
    }

    // Put back everything the trials disturbed, keeping only the focal. The
    // primary model must be built exactly as it would have been with that focal
    // supplied on the command line -- no trial's seed choice, relaxation level
    // or cached geometry may leak into it.
    //
    // The one thing worth keeping is the memoized two-view geometry, and only
    // when the focal did *not* move: those RANSACs are the cost of the candidate
    // scan, and on a dataset that seeds at no level at all (where the scan is
    // exhaustive and the probe learned nothing) throwing them away would make
    // the pipeline pay for the whole ladder twice to reach the same verdict.
    // `probe_focal` is the focal the memoized geometry was computed at (the
    // trials never write the cache), so the cache survives exactly when the
    // committed focal is that same value.
    void restoreAfterBootstrap(double focal, const std::vector<uint32_t>& cams,
                               double probe_focal) {
        for (uint32_t cid : cams) {
            Camera fresh = default_cams_.at(cid);
            fresh.setFocal(focal);
            fresh.k1 = fresh.k2 = fresh.p1 = fresh.p2 = 0;
            default_cams_[cid] = fresh;
        }
        if (std::fabs(focal - probe_focal) > 1e-9) seed_geom_.clear();
        used_seeds_.clear();
        seed_pair_ = nullptr;
        init_relax_ = 0;
        resetModel();
    }

    // `from` is an in/out cursor into the ranked candidate list, so a retry
    // resumes after the seed that just produced too small a model. COLMAP's
    // thresholds (16 deg median angle, 100 inliers) are demanding on purpose:
    // everything is built on the seed. When no pair on the dataset can meet
    // them, halving stepwise is strictly better than not reconstructing (D36).
    //
    // The ladder has two rungs (D48). The first four levels relax the angle and
    // inlier demands with the forward-motion cap in force, which is every
    // dataset that has any sideways pair at all -- they behave exactly as
    // before. The next four repeat the sequence with the cap lifted, for
    // captures where *no* pair is sideways: a dashcam, a dolly shot, a drone
    // flying down a corridor. For those, the cap is not a quality filter but a
    // veto on the whole dataset, and lifting it last means a sideways seed is
    // always preferred when one exists, however weak.
    //   Lifting it is safe only because the criterion the cap approximates --
    // is this pair well enough conditioned to build a model on? -- is measured
    // directly by the seed's median triangulation angle, which is still
    // enforced at every level. A forward pair with 16 deg of median parallax is
    // a better seed than a sideways pair with 2.
    struct InitLevel { double angle; int inliers; bool allow_forward; };

    InitLevel initLevel(int level) const {
        InitLevel l{opt_.init_min_tri_angle_deg, opt_.init_min_inliers, level >= 4};
        for (int r = level % 4; r > 0; r--) {
            l.angle = std::max(l.angle * 0.5, 2.0);
            l.inliers = std::max(l.inliers / 2, 2 * TwoViewOptions().min_num_inliers);
        }
        return l;
    }

    bool initialize(size_t& from) {
        init_tally_ = InitTally();
        // The reached relaxation level is sticky across seed retries (D19's
        // best-of-attempts loop): a fresh attempt resumes at the level that
        // last produced a seed, keeping the `from` cursor's indexing valid.
        const int levels = opt_.init_max_forward_motion >= 1.0 ? 4 : 8;
        for (; init_relax_ < levels; init_relax_++, from = 0) {
            const InitLevel l = initLevel(init_relax_);
            if (init_relax_ && from == 0 && opt_.verbose)
                fprintf(stderr, "[map] no seed at stricter thresholds; relaxing to "
                        "%d inliers / %.0f deg%s\n", l.inliers, l.angle,
                        l.allow_forward ? " / forward motion allowed" : "");
            if (initializeAttempt(from, l.angle, l.inliers, l.allow_forward)) return true;
        }
        return false;
    }

    bool initializeAttempt(size_t& from, double min_ang_deg, int min_inliers,
                           bool allow_forward) {
        // Candidate pairs: verified, non-planar, most inliers first.
        std::vector<const TwoViewMatches*> cand;
        for (const TwoViewMatches& p : db_.pairs)
            if (p.config == (int)TwoViewConfig::Uncalibrated &&
                (int)p.matches.size() >= min_inliers)
                cand.push_back(&p);
        std::sort(cand.begin(), cand.end(),
                  [](auto* a, auto* b) { return a->matches.size() > b->matches.size(); });

        for (size_t ci = from; ci < cand.size(); ci++) {
            const TwoViewMatches* p = cand[ci];
            uint32_t a = p->image1, b = p->image2;
            // Relaxation levels re-scan the (re-sorted) candidate list, so the
            // cursor alone cannot prevent rebuilding a seed a previous attempt
            // already grew and rejected -- burning retry budget on duplicates.
            if (used_seeds_.count({a, b})) continue;
            // A further model must start somewhere no kept model reached, or it
            // would just rebuild the model that is already written (COLMAP
            // FindFirstInitialImage's `num_registrations == 0` rule, D41).
            // Inert while the primary model is being built.
            if (claimed(a) || claimed(b)) continue;
            if (!trySeedPair(*p, min_ang_deg, min_inliers, allow_forward, true)) continue;
            used_seeds_.insert({a, b});
            from = ci + 1;
            return true;
        }
        return false;
    }

    // Build the two-camera seed on one candidate pair and keep it if it clears
    // the level's thresholds; otherwise roll the model back to empty. Split out
    // of the candidate scan so the focal bootstrap can rebuild *the same* pair
    // under a different focal -- comparing hypotheses on one pair is what makes
    // its scores comparable.
    bool trySeedPair(const TwoViewMatches& pm, double min_ang_deg, int min_inliers,
                     bool allow_forward, bool memoize) {
        const TwoViewMatches* p = &pm;
        const uint32_t a = p->image1, b = p->image2;
        // Memoized (D38): relaxation levels and later seed attempts rescan
        // the candidate list, and estimateTwoView (full two-view RANSAC) is by
        // far its cost. The result is identical on every rescan -- cameras are
        // pristine defaults whenever initialize() runs (resetModel precedes
        // it), and the estimator is deterministic. The focal bootstrap passes
        // memoize=false because it is *changing* the cameras between calls,
        // which is exactly the assumption the cache rests on.
        TwoViewGeometry g;
        auto cached = memoize ? seed_geom_.find({a, b}) : seed_geom_.end();
        if (cached != seed_geom_.end()) {
            g = cached->second;
        } else {
            TwoViewOptions tvo;
            tvo.recover_pose = true;
            const Camera& ca = camOf(a);
            const Camera& cb = camOf(b);
            if (ca.wideFov() || cb.wideFov()) {
                // Same reason verification works on bearings (D45): a
                // pinhole seed throws away every wide correspondence, and
                // the seed is what the whole model is built on. The
                // thresholds move from pixels to radians with the focal.
                std::vector<Vec3> b1(p->matches.size()), b2(p->matches.size());
                for (size_t k = 0; k < p->matches.size(); k++) {
                    b1[k] = ca.bearing(kp(a, p->matches[k].idx1));
                    b2[k] = cb.bearing(kp(b, p->matches[k].idx2));
                }
                tvo.ransac.max_error =
                    0.5 * (ca.errRad(tvo.ransac.max_error) + cb.errRad(tvo.ransac.max_error));
                g = estimateTwoViewBearing(b1, b2, tvo);
            } else {
                std::vector<Vec2> q1(p->matches.size()), q2(p->matches.size());
                for (size_t k = 0; k < p->matches.size(); k++) {
                    q1[k] = kp(a, p->matches[k].idx1);
                    q2[k] = kp(b, p->matches[k].idx2);
                }
                tvo.K1 = ca.K();
                tvo.K2 = cb.K();
                g = estimateTwoView(q1, q2, tvo);
            }
            if (memoize) seed_geom_[{a, b}] = g;
        }
        init_tally_.candidates++;
        if (!g.has_pose) { init_tally_.no_pose++; return false; }
        if (g.config != TwoViewConfig::Uncalibrated) { init_tally_.config++; return false; }
        if (g.num_inliers < min_inliers) { init_tally_.few_inliers++; return false; }

        // Near-pure forward motion: the baseline is parallel to the
        // viewing directions, points triangulate on needle-thin cones
        // around the epipole (COLMAP init_max_forward_motion).
        //
        // The veto is about what the *image* covers, not about the motion:
        // a narrow lens pointed along its own baseline sees only the pencil of
        // rays around the epipole. A spherical camera sees the whole sphere, so
        // the points abeam have full parallax and driving straight forward is
        // as well conditioned as any other motion -- the cap would be measuring
        // a degeneracy that is not there (D49).
        double fwd;
        {
            Vec3 c2 = cameraCenter(g.pose).normalized();  // cam1 is at origin
            fwd = std::max(std::fabs(c2.z),
                           std::fabs(c2.x * g.pose.R[6] + c2.y * g.pose.R[7] +
                                     c2.z * g.pose.R[8]));
            init_tally_.best_forward = std::min(init_tally_.best_forward, fwd);
            const bool spherical = camOf(a).isSpherical() && camOf(b).isSpherical();
            if (!allow_forward && !spherical && fwd > opt_.init_max_forward_motion) {
                init_tally_.forward++;
                return false;
            }
        }

        // Set up the two cameras and triangulate the inliers.
        rec_.images[a].pose = {mat3Identity(), {0, 0, 0}};
        rec_.images[a].registered = true;
        rec_.images[b].pose = g.pose;
        rec_.images[b].registered = true;

        std::vector<double> angles;
        int created = 0;
        for (size_t k = 0; k < p->matches.size(); k++) {
            if (!g.inlier_mask[k]) continue;
            uint32_t fa = p->matches[k].idx1, fb = p->matches[k].idx2;
            Vec3 X;
            if (!triangulatePair(a, fa, b, fb, X)) continue;
            double ang = triangulationAngle(X, cameraCenter(rec_.images[a].pose),
                                            cameraCenter(rec_.images[b].pose));
            angles.push_back(ang * 180.0 / M_PI);
            rec_.addPoint3D(X, {{a, fa}, {b, fb}});
            created++;
        }
        double medAng = 0;
        if (!angles.empty()) {
            std::sort(angles.begin(), angles.end());
            medAng = angles[angles.size() / 2];
        }
        init_tally_.best_angle = std::max(init_tally_.best_angle, medAng);
        if (created < std::max(30, min_inliers / 2)) init_tally_.few_points++;
        else if (medAng < min_ang_deg) init_tally_.low_angle++;
        if (created >= std::max(30, min_inliers / 2) && medAng >= min_ang_deg) {
            // The seed cameras are about to be refined by the first global
            // BA; do not let a later focal search overwrite that.
            focal_known_.insert(rec_.images[a].camera_id);
            focal_known_.insert(rec_.images[b].camera_id);
            seed_pair_ = &pm;
            seed_forward_ = fwd;
            if (opt_.verbose)
                fprintf(stderr, "[map] init pair (%u,%u): %d points, median angle %.1f deg, "
                        "baseline %.2f forward\n", a, b, created, medAng, fwd);
            return true;
        }
        // Roll back and let the caller try the next candidate.
        rollbackInit(a, b);
        return false;
    }

    void rollbackInit(uint32_t a, uint32_t b) {
        rec_.points3D.clear();
        rec_.next_point3D_id = 1;
        for (auto& kv : rec_.images) {
            kv.second.registered = false;
            std::fill(kv.second.point3D_ids.begin(), kv.second.point3D_ids.end(), kInvalidPoint3D);
        }
    }

    // ---- registration ----
    // Count 2D-3D correspondences an unregistered image has to the model.
    // Reference implementation: score_cache_ maintains exactly this value
    // incrementally (SSPLAT_SFM_SCORE_CHECK=1 cross-checks them on every ranking).
    int score(uint32_t img) const {
        if (rec_.images.at(img).registered) return -1;
        int n = 0;
        for (uint32_t f = 0; f < feats_[img].count(); f++)
            for (const Correspondence& c : graph_.at(img, f))
                if (rec_.images.at(c.image_id).registered &&
                    rec_.images.at(c.image_id).point3D_ids[c.feature_idx] != kInvalidPoint3D) {
                    n++;
                    break;  // count the feature once
                }
        return n;
    }

    // Incremental next-image scoring (D37). Rescoring every unregistered image
    // after every registration was 1/3 of kitchen279's map time, and almost
    // all of it recomputed unchanged values. Bookkeeping instead:
    // support_[i][f] = number of (registered image, triangulated feature)
    // correspondences of (i, f); score_cache_[i] = #features with support > 0
    // == score(i). attachObservation() updates it when a feature of a
    // registered image joins a 3D point (the only mutation on the pure growth
    // path); everything wholesale (refine, undo, reseed) triggers
    // rebuildScores() instead of being tracked piecemeal.
    void rebuildScores() {
        ProfTimer pt(g_map_prof.choose);
        if (support_.size() != db_.images.size()) {
            support_.resize(db_.images.size());
            for (size_t i = 0; i < db_.images.size(); i++)
                support_[i].assign(feats_[i].count(), 0);
            score_cache_.assign(db_.images.size(), 0);
        } else {
            for (auto& s : support_) std::fill(s.begin(), s.end(), 0);
            std::fill(score_cache_.begin(), score_cache_.end(), 0);
        }
        for (const auto& kv : rec_.images) {
            const Image& im = kv.second;
            if (!im.registered) continue;
            for (uint32_t f = 0; f < (uint32_t)im.point3D_ids.size(); f++)
                if (im.point3D_ids[f] != kInvalidPoint3D) attachObservation(kv.first, f);
        }
    }

    // Feature f of registered image img just joined a 3D point: every
    // correspondence of (img, f) gains one unit of support.
    void attachObservation(uint32_t img, uint32_t f) {
        for (const Correspondence& c : graph_.at(img, f))
            if (++support_[c.image_id][c.feature_idx] == 1) score_cache_[c.image_id]++;
    }

    // Candidates for the next registration, best first: unregistered, still
    // within their trial budget, and seeing at least min_num_pnp_inliers
    // triangulated points.
    std::vector<uint32_t> chooseNextImages() const {
        static const bool score_check = std::getenv("SSPLAT_SFM_SCORE_CHECK") != nullptr;
        std::vector<std::pair<int, uint32_t>> ranked;
        for (uint32_t i = 0; i < db_.images.size(); i++) {
            if (rec_.images.at(i).registered) continue;
            if (reg_trials_[i] >= opt_.max_reg_trials) continue;
            int s = score_cache_[i];
            if (score_check && s != score(i)) {
                fprintf(stderr, "[map] SCORE MISMATCH image %u: cache %d, reference %d\n",
                        i, s, score(i));
                abort();
            }
            if (s >= opt_.min_num_pnp_inliers) ranked.emplace_back(s, i);
        }
        std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) {
            if (a.first != b.first) return a.first > b.first;
            return a.second < b.second;  // stable, index-ordered tie-break
        });
        std::vector<uint32_t> out;
        out.reserve(ranked.size());
        for (const auto& r : ranked) out.push_back(r.second);
        return out;
    }

    bool registerImage(uint32_t img) {
        // Gather 2D-3D correspondences (one 3D point per feature: the first seen).
        std::vector<Vec3> X;
        std::vector<Vec3> br;   // observed unit bearings
        std::vector<uint32_t> feat;
        std::vector<uint64_t> pid;
        for (uint32_t f = 0; f < feats_[img].count(); f++) {
            uint64_t chosen = kInvalidPoint3D;
            for (const Correspondence& c : graph_.at(img, f)) {
                const Image& oi = rec_.images.at(c.image_id);
                if (oi.registered && oi.point3D_ids[c.feature_idx] != kInvalidPoint3D) {
                    chosen = oi.point3D_ids[c.feature_idx];
                    break;
                }
            }
            if (chosen == kInvalidPoint3D) continue;
            X.push_back(rec_.points3D[chosen].xyz);
            br.push_back(bearing(img, f));
            feat.push_back(f);
            pid.push_back(chosen);
        }
        if ((int)X.size() < opt_.min_num_pnp_inliers) return false;

        const uint32_t cid = rec_.images[img].camera_id;
        PnPResult r;
        double swept_f = 0;  // focal chosen by the search, 0 = no search ran
        // A supplied focal does *not* switch the sweep off, though the argument
        // for switching it off is tempting (it is a measurement, and a sweep on
        // one image's correspondences once took a dataset's camera
        // from 549 to 5493 px). Measured, gating on it cost that dataset 75
        // images of the primary model and 39 of total coverage: the sweep also
        // lets a second camera group depart from a focal that was measured
        // over the first. The runaway is caught by sanitizeCameras and by the
        // joint refinement (D45), so the sweep stays.
        if (focal_known_.count(cid) || opt_.focal_search_samples <= 0) {
            r = ransacPnP(X, br, camOf(img).focal(), errPx(img));
        } else {
            // First image of this camera: its focal is still the no-EXIF guess,
            // and P3P consumes *calibrated* bearings, so a wrong focal fails
            // outright rather than degrading. Sweep hypotheses log-uniformly
            // and keep the one with the most inliers; BA refines from there.
            const Camera& c0 = camOf(img);
            const int N = opt_.focal_search_samples;
            for (int s = 0; s < N; s++) {
                double t = N > 1 ? (double)s / (N - 1) : 0.0;
                double ratio = opt_.min_focal_ratio *
                               std::pow(opt_.max_focal_ratio / opt_.min_focal_ratio, t);
                Camera trial = c0;
                trial.setFocal(c0.focal() * ratio);
                trial.k1 = trial.k2 = trial.p1 = trial.p2 = 0;
                std::vector<Vec3> bt(feat.size());
                for (size_t k = 0; k < feat.size(); k++) bt[k] = trial.bearing(kp(img, feat[k]));
                PnPResult t_r = ransacPnP(X, bt, trial.focal(), errPx(img));
                if (t_r.success && t_r.num_inliers > r.num_inliers) { r = t_r; swept_f = trial.focal(); }
            }
        }
        // Acceptance gates (COLMAP's): enough inliers *and* enough of the
        // offered correspondences agreeing. An image whose hundreds of 2D-3D
        // candidates yield a bare-minimum consensus is a misregistration
        // waiting to bend the model (D36).
        if (!r.success || r.num_inliers < opt_.min_num_pnp_inliers) return false;
        if ((double)r.num_inliers < opt_.min_pnp_inlier_ratio * (double)X.size()) return false;

        // Commit the searched focal only after the gates pass, then refine the
        // pose on the inlier set (COLMAP refines every accepted PnP pose; the
        // raw P3P/DLT pose is what the next triangulations would build on).
        // A swept focal is only a coarse grid hypothesis, so it is refined
        // jointly with the pose -- a 20% focal error passes RANSAC fine but
        // locks bad intrinsics into the group, and BA then papers over the
        // mismatch with runaway distortion (D36).
        if (swept_f > 0) {
            const double f0 = camOf(img).focal();
            rec_.cameras[cid].setFocal(swept_f);
            for (size_t k = 0; k < feat.size(); k++) br[k] = bearing(img, feat[k]);
            double fs = 1.0;
            refinePose(X, br, r.inlier_mask, r.pose, &fs);
            double refined = swept_f * fs;
            double ratio = refined / default_cams_.at(cid).focal();
            if (fs != 1.0 && ratio > opt_.min_focal_ratio && ratio < opt_.max_focal_ratio) {
                rec_.cameras[cid].setFocal(refined);
                for (size_t k = 0; k < feat.size(); k++) br[k] = bearing(img, feat[k]);
            }
            if (opt_.verbose && std::fabs(camOf(img).focal() - f0) > 1e-6)
                fprintf(stderr, "[map] camera %u focal %.0f -> %.0f (search+refine, %d inliers)\n",
                        cid, f0, camOf(img).focal(), r.num_inliers);
        } else {
            refinePose(X, br, r.inlier_mask, r.pose);
        }
        // Re-classify against the refined pose; the gates apply to the final
        // consensus, not the RANSAC one.
        double thr = camOf(img).errRad(opt_.max_reproj_error);
        thr *= thr;
        r.num_inliers = 0;
        for (size_t k = 0; k < X.size(); k++) {
            r.inlier_mask[k] = pnpResidualSq(r.pose, X[k], br[k]) < thr;
            r.num_inliers += r.inlier_mask[k] ? 1 : 0;
        }
        if (r.num_inliers < opt_.min_num_pnp_inliers) return false;
        if ((double)r.num_inliers < opt_.min_pnp_inlier_ratio * (double)X.size()) return false;
        focal_known_.insert(cid);

        rec_.images[img].pose = r.pose;
        rec_.images[img].registered = true;

        // Track continuation: attach inlier features to their 3D points.
        for (size_t k = 0; k < feat.size(); k++) {
            if (!r.inlier_mask[k]) continue;
            uint32_t f = feat[k];
            if (rec_.images[img].point3D_ids[f] != kInvalidPoint3D) continue;
            Point3D& pt = rec_.points3D[pid[k]];
            if (reprojErr(img, f, pt.xyz) > errPx(img)) continue;
            pt.track.push_back({img, f});
            rec_.images[img].point3D_ids[f] = pid[k];
            attachObservation(img, f);
        }
        if (opt_.verbose)
            fprintf(stderr, "[map] registered image %u (%d/%zu PnP inliers), %u total\n", img,
                    r.num_inliers, X.size(), rec_.numRegistered());
        return true;
    }

    // Join this image's features to 3D points the model already has, wherever
    // its pose explains them: registerImage's track-continuation step, split
    // out so a pose set from outside (the audit's repair) can use it too.
    // Returns the number of observations attached.
    uint32_t attachExisting(uint32_t img) {
        uint32_t attached = 0;
        Image& im = rec_.images[img];
        for (uint32_t f = 0; f < feats_[img].count(); f++) {
            if (im.point3D_ids[f] != kInvalidPoint3D) continue;
            for (const Correspondence& c : graph_.at(img, f)) {
                if (c.image_id == img) continue;
                const Image& oi = rec_.images.at(c.image_id);
                if (!oi.registered) continue;
                uint64_t pid = oi.point3D_ids[c.feature_idx];
                if (pid == kInvalidPoint3D) continue;
                auto pt = rec_.points3D.find(pid);
                if (pt == rec_.points3D.end()) continue;
                if (reprojErr(img, f, pt->second.xyz) > errPx(img)) continue;
                bool dup = false;  // one observation per image on a track
                for (const TrackElement& e : pt->second.track)
                    if (e.image_id == img) { dup = true; break; }
                if (dup) continue;
                pt->second.track.push_back({img, f});
                im.point3D_ids[f] = pid;
                attachObservation(img, f);
                attached++;
                break;
            }
        }
        return attached;
    }

    // ---- triangulation of new points seen by a freshly registered image ----
    void triangulateForImage(uint32_t img, double err_scale = 1.0) {
        for (uint32_t f = 0; f < feats_[img].count(); f++) {
            if (rec_.images[img].point3D_ids[f] != kInvalidPoint3D) continue;
            // candidate observations: registered, feature not yet on a 3D point
            std::vector<Correspondence> obs;
            for (const Correspondence& c : graph_.at(img, f))
                if (rec_.images.at(c.image_id).registered &&
                    rec_.images[c.image_id].point3D_ids[c.feature_idx] == kInvalidPoint3D)
                    obs.push_back(c);
            if (obs.empty()) continue;

            // Triangulate with the correspondence of maximum parallax.
            Vec3 bestX;
            double bestAng = -1;
            const Correspondence* bestC = nullptr;
            for (const Correspondence& c : obs) {
                Vec3 X;
                if (!triangulatePair(img, f, c.image_id, c.feature_idx, X, err_scale)) continue;
                double ang = triangulationAngle(X, cameraCenter(rec_.images[img].pose),
                                                cameraCenter(rec_.images[c.image_id].pose));
                if (ang > bestAng) { bestAng = ang; bestX = X; bestC = &c; }
            }
            if (!bestC) continue;

            // Build the track: this obs + every candidate that reprojects well.
            std::vector<TrackElement> track;
            track.push_back({img, f});
            for (const Correspondence& c : obs) {
                if (rec_.images[c.image_id].point3D_ids[c.feature_idx] != kInvalidPoint3D) continue;
                if (reprojErr(c.image_id, c.feature_idx, bestX) <= err_scale * errPx(c.image_id))
                    track.push_back({c.image_id, c.feature_idx});
            }
            if (track.size() < 2) continue;
            // Guard against two features from the same image on one track.
            std::sort(track.begin(), track.end(),
                      [](const TrackElement& a, const TrackElement& b) {
                          return a.image_id < b.image_id;
                      });
            track.erase(std::unique(track.begin(), track.end(),
                                    [](const TrackElement& a, const TrackElement& b) {
                                        return a.image_id == b.image_id;
                                    }),
                        track.end());
            if (track.size() < 2) continue;
            rec_.addPoint3D(bestX, track);
            // Every track element is a registered image's feature joining a
            // point. (Also reached from completeAndRetriangulate, where the
            // counts drift against its direct track edits -- harmless, no
            // ranking happens before the post-refine rebuildScores().)
            for (const TrackElement& e : track) attachObservation(e.image_id, e.point2D_idx);
        }
    }

    // ---- global refinement: BA + filtering + image de-registration ----
    //
    // COLMAP's IterativeGlobalRefinement in miniature: bundle-adjust, filter,
    // and repeat while the model keeps changing, then de-register images the
    // filtering hollowed out. De-registration is the mapper's undo: a
    // registration that passed its PnP gates but stopped agreeing with the
    // refined model exits (and may re-register later, better) instead of
    // staying and bending everything around it (D36).
    void globalRefine(bool final_pass) {
        if (rec_.numRegistered() < 2 || rec_.points3D.size() < 10) return;
        int rounds = final_pass ? opt_.ba_max_refinements : 2;
        for (int i = 0; i < rounds; i++) {
            // Observations shredded by the previous round's filtering (or
            // never triangulated because the poses were still rough) get a
            // second chance against the refined geometry -- COLMAP's
            // Retriangulate + CompleteTracks. Without it refinement can only
            // ever LOSE observations, and on sparse match graphs the model
            // starves right after bootstrap (D36).
            if (i > 0 && opt_.retri_scale > 0) {
                ProfTimer pt(g_map_prof.retri);
                completeAndRetriangulate();
            }
            size_t before = countObservations();
            BundleOptions bo;
            bo.real = opt_.ba_real;
            bo.device = opt_.device;
            bo.verbose = false;
            bo.loss = opt_.ba_loss;
            bo.loss_param = (float)(opt_.ba_loss_param * medianPixelScale());
            bo.refine_principal_point = opt_.refine_principal_point || pp_free_;
            bo.pp_min_images = opt_.pp_min_images;
            // Convergence-adaptive iterations (D38): the first round of a
            // growth refine runs to a loose tolerance -- most refines stop
            // there because the model barely changed. A round beyond the
            // first only happens when the previous one moved >ba_refine_change
            // of the observations (retriangulation just rebuilt structure, or
            // filtering shredded it), and those rounds run to the solver's
            // full tolerance: the filters that follow are about to make
            // kill/keep decisions against this geometry, and judging them
            // against a half-converged model is what cost a dataset
            // its tail under D37's fixed cap. Final passes are always tight.
            if (!final_pass && i == 0 && opt_.ba_growth_rtol > 0) {
                bo.rtol = opt_.ba_growth_rtol;
                bo.patience = opt_.ba_growth_patience;
            }
            bo.shared_ctx = &ba_ctx_;
            double cost = runGlobalBA(rec_, bo);
            ProfTimer pt(g_map_prof.filter);
            sanitizeCameras();
            int removedObs = 0, removedPts = 0;
            filterPoints(removedObs, removedPts);
            if (opt_.verbose)
                fprintf(stderr,
                        "[map] global BA (cost %.3e), filtered %d obs / %d points, %zu remain\n",
                        cost, removedObs, removedPts, rec_.points3D.size());
            if (!before || (double)removedObs / (double)before <= opt_.ba_refine_change) break;
        }
        int dropped;
        {
            ProfTimer pt(g_map_prof.filter);
            dropped = filterImages();
        }
        if (dropped && opt_.verbose)
            fprintf(stderr, "[map] de-registered %d image(s) (few points or bogus camera), "
                    "%u remain\n", dropped, rec_.numRegistered());
    }

    // Extend existing tracks to unassigned features that reproject well
    // (transitively: added elements are sources for further completion), then
    // re-run triangulation for every registered image so features whose
    // earlier points were filtered can rebuild them from the refined poses.
    //
    // Everything re-added must clear a *stricter* bar (0.75x) than the filter
    // kills at: without the hysteresis, junk in the borderline band churns --
    // filtered at >4 px, immediately re-created at <=4 px -- and every BA
    // round drags the poses toward it again (D36).
    void completeAndRetriangulate() {
        const double err = opt_.retri_scale * opt_.max_reproj_error;  // x pixel_scale below
        for (auto& kv : rec_.points3D) {
            Point3D& pt = kv.second;
            std::set<uint32_t> in_track;
            for (const TrackElement& e : pt.track) in_track.insert(e.image_id);
            for (size_t ti = 0; ti < pt.track.size(); ti++) {
                const TrackElement e = pt.track[ti];  // copy: track grows below
                for (const Correspondence& c : graph_.at(e.image_id, e.point2D_idx)) {
                    Image& oi = rec_.images.at(c.image_id);
                    if (!oi.registered || in_track.count(c.image_id) ||
                        oi.point3D_ids[c.feature_idx] != kInvalidPoint3D)
                        continue;
                    if (reprojErr(c.image_id, c.feature_idx, pt.xyz) <=
                        camOf(c.image_id).errPx(err)) {
                        pt.track.push_back({c.image_id, c.feature_idx});
                        in_track.insert(c.image_id);
                        oi.point3D_ids[c.feature_idx] = kv.first;
                    }
                }
            }
        }
        for (auto& kv : rec_.images)
            if (kv.second.registered) triangulateForImage(kv.first, opt_.retri_scale);
    }

    size_t countObservations() const {
        size_t n = 0;
        for (const auto& kv : rec_.points3D) n += kv.second.track.size();
        return n;
    }

    // Drop observations that reproject badly and points whose track lost its
    // parallax: a track whose best view pair subtends less than min_tri_angle
    // sits on a near-degenerate cone and feeds PnP unstable geometry (COLMAP
    // filters on both criteria; the old code only checked reprojection).
    void filterPoints(int& removedObs, int& removedPts) {
        std::vector<Vec3> centers(db_.images.size());
        for (uint32_t i = 0; i < db_.images.size(); i++)
            if (rec_.images.at(i).registered) centers[i] = cameraCenter(rec_.images.at(i).pose);
        const double min_ang = opt_.min_tri_angle_deg * M_PI / 180.0;
        std::vector<uint64_t> drop;
        for (auto& kv : rec_.points3D) {
            Point3D& pt = kv.second;
            std::vector<TrackElement> keep;
            for (const TrackElement& e : pt.track) {
                if (reprojErr(e.image_id, e.point2D_idx, pt.xyz) <= errPx(e.image_id))
                    keep.push_back(e);
                else {
                    rec_.images[e.image_id].point3D_ids[e.point2D_idx] = kInvalidPoint3D;
                    removedObs++;
                }
            }
            pt.track = keep;
            bool degenerate = pt.track.size() < 2;
            if (!degenerate) {
                double best = 0;
                for (size_t i = 0; i + 1 < pt.track.size() && best < min_ang; i++)
                    for (size_t j = i + 1; j < pt.track.size() && best < min_ang; j++)
                        best = std::max(best, triangulationAngle(pt.xyz, centers[pt.track[i].image_id],
                                                                 centers[pt.track[j].image_id]));
                degenerate = best < min_ang;
            }
            if (degenerate) {
                for (const TrackElement& e : pt.track) {
                    rec_.images[e.image_id].point3D_ids[e.point2D_idx] = kInvalidPoint3D;
                    removedObs++;
                }
                drop.push_back(kv.first);
            }
        }
        for (uint64_t id : drop) { rec_.points3D.erase(id); removedPts++; }
    }

    // Undo a registration: detach every observation, drop tracks that fall
    // below two views, unregister. reg_trials_ already counted the attempt,
    // so the image can be retried until its budget runs out (D15's ranking
    // naturally re-offers it once more of the scene exists).
    void deregisterImage(uint32_t img) {
        Image& im = rec_.images[img];
        for (uint32_t f = 0; f < (uint32_t)im.point3D_ids.size(); f++) {
            uint64_t id = im.point3D_ids[f];
            if (id == kInvalidPoint3D) continue;
            im.point3D_ids[f] = kInvalidPoint3D;
            auto it = rec_.points3D.find(id);
            if (it == rec_.points3D.end()) continue;
            auto& tr = it->second.track;
            tr.erase(std::remove_if(tr.begin(), tr.end(),
                                    [&](const TrackElement& e) { return e.image_id == img; }),
                     tr.end());
            if (tr.size() < 2) {
                for (const TrackElement& e : tr)
                    rec_.images[e.image_id].point3D_ids[e.point2D_idx] = kInvalidPoint3D;
                rec_.points3D.erase(it);
            }
        }
        im.registered = false;
    }

    // A camera that left the physically plausible regime is pulled back in
    // rather than getting its images de-registered (where COLMAP kills the
    // images, fatal when a camera group covers many of them). If the wild
    // parameters were benign -- high-order distortion terms unconstrained
    // outside the data's radial support -- residuals barely move and
    // everything survives; if they were absorbing bad geometry, the very next
    // reprojection filter exposes exactly the observations that depended on
    // them, and de-registration proceeds from evidence (D36).
    void sanitizeCameras() {
        for (auto& kv : rec_.cameras) {
            Camera& c = kv.second;
            const Camera& d = default_cams_.at(kv.first);
            int fixed = 0;
            double ratio = c.focal() / d.focal();
            if (!(ratio > opt_.min_focal_ratio && ratio < opt_.max_focal_ratio)) {
                c.setFocal(d.focal());
                fixed++;
            }
            for (double* k : {&c.k1, &c.k2, &c.k3, &c.k4, &c.p1, &c.p2, &c.sx1, &c.sy1})
                if (std::fabs(*k) > opt_.max_extra_param) { *k = 0; fixed++; }
            // The mapper never has evidence to move the principal point far
            // from the center (COLMAP does not refine it at all during
            // mapping); a large excursion is BA absorbing something else.
            if (std::fabs(c.cx - d.cx) > 0.2 * c.width) { c.cx = d.cx; fixed++; }
            if (std::fabs(c.cy - d.cy) > 0.2 * c.height) { c.cy = d.cy; fixed++; }
            if (fixed && opt_.verbose)
                fprintf(stderr, "[map] camera %u: reset %d runaway parameter(s)\n",
                        kv.first, fixed);
        }
    }

    // COLMAP's FilterImages, minus the bogus-camera criterion (cameras are
    // sanitized in place instead): de-register images whose observations
    // collapsed under filtering.
    int filterImages() {
        if (rec_.numRegistered() <= 2) return 0;
        std::vector<uint32_t> drop;
        for (auto& kv : rec_.images) {
            Image& im = kv.second;
            if (!im.registered) continue;
            if ((int)im.numPoint3D() < opt_.min_image_points) drop.push_back(kv.first);
        }
        for (uint32_t id : drop) deregisterImage(id);
        if (!drop.empty()) resetOrphanCameras();
        return (int)drop.size();
    }

    // A camera whose registered images all went away starts over from the
    // pristine default: its focal search / BA state was fit to registrations
    // that have been rejected, and a retry must not inherit that.
    void resetOrphanCameras() {
        std::set<uint32_t> used;
        for (const auto& kv : rec_.images)
            if (kv.second.registered) used.insert(kv.second.camera_id);
        for (auto& kv : rec_.cameras)
            if (!used.count(kv.first) && focal_known_.count(kv.first)) {
                kv.second = default_cams_.at(kv.first);
                // A prior is a measurement of the lens, not of the
                // registrations that were just rejected: it stays known, and
                // default_cams_ already holds it.
                if (!opt_.known_focal_cameras.count(kv.first)) focal_known_.erase(kv.first);
            }
    }

    const MatchesDatabase& db_;
    const std::vector<FeatureSet>& feats_;
    MapperOptions opt_;
    std::vector<uint32_t> cam_ids_;   // per image, 1-based camera id
    std::set<uint32_t> focal_known_;  // cameras whose focal is no longer a guess
    bool pp_free_ = false;            // inside polish(): release the principal point
    std::map<uint32_t, Camera> default_cams_;  // pristine per-group defaults
    // Best intrinsics any admitted model has produced for each camera group,
    // with the number of images that constrained them (D45; recordCameras).
    std::map<uint32_t, std::pair<Camera, double>> cam_consensus_;
    bool setup_done_ = false;
    // SSPLAT_SFM_AUDIT_DUMP=1 prints the support ratio of every audited image, which
    // is how the threshold above was chosen against tools/rig_check.py.
    const bool audit_dump_ = std::getenv("SSPLAT_SFM_AUDIT_DUMP") != nullptr;
    mutable double scale_cache_ = 0;  // modelScale(), reset by resetModel()
    int init_relax_ = 0;              // reached seed-threshold relaxation level
    InitTally init_tally_;            // why the last initialize() found nothing
    const TwoViewMatches* seed_pair_ = nullptr;  // pair the last seed was built on
    double seed_forward_ = 0;         // its |baseline . viewing dir|
    std::set<std::pair<uint32_t, uint32_t>> used_seeds_;  // seeds already grown
    std::map<std::pair<uint32_t, uint32_t>, TwoViewGeometry> seed_geom_;  // memoized (D38)
    Reconstruction rec_;
    CorrespondenceGraph graph_;
    std::vector<std::vector<uint16_t>> support_;  // per (image, feature), see rebuildScores
    std::vector<int> score_cache_;                // == score(i) for every image
    std::vector<int> reg_trials_;
    // Per image, how many *kept* models registered it (D41). Empty until the
    // first model is kept, which is what makes the whole multi-model path inert
    // while the primary model is being built.
    std::vector<uint32_t> model_count_;
    // Persistent BA context (D38): device + pipelines survive across the
    // mapper's many global BAs; each solve only creates/frees its own
    // problem-sized buffers. Uninitialized until the first BA initializes it.
    VkContext ba_ctx_;
    std::vector<uint32_t> recent_regs_;  // registered since the last refinement
};

}  // namespace sfm
