// The one option surface of the SfM pipeline: an aggregate of the stage option
// structs, plus a descriptor table that every consumer reads instead of writing
// its own list.
//
// Why a table. The pipeline has ~90 knobs spread over eight option structs, and
// three consumers need the same view of them: the CLI parser, `--help`, and
// (port plan phase 5) the GUI's options editor. Written by hand that is three
// lists that drift; here the table is the source of truth and each consumer is
// one macro expansion over it -- the same shape as the trainer's generated
// SSPLAT_CONFIG_FIELDS / ConfigUI pair.
//
// What is NOT in the table: flags that do not name one scalar field.
// `--camera-model PREFIX=MODEL` and `--focal PREFIX=F` also feed a per-group
// override list, `--no-manage` moves four fields at once, `--check` / `--audit`
// pick what the command does rather than how, and `-o/--output` is positional
// in spirit. Those stay hand-parsed in the CLI, which tries them *before* the
// table, so a hand-parsed name always wins.
//
// Layering: the sub-option structs are untouched, and the pipeline-level knobs
// here are the ones that span stages (`--max-error` is one tolerance measured
// in one frame, D47, but two struct fields; `--device`, `--quiet` and the
// camera settings likewise). finalize() is what fans them out; nothing else
// may, or the CLI and the GUI would have to agree on the fan-out separately.
#pragma once

#include <cmath>
#include <cstdio>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "sfm/core/CameraSetup.h"
#include "sfm/feature/Matcher.h"
#include "sfm/feature/PairSelection.h"
#include "sfm/feature/Sift.h"
#include "sfm/geometry/TwoView.h"
#include "sfm/map/Bottomup.h"
#include "sfm/map/Manager.h"
#include "sfm/map/Mapper.h"
#include "sfm/map/Merge.h"

namespace sfm {

// Which subcommands accept a flag. A name may be reused across commands as long
// as the masks are disjoint -- `--max-error` is the verification/mapping
// tolerance for auto/match/map and the *alignment* tolerance for merge, which
// are genuinely different quantities with the same spelling.
enum CmdMask : uint32_t {
    CMD_AUTO    = 1u << 0,
    CMD_EXTRACT = 1u << 1,
    CMD_MATCH   = 1u << 2,
    CMD_MAP     = 1u << 3,
    CMD_MERGE   = 1u << 4,
    CMD_ALL     = CMD_AUTO | CMD_EXTRACT | CMD_MATCH | CMD_MAP | CMD_MERGE,
};

// Basic is what a first-time user has to think about and what the GUI shows
// unfolded; Advanced is everything else; Alias is a second spelling of a field
// already listed (parsed, but not printed twice or offered to the GUI).
enum class Tier { Basic, Advanced, Alias };

// The aggregate. Sub-option structs are held unchanged so that library callers
// -- the tests, the mapper's own defaults -- are unaffected by anything here.
struct SfmConfig {
    // ---- pipeline-level knobs ----
    // Presets (auto only). Applied before the table's overrides, so an explicit
    // flag always wins; see applyPresets().
    std::string quality = "high";
    std::string data_type = "individual";
    // "auto" resolves at run time: sequential for video, exhaustive below 100
    // images, prefilter at or above. `match` reads it as exhaustive.
    std::string pairs = "auto";
    int overlap = 10;

    // The one geometric tolerance, in extraction pixels (D47): the two-view
    // verifier's inlier radius and the mapper's reprojection cap are the same
    // quantity in the same frame. finalize() writes both.
    double max_error = 3.0;
    int max_image_size = 3200;
    std::string mask_dir;

    // Camera setup. The string forms are what the table and the GUI see; the
    // parsed forms live in `camera` after finalize(). `--camera-model` and
    // `--focal` also accept PREFIX=VALUE, which the CLI routes straight into
    // camera.overrides and never through these.
    std::string camera_mode = "folder";
    std::string camera_model = "opencv";
    double focal = 0;
    // Whether camera_mode is pinned: set by an explicit --camera-mode, and by
    // the internet preset, which chooses per-image cameras deliberately (D20).
    // Unpinned, buildCameras may still switch a Folder default to Image for a
    // set that is plainly a photo collection (D48).
    bool camera_mode_pinned = false;

    // Stage switches that are not a field of any stage's options.
    bool verify = true;                 // match: geometric verification
    bool final_principal_point = true;  // one PP-free global BA at the end (D51)
    bool merge_ba = true;               // merge: bundle-adjust across the seams
    bool in_place = false;              // merge: write back over the input

    // Inputs the `map` subcommand names by flag (auto takes them positionally).
    std::string image_dir;
    std::string feature_dir;
    std::string resume;
    bool check = false;

    // Runtime.
    int threads = 0;           // host worker pools; 0 = hardware_concurrency
    int decode_threads = 0;    // image decode pool; 0 = hardware_concurrency
    int decode_budget_mb = 0;  // 0 = ImageLoadOptions default
    int device = -1;
    bool quiet = false;

    // ---- the stage option structs, unchanged ----
    SiftOptions sift;
    MatchOptions match;
    PairSelectionOptions prefilter;
    TwoViewOptions twoview;
    CameraSetupOptions camera;
    MapperOptions mapper;
    ManagerOptions manager;
    MergeOptions merge;
    BottomUpOptions bup;
    // How the mapper is scheduled: "flat" is one incremental reconstruction of
    // the whole capture, "bottom-up" reconstructs small atoms of the view graph
    // and merges them upwards (D57), and "auto" picks bottom-up at or above
    // bup.min_images, where the flat schedule's whole-model passes start to
    // dominate. "hierarchical" is accepted as the old name for "bottom-up".
    // Flat until the bottom-up schedule has been measured on enough captures.
    std::string mapper_mode = "flat";

    // Fan the pipeline knobs out into the structs above and validate what the
    // table's ranges cannot. Returns an empty string, or the error to print.
    // Call once, after parsing and after the presets.
    std::string finalize(uint32_t cmd);

    // What --pairs names, with "auto" resolving to exhaustive. `auto` alone
    // additionally switches to pair selection above 100 images, which it can
    // only decide once extraction has counted them -- see cmdAuto.
    PairMode pairMode() const;
};

// ---------------------------------------------------------------------------
// The descriptor table
// ---------------------------------------------------------------------------
// F(member, name, cmds, tier, group, lo, hi, choices, help)
//
//   member   path inside SfmConfig; its declared type is the field's type
//            (bool switch, integer, real, string), deduced at each expansion
//   name     the CLI flag without "--"; the GUI's id; the key applyPresets logs
//   cmds     CmdMask bits. A field an editor should offer must include
//            CMD_AUTO, since `auto` is the whole pipeline
//   tier     Tier::Basic is the short list; Tier::Alias is a second spelling
//   group    display group, in the order the rows appear here
//   lo, hi   accepted range; lo >= hi means unbounded
//   choices  "a|b|c" for string fields, "" otherwise
//   help     one sentence, no trailing default (the printer appends it)
//
// A bool field is a switch: `--name` sets it true and `--no-name` false, so a
// field defaulting to true reads as `--no-name` on the command line and is
// printed that way. Keep the name positive.
#define SFM_CONFIG_FIELDS(F)                                                                       \
    /* ---- pipeline ---- */                                                                       \
    F(quality, "quality", CMD_AUTO, Tier::Basic, "pipeline", 0, 0,                                  \
      "low|medium|high|extreme",                                                                    \
      "Working resolution, feature budget and pair-selection breadth")                              \
    F(data_type, "data-type", CMD_AUTO, Tier::Basic, "pipeline", 0, 0,                              \
      "individual|video|internet",                                                                  \
      "What the capture is: individual photos, video frames, or an unordered internet collection")  \
    F(pairs, "pairs", CMD_AUTO | CMD_MATCH, Tier::Basic, "pipeline", 0, 0,                          \
      "auto|exhaustive|sequential|prefilter",                                                       \
      "Which image pairs are matched; auto is sequential for video, exhaustive below 100 images "   \
      "and GPU pair selection at or above")                                                         \
    F(overlap, "overlap", CMD_AUTO | CMD_MATCH, Tier::Advanced, "pipeline", 1, 1000000, "",         \
      "Neighbours each image is paired with under --pairs sequential")                              \
    F(max_error, "max-error", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Advanced, "pipeline",           \
      0.1, 100, "",                                                                                 \
      "Inlier radius for verification and mapping, in pixels of the image SIFT ran on rather "      \
      "than of the source file (D47)")                                                              \
    F(max_image_size, "max-image-size", CMD_AUTO | CMD_EXTRACT, Tier::Advanced, "pipeline",         \
      64, 20000, "",                                                                                \
      "Longest edge SIFT runs on; larger images are downscaled first, and keypoints are still "     \
      "reported in the source image's pixels")                                                      \
    F(mask_dir, "masks", CMD_AUTO | CMD_EXTRACT, Tier::Basic, "pipeline", 0, 0, "",                 \
      "Directory of masks; keypoints on zero (black) pixels are dropped. auto defaults it to "      \
      "`masks` beside the image directory")                                                         \
    F(mask_dir, "mask-dir", CMD_AUTO | CMD_EXTRACT, Tier::Alias, "pipeline", 0, 0, "",              \
      "Alias of --masks")                                                                           \
    /* ---- camera ---- */                                                                          \
    F(camera_mode, "camera-mode", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Basic, "camera", 0, 0,      \
      "single|folder|image",                                                                        \
      "How images are grouped into cameras; every mode splits on image resolution first")           \
    F(camera_model, "camera-model", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Basic, "camera", 0, 0,    \
      "simple-pinhole|pinhole|radial|opencv|full-opencv|opencv-fisheye|thin-prism-fisheye|"         \
      "equirectangular",                                                                            \
      "Distortion model fitted to each camera group; also takes PREFIX=MODEL to set one group")     \
    F(focal, "focal", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Basic, "camera", 0, 10000000, "",       \
      "Starting focal length in pixels, 0 to guess from EXIF or image size; also takes PREFIX=F")   \
    F(camera.exif_focal, "exif-focal", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Advanced, "camera",    \
      0, 0, "",                                                                                     \
      "Use the focal length EXIF recorded when no --focal covers the group")                        \
    F(camera.exif_groups, "exif-groups", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Advanced, "camera",  \
      0, 0, "",                                                                                     \
      "Split camera groups by what EXIF says the body and the focal setting were (D48)")            \
    F(camera.exif_focal_tol, "exif-focal-tol", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Advanced,      \
      "camera", 0.001, 1.0, "",                                                                     \
      "Relative tolerance clustering EXIF focals into one group; must exceed EXIF's 1 mm "          \
      "quantization and stay under a real zoom step")                                               \
    /* ---- features (SIFT) ---- */                                                                 \
    F(sift.max_num_features, "max-features", CMD_AUTO | CMD_EXTRACT, Tier::Advanced, "features",    \
      128, 1000000, "", "Keypoints kept per image, the largest scales first")                       \
    F(sift.num_octaves, "octaves", CMD_AUTO | CMD_EXTRACT, Tier::Advanced, "features", 1, 8, "",    \
      "Scale-space octaves")                                                                        \
    F(sift.peak_threshold, "peak-threshold", CMD_AUTO | CMD_EXTRACT, Tier::Advanced, "features",    \
      0, 1, "", "DoG response a keypoint must reach")                                               \
    F(sift.edge_threshold, "edge-threshold", CMD_AUTO | CMD_EXTRACT, Tier::Advanced, "features",    \
      1, 1000, "", "Principal-curvature ratio above which a keypoint is an edge, not a corner")     \
    F(sift.max_num_orientations, "max-orientations", CMD_AUTO | CMD_EXTRACT, Tier::Advanced,        \
      "features", 1, 8, "", "Descriptors emitted per keypoint when its gradient is ambiguous")      \
    F(sift.profile, "profile", CMD_EXTRACT, Tier::Advanced, "features", 0, 0, "",                   \
      "Print per-stage GPU timings for the extractor")                                              \
    F(sift.spv_path, "spv-path", CMD_EXTRACT, Tier::Advanced, "features", 0, 0, "",                 \
      "Load the SIFT kernels from this SPIR-V file instead of the embedded blob")                   \
    /* ---- matching ---- */                                                                        \
    F(match.max_ratio, "ratio", CMD_AUTO | CMD_MATCH, Tier::Advanced, "matching", 0, 1, "",         \
      "Lowe ratio: a match is kept when the best distance is below this times the second best")     \
    F(match.cross_check, "cross-check", CMD_AUTO | CMD_MATCH, Tier::Advanced, "matching", 0, 0,     \
      "", "Keep only mutual nearest neighbours")                                                    \
    F(match.max_num_matches, "max-matches", CMD_AUTO | CMD_MATCH, Tier::Advanced, "matching",       \
      0, 1000000000, "", "Cap on matches per pair, 0 for no cap")                                   \
    F(verify, "verify", CMD_MATCH, Tier::Advanced, "matching", 0, 0, "",                            \
      "Geometrically verify each pair (F/H RANSAC); off keeps the raw putative matches")            \
    F(twoview.min_num_inliers, "min-inliers", CMD_AUTO | CMD_MATCH, Tier::Advanced, "matching",     \
      4, 1000000, "", "Inliers a verified pair must keep to be recorded")                           \
    F(prefilter.num_features, "prefilter-features", CMD_AUTO | CMD_MATCH, Tier::Advanced,           \
      "matching", 16, 100000, "", "Query descriptors per image in the pair-selection pass")         \
    F(prefilter.train_features, "prefilter-train", CMD_AUTO | CMD_MATCH, Tier::Advanced,            \
      "matching", 0, 1000000, "",                                                                   \
      "Cap on the train side of that pass, 0 to score against every descriptor")                    \
    F(prefilter.num_neighbors, "prefilter-neighbors", CMD_AUTO | CMD_MATCH, Tier::Advanced,         \
      "matching", 1, 100000, "", "Best-scoring partners each image keeps for full matching")        \
    F(prefilter.min_score, "prefilter-min-score", CMD_AUTO | CMD_MATCH, Tier::Advanced,             \
      "matching", 0, 1000000, "", "Mini-matches below which a pair never qualifies")                \
    F(prefilter.ratio, "prefilter-ratio", CMD_AUTO | CMD_MATCH, Tier::Advanced, "matching",         \
      0, 1, "", "Lowe ratio of the scoring pass; looser than the matcher's, as it only ranks")      \
    /* ---- mapper ---- */                                                                          \
    F(mapper.focal_trials, "focal-trials", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",            \
      0, 1000, "",                                                                                  \
      "Trial reconstructions used to pick a focal the motion cannot determine (D48), 0 to skip")    \
    F(mapper.refine_principal_point, "refine-principal-point", CMD_AUTO | CMD_MAP, Tier::Advanced,  \
      "mapper", 0, 0, "",                                                                           \
      "Let BA move cx,cy throughout; off as in COLMAP, where it is nearly a camera rotation (D50)") \
    F(final_principal_point, "final-principal-point", CMD_AUTO | CMD_MAP, Tier::Advanced,           \
      "mapper", 0, 0, "",                                                                           \
      "One principal-point-free global BA on the finished model, for a single camera group (D51)")  \
    F(mapper.pp_min_images, "pp-min-images", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",          \
      2, 1000000, "", "Images a group needs before that final pass runs on it")                     \
    F(mapper.min_tri_angle_deg, "min-tri-angle", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",      \
      0, 90, "", "Triangulation angle a 3D point must subtend to be kept")                          \
    F(mapper.init_min_tri_angle_deg, "init-min-tri-angle", CMD_AUTO | CMD_MAP, Tier::Advanced,      \
      "mapper", 0, 90, "",                                                                          \
      "Median triangulation angle the seed pair must reach; relaxed stepwise if nothing passes")    \
    F(mapper.min_num_pnp_inliers, "min-pnp-inliers", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",  \
      4, 1000000, "", "2D-3D inliers an image needs to register")                                   \
    F(mapper.min_pnp_inlier_ratio, "min-pnp-ratio", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",   \
      0, 1, "", "Inlier fraction an image needs to register, which rejects accidental agreement")   \
    F(mapper.min_image_points, "min-image-points", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",    \
      1, 1000000, "", "Observations below which a registered image is dropped again")               \
    F(mapper.ba_loss, "ba-loss", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 0, 0,                \
      "trivial|huber|cauchy", "Robust loss used by mapping-time bundle adjustment")                 \
    F(mapper.ba_loss_param, "ba-loss-param", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",          \
      0, 1000, "", "Huber delta / Cauchy c, in extraction pixels")                                  \
    F(mapper.ba_solver, "ba-solver", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 0, 0,            \
      "auto|dense|cg",                                                                              \
      "Linear solver for the reduced camera system; auto switches to CG above the size where "      \
      "the dense factorization stops paying")                                                       \
    F(mapper.retri_scale, "retri-scale", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 0, 10, "",   \
      "Retriangulation tolerance as a fraction of --max-error, 0 to skip the pass")                 \
    F(mapper.merge_tracks, "merge-tracks", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 0, 0,      \
      "", "Fuse two 3D points a correspondence says are the same feature")                          \
    F(mapper.rank_by_visibility, "rank-by-visibility", CMD_AUTO | CMD_MAP, Tier::Advanced,          \
      "mapper", 0, 0, "",                                                                           \
      "Rank the next image by how its visible structure spreads over the frame, not by count")      \
    F(mapper.seed_blocking, "seed-blocking", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 0, 0,    \
      "", "A seed retry starts somewhere no earlier attempt reached, instead of rebuilding it")     \
    F(mapper_mode, "mapper", CMD_AUTO | CMD_MAP, Tier::Basic, "mapper", 0, 0,                       \
      "auto|flat|bottom-up",                                                                        \
      "One incremental reconstruction of the whole capture, or small atoms of the view graph "      \
      "reconstructed separately and merged upwards; auto picks the latter for a large capture")     \
    F(bup.min_images, "bup-min-images", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",               \
      2, 1000000, "", "Images at or above which --mapper auto goes bottom-up")                      \
    F(bup.partition.leaf_max_images, "bup-atom-size", CMD_AUTO | CMD_MAP, Tier::Advanced,           \
      "mapper", 8, 100000, "", "Images a view-graph atom is split until it is under")               \
    F(bup.partition.overlap, "bup-overlap", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",           \
      0, 100000, "", "Images each atom borrows from its sibling, which is what a merge aligns on")  \
    F(bup.max_rounds, "bup-rounds", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 1, 100, "",       \
      "Levels of the merge tree before the leftovers are handed to the manage loop")                     \
    F(bup.joint_intrinsics, "bup-joint-intrinsics", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",   \
      0, 0, "", "Bundle-adjust every model in one problem with the intrinsics shared per camera")   \
    F(bup.grow_every, "bup-grow-every", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 0, 100,      \
      "", "Levels between growth passes over the models that did not merge; 0 disables")            \
    F(bup.grow_budget_frac, "bup-grow-budget", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 0,     \
      1000, "", "Images one growth pass may add to a model, as a fraction of what it holds")        \
    F(mapper.ba_growth_ratio, "ba-growth", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper", 1, 100,    \
      "", "Model growth that triggers the next global bundle adjustment")                           \
    F(mapper.ba_growth_rtol, "ba-growth-rtol", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",        \
      0, 1, "", "Relative cost improvement a growth-phase BA stops below, 0 for the solver's")      \
    F(mapper.ba_growth_patience, "ba-growth-patience", CMD_AUTO | CMD_MAP, Tier::Advanced,          \
      "mapper", 1, 1000, "", "Accepted steps below that tolerance before it stops")                 \
    F(mapper.max_num_models, "max-models", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",            \
      1, 100000, "", "Reconstructions a fragmented capture may produce; 1 for a single model")      \
    F(mapper.max_model_overlap, "model-overlap", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",      \
      0, 1000000, "", "Images a further model may take from one already kept -- what a merge "      \
      "later aligns on")                                                                            \
    F(mapper.min_model_size, "min-model-size", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",        \
      2, 1000000, "", "Images a model must reach to be kept at all")                                \
    F(mapper.audit_min_evidence, "audit-evidence", CMD_AUTO | CMD_MAP, Tier::Advanced, "mapper",    \
      0, 1000000, "", "Correspondences an image needs before the audit will judge its pose")        \
    /* ---- manage rounds ---- */                                                                   \
    F(manager.max_rounds, "rounds", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage", 0, 1000, "",      \
      "Merge / grow / prune / reseed rounds run until nothing changes")                             \
    F(manager.do_merge, "merge", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage", 0, 0, "",            \
      "Merge models that share images, on the Sim(3) those shared poses give (D43)")                \
    F(manager.do_grow, "grow", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage", 0, 0, "",              \
      "Register still-unregistered images into the models that exist")                              \
    F(manager.do_reseed, "reseed", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage", 0, 0, "",          \
      "Look for further models among the images nothing covers")                                    \
    F(manager.do_audit, "audit", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage", 0, 0, "",            \
      "Check each pose against the correspondence graph and re-register what a model cannot "       \
      "support (D44)")                                                                              \
    F(manager.do_split, "split", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage", 0, 0, "",            \
      "Break a model its own verified pairs contradict")                                            \
    F(manager.do_duplicate_split, "fold-split", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage",       \
      0, 0, "",                                                                                     \
      "Cut a model that has written two places on top of each other, when the cut severs almost "   \
      "no co-visibility (D46)")                                                                     \
    F(manager.duplicate.max_cut_fraction, "fold-max-cut", CMD_AUTO | CMD_MAP, Tier::Advanced,       \
      "manage", 0, 1, "", "Co-visibility that cut may sever, as a fraction of the model's total")   \
    F(manager.do_joint_ba, "joint-ba", CMD_AUTO | CMD_MAP, Tier::Advanced, "manage", 0, 0, "",      \
      "Refine every component in one problem with intrinsics shared per camera group (D45)")        \
    /* ---- merging ---- */                                                                         \
    F(merge.max_reproj_error, "max-error", CMD_MERGE, Tier::Advanced, "merge", 0.1, 1000, "",       \
      "Alignment inlier threshold in pixels; looser than the mapper's, as the two models were "     \
      "optimized independently")                                                                    \
    F(merge.max_reproj_error, "merge-max-error", CMD_AUTO | CMD_MAP, Tier::Advanced, "merge",       \
      0.1, 1000, "", "Alignment inlier threshold in pixels for the merge step")                     \
    F(merge.min_common_images, "min-common", CMD_MERGE, Tier::Advanced, "merge", 2, 1000000, "",    \
      "Shared images an alignment needs; two determine a Sim(3), three give it a vote")             \
    F(merge.min_common_images, "merge-min-common", CMD_AUTO | CMD_MAP, Tier::Advanced, "merge",     \
      2, 1000000, "", "Shared images an alignment needs in the merge step")                         \
    F(merge.min_inlier_ratio, "min-inlier-ratio", CMD_MERGE, Tier::Advanced, "merge", 0, 1, "",     \
      "Fraction of the shared images an alignment must explain")                                    \
    F(merge.filter_reproj_error, "filter-error", CMD_MERGE, Tier::Advanced, "merge", 0, 1000, "",   \
      "Post-merge observation filtering, in pixels")                                                \
    F(merge.min_tri_angle_deg, "min-tri-angle", CMD_MERGE, Tier::Advanced, "merge", 0, 90, "",      \
      "Post-merge triangulation-angle filtering")                                                   \
    F(merge_ba, "ba", CMD_MERGE, Tier::Advanced, "merge", 0, 0, "",                                 \
      "Bundle-adjust across the seams a merge created; the seam is the one part no BA has seen")    \
    F(in_place, "in-place", CMD_MERGE, Tier::Advanced, "merge", 0, 0, "",                           \
      "Write the merged models back over the input directory")                                      \
    /* ---- inputs ---- */                                                                          \
    F(image_dir, "images", CMD_MAP, Tier::Advanced, "input", 0, 0, "",                              \
      "Image directory, used to put real filenames back into the model")                            \
    F(feature_dir, "features", CMD_MAP, Tier::Advanced, "input", 0, 0, "",                          \
      "Feature directory, if not given positionally")                                               \
    F(resume, "resume", CMD_MAP, Tier::Advanced, "input", 0, 0, "",                                 \
      "Adopt the models in this directory instead of mapping from scratch (D44)")                   \
    F(check, "check", CMD_MAP, Tier::Advanced, "input", 0, 0, "",                                   \
      "With --resume: report how far each model agrees with the two-view geometries it was "        \
      "built from, then exit without writing anything")                                             \
    /* ---- runtime ---- */                                                                         \
    F(threads, "threads", CMD_AUTO | CMD_MATCH | CMD_MAP, Tier::Advanced, "runtime", 0, 4096, "",   \
      "Host worker threads: two-view verification, and the mapper's per-point and per-image "       \
      "passes; 0 is every core, 1 is serial, and results do not depend on the count")               \
    F(decode_threads, "decode-threads", CMD_AUTO | CMD_EXTRACT, Tier::Advanced, "runtime",          \
      0, 4096, "", "Image decode threads; 0 is every core, 1 decodes inline")                       \
    F(decode_budget_mb, "decode-budget", CMD_AUTO | CMD_EXTRACT, Tier::Advanced, "runtime",         \
      0, 1048576, "", "Memory the decode pool may hold in flight, in MB")                           \
    F(device, "device", CMD_ALL, Tier::Advanced, "runtime", -1, 64, "",                             \
      "Vulkan device index; -1 picks the first suitable one")                                       \
    F(quiet, "quiet", CMD_ALL, Tier::Advanced, "runtime", 0, 0, "",                                 \
      "Print only the result lines, not per-stage progress")

// ---------------------------------------------------------------------------
// Consumers
// ---------------------------------------------------------------------------

// Outcome of offering one command-line token to the table.
enum class FieldResult { Unknown, Ok, Error };

// Try to consume argv[i] (and its value, if the field takes one) as a table
// field of `cmd`. `arg` is the whole token, leading dashes included; `seen`
// collects the names that were set, which is what makes an explicit flag beat a
// preset and what tells `map` whether the command line said anything about
// cameras. Returns Unknown if no field of this command has that name.
FieldResult setConfigField(SfmConfig& cfg, uint32_t cmd, const std::string& arg, int argc,
                           char** argv, int& i, std::set<std::string>& seen, std::string& error);

// One field the presets moved, for the "which options a preset changed" report.
struct PresetChange {
    std::string flag, from, to;
};

// Apply --quality / --data-type. Fields in `seen` are left alone, so an
// explicit flag always wins; everything the presets do move is logged.
std::string applyPresets(SfmConfig& cfg, const std::set<std::string>& seen,
                         std::vector<PresetChange>& moved);

// The option block of `--help` for one subcommand: the table's rows for `cmd`,
// grouped, with each field's current value as the printed default.
void printConfigOptions(FILE* out, uint32_t cmd, const SfmConfig& defaults);

// One option line in that same layout, for the flags a CLI parses by hand.
// `flag` carries its own metavar ("-o, --output DIR"); `value` is the bracketed
// note after it, empty for a switch.
void printOptionLine(FILE* out, const std::string& flag, const std::string& value,
                     const std::string& help);

}  // namespace sfm
