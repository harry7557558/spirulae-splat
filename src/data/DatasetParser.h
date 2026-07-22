#pragma once

// DatasetParser -- standalone C++ dataset readers + camera bakes for the
// CLI trainer (app/main.cpp). No Python, no external deps (JSON via
// app/Json.h, PLY reader in NerfstudioParser.cpp, stb handles images in
// the engine's DataManager).
//
// Three formats, mirroring modules/dataparser.py:
//   parse_colmap_dataset      cameras / images / points3D, .bin or .txt
//                             (ColmapParser.cpp)
//   parse_nerfstudio_dataset  transforms.json + PLY point cloud
//                             (NerfstudioParser.cpp)
//   parse_metashape_dataset   Metashape camera-export .xml + .ply (+ optional
//                             .psx project for filename disambiguation)
//                             (MetashapeParser.cpp; XML via app/Xml.h, zips
//                             via external/miniz)
//   parse_dataset             auto-detect: nerfstudio, then COLMAP, then
//                             Metashape (Python probes in the same order)
//
// Both produce a ParsedDataset of PER-INPUT cameras in the raw
// train_frame="points" frame (poses/points as stored; the normalized-frame
// similarity is computed only for the train_frame_scale scalar).
//
// bake_post_split() then expands to the POST-split arrays
// engine_setup_data_manager consumes -- identity (K=1) pass-through, or the
// fisheye/equisolid 5-face / equirectangular 6-face cubemap split when
// warp-to-pinhole is enabled (port of trainer.py _setup_cpp_data_manager).

#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>


// ===========================================================================
// Raw COLMAP records (mirrors colmap_utils.py Camera / Image / Point3D)
// ===========================================================================

struct ColmapCamera {
    int32_t              camera_id = -1;
    std::string          model;       // e.g. "OPENCV", "OPENCV_FISHEYE"
    uint64_t             width  = 0;
    uint64_t             height = 0;
    std::vector<double>  params;      // model-specific, COLMAP order
};

struct ColmapImage {
    int32_t                image_id  = -1;
    std::array<double, 4>  qvec{};    // (w, x, y, z), world->camera
    std::array<double, 3>  tvec{};
    int32_t                camera_id = -1;
    std::string            name;      // path relative to the image dir
};

struct ColmapPoints3D {
    std::vector<float>    xyz;        // [N, 3] flat
    std::vector<uint8_t>  rgb;        // [N, 3] flat
    int64_t num() const { return (int64_t)xyz.size() / 3; }
};

// Low-level readers. Throw std::runtime_error on missing/corrupt files.
std::map<int32_t, ColmapCamera> read_cameras_binary(const std::string& recon_dir);
std::map<int32_t, ColmapImage>  read_images_binary(const std::string& recon_dir);
ColmapPoints3D                  read_points3D_binary(const std::string& recon_dir);

// Text-format equivalents (cameras.txt / images.txt / points3D.txt).
std::map<int32_t, ColmapCamera> read_cameras_text(const std::string& recon_dir);
std::map<int32_t, ColmapImage>  read_images_text(const std::string& recon_dir);
ColmapPoints3D                  read_points3D_text(const std::string& recon_dir);

// PLY point-cloud reader (ascii + binary_little_endian; x/y/z of any float
// or double type, red/green/blue uchar or float). NerfstudioParser.cpp.
ColmapPoints3D read_ply_points(const std::string& path);


// ===========================================================================
// Parser config + baked dataset
// ===========================================================================

struct DatasetParserConfig {
    // COLMAP reconstruction dir relative to the dataset dir. Empty =
    // auto-detect over {sparse/0, colmap/sparse/0, sparse, colmap, .}.
    std::string recon_dir;

    std::string image_dir  = "images";
    std::string mask_dir   = "masks";
    std::string depth_dir  = "depths";
    std::string normal_dir = "normals";

    // When false, frames whose image file is missing on disk are kept (with
    // their computed path) instead of being skipped/rejected. The trainer
    // needs the pixels (default true); the standalone viewer only needs the
    // camera poses, so it parses cameras from a dataset whose images were not
    // provided (e.g. dropping just the sparse/ folder).
    bool require_image_files = true;

    // Fraction of training images held out for validation (linspace-spread,
    // matching get_train_eval_split_fraction). 0 = no validation set.
    float validation_fraction = 0.0f;

    // Train/eval split (dataparser.py eval_mode). The CLI has no eval pass
    // yet, so eval images are simply excluded from the parsed dataset:
    //   "all"      -> every image trains (default).
    //   "fraction" -> linspace-spread ceil(N * train_split_fraction) train.
    //   "interval" -> index % eval_interval == 0 is eval, rest train.
    //   "filename" -> basename contains "train" / "eval".
    std::string eval_mode = "all";
    float       train_split_fraction = 0.9f;
    int         eval_interval = 8;

    // Reject frames whose camera position is more than this many MADs from
    // the geometric median of all camera positions (dataparser.py:319-325).
    // inf = off (default).
    float outlier_threshold = std::numeric_limits<float>::infinity();

    // Divide stored intrinsics by this factor (Mip-NeRF 360 images_2/_4
    // style). 0 = off. TODO: the auto-detect (bool) mode which probes the
    // first image's actual resolution (dataparser.py:363-372).
    float       rescale_camera_to_fit = 0.0f;
    std::string downscale_rounding_mode = "floor";   // floor | ceil | round

    // Metashape inputs (parse_metashape_dataset). Relative paths resolve
    // against the dataset dir; empty = auto-detect a unique candidate in the
    // dataset dir (.psx is optional -- used only to disambiguate camera ->
    // image filename matching).
    std::string metashape_xml;
    std::string metashape_ply;
    std::string metashape_psx;

    // Which <component> group to use when a Metashape export contains several.
    // -1 (default) keeps the historical behavior: train on the largest group.
    // >= 0 selects that component-group index (document order); the viewer uses
    // this to let the user pick a component.
    int metashape_component = -1;
};

// Everything main.cpp needs, PER-INPUT camera (length N), sorted by image
// filename (matching dataparser.py's argsort over fnames).
struct ParsedDataset {
    int64_t num_cameras = 0;

    // CameraModelType as int (Camera.h camera_model_from_name over the
    // nerfstudio-normalized model name).
    std::vector<int32_t>     camera_models;

    // Absolute file paths. mask/depth/normal are all-N with "" for images
    // that have no such auxiliary file, or empty vectors when no image has
    // one (the DataManager convention).
    std::vector<std::string> image_filenames;
    std::vector<std::string> mask_filenames;
    std::vector<std::string> depth_filenames;
    std::vector<std::string> normal_filenames;

    std::vector<int32_t>     widths;
    std::vector<int32_t>     heights;

    // Camera-to-world, [N, 3, 4] flat, OpenGL/nerfstudio convention.
    std::vector<float>       c2w;

    std::vector<float>       intrins;      // [N, 4]  (fx, fy, cx, cy)
    std::vector<float>       dist_coeffs;  // [N, 10] (k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2)

    // Per-INPUT train/val partition (validation_fraction).
    std::vector<int32_t>     train_indices;
    std::vector<int32_t>     val_indices;

    // Seed point cloud in the training frame.
    ColmapPoints3D           points;

    // 1 / scale_factor of the would-be normalized frame (dataparser.py
    // train_frame="points" branch). Computed over ALL frames (before the
    // eval_mode subset is dropped), like the Python dataparser.
    float                    train_frame_scale = 1.0f;

    // Similarity mapping a normalized-frame point into the training frame
    // (dataparser.py:493-537 "train_to_normalized_transform" -- the name is
    // historical; the stored value is inv(T_n_from_train)). The viewer
    // client navigates in the normalized frame; Trainer._render remaps its
    // c2w through this before rendering (trainer.py:557-576). Row-major
    // 4x4; identity when train_frame_scale == 1.
    std::array<float, 16>    train_to_normalized{1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
};

ParsedDataset parse_colmap_dataset(const std::string& dataset_dir,
                                   const DatasetParserConfig& cfg);
ParsedDataset parse_nerfstudio_dataset(const std::string& dataset_dir,
                                       const DatasetParserConfig& cfg);
ParsedDataset parse_metashape_dataset(const std::string& dataset_dir,
                                      const DatasetParserConfig& cfg);

// Nerfstudio back-end over an already-built transforms.json-shaped meta
// (NerfstudioParser.cpp; used by the Metashape front-end).
struct JsonValue;
ParsedDataset parse_nerfstudio_meta(const JsonValue& meta,
                                    const std::string& dataset_dir,
                                    const DatasetParserConfig& cfg);

// Auto-detect (format = "") or dispatch ("colmap" / "nerfstudio" /
// "metashape").
ParsedDataset parse_dataset(const std::string& dataset_dir,
                            const DatasetParserConfig& cfg,
                            const std::string& format);


// ===========================================================================
// POST-split camera bake (trainer.py _setup_cpp_data_manager:257-414)
// ===========================================================================

// Per-input split factor K:
//   FISHEYE / EQUISOLID -> 5 cubemap faces when warp_to_pinhole
//   EQUIRECTANGULAR     -> 6 faces when warp_spherical_to_pinhole,
//                          else K=1 direct-equirect pass-through
//   everything else     -> 1
// K>1 sub-cameras are canonical PINHOLE at out_shape = ceil(sqrt(H*W/K)),
// fx=fy=cx=cy=out_shape/2, zero distortion, cubemap-face rotations matching
// DataManager.cpp's kAxesFisheye5 / kAxesEquirect6 tables.
struct PostSplitCameras {
    int64_t n_post = 0;
    bool any_warp        = false;   // any K > 1
    bool any_fisheye_warp = false;  // fisheye/equisolid split present
                                    // (DataManager synthesizes FOV masks)
    bool direct_equirect = false;   // un-split equirect present

    std::vector<int32_t> K_per_camera;   // [N]  (empty vectors passed to the
    std::vector<int32_t> post_offsets;   // [N]   engine when !any_warp)

    // Engine-convention world-to-camera (R/T inverse + Y/Z flip pre-baked,
    // trainer.py:403-414), POST-split.
    std::vector<float>   viewmats;       // [N_post, 4, 4]
    std::vector<float>   intrins;        // [N_post, 4]
    std::vector<float>   dist_coeffs;    // [N_post, 10]

    // Viewer (engine_viewer_init) arrays, POST-split: camera-to-world in the
    // y/z-flipped form the blit kernel expects (annotation.py:103-106), plus
    // per-post W/H/model.
    std::vector<float>   c2w_flip;       // [N_post, 3, 4]
    std::vector<int32_t> post_widths;    // [N_post]
    std::vector<int32_t> post_heights;   // [N_post]
    std::vector<int32_t> post_models;    // [N_post] CameraModelType (faces = PINHOLE)

    // Per-INPUT copies for the wide-warp kernel; empty when !any_warp.
    std::vector<float>   input_intrins;      // [N, 4]
    std::vector<float>   input_dist_coeffs;  // [N, 10]
};

PostSplitCameras bake_post_split(const ParsedDataset& ds,
                                 bool warp_to_pinhole,
                                 bool warp_spherical_to_pinhole);


// ===========================================================================
// Shared parser internals (DatasetCommon.cpp) -- used by both format
// parsers; exposed here so they stay in one place.
// ===========================================================================
namespace dsparse {

// Normalized-frame scale factor over c2w [N,3,4] (orient="up",
// center="poses", auto-scale). Only the scalar matters for
// train_frame="points". Returns 1/max_abs (the dataparser's scale_factor).
double compute_normalized_scale_factor(const std::vector<float>& c2w, int64_t n);

// Full normalized-frame similarity: writes the row-major 4x4
// T_n_from_camera = scale * [R_align | -R_align @ center] (dataparser.py
// T_n_from_camera, lines 497-499) and returns scale_factor. The viewer
// remap transform is inv(T_n_from_camera @ applied_transform).
double compute_normalized_transform(const std::vector<float>& c2w, int64_t n,
                                    double T_out[16]);

// inv([A|b; 0 1]) for a general invertible 3x3 A (row-major 4x4 in/out).
void invert_affine4x4(const double in[16], double out[16]);

// eval_mode train subset over N sorted frames; identity for "all".
// `names` are image filenames (used by eval_mode="filename").
std::vector<int64_t> train_subset(int64_t n, const std::vector<std::string>& names,
                                  const DatasetParserConfig& cfg);

// validation_fraction partition of 0..N-1 into ds.train_indices/val_indices.
void assign_val_split(ParsedDataset& ds, float validation_fraction);

// Auxiliary mask/depth/normal discovery by filename convention
// (dataparser.py _add_auxiliary_buffers, trimmed candidate list).
std::string find_aux_file(const std::string& aux_dir, const std::string& rel_name,
                          const char* suffix_tag);

// Outlier-frame mask via geometric median of camera positions
// (dataparser.py:138-173, 319-325). Returns keep-flags, all-true when
// threshold is inf. positions = [N, 3].
std::vector<char> outlier_keep_mask(const std::vector<double>& positions,
                                    int64_t n, float threshold);

}  // namespace dsparse
