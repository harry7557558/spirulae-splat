#pragma once

// ColmapParser -- standalone C++ COLMAP reconstruction reader + dataset baker.
//
// C++ port of the COLMAP path through modules/colmap_utils.py +
// modules/dataparser.py (_parse_colmap_data -> _parse_nerfstudio_data), for
// the standalone CLI trainer (app/main.cpp). Produces the flat per-camera
// arrays that engine_setup_data_manager / engine_viewer_init consume
// directly, plus the seed point cloud for splat initialization.
//
// Scope (deliberately narrower than the Python dataparser):
//   - COLMAP binary format only (cameras.bin / images.bin / points3D.bin);
//     text fallback is a TODO. transforms.json / Metashape stay Python-only
//     until this parser grows a JSON dependency.
//   - train_frame = "points" semantics (the current default): poses and
//     points stay in the raw COLMAP-derived frame; the normalized-frame
//     similarity is computed only to derive `train_frame_scale`, which the
//     trainer uses to rescale means_lr / max_world_size / noise_lr.
//   - eval_mode = "all" + validation_fraction split. fraction / filename /
//     interval modes are a TODO (port dataparser.py:51-129).
//   - No warp_to_pinhole camera pre-splitting (K = 1 per camera). The
//     dataset arrays are laid out so adding the split later only means
//     filling K_per_camera / post_offsets and expanding the camera arrays
//     (same algebra as trainer.py _setup_cpp_data_manager:324-401).

#include <array>
#include <cstdint>
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

// Low-level readers. `recon_dir` is the directory holding cameras.bin etc.
// Throw std::runtime_error on missing/corrupt files.
std::map<int32_t, ColmapCamera> read_cameras_binary(const std::string& recon_dir);
std::map<int32_t, ColmapImage>  read_images_binary(const std::string& recon_dir);
ColmapPoints3D                  read_points3D_binary(const std::string& recon_dir);


// ===========================================================================
// Baked dataset
// ===========================================================================

struct ColmapParserConfig {
    // Reconstruction dir relative to the dataset dir. Empty = auto-detect
    // over {sparse/0, colmap/sparse/0, sparse, colmap, .} like
    // dataparser.py:_parse_colmap_data.
    std::string recon_dir;

    std::string image_dir  = "images";
    std::string mask_dir   = "masks";
    std::string depth_dir  = "depths";
    std::string normal_dir = "normals";

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

    // Divide stored intrinsics by this factor (Mip-NeRF 360 images_2/_4
    // style). 0 = off. TODO: port the auto-detect (bool) mode which probes
    // the first image's actual resolution (dataparser.py:363-372).
    float rescale_camera_to_fit = 0.0f;
};

// Everything the CLI trainer needs, in engine-native layout. All per-camera
// arrays are length N (= number of registered COLMAP images), sorted by
// image filename (matching dataparser.py's argsort over fnames).
struct ParsedDataset {
    int64_t num_cameras = 0;

    // CameraModelType as int (Camera.h camera_model_from_name over the
    // nerfstudio-normalized model name, e.g. OPENCV -> PINHOLE).
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

    // Camera-to-world, [N, 3, 4] flat, OpenGL/nerfstudio convention
    // (COLMAP w2c inverted, then Y/Z camera-axis flip). This is what the
    // viewer (engine_viewer_init) consumes.
    std::vector<float>       c2w;

    // World-to-camera, [N, 4, 4] flat, engine/OpenCV convention -- the
    // Y/Z flip + R/T inverse of `c2w`, pre-baked exactly like
    // trainer.py _setup_cpp_data_manager:403-414. This is what
    // engine_setup_data_manager consumes.
    std::vector<float>       viewmats;

    std::vector<float>       intrins;      // [N, 4]  (fx, fy, cx, cy)
    std::vector<float>       dist_coeffs;  // [N, 10] (k1 k2 k3 k4 p1 p2 sx1 sy1 b1 b2)

    std::vector<int32_t>     train_indices;
    std::vector<int32_t>     val_indices;

    // Seed point cloud in the training frame.
    ColmapPoints3D           points;

    // 1 / scale_factor of the would-be normalized frame (dataparser.py
    // train_frame="points" branch). The trainer multiplies means_lr /
    // max_world_size / noise_lr by this so tuning transfers across scenes
    // of different metric scale.
    float                    train_frame_scale = 1.0f;
};

// Parse + bake. Throws std::runtime_error with a user-facing message on
// failure (missing recon dir, missing images, unsupported camera model).
ParsedDataset parse_colmap_dataset(const std::string& dataset_dir,
                                   const ColmapParserConfig& cfg);
