#pragma once

// DataManager — host-side training-data orchestrator.
//
// C++ port of `modules/datamanager.py` + `modules/dataset.py`. Owns the on-disk
// dataset layout, decodes images via stb_image, batches indices by image-shape
// group (so each engine_train_step sees uniform [B,H,W,C] inputs), and feeds
// the resulting batches directly into the engine.
//
// Two caching modes are supported:
//
//   CacheMode::CPU
//       At construction time, every image (RGB + mask + depth + normal as
//       configured) is decoded into a per-image host-pageable byte buffer.
//       next_*_batch() then just gathers the relevant rows into a contiguous
//       batch buffer and hands the engine a host TorchTensorView. Cheap
//       per-step latency, peak RAM ~= the full decoded dataset.
//
//   CacheMode::DISK
//       Nothing is pre-decoded. A prefetch pipeline runs three (or four, when
//       masks are present) independent worker pools — RGB+MASK / DEPTH /
//       NORMAL — each with its own bounded job queue and a bounded ready queue.
//       The pipeline keeps ``prefetch_batches`` decoded batches in flight at a
//       time; next_*_batch() pops the next ready one and blocks only if the
//       pool is empty. This mirrors PyTorch's DataLoader's prefetch_factor +
//       num_workers semantics, with the per-modality split that the user asked
//       for (depth/normal often dominate decode cost on dense-supervision
//       datasets and benefit from dedicated worker threads).
//
// Image format support (stb_image): 8-bit and 16-bit PNG (one channel for
// depth / masks, three or four for RGB / normals), JPEG, BMP, TGA. EXR / DNG
// from the Python path are intentionally dropped — bring them back via a
// pluggable decoder if needed.

#include "Common.cuh"
#include "Tensor.h"

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>


// ===========================================================================
// Public configuration
// ===========================================================================

enum class CacheMode {
    CPU,    // pageable host RAM, all images pre-decoded at construction
    DISK    // async per-modality prefetch pool, nothing pre-decoded
};

struct DataManagerConfig {
    CacheMode cache_mode = CacheMode::CPU;

    bool load_masks   = true;
    bool load_depths  = true;
    bool load_normals = true;

    // Train batch size. Picked per-group at next_train_batch() time — groups
    // with fewer than batch_size images yield a smaller batch.
    int  train_batch_size = 1;

    // Eval batch size. Ignored when val_indices is empty.
    int  val_batch_size   = 1;

    // RNG seed for group sampling and shuffle. 0 lets the implementation pick
    // a nondeterministic value (steady_clock).
    uint64_t seed = 0;

    // Worker count per modality pool (DISK mode only). Defaults to a small
    // fixed number; tune up for fast disks / many CPU cores.
    int workers_rgb    = 16;
    int workers_depth  = 8;
    int workers_normal = 8;

    // How many batches to keep in flight in the prefetch pipeline (DISK only).
    // Hard upper bound on RAM consumed by ready / partially-ready batches.
    int prefetch_batches = 4;

    // Signed boundary offset applied to binarized masks at decode time,
    // expressed as a fraction of sqrt(W*H) of the decoded mask. Positive ->
    // dilate (grow) foreground; negative -> erode (shrink) foreground; zero
    // disables. Implemented via separable Felzenszwalb-Huttenlocher squared
    // Euclidean distance transform on CPU (O(N) per row+col, exact).
    float mask_boundary_offset = 0.0f;

    // When true: cameras with model FISHEYE or EQUISOLID get split into 5
    // cubemap-face pinhole sub-cameras at training time. PINHOLE cameras
    // pass through unchanged. EQUIRECTANGULAR cameras are ALWAYS split into
    // 6 sub-cameras regardless of this flag. The byte->float conversion of
    // the GT image is fused with the warp on GPU, so no full-res
    // intermediate float buffer is allocated. Depth, normal, and PPISP are
    // not supported when any image in the dataset is split.
    bool warp_to_pinhole = false;
};


// ===========================================================================
// Decoded batch
// ===========================================================================

// Element-type tag, matching the engine's set_training_data inputs.
//   RGB:    UINT8 (sRGB 0..255), UINT16 (0..65535), or FLOAT32 (0..1).
//   MASK:   UINT8 (treated as bool; nonzero -> valid pixel).
//   DEPTH:  UINT16 (raw counts, no scaling) or FLOAT32.
//   NORMAL: UINT8 (x/127.5 - 1) or FLOAT32 (already in [-1, 1]).
enum class PixelDType : uint8_t { UINT8 = 1, UINT16 = 2, FLOAT32 = 4 };

inline uint32_t pixel_dtype_size(PixelDType t) { return (uint32_t)t; }


// A single ready-to-consume training batch.
//
// `*_buffer` are host buffers (std::vector for ownership + automatic growth);
// `*_view` is a non-owning TorchTensorView pointing into the corresponding
// buffer in the [B, H, W, C] (or [B, H, W, 1]) channel-last layout expected
// by the engine.
//
// Empty optional modalities (no mask / no depth / no normal) leave the
// corresponding view at all-zero (data_ptr == 0), which the engine interprets
// as "not provided".
struct DecodedBatch {
    // Per-batch metadata. When warp_to_pinhole splits a batch's images, these
    // fields describe the POST-split state -- so `width`, `height`, `num`,
    // `model`, `viewmats`, `intrins`, `dist_coeffs` are what the engine's
    // set_camera_params consumes directly. The INPUT (unwarped) shape is
    // tracked separately below so the GT warp kernel can read the byte
    // staging buffer.
    int32_t                width  = 0;     // post-split W (= input W when K=1)
    int32_t                height = 0;     // post-split H
    int32_t                num    = 0;     // post-split B (= input B * K)
    CameraModelType        model  = (CameraModelType)-1;   // PINHOLE when K>1

    // Selected camera indices, in input batch order (length = input B).
    std::vector<int32_t>   indices;

    // Camera parameters at POST-split resolution, host-contiguous.
    // Length = 16 / 4 / 10 * num (= B*K).
    std::vector<float>     viewmats;
    std::vector<float>     intrins;
    std::vector<float>     dist_coeffs;

    // Per-camera post-split offset for bilagrid / per-image lookup. Length
    // = input B. bilagrid_cam_indices for the batch is built from this:
    // [post_offsets[j] + k for j in 0..B-1 for k in 0..K-1].
    std::vector<int32_t>   post_offsets;

    // Pixel buffers + their dtypes. `rgb_dtype` is the original on-disk dtype
    // (so 16-bit PNGs survive into the upload as uint16); the engine's GPU
    // converter handles the cast.
    std::vector<uint8_t>   rgb_buffer;
    PixelDType             rgb_dtype = PixelDType::UINT8;

    std::vector<uint8_t>   mask_buffer;    // uint8, empty if disabled / none
    std::vector<uint8_t>   depth_buffer;
    PixelDType             depth_dtype = PixelDType::UINT16;

    std::vector<uint8_t>   normal_buffer;
    PixelDType             normal_dtype = PixelDType::UINT8;

    // Per-modality (H, W). Mask / depth / normal images are allowed to be a
    // different resolution from the rendered RGB output — the engine's
    // fused per-pixel loss kernel bilinearly samples from these buffers,
    // so the data manager just exposes the raw on-disk shape and leaves the
    // resize to the loss path. Zero when the modality is disabled / absent.
    int32_t                mask_height   = 0;
    int32_t                mask_width    = 0;
    int32_t                depth_height  = 0;
    int32_t                depth_width   = 0;
    int32_t                normal_height = 0;
    int32_t                normal_width  = 0;

    // ---- Warp metadata (only meaningful when K > 1) -----------------------
    // `input_*` describe the unwarped GT image shape and the input batch;
    // the engine's warp kernel reads `rgb_buffer` / `mask_buffer` at those
    // sizes and writes the warped result at post-split resolution
    // (width / height / num). `K` is the split factor (1 = pass-through,
    // 5 = fisheye, 6 = equirectangular). `input_model` selects the warp
    // kernel (fisheye / equisolid / equirectangular). `axes_dev` is a
    // device pointer to a [K, 3, 3] cubemap axis table -- lifetime is
    // managed by the DataManager.
    int32_t                input_width   = 0;
    int32_t                input_height  = 0;
    int32_t                input_num     = 0;   // = num / K
    int32_t                K             = 1;
    CameraModelType        input_model   = (CameraModelType)-1;
    const float*           axes_dev      = nullptr;

    // Per-INPUT intrins / dist_coeffs (length input_num). Populated when
    // K > 1 and the wide warp kernel (fisheye/equisolid) needs them. The
    // equirectangular path ignores these.
    std::vector<float>     input_intrins;       // 4 * input_num
    std::vector<float>     input_dist_coeffs;   // 10 * input_num
    TorchTensorView        input_intrins_view{0, 0, {}};
    TorchTensorView        input_dist_coeffs_view{0, 0, {}};

    // Engine-facing TorchTensorViews. Lazily filled by `build_views()` once
    // the buffers are populated. The Engine consumes these directly via
    // set_camera_params / set_training_data.
    TorchTensorView        viewmats_view{0, 0, {}};
    TorchTensorView        intrins_view{0, 0, {}};
    TorchTensorView        dist_coeffs_view{0, 0, {}};
    TorchTensorView        rgb_view{0, 0, {}};
    TorchTensorView        mask_view{0, 0, {}};
    TorchTensorView        depth_view{0, 0, {}};
    TorchTensorView        normal_view{0, 0, {}};

    // Populate the *_view fields from the *_buffer fields. Idempotent.
    void build_views();
};


// ===========================================================================
// Internal: per-modality prefetch worker pool (DISK mode)
//
// Forward-declared so DataManager can hold a unique_ptr<Impl>.
// ===========================================================================
class DataManagerImpl;


// ===========================================================================
// DataManager
// ===========================================================================
class DataManager {
public:

    // The constructor takes the parsed-out dataset arrays directly. The
    // dataparser port populates these — DataManager doesn't know anything
    // about transforms.json / COLMAP.
    //
    // All filename vectors are sized to the total dataset; `mask_filenames`
    // / `depth_filenames` / `normal_filenames` may be empty to indicate
    // "no modality of that kind on disk". An individual entry that is the
    // empty string indicates "this image has no mask/depth/normal" and the
    // corresponding batch row is left zero.
    //
    // `train_indices` / `val_indices` partition the dataset. They are NOT
    // required to be disjoint; that's the caller's concern.
    DataManager(
        DataManagerConfig          config,
        // Per-INPUT-camera (length N) ---------------------------------------
        std::vector<int32_t>       camera_models,
        std::vector<std::string>   image_filenames,
        std::vector<std::string>   mask_filenames,
        std::vector<std::string>   depth_filenames,
        std::vector<std::string>   normal_filenames,
        std::vector<int32_t>       widths,
        std::vector<int32_t>       heights,
        // Per-INPUT split factor and post-split offset. K[i] = 1 / 5 / 6
        // depending on the camera model + warp config. Pre-computed on the
        // caller side (typically the Python trainer) so the C++ side only
        // needs to read it. offsets[i] is the starting index in the
        // POST-split camera arrays for input camera i, computed as
        // exclusive-prefix-sum of K. Pass empty vectors (or all-ones K +
        // identity offsets) to disable the warp path.
        std::vector<int32_t>       K_per_camera,        // length N
        std::vector<int32_t>       post_offsets,        // length N
        // Per-POST-split-camera (length N_post = sum(K_per_camera)) -----------
        // When the warp path is inactive, N_post == N and these are just
        // the per-input camera params. When the warp path is active, the
        // caller pre-expands these per cubemap face / equirectangular face
        // using the same algebra as modules/resample.py.
        std::vector<float>         viewmats,            // [N_post, 4, 4]
        std::vector<float>         intrins,             // [N_post, 4]
        std::vector<float>         dist_coeffs,         // [N_post, 10]
        // Per-INPUT intrins/dist_coeffs (length N). Needed by the wide warp
        // kernel (fisheye/equisolid projection). Pass empty vectors when
        // no fisheye-warp camera is in the dataset -- the kernel only reads
        // them on the fisheye path.
        std::vector<float>         input_intrins,       // [N, 4]
        std::vector<float>         input_dist_coeffs,   // [N, 10]
        std::vector<int32_t>       train_indices,
        std::vector<int32_t>       val_indices);

    ~DataManager();

    DataManager(const DataManager&)            = delete;
    DataManager& operator=(const DataManager&) = delete;

    // ---- Batch fetch ------------------------------------------------------

    // Block until the next training batch is ready and return a reference to
    // it. The returned reference is valid only until the next call to
    // `next_train_batch()`.
    const DecodedBatch& next_train_batch();

    // Same for eval. Returns nullptr if no validation indices were configured.
    const DecodedBatch* next_val_batch();

    // ---- Stats ------------------------------------------------------------

    int64_t num_train()     const;
    int64_t num_val()       const;
    bool    has_val()       const;
    CacheMode cache_mode()  const;

    // Whether each modality is present in the dataset AND enabled.
    bool has_masks()        const;
    bool has_depths()       const;
    bool has_normals()      const;

private:
    std::unique_ptr<DataManagerImpl> _impl;
};
