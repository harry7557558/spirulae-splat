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
    int workers_rgb    = 8;
    int workers_depth  = 8;
    int workers_normal = 8;

    // How many batches to keep in flight in the prefetch pipeline (DISK only).
    // Hard upper bound on RAM consumed by ready / partially-ready batches.
    int prefetch_batches = 4;
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
    // Per-batch metadata.
    int32_t                width  = 0;
    int32_t                height = 0;
    int32_t                num    = 0;     // B
    CameraModelType        model  = (CameraModelType)-1;

    // Selected camera indices, in batch order.
    std::vector<int32_t>   indices;

    // Camera parameters, host-contiguous [B,4,4] / [B,4] / [B,10] floats.
    std::vector<float>     viewmats;       // 16 * B
    std::vector<float>     intrins;        // 4 * B
    std::vector<float>     dist_coeffs;    // 10 * B

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
        CameraModelType            model,
        std::vector<std::string>   image_filenames,
        std::vector<std::string>   mask_filenames,
        std::vector<std::string>   depth_filenames,
        std::vector<std::string>   normal_filenames,
        std::vector<int32_t>       widths,
        std::vector<int32_t>       heights,
        std::vector<float>         viewmats,      // [N, 4, 4] row-major
        std::vector<float>         intrins,       // [N, 4]
        std::vector<float>         dist_coeffs,   // [N, 10]
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
