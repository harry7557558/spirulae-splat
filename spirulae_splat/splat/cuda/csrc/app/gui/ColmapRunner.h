#pragma once

// ColmapRunner -- turns raw images or a video into a trainable COLMAP
// dataset by driving the external `colmap` (>= 4.x; the CLI flags follow
// scripts/run_colmap.bash) and, for video, `ffmpeg`, on a worker thread with
// live log streaming and cancellation.
//
// Pipeline:
//   [video] ffmpeg frame extraction (oversampled) -> sharpest-frame
//           selection (FrameSelect, extract_frames.py port); multi-track
//           .insv files split into images/cam<N>/ with one camera per folder
//   [optional] AI masking via the embedded scripts/mask.py (external Python
//           with lang-segment-anything; masks feed COLMAP and the trainer)
//   feature_extractor -> exhaustive / sequential / vocab-tree matcher (the
//           vocab tree is auto-found or downloaded) -> mapper ->
//           [optional] bundle_adjuster refinement
//
// Output layout (what the dataset parsers auto-detect):
//   <workspace>/database.db
//   <workspace>/sparse/0/{cameras,images,points3D}.bin
//   <workspace>/images/            (video input: extracted frames)
//   <workspace>/masks/             (when masking is enabled)
//   <workspace>/ssplat_dataset.json  (records the image dir when the source
//                                     images are referenced in place)
//
// For a folder-of-images input the images are NOT copied; COLMAP indexes
// them where they are (recursively) and the GUI passes the absolute path as
// image_dir (the parsers join dataset_dir / image_dir, absolute wins).

#include <atomic>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

// COLMAP camera models the dataset parser understands (ColmapParser.cpp).
inline const char* kColmapCameraModels[] = {
    "OPENCV", "PINHOLE", "SIMPLE_PINHOLE", "SIMPLE_RADIAL", "RADIAL",
    "FULL_OPENCV", "OPENCV_FISHEYE", "THIN_PRISM_FISHEYE",
    "SIMPLE_RADIAL_FISHEYE", "RADIAL_FISHEYE",
};
inline constexpr int kNumColmapCameraModels = 10;

struct ColmapJob {
    std::string input_path;              // images folder, or a video file
    bool is_video = false;
    std::string workspace;               // output dataset dir (created)
    std::string colmap_exe = "colmap";
    std::string ffmpeg_exe = "ffmpeg";
    std::string python_exe = "python3";  // for the masking script

    // Cameras
    std::string camera_model = "OPENCV"; // ImageReader.camera_model
    int camera_mode = 0;                 // 0 = one shared camera,
                                         // 1 = one per subfolder,
                                         // 2 = one per image

    // Matching
    int quality = 1;                     // 0 = fast, 1 = balanced, 2 = high
    int matcher = 0;                     // 0 auto, 1 exhaustive,
                                         // 2 sequential, 3 vocab tree

    // Video extraction
    float video_fps = 2.0f;              // kept frames per second
    int sharp_window = 3;                // pick sharpest of N candidates (1 = off)
    int max_frames = 100000;

    // Advanced
    int max_num_features = 0;            // 0 = per-quality default
    int max_image_size = 0;              // 0 = COLMAP default
    int seq_overlap = 10;                // sequential matcher overlap
    bool seq_quadratic_overlap = true;   // match frame i with i +- 2^k too
    bool estimate_affine_shape = false;  // + guided feature matching
    bool ba_use_gpu = true;              // Mapper.ba_use_gpu
    bool final_bundle_adjust = true;     // bundle_adjuster refinement pass
    std::string vocab_tree_path;         // "" = auto find / download

    // AI masking (embedded scripts/mask.py; needs Python + lang-sam)
    bool mask_enable = false;
    std::string mask_prompt;             // "people; cars; ..."
    std::string mask_negative_prompt;
    std::string mask_model = "sam2.1_hiera_large";
    int mask_max_image_size = 1600;
};

class ColmapRunner {
public:
    enum class State { Idle, Running, Done, Failed, Cancelled };

    ~ColmapRunner();

    void start(const ColmapJob& job);
    void cancel();

    State state() const { return _state.load(); }
    std::string stage();                 // current pipeline stage label
    std::string error();                 // set when Failed
    std::string dataset_dir();           // valid when Done
    std::string image_dir();             // image_dir to train with ("" = default)
    std::vector<std::string> drain_log();

private:
    void run(ColmapJob job);
    void log(const std::string& line);
    int  exec(const std::vector<std::string>& argv);
    void set_stage(const std::string& s);
    bool check_colmap_version(const ColmapJob& job, std::string& err);
    std::string resolve_vocab_tree(const ColmapJob& job);
    bool run_masking(const ColmapJob& job, const std::string& images,
                     std::string& err);

    std::thread _worker;
    std::atomic<State> _state{State::Idle};
    std::atomic<bool> _cancel{false};
    std::mutex _mu;                      // guards strings + log below
    std::string _stage, _error, _dataset_dir, _image_dir;
    std::vector<std::string> _log;
};

}  // namespace gui
