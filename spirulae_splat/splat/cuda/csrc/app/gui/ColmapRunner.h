#pragma once

// ColmapRunner -- turns raw images or a video into a trainable COLMAP
// dataset by driving the external `colmap` (and, for video, `ffmpeg`) CLIs
// on a worker thread with live log streaming and cancellation.
//
// Output layout (what the dataset parsers auto-detect):
//   <workspace>/database.db
//   <workspace>/sparse/0/{cameras,images,points3D}.bin
//   <workspace>/images/            (video input: extracted frames)
//   <workspace>/ssplat_dataset.json  (records the image dir when the source
//                                     images are referenced in place)
//
// For a folder-of-images input the images are NOT copied; COLMAP indexes
// them where they are and the GUI passes the absolute path as image_dir
// (the parsers join dataset_dir / image_dir, which std::filesystem resolves
// to the absolute path).

#include <atomic>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace gui {

struct ColmapJob {
    std::string input_path;              // images folder, or a video file
    bool is_video = false;
    std::string workspace;               // output dataset dir (created)
    std::string colmap_exe = "colmap";
    std::string ffmpeg_exe = "ffmpeg";
    std::string camera_model = "OPENCV"; // ImageReader.camera_model
    bool single_camera = true;           // share intrinsics across images
    int quality = 1;                     // 0 = fast, 1 = balanced, 2 = high
    int matcher = 0;                     // 0 = auto, 1 = exhaustive, 2 = sequential
    bool use_gpu = true;                 // SIFT extraction/matching on GPU
    float video_fps = 2.0f;              // frames per second to extract
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

    std::thread _worker;
    std::atomic<State> _state{State::Idle};
    std::atomic<bool> _cancel{false};
    std::mutex _mu;                      // guards strings + log below
    std::string _stage, _error, _dataset_dir, _image_dir;
    std::vector<std::string> _log;
};

}  // namespace gui
