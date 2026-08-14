#pragma once

// FeatureWatcher -- the frames of a run, with the features found on them, while
// the extractor is still working through the rest.
//
// It walks the image list in the order the extractor does and waits for each
// one's feature file to appear, so the strip advances at the pace of the stage
// rather than at the pace of the disk. Its own thread because it decodes a
// full-size image per frame, which is not something to do between two ImGui
// calls.
//
// No channel from the child is involved: the feature files ARE the channel,
// and they are written whether or not anybody is watching.

#include "app/gui/FilmStrip.h"

#include <atomic>
#include <string>
#include <thread>

namespace gui {

class FeatureWatcher {
public:
    ~FeatureWatcher();

    // Idempotent for the same pair of folders: starting it again on the folders
    // it is already watching does nothing, which is what lets the screen call
    // it every frame the step is running.
    void start(const std::string& image_dir, const std::string& features_dir,
               FilmStrip* film);
    void stop();

private:
    void run(std::string image_dir, std::string features_dir, FilmStrip* film);

    std::thread _worker;
    std::atomic<bool> _stop{false};
    std::string _image_dir, _features_dir;
};

}  // namespace gui
