#pragma once

// What a dataset run is doing, in the shape the screen draws it.
//
// One object per runner, written by its worker thread and read by the UI
// thread. It replaces a stage string, a progress float and an unclassified
// stream of log lines with the three things the screen actually asks: which
// step is running, how far through it is, and which lines a user who is not
// debugging should be shown.
//
// SfmRunner fills it by parsing the child's stdout; when the SfM module
// becomes a library (docs/notes/sfm-port-plan.md phase 3) only that parser goes.

#include <cstdint>
#include <mutex>
#include <string>
#include <vector>

namespace gui {

// The steps a dataset run goes through, in order. Both engines report through
// the same list; a step a run does not need is never entered.
enum class Stage {
    Frames, Masks, Features, Matching, Mapping, Finishing
};
inline constexpr int kNumStages = 6;

enum class StageStatus { Pending, Running, Done, Skipped, Failed };

struct StageProgress {
    StageStatus status = StageStatus::Pending;
    float   fraction = -1.0f;    // 0..1, or -1 when the step cannot estimate one
    int64_t done = 0, total = 0;
    std::string detail;          // already formatted; the line under the bar
};

// A line of a run's output. `detail` is the stream a developer reads; the rest
// is the handful worth showing without asking for it.
struct RunLine {
    Stage stage = Stage::Frames;
    bool  detail = true;
    std::string text;
};

class RunProgress {
public:
    void reset();

    // Enter a step. Everything before it still Pending becomes Skipped, so a
    // run that masks nothing does not leave "Masks" looking like it is next.
    void enter(Stage s, const std::string& detail);
    void count(Stage s, int64_t done, int64_t total);
    void fraction(Stage s, float f);
    void detail(Stage s, const std::string& text);
    void mark(Stage s, StageStatus st);
    // The run ended: whatever was running takes `st`.
    void finish(StageStatus st);

    void note(const std::string& text, bool detail);
    void note(Stage s, const std::string& text, bool detail);

    std::vector<RunLine> drain();
    StageProgress stage(Stage s) const;
    Stage current() const;
    // The step a job would enter next if it started now -- what decides which
    // half of the form is still editable.
    bool ran(Stage s) const;

private:
    mutable std::mutex _mu;
    StageProgress _st[kNumStages];
    Stage _cur = Stage::Frames;
    std::vector<RunLine> _pending;
};

}  // namespace gui
