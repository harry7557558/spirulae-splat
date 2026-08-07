#pragma once

// Batch training: a list of (dataset, preset, output folder) rows the trainer
// works through one after another, with nobody watching.
//
// The list is data, not a runner. Driving it is a handful of lines in
// GuiApp::frame() -- launch a row, wait for TrainRunner to leave Training,
// record what happened, launch the next -- because that reuses the whole live
// trainer screen: the viewport still shows the run in flight, the metrics
// still plot, the log still scrolls. A second training driver here would be a
// second implementation of the thing app/TrainerCore.h already is.
//
// Two properties an unattended queue has to have, and where they live:
//
//   * a row that fails does not stop the rest. TrainRunner turns any
//     exception into Phase::TrainError, so a failing row is an ordinary
//     transition, not a crash -- BatchJob::status records it and the queue
//     moves on.
//   * problems are reported BEFORE anything starts. batch_check() is the
//     pre-flight: it answers, for every row at once, "is there a reason this
//     is going to fail" -- a dataset folder with no reconstruction in it, a
//     preset file that has gone missing, an option the trainer does not
//     implement, two rows writing to the same place. Fatal issues block the
//     start; the rest are said out loud and then ignored.

#include "config/TrainConfig.h"
#include "i18n/Message.h"

#include <string>
#include <vector>

namespace gui {

// One thing a pre-flight found. `text` is interface copy with an optional
// {0}; `raw` is engine text (a flag name, a parser message) that stays in
// English -- exactly the ui:: / ui::*Raw split. Exactly one of them is set.
struct BatchIssue {
    const spirula::i18n::Msg* text = nullptr;
    std::string arg;      // fills {0} in `text`
    std::string raw;      // the whole line, untranslated
    bool fatal = false;   // blocks the start; otherwise a warning
};

// The line to draw. Already translated and formatted, so it goes on screen
// through a ui::*Raw call.
std::string batch_issue_line(const BatchIssue& issue);

struct BatchJob {
    std::string dataset;
    // "" means the built-in preset `preset_name`; otherwise the preset file at
    // this path, and `preset_name` is what it calls itself (for display).
    std::string preset_path;
    std::string preset_name = "3dgs";
    // Where runs go. "" = <dataset>/outputs, the same default the trainer
    // screen picks when a dataset is opened. Each run still lands in its own
    // timestamped subfolder of this, so two rows sharing one folder is fine.
    std::string output_dir;

    // Per-row overrides of the three flags that get changed often enough that
    // making a whole preset for each combination is the wrong shape of work.
    // "" means "whatever the preset says".
    //
    // Held as TEXT, not as ints: "unset" and "0" are different answers and an
    // integer box cannot show the difference -- 0 is a legal --sh-degree. The
    // pre-flight reports anything that is not a usable number, so a typo is
    // caught before the queue starts rather than silently rounded to a
    // default.
    std::string cap_max_override;
    std::string sh_degree_override;
    std::string iterations_override;

    enum class Status { Pending, Running, Done, Failed, Skipped, Stopped };
    Status status = Status::Pending;

    // Filled as the row runs. `message` is engine text (English).
    std::string message;
    std::string out_dir;      // where it actually wrote
    int steps = 0;
    double seconds = 0.0;

    // From the last batch_check(); empty until one has run.
    std::vector<BatchIssue> issues;
};

// True if any issue found on this row blocks the start.
bool batch_has_error(const BatchJob& job);

// Check one row. `all`/`index` are for the checks that are about the list
// rather than the row -- another row that would do exactly the same work.
// Does no GPU work and touches nothing on disk.
std::vector<BatchIssue> batch_check(const BatchJob& job,
                                    const std::vector<BatchJob>& all,
                                    int index);

// The config this row would train with: the preset, then the row's dataset and
// output folder, then the macro options resolved the way the trainer screen
// resolves them. `preset_base` comes back as the built-in preset name to
// record in the run's config.json. Returns false and fills `error` (English)
// when the preset cannot be read -- the same condition batch_check() reports,
// re-checked here because the file can go away in between.
bool batch_build_config(const BatchJob& job, TrainConfig& cfg,
                        std::string& preset_base, std::string& error);

// The list, kept across sessions in <config_dir>/batch.json. A queue that is
// worth setting up is worth surviving a restart, and losing it to a crash
// three hours into a five-dataset run is the one failure the user cannot
// recover from by trying again. Neither call throws; a list that cannot be
// read comes back empty.
std::vector<BatchJob> load_batch_list();
void save_batch_list(const std::vector<BatchJob>& jobs);

}  // namespace gui
