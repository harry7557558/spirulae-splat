#pragma once

// ConfigUI -- the "Advanced options" editor. All widgets are expanded from
// the SS_CONFIG_FIELDS X-macro (config/TrainConfig.h), so every one of
// the training config fields is editable, with:
//   - grouping under the field table's `section` headings, collapsed by
//     default so novices are not overwhelmed
//   - a detail filter over the field table's `tier`, so the default view is
//     the ~20 flags a first run needs rather than all 178
//   - a search box filtering by flag name and help text, which ignores the
//     detail filter: searching finds a flag wherever it is buried
//   - the field's help text as a hover tooltip (+ the preset default)
//   - modified-from-preset highlighting and right-click "Reset to default"
//   - fields with `choices` as dropdowns, std::optional as auto/override
//
// A row added to the field table appears here automatically -- no GUI change
// needed.

#include "config/TrainConfig.h"

#include <set>
#include <string>

namespace gui {

struct ConfigUIState {
    char search[128] = "";
    bool modified_only = false;
    int  tier = 0;   // rank into kTrainTiers: show this specialist and below

    // Flags the user edited here, by name. Macro options (--quality and
    // friends) leave these alone, so touching a flag by hand takes it out of
    // a macro's reach for good -- see train_resolve_macros(). Reset when the
    // whole config is, i.e. on a preset change.
    std::set<std::string> touched;
};

// Draw the full generated editor. `defaults` is the preset-applied baseline
// used for modified-highlighting and reset. Returns true when any field
// changed this frame; the caller re-resolves the macro options when it does.
bool draw_config_editor(TrainConfig& cfg, const TrainConfig& defaults,
                        ConfigUIState& st);

}  // namespace gui
