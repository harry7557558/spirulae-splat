#pragma once

// ConfigUI -- the "Advanced options" editor. All widgets are expanded from
// the SSPLAT_CONFIG_FIELDS X-macro (config/TrainConfig.h), so every one of
// the training config fields is editable, with:
//   - grouping by config group (Run / Dataset / Data loading / Model /
//     Optimizer), collapsed by default so novices are not overwhelmed
//   - a search box filtering by flag name and help text
//   - the field's help text as a hover tooltip (+ the preset default)
//   - modified-from-preset highlighting and right-click "Reset to default"
//   - fields with `choices` as dropdowns, std::optional as auto/override
//
// A row added to the field table appears here automatically -- no GUI change
// needed.

#include "config/TrainConfig.h"

namespace gui {

struct ConfigUIState {
    char search[128] = "";
    bool modified_only = false;
};

// Draw the full generated editor. `defaults` is the preset-applied baseline
// used for modified-highlighting and reset. Returns true when any field
// changed this frame.
bool draw_config_editor(SsplatConfig& cfg, const SsplatConfig& defaults,
                        ConfigUIState& st);

// Shared helper: hoverable "(?)"-style tooltip with wrapped text.
void help_tooltip_on_hover(const char* text);

}  // namespace gui
