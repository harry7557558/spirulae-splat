#pragma once

// ConfigUI -- the "Advanced options" editor. All widgets are generated from
// the SSPLAT_CONFIG_FIELDS X-macro (generated/cli_config.h), so every one of
// the training config fields is editable, with:
//   - grouping by Python sub-config (Run / Dataset / Data loading / Model /
//     Optimizer), collapsed by default so novices are not overwhelmed
//   - a search box filtering by flag name and help text
//   - the full Python docstring as a hover tooltip (+ the preset default)
//   - modified-from-preset highlighting and right-click "Reset to default"
//   - Literal[...] fields as dropdowns, Optional[...] as auto/override
//
// New config fields added on the Python side appear here automatically after
// re-running generate_cli_config.py -- no GUI change needed.

#include "../generated/cli_config.h"

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
