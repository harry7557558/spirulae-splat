#pragma once

// A saved training preset: the whole training config, under a name the user
// chose, in one JSON file.
//
// The built-in presets (kTrainPresets, config/TrainConfig.h) are code -- a
// handful of curated starting points that ship with the program. This is the
// other half: once someone has tuned the options for their own captures, the
// tuning is worth keeping, sharing, and pointing a batch run at.
//
// The file:
//
//   {
//     "spirula_preset": 1,        // format marker
//     "name": "Indoor handheld",
//     "description": "...",       // may be absent
//     "base_preset": "3dgs",      // the built-in it started from
//     "touched": ["quality", ...],// flags set by hand; see below
//     "config": { <one key per flag, as a run's config.json> }
//   }
//
// Reading is forgiving in the way a settings file has to be: a key that is
// missing keeps the field's default, and a key that is not a training flag is
// ignored. So a preset written by an older build still loads into a newer one
// (new flags come up at their defaults), and a preset from a newer build loads
// into an older one (its unknown flags are dropped).
//
// A run's own config.json loads as a preset too. That file is the same flat
// table under a "preset" key instead of "base_preset", and reusing the
// settings of a run that came out well is the single most likely reason
// anybody wants this feature at all.
//
// `touched` is what stops the macro options from undoing the tuning. --quality
// and friends write the flags they stand for unless the user set those flags
// by hand (train_resolve_macros), and "by hand" is a GUI-session fact that
// would otherwise be lost the moment the preset is saved. A config.json has no
// such list, so one is derived: every flag that differs from the base preset
// counts as deliberate.

#include "config/TrainConfig.h"

#include <set>
#include <string>
#include <vector>

namespace gui {

// What a preset must NOT carry: where the data is and where the run goes. A
// preset answers "how", the trainer screen and a batch row answer "where" --
// so loading one never moves the dataset out from under you, and any preset
// can be paired with any dataset.
#define SS_PRESET_CONTEXT_FIELDS(X) \
    X(data) X(resume) X(output_dir_prefix) X(output_dir_name) \
    /* end */

struct TrainPreset {
    std::string name;           // what the user called it
    std::string description;    // free text, may be empty
    std::string base = "3dgs";  // the built-in preset it started from
    std::string path;           // the file it was read from / written to
    TrainConfig cfg;
    std::set<std::string> touched;
};

// Where presets live by default, created on first call:
// <config_dir>/presets. The picker lists what is in here; Save offers it as
// the default location but any path is allowed.
std::string preset_dir();

// A file name for `name` -- lowercased, spaces and separators folded to '-',
// always ending in ".json". Never empty (an unnameable name becomes
// "preset.json"), and never a path: the caller decides the folder.
std::string preset_file_name(const std::string& name);

// Write. Throws std::runtime_error if the file cannot be created. The context
// fields above are dropped here rather than at the call site, so no caller can
// forget and bake a dataset path into a shared preset.
void save_preset(const TrainPreset& p, const std::string& path);

// Read a preset file, or a run's config.json. Throws std::runtime_error when
// the file is unreadable, is not JSON, or names no training flag at all.
TrainPreset load_preset(const std::string& path);

// Cheap probe for drag-and-drop: would load_preset() accept this file? Never
// throws -- a JSON file that is something else entirely just answers false.
bool is_preset_file(const std::string& path);

// Delete a preset file. Throws std::runtime_error when the file is not a
// preset (nothing else may be deleted through here, whatever path a caller
// hands over) or when the filesystem refuses.
void delete_preset(const std::string& path);

// Every readable *.json in preset_dir(), sorted by name. Files that fail to
// parse are skipped rather than reported: this runs to fill a dropdown.
std::vector<TrainPreset> list_presets();

}  // namespace gui
