// The three consumers of SFM_CONFIG_FIELDS: the flag parser, the preset
// appliers, and the `--help` printer. Each is one macro expansion over the
// table in SfmConfig.h, so a new knob is one row and never a fourth edit.
#include "sfm/SfmConfig.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace sfm {
namespace {

// ---------------------------------------------------------------------------
// Value <-> string, one overload set per field type
// ---------------------------------------------------------------------------

std::string valueString(bool v) { return v ? "on" : "off"; }
std::string valueString(const std::string& v) { return v.empty() ? "none" : v; }

template <class T>
std::enable_if_t<std::is_integral<T>::value && !std::is_same<T, bool>::value, std::string>
valueString(T v) {
    char buf[32];
    if (std::is_signed<T>::value) std::snprintf(buf, sizeof buf, "%lld", (long long)v);
    else std::snprintf(buf, sizeof buf, "%llu", (unsigned long long)v);
    return buf;
}

template <class T>
std::enable_if_t<std::is_floating_point<T>::value, std::string> valueString(T v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%g", (double)v);
    return buf;
}

// ---------------------------------------------------------------------------
// Parsing one token into one field
// ---------------------------------------------------------------------------

std::string rangeText(double lo, double hi) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "expected a value in [%g, %g]", lo, hi);
    return buf;
}

template <class T>
std::enable_if_t<std::is_integral<T>::value && !std::is_same<T, bool>::value, bool>
parseValue(T& out, const std::string& s, double lo, double hi, const char*, std::string& err) {
    long long v = 0;
    try {
        size_t p = 0;
        v = std::stoll(s, &p);
        if (p != s.size()) throw std::invalid_argument("trailing characters");
    } catch (...) {
        err = "expected an integer";
        return false;
    }
    if (std::is_unsigned<T>::value && v < 0) {
        err = "must not be negative";
        return false;
    }
    if (lo < hi && ((double)v < lo || (double)v > hi)) {
        err = rangeText(lo, hi);
        return false;
    }
    out = (T)v;
    return true;
}

template <class T>
std::enable_if_t<std::is_floating_point<T>::value, bool> parseValue(T& out, const std::string& s,
                                                                    double lo, double hi,
                                                                    const char*, std::string& err) {
    double v = 0;
    try {
        size_t p = 0;
        v = std::stod(s, &p);
        if (p != s.size()) throw std::invalid_argument("trailing characters");
    } catch (...) {
        err = "expected a number";
        return false;
    }
    if (lo < hi && (v < lo || v > hi)) {
        err = rangeText(lo, hi);
        return false;
    }
    out = (T)v;
    return true;
}

bool parseValue(std::string& out, const std::string& s, double, double, const char* choices,
                std::string& err) {
    std::string ch = choices ? choices : "";
    if (!ch.empty()) {
        for (size_t pos = 0; pos <= ch.size();) {
            size_t bar = ch.find('|', pos);
            std::string tok = ch.substr(pos, bar == std::string::npos ? bar : bar - pos);
            if (tok == s) { out = s; return true; }
            if (bar == std::string::npos) break;
            pos = bar + 1;
        }
        err = "expected one of: " + ch;
        for (char& c : err) if (c == '|') c = ' ';
        return false;
    }
    out = s;
    return true;
}

// A bool field is a switch, never a value: `--name` / `--no-name`. Keeping it
// valueless is what lets `ssplat-sfm map --no-manage matches.bin feats/` parse
// -- a positional right after a switch must stay a positional.
FieldResult trySetField(bool& out, const std::string& key, const char* name, int, char**, int&,
                        double, double, const char*, std::set<std::string>& seen, std::string&) {
    if (key == name) out = true;
    else if (key.size() > 3 && key.compare(0, 3, "no-") == 0 && key.compare(3, std::string::npos, name) == 0)
        out = false;
    else return FieldResult::Unknown;
    seen.insert(name);
    return FieldResult::Ok;
}

template <class T>
FieldResult trySetField(T& out, const std::string& key, const char* name, int argc, char** argv,
                        int& i, double lo, double hi, const char* choices,
                        std::set<std::string>& seen, std::string& error) {
    if (key != name) return FieldResult::Unknown;
    if (i + 1 >= argc) {
        error = std::string("--") + name + ": missing value";
        return FieldResult::Error;
    }
    std::string why;
    if (!parseValue(out, argv[++i], lo, hi, choices, why)) {
        error = std::string("--") + name + " " + argv[i] + ": " + why;
        return FieldResult::Error;
    }
    seen.insert(name);
    return FieldResult::Ok;
}

// "--Max_Error" and "--max-error" are the same flag; '_' is accepted because
// the GUI and the config table spell ids either way.
std::string normalizeKey(const std::string& arg) {
    std::string s = arg;
    while (!s.empty() && s[0] == '-') s.erase(s.begin());
    for (char& c : s) if (c == '_') c = '-';
    return s;
}

// ---------------------------------------------------------------------------
// Help formatting
// ---------------------------------------------------------------------------

// Long choice lists do not fit the flag column; they move into the help text.
constexpr size_t kInlineChoices = 26;

std::string metavarFor(bool, const char*, const char*) { return ""; }
std::string metavarFor(const std::string&, const char* name, const char* choices) {
    if (choices && *choices && std::string(choices).size() <= kInlineChoices)
        return "{" + std::string(choices) + "}";
    // No metavar column in the table: the flag name says what it takes, and
    // "DIR" reads better than "VALUE" on the handful that take a path.
    std::string n = name;
    if (n.find("dir") != std::string::npos || n == "masks" || n == "images" ||
        n == "features" || n == "resume")
        return "DIR";
    if (n.find("path") != std::string::npos) return "FILE";
    return "VALUE";
}
template <class T>
std::enable_if_t<std::is_integral<T>::value && !std::is_same<T, bool>::value, std::string>
metavarFor(T, const char*, const char*) {
    return "N";
}
template <class T>
std::enable_if_t<std::is_floating_point<T>::value, std::string> metavarFor(T, const char*,
                                                                          const char*) {
    return "X";
}

// A bool is offered in the direction that changes it, which is how it is
// written on a real command line (`--no-verify`, not `--verify 0`).
std::string flagFor(bool v, const char* name) {
    return v ? "--no-" + std::string(name) : "--" + std::string(name);
}
template <class T>
std::string flagFor(const T&, const char* name) {
    return "--" + std::string(name);
}

constexpr int kFlagCol = 42;  // flag column width; the widest flag+metavar fits
constexpr int kWidth = 96;

void printWrapped(FILE* out, const std::string& text, int indent) {
    std::string line;
    size_t pos = 0;
    while (pos < text.size()) {
        size_t sp = text.find(' ', pos);
        std::string word = text.substr(pos, sp == std::string::npos ? sp : sp - pos);
        if (!line.empty() && (int)(line.size() + 1 + word.size()) + indent > kWidth) {
            std::fprintf(out, "%*s%s\n", indent, "", line.c_str());
            line.clear();
        }
        line += (line.empty() ? "" : " ") + word;
        if (sp == std::string::npos) break;
        pos = sp + 1;
    }
    if (!line.empty()) std::fprintf(out, "%*s%s\n", indent, "", line.c_str());
}

void printOption(FILE* out, const std::string& flag, const std::string& metavar,
                 const std::string& value, const std::string& help, const char* choices) {
    std::string left = "  " + flag + (metavar.empty() ? "" : " " + metavar);
    if (value.empty()) {
        std::fprintf(out, "%s\n", left.c_str());
    } else {
        if ((int)left.size() >= kFlagCol) std::fprintf(out, "%s\n%*s", left.c_str(), kFlagCol, "");
        else std::fprintf(out, "%-*s", kFlagCol, left.c_str());
        std::fprintf(out, "[%s]\n", value.c_str());
    }
    std::string h = help;
    std::string ch = choices ? choices : "";
    if (!ch.empty() && ch.size() > kInlineChoices) {
        for (char& c : ch) if (c == '|') c = ' ';
        h += " (one of: " + ch + ")";
    }
    printWrapped(out, h, 6);
}

}  // namespace

// ---------------------------------------------------------------------------
// setConfigField
// ---------------------------------------------------------------------------

FieldResult setConfigField(SfmConfig& cfg, uint32_t cmd, const std::string& arg, int argc,
                           char** argv, int& i, std::set<std::string>& seen, std::string& error) {
    const std::string key = normalizeKey(arg);
    if (key.empty()) return FieldResult::Unknown;
#define SFM_TRY_SET(member, name, cmds, tier, group, lo, hi, choices, help)                        \
    if ((uint32_t)(cmds) & cmd) {                                                                  \
        FieldResult r = trySetField(cfg.member, key, name, argc, argv, i, (double)(lo),            \
                                    (double)(hi), choices, seen, error);                           \
        if (r != FieldResult::Unknown) return r;                                                   \
    }
    SFM_CONFIG_FIELDS(SFM_TRY_SET)
#undef SFM_TRY_SET
    return FieldResult::Unknown;
}

// ---------------------------------------------------------------------------
// Presets
// ---------------------------------------------------------------------------

namespace {
// Move a field unless the command line already claimed it, and say so. The
// report is what tells a user why `--quality low` produced a 1000 px pipeline,
// and it is the same list the GUI highlights as modified.
template <class T, class V>
void presetSet(std::set<std::string> const& seen, std::vector<PresetChange>& moved,
               const char* flag, T& member, V value) {
    if (seen.count(flag)) return;
    T next = (T)value;
    if (next == member) return;
    moved.push_back({flag, valueString(member), valueString(next)});
    member = next;
}
}  // namespace

std::string applyPresets(SfmConfig& cfg, const std::set<std::string>& seen,
                         std::vector<PresetChange>& moved) {
    // Quality mirrors COLMAP's OptionManager::ModifyFor*Quality() as fractions
    // of the extractor's 3200 px default. Pair-selection breadth follows it too
    // (D42): `prefilter-neighbors` is the one selection parameter that trades
    // match time against how much of the view graph verification even sees.
    if (cfg.quality == "low") {
        presetSet(seen, moved, "max-image-size", cfg.max_image_size, 1000);
        presetSet(seen, moved, "max-features", cfg.sift.max_num_features, 2048);
        presetSet(seen, moved, "prefilter-neighbors", cfg.prefilter.num_neighbors, 16);
    } else if (cfg.quality == "medium") {
        presetSet(seen, moved, "max-image-size", cfg.max_image_size, 1600);
        presetSet(seen, moved, "max-features", cfg.sift.max_num_features, 4096);
        presetSet(seen, moved, "prefilter-neighbors", cfg.prefilter.num_neighbors, 24);
    } else if (cfg.quality == "high") {
        presetSet(seen, moved, "max-image-size", cfg.max_image_size, 2400);
        presetSet(seen, moved, "max-features", cfg.sift.max_num_features, 8192);
    } else if (cfg.quality == "extreme") {
        presetSet(seen, moved, "max-image-size", cfg.max_image_size, 3200);
        presetSet(seen, moved, "max-features", cfg.sift.max_num_features, 16384);
        presetSet(seen, moved, "prefilter-neighbors", cfg.prefilter.num_neighbors, 48);
    } else {
        return "unknown --quality '" + cfg.quality + "' (low, medium, high or extreme)";
    }

    if (cfg.data_type == "individual") {
        // Nothing: the defaults are written for a set of individual photos.
    } else if (cfg.data_type == "video") {
        presetSet(seen, moved, "pairs", cfg.pairs, std::string("sequential"));
        // COLMAP ModifyForVideoData: consecutive frames subtend far less angle,
        // so the seed pair cannot be held to the photo-collection threshold.
        presetSet(seen, moved, "init-min-tri-angle", cfg.mapper.init_min_tri_angle_deg,
                  cfg.mapper.init_min_tri_angle_deg / 2);
    } else if (cfg.data_type == "internet") {
        // Two Flickr uploads that happen to be 1024x768 are two cameras that
        // were each downscaled, so grouping on resolution fuses unrelated
        // intrinsics. Give every image its own (D20).
        presetSet(seen, moved, "camera-mode", cfg.camera_mode, std::string("image"));
        cfg.camera_mode_pinned = true;
    } else {
        return "unknown --data-type '" + cfg.data_type + "' (individual, video or internet)";
    }
    return "";
}

// ---------------------------------------------------------------------------
// finalize
// ---------------------------------------------------------------------------

std::string SfmConfig::finalize(uint32_t cmd) {
    if (!parseCamModelName(camera_model, camera.model))
        return "unknown --camera-model '" + camera_model + "'";
    if (!parseCameraMode(camera_mode, camera.mode))
        return "unknown --camera-mode '" + camera_mode + "'";
    camera.mode_explicit = camera_mode_pinned;
    camera.focal = focal;
    // The model a group with no explicit entry falls back to, which is the same
    // question --camera-model answers for the ones that have entries.
    mapper.camera_model = camera.model;

    // One tolerance, two fields (D47).
    twoview.ransac.max_error = max_error;
    mapper.max_reproj_error = max_error;

    sift.device = match.device = prefilter.device = mapper.device = device;
    mapper.threads = threads;
    const bool v = !quiet;
    sift.verbose = mapper.verbose = manager.verbose = merge.verbose = v;

    // Scoring problems are ~1/32 the size of full matching, so the selection
    // pass batches at least as many pairs per submit as the matcher does.
    prefilter.batch_pairs = std::max(prefilter.batch_pairs, match.batch_pairs);

    // Post-merge filtering follows the mapper, because spliced tracks have
    // never seen a bundle adjustment across the seam and should not survive
    // criteria the mapper's own observations would not. `merge` is the
    // exception: there is no mapper in that command, and --filter-error /
    // --min-tri-angle are the flags that set these directly.
    if (cmd != CMD_MERGE) {
        merge.filter_reproj_error = mapper.max_reproj_error;
        merge.min_tri_angle_deg = mapper.min_tri_angle_deg;
        merge.min_image_points = mapper.min_image_points;
    }
    manager.merge = merge;
    manager.merge.verbose = v;
    return "";
}

PairMode SfmConfig::pairMode() const {
    if (pairs == "sequential") return PairMode::Sequential;
    if (pairs == "prefilter") return PairMode::Prefilter;
    return PairMode::Exhaustive;  // "exhaustive", and "auto" until counted
}

// ---------------------------------------------------------------------------
// Help
// ---------------------------------------------------------------------------

void printConfigOptions(FILE* out, uint32_t cmd, const SfmConfig& defaults) {
    const char* cur_group = "";
#define SFM_PRINT_FIELD(member, name, cmds, tier, group, lo, hi, choices, help)                    \
    if (((uint32_t)(cmds) & cmd) && (tier) != Tier::Alias) {                                       \
        if (std::string(cur_group) != group) {                                                     \
            cur_group = group;                                                                     \
            std::fprintf(out, "\n %s:\n", group);                                                  \
        }                                                                                          \
        printOption(out, flagFor(defaults.member, name), metavarFor(defaults.member, name, choices),\
                    valueString(defaults.member), help, choices);                                  \
    }
    SFM_CONFIG_FIELDS(SFM_PRINT_FIELD)
#undef SFM_PRINT_FIELD
}

void printOptionLine(FILE* out, const std::string& flag, const std::string& value,
                     const std::string& help) {
    printOption(out, flag, "", value, help, "");
}

}  // namespace sfm
