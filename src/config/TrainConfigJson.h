#pragma once

// Flat JSON <-> TrainConfig, in one place.
//
// Three things put a training config on disk or read one back, and all three
// have to agree down to how an unset std::optional is spelled:
//
//   * a run's config.json          (app/TrainerCore.cpp save_config_json)
//   * --resume                     (checkpoint/Resume.cpp config_from_json)
//   * the GUI's saved presets      (app/gui/TrainPreset.cpp)
//
// The encoding: one key per flag, spelled exactly as the flag is, flat --
// never nested under the field table's `section`, which is presentational and
// free to be reshuffled. `null` means "unset": an empty std::string, a
// disengaged std::optional. Floats are Python-json compatible, so an infinite
// default serializes as `Infinity` rather than something no parser accepts.
//
// Reading is deliberately forgiving in one direction only: a key that is not
// in the table is ignored, and a field with no key keeps its default -- so a
// file written by an older or newer build still loads. A key that IS in the
// table but holds the wrong shape (a two-element array for a vec3) still
// throws, because that is a corrupt file rather than a version skew.

#include "config/TrainConfig.h"
#include "data/Json.h"

#include <cmath>
#include <cstdio>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace train_json_detail {

// ---- field -> JSON text ---------------------------------------------------

inline std::string emit(bool v) { return v ? "true" : "false"; }
inline std::string emit(int v)  { return std::to_string(v); }
inline std::string emit(float v) {
    if (std::isinf(v)) return v > 0 ? "Infinity" : "-Infinity";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%g", v);
    return buf;
}
inline std::string emit(const std::string& v) {
    if (v.empty()) return "null";
    std::string s = "\"";
    for (char ch : v) { if (ch == '"' || ch == '\\') s += '\\'; s += ch; }
    return s + "\"";
}
template <typename T> std::string emit(const std::optional<T>& v) {
    return v ? emit(*v) : "null";
}
template <typename T, size_t N> std::string emit(const std::array<T, N>& v) {
    std::string s = "[";
    for (size_t i = 0; i < N; i++) s += (i ? ", " : "") + emit(v[i]);
    return s + "]";
}

// ---- JSON value -> field --------------------------------------------------
// A missing or null value leaves the field at its default, except where null
// is the encoding of "unset": an empty std::string, a disengaged optional.

inline void assign(int& out, const JsonValue& v)   { if (!v.is_null()) out = (int)v.as_int(); }
inline void assign(float& out, const JsonValue& v) { if (!v.is_null()) out = (float)v.as_double(); }
inline void assign(bool& out, const JsonValue& v)  { if (!v.is_null()) out = v.as_bool(); }

inline void assign(std::string& out, const JsonValue& v) {
    out = v.is_null() ? std::string() : v.as_string();
}

inline void assign(std::optional<int>& out, const JsonValue& v) {
    if (v.is_null()) out = std::nullopt; else out = (int)v.as_int();
}
inline void assign(std::optional<float>& out, const JsonValue& v) {
    if (v.is_null()) out = std::nullopt; else out = (float)v.as_double();
}
inline void assign(std::optional<bool>& out, const JsonValue& v) {
    if (v.is_null()) out = std::nullopt; else out = v.as_bool();
}

template <typename T, size_t N>
void assign(std::array<T, N>& out, const JsonValue& v) {
    if (v.is_null()) return;
    if (!v.is_array() || v.arr.size() != N)
        throw std::runtime_error("config JSON: expected an array of " +
                                 std::to_string(N) + " numbers");
    for (size_t i = 0; i < N; i++)
        out[i] = (T)(std::is_integral<T>::value ? (double)v.arr[i].as_int()
                                                : v.arr[i].as_double());
}

}  // namespace train_json_detail


// Every flag as (key, JSON value text), in field-table order. The caller
// decides the punctuation and the indent, because the three writers wrap it
// differently -- a run's config.json puts these at the top level, a preset
// file nests them under "config".
inline std::vector<std::pair<const char*, std::string>>
train_config_json_pairs(const TrainConfig& c) {
    std::vector<std::pair<const char*, std::string>> out;
#define SS_JSON_PAIR(type, member, default_, section, tier, choices)          \
    out.emplace_back(#member, train_json_detail::emit(c.member));
    SS_CONFIG_FIELDS(SS_JSON_PAIR)
#undef SS_JSON_PAIR
    return out;
}

// The inverse. `root` is the object holding the flat keys; anything it does
// not name keeps whatever `out` already had, which is what lets a caller pass
// a preset-applied config as the baseline.
inline void train_config_from_json(const JsonValue& root, TrainConfig& out) {
#define SS_JSON_LOAD(type, member, default_, section, tier, choices)          \
    if (const JsonValue* v = root.find(#member))                              \
        train_json_detail::assign(out.member, *v);
    SS_CONFIG_FIELDS(SS_JSON_LOAD)
#undef SS_JSON_LOAD
}

// Does this object name any training flag at all? The probe that tells a
// config.json (or a preset) from some other JSON file that was dropped on the
// window by mistake.
inline bool train_config_json_has_fields(const JsonValue& root) {
    if (!root.is_object()) return false;
    bool any = false;
#define SS_JSON_PROBE(type, member, default_, section, tier, choices)         \
    any = any || root.find(#member) != nullptr;
    SS_CONFIG_FIELDS(SS_JSON_PROBE)
#undef SS_JSON_PROBE
    return any;
}
