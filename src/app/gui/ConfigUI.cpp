// ConfigUI.cpp -- see ConfigUI.h. Every widget below is expanded from the
// SS_CONFIG_FIELDS X-macro; keep this file free of per-field special
// cases (the point is that a new row in the field table shows up here with
// zero GUI work).

#include "app/gui/ConfigUI.h"

#include "app/gui/Ui.h"

#include "i18n/catalog/Gui.h"

#include "imgui.h"
#include "imgui_stdlib.h"

#include <array>
#include <cctype>
#include <cmath>
#include <cstring>
#include <optional>
#include <string>
#include <vector>

namespace msg = spirula::i18n::msg::gui;
using spirula::i18n::Msg;

namespace gui {

namespace {

constexpr float kFieldWidth = 170.0f;
const ImVec4 kModifiedColor(1.0f, 0.72f, 0.25f, 1.0f);

bool icontains(const char* hay, const char* needle) {
    if (!*needle) return true;
    size_t nlen = std::strlen(needle);
    for (const char* h = hay; *h; h++) {
        size_t i = 0;
        while (i < nlen && h[i] &&
               std::tolower((unsigned char)h[i]) == std::tolower((unsigned char)needle[i]))
            i++;
        if (i == nlen) return true;
    }
    return false;
}

// ---- value -> display string (tooltips / reset labels) ---------------------

std::string value_str(bool v)               { return v ? "true" : "false"; }
std::string value_str(int v)                { return std::to_string(v); }
std::string value_str(float v) {
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[32]; std::snprintf(buf, sizeof buf, "%g", v); return buf;
}
std::string value_str(const std::string& v) { return v.empty() ? "none" : v; }
template <typename T> std::string value_str(const std::optional<T>& v) {
    return v ? value_str(*v) : "auto";
}
template <typename T, size_t N> std::string value_str(const std::array<T, N>& v) {
    std::string s;
    for (size_t i = 0; i < N; i++) s += (i ? " " : "") + value_str(v[i]);
    return s;
}

// ---- per-type widgets -------------------------------------------------------

std::vector<std::string> split_choices(const char* choices) {
    std::vector<std::string> out;
    std::string ch = choices;
    size_t pos = 0;
    while (pos <= ch.size()) {
        size_t bar = ch.find('|', pos);
        out.push_back(ch.substr(pos, bar == std::string::npos
                                         ? std::string::npos : bar - pos));
        if (bar == std::string::npos) break;
        pos = bar + 1;
    }
    return out;
}

bool draw_value(bool& v, const char*) {
    return ui::CheckboxRaw("##v", &v);
}

bool draw_value(int& v, const char*) {
    ImGui::SetNextItemWidth(kFieldWidth);
    return ui::InputIntRaw("##v", &v);
}

bool draw_value(float& v, const char*) {
    ImGui::SetNextItemWidth(kFieldWidth);
    return ui::InputFloatRaw("##v", &v, "%g");
}

bool draw_value(std::string& v, const char* choices) {
    std::string ch = choices;
    if (ch.empty() || ch == "none") {
        ImGui::SetNextItemWidth(kFieldWidth);
        return ui::InputTextRaw("##v", &v);
    }
    // Literal[...] -> dropdown; the "none" token maps to the empty string.
    std::vector<std::string> opts = split_choices(choices);
    std::string cur = v.empty() ? "none" : v;
    bool changed = false;
    ImGui::SetNextItemWidth(kFieldWidth);
    // The choices are the literal values `--<flag>` accepts and config.json
    // stores; they are not translated.
    if (ui::BeginComboRaw("##v", cur.c_str())) {
        for (const auto& o : opts) {
            bool sel = (o == cur);
            if (ui::SelectableRaw(o, sel)) {
                v = (o == "none") ? "" : o;
                changed = true;
            }
            if (sel) ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    return changed;
}

bool draw_value(std::optional<bool>& v, const char*) {
    // Literal[True, False, None] -> tri-state dropdown.
    const char* cur = !v.has_value() ? "auto" : (*v ? "true" : "false");
    bool changed = false;
    ImGui::SetNextItemWidth(kFieldWidth);
    if (ui::BeginComboRaw("##v", cur)) {
        if (ui::SelectableRaw("auto",  !v.has_value())) { v = std::nullopt; changed = true; }
        if (ui::SelectableRaw("true",  v.has_value() && *v))  { v = true;  changed = true; }
        if (ui::SelectableRaw("false", v.has_value() && !*v)) { v = false; changed = true; }
        ImGui::EndCombo();
    }
    return changed;
}

template <typename T>
bool draw_value(std::optional<T>& v, const char* choices) {
    bool has = v.has_value();
    bool changed = false;
    if (ui::CheckboxRaw("##has", &has)) {
        v = has ? std::optional<T>(T{}) : std::nullopt;
        changed = true;
    }
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort))
        ui::SetTooltip(msg::cfg_unchecked_is_auto);
    ImGui::SameLine();
    if (v.has_value()) {
        T tmp = *v;
        if (draw_value(tmp, choices)) { v = tmp; changed = true; }
    } else {
        ui::TextDisabled(msg::cfg_auto);
    }
    return changed;
}

template <size_t N>
bool draw_value(std::array<float, N>& v, const char*) {
    ImGui::SetNextItemWidth(kFieldWidth + 90);
    return ImGui::InputScalarN("##v", ImGuiDataType_Float, v.data(), (int)N,
                               nullptr, nullptr, "%g");
}

template <size_t N>
bool draw_value(std::array<int, N>& v, const char*) {
    ImGui::SetNextItemWidth(kFieldWidth + 90);
    return ImGui::InputScalarN("##v", ImGuiDataType_S32, v.data(), (int)N);
}

// ---- one field row -----------------------------------------------------------

template <typename T>
bool field_row(const char* cli_key, T& v, const T& def,
               const char* choices, const char* help) {
    bool modified = !(v == def);
    ImGui::PushID(cli_key);
    bool changed = draw_value(v, choices);
    bool hovered = ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort |
                                        ImGuiHoveredFlags_AllowWhenDisabled);
    ImGui::OpenPopupOnItemClick("ctx", ImGuiPopupFlags_MouseButtonRight);
    ImGui::SameLine();
    // The flag name is the flag name in every language: this row is a direct
    // view of `spirula train --<key>`.
    if (modified) ui::TextColoredRaw(kModifiedColor, std::string(cli_key) + " *");
    else          ui::TextRaw(cli_key);
    hovered |= ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort |
                                    ImGuiHoveredFlags_AllowWhenDisabled);
    ImGui::OpenPopupOnItemClick("ctx", ImGuiPopupFlags_MouseButtonRight);

    if (hovered && ImGui::BeginTooltip()) {
        ImGui::PushTextWrapPos(420.0f);
        ui::TextColoredRaw(kModifiedColor, "--" + std::string(cli_key));
        // The field help is shared with `spirula train --help` and is still
        // English everywhere; translating the 190 of them is the last tier.
        if (*help) ui::TextRaw(help);
        else       ui::Text(msg::cfg_no_description);
        ImGui::Separator();
        ui::Text(msg::cfg_preset_default, {value_str(def)});
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
    if (ImGui::BeginPopup("ctx")) {
        if (ui::MenuItem(msg::cfg_reset_to, {value_str(def)}))
            { v = def; changed = true; }
        ImGui::EndPopup();
    }
    ImGui::PopID();
    return changed;
}

int group_index(const char* g) {
    if (!std::strcmp(g, "trainer"))     return 0;
    if (!std::strcmp(g, "dataparser"))  return 1;
    if (!std::strcmp(g, "datamanager")) return 2;
    if (!std::strcmp(g, "model"))       return 3;
    return 4;  // optimizer
}

const Msg& group_label(const char* g) {
    switch (group_index(g)) {
        case 0: return msg::cfg_group_run;
        case 1: return msg::cfg_group_dataparser;
        case 2: return msg::cfg_group_datamanager;
        case 3: return msg::cfg_group_model;
        default: return msg::cfg_group_optimizer;
    }
}

// Fields owned by dedicated GUI controls; hidden from the generated list.
bool gui_managed(const char* cli_key) {
    return !std::strcmp(cli_key, "data");
}

}  // namespace


bool draw_config_editor(TrainConfig& cfg, const TrainConfig& defaults,
                        ConfigUIState& st) {
    ImGui::SetNextItemWidth(-115);
    ui::InputTextHintBufRaw("##cfgsearch", msg::cfg_search_hint,
                            st.search, sizeof st.search);
    ImGui::SameLine();
    ui::Checkbox(msg::cfg_edited_only, &st.modified_only);
    ui::help_on_hover(msg::cfg_edited_only_help);

    const bool searching = st.search[0] != 0 || st.modified_only;

    auto passes = [&](const char* key, const char* help, bool modified) {
        if (gui_managed(key)) return false;
        if (st.modified_only && !modified) return false;
        if (st.search[0] && !icontains(key, st.search) && !icontains(help, st.search))
            return false;
        return true;
    };

    // Pass 1: per-group visible-field counts (groups are contiguous in the
    // field table, so pass 2 can stream group headers).
    int vis[5] = {0, 0, 0, 0, 0};
#define SS_COUNT(type, member, default_, group, choices, help)                 \
    if (passes(#member, help, !(cfg.member == defaults.member)))               \
        vis[group_index(group)]++;
    SS_CONFIG_FIELDS(SS_COUNT)
#undef SS_COUNT

    // Pass 2: draw.
    bool any_changed = false;
    const char* cur_group = "";
    bool group_open = false;
#define SS_DRAW(type, member, default_, group, choices, help)                  \
    if (std::strcmp(cur_group, group) != 0) {                                  \
        cur_group = group;                                                     \
        if (vis[group_index(group)] == 0) {                                    \
            group_open = false;                                                \
        } else {                                                               \
            if (searching) ImGui::SetNextItemOpen(true);                       \
            group_open = ui::CollapsingHeader(group_label(group));             \
        }                                                                      \
    }                                                                          \
    if (group_open && passes(#member, help, !(cfg.member == defaults.member))) \
        any_changed |= field_row(#member, cfg.member, defaults.member, choices, help);
    SS_CONFIG_FIELDS(SS_DRAW)
#undef SS_DRAW

    return any_changed;
}

}  // namespace gui
