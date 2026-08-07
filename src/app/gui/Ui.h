#pragma once

// ui:: -- the only text-bearing ImGui entry points the GUI is allowed to call.
//
// `Msg` (src/i18n/Message.h) makes an INCOMPLETE translation a compile error.
// It cannot do anything about ImGui::Button("Start"), which compiles perfectly
// well and is simply never translated. This header closes that hole from the
// other side: every wrapper here takes a `const Msg&`, so a bare literal does
// not compile, and tools/check_i18n.sh fails the build if any of the wrapped
// ImGui functions is called anywhere in src/app/gui/ except from this file.
//
// Between the two, "forgot to translate" has no way in:
//   type system  -> a message with a missing language
//   the lint     -> a string that never became a message
//
// The `*Raw` overloads are the deliberate exceptions, and they are named so
// that using one is a visible decision in the diff. They are for text that
// MUST NOT be translated:
//   - paths, filenames, numbers, device names
//   - engine log lines and backend error strings (English, and grepped for)
//   - text already produced by i18n::format(), which is translated
//
// Widget IDs: everything with an ID gets "text###<message name>", so a widget
// keeps its identity when the language changes. Without it, switching to
// Japanese would collapse every open header and reset every scroll position,
// because ImGui derives a widget's ID from its label. Two call sites sharing
// one message therefore share an ID too -- wrap one in ImGui::PushID() exactly
// as you would for two identical literals.

#include "i18n/Message.h"

#include "imgui.h"
#include "imgui_stdlib.h"

#include <cfloat>
#include <string>
#include <vector>

namespace ui {

using ::spirula::i18n::Arg;
using ::spirula::i18n::Msg;
using ::spirula::i18n::format;

namespace detail {

// "text###id" in a small rotating buffer. Valid until this thread has built
// four more labels, which is long enough for any single ImGui call and short
// enough to never grow. ImGui hashes what follows ### and renders what
// precedes it.
inline const char* label(const char* text, const char* id) {
    static thread_local std::string ring[4];
    static thread_local unsigned next = 0;
    std::string& s = ring[next++ & 3u];
    s.assign(text);
    s += "###";
    s += id;
    return s.c_str();
}

inline const char* label(const Msg& m) { return label(m.get(), m.id); }

inline const char* label(const std::string& text, const Msg& m) {
    return label(text.c_str(), m.id);
}

// Msg* list -> the const char*[] ImGui's Combo wants.
inline const std::vector<const char*>& items(
        std::initializer_list<const Msg*> ms) {
    static thread_local std::vector<const char*> v;
    v.clear();
    for (const Msg* m : ms) v.push_back(m->get());
    return v;
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Text
// ---------------------------------------------------------------------------

inline void Text(const Msg& m) { ImGui::TextUnformatted(m.get()); }
inline void Text(const Msg& m, std::initializer_list<Arg> a) {
    ImGui::TextUnformatted(format(m, a).c_str());
}

inline void TextWrapped(const Msg& m) {
    ImGui::PushTextWrapPos();
    ImGui::TextUnformatted(m.get());
    ImGui::PopTextWrapPos();
}
inline void TextWrapped(const Msg& m, std::initializer_list<Arg> a) {
    ImGui::PushTextWrapPos();
    ImGui::TextUnformatted(format(m, a).c_str());
    ImGui::PopTextWrapPos();
}

inline void TextDisabled(const Msg& m) {
    ImGui::TextDisabled("%s", m.get());
}
inline void TextDisabled(const Msg& m, std::initializer_list<Arg> a) {
    ImGui::TextDisabled("%s", format(m, a).c_str());
}

inline void TextColored(const ImVec4& c, const Msg& m) {
    ImGui::TextColored(c, "%s", m.get());
}
inline void TextColored(const ImVec4& c, const Msg& m,
                        std::initializer_list<Arg> a) {
    ImGui::TextColored(c, "%s", format(m, a).c_str());
}

inline void TextColoredWrapped(const ImVec4& c, const Msg& m) {
    ImGui::PushTextWrapPos();
    ImGui::TextColored(c, "%s", m.get());
    ImGui::PopTextWrapPos();
}
inline void TextColoredWrapped(const ImVec4& c, const Msg& m,
                               std::initializer_list<Arg> a) {
    ImGui::PushTextWrapPos();
    ImGui::TextColored(c, "%s", format(m, a).c_str());
    ImGui::PopTextWrapPos();
}

// ---- raw: paths, numbers, log lines, and already-formatted text ----
inline void TextRaw(const char* s)        { ImGui::TextUnformatted(s); }
inline void TextRaw(const std::string& s) { ImGui::TextUnformatted(s.c_str()); }
inline void TextDisabledRaw(const char* s)        { ImGui::TextDisabled("%s", s); }
inline void TextDisabledRaw(const std::string& s) { ImGui::TextDisabled("%s", s.c_str()); }
inline void TextColoredRaw(const ImVec4& c, const char* s) {
    ImGui::TextColored(c, "%s", s);
}
inline void TextColoredRaw(const ImVec4& c, const std::string& s) {
    ImGui::TextColored(c, "%s", s.c_str());
}
inline void TextWrappedRaw(const std::string& s) {
    ImGui::PushTextWrapPos();
    ImGui::TextUnformatted(s.c_str());
    ImGui::PopTextWrapPos();
}
inline void TextColoredWrappedRaw(const ImVec4& c, const std::string& s) {
    ImGui::PushTextWrapPos();
    ImGui::TextColored(c, "%s", s.c_str());
    ImGui::PopTextWrapPos();
}

// ---------------------------------------------------------------------------
// Tooltips
// ---------------------------------------------------------------------------
// The delayed "(?)"-style hover help every option on the dataset and training
// screens carries. Kept here rather than in ConfigUI.h so that the wrapped-text
// tooltip and the widgets it annotates cannot drift apart.

inline void help_on_hover_raw(const char* text) {
    if (!text || !*text) return;
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_DelayShort) &&
        ImGui::BeginTooltip()) {
        ImGui::PushTextWrapPos(420.0f);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}
inline void help_on_hover(const Msg& m) { help_on_hover_raw(m.get()); }
inline void help_on_hover(const Msg& m, std::initializer_list<Arg> a) {
    help_on_hover_raw(format(m, a).c_str());
}

inline void SetTooltip(const Msg& m) { ImGui::SetTooltip("%s", m.get()); }
inline void SetTooltip(const Msg& m, std::initializer_list<Arg> a) {
    ImGui::SetTooltip("%s", format(m, a).c_str());
}
inline void SetTooltipRaw(const char* s) { ImGui::SetTooltip("%s", s); }
inline void SetTooltipRaw(const std::string& s) { ImGui::SetTooltip("%s", s.c_str()); }

// ---------------------------------------------------------------------------
// Buttons and toggles
// ---------------------------------------------------------------------------

inline bool Button(const Msg& m, const ImVec2& size = ImVec2(0, 0)) {
    return ImGui::Button(detail::label(m), size);
}
inline bool Button(const Msg& m, std::initializer_list<Arg> a,
                   const ImVec2& size = ImVec2(0, 0)) {
    return ImGui::Button(detail::label(format(m, a), m), size);
}
inline bool SmallButton(const Msg& m) {
    return ImGui::SmallButton(detail::label(m));
}
// A button whose face is punctuation or a glyph ("...", "<"), not a word.
inline bool ButtonRaw(const char* id, const ImVec2& size = ImVec2(0, 0)) {
    return ImGui::Button(id, size);
}
inline bool Checkbox(const Msg& m, bool* v) {
    return ImGui::Checkbox(detail::label(m), v);
}
inline bool CheckboxRaw(const char* id, bool* v) {
    return ImGui::Checkbox(id, v);
}
inline bool RadioButton(const Msg& m, bool active) {
    return ImGui::RadioButton(detail::label(m), active);
}
inline bool Selectable(const Msg& m, bool selected = false) {
    return ImGui::Selectable(detail::label(m), selected);
}
inline bool SelectableRaw(const char* s, bool selected = false) {
    return ImGui::Selectable(s, selected);
}
inline bool SelectableRaw(const std::string& s, bool selected = false) {
    return ImGui::Selectable(s.c_str(), selected);
}
inline bool SelectableRaw(const char* s, bool selected, ImGuiSelectableFlags f) {
    return ImGui::Selectable(s, selected, f);
}
inline bool RadioButtonRaw(const char* s, bool active) {
    return ImGui::RadioButton(s, active);
}

// ---------------------------------------------------------------------------
// Menus, headers, popups
// ---------------------------------------------------------------------------

inline bool BeginMenu(const Msg& m, bool enabled = true) {
    return ImGui::BeginMenu(detail::label(m), enabled);
}
// For the language menu, whose label is an icon plus a language's own name --
// neither of which is translated, and both of which would be wrong to
// translate. Pass detail::label(text, msg) so it still keeps a stable ID.
inline bool BeginMenuRaw(const char* s, bool enabled = true) {
    return ImGui::BeginMenu(s, enabled);
}
inline bool MenuItem(const Msg& m, const char* shortcut = nullptr,
                     bool selected = false, bool enabled = true) {
    return ImGui::MenuItem(detail::label(m), shortcut, selected, enabled);
}
inline bool MenuItem(const Msg& m, const char* shortcut, bool* selected,
                     bool enabled = true) {
    return ImGui::MenuItem(detail::label(m), shortcut, selected, enabled);
}
inline bool MenuItem(const Msg& m, std::initializer_list<Arg> a) {
    return ImGui::MenuItem(detail::label(format(m, a), m));
}
inline bool CollapsingHeader(const Msg& m, ImGuiTreeNodeFlags flags = 0) {
    return ImGui::CollapsingHeader(detail::label(m), flags);
}
inline bool TreeNode(const Msg& m) { return ImGui::TreeNode(detail::label(m)); }
inline void SeparatorText(const Msg& m) {
    ImGui::SeparatorText(detail::label(m));
}
inline void OpenPopup(const Msg& m) { ImGui::OpenPopup(detail::label(m)); }
inline bool BeginPopupModal(const Msg& m, bool* open = nullptr,
                            ImGuiWindowFlags flags = 0) {
    return ImGui::BeginPopupModal(detail::label(m), open, flags);
}
// A modal whose title the caller already built -- the file dialog is opened
// with a different title depending on what it is picking.
inline bool BeginPopupModalRaw(const char* title, bool* open = nullptr,
                               ImGuiWindowFlags flags = 0) {
    return ImGui::BeginPopupModal(title, open, flags);
}

// ---------------------------------------------------------------------------
// Value editors
// ---------------------------------------------------------------------------
// A combo's items are messages too -- pass them as a braced list of pointers:
//   ui::Combo(msg::quality, &q, {&msg::q_fast, &msg::q_balanced});

inline bool Combo(const Msg& m, int* cur, std::initializer_list<const Msg*> its) {
    const auto& v = detail::items(its);
    return ImGui::Combo(detail::label(m), cur, v.data(), (int)v.size());
}
inline bool ComboRaw(const char* id, int* cur, const char* const items[],
                     int count) {
    return ImGui::Combo(id, cur, items, count);
}
// Unlabelled combo, translated items: the viewport's toolbar is a single row
// of controls with no room for labels, but what they offer is still words.
inline bool ComboRaw(const char* id, int* cur,
                     std::initializer_list<const Msg*> its) {
    const auto& v = detail::items(its);
    return ImGui::Combo(id, cur, v.data(), (int)v.size());
}
inline bool ComboRaw(const char* id, int* cur, const std::vector<const Msg*>& its) {
    static thread_local std::vector<const char*> v;
    v.clear();
    for (const Msg* m : its) v.push_back(m->get());
    return ImGui::Combo(id, cur, v.data(), (int)v.size());
}
inline bool BeginCombo(const Msg& m, const char* preview) {
    return ImGui::BeginCombo(detail::label(m), preview);
}
inline bool BeginComboRaw(const char* id, const char* preview) {
    return ImGui::BeginCombo(id, preview);
}

inline bool SliderInt(const Msg& m, int* v, int lo, int hi) {
    return ImGui::SliderInt(detail::label(m), v, lo, hi);
}
// Unlabelled sliders whose value format string carries the meaning
// ("size x%.2f", "fov %.0f"). The format is a printf pattern ImGui fills in,
// so it cannot be a Msg -- but the words inside it can be, via format().
inline bool SliderIntRaw(const char* id, int* v, int lo, int hi,
                         const char* fmt) {
    return ImGui::SliderInt(id, v, lo, hi, fmt);
}
inline bool SliderFloatRaw(const char* id, float* v, float lo, float hi,
                           const char* fmt, ImGuiSliderFlags flags = 0) {
    return ImGui::SliderFloat(id, v, lo, hi, fmt, flags);
}
inline bool SliderFloat(const Msg& m, float* v, float lo, float hi,
                        const char* fmt = "%.3f") {
    return ImGui::SliderFloat(detail::label(m), v, lo, hi, fmt);
}
inline bool InputInt(const Msg& m, int* v, int step = 0, int step_fast = 0) {
    return ImGui::InputInt(detail::label(m), v, step, step_fast);
}
inline bool InputIntRaw(const char* id, int* v) {
    return ImGui::InputInt(id, v, 0, 0);
}
inline bool InputFloat(const Msg& m, float* v, float step = 0.0f,
                       float step_fast = 0.0f, const char* fmt = "%.3f") {
    return ImGui::InputFloat(detail::label(m), v, step, step_fast, fmt);
}
inline bool InputText(const Msg& m, std::string* v,
                      ImGuiInputTextFlags flags = 0) {
    return ImGui::InputText(detail::label(m), v, flags);
}
inline bool InputTextWithHint(const Msg& m, const Msg& hint, std::string* v,
                              ImGuiInputTextFlags flags = 0) {
    return ImGui::InputTextWithHint(detail::label(m), hint.get(), v, flags);
}

// A field whose CONTENT is English whatever the interface language is,
// because it is not read by a person: a mask prompt goes to a text encoder
// trained on English (src/app/gui/MaskPrompt.h). The label IS translated --
// it tells the user what the field is for -- but `example` is not, because it
// is an example of what to type into an English box, and translating it would
// be advice to type something that works worse.
inline bool InputTextEnglish(const Msg& m, const char* example, std::string* v,
                             ImGuiInputTextFlags flags = 0) {
    return ImGui::InputTextWithHint(detail::label(m), example, v, flags);
}
inline bool InputTextEnglishRaw(const char* id, const char* example,
                                std::string* v,
                                ImGuiInputTextFlags flags = 0) {
    return ImGui::InputTextWithHint(id, example, v, flags);
}

// Unlabelled fields ("##outdir"): no text to translate, but they still go
// through here so the lint does not have to special-case them.
inline bool InputTextRaw(const char* id, std::string* v,
                         ImGuiInputTextFlags flags = 0) {
    return ImGui::InputText(id, v, flags);
}
inline bool InputTextWithHintRaw(const char* id, const Msg& hint,
                                 std::string* v,
                                 ImGuiInputTextFlags flags = 0) {
    return ImGui::InputTextWithHint(id, hint.get(), v, flags);
}
inline bool InputFloatRaw(const char* id, float* v, const char* fmt) {
    return ImGui::InputFloat(id, v, 0.0f, 0.0f, fmt);
}
inline bool InputTextHintBufRaw(const char* id, const Msg& hint, char* buf,
                                size_t buf_size) {
    return ImGui::InputTextWithHint(id, hint.get(), buf, buf_size);
}

// ---------------------------------------------------------------------------
// Tables
// ---------------------------------------------------------------------------
// Only the column headers carry text; BeginTable / TableNextRow /
// TableNextColumn / TableHeadersRow do not, and are called directly.

inline void TableSetupColumn(const Msg& m, ImGuiTableColumnFlags flags = 0,
                             float width = 0.0f) {
    ImGui::TableSetupColumn(detail::label(m), flags, width);
}
// A column whose heading is a glyph or a number ("#"), not a word.
inline void TableSetupColumnRaw(const char* s, ImGuiTableColumnFlags flags = 0,
                                float width = 0.0f) {
    ImGui::TableSetupColumn(s, flags, width);
}

// ---------------------------------------------------------------------------
// Progress
// ---------------------------------------------------------------------------

inline void ProgressBar(float frac, const ImVec2& size, const Msg& m) {
    ImGui::ProgressBar(frac, size, m.get());
}
inline void ProgressBar(float frac, const ImVec2& size, const Msg& m,
                        std::initializer_list<Arg> a) {
    ImGui::ProgressBar(frac, size, format(m, a).c_str());
}
inline void ProgressBarRaw(float frac, const ImVec2& size, const char* overlay) {
    ImGui::ProgressBar(frac, size, overlay);
}
inline void PlotLinesRaw(const char* id, const float* values, int count,
                         const char* overlay, const ImVec2& size) {
    ImGui::PlotLines(id, values, count, 0, overlay, FLT_MAX, FLT_MAX, size);
}

}  // namespace ui
