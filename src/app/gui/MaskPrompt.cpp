// MaskPrompt.cpp -- see MaskPrompt.h.

#include "app/gui/MaskPrompt.h"

#include "app/gui/Ui.h"
#include "i18n/catalog/Dataset.h"

#include <algorithm>
#include <cctype>

namespace dmsg = spirula::i18n::msg::dataset;

namespace gui {

namespace {

// The English terms. Chosen for what SAM 3 actually responds to -- a concrete
// countable noun -- rather than for what reads best in a menu. "my shadow"
// would be the natural English but the model has no idea whose shadow it is,
// so the term is "shadow" and the label carries the nuance.
const std::vector<MaskSubject>& subjects_impl() {
    static const std::vector<MaskSubject> v = {
        {&dmsg::subj_person,          "person"},
        {&dmsg::subj_hand,            "hand"},
        {&dmsg::subj_dog,             "dog"},
        {&dmsg::subj_animal,          "animal"},
        {&dmsg::subj_car,             "car"},
        {&dmsg::subj_bicycle,         "bicycle"},
        {&dmsg::subj_vehicle,         "vehicle"},
        {&dmsg::subj_sky,             "sky"},
        {&dmsg::subj_shadow,          "shadow"},
        {&dmsg::subj_water,           "water"},
        {&dmsg::subj_reflection,      "reflection"},
        {&dmsg::subj_tripod,          "tripod"},
        {&dmsg::subj_fisheye_border,  "black border"},
        {&dmsg::subj_watermark,       "watermark"},
    };
    return v;
}

const std::vector<MaskSubject>& exceptions_impl() {
    static const std::vector<MaskSubject> v = {
        {&dmsg::subj_person_painting, "person in a painting"},
        {&dmsg::subj_statue,          "statue"},
        {&dmsg::subj_mannequin,       "mannequin"},
    };
    return v;
}

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) a++;
    while (b > a && std::isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

std::string lower(std::string s) {
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

std::vector<std::string> split_terms(const std::string& prompt) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= prompt.size()) {
        const size_t sep = prompt.find(';', start);
        const std::string piece =
            trim(prompt.substr(start, sep == std::string::npos
                                          ? std::string::npos : sep - start));
        if (!piece.empty()) out.push_back(piece);
        if (sep == std::string::npos) break;
        start = sep + 1;
    }
    return out;
}

std::string join_terms(const std::vector<std::string>& terms) {
    std::string out;
    for (const std::string& t : terms) {
        if (!out.empty()) out += "; ";
        out += t;
    }
    return out;
}

// One chip. Highlighted when its term is already in the target field, which
// makes the palette double as a readout of what the field says.
bool chip(const MaskSubject& s, std::string& target) {
    const bool on = prompt_has_term(target, s.term);
    if (on) {
        const ImVec4& c = ImGui::GetStyle().Colors[ImGuiCol_ButtonActive];
        ImGui::PushStyleColor(ImGuiCol_Button, c);
    }
    const bool hit = ui::SmallButton(*s.label);
    if (on) ImGui::PopStyleColor();
    if (hit) prompt_toggle_term(target, s.term);
    return hit;
}

// Chips laid out left to right, wrapping at the content edge. ImGui has no
// flow layout, so the wrap is done by hand against the next chip's measured
// width -- SameLine() unconditionally would run off the panel in the
// languages whose words are longest (German, and Portuguese here).
bool chip_row(const std::vector<MaskSubject>& group, std::string& target) {
    const float right = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
    const float spacing = ImGui::GetStyle().ItemSpacing.x;
    bool edited = false;
    for (size_t i = 0; i < group.size(); i++) {
        if (i) {
            const float last = ImGui::GetItemRectMax().x;
            const float w = ImGui::CalcTextSize(group[i].label->get()).x +
                            ImGui::GetStyle().FramePadding.x * 2.0f;
            if (last + spacing + w < right) ImGui::SameLine();
        }
        edited |= chip(group[i], target);
    }
    return edited;
}

}  // namespace

const std::vector<MaskSubject>& mask_subjects() { return subjects_impl(); }
const std::vector<MaskSubject>& mask_exception_subjects() {
    return exceptions_impl();
}

bool prompt_has_term(const std::string& prompt, const char* term) {
    const std::string want = lower(trim(term));
    for (const std::string& t : split_terms(prompt))
        if (lower(t) == want) return true;
    return false;
}

void prompt_toggle_term(std::string& prompt, const char* term) {
    const std::string want = lower(trim(term));
    std::vector<std::string> terms = split_terms(prompt);
    const size_t before = terms.size();
    terms.erase(std::remove_if(terms.begin(), terms.end(),
                               [&](const std::string& t) {
                                   return lower(t) == want;
                               }),
                terms.end());
    if (terms.size() == before) terms.push_back(trim(term));
    prompt = join_terms(terms);
}

bool draw_subject_palette(std::string& prompt, std::string& negative,
                          bool keep_subject) {
    if (!ui::TreeNode(dmsg::mask_subjects)) return false;
    ui::help_on_hover(dmsg::mask_subjects_help);

    // Each group is headed by the label of the field it writes into. That is
    // the whole explanation of where a chip's word goes, and it stays correct
    // when the keep/remove radio flips both labels.
    bool edited = false;
    ui::TextDisabled(keep_subject ? dmsg::mask_what_to_keep
                                  : dmsg::mask_what_to_remove);
    edited |= chip_row(mask_subjects(), prompt);

    ImGui::Spacing();
    ui::TextDisabled(keep_subject ? dmsg::mask_but_remove : dmsg::mask_but_keep);
    edited |= chip_row(mask_exception_subjects(), negative);

    ImGui::TreePop();
    return edited;
}

}  // namespace gui
