// FileDialog.cpp -- see FileDialog.h.

#include "FileDialog.h"

#include "imgui.h"
#include "imgui_stdlib.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

namespace gui {

namespace {

std::string home_dir() {
#ifdef _WIN32
    if (const char* p = std::getenv("USERPROFILE")) return p;
    return "C:\\";
#else
    if (const char* p = std::getenv("HOME")) return p;
    return "/";
#endif
}

std::string lower(std::string s) {
    for (auto& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

}  // namespace

void FileDialog::open(const std::string& title, Mode mode,
                      const std::vector<std::string>& extensions,
                      const std::string& start_dir) {
    _title = title;
    _mode = mode;
    _extensions = extensions;
    std::error_code ec;
    if (!start_dir.empty() && fs::is_directory(start_dir, ec))
        _cwd = fs::absolute(start_dir, ec).string();
    else if (_cwd.empty() || !fs::is_directory(_cwd, ec))
        _cwd = home_dir();
    _selected.clear();
    _result.clear();
    _want_open = true;
    refresh();
}

void FileDialog::refresh() {
    _entries.clear();
    _path_edit = _cwd;
    std::error_code ec;
    for (fs::directory_iterator it(_cwd, fs::directory_options::skip_permission_denied, ec), end;
         !ec && it != end; it.increment(ec)) {
        Entry e;
        e.name = it->path().filename().string();
        if (e.name.empty() || e.name[0] == '.') continue;   // hide dotfiles
        e.is_dir = it->is_directory(ec);
        if (!e.is_dir) {
            if (_mode == Mode::Folder) continue;
            if (!_extensions.empty()) {
                std::string ext = lower(it->path().extension().string());
                if (std::find(_extensions.begin(), _extensions.end(), ext) ==
                    _extensions.end()) continue;
            }
        }
        _entries.push_back(std::move(e));
    }
    std::sort(_entries.begin(), _entries.end(), [](const Entry& a, const Entry& b) {
        if (a.is_dir != b.is_dir) return a.is_dir;
        return lower(a.name) < lower(b.name);
    });
}

bool FileDialog::draw() {
    if (_want_open) {
        ImGui::OpenPopup(_title.c_str());
        _want_open = false;
        _is_open = true;
    }
    if (!_is_open) return false;

    bool confirmed = false;
    ImGui::SetNextWindowSize(ImVec2(680, 480), ImGuiCond_Appearing);
    if (ImGui::BeginPopupModal(_title.c_str(), &_is_open)) {
        std::error_code ec;

        // ---- top bar: up / home / drives / editable path ----
        if (ImGui::Button("Up")) {
            fs::path p(_cwd);
            if (p.has_parent_path() && p.parent_path() != p) {
                _cwd = p.parent_path().string();
                _selected.clear();
                refresh();
            }
        }
        ImGui::SameLine();
        if (ImGui::Button("Home")) {
            _cwd = home_dir();
            _selected.clear();
            refresh();
        }
#ifdef _WIN32
        ImGui::SameLine();
        ImGui::SetNextItemWidth(70);
        if (ImGui::BeginCombo("##drives", "Drive")) {
            for (char d = 'A'; d <= 'Z'; d++) {
                std::string root = std::string(1, d) + ":\\";
                if (fs::exists(root, ec) && ImGui::Selectable(root.c_str())) {
                    _cwd = root;
                    _selected.clear();
                    refresh();
                }
            }
            ImGui::EndCombo();
        }
#endif
        ImGui::SameLine();
        ImGui::SetNextItemWidth(-1);
        if (ImGui::InputText("##path", &_path_edit,
                             ImGuiInputTextFlags_EnterReturnsTrue)) {
            if (fs::is_directory(_path_edit, ec)) {
                _cwd = fs::absolute(_path_edit, ec).string();
                _selected.clear();
                refresh();
            } else if (_mode == Mode::File && fs::is_regular_file(_path_edit, ec)) {
                _result = fs::absolute(_path_edit, ec).string();
                confirmed = true;
            }
        }

        // ---- listing ----
        float footer = ImGui::GetFrameHeightWithSpacing() + 8;
        if (ImGui::BeginChild("##list", ImVec2(0, -footer), ImGuiChildFlags_Borders)) {
            for (const auto& e : _entries) {
                std::string label = (e.is_dir ? "[+] " : "     ") + e.name;
                bool sel = (_selected == e.name);
                if (ImGui::Selectable(label.c_str(), sel,
                                      ImGuiSelectableFlags_AllowDoubleClick)) {
                    _selected = e.name;
                    if (ImGui::IsMouseDoubleClicked(0)) {
                        fs::path full = fs::path(_cwd) / e.name;
                        if (e.is_dir) {
                            _cwd = full.string();
                            _selected.clear();
                            refresh();
                            break;   // _entries invalidated
                        }
                        if (_mode == Mode::File) {
                            _result = full.string();
                            confirmed = true;
                        }
                    }
                }
            }
        }
        ImGui::EndChild();

        // ---- footer ----
        bool have_sel = !_selected.empty();
        if (_mode == Mode::Folder) {
            if (have_sel) {
                if (ImGui::Button("Select Highlighted Folder")) {
                    _result = (fs::path(_cwd) / _selected).string();
                    confirmed = true;
                }
                ImGui::SameLine();
            }
            if (ImGui::Button("Use This Folder")) {
                _result = _cwd;
                confirmed = true;
            }
        } else {
            ImGui::BeginDisabled(!have_sel);
            if (ImGui::Button("Select File") && have_sel) {
                _result = (fs::path(_cwd) / _selected).string();
                confirmed = true;
            }
            ImGui::EndDisabled();
        }
        ImGui::SameLine();
        if (ImGui::Button("Cancel")) _is_open = false;

        if (confirmed) _is_open = false;
        if (!_is_open) ImGui::CloseCurrentPopup();
        ImGui::EndPopup();
    }
    return confirmed;
}

}  // namespace gui
