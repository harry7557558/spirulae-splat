// Fonts.cpp -- see Fonts.h.

#include "app/gui/Fonts.h"

#include "app/gui/AppPaths.h"
#include "core/Env.h"
#include "core/Sha256.h"

#include "app_generated/ui_font.h"    // kUiFont / kUiFontSize
#include "app_generated/cjk_faces.h"  // kCjkFaces, from assets/fonts/cjk_faces.txt

#include "imgui.h"

#include <cstdio>
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

namespace gui {

namespace {

// Base rasterization size. The old built-in font was a 13 px pixel font;
// an outline face at 13 px reads noticeably lighter, so the base goes up.
// DPI is NOT applied here -- GuiMain sets style.FontScaleMain, and since 1.92
// ImGui rasterizes on demand at the scaled size, so the atlas does not have to
// be rebuilt when the window moves to another monitor.
constexpr float kBaseSize = 16.0f;

// Every directory a face may already be in, most specific first.
std::vector<fs::path> search_dirs() {
    std::vector<fs::path> dirs;
    if (const char* d = spirula::env("FONT_DIR")) dirs.emplace_back(d);
    if (!exe_dir().empty()) dirs.push_back(fs::path(exe_dir()) / "fonts");
    dirs.push_back(fs::path(cache_dir()) / "fonts");
    return dirs;
}

bool verify(const fs::path& p, const CjkFace& f) {
    std::error_code ec;
    if (!fs::is_regular_file(p, ec)) return false;
    if (fs::file_size(p, ec) != f.bytes) return false;
    // A font is parsed by stb_truetype, which is exactly the sort of thing
    // that should not be handed unverified bytes off the network.
    return spirula::sha256_file(p.string()) == f.sha256;
}

}  // namespace

const CjkFace* cjk_face_for(spirula::i18n::Lang l) {
    using spirula::i18n::Lang;
    switch (l) {
        case Lang::zh_hans: return find_cjk_face("sc");
        case Lang::zh_hant: return find_cjk_face("tc");
        case Lang::ja:      return find_cjk_face("jp");
        case Lang::ko:      return find_cjk_face("kr");
        default:            return nullptr;
    }
}

const CjkFace* find_cjk_face(const std::string& id) {
    for (const CjkFace& f : kCjkFaces)
        if (id == f.id) return &f;
    return nullptr;
}

std::string cjk_face_path(const CjkFace& f) {
    for (const fs::path& d : search_dirs()) {
        const fs::path p = d / f.file;
        if (verify(p, f)) return p.string();
    }
    return {};
}

std::string cjk_face_download_path(const CjkFace& f) {
    return (fs::path(cache_dir()) / "fonts" / f.file).string();
}

bool FontSet::fetch_enabled() {
#ifdef SS_FONT_CJK_NONE
    return false;
#else
    return true;
#endif
}

void FontSet::ensure() {
    // Deciding what to load hashes up to 8 MB, so it happens when the answer
    // can have changed -- a language switch or a finished download -- and not
    // every frame.
    const spirula::i18n::Lang lang = spirula::i18n::current();
    if (_built && !_dirty && lang == _lang) return;
    _lang = lang;
    _dirty = false;

    const CjkFace* want = cjk_face_for(lang);

    // A face already on disk is worth loading whatever the UI language is:
    // dataset paths and file names are not translated, and a Japanese folder
    // name in an English UI should not be a row of boxes. The current
    // language's face wins when several are present, because Han unification
    // makes that the one whose glyph forms are right.
    std::string chosen, chosen_id;
    if (want) chosen = cjk_face_path(*want);
    if (!chosen.empty()) {
        chosen_id = want->id;
        _missing = nullptr;
    } else {
        _missing = want;                  // null unless the language needs one
        for (const CjkFace& f : kCjkFaces) {
            chosen = cjk_face_path(f);
            if (!chosen.empty()) { chosen_id = f.id; break; }
        }
    }

    if (_built && chosen_id == _loaded_cjk) return;
    _loaded_cjk = chosen_id;

    _cjk_data.clear();
    if (!chosen.empty()) {
        std::ifstream fin(chosen, std::ios::binary | std::ios::ate);
        if (fin) {
            const std::streamoff n = fin.tellg();
            fin.seekg(0);
            _cjk_data.resize((size_t)n);
            fin.read(_cjk_data.data(), n);
            if (!fin) _cjk_data.clear();
        }
        if (_cjk_data.empty())
            std::fprintf(stderr, "warning: could not read %s\n", chosen.c_str());
    }
    rebuild();
    _built = true;
}

void FontSet::rebuild() {
    // ClearFonts(), not Clear(): it drops the font inputs and the cached
    // glyphs but keeps the texture, which is what makes swapping a face
    // between frames safe. Clear() also destroys the texture and is documented
    // as not being callable at this point.
    ImFontAtlas* atlas = ImGui::GetIO().Fonts;
    atlas->ClearFonts();

    // Latin/Cyrillic first: with merged sources ImGui takes the glyph from the
    // first source that has it, and this is the face the UI was designed
    // around. The CJK face carries its own Latin, which would otherwise win
    // and quietly change every English string's appearance.
    ImFontConfig base;
    base.FontDataOwnedByAtlas = false;    // static storage, outlives the atlas
    std::snprintf(base.Name, sizeof base.Name, "SpirulaUI");
    atlas->AddFontFromMemoryTTF(const_cast<unsigned char*>(kUiFont),
                                (int)kUiFontSize, kBaseSize, &base);

    if (!_cjk_data.empty()) {
        ImFontConfig cjk;
        cjk.MergeMode = true;
        cjk.FontDataOwnedByAtlas = false; // owned by _cjk_data, which outlives it
        std::snprintf(cjk.Name, sizeof cjk.Name, "NotoSansCJK-%s",
                      _loaded_cjk.c_str());
        atlas->AddFontFromMemoryTTF(_cjk_data.data(), (int)_cjk_data.size(),
                                    kBaseSize, &cjk);
    }
}

}  // namespace gui
