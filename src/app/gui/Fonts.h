#pragma once

// The glyphs. Which font the UI draws with, and how the CJK faces get onto
// disk.
//
// Two facts shape everything here:
//
//   * ImGui's built-in font is ASCII-only. German umlauts, French accents,
//     Turkish dotless i and Cyrillic were all broken before this file existed,
//     never mind Japanese -- so a Latin/Cyrillic face is EMBEDDED and always
//     loaded. That is assets/fonts/SpirulaUI-Regular.ttf, 59 KB, and it is
//     what makes a default build render twelve of the thirteen languages.
//
//   * A CJK face is 4-8 MB and there are four of them, because Han
//     unification means the shared codepoints render with different default
//     glyph forms per region. A Japanese reader shown the Simplified Chinese
//     face gets kanji in Chinese forms -- legible, and visibly wrong. So the
//     face is per-region, and per-region times 4 is too much to embed. It is
//     downloaded on demand instead (SS_FONT_CJK=fetch, the default) or shipped
//     beside the executable by a regional build (SS_FONT_CJK=sc|tc|jp|kr|all).
//
// A face that is already on disk is loaded whatever the UI language is, not
// only for CJK locales: an English UI still has to draw a dataset path like
// C:\写真\ without turning it into boxes.

#include "i18n/Message.h"

#include <cstdint>
#include <string>
#include <vector>

namespace gui {

// One regional CJK face.
struct CjkFace {
    const char* id;       // "sc" -- also the SS_FONT_CJK value that embeds it
    const char* file;     // basename on disk
    const char* url;
    const char* sha256;   // verified after download and on every load
    uint64_t    bytes;
    const char* label;    // "Simplified Chinese", for the download prompt
};

// The face a language needs, or null for the Latin-script languages.
const CjkFace* cjk_face_for(spirula::i18n::Lang l);
const CjkFace* find_cjk_face(const std::string& id);

// Where that face would live, and whether a verified copy is there. Searched
// in order: $SS_FONT_DIR, <exe dir>/fonts (a regional build ships it there),
// then <cache dir>/fonts (what a download writes). Empty when not found.
std::string cjk_face_path(const CjkFace& f);
std::string cjk_face_download_path(const CjkFace& f);

// The atlas.
//
// Call ensure() once per frame, before ImGui::NewFrame(). It returns
// immediately unless the language changed or invalidate() was called, because
// deciding what to load means hashing a multi-megabyte file and that is not a
// per-frame job.
class FontSet {
public:
    void ensure();

    // Re-examine the disk on the next ensure(). Call when a font download
    // finishes; nothing else changes the answer.
    void invalidate() { _dirty = true; }

    // The face the current language needs and does not have -- i.e. the UI is
    // about to render tofu. Null when all is well. The language picker uses
    // this to offer the download.
    const CjkFace* missing_face() const { return _missing; }

    // A build with SS_FONT_CJK=none has no fetch path, and should say so
    // rather than offer a download that will not happen.
    static bool fetch_enabled();

private:
    void rebuild();

    std::string _loaded_cjk;              // face id currently in the atlas
    const CjkFace* _missing = nullptr;
    std::vector<char> _cjk_data;          // must outlive the atlas
    spirula::i18n::Lang _lang = spirula::i18n::Lang::en;
    bool _built = false;
    bool _dirty = true;
};

}  // namespace gui
