#pragma once

// The glyphs. Which font the UI draws with, and how the CJK faces get onto
// disk.
//
// Three facts shape everything here:
//
//   * ImGui's built-in font is ASCII-only. German umlauts, French accents,
//     Turkish dotless i and Cyrillic were all broken before this file existed,
//     never mind Japanese -- so a Latin/Cyrillic face is EMBEDDED and always
//     loaded. That is assets/fonts/SpirulaUI-Regular.ttf, 59 KB.
//
//   * A full CJK face is 4-8 MB and there are four of them, because Han
//     unification means the shared codepoints render with different default
//     glyph forms per region. A Japanese reader shown the Simplified Chinese
//     face gets kanji in Chinese forms -- legible, and visibly wrong. Four
//     times 4-8 MB is too much to embed.
//
//   * But this program only ever writes ~600 characters per region: its own
//     translations. Subset to those and a regional face is ~110 KB, so all
//     FOUR are embedded (assets/fonts/SpirulaCJK-*.otf, 422 KB together).
//     That is what makes a default build render all thirteen languages, in
//     the right regional forms, with nothing to download -- which matters
//     most in the language picker, the one screen a user who cannot read the
//     current UI language has to be able to read.
//
// The full faces are still fetched on demand, but for a much smaller job than
// they used to do: dataset paths, file names and mask prompts are user data
// and can hold any character at all, and no subset can cover that. Until one
// is fetched, a folder called C:\写真\ may render as boxes; the UI itself
// never does.

#include "i18n/Message.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace gui {

// One regional CJK face, downloadable in full.
struct CjkFace {
    const char* id;       // "sc" -- also the SS_FONT_CJK value that bundles it
    const char* file;     // basename on disk
    const char* url;
    const char* sha256;   // verified after download and on every load
    uint64_t    bytes;
    const char* label;    // "Simplified Chinese", for the download prompt
};

// The embedded counterpart: the same region, cut down to this program's own
// vocabulary. Generated into app_generated/cjk_subsets.h at configure time.
struct CjkSubset {
    const char* id;
    const unsigned char* data;
    size_t size;
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

    // The full face for the current language, when it is not installed. NOT
    // an error state: the UI reads fine without it. It is what the language
    // menu offers so that CJK file names and typed prompts render too.
    const CjkFace* optional_face() const { return _optional; }

    // A build with SS_FONT_CJK=none has no fetch path, and should say so
    // rather than offer a download that will not happen.
    static bool fetch_enabled();

private:
    void rebuild();

    std::string _loaded_cjk;              // full face id currently in the atlas
    std::string _lead_cjk;                // region whose subset is merged first
    const CjkFace* _optional = nullptr;
    std::vector<char> _cjk_data;          // must outlive the atlas
    spirula::i18n::Lang _lang = spirula::i18n::Lang::en;
    bool _built = false;
    bool _dirty = true;
};

}  // namespace gui
