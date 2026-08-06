#pragma once

// The ONE place the language set is written.
//
// Everything else derives from this table: the `Lang` enum, the number of
// slots in a `Msg`, the tag macros in BeginCatalog.h, the language picker, the
// `--lang` parser and the locale detector. Adding a row here breaks the build
// on every message that does not have a translation for it -- which is the
// point. A missing translation is a compile error, not a silent fallback.
//
//   X(id, code, native, english)
//     id       the C++ identifier: `Lang::ja`, and the SS_DEFAULT_LANG value
//     code     what `--lang` / SS_LANG accept and what an OS locale maps to
//     native   the language's own name, for the picker -- NEVER translated
//     english  the English name, for `--help` and diagnostics
//
// Underscores rather than hyphens in the ids (`zh_hans`, not `zh-Hans`)
// because they are identifiers; `parse_lang()` accepts either spelling.
//
// The order is the picker's order: English first, then by how much the UI
// changes shape (CJK, then the Latin-script locales, then Cyrillic).

#define SS_LANGUAGES(X)                                                       \
    X(en,      "en",      "English",     "English")                           \
    X(ja,      "ja",      "日本語",       "Japanese")                          \
    X(zh_hans, "zh-Hans", "简体中文",     "Chinese (Simplified)")               \
    X(zh_hant, "zh-Hant", "繁體中文",     "Chinese (Traditional)")              \
    X(ko,      "ko",      "한국어",       "Korean")                            \
    X(de,      "de",      "Deutsch",     "German")                            \
    X(fr,      "fr",      "Français",    "French")                            \
    X(es,      "es",      "Español",     "Spanish")                           \
    X(pt,      "pt",      "Português",   "Portuguese")                        \
    X(it,      "it",      "Italiano",    "Italian")                           \
    X(nl,      "nl",      "Nederlands",  "Dutch")                             \
    X(ru,      "ru",      "Русский",     "Russian")                           \
    X(tr,      "tr",      "Türkçe",      "Turkish")

// What the language menu is labelled with in the menu bar.
//
// An icon rather than the word "Language", because of who needs that menu:
// someone looking at a UI in a language they cannot read, who therefore
// cannot read the menu that fixes it either. 文A is the widely used mark for
// this -- one Han character next to one Latin letter, which reads as "scripts"
// without belonging to either script. The menu draws it followed by the
// current language's own name, so the menu also says what is selected.
//
// Both characters have to be in the embedded fonts, which is why this lives
// here: tools/make_ui_font.py reads it from this file, so changing it
// regenerates the subsets rather than silently producing a box.
#define SS_LANG_MENU_ICON "文A"

// Which languages need a CJK font to render at all. Used by the font loader
// to decide whether to fetch a face, and by the CMake consistency check that
// refuses SS_DEFAULT_LANG=ja with SS_FONT_CJK=none.
#define SS_LANGUAGES_CJK(X)                                                   \
    X(ja) X(zh_hans) X(zh_hant) X(ko)
