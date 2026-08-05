#pragma once

// Which language the UI is in, and how that is decided.
//
// Resolution order, first hit wins (§6.7 of the localization plan):
//
//   1. `--lang <code>` on the command line
//   2. the SS_LANG environment variable
//   3. the settings file (the GUI's language picker writes it here)
//   4. the OS locale
//   5. SS_DEFAULT_LANG, the compile-time default
//
// Step 5 is a *language*, not "English": a build shipped as
// -DSS_DEFAULT_LANG=ja must come up Japanese in a container with no LANG set,
// which is exactly the case a hard-coded English fallback gets wrong.

#include "i18n/Message.h"

// The compile-time default is an identifier, not a string, so a typo in
// -DSS_DEFAULT_LANG is a compile error naming the bad language rather than a
// binary that silently comes up English.
#ifndef SS_DEFAULT_LANG
#define SS_DEFAULT_LANG en
#endif

namespace spirula {
namespace i18n {

inline constexpr Lang kDefaultLang = Lang::SS_DEFAULT_LANG;

// "zh-Hans" -- the canonical code, what `--lang` prints and what the settings
// file stores.
const char* code(Lang l);
// "简体中文" -- the language's own name, for the picker. Never translated.
const char* native_name(Lang l);
// "Chinese (Simplified)" -- for `--help` and diagnostics, always English.
const char* english_name(Lang l);

// Accepts our own codes and anything an OS is likely to hand us:
// "ja", "ja_JP.UTF-8", "zh-Hans-CN", "zh_TW", "pt-BR", "zh_hant". Returns
// false and leaves `out` alone if nothing matches; "C" and "POSIX" do not
// match (they mean "no preference", not "English").
bool parse_lang(const char* s, Lang* out);

// The OS's idea of the user's language, or kDefaultLang if it has none we
// recognize.
Lang detect_os_lang();

// Runs the resolution chain above and calls set_current(). `cli` and `saved`
// may be null or empty. Returns what it settled on.
Lang init(const char* cli, const char* saved);

// True for the languages that cannot be rendered from the Latin font alone.
bool needs_cjk_font(Lang l);

// Scans argv for `--lang <code>` / `--lang=<code>`, REMOVES the arguments it
// consumes (so no per-tool parser has to know about the flag), and returns the
// value, or nullptr if the flag was absent. An unparseable value is left for
// init() to reject with a message listing the codes.
const char* take_lang_arg(int* argc, char** argv);

// What take_lang_arg() found, for a component that learns the saved preference
// only later: the GUI reads its settings file long after main() has parsed
// argv, and has to re-run the chain without --lang losing to a stored value.
const char* lang_arg();

// One line per language, for `--help`.
std::string language_list();

}  // namespace i18n
}  // namespace spirula
