// Locale.cpp -- see Locale.h.

#include "i18n/Locale.h"

#include "core/Env.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#elif defined(__APPLE__)
#  include <CoreFoundation/CoreFoundation.h>
#endif

namespace spirula {
namespace i18n {

namespace detail {
std::atomic<Lang> g_current{kDefaultLang};
}

void set_current(Lang l) {
    detail::g_current.store(l, std::memory_order_relaxed);
}

const char* code(Lang l) {
    switch (l) {
#define X(id, c, native, english) case Lang::id: return c;
        SS_LANGUAGES(X)
#undef X
    }
    return "en";
}

const char* native_name(Lang l) {
    switch (l) {
#define X(id, c, native, english) case Lang::id: return native;
        SS_LANGUAGES(X)
#undef X
    }
    return "English";
}

const char* english_name(Lang l) {
    switch (l) {
#define X(id, c, native, english) case Lang::id: return english;
        SS_LANGUAGES(X)
#undef X
    }
    return "English";
}

bool needs_cjk_font(Lang l) {
    switch (l) {
#define X(id) case Lang::id: return true;
        SS_LANGUAGES_CJK(X)
#undef X
        default: return false;
    }
}

namespace {

char lower(char c) {
    return (c >= 'A' && c <= 'Z') ? char(c - 'A' + 'a') : c;
}

// "zh_Hans_CN.UTF-8@modifier" -> {"zh", "hans", "cn"}, lowercased, at most
// three subtags (nothing we match on needs a fourth).
struct Subtags {
    std::string t[3];
    int n = 0;
};

Subtags split_locale(const char* s) {
    Subtags out;
    std::string cur;
    for (const char* p = s; *p; p++) {
        // The encoding and the modifier say nothing about the language.
        if (*p == '.' || *p == '@') break;
        if (*p == '-' || *p == '_') {
            if (!cur.empty() && out.n < 3) out.t[out.n++] = cur;
            cur.clear();
            continue;
        }
        cur += lower(*p);
    }
    if (!cur.empty() && out.n < 3) out.t[out.n++] = cur;
    return out;
}

// zh needs a script, and an OS will express it as a script subtag, a region,
// or -- on older Windows -- as the "CHS"/"CHT" pseudo-scripts.
bool resolve_chinese(const Subtags& st, Lang* out) {
    for (int i = 1; i < st.n; i++) {
        const std::string& t = st.t[i];
        if (t == "hans" || t == "chs") { *out = Lang::zh_hans; return true; }
        if (t == "hant" || t == "cht") { *out = Lang::zh_hant; return true; }
        // Region. Macau and Hong Kong are traditional; Singapore and Malaysia
        // follow the mainland.
        if (t == "cn" || t == "sg" || t == "my") { *out = Lang::zh_hans; return true; }
        if (t == "tw" || t == "hk" || t == "mo") { *out = Lang::zh_hant; return true; }
    }
    // A bare "zh" with no script and no region. A Traditional build should
    // stay Traditional; everything else takes the more common script.
    *out = (kDefaultLang == Lang::zh_hant) ? Lang::zh_hant : Lang::zh_hans;
    return true;
}

}  // namespace

bool parse_lang(const char* s, Lang* out) {
    if (!s || !*s) return false;

    const Subtags st = split_locale(s);
    if (st.n == 0) return false;

    // "C" and "POSIX" mean "the user expressed no preference", which is not
    // the same as "the user wants English" -- fall through to the next step of
    // the chain rather than pinning en here.
    if (st.t[0] == "c" || st.t[0] == "posix") return false;

    if (st.t[0] == "zh") return resolve_chinese(st, out);

    // Our own canonical codes, and the identifier spelling ("zh_hans") that
    // SS_DEFAULT_LANG and the settings file use.
    std::string joined = st.t[0];
    for (int i = 1; i < st.n; i++) joined += "-" + st.t[i];

#define X(id, c, native, english)                                             \
    {                                                                         \
        std::string canon(c);                                                 \
        for (auto& ch : canon) ch = lower(ch);                                \
        if (joined == canon || st.t[0] == canon) { *out = Lang::id; return true; } \
    }
    SS_LANGUAGES(X)
#undef X

    return false;
}

Lang detect_os_lang() {
    Lang l = kDefaultLang;

#if defined(_WIN32)
    // The user's UI language, e.g. L"zh-Hans-CN". Present since Vista.
    wchar_t buf[LOCALE_NAME_MAX_LENGTH] = {0};
    if (GetUserDefaultLocaleName(buf, LOCALE_NAME_MAX_LENGTH) > 0) {
        // Locale names are ASCII; a narrowing copy is exact.
        char narrow[LOCALE_NAME_MAX_LENGTH] = {0};
        size_t n = 0;
        for (; n + 1 < sizeof narrow && buf[n]; n++)
            narrow[n] = (buf[n] < 128) ? char(buf[n]) : '?';
        narrow[n] = '\0';
        if (parse_lang(narrow, &l)) return l;
    }
#elif defined(__APPLE__)
    // POSIX env vars are usually absent for an app launched from Finder, so
    // ask CoreFoundation. (CFLocale is C; no Objective-C is involved.)
    if (CFLocaleRef loc = CFLocaleCopyCurrent()) {
        char buf[128] = {0};
        const bool ok = CFStringGetCString(CFLocaleGetIdentifier(loc), buf,
                                           sizeof buf, kCFStringEncodingUTF8);
        CFRelease(loc);
        if (ok && parse_lang(buf, &l)) return l;
    }
#endif

    // POSIX, and the fallback everywhere: LC_ALL beats LC_MESSAGES beats LANG.
    for (const char* var : {"LC_ALL", "LC_MESSAGES", "LANG"})
        if (const char* v = std::getenv(var))
            if (parse_lang(v, &l)) return l;

    return kDefaultLang;
}

Lang init(const char* cli, const char* saved) {
    Lang l = kDefaultLang;

    auto try_set = [&](const char* s, const char* origin) {
        if (!s || !*s) return false;
        if (parse_lang(s, &l)) return true;
        std::fprintf(stderr, "warning: unknown language '%s' (%s); using %s\n"
                             "%s",
                     s, origin, code(kDefaultLang), language_list().c_str());
        return false;
    };

    if (try_set(cli, "--lang")) { set_current(l); return l; }
    if (try_set(env("LANG"), "SS_LANG")) { set_current(l); return l; }
    if (try_set(saved, "settings")) { set_current(l); return l; }

    l = detect_os_lang();
    set_current(l);
    return l;
}

namespace {
const char* g_lang_arg = nullptr;
}

const char* lang_arg() { return g_lang_arg; }

const char* take_lang_arg(int* argc, char** argv) {
    const char* value = nullptr;
    int out = 1;
    for (int i = 1; i < *argc; i++) {
        const char* a = argv[i];
        if (std::strcmp(a, "--lang") == 0 && i + 1 < *argc) {
            value = argv[++i];
            continue;
        }
        if (std::strncmp(a, "--lang=", 7) == 0) {
            value = a + 7;
            continue;
        }
        argv[out++] = argv[i];
    }
    argv[out] = nullptr;
    *argc = out;
    g_lang_arg = value;
    return value;
}

std::string language_list() {
    std::string s;
#define X(id, c, native, english)                                             \
    {                                                                         \
        std::string row = std::string("    ") + c;                            \
        row.resize(14, ' ');                                                  \
        s += row; s += english; s += "\n";                                    \
    }
    SS_LANGUAGES(X)
#undef X
    return s;
}

// ---------------------------------------------------------------------------
// Positional substitution -- see the rule in Message.h.
// ---------------------------------------------------------------------------

Arg::Arg(int v)                : s(std::to_string(v)) {}
Arg::Arg(long v)               : s(std::to_string(v)) {}
Arg::Arg(long long v)          : s(std::to_string(v)) {}
Arg::Arg(unsigned v)           : s(std::to_string(v)) {}
Arg::Arg(unsigned long v)      : s(std::to_string(v)) {}
Arg::Arg(unsigned long long v) : s(std::to_string(v)) {}

Arg::Arg(double v) {
    // std::to_string(double) is always %f, so 0.5 becomes "0.500000" and 1e-7
    // becomes "0.000000". %g is what a UI wants.
    char buf[32];
    std::snprintf(buf, sizeof buf, "%g", v);
    s = buf;
}

std::string format(const char* pattern, std::initializer_list<Arg> args) {
    std::string out;
    if (!pattern) return out;
    for (const char* p = pattern; *p;) {
        if (p[0] == '{' && p[1] == '{') { out += '{'; p += 2; continue; }
        if (p[0] == '}' && p[1] == '}') { out += '}'; p += 2; continue; }
        if (*p != '{') { out += *p++; continue; }

        const char* q = p + 1;
        size_t idx = 0;
        bool digits = false;
        while (*q >= '0' && *q <= '9') { idx = idx * 10 + size_t(*q - '0'); q++; digits = true; }
        if (!digits || *q != '}' || idx >= args.size()) {
            // Not a placeholder, or one with no argument. Copy it through so
            // the mistake is visible as `{2}` in the UI instead of vanishing.
            out += *p++;
            continue;
        }
        out += (args.begin() + idx)->s;
        p = q + 1;
    }
    return out;
}

}  // namespace i18n
}  // namespace spirula
