// Log.cpp -- see Log.h.

#include "sfm/core/Log.h"

#include "i18n/catalog/Sfm.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <mutex>

namespace msg = spirula::i18n::msg::sfm;

namespace sfm {
namespace slog {

namespace {

const spirula::i18n::Msg& tag_msg(Tag t) {
    switch (t) {
        case Tag::Run:     return msg::tag_run;
        case Tag::Extract: return msg::tag_extract;
        case Tag::Match:   return msg::tag_match;
        case Tag::Map:     return msg::tag_map;
        case Tag::Merge:   return msg::tag_merge;
        case Tag::Orient:  return msg::tag_orient;
        case Tag::Device:  return msg::tag_device;
    }
    return msg::tag_run;
}

constexpr int kNumTags = 7;

// The padded prefixes for one language, built once.
struct TagTable {
    spirula::i18n::Lang lang{};
    bool valid = false;
    std::string prefix[kNumTags];
};

std::mutex g_mu;
TagTable g_table;

const TagTable& table() {
    const spirula::i18n::Lang cur = spirula::i18n::current();
    if (g_table.valid && g_table.lang == cur) return g_table;
    int width = 0;
    for (int i = 0; i < kNumTags; i++)
        width = std::max(width, display_width(tag_msg((Tag)i).get()));
    for (int i = 0; i < kNumTags; i++) {
        const char* s = tag_msg((Tag)i).get();
        std::string p = "[";
        p += s;
        p += "]";
        p.append((size_t)std::max(0, width - display_width(s)), ' ');
        p += " ";
        g_table.prefix[i] = std::move(p);
    }
    g_table.lang = cur;
    g_table.valid = true;
    return g_table;
}

void write(std::FILE* f, Tag t, const std::string& body) {
    std::lock_guard<std::mutex> lk(g_mu);
    std::fprintf(f, "%s%s\n", table().prefix[(int)t].c_str(), body.c_str());
    // Unbuffered enough to interleave correctly with the other stream when a
    // parent process is reading both: the GUI shows one terminal, not two.
    std::fflush(f);
}

}  // namespace


int display_width(const char* s) {
    // UTF-8 decode, counting East Asian Wide / Fullwidth as two columns. The
    // ranges are the ones a tag can plausibly land in -- CJK, kana, hangul and
    // fullwidth forms -- not a full UAX #11 table.
    int w = 0;
    for (const unsigned char* p = (const unsigned char*)s; *p;) {
        uint32_t cp = *p;
        int len = 1;
        if (cp >= 0xF0)      { cp &= 0x07; len = 4; }
        else if (cp >= 0xE0) { cp &= 0x0F; len = 3; }
        else if (cp >= 0xC0) { cp &= 0x1F; len = 2; }
        for (int i = 1; i < len; i++) {
            if ((p[i] & 0xC0) != 0x80) { len = i; break; }
            cp = (cp << 6) | (p[i] & 0x3F);
        }
        p += len;
        const bool wide =
            (cp >= 0x1100 && cp <= 0x115F) ||    // hangul jamo
            (cp >= 0x2E80 && cp <= 0x303E) ||    // CJK radicals, punctuation
            (cp >= 0x3041 && cp <= 0x33FF) ||    // kana, hangul compat, CJK compat
            (cp >= 0x3400 && cp <= 0x4DBF) ||    // CJK ext A
            (cp >= 0x4E00 && cp <= 0x9FFF) ||    // CJK unified
            (cp >= 0xA000 && cp <= 0xA4CF) ||    // yi
            (cp >= 0xAC00 && cp <= 0xD7A3) ||    // hangul syllables
            (cp >= 0xF900 && cp <= 0xFAFF) ||    // CJK compat ideographs
            (cp >= 0xFE30 && cp <= 0xFE6F) ||    // CJK compat forms
            (cp >= 0xFF00 && cp <= 0xFF60) ||    // fullwidth forms
            (cp >= 0xFFE0 && cp <= 0xFFE6) ||
            (cp >= 0x20000 && cp <= 0x3FFFD);    // CJK ext B+
        w += wide ? 2 : 1;
    }
    return w;
}

std::string prefix(Tag t) {
    std::lock_guard<std::mutex> lk(g_mu);
    return table().prefix[(int)t];
}

std::string num(double v, int decimals) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", std::max(0, decimals), v);
    return buf;
}

void out(Tag t, const spirula::i18n::Msg& m,
         std::initializer_list<spirula::i18n::Arg> args) {
    write(stdout, t, spirula::i18n::format(m, args));
}

void err(Tag t, const spirula::i18n::Msg& m,
         std::initializer_list<spirula::i18n::Arg> args) {
    write(stderr, t, spirula::i18n::format(m, args));
}

void warn(Tag t, const spirula::i18n::Msg& m,
          std::initializer_list<spirula::i18n::Arg> args) {
    write(stderr, t, std::string(msg::word_warning.get()) + " " +
                         spirula::i18n::format(m, args));
}

void fail(Tag t, const spirula::i18n::Msg& m,
          std::initializer_list<spirula::i18n::Arg> args) {
    write(stderr, t, std::string(msg::word_error.get()) + " " +
                         spirula::i18n::format(m, args));
}

void out_raw(Tag t, const std::string& text) { write(stdout, t, text); }
void err_raw(Tag t, const std::string& text) { write(stderr, t, text); }

}  // namespace slog
}  // namespace sfm
