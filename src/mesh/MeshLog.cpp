// MeshLog.cpp -- see MeshLog.h.

#include "mesh/MeshLog.h"

#include <algorithm>
#include <cstdio>
#include <mutex>

namespace mmsg = spirula::i18n::msg::mesh;

namespace meshing {
namespace mlog {

namespace {

std::mutex g_mu;

void write(std::FILE* f, Stage s, const std::string& body) {
    std::lock_guard<std::mutex> lk(g_mu);
    std::fprintf(f, "%s%s%s\n", tag_prefix().c_str(), stage_label(s).c_str(),
                 body.c_str());
    // Flushed per line so the two streams interleave correctly when a parent
    // process is reading both: the GUI shows one terminal, not two.
    std::fflush(f);
}

}  // namespace

void out(Stage s, const spirula::i18n::Msg& m,
         std::initializer_list<spirula::i18n::Arg> args) {
    write(stdout, s, spirula::i18n::format(m, args));
}

void err(Stage s, const spirula::i18n::Msg& m,
         std::initializer_list<spirula::i18n::Arg> args) {
    write(stderr, s, spirula::i18n::format(m, args));
}

void warn(Stage s, const spirula::i18n::Msg& m,
          std::initializer_list<spirula::i18n::Arg> args) {
    write(stderr, s, std::string(mmsg::word_warning.get()) + " " +
                         spirula::i18n::format(m, args));
}

void fail(Stage s, const spirula::i18n::Msg& m,
          std::initializer_list<spirula::i18n::Arg> args) {
    write(stderr, s, std::string(mmsg::word_error.get()) + " " +
                         spirula::i18n::format(m, args));
}

void out_raw(Stage s, const std::string& text) { write(stdout, s, text); }
void err_raw(Stage s, const std::string& text) { write(stderr, s, text); }

void fail_raw(const std::string& text) {
    std::lock_guard<std::mutex> lk(g_mu);
    std::fprintf(stderr, "%s%s %s\n", tag_prefix().c_str(),
                 mmsg::word_error.get(), text.c_str());
    std::fflush(stderr);
}

std::string num(double v, int decimals) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*f", std::max(0, decimals), v);
    return buf;
}

}  // namespace mlog
}  // namespace meshing
