// Native replacement for slang/build_spirv.py + spirulae_splat/embed_spirv.py.
//
// A single self-contained C++17 host tool (no dependencies beyond the standard
// library) so the Vulkan (SSPLAT_BACKEND=vulkan) build no longer needs Python.
// CMake compiles this once at configure time (try_compile + COPY_FILE) and uses
// it in two modes:
//
//   discover -I <dir>... <source.slang>...
//       Scan the compute shaders and print, one blob per line to stdout:
//           <name>\t<source>\t<entry>\t<defines>\t<dep1> <dep2> ...
//       CMake turns each line into its own slangc custom command (one Ninja
//       edge per blob, so the build's -j governs how many slangc run at once,
//       and Ninja prints "[n/m] SPIR-V <name>" as each finishes). <defines> and
//       <deps> are space-separated (deps = the source's #include closure).
//
//   embed <out.cpp> --list <listfile>
//       Read the compiled .spv blobs named in <listfile>, drop the .noint64
//       variants whose base does not actually declare the Int64 capability
//       (reproducing build_spirv.py's phase-2 gate), verify no capability leaks,
//       and emit the C++ translation unit consumed by VulkanPipelines.cpp.
//
// Blob naming, feature variants, and the capability audit mirror build_spirv.py
// exactly; see that file's history and backend/vulkan/README.md for the why.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

// ---- feature variants (canonical suffix order must match kFeatureSuffixes in
// csrc/backend/vulkan/VulkanPipelines.cpp) --------------------------------
struct Feature { const char* suffix; const char* define; };
const Feature ATOMIC{".atomicadd", "-DSSPLAT_NATIVE_F32_ATOMIC_ADD"};
const Feature INT8{".int8", "-DSSPLAT_NATIVE_INT8"};
const Feature NOINT64{".noint64", "-DSSPLAT_EMULATE_INT64"};

// SPIR-V capability operands (OpCapability = opcode 17).
constexpr uint32_t CAP_INT64 = 11;
constexpr uint32_t CAP_INT8 = 39;
constexpr uint32_t CAP_8BIT_STORAGE = 4437;  // StorageBuffer8BitAccess

bool is_ident(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
           (c >= '0' && c <= '9') || c == '_';
}

std::string read_file(const std::string& path, bool* ok = nullptr) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { if (ok) *ok = false; return {}; }
    if (ok) *ok = true;
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::string basename_of(const std::string& p) {
    size_t s = p.find_last_of("/\\");
    return s == std::string::npos ? p : p.substr(s + 1);
}
std::string dirname_of(const std::string& p) {
    size_t s = p.find_last_of("/\\");
    return s == std::string::npos ? std::string(".") : p.substr(0, s);
}
bool file_exists(const std::string& p) {
    std::ifstream f(p);
    return static_cast<bool>(f);
}

// ---- shader source scanning (mirrors build_spirv.py regexes) -------------

// `[shader("compute")]` ... `void <name>(`, allowing attribute/comment lines in
// between: take the first `void <name>(` after each compute marker.
bool parse_void_name(const std::string& s, size_t from, std::string* name) {
    const std::string kw = "void";
    while (true) {
        size_t v = s.find(kw, from);
        if (v == std::string::npos) return false;
        // must be a whole word
        bool lok = (v == 0) || !is_ident(s[v - 1]);
        size_t p = v + kw.size();
        if (!lok || p >= s.size() || is_ident(s[p])) { from = v + 1; continue; }
        while (p < s.size() && (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' ||
                                s[p] == '\r')) p++;
        size_t ns = p;
        while (p < s.size() && is_ident(s[p])) p++;
        if (p == ns) { from = v + 1; continue; }  // no identifier
        std::string id = s.substr(ns, p - ns);
        while (p < s.size() && (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' ||
                                s[p] == '\r')) p++;
        if (p < s.size() && s[p] == '(') { *name = id; return true; }
        from = v + 1;
    }
}

// Discover entry points: `[shader("compute")]`-marked functions and line-initial
// `DEF_<macro>(<entry>, ...)` instantiations (the literal macro parameter NAME
// is skipped, matching build_spirv.py).
std::vector<std::string> find_entries(const std::string& src) {
    std::vector<std::string> out;
    std::set<std::string> seen;
    auto add = [&](const std::string& n) {
        if (n != "NAME" && seen.insert(n).second) out.push_back(n);
    };
    // [shader("compute")] ... void NAME(
    const std::string marker = "[shader(\"compute\")]";
    size_t pos = 0;
    std::vector<std::string> shader_names;
    while ((pos = src.find(marker, pos)) != std::string::npos) {
        pos += marker.size();
        std::string n;
        if (parse_void_name(src, pos, &n)) shader_names.push_back(n);
    }
    // line-initial DEF_...(entry, ...)
    std::vector<std::string> macro_names;
    size_t line_start = 0;
    while (line_start <= src.size()) {
        size_t nl = src.find('\n', line_start);
        size_t end = (nl == std::string::npos) ? src.size() : nl;
        // line begins with "DEF_"
        if (end - line_start > 4 && src.compare(line_start, 4, "DEF_") == 0) {
            size_t p = line_start + 4;
            while (p < end && is_ident(src[p])) p++;  // rest of macro name
            if (p < end && src[p] == '(') {
                p++;
                while (p < end && (src[p] == ' ' || src[p] == '\t')) p++;
                size_t ns = p;
                while (p < end && is_ident(src[p])) p++;
                if (p > ns) macro_names.push_back(src.substr(ns, p - ns));
            }
        }
        if (nl == std::string::npos) break;
        line_start = nl + 1;
    }
    for (auto& n : shader_names) add(n);
    for (auto& n : macro_names) add(n);
    return out;
}

// Extract `#include "..."` targets (leading whitespace allowed).
std::vector<std::string> parse_includes(const std::string& src) {
    std::vector<std::string> out;
    size_t line_start = 0;
    while (line_start <= src.size()) {
        size_t nl = src.find('\n', line_start);
        size_t end = (nl == std::string::npos) ? src.size() : nl;
        size_t p = line_start;
        while (p < end && (src[p] == ' ' || src[p] == '\t')) p++;
        const std::string inc = "#include";
        if (end - p >= inc.size() && src.compare(p, inc.size(), inc) == 0) {
            p += inc.size();
            while (p < end && (src[p] == ' ' || src[p] == '\t')) p++;
            if (p < end && src[p] == '"') {
                size_t q = src.find('"', p + 1);
                if (q != std::string::npos && q <= end)
                    out.push_back(src.substr(p + 1, q - p - 1));
            }
        }
        if (nl == std::string::npos) break;
        line_start = nl + 1;
    }
    return out;
}

std::string resolve_include(const std::string& inc, const std::string& from_dir,
                            const std::vector<std::string>& incdirs) {
    std::vector<std::string> bases;
    bases.push_back(from_dir);
    for (auto& d : incdirs) bases.push_back(d);
    for (auto& b : bases) {
        std::string cand = b + "/" + inc;
        if (file_exists(cand)) return cand;
    }
    return {};
}

// Transitive #include closure (the source itself first), resolved against the
// including file's directory and the include dirs.
void include_closure(const std::string& path,
                     const std::vector<std::string>& incdirs,
                     std::vector<std::string>& seen) {
    if (std::find(seen.begin(), seen.end(), path) != seen.end()) return;
    bool ok = false;
    std::string src = read_file(path, &ok);
    if (!ok) return;
    seen.push_back(path);
    for (auto& inc : parse_includes(src)) {
        std::string r = resolve_include(inc, dirname_of(path), incdirs);
        if (!r.empty()) include_closure(r, incdirs, seen);
    }
}

// True when a whole-word `token` immediately followed by `(` (whitespace
// permitted) occurs in `text`.
bool calls_token(const std::string& text, const std::string& token) {
    size_t pos = 0;
    while ((pos = text.find(token, pos)) != std::string::npos) {
        bool lok = (pos == 0) || !is_ident(text[pos - 1]);
        size_t p = pos + token.size();
        while (p < text.size() && (text[p] == ' ' || text[p] == '\t' ||
                                   text[p] == '\n' || text[p] == '\r')) p++;
        if (lok && p < text.size() && text[p] == '(') return true;
        pos += token.size();
    }
    return false;
}

// Does the source's include closure call any of `tokens` from `header` (the
// header's own definitions do not count)?
bool closure_calls(const std::vector<std::string>& closure,
                   const std::string& header_base,
                   const std::vector<std::string>& tokens) {
    bool has_header = false;
    for (auto& f : closure)
        if (basename_of(f) == header_base) has_header = true;
    if (!has_header) return false;
    for (auto& f : closure) {
        if (basename_of(f) == header_base) continue;
        std::string t = read_file(f);
        for (auto& tok : tokens)
            if (calls_token(t, tok)) return true;
    }
    return false;
}

bool closure_has(const std::vector<std::string>& closure,
                 const std::string& header_base) {
    for (auto& f : closure)
        if (basename_of(f) == header_base) return true;
    return false;
}

// All subsets of `feats` in canonical order, empty subset first.
std::vector<std::pair<std::string, std::vector<std::string>>>
variant_subsets(const std::vector<Feature>& feats) {
    std::vector<std::pair<std::string, std::vector<std::string>>> out;
    size_t n = feats.size();
    // Ordered by subset size then canonical index order (matches build_spirv.py).
    for (size_t r = 0; r <= n; r++) {
        // iterate combinations of size r preserving index order
        std::vector<size_t> idx(r);
        for (size_t i = 0; i < r; i++) idx[i] = i;
        if (r == 0) { out.push_back({"", {}}); continue; }
        while (true) {
            std::string suffix;
            std::vector<std::string> defs;
            for (size_t i = 0; i < r; i++) {
                suffix += feats[idx[i]].suffix;
                defs.push_back(feats[idx[i]].define);
            }
            out.push_back({suffix, defs});
            // next combination
            long i = (long)r - 1;
            while (i >= 0 && idx[i] == n - r + (size_t)i) i--;
            if (i < 0) break;
            idx[i]++;
            for (size_t j = i + 1; j < r; j++) idx[j] = idx[j - 1] + 1;
        }
    }
    return out;
}

// ---- SPIR-V capability inspection ---------------------------------------
std::set<uint32_t> spirv_capabilities(const std::string& bytes) {
    std::set<uint32_t> caps;
    if (bytes.size() < 24) return caps;
    auto w = [&](size_t off) {
        return (uint32_t)(uint8_t)bytes[off] |
               ((uint32_t)(uint8_t)bytes[off + 1] << 8) |
               ((uint32_t)(uint8_t)bytes[off + 2] << 16) |
               ((uint32_t)(uint8_t)bytes[off + 3] << 24);
    };
    if (w(0) != 0x07230203u) return caps;
    size_t i = 20;
    while (i + 4 <= bytes.size()) {
        uint32_t word = w(i);
        uint32_t opcode = word & 0xffff, count = word >> 16;
        if (count == 0) break;
        if (opcode == 17 && i + 8 <= bytes.size()) {
            caps.insert(w(i + 4));  // OpCapability
        } else if (opcode != 17) {
            break;  // capabilities lead the module
        }
        i += (size_t)count * 4;
    }
    return caps;
}

// ---- discover mode -------------------------------------------------------
int run_discover(const std::vector<std::string>& args) {
    std::vector<std::string> incdirs, sources;
    for (size_t i = 0; i < args.size(); i++) {
        if (args[i] == "-I" && i + 1 < args.size())
            incdirs.push_back(args[++i]);            // "-I" "dir"
        else if (args[i].rfind("-I", 0) == 0)
            incdirs.push_back(args[i].substr(2));    // "-Idir"
        else
            sources.push_back(args[i]);
    }

    // Files pulled in as headers by some source; an entry-less source is only
    // worth a warning when nothing includes it.
    std::set<std::string> included;
    std::map<std::string, std::vector<std::string>> closures;
    for (auto& s : sources) {
        std::vector<std::string> cl;
        include_closure(s, incdirs, cl);
        closures[s] = cl;
        for (auto& f : cl) if (f != s) included.insert(f);
    }

    std::ostringstream out;
    for (auto& source : sources) {
        std::string stem = basename_of(source);
        if (stem.size() > 6 && stem.compare(stem.size() - 6, 6, ".slang") == 0)
            stem = stem.substr(0, stem.size() - 6);
        const auto& closure = closures[source];
        std::vector<std::string> entries = find_entries(read_file(source));
        if (entries.empty() && !included.count(source))
            std::fprintf(stderr, "warning: no compute entry points in %s\n",
                         source.c_str());

        std::vector<Feature> feats;
        if (closure_calls(closure, "atomic_float.slang", {"atomic_add_f32"}))
            feats.push_back(ATOMIC);
        if (closure_calls(closure, "int8_compat.slang",
                          {"u8_load", "s8_load", "u8_store"}))
            feats.push_back(INT8);
        // noint64 candidate: the source can route 64-bit integers through
        // int64_compat.slang. Over-approximates build_spirv.py's compile-then-
        // inspect gate; the embed step drops variants whose base has no Int64.
        bool noint64_cand = closure_has(closure, "int64_compat.slang");

        std::string deps;
        for (size_t i = 0; i < closure.size(); i++)
            deps += (i ? " " : "") + closure[i];

        auto emit = [&](const std::string& entry, const std::string& suffix,
                        const std::vector<std::string>& defs) {
            std::string name = stem + "." + entry + suffix;
            std::string defstr;  // "-" placeholder keeps the field non-empty
            for (size_t i = 0; i < defs.size(); i++)
                defstr += (i ? " " : "") + defs[i];
            if (defstr.empty()) defstr = "-";
            out << name << '\t' << source << '\t' << entry << '\t' << defstr
                << '\t' << deps << '\n';
        };

        auto base_subsets = variant_subsets(feats);
        for (auto& entry : entries) {
            for (auto& [suffix, defs] : base_subsets) emit(entry, suffix, defs);
            if (noint64_cand) {
                for (auto& [suffix, defs] : base_subsets) {
                    std::vector<std::string> d = defs;
                    d.push_back(NOINT64.define);
                    emit(entry, suffix + NOINT64.suffix, d);
                }
            }
        }
    }
    std::fputs(out.str().c_str(), stdout);
    return 0;
}

// ---- embed mode ----------------------------------------------------------
std::string strip_feature_suffixes(std::string name) {
    const char* sfx[] = {NOINT64.suffix, INT8.suffix, ATOMIC.suffix};
    bool changed = true;
    while (changed) {
        changed = false;
        for (const char* s : sfx) {
            size_t n = std::strlen(s);
            if (name.size() > n && name.compare(name.size() - n, n, s) == 0) {
                name.resize(name.size() - n);
                changed = true;
            }
        }
    }
    return name;
}

int run_embed(const std::vector<std::string>& args) {
    std::string out_cpp, listfile;
    for (size_t i = 0; i < args.size(); i++) {
        if (args[i] == "--list" && i + 1 < args.size()) listfile = args[++i];
        else if (out_cpp.empty()) out_cpp = args[i];
    }
    if (out_cpp.empty() || listfile.empty()) {
        std::fprintf(stderr, "embed: usage: embed <out.cpp> --list <file>\n");
        return 2;
    }

    std::vector<std::string> paths;
    {
        std::ifstream lf(listfile);
        std::string line;
        while (std::getline(lf, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) paths.push_back(line);
        }
    }

    // name -> bytes, name -> capabilities
    std::map<std::string, std::string> blob_bytes;
    std::map<std::string, std::set<uint32_t>> blob_caps;
    for (auto& p : paths) {
        std::string name = basename_of(p);
        if (name.size() > 4 && name.compare(name.size() - 4, 4, ".spv") == 0)
            name.resize(name.size() - 4);
        bool ok = false;
        std::string data = read_file(p, &ok);
        if (!ok) {
            std::fprintf(stderr, "embed: missing blob %s\n", p.c_str());
            return 1;
        }
        blob_bytes[name] = data;
        blob_caps[name] = spirv_capabilities(data);
    }

    // Keep every base/feature blob; keep a .noint64 variant only when its base
    // actually declares Int64 (build_spirv.py's phase-2 condition).
    std::vector<std::string> kept;
    for (auto& [name, bytes] : blob_bytes) {
        (void)bytes;
        size_t n = std::strlen(NOINT64.suffix);
        bool is_noint64 = name.size() > n &&
            name.compare(name.size() - n, n, NOINT64.suffix) == 0;
        if (is_noint64) {
            std::string base = strip_feature_suffixes(name);
            auto it = blob_caps.find(base);
            if (it == blob_caps.end() || !it->second.count(CAP_INT64))
                continue;  // base has no Int64 -> no dedicated variant needed
        }
        kept.push_back(name);
    }
    // Sort by the ".spv" filename (not the bare name) so the emitted order is
    // stable and independent of which variant suffixes exist.
    std::sort(kept.begin(), kept.end(),
              [](const std::string& a, const std::string& b) {
                  return a + ".spv" < b + ".spv";
              });

    // Capability audit (mirrors build_spirv.py's final pass).
    bool failed = false;
    for (auto& name : kept) {
        const auto& caps = blob_caps[name];
        size_t n = std::strlen(NOINT64.suffix);
        bool is_noint64 = name.size() > n &&
            name.compare(name.size() - n, n, NOINT64.suffix) == 0;
        if (is_noint64 && caps.count(CAP_INT64)) {
            std::fprintf(stderr, "  CAPABILITY LEAK: %s still declares Int64 "
                "(unguarded 64-bit integer use; route it through "
                "vulkan/int64_compat.slang)\n", name.c_str());
            failed = true;
        }
        if (name.find(INT8.suffix) == std::string::npos &&
            (caps.count(CAP_INT8) || caps.count(CAP_8BIT_STORAGE))) {
            std::fprintf(stderr, "  CAPABILITY LEAK: %s declares 8-bit "
                "capabilities outside the .int8 variant (route byte access "
                "through vulkan/int8_compat.slang)\n", name.c_str());
            failed = true;
        }
    }
    if (failed) return 1;

    // Emit the translation unit (layout matches the former embed_spirv.py).
    std::ostringstream o;
    o << "// GENERATED by slang/spirv_tool.cpp -- DO NOT EDIT.\n"
         "#include <cstddef>\n#include <cstdint>\n\n"
         "namespace backend {\nnamespace vk {\n\n"
         "struct SpirvBlob { const char* name; const uint32_t* data;"
         " size_t size_bytes; };\n\n";
    for (size_t i = 0; i < kept.size(); i++) {
        const std::string& data = blob_bytes[kept[i]];
        o << "static const uint32_t _blob" << i << "[] = {\n";
        size_t words = data.size() / 4;
        char buf[16];
        for (size_t wi = 0; wi < words; wi += 8) {
            o << "   ";
            for (size_t j = wi; j < wi + 8 && j < words; j++) {
                uint32_t v = (uint32_t)(uint8_t)data[j * 4] |
                             ((uint32_t)(uint8_t)data[j * 4 + 1] << 8) |
                             ((uint32_t)(uint8_t)data[j * 4 + 2] << 16) |
                             ((uint32_t)(uint8_t)data[j * 4 + 3] << 24);
                std::snprintf(buf, sizeof(buf), " 0x%08xu,", v);
                o << buf;
            }
            o << "\n";
        }
        o << "};\n";
    }
    o << "\nextern const SpirvBlob g_spirv_blobs[];\n"
         "extern const size_t g_spirv_blob_count;\n"
         "const SpirvBlob g_spirv_blobs[] = {\n";
    for (size_t i = 0; i < kept.size(); i++)
        o << "    {\"" << kept[i] << "\", _blob" << i << ", "
          << blob_bytes[kept[i]].size() << "},\n";
    if (kept.empty())
        o << "    {nullptr, nullptr, 0},  // keep the array non-empty\n";
    o << "};\nconst size_t g_spirv_blob_count = " << kept.size() << ";\n\n"
         "}  // namespace vk\n}  // namespace backend\n";

    // Always (re)write: the embed edge only runs when a .spv changed, so the
    // fresh mtime is what keeps Ninja's output-newer-than-inputs check happy.
    std::string text = o.str();
    std::ofstream f(out_cpp, std::ios::binary | std::ios::trunc);
    if (!f) {
        std::fprintf(stderr, "embed: cannot write %s\n", out_cpp.c_str());
        return 1;
    }
    f << text;
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s discover|embed ...\n", argv[0]);
        return 2;
    }
    std::vector<std::string> args(argv + 2, argv + argc);
    std::string mode = argv[1];
    if (mode == "discover") return run_discover(args);
    if (mode == "embed") return run_embed(args);
    std::fprintf(stderr, "unknown mode: %s\n", mode.c_str());
    return 2;
}
