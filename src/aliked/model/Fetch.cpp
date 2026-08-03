#include "aliked/model/Fetch.h"

#include "nn/core/Error.h"
#include "nn/core/Log.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <vector>

namespace fs = std::filesystem;

namespace aliked {
namespace {

// The same URL / filename / SHA-256 triples COLMAP carries in
// src/colmap/feature/resources.h. Keep them identical: the point of fetching
// from COLMAP's release is that a parity run compares implementations, not
// checkpoints.
const ModelSource kSources[] = {
    {"aliked-n16rot", "aliked-n16rot.onnx",
     "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n16rot.onnx",
     "39c423d0a6f03d39ec89d3d1d61853765c2fb6a8b8381376c703e5758778a547", 2997054ull},
    {"aliked-n32", "aliked-n32.onnx",
     "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-n32.onnx",
     "a077728a02d2de1a775c66df6de8cfeb7c6b51ca57572c64c680131c988c8b3c", 4205634ull},
    {"aliked-lightglue", "aliked-lightglue.onnx",
     "https://github.com/colmap/colmap/releases/download/3.13.0/aliked-lightglue.onnx",
     "b9a5de7204648b18a8cf5dcac819f9d30de1a5961ef03756803c8b86c2dceb8d", 0ull},
};

// ---------------------------------------------------------------------------
// SHA-256 (FIPS 180-4)
//
// Here because the repository has no crypto dependency and should not grow one
// for this. Verifying matters more than usual: these bytes come off the
// network and go straight into a parser.
// ---------------------------------------------------------------------------

struct Sha256 {
    uint32_t h[8] = {0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
                     0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u};
    uint8_t  buf[64] = {};
    size_t   buf_len = 0;
    uint64_t total = 0;

    static uint32_t ror(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

    void block(const uint8_t* p) {
        static const uint32_t K[64] = {
            0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
            0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
            0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
            0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
            0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
            0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
            0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
            0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
            0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
            0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
            0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = ((uint32_t)p[i * 4] << 24) | ((uint32_t)p[i * 4 + 1] << 16) |
                   ((uint32_t)p[i * 4 + 2] << 8) | (uint32_t)p[i * 4 + 3];
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = ror(w[i - 15], 7) ^ ror(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const uint32_t s1 = ror(w[i - 2], 17) ^ ror(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h[0], b = h[1], c = h[2], d = h[3];
        uint32_t e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t S1 = ror(e, 6) ^ ror(e, 11) ^ ror(e, 25);
            const uint32_t ch = (e & f) ^ (~e & g);
            const uint32_t t1 = hh + S1 + ch + K[i] + w[i];
            const uint32_t S0 = ror(a, 2) ^ ror(a, 13) ^ ror(a, 22);
            const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t t2 = S0 + maj;
            hh = g; g = f; f = e; e = d + t1;
            d = c; c = b; b = a; a = t1 + t2;
        }
        h[0] += a; h[1] += b; h[2] += c; h[3] += d;
        h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
    }

    void update(const uint8_t* p, size_t n) {
        total += n;
        while (n) {
            const size_t take = std::min(n, sizeof(buf) - buf_len);
            std::memcpy(buf + buf_len, p, take);
            buf_len += take;
            p += take;
            n -= take;
            if (buf_len == sizeof(buf)) { block(buf); buf_len = 0; }
        }
    }

    std::string hex() {
        const uint64_t bits = total * 8;
        uint8_t pad = 0x80;
        update(&pad, 1);
        pad = 0x00;
        while (buf_len != 56) update(&pad, 1);
        uint8_t len[8];
        for (int i = 0; i < 8; ++i) len[i] = (uint8_t)(bits >> (56 - 8 * i));
        // update() would count these into `total`, but the length is already
        // frozen in `bits`, so the extra count is harmless.
        update(len, 8);

        static const char* kHex = "0123456789abcdef";
        std::string out(64, '0');
        for (int i = 0; i < 8; ++i)
            for (int b = 0; b < 4; ++b) {
                const uint8_t v = (uint8_t)(h[i] >> (24 - 8 * b));
                out[(size_t)(i * 8 + b * 2)] = kHex[v >> 4];
                out[(size_t)(i * 8 + b * 2 + 1)] = kHex[v & 0xF];
            }
        return out;
    }
};

std::string env_str(const char* name) {
    const char* v = std::getenv(name);
    return v ? std::string(v) : std::string();
}

fs::path cache_root() {
#ifdef _WIN32
    fs::path dir = env_str("LOCALAPPDATA");
#else
    fs::path dir = env_str("XDG_CACHE_HOME");
    if (dir.empty()) {
        const std::string home = env_str("HOME");
        if (!home.empty()) dir = fs::path(home) / ".cache";
    }
#endif
    if (dir.empty()) dir = ".";
    return dir / "spirulae-splat";
}

bool have_curl() {
#ifdef _WIN32
    return std::system("curl --version >NUL 2>&1") == 0;
#else
    return std::system("curl --version >/dev/null 2>&1") == 0;
#endif
}

}  // namespace

const ModelSource* find_model_source(const std::string& id) {
    for (const ModelSource& s : kSources)
        if (id == s.id) return &s;
    return nullptr;
}

std::string model_cache_path(const ModelSource& src) {
    return (cache_root() / "models" / src.file).string();
}

std::string sha256_file(const std::string& path) {
    std::ifstream fin(path, std::ios::binary);
    if (!fin) return {};
    Sha256 sha;
    std::vector<char> chunk(1 << 20);
    while (fin) {
        fin.read(chunk.data(), (std::streamsize)chunk.size());
        const std::streamsize got = fin.gcount();
        if (got > 0) sha.update(reinterpret_cast<const uint8_t*>(chunk.data()), (size_t)got);
    }
    return sha.hex();
}

std::string ensure_model(const ModelSource& src) {
    const fs::path dst = model_cache_path(src);
    std::error_code ec;

    if (fs::exists(dst, ec)) {
        const std::string got = sha256_file(dst.string());
        if (got == src.sha256) return dst.string();
        // A cached file that does not hash is a failed or interrupted download
        // from a previous run, not a reason to stop: say so and refetch.
        NN_LOG_WARN("[aliked] cached %s does not match its checksum; re-downloading\n",
                    src.file);
        fs::remove(dst, ec);
    }

    fs::create_directories(dst.parent_path(), ec);
    NN_CHECK(!ec, "cannot create %s: %s", dst.parent_path().string().c_str(),
             ec.message().c_str());

    NN_CHECK(have_curl(),
             "curl was not found, and it is how checkpoints are fetched.\n"
             "  Install curl, or download\n    %s\n  to\n    %s\n  by hand.",
             src.url, dst.string().c_str());

    fs::path part = dst;
    part += ".part";

    NN_LOG_INFO("[aliked] fetching %s (%.1f MB) from %s\n", src.file,
                (double)src.bytes / 1e6, src.url);
    // -C - resumes a partial .part file; -f makes an HTTP error an exit code
    // rather than a saved error page. Downloading into .part and renaming is
    // what keeps an interrupted fetch from ever looking like a complete model.
    std::string cmd = "curl -L -f --progress-bar -C - -o \"" + part.string() + "\" \"" +
                      std::string(src.url) + "\"";
    const int rc = std::system(cmd.c_str());
    if (rc != 0) {
        fs::remove(part, ec);
        nn::fail("downloading %s failed (curl exit %d).\n"
                 "  Fetch it by hand from\n    %s\n  and save it as\n    %s",
                 src.file, rc, src.url, dst.string().c_str());
    }

    const std::string got = sha256_file(part.string());
    if (got != src.sha256) {
        fs::remove(part, ec);
        nn::fail("%s downloaded but its SHA-256 is\n    %s\n  expected\n    %s\n"
                 "  The file was discarded.",
                 src.file, got.c_str(), src.sha256);
    }

    fs::rename(part, dst, ec);
    NN_CHECK(!ec, "cannot move the download into place: %s", ec.message().c_str());
    NN_LOG_INFO("[aliked] saved %s\n", dst.string().c_str());
    return dst.string();
}

std::string resolve_model(const std::string& id_or_path) {
    if (const ModelSource* src = find_model_source(id_or_path)) return ensure_model(*src);

    std::error_code ec;
    NN_CHECK(fs::exists(id_or_path, ec),
             "'%s' is neither a known model id (aliked-n16rot, aliked-n32, "
             "aliked-lightglue) nor a file that exists",
             id_or_path.c_str());
    return id_or_path;
}

}  // namespace aliked
