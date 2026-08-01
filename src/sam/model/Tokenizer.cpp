#include "sam/model/Tokenizer.h"

#include "nn/core/Error.h"
#include "nn/core/Log.h"

#include <algorithm>
#include <climits>
#include <cctype>
#include <istream>

namespace sam {

namespace {

int utf8_len(uint8_t c) {
    if (c < 0x80) return 1;
    if (c < 0xC0) return 1;  // stray continuation byte; consume it alone
    if (c < 0xE0) return 2;
    if (c < 0xF0) return 3;
    return 4;
}

std::string codepoint_to_utf8(int cp) {
    std::string s;
    if (cp < 0x80) {
        s += (char)cp;
    } else if (cp < 0x800) {
        s += (char)(0xC0 | (cp >> 6));
        s += (char)(0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
        s += (char)(0xE0 | (cp >> 12));
        s += (char)(0x80 | ((cp >> 6) & 0x3F));
        s += (char)(0x80 | (cp & 0x3F));
    } else {
        s += (char)(0xF0 | (cp >> 18));
        s += (char)(0x80 | ((cp >> 12) & 0x3F));
        s += (char)(0x80 | ((cp >> 6) & 0x3F));
        s += (char)(0x80 | (cp & 0x3F));
    }
    return s;
}

bool is_letter(const std::string& s, size_t i) {
    uint8_t c = (uint8_t)s[i];
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')) return true;
    return c >= 0xC0;  // any multi-byte start byte counts as a letter
}

// GPT-2 / CLIP byte->unicode table: printable ASCII and Latin-1 map to
// themselves, the remaining 68 bytes to U+0100 upward, so the vocabulary covers
// every possible byte without an <unk>.
std::unordered_map<uint8_t, std::string> build_byte_encoder() {
    std::vector<int> bs;
    for (int i = 33; i <= 126; ++i) bs.push_back(i);
    for (int i = 161; i <= 172; ++i) bs.push_back(i);
    for (int i = 174; i <= 255; ++i) bs.push_back(i);
    std::vector<int> cs(bs.begin(), bs.end());
    int n = 0;
    for (int b = 0; b < 256; ++b)
        if (std::find(bs.begin(), bs.end(), b) == bs.end()) {
            bs.push_back(b);
            cs.push_back(256 + n++);
        }
    std::unordered_map<uint8_t, std::string> enc;
    for (size_t i = 0; i < bs.size(); ++i)
        enc[(uint8_t)bs[i]] = codepoint_to_utf8(cs[i]);
    return enc;
}

// 0x1F cannot appear inside a byte-encoded BPE token, so it is a safe joiner.
std::string merge_key(const std::string& a, const std::string& b) {
    std::string k;
    k.reserve(a.size() + 1 + b.size());
    k += a;
    k += '\x1f';
    k += b;
    return k;
}

std::vector<std::string> utf8_chars(const std::string& s) {
    std::vector<std::string> out;
    size_t i = 0;
    while (i < s.size()) {
        int n = utf8_len((uint8_t)s[i]);
        out.push_back(s.substr(i, (size_t)n));
        i += (size_t)n;
    }
    return out;
}

std::vector<std::string> pretokenize(const std::string& text) {
    std::vector<std::string> tokens;
    size_t i = 0;
    const size_t n = text.size();
    while (i < n) {
        uint8_t c = (uint8_t)text[i];
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') { ++i; continue; }

        if (i + 15 <= n && text.compare(i, 15, "<|startoftext|>") == 0) {
            tokens.push_back("<|startoftext|>");
            i += 15;
            continue;
        }
        if (i + 13 <= n && text.compare(i, 13, "<|endoftext|>") == 0) {
            tokens.push_back("<|endoftext|>");
            i += 13;
            continue;
        }
        // Contractions come before the letter rule: an apostrophe is not a
        // letter, so 's/'t/'re/... would otherwise split wrongly.
        if (c == '\'') {
            if (i + 2 <= n) {
                char c2 = text[i + 1];
                if (c2 == 's' || c2 == 't' || c2 == 'm' || c2 == 'd') {
                    tokens.push_back(text.substr(i, 2));
                    i += 2;
                    continue;
                }
            }
            if (i + 3 <= n) {
                std::string c3 = text.substr(i + 1, 2);
                if (c3 == "re" || c3 == "ve" || c3 == "ll") {
                    tokens.push_back(text.substr(i, 3));
                    i += 3;
                    continue;
                }
            }
        }
        if (is_letter(text, i)) {
            size_t start = i;
            while (i < n && is_letter(text, i)) i += (size_t)utf8_len((uint8_t)text[i]);
            tokens.push_back(text.substr(start, i - start));
            continue;
        }
        if (c >= '0' && c <= '9') {  // one digit at a time, as CLIP's [\p{N}]
            tokens.push_back(text.substr(i, 1));
            ++i;
            continue;
        }
        size_t start = i;
        while (i < n) {
            uint8_t ch = (uint8_t)text[i];
            if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r') break;
            if (is_letter(text, i)) break;
            if (ch >= '0' && ch <= '9') break;
            ++i;
        }
        if (i > start) tokens.push_back(text.substr(start, i - start));
    }
    return tokens;
}

template <typename T>
bool read_pod(std::istream& in, T& v) {
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    return (bool)in;
}

}  // namespace

bool BpeTokenizer::loadFromStream(std::istream& in) {
    uint32_t magic = 0;
    if (!read_pod(in, magic)) return false;
    if (magic != 0x746F6B00u) {
        NN_LOG_ERROR("tokenizer: bad magic 0x%08x (expected 0x746f6b00 \"tok\")\n",
                       magic);
        return false;
    }

    int32_t n_vocab = 0;
    if (!read_pod(in, n_vocab) || n_vocab <= 0 || n_vocab > 1000000) return false;
    encoder_.reserve((size_t)n_vocab);
    for (int32_t i = 0; i < n_vocab; ++i) {
        int32_t len = 0;
        if (!read_pod(in, len) || len < 0 || len > 4096) return false;
        std::string token((size_t)len, '\0');
        in.read(&token[0], len);
        int32_t id = 0;
        if (!read_pod(in, id)) return false;
        encoder_.emplace(std::move(token), id);
    }

    int32_t n_merges = 0;
    if (!read_pod(in, n_merges) || n_merges < 0 || n_merges > 1000000) return false;
    merge_ranks_.reserve((size_t)n_merges);
    for (int32_t i = 0; i < n_merges; ++i) {
        int32_t la = 0, lb = 0;
        if (!read_pod(in, la) || la < 0 || la > 4096) return false;
        std::string a((size_t)la, '\0');
        in.read(&a[0], la);
        if (!read_pod(in, lb) || lb < 0 || lb > 4096) return false;
        std::string b((size_t)lb, '\0');
        in.read(&b[0], lb);
        merge_ranks_.emplace(merge_key(a, b), i);
    }
    if (!in) return false;

    byte_encoder_ = build_byte_encoder();
    NN_LOG_INFO("[model] tokenizer: %zu vocab entries, %zu merges\n", encoder_.size(),
                  merge_ranks_.size());
    return true;
}

std::string BpeTokenizer::bpe(const std::string& token) const {
    auto it = cache_.find(token);
    if (it != cache_.end()) return it->second;

    std::vector<std::string> word = utf8_chars(token);
    if (word.empty()) return "";
    word.back() += "</w>";
    if (word.size() == 1) {
        cache_[token] = word[0];
        return word[0];
    }

    while (true) {
        int best_rank = INT_MAX;
        std::string best_a, best_b;
        for (size_t i = 0; i + 1 < word.size(); ++i) {
            auto mit = merge_ranks_.find(merge_key(word[i], word[i + 1]));
            if (mit != merge_ranks_.end() && mit->second < best_rank) {
                best_rank = mit->second;
                best_a = word[i];
                best_b = word[i + 1];
            }
        }
        if (best_rank == INT_MAX) break;

        const std::string merged = best_a + best_b;
        std::vector<std::string> next;
        next.reserve(word.size());
        for (size_t i = 0; i < word.size();) {
            if (i + 1 < word.size() && word[i] == best_a && word[i + 1] == best_b) {
                next.push_back(merged);
                i += 2;
            } else {
                next.push_back(word[i]);
                ++i;
            }
        }
        word.swap(next);
        if (word.size() == 1) break;
    }

    std::string result;
    for (size_t i = 0; i < word.size(); ++i) {
        if (i) result += ' ';
        result += word[i];
    }
    cache_[token] = result;
    return result;
}

std::vector<int32_t> BpeTokenizer::encode(const std::string& text, int ctx_len) const {
    // Lowercase, then collapse runs of whitespace to a single space.
    std::string clean;
    clean.reserve(text.size());
    bool last_ws = true;
    for (char raw : text) {
        char c = (raw >= 'A' && raw <= 'Z') ? (char)(raw + 32) : raw;
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            if (!last_ws) { clean += ' '; last_ws = true; }
        } else {
            clean += c;
            last_ws = false;
        }
    }
    if (!clean.empty() && clean.back() == ' ') clean.pop_back();

    std::vector<int32_t> ids;
    ids.push_back(sot_);
    for (const auto& word : pretokenize(clean)) {
        std::string encoded;
        for (unsigned char b : word) {
            auto it = byte_encoder_.find(b);
            if (it != byte_encoder_.end()) encoded += it->second;
        }
        const std::string pieces = bpe(encoded);
        size_t start = 0;
        while (start < pieces.size()) {
            size_t end = pieces.find(' ', start);
            if (end == std::string::npos) end = pieces.size();
            auto eit = encoder_.find(pieces.substr(start, end - start));
            // Unknown pieces are dropped, matching CLIP: with the byte encoder
            // above every byte sequence is in the vocabulary, so this only
            // fires on a truncated vocabulary file.
            if (eit != encoder_.end()) ids.push_back(eit->second);
            start = end + 1;
        }
    }
    ids.push_back(eot_);

    if ((int)ids.size() > ctx_len) {
        ids.resize((size_t)ctx_len);
        ids.back() = eot_;
    }
    ids.resize((size_t)ctx_len, 0);
    return ids;
}

}  // namespace sam
