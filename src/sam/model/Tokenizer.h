#pragma once
// CLIP byte-pair-encoding tokenizer, loaded from the vocabulary embedded at the
// end of the checkpoint.
//
// A faithful port of the reference tokenizer, including its two quirks that
// matter for output parity:
//
//   * the byte->unicode table (printable bytes map to themselves, the rest to
//     U+0100..U+0143), so any input byte is representable in the vocab;
//   * a hand-rolled approximation of CLIP's pre-tokenization regex. The real
//     pattern needs Unicode character classes; here any multi-byte UTF-8 start
//     byte counts as a letter, contractions are matched before letters, and
//     digits are split one at a time. For the noun phrases SAM 3 takes as
//     prompts this is exact.

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace sam {

class BpeTokenizer {
public:
    // Reads the "tok\0" block: vocabulary then merge list.
    bool loadFromStream(std::istream& in);
    bool loaded() const { return !encoder_.empty(); }

    // Text -> exactly `ctx_len` ids: [SOT, ..., EOT, 0, 0, ...]. Truncates to
    // ctx_len keeping EOT last, matching the reference.
    std::vector<int32_t> encode(const std::string& text, int ctx_len) const;

    size_t vocabSize() const { return encoder_.size(); }
    size_t mergeCount() const { return merge_ranks_.size(); }

    int sotToken() const { return sot_; }
    int eotToken() const { return eot_; }

private:
    std::string bpe(const std::string& token) const;

    std::unordered_map<std::string, int> encoder_;
    std::unordered_map<std::string, int> merge_ranks_;   // "a\x1fb" -> rank
    std::unordered_map<uint8_t, std::string> byte_encoder_;
    mutable std::unordered_map<std::string, std::string> cache_;
    int sot_ = 49406;
    int eot_ = 49407;
};

}  // namespace sam
