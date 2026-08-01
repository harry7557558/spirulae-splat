#pragma once
// Bit-level readers for video bitstreams.
//
// `BitReader` is the plain MSB-first reader AV1 uses. `RbspReader` adds the
// H.264/H.265 emulation-prevention rule: a 0x03 byte following two zero bytes
// is not payload. Both clamp rather than fault -- a truncated or corrupt NAL
// unit should surface as a parse failure, never as a read past the buffer.

#include <cstdint>
#include <cstring>
#include <vector>

namespace video {

class BitReader {
public:
    BitReader() = default;
    BitReader(const uint8_t* data, size_t size) : data_(data), size_(size) {}

    // Reads `n` (<= 32) bits MSB-first. Past the end reads zeros and latches
    // `overrun`.
    uint32_t u(int n) {
        uint32_t v = 0;
        for (int i = 0; i < n; ++i) v = (v << 1) | bit();
        return v;
    }
    uint32_t bit() {
        const size_t byte = pos_ >> 3;
        if (byte >= size_) {
            overrun_ = true;
            return 0;
        }
        const uint32_t b = (data_[byte] >> (7 - (pos_ & 7))) & 1u;
        ++pos_;
        return b;
    }
    bool flag() { return bit() != 0; }

    // Exp-Golomb, unsigned and signed (H.264/H.265 ue(v) / se(v)).
    uint32_t ue() {
        int zeros = 0;
        while (!overrun_ && bit() == 0 && zeros < 32) ++zeros;
        if (zeros >= 32) {
            overrun_ = true;
            return 0;
        }
        return (1u << zeros) - 1u + (zeros ? u(zeros) : 0u);
    }
    int32_t se() {
        const uint32_t k = ue();
        const int32_t v = (int32_t)((k + 1) >> 1);
        return (k & 1) ? v : -v;
    }

    // AV1 leb128 / uvlc / le(n) / ns(n) / su(n).
    uint64_t leb128(int* bytes_read = nullptr) {
        uint64_t value = 0;
        int i = 0;
        for (; i < 8; ++i) {
            const uint32_t b = u(8);
            value |= (uint64_t)(b & 0x7f) << (i * 7);
            if (!(b & 0x80)) { ++i; break; }
        }
        if (bytes_read) *bytes_read = i;
        return value;
    }
    uint32_t uvlc() {
        int zeros = 0;
        while (!overrun_ && bit() == 0 && zeros < 32) ++zeros;
        if (zeros >= 32) return UINT32_MAX;
        return u(zeros) + (1u << zeros) - 1u;
    }
    uint32_t le(int n) {  // little-endian bytes, byte-aligned
        uint32_t t = 0;
        for (int i = 0; i < n; ++i) t |= u(8) << (i * 8);
        return t;
    }
    uint32_t ns(uint32_t n) {  // AV1 non-symmetric
        if (n <= 1) return 0;
        uint32_t w = 0, x = n;
        while (x) { ++w; x >>= 1; }
        const uint32_t m = (1u << w) - n;
        const uint32_t v = u((int)w - 1);
        if (v < m) return v;
        return (v << 1) - m + bit();
    }
    int32_t su(int n) {
        int32_t value = (int32_t)u(n);
        const int32_t sign_mask = 1 << (n - 1);
        if (value & sign_mask) value -= 2 * sign_mask;
        return value;
    }

    size_t bitPos() const { return pos_; }
    void   seekBits(size_t p) { pos_ = p; }
    void   byteAlign() { pos_ = (pos_ + 7) & ~(size_t)7; }
    size_t sizeBits() const { return size_ * 8; }
    bool   overrun() const { return overrun_; }
    bool   moreRbspData() const {
        // Non-zero bits remain before the rbsp_stop_one_bit's trailing zeros.
        if (pos_ >= size_ * 8) return false;
        size_t last = size_ * 8;
        while (last > pos_) {
            const size_t i = last - 1;
            if ((data_[i >> 3] >> (7 - (i & 7))) & 1u) break;
            --last;
        }
        return last > pos_ + 1;  // strictly before the stop bit
    }

protected:
    const uint8_t* data_ = nullptr;
    size_t size_ = 0;
    size_t pos_ = 0;
    bool   overrun_ = false;
};

// Strips emulation-prevention bytes once, then reads from the copy. Doing it up
// front rather than per-bit keeps the reader itself branch-free and lets the
// H.265 slice header record exact bit offsets into the RBSP, which
// StdVideoDecodeH265PictureInfo::NumBitsForSTRefPicSetInSlice requires.
class RbspReader : public BitReader {
public:
    RbspReader() = default;
    RbspReader(const uint8_t* nal, size_t size, int header_bytes = 0) {
        set(nal, size, header_bytes);
    }

    void set(const uint8_t* nal, size_t size, int header_bytes) {
        rbsp_.clear();
        rbsp_.reserve(size);
        int zeros = 0;
        for (size_t i = (size_t)header_bytes; i < size; ++i) {
            const uint8_t b = nal[i];
            if (zeros >= 2 && b == 0x03) {
                zeros = 0;
                continue;  // emulation prevention byte
            }
            zeros = (b == 0) ? zeros + 1 : 0;
            rbsp_.push_back(b);
        }
        data_ = rbsp_.data();
        size_ = rbsp_.size();
        pos_ = 0;
        overrun_ = false;
    }

    const std::vector<uint8_t>& rbsp() const { return rbsp_; }

private:
    std::vector<uint8_t> rbsp_;
};

}  // namespace video
