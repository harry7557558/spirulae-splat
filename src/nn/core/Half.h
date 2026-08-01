#pragma once
// IEEE-754 binary16 <-> binary32 on the host.
//
// Needed by the weight loader (checkpoints are fp16, and quantized ones
// dequantize to fp16) and by the debug readback path. Round-to-nearest-even,
// with subnormals and infinities handled, so a round-trip through fp16 matches
// what the shader's f16tof32 sees bit for bit.

#include <cstdint>
#include <cstring>

namespace nn {

inline float half_to_float(uint16_t h) {
    const uint32_t sign = (uint32_t)(h & 0x8000u) << 16;
    const uint32_t exp = (h >> 10) & 0x1Fu;
    const uint32_t mant = h & 0x3FFu;
    uint32_t bits;
    if (exp == 0) {
        if (mant == 0) {
            bits = sign;  // +/- zero
        } else {
            // Subnormal: normalize by shifting the mantissa up.
            uint32_t e = 0;
            uint32_t m = mant;
            while ((m & 0x400u) == 0) { m <<= 1; ++e; }
            m &= 0x3FFu;
            bits = sign | ((127 - 15 - e + 1) << 23) | (m << 13);
        }
    } else if (exp == 0x1F) {
        bits = sign | 0x7F800000u | (mant << 13);  // inf / nan
    } else {
        bits = sign | ((exp + (127 - 15)) << 23) | (mant << 13);
    }
    float f;
    std::memcpy(&f, &bits, 4);
    return f;
}

inline uint16_t float_to_half(float f) {
    uint32_t bits;
    std::memcpy(&bits, &f, 4);
    const uint16_t sign = (uint16_t)((bits >> 16) & 0x8000u);
    int32_t exp = (int32_t)((bits >> 23) & 0xFFu) - 127 + 15;
    uint32_t mant = bits & 0x7FFFFFu;

    if (((bits >> 23) & 0xFFu) == 0xFFu)                        // inf / nan
        return (uint16_t)(sign | 0x7C00u | (mant ? 0x200u : 0u));
    if (exp >= 0x1F) return (uint16_t)(sign | 0x7C00u);          // overflow -> inf
    if (exp <= 0) {
        if (exp < -10) return sign;                              // underflow -> 0
        // Subnormal: shift the implicit 1 in and round to nearest even.
        mant |= 0x800000u;
        uint32_t shift = (uint32_t)(14 - exp);
        uint32_t m = mant >> shift;
        uint32_t rem = mant & ((1u << shift) - 1u);
        uint32_t halfway = 1u << (shift - 1);
        if (rem > halfway || (rem == halfway && (m & 1u))) ++m;
        return (uint16_t)(sign | m);
    }
    uint32_t m = mant >> 13;
    uint32_t rem = mant & 0x1FFFu;
    if (rem > 0x1000u || (rem == 0x1000u && (m & 1u))) {
        ++m;
        if (m == 0x400u) { m = 0; ++exp; if (exp >= 0x1F) return (uint16_t)(sign | 0x7C00u); }
    }
    return (uint16_t)(sign | ((uint32_t)exp << 10) | m);
}

}  // namespace nn
