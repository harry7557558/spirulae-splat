#pragma once

// Fast decimal number parsing for the text-format dataset/model readers.
//
// musl's strtod/strtof (the libc Emscripten links) performs arbitrary-
// precision correct rounding (floatscan.c) and costs ~1.5us per call in WASM;
// parsing a 60+ MB ascii PLY or a 100+ MB OBJ makes tens of millions of such
// calls, turning sub-second loads into tens of seconds. Every consumer here
// stores float32, so correct double rounding is wasted effort.
//
// fast_strtod parses the common decimal forms ([+-]ddd[.ddd][eNN]) with an
// int64 mantissa and a pow-10 scale (error <= ~2 ulp of double — invisible in
// float32), and falls back to strtod for anything exotic (hex, inf/nan,
// out-of-range exponents). Interface matches strtod so it's a drop-in.

#include <cstdint>
#include <cstdlib>

namespace fastfloat {

inline double pow10i(int e) {
    static double pos[309], neg[309];
    static bool init = false;
    if (!init) {
        pos[0] = neg[0] = 1.0;
        for (int i = 1; i < 309; i++) {
            pos[i] = pos[i-1] * 10.0;
            neg[i] = neg[i-1] * 0.1;
        }
        init = true;
    }
    return e >= 0 ? pos[e] : neg[-e];
}

}  // namespace fastfloat

inline double fast_strtod(const char* p, char** endptr) {
    const char* s = p;
    while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r') s++;
    bool neg = false;
    if (*s == '-' || *s == '+') { neg = (*s == '-'); s++; }
    if (s[0] == '0' && (s[1] == 'x' || s[1] == 'X'))
        return std::strtod(p, endptr);     // hex float

    uint64_t mant = 0;
    int ndig = 0, frac = 0;
    bool any = false;
    while (*s >= '0' && *s <= '9') {
        if (ndig < 19) { mant = mant * 10 + (uint64_t)(*s - '0'); ndig++; }
        else frac--;                       // dropped integer digit = *10 later
        any = true; s++;
    }
    if (*s == '.') {
        s++;
        while (*s >= '0' && *s <= '9') {
            if (ndig < 19) { mant = mant * 10 + (uint64_t)(*s - '0'); ndig++; frac++; }
            any = true; s++;               // excess fraction digits: truncate
        }
    }
    if (!any) return std::strtod(p, endptr);   // inf / nan / hex / no digits

    int exp10 = 0;
    if (*s == 'e' || *s == 'E') {
        const char* es = s + 1;
        bool eneg = false;
        if (*es == '-' || *es == '+') { eneg = (*es == '-'); es++; }
        if (*es >= '0' && *es <= '9') {
            int ev = 0;
            while (*es >= '0' && *es <= '9' && ev < 100000) { ev = ev * 10 + (*es - '0'); es++; }
            exp10 = eneg ? -ev : ev;
            s = es;
        }
    }
    int e = exp10 - frac;
    if (e < -308 || e > 308) return std::strtod(p, endptr);

    if (endptr) *endptr = (char*)s;
    double v = (double)mant * fastfloat::pow10i(e);
    return neg ? -v : v;
}

inline float fast_strtof(const char* p, char** endptr) {
    return (float)fast_strtod(p, endptr);
}

// Fast [+-]digits integer parse (strtol without base/locale overhead).
// 64-bit accumulator: wasm32 long is 32 bits and e.g. COLMAP point3D ids can
// exceed it.
inline long long fast_strtol(const char* p, char** endptr) {
    const char* s = p;
    while (*s == ' ' || *s == '\t' || *s == '\n' || *s == '\r') s++;
    bool neg = false;
    if (*s == '-' || *s == '+') { neg = (*s == '-'); s++; }
    long long v = 0;
    bool any = false;
    while (*s >= '0' && *s <= '9') { v = v * 10 + (*s - '0'); any = true; s++; }
    if (endptr) *endptr = (char*)(any ? s : p);
    return neg ? -v : v;
}
