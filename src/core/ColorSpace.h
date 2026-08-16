#pragma once

// Gamut matrices and the sRGB transfer, shared by everything that has to move
// pixels between a capture colour space and sRGB: the trainer's seed colours,
// the SfM front end, and the AI models. Header-only so `src/sfm/` and `src/nn/`
// can use it without depending on the app layer.
//
// Matrices are row-major 3x3, source primaries -> Rec.709, chromatically
// adapted (white maps to white).

#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace colorspace {

using Mat3 = std::array<float, 9>;

inline constexpr const char* kGamuts[] = {
    "Rec.709", "ACES2065-1", "ACEScg", "Rec.2020", "AdobeRGB", "DCI-P3",
};

// "" and "none" are Rec.709, i.e. the identity.
inline Mat3 gamut_to_rec709(const std::string& name) {
    if (name.empty() || name == "none" || name == "Rec.709")
        return {1,0,0, 0,1,0, 0,0,1};
    if (name == "ACES2065-1") return {
        2.5247180476f, -1.1325619434f, -0.3921561044f,
        -0.2776344819f, 1.3709123773f, -0.0932778953f,
        -0.0165202369f, -0.1479259606f, 1.1644461975f};
    if (name == "ACEScg") return {
        1.7072552160f, -0.6200352595f, -0.0872199564f,
        -0.1311566587f, 1.1391010566f, -0.0079443978f,
        -0.0245499075f, -0.1248045805f, 1.1493544880f};
    if (name == "Rec.2020") return {
        1.6604910021f, -0.5876411388f, -0.0728498633f,
        -0.1245504745f, 1.1328998971f, -0.0083494226f,
        -0.0181507634f, -0.1005788980f, 1.1187296614f};
    if (name == "AdobeRGB") return {
        1.3983671735f, -0.3983451225f, 0.0000054016f,
        -0.0000103176f, 0.9999916496f, -0.0000039459f,
        -0.0000003709f, -0.0429269510f, 1.0429319656f};
    if (name == "DCI-P3") return {
        1.1548337042f, -0.1451763523f, -0.0096573518f,
        -0.0393300117f, 1.0378282998f, 0.0015017119f,
        -0.0184786235f, -0.0689101110f, 1.0873887345f};
    throw std::runtime_error("unsupported color gamut: " + name);
}

inline Mat3 invert3x3(const Mat3& m) {
    const float a = m[0], b = m[1], c = m[2],
                d = m[3], e = m[4], g = m[5],
                h = m[6], i = m[7], j = m[8];
    const float det = a*(e*j - g*i) - b*(d*j - g*h) + c*(d*i - e*h);
    if (std::abs(det) < 1e-20f) throw std::runtime_error("singular color matrix");
    const float s = 1.0f / det;
    return {(e*j - g*i)*s, (c*i - b*j)*s, (b*g - c*e)*s,
            (g*h - d*j)*s, (a*j - c*h)*s, (c*d - a*g)*s,
            (d*i - e*h)*s, (b*h - a*i)*s, (a*e - b*d)*s};
}

// 0.04045 is the branch point (= 12.92 * 0.0031308); 0.055 is the offset.
inline float srgb_to_linear(float x) {
    return x < 0.04045f ? x * (1.0f / 12.92f)
                        : std::pow((x + 0.055f) * (1.0f / 1.055f), 2.4f);
}

inline float linear_to_srgb(float x) {
    return x < 0.0031308f ? 12.92f * x
                          : 1.055f * std::pow(std::max(x, 0.0f), 1.0f / 2.4f) - 0.055f;
}

inline void apply3x3(const Mat3& m, float v[3]) {
    const float t0 = v[0], t1 = v[1], t2 = v[2];
    for (int r = 0; r < 3; r++)
        v[r] = m[r*3+0]*t0 + m[r*3+1]*t1 + m[r*3+2]*t2;
}

// True when (gamut, is_linear) is already plain sRGB, so callers can skip the
// per-pixel work entirely.
inline bool is_identity(const std::string& gamut, bool is_linear) {
    return !is_linear && (gamut.empty() || gamut == "none" || gamut == "Rec.709");
}

// Interleaved 8-bit RGB, in place. `n` is the pixel count.
//
// Values outside the destination gamut clamp -- these buffers feed feature
// extraction and 8-bit point colours, both of which are 0..255 by definition.
inline void to_srgb_inplace(uint8_t* rgb, size_t n,
                            const std::string& gamut, bool is_linear) {
    if (is_identity(gamut, is_linear)) return;
    const Mat3 m = gamut_to_rec709(gamut);
    for (size_t i = 0; i < n; i++) {
        float v[3];
        for (int c = 0; c < 3; c++) {
            const float x = rgb[3*i + c] * (1.0f / 255.0f);
            v[c] = is_linear ? x : srgb_to_linear(x);
        }
        apply3x3(m, v);
        for (int c = 0; c < 3; c++) {
            const float y = linear_to_srgb(v[c]) * 255.0f + 0.5f;
            rgb[3*i + c] = (uint8_t)(y < 0.0f ? 0.0f : (y > 255.0f ? 255.0f : y));
        }
    }
}

// The inverse of to_srgb_inplace.
inline void from_srgb_inplace(uint8_t* rgb, size_t n,
                              const std::string& gamut, bool is_linear) {
    if (is_identity(gamut, is_linear)) return;
    const Mat3 m = invert3x3(gamut_to_rec709(gamut));
    for (size_t i = 0; i < n; i++) {
        float v[3];
        for (int c = 0; c < 3; c++)
            v[c] = srgb_to_linear(rgb[3*i + c] * (1.0f / 255.0f));
        apply3x3(m, v);
        for (int c = 0; c < 3; c++) {
            const float y = (is_linear ? v[c] : linear_to_srgb(v[c])) * 255.0f + 0.5f;
            rgb[3*i + c] = (uint8_t)(y < 0.0f ? 0.0f : (y > 255.0f ? 255.0f : y));
        }
    }
}

}  // namespace colorspace
