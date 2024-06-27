#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include "config.h"
#include "glm/glm/glm.hpp"

namespace cg = cooperative_groups;


inline __device__ float bessel_j0(const float x) {
    const float A01 = -0.22334141550835468f;
    const float A02 = -0.18054514613169334f;
    const float A03 = 0.047655819492146555f;
    const float A04 = -0.0024383224605644127f;
    const float B01 = -0.2234104067240744f;
    const float B02 = 0.06997488447829185f;
    const float B03 = -0.009545108494228093f;
    const float B04 = 0.0011020589286710896f;
    const float C01 = -0.12645960753014673f;
    const float C02 = -0.046982403187758204f;
    const float D01 = 0.12546821167569544f;
    const float D02 = -0.08179017949118085f;
    if (x < 4.2f) {
        return (1.0f + x * (A01 + x * (A02 + x * (A03 + x * A04)))) / 
               (1.0f + x * (B01 + x * (B02 + x * (B03 + x * B04))));
    } else {
        float inv_x = 1.0f / x;
        return sqrtf(inv_x / M_PI) * 
               ((1.0f + inv_x * (C01 + inv_x * C02)) * cosf(x) + 
                (1.0f + inv_x * (D01 + inv_x * D02)) * sinf(x));
    }
}

inline __device__ float bessel_j1(const float x) {
    const float A11 = 0.5000532505934573f;
    const float A12 = -0.07119189008625842f;
    const float A13 = -0.03516544965310418f;
    const float A14 = 0.005137712441402014f;
    const float B11 = -0.14169754004287138f;
    const float B12 = 0.05321374041943971f;
    const float B13 = -0.006074191793869077f;
    const float B14 = 0.0008890431150018836f;
    const float C11 = -0.37578090667550257f;
    const float C12 = 0.13415846822338878f;
    const float D11 = 0.3771741874195154f;
    const float D12 = 0.08328593955487182f;
    if (x < 4.7f) {
        return (x * (A11 + x * (A12 + x * (A13 + x * A14)))) / 
               (1.0f + x * (B11 + x * (B12 + x * (B13 + x * B14))));
    } else {
        float inv_x = 1.0f / x;
        return -sqrtf(inv_x / M_PI) * 
               ((1.0f + inv_x * (C11 + inv_x * C12)) * cosf(x) - 
                (1.0f + inv_x * (D11 + inv_x * D12)) * sinf(x));
    }
}

inline __device__ float bessel_j2(const float x) {
    const float A22 = 0.12515704623184004f;
    const float A23 = -0.029939064720165425f;
    const float A24 = 0.0010672441356451344f;
    const float B21 = -0.23414231622686957f;
    const float B22 = 0.08321793931659045f;
    const float B23 = -0.012670220025970099f;
    const float B24 = 0.0015767563111494629f;
    const float C21 = 1.874412399273724f;
    const float C22 = -0.8101992991186221f;
    const float C23 = 0.3015954731134034f;
    const float D21 = -1.871274305839778f;
    const float D22 = -0.8861009908575821f;
    if (x < 4.0f) {
        return (x * x * (A22 + x * (A23 + x * A24))) / 
               (1.0f + x * (B21 + x * (B22 + x * (B23 + x * B24))));
    } else {
        float inv_x = 1.0f / x;
        return -sqrtf(inv_x / M_PI) * 
               ((1.0f + inv_x * (C21 + inv_x * (C22 + inv_x * C23))) * cosf(x) + 
                (1.0f + inv_x * (D21 + inv_x * D22)) * sinf(x));
    }
}

inline __device__ float bessel_j3(const float x) {
    const float A33 = 0.020910819133348472f;
    const float A34 = -0.005072094376038419f;
    const float A35 = 0.0002802765938927514f;
    const float B31 = -0.2330616229268751f;
    const float B32 = 0.06455328873550212f;
    const float B33 = -0.008312028977714348f;
    const float B34 = 0.0007466861514973682f;
    const float C31 = -4.376749965939438f;
    const float C32 = -7.327544311795212f;
    const float C33 = 2.8595505732173425f;
    const float D31 = 4.374149309521666f;
    const float D32 = -7.3507982430716545f;
    const float D33 = -3.7324735035522663f;
    if (x < 5.0f) {
        return (x * x * x * (A33 + x * (A34 + x * A35))) / 
               (1.0f + x * (B31 + x * (B32 + x * (B33 + x * B34))));
    } else {
        float inv_x = 1.0f / x;
        return sqrtf(inv_x / M_PI) * 
               ((1.0f + inv_x * (C31 + inv_x * (C32 + inv_x * C33))) * cosf(x) - 
                (1.0f + inv_x * (D31 + inv_x * (D32 + inv_x * D33))) * sinf(x));
    }
}

inline __device__ float bessel_j4(const float x) {
    const float A44 = 0.002644492060608329f;
    const float A45 = -0.0006004507528516955f;
    const float A46 = 0.00003320308950860871f;
    const float B41 = -0.20296731043978247f;
    const float B42 = 0.04338600070178919f;
    const float B43 = -0.0035265908540099847f;
    const float B44 = 0.00013712907221840123f;
    const float B45 = 0.000012746991211123013f;
    const float C41 = 7.881701792737443f;
    const float C42 = -27.37611266073206f;
    const float C43 = -39.62023054032f;
    const float D41 = -7.865445288974054f;
    const float D42 = -27.479039704046176f;
    const float D43 = 49.286435632834696f;
    if (x < 8.0f) {
        return (x * x * x * x * (A44 + x * (A45 + x * A46))) / 
               (1.0f + x * (B41 + x * (B42 + x * (B43 + x * (B44 + x * B45)))));
    } else {
        float inv_x = 1.0f / x;
        return sqrtf(inv_x / M_PI) * 
               ((1.0f + inv_x * (C41 + inv_x * (C42 + inv_x * C43))) * cosf(x) + 
                (1.0f + inv_x * (D41 + inv_x * (D42 + inv_x * D43))) * sinf(x));
    }
}



inline __device__ float bessel_j(
    unsigned m, float x
) {
    // return jnf((int)m, x);
    if (m == 0)
        return bessel_j0(x);
    if (m == 1)
        return bessel_j1(x);
    if (m == 2)
        return bessel_j2(x);
    if (m == 3)
        return bessel_j3(x);
    if (m == 4)
        return bessel_j4(x);
    return 0.0f;
}


inline __device__ glm::vec3 ch_coeffs_to_color(
    const unsigned degree_r,
    const unsigned degree_phi,
    const glm::vec3 *coeffs,
    const glm::vec2 &uv
) {
    const float r = hypot(uv.x, uv.y);
    const float phi = atan2(uv.y, uv.x);
    const float pi = 3.14159265358979;

    assert(degree_phi <= 4);

    int idx = 0;
    glm::vec3 color(0.0f);

    for (unsigned k = 1; k <= degree_r; k++) {
        // m == 0
        float w = bessel_j(0, k*pi*r);
        color += w * coeffs[idx];
        idx += 1;
        // m > 0
        for (unsigned m = 1; m <= degree_phi; m++) {
            float wr = bessel_j(m, k*pi*r);
            float wc = wr * cosf(m*phi);
            float ws = wr * sinf(m*phi);
            color += wc * coeffs[idx];
            idx += 1;
            color += ws * coeffs[idx];
            idx += 1;
        }
    }

    return color;
}


inline __device__ void ch_coeffs_to_color_vjp(
    const unsigned degree_r,
    const unsigned degree_phi,
    const glm::vec3 *coeffs,
    const glm::vec2 &uv,
    const glm::vec3 v_color,
    glm::vec3 &color,
    glm::vec3 *v_coeffs,
    float &v_ch_coeff_abs,
    glm::vec2 &v_uv
) {
    const float r = hypot(uv.x, uv.y);
    const float phi = atan2(uv.y, uv.x);
    const float pi = 3.14159265358979;

    const float inv_r = 1.0f / fmax(r, 1e-6f);
    float v_r = 0.f, v_phi = 0.f;

    assert(degree_phi <= 4);

    int idx = 0;
    color = glm::vec3(0.0f);
    v_ch_coeff_abs = 0.0f;

    // derivative of J:
    // J'(0,x) = -J(1,x)
    // J'(n,x) = J(n-1,x) - n/x J(n,x)
    // J'(n,x) = 1/2 J(n-1,x) - 1/2 J(n+1,x)

    for (unsigned k = 1; k <= degree_r; k++) {
        // m == 0
        float w = bessel_j(0, k*pi*r);
        float v_w = 0.f;
        {
            glm::vec3 coef = coeffs[idx];
            color += w * coef;
            v_w += glm::dot(coef, v_color);
            glm::vec3 v_coeff = w * v_color;
            v_coeffs[idx] = v_coeff;
            v_ch_coeff_abs += abs(v_coeff.x)+abs(v_coeff.y)+abs(v_coeff.z);
        }
        idx += 1;
        float w1 = bessel_j(1, k*pi*r);
        v_r += k*pi * (-w1) * v_w;
        // m > 0
        float prev_wr = w;
        for (unsigned m = 1; m <= degree_phi; m++) {
            float wr = m == 1 ? w1 : bessel_j(m, k*pi*r);
            float cos_m_phi = cosf(m*phi);
            float sin_m_phi = sinf(m*phi);
            float wc = wr * cos_m_phi;
            float ws = wr * sin_m_phi;
            float v_wc = 0.f, v_ws = 0.f;
            {
                glm::vec3 coef = coeffs[idx];
                color += wc * coef;
                v_wc += glm::dot(coef, v_color);
                glm::vec3 v_coeff = wc * v_color;
                v_coeffs[idx] = v_coeff;
                v_ch_coeff_abs += abs(v_coeff.x)+abs(v_coeff.y)+abs(v_coeff.z);
            }
            idx += 1;
            {
                glm::vec3 coef = coeffs[idx];
                color += ws * coef;
                v_ws += glm::dot(coef, v_color);
                glm::vec3 v_coeff = ws * v_color;
                v_coeffs[idx] = v_coeff;
                v_ch_coeff_abs += abs(v_coeff.x)+abs(v_coeff.y)+abs(v_coeff.z);
            }
            idx += 1;
            float v_wr = v_wc * cos_m_phi + v_ws * sin_m_phi;
            v_r += k*pi * (prev_wr - m/(k*pi)*wr*inv_r) * v_wr;
            v_phi += v_wc * wr * m * (-sin_m_phi);
            v_phi += v_ws * wr * m * cos_m_phi;
            prev_wr = wr;
        }
    }

    v_uv.x = v_r * cosf(phi) - v_phi * inv_r * sinf(phi);
    v_uv.y = v_r * sinf(phi) + v_phi * inv_r * cosf(phi);

    v_ch_coeff_abs /= (float)(3*idx);
}


// Optimized for backward of sigmoid of CH color
inline __device__ void ch_coeffs_to_color_sigmoid_vjp(
    const unsigned degree_r,
    const unsigned degree_phi,
    const glm::vec3 *coeffs,
    const glm::vec2 &uv,
    const glm::vec3 &v_color_sigmoid,
    glm::vec3 &color_sigmoid,
    glm::vec3 *v_coeffs,
    float &v_ch_coeff_abs,
    glm::vec2 &v_uv
) {
    const float r = hypot(uv.x, uv.y);
    const float phi = atan2(uv.y, uv.x);
    const float pi = 3.14159265358979;

    const float inv_r = 1.0f / fmax(r, 1e-6f);
    glm::vec3 v_r_vec(0.f), v_phi_vec(0.f);

    assert(degree_phi <= 4);

    int idx = 0;
    glm::vec3 color = glm::vec3(0.0f);

    // derivative of J:
    // J'(0,x) = -J(1,x)
    // J'(n,x) = J(n-1,x) - n/x J(n,x)
    // J'(n,x) = 1/2 J(n-1,x) - 1/2 J(n+1,x)

    for (unsigned k = 1; k <= degree_r; k++) {
        // m == 0
        float w = bessel_j(0, k*pi*r);
        glm::vec3 v_w(0.f);
        {
            glm::vec3 coef = coeffs[idx];
            color += w * coef;
            v_w += coef;
            v_coeffs[idx].x = w;
        }
        idx += 1;
        float w1 = bessel_j(1, k*pi*r);
        v_r_vec += k*pi * (-w1) * v_w;
        // m > 0
        float prev_wr = w;
        for (unsigned m = 1; m <= degree_phi; m++) {
            float wr = m == 1 ? w1 : bessel_j(m, k*pi*r);
            float cos_m_phi = cosf(m*phi);
            float sin_m_phi = sinf(m*phi);
            float wc = wr * cos_m_phi;
            float ws = wr * sin_m_phi;
            glm::vec3 v_wc(0.f), v_ws(0.f);
            {
                glm::vec3 coef = coeffs[idx];
                color += wc * coef;
                v_wc += coef;
                v_coeffs[idx].x = wc;
            }
            idx += 1;
            {
                glm::vec3 coef = coeffs[idx];
                color += ws * coef;
                v_ws += coef;
                v_coeffs[idx].x = ws;
            }
            idx += 1;
            glm::vec3 v_wr = v_wc * cos_m_phi + v_ws * sin_m_phi;
            v_r_vec += k*pi * (prev_wr - m/(k*pi)*wr*inv_r) * v_wr;
            v_phi_vec += v_wc * wr * float(m) * (-sin_m_phi);
            v_phi_vec += v_ws * wr * float(m) * cos_m_phi;
            prev_wr = wr;
        }
    }

    color_sigmoid = 1.0f / (1.0f + glm::exp(-color));
    glm::vec3 v_color = color_sigmoid * (1.0f-color_sigmoid) * v_color_sigmoid;

    v_ch_coeff_abs = 0.0f;
    for (int i = 0; i < idx; i++) {
        glm::vec3 v_coeff = v_coeffs[i].x * v_color;
        v_coeffs[i] = v_coeff;
        #if CH_ABSGRAD_REDUCE == 0
        v_ch_coeff_abs += abs(v_coeff.x) + abs(v_coeff.y) + abs(v_coeff.z);
        #elif CH_ABSGRAD_REDUCE == 1
        v_ch_coeff_abs = max(v_ch_coeff_abs, abs(v_coeff.x));
        v_ch_coeff_abs = max(v_ch_coeff_abs, abs(v_coeff.y));
        v_ch_coeff_abs = max(v_ch_coeff_abs, abs(v_coeff.z));
        #endif
    }
    #if CH_ABSGRAD_REDUCE == 0
    v_ch_coeff_abs /= (float)(3*idx);
    #endif

    float v_r = glm::dot(v_r_vec, v_color);
    float v_phi = glm::dot(v_phi_vec, v_color);
    v_uv.x = v_r * cosf(phi) - v_phi * inv_r * sinf(phi);
    v_uv.y = v_r * sinf(phi) + v_phi * inv_r * cosf(phi);
}
