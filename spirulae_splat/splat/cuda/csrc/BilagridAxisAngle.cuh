#pragma once

#include <cuda_runtime.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626f
#endif


__device__ __forceinline__ float3 cross3(const float3 &a, const float3 &b) {
    return make_float3(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}

__device__ __forceinline__ float dot3(const float3 &a, const float3 &b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

__device__ __forceinline__ float3 mul3(const float3 &a, float s) {
    return make_float3(a.x*s, a.y*s, a.z*s);
}

__device__ __forceinline__ float3 add3(const float3 &a, const float3 &b) {
    return make_float3(a.x+b.x, a.y+b.y, a.z+b.z);
}


// ------------------------- FORWARD -------------------------
inline __device__ float3 axis_angle_rotate_fwd(float3 w, float3 v)
{
    const float eps = 1e-5f;

    float theta2 = dot3(w, w);
    float theta  = sqrtf(theta2 + 1e-20f);

    float c, s, t;
    float3 u;

    if (theta < eps) {
        // Taylor expansions
        float inv6 = 1.f / 6.f;
        c = 1.f - 0.5f * theta2;
        s = theta - inv6 * theta * theta2;
        t = 0.5f * theta2;

        // u ~ w / theta but stable: if theta small, use first-order approx u ~ w / (theta + eps)
        float inv_theta = 1.f / (theta + 1e-20f);
        u = mul3(w, inv_theta);
    } else {
        c = __cosf(theta);
        s = __sinf(theta);
        t = 1.f - c;
        float inv_theta = 1.f / theta;
        u = mul3(w, inv_theta);
    }

    float3 uv  = cross3(u, v);
    float uvd  = dot3(u, v);

    float3 term1 = mul3(v, c);
    float3 term2 = mul3(uv, s);
    float3 term3 = mul3(u, uvd * t);

    return add3(add3(term1, term2), term3);
}

// ------------------------- BACKWARD -------------------------

inline __device__ void axis_angle_rotate_bwd(
    const float3 w,
    const float3 v,
    const float3 grad_out,
    float3& grad_w,
    float3& grad_v)
{
    const float eps = 1e-5f;

    const float theta2 = dot3(w, w);
    const float theta  = sqrtf(theta2 + 1e-20f);

    float A, B, dA_dtheta, dB_dtheta;

    if (theta < eps) {
        // Stable series expansions around theta = 0
        const float t2 = theta2;
        const float t4 = t2 * t2;

        A = 1.0f - t2 * (1.0f / 6.0f) + t4 * (1.0f / 120.0f);
        B = 0.5f - t2 * (1.0f / 24.0f) + t4 * (1.0f / 720.0f);

        dA_dtheta = -theta * (1.0f / 3.0f) + theta * t2 * (1.0f / 30.0f);
        dB_dtheta = -theta * (1.0f / 12.0f) + theta * t2 * (1.0f / 180.0f);
    } else {
        const float s = __sinf(theta);
        const float c = __cosf(theta);

        A = s / theta;
        B = (1.0f - c) / theta2;

        dA_dtheta = (theta * c - s) / theta2;
        dB_dtheta = (theta * s - 2.0f * (1.0f - c)) / (theta2 * theta);
    }

    const float3 c1 = cross3(w, v);     // w x v
    const float3 c2 = cross3(w, c1);    // w x (w x v)

    // TODO: fix nonzero gradient when grad_out is zero

    const float g_c1 = dot3(grad_out, c1);
    const float g_c2 = dot3(grad_out, c2);

    const float inv_theta = 1.0f / (theta + 1e-20f);

    // Chain term from A(theta), B(theta) through theta = ||w||
    const float alpha =
        (dA_dtheta * g_c1 + dB_dtheta * g_c2) * inv_theta;

    // grad_w =
    //   w * alpha
    // + A * (v x grad_out)
    // + B * (c1 x grad_out)
    // + B * [ dot(w, v) grad_out - dot(grad_out, v) w ]
    const float w_dot_v = dot3(w, v);
    const float g_dot_v = dot3(grad_out, v);

    const float3 term_theta = mul3(w, alpha);
    const float3 term_a     = mul3(cross3(v, grad_out), A);
    const float3 term_b     = mul3(cross3(c1, grad_out), B);
    const float3 term_c     = mul3(
        add3(mul3(grad_out, w_dot_v), mul3(w, -g_dot_v)),
        B
    );

    grad_w = add3(term_theta, add3(term_a, add3(term_b, term_c)));

    // grad_v = R(-w) * grad_out
    grad_v = axis_angle_rotate_fwd(mul3(w, -1.0f), grad_out);
}
