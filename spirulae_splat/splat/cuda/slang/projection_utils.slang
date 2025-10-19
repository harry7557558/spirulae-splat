#pragma once

static const float ALPHA_THRESHOLD = 1.0f/255.0f;

enum CameraModelType {
    PINHOLE = 0,
    FISHEYE = 1,
    ORTHO = 2
};

[CudaDeviceExport]
[ForceInline]
[Differentiable]
float3x3 quat_to_rotmat(const float4 quat) {
    float w = quat[0], x = quat[1], y = quat[2], z = quat[3];
    // normalize
    float inv_norm = rsqrt(x * x + y * y + z * z + w * w);
    x *= inv_norm;
    y *= inv_norm;
    z *= inv_norm;
    w *= inv_norm;
    float x2 = x * x, y2 = y * y, z2 = z * z;
    float xy = x * y, xz = x * z, yz = y * z;
    float wx = w * x, wy = w * y, wz = w * z;
    return transpose(float3x3(
        (1.f - 2.f * (y2 + z2)),
        (2.f * (xy + wz)),
        (2.f * (xz - wy)), // 1st col
        (2.f * (xy - wz)),
        (1.f - 2.f * (x2 + z2)),
        (2.f * (yz + wx)), // 2nd col
        (2.f * (xz + wy)),
        (2.f * (yz - wx)),
        (1.f - 2.f * (x2 + y2)) // 3rd col
    ));
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void posW2C(
    const float3x3 R,
    const float3 t,
    const float3 pW, // Input position in world coordinates
    out float3 pC       // Output position in camera coordinates
) {
    pC = mul(R, pW) + t;
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void covarW2C(
    const float3x3 R,
    const float3x3 covarW, // Input covariance matrix in world coordinates
    out float3x3 covarC       // Output covariance matrix in camera coordinates
) {
    covarC = mul(mul(R, covarW), transpose(R));
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void quat_scale_to_covar(
    const float4 quat,
    const float3 scale,
    out float3x3 covar
) {
    float3x3 R = quat_to_rotmat(quat);
    // C = R * S * S * Rt
    float3x3 S = float3x3(scale[0], 0.f, 0.f, 0.f, scale[1], 0.f, 0.f, 0.f, scale[2]);
    float3x3 M = mul(R, S);
    covar = mul(M, transpose(M));
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void quat_scale_to_sqrt_covar(
    const float4 quat,
    const float3 scale,
    out float3x3 M
) {
    float3x3 R = quat_to_rotmat(quat);
    float3x3 S = float3x3(scale[0], 0.f, 0.f, 0.f, scale[1], 0.f, 0.f, 0.f, scale[2]);
    M = mul(R, S);
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
float2x2 inverse(float2x2 m){
    float invdet = 1.0f / (m[0][0] * m[1][1] - m[0][1] * m[1][0]);
    return {
        m[1][1] * invdet,
        -m[0][1] * invdet,
        -m[1][0] * invdet,
        m[0][0] * invdet
    };
}

[ForceInline]
[Differentiable]
float det(float2 u, float2 v) {
    return u.x*v.y - u.y*v.x;
}
