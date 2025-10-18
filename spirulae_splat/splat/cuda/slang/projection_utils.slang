
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
void persp_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    const uint width,
    const uint height,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float tan_fovx = 0.5f * width / fx;
    float tan_fovy = 0.5f * height / fy;
    float lim_x_pos = (width - cx) / fx + 0.3f * tan_fovx;
    float lim_x_neg = cx / fx + 0.3f * tan_fovx;
    float lim_y_pos = (height - cy) / fy + 0.3f * tan_fovy;
    float lim_y_neg = cy / fy + 0.3f * tan_fovy;

    float rz = 1.f / z;
    float rz2 = rz * rz;
    float tx = z * min(lim_x_pos, max(-lim_x_neg, x * rz));
    float ty = z * min(lim_y_pos, max(-lim_y_neg, y * rz));

    float2x3 J = float2x3(
        fx * rz, 0.f, -fx * tx * rz2,  // 1st row
        0.f, fy * rz, -fy * ty * rz2  // 2nd row
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
    mean2d = float2(fx * x * rz + cx, fy * y * rz + cy);
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void persp_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float rz = 1.f / z;
    float rz2 = rz * rz;
    float tx = x, ty = y;

    float2x3 J = float2x3(
        fx * rz, 0.f, -fx * tx * rz2,  // 1st row
        0.f, fy * rz, -fy * ty * rz2  // 2nd row
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
    mean2d = float2(fx * x * rz + cx, fy * y * rz + cy);
}



[CudaDeviceExport]
[ForceInline]
[Differentiable]
void fisheye_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float eps = 0.0000001f;
    float xy_len = length(float2(x, y)) + eps;
    float theta = atan2(xy_len, z + eps);
    mean2d = float2(x * fx * theta / xy_len + cx, y * fy * theta / xy_len + cy);

    float x2 = x * x + eps;
    float y2 = y * y;
    float xy = x * y;
    float x2y2 = x2 + y2;
    float x2y2z2_inv = 1.f / (x2y2 + z * z);

    float b = atan2(xy_len, z) / xy_len / x2y2;
    float a = z * x2y2z2_inv / (x2y2);
    float2x3 J = float2x3(
        fx * (x2 * a + y2 * b), fx * xy * (a - b), -fx * x * x2y2z2_inv,
        fy * xy * (a - b), fy * (y2 * a + x2 * b), -fy * y * x2y2z2_inv
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
}


[CudaDeviceExport]
[ForceInline]
[Differentiable]
void ortho_proj_3dgs(
    // inputs
    const float3 mean3d,
    const float3x3 cov3d,
    const float fx,
    const float fy,
    const float cx,
    const float cy,
    // outputs
    out float2x2 cov2d,
    out float2 mean2d
) {
    float x = mean3d[0], y = mean3d[1], z = mean3d[2];

    float2x3 J = float2x3(
        fx, 0.f, 0.f,
        0.f, fy, 0.f
    );
    cov2d = mul(mul(J, cov3d), transpose(J));;
    mean2d = float2(fx * x + cx, fy * y + cy);
}



[CudaDeviceExport]
[ForceInline]
[Differentiable]
float add_blur(const float eps2d, inout float2x2 covar, out float compensation) {
    float det_orig = covar[0][0] * covar[1][1] - covar[0][1] * covar[1][0];
    covar[0][0] += eps2d;
    covar[1][1] += eps2d;
    float det_blur = covar[0][0] * covar[1][1] - covar[0][1] * covar[1][0];
    compensation = sqrt(max(0.f, det_orig / det_blur));
    return det_blur;
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
