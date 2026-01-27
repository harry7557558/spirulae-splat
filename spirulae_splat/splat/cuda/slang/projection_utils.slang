#pragma once
// ^ not working for [CudaDeviceExport]

#ifndef _PROJECTION_UTILS_SLANG_INC

#define _PROJECTION_UTILS_SLANG_INC


static const float ALPHA_THRESHOLD = 1.0f/255.0f;

enum CameraModelType {
    PINHOLE = 0,
    FISHEYE = 1,
    ORTHO = 2
};

// https://www.shadertoy.com/view/4djSRW

#define fract(x) ((x)-floor(x))

[ForceInline]
float hash14(in float4 p4) {
    p4 = fract(p4 * float4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    return fract((p4.x + p4.y) * (p4.z + p4.w));
}

[ForceInline]
float2 hash24(in float4 p4) {
    p4 = fract(p4 * float4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    return fract((p4.xz + p4.yw) * p4.wx);
}

[ForceInline]
float4 hash44(in float4 p4) {
    p4 = fract(p4 * float4(.1031, .1030, .0973, .1099));
    p4 += dot(p4, p4.wzxy + 33.33);
    return fract((p4.xxyz + p4.yzzw) * p4.zywx);
}

#undef fract

[ForceInline]
float2 box_muller(in float2 uv) {
    return sqrt(-2.0f * log(1.0 - uv.x)) * float2(
        cos(6.28318530718 * uv.y),
        sin(6.28318530718 * uv.y)
    );
}

[ForceInline]
float2 randn2(in float4 seed) {
    float2 val = hash24(seed);
    return box_muller(val);
}

[ForceInline]
float3 randn3(in float4 seed) {
    float4 val = hash44(seed);
    return float3(
        box_muller(val.xz),
        box_muller(val.yw).x
    );
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
float3x3 normalized_quat_to_rotmat(const float4 quat) {
    float w = quat[0], x = quat[1], y = quat[2], z = quat[3];
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
    float3x3 R = normalized_quat_to_rotmat(quat);
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
    float3x3 R = normalized_quat_to_rotmat(quat);
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

[ForceInline]
[Differentiable]
void sh_coeffs_to_color(
    const float3 viewdir,
    const float3 coeffs[16],
    out float3 colors
) {
    colors = 0.2820947917738781f * coeffs[0];

    float norm = sqrt(
        viewdir.x * viewdir.x + viewdir.y * viewdir.y + viewdir.z * viewdir.z
    );
    float x = viewdir.x / norm;
    float y = viewdir.y / norm;
    float z = viewdir.z / norm;

    float fTmp0A = 0.48860251190292f;
    colors += fTmp0A *
            (-y * coeffs[1] +
            z * coeffs[2] -
            x * coeffs[3]);
    float z2 = z * z;

    float fTmp0B = -1.092548430592079f * z;
    float fTmp1A = 0.5462742152960395f;
    float fC1 = x * x - y * y;
    float fS1 = 2.f * x * y;
    float pSH6 = (0.9461746957575601f * z2 - 0.3153915652525201f);
    float pSH7 = fTmp0B * x;
    float pSH5 = fTmp0B * y;
    float pSH8 = fTmp1A * fC1;
    float pSH4 = fTmp1A * fS1;
    colors +=
        pSH4 * coeffs[4] + pSH5 * coeffs[5] +
        pSH6 * coeffs[6] + pSH7 * coeffs[7] +
        pSH8 * coeffs[8];

    float fTmp0C = -2.285228997322329f * z2 + 0.4570457994644658f;
    float fTmp1B = 1.445305721320277f * z;
    float fTmp2A = -0.5900435899266435f;
    float fC2 = x * fC1 - y * fS1;
    float fS2 = x * fS1 + y * fC1;
    float pSH12 = z * (1.865881662950577f * z2 - 1.119528997770346f);
    float pSH13 = fTmp0C * x;
    float pSH11 = fTmp0C * y;
    float pSH14 = fTmp1B * fC1;
    float pSH10 = fTmp1B * fS1;
    float pSH15 = fTmp2A * fC2;
    float pSH9  = fTmp2A * fS2;
    colors += pSH9  * coeffs[9] +
            pSH10 * coeffs[10] +
            pSH11 * coeffs[11] +
            pSH12 * coeffs[12] +
            pSH13 * coeffs[13] +
            pSH14 * coeffs[14] +
            pSH15 * coeffs[15];

    // colors = max(colors+0.5f, 0.0f);
}

[ForceInline]
[Differentiable]
void sv_coeffs_to_color<let num_sv : int>(
    float3 viewdir,
    const float3 sv_sites[num_sv],
    const float3 sv_colors[num_sv],
    out float3 color
) {
    viewdir *= rsqrt(dot(viewdir, viewdir));

    float maxdot = dot(sv_sites[0], viewdir);
    [ForceUnroll] for (int i = 1; i < num_sv; i++) {
        maxdot = dot(sv_sites[i], viewdir);
    }

    float total_w = 0.0;
    float3 total_c = float3(0.0);
    [ForceUnroll] for (int i = 0; i < num_sv; i++) {
        float w = exp(dot(sv_sites[i], viewdir) - maxdot);
        total_w += w;
        total_c += w * sv_colors[i];
    }

    color = total_c / total_w;
}

[ForceInline]
[Differentiable]
float2 camera_distortion(float2 uv, no_diff float[10] dist_coeffs) {
    float k1 = dist_coeffs[0],
          k2 = dist_coeffs[1],
          k3 = dist_coeffs[2],
          k4 = dist_coeffs[3],
          p1 = dist_coeffs[4],
          p2 = dist_coeffs[5],
          sx1 = dist_coeffs[6],
          sy1 = dist_coeffs[7],
          b1 = dist_coeffs[8],
          b2 = dist_coeffs[9];
    float u = uv.x, v = uv.y;
    float r2 = u*u+v*v;
    float radial = 1.0 + r2*(k1 + r2*(k2 + r2*(k3 + r2*k4)));
    float du = 2.0*p1*u*v + p2*(r2+2.0*u*u) + sx1*r2;
    float dv = 2.0*p2*u*v + p1*(r2+2.0*v*v) + sy1*r2;
    uv = uv * radial + float2(du, dv);
    return uv + float2(b1 * uv.x + b2 * uv.y, 0.0);
}

// [ForceInline]
[Differentiable]
float2x2 camera_distortion_jac(float2 uv, no_diff float[10] dist_coeffs) {
    DifferentialPair<float2> dp_uv_x = diffPair(uv);
    bwd_diff(camera_distortion)(dp_uv_x, dist_coeffs, float2(1, 0));
    DifferentialPair<float2> dp_uv_y = diffPair(uv);
    bwd_diff(camera_distortion)(dp_uv_y, dist_coeffs, float2(0, 1));
    return float2x2(dp_uv_x.d, dp_uv_y.d);
}

[ForceInline]
[CudaDeviceExport]
bool is_valid_distortion(no_diff float2 uv, no_diff float[10] dist_coeffs) {
    // TODO: doesn't always work for more than 2 radial coeffs
    float2x2 J = camera_distortion_jac(uv, dist_coeffs);
    return min(determinant(J), min(J[0][0], J[1][1])) > 0.0;
}

// [ForceInline]
[Differentiable]
float valid_distortion_log_barrier(float2 uv, no_diff float[10] dist_coeffs) {
    // TODO: doesn't always work for more than 2 radial coeffs
    float2x2 J = camera_distortion_jac(uv, dist_coeffs);
    float val0 = min(determinant(J), min(J[0][0], J[1][1])); // positive if inside
    return -log(val0);
}

[ForceInline]
[Differentiable]
float2 valid_distortion_log_barrier_grad(float2 uv, no_diff float[10] dist_coeffs) {
    DifferentialPair<float2> dp_uv = diffPair(uv);
    bwd_diff(valid_distortion_log_barrier)(dp_uv, dist_coeffs, 1.0f);
    return dp_uv.d;
}

[ForceInline]
[Differentiable]
float2x2 valid_distortion_log_barrier_hess(float2 uv, no_diff float[10] dist_coeffs) {
    DifferentialPair<float2> dp_uv_x = diffPair(uv);
    bwd_diff(valid_distortion_log_barrier_grad)(dp_uv_x, dist_coeffs, float2(1.0f, 0.0f));
    DifferentialPair<float2> dp_uv_y = diffPair(uv);
    bwd_diff(valid_distortion_log_barrier_grad)(dp_uv_y, dist_coeffs, float2(0.0f, 1.0f));
    return float2x2(dp_uv_x.d, dp_uv_y.d);
}

[ForceInline]
[Differentiable]
bool persp_proj<let assume_valid : bool>(
    float3 p_view, no_diff float4 intrins, no_diff float[10] dist_coeffs,
    out float2 uv
) {
    uv = p_view.xy / p_view.z;
    if (!assume_valid && (p_view.z < 0.0 || !is_valid_distortion(uv, dist_coeffs)))
        return false;
    uv = camera_distortion(uv, dist_coeffs);
    float fx = intrins.x, fy = intrins.y, cx = intrins.z, cy = intrins.w;
    uv = float2(
        fx * uv.x + cx,
        fy * uv.y + cy
    );
    return true;
}

[ForceInline]
[Differentiable]
float2x3 persp_proj_jac(float3 p_view, no_diff float4 intrins, no_diff float[10] dist_coeffs) {
    float2x3 J;
    [ForceUnroll]
    for (int i = 0; i < 2; ++i) {
        float2 seed = 0.0; seed[i] = 1.0;
        DifferentialPair<float3> dp_p_view = diffPair(p_view);
        bwd_diff(persp_proj<true>)(dp_p_view, intrins, dist_coeffs, seed);
        J[i] = dp_p_view.d;
    }
    return J;
}

[ForceInline]
[Differentiable]
float2 fisheye_proj(float3 p_view) {
    float r = length(p_view.xy);
    float theta = atan2(r, p_view.z);
    float k = theta < 1e-3 ? (1.0 - theta*theta/3.0) / p_view.z : theta / r;
    return p_view.xy * k;
}

[ForceInline]
[Differentiable]
bool fisheye_proj<let assume_valid : bool>(
    float3 p_view, no_diff float4 intrins, no_diff float[10] dist_coeffs,
    out float2 uv
) {
    uv = fisheye_proj(p_view);
    if (!assume_valid && !is_valid_distortion(uv, dist_coeffs))
        return false;
    uv = camera_distortion(uv, dist_coeffs);
    float fx = intrins.x, fy = intrins.y, cx = intrins.z, cy = intrins.w;
    uv = float2(
        fx * uv.x + cx,
        fy * uv.y + cy
    );
    return true;
}

[ForceInline]
[Differentiable]
float2x3 fisheye_proj_jac(float3 p_view, no_diff float4 intrins, no_diff float[10] dist_coeffs) {
    float2x3 J;
    [ForceUnroll]
    for (int i = 0; i < 2; ++i) {
        float2 seed = 0.0; seed[i] = 1.0;
        DifferentialPair<float3> dp_p_view = diffPair(p_view);
        bwd_diff(fisheye_proj<true>)(dp_p_view, intrins, dist_coeffs, seed);
        J[i] = dp_p_view.d;
    }
    return J;
}

[ForceInline]
[Differentiable]
float3 fisheye_unproj(float2 uv) {
    float theta = length(uv);
    return float3(
        uv * (theta < 1e-3 ? 1.0-theta*theta/6.0 : sin(theta) / theta),
        cos(theta)
    );
}


[ForceInline]
[Differentiable]
float2 _solve2x2(float2x2 A, float2 b) {
    float inv_det = 1.0 / (A[0][0]*A[1][1] - A[0][1]*A[1][0]);
    return float2(
        ( b.x*A[1][1] - b.y*A[0][1]) * inv_det,
        (-b.x*A[1][0] + b.y*A[0][0]) * inv_det
    );
}

// find uv' s.t. distort(uv') = uv
bool undistort_point(
    no_diff float2 uv, no_diff float[10] dist_coeffs,
    int maxiter, out float2 uv_undist
) {
    float2 q = uv;
    float2 prev_q = q;

    for (int i = 0; i < maxiter; i++) {
        float2 qd = camera_distortion(q, dist_coeffs);
        float2 r = qd - uv;
        float err = dot(r, r);

        // bool use_newton = err < 1e-2;
        bool use_newton = true;

        if (use_newton) {
            float2x2 J = camera_distortion_jac(q, dist_coeffs);
            q -= _solve2x2(J, r);
        } else {
            // float step = 0.5;
            float step = clamp(0.5 + 0.5 / (1.0 + dot(r, r) * 10.0), 0.1, 0.9);
            q = q - step * r;
        }

        prev_q = q;
    }

    uv_undist = q;
    // TODO: better way to tell so
    return is_valid_distortion(uv_undist, dist_coeffs) &&
           length(camera_distortion(uv_undist, dist_coeffs) - uv) < 1e-2;
}

// find uv' that minimizes (distort(uv') - uv)^2 s.t. is_valid_distortion_func(uv') < 0
bool undistort_point_clipped(
    no_diff float2 uv, no_diff float[10] dist_coeffs,
    int maxiter, out float2 q
) {
    q = uv;
    for (int i = 0; i < 4; i++) {
        if (is_valid_distortion(q, dist_coeffs))
            break;
        q *= 0.618f;
    }
    if (!is_valid_distortion(q, dist_coeffs))
        return false;

    float weight_b = 0.0f;
    for (int i = 0; i < maxiter; i++) {
        float2 qd = camera_distortion(q, dist_coeffs);
        float2 r = qd - uv;
        if (length(r) < 1e-4)
            break;

        float2x2 J = camera_distortion_jac(q, dist_coeffs);
        float2 grad_f = mul(transpose(J), r);
        float2x2 hess_f = mul(transpose(J), J);

        float2 grad_b = valid_distortion_log_barrier_grad(q, dist_coeffs);
        float2x2 hess_b = valid_distortion_log_barrier_hess(q, dist_coeffs);

        if (i == 0) {
            weight_b = 10.0f * length(grad_f) / (length(grad_b) + 1e-6f);
        }

        bool valid_found = false;
        for (int j = 0; j < 8; j++) {
            float2 grad = grad_f + weight_b * grad_b;
            float2x2 hess = hess_f + weight_b * hess_b;
            float2 dq = _solve2x2(J, r);
            float2 q1 = q - dq;
            if (is_valid_distortion(q1, dist_coeffs)) {
                valid_found = true;
                q = q1;
                weight_b *= 0.25f;
                break;
            }
            else {
                weight_b *= 1.618f;
            }
        }
        if (!valid_found)
            return false;
    }

    // TODO: better way to tell so
    return true;
    // return length(camera_distortion(q, dist_coeffs) - uv) < 1e-2;
}

[CudaDeviceExport]
[Differentiable]
float2 distort_point(
    no_diff float2 uv, bool is_fisheye,
    no_diff float[10] dist_coeffs
) {
    if (is_fisheye) {
        float r = length(uv);
        float theta = atan(r);
        uv *= (r < 1e-3 ? 1.0-r*r/3.0 : theta/r);
    }
    return camera_distortion(uv, dist_coeffs);
}

[CudaDeviceExport]
[Differentiable]
bool undistort_point(
    no_diff float2 uv, bool is_fisheye,
    no_diff float[10] dist_coeffs,
    out float2 uv_undist
) {
    if (!no_diff undistort_point(uv, dist_coeffs, 8, uv))
        return false;
    float3 raydir = is_fisheye ? fisheye_unproj(uv) : float3(uv, 1.0);
    uv_undist = raydir.xy / max(raydir.z, 1e-12);
    return true;
}

[CudaDeviceExport]
bool unproject_point(
    no_diff float2 uv, bool is_fisheye,
    no_diff float[10] dist_coeffs,
    out float3 raydir
) {
    raydir = float3(0);
    // if (!undistort_point_clipped(uv, dist_coeffs, 8, uv))
    if (!undistort_point(uv, dist_coeffs, 8, uv))
        return false;
    raydir = is_fisheye ? fisheye_unproj(uv) : float3(uv, 1.0);
    return true;
}

[CudaDeviceExport]
[Differentiable]
bool generate_ray(
    no_diff float2 uv, bool is_fisheye,
    no_diff float[10] dist_coeffs,
    out float3 raydir
) {
    if (!no_diff undistort_point(uv, dist_coeffs, 8, uv))
        { raydir = float3(0); return false; }
    raydir = is_fisheye ? fisheye_unproj(uv) : float3(uv, 1.0);
    raydir = normalize(raydir);
    return true;
}

[CudaDeviceExport]
[Differentiable]
float3 transform_ray_o(float3x3 R, float3 t) {
    return -mul(t, R);
}

[CudaDeviceExport]
[Differentiable]
float3 transform_ray_d(float3x3 R, float3 raydir) {
    return mul(raydir, R);
}

[CudaDeviceExport]
float3 undo_transform_ray_d(float3x3 R, float3 raydir) {
    return mul(raydir, transpose(R));
}

[CudaDeviceExport]
void transform_ray_o_vjp(
    float3x3 R, float3 t, float3 v_ray_o,
    out float3x3 v_R, out float3 v_t
) {
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_t = diffPair(t);
    bwd_diff(transform_ray_o)(dp_R, dp_t, v_ray_o);
    v_R = dp_R.d;
    v_t = dp_t.d;
}

[CudaDeviceExport]
void transform_ray_d_vjp(
    float3x3 R, float3 raydir, float3 v_ray_d,
    out float3x3 v_R, out float3 v_raydir
) {
    DifferentialPair<float3x3> dp_R = diffPair(R);
    DifferentialPair<float3> dp_raydir = diffPair(raydir);
    bwd_diff(transform_ray_d)(
        dp_R, dp_raydir,
        v_ray_d
    );
    v_R = dp_R.d;
    v_raydir = dp_raydir.d;
}

[CudaDeviceExport]
[ForceInline]
[Differentiable]
void map_opaque_triangle(
    float3 mean, float4 quat, float3 scale,
    out float3 vert0, out float3 vert1, out float3 vert2
) {
    float sx = exp(scale.x),
          sy = exp(scale.y),
          sz = scale.z - 0.5 * (scale.x + scale.y);
    float3x3 M = normalized_quat_to_rotmat(quat);
    vert0 = mul(M, float3(sx, 0, 0)) + mean;
    vert1 = mul(M, float3(sx * (-0.5 + sz), sy, 0)) + mean;
    vert2 = mul(M, float3(sx * (-0.5 - sz), -sy, 0)) + mean;
}

[ForceInline]
float3x3 outerProduct(float3 u) {
    return float3x3(u.x * u, u.y * u, u.z * u);
}

#endif  // _PROJECTION_UTILS_SLANG_INC
