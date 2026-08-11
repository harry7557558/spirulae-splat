#pragma once
// ^ not working for [CudaDeviceExport]

#ifndef _CAMERA_SOURCE_SLANG_INC

#define _CAMERA_SOURCE_SLANG_INC

// Forward projection of the camera models that no (camera model, distortion
// tier) pair represents. The parser fits those onto a pair and records the
// source here so the re-distort pass can resample from the TRUE source
// projection rather than from an approximation of it -- resampling through the
// fit would be a no-op.
//
// For a COLMAP model, `params` is COLMAP's own parameter array verbatim, so
// the host side copies `cameras.bin` through unchanged. Transcribed from
// colmap/src/colmap/sensor/models.h; the host mirror used by the fitter is
// srccam::project() in data/SourceCamera.cpp and the two must agree exactly or
// the fit is measured against a different camera than the images are resampled
// from.
//
// Kept out of projection_utils.slang deliberately: nothing on the training hot
// path calls this, and that header is included by every kernel TU. That also
// means the two fisheye mappings below are transcribed rather than shared;
// they must stay identical to fisheye_proj / equisolid_proj there.

static const int kSourceMaxParams = 16;

// COLMAP CameraModelId values. Only the ones with no exact tier appear.
static const int kSourceFOV = 7;
static const int kSourceRadTanThinPrismFisheye = 11;
static const int kSourceSimpleDivision = 12;
static const int kSourceDivision = 13;
static const int kSourceEUCM = 16;

// Ours, numbered from 1000 so COLMAP can keep appending to its own enum
// (18, 19, ...) without ever colliding.
static const int kSourceSkewed = 1000;

// kSourceSkewed's base camera, params[13].
static const int kSkewBasePerspective = 0;
static const int kSkewBaseFisheye     = 1;
static const int kSkewBaseEquisolid   = 2;

[ForceInline]
bool source_project_fov(float p[kSourceMaxParams], float3 ray, out float2 uv) {
    uv = float2(0.0f);
    if (ray.z < 1e-12f)
        return false;
    float x = ray.x / ray.z, y = ray.y / ray.z, omega = p[4];
    float r2 = x*x + y*y, om2 = omega * omega;
    const float kEps = 1e-4f;
    float factor;
    if (om2 < kEps) {
        factor = (om2 * r2) / 3.0f - om2 / 12.0f + 1.0f;
    } else if (r2 < kEps) {
        float t = tan(omega * 0.5f);
        factor = (-2.0f * t * (4.0f * r2 * t * t - 3.0f)) / (3.0f * omega);
    } else {
        float r = sqrt(r2);
        factor = atan(r * 2.0f * tan(omega * 0.5f)) / (r * omega);
    }
    uv = float2(p[0] * x * factor + p[2], p[1] * y * factor + p[3]);
    return true;
}

// Fitzgibbon division model. `simple` shares one focal length, which also
// shifts every later parameter down by one.
[ForceInline]
bool source_project_division(float p[kSourceMaxParams], bool simple,
                             float3 ray, out float2 uv) {
    uv = float2(0.0f);
    float fx = p[0], fy = simple ? p[0] : p[1];
    float cx = simple ? p[1] : p[2], cy = simple ? p[2] : p[3];
    float k  = simple ? p[3] : p[4];
    float rho2 = ray.x*ray.x + ray.y*ray.y;
    float disc = ray.z*ray.z - 4.0f * rho2 * k;
    if (disc < 0.0f)
        return false;
    float r = 2.0f / (ray.z + sqrt(disc));
    uv = float2(fx * r * ray.x + cx, fy * r * ray.y + cy);
    return true;
}

[ForceInline]
bool source_project_eucm(float p[kSourceMaxParams], float3 ray, out float2 uv) {
    uv = float2(0.0f);
    if (ray.z < 1e-12f)
        return false;
    float alpha = p[4], beta = p[5];
    float rho2 = beta * (ray.x*ray.x + ray.y*ray.y) + ray.z*ray.z;
    if (rho2 < 0.0f)
        return false;
    float den = alpha * sqrt(rho2) + (1.0f - alpha) * ray.z;
    if (den < 1e-12f)
        return false;
    uv = float2(p[0] * ray.x / den + p[2], p[1] * ray.y / den + p[3]);
    return true;
}

// Aria Fisheye624. Unlike THIN_PRISM_FISHEYE, the tangential and thin-prism
// terms act on the ALREADY radially distorted point.
[ForceInline]
bool source_project_rad_tan_thin_prism(float p[kSourceMaxParams], float3 ray,
                                       out float2 uv) {
    uv = float2(0.0f);
    if (ray.z < 1e-12f)
        return false;
    float xn = ray.x / ray.z, yn = ray.y / ray.z;
    float r = sqrt(xn*xn + yn*yn);
    float uu = xn, vv = yn;
    if (r > 1e-12f) {
        float th = atan(r);
        uu *= th / r; vv *= th / r;
    }
    float th2 = uu*uu + vv*vv, pw = 1.0f, rad = 1.0f;
    [ForceUnroll] for (int i = 0; i < 6; i++) { pw *= th2; rad += p[4 + i] * pw; }
    float x = rad * uu, y = rad * vv;
    float p0 = p[10], p1 = p[11];
    float s0 = p[12], s1 = p[13], s2 = p[14], s3 = p[15];
    float x2 = x*x, y2 = y*y, xy = x*y, r2 = x2 + y2, r4 = r2*r2;
    float xd = x + 2.0f*p1*xy + p0*(r2 + 2.0f*x2) + s0*r2 + s1*r4;
    float yd = y + 2.0f*p0*xy + p1*(r2 + 2.0f*y2) + s2*r2 + s3*r4;
    uv = float2(p[0] * xd + p[2], p[1] * yd + p[3]);
    return true;
}

// Any supported camera plus a sensor skew, which no tier carries: Agisoft
// Metashape's b2. Its b1 is an fx correction and is already folded into
// params[0]; only the off-diagonal term needs re-distorting.
//
//   p[0..3]   fx fy cx cy
//   p[4]      skew -- pixels of u per unit distorted y
//   p[5..12]  coefficients, in DistThinPrism slot order (k1 k2 k3 k4 p1 p2
//             sx1 sy1) or, when p[14] is nonzero, DistRational's (k1..k6 p1 p2)
//   p[13]     base camera, kSkewBase*
//   p[14]     0 = polynomial radial, 1 = rational radial
//
// The two radial forms and the two fisheye mappings repeat projection_utils'
// DistThinPrism / DistRational / fisheye_proj / equisolid_proj; a change there
// has to be made here too.
[ForceInline]
bool source_project_skewed(float p[kSourceMaxParams], float3 ray, out float2 uv) {
    uv = float2(0.0f);
    float2 q;
    int base = (int)p[13];
    if (base == kSkewBasePerspective) {
        if (ray.z < 1e-12f)
            return false;
        q = ray.xy / ray.z;
    } else {
        float r = length(ray.xy);
        float theta = atan2(r, ray.z);
        float k;
        if (base == kSkewBaseFisheye)
            k = theta < 1e-3f ? (1.0f - theta*theta/3.0f) / ray.z : theta / r;
        else
            k = r < 1e-6f ? (1.0f - theta*theta/24.0f) / ray.z
                          : (2.0f * sin(0.5f * theta)) / r;
        q = ray.xy * k;
    }
    float u = q.x, v = q.y, r2 = u*u + v*v;
    float dx, dy;
    if (p[14] != 0.0f) {
        float radial = (1.0f + r2*(p[5] + r2*(p[6] + r2*p[7])))
                     / (1.0f + r2*(p[8] + r2*(p[9] + r2*p[10])));
        dx = u * radial + 2.0f*p[11]*u*v + p[12]*(r2 + 2.0f*u*u);
        dy = v * radial + 2.0f*p[12]*u*v + p[11]*(r2 + 2.0f*v*v);
    } else {
        float radial = 1.0f + r2*(p[5] + r2*(p[6] + r2*(p[7] + r2*p[8])));
        dx = u * radial + 2.0f*p[9]*u*v + p[10]*(r2 + 2.0f*u*u) + p[11]*r2;
        dy = v * radial + 2.0f*p[10]*u*v + p[9]*(r2 + 2.0f*v*v) + p[12]*r2;
    }
    uv = float2(p[0] * dx + p[4] * dy + p[2], p[1] * dy + p[3]);
    return true;
}

// The model's own formula, with no check that it still describes a lens.
[ForceInline]
bool source_project_raw(int model_id, float params[kSourceMaxParams],
                        float3 ray, out float2 uv) {
    if (model_id == kSourceSkewed)
        return source_project_skewed(params, ray, uv);
    if (model_id == kSourceFOV)
        return source_project_fov(params, ray, uv);
    if (model_id == kSourceSimpleDivision)
        return source_project_division(params, true, ray, uv);
    if (model_id == kSourceDivision)
        return source_project_division(params, false, ray, uv);
    if (model_id == kSourceEUCM)
        return source_project_eucm(params, ray, uv);
    if (model_id == kSourceRadTanThinPrismFisheye)
        return source_project_rad_tan_thin_prism(params, ray, uv);
    uv = float2(0.0f);
    return false;
}

// Where the image stops growing outward the model has folded, and resampling
// past that mirrors the image instead of ending it at the lens border. Every
// model here is degree-0 homogeneous in the ray, so the test is the sign of
// the Jacobian taken in the ray's own tangent plane.
[ForceInline]
bool source_unfolded(int model_id, float params[kSourceMaxParams], float3 ray,
                     float2 uv) {
    float3 n = normalize(ray);
    float3 t = abs(n.x) < 0.9f ? float3(1.0f, 0.0f, 0.0f) : float3(0.0f, 1.0f, 0.0f);
    float3 e1 = normalize(cross(t, n));
    float3 e2 = cross(n, e1);
    const float kStep = 1e-3f;  // radians; the models are smooth well below it
    float2 a, b;
    if (!source_project_raw(model_id, params, n + kStep * e1, a)) return false;
    if (!source_project_raw(model_id, params, n + kStep * e2, b)) return false;
    float2 j1 = a - uv, j2 = b - uv;
    return j1.x * j2.y - j1.y * j2.x > 0.0f;
}

// View-space direction (need not be normalized) -> source pixel. false where
// the source camera has no image for that direction.
[CudaDeviceExport]
bool source_project(int model_id, float params[kSourceMaxParams],
                    float3 ray, out float2 uv) {
    if (!source_project_raw(model_id, params, ray, uv))
        return false;
    return source_unfolded(model_id, params, ray, uv);
}

#endif  // _CAMERA_SOURCE_SLANG_INC
