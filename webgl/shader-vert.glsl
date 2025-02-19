#version 300 es
precision highp float;
precision highp int;

#include "shader-utils.glsl"
#line 7

uniform highp usampler2D u_base_texture;
uniform mat4 projection, view;
uniform vec2 focal;
uniform vec2 viewport;
uniform int camera_model;
uniform vec4 distortion;

#define USE_EXACT_DISTORTION 0

uniform ivec2 u_sh_config;
uniform int u_use_aniso;

layout(location = 0) in vec2 vertexPosition;
layout(location = 1) in int index;

uniform highp usampler2D u_sh_texture;

flat out vec4 vColor;
flat out vec3 vPosition;
flat out vec3 vAxesU;
flat out vec3 vAxesV;
flat out vec2 vAnisotropy;
flat out int vIndex;

#define PI 3.14159265


mat3 quat_to_rotmat(vec4 quat) {
    // quat to rotation matrix
    float w = quat.x;
    float x = quat.y;
    float y = quat.z;
    float z = quat.w;

    // glm matrices are column-major
    return mat3(
        1.f - 2.f * (y * y + z * z),  // [0][0]
        2.f * (x * y + w * z),  // [0][1]
        2.f * (x * z - w * y),  // [0][2]
        2.f * (x * y - w * z),  // [1][0]
        1.f - 2.f * (x * x + z * z),  // [1][1]
        2.f * (y * z + w * x),  // [1][2]
        2.f * (x * z + w * y),  // [2][0]
        2.f * (y * z - w * x),  // [2][1]
        1.f - 2.f * (x * x + y * y)  // [2][2]
    );
}


// axis-aligned bounding box for perspective camera
bool project_bound_perspective(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy,
    out vec2 center, out vec2 bound
) {
    // 2d conic coefficients
    // A x^2 + 2B xy + C y^2 + 2D x + 2E y + F = 0
    vec3 V01 = cross(V0, V1);
    vec3 V0T = cross(T, V0);
    vec3 V1T = cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = V0T.z * V0T.x + V1T.z * V1T.x - V01.z * V01.x;
    float E = -V01.z * V01.y + V1T.z * V1T.y + V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (B * B - A * C);
    float V = (A * E - B * D) / (B * B - A * C);
    float S = A * U*U + 2.0*B * U*V + C * V*V + 2.0*D * U + 2.0*E * V + F;

    // bounds
    float W = sqrt(C * S / (B * B - A * C));
    float H = sqrt(A * S / (B * B - A * C));
    center = vec2(fx*U+cx, fx*V+cy);
    bound = vec2(fx*W, fy*H);

    return true;
}

// axis-aligned bounding box for OpenCV camera
bool project_bound_opencv(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy, vec4 dist_coeffs,
    out vec2 center, out vec2 bound
) {
    // 2d conic coefficients
    // A x^2 + 2B xy + C y^2 + 2D x + 2E y + F = 0
    vec3 V01 = cross(V0, V1);
    vec3 V0T = cross(T, V0);
    vec3 V1T = cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = V0T.z * V0T.x + V1T.z * V1T.x - V01.z * V01.x;
    float E = -V01.z * V01.y + V1T.z * V1T.y + V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (B * B - A * C);
    float V = (A * E - B * D) / (B * B - A * C);
    float S = A * U*U + 2.0*B * U*V + C * V*V + 2.0*D * U + 2.0*E * V + F;

    // represent ellipse in principal parametric form
    vec2 p0 = vec2(U, V);
    float delta = sqrt((A-C)*(A-C) + 4.0*B*B);
    float eps = min(1e-8*S, -1e-24);
    float lambda1 = 2.0*S/min(delta-(A+C),eps);
    float lambda2 = 2.0*S/min(-delta-(A+C),eps);
    vec2 v0 = normalize(vec2(-A+C+delta, -2.0*B));
    vec2 v1 = sqrt(lambda1) * v0;
    vec2 v2 = sqrt(lambda2) * vec2(-v0.y, v0.x);

    // rejection
    float fr = opencv_radius(dist_coeffs);

    // brute force sampling on parametric
    // see here for quadratic fit code (not faster): https://github.com/harry7557558/spirulae-splat/blob/f7dc2989e5383bbeb08d5723a4801a80f11242ca/webgl/shader-vert.glsl#L238-L280
    const float N = 16.0;
    float x0 = 1e4, x1 = -1e4, y0 = 1e4, y1 = -1e4;
    for (float i = 0.0; i < N; i++) {
        float t = 2.0*PI*i/N;
        vec2 p = p0+v1*cos(t)+v2*sin(t);
        if (dot(p,p) < fr*fr) {
            vec2 c = distort_opencv(p, dist_coeffs);
            x0 = min(x0, c.x); x1 = max(x1, c.x);
            y0 = min(y0, c.y); y1 = max(y1, c.y);
        }
    }
    if (x0 > x1 || y0 > y1)
        return false;

    // update bounds
    vec2 b0 = vec2(x0, y0);
    vec2 b1 = vec2(x1, y1);
    center = 0.5*(b0+b1) * vec2(fx,fy) + vec2(cx,cy);
    bound = 0.5*(b1-b0) * vec2(fx,fy);
    return true;
}

// axis-aligned bounding box for fisheye camera, by sampling on 2D parametric
bool project_bound_fisheye_2d(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy, vec4 dist_coeffs,
    out vec2 center, out vec2 bound
) {
    // 2d conic coefficients
    // A x^2 + 2B xy + C y^2 + 2D x + 2E y + F = 0
    vec3 V01 = cross(V0, V1);
    vec3 V0T = cross(T, V0);
    vec3 V1T = cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = V0T.z * V0T.x + V1T.z * V1T.x - V01.z * V01.x;
    float E = -V01.z * V01.y + V1T.z * V1T.y + V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (B * B - A * C);
    float V = (A * E - B * D) / (B * B - A * C);
    float S = A * U*U + 2.0*B * U*V + C * V*V + 2.0*D * U + 2.0*E * V + F;

    // represent ellipse in principal parametric form
    vec2 p0 = vec2(U, V);
    float delta = sqrt((A-C)*(A-C) + 4.0*B*B);
    float eps = min(1e-8*S, -1e-24);
    float lambda1 = 2.0*S/min(delta-(A+C),eps);
    float lambda2 = 2.0*S/min(-delta-(A+C),eps);
    vec2 v0 = normalize(vec2(-A+C+delta, -2.0*B));
    vec2 v1 = sqrt(lambda1) * v0;
    vec2 v2 = sqrt(lambda2) * vec2(-v0.y, v0.x);

    // brute force sampling on parametric
    const float N = 16.0;
    float x0 = 1e4, x1 = -1e4, y0 = 1e4, y1 = -1e4;
    for (float i = 0.0; i < N; i++) {
        float t = 2.0*PI*i/N;
        vec2 c = distort_fisheye(p0+v1*cos(t)+v2*sin(t), dist_coeffs);
        x0 = min(x0, c.x); x1 = max(x1, c.x);
        y0 = min(y0, c.y); y1 = max(y1, c.y);
    }

    // update bounds
    vec2 b0 = vec2(x0, y0);
    vec2 b1 = vec2(x1, y1);
    center = 0.5*(b0+b1) * vec2(fx,fy) + vec2(cx,cy);
    bound = 0.5*(b1-b0) * vec2(fx,fy);

    return true;
}

// axis-aligned bounding box for fisheye camera, by sampling on 3D parametric
bool project_bound_fisheye(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy, vec4 dist_coeffs,
    out vec2 center, out vec2 bound
) {
    // sampling
    const float N = 16.0;
    float x0 = 1e4, x1 = -1e4, y0 = 1e4, y1 = -1e4;
    for (float i = 0.0; i < N; i++) {
        float t = 2.0*PI*i/N;
        vec3 p3d = T+V0*cos(t)+V1*sin(t);
        if (p3d.z > 1e-4) {
            vec2 p = p3d.xy / p3d.z;
            vec2 c = distort_fisheye(p, dist_coeffs);
            x0 = min(x0, c.x); x1 = max(x1, c.x);
            y0 = min(y0, c.y); y1 = max(y1, c.y);
        }
    }
    if (x0 > x1 || y0 > y1)
        return false;

    // update bounds
    vec2 b0 = vec2(x0, y0);
    vec2 b1 = vec2(x1, y1);
    center = 0.5*(b0+b1) * vec2(fx,fy) + vec2(cx,cy);
    bound = 0.5*(b1-b0) * vec2(fx,fy);
    return true;
}



#define use_sh bool(u_sh_config.x)
#define sh_degree u_sh_config.y
#define sh_dim (sh_degree * (sh_degree + 2))

ivec3 pix_sh;
uvec4 coeff_sh;

void init_coeffs() {
    int width = textureSize(u_sh_texture, 0).x;
    int pixidx = index * ((sh_dim+1)/2);
    pix_sh = ivec3(pixidx % width, pixidx / width, 0);
}

vec3 next_coeff() {
    uvec2 packed;
    if (pix_sh.z == 0) {
        coeff_sh = texelFetch(u_sh_texture, pix_sh.xy, 0);
        packed = coeff_sh.xy;
        pix_sh.z += 1;
    }
    else {
        packed = coeff_sh.zw;
        pix_sh.x += 1;
        pix_sh.z = 0;
    }
    vec2 c1 = unpackHalf2x16(packed.x);
    vec2 c2 = unpackHalf2x16(packed.y);
    return vec3(c1, c2.x);
}

vec3 sh_coeffs_to_color_fast(
    vec3 base_color,
    int degree,
    vec3 viewdir
) {
    vec3 color = 0.2820947917738781 * base_color;
    if (degree < 1 || !use_sh) {
        return color;
    }

    float norm = length(viewdir);
    float x = viewdir.x / norm;
    float y = viewdir.y / norm;
    float z = viewdir.z / norm;

    float fTmp0A = 0.48860251190292;
    color += fTmp0A * (-y * next_coeff());
    color += fTmp0A * (z * next_coeff());
    color += fTmp0A * (-x * next_coeff());
    if (degree < 2) {
        return color;
    }
    float z2 = z * z;

    float fTmp0B = -1.092548430592079 * z;
    float fTmp1A = 0.5462742152960395;
    float fC1 = x * x - y * y;
    float fS1 = 2. * x * y;
    float pSH6 = (0.9461746957575601 * z2 - 0.3153915652525201);
    float pSH7 = fTmp0B * x;
    float pSH5 = fTmp0B * y;
    float pSH8 = fTmp1A * fC1;
    float pSH4 = fTmp1A * fS1;
    color += pSH4 * next_coeff();
    color += pSH5 * next_coeff();
    color += pSH6 * next_coeff();
    color += pSH7 * next_coeff();
    color += pSH8 * next_coeff();
    if (degree < 3) {
        return color;
    }

    float fTmp0C = -2.285228997322329 * z2 + 0.4570457994644658;
    float fTmp1B = 1.445305721320277 * z;
    float fTmp2A = -0.5900435899266435;
    float fC2 = x * fC1 - y * fS1;
    float fS2 = x * fS1 + y * fC1;
    float pSH12 = z * (1.865881662950577 * z2 - 1.119528997770346);
    float pSH13 = fTmp0C * x;
    float pSH11 = fTmp0C * y;
    float pSH14 = fTmp1B * fC1;
    float pSH10 = fTmp1B * fS1;
    float pSH15 = fTmp2A * fC2;
    float pSH9  = fTmp2A * fS2;
    color += pSH9 * next_coeff();
    color += pSH10 * next_coeff();
    color += pSH11 * next_coeff();
    color += pSH12 * next_coeff();
    color += pSH13 * next_coeff();
    color += pSH14 * next_coeff();
    color += pSH15 * next_coeff();
    if (degree < 4) {
        return color;
    }

    float fTmp0D = z * (-4.683325804901025 * z2 + 2.007139630671868);
    float fTmp1C = 3.31161143515146 * z2 - 0.47308734787878;
    float fTmp2B = -1.770130769779931 * z;
    float fTmp3A = 0.6258357354491763;
    float fC3 = x * fC2 - y * fS2;
    float fS3 = x * fS2 + y * fC2;
    float pSH20 = (1.984313483298443 * z * pSH12 - 1.006230589874905 * pSH6);
    float pSH21 = fTmp0D * x;
    float pSH19 = fTmp0D * y;
    float pSH22 = fTmp1C * fC1;
    float pSH18 = fTmp1C * fS1;
    float pSH23 = fTmp2B * fC2;
    float pSH17 = fTmp2B * fS2;
    float pSH24 = fTmp3A * fC3;
    float pSH16 = fTmp3A * fS3;
    color += pSH16 * next_coeff();
    color += pSH17 * next_coeff();
    color += pSH18 * next_coeff();
    color += pSH19 * next_coeff();
    color += pSH20 * next_coeff();
    color += pSH21 * next_coeff();
    color += pSH22 * next_coeff();
    color += pSH23 * next_coeff();
    color += pSH24 * next_coeff();

    return color;
}


void main () {
    gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
    vPosition = vec3(0);

    if (index == -1) return;

    uvec4 info0 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u, uint(index) >> 10), 0);
    vec3 p_world = uintBitsToFloat(info0.xyz);
    vec4 p_view = view * vec4(p_world, 1);
    vec4 pos2d = projection * p_view;
    if (pos2d.w <= 0.0 || p_view.z <= 0.0)
        return;

    uvec4 info1 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u + 1u, uint(index) >> 10), 0);
    uvec4 info2 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u + 2u, uint(index) >> 10), 0);

    // patch orientation
    vec2 scale = uintBitsToFloat(uvec2(info0.w, info1.x));
    vec4 quat = vec4(unpackHalf2x16(info1.y), unpackHalf2x16(info1.z));
    vec2 anisotropy = uintBitsToFloat(uvec2(info1.w, info2.x));

    mat3 R0 = mat3(view);
    mat3 Rq = quat_to_rotmat(quat);
    mat3 R = R0 * Rq;
    vec3 axis_u = scale.x * R[0];
    vec3 axis_v = scale.y * R[1];

    // distortion
#if USE_EXACT_DISTORTION
    ;
#else  // USE_EXACT_DISTORTION
    vec2 p_proj = p_view.xy/p_view.z;
    vec2 p_dist = p_proj;
    mat2 jac_dist;
    if (camera_model == 0) {  // OPENCV
        float fr = opencv_radius(distortion);
        if (fr > 0.0 && !(dot(p_proj, p_proj) < fr*fr))
            return;
        distort_opencv_with_jac(p_proj, distortion, p_dist, jac_dist);
        axis_u = vec3(jac_dist * vec2(axis_u), axis_u.z);
        axis_v = vec3(jac_dist * vec2(axis_v), axis_v.z);
    }
    else if (camera_model == 1) {  // OPENCV_FISHEYE
        distort_fisheye_with_jac(p_proj, distortion, p_dist, jac_dist);
        float fr = fisheye_radius(0.5*PI, distortion);
        if (!(dot(p_dist,p_dist) < fr*fr))
            return;
        axis_u = distort_jac_3d(p_proj, p_dist, jac_dist, axis_u);
        axis_v = distort_jac_3d(p_proj, p_dist, jac_dist, axis_v);
    }
    p_view.xy = p_dist * p_view.z;
#endif  // USE_EXACT_DISTORTION

    // project to 2d
    float fx = focal.x;
    float fy = focal.y;
    float cx = 0.5 * viewport.x;
    float cy = 0.5 * viewport.y;
    vec2 center;
    vec2 bound = vec2(0.0);
    const float kr = 1.0;

    // compute axis-aligned bounding box
#if USE_EXACT_DISTORTION
    bool intersect;
    if (camera_model == 0)
        intersect = project_bound_opencv(
            p_view.xyz, kr*axis_u, kr*axis_v,
            fx, fy, cx, cy, distortion,
            center, bound
        );
    else if (camera_model == 1)
        intersect = project_bound_fisheye(
            p_view.xyz, kr*axis_u, kr*axis_v,
            fx, fy, cx, cy, distortion,
            center, bound
        );
    else
        intersect = project_bound_perspective(
            p_view.xyz, kr*axis_u, kr*axis_v,
            fx, fy, cx, cy,
            center, bound
        );
#else
    bool intersect = project_bound_perspective(
        p_view.xyz, kr*axis_u, kr*axis_v,
        fx, fy, cx, cy,
        center, bound
    );
#endif

    if (!intersect ||
        center.x+bound.x < 0.0 || center.x-bound.x > viewport.x ||
        center.y+bound.y < 0.0 || center.y-bound.y > viewport.y
    ) return;

    vec4 rgba = vec4(unpackHalf2x16(info2.y), unpackHalf2x16(info2.z));
    vec3 rgb = rgba.xyz;
    if (sh_dim > 0) {
        init_coeffs();
        vec3 viewdir = p_world - inverse(view)[3].xyz;
        rgb = sh_coeffs_to_color_fast(rgb, sh_degree, viewdir) + 0.5;
    }
    rgb = max(rgb, 0.0);

    vColor = vec4(rgb, rgba.w);
    vPosition = p_view.xyz;
    vAxesU = axis_u;
    vAxesV = axis_v;
    vAnisotropy = u_use_aniso == 0 ? vec2(0) : anisotropy;
    vIndex = index;

    gl_Position = vec4(
        vec2(1,-1) * (-1.0 + 2.0 * (center + vertexPosition * bound.xy) / viewport),
        0.0, 1.0);

}
