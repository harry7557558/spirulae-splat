#version 300 es
precision highp float;
precision highp int;

uniform highp usampler2D u_base_texture;
uniform mat4 projection, view;
uniform vec2 focal;
uniform vec2 viewport;
uniform int camera_model;
uniform vec4 distortion;

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


// camera distortion models

float opencv_radius(
    in vec4 dist_coeffs
) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;

    // calculate the extrema of dist(r)
    if (k2 == 0.0) {
        float b = -1.0/(3.0*k1);
        return b > 0.0 ? sqrt(b) : -1.0;
    }
    float disc = 9.0*k1*k1-20.0*k2;
    if (disc <= 0.0) return -1.0;
    disc = sqrt(disc);
    float u1 = (-3.0*k1 + disc) / (10.0*k2);
    float u2 = (-3.0*k1 - disc) / (10.0*k2);
    if (u1 <= 0.0) return u2>0.0 ? sqrt(u2) : -1.0;
    return u2>0.0 ? sqrt(min(u1,u2)) : sqrt(u1);
}

float fisheye_radius(
    in float theta, in vec4 dist_coeffs
) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float theta2 = theta*theta;
    return theta*(1.0 + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));
}

void distort_opencv(
    in vec2 p, in vec4 dist_coeffs,
    out vec2 p_dist, out mat2 jac_dist
) {
    // https://www.desmos.com/calculator/adnh2fzdyi

    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float p1 = dist_coeffs.z;
    float p2 = dist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r_dist = 1.0 + k1*r2 + k2*r2*r2;

    float x_dist = x*r_dist + 2.0*p1*x*y + p2*(r2 + 2.0*x*x);
    float y_dist = y*r_dist + p1*(r2 + 2.0*y*y) + 2.0*p2*x*y;
    p_dist = vec2(x_dist, y_dist);

    float jxx = r_dist + 2.0*( x*x*(k1 + 2.0*k2*r2) + p1*y + 3.0*p2*x );
    float jxy = 2.0*( x*y*(k1 + 2.0*k2*r2) + p1*x + p2*y );
    float jyy = r_dist + 2.0*( y*y*(k1 + 2.0*k2*r2) + 3.0*p1*y + p2*x );
    jac_dist = mat2(jxx, jxy, jxy, jyy);
}

void distort_fisheye(
    in vec2 p, in vec4 dist_coeffs,
    out vec2 p_dist, out mat2 jac_dist
) {
    // https://www.desmos.com/calculator/yg7gqrtgm9

    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r = sqrt(r2);
    float theta = atan(r);
    float theta2 = theta*theta;
    float r_dist = theta*(1.0 + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));
    float s_dist = 1.0 + theta2*(3.0*k1 + theta2*(5.0*k2 + theta2*(7.0*k3 + theta2*9.0*k4)));

    float x_dist = x/r * r_dist;
    float y_dist = y/r * r_dist;
    p_dist = vec2(x_dist, y_dist);

    float jd1 = s_dist / (r2*(r2+1.0));
    float jd2 = r_dist / (r2*r);
    float jxx = jd1 * x*x + jd2 * y*y;
    float jxy = jd1 * x*y - jd2 * x*y;
    float jyy = jd1 * y*y + jd2 * x*x;
    jac_dist = mat2(jxx, jxy, jxy, jyy);
}

vec3 distort_jac_3d(
    in vec2 p2d_u, in vec2 p2d_d,
    in mat2 j2d, in vec3 dp
) {
    // transform dp using the Jacobian of (x,y,z) -> (distort(x/z,y/z)*z,z)
    mat3x2 j3d = j2d * mat3x2(1, 0, 0, 1, -p2d_u.x, -p2d_u.y)
        + mat3x2(0, 0, 0, 0, p2d_d.x, p2d_d.y);
    return vec3(j3d * dp, dp.z);
}


bool project_ellipse_bound(
    vec3 T,
    vec3 V0,
    vec3 V1,
    float fx, float fy, float cx, float cy,
    out vec2 center, out vec3 bound
) {
    // 2d conic coefficients
    vec3 V01 = cross(V0, V1);
    vec3 V0T = cross(T, V0);
    vec3 V1T = cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = 2.0 * V0T.z * V0T.x + 2.0 * V1T.z * V1T.x - 2.0 * V01.z * V01.x;
    float E = -2.0 * V01.z * V01.y + 2.0 * V1T.z * V1T.y + 2.0 * V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (2.0 * (B * B - A * C));
    float V = (A * E - B * D) / (2.0 * (B * B - A * C));
    float S = -(A * U * U + 2.0 * B * U * V + C * V * V + D * U + E * V + F);

    // image transform
    float U_T = fx * U + cx;
    float V_T = fy * V + cy;
    float A_T = A / (fx * fx);
    float B_T = B / (fx * fy);
    float C_T = C / (fy * fy);

    // axis-aligned bounding box
    float W_T = fx * sqrt(C * S / (A * C - B * B));
    float H_T = fy * sqrt(A * S / (A * C - B * B));

    // bounding circle
    float L_T = 0.5 * (A_T + C_T - sqrt((A_T - C_T) * (A_T - C_T) + 4.0 * B_T * B_T));
    float R_T = sqrt(S / L_T);

    // output
    center = vec2(U_T, V_T);
    bound = vec3(W_T, H_T, R_T);
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

    if (index == -1) return;

    uvec4 info0 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u, uint(index) >> 10), 0);
    vec3 p_world = uintBitsToFloat(info0.xyz);
    vec4 p_view = view * vec4(p_world, 1);
    vec4 pos2d = projection * p_view;
    if (pos2d.w < 0.0 || p_view.z < 0.01)
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
    vec2 p_proj = p_view.xy/p_view.z;
    vec2 p_dist = p_proj;
    mat2 jac_dist;
    if (camera_model == 0) {  // OPENCV
        float fr = opencv_radius(distortion);
        if (fr > 0.0 && !(dot(p_proj, p_proj) < fr*fr))
            return;
        distort_opencv(p_proj, distortion, p_dist, jac_dist);
        axis_u = vec3(jac_dist * vec2(axis_u), axis_u.z);
        axis_v = vec3(jac_dist * vec2(axis_v), axis_v.z);
    }
    else if (camera_model == 1) {  // OPENCV_FISHEYE
        distort_fisheye(p_proj, distortion, p_dist, jac_dist);
        float fr = fisheye_radius(1.570796, distortion);
        if (!(dot(p_dist,p_dist) < fr*fr))
            return;
        axis_u = distort_jac_3d(p_proj, p_dist, jac_dist, axis_u);
        axis_v = distort_jac_3d(p_proj, p_dist, jac_dist, axis_v);
        // eliminate splats outside the circle
        // TODO: use an analytical solution
        #if 0
        vec2 u2d = vec2(axis_u)/p_view.z;
        vec2 v2d = vec2(axis_v)/p_view.z;
        const int n = 8;
        for (int i = 0; i < n; i++) {
            float t = 6.283185*float(i)/float(n);
            vec2 p = p_dist + u2d * cos(t) + v2d * sin(t);
            if (!(dot(p,p) < fr*fr)) {
                gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
                return;
            }
        }
        #endif
    }
    p_view.xy = p_dist * p_view.z;

    // project to 2d
    float fx = focal.x;
    float fy = focal.y;
    float cx = 0.5 * viewport.x;
    float cy = 0.5 * viewport.y;
    vec2 center;
    vec3 bound = vec3(0.0);
    const float kr = 1.0;
    bool intersect = project_ellipse_bound(
        p_view.xyz, kr*axis_u, kr*axis_v, fx, fy, cx, cy, center, bound);

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
