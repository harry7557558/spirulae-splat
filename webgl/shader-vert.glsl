#version 300 es
precision highp float;
precision highp int;

uniform highp usampler2D u_base_texture;
uniform mat4 projection, view;
uniform vec2 focal;
uniform vec2 viewport;

in vec2 position;
in int index;

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



const int sh_degree = 3;
const int sh_dim = sh_degree * (sh_degree + 2);

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
    if (degree < 1) {
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

    uvec4 info0 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u, uint(index) >> 10), 0);
    vec3 p_world = uintBitsToFloat(info0.xyz);
    vec4 p_view = view * vec4(p_world, 1);
    vec4 pos2d = projection * p_view;
    if (pos2d.z < 0.02) {
        gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
        return;
    }

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
    vec3 V0 = scale.x * R[0];
    vec3 V1 = scale.y * R[1];

    // project to 2d
    float fx = focal.x;
    float fy = focal.y;
    float cx = 0.5 * viewport.x;
    float cy = 0.5 * viewport.y;
    vec2 center;
    vec3 bound = vec3(0.0);
    const float kr = 1.0;
    bool intersect = project_ellipse_bound(
        p_view.xyz, kr*V0, kr*V1, fx, fy, cx, cy, center, bound);

    if (!intersect ||
        center.x+bound.x < 0.0 || center.x-bound.x > viewport.x ||
        center.y+bound.y < 0.0 || center.y-bound.y > viewport.y) {
        gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
        return;
    }

    vec4 rgba = vec4(unpackHalf2x16(info2.y), unpackHalf2x16(info2.z));
    vec3 rgb = rgba.xyz;
    init_coeffs();
    vec3 viewdir = p_world - inverse(view)[3].xyz;
    rgb = sh_coeffs_to_color_fast(rgb, sh_degree, viewdir) + 0.5;
    rgb = max(rgb, 0.0);

    vColor = vec4(rgb, rgba.w);
    vPosition = p_view.xyz;
    vAxesU = V0;
    vAxesV = V1;
    vAnisotropy = anisotropy;
    vIndex = index;

    gl_Position = vec4(
        vec2(1,-1) * (-1.0 + 2.0 * (center + position * bound.xy) / viewport),
        0.0, 1.0);

}
