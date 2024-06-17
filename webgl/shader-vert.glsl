#version 300 es
precision highp float;
precision highp int;

uniform highp usampler2D u_texture;
uniform mat4 projection, view;
uniform vec2 focal;
uniform vec2 viewport;

in vec2 position;
in int index;

out vec4 vColor;
out vec3 vPosition;
out vec3 vAxesU;
out vec3 vAxesV;


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


void main () {
    uvec4 cen = texelFetch(u_texture,
        ivec2((uint(index) & 0x3ffu) << 1, uint(index) >> 10), 0);
    vec4 p_view = view * vec4(uintBitsToFloat(cen.xyz), 1);
    vec4 pos2d = projection * p_view;
    if (pos2d.z < 0.02) {
        gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
        return;
    }

    // patch orientation
    uvec4 info = texelFetch(u_texture,
        ivec2(((uint(index) & 0x3ffu) << 1) | 1u, uint(index) >> 10), 0);
    vec2 scale = unpackHalf2x16(info.x);
    vec4 quat = vec4(unpackHalf2x16(info.y), unpackHalf2x16(info.z));

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

    vColor = clamp(pos2d.z/pos2d.w+1.0, 0.0, 1.0) * vec4(
        (info.w) & 0xffu, (info.w >> 8) & 0xffu, (info.w >> 16) & 0xffu, (info.w >> 24) & 0xffu) / 255.0;
    vPosition = p_view.xyz;
    vAxesU = V0;
    vAxesV = V1;

    gl_Position = vec4(
        vec2(1,-1) * (-1.0 + 2.0 * (center + position * bound.xy) / viewport),
        0.0, 1.0);

}
