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

layout(location = 0) in vec2 vertexPosition;
layout(location = 1) in int index;

#include "shader-utils-sh.glsl"
#line 20

flat out vec4 vColor;
flat out vec3 vPosition;
flat out vec3 vAxesU;
flat out vec3 vAxesV;
flat out int vIndex;


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

    mat3 R0 = mat3(view);
    mat3 Rq = quat_to_rotmat(quat);
    mat3 R = R0 * Rq;
    vec3 axis_u = scale.x * R[0];
    vec3 axis_v = scale.y * R[1];

    // distortion
    bool useExactDistortion = (camera_model == 1);
    if (!useExactDistortion) {
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
    }

    // project to 2d
    float fx = focal.x;
    float fy = focal.y;
    float cx = 0.5 * viewport.x;
    float cy = 0.5 * viewport.y;
    vec2 center;
    vec2 bound = vec2(0.0);
    const float kr = 1.0;

    // compute axis-aligned bounding box
    bool intersect;
    if (useExactDistortion) {
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
    }
    else {
        intersect = project_bound_perspective(
            p_view.xyz, kr*axis_u, kr*axis_v,
            fx, fy, cx, cy,
            center, bound
        );
    }

    if (!intersect ||
        center.x+bound.x < 0.0 || center.x-bound.x > viewport.x ||
        center.y+bound.y < 0.0 || center.y-bound.y > viewport.y
    ) return;

    vec4 rgba = vec4(unpackHalf2x16(info2.y), unpackHalf2x16(info2.z));
    vec3 rgb = rgba.xyz;
    if (sh_dim > 0) {
        init_coeffs(index);
        vec3 viewdir = p_world - inverse(view)[3].xyz;
        rgb = sh_coeffs_to_color_fast(rgb, sh_degree, viewdir) + 0.5;
    }
    rgb = max(rgb, 0.0);

    vColor = vec4(rgb, rgba.w);
    vPosition = p_view.xyz;
    vAxesU = axis_u;
    vAxesV = axis_v;
    vIndex = index;

    gl_Position = vec4(
        vec2(1,-1) * (-1.0 + 2.0 * (center + vertexPosition * bound.xy) / viewport),
        0.0, 1.0);

}
