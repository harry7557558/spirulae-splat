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

flat out uvec4 vColor0;
flat out uvec4 vColor1;
flat out uvec4 vColor2;

uniform float wh_ratio;


#define USE_EXACT_DISTORTION 1

void main () {
    ivec2 pix_id = ivec2(
        3 * (index % (BBOX_OUTPUT_WIDTH/3)),
        (3*index) / BBOX_OUTPUT_WIDTH
    );
    vec4 xywh = vec4(pix_id.x, pix_id.y, 3, 1) / float(BBOX_OUTPUT_WIDTH);
    gl_Position = vec4(
        -1.0 + 2.0 * (xywh.xy + xywh.zw * (0.5+0.5*vertexPosition)) * vec2(1, wh_ratio),
        0.0, 1.0);

    vColor0 = uvec4(-1, 0, 0, 0);  // id, pos
    vColor1 = uvec4(0, 0, 0, 0);  // uv (6 x half) + rgba (4 bytes)
    vColor2 = uvec4(0, 0, 0, 0);  // xywh bounding box in tiles

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
        init_coeffs(index);
        vec3 viewdir = p_world - inverse(view)[3].xyz;
        rgb = sh_coeffs_to_color_fast(rgb, sh_degree, viewdir) + 0.5;
    }
    rgba.xyz = rgb;

    vColor0.x = uint(index);
    vColor0.y = floatBitsToUint(p_view.x);
    vColor0.z = floatBitsToUint(p_view.y);
    vColor0.w = floatBitsToUint(p_view.z);
    vColor1.x = packHalf2x16(axis_u.xy);
    vColor1.y = packHalf2x16(vec2(axis_u.z, axis_v.x));
    vColor1.z = packHalf2x16(axis_v.yz);
    vColor1.w = pack_rgba(rgba);
    vColor2.xy = uvec2(floor(max(center-bound.xy, 0.0) / float(TILE_SIZE)));
    vColor2.zw = uvec2(ceil(min(center+bound.xy, viewport-1.0) / float(TILE_SIZE)));

}
