#version 300 es
precision highp float;
precision highp int;

#include "shader-utils.glsl"
#line 7

uniform int GS_DIM;

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

flat out uvec4 vColor0;  // position
flat out uvec4 vColor1;  // shape
flat out vec4 vColor2;  // rgba


void main () {
    gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
    vColor0.x = uint(-1);

    if (index == -1) return;

    uvec4 info0 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u, uint(index) >> 10), 0);
    vec3 p_world = uintBitsToFloat(info0.xyz);
    vec4 p_view = view * vec4(p_world, 1);
    vec4 pos2d = projection * p_view;
    if (camera_model != 1 && (pos2d.w <= 0.0 || p_view.z <= 0.0))
        return;

    uvec4 info1 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u + 1u, uint(index) >> 10), 0);
    uvec4 info2 = texelFetch(u_base_texture,
        ivec2((uint(index) & 0x3ffu) * 3u + 2u, uint(index) >> 10), 0);

    // patch orientation
    vec3 scale = vec3(0);
    if (GS_DIM == 2) {
        scale.xy = uintBitsToFloat(uvec2(info0.w, info1.x));
    }
    else if (GS_DIM == 3) {
        scale.xy = unpackHalf2x16(info0.w);
        scale.z = uintBitsToFloat(info1.x);
    }
    vec4 quat = vec4(unpackHalf2x16(info1.y), unpackHalf2x16(info1.z));

    mat3 R0 = mat3(view);
    mat3 Rq = quat_to_rotmat(quat);
    mat3 R = R0 * Rq;
    vec3 axis_u = scale.x * R[0];
    vec3 axis_v = scale.y * R[1];
    vec3 axis_w = scale.z * R[2];

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
            float fr = fisheye_radius(PI, distortion);
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
    mat2 bound = mat2(0.0);  // obb
    
    vec3 V0 = axis_u;
    vec3 V1 = axis_v;

    vec4 rgba = vec4(unpackHalf2x16(info2.y), unpackHalf2x16(info2.z));
    if (GS_DIM == 3) {
        const float min_opac = 0.01;
        float extend = sqrt(2.0*max(log(rgba.w/min_opac), 0.0));
        project_ellipsoid_to_ellipse(
            axis_u, axis_v, axis_w,
            normalize(p_view.xyz),
            V0, V1
        );
        V0 *= extend;
        V1 *= extend;
    }

    // compute axis-aligned bounding box
    bool intersect;
    if (useExactDistortion) {
        vec2 bound2;
        if (camera_model == 0)
            intersect = project_bound_opencv(
                p_view.xyz, V0, V1,
                fx, fy, cx, cy, distortion,
                center, bound2  // TODO
            );
        else if (camera_model == 1)
            intersect = project_bound_fisheye(
                p_view.xyz, V0, V1,
                fx, fy, cx, cy, distortion,
                center, bound
            );
        else
            intersect = project_bound_perspective(
                p_view.xyz, V0, V1,
                fx, fy, cx, cy,
                center, bound
            );
    }
    else {
        intersect = project_bound_perspective(
            p_view.xyz, V0, V1,
            fx, fy, cx, cy,
            center, bound
        );
    }

    vec2 bound_extend = abs(bound[0]) + abs(bound[1]);
    if (!intersect ||
        center.x+bound_extend.x < 0.0 || center.x-bound_extend.x > viewport.x ||
        center.y+bound_extend.y < 0.0 || center.y-bound_extend.y > viewport.y ||
        bound_extend.x > 0.5*viewport.x || bound_extend.y > 0.5*viewport.y
    ) return;

    vec3 rgb = rgba.xyz;
    if (sh_dim > 0) {
        init_coeffs(index);
        vec3 viewdir = p_world - inverse(view)[3].xyz;
        rgb = sh_coeffs_to_color_fast(rgb, sh_degree, viewdir) + 0.5;
    }
    rgb = max(rgb, 0.0);

    vColor0.x = uint(index);
    vColor0.y = floatBitsToUint(p_view.x);
    vColor0.z = floatBitsToUint(p_view.y);
    vColor0.w = floatBitsToUint(p_view.z);
    vColor1.x = packHalf2x16(axis_u.xy);
    vColor1.y = packHalf2x16(vec2(axis_u.z, axis_v.x));
    vColor1.z = packHalf2x16(axis_v.yz);
    vColor1.w = floatBitsToUint(scale.z/(scale.x*scale.y));
    vColor2 = vec4(rgb, rgba.w);

    gl_Position = vec4(
        vec2(1,-1) * (-1.0 + 2.0 * (center + bound * vertexPosition) / viewport),
        0.0, 1.0);

}
