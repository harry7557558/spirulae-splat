#version 300 es
precision highp float;
precision highp int;

#include "shader-utils.glsl"
#line 7

uniform vec2 focal;
uniform vec2 viewport;
uniform int camera_model;
uniform vec4 distortion;
uniform vec4 undistortion;

#include "shader-utils-ch.glsl"
#line 16

flat in vec4 vColor;
flat in vec3 vPosition;
flat in vec3 vAxesU;
flat in vec3 vAxesV;
flat in int vIndex;

out vec4 fragColor;


bool get_intersection(
    vec2 pos_2d,
    out vec3 poi, out vec2 uv
) {
    const float radius = 1.0;
    mat3 A = mat3(
        vAxesU, vAxesV,
        vec3(pos_2d, 1.0)
    );
    if (determinant(A) == 0.0) {
        uv = vec2(-radius);
        return false;
    }
    vec3 uvt = -inverse(A) * vPosition;
    uv = uvt.xy;
    if (length(uv) > radius)
        return false;
    float t = -uvt.z;
    poi = vec3(pos_2d*t, t);
    return t > 0.0;
}


float visibility_kernel(float r2) {
    return r2 > 1.0 ? 0.0 : 1.0-r2;
}

bool get_alpha(
    const vec2 uv,
    const float opac,
    // const vec2 aniso,
    out float alpha
) {
    float r2 = dot(uv, uv);
    float vis = visibility_kernel(r2);
    float t = 0.0; //dot(uv, aniso);
    float m = t<0. ? 1. : t>1. ? 0. :
        t*t*(2.0*t-3.0) + 1.0;
    alpha = opac * vis * m;
    return r2 >= 0. && alpha > 1e-4;
}




void main () {

    if (vPosition == vec3(0))
        discard;

    // uvec4 tex = texelFetch(u_ch_texture, ivec2(gl_FragCoord.xy), 0);
    // if (tex.x != 0u) {
    //     vec2 c1 = unpackHalf2x16(tex.x);
    //     vec2 c2 = unpackHalf2x16(tex.y);
    //     fragColor = vec4(c1, c2.x, 1.0);
    //     return;
    // }

    // fragColor = vColor * vec4(1,1,1,0.2);
    // return;

    vec2 pos2d = vec2(1,-1)*(gl_FragCoord.xy-0.5*viewport)/focal;
    vec2 pos2d_undist = pos2d;

    float fr = -1.0;
    if (camera_model == 0) {
        fr = opencv_radius(distortion);
        // pos2d_undist = undistort_opencv_iterative(pos2d, distortion);
    }
    else if (camera_model == 1) {
        fr = fisheye_radius(0.5*PI, distortion);
        pos2d_undist = undistort_fisheye(pos2d, undistortion);
    }
    const vec3 vignetting_background = vec3(0.1);
    if (fr > 0.0 && length(pos2d) > fr) {
        fragColor = vec4(vignetting_background, 0.0);
        return;
    }

    vec3 poi;
    vec2 uv;
    if (!get_intersection(pos2d_undist, poi, uv))
        discard;
    float alpha;
    if (!get_alpha(uv, vColor.a, alpha))
        discard;
    vec3 color = vColor.rgb;
    if (ch_dim > 0) {
        init_coeffs(vIndex);
        color /= 1.0 + exp(-ch_coeffs_to_color(uv));
    }
    
    if (fr > 0.0) {
        float vignetting = min(length(pos2d)/fr, 1.0);
        vignetting = pow(1.0-pow(vignetting,20.0), 2.0);
        color = mix(vignetting_background, color, vignetting);
    }

    fragColor = vec4(color, alpha);
}
