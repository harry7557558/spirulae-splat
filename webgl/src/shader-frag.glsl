#version 300 es
precision highp float;
precision highp int;

#include "shader-utils.glsl"
#line 7

uniform int GS_DIM;

uniform vec2 focal;
uniform vec2 viewport;
uniform int camera_model;
uniform vec4 distortion;
uniform vec4 undistortion;

#include "shader-utils-ch.glsl"
#line 18

flat in uvec4 vColor0;  // position
flat in uvec4 vColor1;  // shape
flat in vec4 vColor2;  // rgba

out vec4 fragColor;


bool get_intersection(
    vec3 pos, vec3 axis_u, vec3 axis_v,
    vec3 raydir,
    out vec3 poi, out vec2 uv
) {
    const float radius = 1.0;
    mat3 A = mat3(axis_u, axis_v, raydir);
    vec3 uvt = -inverse(A) * pos;
    uv = uvt.xy;
    float t = -uvt.z;
    poi = raydir * t;
    return length(uv) < radius && t > 0.0;
}


float visibility_kernel_2dgs(float r2) {
    return r2 > 1.0 ? 0.0 : 1.0-r2;
}

bool get_alpha_2dgs(
    const vec2 uv,
    const float opac,
    // const vec2 aniso,
    out float alpha
) {
    float r2 = dot(uv, uv);
    float vis = visibility_kernel_2dgs(r2);
    float t = 0.0; //dot(uv, aniso);
    float m = t<0. ? 1. : t>1. ? 0. :
        t*t*(2.0*t-3.0) + 1.0;
    alpha = opac * vis * m;
    return r2 >= 0. && alpha > 1e-4;
}




void main () {
    int index = int(vColor0.x);
    if (index < 0)
        discard;

    vec4 rgba = vColor2;
    // fragColor = rgba; return;

    vec2 pos2d = vec2(1,-1)*(gl_FragCoord.xy-0.5*viewport)/focal;
    const vec3 vignetting_background = vec3(0.1);
    float fr = -1.0;
    if (camera_model == 0)
        fr = opencv_radius(distortion);
    else if (camera_model == 1)
        fr = fisheye_radius(PI, distortion);
    if (fr > 0.0 && length(pos2d) > fr) {
        fragColor = vec4(vignetting_background, 0.0);
        return;
    }

    vec3 color = rgba.rgb;
    float alpha = 0.0;

    // custom 2DGS
    if (GS_DIM == 2) {
        vec3 p_view = vec3(
            uintBitsToFloat(vColor0.y),
            uintBitsToFloat(vColor0.z),
            uintBitsToFloat(vColor0.w)
        );
        mat2x3 axes = mat2x3(
            unpackHalf2x16(vColor1.x),
            unpackHalf2x16(vColor1.y),
            unpackHalf2x16(vColor1.z)
        );

        vec3 viewdir;
        if (camera_model == 0)
            viewdir = vec3(pos2d, 1.0);
        else if (camera_model == 1)
            viewdir = unproject_fisheye(pos2d, undistortion);
        else
            viewdir = vec3(pos2d, 1.0);

        vec3 poi; vec2 uv;
        if (!get_intersection(p_view, axes[0], axes[1], viewdir, poi, uv))
            discard;
        if (!get_alpha_2dgs(uv, rgba.a, alpha))
            discard;
        if (ch_dim > 0) {
            init_coeffs(index);
            color /= 1.0 + exp(-ch_coeffs_to_color(uv));
        }
    }

    // 3DGS
    if (GS_DIM == 3) {
        vec2 t = vec2(
            uintBitsToFloat(vColor0.y),
            uintBitsToFloat(vColor0.z)
        );
        vec3 cov = vec3(
            uintBitsToFloat(vColor1.x),
            uintBitsToFloat(vColor1.y),
            uintBitsToFloat(vColor1.z)
        );
        float dx = t.x - pos2d.x;
        float dy = t.y - pos2d.y;
        float sigma = 0.5 * (cov.x * dx*dx + cov.z * dy*dy) + cov.y * dx*dy;
        alpha = rgba.a * exp(-sigma);
        if (sigma < 0.0 || alpha < 0.01)
            discard;
    }

    if (fr > 0.0) {
        float vignetting = min(length(pos2d)/fr, 1.0);
        vignetting = pow(1.0-pow(vignetting,20.0), 2.0);
        color = mix(vignetting_background, color, vignetting);
    }

    fragColor = vec4(color, alpha);
}
