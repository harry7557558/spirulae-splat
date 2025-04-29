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

uniform highp usampler2D u_psa_texture;
uniform highp usampler2D u_int_texture;

uniform highp usampler2D u_proj_texture;


out vec4 fragColor;


bool get_intersection(
    vec3 pos, vec3 axis_u, vec3 axis_v,
    vec2 pos_2d,
    out vec3 poi, out vec2 uv
) {
    const float radius = 1.0;
    mat3 A = mat3(
        axis_u, axis_v,
        vec3(pos_2d, 1.0)
    );
    if (determinant(A) == 0.0) {
        uv = vec2(-radius);
        return false;
    }
    vec3 uvt = -inverse(A) * pos;
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


#define MAX_N_INT 16
int intersectCount = 0;
uvec2 intersects[MAX_N_INT];  // id, packed rgba

void main () {

    vec2 fragCoord = gl_FragCoord.xy;
    fragCoord.y = viewport.y - fragCoord.y - 1.0;
    ivec2 tile_id = ivec2(fragCoord) / TILE_SIZE;
    ivec2 bounds = ivec2(texelFetch(u_psa_texture, tile_id, 0).xy);
    // fragColor = vec4(1e-3*vec3(float(bounds.y)-float(bounds.x)), 1);

    vec2 pos2d = (fragCoord-0.5*viewport)/focal;
    vec2 pos2d_undist = pos2d;

    float fr = -1.0;
    if (camera_model == 0) {
        fr = opencv_radius(distortion);
        pos2d_undist = undistort_opencv_iterative(pos2d, distortion);
    }
    else if (camera_model == 1) {
        fr = fisheye_radius(0.5*PI, distortion);
        pos2d_undist = undistort_fisheye(pos2d, undistortion);
    }
    const vec3 vignetting_background = vec3(0.1);
    if (fr > 0.0 && length(pos2d) > fr) {
        fragColor = vec4(vignetting_background, 1.0);
        return;
    }

    // normal rasterization

	float T = 1.0;

    for (int i = bounds.y-1; i >= bounds.x; i--) {
        int index = int(texelFetch(u_int_texture, ivec2(i%2048, i/2048), 0).x);

        ivec2 pix_id = ivec2(
            3 * (index % (BBOX_OUTPUT_WIDTH/3)),
            (3*index) / BBOX_OUTPUT_WIDTH
        );
        uvec4 comp0 = texelFetch(u_proj_texture, pix_id+ivec2(0,0), 0);
        uvec4 comp1 = texelFetch(u_proj_texture, pix_id+ivec2(1,0), 0);

        vec3 p_view = vec3(uintBitsToFloat(comp0.y), uintBitsToFloat(comp0.z), uintBitsToFloat(comp0.w));
        mat2x3 axes = mat2x3(
            unpackHalf2x16(comp1.x),
            unpackHalf2x16(comp1.y),
            unpackHalf2x16(comp1.z)
        );
        vec4 rgba = unpack_rgba(comp1.w);

        vec3 poi;
        vec2 uv;
        if (!get_intersection(p_view, axes[0], axes[1], pos2d_undist, poi, uv))
            continue;
        float alpha;
        if (!get_alpha(uv, rgba.w, alpha))
            continue;
        vec3 color = rgba.xyz;
        // if (ch_dim > 0) {
        //     init_coeffs(vIndex);
        //     color /= 1.0 + exp(-ch_coeffs_to_color(uv));
        // }

        if (fr > 0.0) {
            float vignetting = min(length(pos2d)/fr, 1.0);
            vignetting = pow(1.0-pow(vignetting,20.0), 2.0);
            color = mix(vignetting_background, color, vignetting);
        }

		T *= 1.0 - alpha;

        intersects[intersectCount++] = uvec2(
            uint(16777216.0*poi.z),
            pack_rgba(vec4(color, alpha))
        );
        if (intersectCount >= MAX_N_INT)
            break;

		if (T < 0.01)
			break;
    }

    // insertion sort
    for (int i = 1; i < intersectCount; ++i) {
        uvec2 z = intersects[i];
        int j = i-1;
        while (j >= 0 && intersects[j].x > z.x) {
            intersects[j+1] = intersects[j];
            j--;
        }
        intersects[j+1] = z;
    }

    // cumulate color
	vec3 cumColor = vec3(0);
    T = 1.0;
    for (int i = 0; i < intersectCount; ++i) {
        vec4 rgba = unpack_rgba(intersects[i].y);
        cumColor += T * rgba.w * rgba.xyz;
		T *= 1.0 - rgba.w;
    }

    fragColor = vec4(cumColor, 1.0-T);

}
