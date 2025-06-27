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
    vec3 raydir,
    out vec3 poi, out vec2 uv
) {
    const float radius = 1.0;
    mat3 A = mat3(axis_u, axis_v, raydir);
    if (determinant(A) == 0.0) {
        uv = vec2(-radius);
        return false;
    }
    vec3 uvt = -inverse(A) * pos;
    uv = uvt.xy;
    if (length(uv) > radius)
        return false;
    float t = -uvt.z;
    poi = raydir * t;
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

void pps_compare_swap(int i, int j) {
    uvec2 u = intersects[i], v = intersects[j];
    if (u.x > v.x)
        intersects[j] = u, intersects[i] = v;
}

void main () {

    vec2 fragCoord = gl_FragCoord.xy;
    fragCoord.y = viewport.y - fragCoord.y - 1.0;
    ivec2 tile_id = ivec2(fragCoord) / TILE_SIZE;
    ivec2 bounds = ivec2(texelFetch(u_psa_texture, tile_id, 0).xy);
    // fragColor = vec4(1e-3*vec3(float(bounds.y)-float(bounds.x)), 1); return;

    vec2 pos2d = (fragCoord-0.5*viewport)/focal;
    vec3 viewdir;

    float fr = -1.0;
    if (camera_model == 0) {
        fr = opencv_radius(distortion);
        vec2 pos2d_undist = undistort_opencv_iterative(pos2d, distortion);
        viewdir = vec3(pos2d_undist, 1.0);
    }
    else if (camera_model == 1) {
        fr = fisheye_radius(PI, distortion);
        // pos2d_undist = undistort_fisheye(pos2d, undistortion);
        // viewdir = vec3(pos2d_undist, 1.0);
        viewdir = unproject_fisheye(pos2d, undistortion);
    }
    else {
        viewdir = vec3(pos2d, 1.0);
    }
    const vec3 vignetting_background = vec3(0.1);
    if (fr > 0.0 && length(pos2d) > fr) {
        fragColor = vec4(vignetting_background, 1.0);
        return;
    }

    // normal rasterization

	float T = 1.0;

    for (int i = bounds.y-1; i >= bounds.x; i--) {
        uvec4 comp0 = texelFetch(u_int_texture, ivec2((2*i+0)%2048, (2*i+0)/2048), 0);
        uvec4 comp1 = texelFetch(u_int_texture, ivec2((2*i+1)%2048, (2*i+1)/2048), 0);

        vec3 p_view = vec3(uintBitsToFloat(comp0.y), uintBitsToFloat(comp0.z), uintBitsToFloat(comp0.w));
        mat2x3 axes = mat2x3(
            unpackHalf2x16(comp1.x),
            unpackHalf2x16(comp1.y),
            unpackHalf2x16(comp1.z)
        );
        vec4 rgba = unpack_rgba(comp1.w);

        vec3 poi;
        vec2 uv;
        if (!get_intersection(p_view, axes[0], axes[1], viewdir, poi, uv))
            continue;
        float alpha;
        if (!get_alpha(uv, rgba.w, alpha))
            continue;
        vec3 color = rgba.xyz;
        if (ch_dim > 0) {
            init_coeffs(int(comp0.x));
            color /= 1.0 + exp(-ch_coeffs_to_color(uv));
        }

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

#if 0
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
#else
    // https://bertdobbelaere.github.io/sorting_networks.html
    // s = """<paste sorting network>"""
    // for i, j in sum([eval(x) for x in s.strip().split('\n')], []):
    //     print(f"_({i},{j})", end='')
    // print(' ')
    #define _(i,j) pps_compare_swap(i,j);
    if (intersectCount == 2)
        { _(0,1) }
    else if (intersectCount == 3)
        { _(0,2)_(0,1)_(1,2) }
    else if (intersectCount == 4)
        { _(0,2)_(1,3)_(0,1)_(2,3)_(1,2) }
    else if (intersectCount == 5)
        { _(0,3)_(1,4)_(0,2)_(1,3)_(0,1)_(2,4)_(1,2)_(3,4)_(2,3) }
    else if (intersectCount == 6)
        { _(0,5)_(1,3)_(2,4)_(1,2)_(3,4)_(0,3)_(2,5)_(0,1)_(2,3)_(4,5)_(1,2)_(3,4) }
    else if (intersectCount == 7)
        { _(0,6)_(2,3)_(4,5)_(0,2)_(1,4)_(3,6)_(0,1)_(2,5)_(3,4)_(1,2)_(4,6)_(2,3)_(4,5)_(1,2)_(3,4)_(5,6) }
    else if (intersectCount == 8)
        { _(0,2)_(1,3)_(4,6)_(5,7)_(0,4)_(1,5)_(2,6)_(3,7)_(0,1)_(2,3)_(4,5)_(6,7)_(2,4)_(3,5)_(1,4)_(3,6)_(1,2)_(3,4)_(5,6) }
#if MAX_N_INT > 8
    else if (intersectCount == 9)
        { _(0,3)_(1,7)_(2,5)_(4,8)_(0,7)_(2,4)_(3,8)_(5,6)_(0,2)_(1,3)_(4,5)_(7,8)_(1,4)_(3,6)_(5,7)_(0,1)_(2,4)_(3,5)_(6,8)_(2,3)_(4,5)_(6,7)_(1,2)_(3,4)_(5,6) }
    else if (intersectCount == 10)
        { _(0,8)_(1,9)_(2,7)_(3,5)_(4,6)_(0,2)_(1,4)_(5,8)_(7,9)_(0,3)_(2,4)_(5,7)_(6,9)_(0,1)_(3,6)_(8,9)_(1,5)_(2,3)_(4,8)_(6,7)_(1,2)_(3,5)_(4,6)_(7,8)_(2,3)_(4,5)_(6,7)_(3,4)_(5,6) }
    else if (intersectCount == 11)
        { _(0,9)_(1,6)_(2,4)_(3,7)_(5,8)_(0,1)_(3,5)_(4,10)_(6,9)_(7,8)_(1,3)_(2,5)_(4,7)_(8,10)_(0,4)_(1,2)_(3,7)_(5,9)_(6,8)_(0,1)_(2,6)_(4,5)_(7,8)_(9,10)_(2,4)_(3,6)_(5,7)_(8,9)_(1,2)_(3,4)_(5,6)_(7,8)_(2,3)_(4,5)_(6,7) }
    else if (intersectCount == 12)
        { _(0,8)_(1,7)_(2,6)_(3,11)_(4,10)_(5,9)_(0,1)_(2,5)_(3,4)_(6,9)_(7,8)_(10,11)_(0,2)_(1,6)_(5,10)_(9,11)_(0,3)_(1,2)_(4,6)_(5,7)_(8,11)_(9,10)_(1,4)_(3,5)_(6,8)_(7,10)_(1,3)_(2,5)_(6,9)_(8,10)_(2,3)_(4,5)_(6,7)_(8,9)_(4,6)_(5,7)_(3,4)_(5,6)_(7,8) }
#endif  // 8
#if MAX_N_INT > 12
    else if (intersectCount == 13)
        { _(0,11)_(1,7)_(2,4)_(3,5)_(8,9)_(10,12)_(0,2)_(3,6)_(4,12)_(5,7)_(8,10)_(0,8)_(1,3)_(2,5)_(4,9)_(6,11)_(7,12)_(0,1)_(2,10)_(3,8)_(4,6)_(9,11)_(1,3)_(2,4)_(5,10)_(6,8)_(7,9)_(11,12)_(1,2)_(3,4)_(5,8)_(6,9)_(7,10)_(2,3)_(4,7)_(5,6)_(8,11)_(9,10)_(4,5)_(6,7)_(8,9)_(10,11)_(3,4)_(5,6)_(7,8)_(9,10) }
    else if (intersectCount == 14)
        { _(0,1)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(12,13)_(0,2)_(1,3)_(4,8)_(5,9)_(10,12)_(11,13)_(0,4)_(1,2)_(3,7)_(5,8)_(6,10)_(9,13)_(11,12)_(0,6)_(1,5)_(3,9)_(4,10)_(7,13)_(8,12)_(2,10)_(3,11)_(4,6)_(7,9)_(1,3)_(2,8)_(5,11)_(6,7)_(10,12)_(1,4)_(2,6)_(3,5)_(7,11)_(8,10)_(9,12)_(2,4)_(3,6)_(5,8)_(7,10)_(9,11)_(3,4)_(5,6)_(7,8)_(9,10)_(6,7) }
    else if (intersectCount == 15)
        { _(1,2)_(3,10)_(4,14)_(5,8)_(6,13)_(7,12)_(9,11)_(0,14)_(1,5)_(2,8)_(3,7)_(6,9)_(10,12)_(11,13)_(0,7)_(1,6)_(2,9)_(4,10)_(5,11)_(8,13)_(12,14)_(0,6)_(2,4)_(3,5)_(7,11)_(8,10)_(9,12)_(13,14)_(0,3)_(1,2)_(4,7)_(5,9)_(6,8)_(10,11)_(12,13)_(0,1)_(2,3)_(4,6)_(7,9)_(10,12)_(11,13)_(1,2)_(3,5)_(8,10)_(11,12)_(3,4)_(5,6)_(7,8)_(9,10)_(2,3)_(4,5)_(6,7)_(8,9)_(10,11)_(5,6)_(7,8) }
    else if (intersectCount == 16)
        { _(0,13)_(1,12)_(2,15)_(3,14)_(4,8)_(5,6)_(7,11)_(9,10)_(0,5)_(1,7)_(2,9)_(3,4)_(6,13)_(8,14)_(10,15)_(11,12)_(0,1)_(2,3)_(4,5)_(6,8)_(7,9)_(10,11)_(12,13)_(14,15)_(0,2)_(1,3)_(4,10)_(5,11)_(6,7)_(8,9)_(12,14)_(13,15)_(1,2)_(3,12)_(4,6)_(5,7)_(8,10)_(9,11)_(13,14)_(1,4)_(2,6)_(5,8)_(7,10)_(9,13)_(11,14)_(2,4)_(3,6)_(9,12)_(11,13)_(3,5)_(6,8)_(7,9)_(10,12)_(3,4)_(5,6)_(7,8)_(9,10)_(11,12)_(6,7)_(8,9) }
#endif  // 12
    #undef _
#endif

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
