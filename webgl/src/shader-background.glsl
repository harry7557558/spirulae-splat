#version 300 es
precision highp int;
precision highp float;

#include "shader-utils.glsl"
#line 7

uniform vec2 focal;
uniform vec2 viewport;
uniform mat4 projection, view;
uniform int camera_model;
uniform vec4 distortion;
uniform vec4 undistortion;

uniform float sh_degree;
uniform vec3 background_sh[25];

out vec4 fragColor;


void main () {

    vec2 coord0 = vec2(1,-1) * (gl_FragCoord.xy-0.5*viewport) / focal;
    vec2 coord = coord0;
    float vignetting = 0.0;
    if (camera_model == 0) {
        undistort_opencv_iterative(coord0, distortion, coord, vignetting);
    }
    else if (camera_model == 1) {
        // undistort_fisheye_iterative(coord0, distortion, coord, vignetting);
        undistort_fisheye(coord0, distortion, undistortion, coord, vignetting);
    }

    const vec3 vignetting_background = vec3(0.1);
    if (vignetting == 1.0) {
        fragColor = vec4(vignetting_background, 1.0);
        return;
    }

    vec3 dir = normalize(inverse(mat3(view)) * vec3(coord, 1));\
    float x = dir.x, y = dir.y, z = dir.z;
    float xx = x*x, yy = y*y, zz = z*z;

    vec3 color = 0.28209479177387814 * background_sh[0];
    if (sh_degree > 0.5) {
        color += 0.4886025119029199 * y * background_sh[1];
        color += 0.4886025119029199 * z * background_sh[2];
        color += 0.4886025119029199 * x * background_sh[3];
    }
    if (sh_degree > 1.5) {
        color += 1.0925484305920792 * x * y * background_sh[4];
        color += 1.0925484305920792 * y * z * background_sh[5];
        color += (0.9461746957575601 * zz - 0.31539156525251999) * background_sh[6];
        color += 1.0925484305920792 * x * z * background_sh[7];
        color += 0.5462742152960396 * (xx - yy) * background_sh[8];
    }
    if (sh_degree > 2.5) {
        color += 0.5900435899266435 * y * (3. * xx - yy) * background_sh[9];
        color += 2.890611442640554 * x * y * z * background_sh[10];
        color += 0.4570457994644658 * y * (5. * zz - 1.) * background_sh[11];
        color += 0.3731763325901154 * z * (5. * zz - 3.) * background_sh[12];
        color += 0.4570457994644658 * x * (5. * zz - 1.) * background_sh[13];
        color += 1.445305721320277 * z * (xx - yy) * background_sh[14];
        color += 0.5900435899266435 * x * (xx - 3. * yy) * background_sh[15];
    }
    if (sh_degree > 3.5) {
        color += 2.5033429417967046 * x * y * (xx - yy) * background_sh[16];
        color += 1.7701307697799304 * y * z * (3. * xx - yy) * background_sh[17];
        color += 0.9461746957575601 * x * y * (7. * zz - 1.) * background_sh[18];
        color += 0.6690465435572892 * y * z * (7. * zz - 3.) * background_sh[19];
        color += 0.10578554691520431 * (35. * zz * zz - 30. * zz + 3.) * background_sh[20];
        color += 0.6690465435572892 * x * z * (7. * zz - 3.) * background_sh[21];
        color += 0.47308734787878004 * (xx - yy) * (7. * zz - 1.) * background_sh[22];
        color += 1.7701307697799304 * x * z * (xx - 3. * yy) * background_sh[23];
        color += 0.6258357354491761 * (xx * (xx - 3. * yy) - yy * (3. * xx - yy)) * background_sh[24];
    }
    color = max(color+0.5, 0.0);
    if (sh_degree < 0.5)
        color = background_sh[0];

    if (vignetting > 0.0) {
        vignetting = pow(1.0-pow(vignetting,20.0), 2.0);
        color = mix(vignetting_background, color, vignetting);
    }


    // color = vec3(float(iter_count)/16.0);
    // color = vec3(0.5+0.5*(coord-coord0), 0.0);

    // color = 0.5+0.5*dir;
    // float phi = atan(dir.y, dir.x);
    // float theta = atan(length(dir.xy), dir.z);
    // float grid = sin(10.*phi)*sin(10.*theta);
    // color = 0.5+0.5*vec3(clamp(0.1*grid/length(vec2(dFdx(grid),dFdy(grid))),0.,1.));

    fragColor = vec4(color, 1.0);
}
