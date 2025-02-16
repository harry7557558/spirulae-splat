#version 300 es
precision highp float;

uniform vec2 focal;
uniform vec2 viewport;
uniform mat4 projection, view;
uniform int camera_model;
uniform vec4 distortion;

uniform float sh_degree;
uniform vec3 background_sh[25];

out vec4 fragColor;


// camera distortion models

float opencv_radius(
    in vec4 dist_coeffs
) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;

    // calculate the extrema of dist(r)
    if (k2 == 0.0) {
        float b = -1.0/(3.0*k1);
        return b > 0.0 ? sqrt(b) : -1.0;
    }
    float disc = 9.0*k1*k1-20.0*k2;
    if (disc <= 0.0) return -1.0;
    disc = sqrt(disc);
    float u1 = (-3.0*k1 + disc) / (10.0*k2);
    float u2 = (-3.0*k1 - disc) / (10.0*k2);
    if (u1 <= 0.0) return u2>0.0 ? sqrt(u2) : -1.0;
    return u2>0.0 ? sqrt(min(u1,u2)) : sqrt(u1);
}

float fisheye_radius(
    in float theta, in vec4 dist_coeffs
) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float theta2 = theta*theta;
    return theta*(1.0 + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));
}

void distort_opencv(
    in vec2 p, in vec4 dist_coeffs,
    out vec2 p_dist, out mat2 jac_dist
) {
    // https://www.desmos.com/calculator/adnh2fzdyi

    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float p1 = dist_coeffs.z;
    float p2 = dist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r_dist = 1.0 + k1*r2 + k2*r2*r2;

    float x_dist = x*r_dist + 2.0*p1*x*y + p2*(r2 + 2.0*x*x);
    float y_dist = y*r_dist + p1*(r2 + 2.0*y*y) + 2.0*p2*x*y;
    p_dist = vec2(x_dist, y_dist);

    float jxx = r_dist + 2.0*( x*x*(k1 + 2.0*k2*r2) + p1*y + 3.0*p2*x );
    float jxy = 2.0*( x*y*(k1 + 2.0*k2*r2) + p1*x + p2*y );
    float jyy = r_dist + 2.0*( y*y*(k1 + 2.0*k2*r2) + 3.0*p1*y + p2*x );
    jac_dist = mat2(jxx, jxy, jxy, jyy);
}

void distort_fisheye(
    in vec2 p, in vec4 dist_coeffs,
    out vec2 p_dist, out mat2 jac_dist
) {
    // https://www.desmos.com/calculator/yg7gqrtgm9

    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r = sqrt(r2);
    float theta = atan(r);
    float theta2 = theta*theta;
    float r_dist = theta*(1.0 + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));
    float s_dist = 1.0 + theta2*(3.0*k1 + theta2*(5.0*k2 + theta2*(7.0*k3 + theta2*9.0*k4)));

    float x_dist = x/r * r_dist;
    float y_dist = y/r * r_dist;
    p_dist = vec2(x_dist, y_dist);

    float jd1 = s_dist / (r2*(r2+1.0));
    float jd2 = r_dist / (r2*r);
    float jxx = jd1 * x*x + jd2 * y*y;
    float jxy = jd1 * x*y - jd2 * x*y;
    float jyy = jd1 * y*y + jd2 * x*x;
    jac_dist = mat2(jxx, jxy, jxy, jyy);
}

int iter_count = 0;

void undistort_opencv(
    in vec2 p, in vec4 dist_coeffs,
    out vec2 p_ud, out float vignetting
) {
    float fr = opencv_radius(dist_coeffs);
    if (fr > 0.0 && length(p) > fr) {
        vignetting = 1.0;
        return;
    }
    vignetting = fr > 0.0 ? length(p)/fr : 0.0;

    p_ud = p;
    vec2 pd; mat2 jac;
    for (int i = 0; i < 6; i++) {
        iter_count = i+1;
        distort_opencv(p_ud, dist_coeffs, pd, jac);
        vec2 dp = inverse(jac) * (pd-p);
        p_ud -= dp;
        if (length(dp) < 1e-3)
            break;
    }
}

void undistort_fisheye(
    in vec2 p, in vec4 dist_coeffs,
    out vec2 p_ud, out float vignetting
) {
    float fr = fisheye_radius(1.570796, dist_coeffs);
    if (length(p) > fr) {
        vignetting = 1.0;
        return;
    }
    vignetting = length(p)/fr;

    p_ud = p;
    vec2 pd; mat2 jac;
    for (int i = 0; i < 12; i++) {
        iter_count = i+1;
        distort_fisheye(p_ud, dist_coeffs, pd, jac);
        vec2 dp = inverse(jac) * (pd-p);
        p_ud -= dp;
        if (length(dp) < 1e-3)
            break;
    }
}


void main () {

    vec2 coord0 = vec2(1,-1) * (gl_FragCoord.xy-0.5*viewport) / focal;
    vec2 coord = coord0;
    float vignetting = 0.0;
    if (camera_model == 0)
        undistort_opencv(coord0, distortion, coord, vignetting);
    else if (camera_model == 1)
        undistort_fisheye(coord0, distortion, coord, vignetting);

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
