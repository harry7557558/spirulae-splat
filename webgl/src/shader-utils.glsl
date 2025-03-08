#line 2

#define PI 3.14159265


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

vec2 distort_opencv(in vec2 p, in vec4 dist_coeffs) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float p1 = dist_coeffs.z;
    float p2 = dist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r_dist = 1.0 + k1*r2 + k2*r2*r2;

    float x_dist = x*r_dist + 2.0*p1*x*y + p2*(r2 + 2.0*x*x);
    float y_dist = y*r_dist + p1*(r2 + 2.0*y*y) + 2.0*p2*x*y;
    return vec2(x_dist, y_dist);
}

void distort_opencv_with_jac(
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

vec2 distort_fisheye(in vec2 p, in vec4 dist_coeffs) {
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

    return p * (r_dist/r);
}

void distort_fisheye_with_jac(
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

vec3 distort_jac_3d(
    in vec2 p2d_u, in vec2 p2d_d,
    in mat2 j2d, in vec3 dp
) {
    // transform dp using the Jacobian of (x,y,z) -> (distort(x/z,y/z)*z,z)
    mat3x2 j3d = j2d * mat3x2(1, 0, 0, 1, -p2d_u.x, -p2d_u.y)
        + mat3x2(0, 0, 0, 0, p2d_d.x, p2d_d.y);
    return vec3(j3d * dp, dp.z);
}


int iter_count = 0;

void undistort_opencv_iterative(
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
        distort_opencv_with_jac(p_ud, dist_coeffs, pd, jac);
        vec2 dp = inverse(jac) * (pd-p);
        p_ud -= dp;
        if (length(dp) < 1e-3)
            break;
    }
}

vec2 undistort_opencv_iterative(in vec2 p, in vec4 dist_coeffs) {
    vec2 p_ud = p;
    vec2 pd; mat2 jac;
    for (int i = 0; i < 6; i++) {
        iter_count = i+1;
        distort_opencv_with_jac(p_ud, dist_coeffs, pd, jac);
        vec2 dp = inverse(jac) * (pd-p);
        p_ud -= dp;
        if (length(dp) < 1e-3)
            break;
    }
    return p_ud;
}

void undistort_fisheye_iterative(
    in vec2 p, in vec4 dist_coeffs,
    out vec2 p_ud, out float vignetting
) {
    float fr = fisheye_radius(0.5*PI, dist_coeffs);
    if (length(p) > fr) {
        vignetting = 1.0;
        return;
    }
    vignetting = length(p)/fr;

    p_ud = p;
    vec2 pd; mat2 jac;
    for (int i = 0; i < 12; i++) {
        iter_count = i+1;
        distort_fisheye_with_jac(p_ud, dist_coeffs, pd, jac);
        vec2 dp = inverse(jac) * (pd-p);
        p_ud -= dp;
        if (length(dp) < 1e-3)
            break;
    }
}

void undistort_fisheye(
    in vec2 p, in vec4 dist_coeffs, in vec4 undist_coeffs,
    out vec2 p_ud, out float vignetting
) {
    float fr = fisheye_radius(0.5*PI, dist_coeffs);
    if (length(p) > fr) {
        vignetting = 1.0;
        return;
    }
    vignetting = length(p)/fr;

    float l1 = undist_coeffs.x;
    float l2 = undist_coeffs.y;
    float l3 = undist_coeffs.z;
    float l4 = undist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r = sqrt(r2);
    float theta = r*(1.0 + r2*(l1 + r2*(l2 + r2*(l3 + r2*l4))));

    p_ud = p * (tan(theta)/r);
}

vec2 undistort_fisheye(in vec2 p, in vec4 undist_coeffs) {
    float l1 = undist_coeffs.x;
    float l2 = undist_coeffs.y;
    float l3 = undist_coeffs.z;
    float l4 = undist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r = sqrt(r2);
    float theta = r*(1.0 + r2*(l1 + r2*(l2 + r2*(l3 + r2*l4))));

    return p * (tan(theta)/r);
}
