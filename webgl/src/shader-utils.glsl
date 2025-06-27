#line 2

#define PI 3.14159265

// for per pixel sorting
#define TILE_SIZE 12
#define BBOX_OUTPUT_WIDTH 1536

uint pack_rgba(vec4 cf) {
    // 8x3 bit RGB, 6 bit opacity, 2 bit exp
    cf = max(cf, 0.0);
    float sc = max(cf.x, max(cf.y, cf.z));
    uint sb = uint(log2(clamp(2.0*sc, 1.0, 8.0)));
    uvec4 c = uvec4(
        clamp(cf.xyz * exp2(-float(sb)), 0.0, 1.0) * 255.0,
        clamp(cf.w, 0.0, 1.0) * 64.0
    );
    return (c.x<<24) | (c.y<<16) | (c.z<<8) | (c.w<<2) | sb;
}
vec4 unpack_rgba(uint c) {
    uvec4 cf = uvec4(c>>24, (c>>16)&uint(255), (c>>8)&uint(255), c&uint(255));
    vec3 rgb = (vec3(cf.xyz)+0.5) / 255.0 * exp2(float(cf.w & 3u));
    float a = (float(cf.w>>2)+0.5) / 64.0;
    return vec4(rgb, a);
}

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

vec2 project_fisheye(in vec3 p_view, in vec4 dist_coeffs) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float r = length(p_view.xy);
    float theta = atan(r, p_view.z);
    vec2 uv = p_view.xy * theta / r;

    float r2 = dot(uv, uv);
    float radial = 1.0 + r2*(k1 + r2*(k2 + r2*(k3 + r2*k4)));
    return uv * radial;
}

vec3 unproject_fisheye(in vec2 uv, in vec4 undist_coeffs) {

    float l1 = undist_coeffs.x;
    float l2 = undist_coeffs.y;
    float l3 = undist_coeffs.z;
    float l4 = undist_coeffs.w;

    float r2 = dot(uv, uv);
    uv *= 1.0 + r2*(l1 + r2*(l2 + r2*(l3 + r2*l4)));

    float theta = length(uv);
    float theta2 = theta*theta;
    float k = sin(theta)/theta;
    return vec3(uv.x*k, uv.y*k, cos(theta));
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



mat3 quat_to_rotmat(vec4 quat) {
    // quat to rotation matrix
    float w = quat.x;
    float x = quat.y;
    float y = quat.z;
    float z = quat.w;

    // glm matrices are column-major
    return mat3(
        1.f - 2.f * (y * y + z * z),  // [0][0]
        2.f * (x * y + w * z),  // [0][1]
        2.f * (x * z - w * y),  // [0][2]
        2.f * (x * y - w * z),  // [1][0]
        1.f - 2.f * (x * x + z * z),  // [1][1]
        2.f * (y * z + w * x),  // [1][2]
        2.f * (x * z + w * y),  // [2][0]
        2.f * (y * z - w * x),  // [2][1]
        1.f - 2.f * (x * x + y * y)  // [2][2]
    );
}


// axis-aligned bounding box for perspective camera
bool project_bound_perspective(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy,
    out vec2 center, out vec2 bound
) {
    // 2d conic coefficients
    // A x^2 + 2B xy + C y^2 + 2D x + 2E y + F = 0
    vec3 V01 = cross(V0, V1);
    vec3 V0T = cross(T, V0);
    vec3 V1T = cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = V0T.z * V0T.x + V1T.z * V1T.x - V01.z * V01.x;
    float E = -V01.z * V01.y + V1T.z * V1T.y + V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (B * B - A * C);
    float V = (A * E - B * D) / (B * B - A * C);
    float S = A * U*U + 2.0*B * U*V + C * V*V + 2.0*D * U + 2.0*E * V + F;

    // bounds
    float W = sqrt(C * S / (B * B - A * C));
    float H = sqrt(A * S / (B * B - A * C));
    center = vec2(fx*U+cx, fx*V+cy);
    bound = vec2(fx*W, fy*H);

    return true;
}

// axis-aligned bounding box for OpenCV camera
bool project_bound_opencv(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy, vec4 dist_coeffs,
    out vec2 center, out vec2 bound
) {
    // 2d conic coefficients
    // A x^2 + 2B xy + C y^2 + 2D x + 2E y + F = 0
    vec3 V01 = cross(V0, V1);
    vec3 V0T = cross(T, V0);
    vec3 V1T = cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = V0T.z * V0T.x + V1T.z * V1T.x - V01.z * V01.x;
    float E = -V01.z * V01.y + V1T.z * V1T.y + V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (B * B - A * C);
    float V = (A * E - B * D) / (B * B - A * C);
    float S = A * U*U + 2.0*B * U*V + C * V*V + 2.0*D * U + 2.0*E * V + F;

    // represent ellipse in principal parametric form
    vec2 p0 = vec2(U, V);
    float delta = sqrt((A-C)*(A-C) + 4.0*B*B);
    float eps = min(1e-8*S, -1e-24);
    float lambda1 = 2.0*S/min(delta-(A+C),eps);
    float lambda2 = 2.0*S/min(-delta-(A+C),eps);
    vec2 v0 = normalize(vec2(-A+C+delta, -2.0*B));
    vec2 v1 = sqrt(lambda1) * v0;
    vec2 v2 = sqrt(lambda2) * vec2(-v0.y, v0.x);

    // rejection
    float fr = opencv_radius(dist_coeffs);

    // brute force sampling on parametric
    // see here for quadratic fit code (not faster): https://github.com/harry7557558/spirulae-splat/blob/f7dc2989e5383bbeb08d5723a4801a80f11242ca/webgl/shader-vert.glsl#L238-L280
    const float N = 16.0;
    float x0 = 1e4, x1 = -1e4, y0 = 1e4, y1 = -1e4;
    for (float i = 0.0; i < N; i++) {
        float t = 2.0*PI*i/N;
        vec2 p = p0+v1*cos(t)+v2*sin(t);
        if (dot(p,p) < fr*fr) {
            vec2 c = distort_opencv(p, dist_coeffs);
            x0 = min(x0, c.x); x1 = max(x1, c.x);
            y0 = min(y0, c.y); y1 = max(y1, c.y);
        }
    }
    if (x0 > x1 || y0 > y1)
        return false;

    // update bounds
    vec2 b0 = vec2(x0, y0);
    vec2 b1 = vec2(x1, y1);
    center = 0.5*(b0+b1) * vec2(fx,fy) + vec2(cx,cy);
    bound = 0.5*(b1-b0) * vec2(fx,fy);
    return true;
}

// axis-aligned bounding box for fisheye camera, by sampling on 2D parametric
bool project_bound_fisheye_2d(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy, vec4 dist_coeffs,
    out vec2 center, out vec2 bound
) {
    // 2d conic coefficients
    // A x^2 + 2B xy + C y^2 + 2D x + 2E y + F = 0
    vec3 V01 = cross(V0, V1);
    vec3 V0T = cross(T, V0);
    vec3 V1T = cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = V0T.z * V0T.x + V1T.z * V1T.x - V01.z * V01.x;
    float E = -V01.z * V01.y + V1T.z * V1T.y + V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (B * B - A * C);
    float V = (A * E - B * D) / (B * B - A * C);
    float S = A * U*U + 2.0*B * U*V + C * V*V + 2.0*D * U + 2.0*E * V + F;

    // represent ellipse in principal parametric form
    vec2 p0 = vec2(U, V);
    float delta = sqrt((A-C)*(A-C) + 4.0*B*B);
    float eps = min(1e-8*S, -1e-24);
    float lambda1 = 2.0*S/min(delta-(A+C),eps);
    float lambda2 = 2.0*S/min(-delta-(A+C),eps);
    vec2 v0 = normalize(vec2(-A+C+delta, -2.0*B));
    vec2 v1 = sqrt(lambda1) * v0;
    vec2 v2 = sqrt(lambda2) * vec2(-v0.y, v0.x);

    // brute force sampling on parametric
    const float N = 16.0;
    float x0 = 1e4, x1 = -1e4, y0 = 1e4, y1 = -1e4;
    for (float i = 0.0; i < N; i++) {
        float t = 2.0*PI*i/N;
        vec2 c = distort_fisheye(p0+v1*cos(t)+v2*sin(t), dist_coeffs);
        x0 = min(x0, c.x); x1 = max(x1, c.x);
        y0 = min(y0, c.y); y1 = max(y1, c.y);
    }

    // update bounds
    vec2 b0 = vec2(x0, y0);
    vec2 b1 = vec2(x1, y1);
    center = 0.5*(b0+b1) * vec2(fx,fy) + vec2(cx,cy);
    bound = 0.5*(b1-b0) * vec2(fx,fy);

    return true;
}

// axis-aligned bounding box for fisheye camera, by sampling on 3D parametric
bool project_bound_fisheye(
    vec3 T, vec3 V0, vec3 V1,
    float fx, float fy, float cx, float cy, vec4 dist_coeffs,
    out vec2 center, out vec2 bound
) {
    // sampling
    const float N = 16.0;
    float x0 = 1e4, x1 = -1e4, y0 = 1e4, y1 = -1e4;
    for (float i = 0.0; i < N; i++) {
        float t = 2.0*PI*i/N;
        vec3 p3d = T+V0*cos(t)+V1*sin(t);
        vec2 c = project_fisheye(p3d, dist_coeffs);
        x0 = min(x0, c.x); x1 = max(x1, c.x);
        y0 = min(y0, c.y); y1 = max(y1, c.y);
    }
    if (x0 > x1 || y0 > y1)
        return false;

    // update bounds
    vec2 b0 = vec2(x0, y0);
    vec2 b1 = vec2(x1, y1);
    center = 0.5*(b0+b1) * vec2(fx,fy) + vec2(cx,cy);
    bound = 0.5*(b1-b0) * vec2(fx,fy);
    return true;
}
