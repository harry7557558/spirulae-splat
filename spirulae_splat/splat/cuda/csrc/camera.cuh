#include "helpers.cuh"


inline __device__ float fisheye_radius(
    float theta, float4 dist_coeffs
) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float theta2 = theta*theta;
    return theta*(1.0f + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));
}

inline __device__ glm::vec2 distort_fisheye(glm::vec2 p, float4 dist_coeffs) {
    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r = sqrtf(fmaxf(r2, 1e-20f));
    float theta = atanf(r);
    float theta2 = theta*theta;
    float r_dist = theta*(1.0f + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));

    return p * (r_dist/r);
}


// bound of 2d projection of a 3d ellipse
// 3 components for bound are respectively
// xy radius of AABB and radius of bounding circle
inline __device__ bool project_bound_perspective(
    const glm::vec3 T,
    const glm::vec3 V0,
    const glm::vec3 V1,
    float fx, float fy, float cx, float cy,
    float2 &center, float3 &bound
) {
    // 2d conic coefficients
    glm::vec3 V01 = glm::cross(V0, V1);
    glm::vec3 V0T = glm::cross(T, V0);
    glm::vec3 V1T = glm::cross(T, V1);
    float A = V0T.x * V0T.x + V1T.x * V1T.x - V01.x * V01.x;
    float B = -V01.y * V01.x + V1T.y * V1T.x + V0T.y * V0T.x;
    float C = V0T.y * V0T.y + V1T.y * V1T.y - V01.y * V01.y;
    float D = 2.0f * V0T.z * V0T.x + 2.0f * V1T.z * V1T.x - 2.0f * V01.z * V01.x;
    float E = -2.0f * V01.z * V01.y + 2.0f * V1T.z * V1T.y + 2.0f * V0T.z * V0T.y;
    float F = V0T.z * V0T.z + V1T.z * V1T.z - V01.z * V01.z;

    if (!(B * B < A * C))
        return false;

    // translate to origin
    float U = (C * D - B * E) / (2.0f * (B * B - A * C));
    float V = (A * E - B * D) / (2.0f * (B * B - A * C));
    float S = -(A * U * U + 2.0f * B * U * V + C * V * V + D * U + E * V + F);

    // image transform
    float U_T = fx * U + cx;
    float V_T = fy * V + cy;
    float A_T = A / (fx * fx);
    float B_T = B / (fx * fy);
    float C_T = C / (fy * fy);

    // axis-aligned bounding box
    float W_T = fx * sqrt(C * S / (A * C - B * B));
    float H_T = fy * sqrt(A * S / (A * C - B * B));

    // bounding circle
    float L_T = 0.5f * (A_T + C_T - sqrt((A_T - C_T) * (A_T - C_T) + 4.0f * B_T * B_T));
    float R_T = sqrt(S / L_T);

    // output
    center = {U_T, V_T};
    bound = {W_T, H_T, R_T};
    return true;
}


// axis-aligned bounding box for fisheye camera, by sampling on 3D parametric
// gives better results than first-order approximation, especially at large distortion; requires undistorting in rasterize code
inline __device__ bool project_bound_fisheye(
    const glm::vec3 T,
    const glm::vec3 V0,
    const glm::vec3 V1,
    float fx, float fy, float cx, float cy, float4 dist_coeffs,
    float2 &center, float3 &bound
) {
    // exhausive sampling, appears to be the fastest way on GPU
    // see shader-vert.glsl for prototypes of various methods
    const int N = 16;
    float x0 = 1e4f, x1 = -1e4f, y0 = 1e4f, y1 = -1e4f;
    #pragma unroll
    for (int i = 0; i < N; i++) {
        float t = 2.0f*M_PI*(float)i/(float)N;
        glm::vec3 p3d = T+V0*cos(t)+V1*sin(t);
        if (p3d.z > 1e-4f) {
            glm::vec2 p = {p3d.x/p3d.z, p3d.y/p3d.z};
            glm::vec2 c = distort_fisheye(p, dist_coeffs);
            x0 = fmin(x0, c.x); x1 = fmax(x1, c.x);
            y0 = fmin(y0, c.y); y1 = fmax(y1, c.y);
        }
    }
    if (x0 > x1 || y0 > y1)
        return false;

    // update bounds
    center = { 0.5f*(x0+x1) * fx + cx, 0.5f*(y0+y1) * fy + cy };
    bound = { 0.5f*(x1-x0) * fx, 0.5f*(y1-y0) * fy };
    return true;
}



inline __device__ void distort_fisheye_with_jac(
    glm::vec2 p, float4 dist_coeffs,
    glm::vec2 &p_dist, glm::mat2 &jac_dist
) {
    // https://www.desmos.com/calculator/yg7gqrtgm9

    float k1 = dist_coeffs.x;
    float k2 = dist_coeffs.y;
    float k3 = dist_coeffs.z;
    float k4 = dist_coeffs.w;

    float x = p.x, y = p.y;
    float r2 = x*x + y*y;
    float r = sqrtf(r2);
    float theta = atanf(r);
    float theta2 = theta*theta;
    float r_dist = theta*(1.0f + theta2*(k1 + theta2*(k2 + theta2*(k3 + theta2*k4))));
    float s_dist = 1.0f + theta2*(3.0f*k1 + theta2*(5.0f*k2 + theta2*(7.0f*k3 + theta2*9.0f*k4)));

    float x_dist = x/r * r_dist;
    float y_dist = y/r * r_dist;
    p_dist = glm::vec2(x_dist, y_dist);

    float jd1 = s_dist / (r2*(r2+1.0f));
    float jd2 = r_dist / (r2*r);
    float jxx = jd1 * x*x + jd2 * y*y;
    float jxy = jd1 * x*y - jd2 * x*y;
    float jyy = jd1 * y*y + jd2 * x*x;
    jac_dist = glm::mat2(jxx, jxy, jxy, jyy);
}

inline __device__ glm::vec3 distort_jac_3d(
    glm::vec2 p2d_u, glm::vec2 p2d_d,
    glm::mat2 j2d, glm::vec3 dp
) {
    // transform dp using the Jacobian of (x,y,z) -> (distort(x/z,y/z)*z,z)
    glm::mat3x2 j3d = j2d * glm::mat3x2(1.f, 0.f, 0.f, 1.f, -p2d_u.x, -p2d_u.y)
        + glm::mat3x2(0.f, 0.f, 0.f, 0.f, p2d_d.x, p2d_d.y);
    return glm::vec3(j3d * dp, dp.z);
}

// returns true iff is out
inline __device__ bool undistort_fisheye_iterative(
    glm::vec2 p, float4 dist_coeffs,
    glm::vec2 &p_ud, float tol
) {
    float fr = fisheye_radius(0.5*M_PI, dist_coeffs);
    if (glm::length(p) > fr) {
        return true;
    }

    p_ud = p;
    glm::vec2 pd; glm::mat2 jac;
    for (int i = 0; i < 12; i++) {
        distort_fisheye_with_jac(p_ud, dist_coeffs, pd, jac);
        glm::vec2 dp = inverse(jac) * (pd-p);
        p_ud -= dp;
        if (glm::length(dp) < tol)
            return false;
    }
    return true;
}
