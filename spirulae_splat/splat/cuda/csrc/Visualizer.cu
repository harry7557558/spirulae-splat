#include <cuda_runtime.h>
#include <cstdint>

#include "generated/slang.cuh"
namespace SlangProjectionUtils {
#include "generated/set_namespace.cuh"
#include "generated/projection_utils.cuh"
}


#include <gsplat/Common.h>
#include <gsplat/Utils.cuh>

#include "types.cuh"
#include "common.cuh"


inline __device__ float ray_aabb_intersection(
    float3 bc, float3 br,
    float3 ro, float3 rd,
    float t0, float t1
) {
    float3 n = (bc - ro) / rd;
    float3 k = br / float3{fabsf(rd.x), fabsf(rd.y), fabsf(rd.z)};
    float tmin = fmaxf(n.x-k.x, fmaxf(n.y-k.y, n.z-k.z));
    float tmax = fminf(n.x+k.x, fminf(n.y+k.y, n.z+k.z));
    bool intersect = tmin < tmax && tmax > t0 && tmin < t1;
    // if (intersect) t0 = fmaxf(t0, tmin), t1 = fminf(t1, tmax);
    return intersect ? tmin : -1.0f;
}

inline __device__ float ray_linear_swept_sphere_intersection(
    float3 a, float3 b,
    float r0, float r1,
    float3 ro, float3 rd,
    float t0, float t1
) {
    float3 q = b-a;
    float inv_l = rsqrtf(dot(q, q));
    float3 u = q * inv_l;
    float k = (r1 - r0) * inv_l;
    float3 m = ro - a;
    float alpha = dot(u, m), beta = dot(u, rd);
    m = m - u * alpha;
    rd = rd - u * beta;
    float A = dot(rd, rd) - (k * beta) * (k * beta);
    float B = dot(m, rd) - (r0 + k * alpha) * (k * beta);
    float C = dot(m, m) - (r0 + k * alpha) * (r0 + k * alpha);
    float delta = B*B - A*C;
    float D = sqrtf(delta);
    float tmin = ((A > 0.0f ? -D : D) - B) / A;
    float tmax = ((A > 0.0f ? D : -D) - B) / A;
    // float topt = (1.0f / inv_l - alpha) / beta;
    // bool intersect = delta > 0.0f && (tmin - topt) * (tmax - topt) < 0.0f && tmax > t0 && tmin < t1;
    float l = 1.0f / inv_l;
    float zmin = alpha + tmin * beta;
    float zmax = alpha + tmax * beta;
    float t = zmin * (zmin - l) < 0.0f ? tmin : tmax;
    bool intersect = delta > 0.0f && fminf(zmin * (zmin - l), zmax * (zmax - l)) < 0.0f && t > t0 && t < t1;
    return intersect ? t : -1.0f;
}

inline __device__ float ray_quad_intersection(
    float3 p, float3 a, float3 b,
    float3 ro, float3 rd,
    float t0, float t1
) {
    float3 po = p - ro;
    float invdet = 1.0f / dot(rd, cross(a, b));
    float t = dot(po, cross(a, b)) * invdet;
    float u = dot(rd, cross(b, po)) * invdet;
    float v = dot(rd, cross(a, po)) * invdet;
    bool intersect = (t > t0 && t < t1 && fmaxf(fabsf(u), fabsf(v)) < 1.0f);
    return intersect ? t : -1.0f;
}

inline __device__ float ray_pinhole_camera_intersection(
    float4 intrin,
    float width, float height,
    float3x3 R, float3 t,
    float size,
    float3 ro, float3 rd,
    float t0, float t1,
    float3& color
) {
    float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
    float3 corners[4] = {
        {-cx / fx, -cy / fy, 1.0f},
        {(width-cx) / fx, -cy / fy, 1.0f},
        {(width-cx) / fx, (height-cy) / fy, 1.0f},
        {-cx / fx, (height-cy) / fy, 1.0f},
    };
    float sxy = sqrtf(fx * fy / (width * height));
    float sz = rsqrtf(sxy);
    sxy *= sz * size;
    const float r0 = 0.015f*size*sz, r1 = 0.03f*size*sz;
    sz = size / sz;
    float3 pmin = t, pmax = t;
    #pragma unroll
    for (int i = 0; i < 4; ++i) {
        float3 corner = {corners[i].x*sxy, corners[i].y*sxy, sz};
        corner = float3{dot(R[0], corner), dot(R[1], corner), dot(R[2], corner)} + t;
        corners[i] = corner;
        pmin.x = fminf(pmin.x, corner.x), pmin.y = fminf(pmin.y, corner.y), pmin.z = fminf(pmin.z, corner.z);
        pmax.x = fmaxf(pmax.x, corner.x), pmax.y = fmaxf(pmax.y, corner.y), pmax.z = fmaxf(pmax.z, corner.z);
    }
    if (ray_aabb_intersection((pmin+pmax)*0.5f, (pmax-pmin)*0.5f + r1, ro, rd, t0, t1) < t0)
        return -1.0f;
    float tmin = ray_quad_intersection(
        (corners[0]+corners[1]+corners[2]+corners[3])*0.25f,
        (corners[0]-corners[1]-corners[2]+corners[3])*0.25f,
        (corners[0]+corners[1]-corners[2]-corners[3])*0.25f,
        ro, rd, t0, t1
    );
    if (tmin > 0.0f)
        color = float3{0.5f, 0.5f, 0.5f};
    #pragma unroll
    for (int i = 0; i < 4; ++i) {
        float ti = ray_linear_swept_sphere_intersection(
            t, corners[i], r0, r1, ro, rd, t0, tmin > 0.0f ? tmin : t1
        );
        tmin = ti < 0.0f ? tmin : ti;
        color = ti < 0.0f ? color : float3{0.1f, 0.1f, 0.1f};
    }
    #pragma unroll
    for (int i = 0; i < 4; ++i) {
        float ti = ray_linear_swept_sphere_intersection(
            corners[i], corners[(i+1)&3], r1, r1, ro, rd, t0, tmin > 0.0f ? tmin : t1
        );
        tmin = ti < 0.0f ? tmin : ti;
        color = ti < 0.0f ? color : float3{0.1f, 0.1f, 0.1f};
    }
    return tmin;
}



__global__ void blit_train_cameras_kernel(
    TensorView<float, 3> render_rgbs,  // [H, W, 3]
    const TensorView<float, 3> render_depths,  // [H, W, 1]
    const TensorView<float, 3> render_alphas,  // [H, W, 1]
    const bool view_is_fisheye,
    const int image_width,
    const int image_height,
    const float4* __restrict__ view_intrins,  // [4]
    const float* __restrict__ view_viewmat,  // [4, 4]
    const CameraDistortionCoeffsBuffer view_dist_coeffs,
    const int num_cameras,
    const float4* __restrict__ intrins,  // [N, 4]
    const int32_t* __restrict__ widths,  // [N]
    const int32_t* __restrict__ heights,  // [N]
    const float* __restrict__ camera_to_worlds  // [N, 3, 4]
) {
    const int pix_x = blockIdx.x * blockDim.x + threadIdx.x;
    const int pix_y = blockIdx.y * blockDim.y + threadIdx.y;
    if (pix_x >= image_width || pix_y >= image_height)
        return;

    constexpr int MSAA = 2;
    float3 ray_o[MSAA*MSAA], ray_d[MSAA*MSAA];
    bool inside = false;
    {
        float4 intrin = view_intrins[0];
        float fx = intrin.x, fy = intrin.y, cx = intrin.z, cy = intrin.w;
        float3x3 R = {  // row major
            view_viewmat[0], view_viewmat[1], view_viewmat[2],  // 1st row
            view_viewmat[4], view_viewmat[5], view_viewmat[6],  // 2nd row
            view_viewmat[8], view_viewmat[9], view_viewmat[10],  // 3rd row
        };
        float3 t = { view_viewmat[3], view_viewmat[7], view_viewmat[11] };
        CameraDistortionCoeffs dist_coeffs = view_dist_coeffs.load(0);

        #pragma unroll
        for (int i = 0; i < MSAA*MSAA; ++i) {
            const float px = (float)pix_x + (i/MSAA + 0.5f) / (float)MSAA;
            const float py = (float)pix_y + (i%MSAA + 0.5f) / (float)MSAA;
            ray_o[i] = SlangProjectionUtils::transform_ray_o(R, t);
            inside |= SlangProjectionUtils::generate_ray(
                {(px-cx)/fx, (py-cy)/fy},
                view_is_fisheye, dist_coeffs,
                &ray_d[i]
            );
            if (!view_is_fisheye)
                ray_d[i] = ray_d[i] * (1.0f / ray_d[i].z);
            ray_d[i] = SlangProjectionUtils::transform_ray_d(R, ray_d[i]);
        }
    }

    if (!inside)
        return;

    float depth = render_depths.load1(pix_y, pix_x);
    float opac = render_alphas.load1(pix_y, pix_x);
    depth /= fminf(fmaxf(opac, 1e-6f) / 0.5f, 1.0f);

    float3 rgb[MSAA*MSAA];
    float alpha[MSAA*MSAA] = {0.0f};
    float tmax[4] = {depth, depth, depth, depth};

    uint32_t time_mem = 0, time_comp = 0;

    for (int cam_idx = 0; cam_idx < num_cameras; ++cam_idx) {
        auto time_start = clock();
        float4 intrin = intrins[cam_idx];
        float w = (float)widths[cam_idx], h = (float)heights[cam_idx];
        float3x3 R = {  // row major
            camera_to_worlds[0], camera_to_worlds[1], camera_to_worlds[2],  // 1st row
            camera_to_worlds[4], camera_to_worlds[5], camera_to_worlds[6],  // 2nd row
            camera_to_worlds[8], camera_to_worlds[9], camera_to_worlds[10],  // 3rd row
        };
        float3 t = { camera_to_worlds[3], camera_to_worlds[7], camera_to_worlds[11] };
        camera_to_worlds += 12;
        time_mem += clock() - time_start;

        time_start = clock();
        for (int i = 0; i < MSAA*MSAA; ++i) {
            // if (alpha[i] == 1.0f)
            //     continue;
            float ti = ray_pinhole_camera_intersection(
                intrin, w, h, R, t, 0.05f,
                ray_o[i], ray_d[i], 0.0f, tmax[i], rgb[i]
            );
            if (ti > 0.0f)
                alpha[i] = 1.0f, tmax[i] = ti;
        }
        time_comp += clock() - time_start;
    }

    float alpha_final = 0.0f;
    float3 rgb_final = {0.0f, 0.0f, 0.0f};
    #pragma unroll
    for (int i = 0; i < MSAA*MSAA; ++i) {
        float a = alpha[i] / (float)(MSAA*MSAA);
        alpha_final += a;
        rgb_final += rgb[i] * a;
    }

    // float t_intersect = ray_aabb_intersection(
    //     float3{0, 0, 1}, float3{1, 2, 3},
    //     ray_o, ray_d,
    //     0.0f, depth
    // );
    // float t_intersect = ray_linear_swept_sphere_intersection(
    //     float3{-1, -2, -3}, float3{1, 2, 3},
    //     0.1f, 0.5f,
    //     ray_o, ray_d,
    //     0.0f, depth
    // );
    if (alpha_final > 0.0f) {
        float3 rgb = render_rgbs.load3(pix_y, pix_x);
        rgb = rgb * (1.0f - alpha_final) + rgb_final;
        render_rgbs.store3(pix_y, pix_x, rgb);
    }
    // else {
    //     float3 rgb = render_rgbs.load3(pix_y, pix_x);
    //     rgb.x = rgb.y = rgb.z = (float)time_mem / (float)(time_mem + time_comp);
    //     render_rgbs.store3(pix_y, pix_x, rgb);
    // }
}


/*[AutoHeaderGeneratorExport]*/
void blit_train_cameras_tensor(
    at::Tensor render_rgbs,  // [H, W, 3]
    at::Tensor render_depths,  // [H, W, 1]
    at::Tensor render_alphas,  // [H, W, 1]
    const bool view_is_fisheye,
    const at::Tensor view_intrins,  // [4]
    const at::Tensor view_viewmat,  // [4, 4]
    const CameraDistortionCoeffsTensor view_dist_coeffs,
    const at::Tensor intrins,  // [N, 4]
    const at::Tensor widths,  // [N]
    const at::Tensor heights,  // [N]
    const at::Tensor camera_to_worlds  // [N, 3, 4]
) {
    DEVICE_GUARD(render_rgbs);
    CHECK_CUDA(render_rgbs);
    CHECK_CUDA(render_depths);
    CHECK_CUDA(render_alphas);
    CHECK_INPUT(view_intrins);
    CHECK_INPUT(view_viewmat);
    CHECK_INPUT(intrins);
    CHECK_INPUT(widths);
    CHECK_INPUT(heights);
    CHECK_INPUT(camera_to_worlds);

    auto h = render_rgbs.size(-3);
    auto w = render_rgbs.size(-2);
    auto n = intrins.size(0);

    blit_train_cameras_kernel<<<_LAUNCH_ARGS_2D(w, h, 8, 4)>>>(
        tensor2view<float, 3>(render_rgbs),
        tensor2view<float, 3>(render_depths),
        tensor2view<float, 3>(render_alphas),
        view_is_fisheye,
        w, h,
        (float4*)view_intrins.data_ptr<float>(),
        view_viewmat.data_ptr<float>(),
        view_dist_coeffs,
        n,
        (float4*)intrins.data_ptr<float>(),
        widths.data_ptr<int32_t>(),
        heights.data_ptr<int32_t>(),
        camera_to_worlds.data_ptr<float>()
    );
    CHECK_DEVICE_ERROR(cudaGetLastError());
}
