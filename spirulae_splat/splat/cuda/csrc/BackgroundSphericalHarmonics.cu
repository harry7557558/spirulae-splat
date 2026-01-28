#include "BackgroundSphericalHarmonics.cuh"

#include "common.cuh"

#include <gsplat/Common.h>
#include <gsplat/Cameras.cuh>

#include <algorithm>

#include <ATen/ops/empty.h>
#include <ATen/ops/zeros.h>



template<gsplat::CameraModelType CAMERA_MODEL>
__global__ void render_background_sh_forward_kernel(
    const dim3 img_size,
    const float4* intrins,  // fx, fy, cx, cy
    const float* rotation,  // row major 3x3
    const unsigned sh_degree,
    const glm::vec3* __restrict__ sh_coeffs,
    glm::vec3* __restrict__ out_img
) {
    unsigned camera_id = blockIdx.z * blockDim.z + threadIdx.z;

    unsigned i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    int32_t pix_id = (camera_id * img_size.y + i) * img_size.x + j;

    if (i >= img_size.y || j >= img_size.x || camera_id >= img_size.z)
        return;

    intrins += camera_id;
    float fx = intrins[0].x, fy = intrins[0].y, cx = intrins[0].z, cy = intrins[0].w;

    glm::vec2 pos_2d(j+0.5f, i+0.5f);
    CameraRay camera_ray;
    camera_ray.valid_flag = false;

    switch (CAMERA_MODEL) {
    case gsplat::CameraModelType::PINHOLE:
        camera_ray = PerfectPinholeCameraModel(
            { {{img_size.x, img_size.y}, ShutterType::GLOBAL}, {cx, cy}, {fx, fy} }
        ).image_point_to_camera_ray(pos_2d);
        break;
    case gsplat::CameraModelType::FISHEYE:
        camera_ray = OpenCVFisheyeCameraModel(
            { {{img_size.x, img_size.y}, ShutterType::GLOBAL}, {cx, cy}, {fx, fy} }
        ).image_point_to_camera_ray(pos_2d);
        break;
    }
    if (!camera_ray.valid_flag) {
        out_img[pix_id] = {0.0f, 0.0f, 0.0f};
        return;
    }

    float xi = camera_ray.ray_dir.x;
    float yi = -camera_ray.ray_dir.y;
    float zi = -camera_ray.ray_dir.z;
    rotation += 9 * camera_id;
    float x = rotation[0] * xi + rotation[1] * yi + rotation[2] * zi;
    float y = rotation[3] * xi + rotation[4] * yi + rotation[5] * zi;
    float z = rotation[6] * xi + rotation[7] * yi + rotation[8] * zi;

    float xx = x*x, yy = y*y, zz = z*z;

    glm::vec3 color = glm::vec3(0.0f);

    // l0
    color += 0.28209479177387814f * sh_coeffs[0];

    // l1
    if (sh_degree > 1) {
        color += 0.4886025119029199f * y * sh_coeffs[1];
        color += 0.4886025119029199f * z * sh_coeffs[2];
        color += 0.4886025119029199f * x * sh_coeffs[3];
    }

    // l2
    if (sh_degree > 2) {
        color += 1.0925484305920792f * x * y * sh_coeffs[4];
        color += 1.0925484305920792f * y * z * sh_coeffs[5];
        color += (0.9461746957575601f * zz - 0.31539156525251999f) * sh_coeffs[6];
        color += 1.0925484305920792f * x * z * sh_coeffs[7];
        color += 0.5462742152960396f * (xx - yy) * sh_coeffs[8];
    }

    // l3
    if (sh_degree > 3) {
        color += 0.5900435899266435f * y * (3.0f * xx - yy) * sh_coeffs[9];
        color += 2.890611442640554f * x * y * z * sh_coeffs[10];
        color += 0.4570457994644658f * y * (5.0f * zz - 1.0f) * sh_coeffs[11];
        color += 0.3731763325901154f * z * (5.0f * zz - 3.0f) * sh_coeffs[12];
        color += 0.4570457994644658f * x * (5.0f * zz - 1.0f) * sh_coeffs[13];
        color += 1.445305721320277f * z * (xx - yy) * sh_coeffs[14];
        color += 0.5900435899266435f * x * (xx - 3.0f * yy) * sh_coeffs[15];
    }

    // l4
    if (sh_degree > 4) {
        color += 2.5033429417967046f * x * y * (xx - yy) * sh_coeffs[16];
        color += 1.7701307697799304f * y * z * (3.0f * xx - yy) * sh_coeffs[17];
        color += 0.9461746957575601f * x * y * (7.0f * zz - 1.0f) * sh_coeffs[18];
        color += 0.6690465435572892f * y * z * (7.0f * zz - 3.0f) * sh_coeffs[19];
        color += 0.10578554691520431f * (35.0f * zz * zz - 30.0f * zz + 3.0f) * sh_coeffs[20];
        color += 0.6690465435572892f * x * z * (7.0f * zz - 3.0f) * sh_coeffs[21];
        color += 0.47308734787878004f * (xx - yy) * (7.0f * zz - 1.0f) * sh_coeffs[22];
        color += 1.7701307697799304f * x * z * (xx - 3.0f * yy) * sh_coeffs[23];
        color += 0.6258357354491761f * (xx * (xx - 3.0f * yy) - yy * (3.0f * xx - yy)) * sh_coeffs[24];
    }

    color.x = fmaxf(color.x + 0.5f, 0.0f);
    color.y = fmaxf(color.y + 0.5f, 0.0f);
    color.z = fmaxf(color.z + 0.5f, 0.0f);

    out_img[pix_id] = color;
}


template<gsplat::CameraModelType CAMERA_MODEL>
__global__ void render_background_sh_backward_kernel(
    const dim3 img_size,
    const float4* intrins,  // fx, fy, cx, cy
    const float* rotation,  // row major 3x3
    const unsigned sh_degree,
    const glm::vec3* __restrict__ sh_coeffs,
    const glm::vec3* __restrict__ out_color,
    const glm::vec3* __restrict__ v_out_color,
    glm::vec3* __restrict__ v_rotation,
    glm::vec3* __restrict__ v_sh_coeffs
) {
    #if 0
    unsigned camera_id = blockIdx.z * blockDim.z + threadIdx.z;
    unsigned i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    int32_t pix_id = (camera_id * img_size.y + i) * img_size.x + j;
    #else
    unsigned camera_id = blockIdx.y * blockDim.y + threadIdx.y;
    int32_t pix_id = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned j = pix_id % img_size.x, i = pix_id / img_size.x;
    pix_id += camera_id * img_size.y * img_size.x;
    #endif

    bool inside = (i < img_size.y && j < img_size.x && camera_id < img_size.z);

    // backprop color output
    glm::vec3 v_color = glm::vec3(0.0);
    if (inside) {
        glm::vec3 color = out_color[pix_id];
        v_color = v_out_color[pix_id];
        if (color.x == 0.0f || !isfinite(v_color.x)) v_color.x = 0.0f;
        if (color.y == 0.0f || !isfinite(v_color.y)) v_color.y = 0.0f;
        if (color.z == 0.0f || !isfinite(v_color.z)) v_color.z = 0.0f;
    }

    // undistort
    intrins += camera_id;
    float fx = intrins[0].x, fy = intrins[0].y, cx = intrins[0].z, cy = intrins[0].w;

    glm::vec2 pos_2d(j+0.5f, i+0.5f);
    CameraRay camera_ray;
    camera_ray.valid_flag = false;

    switch (CAMERA_MODEL) {
    case gsplat::CameraModelType::PINHOLE:
        camera_ray = PerfectPinholeCameraModel(
            { {{img_size.x, img_size.y}, ShutterType::GLOBAL}, {cx, cy}, {fx, fy} }
        ).image_point_to_camera_ray(pos_2d);
        break;
    case gsplat::CameraModelType::FISHEYE:
        camera_ray = OpenCVFisheyeCameraModel(
            { {{img_size.x, img_size.y}, ShutterType::GLOBAL}, {cx, cy}, {fx, fy} }
        ).image_point_to_camera_ray(pos_2d);
        break;
    }
    if (!camera_ray.valid_flag) {
        inside = false;
        v_color = glm::vec3(0.0);
    }

    // early termination
    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);

    unsigned thread_idx = block.thread_rank();
    // unsigned warp_idx = thread_idx/WARP_SIZE;
    unsigned lane_idx = thread_idx%WARP_SIZE;

    {
        #if 0
        if (__syncthreads_count(inside) == 0)
            return;
        #else
        if (__ballot_sync(~0u, inside) == 0)
            return;
        #endif
    }

    // backward

    float xi = camera_ray.ray_dir.x;
    float yi = -camera_ray.ray_dir.y;
    float zi = -camera_ray.ray_dir.z;
    rotation += 9 * camera_id;
    float x = rotation[0] * xi + rotation[1] * yi + rotation[2] * zi;
    float y = rotation[3] * xi + rotation[4] * yi + rotation[5] * zi;
    float z = rotation[6] * xi + rotation[7] * yi + rotation[8] * zi;

    float xx = x*x, yy = y*y, zz = z*z;

    float v_x = 0.0f, v_y = 0.0f, v_z = 0.0f;
    float v_xx = 0.0f, v_yy = 0.0f, v_zz = 0.0f;

    // __shared__ glm::vec3 atomic_reduce[WARP_SIZE];  // assume WARP_SIZE^2 >= block_size

    glm::vec3 temp3;
    float temp;
    #define _BLOCK_REDUCE_VEC3() \
        warpSum(temp3, warp); \
        if (warp.thread_rank() == 0) \
            atomic_reduce[warp_idx] = temp3; \
        __syncthreads(); \
        temp = 0.0; \
        if (warp_idx < 3 && lane_idx < (blockDim.x*blockDim.y/WARP_SIZE)) \
            temp = atomic_reduce[lane_idx][warp_idx]; \
        warpSum(temp, warp);

    #define _ATOMIC_ADD(address, idx) \
        _BLOCK_REDUCE_VEC3(); \
        if (warp_idx < 3 && lane_idx == 0) \
            atomicAdd((float*)address + (3*idx+warp_idx), temp);

    #undef _ATOMIC_ADD
    #define _ATOMIC_ADD(address, idx) \
        warpSum(temp3, warp); \
        if (lane_idx < 3) { \
            temp = lane_idx == 0 ? temp3.x : lane_idx == 1 ? temp3.y : temp3.z; \
            if (temp != 0.0) atomicAdd((float*)address + (3*idx+lane_idx), temp); \
        }

    // #define _ATOMIC_ADD(address, idx) \
    //     if (inside) { \
    //         atomicAdd((float*)address+3*idx+0, temp3.x); \
    //         atomicAdd((float*)address+3*idx+1, temp3.y); \
    //         atomicAdd((float*)address+3*idx+2, temp3.z); \
    //     }

    // l0
    float v_color_dot_sh_coeff = 0.0f;
    temp3 = 0.28209479177387814f * v_color;
    _ATOMIC_ADD(v_sh_coeffs, 0);

    // l1 - manually calculated
    if (sh_degree > 1) {

        // color += 0.4886025119029199f * y * sh_coeffs[1];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[1]);
        v_y += 0.4886025119029199f * v_color_dot_sh_coeff;
        temp3 = 0.4886025119029199f * y * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 1);

        // color += 0.4886025119029199f * z * sh_coeffs[2];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[2]);
        v_z += 0.4886025119029199f * v_color_dot_sh_coeff;
        temp3 = 0.4886025119029199f * z * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 2);

        // color += 0.4886025119029199f * x * sh_coeffs[3];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[3]);
        v_x += 0.4886025119029199f * v_color_dot_sh_coeff;
        temp3 = 0.4886025119029199f * x * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 3);
    }

    // l2 - manually calculated
    if (sh_degree > 2) {

        // color += 1.0925484305920792f * x * y * sh_coeffs[4];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[4]);
        v_x += 1.0925484305920792f * y * v_color_dot_sh_coeff;
        v_y += 1.0925484305920792f * x * v_color_dot_sh_coeff;
        temp3 = 1.0925484305920792f * x * y * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 4);

        // color += 1.0925484305920792f * y * z * sh_coeffs[5];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[5]);
        v_z += 1.0925484305920792f * y * v_color_dot_sh_coeff;
        v_y += 1.0925484305920792f * z * v_color_dot_sh_coeff;
        temp3 = 1.0925484305920792f * y * z * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 5);

        // color += (0.9461746957575601f * zz - 0.31539156525251999f) * sh_coeffs[6];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[6]);
        v_zz += 0.9461746957575601f * v_color_dot_sh_coeff;
        temp3 = (0.9461746957575601f * zz - 0.31539156525251999f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 6);

        // color += 1.0925484305920792f * x * z * sh_coeffs[7];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[7]);
        v_x += 1.0925484305920792f * z * v_color_dot_sh_coeff;
        v_z += 1.0925484305920792f * x * v_color_dot_sh_coeff;
        temp3 = 1.0925484305920792f * x * z * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 7);

        // color += 0.5462742152960396f * (xx - yy) * sh_coeffs[8];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[8]);
        v_xx += 0.5462742152960396f * v_color_dot_sh_coeff;
        v_yy -= 0.5462742152960396f * v_color_dot_sh_coeff;
        temp3 = 0.5462742152960396f * (xx - yy) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 8);
    }

    // l3 - AI generated, one incorrect line commented
    if (sh_degree > 3) {
        // color += 0.5900435899266435f * y * (3.0f * xx - yy) * sh_coeffs[9];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[9]);
        v_xx += 1.7701307697799305f * y * v_color_dot_sh_coeff;
        v_yy -= 0.5900435899266435f * y * v_color_dot_sh_coeff;
        v_y += 0.5900435899266435f * (3.0f * xx - yy) * v_color_dot_sh_coeff;
        temp3 = 0.5900435899266435f * y * (3.0f * xx - yy) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 9);

        // color += 2.890611442640554f * x * y * z * sh_coeffs[10];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[10]);
        v_x += 2.890611442640554f * y * z * v_color_dot_sh_coeff;
        v_y += 2.890611442640554f * x * z * v_color_dot_sh_coeff;
        v_z += 2.890611442640554f * x * y * v_color_dot_sh_coeff;
        temp3 = 2.890611442640554f * x * y * z * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 10);

        // color += 0.4570457994644658f * y * (5.0f * zz - 1.0f) * sh_coeffs[11];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[11]);
        v_zz += 2.285228997322329f * y * v_color_dot_sh_coeff;
        v_y += 0.4570457994644658f * (5.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        temp3 = 0.4570457994644658f * y * (5.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 11);

        // color += 0.3731763325901154f * z * (5.0f * zz - 3.0f) * sh_coeffs[12];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[12]);
        v_z += 0.3731763325901154f * (5.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_zz += 1.865881662950577f * z * v_color_dot_sh_coeff;
        temp3 = 0.3731763325901154f * z * (5.0f * zz - 3.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 12);

        // color += 0.4570457994644658f * x * (5.0f * zz - 1.0f) * sh_coeffs[13];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[13]);
        v_x += 0.4570457994644658f * (5.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_zz += 2.285228997322329f * x * v_color_dot_sh_coeff;
        temp3 = 0.4570457994644658f * x * (5.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 13);

        // color += 1.445305721320277f * z * (xx - yy) * sh_coeffs[14];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[14]);
        v_xx += 1.445305721320277f * z * v_color_dot_sh_coeff;
        v_yy -= 1.445305721320277f * z * v_color_dot_sh_coeff;
        v_z += 1.445305721320277f * (xx - yy) * v_color_dot_sh_coeff;
        temp3 = 1.445305721320277f * z * (xx - yy) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 14);

        // color += 0.5900435899266435f * x * (xx - 3.0f * yy) * sh_coeffs[15];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[15]);
        // v_xx += 1.1800871798532870f * x * v_color_dot_sh_coeff;
        v_xx += 0.5900435899266435f * x * v_color_dot_sh_coeff;
        v_yy -= 1.7701307697799305f * x * v_color_dot_sh_coeff;
        v_x += 0.5900435899266435f * (xx - 3.0f * yy) * v_color_dot_sh_coeff;
        temp3 = 0.5900435899266435f * x * (xx - 3.0f * yy) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 15);
    }

    // l4 - AI generated, two incorrect lines commented
    if (sh_degree > 4) {
        // color += 2.5033429417967046f * x * y * (xx - yy) * sh_coeffs[16];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[16]);
        v_x += 2.5033429417967046f * y * (xx - yy) * v_color_dot_sh_coeff;
        v_y += 2.5033429417967046f * x * (xx - yy) * v_color_dot_sh_coeff;
        v_xx += 2.5033429417967046f * x * y * v_color_dot_sh_coeff;
        v_yy -= 2.5033429417967046f * x * y * v_color_dot_sh_coeff;
        temp3 = 2.5033429417967046f * x * y * (xx - yy) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 16);

        // color += 1.7701307697799304f * y * z * (3.0f * xx - yy) * sh_coeffs[17];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[17]);
        v_xx += 5.3103923093397912f * y * z * v_color_dot_sh_coeff;
        v_yy -= 1.7701307697799304f * y * z * v_color_dot_sh_coeff;
        v_y += 1.7701307697799304f * z * (3.0f * xx - yy) * v_color_dot_sh_coeff;
        v_z += 1.7701307697799304f * y * (3.0f * xx - yy) * v_color_dot_sh_coeff;
        temp3 = 1.7701307697799304f * y * z * (3.0f * xx - yy) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 17);

        // color += 0.9461746957575601f * x * y * (7.0f * zz - 1.0f) * sh_coeffs[18];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[18]);
        v_x += 0.9461746957575601f * y * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_y += 0.9461746957575601f * x * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_zz += 6.6232228703029207f * x * y * v_color_dot_sh_coeff;
        temp3 = 0.9461746957575601f * x * y * (7.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 18);

        // color += 0.6690465435572892f * y * z * (7.0f * zz - 3.0f) * sh_coeffs[19];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[19]);
        v_y += 0.6690465435572892f * z * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_z += 0.6690465435572892f * y * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_zz += 4.6833258049010244f * y * z * v_color_dot_sh_coeff;
        temp3 = 0.6690465435572892f * y * z * (7.0f * zz - 3.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 19);

        // color += 0.10578554691520431f * (35.0f * zz * zz - 30.0f * zz + 3.0f) * sh_coeffs[20];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[20]);
        v_zz += 0.10578554691520431f * (70.0f * zz - 30.0f) * v_color_dot_sh_coeff;
        temp3 = 0.10578554691520431f * (35.0f * zz * zz - 30.0f * zz + 3.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 20);

        // color += 0.6690465435572892f * x * z * (7.0f * zz - 3.0f) * sh_coeffs[21];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[21]);
        v_x += 0.6690465435572892f * z * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_z += 0.6690465435572892f * x * (7.0f * zz - 3.0f) * v_color_dot_sh_coeff;
        v_zz += 4.6833258049010244f * x * z * v_color_dot_sh_coeff;
        temp3 = 0.6690465435572892f * x * z * (7.0f * zz - 3.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 21);

        // color += 0.47308734787878004f * (xx - yy) * (7.0f * zz - 1.0f) * sh_coeffs[22];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[22]);
        v_xx += 0.47308734787878004f * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_yy -= 0.47308734787878004f * (7.0f * zz - 1.0f) * v_color_dot_sh_coeff;
        v_zz += 3.3116114351514603f * (xx - yy) * v_color_dot_sh_coeff;
        temp3 = 0.47308734787878004f * (xx - yy) * (7.0f * zz - 1.0f) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 22);

        // color += 1.7701307697799304f * x * z * (xx - 3.0f * yy) * sh_coeffs[23];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[23]);
        v_x += 1.7701307697799304f * z * (xx - 3.0f * yy) * v_color_dot_sh_coeff;
        v_z += 1.7701307697799304f * x * (xx - 3.0f * yy) * v_color_dot_sh_coeff;
        v_xx += 1.7701307697799304f * x * z * v_color_dot_sh_coeff;
        v_yy -= 5.3103923093397912f * x * z * v_color_dot_sh_coeff;
        temp3 = 1.7701307697799304f * x * z * (xx - 3.0f * yy) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 23);

        // color += 0.6258357354491761f * (xx * (xx - 3.0f * yy) - yy * (3.0f * xx - yy)) * sh_coeffs[24];
        v_color_dot_sh_coeff = glm::dot(v_color, sh_coeffs[24]);
        // v_xx += 0.6258357354491761f * (4.0f * xx - 6.0f * yy) * v_color_dot_sh_coeff;
        // v_yy += 0.6258357354491761f * (6.0f * yy - 12.0f * xx) * v_color_dot_sh_coeff;
        v_xx += 0.6258357354491761f * (2.0f * xx - 6.0f * yy) * v_color_dot_sh_coeff;
        v_yy += 0.6258357354491761f * (2.0f * yy - 6.0f * xx) * v_color_dot_sh_coeff;
        temp3 = 0.6258357354491761f * (xx * (xx - 3.0f * yy) - yy * (3.0f * xx - yy)) * v_color;
        _ATOMIC_ADD(v_sh_coeffs, 24);
    }

    if (v_rotation == nullptr)
        return;

    v_x += v_xx * 2.0f*x;
    v_y += v_yy * 2.0f*y;
    v_z += v_zz * 2.0f*z;

    glm::vec3 xyz = glm::vec3(x, y, z);
    glm::mat3 dp_dpr = glm::mat3(1.0f) - glm::outerProduct(xyz, xyz);
    glm::vec3 v_p = dp_dpr * glm::vec3(v_x, v_y, v_z);
    v_p *= (inside ? 1.0f : 0.0f);

    v_rotation += 3 * camera_id;
  #if 0
    float tmp[9] = {
        v_p.x * xi, v_p.x * yi, v_p.x * zi,
        v_p.y * xi, v_p.y * yi, v_p.y * zi,
        v_p.z * xi, v_p.z * yi, v_p.z * zi
    };
    #pragma unroll
    for (int i = 0; i < 9; i++) {
        warpSum(tmp[i], warp);
        if (warp.thread_rank() == i)
            atomicAdd((float*)v_rotation + i, tmp[i]);
    }
  #else
    temp3 = glm::vec3(v_p.x * xi, v_p.x * yi, v_p.x * zi);
    _ATOMIC_ADD(v_rotation, 0);
    temp3 = glm::vec3(v_p.y * xi, v_p.y * yi, v_p.y * zi);
    _ATOMIC_ADD(v_rotation, 1);
    temp3 = glm::vec3(v_p.z * xi, v_p.z * yi, v_p.z * zi);
    _ATOMIC_ADD(v_rotation, 2);
  #endif
    #undef _BLOCK_REDUCE_VEC3
    #undef _ATOMIC_ADD
}



/*[AutoHeaderGeneratorExport]*/
at::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    std::string camera_model,
    const at::Tensor &intrins,  // fx, fy, cx, cy
    const at::Tensor &rotation,  // row major 3x3
    const unsigned sh_degree,
    const at::Tensor &sh_coeffs
) {
    DEVICE_GUARD(sh_coeffs);
    CHECK_INPUT(sh_coeffs);
    CHECK_INPUT(intrins);
    CHECK_INPUT(rotation);

    if (intrins.size(-1) != 4) {
        AT_ERROR("intrins must be (..., 4)");
    }
    if (rotation.size(-1) != 3 || rotation.size(-2) != 3) {
        AT_ERROR("rotation must be (..., 3, 3)");
    }
    if (sh_coeffs.size(-2) != sh_degree*sh_degree ||
        sh_coeffs.size(-1) != 3) {
        AT_ERROR("sh_coeffs shape must be (..., sh_regree**2, 3)");
    }

    unsigned b = intrins.numel() / 4;
    if (b != rotation.numel() / 9)
        AT_ERROR("intrins and rotation must have same batch dimension");

    const dim3 img_size = {w, h, b};

    auto options = sh_coeffs.options();
    at::Tensor out_color = at::empty({b, h, w, 3}, options);
    if (b * h * w == 0)
        return out_color;

    if (camera_model == "fisheye") {
        render_background_sh_forward_kernel<gsplat::CameraModelType::FISHEYE>
        <<<_LAUNCH_ARGS_3D(w, h, b, TILE_SIZE, TILE_SIZE, 1)>>>(
            img_size,
            (float4*)intrins.data_ptr<float>(),
            rotation.data_ptr<float>(),
            sh_degree,
            (glm::vec3*)sh_coeffs.contiguous().data_ptr<float>(),
            (glm::vec3*)out_color.contiguous().data_ptr<float>()
        );
    } else {
        render_background_sh_forward_kernel<gsplat::CameraModelType::PINHOLE>
        <<<_LAUNCH_ARGS_3D(w, h, b, TILE_SIZE, TILE_SIZE, 1)>>>(
            img_size,
            (float4*)intrins.data_ptr<float>(),
            rotation.data_ptr<float>(),
            sh_degree,
            (glm::vec3*)sh_coeffs.contiguous().data_ptr<float>(),
            (glm::vec3*)out_color.contiguous().data_ptr<float>()
        );
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    return out_color;
}


/*[AutoHeaderGeneratorExport]*/
std::tuple<
    at::Tensor,  // v_rotation
    at::Tensor  // v_sh_coeffs
> render_background_sh_backward_tensor(
    const unsigned w,
    const unsigned h,
    const std::string camera_model,
    const at::Tensor &intrins,  // fx, fy, cx, cy
    const at::Tensor &rotation,  // row major 3x3
    const unsigned sh_degree,
    const at::Tensor &sh_coeffs,
    const at::Tensor &out_color,
    const at::Tensor &v_out_color
) {
    DEVICE_GUARD(sh_coeffs);
    CHECK_INPUT(sh_coeffs);
    CHECK_INPUT(intrins);
    CHECK_INPUT(rotation);
    CHECK_INPUT(out_color);
    CHECK_INPUT(v_out_color);

    if (intrins.size(-1) != 4) {
        AT_ERROR("intrins must be (..., 4)");
    }
    if (rotation.size(-1) != 3 || rotation.size(-2) != 3) {
        AT_ERROR("rotation must be (..., 3, 3)");
    }
    if (sh_coeffs.size(-2) != sh_degree*sh_degree ||
        sh_coeffs.size(-1) != 3) {
        AT_ERROR("sh_coeffs shape must be (..., sh_regree**2, 3)");
    }
    if (out_color.size(-3) != h ||
        out_color.size(-2) != w ||
        out_color.size(-1) != 3) {
        AT_ERROR("out_color shape must be (..., h, w, 3)");
    }
    if (v_out_color.size(-3) != h ||
        v_out_color.size(-2) != w ||
        v_out_color.size(-1) != 3) {
        AT_ERROR("v_out_color shape must be (... h, w, 3)");
    }
    unsigned b = intrins.numel() / 4;

    // unsigned block_width = TILE_SIZE;
    // unsigned block_width = 32;  // 1024 threads
    const dim3 img_size = {w, h, b};

    auto options = sh_coeffs.options();
    at::Tensor v_rotation = at::zeros({b, 3, 3}, options);
    at::Tensor v_sh_coeffs = at::zeros({sh_degree*sh_degree, 3}, options);
    if (b * h * w == 0)
        return std::make_tuple(v_rotation, v_sh_coeffs);

    if (camera_model == "fisheye") {
        render_background_sh_backward_kernel<gsplat::CameraModelType::FISHEYE>
        // <<<_LAUNCH_ARGS_3D(w, h, b, block_width, block_width, 1)>>>(
        <<<_LAUNCH_ARGS_2D(w*h, b, 256, 1)>>>(
            img_size,
            (float4*)intrins.data_ptr<float>(),
            rotation.data_ptr<float>(),
            sh_degree,
            (glm::vec3*)sh_coeffs.contiguous().data_ptr<float>(),
            (glm::vec3*)out_color.contiguous().data_ptr<float>(),
            (glm::vec3*)v_out_color.contiguous().data_ptr<float>(),
            (glm::vec3*)v_rotation.data_ptr<float>(),
            (glm::vec3*)v_sh_coeffs.data_ptr<float>()
        );
    } else {
        render_background_sh_backward_kernel<gsplat::CameraModelType::PINHOLE>
        // <<<_LAUNCH_ARGS_3D(w, h, b, block_width, block_width, 1)>>>(
        <<<_LAUNCH_ARGS_2D(w*h, b, 256, 1)>>>(
            img_size,
            (float4*)intrins.data_ptr<float>(),
            rotation.data_ptr<float>(),
            sh_degree,
            (glm::vec3*)sh_coeffs.contiguous().data_ptr<float>(),
            (glm::vec3*)out_color.contiguous().data_ptr<float>(),
            (glm::vec3*)v_out_color.contiguous().data_ptr<float>(),
            (glm::vec3*)v_rotation.data_ptr<float>(),
            (glm::vec3*)v_sh_coeffs.data_ptr<float>()
        );
    }
    CHECK_DEVICE_ERROR(cudaGetLastError());

    #undef _TEMP_ARGS

    return std::make_tuple(v_rotation, v_sh_coeffs);
}

