#include "BackgroundSphericalHarmonics.cuh"

#include "common.cuh"
#include "camera.cuh"
#include <algorithm>


enum class CameraType {
    // undistorted vs generic distorted
	Undistorted,
    GenericDistorted,
    // (near-)exact distortion
    OPENCV,
    OPENCV_FISHEYE,
    // approximate distortion
    // same rasterization, distort using Jacobian in projection
    OPENCV_approx,
    OPENCV_FISHEYE_approx,
};



template<CameraType CAMERA_TYPE>
__global__ void render_background_sh_forward_kernel(
    _ARGS_render_background_sh_forward_kernel
) {
    unsigned i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    int32_t pix_id = i * img_size.x + j;

    if (i >= img_size.y || j >= img_size.x) return;

    float fx = intrins.x, fy = intrins.y;
    float cx = intrins.z, cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y)) {
            out_img[pix_id] = {0.0f, 0.0f, 0.0f};
            return;
        }
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }

    float xi = pos_2d.x;
    float yi = -pos_2d.y;
    float zi = -1.0f;
    float xr = rotation[0] * xi + rotation[1] * yi + rotation[2] * zi;
    float yr = rotation[3] * xi + rotation[4] * yi + rotation[5] * zi;
    float zr = rotation[6] * xi + rotation[7] * yi + rotation[8] * zi;
    float norm = sqrtf(fmaxf(xr * xr + yr * yr + zr * zr, 1e-12f));
    float x = isfinite(xr) ? xr / norm : 0.0f;
    float y = isfinite(yr) ? yr / norm : 0.0f;
    float z = isfinite(zr) ? zr / norm : 0.0f;

    float xx = x*x, yy = y*y, zz = z*z;

    glm::vec3 color = glm::vec3(0.0f);
    glm::vec3 *sh_coeffs = (glm::vec3*)sh_coeffs_float3;

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

    out_img[pix_id] = *(float3*)&color;
}


template<CameraType CAMERA_TYPE>
__global__ void render_background_sh_backward_kernel(
    _ARGS_render_background_sh_backward_kernel
) {
    unsigned i = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned j = blockIdx.x * blockDim.x + threadIdx.x;
    int32_t pix_id = i * img_size.x + j;

    bool inside = (i < img_size.y && j < img_size.x);

    unsigned idx = i * img_size.x + j;
    glm::vec3 v_color = glm::vec3(0.0);
    if (inside) {
        glm::vec3 color = ((glm::vec3*)out_color)[idx];
        v_color = ((glm::vec3*)v_out_color)[idx];
        if (color.x == 0.0f || !isfinite(v_color.x)) v_color.x = 0.0f;
        if (color.y == 0.0f || !isfinite(v_color.y)) v_color.y = 0.0f;
        if (color.z == 0.0f || !isfinite(v_color.z)) v_color.z = 0.0f;
        // v_color = glm::clamp(v_color, -glm::vec3(1e4f), glm::vec3(1e4f));
    }

    float fx = intrins.x, fy = intrins.y;
    float cx = intrins.z, cy = intrins.w;

    glm::vec2 pos_2d = { (j + 0.5f - cx) / fx, (i + 0.5f - cy) / fy };
    if (CAMERA_TYPE == CameraType::GenericDistorted && inside) {
        float2 pos_2d_u = undistortion_map[pix_id];
        if (isnan(pos_2d.x+pos_2d.y))
            inside = false;
        else
            pos_2d = { pos_2d_u.x, pos_2d_u.y };
    }
    if (__syncthreads_count(inside) == 0)
        return;

    auto block = cg::this_thread_block();
    cg::thread_block_tile<WARP_SIZE> warp = cg::tiled_partition<WARP_SIZE>(block);

    float xi = pos_2d.x;
    float yi = -pos_2d.y;
    float zi = -1.0f;
    float xr = rotation[0] * xi + rotation[1] * yi + rotation[2] * zi;
    float yr = rotation[3] * xi + rotation[4] * yi + rotation[5] * zi;
    float zr = rotation[6] * xi + rotation[7] * yi + rotation[8] * zi;
    float norm2 = xr * xr + yr * yr + zr * zr;
    float norm = sqrtf(fmaxf(norm2, 1e-12f));
    float x = inside && isfinite(xr) ? xr / norm : 0.0f;
    float y = inside && isfinite(yr) ? yr / norm : 0.0f;
    float z = inside && isfinite(zr) ? zr / norm : 0.0f;

    float xx = x*x, yy = y*y, zz = z*z;

    float v_x = 0.0f, v_y = 0.0f, v_z = 0.0f;
    float v_xx = 0.0f, v_yy = 0.0f, v_zz = 0.0f;

    glm::vec3 *sh_coeffs = (glm::vec3*)sh_coeffs_float3;

    __shared__ glm::vec3 atomic_reduce[WARP_SIZE];  // assume WARP_SIZE^2 >= block_size

    unsigned thread_idx = block.thread_rank();
    unsigned warp_idx = thread_idx/WARP_SIZE;
    unsigned lane_idx = thread_idx%WARP_SIZE;

    glm::vec3 temp3;
    float temp;
    #define _BLOCK_REDUCE_VEC3() \
        warpSum3(temp3, warp); \
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

    v_x += v_xx * 2.0f*x;
    v_y += v_yy * 2.0f*y;
    v_z += v_zz * 2.0f*z;

    glm::vec3 xyz = glm::vec3(x, y, z);
    glm::mat3 dp_dpr = (glm::mat3(1.0f) - glm::outerProduct(xyz, xyz)) / norm;
    glm::vec3 v_p = dp_dpr * glm::vec3(v_x, v_y, v_z);
    v_p *= (inside ? 1.0f : 0.0f);

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



torch::Tensor render_background_sh_forward_tensor(
    const unsigned w,
    const unsigned h,
    std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &rotation,
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs
) {
    DEVICE_GUARD(sh_coeffs);
    CHECK_INPUT(sh_coeffs);
    CHECK_INPUT(rotation);

    if (rotation.numel() != 9) {
        AT_ERROR("rotation must be 3x3");
    }
    if (sh_coeffs.ndimension() != 2 ||
        sh_coeffs.size(0) != sh_degree*sh_degree ||
        sh_coeffs.size(1) != 3) {
        AT_ERROR("sh_coeffs must be (sh_regree**2, 3)");
    }

    const dim3 img_size = {w, h, 1};

    auto options = sh_coeffs.options();
    torch::Tensor out_color = torch::empty({h, w, 3}, options);

    if (camera_model == "") {
        render_background_sh_forward_kernel<CameraType::Undistorted>
        <<<_LAUNCH_ARGS_2D(w, h, TILE_SIZE, TILE_SIZE)>>>(
            img_size,
            tuple2float4(intrins), nullptr,
            rotation.contiguous().data_ptr<float>(),
            sh_degree,
            (float3 *)sh_coeffs.contiguous().data_ptr<float>(),
            (float3 *)out_color.contiguous().data_ptr<float>()
        );
    }

    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);

        render_background_sh_forward_kernel<CameraType::GenericDistorted>
        <<<_LAUNCH_ARGS_2D(w, h, TILE_SIZE, TILE_SIZE)>>>(
            img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            rotation.contiguous().data_ptr<float>(),
            sh_degree,
            (float3 *)sh_coeffs.contiguous().data_ptr<float>(),
            (float3 *)out_color.contiguous().data_ptr<float>()
        );
    }

    return out_color;
}


std::tuple<
    torch::Tensor,  // v_rotation
    torch::Tensor  // v_sh_coeffs
> render_background_sh_backward_tensor(
    const unsigned w,
    const unsigned h,
    const std::string camera_model,
    const std::tuple<float, float, float, float> intrins,
    const std::optional<torch::Tensor> &undistortion_map_,
    const torch::Tensor &rotation,
    const unsigned sh_degree,
    const torch::Tensor &sh_coeffs,
    const torch::Tensor &out_color,
    const torch::Tensor &v_out_color
) {
    DEVICE_GUARD(sh_coeffs);
    CHECK_INPUT(sh_coeffs);
    CHECK_INPUT(rotation);
    CHECK_INPUT(v_out_color);

    if (rotation.numel() != 9) {
        AT_ERROR("rotation must be 3x3");
    }
    if (sh_coeffs.ndimension() != 2 ||
        sh_coeffs.size(0) != sh_degree*sh_degree ||
        sh_coeffs.size(1) != 3) {
        AT_ERROR("sh_coeffs shape must be (sh_regree**2, 3)");
    }
    if (out_color.ndimension() != 3 ||
        out_color.size(0) != h ||
        out_color.size(1) != w ||
        out_color.size(2) != 3) {
        AT_ERROR("out_color shape must be (h, w, 3)");
    }
    if (v_out_color.ndimension() != 3 ||
        v_out_color.size(0) != h ||
        v_out_color.size(1) != w ||
        v_out_color.size(2) != 3) {
        AT_ERROR("v_out_color shape must be (h, w, 3)");
    }

    // unsigned block_width = TILE_SIZE;
    unsigned block_width = 32;  // 1024 threads
    const dim3 img_size = {w, h, 1};

    auto options = sh_coeffs.options();
    torch::Tensor v_rotation = torch::zeros({3, 3}, options);
    torch::Tensor v_sh_coeffs = torch::zeros({sh_degree*sh_degree, 3}, options);

    #define _TEMP_ARGS \
        rotation.contiguous().data_ptr<float>(), \
        sh_degree, \
        (float3 *)sh_coeffs.contiguous().data_ptr<float>(), \
        (float3 *)out_color.contiguous().data_ptr<float>(), \
        (float3 *)v_out_color.contiguous().data_ptr<float>(), \
        (float3 *)v_rotation.contiguous().data_ptr<float>(), \
        (float3 *)v_sh_coeffs.contiguous().data_ptr<float>()

    if (camera_model == "") {
        render_background_sh_backward_kernel<CameraType::Undistorted>
        <<<_LAUNCH_ARGS_2D(w, h, block_width, block_width)>>>(
            img_size,
            tuple2float4(intrins), nullptr,
            _TEMP_ARGS
        );
    }
    else {
        const torch::Tensor& undistortion_map = undistortion_map_.value();
        CHECK_INPUT(undistortion_map);
        render_background_sh_backward_kernel<CameraType::GenericDistorted>
        <<<_LAUNCH_ARGS_2D(w, h, block_width, block_width)>>>(
            img_size,
            tuple2float4(intrins),
            (float2 *)undistortion_map.contiguous().data_ptr<float>(),
            _TEMP_ARGS
        );
    }

    #undef _TEMP_ARGS

    return std::make_tuple(v_rotation, v_sh_coeffs);
}



template __global__ void render_background_sh_forward_kernel<CameraType::Undistorted>(
    _ARGS_render_background_sh_forward_kernel
);
template __global__ void render_background_sh_forward_kernel<CameraType::GenericDistorted>(
    _ARGS_render_background_sh_forward_kernel
);
template __global__ void render_background_sh_backward_kernel<CameraType::Undistorted>(
    _ARGS_render_background_sh_backward_kernel
);
template __global__ void render_background_sh_backward_kernel<CameraType::GenericDistorted>(
    _ARGS_render_background_sh_backward_kernel
);
