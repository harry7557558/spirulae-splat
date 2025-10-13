#include "SphericalHarmonics.cuh"

#include <cooperative_groups.h>

#include "common.cuh"

// TODO: gradient to viewdir; masks;

namespace cg = cooperative_groups;

__device__ __constant__ float SH_C0 = 0.28209479177387814f;
__device__ __constant__ float SH_C1 = 0.4886025119029199f;
__device__ __constant__ float SH_C2[] = {
    1.0925484305920792f,
    -1.0925484305920792f,
    0.31539156525252005f,
    -1.0925484305920792f,
    0.5462742152960396f};
__device__ __constant__ float SH_C3[] = {
    -0.5900435899266435f,
    2.890611442640554f,
    -0.4570457994644658f,
    0.3731763325901154f,
    -0.4570457994644658f,
    1.445305721320277f,
    -0.5900435899266435f};
__device__ __constant__ float SH_C4[] = {
    2.5033429417967046f,
    -1.7701307697799304,
    0.9461746957575601f,
    -0.6690465435572892f,
    0.10578554691520431f,
    -0.6690465435572892f,
    0.47308734787878004f,
    -1.7701307697799304f,
    0.6258357354491761f};


// This function is used in both host and device code
__host__ __device__ unsigned num_sh_bases(const unsigned degree) {
    if (degree == 0)
        return 1;
    if (degree == 1)
        return 4;
    if (degree == 2)
        return 9;
    if (degree == 3)
        return 16;
    return 25;
}



// Evaluate spherical harmonics bases at unit direction for high orders using approach described by
// Efficient Spherical Harmonic Evaluation, Peter-Pike Sloan, JCGT 2013
// See https://jcgt.org/published/0002/02/06/ for reference implementation
__device__ void sh_coeffs_to_color_fast(
    const unsigned degree,
    const float3 &viewdir,
    const float3 &coeffs0,
    const float3 *coeffs,
    float3 &colors
) {
    colors = 0.2820947917738781f * coeffs0;
    if (degree < 1) {
        colors = fmax(colors+0.5f, 0.0f);
        return;
    }

    float norm = sqrt(
        viewdir.x * viewdir.x + viewdir.y * viewdir.y + viewdir.z * viewdir.z
    );
    float x = viewdir.x / norm;
    float y = viewdir.y / norm;
    float z = viewdir.z / norm;

    float fTmp0A = 0.48860251190292f;
    colors += fTmp0A *
            (-y * coeffs[0] +
            z * coeffs[1] -
            x * coeffs[2]);
    if (degree < 2) {
        colors = fmax(colors+0.5f, 0.0f);
        return;
    }
    float z2 = z * z;

    float fTmp0B = -1.092548430592079f * z;
    float fTmp1A = 0.5462742152960395f;
    float fC1 = x * x - y * y;
    float fS1 = 2.f * x * y;
    float pSH6 = (0.9461746957575601f * z2 - 0.3153915652525201f);
    float pSH7 = fTmp0B * x;
    float pSH5 = fTmp0B * y;
    float pSH8 = fTmp1A * fC1;
    float pSH4 = fTmp1A * fS1;
    colors +=
        pSH4 * coeffs[3] + pSH5 * coeffs[4] +
        pSH6 * coeffs[5] + pSH7 * coeffs[6] +
        pSH8 * coeffs[7];
    if (degree < 3) {
        colors = fmax(colors+0.5f, 0.0f);
        return;
    }

    float fTmp0C = -2.285228997322329f * z2 + 0.4570457994644658f;
    float fTmp1B = 1.445305721320277f * z;
    float fTmp2A = -0.5900435899266435f;
    float fC2 = x * fC1 - y * fS1;
    float fS2 = x * fS1 + y * fC1;
    float pSH12 = z * (1.865881662950577f * z2 - 1.119528997770346f);
    float pSH13 = fTmp0C * x;
    float pSH11 = fTmp0C * y;
    float pSH14 = fTmp1B * fC1;
    float pSH10 = fTmp1B * fS1;
    float pSH15 = fTmp2A * fC2;
    float pSH9  = fTmp2A * fS2;
    colors += pSH9  * coeffs[8] +
            pSH10 * coeffs[9] +
            pSH11 * coeffs[10] +
            pSH12 * coeffs[11] +
            pSH13 * coeffs[12] +
            pSH14 * coeffs[13] +
            pSH15 * coeffs[14];
    if (degree < 4) {
        colors = fmax(colors+0.5f, 0.0f);
        return;
    }

    float fTmp0D = z * (-4.683325804901025f * z2 + 2.007139630671868f);
    float fTmp1C = 3.31161143515146f * z2 - 0.47308734787878f;
    float fTmp2B = -1.770130769779931f * z;
    float fTmp3A = 0.6258357354491763f;
    float fC3 = x * fC2 - y * fS2;
    float fS3 = x * fS2 + y * fC2;
    float pSH20 = (1.984313483298443f * z * pSH12 - 1.006230589874905f * pSH6);
    float pSH21 = fTmp0D * x;
    float pSH19 = fTmp0D * y;
    float pSH22 = fTmp1C * fC1;
    float pSH18 = fTmp1C * fS1;
    float pSH23 = fTmp2B * fC2;
    float pSH17 = fTmp2B * fS2;
    float pSH24 = fTmp3A * fC3;
    float pSH16 = fTmp3A * fS3;
    colors += pSH16 * coeffs[15] +
            pSH17 * coeffs[16] +
            pSH18 * coeffs[17] +
            pSH19 * coeffs[18] +
            pSH20 * coeffs[19] +
            pSH21 * coeffs[20] +
            pSH22 * coeffs[21] +
            pSH23 * coeffs[22] +
            pSH24 * coeffs[23];
    colors = fmax(colors+0.5f, 0.0f);
}

__device__ void sh_coeffs_to_color_fast_vjp(
    const unsigned degree,
    const float3 &viewdir,
    const float3 *coeffs,
    const float3 &colors,
    float3 v_colors,
    float3 &v_coeffs0,
    float3 *v_coeffs,
    float3 *v_viewdir
) {
    v_colors.x = (colors.x == 0.0f ? 0.0f : v_colors.x);
    v_colors.y = (colors.y == 0.0f ? 0.0f : v_colors.y);
    v_colors.z = (colors.z == 0.0f ? 0.0f : v_colors.z);

    v_coeffs0 = 0.2820947917738781f * v_colors;
    if (degree < 1) {
        return;
    }
    float inorm = rsqrtf(
        viewdir.x * viewdir.x + viewdir.y * viewdir.y + viewdir.z * viewdir.z
    );
    float x = viewdir.x * inorm;
    float y = viewdir.y * inorm;
    float z = viewdir.z * inorm;
    float v_x = 0.f, v_y = 0.f, v_z = 0.f;


    float fTmp0A = 0.48860251190292f;
    v_coeffs[0] = -fTmp0A * y * v_colors;
    v_coeffs[1] = fTmp0A * z * v_colors;
    v_coeffs[2] = -fTmp0A * x * v_colors;
    if (v_viewdir != nullptr) {
        v_x += -fTmp0A * dot(coeffs[2], v_colors);
        v_y += -fTmp0A * dot(coeffs[0], v_colors);
        v_z += fTmp0A * dot(coeffs[1], v_colors);
    }
    if (degree < 2) {
        if (v_viewdir != nullptr) {
            float3 dir_n = make_float3(x, y, z);
            float3 v_dir_n = make_float3(v_x, v_y, v_z);
            float3 v_d = (v_dir_n - dot(v_dir_n, dir_n) * dir_n) * inorm;
            *v_viewdir = v_d;
        }
        return;
    }

    float z2 = z * z;
    float fTmp0B = -1.092548430592079f * z;
    float fTmp1A = 0.5462742152960395f;
    float fC1 = x * x - y * y;
    float fS1 = 2.f * x * y;
    float pSH6 = (0.9461746957575601f * z2 - 0.3153915652525201f);
    float pSH7 = fTmp0B * x;
    float pSH5 = fTmp0B * y;
    float pSH8 = fTmp1A * fC1;
    float pSH4 = fTmp1A * fS1;
    v_coeffs[3] = pSH4 * v_colors;
    v_coeffs[4] = pSH5 * v_colors;
    v_coeffs[5] = pSH6 * v_colors;
    v_coeffs[6] = pSH7 * v_colors;
    v_coeffs[7] = pSH8 * v_colors;

    float fTmp0B_z, fC1_x, fC1_y, fS1_x, fS1_y, pSH6_z, pSH7_x, pSH7_z, pSH5_y,
        pSH5_z, pSH8_x, pSH8_y, pSH4_x, pSH4_y;
    if (v_viewdir != nullptr) {
        fTmp0B_z = -1.092548430592079f;
        fC1_x = 2.f * x;
        fC1_y = -2.f * y;
        fS1_x = 2.f * y;
        fS1_y = 2.f * x;
        pSH6_z = 2.f * 0.9461746957575601f * z;
        pSH7_x = fTmp0B;
        pSH7_z = fTmp0B_z * x;
        pSH5_y = fTmp0B;
        pSH5_z = fTmp0B_z * y;
        pSH8_x = 0.5462742152960395f * fC1_x;
        pSH8_y = 0.5462742152960395f * fC1_y;
        pSH4_x = 0.5462742152960395f * fS1_x;
        pSH4_y = 0.5462742152960395f * fS1_y;

        v_x += dot(v_colors,
            pSH4_x * coeffs[3] + pSH8_x * coeffs[7] +
                pSH7_x * coeffs[6]);
        v_y += dot(v_colors,
            pSH4_y * coeffs[3] + pSH8_y * coeffs[7] +
                pSH5_y * coeffs[4]);
        v_z += dot(v_colors,
            pSH6_z * coeffs[5] + pSH7_z * coeffs[6] +
                pSH5_z * coeffs[4]);
    }

    if (degree < 3) {
        if (v_viewdir != nullptr) {
            float3 dir_n = make_float3(x, y, z);
            float3 v_dir_n = make_float3(v_x, v_y, v_z);
            float3 v_d = (v_dir_n - dot(v_dir_n, dir_n) * dir_n) * inorm;
            *v_viewdir = v_d;
        }
        return;
    }

    float fTmp0C = -2.285228997322329f * z2 + 0.4570457994644658f;
    float fTmp1B = 1.445305721320277f * z;
    float fTmp2A = -0.5900435899266435f;
    float fC2 = x * fC1 - y * fS1;
    float fS2 = x * fS1 + y * fC1;
    float pSH12 = z * (1.865881662950577f * z2 - 1.119528997770346f);
    float pSH13 = fTmp0C * x;
    float pSH11 = fTmp0C * y;
    float pSH14 = fTmp1B * fC1;
    float pSH10 = fTmp1B * fS1;
    float pSH15 = fTmp2A * fC2;
    float pSH9  = fTmp2A * fS2;
    v_coeffs[8] = pSH9 * v_colors;
    v_coeffs[9] = pSH10 * v_colors;
    v_coeffs[10] = pSH11 * v_colors;
    v_coeffs[11] = pSH12 * v_colors;
    v_coeffs[12] = pSH13 * v_colors;
    v_coeffs[13] = pSH14 * v_colors;
    v_coeffs[14] = pSH15 * v_colors;

    float fTmp0C_z, fTmp1B_z, fC2_x, fC2_y, fS2_x, fS2_y, pSH12_z, pSH13_x,
        pSH13_z, pSH11_y, pSH11_z, pSH14_x, pSH14_y, pSH14_z, pSH10_x, pSH10_y,
        pSH10_z, pSH15_x, pSH15_y, pSH9_x, pSH9_y;
    if (v_viewdir != nullptr) {
        fTmp0C_z = -2.285228997322329f * 2.f * z;
        fTmp1B_z = 1.445305721320277f;
        fC2_x = fC1 + x * fC1_x - y * fS1_x;
        fC2_y = x * fC1_y - fS1 - y * fS1_y;
        fS2_x = fS1 + x * fS1_x + y * fC1_x;
        fS2_y = x * fS1_y + fC1 + y * fC1_y;
        pSH12_z = 3.f * 1.865881662950577f * z2 - 1.119528997770346f;
        pSH13_x = fTmp0C;
        pSH13_z = fTmp0C_z * x;
        pSH11_y = fTmp0C;
        pSH11_z = fTmp0C_z * y;
        pSH14_x = fTmp1B * fC1_x;
        pSH14_y = fTmp1B * fC1_y;
        pSH14_z = fTmp1B_z * fC1;
        pSH10_x = fTmp1B * fS1_x;
        pSH10_y = fTmp1B * fS1_y;
        pSH10_z = fTmp1B_z * fS1;
        pSH15_x = -0.5900435899266435f * fC2_x;
        pSH15_y = -0.5900435899266435f * fC2_y;
        pSH9_x = -0.5900435899266435f * fS2_x;
        pSH9_y = -0.5900435899266435f * fS2_y;

        v_x += dot(v_colors,
            pSH9_x * coeffs[8] + pSH15_x * coeffs[14] +
                pSH10_x * coeffs[9] + pSH14_x * coeffs[13] +
                pSH13_x * coeffs[12]);
        v_y += dot(v_colors,
            pSH9_y * coeffs[8] + pSH15_y * coeffs[14] +
                pSH10_y * coeffs[9] + pSH14_y * coeffs[13] +
                pSH11_y * coeffs[10]);
        v_z += dot(v_colors,
            pSH12_z * coeffs[11] + pSH13_z * coeffs[12] +
                pSH11_z * coeffs[10] + pSH14_z * coeffs[13] +
                pSH10_z * coeffs[9]);
    }

    if (degree < 4) {
        if (v_viewdir != nullptr) {
            float3 dir_n = make_float3(x, y, z);
            float3 v_dir_n = make_float3(v_x, v_y, v_z);
            float3 v_d = (v_dir_n - dot(v_dir_n, dir_n) * dir_n) * inorm;
            *v_viewdir = v_d;
        }
        return;
    }

    float fTmp0D = z * (-4.683325804901025f * z2 + 2.007139630671868f);
    float fTmp1C = 3.31161143515146f * z2 - 0.47308734787878f;
    float fTmp2B = -1.770130769779931f * z;
    float fTmp3A = 0.6258357354491763f;
    float fC3 = x * fC2 - y * fS2;
    float fS3 = x * fS2 + y * fC2;
    float pSH20 = (1.984313483298443f * z * pSH12 + -1.006230589874905f * pSH6);
    float pSH21 = fTmp0D * x;
    float pSH19 = fTmp0D * y;
    float pSH22 = fTmp1C * fC1;
    float pSH18 = fTmp1C * fS1;
    float pSH23 = fTmp2B * fC2;
    float pSH17 = fTmp2B * fS2;
    float pSH24 = fTmp3A * fC3;
    float pSH16 = fTmp3A * fS3;
    v_coeffs[15] = pSH16 * v_colors;
    v_coeffs[16] = pSH17 * v_colors;
    v_coeffs[17] = pSH18 * v_colors;
    v_coeffs[18] = pSH19 * v_colors;
    v_coeffs[19] = pSH20 * v_colors;
    v_coeffs[20] = pSH21 * v_colors;
    v_coeffs[21] = pSH22 * v_colors;
    v_coeffs[22] = pSH23 * v_colors;
    v_coeffs[23] = pSH24 * v_colors;

    float fTmp0D_z, fTmp1C_z, fTmp2B_z, fC3_x, fC3_y, fS3_x, fS3_y, pSH20_z,
        pSH21_x, pSH21_z, pSH19_y, pSH19_z, pSH22_x, pSH22_y, pSH22_z, pSH18_x,
        pSH18_y, pSH18_z, pSH23_x, pSH23_y, pSH23_z, pSH17_x, pSH17_y, pSH17_z,
        pSH24_x, pSH24_y, pSH16_x, pSH16_y;
    if (v_viewdir != nullptr) {
        fTmp0D_z = 3.f * -4.683325804901025f * z2 + 2.007139630671868f;
        fTmp1C_z = 2.f * 3.31161143515146f * z;
        fTmp2B_z = -1.770130769779931f;
        fC3_x = fC2 + x * fC2_x - y * fS2_x;
        fC3_y = x * fC2_y - fS2 - y * fS2_y;
        fS3_x = fS2 + y * fC2_x + x * fS2_x;
        fS3_y = x * fS2_y + fC2 + y * fC2_y;
        pSH20_z = 1.984313483298443f * (pSH12 + z * pSH12_z) +
                  -1.006230589874905f * pSH6_z;
        pSH21_x = fTmp0D;
        pSH21_z = fTmp0D_z * x;
        pSH19_y = fTmp0D;
        pSH19_z = fTmp0D_z * y;
        pSH22_x = fTmp1C * fC1_x;
        pSH22_y = fTmp1C * fC1_y;
        pSH22_z = fTmp1C_z * fC1;
        pSH18_x = fTmp1C * fS1_x;
        pSH18_y = fTmp1C * fS1_y;
        pSH18_z = fTmp1C_z * fS1;
        pSH23_x = fTmp2B * fC2_x;
        pSH23_y = fTmp2B * fC2_y;
        pSH23_z = fTmp2B_z * fC2;
        pSH17_x = fTmp2B * fS2_x;
        pSH17_y = fTmp2B * fS2_y;
        pSH17_z = fTmp2B_z * fS2;
        pSH24_x = 0.6258357354491763f * fC3_x;
        pSH24_y = 0.6258357354491763f * fC3_y;
        pSH16_x = 0.6258357354491763f * fS3_x;
        pSH16_y = 0.6258357354491763f * fS3_y;

        v_x += dot(v_colors,
            pSH16_x * coeffs[15] + pSH24_x * coeffs[23] +
                pSH17_x * coeffs[16] + pSH23_x * coeffs[22] +
                pSH18_x * coeffs[17] + pSH22_x * coeffs[21] +
                pSH21_x * coeffs[20]);
        v_y += dot(v_colors,
            pSH16_y * coeffs[15] + pSH24_y * coeffs[23] +
                pSH17_y * coeffs[16] + pSH23_y * coeffs[22] +
                pSH18_y * coeffs[17] + pSH22_y * coeffs[21] +
                pSH19_y * coeffs[18]);
        v_z += dot(v_colors,
            pSH20_z * coeffs[19] + pSH21_z * coeffs[20] +
                pSH19_z * coeffs[18] + pSH22_z * coeffs[21] +
                pSH18_z * coeffs[17] + pSH23_z * coeffs[22] +
                pSH17_z * coeffs[16]);

        float3 dir_n = make_float3(x, y, z);
        float3 v_dir_n = make_float3(v_x, v_y, v_z);
        float3 v_d = (v_dir_n - dot(v_dir_n, dir_n) * dir_n) * inorm;
        *v_viewdir = v_d;
    }
}

__global__ void compute_sh_forward_kernel(
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    const float3* __restrict__ viewdirs,
    const float3* __restrict__ coeffs0,
    const float3* __restrict__ coeffs,
    float3* __restrict__ colors
) {
    unsigned idx = cg::this_grid().thread_rank();
    if (idx >= num_points) {
        return;
    }
    unsigned num_bases = num_sh_bases(degree);
    unsigned idx_sh0 = idx;
    unsigned idx_sh = (num_bases-1) * idx;
    unsigned idx_col = idx;

    sh_coeffs_to_color_fast(
        degrees_to_use, viewdirs[idx], coeffs0[idx_sh0], &coeffs[idx_sh], colors[idx_col]
    );
}

__global__ void compute_sh_backward_kernel(
    const unsigned num_points,
    const unsigned degree,
    const unsigned degrees_to_use,
    const float3* __restrict__ viewdirs,
    const float3* __restrict__ coeffs,
    const float3* __restrict__ colors,
    const float3* __restrict__ v_colors,
    float3* __restrict__ v_coeffs0,
    float3* __restrict__ v_coeffs,
    float3* __restrict__ v_viewdirs
) {
    unsigned idx = cg::this_grid().thread_rank();
    if (idx >= num_points) {
        return;
    }
    unsigned num_bases = num_sh_bases(degree);
    unsigned idx_sh0 = idx;
    unsigned idx_sh = (num_bases-1) * idx;
    unsigned idx_col = idx;
    
    sh_coeffs_to_color_fast_vjp(
        degrees_to_use,
        viewdirs[idx], &coeffs[idx_sh],
        colors[idx_col], v_colors[idx_col],
        v_coeffs0[idx_sh0], &v_coeffs[idx_sh],
        v_viewdirs ? &v_viewdirs[idx] : nullptr
    );
}


torch::Tensor compute_sh_forward_tensor(
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs0,   // [..., 3]
    torch::Tensor &coeffs   // [..., K, 3]
) {
    DEVICE_GUARD(viewdirs);
    CHECK_INPUT(viewdirs);
    CHECK_INPUT(coeffs0);
    CHECK_INPUT(coeffs);
    unsigned num_bases = num_sh_bases(degree);
    long num_points = viewdirs.numel() / 3;

    if (coeffs0.ndimension() < 2 || coeffs0.size(-1) != 3) {
        AT_ERROR("coeffs0 must have dimensions (..., 3)");
    }
    if (coeffs.ndimension() < 3 || coeffs.size(-1) != 3 || coeffs.size(-2) != num_bases-1) {
        AT_ERROR("coeffs must have dimensions (..., D, 3)");
    }
    torch::Tensor colors = torch::empty_like(coeffs0);

    compute_sh_forward_kernel <<<_LAUNCH_ARGS_1D(num_points, 256)>>>(
        num_points, degree, degrees_to_use,
        (float3*)viewdirs.data_ptr<float>(),
        (float3*)coeffs0.data_ptr<float>(),
        (float3*)coeffs.data_ptr<float>(),
        (float3*)colors.data_ptr<float>()
    );

    return colors;
}



std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
compute_sh_backward_tensor(
    const unsigned degree,
    const unsigned degrees_to_use,
    torch::Tensor &viewdirs,  // [..., 3]
    torch::Tensor &coeffs,  // [..., 3]
    torch::Tensor &colors,  // [..., 3]
    torch::Tensor &v_colors  // [..., 3]
) {
    DEVICE_GUARD(viewdirs);
    CHECK_INPUT(viewdirs);
    CHECK_INPUT(coeffs);
    CHECK_INPUT(colors);
    CHECK_INPUT(v_colors);
    unsigned num_bases = num_sh_bases(degree);
    long num_points = viewdirs.numel() / 3;

    if (viewdirs.ndimension() < 2 || viewdirs.size(-1) != 3) {
        AT_ERROR("viewdirs must have dimensions (..., 3)");
    }
    if (colors.ndimension() < 2 || colors.size(-1) != 3) {
        AT_ERROR("colors must have dimensions (..., 3)");
    }
    if (v_colors.ndimension() < 2 || v_colors.size(-1) != 3) {
        AT_ERROR("v_colors must have dimensions (..., 3)");
    }
    torch::Tensor v_coeffs0 = torch::empty_like(colors);
    torch::Tensor v_coeffs = torch::empty_like(coeffs);
    torch::Tensor v_viewdirs = torch::empty_like(viewdirs);

    compute_sh_backward_kernel<<<_LAUNCH_ARGS_1D(num_points, 256)>>>(
        num_points, degree, degrees_to_use,
        (float3 *)viewdirs.data_ptr<float>(),
        (float3*)coeffs.data_ptr<float>(),
        (float3*)colors.data_ptr<float>(),
        (float3*)v_colors.data_ptr<float>(),
        (float3*)v_coeffs0.data_ptr<float>(),
        (float3*)v_coeffs.data_ptr<float>(),
        (float3*)v_viewdirs.data_ptr<float>()
    );

    return std::make_tuple(v_coeffs0, v_coeffs, v_viewdirs);
}

