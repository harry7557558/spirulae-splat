#pragma once

#include "slang.cuh"

struct CameraDistortion_0
{
    float4  radial_coeffs_0;
    float2  tangential_coeffs_0;
    float2  thin_prism_coeffs_0;
};

inline __device__ CameraDistortion_0 CameraDistortion_x24_syn_dzero_0()
{
    CameraDistortion_0 result_0;
    (&result_0)->radial_coeffs_0 = make_float4 (0.0f);
    float2  _S1 = make_float2 (0.0f);
    (&result_0)->tangential_coeffs_0 = _S1;
    (&result_0)->thin_prism_coeffs_0 = _S1;
    return result_0;
}

inline __device__ CameraDistortion_0 CameraDistortion_x24_syn_dadd_0(CameraDistortion_0 SLANG_anonymous_0_0, CameraDistortion_0 SLANG_anonymous_1_0)
{
    CameraDistortion_0 result_1;
    (&result_1)->radial_coeffs_0 = SLANG_anonymous_0_0.radial_coeffs_0 + SLANG_anonymous_1_0.radial_coeffs_0;
    (&result_1)->tangential_coeffs_0 = SLANG_anonymous_0_0.tangential_coeffs_0 + SLANG_anonymous_1_0.tangential_coeffs_0;
    (&result_1)->thin_prism_coeffs_0 = SLANG_anonymous_0_0.thin_prism_coeffs_0 + SLANG_anonymous_1_0.thin_prism_coeffs_0;
    return result_1;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpy_0, float dOut_0)
{
    DiffPair_float_0 _S2 = *dpx_0;
    float _S3;
    if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
    {
        _S3 = dOut_0;
    }
    else
    {
        if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
        {
            _S3 = 0.0f;
        }
        else
        {
            _S3 = 0.5f * dOut_0;
        }
    }
    dpx_0->primal_0 = _S2.primal_0;
    dpx_0->differential_0 = _S3;
    DiffPair_float_0 _S4 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S2.primal_0))
    {
        _S3 = dOut_0;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_0).primal_0))
        {
            _S3 = 0.0f;
        }
        else
        {
            _S3 = 0.5f * dOut_0;
        }
    }
    dpy_0->primal_0 = _S4.primal_0;
    dpy_0->differential_0 = _S3;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_1, float dOut_1)
{
    float _S5 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_1).primal_0)))))) * dOut_1;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = _S5;
    return;
}

inline __device__ void _d_rsqrt_0(DiffPair_float_0 * dpx_2, float dOut_2)
{
    float _S6 = -0.5f / ((*dpx_2).primal_0 * (F32_sqrt(((*dpx_2).primal_0)))) * dOut_2;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S6;
    return;
}

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_0)
{
    Matrix<float, 3, 3>  result_2;
    int r_0 = int(0);
    for(;;)
    {
        if(r_0 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_0 = int(0);
        for(;;)
        {
            if(c_0 < int(3))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_0)), c_0) = _slang_vector_get_element(x_0.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_2;
}

inline __device__ Matrix<float, 3, 2>  transpose_1(Matrix<float, 2, 3>  x_1)
{
    Matrix<float, 3, 2>  result_3;
    int r_1 = int(0);
    for(;;)
    {
        if(r_1 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_1 = int(0);
        for(;;)
        {
            if(c_1 < int(2))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_3)->rows + (r_1)), c_1) = _slang_vector_get_element(x_1.rows[c_1], r_1);
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_3;
}

inline __device__ Matrix<float, 2, 3>  transpose_2(Matrix<float, 3, 2>  x_2)
{
    Matrix<float, 2, 3>  result_4;
    int r_2 = int(0);
    for(;;)
    {
        if(r_2 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_2 = int(0);
        for(;;)
        {
            if(c_2 < int(3))
            {
            }
            else
            {
                break;
            }
            *_slang_vector_get_element_ptr(((&result_4)->rows + (r_2)), c_2) = _slang_vector_get_element(x_2.rows[c_2], r_2);
            c_2 = c_2 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_4;
}

inline __device__ Matrix<float, 3, 3>  quat_to_rotmat(float4  quat_0)
{
    float x_3 = quat_0.y;
    float inv_norm_0 = (F32_rsqrt((x_3 * x_3 + quat_0.z * quat_0.z + quat_0.w * quat_0.w + quat_0.x * quat_0.x)));
    float x_4 = quat_0.y * inv_norm_0;
    float y_0 = quat_0.z * inv_norm_0;
    float z_0 = quat_0.w * inv_norm_0;
    float w_0 = quat_0.x * inv_norm_0;
    float x2_0 = x_4 * x_4;
    float y2_0 = y_0 * y_0;
    float z2_0 = z_0 * z_0;
    float xy_0 = x_4 * y_0;
    float xz_0 = x_4 * z_0;
    float yz_0 = y_0 * z_0;
    float wx_0 = w_0 * x_4;
    float wy_0 = w_0 * y_0;
    float wz_0 = w_0 * z_0;
    return transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_0 + z2_0), 2.0f * (xy_0 + wz_0), 2.0f * (xz_0 - wy_0), 2.0f * (xy_0 - wz_0), 1.0f - 2.0f * (x2_0 + z2_0), 2.0f * (yz_0 + wx_0), 2.0f * (xz_0 + wy_0), 2.0f * (yz_0 - wx_0), 1.0f - 2.0f * (x2_0 + y2_0)));
}

struct DiffPair_matrixx3Cfloatx2C3x2C3x3E_0
{
    Matrix<float, 3, 3>  primal_0;
    Matrix<float, 3, 3>  differential_0;
};

struct DiffPair_vectorx3Cfloatx2C3x3E_0
{
    float3  primal_0;
    float3  differential_0;
};

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_3)
{
    float _S7 = (*left_0).primal_0.rows[int(0)].x * dOut_3.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_3.x;
    float sum_0 = _S7 + (*left_0).primal_0.rows[int(1)].x * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_3.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_3.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_3.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S8 = (*left_0).primal_0.rows[int(0)].y * dOut_3.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_3.x;
    float sum_2 = _S8 + (*left_0).primal_0.rows[int(1)].y * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_3.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_3.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_3.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S9 = (*left_0).primal_0.rows[int(0)].z * dOut_3.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_3.x;
    float sum_4 = _S9 + (*left_0).primal_0.rows[int(1)].z * dOut_3.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_3.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_3.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_3.z;
    *&((&right_d_result_0)->z) = sum_5;
    left_0->primal_0 = (*left_0).primal_0;
    left_0->differential_0 = left_d_result_0;
    right_0->primal_0 = (*right_0).primal_0;
    right_0->differential_0 = right_d_result_0;
    return;
}

inline __device__ float3  mul_0(Matrix<float, 3, 3>  left_1, float3  right_1)
{
    float3  result_5;
    int i_0 = int(0);
    for(;;)
    {
        if(i_0 < int(3))
        {
        }
        else
        {
            break;
        }
        int j_0 = int(0);
        float sum_6 = 0.0f;
        for(;;)
        {
            if(j_0 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_7 = sum_6 + _slang_vector_get_element(left_1.rows[i_0], j_0) * _slang_vector_get_element(right_1, j_0);
            j_0 = j_0 + int(1);
            sum_6 = sum_7;
        }
        *_slang_vector_get_element_ptr(&result_5, i_0) = sum_6;
        i_0 = i_0 + int(1);
    }
    return result_5;
}

inline __device__ void posW2C(Matrix<float, 3, 3>  R_0, float3  t_0, float3  pW_0, float3  * pC_0)
{
    *pC_0 = mul_0(R_0, pW_0) + t_0;
    return;
}

inline __device__ void mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, Matrix<float, 3, 3>  dOut_4)
{
    Matrix<float, 3, 3>  left_d_result_1;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_1;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_4.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_4.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(0)].x * dOut_4.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_4.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(0)].y * dOut_4.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_4.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(0)].z * dOut_4.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_4.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_4.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_4.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(1)].x * dOut_4.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(1)].y * dOut_4.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_4.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(1)].z * dOut_4.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_4.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_4.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_4.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_4.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_4.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(2)].x * dOut_4.rows[int(2)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_4.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(2)].y * dOut_4.rows[int(2)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(2)].z * dOut_4.rows[int(2)].z;
    left_2->primal_0 = (*left_2).primal_0;
    left_2->differential_0 = left_d_result_1;
    right_2->primal_0 = (*right_2).primal_0;
    right_2->differential_0 = right_d_result_1;
    return;
}

struct DiffPair_matrixx3Cfloatx2C2x2C3x3E_0
{
    Matrix<float, 2, 3>  primal_0;
    Matrix<float, 2, 3>  differential_0;
};

inline __device__ void mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_3, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_3, Matrix<float, 2, 3>  dOut_5)
{
    Matrix<float, 2, 3>  left_d_result_2;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = 0.0f;
    Matrix<float, 3, 3>  right_d_result_2;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = 0.0f;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = 0.0f;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].x * dOut_5.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].y * dOut_5.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_3).primal_0.rows[int(0)].x * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_3).primal_0.rows[int(0)].y * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].z * dOut_5.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_3).primal_0.rows[int(0)].z * dOut_5.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].x * dOut_5.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].y * dOut_5.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_3).primal_0.rows[int(1)].x * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_3).primal_0.rows[int(1)].y * dOut_5.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].z * dOut_5.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_3).primal_0.rows[int(1)].z * dOut_5.rows[int(1)].z;
    left_3->primal_0 = (*left_3).primal_0;
    left_3->differential_0 = left_d_result_2;
    right_3->primal_0 = (*right_3).primal_0;
    right_3->differential_0 = right_d_result_2;
    return;
}

struct DiffPair_matrixx3Cfloatx2C3x2C2x3E_0
{
    Matrix<float, 3, 2>  primal_0;
    Matrix<float, 3, 2>  differential_0;
};

inline __device__ void mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_4, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * right_4, Matrix<float, 2, 2>  dOut_6)
{
    Matrix<float, 2, 3>  left_d_result_3;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = 0.0f;
    Matrix<float, 3, 2>  right_d_result_3;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = 0.0f;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = 0.0f;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = *&(((&left_d_result_3)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_6.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = *&(((&right_d_result_3)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(0)].x * dOut_6.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = *&(((&left_d_result_3)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_6.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = *&(((&right_d_result_3)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(0)].y * dOut_6.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = *&(((&left_d_result_3)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_6.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = *&(((&right_d_result_3)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(0)].z * dOut_6.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = *&(((&left_d_result_3)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_6.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = *&(((&right_d_result_3)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(0)].x * dOut_6.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = *&(((&left_d_result_3)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_6.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = *&(((&right_d_result_3)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(0)].y * dOut_6.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = *&(((&left_d_result_3)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_6.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = *&(((&right_d_result_3)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(0)].z * dOut_6.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = *&(((&left_d_result_3)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_6.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = *&(((&right_d_result_3)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(1)].x * dOut_6.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = *&(((&left_d_result_3)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_6.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = *&(((&right_d_result_3)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(1)].y * dOut_6.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = *&(((&left_d_result_3)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_6.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = *&(((&right_d_result_3)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(1)].z * dOut_6.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = *&(((&left_d_result_3)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_6.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = *&(((&right_d_result_3)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(1)].x * dOut_6.rows[int(1)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = *&(((&left_d_result_3)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_6.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = *&(((&right_d_result_3)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(1)].y * dOut_6.rows[int(1)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = *&(((&left_d_result_3)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_6.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = *&(((&right_d_result_3)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(1)].z * dOut_6.rows[int(1)].y;
    left_4->primal_0 = (*left_4).primal_0;
    left_4->differential_0 = left_d_result_3;
    right_4->primal_0 = (*right_4).primal_0;
    right_4->differential_0 = right_d_result_3;
    return;
}

inline __device__ Matrix<float, 3, 3>  mul_4(Matrix<float, 3, 3>  left_5, Matrix<float, 3, 3>  right_5)
{
    Matrix<float, 3, 3>  result_6;
    int r_3 = int(0);
    for(;;)
    {
        if(r_3 < int(3))
        {
        }
        else
        {
            break;
        }
        int c_3 = int(0);
        for(;;)
        {
            if(c_3 < int(3))
            {
            }
            else
            {
                break;
            }
            int i_1 = int(0);
            float sum_8 = 0.0f;
            for(;;)
            {
                if(i_1 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_9 = sum_8 + _slang_vector_get_element(left_5.rows[r_3], i_1) * _slang_vector_get_element(right_5.rows[i_1], c_3);
                i_1 = i_1 + int(1);
                sum_8 = sum_9;
            }
            *_slang_vector_get_element_ptr(((&result_6)->rows + (r_3)), c_3) = sum_8;
            c_3 = c_3 + int(1);
        }
        r_3 = r_3 + int(1);
    }
    return result_6;
}

inline __device__ Matrix<float, 2, 3>  mul_5(Matrix<float, 2, 3>  left_6, Matrix<float, 3, 3>  right_6)
{
    Matrix<float, 2, 3>  result_7;
    int r_4 = int(0);
    for(;;)
    {
        if(r_4 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_4 = int(0);
        for(;;)
        {
            if(c_4 < int(3))
            {
            }
            else
            {
                break;
            }
            int i_2 = int(0);
            float sum_10 = 0.0f;
            for(;;)
            {
                if(i_2 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_11 = sum_10 + _slang_vector_get_element(left_6.rows[r_4], i_2) * _slang_vector_get_element(right_6.rows[i_2], c_4);
                i_2 = i_2 + int(1);
                sum_10 = sum_11;
            }
            *_slang_vector_get_element_ptr(((&result_7)->rows + (r_4)), c_4) = sum_10;
            c_4 = c_4 + int(1);
        }
        r_4 = r_4 + int(1);
    }
    return result_7;
}

inline __device__ Matrix<float, 2, 2>  mul_6(Matrix<float, 2, 3>  left_7, Matrix<float, 3, 2>  right_7)
{
    Matrix<float, 2, 2>  result_8;
    int r_5 = int(0);
    for(;;)
    {
        if(r_5 < int(2))
        {
        }
        else
        {
            break;
        }
        int c_5 = int(0);
        for(;;)
        {
            if(c_5 < int(2))
            {
            }
            else
            {
                break;
            }
            int i_3 = int(0);
            float sum_12 = 0.0f;
            for(;;)
            {
                if(i_3 < int(3))
                {
                }
                else
                {
                    break;
                }
                float sum_13 = sum_12 + _slang_vector_get_element(left_7.rows[r_5], i_3) * _slang_vector_get_element(right_7.rows[i_3], c_5);
                i_3 = i_3 + int(1);
                sum_12 = sum_13;
            }
            *_slang_vector_get_element_ptr(((&result_8)->rows + (r_5)), c_5) = sum_12;
            c_5 = c_5 + int(1);
        }
        r_5 = r_5 + int(1);
    }
    return result_8;
}

inline __device__ void covarW2C(Matrix<float, 3, 3>  R_1, Matrix<float, 3, 3>  covarW_0, Matrix<float, 3, 3>  * covarC_0)
{
    *covarC_0 = mul_4(mul_4(R_1, covarW_0), transpose_0(R_1));
    return;
}

inline __device__ void quat_scale_to_covar(float4  quat_1, float3  scale_0, Matrix<float, 3, 3>  * covar_0)
{
    float x_5 = quat_1.y;
    float inv_norm_1 = (F32_rsqrt((x_5 * x_5 + quat_1.z * quat_1.z + quat_1.w * quat_1.w + quat_1.x * quat_1.x)));
    float x_6 = quat_1.y * inv_norm_1;
    float y_1 = quat_1.z * inv_norm_1;
    float z_1 = quat_1.w * inv_norm_1;
    float w_1 = quat_1.x * inv_norm_1;
    float x2_1 = x_6 * x_6;
    float y2_1 = y_1 * y_1;
    float z2_1 = z_1 * z_1;
    float xy_1 = x_6 * y_1;
    float xz_1 = x_6 * z_1;
    float yz_1 = y_1 * z_1;
    float wx_1 = w_1 * x_6;
    float wy_1 = w_1 * y_1;
    float wz_1 = w_1 * z_1;
    Matrix<float, 3, 3>  M_0 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_4(M_0, transpose_0(M_0));
    return;
}

inline __device__ void quat_scale_to_sqrt_covar(float4  quat_2, float3  scale_1, Matrix<float, 3, 3>  * M_1)
{
    float x_7 = quat_2.y;
    float inv_norm_2 = (F32_rsqrt((x_7 * x_7 + quat_2.z * quat_2.z + quat_2.w * quat_2.w + quat_2.x * quat_2.x)));
    float x_8 = quat_2.y * inv_norm_2;
    float y_2 = quat_2.z * inv_norm_2;
    float z_2 = quat_2.w * inv_norm_2;
    float w_2 = quat_2.x * inv_norm_2;
    float x2_2 = x_8 * x_8;
    float y2_2 = y_2 * y_2;
    float z2_2 = z_2 * z_2;
    float xy_2 = x_8 * y_2;
    float xz_2 = x_8 * z_2;
    float yz_2 = y_2 * z_2;
    float wx_2 = w_2 * x_8;
    float wy_2 = w_2 * y_2;
    float wz_2 = w_2 * z_2;
    *M_1 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
    return;
}

inline __device__ Matrix<float, 2, 2>  inverse(Matrix<float, 2, 2>  m_0)
{
    float invdet_0 = 1.0f / (m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x);
    return makeMatrix<float, 2, 2> (m_0.rows[int(1)].y * invdet_0, - m_0.rows[int(0)].y * invdet_0, - m_0.rows[int(1)].x * invdet_0, m_0.rows[int(0)].x * invdet_0);
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_3, DiffPair_float_0 * dpy_1, float dOut_7)
{
    DiffPair_float_0 _S10 = *dpx_3;
    float _S11;
    if(((*dpx_3).primal_0) < ((*dpy_1).primal_0))
    {
        _S11 = dOut_7;
    }
    else
    {
        if(((*dpx_3).primal_0) > ((*dpy_1).primal_0))
        {
            _S11 = 0.0f;
        }
        else
        {
            _S11 = 0.5f * dOut_7;
        }
    }
    dpx_3->primal_0 = _S10.primal_0;
    dpx_3->differential_0 = _S11;
    DiffPair_float_0 _S12 = *dpy_1;
    if(((*dpy_1).primal_0) < (_S10.primal_0))
    {
        _S11 = dOut_7;
    }
    else
    {
        if(((*dpy_1).primal_0) > ((*dpx_3).primal_0))
        {
            _S11 = 0.0f;
        }
        else
        {
            _S11 = 0.5f * dOut_7;
        }
    }
    dpy_1->primal_0 = _S12.primal_0;
    dpy_1->differential_0 = _S11;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S13 = float(width_0);
    float _S14 = float(height_0);
    float _S15 = 0.30000001192092896f * (0.5f * _S13 / fx_0);
    float _S16 = 0.30000001192092896f * (0.5f * _S14 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_0 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S13 - cx_0) / fx_0 + _S15), ((F32_max((- (cx_0 / fx_0 + _S15)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S14 - cy_0) / fy_0 + _S16), ((F32_max((- (cy_0 / fy_0 + _S16)), (mean3d_0.y * rz_0))))))) * rz2_0);
    *cov2d_0 = mul_6(mul_5(J_0, cov3d_0), transpose_1(J_0));
    *mean2d_0 = make_float2 (fx_0 * mean3d_0.x * rz_0 + cx_0, fy_0 * mean3d_0.y * rz_0 + cy_0);
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float fx_1, float fy_1, float cx_1, float cy_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float rz_1 = 1.0f / mean3d_1.z;
    float rz2_1 = rz_1 * rz_1;
    Matrix<float, 2, 3>  J_1 = makeMatrix<float, 2, 3> (fx_1 * rz_1, 0.0f, - fx_1 * mean3d_1.x * rz2_1, 0.0f, fy_1 * rz_1, - fy_1 * mean3d_1.y * rz2_1);
    *cov2d_1 = mul_6(mul_5(J_1, cov3d_1), transpose_1(J_1));
    *mean2d_1 = make_float2 (fx_1 * mean3d_1.x * rz_1 + cx_1, fy_1 * mean3d_1.y * rz_1 + cy_1);
    return;
}

inline __device__ CameraDistortion_0 CameraDistortion_x24init_0(float4  radial_coeffs_1, float2  tangential_coeffs_1, float2  thin_prism_coeffs_1)
{
    CameraDistortion_0 _S17;
    (&_S17)->radial_coeffs_0 = radial_coeffs_1;
    (&_S17)->tangential_coeffs_0 = tangential_coeffs_1;
    (&_S17)->thin_prism_coeffs_0 = thin_prism_coeffs_1;
    return _S17;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, float dOut_8)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_2).primal_0.x * dOut_8;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_4).primal_0.x * dOut_8;
    *&((&x_d_result_0)->y) = (*dpy_2).primal_0.y * dOut_8;
    *&((&y_d_result_0)->y) = (*dpx_4).primal_0.y * dOut_8;
    *&((&x_d_result_0)->z) = (*dpy_2).primal_0.z * dOut_8;
    *&((&y_d_result_0)->z) = (*dpx_4).primal_0.z * dOut_8;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = x_d_result_0;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_0;
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void _d_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_5, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpy_3, float dOut_9)
{
    float2  x_d_result_1;
    *&((&x_d_result_1)->x) = (*dpy_3).primal_0.x * dOut_9;
    float2  y_d_result_1;
    *&((&y_d_result_1)->x) = (*dpx_5).primal_0.x * dOut_9;
    *&((&x_d_result_1)->y) = (*dpy_3).primal_0.y * dOut_9;
    *&((&y_d_result_1)->y) = (*dpx_5).primal_0.y * dOut_9;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = x_d_result_1;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = y_d_result_1;
    return;
}

inline __device__ float dot_0(float3  x_9, float3  y_3)
{
    int i_4 = int(0);
    float result_9 = 0.0f;
    for(;;)
    {
        if(i_4 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_10 = result_9 + _slang_vector_get_element(x_9, i_4) * _slang_vector_get_element(y_3, i_4);
        i_4 = i_4 + int(1);
        result_9 = result_10;
    }
    return result_9;
}

inline __device__ float dot_1(float2  x_10, float2  y_4)
{
    int i_5 = int(0);
    float result_11 = 0.0f;
    for(;;)
    {
        if(i_5 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_12 = result_11 + _slang_vector_get_element(x_10, i_5) * _slang_vector_get_element(y_4, i_5);
        i_5 = i_5 + int(1);
        result_11 = result_12;
    }
    return result_11;
}

inline __device__ float length_0(float2  x_11)
{
    return (F32_sqrt((dot_1(x_11, x_11))));
}

inline __device__ float length_1(float3  x_12)
{
    return (F32_sqrt((dot_0(x_12, x_12))));
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_6, float dOut_10)
{
    DiffPair_float_0 _S18 = *dpx_6;
    float _S19 = - (*dpy_4).primal_0 / ((*dpx_6).primal_0 * (*dpx_6).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S19;
    float _S20 = _S18.primal_0 / (_S18.primal_0 * _S18.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S20;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S21, float _S22)
{
    return (F32_atan2((_S21), (_S22)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S23, DiffPair_float_0 * _S24, float _S25)
{
    _d_atan2_0(_S23, _S24, _S25);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S26, float _S27)
{
    _d_sqrt_0(_S26, _S27);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_7, float _s_dOut_0)
{
    float _S28 = (*dpx_7).primal_0.x;
    float _S29 = (*dpx_7).primal_0.y;
    DiffPair_float_0 _S30;
    (&_S30)->primal_0 = _S28 * _S28 + _S29 * _S29;
    (&_S30)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S30, _s_dOut_0);
    float _S31 = (*dpx_7).primal_0.y * _S30.differential_0;
    float _S32 = _S31 + _S31;
    float _S33 = (*dpx_7).primal_0.x * _S30.differential_0;
    float _S34 = _S33 + _S33;
    float2  _S35 = make_float2 (0.0f);
    *&((&_S35)->y) = _S32;
    *&((&_S35)->x) = _S34;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S35;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S36, float _S37)
{
    s_bwd_prop_length_impl_0(_S36, _S37);
    return;
}

inline __device__ void fisheye_proj_3dgs(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float fx_2, float fy_2, float cx_2, float cy_2, float4  radial_coeffs_2, float2  tangential_coeffs_2, float2  thin_prism_coeffs_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    CameraDistortion_0 dist_coeffs_0 = CameraDistortion_x24init_0(radial_coeffs_2, tangential_coeffs_2, thin_prism_coeffs_2);
    float2  _S38 = float2 {mean3d_2.x, mean3d_2.y};
    float r_6 = length_0(_S38);
    float _S39 = mean3d_2.z;
    float theta_0 = (F32_atan2((r_6), (_S39)));
    float k_0;
    if(theta_0 < 0.00100000004749745f)
    {
        k_0 = (1.0f - theta_0 * theta_0 / 3.0f) / _S39;
    }
    else
    {
        k_0 = theta_0 / r_6;
    }
    float2  _S40 = _S38 * make_float2 (k_0);
    float k1_0 = dist_coeffs_0.radial_coeffs_0.x;
    float k2_0 = dist_coeffs_0.radial_coeffs_0.y;
    float k3_0 = dist_coeffs_0.radial_coeffs_0.z;
    float k4_0 = dist_coeffs_0.radial_coeffs_0.w;
    float p1_0 = dist_coeffs_0.tangential_coeffs_0.x;
    float p2_0 = dist_coeffs_0.tangential_coeffs_0.y;
    float sx1_0 = dist_coeffs_0.thin_prism_coeffs_0.x;
    float sy1_0 = dist_coeffs_0.thin_prism_coeffs_0.y;
    float u_0 = _S40.x;
    float v_0 = _S40.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float _S41 = 2.0f * p1_0;
    float _S42 = 2.0f * p2_0;
    float2  _S43 = _S40 * make_float2 (1.0f + r2_0 * (k1_0 + r2_0 * (k2_0 + r2_0 * (k3_0 + r2_0 * k4_0)))) + make_float2 (_S41 * u_0 * v_0 + p2_0 * (r2_0 + 2.0f * u_0 * u_0) + sx1_0 * r2_0, _S42 * u_0 * v_0 + p1_0 * (r2_0 + 2.0f * v_0 * v_0) + sy1_0 * r2_0);
    *mean2d_2 = make_float2 (fx_2 * _S43.x + cx_2, fy_2 * _S43.y + cy_2);
    Matrix<float, 2, 3>  J_2;
    float2  _S44 = make_float2 (0.0f);
    float2  seed_0 = _S44;
    *&((&seed_0)->x) = 1.0f;
    float2  _S45 = seed_0;
    float _S46 = s_primal_ctx_atan2_0(r_6, _S39);
    bool _S47 = _S46 < 0.00100000004749745f;
    float _S48;
    float _S49;
    float _S50;
    if(_S47)
    {
        float _S51 = 1.0f - _S46 * _S46 / 3.0f;
        float _S52 = _S39 * _S39;
        k_0 = _S51 / _S39;
        _S48 = 0.0f;
        _S49 = _S52;
        _S50 = _S51;
    }
    else
    {
        float _S53 = r_6 * r_6;
        k_0 = _S46 / r_6;
        _S48 = _S53;
        _S49 = 0.0f;
        _S50 = 0.0f;
    }
    float2  _S54 = make_float2 (k_0);
    float2  _S55 = _S38 * make_float2 (k_0);
    float u_1 = _S55.x;
    float v_1 = _S55.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float _S56 = k3_0 + r2_1 * k4_0;
    float _S57 = k2_0 + r2_1 * _S56;
    float _S58 = k1_0 + r2_1 * _S57;
    float _S59 = fy_2 * _S45.y;
    float _S60 = fx_2 * _S45.x;
    float2  _S61 = make_float2 (_S60, _S59);
    float2  _S62 = _S55 * _S61;
    float _S63 = p1_0 * _S59;
    float _S64 = p2_0 * _S60;
    float _S65 = _S62.x + _S62.y;
    float _S66 = r2_1 * _S65;
    float _S67 = r2_1 * _S66;
    float _S68 = sy1_0 * _S59 + _S63 + sx1_0 * _S60 + _S64 + _S58 * _S65 + _S57 * _S66 + _S56 * _S67 + k4_0 * (r2_1 * _S67);
    float _S69 = v_1 * _S68;
    float _S70 = u_1 * _S68;
    float2  _S71 = make_float2 (1.0f + r2_1 * _S58) * _S61 + make_float2 (_S42 * (v_1 * _S59) + 2.0f * u_1 * _S64 + 2.0f * (u_1 * _S64) + _S41 * (v_1 * _S60) + _S70 + _S70, 2.0f * v_1 * _S63 + 2.0f * (v_1 * _S63) + _S42 * u_1 * _S59 + _S41 * u_1 * _S60 + _S69 + _S69);
    float2  _S72 = _S38 * _S71;
    float2  _S73 = _S54 * _S71;
    float _S74 = _S72.x + _S72.y;
    if(_S47)
    {
        float _S75 = _S74 / _S49;
        float _S76 = _S50 * - _S75;
        float _S77 = _S46 * (0.3333333432674408f * - (_S39 * _S75));
        k_0 = _S77 + _S77;
        _S48 = _S76;
        _S49 = 0.0f;
    }
    else
    {
        float _S78 = _S74 / _S48;
        float _S79 = _S46 * - _S78;
        k_0 = r_6 * _S78;
        _S48 = 0.0f;
        _S49 = _S79;
    }
    DiffPair_float_0 _S80;
    (&_S80)->primal_0 = r_6;
    (&_S80)->differential_0 = 0.0f;
    DiffPair_float_0 _S81;
    (&_S81)->primal_0 = _S39;
    (&_S81)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S80, &_S81, k_0);
    float _S82 = _S81.differential_0 + _S48;
    float _S83 = _S80.differential_0 + _S49;
    float2  _S84 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S85;
    (&_S85)->primal_0 = _S38;
    (&_S85)->differential_0 = _S84;
    s_bwd_length_impl_0(&_S85, _S83);
    float2  _S86 = _S85.differential_0 + _S73;
    float3  _S87 = make_float3 (_S86.x, _S86.y, _S82);
    J_2[int(0)] = _S87;
    float2  seed_1 = _S44;
    *&((&seed_1)->y) = 1.0f;
    float2  _S88 = seed_1;
    if(_S47)
    {
        float _S89 = 1.0f - _S46 * _S46 / 3.0f;
        float _S90 = _S39 * _S39;
        k_0 = _S89 / _S39;
        _S48 = 0.0f;
        _S49 = _S90;
        _S50 = _S89;
    }
    else
    {
        float _S91 = r_6 * r_6;
        k_0 = _S46 / r_6;
        _S48 = _S91;
        _S49 = 0.0f;
        _S50 = 0.0f;
    }
    float2  _S92 = make_float2 (k_0);
    float2  _S93 = _S38 * make_float2 (k_0);
    float u_2 = _S93.x;
    float v_2 = _S93.y;
    float r2_2 = u_2 * u_2 + v_2 * v_2;
    float _S94 = k3_0 + r2_2 * k4_0;
    float _S95 = k2_0 + r2_2 * _S94;
    float _S96 = k1_0 + r2_2 * _S95;
    float _S97 = fy_2 * _S88.y;
    float _S98 = fx_2 * _S88.x;
    float2  _S99 = make_float2 (_S98, _S97);
    float2  _S100 = _S93 * _S99;
    float _S101 = p1_0 * _S97;
    float _S102 = p2_0 * _S98;
    float _S103 = _S100.x + _S100.y;
    float _S104 = r2_2 * _S103;
    float _S105 = r2_2 * _S104;
    float _S106 = sy1_0 * _S97 + _S101 + sx1_0 * _S98 + _S102 + _S96 * _S103 + _S95 * _S104 + _S94 * _S105 + k4_0 * (r2_2 * _S105);
    float _S107 = v_2 * _S106;
    float _S108 = u_2 * _S106;
    float2  _S109 = make_float2 (1.0f + r2_2 * _S96) * _S99 + make_float2 (_S42 * (v_2 * _S97) + 2.0f * u_2 * _S102 + 2.0f * (u_2 * _S102) + _S41 * (v_2 * _S98) + _S108 + _S108, 2.0f * v_2 * _S101 + 2.0f * (v_2 * _S101) + _S42 * u_2 * _S97 + _S41 * u_2 * _S98 + _S107 + _S107);
    float2  _S110 = _S38 * _S109;
    float2  _S111 = _S92 * _S109;
    float _S112 = _S110.x + _S110.y;
    if(_S47)
    {
        float _S113 = _S112 / _S49;
        float _S114 = _S50 * - _S113;
        float _S115 = _S46 * (0.3333333432674408f * - (_S39 * _S113));
        k_0 = _S115 + _S115;
        _S48 = _S114;
        _S49 = 0.0f;
    }
    else
    {
        float _S116 = _S112 / _S48;
        float _S117 = _S46 * - _S116;
        k_0 = r_6 * _S116;
        _S48 = 0.0f;
        _S49 = _S117;
    }
    DiffPair_float_0 _S118;
    (&_S118)->primal_0 = r_6;
    (&_S118)->differential_0 = 0.0f;
    DiffPair_float_0 _S119;
    (&_S119)->primal_0 = _S39;
    (&_S119)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S118, &_S119, k_0);
    float _S120 = _S119.differential_0 + _S48;
    float _S121 = _S118.differential_0 + _S49;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S122;
    (&_S122)->primal_0 = _S38;
    (&_S122)->differential_0 = _S84;
    s_bwd_length_impl_0(&_S122, _S121);
    float2  _S123 = _S122.differential_0 + _S111;
    float3  _S124 = make_float3 (_S123.x, _S123.y, _S120);
    J_2[int(1)] = _S124;
    *cov2d_2 = mul_6(mul_5(J_2, cov3d_2), transpose_1(J_2));
    return;
}

inline __device__ void ortho_proj_3dgs(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float fx_3, float fy_3, float cx_3, float cy_3, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    Matrix<float, 2, 3>  J_3 = makeMatrix<float, 2, 3> (fx_3, 0.0f, 0.0f, 0.0f, fy_3, 0.0f);
    *cov2d_3 = mul_6(mul_5(J_3, cov3d_3), transpose_1(J_3));
    *mean2d_3 = make_float2 (fx_3 * mean3d_3.x + cx_3, fy_3 * mean3d_3.y + cy_3);
    return;
}

inline __device__ float add_blur(float eps2d_0, Matrix<float, 2, 2>  * covar_1, float * compensation_0)
{
    float det_orig_0 = *&((covar_1->rows + (int(0)))->x) * *&((covar_1->rows + (int(1)))->y) - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *&((covar_1->rows + (int(0)))->x) = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    float _S125 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S125;
    float det_blur_0 = *&((covar_1->rows + (int(0)))->x) * _S125 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_8, float dOut_11)
{
    float _S126 = (F32_exp(((*dpx_8).primal_0))) * dOut_11;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S126;
    return;
}

inline __device__ float3  exp_0(float3  x_13)
{
    float3  result_13;
    int i_6 = int(0);
    for(;;)
    {
        if(i_6 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_13, i_6) = (F32_exp((_slang_vector_get_element(x_13, i_6))));
        i_6 = i_6 + int(1);
    }
    return result_13;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float3  dOut_12)
{
    float3  _S127 = exp_0((*dpx_9).primal_0) * dOut_12;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S127;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_10, float dOut_13)
{
    float _S128 = 1.0f / (*dpx_10).primal_0 * dOut_13;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S128;
    return;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, Matrix<float, 3, 3>  R_2, float3  t_1, float fx_4, float fy_4, float cx_4, float cy_4, float4  radial_coeffs_3, float2  tangential_coeffs_3, float2  thin_prism_coeffs_3, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float * depth_0, float2  * mean2d_4, float3  * conic_0, float * opacity_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_2, mean_0) + t_1;
        float _S129 = mean_c_0.z;
        bool _S130;
        if(_S129 < near_plane_0)
        {
            _S130 = true;
        }
        else
        {
            _S130 = _S129 > far_plane_0;
        }
        if(_S130)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S131 = exp_0(scale_2);
        float x_14 = quat_3.y;
        float inv_norm_3 = (F32_rsqrt((x_14 * x_14 + quat_3.z * quat_3.z + quat_3.w * quat_3.w + quat_3.x * quat_3.x)));
        float x_15 = quat_3.y * inv_norm_3;
        float y_5 = quat_3.z * inv_norm_3;
        float z_3 = quat_3.w * inv_norm_3;
        float w_3 = quat_3.x * inv_norm_3;
        float x2_3 = x_15 * x_15;
        float y2_3 = y_5 * y_5;
        float z2_3 = z_3 * z_3;
        float xy_3 = x_15 * y_5;
        float xz_3 = x_15 * z_3;
        float yz_3 = y_5 * z_3;
        float wx_3 = w_3 * x_15;
        float wy_3 = w_3 * y_5;
        float wz_3 = w_3 * z_3;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S131.x, 0.0f, 0.0f, 0.0f, _S131.y, 0.0f, 0.0f, 0.0f, _S131.z));
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_2, mul_4(M_2, transpose_0(M_2))), transpose_0(R_2));
        Matrix<float, 2, 2>  covar2d_0;
        float _S132 = float(image_width_0);
        float _S133 = float(image_height_0);
        float _S134 = 0.30000001192092896f * (0.5f * _S132 / fx_4);
        float _S135 = 0.30000001192092896f * (0.5f * _S133 / fy_4);
        float rz_2 = 1.0f / mean_c_0.z;
        float rz2_2 = rz_2 * rz_2;
        Matrix<float, 2, 3>  J_4 = makeMatrix<float, 2, 3> (fx_4 * rz_2, 0.0f, - fx_4 * (mean_c_0.z * (F32_min(((_S132 - cx_4) / fx_4 + _S134), ((F32_max((- (cx_4 / fx_4 + _S134)), (mean_c_0.x * rz_2))))))) * rz2_2, 0.0f, fy_4 * rz_2, - fy_4 * (mean_c_0.z * (F32_min(((_S133 - cy_4) / fy_4 + _S135), ((F32_max((- (cy_4 / fy_4 + _S135)), (mean_c_0.y * rz_2))))))) * rz2_2);
        covar2d_0 = mul_6(mul_5(J_4, covar_c_0), transpose_1(J_4));
        *mean2d_4 = make_float2 (fx_4 * mean_c_0.x * rz_2 + cx_4, fy_4 * mean_c_0.y * rz_2 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        *&(((&covar2d_0)->rows + (int(0)))->x) = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        float _S136 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S136;
        float det_blur_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * _S136 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S137 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
        *opacity_0 = 1.0f / (1.0f + (F32_exp((- in_opacity_0))));
        if(antialiased_0)
        {
            *opacity_0 = *opacity_0 * compensation_1;
        }
        if((*opacity_0) < 0.00392156885936856f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float extend_0 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_0 / 0.00392156885936856f)))))))));
        float radius_x_0 = extend_0 * (F32_sqrt((covar2d_0[int(0)].x)));
        float radius_y_0 = extend_0 * (F32_sqrt((covar2d_0[int(1)].y)));
        float xmin_0 = (F32_floor(((*mean2d_4).x - radius_x_0)));
        float xmax_0 = (F32_ceil(((*mean2d_4).x + radius_x_0)));
        float ymin_0 = (F32_floor(((*mean2d_4).y - radius_y_0)));
        float ymax_0 = (F32_ceil(((*mean2d_4).y + radius_y_0)));
        if(xmax_0 <= 0.0f)
        {
            _S130 = true;
        }
        else
        {
            _S130 = xmin_0 >= _S132;
        }
        if(_S130)
        {
            _S130 = true;
        }
        else
        {
            _S130 = ymax_0 <= 0.0f;
        }
        if(_S130)
        {
            _S130 = true;
        }
        else
        {
            _S130 = ymin_0 >= _S133;
        }
        if(_S130)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = _S129;
        *conic_0 = make_float3 (_S137.rows[int(0)].x, _S137.rows[int(0)].y, _S137.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, Matrix<float, 3, 3>  R_3, float3  t_2, float fx_5, float fy_5, float cx_5, float cy_5, float4  radial_coeffs_4, float2  tangential_coeffs_4, float2  thin_prism_coeffs_4, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float * depth_1, float2  * mean2d_5, float3  * conic_1, float * opacity_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_3, mean_1) + t_2;
        float _S138 = mean_c_1.z;
        bool _S139;
        if(_S138 < near_plane_1)
        {
            _S139 = true;
        }
        else
        {
            _S139 = _S138 > far_plane_1;
        }
        if(_S139)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S140 = exp_0(scale_3);
        float x_16 = quat_4.y;
        float inv_norm_4 = (F32_rsqrt((x_16 * x_16 + quat_4.z * quat_4.z + quat_4.w * quat_4.w + quat_4.x * quat_4.x)));
        float x_17 = quat_4.y * inv_norm_4;
        float y_6 = quat_4.z * inv_norm_4;
        float z_4 = quat_4.w * inv_norm_4;
        float w_4 = quat_4.x * inv_norm_4;
        float x2_4 = x_17 * x_17;
        float y2_4 = y_6 * y_6;
        float z2_4 = z_4 * z_4;
        float xy_4 = x_17 * y_6;
        float xz_4 = x_17 * z_4;
        float yz_4 = y_6 * z_4;
        float wx_4 = w_4 * x_17;
        float wy_4 = w_4 * y_6;
        float wz_4 = w_4 * z_4;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S140.x, 0.0f, 0.0f, 0.0f, _S140.y, 0.0f, 0.0f, 0.0f, _S140.z));
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_3, mul_4(M_3, transpose_0(M_3))), transpose_0(R_3));
        Matrix<float, 2, 2>  covar2d_1;
        fisheye_proj_3dgs(mean_c_1, covar_c_1, fx_5, fy_5, cx_5, cy_5, radial_coeffs_4, tangential_coeffs_4, thin_prism_coeffs_4, &covar2d_1, mean2d_5);
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        *&(((&covar2d_1)->rows + (int(0)))->x) = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        float _S141 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S141;
        float det_blur_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * _S141 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S142 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
        *opacity_1 = 1.0f / (1.0f + (F32_exp((- in_opacity_1))));
        if(antialiased_1)
        {
            *opacity_1 = *opacity_1 * compensation_2;
        }
        if((*opacity_1) < 0.00392156885936856f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float extend_1 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_1 / 0.00392156885936856f)))))))));
        float radius_x_1 = extend_1 * (F32_sqrt((covar2d_1[int(0)].x)));
        float radius_y_1 = extend_1 * (F32_sqrt((covar2d_1[int(1)].y)));
        float xmin_1 = (F32_floor(((*mean2d_5).x - radius_x_1)));
        float xmax_1 = (F32_ceil(((*mean2d_5).x + radius_x_1)));
        float ymin_1 = (F32_floor(((*mean2d_5).y - radius_y_1)));
        float ymax_1 = (F32_ceil(((*mean2d_5).y + radius_y_1)));
        if(xmax_1 <= 0.0f)
        {
            _S139 = true;
        }
        else
        {
            _S139 = xmin_1 >= float(image_width_1);
        }
        if(_S139)
        {
            _S139 = true;
        }
        else
        {
            _S139 = ymax_1 <= 0.0f;
        }
        if(_S139)
        {
            _S139 = true;
        }
        else
        {
            _S139 = ymin_1 >= float(image_height_1);
        }
        if(_S139)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = _S138;
        *conic_1 = make_float3 (_S142.rows[int(0)].x, _S142.rows[int(0)].y, _S142.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_6, float fy_6, float cx_6, float cy_6, float4  radial_coeffs_5, float2  tangential_coeffs_5, float2  thin_prism_coeffs_5, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float * depth_2, float2  * mean2d_6, float3  * conic_2, float * opacity_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_4, mean_2) + t_3;
        float _S143 = mean_c_2.z;
        bool _S144;
        if(_S143 < near_plane_2)
        {
            _S144 = true;
        }
        else
        {
            _S144 = _S143 > far_plane_2;
        }
        if(_S144)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S145 = exp_0(scale_4);
        float x_18 = quat_5.y;
        float inv_norm_5 = (F32_rsqrt((x_18 * x_18 + quat_5.z * quat_5.z + quat_5.w * quat_5.w + quat_5.x * quat_5.x)));
        float x_19 = quat_5.y * inv_norm_5;
        float y_7 = quat_5.z * inv_norm_5;
        float z_5 = quat_5.w * inv_norm_5;
        float w_5 = quat_5.x * inv_norm_5;
        float x2_5 = x_19 * x_19;
        float y2_5 = y_7 * y_7;
        float z2_5 = z_5 * z_5;
        float xy_5 = x_19 * y_7;
        float xz_5 = x_19 * z_5;
        float yz_5 = y_7 * z_5;
        float wx_5 = w_5 * x_19;
        float wy_5 = w_5 * y_7;
        float wz_5 = w_5 * z_5;
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_5), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_5), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S145.x, 0.0f, 0.0f, 0.0f, _S145.y, 0.0f, 0.0f, 0.0f, _S145.z));
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_4, mul_4(M_4, transpose_0(M_4))), transpose_0(R_4));
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_5, covar_c_2), transpose_1(J_5));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        *&(((&covar2d_2)->rows + (int(0)))->x) = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        float _S146 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S146;
        float det_blur_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * _S146 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S147 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
        *opacity_2 = 1.0f / (1.0f + (F32_exp((- in_opacity_2))));
        if(antialiased_2)
        {
            *opacity_2 = *opacity_2 * compensation_3;
        }
        if((*opacity_2) < 0.00392156885936856f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float extend_2 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_2 / 0.00392156885936856f)))))))));
        float radius_x_2 = extend_2 * (F32_sqrt((covar2d_2[int(0)].x)));
        float radius_y_2 = extend_2 * (F32_sqrt((covar2d_2[int(1)].y)));
        float xmin_2 = (F32_floor(((*mean2d_6).x - radius_x_2)));
        float xmax_2 = (F32_ceil(((*mean2d_6).x + radius_x_2)));
        float ymin_2 = (F32_floor(((*mean2d_6).y - radius_y_2)));
        float ymax_2 = (F32_ceil(((*mean2d_6).y + radius_y_2)));
        if(xmax_2 <= 0.0f)
        {
            _S144 = true;
        }
        else
        {
            _S144 = xmin_2 >= float(image_width_2);
        }
        if(_S144)
        {
            _S144 = true;
        }
        else
        {
            _S144 = ymax_2 <= 0.0f;
        }
        if(_S144)
        {
            _S144 = true;
        }
        else
        {
            _S144 = ymin_2 >= float(image_height_2);
        }
        if(_S144)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = _S143;
        *conic_2 = make_float3 (_S147.rows[int(0)].x, _S147.rows[int(0)].y, _S147.rows[int(1)].y);
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_7, float fy_7, float cx_7, float cy_7, float4  radial_coeffs_6, float2  tangential_coeffs_6, float2  thin_prism_coeffs_6, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float * depth_3, float2  * mean2d_7, float3  * conic_3, float * opacity_3)
{
    float3  mean_c_3 = mul_0(R_5, mean_3) + t_4;
    float3  _S148 = exp_0(scale_5);
    float x_20 = quat_6.y;
    float inv_norm_6 = (F32_rsqrt((x_20 * x_20 + quat_6.z * quat_6.z + quat_6.w * quat_6.w + quat_6.x * quat_6.x)));
    float x_21 = quat_6.y * inv_norm_6;
    float y_8 = quat_6.z * inv_norm_6;
    float z_6 = quat_6.w * inv_norm_6;
    float w_6 = quat_6.x * inv_norm_6;
    float x2_6 = x_21 * x_21;
    float y2_6 = y_8 * y_8;
    float z2_6 = z_6 * z_6;
    float xy_6 = x_21 * y_8;
    float xz_6 = x_21 * z_6;
    float yz_6 = y_8 * z_6;
    float wx_6 = w_6 * x_21;
    float wy_6 = w_6 * y_8;
    float wz_6 = w_6 * z_6;
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_6), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_6), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S148.x, 0.0f, 0.0f, 0.0f, _S148.y, 0.0f, 0.0f, 0.0f, _S148.z));
    float _S149 = float(image_width_3);
    float _S150 = float(image_height_3);
    float _S151 = 0.30000001192092896f * (0.5f * _S149 / fx_7);
    float _S152 = 0.30000001192092896f * (0.5f * _S150 / fy_7);
    float rz_3 = 1.0f / mean_c_3.z;
    float rz2_3 = rz_3 * rz_3;
    Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, - fx_7 * (mean_c_3.z * (F32_min(((_S149 - cx_7) / fx_7 + _S151), ((F32_max((- (cx_7 / fx_7 + _S151)), (mean_c_3.x * rz_3))))))) * rz2_3, 0.0f, fy_7 * rz_3, - fy_7 * (mean_c_3.z * (F32_min(((_S150 - cy_7) / fy_7 + _S152), ((F32_max((- (cy_7 / fy_7 + _S152)), (mean_c_3.y * rz_3))))))) * rz2_3);
    Matrix<float, 2, 2>  covar2d_3 = mul_6(mul_5(J_6, mul_4(mul_4(R_5, mul_4(M_5, transpose_0(M_5))), transpose_0(R_5))), transpose_1(J_6));
    *mean2d_7 = make_float2 (fx_7 * mean_c_3.x * rz_3 + cx_7, fy_7 * mean_c_3.y * rz_3 + cy_7);
    float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
    *&(((&covar2d_3)->rows + (int(0)))->x) = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
    float _S153 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_3)->rows + (int(1)))->y) = _S153;
    float compensation_4 = (F32_sqrt(((F32_max((0.0f), (det_orig_4 / (*&(((&covar2d_3)->rows + (int(0)))->x) * _S153 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_3.rows[int(0)].x * covar2d_3.rows[int(1)].y - covar2d_3.rows[int(0)].y * covar2d_3.rows[int(1)].x);
    Matrix<float, 2, 2>  _S154 = makeMatrix<float, 2, 2> (covar2d_3.rows[int(1)].y * invdet_4, - covar2d_3.rows[int(0)].y * invdet_4, - covar2d_3.rows[int(1)].x * invdet_4, covar2d_3.rows[int(0)].x * invdet_4);
    *opacity_3 = 1.0f / (1.0f + (F32_exp((- in_opacity_3))));
    if(antialiased_3)
    {
        *opacity_3 = *opacity_3 * compensation_4;
    }
    float extend_3 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_3 / 0.00392156885936856f)))))))));
    float radius_x_3 = extend_3 * (F32_sqrt((covar2d_3[int(0)].x)));
    float radius_y_3 = extend_3 * (F32_sqrt((covar2d_3[int(1)].y)));
    *aabb_xyxy_3 = make_int4 (int((F32_floor(((*mean2d_7).x - radius_x_3)))), int((F32_floor(((*mean2d_7).y - radius_y_3)))), int((F32_ceil(((*mean2d_7).x + radius_x_3)))), int((F32_ceil(((*mean2d_7).y + radius_y_3)))));
    *depth_3 = mean_c_3.z;
    *conic_3 = make_float3 (_S154.rows[int(0)].x, _S154.rows[int(0)].y, _S154.rows[int(1)].y);
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_8, float fy_8, float cx_8, float cy_8, float4  radial_coeffs_7, float2  tangential_coeffs_7, float2  thin_prism_coeffs_7, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float * depth_4, float2  * mean2d_8, float3  * conic_4, float * opacity_4)
{
    float3  mean_c_4 = mul_0(R_6, mean_4) + t_5;
    float3  _S155 = exp_0(scale_6);
    float x_22 = quat_7.y;
    float inv_norm_7 = (F32_rsqrt((x_22 * x_22 + quat_7.z * quat_7.z + quat_7.w * quat_7.w + quat_7.x * quat_7.x)));
    float x_23 = quat_7.y * inv_norm_7;
    float y_9 = quat_7.z * inv_norm_7;
    float z_7 = quat_7.w * inv_norm_7;
    float w_7 = quat_7.x * inv_norm_7;
    float x2_7 = x_23 * x_23;
    float y2_7 = y_9 * y_9;
    float z2_7 = z_7 * z_7;
    float xy_7 = x_23 * y_9;
    float xz_7 = x_23 * z_7;
    float yz_7 = y_9 * z_7;
    float wx_7 = w_7 * x_23;
    float wy_7 = w_7 * y_9;
    float wz_7 = w_7 * z_7;
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_7), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_7), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S155.x, 0.0f, 0.0f, 0.0f, _S155.y, 0.0f, 0.0f, 0.0f, _S155.z));
    Matrix<float, 2, 2>  covar2d_4;
    fisheye_proj_3dgs(mean_c_4, mul_4(mul_4(R_6, mul_4(M_6, transpose_0(M_6))), transpose_0(R_6)), fx_8, fy_8, cx_8, cy_8, radial_coeffs_7, tangential_coeffs_7, thin_prism_coeffs_7, &covar2d_4, mean2d_8);
    float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
    *&(((&covar2d_4)->rows + (int(0)))->x) = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
    float _S156 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_4)->rows + (int(1)))->y) = _S156;
    float compensation_5 = (F32_sqrt(((F32_max((0.0f), (det_orig_5 / (*&(((&covar2d_4)->rows + (int(0)))->x) * _S156 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_4.rows[int(0)].x * covar2d_4.rows[int(1)].y - covar2d_4.rows[int(0)].y * covar2d_4.rows[int(1)].x);
    Matrix<float, 2, 2>  _S157 = makeMatrix<float, 2, 2> (covar2d_4.rows[int(1)].y * invdet_5, - covar2d_4.rows[int(0)].y * invdet_5, - covar2d_4.rows[int(1)].x * invdet_5, covar2d_4.rows[int(0)].x * invdet_5);
    *opacity_4 = 1.0f / (1.0f + (F32_exp((- in_opacity_4))));
    if(antialiased_4)
    {
        *opacity_4 = *opacity_4 * compensation_5;
    }
    float extend_4 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_4 / 0.00392156885936856f)))))))));
    float radius_x_4 = extend_4 * (F32_sqrt((covar2d_4[int(0)].x)));
    float radius_y_4 = extend_4 * (F32_sqrt((covar2d_4[int(1)].y)));
    *aabb_xyxy_4 = make_int4 (int((F32_floor(((*mean2d_8).x - radius_x_4)))), int((F32_floor(((*mean2d_8).y - radius_y_4)))), int((F32_ceil(((*mean2d_8).x + radius_x_4)))), int((F32_ceil(((*mean2d_8).y + radius_y_4)))));
    *depth_4 = mean_c_4.z;
    *conic_4 = make_float3 (_S157.rows[int(0)].x, _S157.rows[int(0)].y, _S157.rows[int(1)].y);
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_9, float fy_9, float cx_9, float cy_9, float4  radial_coeffs_8, float2  tangential_coeffs_8, float2  thin_prism_coeffs_8, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float * depth_5, float2  * mean2d_9, float3  * conic_5, float * opacity_5)
{
    float3  mean_c_5 = mul_0(R_7, mean_5) + t_6;
    float3  _S158 = exp_0(scale_7);
    float x_24 = quat_8.y;
    float inv_norm_8 = (F32_rsqrt((x_24 * x_24 + quat_8.z * quat_8.z + quat_8.w * quat_8.w + quat_8.x * quat_8.x)));
    float x_25 = quat_8.y * inv_norm_8;
    float y_10 = quat_8.z * inv_norm_8;
    float z_8 = quat_8.w * inv_norm_8;
    float w_8 = quat_8.x * inv_norm_8;
    float x2_8 = x_25 * x_25;
    float y2_8 = y_10 * y_10;
    float z2_8 = z_8 * z_8;
    float xy_8 = x_25 * y_10;
    float xz_8 = x_25 * z_8;
    float yz_8 = y_10 * z_8;
    float wx_8 = w_8 * x_25;
    float wy_8 = w_8 * y_10;
    float wz_8 = w_8 * z_8;
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_8), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_8), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S158.x, 0.0f, 0.0f, 0.0f, _S158.y, 0.0f, 0.0f, 0.0f, _S158.z));
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_7, mul_4(mul_4(R_7, mul_4(M_7, transpose_0(M_7))), transpose_0(R_7))), transpose_1(J_7));
    *mean2d_9 = make_float2 (fx_9 * mean_c_5.x + cx_9, fy_9 * mean_c_5.y + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    *&(((&covar2d_5)->rows + (int(0)))->x) = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    float _S159 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S159;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (*&(((&covar2d_5)->rows + (int(0)))->x) * _S159 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S160 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_6, - covar2d_5.rows[int(0)].y * invdet_6, - covar2d_5.rows[int(1)].x * invdet_6, covar2d_5.rows[int(0)].x * invdet_6);
    *opacity_5 = 1.0f / (1.0f + (F32_exp((- in_opacity_5))));
    if(antialiased_5)
    {
        *opacity_5 = *opacity_5 * compensation_6;
    }
    float extend_5 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_5 / 0.00392156885936856f)))))))));
    float radius_x_5 = extend_5 * (F32_sqrt((covar2d_5[int(0)].x)));
    float radius_y_5 = extend_5 * (F32_sqrt((covar2d_5[int(1)].y)));
    *aabb_xyxy_5 = make_int4 (int((F32_floor(((*mean2d_9).x - radius_x_5)))), int((F32_floor(((*mean2d_9).y - radius_y_5)))), int((F32_ceil(((*mean2d_9).x + radius_x_5)))), int((F32_ceil(((*mean2d_9).y + radius_y_5)))));
    *depth_5 = mean_c_5.z;
    *conic_5 = make_float3 (_S160.rows[int(0)].x, _S160.rows[int(0)].y, _S160.rows[int(1)].y);
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S161, float3  _S162)
{
    return mul_0(_S161, _S162);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S163)
{
    return exp_0(_S163);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S164)
{
    return (F32_rsqrt((_S164)));
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S165, Matrix<float, 3, 3>  _S166)
{
    return mul_4(_S165, _S166);
}

inline __device__ float s_primal_ctx_max_0(float _S167, float _S168)
{
    return (F32_max((_S167), (_S168)));
}

inline __device__ float s_primal_ctx_min_0(float _S169, float _S170)
{
    return (F32_min((_S169), (_S170)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_2(Matrix<float, 2, 3>  _S171, Matrix<float, 3, 3>  _S172)
{
    return mul_5(_S171, _S172);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S173, Matrix<float, 3, 2>  _S174)
{
    return mul_6(_S173, _S174);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S175)
{
    return (F32_sqrt((_S175)));
}

inline __device__ float s_primal_ctx_exp_1(float _S176)
{
    return (F32_exp((_S176)));
}

inline __device__ float s_primal_ctx_log_0(float _S177)
{
    return (F32_log((_S177)));
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S178, DiffPair_float_0 * _S179, float _S180)
{
    _d_min_0(_S178, _S179, _S180);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S181, float _S182)
{
    _d_log_0(_S181, _S182);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S183, float _S184)
{
    _d_exp_0(_S183, _S184);
    return;
}

inline __device__ void s_bwd_prop_max_0(DiffPair_float_0 * _S185, DiffPair_float_0 * _S186, float _S187)
{
    _d_max_0(_S185, _S186, _S187);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S188, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S189, Matrix<float, 2, 2>  _S190)
{
    mul_3(_S188, _S189, _S190);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S191, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S192, Matrix<float, 2, 3>  _S193)
{
    mul_2(_S191, _S192, _S193);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S194, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S195, Matrix<float, 3, 3>  _S196)
{
    mul_1(_S194, _S195, _S196);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S197, float _S198)
{
    _d_rsqrt_0(_S197, _S198);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S199, float3  _S200)
{
    _d_exp_vector_0(_S199, _S200);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S201, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S202, float3  _S203)
{
    _d_mul_0(_S201, _S202, _S203);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_10, float fy_10, float cx_10, float cy_10, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, uint image_width_6, uint image_height_6, float v_depth_0, float2  v_mean2d_0, float3  v_conic_0, float v_opacity_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_6 = s_primal_ctx_mul_0(R_8, mean_6) + t_7;
    float3  _S204 = s_primal_ctx_exp_0(scale_8);
    float _S205 = quat_9.y;
    float _S206 = _S205 * _S205 + quat_9.z * quat_9.z + quat_9.w * quat_9.w + quat_9.x * quat_9.x;
    float _S207 = s_primal_ctx_rsqrt_0(_S206);
    float x_26 = quat_9.y * _S207;
    float y_11 = quat_9.z * _S207;
    float z_9 = quat_9.w * _S207;
    float w_9 = quat_9.x * _S207;
    float x2_9 = x_26 * x_26;
    float y2_9 = y_11 * y_11;
    float z2_9 = z_9 * z_9;
    float xy_9 = x_26 * y_11;
    float xz_9 = x_26 * z_9;
    float yz_9 = y_11 * z_9;
    float wx_9 = w_9 * x_26;
    float wy_9 = w_9 * y_11;
    float wz_9 = w_9 * z_9;
    Matrix<float, 3, 3>  _S208 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_9), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_9), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S204.x, 0.0f, 0.0f, 0.0f, _S204.y, 0.0f, 0.0f, 0.0f, _S204.z);
    Matrix<float, 3, 3>  _S209 = s_primal_ctx_mul_1(_S208, S_0);
    Matrix<float, 3, 3>  _S210 = transpose_0(_S209);
    Matrix<float, 3, 3>  _S211 = s_primal_ctx_mul_1(_S209, _S210);
    Matrix<float, 3, 3>  _S212 = s_primal_ctx_mul_1(R_8, _S211);
    Matrix<float, 3, 3>  _S213 = transpose_0(R_8);
    Matrix<float, 3, 3>  _S214 = s_primal_ctx_mul_1(_S212, _S213);
    float _S215 = float(image_width_6);
    float _S216 = float(image_height_6);
    float _S217 = 0.30000001192092896f * (0.5f * _S215 / fx_10);
    float lim_x_pos_0 = (_S215 - cx_10) / fx_10 + _S217;
    float _S218 = 0.30000001192092896f * (0.5f * _S216 / fy_10);
    float lim_y_pos_0 = (_S216 - cy_10) / fy_10 + _S218;
    float rz_4 = 1.0f / mean_c_6.z;
    float _S219 = mean_c_6.z * mean_c_6.z;
    float rz2_4 = rz_4 * rz_4;
    float _S220 = - (cx_10 / fx_10 + _S217);
    float _S221 = mean_c_6.x * rz_4;
    float _S222 = s_primal_ctx_max_0(_S220, _S221);
    float _S223 = s_primal_ctx_min_0(lim_x_pos_0, _S222);
    float _S224 = - (cy_10 / fy_10 + _S218);
    float _S225 = mean_c_6.y * rz_4;
    float _S226 = s_primal_ctx_max_0(_S224, _S225);
    float _S227 = s_primal_ctx_min_0(lim_y_pos_0, _S226);
    float _S228 = - fx_10;
    float _S229 = _S228 * (mean_c_6.z * _S223);
    float _S230 = - fy_10;
    float _S231 = _S230 * (mean_c_6.z * _S227);
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_10 * rz_4, 0.0f, _S229 * rz2_4, 0.0f, fy_10 * rz_4, _S231 * rz2_4);
    Matrix<float, 2, 3>  _S232 = s_primal_ctx_mul_2(J_8, _S214);
    Matrix<float, 3, 2>  _S233 = transpose_1(J_8);
    Matrix<float, 2, 2>  _S234 = s_primal_ctx_mul_3(_S232, _S233);
    float _S235 = fx_10 * mean_c_6.x;
    float _S236 = fy_10 * mean_c_6.y;
    float _S237 = _S234.rows[int(0)].y * _S234.rows[int(1)].x;
    float det_orig_7 = _S234.rows[int(0)].x * _S234.rows[int(1)].y - _S237;
    float _S238 = _S234.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S239 = _S234;
    *&(((&_S239)->rows + (int(0)))->x) = _S238;
    float _S240 = _S234.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S239)->rows + (int(1)))->y) = _S240;
    Matrix<float, 2, 2>  _S241 = _S239;
    Matrix<float, 2, 2>  _S242 = _S239;
    float det_blur_4 = _S238 * _S240 - _S237;
    float _S243 = det_orig_7 / det_blur_4;
    float _S244 = det_blur_4 * det_blur_4;
    float _S245 = s_primal_ctx_max_0(0.0f, _S243);
    float _S246 = s_primal_ctx_sqrt_0(_S245);
    float invdet_7 = 1.0f / det_blur_4;
    float _S247 = - _S234.rows[int(0)].y;
    float _S248 = - _S234.rows[int(1)].x;
    float _S249 = - in_opacity_6;
    float _S250 = 1.0f + s_primal_ctx_exp_1(_S249);
    float _S251 = 1.0f / _S250;
    float _S252 = _S250 * _S250;
    float _S253;
    if(antialiased_6)
    {
        _S253 = _S251 * _S246;
    }
    else
    {
        _S253 = _S251;
    }
    float _S254 = _S253 / 0.00392156885936856f;
    float _S255 = 2.0f * s_primal_ctx_log_0(_S254);
    float _S256 = s_primal_ctx_sqrt_0(_S255);
    float _S257 = _S241.rows[int(0)].x;
    float _S258 = _S242.rows[int(1)].y;
    float2  _S259 = make_float2 (0.0f);
    float2  _S260 = _S259;
    *&((&_S260)->y) = v_conic_0.z;
    float2  _S261 = _S259;
    *&((&_S261)->y) = v_conic_0.y;
    *&((&_S261)->x) = v_conic_0.x;
    DiffPair_float_0 _S262;
    (&_S262)->primal_0 = _S258;
    (&_S262)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S262, 0.0f);
    DiffPair_float_0 _S263;
    (&_S263)->primal_0 = _S257;
    (&_S263)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S263, 0.0f);
    DiffPair_float_0 _S264;
    (&_S264)->primal_0 = 3.32999992370605469f;
    (&_S264)->differential_0 = 0.0f;
    DiffPair_float_0 _S265;
    (&_S265)->primal_0 = _S256;
    (&_S265)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S264, &_S265, 0.0f);
    DiffPair_float_0 _S266;
    (&_S266)->primal_0 = _S255;
    (&_S266)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S266, _S265.differential_0);
    float _S267 = 2.0f * _S266.differential_0;
    DiffPair_float_0 _S268;
    (&_S268)->primal_0 = _S254;
    (&_S268)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S268, _S267);
    float _S269 = v_opacity_0 + 254.9999847412109375f * _S268.differential_0;
    Matrix<float, 2, 2>  _S270 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S271 = _S270;
    _S271[int(1)] = _S260;
    _S271[int(0)] = _S261;
    Matrix<float, 2, 2>  _S272 = _S271;
    float3  _S273 = make_float3 (0.0f, 0.0f, v_depth_0);
    float2  _S274 = make_float2 (0.0f, _S262.differential_0);
    float2  _S275 = make_float2 (_S263.differential_0, 0.0f);
    float _S276;
    if(antialiased_6)
    {
        float _S277 = _S251 * _S269;
        _S253 = _S246 * _S269;
        _S276 = _S277;
    }
    else
    {
        _S253 = _S269;
        _S276 = 0.0f;
    }
    float _S278 = - (_S253 / _S252);
    DiffPair_float_0 _S279;
    (&_S279)->primal_0 = _S249;
    (&_S279)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S279, _S278);
    float _S280 = - _S279.differential_0;
    float _S281 = invdet_7 * _S272.rows[int(1)].y;
    float _S282 = - (invdet_7 * _S272.rows[int(1)].x);
    float _S283 = - (invdet_7 * _S272.rows[int(0)].y);
    float _S284 = invdet_7 * _S272.rows[int(0)].x;
    float _S285 = - ((_S238 * _S272.rows[int(1)].y + _S248 * _S272.rows[int(1)].x + _S247 * _S272.rows[int(0)].y + _S240 * _S272.rows[int(0)].x) / _S244);
    DiffPair_float_0 _S286;
    (&_S286)->primal_0 = _S245;
    (&_S286)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S286, _S276);
    DiffPair_float_0 _S287;
    (&_S287)->primal_0 = 0.0f;
    (&_S287)->differential_0 = 0.0f;
    DiffPair_float_0 _S288;
    (&_S288)->primal_0 = _S243;
    (&_S288)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S287, &_S288, _S286.differential_0);
    float _S289 = _S288.differential_0 / _S244;
    float s_diff_det_orig_T_0 = det_blur_4 * _S289;
    float _S290 = _S285 + det_orig_7 * - _S289;
    float _S291 = - _S290;
    float _S292 = _S238 * _S290;
    float _S293 = _S240 * _S290;
    Matrix<float, 2, 2>  _S294 = _S270;
    _S294[int(1)] = _S274;
    _S294[int(0)] = _S275;
    _S239 = _S294;
    *&(((&_S239)->rows + (int(1)))->y) = 0.0f;
    float _S295 = _S284 + _S292 + _S294.rows[int(1)].y;
    *&(((&_S239)->rows + (int(0)))->x) = 0.0f;
    float _S296 = _S281 + _S293 + _S294.rows[int(0)].x;
    float _S297 = _S291 + - s_diff_det_orig_T_0;
    float _S298 = _S282 + _S234.rows[int(0)].y * _S297;
    float _S299 = _S283 + _S234.rows[int(1)].x * _S297;
    float _S300 = _S234.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S301 = _S295 + _S234.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S302 = _S259;
    *&((&_S302)->x) = _S298;
    *&((&_S302)->y) = _S301;
    float _S303 = _S296 + _S300;
    float2  _S304 = _S259;
    *&((&_S304)->y) = _S299;
    *&((&_S304)->x) = _S303;
    float _S305 = _S236 * v_mean2d_0.y;
    float _S306 = fy_10 * (rz_4 * v_mean2d_0.y);
    float _S307 = _S235 * v_mean2d_0.x;
    float _S308 = fx_10 * (rz_4 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S309 = _S270;
    _S309[int(1)] = _S302;
    _S309[int(0)] = _S304;
    Matrix<float, 2, 2>  _S310 = _S239 + _S309;
    Matrix<float, 2, 3>  _S311 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S312;
    (&_S312)->primal_0 = _S232;
    (&_S312)->differential_0 = _S311;
    Matrix<float, 3, 2>  _S313 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S314;
    (&_S314)->primal_0 = _S233;
    (&_S314)->differential_0 = _S313;
    s_bwd_prop_mul_0(&_S312, &_S314, _S310);
    Matrix<float, 2, 3>  _S315 = transpose_2(_S314.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S316;
    (&_S316)->primal_0 = J_8;
    (&_S316)->differential_0 = _S311;
    Matrix<float, 3, 3>  _S317 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S318;
    (&_S318)->primal_0 = _S214;
    (&_S318)->differential_0 = _S317;
    s_bwd_prop_mul_1(&_S316, &_S318, _S312.differential_0);
    Matrix<float, 2, 3>  _S319 = _S315 + _S316.differential_0;
    float _S320 = _S231 * _S319.rows[int(1)].z;
    float s_diff_ty_T_0 = _S230 * (rz2_4 * _S319.rows[int(1)].z);
    float _S321 = fy_10 * _S319.rows[int(1)].y;
    float _S322 = _S229 * _S319.rows[int(0)].z;
    float s_diff_tx_T_0 = _S228 * (rz2_4 * _S319.rows[int(0)].z);
    float _S323 = fx_10 * _S319.rows[int(0)].x;
    float _S324 = mean_c_6.z * s_diff_ty_T_0;
    float _S325 = _S227 * s_diff_ty_T_0;
    DiffPair_float_0 _S326;
    (&_S326)->primal_0 = lim_y_pos_0;
    (&_S326)->differential_0 = 0.0f;
    DiffPair_float_0 _S327;
    (&_S327)->primal_0 = _S226;
    (&_S327)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S326, &_S327, _S324);
    DiffPair_float_0 _S328;
    (&_S328)->primal_0 = _S224;
    (&_S328)->differential_0 = 0.0f;
    DiffPair_float_0 _S329;
    (&_S329)->primal_0 = _S225;
    (&_S329)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S328, &_S329, _S327.differential_0);
    float _S330 = mean_c_6.y * _S329.differential_0;
    float _S331 = rz_4 * _S329.differential_0;
    float _S332 = mean_c_6.z * s_diff_tx_T_0;
    float _S333 = _S223 * s_diff_tx_T_0;
    DiffPair_float_0 _S334;
    (&_S334)->primal_0 = lim_x_pos_0;
    (&_S334)->differential_0 = 0.0f;
    DiffPair_float_0 _S335;
    (&_S335)->primal_0 = _S222;
    (&_S335)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S334, &_S335, _S332);
    DiffPair_float_0 _S336;
    (&_S336)->primal_0 = _S220;
    (&_S336)->differential_0 = 0.0f;
    DiffPair_float_0 _S337;
    (&_S337)->primal_0 = _S221;
    (&_S337)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S336, &_S337, _S335.differential_0);
    float _S338 = rz_4 * (_S320 + _S322);
    float _S339 = _S325 + _S333 + - ((_S305 + _S307 + _S321 + _S323 + _S330 + mean_c_6.x * _S337.differential_0 + _S338 + _S338) / _S219);
    float _S340 = _S306 + _S331;
    float _S341 = _S308 + rz_4 * _S337.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S342;
    (&_S342)->primal_0 = _S212;
    (&_S342)->differential_0 = _S317;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S343;
    (&_S343)->primal_0 = _S213;
    (&_S343)->differential_0 = _S317;
    s_bwd_prop_mul_2(&_S342, &_S343, _S318.differential_0);
    Matrix<float, 3, 3>  _S344 = transpose_0(_S343.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S345;
    (&_S345)->primal_0 = R_8;
    (&_S345)->differential_0 = _S317;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S346;
    (&_S346)->primal_0 = _S211;
    (&_S346)->differential_0 = _S317;
    s_bwd_prop_mul_2(&_S345, &_S346, _S342.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S347;
    (&_S347)->primal_0 = _S209;
    (&_S347)->differential_0 = _S317;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S348;
    (&_S348)->primal_0 = _S210;
    (&_S348)->differential_0 = _S317;
    s_bwd_prop_mul_2(&_S347, &_S348, _S346.differential_0);
    Matrix<float, 3, 3>  _S349 = _S347.differential_0 + transpose_0(_S348.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S350;
    (&_S350)->primal_0 = _S208;
    (&_S350)->differential_0 = _S317;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S351;
    (&_S351)->primal_0 = S_0;
    (&_S351)->differential_0 = _S317;
    s_bwd_prop_mul_2(&_S350, &_S351, _S349);
    Matrix<float, 3, 3>  _S352 = transpose_0(_S350.differential_0);
    float _S353 = 2.0f * - _S352.rows[int(2)].z;
    float _S354 = 2.0f * _S352.rows[int(2)].y;
    float _S355 = 2.0f * _S352.rows[int(2)].x;
    float _S356 = 2.0f * _S352.rows[int(1)].z;
    float _S357 = 2.0f * - _S352.rows[int(1)].y;
    float _S358 = 2.0f * _S352.rows[int(1)].x;
    float _S359 = 2.0f * _S352.rows[int(0)].z;
    float _S360 = 2.0f * _S352.rows[int(0)].y;
    float _S361 = 2.0f * - _S352.rows[int(0)].x;
    float _S362 = - _S358 + _S360;
    float _S363 = _S355 + - _S359;
    float _S364 = - _S354 + _S356;
    float _S365 = _S354 + _S356;
    float _S366 = _S355 + _S359;
    float _S367 = _S358 + _S360;
    float _S368 = z_9 * (_S357 + _S361);
    float _S369 = y_11 * (_S353 + _S361);
    float _S370 = x_26 * (_S353 + _S357);
    float _S371 = z_9 * _S362 + y_11 * _S363 + x_26 * _S364;
    float _S372 = _S207 * _S371;
    float _S373 = w_9 * _S362 + y_11 * _S365 + x_26 * _S366 + _S368 + _S368;
    float _S374 = _S207 * _S373;
    float _S375 = w_9 * _S363 + z_9 * _S365 + x_26 * _S367 + _S369 + _S369;
    float _S376 = _S207 * _S375;
    float _S377 = w_9 * _S364 + z_9 * _S366 + y_11 * _S367 + _S370 + _S370;
    float _S378 = _S207 * _S377;
    float _S379 = quat_9.x * _S371 + quat_9.w * _S373 + quat_9.z * _S375 + quat_9.y * _S377;
    DiffPair_float_0 _S380;
    (&_S380)->primal_0 = _S206;
    (&_S380)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S380, _S379);
    float _S381 = quat_9.x * _S380.differential_0;
    float _S382 = quat_9.w * _S380.differential_0;
    float _S383 = quat_9.z * _S380.differential_0;
    float _S384 = quat_9.y * _S380.differential_0;
    float _S385 = _S374 + _S382 + _S382;
    float _S386 = _S376 + _S383 + _S383;
    float _S387 = _S378 + _S384 + _S384;
    float _S388 = _S372 + _S381 + _S381;
    float3  _S389 = make_float3 (0.0f);
    float3  _S390 = _S389;
    *&((&_S390)->z) = _S351.differential_0.rows[int(2)].z;
    *&((&_S390)->y) = _S351.differential_0.rows[int(1)].y;
    *&((&_S390)->x) = _S351.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S391;
    (&_S391)->primal_0 = scale_8;
    (&_S391)->differential_0 = _S389;
    s_bwd_prop_exp_1(&_S391, _S390);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S392 = _S391;
    float3  _S393 = _S389;
    *&((&_S393)->z) = _S339;
    *&((&_S393)->y) = _S340;
    *&((&_S393)->x) = _S341;
    float3  _S394 = _S273 + _S393;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S395;
    (&_S395)->primal_0 = R_8;
    (&_S395)->differential_0 = _S317;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S396;
    (&_S396)->primal_0 = mean_6;
    (&_S396)->differential_0 = _S389;
    s_bwd_prop_mul_3(&_S395, &_S396, _S394);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S397 = _S396;
    Matrix<float, 3, 3>  _S398 = _S344 + _S345.differential_0 + _S395.differential_0;
    float4  _S399 = make_float4 (0.0f);
    *&((&_S399)->w) = _S385;
    *&((&_S399)->z) = _S386;
    *&((&_S399)->y) = _S387;
    *&((&_S399)->x) = _S388;
    float4  _S400 = _S399;
    *v_mean_0 = _S397.differential_0;
    *v_quat_0 = _S400;
    *v_scale_0 = _S392.differential_0;
    *v_in_opacity_0 = _S280;
    *v_R_0 = _S398;
    *v_t_0 = _S394;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S401;
    DiffPair_float_0 _S402;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S403;
    DiffPair_float_0 _S404;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S405;
    DiffPair_float_0 _S406;
    DiffPair_float_0 _S407;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S408;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S409;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S410;
};

inline __device__ CameraDistortion_0 s_primal_ctx_CameraDistortion_x24init_0(float4  dpradial_coeffs_0, float2  dptangential_coeffs_0, float2  dpthin_prism_coeffs_0)
{
    CameraDistortion_0 _S411 = { dpradial_coeffs_0, dptangential_coeffs_0, dpthin_prism_coeffs_0 };
    return _S411;
}

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S412, float _S413)
{
    return s_primal_ctx_atan2_0(_S412, _S413);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S414;
    DiffPair_float_0 _S415;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_0)
{
    DiffPair_float_0 _S416 = { 0.0f, 0.0f };
    _s_diff_ctx_0->_S414 = _S416;
    _s_diff_ctx_0->_S415 = _S416;
    (&_s_diff_ctx_0->_S414)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S414)->differential_0 = 0.0f;
    (&_s_diff_ctx_0->_S415)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S415)->differential_0 = 0.0f;
    DiffPair_float_0 _S417 = *dpdpy_0;
    _s_diff_ctx_0->_S414 = *dpdpy_0;
    DiffPair_float_0 _S418 = *dpdpx_0;
    _s_diff_ctx_0->_S415 = *dpdpx_0;
    float _S419 = _S418.primal_0 * _S418.primal_0 + _S417.primal_0 * _S417.primal_0;
    float _S420 = - _S417.primal_0 / _S419 * dpdOut_0;
    float _S421 = _S418.primal_0 / _S419 * dpdOut_0;
    dpdpy_0->primal_0 = _S417.primal_0;
    dpdpy_0->differential_0 = _S421;
    dpdpx_0->primal_0 = _S418.primal_0;
    dpdpx_0->differential_0 = _S420;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S422, DiffPair_float_0 * _S423, float _S424, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_float_0 _S425 = { 0.0f, 0.0f };
    _s_diff_ctx_1->_S401 = _S425;
    _s_diff_ctx_1->_S402 = _S425;
    (&_s_diff_ctx_1->_S401)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S401)->differential_0 = 0.0f;
    (&_s_diff_ctx_1->_S402)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S402)->differential_0 = 0.0f;
    DiffPair_float_0 _S426 = *_S422;
    _s_diff_ctx_1->_S401 = *_S422;
    DiffPair_float_0 _S427 = *_S423;
    _s_diff_ctx_1->_S402 = *_S423;
    DiffPair_float_0 _S428 = _S426;
    DiffPair_float_0 _S429 = _S427;
    s_bwd_prop_d_atan2_Intermediates_0 _S430;
    (&_S430)->_S414 = _S425;
    (&_S430)->_S415 = _S425;
    s_primal_ctx_d_atan2_0(&_S428, &_S429, _S424, &_S430);
    *_S422 = _S428;
    *_S423 = _S429;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S431;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S432;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S433;
    DiffPair_float_0 _S434;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S435;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S436;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S437 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S436 = _S437;
    (&_s_diff_ctx_2->_S436)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S436)->differential_0 = 0.0f;
    DiffPair_float_0 _S438 = *dpdpx_1;
    _s_diff_ctx_2->_S436 = *dpdpx_1;
    float _S439 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S438.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S438.primal_0;
    dpdpx_1->differential_0 = _S439;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S440, float _S441, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S442 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S432 = _S442;
    (&_s_diff_ctx_3->_S432)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S432)->differential_0 = 0.0f;
    DiffPair_float_0 _S443 = *_S440;
    _s_diff_ctx_3->_S432 = *_S440;
    DiffPair_float_0 _S444 = _S443;
    s_bwd_prop_d_sqrt_Intermediates_0 _S445;
    (&_S445)->_S436 = _S442;
    s_primal_ctx_d_sqrt_0(&_S444, _S441, &_S445);
    *_S440 = _S444;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_4)
{
    float2  _S446 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S447 = { _S446, _S446 };
    DiffPair_float_0 _S448 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S449 = { _S448 };
    _s_diff_ctx_4->_S433 = _S447;
    _s_diff_ctx_4->_S434 = _S448;
    _s_diff_ctx_4->_S435 = _S449;
    (&_s_diff_ctx_4->_S433)->primal_0 = _S446;
    (&_s_diff_ctx_4->_S433)->differential_0 = _S446;
    (&_s_diff_ctx_4->_S434)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S434)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S450 = *dpdpx_2;
    _s_diff_ctx_4->_S433 = *dpdpx_2;
    float _S451 = _S450.primal_0.x;
    float _S452 = _S450.primal_0.y;
    DiffPair_float_0 _S453;
    (&_S453)->primal_0 = _S451 * _S451 + _S452 * _S452;
    (&_S453)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S453, dp_s_dOut_0, &_s_diff_ctx_4->_S435);
    _s_diff_ctx_4->_S434 = _S453;
    float _S454 = _S450.primal_0.y * _S453.differential_0;
    float _S455 = _S454 + _S454;
    float _S456 = _S450.primal_0.x * _S453.differential_0;
    float _S457 = _S456 + _S456;
    float2  _S458 = _S446;
    *&((&_S458)->y) = _S455;
    *&((&_S458)->x) = _S457;
    dpdpx_2->primal_0 = _S450.primal_0;
    dpdpx_2->differential_0 = _S458;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S459, float _S460, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_5)
{
    float2  _S461 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S462 = { _S461, _S461 };
    _s_diff_ctx_5->_S431 = _S462;
    (&_s_diff_ctx_5->_S431)->primal_0 = _S461;
    (&_s_diff_ctx_5->_S431)->differential_0 = _S461;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S463 = *_S459;
    _s_diff_ctx_5->_S431 = *_S459;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S464 = _S463;
    DiffPair_float_0 _S465 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S466 = { _S465 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S467;
    (&_S467)->_S433 = _S462;
    (&_S467)->_S434 = _S465;
    (&_S467)->_S435 = _S466;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S464, _S460, &_S467);
    *_S459 = _S464;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, float4  dpradial_coeffs_1, float2  dptangential_coeffs_1, float2  dpthin_prism_coeffs_1, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_6)
{
    DiffPair_float_0 _S468 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S469 = { _S468, _S468 };
    _s_diff_ctx_6->_S403 = _S468;
    _s_diff_ctx_6->_S404 = _S468;
    _s_diff_ctx_6->_S405 = _S469;
    _s_diff_ctx_6->_S406 = _S468;
    _s_diff_ctx_6->_S407 = _S468;
    _s_diff_ctx_6->_S408 = _S469;
    (&_s_diff_ctx_6->_S403)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S403)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S404)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S404)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S406)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S406)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S407)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S407)->differential_0 = 0.0f;
    float2  _S470 = make_float2 (0.0f);
    CameraDistortion_0 _S471 = s_primal_ctx_CameraDistortion_x24init_0(dpradial_coeffs_1, dptangential_coeffs_1, dpthin_prism_coeffs_1);
    float2  _S472 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S473 = length_0(_S472);
    float _S474 = dpmean3d_0.z;
    float _S475 = s_primal_ctx_atan2_0(_S473, _S474);
    float k_1;
    if(_S475 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S475 * _S475 / 3.0f) / _S474;
    }
    else
    {
        k_1 = _S475 / _S473;
    }
    float2  _S476 = _S472 * make_float2 (k_1);
    float k1_1 = _S471.radial_coeffs_0.x;
    float k2_1 = _S471.radial_coeffs_0.y;
    float k3_1 = _S471.radial_coeffs_0.z;
    float k4_1 = _S471.radial_coeffs_0.w;
    float p1_1 = _S471.tangential_coeffs_0.x;
    float p2_1 = _S471.tangential_coeffs_0.y;
    float sx1_1 = _S471.thin_prism_coeffs_0.x;
    float sy1_1 = _S471.thin_prism_coeffs_0.y;
    float u_3 = _S476.x;
    float v_3 = _S476.y;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float _S477 = 2.0f * p1_1;
    float _S478 = 2.0f * p2_1;
    float2  _S479 = _S476 * make_float2 (1.0f + r2_3 * (k1_1 + r2_3 * (k2_1 + r2_3 * (k3_1 + r2_3 * k4_1)))) + make_float2 (_S477 * u_3 * v_3 + p2_1 * (r2_3 + 2.0f * u_3 * u_3) + sx1_1 * r2_3, _S478 * u_3 * v_3 + p1_1 * (r2_3 + 2.0f * v_3 * v_3) + sy1_1 * r2_3);
    float2  _S480 = make_float2 (dpfx_0 * _S479.x + dpcx_0, dpfy_0 * _S479.y + dpcy_0);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (0.0f);
    float _S481 = s_primal_ctx_s_primal_ctx_atan2_0(_S473, _S474);
    bool _S482 = _S481 < 0.00100000004749745f;
    float _S483;
    float _S484;
    float _S485;
    if(_S482)
    {
        float _S486 = 1.0f - _S481 * _S481 / 3.0f;
        float _S487 = _S474 * _S474;
        k_1 = _S486 / _S474;
        _S483 = 0.0f;
        _S484 = _S487;
        _S485 = _S486;
    }
    else
    {
        float _S488 = _S473 * _S473;
        k_1 = _S481 / _S473;
        _S483 = _S488;
        _S484 = 0.0f;
        _S485 = 0.0f;
    }
    float2  _S489 = make_float2 (k_1);
    float2  _S490 = _S472 * make_float2 (k_1);
    float u_4 = _S490.x;
    float v_4 = _S490.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float _S491 = k3_1 + r2_4 * k4_1;
    float _S492 = k2_1 + r2_4 * _S491;
    float _S493 = k1_1 + r2_4 * _S492;
    float2  _S494 = make_float2 (dpfx_0, 0.0f);
    float2  _S495 = _S490 * _S494;
    float _S496 = p2_1 * dpfx_0;
    float _S497 = _S495.x + _S495.y;
    float _S498 = r2_4 * _S497;
    float _S499 = r2_4 * _S498;
    float _S500 = sx1_1 * dpfx_0 + _S496 + _S493 * _S497 + _S492 * _S498 + _S491 * _S499 + k4_1 * (r2_4 * _S499);
    float _S501 = v_4 * _S500;
    float _S502 = u_4 * _S500;
    float2  _S503 = make_float2 (1.0f + r2_4 * _S493) * _S494 + make_float2 (2.0f * u_4 * _S496 + 2.0f * (u_4 * _S496) + _S477 * (v_4 * dpfx_0) + _S502 + _S502, _S477 * u_4 * dpfx_0 + _S501 + _S501);
    float2  _S504 = _S472 * _S503;
    float2  _S505 = _S489 * _S503;
    float _S506 = _S504.x + _S504.y;
    if(_S482)
    {
        float _S507 = _S506 / _S484;
        float _S508 = _S485 * - _S507;
        float _S509 = _S481 * (0.3333333432674408f * - (_S474 * _S507));
        k_1 = _S509 + _S509;
        _S483 = _S508;
        _S484 = 0.0f;
    }
    else
    {
        float _S510 = _S506 / _S483;
        float _S511 = _S481 * - _S510;
        k_1 = _S473 * _S510;
        _S483 = 0.0f;
        _S484 = _S511;
    }
    DiffPair_float_0 _S512;
    (&_S512)->primal_0 = _S473;
    (&_S512)->differential_0 = 0.0f;
    DiffPair_float_0 _S513;
    (&_S513)->primal_0 = _S474;
    (&_S513)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S512, &_S513, k_1, &_s_diff_ctx_6->_S405);
    _s_diff_ctx_6->_S403 = _S512;
    _s_diff_ctx_6->_S404 = _S513;
    float _S514 = _S513.differential_0 + _S483;
    float _S515 = _S512.differential_0 + _S484;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S516 = { _S470, _S470 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S517;
    (&_S517)->primal_0 = _S472;
    (&_S517)->differential_0 = _S470;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S518;
    (&_S518)->_S431 = _S516;
    s_primal_ctx_s_bwd_length_impl_0(&_S517, _S515, &_S518);
    float2  _S519 = _S517.differential_0 + _S505;
    float3  _S520 = make_float3 (_S519.x, _S519.y, _S514);
    Matrix<float, 2, 3>  _S521 = J_9;
    _S521[int(0)] = _S520;
    if(_S482)
    {
        float _S522 = 1.0f - _S481 * _S481 / 3.0f;
        float _S523 = _S474 * _S474;
        k_1 = _S522 / _S474;
        _S483 = 0.0f;
        _S484 = _S523;
        _S485 = _S522;
    }
    else
    {
        float _S524 = _S473 * _S473;
        k_1 = _S481 / _S473;
        _S483 = _S524;
        _S484 = 0.0f;
        _S485 = 0.0f;
    }
    float2  _S525 = make_float2 (k_1);
    float2  _S526 = _S472 * make_float2 (k_1);
    float u_5 = _S526.x;
    float v_5 = _S526.y;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float _S527 = k3_1 + r2_5 * k4_1;
    float _S528 = k2_1 + r2_5 * _S527;
    float _S529 = k1_1 + r2_5 * _S528;
    float2  _S530 = make_float2 (0.0f, dpfy_0);
    float2  _S531 = _S526 * _S530;
    float _S532 = p1_1 * dpfy_0;
    float _S533 = _S531.x + _S531.y;
    float _S534 = r2_5 * _S533;
    float _S535 = r2_5 * _S534;
    float _S536 = sy1_1 * dpfy_0 + _S532 + _S529 * _S533 + _S528 * _S534 + _S527 * _S535 + k4_1 * (r2_5 * _S535);
    float _S537 = v_5 * _S536;
    float _S538 = u_5 * _S536;
    float2  _S539 = make_float2 (1.0f + r2_5 * _S529) * _S530 + make_float2 (_S478 * (v_5 * dpfy_0) + _S538 + _S538, 2.0f * v_5 * _S532 + 2.0f * (v_5 * _S532) + _S478 * u_5 * dpfy_0 + _S537 + _S537);
    float2  _S540 = _S472 * _S539;
    float2  _S541 = _S525 * _S539;
    float _S542 = _S540.x + _S540.y;
    if(_S482)
    {
        float _S543 = _S542 / _S484;
        float _S544 = _S485 * - _S543;
        float _S545 = _S481 * (0.3333333432674408f * - (_S474 * _S543));
        k_1 = _S545 + _S545;
        _S483 = _S544;
        _S484 = 0.0f;
    }
    else
    {
        float _S546 = _S542 / _S483;
        float _S547 = _S481 * - _S546;
        k_1 = _S473 * _S546;
        _S483 = 0.0f;
        _S484 = _S547;
    }
    DiffPair_float_0 _S548;
    (&_S548)->primal_0 = _S473;
    (&_S548)->differential_0 = 0.0f;
    DiffPair_float_0 _S549;
    (&_S549)->primal_0 = _S474;
    (&_S549)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S548, &_S549, k_1, &_s_diff_ctx_6->_S408);
    _s_diff_ctx_6->_S406 = _S548;
    _s_diff_ctx_6->_S407 = _S549;
    float _S550 = _S549.differential_0 + _S483;
    float _S551 = _S548.differential_0 + _S484;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S552;
    (&_S552)->primal_0 = _S472;
    (&_S552)->differential_0 = _S470;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S553;
    (&_S553)->_S431 = _S516;
    s_primal_ctx_s_bwd_length_impl_0(&_S552, _S551, &_S553);
    float2  _S554 = _S552.differential_0 + _S541;
    float3  _S555 = make_float3 (_S554.x, _S554.y, _S550);
    _S521[int(1)] = _S555;
    *dpcov2d_0 = s_primal_ctx_mul_3(s_primal_ctx_mul_2(_S521, dpcov3d_0), transpose_1(_S521));
    *dpmean2d_0 = _S480;
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

struct s_bwd_prop_d_length_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S556;
};

struct DiffPair_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 primal_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 differential_0;
};

struct DiffPair_1
{
    DiffPair_float_0 primal_0;
    DiffPair_float_0 differential_0;
};

inline __device__ void s_bwd_prop_d_sqrt_0(DiffPair_1 * dpdpx_3, DiffPair_float_0 * dpdOut_2, s_bwd_prop_d_sqrt_Intermediates_0 _s_diff_ctx_7)
{
    DiffPair_1 _S557 = *dpdpx_3;
    float _S558 = s_primal_ctx_max_0(1.00000001168609742e-07f, _s_diff_ctx_7._S436.primal_0);
    float _S559 = s_primal_ctx_sqrt_0(_S558);
    float _S560 = 0.5f / _S559 * (*dpdpx_3).differential_0.differential_0;
    float _S561 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S559 * _S559));
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = _S558;
    (&_S562)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S562, _S561);
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = 1.00000001168609742e-07f;
    (&_S563)->differential_0 = 0.0f;
    DiffPair_float_0 _S564;
    (&_S564)->primal_0 = _s_diff_ctx_7._S436.primal_0;
    (&_S564)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S563, &_S564, _S562.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S564.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S560;
    dpdpx_3->primal_0 = _S557.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S565, DiffPair_float_0 * _S566, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _s_diff_ctx_8)
{
    DiffPair_1 _S567 = *_S565;
    DiffPair_float_0 _S568 = _s_diff_ctx_8._S432;
    DiffPair_float_0 _S569 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S570;
    (&_S570)->_S436 = _S569;
    s_primal_ctx_d_sqrt_0(&_S568, (*_S566).primal_0, &_S570);
    DiffPair_float_0 _S571 = { (*_S565).differential_0.primal_0, (*_S565).differential_0.differential_0 };
    DiffPair_1 _S572;
    (&_S572)->primal_0 = _s_diff_ctx_8._S432;
    (&_S572)->differential_0 = _S571;
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = (*_S566).primal_0;
    (&_S573)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_0(&_S572, &_S573, _S570);
    DiffPair_float_0 _S574 = { _S572.differential_0.primal_0, _S572.differential_0.differential_0 };
    _S566->primal_0 = (*_S566).primal_0;
    _S566->differential_0 = _S573.differential_0;
    _S565->primal_0 = _S567.primal_0;
    _S565->differential_0 = _S574;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S575, float _s_dOut_1)
{
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = (*_S575).primal_0;
    (&_S576)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S576, _s_dOut_1);
    _S575->primal_0 = (*_S575).primal_0;
    _S575->differential_0 = _S576.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _s_diff_ctx_9)
{
    DiffPair_0 _S577 = *dpdpx_5;
    float _S578 = _s_diff_ctx_9._S433.primal_0.x;
    float _S579 = _s_diff_ctx_9._S433.primal_0.y;
    float len_0 = _S578 * _S578 + _S579 * _S579;
    DiffPair_float_0 _S580 = { len_0, 0.0f };
    float2  _S581 = make_float2 (0.0f);
    float _S582 = (*dpdpx_5).differential_0.differential_0.x;
    float _S583 = _S582 + _S582;
    float _S584 = _s_diff_ctx_9._S434.differential_0 * _S583;
    float _S585 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S586 = _s_diff_ctx_9._S434.differential_0 * _S585;
    DiffPair_float_0 _S587 = { 0.0f, _s_diff_ctx_9._S433.primal_0.x * _S583 + _s_diff_ctx_9._S433.primal_0.y * _S585 };
    DiffPair_1 _S588;
    (&_S588)->primal_0 = _S580;
    (&_S588)->differential_0 = _S587;
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S589)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S588, &_S589, _s_diff_ctx_9._S435);
    DiffPair_float_0 _S590;
    (&_S590)->primal_0 = len_0;
    (&_S590)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S590, 0.0f);
    float _S591 = _S588.differential_0.primal_0 + _S590.differential_0;
    float _S592 = _s_diff_ctx_9._S433.primal_0.y * _S591;
    float _S593 = _S586 + _S592 + _S592;
    float _S594 = _s_diff_ctx_9._S433.primal_0.x * _S591;
    float _S595 = _S584 + _S594 + _S594;
    float2  _S596 = _S581;
    *&((&_S596)->y) = _S593;
    *&((&_S596)->x) = _S595;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { (*dpdpx_5).differential_0.primal_0 + _S596, _S581 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S589.differential_0;
    dpdpx_5->primal_0 = _S577.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_2)
{
    float _S597 = (*dpdpx_7).primal_0.x;
    float _S598 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S599;
    (&_S599)->primal_0 = _S597 * _S597 + _S598 * _S598;
    (&_S599)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S599, _s_dOut_2);
    float _S600 = (*dpdpx_7).primal_0.y * _S599.differential_0;
    float _S601 = _S600 + _S600;
    float _S602 = (*dpdpx_7).primal_0.x * _S599.differential_0;
    float _S603 = _S602 + _S602;
    float2  _S604 = make_float2 (0.0f);
    *&((&_S604)->y) = _S601;
    *&((&_S604)->x) = _S603;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S604;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S605, DiffPair_float_0 * _S606, s_bwd_prop_s_bwd_length_impl_Intermediates_0 _s_diff_ctx_10)
{
    DiffPair_0 _S607 = *_S605;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S608 = _s_diff_ctx_10._S431;
    float2  _S609 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S610 = { _S609, _S609 };
    DiffPair_float_0 _S611 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S612 = { _S611 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S613;
    (&_S613)->_S433 = _S610;
    (&_S613)->_S434 = _S611;
    (&_S613)->_S435 = _S612;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S608, (*_S606).primal_0, &_S613);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S614 = { (*_S605).differential_0.primal_0, (*_S605).differential_0.differential_0 };
    DiffPair_0 _S615;
    (&_S615)->primal_0 = _s_diff_ctx_10._S431;
    (&_S615)->differential_0 = _S614;
    DiffPair_float_0 _S616;
    (&_S616)->primal_0 = (*_S606).primal_0;
    (&_S616)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S615, &_S616, _S613);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S617;
    (&_S617)->primal_0 = _s_diff_ctx_10._S431.primal_0;
    (&_S617)->differential_0 = _S609;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S617, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S618 = { _S615.differential_0.primal_0 + _S617.differential_0, _S615.differential_0.differential_0 };
    _S606->primal_0 = (*_S606).primal_0;
    _S606->differential_0 = _S616.differential_0;
    _S605->primal_0 = _S607.primal_0;
    _S605->differential_0 = _S618;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 _s_diff_ctx_11)
{
    DiffPair_1 _S619 = *dpdpy_1;
    DiffPair_1 _S620 = *dpdpx_8;
    float _S621 = - _s_diff_ctx_11._S414.primal_0;
    float _S622 = _s_diff_ctx_11._S415.primal_0 * _s_diff_ctx_11._S415.primal_0 + _s_diff_ctx_11._S414.primal_0 * _s_diff_ctx_11._S414.primal_0;
    float _S623 = _S622 * _S622;
    float _S624 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S623;
    float _S625 = _s_diff_ctx_11._S415.primal_0 * - _S624;
    float _S626 = _s_diff_ctx_11._S414.primal_0 * _S625;
    float _S627 = _s_diff_ctx_11._S415.primal_0 * _S625;
    float _S628 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S623;
    float _S629 = _S621 * - _S628;
    float _S630 = _s_diff_ctx_11._S414.primal_0 * _S629;
    float _S631 = _s_diff_ctx_11._S415.primal_0 * _S629;
    DiffPair_float_0 dpdpx_9 = { _S631 + _S631 + ((*dpdpx_8).differential_0.primal_0 + (_S627 + _S627 + _S622 * _S624)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S626 + _S626 + (*dpdpy_1).differential_0.primal_0 + _S630 + _S630 + - (_S622 * _S628), 0.0f };
    float _S632 = _s_diff_ctx_11._S415.primal_0 / _S622 * (*dpdpy_1).differential_0.differential_0 + _S621 / _S622 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S632;
    dpdpy_1->primal_0 = _S619.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S620.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S633, DiffPair_1 * _S634, DiffPair_float_0 * _S635, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _s_diff_ctx_12)
{
    DiffPair_1 _S636 = *_S633;
    DiffPair_1 _S637 = *_S634;
    DiffPair_float_0 _S638 = _s_diff_ctx_12._S401;
    DiffPair_float_0 _S639 = _s_diff_ctx_12._S402;
    DiffPair_float_0 _S640 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S641;
    (&_S641)->_S414 = _S640;
    (&_S641)->_S415 = _S640;
    s_primal_ctx_d_atan2_0(&_S638, &_S639, (*_S635).primal_0, &_S641);
    DiffPair_float_0 _S642 = { (*_S634).differential_0.primal_0, (*_S634).differential_0.differential_0 };
    DiffPair_float_0 _S643 = { (*_S633).differential_0.primal_0, (*_S633).differential_0.differential_0 };
    DiffPair_1 _S644;
    (&_S644)->primal_0 = _s_diff_ctx_12._S401;
    (&_S644)->differential_0 = _S643;
    DiffPair_1 _S645;
    (&_S645)->primal_0 = _s_diff_ctx_12._S402;
    (&_S645)->differential_0 = _S642;
    DiffPair_float_0 _S646;
    (&_S646)->primal_0 = (*_S635).primal_0;
    (&_S646)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_0(&_S644, &_S645, &_S646, _S641);
    DiffPair_float_0 _S647 = { _S645.differential_0.primal_0, _S645.differential_0.differential_0 };
    DiffPair_float_0 _S648 = { _S644.differential_0.primal_0, _S644.differential_0.differential_0 };
    _S635->primal_0 = (*_S635).primal_0;
    _S635->differential_0 = _S646.differential_0;
    _S633->primal_0 = _S636.primal_0;
    _S633->differential_0 = _S648;
    _S634->primal_0 = _S637.primal_0;
    _S634->differential_0 = _S647;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S649, DiffPair_float_0 * _S650, float _s_dOut_3)
{
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = (*_S649).primal_0;
    (&_S651)->differential_0 = 0.0f;
    DiffPair_float_0 _S652;
    (&_S652)->primal_0 = (*_S650).primal_0;
    (&_S652)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S651, &_S652, _s_dOut_3);
    _S650->primal_0 = (*_S650).primal_0;
    _S650->differential_0 = _S652.differential_0;
    _S649->primal_0 = (*_S649).primal_0;
    _S649->differential_0 = _S651.differential_0;
    return;
}

inline __device__ void s_bwd_prop_CameraDistortion_x24init_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_2, CameraDistortion_0 _s_dOut_4)
{
    dpthin_prism_coeffs_2->primal_0 = (*dpthin_prism_coeffs_2).primal_0;
    dpthin_prism_coeffs_2->differential_0 = _s_dOut_4.thin_prism_coeffs_0;
    dptangential_coeffs_2->primal_0 = (*dptangential_coeffs_2).primal_0;
    dptangential_coeffs_2->differential_0 = _s_dOut_4.tangential_coeffs_0;
    dpradial_coeffs_2->primal_0 = (*dpradial_coeffs_2).primal_0;
    dpradial_coeffs_2->differential_0 = _s_dOut_4.radial_coeffs_0;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_3, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _s_diff_ctx_13)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S653 = *dpcov3d_1;
    DiffPair_float_0 _S654 = *dpfx_1;
    DiffPair_float_0 _S655 = *dpfy_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S656 = *dpradial_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S657 = *dptangential_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S658 = *dpthin_prism_coeffs_3;
    float2  _S659 = make_float2 (0.0f);
    CameraDistortion_0 _S660 = s_primal_ctx_CameraDistortion_x24init_0((*dpradial_coeffs_3).primal_0, (*dptangential_coeffs_3).primal_0, (*dpthin_prism_coeffs_3).primal_0);
    float2  _S661 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S662 = length_0(_S661);
    float _S663 = (*dpmean3d_1).primal_0.z;
    float _S664 = s_primal_ctx_atan2_0(_S662, _S663);
    bool _S665 = _S664 < 0.00100000004749745f;
    float k_2;
    float _S666;
    float _S667;
    float _S668;
    if(_S665)
    {
        float _S669 = 1.0f - _S664 * _S664 / 3.0f;
        float _S670 = _S663 * _S663;
        k_2 = _S669 / _S663;
        _S666 = 0.0f;
        _S667 = _S670;
        _S668 = _S669;
    }
    else
    {
        float _S671 = _S662 * _S662;
        k_2 = _S664 / _S662;
        _S666 = _S671;
        _S667 = 0.0f;
        _S668 = 0.0f;
    }
    float2  _S672 = make_float2 (k_2);
    float2  _S673 = _S661 * make_float2 (k_2);
    float k1_2 = _S660.radial_coeffs_0.x;
    float k2_2 = _S660.radial_coeffs_0.y;
    float k3_2 = _S660.radial_coeffs_0.z;
    float k4_2 = _S660.radial_coeffs_0.w;
    float p1_2 = _S660.tangential_coeffs_0.x;
    float p2_2 = _S660.tangential_coeffs_0.y;
    float sx1_2 = _S660.thin_prism_coeffs_0.x;
    float sy1_2 = _S660.thin_prism_coeffs_0.y;
    float u_6 = _S673.x;
    float v_6 = _S673.y;
    float r2_6 = u_6 * u_6 + v_6 * v_6;
    float _S674 = k3_2 + r2_6 * k4_2;
    float _S675 = k2_2 + r2_6 * _S674;
    float _S676 = k1_2 + r2_6 * _S675;
    float radial_0 = 1.0f + r2_6 * _S676;
    float2  _S677 = make_float2 (radial_0);
    float _S678 = 2.0f * p1_2;
    float _S679 = _S678 * u_6;
    float _S680 = 2.0f * u_6;
    float _S681 = r2_6 + _S680 * u_6;
    float _S682 = 2.0f * p2_2;
    float _S683 = _S682 * u_6;
    float _S684 = 2.0f * v_6;
    float _S685 = r2_6 + _S684 * v_6;
    float2  _S686 = _S673 * make_float2 (radial_0) + make_float2 (_S679 * v_6 + p2_2 * _S681 + sx1_2 * r2_6, _S683 * v_6 + p1_2 * _S685 + sy1_2 * r2_6);
    float _S687 = _S686.x;
    float _S688 = _S686.y;
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float _S689 = s_primal_ctx_s_primal_ctx_atan2_0(_S662, _S663);
    bool _S690 = _S689 < 0.00100000004749745f;
    float _S691;
    float _S692;
    float _S693;
    float _S694;
    float _S695;
    if(_S690)
    {
        float _S696 = 1.0f - _S689 * _S689 / 3.0f;
        float _S697 = _S663 * _S663;
        k_2 = _S696 / _S663;
        _S691 = 0.0f;
        _S692 = _S697;
        _S693 = _S696;
        _S694 = 0.0f;
        _S695 = _S697;
    }
    else
    {
        float _S698 = _S662 * _S662;
        k_2 = _S689 / _S662;
        _S691 = _S698;
        _S692 = 0.0f;
        _S693 = 0.0f;
        _S694 = _S698;
        _S695 = 0.0f;
    }
    float2  _S699 = make_float2 (k_2);
    float2  _S700 = make_float2 (k_2);
    float2  _S701 = _S661 * make_float2 (k_2);
    float u_7 = _S701.x;
    float v_7 = _S701.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S702 = k3_2 + r2_7 * k4_2;
    float _S703 = k2_2 + r2_7 * _S702;
    float _S704 = k1_2 + r2_7 * _S703;
    float2  _S705 = make_float2 (1.0f + r2_7 * _S704);
    float _S706 = _S678 * u_7;
    float _S707 = 2.0f * u_7;
    float2  _S708 = make_float2 (_S654.primal_0, 0.0f);
    float2  _S709 = _S701 * _S708;
    float _S710 = p2_2 * _S654.primal_0;
    float _S711 = v_7 * _S654.primal_0;
    float _S712 = _S709.x + _S709.y;
    float _S713 = r2_7 * _S712;
    float _S714 = r2_7 * _S713;
    float _S715 = r2_7 * _S714;
    float _S716 = sx1_2 * _S654.primal_0 + _S710 + _S704 * _S712 + _S703 * _S713 + _S702 * _S714 + k4_2 * _S715;
    float _S717 = v_7 * _S716;
    float _S718 = u_7 * _S716;
    float2  _S719 = _S705 * _S708 + make_float2 (_S707 * _S710 + 2.0f * (u_7 * _S710) + _S678 * _S711 + _S718 + _S718, _S706 * _S654.primal_0 + _S717 + _S717);
    float2  _S720 = _S661 * _S719;
    float2  _S721 = _S700 * _S719;
    float _S722 = _S720.x + _S720.y;
    float k_3;
    float _S723;
    float _S724;
    float _S725;
    float _S726;
    float _S727;
    float _S728;
    float _S729;
    float _S730;
    if(_S690)
    {
        float _S731 = _S722 / _S692;
        float _S732 = _S692 * _S692;
        float _S733 = - _S731;
        float _S734 = _S693 * _S733;
        float _S735 = 0.3333333432674408f * - (_S663 * _S731);
        float _S736 = _S689 * _S735;
        k_2 = _S736 + _S736;
        k_3 = _S734;
        _S723 = 0.0f;
        _S724 = 0.0f;
        _S725 = 0.0f;
        _S726 = 0.0f;
        _S727 = _S735;
        _S728 = _S731;
        _S729 = _S733;
        _S730 = _S732;
    }
    else
    {
        float _S737 = _S722 / _S691;
        float _S738 = _S691 * _S691;
        float _S739 = - _S737;
        float _S740 = _S689 * _S739;
        k_2 = _S662 * _S737;
        k_3 = 0.0f;
        _S723 = _S740;
        _S724 = _S737;
        _S725 = _S739;
        _S726 = _S738;
        _S727 = 0.0f;
        _S728 = 0.0f;
        _S729 = 0.0f;
        _S730 = 0.0f;
    }
    DiffPair_float_0 _S741 = { _S662, 0.0f };
    DiffPair_float_0 _S742 = { _S663, 0.0f };
    float _S743 = _s_diff_ctx_13._S404.differential_0 + k_3;
    float _S744 = _s_diff_ctx_13._S403.differential_0 + _S723;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S745 = { _S659, _S659 };
    s_bwd_prop_d_length_Intermediates_0 _S746;
    (&_S746)->_S556 = _S745;
    (&(&_S746)->_S556)->primal_0 = _S661;
    (&(&_S746)->_S556)->differential_0 = _S659;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S747;
    (&_S747)->primal_0 = _S661;
    (&_S747)->differential_0 = _S659;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S748;
    (&_S748)->_S431 = _S745;
    s_primal_ctx_s_bwd_length_impl_0(&_S747, _S744, &_S748);
    float2  _S749 = _S747.differential_0 + _S721;
    float3  _S750 = make_float3 (_S749.x, _S749.y, _S743);
    Matrix<float, 2, 3>  _S751 = J_10;
    _S751[int(0)] = _S750;
    float _S752;
    float _S753;
    float _S754;
    float _S755;
    if(_S690)
    {
        float _S756 = 1.0f - _S689 * _S689 / 3.0f;
        float _S757 = _S663 * _S663;
        k_3 = _S756 / _S663;
        _S723 = 0.0f;
        _S752 = _S757;
        _S753 = _S756;
        _S754 = 0.0f;
        _S755 = _S757;
    }
    else
    {
        float _S758 = _S662 * _S662;
        k_3 = _S689 / _S662;
        _S723 = _S758;
        _S752 = 0.0f;
        _S753 = 0.0f;
        _S754 = _S758;
        _S755 = 0.0f;
    }
    float2  _S759 = make_float2 (k_3);
    float2  _S760 = make_float2 (k_3);
    float2  _S761 = _S661 * make_float2 (k_3);
    float u_8 = _S761.x;
    float v_8 = _S761.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S762 = k3_2 + r2_8 * k4_2;
    float _S763 = k2_2 + r2_8 * _S762;
    float _S764 = k1_2 + r2_8 * _S763;
    float2  _S765 = make_float2 (1.0f + r2_8 * _S764);
    float _S766 = _S682 * u_8;
    float _S767 = 2.0f * v_8;
    float2  _S768 = make_float2 (0.0f, _S655.primal_0);
    float2  _S769 = _S761 * _S768;
    float _S770 = p1_2 * _S655.primal_0;
    float _S771 = v_8 * _S655.primal_0;
    float _S772 = _S769.x + _S769.y;
    float _S773 = r2_8 * _S772;
    float _S774 = r2_8 * _S773;
    float _S775 = r2_8 * _S774;
    float _S776 = sy1_2 * _S655.primal_0 + _S770 + _S764 * _S772 + _S763 * _S773 + _S762 * _S774 + k4_2 * _S775;
    float _S777 = v_8 * _S776;
    float _S778 = u_8 * _S776;
    float2  _S779 = _S765 * _S768 + make_float2 (_S682 * _S771 + _S778 + _S778, _S767 * _S770 + 2.0f * (v_8 * _S770) + _S766 * _S655.primal_0 + _S777 + _S777);
    float2  _S780 = _S661 * _S779;
    float2  _S781 = _S760 * _S779;
    float _S782 = _S780.x + _S780.y;
    float _S783;
    float _S784;
    float _S785;
    float _S786;
    float _S787;
    float _S788;
    float _S789;
    float _S790;
    float _S791;
    if(_S690)
    {
        float _S792 = _S782 / _S752;
        float _S793 = _S752 * _S752;
        float _S794 = - _S792;
        float _S795 = _S753 * _S794;
        float _S796 = 0.3333333432674408f * - (_S663 * _S792);
        float _S797 = _S689 * _S796;
        k_3 = _S797 + _S797;
        _S783 = _S795;
        _S784 = 0.0f;
        _S785 = 0.0f;
        _S786 = 0.0f;
        _S787 = 0.0f;
        _S788 = _S796;
        _S789 = _S792;
        _S790 = _S794;
        _S791 = _S793;
    }
    else
    {
        float _S798 = _S782 / _S723;
        float _S799 = _S723 * _S723;
        float _S800 = - _S798;
        float _S801 = _S689 * _S800;
        k_3 = _S662 * _S798;
        _S783 = 0.0f;
        _S784 = _S801;
        _S785 = _S798;
        _S786 = _S800;
        _S787 = _S799;
        _S788 = 0.0f;
        _S789 = 0.0f;
        _S790 = 0.0f;
        _S791 = 0.0f;
    }
    float _S802 = _s_diff_ctx_13._S407.differential_0 + _S783;
    float _S803 = _s_diff_ctx_13._S406.differential_0 + _S784;
    s_bwd_prop_d_length_Intermediates_0 _S804;
    (&_S804)->_S556 = _S745;
    (&(&_S804)->_S556)->primal_0 = _S661;
    (&(&_S804)->_S556)->differential_0 = _S659;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S805;
    (&_S805)->primal_0 = _S661;
    (&_S805)->differential_0 = _S659;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S806;
    (&_S806)->_S431 = _S745;
    s_primal_ctx_s_bwd_length_impl_0(&_S805, _S803, &_S806);
    float2  _S807 = _S805.differential_0 + _S781;
    float3  _S808 = make_float3 (_S807.x, _S807.y, _S802);
    _S751[int(1)] = _S808;
    Matrix<float, 3, 2>  _S809 = transpose_1(_S751);
    CameraDistortion_0 _S810 = CameraDistortion_x24_syn_dzero_0();
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S811;
    (&_S811)->primal_0 = s_primal_ctx_mul_2(_S751, _S653.primal_0);
    (&_S811)->differential_0 = J_10;
    Matrix<float, 3, 2>  _S812 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S813;
    (&_S813)->primal_0 = _S809;
    (&_S813)->differential_0 = _S812;
    s_bwd_prop_mul_0(&_S811, &_S813, dpcov2d_1);
    Matrix<float, 2, 3>  _S814 = transpose_2(_S813.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S815;
    (&_S815)->primal_0 = _S751;
    (&_S815)->differential_0 = J_10;
    Matrix<float, 3, 3>  _S816 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S817;
    (&_S817)->primal_0 = _S653.primal_0;
    (&_S817)->differential_0 = _S816;
    s_bwd_prop_mul_1(&_S815, &_S817, _S811.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S818 = _S817;
    Matrix<float, 2, 3>  _S819 = _S814 + _S815.differential_0;
    float2  _S820 = _S659;
    *&((&_S820)->y) = _S819.rows[int(1)].y;
    *&((&_S820)->x) = _S819.rows[int(1)].x;
    float2  _S821 = _S820;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S822 = _S804._S556;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S823;
    (&_S823)->_S431 = _S745;
    s_primal_ctx_s_bwd_length_impl_0(&_S822, _S803, &_S823);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S824 = { _S659, _S820 };
    DiffPair_0 _S825;
    (&_S825)->primal_0 = _S804._S556;
    (&_S825)->differential_0 = _S824;
    DiffPair_float_0 _S826;
    (&_S826)->primal_0 = _S803;
    (&_S826)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_0(&_S825, &_S826, _S823);
    DiffPair_0 _S827 = _S825;
    DiffPair_float_0 _S828 = _S826;
    DiffPair_float_0 _S829 = { 0.0f, _S819.rows[int(1)].z };
    DiffPair_float_0 _S830 = { 0.0f, _S826.differential_0 };
    DiffPair_1 _S831;
    (&_S831)->primal_0 = _S741;
    (&_S831)->differential_0 = _S830;
    DiffPair_1 _S832;
    (&_S832)->primal_0 = _S742;
    (&_S832)->differential_0 = _S829;
    DiffPair_float_0 _S833;
    (&_S833)->primal_0 = k_3;
    (&_S833)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S831, &_S832, &_S833, _s_diff_ctx_13._S408);
    DiffPair_1 _S834 = _S831;
    DiffPair_1 _S835 = _S832;
    DiffPair_float_0 _S836 = _S833;
    if(_S690)
    {
        float _S837 = _S836.differential_0 + _S836.differential_0;
        float _S838 = _S788 * _S837;
        float _S839 = - (0.3333333432674408f * (_S689 * _S837));
        float _S840 = _S790 * _S819.rows[int(1)].z;
        float _S841 = (_S663 * _S839 + - (_S753 * _S819.rows[int(1)].z)) / _S791;
        float _S842 = _S782 * - _S841;
        float _S843 = _S789 * _S839 + _S835.differential_0.primal_0;
        k_3 = _S752 * _S841;
        _S723 = 0.0f;
        _S752 = _S842;
        _S783 = _S840;
        _S784 = _S834.differential_0.primal_0;
        _S785 = _S838;
        _S786 = _S843;
    }
    else
    {
        float _S844 = _S786 * _S828.differential_0;
        float _S845 = (_S662 * _S836.differential_0 + - (_S689 * _S828.differential_0)) / _S787;
        float _S846 = _S782 * - _S845;
        float _S847 = _S785 * _S836.differential_0 + _S834.differential_0.primal_0;
        k_3 = _S723 * _S845;
        _S723 = _S846;
        _S752 = 0.0f;
        _S783 = 0.0f;
        _S784 = _S847;
        _S785 = _S844;
        _S786 = _S835.differential_0.primal_0;
    }
    float2  _S848 = _S760 * _S821;
    float2  _S849 = _S779 * _S821;
    float2  _S850 = _S659;
    *&((&_S850)->y) = k_3;
    *&((&_S850)->x) = k_3;
    float2  _S851 = _S779 * _S850;
    float2  _S852 = _S848 + _S661 * _S850;
    float _S853 = _S852.x;
    float _S854 = _S853 + _S853;
    float _S855 = _S776 * _S854;
    float _S856 = _S852.y + _S852.y;
    float _S857 = _S776 * _S856;
    float _S858 = u_8 * _S854 + v_8 * _S856;
    float _S859 = k4_2 * _S858;
    float _S860 = _S775 * _S858;
    float _S861 = _S774 * _S858;
    float _S862 = _S774 * _S859;
    float _S863 = _S773 * _S858;
    float _S864 = _S762 * _S858 + r2_8 * _S859;
    float _S865 = _S773 * _S864;
    float _S866 = _S772 * _S858;
    float _S867 = _S763 * _S858 + r2_8 * _S864;
    float _S868 = _S772 * _S867;
    float _S869 = _S764 * _S858 + r2_8 * _S867;
    float _S870 = _S682 * _S852.x;
    float _S871 = _S771 * _S852.x;
    float _S872 = v_8 * _S870;
    float _S873 = _S655.primal_0 * _S870;
    float _S874 = _S766 * _S852.y;
    float _S875 = _S655.primal_0 * _S852.y;
    float _S876 = 2.0f * _S852.y;
    float _S877 = _S770 * _S876;
    float _S878 = _S770 * _S852.y;
    float _S879 = _S858 + v_8 * _S876 + _S767 * _S852.y;
    float _S880 = p1_2 * _S879;
    float _S881 = _S655.primal_0 * _S879;
    float _S882 = sy1_2 * _S858;
    float _S883 = _S655.primal_0 * _S858;
    float2  _S884 = _S765 * _S852;
    float2  _S885 = _S768 * _S852;
    float2  _S886 = _S659;
    *&((&_S886)->y) = _S869;
    *&((&_S886)->x) = _S869;
    float _S887 = _S885.x + _S885.y;
    float _S888 = _S866 + r2_8 * _S887;
    float _S889 = _S863 + r2_8 * _S888;
    float _S890 = _S861 + r2_8 * _S889;
    float _S891 = _S862 + _S865 + _S868 + _S764 * _S887 + _S763 * _S888 + _S762 * _S889 + k4_2 * _S890;
    float _S892 = v_8 * _S891;
    float _S893 = u_8 * _S891;
    float2  _S894 = _S851 + _S827.differential_0.primal_0;
    float2  _S895 = _S768 * _S886 + make_float2 (_S855 + _S682 * _S875 + _S893 + _S893, _S857 + _S873 + _S877 + 2.0f * _S878 + _S892 + _S892);
    float _S896 = _S872 + _S874 + _S880 + _S882 + (_S884 + _S761 * _S886).y;
    float _S897 = _S871 + u_8 * _S875;
    float _S898 = _S860 + r2_8 * _S890;
    float2  _S899 = _S661 * _S895;
    float _S900 = _S849.x + _S849.y + _S899.x + _S899.y;
    float2  _S901 = _S759 * _S895 + _S894;
    if(_S690)
    {
        float _S902 = _S900 / _S755;
        float _S903 = _S689 * (0.3333333432674408f * - (_S783 + _S663 * _S902));
        float _S904 = _S663 * _S752 + _S663 * _S752 + _S753 * - _S902 + _S786;
        k_3 = _S903 + _S903 + _S785;
        _S723 = _S904;
        _S752 = _S784;
    }
    else
    {
        float _S905 = _S900 / _S754;
        float _S906 = _S662 * _S723 + _S662 * _S723 + _S689 * - _S905 + _S784;
        k_3 = _S662 * _S905 + _S785;
        _S723 = _S786;
        _S752 = _S906;
    }
    DiffPair_float_0 _S907;
    (&_S907)->primal_0 = _S662;
    (&_S907)->differential_0 = 0.0f;
    DiffPair_float_0 _S908;
    (&_S908)->primal_0 = _S663;
    (&_S908)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S907, &_S908, k_3);
    float _S909 = _S908.differential_0 + _S723;
    float _S910 = _S907.differential_0 + _S752;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S911;
    (&_S911)->primal_0 = _S661;
    (&_S911)->differential_0 = _S659;
    s_bwd_length_impl_0(&_S911, _S910);
    float2  _S912 = _S911.differential_0 + _S901;
    float3  _S913 = make_float3 (_S912.x, _S912.y, _S909);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S914;
    (&_S914)->primal_0 = _S661;
    (&_S914)->differential_0 = _S659;
    s_bwd_length_impl_0(&_S914, 0.0f);
    float3  _S915 = _S913 + make_float3 (_S914.differential_0.x, _S914.differential_0.y, 0.0f);
    float2  _S916 = _S659;
    *&((&_S916)->y) = _S819.rows[int(0)].y;
    *&((&_S916)->x) = _S819.rows[int(0)].x;
    float2  _S917 = _S916;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S918 = _S746._S556;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S919;
    (&_S919)->_S431 = _S745;
    s_primal_ctx_s_bwd_length_impl_0(&_S918, _S744, &_S919);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S920 = { _S659, _S916 };
    DiffPair_0 _S921;
    (&_S921)->primal_0 = _S746._S556;
    (&_S921)->differential_0 = _S920;
    DiffPair_float_0 _S922;
    (&_S922)->primal_0 = _S744;
    (&_S922)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_0(&_S921, &_S922, _S919);
    DiffPair_0 _S923 = _S921;
    DiffPair_float_0 _S924 = _S922;
    DiffPair_float_0 _S925 = { 0.0f, _S819.rows[int(0)].z };
    DiffPair_float_0 _S926 = { 0.0f, _S922.differential_0 };
    DiffPair_1 _S927;
    (&_S927)->primal_0 = _S741;
    (&_S927)->differential_0 = _S926;
    DiffPair_1 _S928;
    (&_S928)->primal_0 = _S742;
    (&_S928)->differential_0 = _S925;
    DiffPair_float_0 _S929;
    (&_S929)->primal_0 = k_2;
    (&_S929)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S927, &_S928, &_S929, _s_diff_ctx_13._S405);
    DiffPair_1 _S930 = _S927;
    DiffPair_1 _S931 = _S928;
    DiffPair_float_0 _S932 = _S929;
    if(_S690)
    {
        float _S933 = _S932.differential_0 + _S932.differential_0;
        float _S934 = _S727 * _S933;
        float _S935 = - (0.3333333432674408f * (_S689 * _S933));
        float _S936 = _S729 * _S819.rows[int(0)].z;
        float _S937 = (_S663 * _S935 + - (_S693 * _S819.rows[int(0)].z)) / _S730;
        float _S938 = _S722 * - _S937;
        float _S939 = _S728 * _S935 + _S931.differential_0.primal_0;
        k_2 = _S692 * _S937;
        _S691 = 0.0f;
        _S692 = _S938;
        k_3 = _S936;
        _S723 = _S930.differential_0.primal_0;
        _S724 = _S934;
        _S725 = _S939;
    }
    else
    {
        float _S940 = _S725 * _S924.differential_0;
        float _S941 = (_S662 * _S932.differential_0 + - (_S689 * _S924.differential_0)) / _S726;
        float _S942 = _S722 * - _S941;
        float _S943 = _S724 * _S932.differential_0 + _S930.differential_0.primal_0;
        k_2 = _S691 * _S941;
        _S691 = _S942;
        _S692 = 0.0f;
        k_3 = 0.0f;
        _S723 = _S943;
        _S724 = _S940;
        _S725 = _S931.differential_0.primal_0;
    }
    float2  _S944 = _S700 * _S917;
    float2  _S945 = _S719 * _S917;
    float2  _S946 = _S659;
    *&((&_S946)->y) = k_2;
    *&((&_S946)->x) = k_2;
    float2  _S947 = _S719 * _S946;
    float2  _S948 = _S944 + _S661 * _S946;
    float _S949 = _S948.x;
    float _S950 = _S949 + _S949;
    float _S951 = _S716 * _S950;
    float _S952 = _S948.y + _S948.y;
    float _S953 = _S716 * _S952;
    float _S954 = u_7 * _S950 + v_7 * _S952;
    float _S955 = k4_2 * _S954;
    float _S956 = _S715 * _S954;
    float _S957 = _S714 * _S954;
    float _S958 = _S714 * _S955;
    float _S959 = _S713 * _S954;
    float _S960 = _S702 * _S954 + r2_7 * _S955;
    float _S961 = _S713 * _S960;
    float _S962 = _S712 * _S954;
    float _S963 = _S703 * _S954 + r2_7 * _S960;
    float _S964 = _S712 * _S963;
    float _S965 = _S704 * _S954 + r2_7 * _S963;
    float _S966 = _S678 * _S948.x;
    float _S967 = _S711 * _S948.x;
    float _S968 = v_7 * _S966;
    float _S969 = _S654.primal_0 * _S966;
    float _S970 = _S706 * _S948.y;
    float _S971 = _S654.primal_0 * _S948.y;
    float _S972 = 2.0f * _S948.x;
    float _S973 = _S710 * _S972;
    float _S974 = _S710 * _S948.x;
    float _S975 = _S954 + u_7 * _S972 + _S707 * _S948.x;
    float _S976 = p2_2 * _S975;
    float _S977 = _S654.primal_0 * _S975;
    float _S978 = sx1_2 * _S954;
    float _S979 = _S654.primal_0 * _S954;
    float2  _S980 = _S705 * _S948;
    float2  _S981 = _S708 * _S948;
    float2  _S982 = _S659;
    *&((&_S982)->y) = _S965;
    *&((&_S982)->x) = _S965;
    float _S983 = _S981.x + _S981.y;
    float _S984 = _S962 + r2_7 * _S983;
    float _S985 = _S959 + r2_7 * _S984;
    float _S986 = _S957 + r2_7 * _S985;
    float _S987 = _S958 + _S961 + _S964 + _S704 * _S983 + _S703 * _S984 + _S702 * _S985 + k4_2 * _S986;
    float _S988 = v_7 * _S987;
    float _S989 = u_7 * _S987;
    float2  _S990 = _S947 + _S923.differential_0.primal_0;
    float _S991 = _S967 + u_7 * _S971;
    float _S992 = _S984 + _S888;
    float _S993 = _S986 + _S890;
    float _S994 = _S968 + _S970 + _S976 + _S978 + (_S980 + _S701 * _S982).x;
    float2  _S995 = _S708 * _S982 + make_float2 (_S951 + _S973 + 2.0f * _S974 + _S678 * _S971 + _S989 + _S989, _S953 + _S969 + _S988 + _S988);
    float _S996 = _S985 + _S889;
    float _S997 = _S956 + r2_7 * _S986 + _S898;
    float2  _S998 = _S661 * _S995;
    float _S999 = _S945.x + _S945.y + _S998.x + _S998.y;
    float2  _S1000 = _S699 * _S995 + _S990;
    if(_S690)
    {
        float _S1001 = _S999 / _S695;
        float _S1002 = _S689 * (0.3333333432674408f * - (k_3 + _S663 * _S1001));
        float _S1003 = _S663 * _S692 + _S663 * _S692 + _S693 * - _S1001 + _S725;
        k_2 = _S1002 + _S1002 + _S724;
        _S691 = _S1003;
        _S692 = _S723;
    }
    else
    {
        float _S1004 = _S999 / _S694;
        float _S1005 = _S662 * _S691 + _S662 * _S691 + _S689 * - _S1004 + _S723;
        k_2 = _S662 * _S1004 + _S724;
        _S691 = _S725;
        _S692 = _S1005;
    }
    DiffPair_float_0 _S1006;
    (&_S1006)->primal_0 = _S662;
    (&_S1006)->differential_0 = 0.0f;
    DiffPair_float_0 _S1007;
    (&_S1007)->primal_0 = _S663;
    (&_S1007)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1006, &_S1007, k_2);
    float _S1008 = _S1007.differential_0 + _S691;
    float _S1009 = _S1006.differential_0 + _S692;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1010;
    (&_S1010)->primal_0 = _S661;
    (&_S1010)->differential_0 = _S659;
    s_bwd_length_impl_0(&_S1010, _S1009);
    float2  _S1011 = _S1010.differential_0 + _S1000;
    float3  _S1012 = make_float3 (_S1011.x, _S1011.y, _S1008);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1013;
    (&_S1013)->primal_0 = _S661;
    (&_S1013)->differential_0 = _S659;
    s_bwd_length_impl_0(&_S1013, 0.0f);
    float _S1014 = _S655.primal_0 * dpmean2d_1.y;
    float _S1015 = _S654.primal_0 * dpmean2d_1.x;
    float2  _S1016 = make_float2 (_S1015, _S1014);
    float2  _S1017 = _S673 * _S1016;
    float _S1018 = p1_2 * _S1014;
    float _S1019 = v_6 * _S1014;
    float _S1020 = p2_2 * _S1015;
    float _S1021 = v_6 * _S1015;
    float _S1022 = _S1017.x + _S1017.y;
    float _S1023 = r2_6 * _S1022;
    float _S1024 = r2_6 * _S1023;
    float _S1025 = r2_6 * _S1024;
    float _S1026 = sy1_2 * _S1014 + _S1018 + sx1_2 * _S1015 + _S1020 + _S676 * _S1022 + _S675 * _S1023 + _S674 * _S1024 + k4_2 * _S1025;
    float _S1027 = v_6 * _S1026;
    float _S1028 = u_6 * _S1026;
    float2  _S1029 = make_float2 (r2_6 * _S1015 + _S979, r2_6 * _S1014 + _S883);
    float2  _S1030 = make_float2 (_S685 * _S1014 + 2.0f * (u_6 * _S1021 + _S991) + _S881, 2.0f * (u_6 * _S1019 + _S897) + _S681 * _S1015 + _S977);
    float4  _S1031 = make_float4 (_S1023 + _S992, _S1024 + _S996, _S1025 + _S993, r2_6 * _S1025 + _S997);
    float3  _S1032 = _S1012 + make_float3 (_S1013.differential_0.x, _S1013.differential_0.y, 0.0f) + _S915;
    float _S1033 = _S687 * dpmean2d_1.x + _S994;
    float _S1034 = _S688 * dpmean2d_1.y + _S896;
    float2  _S1035 = _S677 * _S1016 + make_float2 (_S682 * _S1019 + _S680 * _S1020 + 2.0f * (u_6 * _S1020) + _S678 * _S1021 + _S1028 + _S1028, _S684 * _S1018 + 2.0f * (v_6 * _S1018) + _S683 * _S1014 + _S679 * _S1015 + _S1027 + _S1027);
    CameraDistortion_0 _S1036 = _S810;
    (&_S1036)->thin_prism_coeffs_0 = _S1029;
    (&_S1036)->tangential_coeffs_0 = _S1030;
    (&_S1036)->radial_coeffs_0 = _S1031;
    CameraDistortion_0 _S1037 = CameraDistortion_x24_syn_dadd_0(_S810, _S1036);
    float2  _S1038 = _S661 * _S1035;
    float2  _S1039 = _S672 * _S1035;
    float _S1040 = _S1038.x + _S1038.y;
    if(_S665)
    {
        float _S1041 = _S1040 / _S667;
        float _S1042 = _S668 * - _S1041;
        float _S1043 = _S664 * (0.3333333432674408f * - (_S663 * _S1041));
        k_2 = _S1043 + _S1043;
        _S666 = _S1042;
        _S667 = 0.0f;
    }
    else
    {
        float _S1044 = _S1040 / _S666;
        float _S1045 = _S664 * - _S1044;
        k_2 = _S662 * _S1044;
        _S666 = 0.0f;
        _S667 = _S1045;
    }
    DiffPair_float_0 _S1046;
    (&_S1046)->primal_0 = _S662;
    (&_S1046)->differential_0 = 0.0f;
    DiffPair_float_0 _S1047;
    (&_S1047)->primal_0 = _S663;
    (&_S1047)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1046, &_S1047, k_2);
    float _S1048 = _S1047.differential_0 + _S666;
    float _S1049 = _S1046.differential_0 + _S667;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1050;
    (&_S1050)->primal_0 = _S661;
    (&_S1050)->differential_0 = _S659;
    s_bwd_length_impl_0(&_S1050, _S1049);
    float2  _S1051 = _S1050.differential_0 + _S1039;
    float4  _S1052 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1053;
    (&_S1053)->primal_0 = _S656.primal_0;
    (&_S1053)->differential_0 = _S1052;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1054;
    (&_S1054)->primal_0 = _S657.primal_0;
    (&_S1054)->differential_0 = _S659;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1055;
    (&_S1055)->primal_0 = _S658.primal_0;
    (&_S1055)->differential_0 = _S659;
    s_bwd_prop_CameraDistortion_x24init_0(&_S1053, &_S1054, &_S1055, _S1037);
    dpthin_prism_coeffs_3->primal_0 = (*dpthin_prism_coeffs_3).primal_0;
    dpthin_prism_coeffs_3->differential_0 = _S1055.differential_0;
    dptangential_coeffs_3->primal_0 = (*dptangential_coeffs_3).primal_0;
    dptangential_coeffs_3->differential_0 = _S1054.differential_0;
    dpradial_coeffs_3->primal_0 = (*dpradial_coeffs_3).primal_0;
    dpradial_coeffs_3->differential_0 = _S1053.differential_0;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = dpmean2d_1.y;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = dpmean2d_1.x;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S1034;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S1033;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S818.differential_0;
    float3  _S1056 = _S1032 + make_float3 (_S1051.x, _S1051.y, _S1048);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1056;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_11, float fy_11, float cx_11, float cy_11, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, uint image_width_7, uint image_height_7, float v_depth_1, float2  v_mean2d_1, float3  v_conic_1, float v_opacity_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    Matrix<float, 2, 2>  _S1057 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1058 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1059 = { _S1058, _S1058 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1060 = { _S1058, _S1058, _S1059, _S1058, _S1058, _S1059 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1061;
    (&_S1061)->_S409 = _S1057;
    (&_S1061)->_S410 = _S1060;
    (&_S1061)->_S409 = _S1057;
    float3  mean_c_7 = s_primal_ctx_mul_0(R_9, mean_7) + t_8;
    float3  _S1062 = s_primal_ctx_exp_0(scale_9);
    float _S1063 = quat_10.y;
    float _S1064 = _S1063 * _S1063 + quat_10.z * quat_10.z + quat_10.w * quat_10.w + quat_10.x * quat_10.x;
    float _S1065 = s_primal_ctx_rsqrt_0(_S1064);
    float x_27 = quat_10.y * _S1065;
    float y_12 = quat_10.z * _S1065;
    float z_10 = quat_10.w * _S1065;
    float w_10 = quat_10.x * _S1065;
    float x2_10 = x_27 * x_27;
    float y2_10 = y_12 * y_12;
    float z2_10 = z_10 * z_10;
    float xy_10 = x_27 * y_12;
    float xz_10 = x_27 * z_10;
    float yz_10 = y_12 * z_10;
    float wx_10 = w_10 * x_27;
    float wy_10 = w_10 * y_12;
    float wz_10 = w_10 * z_10;
    Matrix<float, 3, 3>  _S1066 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_10), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_10), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1062.x, 0.0f, 0.0f, 0.0f, _S1062.y, 0.0f, 0.0f, 0.0f, _S1062.z);
    Matrix<float, 3, 3>  _S1067 = s_primal_ctx_mul_1(_S1066, S_1);
    Matrix<float, 3, 3>  _S1068 = transpose_0(_S1067);
    Matrix<float, 3, 3>  _S1069 = s_primal_ctx_mul_1(_S1067, _S1068);
    Matrix<float, 3, 3>  _S1070 = s_primal_ctx_mul_1(R_9, _S1069);
    Matrix<float, 3, 3>  _S1071 = transpose_0(R_9);
    Matrix<float, 3, 3>  _S1072 = s_primal_ctx_mul_1(_S1070, _S1071);
    Matrix<float, 2, 2>  _S1073 = _S1057;
    float2  _S1074 = make_float2 (0.0f);
    float2  _S1075 = _S1074;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_7, _S1072, fx_11, fy_11, cx_11, cy_11, radial_coeffs_10, tangential_coeffs_10, thin_prism_coeffs_10, &_S1073, &_S1075, &(&_S1061)->_S410);
    (&_S1061)->_S409 = _S1073;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1076 = _S1061;
    float _S1077 = _S1061._S409.rows[int(0)].y * _S1061._S409.rows[int(1)].x;
    float det_orig_8 = _S1061._S409.rows[int(0)].x * _S1061._S409.rows[int(1)].y - _S1077;
    float _S1078 = _S1061._S409.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1079 = _S1061._S409;
    *&(((&_S1079)->rows + (int(0)))->x) = _S1078;
    float _S1080 = _S1061._S409.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1079)->rows + (int(1)))->y) = _S1080;
    Matrix<float, 2, 2>  _S1081 = _S1079;
    Matrix<float, 2, 2>  _S1082 = _S1079;
    float det_blur_5 = _S1078 * _S1080 - _S1077;
    float _S1083 = det_orig_8 / det_blur_5;
    float _S1084 = det_blur_5 * det_blur_5;
    float _S1085 = s_primal_ctx_max_0(0.0f, _S1083);
    float _S1086 = s_primal_ctx_sqrt_0(_S1085);
    float invdet_8 = 1.0f / det_blur_5;
    float _S1087 = - _S1061._S409.rows[int(0)].y;
    float _S1088 = - _S1061._S409.rows[int(1)].x;
    float _S1089 = - in_opacity_7;
    float _S1090 = 1.0f + s_primal_ctx_exp_1(_S1089);
    float _S1091 = 1.0f / _S1090;
    float _S1092 = _S1090 * _S1090;
    float _S1093;
    if(antialiased_7)
    {
        _S1093 = _S1091 * _S1086;
    }
    else
    {
        _S1093 = _S1091;
    }
    float _S1094 = _S1093 / 0.00392156885936856f;
    float _S1095 = 2.0f * s_primal_ctx_log_0(_S1094);
    float _S1096 = s_primal_ctx_sqrt_0(_S1095);
    float _S1097 = _S1081.rows[int(0)].x;
    float _S1098 = _S1082.rows[int(1)].y;
    float2  _S1099 = _S1074;
    *&((&_S1099)->y) = v_conic_1.z;
    float2  _S1100 = _S1074;
    *&((&_S1100)->y) = v_conic_1.y;
    *&((&_S1100)->x) = v_conic_1.x;
    DiffPair_float_0 _S1101;
    (&_S1101)->primal_0 = _S1098;
    (&_S1101)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1101, 0.0f);
    DiffPair_float_0 _S1102;
    (&_S1102)->primal_0 = _S1097;
    (&_S1102)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1102, 0.0f);
    DiffPair_float_0 _S1103;
    (&_S1103)->primal_0 = 3.32999992370605469f;
    (&_S1103)->differential_0 = 0.0f;
    DiffPair_float_0 _S1104;
    (&_S1104)->primal_0 = _S1096;
    (&_S1104)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1103, &_S1104, 0.0f);
    DiffPair_float_0 _S1105;
    (&_S1105)->primal_0 = _S1095;
    (&_S1105)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1105, _S1104.differential_0);
    float _S1106 = 2.0f * _S1105.differential_0;
    DiffPair_float_0 _S1107;
    (&_S1107)->primal_0 = _S1094;
    (&_S1107)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1107, _S1106);
    float _S1108 = v_opacity_1 + 254.9999847412109375f * _S1107.differential_0;
    Matrix<float, 2, 2>  _S1109 = _S1057;
    _S1109[int(1)] = _S1099;
    _S1109[int(0)] = _S1100;
    Matrix<float, 2, 2>  _S1110 = _S1109;
    float2  _S1111 = make_float2 (_S1102.differential_0, 0.0f);
    float3  _S1112 = make_float3 (0.0f, 0.0f, v_depth_1);
    float2  _S1113 = make_float2 (0.0f, _S1101.differential_0);
    float _S1114;
    if(antialiased_7)
    {
        float _S1115 = _S1091 * _S1108;
        _S1093 = _S1086 * _S1108;
        _S1114 = _S1115;
    }
    else
    {
        _S1093 = _S1108;
        _S1114 = 0.0f;
    }
    float _S1116 = - (_S1093 / _S1092);
    DiffPair_float_0 _S1117;
    (&_S1117)->primal_0 = _S1089;
    (&_S1117)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1117, _S1116);
    float _S1118 = - _S1117.differential_0;
    float _S1119 = invdet_8 * _S1110.rows[int(1)].y;
    float _S1120 = - (invdet_8 * _S1110.rows[int(1)].x);
    float _S1121 = - (invdet_8 * _S1110.rows[int(0)].y);
    float _S1122 = invdet_8 * _S1110.rows[int(0)].x;
    float _S1123 = - ((_S1078 * _S1110.rows[int(1)].y + _S1088 * _S1110.rows[int(1)].x + _S1087 * _S1110.rows[int(0)].y + _S1080 * _S1110.rows[int(0)].x) / _S1084);
    DiffPair_float_0 _S1124;
    (&_S1124)->primal_0 = _S1085;
    (&_S1124)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1124, _S1114);
    DiffPair_float_0 _S1125;
    (&_S1125)->primal_0 = 0.0f;
    (&_S1125)->differential_0 = 0.0f;
    DiffPair_float_0 _S1126;
    (&_S1126)->primal_0 = _S1083;
    (&_S1126)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1125, &_S1126, _S1124.differential_0);
    float _S1127 = _S1126.differential_0 / _S1084;
    float s_diff_det_orig_T_1 = det_blur_5 * _S1127;
    float _S1128 = _S1123 + det_orig_8 * - _S1127;
    float _S1129 = - _S1128;
    float _S1130 = _S1078 * _S1128;
    float _S1131 = _S1080 * _S1128;
    Matrix<float, 2, 2>  _S1132 = _S1057;
    _S1132[int(1)] = _S1113;
    _S1132[int(0)] = _S1111;
    _S1079 = _S1132;
    *&(((&_S1079)->rows + (int(1)))->y) = 0.0f;
    float _S1133 = _S1122 + _S1130 + _S1132.rows[int(1)].y;
    *&(((&_S1079)->rows + (int(0)))->x) = 0.0f;
    float _S1134 = _S1119 + _S1131 + _S1132.rows[int(0)].x;
    float _S1135 = _S1129 + - s_diff_det_orig_T_1;
    float _S1136 = _S1120 + _S1076._S409.rows[int(0)].y * _S1135;
    float _S1137 = _S1121 + _S1076._S409.rows[int(1)].x * _S1135;
    float _S1138 = _S1076._S409.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S1139 = _S1133 + _S1076._S409.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S1140 = _S1074;
    *&((&_S1140)->x) = _S1136;
    *&((&_S1140)->y) = _S1139;
    float _S1141 = _S1134 + _S1138;
    float2  _S1142 = _S1074;
    *&((&_S1142)->y) = _S1137;
    *&((&_S1142)->x) = _S1141;
    Matrix<float, 2, 2>  _S1143 = _S1057;
    _S1143[int(1)] = _S1140;
    _S1143[int(0)] = _S1142;
    Matrix<float, 2, 2>  _S1144 = _S1079 + _S1143;
    float3  _S1145 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1146;
    (&_S1146)->primal_0 = mean_c_7;
    (&_S1146)->differential_0 = _S1145;
    Matrix<float, 3, 3>  _S1147 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1148;
    (&_S1148)->primal_0 = _S1072;
    (&_S1148)->differential_0 = _S1147;
    DiffPair_float_0 _S1149;
    (&_S1149)->primal_0 = fx_11;
    (&_S1149)->differential_0 = 0.0f;
    DiffPair_float_0 _S1150;
    (&_S1150)->primal_0 = fy_11;
    (&_S1150)->differential_0 = 0.0f;
    DiffPair_float_0 _S1151;
    (&_S1151)->primal_0 = cx_11;
    (&_S1151)->differential_0 = 0.0f;
    DiffPair_float_0 _S1152;
    (&_S1152)->primal_0 = cy_11;
    (&_S1152)->differential_0 = 0.0f;
    float4  _S1153 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1154;
    (&_S1154)->primal_0 = radial_coeffs_10;
    (&_S1154)->differential_0 = _S1153;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1155;
    (&_S1155)->primal_0 = tangential_coeffs_10;
    (&_S1155)->differential_0 = _S1074;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1156;
    (&_S1156)->primal_0 = thin_prism_coeffs_10;
    (&_S1156)->differential_0 = _S1074;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S1146, &_S1148, &_S1149, &_S1150, &_S1151, &_S1152, &_S1154, &_S1155, &_S1156, _S1144, v_mean2d_1, _S1076._S410);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1157;
    (&_S1157)->primal_0 = _S1070;
    (&_S1157)->differential_0 = _S1147;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1158;
    (&_S1158)->primal_0 = _S1071;
    (&_S1158)->differential_0 = _S1147;
    s_bwd_prop_mul_2(&_S1157, &_S1158, _S1148.differential_0);
    Matrix<float, 3, 3>  _S1159 = transpose_0(_S1158.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1160;
    (&_S1160)->primal_0 = R_9;
    (&_S1160)->differential_0 = _S1147;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1161;
    (&_S1161)->primal_0 = _S1069;
    (&_S1161)->differential_0 = _S1147;
    s_bwd_prop_mul_2(&_S1160, &_S1161, _S1157.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1162;
    (&_S1162)->primal_0 = _S1067;
    (&_S1162)->differential_0 = _S1147;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1163;
    (&_S1163)->primal_0 = _S1068;
    (&_S1163)->differential_0 = _S1147;
    s_bwd_prop_mul_2(&_S1162, &_S1163, _S1161.differential_0);
    Matrix<float, 3, 3>  _S1164 = _S1162.differential_0 + transpose_0(_S1163.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1165;
    (&_S1165)->primal_0 = _S1066;
    (&_S1165)->differential_0 = _S1147;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1166;
    (&_S1166)->primal_0 = S_1;
    (&_S1166)->differential_0 = _S1147;
    s_bwd_prop_mul_2(&_S1165, &_S1166, _S1164);
    Matrix<float, 3, 3>  _S1167 = transpose_0(_S1165.differential_0);
    float _S1168 = 2.0f * - _S1167.rows[int(2)].z;
    float _S1169 = 2.0f * _S1167.rows[int(2)].y;
    float _S1170 = 2.0f * _S1167.rows[int(2)].x;
    float _S1171 = 2.0f * _S1167.rows[int(1)].z;
    float _S1172 = 2.0f * - _S1167.rows[int(1)].y;
    float _S1173 = 2.0f * _S1167.rows[int(1)].x;
    float _S1174 = 2.0f * _S1167.rows[int(0)].z;
    float _S1175 = 2.0f * _S1167.rows[int(0)].y;
    float _S1176 = 2.0f * - _S1167.rows[int(0)].x;
    float _S1177 = - _S1173 + _S1175;
    float _S1178 = _S1170 + - _S1174;
    float _S1179 = - _S1169 + _S1171;
    float _S1180 = _S1169 + _S1171;
    float _S1181 = _S1170 + _S1174;
    float _S1182 = _S1173 + _S1175;
    float _S1183 = z_10 * (_S1172 + _S1176);
    float _S1184 = y_12 * (_S1168 + _S1176);
    float _S1185 = x_27 * (_S1168 + _S1172);
    float _S1186 = z_10 * _S1177 + y_12 * _S1178 + x_27 * _S1179;
    float _S1187 = _S1065 * _S1186;
    float _S1188 = w_10 * _S1177 + y_12 * _S1180 + x_27 * _S1181 + _S1183 + _S1183;
    float _S1189 = _S1065 * _S1188;
    float _S1190 = w_10 * _S1178 + z_10 * _S1180 + x_27 * _S1182 + _S1184 + _S1184;
    float _S1191 = _S1065 * _S1190;
    float _S1192 = w_10 * _S1179 + z_10 * _S1181 + y_12 * _S1182 + _S1185 + _S1185;
    float _S1193 = _S1065 * _S1192;
    float _S1194 = quat_10.x * _S1186 + quat_10.w * _S1188 + quat_10.z * _S1190 + quat_10.y * _S1192;
    DiffPair_float_0 _S1195;
    (&_S1195)->primal_0 = _S1064;
    (&_S1195)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1195, _S1194);
    float _S1196 = quat_10.x * _S1195.differential_0;
    float _S1197 = quat_10.w * _S1195.differential_0;
    float _S1198 = quat_10.z * _S1195.differential_0;
    float _S1199 = quat_10.y * _S1195.differential_0;
    float _S1200 = _S1189 + _S1197 + _S1197;
    float _S1201 = _S1191 + _S1198 + _S1198;
    float _S1202 = _S1193 + _S1199 + _S1199;
    float _S1203 = _S1187 + _S1196 + _S1196;
    float3  _S1204 = _S1145;
    *&((&_S1204)->z) = _S1166.differential_0.rows[int(2)].z;
    *&((&_S1204)->y) = _S1166.differential_0.rows[int(1)].y;
    *&((&_S1204)->x) = _S1166.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1205;
    (&_S1205)->primal_0 = scale_9;
    (&_S1205)->differential_0 = _S1145;
    s_bwd_prop_exp_1(&_S1205, _S1204);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1206 = _S1205;
    float3  _S1207 = _S1146.differential_0 + _S1112;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1208;
    (&_S1208)->primal_0 = R_9;
    (&_S1208)->differential_0 = _S1147;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1209;
    (&_S1209)->primal_0 = mean_7;
    (&_S1209)->differential_0 = _S1145;
    s_bwd_prop_mul_3(&_S1208, &_S1209, _S1207);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1210 = _S1209;
    Matrix<float, 3, 3>  _S1211 = _S1159 + _S1160.differential_0 + _S1208.differential_0;
    float4  _S1212 = _S1153;
    *&((&_S1212)->w) = _S1200;
    *&((&_S1212)->z) = _S1201;
    *&((&_S1212)->y) = _S1202;
    *&((&_S1212)->x) = _S1203;
    float4  _S1213 = _S1212;
    *v_mean_1 = _S1210.differential_0;
    *v_quat_1 = _S1213;
    *v_scale_1 = _S1206.differential_0;
    *v_in_opacity_1 = _S1118;
    *v_R_1 = _S1211;
    *v_t_1 = _S1207;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_12, float fy_12, float cx_12, float cy_12, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, uint image_width_8, uint image_height_8, float v_depth_2, float2  v_mean2d_2, float3  v_conic_2, float v_opacity_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    float3  _S1214 = s_primal_ctx_exp_0(scale_10);
    float _S1215 = quat_11.y;
    float _S1216 = _S1215 * _S1215 + quat_11.z * quat_11.z + quat_11.w * quat_11.w + quat_11.x * quat_11.x;
    float _S1217 = s_primal_ctx_rsqrt_0(_S1216);
    float x_28 = quat_11.y * _S1217;
    float y_13 = quat_11.z * _S1217;
    float z_11 = quat_11.w * _S1217;
    float w_11 = quat_11.x * _S1217;
    float x2_11 = x_28 * x_28;
    float y2_11 = y_13 * y_13;
    float z2_11 = z_11 * z_11;
    float xy_11 = x_28 * y_13;
    float xz_11 = x_28 * z_11;
    float yz_11 = y_13 * z_11;
    float wx_11 = w_11 * x_28;
    float wy_11 = w_11 * y_13;
    float wz_11 = w_11 * z_11;
    Matrix<float, 3, 3>  _S1218 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_11), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_11), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S1214.x, 0.0f, 0.0f, 0.0f, _S1214.y, 0.0f, 0.0f, 0.0f, _S1214.z);
    Matrix<float, 3, 3>  _S1219 = s_primal_ctx_mul_1(_S1218, S_2);
    Matrix<float, 3, 3>  _S1220 = transpose_0(_S1219);
    Matrix<float, 3, 3>  _S1221 = s_primal_ctx_mul_1(_S1219, _S1220);
    Matrix<float, 3, 3>  _S1222 = s_primal_ctx_mul_1(R_10, _S1221);
    Matrix<float, 3, 3>  _S1223 = transpose_0(R_10);
    Matrix<float, 3, 3>  _S1224 = s_primal_ctx_mul_1(_S1222, _S1223);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_12, 0.0f, 0.0f, 0.0f, fy_12, 0.0f);
    Matrix<float, 2, 3>  _S1225 = s_primal_ctx_mul_2(J_11, _S1224);
    Matrix<float, 3, 2>  _S1226 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S1227 = s_primal_ctx_mul_3(_S1225, _S1226);
    float _S1228 = _S1227.rows[int(0)].y * _S1227.rows[int(1)].x;
    float det_orig_9 = _S1227.rows[int(0)].x * _S1227.rows[int(1)].y - _S1228;
    float _S1229 = _S1227.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1230 = _S1227;
    *&(((&_S1230)->rows + (int(0)))->x) = _S1229;
    float _S1231 = _S1227.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1230)->rows + (int(1)))->y) = _S1231;
    Matrix<float, 2, 2>  _S1232 = _S1230;
    Matrix<float, 2, 2>  _S1233 = _S1230;
    float det_blur_6 = _S1229 * _S1231 - _S1228;
    float _S1234 = det_orig_9 / det_blur_6;
    float _S1235 = det_blur_6 * det_blur_6;
    float _S1236 = s_primal_ctx_max_0(0.0f, _S1234);
    float _S1237 = s_primal_ctx_sqrt_0(_S1236);
    float invdet_9 = 1.0f / det_blur_6;
    float _S1238 = - _S1227.rows[int(0)].y;
    float _S1239 = - _S1227.rows[int(1)].x;
    float _S1240 = - in_opacity_8;
    float _S1241 = 1.0f + s_primal_ctx_exp_1(_S1240);
    float _S1242 = 1.0f / _S1241;
    float _S1243 = _S1241 * _S1241;
    float _S1244;
    if(antialiased_8)
    {
        _S1244 = _S1242 * _S1237;
    }
    else
    {
        _S1244 = _S1242;
    }
    float _S1245 = _S1244 / 0.00392156885936856f;
    float _S1246 = 2.0f * s_primal_ctx_log_0(_S1245);
    float _S1247 = s_primal_ctx_sqrt_0(_S1246);
    float _S1248 = _S1232.rows[int(0)].x;
    float _S1249 = _S1233.rows[int(1)].y;
    float2  _S1250 = make_float2 (0.0f);
    float2  _S1251 = _S1250;
    *&((&_S1251)->y) = v_conic_2.z;
    float2  _S1252 = _S1250;
    *&((&_S1252)->y) = v_conic_2.y;
    *&((&_S1252)->x) = v_conic_2.x;
    DiffPair_float_0 _S1253;
    (&_S1253)->primal_0 = _S1249;
    (&_S1253)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1253, 0.0f);
    DiffPair_float_0 _S1254;
    (&_S1254)->primal_0 = _S1248;
    (&_S1254)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1254, 0.0f);
    DiffPair_float_0 _S1255;
    (&_S1255)->primal_0 = 3.32999992370605469f;
    (&_S1255)->differential_0 = 0.0f;
    DiffPair_float_0 _S1256;
    (&_S1256)->primal_0 = _S1247;
    (&_S1256)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1255, &_S1256, 0.0f);
    DiffPair_float_0 _S1257;
    (&_S1257)->primal_0 = _S1246;
    (&_S1257)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1257, _S1256.differential_0);
    float _S1258 = 2.0f * _S1257.differential_0;
    DiffPair_float_0 _S1259;
    (&_S1259)->primal_0 = _S1245;
    (&_S1259)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1259, _S1258);
    float _S1260 = v_opacity_2 + 254.9999847412109375f * _S1259.differential_0;
    Matrix<float, 2, 2>  _S1261 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1262 = _S1261;
    _S1262[int(1)] = _S1251;
    _S1262[int(0)] = _S1252;
    Matrix<float, 2, 2>  _S1263 = _S1262;
    float3  _S1264 = make_float3 (0.0f, 0.0f, v_depth_2);
    float2  _S1265 = make_float2 (0.0f, _S1253.differential_0);
    float2  _S1266 = make_float2 (_S1254.differential_0, 0.0f);
    float _S1267;
    if(antialiased_8)
    {
        float _S1268 = _S1242 * _S1260;
        _S1244 = _S1237 * _S1260;
        _S1267 = _S1268;
    }
    else
    {
        _S1244 = _S1260;
        _S1267 = 0.0f;
    }
    float _S1269 = - (_S1244 / _S1243);
    DiffPair_float_0 _S1270;
    (&_S1270)->primal_0 = _S1240;
    (&_S1270)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1270, _S1269);
    float _S1271 = - _S1270.differential_0;
    float _S1272 = invdet_9 * _S1263.rows[int(1)].y;
    float _S1273 = - (invdet_9 * _S1263.rows[int(1)].x);
    float _S1274 = - (invdet_9 * _S1263.rows[int(0)].y);
    float _S1275 = invdet_9 * _S1263.rows[int(0)].x;
    float _S1276 = - ((_S1229 * _S1263.rows[int(1)].y + _S1239 * _S1263.rows[int(1)].x + _S1238 * _S1263.rows[int(0)].y + _S1231 * _S1263.rows[int(0)].x) / _S1235);
    DiffPair_float_0 _S1277;
    (&_S1277)->primal_0 = _S1236;
    (&_S1277)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1277, _S1267);
    DiffPair_float_0 _S1278;
    (&_S1278)->primal_0 = 0.0f;
    (&_S1278)->differential_0 = 0.0f;
    DiffPair_float_0 _S1279;
    (&_S1279)->primal_0 = _S1234;
    (&_S1279)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1278, &_S1279, _S1277.differential_0);
    float _S1280 = _S1279.differential_0 / _S1235;
    float s_diff_det_orig_T_2 = det_blur_6 * _S1280;
    float _S1281 = _S1276 + det_orig_9 * - _S1280;
    float _S1282 = - _S1281;
    float _S1283 = _S1229 * _S1281;
    float _S1284 = _S1231 * _S1281;
    Matrix<float, 2, 2>  _S1285 = _S1261;
    _S1285[int(1)] = _S1265;
    _S1285[int(0)] = _S1266;
    _S1230 = _S1285;
    *&(((&_S1230)->rows + (int(1)))->y) = 0.0f;
    float _S1286 = _S1275 + _S1283 + _S1285.rows[int(1)].y;
    *&(((&_S1230)->rows + (int(0)))->x) = 0.0f;
    float _S1287 = _S1272 + _S1284 + _S1285.rows[int(0)].x;
    float _S1288 = _S1282 + - s_diff_det_orig_T_2;
    float _S1289 = _S1273 + _S1227.rows[int(0)].y * _S1288;
    float _S1290 = _S1274 + _S1227.rows[int(1)].x * _S1288;
    float _S1291 = _S1227.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S1292 = _S1286 + _S1227.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S1293 = _S1250;
    *&((&_S1293)->x) = _S1289;
    *&((&_S1293)->y) = _S1292;
    float _S1294 = _S1287 + _S1291;
    float2  _S1295 = _S1250;
    *&((&_S1295)->y) = _S1290;
    *&((&_S1295)->x) = _S1294;
    float _S1296 = fy_12 * v_mean2d_2.y;
    float _S1297 = fx_12 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S1298 = _S1261;
    _S1298[int(1)] = _S1293;
    _S1298[int(0)] = _S1295;
    Matrix<float, 2, 2>  _S1299 = _S1230 + _S1298;
    Matrix<float, 2, 3>  _S1300 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1301;
    (&_S1301)->primal_0 = _S1225;
    (&_S1301)->differential_0 = _S1300;
    Matrix<float, 3, 2>  _S1302 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1303;
    (&_S1303)->primal_0 = _S1226;
    (&_S1303)->differential_0 = _S1302;
    s_bwd_prop_mul_0(&_S1301, &_S1303, _S1299);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1304;
    (&_S1304)->primal_0 = J_11;
    (&_S1304)->differential_0 = _S1300;
    Matrix<float, 3, 3>  _S1305 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1306;
    (&_S1306)->primal_0 = _S1224;
    (&_S1306)->differential_0 = _S1305;
    s_bwd_prop_mul_1(&_S1304, &_S1306, _S1301.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1307;
    (&_S1307)->primal_0 = _S1222;
    (&_S1307)->differential_0 = _S1305;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1308;
    (&_S1308)->primal_0 = _S1223;
    (&_S1308)->differential_0 = _S1305;
    s_bwd_prop_mul_2(&_S1307, &_S1308, _S1306.differential_0);
    Matrix<float, 3, 3>  _S1309 = transpose_0(_S1308.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1310;
    (&_S1310)->primal_0 = R_10;
    (&_S1310)->differential_0 = _S1305;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1311;
    (&_S1311)->primal_0 = _S1221;
    (&_S1311)->differential_0 = _S1305;
    s_bwd_prop_mul_2(&_S1310, &_S1311, _S1307.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1312;
    (&_S1312)->primal_0 = _S1219;
    (&_S1312)->differential_0 = _S1305;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1313;
    (&_S1313)->primal_0 = _S1220;
    (&_S1313)->differential_0 = _S1305;
    s_bwd_prop_mul_2(&_S1312, &_S1313, _S1311.differential_0);
    Matrix<float, 3, 3>  _S1314 = _S1312.differential_0 + transpose_0(_S1313.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1315;
    (&_S1315)->primal_0 = _S1218;
    (&_S1315)->differential_0 = _S1305;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1316;
    (&_S1316)->primal_0 = S_2;
    (&_S1316)->differential_0 = _S1305;
    s_bwd_prop_mul_2(&_S1315, &_S1316, _S1314);
    Matrix<float, 3, 3>  _S1317 = transpose_0(_S1315.differential_0);
    float _S1318 = 2.0f * - _S1317.rows[int(2)].z;
    float _S1319 = 2.0f * _S1317.rows[int(2)].y;
    float _S1320 = 2.0f * _S1317.rows[int(2)].x;
    float _S1321 = 2.0f * _S1317.rows[int(1)].z;
    float _S1322 = 2.0f * - _S1317.rows[int(1)].y;
    float _S1323 = 2.0f * _S1317.rows[int(1)].x;
    float _S1324 = 2.0f * _S1317.rows[int(0)].z;
    float _S1325 = 2.0f * _S1317.rows[int(0)].y;
    float _S1326 = 2.0f * - _S1317.rows[int(0)].x;
    float _S1327 = - _S1323 + _S1325;
    float _S1328 = _S1320 + - _S1324;
    float _S1329 = - _S1319 + _S1321;
    float _S1330 = _S1319 + _S1321;
    float _S1331 = _S1320 + _S1324;
    float _S1332 = _S1323 + _S1325;
    float _S1333 = z_11 * (_S1322 + _S1326);
    float _S1334 = y_13 * (_S1318 + _S1326);
    float _S1335 = x_28 * (_S1318 + _S1322);
    float _S1336 = z_11 * _S1327 + y_13 * _S1328 + x_28 * _S1329;
    float _S1337 = _S1217 * _S1336;
    float _S1338 = w_11 * _S1327 + y_13 * _S1330 + x_28 * _S1331 + _S1333 + _S1333;
    float _S1339 = _S1217 * _S1338;
    float _S1340 = w_11 * _S1328 + z_11 * _S1330 + x_28 * _S1332 + _S1334 + _S1334;
    float _S1341 = _S1217 * _S1340;
    float _S1342 = w_11 * _S1329 + z_11 * _S1331 + y_13 * _S1332 + _S1335 + _S1335;
    float _S1343 = _S1217 * _S1342;
    float _S1344 = quat_11.x * _S1336 + quat_11.w * _S1338 + quat_11.z * _S1340 + quat_11.y * _S1342;
    DiffPair_float_0 _S1345;
    (&_S1345)->primal_0 = _S1216;
    (&_S1345)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1345, _S1344);
    float _S1346 = quat_11.x * _S1345.differential_0;
    float _S1347 = quat_11.w * _S1345.differential_0;
    float _S1348 = quat_11.z * _S1345.differential_0;
    float _S1349 = quat_11.y * _S1345.differential_0;
    float _S1350 = _S1339 + _S1347 + _S1347;
    float _S1351 = _S1341 + _S1348 + _S1348;
    float _S1352 = _S1343 + _S1349 + _S1349;
    float _S1353 = _S1337 + _S1346 + _S1346;
    float3  _S1354 = make_float3 (0.0f);
    float3  _S1355 = _S1354;
    *&((&_S1355)->z) = _S1316.differential_0.rows[int(2)].z;
    *&((&_S1355)->y) = _S1316.differential_0.rows[int(1)].y;
    *&((&_S1355)->x) = _S1316.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1356;
    (&_S1356)->primal_0 = scale_10;
    (&_S1356)->differential_0 = _S1354;
    s_bwd_prop_exp_1(&_S1356, _S1355);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1357 = _S1356;
    float3  _S1358 = _S1354;
    *&((&_S1358)->y) = _S1296;
    *&((&_S1358)->x) = _S1297;
    float3  _S1359 = _S1264 + _S1358;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1360;
    (&_S1360)->primal_0 = R_10;
    (&_S1360)->differential_0 = _S1305;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1361;
    (&_S1361)->primal_0 = mean_8;
    (&_S1361)->differential_0 = _S1354;
    s_bwd_prop_mul_3(&_S1360, &_S1361, _S1359);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1362 = _S1361;
    Matrix<float, 3, 3>  _S1363 = _S1309 + _S1310.differential_0 + _S1360.differential_0;
    float4  _S1364 = make_float4 (0.0f);
    *&((&_S1364)->w) = _S1350;
    *&((&_S1364)->z) = _S1351;
    *&((&_S1364)->y) = _S1352;
    *&((&_S1364)->x) = _S1353;
    float4  _S1365 = _S1364;
    *v_mean_2 = _S1362.differential_0;
    *v_quat_2 = _S1365;
    *v_scale_2 = _S1357.differential_0;
    *v_in_opacity_2 = _S1271;
    *v_R_2 = _S1363;
    *v_t_2 = _S1359;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_11, float dOut_14)
{
    float _S1366 = _slang_select(((*dpx_11).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_11).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_14;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S1366;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_12, float dOut_15)
{
    float _S1367 = (F32_exp2(((*dpx_12).primal_0))) * 50.693145751953125f * dOut_15;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S1367;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_16)
{
    float _S1368 = dOut_16.y;
    float _S1369 = dOut_16.z;
    float _S1370 = dOut_16.x;
    float _S1371 = (*a_0).primal_0.z * _S1368 + - (*a_0).primal_0.y * _S1369;
    float _S1372 = - (*a_0).primal_0.z * _S1370 + (*a_0).primal_0.x * _S1369;
    float _S1373 = (*a_0).primal_0.y * _S1370 + - (*a_0).primal_0.x * _S1368;
    float3  _S1374 = make_float3 (- (*b_0).primal_0.z * _S1368 + (*b_0).primal_0.y * _S1369, (*b_0).primal_0.z * _S1370 + - (*b_0).primal_0.x * _S1369, - (*b_0).primal_0.y * _S1370 + (*b_0).primal_0.x * _S1368);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S1374;
    float3  _S1375 = make_float3 (_S1371, _S1372, _S1373);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S1375;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S1376 = left_8.y;
    float _S1377 = right_8.z;
    float _S1378 = left_8.z;
    float _S1379 = right_8.y;
    float _S1380 = right_8.x;
    float _S1381 = left_8.x;
    return make_float3 (_S1376 * _S1377 - _S1378 * _S1379, _S1378 * _S1380 - _S1381 * _S1377, _S1381 * _S1379 - _S1376 * _S1380);
}

inline __device__ float3  normalize_0(float3  x_29)
{
    return x_29 / make_float3 (length_1(x_29));
}

inline __device__ float2  normalize_1(float2  x_30)
{
    return x_30 / make_float2 (length_0(x_30));
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_9, float4  quat_12, float3  scale_11, float2  hardness_0, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_13, float fy_13, float cx_13, float cy_13, float4  radial_coeffs_12, float2  tangential_coeffs_12, float2  thin_prism_coeffs_12, uint image_width_9, uint image_height_9, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float * depth_6, float3  * normal_0, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float2  * out_hardness_0)
{
    for(;;)
    {
        float3  mean_c_8 = mul_0(R_11, mean_9) + t_10;
        float _S1382 = mean_c_8.z;
        bool _S1383;
        if(_S1382 < near_plane_6)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = _S1382 > far_plane_6;
        }
        if(_S1383)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S1384 = scale_11.x;
        float sx_0 = (F32_exp((_S1384)));
        float _S1385 = scale_11.y;
        float sy_0 = (F32_exp((_S1385)));
        float sz_0 = scale_11.z - 0.5f * (_S1384 + _S1385);
        float x_31 = quat_12.y;
        float inv_norm_9 = (F32_rsqrt((x_31 * x_31 + quat_12.z * quat_12.z + quat_12.w * quat_12.w + quat_12.x * quat_12.x)));
        float x_32 = quat_12.y * inv_norm_9;
        float y_14 = quat_12.z * inv_norm_9;
        float z_12 = quat_12.w * inv_norm_9;
        float w_12 = quat_12.x * inv_norm_9;
        float x2_12 = x_32 * x_32;
        float y2_12 = y_14 * y_14;
        float z2_12 = z_12 * z_12;
        float xy_12 = x_32 * y_14;
        float xz_12 = x_32 * z_12;
        float yz_12 = y_14 * z_12;
        float wx_12 = w_12 * x_32;
        float wy_12 = w_12 * y_14;
        float wz_12 = w_12 * z_12;
        Matrix<float, 3, 3>  _S1386 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_12), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_12), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12)));
        float3  vert0_0 = mul_0(_S1386, make_float3 (sx_0, 0.0f, 0.0f)) + mean_9;
        float3  vert1_0 = mul_0(_S1386, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_9;
        float3  vert2_0 = mul_0(_S1386, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_9;
        float3  vert0_c_0 = mul_0(R_11, vert0_0) + t_10;
        float3  vert1_c_0 = mul_0(R_11, vert1_0) + t_10;
        float3  vert2_c_0 = mul_0(R_11, vert2_0) + t_10;
        float _S1387 = vert0_c_0.z;
        if(_S1387 < near_plane_6)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = _S1387 > far_plane_6;
        }
        if(_S1383)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = (vert1_c_0.z) < near_plane_6;
        }
        if(_S1383)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = (vert1_c_0.z) > far_plane_6;
        }
        if(_S1383)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = (vert2_c_0.z) < near_plane_6;
        }
        if(_S1383)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = (vert2_c_0.z) > far_plane_6;
        }
        if(_S1383)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S1387);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (vert1_c_0.z);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (vert2_c_0.z);
        float2  _S1388 = make_float2 (fx_13, fy_13);
        float2  _S1389 = make_float2 (cx_13, cy_13);
        *uv0_0 = _S1388 * *uv0_0 + _S1389;
        *uv1_0 = _S1388 * *uv1_0 + _S1389;
        float2  _S1390 = _S1388 * *uv2_0 + _S1389;
        *uv2_0 = _S1390;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S1390 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S1390)));
        float _S1391 = _S1390.x;
        float xmax_3 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S1391))) + offset_0;
        float xmin_3 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S1391))) - offset_0;
        float _S1392 = _S1390.y;
        float ymax_3 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S1392))) + offset_0;
        float ymin_3 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S1392))) - offset_0;
        if(xmax_3 <= 0.0f)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = xmin_3 >= float(image_width_9);
        }
        if(_S1383)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = ymax_3 <= 0.0f;
        }
        if(_S1383)
        {
            _S1383 = true;
        }
        else
        {
            _S1383 = ymin_3 >= float(image_height_9);
        }
        if(_S1383)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_6 = make_int4 (int((F32_floor((xmin_3)))), int((F32_floor((ymin_3)))), int((F32_ceil((xmax_3)))), int((F32_ceil((ymax_3)))));
        *depth_6 = _S1382;
        float3  _S1393 = normalize_0(mul_0(R_11, cross_0(vert1_0 - vert0_0, vert2_0 - vert0_0)));
        *normal_0 = _S1393 * make_float3 (float(- (F32_sign((dot_0(_S1393, mean_c_8))))));
        *out_hardness_0 = hardness_0;
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_10, float4  quat_13, float3  scale_12, float2  hardness_1, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_14, float fy_14, float cx_14, float cy_14, float4  radial_coeffs_13, float2  tangential_coeffs_13, float2  thin_prism_coeffs_13, uint image_width_10, uint image_height_10, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float * depth_7, float3  * normal_1, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float2  * out_hardness_1)
{
    for(;;)
    {
        float3  mean_c_9 = mul_0(R_12, mean_10) + t_11;
        float _S1394 = mean_c_9.z;
        bool _S1395;
        if(_S1394 < near_plane_7)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = _S1394 > far_plane_7;
        }
        if(_S1395)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S1396 = scale_12.x;
        float sx_1 = (F32_exp((_S1396)));
        float _S1397 = scale_12.y;
        float sy_1 = (F32_exp((_S1397)));
        float sz_1 = scale_12.z - 0.5f * (_S1396 + _S1397);
        float x_33 = quat_13.y;
        float inv_norm_10 = (F32_rsqrt((x_33 * x_33 + quat_13.z * quat_13.z + quat_13.w * quat_13.w + quat_13.x * quat_13.x)));
        float x_34 = quat_13.y * inv_norm_10;
        float y_15 = quat_13.z * inv_norm_10;
        float z_13 = quat_13.w * inv_norm_10;
        float w_13 = quat_13.x * inv_norm_10;
        float x2_13 = x_34 * x_34;
        float y2_13 = y_15 * y_15;
        float z2_13 = z_13 * z_13;
        float xy_13 = x_34 * y_15;
        float xz_13 = x_34 * z_13;
        float yz_13 = y_15 * z_13;
        float wx_13 = w_13 * x_34;
        float wy_13 = w_13 * y_15;
        float wz_13 = w_13 * z_13;
        Matrix<float, 3, 3>  _S1398 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_13), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_13), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
        float3  vert0_1 = mul_0(_S1398, make_float3 (sx_1, 0.0f, 0.0f)) + mean_10;
        float3  vert1_1 = mul_0(_S1398, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_10;
        float3  vert2_1 = mul_0(_S1398, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_10;
        float3  vert0_c_1 = mul_0(R_12, vert0_1) + t_11;
        float3  vert1_c_1 = mul_0(R_12, vert1_1) + t_11;
        float3  vert2_c_1 = mul_0(R_12, vert2_1) + t_11;
        float _S1399 = vert0_c_1.z;
        if(_S1399 < near_plane_7)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = _S1399 > far_plane_7;
        }
        if(_S1395)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = (vert1_c_1.z) < near_plane_7;
        }
        if(_S1395)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = (vert1_c_1.z) > far_plane_7;
        }
        if(_S1395)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = (vert2_c_1.z) < near_plane_7;
        }
        if(_S1395)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = (vert2_c_1.z) > far_plane_7;
        }
        if(_S1395)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_1 = CameraDistortion_x24init_0(radial_coeffs_13, tangential_coeffs_13, thin_prism_coeffs_13);
        float2  _S1400 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_7 = length_0(_S1400);
        float theta_1 = (F32_atan2((r_7), (_S1399)));
        float k_4;
        if(theta_1 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_1 * theta_1 / 3.0f) / _S1399;
        }
        else
        {
            k_4 = theta_1 / r_7;
        }
        float2  _S1401 = _S1400 * make_float2 (k_4);
        float k1_3 = dist_coeffs_1.radial_coeffs_0.x;
        float k2_3 = dist_coeffs_1.radial_coeffs_0.y;
        float k3_3 = dist_coeffs_1.radial_coeffs_0.z;
        float k4_3 = dist_coeffs_1.radial_coeffs_0.w;
        float p1_3 = dist_coeffs_1.tangential_coeffs_0.x;
        float p2_3 = dist_coeffs_1.tangential_coeffs_0.y;
        float sx1_3 = dist_coeffs_1.thin_prism_coeffs_0.x;
        float sy1_3 = dist_coeffs_1.thin_prism_coeffs_0.y;
        float u_9 = _S1401.x;
        float v_9 = _S1401.y;
        float r2_9 = u_9 * u_9 + v_9 * v_9;
        float _S1402 = 2.0f * p1_3;
        float _S1403 = 2.0f * p2_3;
        float2  _S1404 = _S1401 * make_float2 (1.0f + r2_9 * (k1_3 + r2_9 * (k2_3 + r2_9 * (k3_3 + r2_9 * k4_3)))) + make_float2 (_S1402 * u_9 * v_9 + p2_3 * (r2_9 + 2.0f * u_9 * u_9) + sx1_3 * r2_9, _S1403 * u_9 * v_9 + p1_3 * (r2_9 + 2.0f * v_9 * v_9) + sy1_3 * r2_9);
        *uv0_1 = make_float2 (fx_14 * _S1404.x + cx_14, fy_14 * _S1404.y + cy_14);
        float2  _S1405 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_8 = length_0(_S1405);
        float _S1406 = vert1_c_1.z;
        float theta_2 = (F32_atan2((r_8), (_S1406)));
        if(theta_2 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_2 * theta_2 / 3.0f) / _S1406;
        }
        else
        {
            k_4 = theta_2 / r_8;
        }
        float2  _S1407 = _S1405 * make_float2 (k_4);
        float u_10 = _S1407.x;
        float v_10 = _S1407.y;
        float r2_10 = u_10 * u_10 + v_10 * v_10;
        float2  _S1408 = _S1407 * make_float2 (1.0f + r2_10 * (k1_3 + r2_10 * (k2_3 + r2_10 * (k3_3 + r2_10 * k4_3)))) + make_float2 (_S1402 * u_10 * v_10 + p2_3 * (r2_10 + 2.0f * u_10 * u_10) + sx1_3 * r2_10, _S1403 * u_10 * v_10 + p1_3 * (r2_10 + 2.0f * v_10 * v_10) + sy1_3 * r2_10);
        *uv1_1 = make_float2 (fx_14 * _S1408.x + cx_14, fy_14 * _S1408.y + cy_14);
        float2  _S1409 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_9 = length_0(_S1409);
        float _S1410 = vert2_c_1.z;
        float theta_3 = (F32_atan2((r_9), (_S1410)));
        if(theta_3 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_3 * theta_3 / 3.0f) / _S1410;
        }
        else
        {
            k_4 = theta_3 / r_9;
        }
        float2  _S1411 = _S1409 * make_float2 (k_4);
        float u_11 = _S1411.x;
        float v_11 = _S1411.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float2  _S1412 = _S1411 * make_float2 (1.0f + r2_11 * (k1_3 + r2_11 * (k2_3 + r2_11 * (k3_3 + r2_11 * k4_3)))) + make_float2 (_S1402 * u_11 * v_11 + p2_3 * (r2_11 + 2.0f * u_11 * u_11) + sx1_3 * r2_11, _S1403 * u_11 * v_11 + p1_3 * (r2_11 + 2.0f * v_11 * v_11) + sy1_3 * r2_11);
        float _S1413 = fx_14 * _S1412.x + cx_14;
        float _S1414 = fy_14 * _S1412.y + cy_14;
        float2  _S1415 = make_float2 (_S1413, _S1414);
        *uv2_1 = _S1415;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S1415 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S1415)));
        float xmax_4 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S1413))) + offset_1;
        float xmin_4 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S1413))) - offset_1;
        float ymax_4 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S1414))) + offset_1;
        float ymin_4 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S1414))) - offset_1;
        if(xmax_4 <= 0.0f)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = xmin_4 >= float(image_width_10);
        }
        if(_S1395)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = ymax_4 <= 0.0f;
        }
        if(_S1395)
        {
            _S1395 = true;
        }
        else
        {
            _S1395 = ymin_4 >= float(image_height_10);
        }
        if(_S1395)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_7 = make_int4 (int((F32_floor((xmin_4)))), int((F32_floor((ymin_4)))), int((F32_ceil((xmax_4)))), int((F32_ceil((ymax_4)))));
        *depth_7 = _S1394;
        float3  _S1416 = normalize_0(mul_0(R_12, cross_0(vert1_1 - vert0_1, vert2_1 - vert0_1)));
        *normal_1 = _S1416 * make_float3 (float(- (F32_sign((dot_0(_S1416, mean_c_9))))));
        *out_hardness_1 = hardness_1;
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_11, float4  quat_14, float3  scale_13, float2  hardness_2, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_15, float fy_15, float cx_15, float cy_15, float4  radial_coeffs_14, float2  tangential_coeffs_14, float2  thin_prism_coeffs_14, uint image_width_11, uint image_height_11, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float * depth_8, float3  * normal_2, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float2  * out_hardness_2)
{
    float3  mean_c_10 = mul_0(R_13, mean_11) + t_12;
    float _S1417 = scale_13.x;
    float sx_2 = (F32_exp((_S1417)));
    float _S1418 = scale_13.y;
    float sy_2 = (F32_exp((_S1418)));
    float sz_2 = scale_13.z - 0.5f * (_S1417 + _S1418);
    float x_35 = quat_14.y;
    float inv_norm_11 = (F32_rsqrt((x_35 * x_35 + quat_14.z * quat_14.z + quat_14.w * quat_14.w + quat_14.x * quat_14.x)));
    float x_36 = quat_14.y * inv_norm_11;
    float y_16 = quat_14.z * inv_norm_11;
    float z_14 = quat_14.w * inv_norm_11;
    float w_14 = quat_14.x * inv_norm_11;
    float x2_14 = x_36 * x_36;
    float y2_14 = y_16 * y_16;
    float z2_14 = z_14 * z_14;
    float xy_14 = x_36 * y_16;
    float xz_14 = x_36 * z_14;
    float yz_14 = y_16 * z_14;
    float wx_14 = w_14 * x_36;
    float wy_14 = w_14 * y_16;
    float wz_14 = w_14 * z_14;
    Matrix<float, 3, 3>  _S1419 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_14), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_14), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    float3  vert0_2 = mul_0(_S1419, make_float3 (sx_2, 0.0f, 0.0f)) + mean_11;
    float3  vert1_2 = mul_0(_S1419, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_11;
    float3  vert2_2 = mul_0(_S1419, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_11;
    float3  vert0_c_2 = mul_0(R_13, vert0_2) + t_12;
    float3  vert1_c_2 = mul_0(R_13, vert1_2) + t_12;
    float3  vert2_c_2 = mul_0(R_13, vert2_2) + t_12;
    *uv0_2 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    *uv1_2 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    *uv2_2 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float2  _S1420 = make_float2 (fx_15, fy_15);
    float2  _S1421 = make_float2 (cx_15, cy_15);
    *uv0_2 = _S1420 * *uv0_2 + _S1421;
    *uv1_2 = _S1420 * *uv1_2 + _S1421;
    float2  _S1422 = _S1420 * *uv2_2 + _S1421;
    *uv2_2 = _S1422;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S1422 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S1422)));
    float _S1423 = _S1422.x;
    float _S1424 = _S1422.y;
    *aabb_xyxy_8 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S1423))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S1424))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S1423))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S1424))) + offset_2)))));
    *depth_8 = mean_c_10.z;
    float3  _S1425 = normalize_0(mul_0(R_13, cross_0(vert1_2 - vert0_2, vert2_2 - vert0_2)));
    *normal_2 = _S1425 * make_float3 (float(- (F32_sign((dot_0(_S1425, mean_c_10))))));
    *out_hardness_2 = hardness_2;
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_12, float4  quat_15, float3  scale_14, float2  hardness_3, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_16, float fy_16, float cx_16, float cy_16, float4  radial_coeffs_15, float2  tangential_coeffs_15, float2  thin_prism_coeffs_15, uint image_width_12, uint image_height_12, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float * depth_9, float3  * normal_3, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float2  * out_hardness_3)
{
    float3  mean_c_11 = mul_0(R_14, mean_12) + t_13;
    float _S1426 = scale_14.x;
    float sx_3 = (F32_exp((_S1426)));
    float _S1427 = scale_14.y;
    float sy_3 = (F32_exp((_S1427)));
    float sz_3 = scale_14.z - 0.5f * (_S1426 + _S1427);
    float x_37 = quat_15.y;
    float inv_norm_12 = (F32_rsqrt((x_37 * x_37 + quat_15.z * quat_15.z + quat_15.w * quat_15.w + quat_15.x * quat_15.x)));
    float x_38 = quat_15.y * inv_norm_12;
    float y_17 = quat_15.z * inv_norm_12;
    float z_15 = quat_15.w * inv_norm_12;
    float w_15 = quat_15.x * inv_norm_12;
    float x2_15 = x_38 * x_38;
    float y2_15 = y_17 * y_17;
    float z2_15 = z_15 * z_15;
    float xy_15 = x_38 * y_17;
    float xz_15 = x_38 * z_15;
    float yz_15 = y_17 * z_15;
    float wx_15 = w_15 * x_38;
    float wy_15 = w_15 * y_17;
    float wz_15 = w_15 * z_15;
    Matrix<float, 3, 3>  _S1428 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_15), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_15), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    float3  vert0_3 = mul_0(_S1428, make_float3 (sx_3, 0.0f, 0.0f)) + mean_12;
    float3  vert1_3 = mul_0(_S1428, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_12;
    float3  vert2_3 = mul_0(_S1428, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_12;
    float3  vert0_c_3 = mul_0(R_14, vert0_3) + t_13;
    float3  vert1_c_3 = mul_0(R_14, vert1_3) + t_13;
    float3  vert2_c_3 = mul_0(R_14, vert2_3) + t_13;
    CameraDistortion_0 dist_coeffs_2 = CameraDistortion_x24init_0(radial_coeffs_15, tangential_coeffs_15, thin_prism_coeffs_15);
    float2  _S1429 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_10 = length_0(_S1429);
    float _S1430 = vert0_c_3.z;
    float theta_4 = (F32_atan2((r_10), (_S1430)));
    float k_5;
    if(theta_4 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_4 * theta_4 / 3.0f) / _S1430;
    }
    else
    {
        k_5 = theta_4 / r_10;
    }
    float2  _S1431 = _S1429 * make_float2 (k_5);
    float k1_4 = dist_coeffs_2.radial_coeffs_0.x;
    float k2_4 = dist_coeffs_2.radial_coeffs_0.y;
    float k3_4 = dist_coeffs_2.radial_coeffs_0.z;
    float k4_4 = dist_coeffs_2.radial_coeffs_0.w;
    float p1_4 = dist_coeffs_2.tangential_coeffs_0.x;
    float p2_4 = dist_coeffs_2.tangential_coeffs_0.y;
    float sx1_4 = dist_coeffs_2.thin_prism_coeffs_0.x;
    float sy1_4 = dist_coeffs_2.thin_prism_coeffs_0.y;
    float u_12 = _S1431.x;
    float v_12 = _S1431.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S1432 = 2.0f * p1_4;
    float _S1433 = 2.0f * p2_4;
    float2  _S1434 = _S1431 * make_float2 (1.0f + r2_12 * (k1_4 + r2_12 * (k2_4 + r2_12 * (k3_4 + r2_12 * k4_4)))) + make_float2 (_S1432 * u_12 * v_12 + p2_4 * (r2_12 + 2.0f * u_12 * u_12) + sx1_4 * r2_12, _S1433 * u_12 * v_12 + p1_4 * (r2_12 + 2.0f * v_12 * v_12) + sy1_4 * r2_12);
    *uv0_3 = make_float2 (fx_16 * _S1434.x + cx_16, fy_16 * _S1434.y + cy_16);
    float2  _S1435 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_11 = length_0(_S1435);
    float _S1436 = vert1_c_3.z;
    float theta_5 = (F32_atan2((r_11), (_S1436)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_5 * theta_5 / 3.0f) / _S1436;
    }
    else
    {
        k_5 = theta_5 / r_11;
    }
    float2  _S1437 = _S1435 * make_float2 (k_5);
    float u_13 = _S1437.x;
    float v_13 = _S1437.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S1438 = _S1437 * make_float2 (1.0f + r2_13 * (k1_4 + r2_13 * (k2_4 + r2_13 * (k3_4 + r2_13 * k4_4)))) + make_float2 (_S1432 * u_13 * v_13 + p2_4 * (r2_13 + 2.0f * u_13 * u_13) + sx1_4 * r2_13, _S1433 * u_13 * v_13 + p1_4 * (r2_13 + 2.0f * v_13 * v_13) + sy1_4 * r2_13);
    *uv1_3 = make_float2 (fx_16 * _S1438.x + cx_16, fy_16 * _S1438.y + cy_16);
    float2  _S1439 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_12 = length_0(_S1439);
    float _S1440 = vert2_c_3.z;
    float theta_6 = (F32_atan2((r_12), (_S1440)));
    if(theta_6 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_6 * theta_6 / 3.0f) / _S1440;
    }
    else
    {
        k_5 = theta_6 / r_12;
    }
    float2  _S1441 = _S1439 * make_float2 (k_5);
    float u_14 = _S1441.x;
    float v_14 = _S1441.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S1442 = _S1441 * make_float2 (1.0f + r2_14 * (k1_4 + r2_14 * (k2_4 + r2_14 * (k3_4 + r2_14 * k4_4)))) + make_float2 (_S1432 * u_14 * v_14 + p2_4 * (r2_14 + 2.0f * u_14 * u_14) + sx1_4 * r2_14, _S1433 * u_14 * v_14 + p1_4 * (r2_14 + 2.0f * v_14 * v_14) + sy1_4 * r2_14);
    float _S1443 = fx_16 * _S1442.x + cx_16;
    float _S1444 = fy_16 * _S1442.y + cy_16;
    float2  _S1445 = make_float2 (_S1443, _S1444);
    *uv2_3 = _S1445;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S1445 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S1445)));
    *aabb_xyxy_9 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S1443))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S1444))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S1443))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S1444))) + offset_3)))));
    *depth_9 = mean_c_11.z;
    float3  _S1446 = normalize_0(mul_0(R_14, cross_0(vert1_3 - vert0_3, vert2_3 - vert0_3)));
    *normal_3 = _S1446 * make_float3 (float(- (F32_sign((dot_0(_S1446, mean_c_11))))));
    *out_hardness_3 = hardness_3;
    return;
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S1447, float3  _S1448)
{
    return cross_0(_S1447, _S1448);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S1449, float3  _S1450)
{
    return dot_0(_S1449, _S1450);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1451, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1452, float _S1453)
{
    _d_dot_0(_S1451, _S1452, _S1453);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_13, float _s_dOut_5)
{
    float _S1454 = (*dpx_13).primal_0.x;
    float _S1455 = (*dpx_13).primal_0.y;
    float _S1456 = (*dpx_13).primal_0.z;
    DiffPair_float_0 _S1457;
    (&_S1457)->primal_0 = _S1454 * _S1454 + _S1455 * _S1455 + _S1456 * _S1456;
    (&_S1457)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1457, _s_dOut_5);
    float _S1458 = (*dpx_13).primal_0.z * _S1457.differential_0;
    float _S1459 = _S1458 + _S1458;
    float _S1460 = (*dpx_13).primal_0.y * _S1457.differential_0;
    float _S1461 = _S1460 + _S1460;
    float _S1462 = (*dpx_13).primal_0.x * _S1457.differential_0;
    float _S1463 = _S1462 + _S1462;
    float3  _S1464 = make_float3 (0.0f);
    *&((&_S1464)->z) = _S1459;
    *&((&_S1464)->y) = _S1461;
    *&((&_S1464)->x) = _S1463;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S1464;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1465, float _S1466)
{
    s_bwd_prop_length_impl_1(_S1465, _S1466);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_14, float3  _s_dOut_6)
{
    float _S1467 = length_1((*dpx_14).primal_0);
    float3  _S1468 = (*dpx_14).primal_0 * _s_dOut_6;
    float3  _S1469 = make_float3 (1.0f / _S1467) * _s_dOut_6;
    float _S1470 = - ((_S1468.x + _S1468.y + _S1468.z) / (_S1467 * _S1467));
    float3  _S1471 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1472;
    (&_S1472)->primal_0 = (*dpx_14).primal_0;
    (&_S1472)->differential_0 = _S1471;
    s_bwd_length_impl_1(&_S1472, _S1470);
    float3  _S1473 = _S1469 + _S1472.differential_0;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S1473;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1474, float3  _S1475)
{
    s_bwd_prop_normalize_impl_0(_S1474, _S1475);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1476, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1477, float3  _S1478)
{
    _d_cross_0(_S1476, _S1477, _S1478);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1479, float _S1480)
{
    _d_exp2_0(_S1479, _S1480);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S1481, float _S1482)
{
    _d_abs_0(_S1481, _S1482);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_13, float4  quat_16, float3  scale_15, float2  hardness_4, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_17, float fy_17, float cx_17, float cy_17, float4  radial_coeffs_16, float2  tangential_coeffs_16, float2  thin_prism_coeffs_16, uint image_width_13, uint image_height_13, float v_depth_3, float3  v_normal_0, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float2  v_out_hardness_0, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float2  * v_hardness_0, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_0(R_15, mean_13) + t_14;
    float _S1483 = scale_15.x;
    float _S1484 = s_primal_ctx_exp_1(_S1483);
    float _S1485 = scale_15.y;
    float _S1486 = s_primal_ctx_exp_1(_S1485);
    float sz_4 = scale_15.z - 0.5f * (_S1483 + _S1485);
    float _S1487 = quat_16.y;
    float _S1488 = _S1487 * _S1487 + quat_16.z * quat_16.z + quat_16.w * quat_16.w + quat_16.x * quat_16.x;
    float _S1489 = s_primal_ctx_rsqrt_0(_S1488);
    float x_39 = quat_16.y * _S1489;
    float y_18 = quat_16.z * _S1489;
    float z_16 = quat_16.w * _S1489;
    float w_16 = quat_16.x * _S1489;
    float x2_16 = x_39 * x_39;
    float y2_16 = y_18 * y_18;
    float z2_16 = z_16 * z_16;
    float xy_16 = x_39 * y_18;
    float xz_16 = x_39 * z_16;
    float yz_16 = y_18 * z_16;
    float wx_16 = w_16 * x_39;
    float wy_16 = w_16 * y_18;
    float wz_16 = w_16 * z_16;
    Matrix<float, 3, 3>  _S1490 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_16), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_16), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    float3  _S1491 = make_float3 (_S1484, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S1490, _S1491) + mean_13;
    float _S1492 = -0.5f + sz_4;
    float3  _S1493 = make_float3 (_S1484 * _S1492, _S1486, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S1490, _S1493) + mean_13;
    float _S1494 = -0.5f - sz_4;
    float3  _S1495 = make_float3 (_S1484 * _S1494, - _S1486, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S1490, _S1495) + mean_13;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_15, vert0_4) + t_14;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_15, vert1_4) + t_14;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_15, vert2_4) + t_14;
    float2  _S1496 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S1497 = vert0_c_4.z;
    float2  _S1498 = make_float2 (_S1497);
    float2  _S1499 = make_float2 (_S1497 * _S1497);
    float2  _S1500 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S1501 = vert1_c_4.z;
    float2  _S1502 = make_float2 (_S1501);
    float2  _S1503 = make_float2 (_S1501 * _S1501);
    float2  _S1504 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S1505 = vert2_c_4.z;
    float2  _S1506 = make_float2 (_S1505);
    float2  _S1507 = make_float2 (_S1505 * _S1505);
    float2  _S1508 = make_float2 (fx_17, fy_17);
    float2  _S1509 = make_float2 (cx_17, cy_17);
    float2  _S1510 = _S1508 * (_S1496 / make_float2 (_S1497)) + _S1509;
    float2  _S1511 = _S1508 * (_S1500 / make_float2 (_S1501)) + _S1509;
    float2  _S1512 = _S1508 * (_S1504 / make_float2 (_S1505)) + _S1509;
    float2  e0_4 = _S1511 - _S1510;
    float2  e1_4 = _S1512 - _S1511;
    float2  e2_0 = _S1510 - _S1512;
    float _S1513 = e0_4.x;
    float _S1514 = e1_4.y;
    float _S1515 = e0_4.y;
    float _S1516 = e1_4.x;
    float _S1517 = _S1513 * _S1514 - _S1515 * _S1516;
    float _S1518 = 1.0f - hardness_4.y;
    float _S1519 = -1.0f / _S1518;
    float _S1520 = _S1518 * _S1518;
    float _S1521 = _S1510.x;
    float _S1522 = _S1511.x;
    float _S1523 = s_primal_ctx_max_0(_S1521, _S1522);
    float _S1524 = _S1512.x;
    float _S1525 = s_primal_ctx_min_0(_S1521, _S1522);
    float _S1526 = _S1510.y;
    float _S1527 = _S1511.y;
    float _S1528 = s_primal_ctx_max_0(_S1526, _S1527);
    float _S1529 = _S1512.y;
    float _S1530 = s_primal_ctx_min_0(_S1526, _S1527);
    float3  _S1531 = vert1_4 - vert0_4;
    float3  _S1532 = vert2_4 - vert0_4;
    float3  _S1533 = s_primal_ctx_cross_0(_S1531, _S1532);
    float3  _S1534 = s_primal_ctx_mul_0(R_15, _S1533);
    float3  _S1535 = normalize_0(_S1534);
    float3  _S1536 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S1535, mean_c_12)))))) * v_normal_0;
    float3  _S1537 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1538;
    (&_S1538)->primal_0 = _S1535;
    (&_S1538)->differential_0 = _S1537;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1539;
    (&_S1539)->primal_0 = mean_c_12;
    (&_S1539)->differential_0 = _S1537;
    s_bwd_prop_dot_0(&_S1538, &_S1539, 0.0f);
    float3  _S1540 = _S1536 + _S1538.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1541;
    (&_S1541)->primal_0 = _S1534;
    (&_S1541)->differential_0 = _S1537;
    s_bwd_normalize_impl_0(&_S1541, _S1540);
    Matrix<float, 3, 3>  _S1542 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1543;
    (&_S1543)->primal_0 = R_15;
    (&_S1543)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1544;
    (&_S1544)->primal_0 = _S1533;
    (&_S1544)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1543, &_S1544, _S1541.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1545;
    (&_S1545)->primal_0 = _S1531;
    (&_S1545)->differential_0 = _S1537;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1546;
    (&_S1546)->primal_0 = _S1532;
    (&_S1546)->differential_0 = _S1537;
    s_bwd_prop_cross_0(&_S1545, &_S1546, _S1544.differential_0);
    float3  _S1547 = - _S1546.differential_0;
    float3  _S1548 = - _S1545.differential_0;
    DiffPair_float_0 _S1549;
    (&_S1549)->primal_0 = _S1530;
    (&_S1549)->differential_0 = 0.0f;
    DiffPair_float_0 _S1550;
    (&_S1550)->primal_0 = _S1529;
    (&_S1550)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1549, &_S1550, 0.0f);
    DiffPair_float_0 _S1551;
    (&_S1551)->primal_0 = _S1526;
    (&_S1551)->differential_0 = 0.0f;
    DiffPair_float_0 _S1552;
    (&_S1552)->primal_0 = _S1527;
    (&_S1552)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1551, &_S1552, _S1549.differential_0);
    DiffPair_float_0 _S1553;
    (&_S1553)->primal_0 = _S1528;
    (&_S1553)->differential_0 = 0.0f;
    DiffPair_float_0 _S1554;
    (&_S1554)->primal_0 = _S1529;
    (&_S1554)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1553, &_S1554, 0.0f);
    float _S1555 = _S1550.differential_0 + _S1554.differential_0;
    DiffPair_float_0 _S1556;
    (&_S1556)->primal_0 = _S1526;
    (&_S1556)->differential_0 = 0.0f;
    DiffPair_float_0 _S1557;
    (&_S1557)->primal_0 = _S1527;
    (&_S1557)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1556, &_S1557, _S1553.differential_0);
    float _S1558 = _S1552.differential_0 + _S1557.differential_0;
    float _S1559 = _S1551.differential_0 + _S1556.differential_0;
    DiffPair_float_0 _S1560;
    (&_S1560)->primal_0 = _S1525;
    (&_S1560)->differential_0 = 0.0f;
    DiffPair_float_0 _S1561;
    (&_S1561)->primal_0 = _S1524;
    (&_S1561)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1560, &_S1561, 0.0f);
    DiffPair_float_0 _S1562;
    (&_S1562)->primal_0 = _S1521;
    (&_S1562)->differential_0 = 0.0f;
    DiffPair_float_0 _S1563;
    (&_S1563)->primal_0 = _S1522;
    (&_S1563)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1562, &_S1563, _S1560.differential_0);
    DiffPair_float_0 _S1564;
    (&_S1564)->primal_0 = _S1523;
    (&_S1564)->differential_0 = 0.0f;
    DiffPair_float_0 _S1565;
    (&_S1565)->primal_0 = _S1524;
    (&_S1565)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1564, &_S1565, 0.0f);
    float _S1566 = _S1561.differential_0 + _S1565.differential_0;
    DiffPair_float_0 _S1567;
    (&_S1567)->primal_0 = _S1521;
    (&_S1567)->differential_0 = 0.0f;
    DiffPair_float_0 _S1568;
    (&_S1568)->primal_0 = _S1522;
    (&_S1568)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1567, &_S1568, _S1564.differential_0);
    float _S1569 = _S1563.differential_0 + _S1568.differential_0;
    float _S1570 = _S1562.differential_0 + _S1567.differential_0;
    DiffPair_float_0 _S1571;
    (&_S1571)->primal_0 = _S1519;
    (&_S1571)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1571, 0.0f);
    float _S1572 = - (-1.0f * - (_S1571.differential_0 / _S1520));
    float2  _S1573 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1574;
    (&_S1574)->primal_0 = e2_0;
    (&_S1574)->differential_0 = _S1573;
    s_bwd_length_impl_0(&_S1574, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1575;
    (&_S1575)->primal_0 = e1_4;
    (&_S1575)->differential_0 = _S1573;
    s_bwd_length_impl_0(&_S1575, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1576;
    (&_S1576)->primal_0 = e0_4;
    (&_S1576)->differential_0 = _S1573;
    s_bwd_length_impl_0(&_S1576, -0.0f);
    DiffPair_float_0 _S1577;
    (&_S1577)->primal_0 = _S1517;
    (&_S1577)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1577, 0.0f);
    float _S1578 = - _S1577.differential_0;
    float2  _S1579 = _S1575.differential_0 + make_float2 (_S1515 * _S1578, _S1513 * _S1577.differential_0);
    float2  _S1580 = _S1576.differential_0 + make_float2 (_S1514 * _S1577.differential_0, _S1516 * _S1578);
    float2  _S1581 = _S1508 * (v_uv2_0 + - _S1574.differential_0 + _S1579 + make_float2 (_S1566, _S1555)) / _S1507;
    float2  _S1582 = _S1504 * - _S1581;
    float2  _S1583 = _S1506 * _S1581;
    float2  _S1584 = _S1508 * (v_uv1_0 + - _S1579 + _S1580 + make_float2 (_S1569, _S1558)) / _S1503;
    float2  _S1585 = _S1500 * - _S1584;
    float2  _S1586 = _S1502 * _S1584;
    float _S1587 = _S1585.x + _S1585.y;
    float2  _S1588 = _S1508 * (v_uv0_0 + _S1574.differential_0 + - _S1580 + make_float2 (_S1570, _S1559)) / _S1499;
    float2  _S1589 = _S1496 * - _S1588;
    float2  _S1590 = _S1498 * _S1588;
    float _S1591 = _S1589.x + _S1589.y;
    float3  s_diff_vert2_c_T_0 = make_float3 (_S1583.x, _S1583.y, _S1582.x + _S1582.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1592;
    (&_S1592)->primal_0 = R_15;
    (&_S1592)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1593;
    (&_S1593)->primal_0 = vert2_4;
    (&_S1593)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1592, &_S1593, s_diff_vert2_c_T_0);
    float3  s_diff_vert1_c_T_0 = make_float3 (_S1586.x, _S1586.y, _S1587);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1594;
    (&_S1594)->primal_0 = R_15;
    (&_S1594)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1595;
    (&_S1595)->primal_0 = vert1_4;
    (&_S1595)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1594, &_S1595, s_diff_vert1_c_T_0);
    float3  s_diff_vert0_c_T_0 = make_float3 (_S1590.x, _S1590.y, _S1591);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1596;
    (&_S1596)->primal_0 = R_15;
    (&_S1596)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1597;
    (&_S1597)->primal_0 = vert0_4;
    (&_S1597)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1596, &_S1597, s_diff_vert0_c_T_0);
    float3  _S1598 = _S1546.differential_0 + _S1593.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1599;
    (&_S1599)->primal_0 = _S1490;
    (&_S1599)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1600;
    (&_S1600)->primal_0 = _S1495;
    (&_S1600)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1599, &_S1600, _S1598);
    float _S1601 = - _S1600.differential_0.y;
    float _S1602 = _S1494 * _S1600.differential_0.x;
    float _S1603 = - (_S1484 * _S1600.differential_0.x);
    float3  _S1604 = _S1545.differential_0 + _S1595.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1605;
    (&_S1605)->primal_0 = _S1490;
    (&_S1605)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1606;
    (&_S1606)->primal_0 = _S1493;
    (&_S1606)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1605, &_S1606, _S1604);
    float _S1607 = _S1484 * _S1606.differential_0.x;
    float _S1608 = _S1492 * _S1606.differential_0.x;
    float3  _S1609 = _S1547 + _S1548 + _S1597.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1610;
    (&_S1610)->primal_0 = _S1490;
    (&_S1610)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1611;
    (&_S1611)->primal_0 = _S1491;
    (&_S1611)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1610, &_S1611, _S1609);
    Matrix<float, 3, 3>  _S1612 = transpose_0(_S1599.differential_0 + _S1605.differential_0 + _S1610.differential_0);
    float _S1613 = 2.0f * - _S1612.rows[int(2)].z;
    float _S1614 = 2.0f * _S1612.rows[int(2)].y;
    float _S1615 = 2.0f * _S1612.rows[int(2)].x;
    float _S1616 = 2.0f * _S1612.rows[int(1)].z;
    float _S1617 = 2.0f * - _S1612.rows[int(1)].y;
    float _S1618 = 2.0f * _S1612.rows[int(1)].x;
    float _S1619 = 2.0f * _S1612.rows[int(0)].z;
    float _S1620 = 2.0f * _S1612.rows[int(0)].y;
    float _S1621 = 2.0f * - _S1612.rows[int(0)].x;
    float _S1622 = - _S1618 + _S1620;
    float _S1623 = _S1615 + - _S1619;
    float _S1624 = - _S1614 + _S1616;
    float _S1625 = _S1614 + _S1616;
    float _S1626 = _S1615 + _S1619;
    float _S1627 = _S1618 + _S1620;
    float _S1628 = z_16 * (_S1617 + _S1621);
    float _S1629 = y_18 * (_S1613 + _S1621);
    float _S1630 = x_39 * (_S1613 + _S1617);
    float _S1631 = z_16 * _S1622 + y_18 * _S1623 + x_39 * _S1624;
    float _S1632 = _S1489 * _S1631;
    float _S1633 = w_16 * _S1622 + y_18 * _S1625 + x_39 * _S1626 + _S1628 + _S1628;
    float _S1634 = _S1489 * _S1633;
    float _S1635 = w_16 * _S1623 + z_16 * _S1625 + x_39 * _S1627 + _S1629 + _S1629;
    float _S1636 = _S1489 * _S1635;
    float _S1637 = w_16 * _S1624 + z_16 * _S1626 + y_18 * _S1627 + _S1630 + _S1630;
    float _S1638 = _S1489 * _S1637;
    float _S1639 = quat_16.x * _S1631 + quat_16.w * _S1633 + quat_16.z * _S1635 + quat_16.y * _S1637;
    DiffPair_float_0 _S1640;
    (&_S1640)->primal_0 = _S1488;
    (&_S1640)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1640, _S1639);
    float _S1641 = quat_16.x * _S1640.differential_0;
    float _S1642 = quat_16.w * _S1640.differential_0;
    float _S1643 = quat_16.z * _S1640.differential_0;
    float _S1644 = quat_16.y * _S1640.differential_0;
    float _S1645 = _S1634 + _S1642 + _S1642;
    float _S1646 = _S1636 + _S1643 + _S1643;
    float _S1647 = _S1638 + _S1644 + _S1644;
    float _S1648 = _S1632 + _S1641 + _S1641;
    float _S1649 = _S1603 + _S1607;
    float _S1650 = 0.5f * - _S1649;
    float _S1651 = _S1601 + _S1606.differential_0.y;
    DiffPair_float_0 _S1652;
    (&_S1652)->primal_0 = _S1485;
    (&_S1652)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1652, _S1651);
    float _S1653 = _S1650 + _S1652.differential_0;
    float _S1654 = _S1602 + _S1608 + _S1611.differential_0.x;
    DiffPair_float_0 _S1655;
    (&_S1655)->primal_0 = _S1483;
    (&_S1655)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1655, _S1654);
    float _S1656 = _S1650 + _S1655.differential_0;
    float3  _S1657 = _S1539.differential_0 + make_float3 (0.0f, 0.0f, v_depth_3);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1658;
    (&_S1658)->primal_0 = R_15;
    (&_S1658)->differential_0 = _S1542;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1659;
    (&_S1659)->primal_0 = mean_13;
    (&_S1659)->differential_0 = _S1537;
    s_bwd_prop_mul_3(&_S1658, &_S1659, _S1657);
    float3  _S1660 = s_diff_vert2_c_T_0 + s_diff_vert1_c_T_0 + s_diff_vert0_c_T_0 + _S1657;
    Matrix<float, 3, 3>  _S1661 = _S1543.differential_0 + _S1592.differential_0 + _S1594.differential_0 + _S1596.differential_0 + _S1658.differential_0;
    float2  _S1662 = v_out_hardness_0 + make_float2 (0.0f, _S1572);
    float3  _S1663 = make_float3 (_S1656, _S1653, _S1649);
    float4  _S1664 = make_float4 (0.0f);
    *&((&_S1664)->w) = _S1645;
    *&((&_S1664)->z) = _S1646;
    *&((&_S1664)->y) = _S1647;
    *&((&_S1664)->x) = _S1648;
    *v_mean_3 = _S1598 + _S1604 + _S1609 + _S1659.differential_0;
    *v_quat_3 = _S1664;
    *v_scale_3 = _S1663;
    *v_hardness_0 = _S1662;
    *v_R_3 = _S1661;
    *v_t_3 = _S1660;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_14, float4  quat_17, float3  scale_16, float2  hardness_5, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_18, float fy_18, float cx_18, float cy_18, float4  radial_coeffs_17, float2  tangential_coeffs_17, float2  thin_prism_coeffs_17, uint image_width_14, uint image_height_14, float v_depth_4, float3  v_normal_1, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float2  v_out_hardness_1, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float2  * v_hardness_1, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float3  mean_c_13 = s_primal_ctx_mul_0(R_16, mean_14) + t_15;
    float _S1665 = scale_16.x;
    float _S1666 = s_primal_ctx_exp_1(_S1665);
    float _S1667 = scale_16.y;
    float _S1668 = s_primal_ctx_exp_1(_S1667);
    float sz_5 = scale_16.z - 0.5f * (_S1665 + _S1667);
    float _S1669 = quat_17.y;
    float _S1670 = _S1669 * _S1669 + quat_17.z * quat_17.z + quat_17.w * quat_17.w + quat_17.x * quat_17.x;
    float _S1671 = s_primal_ctx_rsqrt_0(_S1670);
    float x_40 = quat_17.y * _S1671;
    float y_19 = quat_17.z * _S1671;
    float z_17 = quat_17.w * _S1671;
    float w_17 = quat_17.x * _S1671;
    float x2_17 = x_40 * x_40;
    float y2_17 = y_19 * y_19;
    float z2_17 = z_17 * z_17;
    float xy_17 = x_40 * y_19;
    float xz_17 = x_40 * z_17;
    float yz_17 = y_19 * z_17;
    float wx_17 = w_17 * x_40;
    float wy_17 = w_17 * y_19;
    float wz_17 = w_17 * z_17;
    Matrix<float, 3, 3>  _S1672 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_17), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_17), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    float3  _S1673 = make_float3 (_S1666, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S1672, _S1673) + mean_14;
    float _S1674 = -0.5f + sz_5;
    float3  _S1675 = make_float3 (_S1666 * _S1674, _S1668, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S1672, _S1675) + mean_14;
    float _S1676 = -0.5f - sz_5;
    float3  _S1677 = make_float3 (_S1666 * _S1676, - _S1668, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S1672, _S1677) + mean_14;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_16, vert0_5) + t_15;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_16, vert1_5) + t_15;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_16, vert2_5) + t_15;
    CameraDistortion_0 _S1678 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_17, tangential_coeffs_17, thin_prism_coeffs_17);
    float2  _S1679 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S1680 = length_0(_S1679);
    float _S1681 = vert0_c_5.z;
    float _S1682 = s_primal_ctx_atan2_0(_S1680, _S1681);
    bool _S1683 = _S1682 < 0.00100000004749745f;
    float k_6;
    float _S1684;
    float _S1685;
    float _S1686;
    if(_S1683)
    {
        float _S1687 = 1.0f - _S1682 * _S1682 / 3.0f;
        float _S1688 = _S1681 * _S1681;
        k_6 = _S1687 / _S1681;
        _S1684 = 0.0f;
        _S1685 = _S1688;
        _S1686 = _S1687;
    }
    else
    {
        float _S1689 = _S1680 * _S1680;
        k_6 = _S1682 / _S1680;
        _S1684 = _S1689;
        _S1685 = 0.0f;
        _S1686 = 0.0f;
    }
    float2  _S1690 = make_float2 (k_6);
    float2  _S1691 = _S1679 * make_float2 (k_6);
    float k1_5 = _S1678.radial_coeffs_0.x;
    float k2_5 = _S1678.radial_coeffs_0.y;
    float k3_5 = _S1678.radial_coeffs_0.z;
    float k4_5 = _S1678.radial_coeffs_0.w;
    float p1_5 = _S1678.tangential_coeffs_0.x;
    float p2_5 = _S1678.tangential_coeffs_0.y;
    float sx1_5 = _S1678.thin_prism_coeffs_0.x;
    float sy1_5 = _S1678.thin_prism_coeffs_0.y;
    float u_15 = _S1691.x;
    float v_15 = _S1691.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S1692 = k3_5 + r2_15 * k4_5;
    float _S1693 = k2_5 + r2_15 * _S1692;
    float _S1694 = k1_5 + r2_15 * _S1693;
    float radial_1 = 1.0f + r2_15 * _S1694;
    float2  _S1695 = make_float2 (radial_1);
    float _S1696 = 2.0f * p1_5;
    float _S1697 = _S1696 * u_15;
    float _S1698 = 2.0f * u_15;
    float _S1699 = r2_15 + _S1698 * u_15;
    float _S1700 = 2.0f * p2_5;
    float _S1701 = _S1700 * u_15;
    float _S1702 = 2.0f * v_15;
    float _S1703 = r2_15 + _S1702 * v_15;
    float2  _S1704 = _S1691 * make_float2 (radial_1) + make_float2 (_S1697 * v_15 + p2_5 * _S1699 + sx1_5 * r2_15, _S1701 * v_15 + p1_5 * _S1703 + sy1_5 * r2_15);
    float _S1705 = fx_18 * _S1704.x + cx_18;
    float _S1706 = fy_18 * _S1704.y + cy_18;
    float2  _S1707 = make_float2 (_S1705, _S1706);
    float2  _S1708 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S1709 = length_0(_S1708);
    float _S1710 = vert1_c_5.z;
    float _S1711 = s_primal_ctx_atan2_0(_S1709, _S1710);
    bool _S1712 = _S1711 < 0.00100000004749745f;
    float _S1713;
    float _S1714;
    float _S1715;
    if(_S1712)
    {
        float _S1716 = 1.0f - _S1711 * _S1711 / 3.0f;
        float _S1717 = _S1710 * _S1710;
        k_6 = _S1716 / _S1710;
        _S1713 = 0.0f;
        _S1714 = _S1717;
        _S1715 = _S1716;
    }
    else
    {
        float _S1718 = _S1709 * _S1709;
        k_6 = _S1711 / _S1709;
        _S1713 = _S1718;
        _S1714 = 0.0f;
        _S1715 = 0.0f;
    }
    float2  _S1719 = make_float2 (k_6);
    float2  _S1720 = _S1708 * make_float2 (k_6);
    float u_16 = _S1720.x;
    float v_16 = _S1720.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S1721 = k3_5 + r2_16 * k4_5;
    float _S1722 = k2_5 + r2_16 * _S1721;
    float _S1723 = k1_5 + r2_16 * _S1722;
    float radial_2 = 1.0f + r2_16 * _S1723;
    float2  _S1724 = make_float2 (radial_2);
    float _S1725 = _S1696 * u_16;
    float _S1726 = 2.0f * u_16;
    float _S1727 = r2_16 + _S1726 * u_16;
    float _S1728 = _S1700 * u_16;
    float _S1729 = 2.0f * v_16;
    float _S1730 = r2_16 + _S1729 * v_16;
    float2  _S1731 = _S1720 * make_float2 (radial_2) + make_float2 (_S1725 * v_16 + p2_5 * _S1727 + sx1_5 * r2_16, _S1728 * v_16 + p1_5 * _S1730 + sy1_5 * r2_16);
    float _S1732 = fx_18 * _S1731.x + cx_18;
    float _S1733 = fy_18 * _S1731.y + cy_18;
    float2  _S1734 = make_float2 (_S1732, _S1733);
    float2  _S1735 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S1736 = length_0(_S1735);
    float _S1737 = vert2_c_5.z;
    float _S1738 = s_primal_ctx_atan2_0(_S1736, _S1737);
    bool _S1739 = _S1738 < 0.00100000004749745f;
    float _S1740;
    float _S1741;
    float _S1742;
    if(_S1739)
    {
        float _S1743 = 1.0f - _S1738 * _S1738 / 3.0f;
        float _S1744 = _S1737 * _S1737;
        k_6 = _S1743 / _S1737;
        _S1740 = 0.0f;
        _S1741 = _S1744;
        _S1742 = _S1743;
    }
    else
    {
        float _S1745 = _S1736 * _S1736;
        k_6 = _S1738 / _S1736;
        _S1740 = _S1745;
        _S1741 = 0.0f;
        _S1742 = 0.0f;
    }
    float2  _S1746 = make_float2 (k_6);
    float2  _S1747 = _S1735 * make_float2 (k_6);
    float u_17 = _S1747.x;
    float v_17 = _S1747.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S1748 = k3_5 + r2_17 * k4_5;
    float _S1749 = k2_5 + r2_17 * _S1748;
    float _S1750 = k1_5 + r2_17 * _S1749;
    float radial_3 = 1.0f + r2_17 * _S1750;
    float2  _S1751 = make_float2 (radial_3);
    float _S1752 = _S1696 * u_17;
    float _S1753 = 2.0f * u_17;
    float _S1754 = r2_17 + _S1753 * u_17;
    float _S1755 = _S1700 * u_17;
    float _S1756 = 2.0f * v_17;
    float _S1757 = r2_17 + _S1756 * v_17;
    float2  _S1758 = _S1747 * make_float2 (radial_3) + make_float2 (_S1752 * v_17 + p2_5 * _S1754 + sx1_5 * r2_17, _S1755 * v_17 + p1_5 * _S1757 + sy1_5 * r2_17);
    float _S1759 = fx_18 * _S1758.x + cx_18;
    float _S1760 = fy_18 * _S1758.y + cy_18;
    float2  _S1761 = make_float2 (_S1759, _S1760);
    float2  e0_5 = _S1734 - _S1707;
    float2  e1_5 = _S1761 - _S1734;
    float2  e2_1 = _S1707 - _S1761;
    float _S1762 = e0_5.x;
    float _S1763 = e1_5.y;
    float _S1764 = e0_5.y;
    float _S1765 = e1_5.x;
    float _S1766 = _S1762 * _S1763 - _S1764 * _S1765;
    float _S1767 = 1.0f - hardness_5.y;
    float _S1768 = -1.0f / _S1767;
    float _S1769 = _S1767 * _S1767;
    float _S1770 = s_primal_ctx_max_0(_S1705, _S1732);
    float _S1771 = s_primal_ctx_min_0(_S1705, _S1732);
    float _S1772 = s_primal_ctx_max_0(_S1706, _S1733);
    float _S1773 = s_primal_ctx_min_0(_S1706, _S1733);
    float3  _S1774 = vert1_5 - vert0_5;
    float3  _S1775 = vert2_5 - vert0_5;
    float3  _S1776 = s_primal_ctx_cross_0(_S1774, _S1775);
    float3  _S1777 = s_primal_ctx_mul_0(R_16, _S1776);
    float3  _S1778 = normalize_0(_S1777);
    float3  _S1779 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S1778, mean_c_13)))))) * v_normal_1;
    float3  _S1780 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1781;
    (&_S1781)->primal_0 = _S1778;
    (&_S1781)->differential_0 = _S1780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1782;
    (&_S1782)->primal_0 = mean_c_13;
    (&_S1782)->differential_0 = _S1780;
    s_bwd_prop_dot_0(&_S1781, &_S1782, 0.0f);
    float3  _S1783 = _S1779 + _S1781.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1784;
    (&_S1784)->primal_0 = _S1777;
    (&_S1784)->differential_0 = _S1780;
    s_bwd_normalize_impl_0(&_S1784, _S1783);
    Matrix<float, 3, 3>  _S1785 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1786;
    (&_S1786)->primal_0 = R_16;
    (&_S1786)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1787;
    (&_S1787)->primal_0 = _S1776;
    (&_S1787)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S1786, &_S1787, _S1784.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1788 = _S1786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1789;
    (&_S1789)->primal_0 = _S1774;
    (&_S1789)->differential_0 = _S1780;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1790;
    (&_S1790)->primal_0 = _S1775;
    (&_S1790)->differential_0 = _S1780;
    s_bwd_prop_cross_0(&_S1789, &_S1790, _S1787.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1791 = _S1789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1792 = _S1790;
    float3  _S1793 = - _S1790.differential_0;
    float3  _S1794 = - _S1789.differential_0;
    DiffPair_float_0 _S1795;
    (&_S1795)->primal_0 = _S1773;
    (&_S1795)->differential_0 = 0.0f;
    DiffPair_float_0 _S1796;
    (&_S1796)->primal_0 = _S1760;
    (&_S1796)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1795, &_S1796, 0.0f);
    DiffPair_float_0 _S1797;
    (&_S1797)->primal_0 = _S1706;
    (&_S1797)->differential_0 = 0.0f;
    DiffPair_float_0 _S1798;
    (&_S1798)->primal_0 = _S1733;
    (&_S1798)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1797, &_S1798, _S1795.differential_0);
    DiffPair_float_0 _S1799;
    (&_S1799)->primal_0 = _S1772;
    (&_S1799)->differential_0 = 0.0f;
    DiffPair_float_0 _S1800;
    (&_S1800)->primal_0 = _S1760;
    (&_S1800)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1799, &_S1800, 0.0f);
    DiffPair_float_0 _S1801;
    (&_S1801)->primal_0 = _S1706;
    (&_S1801)->differential_0 = 0.0f;
    DiffPair_float_0 _S1802;
    (&_S1802)->primal_0 = _S1733;
    (&_S1802)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1801, &_S1802, _S1799.differential_0);
    DiffPair_float_0 _S1803;
    (&_S1803)->primal_0 = _S1771;
    (&_S1803)->differential_0 = 0.0f;
    DiffPair_float_0 _S1804;
    (&_S1804)->primal_0 = _S1759;
    (&_S1804)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1803, &_S1804, 0.0f);
    DiffPair_float_0 _S1805;
    (&_S1805)->primal_0 = _S1705;
    (&_S1805)->differential_0 = 0.0f;
    DiffPair_float_0 _S1806;
    (&_S1806)->primal_0 = _S1732;
    (&_S1806)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1805, &_S1806, _S1803.differential_0);
    DiffPair_float_0 _S1807;
    (&_S1807)->primal_0 = _S1770;
    (&_S1807)->differential_0 = 0.0f;
    DiffPair_float_0 _S1808;
    (&_S1808)->primal_0 = _S1759;
    (&_S1808)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1807, &_S1808, 0.0f);
    DiffPair_float_0 _S1809;
    (&_S1809)->primal_0 = _S1705;
    (&_S1809)->differential_0 = 0.0f;
    DiffPair_float_0 _S1810;
    (&_S1810)->primal_0 = _S1732;
    (&_S1810)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S1809, &_S1810, _S1807.differential_0);
    DiffPair_float_0 _S1811;
    (&_S1811)->primal_0 = _S1768;
    (&_S1811)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S1811, 0.0f);
    float _S1812 = - (-1.0f * - (_S1811.differential_0 / _S1769));
    float2  _S1813 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1814;
    (&_S1814)->primal_0 = e2_1;
    (&_S1814)->differential_0 = _S1813;
    s_bwd_length_impl_0(&_S1814, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1815;
    (&_S1815)->primal_0 = e1_5;
    (&_S1815)->differential_0 = _S1813;
    s_bwd_length_impl_0(&_S1815, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1816;
    (&_S1816)->primal_0 = e0_5;
    (&_S1816)->differential_0 = _S1813;
    s_bwd_length_impl_0(&_S1816, -0.0f);
    DiffPair_float_0 _S1817;
    (&_S1817)->primal_0 = _S1766;
    (&_S1817)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S1817, 0.0f);
    float _S1818 = - _S1817.differential_0;
    float2  _S1819 = _S1815.differential_0 + make_float2 (_S1764 * _S1818, _S1762 * _S1817.differential_0);
    float2  _S1820 = _S1816.differential_0 + make_float2 (_S1763 * _S1817.differential_0, _S1765 * _S1818);
    float2  _S1821 = v_uv2_1 + - _S1814.differential_0 + _S1819;
    float _S1822 = fy_18 * (_S1796.differential_0 + _S1800.differential_0 + _S1821.y);
    float _S1823 = fx_18 * (_S1804.differential_0 + _S1808.differential_0 + _S1821.x);
    float2  _S1824 = make_float2 (_S1823, _S1822);
    float2  _S1825 = _S1747 * _S1824;
    float _S1826 = r2_17 * _S1822;
    float _S1827 = p1_5 * _S1822;
    float _S1828 = _S1757 * _S1822;
    float _S1829 = v_17 * _S1822;
    float _S1830 = u_17 * _S1829;
    float _S1831 = r2_17 * _S1823;
    float _S1832 = p2_5 * _S1823;
    float _S1833 = _S1754 * _S1823;
    float _S1834 = v_17 * _S1823;
    float _S1835 = u_17 * _S1834;
    float _S1836 = _S1825.x + _S1825.y;
    float _S1837 = r2_17 * _S1836;
    float _S1838 = r2_17 * _S1837;
    float _S1839 = r2_17 * _S1838;
    float _S1840 = r2_17 * _S1839;
    float _S1841 = sy1_5 * _S1822 + _S1827 + sx1_5 * _S1823 + _S1832 + _S1750 * _S1836 + _S1749 * _S1837 + _S1748 * _S1838 + k4_5 * _S1839;
    float _S1842 = v_17 * _S1841;
    float _S1843 = u_17 * _S1841;
    float2  _S1844 = v_out_hardness_1 + make_float2 (0.0f, _S1812);
    float _S1845 = _S1806.differential_0 + _S1810.differential_0;
    float2  _S1846 = v_uv1_1 + - _S1819 + _S1820;
    float2  _S1847 = v_uv0_1 + _S1814.differential_0 + - _S1820;
    float3  _S1848 = _S1782.differential_0 + make_float3 (0.0f, 0.0f, v_depth_4);
    float3  _S1849 = _S1793 + _S1794;
    float _S1850 = _S1805.differential_0 + _S1809.differential_0;
    float _S1851 = _S1797.differential_0 + _S1801.differential_0;
    float _S1852 = _S1798.differential_0 + _S1802.differential_0;
    float2  _S1853 = _S1751 * _S1824 + make_float2 (_S1700 * _S1829 + _S1753 * _S1832 + 2.0f * (u_17 * _S1832) + _S1696 * _S1834 + _S1843 + _S1843, _S1756 * _S1827 + 2.0f * (v_17 * _S1827) + _S1755 * _S1822 + _S1752 * _S1823 + _S1842 + _S1842);
    float2  _S1854 = _S1735 * _S1853;
    float2  _S1855 = _S1746 * _S1853;
    float _S1856 = _S1854.x + _S1854.y;
    if(_S1739)
    {
        float _S1857 = _S1856 / _S1741;
        float _S1858 = _S1742 * - _S1857;
        float _S1859 = _S1738 * (0.3333333432674408f * - (_S1737 * _S1857));
        k_6 = _S1859 + _S1859;
        _S1740 = _S1858;
        _S1741 = 0.0f;
    }
    else
    {
        float _S1860 = _S1856 / _S1740;
        float _S1861 = _S1738 * - _S1860;
        k_6 = _S1736 * _S1860;
        _S1740 = 0.0f;
        _S1741 = _S1861;
    }
    DiffPair_float_0 _S1862;
    (&_S1862)->primal_0 = _S1736;
    (&_S1862)->differential_0 = 0.0f;
    DiffPair_float_0 _S1863;
    (&_S1863)->primal_0 = _S1737;
    (&_S1863)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1862, &_S1863, k_6);
    float _S1864 = _S1863.differential_0 + _S1740;
    float _S1865 = _S1862.differential_0 + _S1741;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1866;
    (&_S1866)->primal_0 = _S1735;
    (&_S1866)->differential_0 = _S1813;
    s_bwd_length_impl_0(&_S1866, _S1865);
    float2  _S1867 = _S1866.differential_0 + _S1855;
    float _S1868 = fy_18 * (_S1846.y + _S1852);
    float _S1869 = fx_18 * (_S1846.x + _S1845);
    float2  _S1870 = make_float2 (_S1869, _S1868);
    float2  _S1871 = _S1720 * _S1870;
    float _S1872 = p1_5 * _S1868;
    float _S1873 = v_16 * _S1868;
    float _S1874 = p2_5 * _S1869;
    float _S1875 = v_16 * _S1869;
    float _S1876 = _S1871.x + _S1871.y;
    float _S1877 = r2_16 * _S1876;
    float _S1878 = r2_16 * _S1877;
    float _S1879 = r2_16 * _S1878;
    float _S1880 = sy1_5 * _S1868 + _S1872 + sx1_5 * _S1869 + _S1874 + _S1723 * _S1876 + _S1722 * _S1877 + _S1721 * _S1878 + k4_5 * _S1879;
    float _S1881 = v_16 * _S1880;
    float _S1882 = u_16 * _S1880;
    float3  _S1883 = make_float3 (_S1867.x, _S1867.y, _S1864);
    float2  _S1884 = _S1724 * _S1870 + make_float2 (_S1700 * _S1873 + _S1726 * _S1874 + 2.0f * (u_16 * _S1874) + _S1696 * _S1875 + _S1882 + _S1882, _S1729 * _S1872 + 2.0f * (v_16 * _S1872) + _S1728 * _S1868 + _S1725 * _S1869 + _S1881 + _S1881);
    float _S1885 = u_16 * _S1873 + _S1830;
    float _S1886 = u_16 * _S1875 + _S1835;
    float _S1887 = r2_16 * _S1868 + _S1826;
    float _S1888 = r2_16 * _S1879 + _S1840;
    float _S1889 = _S1879 + _S1839;
    float _S1890 = _S1730 * _S1868 + _S1828;
    float _S1891 = _S1878 + _S1838;
    float _S1892 = r2_16 * _S1869 + _S1831;
    float _S1893 = _S1727 * _S1869 + _S1833;
    float _S1894 = _S1877 + _S1837;
    float2  _S1895 = _S1708 * _S1884;
    float2  _S1896 = _S1719 * _S1884;
    float _S1897 = _S1895.x + _S1895.y;
    if(_S1712)
    {
        float _S1898 = _S1897 / _S1714;
        float _S1899 = _S1715 * - _S1898;
        float _S1900 = _S1711 * (0.3333333432674408f * - (_S1710 * _S1898));
        k_6 = _S1900 + _S1900;
        _S1713 = _S1899;
        _S1714 = 0.0f;
    }
    else
    {
        float _S1901 = _S1897 / _S1713;
        float _S1902 = _S1711 * - _S1901;
        k_6 = _S1709 * _S1901;
        _S1713 = 0.0f;
        _S1714 = _S1902;
    }
    DiffPair_float_0 _S1903;
    (&_S1903)->primal_0 = _S1709;
    (&_S1903)->differential_0 = 0.0f;
    DiffPair_float_0 _S1904;
    (&_S1904)->primal_0 = _S1710;
    (&_S1904)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1903, &_S1904, k_6);
    float _S1905 = _S1904.differential_0 + _S1713;
    float _S1906 = _S1903.differential_0 + _S1714;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1907;
    (&_S1907)->primal_0 = _S1708;
    (&_S1907)->differential_0 = _S1813;
    s_bwd_length_impl_0(&_S1907, _S1906);
    float2  _S1908 = _S1907.differential_0 + _S1896;
    float _S1909 = fy_18 * (_S1847.y + _S1851);
    float _S1910 = fx_18 * (_S1847.x + _S1850);
    float2  _S1911 = make_float2 (_S1910, _S1909);
    float2  _S1912 = _S1691 * _S1911;
    float _S1913 = p1_5 * _S1909;
    float _S1914 = v_15 * _S1909;
    float _S1915 = p2_5 * _S1910;
    float _S1916 = v_15 * _S1910;
    float _S1917 = _S1912.x + _S1912.y;
    float _S1918 = r2_15 * _S1917;
    float _S1919 = r2_15 * _S1918;
    float _S1920 = r2_15 * _S1919;
    float _S1921 = sy1_5 * _S1909 + _S1913 + sx1_5 * _S1910 + _S1915 + _S1694 * _S1917 + _S1693 * _S1918 + _S1692 * _S1919 + k4_5 * _S1920;
    float _S1922 = v_15 * _S1921;
    float _S1923 = u_15 * _S1921;
    float2  _S1924 = make_float2 (r2_15 * _S1910 + _S1892, r2_15 * _S1909 + _S1887);
    float2  _S1925 = make_float2 (_S1703 * _S1909 + 2.0f * (u_15 * _S1916 + _S1886) + _S1890, 2.0f * (u_15 * _S1914 + _S1885) + _S1699 * _S1910 + _S1893);
    float4  _S1926 = make_float4 (_S1918 + _S1894, _S1919 + _S1891, _S1920 + _S1889, r2_15 * _S1920 + _S1888);
    float3  _S1927 = make_float3 (_S1908.x, _S1908.y, _S1905);
    float2  _S1928 = _S1695 * _S1911 + make_float2 (_S1700 * _S1914 + _S1698 * _S1915 + 2.0f * (u_15 * _S1915) + _S1696 * _S1916 + _S1923 + _S1923, _S1702 * _S1913 + 2.0f * (v_15 * _S1913) + _S1701 * _S1909 + _S1697 * _S1910 + _S1922 + _S1922);
    CameraDistortion_0 _S1929 = CameraDistortion_x24_syn_dzero_0();
    (&_S1929)->thin_prism_coeffs_0 = _S1924;
    (&_S1929)->tangential_coeffs_0 = _S1925;
    (&_S1929)->radial_coeffs_0 = _S1926;
    float2  _S1930 = _S1679 * _S1928;
    float2  _S1931 = _S1690 * _S1928;
    float _S1932 = _S1930.x + _S1930.y;
    if(_S1683)
    {
        float _S1933 = _S1932 / _S1685;
        float _S1934 = _S1686 * - _S1933;
        float _S1935 = _S1682 * (0.3333333432674408f * - (_S1681 * _S1933));
        k_6 = _S1935 + _S1935;
        _S1684 = _S1934;
        _S1685 = 0.0f;
    }
    else
    {
        float _S1936 = _S1932 / _S1684;
        float _S1937 = _S1682 * - _S1936;
        k_6 = _S1680 * _S1936;
        _S1684 = 0.0f;
        _S1685 = _S1937;
    }
    DiffPair_float_0 _S1938;
    (&_S1938)->primal_0 = _S1680;
    (&_S1938)->differential_0 = 0.0f;
    DiffPair_float_0 _S1939;
    (&_S1939)->primal_0 = _S1681;
    (&_S1939)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1938, &_S1939, k_6);
    float _S1940 = _S1939.differential_0 + _S1684;
    float _S1941 = _S1938.differential_0 + _S1685;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1942;
    (&_S1942)->primal_0 = _S1679;
    (&_S1942)->differential_0 = _S1813;
    s_bwd_length_impl_0(&_S1942, _S1941);
    float2  _S1943 = _S1942.differential_0 + _S1931;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1944;
    (&_S1944)->primal_0 = R_16;
    (&_S1944)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1945;
    (&_S1945)->primal_0 = vert2_5;
    (&_S1945)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S1944, &_S1945, _S1883);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1946;
    (&_S1946)->primal_0 = R_16;
    (&_S1946)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1947;
    (&_S1947)->primal_0 = vert1_5;
    (&_S1947)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S1946, &_S1947, _S1927);
    float3  s_diff_vert0_c_T_1 = make_float3 (_S1943.x, _S1943.y, _S1940);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1948;
    (&_S1948)->primal_0 = R_16;
    (&_S1948)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1949;
    (&_S1949)->primal_0 = vert0_5;
    (&_S1949)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S1948, &_S1949, s_diff_vert0_c_T_1);
    float3  _S1950 = _S1945.differential_0 + _S1792.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1951;
    (&_S1951)->primal_0 = _S1672;
    (&_S1951)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1952;
    (&_S1952)->primal_0 = _S1677;
    (&_S1952)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S1951, &_S1952, _S1950);
    float _S1953 = - _S1952.differential_0.y;
    float _S1954 = _S1676 * _S1952.differential_0.x;
    float _S1955 = - (_S1666 * _S1952.differential_0.x);
    float3  _S1956 = _S1947.differential_0 + _S1791.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1957;
    (&_S1957)->primal_0 = _S1672;
    (&_S1957)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1958;
    (&_S1958)->primal_0 = _S1675;
    (&_S1958)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S1957, &_S1958, _S1956);
    float _S1959 = _S1666 * _S1958.differential_0.x;
    float _S1960 = _S1674 * _S1958.differential_0.x;
    float3  _S1961 = _S1949.differential_0 + _S1849;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1962;
    (&_S1962)->primal_0 = _S1672;
    (&_S1962)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1963;
    (&_S1963)->primal_0 = _S1673;
    (&_S1963)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S1962, &_S1963, _S1961);
    Matrix<float, 3, 3>  _S1964 = transpose_0(_S1951.differential_0 + _S1957.differential_0 + _S1962.differential_0);
    float _S1965 = 2.0f * - _S1964.rows[int(2)].z;
    float _S1966 = 2.0f * _S1964.rows[int(2)].y;
    float _S1967 = 2.0f * _S1964.rows[int(2)].x;
    float _S1968 = 2.0f * _S1964.rows[int(1)].z;
    float _S1969 = 2.0f * - _S1964.rows[int(1)].y;
    float _S1970 = 2.0f * _S1964.rows[int(1)].x;
    float _S1971 = 2.0f * _S1964.rows[int(0)].z;
    float _S1972 = 2.0f * _S1964.rows[int(0)].y;
    float _S1973 = 2.0f * - _S1964.rows[int(0)].x;
    float _S1974 = - _S1970 + _S1972;
    float _S1975 = _S1967 + - _S1971;
    float _S1976 = - _S1966 + _S1968;
    float _S1977 = _S1966 + _S1968;
    float _S1978 = _S1967 + _S1971;
    float _S1979 = _S1970 + _S1972;
    float _S1980 = z_17 * (_S1969 + _S1973);
    float _S1981 = y_19 * (_S1965 + _S1973);
    float _S1982 = x_40 * (_S1965 + _S1969);
    float _S1983 = z_17 * _S1974 + y_19 * _S1975 + x_40 * _S1976;
    float _S1984 = _S1671 * _S1983;
    float _S1985 = w_17 * _S1974 + y_19 * _S1977 + x_40 * _S1978 + _S1980 + _S1980;
    float _S1986 = _S1671 * _S1985;
    float _S1987 = w_17 * _S1975 + z_17 * _S1977 + x_40 * _S1979 + _S1981 + _S1981;
    float _S1988 = _S1671 * _S1987;
    float _S1989 = w_17 * _S1976 + z_17 * _S1978 + y_19 * _S1979 + _S1982 + _S1982;
    float _S1990 = _S1671 * _S1989;
    float _S1991 = quat_17.x * _S1983 + quat_17.w * _S1985 + quat_17.z * _S1987 + quat_17.y * _S1989;
    DiffPair_float_0 _S1992;
    (&_S1992)->primal_0 = _S1670;
    (&_S1992)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1992, _S1991);
    float _S1993 = quat_17.x * _S1992.differential_0;
    float _S1994 = quat_17.w * _S1992.differential_0;
    float _S1995 = quat_17.z * _S1992.differential_0;
    float _S1996 = quat_17.y * _S1992.differential_0;
    float _S1997 = _S1986 + _S1994 + _S1994;
    float _S1998 = _S1988 + _S1995 + _S1995;
    float _S1999 = _S1990 + _S1996 + _S1996;
    float _S2000 = _S1984 + _S1993 + _S1993;
    float _S2001 = _S1955 + _S1959;
    float _S2002 = 0.5f * - _S2001;
    float _S2003 = _S1953 + _S1958.differential_0.y;
    DiffPair_float_0 _S2004;
    (&_S2004)->primal_0 = _S1667;
    (&_S2004)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2004, _S2003);
    float _S2005 = _S2002 + _S2004.differential_0;
    float _S2006 = _S1954 + _S1960 + _S1963.differential_0.x;
    DiffPair_float_0 _S2007;
    (&_S2007)->primal_0 = _S1665;
    (&_S2007)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2007, _S2006);
    float _S2008 = _S2002 + _S2007.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2009;
    (&_S2009)->primal_0 = R_16;
    (&_S2009)->differential_0 = _S1785;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2010;
    (&_S2010)->primal_0 = mean_14;
    (&_S2010)->differential_0 = _S1780;
    s_bwd_prop_mul_3(&_S2009, &_S2010, _S1848);
    float3  _S2011 = _S1883 + _S1927 + s_diff_vert0_c_T_1 + _S1848;
    Matrix<float, 3, 3>  _S2012 = _S1944.differential_0 + _S1946.differential_0 + _S1948.differential_0 + _S2009.differential_0 + _S1788.differential_0;
    float3  _S2013 = make_float3 (_S2008, _S2005, _S2001);
    float4  _S2014 = make_float4 (0.0f);
    *&((&_S2014)->w) = _S1997;
    *&((&_S2014)->z) = _S1998;
    *&((&_S2014)->y) = _S1999;
    *&((&_S2014)->x) = _S2000;
    float4  _S2015 = _S2014;
    float3  _S2016 = _S1950 + _S1956 + _S1961 + _S2010.differential_0;
    *v_mean_4 = _S2016;
    *v_quat_4 = _S2015;
    *v_scale_4 = _S2013;
    *v_hardness_1 = _S1844;
    *v_R_4 = _S2012;
    *v_t_4 = _S2011;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_15, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_17)
{
    DiffPair_float_0 _S2017 = *dpx_15;
    bool _S2018;
    if(((*dpx_15).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S2018 = ((*dpx_15).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S2018 = false;
    }
    float _S2019;
    if(_S2018)
    {
        _S2019 = dOut_17;
    }
    else
    {
        _S2019 = 0.0f;
    }
    dpx_15->primal_0 = _S2017.primal_0;
    dpx_15->differential_0 = _S2019;
    DiffPair_float_0 _S2020 = *dpMin_0;
    if((_S2017.primal_0) < ((*dpMin_0).primal_0))
    {
        _S2019 = dOut_17;
    }
    else
    {
        _S2019 = 0.0f;
    }
    dpMin_0->primal_0 = _S2020.primal_0;
    dpMin_0->differential_0 = _S2019;
    DiffPair_float_0 _S2021 = *dpMax_0;
    if(((*dpx_15).primal_0) > ((*dpMax_0).primal_0))
    {
        _S2019 = dOut_17;
    }
    else
    {
        _S2019 = 0.0f;
    }
    dpMax_0->primal_0 = _S2021.primal_0;
    dpMax_0->differential_0 = _S2019;
    return;
}

inline __device__ float clamp_0(float x_41, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_41), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpy_5, float dOut_18)
{
    if(((*dpx_16).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = 0.0f;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_16).primal_0), ((*dpy_5).primal_0)));
        DiffPair_float_0 _S2022 = *dpx_16;
        float _S2023 = val_0 * (*dpy_5).primal_0 / (*dpx_16).primal_0 * dOut_18;
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = _S2023;
        float _S2024 = val_0 * (F32_log((_S2022.primal_0))) * dOut_18;
        dpy_5->primal_0 = (*dpy_5).primal_0;
        dpy_5->differential_0 = _S2024;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S2025 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S2025))));
    float2  _S2026 = p_0 - v0_0;
    float2  _S2027 = normalize_1(e0_6);
    float2  _S2028 = p_0 - v1_0;
    float2  _S2029 = normalize_1(e1_6);
    float2  _S2030 = p_0 - v2_0;
    float2  _S2031 = normalize_1(e2_2);
    float _S2032 = hardness_6.x;
    float _S2033 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.99500000476837158f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S2026.x * _S2027.y - _S2026.y * _S2027.x)), (se_0 * (_S2028.x * _S2029.y - _S2028.y * _S2029.x))))), (se_0 * (_S2030.x * _S2031.y - _S2030.y * _S2031.x)))) / ((F32_abs((_S2025))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S2033))));
    float _S2034;
    if(a_1 <= 0.0f)
    {
        _S2034 = 0.0f;
    }
    else
    {
        _S2034 = (F32_min(((F32_pow((a_1), (_S2033)))), (0.99900001287460327f)));
    }
    return _S2032 * _S2034;
}

inline __device__ float s_primal_ctx_abs_0(float _S2035)
{
    return (F32_abs((_S2035)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S2036, float _S2037, float _S2038)
{
    return clamp_0(_S2036, _S2037, _S2038);
}

inline __device__ float s_primal_ctx_exp2_0(float _S2039)
{
    return (F32_exp2((_S2039)));
}

inline __device__ float s_primal_ctx_pow_0(float _S2040, float _S2041)
{
    return (F32_pow((_S2040), (_S2041)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S2042, DiffPair_float_0 * _S2043, float _S2044)
{
    _d_pow_0(_S2042, _S2043, _S2044);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S2045, DiffPair_float_0 * _S2046, DiffPair_float_0 * _S2047, float _S2048)
{
    _d_clamp_0(_S2045, _S2046, _S2047, _S2048);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_17, float2  _s_dOut_7)
{
    float _S2049 = length_0((*dpx_17).primal_0);
    float2  _S2050 = (*dpx_17).primal_0 * _s_dOut_7;
    float2  _S2051 = make_float2 (1.0f / _S2049) * _s_dOut_7;
    float _S2052 = - ((_S2050.x + _S2050.y) / (_S2049 * _S2049));
    float2  _S2053 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2054;
    (&_S2054)->primal_0 = (*dpx_17).primal_0;
    (&_S2054)->differential_0 = _S2053;
    s_bwd_length_impl_0(&_S2054, _S2052);
    float2  _S2055 = _S2051 + _S2054.differential_0;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S2055;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2056, float2  _S2057)
{
    s_bwd_prop_normalize_impl_1(_S2056, _S2057);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_8)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S2058 = e0_7.x;
    float _S2059 = e1_7.y;
    float _S2060 = e0_7.y;
    float _S2061 = e1_7.x;
    float _S2062 = _S2058 * _S2059 - _S2060 * _S2061;
    float se_1 = float((F32_sign((_S2062))));
    float2  _S2063 = p_1 - (*dpv0_0).primal_0;
    float2  _S2064 = normalize_1(e0_7);
    float _S2065 = _S2063.x;
    float _S2066 = _S2064.y;
    float _S2067 = _S2063.y;
    float _S2068 = _S2064.x;
    float de0_0 = se_1 * (_S2065 * _S2066 - _S2067 * _S2068);
    float2  _S2069 = p_1 - (*dpv1_0).primal_0;
    float2  _S2070 = normalize_1(e1_7);
    float _S2071 = _S2069.x;
    float _S2072 = _S2070.y;
    float _S2073 = _S2069.y;
    float _S2074 = _S2070.x;
    float de1_0 = se_1 * (_S2071 * _S2072 - _S2073 * _S2074);
    float2  _S2075 = p_1 - (*dpv2_0).primal_0;
    float2  _S2076 = normalize_1(e2_3);
    float _S2077 = _S2075.x;
    float _S2078 = _S2076.y;
    float _S2079 = _S2075.y;
    float _S2080 = _S2076.x;
    float de2_0 = se_1 * (_S2077 * _S2078 - _S2079 * _S2080);
    float _S2081 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S2082 = s_primal_ctx_max_0(_S2081, de2_0);
    float _S2083 = s_primal_ctx_abs_0(_S2062);
    float _S2084 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S2083 / _S2084;
    float _S2085 = _S2084 * _S2084;
    float _S2086 = (*dphardness_0).primal_0.x;
    float _S2087 = (*dphardness_0).primal_0.y;
    float _S2088 = dmax_0 * dmax_0;
    float _S2089 = 1.0f + _S2082 / dmax_0;
    float _S2090 = 1.0f - s_primal_ctx_clamp_0(_S2087, 0.00499999988824129f, 0.99500000476837158f);
    float _S2091 = -1.0f / _S2090;
    float _S2092 = _S2090 * _S2090;
    float _S2093 = 1.0f - s_primal_ctx_exp2_0(_S2091);
    float a_2 = 1.0f - _S2089 * _S2093;
    bool _S2094 = a_2 <= 0.0f;
    float _S2095;
    float _S2096;
    if(_S2094)
    {
        _S2095 = 0.0f;
        _S2096 = 0.0f;
    }
    else
    {
        float _S2097 = s_primal_ctx_pow_0(a_2, _S2090);
        _S2095 = s_primal_ctx_min_0(_S2097, 0.99900001287460327f);
        _S2096 = _S2097;
    }
    float _S2098 = _S2086 * _s_dOut_8;
    float _S2099 = _S2095 * _s_dOut_8;
    if(_S2094)
    {
        _S2095 = 0.0f;
        _S2096 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S2100;
        (&_S2100)->primal_0 = _S2096;
        (&_S2100)->differential_0 = 0.0f;
        DiffPair_float_0 _S2101;
        (&_S2101)->primal_0 = 0.99900001287460327f;
        (&_S2101)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S2100, &_S2101, _S2098);
        DiffPair_float_0 _S2102;
        (&_S2102)->primal_0 = a_2;
        (&_S2102)->differential_0 = 0.0f;
        DiffPair_float_0 _S2103;
        (&_S2103)->primal_0 = _S2090;
        (&_S2103)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S2102, &_S2103, _S2100.differential_0);
        _S2095 = _S2102.differential_0;
        _S2096 = _S2103.differential_0;
    }
    float _S2104 = - _S2095;
    float _S2105 = _S2093 * _S2104;
    float _S2106 = - (_S2089 * _S2104);
    DiffPair_float_0 _S2107;
    (&_S2107)->primal_0 = _S2091;
    (&_S2107)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2107, _S2106);
    float _S2108 = - (-1.0f * - (_S2107.differential_0 / _S2092) + _S2096);
    float _S2109 = _S2105 / _S2088;
    float s_diff_dmax_T_0 = _S2082 * - _S2109;
    float _S2110 = dmax_0 * _S2109;
    DiffPair_float_0 _S2111;
    (&_S2111)->primal_0 = _S2087;
    (&_S2111)->differential_0 = 0.0f;
    DiffPair_float_0 _S2112;
    (&_S2112)->primal_0 = 0.00499999988824129f;
    (&_S2112)->differential_0 = 0.0f;
    DiffPair_float_0 _S2113;
    (&_S2113)->primal_0 = 0.99500000476837158f;
    (&_S2113)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2111, &_S2112, &_S2113, _S2108);
    float _S2114 = s_diff_dmax_T_0 / _S2085;
    float _S2115 = _S2083 * - _S2114;
    float _S2116 = _S2084 * _S2114;
    float2  _S2117 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2118;
    (&_S2118)->primal_0 = e2_3;
    (&_S2118)->differential_0 = _S2117;
    s_bwd_length_impl_0(&_S2118, _S2115);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2119;
    (&_S2119)->primal_0 = e1_7;
    (&_S2119)->differential_0 = _S2117;
    s_bwd_length_impl_0(&_S2119, _S2115);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2120;
    (&_S2120)->primal_0 = e0_7;
    (&_S2120)->differential_0 = _S2117;
    s_bwd_length_impl_0(&_S2120, _S2115);
    DiffPair_float_0 _S2121;
    (&_S2121)->primal_0 = _S2062;
    (&_S2121)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2121, _S2116);
    DiffPair_float_0 _S2122;
    (&_S2122)->primal_0 = _S2081;
    (&_S2122)->differential_0 = 0.0f;
    DiffPair_float_0 _S2123;
    (&_S2123)->primal_0 = de2_0;
    (&_S2123)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2122, &_S2123, _S2110);
    DiffPair_float_0 _S2124;
    (&_S2124)->primal_0 = de0_0;
    (&_S2124)->differential_0 = 0.0f;
    DiffPair_float_0 _S2125;
    (&_S2125)->primal_0 = de1_0;
    (&_S2125)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2124, &_S2125, _S2122.differential_0);
    float _S2126 = se_1 * _S2123.differential_0;
    float _S2127 = - _S2126;
    float _S2128 = _S2080 * _S2127;
    float _S2129 = _S2078 * _S2126;
    float2  _S2130 = make_float2 (_S2079 * _S2127, _S2077 * _S2126);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2131;
    (&_S2131)->primal_0 = e2_3;
    (&_S2131)->differential_0 = _S2117;
    s_bwd_normalize_impl_1(&_S2131, _S2130);
    float2  _S2132 = - make_float2 (_S2129, _S2128);
    float _S2133 = se_1 * _S2125.differential_0;
    float _S2134 = - _S2133;
    float _S2135 = _S2074 * _S2134;
    float _S2136 = _S2072 * _S2133;
    float2  _S2137 = make_float2 (_S2073 * _S2134, _S2071 * _S2133);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2138;
    (&_S2138)->primal_0 = e1_7;
    (&_S2138)->differential_0 = _S2117;
    s_bwd_normalize_impl_1(&_S2138, _S2137);
    float2  _S2139 = - make_float2 (_S2136, _S2135);
    float _S2140 = se_1 * _S2124.differential_0;
    float _S2141 = - _S2140;
    float _S2142 = _S2068 * _S2141;
    float _S2143 = _S2066 * _S2140;
    float2  _S2144 = make_float2 (_S2067 * _S2141, _S2065 * _S2140);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2145;
    (&_S2145)->primal_0 = e0_7;
    (&_S2145)->differential_0 = _S2117;
    s_bwd_normalize_impl_1(&_S2145, _S2144);
    float2  _S2146 = - make_float2 (_S2143, _S2142);
    float _S2147 = - _S2121.differential_0;
    float2  _S2148 = _S2118.differential_0 + _S2131.differential_0;
    float2  _S2149 = - _S2148;
    float2  _S2150 = _S2119.differential_0 + _S2138.differential_0 + make_float2 (_S2060 * _S2147, _S2058 * _S2121.differential_0);
    float2  _S2151 = - _S2150;
    float2  _S2152 = _S2120.differential_0 + _S2145.differential_0 + make_float2 (_S2059 * _S2121.differential_0, _S2061 * _S2147);
    float2  _S2153 = - _S2152;
    float2  _S2154 = make_float2 (_S2099, _S2111.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S2154;
    float2  _S2155 = _S2132 + _S2149 + _S2150;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S2155;
    float2  _S2156 = _S2139 + _S2151 + _S2152;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S2156;
    float2  _S2157 = _S2146 + _S2148 + _S2153;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S2157;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2158, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2159, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2160, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2161, float2  _S2162, float _S2163)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S2158, _S2159, _S2160, _S2161, _S2162, _S2163);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_0, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S2164 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S2164;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S2164;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S2164;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S2164;
    s_bwd_evaluate_alpha_opaque_triangle_fast_0(&dp_v0_0, &dp_v1_0, &dp_v2_0, &dp_hardness_0, p_2, v_alpha_0);
    *v_v0_0 = dp_v0_0.differential_0;
    *v_v1_0 = dp_v2_0.differential_0;
    *v_v2_0 = dp_v1_0.differential_0;
    *v_hardness_2 = dp_hardness_0.differential_0;
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_precise(float2  v0_2, float2  v1_2, float2  v2_2, float2  hardness_8, float2  p_3)
{
    float2  e0_8 = v1_2 - v0_2;
    float2  e1_8 = v2_2 - v1_2;
    float2  e2_4 = v0_2 - v2_2;
    float2  _S2165 = p_3 - v0_2;
    float2  _S2166 = p_3 - v1_2;
    float2  _S2167 = p_3 - v2_2;
    float _S2168 = e0_8.x;
    float _S2169 = e1_8.y;
    float _S2170 = e0_8.y;
    float _S2171 = e1_8.x;
    float _S2172 = _S2168 * _S2169 - _S2170 * _S2171;
    float se_2 = float((F32_sign((_S2172))));
    float _S2173 = hardness_8.x;
    float _S2174 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.99500000476837158f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S2165.x * _S2170 - _S2165.y * _S2168)), (se_2 * (_S2166.x * _S2169 - _S2166.y * _S2171))))), (se_2 * (_S2167.x * e2_4.y - _S2167.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S2165 - e0_8 * make_float2 (clamp_0(dot_1(_S2165, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S2166 - e1_8 * make_float2 (clamp_0(dot_1(_S2166, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S2167 - e2_4 * make_float2 (clamp_0(dot_1(_S2167, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S2172))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S2174))));
    float _S2175;
    if(a_3 <= 0.0f)
    {
        _S2175 = 0.0f;
    }
    else
    {
        _S2175 = (F32_min(((F32_pow((a_3), (_S2174)))), (0.99900001287460327f)));
    }
    return _S2173 * _S2175;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S2176, float2  _S2177)
{
    return dot_1(_S2176, _S2177);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2178, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2179, float _S2180)
{
    _d_dot_1(_S2178, _S2179, _S2180);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_9)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S2181 = p_4 - (*dpv0_1).primal_0;
    float _S2182 = s_primal_ctx_dot_1(_S2181, e0_9);
    float _S2183 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S2184 = _S2182 / _S2183;
    float _S2185 = _S2183 * _S2183;
    float _S2186 = s_primal_ctx_clamp_0(_S2184, 0.0f, 1.0f);
    float2  _S2187 = make_float2 (_S2186);
    float2  _S2188 = _S2181 - e0_9 * make_float2 (_S2186);
    float _S2189 = length_0(_S2188);
    float2  _S2190 = p_4 - (*dpv1_1).primal_0;
    float _S2191 = s_primal_ctx_dot_1(_S2190, e1_9);
    float _S2192 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S2193 = _S2191 / _S2192;
    float _S2194 = _S2192 * _S2192;
    float _S2195 = s_primal_ctx_clamp_0(_S2193, 0.0f, 1.0f);
    float2  _S2196 = make_float2 (_S2195);
    float2  _S2197 = _S2190 - e1_9 * make_float2 (_S2195);
    float _S2198 = length_0(_S2197);
    float2  _S2199 = p_4 - (*dpv2_1).primal_0;
    float _S2200 = s_primal_ctx_dot_1(_S2199, e2_5);
    float _S2201 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S2202 = _S2200 / _S2201;
    float _S2203 = _S2201 * _S2201;
    float _S2204 = s_primal_ctx_clamp_0(_S2202, 0.0f, 1.0f);
    float2  _S2205 = make_float2 (_S2204);
    float2  _S2206 = _S2199 - e2_5 * make_float2 (_S2204);
    float _S2207 = length_0(_S2206);
    float _S2208 = e0_9.x;
    float _S2209 = e1_9.y;
    float _S2210 = e0_9.y;
    float _S2211 = e1_9.x;
    float _S2212 = _S2208 * _S2209 - _S2210 * _S2211;
    float se_3 = float((F32_sign((_S2212))));
    float _S2213 = _S2181.x;
    float _S2214 = _S2181.y;
    float s0_0 = se_3 * (_S2213 * _S2210 - _S2214 * _S2208);
    float _S2215 = _S2190.x;
    float _S2216 = _S2190.y;
    float s1_0 = se_3 * (_S2215 * _S2209 - _S2216 * _S2211);
    float _S2217 = _S2199.x;
    float _S2218 = e2_5.y;
    float _S2219 = _S2199.y;
    float _S2220 = e2_5.x;
    float s2_0 = se_3 * (_S2217 * _S2218 - _S2219 * _S2220);
    float _S2221 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S2221, s2_0)))));
    float _S2222 = s_primal_ctx_min_0(_S2189, _S2198);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S2222, _S2207);
    float _S2223 = s_primal_ctx_abs_0(_S2212);
    float _S2224 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S2223 / _S2224;
    float _S2225 = _S2224 * _S2224;
    float _S2226 = (*dphardness_1).primal_0.x;
    float _S2227 = (*dphardness_1).primal_0.y;
    float _S2228 = dmax_1 * dmax_1;
    float _S2229 = 1.0f + dv_0 / dmax_1;
    float _S2230 = 1.0f - s_primal_ctx_clamp_0(_S2227, 0.00499999988824129f, 0.99500000476837158f);
    float _S2231 = -1.0f / _S2230;
    float _S2232 = _S2230 * _S2230;
    float _S2233 = 1.0f - s_primal_ctx_exp2_0(_S2231);
    float a_4 = 1.0f - _S2229 * _S2233;
    bool _S2234 = a_4 <= 0.0f;
    float _S2235;
    float _S2236;
    if(_S2234)
    {
        _S2235 = 0.0f;
        _S2236 = 0.0f;
    }
    else
    {
        float _S2237 = s_primal_ctx_pow_0(a_4, _S2230);
        _S2235 = s_primal_ctx_min_0(_S2237, 0.99900001287460327f);
        _S2236 = _S2237;
    }
    float _S2238 = _S2226 * _s_dOut_9;
    float _S2239 = _S2235 * _s_dOut_9;
    if(_S2234)
    {
        _S2235 = 0.0f;
        _S2236 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S2240;
        (&_S2240)->primal_0 = _S2236;
        (&_S2240)->differential_0 = 0.0f;
        DiffPair_float_0 _S2241;
        (&_S2241)->primal_0 = 0.99900001287460327f;
        (&_S2241)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S2240, &_S2241, _S2238);
        DiffPair_float_0 _S2242;
        (&_S2242)->primal_0 = a_4;
        (&_S2242)->differential_0 = 0.0f;
        DiffPair_float_0 _S2243;
        (&_S2243)->primal_0 = _S2230;
        (&_S2243)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S2242, &_S2243, _S2240.differential_0);
        _S2235 = _S2242.differential_0;
        _S2236 = _S2243.differential_0;
    }
    float _S2244 = - _S2235;
    float _S2245 = _S2233 * _S2244;
    float _S2246 = - (_S2229 * _S2244);
    DiffPair_float_0 _S2247;
    (&_S2247)->primal_0 = _S2231;
    (&_S2247)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2247, _S2246);
    float _S2248 = - (-1.0f * - (_S2247.differential_0 / _S2232) + _S2236);
    float _S2249 = _S2245 / _S2228;
    float s_diff_dmax_T_1 = dv_0 * - _S2249;
    float s_diff_dv_T_0 = dmax_1 * _S2249;
    DiffPair_float_0 _S2250;
    (&_S2250)->primal_0 = _S2227;
    (&_S2250)->differential_0 = 0.0f;
    DiffPair_float_0 _S2251;
    (&_S2251)->primal_0 = 0.00499999988824129f;
    (&_S2251)->differential_0 = 0.0f;
    DiffPair_float_0 _S2252;
    (&_S2252)->primal_0 = 0.99500000476837158f;
    (&_S2252)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2250, &_S2251, &_S2252, _S2248);
    float _S2253 = s_diff_dmax_T_1 / _S2225;
    float _S2254 = _S2223 * - _S2253;
    float _S2255 = _S2224 * _S2253;
    float2  _S2256 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2257;
    (&_S2257)->primal_0 = e2_5;
    (&_S2257)->differential_0 = _S2256;
    s_bwd_length_impl_0(&_S2257, _S2254);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2258;
    (&_S2258)->primal_0 = e1_9;
    (&_S2258)->differential_0 = _S2256;
    s_bwd_length_impl_0(&_S2258, _S2254);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2259;
    (&_S2259)->primal_0 = e0_9;
    (&_S2259)->differential_0 = _S2256;
    s_bwd_length_impl_0(&_S2259, _S2254);
    DiffPair_float_0 _S2260;
    (&_S2260)->primal_0 = _S2212;
    (&_S2260)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2260, _S2255);
    float _S2261 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S2262;
    (&_S2262)->primal_0 = _S2222;
    (&_S2262)->differential_0 = 0.0f;
    DiffPair_float_0 _S2263;
    (&_S2263)->primal_0 = _S2207;
    (&_S2263)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2262, &_S2263, _S2261);
    DiffPair_float_0 _S2264;
    (&_S2264)->primal_0 = _S2189;
    (&_S2264)->differential_0 = 0.0f;
    DiffPair_float_0 _S2265;
    (&_S2265)->primal_0 = _S2198;
    (&_S2265)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2264, &_S2265, _S2262.differential_0);
    DiffPair_float_0 _S2266;
    (&_S2266)->primal_0 = _S2221;
    (&_S2266)->differential_0 = 0.0f;
    DiffPair_float_0 _S2267;
    (&_S2267)->primal_0 = s2_0;
    (&_S2267)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2266, &_S2267, 0.0f);
    DiffPair_float_0 _S2268;
    (&_S2268)->primal_0 = s0_0;
    (&_S2268)->differential_0 = 0.0f;
    DiffPair_float_0 _S2269;
    (&_S2269)->primal_0 = s1_0;
    (&_S2269)->differential_0 = 0.0f;
    s_bwd_prop_max_0(&_S2268, &_S2269, _S2266.differential_0);
    float _S2270 = se_3 * _S2267.differential_0;
    float _S2271 = - _S2270;
    float _S2272 = _S2219 * _S2271;
    float _S2273 = _S2220 * _S2271;
    float _S2274 = _S2217 * _S2270;
    float _S2275 = _S2218 * _S2270;
    float _S2276 = se_3 * _S2269.differential_0;
    float _S2277 = - _S2276;
    float _S2278 = _S2211 * _S2277;
    float _S2279 = _S2209 * _S2276;
    float _S2280 = se_3 * _S2268.differential_0;
    float _S2281 = - _S2280;
    float _S2282 = _S2208 * _S2281;
    float _S2283 = _S2210 * _S2280;
    float _S2284 = - _S2260.differential_0;
    float _S2285 = _S2216 * _S2277 + _S2210 * _S2284;
    float _S2286 = _S2213 * _S2280 + _S2211 * _S2284;
    float _S2287 = _S2215 * _S2276 + _S2208 * _S2260.differential_0;
    float _S2288 = _S2214 * _S2281 + _S2209 * _S2260.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2289;
    (&_S2289)->primal_0 = _S2206;
    (&_S2289)->differential_0 = _S2256;
    s_bwd_length_impl_0(&_S2289, _S2263.differential_0);
    float2  _S2290 = - _S2289.differential_0;
    float2  _S2291 = e2_5 * _S2290;
    float2  _S2292 = _S2205 * _S2290;
    float _S2293 = _S2291.x + _S2291.y;
    DiffPair_float_0 _S2294;
    (&_S2294)->primal_0 = _S2202;
    (&_S2294)->differential_0 = 0.0f;
    DiffPair_float_0 _S2295;
    (&_S2295)->primal_0 = 0.0f;
    (&_S2295)->differential_0 = 0.0f;
    DiffPair_float_0 _S2296;
    (&_S2296)->primal_0 = 1.0f;
    (&_S2296)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2294, &_S2295, &_S2296, _S2293);
    float _S2297 = _S2294.differential_0 / _S2203;
    float _S2298 = _S2200 * - _S2297;
    float _S2299 = _S2201 * _S2297;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2300;
    (&_S2300)->primal_0 = e2_5;
    (&_S2300)->differential_0 = _S2256;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2301;
    (&_S2301)->primal_0 = e2_5;
    (&_S2301)->differential_0 = _S2256;
    s_bwd_prop_dot_1(&_S2300, &_S2301, _S2298);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2302;
    (&_S2302)->primal_0 = _S2199;
    (&_S2302)->differential_0 = _S2256;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2303;
    (&_S2303)->primal_0 = e2_5;
    (&_S2303)->differential_0 = _S2256;
    s_bwd_prop_dot_1(&_S2302, &_S2303, _S2299);
    float2  _S2304 = - (_S2289.differential_0 + _S2302.differential_0 + make_float2 (_S2275, _S2273));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2305;
    (&_S2305)->primal_0 = _S2197;
    (&_S2305)->differential_0 = _S2256;
    s_bwd_length_impl_0(&_S2305, _S2265.differential_0);
    float2  _S2306 = - _S2305.differential_0;
    float2  _S2307 = e1_9 * _S2306;
    float2  _S2308 = _S2196 * _S2306;
    float _S2309 = _S2307.x + _S2307.y;
    DiffPair_float_0 _S2310;
    (&_S2310)->primal_0 = _S2193;
    (&_S2310)->differential_0 = 0.0f;
    DiffPair_float_0 _S2311;
    (&_S2311)->primal_0 = 0.0f;
    (&_S2311)->differential_0 = 0.0f;
    DiffPair_float_0 _S2312;
    (&_S2312)->primal_0 = 1.0f;
    (&_S2312)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2310, &_S2311, &_S2312, _S2309);
    float _S2313 = _S2310.differential_0 / _S2194;
    float _S2314 = _S2191 * - _S2313;
    float _S2315 = _S2192 * _S2313;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2316;
    (&_S2316)->primal_0 = e1_9;
    (&_S2316)->differential_0 = _S2256;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2317;
    (&_S2317)->primal_0 = e1_9;
    (&_S2317)->differential_0 = _S2256;
    s_bwd_prop_dot_1(&_S2316, &_S2317, _S2314);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2318;
    (&_S2318)->primal_0 = _S2190;
    (&_S2318)->differential_0 = _S2256;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2319;
    (&_S2319)->primal_0 = e1_9;
    (&_S2319)->differential_0 = _S2256;
    s_bwd_prop_dot_1(&_S2318, &_S2319, _S2315);
    float2  _S2320 = - (_S2305.differential_0 + _S2318.differential_0 + make_float2 (_S2279, _S2278));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2321;
    (&_S2321)->primal_0 = _S2188;
    (&_S2321)->differential_0 = _S2256;
    s_bwd_length_impl_0(&_S2321, _S2264.differential_0);
    float2  _S2322 = - _S2321.differential_0;
    float2  _S2323 = e0_9 * _S2322;
    float2  _S2324 = _S2187 * _S2322;
    float _S2325 = _S2323.x + _S2323.y;
    DiffPair_float_0 _S2326;
    (&_S2326)->primal_0 = _S2184;
    (&_S2326)->differential_0 = 0.0f;
    DiffPair_float_0 _S2327;
    (&_S2327)->primal_0 = 0.0f;
    (&_S2327)->differential_0 = 0.0f;
    DiffPair_float_0 _S2328;
    (&_S2328)->primal_0 = 1.0f;
    (&_S2328)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2326, &_S2327, &_S2328, _S2325);
    float _S2329 = _S2326.differential_0 / _S2185;
    float _S2330 = _S2182 * - _S2329;
    float _S2331 = _S2183 * _S2329;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2332;
    (&_S2332)->primal_0 = e0_9;
    (&_S2332)->differential_0 = _S2256;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2333;
    (&_S2333)->primal_0 = e0_9;
    (&_S2333)->differential_0 = _S2256;
    s_bwd_prop_dot_1(&_S2332, &_S2333, _S2330);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2334;
    (&_S2334)->primal_0 = _S2181;
    (&_S2334)->differential_0 = _S2256;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2335;
    (&_S2335)->primal_0 = e0_9;
    (&_S2335)->differential_0 = _S2256;
    s_bwd_prop_dot_1(&_S2334, &_S2335, _S2331);
    float2  _S2336 = - (_S2321.differential_0 + _S2334.differential_0 + make_float2 (_S2283, _S2282));
    float2  _S2337 = _S2257.differential_0 + _S2292 + _S2301.differential_0 + _S2300.differential_0 + _S2303.differential_0 + make_float2 (_S2272, _S2274);
    float2  _S2338 = - _S2337;
    float2  _S2339 = _S2258.differential_0 + _S2308 + _S2317.differential_0 + _S2316.differential_0 + _S2319.differential_0 + make_float2 (_S2285, _S2287);
    float2  _S2340 = - _S2339;
    float2  _S2341 = _S2259.differential_0 + _S2324 + _S2333.differential_0 + _S2332.differential_0 + _S2335.differential_0 + make_float2 (_S2288, _S2286);
    float2  _S2342 = - _S2341;
    float2  _S2343 = make_float2 (_S2239, _S2250.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S2343;
    float2  _S2344 = _S2304 + _S2338 + _S2339;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S2344;
    float2  _S2345 = _S2320 + _S2340 + _S2341;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S2345;
    float2  _S2346 = _S2336 + _S2337 + _S2342;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S2346;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2347, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2348, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2349, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2350, float2  _S2351, float _S2352)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S2347, _S2348, _S2349, _S2350, _S2351, _S2352);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_1, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S2353 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S2353;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S2353;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S2353;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S2353;
    s_bwd_evaluate_alpha_opaque_triangle_precise_0(&dp_v0_1, &dp_v1_1, &dp_v2_1, &dp_hardness_1, p_5, v_alpha_1);
    *v_v0_1 = dp_v0_1.differential_0;
    *v_v1_1 = dp_v2_1.differential_0;
    *v_v2_1 = dp_v1_1.differential_0;
    *v_hardness_3 = dp_hardness_1.differential_0;
    return;
}

