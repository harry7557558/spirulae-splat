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

inline __device__ CameraDistortion_0 CameraDistortion_x24_syn_dadd_0(CameraDistortion_0 * SLANG_anonymous_0_0, CameraDistortion_0 * SLANG_anonymous_1_0)
{
    CameraDistortion_0 result_1;
    (&result_1)->radial_coeffs_0 = SLANG_anonymous_0_0->radial_coeffs_0 + SLANG_anonymous_1_0->radial_coeffs_0;
    (&result_1)->tangential_coeffs_0 = SLANG_anonymous_0_0->tangential_coeffs_0 + SLANG_anonymous_1_0->tangential_coeffs_0;
    (&result_1)->thin_prism_coeffs_0 = SLANG_anonymous_0_0->thin_prism_coeffs_0 + SLANG_anonymous_1_0->thin_prism_coeffs_0;
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
    float _S125 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S125;
    float _S126 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S126;
    float det_blur_0 = _S125 * _S126 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_8, float dOut_11)
{
    float _S127 = (F32_exp(((*dpx_8).primal_0))) * dOut_11;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S127;
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
    float3  _S128 = exp_0((*dpx_9).primal_0) * dOut_12;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S128;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_10, float dOut_13)
{
    float _S129 = 1.0f / (*dpx_10).primal_0 * dOut_13;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S129;
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_5, float3  dOut_14)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_11).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpy_5).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_max_0(&left_dp_0, &right_dp_0, dOut_14.x);
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = left_dp_0.differential_0;
    float3  right_d_result_4;
    *&((&right_d_result_4)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_11).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_5).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_14.y);
    *&((&left_d_result_4)->y) = left_dp_1.differential_0;
    *&((&right_d_result_4)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_11).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_5).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_14.z);
    *&((&left_d_result_4)->z) = left_dp_2.differential_0;
    *&((&right_d_result_4)->z) = right_dp_2.differential_0;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = left_d_result_4;
    dpy_5->primal_0 = (*dpy_5).primal_0;
    dpy_5->differential_0 = right_d_result_4;
    return;
}

inline __device__ float3  max_0(float3  x_14, float3  y_5)
{
    float3  result_14;
    int i_7 = int(0);
    for(;;)
    {
        if(i_7 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_14, i_7) = (F32_max((_slang_vector_get_element(x_14, i_7)), (_slang_vector_get_element(y_5, i_7))));
        i_7 = i_7 + int(1);
    }
    return result_14;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_2, float3  t_1, float fx_4, float fy_4, float cx_4, float cy_4, float4  radial_coeffs_3, float2  tangential_coeffs_3, float2  thin_prism_coeffs_3, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_4, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_2, mean_0) + t_1;
        float _S130 = mean_c_0.z;
        bool _S131;
        if(_S130 < near_plane_0)
        {
            _S131 = true;
        }
        else
        {
            _S131 = _S130 > far_plane_0;
        }
        if(_S131)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S132 = exp_0(scale_2);
        float x_15 = quat_3.y;
        float inv_norm_3 = (F32_rsqrt((x_15 * x_15 + quat_3.z * quat_3.z + quat_3.w * quat_3.w + quat_3.x * quat_3.x)));
        float x_16 = quat_3.y * inv_norm_3;
        float y_6 = quat_3.z * inv_norm_3;
        float z_3 = quat_3.w * inv_norm_3;
        float w_3 = quat_3.x * inv_norm_3;
        float x2_3 = x_16 * x_16;
        float y2_3 = y_6 * y_6;
        float z2_3 = z_3 * z_3;
        float xy_3 = x_16 * y_6;
        float xz_3 = x_16 * z_3;
        float yz_3 = y_6 * z_3;
        float wx_3 = w_3 * x_16;
        float wy_3 = w_3 * y_6;
        float wz_3 = w_3 * z_3;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S132.x, 0.0f, 0.0f, 0.0f, _S132.y, 0.0f, 0.0f, 0.0f, _S132.z));
        Matrix<float, 3, 3>  _S133 = transpose_0(R_2);
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_2, mul_4(M_2, transpose_0(M_2))), _S133);
        Matrix<float, 2, 2>  covar2d_0;
        float _S134 = float(image_width_0);
        float _S135 = float(image_height_0);
        float _S136 = 0.30000001192092896f * (0.5f * _S134 / fx_4);
        float _S137 = 0.30000001192092896f * (0.5f * _S135 / fy_4);
        float rz_2 = 1.0f / mean_c_0.z;
        float rz2_2 = rz_2 * rz_2;
        Matrix<float, 2, 3>  J_4 = makeMatrix<float, 2, 3> (fx_4 * rz_2, 0.0f, - fx_4 * (mean_c_0.z * (F32_min(((_S134 - cx_4) / fx_4 + _S136), ((F32_max((- (cx_4 / fx_4 + _S136)), (mean_c_0.x * rz_2))))))) * rz2_2, 0.0f, fy_4 * rz_2, - fy_4 * (mean_c_0.z * (F32_min(((_S135 - cy_4) / fy_4 + _S137), ((F32_max((- (cy_4 / fy_4 + _S137)), (mean_c_0.y * rz_2))))))) * rz2_2);
        covar2d_0 = mul_6(mul_5(J_4, covar_c_0), transpose_1(J_4));
        *mean2d_4 = make_float2 (fx_4 * mean_c_0.x * rz_2 + cx_4, fy_4 * mean_c_0.y * rz_2 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S138 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S138;
        float _S139 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S139;
        float det_blur_1 = _S138 * _S139 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S140 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S131 = true;
        }
        else
        {
            _S131 = xmin_0 >= _S134;
        }
        if(_S131)
        {
            _S131 = true;
        }
        else
        {
            _S131 = ymax_0 <= 0.0f;
        }
        if(_S131)
        {
            _S131 = true;
        }
        else
        {
            _S131 = ymin_0 >= _S135;
        }
        if(_S131)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = _S130;
        *conic_0 = make_float3 (_S140.rows[int(0)].x, _S140.rows[int(0)].y, _S140.rows[int(1)].y);
        float3  _S141 = mean_0 - - mul_0(_S133, t_1);
        float3  _S142 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S142;
        float _S143 = _S141.x;
        float _S144 = _S141.y;
        float _S145 = _S141.z;
        float norm_0 = (F32_sqrt((_S143 * _S143 + _S144 * _S144 + _S145 * _S145)));
        float x_17 = _S143 / norm_0;
        float y_7 = _S144 / norm_0;
        float z_4 = _S145 / norm_0;
        float3  _S146 = _S142 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_0)[int(1)] + make_float3 (z_4) * (*sh_coeffs_0)[int(2)] - make_float3 (x_17) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S146;
        float z2_4 = z_4 * z_4;
        float fTmp0B_0 = -1.09254848957061768f * z_4;
        float fC1_0 = x_17 * x_17 - y_7 * y_7;
        float fS1_0 = 2.0f * x_17 * y_7;
        float3  _S147 = _S146 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_7) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_4 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_17) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S147;
        float fTmp0C_0 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_4;
        *rgb_0 = max_0(_S147 + (make_float3 (-0.59004360437393188f * (x_17 * fS1_0 + y_7 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_7) * (*sh_coeffs_0)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_4 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_17) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_0 - y_7 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_3, float3  t_2, float fx_5, float fy_5, float cx_5, float cy_5, float4  radial_coeffs_4, float2  tangential_coeffs_4, float2  thin_prism_coeffs_4, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_5, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_3, mean_1) + t_2;
        float _S148 = mean_c_1.z;
        bool _S149;
        if(_S148 < near_plane_1)
        {
            _S149 = true;
        }
        else
        {
            _S149 = _S148 > far_plane_1;
        }
        if(_S149)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S150 = exp_0(scale_3);
        float x_18 = quat_4.y;
        float inv_norm_4 = (F32_rsqrt((x_18 * x_18 + quat_4.z * quat_4.z + quat_4.w * quat_4.w + quat_4.x * quat_4.x)));
        float x_19 = quat_4.y * inv_norm_4;
        float y_8 = quat_4.z * inv_norm_4;
        float z_5 = quat_4.w * inv_norm_4;
        float w_4 = quat_4.x * inv_norm_4;
        float x2_4 = x_19 * x_19;
        float y2_4 = y_8 * y_8;
        float z2_5 = z_5 * z_5;
        float xy_4 = x_19 * y_8;
        float xz_4 = x_19 * z_5;
        float yz_4 = y_8 * z_5;
        float wx_4 = w_4 * x_19;
        float wy_4 = w_4 * y_8;
        float wz_4 = w_4 * z_5;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_5), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_5), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S150.x, 0.0f, 0.0f, 0.0f, _S150.y, 0.0f, 0.0f, 0.0f, _S150.z));
        Matrix<float, 3, 3>  _S151 = transpose_0(R_3);
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_3, mul_4(M_3, transpose_0(M_3))), _S151);
        Matrix<float, 2, 2>  covar2d_1;
        fisheye_proj_3dgs(mean_c_1, covar_c_1, fx_5, fy_5, cx_5, cy_5, radial_coeffs_4, tangential_coeffs_4, thin_prism_coeffs_4, &covar2d_1, mean2d_5);
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S152 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S152;
        float _S153 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S153;
        float det_blur_2 = _S152 * _S153 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S154 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S149 = true;
        }
        else
        {
            _S149 = xmin_1 >= float(image_width_1);
        }
        if(_S149)
        {
            _S149 = true;
        }
        else
        {
            _S149 = ymax_1 <= 0.0f;
        }
        if(_S149)
        {
            _S149 = true;
        }
        else
        {
            _S149 = ymin_1 >= float(image_height_1);
        }
        if(_S149)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = _S148;
        *conic_1 = make_float3 (_S154.rows[int(0)].x, _S154.rows[int(0)].y, _S154.rows[int(1)].y);
        float3  _S155 = mean_1 - - mul_0(_S151, t_2);
        float3  _S156 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S156;
        float _S157 = _S155.x;
        float _S158 = _S155.y;
        float _S159 = _S155.z;
        float norm_1 = (F32_sqrt((_S157 * _S157 + _S158 * _S158 + _S159 * _S159)));
        float x_20 = _S157 / norm_1;
        float y_9 = _S158 / norm_1;
        float z_6 = _S159 / norm_1;
        float3  _S160 = _S156 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_1)[int(1)] + make_float3 (z_6) * (*sh_coeffs_1)[int(2)] - make_float3 (x_20) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S160;
        float z2_6 = z_6 * z_6;
        float fTmp0B_1 = -1.09254848957061768f * z_6;
        float fC1_1 = x_20 * x_20 - y_9 * y_9;
        float fS1_1 = 2.0f * x_20 * y_9;
        float3  _S161 = _S160 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_9) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_6 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_20) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S161;
        float fTmp0C_1 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_6;
        *rgb_1 = max_0(_S161 + (make_float3 (-0.59004360437393188f * (x_20 * fS1_1 + y_9 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_9) * (*sh_coeffs_1)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_6 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_20) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_20 * fC1_1 - y_9 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_6, float fy_6, float cx_6, float cy_6, float4  radial_coeffs_5, float2  tangential_coeffs_5, float2  thin_prism_coeffs_5, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_6, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_4, mean_2) + t_3;
        float _S162 = mean_c_2.z;
        bool _S163;
        if(_S162 < near_plane_2)
        {
            _S163 = true;
        }
        else
        {
            _S163 = _S162 > far_plane_2;
        }
        if(_S163)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S164 = exp_0(scale_4);
        float x_21 = quat_5.y;
        float inv_norm_5 = (F32_rsqrt((x_21 * x_21 + quat_5.z * quat_5.z + quat_5.w * quat_5.w + quat_5.x * quat_5.x)));
        float x_22 = quat_5.y * inv_norm_5;
        float y_10 = quat_5.z * inv_norm_5;
        float z_7 = quat_5.w * inv_norm_5;
        float w_5 = quat_5.x * inv_norm_5;
        float x2_5 = x_22 * x_22;
        float y2_5 = y_10 * y_10;
        float z2_7 = z_7 * z_7;
        float xy_5 = x_22 * y_10;
        float xz_5 = x_22 * z_7;
        float yz_5 = y_10 * z_7;
        float wx_5 = w_5 * x_22;
        float wy_5 = w_5 * y_10;
        float wz_5 = w_5 * z_7;
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_7), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_7), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S164.x, 0.0f, 0.0f, 0.0f, _S164.y, 0.0f, 0.0f, 0.0f, _S164.z));
        Matrix<float, 3, 3>  _S165 = transpose_0(R_4);
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_4, mul_4(M_4, transpose_0(M_4))), _S165);
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_5, covar_c_2), transpose_1(J_5));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S166 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S166;
        float _S167 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S167;
        float det_blur_3 = _S166 * _S167 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S168 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S163 = true;
        }
        else
        {
            _S163 = xmin_2 >= float(image_width_2);
        }
        if(_S163)
        {
            _S163 = true;
        }
        else
        {
            _S163 = ymax_2 <= 0.0f;
        }
        if(_S163)
        {
            _S163 = true;
        }
        else
        {
            _S163 = ymin_2 >= float(image_height_2);
        }
        if(_S163)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = _S162;
        *conic_2 = make_float3 (_S168.rows[int(0)].x, _S168.rows[int(0)].y, _S168.rows[int(1)].y);
        float3  _S169 = mean_2 - - mul_0(_S165, t_3);
        float3  _S170 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S170;
        float _S171 = _S169.x;
        float _S172 = _S169.y;
        float _S173 = _S169.z;
        float norm_2 = (F32_sqrt((_S171 * _S171 + _S172 * _S172 + _S173 * _S173)));
        float x_23 = _S171 / norm_2;
        float y_11 = _S172 / norm_2;
        float z_8 = _S173 / norm_2;
        float3  _S174 = _S170 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_2)[int(1)] + make_float3 (z_8) * (*sh_coeffs_2)[int(2)] - make_float3 (x_23) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S174;
        float z2_8 = z_8 * z_8;
        float fTmp0B_2 = -1.09254848957061768f * z_8;
        float fC1_2 = x_23 * x_23 - y_11 * y_11;
        float fS1_2 = 2.0f * x_23 * y_11;
        float3  _S175 = _S174 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_11) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_8 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_23) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S175;
        float fTmp0C_2 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_8;
        *rgb_2 = max_0(_S175 + (make_float3 (-0.59004360437393188f * (x_23 * fS1_2 + y_11 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_11) * (*sh_coeffs_2)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_8 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_23) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_23 * fC1_2 - y_11 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_7, float fy_7, float cx_7, float cy_7, float4  radial_coeffs_6, float2  tangential_coeffs_6, float2  thin_prism_coeffs_6, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_7, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    float3  mean_c_3 = mul_0(R_5, mean_3) + t_4;
    float3  _S176 = exp_0(scale_5);
    float x_24 = quat_6.y;
    float inv_norm_6 = (F32_rsqrt((x_24 * x_24 + quat_6.z * quat_6.z + quat_6.w * quat_6.w + quat_6.x * quat_6.x)));
    float x_25 = quat_6.y * inv_norm_6;
    float y_12 = quat_6.z * inv_norm_6;
    float z_9 = quat_6.w * inv_norm_6;
    float w_6 = quat_6.x * inv_norm_6;
    float x2_6 = x_25 * x_25;
    float y2_6 = y_12 * y_12;
    float z2_9 = z_9 * z_9;
    float xy_6 = x_25 * y_12;
    float xz_6 = x_25 * z_9;
    float yz_6 = y_12 * z_9;
    float wx_6 = w_6 * x_25;
    float wy_6 = w_6 * y_12;
    float wz_6 = w_6 * z_9;
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_9), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_9), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S176.x, 0.0f, 0.0f, 0.0f, _S176.y, 0.0f, 0.0f, 0.0f, _S176.z));
    Matrix<float, 3, 3>  _S177 = transpose_0(R_5);
    float _S178 = float(image_width_3);
    float _S179 = float(image_height_3);
    float _S180 = 0.30000001192092896f * (0.5f * _S178 / fx_7);
    float _S181 = 0.30000001192092896f * (0.5f * _S179 / fy_7);
    float rz_3 = 1.0f / mean_c_3.z;
    float rz2_3 = rz_3 * rz_3;
    Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, - fx_7 * (mean_c_3.z * (F32_min(((_S178 - cx_7) / fx_7 + _S180), ((F32_max((- (cx_7 / fx_7 + _S180)), (mean_c_3.x * rz_3))))))) * rz2_3, 0.0f, fy_7 * rz_3, - fy_7 * (mean_c_3.z * (F32_min(((_S179 - cy_7) / fy_7 + _S181), ((F32_max((- (cy_7 / fy_7 + _S181)), (mean_c_3.y * rz_3))))))) * rz2_3);
    Matrix<float, 2, 2>  covar2d_3 = mul_6(mul_5(J_6, mul_4(mul_4(R_5, mul_4(M_5, transpose_0(M_5))), _S177)), transpose_1(J_6));
    *mean2d_7 = make_float2 (fx_7 * mean_c_3.x * rz_3 + cx_7, fy_7 * mean_c_3.y * rz_3 + cy_7);
    float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
    float _S182 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_3)->rows + (int(0)))->x) = _S182;
    float _S183 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_3)->rows + (int(1)))->y) = _S183;
    float compensation_4 = (F32_sqrt(((F32_max((0.0f), (det_orig_4 / (_S182 * _S183 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_3.rows[int(0)].x * covar2d_3.rows[int(1)].y - covar2d_3.rows[int(0)].y * covar2d_3.rows[int(1)].x);
    Matrix<float, 2, 2>  _S184 = makeMatrix<float, 2, 2> (covar2d_3.rows[int(1)].y * invdet_4, - covar2d_3.rows[int(0)].y * invdet_4, - covar2d_3.rows[int(1)].x * invdet_4, covar2d_3.rows[int(0)].x * invdet_4);
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
    *conic_3 = make_float3 (_S184.rows[int(0)].x, _S184.rows[int(0)].y, _S184.rows[int(1)].y);
    float3  _S185 = mean_3 - - mul_0(_S177, t_4);
    float3  _S186 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
    *rgb_3 = _S186;
    float _S187 = _S185.x;
    float _S188 = _S185.y;
    float _S189 = _S185.z;
    float norm_3 = (F32_sqrt((_S187 * _S187 + _S188 * _S188 + _S189 * _S189)));
    float x_26 = _S187 / norm_3;
    float y_13 = _S188 / norm_3;
    float z_10 = _S189 / norm_3;
    float3  _S190 = _S186 + make_float3 (0.48860251903533936f) * (make_float3 (- y_13) * (*sh_coeffs_3)[int(1)] + make_float3 (z_10) * (*sh_coeffs_3)[int(2)] - make_float3 (x_26) * (*sh_coeffs_3)[int(3)]);
    *rgb_3 = _S190;
    float z2_10 = z_10 * z_10;
    float fTmp0B_3 = -1.09254848957061768f * z_10;
    float fC1_3 = x_26 * x_26 - y_13 * y_13;
    float fS1_3 = 2.0f * x_26 * y_13;
    float3  _S191 = _S190 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_13) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_10 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_26) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
    *rgb_3 = _S191;
    float fTmp0C_3 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_10;
    *rgb_3 = max_0(_S191 + (make_float3 (-0.59004360437393188f * (x_26 * fS1_3 + y_13 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_13) * (*sh_coeffs_3)[int(11)] + make_float3 (z_10 * (1.86588168144226074f * z2_10 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_26) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_26 * fC1_3 - y_13 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_8, float fy_8, float cx_8, float cy_8, float4  radial_coeffs_7, float2  tangential_coeffs_7, float2  thin_prism_coeffs_7, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_8, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    float3  mean_c_4 = mul_0(R_6, mean_4) + t_5;
    float3  _S192 = exp_0(scale_6);
    float x_27 = quat_7.y;
    float inv_norm_7 = (F32_rsqrt((x_27 * x_27 + quat_7.z * quat_7.z + quat_7.w * quat_7.w + quat_7.x * quat_7.x)));
    float x_28 = quat_7.y * inv_norm_7;
    float y_14 = quat_7.z * inv_norm_7;
    float z_11 = quat_7.w * inv_norm_7;
    float w_7 = quat_7.x * inv_norm_7;
    float x2_7 = x_28 * x_28;
    float y2_7 = y_14 * y_14;
    float z2_11 = z_11 * z_11;
    float xy_7 = x_28 * y_14;
    float xz_7 = x_28 * z_11;
    float yz_7 = y_14 * z_11;
    float wx_7 = w_7 * x_28;
    float wy_7 = w_7 * y_14;
    float wz_7 = w_7 * z_11;
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_11), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_11), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S192.x, 0.0f, 0.0f, 0.0f, _S192.y, 0.0f, 0.0f, 0.0f, _S192.z));
    Matrix<float, 3, 3>  _S193 = transpose_0(R_6);
    Matrix<float, 2, 2>  covar2d_4;
    fisheye_proj_3dgs(mean_c_4, mul_4(mul_4(R_6, mul_4(M_6, transpose_0(M_6))), _S193), fx_8, fy_8, cx_8, cy_8, radial_coeffs_7, tangential_coeffs_7, thin_prism_coeffs_7, &covar2d_4, mean2d_8);
    float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
    float _S194 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_4)->rows + (int(0)))->x) = _S194;
    float _S195 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_4)->rows + (int(1)))->y) = _S195;
    float compensation_5 = (F32_sqrt(((F32_max((0.0f), (det_orig_5 / (_S194 * _S195 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_4.rows[int(0)].x * covar2d_4.rows[int(1)].y - covar2d_4.rows[int(0)].y * covar2d_4.rows[int(1)].x);
    Matrix<float, 2, 2>  _S196 = makeMatrix<float, 2, 2> (covar2d_4.rows[int(1)].y * invdet_5, - covar2d_4.rows[int(0)].y * invdet_5, - covar2d_4.rows[int(1)].x * invdet_5, covar2d_4.rows[int(0)].x * invdet_5);
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
    *conic_4 = make_float3 (_S196.rows[int(0)].x, _S196.rows[int(0)].y, _S196.rows[int(1)].y);
    float3  _S197 = mean_4 - - mul_0(_S193, t_5);
    float3  _S198 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
    *rgb_4 = _S198;
    float _S199 = _S197.x;
    float _S200 = _S197.y;
    float _S201 = _S197.z;
    float norm_4 = (F32_sqrt((_S199 * _S199 + _S200 * _S200 + _S201 * _S201)));
    float x_29 = _S199 / norm_4;
    float y_15 = _S200 / norm_4;
    float z_12 = _S201 / norm_4;
    float3  _S202 = _S198 + make_float3 (0.48860251903533936f) * (make_float3 (- y_15) * (*sh_coeffs_4)[int(1)] + make_float3 (z_12) * (*sh_coeffs_4)[int(2)] - make_float3 (x_29) * (*sh_coeffs_4)[int(3)]);
    *rgb_4 = _S202;
    float z2_12 = z_12 * z_12;
    float fTmp0B_4 = -1.09254848957061768f * z_12;
    float fC1_4 = x_29 * x_29 - y_15 * y_15;
    float fS1_4 = 2.0f * x_29 * y_15;
    float3  _S203 = _S202 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_15) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_29) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
    *rgb_4 = _S203;
    float fTmp0C_4 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_12;
    *rgb_4 = max_0(_S203 + (make_float3 (-0.59004360437393188f * (x_29 * fS1_4 + y_15 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_15) * (*sh_coeffs_4)[int(11)] + make_float3 (z_12 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_29) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_29 * fC1_4 - y_15 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_9, float fy_9, float cx_9, float cy_9, float4  radial_coeffs_8, float2  tangential_coeffs_8, float2  thin_prism_coeffs_8, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_9, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_7, mean_5) + t_6;
    float3  _S204 = exp_0(scale_7);
    float x_30 = quat_8.y;
    float inv_norm_8 = (F32_rsqrt((x_30 * x_30 + quat_8.z * quat_8.z + quat_8.w * quat_8.w + quat_8.x * quat_8.x)));
    float x_31 = quat_8.y * inv_norm_8;
    float y_16 = quat_8.z * inv_norm_8;
    float z_13 = quat_8.w * inv_norm_8;
    float w_8 = quat_8.x * inv_norm_8;
    float x2_8 = x_31 * x_31;
    float y2_8 = y_16 * y_16;
    float z2_13 = z_13 * z_13;
    float xy_8 = x_31 * y_16;
    float xz_8 = x_31 * z_13;
    float yz_8 = y_16 * z_13;
    float wx_8 = w_8 * x_31;
    float wy_8 = w_8 * y_16;
    float wz_8 = w_8 * z_13;
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_13), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_13), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S204.x, 0.0f, 0.0f, 0.0f, _S204.y, 0.0f, 0.0f, 0.0f, _S204.z));
    Matrix<float, 3, 3>  _S205 = transpose_0(R_7);
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_7, mul_4(mul_4(R_7, mul_4(M_7, transpose_0(M_7))), _S205)), transpose_1(J_7));
    *mean2d_9 = make_float2 (fx_9 * mean_c_5.x + cx_9, fy_9 * mean_c_5.y + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S206 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S206;
    float _S207 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S207;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S206 * _S207 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S208 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_6, - covar2d_5.rows[int(0)].y * invdet_6, - covar2d_5.rows[int(1)].x * invdet_6, covar2d_5.rows[int(0)].x * invdet_6);
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
    *conic_5 = make_float3 (_S208.rows[int(0)].x, _S208.rows[int(0)].y, _S208.rows[int(1)].y);
    float3  _S209 = mean_5 - - mul_0(_S205, t_6);
    float3  _S210 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S210;
    float _S211 = _S209.x;
    float _S212 = _S209.y;
    float _S213 = _S209.z;
    float norm_5 = (F32_sqrt((_S211 * _S211 + _S212 * _S212 + _S213 * _S213)));
    float x_32 = _S211 / norm_5;
    float y_17 = _S212 / norm_5;
    float z_14 = _S213 / norm_5;
    float3  _S214 = _S210 + make_float3 (0.48860251903533936f) * (make_float3 (- y_17) * (*sh_coeffs_5)[int(1)] + make_float3 (z_14) * (*sh_coeffs_5)[int(2)] - make_float3 (x_32) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S214;
    float z2_14 = z_14 * z_14;
    float fTmp0B_5 = -1.09254848957061768f * z_14;
    float fC1_5 = x_32 * x_32 - y_17 * y_17;
    float fS1_5 = 2.0f * x_32 * y_17;
    float3  _S215 = _S214 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_17) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_14 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_32) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S215;
    float fTmp0C_5 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_14;
    *rgb_5 = max_0(_S215 + (make_float3 (-0.59004360437393188f * (x_32 * fS1_5 + y_17 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_17) * (*sh_coeffs_5)[int(11)] + make_float3 (z_14 * (1.86588168144226074f * z2_14 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_32) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_32 * fC1_5 - y_17 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S216, float3  _S217)
{
    return mul_0(_S216, _S217);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S218)
{
    return exp_0(_S218);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S219)
{
    return (F32_rsqrt((_S219)));
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S220, Matrix<float, 3, 3>  _S221)
{
    return mul_4(_S220, _S221);
}

inline __device__ float s_primal_ctx_max_0(float _S222, float _S223)
{
    return (F32_max((_S222), (_S223)));
}

inline __device__ float s_primal_ctx_min_0(float _S224, float _S225)
{
    return (F32_min((_S224), (_S225)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_2(Matrix<float, 2, 3>  _S226, Matrix<float, 3, 3>  _S227)
{
    return mul_5(_S226, _S227);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S228, Matrix<float, 3, 2>  _S229)
{
    return mul_6(_S228, _S229);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S230)
{
    return (F32_sqrt((_S230)));
}

inline __device__ float s_primal_ctx_exp_1(float _S231)
{
    return (F32_exp((_S231)));
}

inline __device__ float s_primal_ctx_log_0(float _S232)
{
    return (F32_log((_S232)));
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S233, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S234, float3  _S235)
{
    _d_max_vector_0(_S233, _S234, _S235);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S236, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S237, float3  _S238)
{
    _d_mul_0(_S236, _S237, _S238);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S239, DiffPair_float_0 * _S240, float _S241)
{
    _d_min_0(_S239, _S240, _S241);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S242, float _S243)
{
    _d_log_0(_S242, _S243);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S244, float _S245)
{
    _d_exp_0(_S244, _S245);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S246, DiffPair_float_0 * _S247, float _S248)
{
    _d_max_0(_S246, _S247, _S248);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S249, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S250, Matrix<float, 2, 2>  _S251)
{
    mul_3(_S249, _S250, _S251);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S252, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S253, Matrix<float, 2, 3>  _S254)
{
    mul_2(_S252, _S253, _S254);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S255, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S256, Matrix<float, 3, 3>  _S257)
{
    mul_1(_S255, _S256, _S257);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S258, float _S259)
{
    _d_rsqrt_0(_S258, _S259);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S260, float3  _S261)
{
    _d_exp_vector_0(_S260, _S261);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_10, float fy_10, float cx_10, float cy_10, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, uint image_width_6, uint image_height_6, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_6 = s_primal_ctx_mul_0(R_8, mean_6) + t_7;
    float3  _S262 = s_primal_ctx_exp_0(scale_8);
    float _S263 = quat_9.y;
    float _S264 = _S263 * _S263 + quat_9.z * quat_9.z + quat_9.w * quat_9.w + quat_9.x * quat_9.x;
    float _S265 = s_primal_ctx_rsqrt_0(_S264);
    float x_33 = quat_9.y * _S265;
    float y_18 = quat_9.z * _S265;
    float z_15 = quat_9.w * _S265;
    float w_9 = quat_9.x * _S265;
    float x2_9 = x_33 * x_33;
    float y2_9 = y_18 * y_18;
    float z2_15 = z_15 * z_15;
    float xy_9 = x_33 * y_18;
    float xz_9 = x_33 * z_15;
    float yz_9 = y_18 * z_15;
    float wx_9 = w_9 * x_33;
    float wy_9 = w_9 * y_18;
    float wz_9 = w_9 * z_15;
    Matrix<float, 3, 3>  _S266 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S262.x, 0.0f, 0.0f, 0.0f, _S262.y, 0.0f, 0.0f, 0.0f, _S262.z);
    Matrix<float, 3, 3>  _S267 = s_primal_ctx_mul_1(_S266, S_0);
    Matrix<float, 3, 3>  _S268 = transpose_0(_S267);
    Matrix<float, 3, 3>  _S269 = s_primal_ctx_mul_1(_S267, _S268);
    Matrix<float, 3, 3>  _S270 = s_primal_ctx_mul_1(R_8, _S269);
    Matrix<float, 3, 3>  _S271 = transpose_0(R_8);
    Matrix<float, 3, 3>  _S272 = s_primal_ctx_mul_1(_S270, _S271);
    float _S273 = float(image_width_6);
    float _S274 = float(image_height_6);
    float _S275 = 0.30000001192092896f * (0.5f * _S273 / fx_10);
    float lim_x_pos_0 = (_S273 - cx_10) / fx_10 + _S275;
    float _S276 = 0.30000001192092896f * (0.5f * _S274 / fy_10);
    float lim_y_pos_0 = (_S274 - cy_10) / fy_10 + _S276;
    float rz_4 = 1.0f / mean_c_6.z;
    float _S277 = mean_c_6.z * mean_c_6.z;
    float rz2_4 = rz_4 * rz_4;
    float _S278 = - (cx_10 / fx_10 + _S275);
    float _S279 = mean_c_6.x * rz_4;
    float _S280 = s_primal_ctx_max_0(_S278, _S279);
    float _S281 = s_primal_ctx_min_0(lim_x_pos_0, _S280);
    float _S282 = - (cy_10 / fy_10 + _S276);
    float _S283 = mean_c_6.y * rz_4;
    float _S284 = s_primal_ctx_max_0(_S282, _S283);
    float _S285 = s_primal_ctx_min_0(lim_y_pos_0, _S284);
    float _S286 = - fx_10;
    float _S287 = _S286 * (mean_c_6.z * _S281);
    float _S288 = - fy_10;
    float _S289 = _S288 * (mean_c_6.z * _S285);
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_10 * rz_4, 0.0f, _S287 * rz2_4, 0.0f, fy_10 * rz_4, _S289 * rz2_4);
    Matrix<float, 2, 3>  _S290 = s_primal_ctx_mul_2(J_8, _S272);
    Matrix<float, 3, 2>  _S291 = transpose_1(J_8);
    Matrix<float, 2, 2>  _S292 = s_primal_ctx_mul_3(_S290, _S291);
    float _S293 = fx_10 * mean_c_6.x;
    float _S294 = fy_10 * mean_c_6.y;
    float _S295 = _S292.rows[int(0)].y * _S292.rows[int(1)].x;
    float det_orig_7 = _S292.rows[int(0)].x * _S292.rows[int(1)].y - _S295;
    float _S296 = _S292.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S297 = _S292;
    *&(((&_S297)->rows + (int(0)))->x) = _S296;
    float _S298 = _S292.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S297)->rows + (int(1)))->y) = _S298;
    Matrix<float, 2, 2>  _S299 = _S297;
    Matrix<float, 2, 2>  _S300 = _S297;
    float det_blur_4 = _S296 * _S298 - _S295;
    float _S301 = det_orig_7 / det_blur_4;
    float _S302 = det_blur_4 * det_blur_4;
    float _S303 = s_primal_ctx_max_0(0.0f, _S301);
    float _S304 = s_primal_ctx_sqrt_0(_S303);
    float invdet_7 = 1.0f / det_blur_4;
    float _S305 = - _S292.rows[int(0)].y;
    float _S306 = - _S292.rows[int(1)].x;
    float _S307 = - in_opacity_6;
    float _S308 = 1.0f + s_primal_ctx_exp_1(_S307);
    float _S309 = 1.0f / _S308;
    float _S310 = _S308 * _S308;
    float _S311;
    if(antialiased_6)
    {
        _S311 = _S309 * _S304;
    }
    else
    {
        _S311 = _S309;
    }
    float _S312 = _S311 / 0.00392156885936856f;
    float _S313 = 2.0f * s_primal_ctx_log_0(_S312);
    float _S314 = s_primal_ctx_sqrt_0(_S313);
    float _S315 = _S299.rows[int(0)].x;
    float _S316 = _S300.rows[int(1)].y;
    float3  _S317 = mean_6 - - s_primal_ctx_mul_0(_S271, t_7);
    float _S318 = _S317.x;
    float _S319 = _S317.y;
    float _S320 = _S317.z;
    float _S321 = _S318 * _S318 + _S319 * _S319 + _S320 * _S320;
    float _S322 = s_primal_ctx_sqrt_0(_S321);
    float x_34 = _S318 / _S322;
    float3  _S323 = make_float3 (x_34);
    float _S324 = _S322 * _S322;
    float y_19 = _S319 / _S322;
    float z_16 = _S320 / _S322;
    float3  _S325 = make_float3 (z_16);
    float _S326 = - y_19;
    float3  _S327 = make_float3 (_S326);
    float z2_16 = z_16 * z_16;
    float fTmp0B_6 = -1.09254848957061768f * z_16;
    float fC1_6 = x_34 * x_34 - y_19 * y_19;
    float _S328 = 2.0f * x_34;
    float fS1_6 = _S328 * y_19;
    float pSH6_0 = 0.94617468118667603f * z2_16 - 0.31539157032966614f;
    float3  _S329 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_6 * x_34;
    float3  _S330 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_6 * y_19;
    float3  _S331 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_6;
    float3  _S332 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_6;
    float3  _S333 = make_float3 (pSH4_0);
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_16;
    float _S334 = 1.86588168144226074f * z2_16 - 1.11952900886535645f;
    float pSH12_0 = z_16 * _S334;
    float3  _S335 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_6 * x_34;
    float3  _S336 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_6 * y_19;
    float3  _S337 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_6 * fC1_6;
    float3  _S338 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_6 * fS1_6;
    float3  _S339 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_34 * fC1_6 - y_19 * fS1_6);
    float3  _S340 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_34 * fS1_6 + y_19 * fC1_6);
    float3  _S341 = make_float3 (pSH9_0);
    float3  _S342 = make_float3 (0.0f);
    float3  _S343 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S344;
    (&_S344)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S326) * (*sh_coeffs_6)[int(1)] + make_float3 (z_16) * (*sh_coeffs_6)[int(2)] - make_float3 (x_34) * (*sh_coeffs_6)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_6)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_6)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_6)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_6)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_6)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_6)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_6)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_6)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_6)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_6)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_6)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f);
    (&_S344)->differential_0 = _S343;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S345;
    (&_S345)->primal_0 = _S342;
    (&_S345)->differential_0 = _S343;
    s_bwd_prop_max_0(&_S344, &_S345, v_rgb_0);
    float3  _S346 = _S340 * _S344.differential_0;
    float3  _S347 = (*sh_coeffs_6)[int(15)] * _S344.differential_0;
    float3  _S348 = _S338 * _S344.differential_0;
    float3  _S349 = (*sh_coeffs_6)[int(14)] * _S344.differential_0;
    float3  _S350 = _S336 * _S344.differential_0;
    float3  _S351 = (*sh_coeffs_6)[int(13)] * _S344.differential_0;
    float3  _S352 = _S335 * _S344.differential_0;
    float3  _S353 = (*sh_coeffs_6)[int(12)] * _S344.differential_0;
    float3  _S354 = _S337 * _S344.differential_0;
    float3  _S355 = (*sh_coeffs_6)[int(11)] * _S344.differential_0;
    float3  _S356 = _S339 * _S344.differential_0;
    float3  _S357 = (*sh_coeffs_6)[int(10)] * _S344.differential_0;
    float3  _S358 = _S341 * _S344.differential_0;
    float3  _S359 = (*sh_coeffs_6)[int(9)] * _S344.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S359.x + _S359.y + _S359.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S347.x + _S347.y + _S347.z);
    float _S360 = _S357.x + _S357.y + _S357.z;
    float _S361 = _S349.x + _S349.y + _S349.z;
    float _S362 = _S355.x + _S355.y + _S355.z;
    float _S363 = _S351.x + _S351.y + _S351.z;
    float _S364 = _S353.x + _S353.y + _S353.z;
    float _S365 = - s_diff_fC2_T_0;
    float3  _S366 = _S332 * _S344.differential_0;
    float3  _S367 = (*sh_coeffs_6)[int(8)] * _S344.differential_0;
    float3  _S368 = _S330 * _S344.differential_0;
    float3  _S369 = (*sh_coeffs_6)[int(7)] * _S344.differential_0;
    float3  _S370 = _S329 * _S344.differential_0;
    float3  _S371 = (*sh_coeffs_6)[int(6)] * _S344.differential_0;
    float3  _S372 = _S331 * _S344.differential_0;
    float3  _S373 = (*sh_coeffs_6)[int(5)] * _S344.differential_0;
    float3  _S374 = _S333 * _S344.differential_0;
    float3  _S375 = (*sh_coeffs_6)[int(4)] * _S344.differential_0;
    float _S376 = _S373.x + _S373.y + _S373.z;
    float _S377 = _S369.x + _S369.y + _S369.z;
    float _S378 = fTmp1B_6 * _S360 + x_34 * s_diff_fS2_T_0 + y_19 * _S365 + 0.54627424478530884f * (_S375.x + _S375.y + _S375.z);
    float _S379 = fTmp1B_6 * _S361 + y_19 * s_diff_fS2_T_0 + x_34 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S367.x + _S367.y + _S367.z);
    float _S380 = y_19 * - _S379;
    float _S381 = x_34 * _S379;
    float _S382 = z_16 * (1.86588168144226074f * (z_16 * _S364) + -2.28522896766662598f * (y_19 * _S362 + x_34 * _S363) + 0.94617468118667603f * (_S371.x + _S371.y + _S371.z));
    float3  _S383 = make_float3 (0.48860251903533936f) * _S344.differential_0;
    float3  _S384 = - _S383;
    float3  _S385 = _S323 * _S384;
    float3  _S386 = (*sh_coeffs_6)[int(3)] * _S384;
    float3  _S387 = _S325 * _S383;
    float3  _S388 = (*sh_coeffs_6)[int(2)] * _S383;
    float3  _S389 = _S327 * _S383;
    float3  _S390 = (*sh_coeffs_6)[int(1)] * _S383;
    float _S391 = (_S334 * _S364 + 1.44530570507049561f * (fS1_6 * _S360 + fC1_6 * _S361) + -1.09254848957061768f * (y_19 * _S376 + x_34 * _S377) + _S382 + _S382 + _S388.x + _S388.y + _S388.z) / _S324;
    float _S392 = _S322 * _S391;
    float _S393 = (fTmp0C_6 * _S362 + fC1_6 * s_diff_fS2_T_0 + fS1_6 * _S365 + fTmp0B_6 * _S376 + _S328 * _S378 + _S380 + _S380 + - (_S390.x + _S390.y + _S390.z)) / _S324;
    float _S394 = _S322 * _S393;
    float _S395 = (fTmp0C_6 * _S363 + fS1_6 * s_diff_fS2_T_0 + fC1_6 * s_diff_fC2_T_0 + fTmp0B_6 * _S377 + 2.0f * (y_19 * _S378) + _S381 + _S381 + _S386.x + _S386.y + _S386.z) / _S324;
    float _S396 = _S322 * _S395;
    float _S397 = _S320 * - _S391 + _S319 * - _S393 + _S318 * - _S395;
    DiffPair_float_0 _S398;
    (&_S398)->primal_0 = _S321;
    (&_S398)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S398, _S397);
    float _S399 = _S320 * _S398.differential_0;
    float _S400 = _S319 * _S398.differential_0;
    float _S401 = _S318 * _S398.differential_0;
    float3  _S402 = make_float3 (0.282094806432724f) * _S344.differential_0;
    float3  _S403 = make_float3 (_S396 + _S401 + _S401, _S394 + _S400 + _S400, _S392 + _S399 + _S399);
    float3  _S404 = - - _S403;
    Matrix<float, 3, 3>  _S405 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S406;
    (&_S406)->primal_0 = _S271;
    (&_S406)->differential_0 = _S405;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S407;
    (&_S407)->primal_0 = t_7;
    (&_S407)->differential_0 = _S343;
    s_bwd_prop_mul_0(&_S406, &_S407, _S404);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S408 = _S406;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S409 = _S407;
    float2  _S410 = make_float2 (0.0f);
    float2  _S411 = _S410;
    *&((&_S411)->y) = v_conic_0.z;
    float2  _S412 = _S410;
    *&((&_S412)->y) = v_conic_0.y;
    *&((&_S412)->x) = v_conic_0.x;
    DiffPair_float_0 _S413;
    (&_S413)->primal_0 = _S316;
    (&_S413)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S413, 0.0f);
    DiffPair_float_0 _S414;
    (&_S414)->primal_0 = _S315;
    (&_S414)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S414, 0.0f);
    DiffPair_float_0 _S415;
    (&_S415)->primal_0 = 3.32999992370605469f;
    (&_S415)->differential_0 = 0.0f;
    DiffPair_float_0 _S416;
    (&_S416)->primal_0 = _S314;
    (&_S416)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S415, &_S416, 0.0f);
    DiffPair_float_0 _S417;
    (&_S417)->primal_0 = _S313;
    (&_S417)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S417, _S416.differential_0);
    float _S418 = 2.0f * _S417.differential_0;
    DiffPair_float_0 _S419;
    (&_S419)->primal_0 = _S312;
    (&_S419)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S419, _S418);
    float _S420 = v_opacity_0 + 254.9999847412109375f * _S419.differential_0;
    Matrix<float, 2, 2>  _S421 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S422 = _S421;
    _S422[int(1)] = _S411;
    _S422[int(0)] = _S412;
    Matrix<float, 2, 2>  _S423 = _S422;
    FixedArray<float3 , 16>  _S424;
    _S424[int(0)] = _S343;
    _S424[int(1)] = _S343;
    _S424[int(2)] = _S343;
    _S424[int(3)] = _S343;
    _S424[int(4)] = _S343;
    _S424[int(5)] = _S343;
    _S424[int(6)] = _S343;
    _S424[int(7)] = _S343;
    _S424[int(8)] = _S343;
    _S424[int(9)] = _S343;
    _S424[int(10)] = _S343;
    _S424[int(11)] = _S343;
    _S424[int(12)] = _S343;
    _S424[int(13)] = _S343;
    _S424[int(14)] = _S343;
    _S424[int(15)] = _S343;
    _S424[int(7)] = _S368;
    _S424[int(0)] = _S402;
    _S424[int(1)] = _S389;
    _S424[int(2)] = _S387;
    _S424[int(3)] = _S385;
    _S424[int(4)] = _S374;
    _S424[int(5)] = _S372;
    _S424[int(6)] = _S370;
    _S424[int(15)] = _S346;
    _S424[int(8)] = _S366;
    _S424[int(9)] = _S358;
    _S424[int(10)] = _S356;
    _S424[int(11)] = _S354;
    _S424[int(12)] = _S352;
    _S424[int(13)] = _S350;
    _S424[int(14)] = _S348;
    float3  _S425 = _S424[int(0)];
    float3  _S426 = _S424[int(1)];
    float3  _S427 = _S424[int(2)];
    float3  _S428 = _S424[int(3)];
    float3  _S429 = _S424[int(4)];
    float3  _S430 = _S424[int(5)];
    float3  _S431 = _S424[int(6)];
    float3  _S432 = _S424[int(7)];
    float3  _S433 = _S424[int(8)];
    float3  _S434 = _S424[int(9)];
    float3  _S435 = _S424[int(10)];
    float3  _S436 = _S424[int(11)];
    float3  _S437 = _S424[int(12)];
    float3  _S438 = _S424[int(13)];
    float3  _S439 = _S424[int(14)];
    float3  _S440 = _S424[int(15)];
    float3  _S441 = make_float3 (0.0f, 0.0f, v_depth_0);
    float2  _S442 = make_float2 (0.0f, _S413.differential_0);
    float2  _S443 = make_float2 (_S414.differential_0, 0.0f);
    float _S444;
    if(antialiased_6)
    {
        float _S445 = _S309 * _S420;
        _S311 = _S304 * _S420;
        _S444 = _S445;
    }
    else
    {
        _S311 = _S420;
        _S444 = 0.0f;
    }
    float _S446 = - (_S311 / _S310);
    DiffPair_float_0 _S447;
    (&_S447)->primal_0 = _S307;
    (&_S447)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S447, _S446);
    float _S448 = - _S447.differential_0;
    float _S449 = invdet_7 * _S423.rows[int(1)].y;
    float _S450 = - (invdet_7 * _S423.rows[int(1)].x);
    float _S451 = - (invdet_7 * _S423.rows[int(0)].y);
    float _S452 = invdet_7 * _S423.rows[int(0)].x;
    float _S453 = - ((_S296 * _S423.rows[int(1)].y + _S306 * _S423.rows[int(1)].x + _S305 * _S423.rows[int(0)].y + _S298 * _S423.rows[int(0)].x) / _S302);
    DiffPair_float_0 _S454;
    (&_S454)->primal_0 = _S303;
    (&_S454)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S454, _S444);
    DiffPair_float_0 _S455;
    (&_S455)->primal_0 = 0.0f;
    (&_S455)->differential_0 = 0.0f;
    DiffPair_float_0 _S456;
    (&_S456)->primal_0 = _S301;
    (&_S456)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S455, &_S456, _S454.differential_0);
    float _S457 = _S456.differential_0 / _S302;
    float s_diff_det_orig_T_0 = det_blur_4 * _S457;
    float _S458 = _S453 + det_orig_7 * - _S457;
    float _S459 = - _S458;
    float _S460 = _S296 * _S458;
    float _S461 = _S298 * _S458;
    Matrix<float, 2, 2>  _S462 = _S421;
    _S462[int(1)] = _S442;
    _S462[int(0)] = _S443;
    _S297 = _S462;
    *&(((&_S297)->rows + (int(1)))->y) = 0.0f;
    float _S463 = _S452 + _S460 + _S462.rows[int(1)].y;
    *&(((&_S297)->rows + (int(0)))->x) = 0.0f;
    float _S464 = _S449 + _S461 + _S462.rows[int(0)].x;
    float _S465 = _S459 + - s_diff_det_orig_T_0;
    float _S466 = _S450 + _S292.rows[int(0)].y * _S465;
    float _S467 = _S451 + _S292.rows[int(1)].x * _S465;
    float _S468 = _S292.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S469 = _S463 + _S292.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S470 = _S410;
    *&((&_S470)->x) = _S466;
    *&((&_S470)->y) = _S469;
    float _S471 = _S464 + _S468;
    float2  _S472 = _S410;
    *&((&_S472)->y) = _S467;
    *&((&_S472)->x) = _S471;
    float _S473 = _S294 * v_mean2d_0.y;
    float _S474 = fy_10 * (rz_4 * v_mean2d_0.y);
    float _S475 = _S293 * v_mean2d_0.x;
    float _S476 = fx_10 * (rz_4 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S477 = _S421;
    _S477[int(1)] = _S470;
    _S477[int(0)] = _S472;
    Matrix<float, 2, 2>  _S478 = _S297 + _S477;
    Matrix<float, 2, 3>  _S479 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S480;
    (&_S480)->primal_0 = _S290;
    (&_S480)->differential_0 = _S479;
    Matrix<float, 3, 2>  _S481 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S482;
    (&_S482)->primal_0 = _S291;
    (&_S482)->differential_0 = _S481;
    s_bwd_prop_mul_1(&_S480, &_S482, _S478);
    Matrix<float, 2, 3>  _S483 = transpose_2(_S482.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S484;
    (&_S484)->primal_0 = J_8;
    (&_S484)->differential_0 = _S479;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S485;
    (&_S485)->primal_0 = _S272;
    (&_S485)->differential_0 = _S405;
    s_bwd_prop_mul_2(&_S484, &_S485, _S480.differential_0);
    Matrix<float, 2, 3>  _S486 = _S483 + _S484.differential_0;
    float _S487 = _S289 * _S486.rows[int(1)].z;
    float s_diff_ty_T_0 = _S288 * (rz2_4 * _S486.rows[int(1)].z);
    float _S488 = fy_10 * _S486.rows[int(1)].y;
    float _S489 = _S287 * _S486.rows[int(0)].z;
    float s_diff_tx_T_0 = _S286 * (rz2_4 * _S486.rows[int(0)].z);
    float _S490 = fx_10 * _S486.rows[int(0)].x;
    float _S491 = mean_c_6.z * s_diff_ty_T_0;
    float _S492 = _S285 * s_diff_ty_T_0;
    DiffPair_float_0 _S493;
    (&_S493)->primal_0 = lim_y_pos_0;
    (&_S493)->differential_0 = 0.0f;
    DiffPair_float_0 _S494;
    (&_S494)->primal_0 = _S284;
    (&_S494)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S493, &_S494, _S491);
    DiffPair_float_0 _S495;
    (&_S495)->primal_0 = _S282;
    (&_S495)->differential_0 = 0.0f;
    DiffPair_float_0 _S496;
    (&_S496)->primal_0 = _S283;
    (&_S496)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S495, &_S496, _S494.differential_0);
    float _S497 = mean_c_6.y * _S496.differential_0;
    float _S498 = rz_4 * _S496.differential_0;
    float _S499 = mean_c_6.z * s_diff_tx_T_0;
    float _S500 = _S281 * s_diff_tx_T_0;
    DiffPair_float_0 _S501;
    (&_S501)->primal_0 = lim_x_pos_0;
    (&_S501)->differential_0 = 0.0f;
    DiffPair_float_0 _S502;
    (&_S502)->primal_0 = _S280;
    (&_S502)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S501, &_S502, _S499);
    DiffPair_float_0 _S503;
    (&_S503)->primal_0 = _S278;
    (&_S503)->differential_0 = 0.0f;
    DiffPair_float_0 _S504;
    (&_S504)->primal_0 = _S279;
    (&_S504)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S503, &_S504, _S502.differential_0);
    float _S505 = rz_4 * (_S487 + _S489);
    float _S506 = _S492 + _S500 + - ((_S473 + _S475 + _S488 + _S490 + _S497 + mean_c_6.x * _S504.differential_0 + _S505 + _S505) / _S277);
    float _S507 = _S474 + _S498;
    float _S508 = _S476 + rz_4 * _S504.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S509;
    (&_S509)->primal_0 = _S270;
    (&_S509)->differential_0 = _S405;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S510;
    (&_S510)->primal_0 = _S271;
    (&_S510)->differential_0 = _S405;
    s_bwd_prop_mul_3(&_S509, &_S510, _S485.differential_0);
    Matrix<float, 3, 3>  _S511 = transpose_0(_S510.differential_0 + _S408.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S512;
    (&_S512)->primal_0 = R_8;
    (&_S512)->differential_0 = _S405;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S513;
    (&_S513)->primal_0 = _S269;
    (&_S513)->differential_0 = _S405;
    s_bwd_prop_mul_3(&_S512, &_S513, _S509.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S514;
    (&_S514)->primal_0 = _S267;
    (&_S514)->differential_0 = _S405;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S515;
    (&_S515)->primal_0 = _S268;
    (&_S515)->differential_0 = _S405;
    s_bwd_prop_mul_3(&_S514, &_S515, _S513.differential_0);
    Matrix<float, 3, 3>  _S516 = _S514.differential_0 + transpose_0(_S515.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S517;
    (&_S517)->primal_0 = _S266;
    (&_S517)->differential_0 = _S405;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S518;
    (&_S518)->primal_0 = S_0;
    (&_S518)->differential_0 = _S405;
    s_bwd_prop_mul_3(&_S517, &_S518, _S516);
    Matrix<float, 3, 3>  _S519 = transpose_0(_S517.differential_0);
    float _S520 = 2.0f * - _S519.rows[int(2)].z;
    float _S521 = 2.0f * _S519.rows[int(2)].y;
    float _S522 = 2.0f * _S519.rows[int(2)].x;
    float _S523 = 2.0f * _S519.rows[int(1)].z;
    float _S524 = 2.0f * - _S519.rows[int(1)].y;
    float _S525 = 2.0f * _S519.rows[int(1)].x;
    float _S526 = 2.0f * _S519.rows[int(0)].z;
    float _S527 = 2.0f * _S519.rows[int(0)].y;
    float _S528 = 2.0f * - _S519.rows[int(0)].x;
    float _S529 = - _S525 + _S527;
    float _S530 = _S522 + - _S526;
    float _S531 = - _S521 + _S523;
    float _S532 = _S521 + _S523;
    float _S533 = _S522 + _S526;
    float _S534 = _S525 + _S527;
    float _S535 = z_15 * (_S524 + _S528);
    float _S536 = y_18 * (_S520 + _S528);
    float _S537 = x_33 * (_S520 + _S524);
    float _S538 = z_15 * _S529 + y_18 * _S530 + x_33 * _S531;
    float _S539 = _S265 * _S538;
    float _S540 = w_9 * _S529 + y_18 * _S532 + x_33 * _S533 + _S535 + _S535;
    float _S541 = _S265 * _S540;
    float _S542 = w_9 * _S530 + z_15 * _S532 + x_33 * _S534 + _S536 + _S536;
    float _S543 = _S265 * _S542;
    float _S544 = w_9 * _S531 + z_15 * _S533 + y_18 * _S534 + _S537 + _S537;
    float _S545 = _S265 * _S544;
    float _S546 = quat_9.x * _S538 + quat_9.w * _S540 + quat_9.z * _S542 + quat_9.y * _S544;
    DiffPair_float_0 _S547;
    (&_S547)->primal_0 = _S264;
    (&_S547)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S547, _S546);
    float _S548 = quat_9.x * _S547.differential_0;
    float _S549 = quat_9.w * _S547.differential_0;
    float _S550 = quat_9.z * _S547.differential_0;
    float _S551 = quat_9.y * _S547.differential_0;
    float _S552 = _S541 + _S549 + _S549;
    float _S553 = _S543 + _S550 + _S550;
    float _S554 = _S545 + _S551 + _S551;
    float _S555 = _S539 + _S548 + _S548;
    float3  _S556 = _S343;
    *&((&_S556)->z) = _S518.differential_0.rows[int(2)].z;
    *&((&_S556)->y) = _S518.differential_0.rows[int(1)].y;
    *&((&_S556)->x) = _S518.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S557;
    (&_S557)->primal_0 = scale_8;
    (&_S557)->differential_0 = _S343;
    s_bwd_prop_exp_1(&_S557, _S556);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S558 = _S557;
    float3  _S559 = _S343;
    *&((&_S559)->z) = _S506;
    *&((&_S559)->y) = _S507;
    *&((&_S559)->x) = _S508;
    float3  _S560 = _S441 + _S559;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S561;
    (&_S561)->primal_0 = R_8;
    (&_S561)->differential_0 = _S405;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S562;
    (&_S562)->primal_0 = mean_6;
    (&_S562)->differential_0 = _S343;
    s_bwd_prop_mul_0(&_S561, &_S562, _S560);
    float3  _S563 = _S560 + _S409.differential_0;
    Matrix<float, 3, 3>  _S564 = _S511 + _S512.differential_0 + _S561.differential_0;
    float4  _S565 = make_float4 (0.0f);
    *&((&_S565)->w) = _S552;
    *&((&_S565)->z) = _S553;
    *&((&_S565)->y) = _S554;
    *&((&_S565)->x) = _S555;
    float4  _S566 = _S565;
    float3  _S567 = _S562.differential_0 + _S403;
    *v_mean_0 = _S567;
    *v_quat_0 = _S566;
    *v_scale_0 = _S558.differential_0;
    *v_in_opacity_0 = _S448;
    (*v_sh_coeffs_0)[int(0)] = _S425;
    (*v_sh_coeffs_0)[int(1)] = _S426;
    (*v_sh_coeffs_0)[int(2)] = _S427;
    (*v_sh_coeffs_0)[int(3)] = _S428;
    (*v_sh_coeffs_0)[int(4)] = _S429;
    (*v_sh_coeffs_0)[int(5)] = _S430;
    (*v_sh_coeffs_0)[int(6)] = _S431;
    (*v_sh_coeffs_0)[int(7)] = _S432;
    (*v_sh_coeffs_0)[int(8)] = _S433;
    (*v_sh_coeffs_0)[int(9)] = _S434;
    (*v_sh_coeffs_0)[int(10)] = _S435;
    (*v_sh_coeffs_0)[int(11)] = _S436;
    (*v_sh_coeffs_0)[int(12)] = _S437;
    (*v_sh_coeffs_0)[int(13)] = _S438;
    (*v_sh_coeffs_0)[int(14)] = _S439;
    (*v_sh_coeffs_0)[int(15)] = _S440;
    *v_R_0 = _S564;
    *v_t_0 = _S563;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S568;
    DiffPair_float_0 _S569;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S570;
    DiffPair_float_0 _S571;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S572;
    DiffPair_float_0 _S573;
    DiffPair_float_0 _S574;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S575;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S576;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S577;
};

inline __device__ CameraDistortion_0 s_primal_ctx_CameraDistortion_x24init_0(float4  dpradial_coeffs_0, float2  dptangential_coeffs_0, float2  dpthin_prism_coeffs_0)
{
    CameraDistortion_0 _S578 = { dpradial_coeffs_0, dptangential_coeffs_0, dpthin_prism_coeffs_0 };
    return _S578;
}

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S579, float _S580)
{
    return s_primal_ctx_atan2_0(_S579, _S580);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S581;
    DiffPair_float_0 _S582;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_0)
{
    DiffPair_float_0 _S583 = { 0.0f, 0.0f };
    _s_diff_ctx_0->_S581 = _S583;
    _s_diff_ctx_0->_S582 = _S583;
    (&_s_diff_ctx_0->_S581)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S581)->differential_0 = 0.0f;
    (&_s_diff_ctx_0->_S582)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S582)->differential_0 = 0.0f;
    DiffPair_float_0 _S584 = *dpdpy_0;
    _s_diff_ctx_0->_S581 = *dpdpy_0;
    DiffPair_float_0 _S585 = *dpdpx_0;
    _s_diff_ctx_0->_S582 = *dpdpx_0;
    float _S586 = _S585.primal_0 * _S585.primal_0 + _S584.primal_0 * _S584.primal_0;
    float _S587 = - _S584.primal_0 / _S586 * dpdOut_0;
    float _S588 = _S585.primal_0 / _S586 * dpdOut_0;
    dpdpy_0->primal_0 = _S584.primal_0;
    dpdpy_0->differential_0 = _S588;
    dpdpx_0->primal_0 = _S585.primal_0;
    dpdpx_0->differential_0 = _S587;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S589, DiffPair_float_0 * _S590, float _S591, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_float_0 _S592 = { 0.0f, 0.0f };
    _s_diff_ctx_1->_S568 = _S592;
    _s_diff_ctx_1->_S569 = _S592;
    (&_s_diff_ctx_1->_S568)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S568)->differential_0 = 0.0f;
    (&_s_diff_ctx_1->_S569)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S569)->differential_0 = 0.0f;
    DiffPair_float_0 _S593 = *_S589;
    _s_diff_ctx_1->_S568 = *_S589;
    DiffPair_float_0 _S594 = *_S590;
    _s_diff_ctx_1->_S569 = *_S590;
    DiffPair_float_0 _S595 = _S593;
    DiffPair_float_0 _S596 = _S594;
    s_bwd_prop_d_atan2_Intermediates_0 _S597;
    (&_S597)->_S581 = _S592;
    (&_S597)->_S582 = _S592;
    s_primal_ctx_d_atan2_0(&_S595, &_S596, _S591, &_S597);
    *_S589 = _S595;
    *_S590 = _S596;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S598;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S599;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S600;
    DiffPair_float_0 _S601;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S602;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S603;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S604 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S603 = _S604;
    (&_s_diff_ctx_2->_S603)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S603)->differential_0 = 0.0f;
    DiffPair_float_0 _S605 = *dpdpx_1;
    _s_diff_ctx_2->_S603 = *dpdpx_1;
    float _S606 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S605.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S605.primal_0;
    dpdpx_1->differential_0 = _S606;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S607, float _S608, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S609 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S599 = _S609;
    (&_s_diff_ctx_3->_S599)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S599)->differential_0 = 0.0f;
    DiffPair_float_0 _S610 = *_S607;
    _s_diff_ctx_3->_S599 = *_S607;
    DiffPair_float_0 _S611 = _S610;
    s_bwd_prop_d_sqrt_Intermediates_0 _S612;
    (&_S612)->_S603 = _S609;
    s_primal_ctx_d_sqrt_0(&_S611, _S608, &_S612);
    *_S607 = _S611;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_4)
{
    float2  _S613 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S614 = { _S613, _S613 };
    DiffPair_float_0 _S615 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S616 = { _S615 };
    _s_diff_ctx_4->_S600 = _S614;
    _s_diff_ctx_4->_S601 = _S615;
    _s_diff_ctx_4->_S602 = _S616;
    (&_s_diff_ctx_4->_S600)->primal_0 = _S613;
    (&_s_diff_ctx_4->_S600)->differential_0 = _S613;
    (&_s_diff_ctx_4->_S601)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S601)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S617 = *dpdpx_2;
    _s_diff_ctx_4->_S600 = *dpdpx_2;
    float _S618 = _S617.primal_0.x;
    float _S619 = _S617.primal_0.y;
    DiffPair_float_0 _S620;
    (&_S620)->primal_0 = _S618 * _S618 + _S619 * _S619;
    (&_S620)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S620, dp_s_dOut_0, &_s_diff_ctx_4->_S602);
    _s_diff_ctx_4->_S601 = _S620;
    float _S621 = _S617.primal_0.y * _S620.differential_0;
    float _S622 = _S621 + _S621;
    float _S623 = _S617.primal_0.x * _S620.differential_0;
    float _S624 = _S623 + _S623;
    float2  _S625 = _S613;
    *&((&_S625)->y) = _S622;
    *&((&_S625)->x) = _S624;
    dpdpx_2->primal_0 = _S617.primal_0;
    dpdpx_2->differential_0 = _S625;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S626, float _S627, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_5)
{
    float2  _S628 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S629 = { _S628, _S628 };
    _s_diff_ctx_5->_S598 = _S629;
    (&_s_diff_ctx_5->_S598)->primal_0 = _S628;
    (&_s_diff_ctx_5->_S598)->differential_0 = _S628;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S630 = *_S626;
    _s_diff_ctx_5->_S598 = *_S626;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S631 = _S630;
    DiffPair_float_0 _S632 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S633 = { _S632 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S634;
    (&_S634)->_S600 = _S629;
    (&_S634)->_S601 = _S632;
    (&_S634)->_S602 = _S633;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S631, _S627, &_S634);
    *_S626 = _S631;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, float4  dpradial_coeffs_1, float2  dptangential_coeffs_1, float2  dpthin_prism_coeffs_1, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_6)
{
    DiffPair_float_0 _S635 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S636 = { _S635, _S635 };
    _s_diff_ctx_6->_S570 = _S635;
    _s_diff_ctx_6->_S571 = _S635;
    _s_diff_ctx_6->_S572 = _S636;
    _s_diff_ctx_6->_S573 = _S635;
    _s_diff_ctx_6->_S574 = _S635;
    _s_diff_ctx_6->_S575 = _S636;
    (&_s_diff_ctx_6->_S570)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S570)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S571)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S571)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S573)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S573)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S574)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S574)->differential_0 = 0.0f;
    float2  _S637 = make_float2 (0.0f);
    CameraDistortion_0 _S638 = s_primal_ctx_CameraDistortion_x24init_0(dpradial_coeffs_1, dptangential_coeffs_1, dpthin_prism_coeffs_1);
    float2  _S639 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S640 = length_0(_S639);
    float _S641 = dpmean3d_0.z;
    float _S642 = s_primal_ctx_atan2_0(_S640, _S641);
    float k_1;
    if(_S642 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S642 * _S642 / 3.0f) / _S641;
    }
    else
    {
        k_1 = _S642 / _S640;
    }
    float2  _S643 = _S639 * make_float2 (k_1);
    float k1_1 = _S638.radial_coeffs_0.x;
    float k2_1 = _S638.radial_coeffs_0.y;
    float k3_1 = _S638.radial_coeffs_0.z;
    float k4_1 = _S638.radial_coeffs_0.w;
    float p1_1 = _S638.tangential_coeffs_0.x;
    float p2_1 = _S638.tangential_coeffs_0.y;
    float sx1_1 = _S638.thin_prism_coeffs_0.x;
    float sy1_1 = _S638.thin_prism_coeffs_0.y;
    float u_3 = _S643.x;
    float v_3 = _S643.y;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float _S644 = 2.0f * p1_1;
    float _S645 = 2.0f * p2_1;
    float2  _S646 = _S643 * make_float2 (1.0f + r2_3 * (k1_1 + r2_3 * (k2_1 + r2_3 * (k3_1 + r2_3 * k4_1)))) + make_float2 (_S644 * u_3 * v_3 + p2_1 * (r2_3 + 2.0f * u_3 * u_3) + sx1_1 * r2_3, _S645 * u_3 * v_3 + p1_1 * (r2_3 + 2.0f * v_3 * v_3) + sy1_1 * r2_3);
    float2  _S647 = make_float2 (dpfx_0 * _S646.x + dpcx_0, dpfy_0 * _S646.y + dpcy_0);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (0.0f);
    float _S648 = s_primal_ctx_s_primal_ctx_atan2_0(_S640, _S641);
    bool _S649 = _S648 < 0.00100000004749745f;
    float _S650;
    float _S651;
    float _S652;
    if(_S649)
    {
        float _S653 = 1.0f - _S648 * _S648 / 3.0f;
        float _S654 = _S641 * _S641;
        k_1 = _S653 / _S641;
        _S650 = 0.0f;
        _S651 = _S654;
        _S652 = _S653;
    }
    else
    {
        float _S655 = _S640 * _S640;
        k_1 = _S648 / _S640;
        _S650 = _S655;
        _S651 = 0.0f;
        _S652 = 0.0f;
    }
    float2  _S656 = make_float2 (k_1);
    float2  _S657 = _S639 * make_float2 (k_1);
    float u_4 = _S657.x;
    float v_4 = _S657.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float _S658 = k3_1 + r2_4 * k4_1;
    float _S659 = k2_1 + r2_4 * _S658;
    float _S660 = k1_1 + r2_4 * _S659;
    float2  _S661 = make_float2 (dpfx_0, 0.0f);
    float2  _S662 = _S657 * _S661;
    float _S663 = p2_1 * dpfx_0;
    float _S664 = _S662.x + _S662.y;
    float _S665 = r2_4 * _S664;
    float _S666 = r2_4 * _S665;
    float _S667 = sx1_1 * dpfx_0 + _S663 + _S660 * _S664 + _S659 * _S665 + _S658 * _S666 + k4_1 * (r2_4 * _S666);
    float _S668 = v_4 * _S667;
    float _S669 = u_4 * _S667;
    float2  _S670 = make_float2 (1.0f + r2_4 * _S660) * _S661 + make_float2 (2.0f * u_4 * _S663 + 2.0f * (u_4 * _S663) + _S644 * (v_4 * dpfx_0) + _S669 + _S669, _S644 * u_4 * dpfx_0 + _S668 + _S668);
    float2  _S671 = _S639 * _S670;
    float2  _S672 = _S656 * _S670;
    float _S673 = _S671.x + _S671.y;
    if(_S649)
    {
        float _S674 = _S673 / _S651;
        float _S675 = _S652 * - _S674;
        float _S676 = _S648 * (0.3333333432674408f * - (_S641 * _S674));
        k_1 = _S676 + _S676;
        _S650 = _S675;
        _S651 = 0.0f;
    }
    else
    {
        float _S677 = _S673 / _S650;
        float _S678 = _S648 * - _S677;
        k_1 = _S640 * _S677;
        _S650 = 0.0f;
        _S651 = _S678;
    }
    DiffPair_float_0 _S679;
    (&_S679)->primal_0 = _S640;
    (&_S679)->differential_0 = 0.0f;
    DiffPair_float_0 _S680;
    (&_S680)->primal_0 = _S641;
    (&_S680)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S679, &_S680, k_1, &_s_diff_ctx_6->_S572);
    _s_diff_ctx_6->_S570 = _S679;
    _s_diff_ctx_6->_S571 = _S680;
    float _S681 = _S680.differential_0 + _S650;
    float _S682 = _S679.differential_0 + _S651;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S683;
    (&_S683)->primal_0 = _S639;
    (&_S683)->differential_0 = _S637;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S684 = { _S637, _S637 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S685;
    (&_S685)->_S598 = _S684;
    s_primal_ctx_s_bwd_length_impl_0(&_S683, _S682, &_S685);
    float2  _S686 = _S683.differential_0 + _S672;
    float3  _S687 = make_float3 (_S686.x, _S686.y, _S681);
    Matrix<float, 2, 3>  _S688 = J_9;
    _S688[int(0)] = _S687;
    if(_S649)
    {
        float _S689 = 1.0f - _S648 * _S648 / 3.0f;
        float _S690 = _S641 * _S641;
        k_1 = _S689 / _S641;
        _S650 = 0.0f;
        _S651 = _S690;
        _S652 = _S689;
    }
    else
    {
        float _S691 = _S640 * _S640;
        k_1 = _S648 / _S640;
        _S650 = _S691;
        _S651 = 0.0f;
        _S652 = 0.0f;
    }
    float2  _S692 = make_float2 (k_1);
    float2  _S693 = _S639 * make_float2 (k_1);
    float u_5 = _S693.x;
    float v_5 = _S693.y;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float _S694 = k3_1 + r2_5 * k4_1;
    float _S695 = k2_1 + r2_5 * _S694;
    float _S696 = k1_1 + r2_5 * _S695;
    float2  _S697 = make_float2 (0.0f, dpfy_0);
    float2  _S698 = _S693 * _S697;
    float _S699 = p1_1 * dpfy_0;
    float _S700 = _S698.x + _S698.y;
    float _S701 = r2_5 * _S700;
    float _S702 = r2_5 * _S701;
    float _S703 = sy1_1 * dpfy_0 + _S699 + _S696 * _S700 + _S695 * _S701 + _S694 * _S702 + k4_1 * (r2_5 * _S702);
    float _S704 = v_5 * _S703;
    float _S705 = u_5 * _S703;
    float2  _S706 = make_float2 (1.0f + r2_5 * _S696) * _S697 + make_float2 (_S645 * (v_5 * dpfy_0) + _S705 + _S705, 2.0f * v_5 * _S699 + 2.0f * (v_5 * _S699) + _S645 * u_5 * dpfy_0 + _S704 + _S704);
    float2  _S707 = _S639 * _S706;
    float2  _S708 = _S692 * _S706;
    float _S709 = _S707.x + _S707.y;
    if(_S649)
    {
        float _S710 = _S709 / _S651;
        float _S711 = _S652 * - _S710;
        float _S712 = _S648 * (0.3333333432674408f * - (_S641 * _S710));
        k_1 = _S712 + _S712;
        _S650 = _S711;
        _S651 = 0.0f;
    }
    else
    {
        float _S713 = _S709 / _S650;
        float _S714 = _S648 * - _S713;
        k_1 = _S640 * _S713;
        _S650 = 0.0f;
        _S651 = _S714;
    }
    DiffPair_float_0 _S715;
    (&_S715)->primal_0 = _S640;
    (&_S715)->differential_0 = 0.0f;
    DiffPair_float_0 _S716;
    (&_S716)->primal_0 = _S641;
    (&_S716)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S715, &_S716, k_1, &_s_diff_ctx_6->_S575);
    _s_diff_ctx_6->_S573 = _S715;
    _s_diff_ctx_6->_S574 = _S716;
    float _S717 = _S716.differential_0 + _S650;
    float _S718 = _S715.differential_0 + _S651;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S719;
    (&_S719)->primal_0 = _S639;
    (&_S719)->differential_0 = _S637;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S720;
    (&_S720)->_S598 = _S684;
    s_primal_ctx_s_bwd_length_impl_0(&_S719, _S718, &_S720);
    float2  _S721 = _S719.differential_0 + _S708;
    float3  _S722 = make_float3 (_S721.x, _S721.y, _S717);
    _S688[int(1)] = _S722;
    *dpcov2d_0 = s_primal_ctx_mul_3(s_primal_ctx_mul_2(_S688, dpcov3d_0), transpose_1(_S688));
    *dpmean2d_0 = _S647;
    return;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
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

inline __device__ void s_bwd_prop_d_sqrt_0(DiffPair_1 * dpdpx_3, DiffPair_float_0 * dpdOut_2, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_7)
{
    DiffPair_1 _S723 = *dpdpx_3;
    float _S724 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_7->_S603)->primal_0);
    float _S725 = s_primal_ctx_sqrt_0(_S724);
    float _S726 = 0.5f / _S725 * (*dpdpx_3).differential_0.differential_0;
    float _S727 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S725 * _S725));
    DiffPair_float_0 _S728;
    (&_S728)->primal_0 = _S724;
    (&_S728)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S728, _S727);
    DiffPair_float_0 _S729;
    (&_S729)->primal_0 = 1.00000001168609742e-07f;
    (&_S729)->differential_0 = 0.0f;
    DiffPair_float_0 _S730;
    (&_S730)->primal_0 = (&_s_diff_ctx_7->_S603)->primal_0;
    (&_S730)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S729, &_S730, _S728.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S730.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S726;
    dpdpx_3->primal_0 = _S723.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S731, DiffPair_float_0 * _S732, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_1 _S733 = *_S731;
    DiffPair_float_0 _S734 = _s_diff_ctx_8->_S599;
    DiffPair_float_0 _S735 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S736;
    (&_S736)->_S603 = _S735;
    s_primal_ctx_d_sqrt_0(&_S734, (*_S732).primal_0, &_S736);
    DiffPair_float_0 _S737 = { (*_S731).differential_0.primal_0, (*_S731).differential_0.differential_0 };
    DiffPair_1 _S738;
    (&_S738)->primal_0 = _s_diff_ctx_8->_S599;
    (&_S738)->differential_0 = _S737;
    DiffPair_float_0 _S739;
    (&_S739)->primal_0 = (*_S732).primal_0;
    (&_S739)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S740 = _S736;
    s_bwd_prop_d_sqrt_0(&_S738, &_S739, &_S740);
    DiffPair_float_0 _S741 = { _S738.differential_0.primal_0, _S738.differential_0.differential_0 };
    _S732->primal_0 = (*_S732).primal_0;
    _S732->differential_0 = _S739.differential_0;
    _S731->primal_0 = _S733.primal_0;
    _S731->differential_0 = _S741;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S742, float _s_dOut_1)
{
    DiffPair_float_0 _S743;
    (&_S743)->primal_0 = (*_S742).primal_0;
    (&_S743)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S743, _s_dOut_1);
    _S742->primal_0 = (*_S742).primal_0;
    _S742->differential_0 = _S743.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_9)
{
    DiffPair_0 _S744 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_9->_S600)->primal_0)->x) * *&((&(&_s_diff_ctx_9->_S600)->primal_0)->x) + *&((&(&_s_diff_ctx_9->_S600)->primal_0)->y) * *&((&(&_s_diff_ctx_9->_S600)->primal_0)->y);
    DiffPair_float_0 _S745 = { len_0, 0.0f };
    float2  _S746 = make_float2 (0.0f);
    float _S747 = (*dpdpx_5).differential_0.differential_0.x;
    float _S748 = _S747 + _S747;
    float _S749 = (&_s_diff_ctx_9->_S601)->differential_0 * _S748;
    float _S750 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S751 = (&_s_diff_ctx_9->_S601)->differential_0 * _S750;
    DiffPair_float_0 _S752 = { 0.0f, *&((&(&_s_diff_ctx_9->_S600)->primal_0)->x) * _S748 + *&((&(&_s_diff_ctx_9->_S600)->primal_0)->y) * _S750 };
    DiffPair_1 _S753;
    (&_S753)->primal_0 = _S745;
    (&_S753)->differential_0 = _S752;
    DiffPair_float_0 _S754;
    (&_S754)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S754)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S753, &_S754, &_s_diff_ctx_9->_S602);
    DiffPair_float_0 _S755;
    (&_S755)->primal_0 = len_0;
    (&_S755)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S755, 0.0f);
    float _S756 = _S753.differential_0.primal_0 + _S755.differential_0;
    float _S757 = *&((&(&_s_diff_ctx_9->_S600)->primal_0)->y) * _S756;
    float _S758 = _S751 + _S757 + _S757;
    float _S759 = *&((&(&_s_diff_ctx_9->_S600)->primal_0)->x) * _S756;
    float _S760 = _S749 + _S759 + _S759;
    float2  _S761 = _S746;
    *&((&_S761)->y) = _S758;
    *&((&_S761)->x) = _S760;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S744.differential_0.primal_0 + _S761, _S746 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S754.differential_0;
    dpdpx_5->primal_0 = _S744.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_2)
{
    float _S762 = (*dpdpx_7).primal_0.x;
    float _S763 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S764;
    (&_S764)->primal_0 = _S762 * _S762 + _S763 * _S763;
    (&_S764)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S764, _s_dOut_2);
    float _S765 = (*dpdpx_7).primal_0.y * _S764.differential_0;
    float _S766 = _S765 + _S765;
    float _S767 = (*dpdpx_7).primal_0.x * _S764.differential_0;
    float _S768 = _S767 + _S767;
    float2  _S769 = make_float2 (0.0f);
    *&((&_S769)->y) = _S766;
    *&((&_S769)->x) = _S768;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S769;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S770, DiffPair_float_0 * _S771, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_0 _S772 = *_S770;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S773 = _s_diff_ctx_10->_S598;
    float2  _S774 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S775 = { _S774, _S774 };
    DiffPair_float_0 _S776 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S777 = { _S776 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S778;
    (&_S778)->_S600 = _S775;
    (&_S778)->_S601 = _S776;
    (&_S778)->_S602 = _S777;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S773, (*_S771).primal_0, &_S778);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S779 = { (*_S770).differential_0.primal_0, (*_S770).differential_0.differential_0 };
    DiffPair_0 _S780;
    (&_S780)->primal_0 = _s_diff_ctx_10->_S598;
    (&_S780)->differential_0 = _S779;
    DiffPair_float_0 _S781;
    (&_S781)->primal_0 = (*_S771).primal_0;
    (&_S781)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S782 = _S778;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S780, &_S781, &_S782);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S783;
    (&_S783)->primal_0 = (&_s_diff_ctx_10->_S598)->primal_0;
    (&_S783)->differential_0 = _S774;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S783, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S784 = { _S780.differential_0.primal_0 + _S783.differential_0, _S780.differential_0.differential_0 };
    _S771->primal_0 = (*_S771).primal_0;
    _S771->differential_0 = _S781.differential_0;
    _S770->primal_0 = _S772.primal_0;
    _S770->differential_0 = _S784;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_1 _S785 = *dpdpy_1;
    DiffPair_1 _S786 = *dpdpx_8;
    float _S787 = - (&_s_diff_ctx_11->_S581)->primal_0;
    float _S788 = (&_s_diff_ctx_11->_S582)->primal_0 * (&_s_diff_ctx_11->_S582)->primal_0 + (&_s_diff_ctx_11->_S581)->primal_0 * (&_s_diff_ctx_11->_S581)->primal_0;
    float _S789 = _S788 * _S788;
    float _S790 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S789;
    float _S791 = (&_s_diff_ctx_11->_S582)->primal_0 * - _S790;
    float _S792 = (&_s_diff_ctx_11->_S581)->primal_0 * _S791;
    float _S793 = (&_s_diff_ctx_11->_S582)->primal_0 * _S791;
    float _S794 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S789;
    float _S795 = _S787 * - _S794;
    float _S796 = (&_s_diff_ctx_11->_S581)->primal_0 * _S795;
    float _S797 = (&_s_diff_ctx_11->_S582)->primal_0 * _S795;
    DiffPair_float_0 dpdpx_9 = { _S797 + _S797 + ((*dpdpx_8).differential_0.primal_0 + (_S793 + _S793 + _S788 * _S790)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S792 + _S792 + (*dpdpy_1).differential_0.primal_0 + _S796 + _S796 + - (_S788 * _S794), 0.0f };
    float _S798 = (&_s_diff_ctx_11->_S582)->primal_0 / _S788 * (*dpdpy_1).differential_0.differential_0 + _S787 / _S788 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S798;
    dpdpy_1->primal_0 = _S785.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S786.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S799, DiffPair_1 * _S800, DiffPair_float_0 * _S801, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_1 _S802 = *_S799;
    DiffPair_1 _S803 = *_S800;
    DiffPair_float_0 _S804 = _s_diff_ctx_12->_S568;
    DiffPair_float_0 _S805 = _s_diff_ctx_12->_S569;
    DiffPair_float_0 _S806 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S807;
    (&_S807)->_S581 = _S806;
    (&_S807)->_S582 = _S806;
    s_primal_ctx_d_atan2_0(&_S804, &_S805, (*_S801).primal_0, &_S807);
    DiffPair_float_0 _S808 = { (*_S800).differential_0.primal_0, (*_S800).differential_0.differential_0 };
    DiffPair_float_0 _S809 = { (*_S799).differential_0.primal_0, (*_S799).differential_0.differential_0 };
    DiffPair_1 _S810;
    (&_S810)->primal_0 = _s_diff_ctx_12->_S568;
    (&_S810)->differential_0 = _S809;
    DiffPair_1 _S811;
    (&_S811)->primal_0 = _s_diff_ctx_12->_S569;
    (&_S811)->differential_0 = _S808;
    DiffPair_float_0 _S812;
    (&_S812)->primal_0 = (*_S801).primal_0;
    (&_S812)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S813 = _S807;
    s_bwd_prop_d_atan2_0(&_S810, &_S811, &_S812, &_S813);
    DiffPair_float_0 _S814 = { _S811.differential_0.primal_0, _S811.differential_0.differential_0 };
    DiffPair_float_0 _S815 = { _S810.differential_0.primal_0, _S810.differential_0.differential_0 };
    _S801->primal_0 = (*_S801).primal_0;
    _S801->differential_0 = _S812.differential_0;
    _S799->primal_0 = _S802.primal_0;
    _S799->differential_0 = _S815;
    _S800->primal_0 = _S803.primal_0;
    _S800->differential_0 = _S814;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S816, DiffPair_float_0 * _S817, float _s_dOut_3)
{
    DiffPair_float_0 _S818;
    (&_S818)->primal_0 = (*_S816).primal_0;
    (&_S818)->differential_0 = 0.0f;
    DiffPair_float_0 _S819;
    (&_S819)->primal_0 = (*_S817).primal_0;
    (&_S819)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S818, &_S819, _s_dOut_3);
    _S817->primal_0 = (*_S817).primal_0;
    _S817->differential_0 = _S819.differential_0;
    _S816->primal_0 = (*_S816).primal_0;
    _S816->differential_0 = _S818.differential_0;
    return;
}

inline __device__ void s_bwd_prop_CameraDistortion_x24init_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_2, CameraDistortion_0 * _s_dOut_4)
{
    float2  _S820 = _s_dOut_4->thin_prism_coeffs_0;
    dpthin_prism_coeffs_2->primal_0 = (*dpthin_prism_coeffs_2).primal_0;
    dpthin_prism_coeffs_2->differential_0 = _S820;
    float2  _S821 = _s_dOut_4->tangential_coeffs_0;
    dptangential_coeffs_2->primal_0 = (*dptangential_coeffs_2).primal_0;
    dptangential_coeffs_2->differential_0 = _S821;
    float4  _S822 = _s_dOut_4->radial_coeffs_0;
    dpradial_coeffs_2->primal_0 = (*dpradial_coeffs_2).primal_0;
    dpradial_coeffs_2->differential_0 = _S822;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_3, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S823 = *dpcov3d_1;
    DiffPair_float_0 _S824 = *dpfx_1;
    DiffPair_float_0 _S825 = *dpfy_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S826 = *dpradial_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S827 = *dptangential_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S828 = *dpthin_prism_coeffs_3;
    float2  _S829 = make_float2 (0.0f);
    CameraDistortion_0 _S830 = s_primal_ctx_CameraDistortion_x24init_0((*dpradial_coeffs_3).primal_0, (*dptangential_coeffs_3).primal_0, (*dpthin_prism_coeffs_3).primal_0);
    float2  _S831 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S832 = length_0(_S831);
    float _S833 = (*dpmean3d_1).primal_0.z;
    float _S834 = s_primal_ctx_atan2_0(_S832, _S833);
    bool _S835 = _S834 < 0.00100000004749745f;
    float k_2;
    float _S836;
    float _S837;
    float _S838;
    if(_S835)
    {
        float _S839 = 1.0f - _S834 * _S834 / 3.0f;
        float _S840 = _S833 * _S833;
        k_2 = _S839 / _S833;
        _S836 = 0.0f;
        _S837 = _S840;
        _S838 = _S839;
    }
    else
    {
        float _S841 = _S832 * _S832;
        k_2 = _S834 / _S832;
        _S836 = _S841;
        _S837 = 0.0f;
        _S838 = 0.0f;
    }
    float2  _S842 = make_float2 (k_2);
    float2  _S843 = _S831 * make_float2 (k_2);
    float k1_2 = _S830.radial_coeffs_0.x;
    float k2_2 = _S830.radial_coeffs_0.y;
    float k3_2 = _S830.radial_coeffs_0.z;
    float k4_2 = _S830.radial_coeffs_0.w;
    float p1_2 = _S830.tangential_coeffs_0.x;
    float p2_2 = _S830.tangential_coeffs_0.y;
    float sx1_2 = _S830.thin_prism_coeffs_0.x;
    float sy1_2 = _S830.thin_prism_coeffs_0.y;
    float u_6 = _S843.x;
    float v_6 = _S843.y;
    float r2_6 = u_6 * u_6 + v_6 * v_6;
    float _S844 = k3_2 + r2_6 * k4_2;
    float _S845 = k2_2 + r2_6 * _S844;
    float _S846 = k1_2 + r2_6 * _S845;
    float radial_0 = 1.0f + r2_6 * _S846;
    float2  _S847 = make_float2 (radial_0);
    float _S848 = 2.0f * p1_2;
    float _S849 = _S848 * u_6;
    float _S850 = 2.0f * u_6;
    float _S851 = r2_6 + _S850 * u_6;
    float _S852 = 2.0f * p2_2;
    float _S853 = _S852 * u_6;
    float _S854 = 2.0f * v_6;
    float _S855 = r2_6 + _S854 * v_6;
    float2  _S856 = _S843 * make_float2 (radial_0) + make_float2 (_S849 * v_6 + p2_2 * _S851 + sx1_2 * r2_6, _S853 * v_6 + p1_2 * _S855 + sy1_2 * r2_6);
    float _S857 = _S856.x;
    float _S858 = _S856.y;
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float _S859 = s_primal_ctx_s_primal_ctx_atan2_0(_S832, _S833);
    bool _S860 = _S859 < 0.00100000004749745f;
    float _S861;
    float _S862;
    float _S863;
    if(_S860)
    {
        float _S864 = 1.0f - _S859 * _S859 / 3.0f;
        float _S865 = _S833 * _S833;
        k_2 = _S864 / _S833;
        _S861 = 0.0f;
        _S862 = _S865;
        _S863 = _S864;
    }
    else
    {
        float _S866 = _S832 * _S832;
        k_2 = _S859 / _S832;
        _S861 = _S866;
        _S862 = 0.0f;
        _S863 = 0.0f;
    }
    float2  _S867 = make_float2 (k_2);
    float2  _S868 = _S831 * make_float2 (k_2);
    float u_7 = _S868.x;
    float v_7 = _S868.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S869 = k3_2 + r2_7 * k4_2;
    float _S870 = k2_2 + r2_7 * _S869;
    float _S871 = k1_2 + r2_7 * _S870;
    float2  _S872 = make_float2 (1.0f + r2_7 * _S871);
    float _S873 = _S848 * u_7;
    float _S874 = 2.0f * u_7;
    float2  _S875 = make_float2 (_S824.primal_0, 0.0f);
    float2  _S876 = _S868 * _S875;
    float _S877 = p2_2 * _S824.primal_0;
    float _S878 = v_7 * _S824.primal_0;
    float _S879 = _S876.x + _S876.y;
    float _S880 = r2_7 * _S879;
    float _S881 = r2_7 * _S880;
    float _S882 = r2_7 * _S881;
    float _S883 = sx1_2 * _S824.primal_0 + _S877 + _S871 * _S879 + _S870 * _S880 + _S869 * _S881 + k4_2 * _S882;
    float _S884 = v_7 * _S883;
    float _S885 = u_7 * _S883;
    float2  _S886 = _S872 * _S875 + make_float2 (_S874 * _S877 + 2.0f * (u_7 * _S877) + _S848 * _S878 + _S885 + _S885, _S873 * _S824.primal_0 + _S884 + _S884);
    float2  _S887 = _S831 * _S886;
    float2  _S888 = _S867 * _S886;
    float _S889 = _S887.x + _S887.y;
    float k_3;
    float _S890;
    float _S891;
    float _S892;
    float _S893;
    float _S894;
    float _S895;
    float _S896;
    float _S897;
    if(_S860)
    {
        float _S898 = _S889 / _S862;
        float _S899 = _S862 * _S862;
        float _S900 = - _S898;
        float _S901 = _S863 * _S900;
        float _S902 = 0.3333333432674408f * - (_S833 * _S898);
        float _S903 = _S859 * _S902;
        k_2 = _S903 + _S903;
        k_3 = _S901;
        _S890 = 0.0f;
        _S891 = 0.0f;
        _S892 = 0.0f;
        _S893 = 0.0f;
        _S894 = _S902;
        _S895 = _S898;
        _S896 = _S900;
        _S897 = _S899;
    }
    else
    {
        float _S904 = _S889 / _S861;
        float _S905 = _S861 * _S861;
        float _S906 = - _S904;
        float _S907 = _S859 * _S906;
        k_2 = _S832 * _S904;
        k_3 = 0.0f;
        _S890 = _S907;
        _S891 = _S904;
        _S892 = _S906;
        _S893 = _S905;
        _S894 = 0.0f;
        _S895 = 0.0f;
        _S896 = 0.0f;
        _S897 = 0.0f;
    }
    DiffPair_float_0 _S908 = { _S832, 0.0f };
    DiffPair_float_0 _S909 = { _S833, 0.0f };
    float _S910 = (&_s_diff_ctx_13->_S571)->differential_0 + k_3;
    float _S911 = (&_s_diff_ctx_13->_S570)->differential_0 + _S890;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S912 = { _S831, _S829 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S913;
    (&_S913)->primal_0 = _S831;
    (&_S913)->differential_0 = _S829;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S914 = { _S829, _S829 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S915;
    (&_S915)->_S598 = _S914;
    s_primal_ctx_s_bwd_length_impl_0(&_S913, _S911, &_S915);
    float2  _S916 = _S913.differential_0 + _S888;
    float3  _S917 = make_float3 (_S916.x, _S916.y, _S910);
    Matrix<float, 2, 3>  _S918 = J_10;
    _S918[int(0)] = _S917;
    float _S919;
    float _S920;
    if(_S860)
    {
        float _S921 = 1.0f - _S859 * _S859 / 3.0f;
        float _S922 = _S833 * _S833;
        k_3 = _S921 / _S833;
        _S890 = 0.0f;
        _S919 = _S922;
        _S920 = _S921;
    }
    else
    {
        float _S923 = _S832 * _S832;
        k_3 = _S859 / _S832;
        _S890 = _S923;
        _S919 = 0.0f;
        _S920 = 0.0f;
    }
    float2  _S924 = make_float2 (k_3);
    float2  _S925 = _S831 * make_float2 (k_3);
    float u_8 = _S925.x;
    float v_8 = _S925.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S926 = k3_2 + r2_8 * k4_2;
    float _S927 = k2_2 + r2_8 * _S926;
    float _S928 = k1_2 + r2_8 * _S927;
    float2  _S929 = make_float2 (1.0f + r2_8 * _S928);
    float _S930 = _S852 * u_8;
    float _S931 = 2.0f * v_8;
    float2  _S932 = make_float2 (0.0f, _S825.primal_0);
    float2  _S933 = _S925 * _S932;
    float _S934 = p1_2 * _S825.primal_0;
    float _S935 = v_8 * _S825.primal_0;
    float _S936 = _S933.x + _S933.y;
    float _S937 = r2_8 * _S936;
    float _S938 = r2_8 * _S937;
    float _S939 = r2_8 * _S938;
    float _S940 = sy1_2 * _S825.primal_0 + _S934 + _S928 * _S936 + _S927 * _S937 + _S926 * _S938 + k4_2 * _S939;
    float _S941 = v_8 * _S940;
    float _S942 = u_8 * _S940;
    float2  _S943 = _S929 * _S932 + make_float2 (_S852 * _S935 + _S942 + _S942, _S931 * _S934 + 2.0f * (v_8 * _S934) + _S930 * _S825.primal_0 + _S941 + _S941);
    float2  _S944 = _S831 * _S943;
    float2  _S945 = _S924 * _S943;
    float _S946 = _S944.x + _S944.y;
    float _S947;
    float _S948;
    float _S949;
    float _S950;
    float _S951;
    float _S952;
    float _S953;
    float _S954;
    float _S955;
    if(_S860)
    {
        float _S956 = _S946 / _S919;
        float _S957 = _S919 * _S919;
        float _S958 = - _S956;
        float _S959 = _S920 * _S958;
        float _S960 = 0.3333333432674408f * - (_S833 * _S956);
        float _S961 = _S859 * _S960;
        k_3 = _S961 + _S961;
        _S947 = _S959;
        _S948 = 0.0f;
        _S949 = 0.0f;
        _S950 = 0.0f;
        _S951 = 0.0f;
        _S952 = _S960;
        _S953 = _S956;
        _S954 = _S958;
        _S955 = _S957;
    }
    else
    {
        float _S962 = _S946 / _S890;
        float _S963 = _S890 * _S890;
        float _S964 = - _S962;
        float _S965 = _S859 * _S964;
        k_3 = _S832 * _S962;
        _S947 = 0.0f;
        _S948 = _S965;
        _S949 = _S962;
        _S950 = _S964;
        _S951 = _S963;
        _S952 = 0.0f;
        _S953 = 0.0f;
        _S954 = 0.0f;
        _S955 = 0.0f;
    }
    float _S966 = (&_s_diff_ctx_13->_S574)->differential_0 + _S947;
    float _S967 = (&_s_diff_ctx_13->_S573)->differential_0 + _S948;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S968;
    (&_S968)->primal_0 = _S831;
    (&_S968)->differential_0 = _S829;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S969;
    (&_S969)->_S598 = _S914;
    s_primal_ctx_s_bwd_length_impl_0(&_S968, _S967, &_S969);
    float2  _S970 = _S968.differential_0 + _S945;
    float3  _S971 = make_float3 (_S970.x, _S970.y, _S966);
    _S918[int(1)] = _S971;
    Matrix<float, 3, 2>  _S972 = transpose_1(_S918);
    CameraDistortion_0 _S973 = CameraDistortion_x24_syn_dzero_0();
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S974;
    (&_S974)->primal_0 = s_primal_ctx_mul_2(_S918, _S823.primal_0);
    (&_S974)->differential_0 = J_10;
    Matrix<float, 3, 2>  _S975 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S976;
    (&_S976)->primal_0 = _S972;
    (&_S976)->differential_0 = _S975;
    s_bwd_prop_mul_1(&_S974, &_S976, dpcov2d_1);
    Matrix<float, 2, 3>  _S977 = transpose_2(_S976.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S978;
    (&_S978)->primal_0 = _S918;
    (&_S978)->differential_0 = J_10;
    Matrix<float, 3, 3>  _S979 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S980;
    (&_S980)->primal_0 = _S823.primal_0;
    (&_S980)->differential_0 = _S979;
    s_bwd_prop_mul_2(&_S978, &_S980, _S974.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S981 = _S980;
    Matrix<float, 2, 3>  _S982 = _S977 + _S978.differential_0;
    float2  _S983 = _S829;
    *&((&_S983)->y) = _S982.rows[int(1)].y;
    *&((&_S983)->x) = _S982.rows[int(1)].x;
    float2  _S984 = _S983;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S985 = { _S829, _S983 };
    DiffPair_0 _S986;
    (&_S986)->primal_0 = _S912;
    (&_S986)->differential_0 = _S985;
    DiffPair_float_0 _S987;
    (&_S987)->primal_0 = _S967;
    (&_S987)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S988 = _S969;
    s_bwd_prop_s_bwd_length_impl_0(&_S986, &_S987, &_S988);
    DiffPair_0 _S989 = _S986;
    DiffPair_float_0 _S990 = _S987;
    DiffPair_float_0 _S991 = { 0.0f, _S982.rows[int(1)].z };
    DiffPair_float_0 _S992 = { 0.0f, _S987.differential_0 };
    DiffPair_1 _S993;
    (&_S993)->primal_0 = _S908;
    (&_S993)->differential_0 = _S992;
    DiffPair_1 _S994;
    (&_S994)->primal_0 = _S909;
    (&_S994)->differential_0 = _S991;
    DiffPair_float_0 _S995;
    (&_S995)->primal_0 = k_3;
    (&_S995)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S993, &_S994, &_S995, &_s_diff_ctx_13->_S575);
    DiffPair_1 _S996 = _S993;
    DiffPair_1 _S997 = _S994;
    DiffPair_float_0 _S998 = _S995;
    if(_S860)
    {
        float _S999 = _S998.differential_0 + _S998.differential_0;
        float _S1000 = _S952 * _S999;
        float _S1001 = - (0.3333333432674408f * (_S859 * _S999));
        float _S1002 = _S954 * _S982.rows[int(1)].z;
        float _S1003 = (_S833 * _S1001 + - (_S920 * _S982.rows[int(1)].z)) / _S955;
        float _S1004 = _S946 * - _S1003;
        float _S1005 = _S953 * _S1001 + _S997.differential_0.primal_0;
        k_3 = _S919 * _S1003;
        _S947 = 0.0f;
        _S948 = _S1004;
        _S949 = _S1002;
        _S950 = _S996.differential_0.primal_0;
        _S951 = _S1000;
        _S952 = _S1005;
    }
    else
    {
        float _S1006 = _S950 * _S990.differential_0;
        float _S1007 = (_S832 * _S998.differential_0 + - (_S859 * _S990.differential_0)) / _S951;
        float _S1008 = _S946 * - _S1007;
        float _S1009 = _S949 * _S998.differential_0 + _S996.differential_0.primal_0;
        k_3 = _S890 * _S1007;
        _S947 = _S1008;
        _S948 = 0.0f;
        _S949 = 0.0f;
        _S950 = _S1009;
        _S951 = _S1006;
        _S952 = _S997.differential_0.primal_0;
    }
    float2  _S1010 = _S924 * _S984;
    float2  _S1011 = _S943 * _S984;
    float2  _S1012 = _S829;
    *&((&_S1012)->y) = k_3;
    *&((&_S1012)->x) = k_3;
    float2  _S1013 = _S943 * _S1012;
    float2  _S1014 = _S1010 + _S831 * _S1012;
    float _S1015 = _S1014.x;
    float _S1016 = _S1015 + _S1015;
    float _S1017 = _S940 * _S1016;
    float _S1018 = _S1014.y + _S1014.y;
    float _S1019 = _S940 * _S1018;
    float _S1020 = u_8 * _S1016 + v_8 * _S1018;
    float _S1021 = k4_2 * _S1020;
    float _S1022 = _S939 * _S1020;
    float _S1023 = _S938 * _S1020;
    float _S1024 = _S938 * _S1021;
    float _S1025 = _S937 * _S1020;
    float _S1026 = _S926 * _S1020 + r2_8 * _S1021;
    float _S1027 = _S937 * _S1026;
    float _S1028 = _S936 * _S1020;
    float _S1029 = _S927 * _S1020 + r2_8 * _S1026;
    float _S1030 = _S936 * _S1029;
    float _S1031 = _S928 * _S1020 + r2_8 * _S1029;
    float _S1032 = _S852 * _S1014.x;
    float _S1033 = _S935 * _S1014.x;
    float _S1034 = v_8 * _S1032;
    float _S1035 = _S825.primal_0 * _S1032;
    float _S1036 = _S930 * _S1014.y;
    float _S1037 = _S825.primal_0 * _S1014.y;
    float _S1038 = 2.0f * _S1014.y;
    float _S1039 = _S934 * _S1038;
    float _S1040 = _S934 * _S1014.y;
    float _S1041 = _S1020 + v_8 * _S1038 + _S931 * _S1014.y;
    float _S1042 = p1_2 * _S1041;
    float _S1043 = _S825.primal_0 * _S1041;
    float _S1044 = sy1_2 * _S1020;
    float _S1045 = _S825.primal_0 * _S1020;
    float2  _S1046 = _S929 * _S1014;
    float2  _S1047 = _S932 * _S1014;
    float2  _S1048 = _S829;
    *&((&_S1048)->y) = _S1031;
    *&((&_S1048)->x) = _S1031;
    float _S1049 = _S1047.x + _S1047.y;
    float _S1050 = _S1028 + r2_8 * _S1049;
    float _S1051 = _S1025 + r2_8 * _S1050;
    float _S1052 = _S1023 + r2_8 * _S1051;
    float _S1053 = _S1024 + _S1027 + _S1030 + _S928 * _S1049 + _S927 * _S1050 + _S926 * _S1051 + k4_2 * _S1052;
    float _S1054 = v_8 * _S1053;
    float _S1055 = u_8 * _S1053;
    float2  _S1056 = _S1013 + _S989.differential_0.primal_0;
    float2  _S1057 = _S932 * _S1048 + make_float2 (_S1017 + _S852 * _S1037 + _S1055 + _S1055, _S1019 + _S1035 + _S1039 + 2.0f * _S1040 + _S1054 + _S1054);
    float _S1058 = _S1034 + _S1036 + _S1042 + _S1044 + (_S1046 + _S925 * _S1048).y;
    float _S1059 = _S1033 + u_8 * _S1037;
    float _S1060 = _S1022 + r2_8 * _S1052;
    float2  _S1061 = _S831 * _S1057;
    float _S1062 = _S1011.x + _S1011.y + _S1061.x + _S1061.y;
    float2  _S1063 = _S924 * _S1057 + _S1056;
    if(_S860)
    {
        float _S1064 = _S833 * _S948;
        float _S1065 = _S1062 / _S919;
        float _S1066 = _S859 * (0.3333333432674408f * - (_S949 + _S833 * _S1065));
        float _S1067 = _S1064 + _S1064 + _S920 * - _S1065 + _S952;
        k_3 = _S1066 + _S1066 + _S951;
        _S890 = _S1067;
        _S919 = _S950;
    }
    else
    {
        float _S1068 = _S832 * _S947;
        float _S1069 = _S1062 / _S890;
        float _S1070 = _S1068 + _S1068 + _S859 * - _S1069 + _S950;
        k_3 = _S832 * _S1069 + _S951;
        _S890 = _S952;
        _S919 = _S1070;
    }
    DiffPair_float_0 _S1071;
    (&_S1071)->primal_0 = _S832;
    (&_S1071)->differential_0 = 0.0f;
    DiffPair_float_0 _S1072;
    (&_S1072)->primal_0 = _S833;
    (&_S1072)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1071, &_S1072, k_3);
    float _S1073 = _S1072.differential_0 + _S890;
    float _S1074 = _S1071.differential_0 + _S919;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1075;
    (&_S1075)->primal_0 = _S831;
    (&_S1075)->differential_0 = _S829;
    s_bwd_length_impl_0(&_S1075, _S1074);
    float2  _S1076 = _S1075.differential_0 + _S1063;
    float3  _S1077 = make_float3 (_S1076.x, _S1076.y, _S1073);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1078;
    (&_S1078)->primal_0 = _S831;
    (&_S1078)->differential_0 = _S829;
    s_bwd_length_impl_0(&_S1078, 0.0f);
    float3  _S1079 = _S1077 + make_float3 (_S1078.differential_0.x, _S1078.differential_0.y, 0.0f);
    float2  _S1080 = _S829;
    *&((&_S1080)->y) = _S982.rows[int(0)].y;
    *&((&_S1080)->x) = _S982.rows[int(0)].x;
    float2  _S1081 = _S1080;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1082 = { _S829, _S1080 };
    DiffPair_0 _S1083;
    (&_S1083)->primal_0 = _S912;
    (&_S1083)->differential_0 = _S1082;
    DiffPair_float_0 _S1084;
    (&_S1084)->primal_0 = _S911;
    (&_S1084)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1085 = _S915;
    s_bwd_prop_s_bwd_length_impl_0(&_S1083, &_S1084, &_S1085);
    DiffPair_0 _S1086 = _S1083;
    DiffPair_float_0 _S1087 = _S1084;
    DiffPair_float_0 _S1088 = { 0.0f, _S982.rows[int(0)].z };
    DiffPair_float_0 _S1089 = { 0.0f, _S1084.differential_0 };
    DiffPair_1 _S1090;
    (&_S1090)->primal_0 = _S908;
    (&_S1090)->differential_0 = _S1089;
    DiffPair_1 _S1091;
    (&_S1091)->primal_0 = _S909;
    (&_S1091)->differential_0 = _S1088;
    DiffPair_float_0 _S1092;
    (&_S1092)->primal_0 = k_2;
    (&_S1092)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1090, &_S1091, &_S1092, &_s_diff_ctx_13->_S572);
    DiffPair_1 _S1093 = _S1090;
    DiffPair_1 _S1094 = _S1091;
    DiffPair_float_0 _S1095 = _S1092;
    if(_S860)
    {
        float _S1096 = _S1095.differential_0 + _S1095.differential_0;
        float _S1097 = _S894 * _S1096;
        float _S1098 = - (0.3333333432674408f * (_S859 * _S1096));
        float _S1099 = _S896 * _S982.rows[int(0)].z;
        float _S1100 = (_S833 * _S1098 + - (_S863 * _S982.rows[int(0)].z)) / _S897;
        float _S1101 = _S889 * - _S1100;
        float _S1102 = _S895 * _S1098 + _S1094.differential_0.primal_0;
        k_2 = _S862 * _S1100;
        k_3 = 0.0f;
        _S890 = _S1101;
        _S891 = _S1099;
        _S892 = _S1093.differential_0.primal_0;
        _S893 = _S1097;
        _S894 = _S1102;
    }
    else
    {
        float _S1103 = _S892 * _S1087.differential_0;
        float _S1104 = (_S832 * _S1095.differential_0 + - (_S859 * _S1087.differential_0)) / _S893;
        float _S1105 = _S889 * - _S1104;
        float _S1106 = _S891 * _S1095.differential_0 + _S1093.differential_0.primal_0;
        k_2 = _S861 * _S1104;
        k_3 = _S1105;
        _S890 = 0.0f;
        _S891 = 0.0f;
        _S892 = _S1106;
        _S893 = _S1103;
        _S894 = _S1094.differential_0.primal_0;
    }
    float2  _S1107 = _S867 * _S1081;
    float2  _S1108 = _S886 * _S1081;
    float2  _S1109 = _S829;
    *&((&_S1109)->y) = k_2;
    *&((&_S1109)->x) = k_2;
    float2  _S1110 = _S886 * _S1109;
    float2  _S1111 = _S1107 + _S831 * _S1109;
    float _S1112 = _S1111.x;
    float _S1113 = _S1112 + _S1112;
    float _S1114 = _S883 * _S1113;
    float _S1115 = _S1111.y + _S1111.y;
    float _S1116 = _S883 * _S1115;
    float _S1117 = u_7 * _S1113 + v_7 * _S1115;
    float _S1118 = k4_2 * _S1117;
    float _S1119 = _S882 * _S1117;
    float _S1120 = _S881 * _S1117;
    float _S1121 = _S881 * _S1118;
    float _S1122 = _S880 * _S1117;
    float _S1123 = _S869 * _S1117 + r2_7 * _S1118;
    float _S1124 = _S880 * _S1123;
    float _S1125 = _S879 * _S1117;
    float _S1126 = _S870 * _S1117 + r2_7 * _S1123;
    float _S1127 = _S879 * _S1126;
    float _S1128 = _S871 * _S1117 + r2_7 * _S1126;
    float _S1129 = _S848 * _S1111.x;
    float _S1130 = _S878 * _S1111.x;
    float _S1131 = v_7 * _S1129;
    float _S1132 = _S824.primal_0 * _S1129;
    float _S1133 = _S873 * _S1111.y;
    float _S1134 = _S824.primal_0 * _S1111.y;
    float _S1135 = 2.0f * _S1111.x;
    float _S1136 = _S877 * _S1135;
    float _S1137 = _S877 * _S1111.x;
    float _S1138 = _S1117 + u_7 * _S1135 + _S874 * _S1111.x;
    float _S1139 = p2_2 * _S1138;
    float _S1140 = _S824.primal_0 * _S1138;
    float _S1141 = sx1_2 * _S1117;
    float _S1142 = _S824.primal_0 * _S1117;
    float2  _S1143 = _S872 * _S1111;
    float2  _S1144 = _S875 * _S1111;
    float2  _S1145 = _S829;
    *&((&_S1145)->y) = _S1128;
    *&((&_S1145)->x) = _S1128;
    float _S1146 = _S1144.x + _S1144.y;
    float _S1147 = _S1125 + r2_7 * _S1146;
    float _S1148 = _S1122 + r2_7 * _S1147;
    float _S1149 = _S1120 + r2_7 * _S1148;
    float _S1150 = _S1121 + _S1124 + _S1127 + _S871 * _S1146 + _S870 * _S1147 + _S869 * _S1148 + k4_2 * _S1149;
    float _S1151 = v_7 * _S1150;
    float _S1152 = u_7 * _S1150;
    float2  _S1153 = _S1110 + _S1086.differential_0.primal_0;
    float _S1154 = _S1130 + u_7 * _S1134;
    float _S1155 = _S1147 + _S1050;
    float _S1156 = _S1149 + _S1052;
    float _S1157 = _S1131 + _S1133 + _S1139 + _S1141 + (_S1143 + _S868 * _S1145).x;
    float2  _S1158 = _S875 * _S1145 + make_float2 (_S1114 + _S1136 + 2.0f * _S1137 + _S848 * _S1134 + _S1152 + _S1152, _S1116 + _S1132 + _S1151 + _S1151);
    float _S1159 = _S1148 + _S1051;
    float _S1160 = _S1119 + r2_7 * _S1149 + _S1060;
    float2  _S1161 = _S831 * _S1158;
    float _S1162 = _S1108.x + _S1108.y + _S1161.x + _S1161.y;
    float2  _S1163 = _S867 * _S1158 + _S1153;
    if(_S860)
    {
        float _S1164 = _S833 * _S890;
        float _S1165 = _S1162 / _S862;
        float _S1166 = _S859 * (0.3333333432674408f * - (_S891 + _S833 * _S1165));
        float _S1167 = _S1164 + _S1164 + _S863 * - _S1165 + _S894;
        k_2 = _S1166 + _S1166 + _S893;
        _S861 = _S1167;
        _S862 = _S892;
    }
    else
    {
        float _S1168 = _S832 * k_3;
        float _S1169 = _S1162 / _S861;
        float _S1170 = _S1168 + _S1168 + _S859 * - _S1169 + _S892;
        k_2 = _S832 * _S1169 + _S893;
        _S861 = _S894;
        _S862 = _S1170;
    }
    DiffPair_float_0 _S1171;
    (&_S1171)->primal_0 = _S832;
    (&_S1171)->differential_0 = 0.0f;
    DiffPair_float_0 _S1172;
    (&_S1172)->primal_0 = _S833;
    (&_S1172)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1171, &_S1172, k_2);
    float _S1173 = _S1172.differential_0 + _S861;
    float _S1174 = _S1171.differential_0 + _S862;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1175;
    (&_S1175)->primal_0 = _S831;
    (&_S1175)->differential_0 = _S829;
    s_bwd_length_impl_0(&_S1175, _S1174);
    float2  _S1176 = _S1175.differential_0 + _S1163;
    float3  _S1177 = make_float3 (_S1176.x, _S1176.y, _S1173);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1178;
    (&_S1178)->primal_0 = _S831;
    (&_S1178)->differential_0 = _S829;
    s_bwd_length_impl_0(&_S1178, 0.0f);
    float _S1179 = _S825.primal_0 * dpmean2d_1.y;
    float _S1180 = _S824.primal_0 * dpmean2d_1.x;
    float2  _S1181 = make_float2 (_S1180, _S1179);
    float2  _S1182 = _S843 * _S1181;
    float _S1183 = p1_2 * _S1179;
    float _S1184 = v_6 * _S1179;
    float _S1185 = p2_2 * _S1180;
    float _S1186 = v_6 * _S1180;
    float _S1187 = _S1182.x + _S1182.y;
    float _S1188 = r2_6 * _S1187;
    float _S1189 = r2_6 * _S1188;
    float _S1190 = r2_6 * _S1189;
    float _S1191 = sy1_2 * _S1179 + _S1183 + sx1_2 * _S1180 + _S1185 + _S846 * _S1187 + _S845 * _S1188 + _S844 * _S1189 + k4_2 * _S1190;
    float _S1192 = v_6 * _S1191;
    float _S1193 = u_6 * _S1191;
    float2  _S1194 = make_float2 (r2_6 * _S1180 + _S1142, r2_6 * _S1179 + _S1045);
    float2  _S1195 = make_float2 (_S855 * _S1179 + 2.0f * (u_6 * _S1186 + _S1154) + _S1043, 2.0f * (u_6 * _S1184 + _S1059) + _S851 * _S1180 + _S1140);
    float4  _S1196 = make_float4 (_S1188 + _S1155, _S1189 + _S1159, _S1190 + _S1156, r2_6 * _S1190 + _S1160);
    float3  _S1197 = _S1177 + make_float3 (_S1178.differential_0.x, _S1178.differential_0.y, 0.0f) + _S1079;
    float _S1198 = _S857 * dpmean2d_1.x + _S1157;
    float _S1199 = _S858 * dpmean2d_1.y + _S1058;
    float2  _S1200 = _S847 * _S1181 + make_float2 (_S852 * _S1184 + _S850 * _S1185 + 2.0f * (u_6 * _S1185) + _S848 * _S1186 + _S1193 + _S1193, _S854 * _S1183 + 2.0f * (v_6 * _S1183) + _S853 * _S1179 + _S849 * _S1180 + _S1192 + _S1192);
    CameraDistortion_0 _S1201 = _S973;
    (&_S1201)->thin_prism_coeffs_0 = _S1194;
    (&_S1201)->tangential_coeffs_0 = _S1195;
    (&_S1201)->radial_coeffs_0 = _S1196;
    CameraDistortion_0 _S1202 = _S973;
    CameraDistortion_0 _S1203 = _S1201;
    CameraDistortion_0 _S1204 = CameraDistortion_x24_syn_dadd_0(&_S1202, &_S1203);
    float2  _S1205 = _S831 * _S1200;
    float2  _S1206 = _S842 * _S1200;
    float _S1207 = _S1205.x + _S1205.y;
    if(_S835)
    {
        float _S1208 = _S1207 / _S837;
        float _S1209 = _S838 * - _S1208;
        float _S1210 = _S834 * (0.3333333432674408f * - (_S833 * _S1208));
        k_2 = _S1210 + _S1210;
        _S836 = _S1209;
        _S837 = 0.0f;
    }
    else
    {
        float _S1211 = _S1207 / _S836;
        float _S1212 = _S834 * - _S1211;
        k_2 = _S832 * _S1211;
        _S836 = 0.0f;
        _S837 = _S1212;
    }
    DiffPair_float_0 _S1213;
    (&_S1213)->primal_0 = _S832;
    (&_S1213)->differential_0 = 0.0f;
    DiffPair_float_0 _S1214;
    (&_S1214)->primal_0 = _S833;
    (&_S1214)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1213, &_S1214, k_2);
    float _S1215 = _S1214.differential_0 + _S836;
    float _S1216 = _S1213.differential_0 + _S837;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1217;
    (&_S1217)->primal_0 = _S831;
    (&_S1217)->differential_0 = _S829;
    s_bwd_length_impl_0(&_S1217, _S1216);
    float2  _S1218 = _S1217.differential_0 + _S1206;
    float4  _S1219 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1220;
    (&_S1220)->primal_0 = _S826.primal_0;
    (&_S1220)->differential_0 = _S1219;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1221;
    (&_S1221)->primal_0 = _S827.primal_0;
    (&_S1221)->differential_0 = _S829;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1222;
    (&_S1222)->primal_0 = _S828.primal_0;
    (&_S1222)->differential_0 = _S829;
    CameraDistortion_0 _S1223 = _S1204;
    s_bwd_prop_CameraDistortion_x24init_0(&_S1220, &_S1221, &_S1222, &_S1223);
    dpthin_prism_coeffs_3->primal_0 = (*dpthin_prism_coeffs_3).primal_0;
    dpthin_prism_coeffs_3->differential_0 = _S1222.differential_0;
    dptangential_coeffs_3->primal_0 = (*dptangential_coeffs_3).primal_0;
    dptangential_coeffs_3->differential_0 = _S1221.differential_0;
    dpradial_coeffs_3->primal_0 = (*dpradial_coeffs_3).primal_0;
    dpradial_coeffs_3->differential_0 = _S1220.differential_0;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = dpmean2d_1.y;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = dpmean2d_1.x;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S1199;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S1198;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S981.differential_0;
    float3  _S1224 = _S1197 + make_float3 (_S1218.x, _S1218.y, _S1215);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1224;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_11, float fy_11, float cx_11, float cy_11, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, uint image_width_7, uint image_height_7, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    Matrix<float, 2, 2>  _S1225 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1226 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1227 = { _S1226, _S1226 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1228 = { _S1226, _S1226, _S1227, _S1226, _S1226, _S1227 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1229;
    (&_S1229)->_S576 = _S1225;
    (&_S1229)->_S577 = _S1228;
    float3  mean_c_7 = s_primal_ctx_mul_0(R_9, mean_7) + t_8;
    float3  _S1230 = s_primal_ctx_exp_0(scale_9);
    float _S1231 = quat_10.y;
    float _S1232 = _S1231 * _S1231 + quat_10.z * quat_10.z + quat_10.w * quat_10.w + quat_10.x * quat_10.x;
    float _S1233 = s_primal_ctx_rsqrt_0(_S1232);
    float x_35 = quat_10.y * _S1233;
    float y_20 = quat_10.z * _S1233;
    float z_17 = quat_10.w * _S1233;
    float w_10 = quat_10.x * _S1233;
    float x2_10 = x_35 * x_35;
    float y2_10 = y_20 * y_20;
    float z2_17 = z_17 * z_17;
    float xy_10 = x_35 * y_20;
    float xz_10 = x_35 * z_17;
    float yz_10 = y_20 * z_17;
    float wx_10 = w_10 * x_35;
    float wy_10 = w_10 * y_20;
    float wz_10 = w_10 * z_17;
    Matrix<float, 3, 3>  _S1234 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1230.x, 0.0f, 0.0f, 0.0f, _S1230.y, 0.0f, 0.0f, 0.0f, _S1230.z);
    Matrix<float, 3, 3>  _S1235 = s_primal_ctx_mul_1(_S1234, S_1);
    Matrix<float, 3, 3>  _S1236 = transpose_0(_S1235);
    Matrix<float, 3, 3>  _S1237 = s_primal_ctx_mul_1(_S1235, _S1236);
    Matrix<float, 3, 3>  _S1238 = s_primal_ctx_mul_1(R_9, _S1237);
    Matrix<float, 3, 3>  _S1239 = transpose_0(R_9);
    Matrix<float, 3, 3>  _S1240 = s_primal_ctx_mul_1(_S1238, _S1239);
    Matrix<float, 2, 2>  _S1241 = _S1225;
    float2  _S1242 = make_float2 (0.0f);
    float2  _S1243 = _S1242;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_7, _S1240, fx_11, fy_11, cx_11, cy_11, radial_coeffs_10, tangential_coeffs_10, thin_prism_coeffs_10, &_S1241, &_S1243, &(&_S1229)->_S577);
    (&_S1229)->_S576 = _S1241;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1244 = _S1229;
    float _S1245 = _S1229._S576.rows[int(0)].y * _S1229._S576.rows[int(1)].x;
    float det_orig_8 = _S1229._S576.rows[int(0)].x * _S1229._S576.rows[int(1)].y - _S1245;
    float _S1246 = _S1229._S576.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1247 = _S1229._S576;
    *&(((&_S1247)->rows + (int(0)))->x) = _S1246;
    float _S1248 = _S1229._S576.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1247)->rows + (int(1)))->y) = _S1248;
    Matrix<float, 2, 2>  _S1249 = _S1247;
    Matrix<float, 2, 2>  _S1250 = _S1247;
    float det_blur_5 = _S1246 * _S1248 - _S1245;
    float _S1251 = det_orig_8 / det_blur_5;
    float _S1252 = det_blur_5 * det_blur_5;
    float _S1253 = s_primal_ctx_max_0(0.0f, _S1251);
    float _S1254 = s_primal_ctx_sqrt_0(_S1253);
    float invdet_8 = 1.0f / det_blur_5;
    float _S1255 = - _S1229._S576.rows[int(0)].y;
    float _S1256 = - _S1229._S576.rows[int(1)].x;
    float _S1257 = - in_opacity_7;
    float _S1258 = 1.0f + s_primal_ctx_exp_1(_S1257);
    float _S1259 = 1.0f / _S1258;
    float _S1260 = _S1258 * _S1258;
    float _S1261;
    if(antialiased_7)
    {
        _S1261 = _S1259 * _S1254;
    }
    else
    {
        _S1261 = _S1259;
    }
    float _S1262 = _S1261 / 0.00392156885936856f;
    float _S1263 = 2.0f * s_primal_ctx_log_0(_S1262);
    float _S1264 = s_primal_ctx_sqrt_0(_S1263);
    float _S1265 = _S1249.rows[int(0)].x;
    float _S1266 = _S1250.rows[int(1)].y;
    float3  _S1267 = mean_7 - - s_primal_ctx_mul_0(_S1239, t_8);
    float _S1268 = _S1267.x;
    float _S1269 = _S1267.y;
    float _S1270 = _S1267.z;
    float _S1271 = _S1268 * _S1268 + _S1269 * _S1269 + _S1270 * _S1270;
    float _S1272 = s_primal_ctx_sqrt_0(_S1271);
    float x_36 = _S1268 / _S1272;
    float3  _S1273 = make_float3 (x_36);
    float _S1274 = _S1272 * _S1272;
    float y_21 = _S1269 / _S1272;
    float z_18 = _S1270 / _S1272;
    float3  _S1275 = make_float3 (z_18);
    float _S1276 = - y_21;
    float3  _S1277 = make_float3 (_S1276);
    float z2_18 = z_18 * z_18;
    float fTmp0B_7 = -1.09254848957061768f * z_18;
    float fC1_7 = x_36 * x_36 - y_21 * y_21;
    float _S1278 = 2.0f * x_36;
    float fS1_7 = _S1278 * y_21;
    float pSH6_1 = 0.94617468118667603f * z2_18 - 0.31539157032966614f;
    float3  _S1279 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_7 * x_36;
    float3  _S1280 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_7 * y_21;
    float3  _S1281 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_7;
    float3  _S1282 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_7;
    float3  _S1283 = make_float3 (pSH4_1);
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_18;
    float _S1284 = 1.86588168144226074f * z2_18 - 1.11952900886535645f;
    float pSH12_1 = z_18 * _S1284;
    float3  _S1285 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_7 * x_36;
    float3  _S1286 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_7 * y_21;
    float3  _S1287 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_7 * fC1_7;
    float3  _S1288 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_7 * fS1_7;
    float3  _S1289 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_36 * fC1_7 - y_21 * fS1_7);
    float3  _S1290 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_36 * fS1_7 + y_21 * fC1_7);
    float3  _S1291 = make_float3 (pSH9_1);
    float3  _S1292 = make_float3 (0.0f);
    float3  _S1293 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1294;
    (&_S1294)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1276) * (*sh_coeffs_7)[int(1)] + make_float3 (z_18) * (*sh_coeffs_7)[int(2)] - make_float3 (x_36) * (*sh_coeffs_7)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_7)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_7)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_7)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_7)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_7)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_7)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_7)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_7)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_7)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_7)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_7)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f);
    (&_S1294)->differential_0 = _S1293;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1295;
    (&_S1295)->primal_0 = _S1292;
    (&_S1295)->differential_0 = _S1293;
    s_bwd_prop_max_0(&_S1294, &_S1295, v_rgb_1);
    float3  _S1296 = _S1290 * _S1294.differential_0;
    float3  _S1297 = (*sh_coeffs_7)[int(15)] * _S1294.differential_0;
    float3  _S1298 = _S1288 * _S1294.differential_0;
    float3  _S1299 = (*sh_coeffs_7)[int(14)] * _S1294.differential_0;
    float3  _S1300 = _S1286 * _S1294.differential_0;
    float3  _S1301 = (*sh_coeffs_7)[int(13)] * _S1294.differential_0;
    float3  _S1302 = _S1285 * _S1294.differential_0;
    float3  _S1303 = (*sh_coeffs_7)[int(12)] * _S1294.differential_0;
    float3  _S1304 = _S1287 * _S1294.differential_0;
    float3  _S1305 = (*sh_coeffs_7)[int(11)] * _S1294.differential_0;
    float3  _S1306 = _S1289 * _S1294.differential_0;
    float3  _S1307 = (*sh_coeffs_7)[int(10)] * _S1294.differential_0;
    float3  _S1308 = _S1291 * _S1294.differential_0;
    float3  _S1309 = (*sh_coeffs_7)[int(9)] * _S1294.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1309.x + _S1309.y + _S1309.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1297.x + _S1297.y + _S1297.z);
    float _S1310 = _S1307.x + _S1307.y + _S1307.z;
    float _S1311 = _S1299.x + _S1299.y + _S1299.z;
    float _S1312 = _S1305.x + _S1305.y + _S1305.z;
    float _S1313 = _S1301.x + _S1301.y + _S1301.z;
    float _S1314 = _S1303.x + _S1303.y + _S1303.z;
    float _S1315 = - s_diff_fC2_T_1;
    float3  _S1316 = _S1282 * _S1294.differential_0;
    float3  _S1317 = (*sh_coeffs_7)[int(8)] * _S1294.differential_0;
    float3  _S1318 = _S1280 * _S1294.differential_0;
    float3  _S1319 = (*sh_coeffs_7)[int(7)] * _S1294.differential_0;
    float3  _S1320 = _S1279 * _S1294.differential_0;
    float3  _S1321 = (*sh_coeffs_7)[int(6)] * _S1294.differential_0;
    float3  _S1322 = _S1281 * _S1294.differential_0;
    float3  _S1323 = (*sh_coeffs_7)[int(5)] * _S1294.differential_0;
    float3  _S1324 = _S1283 * _S1294.differential_0;
    float3  _S1325 = (*sh_coeffs_7)[int(4)] * _S1294.differential_0;
    float _S1326 = _S1323.x + _S1323.y + _S1323.z;
    float _S1327 = _S1319.x + _S1319.y + _S1319.z;
    float _S1328 = fTmp1B_7 * _S1310 + x_36 * s_diff_fS2_T_1 + y_21 * _S1315 + 0.54627424478530884f * (_S1325.x + _S1325.y + _S1325.z);
    float _S1329 = fTmp1B_7 * _S1311 + y_21 * s_diff_fS2_T_1 + x_36 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1317.x + _S1317.y + _S1317.z);
    float _S1330 = y_21 * - _S1329;
    float _S1331 = x_36 * _S1329;
    float _S1332 = z_18 * (1.86588168144226074f * (z_18 * _S1314) + -2.28522896766662598f * (y_21 * _S1312 + x_36 * _S1313) + 0.94617468118667603f * (_S1321.x + _S1321.y + _S1321.z));
    float3  _S1333 = make_float3 (0.48860251903533936f) * _S1294.differential_0;
    float3  _S1334 = - _S1333;
    float3  _S1335 = _S1273 * _S1334;
    float3  _S1336 = (*sh_coeffs_7)[int(3)] * _S1334;
    float3  _S1337 = _S1275 * _S1333;
    float3  _S1338 = (*sh_coeffs_7)[int(2)] * _S1333;
    float3  _S1339 = _S1277 * _S1333;
    float3  _S1340 = (*sh_coeffs_7)[int(1)] * _S1333;
    float _S1341 = (_S1284 * _S1314 + 1.44530570507049561f * (fS1_7 * _S1310 + fC1_7 * _S1311) + -1.09254848957061768f * (y_21 * _S1326 + x_36 * _S1327) + _S1332 + _S1332 + _S1338.x + _S1338.y + _S1338.z) / _S1274;
    float _S1342 = _S1272 * _S1341;
    float _S1343 = (fTmp0C_7 * _S1312 + fC1_7 * s_diff_fS2_T_1 + fS1_7 * _S1315 + fTmp0B_7 * _S1326 + _S1278 * _S1328 + _S1330 + _S1330 + - (_S1340.x + _S1340.y + _S1340.z)) / _S1274;
    float _S1344 = _S1272 * _S1343;
    float _S1345 = (fTmp0C_7 * _S1313 + fS1_7 * s_diff_fS2_T_1 + fC1_7 * s_diff_fC2_T_1 + fTmp0B_7 * _S1327 + 2.0f * (y_21 * _S1328) + _S1331 + _S1331 + _S1336.x + _S1336.y + _S1336.z) / _S1274;
    float _S1346 = _S1272 * _S1345;
    float _S1347 = _S1270 * - _S1341 + _S1269 * - _S1343 + _S1268 * - _S1345;
    DiffPair_float_0 _S1348;
    (&_S1348)->primal_0 = _S1271;
    (&_S1348)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1348, _S1347);
    float _S1349 = _S1270 * _S1348.differential_0;
    float _S1350 = _S1269 * _S1348.differential_0;
    float _S1351 = _S1268 * _S1348.differential_0;
    float3  _S1352 = make_float3 (0.282094806432724f) * _S1294.differential_0;
    float3  _S1353 = make_float3 (_S1346 + _S1351 + _S1351, _S1344 + _S1350 + _S1350, _S1342 + _S1349 + _S1349);
    float3  _S1354 = - - _S1353;
    Matrix<float, 3, 3>  _S1355 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1356;
    (&_S1356)->primal_0 = _S1239;
    (&_S1356)->differential_0 = _S1355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1357;
    (&_S1357)->primal_0 = t_8;
    (&_S1357)->differential_0 = _S1293;
    s_bwd_prop_mul_0(&_S1356, &_S1357, _S1354);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1358 = _S1356;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1359 = _S1357;
    float2  _S1360 = _S1242;
    *&((&_S1360)->y) = v_conic_1.z;
    float2  _S1361 = _S1242;
    *&((&_S1361)->y) = v_conic_1.y;
    *&((&_S1361)->x) = v_conic_1.x;
    DiffPair_float_0 _S1362;
    (&_S1362)->primal_0 = _S1266;
    (&_S1362)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1362, 0.0f);
    DiffPair_float_0 _S1363;
    (&_S1363)->primal_0 = _S1265;
    (&_S1363)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1363, 0.0f);
    DiffPair_float_0 _S1364;
    (&_S1364)->primal_0 = 3.32999992370605469f;
    (&_S1364)->differential_0 = 0.0f;
    DiffPair_float_0 _S1365;
    (&_S1365)->primal_0 = _S1264;
    (&_S1365)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1364, &_S1365, 0.0f);
    DiffPair_float_0 _S1366;
    (&_S1366)->primal_0 = _S1263;
    (&_S1366)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1366, _S1365.differential_0);
    float _S1367 = 2.0f * _S1366.differential_0;
    DiffPair_float_0 _S1368;
    (&_S1368)->primal_0 = _S1262;
    (&_S1368)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1368, _S1367);
    float2  _S1369 = make_float2 (_S1363.differential_0, 0.0f);
    float _S1370 = v_opacity_1 + 254.9999847412109375f * _S1368.differential_0;
    Matrix<float, 2, 2>  _S1371 = _S1225;
    _S1371[int(1)] = _S1360;
    _S1371[int(0)] = _S1361;
    Matrix<float, 2, 2>  _S1372 = _S1371;
    FixedArray<float3 , 16>  _S1373;
    _S1373[int(0)] = _S1293;
    _S1373[int(1)] = _S1293;
    _S1373[int(2)] = _S1293;
    _S1373[int(3)] = _S1293;
    _S1373[int(4)] = _S1293;
    _S1373[int(5)] = _S1293;
    _S1373[int(6)] = _S1293;
    _S1373[int(7)] = _S1293;
    _S1373[int(8)] = _S1293;
    _S1373[int(9)] = _S1293;
    _S1373[int(10)] = _S1293;
    _S1373[int(11)] = _S1293;
    _S1373[int(12)] = _S1293;
    _S1373[int(13)] = _S1293;
    _S1373[int(14)] = _S1293;
    _S1373[int(15)] = _S1293;
    _S1373[int(7)] = _S1318;
    _S1373[int(0)] = _S1352;
    _S1373[int(1)] = _S1339;
    _S1373[int(2)] = _S1337;
    _S1373[int(3)] = _S1335;
    _S1373[int(4)] = _S1324;
    _S1373[int(5)] = _S1322;
    _S1373[int(6)] = _S1320;
    _S1373[int(15)] = _S1296;
    _S1373[int(8)] = _S1316;
    _S1373[int(9)] = _S1308;
    _S1373[int(10)] = _S1306;
    _S1373[int(11)] = _S1304;
    _S1373[int(12)] = _S1302;
    _S1373[int(13)] = _S1300;
    _S1373[int(14)] = _S1298;
    float3  _S1374 = _S1373[int(0)];
    float3  _S1375 = _S1373[int(1)];
    float3  _S1376 = _S1373[int(2)];
    float3  _S1377 = _S1373[int(3)];
    float3  _S1378 = _S1373[int(4)];
    float3  _S1379 = _S1373[int(5)];
    float3  _S1380 = _S1373[int(6)];
    float3  _S1381 = _S1373[int(7)];
    float3  _S1382 = _S1373[int(8)];
    float3  _S1383 = _S1373[int(9)];
    float3  _S1384 = _S1373[int(10)];
    float3  _S1385 = _S1373[int(11)];
    float3  _S1386 = _S1373[int(12)];
    float3  _S1387 = _S1373[int(13)];
    float3  _S1388 = _S1373[int(14)];
    float3  _S1389 = _S1373[int(15)];
    float3  _S1390 = make_float3 (0.0f, 0.0f, v_depth_1);
    float2  _S1391 = make_float2 (0.0f, _S1362.differential_0);
    float _S1392;
    if(antialiased_7)
    {
        float _S1393 = _S1259 * _S1370;
        _S1261 = _S1254 * _S1370;
        _S1392 = _S1393;
    }
    else
    {
        _S1261 = _S1370;
        _S1392 = 0.0f;
    }
    float _S1394 = - (_S1261 / _S1260);
    DiffPair_float_0 _S1395;
    (&_S1395)->primal_0 = _S1257;
    (&_S1395)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1395, _S1394);
    float _S1396 = - _S1395.differential_0;
    float _S1397 = invdet_8 * _S1372.rows[int(1)].y;
    float _S1398 = - (invdet_8 * _S1372.rows[int(1)].x);
    float _S1399 = - (invdet_8 * _S1372.rows[int(0)].y);
    float _S1400 = invdet_8 * _S1372.rows[int(0)].x;
    float _S1401 = - ((_S1246 * _S1372.rows[int(1)].y + _S1256 * _S1372.rows[int(1)].x + _S1255 * _S1372.rows[int(0)].y + _S1248 * _S1372.rows[int(0)].x) / _S1252);
    DiffPair_float_0 _S1402;
    (&_S1402)->primal_0 = _S1253;
    (&_S1402)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1402, _S1392);
    DiffPair_float_0 _S1403;
    (&_S1403)->primal_0 = 0.0f;
    (&_S1403)->differential_0 = 0.0f;
    DiffPair_float_0 _S1404;
    (&_S1404)->primal_0 = _S1251;
    (&_S1404)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1403, &_S1404, _S1402.differential_0);
    float _S1405 = _S1404.differential_0 / _S1252;
    float s_diff_det_orig_T_1 = det_blur_5 * _S1405;
    float _S1406 = _S1401 + det_orig_8 * - _S1405;
    float _S1407 = - _S1406;
    float _S1408 = _S1246 * _S1406;
    float _S1409 = _S1248 * _S1406;
    Matrix<float, 2, 2>  _S1410 = _S1225;
    _S1410[int(1)] = _S1391;
    _S1410[int(0)] = _S1369;
    _S1247 = _S1410;
    *&(((&_S1247)->rows + (int(1)))->y) = 0.0f;
    float _S1411 = _S1400 + _S1408 + _S1410.rows[int(1)].y;
    *&(((&_S1247)->rows + (int(0)))->x) = 0.0f;
    float _S1412 = _S1397 + _S1409 + _S1410.rows[int(0)].x;
    float _S1413 = _S1407 + - s_diff_det_orig_T_1;
    float _S1414 = _S1398 + _S1244._S576.rows[int(0)].y * _S1413;
    float _S1415 = _S1399 + _S1244._S576.rows[int(1)].x * _S1413;
    float _S1416 = _S1244._S576.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S1417 = _S1411 + _S1244._S576.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S1418 = _S1242;
    *&((&_S1418)->x) = _S1414;
    *&((&_S1418)->y) = _S1417;
    float _S1419 = _S1412 + _S1416;
    float2  _S1420 = _S1242;
    *&((&_S1420)->y) = _S1415;
    *&((&_S1420)->x) = _S1419;
    Matrix<float, 2, 2>  _S1421 = _S1225;
    _S1421[int(1)] = _S1418;
    _S1421[int(0)] = _S1420;
    Matrix<float, 2, 2>  _S1422 = _S1247 + _S1421;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1423;
    (&_S1423)->primal_0 = mean_c_7;
    (&_S1423)->differential_0 = _S1293;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1424;
    (&_S1424)->primal_0 = _S1240;
    (&_S1424)->differential_0 = _S1355;
    DiffPair_float_0 _S1425;
    (&_S1425)->primal_0 = fx_11;
    (&_S1425)->differential_0 = 0.0f;
    DiffPair_float_0 _S1426;
    (&_S1426)->primal_0 = fy_11;
    (&_S1426)->differential_0 = 0.0f;
    DiffPair_float_0 _S1427;
    (&_S1427)->primal_0 = cx_11;
    (&_S1427)->differential_0 = 0.0f;
    DiffPair_float_0 _S1428;
    (&_S1428)->primal_0 = cy_11;
    (&_S1428)->differential_0 = 0.0f;
    float4  _S1429 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1430;
    (&_S1430)->primal_0 = radial_coeffs_10;
    (&_S1430)->differential_0 = _S1429;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1431;
    (&_S1431)->primal_0 = tangential_coeffs_10;
    (&_S1431)->differential_0 = _S1242;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1432;
    (&_S1432)->primal_0 = thin_prism_coeffs_10;
    (&_S1432)->differential_0 = _S1242;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1433 = _S1244._S577;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S1423, &_S1424, &_S1425, &_S1426, &_S1427, &_S1428, &_S1430, &_S1431, &_S1432, _S1422, v_mean2d_1, &_S1433);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1434;
    (&_S1434)->primal_0 = _S1238;
    (&_S1434)->differential_0 = _S1355;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1435;
    (&_S1435)->primal_0 = _S1239;
    (&_S1435)->differential_0 = _S1355;
    s_bwd_prop_mul_3(&_S1434, &_S1435, _S1424.differential_0);
    Matrix<float, 3, 3>  _S1436 = transpose_0(_S1435.differential_0 + _S1358.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1437;
    (&_S1437)->primal_0 = R_9;
    (&_S1437)->differential_0 = _S1355;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1438;
    (&_S1438)->primal_0 = _S1237;
    (&_S1438)->differential_0 = _S1355;
    s_bwd_prop_mul_3(&_S1437, &_S1438, _S1434.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1439;
    (&_S1439)->primal_0 = _S1235;
    (&_S1439)->differential_0 = _S1355;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1440;
    (&_S1440)->primal_0 = _S1236;
    (&_S1440)->differential_0 = _S1355;
    s_bwd_prop_mul_3(&_S1439, &_S1440, _S1438.differential_0);
    Matrix<float, 3, 3>  _S1441 = _S1439.differential_0 + transpose_0(_S1440.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1442;
    (&_S1442)->primal_0 = _S1234;
    (&_S1442)->differential_0 = _S1355;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1443;
    (&_S1443)->primal_0 = S_1;
    (&_S1443)->differential_0 = _S1355;
    s_bwd_prop_mul_3(&_S1442, &_S1443, _S1441);
    Matrix<float, 3, 3>  _S1444 = transpose_0(_S1442.differential_0);
    float _S1445 = 2.0f * - _S1444.rows[int(2)].z;
    float _S1446 = 2.0f * _S1444.rows[int(2)].y;
    float _S1447 = 2.0f * _S1444.rows[int(2)].x;
    float _S1448 = 2.0f * _S1444.rows[int(1)].z;
    float _S1449 = 2.0f * - _S1444.rows[int(1)].y;
    float _S1450 = 2.0f * _S1444.rows[int(1)].x;
    float _S1451 = 2.0f * _S1444.rows[int(0)].z;
    float _S1452 = 2.0f * _S1444.rows[int(0)].y;
    float _S1453 = 2.0f * - _S1444.rows[int(0)].x;
    float _S1454 = - _S1450 + _S1452;
    float _S1455 = _S1447 + - _S1451;
    float _S1456 = - _S1446 + _S1448;
    float _S1457 = _S1446 + _S1448;
    float _S1458 = _S1447 + _S1451;
    float _S1459 = _S1450 + _S1452;
    float _S1460 = z_17 * (_S1449 + _S1453);
    float _S1461 = y_20 * (_S1445 + _S1453);
    float _S1462 = x_35 * (_S1445 + _S1449);
    float _S1463 = z_17 * _S1454 + y_20 * _S1455 + x_35 * _S1456;
    float _S1464 = _S1233 * _S1463;
    float _S1465 = w_10 * _S1454 + y_20 * _S1457 + x_35 * _S1458 + _S1460 + _S1460;
    float _S1466 = _S1233 * _S1465;
    float _S1467 = w_10 * _S1455 + z_17 * _S1457 + x_35 * _S1459 + _S1461 + _S1461;
    float _S1468 = _S1233 * _S1467;
    float _S1469 = w_10 * _S1456 + z_17 * _S1458 + y_20 * _S1459 + _S1462 + _S1462;
    float _S1470 = _S1233 * _S1469;
    float _S1471 = quat_10.x * _S1463 + quat_10.w * _S1465 + quat_10.z * _S1467 + quat_10.y * _S1469;
    DiffPair_float_0 _S1472;
    (&_S1472)->primal_0 = _S1232;
    (&_S1472)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1472, _S1471);
    float _S1473 = quat_10.x * _S1472.differential_0;
    float _S1474 = quat_10.w * _S1472.differential_0;
    float _S1475 = quat_10.z * _S1472.differential_0;
    float _S1476 = quat_10.y * _S1472.differential_0;
    float _S1477 = _S1466 + _S1474 + _S1474;
    float _S1478 = _S1468 + _S1475 + _S1475;
    float _S1479 = _S1470 + _S1476 + _S1476;
    float _S1480 = _S1464 + _S1473 + _S1473;
    float3  _S1481 = _S1293;
    *&((&_S1481)->z) = _S1443.differential_0.rows[int(2)].z;
    *&((&_S1481)->y) = _S1443.differential_0.rows[int(1)].y;
    *&((&_S1481)->x) = _S1443.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1482;
    (&_S1482)->primal_0 = scale_9;
    (&_S1482)->differential_0 = _S1293;
    s_bwd_prop_exp_1(&_S1482, _S1481);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1483 = _S1482;
    float3  _S1484 = _S1423.differential_0 + _S1390;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1485;
    (&_S1485)->primal_0 = R_9;
    (&_S1485)->differential_0 = _S1355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1486;
    (&_S1486)->primal_0 = mean_7;
    (&_S1486)->differential_0 = _S1293;
    s_bwd_prop_mul_0(&_S1485, &_S1486, _S1484);
    float3  _S1487 = _S1484 + _S1359.differential_0;
    Matrix<float, 3, 3>  _S1488 = _S1436 + _S1437.differential_0 + _S1485.differential_0;
    float4  _S1489 = _S1429;
    *&((&_S1489)->w) = _S1477;
    *&((&_S1489)->z) = _S1478;
    *&((&_S1489)->y) = _S1479;
    *&((&_S1489)->x) = _S1480;
    float4  _S1490 = _S1489;
    float3  _S1491 = _S1486.differential_0 + _S1353;
    *v_mean_1 = _S1491;
    *v_quat_1 = _S1490;
    *v_scale_1 = _S1483.differential_0;
    *v_in_opacity_1 = _S1396;
    (*v_sh_coeffs_1)[int(0)] = _S1374;
    (*v_sh_coeffs_1)[int(1)] = _S1375;
    (*v_sh_coeffs_1)[int(2)] = _S1376;
    (*v_sh_coeffs_1)[int(3)] = _S1377;
    (*v_sh_coeffs_1)[int(4)] = _S1378;
    (*v_sh_coeffs_1)[int(5)] = _S1379;
    (*v_sh_coeffs_1)[int(6)] = _S1380;
    (*v_sh_coeffs_1)[int(7)] = _S1381;
    (*v_sh_coeffs_1)[int(8)] = _S1382;
    (*v_sh_coeffs_1)[int(9)] = _S1383;
    (*v_sh_coeffs_1)[int(10)] = _S1384;
    (*v_sh_coeffs_1)[int(11)] = _S1385;
    (*v_sh_coeffs_1)[int(12)] = _S1386;
    (*v_sh_coeffs_1)[int(13)] = _S1387;
    (*v_sh_coeffs_1)[int(14)] = _S1388;
    (*v_sh_coeffs_1)[int(15)] = _S1389;
    *v_R_1 = _S1488;
    *v_t_1 = _S1487;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_12, float fy_12, float cx_12, float cy_12, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, uint image_width_8, uint image_height_8, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    float3  _S1492 = s_primal_ctx_exp_0(scale_10);
    float _S1493 = quat_11.y;
    float _S1494 = _S1493 * _S1493 + quat_11.z * quat_11.z + quat_11.w * quat_11.w + quat_11.x * quat_11.x;
    float _S1495 = s_primal_ctx_rsqrt_0(_S1494);
    float x_37 = quat_11.y * _S1495;
    float y_22 = quat_11.z * _S1495;
    float z_19 = quat_11.w * _S1495;
    float w_11 = quat_11.x * _S1495;
    float x2_11 = x_37 * x_37;
    float y2_11 = y_22 * y_22;
    float z2_19 = z_19 * z_19;
    float xy_11 = x_37 * y_22;
    float xz_11 = x_37 * z_19;
    float yz_11 = y_22 * z_19;
    float wx_11 = w_11 * x_37;
    float wy_11 = w_11 * y_22;
    float wz_11 = w_11 * z_19;
    Matrix<float, 3, 3>  _S1496 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S1492.x, 0.0f, 0.0f, 0.0f, _S1492.y, 0.0f, 0.0f, 0.0f, _S1492.z);
    Matrix<float, 3, 3>  _S1497 = s_primal_ctx_mul_1(_S1496, S_2);
    Matrix<float, 3, 3>  _S1498 = transpose_0(_S1497);
    Matrix<float, 3, 3>  _S1499 = s_primal_ctx_mul_1(_S1497, _S1498);
    Matrix<float, 3, 3>  _S1500 = s_primal_ctx_mul_1(R_10, _S1499);
    Matrix<float, 3, 3>  _S1501 = transpose_0(R_10);
    Matrix<float, 3, 3>  _S1502 = s_primal_ctx_mul_1(_S1500, _S1501);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_12, 0.0f, 0.0f, 0.0f, fy_12, 0.0f);
    Matrix<float, 2, 3>  _S1503 = s_primal_ctx_mul_2(J_11, _S1502);
    Matrix<float, 3, 2>  _S1504 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S1505 = s_primal_ctx_mul_3(_S1503, _S1504);
    float _S1506 = _S1505.rows[int(0)].y * _S1505.rows[int(1)].x;
    float det_orig_9 = _S1505.rows[int(0)].x * _S1505.rows[int(1)].y - _S1506;
    float _S1507 = _S1505.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1508 = _S1505;
    *&(((&_S1508)->rows + (int(0)))->x) = _S1507;
    float _S1509 = _S1505.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1508)->rows + (int(1)))->y) = _S1509;
    Matrix<float, 2, 2>  _S1510 = _S1508;
    Matrix<float, 2, 2>  _S1511 = _S1508;
    float det_blur_6 = _S1507 * _S1509 - _S1506;
    float _S1512 = det_orig_9 / det_blur_6;
    float _S1513 = det_blur_6 * det_blur_6;
    float _S1514 = s_primal_ctx_max_0(0.0f, _S1512);
    float _S1515 = s_primal_ctx_sqrt_0(_S1514);
    float invdet_9 = 1.0f / det_blur_6;
    float _S1516 = - _S1505.rows[int(0)].y;
    float _S1517 = - _S1505.rows[int(1)].x;
    float _S1518 = - in_opacity_8;
    float _S1519 = 1.0f + s_primal_ctx_exp_1(_S1518);
    float _S1520 = 1.0f / _S1519;
    float _S1521 = _S1519 * _S1519;
    float _S1522;
    if(antialiased_8)
    {
        _S1522 = _S1520 * _S1515;
    }
    else
    {
        _S1522 = _S1520;
    }
    float _S1523 = _S1522 / 0.00392156885936856f;
    float _S1524 = 2.0f * s_primal_ctx_log_0(_S1523);
    float _S1525 = s_primal_ctx_sqrt_0(_S1524);
    float _S1526 = _S1510.rows[int(0)].x;
    float _S1527 = _S1511.rows[int(1)].y;
    float3  _S1528 = mean_8 - - s_primal_ctx_mul_0(_S1501, t_9);
    float _S1529 = _S1528.x;
    float _S1530 = _S1528.y;
    float _S1531 = _S1528.z;
    float _S1532 = _S1529 * _S1529 + _S1530 * _S1530 + _S1531 * _S1531;
    float _S1533 = s_primal_ctx_sqrt_0(_S1532);
    float x_38 = _S1529 / _S1533;
    float3  _S1534 = make_float3 (x_38);
    float _S1535 = _S1533 * _S1533;
    float y_23 = _S1530 / _S1533;
    float z_20 = _S1531 / _S1533;
    float3  _S1536 = make_float3 (z_20);
    float _S1537 = - y_23;
    float3  _S1538 = make_float3 (_S1537);
    float z2_20 = z_20 * z_20;
    float fTmp0B_8 = -1.09254848957061768f * z_20;
    float fC1_8 = x_38 * x_38 - y_23 * y_23;
    float _S1539 = 2.0f * x_38;
    float fS1_8 = _S1539 * y_23;
    float pSH6_2 = 0.94617468118667603f * z2_20 - 0.31539157032966614f;
    float3  _S1540 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_8 * x_38;
    float3  _S1541 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_8 * y_23;
    float3  _S1542 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_8;
    float3  _S1543 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_8;
    float3  _S1544 = make_float3 (pSH4_2);
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_20;
    float _S1545 = 1.86588168144226074f * z2_20 - 1.11952900886535645f;
    float pSH12_2 = z_20 * _S1545;
    float3  _S1546 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_8 * x_38;
    float3  _S1547 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_8 * y_23;
    float3  _S1548 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_8 * fC1_8;
    float3  _S1549 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_8 * fS1_8;
    float3  _S1550 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_38 * fC1_8 - y_23 * fS1_8);
    float3  _S1551 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_38 * fS1_8 + y_23 * fC1_8);
    float3  _S1552 = make_float3 (pSH9_2);
    float3  _S1553 = make_float3 (0.0f);
    float3  _S1554 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1555;
    (&_S1555)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1537) * (*sh_coeffs_8)[int(1)] + make_float3 (z_20) * (*sh_coeffs_8)[int(2)] - make_float3 (x_38) * (*sh_coeffs_8)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_8)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_8)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_8)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_8)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_8)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_8)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_8)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_8)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_8)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_8)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_8)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f);
    (&_S1555)->differential_0 = _S1554;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1556;
    (&_S1556)->primal_0 = _S1553;
    (&_S1556)->differential_0 = _S1554;
    s_bwd_prop_max_0(&_S1555, &_S1556, v_rgb_2);
    float3  _S1557 = _S1551 * _S1555.differential_0;
    float3  _S1558 = (*sh_coeffs_8)[int(15)] * _S1555.differential_0;
    float3  _S1559 = _S1549 * _S1555.differential_0;
    float3  _S1560 = (*sh_coeffs_8)[int(14)] * _S1555.differential_0;
    float3  _S1561 = _S1547 * _S1555.differential_0;
    float3  _S1562 = (*sh_coeffs_8)[int(13)] * _S1555.differential_0;
    float3  _S1563 = _S1546 * _S1555.differential_0;
    float3  _S1564 = (*sh_coeffs_8)[int(12)] * _S1555.differential_0;
    float3  _S1565 = _S1548 * _S1555.differential_0;
    float3  _S1566 = (*sh_coeffs_8)[int(11)] * _S1555.differential_0;
    float3  _S1567 = _S1550 * _S1555.differential_0;
    float3  _S1568 = (*sh_coeffs_8)[int(10)] * _S1555.differential_0;
    float3  _S1569 = _S1552 * _S1555.differential_0;
    float3  _S1570 = (*sh_coeffs_8)[int(9)] * _S1555.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S1570.x + _S1570.y + _S1570.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S1558.x + _S1558.y + _S1558.z);
    float _S1571 = _S1568.x + _S1568.y + _S1568.z;
    float _S1572 = _S1560.x + _S1560.y + _S1560.z;
    float _S1573 = _S1566.x + _S1566.y + _S1566.z;
    float _S1574 = _S1562.x + _S1562.y + _S1562.z;
    float _S1575 = _S1564.x + _S1564.y + _S1564.z;
    float _S1576 = - s_diff_fC2_T_2;
    float3  _S1577 = _S1543 * _S1555.differential_0;
    float3  _S1578 = (*sh_coeffs_8)[int(8)] * _S1555.differential_0;
    float3  _S1579 = _S1541 * _S1555.differential_0;
    float3  _S1580 = (*sh_coeffs_8)[int(7)] * _S1555.differential_0;
    float3  _S1581 = _S1540 * _S1555.differential_0;
    float3  _S1582 = (*sh_coeffs_8)[int(6)] * _S1555.differential_0;
    float3  _S1583 = _S1542 * _S1555.differential_0;
    float3  _S1584 = (*sh_coeffs_8)[int(5)] * _S1555.differential_0;
    float3  _S1585 = _S1544 * _S1555.differential_0;
    float3  _S1586 = (*sh_coeffs_8)[int(4)] * _S1555.differential_0;
    float _S1587 = _S1584.x + _S1584.y + _S1584.z;
    float _S1588 = _S1580.x + _S1580.y + _S1580.z;
    float _S1589 = fTmp1B_8 * _S1571 + x_38 * s_diff_fS2_T_2 + y_23 * _S1576 + 0.54627424478530884f * (_S1586.x + _S1586.y + _S1586.z);
    float _S1590 = fTmp1B_8 * _S1572 + y_23 * s_diff_fS2_T_2 + x_38 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S1578.x + _S1578.y + _S1578.z);
    float _S1591 = y_23 * - _S1590;
    float _S1592 = x_38 * _S1590;
    float _S1593 = z_20 * (1.86588168144226074f * (z_20 * _S1575) + -2.28522896766662598f * (y_23 * _S1573 + x_38 * _S1574) + 0.94617468118667603f * (_S1582.x + _S1582.y + _S1582.z));
    float3  _S1594 = make_float3 (0.48860251903533936f) * _S1555.differential_0;
    float3  _S1595 = - _S1594;
    float3  _S1596 = _S1534 * _S1595;
    float3  _S1597 = (*sh_coeffs_8)[int(3)] * _S1595;
    float3  _S1598 = _S1536 * _S1594;
    float3  _S1599 = (*sh_coeffs_8)[int(2)] * _S1594;
    float3  _S1600 = _S1538 * _S1594;
    float3  _S1601 = (*sh_coeffs_8)[int(1)] * _S1594;
    float _S1602 = (_S1545 * _S1575 + 1.44530570507049561f * (fS1_8 * _S1571 + fC1_8 * _S1572) + -1.09254848957061768f * (y_23 * _S1587 + x_38 * _S1588) + _S1593 + _S1593 + _S1599.x + _S1599.y + _S1599.z) / _S1535;
    float _S1603 = _S1533 * _S1602;
    float _S1604 = (fTmp0C_8 * _S1573 + fC1_8 * s_diff_fS2_T_2 + fS1_8 * _S1576 + fTmp0B_8 * _S1587 + _S1539 * _S1589 + _S1591 + _S1591 + - (_S1601.x + _S1601.y + _S1601.z)) / _S1535;
    float _S1605 = _S1533 * _S1604;
    float _S1606 = (fTmp0C_8 * _S1574 + fS1_8 * s_diff_fS2_T_2 + fC1_8 * s_diff_fC2_T_2 + fTmp0B_8 * _S1588 + 2.0f * (y_23 * _S1589) + _S1592 + _S1592 + _S1597.x + _S1597.y + _S1597.z) / _S1535;
    float _S1607 = _S1533 * _S1606;
    float _S1608 = _S1531 * - _S1602 + _S1530 * - _S1604 + _S1529 * - _S1606;
    DiffPair_float_0 _S1609;
    (&_S1609)->primal_0 = _S1532;
    (&_S1609)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1609, _S1608);
    float _S1610 = _S1531 * _S1609.differential_0;
    float _S1611 = _S1530 * _S1609.differential_0;
    float _S1612 = _S1529 * _S1609.differential_0;
    float3  _S1613 = make_float3 (0.282094806432724f) * _S1555.differential_0;
    float3  _S1614 = make_float3 (_S1607 + _S1612 + _S1612, _S1605 + _S1611 + _S1611, _S1603 + _S1610 + _S1610);
    float3  _S1615 = - - _S1614;
    Matrix<float, 3, 3>  _S1616 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1617;
    (&_S1617)->primal_0 = _S1501;
    (&_S1617)->differential_0 = _S1616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1618;
    (&_S1618)->primal_0 = t_9;
    (&_S1618)->differential_0 = _S1554;
    s_bwd_prop_mul_0(&_S1617, &_S1618, _S1615);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1619 = _S1617;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1620 = _S1618;
    float2  _S1621 = make_float2 (0.0f);
    float2  _S1622 = _S1621;
    *&((&_S1622)->y) = v_conic_2.z;
    float2  _S1623 = _S1621;
    *&((&_S1623)->y) = v_conic_2.y;
    *&((&_S1623)->x) = v_conic_2.x;
    DiffPair_float_0 _S1624;
    (&_S1624)->primal_0 = _S1527;
    (&_S1624)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1624, 0.0f);
    DiffPair_float_0 _S1625;
    (&_S1625)->primal_0 = _S1526;
    (&_S1625)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1625, 0.0f);
    DiffPair_float_0 _S1626;
    (&_S1626)->primal_0 = 3.32999992370605469f;
    (&_S1626)->differential_0 = 0.0f;
    DiffPair_float_0 _S1627;
    (&_S1627)->primal_0 = _S1525;
    (&_S1627)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1626, &_S1627, 0.0f);
    DiffPair_float_0 _S1628;
    (&_S1628)->primal_0 = _S1524;
    (&_S1628)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1628, _S1627.differential_0);
    float _S1629 = 2.0f * _S1628.differential_0;
    DiffPair_float_0 _S1630;
    (&_S1630)->primal_0 = _S1523;
    (&_S1630)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1630, _S1629);
    float _S1631 = v_opacity_2 + 254.9999847412109375f * _S1630.differential_0;
    Matrix<float, 2, 2>  _S1632 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1633 = _S1632;
    _S1633[int(1)] = _S1622;
    _S1633[int(0)] = _S1623;
    Matrix<float, 2, 2>  _S1634 = _S1633;
    FixedArray<float3 , 16>  _S1635;
    _S1635[int(0)] = _S1554;
    _S1635[int(1)] = _S1554;
    _S1635[int(2)] = _S1554;
    _S1635[int(3)] = _S1554;
    _S1635[int(4)] = _S1554;
    _S1635[int(5)] = _S1554;
    _S1635[int(6)] = _S1554;
    _S1635[int(7)] = _S1554;
    _S1635[int(8)] = _S1554;
    _S1635[int(9)] = _S1554;
    _S1635[int(10)] = _S1554;
    _S1635[int(11)] = _S1554;
    _S1635[int(12)] = _S1554;
    _S1635[int(13)] = _S1554;
    _S1635[int(14)] = _S1554;
    _S1635[int(15)] = _S1554;
    _S1635[int(7)] = _S1579;
    _S1635[int(0)] = _S1613;
    _S1635[int(1)] = _S1600;
    _S1635[int(2)] = _S1598;
    _S1635[int(3)] = _S1596;
    _S1635[int(4)] = _S1585;
    _S1635[int(5)] = _S1583;
    _S1635[int(6)] = _S1581;
    _S1635[int(15)] = _S1557;
    _S1635[int(8)] = _S1577;
    _S1635[int(9)] = _S1569;
    _S1635[int(10)] = _S1567;
    _S1635[int(11)] = _S1565;
    _S1635[int(12)] = _S1563;
    _S1635[int(13)] = _S1561;
    _S1635[int(14)] = _S1559;
    float3  _S1636 = _S1635[int(0)];
    float3  _S1637 = _S1635[int(1)];
    float3  _S1638 = _S1635[int(2)];
    float3  _S1639 = _S1635[int(3)];
    float3  _S1640 = _S1635[int(4)];
    float3  _S1641 = _S1635[int(5)];
    float3  _S1642 = _S1635[int(6)];
    float3  _S1643 = _S1635[int(7)];
    float3  _S1644 = _S1635[int(8)];
    float3  _S1645 = _S1635[int(9)];
    float3  _S1646 = _S1635[int(10)];
    float3  _S1647 = _S1635[int(11)];
    float3  _S1648 = _S1635[int(12)];
    float3  _S1649 = _S1635[int(13)];
    float3  _S1650 = _S1635[int(14)];
    float3  _S1651 = _S1635[int(15)];
    float3  _S1652 = make_float3 (0.0f, 0.0f, v_depth_2);
    float2  _S1653 = make_float2 (0.0f, _S1624.differential_0);
    float2  _S1654 = make_float2 (_S1625.differential_0, 0.0f);
    float _S1655;
    if(antialiased_8)
    {
        float _S1656 = _S1520 * _S1631;
        _S1522 = _S1515 * _S1631;
        _S1655 = _S1656;
    }
    else
    {
        _S1522 = _S1631;
        _S1655 = 0.0f;
    }
    float _S1657 = - (_S1522 / _S1521);
    DiffPair_float_0 _S1658;
    (&_S1658)->primal_0 = _S1518;
    (&_S1658)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1658, _S1657);
    float _S1659 = - _S1658.differential_0;
    float _S1660 = invdet_9 * _S1634.rows[int(1)].y;
    float _S1661 = - (invdet_9 * _S1634.rows[int(1)].x);
    float _S1662 = - (invdet_9 * _S1634.rows[int(0)].y);
    float _S1663 = invdet_9 * _S1634.rows[int(0)].x;
    float _S1664 = - ((_S1507 * _S1634.rows[int(1)].y + _S1517 * _S1634.rows[int(1)].x + _S1516 * _S1634.rows[int(0)].y + _S1509 * _S1634.rows[int(0)].x) / _S1513);
    DiffPair_float_0 _S1665;
    (&_S1665)->primal_0 = _S1514;
    (&_S1665)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1665, _S1655);
    DiffPair_float_0 _S1666;
    (&_S1666)->primal_0 = 0.0f;
    (&_S1666)->differential_0 = 0.0f;
    DiffPair_float_0 _S1667;
    (&_S1667)->primal_0 = _S1512;
    (&_S1667)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1666, &_S1667, _S1665.differential_0);
    float _S1668 = _S1667.differential_0 / _S1513;
    float s_diff_det_orig_T_2 = det_blur_6 * _S1668;
    float _S1669 = _S1664 + det_orig_9 * - _S1668;
    float _S1670 = - _S1669;
    float _S1671 = _S1507 * _S1669;
    float _S1672 = _S1509 * _S1669;
    Matrix<float, 2, 2>  _S1673 = _S1632;
    _S1673[int(1)] = _S1653;
    _S1673[int(0)] = _S1654;
    _S1508 = _S1673;
    *&(((&_S1508)->rows + (int(1)))->y) = 0.0f;
    float _S1674 = _S1663 + _S1671 + _S1673.rows[int(1)].y;
    *&(((&_S1508)->rows + (int(0)))->x) = 0.0f;
    float _S1675 = _S1660 + _S1672 + _S1673.rows[int(0)].x;
    float _S1676 = _S1670 + - s_diff_det_orig_T_2;
    float _S1677 = _S1661 + _S1505.rows[int(0)].y * _S1676;
    float _S1678 = _S1662 + _S1505.rows[int(1)].x * _S1676;
    float _S1679 = _S1505.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S1680 = _S1674 + _S1505.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S1681 = _S1621;
    *&((&_S1681)->x) = _S1677;
    *&((&_S1681)->y) = _S1680;
    float _S1682 = _S1675 + _S1679;
    float2  _S1683 = _S1621;
    *&((&_S1683)->y) = _S1678;
    *&((&_S1683)->x) = _S1682;
    float _S1684 = fy_12 * v_mean2d_2.y;
    float _S1685 = fx_12 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S1686 = _S1632;
    _S1686[int(1)] = _S1681;
    _S1686[int(0)] = _S1683;
    Matrix<float, 2, 2>  _S1687 = _S1508 + _S1686;
    Matrix<float, 2, 3>  _S1688 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1689;
    (&_S1689)->primal_0 = _S1503;
    (&_S1689)->differential_0 = _S1688;
    Matrix<float, 3, 2>  _S1690 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1691;
    (&_S1691)->primal_0 = _S1504;
    (&_S1691)->differential_0 = _S1690;
    s_bwd_prop_mul_1(&_S1689, &_S1691, _S1687);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1692;
    (&_S1692)->primal_0 = J_11;
    (&_S1692)->differential_0 = _S1688;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1693;
    (&_S1693)->primal_0 = _S1502;
    (&_S1693)->differential_0 = _S1616;
    s_bwd_prop_mul_2(&_S1692, &_S1693, _S1689.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1694;
    (&_S1694)->primal_0 = _S1500;
    (&_S1694)->differential_0 = _S1616;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1695;
    (&_S1695)->primal_0 = _S1501;
    (&_S1695)->differential_0 = _S1616;
    s_bwd_prop_mul_3(&_S1694, &_S1695, _S1693.differential_0);
    Matrix<float, 3, 3>  _S1696 = transpose_0(_S1695.differential_0 + _S1619.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1697;
    (&_S1697)->primal_0 = R_10;
    (&_S1697)->differential_0 = _S1616;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1698;
    (&_S1698)->primal_0 = _S1499;
    (&_S1698)->differential_0 = _S1616;
    s_bwd_prop_mul_3(&_S1697, &_S1698, _S1694.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1699;
    (&_S1699)->primal_0 = _S1497;
    (&_S1699)->differential_0 = _S1616;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1700;
    (&_S1700)->primal_0 = _S1498;
    (&_S1700)->differential_0 = _S1616;
    s_bwd_prop_mul_3(&_S1699, &_S1700, _S1698.differential_0);
    Matrix<float, 3, 3>  _S1701 = _S1699.differential_0 + transpose_0(_S1700.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1702;
    (&_S1702)->primal_0 = _S1496;
    (&_S1702)->differential_0 = _S1616;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1703;
    (&_S1703)->primal_0 = S_2;
    (&_S1703)->differential_0 = _S1616;
    s_bwd_prop_mul_3(&_S1702, &_S1703, _S1701);
    Matrix<float, 3, 3>  _S1704 = transpose_0(_S1702.differential_0);
    float _S1705 = 2.0f * - _S1704.rows[int(2)].z;
    float _S1706 = 2.0f * _S1704.rows[int(2)].y;
    float _S1707 = 2.0f * _S1704.rows[int(2)].x;
    float _S1708 = 2.0f * _S1704.rows[int(1)].z;
    float _S1709 = 2.0f * - _S1704.rows[int(1)].y;
    float _S1710 = 2.0f * _S1704.rows[int(1)].x;
    float _S1711 = 2.0f * _S1704.rows[int(0)].z;
    float _S1712 = 2.0f * _S1704.rows[int(0)].y;
    float _S1713 = 2.0f * - _S1704.rows[int(0)].x;
    float _S1714 = - _S1710 + _S1712;
    float _S1715 = _S1707 + - _S1711;
    float _S1716 = - _S1706 + _S1708;
    float _S1717 = _S1706 + _S1708;
    float _S1718 = _S1707 + _S1711;
    float _S1719 = _S1710 + _S1712;
    float _S1720 = z_19 * (_S1709 + _S1713);
    float _S1721 = y_22 * (_S1705 + _S1713);
    float _S1722 = x_37 * (_S1705 + _S1709);
    float _S1723 = z_19 * _S1714 + y_22 * _S1715 + x_37 * _S1716;
    float _S1724 = _S1495 * _S1723;
    float _S1725 = w_11 * _S1714 + y_22 * _S1717 + x_37 * _S1718 + _S1720 + _S1720;
    float _S1726 = _S1495 * _S1725;
    float _S1727 = w_11 * _S1715 + z_19 * _S1717 + x_37 * _S1719 + _S1721 + _S1721;
    float _S1728 = _S1495 * _S1727;
    float _S1729 = w_11 * _S1716 + z_19 * _S1718 + y_22 * _S1719 + _S1722 + _S1722;
    float _S1730 = _S1495 * _S1729;
    float _S1731 = quat_11.x * _S1723 + quat_11.w * _S1725 + quat_11.z * _S1727 + quat_11.y * _S1729;
    DiffPair_float_0 _S1732;
    (&_S1732)->primal_0 = _S1494;
    (&_S1732)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1732, _S1731);
    float _S1733 = quat_11.x * _S1732.differential_0;
    float _S1734 = quat_11.w * _S1732.differential_0;
    float _S1735 = quat_11.z * _S1732.differential_0;
    float _S1736 = quat_11.y * _S1732.differential_0;
    float _S1737 = _S1726 + _S1734 + _S1734;
    float _S1738 = _S1728 + _S1735 + _S1735;
    float _S1739 = _S1730 + _S1736 + _S1736;
    float _S1740 = _S1724 + _S1733 + _S1733;
    float3  _S1741 = _S1554;
    *&((&_S1741)->z) = _S1703.differential_0.rows[int(2)].z;
    *&((&_S1741)->y) = _S1703.differential_0.rows[int(1)].y;
    *&((&_S1741)->x) = _S1703.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1742;
    (&_S1742)->primal_0 = scale_10;
    (&_S1742)->differential_0 = _S1554;
    s_bwd_prop_exp_1(&_S1742, _S1741);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1743 = _S1742;
    float3  _S1744 = _S1554;
    *&((&_S1744)->y) = _S1684;
    *&((&_S1744)->x) = _S1685;
    float3  _S1745 = _S1652 + _S1744;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1746;
    (&_S1746)->primal_0 = R_10;
    (&_S1746)->differential_0 = _S1616;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1747;
    (&_S1747)->primal_0 = mean_8;
    (&_S1747)->differential_0 = _S1554;
    s_bwd_prop_mul_0(&_S1746, &_S1747, _S1745);
    float3  _S1748 = _S1745 + _S1620.differential_0;
    Matrix<float, 3, 3>  _S1749 = _S1696 + _S1697.differential_0 + _S1746.differential_0;
    float4  _S1750 = make_float4 (0.0f);
    *&((&_S1750)->w) = _S1737;
    *&((&_S1750)->z) = _S1738;
    *&((&_S1750)->y) = _S1739;
    *&((&_S1750)->x) = _S1740;
    float4  _S1751 = _S1750;
    float3  _S1752 = _S1747.differential_0 + _S1614;
    *v_mean_2 = _S1752;
    *v_quat_2 = _S1751;
    *v_scale_2 = _S1743.differential_0;
    *v_in_opacity_2 = _S1659;
    (*v_sh_coeffs_2)[int(0)] = _S1636;
    (*v_sh_coeffs_2)[int(1)] = _S1637;
    (*v_sh_coeffs_2)[int(2)] = _S1638;
    (*v_sh_coeffs_2)[int(3)] = _S1639;
    (*v_sh_coeffs_2)[int(4)] = _S1640;
    (*v_sh_coeffs_2)[int(5)] = _S1641;
    (*v_sh_coeffs_2)[int(6)] = _S1642;
    (*v_sh_coeffs_2)[int(7)] = _S1643;
    (*v_sh_coeffs_2)[int(8)] = _S1644;
    (*v_sh_coeffs_2)[int(9)] = _S1645;
    (*v_sh_coeffs_2)[int(10)] = _S1646;
    (*v_sh_coeffs_2)[int(11)] = _S1647;
    (*v_sh_coeffs_2)[int(12)] = _S1648;
    (*v_sh_coeffs_2)[int(13)] = _S1649;
    (*v_sh_coeffs_2)[int(14)] = _S1650;
    (*v_sh_coeffs_2)[int(15)] = _S1651;
    *v_R_2 = _S1749;
    *v_t_2 = _S1748;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_12, float dOut_15)
{
    float _S1753 = _slang_select(((*dpx_12).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_12).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_15;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S1753;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_13, float dOut_16)
{
    float _S1754 = (F32_exp2(((*dpx_13).primal_0))) * 50.693145751953125f * dOut_16;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S1754;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_17)
{
    float _S1755 = dOut_17.y;
    float _S1756 = dOut_17.z;
    float _S1757 = dOut_17.x;
    float _S1758 = (*a_0).primal_0.z * _S1755 + - (*a_0).primal_0.y * _S1756;
    float _S1759 = - (*a_0).primal_0.z * _S1757 + (*a_0).primal_0.x * _S1756;
    float _S1760 = (*a_0).primal_0.y * _S1757 + - (*a_0).primal_0.x * _S1755;
    float3  _S1761 = make_float3 (- (*b_0).primal_0.z * _S1755 + (*b_0).primal_0.y * _S1756, (*b_0).primal_0.z * _S1757 + - (*b_0).primal_0.x * _S1756, - (*b_0).primal_0.y * _S1757 + (*b_0).primal_0.x * _S1755);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S1761;
    float3  _S1762 = make_float3 (_S1758, _S1759, _S1760);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S1762;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S1763 = left_8.y;
    float _S1764 = right_8.z;
    float _S1765 = left_8.z;
    float _S1766 = right_8.y;
    float _S1767 = right_8.x;
    float _S1768 = left_8.x;
    return make_float3 (_S1763 * _S1764 - _S1765 * _S1766, _S1765 * _S1767 - _S1768 * _S1764, _S1768 * _S1766 - _S1763 * _S1767);
}

inline __device__ float3  normalize_0(float3  x_39)
{
    return x_39 / make_float3 (length_1(x_39));
}

inline __device__ float2  normalize_1(float2  x_40)
{
    return x_40 / make_float2 (length_0(x_40));
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_9, float4  quat_12, float3  scale_11, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_13, float fy_13, float cx_13, float cy_13, float4  radial_coeffs_12, float2  tangential_coeffs_12, float2  thin_prism_coeffs_12, uint image_width_9, uint image_height_9, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float * depth_6, float2  * out_hardness_0, float3  * rgb_6, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_8 = mul_0(R_11, mean_9) + t_10;
        float _S1769 = mean_c_8.z;
        bool _S1770;
        if(_S1769 < near_plane_6)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = _S1769 > far_plane_6;
        }
        if(_S1770)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S1771 = scale_11.x;
        float sx_0 = (F32_exp((_S1771)));
        float _S1772 = scale_11.y;
        float sy_0 = (F32_exp((_S1772)));
        float sz_0 = scale_11.z - 0.5f * (_S1771 + _S1772);
        float x_41 = quat_12.y;
        float inv_norm_9 = (F32_rsqrt((x_41 * x_41 + quat_12.z * quat_12.z + quat_12.w * quat_12.w + quat_12.x * quat_12.x)));
        float x_42 = quat_12.y * inv_norm_9;
        float y_24 = quat_12.z * inv_norm_9;
        float z_21 = quat_12.w * inv_norm_9;
        float w_12 = quat_12.x * inv_norm_9;
        float x2_12 = x_42 * x_42;
        float y2_12 = y_24 * y_24;
        float z2_21 = z_21 * z_21;
        float xy_12 = x_42 * y_24;
        float xz_12 = x_42 * z_21;
        float yz_12 = y_24 * z_21;
        float wx_12 = w_12 * x_42;
        float wy_12 = w_12 * y_24;
        float wz_12 = w_12 * z_21;
        Matrix<float, 3, 3>  _S1773 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12)));
        float3  vert0_0 = mul_0(_S1773, make_float3 (sx_0, 0.0f, 0.0f)) + mean_9;
        float3  vert1_0 = mul_0(_S1773, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_9;
        float3  vert2_0 = mul_0(_S1773, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_9;
        float3  vert0_c_0 = mul_0(R_11, vert0_0) + t_10;
        float3  vert1_c_0 = mul_0(R_11, vert1_0) + t_10;
        float3  vert2_c_0 = mul_0(R_11, vert2_0) + t_10;
        float _S1774 = vert0_c_0.z;
        if(_S1774 < near_plane_6)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = _S1774 > far_plane_6;
        }
        if(_S1770)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = (vert1_c_0.z) < near_plane_6;
        }
        if(_S1770)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = (vert1_c_0.z) > far_plane_6;
        }
        if(_S1770)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = (vert2_c_0.z) < near_plane_6;
        }
        if(_S1770)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = (vert2_c_0.z) > far_plane_6;
        }
        if(_S1770)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S1774);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (vert1_c_0.z);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (vert2_c_0.z);
        float2  _S1775 = make_float2 (fx_13, fy_13);
        float2  _S1776 = make_float2 (cx_13, cy_13);
        *uv0_0 = _S1775 * *uv0_0 + _S1776;
        *uv1_0 = _S1775 * *uv1_0 + _S1776;
        float2  _S1777 = _S1775 * *uv2_0 + _S1776;
        *uv2_0 = _S1777;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S1777 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S1777)));
        float _S1778 = _S1777.x;
        float xmax_3 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S1778))) + offset_0;
        float xmin_3 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S1778))) - offset_0;
        float _S1779 = _S1777.y;
        float ymax_3 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S1779))) + offset_0;
        float ymin_3 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S1779))) - offset_0;
        if(xmax_3 <= 0.0f)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = xmin_3 >= float(image_width_9);
        }
        if(_S1770)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = ymax_3 <= 0.0f;
        }
        if(_S1770)
        {
            _S1770 = true;
        }
        else
        {
            _S1770 = ymin_3 >= float(image_height_9);
        }
        if(_S1770)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_6 = make_int4 (int((F32_floor((xmin_3)))), int((F32_floor((ymin_3)))), int((F32_ceil((xmax_3)))), int((F32_ceil((ymax_3)))));
        *depth_6 = _S1769;
        *out_hardness_0 = hardness_0;
        float3  _S1780 = mean_9 - - mul_0(transpose_0(R_11), t_10);
        float3  _S1781 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
        *rgb_6 = _S1781;
        float _S1782 = _S1780.x;
        float _S1783 = _S1780.y;
        float _S1784 = _S1780.z;
        float norm_6 = (F32_sqrt((_S1782 * _S1782 + _S1783 * _S1783 + _S1784 * _S1784)));
        float x_43 = _S1782 / norm_6;
        float y_25 = _S1783 / norm_6;
        float z_22 = _S1784 / norm_6;
        float3  _S1785 = _S1781 + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_9)[int(1)] + make_float3 (z_22) * (*sh_coeffs_9)[int(2)] - make_float3 (x_43) * (*sh_coeffs_9)[int(3)]);
        *rgb_6 = _S1785;
        float z2_22 = z_22 * z_22;
        float fTmp0B_9 = -1.09254848957061768f * z_22;
        float fC1_9 = x_43 * x_43 - y_25 * y_25;
        float fS1_9 = 2.0f * x_43 * y_25;
        float3  _S1786 = _S1785 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_25) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_22 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_43) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
        *rgb_6 = _S1786;
        float fTmp0C_9 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
        float fTmp1B_9 = 1.44530570507049561f * z_22;
        *rgb_6 = max_0(_S1786 + (make_float3 (-0.59004360437393188f * (x_43 * fS1_9 + y_25 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_25) * (*sh_coeffs_9)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_22 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_43) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_43 * fC1_9 - y_25 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        float3  _S1787 = normalize_0(mul_0(R_11, cross_0(vert1_0 - vert0_0, vert2_0 - vert0_0)));
        *normal_0 = _S1787 * make_float3 (float(- (F32_sign((dot_0(_S1787, mean_c_8))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_10, float4  quat_13, float3  scale_12, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_14, float fy_14, float cx_14, float cy_14, float4  radial_coeffs_13, float2  tangential_coeffs_13, float2  thin_prism_coeffs_13, uint image_width_10, uint image_height_10, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float * depth_7, float2  * out_hardness_1, float3  * rgb_7, float3  * normal_1)
{
    for(;;)
    {
        float3  mean_c_9 = mul_0(R_12, mean_10) + t_11;
        float _S1788 = mean_c_9.z;
        bool _S1789;
        if(_S1788 < near_plane_7)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = _S1788 > far_plane_7;
        }
        if(_S1789)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S1790 = scale_12.x;
        float sx_1 = (F32_exp((_S1790)));
        float _S1791 = scale_12.y;
        float sy_1 = (F32_exp((_S1791)));
        float sz_1 = scale_12.z - 0.5f * (_S1790 + _S1791);
        float x_44 = quat_13.y;
        float inv_norm_10 = (F32_rsqrt((x_44 * x_44 + quat_13.z * quat_13.z + quat_13.w * quat_13.w + quat_13.x * quat_13.x)));
        float x_45 = quat_13.y * inv_norm_10;
        float y_26 = quat_13.z * inv_norm_10;
        float z_23 = quat_13.w * inv_norm_10;
        float w_13 = quat_13.x * inv_norm_10;
        float x2_13 = x_45 * x_45;
        float y2_13 = y_26 * y_26;
        float z2_23 = z_23 * z_23;
        float xy_13 = x_45 * y_26;
        float xz_13 = x_45 * z_23;
        float yz_13 = y_26 * z_23;
        float wx_13 = w_13 * x_45;
        float wy_13 = w_13 * y_26;
        float wz_13 = w_13 * z_23;
        Matrix<float, 3, 3>  _S1792 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
        float3  vert0_1 = mul_0(_S1792, make_float3 (sx_1, 0.0f, 0.0f)) + mean_10;
        float3  vert1_1 = mul_0(_S1792, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_10;
        float3  vert2_1 = mul_0(_S1792, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_10;
        float3  vert0_c_1 = mul_0(R_12, vert0_1) + t_11;
        float3  vert1_c_1 = mul_0(R_12, vert1_1) + t_11;
        float3  vert2_c_1 = mul_0(R_12, vert2_1) + t_11;
        float _S1793 = vert0_c_1.z;
        if(_S1793 < near_plane_7)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = _S1793 > far_plane_7;
        }
        if(_S1789)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = (vert1_c_1.z) < near_plane_7;
        }
        if(_S1789)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = (vert1_c_1.z) > far_plane_7;
        }
        if(_S1789)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = (vert2_c_1.z) < near_plane_7;
        }
        if(_S1789)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = (vert2_c_1.z) > far_plane_7;
        }
        if(_S1789)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_1 = CameraDistortion_x24init_0(radial_coeffs_13, tangential_coeffs_13, thin_prism_coeffs_13);
        float2  _S1794 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_7 = length_0(_S1794);
        float theta_1 = (F32_atan2((r_7), (_S1793)));
        float k_4;
        if(theta_1 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_1 * theta_1 / 3.0f) / _S1793;
        }
        else
        {
            k_4 = theta_1 / r_7;
        }
        float2  _S1795 = _S1794 * make_float2 (k_4);
        float k1_3 = dist_coeffs_1.radial_coeffs_0.x;
        float k2_3 = dist_coeffs_1.radial_coeffs_0.y;
        float k3_3 = dist_coeffs_1.radial_coeffs_0.z;
        float k4_3 = dist_coeffs_1.radial_coeffs_0.w;
        float p1_3 = dist_coeffs_1.tangential_coeffs_0.x;
        float p2_3 = dist_coeffs_1.tangential_coeffs_0.y;
        float sx1_3 = dist_coeffs_1.thin_prism_coeffs_0.x;
        float sy1_3 = dist_coeffs_1.thin_prism_coeffs_0.y;
        float u_9 = _S1795.x;
        float v_9 = _S1795.y;
        float r2_9 = u_9 * u_9 + v_9 * v_9;
        float _S1796 = 2.0f * p1_3;
        float _S1797 = 2.0f * p2_3;
        float2  _S1798 = _S1795 * make_float2 (1.0f + r2_9 * (k1_3 + r2_9 * (k2_3 + r2_9 * (k3_3 + r2_9 * k4_3)))) + make_float2 (_S1796 * u_9 * v_9 + p2_3 * (r2_9 + 2.0f * u_9 * u_9) + sx1_3 * r2_9, _S1797 * u_9 * v_9 + p1_3 * (r2_9 + 2.0f * v_9 * v_9) + sy1_3 * r2_9);
        *uv0_1 = make_float2 (fx_14 * _S1798.x + cx_14, fy_14 * _S1798.y + cy_14);
        float2  _S1799 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_8 = length_0(_S1799);
        float _S1800 = vert1_c_1.z;
        float theta_2 = (F32_atan2((r_8), (_S1800)));
        if(theta_2 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_2 * theta_2 / 3.0f) / _S1800;
        }
        else
        {
            k_4 = theta_2 / r_8;
        }
        float2  _S1801 = _S1799 * make_float2 (k_4);
        float u_10 = _S1801.x;
        float v_10 = _S1801.y;
        float r2_10 = u_10 * u_10 + v_10 * v_10;
        float2  _S1802 = _S1801 * make_float2 (1.0f + r2_10 * (k1_3 + r2_10 * (k2_3 + r2_10 * (k3_3 + r2_10 * k4_3)))) + make_float2 (_S1796 * u_10 * v_10 + p2_3 * (r2_10 + 2.0f * u_10 * u_10) + sx1_3 * r2_10, _S1797 * u_10 * v_10 + p1_3 * (r2_10 + 2.0f * v_10 * v_10) + sy1_3 * r2_10);
        *uv1_1 = make_float2 (fx_14 * _S1802.x + cx_14, fy_14 * _S1802.y + cy_14);
        float2  _S1803 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_9 = length_0(_S1803);
        float _S1804 = vert2_c_1.z;
        float theta_3 = (F32_atan2((r_9), (_S1804)));
        if(theta_3 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_3 * theta_3 / 3.0f) / _S1804;
        }
        else
        {
            k_4 = theta_3 / r_9;
        }
        float2  _S1805 = _S1803 * make_float2 (k_4);
        float u_11 = _S1805.x;
        float v_11 = _S1805.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float2  _S1806 = _S1805 * make_float2 (1.0f + r2_11 * (k1_3 + r2_11 * (k2_3 + r2_11 * (k3_3 + r2_11 * k4_3)))) + make_float2 (_S1796 * u_11 * v_11 + p2_3 * (r2_11 + 2.0f * u_11 * u_11) + sx1_3 * r2_11, _S1797 * u_11 * v_11 + p1_3 * (r2_11 + 2.0f * v_11 * v_11) + sy1_3 * r2_11);
        float _S1807 = fx_14 * _S1806.x + cx_14;
        float _S1808 = fy_14 * _S1806.y + cy_14;
        float2  _S1809 = make_float2 (_S1807, _S1808);
        *uv2_1 = _S1809;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S1809 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S1809)));
        float xmax_4 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S1807))) + offset_1;
        float xmin_4 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S1807))) - offset_1;
        float ymax_4 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S1808))) + offset_1;
        float ymin_4 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S1808))) - offset_1;
        if(xmax_4 <= 0.0f)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = xmin_4 >= float(image_width_10);
        }
        if(_S1789)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = ymax_4 <= 0.0f;
        }
        if(_S1789)
        {
            _S1789 = true;
        }
        else
        {
            _S1789 = ymin_4 >= float(image_height_10);
        }
        if(_S1789)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_7 = make_int4 (int((F32_floor((xmin_4)))), int((F32_floor((ymin_4)))), int((F32_ceil((xmax_4)))), int((F32_ceil((ymax_4)))));
        *depth_7 = _S1788;
        *out_hardness_1 = hardness_1;
        float3  _S1810 = mean_10 - - mul_0(transpose_0(R_12), t_11);
        float3  _S1811 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)];
        *rgb_7 = _S1811;
        float _S1812 = _S1810.x;
        float _S1813 = _S1810.y;
        float _S1814 = _S1810.z;
        float norm_7 = (F32_sqrt((_S1812 * _S1812 + _S1813 * _S1813 + _S1814 * _S1814)));
        float x_46 = _S1812 / norm_7;
        float y_27 = _S1813 / norm_7;
        float z_24 = _S1814 / norm_7;
        float3  _S1815 = _S1811 + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_10)[int(1)] + make_float3 (z_24) * (*sh_coeffs_10)[int(2)] - make_float3 (x_46) * (*sh_coeffs_10)[int(3)]);
        *rgb_7 = _S1815;
        float z2_24 = z_24 * z_24;
        float fTmp0B_10 = -1.09254848957061768f * z_24;
        float fC1_10 = x_46 * x_46 - y_27 * y_27;
        float fS1_10 = 2.0f * x_46 * y_27;
        float3  _S1816 = _S1815 + (make_float3 (0.54627424478530884f * fS1_10) * (*sh_coeffs_10)[int(4)] + make_float3 (fTmp0B_10 * y_27) * (*sh_coeffs_10)[int(5)] + make_float3 (0.94617468118667603f * z2_24 - 0.31539157032966614f) * (*sh_coeffs_10)[int(6)] + make_float3 (fTmp0B_10 * x_46) * (*sh_coeffs_10)[int(7)] + make_float3 (0.54627424478530884f * fC1_10) * (*sh_coeffs_10)[int(8)]);
        *rgb_7 = _S1816;
        float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
        float fTmp1B_10 = 1.44530570507049561f * z_24;
        *rgb_7 = max_0(_S1816 + (make_float3 (-0.59004360437393188f * (x_46 * fS1_10 + y_27 * fC1_10)) * (*sh_coeffs_10)[int(9)] + make_float3 (fTmp1B_10 * fS1_10) * (*sh_coeffs_10)[int(10)] + make_float3 (fTmp0C_10 * y_27) * (*sh_coeffs_10)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_24 - 1.11952900886535645f)) * (*sh_coeffs_10)[int(12)] + make_float3 (fTmp0C_10 * x_46) * (*sh_coeffs_10)[int(13)] + make_float3 (fTmp1B_10 * fC1_10) * (*sh_coeffs_10)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_10 - y_27 * fS1_10)) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        float3  _S1817 = normalize_0(mul_0(R_12, cross_0(vert1_1 - vert0_1, vert2_1 - vert0_1)));
        *normal_1 = _S1817 * make_float3 (float(- (F32_sign((dot_0(_S1817, mean_c_9))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_11, float4  quat_14, float3  scale_13, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_15, float fy_15, float cx_15, float cy_15, float4  radial_coeffs_14, float2  tangential_coeffs_14, float2  thin_prism_coeffs_14, uint image_width_11, uint image_height_11, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float * depth_8, float2  * out_hardness_2, float3  * rgb_8, float3  * normal_2)
{
    float3  mean_c_10 = mul_0(R_13, mean_11) + t_12;
    float _S1818 = scale_13.x;
    float sx_2 = (F32_exp((_S1818)));
    float _S1819 = scale_13.y;
    float sy_2 = (F32_exp((_S1819)));
    float sz_2 = scale_13.z - 0.5f * (_S1818 + _S1819);
    float x_47 = quat_14.y;
    float inv_norm_11 = (F32_rsqrt((x_47 * x_47 + quat_14.z * quat_14.z + quat_14.w * quat_14.w + quat_14.x * quat_14.x)));
    float x_48 = quat_14.y * inv_norm_11;
    float y_28 = quat_14.z * inv_norm_11;
    float z_25 = quat_14.w * inv_norm_11;
    float w_14 = quat_14.x * inv_norm_11;
    float x2_14 = x_48 * x_48;
    float y2_14 = y_28 * y_28;
    float z2_25 = z_25 * z_25;
    float xy_14 = x_48 * y_28;
    float xz_14 = x_48 * z_25;
    float yz_14 = y_28 * z_25;
    float wx_14 = w_14 * x_48;
    float wy_14 = w_14 * y_28;
    float wz_14 = w_14 * z_25;
    Matrix<float, 3, 3>  _S1820 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    float3  vert0_2 = mul_0(_S1820, make_float3 (sx_2, 0.0f, 0.0f)) + mean_11;
    float3  vert1_2 = mul_0(_S1820, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_11;
    float3  vert2_2 = mul_0(_S1820, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_11;
    float3  vert0_c_2 = mul_0(R_13, vert0_2) + t_12;
    float3  vert1_c_2 = mul_0(R_13, vert1_2) + t_12;
    float3  vert2_c_2 = mul_0(R_13, vert2_2) + t_12;
    *uv0_2 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    *uv1_2 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    *uv2_2 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float2  _S1821 = make_float2 (fx_15, fy_15);
    float2  _S1822 = make_float2 (cx_15, cy_15);
    *uv0_2 = _S1821 * *uv0_2 + _S1822;
    *uv1_2 = _S1821 * *uv1_2 + _S1822;
    float2  _S1823 = _S1821 * *uv2_2 + _S1822;
    *uv2_2 = _S1823;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S1823 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S1823)));
    float _S1824 = _S1823.x;
    float _S1825 = _S1823.y;
    *aabb_xyxy_8 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S1824))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S1825))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S1824))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S1825))) + offset_2)))));
    *depth_8 = mean_c_10.z;
    *out_hardness_2 = hardness_2;
    float3  _S1826 = mean_11 - - mul_0(transpose_0(R_13), t_12);
    float3  _S1827 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)];
    *rgb_8 = _S1827;
    float _S1828 = _S1826.x;
    float _S1829 = _S1826.y;
    float _S1830 = _S1826.z;
    float norm_8 = (F32_sqrt((_S1828 * _S1828 + _S1829 * _S1829 + _S1830 * _S1830)));
    float x_49 = _S1828 / norm_8;
    float y_29 = _S1829 / norm_8;
    float z_26 = _S1830 / norm_8;
    float3  _S1831 = _S1827 + make_float3 (0.48860251903533936f) * (make_float3 (- y_29) * (*sh_coeffs_11)[int(1)] + make_float3 (z_26) * (*sh_coeffs_11)[int(2)] - make_float3 (x_49) * (*sh_coeffs_11)[int(3)]);
    *rgb_8 = _S1831;
    float z2_26 = z_26 * z_26;
    float fTmp0B_11 = -1.09254848957061768f * z_26;
    float fC1_11 = x_49 * x_49 - y_29 * y_29;
    float fS1_11 = 2.0f * x_49 * y_29;
    float3  _S1832 = _S1831 + (make_float3 (0.54627424478530884f * fS1_11) * (*sh_coeffs_11)[int(4)] + make_float3 (fTmp0B_11 * y_29) * (*sh_coeffs_11)[int(5)] + make_float3 (0.94617468118667603f * z2_26 - 0.31539157032966614f) * (*sh_coeffs_11)[int(6)] + make_float3 (fTmp0B_11 * x_49) * (*sh_coeffs_11)[int(7)] + make_float3 (0.54627424478530884f * fC1_11) * (*sh_coeffs_11)[int(8)]);
    *rgb_8 = _S1832;
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_26;
    *rgb_8 = max_0(_S1832 + (make_float3 (-0.59004360437393188f * (x_49 * fS1_11 + y_29 * fC1_11)) * (*sh_coeffs_11)[int(9)] + make_float3 (fTmp1B_11 * fS1_11) * (*sh_coeffs_11)[int(10)] + make_float3 (fTmp0C_11 * y_29) * (*sh_coeffs_11)[int(11)] + make_float3 (z_26 * (1.86588168144226074f * z2_26 - 1.11952900886535645f)) * (*sh_coeffs_11)[int(12)] + make_float3 (fTmp0C_11 * x_49) * (*sh_coeffs_11)[int(13)] + make_float3 (fTmp1B_11 * fC1_11) * (*sh_coeffs_11)[int(14)] + make_float3 (-0.59004360437393188f * (x_49 * fC1_11 - y_29 * fS1_11)) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    float3  _S1833 = normalize_0(mul_0(R_13, cross_0(vert1_2 - vert0_2, vert2_2 - vert0_2)));
    *normal_2 = _S1833 * make_float3 (float(- (F32_sign((dot_0(_S1833, mean_c_10))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_12, float4  quat_15, float3  scale_14, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_16, float fy_16, float cx_16, float cy_16, float4  radial_coeffs_15, float2  tangential_coeffs_15, float2  thin_prism_coeffs_15, uint image_width_12, uint image_height_12, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float * depth_9, float2  * out_hardness_3, float3  * rgb_9, float3  * normal_3)
{
    float3  mean_c_11 = mul_0(R_14, mean_12) + t_13;
    float _S1834 = scale_14.x;
    float sx_3 = (F32_exp((_S1834)));
    float _S1835 = scale_14.y;
    float sy_3 = (F32_exp((_S1835)));
    float sz_3 = scale_14.z - 0.5f * (_S1834 + _S1835);
    float x_50 = quat_15.y;
    float inv_norm_12 = (F32_rsqrt((x_50 * x_50 + quat_15.z * quat_15.z + quat_15.w * quat_15.w + quat_15.x * quat_15.x)));
    float x_51 = quat_15.y * inv_norm_12;
    float y_30 = quat_15.z * inv_norm_12;
    float z_27 = quat_15.w * inv_norm_12;
    float w_15 = quat_15.x * inv_norm_12;
    float x2_15 = x_51 * x_51;
    float y2_15 = y_30 * y_30;
    float z2_27 = z_27 * z_27;
    float xy_15 = x_51 * y_30;
    float xz_15 = x_51 * z_27;
    float yz_15 = y_30 * z_27;
    float wx_15 = w_15 * x_51;
    float wy_15 = w_15 * y_30;
    float wz_15 = w_15 * z_27;
    Matrix<float, 3, 3>  _S1836 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    float3  vert0_3 = mul_0(_S1836, make_float3 (sx_3, 0.0f, 0.0f)) + mean_12;
    float3  vert1_3 = mul_0(_S1836, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_12;
    float3  vert2_3 = mul_0(_S1836, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_12;
    float3  vert0_c_3 = mul_0(R_14, vert0_3) + t_13;
    float3  vert1_c_3 = mul_0(R_14, vert1_3) + t_13;
    float3  vert2_c_3 = mul_0(R_14, vert2_3) + t_13;
    CameraDistortion_0 dist_coeffs_2 = CameraDistortion_x24init_0(radial_coeffs_15, tangential_coeffs_15, thin_prism_coeffs_15);
    float2  _S1837 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_10 = length_0(_S1837);
    float _S1838 = vert0_c_3.z;
    float theta_4 = (F32_atan2((r_10), (_S1838)));
    float k_5;
    if(theta_4 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_4 * theta_4 / 3.0f) / _S1838;
    }
    else
    {
        k_5 = theta_4 / r_10;
    }
    float2  _S1839 = _S1837 * make_float2 (k_5);
    float k1_4 = dist_coeffs_2.radial_coeffs_0.x;
    float k2_4 = dist_coeffs_2.radial_coeffs_0.y;
    float k3_4 = dist_coeffs_2.radial_coeffs_0.z;
    float k4_4 = dist_coeffs_2.radial_coeffs_0.w;
    float p1_4 = dist_coeffs_2.tangential_coeffs_0.x;
    float p2_4 = dist_coeffs_2.tangential_coeffs_0.y;
    float sx1_4 = dist_coeffs_2.thin_prism_coeffs_0.x;
    float sy1_4 = dist_coeffs_2.thin_prism_coeffs_0.y;
    float u_12 = _S1839.x;
    float v_12 = _S1839.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S1840 = 2.0f * p1_4;
    float _S1841 = 2.0f * p2_4;
    float2  _S1842 = _S1839 * make_float2 (1.0f + r2_12 * (k1_4 + r2_12 * (k2_4 + r2_12 * (k3_4 + r2_12 * k4_4)))) + make_float2 (_S1840 * u_12 * v_12 + p2_4 * (r2_12 + 2.0f * u_12 * u_12) + sx1_4 * r2_12, _S1841 * u_12 * v_12 + p1_4 * (r2_12 + 2.0f * v_12 * v_12) + sy1_4 * r2_12);
    *uv0_3 = make_float2 (fx_16 * _S1842.x + cx_16, fy_16 * _S1842.y + cy_16);
    float2  _S1843 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_11 = length_0(_S1843);
    float _S1844 = vert1_c_3.z;
    float theta_5 = (F32_atan2((r_11), (_S1844)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_5 * theta_5 / 3.0f) / _S1844;
    }
    else
    {
        k_5 = theta_5 / r_11;
    }
    float2  _S1845 = _S1843 * make_float2 (k_5);
    float u_13 = _S1845.x;
    float v_13 = _S1845.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S1846 = _S1845 * make_float2 (1.0f + r2_13 * (k1_4 + r2_13 * (k2_4 + r2_13 * (k3_4 + r2_13 * k4_4)))) + make_float2 (_S1840 * u_13 * v_13 + p2_4 * (r2_13 + 2.0f * u_13 * u_13) + sx1_4 * r2_13, _S1841 * u_13 * v_13 + p1_4 * (r2_13 + 2.0f * v_13 * v_13) + sy1_4 * r2_13);
    *uv1_3 = make_float2 (fx_16 * _S1846.x + cx_16, fy_16 * _S1846.y + cy_16);
    float2  _S1847 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_12 = length_0(_S1847);
    float _S1848 = vert2_c_3.z;
    float theta_6 = (F32_atan2((r_12), (_S1848)));
    if(theta_6 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_6 * theta_6 / 3.0f) / _S1848;
    }
    else
    {
        k_5 = theta_6 / r_12;
    }
    float2  _S1849 = _S1847 * make_float2 (k_5);
    float u_14 = _S1849.x;
    float v_14 = _S1849.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S1850 = _S1849 * make_float2 (1.0f + r2_14 * (k1_4 + r2_14 * (k2_4 + r2_14 * (k3_4 + r2_14 * k4_4)))) + make_float2 (_S1840 * u_14 * v_14 + p2_4 * (r2_14 + 2.0f * u_14 * u_14) + sx1_4 * r2_14, _S1841 * u_14 * v_14 + p1_4 * (r2_14 + 2.0f * v_14 * v_14) + sy1_4 * r2_14);
    float _S1851 = fx_16 * _S1850.x + cx_16;
    float _S1852 = fy_16 * _S1850.y + cy_16;
    float2  _S1853 = make_float2 (_S1851, _S1852);
    *uv2_3 = _S1853;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S1853 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S1853)));
    *aabb_xyxy_9 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S1851))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S1852))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S1851))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S1852))) + offset_3)))));
    *depth_9 = mean_c_11.z;
    *out_hardness_3 = hardness_3;
    float3  _S1854 = mean_12 - - mul_0(transpose_0(R_14), t_13);
    float3  _S1855 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)];
    *rgb_9 = _S1855;
    float _S1856 = _S1854.x;
    float _S1857 = _S1854.y;
    float _S1858 = _S1854.z;
    float norm_9 = (F32_sqrt((_S1856 * _S1856 + _S1857 * _S1857 + _S1858 * _S1858)));
    float x_52 = _S1856 / norm_9;
    float y_31 = _S1857 / norm_9;
    float z_28 = _S1858 / norm_9;
    float3  _S1859 = _S1855 + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * (*sh_coeffs_12)[int(1)] + make_float3 (z_28) * (*sh_coeffs_12)[int(2)] - make_float3 (x_52) * (*sh_coeffs_12)[int(3)]);
    *rgb_9 = _S1859;
    float z2_28 = z_28 * z_28;
    float fTmp0B_12 = -1.09254848957061768f * z_28;
    float fC1_12 = x_52 * x_52 - y_31 * y_31;
    float fS1_12 = 2.0f * x_52 * y_31;
    float3  _S1860 = _S1859 + (make_float3 (0.54627424478530884f * fS1_12) * (*sh_coeffs_12)[int(4)] + make_float3 (fTmp0B_12 * y_31) * (*sh_coeffs_12)[int(5)] + make_float3 (0.94617468118667603f * z2_28 - 0.31539157032966614f) * (*sh_coeffs_12)[int(6)] + make_float3 (fTmp0B_12 * x_52) * (*sh_coeffs_12)[int(7)] + make_float3 (0.54627424478530884f * fC1_12) * (*sh_coeffs_12)[int(8)]);
    *rgb_9 = _S1860;
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_28;
    *rgb_9 = max_0(_S1860 + (make_float3 (-0.59004360437393188f * (x_52 * fS1_12 + y_31 * fC1_12)) * (*sh_coeffs_12)[int(9)] + make_float3 (fTmp1B_12 * fS1_12) * (*sh_coeffs_12)[int(10)] + make_float3 (fTmp0C_12 * y_31) * (*sh_coeffs_12)[int(11)] + make_float3 (z_28 * (1.86588168144226074f * z2_28 - 1.11952900886535645f)) * (*sh_coeffs_12)[int(12)] + make_float3 (fTmp0C_12 * x_52) * (*sh_coeffs_12)[int(13)] + make_float3 (fTmp1B_12 * fC1_12) * (*sh_coeffs_12)[int(14)] + make_float3 (-0.59004360437393188f * (x_52 * fC1_12 - y_31 * fS1_12)) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    float3  _S1861 = normalize_0(mul_0(R_14, cross_0(vert1_3 - vert0_3, vert2_3 - vert0_3)));
    *normal_3 = _S1861 * make_float3 (float(- (F32_sign((dot_0(_S1861, mean_c_11))))));
    return;
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S1862, float3  _S1863)
{
    return cross_0(_S1862, _S1863);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S1864, float3  _S1865)
{
    return dot_0(_S1864, _S1865);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1866, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1867, float _S1868)
{
    _d_dot_0(_S1866, _S1867, _S1868);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_14, float _s_dOut_5)
{
    float _S1869 = (*dpx_14).primal_0.x;
    float _S1870 = (*dpx_14).primal_0.y;
    float _S1871 = (*dpx_14).primal_0.z;
    DiffPair_float_0 _S1872;
    (&_S1872)->primal_0 = _S1869 * _S1869 + _S1870 * _S1870 + _S1871 * _S1871;
    (&_S1872)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1872, _s_dOut_5);
    float _S1873 = (*dpx_14).primal_0.z * _S1872.differential_0;
    float _S1874 = _S1873 + _S1873;
    float _S1875 = (*dpx_14).primal_0.y * _S1872.differential_0;
    float _S1876 = _S1875 + _S1875;
    float _S1877 = (*dpx_14).primal_0.x * _S1872.differential_0;
    float _S1878 = _S1877 + _S1877;
    float3  _S1879 = make_float3 (0.0f);
    *&((&_S1879)->z) = _S1874;
    *&((&_S1879)->y) = _S1876;
    *&((&_S1879)->x) = _S1878;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S1879;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1880, float _S1881)
{
    s_bwd_prop_length_impl_1(_S1880, _S1881);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float3  _s_dOut_6)
{
    float _S1882 = length_1((*dpx_15).primal_0);
    float3  _S1883 = (*dpx_15).primal_0 * _s_dOut_6;
    float3  _S1884 = make_float3 (1.0f / _S1882) * _s_dOut_6;
    float _S1885 = - ((_S1883.x + _S1883.y + _S1883.z) / (_S1882 * _S1882));
    float3  _S1886 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1887;
    (&_S1887)->primal_0 = (*dpx_15).primal_0;
    (&_S1887)->differential_0 = _S1886;
    s_bwd_length_impl_1(&_S1887, _S1885);
    float3  _S1888 = _S1884 + _S1887.differential_0;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S1888;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1889, float3  _S1890)
{
    s_bwd_prop_normalize_impl_0(_S1889, _S1890);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1891, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1892, float3  _S1893)
{
    _d_cross_0(_S1891, _S1892, _S1893);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1894, float _S1895)
{
    _d_exp2_0(_S1894, _S1895);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S1896, float _S1897)
{
    _d_abs_0(_S1896, _S1897);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_13, float4  quat_16, float3  scale_15, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_17, float fy_17, float cx_17, float cy_17, float4  radial_coeffs_16, float2  tangential_coeffs_16, float2  thin_prism_coeffs_16, uint image_width_13, uint image_height_13, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float v_depth_3, float2  v_out_hardness_0, float3  v_rgb_3, float3  v_normal_0, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_0(R_15, mean_13) + t_14;
    float _S1898 = scale_15.x;
    float _S1899 = s_primal_ctx_exp_1(_S1898);
    float _S1900 = scale_15.y;
    float _S1901 = s_primal_ctx_exp_1(_S1900);
    float sz_4 = scale_15.z - 0.5f * (_S1898 + _S1900);
    float _S1902 = quat_16.y;
    float _S1903 = _S1902 * _S1902 + quat_16.z * quat_16.z + quat_16.w * quat_16.w + quat_16.x * quat_16.x;
    float _S1904 = s_primal_ctx_rsqrt_0(_S1903);
    float x_53 = quat_16.y * _S1904;
    float y_32 = quat_16.z * _S1904;
    float z_29 = quat_16.w * _S1904;
    float w_16 = quat_16.x * _S1904;
    float x2_16 = x_53 * x_53;
    float y2_16 = y_32 * y_32;
    float z2_29 = z_29 * z_29;
    float xy_16 = x_53 * y_32;
    float xz_16 = x_53 * z_29;
    float yz_16 = y_32 * z_29;
    float wx_16 = w_16 * x_53;
    float wy_16 = w_16 * y_32;
    float wz_16 = w_16 * z_29;
    Matrix<float, 3, 3>  _S1905 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    float3  _S1906 = make_float3 (_S1899, 0.0f, 0.0f);
    float3  vert0_4 = s_primal_ctx_mul_0(_S1905, _S1906) + mean_13;
    float _S1907 = -0.5f + sz_4;
    float3  _S1908 = make_float3 (_S1899 * _S1907, _S1901, 0.0f);
    float3  vert1_4 = s_primal_ctx_mul_0(_S1905, _S1908) + mean_13;
    float _S1909 = -0.5f - sz_4;
    float3  _S1910 = make_float3 (_S1899 * _S1909, - _S1901, 0.0f);
    float3  vert2_4 = s_primal_ctx_mul_0(_S1905, _S1910) + mean_13;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_15, vert0_4) + t_14;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_15, vert1_4) + t_14;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_15, vert2_4) + t_14;
    float2  _S1911 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S1912 = vert0_c_4.z;
    float2  _S1913 = make_float2 (_S1912);
    float2  _S1914 = make_float2 (_S1912 * _S1912);
    float2  _S1915 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S1916 = vert1_c_4.z;
    float2  _S1917 = make_float2 (_S1916);
    float2  _S1918 = make_float2 (_S1916 * _S1916);
    float2  _S1919 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S1920 = vert2_c_4.z;
    float2  _S1921 = make_float2 (_S1920);
    float2  _S1922 = make_float2 (_S1920 * _S1920);
    float2  _S1923 = make_float2 (fx_17, fy_17);
    float2  _S1924 = make_float2 (cx_17, cy_17);
    float2  _S1925 = _S1923 * (_S1911 / make_float2 (_S1912)) + _S1924;
    float2  _S1926 = _S1923 * (_S1915 / make_float2 (_S1916)) + _S1924;
    float2  _S1927 = _S1923 * (_S1919 / make_float2 (_S1920)) + _S1924;
    float2  e0_4 = _S1926 - _S1925;
    float2  e1_4 = _S1927 - _S1926;
    float2  e2_0 = _S1925 - _S1927;
    float _S1928 = e0_4.x;
    float _S1929 = e1_4.y;
    float _S1930 = e0_4.y;
    float _S1931 = e1_4.x;
    float _S1932 = _S1928 * _S1929 - _S1930 * _S1931;
    float _S1933 = 1.0f - hardness_4.y;
    float _S1934 = -1.0f / _S1933;
    float _S1935 = _S1933 * _S1933;
    float _S1936 = _S1925.x;
    float _S1937 = _S1926.x;
    float _S1938 = s_primal_ctx_max_0(_S1936, _S1937);
    float _S1939 = _S1927.x;
    float _S1940 = s_primal_ctx_min_0(_S1936, _S1937);
    float _S1941 = _S1925.y;
    float _S1942 = _S1926.y;
    float _S1943 = s_primal_ctx_max_0(_S1941, _S1942);
    float _S1944 = _S1927.y;
    float _S1945 = s_primal_ctx_min_0(_S1941, _S1942);
    Matrix<float, 3, 3>  _S1946 = transpose_0(R_15);
    float3  _S1947 = mean_13 - - s_primal_ctx_mul_0(_S1946, t_14);
    float _S1948 = _S1947.x;
    float _S1949 = _S1947.y;
    float _S1950 = _S1947.z;
    float _S1951 = _S1948 * _S1948 + _S1949 * _S1949 + _S1950 * _S1950;
    float _S1952 = s_primal_ctx_sqrt_0(_S1951);
    float x_54 = _S1948 / _S1952;
    float3  _S1953 = make_float3 (x_54);
    float _S1954 = _S1952 * _S1952;
    float y_33 = _S1949 / _S1952;
    float z_30 = _S1950 / _S1952;
    float3  _S1955 = make_float3 (z_30);
    float _S1956 = - y_33;
    float3  _S1957 = make_float3 (_S1956);
    float z2_30 = z_30 * z_30;
    float fTmp0B_13 = -1.09254848957061768f * z_30;
    float fC1_13 = x_54 * x_54 - y_33 * y_33;
    float _S1958 = 2.0f * x_54;
    float fS1_13 = _S1958 * y_33;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S1959 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_54;
    float3  _S1960 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_33;
    float3  _S1961 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S1962 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S1963 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_30;
    float _S1964 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_30 * _S1964;
    float3  _S1965 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_54;
    float3  _S1966 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_33;
    float3  _S1967 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S1968 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S1969 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_54 * fC1_13 - y_33 * fS1_13);
    float3  _S1970 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_54 * fS1_13 + y_33 * fC1_13);
    float3  _S1971 = make_float3 (pSH9_3);
    float3  _S1972 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1956) * (*sh_coeffs_13)[int(1)] + make_float3 (z_30) * (*sh_coeffs_13)[int(2)] - make_float3 (x_54) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    float3  _S1973 = make_float3 (0.0f);
    float3  _S1974 = vert1_4 - vert0_4;
    float3  _S1975 = vert2_4 - vert0_4;
    float3  _S1976 = s_primal_ctx_cross_0(_S1974, _S1975);
    float3  _S1977 = s_primal_ctx_mul_0(R_15, _S1976);
    float3  _S1978 = normalize_0(_S1977);
    float3  _S1979 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S1978, mean_c_12)))))) * v_normal_0;
    float3  _S1980 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1981;
    (&_S1981)->primal_0 = _S1978;
    (&_S1981)->differential_0 = _S1980;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1982;
    (&_S1982)->primal_0 = mean_c_12;
    (&_S1982)->differential_0 = _S1980;
    s_bwd_prop_dot_0(&_S1981, &_S1982, 0.0f);
    float3  _S1983 = _S1979 + _S1981.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1984;
    (&_S1984)->primal_0 = _S1977;
    (&_S1984)->differential_0 = _S1980;
    s_bwd_normalize_impl_0(&_S1984, _S1983);
    Matrix<float, 3, 3>  _S1985 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1986;
    (&_S1986)->primal_0 = R_15;
    (&_S1986)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1987;
    (&_S1987)->primal_0 = _S1976;
    (&_S1987)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S1986, &_S1987, _S1984.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1988;
    (&_S1988)->primal_0 = _S1974;
    (&_S1988)->differential_0 = _S1980;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1989;
    (&_S1989)->primal_0 = _S1975;
    (&_S1989)->differential_0 = _S1980;
    s_bwd_prop_cross_0(&_S1988, &_S1989, _S1987.differential_0);
    float3  _S1990 = - _S1989.differential_0;
    float3  _S1991 = - _S1988.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1992;
    (&_S1992)->primal_0 = _S1972;
    (&_S1992)->differential_0 = _S1980;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1993;
    (&_S1993)->primal_0 = _S1973;
    (&_S1993)->differential_0 = _S1980;
    s_bwd_prop_max_0(&_S1992, &_S1993, v_rgb_3);
    float3  _S1994 = _S1970 * _S1992.differential_0;
    float3  _S1995 = (*sh_coeffs_13)[int(15)] * _S1992.differential_0;
    float3  _S1996 = _S1968 * _S1992.differential_0;
    float3  _S1997 = (*sh_coeffs_13)[int(14)] * _S1992.differential_0;
    float3  _S1998 = _S1966 * _S1992.differential_0;
    float3  _S1999 = (*sh_coeffs_13)[int(13)] * _S1992.differential_0;
    float3  _S2000 = _S1965 * _S1992.differential_0;
    float3  _S2001 = (*sh_coeffs_13)[int(12)] * _S1992.differential_0;
    float3  _S2002 = _S1967 * _S1992.differential_0;
    float3  _S2003 = (*sh_coeffs_13)[int(11)] * _S1992.differential_0;
    float3  _S2004 = _S1969 * _S1992.differential_0;
    float3  _S2005 = (*sh_coeffs_13)[int(10)] * _S1992.differential_0;
    float3  _S2006 = _S1971 * _S1992.differential_0;
    float3  _S2007 = (*sh_coeffs_13)[int(9)] * _S1992.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2007.x + _S2007.y + _S2007.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S1995.x + _S1995.y + _S1995.z);
    float _S2008 = _S2005.x + _S2005.y + _S2005.z;
    float _S2009 = _S1997.x + _S1997.y + _S1997.z;
    float _S2010 = _S2003.x + _S2003.y + _S2003.z;
    float _S2011 = _S1999.x + _S1999.y + _S1999.z;
    float _S2012 = _S2001.x + _S2001.y + _S2001.z;
    float _S2013 = - s_diff_fC2_T_3;
    float3  _S2014 = _S1962 * _S1992.differential_0;
    float3  _S2015 = (*sh_coeffs_13)[int(8)] * _S1992.differential_0;
    float3  _S2016 = _S1960 * _S1992.differential_0;
    float3  _S2017 = (*sh_coeffs_13)[int(7)] * _S1992.differential_0;
    float3  _S2018 = _S1959 * _S1992.differential_0;
    float3  _S2019 = (*sh_coeffs_13)[int(6)] * _S1992.differential_0;
    float3  _S2020 = _S1961 * _S1992.differential_0;
    float3  _S2021 = (*sh_coeffs_13)[int(5)] * _S1992.differential_0;
    float3  _S2022 = _S1963 * _S1992.differential_0;
    float3  _S2023 = (*sh_coeffs_13)[int(4)] * _S1992.differential_0;
    float _S2024 = _S2021.x + _S2021.y + _S2021.z;
    float _S2025 = _S2017.x + _S2017.y + _S2017.z;
    float _S2026 = fTmp1B_13 * _S2008 + x_54 * s_diff_fS2_T_3 + y_33 * _S2013 + 0.54627424478530884f * (_S2023.x + _S2023.y + _S2023.z);
    float _S2027 = fTmp1B_13 * _S2009 + y_33 * s_diff_fS2_T_3 + x_54 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2015.x + _S2015.y + _S2015.z);
    float _S2028 = y_33 * - _S2027;
    float _S2029 = x_54 * _S2027;
    float _S2030 = z_30 * (1.86588168144226074f * (z_30 * _S2012) + -2.28522896766662598f * (y_33 * _S2010 + x_54 * _S2011) + 0.94617468118667603f * (_S2019.x + _S2019.y + _S2019.z));
    float3  _S2031 = make_float3 (0.48860251903533936f) * _S1992.differential_0;
    float3  _S2032 = - _S2031;
    float3  _S2033 = _S1953 * _S2032;
    float3  _S2034 = (*sh_coeffs_13)[int(3)] * _S2032;
    float3  _S2035 = _S1955 * _S2031;
    float3  _S2036 = (*sh_coeffs_13)[int(2)] * _S2031;
    float3  _S2037 = _S1957 * _S2031;
    float3  _S2038 = (*sh_coeffs_13)[int(1)] * _S2031;
    float _S2039 = (_S1964 * _S2012 + 1.44530570507049561f * (fS1_13 * _S2008 + fC1_13 * _S2009) + -1.09254848957061768f * (y_33 * _S2024 + x_54 * _S2025) + _S2030 + _S2030 + _S2036.x + _S2036.y + _S2036.z) / _S1954;
    float _S2040 = _S1952 * _S2039;
    float _S2041 = (fTmp0C_13 * _S2010 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2013 + fTmp0B_13 * _S2024 + _S1958 * _S2026 + _S2028 + _S2028 + - (_S2038.x + _S2038.y + _S2038.z)) / _S1954;
    float _S2042 = _S1952 * _S2041;
    float _S2043 = (fTmp0C_13 * _S2011 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2025 + 2.0f * (y_33 * _S2026) + _S2029 + _S2029 + _S2034.x + _S2034.y + _S2034.z) / _S1954;
    float _S2044 = _S1952 * _S2043;
    float _S2045 = _S1950 * - _S2039 + _S1949 * - _S2041 + _S1948 * - _S2043;
    DiffPair_float_0 _S2046;
    (&_S2046)->primal_0 = _S1951;
    (&_S2046)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2046, _S2045);
    float _S2047 = _S1950 * _S2046.differential_0;
    float _S2048 = _S1949 * _S2046.differential_0;
    float _S2049 = _S1948 * _S2046.differential_0;
    float3  _S2050 = make_float3 (0.282094806432724f) * _S1992.differential_0;
    float3  _S2051 = make_float3 (_S2044 + _S2049 + _S2049, _S2042 + _S2048 + _S2048, _S2040 + _S2047 + _S2047);
    float3  _S2052 = - - _S2051;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2053;
    (&_S2053)->primal_0 = _S1946;
    (&_S2053)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2054;
    (&_S2054)->primal_0 = t_14;
    (&_S2054)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2053, &_S2054, _S2052);
    Matrix<float, 3, 3>  _S2055 = transpose_0(_S2053.differential_0);
    DiffPair_float_0 _S2056;
    (&_S2056)->primal_0 = _S1945;
    (&_S2056)->differential_0 = 0.0f;
    DiffPair_float_0 _S2057;
    (&_S2057)->primal_0 = _S1944;
    (&_S2057)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2056, &_S2057, 0.0f);
    DiffPair_float_0 _S2058;
    (&_S2058)->primal_0 = _S1941;
    (&_S2058)->differential_0 = 0.0f;
    DiffPair_float_0 _S2059;
    (&_S2059)->primal_0 = _S1942;
    (&_S2059)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2058, &_S2059, _S2056.differential_0);
    DiffPair_float_0 _S2060;
    (&_S2060)->primal_0 = _S1943;
    (&_S2060)->differential_0 = 0.0f;
    DiffPair_float_0 _S2061;
    (&_S2061)->primal_0 = _S1944;
    (&_S2061)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2060, &_S2061, 0.0f);
    float _S2062 = _S2057.differential_0 + _S2061.differential_0;
    DiffPair_float_0 _S2063;
    (&_S2063)->primal_0 = _S1941;
    (&_S2063)->differential_0 = 0.0f;
    DiffPair_float_0 _S2064;
    (&_S2064)->primal_0 = _S1942;
    (&_S2064)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2063, &_S2064, _S2060.differential_0);
    float _S2065 = _S2059.differential_0 + _S2064.differential_0;
    float _S2066 = _S2058.differential_0 + _S2063.differential_0;
    DiffPair_float_0 _S2067;
    (&_S2067)->primal_0 = _S1940;
    (&_S2067)->differential_0 = 0.0f;
    DiffPair_float_0 _S2068;
    (&_S2068)->primal_0 = _S1939;
    (&_S2068)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2067, &_S2068, 0.0f);
    DiffPair_float_0 _S2069;
    (&_S2069)->primal_0 = _S1936;
    (&_S2069)->differential_0 = 0.0f;
    DiffPair_float_0 _S2070;
    (&_S2070)->primal_0 = _S1937;
    (&_S2070)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2069, &_S2070, _S2067.differential_0);
    DiffPair_float_0 _S2071;
    (&_S2071)->primal_0 = _S1938;
    (&_S2071)->differential_0 = 0.0f;
    DiffPair_float_0 _S2072;
    (&_S2072)->primal_0 = _S1939;
    (&_S2072)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2071, &_S2072, 0.0f);
    float _S2073 = _S2068.differential_0 + _S2072.differential_0;
    DiffPair_float_0 _S2074;
    (&_S2074)->primal_0 = _S1936;
    (&_S2074)->differential_0 = 0.0f;
    DiffPair_float_0 _S2075;
    (&_S2075)->primal_0 = _S1937;
    (&_S2075)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2074, &_S2075, _S2071.differential_0);
    float _S2076 = _S2070.differential_0 + _S2075.differential_0;
    float _S2077 = _S2069.differential_0 + _S2074.differential_0;
    DiffPair_float_0 _S2078;
    (&_S2078)->primal_0 = _S1934;
    (&_S2078)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2078, 0.0f);
    float _S2079 = - (-1.0f * - (_S2078.differential_0 / _S1935));
    float2  _S2080 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2081;
    (&_S2081)->primal_0 = e2_0;
    (&_S2081)->differential_0 = _S2080;
    s_bwd_length_impl_0(&_S2081, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2082;
    (&_S2082)->primal_0 = e1_4;
    (&_S2082)->differential_0 = _S2080;
    s_bwd_length_impl_0(&_S2082, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2083;
    (&_S2083)->primal_0 = e0_4;
    (&_S2083)->differential_0 = _S2080;
    s_bwd_length_impl_0(&_S2083, -0.0f);
    DiffPair_float_0 _S2084;
    (&_S2084)->primal_0 = _S1932;
    (&_S2084)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2084, 0.0f);
    float _S2085 = - _S2084.differential_0;
    float2  _S2086 = _S2082.differential_0 + make_float2 (_S1930 * _S2085, _S1928 * _S2084.differential_0);
    float2  _S2087 = _S2083.differential_0 + make_float2 (_S1929 * _S2084.differential_0, _S1931 * _S2085);
    float2  _S2088 = _S1923 * (v_uv2_0 + - _S2081.differential_0 + _S2086 + make_float2 (_S2073, _S2062)) / _S1922;
    float2  _S2089 = _S1919 * - _S2088;
    float2  _S2090 = _S1921 * _S2088;
    float2  _S2091 = _S1923 * (v_uv1_0 + - _S2086 + _S2087 + make_float2 (_S2076, _S2065)) / _S1918;
    float2  _S2092 = _S1915 * - _S2091;
    float2  _S2093 = _S1917 * _S2091;
    float _S2094 = _S2092.x + _S2092.y;
    float2  _S2095 = _S1923 * (v_uv0_0 + _S2081.differential_0 + - _S2087 + make_float2 (_S2077, _S2066)) / _S1914;
    float2  _S2096 = _S1911 * - _S2095;
    float2  _S2097 = _S1913 * _S2095;
    float _S2098 = _S2096.x + _S2096.y;
    float3  s_diff_vert2_c_T_0 = make_float3 (_S2090.x, _S2090.y, _S2089.x + _S2089.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2099;
    (&_S2099)->primal_0 = R_15;
    (&_S2099)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2100;
    (&_S2100)->primal_0 = vert2_4;
    (&_S2100)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2099, &_S2100, s_diff_vert2_c_T_0);
    float3  s_diff_vert1_c_T_0 = make_float3 (_S2093.x, _S2093.y, _S2094);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2101;
    (&_S2101)->primal_0 = R_15;
    (&_S2101)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2102;
    (&_S2102)->primal_0 = vert1_4;
    (&_S2102)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2101, &_S2102, s_diff_vert1_c_T_0);
    float3  s_diff_vert0_c_T_0 = make_float3 (_S2097.x, _S2097.y, _S2098);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2103;
    (&_S2103)->primal_0 = R_15;
    (&_S2103)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2104;
    (&_S2104)->primal_0 = vert0_4;
    (&_S2104)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2103, &_S2104, s_diff_vert0_c_T_0);
    float3  _S2105 = _S1989.differential_0 + _S2100.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2106;
    (&_S2106)->primal_0 = _S1905;
    (&_S2106)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2107;
    (&_S2107)->primal_0 = _S1910;
    (&_S2107)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2106, &_S2107, _S2105);
    float _S2108 = - _S2107.differential_0.y;
    float _S2109 = _S1909 * _S2107.differential_0.x;
    float _S2110 = - (_S1899 * _S2107.differential_0.x);
    float3  _S2111 = _S1988.differential_0 + _S2102.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2112;
    (&_S2112)->primal_0 = _S1905;
    (&_S2112)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2113;
    (&_S2113)->primal_0 = _S1908;
    (&_S2113)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2112, &_S2113, _S2111);
    float _S2114 = _S1899 * _S2113.differential_0.x;
    float _S2115 = _S1907 * _S2113.differential_0.x;
    float3  _S2116 = _S1990 + _S1991 + _S2104.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2117;
    (&_S2117)->primal_0 = _S1905;
    (&_S2117)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2118;
    (&_S2118)->primal_0 = _S1906;
    (&_S2118)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2117, &_S2118, _S2116);
    Matrix<float, 3, 3>  _S2119 = transpose_0(_S2106.differential_0 + _S2112.differential_0 + _S2117.differential_0);
    float _S2120 = 2.0f * - _S2119.rows[int(2)].z;
    float _S2121 = 2.0f * _S2119.rows[int(2)].y;
    float _S2122 = 2.0f * _S2119.rows[int(2)].x;
    float _S2123 = 2.0f * _S2119.rows[int(1)].z;
    float _S2124 = 2.0f * - _S2119.rows[int(1)].y;
    float _S2125 = 2.0f * _S2119.rows[int(1)].x;
    float _S2126 = 2.0f * _S2119.rows[int(0)].z;
    float _S2127 = 2.0f * _S2119.rows[int(0)].y;
    float _S2128 = 2.0f * - _S2119.rows[int(0)].x;
    float _S2129 = - _S2125 + _S2127;
    float _S2130 = _S2122 + - _S2126;
    float _S2131 = - _S2121 + _S2123;
    float _S2132 = _S2121 + _S2123;
    float _S2133 = _S2122 + _S2126;
    float _S2134 = _S2125 + _S2127;
    float _S2135 = z_29 * (_S2124 + _S2128);
    float _S2136 = y_32 * (_S2120 + _S2128);
    float _S2137 = x_53 * (_S2120 + _S2124);
    float _S2138 = z_29 * _S2129 + y_32 * _S2130 + x_53 * _S2131;
    float _S2139 = _S1904 * _S2138;
    float _S2140 = w_16 * _S2129 + y_32 * _S2132 + x_53 * _S2133 + _S2135 + _S2135;
    float _S2141 = _S1904 * _S2140;
    float _S2142 = w_16 * _S2130 + z_29 * _S2132 + x_53 * _S2134 + _S2136 + _S2136;
    float _S2143 = _S1904 * _S2142;
    float _S2144 = w_16 * _S2131 + z_29 * _S2133 + y_32 * _S2134 + _S2137 + _S2137;
    float _S2145 = _S1904 * _S2144;
    float _S2146 = quat_16.x * _S2138 + quat_16.w * _S2140 + quat_16.z * _S2142 + quat_16.y * _S2144;
    DiffPair_float_0 _S2147;
    (&_S2147)->primal_0 = _S1903;
    (&_S2147)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S2147, _S2146);
    float _S2148 = quat_16.x * _S2147.differential_0;
    float _S2149 = quat_16.w * _S2147.differential_0;
    float _S2150 = quat_16.z * _S2147.differential_0;
    float _S2151 = quat_16.y * _S2147.differential_0;
    float _S2152 = _S2141 + _S2149 + _S2149;
    float _S2153 = _S2143 + _S2150 + _S2150;
    float _S2154 = _S2145 + _S2151 + _S2151;
    float _S2155 = _S2139 + _S2148 + _S2148;
    float _S2156 = _S2110 + _S2114;
    float _S2157 = 0.5f * - _S2156;
    float _S2158 = _S2108 + _S2113.differential_0.y;
    DiffPair_float_0 _S2159;
    (&_S2159)->primal_0 = _S1900;
    (&_S2159)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2159, _S2158);
    float _S2160 = _S2157 + _S2159.differential_0;
    float _S2161 = _S2109 + _S2115 + _S2118.differential_0.x;
    DiffPair_float_0 _S2162;
    (&_S2162)->primal_0 = _S1898;
    (&_S2162)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2162, _S2161);
    float _S2163 = _S2157 + _S2162.differential_0;
    float3  _S2164 = _S1982.differential_0 + make_float3 (0.0f, 0.0f, v_depth_3);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2165;
    (&_S2165)->primal_0 = R_15;
    (&_S2165)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2166;
    (&_S2166)->primal_0 = mean_13;
    (&_S2166)->differential_0 = _S1980;
    s_bwd_prop_mul_0(&_S2165, &_S2166, _S2164);
    float3  _S2167 = _S2054.differential_0 + s_diff_vert2_c_T_0 + s_diff_vert1_c_T_0 + s_diff_vert0_c_T_0 + _S2164;
    Matrix<float, 3, 3>  _S2168 = _S1986.differential_0 + _S2055 + _S2099.differential_0 + _S2101.differential_0 + _S2103.differential_0 + _S2165.differential_0;
    FixedArray<float3 , 16>  _S2169;
    _S2169[int(0)] = _S1980;
    _S2169[int(1)] = _S1980;
    _S2169[int(2)] = _S1980;
    _S2169[int(3)] = _S1980;
    _S2169[int(4)] = _S1980;
    _S2169[int(5)] = _S1980;
    _S2169[int(6)] = _S1980;
    _S2169[int(7)] = _S1980;
    _S2169[int(8)] = _S1980;
    _S2169[int(9)] = _S1980;
    _S2169[int(10)] = _S1980;
    _S2169[int(11)] = _S1980;
    _S2169[int(12)] = _S1980;
    _S2169[int(13)] = _S1980;
    _S2169[int(14)] = _S1980;
    _S2169[int(15)] = _S1980;
    _S2169[int(15)] = _S1994;
    _S2169[int(14)] = _S1996;
    _S2169[int(13)] = _S1998;
    _S2169[int(12)] = _S2000;
    _S2169[int(11)] = _S2002;
    _S2169[int(10)] = _S2004;
    _S2169[int(9)] = _S2006;
    _S2169[int(8)] = _S2014;
    _S2169[int(7)] = _S2016;
    _S2169[int(6)] = _S2018;
    _S2169[int(5)] = _S2020;
    _S2169[int(4)] = _S2022;
    _S2169[int(3)] = _S2033;
    _S2169[int(2)] = _S2035;
    _S2169[int(1)] = _S2037;
    _S2169[int(0)] = _S2050;
    float2  _S2170 = v_out_hardness_0 + make_float2 (0.0f, _S2079);
    float3  _S2171 = make_float3 (_S2163, _S2160, _S2156);
    float4  _S2172 = make_float4 (0.0f);
    *&((&_S2172)->w) = _S2152;
    *&((&_S2172)->z) = _S2153;
    *&((&_S2172)->y) = _S2154;
    *&((&_S2172)->x) = _S2155;
    *v_mean_3 = _S2051 + _S2105 + _S2111 + _S2116 + _S2166.differential_0;
    *v_quat_3 = _S2172;
    *v_scale_3 = _S2171;
    *v_hardness_0 = _S2170;
    *v_sh_coeffs_3 = _S2169;
    *v_R_3 = _S2168;
    *v_t_3 = _S2167;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_14, float4  quat_17, float3  scale_16, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_18, float fy_18, float cx_18, float cy_18, float4  radial_coeffs_17, float2  tangential_coeffs_17, float2  thin_prism_coeffs_17, uint image_width_14, uint image_height_14, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float v_depth_4, float2  v_out_hardness_1, float3  v_rgb_4, float3  v_normal_1, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float3  mean_c_13 = s_primal_ctx_mul_0(R_16, mean_14) + t_15;
    float _S2173 = scale_16.x;
    float _S2174 = s_primal_ctx_exp_1(_S2173);
    float _S2175 = scale_16.y;
    float _S2176 = s_primal_ctx_exp_1(_S2175);
    float sz_5 = scale_16.z - 0.5f * (_S2173 + _S2175);
    float _S2177 = quat_17.y;
    float _S2178 = _S2177 * _S2177 + quat_17.z * quat_17.z + quat_17.w * quat_17.w + quat_17.x * quat_17.x;
    float _S2179 = s_primal_ctx_rsqrt_0(_S2178);
    float x_55 = quat_17.y * _S2179;
    float y_34 = quat_17.z * _S2179;
    float z_31 = quat_17.w * _S2179;
    float w_17 = quat_17.x * _S2179;
    float x2_17 = x_55 * x_55;
    float y2_17 = y_34 * y_34;
    float z2_31 = z_31 * z_31;
    float xy_17 = x_55 * y_34;
    float xz_17 = x_55 * z_31;
    float yz_17 = y_34 * z_31;
    float wx_17 = w_17 * x_55;
    float wy_17 = w_17 * y_34;
    float wz_17 = w_17 * z_31;
    Matrix<float, 3, 3>  _S2180 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    float3  _S2181 = make_float3 (_S2174, 0.0f, 0.0f);
    float3  vert0_5 = s_primal_ctx_mul_0(_S2180, _S2181) + mean_14;
    float _S2182 = -0.5f + sz_5;
    float3  _S2183 = make_float3 (_S2174 * _S2182, _S2176, 0.0f);
    float3  vert1_5 = s_primal_ctx_mul_0(_S2180, _S2183) + mean_14;
    float _S2184 = -0.5f - sz_5;
    float3  _S2185 = make_float3 (_S2174 * _S2184, - _S2176, 0.0f);
    float3  vert2_5 = s_primal_ctx_mul_0(_S2180, _S2185) + mean_14;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_16, vert0_5) + t_15;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_16, vert1_5) + t_15;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_16, vert2_5) + t_15;
    CameraDistortion_0 _S2186 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_17, tangential_coeffs_17, thin_prism_coeffs_17);
    float2  _S2187 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S2188 = length_0(_S2187);
    float _S2189 = vert0_c_5.z;
    float _S2190 = s_primal_ctx_atan2_0(_S2188, _S2189);
    bool _S2191 = _S2190 < 0.00100000004749745f;
    float k_6;
    float _S2192;
    float _S2193;
    float _S2194;
    if(_S2191)
    {
        float _S2195 = 1.0f - _S2190 * _S2190 / 3.0f;
        float _S2196 = _S2189 * _S2189;
        k_6 = _S2195 / _S2189;
        _S2192 = 0.0f;
        _S2193 = _S2196;
        _S2194 = _S2195;
    }
    else
    {
        float _S2197 = _S2188 * _S2188;
        k_6 = _S2190 / _S2188;
        _S2192 = _S2197;
        _S2193 = 0.0f;
        _S2194 = 0.0f;
    }
    float2  _S2198 = make_float2 (k_6);
    float2  _S2199 = _S2187 * make_float2 (k_6);
    float k1_5 = _S2186.radial_coeffs_0.x;
    float k2_5 = _S2186.radial_coeffs_0.y;
    float k3_5 = _S2186.radial_coeffs_0.z;
    float k4_5 = _S2186.radial_coeffs_0.w;
    float p1_5 = _S2186.tangential_coeffs_0.x;
    float p2_5 = _S2186.tangential_coeffs_0.y;
    float sx1_5 = _S2186.thin_prism_coeffs_0.x;
    float sy1_5 = _S2186.thin_prism_coeffs_0.y;
    float u_15 = _S2199.x;
    float v_15 = _S2199.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S2200 = k3_5 + r2_15 * k4_5;
    float _S2201 = k2_5 + r2_15 * _S2200;
    float _S2202 = k1_5 + r2_15 * _S2201;
    float radial_1 = 1.0f + r2_15 * _S2202;
    float2  _S2203 = make_float2 (radial_1);
    float _S2204 = 2.0f * p1_5;
    float _S2205 = _S2204 * u_15;
    float _S2206 = 2.0f * u_15;
    float _S2207 = r2_15 + _S2206 * u_15;
    float _S2208 = 2.0f * p2_5;
    float _S2209 = _S2208 * u_15;
    float _S2210 = 2.0f * v_15;
    float _S2211 = r2_15 + _S2210 * v_15;
    float2  _S2212 = _S2199 * make_float2 (radial_1) + make_float2 (_S2205 * v_15 + p2_5 * _S2207 + sx1_5 * r2_15, _S2209 * v_15 + p1_5 * _S2211 + sy1_5 * r2_15);
    float _S2213 = fx_18 * _S2212.x + cx_18;
    float _S2214 = fy_18 * _S2212.y + cy_18;
    float2  _S2215 = make_float2 (_S2213, _S2214);
    float2  _S2216 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S2217 = length_0(_S2216);
    float _S2218 = vert1_c_5.z;
    float _S2219 = s_primal_ctx_atan2_0(_S2217, _S2218);
    bool _S2220 = _S2219 < 0.00100000004749745f;
    float _S2221;
    float _S2222;
    float _S2223;
    if(_S2220)
    {
        float _S2224 = 1.0f - _S2219 * _S2219 / 3.0f;
        float _S2225 = _S2218 * _S2218;
        k_6 = _S2224 / _S2218;
        _S2221 = 0.0f;
        _S2222 = _S2225;
        _S2223 = _S2224;
    }
    else
    {
        float _S2226 = _S2217 * _S2217;
        k_6 = _S2219 / _S2217;
        _S2221 = _S2226;
        _S2222 = 0.0f;
        _S2223 = 0.0f;
    }
    float2  _S2227 = make_float2 (k_6);
    float2  _S2228 = _S2216 * make_float2 (k_6);
    float u_16 = _S2228.x;
    float v_16 = _S2228.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S2229 = k3_5 + r2_16 * k4_5;
    float _S2230 = k2_5 + r2_16 * _S2229;
    float _S2231 = k1_5 + r2_16 * _S2230;
    float radial_2 = 1.0f + r2_16 * _S2231;
    float2  _S2232 = make_float2 (radial_2);
    float _S2233 = _S2204 * u_16;
    float _S2234 = 2.0f * u_16;
    float _S2235 = r2_16 + _S2234 * u_16;
    float _S2236 = _S2208 * u_16;
    float _S2237 = 2.0f * v_16;
    float _S2238 = r2_16 + _S2237 * v_16;
    float2  _S2239 = _S2228 * make_float2 (radial_2) + make_float2 (_S2233 * v_16 + p2_5 * _S2235 + sx1_5 * r2_16, _S2236 * v_16 + p1_5 * _S2238 + sy1_5 * r2_16);
    float _S2240 = fx_18 * _S2239.x + cx_18;
    float _S2241 = fy_18 * _S2239.y + cy_18;
    float2  _S2242 = make_float2 (_S2240, _S2241);
    float2  _S2243 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S2244 = length_0(_S2243);
    float _S2245 = vert2_c_5.z;
    float _S2246 = s_primal_ctx_atan2_0(_S2244, _S2245);
    bool _S2247 = _S2246 < 0.00100000004749745f;
    float _S2248;
    float _S2249;
    float _S2250;
    if(_S2247)
    {
        float _S2251 = 1.0f - _S2246 * _S2246 / 3.0f;
        float _S2252 = _S2245 * _S2245;
        k_6 = _S2251 / _S2245;
        _S2248 = 0.0f;
        _S2249 = _S2252;
        _S2250 = _S2251;
    }
    else
    {
        float _S2253 = _S2244 * _S2244;
        k_6 = _S2246 / _S2244;
        _S2248 = _S2253;
        _S2249 = 0.0f;
        _S2250 = 0.0f;
    }
    float2  _S2254 = make_float2 (k_6);
    float2  _S2255 = _S2243 * make_float2 (k_6);
    float u_17 = _S2255.x;
    float v_17 = _S2255.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S2256 = k3_5 + r2_17 * k4_5;
    float _S2257 = k2_5 + r2_17 * _S2256;
    float _S2258 = k1_5 + r2_17 * _S2257;
    float radial_3 = 1.0f + r2_17 * _S2258;
    float2  _S2259 = make_float2 (radial_3);
    float _S2260 = _S2204 * u_17;
    float _S2261 = 2.0f * u_17;
    float _S2262 = r2_17 + _S2261 * u_17;
    float _S2263 = _S2208 * u_17;
    float _S2264 = 2.0f * v_17;
    float _S2265 = r2_17 + _S2264 * v_17;
    float2  _S2266 = _S2255 * make_float2 (radial_3) + make_float2 (_S2260 * v_17 + p2_5 * _S2262 + sx1_5 * r2_17, _S2263 * v_17 + p1_5 * _S2265 + sy1_5 * r2_17);
    float _S2267 = fx_18 * _S2266.x + cx_18;
    float _S2268 = fy_18 * _S2266.y + cy_18;
    float2  _S2269 = make_float2 (_S2267, _S2268);
    float2  e0_5 = _S2242 - _S2215;
    float2  e1_5 = _S2269 - _S2242;
    float2  e2_1 = _S2215 - _S2269;
    float _S2270 = e0_5.x;
    float _S2271 = e1_5.y;
    float _S2272 = e0_5.y;
    float _S2273 = e1_5.x;
    float _S2274 = _S2270 * _S2271 - _S2272 * _S2273;
    float _S2275 = 1.0f - hardness_5.y;
    float _S2276 = -1.0f / _S2275;
    float _S2277 = _S2275 * _S2275;
    float _S2278 = s_primal_ctx_max_0(_S2213, _S2240);
    float _S2279 = s_primal_ctx_min_0(_S2213, _S2240);
    float _S2280 = s_primal_ctx_max_0(_S2214, _S2241);
    float _S2281 = s_primal_ctx_min_0(_S2214, _S2241);
    Matrix<float, 3, 3>  _S2282 = transpose_0(R_16);
    float3  _S2283 = mean_14 - - s_primal_ctx_mul_0(_S2282, t_15);
    float _S2284 = _S2283.x;
    float _S2285 = _S2283.y;
    float _S2286 = _S2283.z;
    float _S2287 = _S2284 * _S2284 + _S2285 * _S2285 + _S2286 * _S2286;
    float _S2288 = s_primal_ctx_sqrt_0(_S2287);
    float x_56 = _S2284 / _S2288;
    float3  _S2289 = make_float3 (x_56);
    float _S2290 = _S2288 * _S2288;
    float y_35 = _S2285 / _S2288;
    float z_32 = _S2286 / _S2288;
    float3  _S2291 = make_float3 (z_32);
    float _S2292 = - y_35;
    float3  _S2293 = make_float3 (_S2292);
    float z2_32 = z_32 * z_32;
    float fTmp0B_14 = -1.09254848957061768f * z_32;
    float fC1_14 = x_56 * x_56 - y_35 * y_35;
    float _S2294 = 2.0f * x_56;
    float fS1_14 = _S2294 * y_35;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2295 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_56;
    float3  _S2296 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_35;
    float3  _S2297 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2298 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2299 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_32;
    float _S2300 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_32 * _S2300;
    float3  _S2301 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_56;
    float3  _S2302 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_35;
    float3  _S2303 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2304 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2305 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_56 * fC1_14 - y_35 * fS1_14);
    float3  _S2306 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_56 * fS1_14 + y_35 * fC1_14);
    float3  _S2307 = make_float3 (pSH9_4);
    float3  _S2308 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2292) * (*sh_coeffs_14)[int(1)] + make_float3 (z_32) * (*sh_coeffs_14)[int(2)] - make_float3 (x_56) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    float3  _S2309 = make_float3 (0.0f);
    float3  _S2310 = vert1_5 - vert0_5;
    float3  _S2311 = vert2_5 - vert0_5;
    float3  _S2312 = s_primal_ctx_cross_0(_S2310, _S2311);
    float3  _S2313 = s_primal_ctx_mul_0(R_16, _S2312);
    float3  _S2314 = normalize_0(_S2313);
    float3  _S2315 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S2314, mean_c_13)))))) * v_normal_1;
    float3  _S2316 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2317;
    (&_S2317)->primal_0 = _S2314;
    (&_S2317)->differential_0 = _S2316;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2318;
    (&_S2318)->primal_0 = mean_c_13;
    (&_S2318)->differential_0 = _S2316;
    s_bwd_prop_dot_0(&_S2317, &_S2318, 0.0f);
    float3  _S2319 = _S2315 + _S2317.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2320;
    (&_S2320)->primal_0 = _S2313;
    (&_S2320)->differential_0 = _S2316;
    s_bwd_normalize_impl_0(&_S2320, _S2319);
    Matrix<float, 3, 3>  _S2321 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2322;
    (&_S2322)->primal_0 = R_16;
    (&_S2322)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2323;
    (&_S2323)->primal_0 = _S2312;
    (&_S2323)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2322, &_S2323, _S2320.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2324;
    (&_S2324)->primal_0 = _S2310;
    (&_S2324)->differential_0 = _S2316;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2325;
    (&_S2325)->primal_0 = _S2311;
    (&_S2325)->differential_0 = _S2316;
    s_bwd_prop_cross_0(&_S2324, &_S2325, _S2323.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2326 = _S2324;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2327 = _S2325;
    float3  _S2328 = - _S2325.differential_0;
    float3  _S2329 = - _S2324.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2330;
    (&_S2330)->primal_0 = _S2308;
    (&_S2330)->differential_0 = _S2316;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2331;
    (&_S2331)->primal_0 = _S2309;
    (&_S2331)->differential_0 = _S2316;
    s_bwd_prop_max_0(&_S2330, &_S2331, v_rgb_4);
    float3  _S2332 = _S2306 * _S2330.differential_0;
    float3  _S2333 = (*sh_coeffs_14)[int(15)] * _S2330.differential_0;
    float3  _S2334 = _S2304 * _S2330.differential_0;
    float3  _S2335 = (*sh_coeffs_14)[int(14)] * _S2330.differential_0;
    float3  _S2336 = _S2302 * _S2330.differential_0;
    float3  _S2337 = (*sh_coeffs_14)[int(13)] * _S2330.differential_0;
    float3  _S2338 = _S2301 * _S2330.differential_0;
    float3  _S2339 = (*sh_coeffs_14)[int(12)] * _S2330.differential_0;
    float3  _S2340 = _S2303 * _S2330.differential_0;
    float3  _S2341 = (*sh_coeffs_14)[int(11)] * _S2330.differential_0;
    float3  _S2342 = _S2305 * _S2330.differential_0;
    float3  _S2343 = (*sh_coeffs_14)[int(10)] * _S2330.differential_0;
    float3  _S2344 = _S2307 * _S2330.differential_0;
    float3  _S2345 = (*sh_coeffs_14)[int(9)] * _S2330.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2345.x + _S2345.y + _S2345.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2333.x + _S2333.y + _S2333.z);
    float _S2346 = _S2343.x + _S2343.y + _S2343.z;
    float _S2347 = _S2335.x + _S2335.y + _S2335.z;
    float _S2348 = _S2341.x + _S2341.y + _S2341.z;
    float _S2349 = _S2337.x + _S2337.y + _S2337.z;
    float _S2350 = _S2339.x + _S2339.y + _S2339.z;
    float _S2351 = - s_diff_fC2_T_4;
    float3  _S2352 = _S2298 * _S2330.differential_0;
    float3  _S2353 = (*sh_coeffs_14)[int(8)] * _S2330.differential_0;
    float3  _S2354 = _S2296 * _S2330.differential_0;
    float3  _S2355 = (*sh_coeffs_14)[int(7)] * _S2330.differential_0;
    float3  _S2356 = _S2295 * _S2330.differential_0;
    float3  _S2357 = (*sh_coeffs_14)[int(6)] * _S2330.differential_0;
    float3  _S2358 = _S2297 * _S2330.differential_0;
    float3  _S2359 = (*sh_coeffs_14)[int(5)] * _S2330.differential_0;
    float3  _S2360 = _S2299 * _S2330.differential_0;
    float3  _S2361 = (*sh_coeffs_14)[int(4)] * _S2330.differential_0;
    float _S2362 = _S2359.x + _S2359.y + _S2359.z;
    float _S2363 = _S2355.x + _S2355.y + _S2355.z;
    float _S2364 = fTmp1B_14 * _S2346 + x_56 * s_diff_fS2_T_4 + y_35 * _S2351 + 0.54627424478530884f * (_S2361.x + _S2361.y + _S2361.z);
    float _S2365 = fTmp1B_14 * _S2347 + y_35 * s_diff_fS2_T_4 + x_56 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2353.x + _S2353.y + _S2353.z);
    float _S2366 = y_35 * - _S2365;
    float _S2367 = x_56 * _S2365;
    float _S2368 = z_32 * (1.86588168144226074f * (z_32 * _S2350) + -2.28522896766662598f * (y_35 * _S2348 + x_56 * _S2349) + 0.94617468118667603f * (_S2357.x + _S2357.y + _S2357.z));
    float3  _S2369 = make_float3 (0.48860251903533936f) * _S2330.differential_0;
    float3  _S2370 = - _S2369;
    float3  _S2371 = _S2289 * _S2370;
    float3  _S2372 = (*sh_coeffs_14)[int(3)] * _S2370;
    float3  _S2373 = _S2291 * _S2369;
    float3  _S2374 = (*sh_coeffs_14)[int(2)] * _S2369;
    float3  _S2375 = _S2293 * _S2369;
    float3  _S2376 = (*sh_coeffs_14)[int(1)] * _S2369;
    float _S2377 = (_S2300 * _S2350 + 1.44530570507049561f * (fS1_14 * _S2346 + fC1_14 * _S2347) + -1.09254848957061768f * (y_35 * _S2362 + x_56 * _S2363) + _S2368 + _S2368 + _S2374.x + _S2374.y + _S2374.z) / _S2290;
    float _S2378 = _S2288 * _S2377;
    float _S2379 = (fTmp0C_14 * _S2348 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2351 + fTmp0B_14 * _S2362 + _S2294 * _S2364 + _S2366 + _S2366 + - (_S2376.x + _S2376.y + _S2376.z)) / _S2290;
    float _S2380 = _S2288 * _S2379;
    float _S2381 = (fTmp0C_14 * _S2349 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2363 + 2.0f * (y_35 * _S2364) + _S2367 + _S2367 + _S2372.x + _S2372.y + _S2372.z) / _S2290;
    float _S2382 = _S2288 * _S2381;
    float _S2383 = _S2286 * - _S2377 + _S2285 * - _S2379 + _S2284 * - _S2381;
    DiffPair_float_0 _S2384;
    (&_S2384)->primal_0 = _S2287;
    (&_S2384)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2384, _S2383);
    float _S2385 = _S2286 * _S2384.differential_0;
    float _S2386 = _S2285 * _S2384.differential_0;
    float _S2387 = _S2284 * _S2384.differential_0;
    float3  _S2388 = make_float3 (0.282094806432724f) * _S2330.differential_0;
    float3  _S2389 = make_float3 (_S2382 + _S2387 + _S2387, _S2380 + _S2386 + _S2386, _S2378 + _S2385 + _S2385);
    float3  _S2390 = - - _S2389;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2391;
    (&_S2391)->primal_0 = _S2282;
    (&_S2391)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2392;
    (&_S2392)->primal_0 = t_15;
    (&_S2392)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2391, &_S2392, _S2390);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2393 = _S2392;
    Matrix<float, 3, 3>  _S2394 = transpose_0(_S2391.differential_0);
    DiffPair_float_0 _S2395;
    (&_S2395)->primal_0 = _S2281;
    (&_S2395)->differential_0 = 0.0f;
    DiffPair_float_0 _S2396;
    (&_S2396)->primal_0 = _S2268;
    (&_S2396)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2395, &_S2396, 0.0f);
    DiffPair_float_0 _S2397;
    (&_S2397)->primal_0 = _S2214;
    (&_S2397)->differential_0 = 0.0f;
    DiffPair_float_0 _S2398;
    (&_S2398)->primal_0 = _S2241;
    (&_S2398)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2397, &_S2398, _S2395.differential_0);
    DiffPair_float_0 _S2399;
    (&_S2399)->primal_0 = _S2280;
    (&_S2399)->differential_0 = 0.0f;
    DiffPair_float_0 _S2400;
    (&_S2400)->primal_0 = _S2268;
    (&_S2400)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2399, &_S2400, 0.0f);
    DiffPair_float_0 _S2401;
    (&_S2401)->primal_0 = _S2214;
    (&_S2401)->differential_0 = 0.0f;
    DiffPair_float_0 _S2402;
    (&_S2402)->primal_0 = _S2241;
    (&_S2402)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2401, &_S2402, _S2399.differential_0);
    DiffPair_float_0 _S2403;
    (&_S2403)->primal_0 = _S2279;
    (&_S2403)->differential_0 = 0.0f;
    DiffPair_float_0 _S2404;
    (&_S2404)->primal_0 = _S2267;
    (&_S2404)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2403, &_S2404, 0.0f);
    DiffPair_float_0 _S2405;
    (&_S2405)->primal_0 = _S2213;
    (&_S2405)->differential_0 = 0.0f;
    DiffPair_float_0 _S2406;
    (&_S2406)->primal_0 = _S2240;
    (&_S2406)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2405, &_S2406, _S2403.differential_0);
    DiffPair_float_0 _S2407;
    (&_S2407)->primal_0 = _S2278;
    (&_S2407)->differential_0 = 0.0f;
    DiffPair_float_0 _S2408;
    (&_S2408)->primal_0 = _S2267;
    (&_S2408)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2407, &_S2408, 0.0f);
    DiffPair_float_0 _S2409;
    (&_S2409)->primal_0 = _S2213;
    (&_S2409)->differential_0 = 0.0f;
    DiffPair_float_0 _S2410;
    (&_S2410)->primal_0 = _S2240;
    (&_S2410)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2409, &_S2410, _S2407.differential_0);
    DiffPair_float_0 _S2411;
    (&_S2411)->primal_0 = _S2276;
    (&_S2411)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2411, 0.0f);
    float _S2412 = - (-1.0f * - (_S2411.differential_0 / _S2277));
    float2  _S2413 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2414;
    (&_S2414)->primal_0 = e2_1;
    (&_S2414)->differential_0 = _S2413;
    s_bwd_length_impl_0(&_S2414, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2415;
    (&_S2415)->primal_0 = e1_5;
    (&_S2415)->differential_0 = _S2413;
    s_bwd_length_impl_0(&_S2415, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2416;
    (&_S2416)->primal_0 = e0_5;
    (&_S2416)->differential_0 = _S2413;
    s_bwd_length_impl_0(&_S2416, -0.0f);
    DiffPair_float_0 _S2417;
    (&_S2417)->primal_0 = _S2274;
    (&_S2417)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2417, 0.0f);
    float _S2418 = - _S2417.differential_0;
    float2  _S2419 = _S2415.differential_0 + make_float2 (_S2272 * _S2418, _S2270 * _S2417.differential_0);
    float2  _S2420 = _S2416.differential_0 + make_float2 (_S2271 * _S2417.differential_0, _S2273 * _S2418);
    float2  _S2421 = v_uv2_1 + - _S2414.differential_0 + _S2419;
    float _S2422 = fy_18 * (_S2396.differential_0 + _S2400.differential_0 + _S2421.y);
    float _S2423 = fx_18 * (_S2404.differential_0 + _S2408.differential_0 + _S2421.x);
    float2  _S2424 = make_float2 (_S2423, _S2422);
    float2  _S2425 = _S2255 * _S2424;
    float2  _S2426 = _S2259 * _S2424;
    float _S2427 = r2_17 * _S2422;
    float _S2428 = p1_5 * _S2422;
    float _S2429 = _S2265 * _S2422;
    float _S2430 = v_17 * _S2422;
    float _S2431 = u_17 * _S2430;
    float _S2432 = r2_17 * _S2423;
    float _S2433 = p2_5 * _S2423;
    float _S2434 = _S2262 * _S2423;
    float _S2435 = v_17 * _S2423;
    float _S2436 = u_17 * _S2435;
    float _S2437 = _S2425.x + _S2425.y;
    float _S2438 = r2_17 * _S2437;
    float _S2439 = r2_17 * _S2438;
    float _S2440 = r2_17 * _S2439;
    float _S2441 = r2_17 * _S2440;
    float _S2442 = sy1_5 * _S2422 + _S2428 + sx1_5 * _S2423 + _S2433 + _S2258 * _S2437 + _S2257 * _S2438 + _S2256 * _S2439 + k4_5 * _S2440;
    float _S2443 = v_17 * _S2442;
    float _S2444 = u_17 * _S2442;
    float _S2445 = _S2264 * _S2428 + 2.0f * (v_17 * _S2428) + _S2263 * _S2422 + _S2260 * _S2423 + _S2443 + _S2443;
    float _S2446 = _S2208 * _S2430 + _S2261 * _S2433 + 2.0f * (u_17 * _S2433) + _S2204 * _S2435 + _S2444 + _S2444;
    float2  _S2447 = v_uv0_1 + _S2414.differential_0 + - _S2420;
    float2  _S2448 = v_out_hardness_1 + make_float2 (0.0f, _S2412);
    float _S2449 = _S2406.differential_0 + _S2410.differential_0;
    float2  _S2450 = v_uv1_1 + - _S2419 + _S2420;
    float3  _S2451 = _S2318.differential_0 + make_float3 (0.0f, 0.0f, v_depth_4);
    Matrix<float, 3, 3>  _S2452 = _S2322.differential_0 + _S2394;
    float3  _S2453 = _S2328 + _S2329;
    FixedArray<float3 , 16>  _S2454;
    _S2454[int(0)] = _S2316;
    _S2454[int(1)] = _S2316;
    _S2454[int(2)] = _S2316;
    _S2454[int(3)] = _S2316;
    _S2454[int(4)] = _S2316;
    _S2454[int(5)] = _S2316;
    _S2454[int(6)] = _S2316;
    _S2454[int(7)] = _S2316;
    _S2454[int(8)] = _S2316;
    _S2454[int(9)] = _S2316;
    _S2454[int(10)] = _S2316;
    _S2454[int(11)] = _S2316;
    _S2454[int(12)] = _S2316;
    _S2454[int(13)] = _S2316;
    _S2454[int(14)] = _S2316;
    _S2454[int(15)] = _S2316;
    _S2454[int(7)] = _S2354;
    _S2454[int(0)] = _S2388;
    _S2454[int(1)] = _S2375;
    _S2454[int(2)] = _S2373;
    _S2454[int(3)] = _S2371;
    _S2454[int(4)] = _S2360;
    _S2454[int(5)] = _S2358;
    _S2454[int(6)] = _S2356;
    _S2454[int(15)] = _S2332;
    _S2454[int(8)] = _S2352;
    _S2454[int(9)] = _S2344;
    _S2454[int(10)] = _S2342;
    _S2454[int(11)] = _S2340;
    _S2454[int(12)] = _S2338;
    _S2454[int(13)] = _S2336;
    _S2454[int(14)] = _S2334;
    float3  _S2455 = _S2454[int(0)];
    float3  _S2456 = _S2454[int(1)];
    float3  _S2457 = _S2454[int(2)];
    float3  _S2458 = _S2454[int(3)];
    float3  _S2459 = _S2454[int(4)];
    float3  _S2460 = _S2454[int(5)];
    float3  _S2461 = _S2454[int(6)];
    float3  _S2462 = _S2454[int(7)];
    float3  _S2463 = _S2454[int(8)];
    float3  _S2464 = _S2454[int(9)];
    float3  _S2465 = _S2454[int(10)];
    float3  _S2466 = _S2454[int(11)];
    float3  _S2467 = _S2454[int(12)];
    float3  _S2468 = _S2454[int(13)];
    float3  _S2469 = _S2454[int(14)];
    float3  _S2470 = _S2454[int(15)];
    float _S2471 = _S2405.differential_0 + _S2409.differential_0;
    float _S2472 = _S2397.differential_0 + _S2401.differential_0;
    float _S2473 = _S2398.differential_0 + _S2402.differential_0;
    float2  _S2474 = _S2426 + make_float2 (_S2446, _S2445);
    float2  _S2475 = _S2243 * _S2474;
    float2  _S2476 = _S2254 * _S2474;
    float _S2477 = _S2475.x + _S2475.y;
    if(_S2247)
    {
        float _S2478 = _S2477 / _S2249;
        float _S2479 = _S2250 * - _S2478;
        float _S2480 = _S2246 * (0.3333333432674408f * - (_S2245 * _S2478));
        k_6 = _S2480 + _S2480;
        _S2248 = _S2479;
        _S2249 = 0.0f;
    }
    else
    {
        float _S2481 = _S2477 / _S2248;
        float _S2482 = _S2246 * - _S2481;
        k_6 = _S2244 * _S2481;
        _S2248 = 0.0f;
        _S2249 = _S2482;
    }
    DiffPair_float_0 _S2483;
    (&_S2483)->primal_0 = _S2244;
    (&_S2483)->differential_0 = 0.0f;
    DiffPair_float_0 _S2484;
    (&_S2484)->primal_0 = _S2245;
    (&_S2484)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2483, &_S2484, k_6);
    float _S2485 = _S2484.differential_0 + _S2248;
    float _S2486 = _S2483.differential_0 + _S2249;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2487;
    (&_S2487)->primal_0 = _S2243;
    (&_S2487)->differential_0 = _S2413;
    s_bwd_length_impl_0(&_S2487, _S2486);
    float2  _S2488 = _S2487.differential_0 + _S2476;
    float _S2489 = fy_18 * (_S2450.y + _S2473);
    float _S2490 = fx_18 * (_S2450.x + _S2449);
    float2  _S2491 = make_float2 (_S2490, _S2489);
    float2  _S2492 = _S2228 * _S2491;
    float _S2493 = p1_5 * _S2489;
    float _S2494 = v_16 * _S2489;
    float _S2495 = p2_5 * _S2490;
    float _S2496 = v_16 * _S2490;
    float _S2497 = _S2492.x + _S2492.y;
    float _S2498 = r2_16 * _S2497;
    float _S2499 = r2_16 * _S2498;
    float _S2500 = r2_16 * _S2499;
    float _S2501 = sy1_5 * _S2489 + _S2493 + sx1_5 * _S2490 + _S2495 + _S2231 * _S2497 + _S2230 * _S2498 + _S2229 * _S2499 + k4_5 * _S2500;
    float _S2502 = v_16 * _S2501;
    float _S2503 = u_16 * _S2501;
    float3  _S2504 = make_float3 (_S2488.x, _S2488.y, _S2485);
    float2  _S2505 = _S2232 * _S2491 + make_float2 (_S2208 * _S2494 + _S2234 * _S2495 + 2.0f * (u_16 * _S2495) + _S2204 * _S2496 + _S2503 + _S2503, _S2237 * _S2493 + 2.0f * (v_16 * _S2493) + _S2236 * _S2489 + _S2233 * _S2490 + _S2502 + _S2502);
    float _S2506 = u_16 * _S2494 + _S2431;
    float _S2507 = u_16 * _S2496 + _S2436;
    float _S2508 = r2_16 * _S2489 + _S2427;
    float _S2509 = r2_16 * _S2500 + _S2441;
    float _S2510 = _S2500 + _S2440;
    float _S2511 = _S2238 * _S2489 + _S2429;
    float _S2512 = _S2499 + _S2439;
    float _S2513 = r2_16 * _S2490 + _S2432;
    float _S2514 = _S2235 * _S2490 + _S2434;
    float _S2515 = _S2498 + _S2438;
    float2  _S2516 = _S2216 * _S2505;
    float2  _S2517 = _S2227 * _S2505;
    float _S2518 = _S2516.x + _S2516.y;
    if(_S2220)
    {
        float _S2519 = _S2518 / _S2222;
        float _S2520 = _S2223 * - _S2519;
        float _S2521 = _S2219 * (0.3333333432674408f * - (_S2218 * _S2519));
        k_6 = _S2521 + _S2521;
        _S2221 = _S2520;
        _S2222 = 0.0f;
    }
    else
    {
        float _S2522 = _S2518 / _S2221;
        float _S2523 = _S2219 * - _S2522;
        k_6 = _S2217 * _S2522;
        _S2221 = 0.0f;
        _S2222 = _S2523;
    }
    DiffPair_float_0 _S2524;
    (&_S2524)->primal_0 = _S2217;
    (&_S2524)->differential_0 = 0.0f;
    DiffPair_float_0 _S2525;
    (&_S2525)->primal_0 = _S2218;
    (&_S2525)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2524, &_S2525, k_6);
    float _S2526 = _S2525.differential_0 + _S2221;
    float _S2527 = _S2524.differential_0 + _S2222;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2528;
    (&_S2528)->primal_0 = _S2216;
    (&_S2528)->differential_0 = _S2413;
    s_bwd_length_impl_0(&_S2528, _S2527);
    float2  _S2529 = _S2528.differential_0 + _S2517;
    float _S2530 = fy_18 * (_S2447.y + _S2472);
    float _S2531 = fx_18 * (_S2447.x + _S2471);
    float2  _S2532 = make_float2 (_S2531, _S2530);
    float2  _S2533 = _S2199 * _S2532;
    float _S2534 = p1_5 * _S2530;
    float _S2535 = v_15 * _S2530;
    float _S2536 = p2_5 * _S2531;
    float _S2537 = v_15 * _S2531;
    float _S2538 = _S2533.x + _S2533.y;
    float _S2539 = r2_15 * _S2538;
    float _S2540 = r2_15 * _S2539;
    float _S2541 = r2_15 * _S2540;
    float _S2542 = sy1_5 * _S2530 + _S2534 + sx1_5 * _S2531 + _S2536 + _S2202 * _S2538 + _S2201 * _S2539 + _S2200 * _S2540 + k4_5 * _S2541;
    float _S2543 = v_15 * _S2542;
    float _S2544 = u_15 * _S2542;
    float2  _S2545 = make_float2 (r2_15 * _S2531 + _S2513, r2_15 * _S2530 + _S2508);
    float2  _S2546 = make_float2 (_S2211 * _S2530 + 2.0f * (u_15 * _S2537 + _S2507) + _S2511, 2.0f * (u_15 * _S2535 + _S2506) + _S2207 * _S2531 + _S2514);
    float4  _S2547 = make_float4 (_S2539 + _S2515, _S2540 + _S2512, _S2541 + _S2510, r2_15 * _S2541 + _S2509);
    float3  _S2548 = make_float3 (_S2529.x, _S2529.y, _S2526);
    float2  _S2549 = _S2203 * _S2532 + make_float2 (_S2208 * _S2535 + _S2206 * _S2536 + 2.0f * (u_15 * _S2536) + _S2204 * _S2537 + _S2544 + _S2544, _S2210 * _S2534 + 2.0f * (v_15 * _S2534) + _S2209 * _S2530 + _S2205 * _S2531 + _S2543 + _S2543);
    CameraDistortion_0 _S2550 = CameraDistortion_x24_syn_dzero_0();
    (&_S2550)->thin_prism_coeffs_0 = _S2545;
    (&_S2550)->tangential_coeffs_0 = _S2546;
    (&_S2550)->radial_coeffs_0 = _S2547;
    float2  _S2551 = _S2187 * _S2549;
    float2  _S2552 = _S2198 * _S2549;
    float _S2553 = _S2551.x + _S2551.y;
    if(_S2191)
    {
        float _S2554 = _S2553 / _S2193;
        float _S2555 = _S2194 * - _S2554;
        float _S2556 = _S2190 * (0.3333333432674408f * - (_S2189 * _S2554));
        k_6 = _S2556 + _S2556;
        _S2192 = _S2555;
        _S2193 = 0.0f;
    }
    else
    {
        float _S2557 = _S2553 / _S2192;
        float _S2558 = _S2190 * - _S2557;
        k_6 = _S2188 * _S2557;
        _S2192 = 0.0f;
        _S2193 = _S2558;
    }
    DiffPair_float_0 _S2559;
    (&_S2559)->primal_0 = _S2188;
    (&_S2559)->differential_0 = 0.0f;
    DiffPair_float_0 _S2560;
    (&_S2560)->primal_0 = _S2189;
    (&_S2560)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2559, &_S2560, k_6);
    float _S2561 = _S2560.differential_0 + _S2192;
    float _S2562 = _S2559.differential_0 + _S2193;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2563;
    (&_S2563)->primal_0 = _S2187;
    (&_S2563)->differential_0 = _S2413;
    s_bwd_length_impl_0(&_S2563, _S2562);
    float2  _S2564 = _S2563.differential_0 + _S2552;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2565;
    (&_S2565)->primal_0 = R_16;
    (&_S2565)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2566;
    (&_S2566)->primal_0 = vert2_5;
    (&_S2566)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2565, &_S2566, _S2504);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2567;
    (&_S2567)->primal_0 = R_16;
    (&_S2567)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2568;
    (&_S2568)->primal_0 = vert1_5;
    (&_S2568)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2567, &_S2568, _S2548);
    float3  s_diff_vert0_c_T_1 = make_float3 (_S2564.x, _S2564.y, _S2561);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2569;
    (&_S2569)->primal_0 = R_16;
    (&_S2569)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2570;
    (&_S2570)->primal_0 = vert0_5;
    (&_S2570)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2569, &_S2570, s_diff_vert0_c_T_1);
    float3  _S2571 = _S2566.differential_0 + _S2327.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2572;
    (&_S2572)->primal_0 = _S2180;
    (&_S2572)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2573;
    (&_S2573)->primal_0 = _S2185;
    (&_S2573)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2572, &_S2573, _S2571);
    float _S2574 = - _S2573.differential_0.y;
    float _S2575 = _S2184 * _S2573.differential_0.x;
    float _S2576 = - (_S2174 * _S2573.differential_0.x);
    float3  _S2577 = _S2568.differential_0 + _S2326.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2578;
    (&_S2578)->primal_0 = _S2180;
    (&_S2578)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2579;
    (&_S2579)->primal_0 = _S2183;
    (&_S2579)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2578, &_S2579, _S2577);
    float _S2580 = _S2174 * _S2579.differential_0.x;
    float _S2581 = _S2182 * _S2579.differential_0.x;
    float3  _S2582 = _S2570.differential_0 + _S2453;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2583;
    (&_S2583)->primal_0 = _S2180;
    (&_S2583)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2584;
    (&_S2584)->primal_0 = _S2181;
    (&_S2584)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2583, &_S2584, _S2582);
    Matrix<float, 3, 3>  _S2585 = transpose_0(_S2572.differential_0 + _S2578.differential_0 + _S2583.differential_0);
    float _S2586 = 2.0f * - _S2585.rows[int(2)].z;
    float _S2587 = 2.0f * _S2585.rows[int(2)].y;
    float _S2588 = 2.0f * _S2585.rows[int(2)].x;
    float _S2589 = 2.0f * _S2585.rows[int(1)].z;
    float _S2590 = 2.0f * - _S2585.rows[int(1)].y;
    float _S2591 = 2.0f * _S2585.rows[int(1)].x;
    float _S2592 = 2.0f * _S2585.rows[int(0)].z;
    float _S2593 = 2.0f * _S2585.rows[int(0)].y;
    float _S2594 = 2.0f * - _S2585.rows[int(0)].x;
    float _S2595 = - _S2591 + _S2593;
    float _S2596 = _S2588 + - _S2592;
    float _S2597 = - _S2587 + _S2589;
    float _S2598 = _S2587 + _S2589;
    float _S2599 = _S2588 + _S2592;
    float _S2600 = _S2591 + _S2593;
    float _S2601 = z_31 * (_S2590 + _S2594);
    float _S2602 = y_34 * (_S2586 + _S2594);
    float _S2603 = x_55 * (_S2586 + _S2590);
    float _S2604 = z_31 * _S2595 + y_34 * _S2596 + x_55 * _S2597;
    float _S2605 = _S2179 * _S2604;
    float _S2606 = w_17 * _S2595 + y_34 * _S2598 + x_55 * _S2599 + _S2601 + _S2601;
    float _S2607 = _S2179 * _S2606;
    float _S2608 = w_17 * _S2596 + z_31 * _S2598 + x_55 * _S2600 + _S2602 + _S2602;
    float _S2609 = _S2179 * _S2608;
    float _S2610 = w_17 * _S2597 + z_31 * _S2599 + y_34 * _S2600 + _S2603 + _S2603;
    float _S2611 = _S2179 * _S2610;
    float _S2612 = quat_17.x * _S2604 + quat_17.w * _S2606 + quat_17.z * _S2608 + quat_17.y * _S2610;
    DiffPair_float_0 _S2613;
    (&_S2613)->primal_0 = _S2178;
    (&_S2613)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S2613, _S2612);
    float _S2614 = quat_17.x * _S2613.differential_0;
    float _S2615 = quat_17.w * _S2613.differential_0;
    float _S2616 = quat_17.z * _S2613.differential_0;
    float _S2617 = quat_17.y * _S2613.differential_0;
    float _S2618 = _S2607 + _S2615 + _S2615;
    float _S2619 = _S2609 + _S2616 + _S2616;
    float _S2620 = _S2611 + _S2617 + _S2617;
    float _S2621 = _S2605 + _S2614 + _S2614;
    float _S2622 = _S2576 + _S2580;
    float _S2623 = 0.5f * - _S2622;
    float _S2624 = _S2574 + _S2579.differential_0.y;
    DiffPair_float_0 _S2625;
    (&_S2625)->primal_0 = _S2175;
    (&_S2625)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2625, _S2624);
    float _S2626 = _S2623 + _S2625.differential_0;
    float _S2627 = _S2575 + _S2581 + _S2584.differential_0.x;
    DiffPair_float_0 _S2628;
    (&_S2628)->primal_0 = _S2173;
    (&_S2628)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2628, _S2627);
    float _S2629 = _S2623 + _S2628.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2630;
    (&_S2630)->primal_0 = R_16;
    (&_S2630)->differential_0 = _S2321;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2631;
    (&_S2631)->primal_0 = mean_14;
    (&_S2631)->differential_0 = _S2316;
    s_bwd_prop_mul_0(&_S2630, &_S2631, _S2451);
    float3  _S2632 = _S2504 + _S2548 + s_diff_vert0_c_T_1 + _S2451 + _S2393.differential_0;
    Matrix<float, 3, 3>  _S2633 = _S2565.differential_0 + _S2567.differential_0 + _S2569.differential_0 + _S2630.differential_0 + _S2452;
    float3  _S2634 = make_float3 (_S2629, _S2626, _S2622);
    float4  _S2635 = make_float4 (0.0f);
    *&((&_S2635)->w) = _S2618;
    *&((&_S2635)->z) = _S2619;
    *&((&_S2635)->y) = _S2620;
    *&((&_S2635)->x) = _S2621;
    float4  _S2636 = _S2635;
    float3  _S2637 = _S2571 + _S2577 + _S2582 + _S2631.differential_0 + _S2389;
    *v_mean_4 = _S2637;
    *v_quat_4 = _S2636;
    *v_scale_4 = _S2634;
    *v_hardness_1 = _S2448;
    (*v_sh_coeffs_4)[int(0)] = _S2455;
    (*v_sh_coeffs_4)[int(1)] = _S2456;
    (*v_sh_coeffs_4)[int(2)] = _S2457;
    (*v_sh_coeffs_4)[int(3)] = _S2458;
    (*v_sh_coeffs_4)[int(4)] = _S2459;
    (*v_sh_coeffs_4)[int(5)] = _S2460;
    (*v_sh_coeffs_4)[int(6)] = _S2461;
    (*v_sh_coeffs_4)[int(7)] = _S2462;
    (*v_sh_coeffs_4)[int(8)] = _S2463;
    (*v_sh_coeffs_4)[int(9)] = _S2464;
    (*v_sh_coeffs_4)[int(10)] = _S2465;
    (*v_sh_coeffs_4)[int(11)] = _S2466;
    (*v_sh_coeffs_4)[int(12)] = _S2467;
    (*v_sh_coeffs_4)[int(13)] = _S2468;
    (*v_sh_coeffs_4)[int(14)] = _S2469;
    (*v_sh_coeffs_4)[int(15)] = _S2470;
    *v_R_4 = _S2633;
    *v_t_4 = _S2632;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_18)
{
    DiffPair_float_0 _S2638 = *dpx_16;
    bool _S2639;
    if(((*dpx_16).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S2639 = ((*dpx_16).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S2639 = false;
    }
    float _S2640;
    if(_S2639)
    {
        _S2640 = dOut_18;
    }
    else
    {
        _S2640 = 0.0f;
    }
    dpx_16->primal_0 = _S2638.primal_0;
    dpx_16->differential_0 = _S2640;
    DiffPair_float_0 _S2641 = *dpMin_0;
    if((_S2638.primal_0) < ((*dpMin_0).primal_0))
    {
        _S2640 = dOut_18;
    }
    else
    {
        _S2640 = 0.0f;
    }
    dpMin_0->primal_0 = _S2641.primal_0;
    dpMin_0->differential_0 = _S2640;
    DiffPair_float_0 _S2642 = *dpMax_0;
    if(((*dpx_16).primal_0) > ((*dpMax_0).primal_0))
    {
        _S2640 = dOut_18;
    }
    else
    {
        _S2640 = 0.0f;
    }
    dpMax_0->primal_0 = _S2642.primal_0;
    dpMax_0->differential_0 = _S2640;
    return;
}

inline __device__ float clamp_0(float x_57, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_57), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_17, DiffPair_float_0 * dpy_6, float dOut_19)
{
    if(((*dpx_17).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = 0.0f;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_17).primal_0), ((*dpy_6).primal_0)));
        DiffPair_float_0 _S2643 = *dpx_17;
        float _S2644 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_19;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S2644;
        float _S2645 = val_0 * (F32_log((_S2643.primal_0))) * dOut_19;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S2645;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S2646 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S2646))));
    float2  _S2647 = p_0 - v0_0;
    float2  _S2648 = normalize_1(e0_6);
    float2  _S2649 = p_0 - v1_0;
    float2  _S2650 = normalize_1(e1_6);
    float2  _S2651 = p_0 - v2_0;
    float2  _S2652 = normalize_1(e2_2);
    float _S2653 = hardness_6.x;
    float _S2654 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.99500000476837158f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S2647.x * _S2648.y - _S2647.y * _S2648.x)), (se_0 * (_S2649.x * _S2650.y - _S2649.y * _S2650.x))))), (se_0 * (_S2651.x * _S2652.y - _S2651.y * _S2652.x)))) / ((F32_abs((_S2646))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S2654))));
    float _S2655;
    if(a_1 <= 0.0f)
    {
        _S2655 = 0.0f;
    }
    else
    {
        _S2655 = (F32_min(((F32_pow((a_1), (_S2654)))), (0.99900001287460327f)));
    }
    return _S2653 * _S2655;
}

inline __device__ float s_primal_ctx_abs_0(float _S2656)
{
    return (F32_abs((_S2656)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S2657, float _S2658, float _S2659)
{
    return clamp_0(_S2657, _S2658, _S2659);
}

inline __device__ float s_primal_ctx_exp2_0(float _S2660)
{
    return (F32_exp2((_S2660)));
}

inline __device__ float s_primal_ctx_pow_0(float _S2661, float _S2662)
{
    return (F32_pow((_S2661), (_S2662)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S2663, DiffPair_float_0 * _S2664, float _S2665)
{
    _d_pow_0(_S2663, _S2664, _S2665);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S2666, DiffPair_float_0 * _S2667, DiffPair_float_0 * _S2668, float _S2669)
{
    _d_clamp_0(_S2666, _S2667, _S2668, _S2669);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_7)
{
    float _S2670 = length_0((*dpx_18).primal_0);
    float2  _S2671 = (*dpx_18).primal_0 * _s_dOut_7;
    float2  _S2672 = make_float2 (1.0f / _S2670) * _s_dOut_7;
    float _S2673 = - ((_S2671.x + _S2671.y) / (_S2670 * _S2670));
    float2  _S2674 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2675;
    (&_S2675)->primal_0 = (*dpx_18).primal_0;
    (&_S2675)->differential_0 = _S2674;
    s_bwd_length_impl_0(&_S2675, _S2673);
    float2  _S2676 = _S2672 + _S2675.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S2676;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2677, float2  _S2678)
{
    s_bwd_prop_normalize_impl_1(_S2677, _S2678);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_8)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S2679 = e0_7.x;
    float _S2680 = e1_7.y;
    float _S2681 = e0_7.y;
    float _S2682 = e1_7.x;
    float _S2683 = _S2679 * _S2680 - _S2681 * _S2682;
    float se_1 = float((F32_sign((_S2683))));
    float2  _S2684 = p_1 - (*dpv0_0).primal_0;
    float2  _S2685 = normalize_1(e0_7);
    float _S2686 = _S2684.x;
    float _S2687 = _S2685.y;
    float _S2688 = _S2684.y;
    float _S2689 = _S2685.x;
    float de0_0 = se_1 * (_S2686 * _S2687 - _S2688 * _S2689);
    float2  _S2690 = p_1 - (*dpv1_0).primal_0;
    float2  _S2691 = normalize_1(e1_7);
    float _S2692 = _S2690.x;
    float _S2693 = _S2691.y;
    float _S2694 = _S2690.y;
    float _S2695 = _S2691.x;
    float de1_0 = se_1 * (_S2692 * _S2693 - _S2694 * _S2695);
    float2  _S2696 = p_1 - (*dpv2_0).primal_0;
    float2  _S2697 = normalize_1(e2_3);
    float _S2698 = _S2696.x;
    float _S2699 = _S2697.y;
    float _S2700 = _S2696.y;
    float _S2701 = _S2697.x;
    float de2_0 = se_1 * (_S2698 * _S2699 - _S2700 * _S2701);
    float _S2702 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S2703 = s_primal_ctx_max_0(_S2702, de2_0);
    float _S2704 = s_primal_ctx_abs_0(_S2683);
    float _S2705 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S2704 / _S2705;
    float _S2706 = _S2705 * _S2705;
    float _S2707 = (*dphardness_0).primal_0.x;
    float _S2708 = (*dphardness_0).primal_0.y;
    float _S2709 = dmax_0 * dmax_0;
    float _S2710 = 1.0f + _S2703 / dmax_0;
    float _S2711 = 1.0f - s_primal_ctx_clamp_0(_S2708, 0.00499999988824129f, 0.99500000476837158f);
    float _S2712 = -1.0f / _S2711;
    float _S2713 = _S2711 * _S2711;
    float _S2714 = 1.0f - s_primal_ctx_exp2_0(_S2712);
    float a_2 = 1.0f - _S2710 * _S2714;
    bool _S2715 = a_2 <= 0.0f;
    float _S2716;
    float _S2717;
    if(_S2715)
    {
        _S2716 = 0.0f;
        _S2717 = 0.0f;
    }
    else
    {
        float _S2718 = s_primal_ctx_pow_0(a_2, _S2711);
        _S2716 = s_primal_ctx_min_0(_S2718, 0.99900001287460327f);
        _S2717 = _S2718;
    }
    float _S2719 = _S2707 * _s_dOut_8;
    float _S2720 = _S2716 * _s_dOut_8;
    if(_S2715)
    {
        _S2716 = 0.0f;
        _S2717 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S2721;
        (&_S2721)->primal_0 = _S2717;
        (&_S2721)->differential_0 = 0.0f;
        DiffPair_float_0 _S2722;
        (&_S2722)->primal_0 = 0.99900001287460327f;
        (&_S2722)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S2721, &_S2722, _S2719);
        DiffPair_float_0 _S2723;
        (&_S2723)->primal_0 = a_2;
        (&_S2723)->differential_0 = 0.0f;
        DiffPair_float_0 _S2724;
        (&_S2724)->primal_0 = _S2711;
        (&_S2724)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S2723, &_S2724, _S2721.differential_0);
        _S2716 = _S2723.differential_0;
        _S2717 = _S2724.differential_0;
    }
    float _S2725 = - _S2716;
    float _S2726 = _S2714 * _S2725;
    float _S2727 = - (_S2710 * _S2725);
    DiffPair_float_0 _S2728;
    (&_S2728)->primal_0 = _S2712;
    (&_S2728)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2728, _S2727);
    float _S2729 = - (-1.0f * - (_S2728.differential_0 / _S2713) + _S2717);
    float _S2730 = _S2726 / _S2709;
    float s_diff_dmax_T_0 = _S2703 * - _S2730;
    float _S2731 = dmax_0 * _S2730;
    DiffPair_float_0 _S2732;
    (&_S2732)->primal_0 = _S2708;
    (&_S2732)->differential_0 = 0.0f;
    DiffPair_float_0 _S2733;
    (&_S2733)->primal_0 = 0.00499999988824129f;
    (&_S2733)->differential_0 = 0.0f;
    DiffPair_float_0 _S2734;
    (&_S2734)->primal_0 = 0.99500000476837158f;
    (&_S2734)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2732, &_S2733, &_S2734, _S2729);
    float _S2735 = s_diff_dmax_T_0 / _S2706;
    float _S2736 = _S2704 * - _S2735;
    float _S2737 = _S2705 * _S2735;
    float2  _S2738 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2739;
    (&_S2739)->primal_0 = e2_3;
    (&_S2739)->differential_0 = _S2738;
    s_bwd_length_impl_0(&_S2739, _S2736);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2740;
    (&_S2740)->primal_0 = e1_7;
    (&_S2740)->differential_0 = _S2738;
    s_bwd_length_impl_0(&_S2740, _S2736);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2741;
    (&_S2741)->primal_0 = e0_7;
    (&_S2741)->differential_0 = _S2738;
    s_bwd_length_impl_0(&_S2741, _S2736);
    DiffPair_float_0 _S2742;
    (&_S2742)->primal_0 = _S2683;
    (&_S2742)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2742, _S2737);
    DiffPair_float_0 _S2743;
    (&_S2743)->primal_0 = _S2702;
    (&_S2743)->differential_0 = 0.0f;
    DiffPair_float_0 _S2744;
    (&_S2744)->primal_0 = de2_0;
    (&_S2744)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2743, &_S2744, _S2731);
    DiffPair_float_0 _S2745;
    (&_S2745)->primal_0 = de0_0;
    (&_S2745)->differential_0 = 0.0f;
    DiffPair_float_0 _S2746;
    (&_S2746)->primal_0 = de1_0;
    (&_S2746)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2745, &_S2746, _S2743.differential_0);
    float _S2747 = se_1 * _S2744.differential_0;
    float _S2748 = - _S2747;
    float _S2749 = _S2701 * _S2748;
    float _S2750 = _S2699 * _S2747;
    float2  _S2751 = make_float2 (_S2700 * _S2748, _S2698 * _S2747);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2752;
    (&_S2752)->primal_0 = e2_3;
    (&_S2752)->differential_0 = _S2738;
    s_bwd_normalize_impl_1(&_S2752, _S2751);
    float2  _S2753 = - make_float2 (_S2750, _S2749);
    float _S2754 = se_1 * _S2746.differential_0;
    float _S2755 = - _S2754;
    float _S2756 = _S2695 * _S2755;
    float _S2757 = _S2693 * _S2754;
    float2  _S2758 = make_float2 (_S2694 * _S2755, _S2692 * _S2754);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2759;
    (&_S2759)->primal_0 = e1_7;
    (&_S2759)->differential_0 = _S2738;
    s_bwd_normalize_impl_1(&_S2759, _S2758);
    float2  _S2760 = - make_float2 (_S2757, _S2756);
    float _S2761 = se_1 * _S2745.differential_0;
    float _S2762 = - _S2761;
    float _S2763 = _S2689 * _S2762;
    float _S2764 = _S2687 * _S2761;
    float2  _S2765 = make_float2 (_S2688 * _S2762, _S2686 * _S2761);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2766;
    (&_S2766)->primal_0 = e0_7;
    (&_S2766)->differential_0 = _S2738;
    s_bwd_normalize_impl_1(&_S2766, _S2765);
    float2  _S2767 = - make_float2 (_S2764, _S2763);
    float _S2768 = - _S2742.differential_0;
    float2  _S2769 = _S2739.differential_0 + _S2752.differential_0;
    float2  _S2770 = - _S2769;
    float2  _S2771 = _S2740.differential_0 + _S2759.differential_0 + make_float2 (_S2681 * _S2768, _S2679 * _S2742.differential_0);
    float2  _S2772 = - _S2771;
    float2  _S2773 = _S2741.differential_0 + _S2766.differential_0 + make_float2 (_S2680 * _S2742.differential_0, _S2682 * _S2768);
    float2  _S2774 = - _S2773;
    float2  _S2775 = make_float2 (_S2720, _S2732.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S2775;
    float2  _S2776 = _S2753 + _S2770 + _S2771;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S2776;
    float2  _S2777 = _S2760 + _S2772 + _S2773;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S2777;
    float2  _S2778 = _S2767 + _S2769 + _S2774;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S2778;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2779, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2780, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2781, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2782, float2  _S2783, float _S2784)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S2779, _S2780, _S2781, _S2782, _S2783, _S2784);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_0, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S2785 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S2785;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S2785;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S2785;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S2785;
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
    float2  _S2786 = p_3 - v0_2;
    float2  _S2787 = p_3 - v1_2;
    float2  _S2788 = p_3 - v2_2;
    float _S2789 = e0_8.x;
    float _S2790 = e1_8.y;
    float _S2791 = e0_8.y;
    float _S2792 = e1_8.x;
    float _S2793 = _S2789 * _S2790 - _S2791 * _S2792;
    float se_2 = float((F32_sign((_S2793))));
    float _S2794 = hardness_8.x;
    float _S2795 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.99500000476837158f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S2786.x * _S2791 - _S2786.y * _S2789)), (se_2 * (_S2787.x * _S2790 - _S2787.y * _S2792))))), (se_2 * (_S2788.x * e2_4.y - _S2788.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S2786 - e0_8 * make_float2 (clamp_0(dot_1(_S2786, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S2787 - e1_8 * make_float2 (clamp_0(dot_1(_S2787, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S2788 - e2_4 * make_float2 (clamp_0(dot_1(_S2788, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S2793))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S2795))));
    float _S2796;
    if(a_3 <= 0.0f)
    {
        _S2796 = 0.0f;
    }
    else
    {
        _S2796 = (F32_min(((F32_pow((a_3), (_S2795)))), (0.99900001287460327f)));
    }
    return _S2794 * _S2796;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S2797, float2  _S2798)
{
    return dot_1(_S2797, _S2798);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2799, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2800, float _S2801)
{
    _d_dot_1(_S2799, _S2800, _S2801);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_9)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S2802 = p_4 - (*dpv0_1).primal_0;
    float _S2803 = s_primal_ctx_dot_1(_S2802, e0_9);
    float _S2804 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S2805 = _S2803 / _S2804;
    float _S2806 = _S2804 * _S2804;
    float _S2807 = s_primal_ctx_clamp_0(_S2805, 0.0f, 1.0f);
    float2  _S2808 = make_float2 (_S2807);
    float2  _S2809 = _S2802 - e0_9 * make_float2 (_S2807);
    float _S2810 = length_0(_S2809);
    float2  _S2811 = p_4 - (*dpv1_1).primal_0;
    float _S2812 = s_primal_ctx_dot_1(_S2811, e1_9);
    float _S2813 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S2814 = _S2812 / _S2813;
    float _S2815 = _S2813 * _S2813;
    float _S2816 = s_primal_ctx_clamp_0(_S2814, 0.0f, 1.0f);
    float2  _S2817 = make_float2 (_S2816);
    float2  _S2818 = _S2811 - e1_9 * make_float2 (_S2816);
    float _S2819 = length_0(_S2818);
    float2  _S2820 = p_4 - (*dpv2_1).primal_0;
    float _S2821 = s_primal_ctx_dot_1(_S2820, e2_5);
    float _S2822 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S2823 = _S2821 / _S2822;
    float _S2824 = _S2822 * _S2822;
    float _S2825 = s_primal_ctx_clamp_0(_S2823, 0.0f, 1.0f);
    float2  _S2826 = make_float2 (_S2825);
    float2  _S2827 = _S2820 - e2_5 * make_float2 (_S2825);
    float _S2828 = length_0(_S2827);
    float _S2829 = e0_9.x;
    float _S2830 = e1_9.y;
    float _S2831 = e0_9.y;
    float _S2832 = e1_9.x;
    float _S2833 = _S2829 * _S2830 - _S2831 * _S2832;
    float se_3 = float((F32_sign((_S2833))));
    float _S2834 = _S2802.x;
    float _S2835 = _S2802.y;
    float s0_0 = se_3 * (_S2834 * _S2831 - _S2835 * _S2829);
    float _S2836 = _S2811.x;
    float _S2837 = _S2811.y;
    float s1_0 = se_3 * (_S2836 * _S2830 - _S2837 * _S2832);
    float _S2838 = _S2820.x;
    float _S2839 = e2_5.y;
    float _S2840 = _S2820.y;
    float _S2841 = e2_5.x;
    float s2_0 = se_3 * (_S2838 * _S2839 - _S2840 * _S2841);
    float _S2842 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S2842, s2_0)))));
    float _S2843 = s_primal_ctx_min_0(_S2810, _S2819);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S2843, _S2828);
    float _S2844 = s_primal_ctx_abs_0(_S2833);
    float _S2845 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S2844 / _S2845;
    float _S2846 = _S2845 * _S2845;
    float _S2847 = (*dphardness_1).primal_0.x;
    float _S2848 = (*dphardness_1).primal_0.y;
    float _S2849 = dmax_1 * dmax_1;
    float _S2850 = 1.0f + dv_0 / dmax_1;
    float _S2851 = 1.0f - s_primal_ctx_clamp_0(_S2848, 0.00499999988824129f, 0.99500000476837158f);
    float _S2852 = -1.0f / _S2851;
    float _S2853 = _S2851 * _S2851;
    float _S2854 = 1.0f - s_primal_ctx_exp2_0(_S2852);
    float a_4 = 1.0f - _S2850 * _S2854;
    bool _S2855 = a_4 <= 0.0f;
    float _S2856;
    float _S2857;
    if(_S2855)
    {
        _S2856 = 0.0f;
        _S2857 = 0.0f;
    }
    else
    {
        float _S2858 = s_primal_ctx_pow_0(a_4, _S2851);
        _S2856 = s_primal_ctx_min_0(_S2858, 0.99900001287460327f);
        _S2857 = _S2858;
    }
    float _S2859 = _S2847 * _s_dOut_9;
    float _S2860 = _S2856 * _s_dOut_9;
    if(_S2855)
    {
        _S2856 = 0.0f;
        _S2857 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S2861;
        (&_S2861)->primal_0 = _S2857;
        (&_S2861)->differential_0 = 0.0f;
        DiffPair_float_0 _S2862;
        (&_S2862)->primal_0 = 0.99900001287460327f;
        (&_S2862)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S2861, &_S2862, _S2859);
        DiffPair_float_0 _S2863;
        (&_S2863)->primal_0 = a_4;
        (&_S2863)->differential_0 = 0.0f;
        DiffPair_float_0 _S2864;
        (&_S2864)->primal_0 = _S2851;
        (&_S2864)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S2863, &_S2864, _S2861.differential_0);
        _S2856 = _S2863.differential_0;
        _S2857 = _S2864.differential_0;
    }
    float _S2865 = - _S2856;
    float _S2866 = _S2854 * _S2865;
    float _S2867 = - (_S2850 * _S2865);
    DiffPair_float_0 _S2868;
    (&_S2868)->primal_0 = _S2852;
    (&_S2868)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2868, _S2867);
    float _S2869 = - (-1.0f * - (_S2868.differential_0 / _S2853) + _S2857);
    float _S2870 = _S2866 / _S2849;
    float s_diff_dmax_T_1 = dv_0 * - _S2870;
    float s_diff_dv_T_0 = dmax_1 * _S2870;
    DiffPair_float_0 _S2871;
    (&_S2871)->primal_0 = _S2848;
    (&_S2871)->differential_0 = 0.0f;
    DiffPair_float_0 _S2872;
    (&_S2872)->primal_0 = 0.00499999988824129f;
    (&_S2872)->differential_0 = 0.0f;
    DiffPair_float_0 _S2873;
    (&_S2873)->primal_0 = 0.99500000476837158f;
    (&_S2873)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2871, &_S2872, &_S2873, _S2869);
    float _S2874 = s_diff_dmax_T_1 / _S2846;
    float _S2875 = _S2844 * - _S2874;
    float _S2876 = _S2845 * _S2874;
    float2  _S2877 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2878;
    (&_S2878)->primal_0 = e2_5;
    (&_S2878)->differential_0 = _S2877;
    s_bwd_length_impl_0(&_S2878, _S2875);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2879;
    (&_S2879)->primal_0 = e1_9;
    (&_S2879)->differential_0 = _S2877;
    s_bwd_length_impl_0(&_S2879, _S2875);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2880;
    (&_S2880)->primal_0 = e0_9;
    (&_S2880)->differential_0 = _S2877;
    s_bwd_length_impl_0(&_S2880, _S2875);
    DiffPair_float_0 _S2881;
    (&_S2881)->primal_0 = _S2833;
    (&_S2881)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2881, _S2876);
    float _S2882 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S2883;
    (&_S2883)->primal_0 = _S2843;
    (&_S2883)->differential_0 = 0.0f;
    DiffPair_float_0 _S2884;
    (&_S2884)->primal_0 = _S2828;
    (&_S2884)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2883, &_S2884, _S2882);
    DiffPair_float_0 _S2885;
    (&_S2885)->primal_0 = _S2810;
    (&_S2885)->differential_0 = 0.0f;
    DiffPair_float_0 _S2886;
    (&_S2886)->primal_0 = _S2819;
    (&_S2886)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2885, &_S2886, _S2883.differential_0);
    DiffPair_float_0 _S2887;
    (&_S2887)->primal_0 = _S2842;
    (&_S2887)->differential_0 = 0.0f;
    DiffPair_float_0 _S2888;
    (&_S2888)->primal_0 = s2_0;
    (&_S2888)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2887, &_S2888, 0.0f);
    DiffPair_float_0 _S2889;
    (&_S2889)->primal_0 = s0_0;
    (&_S2889)->differential_0 = 0.0f;
    DiffPair_float_0 _S2890;
    (&_S2890)->primal_0 = s1_0;
    (&_S2890)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2889, &_S2890, _S2887.differential_0);
    float _S2891 = se_3 * _S2888.differential_0;
    float _S2892 = - _S2891;
    float _S2893 = _S2840 * _S2892;
    float _S2894 = _S2841 * _S2892;
    float _S2895 = _S2838 * _S2891;
    float _S2896 = _S2839 * _S2891;
    float _S2897 = se_3 * _S2890.differential_0;
    float _S2898 = - _S2897;
    float _S2899 = _S2832 * _S2898;
    float _S2900 = _S2830 * _S2897;
    float _S2901 = se_3 * _S2889.differential_0;
    float _S2902 = - _S2901;
    float _S2903 = _S2829 * _S2902;
    float _S2904 = _S2831 * _S2901;
    float _S2905 = - _S2881.differential_0;
    float _S2906 = _S2837 * _S2898 + _S2831 * _S2905;
    float _S2907 = _S2834 * _S2901 + _S2832 * _S2905;
    float _S2908 = _S2836 * _S2897 + _S2829 * _S2881.differential_0;
    float _S2909 = _S2835 * _S2902 + _S2830 * _S2881.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2910;
    (&_S2910)->primal_0 = _S2827;
    (&_S2910)->differential_0 = _S2877;
    s_bwd_length_impl_0(&_S2910, _S2884.differential_0);
    float2  _S2911 = - _S2910.differential_0;
    float2  _S2912 = e2_5 * _S2911;
    float2  _S2913 = _S2826 * _S2911;
    float _S2914 = _S2912.x + _S2912.y;
    DiffPair_float_0 _S2915;
    (&_S2915)->primal_0 = _S2823;
    (&_S2915)->differential_0 = 0.0f;
    DiffPair_float_0 _S2916;
    (&_S2916)->primal_0 = 0.0f;
    (&_S2916)->differential_0 = 0.0f;
    DiffPair_float_0 _S2917;
    (&_S2917)->primal_0 = 1.0f;
    (&_S2917)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2915, &_S2916, &_S2917, _S2914);
    float _S2918 = _S2915.differential_0 / _S2824;
    float _S2919 = _S2821 * - _S2918;
    float _S2920 = _S2822 * _S2918;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2921;
    (&_S2921)->primal_0 = e2_5;
    (&_S2921)->differential_0 = _S2877;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2922;
    (&_S2922)->primal_0 = e2_5;
    (&_S2922)->differential_0 = _S2877;
    s_bwd_prop_dot_1(&_S2921, &_S2922, _S2919);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2923;
    (&_S2923)->primal_0 = _S2820;
    (&_S2923)->differential_0 = _S2877;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2924;
    (&_S2924)->primal_0 = e2_5;
    (&_S2924)->differential_0 = _S2877;
    s_bwd_prop_dot_1(&_S2923, &_S2924, _S2920);
    float2  _S2925 = - (_S2910.differential_0 + _S2923.differential_0 + make_float2 (_S2896, _S2894));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2926;
    (&_S2926)->primal_0 = _S2818;
    (&_S2926)->differential_0 = _S2877;
    s_bwd_length_impl_0(&_S2926, _S2886.differential_0);
    float2  _S2927 = - _S2926.differential_0;
    float2  _S2928 = e1_9 * _S2927;
    float2  _S2929 = _S2817 * _S2927;
    float _S2930 = _S2928.x + _S2928.y;
    DiffPair_float_0 _S2931;
    (&_S2931)->primal_0 = _S2814;
    (&_S2931)->differential_0 = 0.0f;
    DiffPair_float_0 _S2932;
    (&_S2932)->primal_0 = 0.0f;
    (&_S2932)->differential_0 = 0.0f;
    DiffPair_float_0 _S2933;
    (&_S2933)->primal_0 = 1.0f;
    (&_S2933)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2931, &_S2932, &_S2933, _S2930);
    float _S2934 = _S2931.differential_0 / _S2815;
    float _S2935 = _S2812 * - _S2934;
    float _S2936 = _S2813 * _S2934;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2937;
    (&_S2937)->primal_0 = e1_9;
    (&_S2937)->differential_0 = _S2877;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2938;
    (&_S2938)->primal_0 = e1_9;
    (&_S2938)->differential_0 = _S2877;
    s_bwd_prop_dot_1(&_S2937, &_S2938, _S2935);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2939;
    (&_S2939)->primal_0 = _S2811;
    (&_S2939)->differential_0 = _S2877;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2940;
    (&_S2940)->primal_0 = e1_9;
    (&_S2940)->differential_0 = _S2877;
    s_bwd_prop_dot_1(&_S2939, &_S2940, _S2936);
    float2  _S2941 = - (_S2926.differential_0 + _S2939.differential_0 + make_float2 (_S2900, _S2899));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2942;
    (&_S2942)->primal_0 = _S2809;
    (&_S2942)->differential_0 = _S2877;
    s_bwd_length_impl_0(&_S2942, _S2885.differential_0);
    float2  _S2943 = - _S2942.differential_0;
    float2  _S2944 = e0_9 * _S2943;
    float2  _S2945 = _S2808 * _S2943;
    float _S2946 = _S2944.x + _S2944.y;
    DiffPair_float_0 _S2947;
    (&_S2947)->primal_0 = _S2805;
    (&_S2947)->differential_0 = 0.0f;
    DiffPair_float_0 _S2948;
    (&_S2948)->primal_0 = 0.0f;
    (&_S2948)->differential_0 = 0.0f;
    DiffPair_float_0 _S2949;
    (&_S2949)->primal_0 = 1.0f;
    (&_S2949)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2947, &_S2948, &_S2949, _S2946);
    float _S2950 = _S2947.differential_0 / _S2806;
    float _S2951 = _S2803 * - _S2950;
    float _S2952 = _S2804 * _S2950;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2953;
    (&_S2953)->primal_0 = e0_9;
    (&_S2953)->differential_0 = _S2877;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2954;
    (&_S2954)->primal_0 = e0_9;
    (&_S2954)->differential_0 = _S2877;
    s_bwd_prop_dot_1(&_S2953, &_S2954, _S2951);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2955;
    (&_S2955)->primal_0 = _S2802;
    (&_S2955)->differential_0 = _S2877;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2956;
    (&_S2956)->primal_0 = e0_9;
    (&_S2956)->differential_0 = _S2877;
    s_bwd_prop_dot_1(&_S2955, &_S2956, _S2952);
    float2  _S2957 = - (_S2942.differential_0 + _S2955.differential_0 + make_float2 (_S2904, _S2903));
    float2  _S2958 = _S2878.differential_0 + _S2913 + _S2922.differential_0 + _S2921.differential_0 + _S2924.differential_0 + make_float2 (_S2893, _S2895);
    float2  _S2959 = - _S2958;
    float2  _S2960 = _S2879.differential_0 + _S2929 + _S2938.differential_0 + _S2937.differential_0 + _S2940.differential_0 + make_float2 (_S2906, _S2908);
    float2  _S2961 = - _S2960;
    float2  _S2962 = _S2880.differential_0 + _S2945 + _S2954.differential_0 + _S2953.differential_0 + _S2956.differential_0 + make_float2 (_S2909, _S2907);
    float2  _S2963 = - _S2962;
    float2  _S2964 = make_float2 (_S2860, _S2871.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S2964;
    float2  _S2965 = _S2925 + _S2959 + _S2960;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S2965;
    float2  _S2966 = _S2941 + _S2961 + _S2962;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S2966;
    float2  _S2967 = _S2957 + _S2958 + _S2963;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S2967;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2968, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2969, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2970, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2971, float2  _S2972, float _S2973)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S2968, _S2969, _S2970, _S2971, _S2972, _S2973);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_1, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S2974 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S2974;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S2974;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S2974;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S2974;
    s_bwd_evaluate_alpha_opaque_triangle_precise_0(&dp_v0_1, &dp_v1_1, &dp_v2_1, &dp_hardness_1, p_5, v_alpha_1);
    *v_v0_1 = dp_v0_1.differential_0;
    *v_v1_1 = dp_v2_1.differential_0;
    *v_v2_1 = dp_v1_1.differential_0;
    *v_hardness_3 = dp_hardness_1.differential_0;
    return;
}

