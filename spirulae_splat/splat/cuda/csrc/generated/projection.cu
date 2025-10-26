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

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_8, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_8, float3  dOut_7)
{
    float _S10 = (*right_8).primal_0.rows[int(0)].x * dOut_7.x;
    Matrix<float, 3, 3>  right_d_result_4;
    *&(((&right_d_result_4)->rows + (int(0)))->x) = (*left_8).primal_0.x * dOut_7.x;
    float sum_14 = _S10 + (*right_8).primal_0.rows[int(0)].y * dOut_7.y;
    *&(((&right_d_result_4)->rows + (int(0)))->y) = (*left_8).primal_0.x * dOut_7.y;
    float sum_15 = sum_14 + (*right_8).primal_0.rows[int(0)].z * dOut_7.z;
    *&(((&right_d_result_4)->rows + (int(0)))->z) = (*left_8).primal_0.x * dOut_7.z;
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = sum_15;
    float _S11 = (*right_8).primal_0.rows[int(1)].x * dOut_7.x;
    *&(((&right_d_result_4)->rows + (int(1)))->x) = (*left_8).primal_0.y * dOut_7.x;
    float sum_16 = _S11 + (*right_8).primal_0.rows[int(1)].y * dOut_7.y;
    *&(((&right_d_result_4)->rows + (int(1)))->y) = (*left_8).primal_0.y * dOut_7.y;
    float sum_17 = sum_16 + (*right_8).primal_0.rows[int(1)].z * dOut_7.z;
    *&(((&right_d_result_4)->rows + (int(1)))->z) = (*left_8).primal_0.y * dOut_7.z;
    *&((&left_d_result_4)->y) = sum_17;
    float _S12 = (*right_8).primal_0.rows[int(2)].x * dOut_7.x;
    *&(((&right_d_result_4)->rows + (int(2)))->x) = (*left_8).primal_0.z * dOut_7.x;
    float sum_18 = _S12 + (*right_8).primal_0.rows[int(2)].y * dOut_7.y;
    *&(((&right_d_result_4)->rows + (int(2)))->y) = (*left_8).primal_0.z * dOut_7.y;
    float sum_19 = sum_18 + (*right_8).primal_0.rows[int(2)].z * dOut_7.z;
    *&(((&right_d_result_4)->rows + (int(2)))->z) = (*left_8).primal_0.z * dOut_7.z;
    *&((&left_d_result_4)->z) = sum_19;
    left_8->primal_0 = (*left_8).primal_0;
    left_8->differential_0 = left_d_result_4;
    right_8->primal_0 = (*right_8).primal_0;
    right_8->differential_0 = right_d_result_4;
    return;
}

inline __device__ float3  mul_7(float3  left_9, Matrix<float, 3, 3>  right_9)
{
    float3  result_9;
    int j_1 = int(0);
    for(;;)
    {
        if(j_1 < int(3))
        {
        }
        else
        {
            break;
        }
        int i_4 = int(0);
        float sum_20 = 0.0f;
        for(;;)
        {
            if(i_4 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_21 = sum_20 + _slang_vector_get_element(left_9, i_4) * _slang_vector_get_element(right_9.rows[i_4], j_1);
            i_4 = i_4 + int(1);
            sum_20 = sum_21;
        }
        *_slang_vector_get_element_ptr(&result_9, j_1) = sum_20;
        j_1 = j_1 + int(1);
    }
    return result_9;
}

inline __device__ float2  undistort_point_0(float2  uv_0, CameraDistortion_0 * dist_coeffs_0, int maxiter_0)
{
    int i_5 = int(0);
    float2  q_0 = uv_0;
    for(;;)
    {
        if(i_5 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float k4_0 = dist_coeffs_0->radial_coeffs_0.w;
        float p1_0 = dist_coeffs_0->tangential_coeffs_0.x;
        float p2_0 = dist_coeffs_0->tangential_coeffs_0.y;
        float sx1_0 = dist_coeffs_0->thin_prism_coeffs_0.x;
        float sy1_0 = dist_coeffs_0->thin_prism_coeffs_0.y;
        float u_0 = q_0.x;
        float v_0 = q_0.y;
        float r2_0 = u_0 * u_0 + v_0 * v_0;
        float _S13 = dist_coeffs_0->radial_coeffs_0.z + r2_0 * k4_0;
        float _S14 = dist_coeffs_0->radial_coeffs_0.y + r2_0 * _S13;
        float _S15 = dist_coeffs_0->radial_coeffs_0.x + r2_0 * _S14;
        float radial_0 = 1.0f + r2_0 * _S15;
        float _S16 = 2.0f * p1_0;
        float _S17 = _S16 * u_0;
        float _S18 = 2.0f * u_0;
        float _S19 = 2.0f * p2_0;
        float _S20 = _S19 * u_0;
        float _S21 = 2.0f * v_0;
        float2  _S22 = q_0 * make_float2 (radial_0) + make_float2 (_S17 * v_0 + p2_0 * (r2_0 + _S18 * u_0) + sx1_0 * r2_0, _S20 * v_0 + p1_0 * (r2_0 + _S21 * v_0) + sy1_0 * r2_0);
        float2  _S23 = make_float2 (0.0f);
        float2  seed_0 = _S23;
        *&((&seed_0)->x) = 1.0f;
        float2  _S24 = make_float2 (radial_0);
        float2  _S25 = q_0 * seed_0;
        float _S26 = p1_0 * seed_0.y;
        float _S27 = p2_0 * seed_0.x;
        float _S28 = _S25.x + _S25.y;
        float _S29 = r2_0 * _S28;
        float _S30 = r2_0 * _S29;
        float _S31 = sy1_0 * seed_0.y + _S26 + sx1_0 * seed_0.x + _S27 + _S15 * _S28 + _S14 * _S29 + _S13 * _S30 + k4_0 * (r2_0 * _S30);
        float _S32 = v_0 * _S31;
        float _S33 = u_0 * _S31;
        Matrix<float, 2, 2>  J_0;
        J_0[int(0)] = _S24 * seed_0 + make_float2 (_S19 * (v_0 * seed_0.y) + _S18 * _S27 + 2.0f * (u_0 * _S27) + _S16 * (v_0 * seed_0.x) + _S33 + _S33, _S21 * _S26 + 2.0f * (v_0 * _S26) + _S20 * seed_0.y + _S17 * seed_0.x + _S32 + _S32);
        float2  seed_1 = _S23;
        *&((&seed_1)->y) = 1.0f;
        float2  _S34 = q_0 * seed_1;
        float _S35 = p1_0 * seed_1.y;
        float _S36 = p2_0 * seed_1.x;
        float _S37 = _S34.x + _S34.y;
        float _S38 = r2_0 * _S37;
        float _S39 = r2_0 * _S38;
        float _S40 = sy1_0 * seed_1.y + _S35 + sx1_0 * seed_1.x + _S36 + _S15 * _S37 + _S14 * _S38 + _S13 * _S39 + k4_0 * (r2_0 * _S39);
        float _S41 = v_0 * _S40;
        float _S42 = u_0 * _S40;
        J_0[int(1)] = _S24 * seed_1 + make_float2 (_S19 * (v_0 * seed_1.y) + _S18 * _S36 + 2.0f * (u_0 * _S36) + _S16 * (v_0 * seed_1.x) + _S42 + _S42, _S21 * _S35 + 2.0f * (v_0 * _S35) + _S20 * seed_1.y + _S17 * seed_1.x + _S41 + _S41);
        float2  _S43 = _S22 - uv_0;
        float inv_det_0 = 1.0f / (J_0.rows[int(0)].x * J_0.rows[int(1)].y - J_0.rows[int(0)].y * J_0.rows[int(1)].x);
        float _S44 = _S43.x;
        float _S45 = _S43.y;
        float2  q_1 = q_0 - make_float2 ((_S44 * J_0.rows[int(1)].y - _S45 * J_0.rows[int(0)].y) * inv_det_0, (- _S44 * J_0.rows[int(1)].x + _S45 * J_0.rows[int(0)].x) * inv_det_0);
        i_5 = i_5 + int(1);
        q_0 = q_1;
    }
    return q_0;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_1, float dOut_8)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_1).primal_0.x * dOut_8;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_3).primal_0.x * dOut_8;
    *&((&x_d_result_0)->y) = (*dpy_1).primal_0.y * dOut_8;
    *&((&y_d_result_0)->y) = (*dpx_3).primal_0.y * dOut_8;
    *&((&x_d_result_0)->z) = (*dpy_1).primal_0.z * dOut_8;
    *&((&y_d_result_0)->z) = (*dpx_3).primal_0.z * dOut_8;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = x_d_result_0;
    dpy_1->primal_0 = (*dpy_1).primal_0;
    dpy_1->differential_0 = y_d_result_0;
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void _d_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_4, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpy_2, float dOut_9)
{
    float2  x_d_result_1;
    *&((&x_d_result_1)->x) = (*dpy_2).primal_0.x * dOut_9;
    float2  y_d_result_1;
    *&((&y_d_result_1)->x) = (*dpx_4).primal_0.x * dOut_9;
    *&((&x_d_result_1)->y) = (*dpy_2).primal_0.y * dOut_9;
    *&((&y_d_result_1)->y) = (*dpx_4).primal_0.y * dOut_9;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = x_d_result_1;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_1;
    return;
}

inline __device__ float dot_0(float3  x_9, float3  y_3)
{
    int i_6 = int(0);
    float result_10 = 0.0f;
    for(;;)
    {
        if(i_6 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_11 = result_10 + _slang_vector_get_element(x_9, i_6) * _slang_vector_get_element(y_3, i_6);
        i_6 = i_6 + int(1);
        result_10 = result_11;
    }
    return result_10;
}

inline __device__ float dot_1(float2  x_10, float2  y_4)
{
    int i_7 = int(0);
    float result_12 = 0.0f;
    for(;;)
    {
        if(i_7 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_13 = result_12 + _slang_vector_get_element(x_10, i_7) * _slang_vector_get_element(y_4, i_7);
        i_7 = i_7 + int(1);
        result_12 = result_13;
    }
    return result_12;
}

inline __device__ float length_0(float2  x_11)
{
    return (F32_sqrt((dot_1(x_11, x_11))));
}

inline __device__ float length_1(float3  x_12)
{
    return (F32_sqrt((dot_0(x_12, x_12))));
}

inline __device__ float3  normalize_0(float3  x_13)
{
    return x_13 / make_float3 (length_1(x_13));
}

inline __device__ float2  normalize_1(float2  x_14)
{
    return x_14 / make_float2 (length_0(x_14));
}

inline __device__ void generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_1, bool is_fisheye_0, float4  radial_coeffs_1, float2  tangential_coeffs_1, float2  thin_prism_coeffs_1, float3  * ray_o_0, float3  * ray_d_0)
{
    *ray_o_0 = - mul_7(t_1, R_2);
    CameraDistortion_0 _S46;
    (&_S46)->radial_coeffs_0 = radial_coeffs_1;
    (&_S46)->tangential_coeffs_0 = tangential_coeffs_1;
    (&_S46)->thin_prism_coeffs_0 = thin_prism_coeffs_1;
    float2  _S47 = undistort_point_0(uv_1, &_S46, int(8));
    float3  raydir_0;
    if(is_fisheye_0)
    {
        float theta_0 = length_0(_S47);
        float _S48;
        if(theta_0 < 0.00100000004749745f)
        {
            _S48 = 1.0f - theta_0 * theta_0 / 6.0f;
        }
        else
        {
            _S48 = (F32_sin((theta_0))) / theta_0;
        }
        float3  _S49 = make_float3 ((_S47 * make_float2 (_S48)).x, (_S47 * make_float2 (_S48)).y, (F32_cos((theta_0))));
        raydir_0 = _S49;
    }
    else
    {
        raydir_0 = make_float3 (_S47.x, _S47.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_7(raydir_0, R_2));
    return;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S50;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S51, Matrix<float, 3, 3>  _S52)
{
    return mul_7(_S51, _S52);
}

inline __device__ float s_primal_ctx_sin_0(float _S53)
{
    return (F32_sin((_S53)));
}

inline __device__ float s_primal_ctx_cos_0(float _S54)
{
    return (F32_cos((_S54)));
}

inline __device__ void s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_2, bool is_fisheye_1, float4  radial_coeffs_2, float2  tangential_coeffs_2, float2  thin_prism_coeffs_2, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S50 = make_float2 (0.0f);
    float3  _S55 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    CameraDistortion_0 _S56;
    (&_S56)->radial_coeffs_0 = radial_coeffs_2;
    (&_S56)->tangential_coeffs_0 = tangential_coeffs_2;
    (&_S56)->thin_prism_coeffs_0 = thin_prism_coeffs_2;
    float2  _S57 = undistort_point_0(uv_2, &_S56, int(8));
    _s_diff_ctx_0->_S50 = _S57;
    float3  raydir_1;
    if(is_fisheye_1)
    {
        float _S58 = length_0(_S57);
        float _S59;
        if(_S58 < 0.00100000004749745f)
        {
            _S59 = 1.0f - _S58 * _S58 / 6.0f;
        }
        else
        {
            _S59 = s_primal_ctx_sin_0(_S58) / _S58;
        }
        float3  _S60 = make_float3 ((_S57 * make_float2 (_S59)).x, (_S57 * make_float2 (_S59)).y, s_primal_ctx_cos_0(_S58));
        raydir_1 = _S60;
    }
    else
    {
        raydir_1 = make_float3 (_S57.x, _S57.y, 1.0f);
    }
    float3  _S61 = normalize_0(s_primal_ctx_mul_0(raydir_1, dpR_0));
    *dpray_o_0 = _S55;
    *dpray_d_0 = _S61;
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S62, float _S63)
{
    _d_sqrt_0(_S62, _S63);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float _s_dOut_0)
{
    float _S64 = (*dpx_5).primal_0.x;
    float _S65 = (*dpx_5).primal_0.y;
    float _S66 = (*dpx_5).primal_0.z;
    DiffPair_float_0 _S67;
    (&_S67)->primal_0 = _S64 * _S64 + _S65 * _S65 + _S66 * _S66;
    (&_S67)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S67, _s_dOut_0);
    float _S68 = (*dpx_5).primal_0.z * _S67.differential_0;
    float _S69 = _S68 + _S68;
    float _S70 = (*dpx_5).primal_0.y * _S67.differential_0;
    float _S71 = _S70 + _S70;
    float _S72 = (*dpx_5).primal_0.x * _S67.differential_0;
    float _S73 = _S72 + _S72;
    float3  _S74 = make_float3 (0.0f);
    *&((&_S74)->z) = _S69;
    *&((&_S74)->y) = _S71;
    *&((&_S74)->x) = _S73;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S74;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S75, float _S76)
{
    s_bwd_prop_length_impl_0(_S75, _S76);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float3  _s_dOut_1)
{
    float _S77 = length_1((*dpx_6).primal_0);
    float3  _S78 = (*dpx_6).primal_0 * _s_dOut_1;
    float3  _S79 = make_float3 (1.0f / _S77) * _s_dOut_1;
    float _S80 = - ((_S78.x + _S78.y + _S78.z) / (_S77 * _S77));
    float3  _S81 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S82;
    (&_S82)->primal_0 = (*dpx_6).primal_0;
    (&_S82)->differential_0 = _S81;
    s_bwd_length_impl_0(&_S82, _S80);
    float3  _S83 = _S79 + _S82.differential_0;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S83;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S84, float3  _S85)
{
    s_bwd_prop_normalize_impl_0(_S84, _S85);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S86, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S87, float3  _S88)
{
    _d_mul_1(_S86, _S87, _S88);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_3, bool is_fisheye_2, float4  radial_coeffs_3, float2  tangential_coeffs_3, float2  thin_prism_coeffs_3, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S89 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S90 = *dpt_1;
    float3  raydir_2;
    if(is_fisheye_2)
    {
        float _S91 = length_0(_s_diff_ctx_1->_S50);
        float _S92;
        if(_S91 < 0.00100000004749745f)
        {
            _S92 = 1.0f - _S91 * _S91 / 6.0f;
        }
        else
        {
            _S92 = s_primal_ctx_sin_0(_S91) / _S91;
        }
        float3  _S93 = make_float3 ((_s_diff_ctx_1->_S50 * make_float2 (_S92)).x, (_s_diff_ctx_1->_S50 * make_float2 (_S92)).y, s_primal_ctx_cos_0(_S91));
        raydir_2 = _S93;
    }
    else
    {
        raydir_2 = make_float3 (_s_diff_ctx_1->_S50.x, _s_diff_ctx_1->_S50.y, 1.0f);
    }
    float3  _S94 = s_primal_ctx_mul_0(raydir_2, _S89.primal_0);
    float3  _S95 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S96;
    (&_S96)->primal_0 = _S94;
    (&_S96)->differential_0 = _S95;
    s_bwd_normalize_impl_0(&_S96, dpray_d_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S97;
    (&_S97)->primal_0 = raydir_2;
    (&_S97)->differential_0 = _S95;
    Matrix<float, 3, 3>  _S98 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S99;
    (&_S99)->primal_0 = _S89.primal_0;
    (&_S99)->differential_0 = _S98;
    s_bwd_prop_mul_0(&_S97, &_S99, _S96.differential_0);
    float3  _S100 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S101;
    (&_S101)->primal_0 = _S90.primal_0;
    (&_S101)->differential_0 = _S95;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S102;
    (&_S102)->primal_0 = _S89.primal_0;
    (&_S102)->differential_0 = _S98;
    s_bwd_prop_mul_0(&_S101, &_S102, _S100);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S101.differential_0;
    Matrix<float, 3, 3>  _S103 = _S102.differential_0 + _S99.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S103;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S104, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S105, float2  _S106, bool _S107, float4  _S108, float2  _S109, float2  _S110, float3  _S111, float3  _S112)
{
    float3  _S113;
    float3  _S114;
    s_bwd_prop_generate_ray_Intermediates_0 _S115;
    s_primal_ctx_generate_ray_0((*_S104).primal_0, (*_S105).primal_0, _S106, _S107, _S108, _S109, _S110, &_S113, &_S114, &_S115);
    s_bwd_prop_generate_ray_Intermediates_0 _S116 = _S115;
    s_bwd_prop_generate_ray_0(_S104, _S105, _S106, _S107, _S108, _S109, _S110, _S111, _S112, &_S116);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_4, bool is_fisheye_3, float4  radial_coeffs_4, float2  tangential_coeffs_4, float2  thin_prism_coeffs_4, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S117 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S117;
    float3  _S118 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S118;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_4, is_fisheye_3, radial_coeffs_4, tangential_coeffs_4, thin_prism_coeffs_4, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_7, DiffPair_float_0 * dpy_3, float dOut_10)
{
    DiffPair_float_0 _S119 = *dpx_7;
    float _S120;
    if(((*dpx_7).primal_0) < ((*dpy_3).primal_0))
    {
        _S120 = dOut_10;
    }
    else
    {
        if(((*dpx_7).primal_0) > ((*dpy_3).primal_0))
        {
            _S120 = 0.0f;
        }
        else
        {
            _S120 = 0.5f * dOut_10;
        }
    }
    dpx_7->primal_0 = _S119.primal_0;
    dpx_7->differential_0 = _S120;
    DiffPair_float_0 _S121 = *dpy_3;
    if(((*dpy_3).primal_0) < (_S119.primal_0))
    {
        _S120 = dOut_10;
    }
    else
    {
        if(((*dpy_3).primal_0) > ((*dpx_7).primal_0))
        {
            _S120 = 0.0f;
        }
        else
        {
            _S120 = 0.5f * dOut_10;
        }
    }
    dpy_3->primal_0 = _S121.primal_0;
    dpy_3->differential_0 = _S120;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S122 = float(width_0);
    float _S123 = float(height_0);
    float _S124 = 0.30000001192092896f * (0.5f * _S122 / fx_0);
    float _S125 = 0.30000001192092896f * (0.5f * _S123 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_1 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S122 - cx_0) / fx_0 + _S124), ((F32_max((- (cx_0 / fx_0 + _S124)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S123 - cy_0) / fy_0 + _S125), ((F32_max((- (cy_0 / fy_0 + _S125)), (mean3d_0.y * rz_0))))))) * rz2_0);
    *cov2d_0 = mul_6(mul_5(J_1, cov3d_0), transpose_1(J_1));
    *mean2d_0 = make_float2 (fx_0 * mean3d_0.x * rz_0 + cx_0, fy_0 * mean3d_0.y * rz_0 + cy_0);
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float fx_1, float fy_1, float cx_1, float cy_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float rz_1 = 1.0f / mean3d_1.z;
    float rz2_1 = rz_1 * rz_1;
    Matrix<float, 2, 3>  J_2 = makeMatrix<float, 2, 3> (fx_1 * rz_1, 0.0f, - fx_1 * mean3d_1.x * rz2_1, 0.0f, fy_1 * rz_1, - fy_1 * mean3d_1.y * rz2_1);
    *cov2d_1 = mul_6(mul_5(J_2, cov3d_1), transpose_1(J_2));
    *mean2d_1 = make_float2 (fx_1 * mean3d_1.x * rz_1 + cx_1, fy_1 * mean3d_1.y * rz_1 + cy_1);
    return;
}

inline __device__ CameraDistortion_0 CameraDistortion_x24init_0(float4  radial_coeffs_5, float2  tangential_coeffs_5, float2  thin_prism_coeffs_5)
{
    CameraDistortion_0 _S126;
    (&_S126)->radial_coeffs_0 = radial_coeffs_5;
    (&_S126)->tangential_coeffs_0 = tangential_coeffs_5;
    (&_S126)->thin_prism_coeffs_0 = thin_prism_coeffs_5;
    return _S126;
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_8, float dOut_11)
{
    DiffPair_float_0 _S127 = *dpx_8;
    float _S128 = - (*dpy_4).primal_0 / ((*dpx_8).primal_0 * (*dpx_8).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_11;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S128;
    float _S129 = _S127.primal_0 / (_S127.primal_0 * _S127.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_11;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S129;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S130, float _S131)
{
    return (F32_atan2((_S130), (_S131)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S132, DiffPair_float_0 * _S133, float _S134)
{
    _d_atan2_0(_S132, _S133, _S134);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_9, float _s_dOut_2)
{
    float _S135 = (*dpx_9).primal_0.x;
    float _S136 = (*dpx_9).primal_0.y;
    DiffPair_float_0 _S137;
    (&_S137)->primal_0 = _S135 * _S135 + _S136 * _S136;
    (&_S137)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S137, _s_dOut_2);
    float _S138 = (*dpx_9).primal_0.y * _S137.differential_0;
    float _S139 = _S138 + _S138;
    float _S140 = (*dpx_9).primal_0.x * _S137.differential_0;
    float _S141 = _S140 + _S140;
    float2  _S142 = make_float2 (0.0f);
    *&((&_S142)->y) = _S139;
    *&((&_S142)->x) = _S141;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S142;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S143, float _S144)
{
    s_bwd_prop_length_impl_1(_S143, _S144);
    return;
}

inline __device__ void fisheye_proj_3dgs(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float fx_2, float fy_2, float cx_2, float cy_2, float4  radial_coeffs_6, float2  tangential_coeffs_6, float2  thin_prism_coeffs_6, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    CameraDistortion_0 dist_coeffs_1 = CameraDistortion_x24init_0(radial_coeffs_6, tangential_coeffs_6, thin_prism_coeffs_6);
    float2  _S145 = float2 {mean3d_2.x, mean3d_2.y};
    float r_6 = length_0(_S145);
    float _S146 = mean3d_2.z;
    float theta_1 = (F32_atan2((r_6), (_S146)));
    float k_0;
    if(theta_1 < 0.00100000004749745f)
    {
        k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S146;
    }
    else
    {
        k_0 = theta_1 / r_6;
    }
    float2  _S147 = _S145 * make_float2 (k_0);
    float k1_0 = dist_coeffs_1.radial_coeffs_0.x;
    float k2_0 = dist_coeffs_1.radial_coeffs_0.y;
    float k3_0 = dist_coeffs_1.radial_coeffs_0.z;
    float k4_1 = dist_coeffs_1.radial_coeffs_0.w;
    float p1_1 = dist_coeffs_1.tangential_coeffs_0.x;
    float p2_1 = dist_coeffs_1.tangential_coeffs_0.y;
    float sx1_1 = dist_coeffs_1.thin_prism_coeffs_0.x;
    float sy1_1 = dist_coeffs_1.thin_prism_coeffs_0.y;
    float u_1 = _S147.x;
    float v_1 = _S147.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float _S148 = 2.0f * p1_1;
    float _S149 = 2.0f * p2_1;
    float2  _S150 = _S147 * make_float2 (1.0f + r2_1 * (k1_0 + r2_1 * (k2_0 + r2_1 * (k3_0 + r2_1 * k4_1)))) + make_float2 (_S148 * u_1 * v_1 + p2_1 * (r2_1 + 2.0f * u_1 * u_1) + sx1_1 * r2_1, _S149 * u_1 * v_1 + p1_1 * (r2_1 + 2.0f * v_1 * v_1) + sy1_1 * r2_1);
    *mean2d_2 = make_float2 (fx_2 * _S150.x + cx_2, fy_2 * _S150.y + cy_2);
    Matrix<float, 2, 3>  J_3;
    float2  _S151 = make_float2 (0.0f);
    float2  seed_2 = _S151;
    *&((&seed_2)->x) = 1.0f;
    float2  _S152 = seed_2;
    float _S153 = s_primal_ctx_atan2_0(r_6, _S146);
    bool _S154 = _S153 < 0.00100000004749745f;
    float _S155;
    float _S156;
    float _S157;
    if(_S154)
    {
        float _S158 = 1.0f - _S153 * _S153 / 3.0f;
        float _S159 = _S146 * _S146;
        k_0 = _S158 / _S146;
        _S155 = 0.0f;
        _S156 = _S159;
        _S157 = _S158;
    }
    else
    {
        float _S160 = r_6 * r_6;
        k_0 = _S153 / r_6;
        _S155 = _S160;
        _S156 = 0.0f;
        _S157 = 0.0f;
    }
    float2  _S161 = make_float2 (k_0);
    float2  _S162 = _S145 * make_float2 (k_0);
    float u_2 = _S162.x;
    float v_2 = _S162.y;
    float r2_2 = u_2 * u_2 + v_2 * v_2;
    float _S163 = k3_0 + r2_2 * k4_1;
    float _S164 = k2_0 + r2_2 * _S163;
    float _S165 = k1_0 + r2_2 * _S164;
    float _S166 = fy_2 * _S152.y;
    float _S167 = fx_2 * _S152.x;
    float2  _S168 = make_float2 (_S167, _S166);
    float2  _S169 = _S162 * _S168;
    float _S170 = p1_1 * _S166;
    float _S171 = p2_1 * _S167;
    float _S172 = _S169.x + _S169.y;
    float _S173 = r2_2 * _S172;
    float _S174 = r2_2 * _S173;
    float _S175 = sy1_1 * _S166 + _S170 + sx1_1 * _S167 + _S171 + _S165 * _S172 + _S164 * _S173 + _S163 * _S174 + k4_1 * (r2_2 * _S174);
    float _S176 = v_2 * _S175;
    float _S177 = u_2 * _S175;
    float2  _S178 = make_float2 (1.0f + r2_2 * _S165) * _S168 + make_float2 (_S149 * (v_2 * _S166) + 2.0f * u_2 * _S171 + 2.0f * (u_2 * _S171) + _S148 * (v_2 * _S167) + _S177 + _S177, 2.0f * v_2 * _S170 + 2.0f * (v_2 * _S170) + _S149 * u_2 * _S166 + _S148 * u_2 * _S167 + _S176 + _S176);
    float2  _S179 = _S145 * _S178;
    float2  _S180 = _S161 * _S178;
    float _S181 = _S179.x + _S179.y;
    if(_S154)
    {
        float _S182 = _S181 / _S156;
        float _S183 = _S157 * - _S182;
        float _S184 = _S153 * (0.3333333432674408f * - (_S146 * _S182));
        k_0 = _S184 + _S184;
        _S155 = _S183;
        _S156 = 0.0f;
    }
    else
    {
        float _S185 = _S181 / _S155;
        float _S186 = _S153 * - _S185;
        k_0 = r_6 * _S185;
        _S155 = 0.0f;
        _S156 = _S186;
    }
    DiffPair_float_0 _S187;
    (&_S187)->primal_0 = r_6;
    (&_S187)->differential_0 = 0.0f;
    DiffPair_float_0 _S188;
    (&_S188)->primal_0 = _S146;
    (&_S188)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S187, &_S188, k_0);
    float _S189 = _S188.differential_0 + _S155;
    float _S190 = _S187.differential_0 + _S156;
    float2  _S191 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S192;
    (&_S192)->primal_0 = _S145;
    (&_S192)->differential_0 = _S191;
    s_bwd_length_impl_1(&_S192, _S190);
    float2  _S193 = _S192.differential_0 + _S180;
    float3  _S194 = make_float3 (_S193.x, _S193.y, _S189);
    J_3[int(0)] = _S194;
    float2  seed_3 = _S151;
    *&((&seed_3)->y) = 1.0f;
    float2  _S195 = seed_3;
    if(_S154)
    {
        float _S196 = 1.0f - _S153 * _S153 / 3.0f;
        float _S197 = _S146 * _S146;
        k_0 = _S196 / _S146;
        _S155 = 0.0f;
        _S156 = _S197;
        _S157 = _S196;
    }
    else
    {
        float _S198 = r_6 * r_6;
        k_0 = _S153 / r_6;
        _S155 = _S198;
        _S156 = 0.0f;
        _S157 = 0.0f;
    }
    float2  _S199 = make_float2 (k_0);
    float2  _S200 = _S145 * make_float2 (k_0);
    float u_3 = _S200.x;
    float v_3 = _S200.y;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float _S201 = k3_0 + r2_3 * k4_1;
    float _S202 = k2_0 + r2_3 * _S201;
    float _S203 = k1_0 + r2_3 * _S202;
    float _S204 = fy_2 * _S195.y;
    float _S205 = fx_2 * _S195.x;
    float2  _S206 = make_float2 (_S205, _S204);
    float2  _S207 = _S200 * _S206;
    float _S208 = p1_1 * _S204;
    float _S209 = p2_1 * _S205;
    float _S210 = _S207.x + _S207.y;
    float _S211 = r2_3 * _S210;
    float _S212 = r2_3 * _S211;
    float _S213 = sy1_1 * _S204 + _S208 + sx1_1 * _S205 + _S209 + _S203 * _S210 + _S202 * _S211 + _S201 * _S212 + k4_1 * (r2_3 * _S212);
    float _S214 = v_3 * _S213;
    float _S215 = u_3 * _S213;
    float2  _S216 = make_float2 (1.0f + r2_3 * _S203) * _S206 + make_float2 (_S149 * (v_3 * _S204) + 2.0f * u_3 * _S209 + 2.0f * (u_3 * _S209) + _S148 * (v_3 * _S205) + _S215 + _S215, 2.0f * v_3 * _S208 + 2.0f * (v_3 * _S208) + _S149 * u_3 * _S204 + _S148 * u_3 * _S205 + _S214 + _S214);
    float2  _S217 = _S145 * _S216;
    float2  _S218 = _S199 * _S216;
    float _S219 = _S217.x + _S217.y;
    if(_S154)
    {
        float _S220 = _S219 / _S156;
        float _S221 = _S157 * - _S220;
        float _S222 = _S153 * (0.3333333432674408f * - (_S146 * _S220));
        k_0 = _S222 + _S222;
        _S155 = _S221;
        _S156 = 0.0f;
    }
    else
    {
        float _S223 = _S219 / _S155;
        float _S224 = _S153 * - _S223;
        k_0 = r_6 * _S223;
        _S155 = 0.0f;
        _S156 = _S224;
    }
    DiffPair_float_0 _S225;
    (&_S225)->primal_0 = r_6;
    (&_S225)->differential_0 = 0.0f;
    DiffPair_float_0 _S226;
    (&_S226)->primal_0 = _S146;
    (&_S226)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S225, &_S226, k_0);
    float _S227 = _S226.differential_0 + _S155;
    float _S228 = _S225.differential_0 + _S156;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S229;
    (&_S229)->primal_0 = _S145;
    (&_S229)->differential_0 = _S191;
    s_bwd_length_impl_1(&_S229, _S228);
    float2  _S230 = _S229.differential_0 + _S218;
    float3  _S231 = make_float3 (_S230.x, _S230.y, _S227);
    J_3[int(1)] = _S231;
    *cov2d_2 = mul_6(mul_5(J_3, cov3d_2), transpose_1(J_3));
    return;
}

inline __device__ void ortho_proj_3dgs(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float fx_3, float fy_3, float cx_3, float cy_3, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    Matrix<float, 2, 3>  J_4 = makeMatrix<float, 2, 3> (fx_3, 0.0f, 0.0f, 0.0f, fy_3, 0.0f);
    *cov2d_3 = mul_6(mul_5(J_4, cov3d_3), transpose_1(J_4));
    *mean2d_3 = make_float2 (fx_3 * mean3d_3.x + cx_3, fy_3 * mean3d_3.y + cy_3);
    return;
}

inline __device__ float add_blur(float eps2d_0, Matrix<float, 2, 2>  * covar_1, float * compensation_0)
{
    float det_orig_0 = *&((covar_1->rows + (int(0)))->x) * *&((covar_1->rows + (int(1)))->y) - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    float _S232 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S232;
    float _S233 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S233;
    float det_blur_0 = _S232 * _S233 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_10, float dOut_12)
{
    float _S234 = (F32_exp(((*dpx_10).primal_0))) * dOut_12;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S234;
    return;
}

inline __device__ float3  exp_0(float3  x_15)
{
    float3  result_14;
    int i_8 = int(0);
    for(;;)
    {
        if(i_8 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_14, i_8) = (F32_exp((_slang_vector_get_element(x_15, i_8))));
        i_8 = i_8 + int(1);
    }
    return result_14;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, float3  dOut_13)
{
    float3  _S235 = exp_0((*dpx_11).primal_0) * dOut_13;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S235;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_12, float dOut_14)
{
    float _S236 = 1.0f / (*dpx_12).primal_0 * dOut_14;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S236;
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_13, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_5, float3  dOut_15)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_13).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpy_5).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_max_0(&left_dp_0, &right_dp_0, dOut_15.x);
    float3  left_d_result_5;
    *&((&left_d_result_5)->x) = left_dp_0.differential_0;
    float3  right_d_result_5;
    *&((&right_d_result_5)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_13).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_5).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_15.y);
    *&((&left_d_result_5)->y) = left_dp_1.differential_0;
    *&((&right_d_result_5)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_13).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_5).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_15.z);
    *&((&left_d_result_5)->z) = left_dp_2.differential_0;
    *&((&right_d_result_5)->z) = right_dp_2.differential_0;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = left_d_result_5;
    dpy_5->primal_0 = (*dpy_5).primal_0;
    dpy_5->differential_0 = right_d_result_5;
    return;
}

inline __device__ float3  max_0(float3  x_16, float3  y_5)
{
    float3  result_15;
    int i_9 = int(0);
    for(;;)
    {
        if(i_9 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_15, i_9) = (F32_max((_slang_vector_get_element(x_16, i_9)), (_slang_vector_get_element(y_5, i_9))));
        i_9 = i_9 + int(1);
    }
    return result_15;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_4, float fy_4, float cx_4, float cy_4, float4  radial_coeffs_7, float2  tangential_coeffs_7, float2  thin_prism_coeffs_7, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_4, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_4, mean_0) + t_3;
        float _S237 = mean_c_0.z;
        bool _S238;
        if(_S237 < near_plane_0)
        {
            _S238 = true;
        }
        else
        {
            _S238 = _S237 > far_plane_0;
        }
        if(_S238)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S239 = exp_0(scale_2);
        float x_17 = quat_3.y;
        float inv_norm_3 = (F32_rsqrt((x_17 * x_17 + quat_3.z * quat_3.z + quat_3.w * quat_3.w + quat_3.x * quat_3.x)));
        float x_18 = quat_3.y * inv_norm_3;
        float y_6 = quat_3.z * inv_norm_3;
        float z_3 = quat_3.w * inv_norm_3;
        float w_3 = quat_3.x * inv_norm_3;
        float x2_3 = x_18 * x_18;
        float y2_3 = y_6 * y_6;
        float z2_3 = z_3 * z_3;
        float xy_3 = x_18 * y_6;
        float xz_3 = x_18 * z_3;
        float yz_3 = y_6 * z_3;
        float wx_3 = w_3 * x_18;
        float wy_3 = w_3 * y_6;
        float wz_3 = w_3 * z_3;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S239.x, 0.0f, 0.0f, 0.0f, _S239.y, 0.0f, 0.0f, 0.0f, _S239.z));
        Matrix<float, 3, 3>  _S240 = transpose_0(R_4);
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_4, mul_4(M_2, transpose_0(M_2))), _S240);
        Matrix<float, 2, 2>  covar2d_0;
        float _S241 = float(image_width_0);
        float _S242 = float(image_height_0);
        float _S243 = 0.30000001192092896f * (0.5f * _S241 / fx_4);
        float _S244 = 0.30000001192092896f * (0.5f * _S242 / fy_4);
        float rz_2 = 1.0f / mean_c_0.z;
        float rz2_2 = rz_2 * rz_2;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_4 * rz_2, 0.0f, - fx_4 * (mean_c_0.z * (F32_min(((_S241 - cx_4) / fx_4 + _S243), ((F32_max((- (cx_4 / fx_4 + _S243)), (mean_c_0.x * rz_2))))))) * rz2_2, 0.0f, fy_4 * rz_2, - fy_4 * (mean_c_0.z * (F32_min(((_S242 - cy_4) / fy_4 + _S244), ((F32_max((- (cy_4 / fy_4 + _S244)), (mean_c_0.y * rz_2))))))) * rz2_2);
        covar2d_0 = mul_6(mul_5(J_5, covar_c_0), transpose_1(J_5));
        *mean2d_4 = make_float2 (fx_4 * mean_c_0.x * rz_2 + cx_4, fy_4 * mean_c_0.y * rz_2 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S245 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S245;
        float _S246 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S246;
        float det_blur_1 = _S245 * _S246 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S247 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S238 = true;
        }
        else
        {
            _S238 = xmin_0 >= _S241;
        }
        if(_S238)
        {
            _S238 = true;
        }
        else
        {
            _S238 = ymax_0 <= 0.0f;
        }
        if(_S238)
        {
            _S238 = true;
        }
        else
        {
            _S238 = ymin_0 >= _S242;
        }
        if(_S238)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = length_1(mean_c_0);
        *conic_0 = make_float3 (_S247.rows[int(0)].x, _S247.rows[int(0)].y, _S247.rows[int(1)].y);
        float3  _S248 = mean_0 - - mul_0(_S240, t_3);
        float3  _S249 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S249;
        float _S250 = _S248.x;
        float _S251 = _S248.y;
        float _S252 = _S248.z;
        float norm_0 = (F32_sqrt((_S250 * _S250 + _S251 * _S251 + _S252 * _S252)));
        float x_19 = _S250 / norm_0;
        float y_7 = _S251 / norm_0;
        float z_4 = _S252 / norm_0;
        float3  _S253 = _S249 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_0)[int(1)] + make_float3 (z_4) * (*sh_coeffs_0)[int(2)] - make_float3 (x_19) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S253;
        float z2_4 = z_4 * z_4;
        float fTmp0B_0 = -1.09254848957061768f * z_4;
        float fC1_0 = x_19 * x_19 - y_7 * y_7;
        float fS1_0 = 2.0f * x_19 * y_7;
        float3  _S254 = _S253 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_7) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_4 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_19) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S254;
        float fTmp0C_0 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_4;
        *rgb_0 = max_0(_S254 + (make_float3 (-0.59004360437393188f * (x_19 * fS1_0 + y_7 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_7) * (*sh_coeffs_0)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_4 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_19) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_19 * fC1_0 - y_7 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_5, float fy_5, float cx_5, float cy_5, float4  radial_coeffs_8, float2  tangential_coeffs_8, float2  thin_prism_coeffs_8, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_5, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_5, mean_1) + t_4;
        float _S255 = length_1(mean_c_1);
        bool _S256;
        if(_S255 < near_plane_1)
        {
            _S256 = true;
        }
        else
        {
            _S256 = _S255 > far_plane_1;
        }
        if(_S256)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S257 = exp_0(scale_3);
        float x_20 = quat_4.y;
        float inv_norm_4 = (F32_rsqrt((x_20 * x_20 + quat_4.z * quat_4.z + quat_4.w * quat_4.w + quat_4.x * quat_4.x)));
        float x_21 = quat_4.y * inv_norm_4;
        float y_8 = quat_4.z * inv_norm_4;
        float z_5 = quat_4.w * inv_norm_4;
        float w_4 = quat_4.x * inv_norm_4;
        float x2_4 = x_21 * x_21;
        float y2_4 = y_8 * y_8;
        float z2_5 = z_5 * z_5;
        float xy_4 = x_21 * y_8;
        float xz_4 = x_21 * z_5;
        float yz_4 = y_8 * z_5;
        float wx_4 = w_4 * x_21;
        float wy_4 = w_4 * y_8;
        float wz_4 = w_4 * z_5;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_5), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_5), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S257.x, 0.0f, 0.0f, 0.0f, _S257.y, 0.0f, 0.0f, 0.0f, _S257.z));
        Matrix<float, 3, 3>  _S258 = transpose_0(R_5);
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_5, mul_4(M_3, transpose_0(M_3))), _S258);
        Matrix<float, 2, 2>  covar2d_1;
        fisheye_proj_3dgs(mean_c_1, covar_c_1, fx_5, fy_5, cx_5, cy_5, radial_coeffs_8, tangential_coeffs_8, thin_prism_coeffs_8, &covar2d_1, mean2d_5);
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S259 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S259;
        float _S260 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S260;
        float det_blur_2 = _S259 * _S260 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S261 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S256 = true;
        }
        else
        {
            _S256 = xmin_1 >= float(image_width_1);
        }
        if(_S256)
        {
            _S256 = true;
        }
        else
        {
            _S256 = ymax_1 <= 0.0f;
        }
        if(_S256)
        {
            _S256 = true;
        }
        else
        {
            _S256 = ymin_1 >= float(image_height_1);
        }
        if(_S256)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = _S255;
        *conic_1 = make_float3 (_S261.rows[int(0)].x, _S261.rows[int(0)].y, _S261.rows[int(1)].y);
        float3  _S262 = mean_1 - - mul_0(_S258, t_4);
        float3  _S263 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S263;
        float _S264 = _S262.x;
        float _S265 = _S262.y;
        float _S266 = _S262.z;
        float norm_1 = (F32_sqrt((_S264 * _S264 + _S265 * _S265 + _S266 * _S266)));
        float x_22 = _S264 / norm_1;
        float y_9 = _S265 / norm_1;
        float z_6 = _S266 / norm_1;
        float3  _S267 = _S263 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_1)[int(1)] + make_float3 (z_6) * (*sh_coeffs_1)[int(2)] - make_float3 (x_22) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S267;
        float z2_6 = z_6 * z_6;
        float fTmp0B_1 = -1.09254848957061768f * z_6;
        float fC1_1 = x_22 * x_22 - y_9 * y_9;
        float fS1_1 = 2.0f * x_22 * y_9;
        float3  _S268 = _S267 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_9) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_6 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_22) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S268;
        float fTmp0C_1 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_6;
        *rgb_1 = max_0(_S268 + (make_float3 (-0.59004360437393188f * (x_22 * fS1_1 + y_9 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_9) * (*sh_coeffs_1)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_6 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_22) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_22 * fC1_1 - y_9 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_6, float fy_6, float cx_6, float cy_6, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_6, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_6, mean_2) + t_5;
        float _S269 = mean_c_2.z;
        bool _S270;
        if(_S269 < near_plane_2)
        {
            _S270 = true;
        }
        else
        {
            _S270 = _S269 > far_plane_2;
        }
        if(_S270)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S271 = exp_0(scale_4);
        float x_23 = quat_5.y;
        float inv_norm_5 = (F32_rsqrt((x_23 * x_23 + quat_5.z * quat_5.z + quat_5.w * quat_5.w + quat_5.x * quat_5.x)));
        float x_24 = quat_5.y * inv_norm_5;
        float y_10 = quat_5.z * inv_norm_5;
        float z_7 = quat_5.w * inv_norm_5;
        float w_5 = quat_5.x * inv_norm_5;
        float x2_5 = x_24 * x_24;
        float y2_5 = y_10 * y_10;
        float z2_7 = z_7 * z_7;
        float xy_5 = x_24 * y_10;
        float xz_5 = x_24 * z_7;
        float yz_5 = y_10 * z_7;
        float wx_5 = w_5 * x_24;
        float wy_5 = w_5 * y_10;
        float wz_5 = w_5 * z_7;
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_7), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_7), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S271.x, 0.0f, 0.0f, 0.0f, _S271.y, 0.0f, 0.0f, 0.0f, _S271.z));
        Matrix<float, 3, 3>  _S272 = transpose_0(R_6);
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_6, mul_4(M_4, transpose_0(M_4))), _S272);
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_6, covar_c_2), transpose_1(J_6));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S273 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S273;
        float _S274 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S274;
        float det_blur_3 = _S273 * _S274 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S275 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S270 = true;
        }
        else
        {
            _S270 = xmin_2 >= float(image_width_2);
        }
        if(_S270)
        {
            _S270 = true;
        }
        else
        {
            _S270 = ymax_2 <= 0.0f;
        }
        if(_S270)
        {
            _S270 = true;
        }
        else
        {
            _S270 = ymin_2 >= float(image_height_2);
        }
        if(_S270)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = length_1(mean_c_2);
        *conic_2 = make_float3 (_S275.rows[int(0)].x, _S275.rows[int(0)].y, _S275.rows[int(1)].y);
        float3  _S276 = mean_2 - - mul_0(_S272, t_5);
        float3  _S277 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S277;
        float _S278 = _S276.x;
        float _S279 = _S276.y;
        float _S280 = _S276.z;
        float norm_2 = (F32_sqrt((_S278 * _S278 + _S279 * _S279 + _S280 * _S280)));
        float x_25 = _S278 / norm_2;
        float y_11 = _S279 / norm_2;
        float z_8 = _S280 / norm_2;
        float3  _S281 = _S277 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_2)[int(1)] + make_float3 (z_8) * (*sh_coeffs_2)[int(2)] - make_float3 (x_25) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S281;
        float z2_8 = z_8 * z_8;
        float fTmp0B_2 = -1.09254848957061768f * z_8;
        float fC1_2 = x_25 * x_25 - y_11 * y_11;
        float fS1_2 = 2.0f * x_25 * y_11;
        float3  _S282 = _S281 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_11) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_8 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_25) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S282;
        float fTmp0C_2 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_8;
        *rgb_2 = max_0(_S282 + (make_float3 (-0.59004360437393188f * (x_25 * fS1_2 + y_11 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_11) * (*sh_coeffs_2)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_8 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_25) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_25 * fC1_2 - y_11 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_7, float fy_7, float cx_7, float cy_7, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_7, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    float3  mean_c_3 = mul_0(R_7, mean_3) + t_6;
    float3  _S283 = exp_0(scale_5);
    float x_26 = quat_6.y;
    float inv_norm_6 = (F32_rsqrt((x_26 * x_26 + quat_6.z * quat_6.z + quat_6.w * quat_6.w + quat_6.x * quat_6.x)));
    float x_27 = quat_6.y * inv_norm_6;
    float y_12 = quat_6.z * inv_norm_6;
    float z_9 = quat_6.w * inv_norm_6;
    float w_6 = quat_6.x * inv_norm_6;
    float x2_6 = x_27 * x_27;
    float y2_6 = y_12 * y_12;
    float z2_9 = z_9 * z_9;
    float xy_6 = x_27 * y_12;
    float xz_6 = x_27 * z_9;
    float yz_6 = y_12 * z_9;
    float wx_6 = w_6 * x_27;
    float wy_6 = w_6 * y_12;
    float wz_6 = w_6 * z_9;
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_9), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_9), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S283.x, 0.0f, 0.0f, 0.0f, _S283.y, 0.0f, 0.0f, 0.0f, _S283.z));
    Matrix<float, 3, 3>  _S284 = transpose_0(R_7);
    float _S285 = float(image_width_3);
    float _S286 = float(image_height_3);
    float _S287 = 0.30000001192092896f * (0.5f * _S285 / fx_7);
    float _S288 = 0.30000001192092896f * (0.5f * _S286 / fy_7);
    float rz_3 = 1.0f / mean_c_3.z;
    float rz2_3 = rz_3 * rz_3;
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, - fx_7 * (mean_c_3.z * (F32_min(((_S285 - cx_7) / fx_7 + _S287), ((F32_max((- (cx_7 / fx_7 + _S287)), (mean_c_3.x * rz_3))))))) * rz2_3, 0.0f, fy_7 * rz_3, - fy_7 * (mean_c_3.z * (F32_min(((_S286 - cy_7) / fy_7 + _S288), ((F32_max((- (cy_7 / fy_7 + _S288)), (mean_c_3.y * rz_3))))))) * rz2_3);
    Matrix<float, 2, 2>  covar2d_3 = mul_6(mul_5(J_7, mul_4(mul_4(R_7, mul_4(M_5, transpose_0(M_5))), _S284)), transpose_1(J_7));
    *mean2d_7 = make_float2 (fx_7 * mean_c_3.x * rz_3 + cx_7, fy_7 * mean_c_3.y * rz_3 + cy_7);
    float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
    float _S289 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_3)->rows + (int(0)))->x) = _S289;
    float _S290 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_3)->rows + (int(1)))->y) = _S290;
    float compensation_4 = (F32_sqrt(((F32_max((0.0f), (det_orig_4 / (_S289 * _S290 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_3.rows[int(0)].x * covar2d_3.rows[int(1)].y - covar2d_3.rows[int(0)].y * covar2d_3.rows[int(1)].x);
    Matrix<float, 2, 2>  _S291 = makeMatrix<float, 2, 2> (covar2d_3.rows[int(1)].y * invdet_4, - covar2d_3.rows[int(0)].y * invdet_4, - covar2d_3.rows[int(1)].x * invdet_4, covar2d_3.rows[int(0)].x * invdet_4);
    *opacity_3 = 1.0f / (1.0f + (F32_exp((- in_opacity_3))));
    if(antialiased_3)
    {
        *opacity_3 = *opacity_3 * compensation_4;
    }
    float extend_3 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_3 / 0.00392156885936856f)))))))));
    float radius_x_3 = extend_3 * (F32_sqrt((covar2d_3[int(0)].x)));
    float radius_y_3 = extend_3 * (F32_sqrt((covar2d_3[int(1)].y)));
    *aabb_xyxy_3 = make_int4 (int((F32_floor(((*mean2d_7).x - radius_x_3)))), int((F32_floor(((*mean2d_7).y - radius_y_3)))), int((F32_ceil(((*mean2d_7).x + radius_x_3)))), int((F32_ceil(((*mean2d_7).y + radius_y_3)))));
    *depth_3 = length_1(mean_c_3);
    *conic_3 = make_float3 (_S291.rows[int(0)].x, _S291.rows[int(0)].y, _S291.rows[int(1)].y);
    float3  _S292 = mean_3 - - mul_0(_S284, t_6);
    float3  _S293 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
    *rgb_3 = _S293;
    float _S294 = _S292.x;
    float _S295 = _S292.y;
    float _S296 = _S292.z;
    float norm_3 = (F32_sqrt((_S294 * _S294 + _S295 * _S295 + _S296 * _S296)));
    float x_28 = _S294 / norm_3;
    float y_13 = _S295 / norm_3;
    float z_10 = _S296 / norm_3;
    float3  _S297 = _S293 + make_float3 (0.48860251903533936f) * (make_float3 (- y_13) * (*sh_coeffs_3)[int(1)] + make_float3 (z_10) * (*sh_coeffs_3)[int(2)] - make_float3 (x_28) * (*sh_coeffs_3)[int(3)]);
    *rgb_3 = _S297;
    float z2_10 = z_10 * z_10;
    float fTmp0B_3 = -1.09254848957061768f * z_10;
    float fC1_3 = x_28 * x_28 - y_13 * y_13;
    float fS1_3 = 2.0f * x_28 * y_13;
    float3  _S298 = _S297 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_13) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_10 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_28) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
    *rgb_3 = _S298;
    float fTmp0C_3 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
    float fTmp1B_3 = 1.44530570507049561f * z_10;
    *rgb_3 = max_0(_S298 + (make_float3 (-0.59004360437393188f * (x_28 * fS1_3 + y_13 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_13) * (*sh_coeffs_3)[int(11)] + make_float3 (z_10 * (1.86588168144226074f * z2_10 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_28) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_28 * fC1_3 - y_13 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_8, float fy_8, float cx_8, float cy_8, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_8, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    float3  mean_c_4 = mul_0(R_8, mean_4) + t_7;
    float3  _S299 = exp_0(scale_6);
    float x_29 = quat_7.y;
    float inv_norm_7 = (F32_rsqrt((x_29 * x_29 + quat_7.z * quat_7.z + quat_7.w * quat_7.w + quat_7.x * quat_7.x)));
    float x_30 = quat_7.y * inv_norm_7;
    float y_14 = quat_7.z * inv_norm_7;
    float z_11 = quat_7.w * inv_norm_7;
    float w_7 = quat_7.x * inv_norm_7;
    float x2_7 = x_30 * x_30;
    float y2_7 = y_14 * y_14;
    float z2_11 = z_11 * z_11;
    float xy_7 = x_30 * y_14;
    float xz_7 = x_30 * z_11;
    float yz_7 = y_14 * z_11;
    float wx_7 = w_7 * x_30;
    float wy_7 = w_7 * y_14;
    float wz_7 = w_7 * z_11;
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_11), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_11), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S299.x, 0.0f, 0.0f, 0.0f, _S299.y, 0.0f, 0.0f, 0.0f, _S299.z));
    Matrix<float, 3, 3>  _S300 = transpose_0(R_8);
    Matrix<float, 2, 2>  covar2d_4;
    fisheye_proj_3dgs(mean_c_4, mul_4(mul_4(R_8, mul_4(M_6, transpose_0(M_6))), _S300), fx_8, fy_8, cx_8, cy_8, radial_coeffs_11, tangential_coeffs_11, thin_prism_coeffs_11, &covar2d_4, mean2d_8);
    float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
    float _S301 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_4)->rows + (int(0)))->x) = _S301;
    float _S302 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_4)->rows + (int(1)))->y) = _S302;
    float compensation_5 = (F32_sqrt(((F32_max((0.0f), (det_orig_5 / (_S301 * _S302 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_4.rows[int(0)].x * covar2d_4.rows[int(1)].y - covar2d_4.rows[int(0)].y * covar2d_4.rows[int(1)].x);
    Matrix<float, 2, 2>  _S303 = makeMatrix<float, 2, 2> (covar2d_4.rows[int(1)].y * invdet_5, - covar2d_4.rows[int(0)].y * invdet_5, - covar2d_4.rows[int(1)].x * invdet_5, covar2d_4.rows[int(0)].x * invdet_5);
    *opacity_4 = 1.0f / (1.0f + (F32_exp((- in_opacity_4))));
    if(antialiased_4)
    {
        *opacity_4 = *opacity_4 * compensation_5;
    }
    float extend_4 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_4 / 0.00392156885936856f)))))))));
    float radius_x_4 = extend_4 * (F32_sqrt((covar2d_4[int(0)].x)));
    float radius_y_4 = extend_4 * (F32_sqrt((covar2d_4[int(1)].y)));
    *aabb_xyxy_4 = make_int4 (int((F32_floor(((*mean2d_8).x - radius_x_4)))), int((F32_floor(((*mean2d_8).y - radius_y_4)))), int((F32_ceil(((*mean2d_8).x + radius_x_4)))), int((F32_ceil(((*mean2d_8).y + radius_y_4)))));
    *depth_4 = length_1(mean_c_4);
    *conic_4 = make_float3 (_S303.rows[int(0)].x, _S303.rows[int(0)].y, _S303.rows[int(1)].y);
    float3  _S304 = mean_4 - - mul_0(_S300, t_7);
    float3  _S305 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
    *rgb_4 = _S305;
    float _S306 = _S304.x;
    float _S307 = _S304.y;
    float _S308 = _S304.z;
    float norm_4 = (F32_sqrt((_S306 * _S306 + _S307 * _S307 + _S308 * _S308)));
    float x_31 = _S306 / norm_4;
    float y_15 = _S307 / norm_4;
    float z_12 = _S308 / norm_4;
    float3  _S309 = _S305 + make_float3 (0.48860251903533936f) * (make_float3 (- y_15) * (*sh_coeffs_4)[int(1)] + make_float3 (z_12) * (*sh_coeffs_4)[int(2)] - make_float3 (x_31) * (*sh_coeffs_4)[int(3)]);
    *rgb_4 = _S309;
    float z2_12 = z_12 * z_12;
    float fTmp0B_4 = -1.09254848957061768f * z_12;
    float fC1_4 = x_31 * x_31 - y_15 * y_15;
    float fS1_4 = 2.0f * x_31 * y_15;
    float3  _S310 = _S309 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_15) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_31) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
    *rgb_4 = _S310;
    float fTmp0C_4 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
    float fTmp1B_4 = 1.44530570507049561f * z_12;
    *rgb_4 = max_0(_S310 + (make_float3 (-0.59004360437393188f * (x_31 * fS1_4 + y_15 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_15) * (*sh_coeffs_4)[int(11)] + make_float3 (z_12 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_31) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_31 * fC1_4 - y_15 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_9, float fy_9, float cx_9, float cy_9, float4  radial_coeffs_12, float2  tangential_coeffs_12, float2  thin_prism_coeffs_12, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_9, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_9, mean_5) + t_8;
    float3  _S311 = exp_0(scale_7);
    float x_32 = quat_8.y;
    float inv_norm_8 = (F32_rsqrt((x_32 * x_32 + quat_8.z * quat_8.z + quat_8.w * quat_8.w + quat_8.x * quat_8.x)));
    float x_33 = quat_8.y * inv_norm_8;
    float y_16 = quat_8.z * inv_norm_8;
    float z_13 = quat_8.w * inv_norm_8;
    float w_8 = quat_8.x * inv_norm_8;
    float x2_8 = x_33 * x_33;
    float y2_8 = y_16 * y_16;
    float z2_13 = z_13 * z_13;
    float xy_8 = x_33 * y_16;
    float xz_8 = x_33 * z_13;
    float yz_8 = y_16 * z_13;
    float wx_8 = w_8 * x_33;
    float wy_8 = w_8 * y_16;
    float wz_8 = w_8 * z_13;
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_13), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_13), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S311.x, 0.0f, 0.0f, 0.0f, _S311.y, 0.0f, 0.0f, 0.0f, _S311.z));
    Matrix<float, 3, 3>  _S312 = transpose_0(R_9);
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_8, mul_4(mul_4(R_9, mul_4(M_7, transpose_0(M_7))), _S312)), transpose_1(J_8));
    *mean2d_9 = make_float2 (fx_9 * mean_c_5.x + cx_9, fy_9 * mean_c_5.y + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S313 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S313;
    float _S314 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S314;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S313 * _S314 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S315 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_6, - covar2d_5.rows[int(0)].y * invdet_6, - covar2d_5.rows[int(1)].x * invdet_6, covar2d_5.rows[int(0)].x * invdet_6);
    *opacity_5 = 1.0f / (1.0f + (F32_exp((- in_opacity_5))));
    if(antialiased_5)
    {
        *opacity_5 = *opacity_5 * compensation_6;
    }
    float extend_5 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_5 / 0.00392156885936856f)))))))));
    float radius_x_5 = extend_5 * (F32_sqrt((covar2d_5[int(0)].x)));
    float radius_y_5 = extend_5 * (F32_sqrt((covar2d_5[int(1)].y)));
    *aabb_xyxy_5 = make_int4 (int((F32_floor(((*mean2d_9).x - radius_x_5)))), int((F32_floor(((*mean2d_9).y - radius_y_5)))), int((F32_ceil(((*mean2d_9).x + radius_x_5)))), int((F32_ceil(((*mean2d_9).y + radius_y_5)))));
    *depth_5 = length_1(mean_c_5);
    *conic_5 = make_float3 (_S315.rows[int(0)].x, _S315.rows[int(0)].y, _S315.rows[int(1)].y);
    float3  _S316 = mean_5 - - mul_0(_S312, t_8);
    float3  _S317 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S317;
    float _S318 = _S316.x;
    float _S319 = _S316.y;
    float _S320 = _S316.z;
    float norm_5 = (F32_sqrt((_S318 * _S318 + _S319 * _S319 + _S320 * _S320)));
    float x_34 = _S318 / norm_5;
    float y_17 = _S319 / norm_5;
    float z_14 = _S320 / norm_5;
    float3  _S321 = _S317 + make_float3 (0.48860251903533936f) * (make_float3 (- y_17) * (*sh_coeffs_5)[int(1)] + make_float3 (z_14) * (*sh_coeffs_5)[int(2)] - make_float3 (x_34) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S321;
    float z2_14 = z_14 * z_14;
    float fTmp0B_5 = -1.09254848957061768f * z_14;
    float fC1_5 = x_34 * x_34 - y_17 * y_17;
    float fS1_5 = 2.0f * x_34 * y_17;
    float3  _S322 = _S321 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_17) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_14 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_34) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S322;
    float fTmp0C_5 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_14;
    *rgb_5 = max_0(_S322 + (make_float3 (-0.59004360437393188f * (x_34 * fS1_5 + y_17 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_17) * (*sh_coeffs_5)[int(11)] + make_float3 (z_14 * (1.86588168144226074f * z2_14 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_34) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_34 * fC1_5 - y_17 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S323, float3  _S324)
{
    return mul_0(_S323, _S324);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S325)
{
    return exp_0(_S325);
}

inline __device__ float s_primal_ctx_rsqrt_0(float _S326)
{
    return (F32_rsqrt((_S326)));
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S327, Matrix<float, 3, 3>  _S328)
{
    return mul_4(_S327, _S328);
}

inline __device__ float s_primal_ctx_max_0(float _S329, float _S330)
{
    return (F32_max((_S329), (_S330)));
}

inline __device__ float s_primal_ctx_min_0(float _S331, float _S332)
{
    return (F32_min((_S331), (_S332)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S333, Matrix<float, 3, 3>  _S334)
{
    return mul_5(_S333, _S334);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_4(Matrix<float, 2, 3>  _S335, Matrix<float, 3, 2>  _S336)
{
    return mul_6(_S335, _S336);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S337)
{
    return (F32_sqrt((_S337)));
}

inline __device__ float s_primal_ctx_exp_1(float _S338)
{
    return (F32_exp((_S338)));
}

inline __device__ float s_primal_ctx_log_0(float _S339)
{
    return (F32_log((_S339)));
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S340, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S341, float3  _S342)
{
    _d_max_vector_0(_S340, _S341, _S342);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S343, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S344, float3  _S345)
{
    _d_mul_0(_S343, _S344, _S345);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S346, DiffPair_float_0 * _S347, float _S348)
{
    _d_min_0(_S346, _S347, _S348);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S349, float _S350)
{
    _d_log_0(_S349, _S350);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S351, float _S352)
{
    _d_exp_0(_S351, _S352);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S353, DiffPair_float_0 * _S354, float _S355)
{
    _d_max_0(_S353, _S354, _S355);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S356, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S357, Matrix<float, 2, 2>  _S358)
{
    mul_3(_S356, _S357, _S358);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S359, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S360, Matrix<float, 2, 3>  _S361)
{
    mul_2(_S359, _S360, _S361);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S362, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S363, Matrix<float, 3, 3>  _S364)
{
    mul_1(_S362, _S363, _S364);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S365, float _S366)
{
    _d_rsqrt_0(_S365, _S366);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S367, float3  _S368)
{
    _d_exp_vector_0(_S367, _S368);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_10, float fy_10, float cx_10, float cy_10, float4  radial_coeffs_13, float2  tangential_coeffs_13, float2  thin_prism_coeffs_13, uint image_width_6, uint image_height_6, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_6 = s_primal_ctx_mul_1(R_10, mean_6) + t_9;
    float3  _S369 = s_primal_ctx_exp_0(scale_8);
    float _S370 = quat_9.y;
    float _S371 = _S370 * _S370 + quat_9.z * quat_9.z + quat_9.w * quat_9.w + quat_9.x * quat_9.x;
    float _S372 = s_primal_ctx_rsqrt_0(_S371);
    float x_35 = quat_9.y * _S372;
    float y_18 = quat_9.z * _S372;
    float z_15 = quat_9.w * _S372;
    float w_9 = quat_9.x * _S372;
    float x2_9 = x_35 * x_35;
    float y2_9 = y_18 * y_18;
    float z2_15 = z_15 * z_15;
    float xy_9 = x_35 * y_18;
    float xz_9 = x_35 * z_15;
    float yz_9 = y_18 * z_15;
    float wx_9 = w_9 * x_35;
    float wy_9 = w_9 * y_18;
    float wz_9 = w_9 * z_15;
    Matrix<float, 3, 3>  _S373 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S369.x, 0.0f, 0.0f, 0.0f, _S369.y, 0.0f, 0.0f, 0.0f, _S369.z);
    Matrix<float, 3, 3>  _S374 = s_primal_ctx_mul_2(_S373, S_0);
    Matrix<float, 3, 3>  _S375 = transpose_0(_S374);
    Matrix<float, 3, 3>  _S376 = s_primal_ctx_mul_2(_S374, _S375);
    Matrix<float, 3, 3>  _S377 = s_primal_ctx_mul_2(R_10, _S376);
    Matrix<float, 3, 3>  _S378 = transpose_0(R_10);
    Matrix<float, 3, 3>  _S379 = s_primal_ctx_mul_2(_S377, _S378);
    float _S380 = float(image_width_6);
    float _S381 = float(image_height_6);
    float _S382 = 0.30000001192092896f * (0.5f * _S380 / fx_10);
    float lim_x_pos_0 = (_S380 - cx_10) / fx_10 + _S382;
    float _S383 = 0.30000001192092896f * (0.5f * _S381 / fy_10);
    float lim_y_pos_0 = (_S381 - cy_10) / fy_10 + _S383;
    float rz_4 = 1.0f / mean_c_6.z;
    float _S384 = mean_c_6.z * mean_c_6.z;
    float rz2_4 = rz_4 * rz_4;
    float _S385 = - (cx_10 / fx_10 + _S382);
    float _S386 = mean_c_6.x * rz_4;
    float _S387 = s_primal_ctx_max_0(_S385, _S386);
    float _S388 = s_primal_ctx_min_0(lim_x_pos_0, _S387);
    float _S389 = - (cy_10 / fy_10 + _S383);
    float _S390 = mean_c_6.y * rz_4;
    float _S391 = s_primal_ctx_max_0(_S389, _S390);
    float _S392 = s_primal_ctx_min_0(lim_y_pos_0, _S391);
    float _S393 = - fx_10;
    float _S394 = _S393 * (mean_c_6.z * _S388);
    float _S395 = - fy_10;
    float _S396 = _S395 * (mean_c_6.z * _S392);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_10 * rz_4, 0.0f, _S394 * rz2_4, 0.0f, fy_10 * rz_4, _S396 * rz2_4);
    Matrix<float, 2, 3>  _S397 = s_primal_ctx_mul_3(J_9, _S379);
    Matrix<float, 3, 2>  _S398 = transpose_1(J_9);
    Matrix<float, 2, 2>  _S399 = s_primal_ctx_mul_4(_S397, _S398);
    float _S400 = fx_10 * mean_c_6.x;
    float _S401 = fy_10 * mean_c_6.y;
    float _S402 = _S399.rows[int(0)].y * _S399.rows[int(1)].x;
    float det_orig_7 = _S399.rows[int(0)].x * _S399.rows[int(1)].y - _S402;
    float _S403 = _S399.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S404 = _S399;
    *&(((&_S404)->rows + (int(0)))->x) = _S403;
    float _S405 = _S399.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S404)->rows + (int(1)))->y) = _S405;
    Matrix<float, 2, 2>  _S406 = _S404;
    Matrix<float, 2, 2>  _S407 = _S404;
    float det_blur_4 = _S403 * _S405 - _S402;
    float _S408 = det_orig_7 / det_blur_4;
    float _S409 = det_blur_4 * det_blur_4;
    float _S410 = s_primal_ctx_max_0(0.0f, _S408);
    float _S411 = s_primal_ctx_sqrt_0(_S410);
    float invdet_7 = 1.0f / det_blur_4;
    float _S412 = - _S399.rows[int(0)].y;
    float _S413 = - _S399.rows[int(1)].x;
    float _S414 = - in_opacity_6;
    float _S415 = 1.0f + s_primal_ctx_exp_1(_S414);
    float _S416 = 1.0f / _S415;
    float _S417 = _S415 * _S415;
    float _S418;
    if(antialiased_6)
    {
        _S418 = _S416 * _S411;
    }
    else
    {
        _S418 = _S416;
    }
    float _S419 = _S418 / 0.00392156885936856f;
    float _S420 = 2.0f * s_primal_ctx_log_0(_S419);
    float _S421 = s_primal_ctx_sqrt_0(_S420);
    float _S422 = _S406.rows[int(0)].x;
    float _S423 = _S407.rows[int(1)].y;
    float3  _S424 = mean_6 - - s_primal_ctx_mul_1(_S378, t_9);
    float _S425 = _S424.x;
    float _S426 = _S424.y;
    float _S427 = _S424.z;
    float _S428 = _S425 * _S425 + _S426 * _S426 + _S427 * _S427;
    float _S429 = s_primal_ctx_sqrt_0(_S428);
    float x_36 = _S425 / _S429;
    float3  _S430 = make_float3 (x_36);
    float _S431 = _S429 * _S429;
    float y_19 = _S426 / _S429;
    float z_16 = _S427 / _S429;
    float3  _S432 = make_float3 (z_16);
    float _S433 = - y_19;
    float3  _S434 = make_float3 (_S433);
    float z2_16 = z_16 * z_16;
    float fTmp0B_6 = -1.09254848957061768f * z_16;
    float fC1_6 = x_36 * x_36 - y_19 * y_19;
    float _S435 = 2.0f * x_36;
    float fS1_6 = _S435 * y_19;
    float pSH6_0 = 0.94617468118667603f * z2_16 - 0.31539157032966614f;
    float3  _S436 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_6 * x_36;
    float3  _S437 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_6 * y_19;
    float3  _S438 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_6;
    float3  _S439 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_6;
    float3  _S440 = make_float3 (pSH4_0);
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_16;
    float _S441 = 1.86588168144226074f * z2_16 - 1.11952900886535645f;
    float pSH12_0 = z_16 * _S441;
    float3  _S442 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_6 * x_36;
    float3  _S443 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_6 * y_19;
    float3  _S444 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_6 * fC1_6;
    float3  _S445 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_6 * fS1_6;
    float3  _S446 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_36 * fC1_6 - y_19 * fS1_6);
    float3  _S447 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_36 * fS1_6 + y_19 * fC1_6);
    float3  _S448 = make_float3 (pSH9_0);
    float3  _S449 = make_float3 (0.0f);
    float3  _S450 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S451;
    (&_S451)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S433) * (*sh_coeffs_6)[int(1)] + make_float3 (z_16) * (*sh_coeffs_6)[int(2)] - make_float3 (x_36) * (*sh_coeffs_6)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_6)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_6)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_6)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_6)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_6)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_6)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_6)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_6)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_6)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_6)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_6)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f);
    (&_S451)->differential_0 = _S450;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S452;
    (&_S452)->primal_0 = _S449;
    (&_S452)->differential_0 = _S450;
    s_bwd_prop_max_0(&_S451, &_S452, v_rgb_0);
    float3  _S453 = _S447 * _S451.differential_0;
    float3  _S454 = (*sh_coeffs_6)[int(15)] * _S451.differential_0;
    float3  _S455 = _S445 * _S451.differential_0;
    float3  _S456 = (*sh_coeffs_6)[int(14)] * _S451.differential_0;
    float3  _S457 = _S443 * _S451.differential_0;
    float3  _S458 = (*sh_coeffs_6)[int(13)] * _S451.differential_0;
    float3  _S459 = _S442 * _S451.differential_0;
    float3  _S460 = (*sh_coeffs_6)[int(12)] * _S451.differential_0;
    float3  _S461 = _S444 * _S451.differential_0;
    float3  _S462 = (*sh_coeffs_6)[int(11)] * _S451.differential_0;
    float3  _S463 = _S446 * _S451.differential_0;
    float3  _S464 = (*sh_coeffs_6)[int(10)] * _S451.differential_0;
    float3  _S465 = _S448 * _S451.differential_0;
    float3  _S466 = (*sh_coeffs_6)[int(9)] * _S451.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S466.x + _S466.y + _S466.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S454.x + _S454.y + _S454.z);
    float _S467 = _S464.x + _S464.y + _S464.z;
    float _S468 = _S456.x + _S456.y + _S456.z;
    float _S469 = _S462.x + _S462.y + _S462.z;
    float _S470 = _S458.x + _S458.y + _S458.z;
    float _S471 = _S460.x + _S460.y + _S460.z;
    float _S472 = - s_diff_fC2_T_0;
    float3  _S473 = _S439 * _S451.differential_0;
    float3  _S474 = (*sh_coeffs_6)[int(8)] * _S451.differential_0;
    float3  _S475 = _S437 * _S451.differential_0;
    float3  _S476 = (*sh_coeffs_6)[int(7)] * _S451.differential_0;
    float3  _S477 = _S436 * _S451.differential_0;
    float3  _S478 = (*sh_coeffs_6)[int(6)] * _S451.differential_0;
    float3  _S479 = _S438 * _S451.differential_0;
    float3  _S480 = (*sh_coeffs_6)[int(5)] * _S451.differential_0;
    float3  _S481 = _S440 * _S451.differential_0;
    float3  _S482 = (*sh_coeffs_6)[int(4)] * _S451.differential_0;
    float _S483 = _S480.x + _S480.y + _S480.z;
    float _S484 = _S476.x + _S476.y + _S476.z;
    float _S485 = fTmp1B_6 * _S467 + x_36 * s_diff_fS2_T_0 + y_19 * _S472 + 0.54627424478530884f * (_S482.x + _S482.y + _S482.z);
    float _S486 = fTmp1B_6 * _S468 + y_19 * s_diff_fS2_T_0 + x_36 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S474.x + _S474.y + _S474.z);
    float _S487 = y_19 * - _S486;
    float _S488 = x_36 * _S486;
    float _S489 = z_16 * (1.86588168144226074f * (z_16 * _S471) + -2.28522896766662598f * (y_19 * _S469 + x_36 * _S470) + 0.94617468118667603f * (_S478.x + _S478.y + _S478.z));
    float3  _S490 = make_float3 (0.48860251903533936f) * _S451.differential_0;
    float3  _S491 = - _S490;
    float3  _S492 = _S430 * _S491;
    float3  _S493 = (*sh_coeffs_6)[int(3)] * _S491;
    float3  _S494 = _S432 * _S490;
    float3  _S495 = (*sh_coeffs_6)[int(2)] * _S490;
    float3  _S496 = _S434 * _S490;
    float3  _S497 = (*sh_coeffs_6)[int(1)] * _S490;
    float _S498 = (_S441 * _S471 + 1.44530570507049561f * (fS1_6 * _S467 + fC1_6 * _S468) + -1.09254848957061768f * (y_19 * _S483 + x_36 * _S484) + _S489 + _S489 + _S495.x + _S495.y + _S495.z) / _S431;
    float _S499 = _S429 * _S498;
    float _S500 = (fTmp0C_6 * _S469 + fC1_6 * s_diff_fS2_T_0 + fS1_6 * _S472 + fTmp0B_6 * _S483 + _S435 * _S485 + _S487 + _S487 + - (_S497.x + _S497.y + _S497.z)) / _S431;
    float _S501 = _S429 * _S500;
    float _S502 = (fTmp0C_6 * _S470 + fS1_6 * s_diff_fS2_T_0 + fC1_6 * s_diff_fC2_T_0 + fTmp0B_6 * _S484 + 2.0f * (y_19 * _S485) + _S488 + _S488 + _S493.x + _S493.y + _S493.z) / _S431;
    float _S503 = _S429 * _S502;
    float _S504 = _S427 * - _S498 + _S426 * - _S500 + _S425 * - _S502;
    DiffPair_float_0 _S505;
    (&_S505)->primal_0 = _S428;
    (&_S505)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S505, _S504);
    float _S506 = _S427 * _S505.differential_0;
    float _S507 = _S426 * _S505.differential_0;
    float _S508 = _S425 * _S505.differential_0;
    float3  _S509 = make_float3 (0.282094806432724f) * _S451.differential_0;
    float3  _S510 = make_float3 (_S503 + _S508 + _S508, _S501 + _S507 + _S507, _S499 + _S506 + _S506);
    float3  _S511 = - - _S510;
    Matrix<float, 3, 3>  _S512 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S513;
    (&_S513)->primal_0 = _S378;
    (&_S513)->differential_0 = _S512;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S514;
    (&_S514)->primal_0 = t_9;
    (&_S514)->differential_0 = _S450;
    s_bwd_prop_mul_1(&_S513, &_S514, _S511);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S515 = _S513;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S516 = _S514;
    float2  _S517 = make_float2 (0.0f);
    float2  _S518 = _S517;
    *&((&_S518)->y) = v_conic_0.z;
    float2  _S519 = _S517;
    *&((&_S519)->y) = v_conic_0.y;
    *&((&_S519)->x) = v_conic_0.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S520;
    (&_S520)->primal_0 = mean_c_6;
    (&_S520)->differential_0 = _S450;
    s_bwd_length_impl_0(&_S520, v_depth_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S521 = _S520;
    DiffPair_float_0 _S522;
    (&_S522)->primal_0 = _S423;
    (&_S522)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S522, 0.0f);
    DiffPair_float_0 _S523;
    (&_S523)->primal_0 = _S422;
    (&_S523)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S523, 0.0f);
    DiffPair_float_0 _S524;
    (&_S524)->primal_0 = 3.32999992370605469f;
    (&_S524)->differential_0 = 0.0f;
    DiffPair_float_0 _S525;
    (&_S525)->primal_0 = _S421;
    (&_S525)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S524, &_S525, 0.0f);
    DiffPair_float_0 _S526;
    (&_S526)->primal_0 = _S420;
    (&_S526)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S526, _S525.differential_0);
    float _S527 = 2.0f * _S526.differential_0;
    DiffPair_float_0 _S528;
    (&_S528)->primal_0 = _S419;
    (&_S528)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S528, _S527);
    float _S529 = v_opacity_0 + 254.9999847412109375f * _S528.differential_0;
    Matrix<float, 2, 2>  _S530 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S531 = _S530;
    _S531[int(1)] = _S518;
    _S531[int(0)] = _S519;
    Matrix<float, 2, 2>  _S532 = _S531;
    FixedArray<float3 , 16>  _S533;
    _S533[int(0)] = _S450;
    _S533[int(1)] = _S450;
    _S533[int(2)] = _S450;
    _S533[int(3)] = _S450;
    _S533[int(4)] = _S450;
    _S533[int(5)] = _S450;
    _S533[int(6)] = _S450;
    _S533[int(7)] = _S450;
    _S533[int(8)] = _S450;
    _S533[int(9)] = _S450;
    _S533[int(10)] = _S450;
    _S533[int(11)] = _S450;
    _S533[int(12)] = _S450;
    _S533[int(13)] = _S450;
    _S533[int(14)] = _S450;
    _S533[int(15)] = _S450;
    _S533[int(7)] = _S475;
    _S533[int(0)] = _S509;
    _S533[int(1)] = _S496;
    _S533[int(2)] = _S494;
    _S533[int(3)] = _S492;
    _S533[int(4)] = _S481;
    _S533[int(5)] = _S479;
    _S533[int(6)] = _S477;
    _S533[int(15)] = _S453;
    _S533[int(8)] = _S473;
    _S533[int(9)] = _S465;
    _S533[int(10)] = _S463;
    _S533[int(11)] = _S461;
    _S533[int(12)] = _S459;
    _S533[int(13)] = _S457;
    _S533[int(14)] = _S455;
    float3  _S534 = _S533[int(0)];
    float3  _S535 = _S533[int(1)];
    float3  _S536 = _S533[int(2)];
    float3  _S537 = _S533[int(3)];
    float3  _S538 = _S533[int(4)];
    float3  _S539 = _S533[int(5)];
    float3  _S540 = _S533[int(6)];
    float3  _S541 = _S533[int(7)];
    float3  _S542 = _S533[int(8)];
    float3  _S543 = _S533[int(9)];
    float3  _S544 = _S533[int(10)];
    float3  _S545 = _S533[int(11)];
    float3  _S546 = _S533[int(12)];
    float3  _S547 = _S533[int(13)];
    float3  _S548 = _S533[int(14)];
    float3  _S549 = _S533[int(15)];
    float2  _S550 = make_float2 (0.0f, _S522.differential_0);
    float2  _S551 = make_float2 (_S523.differential_0, 0.0f);
    float _S552;
    if(antialiased_6)
    {
        float _S553 = _S416 * _S529;
        _S418 = _S411 * _S529;
        _S552 = _S553;
    }
    else
    {
        _S418 = _S529;
        _S552 = 0.0f;
    }
    float _S554 = - (_S418 / _S417);
    DiffPair_float_0 _S555;
    (&_S555)->primal_0 = _S414;
    (&_S555)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S555, _S554);
    float _S556 = - _S555.differential_0;
    float _S557 = invdet_7 * _S532.rows[int(1)].y;
    float _S558 = - (invdet_7 * _S532.rows[int(1)].x);
    float _S559 = - (invdet_7 * _S532.rows[int(0)].y);
    float _S560 = invdet_7 * _S532.rows[int(0)].x;
    float _S561 = - ((_S403 * _S532.rows[int(1)].y + _S413 * _S532.rows[int(1)].x + _S412 * _S532.rows[int(0)].y + _S405 * _S532.rows[int(0)].x) / _S409);
    DiffPair_float_0 _S562;
    (&_S562)->primal_0 = _S410;
    (&_S562)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S562, _S552);
    DiffPair_float_0 _S563;
    (&_S563)->primal_0 = 0.0f;
    (&_S563)->differential_0 = 0.0f;
    DiffPair_float_0 _S564;
    (&_S564)->primal_0 = _S408;
    (&_S564)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S563, &_S564, _S562.differential_0);
    float _S565 = _S564.differential_0 / _S409;
    float s_diff_det_orig_T_0 = det_blur_4 * _S565;
    float _S566 = _S561 + det_orig_7 * - _S565;
    float _S567 = - _S566;
    float _S568 = _S403 * _S566;
    float _S569 = _S405 * _S566;
    Matrix<float, 2, 2>  _S570 = _S530;
    _S570[int(1)] = _S550;
    _S570[int(0)] = _S551;
    _S404 = _S570;
    *&(((&_S404)->rows + (int(1)))->y) = 0.0f;
    float _S571 = _S560 + _S568 + _S570.rows[int(1)].y;
    *&(((&_S404)->rows + (int(0)))->x) = 0.0f;
    float _S572 = _S557 + _S569 + _S570.rows[int(0)].x;
    float _S573 = _S567 + - s_diff_det_orig_T_0;
    float _S574 = _S558 + _S399.rows[int(0)].y * _S573;
    float _S575 = _S559 + _S399.rows[int(1)].x * _S573;
    float _S576 = _S399.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S577 = _S571 + _S399.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S578 = _S517;
    *&((&_S578)->x) = _S574;
    *&((&_S578)->y) = _S577;
    float _S579 = _S572 + _S576;
    float2  _S580 = _S517;
    *&((&_S580)->y) = _S575;
    *&((&_S580)->x) = _S579;
    float _S581 = _S401 * v_mean2d_0.y;
    float _S582 = fy_10 * (rz_4 * v_mean2d_0.y);
    float _S583 = _S400 * v_mean2d_0.x;
    float _S584 = fx_10 * (rz_4 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S585 = _S530;
    _S585[int(1)] = _S578;
    _S585[int(0)] = _S580;
    Matrix<float, 2, 2>  _S586 = _S404 + _S585;
    Matrix<float, 2, 3>  _S587 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S588;
    (&_S588)->primal_0 = _S397;
    (&_S588)->differential_0 = _S587;
    Matrix<float, 3, 2>  _S589 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S590;
    (&_S590)->primal_0 = _S398;
    (&_S590)->differential_0 = _S589;
    s_bwd_prop_mul_2(&_S588, &_S590, _S586);
    Matrix<float, 2, 3>  _S591 = transpose_2(_S590.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S592;
    (&_S592)->primal_0 = J_9;
    (&_S592)->differential_0 = _S587;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S593;
    (&_S593)->primal_0 = _S379;
    (&_S593)->differential_0 = _S512;
    s_bwd_prop_mul_3(&_S592, &_S593, _S588.differential_0);
    Matrix<float, 2, 3>  _S594 = _S591 + _S592.differential_0;
    float _S595 = _S396 * _S594.rows[int(1)].z;
    float s_diff_ty_T_0 = _S395 * (rz2_4 * _S594.rows[int(1)].z);
    float _S596 = fy_10 * _S594.rows[int(1)].y;
    float _S597 = _S394 * _S594.rows[int(0)].z;
    float s_diff_tx_T_0 = _S393 * (rz2_4 * _S594.rows[int(0)].z);
    float _S598 = fx_10 * _S594.rows[int(0)].x;
    float _S599 = mean_c_6.z * s_diff_ty_T_0;
    float _S600 = _S392 * s_diff_ty_T_0;
    DiffPair_float_0 _S601;
    (&_S601)->primal_0 = lim_y_pos_0;
    (&_S601)->differential_0 = 0.0f;
    DiffPair_float_0 _S602;
    (&_S602)->primal_0 = _S391;
    (&_S602)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S601, &_S602, _S599);
    DiffPair_float_0 _S603;
    (&_S603)->primal_0 = _S389;
    (&_S603)->differential_0 = 0.0f;
    DiffPair_float_0 _S604;
    (&_S604)->primal_0 = _S390;
    (&_S604)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S603, &_S604, _S602.differential_0);
    float _S605 = mean_c_6.y * _S604.differential_0;
    float _S606 = rz_4 * _S604.differential_0;
    float _S607 = mean_c_6.z * s_diff_tx_T_0;
    float _S608 = _S388 * s_diff_tx_T_0;
    DiffPair_float_0 _S609;
    (&_S609)->primal_0 = lim_x_pos_0;
    (&_S609)->differential_0 = 0.0f;
    DiffPair_float_0 _S610;
    (&_S610)->primal_0 = _S387;
    (&_S610)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S609, &_S610, _S607);
    DiffPair_float_0 _S611;
    (&_S611)->primal_0 = _S385;
    (&_S611)->differential_0 = 0.0f;
    DiffPair_float_0 _S612;
    (&_S612)->primal_0 = _S386;
    (&_S612)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S611, &_S612, _S610.differential_0);
    float _S613 = rz_4 * (_S595 + _S597);
    float _S614 = _S600 + _S608 + - ((_S581 + _S583 + _S596 + _S598 + _S605 + mean_c_6.x * _S612.differential_0 + _S613 + _S613) / _S384);
    float _S615 = _S582 + _S606;
    float _S616 = _S584 + rz_4 * _S612.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S617;
    (&_S617)->primal_0 = _S377;
    (&_S617)->differential_0 = _S512;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S618;
    (&_S618)->primal_0 = _S378;
    (&_S618)->differential_0 = _S512;
    s_bwd_prop_mul_4(&_S617, &_S618, _S593.differential_0);
    Matrix<float, 3, 3>  _S619 = transpose_0(_S618.differential_0 + _S515.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S620;
    (&_S620)->primal_0 = R_10;
    (&_S620)->differential_0 = _S512;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S621;
    (&_S621)->primal_0 = _S376;
    (&_S621)->differential_0 = _S512;
    s_bwd_prop_mul_4(&_S620, &_S621, _S617.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S622;
    (&_S622)->primal_0 = _S374;
    (&_S622)->differential_0 = _S512;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S623;
    (&_S623)->primal_0 = _S375;
    (&_S623)->differential_0 = _S512;
    s_bwd_prop_mul_4(&_S622, &_S623, _S621.differential_0);
    Matrix<float, 3, 3>  _S624 = _S622.differential_0 + transpose_0(_S623.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S625;
    (&_S625)->primal_0 = _S373;
    (&_S625)->differential_0 = _S512;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S626;
    (&_S626)->primal_0 = S_0;
    (&_S626)->differential_0 = _S512;
    s_bwd_prop_mul_4(&_S625, &_S626, _S624);
    Matrix<float, 3, 3>  _S627 = transpose_0(_S625.differential_0);
    float _S628 = 2.0f * - _S627.rows[int(2)].z;
    float _S629 = 2.0f * _S627.rows[int(2)].y;
    float _S630 = 2.0f * _S627.rows[int(2)].x;
    float _S631 = 2.0f * _S627.rows[int(1)].z;
    float _S632 = 2.0f * - _S627.rows[int(1)].y;
    float _S633 = 2.0f * _S627.rows[int(1)].x;
    float _S634 = 2.0f * _S627.rows[int(0)].z;
    float _S635 = 2.0f * _S627.rows[int(0)].y;
    float _S636 = 2.0f * - _S627.rows[int(0)].x;
    float _S637 = - _S633 + _S635;
    float _S638 = _S630 + - _S634;
    float _S639 = - _S629 + _S631;
    float _S640 = _S629 + _S631;
    float _S641 = _S630 + _S634;
    float _S642 = _S633 + _S635;
    float _S643 = z_15 * (_S632 + _S636);
    float _S644 = y_18 * (_S628 + _S636);
    float _S645 = x_35 * (_S628 + _S632);
    float _S646 = z_15 * _S637 + y_18 * _S638 + x_35 * _S639;
    float _S647 = _S372 * _S646;
    float _S648 = w_9 * _S637 + y_18 * _S640 + x_35 * _S641 + _S643 + _S643;
    float _S649 = _S372 * _S648;
    float _S650 = w_9 * _S638 + z_15 * _S640 + x_35 * _S642 + _S644 + _S644;
    float _S651 = _S372 * _S650;
    float _S652 = w_9 * _S639 + z_15 * _S641 + y_18 * _S642 + _S645 + _S645;
    float _S653 = _S372 * _S652;
    float _S654 = quat_9.x * _S646 + quat_9.w * _S648 + quat_9.z * _S650 + quat_9.y * _S652;
    DiffPair_float_0 _S655;
    (&_S655)->primal_0 = _S371;
    (&_S655)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S655, _S654);
    float _S656 = quat_9.x * _S655.differential_0;
    float _S657 = quat_9.w * _S655.differential_0;
    float _S658 = quat_9.z * _S655.differential_0;
    float _S659 = quat_9.y * _S655.differential_0;
    float _S660 = _S649 + _S657 + _S657;
    float _S661 = _S651 + _S658 + _S658;
    float _S662 = _S653 + _S659 + _S659;
    float _S663 = _S647 + _S656 + _S656;
    float3  _S664 = _S450;
    *&((&_S664)->z) = _S626.differential_0.rows[int(2)].z;
    *&((&_S664)->y) = _S626.differential_0.rows[int(1)].y;
    *&((&_S664)->x) = _S626.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S665;
    (&_S665)->primal_0 = scale_8;
    (&_S665)->differential_0 = _S450;
    s_bwd_prop_exp_1(&_S665, _S664);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S666 = _S665;
    float3  _S667 = _S450;
    *&((&_S667)->z) = _S614;
    *&((&_S667)->y) = _S615;
    *&((&_S667)->x) = _S616;
    float3  _S668 = _S521.differential_0 + _S667;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S669;
    (&_S669)->primal_0 = R_10;
    (&_S669)->differential_0 = _S512;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S670;
    (&_S670)->primal_0 = mean_6;
    (&_S670)->differential_0 = _S450;
    s_bwd_prop_mul_1(&_S669, &_S670, _S668);
    float3  _S671 = _S668 + _S516.differential_0;
    Matrix<float, 3, 3>  _S672 = _S619 + _S620.differential_0 + _S669.differential_0;
    float4  _S673 = make_float4 (0.0f);
    *&((&_S673)->w) = _S660;
    *&((&_S673)->z) = _S661;
    *&((&_S673)->y) = _S662;
    *&((&_S673)->x) = _S663;
    float4  _S674 = _S673;
    float3  _S675 = _S670.differential_0 + _S510;
    *v_mean_0 = _S675;
    *v_quat_0 = _S674;
    *v_scale_0 = _S666.differential_0;
    *v_in_opacity_0 = _S556;
    (*v_sh_coeffs_0)[int(0)] = _S534;
    (*v_sh_coeffs_0)[int(1)] = _S535;
    (*v_sh_coeffs_0)[int(2)] = _S536;
    (*v_sh_coeffs_0)[int(3)] = _S537;
    (*v_sh_coeffs_0)[int(4)] = _S538;
    (*v_sh_coeffs_0)[int(5)] = _S539;
    (*v_sh_coeffs_0)[int(6)] = _S540;
    (*v_sh_coeffs_0)[int(7)] = _S541;
    (*v_sh_coeffs_0)[int(8)] = _S542;
    (*v_sh_coeffs_0)[int(9)] = _S543;
    (*v_sh_coeffs_0)[int(10)] = _S544;
    (*v_sh_coeffs_0)[int(11)] = _S545;
    (*v_sh_coeffs_0)[int(12)] = _S546;
    (*v_sh_coeffs_0)[int(13)] = _S547;
    (*v_sh_coeffs_0)[int(14)] = _S548;
    (*v_sh_coeffs_0)[int(15)] = _S549;
    *v_R_1 = _S672;
    *v_t_1 = _S671;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S676;
    DiffPair_float_0 _S677;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S678;
    DiffPair_float_0 _S679;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S680;
    DiffPair_float_0 _S681;
    DiffPair_float_0 _S682;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S683;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S684;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S685;
};

inline __device__ CameraDistortion_0 s_primal_ctx_CameraDistortion_x24init_0(float4  dpradial_coeffs_0, float2  dptangential_coeffs_0, float2  dpthin_prism_coeffs_0)
{
    CameraDistortion_0 _S686 = { dpradial_coeffs_0, dptangential_coeffs_0, dpthin_prism_coeffs_0 };
    return _S686;
}

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S687, float _S688)
{
    return s_primal_ctx_atan2_0(_S687, _S688);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S689;
    DiffPair_float_0 _S690;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S691 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S689 = _S691;
    _s_diff_ctx_2->_S690 = _S691;
    (&_s_diff_ctx_2->_S689)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S689)->differential_0 = 0.0f;
    (&_s_diff_ctx_2->_S690)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S690)->differential_0 = 0.0f;
    DiffPair_float_0 _S692 = *dpdpy_0;
    _s_diff_ctx_2->_S689 = *dpdpy_0;
    DiffPair_float_0 _S693 = *dpdpx_0;
    _s_diff_ctx_2->_S690 = *dpdpx_0;
    float _S694 = _S693.primal_0 * _S693.primal_0 + _S692.primal_0 * _S692.primal_0;
    float _S695 = - _S692.primal_0 / _S694 * dpdOut_0;
    float _S696 = _S693.primal_0 / _S694 * dpdOut_0;
    dpdpy_0->primal_0 = _S692.primal_0;
    dpdpy_0->differential_0 = _S696;
    dpdpx_0->primal_0 = _S693.primal_0;
    dpdpx_0->differential_0 = _S695;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S697, DiffPair_float_0 * _S698, float _S699, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S700 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S676 = _S700;
    _s_diff_ctx_3->_S677 = _S700;
    (&_s_diff_ctx_3->_S676)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S676)->differential_0 = 0.0f;
    (&_s_diff_ctx_3->_S677)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S677)->differential_0 = 0.0f;
    DiffPair_float_0 _S701 = *_S697;
    _s_diff_ctx_3->_S676 = *_S697;
    DiffPair_float_0 _S702 = *_S698;
    _s_diff_ctx_3->_S677 = *_S698;
    DiffPair_float_0 _S703 = _S701;
    DiffPair_float_0 _S704 = _S702;
    s_bwd_prop_d_atan2_Intermediates_0 _S705;
    (&_S705)->_S689 = _S700;
    (&_S705)->_S690 = _S700;
    s_primal_ctx_d_atan2_0(&_S703, &_S704, _S699, &_S705);
    *_S697 = _S703;
    *_S698 = _S704;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S706;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S707;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S708;
    DiffPair_float_0 _S709;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S710;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S711;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S712 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S711 = _S712;
    (&_s_diff_ctx_4->_S711)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S711)->differential_0 = 0.0f;
    DiffPair_float_0 _S713 = *dpdpx_1;
    _s_diff_ctx_4->_S711 = *dpdpx_1;
    float _S714 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S713.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S713.primal_0;
    dpdpx_1->differential_0 = _S714;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S715, float _S716, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S717 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S707 = _S717;
    (&_s_diff_ctx_5->_S707)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S707)->differential_0 = 0.0f;
    DiffPair_float_0 _S718 = *_S715;
    _s_diff_ctx_5->_S707 = *_S715;
    DiffPair_float_0 _S719 = _S718;
    s_bwd_prop_d_sqrt_Intermediates_0 _S720;
    (&_S720)->_S711 = _S717;
    s_primal_ctx_d_sqrt_0(&_S719, _S716, &_S720);
    *_S715 = _S719;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_6)
{
    float2  _S721 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S722 = { _S721, _S721 };
    DiffPair_float_0 _S723 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S724 = { _S723 };
    _s_diff_ctx_6->_S708 = _S722;
    _s_diff_ctx_6->_S709 = _S723;
    _s_diff_ctx_6->_S710 = _S724;
    (&_s_diff_ctx_6->_S708)->primal_0 = _S721;
    (&_s_diff_ctx_6->_S708)->differential_0 = _S721;
    (&_s_diff_ctx_6->_S709)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S709)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S725 = *dpdpx_2;
    _s_diff_ctx_6->_S708 = *dpdpx_2;
    float _S726 = _S725.primal_0.x;
    float _S727 = _S725.primal_0.y;
    DiffPair_float_0 _S728;
    (&_S728)->primal_0 = _S726 * _S726 + _S727 * _S727;
    (&_S728)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S728, dp_s_dOut_0, &_s_diff_ctx_6->_S710);
    _s_diff_ctx_6->_S709 = _S728;
    float _S729 = _S725.primal_0.y * _S728.differential_0;
    float _S730 = _S729 + _S729;
    float _S731 = _S725.primal_0.x * _S728.differential_0;
    float _S732 = _S731 + _S731;
    float2  _S733 = _S721;
    *&((&_S733)->y) = _S730;
    *&((&_S733)->x) = _S732;
    dpdpx_2->primal_0 = _S725.primal_0;
    dpdpx_2->differential_0 = _S733;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S734, float _S735, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_7)
{
    float2  _S736 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S737 = { _S736, _S736 };
    _s_diff_ctx_7->_S706 = _S737;
    (&_s_diff_ctx_7->_S706)->primal_0 = _S736;
    (&_s_diff_ctx_7->_S706)->differential_0 = _S736;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S738 = *_S734;
    _s_diff_ctx_7->_S706 = *_S734;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S739 = _S738;
    DiffPair_float_0 _S740 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S741 = { _S740 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S742;
    (&_S742)->_S708 = _S737;
    (&_S742)->_S709 = _S740;
    (&_S742)->_S710 = _S741;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S739, _S735, &_S742);
    *_S734 = _S739;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, float4  dpradial_coeffs_1, float2  dptangential_coeffs_1, float2  dpthin_prism_coeffs_1, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_float_0 _S743 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S744 = { _S743, _S743 };
    _s_diff_ctx_8->_S678 = _S743;
    _s_diff_ctx_8->_S679 = _S743;
    _s_diff_ctx_8->_S680 = _S744;
    _s_diff_ctx_8->_S681 = _S743;
    _s_diff_ctx_8->_S682 = _S743;
    _s_diff_ctx_8->_S683 = _S744;
    (&_s_diff_ctx_8->_S678)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S678)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S679)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S679)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S681)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S681)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S682)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S682)->differential_0 = 0.0f;
    float2  _S745 = make_float2 (0.0f);
    CameraDistortion_0 _S746 = s_primal_ctx_CameraDistortion_x24init_0(dpradial_coeffs_1, dptangential_coeffs_1, dpthin_prism_coeffs_1);
    float2  _S747 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S748 = length_0(_S747);
    float _S749 = dpmean3d_0.z;
    float _S750 = s_primal_ctx_atan2_0(_S748, _S749);
    float k_1;
    if(_S750 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S750 * _S750 / 3.0f) / _S749;
    }
    else
    {
        k_1 = _S750 / _S748;
    }
    float2  _S751 = _S747 * make_float2 (k_1);
    float k1_1 = _S746.radial_coeffs_0.x;
    float k2_1 = _S746.radial_coeffs_0.y;
    float k3_1 = _S746.radial_coeffs_0.z;
    float k4_2 = _S746.radial_coeffs_0.w;
    float p1_2 = _S746.tangential_coeffs_0.x;
    float p2_2 = _S746.tangential_coeffs_0.y;
    float sx1_2 = _S746.thin_prism_coeffs_0.x;
    float sy1_2 = _S746.thin_prism_coeffs_0.y;
    float u_4 = _S751.x;
    float v_4 = _S751.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float _S752 = 2.0f * p1_2;
    float _S753 = 2.0f * p2_2;
    float2  _S754 = _S751 * make_float2 (1.0f + r2_4 * (k1_1 + r2_4 * (k2_1 + r2_4 * (k3_1 + r2_4 * k4_2)))) + make_float2 (_S752 * u_4 * v_4 + p2_2 * (r2_4 + 2.0f * u_4 * u_4) + sx1_2 * r2_4, _S753 * u_4 * v_4 + p1_2 * (r2_4 + 2.0f * v_4 * v_4) + sy1_2 * r2_4);
    float2  _S755 = make_float2 (dpfx_0 * _S754.x + dpcx_0, dpfy_0 * _S754.y + dpcy_0);
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float _S756 = s_primal_ctx_s_primal_ctx_atan2_0(_S748, _S749);
    bool _S757 = _S756 < 0.00100000004749745f;
    float _S758;
    float _S759;
    float _S760;
    if(_S757)
    {
        float _S761 = 1.0f - _S756 * _S756 / 3.0f;
        float _S762 = _S749 * _S749;
        k_1 = _S761 / _S749;
        _S758 = 0.0f;
        _S759 = _S762;
        _S760 = _S761;
    }
    else
    {
        float _S763 = _S748 * _S748;
        k_1 = _S756 / _S748;
        _S758 = _S763;
        _S759 = 0.0f;
        _S760 = 0.0f;
    }
    float2  _S764 = make_float2 (k_1);
    float2  _S765 = _S747 * make_float2 (k_1);
    float u_5 = _S765.x;
    float v_5 = _S765.y;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float _S766 = k3_1 + r2_5 * k4_2;
    float _S767 = k2_1 + r2_5 * _S766;
    float _S768 = k1_1 + r2_5 * _S767;
    float2  _S769 = make_float2 (dpfx_0, 0.0f);
    float2  _S770 = _S765 * _S769;
    float _S771 = p2_2 * dpfx_0;
    float _S772 = _S770.x + _S770.y;
    float _S773 = r2_5 * _S772;
    float _S774 = r2_5 * _S773;
    float _S775 = sx1_2 * dpfx_0 + _S771 + _S768 * _S772 + _S767 * _S773 + _S766 * _S774 + k4_2 * (r2_5 * _S774);
    float _S776 = v_5 * _S775;
    float _S777 = u_5 * _S775;
    float2  _S778 = make_float2 (1.0f + r2_5 * _S768) * _S769 + make_float2 (2.0f * u_5 * _S771 + 2.0f * (u_5 * _S771) + _S752 * (v_5 * dpfx_0) + _S777 + _S777, _S752 * u_5 * dpfx_0 + _S776 + _S776);
    float2  _S779 = _S747 * _S778;
    float2  _S780 = _S764 * _S778;
    float _S781 = _S779.x + _S779.y;
    if(_S757)
    {
        float _S782 = _S781 / _S759;
        float _S783 = _S760 * - _S782;
        float _S784 = _S756 * (0.3333333432674408f * - (_S749 * _S782));
        k_1 = _S784 + _S784;
        _S758 = _S783;
        _S759 = 0.0f;
    }
    else
    {
        float _S785 = _S781 / _S758;
        float _S786 = _S756 * - _S785;
        k_1 = _S748 * _S785;
        _S758 = 0.0f;
        _S759 = _S786;
    }
    DiffPair_float_0 _S787;
    (&_S787)->primal_0 = _S748;
    (&_S787)->differential_0 = 0.0f;
    DiffPair_float_0 _S788;
    (&_S788)->primal_0 = _S749;
    (&_S788)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S787, &_S788, k_1, &_s_diff_ctx_8->_S680);
    _s_diff_ctx_8->_S678 = _S787;
    _s_diff_ctx_8->_S679 = _S788;
    float _S789 = _S788.differential_0 + _S758;
    float _S790 = _S787.differential_0 + _S759;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S791;
    (&_S791)->primal_0 = _S747;
    (&_S791)->differential_0 = _S745;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S792 = { _S745, _S745 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S793;
    (&_S793)->_S706 = _S792;
    s_primal_ctx_s_bwd_length_impl_0(&_S791, _S790, &_S793);
    float2  _S794 = _S791.differential_0 + _S780;
    float3  _S795 = make_float3 (_S794.x, _S794.y, _S789);
    Matrix<float, 2, 3>  _S796 = J_10;
    _S796[int(0)] = _S795;
    if(_S757)
    {
        float _S797 = 1.0f - _S756 * _S756 / 3.0f;
        float _S798 = _S749 * _S749;
        k_1 = _S797 / _S749;
        _S758 = 0.0f;
        _S759 = _S798;
        _S760 = _S797;
    }
    else
    {
        float _S799 = _S748 * _S748;
        k_1 = _S756 / _S748;
        _S758 = _S799;
        _S759 = 0.0f;
        _S760 = 0.0f;
    }
    float2  _S800 = make_float2 (k_1);
    float2  _S801 = _S747 * make_float2 (k_1);
    float u_6 = _S801.x;
    float v_6 = _S801.y;
    float r2_6 = u_6 * u_6 + v_6 * v_6;
    float _S802 = k3_1 + r2_6 * k4_2;
    float _S803 = k2_1 + r2_6 * _S802;
    float _S804 = k1_1 + r2_6 * _S803;
    float2  _S805 = make_float2 (0.0f, dpfy_0);
    float2  _S806 = _S801 * _S805;
    float _S807 = p1_2 * dpfy_0;
    float _S808 = _S806.x + _S806.y;
    float _S809 = r2_6 * _S808;
    float _S810 = r2_6 * _S809;
    float _S811 = sy1_2 * dpfy_0 + _S807 + _S804 * _S808 + _S803 * _S809 + _S802 * _S810 + k4_2 * (r2_6 * _S810);
    float _S812 = v_6 * _S811;
    float _S813 = u_6 * _S811;
    float2  _S814 = make_float2 (1.0f + r2_6 * _S804) * _S805 + make_float2 (_S753 * (v_6 * dpfy_0) + _S813 + _S813, 2.0f * v_6 * _S807 + 2.0f * (v_6 * _S807) + _S753 * u_6 * dpfy_0 + _S812 + _S812);
    float2  _S815 = _S747 * _S814;
    float2  _S816 = _S800 * _S814;
    float _S817 = _S815.x + _S815.y;
    if(_S757)
    {
        float _S818 = _S817 / _S759;
        float _S819 = _S760 * - _S818;
        float _S820 = _S756 * (0.3333333432674408f * - (_S749 * _S818));
        k_1 = _S820 + _S820;
        _S758 = _S819;
        _S759 = 0.0f;
    }
    else
    {
        float _S821 = _S817 / _S758;
        float _S822 = _S756 * - _S821;
        k_1 = _S748 * _S821;
        _S758 = 0.0f;
        _S759 = _S822;
    }
    DiffPair_float_0 _S823;
    (&_S823)->primal_0 = _S748;
    (&_S823)->differential_0 = 0.0f;
    DiffPair_float_0 _S824;
    (&_S824)->primal_0 = _S749;
    (&_S824)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S823, &_S824, k_1, &_s_diff_ctx_8->_S683);
    _s_diff_ctx_8->_S681 = _S823;
    _s_diff_ctx_8->_S682 = _S824;
    float _S825 = _S824.differential_0 + _S758;
    float _S826 = _S823.differential_0 + _S759;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S827;
    (&_S827)->primal_0 = _S747;
    (&_S827)->differential_0 = _S745;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S828;
    (&_S828)->_S706 = _S792;
    s_primal_ctx_s_bwd_length_impl_0(&_S827, _S826, &_S828);
    float2  _S829 = _S827.differential_0 + _S816;
    float3  _S830 = make_float3 (_S829.x, _S829.y, _S825);
    _S796[int(1)] = _S830;
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S796, dpcov3d_0), transpose_1(_S796));
    *dpmean2d_0 = _S755;
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

inline __device__ void s_bwd_prop_d_sqrt_0(DiffPair_1 * dpdpx_3, DiffPair_float_0 * dpdOut_2, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_9)
{
    DiffPair_1 _S831 = *dpdpx_3;
    float _S832 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_9->_S711)->primal_0);
    float _S833 = s_primal_ctx_sqrt_0(_S832);
    float _S834 = 0.5f / _S833 * (*dpdpx_3).differential_0.differential_0;
    float _S835 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S833 * _S833));
    DiffPair_float_0 _S836;
    (&_S836)->primal_0 = _S832;
    (&_S836)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S836, _S835);
    DiffPair_float_0 _S837;
    (&_S837)->primal_0 = 1.00000001168609742e-07f;
    (&_S837)->differential_0 = 0.0f;
    DiffPair_float_0 _S838;
    (&_S838)->primal_0 = (&_s_diff_ctx_9->_S711)->primal_0;
    (&_S838)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S837, &_S838, _S836.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S838.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S834;
    dpdpx_3->primal_0 = _S831.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S839, DiffPair_float_0 * _S840, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_1 _S841 = *_S839;
    DiffPair_float_0 _S842 = _s_diff_ctx_10->_S707;
    DiffPair_float_0 _S843 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S844;
    (&_S844)->_S711 = _S843;
    s_primal_ctx_d_sqrt_0(&_S842, (*_S840).primal_0, &_S844);
    DiffPair_float_0 _S845 = { (*_S839).differential_0.primal_0, (*_S839).differential_0.differential_0 };
    DiffPair_1 _S846;
    (&_S846)->primal_0 = _s_diff_ctx_10->_S707;
    (&_S846)->differential_0 = _S845;
    DiffPair_float_0 _S847;
    (&_S847)->primal_0 = (*_S840).primal_0;
    (&_S847)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S848 = _S844;
    s_bwd_prop_d_sqrt_0(&_S846, &_S847, &_S848);
    DiffPair_float_0 _S849 = { _S846.differential_0.primal_0, _S846.differential_0.differential_0 };
    _S840->primal_0 = (*_S840).primal_0;
    _S840->differential_0 = _S847.differential_0;
    _S839->primal_0 = _S841.primal_0;
    _S839->differential_0 = _S849;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S850, float _s_dOut_3)
{
    DiffPair_float_0 _S851;
    (&_S851)->primal_0 = (*_S850).primal_0;
    (&_S851)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S851, _s_dOut_3);
    _S850->primal_0 = (*_S850).primal_0;
    _S850->differential_0 = _S851.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_0 _S852 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_11->_S708)->primal_0)->x) * *&((&(&_s_diff_ctx_11->_S708)->primal_0)->x) + *&((&(&_s_diff_ctx_11->_S708)->primal_0)->y) * *&((&(&_s_diff_ctx_11->_S708)->primal_0)->y);
    DiffPair_float_0 _S853 = { len_0, 0.0f };
    float2  _S854 = make_float2 (0.0f);
    float _S855 = (*dpdpx_5).differential_0.differential_0.x;
    float _S856 = _S855 + _S855;
    float _S857 = (&_s_diff_ctx_11->_S709)->differential_0 * _S856;
    float _S858 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S859 = (&_s_diff_ctx_11->_S709)->differential_0 * _S858;
    DiffPair_float_0 _S860 = { 0.0f, *&((&(&_s_diff_ctx_11->_S708)->primal_0)->x) * _S856 + *&((&(&_s_diff_ctx_11->_S708)->primal_0)->y) * _S858 };
    DiffPair_1 _S861;
    (&_S861)->primal_0 = _S853;
    (&_S861)->differential_0 = _S860;
    DiffPair_float_0 _S862;
    (&_S862)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S862)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S861, &_S862, &_s_diff_ctx_11->_S710);
    DiffPair_float_0 _S863;
    (&_S863)->primal_0 = len_0;
    (&_S863)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S863, 0.0f);
    float _S864 = _S861.differential_0.primal_0 + _S863.differential_0;
    float _S865 = *&((&(&_s_diff_ctx_11->_S708)->primal_0)->y) * _S864;
    float _S866 = _S859 + _S865 + _S865;
    float _S867 = *&((&(&_s_diff_ctx_11->_S708)->primal_0)->x) * _S864;
    float _S868 = _S857 + _S867 + _S867;
    float2  _S869 = _S854;
    *&((&_S869)->y) = _S866;
    *&((&_S869)->x) = _S868;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S852.differential_0.primal_0 + _S869, _S854 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S862.differential_0;
    dpdpx_5->primal_0 = _S852.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S870 = (*dpdpx_7).primal_0.x;
    float _S871 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S872;
    (&_S872)->primal_0 = _S870 * _S870 + _S871 * _S871;
    (&_S872)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S872, _s_dOut_4);
    float _S873 = (*dpdpx_7).primal_0.y * _S872.differential_0;
    float _S874 = _S873 + _S873;
    float _S875 = (*dpdpx_7).primal_0.x * _S872.differential_0;
    float _S876 = _S875 + _S875;
    float2  _S877 = make_float2 (0.0f);
    *&((&_S877)->y) = _S874;
    *&((&_S877)->x) = _S876;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S877;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S878, DiffPair_float_0 * _S879, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_0 _S880 = *_S878;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S881 = _s_diff_ctx_12->_S706;
    float2  _S882 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S883 = { _S882, _S882 };
    DiffPair_float_0 _S884 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S885 = { _S884 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S886;
    (&_S886)->_S708 = _S883;
    (&_S886)->_S709 = _S884;
    (&_S886)->_S710 = _S885;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S881, (*_S879).primal_0, &_S886);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S887 = { (*_S878).differential_0.primal_0, (*_S878).differential_0.differential_0 };
    DiffPair_0 _S888;
    (&_S888)->primal_0 = _s_diff_ctx_12->_S706;
    (&_S888)->differential_0 = _S887;
    DiffPair_float_0 _S889;
    (&_S889)->primal_0 = (*_S879).primal_0;
    (&_S889)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S890 = _S886;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S888, &_S889, &_S890);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S891;
    (&_S891)->primal_0 = (&_s_diff_ctx_12->_S706)->primal_0;
    (&_S891)->differential_0 = _S882;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S891, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S892 = { _S888.differential_0.primal_0 + _S891.differential_0, _S888.differential_0.differential_0 };
    _S879->primal_0 = (*_S879).primal_0;
    _S879->differential_0 = _S889.differential_0;
    _S878->primal_0 = _S880.primal_0;
    _S878->differential_0 = _S892;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_1 _S893 = *dpdpy_1;
    DiffPair_1 _S894 = *dpdpx_8;
    float _S895 = - (&_s_diff_ctx_13->_S689)->primal_0;
    float _S896 = (&_s_diff_ctx_13->_S690)->primal_0 * (&_s_diff_ctx_13->_S690)->primal_0 + (&_s_diff_ctx_13->_S689)->primal_0 * (&_s_diff_ctx_13->_S689)->primal_0;
    float _S897 = _S896 * _S896;
    float _S898 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S897;
    float _S899 = (&_s_diff_ctx_13->_S690)->primal_0 * - _S898;
    float _S900 = (&_s_diff_ctx_13->_S689)->primal_0 * _S899;
    float _S901 = (&_s_diff_ctx_13->_S690)->primal_0 * _S899;
    float _S902 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S897;
    float _S903 = _S895 * - _S902;
    float _S904 = (&_s_diff_ctx_13->_S689)->primal_0 * _S903;
    float _S905 = (&_s_diff_ctx_13->_S690)->primal_0 * _S903;
    DiffPair_float_0 dpdpx_9 = { _S905 + _S905 + ((*dpdpx_8).differential_0.primal_0 + (_S901 + _S901 + _S896 * _S898)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S900 + _S900 + (*dpdpy_1).differential_0.primal_0 + _S904 + _S904 + - (_S896 * _S902), 0.0f };
    float _S906 = (&_s_diff_ctx_13->_S690)->primal_0 / _S896 * (*dpdpy_1).differential_0.differential_0 + _S895 / _S896 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S906;
    dpdpy_1->primal_0 = _S893.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S894.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S907, DiffPair_1 * _S908, DiffPair_float_0 * _S909, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_1 _S910 = *_S907;
    DiffPair_1 _S911 = *_S908;
    DiffPair_float_0 _S912 = _s_diff_ctx_14->_S676;
    DiffPair_float_0 _S913 = _s_diff_ctx_14->_S677;
    DiffPair_float_0 _S914 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S915;
    (&_S915)->_S689 = _S914;
    (&_S915)->_S690 = _S914;
    s_primal_ctx_d_atan2_0(&_S912, &_S913, (*_S909).primal_0, &_S915);
    DiffPair_float_0 _S916 = { (*_S908).differential_0.primal_0, (*_S908).differential_0.differential_0 };
    DiffPair_float_0 _S917 = { (*_S907).differential_0.primal_0, (*_S907).differential_0.differential_0 };
    DiffPair_1 _S918;
    (&_S918)->primal_0 = _s_diff_ctx_14->_S676;
    (&_S918)->differential_0 = _S917;
    DiffPair_1 _S919;
    (&_S919)->primal_0 = _s_diff_ctx_14->_S677;
    (&_S919)->differential_0 = _S916;
    DiffPair_float_0 _S920;
    (&_S920)->primal_0 = (*_S909).primal_0;
    (&_S920)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S921 = _S915;
    s_bwd_prop_d_atan2_0(&_S918, &_S919, &_S920, &_S921);
    DiffPair_float_0 _S922 = { _S919.differential_0.primal_0, _S919.differential_0.differential_0 };
    DiffPair_float_0 _S923 = { _S918.differential_0.primal_0, _S918.differential_0.differential_0 };
    _S909->primal_0 = (*_S909).primal_0;
    _S909->differential_0 = _S920.differential_0;
    _S907->primal_0 = _S910.primal_0;
    _S907->differential_0 = _S923;
    _S908->primal_0 = _S911.primal_0;
    _S908->differential_0 = _S922;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S924, DiffPair_float_0 * _S925, float _s_dOut_5)
{
    DiffPair_float_0 _S926;
    (&_S926)->primal_0 = (*_S924).primal_0;
    (&_S926)->differential_0 = 0.0f;
    DiffPair_float_0 _S927;
    (&_S927)->primal_0 = (*_S925).primal_0;
    (&_S927)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S926, &_S927, _s_dOut_5);
    _S925->primal_0 = (*_S925).primal_0;
    _S925->differential_0 = _S927.differential_0;
    _S924->primal_0 = (*_S924).primal_0;
    _S924->differential_0 = _S926.differential_0;
    return;
}

inline __device__ void s_bwd_prop_CameraDistortion_x24init_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_2, CameraDistortion_0 * _s_dOut_6)
{
    float2  _S928 = _s_dOut_6->thin_prism_coeffs_0;
    dpthin_prism_coeffs_2->primal_0 = (*dpthin_prism_coeffs_2).primal_0;
    dpthin_prism_coeffs_2->differential_0 = _S928;
    float2  _S929 = _s_dOut_6->tangential_coeffs_0;
    dptangential_coeffs_2->primal_0 = (*dptangential_coeffs_2).primal_0;
    dptangential_coeffs_2->differential_0 = _S929;
    float4  _S930 = _s_dOut_6->radial_coeffs_0;
    dpradial_coeffs_2->primal_0 = (*dpradial_coeffs_2).primal_0;
    dpradial_coeffs_2->differential_0 = _S930;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_3, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S931 = *dpcov3d_1;
    DiffPair_float_0 _S932 = *dpfx_1;
    DiffPair_float_0 _S933 = *dpfy_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S934 = *dpradial_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S935 = *dptangential_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S936 = *dpthin_prism_coeffs_3;
    float2  _S937 = make_float2 (0.0f);
    CameraDistortion_0 _S938 = s_primal_ctx_CameraDistortion_x24init_0((*dpradial_coeffs_3).primal_0, (*dptangential_coeffs_3).primal_0, (*dpthin_prism_coeffs_3).primal_0);
    float2  _S939 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S940 = length_0(_S939);
    float _S941 = (*dpmean3d_1).primal_0.z;
    float _S942 = s_primal_ctx_atan2_0(_S940, _S941);
    bool _S943 = _S942 < 0.00100000004749745f;
    float k_2;
    float _S944;
    float _S945;
    float _S946;
    if(_S943)
    {
        float _S947 = 1.0f - _S942 * _S942 / 3.0f;
        float _S948 = _S941 * _S941;
        k_2 = _S947 / _S941;
        _S944 = 0.0f;
        _S945 = _S948;
        _S946 = _S947;
    }
    else
    {
        float _S949 = _S940 * _S940;
        k_2 = _S942 / _S940;
        _S944 = _S949;
        _S945 = 0.0f;
        _S946 = 0.0f;
    }
    float2  _S950 = make_float2 (k_2);
    float2  _S951 = _S939 * make_float2 (k_2);
    float k1_2 = _S938.radial_coeffs_0.x;
    float k2_2 = _S938.radial_coeffs_0.y;
    float k3_2 = _S938.radial_coeffs_0.z;
    float k4_3 = _S938.radial_coeffs_0.w;
    float p1_3 = _S938.tangential_coeffs_0.x;
    float p2_3 = _S938.tangential_coeffs_0.y;
    float sx1_3 = _S938.thin_prism_coeffs_0.x;
    float sy1_3 = _S938.thin_prism_coeffs_0.y;
    float u_7 = _S951.x;
    float v_7 = _S951.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S952 = k3_2 + r2_7 * k4_3;
    float _S953 = k2_2 + r2_7 * _S952;
    float _S954 = k1_2 + r2_7 * _S953;
    float radial_1 = 1.0f + r2_7 * _S954;
    float2  _S955 = make_float2 (radial_1);
    float _S956 = 2.0f * p1_3;
    float _S957 = _S956 * u_7;
    float _S958 = 2.0f * u_7;
    float _S959 = r2_7 + _S958 * u_7;
    float _S960 = 2.0f * p2_3;
    float _S961 = _S960 * u_7;
    float _S962 = 2.0f * v_7;
    float _S963 = r2_7 + _S962 * v_7;
    float2  _S964 = _S951 * make_float2 (radial_1) + make_float2 (_S957 * v_7 + p2_3 * _S959 + sx1_3 * r2_7, _S961 * v_7 + p1_3 * _S963 + sy1_3 * r2_7);
    float _S965 = _S964.x;
    float _S966 = _S964.y;
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (0.0f);
    float _S967 = s_primal_ctx_s_primal_ctx_atan2_0(_S940, _S941);
    bool _S968 = _S967 < 0.00100000004749745f;
    float _S969;
    float _S970;
    float _S971;
    if(_S968)
    {
        float _S972 = 1.0f - _S967 * _S967 / 3.0f;
        float _S973 = _S941 * _S941;
        k_2 = _S972 / _S941;
        _S969 = 0.0f;
        _S970 = _S973;
        _S971 = _S972;
    }
    else
    {
        float _S974 = _S940 * _S940;
        k_2 = _S967 / _S940;
        _S969 = _S974;
        _S970 = 0.0f;
        _S971 = 0.0f;
    }
    float2  _S975 = make_float2 (k_2);
    float2  _S976 = _S939 * make_float2 (k_2);
    float u_8 = _S976.x;
    float v_8 = _S976.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S977 = k3_2 + r2_8 * k4_3;
    float _S978 = k2_2 + r2_8 * _S977;
    float _S979 = k1_2 + r2_8 * _S978;
    float2  _S980 = make_float2 (1.0f + r2_8 * _S979);
    float _S981 = _S956 * u_8;
    float _S982 = 2.0f * u_8;
    float2  _S983 = make_float2 (_S932.primal_0, 0.0f);
    float2  _S984 = _S976 * _S983;
    float _S985 = p2_3 * _S932.primal_0;
    float _S986 = v_8 * _S932.primal_0;
    float _S987 = _S984.x + _S984.y;
    float _S988 = r2_8 * _S987;
    float _S989 = r2_8 * _S988;
    float _S990 = r2_8 * _S989;
    float _S991 = sx1_3 * _S932.primal_0 + _S985 + _S979 * _S987 + _S978 * _S988 + _S977 * _S989 + k4_3 * _S990;
    float _S992 = v_8 * _S991;
    float _S993 = u_8 * _S991;
    float2  _S994 = _S980 * _S983 + make_float2 (_S982 * _S985 + 2.0f * (u_8 * _S985) + _S956 * _S986 + _S993 + _S993, _S981 * _S932.primal_0 + _S992 + _S992);
    float2  _S995 = _S939 * _S994;
    float2  _S996 = _S975 * _S994;
    float _S997 = _S995.x + _S995.y;
    float k_3;
    float _S998;
    float _S999;
    float _S1000;
    float _S1001;
    float _S1002;
    float _S1003;
    float _S1004;
    float _S1005;
    if(_S968)
    {
        float _S1006 = _S997 / _S970;
        float _S1007 = _S970 * _S970;
        float _S1008 = - _S1006;
        float _S1009 = _S971 * _S1008;
        float _S1010 = 0.3333333432674408f * - (_S941 * _S1006);
        float _S1011 = _S967 * _S1010;
        k_2 = _S1011 + _S1011;
        k_3 = _S1009;
        _S998 = 0.0f;
        _S999 = 0.0f;
        _S1000 = 0.0f;
        _S1001 = 0.0f;
        _S1002 = _S1010;
        _S1003 = _S1006;
        _S1004 = _S1008;
        _S1005 = _S1007;
    }
    else
    {
        float _S1012 = _S997 / _S969;
        float _S1013 = _S969 * _S969;
        float _S1014 = - _S1012;
        float _S1015 = _S967 * _S1014;
        k_2 = _S940 * _S1012;
        k_3 = 0.0f;
        _S998 = _S1015;
        _S999 = _S1012;
        _S1000 = _S1014;
        _S1001 = _S1013;
        _S1002 = 0.0f;
        _S1003 = 0.0f;
        _S1004 = 0.0f;
        _S1005 = 0.0f;
    }
    DiffPair_float_0 _S1016 = { _S940, 0.0f };
    DiffPair_float_0 _S1017 = { _S941, 0.0f };
    float _S1018 = (&_s_diff_ctx_15->_S679)->differential_0 + k_3;
    float _S1019 = (&_s_diff_ctx_15->_S678)->differential_0 + _S998;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1020 = { _S939, _S937 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1021;
    (&_S1021)->primal_0 = _S939;
    (&_S1021)->differential_0 = _S937;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1022 = { _S937, _S937 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1023;
    (&_S1023)->_S706 = _S1022;
    s_primal_ctx_s_bwd_length_impl_0(&_S1021, _S1019, &_S1023);
    float2  _S1024 = _S1021.differential_0 + _S996;
    float3  _S1025 = make_float3 (_S1024.x, _S1024.y, _S1018);
    Matrix<float, 2, 3>  _S1026 = J_11;
    _S1026[int(0)] = _S1025;
    float _S1027;
    float _S1028;
    if(_S968)
    {
        float _S1029 = 1.0f - _S967 * _S967 / 3.0f;
        float _S1030 = _S941 * _S941;
        k_3 = _S1029 / _S941;
        _S998 = 0.0f;
        _S1027 = _S1030;
        _S1028 = _S1029;
    }
    else
    {
        float _S1031 = _S940 * _S940;
        k_3 = _S967 / _S940;
        _S998 = _S1031;
        _S1027 = 0.0f;
        _S1028 = 0.0f;
    }
    float2  _S1032 = make_float2 (k_3);
    float2  _S1033 = _S939 * make_float2 (k_3);
    float u_9 = _S1033.x;
    float v_9 = _S1033.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S1034 = k3_2 + r2_9 * k4_3;
    float _S1035 = k2_2 + r2_9 * _S1034;
    float _S1036 = k1_2 + r2_9 * _S1035;
    float2  _S1037 = make_float2 (1.0f + r2_9 * _S1036);
    float _S1038 = _S960 * u_9;
    float _S1039 = 2.0f * v_9;
    float2  _S1040 = make_float2 (0.0f, _S933.primal_0);
    float2  _S1041 = _S1033 * _S1040;
    float _S1042 = p1_3 * _S933.primal_0;
    float _S1043 = v_9 * _S933.primal_0;
    float _S1044 = _S1041.x + _S1041.y;
    float _S1045 = r2_9 * _S1044;
    float _S1046 = r2_9 * _S1045;
    float _S1047 = r2_9 * _S1046;
    float _S1048 = sy1_3 * _S933.primal_0 + _S1042 + _S1036 * _S1044 + _S1035 * _S1045 + _S1034 * _S1046 + k4_3 * _S1047;
    float _S1049 = v_9 * _S1048;
    float _S1050 = u_9 * _S1048;
    float2  _S1051 = _S1037 * _S1040 + make_float2 (_S960 * _S1043 + _S1050 + _S1050, _S1039 * _S1042 + 2.0f * (v_9 * _S1042) + _S1038 * _S933.primal_0 + _S1049 + _S1049);
    float2  _S1052 = _S939 * _S1051;
    float2  _S1053 = _S1032 * _S1051;
    float _S1054 = _S1052.x + _S1052.y;
    float _S1055;
    float _S1056;
    float _S1057;
    float _S1058;
    float _S1059;
    float _S1060;
    float _S1061;
    float _S1062;
    float _S1063;
    if(_S968)
    {
        float _S1064 = _S1054 / _S1027;
        float _S1065 = _S1027 * _S1027;
        float _S1066 = - _S1064;
        float _S1067 = _S1028 * _S1066;
        float _S1068 = 0.3333333432674408f * - (_S941 * _S1064);
        float _S1069 = _S967 * _S1068;
        k_3 = _S1069 + _S1069;
        _S1055 = _S1067;
        _S1056 = 0.0f;
        _S1057 = 0.0f;
        _S1058 = 0.0f;
        _S1059 = 0.0f;
        _S1060 = _S1068;
        _S1061 = _S1064;
        _S1062 = _S1066;
        _S1063 = _S1065;
    }
    else
    {
        float _S1070 = _S1054 / _S998;
        float _S1071 = _S998 * _S998;
        float _S1072 = - _S1070;
        float _S1073 = _S967 * _S1072;
        k_3 = _S940 * _S1070;
        _S1055 = 0.0f;
        _S1056 = _S1073;
        _S1057 = _S1070;
        _S1058 = _S1072;
        _S1059 = _S1071;
        _S1060 = 0.0f;
        _S1061 = 0.0f;
        _S1062 = 0.0f;
        _S1063 = 0.0f;
    }
    float _S1074 = (&_s_diff_ctx_15->_S682)->differential_0 + _S1055;
    float _S1075 = (&_s_diff_ctx_15->_S681)->differential_0 + _S1056;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1076;
    (&_S1076)->primal_0 = _S939;
    (&_S1076)->differential_0 = _S937;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1077;
    (&_S1077)->_S706 = _S1022;
    s_primal_ctx_s_bwd_length_impl_0(&_S1076, _S1075, &_S1077);
    float2  _S1078 = _S1076.differential_0 + _S1053;
    float3  _S1079 = make_float3 (_S1078.x, _S1078.y, _S1074);
    _S1026[int(1)] = _S1079;
    Matrix<float, 3, 2>  _S1080 = transpose_1(_S1026);
    CameraDistortion_0 _S1081 = CameraDistortion_x24_syn_dzero_0();
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1082;
    (&_S1082)->primal_0 = s_primal_ctx_mul_3(_S1026, _S931.primal_0);
    (&_S1082)->differential_0 = J_11;
    Matrix<float, 3, 2>  _S1083 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1084;
    (&_S1084)->primal_0 = _S1080;
    (&_S1084)->differential_0 = _S1083;
    s_bwd_prop_mul_2(&_S1082, &_S1084, dpcov2d_1);
    Matrix<float, 2, 3>  _S1085 = transpose_2(_S1084.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1086;
    (&_S1086)->primal_0 = _S1026;
    (&_S1086)->differential_0 = J_11;
    Matrix<float, 3, 3>  _S1087 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1088;
    (&_S1088)->primal_0 = _S931.primal_0;
    (&_S1088)->differential_0 = _S1087;
    s_bwd_prop_mul_3(&_S1086, &_S1088, _S1082.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1089 = _S1088;
    Matrix<float, 2, 3>  _S1090 = _S1085 + _S1086.differential_0;
    float2  _S1091 = _S937;
    *&((&_S1091)->y) = _S1090.rows[int(1)].y;
    *&((&_S1091)->x) = _S1090.rows[int(1)].x;
    float2  _S1092 = _S1091;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1093 = { _S937, _S1091 };
    DiffPair_0 _S1094;
    (&_S1094)->primal_0 = _S1020;
    (&_S1094)->differential_0 = _S1093;
    DiffPair_float_0 _S1095;
    (&_S1095)->primal_0 = _S1075;
    (&_S1095)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1096 = _S1077;
    s_bwd_prop_s_bwd_length_impl_0(&_S1094, &_S1095, &_S1096);
    DiffPair_0 _S1097 = _S1094;
    DiffPair_float_0 _S1098 = _S1095;
    DiffPair_float_0 _S1099 = { 0.0f, _S1090.rows[int(1)].z };
    DiffPair_float_0 _S1100 = { 0.0f, _S1095.differential_0 };
    DiffPair_1 _S1101;
    (&_S1101)->primal_0 = _S1016;
    (&_S1101)->differential_0 = _S1100;
    DiffPair_1 _S1102;
    (&_S1102)->primal_0 = _S1017;
    (&_S1102)->differential_0 = _S1099;
    DiffPair_float_0 _S1103;
    (&_S1103)->primal_0 = k_3;
    (&_S1103)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1101, &_S1102, &_S1103, &_s_diff_ctx_15->_S683);
    DiffPair_1 _S1104 = _S1101;
    DiffPair_1 _S1105 = _S1102;
    DiffPair_float_0 _S1106 = _S1103;
    if(_S968)
    {
        float _S1107 = _S1106.differential_0 + _S1106.differential_0;
        float _S1108 = _S1060 * _S1107;
        float _S1109 = - (0.3333333432674408f * (_S967 * _S1107));
        float _S1110 = _S1062 * _S1090.rows[int(1)].z;
        float _S1111 = (_S941 * _S1109 + - (_S1028 * _S1090.rows[int(1)].z)) / _S1063;
        float _S1112 = _S1054 * - _S1111;
        float _S1113 = _S1061 * _S1109 + _S1105.differential_0.primal_0;
        k_3 = _S1027 * _S1111;
        _S1055 = 0.0f;
        _S1056 = _S1112;
        _S1057 = _S1110;
        _S1058 = _S1104.differential_0.primal_0;
        _S1059 = _S1108;
        _S1060 = _S1113;
    }
    else
    {
        float _S1114 = _S1058 * _S1098.differential_0;
        float _S1115 = (_S940 * _S1106.differential_0 + - (_S967 * _S1098.differential_0)) / _S1059;
        float _S1116 = _S1054 * - _S1115;
        float _S1117 = _S1057 * _S1106.differential_0 + _S1104.differential_0.primal_0;
        k_3 = _S998 * _S1115;
        _S1055 = _S1116;
        _S1056 = 0.0f;
        _S1057 = 0.0f;
        _S1058 = _S1117;
        _S1059 = _S1114;
        _S1060 = _S1105.differential_0.primal_0;
    }
    float2  _S1118 = _S1032 * _S1092;
    float2  _S1119 = _S1051 * _S1092;
    float2  _S1120 = _S937;
    *&((&_S1120)->y) = k_3;
    *&((&_S1120)->x) = k_3;
    float2  _S1121 = _S1051 * _S1120;
    float2  _S1122 = _S1118 + _S939 * _S1120;
    float _S1123 = _S1122.x;
    float _S1124 = _S1123 + _S1123;
    float _S1125 = _S1048 * _S1124;
    float _S1126 = _S1122.y + _S1122.y;
    float _S1127 = _S1048 * _S1126;
    float _S1128 = u_9 * _S1124 + v_9 * _S1126;
    float _S1129 = k4_3 * _S1128;
    float _S1130 = _S1047 * _S1128;
    float _S1131 = _S1046 * _S1128;
    float _S1132 = _S1046 * _S1129;
    float _S1133 = _S1045 * _S1128;
    float _S1134 = _S1034 * _S1128 + r2_9 * _S1129;
    float _S1135 = _S1045 * _S1134;
    float _S1136 = _S1044 * _S1128;
    float _S1137 = _S1035 * _S1128 + r2_9 * _S1134;
    float _S1138 = _S1044 * _S1137;
    float _S1139 = _S1036 * _S1128 + r2_9 * _S1137;
    float _S1140 = _S960 * _S1122.x;
    float _S1141 = _S1043 * _S1122.x;
    float _S1142 = v_9 * _S1140;
    float _S1143 = _S933.primal_0 * _S1140;
    float _S1144 = _S1038 * _S1122.y;
    float _S1145 = _S933.primal_0 * _S1122.y;
    float _S1146 = 2.0f * _S1122.y;
    float _S1147 = _S1042 * _S1146;
    float _S1148 = _S1042 * _S1122.y;
    float _S1149 = _S1128 + v_9 * _S1146 + _S1039 * _S1122.y;
    float _S1150 = p1_3 * _S1149;
    float _S1151 = _S933.primal_0 * _S1149;
    float _S1152 = sy1_3 * _S1128;
    float _S1153 = _S933.primal_0 * _S1128;
    float2  _S1154 = _S1037 * _S1122;
    float2  _S1155 = _S1040 * _S1122;
    float2  _S1156 = _S937;
    *&((&_S1156)->y) = _S1139;
    *&((&_S1156)->x) = _S1139;
    float _S1157 = _S1155.x + _S1155.y;
    float _S1158 = _S1136 + r2_9 * _S1157;
    float _S1159 = _S1133 + r2_9 * _S1158;
    float _S1160 = _S1131 + r2_9 * _S1159;
    float _S1161 = _S1132 + _S1135 + _S1138 + _S1036 * _S1157 + _S1035 * _S1158 + _S1034 * _S1159 + k4_3 * _S1160;
    float _S1162 = v_9 * _S1161;
    float _S1163 = u_9 * _S1161;
    float2  _S1164 = _S1121 + _S1097.differential_0.primal_0;
    float2  _S1165 = _S1040 * _S1156 + make_float2 (_S1125 + _S960 * _S1145 + _S1163 + _S1163, _S1127 + _S1143 + _S1147 + 2.0f * _S1148 + _S1162 + _S1162);
    float _S1166 = _S1142 + _S1144 + _S1150 + _S1152 + (_S1154 + _S1033 * _S1156).y;
    float _S1167 = _S1141 + u_9 * _S1145;
    float _S1168 = _S1130 + r2_9 * _S1160;
    float2  _S1169 = _S939 * _S1165;
    float _S1170 = _S1119.x + _S1119.y + _S1169.x + _S1169.y;
    float2  _S1171 = _S1032 * _S1165 + _S1164;
    if(_S968)
    {
        float _S1172 = _S941 * _S1056;
        float _S1173 = _S1170 / _S1027;
        float _S1174 = _S967 * (0.3333333432674408f * - (_S1057 + _S941 * _S1173));
        float _S1175 = _S1172 + _S1172 + _S1028 * - _S1173 + _S1060;
        k_3 = _S1174 + _S1174 + _S1059;
        _S998 = _S1175;
        _S1027 = _S1058;
    }
    else
    {
        float _S1176 = _S940 * _S1055;
        float _S1177 = _S1170 / _S998;
        float _S1178 = _S1176 + _S1176 + _S967 * - _S1177 + _S1058;
        k_3 = _S940 * _S1177 + _S1059;
        _S998 = _S1060;
        _S1027 = _S1178;
    }
    DiffPair_float_0 _S1179;
    (&_S1179)->primal_0 = _S940;
    (&_S1179)->differential_0 = 0.0f;
    DiffPair_float_0 _S1180;
    (&_S1180)->primal_0 = _S941;
    (&_S1180)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1179, &_S1180, k_3);
    float _S1181 = _S1180.differential_0 + _S998;
    float _S1182 = _S1179.differential_0 + _S1027;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1183;
    (&_S1183)->primal_0 = _S939;
    (&_S1183)->differential_0 = _S937;
    s_bwd_length_impl_1(&_S1183, _S1182);
    float2  _S1184 = _S1183.differential_0 + _S1171;
    float3  _S1185 = make_float3 (_S1184.x, _S1184.y, _S1181);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1186;
    (&_S1186)->primal_0 = _S939;
    (&_S1186)->differential_0 = _S937;
    s_bwd_length_impl_1(&_S1186, 0.0f);
    float3  _S1187 = _S1185 + make_float3 (_S1186.differential_0.x, _S1186.differential_0.y, 0.0f);
    float2  _S1188 = _S937;
    *&((&_S1188)->y) = _S1090.rows[int(0)].y;
    *&((&_S1188)->x) = _S1090.rows[int(0)].x;
    float2  _S1189 = _S1188;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1190 = { _S937, _S1188 };
    DiffPair_0 _S1191;
    (&_S1191)->primal_0 = _S1020;
    (&_S1191)->differential_0 = _S1190;
    DiffPair_float_0 _S1192;
    (&_S1192)->primal_0 = _S1019;
    (&_S1192)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1193 = _S1023;
    s_bwd_prop_s_bwd_length_impl_0(&_S1191, &_S1192, &_S1193);
    DiffPair_0 _S1194 = _S1191;
    DiffPair_float_0 _S1195 = _S1192;
    DiffPair_float_0 _S1196 = { 0.0f, _S1090.rows[int(0)].z };
    DiffPair_float_0 _S1197 = { 0.0f, _S1192.differential_0 };
    DiffPair_1 _S1198;
    (&_S1198)->primal_0 = _S1016;
    (&_S1198)->differential_0 = _S1197;
    DiffPair_1 _S1199;
    (&_S1199)->primal_0 = _S1017;
    (&_S1199)->differential_0 = _S1196;
    DiffPair_float_0 _S1200;
    (&_S1200)->primal_0 = k_2;
    (&_S1200)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1198, &_S1199, &_S1200, &_s_diff_ctx_15->_S680);
    DiffPair_1 _S1201 = _S1198;
    DiffPair_1 _S1202 = _S1199;
    DiffPair_float_0 _S1203 = _S1200;
    if(_S968)
    {
        float _S1204 = _S1203.differential_0 + _S1203.differential_0;
        float _S1205 = _S1002 * _S1204;
        float _S1206 = - (0.3333333432674408f * (_S967 * _S1204));
        float _S1207 = _S1004 * _S1090.rows[int(0)].z;
        float _S1208 = (_S941 * _S1206 + - (_S971 * _S1090.rows[int(0)].z)) / _S1005;
        float _S1209 = _S997 * - _S1208;
        float _S1210 = _S1003 * _S1206 + _S1202.differential_0.primal_0;
        k_2 = _S970 * _S1208;
        k_3 = 0.0f;
        _S998 = _S1209;
        _S999 = _S1207;
        _S1000 = _S1201.differential_0.primal_0;
        _S1001 = _S1205;
        _S1002 = _S1210;
    }
    else
    {
        float _S1211 = _S1000 * _S1195.differential_0;
        float _S1212 = (_S940 * _S1203.differential_0 + - (_S967 * _S1195.differential_0)) / _S1001;
        float _S1213 = _S997 * - _S1212;
        float _S1214 = _S999 * _S1203.differential_0 + _S1201.differential_0.primal_0;
        k_2 = _S969 * _S1212;
        k_3 = _S1213;
        _S998 = 0.0f;
        _S999 = 0.0f;
        _S1000 = _S1214;
        _S1001 = _S1211;
        _S1002 = _S1202.differential_0.primal_0;
    }
    float2  _S1215 = _S975 * _S1189;
    float2  _S1216 = _S994 * _S1189;
    float2  _S1217 = _S937;
    *&((&_S1217)->y) = k_2;
    *&((&_S1217)->x) = k_2;
    float2  _S1218 = _S994 * _S1217;
    float2  _S1219 = _S1215 + _S939 * _S1217;
    float _S1220 = _S1219.x;
    float _S1221 = _S1220 + _S1220;
    float _S1222 = _S991 * _S1221;
    float _S1223 = _S1219.y + _S1219.y;
    float _S1224 = _S991 * _S1223;
    float _S1225 = u_8 * _S1221 + v_8 * _S1223;
    float _S1226 = k4_3 * _S1225;
    float _S1227 = _S990 * _S1225;
    float _S1228 = _S989 * _S1225;
    float _S1229 = _S989 * _S1226;
    float _S1230 = _S988 * _S1225;
    float _S1231 = _S977 * _S1225 + r2_8 * _S1226;
    float _S1232 = _S988 * _S1231;
    float _S1233 = _S987 * _S1225;
    float _S1234 = _S978 * _S1225 + r2_8 * _S1231;
    float _S1235 = _S987 * _S1234;
    float _S1236 = _S979 * _S1225 + r2_8 * _S1234;
    float _S1237 = _S956 * _S1219.x;
    float _S1238 = _S986 * _S1219.x;
    float _S1239 = v_8 * _S1237;
    float _S1240 = _S932.primal_0 * _S1237;
    float _S1241 = _S981 * _S1219.y;
    float _S1242 = _S932.primal_0 * _S1219.y;
    float _S1243 = 2.0f * _S1219.x;
    float _S1244 = _S985 * _S1243;
    float _S1245 = _S985 * _S1219.x;
    float _S1246 = _S1225 + u_8 * _S1243 + _S982 * _S1219.x;
    float _S1247 = p2_3 * _S1246;
    float _S1248 = _S932.primal_0 * _S1246;
    float _S1249 = sx1_3 * _S1225;
    float _S1250 = _S932.primal_0 * _S1225;
    float2  _S1251 = _S980 * _S1219;
    float2  _S1252 = _S983 * _S1219;
    float2  _S1253 = _S937;
    *&((&_S1253)->y) = _S1236;
    *&((&_S1253)->x) = _S1236;
    float _S1254 = _S1252.x + _S1252.y;
    float _S1255 = _S1233 + r2_8 * _S1254;
    float _S1256 = _S1230 + r2_8 * _S1255;
    float _S1257 = _S1228 + r2_8 * _S1256;
    float _S1258 = _S1229 + _S1232 + _S1235 + _S979 * _S1254 + _S978 * _S1255 + _S977 * _S1256 + k4_3 * _S1257;
    float _S1259 = v_8 * _S1258;
    float _S1260 = u_8 * _S1258;
    float2  _S1261 = _S1218 + _S1194.differential_0.primal_0;
    float _S1262 = _S1238 + u_8 * _S1242;
    float _S1263 = _S1255 + _S1158;
    float _S1264 = _S1257 + _S1160;
    float _S1265 = _S1239 + _S1241 + _S1247 + _S1249 + (_S1251 + _S976 * _S1253).x;
    float2  _S1266 = _S983 * _S1253 + make_float2 (_S1222 + _S1244 + 2.0f * _S1245 + _S956 * _S1242 + _S1260 + _S1260, _S1224 + _S1240 + _S1259 + _S1259);
    float _S1267 = _S1256 + _S1159;
    float _S1268 = _S1227 + r2_8 * _S1257 + _S1168;
    float2  _S1269 = _S939 * _S1266;
    float _S1270 = _S1216.x + _S1216.y + _S1269.x + _S1269.y;
    float2  _S1271 = _S975 * _S1266 + _S1261;
    if(_S968)
    {
        float _S1272 = _S941 * _S998;
        float _S1273 = _S1270 / _S970;
        float _S1274 = _S967 * (0.3333333432674408f * - (_S999 + _S941 * _S1273));
        float _S1275 = _S1272 + _S1272 + _S971 * - _S1273 + _S1002;
        k_2 = _S1274 + _S1274 + _S1001;
        _S969 = _S1275;
        _S970 = _S1000;
    }
    else
    {
        float _S1276 = _S940 * k_3;
        float _S1277 = _S1270 / _S969;
        float _S1278 = _S1276 + _S1276 + _S967 * - _S1277 + _S1000;
        k_2 = _S940 * _S1277 + _S1001;
        _S969 = _S1002;
        _S970 = _S1278;
    }
    DiffPair_float_0 _S1279;
    (&_S1279)->primal_0 = _S940;
    (&_S1279)->differential_0 = 0.0f;
    DiffPair_float_0 _S1280;
    (&_S1280)->primal_0 = _S941;
    (&_S1280)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1279, &_S1280, k_2);
    float _S1281 = _S1280.differential_0 + _S969;
    float _S1282 = _S1279.differential_0 + _S970;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1283;
    (&_S1283)->primal_0 = _S939;
    (&_S1283)->differential_0 = _S937;
    s_bwd_length_impl_1(&_S1283, _S1282);
    float2  _S1284 = _S1283.differential_0 + _S1271;
    float3  _S1285 = make_float3 (_S1284.x, _S1284.y, _S1281);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1286;
    (&_S1286)->primal_0 = _S939;
    (&_S1286)->differential_0 = _S937;
    s_bwd_length_impl_1(&_S1286, 0.0f);
    float _S1287 = _S933.primal_0 * dpmean2d_1.y;
    float _S1288 = _S932.primal_0 * dpmean2d_1.x;
    float2  _S1289 = make_float2 (_S1288, _S1287);
    float2  _S1290 = _S951 * _S1289;
    float _S1291 = p1_3 * _S1287;
    float _S1292 = v_7 * _S1287;
    float _S1293 = p2_3 * _S1288;
    float _S1294 = v_7 * _S1288;
    float _S1295 = _S1290.x + _S1290.y;
    float _S1296 = r2_7 * _S1295;
    float _S1297 = r2_7 * _S1296;
    float _S1298 = r2_7 * _S1297;
    float _S1299 = sy1_3 * _S1287 + _S1291 + sx1_3 * _S1288 + _S1293 + _S954 * _S1295 + _S953 * _S1296 + _S952 * _S1297 + k4_3 * _S1298;
    float _S1300 = v_7 * _S1299;
    float _S1301 = u_7 * _S1299;
    float2  _S1302 = make_float2 (r2_7 * _S1288 + _S1250, r2_7 * _S1287 + _S1153);
    float2  _S1303 = make_float2 (_S963 * _S1287 + 2.0f * (u_7 * _S1294 + _S1262) + _S1151, 2.0f * (u_7 * _S1292 + _S1167) + _S959 * _S1288 + _S1248);
    float4  _S1304 = make_float4 (_S1296 + _S1263, _S1297 + _S1267, _S1298 + _S1264, r2_7 * _S1298 + _S1268);
    float3  _S1305 = _S1285 + make_float3 (_S1286.differential_0.x, _S1286.differential_0.y, 0.0f) + _S1187;
    float _S1306 = _S965 * dpmean2d_1.x + _S1265;
    float _S1307 = _S966 * dpmean2d_1.y + _S1166;
    float2  _S1308 = _S955 * _S1289 + make_float2 (_S960 * _S1292 + _S958 * _S1293 + 2.0f * (u_7 * _S1293) + _S956 * _S1294 + _S1301 + _S1301, _S962 * _S1291 + 2.0f * (v_7 * _S1291) + _S961 * _S1287 + _S957 * _S1288 + _S1300 + _S1300);
    CameraDistortion_0 _S1309 = _S1081;
    (&_S1309)->thin_prism_coeffs_0 = _S1302;
    (&_S1309)->tangential_coeffs_0 = _S1303;
    (&_S1309)->radial_coeffs_0 = _S1304;
    CameraDistortion_0 _S1310 = _S1081;
    CameraDistortion_0 _S1311 = _S1309;
    CameraDistortion_0 _S1312 = CameraDistortion_x24_syn_dadd_0(&_S1310, &_S1311);
    float2  _S1313 = _S939 * _S1308;
    float2  _S1314 = _S950 * _S1308;
    float _S1315 = _S1313.x + _S1313.y;
    if(_S943)
    {
        float _S1316 = _S1315 / _S945;
        float _S1317 = _S946 * - _S1316;
        float _S1318 = _S942 * (0.3333333432674408f * - (_S941 * _S1316));
        k_2 = _S1318 + _S1318;
        _S944 = _S1317;
        _S945 = 0.0f;
    }
    else
    {
        float _S1319 = _S1315 / _S944;
        float _S1320 = _S942 * - _S1319;
        k_2 = _S940 * _S1319;
        _S944 = 0.0f;
        _S945 = _S1320;
    }
    DiffPair_float_0 _S1321;
    (&_S1321)->primal_0 = _S940;
    (&_S1321)->differential_0 = 0.0f;
    DiffPair_float_0 _S1322;
    (&_S1322)->primal_0 = _S941;
    (&_S1322)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1321, &_S1322, k_2);
    float _S1323 = _S1322.differential_0 + _S944;
    float _S1324 = _S1321.differential_0 + _S945;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1325;
    (&_S1325)->primal_0 = _S939;
    (&_S1325)->differential_0 = _S937;
    s_bwd_length_impl_1(&_S1325, _S1324);
    float2  _S1326 = _S1325.differential_0 + _S1314;
    float4  _S1327 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1328;
    (&_S1328)->primal_0 = _S934.primal_0;
    (&_S1328)->differential_0 = _S1327;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1329;
    (&_S1329)->primal_0 = _S935.primal_0;
    (&_S1329)->differential_0 = _S937;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1330;
    (&_S1330)->primal_0 = _S936.primal_0;
    (&_S1330)->differential_0 = _S937;
    CameraDistortion_0 _S1331 = _S1312;
    s_bwd_prop_CameraDistortion_x24init_0(&_S1328, &_S1329, &_S1330, &_S1331);
    dpthin_prism_coeffs_3->primal_0 = (*dpthin_prism_coeffs_3).primal_0;
    dpthin_prism_coeffs_3->differential_0 = _S1330.differential_0;
    dptangential_coeffs_3->primal_0 = (*dptangential_coeffs_3).primal_0;
    dptangential_coeffs_3->differential_0 = _S1329.differential_0;
    dpradial_coeffs_3->primal_0 = (*dpradial_coeffs_3).primal_0;
    dpradial_coeffs_3->differential_0 = _S1328.differential_0;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = dpmean2d_1.y;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = dpmean2d_1.x;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S1307;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S1306;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1089.differential_0;
    float3  _S1332 = _S1305 + make_float3 (_S1326.x, _S1326.y, _S1323);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1332;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_11, float fy_11, float cx_11, float cy_11, float4  radial_coeffs_14, float2  tangential_coeffs_14, float2  thin_prism_coeffs_14, uint image_width_7, uint image_height_7, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1333 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1334 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1335 = { _S1334, _S1334 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1336 = { _S1334, _S1334, _S1335, _S1334, _S1334, _S1335 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1337;
    (&_S1337)->_S684 = _S1333;
    (&_S1337)->_S685 = _S1336;
    float3  mean_c_7 = s_primal_ctx_mul_1(R_11, mean_7) + t_10;
    float3  _S1338 = s_primal_ctx_exp_0(scale_9);
    float _S1339 = quat_10.y;
    float _S1340 = _S1339 * _S1339 + quat_10.z * quat_10.z + quat_10.w * quat_10.w + quat_10.x * quat_10.x;
    float _S1341 = s_primal_ctx_rsqrt_0(_S1340);
    float x_37 = quat_10.y * _S1341;
    float y_20 = quat_10.z * _S1341;
    float z_17 = quat_10.w * _S1341;
    float w_10 = quat_10.x * _S1341;
    float x2_10 = x_37 * x_37;
    float y2_10 = y_20 * y_20;
    float z2_17 = z_17 * z_17;
    float xy_10 = x_37 * y_20;
    float xz_10 = x_37 * z_17;
    float yz_10 = y_20 * z_17;
    float wx_10 = w_10 * x_37;
    float wy_10 = w_10 * y_20;
    float wz_10 = w_10 * z_17;
    Matrix<float, 3, 3>  _S1342 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1338.x, 0.0f, 0.0f, 0.0f, _S1338.y, 0.0f, 0.0f, 0.0f, _S1338.z);
    Matrix<float, 3, 3>  _S1343 = s_primal_ctx_mul_2(_S1342, S_1);
    Matrix<float, 3, 3>  _S1344 = transpose_0(_S1343);
    Matrix<float, 3, 3>  _S1345 = s_primal_ctx_mul_2(_S1343, _S1344);
    Matrix<float, 3, 3>  _S1346 = s_primal_ctx_mul_2(R_11, _S1345);
    Matrix<float, 3, 3>  _S1347 = transpose_0(R_11);
    Matrix<float, 3, 3>  _S1348 = s_primal_ctx_mul_2(_S1346, _S1347);
    Matrix<float, 2, 2>  _S1349 = _S1333;
    float2  _S1350 = make_float2 (0.0f);
    float2  _S1351 = _S1350;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_7, _S1348, fx_11, fy_11, cx_11, cy_11, radial_coeffs_14, tangential_coeffs_14, thin_prism_coeffs_14, &_S1349, &_S1351, &(&_S1337)->_S685);
    (&_S1337)->_S684 = _S1349;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1352 = _S1337;
    float _S1353 = _S1337._S684.rows[int(0)].y * _S1337._S684.rows[int(1)].x;
    float det_orig_8 = _S1337._S684.rows[int(0)].x * _S1337._S684.rows[int(1)].y - _S1353;
    float _S1354 = _S1337._S684.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1355 = _S1337._S684;
    *&(((&_S1355)->rows + (int(0)))->x) = _S1354;
    float _S1356 = _S1337._S684.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1355)->rows + (int(1)))->y) = _S1356;
    Matrix<float, 2, 2>  _S1357 = _S1355;
    Matrix<float, 2, 2>  _S1358 = _S1355;
    float det_blur_5 = _S1354 * _S1356 - _S1353;
    float _S1359 = det_orig_8 / det_blur_5;
    float _S1360 = det_blur_5 * det_blur_5;
    float _S1361 = s_primal_ctx_max_0(0.0f, _S1359);
    float _S1362 = s_primal_ctx_sqrt_0(_S1361);
    float invdet_8 = 1.0f / det_blur_5;
    float _S1363 = - _S1337._S684.rows[int(0)].y;
    float _S1364 = - _S1337._S684.rows[int(1)].x;
    float _S1365 = - in_opacity_7;
    float _S1366 = 1.0f + s_primal_ctx_exp_1(_S1365);
    float _S1367 = 1.0f / _S1366;
    float _S1368 = _S1366 * _S1366;
    float _S1369;
    if(antialiased_7)
    {
        _S1369 = _S1367 * _S1362;
    }
    else
    {
        _S1369 = _S1367;
    }
    float _S1370 = _S1369 / 0.00392156885936856f;
    float _S1371 = 2.0f * s_primal_ctx_log_0(_S1370);
    float _S1372 = s_primal_ctx_sqrt_0(_S1371);
    float _S1373 = _S1357.rows[int(0)].x;
    float _S1374 = _S1358.rows[int(1)].y;
    float3  _S1375 = mean_7 - - s_primal_ctx_mul_1(_S1347, t_10);
    float _S1376 = _S1375.x;
    float _S1377 = _S1375.y;
    float _S1378 = _S1375.z;
    float _S1379 = _S1376 * _S1376 + _S1377 * _S1377 + _S1378 * _S1378;
    float _S1380 = s_primal_ctx_sqrt_0(_S1379);
    float x_38 = _S1376 / _S1380;
    float3  _S1381 = make_float3 (x_38);
    float _S1382 = _S1380 * _S1380;
    float y_21 = _S1377 / _S1380;
    float z_18 = _S1378 / _S1380;
    float3  _S1383 = make_float3 (z_18);
    float _S1384 = - y_21;
    float3  _S1385 = make_float3 (_S1384);
    float z2_18 = z_18 * z_18;
    float fTmp0B_7 = -1.09254848957061768f * z_18;
    float fC1_7 = x_38 * x_38 - y_21 * y_21;
    float _S1386 = 2.0f * x_38;
    float fS1_7 = _S1386 * y_21;
    float pSH6_1 = 0.94617468118667603f * z2_18 - 0.31539157032966614f;
    float3  _S1387 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_7 * x_38;
    float3  _S1388 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_7 * y_21;
    float3  _S1389 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_7;
    float3  _S1390 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_7;
    float3  _S1391 = make_float3 (pSH4_1);
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_18;
    float _S1392 = 1.86588168144226074f * z2_18 - 1.11952900886535645f;
    float pSH12_1 = z_18 * _S1392;
    float3  _S1393 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_7 * x_38;
    float3  _S1394 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_7 * y_21;
    float3  _S1395 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_7 * fC1_7;
    float3  _S1396 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_7 * fS1_7;
    float3  _S1397 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_38 * fC1_7 - y_21 * fS1_7);
    float3  _S1398 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_38 * fS1_7 + y_21 * fC1_7);
    float3  _S1399 = make_float3 (pSH9_1);
    float3  _S1400 = make_float3 (0.0f);
    float3  _S1401 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1402;
    (&_S1402)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1384) * (*sh_coeffs_7)[int(1)] + make_float3 (z_18) * (*sh_coeffs_7)[int(2)] - make_float3 (x_38) * (*sh_coeffs_7)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_7)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_7)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_7)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_7)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_7)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_7)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_7)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_7)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_7)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_7)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_7)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f);
    (&_S1402)->differential_0 = _S1401;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1403;
    (&_S1403)->primal_0 = _S1400;
    (&_S1403)->differential_0 = _S1401;
    s_bwd_prop_max_0(&_S1402, &_S1403, v_rgb_1);
    float3  _S1404 = _S1398 * _S1402.differential_0;
    float3  _S1405 = (*sh_coeffs_7)[int(15)] * _S1402.differential_0;
    float3  _S1406 = _S1396 * _S1402.differential_0;
    float3  _S1407 = (*sh_coeffs_7)[int(14)] * _S1402.differential_0;
    float3  _S1408 = _S1394 * _S1402.differential_0;
    float3  _S1409 = (*sh_coeffs_7)[int(13)] * _S1402.differential_0;
    float3  _S1410 = _S1393 * _S1402.differential_0;
    float3  _S1411 = (*sh_coeffs_7)[int(12)] * _S1402.differential_0;
    float3  _S1412 = _S1395 * _S1402.differential_0;
    float3  _S1413 = (*sh_coeffs_7)[int(11)] * _S1402.differential_0;
    float3  _S1414 = _S1397 * _S1402.differential_0;
    float3  _S1415 = (*sh_coeffs_7)[int(10)] * _S1402.differential_0;
    float3  _S1416 = _S1399 * _S1402.differential_0;
    float3  _S1417 = (*sh_coeffs_7)[int(9)] * _S1402.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1417.x + _S1417.y + _S1417.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1405.x + _S1405.y + _S1405.z);
    float _S1418 = _S1415.x + _S1415.y + _S1415.z;
    float _S1419 = _S1407.x + _S1407.y + _S1407.z;
    float _S1420 = _S1413.x + _S1413.y + _S1413.z;
    float _S1421 = _S1409.x + _S1409.y + _S1409.z;
    float _S1422 = _S1411.x + _S1411.y + _S1411.z;
    float _S1423 = - s_diff_fC2_T_1;
    float3  _S1424 = _S1390 * _S1402.differential_0;
    float3  _S1425 = (*sh_coeffs_7)[int(8)] * _S1402.differential_0;
    float3  _S1426 = _S1388 * _S1402.differential_0;
    float3  _S1427 = (*sh_coeffs_7)[int(7)] * _S1402.differential_0;
    float3  _S1428 = _S1387 * _S1402.differential_0;
    float3  _S1429 = (*sh_coeffs_7)[int(6)] * _S1402.differential_0;
    float3  _S1430 = _S1389 * _S1402.differential_0;
    float3  _S1431 = (*sh_coeffs_7)[int(5)] * _S1402.differential_0;
    float3  _S1432 = _S1391 * _S1402.differential_0;
    float3  _S1433 = (*sh_coeffs_7)[int(4)] * _S1402.differential_0;
    float _S1434 = _S1431.x + _S1431.y + _S1431.z;
    float _S1435 = _S1427.x + _S1427.y + _S1427.z;
    float _S1436 = fTmp1B_7 * _S1418 + x_38 * s_diff_fS2_T_1 + y_21 * _S1423 + 0.54627424478530884f * (_S1433.x + _S1433.y + _S1433.z);
    float _S1437 = fTmp1B_7 * _S1419 + y_21 * s_diff_fS2_T_1 + x_38 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1425.x + _S1425.y + _S1425.z);
    float _S1438 = y_21 * - _S1437;
    float _S1439 = x_38 * _S1437;
    float _S1440 = z_18 * (1.86588168144226074f * (z_18 * _S1422) + -2.28522896766662598f * (y_21 * _S1420 + x_38 * _S1421) + 0.94617468118667603f * (_S1429.x + _S1429.y + _S1429.z));
    float3  _S1441 = make_float3 (0.48860251903533936f) * _S1402.differential_0;
    float3  _S1442 = - _S1441;
    float3  _S1443 = _S1381 * _S1442;
    float3  _S1444 = (*sh_coeffs_7)[int(3)] * _S1442;
    float3  _S1445 = _S1383 * _S1441;
    float3  _S1446 = (*sh_coeffs_7)[int(2)] * _S1441;
    float3  _S1447 = _S1385 * _S1441;
    float3  _S1448 = (*sh_coeffs_7)[int(1)] * _S1441;
    float _S1449 = (_S1392 * _S1422 + 1.44530570507049561f * (fS1_7 * _S1418 + fC1_7 * _S1419) + -1.09254848957061768f * (y_21 * _S1434 + x_38 * _S1435) + _S1440 + _S1440 + _S1446.x + _S1446.y + _S1446.z) / _S1382;
    float _S1450 = _S1380 * _S1449;
    float _S1451 = (fTmp0C_7 * _S1420 + fC1_7 * s_diff_fS2_T_1 + fS1_7 * _S1423 + fTmp0B_7 * _S1434 + _S1386 * _S1436 + _S1438 + _S1438 + - (_S1448.x + _S1448.y + _S1448.z)) / _S1382;
    float _S1452 = _S1380 * _S1451;
    float _S1453 = (fTmp0C_7 * _S1421 + fS1_7 * s_diff_fS2_T_1 + fC1_7 * s_diff_fC2_T_1 + fTmp0B_7 * _S1435 + 2.0f * (y_21 * _S1436) + _S1439 + _S1439 + _S1444.x + _S1444.y + _S1444.z) / _S1382;
    float _S1454 = _S1380 * _S1453;
    float _S1455 = _S1378 * - _S1449 + _S1377 * - _S1451 + _S1376 * - _S1453;
    DiffPair_float_0 _S1456;
    (&_S1456)->primal_0 = _S1379;
    (&_S1456)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1456, _S1455);
    float _S1457 = _S1378 * _S1456.differential_0;
    float _S1458 = _S1377 * _S1456.differential_0;
    float _S1459 = _S1376 * _S1456.differential_0;
    float3  _S1460 = make_float3 (0.282094806432724f) * _S1402.differential_0;
    float3  _S1461 = make_float3 (_S1454 + _S1459 + _S1459, _S1452 + _S1458 + _S1458, _S1450 + _S1457 + _S1457);
    float3  _S1462 = - - _S1461;
    Matrix<float, 3, 3>  _S1463 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1464;
    (&_S1464)->primal_0 = _S1347;
    (&_S1464)->differential_0 = _S1463;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1465;
    (&_S1465)->primal_0 = t_10;
    (&_S1465)->differential_0 = _S1401;
    s_bwd_prop_mul_1(&_S1464, &_S1465, _S1462);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1466 = _S1464;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1467 = _S1465;
    float2  _S1468 = _S1350;
    *&((&_S1468)->y) = v_conic_1.z;
    float2  _S1469 = _S1350;
    *&((&_S1469)->y) = v_conic_1.y;
    *&((&_S1469)->x) = v_conic_1.x;
    DiffPair_float_0 _S1470;
    (&_S1470)->primal_0 = _S1374;
    (&_S1470)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1470, 0.0f);
    DiffPair_float_0 _S1471;
    (&_S1471)->primal_0 = _S1373;
    (&_S1471)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1471, 0.0f);
    DiffPair_float_0 _S1472;
    (&_S1472)->primal_0 = 3.32999992370605469f;
    (&_S1472)->differential_0 = 0.0f;
    DiffPair_float_0 _S1473;
    (&_S1473)->primal_0 = _S1372;
    (&_S1473)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1472, &_S1473, 0.0f);
    DiffPair_float_0 _S1474;
    (&_S1474)->primal_0 = _S1371;
    (&_S1474)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1474, _S1473.differential_0);
    float _S1475 = 2.0f * _S1474.differential_0;
    DiffPair_float_0 _S1476;
    (&_S1476)->primal_0 = _S1370;
    (&_S1476)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1476, _S1475);
    float2  _S1477 = make_float2 (_S1471.differential_0, 0.0f);
    float _S1478 = v_opacity_1 + 254.9999847412109375f * _S1476.differential_0;
    FixedArray<float3 , 16>  _S1479;
    _S1479[int(0)] = _S1401;
    _S1479[int(1)] = _S1401;
    _S1479[int(2)] = _S1401;
    _S1479[int(3)] = _S1401;
    _S1479[int(4)] = _S1401;
    _S1479[int(5)] = _S1401;
    _S1479[int(6)] = _S1401;
    _S1479[int(7)] = _S1401;
    _S1479[int(8)] = _S1401;
    _S1479[int(9)] = _S1401;
    _S1479[int(10)] = _S1401;
    _S1479[int(11)] = _S1401;
    _S1479[int(12)] = _S1401;
    _S1479[int(13)] = _S1401;
    _S1479[int(14)] = _S1401;
    _S1479[int(15)] = _S1401;
    _S1479[int(7)] = _S1426;
    _S1479[int(0)] = _S1460;
    _S1479[int(1)] = _S1447;
    _S1479[int(2)] = _S1445;
    _S1479[int(3)] = _S1443;
    _S1479[int(4)] = _S1432;
    _S1479[int(5)] = _S1430;
    _S1479[int(6)] = _S1428;
    _S1479[int(15)] = _S1404;
    _S1479[int(8)] = _S1424;
    _S1479[int(9)] = _S1416;
    _S1479[int(10)] = _S1414;
    _S1479[int(11)] = _S1412;
    _S1479[int(12)] = _S1410;
    _S1479[int(13)] = _S1408;
    _S1479[int(14)] = _S1406;
    float3  _S1480 = _S1479[int(0)];
    float3  _S1481 = _S1479[int(1)];
    float3  _S1482 = _S1479[int(2)];
    float3  _S1483 = _S1479[int(3)];
    float3  _S1484 = _S1479[int(4)];
    float3  _S1485 = _S1479[int(5)];
    float3  _S1486 = _S1479[int(6)];
    float3  _S1487 = _S1479[int(7)];
    float3  _S1488 = _S1479[int(8)];
    float3  _S1489 = _S1479[int(9)];
    float3  _S1490 = _S1479[int(10)];
    float3  _S1491 = _S1479[int(11)];
    float3  _S1492 = _S1479[int(12)];
    float3  _S1493 = _S1479[int(13)];
    float3  _S1494 = _S1479[int(14)];
    float3  _S1495 = _S1479[int(15)];
    Matrix<float, 2, 2>  _S1496 = _S1333;
    _S1496[int(1)] = _S1468;
    _S1496[int(0)] = _S1469;
    Matrix<float, 2, 2>  _S1497 = _S1496;
    float2  _S1498 = make_float2 (0.0f, _S1470.differential_0);
    float _S1499;
    if(antialiased_7)
    {
        float _S1500 = _S1367 * _S1478;
        _S1369 = _S1362 * _S1478;
        _S1499 = _S1500;
    }
    else
    {
        _S1369 = _S1478;
        _S1499 = 0.0f;
    }
    float _S1501 = - (_S1369 / _S1368);
    DiffPair_float_0 _S1502;
    (&_S1502)->primal_0 = _S1365;
    (&_S1502)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1502, _S1501);
    float _S1503 = - _S1502.differential_0;
    float _S1504 = invdet_8 * _S1497.rows[int(1)].y;
    float _S1505 = - (invdet_8 * _S1497.rows[int(1)].x);
    float _S1506 = - (invdet_8 * _S1497.rows[int(0)].y);
    float _S1507 = invdet_8 * _S1497.rows[int(0)].x;
    float _S1508 = - ((_S1354 * _S1497.rows[int(1)].y + _S1364 * _S1497.rows[int(1)].x + _S1363 * _S1497.rows[int(0)].y + _S1356 * _S1497.rows[int(0)].x) / _S1360);
    DiffPair_float_0 _S1509;
    (&_S1509)->primal_0 = _S1361;
    (&_S1509)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1509, _S1499);
    DiffPair_float_0 _S1510;
    (&_S1510)->primal_0 = 0.0f;
    (&_S1510)->differential_0 = 0.0f;
    DiffPair_float_0 _S1511;
    (&_S1511)->primal_0 = _S1359;
    (&_S1511)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1510, &_S1511, _S1509.differential_0);
    float _S1512 = _S1511.differential_0 / _S1360;
    float s_diff_det_orig_T_1 = det_blur_5 * _S1512;
    float _S1513 = _S1508 + det_orig_8 * - _S1512;
    float _S1514 = - _S1513;
    float _S1515 = _S1354 * _S1513;
    float _S1516 = _S1356 * _S1513;
    Matrix<float, 2, 2>  _S1517 = _S1333;
    _S1517[int(1)] = _S1498;
    _S1517[int(0)] = _S1477;
    _S1355 = _S1517;
    *&(((&_S1355)->rows + (int(1)))->y) = 0.0f;
    float _S1518 = _S1507 + _S1515 + _S1517.rows[int(1)].y;
    *&(((&_S1355)->rows + (int(0)))->x) = 0.0f;
    float _S1519 = _S1504 + _S1516 + _S1517.rows[int(0)].x;
    float _S1520 = _S1514 + - s_diff_det_orig_T_1;
    float _S1521 = _S1505 + _S1352._S684.rows[int(0)].y * _S1520;
    float _S1522 = _S1506 + _S1352._S684.rows[int(1)].x * _S1520;
    float _S1523 = _S1352._S684.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S1524 = _S1518 + _S1352._S684.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S1525 = _S1350;
    *&((&_S1525)->x) = _S1521;
    *&((&_S1525)->y) = _S1524;
    float _S1526 = _S1519 + _S1523;
    float2  _S1527 = _S1350;
    *&((&_S1527)->y) = _S1522;
    *&((&_S1527)->x) = _S1526;
    Matrix<float, 2, 2>  _S1528 = _S1333;
    _S1528[int(1)] = _S1525;
    _S1528[int(0)] = _S1527;
    Matrix<float, 2, 2>  _S1529 = _S1355 + _S1528;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1530;
    (&_S1530)->primal_0 = mean_c_7;
    (&_S1530)->differential_0 = _S1401;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1531;
    (&_S1531)->primal_0 = _S1348;
    (&_S1531)->differential_0 = _S1463;
    DiffPair_float_0 _S1532;
    (&_S1532)->primal_0 = fx_11;
    (&_S1532)->differential_0 = 0.0f;
    DiffPair_float_0 _S1533;
    (&_S1533)->primal_0 = fy_11;
    (&_S1533)->differential_0 = 0.0f;
    DiffPair_float_0 _S1534;
    (&_S1534)->primal_0 = cx_11;
    (&_S1534)->differential_0 = 0.0f;
    DiffPair_float_0 _S1535;
    (&_S1535)->primal_0 = cy_11;
    (&_S1535)->differential_0 = 0.0f;
    float4  _S1536 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1537;
    (&_S1537)->primal_0 = radial_coeffs_14;
    (&_S1537)->differential_0 = _S1536;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1538;
    (&_S1538)->primal_0 = tangential_coeffs_14;
    (&_S1538)->differential_0 = _S1350;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1539;
    (&_S1539)->primal_0 = thin_prism_coeffs_14;
    (&_S1539)->differential_0 = _S1350;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1540 = _S1352._S685;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S1530, &_S1531, &_S1532, &_S1533, &_S1534, &_S1535, &_S1537, &_S1538, &_S1539, _S1529, v_mean2d_1, &_S1540);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1541;
    (&_S1541)->primal_0 = _S1346;
    (&_S1541)->differential_0 = _S1463;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1542;
    (&_S1542)->primal_0 = _S1347;
    (&_S1542)->differential_0 = _S1463;
    s_bwd_prop_mul_4(&_S1541, &_S1542, _S1531.differential_0);
    Matrix<float, 3, 3>  _S1543 = transpose_0(_S1542.differential_0 + _S1466.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1544;
    (&_S1544)->primal_0 = R_11;
    (&_S1544)->differential_0 = _S1463;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1545;
    (&_S1545)->primal_0 = _S1345;
    (&_S1545)->differential_0 = _S1463;
    s_bwd_prop_mul_4(&_S1544, &_S1545, _S1541.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1546;
    (&_S1546)->primal_0 = _S1343;
    (&_S1546)->differential_0 = _S1463;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1547;
    (&_S1547)->primal_0 = _S1344;
    (&_S1547)->differential_0 = _S1463;
    s_bwd_prop_mul_4(&_S1546, &_S1547, _S1545.differential_0);
    Matrix<float, 3, 3>  _S1548 = _S1546.differential_0 + transpose_0(_S1547.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1549;
    (&_S1549)->primal_0 = _S1342;
    (&_S1549)->differential_0 = _S1463;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1550;
    (&_S1550)->primal_0 = S_1;
    (&_S1550)->differential_0 = _S1463;
    s_bwd_prop_mul_4(&_S1549, &_S1550, _S1548);
    Matrix<float, 3, 3>  _S1551 = transpose_0(_S1549.differential_0);
    float _S1552 = 2.0f * - _S1551.rows[int(2)].z;
    float _S1553 = 2.0f * _S1551.rows[int(2)].y;
    float _S1554 = 2.0f * _S1551.rows[int(2)].x;
    float _S1555 = 2.0f * _S1551.rows[int(1)].z;
    float _S1556 = 2.0f * - _S1551.rows[int(1)].y;
    float _S1557 = 2.0f * _S1551.rows[int(1)].x;
    float _S1558 = 2.0f * _S1551.rows[int(0)].z;
    float _S1559 = 2.0f * _S1551.rows[int(0)].y;
    float _S1560 = 2.0f * - _S1551.rows[int(0)].x;
    float _S1561 = - _S1557 + _S1559;
    float _S1562 = _S1554 + - _S1558;
    float _S1563 = - _S1553 + _S1555;
    float _S1564 = _S1553 + _S1555;
    float _S1565 = _S1554 + _S1558;
    float _S1566 = _S1557 + _S1559;
    float _S1567 = z_17 * (_S1556 + _S1560);
    float _S1568 = y_20 * (_S1552 + _S1560);
    float _S1569 = x_37 * (_S1552 + _S1556);
    float _S1570 = z_17 * _S1561 + y_20 * _S1562 + x_37 * _S1563;
    float _S1571 = _S1341 * _S1570;
    float _S1572 = w_10 * _S1561 + y_20 * _S1564 + x_37 * _S1565 + _S1567 + _S1567;
    float _S1573 = _S1341 * _S1572;
    float _S1574 = w_10 * _S1562 + z_17 * _S1564 + x_37 * _S1566 + _S1568 + _S1568;
    float _S1575 = _S1341 * _S1574;
    float _S1576 = w_10 * _S1563 + z_17 * _S1565 + y_20 * _S1566 + _S1569 + _S1569;
    float _S1577 = _S1341 * _S1576;
    float _S1578 = quat_10.x * _S1570 + quat_10.w * _S1572 + quat_10.z * _S1574 + quat_10.y * _S1576;
    DiffPair_float_0 _S1579;
    (&_S1579)->primal_0 = _S1340;
    (&_S1579)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1579, _S1578);
    float _S1580 = quat_10.x * _S1579.differential_0;
    float _S1581 = quat_10.w * _S1579.differential_0;
    float _S1582 = quat_10.z * _S1579.differential_0;
    float _S1583 = quat_10.y * _S1579.differential_0;
    float _S1584 = _S1573 + _S1581 + _S1581;
    float _S1585 = _S1575 + _S1582 + _S1582;
    float _S1586 = _S1577 + _S1583 + _S1583;
    float _S1587 = _S1571 + _S1580 + _S1580;
    float3  _S1588 = _S1401;
    *&((&_S1588)->z) = _S1550.differential_0.rows[int(2)].z;
    *&((&_S1588)->y) = _S1550.differential_0.rows[int(1)].y;
    *&((&_S1588)->x) = _S1550.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1589;
    (&_S1589)->primal_0 = scale_9;
    (&_S1589)->differential_0 = _S1401;
    s_bwd_prop_exp_1(&_S1589, _S1588);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1590 = _S1589;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1591;
    (&_S1591)->primal_0 = mean_c_7;
    (&_S1591)->differential_0 = _S1401;
    s_bwd_length_impl_0(&_S1591, v_depth_1);
    float3  _S1592 = _S1530.differential_0 + _S1591.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1593;
    (&_S1593)->primal_0 = R_11;
    (&_S1593)->differential_0 = _S1463;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1594;
    (&_S1594)->primal_0 = mean_7;
    (&_S1594)->differential_0 = _S1401;
    s_bwd_prop_mul_1(&_S1593, &_S1594, _S1592);
    float3  _S1595 = _S1592 + _S1467.differential_0;
    Matrix<float, 3, 3>  _S1596 = _S1543 + _S1544.differential_0 + _S1593.differential_0;
    float4  _S1597 = _S1536;
    *&((&_S1597)->w) = _S1584;
    *&((&_S1597)->z) = _S1585;
    *&((&_S1597)->y) = _S1586;
    *&((&_S1597)->x) = _S1587;
    float4  _S1598 = _S1597;
    float3  _S1599 = _S1594.differential_0 + _S1461;
    *v_mean_1 = _S1599;
    *v_quat_1 = _S1598;
    *v_scale_1 = _S1590.differential_0;
    *v_in_opacity_1 = _S1503;
    (*v_sh_coeffs_1)[int(0)] = _S1480;
    (*v_sh_coeffs_1)[int(1)] = _S1481;
    (*v_sh_coeffs_1)[int(2)] = _S1482;
    (*v_sh_coeffs_1)[int(3)] = _S1483;
    (*v_sh_coeffs_1)[int(4)] = _S1484;
    (*v_sh_coeffs_1)[int(5)] = _S1485;
    (*v_sh_coeffs_1)[int(6)] = _S1486;
    (*v_sh_coeffs_1)[int(7)] = _S1487;
    (*v_sh_coeffs_1)[int(8)] = _S1488;
    (*v_sh_coeffs_1)[int(9)] = _S1489;
    (*v_sh_coeffs_1)[int(10)] = _S1490;
    (*v_sh_coeffs_1)[int(11)] = _S1491;
    (*v_sh_coeffs_1)[int(12)] = _S1492;
    (*v_sh_coeffs_1)[int(13)] = _S1493;
    (*v_sh_coeffs_1)[int(14)] = _S1494;
    (*v_sh_coeffs_1)[int(15)] = _S1495;
    *v_R_2 = _S1596;
    *v_t_2 = _S1595;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_12, float fy_12, float cx_12, float cy_12, float4  radial_coeffs_15, float2  tangential_coeffs_15, float2  thin_prism_coeffs_15, uint image_width_8, uint image_height_8, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_8 = s_primal_ctx_mul_1(R_12, mean_8) + t_11;
    float3  _S1600 = s_primal_ctx_exp_0(scale_10);
    float _S1601 = quat_11.y;
    float _S1602 = _S1601 * _S1601 + quat_11.z * quat_11.z + quat_11.w * quat_11.w + quat_11.x * quat_11.x;
    float _S1603 = s_primal_ctx_rsqrt_0(_S1602);
    float x_39 = quat_11.y * _S1603;
    float y_22 = quat_11.z * _S1603;
    float z_19 = quat_11.w * _S1603;
    float w_11 = quat_11.x * _S1603;
    float x2_11 = x_39 * x_39;
    float y2_11 = y_22 * y_22;
    float z2_19 = z_19 * z_19;
    float xy_11 = x_39 * y_22;
    float xz_11 = x_39 * z_19;
    float yz_11 = y_22 * z_19;
    float wx_11 = w_11 * x_39;
    float wy_11 = w_11 * y_22;
    float wz_11 = w_11 * z_19;
    Matrix<float, 3, 3>  _S1604 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S1600.x, 0.0f, 0.0f, 0.0f, _S1600.y, 0.0f, 0.0f, 0.0f, _S1600.z);
    Matrix<float, 3, 3>  _S1605 = s_primal_ctx_mul_2(_S1604, S_2);
    Matrix<float, 3, 3>  _S1606 = transpose_0(_S1605);
    Matrix<float, 3, 3>  _S1607 = s_primal_ctx_mul_2(_S1605, _S1606);
    Matrix<float, 3, 3>  _S1608 = s_primal_ctx_mul_2(R_12, _S1607);
    Matrix<float, 3, 3>  _S1609 = transpose_0(R_12);
    Matrix<float, 3, 3>  _S1610 = s_primal_ctx_mul_2(_S1608, _S1609);
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (fx_12, 0.0f, 0.0f, 0.0f, fy_12, 0.0f);
    Matrix<float, 2, 3>  _S1611 = s_primal_ctx_mul_3(J_12, _S1610);
    Matrix<float, 3, 2>  _S1612 = transpose_1(J_12);
    Matrix<float, 2, 2>  _S1613 = s_primal_ctx_mul_4(_S1611, _S1612);
    float _S1614 = _S1613.rows[int(0)].y * _S1613.rows[int(1)].x;
    float det_orig_9 = _S1613.rows[int(0)].x * _S1613.rows[int(1)].y - _S1614;
    float _S1615 = _S1613.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1616 = _S1613;
    *&(((&_S1616)->rows + (int(0)))->x) = _S1615;
    float _S1617 = _S1613.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1616)->rows + (int(1)))->y) = _S1617;
    Matrix<float, 2, 2>  _S1618 = _S1616;
    Matrix<float, 2, 2>  _S1619 = _S1616;
    float det_blur_6 = _S1615 * _S1617 - _S1614;
    float _S1620 = det_orig_9 / det_blur_6;
    float _S1621 = det_blur_6 * det_blur_6;
    float _S1622 = s_primal_ctx_max_0(0.0f, _S1620);
    float _S1623 = s_primal_ctx_sqrt_0(_S1622);
    float invdet_9 = 1.0f / det_blur_6;
    float _S1624 = - _S1613.rows[int(0)].y;
    float _S1625 = - _S1613.rows[int(1)].x;
    float _S1626 = - in_opacity_8;
    float _S1627 = 1.0f + s_primal_ctx_exp_1(_S1626);
    float _S1628 = 1.0f / _S1627;
    float _S1629 = _S1627 * _S1627;
    float _S1630;
    if(antialiased_8)
    {
        _S1630 = _S1628 * _S1623;
    }
    else
    {
        _S1630 = _S1628;
    }
    float _S1631 = _S1630 / 0.00392156885936856f;
    float _S1632 = 2.0f * s_primal_ctx_log_0(_S1631);
    float _S1633 = s_primal_ctx_sqrt_0(_S1632);
    float _S1634 = _S1618.rows[int(0)].x;
    float _S1635 = _S1619.rows[int(1)].y;
    float3  _S1636 = mean_8 - - s_primal_ctx_mul_1(_S1609, t_11);
    float _S1637 = _S1636.x;
    float _S1638 = _S1636.y;
    float _S1639 = _S1636.z;
    float _S1640 = _S1637 * _S1637 + _S1638 * _S1638 + _S1639 * _S1639;
    float _S1641 = s_primal_ctx_sqrt_0(_S1640);
    float x_40 = _S1637 / _S1641;
    float3  _S1642 = make_float3 (x_40);
    float _S1643 = _S1641 * _S1641;
    float y_23 = _S1638 / _S1641;
    float z_20 = _S1639 / _S1641;
    float3  _S1644 = make_float3 (z_20);
    float _S1645 = - y_23;
    float3  _S1646 = make_float3 (_S1645);
    float z2_20 = z_20 * z_20;
    float fTmp0B_8 = -1.09254848957061768f * z_20;
    float fC1_8 = x_40 * x_40 - y_23 * y_23;
    float _S1647 = 2.0f * x_40;
    float fS1_8 = _S1647 * y_23;
    float pSH6_2 = 0.94617468118667603f * z2_20 - 0.31539157032966614f;
    float3  _S1648 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_8 * x_40;
    float3  _S1649 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_8 * y_23;
    float3  _S1650 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_8;
    float3  _S1651 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_8;
    float3  _S1652 = make_float3 (pSH4_2);
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_20;
    float _S1653 = 1.86588168144226074f * z2_20 - 1.11952900886535645f;
    float pSH12_2 = z_20 * _S1653;
    float3  _S1654 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_8 * x_40;
    float3  _S1655 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_8 * y_23;
    float3  _S1656 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_8 * fC1_8;
    float3  _S1657 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_8 * fS1_8;
    float3  _S1658 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_40 * fC1_8 - y_23 * fS1_8);
    float3  _S1659 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_40 * fS1_8 + y_23 * fC1_8);
    float3  _S1660 = make_float3 (pSH9_2);
    float3  _S1661 = make_float3 (0.0f);
    float3  _S1662 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1663;
    (&_S1663)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1645) * (*sh_coeffs_8)[int(1)] + make_float3 (z_20) * (*sh_coeffs_8)[int(2)] - make_float3 (x_40) * (*sh_coeffs_8)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_8)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_8)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_8)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_8)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_8)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_8)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_8)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_8)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_8)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_8)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_8)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f);
    (&_S1663)->differential_0 = _S1662;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1664;
    (&_S1664)->primal_0 = _S1661;
    (&_S1664)->differential_0 = _S1662;
    s_bwd_prop_max_0(&_S1663, &_S1664, v_rgb_2);
    float3  _S1665 = _S1659 * _S1663.differential_0;
    float3  _S1666 = (*sh_coeffs_8)[int(15)] * _S1663.differential_0;
    float3  _S1667 = _S1657 * _S1663.differential_0;
    float3  _S1668 = (*sh_coeffs_8)[int(14)] * _S1663.differential_0;
    float3  _S1669 = _S1655 * _S1663.differential_0;
    float3  _S1670 = (*sh_coeffs_8)[int(13)] * _S1663.differential_0;
    float3  _S1671 = _S1654 * _S1663.differential_0;
    float3  _S1672 = (*sh_coeffs_8)[int(12)] * _S1663.differential_0;
    float3  _S1673 = _S1656 * _S1663.differential_0;
    float3  _S1674 = (*sh_coeffs_8)[int(11)] * _S1663.differential_0;
    float3  _S1675 = _S1658 * _S1663.differential_0;
    float3  _S1676 = (*sh_coeffs_8)[int(10)] * _S1663.differential_0;
    float3  _S1677 = _S1660 * _S1663.differential_0;
    float3  _S1678 = (*sh_coeffs_8)[int(9)] * _S1663.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S1678.x + _S1678.y + _S1678.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S1666.x + _S1666.y + _S1666.z);
    float _S1679 = _S1676.x + _S1676.y + _S1676.z;
    float _S1680 = _S1668.x + _S1668.y + _S1668.z;
    float _S1681 = _S1674.x + _S1674.y + _S1674.z;
    float _S1682 = _S1670.x + _S1670.y + _S1670.z;
    float _S1683 = _S1672.x + _S1672.y + _S1672.z;
    float _S1684 = - s_diff_fC2_T_2;
    float3  _S1685 = _S1651 * _S1663.differential_0;
    float3  _S1686 = (*sh_coeffs_8)[int(8)] * _S1663.differential_0;
    float3  _S1687 = _S1649 * _S1663.differential_0;
    float3  _S1688 = (*sh_coeffs_8)[int(7)] * _S1663.differential_0;
    float3  _S1689 = _S1648 * _S1663.differential_0;
    float3  _S1690 = (*sh_coeffs_8)[int(6)] * _S1663.differential_0;
    float3  _S1691 = _S1650 * _S1663.differential_0;
    float3  _S1692 = (*sh_coeffs_8)[int(5)] * _S1663.differential_0;
    float3  _S1693 = _S1652 * _S1663.differential_0;
    float3  _S1694 = (*sh_coeffs_8)[int(4)] * _S1663.differential_0;
    float _S1695 = _S1692.x + _S1692.y + _S1692.z;
    float _S1696 = _S1688.x + _S1688.y + _S1688.z;
    float _S1697 = fTmp1B_8 * _S1679 + x_40 * s_diff_fS2_T_2 + y_23 * _S1684 + 0.54627424478530884f * (_S1694.x + _S1694.y + _S1694.z);
    float _S1698 = fTmp1B_8 * _S1680 + y_23 * s_diff_fS2_T_2 + x_40 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S1686.x + _S1686.y + _S1686.z);
    float _S1699 = y_23 * - _S1698;
    float _S1700 = x_40 * _S1698;
    float _S1701 = z_20 * (1.86588168144226074f * (z_20 * _S1683) + -2.28522896766662598f * (y_23 * _S1681 + x_40 * _S1682) + 0.94617468118667603f * (_S1690.x + _S1690.y + _S1690.z));
    float3  _S1702 = make_float3 (0.48860251903533936f) * _S1663.differential_0;
    float3  _S1703 = - _S1702;
    float3  _S1704 = _S1642 * _S1703;
    float3  _S1705 = (*sh_coeffs_8)[int(3)] * _S1703;
    float3  _S1706 = _S1644 * _S1702;
    float3  _S1707 = (*sh_coeffs_8)[int(2)] * _S1702;
    float3  _S1708 = _S1646 * _S1702;
    float3  _S1709 = (*sh_coeffs_8)[int(1)] * _S1702;
    float _S1710 = (_S1653 * _S1683 + 1.44530570507049561f * (fS1_8 * _S1679 + fC1_8 * _S1680) + -1.09254848957061768f * (y_23 * _S1695 + x_40 * _S1696) + _S1701 + _S1701 + _S1707.x + _S1707.y + _S1707.z) / _S1643;
    float _S1711 = _S1641 * _S1710;
    float _S1712 = (fTmp0C_8 * _S1681 + fC1_8 * s_diff_fS2_T_2 + fS1_8 * _S1684 + fTmp0B_8 * _S1695 + _S1647 * _S1697 + _S1699 + _S1699 + - (_S1709.x + _S1709.y + _S1709.z)) / _S1643;
    float _S1713 = _S1641 * _S1712;
    float _S1714 = (fTmp0C_8 * _S1682 + fS1_8 * s_diff_fS2_T_2 + fC1_8 * s_diff_fC2_T_2 + fTmp0B_8 * _S1696 + 2.0f * (y_23 * _S1697) + _S1700 + _S1700 + _S1705.x + _S1705.y + _S1705.z) / _S1643;
    float _S1715 = _S1641 * _S1714;
    float _S1716 = _S1639 * - _S1710 + _S1638 * - _S1712 + _S1637 * - _S1714;
    DiffPair_float_0 _S1717;
    (&_S1717)->primal_0 = _S1640;
    (&_S1717)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1717, _S1716);
    float _S1718 = _S1639 * _S1717.differential_0;
    float _S1719 = _S1638 * _S1717.differential_0;
    float _S1720 = _S1637 * _S1717.differential_0;
    float3  _S1721 = make_float3 (0.282094806432724f) * _S1663.differential_0;
    float3  _S1722 = make_float3 (_S1715 + _S1720 + _S1720, _S1713 + _S1719 + _S1719, _S1711 + _S1718 + _S1718);
    float3  _S1723 = - - _S1722;
    Matrix<float, 3, 3>  _S1724 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1725;
    (&_S1725)->primal_0 = _S1609;
    (&_S1725)->differential_0 = _S1724;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1726;
    (&_S1726)->primal_0 = t_11;
    (&_S1726)->differential_0 = _S1662;
    s_bwd_prop_mul_1(&_S1725, &_S1726, _S1723);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1727 = _S1725;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1728 = _S1726;
    float2  _S1729 = make_float2 (0.0f);
    float2  _S1730 = _S1729;
    *&((&_S1730)->y) = v_conic_2.z;
    float2  _S1731 = _S1729;
    *&((&_S1731)->y) = v_conic_2.y;
    *&((&_S1731)->x) = v_conic_2.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1732;
    (&_S1732)->primal_0 = mean_c_8;
    (&_S1732)->differential_0 = _S1662;
    s_bwd_length_impl_0(&_S1732, v_depth_2);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1733 = _S1732;
    DiffPair_float_0 _S1734;
    (&_S1734)->primal_0 = _S1635;
    (&_S1734)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1734, 0.0f);
    DiffPair_float_0 _S1735;
    (&_S1735)->primal_0 = _S1634;
    (&_S1735)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1735, 0.0f);
    DiffPair_float_0 _S1736;
    (&_S1736)->primal_0 = 3.32999992370605469f;
    (&_S1736)->differential_0 = 0.0f;
    DiffPair_float_0 _S1737;
    (&_S1737)->primal_0 = _S1633;
    (&_S1737)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1736, &_S1737, 0.0f);
    DiffPair_float_0 _S1738;
    (&_S1738)->primal_0 = _S1632;
    (&_S1738)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1738, _S1737.differential_0);
    float _S1739 = 2.0f * _S1738.differential_0;
    DiffPair_float_0 _S1740;
    (&_S1740)->primal_0 = _S1631;
    (&_S1740)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1740, _S1739);
    float _S1741 = v_opacity_2 + 254.9999847412109375f * _S1740.differential_0;
    Matrix<float, 2, 2>  _S1742 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1743 = _S1742;
    _S1743[int(1)] = _S1730;
    _S1743[int(0)] = _S1731;
    Matrix<float, 2, 2>  _S1744 = _S1743;
    FixedArray<float3 , 16>  _S1745;
    _S1745[int(0)] = _S1662;
    _S1745[int(1)] = _S1662;
    _S1745[int(2)] = _S1662;
    _S1745[int(3)] = _S1662;
    _S1745[int(4)] = _S1662;
    _S1745[int(5)] = _S1662;
    _S1745[int(6)] = _S1662;
    _S1745[int(7)] = _S1662;
    _S1745[int(8)] = _S1662;
    _S1745[int(9)] = _S1662;
    _S1745[int(10)] = _S1662;
    _S1745[int(11)] = _S1662;
    _S1745[int(12)] = _S1662;
    _S1745[int(13)] = _S1662;
    _S1745[int(14)] = _S1662;
    _S1745[int(15)] = _S1662;
    _S1745[int(7)] = _S1687;
    _S1745[int(0)] = _S1721;
    _S1745[int(1)] = _S1708;
    _S1745[int(2)] = _S1706;
    _S1745[int(3)] = _S1704;
    _S1745[int(4)] = _S1693;
    _S1745[int(5)] = _S1691;
    _S1745[int(6)] = _S1689;
    _S1745[int(15)] = _S1665;
    _S1745[int(8)] = _S1685;
    _S1745[int(9)] = _S1677;
    _S1745[int(10)] = _S1675;
    _S1745[int(11)] = _S1673;
    _S1745[int(12)] = _S1671;
    _S1745[int(13)] = _S1669;
    _S1745[int(14)] = _S1667;
    float3  _S1746 = _S1745[int(0)];
    float3  _S1747 = _S1745[int(1)];
    float3  _S1748 = _S1745[int(2)];
    float3  _S1749 = _S1745[int(3)];
    float3  _S1750 = _S1745[int(4)];
    float3  _S1751 = _S1745[int(5)];
    float3  _S1752 = _S1745[int(6)];
    float3  _S1753 = _S1745[int(7)];
    float3  _S1754 = _S1745[int(8)];
    float3  _S1755 = _S1745[int(9)];
    float3  _S1756 = _S1745[int(10)];
    float3  _S1757 = _S1745[int(11)];
    float3  _S1758 = _S1745[int(12)];
    float3  _S1759 = _S1745[int(13)];
    float3  _S1760 = _S1745[int(14)];
    float3  _S1761 = _S1745[int(15)];
    float2  _S1762 = make_float2 (0.0f, _S1734.differential_0);
    float2  _S1763 = make_float2 (_S1735.differential_0, 0.0f);
    float _S1764;
    if(antialiased_8)
    {
        float _S1765 = _S1628 * _S1741;
        _S1630 = _S1623 * _S1741;
        _S1764 = _S1765;
    }
    else
    {
        _S1630 = _S1741;
        _S1764 = 0.0f;
    }
    float _S1766 = - (_S1630 / _S1629);
    DiffPair_float_0 _S1767;
    (&_S1767)->primal_0 = _S1626;
    (&_S1767)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1767, _S1766);
    float _S1768 = - _S1767.differential_0;
    float _S1769 = invdet_9 * _S1744.rows[int(1)].y;
    float _S1770 = - (invdet_9 * _S1744.rows[int(1)].x);
    float _S1771 = - (invdet_9 * _S1744.rows[int(0)].y);
    float _S1772 = invdet_9 * _S1744.rows[int(0)].x;
    float _S1773 = - ((_S1615 * _S1744.rows[int(1)].y + _S1625 * _S1744.rows[int(1)].x + _S1624 * _S1744.rows[int(0)].y + _S1617 * _S1744.rows[int(0)].x) / _S1621);
    DiffPair_float_0 _S1774;
    (&_S1774)->primal_0 = _S1622;
    (&_S1774)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1774, _S1764);
    DiffPair_float_0 _S1775;
    (&_S1775)->primal_0 = 0.0f;
    (&_S1775)->differential_0 = 0.0f;
    DiffPair_float_0 _S1776;
    (&_S1776)->primal_0 = _S1620;
    (&_S1776)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1775, &_S1776, _S1774.differential_0);
    float _S1777 = _S1776.differential_0 / _S1621;
    float s_diff_det_orig_T_2 = det_blur_6 * _S1777;
    float _S1778 = _S1773 + det_orig_9 * - _S1777;
    float _S1779 = - _S1778;
    float _S1780 = _S1615 * _S1778;
    float _S1781 = _S1617 * _S1778;
    Matrix<float, 2, 2>  _S1782 = _S1742;
    _S1782[int(1)] = _S1762;
    _S1782[int(0)] = _S1763;
    _S1616 = _S1782;
    *&(((&_S1616)->rows + (int(1)))->y) = 0.0f;
    float _S1783 = _S1772 + _S1780 + _S1782.rows[int(1)].y;
    *&(((&_S1616)->rows + (int(0)))->x) = 0.0f;
    float _S1784 = _S1769 + _S1781 + _S1782.rows[int(0)].x;
    float _S1785 = _S1779 + - s_diff_det_orig_T_2;
    float _S1786 = _S1770 + _S1613.rows[int(0)].y * _S1785;
    float _S1787 = _S1771 + _S1613.rows[int(1)].x * _S1785;
    float _S1788 = _S1613.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S1789 = _S1783 + _S1613.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S1790 = _S1729;
    *&((&_S1790)->x) = _S1786;
    *&((&_S1790)->y) = _S1789;
    float _S1791 = _S1784 + _S1788;
    float2  _S1792 = _S1729;
    *&((&_S1792)->y) = _S1787;
    *&((&_S1792)->x) = _S1791;
    float _S1793 = fy_12 * v_mean2d_2.y;
    float _S1794 = fx_12 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S1795 = _S1742;
    _S1795[int(1)] = _S1790;
    _S1795[int(0)] = _S1792;
    Matrix<float, 2, 2>  _S1796 = _S1616 + _S1795;
    Matrix<float, 2, 3>  _S1797 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1798;
    (&_S1798)->primal_0 = _S1611;
    (&_S1798)->differential_0 = _S1797;
    Matrix<float, 3, 2>  _S1799 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1800;
    (&_S1800)->primal_0 = _S1612;
    (&_S1800)->differential_0 = _S1799;
    s_bwd_prop_mul_2(&_S1798, &_S1800, _S1796);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1801;
    (&_S1801)->primal_0 = J_12;
    (&_S1801)->differential_0 = _S1797;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1802;
    (&_S1802)->primal_0 = _S1610;
    (&_S1802)->differential_0 = _S1724;
    s_bwd_prop_mul_3(&_S1801, &_S1802, _S1798.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1803;
    (&_S1803)->primal_0 = _S1608;
    (&_S1803)->differential_0 = _S1724;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1804;
    (&_S1804)->primal_0 = _S1609;
    (&_S1804)->differential_0 = _S1724;
    s_bwd_prop_mul_4(&_S1803, &_S1804, _S1802.differential_0);
    Matrix<float, 3, 3>  _S1805 = transpose_0(_S1804.differential_0 + _S1727.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1806;
    (&_S1806)->primal_0 = R_12;
    (&_S1806)->differential_0 = _S1724;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1807;
    (&_S1807)->primal_0 = _S1607;
    (&_S1807)->differential_0 = _S1724;
    s_bwd_prop_mul_4(&_S1806, &_S1807, _S1803.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1808;
    (&_S1808)->primal_0 = _S1605;
    (&_S1808)->differential_0 = _S1724;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1809;
    (&_S1809)->primal_0 = _S1606;
    (&_S1809)->differential_0 = _S1724;
    s_bwd_prop_mul_4(&_S1808, &_S1809, _S1807.differential_0);
    Matrix<float, 3, 3>  _S1810 = _S1808.differential_0 + transpose_0(_S1809.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1811;
    (&_S1811)->primal_0 = _S1604;
    (&_S1811)->differential_0 = _S1724;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1812;
    (&_S1812)->primal_0 = S_2;
    (&_S1812)->differential_0 = _S1724;
    s_bwd_prop_mul_4(&_S1811, &_S1812, _S1810);
    Matrix<float, 3, 3>  _S1813 = transpose_0(_S1811.differential_0);
    float _S1814 = 2.0f * - _S1813.rows[int(2)].z;
    float _S1815 = 2.0f * _S1813.rows[int(2)].y;
    float _S1816 = 2.0f * _S1813.rows[int(2)].x;
    float _S1817 = 2.0f * _S1813.rows[int(1)].z;
    float _S1818 = 2.0f * - _S1813.rows[int(1)].y;
    float _S1819 = 2.0f * _S1813.rows[int(1)].x;
    float _S1820 = 2.0f * _S1813.rows[int(0)].z;
    float _S1821 = 2.0f * _S1813.rows[int(0)].y;
    float _S1822 = 2.0f * - _S1813.rows[int(0)].x;
    float _S1823 = - _S1819 + _S1821;
    float _S1824 = _S1816 + - _S1820;
    float _S1825 = - _S1815 + _S1817;
    float _S1826 = _S1815 + _S1817;
    float _S1827 = _S1816 + _S1820;
    float _S1828 = _S1819 + _S1821;
    float _S1829 = z_19 * (_S1818 + _S1822);
    float _S1830 = y_22 * (_S1814 + _S1822);
    float _S1831 = x_39 * (_S1814 + _S1818);
    float _S1832 = z_19 * _S1823 + y_22 * _S1824 + x_39 * _S1825;
    float _S1833 = _S1603 * _S1832;
    float _S1834 = w_11 * _S1823 + y_22 * _S1826 + x_39 * _S1827 + _S1829 + _S1829;
    float _S1835 = _S1603 * _S1834;
    float _S1836 = w_11 * _S1824 + z_19 * _S1826 + x_39 * _S1828 + _S1830 + _S1830;
    float _S1837 = _S1603 * _S1836;
    float _S1838 = w_11 * _S1825 + z_19 * _S1827 + y_22 * _S1828 + _S1831 + _S1831;
    float _S1839 = _S1603 * _S1838;
    float _S1840 = quat_11.x * _S1832 + quat_11.w * _S1834 + quat_11.z * _S1836 + quat_11.y * _S1838;
    DiffPair_float_0 _S1841;
    (&_S1841)->primal_0 = _S1602;
    (&_S1841)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1841, _S1840);
    float _S1842 = quat_11.x * _S1841.differential_0;
    float _S1843 = quat_11.w * _S1841.differential_0;
    float _S1844 = quat_11.z * _S1841.differential_0;
    float _S1845 = quat_11.y * _S1841.differential_0;
    float _S1846 = _S1835 + _S1843 + _S1843;
    float _S1847 = _S1837 + _S1844 + _S1844;
    float _S1848 = _S1839 + _S1845 + _S1845;
    float _S1849 = _S1833 + _S1842 + _S1842;
    float3  _S1850 = _S1662;
    *&((&_S1850)->z) = _S1812.differential_0.rows[int(2)].z;
    *&((&_S1850)->y) = _S1812.differential_0.rows[int(1)].y;
    *&((&_S1850)->x) = _S1812.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1851;
    (&_S1851)->primal_0 = scale_10;
    (&_S1851)->differential_0 = _S1662;
    s_bwd_prop_exp_1(&_S1851, _S1850);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1852 = _S1851;
    float3  _S1853 = _S1662;
    *&((&_S1853)->y) = _S1793;
    *&((&_S1853)->x) = _S1794;
    float3  _S1854 = _S1733.differential_0 + _S1853;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1855;
    (&_S1855)->primal_0 = R_12;
    (&_S1855)->differential_0 = _S1724;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1856;
    (&_S1856)->primal_0 = mean_8;
    (&_S1856)->differential_0 = _S1662;
    s_bwd_prop_mul_1(&_S1855, &_S1856, _S1854);
    float3  _S1857 = _S1854 + _S1728.differential_0;
    Matrix<float, 3, 3>  _S1858 = _S1805 + _S1806.differential_0 + _S1855.differential_0;
    float4  _S1859 = make_float4 (0.0f);
    *&((&_S1859)->w) = _S1846;
    *&((&_S1859)->z) = _S1847;
    *&((&_S1859)->y) = _S1848;
    *&((&_S1859)->x) = _S1849;
    float4  _S1860 = _S1859;
    float3  _S1861 = _S1856.differential_0 + _S1722;
    *v_mean_2 = _S1861;
    *v_quat_2 = _S1860;
    *v_scale_2 = _S1852.differential_0;
    *v_in_opacity_2 = _S1768;
    (*v_sh_coeffs_2)[int(0)] = _S1746;
    (*v_sh_coeffs_2)[int(1)] = _S1747;
    (*v_sh_coeffs_2)[int(2)] = _S1748;
    (*v_sh_coeffs_2)[int(3)] = _S1749;
    (*v_sh_coeffs_2)[int(4)] = _S1750;
    (*v_sh_coeffs_2)[int(5)] = _S1751;
    (*v_sh_coeffs_2)[int(6)] = _S1752;
    (*v_sh_coeffs_2)[int(7)] = _S1753;
    (*v_sh_coeffs_2)[int(8)] = _S1754;
    (*v_sh_coeffs_2)[int(9)] = _S1755;
    (*v_sh_coeffs_2)[int(10)] = _S1756;
    (*v_sh_coeffs_2)[int(11)] = _S1757;
    (*v_sh_coeffs_2)[int(12)] = _S1758;
    (*v_sh_coeffs_2)[int(13)] = _S1759;
    (*v_sh_coeffs_2)[int(14)] = _S1760;
    (*v_sh_coeffs_2)[int(15)] = _S1761;
    *v_R_3 = _S1858;
    *v_t_3 = _S1857;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_16)
{
    float _S1862 = dOut_16.y;
    float _S1863 = dOut_16.z;
    float _S1864 = dOut_16.x;
    float _S1865 = (*a_0).primal_0.z * _S1862 + - (*a_0).primal_0.y * _S1863;
    float _S1866 = - (*a_0).primal_0.z * _S1864 + (*a_0).primal_0.x * _S1863;
    float _S1867 = (*a_0).primal_0.y * _S1864 + - (*a_0).primal_0.x * _S1862;
    float3  _S1868 = make_float3 (- (*b_0).primal_0.z * _S1862 + (*b_0).primal_0.y * _S1863, (*b_0).primal_0.z * _S1864 + - (*b_0).primal_0.x * _S1863, - (*b_0).primal_0.y * _S1864 + (*b_0).primal_0.x * _S1862);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S1868;
    float3  _S1869 = make_float3 (_S1865, _S1866, _S1867);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S1869;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S1870 = left_10.y;
    float _S1871 = right_10.z;
    float _S1872 = left_10.z;
    float _S1873 = right_10.y;
    float _S1874 = right_10.x;
    float _S1875 = left_10.x;
    return make_float3 (_S1870 * _S1871 - _S1872 * _S1873, _S1872 * _S1874 - _S1875 * _S1871, _S1875 * _S1873 - _S1870 * _S1874);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_9, float4  quat_12, float3  scale_11, float opacity_6, float3  ray_o_1, float3  ray_d_1)
{
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
    Matrix<float, 3, 3>  iscl_rot_0 = mul_4(makeMatrix<float, 3, 3> ((F32_exp((- scale_11.x))), 0.0f, 0.0f, 0.0f, (F32_exp((- scale_11.y))), 0.0f, 0.0f, 0.0f, (F32_exp((- scale_11.z)))), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12)))));
    float3  gcrod_0 = cross_0(normalize_0(mul_0(iscl_rot_0, ray_d_1)), mul_0(iscl_rot_0, ray_o_1 - mean_9));
    return (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0)))) / (1.0f + (F32_exp((- opacity_6))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S1876, float3  _S1877)
{
    return cross_0(_S1876, _S1877);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S1878, float3  _S1879)
{
    return dot_0(_S1878, _S1879);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1880, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1881, float _S1882)
{
    _d_dot_0(_S1880, _S1881, _S1882);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1883, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1884, float3  _S1885)
{
    _d_cross_0(_S1883, _S1884, _S1885);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_7)
{
    float _S1886 = (*dpquat_0).primal_0.y;
    float _S1887 = _S1886 * _S1886 + (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z + (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w + (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.x;
    float _S1888 = s_primal_ctx_rsqrt_0(_S1887);
    float x_43 = (*dpquat_0).primal_0.y * _S1888;
    float y_25 = (*dpquat_0).primal_0.z * _S1888;
    float z_22 = (*dpquat_0).primal_0.w * _S1888;
    float w_13 = (*dpquat_0).primal_0.x * _S1888;
    float x2_13 = x_43 * x_43;
    float y2_13 = y_25 * y_25;
    float z2_22 = z_22 * z_22;
    float xy_13 = x_43 * y_25;
    float xz_13 = x_43 * z_22;
    float yz_13 = y_25 * z_22;
    float wx_13 = w_13 * x_43;
    float wy_13 = w_13 * y_25;
    float wz_13 = w_13 * z_22;
    float _S1889 = - (*dpscale_0).primal_0.x;
    float _S1890 = - (*dpscale_0).primal_0.y;
    float _S1891 = - (*dpscale_0).primal_0.z;
    Matrix<float, 3, 3>  S_3 = makeMatrix<float, 3, 3> (s_primal_ctx_exp_1(_S1889), 0.0f, 0.0f, 0.0f, s_primal_ctx_exp_1(_S1890), 0.0f, 0.0f, 0.0f, s_primal_ctx_exp_1(_S1891));
    Matrix<float, 3, 3>  _S1892 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_22), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_22), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13))));
    Matrix<float, 3, 3>  _S1893 = s_primal_ctx_mul_2(S_3, _S1892);
    float3  _S1894 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S1895 = s_primal_ctx_mul_1(_S1893, _S1894);
    float3  _S1896 = s_primal_ctx_mul_1(_S1893, (*dpray_d_2).primal_0);
    float3  _S1897 = normalize_0(_S1896);
    float3  _S1898 = s_primal_ctx_cross_0(_S1897, _S1895);
    float _S1899 = -0.5f * s_primal_ctx_dot_0(_S1898, _S1898);
    float _S1900 = - (*dpopacity_0).primal_0;
    float _S1901 = 1.0f + s_primal_ctx_exp_1(_S1900);
    float _S1902 = _s_dOut_7 / (_S1901 * _S1901);
    float _S1903 = s_primal_ctx_exp_1(_S1899) * - _S1902;
    float _S1904 = _S1901 * _S1902;
    DiffPair_float_0 _S1905;
    (&_S1905)->primal_0 = _S1900;
    (&_S1905)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1905, _S1903);
    float _S1906 = - _S1905.differential_0;
    DiffPair_float_0 _S1907;
    (&_S1907)->primal_0 = _S1899;
    (&_S1907)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1907, _S1904);
    float _S1908 = -0.5f * _S1907.differential_0;
    float3  _S1909 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1910;
    (&_S1910)->primal_0 = _S1898;
    (&_S1910)->differential_0 = _S1909;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1911;
    (&_S1911)->primal_0 = _S1898;
    (&_S1911)->differential_0 = _S1909;
    s_bwd_prop_dot_0(&_S1910, &_S1911, _S1908);
    float3  _S1912 = _S1911.differential_0 + _S1910.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1913;
    (&_S1913)->primal_0 = _S1897;
    (&_S1913)->differential_0 = _S1909;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1914;
    (&_S1914)->primal_0 = _S1895;
    (&_S1914)->differential_0 = _S1909;
    s_bwd_prop_cross_0(&_S1913, &_S1914, _S1912);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1915;
    (&_S1915)->primal_0 = _S1896;
    (&_S1915)->differential_0 = _S1909;
    s_bwd_normalize_impl_0(&_S1915, _S1913.differential_0);
    Matrix<float, 3, 3>  _S1916 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1917;
    (&_S1917)->primal_0 = _S1893;
    (&_S1917)->differential_0 = _S1916;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1918;
    (&_S1918)->primal_0 = (*dpray_d_2).primal_0;
    (&_S1918)->differential_0 = _S1909;
    s_bwd_prop_mul_1(&_S1917, &_S1918, _S1915.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1919;
    (&_S1919)->primal_0 = _S1893;
    (&_S1919)->differential_0 = _S1916;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1920;
    (&_S1920)->primal_0 = _S1894;
    (&_S1920)->differential_0 = _S1909;
    s_bwd_prop_mul_1(&_S1919, &_S1920, _S1914.differential_0);
    float3  _S1921 = - _S1920.differential_0;
    Matrix<float, 3, 3>  _S1922 = _S1917.differential_0 + _S1919.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1923;
    (&_S1923)->primal_0 = S_3;
    (&_S1923)->differential_0 = _S1916;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1924;
    (&_S1924)->primal_0 = _S1892;
    (&_S1924)->differential_0 = _S1916;
    s_bwd_prop_mul_4(&_S1923, &_S1924, _S1922);
    Matrix<float, 3, 3>  _S1925 = transpose_0(_S1924.differential_0);
    DiffPair_float_0 _S1926;
    (&_S1926)->primal_0 = _S1891;
    (&_S1926)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1926, _S1923.differential_0.rows[int(2)].z);
    float _S1927 = - _S1926.differential_0;
    DiffPair_float_0 _S1928;
    (&_S1928)->primal_0 = _S1890;
    (&_S1928)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1928, _S1923.differential_0.rows[int(1)].y);
    float _S1929 = - _S1928.differential_0;
    DiffPair_float_0 _S1930;
    (&_S1930)->primal_0 = _S1889;
    (&_S1930)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1930, _S1923.differential_0.rows[int(0)].x);
    float _S1931 = - _S1930.differential_0;
    Matrix<float, 3, 3>  _S1932 = transpose_0(_S1925);
    float _S1933 = 2.0f * - _S1932.rows[int(2)].z;
    float _S1934 = 2.0f * _S1932.rows[int(2)].y;
    float _S1935 = 2.0f * _S1932.rows[int(2)].x;
    float _S1936 = 2.0f * _S1932.rows[int(1)].z;
    float _S1937 = 2.0f * - _S1932.rows[int(1)].y;
    float _S1938 = 2.0f * _S1932.rows[int(1)].x;
    float _S1939 = 2.0f * _S1932.rows[int(0)].z;
    float _S1940 = 2.0f * _S1932.rows[int(0)].y;
    float _S1941 = 2.0f * - _S1932.rows[int(0)].x;
    float _S1942 = - _S1938 + _S1940;
    float _S1943 = _S1935 + - _S1939;
    float _S1944 = - _S1934 + _S1936;
    float _S1945 = _S1934 + _S1936;
    float _S1946 = _S1935 + _S1939;
    float _S1947 = _S1938 + _S1940;
    float _S1948 = z_22 * (_S1937 + _S1941);
    float _S1949 = y_25 * (_S1933 + _S1941);
    float _S1950 = x_43 * (_S1933 + _S1937);
    float _S1951 = z_22 * _S1942 + y_25 * _S1943 + x_43 * _S1944;
    float _S1952 = _S1888 * _S1951;
    float _S1953 = w_13 * _S1942 + y_25 * _S1945 + x_43 * _S1946 + _S1948 + _S1948;
    float _S1954 = _S1888 * _S1953;
    float _S1955 = w_13 * _S1943 + z_22 * _S1945 + x_43 * _S1947 + _S1949 + _S1949;
    float _S1956 = _S1888 * _S1955;
    float _S1957 = w_13 * _S1944 + z_22 * _S1946 + y_25 * _S1947 + _S1950 + _S1950;
    float _S1958 = _S1888 * _S1957;
    float _S1959 = (*dpquat_0).primal_0.x * _S1951 + (*dpquat_0).primal_0.w * _S1953 + (*dpquat_0).primal_0.z * _S1955 + (*dpquat_0).primal_0.y * _S1957;
    DiffPair_float_0 _S1960;
    (&_S1960)->primal_0 = _S1887;
    (&_S1960)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1960, _S1959);
    float _S1961 = (*dpquat_0).primal_0.x * _S1960.differential_0;
    float _S1962 = (*dpquat_0).primal_0.w * _S1960.differential_0;
    float _S1963 = (*dpquat_0).primal_0.z * _S1960.differential_0;
    float _S1964 = (*dpquat_0).primal_0.y * _S1960.differential_0;
    float _S1965 = _S1954 + _S1962 + _S1962;
    float _S1966 = _S1956 + _S1963 + _S1963;
    float _S1967 = _S1958 + _S1964 + _S1964;
    float _S1968 = _S1952 + _S1961 + _S1961;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S1918.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S1920.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S1906;
    float3  _S1969 = make_float3 (_S1931, _S1929, _S1927);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S1969;
    float4  _S1970 = make_float4 (0.0f);
    *&((&_S1970)->w) = _S1965;
    *&((&_S1970)->z) = _S1966;
    *&((&_S1970)->y) = _S1967;
    *&((&_S1970)->x) = _S1968;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S1970;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S1921;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1971, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1972, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1973, DiffPair_float_0 * _S1974, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1975, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1976, float _S1977)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S1971, _S1972, _S1973, _S1974, _S1975, _S1976, _S1977);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_10, float4  quat_13, float3  scale_12, float opacity_7, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_opacity_3, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S1978 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_10;
    (&dp_mean_0)->differential_0 = _S1978;
    float4  _S1979 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_13;
    (&dp_quat_0)->differential_0 = _S1979;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_12;
    (&dp_scale_0)->differential_0 = _S1978;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_7;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S1978;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S1978;
    s_bwd_evaluate_alpha_3dgs_0(&dp_mean_0, &dp_quat_0, &dp_scale_0, &dp_opacity_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_mean_3 = dp_mean_0.differential_0;
    *v_quat_3 = dp_quat_0.differential_0;
    *v_scale_3 = dp_scale_0.differential_0;
    *v_opacity_3 = dp_opacity_0.differential_0;
    *v_ray_o_1 = dp_ray_o_0.differential_0;
    *v_ray_d_1 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_3dgs(float3  mean_11, float4  quat_14, float3  scale_13, float opacity_8, float3  rgb_6, float3  ray_o_3, float3  ray_d_3, float3  * out_rgb_0, float * depth_6)
{
    *out_rgb_0 = rgb_6;
    *depth_6 = length_1(mean_11 - ray_o_3);
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S1980 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1981;
    (&_S1981)->primal_0 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    (&_S1981)->differential_0 = _S1980;
    s_bwd_length_impl_0(&_S1981, dpdepth_0);
    float3  _S1982 = - _S1981.differential_0;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S1980;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S1982;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    dpscale_1->primal_0 = (*dpscale_1).primal_0;
    dpscale_1->differential_0 = _S1980;
    float4  _S1983 = make_float4 (0.0f);
    dpquat_1->primal_0 = (*dpquat_1).primal_0;
    dpquat_1->differential_0 = _S1983;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S1981.differential_0;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1984, DiffPair_vectorx3Cfloatx2C4x3E_0 * _S1985, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1986, DiffPair_float_0 * _S1987, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1988, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1989, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1990, float3  _S1991, float _S1992)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S1984, _S1985, _S1986, _S1987, _S1988, _S1989, _S1990, _S1991, _S1992);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_12, float4  quat_15, float3  scale_14, float opacity_9, float3  rgb_7, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_3, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_opacity_4, float3  * v_rgb_3, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S1993 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_12;
    (&dp_mean_1)->differential_0 = _S1993;
    float4  _S1994 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_1;
    (&dp_quat_1)->primal_0 = quat_15;
    (&dp_quat_1)->differential_0 = _S1994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_1;
    (&dp_scale_1)->primal_0 = scale_14;
    (&dp_scale_1)->differential_0 = _S1993;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_9;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_7;
    (&dp_rgb_0)->differential_0 = _S1993;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S1993;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S1993;
    s_bwd_evaluate_color_3dgs_0(&dp_mean_1, &dp_quat_1, &dp_scale_1, &dp_opacity_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_3);
    *v_mean_4 = dp_mean_1.differential_0;
    *v_quat_4 = dp_quat_1.differential_0;
    *v_scale_4 = dp_scale_1.differential_0;
    *v_opacity_4 = dp_opacity_1.differential_0;
    *v_rgb_3 = dp_rgb_0.differential_0;
    *v_ray_o_2 = dp_ray_o_1.differential_0;
    *v_ray_d_2 = dp_ray_d_1.differential_0;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_14, float dOut_17)
{
    float _S1995 = _slang_select(((*dpx_14).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_14).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S1995;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_15, float dOut_18)
{
    float _S1996 = (F32_exp2(((*dpx_15).primal_0))) * 50.693145751953125f * dOut_18;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S1996;
    return;
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_13, float4  quat_16, float3  scale_15, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_9, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_13, float fy_13, float cx_13, float cy_13, float4  radial_coeffs_16, float2  tangential_coeffs_16, float2  thin_prism_coeffs_16, uint image_width_9, uint image_height_9, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_7, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_8, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_9 = mul_0(R_13, mean_13) + t_12;
        float _S1997 = mean_c_9.z;
        bool _S1998;
        if(_S1997 < near_plane_6)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = _S1997 > far_plane_6;
        }
        if(_S1998)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S1999 = scale_15.x;
        float sx_0 = (F32_exp((_S1999)));
        float _S2000 = scale_15.y;
        float sy_0 = (F32_exp((_S2000)));
        float sz_0 = scale_15.z - 0.5f * (_S1999 + _S2000);
        float x_44 = quat_16.y;
        float inv_norm_10 = (F32_rsqrt((x_44 * x_44 + quat_16.z * quat_16.z + quat_16.w * quat_16.w + quat_16.x * quat_16.x)));
        float x_45 = quat_16.y * inv_norm_10;
        float y_26 = quat_16.z * inv_norm_10;
        float z_23 = quat_16.w * inv_norm_10;
        float w_14 = quat_16.x * inv_norm_10;
        float x2_14 = x_45 * x_45;
        float y2_14 = y_26 * y_26;
        float z2_23 = z_23 * z_23;
        float xy_14 = x_45 * y_26;
        float xz_14 = x_45 * z_23;
        float yz_14 = y_26 * z_23;
        float wx_14 = w_14 * x_45;
        float wy_14 = w_14 * y_26;
        float wz_14 = w_14 * z_23;
        Matrix<float, 3, 3>  _S2001 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_23), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_23), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
        float3  vert0_c_0 = mul_0(R_13, mul_0(_S2001, make_float3 (sx_0, 0.0f, 0.0f)) + mean_13) + t_12;
        float3  vert1_c_0 = mul_0(R_13, mul_0(_S2001, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_13) + t_12;
        float3  vert2_c_0 = mul_0(R_13, mul_0(_S2001, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_13) + t_12;
        float _S2002 = vert0_c_0.z;
        float _S2003 = vert1_c_0.z;
        float _S2004 = vert2_c_0.z;
        if(_S2002 < near_plane_6)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = _S2002 > far_plane_6;
        }
        if(_S1998)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = _S2003 < near_plane_6;
        }
        if(_S1998)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = _S2003 > far_plane_6;
        }
        if(_S1998)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = _S2004 < near_plane_6;
        }
        if(_S1998)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = _S2004 > far_plane_6;
        }
        if(_S1998)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S2002);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S2003);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S2004);
        float2  _S2005 = make_float2 (fx_13, fy_13);
        float2  _S2006 = make_float2 (cx_13, cy_13);
        *uv0_0 = _S2005 * *uv0_0 + _S2006;
        *uv1_0 = _S2005 * *uv1_0 + _S2006;
        float2  _S2007 = _S2005 * *uv2_0 + _S2006;
        *uv2_0 = _S2007;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S2007 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S2007)));
        float _S2008 = _S2007.x;
        float xmax_3 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S2008))) + offset_0;
        float xmin_3 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S2008))) - offset_0;
        float _S2009 = _S2007.y;
        float ymax_3 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S2009))) + offset_0;
        float ymin_3 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S2009))) - offset_0;
        if(xmax_3 <= 0.0f)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = xmin_3 >= float(image_width_9);
        }
        if(_S1998)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = ymax_3 <= 0.0f;
        }
        if(_S1998)
        {
            _S1998 = true;
        }
        else
        {
            _S1998 = ymin_3 >= float(image_height_9);
        }
        if(_S1998)
        {
            _S1998 = true;
        }
        else
        {
            if(_S1997 <= 0.0f)
            {
                if(xmin_3 <= 0.0f)
                {
                    _S1998 = xmax_3 >= float(image_width_9);
                }
                else
                {
                    _S1998 = false;
                }
                if(_S1998)
                {
                    _S1998 = true;
                }
                else
                {
                    if(ymin_3 <= 0.0f)
                    {
                        _S1998 = ymax_3 >= float(image_width_9);
                    }
                    else
                    {
                        _S1998 = false;
                    }
                }
            }
            else
            {
                _S1998 = false;
            }
        }
        if(_S1998)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_6 = make_int4 (int((F32_floor((xmin_3)))), int((F32_floor((ymin_3)))), int((F32_ceil((xmax_3)))), int((F32_ceil((ymax_3)))));
        *depth_7 = make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0));
        *out_hardness_0 = hardness_0;
        float3  _S2010 = mean_13 - - mul_0(transpose_0(R_13), t_12);
        float _S2011 = _S2010.x;
        float _S2012 = _S2010.y;
        float _S2013 = _S2010.z;
        float norm_6 = (F32_sqrt((_S2011 * _S2011 + _S2012 * _S2012 + _S2013 * _S2013)));
        float x_46 = _S2011 / norm_6;
        float y_27 = _S2012 / norm_6;
        float z_24 = _S2013 / norm_6;
        float z2_24 = z_24 * z_24;
        float fTmp0B_9 = -1.09254848957061768f * z_24;
        float fC1_9 = x_46 * x_46 - y_27 * y_27;
        float fS1_9 = 2.0f * x_46 * y_27;
        float fTmp0C_9 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
        float fTmp1B_9 = 1.44530570507049561f * z_24;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_9)[int(1)] + make_float3 (z_24) * (*sh_coeffs_9)[int(2)] - make_float3 (x_46) * (*sh_coeffs_9)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_27) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_24 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_46) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_9 + y_27 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_27) * (*sh_coeffs_9)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_24 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_46) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_9 - y_27 * fS1_9)) * (*sh_coeffs_9)[int(15)]);
        float3  _S2014 = make_float3 (0.0f);
        (*rgb_8)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S2014);
        float3  _S2015 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S2016 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_8)[int(1)] = max_0(_S2015 + _S2016 + make_float3 (0.5f), _S2014);
        (*rgb_8)[int(2)] = max_0(_S2015 - _S2016 + make_float3 (0.5f), _S2014);
        float3  _S2017 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S2017 * make_float3 (float(- (F32_sign((dot_0(_S2017, mean_c_9))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_14, float4  quat_17, float3  scale_16, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_10, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_14, float fy_14, float cx_14, float cy_14, float4  radial_coeffs_17, float2  tangential_coeffs_17, float2  thin_prism_coeffs_17, uint image_width_10, uint image_height_10, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_8, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_9, float3  * normal_1)
{
    for(;;)
    {
        float3  mean_c_10 = mul_0(R_14, mean_14) + t_13;
        float _S2018 = length_1(mean_c_10);
        bool _S2019;
        if(_S2018 < near_plane_7)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = _S2018 > far_plane_7;
        }
        if(_S2019)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S2020 = scale_16.x;
        float sx_1 = (F32_exp((_S2020)));
        float _S2021 = scale_16.y;
        float sy_1 = (F32_exp((_S2021)));
        float sz_1 = scale_16.z - 0.5f * (_S2020 + _S2021);
        float x_47 = quat_17.y;
        float inv_norm_11 = (F32_rsqrt((x_47 * x_47 + quat_17.z * quat_17.z + quat_17.w * quat_17.w + quat_17.x * quat_17.x)));
        float x_48 = quat_17.y * inv_norm_11;
        float y_28 = quat_17.z * inv_norm_11;
        float z_25 = quat_17.w * inv_norm_11;
        float w_15 = quat_17.x * inv_norm_11;
        float x2_15 = x_48 * x_48;
        float y2_15 = y_28 * y_28;
        float z2_25 = z_25 * z_25;
        float xy_15 = x_48 * y_28;
        float xz_15 = x_48 * z_25;
        float yz_15 = y_28 * z_25;
        float wx_15 = w_15 * x_48;
        float wy_15 = w_15 * y_28;
        float wz_15 = w_15 * z_25;
        Matrix<float, 3, 3>  _S2022 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_25), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_25), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
        float3  vert0_c_1 = mul_0(R_14, mul_0(_S2022, make_float3 (sx_1, 0.0f, 0.0f)) + mean_14) + t_13;
        float3  vert1_c_1 = mul_0(R_14, mul_0(_S2022, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_14) + t_13;
        float3  vert2_c_1 = mul_0(R_14, mul_0(_S2022, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_14) + t_13;
        float _S2023 = length_1(vert0_c_1);
        float _S2024 = length_1(vert1_c_1);
        float _S2025 = length_1(vert2_c_1);
        if(_S2023 < near_plane_7)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = _S2023 > far_plane_7;
        }
        if(_S2019)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = _S2024 < near_plane_7;
        }
        if(_S2019)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = _S2024 > far_plane_7;
        }
        if(_S2019)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = _S2025 < near_plane_7;
        }
        if(_S2019)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = _S2025 > far_plane_7;
        }
        if(_S2019)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_2 = CameraDistortion_x24init_0(radial_coeffs_17, tangential_coeffs_17, thin_prism_coeffs_17);
        float2  _S2026 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_7 = length_0(_S2026);
        float _S2027 = vert0_c_1.z;
        float theta_2 = (F32_atan2((r_7), (_S2027)));
        float k_4;
        if(theta_2 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_2 * theta_2 / 3.0f) / _S2027;
        }
        else
        {
            k_4 = theta_2 / r_7;
        }
        float2  _S2028 = _S2026 * make_float2 (k_4);
        float k1_3 = dist_coeffs_2.radial_coeffs_0.x;
        float k2_3 = dist_coeffs_2.radial_coeffs_0.y;
        float k3_3 = dist_coeffs_2.radial_coeffs_0.z;
        float k4_4 = dist_coeffs_2.radial_coeffs_0.w;
        float p1_4 = dist_coeffs_2.tangential_coeffs_0.x;
        float p2_4 = dist_coeffs_2.tangential_coeffs_0.y;
        float sx1_4 = dist_coeffs_2.thin_prism_coeffs_0.x;
        float sy1_4 = dist_coeffs_2.thin_prism_coeffs_0.y;
        float u_10 = _S2028.x;
        float v_10 = _S2028.y;
        float r2_10 = u_10 * u_10 + v_10 * v_10;
        float _S2029 = 2.0f * p1_4;
        float _S2030 = 2.0f * p2_4;
        float2  _S2031 = _S2028 * make_float2 (1.0f + r2_10 * (k1_3 + r2_10 * (k2_3 + r2_10 * (k3_3 + r2_10 * k4_4)))) + make_float2 (_S2029 * u_10 * v_10 + p2_4 * (r2_10 + 2.0f * u_10 * u_10) + sx1_4 * r2_10, _S2030 * u_10 * v_10 + p1_4 * (r2_10 + 2.0f * v_10 * v_10) + sy1_4 * r2_10);
        *uv0_1 = make_float2 (fx_14 * _S2031.x + cx_14, fy_14 * _S2031.y + cy_14);
        float2  _S2032 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_8 = length_0(_S2032);
        float _S2033 = vert1_c_1.z;
        float theta_3 = (F32_atan2((r_8), (_S2033)));
        if(theta_3 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_3 * theta_3 / 3.0f) / _S2033;
        }
        else
        {
            k_4 = theta_3 / r_8;
        }
        float2  _S2034 = _S2032 * make_float2 (k_4);
        float u_11 = _S2034.x;
        float v_11 = _S2034.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float2  _S2035 = _S2034 * make_float2 (1.0f + r2_11 * (k1_3 + r2_11 * (k2_3 + r2_11 * (k3_3 + r2_11 * k4_4)))) + make_float2 (_S2029 * u_11 * v_11 + p2_4 * (r2_11 + 2.0f * u_11 * u_11) + sx1_4 * r2_11, _S2030 * u_11 * v_11 + p1_4 * (r2_11 + 2.0f * v_11 * v_11) + sy1_4 * r2_11);
        *uv1_1 = make_float2 (fx_14 * _S2035.x + cx_14, fy_14 * _S2035.y + cy_14);
        float2  _S2036 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_9 = length_0(_S2036);
        float _S2037 = vert2_c_1.z;
        float theta_4 = (F32_atan2((r_9), (_S2037)));
        if(theta_4 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_4 * theta_4 / 3.0f) / _S2037;
        }
        else
        {
            k_4 = theta_4 / r_9;
        }
        float2  _S2038 = _S2036 * make_float2 (k_4);
        float u_12 = _S2038.x;
        float v_12 = _S2038.y;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float2  _S2039 = _S2038 * make_float2 (1.0f + r2_12 * (k1_3 + r2_12 * (k2_3 + r2_12 * (k3_3 + r2_12 * k4_4)))) + make_float2 (_S2029 * u_12 * v_12 + p2_4 * (r2_12 + 2.0f * u_12 * u_12) + sx1_4 * r2_12, _S2030 * u_12 * v_12 + p1_4 * (r2_12 + 2.0f * v_12 * v_12) + sy1_4 * r2_12);
        float _S2040 = fx_14 * _S2039.x + cx_14;
        float _S2041 = fy_14 * _S2039.y + cy_14;
        float2  _S2042 = make_float2 (_S2040, _S2041);
        *uv2_1 = _S2042;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S2042 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S2042)));
        float xmax_4 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S2040))) + offset_1;
        float xmin_4 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S2040))) - offset_1;
        float ymax_4 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S2041))) + offset_1;
        float ymin_4 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S2041))) - offset_1;
        if(xmax_4 <= 0.0f)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = xmin_4 >= float(image_width_10);
        }
        if(_S2019)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = ymax_4 <= 0.0f;
        }
        if(_S2019)
        {
            _S2019 = true;
        }
        else
        {
            _S2019 = ymin_4 >= float(image_height_10);
        }
        if(_S2019)
        {
            _S2019 = true;
        }
        else
        {
            if((mean_c_10.z) <= 0.0f)
            {
                if(xmin_4 <= 0.0f)
                {
                    _S2019 = xmax_4 >= float(image_width_10);
                }
                else
                {
                    _S2019 = false;
                }
                if(_S2019)
                {
                    _S2019 = true;
                }
                else
                {
                    if(ymin_4 <= 0.0f)
                    {
                        _S2019 = ymax_4 >= float(image_width_10);
                    }
                    else
                    {
                        _S2019 = false;
                    }
                }
            }
            else
            {
                _S2019 = false;
            }
        }
        if(_S2019)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_7 = make_int4 (int((F32_floor((xmin_4)))), int((F32_floor((ymin_4)))), int((F32_ceil((xmax_4)))), int((F32_ceil((ymax_4)))));
        *depth_8 = make_float3 (_S2023, _S2024, _S2025);
        *out_hardness_1 = hardness_1;
        float3  _S2043 = mean_14 - - mul_0(transpose_0(R_14), t_13);
        float _S2044 = _S2043.x;
        float _S2045 = _S2043.y;
        float _S2046 = _S2043.z;
        float norm_7 = (F32_sqrt((_S2044 * _S2044 + _S2045 * _S2045 + _S2046 * _S2046)));
        float x_49 = _S2044 / norm_7;
        float y_29 = _S2045 / norm_7;
        float z_26 = _S2046 / norm_7;
        float z2_26 = z_26 * z_26;
        float fTmp0B_10 = -1.09254848957061768f * z_26;
        float fC1_10 = x_49 * x_49 - y_29 * y_29;
        float fS1_10 = 2.0f * x_49 * y_29;
        float fTmp0C_10 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
        float fTmp1B_10 = 1.44530570507049561f * z_26;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_29) * (*sh_coeffs_10)[int(1)] + make_float3 (z_26) * (*sh_coeffs_10)[int(2)] - make_float3 (x_49) * (*sh_coeffs_10)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_10) * (*sh_coeffs_10)[int(4)] + make_float3 (fTmp0B_10 * y_29) * (*sh_coeffs_10)[int(5)] + make_float3 (0.94617468118667603f * z2_26 - 0.31539157032966614f) * (*sh_coeffs_10)[int(6)] + make_float3 (fTmp0B_10 * x_49) * (*sh_coeffs_10)[int(7)] + make_float3 (0.54627424478530884f * fC1_10) * (*sh_coeffs_10)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_49 * fS1_10 + y_29 * fC1_10)) * (*sh_coeffs_10)[int(9)] + make_float3 (fTmp1B_10 * fS1_10) * (*sh_coeffs_10)[int(10)] + make_float3 (fTmp0C_10 * y_29) * (*sh_coeffs_10)[int(11)] + make_float3 (z_26 * (1.86588168144226074f * z2_26 - 1.11952900886535645f)) * (*sh_coeffs_10)[int(12)] + make_float3 (fTmp0C_10 * x_49) * (*sh_coeffs_10)[int(13)] + make_float3 (fTmp1B_10 * fC1_10) * (*sh_coeffs_10)[int(14)] + make_float3 (-0.59004360437393188f * (x_49 * fC1_10 - y_29 * fS1_10)) * (*sh_coeffs_10)[int(15)]);
        float3  _S2047 = make_float3 (0.0f);
        (*rgb_9)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S2047);
        float3  _S2048 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S2049 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_9)[int(1)] = max_0(_S2048 + _S2049 + make_float3 (0.5f), _S2047);
        (*rgb_9)[int(2)] = max_0(_S2048 - _S2049 + make_float3 (0.5f), _S2047);
        float3  _S2050 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S2050 * make_float3 (float(- (F32_sign((dot_0(_S2050, mean_c_10))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_15, float4  quat_18, float3  scale_17, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_11, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_15, float fy_15, float cx_15, float cy_15, float4  radial_coeffs_18, float2  tangential_coeffs_18, float2  thin_prism_coeffs_18, uint image_width_11, uint image_height_11, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_9, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_10, float3  * normal_2)
{
    float3  mean_c_11 = mul_0(R_15, mean_15) + t_14;
    float _S2051 = scale_17.x;
    float sx_2 = (F32_exp((_S2051)));
    float _S2052 = scale_17.y;
    float sy_2 = (F32_exp((_S2052)));
    float sz_2 = scale_17.z - 0.5f * (_S2051 + _S2052);
    float x_50 = quat_18.y;
    float inv_norm_12 = (F32_rsqrt((x_50 * x_50 + quat_18.z * quat_18.z + quat_18.w * quat_18.w + quat_18.x * quat_18.x)));
    float x_51 = quat_18.y * inv_norm_12;
    float y_30 = quat_18.z * inv_norm_12;
    float z_27 = quat_18.w * inv_norm_12;
    float w_16 = quat_18.x * inv_norm_12;
    float x2_16 = x_51 * x_51;
    float y2_16 = y_30 * y_30;
    float z2_27 = z_27 * z_27;
    float xy_16 = x_51 * y_30;
    float xz_16 = x_51 * z_27;
    float yz_16 = y_30 * z_27;
    float wx_16 = w_16 * x_51;
    float wy_16 = w_16 * y_30;
    float wz_16 = w_16 * z_27;
    Matrix<float, 3, 3>  _S2053 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_27), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_27), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    float3  vert0_c_2 = mul_0(R_15, mul_0(_S2053, make_float3 (sx_2, 0.0f, 0.0f)) + mean_15) + t_14;
    float3  vert1_c_2 = mul_0(R_15, mul_0(_S2053, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_15) + t_14;
    float3  vert2_c_2 = mul_0(R_15, mul_0(_S2053, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_15) + t_14;
    *uv0_2 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    *uv1_2 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    *uv2_2 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float2  _S2054 = make_float2 (fx_15, fy_15);
    float2  _S2055 = make_float2 (cx_15, cy_15);
    *uv0_2 = _S2054 * *uv0_2 + _S2055;
    *uv1_2 = _S2054 * *uv1_2 + _S2055;
    float2  _S2056 = _S2054 * *uv2_2 + _S2055;
    *uv2_2 = _S2056;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S2056 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S2056)));
    float _S2057 = _S2056.x;
    float _S2058 = _S2056.y;
    *aabb_xyxy_8 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S2057))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S2058))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S2057))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S2058))) + offset_2)))));
    *depth_9 = make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2));
    *out_hardness_2 = hardness_2;
    float3  _S2059 = mean_15 - - mul_0(transpose_0(R_15), t_14);
    float _S2060 = _S2059.x;
    float _S2061 = _S2059.y;
    float _S2062 = _S2059.z;
    float norm_8 = (F32_sqrt((_S2060 * _S2060 + _S2061 * _S2061 + _S2062 * _S2062)));
    float x_52 = _S2060 / norm_8;
    float y_31 = _S2061 / norm_8;
    float z_28 = _S2062 / norm_8;
    float z2_28 = z_28 * z_28;
    float fTmp0B_11 = -1.09254848957061768f * z_28;
    float fC1_11 = x_52 * x_52 - y_31 * y_31;
    float fS1_11 = 2.0f * x_52 * y_31;
    float fTmp0C_11 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_28;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * (*sh_coeffs_11)[int(1)] + make_float3 (z_28) * (*sh_coeffs_11)[int(2)] - make_float3 (x_52) * (*sh_coeffs_11)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_11) * (*sh_coeffs_11)[int(4)] + make_float3 (fTmp0B_11 * y_31) * (*sh_coeffs_11)[int(5)] + make_float3 (0.94617468118667603f * z2_28 - 0.31539157032966614f) * (*sh_coeffs_11)[int(6)] + make_float3 (fTmp0B_11 * x_52) * (*sh_coeffs_11)[int(7)] + make_float3 (0.54627424478530884f * fC1_11) * (*sh_coeffs_11)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_52 * fS1_11 + y_31 * fC1_11)) * (*sh_coeffs_11)[int(9)] + make_float3 (fTmp1B_11 * fS1_11) * (*sh_coeffs_11)[int(10)] + make_float3 (fTmp0C_11 * y_31) * (*sh_coeffs_11)[int(11)] + make_float3 (z_28 * (1.86588168144226074f * z2_28 - 1.11952900886535645f)) * (*sh_coeffs_11)[int(12)] + make_float3 (fTmp0C_11 * x_52) * (*sh_coeffs_11)[int(13)] + make_float3 (fTmp1B_11 * fC1_11) * (*sh_coeffs_11)[int(14)] + make_float3 (-0.59004360437393188f * (x_52 * fC1_11 - y_31 * fS1_11)) * (*sh_coeffs_11)[int(15)]);
    float3  _S2063 = make_float3 (0.0f);
    (*rgb_10)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S2063);
    float3  _S2064 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S2065 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_10)[int(1)] = max_0(_S2064 + _S2065 + make_float3 (0.5f), _S2063);
    (*rgb_10)[int(2)] = max_0(_S2064 - _S2065 + make_float3 (0.5f), _S2063);
    float3  _S2066 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S2066 * make_float3 (float(- (F32_sign((dot_0(_S2066, mean_c_11))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_16, float4  quat_19, float3  scale_18, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_12, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_16, float fy_16, float cx_16, float cy_16, float4  radial_coeffs_19, float2  tangential_coeffs_19, float2  thin_prism_coeffs_19, uint image_width_12, uint image_height_12, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_10, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_11, float3  * normal_3)
{
    float3  mean_c_12 = mul_0(R_16, mean_16) + t_15;
    float _S2067 = scale_18.x;
    float sx_3 = (F32_exp((_S2067)));
    float _S2068 = scale_18.y;
    float sy_3 = (F32_exp((_S2068)));
    float sz_3 = scale_18.z - 0.5f * (_S2067 + _S2068);
    float x_53 = quat_19.y;
    float inv_norm_13 = (F32_rsqrt((x_53 * x_53 + quat_19.z * quat_19.z + quat_19.w * quat_19.w + quat_19.x * quat_19.x)));
    float x_54 = quat_19.y * inv_norm_13;
    float y_32 = quat_19.z * inv_norm_13;
    float z_29 = quat_19.w * inv_norm_13;
    float w_17 = quat_19.x * inv_norm_13;
    float x2_17 = x_54 * x_54;
    float y2_17 = y_32 * y_32;
    float z2_29 = z_29 * z_29;
    float xy_17 = x_54 * y_32;
    float xz_17 = x_54 * z_29;
    float yz_17 = y_32 * z_29;
    float wx_17 = w_17 * x_54;
    float wy_17 = w_17 * y_32;
    float wz_17 = w_17 * z_29;
    Matrix<float, 3, 3>  _S2069 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_29), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_29), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    float3  vert0_c_3 = mul_0(R_16, mul_0(_S2069, make_float3 (sx_3, 0.0f, 0.0f)) + mean_16) + t_15;
    float3  vert1_c_3 = mul_0(R_16, mul_0(_S2069, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_16) + t_15;
    float3  vert2_c_3 = mul_0(R_16, mul_0(_S2069, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_16) + t_15;
    CameraDistortion_0 dist_coeffs_3 = CameraDistortion_x24init_0(radial_coeffs_19, tangential_coeffs_19, thin_prism_coeffs_19);
    float2  _S2070 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_10 = length_0(_S2070);
    float _S2071 = vert0_c_3.z;
    float theta_5 = (F32_atan2((r_10), (_S2071)));
    float k_5;
    if(theta_5 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_5 * theta_5 / 3.0f) / _S2071;
    }
    else
    {
        k_5 = theta_5 / r_10;
    }
    float2  _S2072 = _S2070 * make_float2 (k_5);
    float k1_4 = dist_coeffs_3.radial_coeffs_0.x;
    float k2_4 = dist_coeffs_3.radial_coeffs_0.y;
    float k3_4 = dist_coeffs_3.radial_coeffs_0.z;
    float k4_5 = dist_coeffs_3.radial_coeffs_0.w;
    float p1_5 = dist_coeffs_3.tangential_coeffs_0.x;
    float p2_5 = dist_coeffs_3.tangential_coeffs_0.y;
    float sx1_5 = dist_coeffs_3.thin_prism_coeffs_0.x;
    float sy1_5 = dist_coeffs_3.thin_prism_coeffs_0.y;
    float u_13 = _S2072.x;
    float v_13 = _S2072.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float _S2073 = 2.0f * p1_5;
    float _S2074 = 2.0f * p2_5;
    float2  _S2075 = _S2072 * make_float2 (1.0f + r2_13 * (k1_4 + r2_13 * (k2_4 + r2_13 * (k3_4 + r2_13 * k4_5)))) + make_float2 (_S2073 * u_13 * v_13 + p2_5 * (r2_13 + 2.0f * u_13 * u_13) + sx1_5 * r2_13, _S2074 * u_13 * v_13 + p1_5 * (r2_13 + 2.0f * v_13 * v_13) + sy1_5 * r2_13);
    *uv0_3 = make_float2 (fx_16 * _S2075.x + cx_16, fy_16 * _S2075.y + cy_16);
    float2  _S2076 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_11 = length_0(_S2076);
    float _S2077 = vert1_c_3.z;
    float theta_6 = (F32_atan2((r_11), (_S2077)));
    if(theta_6 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_6 * theta_6 / 3.0f) / _S2077;
    }
    else
    {
        k_5 = theta_6 / r_11;
    }
    float2  _S2078 = _S2076 * make_float2 (k_5);
    float u_14 = _S2078.x;
    float v_14 = _S2078.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S2079 = _S2078 * make_float2 (1.0f + r2_14 * (k1_4 + r2_14 * (k2_4 + r2_14 * (k3_4 + r2_14 * k4_5)))) + make_float2 (_S2073 * u_14 * v_14 + p2_5 * (r2_14 + 2.0f * u_14 * u_14) + sx1_5 * r2_14, _S2074 * u_14 * v_14 + p1_5 * (r2_14 + 2.0f * v_14 * v_14) + sy1_5 * r2_14);
    *uv1_3 = make_float2 (fx_16 * _S2079.x + cx_16, fy_16 * _S2079.y + cy_16);
    float2  _S2080 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_12 = length_0(_S2080);
    float _S2081 = vert2_c_3.z;
    float theta_7 = (F32_atan2((r_12), (_S2081)));
    if(theta_7 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_7 * theta_7 / 3.0f) / _S2081;
    }
    else
    {
        k_5 = theta_7 / r_12;
    }
    float2  _S2082 = _S2080 * make_float2 (k_5);
    float u_15 = _S2082.x;
    float v_15 = _S2082.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S2083 = _S2082 * make_float2 (1.0f + r2_15 * (k1_4 + r2_15 * (k2_4 + r2_15 * (k3_4 + r2_15 * k4_5)))) + make_float2 (_S2073 * u_15 * v_15 + p2_5 * (r2_15 + 2.0f * u_15 * u_15) + sx1_5 * r2_15, _S2074 * u_15 * v_15 + p1_5 * (r2_15 + 2.0f * v_15 * v_15) + sy1_5 * r2_15);
    float _S2084 = fx_16 * _S2083.x + cx_16;
    float _S2085 = fy_16 * _S2083.y + cy_16;
    float2  _S2086 = make_float2 (_S2084, _S2085);
    *uv2_3 = _S2086;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S2086 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S2086)));
    *aabb_xyxy_9 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S2084))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S2085))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S2084))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S2085))) + offset_3)))));
    *depth_10 = make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3));
    *out_hardness_3 = hardness_3;
    float3  _S2087 = mean_16 - - mul_0(transpose_0(R_16), t_15);
    float _S2088 = _S2087.x;
    float _S2089 = _S2087.y;
    float _S2090 = _S2087.z;
    float norm_9 = (F32_sqrt((_S2088 * _S2088 + _S2089 * _S2089 + _S2090 * _S2090)));
    float x_55 = _S2088 / norm_9;
    float y_33 = _S2089 / norm_9;
    float z_30 = _S2090 / norm_9;
    float z2_30 = z_30 * z_30;
    float fTmp0B_12 = -1.09254848957061768f * z_30;
    float fC1_12 = x_55 * x_55 - y_33 * y_33;
    float fS1_12 = 2.0f * x_55 * y_33;
    float fTmp0C_12 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_30;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_33) * (*sh_coeffs_12)[int(1)] + make_float3 (z_30) * (*sh_coeffs_12)[int(2)] - make_float3 (x_55) * (*sh_coeffs_12)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_12) * (*sh_coeffs_12)[int(4)] + make_float3 (fTmp0B_12 * y_33) * (*sh_coeffs_12)[int(5)] + make_float3 (0.94617468118667603f * z2_30 - 0.31539157032966614f) * (*sh_coeffs_12)[int(6)] + make_float3 (fTmp0B_12 * x_55) * (*sh_coeffs_12)[int(7)] + make_float3 (0.54627424478530884f * fC1_12) * (*sh_coeffs_12)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_55 * fS1_12 + y_33 * fC1_12)) * (*sh_coeffs_12)[int(9)] + make_float3 (fTmp1B_12 * fS1_12) * (*sh_coeffs_12)[int(10)] + make_float3 (fTmp0C_12 * y_33) * (*sh_coeffs_12)[int(11)] + make_float3 (z_30 * (1.86588168144226074f * z2_30 - 1.11952900886535645f)) * (*sh_coeffs_12)[int(12)] + make_float3 (fTmp0C_12 * x_55) * (*sh_coeffs_12)[int(13)] + make_float3 (fTmp1B_12 * fC1_12) * (*sh_coeffs_12)[int(14)] + make_float3 (-0.59004360437393188f * (x_55 * fC1_12 - y_33 * fS1_12)) * (*sh_coeffs_12)[int(15)]);
    float3  _S2091 = make_float3 (0.0f);
    (*rgb_11)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S2091);
    float3  _S2092 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S2093 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_11)[int(1)] = max_0(_S2092 + _S2093 + make_float3 (0.5f), _S2091);
    (*rgb_11)[int(2)] = max_0(_S2092 - _S2093 + make_float3 (0.5f), _S2091);
    float3  _S2094 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S2094 * make_float3 (float(- (F32_sign((dot_0(_S2094, mean_c_12))))));
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S2095, float _S2096)
{
    _d_exp2_0(_S2095, _S2096);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S2097, float _S2098)
{
    _d_abs_0(_S2097, _S2098);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_17, float4  quat_20, float3  scale_19, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_13, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_17, float fy_17, float cx_17, float cy_17, float4  radial_coeffs_20, float2  tangential_coeffs_20, float2  thin_prism_coeffs_20, uint image_width_13, uint image_height_13, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_4, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_4, float3  v_normal_0, float3  * v_mean_5, float4  * v_quat_5, float3  * v_scale_5, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_3, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_17) + t_16;
    float _S2099 = scale_19.x;
    float _S2100 = s_primal_ctx_exp_1(_S2099);
    float _S2101 = scale_19.y;
    float _S2102 = s_primal_ctx_exp_1(_S2101);
    float sz_4 = scale_19.z - 0.5f * (_S2099 + _S2101);
    float _S2103 = quat_20.y;
    float _S2104 = _S2103 * _S2103 + quat_20.z * quat_20.z + quat_20.w * quat_20.w + quat_20.x * quat_20.x;
    float _S2105 = s_primal_ctx_rsqrt_0(_S2104);
    float x_56 = quat_20.y * _S2105;
    float y_34 = quat_20.z * _S2105;
    float z_31 = quat_20.w * _S2105;
    float w_18 = quat_20.x * _S2105;
    float x2_18 = x_56 * x_56;
    float y2_18 = y_34 * y_34;
    float z2_31 = z_31 * z_31;
    float xy_18 = x_56 * y_34;
    float xz_18 = x_56 * z_31;
    float yz_18 = y_34 * z_31;
    float wx_18 = w_18 * x_56;
    float wy_18 = w_18 * y_34;
    float wz_18 = w_18 * z_31;
    Matrix<float, 3, 3>  _S2106 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_18 + z2_31), 2.0f * (xy_18 + wz_18), 2.0f * (xz_18 - wy_18), 2.0f * (xy_18 - wz_18), 1.0f - 2.0f * (x2_18 + z2_31), 2.0f * (yz_18 + wx_18), 2.0f * (xz_18 + wy_18), 2.0f * (yz_18 - wx_18), 1.0f - 2.0f * (x2_18 + y2_18)));
    float3  _S2107 = make_float3 (_S2100, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_1(_S2106, _S2107) + mean_17;
    float _S2108 = -0.5f + sz_4;
    float3  _S2109 = make_float3 (_S2100 * _S2108, _S2102, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_1(_S2106, _S2109) + mean_17;
    float _S2110 = -0.5f - sz_4;
    float3  _S2111 = make_float3 (_S2100 * _S2110, - _S2102, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_1(_S2106, _S2111) + mean_17;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_17, vert0_0) + t_16;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_17, vert1_0) + t_16;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_17, vert2_0) + t_16;
    float2  _S2112 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S2113 = vert0_c_4.z;
    float2  _S2114 = make_float2 (_S2113);
    float2  _S2115 = make_float2 (_S2113 * _S2113);
    float2  _S2116 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S2117 = vert1_c_4.z;
    float2  _S2118 = make_float2 (_S2117);
    float2  _S2119 = make_float2 (_S2117 * _S2117);
    float2  _S2120 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S2121 = vert2_c_4.z;
    float2  _S2122 = make_float2 (_S2121);
    float2  _S2123 = make_float2 (_S2121 * _S2121);
    float2  _S2124 = make_float2 (fx_17, fy_17);
    float2  _S2125 = make_float2 (cx_17, cy_17);
    float2  _S2126 = _S2124 * (_S2112 / make_float2 (_S2113)) + _S2125;
    float2  _S2127 = _S2124 * (_S2116 / make_float2 (_S2117)) + _S2125;
    float2  _S2128 = _S2124 * (_S2120 / make_float2 (_S2121)) + _S2125;
    float2  e0_4 = _S2127 - _S2126;
    float2  e1_4 = _S2128 - _S2127;
    float2  e2_0 = _S2126 - _S2128;
    float _S2129 = e0_4.x;
    float _S2130 = e1_4.y;
    float _S2131 = e0_4.y;
    float _S2132 = e1_4.x;
    float _S2133 = _S2129 * _S2130 - _S2131 * _S2132;
    float _S2134 = 1.0f - hardness_4.y;
    float _S2135 = -1.0f / _S2134;
    float _S2136 = _S2134 * _S2134;
    float _S2137 = _S2126.x;
    float _S2138 = _S2127.x;
    float _S2139 = s_primal_ctx_max_0(_S2137, _S2138);
    float _S2140 = _S2128.x;
    float _S2141 = s_primal_ctx_min_0(_S2137, _S2138);
    float _S2142 = _S2126.y;
    float _S2143 = _S2127.y;
    float _S2144 = s_primal_ctx_max_0(_S2142, _S2143);
    float _S2145 = _S2128.y;
    float _S2146 = s_primal_ctx_min_0(_S2142, _S2143);
    Matrix<float, 3, 3>  _S2147 = transpose_0(R_17);
    float3  _S2148 = mean_17 - - s_primal_ctx_mul_1(_S2147, t_16);
    float _S2149 = _S2148.x;
    float _S2150 = _S2148.y;
    float _S2151 = _S2148.z;
    float _S2152 = _S2149 * _S2149 + _S2150 * _S2150 + _S2151 * _S2151;
    float _S2153 = s_primal_ctx_sqrt_0(_S2152);
    float x_57 = _S2149 / _S2153;
    float3  _S2154 = make_float3 (x_57);
    float _S2155 = _S2153 * _S2153;
    float y_35 = _S2150 / _S2153;
    float z_32 = _S2151 / _S2153;
    float3  _S2156 = make_float3 (z_32);
    float _S2157 = - y_35;
    float3  _S2158 = make_float3 (_S2157);
    float z2_32 = z_32 * z_32;
    float fTmp0B_13 = -1.09254848957061768f * z_32;
    float fC1_13 = x_57 * x_57 - y_35 * y_35;
    float _S2159 = 2.0f * x_57;
    float fS1_13 = _S2159 * y_35;
    float pSH6_3 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2160 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_57;
    float3  _S2161 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_35;
    float3  _S2162 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S2163 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S2164 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_32;
    float _S2165 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_3 = z_32 * _S2165;
    float3  _S2166 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_57;
    float3  _S2167 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_35;
    float3  _S2168 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S2169 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S2170 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_57 * fC1_13 - y_35 * fS1_13);
    float3  _S2171 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_57 * fS1_13 + y_35 * fC1_13);
    float3  _S2172 = make_float3 (pSH9_3);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2157) * (*sh_coeffs_13)[int(1)] + make_float3 (z_32) * (*sh_coeffs_13)[int(2)] - make_float3 (x_57) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]);
    float3  _S2173 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S2174 = make_float3 (0.0f);
    float3  _S2175 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S2176 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S2177 = make_float3 (_S2176);
    float3  _S2178 = (*ch_coeffs_4)[int(1)] * make_float3 (_S2176);
    float3  _S2179 = _S2175 + _S2178 + make_float3 (0.5f);
    float3  _S2180 = _S2175 - _S2178 + make_float3 (0.5f);
    float3  _S2181 = vert1_c_4 - vert0_c_4;
    float3  _S2182 = vert2_c_4 - vert0_c_4;
    float3  _S2183 = s_primal_ctx_cross_0(_S2181, _S2182);
    float3  _S2184 = normalize_0(_S2183);
    float3  _S2185 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S2184, mean_c_13)))))) * v_normal_0;
    float3  _S2186 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2187;
    (&_S2187)->primal_0 = _S2184;
    (&_S2187)->differential_0 = _S2186;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2188;
    (&_S2188)->primal_0 = mean_c_13;
    (&_S2188)->differential_0 = _S2186;
    s_bwd_prop_dot_0(&_S2187, &_S2188, 0.0f);
    float3  _S2189 = _S2185 + _S2187.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2190;
    (&_S2190)->primal_0 = _S2183;
    (&_S2190)->differential_0 = _S2186;
    s_bwd_normalize_impl_0(&_S2190, _S2189);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2191;
    (&_S2191)->primal_0 = _S2181;
    (&_S2191)->differential_0 = _S2186;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2192;
    (&_S2192)->primal_0 = _S2182;
    (&_S2192)->differential_0 = _S2186;
    s_bwd_prop_cross_0(&_S2191, &_S2192, _S2190.differential_0);
    float3  _S2193 = - _S2192.differential_0;
    float3  _S2194 = - _S2191.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2195;
    (&_S2195)->primal_0 = _S2180;
    (&_S2195)->differential_0 = _S2186;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2196;
    (&_S2196)->primal_0 = _S2174;
    (&_S2196)->differential_0 = _S2186;
    s_bwd_prop_max_0(&_S2195, &_S2196, (*v_rgb_4)[int(2)]);
    float3  _S2197 = - _S2195.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2198;
    (&_S2198)->primal_0 = _S2179;
    (&_S2198)->differential_0 = _S2186;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2199;
    (&_S2199)->primal_0 = _S2174;
    (&_S2199)->differential_0 = _S2186;
    s_bwd_prop_max_0(&_S2198, &_S2199, (*v_rgb_4)[int(1)]);
    float3  _S2200 = _S2177 * (_S2197 + _S2198.differential_0);
    float3  _S2201 = _S2195.differential_0 + _S2198.differential_0;
    float3  _S2202 = make_float3 (0.5f) * - _S2201;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2203;
    (&_S2203)->primal_0 = _S2173;
    (&_S2203)->differential_0 = _S2186;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2204;
    (&_S2204)->primal_0 = _S2174;
    (&_S2204)->differential_0 = _S2186;
    s_bwd_prop_max_0(&_S2203, &_S2204, (*v_rgb_4)[int(0)]);
    float3  _S2205 = _S2202 + _S2203.differential_0;
    float3  _S2206 = _S2201 + _S2203.differential_0;
    float3  _S2207 = _S2171 * _S2206;
    float3  _S2208 = (*sh_coeffs_13)[int(15)] * _S2206;
    float3  _S2209 = _S2169 * _S2206;
    float3  _S2210 = (*sh_coeffs_13)[int(14)] * _S2206;
    float3  _S2211 = _S2167 * _S2206;
    float3  _S2212 = (*sh_coeffs_13)[int(13)] * _S2206;
    float3  _S2213 = _S2166 * _S2206;
    float3  _S2214 = (*sh_coeffs_13)[int(12)] * _S2206;
    float3  _S2215 = _S2168 * _S2206;
    float3  _S2216 = (*sh_coeffs_13)[int(11)] * _S2206;
    float3  _S2217 = _S2170 * _S2206;
    float3  _S2218 = (*sh_coeffs_13)[int(10)] * _S2206;
    float3  _S2219 = _S2172 * _S2206;
    float3  _S2220 = (*sh_coeffs_13)[int(9)] * _S2206;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2220.x + _S2220.y + _S2220.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S2208.x + _S2208.y + _S2208.z);
    float _S2221 = _S2218.x + _S2218.y + _S2218.z;
    float _S2222 = _S2210.x + _S2210.y + _S2210.z;
    float _S2223 = _S2216.x + _S2216.y + _S2216.z;
    float _S2224 = _S2212.x + _S2212.y + _S2212.z;
    float _S2225 = _S2214.x + _S2214.y + _S2214.z;
    float _S2226 = - s_diff_fC2_T_3;
    float3  _S2227 = _S2163 * _S2206;
    float3  _S2228 = (*sh_coeffs_13)[int(8)] * _S2206;
    float3  _S2229 = _S2161 * _S2206;
    float3  _S2230 = (*sh_coeffs_13)[int(7)] * _S2206;
    float3  _S2231 = _S2160 * _S2206;
    float3  _S2232 = (*sh_coeffs_13)[int(6)] * _S2206;
    float3  _S2233 = _S2162 * _S2206;
    float3  _S2234 = (*sh_coeffs_13)[int(5)] * _S2206;
    float3  _S2235 = _S2164 * _S2206;
    float3  _S2236 = (*sh_coeffs_13)[int(4)] * _S2206;
    float _S2237 = _S2234.x + _S2234.y + _S2234.z;
    float _S2238 = _S2230.x + _S2230.y + _S2230.z;
    float _S2239 = fTmp1B_13 * _S2221 + x_57 * s_diff_fS2_T_3 + y_35 * _S2226 + 0.54627424478530884f * (_S2236.x + _S2236.y + _S2236.z);
    float _S2240 = fTmp1B_13 * _S2222 + y_35 * s_diff_fS2_T_3 + x_57 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2228.x + _S2228.y + _S2228.z);
    float _S2241 = y_35 * - _S2240;
    float _S2242 = x_57 * _S2240;
    float _S2243 = z_32 * (1.86588168144226074f * (z_32 * _S2225) + -2.28522896766662598f * (y_35 * _S2223 + x_57 * _S2224) + 0.94617468118667603f * (_S2232.x + _S2232.y + _S2232.z));
    float3  _S2244 = make_float3 (0.48860251903533936f) * _S2206;
    float3  _S2245 = - _S2244;
    float3  _S2246 = _S2154 * _S2245;
    float3  _S2247 = (*sh_coeffs_13)[int(3)] * _S2245;
    float3  _S2248 = _S2156 * _S2244;
    float3  _S2249 = (*sh_coeffs_13)[int(2)] * _S2244;
    float3  _S2250 = _S2158 * _S2244;
    float3  _S2251 = (*sh_coeffs_13)[int(1)] * _S2244;
    float _S2252 = (_S2165 * _S2225 + 1.44530570507049561f * (fS1_13 * _S2221 + fC1_13 * _S2222) + -1.09254848957061768f * (y_35 * _S2237 + x_57 * _S2238) + _S2243 + _S2243 + _S2249.x + _S2249.y + _S2249.z) / _S2155;
    float _S2253 = _S2153 * _S2252;
    float _S2254 = (fTmp0C_13 * _S2223 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2226 + fTmp0B_13 * _S2237 + _S2159 * _S2239 + _S2241 + _S2241 + - (_S2251.x + _S2251.y + _S2251.z)) / _S2155;
    float _S2255 = _S2153 * _S2254;
    float _S2256 = (fTmp0C_13 * _S2224 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2238 + 2.0f * (y_35 * _S2239) + _S2242 + _S2242 + _S2247.x + _S2247.y + _S2247.z) / _S2155;
    float _S2257 = _S2153 * _S2256;
    float _S2258 = _S2151 * - _S2252 + _S2150 * - _S2254 + _S2149 * - _S2256;
    DiffPair_float_0 _S2259;
    (&_S2259)->primal_0 = _S2152;
    (&_S2259)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2259, _S2258);
    float _S2260 = _S2151 * _S2259.differential_0;
    float _S2261 = _S2150 * _S2259.differential_0;
    float _S2262 = _S2149 * _S2259.differential_0;
    float3  _S2263 = make_float3 (0.282094806432724f) * _S2206;
    float3  _S2264 = make_float3 (_S2257 + _S2262 + _S2262, _S2255 + _S2261 + _S2261, _S2253 + _S2260 + _S2260);
    float3  _S2265 = - - _S2264;
    Matrix<float, 3, 3>  _S2266 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2267;
    (&_S2267)->primal_0 = _S2147;
    (&_S2267)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2268;
    (&_S2268)->primal_0 = t_16;
    (&_S2268)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2267, &_S2268, _S2265);
    Matrix<float, 3, 3>  _S2269 = transpose_0(_S2267.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2270;
    (&_S2270)->primal_0 = vert2_c_4;
    (&_S2270)->differential_0 = _S2186;
    s_bwd_length_impl_0(&_S2270, v_depth_4.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2271;
    (&_S2271)->primal_0 = vert1_c_4;
    (&_S2271)->differential_0 = _S2186;
    s_bwd_length_impl_0(&_S2271, v_depth_4.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2272;
    (&_S2272)->primal_0 = vert0_c_4;
    (&_S2272)->differential_0 = _S2186;
    s_bwd_length_impl_0(&_S2272, v_depth_4.x);
    DiffPair_float_0 _S2273;
    (&_S2273)->primal_0 = _S2146;
    (&_S2273)->differential_0 = 0.0f;
    DiffPair_float_0 _S2274;
    (&_S2274)->primal_0 = _S2145;
    (&_S2274)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2273, &_S2274, 0.0f);
    DiffPair_float_0 _S2275;
    (&_S2275)->primal_0 = _S2142;
    (&_S2275)->differential_0 = 0.0f;
    DiffPair_float_0 _S2276;
    (&_S2276)->primal_0 = _S2143;
    (&_S2276)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2275, &_S2276, _S2273.differential_0);
    DiffPair_float_0 _S2277;
    (&_S2277)->primal_0 = _S2144;
    (&_S2277)->differential_0 = 0.0f;
    DiffPair_float_0 _S2278;
    (&_S2278)->primal_0 = _S2145;
    (&_S2278)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2277, &_S2278, 0.0f);
    float _S2279 = _S2274.differential_0 + _S2278.differential_0;
    DiffPair_float_0 _S2280;
    (&_S2280)->primal_0 = _S2142;
    (&_S2280)->differential_0 = 0.0f;
    DiffPair_float_0 _S2281;
    (&_S2281)->primal_0 = _S2143;
    (&_S2281)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2280, &_S2281, _S2277.differential_0);
    float _S2282 = _S2276.differential_0 + _S2281.differential_0;
    float _S2283 = _S2275.differential_0 + _S2280.differential_0;
    DiffPair_float_0 _S2284;
    (&_S2284)->primal_0 = _S2141;
    (&_S2284)->differential_0 = 0.0f;
    DiffPair_float_0 _S2285;
    (&_S2285)->primal_0 = _S2140;
    (&_S2285)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2284, &_S2285, 0.0f);
    DiffPair_float_0 _S2286;
    (&_S2286)->primal_0 = _S2137;
    (&_S2286)->differential_0 = 0.0f;
    DiffPair_float_0 _S2287;
    (&_S2287)->primal_0 = _S2138;
    (&_S2287)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2286, &_S2287, _S2284.differential_0);
    DiffPair_float_0 _S2288;
    (&_S2288)->primal_0 = _S2139;
    (&_S2288)->differential_0 = 0.0f;
    DiffPair_float_0 _S2289;
    (&_S2289)->primal_0 = _S2140;
    (&_S2289)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2288, &_S2289, 0.0f);
    float _S2290 = _S2285.differential_0 + _S2289.differential_0;
    DiffPair_float_0 _S2291;
    (&_S2291)->primal_0 = _S2137;
    (&_S2291)->differential_0 = 0.0f;
    DiffPair_float_0 _S2292;
    (&_S2292)->primal_0 = _S2138;
    (&_S2292)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2291, &_S2292, _S2288.differential_0);
    float _S2293 = _S2287.differential_0 + _S2292.differential_0;
    float _S2294 = _S2286.differential_0 + _S2291.differential_0;
    DiffPair_float_0 _S2295;
    (&_S2295)->primal_0 = _S2135;
    (&_S2295)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2295, 0.0f);
    float _S2296 = - (-1.0f * - (_S2295.differential_0 / _S2136));
    float2  _S2297 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2298;
    (&_S2298)->primal_0 = e2_0;
    (&_S2298)->differential_0 = _S2297;
    s_bwd_length_impl_1(&_S2298, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2299;
    (&_S2299)->primal_0 = e1_4;
    (&_S2299)->differential_0 = _S2297;
    s_bwd_length_impl_1(&_S2299, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2300;
    (&_S2300)->primal_0 = e0_4;
    (&_S2300)->differential_0 = _S2297;
    s_bwd_length_impl_1(&_S2300, -0.0f);
    DiffPair_float_0 _S2301;
    (&_S2301)->primal_0 = _S2133;
    (&_S2301)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2301, 0.0f);
    float _S2302 = - _S2301.differential_0;
    float2  _S2303 = _S2299.differential_0 + make_float2 (_S2131 * _S2302, _S2129 * _S2301.differential_0);
    float2  _S2304 = _S2300.differential_0 + make_float2 (_S2130 * _S2301.differential_0, _S2132 * _S2302);
    float2  _S2305 = _S2124 * (v_uv2_0 + - _S2298.differential_0 + _S2303 + make_float2 (_S2290, _S2279)) / _S2123;
    float2  _S2306 = _S2120 * - _S2305;
    float2  _S2307 = _S2122 * _S2305;
    float2  _S2308 = _S2124 * (v_uv1_0 + - _S2303 + _S2304 + make_float2 (_S2293, _S2282)) / _S2119;
    float2  _S2309 = _S2116 * - _S2308;
    float2  _S2310 = _S2118 * _S2308;
    float _S2311 = _S2309.x + _S2309.y;
    float2  _S2312 = _S2124 * (v_uv0_0 + _S2298.differential_0 + - _S2304 + make_float2 (_S2294, _S2283)) / _S2115;
    float2  _S2313 = _S2112 * - _S2312;
    float2  _S2314 = _S2114 * _S2312;
    float _S2315 = _S2313.x + _S2313.y;
    float3  _S2316 = _S2192.differential_0 + _S2270.differential_0 + make_float3 (_S2307.x, _S2307.y, _S2306.x + _S2306.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2317;
    (&_S2317)->primal_0 = R_17;
    (&_S2317)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2318;
    (&_S2318)->primal_0 = vert2_0;
    (&_S2318)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2317, &_S2318, _S2316);
    float3  _S2319 = _S2191.differential_0 + _S2271.differential_0 + make_float3 (_S2310.x, _S2310.y, _S2311);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2320;
    (&_S2320)->primal_0 = R_17;
    (&_S2320)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2321;
    (&_S2321)->primal_0 = vert1_0;
    (&_S2321)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2320, &_S2321, _S2319);
    float3  _S2322 = _S2193 + _S2194 + _S2272.differential_0 + make_float3 (_S2314.x, _S2314.y, _S2315);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2323;
    (&_S2323)->primal_0 = R_17;
    (&_S2323)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2324;
    (&_S2324)->primal_0 = vert0_0;
    (&_S2324)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2323, &_S2324, _S2322);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2325;
    (&_S2325)->primal_0 = _S2106;
    (&_S2325)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2326;
    (&_S2326)->primal_0 = _S2111;
    (&_S2326)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2325, &_S2326, _S2318.differential_0);
    float _S2327 = - _S2326.differential_0.y;
    float _S2328 = _S2110 * _S2326.differential_0.x;
    float _S2329 = - (_S2100 * _S2326.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2330;
    (&_S2330)->primal_0 = _S2106;
    (&_S2330)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2331;
    (&_S2331)->primal_0 = _S2109;
    (&_S2331)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2330, &_S2331, _S2321.differential_0);
    float _S2332 = _S2100 * _S2331.differential_0.x;
    float _S2333 = _S2108 * _S2331.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2334;
    (&_S2334)->primal_0 = _S2106;
    (&_S2334)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2335;
    (&_S2335)->primal_0 = _S2107;
    (&_S2335)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2334, &_S2335, _S2324.differential_0);
    Matrix<float, 3, 3>  _S2336 = transpose_0(_S2325.differential_0 + _S2330.differential_0 + _S2334.differential_0);
    float _S2337 = 2.0f * - _S2336.rows[int(2)].z;
    float _S2338 = 2.0f * _S2336.rows[int(2)].y;
    float _S2339 = 2.0f * _S2336.rows[int(2)].x;
    float _S2340 = 2.0f * _S2336.rows[int(1)].z;
    float _S2341 = 2.0f * - _S2336.rows[int(1)].y;
    float _S2342 = 2.0f * _S2336.rows[int(1)].x;
    float _S2343 = 2.0f * _S2336.rows[int(0)].z;
    float _S2344 = 2.0f * _S2336.rows[int(0)].y;
    float _S2345 = 2.0f * - _S2336.rows[int(0)].x;
    float _S2346 = - _S2342 + _S2344;
    float _S2347 = _S2339 + - _S2343;
    float _S2348 = - _S2338 + _S2340;
    float _S2349 = _S2338 + _S2340;
    float _S2350 = _S2339 + _S2343;
    float _S2351 = _S2342 + _S2344;
    float _S2352 = z_31 * (_S2341 + _S2345);
    float _S2353 = y_34 * (_S2337 + _S2345);
    float _S2354 = x_56 * (_S2337 + _S2341);
    float _S2355 = z_31 * _S2346 + y_34 * _S2347 + x_56 * _S2348;
    float _S2356 = _S2105 * _S2355;
    float _S2357 = w_18 * _S2346 + y_34 * _S2349 + x_56 * _S2350 + _S2352 + _S2352;
    float _S2358 = _S2105 * _S2357;
    float _S2359 = w_18 * _S2347 + z_31 * _S2349 + x_56 * _S2351 + _S2353 + _S2353;
    float _S2360 = _S2105 * _S2359;
    float _S2361 = w_18 * _S2348 + z_31 * _S2350 + y_34 * _S2351 + _S2354 + _S2354;
    float _S2362 = _S2105 * _S2361;
    float _S2363 = quat_20.x * _S2355 + quat_20.w * _S2357 + quat_20.z * _S2359 + quat_20.y * _S2361;
    DiffPair_float_0 _S2364;
    (&_S2364)->primal_0 = _S2104;
    (&_S2364)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S2364, _S2363);
    float _S2365 = quat_20.x * _S2364.differential_0;
    float _S2366 = quat_20.w * _S2364.differential_0;
    float _S2367 = quat_20.z * _S2364.differential_0;
    float _S2368 = quat_20.y * _S2364.differential_0;
    float _S2369 = _S2358 + _S2366 + _S2366;
    float _S2370 = _S2360 + _S2367 + _S2367;
    float _S2371 = _S2362 + _S2368 + _S2368;
    float _S2372 = _S2356 + _S2365 + _S2365;
    float _S2373 = _S2329 + _S2332;
    float _S2374 = 0.5f * - _S2373;
    float _S2375 = _S2327 + _S2331.differential_0.y;
    DiffPair_float_0 _S2376;
    (&_S2376)->primal_0 = _S2101;
    (&_S2376)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2376, _S2375);
    float _S2377 = _S2374 + _S2376.differential_0;
    float _S2378 = _S2328 + _S2333 + _S2335.differential_0.x;
    DiffPair_float_0 _S2379;
    (&_S2379)->primal_0 = _S2099;
    (&_S2379)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2379, _S2378);
    float _S2380 = _S2374 + _S2379.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2381;
    (&_S2381)->primal_0 = R_17;
    (&_S2381)->differential_0 = _S2266;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2382;
    (&_S2382)->primal_0 = mean_17;
    (&_S2382)->differential_0 = _S2186;
    s_bwd_prop_mul_1(&_S2381, &_S2382, _S2188.differential_0);
    float3  _S2383 = _S2268.differential_0 + _S2316 + _S2319 + _S2322 + _S2188.differential_0;
    Matrix<float, 3, 3>  _S2384 = _S2269 + _S2317.differential_0 + _S2320.differential_0 + _S2323.differential_0 + _S2381.differential_0;
    FixedArray<float3 , 2>  _S2385;
    _S2385[int(0)] = _S2186;
    _S2385[int(1)] = _S2186;
    _S2385[int(1)] = _S2200;
    _S2385[int(0)] = _S2205;
    FixedArray<float3 , 16>  _S2386;
    _S2386[int(0)] = _S2186;
    _S2386[int(1)] = _S2186;
    _S2386[int(2)] = _S2186;
    _S2386[int(3)] = _S2186;
    _S2386[int(4)] = _S2186;
    _S2386[int(5)] = _S2186;
    _S2386[int(6)] = _S2186;
    _S2386[int(7)] = _S2186;
    _S2386[int(8)] = _S2186;
    _S2386[int(9)] = _S2186;
    _S2386[int(10)] = _S2186;
    _S2386[int(11)] = _S2186;
    _S2386[int(12)] = _S2186;
    _S2386[int(13)] = _S2186;
    _S2386[int(14)] = _S2186;
    _S2386[int(15)] = _S2186;
    _S2386[int(15)] = _S2207;
    _S2386[int(14)] = _S2209;
    _S2386[int(13)] = _S2211;
    _S2386[int(12)] = _S2213;
    _S2386[int(11)] = _S2215;
    _S2386[int(10)] = _S2217;
    _S2386[int(9)] = _S2219;
    _S2386[int(8)] = _S2227;
    _S2386[int(7)] = _S2229;
    _S2386[int(6)] = _S2231;
    _S2386[int(5)] = _S2233;
    _S2386[int(4)] = _S2235;
    _S2386[int(3)] = _S2246;
    _S2386[int(2)] = _S2248;
    _S2386[int(1)] = _S2250;
    _S2386[int(0)] = _S2263;
    float2  _S2387 = v_out_hardness_0 + make_float2 (0.0f, _S2296);
    float3  _S2388 = make_float3 (_S2380, _S2377, _S2373);
    float4  _S2389 = make_float4 (0.0f);
    *&((&_S2389)->w) = _S2369;
    *&((&_S2389)->z) = _S2370;
    *&((&_S2389)->y) = _S2371;
    *&((&_S2389)->x) = _S2372;
    *v_mean_5 = _S2264 + _S2318.differential_0 + _S2321.differential_0 + _S2324.differential_0 + _S2382.differential_0;
    *v_quat_5 = _S2389;
    *v_scale_5 = _S2388;
    *v_hardness_0 = _S2387;
    *v_sh_coeffs_3 = _S2386;
    *v_ch_coeffs_0 = _S2385;
    *v_R_4 = _S2384;
    *v_t_4 = _S2383;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_18, float4  quat_21, float3  scale_20, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_14, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_18, float fy_18, float cx_18, float cy_18, float4  radial_coeffs_21, float2  tangential_coeffs_21, float2  thin_prism_coeffs_21, uint image_width_14, uint image_height_14, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_5, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_5, float3  v_normal_1, float3  * v_mean_6, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_4, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_18) + t_17;
    float _S2390 = scale_20.x;
    float _S2391 = s_primal_ctx_exp_1(_S2390);
    float _S2392 = scale_20.y;
    float _S2393 = s_primal_ctx_exp_1(_S2392);
    float sz_5 = scale_20.z - 0.5f * (_S2390 + _S2392);
    float _S2394 = quat_21.y;
    float _S2395 = _S2394 * _S2394 + quat_21.z * quat_21.z + quat_21.w * quat_21.w + quat_21.x * quat_21.x;
    float _S2396 = s_primal_ctx_rsqrt_0(_S2395);
    float x_58 = quat_21.y * _S2396;
    float y_36 = quat_21.z * _S2396;
    float z_33 = quat_21.w * _S2396;
    float w_19 = quat_21.x * _S2396;
    float x2_19 = x_58 * x_58;
    float y2_19 = y_36 * y_36;
    float z2_33 = z_33 * z_33;
    float xy_19 = x_58 * y_36;
    float xz_19 = x_58 * z_33;
    float yz_19 = y_36 * z_33;
    float wx_19 = w_19 * x_58;
    float wy_19 = w_19 * y_36;
    float wz_19 = w_19 * z_33;
    Matrix<float, 3, 3>  _S2397 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_33), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_33), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19)));
    float3  _S2398 = make_float3 (_S2391, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S2397, _S2398) + mean_18;
    float _S2399 = -0.5f + sz_5;
    float3  _S2400 = make_float3 (_S2391 * _S2399, _S2393, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S2397, _S2400) + mean_18;
    float _S2401 = -0.5f - sz_5;
    float3  _S2402 = make_float3 (_S2391 * _S2401, - _S2393, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S2397, _S2402) + mean_18;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_18, vert0_1) + t_17;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_18, vert1_1) + t_17;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_18, vert2_1) + t_17;
    CameraDistortion_0 _S2403 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_21, tangential_coeffs_21, thin_prism_coeffs_21);
    float2  _S2404 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S2405 = length_0(_S2404);
    float _S2406 = vert0_c_5.z;
    float _S2407 = s_primal_ctx_atan2_0(_S2405, _S2406);
    bool _S2408 = _S2407 < 0.00100000004749745f;
    float k_6;
    float _S2409;
    float _S2410;
    float _S2411;
    if(_S2408)
    {
        float _S2412 = 1.0f - _S2407 * _S2407 / 3.0f;
        float _S2413 = _S2406 * _S2406;
        k_6 = _S2412 / _S2406;
        _S2409 = 0.0f;
        _S2410 = _S2413;
        _S2411 = _S2412;
    }
    else
    {
        float _S2414 = _S2405 * _S2405;
        k_6 = _S2407 / _S2405;
        _S2409 = _S2414;
        _S2410 = 0.0f;
        _S2411 = 0.0f;
    }
    float2  _S2415 = make_float2 (k_6);
    float2  _S2416 = _S2404 * make_float2 (k_6);
    float k1_5 = _S2403.radial_coeffs_0.x;
    float k2_5 = _S2403.radial_coeffs_0.y;
    float k3_5 = _S2403.radial_coeffs_0.z;
    float k4_6 = _S2403.radial_coeffs_0.w;
    float p1_6 = _S2403.tangential_coeffs_0.x;
    float p2_6 = _S2403.tangential_coeffs_0.y;
    float sx1_6 = _S2403.thin_prism_coeffs_0.x;
    float sy1_6 = _S2403.thin_prism_coeffs_0.y;
    float u_16 = _S2416.x;
    float v_16 = _S2416.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S2417 = k3_5 + r2_16 * k4_6;
    float _S2418 = k2_5 + r2_16 * _S2417;
    float _S2419 = k1_5 + r2_16 * _S2418;
    float radial_2 = 1.0f + r2_16 * _S2419;
    float2  _S2420 = make_float2 (radial_2);
    float _S2421 = 2.0f * p1_6;
    float _S2422 = _S2421 * u_16;
    float _S2423 = 2.0f * u_16;
    float _S2424 = r2_16 + _S2423 * u_16;
    float _S2425 = 2.0f * p2_6;
    float _S2426 = _S2425 * u_16;
    float _S2427 = 2.0f * v_16;
    float _S2428 = r2_16 + _S2427 * v_16;
    float2  _S2429 = _S2416 * make_float2 (radial_2) + make_float2 (_S2422 * v_16 + p2_6 * _S2424 + sx1_6 * r2_16, _S2426 * v_16 + p1_6 * _S2428 + sy1_6 * r2_16);
    float _S2430 = fx_18 * _S2429.x + cx_18;
    float _S2431 = fy_18 * _S2429.y + cy_18;
    float2  _S2432 = make_float2 (_S2430, _S2431);
    float2  _S2433 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S2434 = length_0(_S2433);
    float _S2435 = vert1_c_5.z;
    float _S2436 = s_primal_ctx_atan2_0(_S2434, _S2435);
    bool _S2437 = _S2436 < 0.00100000004749745f;
    float _S2438;
    float _S2439;
    float _S2440;
    if(_S2437)
    {
        float _S2441 = 1.0f - _S2436 * _S2436 / 3.0f;
        float _S2442 = _S2435 * _S2435;
        k_6 = _S2441 / _S2435;
        _S2438 = 0.0f;
        _S2439 = _S2442;
        _S2440 = _S2441;
    }
    else
    {
        float _S2443 = _S2434 * _S2434;
        k_6 = _S2436 / _S2434;
        _S2438 = _S2443;
        _S2439 = 0.0f;
        _S2440 = 0.0f;
    }
    float2  _S2444 = make_float2 (k_6);
    float2  _S2445 = _S2433 * make_float2 (k_6);
    float u_17 = _S2445.x;
    float v_17 = _S2445.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S2446 = k3_5 + r2_17 * k4_6;
    float _S2447 = k2_5 + r2_17 * _S2446;
    float _S2448 = k1_5 + r2_17 * _S2447;
    float radial_3 = 1.0f + r2_17 * _S2448;
    float2  _S2449 = make_float2 (radial_3);
    float _S2450 = _S2421 * u_17;
    float _S2451 = 2.0f * u_17;
    float _S2452 = r2_17 + _S2451 * u_17;
    float _S2453 = _S2425 * u_17;
    float _S2454 = 2.0f * v_17;
    float _S2455 = r2_17 + _S2454 * v_17;
    float2  _S2456 = _S2445 * make_float2 (radial_3) + make_float2 (_S2450 * v_17 + p2_6 * _S2452 + sx1_6 * r2_17, _S2453 * v_17 + p1_6 * _S2455 + sy1_6 * r2_17);
    float _S2457 = fx_18 * _S2456.x + cx_18;
    float _S2458 = fy_18 * _S2456.y + cy_18;
    float2  _S2459 = make_float2 (_S2457, _S2458);
    float2  _S2460 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S2461 = length_0(_S2460);
    float _S2462 = vert2_c_5.z;
    float _S2463 = s_primal_ctx_atan2_0(_S2461, _S2462);
    bool _S2464 = _S2463 < 0.00100000004749745f;
    float _S2465;
    float _S2466;
    float _S2467;
    if(_S2464)
    {
        float _S2468 = 1.0f - _S2463 * _S2463 / 3.0f;
        float _S2469 = _S2462 * _S2462;
        k_6 = _S2468 / _S2462;
        _S2465 = 0.0f;
        _S2466 = _S2469;
        _S2467 = _S2468;
    }
    else
    {
        float _S2470 = _S2461 * _S2461;
        k_6 = _S2463 / _S2461;
        _S2465 = _S2470;
        _S2466 = 0.0f;
        _S2467 = 0.0f;
    }
    float2  _S2471 = make_float2 (k_6);
    float2  _S2472 = _S2460 * make_float2 (k_6);
    float u_18 = _S2472.x;
    float v_18 = _S2472.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S2473 = k3_5 + r2_18 * k4_6;
    float _S2474 = k2_5 + r2_18 * _S2473;
    float _S2475 = k1_5 + r2_18 * _S2474;
    float radial_4 = 1.0f + r2_18 * _S2475;
    float2  _S2476 = make_float2 (radial_4);
    float _S2477 = _S2421 * u_18;
    float _S2478 = 2.0f * u_18;
    float _S2479 = r2_18 + _S2478 * u_18;
    float _S2480 = _S2425 * u_18;
    float _S2481 = 2.0f * v_18;
    float _S2482 = r2_18 + _S2481 * v_18;
    float2  _S2483 = _S2472 * make_float2 (radial_4) + make_float2 (_S2477 * v_18 + p2_6 * _S2479 + sx1_6 * r2_18, _S2480 * v_18 + p1_6 * _S2482 + sy1_6 * r2_18);
    float _S2484 = fx_18 * _S2483.x + cx_18;
    float _S2485 = fy_18 * _S2483.y + cy_18;
    float2  _S2486 = make_float2 (_S2484, _S2485);
    float2  e0_5 = _S2459 - _S2432;
    float2  e1_5 = _S2486 - _S2459;
    float2  e2_1 = _S2432 - _S2486;
    float _S2487 = e0_5.x;
    float _S2488 = e1_5.y;
    float _S2489 = e0_5.y;
    float _S2490 = e1_5.x;
    float _S2491 = _S2487 * _S2488 - _S2489 * _S2490;
    float _S2492 = 1.0f - hardness_5.y;
    float _S2493 = -1.0f / _S2492;
    float _S2494 = _S2492 * _S2492;
    float _S2495 = s_primal_ctx_max_0(_S2430, _S2457);
    float _S2496 = s_primal_ctx_min_0(_S2430, _S2457);
    float _S2497 = s_primal_ctx_max_0(_S2431, _S2458);
    float _S2498 = s_primal_ctx_min_0(_S2431, _S2458);
    Matrix<float, 3, 3>  _S2499 = transpose_0(R_18);
    float3  _S2500 = mean_18 - - s_primal_ctx_mul_1(_S2499, t_17);
    float _S2501 = _S2500.x;
    float _S2502 = _S2500.y;
    float _S2503 = _S2500.z;
    float _S2504 = _S2501 * _S2501 + _S2502 * _S2502 + _S2503 * _S2503;
    float _S2505 = s_primal_ctx_sqrt_0(_S2504);
    float x_59 = _S2501 / _S2505;
    float3  _S2506 = make_float3 (x_59);
    float _S2507 = _S2505 * _S2505;
    float y_37 = _S2502 / _S2505;
    float z_34 = _S2503 / _S2505;
    float3  _S2508 = make_float3 (z_34);
    float _S2509 = - y_37;
    float3  _S2510 = make_float3 (_S2509);
    float z2_34 = z_34 * z_34;
    float fTmp0B_14 = -1.09254848957061768f * z_34;
    float fC1_14 = x_59 * x_59 - y_37 * y_37;
    float _S2511 = 2.0f * x_59;
    float fS1_14 = _S2511 * y_37;
    float pSH6_4 = 0.94617468118667603f * z2_34 - 0.31539157032966614f;
    float3  _S2512 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_59;
    float3  _S2513 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_37;
    float3  _S2514 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2515 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2516 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_34 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_34;
    float _S2517 = 1.86588168144226074f * z2_34 - 1.11952900886535645f;
    float pSH12_4 = z_34 * _S2517;
    float3  _S2518 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_59;
    float3  _S2519 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_37;
    float3  _S2520 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2521 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2522 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_59 * fC1_14 - y_37 * fS1_14);
    float3  _S2523 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_59 * fS1_14 + y_37 * fC1_14);
    float3  _S2524 = make_float3 (pSH9_4);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2509) * (*sh_coeffs_14)[int(1)] + make_float3 (z_34) * (*sh_coeffs_14)[int(2)] - make_float3 (x_59) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]);
    float3  _S2525 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S2526 = make_float3 (0.0f);
    float3  _S2527 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S2528 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S2529 = make_float3 (_S2528);
    float3  _S2530 = (*ch_coeffs_5)[int(1)] * make_float3 (_S2528);
    float3  _S2531 = _S2527 + _S2530 + make_float3 (0.5f);
    float3  _S2532 = _S2527 - _S2530 + make_float3 (0.5f);
    float3  _S2533 = vert1_c_5 - vert0_c_5;
    float3  _S2534 = vert2_c_5 - vert0_c_5;
    float3  _S2535 = s_primal_ctx_cross_0(_S2533, _S2534);
    float3  _S2536 = normalize_0(_S2535);
    float3  _S2537 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S2536, mean_c_14)))))) * v_normal_1;
    float3  _S2538 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2539;
    (&_S2539)->primal_0 = _S2536;
    (&_S2539)->differential_0 = _S2538;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2540;
    (&_S2540)->primal_0 = mean_c_14;
    (&_S2540)->differential_0 = _S2538;
    s_bwd_prop_dot_0(&_S2539, &_S2540, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2541 = _S2540;
    float3  _S2542 = _S2537 + _S2539.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2543;
    (&_S2543)->primal_0 = _S2535;
    (&_S2543)->differential_0 = _S2538;
    s_bwd_normalize_impl_0(&_S2543, _S2542);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2544;
    (&_S2544)->primal_0 = _S2533;
    (&_S2544)->differential_0 = _S2538;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2545;
    (&_S2545)->primal_0 = _S2534;
    (&_S2545)->differential_0 = _S2538;
    s_bwd_prop_cross_0(&_S2544, &_S2545, _S2543.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2546 = _S2544;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2547 = _S2545;
    float3  _S2548 = - _S2545.differential_0;
    float3  _S2549 = - _S2544.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2550;
    (&_S2550)->primal_0 = _S2532;
    (&_S2550)->differential_0 = _S2538;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2551;
    (&_S2551)->primal_0 = _S2526;
    (&_S2551)->differential_0 = _S2538;
    s_bwd_prop_max_0(&_S2550, &_S2551, (*v_rgb_5)[int(2)]);
    float3  _S2552 = - _S2550.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2553;
    (&_S2553)->primal_0 = _S2531;
    (&_S2553)->differential_0 = _S2538;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2554;
    (&_S2554)->primal_0 = _S2526;
    (&_S2554)->differential_0 = _S2538;
    s_bwd_prop_max_0(&_S2553, &_S2554, (*v_rgb_5)[int(1)]);
    float3  _S2555 = _S2529 * (_S2552 + _S2553.differential_0);
    float3  _S2556 = _S2550.differential_0 + _S2553.differential_0;
    float3  _S2557 = make_float3 (0.5f) * - _S2556;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2558;
    (&_S2558)->primal_0 = _S2525;
    (&_S2558)->differential_0 = _S2538;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2559;
    (&_S2559)->primal_0 = _S2526;
    (&_S2559)->differential_0 = _S2538;
    s_bwd_prop_max_0(&_S2558, &_S2559, (*v_rgb_5)[int(0)]);
    float3  _S2560 = _S2557 + _S2558.differential_0;
    float3  _S2561 = _S2556 + _S2558.differential_0;
    float3  _S2562 = _S2523 * _S2561;
    float3  _S2563 = (*sh_coeffs_14)[int(15)] * _S2561;
    float3  _S2564 = _S2521 * _S2561;
    float3  _S2565 = (*sh_coeffs_14)[int(14)] * _S2561;
    float3  _S2566 = _S2519 * _S2561;
    float3  _S2567 = (*sh_coeffs_14)[int(13)] * _S2561;
    float3  _S2568 = _S2518 * _S2561;
    float3  _S2569 = (*sh_coeffs_14)[int(12)] * _S2561;
    float3  _S2570 = _S2520 * _S2561;
    float3  _S2571 = (*sh_coeffs_14)[int(11)] * _S2561;
    float3  _S2572 = _S2522 * _S2561;
    float3  _S2573 = (*sh_coeffs_14)[int(10)] * _S2561;
    float3  _S2574 = _S2524 * _S2561;
    float3  _S2575 = (*sh_coeffs_14)[int(9)] * _S2561;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2575.x + _S2575.y + _S2575.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2563.x + _S2563.y + _S2563.z);
    float _S2576 = _S2573.x + _S2573.y + _S2573.z;
    float _S2577 = _S2565.x + _S2565.y + _S2565.z;
    float _S2578 = _S2571.x + _S2571.y + _S2571.z;
    float _S2579 = _S2567.x + _S2567.y + _S2567.z;
    float _S2580 = _S2569.x + _S2569.y + _S2569.z;
    float _S2581 = - s_diff_fC2_T_4;
    float3  _S2582 = _S2515 * _S2561;
    float3  _S2583 = (*sh_coeffs_14)[int(8)] * _S2561;
    float3  _S2584 = _S2513 * _S2561;
    float3  _S2585 = (*sh_coeffs_14)[int(7)] * _S2561;
    float3  _S2586 = _S2512 * _S2561;
    float3  _S2587 = (*sh_coeffs_14)[int(6)] * _S2561;
    float3  _S2588 = _S2514 * _S2561;
    float3  _S2589 = (*sh_coeffs_14)[int(5)] * _S2561;
    float3  _S2590 = _S2516 * _S2561;
    float3  _S2591 = (*sh_coeffs_14)[int(4)] * _S2561;
    float _S2592 = _S2589.x + _S2589.y + _S2589.z;
    float _S2593 = _S2585.x + _S2585.y + _S2585.z;
    float _S2594 = fTmp1B_14 * _S2576 + x_59 * s_diff_fS2_T_4 + y_37 * _S2581 + 0.54627424478530884f * (_S2591.x + _S2591.y + _S2591.z);
    float _S2595 = fTmp1B_14 * _S2577 + y_37 * s_diff_fS2_T_4 + x_59 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2583.x + _S2583.y + _S2583.z);
    float _S2596 = y_37 * - _S2595;
    float _S2597 = x_59 * _S2595;
    float _S2598 = z_34 * (1.86588168144226074f * (z_34 * _S2580) + -2.28522896766662598f * (y_37 * _S2578 + x_59 * _S2579) + 0.94617468118667603f * (_S2587.x + _S2587.y + _S2587.z));
    float3  _S2599 = make_float3 (0.48860251903533936f) * _S2561;
    float3  _S2600 = - _S2599;
    float3  _S2601 = _S2506 * _S2600;
    float3  _S2602 = (*sh_coeffs_14)[int(3)] * _S2600;
    float3  _S2603 = _S2508 * _S2599;
    float3  _S2604 = (*sh_coeffs_14)[int(2)] * _S2599;
    float3  _S2605 = _S2510 * _S2599;
    float3  _S2606 = (*sh_coeffs_14)[int(1)] * _S2599;
    float _S2607 = (_S2517 * _S2580 + 1.44530570507049561f * (fS1_14 * _S2576 + fC1_14 * _S2577) + -1.09254848957061768f * (y_37 * _S2592 + x_59 * _S2593) + _S2598 + _S2598 + _S2604.x + _S2604.y + _S2604.z) / _S2507;
    float _S2608 = _S2505 * _S2607;
    float _S2609 = (fTmp0C_14 * _S2578 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2581 + fTmp0B_14 * _S2592 + _S2511 * _S2594 + _S2596 + _S2596 + - (_S2606.x + _S2606.y + _S2606.z)) / _S2507;
    float _S2610 = _S2505 * _S2609;
    float _S2611 = (fTmp0C_14 * _S2579 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2593 + 2.0f * (y_37 * _S2594) + _S2597 + _S2597 + _S2602.x + _S2602.y + _S2602.z) / _S2507;
    float _S2612 = _S2505 * _S2611;
    float _S2613 = _S2503 * - _S2607 + _S2502 * - _S2609 + _S2501 * - _S2611;
    DiffPair_float_0 _S2614;
    (&_S2614)->primal_0 = _S2504;
    (&_S2614)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2614, _S2613);
    float _S2615 = _S2503 * _S2614.differential_0;
    float _S2616 = _S2502 * _S2614.differential_0;
    float _S2617 = _S2501 * _S2614.differential_0;
    float3  _S2618 = make_float3 (0.282094806432724f) * _S2561;
    float3  _S2619 = make_float3 (_S2612 + _S2617 + _S2617, _S2610 + _S2616 + _S2616, _S2608 + _S2615 + _S2615);
    float3  _S2620 = - - _S2619;
    Matrix<float, 3, 3>  _S2621 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2622;
    (&_S2622)->primal_0 = _S2499;
    (&_S2622)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2623;
    (&_S2623)->primal_0 = t_17;
    (&_S2623)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2622, &_S2623, _S2620);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2624 = _S2623;
    Matrix<float, 3, 3>  _S2625 = transpose_0(_S2622.differential_0);
    DiffPair_float_0 _S2626;
    (&_S2626)->primal_0 = _S2498;
    (&_S2626)->differential_0 = 0.0f;
    DiffPair_float_0 _S2627;
    (&_S2627)->primal_0 = _S2485;
    (&_S2627)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2626, &_S2627, 0.0f);
    DiffPair_float_0 _S2628;
    (&_S2628)->primal_0 = _S2431;
    (&_S2628)->differential_0 = 0.0f;
    DiffPair_float_0 _S2629;
    (&_S2629)->primal_0 = _S2458;
    (&_S2629)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2628, &_S2629, _S2626.differential_0);
    DiffPair_float_0 _S2630;
    (&_S2630)->primal_0 = _S2497;
    (&_S2630)->differential_0 = 0.0f;
    DiffPair_float_0 _S2631;
    (&_S2631)->primal_0 = _S2485;
    (&_S2631)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2630, &_S2631, 0.0f);
    DiffPair_float_0 _S2632;
    (&_S2632)->primal_0 = _S2431;
    (&_S2632)->differential_0 = 0.0f;
    DiffPair_float_0 _S2633;
    (&_S2633)->primal_0 = _S2458;
    (&_S2633)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2632, &_S2633, _S2630.differential_0);
    DiffPair_float_0 _S2634;
    (&_S2634)->primal_0 = _S2496;
    (&_S2634)->differential_0 = 0.0f;
    DiffPair_float_0 _S2635;
    (&_S2635)->primal_0 = _S2484;
    (&_S2635)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2634, &_S2635, 0.0f);
    DiffPair_float_0 _S2636;
    (&_S2636)->primal_0 = _S2430;
    (&_S2636)->differential_0 = 0.0f;
    DiffPair_float_0 _S2637;
    (&_S2637)->primal_0 = _S2457;
    (&_S2637)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2636, &_S2637, _S2634.differential_0);
    DiffPair_float_0 _S2638;
    (&_S2638)->primal_0 = _S2495;
    (&_S2638)->differential_0 = 0.0f;
    DiffPair_float_0 _S2639;
    (&_S2639)->primal_0 = _S2484;
    (&_S2639)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2638, &_S2639, 0.0f);
    DiffPair_float_0 _S2640;
    (&_S2640)->primal_0 = _S2430;
    (&_S2640)->differential_0 = 0.0f;
    DiffPair_float_0 _S2641;
    (&_S2641)->primal_0 = _S2457;
    (&_S2641)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2640, &_S2641, _S2638.differential_0);
    DiffPair_float_0 _S2642;
    (&_S2642)->primal_0 = _S2493;
    (&_S2642)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2642, 0.0f);
    float _S2643 = - (-1.0f * - (_S2642.differential_0 / _S2494));
    float2  _S2644 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2645;
    (&_S2645)->primal_0 = e2_1;
    (&_S2645)->differential_0 = _S2644;
    s_bwd_length_impl_1(&_S2645, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2646;
    (&_S2646)->primal_0 = e1_5;
    (&_S2646)->differential_0 = _S2644;
    s_bwd_length_impl_1(&_S2646, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2647;
    (&_S2647)->primal_0 = e0_5;
    (&_S2647)->differential_0 = _S2644;
    s_bwd_length_impl_1(&_S2647, -0.0f);
    DiffPair_float_0 _S2648;
    (&_S2648)->primal_0 = _S2491;
    (&_S2648)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2648, 0.0f);
    float _S2649 = - _S2648.differential_0;
    float2  _S2650 = _S2646.differential_0 + make_float2 (_S2489 * _S2649, _S2487 * _S2648.differential_0);
    float2  _S2651 = _S2647.differential_0 + make_float2 (_S2488 * _S2648.differential_0, _S2490 * _S2649);
    float2  _S2652 = v_uv2_1 + - _S2645.differential_0 + _S2650;
    float _S2653 = fy_18 * (_S2627.differential_0 + _S2631.differential_0 + _S2652.y);
    float _S2654 = fx_18 * (_S2635.differential_0 + _S2639.differential_0 + _S2652.x);
    float2  _S2655 = make_float2 (_S2654, _S2653);
    float2  _S2656 = _S2472 * _S2655;
    float2  _S2657 = _S2476 * _S2655;
    float _S2658 = r2_18 * _S2653;
    float _S2659 = p1_6 * _S2653;
    float _S2660 = _S2482 * _S2653;
    float _S2661 = v_18 * _S2653;
    float _S2662 = u_18 * _S2661;
    float _S2663 = r2_18 * _S2654;
    float _S2664 = p2_6 * _S2654;
    float _S2665 = _S2479 * _S2654;
    float _S2666 = v_18 * _S2654;
    float _S2667 = u_18 * _S2666;
    float _S2668 = _S2656.x + _S2656.y;
    float _S2669 = r2_18 * _S2668;
    float _S2670 = r2_18 * _S2669;
    float _S2671 = r2_18 * _S2670;
    float _S2672 = r2_18 * _S2671;
    float _S2673 = sy1_6 * _S2653 + _S2659 + sx1_6 * _S2654 + _S2664 + _S2475 * _S2668 + _S2474 * _S2669 + _S2473 * _S2670 + k4_6 * _S2671;
    float _S2674 = v_18 * _S2673;
    float _S2675 = u_18 * _S2673;
    float _S2676 = _S2481 * _S2659 + 2.0f * (v_18 * _S2659) + _S2480 * _S2653 + _S2477 * _S2654 + _S2674 + _S2674;
    float _S2677 = _S2425 * _S2661 + _S2478 * _S2664 + 2.0f * (u_18 * _S2664) + _S2421 * _S2666 + _S2675 + _S2675;
    float2  _S2678 = v_uv0_1 + _S2645.differential_0 + - _S2651;
    float2  _S2679 = v_out_hardness_1 + make_float2 (0.0f, _S2643);
    float _S2680 = _S2637.differential_0 + _S2641.differential_0;
    float2  _S2681 = v_uv1_1 + - _S2650 + _S2651;
    float3  _S2682 = _S2548 + _S2549;
    FixedArray<float3 , 2>  _S2683;
    _S2683[int(0)] = _S2538;
    _S2683[int(1)] = _S2538;
    _S2683[int(1)] = _S2555;
    _S2683[int(0)] = _S2560;
    float3  _S2684 = _S2683[int(0)];
    float3  _S2685 = _S2683[int(1)];
    FixedArray<float3 , 16>  _S2686;
    _S2686[int(0)] = _S2538;
    _S2686[int(1)] = _S2538;
    _S2686[int(2)] = _S2538;
    _S2686[int(3)] = _S2538;
    _S2686[int(4)] = _S2538;
    _S2686[int(5)] = _S2538;
    _S2686[int(6)] = _S2538;
    _S2686[int(7)] = _S2538;
    _S2686[int(8)] = _S2538;
    _S2686[int(9)] = _S2538;
    _S2686[int(10)] = _S2538;
    _S2686[int(11)] = _S2538;
    _S2686[int(12)] = _S2538;
    _S2686[int(13)] = _S2538;
    _S2686[int(14)] = _S2538;
    _S2686[int(15)] = _S2538;
    _S2686[int(7)] = _S2584;
    _S2686[int(0)] = _S2618;
    _S2686[int(1)] = _S2605;
    _S2686[int(2)] = _S2603;
    _S2686[int(3)] = _S2601;
    _S2686[int(4)] = _S2590;
    _S2686[int(5)] = _S2588;
    _S2686[int(6)] = _S2586;
    _S2686[int(15)] = _S2562;
    _S2686[int(8)] = _S2582;
    _S2686[int(9)] = _S2574;
    _S2686[int(10)] = _S2572;
    _S2686[int(11)] = _S2570;
    _S2686[int(12)] = _S2568;
    _S2686[int(13)] = _S2566;
    _S2686[int(14)] = _S2564;
    float3  _S2687 = _S2686[int(0)];
    float3  _S2688 = _S2686[int(1)];
    float3  _S2689 = _S2686[int(2)];
    float3  _S2690 = _S2686[int(3)];
    float3  _S2691 = _S2686[int(4)];
    float3  _S2692 = _S2686[int(5)];
    float3  _S2693 = _S2686[int(6)];
    float3  _S2694 = _S2686[int(7)];
    float3  _S2695 = _S2686[int(8)];
    float3  _S2696 = _S2686[int(9)];
    float3  _S2697 = _S2686[int(10)];
    float3  _S2698 = _S2686[int(11)];
    float3  _S2699 = _S2686[int(12)];
    float3  _S2700 = _S2686[int(13)];
    float3  _S2701 = _S2686[int(14)];
    float3  _S2702 = _S2686[int(15)];
    float _S2703 = _S2636.differential_0 + _S2640.differential_0;
    float _S2704 = _S2628.differential_0 + _S2632.differential_0;
    float _S2705 = _S2629.differential_0 + _S2633.differential_0;
    float2  _S2706 = _S2657 + make_float2 (_S2677, _S2676);
    float2  _S2707 = _S2460 * _S2706;
    float2  _S2708 = _S2471 * _S2706;
    float _S2709 = _S2707.x + _S2707.y;
    if(_S2464)
    {
        float _S2710 = _S2709 / _S2466;
        float _S2711 = _S2467 * - _S2710;
        float _S2712 = _S2463 * (0.3333333432674408f * - (_S2462 * _S2710));
        k_6 = _S2712 + _S2712;
        _S2465 = _S2711;
        _S2466 = 0.0f;
    }
    else
    {
        float _S2713 = _S2709 / _S2465;
        float _S2714 = _S2463 * - _S2713;
        k_6 = _S2461 * _S2713;
        _S2465 = 0.0f;
        _S2466 = _S2714;
    }
    DiffPair_float_0 _S2715;
    (&_S2715)->primal_0 = _S2461;
    (&_S2715)->differential_0 = 0.0f;
    DiffPair_float_0 _S2716;
    (&_S2716)->primal_0 = _S2462;
    (&_S2716)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2715, &_S2716, k_6);
    float _S2717 = _S2716.differential_0 + _S2465;
    float _S2718 = _S2715.differential_0 + _S2466;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2719;
    (&_S2719)->primal_0 = _S2460;
    (&_S2719)->differential_0 = _S2644;
    s_bwd_length_impl_1(&_S2719, _S2718);
    float2  _S2720 = _S2719.differential_0 + _S2708;
    float _S2721 = fy_18 * (_S2681.y + _S2705);
    float _S2722 = fx_18 * (_S2681.x + _S2680);
    float2  _S2723 = make_float2 (_S2722, _S2721);
    float2  _S2724 = _S2445 * _S2723;
    float _S2725 = p1_6 * _S2721;
    float _S2726 = v_17 * _S2721;
    float _S2727 = p2_6 * _S2722;
    float _S2728 = v_17 * _S2722;
    float _S2729 = _S2724.x + _S2724.y;
    float _S2730 = r2_17 * _S2729;
    float _S2731 = r2_17 * _S2730;
    float _S2732 = r2_17 * _S2731;
    float _S2733 = sy1_6 * _S2721 + _S2725 + sx1_6 * _S2722 + _S2727 + _S2448 * _S2729 + _S2447 * _S2730 + _S2446 * _S2731 + k4_6 * _S2732;
    float _S2734 = v_17 * _S2733;
    float _S2735 = u_17 * _S2733;
    float3  _S2736 = _S2547.differential_0 + make_float3 (_S2720.x, _S2720.y, _S2717);
    float2  _S2737 = _S2449 * _S2723 + make_float2 (_S2425 * _S2726 + _S2451 * _S2727 + 2.0f * (u_17 * _S2727) + _S2421 * _S2728 + _S2735 + _S2735, _S2454 * _S2725 + 2.0f * (v_17 * _S2725) + _S2453 * _S2721 + _S2450 * _S2722 + _S2734 + _S2734);
    float _S2738 = u_17 * _S2726 + _S2662;
    float _S2739 = u_17 * _S2728 + _S2667;
    float _S2740 = r2_17 * _S2721 + _S2658;
    float _S2741 = r2_17 * _S2732 + _S2672;
    float _S2742 = _S2732 + _S2671;
    float _S2743 = _S2455 * _S2721 + _S2660;
    float _S2744 = _S2731 + _S2670;
    float _S2745 = r2_17 * _S2722 + _S2663;
    float _S2746 = _S2452 * _S2722 + _S2665;
    float _S2747 = _S2730 + _S2669;
    float2  _S2748 = _S2433 * _S2737;
    float2  _S2749 = _S2444 * _S2737;
    float _S2750 = _S2748.x + _S2748.y;
    if(_S2437)
    {
        float _S2751 = _S2750 / _S2439;
        float _S2752 = _S2440 * - _S2751;
        float _S2753 = _S2436 * (0.3333333432674408f * - (_S2435 * _S2751));
        k_6 = _S2753 + _S2753;
        _S2438 = _S2752;
        _S2439 = 0.0f;
    }
    else
    {
        float _S2754 = _S2750 / _S2438;
        float _S2755 = _S2436 * - _S2754;
        k_6 = _S2434 * _S2754;
        _S2438 = 0.0f;
        _S2439 = _S2755;
    }
    DiffPair_float_0 _S2756;
    (&_S2756)->primal_0 = _S2434;
    (&_S2756)->differential_0 = 0.0f;
    DiffPair_float_0 _S2757;
    (&_S2757)->primal_0 = _S2435;
    (&_S2757)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2756, &_S2757, k_6);
    float _S2758 = _S2757.differential_0 + _S2438;
    float _S2759 = _S2756.differential_0 + _S2439;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2760;
    (&_S2760)->primal_0 = _S2433;
    (&_S2760)->differential_0 = _S2644;
    s_bwd_length_impl_1(&_S2760, _S2759);
    float2  _S2761 = _S2760.differential_0 + _S2749;
    float _S2762 = fy_18 * (_S2678.y + _S2704);
    float _S2763 = fx_18 * (_S2678.x + _S2703);
    float2  _S2764 = make_float2 (_S2763, _S2762);
    float2  _S2765 = _S2416 * _S2764;
    float _S2766 = p1_6 * _S2762;
    float _S2767 = v_16 * _S2762;
    float _S2768 = p2_6 * _S2763;
    float _S2769 = v_16 * _S2763;
    float _S2770 = _S2765.x + _S2765.y;
    float _S2771 = r2_16 * _S2770;
    float _S2772 = r2_16 * _S2771;
    float _S2773 = r2_16 * _S2772;
    float _S2774 = sy1_6 * _S2762 + _S2766 + sx1_6 * _S2763 + _S2768 + _S2419 * _S2770 + _S2418 * _S2771 + _S2417 * _S2772 + k4_6 * _S2773;
    float _S2775 = v_16 * _S2774;
    float _S2776 = u_16 * _S2774;
    float2  _S2777 = make_float2 (r2_16 * _S2763 + _S2745, r2_16 * _S2762 + _S2740);
    float2  _S2778 = make_float2 (_S2428 * _S2762 + 2.0f * (u_16 * _S2769 + _S2739) + _S2743, 2.0f * (u_16 * _S2767 + _S2738) + _S2424 * _S2763 + _S2746);
    float4  _S2779 = make_float4 (_S2771 + _S2747, _S2772 + _S2744, _S2773 + _S2742, r2_16 * _S2773 + _S2741);
    float3  _S2780 = _S2546.differential_0 + make_float3 (_S2761.x, _S2761.y, _S2758);
    float2  _S2781 = _S2420 * _S2764 + make_float2 (_S2425 * _S2767 + _S2423 * _S2768 + 2.0f * (u_16 * _S2768) + _S2421 * _S2769 + _S2776 + _S2776, _S2427 * _S2766 + 2.0f * (v_16 * _S2766) + _S2426 * _S2762 + _S2422 * _S2763 + _S2775 + _S2775);
    CameraDistortion_0 _S2782 = CameraDistortion_x24_syn_dzero_0();
    (&_S2782)->thin_prism_coeffs_0 = _S2777;
    (&_S2782)->tangential_coeffs_0 = _S2778;
    (&_S2782)->radial_coeffs_0 = _S2779;
    float2  _S2783 = _S2404 * _S2781;
    float2  _S2784 = _S2415 * _S2781;
    float _S2785 = _S2783.x + _S2783.y;
    if(_S2408)
    {
        float _S2786 = _S2785 / _S2410;
        float _S2787 = _S2411 * - _S2786;
        float _S2788 = _S2407 * (0.3333333432674408f * - (_S2406 * _S2786));
        k_6 = _S2788 + _S2788;
        _S2409 = _S2787;
        _S2410 = 0.0f;
    }
    else
    {
        float _S2789 = _S2785 / _S2409;
        float _S2790 = _S2407 * - _S2789;
        k_6 = _S2405 * _S2789;
        _S2409 = 0.0f;
        _S2410 = _S2790;
    }
    DiffPair_float_0 _S2791;
    (&_S2791)->primal_0 = _S2405;
    (&_S2791)->differential_0 = 0.0f;
    DiffPair_float_0 _S2792;
    (&_S2792)->primal_0 = _S2406;
    (&_S2792)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2791, &_S2792, k_6);
    float _S2793 = _S2792.differential_0 + _S2409;
    float _S2794 = _S2791.differential_0 + _S2410;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2795;
    (&_S2795)->primal_0 = _S2404;
    (&_S2795)->differential_0 = _S2644;
    s_bwd_length_impl_1(&_S2795, _S2794);
    float2  _S2796 = _S2795.differential_0 + _S2784;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2797;
    (&_S2797)->primal_0 = vert2_c_5;
    (&_S2797)->differential_0 = _S2538;
    s_bwd_length_impl_0(&_S2797, v_depth_5.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2798;
    (&_S2798)->primal_0 = vert1_c_5;
    (&_S2798)->differential_0 = _S2538;
    s_bwd_length_impl_0(&_S2798, v_depth_5.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2799;
    (&_S2799)->primal_0 = vert0_c_5;
    (&_S2799)->differential_0 = _S2538;
    s_bwd_length_impl_0(&_S2799, v_depth_5.x);
    float3  _S2800 = _S2797.differential_0 + _S2736;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2801;
    (&_S2801)->primal_0 = R_18;
    (&_S2801)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2802;
    (&_S2802)->primal_0 = vert2_1;
    (&_S2802)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2801, &_S2802, _S2800);
    float3  _S2803 = _S2798.differential_0 + _S2780;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2804;
    (&_S2804)->primal_0 = R_18;
    (&_S2804)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2805;
    (&_S2805)->primal_0 = vert1_1;
    (&_S2805)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2804, &_S2805, _S2803);
    float3  _S2806 = _S2799.differential_0 + _S2682 + make_float3 (_S2796.x, _S2796.y, _S2793);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2807;
    (&_S2807)->primal_0 = R_18;
    (&_S2807)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2808;
    (&_S2808)->primal_0 = vert0_1;
    (&_S2808)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2807, &_S2808, _S2806);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2809;
    (&_S2809)->primal_0 = _S2397;
    (&_S2809)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2810;
    (&_S2810)->primal_0 = _S2402;
    (&_S2810)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2809, &_S2810, _S2802.differential_0);
    float _S2811 = - _S2810.differential_0.y;
    float _S2812 = _S2401 * _S2810.differential_0.x;
    float _S2813 = - (_S2391 * _S2810.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2814;
    (&_S2814)->primal_0 = _S2397;
    (&_S2814)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2815;
    (&_S2815)->primal_0 = _S2400;
    (&_S2815)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2814, &_S2815, _S2805.differential_0);
    float _S2816 = _S2391 * _S2815.differential_0.x;
    float _S2817 = _S2399 * _S2815.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2818;
    (&_S2818)->primal_0 = _S2397;
    (&_S2818)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2819;
    (&_S2819)->primal_0 = _S2398;
    (&_S2819)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2818, &_S2819, _S2808.differential_0);
    Matrix<float, 3, 3>  _S2820 = transpose_0(_S2809.differential_0 + _S2814.differential_0 + _S2818.differential_0);
    float _S2821 = 2.0f * - _S2820.rows[int(2)].z;
    float _S2822 = 2.0f * _S2820.rows[int(2)].y;
    float _S2823 = 2.0f * _S2820.rows[int(2)].x;
    float _S2824 = 2.0f * _S2820.rows[int(1)].z;
    float _S2825 = 2.0f * - _S2820.rows[int(1)].y;
    float _S2826 = 2.0f * _S2820.rows[int(1)].x;
    float _S2827 = 2.0f * _S2820.rows[int(0)].z;
    float _S2828 = 2.0f * _S2820.rows[int(0)].y;
    float _S2829 = 2.0f * - _S2820.rows[int(0)].x;
    float _S2830 = - _S2826 + _S2828;
    float _S2831 = _S2823 + - _S2827;
    float _S2832 = - _S2822 + _S2824;
    float _S2833 = _S2822 + _S2824;
    float _S2834 = _S2823 + _S2827;
    float _S2835 = _S2826 + _S2828;
    float _S2836 = z_33 * (_S2825 + _S2829);
    float _S2837 = y_36 * (_S2821 + _S2829);
    float _S2838 = x_58 * (_S2821 + _S2825);
    float _S2839 = z_33 * _S2830 + y_36 * _S2831 + x_58 * _S2832;
    float _S2840 = _S2396 * _S2839;
    float _S2841 = w_19 * _S2830 + y_36 * _S2833 + x_58 * _S2834 + _S2836 + _S2836;
    float _S2842 = _S2396 * _S2841;
    float _S2843 = w_19 * _S2831 + z_33 * _S2833 + x_58 * _S2835 + _S2837 + _S2837;
    float _S2844 = _S2396 * _S2843;
    float _S2845 = w_19 * _S2832 + z_33 * _S2834 + y_36 * _S2835 + _S2838 + _S2838;
    float _S2846 = _S2396 * _S2845;
    float _S2847 = quat_21.x * _S2839 + quat_21.w * _S2841 + quat_21.z * _S2843 + quat_21.y * _S2845;
    DiffPair_float_0 _S2848;
    (&_S2848)->primal_0 = _S2395;
    (&_S2848)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S2848, _S2847);
    float _S2849 = quat_21.x * _S2848.differential_0;
    float _S2850 = quat_21.w * _S2848.differential_0;
    float _S2851 = quat_21.z * _S2848.differential_0;
    float _S2852 = quat_21.y * _S2848.differential_0;
    float _S2853 = _S2842 + _S2850 + _S2850;
    float _S2854 = _S2844 + _S2851 + _S2851;
    float _S2855 = _S2846 + _S2852 + _S2852;
    float _S2856 = _S2840 + _S2849 + _S2849;
    float _S2857 = _S2813 + _S2816;
    float _S2858 = 0.5f * - _S2857;
    float _S2859 = _S2811 + _S2815.differential_0.y;
    DiffPair_float_0 _S2860;
    (&_S2860)->primal_0 = _S2392;
    (&_S2860)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2860, _S2859);
    float _S2861 = _S2858 + _S2860.differential_0;
    float _S2862 = _S2812 + _S2817 + _S2819.differential_0.x;
    DiffPair_float_0 _S2863;
    (&_S2863)->primal_0 = _S2390;
    (&_S2863)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2863, _S2862);
    float _S2864 = _S2858 + _S2863.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2865;
    (&_S2865)->primal_0 = mean_c_14;
    (&_S2865)->differential_0 = _S2538;
    s_bwd_length_impl_0(&_S2865, 0.0f);
    float3  _S2866 = _S2865.differential_0 + _S2541.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2867;
    (&_S2867)->primal_0 = R_18;
    (&_S2867)->differential_0 = _S2621;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2868;
    (&_S2868)->primal_0 = mean_18;
    (&_S2868)->differential_0 = _S2538;
    s_bwd_prop_mul_1(&_S2867, &_S2868, _S2866);
    float3  _S2869 = _S2800 + _S2803 + _S2806 + _S2866 + _S2624.differential_0;
    Matrix<float, 3, 3>  _S2870 = _S2801.differential_0 + _S2804.differential_0 + _S2807.differential_0 + _S2867.differential_0 + _S2625;
    float3  _S2871 = make_float3 (_S2864, _S2861, _S2857);
    float4  _S2872 = make_float4 (0.0f);
    *&((&_S2872)->w) = _S2853;
    *&((&_S2872)->z) = _S2854;
    *&((&_S2872)->y) = _S2855;
    *&((&_S2872)->x) = _S2856;
    float4  _S2873 = _S2872;
    float3  _S2874 = _S2802.differential_0 + _S2805.differential_0 + _S2808.differential_0 + _S2868.differential_0 + _S2619;
    *v_mean_6 = _S2874;
    *v_quat_6 = _S2873;
    *v_scale_6 = _S2871;
    *v_hardness_1 = _S2679;
    (*v_sh_coeffs_4)[int(0)] = _S2687;
    (*v_sh_coeffs_4)[int(1)] = _S2688;
    (*v_sh_coeffs_4)[int(2)] = _S2689;
    (*v_sh_coeffs_4)[int(3)] = _S2690;
    (*v_sh_coeffs_4)[int(4)] = _S2691;
    (*v_sh_coeffs_4)[int(5)] = _S2692;
    (*v_sh_coeffs_4)[int(6)] = _S2693;
    (*v_sh_coeffs_4)[int(7)] = _S2694;
    (*v_sh_coeffs_4)[int(8)] = _S2695;
    (*v_sh_coeffs_4)[int(9)] = _S2696;
    (*v_sh_coeffs_4)[int(10)] = _S2697;
    (*v_sh_coeffs_4)[int(11)] = _S2698;
    (*v_sh_coeffs_4)[int(12)] = _S2699;
    (*v_sh_coeffs_4)[int(13)] = _S2700;
    (*v_sh_coeffs_4)[int(14)] = _S2701;
    (*v_sh_coeffs_4)[int(15)] = _S2702;
    (*v_ch_coeffs_1)[int(0)] = _S2684;
    (*v_ch_coeffs_1)[int(1)] = _S2685;
    *v_R_5 = _S2870;
    *v_t_5 = _S2869;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_19)
{
    DiffPair_float_0 _S2875 = *dpx_16;
    bool _S2876;
    if(((*dpx_16).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S2876 = ((*dpx_16).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S2876 = false;
    }
    float _S2877;
    if(_S2876)
    {
        _S2877 = dOut_19;
    }
    else
    {
        _S2877 = 0.0f;
    }
    dpx_16->primal_0 = _S2875.primal_0;
    dpx_16->differential_0 = _S2877;
    DiffPair_float_0 _S2878 = *dpMin_0;
    if((_S2875.primal_0) < ((*dpMin_0).primal_0))
    {
        _S2877 = dOut_19;
    }
    else
    {
        _S2877 = 0.0f;
    }
    dpMin_0->primal_0 = _S2878.primal_0;
    dpMin_0->differential_0 = _S2877;
    DiffPair_float_0 _S2879 = *dpMax_0;
    if(((*dpx_16).primal_0) > ((*dpMax_0).primal_0))
    {
        _S2877 = dOut_19;
    }
    else
    {
        _S2877 = 0.0f;
    }
    dpMax_0->primal_0 = _S2879.primal_0;
    dpMax_0->differential_0 = _S2877;
    return;
}

inline __device__ float clamp_0(float x_60, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_60), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_17, DiffPair_float_0 * dpy_6, float dOut_20)
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
        DiffPair_float_0 _S2880 = *dpx_17;
        float _S2881 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S2881;
        float _S2882 = val_0 * (F32_log((_S2880.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S2882;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S2883 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S2883))));
    float2  _S2884 = p_0 - v0_0;
    float2  _S2885 = normalize_1(e0_6);
    float2  _S2886 = p_0 - v1_0;
    float2  _S2887 = normalize_1(e1_6);
    float2  _S2888 = p_0 - v2_0;
    float2  _S2889 = normalize_1(e2_2);
    float _S2890 = hardness_6.x;
    float _S2891 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S2884.x * _S2885.y - _S2884.y * _S2885.x)), (se_0 * (_S2886.x * _S2887.y - _S2886.y * _S2887.x))))), (se_0 * (_S2888.x * _S2889.y - _S2888.y * _S2889.x)))) / ((F32_abs((_S2883))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S2891))));
    float _S2892;
    if(a_1 <= 0.0f)
    {
        _S2892 = 0.0f;
    }
    else
    {
        _S2892 = (F32_min(((F32_pow((a_1), (_S2891)))), (0.99900001287460327f)));
    }
    return _S2890 * _S2892;
}

inline __device__ float s_primal_ctx_abs_0(float _S2893)
{
    return (F32_abs((_S2893)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S2894, float _S2895, float _S2896)
{
    return clamp_0(_S2894, _S2895, _S2896);
}

inline __device__ float s_primal_ctx_exp2_0(float _S2897)
{
    return (F32_exp2((_S2897)));
}

inline __device__ float s_primal_ctx_pow_0(float _S2898, float _S2899)
{
    return (F32_pow((_S2898), (_S2899)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S2900, DiffPair_float_0 * _S2901, float _S2902)
{
    _d_pow_0(_S2900, _S2901, _S2902);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S2903, DiffPair_float_0 * _S2904, DiffPair_float_0 * _S2905, float _S2906)
{
    _d_clamp_0(_S2903, _S2904, _S2905, _S2906);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_8)
{
    float _S2907 = length_0((*dpx_18).primal_0);
    float2  _S2908 = (*dpx_18).primal_0 * _s_dOut_8;
    float2  _S2909 = make_float2 (1.0f / _S2907) * _s_dOut_8;
    float _S2910 = - ((_S2908.x + _S2908.y) / (_S2907 * _S2907));
    float2  _S2911 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2912;
    (&_S2912)->primal_0 = (*dpx_18).primal_0;
    (&_S2912)->differential_0 = _S2911;
    s_bwd_length_impl_1(&_S2912, _S2910);
    float2  _S2913 = _S2909 + _S2912.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S2913;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2914, float2  _S2915)
{
    s_bwd_prop_normalize_impl_1(_S2914, _S2915);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_9)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S2916 = e0_7.x;
    float _S2917 = e1_7.y;
    float _S2918 = e0_7.y;
    float _S2919 = e1_7.x;
    float _S2920 = _S2916 * _S2917 - _S2918 * _S2919;
    float se_1 = float((F32_sign((_S2920))));
    float2  _S2921 = p_1 - (*dpv0_0).primal_0;
    float2  _S2922 = normalize_1(e0_7);
    float _S2923 = _S2921.x;
    float _S2924 = _S2922.y;
    float _S2925 = _S2921.y;
    float _S2926 = _S2922.x;
    float de0_0 = se_1 * (_S2923 * _S2924 - _S2925 * _S2926);
    float2  _S2927 = p_1 - (*dpv1_0).primal_0;
    float2  _S2928 = normalize_1(e1_7);
    float _S2929 = _S2927.x;
    float _S2930 = _S2928.y;
    float _S2931 = _S2927.y;
    float _S2932 = _S2928.x;
    float de1_0 = se_1 * (_S2929 * _S2930 - _S2931 * _S2932);
    float2  _S2933 = p_1 - (*dpv2_0).primal_0;
    float2  _S2934 = normalize_1(e2_3);
    float _S2935 = _S2933.x;
    float _S2936 = _S2934.y;
    float _S2937 = _S2933.y;
    float _S2938 = _S2934.x;
    float de2_0 = se_1 * (_S2935 * _S2936 - _S2937 * _S2938);
    float _S2939 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S2940 = s_primal_ctx_max_0(_S2939, de2_0);
    float _S2941 = s_primal_ctx_abs_0(_S2920);
    float _S2942 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S2941 / _S2942;
    float _S2943 = _S2942 * _S2942;
    float _S2944 = (*dphardness_0).primal_0.x;
    float _S2945 = (*dphardness_0).primal_0.y;
    float _S2946 = dmax_0 * dmax_0;
    float _S2947 = 1.0f + _S2940 / dmax_0;
    float _S2948 = 1.0f - s_primal_ctx_clamp_0(_S2945, 0.00499999988824129f, 0.98000001907348633f);
    float _S2949 = -1.0f / _S2948;
    float _S2950 = _S2948 * _S2948;
    float _S2951 = 1.0f - s_primal_ctx_exp2_0(_S2949);
    float a_2 = 1.0f - _S2947 * _S2951;
    bool _S2952 = a_2 <= 0.0f;
    float _S2953;
    float _S2954;
    if(_S2952)
    {
        _S2953 = 0.0f;
        _S2954 = 0.0f;
    }
    else
    {
        float _S2955 = s_primal_ctx_pow_0(a_2, _S2948);
        _S2953 = s_primal_ctx_min_0(_S2955, 0.99900001287460327f);
        _S2954 = _S2955;
    }
    float _S2956 = _S2944 * _s_dOut_9;
    float _S2957 = _S2953 * _s_dOut_9;
    if(_S2952)
    {
        _S2953 = 0.0f;
        _S2954 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S2958;
        (&_S2958)->primal_0 = _S2954;
        (&_S2958)->differential_0 = 0.0f;
        DiffPair_float_0 _S2959;
        (&_S2959)->primal_0 = 0.99900001287460327f;
        (&_S2959)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S2958, &_S2959, _S2956);
        DiffPair_float_0 _S2960;
        (&_S2960)->primal_0 = a_2;
        (&_S2960)->differential_0 = 0.0f;
        DiffPair_float_0 _S2961;
        (&_S2961)->primal_0 = _S2948;
        (&_S2961)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S2960, &_S2961, _S2958.differential_0);
        _S2953 = _S2960.differential_0;
        _S2954 = _S2961.differential_0;
    }
    float _S2962 = - _S2953;
    float _S2963 = _S2951 * _S2962;
    float _S2964 = - (_S2947 * _S2962);
    DiffPair_float_0 _S2965;
    (&_S2965)->primal_0 = _S2949;
    (&_S2965)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2965, _S2964);
    float _S2966 = - (-1.0f * - (_S2965.differential_0 / _S2950) + _S2954);
    float _S2967 = _S2963 / _S2946;
    float s_diff_dmax_T_0 = _S2940 * - _S2967;
    float _S2968 = dmax_0 * _S2967;
    DiffPair_float_0 _S2969;
    (&_S2969)->primal_0 = _S2945;
    (&_S2969)->differential_0 = 0.0f;
    DiffPair_float_0 _S2970;
    (&_S2970)->primal_0 = 0.00499999988824129f;
    (&_S2970)->differential_0 = 0.0f;
    DiffPair_float_0 _S2971;
    (&_S2971)->primal_0 = 0.98000001907348633f;
    (&_S2971)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2969, &_S2970, &_S2971, _S2966);
    float _S2972 = s_diff_dmax_T_0 / _S2943;
    float _S2973 = _S2941 * - _S2972;
    float _S2974 = _S2942 * _S2972;
    float2  _S2975 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2976;
    (&_S2976)->primal_0 = e2_3;
    (&_S2976)->differential_0 = _S2975;
    s_bwd_length_impl_1(&_S2976, _S2973);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2977;
    (&_S2977)->primal_0 = e1_7;
    (&_S2977)->differential_0 = _S2975;
    s_bwd_length_impl_1(&_S2977, _S2973);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2978;
    (&_S2978)->primal_0 = e0_7;
    (&_S2978)->differential_0 = _S2975;
    s_bwd_length_impl_1(&_S2978, _S2973);
    DiffPair_float_0 _S2979;
    (&_S2979)->primal_0 = _S2920;
    (&_S2979)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2979, _S2974);
    DiffPair_float_0 _S2980;
    (&_S2980)->primal_0 = _S2939;
    (&_S2980)->differential_0 = 0.0f;
    DiffPair_float_0 _S2981;
    (&_S2981)->primal_0 = de2_0;
    (&_S2981)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2980, &_S2981, _S2968);
    DiffPair_float_0 _S2982;
    (&_S2982)->primal_0 = de0_0;
    (&_S2982)->differential_0 = 0.0f;
    DiffPair_float_0 _S2983;
    (&_S2983)->primal_0 = de1_0;
    (&_S2983)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2982, &_S2983, _S2980.differential_0);
    float _S2984 = se_1 * _S2981.differential_0;
    float _S2985 = - _S2984;
    float _S2986 = _S2938 * _S2985;
    float _S2987 = _S2936 * _S2984;
    float2  _S2988 = make_float2 (_S2937 * _S2985, _S2935 * _S2984);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2989;
    (&_S2989)->primal_0 = e2_3;
    (&_S2989)->differential_0 = _S2975;
    s_bwd_normalize_impl_1(&_S2989, _S2988);
    float2  _S2990 = - make_float2 (_S2987, _S2986);
    float _S2991 = se_1 * _S2983.differential_0;
    float _S2992 = - _S2991;
    float _S2993 = _S2932 * _S2992;
    float _S2994 = _S2930 * _S2991;
    float2  _S2995 = make_float2 (_S2931 * _S2992, _S2929 * _S2991);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2996;
    (&_S2996)->primal_0 = e1_7;
    (&_S2996)->differential_0 = _S2975;
    s_bwd_normalize_impl_1(&_S2996, _S2995);
    float2  _S2997 = - make_float2 (_S2994, _S2993);
    float _S2998 = se_1 * _S2982.differential_0;
    float _S2999 = - _S2998;
    float _S3000 = _S2926 * _S2999;
    float _S3001 = _S2924 * _S2998;
    float2  _S3002 = make_float2 (_S2925 * _S2999, _S2923 * _S2998);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3003;
    (&_S3003)->primal_0 = e0_7;
    (&_S3003)->differential_0 = _S2975;
    s_bwd_normalize_impl_1(&_S3003, _S3002);
    float2  _S3004 = - make_float2 (_S3001, _S3000);
    float _S3005 = - _S2979.differential_0;
    float2  _S3006 = _S2976.differential_0 + _S2989.differential_0;
    float2  _S3007 = - _S3006;
    float2  _S3008 = _S2977.differential_0 + _S2996.differential_0 + make_float2 (_S2918 * _S3005, _S2916 * _S2979.differential_0);
    float2  _S3009 = - _S3008;
    float2  _S3010 = _S2978.differential_0 + _S3003.differential_0 + make_float2 (_S2917 * _S2979.differential_0, _S2919 * _S3005);
    float2  _S3011 = - _S3010;
    float2  _S3012 = make_float2 (_S2957, _S2969.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S3012;
    float2  _S3013 = _S2990 + _S3007 + _S3008;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S3013;
    float2  _S3014 = _S2997 + _S3009 + _S3010;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S3014;
    float2  _S3015 = _S3004 + _S3006 + _S3011;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S3015;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3016, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3017, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3018, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3019, float2  _S3020, float _S3021)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S3016, _S3017, _S3018, _S3019, _S3020, _S3021);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S3022 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S3022;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S3022;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S3022;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S3022;
    s_bwd_evaluate_alpha_opaque_triangle_fast_0(&dp_v0_0, &dp_v1_0, &dp_v2_0, &dp_hardness_0, p_2, v_alpha_1);
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
    float2  _S3023 = p_3 - v0_2;
    float2  _S3024 = p_3 - v1_2;
    float2  _S3025 = p_3 - v2_2;
    float _S3026 = e0_8.x;
    float _S3027 = e1_8.y;
    float _S3028 = e0_8.y;
    float _S3029 = e1_8.x;
    float _S3030 = _S3026 * _S3027 - _S3028 * _S3029;
    float se_2 = float((F32_sign((_S3030))));
    float _S3031 = hardness_8.x;
    float _S3032 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S3023.x * _S3028 - _S3023.y * _S3026)), (se_2 * (_S3024.x * _S3027 - _S3024.y * _S3029))))), (se_2 * (_S3025.x * e2_4.y - _S3025.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S3023 - e0_8 * make_float2 (clamp_0(dot_1(_S3023, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S3024 - e1_8 * make_float2 (clamp_0(dot_1(_S3024, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S3025 - e2_4 * make_float2 (clamp_0(dot_1(_S3025, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S3030))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S3032))));
    float _S3033;
    if(a_3 <= 0.0f)
    {
        _S3033 = 0.0f;
    }
    else
    {
        _S3033 = (F32_min(((F32_pow((a_3), (_S3032)))), (0.99900001287460327f)));
    }
    return _S3031 * _S3033;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S3034, float2  _S3035)
{
    return dot_1(_S3034, _S3035);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3036, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3037, float _S3038)
{
    _d_dot_1(_S3036, _S3037, _S3038);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_10)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S3039 = p_4 - (*dpv0_1).primal_0;
    float _S3040 = s_primal_ctx_dot_1(_S3039, e0_9);
    float _S3041 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S3042 = _S3040 / _S3041;
    float _S3043 = _S3041 * _S3041;
    float _S3044 = s_primal_ctx_clamp_0(_S3042, 0.0f, 1.0f);
    float2  _S3045 = make_float2 (_S3044);
    float2  _S3046 = _S3039 - e0_9 * make_float2 (_S3044);
    float _S3047 = length_0(_S3046);
    float2  _S3048 = p_4 - (*dpv1_1).primal_0;
    float _S3049 = s_primal_ctx_dot_1(_S3048, e1_9);
    float _S3050 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S3051 = _S3049 / _S3050;
    float _S3052 = _S3050 * _S3050;
    float _S3053 = s_primal_ctx_clamp_0(_S3051, 0.0f, 1.0f);
    float2  _S3054 = make_float2 (_S3053);
    float2  _S3055 = _S3048 - e1_9 * make_float2 (_S3053);
    float _S3056 = length_0(_S3055);
    float2  _S3057 = p_4 - (*dpv2_1).primal_0;
    float _S3058 = s_primal_ctx_dot_1(_S3057, e2_5);
    float _S3059 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S3060 = _S3058 / _S3059;
    float _S3061 = _S3059 * _S3059;
    float _S3062 = s_primal_ctx_clamp_0(_S3060, 0.0f, 1.0f);
    float2  _S3063 = make_float2 (_S3062);
    float2  _S3064 = _S3057 - e2_5 * make_float2 (_S3062);
    float _S3065 = length_0(_S3064);
    float _S3066 = e0_9.x;
    float _S3067 = e1_9.y;
    float _S3068 = e0_9.y;
    float _S3069 = e1_9.x;
    float _S3070 = _S3066 * _S3067 - _S3068 * _S3069;
    float se_3 = float((F32_sign((_S3070))));
    float _S3071 = _S3039.x;
    float _S3072 = _S3039.y;
    float s0_0 = se_3 * (_S3071 * _S3068 - _S3072 * _S3066);
    float _S3073 = _S3048.x;
    float _S3074 = _S3048.y;
    float s1_0 = se_3 * (_S3073 * _S3067 - _S3074 * _S3069);
    float _S3075 = _S3057.x;
    float _S3076 = e2_5.y;
    float _S3077 = _S3057.y;
    float _S3078 = e2_5.x;
    float s2_0 = se_3 * (_S3075 * _S3076 - _S3077 * _S3078);
    float _S3079 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S3079, s2_0)))));
    float _S3080 = s_primal_ctx_min_0(_S3047, _S3056);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S3080, _S3065);
    float _S3081 = s_primal_ctx_abs_0(_S3070);
    float _S3082 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S3081 / _S3082;
    float _S3083 = _S3082 * _S3082;
    float _S3084 = (*dphardness_1).primal_0.x;
    float _S3085 = (*dphardness_1).primal_0.y;
    float _S3086 = dmax_1 * dmax_1;
    float _S3087 = 1.0f + dv_0 / dmax_1;
    float _S3088 = 1.0f - s_primal_ctx_clamp_0(_S3085, 0.00499999988824129f, 0.98000001907348633f);
    float _S3089 = -1.0f / _S3088;
    float _S3090 = _S3088 * _S3088;
    float _S3091 = 1.0f - s_primal_ctx_exp2_0(_S3089);
    float a_4 = 1.0f - _S3087 * _S3091;
    bool _S3092 = a_4 <= 0.0f;
    float _S3093;
    float _S3094;
    if(_S3092)
    {
        _S3093 = 0.0f;
        _S3094 = 0.0f;
    }
    else
    {
        float _S3095 = s_primal_ctx_pow_0(a_4, _S3088);
        _S3093 = s_primal_ctx_min_0(_S3095, 0.99900001287460327f);
        _S3094 = _S3095;
    }
    float _S3096 = _S3084 * _s_dOut_10;
    float _S3097 = _S3093 * _s_dOut_10;
    if(_S3092)
    {
        _S3093 = 0.0f;
        _S3094 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S3098;
        (&_S3098)->primal_0 = _S3094;
        (&_S3098)->differential_0 = 0.0f;
        DiffPair_float_0 _S3099;
        (&_S3099)->primal_0 = 0.99900001287460327f;
        (&_S3099)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S3098, &_S3099, _S3096);
        DiffPair_float_0 _S3100;
        (&_S3100)->primal_0 = a_4;
        (&_S3100)->differential_0 = 0.0f;
        DiffPair_float_0 _S3101;
        (&_S3101)->primal_0 = _S3088;
        (&_S3101)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S3100, &_S3101, _S3098.differential_0);
        _S3093 = _S3100.differential_0;
        _S3094 = _S3101.differential_0;
    }
    float _S3102 = - _S3093;
    float _S3103 = _S3091 * _S3102;
    float _S3104 = - (_S3087 * _S3102);
    DiffPair_float_0 _S3105;
    (&_S3105)->primal_0 = _S3089;
    (&_S3105)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3105, _S3104);
    float _S3106 = - (-1.0f * - (_S3105.differential_0 / _S3090) + _S3094);
    float _S3107 = _S3103 / _S3086;
    float s_diff_dmax_T_1 = dv_0 * - _S3107;
    float s_diff_dv_T_0 = dmax_1 * _S3107;
    DiffPair_float_0 _S3108;
    (&_S3108)->primal_0 = _S3085;
    (&_S3108)->differential_0 = 0.0f;
    DiffPair_float_0 _S3109;
    (&_S3109)->primal_0 = 0.00499999988824129f;
    (&_S3109)->differential_0 = 0.0f;
    DiffPair_float_0 _S3110;
    (&_S3110)->primal_0 = 0.98000001907348633f;
    (&_S3110)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3108, &_S3109, &_S3110, _S3106);
    float _S3111 = s_diff_dmax_T_1 / _S3083;
    float _S3112 = _S3081 * - _S3111;
    float _S3113 = _S3082 * _S3111;
    float2  _S3114 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3115;
    (&_S3115)->primal_0 = e2_5;
    (&_S3115)->differential_0 = _S3114;
    s_bwd_length_impl_1(&_S3115, _S3112);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3116;
    (&_S3116)->primal_0 = e1_9;
    (&_S3116)->differential_0 = _S3114;
    s_bwd_length_impl_1(&_S3116, _S3112);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3117;
    (&_S3117)->primal_0 = e0_9;
    (&_S3117)->differential_0 = _S3114;
    s_bwd_length_impl_1(&_S3117, _S3112);
    DiffPair_float_0 _S3118;
    (&_S3118)->primal_0 = _S3070;
    (&_S3118)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3118, _S3113);
    float _S3119 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S3120;
    (&_S3120)->primal_0 = _S3080;
    (&_S3120)->differential_0 = 0.0f;
    DiffPair_float_0 _S3121;
    (&_S3121)->primal_0 = _S3065;
    (&_S3121)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3120, &_S3121, _S3119);
    DiffPair_float_0 _S3122;
    (&_S3122)->primal_0 = _S3047;
    (&_S3122)->differential_0 = 0.0f;
    DiffPair_float_0 _S3123;
    (&_S3123)->primal_0 = _S3056;
    (&_S3123)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3122, &_S3123, _S3120.differential_0);
    DiffPair_float_0 _S3124;
    (&_S3124)->primal_0 = _S3079;
    (&_S3124)->differential_0 = 0.0f;
    DiffPair_float_0 _S3125;
    (&_S3125)->primal_0 = s2_0;
    (&_S3125)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3124, &_S3125, 0.0f);
    DiffPair_float_0 _S3126;
    (&_S3126)->primal_0 = s0_0;
    (&_S3126)->differential_0 = 0.0f;
    DiffPair_float_0 _S3127;
    (&_S3127)->primal_0 = s1_0;
    (&_S3127)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3126, &_S3127, _S3124.differential_0);
    float _S3128 = se_3 * _S3125.differential_0;
    float _S3129 = - _S3128;
    float _S3130 = _S3077 * _S3129;
    float _S3131 = _S3078 * _S3129;
    float _S3132 = _S3075 * _S3128;
    float _S3133 = _S3076 * _S3128;
    float _S3134 = se_3 * _S3127.differential_0;
    float _S3135 = - _S3134;
    float _S3136 = _S3069 * _S3135;
    float _S3137 = _S3067 * _S3134;
    float _S3138 = se_3 * _S3126.differential_0;
    float _S3139 = - _S3138;
    float _S3140 = _S3066 * _S3139;
    float _S3141 = _S3068 * _S3138;
    float _S3142 = - _S3118.differential_0;
    float _S3143 = _S3074 * _S3135 + _S3068 * _S3142;
    float _S3144 = _S3071 * _S3138 + _S3069 * _S3142;
    float _S3145 = _S3073 * _S3134 + _S3066 * _S3118.differential_0;
    float _S3146 = _S3072 * _S3139 + _S3067 * _S3118.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3147;
    (&_S3147)->primal_0 = _S3064;
    (&_S3147)->differential_0 = _S3114;
    s_bwd_length_impl_1(&_S3147, _S3121.differential_0);
    float2  _S3148 = - _S3147.differential_0;
    float2  _S3149 = e2_5 * _S3148;
    float2  _S3150 = _S3063 * _S3148;
    float _S3151 = _S3149.x + _S3149.y;
    DiffPair_float_0 _S3152;
    (&_S3152)->primal_0 = _S3060;
    (&_S3152)->differential_0 = 0.0f;
    DiffPair_float_0 _S3153;
    (&_S3153)->primal_0 = 0.0f;
    (&_S3153)->differential_0 = 0.0f;
    DiffPair_float_0 _S3154;
    (&_S3154)->primal_0 = 1.0f;
    (&_S3154)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3152, &_S3153, &_S3154, _S3151);
    float _S3155 = _S3152.differential_0 / _S3061;
    float _S3156 = _S3058 * - _S3155;
    float _S3157 = _S3059 * _S3155;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3158;
    (&_S3158)->primal_0 = e2_5;
    (&_S3158)->differential_0 = _S3114;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3159;
    (&_S3159)->primal_0 = e2_5;
    (&_S3159)->differential_0 = _S3114;
    s_bwd_prop_dot_1(&_S3158, &_S3159, _S3156);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3160;
    (&_S3160)->primal_0 = _S3057;
    (&_S3160)->differential_0 = _S3114;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3161;
    (&_S3161)->primal_0 = e2_5;
    (&_S3161)->differential_0 = _S3114;
    s_bwd_prop_dot_1(&_S3160, &_S3161, _S3157);
    float2  _S3162 = - (_S3147.differential_0 + _S3160.differential_0 + make_float2 (_S3133, _S3131));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3163;
    (&_S3163)->primal_0 = _S3055;
    (&_S3163)->differential_0 = _S3114;
    s_bwd_length_impl_1(&_S3163, _S3123.differential_0);
    float2  _S3164 = - _S3163.differential_0;
    float2  _S3165 = e1_9 * _S3164;
    float2  _S3166 = _S3054 * _S3164;
    float _S3167 = _S3165.x + _S3165.y;
    DiffPair_float_0 _S3168;
    (&_S3168)->primal_0 = _S3051;
    (&_S3168)->differential_0 = 0.0f;
    DiffPair_float_0 _S3169;
    (&_S3169)->primal_0 = 0.0f;
    (&_S3169)->differential_0 = 0.0f;
    DiffPair_float_0 _S3170;
    (&_S3170)->primal_0 = 1.0f;
    (&_S3170)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3168, &_S3169, &_S3170, _S3167);
    float _S3171 = _S3168.differential_0 / _S3052;
    float _S3172 = _S3049 * - _S3171;
    float _S3173 = _S3050 * _S3171;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3174;
    (&_S3174)->primal_0 = e1_9;
    (&_S3174)->differential_0 = _S3114;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3175;
    (&_S3175)->primal_0 = e1_9;
    (&_S3175)->differential_0 = _S3114;
    s_bwd_prop_dot_1(&_S3174, &_S3175, _S3172);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3176;
    (&_S3176)->primal_0 = _S3048;
    (&_S3176)->differential_0 = _S3114;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3177;
    (&_S3177)->primal_0 = e1_9;
    (&_S3177)->differential_0 = _S3114;
    s_bwd_prop_dot_1(&_S3176, &_S3177, _S3173);
    float2  _S3178 = - (_S3163.differential_0 + _S3176.differential_0 + make_float2 (_S3137, _S3136));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3179;
    (&_S3179)->primal_0 = _S3046;
    (&_S3179)->differential_0 = _S3114;
    s_bwd_length_impl_1(&_S3179, _S3122.differential_0);
    float2  _S3180 = - _S3179.differential_0;
    float2  _S3181 = e0_9 * _S3180;
    float2  _S3182 = _S3045 * _S3180;
    float _S3183 = _S3181.x + _S3181.y;
    DiffPair_float_0 _S3184;
    (&_S3184)->primal_0 = _S3042;
    (&_S3184)->differential_0 = 0.0f;
    DiffPair_float_0 _S3185;
    (&_S3185)->primal_0 = 0.0f;
    (&_S3185)->differential_0 = 0.0f;
    DiffPair_float_0 _S3186;
    (&_S3186)->primal_0 = 1.0f;
    (&_S3186)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3184, &_S3185, &_S3186, _S3183);
    float _S3187 = _S3184.differential_0 / _S3043;
    float _S3188 = _S3040 * - _S3187;
    float _S3189 = _S3041 * _S3187;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3190;
    (&_S3190)->primal_0 = e0_9;
    (&_S3190)->differential_0 = _S3114;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3191;
    (&_S3191)->primal_0 = e0_9;
    (&_S3191)->differential_0 = _S3114;
    s_bwd_prop_dot_1(&_S3190, &_S3191, _S3188);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3192;
    (&_S3192)->primal_0 = _S3039;
    (&_S3192)->differential_0 = _S3114;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3193;
    (&_S3193)->primal_0 = e0_9;
    (&_S3193)->differential_0 = _S3114;
    s_bwd_prop_dot_1(&_S3192, &_S3193, _S3189);
    float2  _S3194 = - (_S3179.differential_0 + _S3192.differential_0 + make_float2 (_S3141, _S3140));
    float2  _S3195 = _S3115.differential_0 + _S3150 + _S3159.differential_0 + _S3158.differential_0 + _S3161.differential_0 + make_float2 (_S3130, _S3132);
    float2  _S3196 = - _S3195;
    float2  _S3197 = _S3116.differential_0 + _S3166 + _S3175.differential_0 + _S3174.differential_0 + _S3177.differential_0 + make_float2 (_S3143, _S3145);
    float2  _S3198 = - _S3197;
    float2  _S3199 = _S3117.differential_0 + _S3182 + _S3191.differential_0 + _S3190.differential_0 + _S3193.differential_0 + make_float2 (_S3146, _S3144);
    float2  _S3200 = - _S3199;
    float2  _S3201 = make_float2 (_S3097, _S3108.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S3201;
    float2  _S3202 = _S3162 + _S3196 + _S3197;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S3202;
    float2  _S3203 = _S3178 + _S3198 + _S3199;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S3203;
    float2  _S3204 = _S3194 + _S3195 + _S3200;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S3204;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3205, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3206, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3207, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3208, float2  _S3209, float _S3210)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S3205, _S3206, _S3207, _S3208, _S3209, _S3210);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S3211 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S3211;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S3211;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S3211;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S3211;
    s_bwd_evaluate_alpha_opaque_triangle_precise_0(&dp_v0_1, &dp_v1_1, &dp_v2_1, &dp_hardness_1, p_5, v_alpha_2);
    *v_v0_1 = dp_v0_1.differential_0;
    *v_v1_1 = dp_v2_1.differential_0;
    *v_v2_1 = dp_v1_1.differential_0;
    *v_hardness_3 = dp_hardness_1.differential_0;
    return;
}

inline __device__ void evaluate_color_opaque_triangle(float2  v0_4, float2  v1_4, float2  v2_4, FixedArray<float3 , 3>  * colors_0, float3  depths_0, float2  p_6, float3  * color_6, float * depth_11)
{
    *color_6 = ((*colors_0)[int(0)] + (*colors_0)[int(1)] + (*colors_0)[int(2)]) / make_float3 (3.0f);
    *depth_11 = (depths_0.x + depths_0.y + depths_0.z) / 3.0f;
    return;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_2, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpcolors_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepths_0, float2  p_7, float3  dpcolor_0, float dpdepth_1)
{
    float _S3212 = 0.3333333432674408f * dpdepth_1;
    float3  _S3213 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S3214 = make_float3 (0.0f);
    float3  _S3215 = _S3214;
    *&((&_S3215)->z) = _S3212;
    *&((&_S3215)->y) = _S3212;
    *&((&_S3215)->x) = _S3212;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S3215;
    FixedArray<float3 , 3>  _S3216;
    _S3216[int(0)] = _S3214;
    _S3216[int(1)] = _S3214;
    _S3216[int(2)] = _S3214;
    _S3216[int(2)] = _S3213;
    _S3216[int(1)] = _S3213;
    _S3216[int(0)] = _S3213;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S3216;
    float2  _S3217 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S3217;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S3217;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S3217;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3218, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3219, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3220, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S3221, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3222, float2  _S3223, float3  _S3224, float _S3225)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S3218, _S3219, _S3220, _S3221, _S3222, _S3223, _S3224, _S3225);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_8, float3  v_color_0, float v_depth_6, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S3226 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S3226;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S3226;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S3226;
    float3  _S3227 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S3228 = { _S3227, _S3227, _S3227 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S3228;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S3227;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_8, v_color_0, v_depth_6);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

