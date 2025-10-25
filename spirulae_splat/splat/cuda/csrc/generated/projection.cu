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
        *depth_0 = length_1(mean_c_0);
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
        float _S148 = length_1(mean_c_1);
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
        *depth_2 = length_1(mean_c_2);
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
    *depth_3 = length_1(mean_c_3);
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
    *depth_4 = length_1(mean_c_4);
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
    *depth_5 = length_1(mean_c_5);
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

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_12, float _s_dOut_1)
{
    float _S239 = (*dpx_12).primal_0.x;
    float _S240 = (*dpx_12).primal_0.y;
    float _S241 = (*dpx_12).primal_0.z;
    DiffPair_float_0 _S242;
    (&_S242)->primal_0 = _S239 * _S239 + _S240 * _S240 + _S241 * _S241;
    (&_S242)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S242, _s_dOut_1);
    float _S243 = (*dpx_12).primal_0.z * _S242.differential_0;
    float _S244 = _S243 + _S243;
    float _S245 = (*dpx_12).primal_0.y * _S242.differential_0;
    float _S246 = _S245 + _S245;
    float _S247 = (*dpx_12).primal_0.x * _S242.differential_0;
    float _S248 = _S247 + _S247;
    float3  _S249 = make_float3 (0.0f);
    *&((&_S249)->z) = _S244;
    *&((&_S249)->y) = _S246;
    *&((&_S249)->x) = _S248;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S249;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S250, float _S251)
{
    s_bwd_prop_length_impl_1(_S250, _S251);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S252, DiffPair_float_0 * _S253, float _S254)
{
    _d_min_0(_S252, _S253, _S254);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S255, float _S256)
{
    _d_log_0(_S255, _S256);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S257, float _S258)
{
    _d_exp_0(_S257, _S258);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S259, DiffPair_float_0 * _S260, float _S261)
{
    _d_max_0(_S259, _S260, _S261);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S262, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S263, Matrix<float, 2, 2>  _S264)
{
    mul_3(_S262, _S263, _S264);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S265, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S266, Matrix<float, 2, 3>  _S267)
{
    mul_2(_S265, _S266, _S267);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S268, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S269, Matrix<float, 3, 3>  _S270)
{
    mul_1(_S268, _S269, _S270);
    return;
}

inline __device__ void s_bwd_prop_rsqrt_0(DiffPair_float_0 * _S271, float _S272)
{
    _d_rsqrt_0(_S271, _S272);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S273, float3  _S274)
{
    _d_exp_vector_0(_S273, _S274);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_10, float fy_10, float cx_10, float cy_10, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, uint image_width_6, uint image_height_6, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    float3  mean_c_6 = s_primal_ctx_mul_0(R_8, mean_6) + t_7;
    float3  _S275 = s_primal_ctx_exp_0(scale_8);
    float _S276 = quat_9.y;
    float _S277 = _S276 * _S276 + quat_9.z * quat_9.z + quat_9.w * quat_9.w + quat_9.x * quat_9.x;
    float _S278 = s_primal_ctx_rsqrt_0(_S277);
    float x_33 = quat_9.y * _S278;
    float y_18 = quat_9.z * _S278;
    float z_15 = quat_9.w * _S278;
    float w_9 = quat_9.x * _S278;
    float x2_9 = x_33 * x_33;
    float y2_9 = y_18 * y_18;
    float z2_15 = z_15 * z_15;
    float xy_9 = x_33 * y_18;
    float xz_9 = x_33 * z_15;
    float yz_9 = y_18 * z_15;
    float wx_9 = w_9 * x_33;
    float wy_9 = w_9 * y_18;
    float wz_9 = w_9 * z_15;
    Matrix<float, 3, 3>  _S279 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S275.x, 0.0f, 0.0f, 0.0f, _S275.y, 0.0f, 0.0f, 0.0f, _S275.z);
    Matrix<float, 3, 3>  _S280 = s_primal_ctx_mul_1(_S279, S_0);
    Matrix<float, 3, 3>  _S281 = transpose_0(_S280);
    Matrix<float, 3, 3>  _S282 = s_primal_ctx_mul_1(_S280, _S281);
    Matrix<float, 3, 3>  _S283 = s_primal_ctx_mul_1(R_8, _S282);
    Matrix<float, 3, 3>  _S284 = transpose_0(R_8);
    Matrix<float, 3, 3>  _S285 = s_primal_ctx_mul_1(_S283, _S284);
    float _S286 = float(image_width_6);
    float _S287 = float(image_height_6);
    float _S288 = 0.30000001192092896f * (0.5f * _S286 / fx_10);
    float lim_x_pos_0 = (_S286 - cx_10) / fx_10 + _S288;
    float _S289 = 0.30000001192092896f * (0.5f * _S287 / fy_10);
    float lim_y_pos_0 = (_S287 - cy_10) / fy_10 + _S289;
    float rz_4 = 1.0f / mean_c_6.z;
    float _S290 = mean_c_6.z * mean_c_6.z;
    float rz2_4 = rz_4 * rz_4;
    float _S291 = - (cx_10 / fx_10 + _S288);
    float _S292 = mean_c_6.x * rz_4;
    float _S293 = s_primal_ctx_max_0(_S291, _S292);
    float _S294 = s_primal_ctx_min_0(lim_x_pos_0, _S293);
    float _S295 = - (cy_10 / fy_10 + _S289);
    float _S296 = mean_c_6.y * rz_4;
    float _S297 = s_primal_ctx_max_0(_S295, _S296);
    float _S298 = s_primal_ctx_min_0(lim_y_pos_0, _S297);
    float _S299 = - fx_10;
    float _S300 = _S299 * (mean_c_6.z * _S294);
    float _S301 = - fy_10;
    float _S302 = _S301 * (mean_c_6.z * _S298);
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_10 * rz_4, 0.0f, _S300 * rz2_4, 0.0f, fy_10 * rz_4, _S302 * rz2_4);
    Matrix<float, 2, 3>  _S303 = s_primal_ctx_mul_2(J_8, _S285);
    Matrix<float, 3, 2>  _S304 = transpose_1(J_8);
    Matrix<float, 2, 2>  _S305 = s_primal_ctx_mul_3(_S303, _S304);
    float _S306 = fx_10 * mean_c_6.x;
    float _S307 = fy_10 * mean_c_6.y;
    float _S308 = _S305.rows[int(0)].y * _S305.rows[int(1)].x;
    float det_orig_7 = _S305.rows[int(0)].x * _S305.rows[int(1)].y - _S308;
    float _S309 = _S305.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S310 = _S305;
    *&(((&_S310)->rows + (int(0)))->x) = _S309;
    float _S311 = _S305.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S310)->rows + (int(1)))->y) = _S311;
    Matrix<float, 2, 2>  _S312 = _S310;
    Matrix<float, 2, 2>  _S313 = _S310;
    float det_blur_4 = _S309 * _S311 - _S308;
    float _S314 = det_orig_7 / det_blur_4;
    float _S315 = det_blur_4 * det_blur_4;
    float _S316 = s_primal_ctx_max_0(0.0f, _S314);
    float _S317 = s_primal_ctx_sqrt_0(_S316);
    float invdet_7 = 1.0f / det_blur_4;
    float _S318 = - _S305.rows[int(0)].y;
    float _S319 = - _S305.rows[int(1)].x;
    float _S320 = - in_opacity_6;
    float _S321 = 1.0f + s_primal_ctx_exp_1(_S320);
    float _S322 = 1.0f / _S321;
    float _S323 = _S321 * _S321;
    float _S324;
    if(antialiased_6)
    {
        _S324 = _S322 * _S317;
    }
    else
    {
        _S324 = _S322;
    }
    float _S325 = _S324 / 0.00392156885936856f;
    float _S326 = 2.0f * s_primal_ctx_log_0(_S325);
    float _S327 = s_primal_ctx_sqrt_0(_S326);
    float _S328 = _S312.rows[int(0)].x;
    float _S329 = _S313.rows[int(1)].y;
    float3  _S330 = mean_6 - - s_primal_ctx_mul_0(_S284, t_7);
    float _S331 = _S330.x;
    float _S332 = _S330.y;
    float _S333 = _S330.z;
    float _S334 = _S331 * _S331 + _S332 * _S332 + _S333 * _S333;
    float _S335 = s_primal_ctx_sqrt_0(_S334);
    float x_34 = _S331 / _S335;
    float3  _S336 = make_float3 (x_34);
    float _S337 = _S335 * _S335;
    float y_19 = _S332 / _S335;
    float z_16 = _S333 / _S335;
    float3  _S338 = make_float3 (z_16);
    float _S339 = - y_19;
    float3  _S340 = make_float3 (_S339);
    float z2_16 = z_16 * z_16;
    float fTmp0B_6 = -1.09254848957061768f * z_16;
    float fC1_6 = x_34 * x_34 - y_19 * y_19;
    float _S341 = 2.0f * x_34;
    float fS1_6 = _S341 * y_19;
    float pSH6_0 = 0.94617468118667603f * z2_16 - 0.31539157032966614f;
    float3  _S342 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_6 * x_34;
    float3  _S343 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_6 * y_19;
    float3  _S344 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_6;
    float3  _S345 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_6;
    float3  _S346 = make_float3 (pSH4_0);
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_16;
    float _S347 = 1.86588168144226074f * z2_16 - 1.11952900886535645f;
    float pSH12_0 = z_16 * _S347;
    float3  _S348 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_6 * x_34;
    float3  _S349 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_6 * y_19;
    float3  _S350 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_6 * fC1_6;
    float3  _S351 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_6 * fS1_6;
    float3  _S352 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_34 * fC1_6 - y_19 * fS1_6);
    float3  _S353 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_34 * fS1_6 + y_19 * fC1_6);
    float3  _S354 = make_float3 (pSH9_0);
    float3  _S355 = make_float3 (0.0f);
    float3  _S356 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S357;
    (&_S357)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S339) * (*sh_coeffs_6)[int(1)] + make_float3 (z_16) * (*sh_coeffs_6)[int(2)] - make_float3 (x_34) * (*sh_coeffs_6)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_6)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_6)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_6)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_6)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_6)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_6)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_6)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_6)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_6)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_6)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_6)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f);
    (&_S357)->differential_0 = _S356;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S358;
    (&_S358)->primal_0 = _S355;
    (&_S358)->differential_0 = _S356;
    s_bwd_prop_max_0(&_S357, &_S358, v_rgb_0);
    float3  _S359 = _S353 * _S357.differential_0;
    float3  _S360 = (*sh_coeffs_6)[int(15)] * _S357.differential_0;
    float3  _S361 = _S351 * _S357.differential_0;
    float3  _S362 = (*sh_coeffs_6)[int(14)] * _S357.differential_0;
    float3  _S363 = _S349 * _S357.differential_0;
    float3  _S364 = (*sh_coeffs_6)[int(13)] * _S357.differential_0;
    float3  _S365 = _S348 * _S357.differential_0;
    float3  _S366 = (*sh_coeffs_6)[int(12)] * _S357.differential_0;
    float3  _S367 = _S350 * _S357.differential_0;
    float3  _S368 = (*sh_coeffs_6)[int(11)] * _S357.differential_0;
    float3  _S369 = _S352 * _S357.differential_0;
    float3  _S370 = (*sh_coeffs_6)[int(10)] * _S357.differential_0;
    float3  _S371 = _S354 * _S357.differential_0;
    float3  _S372 = (*sh_coeffs_6)[int(9)] * _S357.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S372.x + _S372.y + _S372.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S360.x + _S360.y + _S360.z);
    float _S373 = _S370.x + _S370.y + _S370.z;
    float _S374 = _S362.x + _S362.y + _S362.z;
    float _S375 = _S368.x + _S368.y + _S368.z;
    float _S376 = _S364.x + _S364.y + _S364.z;
    float _S377 = _S366.x + _S366.y + _S366.z;
    float _S378 = - s_diff_fC2_T_0;
    float3  _S379 = _S345 * _S357.differential_0;
    float3  _S380 = (*sh_coeffs_6)[int(8)] * _S357.differential_0;
    float3  _S381 = _S343 * _S357.differential_0;
    float3  _S382 = (*sh_coeffs_6)[int(7)] * _S357.differential_0;
    float3  _S383 = _S342 * _S357.differential_0;
    float3  _S384 = (*sh_coeffs_6)[int(6)] * _S357.differential_0;
    float3  _S385 = _S344 * _S357.differential_0;
    float3  _S386 = (*sh_coeffs_6)[int(5)] * _S357.differential_0;
    float3  _S387 = _S346 * _S357.differential_0;
    float3  _S388 = (*sh_coeffs_6)[int(4)] * _S357.differential_0;
    float _S389 = _S386.x + _S386.y + _S386.z;
    float _S390 = _S382.x + _S382.y + _S382.z;
    float _S391 = fTmp1B_6 * _S373 + x_34 * s_diff_fS2_T_0 + y_19 * _S378 + 0.54627424478530884f * (_S388.x + _S388.y + _S388.z);
    float _S392 = fTmp1B_6 * _S374 + y_19 * s_diff_fS2_T_0 + x_34 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S380.x + _S380.y + _S380.z);
    float _S393 = y_19 * - _S392;
    float _S394 = x_34 * _S392;
    float _S395 = z_16 * (1.86588168144226074f * (z_16 * _S377) + -2.28522896766662598f * (y_19 * _S375 + x_34 * _S376) + 0.94617468118667603f * (_S384.x + _S384.y + _S384.z));
    float3  _S396 = make_float3 (0.48860251903533936f) * _S357.differential_0;
    float3  _S397 = - _S396;
    float3  _S398 = _S336 * _S397;
    float3  _S399 = (*sh_coeffs_6)[int(3)] * _S397;
    float3  _S400 = _S338 * _S396;
    float3  _S401 = (*sh_coeffs_6)[int(2)] * _S396;
    float3  _S402 = _S340 * _S396;
    float3  _S403 = (*sh_coeffs_6)[int(1)] * _S396;
    float _S404 = (_S347 * _S377 + 1.44530570507049561f * (fS1_6 * _S373 + fC1_6 * _S374) + -1.09254848957061768f * (y_19 * _S389 + x_34 * _S390) + _S395 + _S395 + _S401.x + _S401.y + _S401.z) / _S337;
    float _S405 = _S335 * _S404;
    float _S406 = (fTmp0C_6 * _S375 + fC1_6 * s_diff_fS2_T_0 + fS1_6 * _S378 + fTmp0B_6 * _S389 + _S341 * _S391 + _S393 + _S393 + - (_S403.x + _S403.y + _S403.z)) / _S337;
    float _S407 = _S335 * _S406;
    float _S408 = (fTmp0C_6 * _S376 + fS1_6 * s_diff_fS2_T_0 + fC1_6 * s_diff_fC2_T_0 + fTmp0B_6 * _S390 + 2.0f * (y_19 * _S391) + _S394 + _S394 + _S399.x + _S399.y + _S399.z) / _S337;
    float _S409 = _S335 * _S408;
    float _S410 = _S333 * - _S404 + _S332 * - _S406 + _S331 * - _S408;
    DiffPair_float_0 _S411;
    (&_S411)->primal_0 = _S334;
    (&_S411)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S411, _S410);
    float _S412 = _S333 * _S411.differential_0;
    float _S413 = _S332 * _S411.differential_0;
    float _S414 = _S331 * _S411.differential_0;
    float3  _S415 = make_float3 (0.282094806432724f) * _S357.differential_0;
    float3  _S416 = make_float3 (_S409 + _S414 + _S414, _S407 + _S413 + _S413, _S405 + _S412 + _S412);
    float3  _S417 = - - _S416;
    Matrix<float, 3, 3>  _S418 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S419;
    (&_S419)->primal_0 = _S284;
    (&_S419)->differential_0 = _S418;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S420;
    (&_S420)->primal_0 = t_7;
    (&_S420)->differential_0 = _S356;
    s_bwd_prop_mul_0(&_S419, &_S420, _S417);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S421 = _S419;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S422 = _S420;
    float2  _S423 = make_float2 (0.0f);
    float2  _S424 = _S423;
    *&((&_S424)->y) = v_conic_0.z;
    float2  _S425 = _S423;
    *&((&_S425)->y) = v_conic_0.y;
    *&((&_S425)->x) = v_conic_0.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S426;
    (&_S426)->primal_0 = mean_c_6;
    (&_S426)->differential_0 = _S356;
    s_bwd_length_impl_1(&_S426, v_depth_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S427 = _S426;
    DiffPair_float_0 _S428;
    (&_S428)->primal_0 = _S329;
    (&_S428)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S428, 0.0f);
    DiffPair_float_0 _S429;
    (&_S429)->primal_0 = _S328;
    (&_S429)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S429, 0.0f);
    DiffPair_float_0 _S430;
    (&_S430)->primal_0 = 3.32999992370605469f;
    (&_S430)->differential_0 = 0.0f;
    DiffPair_float_0 _S431;
    (&_S431)->primal_0 = _S327;
    (&_S431)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S430, &_S431, 0.0f);
    DiffPair_float_0 _S432;
    (&_S432)->primal_0 = _S326;
    (&_S432)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S432, _S431.differential_0);
    float _S433 = 2.0f * _S432.differential_0;
    DiffPair_float_0 _S434;
    (&_S434)->primal_0 = _S325;
    (&_S434)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S434, _S433);
    float _S435 = v_opacity_0 + 254.9999847412109375f * _S434.differential_0;
    Matrix<float, 2, 2>  _S436 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S437 = _S436;
    _S437[int(1)] = _S424;
    _S437[int(0)] = _S425;
    Matrix<float, 2, 2>  _S438 = _S437;
    FixedArray<float3 , 16>  _S439;
    _S439[int(0)] = _S356;
    _S439[int(1)] = _S356;
    _S439[int(2)] = _S356;
    _S439[int(3)] = _S356;
    _S439[int(4)] = _S356;
    _S439[int(5)] = _S356;
    _S439[int(6)] = _S356;
    _S439[int(7)] = _S356;
    _S439[int(8)] = _S356;
    _S439[int(9)] = _S356;
    _S439[int(10)] = _S356;
    _S439[int(11)] = _S356;
    _S439[int(12)] = _S356;
    _S439[int(13)] = _S356;
    _S439[int(14)] = _S356;
    _S439[int(15)] = _S356;
    _S439[int(7)] = _S381;
    _S439[int(0)] = _S415;
    _S439[int(1)] = _S402;
    _S439[int(2)] = _S400;
    _S439[int(3)] = _S398;
    _S439[int(4)] = _S387;
    _S439[int(5)] = _S385;
    _S439[int(6)] = _S383;
    _S439[int(15)] = _S359;
    _S439[int(8)] = _S379;
    _S439[int(9)] = _S371;
    _S439[int(10)] = _S369;
    _S439[int(11)] = _S367;
    _S439[int(12)] = _S365;
    _S439[int(13)] = _S363;
    _S439[int(14)] = _S361;
    float3  _S440 = _S439[int(0)];
    float3  _S441 = _S439[int(1)];
    float3  _S442 = _S439[int(2)];
    float3  _S443 = _S439[int(3)];
    float3  _S444 = _S439[int(4)];
    float3  _S445 = _S439[int(5)];
    float3  _S446 = _S439[int(6)];
    float3  _S447 = _S439[int(7)];
    float3  _S448 = _S439[int(8)];
    float3  _S449 = _S439[int(9)];
    float3  _S450 = _S439[int(10)];
    float3  _S451 = _S439[int(11)];
    float3  _S452 = _S439[int(12)];
    float3  _S453 = _S439[int(13)];
    float3  _S454 = _S439[int(14)];
    float3  _S455 = _S439[int(15)];
    float2  _S456 = make_float2 (0.0f, _S428.differential_0);
    float2  _S457 = make_float2 (_S429.differential_0, 0.0f);
    float _S458;
    if(antialiased_6)
    {
        float _S459 = _S322 * _S435;
        _S324 = _S317 * _S435;
        _S458 = _S459;
    }
    else
    {
        _S324 = _S435;
        _S458 = 0.0f;
    }
    float _S460 = - (_S324 / _S323);
    DiffPair_float_0 _S461;
    (&_S461)->primal_0 = _S320;
    (&_S461)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S461, _S460);
    float _S462 = - _S461.differential_0;
    float _S463 = invdet_7 * _S438.rows[int(1)].y;
    float _S464 = - (invdet_7 * _S438.rows[int(1)].x);
    float _S465 = - (invdet_7 * _S438.rows[int(0)].y);
    float _S466 = invdet_7 * _S438.rows[int(0)].x;
    float _S467 = - ((_S309 * _S438.rows[int(1)].y + _S319 * _S438.rows[int(1)].x + _S318 * _S438.rows[int(0)].y + _S311 * _S438.rows[int(0)].x) / _S315);
    DiffPair_float_0 _S468;
    (&_S468)->primal_0 = _S316;
    (&_S468)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S468, _S458);
    DiffPair_float_0 _S469;
    (&_S469)->primal_0 = 0.0f;
    (&_S469)->differential_0 = 0.0f;
    DiffPair_float_0 _S470;
    (&_S470)->primal_0 = _S314;
    (&_S470)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S469, &_S470, _S468.differential_0);
    float _S471 = _S470.differential_0 / _S315;
    float s_diff_det_orig_T_0 = det_blur_4 * _S471;
    float _S472 = _S467 + det_orig_7 * - _S471;
    float _S473 = - _S472;
    float _S474 = _S309 * _S472;
    float _S475 = _S311 * _S472;
    Matrix<float, 2, 2>  _S476 = _S436;
    _S476[int(1)] = _S456;
    _S476[int(0)] = _S457;
    _S310 = _S476;
    *&(((&_S310)->rows + (int(1)))->y) = 0.0f;
    float _S477 = _S466 + _S474 + _S476.rows[int(1)].y;
    *&(((&_S310)->rows + (int(0)))->x) = 0.0f;
    float _S478 = _S463 + _S475 + _S476.rows[int(0)].x;
    float _S479 = _S473 + - s_diff_det_orig_T_0;
    float _S480 = _S464 + _S305.rows[int(0)].y * _S479;
    float _S481 = _S465 + _S305.rows[int(1)].x * _S479;
    float _S482 = _S305.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S483 = _S477 + _S305.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S484 = _S423;
    *&((&_S484)->x) = _S480;
    *&((&_S484)->y) = _S483;
    float _S485 = _S478 + _S482;
    float2  _S486 = _S423;
    *&((&_S486)->y) = _S481;
    *&((&_S486)->x) = _S485;
    float _S487 = _S307 * v_mean2d_0.y;
    float _S488 = fy_10 * (rz_4 * v_mean2d_0.y);
    float _S489 = _S306 * v_mean2d_0.x;
    float _S490 = fx_10 * (rz_4 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S491 = _S436;
    _S491[int(1)] = _S484;
    _S491[int(0)] = _S486;
    Matrix<float, 2, 2>  _S492 = _S310 + _S491;
    Matrix<float, 2, 3>  _S493 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S494;
    (&_S494)->primal_0 = _S303;
    (&_S494)->differential_0 = _S493;
    Matrix<float, 3, 2>  _S495 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S496;
    (&_S496)->primal_0 = _S304;
    (&_S496)->differential_0 = _S495;
    s_bwd_prop_mul_1(&_S494, &_S496, _S492);
    Matrix<float, 2, 3>  _S497 = transpose_2(_S496.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S498;
    (&_S498)->primal_0 = J_8;
    (&_S498)->differential_0 = _S493;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S499;
    (&_S499)->primal_0 = _S285;
    (&_S499)->differential_0 = _S418;
    s_bwd_prop_mul_2(&_S498, &_S499, _S494.differential_0);
    Matrix<float, 2, 3>  _S500 = _S497 + _S498.differential_0;
    float _S501 = _S302 * _S500.rows[int(1)].z;
    float s_diff_ty_T_0 = _S301 * (rz2_4 * _S500.rows[int(1)].z);
    float _S502 = fy_10 * _S500.rows[int(1)].y;
    float _S503 = _S300 * _S500.rows[int(0)].z;
    float s_diff_tx_T_0 = _S299 * (rz2_4 * _S500.rows[int(0)].z);
    float _S504 = fx_10 * _S500.rows[int(0)].x;
    float _S505 = mean_c_6.z * s_diff_ty_T_0;
    float _S506 = _S298 * s_diff_ty_T_0;
    DiffPair_float_0 _S507;
    (&_S507)->primal_0 = lim_y_pos_0;
    (&_S507)->differential_0 = 0.0f;
    DiffPair_float_0 _S508;
    (&_S508)->primal_0 = _S297;
    (&_S508)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S507, &_S508, _S505);
    DiffPair_float_0 _S509;
    (&_S509)->primal_0 = _S295;
    (&_S509)->differential_0 = 0.0f;
    DiffPair_float_0 _S510;
    (&_S510)->primal_0 = _S296;
    (&_S510)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S509, &_S510, _S508.differential_0);
    float _S511 = mean_c_6.y * _S510.differential_0;
    float _S512 = rz_4 * _S510.differential_0;
    float _S513 = mean_c_6.z * s_diff_tx_T_0;
    float _S514 = _S294 * s_diff_tx_T_0;
    DiffPair_float_0 _S515;
    (&_S515)->primal_0 = lim_x_pos_0;
    (&_S515)->differential_0 = 0.0f;
    DiffPair_float_0 _S516;
    (&_S516)->primal_0 = _S293;
    (&_S516)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S515, &_S516, _S513);
    DiffPair_float_0 _S517;
    (&_S517)->primal_0 = _S291;
    (&_S517)->differential_0 = 0.0f;
    DiffPair_float_0 _S518;
    (&_S518)->primal_0 = _S292;
    (&_S518)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S517, &_S518, _S516.differential_0);
    float _S519 = rz_4 * (_S501 + _S503);
    float _S520 = _S506 + _S514 + - ((_S487 + _S489 + _S502 + _S504 + _S511 + mean_c_6.x * _S518.differential_0 + _S519 + _S519) / _S290);
    float _S521 = _S488 + _S512;
    float _S522 = _S490 + rz_4 * _S518.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S523;
    (&_S523)->primal_0 = _S283;
    (&_S523)->differential_0 = _S418;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S524;
    (&_S524)->primal_0 = _S284;
    (&_S524)->differential_0 = _S418;
    s_bwd_prop_mul_3(&_S523, &_S524, _S499.differential_0);
    Matrix<float, 3, 3>  _S525 = transpose_0(_S524.differential_0 + _S421.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S526;
    (&_S526)->primal_0 = R_8;
    (&_S526)->differential_0 = _S418;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S527;
    (&_S527)->primal_0 = _S282;
    (&_S527)->differential_0 = _S418;
    s_bwd_prop_mul_3(&_S526, &_S527, _S523.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S528;
    (&_S528)->primal_0 = _S280;
    (&_S528)->differential_0 = _S418;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S529;
    (&_S529)->primal_0 = _S281;
    (&_S529)->differential_0 = _S418;
    s_bwd_prop_mul_3(&_S528, &_S529, _S527.differential_0);
    Matrix<float, 3, 3>  _S530 = _S528.differential_0 + transpose_0(_S529.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S531;
    (&_S531)->primal_0 = _S279;
    (&_S531)->differential_0 = _S418;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S532;
    (&_S532)->primal_0 = S_0;
    (&_S532)->differential_0 = _S418;
    s_bwd_prop_mul_3(&_S531, &_S532, _S530);
    Matrix<float, 3, 3>  _S533 = transpose_0(_S531.differential_0);
    float _S534 = 2.0f * - _S533.rows[int(2)].z;
    float _S535 = 2.0f * _S533.rows[int(2)].y;
    float _S536 = 2.0f * _S533.rows[int(2)].x;
    float _S537 = 2.0f * _S533.rows[int(1)].z;
    float _S538 = 2.0f * - _S533.rows[int(1)].y;
    float _S539 = 2.0f * _S533.rows[int(1)].x;
    float _S540 = 2.0f * _S533.rows[int(0)].z;
    float _S541 = 2.0f * _S533.rows[int(0)].y;
    float _S542 = 2.0f * - _S533.rows[int(0)].x;
    float _S543 = - _S539 + _S541;
    float _S544 = _S536 + - _S540;
    float _S545 = - _S535 + _S537;
    float _S546 = _S535 + _S537;
    float _S547 = _S536 + _S540;
    float _S548 = _S539 + _S541;
    float _S549 = z_15 * (_S538 + _S542);
    float _S550 = y_18 * (_S534 + _S542);
    float _S551 = x_33 * (_S534 + _S538);
    float _S552 = z_15 * _S543 + y_18 * _S544 + x_33 * _S545;
    float _S553 = _S278 * _S552;
    float _S554 = w_9 * _S543 + y_18 * _S546 + x_33 * _S547 + _S549 + _S549;
    float _S555 = _S278 * _S554;
    float _S556 = w_9 * _S544 + z_15 * _S546 + x_33 * _S548 + _S550 + _S550;
    float _S557 = _S278 * _S556;
    float _S558 = w_9 * _S545 + z_15 * _S547 + y_18 * _S548 + _S551 + _S551;
    float _S559 = _S278 * _S558;
    float _S560 = quat_9.x * _S552 + quat_9.w * _S554 + quat_9.z * _S556 + quat_9.y * _S558;
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = _S277;
    (&_S561)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S561, _S560);
    float _S562 = quat_9.x * _S561.differential_0;
    float _S563 = quat_9.w * _S561.differential_0;
    float _S564 = quat_9.z * _S561.differential_0;
    float _S565 = quat_9.y * _S561.differential_0;
    float _S566 = _S555 + _S563 + _S563;
    float _S567 = _S557 + _S564 + _S564;
    float _S568 = _S559 + _S565 + _S565;
    float _S569 = _S553 + _S562 + _S562;
    float3  _S570 = _S356;
    *&((&_S570)->z) = _S532.differential_0.rows[int(2)].z;
    *&((&_S570)->y) = _S532.differential_0.rows[int(1)].y;
    *&((&_S570)->x) = _S532.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S571;
    (&_S571)->primal_0 = scale_8;
    (&_S571)->differential_0 = _S356;
    s_bwd_prop_exp_1(&_S571, _S570);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S572 = _S571;
    float3  _S573 = _S356;
    *&((&_S573)->z) = _S520;
    *&((&_S573)->y) = _S521;
    *&((&_S573)->x) = _S522;
    float3  _S574 = _S427.differential_0 + _S573;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S575;
    (&_S575)->primal_0 = R_8;
    (&_S575)->differential_0 = _S418;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S576;
    (&_S576)->primal_0 = mean_6;
    (&_S576)->differential_0 = _S356;
    s_bwd_prop_mul_0(&_S575, &_S576, _S574);
    float3  _S577 = _S574 + _S422.differential_0;
    Matrix<float, 3, 3>  _S578 = _S525 + _S526.differential_0 + _S575.differential_0;
    float4  _S579 = make_float4 (0.0f);
    *&((&_S579)->w) = _S566;
    *&((&_S579)->z) = _S567;
    *&((&_S579)->y) = _S568;
    *&((&_S579)->x) = _S569;
    float4  _S580 = _S579;
    float3  _S581 = _S576.differential_0 + _S416;
    *v_mean_0 = _S581;
    *v_quat_0 = _S580;
    *v_scale_0 = _S572.differential_0;
    *v_in_opacity_0 = _S462;
    (*v_sh_coeffs_0)[int(0)] = _S440;
    (*v_sh_coeffs_0)[int(1)] = _S441;
    (*v_sh_coeffs_0)[int(2)] = _S442;
    (*v_sh_coeffs_0)[int(3)] = _S443;
    (*v_sh_coeffs_0)[int(4)] = _S444;
    (*v_sh_coeffs_0)[int(5)] = _S445;
    (*v_sh_coeffs_0)[int(6)] = _S446;
    (*v_sh_coeffs_0)[int(7)] = _S447;
    (*v_sh_coeffs_0)[int(8)] = _S448;
    (*v_sh_coeffs_0)[int(9)] = _S449;
    (*v_sh_coeffs_0)[int(10)] = _S450;
    (*v_sh_coeffs_0)[int(11)] = _S451;
    (*v_sh_coeffs_0)[int(12)] = _S452;
    (*v_sh_coeffs_0)[int(13)] = _S453;
    (*v_sh_coeffs_0)[int(14)] = _S454;
    (*v_sh_coeffs_0)[int(15)] = _S455;
    *v_R_0 = _S578;
    *v_t_0 = _S577;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S582;
    DiffPair_float_0 _S583;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S584;
    DiffPair_float_0 _S585;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S586;
    DiffPair_float_0 _S587;
    DiffPair_float_0 _S588;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S589;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S590;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S591;
};

inline __device__ CameraDistortion_0 s_primal_ctx_CameraDistortion_x24init_0(float4  dpradial_coeffs_0, float2  dptangential_coeffs_0, float2  dpthin_prism_coeffs_0)
{
    CameraDistortion_0 _S592 = { dpradial_coeffs_0, dptangential_coeffs_0, dpthin_prism_coeffs_0 };
    return _S592;
}

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S593, float _S594)
{
    return s_primal_ctx_atan2_0(_S593, _S594);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S595;
    DiffPair_float_0 _S596;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_0)
{
    DiffPair_float_0 _S597 = { 0.0f, 0.0f };
    _s_diff_ctx_0->_S595 = _S597;
    _s_diff_ctx_0->_S596 = _S597;
    (&_s_diff_ctx_0->_S595)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S595)->differential_0 = 0.0f;
    (&_s_diff_ctx_0->_S596)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S596)->differential_0 = 0.0f;
    DiffPair_float_0 _S598 = *dpdpy_0;
    _s_diff_ctx_0->_S595 = *dpdpy_0;
    DiffPair_float_0 _S599 = *dpdpx_0;
    _s_diff_ctx_0->_S596 = *dpdpx_0;
    float _S600 = _S599.primal_0 * _S599.primal_0 + _S598.primal_0 * _S598.primal_0;
    float _S601 = - _S598.primal_0 / _S600 * dpdOut_0;
    float _S602 = _S599.primal_0 / _S600 * dpdOut_0;
    dpdpy_0->primal_0 = _S598.primal_0;
    dpdpy_0->differential_0 = _S602;
    dpdpx_0->primal_0 = _S599.primal_0;
    dpdpx_0->differential_0 = _S601;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S603, DiffPair_float_0 * _S604, float _S605, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_float_0 _S606 = { 0.0f, 0.0f };
    _s_diff_ctx_1->_S582 = _S606;
    _s_diff_ctx_1->_S583 = _S606;
    (&_s_diff_ctx_1->_S582)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S582)->differential_0 = 0.0f;
    (&_s_diff_ctx_1->_S583)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S583)->differential_0 = 0.0f;
    DiffPair_float_0 _S607 = *_S603;
    _s_diff_ctx_1->_S582 = *_S603;
    DiffPair_float_0 _S608 = *_S604;
    _s_diff_ctx_1->_S583 = *_S604;
    DiffPair_float_0 _S609 = _S607;
    DiffPair_float_0 _S610 = _S608;
    s_bwd_prop_d_atan2_Intermediates_0 _S611;
    (&_S611)->_S595 = _S606;
    (&_S611)->_S596 = _S606;
    s_primal_ctx_d_atan2_0(&_S609, &_S610, _S605, &_S611);
    *_S603 = _S609;
    *_S604 = _S610;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S612;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S613;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S614;
    DiffPair_float_0 _S615;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S616;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S617;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S618 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S617 = _S618;
    (&_s_diff_ctx_2->_S617)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S617)->differential_0 = 0.0f;
    DiffPair_float_0 _S619 = *dpdpx_1;
    _s_diff_ctx_2->_S617 = *dpdpx_1;
    float _S620 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S619.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S619.primal_0;
    dpdpx_1->differential_0 = _S620;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S621, float _S622, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S623 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S613 = _S623;
    (&_s_diff_ctx_3->_S613)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S613)->differential_0 = 0.0f;
    DiffPair_float_0 _S624 = *_S621;
    _s_diff_ctx_3->_S613 = *_S621;
    DiffPair_float_0 _S625 = _S624;
    s_bwd_prop_d_sqrt_Intermediates_0 _S626;
    (&_S626)->_S617 = _S623;
    s_primal_ctx_d_sqrt_0(&_S625, _S622, &_S626);
    *_S621 = _S625;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_4)
{
    float2  _S627 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S628 = { _S627, _S627 };
    DiffPair_float_0 _S629 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S630 = { _S629 };
    _s_diff_ctx_4->_S614 = _S628;
    _s_diff_ctx_4->_S615 = _S629;
    _s_diff_ctx_4->_S616 = _S630;
    (&_s_diff_ctx_4->_S614)->primal_0 = _S627;
    (&_s_diff_ctx_4->_S614)->differential_0 = _S627;
    (&_s_diff_ctx_4->_S615)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S615)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S631 = *dpdpx_2;
    _s_diff_ctx_4->_S614 = *dpdpx_2;
    float _S632 = _S631.primal_0.x;
    float _S633 = _S631.primal_0.y;
    DiffPair_float_0 _S634;
    (&_S634)->primal_0 = _S632 * _S632 + _S633 * _S633;
    (&_S634)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S634, dp_s_dOut_0, &_s_diff_ctx_4->_S616);
    _s_diff_ctx_4->_S615 = _S634;
    float _S635 = _S631.primal_0.y * _S634.differential_0;
    float _S636 = _S635 + _S635;
    float _S637 = _S631.primal_0.x * _S634.differential_0;
    float _S638 = _S637 + _S637;
    float2  _S639 = _S627;
    *&((&_S639)->y) = _S636;
    *&((&_S639)->x) = _S638;
    dpdpx_2->primal_0 = _S631.primal_0;
    dpdpx_2->differential_0 = _S639;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S640, float _S641, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_5)
{
    float2  _S642 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S643 = { _S642, _S642 };
    _s_diff_ctx_5->_S612 = _S643;
    (&_s_diff_ctx_5->_S612)->primal_0 = _S642;
    (&_s_diff_ctx_5->_S612)->differential_0 = _S642;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S644 = *_S640;
    _s_diff_ctx_5->_S612 = *_S640;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S645 = _S644;
    DiffPair_float_0 _S646 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S647 = { _S646 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S648;
    (&_S648)->_S614 = _S643;
    (&_S648)->_S615 = _S646;
    (&_S648)->_S616 = _S647;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S645, _S641, &_S648);
    *_S640 = _S645;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, float4  dpradial_coeffs_1, float2  dptangential_coeffs_1, float2  dpthin_prism_coeffs_1, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_6)
{
    DiffPair_float_0 _S649 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S650 = { _S649, _S649 };
    _s_diff_ctx_6->_S584 = _S649;
    _s_diff_ctx_6->_S585 = _S649;
    _s_diff_ctx_6->_S586 = _S650;
    _s_diff_ctx_6->_S587 = _S649;
    _s_diff_ctx_6->_S588 = _S649;
    _s_diff_ctx_6->_S589 = _S650;
    (&_s_diff_ctx_6->_S584)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S584)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S585)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S585)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S587)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S587)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S588)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S588)->differential_0 = 0.0f;
    float2  _S651 = make_float2 (0.0f);
    CameraDistortion_0 _S652 = s_primal_ctx_CameraDistortion_x24init_0(dpradial_coeffs_1, dptangential_coeffs_1, dpthin_prism_coeffs_1);
    float2  _S653 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S654 = length_0(_S653);
    float _S655 = dpmean3d_0.z;
    float _S656 = s_primal_ctx_atan2_0(_S654, _S655);
    float k_1;
    if(_S656 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S656 * _S656 / 3.0f) / _S655;
    }
    else
    {
        k_1 = _S656 / _S654;
    }
    float2  _S657 = _S653 * make_float2 (k_1);
    float k1_1 = _S652.radial_coeffs_0.x;
    float k2_1 = _S652.radial_coeffs_0.y;
    float k3_1 = _S652.radial_coeffs_0.z;
    float k4_1 = _S652.radial_coeffs_0.w;
    float p1_1 = _S652.tangential_coeffs_0.x;
    float p2_1 = _S652.tangential_coeffs_0.y;
    float sx1_1 = _S652.thin_prism_coeffs_0.x;
    float sy1_1 = _S652.thin_prism_coeffs_0.y;
    float u_3 = _S657.x;
    float v_3 = _S657.y;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float _S658 = 2.0f * p1_1;
    float _S659 = 2.0f * p2_1;
    float2  _S660 = _S657 * make_float2 (1.0f + r2_3 * (k1_1 + r2_3 * (k2_1 + r2_3 * (k3_1 + r2_3 * k4_1)))) + make_float2 (_S658 * u_3 * v_3 + p2_1 * (r2_3 + 2.0f * u_3 * u_3) + sx1_1 * r2_3, _S659 * u_3 * v_3 + p1_1 * (r2_3 + 2.0f * v_3 * v_3) + sy1_1 * r2_3);
    float2  _S661 = make_float2 (dpfx_0 * _S660.x + dpcx_0, dpfy_0 * _S660.y + dpcy_0);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (0.0f);
    float _S662 = s_primal_ctx_s_primal_ctx_atan2_0(_S654, _S655);
    bool _S663 = _S662 < 0.00100000004749745f;
    float _S664;
    float _S665;
    float _S666;
    if(_S663)
    {
        float _S667 = 1.0f - _S662 * _S662 / 3.0f;
        float _S668 = _S655 * _S655;
        k_1 = _S667 / _S655;
        _S664 = 0.0f;
        _S665 = _S668;
        _S666 = _S667;
    }
    else
    {
        float _S669 = _S654 * _S654;
        k_1 = _S662 / _S654;
        _S664 = _S669;
        _S665 = 0.0f;
        _S666 = 0.0f;
    }
    float2  _S670 = make_float2 (k_1);
    float2  _S671 = _S653 * make_float2 (k_1);
    float u_4 = _S671.x;
    float v_4 = _S671.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float _S672 = k3_1 + r2_4 * k4_1;
    float _S673 = k2_1 + r2_4 * _S672;
    float _S674 = k1_1 + r2_4 * _S673;
    float2  _S675 = make_float2 (dpfx_0, 0.0f);
    float2  _S676 = _S671 * _S675;
    float _S677 = p2_1 * dpfx_0;
    float _S678 = _S676.x + _S676.y;
    float _S679 = r2_4 * _S678;
    float _S680 = r2_4 * _S679;
    float _S681 = sx1_1 * dpfx_0 + _S677 + _S674 * _S678 + _S673 * _S679 + _S672 * _S680 + k4_1 * (r2_4 * _S680);
    float _S682 = v_4 * _S681;
    float _S683 = u_4 * _S681;
    float2  _S684 = make_float2 (1.0f + r2_4 * _S674) * _S675 + make_float2 (2.0f * u_4 * _S677 + 2.0f * (u_4 * _S677) + _S658 * (v_4 * dpfx_0) + _S683 + _S683, _S658 * u_4 * dpfx_0 + _S682 + _S682);
    float2  _S685 = _S653 * _S684;
    float2  _S686 = _S670 * _S684;
    float _S687 = _S685.x + _S685.y;
    if(_S663)
    {
        float _S688 = _S687 / _S665;
        float _S689 = _S666 * - _S688;
        float _S690 = _S662 * (0.3333333432674408f * - (_S655 * _S688));
        k_1 = _S690 + _S690;
        _S664 = _S689;
        _S665 = 0.0f;
    }
    else
    {
        float _S691 = _S687 / _S664;
        float _S692 = _S662 * - _S691;
        k_1 = _S654 * _S691;
        _S664 = 0.0f;
        _S665 = _S692;
    }
    DiffPair_float_0 _S693;
    (&_S693)->primal_0 = _S654;
    (&_S693)->differential_0 = 0.0f;
    DiffPair_float_0 _S694;
    (&_S694)->primal_0 = _S655;
    (&_S694)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S693, &_S694, k_1, &_s_diff_ctx_6->_S586);
    _s_diff_ctx_6->_S584 = _S693;
    _s_diff_ctx_6->_S585 = _S694;
    float _S695 = _S694.differential_0 + _S664;
    float _S696 = _S693.differential_0 + _S665;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S697;
    (&_S697)->primal_0 = _S653;
    (&_S697)->differential_0 = _S651;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S698 = { _S651, _S651 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S699;
    (&_S699)->_S612 = _S698;
    s_primal_ctx_s_bwd_length_impl_0(&_S697, _S696, &_S699);
    float2  _S700 = _S697.differential_0 + _S686;
    float3  _S701 = make_float3 (_S700.x, _S700.y, _S695);
    Matrix<float, 2, 3>  _S702 = J_9;
    _S702[int(0)] = _S701;
    if(_S663)
    {
        float _S703 = 1.0f - _S662 * _S662 / 3.0f;
        float _S704 = _S655 * _S655;
        k_1 = _S703 / _S655;
        _S664 = 0.0f;
        _S665 = _S704;
        _S666 = _S703;
    }
    else
    {
        float _S705 = _S654 * _S654;
        k_1 = _S662 / _S654;
        _S664 = _S705;
        _S665 = 0.0f;
        _S666 = 0.0f;
    }
    float2  _S706 = make_float2 (k_1);
    float2  _S707 = _S653 * make_float2 (k_1);
    float u_5 = _S707.x;
    float v_5 = _S707.y;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float _S708 = k3_1 + r2_5 * k4_1;
    float _S709 = k2_1 + r2_5 * _S708;
    float _S710 = k1_1 + r2_5 * _S709;
    float2  _S711 = make_float2 (0.0f, dpfy_0);
    float2  _S712 = _S707 * _S711;
    float _S713 = p1_1 * dpfy_0;
    float _S714 = _S712.x + _S712.y;
    float _S715 = r2_5 * _S714;
    float _S716 = r2_5 * _S715;
    float _S717 = sy1_1 * dpfy_0 + _S713 + _S710 * _S714 + _S709 * _S715 + _S708 * _S716 + k4_1 * (r2_5 * _S716);
    float _S718 = v_5 * _S717;
    float _S719 = u_5 * _S717;
    float2  _S720 = make_float2 (1.0f + r2_5 * _S710) * _S711 + make_float2 (_S659 * (v_5 * dpfy_0) + _S719 + _S719, 2.0f * v_5 * _S713 + 2.0f * (v_5 * _S713) + _S659 * u_5 * dpfy_0 + _S718 + _S718);
    float2  _S721 = _S653 * _S720;
    float2  _S722 = _S706 * _S720;
    float _S723 = _S721.x + _S721.y;
    if(_S663)
    {
        float _S724 = _S723 / _S665;
        float _S725 = _S666 * - _S724;
        float _S726 = _S662 * (0.3333333432674408f * - (_S655 * _S724));
        k_1 = _S726 + _S726;
        _S664 = _S725;
        _S665 = 0.0f;
    }
    else
    {
        float _S727 = _S723 / _S664;
        float _S728 = _S662 * - _S727;
        k_1 = _S654 * _S727;
        _S664 = 0.0f;
        _S665 = _S728;
    }
    DiffPair_float_0 _S729;
    (&_S729)->primal_0 = _S654;
    (&_S729)->differential_0 = 0.0f;
    DiffPair_float_0 _S730;
    (&_S730)->primal_0 = _S655;
    (&_S730)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S729, &_S730, k_1, &_s_diff_ctx_6->_S589);
    _s_diff_ctx_6->_S587 = _S729;
    _s_diff_ctx_6->_S588 = _S730;
    float _S731 = _S730.differential_0 + _S664;
    float _S732 = _S729.differential_0 + _S665;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S733;
    (&_S733)->primal_0 = _S653;
    (&_S733)->differential_0 = _S651;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S734;
    (&_S734)->_S612 = _S698;
    s_primal_ctx_s_bwd_length_impl_0(&_S733, _S732, &_S734);
    float2  _S735 = _S733.differential_0 + _S722;
    float3  _S736 = make_float3 (_S735.x, _S735.y, _S731);
    _S702[int(1)] = _S736;
    *dpcov2d_0 = s_primal_ctx_mul_3(s_primal_ctx_mul_2(_S702, dpcov3d_0), transpose_1(_S702));
    *dpmean2d_0 = _S661;
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
    DiffPair_1 _S737 = *dpdpx_3;
    float _S738 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_7->_S617)->primal_0);
    float _S739 = s_primal_ctx_sqrt_0(_S738);
    float _S740 = 0.5f / _S739 * (*dpdpx_3).differential_0.differential_0;
    float _S741 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S739 * _S739));
    DiffPair_float_0 _S742;
    (&_S742)->primal_0 = _S738;
    (&_S742)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S742, _S741);
    DiffPair_float_0 _S743;
    (&_S743)->primal_0 = 1.00000001168609742e-07f;
    (&_S743)->differential_0 = 0.0f;
    DiffPair_float_0 _S744;
    (&_S744)->primal_0 = (&_s_diff_ctx_7->_S617)->primal_0;
    (&_S744)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S743, &_S744, _S742.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S744.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S740;
    dpdpx_3->primal_0 = _S737.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S745, DiffPair_float_0 * _S746, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_1 _S747 = *_S745;
    DiffPair_float_0 _S748 = _s_diff_ctx_8->_S613;
    DiffPair_float_0 _S749 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S750;
    (&_S750)->_S617 = _S749;
    s_primal_ctx_d_sqrt_0(&_S748, (*_S746).primal_0, &_S750);
    DiffPair_float_0 _S751 = { (*_S745).differential_0.primal_0, (*_S745).differential_0.differential_0 };
    DiffPair_1 _S752;
    (&_S752)->primal_0 = _s_diff_ctx_8->_S613;
    (&_S752)->differential_0 = _S751;
    DiffPair_float_0 _S753;
    (&_S753)->primal_0 = (*_S746).primal_0;
    (&_S753)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S754 = _S750;
    s_bwd_prop_d_sqrt_0(&_S752, &_S753, &_S754);
    DiffPair_float_0 _S755 = { _S752.differential_0.primal_0, _S752.differential_0.differential_0 };
    _S746->primal_0 = (*_S746).primal_0;
    _S746->differential_0 = _S753.differential_0;
    _S745->primal_0 = _S747.primal_0;
    _S745->differential_0 = _S755;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S756, float _s_dOut_2)
{
    DiffPair_float_0 _S757;
    (&_S757)->primal_0 = (*_S756).primal_0;
    (&_S757)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S757, _s_dOut_2);
    _S756->primal_0 = (*_S756).primal_0;
    _S756->differential_0 = _S757.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_9)
{
    DiffPair_0 _S758 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_9->_S614)->primal_0)->x) * *&((&(&_s_diff_ctx_9->_S614)->primal_0)->x) + *&((&(&_s_diff_ctx_9->_S614)->primal_0)->y) * *&((&(&_s_diff_ctx_9->_S614)->primal_0)->y);
    DiffPair_float_0 _S759 = { len_0, 0.0f };
    float2  _S760 = make_float2 (0.0f);
    float _S761 = (*dpdpx_5).differential_0.differential_0.x;
    float _S762 = _S761 + _S761;
    float _S763 = (&_s_diff_ctx_9->_S615)->differential_0 * _S762;
    float _S764 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S765 = (&_s_diff_ctx_9->_S615)->differential_0 * _S764;
    DiffPair_float_0 _S766 = { 0.0f, *&((&(&_s_diff_ctx_9->_S614)->primal_0)->x) * _S762 + *&((&(&_s_diff_ctx_9->_S614)->primal_0)->y) * _S764 };
    DiffPair_1 _S767;
    (&_S767)->primal_0 = _S759;
    (&_S767)->differential_0 = _S766;
    DiffPair_float_0 _S768;
    (&_S768)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S768)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S767, &_S768, &_s_diff_ctx_9->_S616);
    DiffPair_float_0 _S769;
    (&_S769)->primal_0 = len_0;
    (&_S769)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S769, 0.0f);
    float _S770 = _S767.differential_0.primal_0 + _S769.differential_0;
    float _S771 = *&((&(&_s_diff_ctx_9->_S614)->primal_0)->y) * _S770;
    float _S772 = _S765 + _S771 + _S771;
    float _S773 = *&((&(&_s_diff_ctx_9->_S614)->primal_0)->x) * _S770;
    float _S774 = _S763 + _S773 + _S773;
    float2  _S775 = _S760;
    *&((&_S775)->y) = _S772;
    *&((&_S775)->x) = _S774;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S758.differential_0.primal_0 + _S775, _S760 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S768.differential_0;
    dpdpx_5->primal_0 = _S758.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_3)
{
    float _S776 = (*dpdpx_7).primal_0.x;
    float _S777 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S778;
    (&_S778)->primal_0 = _S776 * _S776 + _S777 * _S777;
    (&_S778)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S778, _s_dOut_3);
    float _S779 = (*dpdpx_7).primal_0.y * _S778.differential_0;
    float _S780 = _S779 + _S779;
    float _S781 = (*dpdpx_7).primal_0.x * _S778.differential_0;
    float _S782 = _S781 + _S781;
    float2  _S783 = make_float2 (0.0f);
    *&((&_S783)->y) = _S780;
    *&((&_S783)->x) = _S782;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S783;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S784, DiffPair_float_0 * _S785, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_0 _S786 = *_S784;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S787 = _s_diff_ctx_10->_S612;
    float2  _S788 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S789 = { _S788, _S788 };
    DiffPair_float_0 _S790 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S791 = { _S790 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S792;
    (&_S792)->_S614 = _S789;
    (&_S792)->_S615 = _S790;
    (&_S792)->_S616 = _S791;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S787, (*_S785).primal_0, &_S792);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S793 = { (*_S784).differential_0.primal_0, (*_S784).differential_0.differential_0 };
    DiffPair_0 _S794;
    (&_S794)->primal_0 = _s_diff_ctx_10->_S612;
    (&_S794)->differential_0 = _S793;
    DiffPair_float_0 _S795;
    (&_S795)->primal_0 = (*_S785).primal_0;
    (&_S795)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S796 = _S792;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S794, &_S795, &_S796);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S797;
    (&_S797)->primal_0 = (&_s_diff_ctx_10->_S612)->primal_0;
    (&_S797)->differential_0 = _S788;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S797, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S798 = { _S794.differential_0.primal_0 + _S797.differential_0, _S794.differential_0.differential_0 };
    _S785->primal_0 = (*_S785).primal_0;
    _S785->differential_0 = _S795.differential_0;
    _S784->primal_0 = _S786.primal_0;
    _S784->differential_0 = _S798;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_1 _S799 = *dpdpy_1;
    DiffPair_1 _S800 = *dpdpx_8;
    float _S801 = - (&_s_diff_ctx_11->_S595)->primal_0;
    float _S802 = (&_s_diff_ctx_11->_S596)->primal_0 * (&_s_diff_ctx_11->_S596)->primal_0 + (&_s_diff_ctx_11->_S595)->primal_0 * (&_s_diff_ctx_11->_S595)->primal_0;
    float _S803 = _S802 * _S802;
    float _S804 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S803;
    float _S805 = (&_s_diff_ctx_11->_S596)->primal_0 * - _S804;
    float _S806 = (&_s_diff_ctx_11->_S595)->primal_0 * _S805;
    float _S807 = (&_s_diff_ctx_11->_S596)->primal_0 * _S805;
    float _S808 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S803;
    float _S809 = _S801 * - _S808;
    float _S810 = (&_s_diff_ctx_11->_S595)->primal_0 * _S809;
    float _S811 = (&_s_diff_ctx_11->_S596)->primal_0 * _S809;
    DiffPair_float_0 dpdpx_9 = { _S811 + _S811 + ((*dpdpx_8).differential_0.primal_0 + (_S807 + _S807 + _S802 * _S804)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S806 + _S806 + (*dpdpy_1).differential_0.primal_0 + _S810 + _S810 + - (_S802 * _S808), 0.0f };
    float _S812 = (&_s_diff_ctx_11->_S596)->primal_0 / _S802 * (*dpdpy_1).differential_0.differential_0 + _S801 / _S802 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S812;
    dpdpy_1->primal_0 = _S799.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S800.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S813, DiffPair_1 * _S814, DiffPair_float_0 * _S815, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_1 _S816 = *_S813;
    DiffPair_1 _S817 = *_S814;
    DiffPair_float_0 _S818 = _s_diff_ctx_12->_S582;
    DiffPair_float_0 _S819 = _s_diff_ctx_12->_S583;
    DiffPair_float_0 _S820 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S821;
    (&_S821)->_S595 = _S820;
    (&_S821)->_S596 = _S820;
    s_primal_ctx_d_atan2_0(&_S818, &_S819, (*_S815).primal_0, &_S821);
    DiffPair_float_0 _S822 = { (*_S814).differential_0.primal_0, (*_S814).differential_0.differential_0 };
    DiffPair_float_0 _S823 = { (*_S813).differential_0.primal_0, (*_S813).differential_0.differential_0 };
    DiffPair_1 _S824;
    (&_S824)->primal_0 = _s_diff_ctx_12->_S582;
    (&_S824)->differential_0 = _S823;
    DiffPair_1 _S825;
    (&_S825)->primal_0 = _s_diff_ctx_12->_S583;
    (&_S825)->differential_0 = _S822;
    DiffPair_float_0 _S826;
    (&_S826)->primal_0 = (*_S815).primal_0;
    (&_S826)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S827 = _S821;
    s_bwd_prop_d_atan2_0(&_S824, &_S825, &_S826, &_S827);
    DiffPair_float_0 _S828 = { _S825.differential_0.primal_0, _S825.differential_0.differential_0 };
    DiffPair_float_0 _S829 = { _S824.differential_0.primal_0, _S824.differential_0.differential_0 };
    _S815->primal_0 = (*_S815).primal_0;
    _S815->differential_0 = _S826.differential_0;
    _S813->primal_0 = _S816.primal_0;
    _S813->differential_0 = _S829;
    _S814->primal_0 = _S817.primal_0;
    _S814->differential_0 = _S828;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S830, DiffPair_float_0 * _S831, float _s_dOut_4)
{
    DiffPair_float_0 _S832;
    (&_S832)->primal_0 = (*_S830).primal_0;
    (&_S832)->differential_0 = 0.0f;
    DiffPair_float_0 _S833;
    (&_S833)->primal_0 = (*_S831).primal_0;
    (&_S833)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S832, &_S833, _s_dOut_4);
    _S831->primal_0 = (*_S831).primal_0;
    _S831->differential_0 = _S833.differential_0;
    _S830->primal_0 = (*_S830).primal_0;
    _S830->differential_0 = _S832.differential_0;
    return;
}

inline __device__ void s_bwd_prop_CameraDistortion_x24init_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_2, CameraDistortion_0 * _s_dOut_5)
{
    float2  _S834 = _s_dOut_5->thin_prism_coeffs_0;
    dpthin_prism_coeffs_2->primal_0 = (*dpthin_prism_coeffs_2).primal_0;
    dpthin_prism_coeffs_2->differential_0 = _S834;
    float2  _S835 = _s_dOut_5->tangential_coeffs_0;
    dptangential_coeffs_2->primal_0 = (*dptangential_coeffs_2).primal_0;
    dptangential_coeffs_2->differential_0 = _S835;
    float4  _S836 = _s_dOut_5->radial_coeffs_0;
    dpradial_coeffs_2->primal_0 = (*dpradial_coeffs_2).primal_0;
    dpradial_coeffs_2->differential_0 = _S836;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_3, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S837 = *dpcov3d_1;
    DiffPair_float_0 _S838 = *dpfx_1;
    DiffPair_float_0 _S839 = *dpfy_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S840 = *dpradial_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S841 = *dptangential_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S842 = *dpthin_prism_coeffs_3;
    float2  _S843 = make_float2 (0.0f);
    CameraDistortion_0 _S844 = s_primal_ctx_CameraDistortion_x24init_0((*dpradial_coeffs_3).primal_0, (*dptangential_coeffs_3).primal_0, (*dpthin_prism_coeffs_3).primal_0);
    float2  _S845 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S846 = length_0(_S845);
    float _S847 = (*dpmean3d_1).primal_0.z;
    float _S848 = s_primal_ctx_atan2_0(_S846, _S847);
    bool _S849 = _S848 < 0.00100000004749745f;
    float k_2;
    float _S850;
    float _S851;
    float _S852;
    if(_S849)
    {
        float _S853 = 1.0f - _S848 * _S848 / 3.0f;
        float _S854 = _S847 * _S847;
        k_2 = _S853 / _S847;
        _S850 = 0.0f;
        _S851 = _S854;
        _S852 = _S853;
    }
    else
    {
        float _S855 = _S846 * _S846;
        k_2 = _S848 / _S846;
        _S850 = _S855;
        _S851 = 0.0f;
        _S852 = 0.0f;
    }
    float2  _S856 = make_float2 (k_2);
    float2  _S857 = _S845 * make_float2 (k_2);
    float k1_2 = _S844.radial_coeffs_0.x;
    float k2_2 = _S844.radial_coeffs_0.y;
    float k3_2 = _S844.radial_coeffs_0.z;
    float k4_2 = _S844.radial_coeffs_0.w;
    float p1_2 = _S844.tangential_coeffs_0.x;
    float p2_2 = _S844.tangential_coeffs_0.y;
    float sx1_2 = _S844.thin_prism_coeffs_0.x;
    float sy1_2 = _S844.thin_prism_coeffs_0.y;
    float u_6 = _S857.x;
    float v_6 = _S857.y;
    float r2_6 = u_6 * u_6 + v_6 * v_6;
    float _S858 = k3_2 + r2_6 * k4_2;
    float _S859 = k2_2 + r2_6 * _S858;
    float _S860 = k1_2 + r2_6 * _S859;
    float radial_0 = 1.0f + r2_6 * _S860;
    float2  _S861 = make_float2 (radial_0);
    float _S862 = 2.0f * p1_2;
    float _S863 = _S862 * u_6;
    float _S864 = 2.0f * u_6;
    float _S865 = r2_6 + _S864 * u_6;
    float _S866 = 2.0f * p2_2;
    float _S867 = _S866 * u_6;
    float _S868 = 2.0f * v_6;
    float _S869 = r2_6 + _S868 * v_6;
    float2  _S870 = _S857 * make_float2 (radial_0) + make_float2 (_S863 * v_6 + p2_2 * _S865 + sx1_2 * r2_6, _S867 * v_6 + p1_2 * _S869 + sy1_2 * r2_6);
    float _S871 = _S870.x;
    float _S872 = _S870.y;
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float _S873 = s_primal_ctx_s_primal_ctx_atan2_0(_S846, _S847);
    bool _S874 = _S873 < 0.00100000004749745f;
    float _S875;
    float _S876;
    float _S877;
    if(_S874)
    {
        float _S878 = 1.0f - _S873 * _S873 / 3.0f;
        float _S879 = _S847 * _S847;
        k_2 = _S878 / _S847;
        _S875 = 0.0f;
        _S876 = _S879;
        _S877 = _S878;
    }
    else
    {
        float _S880 = _S846 * _S846;
        k_2 = _S873 / _S846;
        _S875 = _S880;
        _S876 = 0.0f;
        _S877 = 0.0f;
    }
    float2  _S881 = make_float2 (k_2);
    float2  _S882 = _S845 * make_float2 (k_2);
    float u_7 = _S882.x;
    float v_7 = _S882.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S883 = k3_2 + r2_7 * k4_2;
    float _S884 = k2_2 + r2_7 * _S883;
    float _S885 = k1_2 + r2_7 * _S884;
    float2  _S886 = make_float2 (1.0f + r2_7 * _S885);
    float _S887 = _S862 * u_7;
    float _S888 = 2.0f * u_7;
    float2  _S889 = make_float2 (_S838.primal_0, 0.0f);
    float2  _S890 = _S882 * _S889;
    float _S891 = p2_2 * _S838.primal_0;
    float _S892 = v_7 * _S838.primal_0;
    float _S893 = _S890.x + _S890.y;
    float _S894 = r2_7 * _S893;
    float _S895 = r2_7 * _S894;
    float _S896 = r2_7 * _S895;
    float _S897 = sx1_2 * _S838.primal_0 + _S891 + _S885 * _S893 + _S884 * _S894 + _S883 * _S895 + k4_2 * _S896;
    float _S898 = v_7 * _S897;
    float _S899 = u_7 * _S897;
    float2  _S900 = _S886 * _S889 + make_float2 (_S888 * _S891 + 2.0f * (u_7 * _S891) + _S862 * _S892 + _S899 + _S899, _S887 * _S838.primal_0 + _S898 + _S898);
    float2  _S901 = _S845 * _S900;
    float2  _S902 = _S881 * _S900;
    float _S903 = _S901.x + _S901.y;
    float k_3;
    float _S904;
    float _S905;
    float _S906;
    float _S907;
    float _S908;
    float _S909;
    float _S910;
    float _S911;
    if(_S874)
    {
        float _S912 = _S903 / _S876;
        float _S913 = _S876 * _S876;
        float _S914 = - _S912;
        float _S915 = _S877 * _S914;
        float _S916 = 0.3333333432674408f * - (_S847 * _S912);
        float _S917 = _S873 * _S916;
        k_2 = _S917 + _S917;
        k_3 = _S915;
        _S904 = 0.0f;
        _S905 = 0.0f;
        _S906 = 0.0f;
        _S907 = 0.0f;
        _S908 = _S916;
        _S909 = _S912;
        _S910 = _S914;
        _S911 = _S913;
    }
    else
    {
        float _S918 = _S903 / _S875;
        float _S919 = _S875 * _S875;
        float _S920 = - _S918;
        float _S921 = _S873 * _S920;
        k_2 = _S846 * _S918;
        k_3 = 0.0f;
        _S904 = _S921;
        _S905 = _S918;
        _S906 = _S920;
        _S907 = _S919;
        _S908 = 0.0f;
        _S909 = 0.0f;
        _S910 = 0.0f;
        _S911 = 0.0f;
    }
    DiffPair_float_0 _S922 = { _S846, 0.0f };
    DiffPair_float_0 _S923 = { _S847, 0.0f };
    float _S924 = (&_s_diff_ctx_13->_S585)->differential_0 + k_3;
    float _S925 = (&_s_diff_ctx_13->_S584)->differential_0 + _S904;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S926 = { _S845, _S843 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S927;
    (&_S927)->primal_0 = _S845;
    (&_S927)->differential_0 = _S843;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S928 = { _S843, _S843 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S929;
    (&_S929)->_S612 = _S928;
    s_primal_ctx_s_bwd_length_impl_0(&_S927, _S925, &_S929);
    float2  _S930 = _S927.differential_0 + _S902;
    float3  _S931 = make_float3 (_S930.x, _S930.y, _S924);
    Matrix<float, 2, 3>  _S932 = J_10;
    _S932[int(0)] = _S931;
    float _S933;
    float _S934;
    if(_S874)
    {
        float _S935 = 1.0f - _S873 * _S873 / 3.0f;
        float _S936 = _S847 * _S847;
        k_3 = _S935 / _S847;
        _S904 = 0.0f;
        _S933 = _S936;
        _S934 = _S935;
    }
    else
    {
        float _S937 = _S846 * _S846;
        k_3 = _S873 / _S846;
        _S904 = _S937;
        _S933 = 0.0f;
        _S934 = 0.0f;
    }
    float2  _S938 = make_float2 (k_3);
    float2  _S939 = _S845 * make_float2 (k_3);
    float u_8 = _S939.x;
    float v_8 = _S939.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S940 = k3_2 + r2_8 * k4_2;
    float _S941 = k2_2 + r2_8 * _S940;
    float _S942 = k1_2 + r2_8 * _S941;
    float2  _S943 = make_float2 (1.0f + r2_8 * _S942);
    float _S944 = _S866 * u_8;
    float _S945 = 2.0f * v_8;
    float2  _S946 = make_float2 (0.0f, _S839.primal_0);
    float2  _S947 = _S939 * _S946;
    float _S948 = p1_2 * _S839.primal_0;
    float _S949 = v_8 * _S839.primal_0;
    float _S950 = _S947.x + _S947.y;
    float _S951 = r2_8 * _S950;
    float _S952 = r2_8 * _S951;
    float _S953 = r2_8 * _S952;
    float _S954 = sy1_2 * _S839.primal_0 + _S948 + _S942 * _S950 + _S941 * _S951 + _S940 * _S952 + k4_2 * _S953;
    float _S955 = v_8 * _S954;
    float _S956 = u_8 * _S954;
    float2  _S957 = _S943 * _S946 + make_float2 (_S866 * _S949 + _S956 + _S956, _S945 * _S948 + 2.0f * (v_8 * _S948) + _S944 * _S839.primal_0 + _S955 + _S955);
    float2  _S958 = _S845 * _S957;
    float2  _S959 = _S938 * _S957;
    float _S960 = _S958.x + _S958.y;
    float _S961;
    float _S962;
    float _S963;
    float _S964;
    float _S965;
    float _S966;
    float _S967;
    float _S968;
    float _S969;
    if(_S874)
    {
        float _S970 = _S960 / _S933;
        float _S971 = _S933 * _S933;
        float _S972 = - _S970;
        float _S973 = _S934 * _S972;
        float _S974 = 0.3333333432674408f * - (_S847 * _S970);
        float _S975 = _S873 * _S974;
        k_3 = _S975 + _S975;
        _S961 = _S973;
        _S962 = 0.0f;
        _S963 = 0.0f;
        _S964 = 0.0f;
        _S965 = 0.0f;
        _S966 = _S974;
        _S967 = _S970;
        _S968 = _S972;
        _S969 = _S971;
    }
    else
    {
        float _S976 = _S960 / _S904;
        float _S977 = _S904 * _S904;
        float _S978 = - _S976;
        float _S979 = _S873 * _S978;
        k_3 = _S846 * _S976;
        _S961 = 0.0f;
        _S962 = _S979;
        _S963 = _S976;
        _S964 = _S978;
        _S965 = _S977;
        _S966 = 0.0f;
        _S967 = 0.0f;
        _S968 = 0.0f;
        _S969 = 0.0f;
    }
    float _S980 = (&_s_diff_ctx_13->_S588)->differential_0 + _S961;
    float _S981 = (&_s_diff_ctx_13->_S587)->differential_0 + _S962;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S982;
    (&_S982)->primal_0 = _S845;
    (&_S982)->differential_0 = _S843;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S983;
    (&_S983)->_S612 = _S928;
    s_primal_ctx_s_bwd_length_impl_0(&_S982, _S981, &_S983);
    float2  _S984 = _S982.differential_0 + _S959;
    float3  _S985 = make_float3 (_S984.x, _S984.y, _S980);
    _S932[int(1)] = _S985;
    Matrix<float, 3, 2>  _S986 = transpose_1(_S932);
    CameraDistortion_0 _S987 = CameraDistortion_x24_syn_dzero_0();
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S988;
    (&_S988)->primal_0 = s_primal_ctx_mul_2(_S932, _S837.primal_0);
    (&_S988)->differential_0 = J_10;
    Matrix<float, 3, 2>  _S989 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S990;
    (&_S990)->primal_0 = _S986;
    (&_S990)->differential_0 = _S989;
    s_bwd_prop_mul_1(&_S988, &_S990, dpcov2d_1);
    Matrix<float, 2, 3>  _S991 = transpose_2(_S990.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S992;
    (&_S992)->primal_0 = _S932;
    (&_S992)->differential_0 = J_10;
    Matrix<float, 3, 3>  _S993 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S994;
    (&_S994)->primal_0 = _S837.primal_0;
    (&_S994)->differential_0 = _S993;
    s_bwd_prop_mul_2(&_S992, &_S994, _S988.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S995 = _S994;
    Matrix<float, 2, 3>  _S996 = _S991 + _S992.differential_0;
    float2  _S997 = _S843;
    *&((&_S997)->y) = _S996.rows[int(1)].y;
    *&((&_S997)->x) = _S996.rows[int(1)].x;
    float2  _S998 = _S997;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S999 = { _S843, _S997 };
    DiffPair_0 _S1000;
    (&_S1000)->primal_0 = _S926;
    (&_S1000)->differential_0 = _S999;
    DiffPair_float_0 _S1001;
    (&_S1001)->primal_0 = _S981;
    (&_S1001)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1002 = _S983;
    s_bwd_prop_s_bwd_length_impl_0(&_S1000, &_S1001, &_S1002);
    DiffPair_0 _S1003 = _S1000;
    DiffPair_float_0 _S1004 = _S1001;
    DiffPair_float_0 _S1005 = { 0.0f, _S996.rows[int(1)].z };
    DiffPair_float_0 _S1006 = { 0.0f, _S1001.differential_0 };
    DiffPair_1 _S1007;
    (&_S1007)->primal_0 = _S922;
    (&_S1007)->differential_0 = _S1006;
    DiffPair_1 _S1008;
    (&_S1008)->primal_0 = _S923;
    (&_S1008)->differential_0 = _S1005;
    DiffPair_float_0 _S1009;
    (&_S1009)->primal_0 = k_3;
    (&_S1009)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1007, &_S1008, &_S1009, &_s_diff_ctx_13->_S589);
    DiffPair_1 _S1010 = _S1007;
    DiffPair_1 _S1011 = _S1008;
    DiffPair_float_0 _S1012 = _S1009;
    if(_S874)
    {
        float _S1013 = _S1012.differential_0 + _S1012.differential_0;
        float _S1014 = _S966 * _S1013;
        float _S1015 = - (0.3333333432674408f * (_S873 * _S1013));
        float _S1016 = _S968 * _S996.rows[int(1)].z;
        float _S1017 = (_S847 * _S1015 + - (_S934 * _S996.rows[int(1)].z)) / _S969;
        float _S1018 = _S960 * - _S1017;
        float _S1019 = _S967 * _S1015 + _S1011.differential_0.primal_0;
        k_3 = _S933 * _S1017;
        _S961 = 0.0f;
        _S962 = _S1018;
        _S963 = _S1016;
        _S964 = _S1010.differential_0.primal_0;
        _S965 = _S1014;
        _S966 = _S1019;
    }
    else
    {
        float _S1020 = _S964 * _S1004.differential_0;
        float _S1021 = (_S846 * _S1012.differential_0 + - (_S873 * _S1004.differential_0)) / _S965;
        float _S1022 = _S960 * - _S1021;
        float _S1023 = _S963 * _S1012.differential_0 + _S1010.differential_0.primal_0;
        k_3 = _S904 * _S1021;
        _S961 = _S1022;
        _S962 = 0.0f;
        _S963 = 0.0f;
        _S964 = _S1023;
        _S965 = _S1020;
        _S966 = _S1011.differential_0.primal_0;
    }
    float2  _S1024 = _S938 * _S998;
    float2  _S1025 = _S957 * _S998;
    float2  _S1026 = _S843;
    *&((&_S1026)->y) = k_3;
    *&((&_S1026)->x) = k_3;
    float2  _S1027 = _S957 * _S1026;
    float2  _S1028 = _S1024 + _S845 * _S1026;
    float _S1029 = _S1028.x;
    float _S1030 = _S1029 + _S1029;
    float _S1031 = _S954 * _S1030;
    float _S1032 = _S1028.y + _S1028.y;
    float _S1033 = _S954 * _S1032;
    float _S1034 = u_8 * _S1030 + v_8 * _S1032;
    float _S1035 = k4_2 * _S1034;
    float _S1036 = _S953 * _S1034;
    float _S1037 = _S952 * _S1034;
    float _S1038 = _S952 * _S1035;
    float _S1039 = _S951 * _S1034;
    float _S1040 = _S940 * _S1034 + r2_8 * _S1035;
    float _S1041 = _S951 * _S1040;
    float _S1042 = _S950 * _S1034;
    float _S1043 = _S941 * _S1034 + r2_8 * _S1040;
    float _S1044 = _S950 * _S1043;
    float _S1045 = _S942 * _S1034 + r2_8 * _S1043;
    float _S1046 = _S866 * _S1028.x;
    float _S1047 = _S949 * _S1028.x;
    float _S1048 = v_8 * _S1046;
    float _S1049 = _S839.primal_0 * _S1046;
    float _S1050 = _S944 * _S1028.y;
    float _S1051 = _S839.primal_0 * _S1028.y;
    float _S1052 = 2.0f * _S1028.y;
    float _S1053 = _S948 * _S1052;
    float _S1054 = _S948 * _S1028.y;
    float _S1055 = _S1034 + v_8 * _S1052 + _S945 * _S1028.y;
    float _S1056 = p1_2 * _S1055;
    float _S1057 = _S839.primal_0 * _S1055;
    float _S1058 = sy1_2 * _S1034;
    float _S1059 = _S839.primal_0 * _S1034;
    float2  _S1060 = _S943 * _S1028;
    float2  _S1061 = _S946 * _S1028;
    float2  _S1062 = _S843;
    *&((&_S1062)->y) = _S1045;
    *&((&_S1062)->x) = _S1045;
    float _S1063 = _S1061.x + _S1061.y;
    float _S1064 = _S1042 + r2_8 * _S1063;
    float _S1065 = _S1039 + r2_8 * _S1064;
    float _S1066 = _S1037 + r2_8 * _S1065;
    float _S1067 = _S1038 + _S1041 + _S1044 + _S942 * _S1063 + _S941 * _S1064 + _S940 * _S1065 + k4_2 * _S1066;
    float _S1068 = v_8 * _S1067;
    float _S1069 = u_8 * _S1067;
    float2  _S1070 = _S1027 + _S1003.differential_0.primal_0;
    float2  _S1071 = _S946 * _S1062 + make_float2 (_S1031 + _S866 * _S1051 + _S1069 + _S1069, _S1033 + _S1049 + _S1053 + 2.0f * _S1054 + _S1068 + _S1068);
    float _S1072 = _S1048 + _S1050 + _S1056 + _S1058 + (_S1060 + _S939 * _S1062).y;
    float _S1073 = _S1047 + u_8 * _S1051;
    float _S1074 = _S1036 + r2_8 * _S1066;
    float2  _S1075 = _S845 * _S1071;
    float _S1076 = _S1025.x + _S1025.y + _S1075.x + _S1075.y;
    float2  _S1077 = _S938 * _S1071 + _S1070;
    if(_S874)
    {
        float _S1078 = _S847 * _S962;
        float _S1079 = _S1076 / _S933;
        float _S1080 = _S873 * (0.3333333432674408f * - (_S963 + _S847 * _S1079));
        float _S1081 = _S1078 + _S1078 + _S934 * - _S1079 + _S966;
        k_3 = _S1080 + _S1080 + _S965;
        _S904 = _S1081;
        _S933 = _S964;
    }
    else
    {
        float _S1082 = _S846 * _S961;
        float _S1083 = _S1076 / _S904;
        float _S1084 = _S1082 + _S1082 + _S873 * - _S1083 + _S964;
        k_3 = _S846 * _S1083 + _S965;
        _S904 = _S966;
        _S933 = _S1084;
    }
    DiffPair_float_0 _S1085;
    (&_S1085)->primal_0 = _S846;
    (&_S1085)->differential_0 = 0.0f;
    DiffPair_float_0 _S1086;
    (&_S1086)->primal_0 = _S847;
    (&_S1086)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1085, &_S1086, k_3);
    float _S1087 = _S1086.differential_0 + _S904;
    float _S1088 = _S1085.differential_0 + _S933;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1089;
    (&_S1089)->primal_0 = _S845;
    (&_S1089)->differential_0 = _S843;
    s_bwd_length_impl_0(&_S1089, _S1088);
    float2  _S1090 = _S1089.differential_0 + _S1077;
    float3  _S1091 = make_float3 (_S1090.x, _S1090.y, _S1087);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1092;
    (&_S1092)->primal_0 = _S845;
    (&_S1092)->differential_0 = _S843;
    s_bwd_length_impl_0(&_S1092, 0.0f);
    float3  _S1093 = _S1091 + make_float3 (_S1092.differential_0.x, _S1092.differential_0.y, 0.0f);
    float2  _S1094 = _S843;
    *&((&_S1094)->y) = _S996.rows[int(0)].y;
    *&((&_S1094)->x) = _S996.rows[int(0)].x;
    float2  _S1095 = _S1094;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1096 = { _S843, _S1094 };
    DiffPair_0 _S1097;
    (&_S1097)->primal_0 = _S926;
    (&_S1097)->differential_0 = _S1096;
    DiffPair_float_0 _S1098;
    (&_S1098)->primal_0 = _S925;
    (&_S1098)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1099 = _S929;
    s_bwd_prop_s_bwd_length_impl_0(&_S1097, &_S1098, &_S1099);
    DiffPair_0 _S1100 = _S1097;
    DiffPair_float_0 _S1101 = _S1098;
    DiffPair_float_0 _S1102 = { 0.0f, _S996.rows[int(0)].z };
    DiffPair_float_0 _S1103 = { 0.0f, _S1098.differential_0 };
    DiffPair_1 _S1104;
    (&_S1104)->primal_0 = _S922;
    (&_S1104)->differential_0 = _S1103;
    DiffPair_1 _S1105;
    (&_S1105)->primal_0 = _S923;
    (&_S1105)->differential_0 = _S1102;
    DiffPair_float_0 _S1106;
    (&_S1106)->primal_0 = k_2;
    (&_S1106)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1104, &_S1105, &_S1106, &_s_diff_ctx_13->_S586);
    DiffPair_1 _S1107 = _S1104;
    DiffPair_1 _S1108 = _S1105;
    DiffPair_float_0 _S1109 = _S1106;
    if(_S874)
    {
        float _S1110 = _S1109.differential_0 + _S1109.differential_0;
        float _S1111 = _S908 * _S1110;
        float _S1112 = - (0.3333333432674408f * (_S873 * _S1110));
        float _S1113 = _S910 * _S996.rows[int(0)].z;
        float _S1114 = (_S847 * _S1112 + - (_S877 * _S996.rows[int(0)].z)) / _S911;
        float _S1115 = _S903 * - _S1114;
        float _S1116 = _S909 * _S1112 + _S1108.differential_0.primal_0;
        k_2 = _S876 * _S1114;
        k_3 = 0.0f;
        _S904 = _S1115;
        _S905 = _S1113;
        _S906 = _S1107.differential_0.primal_0;
        _S907 = _S1111;
        _S908 = _S1116;
    }
    else
    {
        float _S1117 = _S906 * _S1101.differential_0;
        float _S1118 = (_S846 * _S1109.differential_0 + - (_S873 * _S1101.differential_0)) / _S907;
        float _S1119 = _S903 * - _S1118;
        float _S1120 = _S905 * _S1109.differential_0 + _S1107.differential_0.primal_0;
        k_2 = _S875 * _S1118;
        k_3 = _S1119;
        _S904 = 0.0f;
        _S905 = 0.0f;
        _S906 = _S1120;
        _S907 = _S1117;
        _S908 = _S1108.differential_0.primal_0;
    }
    float2  _S1121 = _S881 * _S1095;
    float2  _S1122 = _S900 * _S1095;
    float2  _S1123 = _S843;
    *&((&_S1123)->y) = k_2;
    *&((&_S1123)->x) = k_2;
    float2  _S1124 = _S900 * _S1123;
    float2  _S1125 = _S1121 + _S845 * _S1123;
    float _S1126 = _S1125.x;
    float _S1127 = _S1126 + _S1126;
    float _S1128 = _S897 * _S1127;
    float _S1129 = _S1125.y + _S1125.y;
    float _S1130 = _S897 * _S1129;
    float _S1131 = u_7 * _S1127 + v_7 * _S1129;
    float _S1132 = k4_2 * _S1131;
    float _S1133 = _S896 * _S1131;
    float _S1134 = _S895 * _S1131;
    float _S1135 = _S895 * _S1132;
    float _S1136 = _S894 * _S1131;
    float _S1137 = _S883 * _S1131 + r2_7 * _S1132;
    float _S1138 = _S894 * _S1137;
    float _S1139 = _S893 * _S1131;
    float _S1140 = _S884 * _S1131 + r2_7 * _S1137;
    float _S1141 = _S893 * _S1140;
    float _S1142 = _S885 * _S1131 + r2_7 * _S1140;
    float _S1143 = _S862 * _S1125.x;
    float _S1144 = _S892 * _S1125.x;
    float _S1145 = v_7 * _S1143;
    float _S1146 = _S838.primal_0 * _S1143;
    float _S1147 = _S887 * _S1125.y;
    float _S1148 = _S838.primal_0 * _S1125.y;
    float _S1149 = 2.0f * _S1125.x;
    float _S1150 = _S891 * _S1149;
    float _S1151 = _S891 * _S1125.x;
    float _S1152 = _S1131 + u_7 * _S1149 + _S888 * _S1125.x;
    float _S1153 = p2_2 * _S1152;
    float _S1154 = _S838.primal_0 * _S1152;
    float _S1155 = sx1_2 * _S1131;
    float _S1156 = _S838.primal_0 * _S1131;
    float2  _S1157 = _S886 * _S1125;
    float2  _S1158 = _S889 * _S1125;
    float2  _S1159 = _S843;
    *&((&_S1159)->y) = _S1142;
    *&((&_S1159)->x) = _S1142;
    float _S1160 = _S1158.x + _S1158.y;
    float _S1161 = _S1139 + r2_7 * _S1160;
    float _S1162 = _S1136 + r2_7 * _S1161;
    float _S1163 = _S1134 + r2_7 * _S1162;
    float _S1164 = _S1135 + _S1138 + _S1141 + _S885 * _S1160 + _S884 * _S1161 + _S883 * _S1162 + k4_2 * _S1163;
    float _S1165 = v_7 * _S1164;
    float _S1166 = u_7 * _S1164;
    float2  _S1167 = _S1124 + _S1100.differential_0.primal_0;
    float _S1168 = _S1144 + u_7 * _S1148;
    float _S1169 = _S1161 + _S1064;
    float _S1170 = _S1163 + _S1066;
    float _S1171 = _S1145 + _S1147 + _S1153 + _S1155 + (_S1157 + _S882 * _S1159).x;
    float2  _S1172 = _S889 * _S1159 + make_float2 (_S1128 + _S1150 + 2.0f * _S1151 + _S862 * _S1148 + _S1166 + _S1166, _S1130 + _S1146 + _S1165 + _S1165);
    float _S1173 = _S1162 + _S1065;
    float _S1174 = _S1133 + r2_7 * _S1163 + _S1074;
    float2  _S1175 = _S845 * _S1172;
    float _S1176 = _S1122.x + _S1122.y + _S1175.x + _S1175.y;
    float2  _S1177 = _S881 * _S1172 + _S1167;
    if(_S874)
    {
        float _S1178 = _S847 * _S904;
        float _S1179 = _S1176 / _S876;
        float _S1180 = _S873 * (0.3333333432674408f * - (_S905 + _S847 * _S1179));
        float _S1181 = _S1178 + _S1178 + _S877 * - _S1179 + _S908;
        k_2 = _S1180 + _S1180 + _S907;
        _S875 = _S1181;
        _S876 = _S906;
    }
    else
    {
        float _S1182 = _S846 * k_3;
        float _S1183 = _S1176 / _S875;
        float _S1184 = _S1182 + _S1182 + _S873 * - _S1183 + _S906;
        k_2 = _S846 * _S1183 + _S907;
        _S875 = _S908;
        _S876 = _S1184;
    }
    DiffPair_float_0 _S1185;
    (&_S1185)->primal_0 = _S846;
    (&_S1185)->differential_0 = 0.0f;
    DiffPair_float_0 _S1186;
    (&_S1186)->primal_0 = _S847;
    (&_S1186)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1185, &_S1186, k_2);
    float _S1187 = _S1186.differential_0 + _S875;
    float _S1188 = _S1185.differential_0 + _S876;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1189;
    (&_S1189)->primal_0 = _S845;
    (&_S1189)->differential_0 = _S843;
    s_bwd_length_impl_0(&_S1189, _S1188);
    float2  _S1190 = _S1189.differential_0 + _S1177;
    float3  _S1191 = make_float3 (_S1190.x, _S1190.y, _S1187);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1192;
    (&_S1192)->primal_0 = _S845;
    (&_S1192)->differential_0 = _S843;
    s_bwd_length_impl_0(&_S1192, 0.0f);
    float _S1193 = _S839.primal_0 * dpmean2d_1.y;
    float _S1194 = _S838.primal_0 * dpmean2d_1.x;
    float2  _S1195 = make_float2 (_S1194, _S1193);
    float2  _S1196 = _S857 * _S1195;
    float _S1197 = p1_2 * _S1193;
    float _S1198 = v_6 * _S1193;
    float _S1199 = p2_2 * _S1194;
    float _S1200 = v_6 * _S1194;
    float _S1201 = _S1196.x + _S1196.y;
    float _S1202 = r2_6 * _S1201;
    float _S1203 = r2_6 * _S1202;
    float _S1204 = r2_6 * _S1203;
    float _S1205 = sy1_2 * _S1193 + _S1197 + sx1_2 * _S1194 + _S1199 + _S860 * _S1201 + _S859 * _S1202 + _S858 * _S1203 + k4_2 * _S1204;
    float _S1206 = v_6 * _S1205;
    float _S1207 = u_6 * _S1205;
    float2  _S1208 = make_float2 (r2_6 * _S1194 + _S1156, r2_6 * _S1193 + _S1059);
    float2  _S1209 = make_float2 (_S869 * _S1193 + 2.0f * (u_6 * _S1200 + _S1168) + _S1057, 2.0f * (u_6 * _S1198 + _S1073) + _S865 * _S1194 + _S1154);
    float4  _S1210 = make_float4 (_S1202 + _S1169, _S1203 + _S1173, _S1204 + _S1170, r2_6 * _S1204 + _S1174);
    float3  _S1211 = _S1191 + make_float3 (_S1192.differential_0.x, _S1192.differential_0.y, 0.0f) + _S1093;
    float _S1212 = _S871 * dpmean2d_1.x + _S1171;
    float _S1213 = _S872 * dpmean2d_1.y + _S1072;
    float2  _S1214 = _S861 * _S1195 + make_float2 (_S866 * _S1198 + _S864 * _S1199 + 2.0f * (u_6 * _S1199) + _S862 * _S1200 + _S1207 + _S1207, _S868 * _S1197 + 2.0f * (v_6 * _S1197) + _S867 * _S1193 + _S863 * _S1194 + _S1206 + _S1206);
    CameraDistortion_0 _S1215 = _S987;
    (&_S1215)->thin_prism_coeffs_0 = _S1208;
    (&_S1215)->tangential_coeffs_0 = _S1209;
    (&_S1215)->radial_coeffs_0 = _S1210;
    CameraDistortion_0 _S1216 = _S987;
    CameraDistortion_0 _S1217 = _S1215;
    CameraDistortion_0 _S1218 = CameraDistortion_x24_syn_dadd_0(&_S1216, &_S1217);
    float2  _S1219 = _S845 * _S1214;
    float2  _S1220 = _S856 * _S1214;
    float _S1221 = _S1219.x + _S1219.y;
    if(_S849)
    {
        float _S1222 = _S1221 / _S851;
        float _S1223 = _S852 * - _S1222;
        float _S1224 = _S848 * (0.3333333432674408f * - (_S847 * _S1222));
        k_2 = _S1224 + _S1224;
        _S850 = _S1223;
        _S851 = 0.0f;
    }
    else
    {
        float _S1225 = _S1221 / _S850;
        float _S1226 = _S848 * - _S1225;
        k_2 = _S846 * _S1225;
        _S850 = 0.0f;
        _S851 = _S1226;
    }
    DiffPair_float_0 _S1227;
    (&_S1227)->primal_0 = _S846;
    (&_S1227)->differential_0 = 0.0f;
    DiffPair_float_0 _S1228;
    (&_S1228)->primal_0 = _S847;
    (&_S1228)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1227, &_S1228, k_2);
    float _S1229 = _S1228.differential_0 + _S850;
    float _S1230 = _S1227.differential_0 + _S851;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1231;
    (&_S1231)->primal_0 = _S845;
    (&_S1231)->differential_0 = _S843;
    s_bwd_length_impl_0(&_S1231, _S1230);
    float2  _S1232 = _S1231.differential_0 + _S1220;
    float4  _S1233 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1234;
    (&_S1234)->primal_0 = _S840.primal_0;
    (&_S1234)->differential_0 = _S1233;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1235;
    (&_S1235)->primal_0 = _S841.primal_0;
    (&_S1235)->differential_0 = _S843;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1236;
    (&_S1236)->primal_0 = _S842.primal_0;
    (&_S1236)->differential_0 = _S843;
    CameraDistortion_0 _S1237 = _S1218;
    s_bwd_prop_CameraDistortion_x24init_0(&_S1234, &_S1235, &_S1236, &_S1237);
    dpthin_prism_coeffs_3->primal_0 = (*dpthin_prism_coeffs_3).primal_0;
    dpthin_prism_coeffs_3->differential_0 = _S1236.differential_0;
    dptangential_coeffs_3->primal_0 = (*dptangential_coeffs_3).primal_0;
    dptangential_coeffs_3->differential_0 = _S1235.differential_0;
    dpradial_coeffs_3->primal_0 = (*dpradial_coeffs_3).primal_0;
    dpradial_coeffs_3->differential_0 = _S1234.differential_0;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = dpmean2d_1.y;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = dpmean2d_1.x;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S1213;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S1212;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S995.differential_0;
    float3  _S1238 = _S1211 + make_float3 (_S1232.x, _S1232.y, _S1229);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1238;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_11, float fy_11, float cx_11, float cy_11, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, uint image_width_7, uint image_height_7, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    Matrix<float, 2, 2>  _S1239 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1240 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1241 = { _S1240, _S1240 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1242 = { _S1240, _S1240, _S1241, _S1240, _S1240, _S1241 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1243;
    (&_S1243)->_S590 = _S1239;
    (&_S1243)->_S591 = _S1242;
    float3  mean_c_7 = s_primal_ctx_mul_0(R_9, mean_7) + t_8;
    float3  _S1244 = s_primal_ctx_exp_0(scale_9);
    float _S1245 = quat_10.y;
    float _S1246 = _S1245 * _S1245 + quat_10.z * quat_10.z + quat_10.w * quat_10.w + quat_10.x * quat_10.x;
    float _S1247 = s_primal_ctx_rsqrt_0(_S1246);
    float x_35 = quat_10.y * _S1247;
    float y_20 = quat_10.z * _S1247;
    float z_17 = quat_10.w * _S1247;
    float w_10 = quat_10.x * _S1247;
    float x2_10 = x_35 * x_35;
    float y2_10 = y_20 * y_20;
    float z2_17 = z_17 * z_17;
    float xy_10 = x_35 * y_20;
    float xz_10 = x_35 * z_17;
    float yz_10 = y_20 * z_17;
    float wx_10 = w_10 * x_35;
    float wy_10 = w_10 * y_20;
    float wz_10 = w_10 * z_17;
    Matrix<float, 3, 3>  _S1248 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1244.x, 0.0f, 0.0f, 0.0f, _S1244.y, 0.0f, 0.0f, 0.0f, _S1244.z);
    Matrix<float, 3, 3>  _S1249 = s_primal_ctx_mul_1(_S1248, S_1);
    Matrix<float, 3, 3>  _S1250 = transpose_0(_S1249);
    Matrix<float, 3, 3>  _S1251 = s_primal_ctx_mul_1(_S1249, _S1250);
    Matrix<float, 3, 3>  _S1252 = s_primal_ctx_mul_1(R_9, _S1251);
    Matrix<float, 3, 3>  _S1253 = transpose_0(R_9);
    Matrix<float, 3, 3>  _S1254 = s_primal_ctx_mul_1(_S1252, _S1253);
    Matrix<float, 2, 2>  _S1255 = _S1239;
    float2  _S1256 = make_float2 (0.0f);
    float2  _S1257 = _S1256;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_7, _S1254, fx_11, fy_11, cx_11, cy_11, radial_coeffs_10, tangential_coeffs_10, thin_prism_coeffs_10, &_S1255, &_S1257, &(&_S1243)->_S591);
    (&_S1243)->_S590 = _S1255;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1258 = _S1243;
    float _S1259 = _S1243._S590.rows[int(0)].y * _S1243._S590.rows[int(1)].x;
    float det_orig_8 = _S1243._S590.rows[int(0)].x * _S1243._S590.rows[int(1)].y - _S1259;
    float _S1260 = _S1243._S590.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1261 = _S1243._S590;
    *&(((&_S1261)->rows + (int(0)))->x) = _S1260;
    float _S1262 = _S1243._S590.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1261)->rows + (int(1)))->y) = _S1262;
    Matrix<float, 2, 2>  _S1263 = _S1261;
    Matrix<float, 2, 2>  _S1264 = _S1261;
    float det_blur_5 = _S1260 * _S1262 - _S1259;
    float _S1265 = det_orig_8 / det_blur_5;
    float _S1266 = det_blur_5 * det_blur_5;
    float _S1267 = s_primal_ctx_max_0(0.0f, _S1265);
    float _S1268 = s_primal_ctx_sqrt_0(_S1267);
    float invdet_8 = 1.0f / det_blur_5;
    float _S1269 = - _S1243._S590.rows[int(0)].y;
    float _S1270 = - _S1243._S590.rows[int(1)].x;
    float _S1271 = - in_opacity_7;
    float _S1272 = 1.0f + s_primal_ctx_exp_1(_S1271);
    float _S1273 = 1.0f / _S1272;
    float _S1274 = _S1272 * _S1272;
    float _S1275;
    if(antialiased_7)
    {
        _S1275 = _S1273 * _S1268;
    }
    else
    {
        _S1275 = _S1273;
    }
    float _S1276 = _S1275 / 0.00392156885936856f;
    float _S1277 = 2.0f * s_primal_ctx_log_0(_S1276);
    float _S1278 = s_primal_ctx_sqrt_0(_S1277);
    float _S1279 = _S1263.rows[int(0)].x;
    float _S1280 = _S1264.rows[int(1)].y;
    float3  _S1281 = mean_7 - - s_primal_ctx_mul_0(_S1253, t_8);
    float _S1282 = _S1281.x;
    float _S1283 = _S1281.y;
    float _S1284 = _S1281.z;
    float _S1285 = _S1282 * _S1282 + _S1283 * _S1283 + _S1284 * _S1284;
    float _S1286 = s_primal_ctx_sqrt_0(_S1285);
    float x_36 = _S1282 / _S1286;
    float3  _S1287 = make_float3 (x_36);
    float _S1288 = _S1286 * _S1286;
    float y_21 = _S1283 / _S1286;
    float z_18 = _S1284 / _S1286;
    float3  _S1289 = make_float3 (z_18);
    float _S1290 = - y_21;
    float3  _S1291 = make_float3 (_S1290);
    float z2_18 = z_18 * z_18;
    float fTmp0B_7 = -1.09254848957061768f * z_18;
    float fC1_7 = x_36 * x_36 - y_21 * y_21;
    float _S1292 = 2.0f * x_36;
    float fS1_7 = _S1292 * y_21;
    float pSH6_1 = 0.94617468118667603f * z2_18 - 0.31539157032966614f;
    float3  _S1293 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_7 * x_36;
    float3  _S1294 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_7 * y_21;
    float3  _S1295 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_7;
    float3  _S1296 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_7;
    float3  _S1297 = make_float3 (pSH4_1);
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_18;
    float _S1298 = 1.86588168144226074f * z2_18 - 1.11952900886535645f;
    float pSH12_1 = z_18 * _S1298;
    float3  _S1299 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_7 * x_36;
    float3  _S1300 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_7 * y_21;
    float3  _S1301 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_7 * fC1_7;
    float3  _S1302 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_7 * fS1_7;
    float3  _S1303 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_36 * fC1_7 - y_21 * fS1_7);
    float3  _S1304 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_36 * fS1_7 + y_21 * fC1_7);
    float3  _S1305 = make_float3 (pSH9_1);
    float3  _S1306 = make_float3 (0.0f);
    float3  _S1307 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1308;
    (&_S1308)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1290) * (*sh_coeffs_7)[int(1)] + make_float3 (z_18) * (*sh_coeffs_7)[int(2)] - make_float3 (x_36) * (*sh_coeffs_7)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_7)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_7)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_7)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_7)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_7)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_7)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_7)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_7)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_7)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_7)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_7)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f);
    (&_S1308)->differential_0 = _S1307;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1309;
    (&_S1309)->primal_0 = _S1306;
    (&_S1309)->differential_0 = _S1307;
    s_bwd_prop_max_0(&_S1308, &_S1309, v_rgb_1);
    float3  _S1310 = _S1304 * _S1308.differential_0;
    float3  _S1311 = (*sh_coeffs_7)[int(15)] * _S1308.differential_0;
    float3  _S1312 = _S1302 * _S1308.differential_0;
    float3  _S1313 = (*sh_coeffs_7)[int(14)] * _S1308.differential_0;
    float3  _S1314 = _S1300 * _S1308.differential_0;
    float3  _S1315 = (*sh_coeffs_7)[int(13)] * _S1308.differential_0;
    float3  _S1316 = _S1299 * _S1308.differential_0;
    float3  _S1317 = (*sh_coeffs_7)[int(12)] * _S1308.differential_0;
    float3  _S1318 = _S1301 * _S1308.differential_0;
    float3  _S1319 = (*sh_coeffs_7)[int(11)] * _S1308.differential_0;
    float3  _S1320 = _S1303 * _S1308.differential_0;
    float3  _S1321 = (*sh_coeffs_7)[int(10)] * _S1308.differential_0;
    float3  _S1322 = _S1305 * _S1308.differential_0;
    float3  _S1323 = (*sh_coeffs_7)[int(9)] * _S1308.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1323.x + _S1323.y + _S1323.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1311.x + _S1311.y + _S1311.z);
    float _S1324 = _S1321.x + _S1321.y + _S1321.z;
    float _S1325 = _S1313.x + _S1313.y + _S1313.z;
    float _S1326 = _S1319.x + _S1319.y + _S1319.z;
    float _S1327 = _S1315.x + _S1315.y + _S1315.z;
    float _S1328 = _S1317.x + _S1317.y + _S1317.z;
    float _S1329 = - s_diff_fC2_T_1;
    float3  _S1330 = _S1296 * _S1308.differential_0;
    float3  _S1331 = (*sh_coeffs_7)[int(8)] * _S1308.differential_0;
    float3  _S1332 = _S1294 * _S1308.differential_0;
    float3  _S1333 = (*sh_coeffs_7)[int(7)] * _S1308.differential_0;
    float3  _S1334 = _S1293 * _S1308.differential_0;
    float3  _S1335 = (*sh_coeffs_7)[int(6)] * _S1308.differential_0;
    float3  _S1336 = _S1295 * _S1308.differential_0;
    float3  _S1337 = (*sh_coeffs_7)[int(5)] * _S1308.differential_0;
    float3  _S1338 = _S1297 * _S1308.differential_0;
    float3  _S1339 = (*sh_coeffs_7)[int(4)] * _S1308.differential_0;
    float _S1340 = _S1337.x + _S1337.y + _S1337.z;
    float _S1341 = _S1333.x + _S1333.y + _S1333.z;
    float _S1342 = fTmp1B_7 * _S1324 + x_36 * s_diff_fS2_T_1 + y_21 * _S1329 + 0.54627424478530884f * (_S1339.x + _S1339.y + _S1339.z);
    float _S1343 = fTmp1B_7 * _S1325 + y_21 * s_diff_fS2_T_1 + x_36 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1331.x + _S1331.y + _S1331.z);
    float _S1344 = y_21 * - _S1343;
    float _S1345 = x_36 * _S1343;
    float _S1346 = z_18 * (1.86588168144226074f * (z_18 * _S1328) + -2.28522896766662598f * (y_21 * _S1326 + x_36 * _S1327) + 0.94617468118667603f * (_S1335.x + _S1335.y + _S1335.z));
    float3  _S1347 = make_float3 (0.48860251903533936f) * _S1308.differential_0;
    float3  _S1348 = - _S1347;
    float3  _S1349 = _S1287 * _S1348;
    float3  _S1350 = (*sh_coeffs_7)[int(3)] * _S1348;
    float3  _S1351 = _S1289 * _S1347;
    float3  _S1352 = (*sh_coeffs_7)[int(2)] * _S1347;
    float3  _S1353 = _S1291 * _S1347;
    float3  _S1354 = (*sh_coeffs_7)[int(1)] * _S1347;
    float _S1355 = (_S1298 * _S1328 + 1.44530570507049561f * (fS1_7 * _S1324 + fC1_7 * _S1325) + -1.09254848957061768f * (y_21 * _S1340 + x_36 * _S1341) + _S1346 + _S1346 + _S1352.x + _S1352.y + _S1352.z) / _S1288;
    float _S1356 = _S1286 * _S1355;
    float _S1357 = (fTmp0C_7 * _S1326 + fC1_7 * s_diff_fS2_T_1 + fS1_7 * _S1329 + fTmp0B_7 * _S1340 + _S1292 * _S1342 + _S1344 + _S1344 + - (_S1354.x + _S1354.y + _S1354.z)) / _S1288;
    float _S1358 = _S1286 * _S1357;
    float _S1359 = (fTmp0C_7 * _S1327 + fS1_7 * s_diff_fS2_T_1 + fC1_7 * s_diff_fC2_T_1 + fTmp0B_7 * _S1341 + 2.0f * (y_21 * _S1342) + _S1345 + _S1345 + _S1350.x + _S1350.y + _S1350.z) / _S1288;
    float _S1360 = _S1286 * _S1359;
    float _S1361 = _S1284 * - _S1355 + _S1283 * - _S1357 + _S1282 * - _S1359;
    DiffPair_float_0 _S1362;
    (&_S1362)->primal_0 = _S1285;
    (&_S1362)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1362, _S1361);
    float _S1363 = _S1284 * _S1362.differential_0;
    float _S1364 = _S1283 * _S1362.differential_0;
    float _S1365 = _S1282 * _S1362.differential_0;
    float3  _S1366 = make_float3 (0.282094806432724f) * _S1308.differential_0;
    float3  _S1367 = make_float3 (_S1360 + _S1365 + _S1365, _S1358 + _S1364 + _S1364, _S1356 + _S1363 + _S1363);
    float3  _S1368 = - - _S1367;
    Matrix<float, 3, 3>  _S1369 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1370;
    (&_S1370)->primal_0 = _S1253;
    (&_S1370)->differential_0 = _S1369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1371;
    (&_S1371)->primal_0 = t_8;
    (&_S1371)->differential_0 = _S1307;
    s_bwd_prop_mul_0(&_S1370, &_S1371, _S1368);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1372 = _S1370;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1373 = _S1371;
    float2  _S1374 = _S1256;
    *&((&_S1374)->y) = v_conic_1.z;
    float2  _S1375 = _S1256;
    *&((&_S1375)->y) = v_conic_1.y;
    *&((&_S1375)->x) = v_conic_1.x;
    DiffPair_float_0 _S1376;
    (&_S1376)->primal_0 = _S1280;
    (&_S1376)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1376, 0.0f);
    DiffPair_float_0 _S1377;
    (&_S1377)->primal_0 = _S1279;
    (&_S1377)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1377, 0.0f);
    DiffPair_float_0 _S1378;
    (&_S1378)->primal_0 = 3.32999992370605469f;
    (&_S1378)->differential_0 = 0.0f;
    DiffPair_float_0 _S1379;
    (&_S1379)->primal_0 = _S1278;
    (&_S1379)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1378, &_S1379, 0.0f);
    DiffPair_float_0 _S1380;
    (&_S1380)->primal_0 = _S1277;
    (&_S1380)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1380, _S1379.differential_0);
    float _S1381 = 2.0f * _S1380.differential_0;
    DiffPair_float_0 _S1382;
    (&_S1382)->primal_0 = _S1276;
    (&_S1382)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1382, _S1381);
    float2  _S1383 = make_float2 (_S1377.differential_0, 0.0f);
    float _S1384 = v_opacity_1 + 254.9999847412109375f * _S1382.differential_0;
    FixedArray<float3 , 16>  _S1385;
    _S1385[int(0)] = _S1307;
    _S1385[int(1)] = _S1307;
    _S1385[int(2)] = _S1307;
    _S1385[int(3)] = _S1307;
    _S1385[int(4)] = _S1307;
    _S1385[int(5)] = _S1307;
    _S1385[int(6)] = _S1307;
    _S1385[int(7)] = _S1307;
    _S1385[int(8)] = _S1307;
    _S1385[int(9)] = _S1307;
    _S1385[int(10)] = _S1307;
    _S1385[int(11)] = _S1307;
    _S1385[int(12)] = _S1307;
    _S1385[int(13)] = _S1307;
    _S1385[int(14)] = _S1307;
    _S1385[int(15)] = _S1307;
    _S1385[int(7)] = _S1332;
    _S1385[int(0)] = _S1366;
    _S1385[int(1)] = _S1353;
    _S1385[int(2)] = _S1351;
    _S1385[int(3)] = _S1349;
    _S1385[int(4)] = _S1338;
    _S1385[int(5)] = _S1336;
    _S1385[int(6)] = _S1334;
    _S1385[int(15)] = _S1310;
    _S1385[int(8)] = _S1330;
    _S1385[int(9)] = _S1322;
    _S1385[int(10)] = _S1320;
    _S1385[int(11)] = _S1318;
    _S1385[int(12)] = _S1316;
    _S1385[int(13)] = _S1314;
    _S1385[int(14)] = _S1312;
    float3  _S1386 = _S1385[int(0)];
    float3  _S1387 = _S1385[int(1)];
    float3  _S1388 = _S1385[int(2)];
    float3  _S1389 = _S1385[int(3)];
    float3  _S1390 = _S1385[int(4)];
    float3  _S1391 = _S1385[int(5)];
    float3  _S1392 = _S1385[int(6)];
    float3  _S1393 = _S1385[int(7)];
    float3  _S1394 = _S1385[int(8)];
    float3  _S1395 = _S1385[int(9)];
    float3  _S1396 = _S1385[int(10)];
    float3  _S1397 = _S1385[int(11)];
    float3  _S1398 = _S1385[int(12)];
    float3  _S1399 = _S1385[int(13)];
    float3  _S1400 = _S1385[int(14)];
    float3  _S1401 = _S1385[int(15)];
    Matrix<float, 2, 2>  _S1402 = _S1239;
    _S1402[int(1)] = _S1374;
    _S1402[int(0)] = _S1375;
    Matrix<float, 2, 2>  _S1403 = _S1402;
    float2  _S1404 = make_float2 (0.0f, _S1376.differential_0);
    float _S1405;
    if(antialiased_7)
    {
        float _S1406 = _S1273 * _S1384;
        _S1275 = _S1268 * _S1384;
        _S1405 = _S1406;
    }
    else
    {
        _S1275 = _S1384;
        _S1405 = 0.0f;
    }
    float _S1407 = - (_S1275 / _S1274);
    DiffPair_float_0 _S1408;
    (&_S1408)->primal_0 = _S1271;
    (&_S1408)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1408, _S1407);
    float _S1409 = - _S1408.differential_0;
    float _S1410 = invdet_8 * _S1403.rows[int(1)].y;
    float _S1411 = - (invdet_8 * _S1403.rows[int(1)].x);
    float _S1412 = - (invdet_8 * _S1403.rows[int(0)].y);
    float _S1413 = invdet_8 * _S1403.rows[int(0)].x;
    float _S1414 = - ((_S1260 * _S1403.rows[int(1)].y + _S1270 * _S1403.rows[int(1)].x + _S1269 * _S1403.rows[int(0)].y + _S1262 * _S1403.rows[int(0)].x) / _S1266);
    DiffPair_float_0 _S1415;
    (&_S1415)->primal_0 = _S1267;
    (&_S1415)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1415, _S1405);
    DiffPair_float_0 _S1416;
    (&_S1416)->primal_0 = 0.0f;
    (&_S1416)->differential_0 = 0.0f;
    DiffPair_float_0 _S1417;
    (&_S1417)->primal_0 = _S1265;
    (&_S1417)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1416, &_S1417, _S1415.differential_0);
    float _S1418 = _S1417.differential_0 / _S1266;
    float s_diff_det_orig_T_1 = det_blur_5 * _S1418;
    float _S1419 = _S1414 + det_orig_8 * - _S1418;
    float _S1420 = - _S1419;
    float _S1421 = _S1260 * _S1419;
    float _S1422 = _S1262 * _S1419;
    Matrix<float, 2, 2>  _S1423 = _S1239;
    _S1423[int(1)] = _S1404;
    _S1423[int(0)] = _S1383;
    _S1261 = _S1423;
    *&(((&_S1261)->rows + (int(1)))->y) = 0.0f;
    float _S1424 = _S1413 + _S1421 + _S1423.rows[int(1)].y;
    *&(((&_S1261)->rows + (int(0)))->x) = 0.0f;
    float _S1425 = _S1410 + _S1422 + _S1423.rows[int(0)].x;
    float _S1426 = _S1420 + - s_diff_det_orig_T_1;
    float _S1427 = _S1411 + _S1258._S590.rows[int(0)].y * _S1426;
    float _S1428 = _S1412 + _S1258._S590.rows[int(1)].x * _S1426;
    float _S1429 = _S1258._S590.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S1430 = _S1424 + _S1258._S590.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S1431 = _S1256;
    *&((&_S1431)->x) = _S1427;
    *&((&_S1431)->y) = _S1430;
    float _S1432 = _S1425 + _S1429;
    float2  _S1433 = _S1256;
    *&((&_S1433)->y) = _S1428;
    *&((&_S1433)->x) = _S1432;
    Matrix<float, 2, 2>  _S1434 = _S1239;
    _S1434[int(1)] = _S1431;
    _S1434[int(0)] = _S1433;
    Matrix<float, 2, 2>  _S1435 = _S1261 + _S1434;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1436;
    (&_S1436)->primal_0 = mean_c_7;
    (&_S1436)->differential_0 = _S1307;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1437;
    (&_S1437)->primal_0 = _S1254;
    (&_S1437)->differential_0 = _S1369;
    DiffPair_float_0 _S1438;
    (&_S1438)->primal_0 = fx_11;
    (&_S1438)->differential_0 = 0.0f;
    DiffPair_float_0 _S1439;
    (&_S1439)->primal_0 = fy_11;
    (&_S1439)->differential_0 = 0.0f;
    DiffPair_float_0 _S1440;
    (&_S1440)->primal_0 = cx_11;
    (&_S1440)->differential_0 = 0.0f;
    DiffPair_float_0 _S1441;
    (&_S1441)->primal_0 = cy_11;
    (&_S1441)->differential_0 = 0.0f;
    float4  _S1442 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1443;
    (&_S1443)->primal_0 = radial_coeffs_10;
    (&_S1443)->differential_0 = _S1442;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1444;
    (&_S1444)->primal_0 = tangential_coeffs_10;
    (&_S1444)->differential_0 = _S1256;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1445;
    (&_S1445)->primal_0 = thin_prism_coeffs_10;
    (&_S1445)->differential_0 = _S1256;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1446 = _S1258._S591;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S1436, &_S1437, &_S1438, &_S1439, &_S1440, &_S1441, &_S1443, &_S1444, &_S1445, _S1435, v_mean2d_1, &_S1446);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1447;
    (&_S1447)->primal_0 = _S1252;
    (&_S1447)->differential_0 = _S1369;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1448;
    (&_S1448)->primal_0 = _S1253;
    (&_S1448)->differential_0 = _S1369;
    s_bwd_prop_mul_3(&_S1447, &_S1448, _S1437.differential_0);
    Matrix<float, 3, 3>  _S1449 = transpose_0(_S1448.differential_0 + _S1372.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1450;
    (&_S1450)->primal_0 = R_9;
    (&_S1450)->differential_0 = _S1369;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1451;
    (&_S1451)->primal_0 = _S1251;
    (&_S1451)->differential_0 = _S1369;
    s_bwd_prop_mul_3(&_S1450, &_S1451, _S1447.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1452;
    (&_S1452)->primal_0 = _S1249;
    (&_S1452)->differential_0 = _S1369;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1453;
    (&_S1453)->primal_0 = _S1250;
    (&_S1453)->differential_0 = _S1369;
    s_bwd_prop_mul_3(&_S1452, &_S1453, _S1451.differential_0);
    Matrix<float, 3, 3>  _S1454 = _S1452.differential_0 + transpose_0(_S1453.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1455;
    (&_S1455)->primal_0 = _S1248;
    (&_S1455)->differential_0 = _S1369;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1456;
    (&_S1456)->primal_0 = S_1;
    (&_S1456)->differential_0 = _S1369;
    s_bwd_prop_mul_3(&_S1455, &_S1456, _S1454);
    Matrix<float, 3, 3>  _S1457 = transpose_0(_S1455.differential_0);
    float _S1458 = 2.0f * - _S1457.rows[int(2)].z;
    float _S1459 = 2.0f * _S1457.rows[int(2)].y;
    float _S1460 = 2.0f * _S1457.rows[int(2)].x;
    float _S1461 = 2.0f * _S1457.rows[int(1)].z;
    float _S1462 = 2.0f * - _S1457.rows[int(1)].y;
    float _S1463 = 2.0f * _S1457.rows[int(1)].x;
    float _S1464 = 2.0f * _S1457.rows[int(0)].z;
    float _S1465 = 2.0f * _S1457.rows[int(0)].y;
    float _S1466 = 2.0f * - _S1457.rows[int(0)].x;
    float _S1467 = - _S1463 + _S1465;
    float _S1468 = _S1460 + - _S1464;
    float _S1469 = - _S1459 + _S1461;
    float _S1470 = _S1459 + _S1461;
    float _S1471 = _S1460 + _S1464;
    float _S1472 = _S1463 + _S1465;
    float _S1473 = z_17 * (_S1462 + _S1466);
    float _S1474 = y_20 * (_S1458 + _S1466);
    float _S1475 = x_35 * (_S1458 + _S1462);
    float _S1476 = z_17 * _S1467 + y_20 * _S1468 + x_35 * _S1469;
    float _S1477 = _S1247 * _S1476;
    float _S1478 = w_10 * _S1467 + y_20 * _S1470 + x_35 * _S1471 + _S1473 + _S1473;
    float _S1479 = _S1247 * _S1478;
    float _S1480 = w_10 * _S1468 + z_17 * _S1470 + x_35 * _S1472 + _S1474 + _S1474;
    float _S1481 = _S1247 * _S1480;
    float _S1482 = w_10 * _S1469 + z_17 * _S1471 + y_20 * _S1472 + _S1475 + _S1475;
    float _S1483 = _S1247 * _S1482;
    float _S1484 = quat_10.x * _S1476 + quat_10.w * _S1478 + quat_10.z * _S1480 + quat_10.y * _S1482;
    DiffPair_float_0 _S1485;
    (&_S1485)->primal_0 = _S1246;
    (&_S1485)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1485, _S1484);
    float _S1486 = quat_10.x * _S1485.differential_0;
    float _S1487 = quat_10.w * _S1485.differential_0;
    float _S1488 = quat_10.z * _S1485.differential_0;
    float _S1489 = quat_10.y * _S1485.differential_0;
    float _S1490 = _S1479 + _S1487 + _S1487;
    float _S1491 = _S1481 + _S1488 + _S1488;
    float _S1492 = _S1483 + _S1489 + _S1489;
    float _S1493 = _S1477 + _S1486 + _S1486;
    float3  _S1494 = _S1307;
    *&((&_S1494)->z) = _S1456.differential_0.rows[int(2)].z;
    *&((&_S1494)->y) = _S1456.differential_0.rows[int(1)].y;
    *&((&_S1494)->x) = _S1456.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1495;
    (&_S1495)->primal_0 = scale_9;
    (&_S1495)->differential_0 = _S1307;
    s_bwd_prop_exp_1(&_S1495, _S1494);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1496 = _S1495;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1497;
    (&_S1497)->primal_0 = mean_c_7;
    (&_S1497)->differential_0 = _S1307;
    s_bwd_length_impl_1(&_S1497, v_depth_1);
    float3  _S1498 = _S1436.differential_0 + _S1497.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1499;
    (&_S1499)->primal_0 = R_9;
    (&_S1499)->differential_0 = _S1369;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1500;
    (&_S1500)->primal_0 = mean_7;
    (&_S1500)->differential_0 = _S1307;
    s_bwd_prop_mul_0(&_S1499, &_S1500, _S1498);
    float3  _S1501 = _S1498 + _S1373.differential_0;
    Matrix<float, 3, 3>  _S1502 = _S1449 + _S1450.differential_0 + _S1499.differential_0;
    float4  _S1503 = _S1442;
    *&((&_S1503)->w) = _S1490;
    *&((&_S1503)->z) = _S1491;
    *&((&_S1503)->y) = _S1492;
    *&((&_S1503)->x) = _S1493;
    float4  _S1504 = _S1503;
    float3  _S1505 = _S1500.differential_0 + _S1367;
    *v_mean_1 = _S1505;
    *v_quat_1 = _S1504;
    *v_scale_1 = _S1496.differential_0;
    *v_in_opacity_1 = _S1409;
    (*v_sh_coeffs_1)[int(0)] = _S1386;
    (*v_sh_coeffs_1)[int(1)] = _S1387;
    (*v_sh_coeffs_1)[int(2)] = _S1388;
    (*v_sh_coeffs_1)[int(3)] = _S1389;
    (*v_sh_coeffs_1)[int(4)] = _S1390;
    (*v_sh_coeffs_1)[int(5)] = _S1391;
    (*v_sh_coeffs_1)[int(6)] = _S1392;
    (*v_sh_coeffs_1)[int(7)] = _S1393;
    (*v_sh_coeffs_1)[int(8)] = _S1394;
    (*v_sh_coeffs_1)[int(9)] = _S1395;
    (*v_sh_coeffs_1)[int(10)] = _S1396;
    (*v_sh_coeffs_1)[int(11)] = _S1397;
    (*v_sh_coeffs_1)[int(12)] = _S1398;
    (*v_sh_coeffs_1)[int(13)] = _S1399;
    (*v_sh_coeffs_1)[int(14)] = _S1400;
    (*v_sh_coeffs_1)[int(15)] = _S1401;
    *v_R_1 = _S1502;
    *v_t_1 = _S1501;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_12, float fy_12, float cx_12, float cy_12, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, uint image_width_8, uint image_height_8, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    float3  mean_c_8 = s_primal_ctx_mul_0(R_10, mean_8) + t_9;
    float3  _S1506 = s_primal_ctx_exp_0(scale_10);
    float _S1507 = quat_11.y;
    float _S1508 = _S1507 * _S1507 + quat_11.z * quat_11.z + quat_11.w * quat_11.w + quat_11.x * quat_11.x;
    float _S1509 = s_primal_ctx_rsqrt_0(_S1508);
    float x_37 = quat_11.y * _S1509;
    float y_22 = quat_11.z * _S1509;
    float z_19 = quat_11.w * _S1509;
    float w_11 = quat_11.x * _S1509;
    float x2_11 = x_37 * x_37;
    float y2_11 = y_22 * y_22;
    float z2_19 = z_19 * z_19;
    float xy_11 = x_37 * y_22;
    float xz_11 = x_37 * z_19;
    float yz_11 = y_22 * z_19;
    float wx_11 = w_11 * x_37;
    float wy_11 = w_11 * y_22;
    float wz_11 = w_11 * z_19;
    Matrix<float, 3, 3>  _S1510 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S1506.x, 0.0f, 0.0f, 0.0f, _S1506.y, 0.0f, 0.0f, 0.0f, _S1506.z);
    Matrix<float, 3, 3>  _S1511 = s_primal_ctx_mul_1(_S1510, S_2);
    Matrix<float, 3, 3>  _S1512 = transpose_0(_S1511);
    Matrix<float, 3, 3>  _S1513 = s_primal_ctx_mul_1(_S1511, _S1512);
    Matrix<float, 3, 3>  _S1514 = s_primal_ctx_mul_1(R_10, _S1513);
    Matrix<float, 3, 3>  _S1515 = transpose_0(R_10);
    Matrix<float, 3, 3>  _S1516 = s_primal_ctx_mul_1(_S1514, _S1515);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_12, 0.0f, 0.0f, 0.0f, fy_12, 0.0f);
    Matrix<float, 2, 3>  _S1517 = s_primal_ctx_mul_2(J_11, _S1516);
    Matrix<float, 3, 2>  _S1518 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S1519 = s_primal_ctx_mul_3(_S1517, _S1518);
    float _S1520 = _S1519.rows[int(0)].y * _S1519.rows[int(1)].x;
    float det_orig_9 = _S1519.rows[int(0)].x * _S1519.rows[int(1)].y - _S1520;
    float _S1521 = _S1519.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1522 = _S1519;
    *&(((&_S1522)->rows + (int(0)))->x) = _S1521;
    float _S1523 = _S1519.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1522)->rows + (int(1)))->y) = _S1523;
    Matrix<float, 2, 2>  _S1524 = _S1522;
    Matrix<float, 2, 2>  _S1525 = _S1522;
    float det_blur_6 = _S1521 * _S1523 - _S1520;
    float _S1526 = det_orig_9 / det_blur_6;
    float _S1527 = det_blur_6 * det_blur_6;
    float _S1528 = s_primal_ctx_max_0(0.0f, _S1526);
    float _S1529 = s_primal_ctx_sqrt_0(_S1528);
    float invdet_9 = 1.0f / det_blur_6;
    float _S1530 = - _S1519.rows[int(0)].y;
    float _S1531 = - _S1519.rows[int(1)].x;
    float _S1532 = - in_opacity_8;
    float _S1533 = 1.0f + s_primal_ctx_exp_1(_S1532);
    float _S1534 = 1.0f / _S1533;
    float _S1535 = _S1533 * _S1533;
    float _S1536;
    if(antialiased_8)
    {
        _S1536 = _S1534 * _S1529;
    }
    else
    {
        _S1536 = _S1534;
    }
    float _S1537 = _S1536 / 0.00392156885936856f;
    float _S1538 = 2.0f * s_primal_ctx_log_0(_S1537);
    float _S1539 = s_primal_ctx_sqrt_0(_S1538);
    float _S1540 = _S1524.rows[int(0)].x;
    float _S1541 = _S1525.rows[int(1)].y;
    float3  _S1542 = mean_8 - - s_primal_ctx_mul_0(_S1515, t_9);
    float _S1543 = _S1542.x;
    float _S1544 = _S1542.y;
    float _S1545 = _S1542.z;
    float _S1546 = _S1543 * _S1543 + _S1544 * _S1544 + _S1545 * _S1545;
    float _S1547 = s_primal_ctx_sqrt_0(_S1546);
    float x_38 = _S1543 / _S1547;
    float3  _S1548 = make_float3 (x_38);
    float _S1549 = _S1547 * _S1547;
    float y_23 = _S1544 / _S1547;
    float z_20 = _S1545 / _S1547;
    float3  _S1550 = make_float3 (z_20);
    float _S1551 = - y_23;
    float3  _S1552 = make_float3 (_S1551);
    float z2_20 = z_20 * z_20;
    float fTmp0B_8 = -1.09254848957061768f * z_20;
    float fC1_8 = x_38 * x_38 - y_23 * y_23;
    float _S1553 = 2.0f * x_38;
    float fS1_8 = _S1553 * y_23;
    float pSH6_2 = 0.94617468118667603f * z2_20 - 0.31539157032966614f;
    float3  _S1554 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_8 * x_38;
    float3  _S1555 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_8 * y_23;
    float3  _S1556 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_8;
    float3  _S1557 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_8;
    float3  _S1558 = make_float3 (pSH4_2);
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_20;
    float _S1559 = 1.86588168144226074f * z2_20 - 1.11952900886535645f;
    float pSH12_2 = z_20 * _S1559;
    float3  _S1560 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_8 * x_38;
    float3  _S1561 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_8 * y_23;
    float3  _S1562 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_8 * fC1_8;
    float3  _S1563 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_8 * fS1_8;
    float3  _S1564 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_38 * fC1_8 - y_23 * fS1_8);
    float3  _S1565 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_38 * fS1_8 + y_23 * fC1_8);
    float3  _S1566 = make_float3 (pSH9_2);
    float3  _S1567 = make_float3 (0.0f);
    float3  _S1568 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1569;
    (&_S1569)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1551) * (*sh_coeffs_8)[int(1)] + make_float3 (z_20) * (*sh_coeffs_8)[int(2)] - make_float3 (x_38) * (*sh_coeffs_8)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_8)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_8)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_8)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_8)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_8)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_8)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_8)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_8)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_8)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_8)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_8)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f);
    (&_S1569)->differential_0 = _S1568;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1570;
    (&_S1570)->primal_0 = _S1567;
    (&_S1570)->differential_0 = _S1568;
    s_bwd_prop_max_0(&_S1569, &_S1570, v_rgb_2);
    float3  _S1571 = _S1565 * _S1569.differential_0;
    float3  _S1572 = (*sh_coeffs_8)[int(15)] * _S1569.differential_0;
    float3  _S1573 = _S1563 * _S1569.differential_0;
    float3  _S1574 = (*sh_coeffs_8)[int(14)] * _S1569.differential_0;
    float3  _S1575 = _S1561 * _S1569.differential_0;
    float3  _S1576 = (*sh_coeffs_8)[int(13)] * _S1569.differential_0;
    float3  _S1577 = _S1560 * _S1569.differential_0;
    float3  _S1578 = (*sh_coeffs_8)[int(12)] * _S1569.differential_0;
    float3  _S1579 = _S1562 * _S1569.differential_0;
    float3  _S1580 = (*sh_coeffs_8)[int(11)] * _S1569.differential_0;
    float3  _S1581 = _S1564 * _S1569.differential_0;
    float3  _S1582 = (*sh_coeffs_8)[int(10)] * _S1569.differential_0;
    float3  _S1583 = _S1566 * _S1569.differential_0;
    float3  _S1584 = (*sh_coeffs_8)[int(9)] * _S1569.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S1584.x + _S1584.y + _S1584.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S1572.x + _S1572.y + _S1572.z);
    float _S1585 = _S1582.x + _S1582.y + _S1582.z;
    float _S1586 = _S1574.x + _S1574.y + _S1574.z;
    float _S1587 = _S1580.x + _S1580.y + _S1580.z;
    float _S1588 = _S1576.x + _S1576.y + _S1576.z;
    float _S1589 = _S1578.x + _S1578.y + _S1578.z;
    float _S1590 = - s_diff_fC2_T_2;
    float3  _S1591 = _S1557 * _S1569.differential_0;
    float3  _S1592 = (*sh_coeffs_8)[int(8)] * _S1569.differential_0;
    float3  _S1593 = _S1555 * _S1569.differential_0;
    float3  _S1594 = (*sh_coeffs_8)[int(7)] * _S1569.differential_0;
    float3  _S1595 = _S1554 * _S1569.differential_0;
    float3  _S1596 = (*sh_coeffs_8)[int(6)] * _S1569.differential_0;
    float3  _S1597 = _S1556 * _S1569.differential_0;
    float3  _S1598 = (*sh_coeffs_8)[int(5)] * _S1569.differential_0;
    float3  _S1599 = _S1558 * _S1569.differential_0;
    float3  _S1600 = (*sh_coeffs_8)[int(4)] * _S1569.differential_0;
    float _S1601 = _S1598.x + _S1598.y + _S1598.z;
    float _S1602 = _S1594.x + _S1594.y + _S1594.z;
    float _S1603 = fTmp1B_8 * _S1585 + x_38 * s_diff_fS2_T_2 + y_23 * _S1590 + 0.54627424478530884f * (_S1600.x + _S1600.y + _S1600.z);
    float _S1604 = fTmp1B_8 * _S1586 + y_23 * s_diff_fS2_T_2 + x_38 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S1592.x + _S1592.y + _S1592.z);
    float _S1605 = y_23 * - _S1604;
    float _S1606 = x_38 * _S1604;
    float _S1607 = z_20 * (1.86588168144226074f * (z_20 * _S1589) + -2.28522896766662598f * (y_23 * _S1587 + x_38 * _S1588) + 0.94617468118667603f * (_S1596.x + _S1596.y + _S1596.z));
    float3  _S1608 = make_float3 (0.48860251903533936f) * _S1569.differential_0;
    float3  _S1609 = - _S1608;
    float3  _S1610 = _S1548 * _S1609;
    float3  _S1611 = (*sh_coeffs_8)[int(3)] * _S1609;
    float3  _S1612 = _S1550 * _S1608;
    float3  _S1613 = (*sh_coeffs_8)[int(2)] * _S1608;
    float3  _S1614 = _S1552 * _S1608;
    float3  _S1615 = (*sh_coeffs_8)[int(1)] * _S1608;
    float _S1616 = (_S1559 * _S1589 + 1.44530570507049561f * (fS1_8 * _S1585 + fC1_8 * _S1586) + -1.09254848957061768f * (y_23 * _S1601 + x_38 * _S1602) + _S1607 + _S1607 + _S1613.x + _S1613.y + _S1613.z) / _S1549;
    float _S1617 = _S1547 * _S1616;
    float _S1618 = (fTmp0C_8 * _S1587 + fC1_8 * s_diff_fS2_T_2 + fS1_8 * _S1590 + fTmp0B_8 * _S1601 + _S1553 * _S1603 + _S1605 + _S1605 + - (_S1615.x + _S1615.y + _S1615.z)) / _S1549;
    float _S1619 = _S1547 * _S1618;
    float _S1620 = (fTmp0C_8 * _S1588 + fS1_8 * s_diff_fS2_T_2 + fC1_8 * s_diff_fC2_T_2 + fTmp0B_8 * _S1602 + 2.0f * (y_23 * _S1603) + _S1606 + _S1606 + _S1611.x + _S1611.y + _S1611.z) / _S1549;
    float _S1621 = _S1547 * _S1620;
    float _S1622 = _S1545 * - _S1616 + _S1544 * - _S1618 + _S1543 * - _S1620;
    DiffPair_float_0 _S1623;
    (&_S1623)->primal_0 = _S1546;
    (&_S1623)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1623, _S1622);
    float _S1624 = _S1545 * _S1623.differential_0;
    float _S1625 = _S1544 * _S1623.differential_0;
    float _S1626 = _S1543 * _S1623.differential_0;
    float3  _S1627 = make_float3 (0.282094806432724f) * _S1569.differential_0;
    float3  _S1628 = make_float3 (_S1621 + _S1626 + _S1626, _S1619 + _S1625 + _S1625, _S1617 + _S1624 + _S1624);
    float3  _S1629 = - - _S1628;
    Matrix<float, 3, 3>  _S1630 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1631;
    (&_S1631)->primal_0 = _S1515;
    (&_S1631)->differential_0 = _S1630;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1632;
    (&_S1632)->primal_0 = t_9;
    (&_S1632)->differential_0 = _S1568;
    s_bwd_prop_mul_0(&_S1631, &_S1632, _S1629);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1633 = _S1631;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1634 = _S1632;
    float2  _S1635 = make_float2 (0.0f);
    float2  _S1636 = _S1635;
    *&((&_S1636)->y) = v_conic_2.z;
    float2  _S1637 = _S1635;
    *&((&_S1637)->y) = v_conic_2.y;
    *&((&_S1637)->x) = v_conic_2.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1638;
    (&_S1638)->primal_0 = mean_c_8;
    (&_S1638)->differential_0 = _S1568;
    s_bwd_length_impl_1(&_S1638, v_depth_2);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1639 = _S1638;
    DiffPair_float_0 _S1640;
    (&_S1640)->primal_0 = _S1541;
    (&_S1640)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1640, 0.0f);
    DiffPair_float_0 _S1641;
    (&_S1641)->primal_0 = _S1540;
    (&_S1641)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1641, 0.0f);
    DiffPair_float_0 _S1642;
    (&_S1642)->primal_0 = 3.32999992370605469f;
    (&_S1642)->differential_0 = 0.0f;
    DiffPair_float_0 _S1643;
    (&_S1643)->primal_0 = _S1539;
    (&_S1643)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1642, &_S1643, 0.0f);
    DiffPair_float_0 _S1644;
    (&_S1644)->primal_0 = _S1538;
    (&_S1644)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1644, _S1643.differential_0);
    float _S1645 = 2.0f * _S1644.differential_0;
    DiffPair_float_0 _S1646;
    (&_S1646)->primal_0 = _S1537;
    (&_S1646)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1646, _S1645);
    float _S1647 = v_opacity_2 + 254.9999847412109375f * _S1646.differential_0;
    Matrix<float, 2, 2>  _S1648 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1649 = _S1648;
    _S1649[int(1)] = _S1636;
    _S1649[int(0)] = _S1637;
    Matrix<float, 2, 2>  _S1650 = _S1649;
    FixedArray<float3 , 16>  _S1651;
    _S1651[int(0)] = _S1568;
    _S1651[int(1)] = _S1568;
    _S1651[int(2)] = _S1568;
    _S1651[int(3)] = _S1568;
    _S1651[int(4)] = _S1568;
    _S1651[int(5)] = _S1568;
    _S1651[int(6)] = _S1568;
    _S1651[int(7)] = _S1568;
    _S1651[int(8)] = _S1568;
    _S1651[int(9)] = _S1568;
    _S1651[int(10)] = _S1568;
    _S1651[int(11)] = _S1568;
    _S1651[int(12)] = _S1568;
    _S1651[int(13)] = _S1568;
    _S1651[int(14)] = _S1568;
    _S1651[int(15)] = _S1568;
    _S1651[int(7)] = _S1593;
    _S1651[int(0)] = _S1627;
    _S1651[int(1)] = _S1614;
    _S1651[int(2)] = _S1612;
    _S1651[int(3)] = _S1610;
    _S1651[int(4)] = _S1599;
    _S1651[int(5)] = _S1597;
    _S1651[int(6)] = _S1595;
    _S1651[int(15)] = _S1571;
    _S1651[int(8)] = _S1591;
    _S1651[int(9)] = _S1583;
    _S1651[int(10)] = _S1581;
    _S1651[int(11)] = _S1579;
    _S1651[int(12)] = _S1577;
    _S1651[int(13)] = _S1575;
    _S1651[int(14)] = _S1573;
    float3  _S1652 = _S1651[int(0)];
    float3  _S1653 = _S1651[int(1)];
    float3  _S1654 = _S1651[int(2)];
    float3  _S1655 = _S1651[int(3)];
    float3  _S1656 = _S1651[int(4)];
    float3  _S1657 = _S1651[int(5)];
    float3  _S1658 = _S1651[int(6)];
    float3  _S1659 = _S1651[int(7)];
    float3  _S1660 = _S1651[int(8)];
    float3  _S1661 = _S1651[int(9)];
    float3  _S1662 = _S1651[int(10)];
    float3  _S1663 = _S1651[int(11)];
    float3  _S1664 = _S1651[int(12)];
    float3  _S1665 = _S1651[int(13)];
    float3  _S1666 = _S1651[int(14)];
    float3  _S1667 = _S1651[int(15)];
    float2  _S1668 = make_float2 (0.0f, _S1640.differential_0);
    float2  _S1669 = make_float2 (_S1641.differential_0, 0.0f);
    float _S1670;
    if(antialiased_8)
    {
        float _S1671 = _S1534 * _S1647;
        _S1536 = _S1529 * _S1647;
        _S1670 = _S1671;
    }
    else
    {
        _S1536 = _S1647;
        _S1670 = 0.0f;
    }
    float _S1672 = - (_S1536 / _S1535);
    DiffPair_float_0 _S1673;
    (&_S1673)->primal_0 = _S1532;
    (&_S1673)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1673, _S1672);
    float _S1674 = - _S1673.differential_0;
    float _S1675 = invdet_9 * _S1650.rows[int(1)].y;
    float _S1676 = - (invdet_9 * _S1650.rows[int(1)].x);
    float _S1677 = - (invdet_9 * _S1650.rows[int(0)].y);
    float _S1678 = invdet_9 * _S1650.rows[int(0)].x;
    float _S1679 = - ((_S1521 * _S1650.rows[int(1)].y + _S1531 * _S1650.rows[int(1)].x + _S1530 * _S1650.rows[int(0)].y + _S1523 * _S1650.rows[int(0)].x) / _S1527);
    DiffPair_float_0 _S1680;
    (&_S1680)->primal_0 = _S1528;
    (&_S1680)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1680, _S1670);
    DiffPair_float_0 _S1681;
    (&_S1681)->primal_0 = 0.0f;
    (&_S1681)->differential_0 = 0.0f;
    DiffPair_float_0 _S1682;
    (&_S1682)->primal_0 = _S1526;
    (&_S1682)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1681, &_S1682, _S1680.differential_0);
    float _S1683 = _S1682.differential_0 / _S1527;
    float s_diff_det_orig_T_2 = det_blur_6 * _S1683;
    float _S1684 = _S1679 + det_orig_9 * - _S1683;
    float _S1685 = - _S1684;
    float _S1686 = _S1521 * _S1684;
    float _S1687 = _S1523 * _S1684;
    Matrix<float, 2, 2>  _S1688 = _S1648;
    _S1688[int(1)] = _S1668;
    _S1688[int(0)] = _S1669;
    _S1522 = _S1688;
    *&(((&_S1522)->rows + (int(1)))->y) = 0.0f;
    float _S1689 = _S1678 + _S1686 + _S1688.rows[int(1)].y;
    *&(((&_S1522)->rows + (int(0)))->x) = 0.0f;
    float _S1690 = _S1675 + _S1687 + _S1688.rows[int(0)].x;
    float _S1691 = _S1685 + - s_diff_det_orig_T_2;
    float _S1692 = _S1676 + _S1519.rows[int(0)].y * _S1691;
    float _S1693 = _S1677 + _S1519.rows[int(1)].x * _S1691;
    float _S1694 = _S1519.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S1695 = _S1689 + _S1519.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S1696 = _S1635;
    *&((&_S1696)->x) = _S1692;
    *&((&_S1696)->y) = _S1695;
    float _S1697 = _S1690 + _S1694;
    float2  _S1698 = _S1635;
    *&((&_S1698)->y) = _S1693;
    *&((&_S1698)->x) = _S1697;
    float _S1699 = fy_12 * v_mean2d_2.y;
    float _S1700 = fx_12 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S1701 = _S1648;
    _S1701[int(1)] = _S1696;
    _S1701[int(0)] = _S1698;
    Matrix<float, 2, 2>  _S1702 = _S1522 + _S1701;
    Matrix<float, 2, 3>  _S1703 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1704;
    (&_S1704)->primal_0 = _S1517;
    (&_S1704)->differential_0 = _S1703;
    Matrix<float, 3, 2>  _S1705 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1706;
    (&_S1706)->primal_0 = _S1518;
    (&_S1706)->differential_0 = _S1705;
    s_bwd_prop_mul_1(&_S1704, &_S1706, _S1702);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1707;
    (&_S1707)->primal_0 = J_11;
    (&_S1707)->differential_0 = _S1703;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1708;
    (&_S1708)->primal_0 = _S1516;
    (&_S1708)->differential_0 = _S1630;
    s_bwd_prop_mul_2(&_S1707, &_S1708, _S1704.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1709;
    (&_S1709)->primal_0 = _S1514;
    (&_S1709)->differential_0 = _S1630;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1710;
    (&_S1710)->primal_0 = _S1515;
    (&_S1710)->differential_0 = _S1630;
    s_bwd_prop_mul_3(&_S1709, &_S1710, _S1708.differential_0);
    Matrix<float, 3, 3>  _S1711 = transpose_0(_S1710.differential_0 + _S1633.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1712;
    (&_S1712)->primal_0 = R_10;
    (&_S1712)->differential_0 = _S1630;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1713;
    (&_S1713)->primal_0 = _S1513;
    (&_S1713)->differential_0 = _S1630;
    s_bwd_prop_mul_3(&_S1712, &_S1713, _S1709.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1714;
    (&_S1714)->primal_0 = _S1511;
    (&_S1714)->differential_0 = _S1630;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1715;
    (&_S1715)->primal_0 = _S1512;
    (&_S1715)->differential_0 = _S1630;
    s_bwd_prop_mul_3(&_S1714, &_S1715, _S1713.differential_0);
    Matrix<float, 3, 3>  _S1716 = _S1714.differential_0 + transpose_0(_S1715.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1717;
    (&_S1717)->primal_0 = _S1510;
    (&_S1717)->differential_0 = _S1630;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1718;
    (&_S1718)->primal_0 = S_2;
    (&_S1718)->differential_0 = _S1630;
    s_bwd_prop_mul_3(&_S1717, &_S1718, _S1716);
    Matrix<float, 3, 3>  _S1719 = transpose_0(_S1717.differential_0);
    float _S1720 = 2.0f * - _S1719.rows[int(2)].z;
    float _S1721 = 2.0f * _S1719.rows[int(2)].y;
    float _S1722 = 2.0f * _S1719.rows[int(2)].x;
    float _S1723 = 2.0f * _S1719.rows[int(1)].z;
    float _S1724 = 2.0f * - _S1719.rows[int(1)].y;
    float _S1725 = 2.0f * _S1719.rows[int(1)].x;
    float _S1726 = 2.0f * _S1719.rows[int(0)].z;
    float _S1727 = 2.0f * _S1719.rows[int(0)].y;
    float _S1728 = 2.0f * - _S1719.rows[int(0)].x;
    float _S1729 = - _S1725 + _S1727;
    float _S1730 = _S1722 + - _S1726;
    float _S1731 = - _S1721 + _S1723;
    float _S1732 = _S1721 + _S1723;
    float _S1733 = _S1722 + _S1726;
    float _S1734 = _S1725 + _S1727;
    float _S1735 = z_19 * (_S1724 + _S1728);
    float _S1736 = y_22 * (_S1720 + _S1728);
    float _S1737 = x_37 * (_S1720 + _S1724);
    float _S1738 = z_19 * _S1729 + y_22 * _S1730 + x_37 * _S1731;
    float _S1739 = _S1509 * _S1738;
    float _S1740 = w_11 * _S1729 + y_22 * _S1732 + x_37 * _S1733 + _S1735 + _S1735;
    float _S1741 = _S1509 * _S1740;
    float _S1742 = w_11 * _S1730 + z_19 * _S1732 + x_37 * _S1734 + _S1736 + _S1736;
    float _S1743 = _S1509 * _S1742;
    float _S1744 = w_11 * _S1731 + z_19 * _S1733 + y_22 * _S1734 + _S1737 + _S1737;
    float _S1745 = _S1509 * _S1744;
    float _S1746 = quat_11.x * _S1738 + quat_11.w * _S1740 + quat_11.z * _S1742 + quat_11.y * _S1744;
    DiffPair_float_0 _S1747;
    (&_S1747)->primal_0 = _S1508;
    (&_S1747)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S1747, _S1746);
    float _S1748 = quat_11.x * _S1747.differential_0;
    float _S1749 = quat_11.w * _S1747.differential_0;
    float _S1750 = quat_11.z * _S1747.differential_0;
    float _S1751 = quat_11.y * _S1747.differential_0;
    float _S1752 = _S1741 + _S1749 + _S1749;
    float _S1753 = _S1743 + _S1750 + _S1750;
    float _S1754 = _S1745 + _S1751 + _S1751;
    float _S1755 = _S1739 + _S1748 + _S1748;
    float3  _S1756 = _S1568;
    *&((&_S1756)->z) = _S1718.differential_0.rows[int(2)].z;
    *&((&_S1756)->y) = _S1718.differential_0.rows[int(1)].y;
    *&((&_S1756)->x) = _S1718.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1757;
    (&_S1757)->primal_0 = scale_10;
    (&_S1757)->differential_0 = _S1568;
    s_bwd_prop_exp_1(&_S1757, _S1756);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1758 = _S1757;
    float3  _S1759 = _S1568;
    *&((&_S1759)->y) = _S1699;
    *&((&_S1759)->x) = _S1700;
    float3  _S1760 = _S1639.differential_0 + _S1759;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1761;
    (&_S1761)->primal_0 = R_10;
    (&_S1761)->differential_0 = _S1630;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1762;
    (&_S1762)->primal_0 = mean_8;
    (&_S1762)->differential_0 = _S1568;
    s_bwd_prop_mul_0(&_S1761, &_S1762, _S1760);
    float3  _S1763 = _S1760 + _S1634.differential_0;
    Matrix<float, 3, 3>  _S1764 = _S1711 + _S1712.differential_0 + _S1761.differential_0;
    float4  _S1765 = make_float4 (0.0f);
    *&((&_S1765)->w) = _S1752;
    *&((&_S1765)->z) = _S1753;
    *&((&_S1765)->y) = _S1754;
    *&((&_S1765)->x) = _S1755;
    float4  _S1766 = _S1765;
    float3  _S1767 = _S1762.differential_0 + _S1628;
    *v_mean_2 = _S1767;
    *v_quat_2 = _S1766;
    *v_scale_2 = _S1758.differential_0;
    *v_in_opacity_2 = _S1674;
    (*v_sh_coeffs_2)[int(0)] = _S1652;
    (*v_sh_coeffs_2)[int(1)] = _S1653;
    (*v_sh_coeffs_2)[int(2)] = _S1654;
    (*v_sh_coeffs_2)[int(3)] = _S1655;
    (*v_sh_coeffs_2)[int(4)] = _S1656;
    (*v_sh_coeffs_2)[int(5)] = _S1657;
    (*v_sh_coeffs_2)[int(6)] = _S1658;
    (*v_sh_coeffs_2)[int(7)] = _S1659;
    (*v_sh_coeffs_2)[int(8)] = _S1660;
    (*v_sh_coeffs_2)[int(9)] = _S1661;
    (*v_sh_coeffs_2)[int(10)] = _S1662;
    (*v_sh_coeffs_2)[int(11)] = _S1663;
    (*v_sh_coeffs_2)[int(12)] = _S1664;
    (*v_sh_coeffs_2)[int(13)] = _S1665;
    (*v_sh_coeffs_2)[int(14)] = _S1666;
    (*v_sh_coeffs_2)[int(15)] = _S1667;
    *v_R_2 = _S1764;
    *v_t_2 = _S1763;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_13, float dOut_15)
{
    float _S1768 = _slang_select(((*dpx_13).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_13).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_15;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S1768;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_14, float dOut_16)
{
    float _S1769 = (F32_exp2(((*dpx_14).primal_0))) * 50.693145751953125f * dOut_16;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S1769;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_17)
{
    float _S1770 = dOut_17.y;
    float _S1771 = dOut_17.z;
    float _S1772 = dOut_17.x;
    float _S1773 = (*a_0).primal_0.z * _S1770 + - (*a_0).primal_0.y * _S1771;
    float _S1774 = - (*a_0).primal_0.z * _S1772 + (*a_0).primal_0.x * _S1771;
    float _S1775 = (*a_0).primal_0.y * _S1772 + - (*a_0).primal_0.x * _S1770;
    float3  _S1776 = make_float3 (- (*b_0).primal_0.z * _S1770 + (*b_0).primal_0.y * _S1771, (*b_0).primal_0.z * _S1772 + - (*b_0).primal_0.x * _S1771, - (*b_0).primal_0.y * _S1772 + (*b_0).primal_0.x * _S1770);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S1776;
    float3  _S1777 = make_float3 (_S1773, _S1774, _S1775);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S1777;
    return;
}

inline __device__ float3  cross_0(float3  left_8, float3  right_8)
{
    float _S1778 = left_8.y;
    float _S1779 = right_8.z;
    float _S1780 = left_8.z;
    float _S1781 = right_8.y;
    float _S1782 = right_8.x;
    float _S1783 = left_8.x;
    return make_float3 (_S1778 * _S1779 - _S1780 * _S1781, _S1780 * _S1782 - _S1783 * _S1779, _S1783 * _S1781 - _S1778 * _S1782);
}

inline __device__ float3  normalize_0(float3  x_39)
{
    return x_39 / make_float3 (length_1(x_39));
}

inline __device__ float2  normalize_1(float2  x_40)
{
    return x_40 / make_float2 (length_0(x_40));
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_9, float4  quat_12, float3  scale_11, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_9, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_13, float fy_13, float cx_13, float cy_13, float4  radial_coeffs_12, float2  tangential_coeffs_12, float2  thin_prism_coeffs_12, uint image_width_9, uint image_height_9, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_6, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_6, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_9 = mul_0(R_11, mean_9) + t_10;
        float _S1784 = mean_c_9.z;
        bool _S1785;
        if(_S1784 < near_plane_6)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = _S1784 > far_plane_6;
        }
        if(_S1785)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S1786 = scale_11.x;
        float sx_0 = (F32_exp((_S1786)));
        float _S1787 = scale_11.y;
        float sy_0 = (F32_exp((_S1787)));
        float sz_0 = scale_11.z - 0.5f * (_S1786 + _S1787);
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
        Matrix<float, 3, 3>  _S1788 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12)));
        float3  vert0_c_0 = mul_0(R_11, mul_0(_S1788, make_float3 (sx_0, 0.0f, 0.0f)) + mean_9) + t_10;
        float3  vert1_c_0 = mul_0(R_11, mul_0(_S1788, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_9) + t_10;
        float3  vert2_c_0 = mul_0(R_11, mul_0(_S1788, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_9) + t_10;
        float _S1789 = vert0_c_0.z;
        float _S1790 = vert1_c_0.z;
        float _S1791 = vert2_c_0.z;
        if(_S1789 < near_plane_6)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = _S1789 > far_plane_6;
        }
        if(_S1785)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = _S1790 < near_plane_6;
        }
        if(_S1785)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = _S1790 > far_plane_6;
        }
        if(_S1785)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = _S1791 < near_plane_6;
        }
        if(_S1785)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = _S1791 > far_plane_6;
        }
        if(_S1785)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S1789);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S1790);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S1791);
        float2  _S1792 = make_float2 (fx_13, fy_13);
        float2  _S1793 = make_float2 (cx_13, cy_13);
        *uv0_0 = _S1792 * *uv0_0 + _S1793;
        *uv1_0 = _S1792 * *uv1_0 + _S1793;
        float2  _S1794 = _S1792 * *uv2_0 + _S1793;
        *uv2_0 = _S1794;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S1794 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S1794)));
        float _S1795 = _S1794.x;
        float xmax_3 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S1795))) + offset_0;
        float xmin_3 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S1795))) - offset_0;
        float _S1796 = _S1794.y;
        float ymax_3 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S1796))) + offset_0;
        float ymin_3 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S1796))) - offset_0;
        if(xmax_3 <= 0.0f)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = xmin_3 >= float(image_width_9);
        }
        if(_S1785)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = ymax_3 <= 0.0f;
        }
        if(_S1785)
        {
            _S1785 = true;
        }
        else
        {
            _S1785 = ymin_3 >= float(image_height_9);
        }
        if(_S1785)
        {
            _S1785 = true;
        }
        else
        {
            if(_S1784 <= 0.0f)
            {
                if(xmin_3 <= 0.0f)
                {
                    _S1785 = xmax_3 >= float(image_width_9);
                }
                else
                {
                    _S1785 = false;
                }
                if(_S1785)
                {
                    _S1785 = true;
                }
                else
                {
                    if(ymin_3 <= 0.0f)
                    {
                        _S1785 = ymax_3 >= float(image_width_9);
                    }
                    else
                    {
                        _S1785 = false;
                    }
                }
            }
            else
            {
                _S1785 = false;
            }
        }
        if(_S1785)
        {
            *aabb_xyxy_6 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_6 = make_int4 (int((F32_floor((xmin_3)))), int((F32_floor((ymin_3)))), int((F32_ceil((xmax_3)))), int((F32_ceil((ymax_3)))));
        *depth_6 = make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0));
        *out_hardness_0 = hardness_0;
        float3  _S1797 = mean_9 - - mul_0(transpose_0(R_11), t_10);
        float _S1798 = _S1797.x;
        float _S1799 = _S1797.y;
        float _S1800 = _S1797.z;
        float norm_6 = (F32_sqrt((_S1798 * _S1798 + _S1799 * _S1799 + _S1800 * _S1800)));
        float x_43 = _S1798 / norm_6;
        float y_25 = _S1799 / norm_6;
        float z_22 = _S1800 / norm_6;
        float z2_22 = z_22 * z_22;
        float fTmp0B_9 = -1.09254848957061768f * z_22;
        float fC1_9 = x_43 * x_43 - y_25 * y_25;
        float fS1_9 = 2.0f * x_43 * y_25;
        float fTmp0C_9 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
        float fTmp1B_9 = 1.44530570507049561f * z_22;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_9)[int(1)] + make_float3 (z_22) * (*sh_coeffs_9)[int(2)] - make_float3 (x_43) * (*sh_coeffs_9)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_25) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_22 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_43) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_43 * fS1_9 + y_25 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_25) * (*sh_coeffs_9)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_22 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_43) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_43 * fC1_9 - y_25 * fS1_9)) * (*sh_coeffs_9)[int(15)]);
        float3  _S1801 = make_float3 (0.0f);
        (*rgb_6)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S1801);
        float3  _S1802 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S1803 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_6)[int(1)] = max_0(_S1802 + _S1803 + make_float3 (0.5f), _S1801);
        (*rgb_6)[int(2)] = max_0(_S1802 - _S1803 + make_float3 (0.5f), _S1801);
        float3  _S1804 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S1804 * make_float3 (float(- (F32_sign((dot_0(_S1804, mean_c_9))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_10, float4  quat_13, float3  scale_12, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_10, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_14, float fy_14, float cx_14, float cy_14, float4  radial_coeffs_13, float2  tangential_coeffs_13, float2  thin_prism_coeffs_13, uint image_width_10, uint image_height_10, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_7, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_7, float3  * normal_1)
{
    for(;;)
    {
        float3  mean_c_10 = mul_0(R_12, mean_10) + t_11;
        float _S1805 = length_1(mean_c_10);
        bool _S1806;
        if(_S1805 < near_plane_7)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = _S1805 > far_plane_7;
        }
        if(_S1806)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S1807 = scale_12.x;
        float sx_1 = (F32_exp((_S1807)));
        float _S1808 = scale_12.y;
        float sy_1 = (F32_exp((_S1808)));
        float sz_1 = scale_12.z - 0.5f * (_S1807 + _S1808);
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
        Matrix<float, 3, 3>  _S1809 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
        float3  vert0_c_1 = mul_0(R_12, mul_0(_S1809, make_float3 (sx_1, 0.0f, 0.0f)) + mean_10) + t_11;
        float3  vert1_c_1 = mul_0(R_12, mul_0(_S1809, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_10) + t_11;
        float3  vert2_c_1 = mul_0(R_12, mul_0(_S1809, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_10) + t_11;
        float _S1810 = length_1(vert0_c_1);
        float _S1811 = length_1(vert1_c_1);
        float _S1812 = length_1(vert2_c_1);
        if(_S1810 < near_plane_7)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = _S1810 > far_plane_7;
        }
        if(_S1806)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = _S1811 < near_plane_7;
        }
        if(_S1806)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = _S1811 > far_plane_7;
        }
        if(_S1806)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = _S1812 < near_plane_7;
        }
        if(_S1806)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = _S1812 > far_plane_7;
        }
        if(_S1806)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_1 = CameraDistortion_x24init_0(radial_coeffs_13, tangential_coeffs_13, thin_prism_coeffs_13);
        float2  _S1813 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_7 = length_0(_S1813);
        float _S1814 = vert0_c_1.z;
        float theta_1 = (F32_atan2((r_7), (_S1814)));
        float k_4;
        if(theta_1 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_1 * theta_1 / 3.0f) / _S1814;
        }
        else
        {
            k_4 = theta_1 / r_7;
        }
        float2  _S1815 = _S1813 * make_float2 (k_4);
        float k1_3 = dist_coeffs_1.radial_coeffs_0.x;
        float k2_3 = dist_coeffs_1.radial_coeffs_0.y;
        float k3_3 = dist_coeffs_1.radial_coeffs_0.z;
        float k4_3 = dist_coeffs_1.radial_coeffs_0.w;
        float p1_3 = dist_coeffs_1.tangential_coeffs_0.x;
        float p2_3 = dist_coeffs_1.tangential_coeffs_0.y;
        float sx1_3 = dist_coeffs_1.thin_prism_coeffs_0.x;
        float sy1_3 = dist_coeffs_1.thin_prism_coeffs_0.y;
        float u_9 = _S1815.x;
        float v_9 = _S1815.y;
        float r2_9 = u_9 * u_9 + v_9 * v_9;
        float _S1816 = 2.0f * p1_3;
        float _S1817 = 2.0f * p2_3;
        float2  _S1818 = _S1815 * make_float2 (1.0f + r2_9 * (k1_3 + r2_9 * (k2_3 + r2_9 * (k3_3 + r2_9 * k4_3)))) + make_float2 (_S1816 * u_9 * v_9 + p2_3 * (r2_9 + 2.0f * u_9 * u_9) + sx1_3 * r2_9, _S1817 * u_9 * v_9 + p1_3 * (r2_9 + 2.0f * v_9 * v_9) + sy1_3 * r2_9);
        *uv0_1 = make_float2 (fx_14 * _S1818.x + cx_14, fy_14 * _S1818.y + cy_14);
        float2  _S1819 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_8 = length_0(_S1819);
        float _S1820 = vert1_c_1.z;
        float theta_2 = (F32_atan2((r_8), (_S1820)));
        if(theta_2 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_2 * theta_2 / 3.0f) / _S1820;
        }
        else
        {
            k_4 = theta_2 / r_8;
        }
        float2  _S1821 = _S1819 * make_float2 (k_4);
        float u_10 = _S1821.x;
        float v_10 = _S1821.y;
        float r2_10 = u_10 * u_10 + v_10 * v_10;
        float2  _S1822 = _S1821 * make_float2 (1.0f + r2_10 * (k1_3 + r2_10 * (k2_3 + r2_10 * (k3_3 + r2_10 * k4_3)))) + make_float2 (_S1816 * u_10 * v_10 + p2_3 * (r2_10 + 2.0f * u_10 * u_10) + sx1_3 * r2_10, _S1817 * u_10 * v_10 + p1_3 * (r2_10 + 2.0f * v_10 * v_10) + sy1_3 * r2_10);
        *uv1_1 = make_float2 (fx_14 * _S1822.x + cx_14, fy_14 * _S1822.y + cy_14);
        float2  _S1823 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_9 = length_0(_S1823);
        float _S1824 = vert2_c_1.z;
        float theta_3 = (F32_atan2((r_9), (_S1824)));
        if(theta_3 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_3 * theta_3 / 3.0f) / _S1824;
        }
        else
        {
            k_4 = theta_3 / r_9;
        }
        float2  _S1825 = _S1823 * make_float2 (k_4);
        float u_11 = _S1825.x;
        float v_11 = _S1825.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float2  _S1826 = _S1825 * make_float2 (1.0f + r2_11 * (k1_3 + r2_11 * (k2_3 + r2_11 * (k3_3 + r2_11 * k4_3)))) + make_float2 (_S1816 * u_11 * v_11 + p2_3 * (r2_11 + 2.0f * u_11 * u_11) + sx1_3 * r2_11, _S1817 * u_11 * v_11 + p1_3 * (r2_11 + 2.0f * v_11 * v_11) + sy1_3 * r2_11);
        float _S1827 = fx_14 * _S1826.x + cx_14;
        float _S1828 = fy_14 * _S1826.y + cy_14;
        float2  _S1829 = make_float2 (_S1827, _S1828);
        *uv2_1 = _S1829;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S1829 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S1829)));
        float xmax_4 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S1827))) + offset_1;
        float xmin_4 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S1827))) - offset_1;
        float ymax_4 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S1828))) + offset_1;
        float ymin_4 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S1828))) - offset_1;
        if(xmax_4 <= 0.0f)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = xmin_4 >= float(image_width_10);
        }
        if(_S1806)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = ymax_4 <= 0.0f;
        }
        if(_S1806)
        {
            _S1806 = true;
        }
        else
        {
            _S1806 = ymin_4 >= float(image_height_10);
        }
        if(_S1806)
        {
            _S1806 = true;
        }
        else
        {
            if((mean_c_10.z) <= 0.0f)
            {
                if(xmin_4 <= 0.0f)
                {
                    _S1806 = xmax_4 >= float(image_width_10);
                }
                else
                {
                    _S1806 = false;
                }
                if(_S1806)
                {
                    _S1806 = true;
                }
                else
                {
                    if(ymin_4 <= 0.0f)
                    {
                        _S1806 = ymax_4 >= float(image_width_10);
                    }
                    else
                    {
                        _S1806 = false;
                    }
                }
            }
            else
            {
                _S1806 = false;
            }
        }
        if(_S1806)
        {
            *aabb_xyxy_7 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_7 = make_int4 (int((F32_floor((xmin_4)))), int((F32_floor((ymin_4)))), int((F32_ceil((xmax_4)))), int((F32_ceil((ymax_4)))));
        *depth_7 = make_float3 (_S1810, _S1811, _S1812);
        *out_hardness_1 = hardness_1;
        float3  _S1830 = mean_10 - - mul_0(transpose_0(R_12), t_11);
        float _S1831 = _S1830.x;
        float _S1832 = _S1830.y;
        float _S1833 = _S1830.z;
        float norm_7 = (F32_sqrt((_S1831 * _S1831 + _S1832 * _S1832 + _S1833 * _S1833)));
        float x_46 = _S1831 / norm_7;
        float y_27 = _S1832 / norm_7;
        float z_24 = _S1833 / norm_7;
        float z2_24 = z_24 * z_24;
        float fTmp0B_10 = -1.09254848957061768f * z_24;
        float fC1_10 = x_46 * x_46 - y_27 * y_27;
        float fS1_10 = 2.0f * x_46 * y_27;
        float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
        float fTmp1B_10 = 1.44530570507049561f * z_24;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_10)[int(1)] + make_float3 (z_24) * (*sh_coeffs_10)[int(2)] - make_float3 (x_46) * (*sh_coeffs_10)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_10) * (*sh_coeffs_10)[int(4)] + make_float3 (fTmp0B_10 * y_27) * (*sh_coeffs_10)[int(5)] + make_float3 (0.94617468118667603f * z2_24 - 0.31539157032966614f) * (*sh_coeffs_10)[int(6)] + make_float3 (fTmp0B_10 * x_46) * (*sh_coeffs_10)[int(7)] + make_float3 (0.54627424478530884f * fC1_10) * (*sh_coeffs_10)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_10 + y_27 * fC1_10)) * (*sh_coeffs_10)[int(9)] + make_float3 (fTmp1B_10 * fS1_10) * (*sh_coeffs_10)[int(10)] + make_float3 (fTmp0C_10 * y_27) * (*sh_coeffs_10)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_24 - 1.11952900886535645f)) * (*sh_coeffs_10)[int(12)] + make_float3 (fTmp0C_10 * x_46) * (*sh_coeffs_10)[int(13)] + make_float3 (fTmp1B_10 * fC1_10) * (*sh_coeffs_10)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_10 - y_27 * fS1_10)) * (*sh_coeffs_10)[int(15)]);
        float3  _S1834 = make_float3 (0.0f);
        (*rgb_7)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S1834);
        float3  _S1835 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S1836 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_7)[int(1)] = max_0(_S1835 + _S1836 + make_float3 (0.5f), _S1834);
        (*rgb_7)[int(2)] = max_0(_S1835 - _S1836 + make_float3 (0.5f), _S1834);
        float3  _S1837 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S1837 * make_float3 (float(- (F32_sign((dot_0(_S1837, mean_c_10))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_11, float4  quat_14, float3  scale_13, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_11, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_15, float fy_15, float cx_15, float cy_15, float4  radial_coeffs_14, float2  tangential_coeffs_14, float2  thin_prism_coeffs_14, uint image_width_11, uint image_height_11, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_8, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_8, float3  * normal_2)
{
    float3  mean_c_11 = mul_0(R_13, mean_11) + t_12;
    float _S1838 = scale_13.x;
    float sx_2 = (F32_exp((_S1838)));
    float _S1839 = scale_13.y;
    float sy_2 = (F32_exp((_S1839)));
    float sz_2 = scale_13.z - 0.5f * (_S1838 + _S1839);
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
    Matrix<float, 3, 3>  _S1840 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    float3  vert0_c_2 = mul_0(R_13, mul_0(_S1840, make_float3 (sx_2, 0.0f, 0.0f)) + mean_11) + t_12;
    float3  vert1_c_2 = mul_0(R_13, mul_0(_S1840, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_11) + t_12;
    float3  vert2_c_2 = mul_0(R_13, mul_0(_S1840, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_11) + t_12;
    *uv0_2 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    *uv1_2 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    *uv2_2 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float2  _S1841 = make_float2 (fx_15, fy_15);
    float2  _S1842 = make_float2 (cx_15, cy_15);
    *uv0_2 = _S1841 * *uv0_2 + _S1842;
    *uv1_2 = _S1841 * *uv1_2 + _S1842;
    float2  _S1843 = _S1841 * *uv2_2 + _S1842;
    *uv2_2 = _S1843;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S1843 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S1843)));
    float _S1844 = _S1843.x;
    float _S1845 = _S1843.y;
    *aabb_xyxy_8 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S1844))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S1845))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S1844))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S1845))) + offset_2)))));
    *depth_8 = make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2));
    *out_hardness_2 = hardness_2;
    float3  _S1846 = mean_11 - - mul_0(transpose_0(R_13), t_12);
    float _S1847 = _S1846.x;
    float _S1848 = _S1846.y;
    float _S1849 = _S1846.z;
    float norm_8 = (F32_sqrt((_S1847 * _S1847 + _S1848 * _S1848 + _S1849 * _S1849)));
    float x_49 = _S1847 / norm_8;
    float y_29 = _S1848 / norm_8;
    float z_26 = _S1849 / norm_8;
    float z2_26 = z_26 * z_26;
    float fTmp0B_11 = -1.09254848957061768f * z_26;
    float fC1_11 = x_49 * x_49 - y_29 * y_29;
    float fS1_11 = 2.0f * x_49 * y_29;
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_26;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_29) * (*sh_coeffs_11)[int(1)] + make_float3 (z_26) * (*sh_coeffs_11)[int(2)] - make_float3 (x_49) * (*sh_coeffs_11)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_11) * (*sh_coeffs_11)[int(4)] + make_float3 (fTmp0B_11 * y_29) * (*sh_coeffs_11)[int(5)] + make_float3 (0.94617468118667603f * z2_26 - 0.31539157032966614f) * (*sh_coeffs_11)[int(6)] + make_float3 (fTmp0B_11 * x_49) * (*sh_coeffs_11)[int(7)] + make_float3 (0.54627424478530884f * fC1_11) * (*sh_coeffs_11)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_49 * fS1_11 + y_29 * fC1_11)) * (*sh_coeffs_11)[int(9)] + make_float3 (fTmp1B_11 * fS1_11) * (*sh_coeffs_11)[int(10)] + make_float3 (fTmp0C_11 * y_29) * (*sh_coeffs_11)[int(11)] + make_float3 (z_26 * (1.86588168144226074f * z2_26 - 1.11952900886535645f)) * (*sh_coeffs_11)[int(12)] + make_float3 (fTmp0C_11 * x_49) * (*sh_coeffs_11)[int(13)] + make_float3 (fTmp1B_11 * fC1_11) * (*sh_coeffs_11)[int(14)] + make_float3 (-0.59004360437393188f * (x_49 * fC1_11 - y_29 * fS1_11)) * (*sh_coeffs_11)[int(15)]);
    float3  _S1850 = make_float3 (0.0f);
    (*rgb_8)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S1850);
    float3  _S1851 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S1852 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_8)[int(1)] = max_0(_S1851 + _S1852 + make_float3 (0.5f), _S1850);
    (*rgb_8)[int(2)] = max_0(_S1851 - _S1852 + make_float3 (0.5f), _S1850);
    float3  _S1853 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S1853 * make_float3 (float(- (F32_sign((dot_0(_S1853, mean_c_11))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_12, float4  quat_15, float3  scale_14, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_12, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_16, float fy_16, float cx_16, float cy_16, float4  radial_coeffs_15, float2  tangential_coeffs_15, float2  thin_prism_coeffs_15, uint image_width_12, uint image_height_12, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_9, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_9, float3  * normal_3)
{
    float3  mean_c_12 = mul_0(R_14, mean_12) + t_13;
    float _S1854 = scale_14.x;
    float sx_3 = (F32_exp((_S1854)));
    float _S1855 = scale_14.y;
    float sy_3 = (F32_exp((_S1855)));
    float sz_3 = scale_14.z - 0.5f * (_S1854 + _S1855);
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
    Matrix<float, 3, 3>  _S1856 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    float3  vert0_c_3 = mul_0(R_14, mul_0(_S1856, make_float3 (sx_3, 0.0f, 0.0f)) + mean_12) + t_13;
    float3  vert1_c_3 = mul_0(R_14, mul_0(_S1856, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_12) + t_13;
    float3  vert2_c_3 = mul_0(R_14, mul_0(_S1856, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_12) + t_13;
    CameraDistortion_0 dist_coeffs_2 = CameraDistortion_x24init_0(radial_coeffs_15, tangential_coeffs_15, thin_prism_coeffs_15);
    float2  _S1857 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_10 = length_0(_S1857);
    float _S1858 = vert0_c_3.z;
    float theta_4 = (F32_atan2((r_10), (_S1858)));
    float k_5;
    if(theta_4 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_4 * theta_4 / 3.0f) / _S1858;
    }
    else
    {
        k_5 = theta_4 / r_10;
    }
    float2  _S1859 = _S1857 * make_float2 (k_5);
    float k1_4 = dist_coeffs_2.radial_coeffs_0.x;
    float k2_4 = dist_coeffs_2.radial_coeffs_0.y;
    float k3_4 = dist_coeffs_2.radial_coeffs_0.z;
    float k4_4 = dist_coeffs_2.radial_coeffs_0.w;
    float p1_4 = dist_coeffs_2.tangential_coeffs_0.x;
    float p2_4 = dist_coeffs_2.tangential_coeffs_0.y;
    float sx1_4 = dist_coeffs_2.thin_prism_coeffs_0.x;
    float sy1_4 = dist_coeffs_2.thin_prism_coeffs_0.y;
    float u_12 = _S1859.x;
    float v_12 = _S1859.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float _S1860 = 2.0f * p1_4;
    float _S1861 = 2.0f * p2_4;
    float2  _S1862 = _S1859 * make_float2 (1.0f + r2_12 * (k1_4 + r2_12 * (k2_4 + r2_12 * (k3_4 + r2_12 * k4_4)))) + make_float2 (_S1860 * u_12 * v_12 + p2_4 * (r2_12 + 2.0f * u_12 * u_12) + sx1_4 * r2_12, _S1861 * u_12 * v_12 + p1_4 * (r2_12 + 2.0f * v_12 * v_12) + sy1_4 * r2_12);
    *uv0_3 = make_float2 (fx_16 * _S1862.x + cx_16, fy_16 * _S1862.y + cy_16);
    float2  _S1863 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_11 = length_0(_S1863);
    float _S1864 = vert1_c_3.z;
    float theta_5 = (F32_atan2((r_11), (_S1864)));
    if(theta_5 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_5 * theta_5 / 3.0f) / _S1864;
    }
    else
    {
        k_5 = theta_5 / r_11;
    }
    float2  _S1865 = _S1863 * make_float2 (k_5);
    float u_13 = _S1865.x;
    float v_13 = _S1865.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S1866 = _S1865 * make_float2 (1.0f + r2_13 * (k1_4 + r2_13 * (k2_4 + r2_13 * (k3_4 + r2_13 * k4_4)))) + make_float2 (_S1860 * u_13 * v_13 + p2_4 * (r2_13 + 2.0f * u_13 * u_13) + sx1_4 * r2_13, _S1861 * u_13 * v_13 + p1_4 * (r2_13 + 2.0f * v_13 * v_13) + sy1_4 * r2_13);
    *uv1_3 = make_float2 (fx_16 * _S1866.x + cx_16, fy_16 * _S1866.y + cy_16);
    float2  _S1867 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_12 = length_0(_S1867);
    float _S1868 = vert2_c_3.z;
    float theta_6 = (F32_atan2((r_12), (_S1868)));
    if(theta_6 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_6 * theta_6 / 3.0f) / _S1868;
    }
    else
    {
        k_5 = theta_6 / r_12;
    }
    float2  _S1869 = _S1867 * make_float2 (k_5);
    float u_14 = _S1869.x;
    float v_14 = _S1869.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S1870 = _S1869 * make_float2 (1.0f + r2_14 * (k1_4 + r2_14 * (k2_4 + r2_14 * (k3_4 + r2_14 * k4_4)))) + make_float2 (_S1860 * u_14 * v_14 + p2_4 * (r2_14 + 2.0f * u_14 * u_14) + sx1_4 * r2_14, _S1861 * u_14 * v_14 + p1_4 * (r2_14 + 2.0f * v_14 * v_14) + sy1_4 * r2_14);
    float _S1871 = fx_16 * _S1870.x + cx_16;
    float _S1872 = fy_16 * _S1870.y + cy_16;
    float2  _S1873 = make_float2 (_S1871, _S1872);
    *uv2_3 = _S1873;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S1873 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S1873)));
    *aabb_xyxy_9 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S1871))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S1872))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S1871))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S1872))) + offset_3)))));
    *depth_9 = make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3));
    *out_hardness_3 = hardness_3;
    float3  _S1874 = mean_12 - - mul_0(transpose_0(R_14), t_13);
    float _S1875 = _S1874.x;
    float _S1876 = _S1874.y;
    float _S1877 = _S1874.z;
    float norm_9 = (F32_sqrt((_S1875 * _S1875 + _S1876 * _S1876 + _S1877 * _S1877)));
    float x_52 = _S1875 / norm_9;
    float y_31 = _S1876 / norm_9;
    float z_28 = _S1877 / norm_9;
    float z2_28 = z_28 * z_28;
    float fTmp0B_12 = -1.09254848957061768f * z_28;
    float fC1_12 = x_52 * x_52 - y_31 * y_31;
    float fS1_12 = 2.0f * x_52 * y_31;
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_28;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * (*sh_coeffs_12)[int(1)] + make_float3 (z_28) * (*sh_coeffs_12)[int(2)] - make_float3 (x_52) * (*sh_coeffs_12)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_12) * (*sh_coeffs_12)[int(4)] + make_float3 (fTmp0B_12 * y_31) * (*sh_coeffs_12)[int(5)] + make_float3 (0.94617468118667603f * z2_28 - 0.31539157032966614f) * (*sh_coeffs_12)[int(6)] + make_float3 (fTmp0B_12 * x_52) * (*sh_coeffs_12)[int(7)] + make_float3 (0.54627424478530884f * fC1_12) * (*sh_coeffs_12)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_52 * fS1_12 + y_31 * fC1_12)) * (*sh_coeffs_12)[int(9)] + make_float3 (fTmp1B_12 * fS1_12) * (*sh_coeffs_12)[int(10)] + make_float3 (fTmp0C_12 * y_31) * (*sh_coeffs_12)[int(11)] + make_float3 (z_28 * (1.86588168144226074f * z2_28 - 1.11952900886535645f)) * (*sh_coeffs_12)[int(12)] + make_float3 (fTmp0C_12 * x_52) * (*sh_coeffs_12)[int(13)] + make_float3 (fTmp1B_12 * fC1_12) * (*sh_coeffs_12)[int(14)] + make_float3 (-0.59004360437393188f * (x_52 * fC1_12 - y_31 * fS1_12)) * (*sh_coeffs_12)[int(15)]);
    float3  _S1878 = make_float3 (0.0f);
    (*rgb_9)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S1878);
    float3  _S1879 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S1880 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_9)[int(1)] = max_0(_S1879 + _S1880 + make_float3 (0.5f), _S1878);
    (*rgb_9)[int(2)] = max_0(_S1879 - _S1880 + make_float3 (0.5f), _S1878);
    float3  _S1881 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S1881 * make_float3 (float(- (F32_sign((dot_0(_S1881, mean_c_12))))));
    return;
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S1882, float3  _S1883)
{
    return cross_0(_S1882, _S1883);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S1884, float3  _S1885)
{
    return dot_0(_S1884, _S1885);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1886, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1887, float _S1888)
{
    _d_dot_0(_S1886, _S1887, _S1888);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float3  _s_dOut_6)
{
    float _S1889 = length_1((*dpx_15).primal_0);
    float3  _S1890 = (*dpx_15).primal_0 * _s_dOut_6;
    float3  _S1891 = make_float3 (1.0f / _S1889) * _s_dOut_6;
    float _S1892 = - ((_S1890.x + _S1890.y + _S1890.z) / (_S1889 * _S1889));
    float3  _S1893 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1894;
    (&_S1894)->primal_0 = (*dpx_15).primal_0;
    (&_S1894)->differential_0 = _S1893;
    s_bwd_length_impl_1(&_S1894, _S1892);
    float3  _S1895 = _S1891 + _S1894.differential_0;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S1895;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1896, float3  _S1897)
{
    s_bwd_prop_normalize_impl_0(_S1896, _S1897);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1898, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1899, float3  _S1900)
{
    _d_cross_0(_S1898, _S1899, _S1900);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S1901, float _S1902)
{
    _d_exp2_0(_S1901, _S1902);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S1903, float _S1904)
{
    _d_abs_0(_S1903, _S1904);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_13, float4  quat_16, float3  scale_15, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_13, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_17, float fy_17, float cx_17, float cy_17, float4  radial_coeffs_16, float2  tangential_coeffs_16, float2  thin_prism_coeffs_16, uint image_width_13, uint image_height_13, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_3, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_3, float3  v_normal_0, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_3, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_13 = s_primal_ctx_mul_0(R_15, mean_13) + t_14;
    float _S1905 = scale_15.x;
    float _S1906 = s_primal_ctx_exp_1(_S1905);
    float _S1907 = scale_15.y;
    float _S1908 = s_primal_ctx_exp_1(_S1907);
    float sz_4 = scale_15.z - 0.5f * (_S1905 + _S1907);
    float _S1909 = quat_16.y;
    float _S1910 = _S1909 * _S1909 + quat_16.z * quat_16.z + quat_16.w * quat_16.w + quat_16.x * quat_16.x;
    float _S1911 = s_primal_ctx_rsqrt_0(_S1910);
    float x_53 = quat_16.y * _S1911;
    float y_32 = quat_16.z * _S1911;
    float z_29 = quat_16.w * _S1911;
    float w_16 = quat_16.x * _S1911;
    float x2_16 = x_53 * x_53;
    float y2_16 = y_32 * y_32;
    float z2_29 = z_29 * z_29;
    float xy_16 = x_53 * y_32;
    float xz_16 = x_53 * z_29;
    float yz_16 = y_32 * z_29;
    float wx_16 = w_16 * x_53;
    float wy_16 = w_16 * y_32;
    float wz_16 = w_16 * z_29;
    Matrix<float, 3, 3>  _S1912 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    float3  _S1913 = make_float3 (_S1906, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_0(_S1912, _S1913) + mean_13;
    float _S1914 = -0.5f + sz_4;
    float3  _S1915 = make_float3 (_S1906 * _S1914, _S1908, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_0(_S1912, _S1915) + mean_13;
    float _S1916 = -0.5f - sz_4;
    float3  _S1917 = make_float3 (_S1906 * _S1916, - _S1908, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_0(_S1912, _S1917) + mean_13;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_15, vert0_0) + t_14;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_15, vert1_0) + t_14;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_15, vert2_0) + t_14;
    float2  _S1918 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S1919 = vert0_c_4.z;
    float2  _S1920 = make_float2 (_S1919);
    float2  _S1921 = make_float2 (_S1919 * _S1919);
    float2  _S1922 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S1923 = vert1_c_4.z;
    float2  _S1924 = make_float2 (_S1923);
    float2  _S1925 = make_float2 (_S1923 * _S1923);
    float2  _S1926 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S1927 = vert2_c_4.z;
    float2  _S1928 = make_float2 (_S1927);
    float2  _S1929 = make_float2 (_S1927 * _S1927);
    float2  _S1930 = make_float2 (fx_17, fy_17);
    float2  _S1931 = make_float2 (cx_17, cy_17);
    float2  _S1932 = _S1930 * (_S1918 / make_float2 (_S1919)) + _S1931;
    float2  _S1933 = _S1930 * (_S1922 / make_float2 (_S1923)) + _S1931;
    float2  _S1934 = _S1930 * (_S1926 / make_float2 (_S1927)) + _S1931;
    float2  e0_4 = _S1933 - _S1932;
    float2  e1_4 = _S1934 - _S1933;
    float2  e2_0 = _S1932 - _S1934;
    float _S1935 = e0_4.x;
    float _S1936 = e1_4.y;
    float _S1937 = e0_4.y;
    float _S1938 = e1_4.x;
    float _S1939 = _S1935 * _S1936 - _S1937 * _S1938;
    float _S1940 = 1.0f - hardness_4.y;
    float _S1941 = -1.0f / _S1940;
    float _S1942 = _S1940 * _S1940;
    float _S1943 = _S1932.x;
    float _S1944 = _S1933.x;
    float _S1945 = s_primal_ctx_max_0(_S1943, _S1944);
    float _S1946 = _S1934.x;
    float _S1947 = s_primal_ctx_min_0(_S1943, _S1944);
    float _S1948 = _S1932.y;
    float _S1949 = _S1933.y;
    float _S1950 = s_primal_ctx_max_0(_S1948, _S1949);
    float _S1951 = _S1934.y;
    float _S1952 = s_primal_ctx_min_0(_S1948, _S1949);
    Matrix<float, 3, 3>  _S1953 = transpose_0(R_15);
    float3  _S1954 = mean_13 - - s_primal_ctx_mul_0(_S1953, t_14);
    float _S1955 = _S1954.x;
    float _S1956 = _S1954.y;
    float _S1957 = _S1954.z;
    float _S1958 = _S1955 * _S1955 + _S1956 * _S1956 + _S1957 * _S1957;
    float _S1959 = s_primal_ctx_sqrt_0(_S1958);
    float x_54 = _S1955 / _S1959;
    float3  _S1960 = make_float3 (x_54);
    float _S1961 = _S1959 * _S1959;
    float y_33 = _S1956 / _S1959;
    float z_30 = _S1957 / _S1959;
    float3  _S1962 = make_float3 (z_30);
    float _S1963 = - y_33;
    float3  _S1964 = make_float3 (_S1963);
    float z2_30 = z_30 * z_30;
    float fTmp0B_13 = -1.09254848957061768f * z_30;
    float fC1_13 = x_54 * x_54 - y_33 * y_33;
    float _S1965 = 2.0f * x_54;
    float fS1_13 = _S1965 * y_33;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S1966 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_54;
    float3  _S1967 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_33;
    float3  _S1968 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S1969 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S1970 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_30;
    float _S1971 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_30 * _S1971;
    float3  _S1972 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_54;
    float3  _S1973 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_33;
    float3  _S1974 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S1975 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S1976 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_54 * fC1_13 - y_33 * fS1_13);
    float3  _S1977 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_54 * fS1_13 + y_33 * fC1_13);
    float3  _S1978 = make_float3 (pSH9_3);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1963) * (*sh_coeffs_13)[int(1)] + make_float3 (z_30) * (*sh_coeffs_13)[int(2)] - make_float3 (x_54) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]);
    float3  _S1979 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S1980 = make_float3 (0.0f);
    float3  _S1981 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S1982 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S1983 = make_float3 (_S1982);
    float3  _S1984 = (*ch_coeffs_4)[int(1)] * make_float3 (_S1982);
    float3  _S1985 = _S1981 + _S1984 + make_float3 (0.5f);
    float3  _S1986 = _S1981 - _S1984 + make_float3 (0.5f);
    float3  _S1987 = vert1_c_4 - vert0_c_4;
    float3  _S1988 = vert2_c_4 - vert0_c_4;
    float3  _S1989 = s_primal_ctx_cross_0(_S1987, _S1988);
    float3  _S1990 = normalize_0(_S1989);
    float3  _S1991 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S1990, mean_c_13)))))) * v_normal_0;
    float3  _S1992 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1993;
    (&_S1993)->primal_0 = _S1990;
    (&_S1993)->differential_0 = _S1992;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1994;
    (&_S1994)->primal_0 = mean_c_13;
    (&_S1994)->differential_0 = _S1992;
    s_bwd_prop_dot_0(&_S1993, &_S1994, 0.0f);
    float3  _S1995 = _S1991 + _S1993.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1996;
    (&_S1996)->primal_0 = _S1989;
    (&_S1996)->differential_0 = _S1992;
    s_bwd_normalize_impl_0(&_S1996, _S1995);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1997;
    (&_S1997)->primal_0 = _S1987;
    (&_S1997)->differential_0 = _S1992;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1998;
    (&_S1998)->primal_0 = _S1988;
    (&_S1998)->differential_0 = _S1992;
    s_bwd_prop_cross_0(&_S1997, &_S1998, _S1996.differential_0);
    float3  _S1999 = - _S1998.differential_0;
    float3  _S2000 = - _S1997.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2001;
    (&_S2001)->primal_0 = _S1986;
    (&_S2001)->differential_0 = _S1992;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2002;
    (&_S2002)->primal_0 = _S1980;
    (&_S2002)->differential_0 = _S1992;
    s_bwd_prop_max_0(&_S2001, &_S2002, (*v_rgb_3)[int(2)]);
    float3  _S2003 = - _S2001.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2004;
    (&_S2004)->primal_0 = _S1985;
    (&_S2004)->differential_0 = _S1992;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2005;
    (&_S2005)->primal_0 = _S1980;
    (&_S2005)->differential_0 = _S1992;
    s_bwd_prop_max_0(&_S2004, &_S2005, (*v_rgb_3)[int(1)]);
    float3  _S2006 = _S1983 * (_S2003 + _S2004.differential_0);
    float3  _S2007 = _S2001.differential_0 + _S2004.differential_0;
    float3  _S2008 = make_float3 (0.5f) * - _S2007;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2009;
    (&_S2009)->primal_0 = _S1979;
    (&_S2009)->differential_0 = _S1992;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2010;
    (&_S2010)->primal_0 = _S1980;
    (&_S2010)->differential_0 = _S1992;
    s_bwd_prop_max_0(&_S2009, &_S2010, (*v_rgb_3)[int(0)]);
    float3  _S2011 = _S2008 + _S2009.differential_0;
    float3  _S2012 = _S2007 + _S2009.differential_0;
    float3  _S2013 = _S1977 * _S2012;
    float3  _S2014 = (*sh_coeffs_13)[int(15)] * _S2012;
    float3  _S2015 = _S1975 * _S2012;
    float3  _S2016 = (*sh_coeffs_13)[int(14)] * _S2012;
    float3  _S2017 = _S1973 * _S2012;
    float3  _S2018 = (*sh_coeffs_13)[int(13)] * _S2012;
    float3  _S2019 = _S1972 * _S2012;
    float3  _S2020 = (*sh_coeffs_13)[int(12)] * _S2012;
    float3  _S2021 = _S1974 * _S2012;
    float3  _S2022 = (*sh_coeffs_13)[int(11)] * _S2012;
    float3  _S2023 = _S1976 * _S2012;
    float3  _S2024 = (*sh_coeffs_13)[int(10)] * _S2012;
    float3  _S2025 = _S1978 * _S2012;
    float3  _S2026 = (*sh_coeffs_13)[int(9)] * _S2012;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2026.x + _S2026.y + _S2026.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S2014.x + _S2014.y + _S2014.z);
    float _S2027 = _S2024.x + _S2024.y + _S2024.z;
    float _S2028 = _S2016.x + _S2016.y + _S2016.z;
    float _S2029 = _S2022.x + _S2022.y + _S2022.z;
    float _S2030 = _S2018.x + _S2018.y + _S2018.z;
    float _S2031 = _S2020.x + _S2020.y + _S2020.z;
    float _S2032 = - s_diff_fC2_T_3;
    float3  _S2033 = _S1969 * _S2012;
    float3  _S2034 = (*sh_coeffs_13)[int(8)] * _S2012;
    float3  _S2035 = _S1967 * _S2012;
    float3  _S2036 = (*sh_coeffs_13)[int(7)] * _S2012;
    float3  _S2037 = _S1966 * _S2012;
    float3  _S2038 = (*sh_coeffs_13)[int(6)] * _S2012;
    float3  _S2039 = _S1968 * _S2012;
    float3  _S2040 = (*sh_coeffs_13)[int(5)] * _S2012;
    float3  _S2041 = _S1970 * _S2012;
    float3  _S2042 = (*sh_coeffs_13)[int(4)] * _S2012;
    float _S2043 = _S2040.x + _S2040.y + _S2040.z;
    float _S2044 = _S2036.x + _S2036.y + _S2036.z;
    float _S2045 = fTmp1B_13 * _S2027 + x_54 * s_diff_fS2_T_3 + y_33 * _S2032 + 0.54627424478530884f * (_S2042.x + _S2042.y + _S2042.z);
    float _S2046 = fTmp1B_13 * _S2028 + y_33 * s_diff_fS2_T_3 + x_54 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2034.x + _S2034.y + _S2034.z);
    float _S2047 = y_33 * - _S2046;
    float _S2048 = x_54 * _S2046;
    float _S2049 = z_30 * (1.86588168144226074f * (z_30 * _S2031) + -2.28522896766662598f * (y_33 * _S2029 + x_54 * _S2030) + 0.94617468118667603f * (_S2038.x + _S2038.y + _S2038.z));
    float3  _S2050 = make_float3 (0.48860251903533936f) * _S2012;
    float3  _S2051 = - _S2050;
    float3  _S2052 = _S1960 * _S2051;
    float3  _S2053 = (*sh_coeffs_13)[int(3)] * _S2051;
    float3  _S2054 = _S1962 * _S2050;
    float3  _S2055 = (*sh_coeffs_13)[int(2)] * _S2050;
    float3  _S2056 = _S1964 * _S2050;
    float3  _S2057 = (*sh_coeffs_13)[int(1)] * _S2050;
    float _S2058 = (_S1971 * _S2031 + 1.44530570507049561f * (fS1_13 * _S2027 + fC1_13 * _S2028) + -1.09254848957061768f * (y_33 * _S2043 + x_54 * _S2044) + _S2049 + _S2049 + _S2055.x + _S2055.y + _S2055.z) / _S1961;
    float _S2059 = _S1959 * _S2058;
    float _S2060 = (fTmp0C_13 * _S2029 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2032 + fTmp0B_13 * _S2043 + _S1965 * _S2045 + _S2047 + _S2047 + - (_S2057.x + _S2057.y + _S2057.z)) / _S1961;
    float _S2061 = _S1959 * _S2060;
    float _S2062 = (fTmp0C_13 * _S2030 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2044 + 2.0f * (y_33 * _S2045) + _S2048 + _S2048 + _S2053.x + _S2053.y + _S2053.z) / _S1961;
    float _S2063 = _S1959 * _S2062;
    float _S2064 = _S1957 * - _S2058 + _S1956 * - _S2060 + _S1955 * - _S2062;
    DiffPair_float_0 _S2065;
    (&_S2065)->primal_0 = _S1958;
    (&_S2065)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2065, _S2064);
    float _S2066 = _S1957 * _S2065.differential_0;
    float _S2067 = _S1956 * _S2065.differential_0;
    float _S2068 = _S1955 * _S2065.differential_0;
    float3  _S2069 = make_float3 (0.282094806432724f) * _S2012;
    float3  _S2070 = make_float3 (_S2063 + _S2068 + _S2068, _S2061 + _S2067 + _S2067, _S2059 + _S2066 + _S2066);
    float3  _S2071 = - - _S2070;
    Matrix<float, 3, 3>  _S2072 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2073;
    (&_S2073)->primal_0 = _S1953;
    (&_S2073)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2074;
    (&_S2074)->primal_0 = t_14;
    (&_S2074)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2073, &_S2074, _S2071);
    Matrix<float, 3, 3>  _S2075 = transpose_0(_S2073.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2076;
    (&_S2076)->primal_0 = vert2_c_4;
    (&_S2076)->differential_0 = _S1992;
    s_bwd_length_impl_1(&_S2076, v_depth_3.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2077;
    (&_S2077)->primal_0 = vert1_c_4;
    (&_S2077)->differential_0 = _S1992;
    s_bwd_length_impl_1(&_S2077, v_depth_3.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2078;
    (&_S2078)->primal_0 = vert0_c_4;
    (&_S2078)->differential_0 = _S1992;
    s_bwd_length_impl_1(&_S2078, v_depth_3.x);
    DiffPair_float_0 _S2079;
    (&_S2079)->primal_0 = _S1952;
    (&_S2079)->differential_0 = 0.0f;
    DiffPair_float_0 _S2080;
    (&_S2080)->primal_0 = _S1951;
    (&_S2080)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2079, &_S2080, 0.0f);
    DiffPair_float_0 _S2081;
    (&_S2081)->primal_0 = _S1948;
    (&_S2081)->differential_0 = 0.0f;
    DiffPair_float_0 _S2082;
    (&_S2082)->primal_0 = _S1949;
    (&_S2082)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2081, &_S2082, _S2079.differential_0);
    DiffPair_float_0 _S2083;
    (&_S2083)->primal_0 = _S1950;
    (&_S2083)->differential_0 = 0.0f;
    DiffPair_float_0 _S2084;
    (&_S2084)->primal_0 = _S1951;
    (&_S2084)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2083, &_S2084, 0.0f);
    float _S2085 = _S2080.differential_0 + _S2084.differential_0;
    DiffPair_float_0 _S2086;
    (&_S2086)->primal_0 = _S1948;
    (&_S2086)->differential_0 = 0.0f;
    DiffPair_float_0 _S2087;
    (&_S2087)->primal_0 = _S1949;
    (&_S2087)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2086, &_S2087, _S2083.differential_0);
    float _S2088 = _S2082.differential_0 + _S2087.differential_0;
    float _S2089 = _S2081.differential_0 + _S2086.differential_0;
    DiffPair_float_0 _S2090;
    (&_S2090)->primal_0 = _S1947;
    (&_S2090)->differential_0 = 0.0f;
    DiffPair_float_0 _S2091;
    (&_S2091)->primal_0 = _S1946;
    (&_S2091)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2090, &_S2091, 0.0f);
    DiffPair_float_0 _S2092;
    (&_S2092)->primal_0 = _S1943;
    (&_S2092)->differential_0 = 0.0f;
    DiffPair_float_0 _S2093;
    (&_S2093)->primal_0 = _S1944;
    (&_S2093)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2092, &_S2093, _S2090.differential_0);
    DiffPair_float_0 _S2094;
    (&_S2094)->primal_0 = _S1945;
    (&_S2094)->differential_0 = 0.0f;
    DiffPair_float_0 _S2095;
    (&_S2095)->primal_0 = _S1946;
    (&_S2095)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2094, &_S2095, 0.0f);
    float _S2096 = _S2091.differential_0 + _S2095.differential_0;
    DiffPair_float_0 _S2097;
    (&_S2097)->primal_0 = _S1943;
    (&_S2097)->differential_0 = 0.0f;
    DiffPair_float_0 _S2098;
    (&_S2098)->primal_0 = _S1944;
    (&_S2098)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2097, &_S2098, _S2094.differential_0);
    float _S2099 = _S2093.differential_0 + _S2098.differential_0;
    float _S2100 = _S2092.differential_0 + _S2097.differential_0;
    DiffPair_float_0 _S2101;
    (&_S2101)->primal_0 = _S1941;
    (&_S2101)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2101, 0.0f);
    float _S2102 = - (-1.0f * - (_S2101.differential_0 / _S1942));
    float2  _S2103 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2104;
    (&_S2104)->primal_0 = e2_0;
    (&_S2104)->differential_0 = _S2103;
    s_bwd_length_impl_0(&_S2104, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2105;
    (&_S2105)->primal_0 = e1_4;
    (&_S2105)->differential_0 = _S2103;
    s_bwd_length_impl_0(&_S2105, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2106;
    (&_S2106)->primal_0 = e0_4;
    (&_S2106)->differential_0 = _S2103;
    s_bwd_length_impl_0(&_S2106, -0.0f);
    DiffPair_float_0 _S2107;
    (&_S2107)->primal_0 = _S1939;
    (&_S2107)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2107, 0.0f);
    float _S2108 = - _S2107.differential_0;
    float2  _S2109 = _S2105.differential_0 + make_float2 (_S1937 * _S2108, _S1935 * _S2107.differential_0);
    float2  _S2110 = _S2106.differential_0 + make_float2 (_S1936 * _S2107.differential_0, _S1938 * _S2108);
    float2  _S2111 = _S1930 * (v_uv2_0 + - _S2104.differential_0 + _S2109 + make_float2 (_S2096, _S2085)) / _S1929;
    float2  _S2112 = _S1926 * - _S2111;
    float2  _S2113 = _S1928 * _S2111;
    float2  _S2114 = _S1930 * (v_uv1_0 + - _S2109 + _S2110 + make_float2 (_S2099, _S2088)) / _S1925;
    float2  _S2115 = _S1922 * - _S2114;
    float2  _S2116 = _S1924 * _S2114;
    float _S2117 = _S2115.x + _S2115.y;
    float2  _S2118 = _S1930 * (v_uv0_0 + _S2104.differential_0 + - _S2110 + make_float2 (_S2100, _S2089)) / _S1921;
    float2  _S2119 = _S1918 * - _S2118;
    float2  _S2120 = _S1920 * _S2118;
    float _S2121 = _S2119.x + _S2119.y;
    float3  _S2122 = _S1998.differential_0 + _S2076.differential_0 + make_float3 (_S2113.x, _S2113.y, _S2112.x + _S2112.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2123;
    (&_S2123)->primal_0 = R_15;
    (&_S2123)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2124;
    (&_S2124)->primal_0 = vert2_0;
    (&_S2124)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2123, &_S2124, _S2122);
    float3  _S2125 = _S1997.differential_0 + _S2077.differential_0 + make_float3 (_S2116.x, _S2116.y, _S2117);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2126;
    (&_S2126)->primal_0 = R_15;
    (&_S2126)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2127;
    (&_S2127)->primal_0 = vert1_0;
    (&_S2127)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2126, &_S2127, _S2125);
    float3  _S2128 = _S1999 + _S2000 + _S2078.differential_0 + make_float3 (_S2120.x, _S2120.y, _S2121);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2129;
    (&_S2129)->primal_0 = R_15;
    (&_S2129)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2130;
    (&_S2130)->primal_0 = vert0_0;
    (&_S2130)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2129, &_S2130, _S2128);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2131;
    (&_S2131)->primal_0 = _S1912;
    (&_S2131)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2132;
    (&_S2132)->primal_0 = _S1917;
    (&_S2132)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2131, &_S2132, _S2124.differential_0);
    float _S2133 = - _S2132.differential_0.y;
    float _S2134 = _S1916 * _S2132.differential_0.x;
    float _S2135 = - (_S1906 * _S2132.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2136;
    (&_S2136)->primal_0 = _S1912;
    (&_S2136)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2137;
    (&_S2137)->primal_0 = _S1915;
    (&_S2137)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2136, &_S2137, _S2127.differential_0);
    float _S2138 = _S1906 * _S2137.differential_0.x;
    float _S2139 = _S1914 * _S2137.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2140;
    (&_S2140)->primal_0 = _S1912;
    (&_S2140)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2141;
    (&_S2141)->primal_0 = _S1913;
    (&_S2141)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2140, &_S2141, _S2130.differential_0);
    Matrix<float, 3, 3>  _S2142 = transpose_0(_S2131.differential_0 + _S2136.differential_0 + _S2140.differential_0);
    float _S2143 = 2.0f * - _S2142.rows[int(2)].z;
    float _S2144 = 2.0f * _S2142.rows[int(2)].y;
    float _S2145 = 2.0f * _S2142.rows[int(2)].x;
    float _S2146 = 2.0f * _S2142.rows[int(1)].z;
    float _S2147 = 2.0f * - _S2142.rows[int(1)].y;
    float _S2148 = 2.0f * _S2142.rows[int(1)].x;
    float _S2149 = 2.0f * _S2142.rows[int(0)].z;
    float _S2150 = 2.0f * _S2142.rows[int(0)].y;
    float _S2151 = 2.0f * - _S2142.rows[int(0)].x;
    float _S2152 = - _S2148 + _S2150;
    float _S2153 = _S2145 + - _S2149;
    float _S2154 = - _S2144 + _S2146;
    float _S2155 = _S2144 + _S2146;
    float _S2156 = _S2145 + _S2149;
    float _S2157 = _S2148 + _S2150;
    float _S2158 = z_29 * (_S2147 + _S2151);
    float _S2159 = y_32 * (_S2143 + _S2151);
    float _S2160 = x_53 * (_S2143 + _S2147);
    float _S2161 = z_29 * _S2152 + y_32 * _S2153 + x_53 * _S2154;
    float _S2162 = _S1911 * _S2161;
    float _S2163 = w_16 * _S2152 + y_32 * _S2155 + x_53 * _S2156 + _S2158 + _S2158;
    float _S2164 = _S1911 * _S2163;
    float _S2165 = w_16 * _S2153 + z_29 * _S2155 + x_53 * _S2157 + _S2159 + _S2159;
    float _S2166 = _S1911 * _S2165;
    float _S2167 = w_16 * _S2154 + z_29 * _S2156 + y_32 * _S2157 + _S2160 + _S2160;
    float _S2168 = _S1911 * _S2167;
    float _S2169 = quat_16.x * _S2161 + quat_16.w * _S2163 + quat_16.z * _S2165 + quat_16.y * _S2167;
    DiffPair_float_0 _S2170;
    (&_S2170)->primal_0 = _S1910;
    (&_S2170)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S2170, _S2169);
    float _S2171 = quat_16.x * _S2170.differential_0;
    float _S2172 = quat_16.w * _S2170.differential_0;
    float _S2173 = quat_16.z * _S2170.differential_0;
    float _S2174 = quat_16.y * _S2170.differential_0;
    float _S2175 = _S2164 + _S2172 + _S2172;
    float _S2176 = _S2166 + _S2173 + _S2173;
    float _S2177 = _S2168 + _S2174 + _S2174;
    float _S2178 = _S2162 + _S2171 + _S2171;
    float _S2179 = _S2135 + _S2138;
    float _S2180 = 0.5f * - _S2179;
    float _S2181 = _S2133 + _S2137.differential_0.y;
    DiffPair_float_0 _S2182;
    (&_S2182)->primal_0 = _S1907;
    (&_S2182)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2182, _S2181);
    float _S2183 = _S2180 + _S2182.differential_0;
    float _S2184 = _S2134 + _S2139 + _S2141.differential_0.x;
    DiffPair_float_0 _S2185;
    (&_S2185)->primal_0 = _S1905;
    (&_S2185)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2185, _S2184);
    float _S2186 = _S2180 + _S2185.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2187;
    (&_S2187)->primal_0 = R_15;
    (&_S2187)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2188;
    (&_S2188)->primal_0 = mean_13;
    (&_S2188)->differential_0 = _S1992;
    s_bwd_prop_mul_0(&_S2187, &_S2188, _S1994.differential_0);
    float3  _S2189 = _S2074.differential_0 + _S2122 + _S2125 + _S2128 + _S1994.differential_0;
    Matrix<float, 3, 3>  _S2190 = _S2075 + _S2123.differential_0 + _S2126.differential_0 + _S2129.differential_0 + _S2187.differential_0;
    FixedArray<float3 , 2>  _S2191;
    _S2191[int(0)] = _S1992;
    _S2191[int(1)] = _S1992;
    _S2191[int(1)] = _S2006;
    _S2191[int(0)] = _S2011;
    FixedArray<float3 , 16>  _S2192;
    _S2192[int(0)] = _S1992;
    _S2192[int(1)] = _S1992;
    _S2192[int(2)] = _S1992;
    _S2192[int(3)] = _S1992;
    _S2192[int(4)] = _S1992;
    _S2192[int(5)] = _S1992;
    _S2192[int(6)] = _S1992;
    _S2192[int(7)] = _S1992;
    _S2192[int(8)] = _S1992;
    _S2192[int(9)] = _S1992;
    _S2192[int(10)] = _S1992;
    _S2192[int(11)] = _S1992;
    _S2192[int(12)] = _S1992;
    _S2192[int(13)] = _S1992;
    _S2192[int(14)] = _S1992;
    _S2192[int(15)] = _S1992;
    _S2192[int(15)] = _S2013;
    _S2192[int(14)] = _S2015;
    _S2192[int(13)] = _S2017;
    _S2192[int(12)] = _S2019;
    _S2192[int(11)] = _S2021;
    _S2192[int(10)] = _S2023;
    _S2192[int(9)] = _S2025;
    _S2192[int(8)] = _S2033;
    _S2192[int(7)] = _S2035;
    _S2192[int(6)] = _S2037;
    _S2192[int(5)] = _S2039;
    _S2192[int(4)] = _S2041;
    _S2192[int(3)] = _S2052;
    _S2192[int(2)] = _S2054;
    _S2192[int(1)] = _S2056;
    _S2192[int(0)] = _S2069;
    float2  _S2193 = v_out_hardness_0 + make_float2 (0.0f, _S2102);
    float3  _S2194 = make_float3 (_S2186, _S2183, _S2179);
    float4  _S2195 = make_float4 (0.0f);
    *&((&_S2195)->w) = _S2175;
    *&((&_S2195)->z) = _S2176;
    *&((&_S2195)->y) = _S2177;
    *&((&_S2195)->x) = _S2178;
    *v_mean_3 = _S2070 + _S2124.differential_0 + _S2127.differential_0 + _S2130.differential_0 + _S2188.differential_0;
    *v_quat_3 = _S2195;
    *v_scale_3 = _S2194;
    *v_hardness_0 = _S2193;
    *v_sh_coeffs_3 = _S2192;
    *v_ch_coeffs_0 = _S2191;
    *v_R_3 = _S2190;
    *v_t_3 = _S2189;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_14, float4  quat_17, float3  scale_16, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_14, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_18, float fy_18, float cx_18, float cy_18, float4  radial_coeffs_17, float2  tangential_coeffs_17, float2  thin_prism_coeffs_17, uint image_width_14, uint image_height_14, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_4, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_4, float3  v_normal_1, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_4, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float3  mean_c_14 = s_primal_ctx_mul_0(R_16, mean_14) + t_15;
    float _S2196 = scale_16.x;
    float _S2197 = s_primal_ctx_exp_1(_S2196);
    float _S2198 = scale_16.y;
    float _S2199 = s_primal_ctx_exp_1(_S2198);
    float sz_5 = scale_16.z - 0.5f * (_S2196 + _S2198);
    float _S2200 = quat_17.y;
    float _S2201 = _S2200 * _S2200 + quat_17.z * quat_17.z + quat_17.w * quat_17.w + quat_17.x * quat_17.x;
    float _S2202 = s_primal_ctx_rsqrt_0(_S2201);
    float x_55 = quat_17.y * _S2202;
    float y_34 = quat_17.z * _S2202;
    float z_31 = quat_17.w * _S2202;
    float w_17 = quat_17.x * _S2202;
    float x2_17 = x_55 * x_55;
    float y2_17 = y_34 * y_34;
    float z2_31 = z_31 * z_31;
    float xy_17 = x_55 * y_34;
    float xz_17 = x_55 * z_31;
    float yz_17 = y_34 * z_31;
    float wx_17 = w_17 * x_55;
    float wy_17 = w_17 * y_34;
    float wz_17 = w_17 * z_31;
    Matrix<float, 3, 3>  _S2203 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    float3  _S2204 = make_float3 (_S2197, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_0(_S2203, _S2204) + mean_14;
    float _S2205 = -0.5f + sz_5;
    float3  _S2206 = make_float3 (_S2197 * _S2205, _S2199, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_0(_S2203, _S2206) + mean_14;
    float _S2207 = -0.5f - sz_5;
    float3  _S2208 = make_float3 (_S2197 * _S2207, - _S2199, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_0(_S2203, _S2208) + mean_14;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_16, vert0_1) + t_15;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_16, vert1_1) + t_15;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_16, vert2_1) + t_15;
    CameraDistortion_0 _S2209 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_17, tangential_coeffs_17, thin_prism_coeffs_17);
    float2  _S2210 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S2211 = length_0(_S2210);
    float _S2212 = vert0_c_5.z;
    float _S2213 = s_primal_ctx_atan2_0(_S2211, _S2212);
    bool _S2214 = _S2213 < 0.00100000004749745f;
    float k_6;
    float _S2215;
    float _S2216;
    float _S2217;
    if(_S2214)
    {
        float _S2218 = 1.0f - _S2213 * _S2213 / 3.0f;
        float _S2219 = _S2212 * _S2212;
        k_6 = _S2218 / _S2212;
        _S2215 = 0.0f;
        _S2216 = _S2219;
        _S2217 = _S2218;
    }
    else
    {
        float _S2220 = _S2211 * _S2211;
        k_6 = _S2213 / _S2211;
        _S2215 = _S2220;
        _S2216 = 0.0f;
        _S2217 = 0.0f;
    }
    float2  _S2221 = make_float2 (k_6);
    float2  _S2222 = _S2210 * make_float2 (k_6);
    float k1_5 = _S2209.radial_coeffs_0.x;
    float k2_5 = _S2209.radial_coeffs_0.y;
    float k3_5 = _S2209.radial_coeffs_0.z;
    float k4_5 = _S2209.radial_coeffs_0.w;
    float p1_5 = _S2209.tangential_coeffs_0.x;
    float p2_5 = _S2209.tangential_coeffs_0.y;
    float sx1_5 = _S2209.thin_prism_coeffs_0.x;
    float sy1_5 = _S2209.thin_prism_coeffs_0.y;
    float u_15 = _S2222.x;
    float v_15 = _S2222.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float _S2223 = k3_5 + r2_15 * k4_5;
    float _S2224 = k2_5 + r2_15 * _S2223;
    float _S2225 = k1_5 + r2_15 * _S2224;
    float radial_1 = 1.0f + r2_15 * _S2225;
    float2  _S2226 = make_float2 (radial_1);
    float _S2227 = 2.0f * p1_5;
    float _S2228 = _S2227 * u_15;
    float _S2229 = 2.0f * u_15;
    float _S2230 = r2_15 + _S2229 * u_15;
    float _S2231 = 2.0f * p2_5;
    float _S2232 = _S2231 * u_15;
    float _S2233 = 2.0f * v_15;
    float _S2234 = r2_15 + _S2233 * v_15;
    float2  _S2235 = _S2222 * make_float2 (radial_1) + make_float2 (_S2228 * v_15 + p2_5 * _S2230 + sx1_5 * r2_15, _S2232 * v_15 + p1_5 * _S2234 + sy1_5 * r2_15);
    float _S2236 = fx_18 * _S2235.x + cx_18;
    float _S2237 = fy_18 * _S2235.y + cy_18;
    float2  _S2238 = make_float2 (_S2236, _S2237);
    float2  _S2239 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S2240 = length_0(_S2239);
    float _S2241 = vert1_c_5.z;
    float _S2242 = s_primal_ctx_atan2_0(_S2240, _S2241);
    bool _S2243 = _S2242 < 0.00100000004749745f;
    float _S2244;
    float _S2245;
    float _S2246;
    if(_S2243)
    {
        float _S2247 = 1.0f - _S2242 * _S2242 / 3.0f;
        float _S2248 = _S2241 * _S2241;
        k_6 = _S2247 / _S2241;
        _S2244 = 0.0f;
        _S2245 = _S2248;
        _S2246 = _S2247;
    }
    else
    {
        float _S2249 = _S2240 * _S2240;
        k_6 = _S2242 / _S2240;
        _S2244 = _S2249;
        _S2245 = 0.0f;
        _S2246 = 0.0f;
    }
    float2  _S2250 = make_float2 (k_6);
    float2  _S2251 = _S2239 * make_float2 (k_6);
    float u_16 = _S2251.x;
    float v_16 = _S2251.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S2252 = k3_5 + r2_16 * k4_5;
    float _S2253 = k2_5 + r2_16 * _S2252;
    float _S2254 = k1_5 + r2_16 * _S2253;
    float radial_2 = 1.0f + r2_16 * _S2254;
    float2  _S2255 = make_float2 (radial_2);
    float _S2256 = _S2227 * u_16;
    float _S2257 = 2.0f * u_16;
    float _S2258 = r2_16 + _S2257 * u_16;
    float _S2259 = _S2231 * u_16;
    float _S2260 = 2.0f * v_16;
    float _S2261 = r2_16 + _S2260 * v_16;
    float2  _S2262 = _S2251 * make_float2 (radial_2) + make_float2 (_S2256 * v_16 + p2_5 * _S2258 + sx1_5 * r2_16, _S2259 * v_16 + p1_5 * _S2261 + sy1_5 * r2_16);
    float _S2263 = fx_18 * _S2262.x + cx_18;
    float _S2264 = fy_18 * _S2262.y + cy_18;
    float2  _S2265 = make_float2 (_S2263, _S2264);
    float2  _S2266 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S2267 = length_0(_S2266);
    float _S2268 = vert2_c_5.z;
    float _S2269 = s_primal_ctx_atan2_0(_S2267, _S2268);
    bool _S2270 = _S2269 < 0.00100000004749745f;
    float _S2271;
    float _S2272;
    float _S2273;
    if(_S2270)
    {
        float _S2274 = 1.0f - _S2269 * _S2269 / 3.0f;
        float _S2275 = _S2268 * _S2268;
        k_6 = _S2274 / _S2268;
        _S2271 = 0.0f;
        _S2272 = _S2275;
        _S2273 = _S2274;
    }
    else
    {
        float _S2276 = _S2267 * _S2267;
        k_6 = _S2269 / _S2267;
        _S2271 = _S2276;
        _S2272 = 0.0f;
        _S2273 = 0.0f;
    }
    float2  _S2277 = make_float2 (k_6);
    float2  _S2278 = _S2266 * make_float2 (k_6);
    float u_17 = _S2278.x;
    float v_17 = _S2278.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S2279 = k3_5 + r2_17 * k4_5;
    float _S2280 = k2_5 + r2_17 * _S2279;
    float _S2281 = k1_5 + r2_17 * _S2280;
    float radial_3 = 1.0f + r2_17 * _S2281;
    float2  _S2282 = make_float2 (radial_3);
    float _S2283 = _S2227 * u_17;
    float _S2284 = 2.0f * u_17;
    float _S2285 = r2_17 + _S2284 * u_17;
    float _S2286 = _S2231 * u_17;
    float _S2287 = 2.0f * v_17;
    float _S2288 = r2_17 + _S2287 * v_17;
    float2  _S2289 = _S2278 * make_float2 (radial_3) + make_float2 (_S2283 * v_17 + p2_5 * _S2285 + sx1_5 * r2_17, _S2286 * v_17 + p1_5 * _S2288 + sy1_5 * r2_17);
    float _S2290 = fx_18 * _S2289.x + cx_18;
    float _S2291 = fy_18 * _S2289.y + cy_18;
    float2  _S2292 = make_float2 (_S2290, _S2291);
    float2  e0_5 = _S2265 - _S2238;
    float2  e1_5 = _S2292 - _S2265;
    float2  e2_1 = _S2238 - _S2292;
    float _S2293 = e0_5.x;
    float _S2294 = e1_5.y;
    float _S2295 = e0_5.y;
    float _S2296 = e1_5.x;
    float _S2297 = _S2293 * _S2294 - _S2295 * _S2296;
    float _S2298 = 1.0f - hardness_5.y;
    float _S2299 = -1.0f / _S2298;
    float _S2300 = _S2298 * _S2298;
    float _S2301 = s_primal_ctx_max_0(_S2236, _S2263);
    float _S2302 = s_primal_ctx_min_0(_S2236, _S2263);
    float _S2303 = s_primal_ctx_max_0(_S2237, _S2264);
    float _S2304 = s_primal_ctx_min_0(_S2237, _S2264);
    Matrix<float, 3, 3>  _S2305 = transpose_0(R_16);
    float3  _S2306 = mean_14 - - s_primal_ctx_mul_0(_S2305, t_15);
    float _S2307 = _S2306.x;
    float _S2308 = _S2306.y;
    float _S2309 = _S2306.z;
    float _S2310 = _S2307 * _S2307 + _S2308 * _S2308 + _S2309 * _S2309;
    float _S2311 = s_primal_ctx_sqrt_0(_S2310);
    float x_56 = _S2307 / _S2311;
    float3  _S2312 = make_float3 (x_56);
    float _S2313 = _S2311 * _S2311;
    float y_35 = _S2308 / _S2311;
    float z_32 = _S2309 / _S2311;
    float3  _S2314 = make_float3 (z_32);
    float _S2315 = - y_35;
    float3  _S2316 = make_float3 (_S2315);
    float z2_32 = z_32 * z_32;
    float fTmp0B_14 = -1.09254848957061768f * z_32;
    float fC1_14 = x_56 * x_56 - y_35 * y_35;
    float _S2317 = 2.0f * x_56;
    float fS1_14 = _S2317 * y_35;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2318 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_56;
    float3  _S2319 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_35;
    float3  _S2320 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2321 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2322 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_32;
    float _S2323 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_32 * _S2323;
    float3  _S2324 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_56;
    float3  _S2325 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_35;
    float3  _S2326 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2327 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2328 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_56 * fC1_14 - y_35 * fS1_14);
    float3  _S2329 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_56 * fS1_14 + y_35 * fC1_14);
    float3  _S2330 = make_float3 (pSH9_4);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2315) * (*sh_coeffs_14)[int(1)] + make_float3 (z_32) * (*sh_coeffs_14)[int(2)] - make_float3 (x_56) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]);
    float3  _S2331 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S2332 = make_float3 (0.0f);
    float3  _S2333 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S2334 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S2335 = make_float3 (_S2334);
    float3  _S2336 = (*ch_coeffs_5)[int(1)] * make_float3 (_S2334);
    float3  _S2337 = _S2333 + _S2336 + make_float3 (0.5f);
    float3  _S2338 = _S2333 - _S2336 + make_float3 (0.5f);
    float3  _S2339 = vert1_c_5 - vert0_c_5;
    float3  _S2340 = vert2_c_5 - vert0_c_5;
    float3  _S2341 = s_primal_ctx_cross_0(_S2339, _S2340);
    float3  _S2342 = normalize_0(_S2341);
    float3  _S2343 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S2342, mean_c_14)))))) * v_normal_1;
    float3  _S2344 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2345;
    (&_S2345)->primal_0 = _S2342;
    (&_S2345)->differential_0 = _S2344;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2346;
    (&_S2346)->primal_0 = mean_c_14;
    (&_S2346)->differential_0 = _S2344;
    s_bwd_prop_dot_0(&_S2345, &_S2346, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2347 = _S2346;
    float3  _S2348 = _S2343 + _S2345.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2349;
    (&_S2349)->primal_0 = _S2341;
    (&_S2349)->differential_0 = _S2344;
    s_bwd_normalize_impl_0(&_S2349, _S2348);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2350;
    (&_S2350)->primal_0 = _S2339;
    (&_S2350)->differential_0 = _S2344;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2351;
    (&_S2351)->primal_0 = _S2340;
    (&_S2351)->differential_0 = _S2344;
    s_bwd_prop_cross_0(&_S2350, &_S2351, _S2349.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2352 = _S2350;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2353 = _S2351;
    float3  _S2354 = - _S2351.differential_0;
    float3  _S2355 = - _S2350.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2356;
    (&_S2356)->primal_0 = _S2338;
    (&_S2356)->differential_0 = _S2344;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2357;
    (&_S2357)->primal_0 = _S2332;
    (&_S2357)->differential_0 = _S2344;
    s_bwd_prop_max_0(&_S2356, &_S2357, (*v_rgb_4)[int(2)]);
    float3  _S2358 = - _S2356.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2359;
    (&_S2359)->primal_0 = _S2337;
    (&_S2359)->differential_0 = _S2344;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2360;
    (&_S2360)->primal_0 = _S2332;
    (&_S2360)->differential_0 = _S2344;
    s_bwd_prop_max_0(&_S2359, &_S2360, (*v_rgb_4)[int(1)]);
    float3  _S2361 = _S2335 * (_S2358 + _S2359.differential_0);
    float3  _S2362 = _S2356.differential_0 + _S2359.differential_0;
    float3  _S2363 = make_float3 (0.5f) * - _S2362;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2364;
    (&_S2364)->primal_0 = _S2331;
    (&_S2364)->differential_0 = _S2344;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2365;
    (&_S2365)->primal_0 = _S2332;
    (&_S2365)->differential_0 = _S2344;
    s_bwd_prop_max_0(&_S2364, &_S2365, (*v_rgb_4)[int(0)]);
    float3  _S2366 = _S2363 + _S2364.differential_0;
    float3  _S2367 = _S2362 + _S2364.differential_0;
    float3  _S2368 = _S2329 * _S2367;
    float3  _S2369 = (*sh_coeffs_14)[int(15)] * _S2367;
    float3  _S2370 = _S2327 * _S2367;
    float3  _S2371 = (*sh_coeffs_14)[int(14)] * _S2367;
    float3  _S2372 = _S2325 * _S2367;
    float3  _S2373 = (*sh_coeffs_14)[int(13)] * _S2367;
    float3  _S2374 = _S2324 * _S2367;
    float3  _S2375 = (*sh_coeffs_14)[int(12)] * _S2367;
    float3  _S2376 = _S2326 * _S2367;
    float3  _S2377 = (*sh_coeffs_14)[int(11)] * _S2367;
    float3  _S2378 = _S2328 * _S2367;
    float3  _S2379 = (*sh_coeffs_14)[int(10)] * _S2367;
    float3  _S2380 = _S2330 * _S2367;
    float3  _S2381 = (*sh_coeffs_14)[int(9)] * _S2367;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2381.x + _S2381.y + _S2381.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2369.x + _S2369.y + _S2369.z);
    float _S2382 = _S2379.x + _S2379.y + _S2379.z;
    float _S2383 = _S2371.x + _S2371.y + _S2371.z;
    float _S2384 = _S2377.x + _S2377.y + _S2377.z;
    float _S2385 = _S2373.x + _S2373.y + _S2373.z;
    float _S2386 = _S2375.x + _S2375.y + _S2375.z;
    float _S2387 = - s_diff_fC2_T_4;
    float3  _S2388 = _S2321 * _S2367;
    float3  _S2389 = (*sh_coeffs_14)[int(8)] * _S2367;
    float3  _S2390 = _S2319 * _S2367;
    float3  _S2391 = (*sh_coeffs_14)[int(7)] * _S2367;
    float3  _S2392 = _S2318 * _S2367;
    float3  _S2393 = (*sh_coeffs_14)[int(6)] * _S2367;
    float3  _S2394 = _S2320 * _S2367;
    float3  _S2395 = (*sh_coeffs_14)[int(5)] * _S2367;
    float3  _S2396 = _S2322 * _S2367;
    float3  _S2397 = (*sh_coeffs_14)[int(4)] * _S2367;
    float _S2398 = _S2395.x + _S2395.y + _S2395.z;
    float _S2399 = _S2391.x + _S2391.y + _S2391.z;
    float _S2400 = fTmp1B_14 * _S2382 + x_56 * s_diff_fS2_T_4 + y_35 * _S2387 + 0.54627424478530884f * (_S2397.x + _S2397.y + _S2397.z);
    float _S2401 = fTmp1B_14 * _S2383 + y_35 * s_diff_fS2_T_4 + x_56 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2389.x + _S2389.y + _S2389.z);
    float _S2402 = y_35 * - _S2401;
    float _S2403 = x_56 * _S2401;
    float _S2404 = z_32 * (1.86588168144226074f * (z_32 * _S2386) + -2.28522896766662598f * (y_35 * _S2384 + x_56 * _S2385) + 0.94617468118667603f * (_S2393.x + _S2393.y + _S2393.z));
    float3  _S2405 = make_float3 (0.48860251903533936f) * _S2367;
    float3  _S2406 = - _S2405;
    float3  _S2407 = _S2312 * _S2406;
    float3  _S2408 = (*sh_coeffs_14)[int(3)] * _S2406;
    float3  _S2409 = _S2314 * _S2405;
    float3  _S2410 = (*sh_coeffs_14)[int(2)] * _S2405;
    float3  _S2411 = _S2316 * _S2405;
    float3  _S2412 = (*sh_coeffs_14)[int(1)] * _S2405;
    float _S2413 = (_S2323 * _S2386 + 1.44530570507049561f * (fS1_14 * _S2382 + fC1_14 * _S2383) + -1.09254848957061768f * (y_35 * _S2398 + x_56 * _S2399) + _S2404 + _S2404 + _S2410.x + _S2410.y + _S2410.z) / _S2313;
    float _S2414 = _S2311 * _S2413;
    float _S2415 = (fTmp0C_14 * _S2384 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2387 + fTmp0B_14 * _S2398 + _S2317 * _S2400 + _S2402 + _S2402 + - (_S2412.x + _S2412.y + _S2412.z)) / _S2313;
    float _S2416 = _S2311 * _S2415;
    float _S2417 = (fTmp0C_14 * _S2385 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2399 + 2.0f * (y_35 * _S2400) + _S2403 + _S2403 + _S2408.x + _S2408.y + _S2408.z) / _S2313;
    float _S2418 = _S2311 * _S2417;
    float _S2419 = _S2309 * - _S2413 + _S2308 * - _S2415 + _S2307 * - _S2417;
    DiffPair_float_0 _S2420;
    (&_S2420)->primal_0 = _S2310;
    (&_S2420)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2420, _S2419);
    float _S2421 = _S2309 * _S2420.differential_0;
    float _S2422 = _S2308 * _S2420.differential_0;
    float _S2423 = _S2307 * _S2420.differential_0;
    float3  _S2424 = make_float3 (0.282094806432724f) * _S2367;
    float3  _S2425 = make_float3 (_S2418 + _S2423 + _S2423, _S2416 + _S2422 + _S2422, _S2414 + _S2421 + _S2421);
    float3  _S2426 = - - _S2425;
    Matrix<float, 3, 3>  _S2427 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2428;
    (&_S2428)->primal_0 = _S2305;
    (&_S2428)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2429;
    (&_S2429)->primal_0 = t_15;
    (&_S2429)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2428, &_S2429, _S2426);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2430 = _S2429;
    Matrix<float, 3, 3>  _S2431 = transpose_0(_S2428.differential_0);
    DiffPair_float_0 _S2432;
    (&_S2432)->primal_0 = _S2304;
    (&_S2432)->differential_0 = 0.0f;
    DiffPair_float_0 _S2433;
    (&_S2433)->primal_0 = _S2291;
    (&_S2433)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2432, &_S2433, 0.0f);
    DiffPair_float_0 _S2434;
    (&_S2434)->primal_0 = _S2237;
    (&_S2434)->differential_0 = 0.0f;
    DiffPair_float_0 _S2435;
    (&_S2435)->primal_0 = _S2264;
    (&_S2435)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2434, &_S2435, _S2432.differential_0);
    DiffPair_float_0 _S2436;
    (&_S2436)->primal_0 = _S2303;
    (&_S2436)->differential_0 = 0.0f;
    DiffPair_float_0 _S2437;
    (&_S2437)->primal_0 = _S2291;
    (&_S2437)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2436, &_S2437, 0.0f);
    DiffPair_float_0 _S2438;
    (&_S2438)->primal_0 = _S2237;
    (&_S2438)->differential_0 = 0.0f;
    DiffPair_float_0 _S2439;
    (&_S2439)->primal_0 = _S2264;
    (&_S2439)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2438, &_S2439, _S2436.differential_0);
    DiffPair_float_0 _S2440;
    (&_S2440)->primal_0 = _S2302;
    (&_S2440)->differential_0 = 0.0f;
    DiffPair_float_0 _S2441;
    (&_S2441)->primal_0 = _S2290;
    (&_S2441)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2440, &_S2441, 0.0f);
    DiffPair_float_0 _S2442;
    (&_S2442)->primal_0 = _S2236;
    (&_S2442)->differential_0 = 0.0f;
    DiffPair_float_0 _S2443;
    (&_S2443)->primal_0 = _S2263;
    (&_S2443)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2442, &_S2443, _S2440.differential_0);
    DiffPair_float_0 _S2444;
    (&_S2444)->primal_0 = _S2301;
    (&_S2444)->differential_0 = 0.0f;
    DiffPair_float_0 _S2445;
    (&_S2445)->primal_0 = _S2290;
    (&_S2445)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2444, &_S2445, 0.0f);
    DiffPair_float_0 _S2446;
    (&_S2446)->primal_0 = _S2236;
    (&_S2446)->differential_0 = 0.0f;
    DiffPair_float_0 _S2447;
    (&_S2447)->primal_0 = _S2263;
    (&_S2447)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2446, &_S2447, _S2444.differential_0);
    DiffPair_float_0 _S2448;
    (&_S2448)->primal_0 = _S2299;
    (&_S2448)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2448, 0.0f);
    float _S2449 = - (-1.0f * - (_S2448.differential_0 / _S2300));
    float2  _S2450 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2451;
    (&_S2451)->primal_0 = e2_1;
    (&_S2451)->differential_0 = _S2450;
    s_bwd_length_impl_0(&_S2451, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2452;
    (&_S2452)->primal_0 = e1_5;
    (&_S2452)->differential_0 = _S2450;
    s_bwd_length_impl_0(&_S2452, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2453;
    (&_S2453)->primal_0 = e0_5;
    (&_S2453)->differential_0 = _S2450;
    s_bwd_length_impl_0(&_S2453, -0.0f);
    DiffPair_float_0 _S2454;
    (&_S2454)->primal_0 = _S2297;
    (&_S2454)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2454, 0.0f);
    float _S2455 = - _S2454.differential_0;
    float2  _S2456 = _S2452.differential_0 + make_float2 (_S2295 * _S2455, _S2293 * _S2454.differential_0);
    float2  _S2457 = _S2453.differential_0 + make_float2 (_S2294 * _S2454.differential_0, _S2296 * _S2455);
    float2  _S2458 = v_uv2_1 + - _S2451.differential_0 + _S2456;
    float _S2459 = fy_18 * (_S2433.differential_0 + _S2437.differential_0 + _S2458.y);
    float _S2460 = fx_18 * (_S2441.differential_0 + _S2445.differential_0 + _S2458.x);
    float2  _S2461 = make_float2 (_S2460, _S2459);
    float2  _S2462 = _S2278 * _S2461;
    float2  _S2463 = _S2282 * _S2461;
    float _S2464 = r2_17 * _S2459;
    float _S2465 = p1_5 * _S2459;
    float _S2466 = _S2288 * _S2459;
    float _S2467 = v_17 * _S2459;
    float _S2468 = u_17 * _S2467;
    float _S2469 = r2_17 * _S2460;
    float _S2470 = p2_5 * _S2460;
    float _S2471 = _S2285 * _S2460;
    float _S2472 = v_17 * _S2460;
    float _S2473 = u_17 * _S2472;
    float _S2474 = _S2462.x + _S2462.y;
    float _S2475 = r2_17 * _S2474;
    float _S2476 = r2_17 * _S2475;
    float _S2477 = r2_17 * _S2476;
    float _S2478 = r2_17 * _S2477;
    float _S2479 = sy1_5 * _S2459 + _S2465 + sx1_5 * _S2460 + _S2470 + _S2281 * _S2474 + _S2280 * _S2475 + _S2279 * _S2476 + k4_5 * _S2477;
    float _S2480 = v_17 * _S2479;
    float _S2481 = u_17 * _S2479;
    float _S2482 = _S2287 * _S2465 + 2.0f * (v_17 * _S2465) + _S2286 * _S2459 + _S2283 * _S2460 + _S2480 + _S2480;
    float _S2483 = _S2231 * _S2467 + _S2284 * _S2470 + 2.0f * (u_17 * _S2470) + _S2227 * _S2472 + _S2481 + _S2481;
    float2  _S2484 = v_uv0_1 + _S2451.differential_0 + - _S2457;
    float2  _S2485 = v_out_hardness_1 + make_float2 (0.0f, _S2449);
    float _S2486 = _S2443.differential_0 + _S2447.differential_0;
    float2  _S2487 = v_uv1_1 + - _S2456 + _S2457;
    float3  _S2488 = _S2354 + _S2355;
    FixedArray<float3 , 2>  _S2489;
    _S2489[int(0)] = _S2344;
    _S2489[int(1)] = _S2344;
    _S2489[int(1)] = _S2361;
    _S2489[int(0)] = _S2366;
    float3  _S2490 = _S2489[int(0)];
    float3  _S2491 = _S2489[int(1)];
    FixedArray<float3 , 16>  _S2492;
    _S2492[int(0)] = _S2344;
    _S2492[int(1)] = _S2344;
    _S2492[int(2)] = _S2344;
    _S2492[int(3)] = _S2344;
    _S2492[int(4)] = _S2344;
    _S2492[int(5)] = _S2344;
    _S2492[int(6)] = _S2344;
    _S2492[int(7)] = _S2344;
    _S2492[int(8)] = _S2344;
    _S2492[int(9)] = _S2344;
    _S2492[int(10)] = _S2344;
    _S2492[int(11)] = _S2344;
    _S2492[int(12)] = _S2344;
    _S2492[int(13)] = _S2344;
    _S2492[int(14)] = _S2344;
    _S2492[int(15)] = _S2344;
    _S2492[int(7)] = _S2390;
    _S2492[int(0)] = _S2424;
    _S2492[int(1)] = _S2411;
    _S2492[int(2)] = _S2409;
    _S2492[int(3)] = _S2407;
    _S2492[int(4)] = _S2396;
    _S2492[int(5)] = _S2394;
    _S2492[int(6)] = _S2392;
    _S2492[int(15)] = _S2368;
    _S2492[int(8)] = _S2388;
    _S2492[int(9)] = _S2380;
    _S2492[int(10)] = _S2378;
    _S2492[int(11)] = _S2376;
    _S2492[int(12)] = _S2374;
    _S2492[int(13)] = _S2372;
    _S2492[int(14)] = _S2370;
    float3  _S2493 = _S2492[int(0)];
    float3  _S2494 = _S2492[int(1)];
    float3  _S2495 = _S2492[int(2)];
    float3  _S2496 = _S2492[int(3)];
    float3  _S2497 = _S2492[int(4)];
    float3  _S2498 = _S2492[int(5)];
    float3  _S2499 = _S2492[int(6)];
    float3  _S2500 = _S2492[int(7)];
    float3  _S2501 = _S2492[int(8)];
    float3  _S2502 = _S2492[int(9)];
    float3  _S2503 = _S2492[int(10)];
    float3  _S2504 = _S2492[int(11)];
    float3  _S2505 = _S2492[int(12)];
    float3  _S2506 = _S2492[int(13)];
    float3  _S2507 = _S2492[int(14)];
    float3  _S2508 = _S2492[int(15)];
    float _S2509 = _S2442.differential_0 + _S2446.differential_0;
    float _S2510 = _S2434.differential_0 + _S2438.differential_0;
    float _S2511 = _S2435.differential_0 + _S2439.differential_0;
    float2  _S2512 = _S2463 + make_float2 (_S2483, _S2482);
    float2  _S2513 = _S2266 * _S2512;
    float2  _S2514 = _S2277 * _S2512;
    float _S2515 = _S2513.x + _S2513.y;
    if(_S2270)
    {
        float _S2516 = _S2515 / _S2272;
        float _S2517 = _S2273 * - _S2516;
        float _S2518 = _S2269 * (0.3333333432674408f * - (_S2268 * _S2516));
        k_6 = _S2518 + _S2518;
        _S2271 = _S2517;
        _S2272 = 0.0f;
    }
    else
    {
        float _S2519 = _S2515 / _S2271;
        float _S2520 = _S2269 * - _S2519;
        k_6 = _S2267 * _S2519;
        _S2271 = 0.0f;
        _S2272 = _S2520;
    }
    DiffPair_float_0 _S2521;
    (&_S2521)->primal_0 = _S2267;
    (&_S2521)->differential_0 = 0.0f;
    DiffPair_float_0 _S2522;
    (&_S2522)->primal_0 = _S2268;
    (&_S2522)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2521, &_S2522, k_6);
    float _S2523 = _S2522.differential_0 + _S2271;
    float _S2524 = _S2521.differential_0 + _S2272;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2525;
    (&_S2525)->primal_0 = _S2266;
    (&_S2525)->differential_0 = _S2450;
    s_bwd_length_impl_0(&_S2525, _S2524);
    float2  _S2526 = _S2525.differential_0 + _S2514;
    float _S2527 = fy_18 * (_S2487.y + _S2511);
    float _S2528 = fx_18 * (_S2487.x + _S2486);
    float2  _S2529 = make_float2 (_S2528, _S2527);
    float2  _S2530 = _S2251 * _S2529;
    float _S2531 = p1_5 * _S2527;
    float _S2532 = v_16 * _S2527;
    float _S2533 = p2_5 * _S2528;
    float _S2534 = v_16 * _S2528;
    float _S2535 = _S2530.x + _S2530.y;
    float _S2536 = r2_16 * _S2535;
    float _S2537 = r2_16 * _S2536;
    float _S2538 = r2_16 * _S2537;
    float _S2539 = sy1_5 * _S2527 + _S2531 + sx1_5 * _S2528 + _S2533 + _S2254 * _S2535 + _S2253 * _S2536 + _S2252 * _S2537 + k4_5 * _S2538;
    float _S2540 = v_16 * _S2539;
    float _S2541 = u_16 * _S2539;
    float3  _S2542 = _S2353.differential_0 + make_float3 (_S2526.x, _S2526.y, _S2523);
    float2  _S2543 = _S2255 * _S2529 + make_float2 (_S2231 * _S2532 + _S2257 * _S2533 + 2.0f * (u_16 * _S2533) + _S2227 * _S2534 + _S2541 + _S2541, _S2260 * _S2531 + 2.0f * (v_16 * _S2531) + _S2259 * _S2527 + _S2256 * _S2528 + _S2540 + _S2540);
    float _S2544 = u_16 * _S2532 + _S2468;
    float _S2545 = u_16 * _S2534 + _S2473;
    float _S2546 = r2_16 * _S2527 + _S2464;
    float _S2547 = r2_16 * _S2538 + _S2478;
    float _S2548 = _S2538 + _S2477;
    float _S2549 = _S2261 * _S2527 + _S2466;
    float _S2550 = _S2537 + _S2476;
    float _S2551 = r2_16 * _S2528 + _S2469;
    float _S2552 = _S2258 * _S2528 + _S2471;
    float _S2553 = _S2536 + _S2475;
    float2  _S2554 = _S2239 * _S2543;
    float2  _S2555 = _S2250 * _S2543;
    float _S2556 = _S2554.x + _S2554.y;
    if(_S2243)
    {
        float _S2557 = _S2556 / _S2245;
        float _S2558 = _S2246 * - _S2557;
        float _S2559 = _S2242 * (0.3333333432674408f * - (_S2241 * _S2557));
        k_6 = _S2559 + _S2559;
        _S2244 = _S2558;
        _S2245 = 0.0f;
    }
    else
    {
        float _S2560 = _S2556 / _S2244;
        float _S2561 = _S2242 * - _S2560;
        k_6 = _S2240 * _S2560;
        _S2244 = 0.0f;
        _S2245 = _S2561;
    }
    DiffPair_float_0 _S2562;
    (&_S2562)->primal_0 = _S2240;
    (&_S2562)->differential_0 = 0.0f;
    DiffPair_float_0 _S2563;
    (&_S2563)->primal_0 = _S2241;
    (&_S2563)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2562, &_S2563, k_6);
    float _S2564 = _S2563.differential_0 + _S2244;
    float _S2565 = _S2562.differential_0 + _S2245;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2566;
    (&_S2566)->primal_0 = _S2239;
    (&_S2566)->differential_0 = _S2450;
    s_bwd_length_impl_0(&_S2566, _S2565);
    float2  _S2567 = _S2566.differential_0 + _S2555;
    float _S2568 = fy_18 * (_S2484.y + _S2510);
    float _S2569 = fx_18 * (_S2484.x + _S2509);
    float2  _S2570 = make_float2 (_S2569, _S2568);
    float2  _S2571 = _S2222 * _S2570;
    float _S2572 = p1_5 * _S2568;
    float _S2573 = v_15 * _S2568;
    float _S2574 = p2_5 * _S2569;
    float _S2575 = v_15 * _S2569;
    float _S2576 = _S2571.x + _S2571.y;
    float _S2577 = r2_15 * _S2576;
    float _S2578 = r2_15 * _S2577;
    float _S2579 = r2_15 * _S2578;
    float _S2580 = sy1_5 * _S2568 + _S2572 + sx1_5 * _S2569 + _S2574 + _S2225 * _S2576 + _S2224 * _S2577 + _S2223 * _S2578 + k4_5 * _S2579;
    float _S2581 = v_15 * _S2580;
    float _S2582 = u_15 * _S2580;
    float2  _S2583 = make_float2 (r2_15 * _S2569 + _S2551, r2_15 * _S2568 + _S2546);
    float2  _S2584 = make_float2 (_S2234 * _S2568 + 2.0f * (u_15 * _S2575 + _S2545) + _S2549, 2.0f * (u_15 * _S2573 + _S2544) + _S2230 * _S2569 + _S2552);
    float4  _S2585 = make_float4 (_S2577 + _S2553, _S2578 + _S2550, _S2579 + _S2548, r2_15 * _S2579 + _S2547);
    float3  _S2586 = _S2352.differential_0 + make_float3 (_S2567.x, _S2567.y, _S2564);
    float2  _S2587 = _S2226 * _S2570 + make_float2 (_S2231 * _S2573 + _S2229 * _S2574 + 2.0f * (u_15 * _S2574) + _S2227 * _S2575 + _S2582 + _S2582, _S2233 * _S2572 + 2.0f * (v_15 * _S2572) + _S2232 * _S2568 + _S2228 * _S2569 + _S2581 + _S2581);
    CameraDistortion_0 _S2588 = CameraDistortion_x24_syn_dzero_0();
    (&_S2588)->thin_prism_coeffs_0 = _S2583;
    (&_S2588)->tangential_coeffs_0 = _S2584;
    (&_S2588)->radial_coeffs_0 = _S2585;
    float2  _S2589 = _S2210 * _S2587;
    float2  _S2590 = _S2221 * _S2587;
    float _S2591 = _S2589.x + _S2589.y;
    if(_S2214)
    {
        float _S2592 = _S2591 / _S2216;
        float _S2593 = _S2217 * - _S2592;
        float _S2594 = _S2213 * (0.3333333432674408f * - (_S2212 * _S2592));
        k_6 = _S2594 + _S2594;
        _S2215 = _S2593;
        _S2216 = 0.0f;
    }
    else
    {
        float _S2595 = _S2591 / _S2215;
        float _S2596 = _S2213 * - _S2595;
        k_6 = _S2211 * _S2595;
        _S2215 = 0.0f;
        _S2216 = _S2596;
    }
    DiffPair_float_0 _S2597;
    (&_S2597)->primal_0 = _S2211;
    (&_S2597)->differential_0 = 0.0f;
    DiffPair_float_0 _S2598;
    (&_S2598)->primal_0 = _S2212;
    (&_S2598)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S2597, &_S2598, k_6);
    float _S2599 = _S2598.differential_0 + _S2215;
    float _S2600 = _S2597.differential_0 + _S2216;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2601;
    (&_S2601)->primal_0 = _S2210;
    (&_S2601)->differential_0 = _S2450;
    s_bwd_length_impl_0(&_S2601, _S2600);
    float2  _S2602 = _S2601.differential_0 + _S2590;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2603;
    (&_S2603)->primal_0 = vert2_c_5;
    (&_S2603)->differential_0 = _S2344;
    s_bwd_length_impl_1(&_S2603, v_depth_4.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2604;
    (&_S2604)->primal_0 = vert1_c_5;
    (&_S2604)->differential_0 = _S2344;
    s_bwd_length_impl_1(&_S2604, v_depth_4.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2605;
    (&_S2605)->primal_0 = vert0_c_5;
    (&_S2605)->differential_0 = _S2344;
    s_bwd_length_impl_1(&_S2605, v_depth_4.x);
    float3  _S2606 = _S2603.differential_0 + _S2542;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2607;
    (&_S2607)->primal_0 = R_16;
    (&_S2607)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2608;
    (&_S2608)->primal_0 = vert2_1;
    (&_S2608)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2607, &_S2608, _S2606);
    float3  _S2609 = _S2604.differential_0 + _S2586;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2610;
    (&_S2610)->primal_0 = R_16;
    (&_S2610)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2611;
    (&_S2611)->primal_0 = vert1_1;
    (&_S2611)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2610, &_S2611, _S2609);
    float3  _S2612 = _S2605.differential_0 + _S2488 + make_float3 (_S2602.x, _S2602.y, _S2599);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2613;
    (&_S2613)->primal_0 = R_16;
    (&_S2613)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2614;
    (&_S2614)->primal_0 = vert0_1;
    (&_S2614)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2613, &_S2614, _S2612);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2615;
    (&_S2615)->primal_0 = _S2203;
    (&_S2615)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2616;
    (&_S2616)->primal_0 = _S2208;
    (&_S2616)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2615, &_S2616, _S2608.differential_0);
    float _S2617 = - _S2616.differential_0.y;
    float _S2618 = _S2207 * _S2616.differential_0.x;
    float _S2619 = - (_S2197 * _S2616.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2620;
    (&_S2620)->primal_0 = _S2203;
    (&_S2620)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2621;
    (&_S2621)->primal_0 = _S2206;
    (&_S2621)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2620, &_S2621, _S2611.differential_0);
    float _S2622 = _S2197 * _S2621.differential_0.x;
    float _S2623 = _S2205 * _S2621.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2624;
    (&_S2624)->primal_0 = _S2203;
    (&_S2624)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2625;
    (&_S2625)->primal_0 = _S2204;
    (&_S2625)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2624, &_S2625, _S2614.differential_0);
    Matrix<float, 3, 3>  _S2626 = transpose_0(_S2615.differential_0 + _S2620.differential_0 + _S2624.differential_0);
    float _S2627 = 2.0f * - _S2626.rows[int(2)].z;
    float _S2628 = 2.0f * _S2626.rows[int(2)].y;
    float _S2629 = 2.0f * _S2626.rows[int(2)].x;
    float _S2630 = 2.0f * _S2626.rows[int(1)].z;
    float _S2631 = 2.0f * - _S2626.rows[int(1)].y;
    float _S2632 = 2.0f * _S2626.rows[int(1)].x;
    float _S2633 = 2.0f * _S2626.rows[int(0)].z;
    float _S2634 = 2.0f * _S2626.rows[int(0)].y;
    float _S2635 = 2.0f * - _S2626.rows[int(0)].x;
    float _S2636 = - _S2632 + _S2634;
    float _S2637 = _S2629 + - _S2633;
    float _S2638 = - _S2628 + _S2630;
    float _S2639 = _S2628 + _S2630;
    float _S2640 = _S2629 + _S2633;
    float _S2641 = _S2632 + _S2634;
    float _S2642 = z_31 * (_S2631 + _S2635);
    float _S2643 = y_34 * (_S2627 + _S2635);
    float _S2644 = x_55 * (_S2627 + _S2631);
    float _S2645 = z_31 * _S2636 + y_34 * _S2637 + x_55 * _S2638;
    float _S2646 = _S2202 * _S2645;
    float _S2647 = w_17 * _S2636 + y_34 * _S2639 + x_55 * _S2640 + _S2642 + _S2642;
    float _S2648 = _S2202 * _S2647;
    float _S2649 = w_17 * _S2637 + z_31 * _S2639 + x_55 * _S2641 + _S2643 + _S2643;
    float _S2650 = _S2202 * _S2649;
    float _S2651 = w_17 * _S2638 + z_31 * _S2640 + y_34 * _S2641 + _S2644 + _S2644;
    float _S2652 = _S2202 * _S2651;
    float _S2653 = quat_17.x * _S2645 + quat_17.w * _S2647 + quat_17.z * _S2649 + quat_17.y * _S2651;
    DiffPair_float_0 _S2654;
    (&_S2654)->primal_0 = _S2201;
    (&_S2654)->differential_0 = 0.0f;
    s_bwd_prop_rsqrt_0(&_S2654, _S2653);
    float _S2655 = quat_17.x * _S2654.differential_0;
    float _S2656 = quat_17.w * _S2654.differential_0;
    float _S2657 = quat_17.z * _S2654.differential_0;
    float _S2658 = quat_17.y * _S2654.differential_0;
    float _S2659 = _S2648 + _S2656 + _S2656;
    float _S2660 = _S2650 + _S2657 + _S2657;
    float _S2661 = _S2652 + _S2658 + _S2658;
    float _S2662 = _S2646 + _S2655 + _S2655;
    float _S2663 = _S2619 + _S2622;
    float _S2664 = 0.5f * - _S2663;
    float _S2665 = _S2617 + _S2621.differential_0.y;
    DiffPair_float_0 _S2666;
    (&_S2666)->primal_0 = _S2198;
    (&_S2666)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2666, _S2665);
    float _S2667 = _S2664 + _S2666.differential_0;
    float _S2668 = _S2618 + _S2623 + _S2625.differential_0.x;
    DiffPair_float_0 _S2669;
    (&_S2669)->primal_0 = _S2196;
    (&_S2669)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2669, _S2668);
    float _S2670 = _S2664 + _S2669.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2671;
    (&_S2671)->primal_0 = mean_c_14;
    (&_S2671)->differential_0 = _S2344;
    s_bwd_length_impl_1(&_S2671, 0.0f);
    float3  _S2672 = _S2671.differential_0 + _S2347.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2673;
    (&_S2673)->primal_0 = R_16;
    (&_S2673)->differential_0 = _S2427;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2674;
    (&_S2674)->primal_0 = mean_14;
    (&_S2674)->differential_0 = _S2344;
    s_bwd_prop_mul_0(&_S2673, &_S2674, _S2672);
    float3  _S2675 = _S2606 + _S2609 + _S2612 + _S2672 + _S2430.differential_0;
    Matrix<float, 3, 3>  _S2676 = _S2607.differential_0 + _S2610.differential_0 + _S2613.differential_0 + _S2673.differential_0 + _S2431;
    float3  _S2677 = make_float3 (_S2670, _S2667, _S2663);
    float4  _S2678 = make_float4 (0.0f);
    *&((&_S2678)->w) = _S2659;
    *&((&_S2678)->z) = _S2660;
    *&((&_S2678)->y) = _S2661;
    *&((&_S2678)->x) = _S2662;
    float4  _S2679 = _S2678;
    float3  _S2680 = _S2608.differential_0 + _S2611.differential_0 + _S2614.differential_0 + _S2674.differential_0 + _S2425;
    *v_mean_4 = _S2680;
    *v_quat_4 = _S2679;
    *v_scale_4 = _S2677;
    *v_hardness_1 = _S2485;
    (*v_sh_coeffs_4)[int(0)] = _S2493;
    (*v_sh_coeffs_4)[int(1)] = _S2494;
    (*v_sh_coeffs_4)[int(2)] = _S2495;
    (*v_sh_coeffs_4)[int(3)] = _S2496;
    (*v_sh_coeffs_4)[int(4)] = _S2497;
    (*v_sh_coeffs_4)[int(5)] = _S2498;
    (*v_sh_coeffs_4)[int(6)] = _S2499;
    (*v_sh_coeffs_4)[int(7)] = _S2500;
    (*v_sh_coeffs_4)[int(8)] = _S2501;
    (*v_sh_coeffs_4)[int(9)] = _S2502;
    (*v_sh_coeffs_4)[int(10)] = _S2503;
    (*v_sh_coeffs_4)[int(11)] = _S2504;
    (*v_sh_coeffs_4)[int(12)] = _S2505;
    (*v_sh_coeffs_4)[int(13)] = _S2506;
    (*v_sh_coeffs_4)[int(14)] = _S2507;
    (*v_sh_coeffs_4)[int(15)] = _S2508;
    (*v_ch_coeffs_1)[int(0)] = _S2490;
    (*v_ch_coeffs_1)[int(1)] = _S2491;
    *v_R_4 = _S2676;
    *v_t_4 = _S2675;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_18)
{
    DiffPair_float_0 _S2681 = *dpx_16;
    bool _S2682;
    if(((*dpx_16).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S2682 = ((*dpx_16).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S2682 = false;
    }
    float _S2683;
    if(_S2682)
    {
        _S2683 = dOut_18;
    }
    else
    {
        _S2683 = 0.0f;
    }
    dpx_16->primal_0 = _S2681.primal_0;
    dpx_16->differential_0 = _S2683;
    DiffPair_float_0 _S2684 = *dpMin_0;
    if((_S2681.primal_0) < ((*dpMin_0).primal_0))
    {
        _S2683 = dOut_18;
    }
    else
    {
        _S2683 = 0.0f;
    }
    dpMin_0->primal_0 = _S2684.primal_0;
    dpMin_0->differential_0 = _S2683;
    DiffPair_float_0 _S2685 = *dpMax_0;
    if(((*dpx_16).primal_0) > ((*dpMax_0).primal_0))
    {
        _S2683 = dOut_18;
    }
    else
    {
        _S2683 = 0.0f;
    }
    dpMax_0->primal_0 = _S2685.primal_0;
    dpMax_0->differential_0 = _S2683;
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
        DiffPair_float_0 _S2686 = *dpx_17;
        float _S2687 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_19;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S2687;
        float _S2688 = val_0 * (F32_log((_S2686.primal_0))) * dOut_19;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S2688;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S2689 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S2689))));
    float2  _S2690 = p_0 - v0_0;
    float2  _S2691 = normalize_1(e0_6);
    float2  _S2692 = p_0 - v1_0;
    float2  _S2693 = normalize_1(e1_6);
    float2  _S2694 = p_0 - v2_0;
    float2  _S2695 = normalize_1(e2_2);
    float _S2696 = hardness_6.x;
    float _S2697 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S2690.x * _S2691.y - _S2690.y * _S2691.x)), (se_0 * (_S2692.x * _S2693.y - _S2692.y * _S2693.x))))), (se_0 * (_S2694.x * _S2695.y - _S2694.y * _S2695.x)))) / ((F32_abs((_S2689))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S2697))));
    float _S2698;
    if(a_1 <= 0.0f)
    {
        _S2698 = 0.0f;
    }
    else
    {
        _S2698 = (F32_min(((F32_pow((a_1), (_S2697)))), (0.99900001287460327f)));
    }
    return _S2696 * _S2698;
}

inline __device__ float s_primal_ctx_abs_0(float _S2699)
{
    return (F32_abs((_S2699)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S2700, float _S2701, float _S2702)
{
    return clamp_0(_S2700, _S2701, _S2702);
}

inline __device__ float s_primal_ctx_exp2_0(float _S2703)
{
    return (F32_exp2((_S2703)));
}

inline __device__ float s_primal_ctx_pow_0(float _S2704, float _S2705)
{
    return (F32_pow((_S2704), (_S2705)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S2706, DiffPair_float_0 * _S2707, float _S2708)
{
    _d_pow_0(_S2706, _S2707, _S2708);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S2709, DiffPair_float_0 * _S2710, DiffPair_float_0 * _S2711, float _S2712)
{
    _d_clamp_0(_S2709, _S2710, _S2711, _S2712);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_7)
{
    float _S2713 = length_0((*dpx_18).primal_0);
    float2  _S2714 = (*dpx_18).primal_0 * _s_dOut_7;
    float2  _S2715 = make_float2 (1.0f / _S2713) * _s_dOut_7;
    float _S2716 = - ((_S2714.x + _S2714.y) / (_S2713 * _S2713));
    float2  _S2717 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2718;
    (&_S2718)->primal_0 = (*dpx_18).primal_0;
    (&_S2718)->differential_0 = _S2717;
    s_bwd_length_impl_0(&_S2718, _S2716);
    float2  _S2719 = _S2715 + _S2718.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S2719;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2720, float2  _S2721)
{
    s_bwd_prop_normalize_impl_1(_S2720, _S2721);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_8)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S2722 = e0_7.x;
    float _S2723 = e1_7.y;
    float _S2724 = e0_7.y;
    float _S2725 = e1_7.x;
    float _S2726 = _S2722 * _S2723 - _S2724 * _S2725;
    float se_1 = float((F32_sign((_S2726))));
    float2  _S2727 = p_1 - (*dpv0_0).primal_0;
    float2  _S2728 = normalize_1(e0_7);
    float _S2729 = _S2727.x;
    float _S2730 = _S2728.y;
    float _S2731 = _S2727.y;
    float _S2732 = _S2728.x;
    float de0_0 = se_1 * (_S2729 * _S2730 - _S2731 * _S2732);
    float2  _S2733 = p_1 - (*dpv1_0).primal_0;
    float2  _S2734 = normalize_1(e1_7);
    float _S2735 = _S2733.x;
    float _S2736 = _S2734.y;
    float _S2737 = _S2733.y;
    float _S2738 = _S2734.x;
    float de1_0 = se_1 * (_S2735 * _S2736 - _S2737 * _S2738);
    float2  _S2739 = p_1 - (*dpv2_0).primal_0;
    float2  _S2740 = normalize_1(e2_3);
    float _S2741 = _S2739.x;
    float _S2742 = _S2740.y;
    float _S2743 = _S2739.y;
    float _S2744 = _S2740.x;
    float de2_0 = se_1 * (_S2741 * _S2742 - _S2743 * _S2744);
    float _S2745 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S2746 = s_primal_ctx_max_0(_S2745, de2_0);
    float _S2747 = s_primal_ctx_abs_0(_S2726);
    float _S2748 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S2747 / _S2748;
    float _S2749 = _S2748 * _S2748;
    float _S2750 = (*dphardness_0).primal_0.x;
    float _S2751 = (*dphardness_0).primal_0.y;
    float _S2752 = dmax_0 * dmax_0;
    float _S2753 = 1.0f + _S2746 / dmax_0;
    float _S2754 = 1.0f - s_primal_ctx_clamp_0(_S2751, 0.00499999988824129f, 0.98000001907348633f);
    float _S2755 = -1.0f / _S2754;
    float _S2756 = _S2754 * _S2754;
    float _S2757 = 1.0f - s_primal_ctx_exp2_0(_S2755);
    float a_2 = 1.0f - _S2753 * _S2757;
    bool _S2758 = a_2 <= 0.0f;
    float _S2759;
    float _S2760;
    if(_S2758)
    {
        _S2759 = 0.0f;
        _S2760 = 0.0f;
    }
    else
    {
        float _S2761 = s_primal_ctx_pow_0(a_2, _S2754);
        _S2759 = s_primal_ctx_min_0(_S2761, 0.99900001287460327f);
        _S2760 = _S2761;
    }
    float _S2762 = _S2750 * _s_dOut_8;
    float _S2763 = _S2759 * _s_dOut_8;
    if(_S2758)
    {
        _S2759 = 0.0f;
        _S2760 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S2764;
        (&_S2764)->primal_0 = _S2760;
        (&_S2764)->differential_0 = 0.0f;
        DiffPair_float_0 _S2765;
        (&_S2765)->primal_0 = 0.99900001287460327f;
        (&_S2765)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S2764, &_S2765, _S2762);
        DiffPair_float_0 _S2766;
        (&_S2766)->primal_0 = a_2;
        (&_S2766)->differential_0 = 0.0f;
        DiffPair_float_0 _S2767;
        (&_S2767)->primal_0 = _S2754;
        (&_S2767)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S2766, &_S2767, _S2764.differential_0);
        _S2759 = _S2766.differential_0;
        _S2760 = _S2767.differential_0;
    }
    float _S2768 = - _S2759;
    float _S2769 = _S2757 * _S2768;
    float _S2770 = - (_S2753 * _S2768);
    DiffPair_float_0 _S2771;
    (&_S2771)->primal_0 = _S2755;
    (&_S2771)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2771, _S2770);
    float _S2772 = - (-1.0f * - (_S2771.differential_0 / _S2756) + _S2760);
    float _S2773 = _S2769 / _S2752;
    float s_diff_dmax_T_0 = _S2746 * - _S2773;
    float _S2774 = dmax_0 * _S2773;
    DiffPair_float_0 _S2775;
    (&_S2775)->primal_0 = _S2751;
    (&_S2775)->differential_0 = 0.0f;
    DiffPair_float_0 _S2776;
    (&_S2776)->primal_0 = 0.00499999988824129f;
    (&_S2776)->differential_0 = 0.0f;
    DiffPair_float_0 _S2777;
    (&_S2777)->primal_0 = 0.98000001907348633f;
    (&_S2777)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2775, &_S2776, &_S2777, _S2772);
    float _S2778 = s_diff_dmax_T_0 / _S2749;
    float _S2779 = _S2747 * - _S2778;
    float _S2780 = _S2748 * _S2778;
    float2  _S2781 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2782;
    (&_S2782)->primal_0 = e2_3;
    (&_S2782)->differential_0 = _S2781;
    s_bwd_length_impl_0(&_S2782, _S2779);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2783;
    (&_S2783)->primal_0 = e1_7;
    (&_S2783)->differential_0 = _S2781;
    s_bwd_length_impl_0(&_S2783, _S2779);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2784;
    (&_S2784)->primal_0 = e0_7;
    (&_S2784)->differential_0 = _S2781;
    s_bwd_length_impl_0(&_S2784, _S2779);
    DiffPair_float_0 _S2785;
    (&_S2785)->primal_0 = _S2726;
    (&_S2785)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2785, _S2780);
    DiffPair_float_0 _S2786;
    (&_S2786)->primal_0 = _S2745;
    (&_S2786)->differential_0 = 0.0f;
    DiffPair_float_0 _S2787;
    (&_S2787)->primal_0 = de2_0;
    (&_S2787)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2786, &_S2787, _S2774);
    DiffPair_float_0 _S2788;
    (&_S2788)->primal_0 = de0_0;
    (&_S2788)->differential_0 = 0.0f;
    DiffPair_float_0 _S2789;
    (&_S2789)->primal_0 = de1_0;
    (&_S2789)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2788, &_S2789, _S2786.differential_0);
    float _S2790 = se_1 * _S2787.differential_0;
    float _S2791 = - _S2790;
    float _S2792 = _S2744 * _S2791;
    float _S2793 = _S2742 * _S2790;
    float2  _S2794 = make_float2 (_S2743 * _S2791, _S2741 * _S2790);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2795;
    (&_S2795)->primal_0 = e2_3;
    (&_S2795)->differential_0 = _S2781;
    s_bwd_normalize_impl_1(&_S2795, _S2794);
    float2  _S2796 = - make_float2 (_S2793, _S2792);
    float _S2797 = se_1 * _S2789.differential_0;
    float _S2798 = - _S2797;
    float _S2799 = _S2738 * _S2798;
    float _S2800 = _S2736 * _S2797;
    float2  _S2801 = make_float2 (_S2737 * _S2798, _S2735 * _S2797);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2802;
    (&_S2802)->primal_0 = e1_7;
    (&_S2802)->differential_0 = _S2781;
    s_bwd_normalize_impl_1(&_S2802, _S2801);
    float2  _S2803 = - make_float2 (_S2800, _S2799);
    float _S2804 = se_1 * _S2788.differential_0;
    float _S2805 = - _S2804;
    float _S2806 = _S2732 * _S2805;
    float _S2807 = _S2730 * _S2804;
    float2  _S2808 = make_float2 (_S2731 * _S2805, _S2729 * _S2804);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2809;
    (&_S2809)->primal_0 = e0_7;
    (&_S2809)->differential_0 = _S2781;
    s_bwd_normalize_impl_1(&_S2809, _S2808);
    float2  _S2810 = - make_float2 (_S2807, _S2806);
    float _S2811 = - _S2785.differential_0;
    float2  _S2812 = _S2782.differential_0 + _S2795.differential_0;
    float2  _S2813 = - _S2812;
    float2  _S2814 = _S2783.differential_0 + _S2802.differential_0 + make_float2 (_S2724 * _S2811, _S2722 * _S2785.differential_0);
    float2  _S2815 = - _S2814;
    float2  _S2816 = _S2784.differential_0 + _S2809.differential_0 + make_float2 (_S2723 * _S2785.differential_0, _S2725 * _S2811);
    float2  _S2817 = - _S2816;
    float2  _S2818 = make_float2 (_S2763, _S2775.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S2818;
    float2  _S2819 = _S2796 + _S2813 + _S2814;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S2819;
    float2  _S2820 = _S2803 + _S2815 + _S2816;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S2820;
    float2  _S2821 = _S2810 + _S2812 + _S2817;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S2821;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2822, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2823, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2824, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2825, float2  _S2826, float _S2827)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S2822, _S2823, _S2824, _S2825, _S2826, _S2827);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_0, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S2828 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S2828;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S2828;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S2828;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S2828;
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
    float2  _S2829 = p_3 - v0_2;
    float2  _S2830 = p_3 - v1_2;
    float2  _S2831 = p_3 - v2_2;
    float _S2832 = e0_8.x;
    float _S2833 = e1_8.y;
    float _S2834 = e0_8.y;
    float _S2835 = e1_8.x;
    float _S2836 = _S2832 * _S2833 - _S2834 * _S2835;
    float se_2 = float((F32_sign((_S2836))));
    float _S2837 = hardness_8.x;
    float _S2838 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S2829.x * _S2834 - _S2829.y * _S2832)), (se_2 * (_S2830.x * _S2833 - _S2830.y * _S2835))))), (se_2 * (_S2831.x * e2_4.y - _S2831.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S2829 - e0_8 * make_float2 (clamp_0(dot_1(_S2829, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S2830 - e1_8 * make_float2 (clamp_0(dot_1(_S2830, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S2831 - e2_4 * make_float2 (clamp_0(dot_1(_S2831, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S2836))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S2838))));
    float _S2839;
    if(a_3 <= 0.0f)
    {
        _S2839 = 0.0f;
    }
    else
    {
        _S2839 = (F32_min(((F32_pow((a_3), (_S2838)))), (0.99900001287460327f)));
    }
    return _S2837 * _S2839;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S2840, float2  _S2841)
{
    return dot_1(_S2840, _S2841);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2842, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S2843, float _S2844)
{
    _d_dot_1(_S2842, _S2843, _S2844);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_9)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S2845 = p_4 - (*dpv0_1).primal_0;
    float _S2846 = s_primal_ctx_dot_1(_S2845, e0_9);
    float _S2847 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S2848 = _S2846 / _S2847;
    float _S2849 = _S2847 * _S2847;
    float _S2850 = s_primal_ctx_clamp_0(_S2848, 0.0f, 1.0f);
    float2  _S2851 = make_float2 (_S2850);
    float2  _S2852 = _S2845 - e0_9 * make_float2 (_S2850);
    float _S2853 = length_0(_S2852);
    float2  _S2854 = p_4 - (*dpv1_1).primal_0;
    float _S2855 = s_primal_ctx_dot_1(_S2854, e1_9);
    float _S2856 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S2857 = _S2855 / _S2856;
    float _S2858 = _S2856 * _S2856;
    float _S2859 = s_primal_ctx_clamp_0(_S2857, 0.0f, 1.0f);
    float2  _S2860 = make_float2 (_S2859);
    float2  _S2861 = _S2854 - e1_9 * make_float2 (_S2859);
    float _S2862 = length_0(_S2861);
    float2  _S2863 = p_4 - (*dpv2_1).primal_0;
    float _S2864 = s_primal_ctx_dot_1(_S2863, e2_5);
    float _S2865 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S2866 = _S2864 / _S2865;
    float _S2867 = _S2865 * _S2865;
    float _S2868 = s_primal_ctx_clamp_0(_S2866, 0.0f, 1.0f);
    float2  _S2869 = make_float2 (_S2868);
    float2  _S2870 = _S2863 - e2_5 * make_float2 (_S2868);
    float _S2871 = length_0(_S2870);
    float _S2872 = e0_9.x;
    float _S2873 = e1_9.y;
    float _S2874 = e0_9.y;
    float _S2875 = e1_9.x;
    float _S2876 = _S2872 * _S2873 - _S2874 * _S2875;
    float se_3 = float((F32_sign((_S2876))));
    float _S2877 = _S2845.x;
    float _S2878 = _S2845.y;
    float s0_0 = se_3 * (_S2877 * _S2874 - _S2878 * _S2872);
    float _S2879 = _S2854.x;
    float _S2880 = _S2854.y;
    float s1_0 = se_3 * (_S2879 * _S2873 - _S2880 * _S2875);
    float _S2881 = _S2863.x;
    float _S2882 = e2_5.y;
    float _S2883 = _S2863.y;
    float _S2884 = e2_5.x;
    float s2_0 = se_3 * (_S2881 * _S2882 - _S2883 * _S2884);
    float _S2885 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S2885, s2_0)))));
    float _S2886 = s_primal_ctx_min_0(_S2853, _S2862);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S2886, _S2871);
    float _S2887 = s_primal_ctx_abs_0(_S2876);
    float _S2888 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S2887 / _S2888;
    float _S2889 = _S2888 * _S2888;
    float _S2890 = (*dphardness_1).primal_0.x;
    float _S2891 = (*dphardness_1).primal_0.y;
    float _S2892 = dmax_1 * dmax_1;
    float _S2893 = 1.0f + dv_0 / dmax_1;
    float _S2894 = 1.0f - s_primal_ctx_clamp_0(_S2891, 0.00499999988824129f, 0.98000001907348633f);
    float _S2895 = -1.0f / _S2894;
    float _S2896 = _S2894 * _S2894;
    float _S2897 = 1.0f - s_primal_ctx_exp2_0(_S2895);
    float a_4 = 1.0f - _S2893 * _S2897;
    bool _S2898 = a_4 <= 0.0f;
    float _S2899;
    float _S2900;
    if(_S2898)
    {
        _S2899 = 0.0f;
        _S2900 = 0.0f;
    }
    else
    {
        float _S2901 = s_primal_ctx_pow_0(a_4, _S2894);
        _S2899 = s_primal_ctx_min_0(_S2901, 0.99900001287460327f);
        _S2900 = _S2901;
    }
    float _S2902 = _S2890 * _s_dOut_9;
    float _S2903 = _S2899 * _s_dOut_9;
    if(_S2898)
    {
        _S2899 = 0.0f;
        _S2900 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S2904;
        (&_S2904)->primal_0 = _S2900;
        (&_S2904)->differential_0 = 0.0f;
        DiffPair_float_0 _S2905;
        (&_S2905)->primal_0 = 0.99900001287460327f;
        (&_S2905)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S2904, &_S2905, _S2902);
        DiffPair_float_0 _S2906;
        (&_S2906)->primal_0 = a_4;
        (&_S2906)->differential_0 = 0.0f;
        DiffPair_float_0 _S2907;
        (&_S2907)->primal_0 = _S2894;
        (&_S2907)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S2906, &_S2907, _S2904.differential_0);
        _S2899 = _S2906.differential_0;
        _S2900 = _S2907.differential_0;
    }
    float _S2908 = - _S2899;
    float _S2909 = _S2897 * _S2908;
    float _S2910 = - (_S2893 * _S2908);
    DiffPair_float_0 _S2911;
    (&_S2911)->primal_0 = _S2895;
    (&_S2911)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2911, _S2910);
    float _S2912 = - (-1.0f * - (_S2911.differential_0 / _S2896) + _S2900);
    float _S2913 = _S2909 / _S2892;
    float s_diff_dmax_T_1 = dv_0 * - _S2913;
    float s_diff_dv_T_0 = dmax_1 * _S2913;
    DiffPair_float_0 _S2914;
    (&_S2914)->primal_0 = _S2891;
    (&_S2914)->differential_0 = 0.0f;
    DiffPair_float_0 _S2915;
    (&_S2915)->primal_0 = 0.00499999988824129f;
    (&_S2915)->differential_0 = 0.0f;
    DiffPair_float_0 _S2916;
    (&_S2916)->primal_0 = 0.98000001907348633f;
    (&_S2916)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2914, &_S2915, &_S2916, _S2912);
    float _S2917 = s_diff_dmax_T_1 / _S2889;
    float _S2918 = _S2887 * - _S2917;
    float _S2919 = _S2888 * _S2917;
    float2  _S2920 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2921;
    (&_S2921)->primal_0 = e2_5;
    (&_S2921)->differential_0 = _S2920;
    s_bwd_length_impl_0(&_S2921, _S2918);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2922;
    (&_S2922)->primal_0 = e1_9;
    (&_S2922)->differential_0 = _S2920;
    s_bwd_length_impl_0(&_S2922, _S2918);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2923;
    (&_S2923)->primal_0 = e0_9;
    (&_S2923)->differential_0 = _S2920;
    s_bwd_length_impl_0(&_S2923, _S2918);
    DiffPair_float_0 _S2924;
    (&_S2924)->primal_0 = _S2876;
    (&_S2924)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2924, _S2919);
    float _S2925 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S2926;
    (&_S2926)->primal_0 = _S2886;
    (&_S2926)->differential_0 = 0.0f;
    DiffPair_float_0 _S2927;
    (&_S2927)->primal_0 = _S2871;
    (&_S2927)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2926, &_S2927, _S2925);
    DiffPair_float_0 _S2928;
    (&_S2928)->primal_0 = _S2853;
    (&_S2928)->differential_0 = 0.0f;
    DiffPair_float_0 _S2929;
    (&_S2929)->primal_0 = _S2862;
    (&_S2929)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2928, &_S2929, _S2926.differential_0);
    DiffPair_float_0 _S2930;
    (&_S2930)->primal_0 = _S2885;
    (&_S2930)->differential_0 = 0.0f;
    DiffPair_float_0 _S2931;
    (&_S2931)->primal_0 = s2_0;
    (&_S2931)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2930, &_S2931, 0.0f);
    DiffPair_float_0 _S2932;
    (&_S2932)->primal_0 = s0_0;
    (&_S2932)->differential_0 = 0.0f;
    DiffPair_float_0 _S2933;
    (&_S2933)->primal_0 = s1_0;
    (&_S2933)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2932, &_S2933, _S2930.differential_0);
    float _S2934 = se_3 * _S2931.differential_0;
    float _S2935 = - _S2934;
    float _S2936 = _S2883 * _S2935;
    float _S2937 = _S2884 * _S2935;
    float _S2938 = _S2881 * _S2934;
    float _S2939 = _S2882 * _S2934;
    float _S2940 = se_3 * _S2933.differential_0;
    float _S2941 = - _S2940;
    float _S2942 = _S2875 * _S2941;
    float _S2943 = _S2873 * _S2940;
    float _S2944 = se_3 * _S2932.differential_0;
    float _S2945 = - _S2944;
    float _S2946 = _S2872 * _S2945;
    float _S2947 = _S2874 * _S2944;
    float _S2948 = - _S2924.differential_0;
    float _S2949 = _S2880 * _S2941 + _S2874 * _S2948;
    float _S2950 = _S2877 * _S2944 + _S2875 * _S2948;
    float _S2951 = _S2879 * _S2940 + _S2872 * _S2924.differential_0;
    float _S2952 = _S2878 * _S2945 + _S2873 * _S2924.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2953;
    (&_S2953)->primal_0 = _S2870;
    (&_S2953)->differential_0 = _S2920;
    s_bwd_length_impl_0(&_S2953, _S2927.differential_0);
    float2  _S2954 = - _S2953.differential_0;
    float2  _S2955 = e2_5 * _S2954;
    float2  _S2956 = _S2869 * _S2954;
    float _S2957 = _S2955.x + _S2955.y;
    DiffPair_float_0 _S2958;
    (&_S2958)->primal_0 = _S2866;
    (&_S2958)->differential_0 = 0.0f;
    DiffPair_float_0 _S2959;
    (&_S2959)->primal_0 = 0.0f;
    (&_S2959)->differential_0 = 0.0f;
    DiffPair_float_0 _S2960;
    (&_S2960)->primal_0 = 1.0f;
    (&_S2960)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2958, &_S2959, &_S2960, _S2957);
    float _S2961 = _S2958.differential_0 / _S2867;
    float _S2962 = _S2864 * - _S2961;
    float _S2963 = _S2865 * _S2961;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2964;
    (&_S2964)->primal_0 = e2_5;
    (&_S2964)->differential_0 = _S2920;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2965;
    (&_S2965)->primal_0 = e2_5;
    (&_S2965)->differential_0 = _S2920;
    s_bwd_prop_dot_1(&_S2964, &_S2965, _S2962);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2966;
    (&_S2966)->primal_0 = _S2863;
    (&_S2966)->differential_0 = _S2920;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2967;
    (&_S2967)->primal_0 = e2_5;
    (&_S2967)->differential_0 = _S2920;
    s_bwd_prop_dot_1(&_S2966, &_S2967, _S2963);
    float2  _S2968 = - (_S2953.differential_0 + _S2966.differential_0 + make_float2 (_S2939, _S2937));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2969;
    (&_S2969)->primal_0 = _S2861;
    (&_S2969)->differential_0 = _S2920;
    s_bwd_length_impl_0(&_S2969, _S2929.differential_0);
    float2  _S2970 = - _S2969.differential_0;
    float2  _S2971 = e1_9 * _S2970;
    float2  _S2972 = _S2860 * _S2970;
    float _S2973 = _S2971.x + _S2971.y;
    DiffPair_float_0 _S2974;
    (&_S2974)->primal_0 = _S2857;
    (&_S2974)->differential_0 = 0.0f;
    DiffPair_float_0 _S2975;
    (&_S2975)->primal_0 = 0.0f;
    (&_S2975)->differential_0 = 0.0f;
    DiffPair_float_0 _S2976;
    (&_S2976)->primal_0 = 1.0f;
    (&_S2976)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2974, &_S2975, &_S2976, _S2973);
    float _S2977 = _S2974.differential_0 / _S2858;
    float _S2978 = _S2855 * - _S2977;
    float _S2979 = _S2856 * _S2977;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2980;
    (&_S2980)->primal_0 = e1_9;
    (&_S2980)->differential_0 = _S2920;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2981;
    (&_S2981)->primal_0 = e1_9;
    (&_S2981)->differential_0 = _S2920;
    s_bwd_prop_dot_1(&_S2980, &_S2981, _S2978);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2982;
    (&_S2982)->primal_0 = _S2854;
    (&_S2982)->differential_0 = _S2920;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2983;
    (&_S2983)->primal_0 = e1_9;
    (&_S2983)->differential_0 = _S2920;
    s_bwd_prop_dot_1(&_S2982, &_S2983, _S2979);
    float2  _S2984 = - (_S2969.differential_0 + _S2982.differential_0 + make_float2 (_S2943, _S2942));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2985;
    (&_S2985)->primal_0 = _S2852;
    (&_S2985)->differential_0 = _S2920;
    s_bwd_length_impl_0(&_S2985, _S2928.differential_0);
    float2  _S2986 = - _S2985.differential_0;
    float2  _S2987 = e0_9 * _S2986;
    float2  _S2988 = _S2851 * _S2986;
    float _S2989 = _S2987.x + _S2987.y;
    DiffPair_float_0 _S2990;
    (&_S2990)->primal_0 = _S2848;
    (&_S2990)->differential_0 = 0.0f;
    DiffPair_float_0 _S2991;
    (&_S2991)->primal_0 = 0.0f;
    (&_S2991)->differential_0 = 0.0f;
    DiffPair_float_0 _S2992;
    (&_S2992)->primal_0 = 1.0f;
    (&_S2992)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S2990, &_S2991, &_S2992, _S2989);
    float _S2993 = _S2990.differential_0 / _S2849;
    float _S2994 = _S2846 * - _S2993;
    float _S2995 = _S2847 * _S2993;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2996;
    (&_S2996)->primal_0 = e0_9;
    (&_S2996)->differential_0 = _S2920;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2997;
    (&_S2997)->primal_0 = e0_9;
    (&_S2997)->differential_0 = _S2920;
    s_bwd_prop_dot_1(&_S2996, &_S2997, _S2994);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2998;
    (&_S2998)->primal_0 = _S2845;
    (&_S2998)->differential_0 = _S2920;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2999;
    (&_S2999)->primal_0 = e0_9;
    (&_S2999)->differential_0 = _S2920;
    s_bwd_prop_dot_1(&_S2998, &_S2999, _S2995);
    float2  _S3000 = - (_S2985.differential_0 + _S2998.differential_0 + make_float2 (_S2947, _S2946));
    float2  _S3001 = _S2921.differential_0 + _S2956 + _S2965.differential_0 + _S2964.differential_0 + _S2967.differential_0 + make_float2 (_S2936, _S2938);
    float2  _S3002 = - _S3001;
    float2  _S3003 = _S2922.differential_0 + _S2972 + _S2981.differential_0 + _S2980.differential_0 + _S2983.differential_0 + make_float2 (_S2949, _S2951);
    float2  _S3004 = - _S3003;
    float2  _S3005 = _S2923.differential_0 + _S2988 + _S2997.differential_0 + _S2996.differential_0 + _S2999.differential_0 + make_float2 (_S2952, _S2950);
    float2  _S3006 = - _S3005;
    float2  _S3007 = make_float2 (_S2903, _S2914.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S3007;
    float2  _S3008 = _S2968 + _S3002 + _S3003;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S3008;
    float2  _S3009 = _S2984 + _S3004 + _S3005;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S3009;
    float2  _S3010 = _S3000 + _S3001 + _S3006;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S3010;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3011, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3012, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3013, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3014, float2  _S3015, float _S3016)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S3011, _S3012, _S3013, _S3014, _S3015, _S3016);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_1, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S3017 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S3017;
    s_bwd_evaluate_alpha_opaque_triangle_precise_0(&dp_v0_1, &dp_v1_1, &dp_v2_1, &dp_hardness_1, p_5, v_alpha_1);
    *v_v0_1 = dp_v0_1.differential_0;
    *v_v1_1 = dp_v2_1.differential_0;
    *v_v2_1 = dp_v1_1.differential_0;
    *v_hardness_3 = dp_hardness_1.differential_0;
    return;
}

inline __device__ void evaluate_color_opaque_triangle(float2  v0_4, float2  v1_4, float2  v2_4, FixedArray<float3 , 3>  * colors_0, float3  depths_0, float2  p_6, float3  * color_6, float * depth_10)
{
    *color_6 = ((*colors_0)[int(0)] + (*colors_0)[int(1)] + (*colors_0)[int(2)]) / make_float3 (3.0f);
    *depth_10 = (depths_0.x + depths_0.y + depths_0.z) / 3.0f;
    return;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_2, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpcolors_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepths_0, float2  p_7, float3  dpcolor_0, float dpdepth_0)
{
    float _S3018 = 0.3333333432674408f * dpdepth_0;
    float3  _S3019 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S3020 = make_float3 (0.0f);
    float3  _S3021 = _S3020;
    *&((&_S3021)->z) = _S3018;
    *&((&_S3021)->y) = _S3018;
    *&((&_S3021)->x) = _S3018;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S3021;
    FixedArray<float3 , 3>  _S3022;
    _S3022[int(0)] = _S3020;
    _S3022[int(1)] = _S3020;
    _S3022[int(2)] = _S3020;
    _S3022[int(2)] = _S3019;
    _S3022[int(1)] = _S3019;
    _S3022[int(0)] = _S3019;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S3022;
    float2  _S3023 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S3023;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S3023;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S3023;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3024, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3025, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3026, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S3027, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3028, float2  _S3029, float3  _S3030, float _S3031)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S3024, _S3025, _S3026, _S3027, _S3028, _S3029, _S3030, _S3031);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_8, float3  v_color_0, float v_depth_5, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S3032 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S3032;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S3032;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S3032;
    float3  _S3033 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S3034 = { _S3033, _S3033, _S3033 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S3034;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S3033;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_8, v_color_0, v_depth_5);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

