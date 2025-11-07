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

inline __device__ Matrix<float, 3, 3>  normalized_quat_to_rotmat(float4  quat_0)
{
    float x_3 = quat_0.y;
    float x2_0 = x_3 * x_3;
    float y2_0 = quat_0.z * quat_0.z;
    float z2_0 = quat_0.w * quat_0.w;
    float xy_0 = quat_0.y * quat_0.z;
    float xz_0 = quat_0.y * quat_0.w;
    float yz_0 = quat_0.z * quat_0.w;
    float wx_0 = quat_0.x * quat_0.y;
    float wy_0 = quat_0.x * quat_0.z;
    float wz_0 = quat_0.x * quat_0.w;
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

inline __device__ void _d_mul_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * right_0, float3  dOut_0)
{
    float _S2 = (*left_0).primal_0.rows[int(0)].x * dOut_0.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_0.x;
    float sum_0 = _S2 + (*left_0).primal_0.rows[int(1)].x * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_0.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_0.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S3 = (*left_0).primal_0.rows[int(0)].y * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_0.x;
    float sum_2 = _S3 + (*left_0).primal_0.rows[int(1)].y * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_0.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_0.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S4 = (*left_0).primal_0.rows[int(0)].z * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_0.x;
    float sum_4 = _S4 + (*left_0).primal_0.rows[int(1)].z * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->z) = (*right_0).primal_0.z * dOut_0.y;
    float sum_5 = sum_4 + (*left_0).primal_0.rows[int(2)].z * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->z) = (*right_0).primal_0.z * dOut_0.z;
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

inline __device__ void mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * left_2, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_2, Matrix<float, 3, 3>  dOut_1)
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
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_1.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(0)].x * dOut_1.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_1.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(0)].y * dOut_1.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_1.rows[int(0)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(0)].z * dOut_1.rows[int(0)].x;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_1.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(0)].x * dOut_1.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_1.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(0)].y * dOut_1.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_1.rows[int(0)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(0)].z * dOut_1.rows[int(0)].y;
    *&(((&left_d_result_1)->rows + (int(0)))->x) = *&(((&left_d_result_1)->rows + (int(0)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_1.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(0)].x * dOut_1.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(0)))->y) = *&(((&left_d_result_1)->rows + (int(0)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_1.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(0)].y * dOut_1.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(0)))->z) = *&(((&left_d_result_1)->rows + (int(0)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_1.rows[int(0)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(0)].z * dOut_1.rows[int(0)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_1.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(1)].x * dOut_1.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_1.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(1)].y * dOut_1.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_1.rows[int(1)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(1)].z * dOut_1.rows[int(1)].x;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_1.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(1)].x * dOut_1.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_1.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(1)].y * dOut_1.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_1.rows[int(1)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(1)].z * dOut_1.rows[int(1)].y;
    *&(((&left_d_result_1)->rows + (int(1)))->x) = *&(((&left_d_result_1)->rows + (int(1)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_1.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(1)].x * dOut_1.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->y) = *&(((&left_d_result_1)->rows + (int(1)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_1.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(1)].y * dOut_1.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(1)))->z) = *&(((&left_d_result_1)->rows + (int(1)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_1.rows[int(1)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(1)].z * dOut_1.rows[int(1)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].x * dOut_1.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(0)))->x) = *&(((&right_d_result_1)->rows + (int(0)))->x) + (*left_2).primal_0.rows[int(2)].x * dOut_1.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].x * dOut_1.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(1)))->x) = *&(((&right_d_result_1)->rows + (int(1)))->x) + (*left_2).primal_0.rows[int(2)].y * dOut_1.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].x * dOut_1.rows[int(2)].x;
    *&(((&right_d_result_1)->rows + (int(2)))->x) = *&(((&right_d_result_1)->rows + (int(2)))->x) + (*left_2).primal_0.rows[int(2)].z * dOut_1.rows[int(2)].x;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].y * dOut_1.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(0)))->y) = *&(((&right_d_result_1)->rows + (int(0)))->y) + (*left_2).primal_0.rows[int(2)].x * dOut_1.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].y * dOut_1.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(1)))->y) = *&(((&right_d_result_1)->rows + (int(1)))->y) + (*left_2).primal_0.rows[int(2)].y * dOut_1.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].y * dOut_1.rows[int(2)].y;
    *&(((&right_d_result_1)->rows + (int(2)))->y) = *&(((&right_d_result_1)->rows + (int(2)))->y) + (*left_2).primal_0.rows[int(2)].z * dOut_1.rows[int(2)].y;
    *&(((&left_d_result_1)->rows + (int(2)))->x) = *&(((&left_d_result_1)->rows + (int(2)))->x) + (*right_2).primal_0.rows[int(0)].z * dOut_1.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(0)))->z) = *&(((&right_d_result_1)->rows + (int(0)))->z) + (*left_2).primal_0.rows[int(2)].x * dOut_1.rows[int(2)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->y) = *&(((&left_d_result_1)->rows + (int(2)))->y) + (*right_2).primal_0.rows[int(1)].z * dOut_1.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(1)))->z) = *&(((&right_d_result_1)->rows + (int(1)))->z) + (*left_2).primal_0.rows[int(2)].y * dOut_1.rows[int(2)].z;
    *&(((&left_d_result_1)->rows + (int(2)))->z) = *&(((&left_d_result_1)->rows + (int(2)))->z) + (*right_2).primal_0.rows[int(2)].z * dOut_1.rows[int(2)].z;
    *&(((&right_d_result_1)->rows + (int(2)))->z) = *&(((&right_d_result_1)->rows + (int(2)))->z) + (*left_2).primal_0.rows[int(2)].z * dOut_1.rows[int(2)].z;
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

inline __device__ void mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_3, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_3, Matrix<float, 2, 3>  dOut_2)
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
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].x * dOut_2.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_3).primal_0.rows[int(0)].x * dOut_2.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].x * dOut_2.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_3).primal_0.rows[int(0)].y * dOut_2.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].x * dOut_2.rows[int(0)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_3).primal_0.rows[int(0)].z * dOut_2.rows[int(0)].x;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].y * dOut_2.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_3).primal_0.rows[int(0)].x * dOut_2.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].y * dOut_2.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_3).primal_0.rows[int(0)].y * dOut_2.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].y * dOut_2.rows[int(0)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_3).primal_0.rows[int(0)].z * dOut_2.rows[int(0)].y;
    *&(((&left_d_result_2)->rows + (int(0)))->x) = *&(((&left_d_result_2)->rows + (int(0)))->x) + (*right_3).primal_0.rows[int(0)].z * dOut_2.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_3).primal_0.rows[int(0)].x * dOut_2.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->y) = *&(((&left_d_result_2)->rows + (int(0)))->y) + (*right_3).primal_0.rows[int(1)].z * dOut_2.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_3).primal_0.rows[int(0)].y * dOut_2.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(0)))->z) = *&(((&left_d_result_2)->rows + (int(0)))->z) + (*right_3).primal_0.rows[int(2)].z * dOut_2.rows[int(0)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_3).primal_0.rows[int(0)].z * dOut_2.rows[int(0)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].x * dOut_2.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(0)))->x) = *&(((&right_d_result_2)->rows + (int(0)))->x) + (*left_3).primal_0.rows[int(1)].x * dOut_2.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].x * dOut_2.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(1)))->x) = *&(((&right_d_result_2)->rows + (int(1)))->x) + (*left_3).primal_0.rows[int(1)].y * dOut_2.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].x * dOut_2.rows[int(1)].x;
    *&(((&right_d_result_2)->rows + (int(2)))->x) = *&(((&right_d_result_2)->rows + (int(2)))->x) + (*left_3).primal_0.rows[int(1)].z * dOut_2.rows[int(1)].x;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].y * dOut_2.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(0)))->y) = *&(((&right_d_result_2)->rows + (int(0)))->y) + (*left_3).primal_0.rows[int(1)].x * dOut_2.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].y * dOut_2.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(1)))->y) = *&(((&right_d_result_2)->rows + (int(1)))->y) + (*left_3).primal_0.rows[int(1)].y * dOut_2.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].y * dOut_2.rows[int(1)].y;
    *&(((&right_d_result_2)->rows + (int(2)))->y) = *&(((&right_d_result_2)->rows + (int(2)))->y) + (*left_3).primal_0.rows[int(1)].z * dOut_2.rows[int(1)].y;
    *&(((&left_d_result_2)->rows + (int(1)))->x) = *&(((&left_d_result_2)->rows + (int(1)))->x) + (*right_3).primal_0.rows[int(0)].z * dOut_2.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(0)))->z) = *&(((&right_d_result_2)->rows + (int(0)))->z) + (*left_3).primal_0.rows[int(1)].x * dOut_2.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->y) = *&(((&left_d_result_2)->rows + (int(1)))->y) + (*right_3).primal_0.rows[int(1)].z * dOut_2.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(1)))->z) = *&(((&right_d_result_2)->rows + (int(1)))->z) + (*left_3).primal_0.rows[int(1)].y * dOut_2.rows[int(1)].z;
    *&(((&left_d_result_2)->rows + (int(1)))->z) = *&(((&left_d_result_2)->rows + (int(1)))->z) + (*right_3).primal_0.rows[int(2)].z * dOut_2.rows[int(1)].z;
    *&(((&right_d_result_2)->rows + (int(2)))->z) = *&(((&right_d_result_2)->rows + (int(2)))->z) + (*left_3).primal_0.rows[int(1)].z * dOut_2.rows[int(1)].z;
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

inline __device__ void mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * left_4, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * right_4, Matrix<float, 2, 2>  dOut_3)
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
    *&(((&left_d_result_3)->rows + (int(0)))->x) = *&(((&left_d_result_3)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_3.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = *&(((&right_d_result_3)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(0)].x * dOut_3.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = *&(((&left_d_result_3)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_3.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = *&(((&right_d_result_3)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(0)].y * dOut_3.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = *&(((&left_d_result_3)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_3.rows[int(0)].x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = *&(((&right_d_result_3)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(0)].z * dOut_3.rows[int(0)].x;
    *&(((&left_d_result_3)->rows + (int(0)))->x) = *&(((&left_d_result_3)->rows + (int(0)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_3.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = *&(((&right_d_result_3)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(0)].x * dOut_3.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(0)))->y) = *&(((&left_d_result_3)->rows + (int(0)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_3.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = *&(((&right_d_result_3)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(0)].y * dOut_3.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(0)))->z) = *&(((&left_d_result_3)->rows + (int(0)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_3.rows[int(0)].y;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = *&(((&right_d_result_3)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(0)].z * dOut_3.rows[int(0)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = *&(((&left_d_result_3)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].x * dOut_3.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(0)))->x) = *&(((&right_d_result_3)->rows + (int(0)))->x) + (*left_4).primal_0.rows[int(1)].x * dOut_3.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = *&(((&left_d_result_3)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].x * dOut_3.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(1)))->x) = *&(((&right_d_result_3)->rows + (int(1)))->x) + (*left_4).primal_0.rows[int(1)].y * dOut_3.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = *&(((&left_d_result_3)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].x * dOut_3.rows[int(1)].x;
    *&(((&right_d_result_3)->rows + (int(2)))->x) = *&(((&right_d_result_3)->rows + (int(2)))->x) + (*left_4).primal_0.rows[int(1)].z * dOut_3.rows[int(1)].x;
    *&(((&left_d_result_3)->rows + (int(1)))->x) = *&(((&left_d_result_3)->rows + (int(1)))->x) + (*right_4).primal_0.rows[int(0)].y * dOut_3.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(0)))->y) = *&(((&right_d_result_3)->rows + (int(0)))->y) + (*left_4).primal_0.rows[int(1)].x * dOut_3.rows[int(1)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->y) = *&(((&left_d_result_3)->rows + (int(1)))->y) + (*right_4).primal_0.rows[int(1)].y * dOut_3.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(1)))->y) = *&(((&right_d_result_3)->rows + (int(1)))->y) + (*left_4).primal_0.rows[int(1)].y * dOut_3.rows[int(1)].y;
    *&(((&left_d_result_3)->rows + (int(1)))->z) = *&(((&left_d_result_3)->rows + (int(1)))->z) + (*right_4).primal_0.rows[int(2)].y * dOut_3.rows[int(1)].y;
    *&(((&right_d_result_3)->rows + (int(2)))->y) = *&(((&right_d_result_3)->rows + (int(2)))->y) + (*left_4).primal_0.rows[int(1)].z * dOut_3.rows[int(1)].y;
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
    float x_4 = quat_1.y;
    float x2_1 = x_4 * x_4;
    float y2_1 = quat_1.z * quat_1.z;
    float z2_1 = quat_1.w * quat_1.w;
    float xy_1 = quat_1.y * quat_1.z;
    float xz_1 = quat_1.y * quat_1.w;
    float yz_1 = quat_1.z * quat_1.w;
    float wx_1 = quat_1.x * quat_1.y;
    float wy_1 = quat_1.x * quat_1.z;
    float wz_1 = quat_1.x * quat_1.w;
    Matrix<float, 3, 3>  M_0 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_1 + z2_1), 2.0f * (xy_1 + wz_1), 2.0f * (xz_1 - wy_1), 2.0f * (xy_1 - wz_1), 1.0f - 2.0f * (x2_1 + z2_1), 2.0f * (yz_1 + wx_1), 2.0f * (xz_1 + wy_1), 2.0f * (yz_1 - wx_1), 1.0f - 2.0f * (x2_1 + y2_1))), makeMatrix<float, 3, 3> (scale_0.x, 0.0f, 0.0f, 0.0f, scale_0.y, 0.0f, 0.0f, 0.0f, scale_0.z));
    *covar_0 = mul_4(M_0, transpose_0(M_0));
    return;
}

inline __device__ void quat_scale_to_sqrt_covar(float4  quat_2, float3  scale_1, Matrix<float, 3, 3>  * M_1)
{
    float x_5 = quat_2.y;
    float x2_2 = x_5 * x_5;
    float y2_2 = quat_2.z * quat_2.z;
    float z2_2 = quat_2.w * quat_2.w;
    float xy_2 = quat_2.y * quat_2.z;
    float xz_2 = quat_2.y * quat_2.w;
    float yz_2 = quat_2.z * quat_2.w;
    float wx_2 = quat_2.x * quat_2.y;
    float wy_2 = quat_2.x * quat_2.z;
    float wz_2 = quat_2.x * quat_2.w;
    *M_1 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_2 + z2_2), 2.0f * (xy_2 + wz_2), 2.0f * (xz_2 - wy_2), 2.0f * (xy_2 - wz_2), 1.0f - 2.0f * (x2_2 + z2_2), 2.0f * (yz_2 + wx_2), 2.0f * (xz_2 + wy_2), 2.0f * (yz_2 - wx_2), 1.0f - 2.0f * (x2_2 + y2_2))), makeMatrix<float, 3, 3> (scale_1.x, 0.0f, 0.0f, 0.0f, scale_1.y, 0.0f, 0.0f, 0.0f, scale_1.z));
    return;
}

inline __device__ Matrix<float, 2, 2>  inverse(Matrix<float, 2, 2>  m_0)
{
    float invdet_0 = 1.0f / (m_0.rows[int(0)].x * m_0.rows[int(1)].y - m_0.rows[int(0)].y * m_0.rows[int(1)].x);
    return makeMatrix<float, 2, 2> (m_0.rows[int(1)].y * invdet_0, - m_0.rows[int(0)].y * invdet_0, - m_0.rows[int(1)].x * invdet_0, m_0.rows[int(0)].x * invdet_0);
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpy_0, float dOut_4)
{
    DiffPair_float_0 _S5 = *dpx_0;
    float _S6;
    if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
    {
        _S6 = dOut_4;
    }
    else
    {
        if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
        {
            _S6 = 0.0f;
        }
        else
        {
            _S6 = 0.5f * dOut_4;
        }
    }
    dpx_0->primal_0 = _S5.primal_0;
    dpx_0->differential_0 = _S6;
    DiffPair_float_0 _S7 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S5.primal_0))
    {
        _S6 = dOut_4;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_0).primal_0))
        {
            _S6 = 0.0f;
        }
        else
        {
            _S6 = 0.5f * dOut_4;
        }
    }
    dpy_0->primal_0 = _S7.primal_0;
    dpy_0->differential_0 = _S6;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_1, float dOut_5)
{
    float _S8 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_1).primal_0)))))) * dOut_5;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = _S8;
    return;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_1, float dOut_6)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_1).primal_0.x * dOut_6;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_2).primal_0.x * dOut_6;
    *&((&x_d_result_0)->y) = (*dpy_1).primal_0.y * dOut_6;
    *&((&y_d_result_0)->y) = (*dpx_2).primal_0.y * dOut_6;
    *&((&x_d_result_0)->z) = (*dpy_1).primal_0.z * dOut_6;
    *&((&y_d_result_0)->z) = (*dpx_2).primal_0.z * dOut_6;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = x_d_result_0;
    dpy_1->primal_0 = (*dpy_1).primal_0;
    dpy_1->differential_0 = y_d_result_0;
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void _d_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpy_2, float dOut_7)
{
    float2  x_d_result_1;
    *&((&x_d_result_1)->x) = (*dpy_2).primal_0.x * dOut_7;
    float2  y_d_result_1;
    *&((&y_d_result_1)->x) = (*dpx_3).primal_0.x * dOut_7;
    *&((&x_d_result_1)->y) = (*dpy_2).primal_0.y * dOut_7;
    *&((&y_d_result_1)->y) = (*dpx_3).primal_0.y * dOut_7;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = x_d_result_1;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_1;
    return;
}

inline __device__ float dot_0(float3  x_6, float3  y_0)
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
        float result_10 = result_9 + _slang_vector_get_element(x_6, i_4) * _slang_vector_get_element(y_0, i_4);
        i_4 = i_4 + int(1);
        result_9 = result_10;
    }
    return result_9;
}

inline __device__ float dot_1(float2  x_7, float2  y_1)
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
        float result_12 = result_11 + _slang_vector_get_element(x_7, i_5) * _slang_vector_get_element(y_1, i_5);
        i_5 = i_5 + int(1);
        result_11 = result_12;
    }
    return result_11;
}

inline __device__ float length_0(float2  x_8)
{
    return (F32_sqrt((dot_1(x_8, x_8))));
}

inline __device__ float length_1(float3  x_9)
{
    return (F32_sqrt((dot_0(x_9, x_9))));
}

inline __device__ CameraDistortion_0 CameraDistortion_x24init_0(float4  radial_coeffs_1, float2  tangential_coeffs_1, float2  thin_prism_coeffs_1)
{
    CameraDistortion_0 _S9;
    (&_S9)->radial_coeffs_0 = radial_coeffs_1;
    (&_S9)->tangential_coeffs_0 = tangential_coeffs_1;
    (&_S9)->thin_prism_coeffs_0 = thin_prism_coeffs_1;
    return _S9;
}

inline __device__ float2  distort_point(float2  uv_0, bool is_fisheye_0, float4  radial_coeffs_2, float2  tangential_coeffs_2, float2  thin_prism_coeffs_2)
{
    float2  _S10;
    if(is_fisheye_0)
    {
        float r_6 = length_0(uv_0);
        float theta_0 = (F32_atan((r_6)));
        float _S11;
        if(r_6 < 0.00100000004749745f)
        {
            _S11 = 1.0f - r_6 * r_6 / 3.0f;
        }
        else
        {
            _S11 = theta_0 / r_6;
        }
        _S10 = uv_0 * make_float2 (_S11);
    }
    else
    {
        _S10 = uv_0;
    }
    CameraDistortion_0 _S12 = CameraDistortion_x24init_0(radial_coeffs_2, tangential_coeffs_2, thin_prism_coeffs_2);
    float p1_0 = _S12.tangential_coeffs_0.x;
    float p2_0 = _S12.tangential_coeffs_0.y;
    float u_0 = _S10.x;
    float v_0 = _S10.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    return _S10 * make_float2 (1.0f + r2_0 * (_S12.radial_coeffs_0.x + r2_0 * (_S12.radial_coeffs_0.y + r2_0 * (_S12.radial_coeffs_0.z + r2_0 * _S12.radial_coeffs_0.w)))) + make_float2 (2.0f * p1_0 * u_0 * v_0 + p2_0 * (r2_0 + 2.0f * u_0 * u_0) + _S12.thin_prism_coeffs_0.x * r2_0, 2.0f * p2_0 * u_0 * v_0 + p1_0 * (r2_0 + 2.0f * v_0 * v_0) + _S12.thin_prism_coeffs_0.y * r2_0);
}

inline __device__ float2  undistort_point_0(float2  uv_1, CameraDistortion_0 * dist_coeffs_0, int maxiter_0)
{
    int i_6 = int(0);
    float2  q_0 = uv_1;
    for(;;)
    {
        if(i_6 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float k4_0 = dist_coeffs_0->radial_coeffs_0.w;
        float p1_1 = dist_coeffs_0->tangential_coeffs_0.x;
        float p2_1 = dist_coeffs_0->tangential_coeffs_0.y;
        float sx1_0 = dist_coeffs_0->thin_prism_coeffs_0.x;
        float sy1_0 = dist_coeffs_0->thin_prism_coeffs_0.y;
        float u_1 = q_0.x;
        float v_1 = q_0.y;
        float r2_1 = u_1 * u_1 + v_1 * v_1;
        float _S13 = dist_coeffs_0->radial_coeffs_0.z + r2_1 * k4_0;
        float _S14 = dist_coeffs_0->radial_coeffs_0.y + r2_1 * _S13;
        float _S15 = dist_coeffs_0->radial_coeffs_0.x + r2_1 * _S14;
        float radial_0 = 1.0f + r2_1 * _S15;
        float _S16 = 2.0f * p1_1;
        float _S17 = _S16 * u_1;
        float _S18 = 2.0f * u_1;
        float _S19 = 2.0f * p2_1;
        float _S20 = _S19 * u_1;
        float _S21 = 2.0f * v_1;
        float2  _S22 = q_0 * make_float2 (radial_0) + make_float2 (_S17 * v_1 + p2_1 * (r2_1 + _S18 * u_1) + sx1_0 * r2_1, _S20 * v_1 + p1_1 * (r2_1 + _S21 * v_1) + sy1_0 * r2_1);
        float2  _S23 = make_float2 (0.0f);
        float2  seed_0 = _S23;
        *&((&seed_0)->x) = 1.0f;
        float2  _S24 = make_float2 (radial_0);
        float2  _S25 = q_0 * seed_0;
        float _S26 = p1_1 * seed_0.y;
        float _S27 = p2_1 * seed_0.x;
        float _S28 = _S25.x + _S25.y;
        float _S29 = r2_1 * _S28;
        float _S30 = r2_1 * _S29;
        float _S31 = sy1_0 * seed_0.y + _S26 + sx1_0 * seed_0.x + _S27 + _S15 * _S28 + _S14 * _S29 + _S13 * _S30 + k4_0 * (r2_1 * _S30);
        float _S32 = v_1 * _S31;
        float _S33 = u_1 * _S31;
        Matrix<float, 2, 2>  J_0;
        J_0[int(0)] = _S24 * seed_0 + make_float2 (_S19 * (v_1 * seed_0.y) + _S18 * _S27 + 2.0f * (u_1 * _S27) + _S16 * (v_1 * seed_0.x) + _S33 + _S33, _S21 * _S26 + 2.0f * (v_1 * _S26) + _S20 * seed_0.y + _S17 * seed_0.x + _S32 + _S32);
        float2  seed_1 = _S23;
        *&((&seed_1)->y) = 1.0f;
        float2  _S34 = q_0 * seed_1;
        float _S35 = p1_1 * seed_1.y;
        float _S36 = p2_1 * seed_1.x;
        float _S37 = _S34.x + _S34.y;
        float _S38 = r2_1 * _S37;
        float _S39 = r2_1 * _S38;
        float _S40 = sy1_0 * seed_1.y + _S35 + sx1_0 * seed_1.x + _S36 + _S15 * _S37 + _S14 * _S38 + _S13 * _S39 + k4_0 * (r2_1 * _S39);
        float _S41 = v_1 * _S40;
        float _S42 = u_1 * _S40;
        J_0[int(1)] = _S24 * seed_1 + make_float2 (_S19 * (v_1 * seed_1.y) + _S18 * _S36 + 2.0f * (u_1 * _S36) + _S16 * (v_1 * seed_1.x) + _S42 + _S42, _S21 * _S35 + 2.0f * (v_1 * _S35) + _S20 * seed_1.y + _S17 * seed_1.x + _S41 + _S41);
        float2  _S43 = _S22 - uv_1;
        float inv_det_0 = 1.0f / (J_0.rows[int(0)].x * J_0.rows[int(1)].y - J_0.rows[int(0)].y * J_0.rows[int(1)].x);
        float _S44 = _S43.x;
        float _S45 = _S43.y;
        float2  q_1 = q_0 - make_float2 ((_S44 * J_0.rows[int(1)].y - _S45 * J_0.rows[int(0)].y) * inv_det_0, (- _S44 * J_0.rows[int(1)].x + _S45 * J_0.rows[int(0)].x) * inv_det_0);
        i_6 = i_6 + int(1);
        q_0 = q_1;
    }
    return q_0;
}

inline __device__ float2  undistort_point(float2  uv_2, bool is_fisheye_1, float4  radial_coeffs_3, float2  tangential_coeffs_3, float2  thin_prism_coeffs_3)
{
    CameraDistortion_0 _S46;
    (&_S46)->radial_coeffs_0 = radial_coeffs_3;
    (&_S46)->tangential_coeffs_0 = tangential_coeffs_3;
    (&_S46)->thin_prism_coeffs_0 = thin_prism_coeffs_3;
    float2  _S47 = undistort_point_0(uv_2, &_S46, int(8));
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float theta_1 = length_0(_S47);
        float _S48;
        if(theta_1 < 0.00100000004749745f)
        {
            _S48 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S48 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S49 = make_float3 ((_S47 * make_float2 (_S48)).x, (_S47 * make_float2 (_S48)).y, (F32_cos((theta_1))));
        raydir_0 = _S49;
    }
    else
    {
        raydir_0 = make_float3 (_S47.x, _S47.y, 1.0f);
    }
    return float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
}

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_8, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_8, float3  dOut_8)
{
    float _S50 = (*right_8).primal_0.rows[int(0)].x * dOut_8.x;
    Matrix<float, 3, 3>  right_d_result_4;
    *&(((&right_d_result_4)->rows + (int(0)))->x) = (*left_8).primal_0.x * dOut_8.x;
    float sum_14 = _S50 + (*right_8).primal_0.rows[int(0)].y * dOut_8.y;
    *&(((&right_d_result_4)->rows + (int(0)))->y) = (*left_8).primal_0.x * dOut_8.y;
    float sum_15 = sum_14 + (*right_8).primal_0.rows[int(0)].z * dOut_8.z;
    *&(((&right_d_result_4)->rows + (int(0)))->z) = (*left_8).primal_0.x * dOut_8.z;
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = sum_15;
    float _S51 = (*right_8).primal_0.rows[int(1)].x * dOut_8.x;
    *&(((&right_d_result_4)->rows + (int(1)))->x) = (*left_8).primal_0.y * dOut_8.x;
    float sum_16 = _S51 + (*right_8).primal_0.rows[int(1)].y * dOut_8.y;
    *&(((&right_d_result_4)->rows + (int(1)))->y) = (*left_8).primal_0.y * dOut_8.y;
    float sum_17 = sum_16 + (*right_8).primal_0.rows[int(1)].z * dOut_8.z;
    *&(((&right_d_result_4)->rows + (int(1)))->z) = (*left_8).primal_0.y * dOut_8.z;
    *&((&left_d_result_4)->y) = sum_17;
    float _S52 = (*right_8).primal_0.rows[int(2)].x * dOut_8.x;
    *&(((&right_d_result_4)->rows + (int(2)))->x) = (*left_8).primal_0.z * dOut_8.x;
    float sum_18 = _S52 + (*right_8).primal_0.rows[int(2)].y * dOut_8.y;
    *&(((&right_d_result_4)->rows + (int(2)))->y) = (*left_8).primal_0.z * dOut_8.y;
    float sum_19 = sum_18 + (*right_8).primal_0.rows[int(2)].z * dOut_8.z;
    *&(((&right_d_result_4)->rows + (int(2)))->z) = (*left_8).primal_0.z * dOut_8.z;
    *&((&left_d_result_4)->z) = sum_19;
    left_8->primal_0 = (*left_8).primal_0;
    left_8->differential_0 = left_d_result_4;
    right_8->primal_0 = (*right_8).primal_0;
    right_8->differential_0 = right_d_result_4;
    return;
}

inline __device__ float3  mul_7(float3  left_9, Matrix<float, 3, 3>  right_9)
{
    float3  result_13;
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
        int i_7 = int(0);
        float sum_20 = 0.0f;
        for(;;)
        {
            if(i_7 < int(3))
            {
            }
            else
            {
                break;
            }
            float sum_21 = sum_20 + _slang_vector_get_element(left_9, i_7) * _slang_vector_get_element(right_9.rows[i_7], j_1);
            i_7 = i_7 + int(1);
            sum_20 = sum_21;
        }
        *_slang_vector_get_element_ptr(&result_13, j_1) = sum_20;
        j_1 = j_1 + int(1);
    }
    return result_13;
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float2  normalize_1(float2  x_11)
{
    return x_11 / make_float2 (length_0(x_11));
}

inline __device__ void generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_3, bool is_fisheye_2, float4  radial_coeffs_4, float2  tangential_coeffs_4, float2  thin_prism_coeffs_4, float3  * ray_o_0, float3  * ray_d_0)
{
    *ray_o_0 = - mul_7(t_1, R_2);
    CameraDistortion_0 _S53;
    (&_S53)->radial_coeffs_0 = radial_coeffs_4;
    (&_S53)->tangential_coeffs_0 = tangential_coeffs_4;
    (&_S53)->thin_prism_coeffs_0 = thin_prism_coeffs_4;
    float2  _S54 = undistort_point_0(uv_3, &_S53, int(8));
    float3  raydir_1;
    if(is_fisheye_2)
    {
        float theta_2 = length_0(_S54);
        float _S55;
        if(theta_2 < 0.00100000004749745f)
        {
            _S55 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S55 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S56 = make_float3 ((_S54 * make_float2 (_S55)).x, (_S54 * make_float2 (_S55)).y, (F32_cos((theta_2))));
        raydir_1 = _S56;
    }
    else
    {
        raydir_1 = make_float3 (_S54.x, _S54.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_7(raydir_1, R_2));
    return;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S57;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S58, Matrix<float, 3, 3>  _S59)
{
    return mul_7(_S58, _S59);
}

inline __device__ float s_primal_ctx_sin_0(float _S60)
{
    return (F32_sin((_S60)));
}

inline __device__ float s_primal_ctx_cos_0(float _S61)
{
    return (F32_cos((_S61)));
}

inline __device__ void s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_4, bool is_fisheye_3, float4  radial_coeffs_5, float2  tangential_coeffs_5, float2  thin_prism_coeffs_5, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S57 = make_float2 (0.0f);
    float3  _S62 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    CameraDistortion_0 _S63;
    (&_S63)->radial_coeffs_0 = radial_coeffs_5;
    (&_S63)->tangential_coeffs_0 = tangential_coeffs_5;
    (&_S63)->thin_prism_coeffs_0 = thin_prism_coeffs_5;
    float2  _S64 = undistort_point_0(uv_4, &_S63, int(8));
    _s_diff_ctx_0->_S57 = _S64;
    float3  raydir_2;
    if(is_fisheye_3)
    {
        float _S65 = length_0(_S64);
        float _S66;
        if(_S65 < 0.00100000004749745f)
        {
            _S66 = 1.0f - _S65 * _S65 / 6.0f;
        }
        else
        {
            _S66 = s_primal_ctx_sin_0(_S65) / _S65;
        }
        float3  _S67 = make_float3 ((_S64 * make_float2 (_S66)).x, (_S64 * make_float2 (_S66)).y, s_primal_ctx_cos_0(_S65));
        raydir_2 = _S67;
    }
    else
    {
        raydir_2 = make_float3 (_S64.x, _S64.y, 1.0f);
    }
    float3  _S68 = normalize_0(s_primal_ctx_mul_0(raydir_2, dpR_0));
    *dpray_o_0 = _S62;
    *dpray_d_0 = _S68;
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S69, float _S70)
{
    _d_sqrt_0(_S69, _S70);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_4, float _s_dOut_0)
{
    float _S71 = (*dpx_4).primal_0.x;
    float _S72 = (*dpx_4).primal_0.y;
    float _S73 = (*dpx_4).primal_0.z;
    DiffPair_float_0 _S74;
    (&_S74)->primal_0 = _S71 * _S71 + _S72 * _S72 + _S73 * _S73;
    (&_S74)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S74, _s_dOut_0);
    float _S75 = (*dpx_4).primal_0.z * _S74.differential_0;
    float _S76 = _S75 + _S75;
    float _S77 = (*dpx_4).primal_0.y * _S74.differential_0;
    float _S78 = _S77 + _S77;
    float _S79 = (*dpx_4).primal_0.x * _S74.differential_0;
    float _S80 = _S79 + _S79;
    float3  _S81 = make_float3 (0.0f);
    *&((&_S81)->z) = _S76;
    *&((&_S81)->y) = _S78;
    *&((&_S81)->x) = _S80;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S81;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S82, float _S83)
{
    s_bwd_prop_length_impl_0(_S82, _S83);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float3  _s_dOut_1)
{
    float _S84 = length_1((*dpx_5).primal_0);
    float3  _S85 = (*dpx_5).primal_0 * _s_dOut_1;
    float3  _S86 = make_float3 (1.0f / _S84) * _s_dOut_1;
    float _S87 = - ((_S85.x + _S85.y + _S85.z) / (_S84 * _S84));
    float3  _S88 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S89;
    (&_S89)->primal_0 = (*dpx_5).primal_0;
    (&_S89)->differential_0 = _S88;
    s_bwd_length_impl_0(&_S89, _S87);
    float3  _S90 = _S86 + _S89.differential_0;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S90;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S91, float3  _S92)
{
    s_bwd_prop_normalize_impl_0(_S91, _S92);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S93, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S94, float3  _S95)
{
    _d_mul_1(_S93, _S94, _S95);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_5, bool is_fisheye_4, float4  radial_coeffs_6, float2  tangential_coeffs_6, float2  thin_prism_coeffs_6, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S96 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S97 = *dpt_1;
    float3  raydir_3;
    if(is_fisheye_4)
    {
        float _S98 = length_0(_s_diff_ctx_1->_S57);
        float _S99;
        if(_S98 < 0.00100000004749745f)
        {
            _S99 = 1.0f - _S98 * _S98 / 6.0f;
        }
        else
        {
            _S99 = s_primal_ctx_sin_0(_S98) / _S98;
        }
        float3  _S100 = make_float3 ((_s_diff_ctx_1->_S57 * make_float2 (_S99)).x, (_s_diff_ctx_1->_S57 * make_float2 (_S99)).y, s_primal_ctx_cos_0(_S98));
        raydir_3 = _S100;
    }
    else
    {
        raydir_3 = make_float3 (_s_diff_ctx_1->_S57.x, _s_diff_ctx_1->_S57.y, 1.0f);
    }
    float3  _S101 = s_primal_ctx_mul_0(raydir_3, _S96.primal_0);
    float3  _S102 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S103;
    (&_S103)->primal_0 = _S101;
    (&_S103)->differential_0 = _S102;
    s_bwd_normalize_impl_0(&_S103, dpray_d_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S104;
    (&_S104)->primal_0 = raydir_3;
    (&_S104)->differential_0 = _S102;
    Matrix<float, 3, 3>  _S105 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S106;
    (&_S106)->primal_0 = _S96.primal_0;
    (&_S106)->differential_0 = _S105;
    s_bwd_prop_mul_0(&_S104, &_S106, _S103.differential_0);
    float3  _S107 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S108;
    (&_S108)->primal_0 = _S97.primal_0;
    (&_S108)->differential_0 = _S102;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S109;
    (&_S109)->primal_0 = _S96.primal_0;
    (&_S109)->differential_0 = _S105;
    s_bwd_prop_mul_0(&_S108, &_S109, _S107);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S108.differential_0;
    Matrix<float, 3, 3>  _S110 = _S109.differential_0 + _S106.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S110;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S111, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S112, float2  _S113, bool _S114, float4  _S115, float2  _S116, float2  _S117, float3  _S118, float3  _S119)
{
    float3  _S120;
    float3  _S121;
    s_bwd_prop_generate_ray_Intermediates_0 _S122;
    s_primal_ctx_generate_ray_0((*_S111).primal_0, (*_S112).primal_0, _S113, _S114, _S115, _S116, _S117, &_S120, &_S121, &_S122);
    s_bwd_prop_generate_ray_Intermediates_0 _S123 = _S122;
    s_bwd_prop_generate_ray_0(_S111, _S112, _S113, _S114, _S115, _S116, _S117, _S118, _S119, &_S123);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_6, bool is_fisheye_5, float4  radial_coeffs_7, float2  tangential_coeffs_7, float2  thin_prism_coeffs_7, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S124 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S124;
    float3  _S125 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S125;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_6, is_fisheye_5, radial_coeffs_7, tangential_coeffs_7, thin_prism_coeffs_7, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_6, DiffPair_float_0 * dpy_3, float dOut_9)
{
    DiffPair_float_0 _S126 = *dpx_6;
    float _S127;
    if(((*dpx_6).primal_0) < ((*dpy_3).primal_0))
    {
        _S127 = dOut_9;
    }
    else
    {
        if(((*dpx_6).primal_0) > ((*dpy_3).primal_0))
        {
            _S127 = 0.0f;
        }
        else
        {
            _S127 = 0.5f * dOut_9;
        }
    }
    dpx_6->primal_0 = _S126.primal_0;
    dpx_6->differential_0 = _S127;
    DiffPair_float_0 _S128 = *dpy_3;
    if(((*dpy_3).primal_0) < (_S126.primal_0))
    {
        _S127 = dOut_9;
    }
    else
    {
        if(((*dpy_3).primal_0) > ((*dpx_6).primal_0))
        {
            _S127 = 0.0f;
        }
        else
        {
            _S127 = 0.5f * dOut_9;
        }
    }
    dpy_3->primal_0 = _S128.primal_0;
    dpy_3->differential_0 = _S127;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S129 = float(width_0);
    float _S130 = float(height_0);
    float _S131 = 0.30000001192092896f * (0.5f * _S129 / fx_0);
    float _S132 = 0.30000001192092896f * (0.5f * _S130 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_1 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S129 - cx_0) / fx_0 + _S131), ((F32_max((- (cx_0 / fx_0 + _S131)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S130 - cy_0) / fy_0 + _S132), ((F32_max((- (cy_0 / fy_0 + _S132)), (mean3d_0.y * rz_0))))))) * rz2_0);
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

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_7, float dOut_10)
{
    DiffPair_float_0 _S133 = *dpx_7;
    float _S134 = - (*dpy_4).primal_0 / ((*dpx_7).primal_0 * (*dpx_7).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S134;
    float _S135 = _S133.primal_0 / (_S133.primal_0 * _S133.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S135;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S136, float _S137)
{
    return (F32_atan2((_S136), (_S137)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S138, DiffPair_float_0 * _S139, float _S140)
{
    _d_atan2_0(_S138, _S139, _S140);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_8, float _s_dOut_2)
{
    float _S141 = (*dpx_8).primal_0.x;
    float _S142 = (*dpx_8).primal_0.y;
    DiffPair_float_0 _S143;
    (&_S143)->primal_0 = _S141 * _S141 + _S142 * _S142;
    (&_S143)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S143, _s_dOut_2);
    float _S144 = (*dpx_8).primal_0.y * _S143.differential_0;
    float _S145 = _S144 + _S144;
    float _S146 = (*dpx_8).primal_0.x * _S143.differential_0;
    float _S147 = _S146 + _S146;
    float2  _S148 = make_float2 (0.0f);
    *&((&_S148)->y) = _S145;
    *&((&_S148)->x) = _S147;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S148;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S149, float _S150)
{
    s_bwd_prop_length_impl_1(_S149, _S150);
    return;
}

inline __device__ void fisheye_proj_3dgs(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float fx_2, float fy_2, float cx_2, float cy_2, float4  radial_coeffs_8, float2  tangential_coeffs_8, float2  thin_prism_coeffs_8, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    CameraDistortion_0 dist_coeffs_1 = CameraDistortion_x24init_0(radial_coeffs_8, tangential_coeffs_8, thin_prism_coeffs_8);
    float2  _S151 = float2 {mean3d_2.x, mean3d_2.y};
    float r_7 = length_0(_S151);
    float _S152 = mean3d_2.z;
    float theta_3 = (F32_atan2((r_7), (_S152)));
    float k_0;
    if(theta_3 < 0.00100000004749745f)
    {
        k_0 = (1.0f - theta_3 * theta_3 / 3.0f) / _S152;
    }
    else
    {
        k_0 = theta_3 / r_7;
    }
    float2  _S153 = _S151 * make_float2 (k_0);
    float k1_0 = dist_coeffs_1.radial_coeffs_0.x;
    float k2_0 = dist_coeffs_1.radial_coeffs_0.y;
    float k3_0 = dist_coeffs_1.radial_coeffs_0.z;
    float k4_1 = dist_coeffs_1.radial_coeffs_0.w;
    float p1_2 = dist_coeffs_1.tangential_coeffs_0.x;
    float p2_2 = dist_coeffs_1.tangential_coeffs_0.y;
    float sx1_1 = dist_coeffs_1.thin_prism_coeffs_0.x;
    float sy1_1 = dist_coeffs_1.thin_prism_coeffs_0.y;
    float u_2 = _S153.x;
    float v_2 = _S153.y;
    float r2_2 = u_2 * u_2 + v_2 * v_2;
    float _S154 = 2.0f * p1_2;
    float _S155 = 2.0f * p2_2;
    float2  _S156 = _S153 * make_float2 (1.0f + r2_2 * (k1_0 + r2_2 * (k2_0 + r2_2 * (k3_0 + r2_2 * k4_1)))) + make_float2 (_S154 * u_2 * v_2 + p2_2 * (r2_2 + 2.0f * u_2 * u_2) + sx1_1 * r2_2, _S155 * u_2 * v_2 + p1_2 * (r2_2 + 2.0f * v_2 * v_2) + sy1_1 * r2_2);
    *mean2d_2 = make_float2 (fx_2 * _S156.x + cx_2, fy_2 * _S156.y + cy_2);
    Matrix<float, 2, 3>  J_3;
    float2  _S157 = make_float2 (0.0f);
    float2  seed_2 = _S157;
    *&((&seed_2)->x) = 1.0f;
    float2  _S158 = seed_2;
    float _S159 = s_primal_ctx_atan2_0(r_7, _S152);
    bool _S160 = _S159 < 0.00100000004749745f;
    float _S161;
    float _S162;
    float _S163;
    if(_S160)
    {
        float _S164 = 1.0f - _S159 * _S159 / 3.0f;
        float _S165 = _S152 * _S152;
        k_0 = _S164 / _S152;
        _S161 = 0.0f;
        _S162 = _S165;
        _S163 = _S164;
    }
    else
    {
        float _S166 = r_7 * r_7;
        k_0 = _S159 / r_7;
        _S161 = _S166;
        _S162 = 0.0f;
        _S163 = 0.0f;
    }
    float2  _S167 = make_float2 (k_0);
    float2  _S168 = _S151 * make_float2 (k_0);
    float u_3 = _S168.x;
    float v_3 = _S168.y;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float _S169 = k3_0 + r2_3 * k4_1;
    float _S170 = k2_0 + r2_3 * _S169;
    float _S171 = k1_0 + r2_3 * _S170;
    float _S172 = fy_2 * _S158.y;
    float _S173 = fx_2 * _S158.x;
    float2  _S174 = make_float2 (_S173, _S172);
    float2  _S175 = _S168 * _S174;
    float _S176 = p1_2 * _S172;
    float _S177 = p2_2 * _S173;
    float _S178 = _S175.x + _S175.y;
    float _S179 = r2_3 * _S178;
    float _S180 = r2_3 * _S179;
    float _S181 = sy1_1 * _S172 + _S176 + sx1_1 * _S173 + _S177 + _S171 * _S178 + _S170 * _S179 + _S169 * _S180 + k4_1 * (r2_3 * _S180);
    float _S182 = v_3 * _S181;
    float _S183 = u_3 * _S181;
    float2  _S184 = make_float2 (1.0f + r2_3 * _S171) * _S174 + make_float2 (_S155 * (v_3 * _S172) + 2.0f * u_3 * _S177 + 2.0f * (u_3 * _S177) + _S154 * (v_3 * _S173) + _S183 + _S183, 2.0f * v_3 * _S176 + 2.0f * (v_3 * _S176) + _S155 * u_3 * _S172 + _S154 * u_3 * _S173 + _S182 + _S182);
    float2  _S185 = _S151 * _S184;
    float2  _S186 = _S167 * _S184;
    float _S187 = _S185.x + _S185.y;
    if(_S160)
    {
        float _S188 = _S187 / _S162;
        float _S189 = _S163 * - _S188;
        float _S190 = _S159 * (0.3333333432674408f * - (_S152 * _S188));
        k_0 = _S190 + _S190;
        _S161 = _S189;
        _S162 = 0.0f;
    }
    else
    {
        float _S191 = _S187 / _S161;
        float _S192 = _S159 * - _S191;
        k_0 = r_7 * _S191;
        _S161 = 0.0f;
        _S162 = _S192;
    }
    DiffPair_float_0 _S193;
    (&_S193)->primal_0 = r_7;
    (&_S193)->differential_0 = 0.0f;
    DiffPair_float_0 _S194;
    (&_S194)->primal_0 = _S152;
    (&_S194)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S193, &_S194, k_0);
    float _S195 = _S194.differential_0 + _S161;
    float _S196 = _S193.differential_0 + _S162;
    float2  _S197 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S198;
    (&_S198)->primal_0 = _S151;
    (&_S198)->differential_0 = _S197;
    s_bwd_length_impl_1(&_S198, _S196);
    float2  _S199 = _S198.differential_0 + _S186;
    float3  _S200 = make_float3 (_S199.x, _S199.y, _S195);
    J_3[int(0)] = _S200;
    float2  seed_3 = _S157;
    *&((&seed_3)->y) = 1.0f;
    float2  _S201 = seed_3;
    if(_S160)
    {
        float _S202 = 1.0f - _S159 * _S159 / 3.0f;
        float _S203 = _S152 * _S152;
        k_0 = _S202 / _S152;
        _S161 = 0.0f;
        _S162 = _S203;
        _S163 = _S202;
    }
    else
    {
        float _S204 = r_7 * r_7;
        k_0 = _S159 / r_7;
        _S161 = _S204;
        _S162 = 0.0f;
        _S163 = 0.0f;
    }
    float2  _S205 = make_float2 (k_0);
    float2  _S206 = _S151 * make_float2 (k_0);
    float u_4 = _S206.x;
    float v_4 = _S206.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float _S207 = k3_0 + r2_4 * k4_1;
    float _S208 = k2_0 + r2_4 * _S207;
    float _S209 = k1_0 + r2_4 * _S208;
    float _S210 = fy_2 * _S201.y;
    float _S211 = fx_2 * _S201.x;
    float2  _S212 = make_float2 (_S211, _S210);
    float2  _S213 = _S206 * _S212;
    float _S214 = p1_2 * _S210;
    float _S215 = p2_2 * _S211;
    float _S216 = _S213.x + _S213.y;
    float _S217 = r2_4 * _S216;
    float _S218 = r2_4 * _S217;
    float _S219 = sy1_1 * _S210 + _S214 + sx1_1 * _S211 + _S215 + _S209 * _S216 + _S208 * _S217 + _S207 * _S218 + k4_1 * (r2_4 * _S218);
    float _S220 = v_4 * _S219;
    float _S221 = u_4 * _S219;
    float2  _S222 = make_float2 (1.0f + r2_4 * _S209) * _S212 + make_float2 (_S155 * (v_4 * _S210) + 2.0f * u_4 * _S215 + 2.0f * (u_4 * _S215) + _S154 * (v_4 * _S211) + _S221 + _S221, 2.0f * v_4 * _S214 + 2.0f * (v_4 * _S214) + _S155 * u_4 * _S210 + _S154 * u_4 * _S211 + _S220 + _S220);
    float2  _S223 = _S151 * _S222;
    float2  _S224 = _S205 * _S222;
    float _S225 = _S223.x + _S223.y;
    if(_S160)
    {
        float _S226 = _S225 / _S162;
        float _S227 = _S163 * - _S226;
        float _S228 = _S159 * (0.3333333432674408f * - (_S152 * _S226));
        k_0 = _S228 + _S228;
        _S161 = _S227;
        _S162 = 0.0f;
    }
    else
    {
        float _S229 = _S225 / _S161;
        float _S230 = _S159 * - _S229;
        k_0 = r_7 * _S229;
        _S161 = 0.0f;
        _S162 = _S230;
    }
    DiffPair_float_0 _S231;
    (&_S231)->primal_0 = r_7;
    (&_S231)->differential_0 = 0.0f;
    DiffPair_float_0 _S232;
    (&_S232)->primal_0 = _S152;
    (&_S232)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S231, &_S232, k_0);
    float _S233 = _S232.differential_0 + _S161;
    float _S234 = _S231.differential_0 + _S162;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S235;
    (&_S235)->primal_0 = _S151;
    (&_S235)->differential_0 = _S197;
    s_bwd_length_impl_1(&_S235, _S234);
    float2  _S236 = _S235.differential_0 + _S224;
    float3  _S237 = make_float3 (_S236.x, _S236.y, _S233);
    J_3[int(1)] = _S237;
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
    float _S238 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S238;
    float _S239 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S239;
    float det_blur_0 = _S238 * _S239 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_9, float dOut_11)
{
    float _S240 = (F32_exp(((*dpx_9).primal_0))) * dOut_11;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S240;
    return;
}

inline __device__ float3  exp_0(float3  x_12)
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
        *_slang_vector_get_element_ptr(&result_14, i_8) = (F32_exp((_slang_vector_get_element(x_12, i_8))));
        i_8 = i_8 + int(1);
    }
    return result_14;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_10, float3  dOut_12)
{
    float3  _S241 = exp_0((*dpx_10).primal_0) * dOut_12;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S241;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_11, float dOut_13)
{
    float _S242 = 1.0f / (*dpx_11).primal_0 * dOut_13;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S242;
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_12, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_5, float3  dOut_14)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_12).primal_0.x;
    (&left_dp_0)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_0;
    (&right_dp_0)->primal_0 = (*dpy_5).primal_0.x;
    (&right_dp_0)->differential_0 = 0.0f;
    _d_max_0(&left_dp_0, &right_dp_0, dOut_14.x);
    float3  left_d_result_5;
    *&((&left_d_result_5)->x) = left_dp_0.differential_0;
    float3  right_d_result_5;
    *&((&right_d_result_5)->x) = right_dp_0.differential_0;
    DiffPair_float_0 left_dp_1;
    (&left_dp_1)->primal_0 = (*dpx_12).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_5).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_14.y);
    *&((&left_d_result_5)->y) = left_dp_1.differential_0;
    *&((&right_d_result_5)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_12).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_5).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_14.z);
    *&((&left_d_result_5)->z) = left_dp_2.differential_0;
    *&((&right_d_result_5)->z) = right_dp_2.differential_0;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = left_d_result_5;
    dpy_5->primal_0 = (*dpy_5).primal_0;
    dpy_5->differential_0 = right_d_result_5;
    return;
}

inline __device__ float3  max_0(float3  x_13, float3  y_2)
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
        *_slang_vector_get_element_ptr(&result_15, i_9) = (F32_max((_slang_vector_get_element(x_13, i_9)), (_slang_vector_get_element(y_2, i_9))));
        i_9 = i_9 + int(1);
    }
    return result_15;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_4, float fy_4, float cx_4, float cy_4, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_4, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_4, mean_0) + t_3;
        float _S243 = mean_c_0.z;
        bool _S244;
        if(_S243 < near_plane_0)
        {
            _S244 = true;
        }
        else
        {
            _S244 = _S243 > far_plane_0;
        }
        if(_S244)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S245 = exp_0(scale_2);
        float x_14 = quat_3.y;
        float x2_3 = x_14 * x_14;
        float y2_3 = quat_3.z * quat_3.z;
        float z2_3 = quat_3.w * quat_3.w;
        float xy_3 = quat_3.y * quat_3.z;
        float xz_3 = quat_3.y * quat_3.w;
        float yz_3 = quat_3.z * quat_3.w;
        float wx_3 = quat_3.x * quat_3.y;
        float wy_3 = quat_3.x * quat_3.z;
        float wz_3 = quat_3.x * quat_3.w;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S245.x, 0.0f, 0.0f, 0.0f, _S245.y, 0.0f, 0.0f, 0.0f, _S245.z));
        Matrix<float, 3, 3>  _S246 = transpose_0(R_4);
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_4, mul_4(M_2, transpose_0(M_2))), _S246);
        Matrix<float, 2, 2>  covar2d_0;
        float _S247 = float(image_width_0);
        float _S248 = float(image_height_0);
        float _S249 = 0.30000001192092896f * (0.5f * _S247 / fx_4);
        float _S250 = 0.30000001192092896f * (0.5f * _S248 / fy_4);
        float rz_2 = 1.0f / mean_c_0.z;
        float rz2_2 = rz_2 * rz_2;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_4 * rz_2, 0.0f, - fx_4 * (mean_c_0.z * (F32_min(((_S247 - cx_4) / fx_4 + _S249), ((F32_max((- (cx_4 / fx_4 + _S249)), (mean_c_0.x * rz_2))))))) * rz2_2, 0.0f, fy_4 * rz_2, - fy_4 * (mean_c_0.z * (F32_min(((_S248 - cy_4) / fy_4 + _S250), ((F32_max((- (cy_4 / fy_4 + _S250)), (mean_c_0.y * rz_2))))))) * rz2_2);
        covar2d_0 = mul_6(mul_5(J_5, covar_c_0), transpose_1(J_5));
        *mean2d_4 = make_float2 (fx_4 * mean_c_0.x * rz_2 + cx_4, fy_4 * mean_c_0.y * rz_2 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S251 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S251;
        float _S252 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S252;
        float det_blur_1 = _S251 * _S252 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S253 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S244 = true;
        }
        else
        {
            _S244 = xmin_0 >= _S247;
        }
        if(_S244)
        {
            _S244 = true;
        }
        else
        {
            _S244 = ymax_0 <= 0.0f;
        }
        if(_S244)
        {
            _S244 = true;
        }
        else
        {
            _S244 = ymin_0 >= _S248;
        }
        if(_S244)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = 0.5f * (F32_log((dot_0(mean_c_0, mean_c_0) + 9.99999997475242708e-07f)));
        *conic_0 = make_float3 (_S253.rows[int(0)].x, _S253.rows[int(0)].y, _S253.rows[int(1)].y);
        float3  _S254 = mean_0 - - mul_0(_S246, t_3);
        float3  _S255 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S255;
        float _S256 = _S254.x;
        float _S257 = _S254.y;
        float _S258 = _S254.z;
        float norm_0 = (F32_sqrt((_S256 * _S256 + _S257 * _S257 + _S258 * _S258)));
        float x_15 = _S256 / norm_0;
        float y_3 = _S257 / norm_0;
        float z_0 = _S258 / norm_0;
        float3  _S259 = _S255 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_15) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S259;
        float z2_4 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_15 * x_15 - y_3 * y_3;
        float fS1_0 = 2.0f * x_15 * y_3;
        float3  _S260 = _S259 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_4 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_15) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S260;
        float fTmp0C_0 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgb_0 = max_0(_S260 + (make_float3 (-0.59004360437393188f * (x_15 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_4 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_15) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_5, float fy_5, float cx_5, float cy_5, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_5, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_5, mean_1) + t_4;
        float _S261 = length_1(mean_c_1);
        bool _S262;
        if(_S261 < near_plane_1)
        {
            _S262 = true;
        }
        else
        {
            _S262 = _S261 > far_plane_1;
        }
        if(_S262)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S263 = exp_0(scale_3);
        float x_16 = quat_4.y;
        float x2_4 = x_16 * x_16;
        float y2_4 = quat_4.z * quat_4.z;
        float z2_5 = quat_4.w * quat_4.w;
        float xy_4 = quat_4.y * quat_4.z;
        float xz_4 = quat_4.y * quat_4.w;
        float yz_4 = quat_4.z * quat_4.w;
        float wx_4 = quat_4.x * quat_4.y;
        float wy_4 = quat_4.x * quat_4.z;
        float wz_4 = quat_4.x * quat_4.w;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_5), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_5), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S263.x, 0.0f, 0.0f, 0.0f, _S263.y, 0.0f, 0.0f, 0.0f, _S263.z));
        Matrix<float, 3, 3>  _S264 = transpose_0(R_5);
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_5, mul_4(M_3, transpose_0(M_3))), _S264);
        Matrix<float, 2, 2>  covar2d_1;
        fisheye_proj_3dgs(mean_c_1, covar_c_1, fx_5, fy_5, cx_5, cy_5, radial_coeffs_10, tangential_coeffs_10, thin_prism_coeffs_10, &covar2d_1, mean2d_5);
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S265 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S265;
        float _S266 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S266;
        float det_blur_2 = _S265 * _S266 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S267 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S262 = true;
        }
        else
        {
            _S262 = xmin_1 >= float(image_width_1);
        }
        if(_S262)
        {
            _S262 = true;
        }
        else
        {
            _S262 = ymax_1 <= 0.0f;
        }
        if(_S262)
        {
            _S262 = true;
        }
        else
        {
            _S262 = ymin_1 >= float(image_height_1);
        }
        if(_S262)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = 0.5f * (F32_log((dot_0(mean_c_1, mean_c_1) + 9.99999997475242708e-07f)));
        *conic_1 = make_float3 (_S267.rows[int(0)].x, _S267.rows[int(0)].y, _S267.rows[int(1)].y);
        float3  _S268 = mean_1 - - mul_0(_S264, t_4);
        float3  _S269 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S269;
        float _S270 = _S268.x;
        float _S271 = _S268.y;
        float _S272 = _S268.z;
        float norm_1 = (F32_sqrt((_S270 * _S270 + _S271 * _S271 + _S272 * _S272)));
        float x_17 = _S270 / norm_1;
        float y_4 = _S271 / norm_1;
        float z_1 = _S272 / norm_1;
        float3  _S273 = _S269 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_1)[int(1)] + make_float3 (z_1) * (*sh_coeffs_1)[int(2)] - make_float3 (x_17) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S273;
        float z2_6 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_17 * x_17 - y_4 * y_4;
        float fS1_1 = 2.0f * x_17 * y_4;
        float3  _S274 = _S273 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_4) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_6 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_17) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S274;
        float fTmp0C_1 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgb_1 = max_0(_S274 + (make_float3 (-0.59004360437393188f * (x_17 * fS1_1 + y_4 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_4) * (*sh_coeffs_1)[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_6 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_17) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_1 - y_4 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_6, float fy_6, float cx_6, float cy_6, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_6, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_6, mean_2) + t_5;
        float _S275 = mean_c_2.z;
        bool _S276;
        if(_S275 < near_plane_2)
        {
            _S276 = true;
        }
        else
        {
            _S276 = _S275 > far_plane_2;
        }
        if(_S276)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S277 = exp_0(scale_4);
        float x_18 = quat_5.y;
        float x2_5 = x_18 * x_18;
        float y2_5 = quat_5.z * quat_5.z;
        float z2_7 = quat_5.w * quat_5.w;
        float xy_5 = quat_5.y * quat_5.z;
        float xz_5 = quat_5.y * quat_5.w;
        float yz_5 = quat_5.z * quat_5.w;
        float wx_5 = quat_5.x * quat_5.y;
        float wy_5 = quat_5.x * quat_5.z;
        float wz_5 = quat_5.x * quat_5.w;
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_7), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_7), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S277.x, 0.0f, 0.0f, 0.0f, _S277.y, 0.0f, 0.0f, 0.0f, _S277.z));
        Matrix<float, 3, 3>  _S278 = transpose_0(R_6);
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_6, mul_4(M_4, transpose_0(M_4))), _S278);
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_6, covar_c_2), transpose_1(J_6));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S279 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S279;
        float _S280 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S280;
        float det_blur_3 = _S279 * _S280 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S281 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S276 = true;
        }
        else
        {
            _S276 = xmin_2 >= float(image_width_2);
        }
        if(_S276)
        {
            _S276 = true;
        }
        else
        {
            _S276 = ymax_2 <= 0.0f;
        }
        if(_S276)
        {
            _S276 = true;
        }
        else
        {
            _S276 = ymin_2 >= float(image_height_2);
        }
        if(_S276)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = 0.5f * (F32_log((dot_0(mean_c_2, mean_c_2) + 9.99999997475242708e-07f)));
        *conic_2 = make_float3 (_S281.rows[int(0)].x, _S281.rows[int(0)].y, _S281.rows[int(1)].y);
        float3  _S282 = mean_2 - - mul_0(_S278, t_5);
        float3  _S283 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S283;
        float _S284 = _S282.x;
        float _S285 = _S282.y;
        float _S286 = _S282.z;
        float norm_2 = (F32_sqrt((_S284 * _S284 + _S285 * _S285 + _S286 * _S286)));
        float x_19 = _S284 / norm_2;
        float y_5 = _S285 / norm_2;
        float z_2 = _S286 / norm_2;
        float3  _S287 = _S283 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_2)[int(1)] + make_float3 (z_2) * (*sh_coeffs_2)[int(2)] - make_float3 (x_19) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S287;
        float z2_8 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_19 * x_19 - y_5 * y_5;
        float fS1_2 = 2.0f * x_19 * y_5;
        float3  _S288 = _S287 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_5) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_8 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_19) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S288;
        float fTmp0C_2 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        *rgb_2 = max_0(_S288 + (make_float3 (-0.59004360437393188f * (x_19 * fS1_2 + y_5 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_5) * (*sh_coeffs_2)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_8 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_19) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_19 * fC1_2 - y_5 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_persp(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_7, float fy_7, float cx_7, float cy_7, float4  radial_coeffs_12, float2  tangential_coeffs_12, float2  thin_prism_coeffs_12, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_7, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_0(R_7, mean_3) + t_6;
        float _S289 = mean_c_3.z;
        bool _S290;
        if(_S289 < near_plane_3)
        {
            _S290 = true;
        }
        else
        {
            _S290 = _S289 > far_plane_3;
        }
        if(_S290)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S291 = exp_0(scale_5);
        float x_20 = quat_6.y;
        float x2_6 = x_20 * x_20;
        float y2_6 = quat_6.z * quat_6.z;
        float z2_9 = quat_6.w * quat_6.w;
        float xy_6 = quat_6.y * quat_6.z;
        float xz_6 = quat_6.y * quat_6.w;
        float yz_6 = quat_6.z * quat_6.w;
        float wx_6 = quat_6.x * quat_6.y;
        float wy_6 = quat_6.x * quat_6.z;
        float wz_6 = quat_6.x * quat_6.w;
        Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_9), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_9), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S291.x, 0.0f, 0.0f, 0.0f, _S291.y, 0.0f, 0.0f, 0.0f, _S291.z));
        Matrix<float, 3, 3>  _S292 = transpose_0(R_7);
        Matrix<float, 3, 3>  covar_c_3 = mul_4(mul_4(R_7, mul_4(M_5, transpose_0(M_5))), _S292);
        Matrix<float, 2, 2>  covar2d_3;
        float _S293 = float(image_width_3);
        float _S294 = float(image_height_3);
        float _S295 = 0.30000001192092896f * (0.5f * _S293 / fx_7);
        float _S296 = 0.30000001192092896f * (0.5f * _S294 / fy_7);
        float rz_3 = 1.0f / mean_c_3.z;
        float rz2_3 = rz_3 * rz_3;
        Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, - fx_7 * (mean_c_3.z * (F32_min(((_S293 - cx_7) / fx_7 + _S295), ((F32_max((- (cx_7 / fx_7 + _S295)), (mean_c_3.x * rz_3))))))) * rz2_3, 0.0f, fy_7 * rz_3, - fy_7 * (mean_c_3.z * (F32_min(((_S294 - cy_7) / fy_7 + _S296), ((F32_max((- (cy_7 / fy_7 + _S296)), (mean_c_3.y * rz_3))))))) * rz2_3);
        covar2d_3 = mul_6(mul_5(J_7, covar_c_3), transpose_1(J_7));
        *mean2d_7 = make_float2 (fx_7 * mean_c_3.x * rz_3 + cx_7, fy_7 * mean_c_3.y * rz_3 + cy_7);
        float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S297 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S297;
        float _S298 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S298;
        float det_blur_4 = _S297 * _S298 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float compensation_4 = (F32_sqrt(((F32_max((0.0f), (det_orig_4 / det_blur_4))))));
        if(det_blur_4 <= 0.0f)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *opacity_3 = 1.0f / (1.0f + (F32_exp((- in_opacity_3))));
        if(antialiased_3)
        {
            *opacity_3 = *opacity_3 * compensation_4;
        }
        if((*opacity_3) < 0.00392156885936856f)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float extend_3 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_3 / 0.00392156885936856f)))))))));
        float radius_x_3 = extend_3 * (F32_sqrt((covar2d_3[int(0)].x)));
        float radius_y_3 = extend_3 * (F32_sqrt((covar2d_3[int(1)].y)));
        float xmin_3 = (F32_floor(((*mean2d_7).x - radius_x_3)));
        float xmax_3 = (F32_ceil(((*mean2d_7).x + radius_x_3)));
        float ymin_3 = (F32_floor(((*mean2d_7).y - radius_y_3)));
        float ymax_3 = (F32_ceil(((*mean2d_7).y + radius_y_3)));
        if(xmax_3 <= 0.0f)
        {
            _S290 = true;
        }
        else
        {
            _S290 = xmin_3 >= _S293;
        }
        if(_S290)
        {
            _S290 = true;
        }
        else
        {
            _S290 = ymax_3 <= 0.0f;
        }
        if(_S290)
        {
            _S290 = true;
        }
        else
        {
            _S290 = ymin_3 >= _S294;
        }
        if(_S290)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_3 = make_int4 (int(xmin_3), int(ymin_3), int(xmax_3), int(ymax_3));
        *depth_3 = 0.5f * (F32_log((dot_0(mean_c_3, mean_c_3) + 9.99999997475242708e-07f)));
        *conic_3 = exp_0(- scale_5);
        float3  _S299 = mean_3 - - mul_0(_S292, t_6);
        float3  _S300 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
        *rgb_3 = _S300;
        float _S301 = _S299.x;
        float _S302 = _S299.y;
        float _S303 = _S299.z;
        float norm_3 = (F32_sqrt((_S301 * _S301 + _S302 * _S302 + _S303 * _S303)));
        float x_21 = _S301 / norm_3;
        float y_6 = _S302 / norm_3;
        float z_3 = _S303 / norm_3;
        float3  _S304 = _S300 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_3)[int(1)] + make_float3 (z_3) * (*sh_coeffs_3)[int(2)] - make_float3 (x_21) * (*sh_coeffs_3)[int(3)]);
        *rgb_3 = _S304;
        float z2_10 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_21 * x_21 - y_6 * y_6;
        float fS1_3 = 2.0f * x_21 * y_6;
        float3  _S305 = _S304 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_6) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_10 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_21) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
        *rgb_3 = _S305;
        float fTmp0C_3 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        *rgb_3 = max_0(_S305 + (make_float3 (-0.59004360437393188f * (x_21 * fS1_3 + y_6 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_6) * (*sh_coeffs_3)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_10 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_21) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_21 * fC1_3 - y_6 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_fisheye(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_8, float fy_8, float cx_8, float cy_8, float4  radial_coeffs_13, float2  tangential_coeffs_13, float2  thin_prism_coeffs_13, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_8, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    for(;;)
    {
        float3  mean_c_4 = mul_0(R_8, mean_4) + t_7;
        float _S306 = length_1(mean_c_4);
        bool _S307;
        if(_S306 < near_plane_4)
        {
            _S307 = true;
        }
        else
        {
            _S307 = _S306 > far_plane_4;
        }
        if(_S307)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S308 = exp_0(scale_6);
        float x_22 = quat_7.y;
        float x2_7 = x_22 * x_22;
        float y2_7 = quat_7.z * quat_7.z;
        float z2_11 = quat_7.w * quat_7.w;
        float xy_7 = quat_7.y * quat_7.z;
        float xz_7 = quat_7.y * quat_7.w;
        float yz_7 = quat_7.z * quat_7.w;
        float wx_7 = quat_7.x * quat_7.y;
        float wy_7 = quat_7.x * quat_7.z;
        float wz_7 = quat_7.x * quat_7.w;
        Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_11), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_11), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S308.x, 0.0f, 0.0f, 0.0f, _S308.y, 0.0f, 0.0f, 0.0f, _S308.z));
        Matrix<float, 3, 3>  _S309 = transpose_0(R_8);
        Matrix<float, 3, 3>  covar_c_4 = mul_4(mul_4(R_8, mul_4(M_6, transpose_0(M_6))), _S309);
        Matrix<float, 2, 2>  covar2d_4;
        fisheye_proj_3dgs(mean_c_4, covar_c_4, fx_8, fy_8, cx_8, cy_8, radial_coeffs_13, tangential_coeffs_13, thin_prism_coeffs_13, &covar2d_4, mean2d_8);
        float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S310 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S310;
        float _S311 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S311;
        float det_blur_5 = _S310 * _S311 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float compensation_5 = (F32_sqrt(((F32_max((0.0f), (det_orig_5 / det_blur_5))))));
        if(det_blur_5 <= 0.0f)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *opacity_4 = 1.0f / (1.0f + (F32_exp((- in_opacity_4))));
        if(antialiased_4)
        {
            *opacity_4 = *opacity_4 * compensation_5;
        }
        if((*opacity_4) < 0.00392156885936856f)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float extend_4 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_4 / 0.00392156885936856f)))))))));
        float radius_x_4 = extend_4 * (F32_sqrt((covar2d_4[int(0)].x)));
        float radius_y_4 = extend_4 * (F32_sqrt((covar2d_4[int(1)].y)));
        float xmin_4 = (F32_floor(((*mean2d_8).x - radius_x_4)));
        float xmax_4 = (F32_ceil(((*mean2d_8).x + radius_x_4)));
        float ymin_4 = (F32_floor(((*mean2d_8).y - radius_y_4)));
        float ymax_4 = (F32_ceil(((*mean2d_8).y + radius_y_4)));
        if(xmax_4 <= 0.0f)
        {
            _S307 = true;
        }
        else
        {
            _S307 = xmin_4 >= float(image_width_4);
        }
        if(_S307)
        {
            _S307 = true;
        }
        else
        {
            _S307 = ymax_4 <= 0.0f;
        }
        if(_S307)
        {
            _S307 = true;
        }
        else
        {
            _S307 = ymin_4 >= float(image_height_4);
        }
        if(_S307)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_4 = make_int4 (int(xmin_4), int(ymin_4), int(xmax_4), int(ymax_4));
        *depth_4 = 0.5f * (F32_log((dot_0(mean_c_4, mean_c_4) + 9.99999997475242708e-07f)));
        *conic_4 = exp_0(- scale_6);
        float3  _S312 = mean_4 - - mul_0(_S309, t_7);
        float3  _S313 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
        *rgb_4 = _S313;
        float _S314 = _S312.x;
        float _S315 = _S312.y;
        float _S316 = _S312.z;
        float norm_4 = (F32_sqrt((_S314 * _S314 + _S315 * _S315 + _S316 * _S316)));
        float x_23 = _S314 / norm_4;
        float y_7 = _S315 / norm_4;
        float z_4 = _S316 / norm_4;
        float3  _S317 = _S313 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_4)[int(1)] + make_float3 (z_4) * (*sh_coeffs_4)[int(2)] - make_float3 (x_23) * (*sh_coeffs_4)[int(3)]);
        *rgb_4 = _S317;
        float z2_12 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_23 * x_23 - y_7 * y_7;
        float fS1_4 = 2.0f * x_23 * y_7;
        float3  _S318 = _S317 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_7) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_23) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
        *rgb_4 = _S318;
        float fTmp0C_4 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        *rgb_4 = max_0(_S318 + (make_float3 (-0.59004360437393188f * (x_23 * fS1_4 + y_7 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_7) * (*sh_coeffs_4)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_23) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_23 * fC1_4 - y_7 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_9, float fy_9, float cx_9, float cy_9, float4  radial_coeffs_14, float2  tangential_coeffs_14, float2  thin_prism_coeffs_14, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_9, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_9, mean_5) + t_8;
    float3  _S319 = exp_0(scale_7);
    float x_24 = quat_8.y;
    float x2_8 = x_24 * x_24;
    float y2_8 = quat_8.z * quat_8.z;
    float z2_13 = quat_8.w * quat_8.w;
    float xy_8 = quat_8.y * quat_8.z;
    float xz_8 = quat_8.y * quat_8.w;
    float yz_8 = quat_8.z * quat_8.w;
    float wx_8 = quat_8.x * quat_8.y;
    float wy_8 = quat_8.x * quat_8.z;
    float wz_8 = quat_8.x * quat_8.w;
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_13), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_13), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S319.x, 0.0f, 0.0f, 0.0f, _S319.y, 0.0f, 0.0f, 0.0f, _S319.z));
    Matrix<float, 3, 3>  _S320 = transpose_0(R_9);
    float _S321 = float(image_width_5);
    float _S322 = float(image_height_5);
    float _S323 = 0.30000001192092896f * (0.5f * _S321 / fx_9);
    float _S324 = 0.30000001192092896f * (0.5f * _S322 / fy_9);
    float rz_4 = 1.0f / mean_c_5.z;
    float rz2_4 = rz_4 * rz_4;
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_9 * rz_4, 0.0f, - fx_9 * (mean_c_5.z * (F32_min(((_S321 - cx_9) / fx_9 + _S323), ((F32_max((- (cx_9 / fx_9 + _S323)), (mean_c_5.x * rz_4))))))) * rz2_4, 0.0f, fy_9 * rz_4, - fy_9 * (mean_c_5.z * (F32_min(((_S322 - cy_9) / fy_9 + _S324), ((F32_max((- (cy_9 / fy_9 + _S324)), (mean_c_5.y * rz_4))))))) * rz2_4);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_8, mul_4(mul_4(R_9, mul_4(M_7, transpose_0(M_7))), _S320)), transpose_1(J_8));
    *mean2d_9 = make_float2 (fx_9 * mean_c_5.x * rz_4 + cx_9, fy_9 * mean_c_5.y * rz_4 + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S325 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S325;
    float _S326 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S326;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S325 * _S326 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S327 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_4, - covar2d_5.rows[int(0)].y * invdet_4, - covar2d_5.rows[int(1)].x * invdet_4, covar2d_5.rows[int(0)].x * invdet_4);
    *opacity_5 = 1.0f / (1.0f + (F32_exp((- in_opacity_5))));
    if(antialiased_5)
    {
        *opacity_5 = *opacity_5 * compensation_6;
    }
    float extend_5 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_5 / 0.00392156885936856f)))))))));
    float radius_x_5 = extend_5 * (F32_sqrt((covar2d_5[int(0)].x)));
    float radius_y_5 = extend_5 * (F32_sqrt((covar2d_5[int(1)].y)));
    *aabb_xyxy_5 = make_int4 (int((F32_floor(((*mean2d_9).x - radius_x_5)))), int((F32_floor(((*mean2d_9).y - radius_y_5)))), int((F32_ceil(((*mean2d_9).x + radius_x_5)))), int((F32_ceil(((*mean2d_9).y + radius_y_5)))));
    *depth_5 = 0.5f * (F32_log((dot_0(mean_c_5, mean_c_5) + 9.99999997475242708e-07f)));
    *conic_5 = make_float3 (_S327.rows[int(0)].x, _S327.rows[int(0)].y, _S327.rows[int(1)].y);
    float3  _S328 = mean_5 - - mul_0(_S320, t_8);
    float3  _S329 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S329;
    float _S330 = _S328.x;
    float _S331 = _S328.y;
    float _S332 = _S328.z;
    float norm_5 = (F32_sqrt((_S330 * _S330 + _S331 * _S331 + _S332 * _S332)));
    float x_25 = _S330 / norm_5;
    float y_8 = _S331 / norm_5;
    float z_5 = _S332 / norm_5;
    float3  _S333 = _S329 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_5)[int(1)] + make_float3 (z_5) * (*sh_coeffs_5)[int(2)] - make_float3 (x_25) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S333;
    float z2_14 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_25 * x_25 - y_8 * y_8;
    float fS1_5 = 2.0f * x_25 * y_8;
    float3  _S334 = _S333 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_8) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_14 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_25) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S334;
    float fTmp0C_5 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    *rgb_5 = max_0(_S334 + (make_float3 (-0.59004360437393188f * (x_25 * fS1_5 + y_8 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_8) * (*sh_coeffs_5)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_14 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_25) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_25 * fC1_5 - y_8 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_10, float fy_10, float cx_10, float cy_10, float4  radial_coeffs_15, float2  tangential_coeffs_15, float2  thin_prism_coeffs_15, uint image_width_6, uint image_height_6, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * mean2d_10, float * depth_6, float3  * conic_6, float * opacity_6, float3  * rgb_6)
{
    float3  mean_c_6 = mul_0(R_10, mean_6) + t_9;
    float3  _S335 = exp_0(scale_8);
    float x_26 = quat_9.y;
    float x2_9 = x_26 * x_26;
    float y2_9 = quat_9.z * quat_9.z;
    float z2_15 = quat_9.w * quat_9.w;
    float xy_9 = quat_9.y * quat_9.z;
    float xz_9 = quat_9.y * quat_9.w;
    float yz_9 = quat_9.z * quat_9.w;
    float wx_9 = quat_9.x * quat_9.y;
    float wy_9 = quat_9.x * quat_9.z;
    float wz_9 = quat_9.x * quat_9.w;
    Matrix<float, 3, 3>  M_8 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S335.x, 0.0f, 0.0f, 0.0f, _S335.y, 0.0f, 0.0f, 0.0f, _S335.z));
    Matrix<float, 3, 3>  _S336 = transpose_0(R_10);
    Matrix<float, 2, 2>  covar2d_6;
    fisheye_proj_3dgs(mean_c_6, mul_4(mul_4(R_10, mul_4(M_8, transpose_0(M_8))), _S336), fx_10, fy_10, cx_10, cy_10, radial_coeffs_15, tangential_coeffs_15, thin_prism_coeffs_15, &covar2d_6, mean2d_10);
    float det_orig_7 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
    float _S337 = *&(((&covar2d_6)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(0)))->x) = _S337;
    float _S338 = *&(((&covar2d_6)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(1)))->y) = _S338;
    float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / (_S337 * _S338 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
    Matrix<float, 2, 2>  _S339 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_5, - covar2d_6.rows[int(0)].y * invdet_5, - covar2d_6.rows[int(1)].x * invdet_5, covar2d_6.rows[int(0)].x * invdet_5);
    *opacity_6 = 1.0f / (1.0f + (F32_exp((- in_opacity_6))));
    if(antialiased_6)
    {
        *opacity_6 = *opacity_6 * compensation_7;
    }
    float extend_6 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_6 / 0.00392156885936856f)))))))));
    float radius_x_6 = extend_6 * (F32_sqrt((covar2d_6[int(0)].x)));
    float radius_y_6 = extend_6 * (F32_sqrt((covar2d_6[int(1)].y)));
    *aabb_xyxy_6 = make_int4 (int((F32_floor(((*mean2d_10).x - radius_x_6)))), int((F32_floor(((*mean2d_10).y - radius_y_6)))), int((F32_ceil(((*mean2d_10).x + radius_x_6)))), int((F32_ceil(((*mean2d_10).y + radius_y_6)))));
    *depth_6 = 0.5f * (F32_log((dot_0(mean_c_6, mean_c_6) + 9.99999997475242708e-07f)));
    *conic_6 = make_float3 (_S339.rows[int(0)].x, _S339.rows[int(0)].y, _S339.rows[int(1)].y);
    float3  _S340 = mean_6 - - mul_0(_S336, t_9);
    float3  _S341 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)];
    *rgb_6 = _S341;
    float _S342 = _S340.x;
    float _S343 = _S340.y;
    float _S344 = _S340.z;
    float norm_6 = (F32_sqrt((_S342 * _S342 + _S343 * _S343 + _S344 * _S344)));
    float x_27 = _S342 / norm_6;
    float y_9 = _S343 / norm_6;
    float z_6 = _S344 / norm_6;
    float3  _S345 = _S341 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_6)[int(1)] + make_float3 (z_6) * (*sh_coeffs_6)[int(2)] - make_float3 (x_27) * (*sh_coeffs_6)[int(3)]);
    *rgb_6 = _S345;
    float z2_16 = z_6 * z_6;
    float fTmp0B_6 = -1.09254848957061768f * z_6;
    float fC1_6 = x_27 * x_27 - y_9 * y_9;
    float fS1_6 = 2.0f * x_27 * y_9;
    float3  _S346 = _S345 + (make_float3 (0.54627424478530884f * fS1_6) * (*sh_coeffs_6)[int(4)] + make_float3 (fTmp0B_6 * y_9) * (*sh_coeffs_6)[int(5)] + make_float3 (0.94617468118667603f * z2_16 - 0.31539157032966614f) * (*sh_coeffs_6)[int(6)] + make_float3 (fTmp0B_6 * x_27) * (*sh_coeffs_6)[int(7)] + make_float3 (0.54627424478530884f * fC1_6) * (*sh_coeffs_6)[int(8)]);
    *rgb_6 = _S346;
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_6;
    *rgb_6 = max_0(_S346 + (make_float3 (-0.59004360437393188f * (x_27 * fS1_6 + y_9 * fC1_6)) * (*sh_coeffs_6)[int(9)] + make_float3 (fTmp1B_6 * fS1_6) * (*sh_coeffs_6)[int(10)] + make_float3 (fTmp0C_6 * y_9) * (*sh_coeffs_6)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_16 - 1.11952900886535645f)) * (*sh_coeffs_6)[int(12)] + make_float3 (fTmp0C_6 * x_27) * (*sh_coeffs_6)[int(13)] + make_float3 (fTmp1B_6 * fC1_6) * (*sh_coeffs_6)[int(14)] + make_float3 (-0.59004360437393188f * (x_27 * fC1_6 - y_9 * fS1_6)) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_11, float fy_11, float cx_11, float cy_11, float4  radial_coeffs_16, float2  tangential_coeffs_16, float2  thin_prism_coeffs_16, uint image_width_7, uint image_height_7, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * mean2d_11, float * depth_7, float3  * conic_7, float * opacity_7, float3  * rgb_7)
{
    float3  mean_c_7 = mul_0(R_11, mean_7) + t_10;
    float3  _S347 = exp_0(scale_9);
    float x_28 = quat_10.y;
    float x2_10 = x_28 * x_28;
    float y2_10 = quat_10.z * quat_10.z;
    float z2_17 = quat_10.w * quat_10.w;
    float xy_10 = quat_10.y * quat_10.z;
    float xz_10 = quat_10.y * quat_10.w;
    float yz_10 = quat_10.z * quat_10.w;
    float wx_10 = quat_10.x * quat_10.y;
    float wy_10 = quat_10.x * quat_10.z;
    float wz_10 = quat_10.x * quat_10.w;
    Matrix<float, 3, 3>  M_9 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S347.x, 0.0f, 0.0f, 0.0f, _S347.y, 0.0f, 0.0f, 0.0f, _S347.z));
    Matrix<float, 3, 3>  _S348 = transpose_0(R_11);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_11, 0.0f, 0.0f, 0.0f, fy_11, 0.0f);
    Matrix<float, 2, 2>  covar2d_7 = mul_6(mul_5(J_9, mul_4(mul_4(R_11, mul_4(M_9, transpose_0(M_9))), _S348)), transpose_1(J_9));
    *mean2d_11 = make_float2 (fx_11 * mean_c_7.x + cx_11, fy_11 * mean_c_7.y + cy_11);
    float det_orig_8 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
    float _S349 = *&(((&covar2d_7)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(0)))->x) = _S349;
    float _S350 = *&(((&covar2d_7)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(1)))->y) = _S350;
    float compensation_8 = (F32_sqrt(((F32_max((0.0f), (det_orig_8 / (_S349 * _S350 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
    Matrix<float, 2, 2>  _S351 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_6, - covar2d_7.rows[int(0)].y * invdet_6, - covar2d_7.rows[int(1)].x * invdet_6, covar2d_7.rows[int(0)].x * invdet_6);
    *opacity_7 = 1.0f / (1.0f + (F32_exp((- in_opacity_7))));
    if(antialiased_7)
    {
        *opacity_7 = *opacity_7 * compensation_8;
    }
    float extend_7 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_7 / 0.00392156885936856f)))))))));
    float radius_x_7 = extend_7 * (F32_sqrt((covar2d_7[int(0)].x)));
    float radius_y_7 = extend_7 * (F32_sqrt((covar2d_7[int(1)].y)));
    *aabb_xyxy_7 = make_int4 (int((F32_floor(((*mean2d_11).x - radius_x_7)))), int((F32_floor(((*mean2d_11).y - radius_y_7)))), int((F32_ceil(((*mean2d_11).x + radius_x_7)))), int((F32_ceil(((*mean2d_11).y + radius_y_7)))));
    *depth_7 = 0.5f * (F32_log((dot_0(mean_c_7, mean_c_7) + 9.99999997475242708e-07f)));
    *conic_7 = make_float3 (_S351.rows[int(0)].x, _S351.rows[int(0)].y, _S351.rows[int(1)].y);
    float3  _S352 = mean_7 - - mul_0(_S348, t_10);
    float3  _S353 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)];
    *rgb_7 = _S353;
    float _S354 = _S352.x;
    float _S355 = _S352.y;
    float _S356 = _S352.z;
    float norm_7 = (F32_sqrt((_S354 * _S354 + _S355 * _S355 + _S356 * _S356)));
    float x_29 = _S354 / norm_7;
    float y_10 = _S355 / norm_7;
    float z_7 = _S356 / norm_7;
    float3  _S357 = _S353 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * (*sh_coeffs_7)[int(1)] + make_float3 (z_7) * (*sh_coeffs_7)[int(2)] - make_float3 (x_29) * (*sh_coeffs_7)[int(3)]);
    *rgb_7 = _S357;
    float z2_18 = z_7 * z_7;
    float fTmp0B_7 = -1.09254848957061768f * z_7;
    float fC1_7 = x_29 * x_29 - y_10 * y_10;
    float fS1_7 = 2.0f * x_29 * y_10;
    float3  _S358 = _S357 + (make_float3 (0.54627424478530884f * fS1_7) * (*sh_coeffs_7)[int(4)] + make_float3 (fTmp0B_7 * y_10) * (*sh_coeffs_7)[int(5)] + make_float3 (0.94617468118667603f * z2_18 - 0.31539157032966614f) * (*sh_coeffs_7)[int(6)] + make_float3 (fTmp0B_7 * x_29) * (*sh_coeffs_7)[int(7)] + make_float3 (0.54627424478530884f * fC1_7) * (*sh_coeffs_7)[int(8)]);
    *rgb_7 = _S358;
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_7;
    *rgb_7 = max_0(_S358 + (make_float3 (-0.59004360437393188f * (x_29 * fS1_7 + y_10 * fC1_7)) * (*sh_coeffs_7)[int(9)] + make_float3 (fTmp1B_7 * fS1_7) * (*sh_coeffs_7)[int(10)] + make_float3 (fTmp0C_7 * y_10) * (*sh_coeffs_7)[int(11)] + make_float3 (z_7 * (1.86588168144226074f * z2_18 - 1.11952900886535645f)) * (*sh_coeffs_7)[int(12)] + make_float3 (fTmp0C_7 * x_29) * (*sh_coeffs_7)[int(13)] + make_float3 (fTmp1B_7 * fC1_7) * (*sh_coeffs_7)[int(14)] + make_float3 (-0.59004360437393188f * (x_29 * fC1_7 - y_10 * fS1_7)) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_persp_differentiable(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_12, float fy_12, float cx_12, float cy_12, float4  radial_coeffs_17, float2  tangential_coeffs_17, float2  thin_prism_coeffs_17, uint image_width_8, uint image_height_8, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * mean2d_12, float * depth_8, float3  * conic_8, float * opacity_8, float3  * rgb_8)
{
    float3  mean_c_8 = mul_0(R_12, mean_8) + t_11;
    float3  _S359 = exp_0(scale_10);
    float x_30 = quat_11.y;
    float x2_11 = x_30 * x_30;
    float y2_11 = quat_11.z * quat_11.z;
    float z2_19 = quat_11.w * quat_11.w;
    float xy_11 = quat_11.y * quat_11.z;
    float xz_11 = quat_11.y * quat_11.w;
    float yz_11 = quat_11.z * quat_11.w;
    float wx_11 = quat_11.x * quat_11.y;
    float wy_11 = quat_11.x * quat_11.z;
    float wz_11 = quat_11.x * quat_11.w;
    Matrix<float, 3, 3>  M_10 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))), makeMatrix<float, 3, 3> (_S359.x, 0.0f, 0.0f, 0.0f, _S359.y, 0.0f, 0.0f, 0.0f, _S359.z));
    Matrix<float, 3, 3>  _S360 = transpose_0(R_12);
    float _S361 = float(image_width_8);
    float _S362 = float(image_height_8);
    float _S363 = 0.30000001192092896f * (0.5f * _S361 / fx_12);
    float _S364 = 0.30000001192092896f * (0.5f * _S362 / fy_12);
    float rz_5 = 1.0f / mean_c_8.z;
    float rz2_5 = rz_5 * rz_5;
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (fx_12 * rz_5, 0.0f, - fx_12 * (mean_c_8.z * (F32_min(((_S361 - cx_12) / fx_12 + _S363), ((F32_max((- (cx_12 / fx_12 + _S363)), (mean_c_8.x * rz_5))))))) * rz2_5, 0.0f, fy_12 * rz_5, - fy_12 * (mean_c_8.z * (F32_min(((_S362 - cy_12) / fy_12 + _S364), ((F32_max((- (cy_12 / fy_12 + _S364)), (mean_c_8.y * rz_5))))))) * rz2_5);
    Matrix<float, 2, 2>  covar2d_8 = mul_6(mul_5(J_10, mul_4(mul_4(R_12, mul_4(M_10, transpose_0(M_10))), _S360)), transpose_1(J_10));
    *mean2d_12 = make_float2 (fx_12 * mean_c_8.x * rz_5 + cx_12, fy_12 * mean_c_8.y * rz_5 + cy_12);
    float det_orig_9 = *&(((&covar2d_8)->rows + (int(0)))->x) * *&(((&covar2d_8)->rows + (int(1)))->y) - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x);
    float _S365 = *&(((&covar2d_8)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(0)))->x) = _S365;
    float _S366 = *&(((&covar2d_8)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(1)))->y) = _S366;
    float compensation_9 = (F32_sqrt(((F32_max((0.0f), (det_orig_9 / (_S365 * _S366 - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x))))))));
    *opacity_8 = 1.0f / (1.0f + (F32_exp((- in_opacity_8))));
    if(antialiased_8)
    {
        *opacity_8 = *opacity_8 * compensation_9;
    }
    float extend_8 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_8 / 0.00392156885936856f)))))))));
    float radius_x_8 = extend_8 * (F32_sqrt((covar2d_8[int(0)].x)));
    float radius_y_8 = extend_8 * (F32_sqrt((covar2d_8[int(1)].y)));
    *aabb_xyxy_8 = make_int4 (int((F32_floor(((*mean2d_12).x - radius_x_8)))), int((F32_floor(((*mean2d_12).y - radius_y_8)))), int((F32_ceil(((*mean2d_12).x + radius_x_8)))), int((F32_ceil(((*mean2d_12).y + radius_y_8)))));
    *depth_8 = 0.5f * (F32_log((dot_0(mean_c_8, mean_c_8) + 9.99999997475242708e-07f)));
    *conic_8 = exp_0(- scale_10);
    float3  _S367 = mean_8 - - mul_0(_S360, t_11);
    float3  _S368 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)];
    *rgb_8 = _S368;
    float _S369 = _S367.x;
    float _S370 = _S367.y;
    float _S371 = _S367.z;
    float norm_8 = (F32_sqrt((_S369 * _S369 + _S370 * _S370 + _S371 * _S371)));
    float x_31 = _S369 / norm_8;
    float y_11 = _S370 / norm_8;
    float z_8 = _S371 / norm_8;
    float3  _S372 = _S368 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_8)[int(1)] + make_float3 (z_8) * (*sh_coeffs_8)[int(2)] - make_float3 (x_31) * (*sh_coeffs_8)[int(3)]);
    *rgb_8 = _S372;
    float z2_20 = z_8 * z_8;
    float fTmp0B_8 = -1.09254848957061768f * z_8;
    float fC1_8 = x_31 * x_31 - y_11 * y_11;
    float fS1_8 = 2.0f * x_31 * y_11;
    float3  _S373 = _S372 + (make_float3 (0.54627424478530884f * fS1_8) * (*sh_coeffs_8)[int(4)] + make_float3 (fTmp0B_8 * y_11) * (*sh_coeffs_8)[int(5)] + make_float3 (0.94617468118667603f * z2_20 - 0.31539157032966614f) * (*sh_coeffs_8)[int(6)] + make_float3 (fTmp0B_8 * x_31) * (*sh_coeffs_8)[int(7)] + make_float3 (0.54627424478530884f * fC1_8) * (*sh_coeffs_8)[int(8)]);
    *rgb_8 = _S373;
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_8;
    *rgb_8 = max_0(_S373 + (make_float3 (-0.59004360437393188f * (x_31 * fS1_8 + y_11 * fC1_8)) * (*sh_coeffs_8)[int(9)] + make_float3 (fTmp1B_8 * fS1_8) * (*sh_coeffs_8)[int(10)] + make_float3 (fTmp0C_8 * y_11) * (*sh_coeffs_8)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_20 - 1.11952900886535645f)) * (*sh_coeffs_8)[int(12)] + make_float3 (fTmp0C_8 * x_31) * (*sh_coeffs_8)[int(13)] + make_float3 (fTmp1B_8 * fC1_8) * (*sh_coeffs_8)[int(14)] + make_float3 (-0.59004360437393188f * (x_31 * fC1_8 - y_11 * fS1_8)) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_fisheye_differentiable(bool antialiased_9, float3  mean_9, float4  quat_12, float3  scale_11, float in_opacity_9, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_13, float fy_13, float cx_13, float cy_13, float4  radial_coeffs_18, float2  tangential_coeffs_18, float2  thin_prism_coeffs_18, uint image_width_9, uint image_height_9, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * mean2d_13, float * depth_9, float3  * conic_9, float * opacity_9, float3  * rgb_9)
{
    float3  mean_c_9 = mul_0(R_13, mean_9) + t_12;
    float3  _S374 = exp_0(scale_11);
    float x_32 = quat_12.y;
    float x2_12 = x_32 * x_32;
    float y2_12 = quat_12.z * quat_12.z;
    float z2_21 = quat_12.w * quat_12.w;
    float xy_12 = quat_12.y * quat_12.z;
    float xz_12 = quat_12.y * quat_12.w;
    float yz_12 = quat_12.z * quat_12.w;
    float wx_12 = quat_12.x * quat_12.y;
    float wy_12 = quat_12.x * quat_12.z;
    float wz_12 = quat_12.x * quat_12.w;
    Matrix<float, 3, 3>  M_11 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))), makeMatrix<float, 3, 3> (_S374.x, 0.0f, 0.0f, 0.0f, _S374.y, 0.0f, 0.0f, 0.0f, _S374.z));
    Matrix<float, 3, 3>  _S375 = transpose_0(R_13);
    Matrix<float, 2, 2>  covar2d_9;
    fisheye_proj_3dgs(mean_c_9, mul_4(mul_4(R_13, mul_4(M_11, transpose_0(M_11))), _S375), fx_13, fy_13, cx_13, cy_13, radial_coeffs_18, tangential_coeffs_18, thin_prism_coeffs_18, &covar2d_9, mean2d_13);
    float det_orig_10 = *&(((&covar2d_9)->rows + (int(0)))->x) * *&(((&covar2d_9)->rows + (int(1)))->y) - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x);
    float _S376 = *&(((&covar2d_9)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(0)))->x) = _S376;
    float _S377 = *&(((&covar2d_9)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(1)))->y) = _S377;
    float compensation_10 = (F32_sqrt(((F32_max((0.0f), (det_orig_10 / (_S376 * _S377 - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x))))))));
    *opacity_9 = 1.0f / (1.0f + (F32_exp((- in_opacity_9))));
    if(antialiased_9)
    {
        *opacity_9 = *opacity_9 * compensation_10;
    }
    float extend_9 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_9 / 0.00392156885936856f)))))))));
    float radius_x_9 = extend_9 * (F32_sqrt((covar2d_9[int(0)].x)));
    float radius_y_9 = extend_9 * (F32_sqrt((covar2d_9[int(1)].y)));
    *aabb_xyxy_9 = make_int4 (int((F32_floor(((*mean2d_13).x - radius_x_9)))), int((F32_floor(((*mean2d_13).y - radius_y_9)))), int((F32_ceil(((*mean2d_13).x + radius_x_9)))), int((F32_ceil(((*mean2d_13).y + radius_y_9)))));
    *depth_9 = 0.5f * (F32_log((dot_0(mean_c_9, mean_c_9) + 9.99999997475242708e-07f)));
    *conic_9 = exp_0(- scale_11);
    float3  _S378 = mean_9 - - mul_0(_S375, t_12);
    float3  _S379 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
    *rgb_9 = _S379;
    float _S380 = _S378.x;
    float _S381 = _S378.y;
    float _S382 = _S378.z;
    float norm_9 = (F32_sqrt((_S380 * _S380 + _S381 * _S381 + _S382 * _S382)));
    float x_33 = _S380 / norm_9;
    float y_12 = _S381 / norm_9;
    float z_9 = _S382 / norm_9;
    float3  _S383 = _S379 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * (*sh_coeffs_9)[int(1)] + make_float3 (z_9) * (*sh_coeffs_9)[int(2)] - make_float3 (x_33) * (*sh_coeffs_9)[int(3)]);
    *rgb_9 = _S383;
    float z2_22 = z_9 * z_9;
    float fTmp0B_9 = -1.09254848957061768f * z_9;
    float fC1_9 = x_33 * x_33 - y_12 * y_12;
    float fS1_9 = 2.0f * x_33 * y_12;
    float3  _S384 = _S383 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_12) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_22 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_33) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
    *rgb_9 = _S384;
    float fTmp0C_9 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_9;
    *rgb_9 = max_0(_S384 + (make_float3 (-0.59004360437393188f * (x_33 * fS1_9 + y_12 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_12) * (*sh_coeffs_9)[int(11)] + make_float3 (z_9 * (1.86588168144226074f * z2_22 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_33) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_33 * fC1_9 - y_12 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S385, float3  _S386)
{
    return mul_0(_S385, _S386);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S387)
{
    return exp_0(_S387);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S388, Matrix<float, 3, 3>  _S389)
{
    return mul_4(_S388, _S389);
}

inline __device__ float s_primal_ctx_max_0(float _S390, float _S391)
{
    return (F32_max((_S390), (_S391)));
}

inline __device__ float s_primal_ctx_min_0(float _S392, float _S393)
{
    return (F32_min((_S392), (_S393)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S394, Matrix<float, 3, 3>  _S395)
{
    return mul_5(_S394, _S395);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_4(Matrix<float, 2, 3>  _S396, Matrix<float, 3, 2>  _S397)
{
    return mul_6(_S396, _S397);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S398)
{
    return (F32_sqrt((_S398)));
}

inline __device__ float s_primal_ctx_exp_1(float _S399)
{
    return (F32_exp((_S399)));
}

inline __device__ float s_primal_ctx_log_0(float _S400)
{
    return (F32_log((_S400)));
}

inline __device__ float s_primal_ctx_dot_0(float3  _S401, float3  _S402)
{
    return dot_0(_S401, _S402);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S403, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S404, float3  _S405)
{
    _d_max_vector_0(_S403, _S404, _S405);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S406, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S407, float3  _S408)
{
    _d_mul_0(_S406, _S407, _S408);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S409, float _S410)
{
    _d_log_0(_S409, _S410);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S411, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S412, float _S413)
{
    _d_dot_0(_S411, _S412, _S413);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S414, DiffPair_float_0 * _S415, float _S416)
{
    _d_min_0(_S414, _S415, _S416);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S417, float _S418)
{
    _d_exp_0(_S417, _S418);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S419, DiffPair_float_0 * _S420, float _S421)
{
    _d_max_0(_S419, _S420, _S421);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S422, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S423, Matrix<float, 2, 2>  _S424)
{
    mul_3(_S422, _S423, _S424);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S425, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S426, Matrix<float, 2, 3>  _S427)
{
    mul_2(_S425, _S426, _S427);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S428, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S429, Matrix<float, 3, 3>  _S430)
{
    mul_1(_S428, _S429, _S430);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S431, float3  _S432)
{
    _d_exp_vector_0(_S431, _S432);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_14, float fy_14, float cx_14, float cy_14, float4  radial_coeffs_19, float2  tangential_coeffs_19, float2  thin_prism_coeffs_19, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_10 = s_primal_ctx_mul_1(R_14, mean_10) + t_13;
    float3  _S433 = s_primal_ctx_exp_0(scale_12);
    float _S434 = quat_13.y;
    float x2_13 = _S434 * _S434;
    float y2_13 = quat_13.z * quat_13.z;
    float z2_23 = quat_13.w * quat_13.w;
    float xy_13 = quat_13.y * quat_13.z;
    float xz_13 = quat_13.y * quat_13.w;
    float yz_13 = quat_13.z * quat_13.w;
    float wx_13 = quat_13.x * quat_13.y;
    float wy_13 = quat_13.x * quat_13.z;
    float wz_13 = quat_13.x * quat_13.w;
    Matrix<float, 3, 3>  _S435 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S433.x, 0.0f, 0.0f, 0.0f, _S433.y, 0.0f, 0.0f, 0.0f, _S433.z);
    Matrix<float, 3, 3>  _S436 = s_primal_ctx_mul_2(_S435, S_0);
    Matrix<float, 3, 3>  _S437 = transpose_0(_S436);
    Matrix<float, 3, 3>  _S438 = s_primal_ctx_mul_2(_S436, _S437);
    Matrix<float, 3, 3>  _S439 = s_primal_ctx_mul_2(R_14, _S438);
    Matrix<float, 3, 3>  _S440 = transpose_0(R_14);
    Matrix<float, 3, 3>  _S441 = s_primal_ctx_mul_2(_S439, _S440);
    float _S442 = float(image_width_10);
    float _S443 = float(image_height_10);
    float _S444 = 0.30000001192092896f * (0.5f * _S442 / fx_14);
    float lim_x_pos_0 = (_S442 - cx_14) / fx_14 + _S444;
    float _S445 = 0.30000001192092896f * (0.5f * _S443 / fy_14);
    float lim_y_pos_0 = (_S443 - cy_14) / fy_14 + _S445;
    float rz_6 = 1.0f / mean_c_10.z;
    float _S446 = mean_c_10.z * mean_c_10.z;
    float rz2_6 = rz_6 * rz_6;
    float _S447 = - (cx_14 / fx_14 + _S444);
    float _S448 = mean_c_10.x * rz_6;
    float _S449 = s_primal_ctx_max_0(_S447, _S448);
    float _S450 = s_primal_ctx_min_0(lim_x_pos_0, _S449);
    float _S451 = - (cy_14 / fy_14 + _S445);
    float _S452 = mean_c_10.y * rz_6;
    float _S453 = s_primal_ctx_max_0(_S451, _S452);
    float _S454 = s_primal_ctx_min_0(lim_y_pos_0, _S453);
    float _S455 = - fx_14;
    float _S456 = _S455 * (mean_c_10.z * _S450);
    float _S457 = - fy_14;
    float _S458 = _S457 * (mean_c_10.z * _S454);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_14 * rz_6, 0.0f, _S456 * rz2_6, 0.0f, fy_14 * rz_6, _S458 * rz2_6);
    Matrix<float, 2, 3>  _S459 = s_primal_ctx_mul_3(J_11, _S441);
    Matrix<float, 3, 2>  _S460 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S461 = s_primal_ctx_mul_4(_S459, _S460);
    float _S462 = fx_14 * mean_c_10.x;
    float _S463 = fy_14 * mean_c_10.y;
    float _S464 = _S461.rows[int(0)].y * _S461.rows[int(1)].x;
    float det_orig_11 = _S461.rows[int(0)].x * _S461.rows[int(1)].y - _S464;
    float _S465 = _S461.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S466 = _S461;
    *&(((&_S466)->rows + (int(0)))->x) = _S465;
    float _S467 = _S461.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S466)->rows + (int(1)))->y) = _S467;
    Matrix<float, 2, 2>  _S468 = _S466;
    Matrix<float, 2, 2>  _S469 = _S466;
    float det_blur_6 = _S465 * _S467 - _S464;
    float _S470 = det_orig_11 / det_blur_6;
    float _S471 = det_blur_6 * det_blur_6;
    float _S472 = s_primal_ctx_max_0(0.0f, _S470);
    float _S473 = s_primal_ctx_sqrt_0(_S472);
    float invdet_7 = 1.0f / det_blur_6;
    float _S474 = - _S461.rows[int(0)].y;
    float _S475 = - _S461.rows[int(1)].x;
    float _S476 = - in_opacity_10;
    float _S477 = 1.0f + s_primal_ctx_exp_1(_S476);
    float _S478 = 1.0f / _S477;
    float _S479 = _S477 * _S477;
    float _S480;
    if(antialiased_10)
    {
        _S480 = _S478 * _S473;
    }
    else
    {
        _S480 = _S478;
    }
    float _S481 = _S480 / 0.00392156885936856f;
    float _S482 = 2.0f * s_primal_ctx_log_0(_S481);
    float _S483 = s_primal_ctx_sqrt_0(_S482);
    float _S484 = _S468.rows[int(0)].x;
    float _S485 = _S469.rows[int(1)].y;
    float _S486 = s_primal_ctx_dot_0(mean_c_10, mean_c_10) + 9.99999997475242708e-07f;
    float3  _S487 = mean_10 - - s_primal_ctx_mul_1(_S440, t_13);
    float _S488 = _S487.x;
    float _S489 = _S487.y;
    float _S490 = _S487.z;
    float _S491 = _S488 * _S488 + _S489 * _S489 + _S490 * _S490;
    float _S492 = s_primal_ctx_sqrt_0(_S491);
    float x_34 = _S488 / _S492;
    float3  _S493 = make_float3 (x_34);
    float _S494 = _S492 * _S492;
    float y_13 = _S489 / _S492;
    float z_10 = _S490 / _S492;
    float3  _S495 = make_float3 (z_10);
    float _S496 = - y_13;
    float3  _S497 = make_float3 (_S496);
    float z2_24 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_34 * x_34 - y_13 * y_13;
    float _S498 = 2.0f * x_34;
    float fS1_10 = _S498 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_24 - 0.31539157032966614f;
    float3  _S499 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_34;
    float3  _S500 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S501 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S502 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S503 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S504 = 1.86588168144226074f * z2_24 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S504;
    float3  _S505 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_34;
    float3  _S506 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S507 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S508 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S509 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_34 * fC1_10 - y_13 * fS1_10);
    float3  _S510 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_34 * fS1_10 + y_13 * fC1_10);
    float3  _S511 = make_float3 (pSH9_0);
    float3  _S512 = make_float3 (0.0f);
    float3  _S513 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S514;
    (&_S514)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S496) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_34) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S514)->differential_0 = _S513;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S515;
    (&_S515)->primal_0 = _S512;
    (&_S515)->differential_0 = _S513;
    s_bwd_prop_max_0(&_S514, &_S515, v_rgb_0);
    float3  _S516 = _S510 * _S514.differential_0;
    float3  _S517 = (*sh_coeffs_10)[int(15)] * _S514.differential_0;
    float3  _S518 = _S508 * _S514.differential_0;
    float3  _S519 = (*sh_coeffs_10)[int(14)] * _S514.differential_0;
    float3  _S520 = _S506 * _S514.differential_0;
    float3  _S521 = (*sh_coeffs_10)[int(13)] * _S514.differential_0;
    float3  _S522 = _S505 * _S514.differential_0;
    float3  _S523 = (*sh_coeffs_10)[int(12)] * _S514.differential_0;
    float3  _S524 = _S507 * _S514.differential_0;
    float3  _S525 = (*sh_coeffs_10)[int(11)] * _S514.differential_0;
    float3  _S526 = _S509 * _S514.differential_0;
    float3  _S527 = (*sh_coeffs_10)[int(10)] * _S514.differential_0;
    float3  _S528 = _S511 * _S514.differential_0;
    float3  _S529 = (*sh_coeffs_10)[int(9)] * _S514.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S529.x + _S529.y + _S529.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S517.x + _S517.y + _S517.z);
    float _S530 = _S527.x + _S527.y + _S527.z;
    float _S531 = _S519.x + _S519.y + _S519.z;
    float _S532 = _S525.x + _S525.y + _S525.z;
    float _S533 = _S521.x + _S521.y + _S521.z;
    float _S534 = _S523.x + _S523.y + _S523.z;
    float _S535 = - s_diff_fC2_T_0;
    float3  _S536 = _S502 * _S514.differential_0;
    float3  _S537 = (*sh_coeffs_10)[int(8)] * _S514.differential_0;
    float3  _S538 = _S500 * _S514.differential_0;
    float3  _S539 = (*sh_coeffs_10)[int(7)] * _S514.differential_0;
    float3  _S540 = _S499 * _S514.differential_0;
    float3  _S541 = (*sh_coeffs_10)[int(6)] * _S514.differential_0;
    float3  _S542 = _S501 * _S514.differential_0;
    float3  _S543 = (*sh_coeffs_10)[int(5)] * _S514.differential_0;
    float3  _S544 = _S503 * _S514.differential_0;
    float3  _S545 = (*sh_coeffs_10)[int(4)] * _S514.differential_0;
    float _S546 = _S543.x + _S543.y + _S543.z;
    float _S547 = _S539.x + _S539.y + _S539.z;
    float _S548 = fTmp1B_10 * _S530 + x_34 * s_diff_fS2_T_0 + y_13 * _S535 + 0.54627424478530884f * (_S545.x + _S545.y + _S545.z);
    float _S549 = fTmp1B_10 * _S531 + y_13 * s_diff_fS2_T_0 + x_34 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S537.x + _S537.y + _S537.z);
    float _S550 = y_13 * - _S549;
    float _S551 = x_34 * _S549;
    float _S552 = z_10 * (1.86588168144226074f * (z_10 * _S534) + -2.28522896766662598f * (y_13 * _S532 + x_34 * _S533) + 0.94617468118667603f * (_S541.x + _S541.y + _S541.z));
    float3  _S553 = make_float3 (0.48860251903533936f) * _S514.differential_0;
    float3  _S554 = - _S553;
    float3  _S555 = _S493 * _S554;
    float3  _S556 = (*sh_coeffs_10)[int(3)] * _S554;
    float3  _S557 = _S495 * _S553;
    float3  _S558 = (*sh_coeffs_10)[int(2)] * _S553;
    float3  _S559 = _S497 * _S553;
    float3  _S560 = (*sh_coeffs_10)[int(1)] * _S553;
    float _S561 = (_S504 * _S534 + 1.44530570507049561f * (fS1_10 * _S530 + fC1_10 * _S531) + -1.09254848957061768f * (y_13 * _S546 + x_34 * _S547) + _S552 + _S552 + _S558.x + _S558.y + _S558.z) / _S494;
    float _S562 = _S492 * _S561;
    float _S563 = (fTmp0C_10 * _S532 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S535 + fTmp0B_10 * _S546 + _S498 * _S548 + _S550 + _S550 + - (_S560.x + _S560.y + _S560.z)) / _S494;
    float _S564 = _S492 * _S563;
    float _S565 = (fTmp0C_10 * _S533 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S547 + 2.0f * (y_13 * _S548) + _S551 + _S551 + _S556.x + _S556.y + _S556.z) / _S494;
    float _S566 = _S492 * _S565;
    float _S567 = _S490 * - _S561 + _S489 * - _S563 + _S488 * - _S565;
    DiffPair_float_0 _S568;
    (&_S568)->primal_0 = _S491;
    (&_S568)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S568, _S567);
    float _S569 = _S490 * _S568.differential_0;
    float _S570 = _S489 * _S568.differential_0;
    float _S571 = _S488 * _S568.differential_0;
    float3  _S572 = make_float3 (0.282094806432724f) * _S514.differential_0;
    float3  _S573 = make_float3 (_S566 + _S571 + _S571, _S564 + _S570 + _S570, _S562 + _S569 + _S569);
    float3  _S574 = - - _S573;
    Matrix<float, 3, 3>  _S575 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S576;
    (&_S576)->primal_0 = _S440;
    (&_S576)->differential_0 = _S575;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S577;
    (&_S577)->primal_0 = t_13;
    (&_S577)->differential_0 = _S513;
    s_bwd_prop_mul_1(&_S576, &_S577, _S574);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S578 = _S576;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S579 = _S577;
    float2  _S580 = make_float2 (0.0f);
    float2  _S581 = _S580;
    *&((&_S581)->y) = v_conic_0.z;
    float2  _S582 = _S580;
    *&((&_S582)->y) = v_conic_0.y;
    *&((&_S582)->x) = v_conic_0.x;
    float _S583 = 0.5f * v_depth_0;
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S486;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S584, _S583);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S585;
    (&_S585)->primal_0 = mean_c_10;
    (&_S585)->differential_0 = _S513;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S586;
    (&_S586)->primal_0 = mean_c_10;
    (&_S586)->differential_0 = _S513;
    s_bwd_prop_dot_0(&_S585, &_S586, _S584.differential_0);
    DiffPair_float_0 _S587;
    (&_S587)->primal_0 = _S485;
    (&_S587)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S587, 0.0f);
    DiffPair_float_0 _S588;
    (&_S588)->primal_0 = _S484;
    (&_S588)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S588, 0.0f);
    DiffPair_float_0 _S589;
    (&_S589)->primal_0 = 3.32999992370605469f;
    (&_S589)->differential_0 = 0.0f;
    DiffPair_float_0 _S590;
    (&_S590)->primal_0 = _S483;
    (&_S590)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S589, &_S590, 0.0f);
    DiffPair_float_0 _S591;
    (&_S591)->primal_0 = _S482;
    (&_S591)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S591, _S590.differential_0);
    float _S592 = 2.0f * _S591.differential_0;
    DiffPair_float_0 _S593;
    (&_S593)->primal_0 = _S481;
    (&_S593)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S593, _S592);
    float _S594 = v_opacity_0 + 254.9999847412109375f * _S593.differential_0;
    Matrix<float, 2, 2>  _S595 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S596 = _S595;
    _S596[int(1)] = _S581;
    _S596[int(0)] = _S582;
    Matrix<float, 2, 2>  _S597 = _S596;
    FixedArray<float3 , 16>  _S598;
    _S598[int(0)] = _S513;
    _S598[int(1)] = _S513;
    _S598[int(2)] = _S513;
    _S598[int(3)] = _S513;
    _S598[int(4)] = _S513;
    _S598[int(5)] = _S513;
    _S598[int(6)] = _S513;
    _S598[int(7)] = _S513;
    _S598[int(8)] = _S513;
    _S598[int(9)] = _S513;
    _S598[int(10)] = _S513;
    _S598[int(11)] = _S513;
    _S598[int(12)] = _S513;
    _S598[int(13)] = _S513;
    _S598[int(14)] = _S513;
    _S598[int(15)] = _S513;
    _S598[int(7)] = _S538;
    _S598[int(0)] = _S572;
    _S598[int(1)] = _S559;
    _S598[int(2)] = _S557;
    _S598[int(3)] = _S555;
    _S598[int(4)] = _S544;
    _S598[int(5)] = _S542;
    _S598[int(6)] = _S540;
    _S598[int(15)] = _S516;
    _S598[int(8)] = _S536;
    _S598[int(9)] = _S528;
    _S598[int(10)] = _S526;
    _S598[int(11)] = _S524;
    _S598[int(12)] = _S522;
    _S598[int(13)] = _S520;
    _S598[int(14)] = _S518;
    float3  _S599 = _S598[int(0)];
    float3  _S600 = _S598[int(1)];
    float3  _S601 = _S598[int(2)];
    float3  _S602 = _S598[int(3)];
    float3  _S603 = _S598[int(4)];
    float3  _S604 = _S598[int(5)];
    float3  _S605 = _S598[int(6)];
    float3  _S606 = _S598[int(7)];
    float3  _S607 = _S598[int(8)];
    float3  _S608 = _S598[int(9)];
    float3  _S609 = _S598[int(10)];
    float3  _S610 = _S598[int(11)];
    float3  _S611 = _S598[int(12)];
    float3  _S612 = _S598[int(13)];
    float3  _S613 = _S598[int(14)];
    float3  _S614 = _S598[int(15)];
    float3  _S615 = _S586.differential_0 + _S585.differential_0;
    float2  _S616 = make_float2 (0.0f, _S587.differential_0);
    float2  _S617 = make_float2 (_S588.differential_0, 0.0f);
    float _S618;
    if(antialiased_10)
    {
        float _S619 = _S478 * _S594;
        _S480 = _S473 * _S594;
        _S618 = _S619;
    }
    else
    {
        _S480 = _S594;
        _S618 = 0.0f;
    }
    float _S620 = - (_S480 / _S479);
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = _S476;
    (&_S621)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S621, _S620);
    float _S622 = - _S621.differential_0;
    float _S623 = invdet_7 * _S597.rows[int(1)].y;
    float _S624 = - (invdet_7 * _S597.rows[int(1)].x);
    float _S625 = - (invdet_7 * _S597.rows[int(0)].y);
    float _S626 = invdet_7 * _S597.rows[int(0)].x;
    float _S627 = - ((_S465 * _S597.rows[int(1)].y + _S475 * _S597.rows[int(1)].x + _S474 * _S597.rows[int(0)].y + _S467 * _S597.rows[int(0)].x) / _S471);
    DiffPair_float_0 _S628;
    (&_S628)->primal_0 = _S472;
    (&_S628)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S628, _S618);
    DiffPair_float_0 _S629;
    (&_S629)->primal_0 = 0.0f;
    (&_S629)->differential_0 = 0.0f;
    DiffPair_float_0 _S630;
    (&_S630)->primal_0 = _S470;
    (&_S630)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S629, &_S630, _S628.differential_0);
    float _S631 = _S630.differential_0 / _S471;
    float s_diff_det_orig_T_0 = det_blur_6 * _S631;
    float _S632 = _S627 + det_orig_11 * - _S631;
    float _S633 = - _S632;
    float _S634 = _S465 * _S632;
    float _S635 = _S467 * _S632;
    Matrix<float, 2, 2>  _S636 = _S595;
    _S636[int(1)] = _S616;
    _S636[int(0)] = _S617;
    _S466 = _S636;
    *&(((&_S466)->rows + (int(1)))->y) = 0.0f;
    float _S637 = _S626 + _S634 + _S636.rows[int(1)].y;
    *&(((&_S466)->rows + (int(0)))->x) = 0.0f;
    float _S638 = _S623 + _S635 + _S636.rows[int(0)].x;
    float _S639 = _S633 + - s_diff_det_orig_T_0;
    float _S640 = _S624 + _S461.rows[int(0)].y * _S639;
    float _S641 = _S625 + _S461.rows[int(1)].x * _S639;
    float _S642 = _S461.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S643 = _S637 + _S461.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S644 = _S580;
    *&((&_S644)->x) = _S640;
    *&((&_S644)->y) = _S643;
    float _S645 = _S638 + _S642;
    float2  _S646 = _S580;
    *&((&_S646)->y) = _S641;
    *&((&_S646)->x) = _S645;
    float _S647 = _S463 * v_mean2d_0.y;
    float _S648 = fy_14 * (rz_6 * v_mean2d_0.y);
    float _S649 = _S462 * v_mean2d_0.x;
    float _S650 = fx_14 * (rz_6 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S651 = _S595;
    _S651[int(1)] = _S644;
    _S651[int(0)] = _S646;
    Matrix<float, 2, 2>  _S652 = _S466 + _S651;
    Matrix<float, 2, 3>  _S653 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S654;
    (&_S654)->primal_0 = _S459;
    (&_S654)->differential_0 = _S653;
    Matrix<float, 3, 2>  _S655 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S656;
    (&_S656)->primal_0 = _S460;
    (&_S656)->differential_0 = _S655;
    s_bwd_prop_mul_2(&_S654, &_S656, _S652);
    Matrix<float, 2, 3>  _S657 = transpose_2(_S656.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S658;
    (&_S658)->primal_0 = J_11;
    (&_S658)->differential_0 = _S653;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S659;
    (&_S659)->primal_0 = _S441;
    (&_S659)->differential_0 = _S575;
    s_bwd_prop_mul_3(&_S658, &_S659, _S654.differential_0);
    Matrix<float, 2, 3>  _S660 = _S657 + _S658.differential_0;
    float _S661 = _S458 * _S660.rows[int(1)].z;
    float s_diff_ty_T_0 = _S457 * (rz2_6 * _S660.rows[int(1)].z);
    float _S662 = fy_14 * _S660.rows[int(1)].y;
    float _S663 = _S456 * _S660.rows[int(0)].z;
    float s_diff_tx_T_0 = _S455 * (rz2_6 * _S660.rows[int(0)].z);
    float _S664 = fx_14 * _S660.rows[int(0)].x;
    float _S665 = mean_c_10.z * s_diff_ty_T_0;
    float _S666 = _S454 * s_diff_ty_T_0;
    DiffPair_float_0 _S667;
    (&_S667)->primal_0 = lim_y_pos_0;
    (&_S667)->differential_0 = 0.0f;
    DiffPair_float_0 _S668;
    (&_S668)->primal_0 = _S453;
    (&_S668)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S667, &_S668, _S665);
    DiffPair_float_0 _S669;
    (&_S669)->primal_0 = _S451;
    (&_S669)->differential_0 = 0.0f;
    DiffPair_float_0 _S670;
    (&_S670)->primal_0 = _S452;
    (&_S670)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S669, &_S670, _S668.differential_0);
    float _S671 = mean_c_10.y * _S670.differential_0;
    float _S672 = rz_6 * _S670.differential_0;
    float _S673 = mean_c_10.z * s_diff_tx_T_0;
    float _S674 = _S450 * s_diff_tx_T_0;
    DiffPair_float_0 _S675;
    (&_S675)->primal_0 = lim_x_pos_0;
    (&_S675)->differential_0 = 0.0f;
    DiffPair_float_0 _S676;
    (&_S676)->primal_0 = _S449;
    (&_S676)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S675, &_S676, _S673);
    DiffPair_float_0 _S677;
    (&_S677)->primal_0 = _S447;
    (&_S677)->differential_0 = 0.0f;
    DiffPair_float_0 _S678;
    (&_S678)->primal_0 = _S448;
    (&_S678)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S677, &_S678, _S676.differential_0);
    float _S679 = rz_6 * (_S661 + _S663);
    float _S680 = _S666 + _S674 + - ((_S647 + _S649 + _S662 + _S664 + _S671 + mean_c_10.x * _S678.differential_0 + _S679 + _S679) / _S446);
    float _S681 = _S648 + _S672;
    float _S682 = _S650 + rz_6 * _S678.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S683;
    (&_S683)->primal_0 = _S439;
    (&_S683)->differential_0 = _S575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S684;
    (&_S684)->primal_0 = _S440;
    (&_S684)->differential_0 = _S575;
    s_bwd_prop_mul_4(&_S683, &_S684, _S659.differential_0);
    Matrix<float, 3, 3>  _S685 = transpose_0(_S684.differential_0 + _S578.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S686;
    (&_S686)->primal_0 = R_14;
    (&_S686)->differential_0 = _S575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S687;
    (&_S687)->primal_0 = _S438;
    (&_S687)->differential_0 = _S575;
    s_bwd_prop_mul_4(&_S686, &_S687, _S683.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S688;
    (&_S688)->primal_0 = _S436;
    (&_S688)->differential_0 = _S575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S689;
    (&_S689)->primal_0 = _S437;
    (&_S689)->differential_0 = _S575;
    s_bwd_prop_mul_4(&_S688, &_S689, _S687.differential_0);
    Matrix<float, 3, 3>  _S690 = _S688.differential_0 + transpose_0(_S689.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S691;
    (&_S691)->primal_0 = _S435;
    (&_S691)->differential_0 = _S575;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S692;
    (&_S692)->primal_0 = S_0;
    (&_S692)->differential_0 = _S575;
    s_bwd_prop_mul_4(&_S691, &_S692, _S690);
    Matrix<float, 3, 3>  _S693 = transpose_0(_S691.differential_0);
    float _S694 = 2.0f * - _S693.rows[int(2)].z;
    float _S695 = 2.0f * _S693.rows[int(2)].y;
    float _S696 = 2.0f * _S693.rows[int(2)].x;
    float _S697 = 2.0f * _S693.rows[int(1)].z;
    float _S698 = 2.0f * - _S693.rows[int(1)].y;
    float _S699 = 2.0f * _S693.rows[int(1)].x;
    float _S700 = 2.0f * _S693.rows[int(0)].z;
    float _S701 = 2.0f * _S693.rows[int(0)].y;
    float _S702 = 2.0f * - _S693.rows[int(0)].x;
    float _S703 = - _S699 + _S701;
    float _S704 = _S696 + - _S700;
    float _S705 = - _S695 + _S697;
    float _S706 = _S695 + _S697;
    float _S707 = _S696 + _S700;
    float _S708 = _S699 + _S701;
    float _S709 = quat_13.w * (_S698 + _S702);
    float _S710 = quat_13.z * (_S694 + _S702);
    float _S711 = quat_13.y * (_S694 + _S698);
    float _S712 = quat_13.x * _S703 + quat_13.z * _S706 + quat_13.y * _S707 + _S709 + _S709;
    float _S713 = quat_13.x * _S704 + quat_13.w * _S706 + quat_13.y * _S708 + _S710 + _S710;
    float _S714 = quat_13.x * _S705 + quat_13.w * _S707 + quat_13.z * _S708 + _S711 + _S711;
    float _S715 = quat_13.w * _S703 + quat_13.z * _S704 + quat_13.y * _S705;
    float3  _S716 = _S513;
    *&((&_S716)->z) = _S692.differential_0.rows[int(2)].z;
    *&((&_S716)->y) = _S692.differential_0.rows[int(1)].y;
    *&((&_S716)->x) = _S692.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S717;
    (&_S717)->primal_0 = scale_12;
    (&_S717)->differential_0 = _S513;
    s_bwd_prop_exp_1(&_S717, _S716);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S718 = _S717;
    float3  _S719 = _S513;
    *&((&_S719)->z) = _S680;
    *&((&_S719)->y) = _S681;
    *&((&_S719)->x) = _S682;
    float3  _S720 = _S615 + _S719;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S721;
    (&_S721)->primal_0 = R_14;
    (&_S721)->differential_0 = _S575;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S722;
    (&_S722)->primal_0 = mean_10;
    (&_S722)->differential_0 = _S513;
    s_bwd_prop_mul_1(&_S721, &_S722, _S720);
    float3  _S723 = _S720 + _S579.differential_0;
    Matrix<float, 3, 3>  _S724 = _S685 + _S686.differential_0 + _S721.differential_0;
    float4  _S725 = make_float4 (0.0f);
    *&((&_S725)->w) = _S712;
    *&((&_S725)->z) = _S713;
    *&((&_S725)->y) = _S714;
    *&((&_S725)->x) = _S715;
    float4  _S726 = _S725;
    float3  _S727 = _S722.differential_0 + _S573;
    *v_mean_0 = _S727;
    *v_quat_0 = _S726;
    *v_scale_0 = _S718.differential_0;
    *v_in_opacity_0 = _S622;
    (*v_sh_coeffs_0)[int(0)] = _S599;
    (*v_sh_coeffs_0)[int(1)] = _S600;
    (*v_sh_coeffs_0)[int(2)] = _S601;
    (*v_sh_coeffs_0)[int(3)] = _S602;
    (*v_sh_coeffs_0)[int(4)] = _S603;
    (*v_sh_coeffs_0)[int(5)] = _S604;
    (*v_sh_coeffs_0)[int(6)] = _S605;
    (*v_sh_coeffs_0)[int(7)] = _S606;
    (*v_sh_coeffs_0)[int(8)] = _S607;
    (*v_sh_coeffs_0)[int(9)] = _S608;
    (*v_sh_coeffs_0)[int(10)] = _S609;
    (*v_sh_coeffs_0)[int(11)] = _S610;
    (*v_sh_coeffs_0)[int(12)] = _S611;
    (*v_sh_coeffs_0)[int(13)] = _S612;
    (*v_sh_coeffs_0)[int(14)] = _S613;
    (*v_sh_coeffs_0)[int(15)] = _S614;
    *v_R_1 = _S724;
    *v_t_1 = _S723;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S728;
    DiffPair_float_0 _S729;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S730;
    DiffPair_float_0 _S731;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S732;
    DiffPair_float_0 _S733;
    DiffPair_float_0 _S734;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S735;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S736;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S737;
};

inline __device__ CameraDistortion_0 s_primal_ctx_CameraDistortion_x24init_0(float4  dpradial_coeffs_0, float2  dptangential_coeffs_0, float2  dpthin_prism_coeffs_0)
{
    CameraDistortion_0 _S738 = { dpradial_coeffs_0, dptangential_coeffs_0, dpthin_prism_coeffs_0 };
    return _S738;
}

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S739, float _S740)
{
    return s_primal_ctx_atan2_0(_S739, _S740);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S741;
    DiffPair_float_0 _S742;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S743 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S741 = _S743;
    _s_diff_ctx_2->_S742 = _S743;
    (&_s_diff_ctx_2->_S741)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S741)->differential_0 = 0.0f;
    (&_s_diff_ctx_2->_S742)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S742)->differential_0 = 0.0f;
    DiffPair_float_0 _S744 = *dpdpy_0;
    _s_diff_ctx_2->_S741 = *dpdpy_0;
    DiffPair_float_0 _S745 = *dpdpx_0;
    _s_diff_ctx_2->_S742 = *dpdpx_0;
    float _S746 = _S745.primal_0 * _S745.primal_0 + _S744.primal_0 * _S744.primal_0;
    float _S747 = - _S744.primal_0 / _S746 * dpdOut_0;
    float _S748 = _S745.primal_0 / _S746 * dpdOut_0;
    dpdpy_0->primal_0 = _S744.primal_0;
    dpdpy_0->differential_0 = _S748;
    dpdpx_0->primal_0 = _S745.primal_0;
    dpdpx_0->differential_0 = _S747;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S749, DiffPair_float_0 * _S750, float _S751, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S752 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S728 = _S752;
    _s_diff_ctx_3->_S729 = _S752;
    (&_s_diff_ctx_3->_S728)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S728)->differential_0 = 0.0f;
    (&_s_diff_ctx_3->_S729)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S729)->differential_0 = 0.0f;
    DiffPair_float_0 _S753 = *_S749;
    _s_diff_ctx_3->_S728 = *_S749;
    DiffPair_float_0 _S754 = *_S750;
    _s_diff_ctx_3->_S729 = *_S750;
    DiffPair_float_0 _S755 = _S753;
    DiffPair_float_0 _S756 = _S754;
    s_bwd_prop_d_atan2_Intermediates_0 _S757;
    (&_S757)->_S741 = _S752;
    (&_S757)->_S742 = _S752;
    s_primal_ctx_d_atan2_0(&_S755, &_S756, _S751, &_S757);
    *_S749 = _S755;
    *_S750 = _S756;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S758;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S759;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S760;
    DiffPair_float_0 _S761;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S762;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S763;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S764 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S763 = _S764;
    (&_s_diff_ctx_4->_S763)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S763)->differential_0 = 0.0f;
    DiffPair_float_0 _S765 = *dpdpx_1;
    _s_diff_ctx_4->_S763 = *dpdpx_1;
    float _S766 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S765.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S765.primal_0;
    dpdpx_1->differential_0 = _S766;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S767, float _S768, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S769 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S759 = _S769;
    (&_s_diff_ctx_5->_S759)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S759)->differential_0 = 0.0f;
    DiffPair_float_0 _S770 = *_S767;
    _s_diff_ctx_5->_S759 = *_S767;
    DiffPair_float_0 _S771 = _S770;
    s_bwd_prop_d_sqrt_Intermediates_0 _S772;
    (&_S772)->_S763 = _S769;
    s_primal_ctx_d_sqrt_0(&_S771, _S768, &_S772);
    *_S767 = _S771;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_6)
{
    float2  _S773 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S774 = { _S773, _S773 };
    DiffPair_float_0 _S775 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S776 = { _S775 };
    _s_diff_ctx_6->_S760 = _S774;
    _s_diff_ctx_6->_S761 = _S775;
    _s_diff_ctx_6->_S762 = _S776;
    (&_s_diff_ctx_6->_S760)->primal_0 = _S773;
    (&_s_diff_ctx_6->_S760)->differential_0 = _S773;
    (&_s_diff_ctx_6->_S761)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S761)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S777 = *dpdpx_2;
    _s_diff_ctx_6->_S760 = *dpdpx_2;
    float _S778 = _S777.primal_0.x;
    float _S779 = _S777.primal_0.y;
    DiffPair_float_0 _S780;
    (&_S780)->primal_0 = _S778 * _S778 + _S779 * _S779;
    (&_S780)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S780, dp_s_dOut_0, &_s_diff_ctx_6->_S762);
    _s_diff_ctx_6->_S761 = _S780;
    float _S781 = _S777.primal_0.y * _S780.differential_0;
    float _S782 = _S781 + _S781;
    float _S783 = _S777.primal_0.x * _S780.differential_0;
    float _S784 = _S783 + _S783;
    float2  _S785 = _S773;
    *&((&_S785)->y) = _S782;
    *&((&_S785)->x) = _S784;
    dpdpx_2->primal_0 = _S777.primal_0;
    dpdpx_2->differential_0 = _S785;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S786, float _S787, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_7)
{
    float2  _S788 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S789 = { _S788, _S788 };
    _s_diff_ctx_7->_S758 = _S789;
    (&_s_diff_ctx_7->_S758)->primal_0 = _S788;
    (&_s_diff_ctx_7->_S758)->differential_0 = _S788;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S790 = *_S786;
    _s_diff_ctx_7->_S758 = *_S786;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S791 = _S790;
    DiffPair_float_0 _S792 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S793 = { _S792 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S794;
    (&_S794)->_S760 = _S789;
    (&_S794)->_S761 = _S792;
    (&_S794)->_S762 = _S793;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S791, _S787, &_S794);
    *_S786 = _S791;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, float4  dpradial_coeffs_1, float2  dptangential_coeffs_1, float2  dpthin_prism_coeffs_1, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_float_0 _S795 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S796 = { _S795, _S795 };
    _s_diff_ctx_8->_S730 = _S795;
    _s_diff_ctx_8->_S731 = _S795;
    _s_diff_ctx_8->_S732 = _S796;
    _s_diff_ctx_8->_S733 = _S795;
    _s_diff_ctx_8->_S734 = _S795;
    _s_diff_ctx_8->_S735 = _S796;
    (&_s_diff_ctx_8->_S730)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S730)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S731)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S731)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S733)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S733)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S734)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S734)->differential_0 = 0.0f;
    float2  _S797 = make_float2 (0.0f);
    CameraDistortion_0 _S798 = s_primal_ctx_CameraDistortion_x24init_0(dpradial_coeffs_1, dptangential_coeffs_1, dpthin_prism_coeffs_1);
    float2  _S799 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S800 = length_0(_S799);
    float _S801 = dpmean3d_0.z;
    float _S802 = s_primal_ctx_atan2_0(_S800, _S801);
    float k_1;
    if(_S802 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S802 * _S802 / 3.0f) / _S801;
    }
    else
    {
        k_1 = _S802 / _S800;
    }
    float2  _S803 = _S799 * make_float2 (k_1);
    float k1_1 = _S798.radial_coeffs_0.x;
    float k2_1 = _S798.radial_coeffs_0.y;
    float k3_1 = _S798.radial_coeffs_0.z;
    float k4_2 = _S798.radial_coeffs_0.w;
    float p1_3 = _S798.tangential_coeffs_0.x;
    float p2_3 = _S798.tangential_coeffs_0.y;
    float sx1_2 = _S798.thin_prism_coeffs_0.x;
    float sy1_2 = _S798.thin_prism_coeffs_0.y;
    float u_5 = _S803.x;
    float v_5 = _S803.y;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float _S804 = 2.0f * p1_3;
    float _S805 = 2.0f * p2_3;
    float2  _S806 = _S803 * make_float2 (1.0f + r2_5 * (k1_1 + r2_5 * (k2_1 + r2_5 * (k3_1 + r2_5 * k4_2)))) + make_float2 (_S804 * u_5 * v_5 + p2_3 * (r2_5 + 2.0f * u_5 * u_5) + sx1_2 * r2_5, _S805 * u_5 * v_5 + p1_3 * (r2_5 + 2.0f * v_5 * v_5) + sy1_2 * r2_5);
    float2  _S807 = make_float2 (dpfx_0 * _S806.x + dpcx_0, dpfy_0 * _S806.y + dpcy_0);
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (0.0f);
    float _S808 = s_primal_ctx_s_primal_ctx_atan2_0(_S800, _S801);
    bool _S809 = _S808 < 0.00100000004749745f;
    float _S810;
    float _S811;
    float _S812;
    if(_S809)
    {
        float _S813 = 1.0f - _S808 * _S808 / 3.0f;
        float _S814 = _S801 * _S801;
        k_1 = _S813 / _S801;
        _S810 = 0.0f;
        _S811 = _S814;
        _S812 = _S813;
    }
    else
    {
        float _S815 = _S800 * _S800;
        k_1 = _S808 / _S800;
        _S810 = _S815;
        _S811 = 0.0f;
        _S812 = 0.0f;
    }
    float2  _S816 = make_float2 (k_1);
    float2  _S817 = _S799 * make_float2 (k_1);
    float u_6 = _S817.x;
    float v_6 = _S817.y;
    float r2_6 = u_6 * u_6 + v_6 * v_6;
    float _S818 = k3_1 + r2_6 * k4_2;
    float _S819 = k2_1 + r2_6 * _S818;
    float _S820 = k1_1 + r2_6 * _S819;
    float2  _S821 = make_float2 (dpfx_0, 0.0f);
    float2  _S822 = _S817 * _S821;
    float _S823 = p2_3 * dpfx_0;
    float _S824 = _S822.x + _S822.y;
    float _S825 = r2_6 * _S824;
    float _S826 = r2_6 * _S825;
    float _S827 = sx1_2 * dpfx_0 + _S823 + _S820 * _S824 + _S819 * _S825 + _S818 * _S826 + k4_2 * (r2_6 * _S826);
    float _S828 = v_6 * _S827;
    float _S829 = u_6 * _S827;
    float2  _S830 = make_float2 (1.0f + r2_6 * _S820) * _S821 + make_float2 (2.0f * u_6 * _S823 + 2.0f * (u_6 * _S823) + _S804 * (v_6 * dpfx_0) + _S829 + _S829, _S804 * u_6 * dpfx_0 + _S828 + _S828);
    float2  _S831 = _S799 * _S830;
    float2  _S832 = _S816 * _S830;
    float _S833 = _S831.x + _S831.y;
    if(_S809)
    {
        float _S834 = _S833 / _S811;
        float _S835 = _S812 * - _S834;
        float _S836 = _S808 * (0.3333333432674408f * - (_S801 * _S834));
        k_1 = _S836 + _S836;
        _S810 = _S835;
        _S811 = 0.0f;
    }
    else
    {
        float _S837 = _S833 / _S810;
        float _S838 = _S808 * - _S837;
        k_1 = _S800 * _S837;
        _S810 = 0.0f;
        _S811 = _S838;
    }
    DiffPair_float_0 _S839;
    (&_S839)->primal_0 = _S800;
    (&_S839)->differential_0 = 0.0f;
    DiffPair_float_0 _S840;
    (&_S840)->primal_0 = _S801;
    (&_S840)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S839, &_S840, k_1, &_s_diff_ctx_8->_S732);
    _s_diff_ctx_8->_S730 = _S839;
    _s_diff_ctx_8->_S731 = _S840;
    float _S841 = _S840.differential_0 + _S810;
    float _S842 = _S839.differential_0 + _S811;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S843;
    (&_S843)->primal_0 = _S799;
    (&_S843)->differential_0 = _S797;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S844 = { _S797, _S797 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S845;
    (&_S845)->_S758 = _S844;
    s_primal_ctx_s_bwd_length_impl_0(&_S843, _S842, &_S845);
    float2  _S846 = _S843.differential_0 + _S832;
    float3  _S847 = make_float3 (_S846.x, _S846.y, _S841);
    Matrix<float, 2, 3>  _S848 = J_12;
    _S848[int(0)] = _S847;
    if(_S809)
    {
        float _S849 = 1.0f - _S808 * _S808 / 3.0f;
        float _S850 = _S801 * _S801;
        k_1 = _S849 / _S801;
        _S810 = 0.0f;
        _S811 = _S850;
        _S812 = _S849;
    }
    else
    {
        float _S851 = _S800 * _S800;
        k_1 = _S808 / _S800;
        _S810 = _S851;
        _S811 = 0.0f;
        _S812 = 0.0f;
    }
    float2  _S852 = make_float2 (k_1);
    float2  _S853 = _S799 * make_float2 (k_1);
    float u_7 = _S853.x;
    float v_7 = _S853.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S854 = k3_1 + r2_7 * k4_2;
    float _S855 = k2_1 + r2_7 * _S854;
    float _S856 = k1_1 + r2_7 * _S855;
    float2  _S857 = make_float2 (0.0f, dpfy_0);
    float2  _S858 = _S853 * _S857;
    float _S859 = p1_3 * dpfy_0;
    float _S860 = _S858.x + _S858.y;
    float _S861 = r2_7 * _S860;
    float _S862 = r2_7 * _S861;
    float _S863 = sy1_2 * dpfy_0 + _S859 + _S856 * _S860 + _S855 * _S861 + _S854 * _S862 + k4_2 * (r2_7 * _S862);
    float _S864 = v_7 * _S863;
    float _S865 = u_7 * _S863;
    float2  _S866 = make_float2 (1.0f + r2_7 * _S856) * _S857 + make_float2 (_S805 * (v_7 * dpfy_0) + _S865 + _S865, 2.0f * v_7 * _S859 + 2.0f * (v_7 * _S859) + _S805 * u_7 * dpfy_0 + _S864 + _S864);
    float2  _S867 = _S799 * _S866;
    float2  _S868 = _S852 * _S866;
    float _S869 = _S867.x + _S867.y;
    if(_S809)
    {
        float _S870 = _S869 / _S811;
        float _S871 = _S812 * - _S870;
        float _S872 = _S808 * (0.3333333432674408f * - (_S801 * _S870));
        k_1 = _S872 + _S872;
        _S810 = _S871;
        _S811 = 0.0f;
    }
    else
    {
        float _S873 = _S869 / _S810;
        float _S874 = _S808 * - _S873;
        k_1 = _S800 * _S873;
        _S810 = 0.0f;
        _S811 = _S874;
    }
    DiffPair_float_0 _S875;
    (&_S875)->primal_0 = _S800;
    (&_S875)->differential_0 = 0.0f;
    DiffPair_float_0 _S876;
    (&_S876)->primal_0 = _S801;
    (&_S876)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S875, &_S876, k_1, &_s_diff_ctx_8->_S735);
    _s_diff_ctx_8->_S733 = _S875;
    _s_diff_ctx_8->_S734 = _S876;
    float _S877 = _S876.differential_0 + _S810;
    float _S878 = _S875.differential_0 + _S811;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S879;
    (&_S879)->primal_0 = _S799;
    (&_S879)->differential_0 = _S797;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S880;
    (&_S880)->_S758 = _S844;
    s_primal_ctx_s_bwd_length_impl_0(&_S879, _S878, &_S880);
    float2  _S881 = _S879.differential_0 + _S868;
    float3  _S882 = make_float3 (_S881.x, _S881.y, _S877);
    _S848[int(1)] = _S882;
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S848, dpcov3d_0), transpose_1(_S848));
    *dpmean2d_0 = _S807;
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
    DiffPair_1 _S883 = *dpdpx_3;
    float _S884 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_9->_S763)->primal_0);
    float _S885 = s_primal_ctx_sqrt_0(_S884);
    float _S886 = 0.5f / _S885 * (*dpdpx_3).differential_0.differential_0;
    float _S887 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S885 * _S885));
    DiffPair_float_0 _S888;
    (&_S888)->primal_0 = _S884;
    (&_S888)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S888, _S887);
    DiffPair_float_0 _S889;
    (&_S889)->primal_0 = 1.00000001168609742e-07f;
    (&_S889)->differential_0 = 0.0f;
    DiffPair_float_0 _S890;
    (&_S890)->primal_0 = (&_s_diff_ctx_9->_S763)->primal_0;
    (&_S890)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S889, &_S890, _S888.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S890.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S886;
    dpdpx_3->primal_0 = _S883.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S891, DiffPair_float_0 * _S892, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_1 _S893 = *_S891;
    DiffPair_float_0 _S894 = _s_diff_ctx_10->_S759;
    DiffPair_float_0 _S895 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S896;
    (&_S896)->_S763 = _S895;
    s_primal_ctx_d_sqrt_0(&_S894, (*_S892).primal_0, &_S896);
    DiffPair_float_0 _S897 = { (*_S891).differential_0.primal_0, (*_S891).differential_0.differential_0 };
    DiffPair_1 _S898;
    (&_S898)->primal_0 = _s_diff_ctx_10->_S759;
    (&_S898)->differential_0 = _S897;
    DiffPair_float_0 _S899;
    (&_S899)->primal_0 = (*_S892).primal_0;
    (&_S899)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S900 = _S896;
    s_bwd_prop_d_sqrt_0(&_S898, &_S899, &_S900);
    DiffPair_float_0 _S901 = { _S898.differential_0.primal_0, _S898.differential_0.differential_0 };
    _S892->primal_0 = (*_S892).primal_0;
    _S892->differential_0 = _S899.differential_0;
    _S891->primal_0 = _S893.primal_0;
    _S891->differential_0 = _S901;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S902, float _s_dOut_3)
{
    DiffPair_float_0 _S903;
    (&_S903)->primal_0 = (*_S902).primal_0;
    (&_S903)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S903, _s_dOut_3);
    _S902->primal_0 = (*_S902).primal_0;
    _S902->differential_0 = _S903.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_0 _S904 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_11->_S760)->primal_0)->x) * *&((&(&_s_diff_ctx_11->_S760)->primal_0)->x) + *&((&(&_s_diff_ctx_11->_S760)->primal_0)->y) * *&((&(&_s_diff_ctx_11->_S760)->primal_0)->y);
    DiffPair_float_0 _S905 = { len_0, 0.0f };
    float2  _S906 = make_float2 (0.0f);
    float _S907 = (*dpdpx_5).differential_0.differential_0.x;
    float _S908 = _S907 + _S907;
    float _S909 = (&_s_diff_ctx_11->_S761)->differential_0 * _S908;
    float _S910 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S911 = (&_s_diff_ctx_11->_S761)->differential_0 * _S910;
    DiffPair_float_0 _S912 = { 0.0f, *&((&(&_s_diff_ctx_11->_S760)->primal_0)->x) * _S908 + *&((&(&_s_diff_ctx_11->_S760)->primal_0)->y) * _S910 };
    DiffPair_1 _S913;
    (&_S913)->primal_0 = _S905;
    (&_S913)->differential_0 = _S912;
    DiffPair_float_0 _S914;
    (&_S914)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S914)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S913, &_S914, &_s_diff_ctx_11->_S762);
    DiffPair_float_0 _S915;
    (&_S915)->primal_0 = len_0;
    (&_S915)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S915, 0.0f);
    float _S916 = _S913.differential_0.primal_0 + _S915.differential_0;
    float _S917 = *&((&(&_s_diff_ctx_11->_S760)->primal_0)->y) * _S916;
    float _S918 = _S911 + _S917 + _S917;
    float _S919 = *&((&(&_s_diff_ctx_11->_S760)->primal_0)->x) * _S916;
    float _S920 = _S909 + _S919 + _S919;
    float2  _S921 = _S906;
    *&((&_S921)->y) = _S918;
    *&((&_S921)->x) = _S920;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S904.differential_0.primal_0 + _S921, _S906 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S914.differential_0;
    dpdpx_5->primal_0 = _S904.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S922 = (*dpdpx_7).primal_0.x;
    float _S923 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S924;
    (&_S924)->primal_0 = _S922 * _S922 + _S923 * _S923;
    (&_S924)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S924, _s_dOut_4);
    float _S925 = (*dpdpx_7).primal_0.y * _S924.differential_0;
    float _S926 = _S925 + _S925;
    float _S927 = (*dpdpx_7).primal_0.x * _S924.differential_0;
    float _S928 = _S927 + _S927;
    float2  _S929 = make_float2 (0.0f);
    *&((&_S929)->y) = _S926;
    *&((&_S929)->x) = _S928;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S929;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S930, DiffPair_float_0 * _S931, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_0 _S932 = *_S930;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S933 = _s_diff_ctx_12->_S758;
    float2  _S934 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S935 = { _S934, _S934 };
    DiffPair_float_0 _S936 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S937 = { _S936 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S938;
    (&_S938)->_S760 = _S935;
    (&_S938)->_S761 = _S936;
    (&_S938)->_S762 = _S937;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S933, (*_S931).primal_0, &_S938);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S939 = { (*_S930).differential_0.primal_0, (*_S930).differential_0.differential_0 };
    DiffPair_0 _S940;
    (&_S940)->primal_0 = _s_diff_ctx_12->_S758;
    (&_S940)->differential_0 = _S939;
    DiffPair_float_0 _S941;
    (&_S941)->primal_0 = (*_S931).primal_0;
    (&_S941)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S942 = _S938;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S940, &_S941, &_S942);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S943;
    (&_S943)->primal_0 = (&_s_diff_ctx_12->_S758)->primal_0;
    (&_S943)->differential_0 = _S934;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S943, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S944 = { _S940.differential_0.primal_0 + _S943.differential_0, _S940.differential_0.differential_0 };
    _S931->primal_0 = (*_S931).primal_0;
    _S931->differential_0 = _S941.differential_0;
    _S930->primal_0 = _S932.primal_0;
    _S930->differential_0 = _S944;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_1 _S945 = *dpdpy_1;
    DiffPair_1 _S946 = *dpdpx_8;
    float _S947 = - (&_s_diff_ctx_13->_S741)->primal_0;
    float _S948 = (&_s_diff_ctx_13->_S742)->primal_0 * (&_s_diff_ctx_13->_S742)->primal_0 + (&_s_diff_ctx_13->_S741)->primal_0 * (&_s_diff_ctx_13->_S741)->primal_0;
    float _S949 = _S948 * _S948;
    float _S950 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S949;
    float _S951 = (&_s_diff_ctx_13->_S742)->primal_0 * - _S950;
    float _S952 = (&_s_diff_ctx_13->_S741)->primal_0 * _S951;
    float _S953 = (&_s_diff_ctx_13->_S742)->primal_0 * _S951;
    float _S954 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S949;
    float _S955 = _S947 * - _S954;
    float _S956 = (&_s_diff_ctx_13->_S741)->primal_0 * _S955;
    float _S957 = (&_s_diff_ctx_13->_S742)->primal_0 * _S955;
    DiffPair_float_0 dpdpx_9 = { _S957 + _S957 + ((*dpdpx_8).differential_0.primal_0 + (_S953 + _S953 + _S948 * _S950)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S952 + _S952 + (*dpdpy_1).differential_0.primal_0 + _S956 + _S956 + - (_S948 * _S954), 0.0f };
    float _S958 = (&_s_diff_ctx_13->_S742)->primal_0 / _S948 * (*dpdpy_1).differential_0.differential_0 + _S947 / _S948 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S958;
    dpdpy_1->primal_0 = _S945.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S946.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S959, DiffPair_1 * _S960, DiffPair_float_0 * _S961, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_1 _S962 = *_S959;
    DiffPair_1 _S963 = *_S960;
    DiffPair_float_0 _S964 = _s_diff_ctx_14->_S728;
    DiffPair_float_0 _S965 = _s_diff_ctx_14->_S729;
    DiffPair_float_0 _S966 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S967;
    (&_S967)->_S741 = _S966;
    (&_S967)->_S742 = _S966;
    s_primal_ctx_d_atan2_0(&_S964, &_S965, (*_S961).primal_0, &_S967);
    DiffPair_float_0 _S968 = { (*_S960).differential_0.primal_0, (*_S960).differential_0.differential_0 };
    DiffPair_float_0 _S969 = { (*_S959).differential_0.primal_0, (*_S959).differential_0.differential_0 };
    DiffPair_1 _S970;
    (&_S970)->primal_0 = _s_diff_ctx_14->_S728;
    (&_S970)->differential_0 = _S969;
    DiffPair_1 _S971;
    (&_S971)->primal_0 = _s_diff_ctx_14->_S729;
    (&_S971)->differential_0 = _S968;
    DiffPair_float_0 _S972;
    (&_S972)->primal_0 = (*_S961).primal_0;
    (&_S972)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S973 = _S967;
    s_bwd_prop_d_atan2_0(&_S970, &_S971, &_S972, &_S973);
    DiffPair_float_0 _S974 = { _S971.differential_0.primal_0, _S971.differential_0.differential_0 };
    DiffPair_float_0 _S975 = { _S970.differential_0.primal_0, _S970.differential_0.differential_0 };
    _S961->primal_0 = (*_S961).primal_0;
    _S961->differential_0 = _S972.differential_0;
    _S959->primal_0 = _S962.primal_0;
    _S959->differential_0 = _S975;
    _S960->primal_0 = _S963.primal_0;
    _S960->differential_0 = _S974;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S976, DiffPair_float_0 * _S977, float _s_dOut_5)
{
    DiffPair_float_0 _S978;
    (&_S978)->primal_0 = (*_S976).primal_0;
    (&_S978)->differential_0 = 0.0f;
    DiffPair_float_0 _S979;
    (&_S979)->primal_0 = (*_S977).primal_0;
    (&_S979)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S978, &_S979, _s_dOut_5);
    _S977->primal_0 = (*_S977).primal_0;
    _S977->differential_0 = _S979.differential_0;
    _S976->primal_0 = (*_S976).primal_0;
    _S976->differential_0 = _S978.differential_0;
    return;
}

inline __device__ void s_bwd_prop_CameraDistortion_x24init_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_2, CameraDistortion_0 * _s_dOut_6)
{
    float2  _S980 = _s_dOut_6->thin_prism_coeffs_0;
    dpthin_prism_coeffs_2->primal_0 = (*dpthin_prism_coeffs_2).primal_0;
    dpthin_prism_coeffs_2->differential_0 = _S980;
    float2  _S981 = _s_dOut_6->tangential_coeffs_0;
    dptangential_coeffs_2->primal_0 = (*dptangential_coeffs_2).primal_0;
    dptangential_coeffs_2->differential_0 = _S981;
    float4  _S982 = _s_dOut_6->radial_coeffs_0;
    dpradial_coeffs_2->primal_0 = (*dpradial_coeffs_2).primal_0;
    dpradial_coeffs_2->differential_0 = _S982;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_3, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S983 = *dpcov3d_1;
    DiffPair_float_0 _S984 = *dpfx_1;
    DiffPair_float_0 _S985 = *dpfy_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S986 = *dpradial_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S987 = *dptangential_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S988 = *dpthin_prism_coeffs_3;
    float2  _S989 = make_float2 (0.0f);
    CameraDistortion_0 _S990 = s_primal_ctx_CameraDistortion_x24init_0((*dpradial_coeffs_3).primal_0, (*dptangential_coeffs_3).primal_0, (*dpthin_prism_coeffs_3).primal_0);
    float2  _S991 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S992 = length_0(_S991);
    float _S993 = (*dpmean3d_1).primal_0.z;
    float _S994 = s_primal_ctx_atan2_0(_S992, _S993);
    bool _S995 = _S994 < 0.00100000004749745f;
    float k_2;
    float _S996;
    float _S997;
    float _S998;
    if(_S995)
    {
        float _S999 = 1.0f - _S994 * _S994 / 3.0f;
        float _S1000 = _S993 * _S993;
        k_2 = _S999 / _S993;
        _S996 = 0.0f;
        _S997 = _S1000;
        _S998 = _S999;
    }
    else
    {
        float _S1001 = _S992 * _S992;
        k_2 = _S994 / _S992;
        _S996 = _S1001;
        _S997 = 0.0f;
        _S998 = 0.0f;
    }
    float2  _S1002 = make_float2 (k_2);
    float2  _S1003 = _S991 * make_float2 (k_2);
    float k1_2 = _S990.radial_coeffs_0.x;
    float k2_2 = _S990.radial_coeffs_0.y;
    float k3_2 = _S990.radial_coeffs_0.z;
    float k4_3 = _S990.radial_coeffs_0.w;
    float p1_4 = _S990.tangential_coeffs_0.x;
    float p2_4 = _S990.tangential_coeffs_0.y;
    float sx1_3 = _S990.thin_prism_coeffs_0.x;
    float sy1_3 = _S990.thin_prism_coeffs_0.y;
    float u_8 = _S1003.x;
    float v_8 = _S1003.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S1004 = k3_2 + r2_8 * k4_3;
    float _S1005 = k2_2 + r2_8 * _S1004;
    float _S1006 = k1_2 + r2_8 * _S1005;
    float radial_1 = 1.0f + r2_8 * _S1006;
    float2  _S1007 = make_float2 (radial_1);
    float _S1008 = 2.0f * p1_4;
    float _S1009 = _S1008 * u_8;
    float _S1010 = 2.0f * u_8;
    float _S1011 = r2_8 + _S1010 * u_8;
    float _S1012 = 2.0f * p2_4;
    float _S1013 = _S1012 * u_8;
    float _S1014 = 2.0f * v_8;
    float _S1015 = r2_8 + _S1014 * v_8;
    float2  _S1016 = _S1003 * make_float2 (radial_1) + make_float2 (_S1009 * v_8 + p2_4 * _S1011 + sx1_3 * r2_8, _S1013 * v_8 + p1_4 * _S1015 + sy1_3 * r2_8);
    float _S1017 = _S1016.x;
    float _S1018 = _S1016.y;
    Matrix<float, 2, 3>  J_13 = makeMatrix<float, 2, 3> (0.0f);
    float _S1019 = s_primal_ctx_s_primal_ctx_atan2_0(_S992, _S993);
    bool _S1020 = _S1019 < 0.00100000004749745f;
    float _S1021;
    float _S1022;
    float _S1023;
    if(_S1020)
    {
        float _S1024 = 1.0f - _S1019 * _S1019 / 3.0f;
        float _S1025 = _S993 * _S993;
        k_2 = _S1024 / _S993;
        _S1021 = 0.0f;
        _S1022 = _S1025;
        _S1023 = _S1024;
    }
    else
    {
        float _S1026 = _S992 * _S992;
        k_2 = _S1019 / _S992;
        _S1021 = _S1026;
        _S1022 = 0.0f;
        _S1023 = 0.0f;
    }
    float2  _S1027 = make_float2 (k_2);
    float2  _S1028 = _S991 * make_float2 (k_2);
    float u_9 = _S1028.x;
    float v_9 = _S1028.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S1029 = k3_2 + r2_9 * k4_3;
    float _S1030 = k2_2 + r2_9 * _S1029;
    float _S1031 = k1_2 + r2_9 * _S1030;
    float2  _S1032 = make_float2 (1.0f + r2_9 * _S1031);
    float _S1033 = _S1008 * u_9;
    float _S1034 = 2.0f * u_9;
    float2  _S1035 = make_float2 (_S984.primal_0, 0.0f);
    float2  _S1036 = _S1028 * _S1035;
    float _S1037 = p2_4 * _S984.primal_0;
    float _S1038 = v_9 * _S984.primal_0;
    float _S1039 = _S1036.x + _S1036.y;
    float _S1040 = r2_9 * _S1039;
    float _S1041 = r2_9 * _S1040;
    float _S1042 = r2_9 * _S1041;
    float _S1043 = sx1_3 * _S984.primal_0 + _S1037 + _S1031 * _S1039 + _S1030 * _S1040 + _S1029 * _S1041 + k4_3 * _S1042;
    float _S1044 = v_9 * _S1043;
    float _S1045 = u_9 * _S1043;
    float2  _S1046 = _S1032 * _S1035 + make_float2 (_S1034 * _S1037 + 2.0f * (u_9 * _S1037) + _S1008 * _S1038 + _S1045 + _S1045, _S1033 * _S984.primal_0 + _S1044 + _S1044);
    float2  _S1047 = _S991 * _S1046;
    float2  _S1048 = _S1027 * _S1046;
    float _S1049 = _S1047.x + _S1047.y;
    float k_3;
    float _S1050;
    float _S1051;
    float _S1052;
    float _S1053;
    float _S1054;
    float _S1055;
    float _S1056;
    float _S1057;
    if(_S1020)
    {
        float _S1058 = _S1049 / _S1022;
        float _S1059 = _S1022 * _S1022;
        float _S1060 = - _S1058;
        float _S1061 = _S1023 * _S1060;
        float _S1062 = 0.3333333432674408f * - (_S993 * _S1058);
        float _S1063 = _S1019 * _S1062;
        k_2 = _S1063 + _S1063;
        k_3 = _S1061;
        _S1050 = 0.0f;
        _S1051 = 0.0f;
        _S1052 = 0.0f;
        _S1053 = 0.0f;
        _S1054 = _S1062;
        _S1055 = _S1058;
        _S1056 = _S1060;
        _S1057 = _S1059;
    }
    else
    {
        float _S1064 = _S1049 / _S1021;
        float _S1065 = _S1021 * _S1021;
        float _S1066 = - _S1064;
        float _S1067 = _S1019 * _S1066;
        k_2 = _S992 * _S1064;
        k_3 = 0.0f;
        _S1050 = _S1067;
        _S1051 = _S1064;
        _S1052 = _S1066;
        _S1053 = _S1065;
        _S1054 = 0.0f;
        _S1055 = 0.0f;
        _S1056 = 0.0f;
        _S1057 = 0.0f;
    }
    DiffPair_float_0 _S1068 = { _S992, 0.0f };
    DiffPair_float_0 _S1069 = { _S993, 0.0f };
    float _S1070 = (&_s_diff_ctx_15->_S731)->differential_0 + k_3;
    float _S1071 = (&_s_diff_ctx_15->_S730)->differential_0 + _S1050;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1072 = { _S991, _S989 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1073;
    (&_S1073)->primal_0 = _S991;
    (&_S1073)->differential_0 = _S989;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1074 = { _S989, _S989 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1075;
    (&_S1075)->_S758 = _S1074;
    s_primal_ctx_s_bwd_length_impl_0(&_S1073, _S1071, &_S1075);
    float2  _S1076 = _S1073.differential_0 + _S1048;
    float3  _S1077 = make_float3 (_S1076.x, _S1076.y, _S1070);
    Matrix<float, 2, 3>  _S1078 = J_13;
    _S1078[int(0)] = _S1077;
    float _S1079;
    float _S1080;
    if(_S1020)
    {
        float _S1081 = 1.0f - _S1019 * _S1019 / 3.0f;
        float _S1082 = _S993 * _S993;
        k_3 = _S1081 / _S993;
        _S1050 = 0.0f;
        _S1079 = _S1082;
        _S1080 = _S1081;
    }
    else
    {
        float _S1083 = _S992 * _S992;
        k_3 = _S1019 / _S992;
        _S1050 = _S1083;
        _S1079 = 0.0f;
        _S1080 = 0.0f;
    }
    float2  _S1084 = make_float2 (k_3);
    float2  _S1085 = _S991 * make_float2 (k_3);
    float u_10 = _S1085.x;
    float v_10 = _S1085.y;
    float r2_10 = u_10 * u_10 + v_10 * v_10;
    float _S1086 = k3_2 + r2_10 * k4_3;
    float _S1087 = k2_2 + r2_10 * _S1086;
    float _S1088 = k1_2 + r2_10 * _S1087;
    float2  _S1089 = make_float2 (1.0f + r2_10 * _S1088);
    float _S1090 = _S1012 * u_10;
    float _S1091 = 2.0f * v_10;
    float2  _S1092 = make_float2 (0.0f, _S985.primal_0);
    float2  _S1093 = _S1085 * _S1092;
    float _S1094 = p1_4 * _S985.primal_0;
    float _S1095 = v_10 * _S985.primal_0;
    float _S1096 = _S1093.x + _S1093.y;
    float _S1097 = r2_10 * _S1096;
    float _S1098 = r2_10 * _S1097;
    float _S1099 = r2_10 * _S1098;
    float _S1100 = sy1_3 * _S985.primal_0 + _S1094 + _S1088 * _S1096 + _S1087 * _S1097 + _S1086 * _S1098 + k4_3 * _S1099;
    float _S1101 = v_10 * _S1100;
    float _S1102 = u_10 * _S1100;
    float2  _S1103 = _S1089 * _S1092 + make_float2 (_S1012 * _S1095 + _S1102 + _S1102, _S1091 * _S1094 + 2.0f * (v_10 * _S1094) + _S1090 * _S985.primal_0 + _S1101 + _S1101);
    float2  _S1104 = _S991 * _S1103;
    float2  _S1105 = _S1084 * _S1103;
    float _S1106 = _S1104.x + _S1104.y;
    float _S1107;
    float _S1108;
    float _S1109;
    float _S1110;
    float _S1111;
    float _S1112;
    float _S1113;
    float _S1114;
    float _S1115;
    if(_S1020)
    {
        float _S1116 = _S1106 / _S1079;
        float _S1117 = _S1079 * _S1079;
        float _S1118 = - _S1116;
        float _S1119 = _S1080 * _S1118;
        float _S1120 = 0.3333333432674408f * - (_S993 * _S1116);
        float _S1121 = _S1019 * _S1120;
        k_3 = _S1121 + _S1121;
        _S1107 = _S1119;
        _S1108 = 0.0f;
        _S1109 = 0.0f;
        _S1110 = 0.0f;
        _S1111 = 0.0f;
        _S1112 = _S1120;
        _S1113 = _S1116;
        _S1114 = _S1118;
        _S1115 = _S1117;
    }
    else
    {
        float _S1122 = _S1106 / _S1050;
        float _S1123 = _S1050 * _S1050;
        float _S1124 = - _S1122;
        float _S1125 = _S1019 * _S1124;
        k_3 = _S992 * _S1122;
        _S1107 = 0.0f;
        _S1108 = _S1125;
        _S1109 = _S1122;
        _S1110 = _S1124;
        _S1111 = _S1123;
        _S1112 = 0.0f;
        _S1113 = 0.0f;
        _S1114 = 0.0f;
        _S1115 = 0.0f;
    }
    float _S1126 = (&_s_diff_ctx_15->_S734)->differential_0 + _S1107;
    float _S1127 = (&_s_diff_ctx_15->_S733)->differential_0 + _S1108;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1128;
    (&_S1128)->primal_0 = _S991;
    (&_S1128)->differential_0 = _S989;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1129;
    (&_S1129)->_S758 = _S1074;
    s_primal_ctx_s_bwd_length_impl_0(&_S1128, _S1127, &_S1129);
    float2  _S1130 = _S1128.differential_0 + _S1105;
    float3  _S1131 = make_float3 (_S1130.x, _S1130.y, _S1126);
    _S1078[int(1)] = _S1131;
    Matrix<float, 3, 2>  _S1132 = transpose_1(_S1078);
    CameraDistortion_0 _S1133 = CameraDistortion_x24_syn_dzero_0();
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1134;
    (&_S1134)->primal_0 = s_primal_ctx_mul_3(_S1078, _S983.primal_0);
    (&_S1134)->differential_0 = J_13;
    Matrix<float, 3, 2>  _S1135 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1136;
    (&_S1136)->primal_0 = _S1132;
    (&_S1136)->differential_0 = _S1135;
    s_bwd_prop_mul_2(&_S1134, &_S1136, dpcov2d_1);
    Matrix<float, 2, 3>  _S1137 = transpose_2(_S1136.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1138;
    (&_S1138)->primal_0 = _S1078;
    (&_S1138)->differential_0 = J_13;
    Matrix<float, 3, 3>  _S1139 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1140;
    (&_S1140)->primal_0 = _S983.primal_0;
    (&_S1140)->differential_0 = _S1139;
    s_bwd_prop_mul_3(&_S1138, &_S1140, _S1134.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1141 = _S1140;
    Matrix<float, 2, 3>  _S1142 = _S1137 + _S1138.differential_0;
    float2  _S1143 = _S989;
    *&((&_S1143)->y) = _S1142.rows[int(1)].y;
    *&((&_S1143)->x) = _S1142.rows[int(1)].x;
    float2  _S1144 = _S1143;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1145 = { _S989, _S1143 };
    DiffPair_0 _S1146;
    (&_S1146)->primal_0 = _S1072;
    (&_S1146)->differential_0 = _S1145;
    DiffPair_float_0 _S1147;
    (&_S1147)->primal_0 = _S1127;
    (&_S1147)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1148 = _S1129;
    s_bwd_prop_s_bwd_length_impl_0(&_S1146, &_S1147, &_S1148);
    DiffPair_0 _S1149 = _S1146;
    DiffPair_float_0 _S1150 = _S1147;
    DiffPair_float_0 _S1151 = { 0.0f, _S1142.rows[int(1)].z };
    DiffPair_float_0 _S1152 = { 0.0f, _S1147.differential_0 };
    DiffPair_1 _S1153;
    (&_S1153)->primal_0 = _S1068;
    (&_S1153)->differential_0 = _S1152;
    DiffPair_1 _S1154;
    (&_S1154)->primal_0 = _S1069;
    (&_S1154)->differential_0 = _S1151;
    DiffPair_float_0 _S1155;
    (&_S1155)->primal_0 = k_3;
    (&_S1155)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1153, &_S1154, &_S1155, &_s_diff_ctx_15->_S735);
    DiffPair_1 _S1156 = _S1153;
    DiffPair_1 _S1157 = _S1154;
    DiffPair_float_0 _S1158 = _S1155;
    if(_S1020)
    {
        float _S1159 = _S1158.differential_0 + _S1158.differential_0;
        float _S1160 = _S1112 * _S1159;
        float _S1161 = - (0.3333333432674408f * (_S1019 * _S1159));
        float _S1162 = _S1114 * _S1142.rows[int(1)].z;
        float _S1163 = (_S993 * _S1161 + - (_S1080 * _S1142.rows[int(1)].z)) / _S1115;
        float _S1164 = _S1106 * - _S1163;
        float _S1165 = _S1113 * _S1161 + _S1157.differential_0.primal_0;
        k_3 = _S1079 * _S1163;
        _S1107 = 0.0f;
        _S1108 = _S1164;
        _S1109 = _S1162;
        _S1110 = _S1156.differential_0.primal_0;
        _S1111 = _S1160;
        _S1112 = _S1165;
    }
    else
    {
        float _S1166 = _S1110 * _S1150.differential_0;
        float _S1167 = (_S992 * _S1158.differential_0 + - (_S1019 * _S1150.differential_0)) / _S1111;
        float _S1168 = _S1106 * - _S1167;
        float _S1169 = _S1109 * _S1158.differential_0 + _S1156.differential_0.primal_0;
        k_3 = _S1050 * _S1167;
        _S1107 = _S1168;
        _S1108 = 0.0f;
        _S1109 = 0.0f;
        _S1110 = _S1169;
        _S1111 = _S1166;
        _S1112 = _S1157.differential_0.primal_0;
    }
    float2  _S1170 = _S1084 * _S1144;
    float2  _S1171 = _S1103 * _S1144;
    float2  _S1172 = _S989;
    *&((&_S1172)->y) = k_3;
    *&((&_S1172)->x) = k_3;
    float2  _S1173 = _S1103 * _S1172;
    float2  _S1174 = _S1170 + _S991 * _S1172;
    float _S1175 = _S1174.x;
    float _S1176 = _S1175 + _S1175;
    float _S1177 = _S1100 * _S1176;
    float _S1178 = _S1174.y + _S1174.y;
    float _S1179 = _S1100 * _S1178;
    float _S1180 = u_10 * _S1176 + v_10 * _S1178;
    float _S1181 = k4_3 * _S1180;
    float _S1182 = _S1099 * _S1180;
    float _S1183 = _S1098 * _S1180;
    float _S1184 = _S1098 * _S1181;
    float _S1185 = _S1097 * _S1180;
    float _S1186 = _S1086 * _S1180 + r2_10 * _S1181;
    float _S1187 = _S1097 * _S1186;
    float _S1188 = _S1096 * _S1180;
    float _S1189 = _S1087 * _S1180 + r2_10 * _S1186;
    float _S1190 = _S1096 * _S1189;
    float _S1191 = _S1088 * _S1180 + r2_10 * _S1189;
    float _S1192 = _S1012 * _S1174.x;
    float _S1193 = _S1095 * _S1174.x;
    float _S1194 = v_10 * _S1192;
    float _S1195 = _S985.primal_0 * _S1192;
    float _S1196 = _S1090 * _S1174.y;
    float _S1197 = _S985.primal_0 * _S1174.y;
    float _S1198 = 2.0f * _S1174.y;
    float _S1199 = _S1094 * _S1198;
    float _S1200 = _S1094 * _S1174.y;
    float _S1201 = _S1180 + v_10 * _S1198 + _S1091 * _S1174.y;
    float _S1202 = p1_4 * _S1201;
    float _S1203 = _S985.primal_0 * _S1201;
    float _S1204 = sy1_3 * _S1180;
    float _S1205 = _S985.primal_0 * _S1180;
    float2  _S1206 = _S1089 * _S1174;
    float2  _S1207 = _S1092 * _S1174;
    float2  _S1208 = _S989;
    *&((&_S1208)->y) = _S1191;
    *&((&_S1208)->x) = _S1191;
    float _S1209 = _S1207.x + _S1207.y;
    float _S1210 = _S1188 + r2_10 * _S1209;
    float _S1211 = _S1185 + r2_10 * _S1210;
    float _S1212 = _S1183 + r2_10 * _S1211;
    float _S1213 = _S1184 + _S1187 + _S1190 + _S1088 * _S1209 + _S1087 * _S1210 + _S1086 * _S1211 + k4_3 * _S1212;
    float _S1214 = v_10 * _S1213;
    float _S1215 = u_10 * _S1213;
    float2  _S1216 = _S1173 + _S1149.differential_0.primal_0;
    float2  _S1217 = _S1092 * _S1208 + make_float2 (_S1177 + _S1012 * _S1197 + _S1215 + _S1215, _S1179 + _S1195 + _S1199 + 2.0f * _S1200 + _S1214 + _S1214);
    float _S1218 = _S1194 + _S1196 + _S1202 + _S1204 + (_S1206 + _S1085 * _S1208).y;
    float _S1219 = _S1193 + u_10 * _S1197;
    float _S1220 = _S1182 + r2_10 * _S1212;
    float2  _S1221 = _S991 * _S1217;
    float _S1222 = _S1171.x + _S1171.y + _S1221.x + _S1221.y;
    float2  _S1223 = _S1084 * _S1217 + _S1216;
    if(_S1020)
    {
        float _S1224 = _S993 * _S1108;
        float _S1225 = _S1222 / _S1079;
        float _S1226 = _S1019 * (0.3333333432674408f * - (_S1109 + _S993 * _S1225));
        float _S1227 = _S1224 + _S1224 + _S1080 * - _S1225 + _S1112;
        k_3 = _S1226 + _S1226 + _S1111;
        _S1050 = _S1227;
        _S1079 = _S1110;
    }
    else
    {
        float _S1228 = _S992 * _S1107;
        float _S1229 = _S1222 / _S1050;
        float _S1230 = _S1228 + _S1228 + _S1019 * - _S1229 + _S1110;
        k_3 = _S992 * _S1229 + _S1111;
        _S1050 = _S1112;
        _S1079 = _S1230;
    }
    DiffPair_float_0 _S1231;
    (&_S1231)->primal_0 = _S992;
    (&_S1231)->differential_0 = 0.0f;
    DiffPair_float_0 _S1232;
    (&_S1232)->primal_0 = _S993;
    (&_S1232)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1231, &_S1232, k_3);
    float _S1233 = _S1232.differential_0 + _S1050;
    float _S1234 = _S1231.differential_0 + _S1079;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1235;
    (&_S1235)->primal_0 = _S991;
    (&_S1235)->differential_0 = _S989;
    s_bwd_length_impl_1(&_S1235, _S1234);
    float2  _S1236 = _S1235.differential_0 + _S1223;
    float3  _S1237 = make_float3 (_S1236.x, _S1236.y, _S1233);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1238;
    (&_S1238)->primal_0 = _S991;
    (&_S1238)->differential_0 = _S989;
    s_bwd_length_impl_1(&_S1238, 0.0f);
    float3  _S1239 = _S1237 + make_float3 (_S1238.differential_0.x, _S1238.differential_0.y, 0.0f);
    float2  _S1240 = _S989;
    *&((&_S1240)->y) = _S1142.rows[int(0)].y;
    *&((&_S1240)->x) = _S1142.rows[int(0)].x;
    float2  _S1241 = _S1240;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1242 = { _S989, _S1240 };
    DiffPair_0 _S1243;
    (&_S1243)->primal_0 = _S1072;
    (&_S1243)->differential_0 = _S1242;
    DiffPair_float_0 _S1244;
    (&_S1244)->primal_0 = _S1071;
    (&_S1244)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1245 = _S1075;
    s_bwd_prop_s_bwd_length_impl_0(&_S1243, &_S1244, &_S1245);
    DiffPair_0 _S1246 = _S1243;
    DiffPair_float_0 _S1247 = _S1244;
    DiffPair_float_0 _S1248 = { 0.0f, _S1142.rows[int(0)].z };
    DiffPair_float_0 _S1249 = { 0.0f, _S1244.differential_0 };
    DiffPair_1 _S1250;
    (&_S1250)->primal_0 = _S1068;
    (&_S1250)->differential_0 = _S1249;
    DiffPair_1 _S1251;
    (&_S1251)->primal_0 = _S1069;
    (&_S1251)->differential_0 = _S1248;
    DiffPair_float_0 _S1252;
    (&_S1252)->primal_0 = k_2;
    (&_S1252)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1250, &_S1251, &_S1252, &_s_diff_ctx_15->_S732);
    DiffPair_1 _S1253 = _S1250;
    DiffPair_1 _S1254 = _S1251;
    DiffPair_float_0 _S1255 = _S1252;
    if(_S1020)
    {
        float _S1256 = _S1255.differential_0 + _S1255.differential_0;
        float _S1257 = _S1054 * _S1256;
        float _S1258 = - (0.3333333432674408f * (_S1019 * _S1256));
        float _S1259 = _S1056 * _S1142.rows[int(0)].z;
        float _S1260 = (_S993 * _S1258 + - (_S1023 * _S1142.rows[int(0)].z)) / _S1057;
        float _S1261 = _S1049 * - _S1260;
        float _S1262 = _S1055 * _S1258 + _S1254.differential_0.primal_0;
        k_2 = _S1022 * _S1260;
        k_3 = 0.0f;
        _S1050 = _S1261;
        _S1051 = _S1259;
        _S1052 = _S1253.differential_0.primal_0;
        _S1053 = _S1257;
        _S1054 = _S1262;
    }
    else
    {
        float _S1263 = _S1052 * _S1247.differential_0;
        float _S1264 = (_S992 * _S1255.differential_0 + - (_S1019 * _S1247.differential_0)) / _S1053;
        float _S1265 = _S1049 * - _S1264;
        float _S1266 = _S1051 * _S1255.differential_0 + _S1253.differential_0.primal_0;
        k_2 = _S1021 * _S1264;
        k_3 = _S1265;
        _S1050 = 0.0f;
        _S1051 = 0.0f;
        _S1052 = _S1266;
        _S1053 = _S1263;
        _S1054 = _S1254.differential_0.primal_0;
    }
    float2  _S1267 = _S1027 * _S1241;
    float2  _S1268 = _S1046 * _S1241;
    float2  _S1269 = _S989;
    *&((&_S1269)->y) = k_2;
    *&((&_S1269)->x) = k_2;
    float2  _S1270 = _S1046 * _S1269;
    float2  _S1271 = _S1267 + _S991 * _S1269;
    float _S1272 = _S1271.x;
    float _S1273 = _S1272 + _S1272;
    float _S1274 = _S1043 * _S1273;
    float _S1275 = _S1271.y + _S1271.y;
    float _S1276 = _S1043 * _S1275;
    float _S1277 = u_9 * _S1273 + v_9 * _S1275;
    float _S1278 = k4_3 * _S1277;
    float _S1279 = _S1042 * _S1277;
    float _S1280 = _S1041 * _S1277;
    float _S1281 = _S1041 * _S1278;
    float _S1282 = _S1040 * _S1277;
    float _S1283 = _S1029 * _S1277 + r2_9 * _S1278;
    float _S1284 = _S1040 * _S1283;
    float _S1285 = _S1039 * _S1277;
    float _S1286 = _S1030 * _S1277 + r2_9 * _S1283;
    float _S1287 = _S1039 * _S1286;
    float _S1288 = _S1031 * _S1277 + r2_9 * _S1286;
    float _S1289 = _S1008 * _S1271.x;
    float _S1290 = _S1038 * _S1271.x;
    float _S1291 = v_9 * _S1289;
    float _S1292 = _S984.primal_0 * _S1289;
    float _S1293 = _S1033 * _S1271.y;
    float _S1294 = _S984.primal_0 * _S1271.y;
    float _S1295 = 2.0f * _S1271.x;
    float _S1296 = _S1037 * _S1295;
    float _S1297 = _S1037 * _S1271.x;
    float _S1298 = _S1277 + u_9 * _S1295 + _S1034 * _S1271.x;
    float _S1299 = p2_4 * _S1298;
    float _S1300 = _S984.primal_0 * _S1298;
    float _S1301 = sx1_3 * _S1277;
    float _S1302 = _S984.primal_0 * _S1277;
    float2  _S1303 = _S1032 * _S1271;
    float2  _S1304 = _S1035 * _S1271;
    float2  _S1305 = _S989;
    *&((&_S1305)->y) = _S1288;
    *&((&_S1305)->x) = _S1288;
    float _S1306 = _S1304.x + _S1304.y;
    float _S1307 = _S1285 + r2_9 * _S1306;
    float _S1308 = _S1282 + r2_9 * _S1307;
    float _S1309 = _S1280 + r2_9 * _S1308;
    float _S1310 = _S1281 + _S1284 + _S1287 + _S1031 * _S1306 + _S1030 * _S1307 + _S1029 * _S1308 + k4_3 * _S1309;
    float _S1311 = v_9 * _S1310;
    float _S1312 = u_9 * _S1310;
    float2  _S1313 = _S1270 + _S1246.differential_0.primal_0;
    float _S1314 = _S1290 + u_9 * _S1294;
    float _S1315 = _S1307 + _S1210;
    float _S1316 = _S1309 + _S1212;
    float _S1317 = _S1291 + _S1293 + _S1299 + _S1301 + (_S1303 + _S1028 * _S1305).x;
    float2  _S1318 = _S1035 * _S1305 + make_float2 (_S1274 + _S1296 + 2.0f * _S1297 + _S1008 * _S1294 + _S1312 + _S1312, _S1276 + _S1292 + _S1311 + _S1311);
    float _S1319 = _S1308 + _S1211;
    float _S1320 = _S1279 + r2_9 * _S1309 + _S1220;
    float2  _S1321 = _S991 * _S1318;
    float _S1322 = _S1268.x + _S1268.y + _S1321.x + _S1321.y;
    float2  _S1323 = _S1027 * _S1318 + _S1313;
    if(_S1020)
    {
        float _S1324 = _S993 * _S1050;
        float _S1325 = _S1322 / _S1022;
        float _S1326 = _S1019 * (0.3333333432674408f * - (_S1051 + _S993 * _S1325));
        float _S1327 = _S1324 + _S1324 + _S1023 * - _S1325 + _S1054;
        k_2 = _S1326 + _S1326 + _S1053;
        _S1021 = _S1327;
        _S1022 = _S1052;
    }
    else
    {
        float _S1328 = _S992 * k_3;
        float _S1329 = _S1322 / _S1021;
        float _S1330 = _S1328 + _S1328 + _S1019 * - _S1329 + _S1052;
        k_2 = _S992 * _S1329 + _S1053;
        _S1021 = _S1054;
        _S1022 = _S1330;
    }
    DiffPair_float_0 _S1331;
    (&_S1331)->primal_0 = _S992;
    (&_S1331)->differential_0 = 0.0f;
    DiffPair_float_0 _S1332;
    (&_S1332)->primal_0 = _S993;
    (&_S1332)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1331, &_S1332, k_2);
    float _S1333 = _S1332.differential_0 + _S1021;
    float _S1334 = _S1331.differential_0 + _S1022;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1335;
    (&_S1335)->primal_0 = _S991;
    (&_S1335)->differential_0 = _S989;
    s_bwd_length_impl_1(&_S1335, _S1334);
    float2  _S1336 = _S1335.differential_0 + _S1323;
    float3  _S1337 = make_float3 (_S1336.x, _S1336.y, _S1333);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1338;
    (&_S1338)->primal_0 = _S991;
    (&_S1338)->differential_0 = _S989;
    s_bwd_length_impl_1(&_S1338, 0.0f);
    float _S1339 = _S985.primal_0 * dpmean2d_1.y;
    float _S1340 = _S984.primal_0 * dpmean2d_1.x;
    float2  _S1341 = make_float2 (_S1340, _S1339);
    float2  _S1342 = _S1003 * _S1341;
    float _S1343 = p1_4 * _S1339;
    float _S1344 = v_8 * _S1339;
    float _S1345 = p2_4 * _S1340;
    float _S1346 = v_8 * _S1340;
    float _S1347 = _S1342.x + _S1342.y;
    float _S1348 = r2_8 * _S1347;
    float _S1349 = r2_8 * _S1348;
    float _S1350 = r2_8 * _S1349;
    float _S1351 = sy1_3 * _S1339 + _S1343 + sx1_3 * _S1340 + _S1345 + _S1006 * _S1347 + _S1005 * _S1348 + _S1004 * _S1349 + k4_3 * _S1350;
    float _S1352 = v_8 * _S1351;
    float _S1353 = u_8 * _S1351;
    float2  _S1354 = make_float2 (r2_8 * _S1340 + _S1302, r2_8 * _S1339 + _S1205);
    float2  _S1355 = make_float2 (_S1015 * _S1339 + 2.0f * (u_8 * _S1346 + _S1314) + _S1203, 2.0f * (u_8 * _S1344 + _S1219) + _S1011 * _S1340 + _S1300);
    float4  _S1356 = make_float4 (_S1348 + _S1315, _S1349 + _S1319, _S1350 + _S1316, r2_8 * _S1350 + _S1320);
    float3  _S1357 = _S1337 + make_float3 (_S1338.differential_0.x, _S1338.differential_0.y, 0.0f) + _S1239;
    float _S1358 = _S1017 * dpmean2d_1.x + _S1317;
    float _S1359 = _S1018 * dpmean2d_1.y + _S1218;
    float2  _S1360 = _S1007 * _S1341 + make_float2 (_S1012 * _S1344 + _S1010 * _S1345 + 2.0f * (u_8 * _S1345) + _S1008 * _S1346 + _S1353 + _S1353, _S1014 * _S1343 + 2.0f * (v_8 * _S1343) + _S1013 * _S1339 + _S1009 * _S1340 + _S1352 + _S1352);
    CameraDistortion_0 _S1361 = _S1133;
    (&_S1361)->thin_prism_coeffs_0 = _S1354;
    (&_S1361)->tangential_coeffs_0 = _S1355;
    (&_S1361)->radial_coeffs_0 = _S1356;
    CameraDistortion_0 _S1362 = _S1133;
    CameraDistortion_0 _S1363 = _S1361;
    CameraDistortion_0 _S1364 = CameraDistortion_x24_syn_dadd_0(&_S1362, &_S1363);
    float2  _S1365 = _S991 * _S1360;
    float2  _S1366 = _S1002 * _S1360;
    float _S1367 = _S1365.x + _S1365.y;
    if(_S995)
    {
        float _S1368 = _S1367 / _S997;
        float _S1369 = _S998 * - _S1368;
        float _S1370 = _S994 * (0.3333333432674408f * - (_S993 * _S1368));
        k_2 = _S1370 + _S1370;
        _S996 = _S1369;
        _S997 = 0.0f;
    }
    else
    {
        float _S1371 = _S1367 / _S996;
        float _S1372 = _S994 * - _S1371;
        k_2 = _S992 * _S1371;
        _S996 = 0.0f;
        _S997 = _S1372;
    }
    DiffPair_float_0 _S1373;
    (&_S1373)->primal_0 = _S992;
    (&_S1373)->differential_0 = 0.0f;
    DiffPair_float_0 _S1374;
    (&_S1374)->primal_0 = _S993;
    (&_S1374)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1373, &_S1374, k_2);
    float _S1375 = _S1374.differential_0 + _S996;
    float _S1376 = _S1373.differential_0 + _S997;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1377;
    (&_S1377)->primal_0 = _S991;
    (&_S1377)->differential_0 = _S989;
    s_bwd_length_impl_1(&_S1377, _S1376);
    float2  _S1378 = _S1377.differential_0 + _S1366;
    float4  _S1379 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1380;
    (&_S1380)->primal_0 = _S986.primal_0;
    (&_S1380)->differential_0 = _S1379;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1381;
    (&_S1381)->primal_0 = _S987.primal_0;
    (&_S1381)->differential_0 = _S989;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1382;
    (&_S1382)->primal_0 = _S988.primal_0;
    (&_S1382)->differential_0 = _S989;
    CameraDistortion_0 _S1383 = _S1364;
    s_bwd_prop_CameraDistortion_x24init_0(&_S1380, &_S1381, &_S1382, &_S1383);
    dpthin_prism_coeffs_3->primal_0 = (*dpthin_prism_coeffs_3).primal_0;
    dpthin_prism_coeffs_3->differential_0 = _S1382.differential_0;
    dptangential_coeffs_3->primal_0 = (*dptangential_coeffs_3).primal_0;
    dptangential_coeffs_3->differential_0 = _S1381.differential_0;
    dpradial_coeffs_3->primal_0 = (*dpradial_coeffs_3).primal_0;
    dpradial_coeffs_3->differential_0 = _S1380.differential_0;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = dpmean2d_1.y;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = dpmean2d_1.x;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S1359;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S1358;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1141.differential_0;
    float3  _S1384 = _S1357 + make_float3 (_S1378.x, _S1378.y, _S1375);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1384;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_15, float fy_15, float cx_15, float cy_15, float4  radial_coeffs_20, float2  tangential_coeffs_20, float2  thin_prism_coeffs_20, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1385 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1386 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1387 = { _S1386, _S1386 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1388 = { _S1386, _S1386, _S1387, _S1386, _S1386, _S1387 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1389;
    (&_S1389)->_S736 = _S1385;
    (&_S1389)->_S737 = _S1388;
    float3  mean_c_11 = s_primal_ctx_mul_1(R_15, mean_11) + t_14;
    float3  _S1390 = s_primal_ctx_exp_0(scale_13);
    float _S1391 = quat_14.y;
    float x2_14 = _S1391 * _S1391;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_25 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S1392 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1390.x, 0.0f, 0.0f, 0.0f, _S1390.y, 0.0f, 0.0f, 0.0f, _S1390.z);
    Matrix<float, 3, 3>  _S1393 = s_primal_ctx_mul_2(_S1392, S_1);
    Matrix<float, 3, 3>  _S1394 = transpose_0(_S1393);
    Matrix<float, 3, 3>  _S1395 = s_primal_ctx_mul_2(_S1393, _S1394);
    Matrix<float, 3, 3>  _S1396 = s_primal_ctx_mul_2(R_15, _S1395);
    Matrix<float, 3, 3>  _S1397 = transpose_0(R_15);
    Matrix<float, 3, 3>  _S1398 = s_primal_ctx_mul_2(_S1396, _S1397);
    Matrix<float, 2, 2>  _S1399 = _S1385;
    float2  _S1400 = make_float2 (0.0f);
    float2  _S1401 = _S1400;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1398, fx_15, fy_15, cx_15, cy_15, radial_coeffs_20, tangential_coeffs_20, thin_prism_coeffs_20, &_S1399, &_S1401, &(&_S1389)->_S737);
    (&_S1389)->_S736 = _S1399;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1402 = _S1389;
    float _S1403 = _S1389._S736.rows[int(0)].y * _S1389._S736.rows[int(1)].x;
    float det_orig_12 = _S1389._S736.rows[int(0)].x * _S1389._S736.rows[int(1)].y - _S1403;
    float _S1404 = _S1389._S736.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1405 = _S1389._S736;
    *&(((&_S1405)->rows + (int(0)))->x) = _S1404;
    float _S1406 = _S1389._S736.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1405)->rows + (int(1)))->y) = _S1406;
    Matrix<float, 2, 2>  _S1407 = _S1405;
    Matrix<float, 2, 2>  _S1408 = _S1405;
    float det_blur_7 = _S1404 * _S1406 - _S1403;
    float _S1409 = det_orig_12 / det_blur_7;
    float _S1410 = det_blur_7 * det_blur_7;
    float _S1411 = s_primal_ctx_max_0(0.0f, _S1409);
    float _S1412 = s_primal_ctx_sqrt_0(_S1411);
    float invdet_8 = 1.0f / det_blur_7;
    float _S1413 = - _S1389._S736.rows[int(0)].y;
    float _S1414 = - _S1389._S736.rows[int(1)].x;
    float _S1415 = - in_opacity_11;
    float _S1416 = 1.0f + s_primal_ctx_exp_1(_S1415);
    float _S1417 = 1.0f / _S1416;
    float _S1418 = _S1416 * _S1416;
    float _S1419;
    if(antialiased_11)
    {
        _S1419 = _S1417 * _S1412;
    }
    else
    {
        _S1419 = _S1417;
    }
    float _S1420 = _S1419 / 0.00392156885936856f;
    float _S1421 = 2.0f * s_primal_ctx_log_0(_S1420);
    float _S1422 = s_primal_ctx_sqrt_0(_S1421);
    float _S1423 = _S1407.rows[int(0)].x;
    float _S1424 = _S1408.rows[int(1)].y;
    float _S1425 = s_primal_ctx_dot_0(mean_c_11, mean_c_11) + 9.99999997475242708e-07f;
    float3  _S1426 = mean_11 - - s_primal_ctx_mul_1(_S1397, t_14);
    float _S1427 = _S1426.x;
    float _S1428 = _S1426.y;
    float _S1429 = _S1426.z;
    float _S1430 = _S1427 * _S1427 + _S1428 * _S1428 + _S1429 * _S1429;
    float _S1431 = s_primal_ctx_sqrt_0(_S1430);
    float x_35 = _S1427 / _S1431;
    float3  _S1432 = make_float3 (x_35);
    float _S1433 = _S1431 * _S1431;
    float y_14 = _S1428 / _S1431;
    float z_11 = _S1429 / _S1431;
    float3  _S1434 = make_float3 (z_11);
    float _S1435 = - y_14;
    float3  _S1436 = make_float3 (_S1435);
    float z2_26 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_35 * x_35 - y_14 * y_14;
    float _S1437 = 2.0f * x_35;
    float fS1_11 = _S1437 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_26 - 0.31539157032966614f;
    float3  _S1438 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_35;
    float3  _S1439 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S1440 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S1441 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S1442 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S1443 = 1.86588168144226074f * z2_26 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S1443;
    float3  _S1444 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_35;
    float3  _S1445 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S1446 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S1447 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S1448 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_35 * fC1_11 - y_14 * fS1_11);
    float3  _S1449 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_35 * fS1_11 + y_14 * fC1_11);
    float3  _S1450 = make_float3 (pSH9_1);
    float3  _S1451 = make_float3 (0.0f);
    float3  _S1452 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1453;
    (&_S1453)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1435) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_35) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S1453)->differential_0 = _S1452;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1454;
    (&_S1454)->primal_0 = _S1451;
    (&_S1454)->differential_0 = _S1452;
    s_bwd_prop_max_0(&_S1453, &_S1454, v_rgb_1);
    float3  _S1455 = _S1449 * _S1453.differential_0;
    float3  _S1456 = (*sh_coeffs_11)[int(15)] * _S1453.differential_0;
    float3  _S1457 = _S1447 * _S1453.differential_0;
    float3  _S1458 = (*sh_coeffs_11)[int(14)] * _S1453.differential_0;
    float3  _S1459 = _S1445 * _S1453.differential_0;
    float3  _S1460 = (*sh_coeffs_11)[int(13)] * _S1453.differential_0;
    float3  _S1461 = _S1444 * _S1453.differential_0;
    float3  _S1462 = (*sh_coeffs_11)[int(12)] * _S1453.differential_0;
    float3  _S1463 = _S1446 * _S1453.differential_0;
    float3  _S1464 = (*sh_coeffs_11)[int(11)] * _S1453.differential_0;
    float3  _S1465 = _S1448 * _S1453.differential_0;
    float3  _S1466 = (*sh_coeffs_11)[int(10)] * _S1453.differential_0;
    float3  _S1467 = _S1450 * _S1453.differential_0;
    float3  _S1468 = (*sh_coeffs_11)[int(9)] * _S1453.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1468.x + _S1468.y + _S1468.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1456.x + _S1456.y + _S1456.z);
    float _S1469 = _S1466.x + _S1466.y + _S1466.z;
    float _S1470 = _S1458.x + _S1458.y + _S1458.z;
    float _S1471 = _S1464.x + _S1464.y + _S1464.z;
    float _S1472 = _S1460.x + _S1460.y + _S1460.z;
    float _S1473 = _S1462.x + _S1462.y + _S1462.z;
    float _S1474 = - s_diff_fC2_T_1;
    float3  _S1475 = _S1441 * _S1453.differential_0;
    float3  _S1476 = (*sh_coeffs_11)[int(8)] * _S1453.differential_0;
    float3  _S1477 = _S1439 * _S1453.differential_0;
    float3  _S1478 = (*sh_coeffs_11)[int(7)] * _S1453.differential_0;
    float3  _S1479 = _S1438 * _S1453.differential_0;
    float3  _S1480 = (*sh_coeffs_11)[int(6)] * _S1453.differential_0;
    float3  _S1481 = _S1440 * _S1453.differential_0;
    float3  _S1482 = (*sh_coeffs_11)[int(5)] * _S1453.differential_0;
    float3  _S1483 = _S1442 * _S1453.differential_0;
    float3  _S1484 = (*sh_coeffs_11)[int(4)] * _S1453.differential_0;
    float _S1485 = _S1482.x + _S1482.y + _S1482.z;
    float _S1486 = _S1478.x + _S1478.y + _S1478.z;
    float _S1487 = fTmp1B_11 * _S1469 + x_35 * s_diff_fS2_T_1 + y_14 * _S1474 + 0.54627424478530884f * (_S1484.x + _S1484.y + _S1484.z);
    float _S1488 = fTmp1B_11 * _S1470 + y_14 * s_diff_fS2_T_1 + x_35 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1476.x + _S1476.y + _S1476.z);
    float _S1489 = y_14 * - _S1488;
    float _S1490 = x_35 * _S1488;
    float _S1491 = z_11 * (1.86588168144226074f * (z_11 * _S1473) + -2.28522896766662598f * (y_14 * _S1471 + x_35 * _S1472) + 0.94617468118667603f * (_S1480.x + _S1480.y + _S1480.z));
    float3  _S1492 = make_float3 (0.48860251903533936f) * _S1453.differential_0;
    float3  _S1493 = - _S1492;
    float3  _S1494 = _S1432 * _S1493;
    float3  _S1495 = (*sh_coeffs_11)[int(3)] * _S1493;
    float3  _S1496 = _S1434 * _S1492;
    float3  _S1497 = (*sh_coeffs_11)[int(2)] * _S1492;
    float3  _S1498 = _S1436 * _S1492;
    float3  _S1499 = (*sh_coeffs_11)[int(1)] * _S1492;
    float _S1500 = (_S1443 * _S1473 + 1.44530570507049561f * (fS1_11 * _S1469 + fC1_11 * _S1470) + -1.09254848957061768f * (y_14 * _S1485 + x_35 * _S1486) + _S1491 + _S1491 + _S1497.x + _S1497.y + _S1497.z) / _S1433;
    float _S1501 = _S1431 * _S1500;
    float _S1502 = (fTmp0C_11 * _S1471 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S1474 + fTmp0B_11 * _S1485 + _S1437 * _S1487 + _S1489 + _S1489 + - (_S1499.x + _S1499.y + _S1499.z)) / _S1433;
    float _S1503 = _S1431 * _S1502;
    float _S1504 = (fTmp0C_11 * _S1472 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S1486 + 2.0f * (y_14 * _S1487) + _S1490 + _S1490 + _S1495.x + _S1495.y + _S1495.z) / _S1433;
    float _S1505 = _S1431 * _S1504;
    float _S1506 = _S1429 * - _S1500 + _S1428 * - _S1502 + _S1427 * - _S1504;
    DiffPair_float_0 _S1507;
    (&_S1507)->primal_0 = _S1430;
    (&_S1507)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1507, _S1506);
    float _S1508 = _S1429 * _S1507.differential_0;
    float _S1509 = _S1428 * _S1507.differential_0;
    float _S1510 = _S1427 * _S1507.differential_0;
    float3  _S1511 = make_float3 (0.282094806432724f) * _S1453.differential_0;
    float3  _S1512 = make_float3 (_S1505 + _S1510 + _S1510, _S1503 + _S1509 + _S1509, _S1501 + _S1508 + _S1508);
    float3  _S1513 = - - _S1512;
    Matrix<float, 3, 3>  _S1514 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1515;
    (&_S1515)->primal_0 = _S1397;
    (&_S1515)->differential_0 = _S1514;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1516;
    (&_S1516)->primal_0 = t_14;
    (&_S1516)->differential_0 = _S1452;
    s_bwd_prop_mul_1(&_S1515, &_S1516, _S1513);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1517 = _S1515;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1518 = _S1516;
    float2  _S1519 = _S1400;
    *&((&_S1519)->y) = v_conic_1.z;
    float2  _S1520 = _S1400;
    *&((&_S1520)->y) = v_conic_1.y;
    *&((&_S1520)->x) = v_conic_1.x;
    float _S1521 = 0.5f * v_depth_1;
    DiffPair_float_0 _S1522;
    (&_S1522)->primal_0 = _S1425;
    (&_S1522)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1522, _S1521);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1523;
    (&_S1523)->primal_0 = mean_c_11;
    (&_S1523)->differential_0 = _S1452;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1524;
    (&_S1524)->primal_0 = mean_c_11;
    (&_S1524)->differential_0 = _S1452;
    s_bwd_prop_dot_0(&_S1523, &_S1524, _S1522.differential_0);
    DiffPair_float_0 _S1525;
    (&_S1525)->primal_0 = _S1424;
    (&_S1525)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1525, 0.0f);
    DiffPair_float_0 _S1526;
    (&_S1526)->primal_0 = _S1423;
    (&_S1526)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1526, 0.0f);
    DiffPair_float_0 _S1527;
    (&_S1527)->primal_0 = 3.32999992370605469f;
    (&_S1527)->differential_0 = 0.0f;
    DiffPair_float_0 _S1528;
    (&_S1528)->primal_0 = _S1422;
    (&_S1528)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1527, &_S1528, 0.0f);
    DiffPair_float_0 _S1529;
    (&_S1529)->primal_0 = _S1421;
    (&_S1529)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1529, _S1528.differential_0);
    float _S1530 = 2.0f * _S1529.differential_0;
    DiffPair_float_0 _S1531;
    (&_S1531)->primal_0 = _S1420;
    (&_S1531)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1531, _S1530);
    float2  _S1532 = make_float2 (_S1526.differential_0, 0.0f);
    float _S1533 = v_opacity_1 + 254.9999847412109375f * _S1531.differential_0;
    Matrix<float, 2, 2>  _S1534 = _S1385;
    _S1534[int(1)] = _S1519;
    _S1534[int(0)] = _S1520;
    Matrix<float, 2, 2>  _S1535 = _S1534;
    FixedArray<float3 , 16>  _S1536;
    _S1536[int(0)] = _S1452;
    _S1536[int(1)] = _S1452;
    _S1536[int(2)] = _S1452;
    _S1536[int(3)] = _S1452;
    _S1536[int(4)] = _S1452;
    _S1536[int(5)] = _S1452;
    _S1536[int(6)] = _S1452;
    _S1536[int(7)] = _S1452;
    _S1536[int(8)] = _S1452;
    _S1536[int(9)] = _S1452;
    _S1536[int(10)] = _S1452;
    _S1536[int(11)] = _S1452;
    _S1536[int(12)] = _S1452;
    _S1536[int(13)] = _S1452;
    _S1536[int(14)] = _S1452;
    _S1536[int(15)] = _S1452;
    _S1536[int(7)] = _S1477;
    _S1536[int(0)] = _S1511;
    _S1536[int(1)] = _S1498;
    _S1536[int(2)] = _S1496;
    _S1536[int(3)] = _S1494;
    _S1536[int(4)] = _S1483;
    _S1536[int(5)] = _S1481;
    _S1536[int(6)] = _S1479;
    _S1536[int(15)] = _S1455;
    _S1536[int(8)] = _S1475;
    _S1536[int(9)] = _S1467;
    _S1536[int(10)] = _S1465;
    _S1536[int(11)] = _S1463;
    _S1536[int(12)] = _S1461;
    _S1536[int(13)] = _S1459;
    _S1536[int(14)] = _S1457;
    float3  _S1537 = _S1536[int(0)];
    float3  _S1538 = _S1536[int(1)];
    float3  _S1539 = _S1536[int(2)];
    float3  _S1540 = _S1536[int(3)];
    float3  _S1541 = _S1536[int(4)];
    float3  _S1542 = _S1536[int(5)];
    float3  _S1543 = _S1536[int(6)];
    float3  _S1544 = _S1536[int(7)];
    float3  _S1545 = _S1536[int(8)];
    float3  _S1546 = _S1536[int(9)];
    float3  _S1547 = _S1536[int(10)];
    float3  _S1548 = _S1536[int(11)];
    float3  _S1549 = _S1536[int(12)];
    float3  _S1550 = _S1536[int(13)];
    float3  _S1551 = _S1536[int(14)];
    float3  _S1552 = _S1536[int(15)];
    float3  _S1553 = _S1524.differential_0 + _S1523.differential_0;
    float2  _S1554 = make_float2 (0.0f, _S1525.differential_0);
    float _S1555;
    if(antialiased_11)
    {
        float _S1556 = _S1417 * _S1533;
        _S1419 = _S1412 * _S1533;
        _S1555 = _S1556;
    }
    else
    {
        _S1419 = _S1533;
        _S1555 = 0.0f;
    }
    float _S1557 = - (_S1419 / _S1418);
    DiffPair_float_0 _S1558;
    (&_S1558)->primal_0 = _S1415;
    (&_S1558)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1558, _S1557);
    float _S1559 = - _S1558.differential_0;
    float _S1560 = invdet_8 * _S1535.rows[int(1)].y;
    float _S1561 = - (invdet_8 * _S1535.rows[int(1)].x);
    float _S1562 = - (invdet_8 * _S1535.rows[int(0)].y);
    float _S1563 = invdet_8 * _S1535.rows[int(0)].x;
    float _S1564 = - ((_S1404 * _S1535.rows[int(1)].y + _S1414 * _S1535.rows[int(1)].x + _S1413 * _S1535.rows[int(0)].y + _S1406 * _S1535.rows[int(0)].x) / _S1410);
    DiffPair_float_0 _S1565;
    (&_S1565)->primal_0 = _S1411;
    (&_S1565)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1565, _S1555);
    DiffPair_float_0 _S1566;
    (&_S1566)->primal_0 = 0.0f;
    (&_S1566)->differential_0 = 0.0f;
    DiffPair_float_0 _S1567;
    (&_S1567)->primal_0 = _S1409;
    (&_S1567)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1566, &_S1567, _S1565.differential_0);
    float _S1568 = _S1567.differential_0 / _S1410;
    float s_diff_det_orig_T_1 = det_blur_7 * _S1568;
    float _S1569 = _S1564 + det_orig_12 * - _S1568;
    float _S1570 = - _S1569;
    float _S1571 = _S1404 * _S1569;
    float _S1572 = _S1406 * _S1569;
    Matrix<float, 2, 2>  _S1573 = _S1385;
    _S1573[int(1)] = _S1554;
    _S1573[int(0)] = _S1532;
    _S1405 = _S1573;
    *&(((&_S1405)->rows + (int(1)))->y) = 0.0f;
    float _S1574 = _S1563 + _S1571 + _S1573.rows[int(1)].y;
    *&(((&_S1405)->rows + (int(0)))->x) = 0.0f;
    float _S1575 = _S1560 + _S1572 + _S1573.rows[int(0)].x;
    float _S1576 = _S1570 + - s_diff_det_orig_T_1;
    float _S1577 = _S1561 + _S1402._S736.rows[int(0)].y * _S1576;
    float _S1578 = _S1562 + _S1402._S736.rows[int(1)].x * _S1576;
    float _S1579 = _S1402._S736.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S1580 = _S1574 + _S1402._S736.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S1581 = _S1400;
    *&((&_S1581)->x) = _S1577;
    *&((&_S1581)->y) = _S1580;
    float _S1582 = _S1575 + _S1579;
    float2  _S1583 = _S1400;
    *&((&_S1583)->y) = _S1578;
    *&((&_S1583)->x) = _S1582;
    Matrix<float, 2, 2>  _S1584 = _S1385;
    _S1584[int(1)] = _S1581;
    _S1584[int(0)] = _S1583;
    Matrix<float, 2, 2>  _S1585 = _S1405 + _S1584;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1586;
    (&_S1586)->primal_0 = mean_c_11;
    (&_S1586)->differential_0 = _S1452;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1587;
    (&_S1587)->primal_0 = _S1398;
    (&_S1587)->differential_0 = _S1514;
    DiffPair_float_0 _S1588;
    (&_S1588)->primal_0 = fx_15;
    (&_S1588)->differential_0 = 0.0f;
    DiffPair_float_0 _S1589;
    (&_S1589)->primal_0 = fy_15;
    (&_S1589)->differential_0 = 0.0f;
    DiffPair_float_0 _S1590;
    (&_S1590)->primal_0 = cx_15;
    (&_S1590)->differential_0 = 0.0f;
    DiffPair_float_0 _S1591;
    (&_S1591)->primal_0 = cy_15;
    (&_S1591)->differential_0 = 0.0f;
    float4  _S1592 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1593;
    (&_S1593)->primal_0 = radial_coeffs_20;
    (&_S1593)->differential_0 = _S1592;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1594;
    (&_S1594)->primal_0 = tangential_coeffs_20;
    (&_S1594)->differential_0 = _S1400;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1595;
    (&_S1595)->primal_0 = thin_prism_coeffs_20;
    (&_S1595)->differential_0 = _S1400;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1596 = _S1402._S737;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S1586, &_S1587, &_S1588, &_S1589, &_S1590, &_S1591, &_S1593, &_S1594, &_S1595, _S1585, v_mean2d_1, &_S1596);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1597;
    (&_S1597)->primal_0 = _S1396;
    (&_S1597)->differential_0 = _S1514;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1598;
    (&_S1598)->primal_0 = _S1397;
    (&_S1598)->differential_0 = _S1514;
    s_bwd_prop_mul_4(&_S1597, &_S1598, _S1587.differential_0);
    Matrix<float, 3, 3>  _S1599 = transpose_0(_S1598.differential_0 + _S1517.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1600;
    (&_S1600)->primal_0 = R_15;
    (&_S1600)->differential_0 = _S1514;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1601;
    (&_S1601)->primal_0 = _S1395;
    (&_S1601)->differential_0 = _S1514;
    s_bwd_prop_mul_4(&_S1600, &_S1601, _S1597.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1602;
    (&_S1602)->primal_0 = _S1393;
    (&_S1602)->differential_0 = _S1514;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1603;
    (&_S1603)->primal_0 = _S1394;
    (&_S1603)->differential_0 = _S1514;
    s_bwd_prop_mul_4(&_S1602, &_S1603, _S1601.differential_0);
    Matrix<float, 3, 3>  _S1604 = _S1602.differential_0 + transpose_0(_S1603.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1605;
    (&_S1605)->primal_0 = _S1392;
    (&_S1605)->differential_0 = _S1514;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1606;
    (&_S1606)->primal_0 = S_1;
    (&_S1606)->differential_0 = _S1514;
    s_bwd_prop_mul_4(&_S1605, &_S1606, _S1604);
    Matrix<float, 3, 3>  _S1607 = transpose_0(_S1605.differential_0);
    float _S1608 = 2.0f * - _S1607.rows[int(2)].z;
    float _S1609 = 2.0f * _S1607.rows[int(2)].y;
    float _S1610 = 2.0f * _S1607.rows[int(2)].x;
    float _S1611 = 2.0f * _S1607.rows[int(1)].z;
    float _S1612 = 2.0f * - _S1607.rows[int(1)].y;
    float _S1613 = 2.0f * _S1607.rows[int(1)].x;
    float _S1614 = 2.0f * _S1607.rows[int(0)].z;
    float _S1615 = 2.0f * _S1607.rows[int(0)].y;
    float _S1616 = 2.0f * - _S1607.rows[int(0)].x;
    float _S1617 = - _S1613 + _S1615;
    float _S1618 = _S1610 + - _S1614;
    float _S1619 = - _S1609 + _S1611;
    float _S1620 = _S1609 + _S1611;
    float _S1621 = _S1610 + _S1614;
    float _S1622 = _S1613 + _S1615;
    float _S1623 = quat_14.w * (_S1612 + _S1616);
    float _S1624 = quat_14.z * (_S1608 + _S1616);
    float _S1625 = quat_14.y * (_S1608 + _S1612);
    float _S1626 = quat_14.x * _S1617 + quat_14.z * _S1620 + quat_14.y * _S1621 + _S1623 + _S1623;
    float _S1627 = quat_14.x * _S1618 + quat_14.w * _S1620 + quat_14.y * _S1622 + _S1624 + _S1624;
    float _S1628 = quat_14.x * _S1619 + quat_14.w * _S1621 + quat_14.z * _S1622 + _S1625 + _S1625;
    float _S1629 = quat_14.w * _S1617 + quat_14.z * _S1618 + quat_14.y * _S1619;
    float3  _S1630 = _S1452;
    *&((&_S1630)->z) = _S1606.differential_0.rows[int(2)].z;
    *&((&_S1630)->y) = _S1606.differential_0.rows[int(1)].y;
    *&((&_S1630)->x) = _S1606.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1631;
    (&_S1631)->primal_0 = scale_13;
    (&_S1631)->differential_0 = _S1452;
    s_bwd_prop_exp_1(&_S1631, _S1630);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1632 = _S1631;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1633;
    (&_S1633)->primal_0 = mean_c_11;
    (&_S1633)->differential_0 = _S1452;
    s_bwd_length_impl_0(&_S1633, 0.0f);
    float3  _S1634 = _S1586.differential_0 + _S1633.differential_0 + _S1553;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1635;
    (&_S1635)->primal_0 = R_15;
    (&_S1635)->differential_0 = _S1514;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1636;
    (&_S1636)->primal_0 = mean_11;
    (&_S1636)->differential_0 = _S1452;
    s_bwd_prop_mul_1(&_S1635, &_S1636, _S1634);
    float3  _S1637 = _S1634 + _S1518.differential_0;
    Matrix<float, 3, 3>  _S1638 = _S1599 + _S1600.differential_0 + _S1635.differential_0;
    float4  _S1639 = _S1592;
    *&((&_S1639)->w) = _S1626;
    *&((&_S1639)->z) = _S1627;
    *&((&_S1639)->y) = _S1628;
    *&((&_S1639)->x) = _S1629;
    float4  _S1640 = _S1639;
    float3  _S1641 = _S1636.differential_0 + _S1512;
    *v_mean_1 = _S1641;
    *v_quat_1 = _S1640;
    *v_scale_1 = _S1632.differential_0;
    *v_in_opacity_1 = _S1559;
    (*v_sh_coeffs_1)[int(0)] = _S1537;
    (*v_sh_coeffs_1)[int(1)] = _S1538;
    (*v_sh_coeffs_1)[int(2)] = _S1539;
    (*v_sh_coeffs_1)[int(3)] = _S1540;
    (*v_sh_coeffs_1)[int(4)] = _S1541;
    (*v_sh_coeffs_1)[int(5)] = _S1542;
    (*v_sh_coeffs_1)[int(6)] = _S1543;
    (*v_sh_coeffs_1)[int(7)] = _S1544;
    (*v_sh_coeffs_1)[int(8)] = _S1545;
    (*v_sh_coeffs_1)[int(9)] = _S1546;
    (*v_sh_coeffs_1)[int(10)] = _S1547;
    (*v_sh_coeffs_1)[int(11)] = _S1548;
    (*v_sh_coeffs_1)[int(12)] = _S1549;
    (*v_sh_coeffs_1)[int(13)] = _S1550;
    (*v_sh_coeffs_1)[int(14)] = _S1551;
    (*v_sh_coeffs_1)[int(15)] = _S1552;
    *v_R_2 = _S1638;
    *v_t_2 = _S1637;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_16, float fy_16, float cx_16, float cy_16, float4  radial_coeffs_21, float2  tangential_coeffs_21, float2  thin_prism_coeffs_21, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_1(R_16, mean_12) + t_15;
    float3  _S1642 = s_primal_ctx_exp_0(scale_14);
    float _S1643 = quat_15.y;
    float x2_15 = _S1643 * _S1643;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_27 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S1644 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S1642.x, 0.0f, 0.0f, 0.0f, _S1642.y, 0.0f, 0.0f, 0.0f, _S1642.z);
    Matrix<float, 3, 3>  _S1645 = s_primal_ctx_mul_2(_S1644, S_2);
    Matrix<float, 3, 3>  _S1646 = transpose_0(_S1645);
    Matrix<float, 3, 3>  _S1647 = s_primal_ctx_mul_2(_S1645, _S1646);
    Matrix<float, 3, 3>  _S1648 = s_primal_ctx_mul_2(R_16, _S1647);
    Matrix<float, 3, 3>  _S1649 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S1650 = s_primal_ctx_mul_2(_S1648, _S1649);
    Matrix<float, 2, 3>  J_14 = makeMatrix<float, 2, 3> (fx_16, 0.0f, 0.0f, 0.0f, fy_16, 0.0f);
    Matrix<float, 2, 3>  _S1651 = s_primal_ctx_mul_3(J_14, _S1650);
    Matrix<float, 3, 2>  _S1652 = transpose_1(J_14);
    Matrix<float, 2, 2>  _S1653 = s_primal_ctx_mul_4(_S1651, _S1652);
    float _S1654 = _S1653.rows[int(0)].y * _S1653.rows[int(1)].x;
    float det_orig_13 = _S1653.rows[int(0)].x * _S1653.rows[int(1)].y - _S1654;
    float _S1655 = _S1653.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1656 = _S1653;
    *&(((&_S1656)->rows + (int(0)))->x) = _S1655;
    float _S1657 = _S1653.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1656)->rows + (int(1)))->y) = _S1657;
    Matrix<float, 2, 2>  _S1658 = _S1656;
    Matrix<float, 2, 2>  _S1659 = _S1656;
    float det_blur_8 = _S1655 * _S1657 - _S1654;
    float _S1660 = det_orig_13 / det_blur_8;
    float _S1661 = det_blur_8 * det_blur_8;
    float _S1662 = s_primal_ctx_max_0(0.0f, _S1660);
    float _S1663 = s_primal_ctx_sqrt_0(_S1662);
    float invdet_9 = 1.0f / det_blur_8;
    float _S1664 = - _S1653.rows[int(0)].y;
    float _S1665 = - _S1653.rows[int(1)].x;
    float _S1666 = - in_opacity_12;
    float _S1667 = 1.0f + s_primal_ctx_exp_1(_S1666);
    float _S1668 = 1.0f / _S1667;
    float _S1669 = _S1667 * _S1667;
    float _S1670;
    if(antialiased_12)
    {
        _S1670 = _S1668 * _S1663;
    }
    else
    {
        _S1670 = _S1668;
    }
    float _S1671 = _S1670 / 0.00392156885936856f;
    float _S1672 = 2.0f * s_primal_ctx_log_0(_S1671);
    float _S1673 = s_primal_ctx_sqrt_0(_S1672);
    float _S1674 = _S1658.rows[int(0)].x;
    float _S1675 = _S1659.rows[int(1)].y;
    float _S1676 = s_primal_ctx_dot_0(mean_c_12, mean_c_12) + 9.99999997475242708e-07f;
    float3  _S1677 = mean_12 - - s_primal_ctx_mul_1(_S1649, t_15);
    float _S1678 = _S1677.x;
    float _S1679 = _S1677.y;
    float _S1680 = _S1677.z;
    float _S1681 = _S1678 * _S1678 + _S1679 * _S1679 + _S1680 * _S1680;
    float _S1682 = s_primal_ctx_sqrt_0(_S1681);
    float x_36 = _S1678 / _S1682;
    float3  _S1683 = make_float3 (x_36);
    float _S1684 = _S1682 * _S1682;
    float y_15 = _S1679 / _S1682;
    float z_12 = _S1680 / _S1682;
    float3  _S1685 = make_float3 (z_12);
    float _S1686 = - y_15;
    float3  _S1687 = make_float3 (_S1686);
    float z2_28 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_36 * x_36 - y_15 * y_15;
    float _S1688 = 2.0f * x_36;
    float fS1_12 = _S1688 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_28 - 0.31539157032966614f;
    float3  _S1689 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_36;
    float3  _S1690 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S1691 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S1692 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S1693 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S1694 = 1.86588168144226074f * z2_28 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S1694;
    float3  _S1695 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_36;
    float3  _S1696 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S1697 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S1698 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S1699 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_36 * fC1_12 - y_15 * fS1_12);
    float3  _S1700 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_36 * fS1_12 + y_15 * fC1_12);
    float3  _S1701 = make_float3 (pSH9_2);
    float3  _S1702 = make_float3 (0.0f);
    float3  _S1703 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1704;
    (&_S1704)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1686) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_36) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S1704)->differential_0 = _S1703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1705;
    (&_S1705)->primal_0 = _S1702;
    (&_S1705)->differential_0 = _S1703;
    s_bwd_prop_max_0(&_S1704, &_S1705, v_rgb_2);
    float3  _S1706 = _S1700 * _S1704.differential_0;
    float3  _S1707 = (*sh_coeffs_12)[int(15)] * _S1704.differential_0;
    float3  _S1708 = _S1698 * _S1704.differential_0;
    float3  _S1709 = (*sh_coeffs_12)[int(14)] * _S1704.differential_0;
    float3  _S1710 = _S1696 * _S1704.differential_0;
    float3  _S1711 = (*sh_coeffs_12)[int(13)] * _S1704.differential_0;
    float3  _S1712 = _S1695 * _S1704.differential_0;
    float3  _S1713 = (*sh_coeffs_12)[int(12)] * _S1704.differential_0;
    float3  _S1714 = _S1697 * _S1704.differential_0;
    float3  _S1715 = (*sh_coeffs_12)[int(11)] * _S1704.differential_0;
    float3  _S1716 = _S1699 * _S1704.differential_0;
    float3  _S1717 = (*sh_coeffs_12)[int(10)] * _S1704.differential_0;
    float3  _S1718 = _S1701 * _S1704.differential_0;
    float3  _S1719 = (*sh_coeffs_12)[int(9)] * _S1704.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S1719.x + _S1719.y + _S1719.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S1707.x + _S1707.y + _S1707.z);
    float _S1720 = _S1717.x + _S1717.y + _S1717.z;
    float _S1721 = _S1709.x + _S1709.y + _S1709.z;
    float _S1722 = _S1715.x + _S1715.y + _S1715.z;
    float _S1723 = _S1711.x + _S1711.y + _S1711.z;
    float _S1724 = _S1713.x + _S1713.y + _S1713.z;
    float _S1725 = - s_diff_fC2_T_2;
    float3  _S1726 = _S1692 * _S1704.differential_0;
    float3  _S1727 = (*sh_coeffs_12)[int(8)] * _S1704.differential_0;
    float3  _S1728 = _S1690 * _S1704.differential_0;
    float3  _S1729 = (*sh_coeffs_12)[int(7)] * _S1704.differential_0;
    float3  _S1730 = _S1689 * _S1704.differential_0;
    float3  _S1731 = (*sh_coeffs_12)[int(6)] * _S1704.differential_0;
    float3  _S1732 = _S1691 * _S1704.differential_0;
    float3  _S1733 = (*sh_coeffs_12)[int(5)] * _S1704.differential_0;
    float3  _S1734 = _S1693 * _S1704.differential_0;
    float3  _S1735 = (*sh_coeffs_12)[int(4)] * _S1704.differential_0;
    float _S1736 = _S1733.x + _S1733.y + _S1733.z;
    float _S1737 = _S1729.x + _S1729.y + _S1729.z;
    float _S1738 = fTmp1B_12 * _S1720 + x_36 * s_diff_fS2_T_2 + y_15 * _S1725 + 0.54627424478530884f * (_S1735.x + _S1735.y + _S1735.z);
    float _S1739 = fTmp1B_12 * _S1721 + y_15 * s_diff_fS2_T_2 + x_36 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S1727.x + _S1727.y + _S1727.z);
    float _S1740 = y_15 * - _S1739;
    float _S1741 = x_36 * _S1739;
    float _S1742 = z_12 * (1.86588168144226074f * (z_12 * _S1724) + -2.28522896766662598f * (y_15 * _S1722 + x_36 * _S1723) + 0.94617468118667603f * (_S1731.x + _S1731.y + _S1731.z));
    float3  _S1743 = make_float3 (0.48860251903533936f) * _S1704.differential_0;
    float3  _S1744 = - _S1743;
    float3  _S1745 = _S1683 * _S1744;
    float3  _S1746 = (*sh_coeffs_12)[int(3)] * _S1744;
    float3  _S1747 = _S1685 * _S1743;
    float3  _S1748 = (*sh_coeffs_12)[int(2)] * _S1743;
    float3  _S1749 = _S1687 * _S1743;
    float3  _S1750 = (*sh_coeffs_12)[int(1)] * _S1743;
    float _S1751 = (_S1694 * _S1724 + 1.44530570507049561f * (fS1_12 * _S1720 + fC1_12 * _S1721) + -1.09254848957061768f * (y_15 * _S1736 + x_36 * _S1737) + _S1742 + _S1742 + _S1748.x + _S1748.y + _S1748.z) / _S1684;
    float _S1752 = _S1682 * _S1751;
    float _S1753 = (fTmp0C_12 * _S1722 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S1725 + fTmp0B_12 * _S1736 + _S1688 * _S1738 + _S1740 + _S1740 + - (_S1750.x + _S1750.y + _S1750.z)) / _S1684;
    float _S1754 = _S1682 * _S1753;
    float _S1755 = (fTmp0C_12 * _S1723 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S1737 + 2.0f * (y_15 * _S1738) + _S1741 + _S1741 + _S1746.x + _S1746.y + _S1746.z) / _S1684;
    float _S1756 = _S1682 * _S1755;
    float _S1757 = _S1680 * - _S1751 + _S1679 * - _S1753 + _S1678 * - _S1755;
    DiffPair_float_0 _S1758;
    (&_S1758)->primal_0 = _S1681;
    (&_S1758)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1758, _S1757);
    float _S1759 = _S1680 * _S1758.differential_0;
    float _S1760 = _S1679 * _S1758.differential_0;
    float _S1761 = _S1678 * _S1758.differential_0;
    float3  _S1762 = make_float3 (0.282094806432724f) * _S1704.differential_0;
    float3  _S1763 = make_float3 (_S1756 + _S1761 + _S1761, _S1754 + _S1760 + _S1760, _S1752 + _S1759 + _S1759);
    float3  _S1764 = - - _S1763;
    Matrix<float, 3, 3>  _S1765 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1766;
    (&_S1766)->primal_0 = _S1649;
    (&_S1766)->differential_0 = _S1765;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1767;
    (&_S1767)->primal_0 = t_15;
    (&_S1767)->differential_0 = _S1703;
    s_bwd_prop_mul_1(&_S1766, &_S1767, _S1764);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1768 = _S1766;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1769 = _S1767;
    float2  _S1770 = make_float2 (0.0f);
    float2  _S1771 = _S1770;
    *&((&_S1771)->y) = v_conic_2.z;
    float2  _S1772 = _S1770;
    *&((&_S1772)->y) = v_conic_2.y;
    *&((&_S1772)->x) = v_conic_2.x;
    float _S1773 = 0.5f * v_depth_2;
    DiffPair_float_0 _S1774;
    (&_S1774)->primal_0 = _S1676;
    (&_S1774)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1774, _S1773);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1775;
    (&_S1775)->primal_0 = mean_c_12;
    (&_S1775)->differential_0 = _S1703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1776;
    (&_S1776)->primal_0 = mean_c_12;
    (&_S1776)->differential_0 = _S1703;
    s_bwd_prop_dot_0(&_S1775, &_S1776, _S1774.differential_0);
    DiffPair_float_0 _S1777;
    (&_S1777)->primal_0 = _S1675;
    (&_S1777)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1777, 0.0f);
    DiffPair_float_0 _S1778;
    (&_S1778)->primal_0 = _S1674;
    (&_S1778)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1778, 0.0f);
    DiffPair_float_0 _S1779;
    (&_S1779)->primal_0 = 3.32999992370605469f;
    (&_S1779)->differential_0 = 0.0f;
    DiffPair_float_0 _S1780;
    (&_S1780)->primal_0 = _S1673;
    (&_S1780)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1779, &_S1780, 0.0f);
    DiffPair_float_0 _S1781;
    (&_S1781)->primal_0 = _S1672;
    (&_S1781)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1781, _S1780.differential_0);
    float _S1782 = 2.0f * _S1781.differential_0;
    DiffPair_float_0 _S1783;
    (&_S1783)->primal_0 = _S1671;
    (&_S1783)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1783, _S1782);
    float _S1784 = v_opacity_2 + 254.9999847412109375f * _S1783.differential_0;
    Matrix<float, 2, 2>  _S1785 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1786 = _S1785;
    _S1786[int(1)] = _S1771;
    _S1786[int(0)] = _S1772;
    Matrix<float, 2, 2>  _S1787 = _S1786;
    FixedArray<float3 , 16>  _S1788;
    _S1788[int(0)] = _S1703;
    _S1788[int(1)] = _S1703;
    _S1788[int(2)] = _S1703;
    _S1788[int(3)] = _S1703;
    _S1788[int(4)] = _S1703;
    _S1788[int(5)] = _S1703;
    _S1788[int(6)] = _S1703;
    _S1788[int(7)] = _S1703;
    _S1788[int(8)] = _S1703;
    _S1788[int(9)] = _S1703;
    _S1788[int(10)] = _S1703;
    _S1788[int(11)] = _S1703;
    _S1788[int(12)] = _S1703;
    _S1788[int(13)] = _S1703;
    _S1788[int(14)] = _S1703;
    _S1788[int(15)] = _S1703;
    _S1788[int(7)] = _S1728;
    _S1788[int(0)] = _S1762;
    _S1788[int(1)] = _S1749;
    _S1788[int(2)] = _S1747;
    _S1788[int(3)] = _S1745;
    _S1788[int(4)] = _S1734;
    _S1788[int(5)] = _S1732;
    _S1788[int(6)] = _S1730;
    _S1788[int(15)] = _S1706;
    _S1788[int(8)] = _S1726;
    _S1788[int(9)] = _S1718;
    _S1788[int(10)] = _S1716;
    _S1788[int(11)] = _S1714;
    _S1788[int(12)] = _S1712;
    _S1788[int(13)] = _S1710;
    _S1788[int(14)] = _S1708;
    float3  _S1789 = _S1788[int(0)];
    float3  _S1790 = _S1788[int(1)];
    float3  _S1791 = _S1788[int(2)];
    float3  _S1792 = _S1788[int(3)];
    float3  _S1793 = _S1788[int(4)];
    float3  _S1794 = _S1788[int(5)];
    float3  _S1795 = _S1788[int(6)];
    float3  _S1796 = _S1788[int(7)];
    float3  _S1797 = _S1788[int(8)];
    float3  _S1798 = _S1788[int(9)];
    float3  _S1799 = _S1788[int(10)];
    float3  _S1800 = _S1788[int(11)];
    float3  _S1801 = _S1788[int(12)];
    float3  _S1802 = _S1788[int(13)];
    float3  _S1803 = _S1788[int(14)];
    float3  _S1804 = _S1788[int(15)];
    float3  _S1805 = _S1776.differential_0 + _S1775.differential_0;
    float2  _S1806 = make_float2 (0.0f, _S1777.differential_0);
    float2  _S1807 = make_float2 (_S1778.differential_0, 0.0f);
    float _S1808;
    if(antialiased_12)
    {
        float _S1809 = _S1668 * _S1784;
        _S1670 = _S1663 * _S1784;
        _S1808 = _S1809;
    }
    else
    {
        _S1670 = _S1784;
        _S1808 = 0.0f;
    }
    float _S1810 = - (_S1670 / _S1669);
    DiffPair_float_0 _S1811;
    (&_S1811)->primal_0 = _S1666;
    (&_S1811)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1811, _S1810);
    float _S1812 = - _S1811.differential_0;
    float _S1813 = invdet_9 * _S1787.rows[int(1)].y;
    float _S1814 = - (invdet_9 * _S1787.rows[int(1)].x);
    float _S1815 = - (invdet_9 * _S1787.rows[int(0)].y);
    float _S1816 = invdet_9 * _S1787.rows[int(0)].x;
    float _S1817 = - ((_S1655 * _S1787.rows[int(1)].y + _S1665 * _S1787.rows[int(1)].x + _S1664 * _S1787.rows[int(0)].y + _S1657 * _S1787.rows[int(0)].x) / _S1661);
    DiffPair_float_0 _S1818;
    (&_S1818)->primal_0 = _S1662;
    (&_S1818)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1818, _S1808);
    DiffPair_float_0 _S1819;
    (&_S1819)->primal_0 = 0.0f;
    (&_S1819)->differential_0 = 0.0f;
    DiffPair_float_0 _S1820;
    (&_S1820)->primal_0 = _S1660;
    (&_S1820)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1819, &_S1820, _S1818.differential_0);
    float _S1821 = _S1820.differential_0 / _S1661;
    float s_diff_det_orig_T_2 = det_blur_8 * _S1821;
    float _S1822 = _S1817 + det_orig_13 * - _S1821;
    float _S1823 = - _S1822;
    float _S1824 = _S1655 * _S1822;
    float _S1825 = _S1657 * _S1822;
    Matrix<float, 2, 2>  _S1826 = _S1785;
    _S1826[int(1)] = _S1806;
    _S1826[int(0)] = _S1807;
    _S1656 = _S1826;
    *&(((&_S1656)->rows + (int(1)))->y) = 0.0f;
    float _S1827 = _S1816 + _S1824 + _S1826.rows[int(1)].y;
    *&(((&_S1656)->rows + (int(0)))->x) = 0.0f;
    float _S1828 = _S1813 + _S1825 + _S1826.rows[int(0)].x;
    float _S1829 = _S1823 + - s_diff_det_orig_T_2;
    float _S1830 = _S1814 + _S1653.rows[int(0)].y * _S1829;
    float _S1831 = _S1815 + _S1653.rows[int(1)].x * _S1829;
    float _S1832 = _S1653.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S1833 = _S1827 + _S1653.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S1834 = _S1770;
    *&((&_S1834)->x) = _S1830;
    *&((&_S1834)->y) = _S1833;
    float _S1835 = _S1828 + _S1832;
    float2  _S1836 = _S1770;
    *&((&_S1836)->y) = _S1831;
    *&((&_S1836)->x) = _S1835;
    float _S1837 = fy_16 * v_mean2d_2.y;
    float _S1838 = fx_16 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S1839 = _S1785;
    _S1839[int(1)] = _S1834;
    _S1839[int(0)] = _S1836;
    Matrix<float, 2, 2>  _S1840 = _S1656 + _S1839;
    Matrix<float, 2, 3>  _S1841 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1842;
    (&_S1842)->primal_0 = _S1651;
    (&_S1842)->differential_0 = _S1841;
    Matrix<float, 3, 2>  _S1843 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1844;
    (&_S1844)->primal_0 = _S1652;
    (&_S1844)->differential_0 = _S1843;
    s_bwd_prop_mul_2(&_S1842, &_S1844, _S1840);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1845;
    (&_S1845)->primal_0 = J_14;
    (&_S1845)->differential_0 = _S1841;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1846;
    (&_S1846)->primal_0 = _S1650;
    (&_S1846)->differential_0 = _S1765;
    s_bwd_prop_mul_3(&_S1845, &_S1846, _S1842.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1847;
    (&_S1847)->primal_0 = _S1648;
    (&_S1847)->differential_0 = _S1765;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1848;
    (&_S1848)->primal_0 = _S1649;
    (&_S1848)->differential_0 = _S1765;
    s_bwd_prop_mul_4(&_S1847, &_S1848, _S1846.differential_0);
    Matrix<float, 3, 3>  _S1849 = transpose_0(_S1848.differential_0 + _S1768.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1850;
    (&_S1850)->primal_0 = R_16;
    (&_S1850)->differential_0 = _S1765;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1851;
    (&_S1851)->primal_0 = _S1647;
    (&_S1851)->differential_0 = _S1765;
    s_bwd_prop_mul_4(&_S1850, &_S1851, _S1847.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1852;
    (&_S1852)->primal_0 = _S1645;
    (&_S1852)->differential_0 = _S1765;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1853;
    (&_S1853)->primal_0 = _S1646;
    (&_S1853)->differential_0 = _S1765;
    s_bwd_prop_mul_4(&_S1852, &_S1853, _S1851.differential_0);
    Matrix<float, 3, 3>  _S1854 = _S1852.differential_0 + transpose_0(_S1853.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1855;
    (&_S1855)->primal_0 = _S1644;
    (&_S1855)->differential_0 = _S1765;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1856;
    (&_S1856)->primal_0 = S_2;
    (&_S1856)->differential_0 = _S1765;
    s_bwd_prop_mul_4(&_S1855, &_S1856, _S1854);
    Matrix<float, 3, 3>  _S1857 = transpose_0(_S1855.differential_0);
    float _S1858 = 2.0f * - _S1857.rows[int(2)].z;
    float _S1859 = 2.0f * _S1857.rows[int(2)].y;
    float _S1860 = 2.0f * _S1857.rows[int(2)].x;
    float _S1861 = 2.0f * _S1857.rows[int(1)].z;
    float _S1862 = 2.0f * - _S1857.rows[int(1)].y;
    float _S1863 = 2.0f * _S1857.rows[int(1)].x;
    float _S1864 = 2.0f * _S1857.rows[int(0)].z;
    float _S1865 = 2.0f * _S1857.rows[int(0)].y;
    float _S1866 = 2.0f * - _S1857.rows[int(0)].x;
    float _S1867 = - _S1863 + _S1865;
    float _S1868 = _S1860 + - _S1864;
    float _S1869 = - _S1859 + _S1861;
    float _S1870 = _S1859 + _S1861;
    float _S1871 = _S1860 + _S1864;
    float _S1872 = _S1863 + _S1865;
    float _S1873 = quat_15.w * (_S1862 + _S1866);
    float _S1874 = quat_15.z * (_S1858 + _S1866);
    float _S1875 = quat_15.y * (_S1858 + _S1862);
    float _S1876 = quat_15.x * _S1867 + quat_15.z * _S1870 + quat_15.y * _S1871 + _S1873 + _S1873;
    float _S1877 = quat_15.x * _S1868 + quat_15.w * _S1870 + quat_15.y * _S1872 + _S1874 + _S1874;
    float _S1878 = quat_15.x * _S1869 + quat_15.w * _S1871 + quat_15.z * _S1872 + _S1875 + _S1875;
    float _S1879 = quat_15.w * _S1867 + quat_15.z * _S1868 + quat_15.y * _S1869;
    float3  _S1880 = _S1703;
    *&((&_S1880)->z) = _S1856.differential_0.rows[int(2)].z;
    *&((&_S1880)->y) = _S1856.differential_0.rows[int(1)].y;
    *&((&_S1880)->x) = _S1856.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1881;
    (&_S1881)->primal_0 = scale_14;
    (&_S1881)->differential_0 = _S1703;
    s_bwd_prop_exp_1(&_S1881, _S1880);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1882 = _S1881;
    float3  _S1883 = _S1703;
    *&((&_S1883)->y) = _S1837;
    *&((&_S1883)->x) = _S1838;
    float3  _S1884 = _S1805 + _S1883;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1885;
    (&_S1885)->primal_0 = R_16;
    (&_S1885)->differential_0 = _S1765;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1886;
    (&_S1886)->primal_0 = mean_12;
    (&_S1886)->differential_0 = _S1703;
    s_bwd_prop_mul_1(&_S1885, &_S1886, _S1884);
    float3  _S1887 = _S1884 + _S1769.differential_0;
    Matrix<float, 3, 3>  _S1888 = _S1849 + _S1850.differential_0 + _S1885.differential_0;
    float4  _S1889 = make_float4 (0.0f);
    *&((&_S1889)->w) = _S1876;
    *&((&_S1889)->z) = _S1877;
    *&((&_S1889)->y) = _S1878;
    *&((&_S1889)->x) = _S1879;
    float4  _S1890 = _S1889;
    float3  _S1891 = _S1886.differential_0 + _S1763;
    *v_mean_2 = _S1891;
    *v_quat_2 = _S1890;
    *v_scale_2 = _S1882.differential_0;
    *v_in_opacity_2 = _S1812;
    (*v_sh_coeffs_2)[int(0)] = _S1789;
    (*v_sh_coeffs_2)[int(1)] = _S1790;
    (*v_sh_coeffs_2)[int(2)] = _S1791;
    (*v_sh_coeffs_2)[int(3)] = _S1792;
    (*v_sh_coeffs_2)[int(4)] = _S1793;
    (*v_sh_coeffs_2)[int(5)] = _S1794;
    (*v_sh_coeffs_2)[int(6)] = _S1795;
    (*v_sh_coeffs_2)[int(7)] = _S1796;
    (*v_sh_coeffs_2)[int(8)] = _S1797;
    (*v_sh_coeffs_2)[int(9)] = _S1798;
    (*v_sh_coeffs_2)[int(10)] = _S1799;
    (*v_sh_coeffs_2)[int(11)] = _S1800;
    (*v_sh_coeffs_2)[int(12)] = _S1801;
    (*v_sh_coeffs_2)[int(13)] = _S1802;
    (*v_sh_coeffs_2)[int(14)] = _S1803;
    (*v_sh_coeffs_2)[int(15)] = _S1804;
    *v_R_3 = _S1888;
    *v_t_3 = _S1887;
    return;
}

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_17, float fy_17, float cx_17, float cy_17, float4  radial_coeffs_22, float2  tangential_coeffs_22, float2  thin_prism_coeffs_22, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_13) + t_16;
    float3  _S1892 = s_primal_ctx_exp_0(scale_15);
    float _S1893 = quat_16.y;
    float x2_16 = _S1893 * _S1893;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_29 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S1894 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    Matrix<float, 3, 3>  S_3 = makeMatrix<float, 3, 3> (_S1892.x, 0.0f, 0.0f, 0.0f, _S1892.y, 0.0f, 0.0f, 0.0f, _S1892.z);
    Matrix<float, 3, 3>  _S1895 = s_primal_ctx_mul_2(_S1894, S_3);
    Matrix<float, 3, 3>  _S1896 = transpose_0(_S1895);
    Matrix<float, 3, 3>  _S1897 = s_primal_ctx_mul_2(_S1895, _S1896);
    Matrix<float, 3, 3>  _S1898 = s_primal_ctx_mul_2(R_17, _S1897);
    Matrix<float, 3, 3>  _S1899 = transpose_0(R_17);
    Matrix<float, 3, 3>  _S1900 = s_primal_ctx_mul_2(_S1898, _S1899);
    float _S1901 = float(image_width_13);
    float _S1902 = float(image_height_13);
    float _S1903 = 0.30000001192092896f * (0.5f * _S1901 / fx_17);
    float lim_x_pos_1 = (_S1901 - cx_17) / fx_17 + _S1903;
    float _S1904 = 0.30000001192092896f * (0.5f * _S1902 / fy_17);
    float lim_y_pos_1 = (_S1902 - cy_17) / fy_17 + _S1904;
    float rz_7 = 1.0f / mean_c_13.z;
    float _S1905 = mean_c_13.z * mean_c_13.z;
    float rz2_7 = rz_7 * rz_7;
    float _S1906 = - (cx_17 / fx_17 + _S1903);
    float _S1907 = mean_c_13.x * rz_7;
    float _S1908 = s_primal_ctx_max_0(_S1906, _S1907);
    float _S1909 = s_primal_ctx_min_0(lim_x_pos_1, _S1908);
    float _S1910 = - (cy_17 / fy_17 + _S1904);
    float _S1911 = mean_c_13.y * rz_7;
    float _S1912 = s_primal_ctx_max_0(_S1910, _S1911);
    float _S1913 = s_primal_ctx_min_0(lim_y_pos_1, _S1912);
    float _S1914 = - fx_17;
    float _S1915 = _S1914 * (mean_c_13.z * _S1909);
    float _S1916 = - fy_17;
    float _S1917 = _S1916 * (mean_c_13.z * _S1913);
    Matrix<float, 2, 3>  J_15 = makeMatrix<float, 2, 3> (fx_17 * rz_7, 0.0f, _S1915 * rz2_7, 0.0f, fy_17 * rz_7, _S1917 * rz2_7);
    Matrix<float, 2, 3>  _S1918 = s_primal_ctx_mul_3(J_15, _S1900);
    Matrix<float, 3, 2>  _S1919 = transpose_1(J_15);
    Matrix<float, 2, 2>  _S1920 = s_primal_ctx_mul_4(_S1918, _S1919);
    float _S1921 = fx_17 * mean_c_13.x;
    float _S1922 = fy_17 * mean_c_13.y;
    float _S1923 = _S1920.rows[int(0)].y * _S1920.rows[int(1)].x;
    float det_orig_14 = _S1920.rows[int(0)].x * _S1920.rows[int(1)].y - _S1923;
    float _S1924 = _S1920.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1925 = _S1920;
    *&(((&_S1925)->rows + (int(0)))->x) = _S1924;
    float _S1926 = _S1920.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1925)->rows + (int(1)))->y) = _S1926;
    Matrix<float, 2, 2>  _S1927 = _S1925;
    Matrix<float, 2, 2>  _S1928 = _S1925;
    float det_blur_9 = _S1924 * _S1926 - _S1923;
    float _S1929 = det_orig_14 / det_blur_9;
    float _S1930 = det_blur_9 * det_blur_9;
    float _S1931 = s_primal_ctx_max_0(0.0f, _S1929);
    float _S1932 = s_primal_ctx_sqrt_0(_S1931);
    float _S1933 = - in_opacity_13;
    float _S1934 = 1.0f + s_primal_ctx_exp_1(_S1933);
    float _S1935 = 1.0f / _S1934;
    float _S1936 = _S1934 * _S1934;
    float _S1937;
    if(antialiased_13)
    {
        _S1937 = _S1935 * _S1932;
    }
    else
    {
        _S1937 = _S1935;
    }
    float _S1938 = _S1937 / 0.00392156885936856f;
    float _S1939 = 2.0f * s_primal_ctx_log_0(_S1938);
    float _S1940 = s_primal_ctx_sqrt_0(_S1939);
    float _S1941 = _S1927.rows[int(0)].x;
    float _S1942 = _S1928.rows[int(1)].y;
    float _S1943 = s_primal_ctx_dot_0(mean_c_13, mean_c_13) + 9.99999997475242708e-07f;
    float3  _S1944 = - scale_15;
    float3  _S1945 = mean_13 - - s_primal_ctx_mul_1(_S1899, t_16);
    float _S1946 = _S1945.x;
    float _S1947 = _S1945.y;
    float _S1948 = _S1945.z;
    float _S1949 = _S1946 * _S1946 + _S1947 * _S1947 + _S1948 * _S1948;
    float _S1950 = s_primal_ctx_sqrt_0(_S1949);
    float x_37 = _S1946 / _S1950;
    float3  _S1951 = make_float3 (x_37);
    float _S1952 = _S1950 * _S1950;
    float y_16 = _S1947 / _S1950;
    float z_13 = _S1948 / _S1950;
    float3  _S1953 = make_float3 (z_13);
    float _S1954 = - y_16;
    float3  _S1955 = make_float3 (_S1954);
    float z2_30 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_37 * x_37 - y_16 * y_16;
    float _S1956 = 2.0f * x_37;
    float fS1_13 = _S1956 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S1957 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_37;
    float3  _S1958 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S1959 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S1960 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S1961 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S1962 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S1962;
    float3  _S1963 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_37;
    float3  _S1964 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S1965 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S1966 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S1967 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_37 * fC1_13 - y_16 * fS1_13);
    float3  _S1968 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_37 * fS1_13 + y_16 * fC1_13);
    float3  _S1969 = make_float3 (pSH9_3);
    float3  _S1970 = make_float3 (0.0f);
    float3  _S1971 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1972;
    (&_S1972)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1954) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_37) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S1972)->differential_0 = _S1971;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1973;
    (&_S1973)->primal_0 = _S1970;
    (&_S1973)->differential_0 = _S1971;
    s_bwd_prop_max_0(&_S1972, &_S1973, v_rgb_3);
    float3  _S1974 = _S1968 * _S1972.differential_0;
    float3  _S1975 = (*sh_coeffs_13)[int(15)] * _S1972.differential_0;
    float3  _S1976 = _S1966 * _S1972.differential_0;
    float3  _S1977 = (*sh_coeffs_13)[int(14)] * _S1972.differential_0;
    float3  _S1978 = _S1964 * _S1972.differential_0;
    float3  _S1979 = (*sh_coeffs_13)[int(13)] * _S1972.differential_0;
    float3  _S1980 = _S1963 * _S1972.differential_0;
    float3  _S1981 = (*sh_coeffs_13)[int(12)] * _S1972.differential_0;
    float3  _S1982 = _S1965 * _S1972.differential_0;
    float3  _S1983 = (*sh_coeffs_13)[int(11)] * _S1972.differential_0;
    float3  _S1984 = _S1967 * _S1972.differential_0;
    float3  _S1985 = (*sh_coeffs_13)[int(10)] * _S1972.differential_0;
    float3  _S1986 = _S1969 * _S1972.differential_0;
    float3  _S1987 = (*sh_coeffs_13)[int(9)] * _S1972.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S1987.x + _S1987.y + _S1987.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S1975.x + _S1975.y + _S1975.z);
    float _S1988 = _S1985.x + _S1985.y + _S1985.z;
    float _S1989 = _S1977.x + _S1977.y + _S1977.z;
    float _S1990 = _S1983.x + _S1983.y + _S1983.z;
    float _S1991 = _S1979.x + _S1979.y + _S1979.z;
    float _S1992 = _S1981.x + _S1981.y + _S1981.z;
    float _S1993 = - s_diff_fC2_T_3;
    float3  _S1994 = _S1960 * _S1972.differential_0;
    float3  _S1995 = (*sh_coeffs_13)[int(8)] * _S1972.differential_0;
    float3  _S1996 = _S1958 * _S1972.differential_0;
    float3  _S1997 = (*sh_coeffs_13)[int(7)] * _S1972.differential_0;
    float3  _S1998 = _S1957 * _S1972.differential_0;
    float3  _S1999 = (*sh_coeffs_13)[int(6)] * _S1972.differential_0;
    float3  _S2000 = _S1959 * _S1972.differential_0;
    float3  _S2001 = (*sh_coeffs_13)[int(5)] * _S1972.differential_0;
    float3  _S2002 = _S1961 * _S1972.differential_0;
    float3  _S2003 = (*sh_coeffs_13)[int(4)] * _S1972.differential_0;
    float _S2004 = _S2001.x + _S2001.y + _S2001.z;
    float _S2005 = _S1997.x + _S1997.y + _S1997.z;
    float _S2006 = fTmp1B_13 * _S1988 + x_37 * s_diff_fS2_T_3 + y_16 * _S1993 + 0.54627424478530884f * (_S2003.x + _S2003.y + _S2003.z);
    float _S2007 = fTmp1B_13 * _S1989 + y_16 * s_diff_fS2_T_3 + x_37 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S1995.x + _S1995.y + _S1995.z);
    float _S2008 = y_16 * - _S2007;
    float _S2009 = x_37 * _S2007;
    float _S2010 = z_13 * (1.86588168144226074f * (z_13 * _S1992) + -2.28522896766662598f * (y_16 * _S1990 + x_37 * _S1991) + 0.94617468118667603f * (_S1999.x + _S1999.y + _S1999.z));
    float3  _S2011 = make_float3 (0.48860251903533936f) * _S1972.differential_0;
    float3  _S2012 = - _S2011;
    float3  _S2013 = _S1951 * _S2012;
    float3  _S2014 = (*sh_coeffs_13)[int(3)] * _S2012;
    float3  _S2015 = _S1953 * _S2011;
    float3  _S2016 = (*sh_coeffs_13)[int(2)] * _S2011;
    float3  _S2017 = _S1955 * _S2011;
    float3  _S2018 = (*sh_coeffs_13)[int(1)] * _S2011;
    float _S2019 = (_S1962 * _S1992 + 1.44530570507049561f * (fS1_13 * _S1988 + fC1_13 * _S1989) + -1.09254848957061768f * (y_16 * _S2004 + x_37 * _S2005) + _S2010 + _S2010 + _S2016.x + _S2016.y + _S2016.z) / _S1952;
    float _S2020 = _S1950 * _S2019;
    float _S2021 = (fTmp0C_13 * _S1990 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S1993 + fTmp0B_13 * _S2004 + _S1956 * _S2006 + _S2008 + _S2008 + - (_S2018.x + _S2018.y + _S2018.z)) / _S1952;
    float _S2022 = _S1950 * _S2021;
    float _S2023 = (fTmp0C_13 * _S1991 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2005 + 2.0f * (y_16 * _S2006) + _S2009 + _S2009 + _S2014.x + _S2014.y + _S2014.z) / _S1952;
    float _S2024 = _S1950 * _S2023;
    float _S2025 = _S1948 * - _S2019 + _S1947 * - _S2021 + _S1946 * - _S2023;
    DiffPair_float_0 _S2026;
    (&_S2026)->primal_0 = _S1949;
    (&_S2026)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2026, _S2025);
    float _S2027 = _S1948 * _S2026.differential_0;
    float _S2028 = _S1947 * _S2026.differential_0;
    float _S2029 = _S1946 * _S2026.differential_0;
    float3  _S2030 = make_float3 (0.282094806432724f) * _S1972.differential_0;
    float3  _S2031 = make_float3 (_S2024 + _S2029 + _S2029, _S2022 + _S2028 + _S2028, _S2020 + _S2027 + _S2027);
    float3  _S2032 = - - _S2031;
    Matrix<float, 3, 3>  _S2033 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2034;
    (&_S2034)->primal_0 = _S1899;
    (&_S2034)->differential_0 = _S2033;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2035;
    (&_S2035)->primal_0 = t_16;
    (&_S2035)->differential_0 = _S1971;
    s_bwd_prop_mul_1(&_S2034, &_S2035, _S2032);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2036 = _S2034;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2037 = _S2035;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2038;
    (&_S2038)->primal_0 = _S1944;
    (&_S2038)->differential_0 = _S1971;
    s_bwd_prop_exp_1(&_S2038, v_conic_3);
    float3  _S2039 = - _S2038.differential_0;
    float _S2040 = 0.5f * v_depth_3;
    DiffPair_float_0 _S2041;
    (&_S2041)->primal_0 = _S1943;
    (&_S2041)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2041, _S2040);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2042;
    (&_S2042)->primal_0 = mean_c_13;
    (&_S2042)->differential_0 = _S1971;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2043;
    (&_S2043)->primal_0 = mean_c_13;
    (&_S2043)->differential_0 = _S1971;
    s_bwd_prop_dot_0(&_S2042, &_S2043, _S2041.differential_0);
    DiffPair_float_0 _S2044;
    (&_S2044)->primal_0 = _S1942;
    (&_S2044)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2044, 0.0f);
    DiffPair_float_0 _S2045;
    (&_S2045)->primal_0 = _S1941;
    (&_S2045)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2045, 0.0f);
    DiffPair_float_0 _S2046;
    (&_S2046)->primal_0 = 3.32999992370605469f;
    (&_S2046)->differential_0 = 0.0f;
    DiffPair_float_0 _S2047;
    (&_S2047)->primal_0 = _S1940;
    (&_S2047)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2046, &_S2047, 0.0f);
    DiffPair_float_0 _S2048;
    (&_S2048)->primal_0 = _S1939;
    (&_S2048)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2048, _S2047.differential_0);
    float _S2049 = 2.0f * _S2048.differential_0;
    DiffPair_float_0 _S2050;
    (&_S2050)->primal_0 = _S1938;
    (&_S2050)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2050, _S2049);
    float _S2051 = v_opacity_3 + 254.9999847412109375f * _S2050.differential_0;
    FixedArray<float3 , 16>  _S2052;
    _S2052[int(0)] = _S1971;
    _S2052[int(1)] = _S1971;
    _S2052[int(2)] = _S1971;
    _S2052[int(3)] = _S1971;
    _S2052[int(4)] = _S1971;
    _S2052[int(5)] = _S1971;
    _S2052[int(6)] = _S1971;
    _S2052[int(7)] = _S1971;
    _S2052[int(8)] = _S1971;
    _S2052[int(9)] = _S1971;
    _S2052[int(10)] = _S1971;
    _S2052[int(11)] = _S1971;
    _S2052[int(12)] = _S1971;
    _S2052[int(13)] = _S1971;
    _S2052[int(14)] = _S1971;
    _S2052[int(15)] = _S1971;
    _S2052[int(7)] = _S1996;
    _S2052[int(0)] = _S2030;
    _S2052[int(1)] = _S2017;
    _S2052[int(2)] = _S2015;
    _S2052[int(3)] = _S2013;
    _S2052[int(4)] = _S2002;
    _S2052[int(5)] = _S2000;
    _S2052[int(6)] = _S1998;
    _S2052[int(15)] = _S1974;
    _S2052[int(8)] = _S1994;
    _S2052[int(9)] = _S1986;
    _S2052[int(10)] = _S1984;
    _S2052[int(11)] = _S1982;
    _S2052[int(12)] = _S1980;
    _S2052[int(13)] = _S1978;
    _S2052[int(14)] = _S1976;
    float3  _S2053 = _S2052[int(0)];
    float3  _S2054 = _S2052[int(1)];
    float3  _S2055 = _S2052[int(2)];
    float3  _S2056 = _S2052[int(3)];
    float3  _S2057 = _S2052[int(4)];
    float3  _S2058 = _S2052[int(5)];
    float3  _S2059 = _S2052[int(6)];
    float3  _S2060 = _S2052[int(7)];
    float3  _S2061 = _S2052[int(8)];
    float3  _S2062 = _S2052[int(9)];
    float3  _S2063 = _S2052[int(10)];
    float3  _S2064 = _S2052[int(11)];
    float3  _S2065 = _S2052[int(12)];
    float3  _S2066 = _S2052[int(13)];
    float3  _S2067 = _S2052[int(14)];
    float3  _S2068 = _S2052[int(15)];
    float3  _S2069 = _S2043.differential_0 + _S2042.differential_0;
    float2  _S2070 = make_float2 (0.0f, _S2044.differential_0);
    float2  _S2071 = make_float2 (_S2045.differential_0, 0.0f);
    float _S2072;
    if(antialiased_13)
    {
        float _S2073 = _S1935 * _S2051;
        _S1937 = _S1932 * _S2051;
        _S2072 = _S2073;
    }
    else
    {
        _S1937 = _S2051;
        _S2072 = 0.0f;
    }
    float _S2074 = - (_S1937 / _S1936);
    DiffPair_float_0 _S2075;
    (&_S2075)->primal_0 = _S1933;
    (&_S2075)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2075, _S2074);
    float _S2076 = - _S2075.differential_0;
    DiffPair_float_0 _S2077;
    (&_S2077)->primal_0 = _S1931;
    (&_S2077)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2077, _S2072);
    DiffPair_float_0 _S2078;
    (&_S2078)->primal_0 = 0.0f;
    (&_S2078)->differential_0 = 0.0f;
    DiffPair_float_0 _S2079;
    (&_S2079)->primal_0 = _S1929;
    (&_S2079)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2078, &_S2079, _S2077.differential_0);
    float _S2080 = _S2079.differential_0 / _S1930;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2080;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2080;
    float _S2081 = - s_diff_det_blur_T_0;
    float _S2082 = _S1924 * s_diff_det_blur_T_0;
    float _S2083 = _S1926 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2084 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2085 = _S2084;
    _S2085[int(1)] = _S2070;
    _S2085[int(0)] = _S2071;
    _S1925 = _S2085;
    *&(((&_S1925)->rows + (int(1)))->y) = 0.0f;
    float _S2086 = _S2082 + _S2085.rows[int(1)].y;
    *&(((&_S1925)->rows + (int(0)))->x) = 0.0f;
    float _S2087 = _S2083 + _S2085.rows[int(0)].x;
    float _S2088 = _S2081 + - s_diff_det_orig_T_3;
    float _S2089 = _S1920.rows[int(0)].y * _S2088;
    float _S2090 = _S1920.rows[int(1)].x * _S2088;
    float _S2091 = _S1920.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2092 = _S2086 + _S1920.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2093 = make_float2 (0.0f);
    float2  _S2094 = _S2093;
    *&((&_S2094)->x) = _S2089;
    *&((&_S2094)->y) = _S2092;
    float _S2095 = _S2087 + _S2091;
    float2  _S2096 = _S2093;
    *&((&_S2096)->y) = _S2090;
    *&((&_S2096)->x) = _S2095;
    float _S2097 = _S1922 * v_mean2d_3.y;
    float _S2098 = fy_17 * (rz_7 * v_mean2d_3.y);
    float _S2099 = _S1921 * v_mean2d_3.x;
    float _S2100 = fx_17 * (rz_7 * v_mean2d_3.x);
    Matrix<float, 2, 2>  _S2101 = _S2084;
    _S2101[int(1)] = _S2094;
    _S2101[int(0)] = _S2096;
    Matrix<float, 2, 2>  _S2102 = _S1925 + _S2101;
    Matrix<float, 2, 3>  _S2103 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2104;
    (&_S2104)->primal_0 = _S1918;
    (&_S2104)->differential_0 = _S2103;
    Matrix<float, 3, 2>  _S2105 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2106;
    (&_S2106)->primal_0 = _S1919;
    (&_S2106)->differential_0 = _S2105;
    s_bwd_prop_mul_2(&_S2104, &_S2106, _S2102);
    Matrix<float, 2, 3>  _S2107 = transpose_2(_S2106.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2108;
    (&_S2108)->primal_0 = J_15;
    (&_S2108)->differential_0 = _S2103;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2109;
    (&_S2109)->primal_0 = _S1900;
    (&_S2109)->differential_0 = _S2033;
    s_bwd_prop_mul_3(&_S2108, &_S2109, _S2104.differential_0);
    Matrix<float, 2, 3>  _S2110 = _S2107 + _S2108.differential_0;
    float _S2111 = _S1917 * _S2110.rows[int(1)].z;
    float s_diff_ty_T_1 = _S1916 * (rz2_7 * _S2110.rows[int(1)].z);
    float _S2112 = fy_17 * _S2110.rows[int(1)].y;
    float _S2113 = _S1915 * _S2110.rows[int(0)].z;
    float s_diff_tx_T_1 = _S1914 * (rz2_7 * _S2110.rows[int(0)].z);
    float _S2114 = fx_17 * _S2110.rows[int(0)].x;
    float _S2115 = mean_c_13.z * s_diff_ty_T_1;
    float _S2116 = _S1913 * s_diff_ty_T_1;
    DiffPair_float_0 _S2117;
    (&_S2117)->primal_0 = lim_y_pos_1;
    (&_S2117)->differential_0 = 0.0f;
    DiffPair_float_0 _S2118;
    (&_S2118)->primal_0 = _S1912;
    (&_S2118)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2117, &_S2118, _S2115);
    DiffPair_float_0 _S2119;
    (&_S2119)->primal_0 = _S1910;
    (&_S2119)->differential_0 = 0.0f;
    DiffPair_float_0 _S2120;
    (&_S2120)->primal_0 = _S1911;
    (&_S2120)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2119, &_S2120, _S2118.differential_0);
    float _S2121 = mean_c_13.y * _S2120.differential_0;
    float _S2122 = rz_7 * _S2120.differential_0;
    float _S2123 = mean_c_13.z * s_diff_tx_T_1;
    float _S2124 = _S1909 * s_diff_tx_T_1;
    DiffPair_float_0 _S2125;
    (&_S2125)->primal_0 = lim_x_pos_1;
    (&_S2125)->differential_0 = 0.0f;
    DiffPair_float_0 _S2126;
    (&_S2126)->primal_0 = _S1908;
    (&_S2126)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2125, &_S2126, _S2123);
    DiffPair_float_0 _S2127;
    (&_S2127)->primal_0 = _S1906;
    (&_S2127)->differential_0 = 0.0f;
    DiffPair_float_0 _S2128;
    (&_S2128)->primal_0 = _S1907;
    (&_S2128)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2127, &_S2128, _S2126.differential_0);
    float _S2129 = rz_7 * (_S2111 + _S2113);
    float _S2130 = _S2116 + _S2124 + - ((_S2097 + _S2099 + _S2112 + _S2114 + _S2121 + mean_c_13.x * _S2128.differential_0 + _S2129 + _S2129) / _S1905);
    float _S2131 = _S2098 + _S2122;
    float _S2132 = _S2100 + rz_7 * _S2128.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2133;
    (&_S2133)->primal_0 = _S1898;
    (&_S2133)->differential_0 = _S2033;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2134;
    (&_S2134)->primal_0 = _S1899;
    (&_S2134)->differential_0 = _S2033;
    s_bwd_prop_mul_4(&_S2133, &_S2134, _S2109.differential_0);
    Matrix<float, 3, 3>  _S2135 = transpose_0(_S2134.differential_0 + _S2036.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2136;
    (&_S2136)->primal_0 = R_17;
    (&_S2136)->differential_0 = _S2033;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2137;
    (&_S2137)->primal_0 = _S1897;
    (&_S2137)->differential_0 = _S2033;
    s_bwd_prop_mul_4(&_S2136, &_S2137, _S2133.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2138;
    (&_S2138)->primal_0 = _S1895;
    (&_S2138)->differential_0 = _S2033;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2139;
    (&_S2139)->primal_0 = _S1896;
    (&_S2139)->differential_0 = _S2033;
    s_bwd_prop_mul_4(&_S2138, &_S2139, _S2137.differential_0);
    Matrix<float, 3, 3>  _S2140 = _S2138.differential_0 + transpose_0(_S2139.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2141;
    (&_S2141)->primal_0 = _S1894;
    (&_S2141)->differential_0 = _S2033;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2142;
    (&_S2142)->primal_0 = S_3;
    (&_S2142)->differential_0 = _S2033;
    s_bwd_prop_mul_4(&_S2141, &_S2142, _S2140);
    Matrix<float, 3, 3>  _S2143 = transpose_0(_S2141.differential_0);
    float _S2144 = 2.0f * - _S2143.rows[int(2)].z;
    float _S2145 = 2.0f * _S2143.rows[int(2)].y;
    float _S2146 = 2.0f * _S2143.rows[int(2)].x;
    float _S2147 = 2.0f * _S2143.rows[int(1)].z;
    float _S2148 = 2.0f * - _S2143.rows[int(1)].y;
    float _S2149 = 2.0f * _S2143.rows[int(1)].x;
    float _S2150 = 2.0f * _S2143.rows[int(0)].z;
    float _S2151 = 2.0f * _S2143.rows[int(0)].y;
    float _S2152 = 2.0f * - _S2143.rows[int(0)].x;
    float _S2153 = - _S2149 + _S2151;
    float _S2154 = _S2146 + - _S2150;
    float _S2155 = - _S2145 + _S2147;
    float _S2156 = _S2145 + _S2147;
    float _S2157 = _S2146 + _S2150;
    float _S2158 = _S2149 + _S2151;
    float _S2159 = quat_16.w * (_S2148 + _S2152);
    float _S2160 = quat_16.z * (_S2144 + _S2152);
    float _S2161 = quat_16.y * (_S2144 + _S2148);
    float _S2162 = quat_16.x * _S2153 + quat_16.z * _S2156 + quat_16.y * _S2157 + _S2159 + _S2159;
    float _S2163 = quat_16.x * _S2154 + quat_16.w * _S2156 + quat_16.y * _S2158 + _S2160 + _S2160;
    float _S2164 = quat_16.x * _S2155 + quat_16.w * _S2157 + quat_16.z * _S2158 + _S2161 + _S2161;
    float _S2165 = quat_16.w * _S2153 + quat_16.z * _S2154 + quat_16.y * _S2155;
    float3  _S2166 = _S1971;
    *&((&_S2166)->z) = _S2142.differential_0.rows[int(2)].z;
    *&((&_S2166)->y) = _S2142.differential_0.rows[int(1)].y;
    *&((&_S2166)->x) = _S2142.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2167;
    (&_S2167)->primal_0 = scale_15;
    (&_S2167)->differential_0 = _S1971;
    s_bwd_prop_exp_1(&_S2167, _S2166);
    float3  _S2168 = _S1971;
    *&((&_S2168)->z) = _S2130;
    *&((&_S2168)->y) = _S2131;
    *&((&_S2168)->x) = _S2132;
    float3  _S2169 = _S2069 + _S2168;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2170;
    (&_S2170)->primal_0 = R_17;
    (&_S2170)->differential_0 = _S2033;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2171;
    (&_S2171)->primal_0 = mean_13;
    (&_S2171)->differential_0 = _S1971;
    s_bwd_prop_mul_1(&_S2170, &_S2171, _S2169);
    float3  _S2172 = _S2169 + _S2037.differential_0;
    Matrix<float, 3, 3>  _S2173 = _S2135 + _S2136.differential_0 + _S2170.differential_0;
    float3  _S2174 = _S2167.differential_0 + _S2039;
    float4  _S2175 = make_float4 (0.0f);
    *&((&_S2175)->w) = _S2162;
    *&((&_S2175)->z) = _S2163;
    *&((&_S2175)->y) = _S2164;
    *&((&_S2175)->x) = _S2165;
    float4  _S2176 = _S2175;
    float3  _S2177 = _S2171.differential_0 + _S2031;
    *v_mean_3 = _S2177;
    *v_quat_3 = _S2176;
    *v_scale_3 = _S2174;
    *v_in_opacity_3 = _S2076;
    (*v_sh_coeffs_3)[int(0)] = _S2053;
    (*v_sh_coeffs_3)[int(1)] = _S2054;
    (*v_sh_coeffs_3)[int(2)] = _S2055;
    (*v_sh_coeffs_3)[int(3)] = _S2056;
    (*v_sh_coeffs_3)[int(4)] = _S2057;
    (*v_sh_coeffs_3)[int(5)] = _S2058;
    (*v_sh_coeffs_3)[int(6)] = _S2059;
    (*v_sh_coeffs_3)[int(7)] = _S2060;
    (*v_sh_coeffs_3)[int(8)] = _S2061;
    (*v_sh_coeffs_3)[int(9)] = _S2062;
    (*v_sh_coeffs_3)[int(10)] = _S2063;
    (*v_sh_coeffs_3)[int(11)] = _S2064;
    (*v_sh_coeffs_3)[int(12)] = _S2065;
    (*v_sh_coeffs_3)[int(13)] = _S2066;
    (*v_sh_coeffs_3)[int(14)] = _S2067;
    (*v_sh_coeffs_3)[int(15)] = _S2068;
    *v_R_4 = _S2173;
    *v_t_4 = _S2172;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2178;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2179;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_18, float fy_18, float cx_18, float cy_18, float4  radial_coeffs_23, float2  tangential_coeffs_23, float2  thin_prism_coeffs_23, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2180 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S2181 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S2182 = { _S2181, _S2181 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2183 = { _S2181, _S2181, _S2182, _S2181, _S2181, _S2182 };
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2184;
    (&_S2184)->_S2178 = _S2180;
    (&_S2184)->_S2179 = _S2183;
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_14) + t_17;
    float3  _S2185 = s_primal_ctx_exp_0(scale_16);
    float _S2186 = quat_17.y;
    float x2_17 = _S2186 * _S2186;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_31 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2187 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    Matrix<float, 3, 3>  S_4 = makeMatrix<float, 3, 3> (_S2185.x, 0.0f, 0.0f, 0.0f, _S2185.y, 0.0f, 0.0f, 0.0f, _S2185.z);
    Matrix<float, 3, 3>  _S2188 = s_primal_ctx_mul_2(_S2187, S_4);
    Matrix<float, 3, 3>  _S2189 = transpose_0(_S2188);
    Matrix<float, 3, 3>  _S2190 = s_primal_ctx_mul_2(_S2188, _S2189);
    Matrix<float, 3, 3>  _S2191 = s_primal_ctx_mul_2(R_18, _S2190);
    Matrix<float, 3, 3>  _S2192 = transpose_0(R_18);
    Matrix<float, 3, 3>  _S2193 = s_primal_ctx_mul_2(_S2191, _S2192);
    Matrix<float, 2, 2>  _S2194 = _S2180;
    float2  _S2195 = make_float2 (0.0f);
    float2  _S2196 = _S2195;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_14, _S2193, fx_18, fy_18, cx_18, cy_18, radial_coeffs_23, tangential_coeffs_23, thin_prism_coeffs_23, &_S2194, &_S2196, &(&_S2184)->_S2179);
    (&_S2184)->_S2178 = _S2194;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2197 = _S2184;
    float _S2198 = _S2184._S2178.rows[int(0)].y * _S2184._S2178.rows[int(1)].x;
    float det_orig_15 = _S2184._S2178.rows[int(0)].x * _S2184._S2178.rows[int(1)].y - _S2198;
    float _S2199 = _S2184._S2178.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2200 = _S2184._S2178;
    *&(((&_S2200)->rows + (int(0)))->x) = _S2199;
    float _S2201 = _S2184._S2178.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2200)->rows + (int(1)))->y) = _S2201;
    Matrix<float, 2, 2>  _S2202 = _S2200;
    Matrix<float, 2, 2>  _S2203 = _S2200;
    float det_blur_10 = _S2199 * _S2201 - _S2198;
    float _S2204 = det_orig_15 / det_blur_10;
    float _S2205 = det_blur_10 * det_blur_10;
    float _S2206 = s_primal_ctx_max_0(0.0f, _S2204);
    float _S2207 = s_primal_ctx_sqrt_0(_S2206);
    float _S2208 = - in_opacity_14;
    float _S2209 = 1.0f + s_primal_ctx_exp_1(_S2208);
    float _S2210 = 1.0f / _S2209;
    float _S2211 = _S2209 * _S2209;
    float _S2212;
    if(antialiased_14)
    {
        _S2212 = _S2210 * _S2207;
    }
    else
    {
        _S2212 = _S2210;
    }
    float _S2213 = _S2212 / 0.00392156885936856f;
    float _S2214 = 2.0f * s_primal_ctx_log_0(_S2213);
    float _S2215 = s_primal_ctx_sqrt_0(_S2214);
    float _S2216 = _S2202.rows[int(0)].x;
    float _S2217 = _S2203.rows[int(1)].y;
    float _S2218 = s_primal_ctx_dot_0(mean_c_14, mean_c_14) + 9.99999997475242708e-07f;
    float3  _S2219 = - scale_16;
    float3  _S2220 = mean_14 - - s_primal_ctx_mul_1(_S2192, t_17);
    float _S2221 = _S2220.x;
    float _S2222 = _S2220.y;
    float _S2223 = _S2220.z;
    float _S2224 = _S2221 * _S2221 + _S2222 * _S2222 + _S2223 * _S2223;
    float _S2225 = s_primal_ctx_sqrt_0(_S2224);
    float x_38 = _S2221 / _S2225;
    float3  _S2226 = make_float3 (x_38);
    float _S2227 = _S2225 * _S2225;
    float y_17 = _S2222 / _S2225;
    float z_14 = _S2223 / _S2225;
    float3  _S2228 = make_float3 (z_14);
    float _S2229 = - y_17;
    float3  _S2230 = make_float3 (_S2229);
    float z2_32 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_38 * x_38 - y_17 * y_17;
    float _S2231 = 2.0f * x_38;
    float fS1_14 = _S2231 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2232 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_38;
    float3  _S2233 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2234 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2235 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2236 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2237 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2237;
    float3  _S2238 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_38;
    float3  _S2239 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2240 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2241 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2242 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_38 * fC1_14 - y_17 * fS1_14);
    float3  _S2243 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_38 * fS1_14 + y_17 * fC1_14);
    float3  _S2244 = make_float3 (pSH9_4);
    float3  _S2245 = make_float3 (0.0f);
    float3  _S2246 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2247;
    (&_S2247)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2229) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_38) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2247)->differential_0 = _S2246;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2248;
    (&_S2248)->primal_0 = _S2245;
    (&_S2248)->differential_0 = _S2246;
    s_bwd_prop_max_0(&_S2247, &_S2248, v_rgb_4);
    float3  _S2249 = _S2243 * _S2247.differential_0;
    float3  _S2250 = (*sh_coeffs_14)[int(15)] * _S2247.differential_0;
    float3  _S2251 = _S2241 * _S2247.differential_0;
    float3  _S2252 = (*sh_coeffs_14)[int(14)] * _S2247.differential_0;
    float3  _S2253 = _S2239 * _S2247.differential_0;
    float3  _S2254 = (*sh_coeffs_14)[int(13)] * _S2247.differential_0;
    float3  _S2255 = _S2238 * _S2247.differential_0;
    float3  _S2256 = (*sh_coeffs_14)[int(12)] * _S2247.differential_0;
    float3  _S2257 = _S2240 * _S2247.differential_0;
    float3  _S2258 = (*sh_coeffs_14)[int(11)] * _S2247.differential_0;
    float3  _S2259 = _S2242 * _S2247.differential_0;
    float3  _S2260 = (*sh_coeffs_14)[int(10)] * _S2247.differential_0;
    float3  _S2261 = _S2244 * _S2247.differential_0;
    float3  _S2262 = (*sh_coeffs_14)[int(9)] * _S2247.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2262.x + _S2262.y + _S2262.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2250.x + _S2250.y + _S2250.z);
    float _S2263 = _S2260.x + _S2260.y + _S2260.z;
    float _S2264 = _S2252.x + _S2252.y + _S2252.z;
    float _S2265 = _S2258.x + _S2258.y + _S2258.z;
    float _S2266 = _S2254.x + _S2254.y + _S2254.z;
    float _S2267 = _S2256.x + _S2256.y + _S2256.z;
    float _S2268 = - s_diff_fC2_T_4;
    float3  _S2269 = _S2235 * _S2247.differential_0;
    float3  _S2270 = (*sh_coeffs_14)[int(8)] * _S2247.differential_0;
    float3  _S2271 = _S2233 * _S2247.differential_0;
    float3  _S2272 = (*sh_coeffs_14)[int(7)] * _S2247.differential_0;
    float3  _S2273 = _S2232 * _S2247.differential_0;
    float3  _S2274 = (*sh_coeffs_14)[int(6)] * _S2247.differential_0;
    float3  _S2275 = _S2234 * _S2247.differential_0;
    float3  _S2276 = (*sh_coeffs_14)[int(5)] * _S2247.differential_0;
    float3  _S2277 = _S2236 * _S2247.differential_0;
    float3  _S2278 = (*sh_coeffs_14)[int(4)] * _S2247.differential_0;
    float _S2279 = _S2276.x + _S2276.y + _S2276.z;
    float _S2280 = _S2272.x + _S2272.y + _S2272.z;
    float _S2281 = fTmp1B_14 * _S2263 + x_38 * s_diff_fS2_T_4 + y_17 * _S2268 + 0.54627424478530884f * (_S2278.x + _S2278.y + _S2278.z);
    float _S2282 = fTmp1B_14 * _S2264 + y_17 * s_diff_fS2_T_4 + x_38 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2270.x + _S2270.y + _S2270.z);
    float _S2283 = y_17 * - _S2282;
    float _S2284 = x_38 * _S2282;
    float _S2285 = z_14 * (1.86588168144226074f * (z_14 * _S2267) + -2.28522896766662598f * (y_17 * _S2265 + x_38 * _S2266) + 0.94617468118667603f * (_S2274.x + _S2274.y + _S2274.z));
    float3  _S2286 = make_float3 (0.48860251903533936f) * _S2247.differential_0;
    float3  _S2287 = - _S2286;
    float3  _S2288 = _S2226 * _S2287;
    float3  _S2289 = (*sh_coeffs_14)[int(3)] * _S2287;
    float3  _S2290 = _S2228 * _S2286;
    float3  _S2291 = (*sh_coeffs_14)[int(2)] * _S2286;
    float3  _S2292 = _S2230 * _S2286;
    float3  _S2293 = (*sh_coeffs_14)[int(1)] * _S2286;
    float _S2294 = (_S2237 * _S2267 + 1.44530570507049561f * (fS1_14 * _S2263 + fC1_14 * _S2264) + -1.09254848957061768f * (y_17 * _S2279 + x_38 * _S2280) + _S2285 + _S2285 + _S2291.x + _S2291.y + _S2291.z) / _S2227;
    float _S2295 = _S2225 * _S2294;
    float _S2296 = (fTmp0C_14 * _S2265 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2268 + fTmp0B_14 * _S2279 + _S2231 * _S2281 + _S2283 + _S2283 + - (_S2293.x + _S2293.y + _S2293.z)) / _S2227;
    float _S2297 = _S2225 * _S2296;
    float _S2298 = (fTmp0C_14 * _S2266 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2280 + 2.0f * (y_17 * _S2281) + _S2284 + _S2284 + _S2289.x + _S2289.y + _S2289.z) / _S2227;
    float _S2299 = _S2225 * _S2298;
    float _S2300 = _S2223 * - _S2294 + _S2222 * - _S2296 + _S2221 * - _S2298;
    DiffPair_float_0 _S2301;
    (&_S2301)->primal_0 = _S2224;
    (&_S2301)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2301, _S2300);
    float _S2302 = _S2223 * _S2301.differential_0;
    float _S2303 = _S2222 * _S2301.differential_0;
    float _S2304 = _S2221 * _S2301.differential_0;
    float3  _S2305 = make_float3 (0.282094806432724f) * _S2247.differential_0;
    float3  _S2306 = make_float3 (_S2299 + _S2304 + _S2304, _S2297 + _S2303 + _S2303, _S2295 + _S2302 + _S2302);
    float3  _S2307 = - - _S2306;
    Matrix<float, 3, 3>  _S2308 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2309;
    (&_S2309)->primal_0 = _S2192;
    (&_S2309)->differential_0 = _S2308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2310;
    (&_S2310)->primal_0 = t_17;
    (&_S2310)->differential_0 = _S2246;
    s_bwd_prop_mul_1(&_S2309, &_S2310, _S2307);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2311 = _S2309;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2312 = _S2310;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2313;
    (&_S2313)->primal_0 = _S2219;
    (&_S2313)->differential_0 = _S2246;
    s_bwd_prop_exp_1(&_S2313, v_conic_4);
    float3  _S2314 = - _S2313.differential_0;
    float _S2315 = 0.5f * v_depth_4;
    DiffPair_float_0 _S2316;
    (&_S2316)->primal_0 = _S2218;
    (&_S2316)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2316, _S2315);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2317;
    (&_S2317)->primal_0 = mean_c_14;
    (&_S2317)->differential_0 = _S2246;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2318;
    (&_S2318)->primal_0 = mean_c_14;
    (&_S2318)->differential_0 = _S2246;
    s_bwd_prop_dot_0(&_S2317, &_S2318, _S2316.differential_0);
    DiffPair_float_0 _S2319;
    (&_S2319)->primal_0 = _S2217;
    (&_S2319)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2319, 0.0f);
    DiffPair_float_0 _S2320;
    (&_S2320)->primal_0 = _S2216;
    (&_S2320)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2320, 0.0f);
    DiffPair_float_0 _S2321;
    (&_S2321)->primal_0 = 3.32999992370605469f;
    (&_S2321)->differential_0 = 0.0f;
    DiffPair_float_0 _S2322;
    (&_S2322)->primal_0 = _S2215;
    (&_S2322)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2321, &_S2322, 0.0f);
    DiffPair_float_0 _S2323;
    (&_S2323)->primal_0 = _S2214;
    (&_S2323)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2323, _S2322.differential_0);
    float _S2324 = 2.0f * _S2323.differential_0;
    DiffPair_float_0 _S2325;
    (&_S2325)->primal_0 = _S2213;
    (&_S2325)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2325, _S2324);
    float2  _S2326 = make_float2 (_S2320.differential_0, 0.0f);
    float _S2327 = v_opacity_4 + 254.9999847412109375f * _S2325.differential_0;
    FixedArray<float3 , 16>  _S2328;
    _S2328[int(0)] = _S2246;
    _S2328[int(1)] = _S2246;
    _S2328[int(2)] = _S2246;
    _S2328[int(3)] = _S2246;
    _S2328[int(4)] = _S2246;
    _S2328[int(5)] = _S2246;
    _S2328[int(6)] = _S2246;
    _S2328[int(7)] = _S2246;
    _S2328[int(8)] = _S2246;
    _S2328[int(9)] = _S2246;
    _S2328[int(10)] = _S2246;
    _S2328[int(11)] = _S2246;
    _S2328[int(12)] = _S2246;
    _S2328[int(13)] = _S2246;
    _S2328[int(14)] = _S2246;
    _S2328[int(15)] = _S2246;
    _S2328[int(7)] = _S2271;
    _S2328[int(0)] = _S2305;
    _S2328[int(1)] = _S2292;
    _S2328[int(2)] = _S2290;
    _S2328[int(3)] = _S2288;
    _S2328[int(4)] = _S2277;
    _S2328[int(5)] = _S2275;
    _S2328[int(6)] = _S2273;
    _S2328[int(15)] = _S2249;
    _S2328[int(8)] = _S2269;
    _S2328[int(9)] = _S2261;
    _S2328[int(10)] = _S2259;
    _S2328[int(11)] = _S2257;
    _S2328[int(12)] = _S2255;
    _S2328[int(13)] = _S2253;
    _S2328[int(14)] = _S2251;
    float3  _S2329 = _S2328[int(0)];
    float3  _S2330 = _S2328[int(1)];
    float3  _S2331 = _S2328[int(2)];
    float3  _S2332 = _S2328[int(3)];
    float3  _S2333 = _S2328[int(4)];
    float3  _S2334 = _S2328[int(5)];
    float3  _S2335 = _S2328[int(6)];
    float3  _S2336 = _S2328[int(7)];
    float3  _S2337 = _S2328[int(8)];
    float3  _S2338 = _S2328[int(9)];
    float3  _S2339 = _S2328[int(10)];
    float3  _S2340 = _S2328[int(11)];
    float3  _S2341 = _S2328[int(12)];
    float3  _S2342 = _S2328[int(13)];
    float3  _S2343 = _S2328[int(14)];
    float3  _S2344 = _S2328[int(15)];
    float3  _S2345 = _S2318.differential_0 + _S2317.differential_0;
    float2  _S2346 = make_float2 (0.0f, _S2319.differential_0);
    float _S2347;
    if(antialiased_14)
    {
        float _S2348 = _S2210 * _S2327;
        _S2212 = _S2207 * _S2327;
        _S2347 = _S2348;
    }
    else
    {
        _S2212 = _S2327;
        _S2347 = 0.0f;
    }
    float _S2349 = - (_S2212 / _S2211);
    DiffPair_float_0 _S2350;
    (&_S2350)->primal_0 = _S2208;
    (&_S2350)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2350, _S2349);
    float _S2351 = - _S2350.differential_0;
    DiffPair_float_0 _S2352;
    (&_S2352)->primal_0 = _S2206;
    (&_S2352)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2352, _S2347);
    DiffPair_float_0 _S2353;
    (&_S2353)->primal_0 = 0.0f;
    (&_S2353)->differential_0 = 0.0f;
    DiffPair_float_0 _S2354;
    (&_S2354)->primal_0 = _S2204;
    (&_S2354)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2353, &_S2354, _S2352.differential_0);
    float _S2355 = _S2354.differential_0 / _S2205;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2355;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2355;
    float _S2356 = - s_diff_det_blur_T_1;
    float _S2357 = _S2199 * s_diff_det_blur_T_1;
    float _S2358 = _S2201 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2359 = _S2180;
    _S2359[int(1)] = _S2346;
    _S2359[int(0)] = _S2326;
    _S2200 = _S2359;
    *&(((&_S2200)->rows + (int(1)))->y) = 0.0f;
    float _S2360 = _S2357 + _S2359.rows[int(1)].y;
    *&(((&_S2200)->rows + (int(0)))->x) = 0.0f;
    float _S2361 = _S2358 + _S2359.rows[int(0)].x;
    float _S2362 = _S2356 + - s_diff_det_orig_T_4;
    float _S2363 = _S2197._S2178.rows[int(0)].y * _S2362;
    float _S2364 = _S2197._S2178.rows[int(1)].x * _S2362;
    float _S2365 = _S2197._S2178.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2366 = _S2360 + _S2197._S2178.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2367 = _S2195;
    *&((&_S2367)->x) = _S2363;
    *&((&_S2367)->y) = _S2366;
    float _S2368 = _S2361 + _S2365;
    float2  _S2369 = _S2195;
    *&((&_S2369)->y) = _S2364;
    *&((&_S2369)->x) = _S2368;
    Matrix<float, 2, 2>  _S2370 = _S2180;
    _S2370[int(1)] = _S2367;
    _S2370[int(0)] = _S2369;
    Matrix<float, 2, 2>  _S2371 = _S2200 + _S2370;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2372;
    (&_S2372)->primal_0 = mean_c_14;
    (&_S2372)->differential_0 = _S2246;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2373;
    (&_S2373)->primal_0 = _S2193;
    (&_S2373)->differential_0 = _S2308;
    DiffPair_float_0 _S2374;
    (&_S2374)->primal_0 = fx_18;
    (&_S2374)->differential_0 = 0.0f;
    DiffPair_float_0 _S2375;
    (&_S2375)->primal_0 = fy_18;
    (&_S2375)->differential_0 = 0.0f;
    DiffPair_float_0 _S2376;
    (&_S2376)->primal_0 = cx_18;
    (&_S2376)->differential_0 = 0.0f;
    DiffPair_float_0 _S2377;
    (&_S2377)->primal_0 = cy_18;
    (&_S2377)->differential_0 = 0.0f;
    float4  _S2378 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2379;
    (&_S2379)->primal_0 = radial_coeffs_23;
    (&_S2379)->differential_0 = _S2378;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2380;
    (&_S2380)->primal_0 = tangential_coeffs_23;
    (&_S2380)->differential_0 = _S2195;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2381;
    (&_S2381)->primal_0 = thin_prism_coeffs_23;
    (&_S2381)->differential_0 = _S2195;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2382 = _S2197._S2179;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2372, &_S2373, &_S2374, &_S2375, &_S2376, &_S2377, &_S2379, &_S2380, &_S2381, _S2371, v_mean2d_4, &_S2382);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2383;
    (&_S2383)->primal_0 = _S2191;
    (&_S2383)->differential_0 = _S2308;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2384;
    (&_S2384)->primal_0 = _S2192;
    (&_S2384)->differential_0 = _S2308;
    s_bwd_prop_mul_4(&_S2383, &_S2384, _S2373.differential_0);
    Matrix<float, 3, 3>  _S2385 = transpose_0(_S2384.differential_0 + _S2311.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2386;
    (&_S2386)->primal_0 = R_18;
    (&_S2386)->differential_0 = _S2308;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2387;
    (&_S2387)->primal_0 = _S2190;
    (&_S2387)->differential_0 = _S2308;
    s_bwd_prop_mul_4(&_S2386, &_S2387, _S2383.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2388;
    (&_S2388)->primal_0 = _S2188;
    (&_S2388)->differential_0 = _S2308;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2389;
    (&_S2389)->primal_0 = _S2189;
    (&_S2389)->differential_0 = _S2308;
    s_bwd_prop_mul_4(&_S2388, &_S2389, _S2387.differential_0);
    Matrix<float, 3, 3>  _S2390 = _S2388.differential_0 + transpose_0(_S2389.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2391;
    (&_S2391)->primal_0 = _S2187;
    (&_S2391)->differential_0 = _S2308;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2392;
    (&_S2392)->primal_0 = S_4;
    (&_S2392)->differential_0 = _S2308;
    s_bwd_prop_mul_4(&_S2391, &_S2392, _S2390);
    Matrix<float, 3, 3>  _S2393 = transpose_0(_S2391.differential_0);
    float _S2394 = 2.0f * - _S2393.rows[int(2)].z;
    float _S2395 = 2.0f * _S2393.rows[int(2)].y;
    float _S2396 = 2.0f * _S2393.rows[int(2)].x;
    float _S2397 = 2.0f * _S2393.rows[int(1)].z;
    float _S2398 = 2.0f * - _S2393.rows[int(1)].y;
    float _S2399 = 2.0f * _S2393.rows[int(1)].x;
    float _S2400 = 2.0f * _S2393.rows[int(0)].z;
    float _S2401 = 2.0f * _S2393.rows[int(0)].y;
    float _S2402 = 2.0f * - _S2393.rows[int(0)].x;
    float _S2403 = - _S2399 + _S2401;
    float _S2404 = _S2396 + - _S2400;
    float _S2405 = - _S2395 + _S2397;
    float _S2406 = _S2395 + _S2397;
    float _S2407 = _S2396 + _S2400;
    float _S2408 = _S2399 + _S2401;
    float _S2409 = quat_17.w * (_S2398 + _S2402);
    float _S2410 = quat_17.z * (_S2394 + _S2402);
    float _S2411 = quat_17.y * (_S2394 + _S2398);
    float _S2412 = quat_17.x * _S2403 + quat_17.z * _S2406 + quat_17.y * _S2407 + _S2409 + _S2409;
    float _S2413 = quat_17.x * _S2404 + quat_17.w * _S2406 + quat_17.y * _S2408 + _S2410 + _S2410;
    float _S2414 = quat_17.x * _S2405 + quat_17.w * _S2407 + quat_17.z * _S2408 + _S2411 + _S2411;
    float _S2415 = quat_17.w * _S2403 + quat_17.z * _S2404 + quat_17.y * _S2405;
    float3  _S2416 = _S2246;
    *&((&_S2416)->z) = _S2392.differential_0.rows[int(2)].z;
    *&((&_S2416)->y) = _S2392.differential_0.rows[int(1)].y;
    *&((&_S2416)->x) = _S2392.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2417;
    (&_S2417)->primal_0 = scale_16;
    (&_S2417)->differential_0 = _S2246;
    s_bwd_prop_exp_1(&_S2417, _S2416);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2418;
    (&_S2418)->primal_0 = mean_c_14;
    (&_S2418)->differential_0 = _S2246;
    s_bwd_length_impl_0(&_S2418, 0.0f);
    float3  _S2419 = _S2372.differential_0 + _S2418.differential_0 + _S2345;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2420;
    (&_S2420)->primal_0 = R_18;
    (&_S2420)->differential_0 = _S2308;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2421;
    (&_S2421)->primal_0 = mean_14;
    (&_S2421)->differential_0 = _S2246;
    s_bwd_prop_mul_1(&_S2420, &_S2421, _S2419);
    float3  _S2422 = _S2419 + _S2312.differential_0;
    Matrix<float, 3, 3>  _S2423 = _S2385 + _S2386.differential_0 + _S2420.differential_0;
    float3  _S2424 = _S2417.differential_0 + _S2314;
    float4  _S2425 = _S2378;
    *&((&_S2425)->w) = _S2412;
    *&((&_S2425)->z) = _S2413;
    *&((&_S2425)->y) = _S2414;
    *&((&_S2425)->x) = _S2415;
    float4  _S2426 = _S2425;
    float3  _S2427 = _S2421.differential_0 + _S2306;
    *v_mean_4 = _S2427;
    *v_quat_4 = _S2426;
    *v_scale_4 = _S2424;
    *v_in_opacity_4 = _S2351;
    (*v_sh_coeffs_4)[int(0)] = _S2329;
    (*v_sh_coeffs_4)[int(1)] = _S2330;
    (*v_sh_coeffs_4)[int(2)] = _S2331;
    (*v_sh_coeffs_4)[int(3)] = _S2332;
    (*v_sh_coeffs_4)[int(4)] = _S2333;
    (*v_sh_coeffs_4)[int(5)] = _S2334;
    (*v_sh_coeffs_4)[int(6)] = _S2335;
    (*v_sh_coeffs_4)[int(7)] = _S2336;
    (*v_sh_coeffs_4)[int(8)] = _S2337;
    (*v_sh_coeffs_4)[int(9)] = _S2338;
    (*v_sh_coeffs_4)[int(10)] = _S2339;
    (*v_sh_coeffs_4)[int(11)] = _S2340;
    (*v_sh_coeffs_4)[int(12)] = _S2341;
    (*v_sh_coeffs_4)[int(13)] = _S2342;
    (*v_sh_coeffs_4)[int(14)] = _S2343;
    (*v_sh_coeffs_4)[int(15)] = _S2344;
    *v_R_5 = _S2423;
    *v_t_5 = _S2422;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_18, float3  scale_17)
{
    float x_39 = quat_18.y;
    float x2_18 = x_39 * x_39;
    float y2_18 = quat_18.z * quat_18.z;
    float z2_33 = quat_18.w * quat_18.w;
    float xy_18 = quat_18.y * quat_18.z;
    float xz_18 = quat_18.y * quat_18.w;
    float yz_18 = quat_18.z * quat_18.w;
    float wx_18 = quat_18.x * quat_18.y;
    float wy_18 = quat_18.x * quat_18.z;
    float wz_18 = quat_18.x * quat_18.w;
    return mul_4(makeMatrix<float, 3, 3> (scale_17.x, 0.0f, 0.0f, 0.0f, scale_17.y, 0.0f, 0.0f, 0.0f, scale_17.z), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_18 + z2_33), 2.0f * (xy_18 + wz_18), 2.0f * (xz_18 - wy_18), 2.0f * (xy_18 - wz_18), 1.0f - 2.0f * (x2_18 + z2_33), 2.0f * (yz_18 + wx_18), 2.0f * (xz_18 + wy_18), 2.0f * (yz_18 - wx_18), 1.0f - 2.0f * (x2_18 + y2_18)))));
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_7)
{
    float _S2428 = (*dpquat_0).primal_0.y;
    float x2_19 = _S2428 * _S2428;
    float y2_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_34 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2429 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_34), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_34), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19))));
    Matrix<float, 3, 3>  _S2430 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2431;
    (&_S2431)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2431)->differential_0 = _S2430;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2432;
    (&_S2432)->primal_0 = _S2429;
    (&_S2432)->differential_0 = _S2430;
    s_bwd_prop_mul_4(&_S2431, &_S2432, _s_dOut_7);
    Matrix<float, 3, 3>  _S2433 = transpose_0(transpose_0(_S2432.differential_0));
    float _S2434 = 2.0f * - _S2433.rows[int(2)].z;
    float _S2435 = 2.0f * _S2433.rows[int(2)].y;
    float _S2436 = 2.0f * _S2433.rows[int(2)].x;
    float _S2437 = 2.0f * _S2433.rows[int(1)].z;
    float _S2438 = 2.0f * - _S2433.rows[int(1)].y;
    float _S2439 = 2.0f * _S2433.rows[int(1)].x;
    float _S2440 = 2.0f * _S2433.rows[int(0)].z;
    float _S2441 = 2.0f * _S2433.rows[int(0)].y;
    float _S2442 = 2.0f * - _S2433.rows[int(0)].x;
    float _S2443 = - _S2439 + _S2441;
    float _S2444 = _S2436 + - _S2440;
    float _S2445 = - _S2435 + _S2437;
    float _S2446 = _S2435 + _S2437;
    float _S2447 = _S2436 + _S2440;
    float _S2448 = _S2439 + _S2441;
    float _S2449 = (*dpquat_0).primal_0.w * (_S2438 + _S2442);
    float _S2450 = (*dpquat_0).primal_0.z * (_S2434 + _S2442);
    float _S2451 = (*dpquat_0).primal_0.y * (_S2434 + _S2438);
    float _S2452 = (*dpquat_0).primal_0.x * _S2443 + (*dpquat_0).primal_0.z * _S2446 + (*dpquat_0).primal_0.y * _S2447 + _S2449 + _S2449;
    float _S2453 = (*dpquat_0).primal_0.x * _S2444 + (*dpquat_0).primal_0.w * _S2446 + (*dpquat_0).primal_0.y * _S2448 + _S2450 + _S2450;
    float _S2454 = (*dpquat_0).primal_0.x * _S2445 + (*dpquat_0).primal_0.w * _S2447 + (*dpquat_0).primal_0.z * _S2448 + _S2451 + _S2451;
    float _S2455 = (*dpquat_0).primal_0.w * _S2443 + (*dpquat_0).primal_0.z * _S2444 + (*dpquat_0).primal_0.y * _S2445;
    float3  _S2456 = make_float3 (_S2431.differential_0.rows[int(0)].x, _S2431.differential_0.rows[int(1)].y, _S2431.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S2456;
    float4  _S2457 = make_float4 (0.0f);
    *&((&_S2457)->w) = _S2452;
    *&((&_S2457)->z) = _S2453;
    *&((&_S2457)->y) = _S2454;
    *&((&_S2457)->x) = _S2455;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S2457;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S2458, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2459, Matrix<float, 3, 3>  _S2460)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S2458, _S2459, _S2460);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_19, float3  scale_18, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S2461 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_19;
    (&dp_quat_0)->differential_0 = _S2461;
    float3  _S2462 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_18;
    (&dp_scale_0)->differential_0 = _S2462;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_15)
{
    float _S2463 = dOut_15.y;
    float _S2464 = dOut_15.z;
    float _S2465 = dOut_15.x;
    float _S2466 = (*a_0).primal_0.z * _S2463 + - (*a_0).primal_0.y * _S2464;
    float _S2467 = - (*a_0).primal_0.z * _S2465 + (*a_0).primal_0.x * _S2464;
    float _S2468 = (*a_0).primal_0.y * _S2465 + - (*a_0).primal_0.x * _S2463;
    float3  _S2469 = make_float3 (- (*b_0).primal_0.z * _S2463 + (*b_0).primal_0.y * _S2464, (*b_0).primal_0.z * _S2465 + - (*b_0).primal_0.x * _S2464, - (*b_0).primal_0.y * _S2465 + (*b_0).primal_0.x * _S2463);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S2469;
    float3  _S2470 = make_float3 (_S2466, _S2467, _S2468);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S2470;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S2471 = left_10.y;
    float _S2472 = right_10.z;
    float _S2473 = left_10.z;
    float _S2474 = right_10.y;
    float _S2475 = right_10.x;
    float _S2476 = left_10.x;
    return make_float3 (_S2471 * _S2472 - _S2473 * _S2474, _S2473 * _S2475 - _S2476 * _S2472, _S2476 * _S2474 - _S2471 * _S2475);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_15, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_1, float3  ray_d_1)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_1);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_1 - mean_15));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S2477, float3  _S2478)
{
    return cross_0(_S2477, _S2478);
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2479, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2480, float3  _S2481)
{
    _d_cross_0(_S2479, _S2480, _S2481);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_8)
{
    float3  _S2482 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S2483 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, _S2482);
    float3  _S2484 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, (*dpray_d_2).primal_0);
    float3  _S2485 = s_primal_ctx_cross_0(_S2484, _S2483);
    float _S2486 = -0.5f * s_primal_ctx_dot_0(_S2485, _S2485);
    float _S2487 = s_primal_ctx_dot_0(_S2484, _S2484);
    float _S2488 = _S2486 / _S2487;
    float _S2489 = _S2487 * _S2487;
    float _S2490 = (*dpopacity_0).primal_0 * _s_dOut_8;
    float _S2491 = s_primal_ctx_exp_1(_S2488) * _s_dOut_8;
    DiffPair_float_0 _S2492;
    (&_S2492)->primal_0 = _S2488;
    (&_S2492)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2492, _S2490);
    float _S2493 = _S2492.differential_0 / _S2489;
    float _S2494 = _S2486 * - _S2493;
    float _S2495 = _S2487 * _S2493;
    float3  _S2496 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2497;
    (&_S2497)->primal_0 = _S2484;
    (&_S2497)->differential_0 = _S2496;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2498;
    (&_S2498)->primal_0 = _S2484;
    (&_S2498)->differential_0 = _S2496;
    s_bwd_prop_dot_0(&_S2497, &_S2498, _S2494);
    float _S2499 = -0.5f * _S2495;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2500;
    (&_S2500)->primal_0 = _S2485;
    (&_S2500)->differential_0 = _S2496;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2501;
    (&_S2501)->primal_0 = _S2485;
    (&_S2501)->differential_0 = _S2496;
    s_bwd_prop_dot_0(&_S2500, &_S2501, _S2499);
    float3  _S2502 = _S2501.differential_0 + _S2500.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2503;
    (&_S2503)->primal_0 = _S2484;
    (&_S2503)->differential_0 = _S2496;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2504;
    (&_S2504)->primal_0 = _S2483;
    (&_S2504)->differential_0 = _S2496;
    s_bwd_prop_cross_0(&_S2503, &_S2504, _S2502);
    float3  _S2505 = _S2498.differential_0 + _S2497.differential_0 + _S2503.differential_0;
    Matrix<float, 3, 3>  _S2506 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2507;
    (&_S2507)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2507)->differential_0 = _S2506;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2508;
    (&_S2508)->primal_0 = (*dpray_d_2).primal_0;
    (&_S2508)->differential_0 = _S2496;
    s_bwd_prop_mul_1(&_S2507, &_S2508, _S2505);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2509;
    (&_S2509)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2509)->differential_0 = _S2506;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2510;
    (&_S2510)->primal_0 = _S2482;
    (&_S2510)->differential_0 = _S2496;
    s_bwd_prop_mul_1(&_S2509, &_S2510, _S2504.differential_0);
    float3  _S2511 = - _S2510.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S2508.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S2510.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S2491;
    Matrix<float, 3, 3>  _S2512 = _S2507.differential_0 + _S2509.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S2512;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S2511;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2513, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S2514, DiffPair_float_0 * _S2515, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2516, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2517, float _S2518)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S2513, _S2514, _S2515, _S2516, _S2517, _S2518);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S2519 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_16;
    (&dp_mean_0)->differential_0 = _S2519;
    Matrix<float, 3, 3>  _S2520 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S2520;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S2519;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S2519;
    s_bwd_evaluate_alpha_3dgs_0(&dp_mean_0, &dp_iscl_rot_0, &dp_opacity_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_mean_5 = dp_mean_0.differential_0;
    *v_iscl_rot_1 = dp_iscl_rot_0.differential_0;
    *v_opacity_5 = dp_opacity_0.differential_0;
    *v_ray_o_1 = dp_ray_o_0.differential_0;
    *v_ray_d_1 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_3dgs(float3  mean_17, Matrix<float, 3, 3>  iscl_rot_2, float opacity_12, float3  rgb_10, float3  ray_o_3, float3  ray_d_3, float3  * out_rgb_0, float * depth_10)
{
    *out_rgb_0 = rgb_10;
    float3  _S2521 = mean_17 - ray_o_3;
    *depth_10 = 0.5f * (F32_log((dot_0(_S2521, _S2521) + 9.99999997475242708e-07f)));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S2522 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    float _S2523 = 0.5f * dpdepth_0;
    DiffPair_float_0 _S2524;
    (&_S2524)->primal_0 = s_primal_ctx_dot_0(_S2522, _S2522) + 9.99999997475242708e-07f;
    (&_S2524)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2524, _S2523);
    float3  _S2525 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2526;
    (&_S2526)->primal_0 = _S2522;
    (&_S2526)->differential_0 = _S2525;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2527;
    (&_S2527)->primal_0 = _S2522;
    (&_S2527)->differential_0 = _S2525;
    s_bwd_prop_dot_0(&_S2526, &_S2527, _S2524.differential_0);
    float3  _S2528 = _S2527.differential_0 + _S2526.differential_0;
    float3  _S2529 = - _S2528;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S2525;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S2529;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S2530 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S2530;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S2528;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2531, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S2532, DiffPair_float_0 * _S2533, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2534, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2535, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2536, float3  _S2537, float _S2538)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S2531, _S2532, _S2533, _S2534, _S2535, _S2536, _S2537, _S2538);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S2539 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_18;
    (&dp_mean_1)->differential_0 = _S2539;
    Matrix<float, 3, 3>  _S2540 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S2540;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S2539;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S2539;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S2539;
    s_bwd_evaluate_color_3dgs_0(&dp_mean_1, &dp_iscl_rot_1, &dp_opacity_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_5);
    *v_mean_6 = dp_mean_1.differential_0;
    *v_iscl_rot_2 = dp_iscl_rot_1.differential_0;
    *v_opacity_6 = dp_opacity_1.differential_0;
    *v_rgb_5 = dp_rgb_0.differential_0;
    *v_ray_o_2 = dp_ray_o_1.differential_0;
    *v_ray_d_2 = dp_ray_d_1.differential_0;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_13, float dOut_16)
{
    float _S2541 = _slang_select(((*dpx_13).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_13).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_16;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S2541;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_14, float dOut_17)
{
    float _S2542 = (F32_exp2(((*dpx_14).primal_0))) * 50.693145751953125f * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S2542;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float3  dOut_18)
{
    float3  _S2543 = make_float3 (1.0f) / (*dpx_15).primal_0 * dOut_18;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S2543;
    return;
}

inline __device__ float3  log_0(float3  x_40)
{
    float3  result_16;
    int i_10 = int(0);
    for(;;)
    {
        if(i_10 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_16, i_10) = (F32_log((_slang_vector_get_element(x_40, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_16;
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_19, float4  quat_20, float3  scale_19, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_19, float3  t_18, float fx_19, float fy_19, float cx_19, float cy_19, float4  radial_coeffs_24, float2  tangential_coeffs_24, float2  thin_prism_coeffs_24, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_19, mean_19) + t_18;
        float _S2544 = mean_c_15.z;
        bool _S2545;
        if(_S2544 < near_plane_10)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = _S2544 > far_plane_10;
        }
        if(_S2545)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S2546 = scale_19.x;
        float sx_0 = (F32_exp((_S2546)));
        float _S2547 = scale_19.y;
        float sy_0 = (F32_exp((_S2547)));
        float sz_0 = scale_19.z - 0.5f * (_S2546 + _S2547);
        float x_41 = quat_20.y;
        float x2_20 = x_41 * x_41;
        float y2_20 = quat_20.z * quat_20.z;
        float z2_35 = quat_20.w * quat_20.w;
        float xy_20 = quat_20.y * quat_20.z;
        float xz_20 = quat_20.y * quat_20.w;
        float yz_20 = quat_20.z * quat_20.w;
        float wx_20 = quat_20.x * quat_20.y;
        float wy_20 = quat_20.x * quat_20.z;
        float wz_20 = quat_20.x * quat_20.w;
        Matrix<float, 3, 3>  _S2548 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20)));
        float3  vert0_c_0 = mul_0(R_19, mul_0(_S2548, make_float3 (sx_0, 0.0f, 0.0f)) + mean_19) + t_18;
        float3  vert1_c_0 = mul_0(R_19, mul_0(_S2548, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_19) + t_18;
        float3  vert2_c_0 = mul_0(R_19, mul_0(_S2548, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_19) + t_18;
        float _S2549 = vert0_c_0.z;
        float _S2550 = vert1_c_0.z;
        float _S2551 = vert2_c_0.z;
        if(_S2549 < near_plane_10)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = _S2549 > far_plane_10;
        }
        if(_S2545)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = _S2550 < near_plane_10;
        }
        if(_S2545)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = _S2550 > far_plane_10;
        }
        if(_S2545)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = _S2551 < near_plane_10;
        }
        if(_S2545)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = _S2551 > far_plane_10;
        }
        if(_S2545)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S2549);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S2550);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S2551);
        float2  _S2552 = make_float2 (fx_19, fy_19);
        float2  _S2553 = make_float2 (cx_19, cy_19);
        *uv0_0 = _S2552 * *uv0_0 + _S2553;
        *uv1_0 = _S2552 * *uv1_0 + _S2553;
        float2  _S2554 = _S2552 * *uv2_0 + _S2553;
        *uv2_0 = _S2554;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S2554 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S2554)));
        float _S2555 = _S2554.x;
        float xmax_5 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S2555))) + offset_0;
        float xmin_5 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S2555))) - offset_0;
        float _S2556 = _S2554.y;
        float ymax_5 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S2556))) + offset_0;
        float ymin_5 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S2556))) - offset_0;
        if(xmax_5 <= 0.0f)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = xmin_5 >= float(image_width_15);
        }
        if(_S2545)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = ymax_5 <= 0.0f;
        }
        if(_S2545)
        {
            _S2545 = true;
        }
        else
        {
            _S2545 = ymin_5 >= float(image_height_15);
        }
        if(_S2545)
        {
            _S2545 = true;
        }
        else
        {
            if(_S2544 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S2545 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S2545 = false;
                }
                if(_S2545)
                {
                    _S2545 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S2545 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S2545 = false;
                    }
                }
            }
            else
            {
                _S2545 = false;
            }
        }
        if(_S2545)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S2557 = mean_19 - - mul_0(transpose_0(R_19), t_18);
        float _S2558 = _S2557.x;
        float _S2559 = _S2557.y;
        float _S2560 = _S2557.z;
        float norm_10 = (F32_sqrt((_S2558 * _S2558 + _S2559 * _S2559 + _S2560 * _S2560)));
        float x_42 = _S2558 / norm_10;
        float y_18 = _S2559 / norm_10;
        float z_15 = _S2560 / norm_10;
        float z2_36 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_42 * x_42 - y_18 * y_18;
        float fS1_15 = 2.0f * x_42 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_36 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_42) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_36 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_42) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_42 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_36 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_42) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_42 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S2561 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S2561);
        float3  _S2562 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S2563 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S2562 + _S2563 + make_float3 (0.5f), _S2561);
        (*rgb_12)[int(2)] = max_0(_S2562 - _S2563 + make_float3 (0.5f), _S2561);
        float3  _S2564 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S2564 * make_float3 (float(- (F32_sign((dot_0(_S2564, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_20, float3  t_19, float fx_20, float fy_20, float cx_20, float cy_20, float4  radial_coeffs_25, float2  tangential_coeffs_25, float2  thin_prism_coeffs_25, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_20, mean_20) + t_19;
        float _S2565 = length_1(mean_c_16);
        bool _S2566;
        if(_S2565 < near_plane_11)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = _S2565 > far_plane_11;
        }
        if(_S2566)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S2567 = scale_20.x;
        float sx_1 = (F32_exp((_S2567)));
        float _S2568 = scale_20.y;
        float sy_1 = (F32_exp((_S2568)));
        float sz_1 = scale_20.z - 0.5f * (_S2567 + _S2568);
        float x_43 = quat_21.y;
        float x2_21 = x_43 * x_43;
        float y2_21 = quat_21.z * quat_21.z;
        float z2_37 = quat_21.w * quat_21.w;
        float xy_21 = quat_21.y * quat_21.z;
        float xz_21 = quat_21.y * quat_21.w;
        float yz_21 = quat_21.z * quat_21.w;
        float wx_21 = quat_21.x * quat_21.y;
        float wy_21 = quat_21.x * quat_21.z;
        float wz_21 = quat_21.x * quat_21.w;
        Matrix<float, 3, 3>  _S2569 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_37), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_37), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_1 = mul_0(R_20, mul_0(_S2569, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_19;
        float3  vert1_c_1 = mul_0(R_20, mul_0(_S2569, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_19;
        float3  vert2_c_1 = mul_0(R_20, mul_0(_S2569, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_19;
        float _S2570 = length_1(vert0_c_1);
        float _S2571 = length_1(vert1_c_1);
        float _S2572 = length_1(vert2_c_1);
        if(_S2570 < near_plane_11)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = _S2570 > far_plane_11;
        }
        if(_S2566)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = _S2571 < near_plane_11;
        }
        if(_S2566)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = _S2571 > far_plane_11;
        }
        if(_S2566)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = _S2572 < near_plane_11;
        }
        if(_S2566)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = _S2572 > far_plane_11;
        }
        if(_S2566)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_2 = CameraDistortion_x24init_0(radial_coeffs_25, tangential_coeffs_25, thin_prism_coeffs_25);
        float2  _S2573 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_8 = length_0(_S2573);
        float _S2574 = vert0_c_1.z;
        float theta_4 = (F32_atan2((r_8), (_S2574)));
        float k_4;
        if(theta_4 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_4 * theta_4 / 3.0f) / _S2574;
        }
        else
        {
            k_4 = theta_4 / r_8;
        }
        float2  _S2575 = _S2573 * make_float2 (k_4);
        float k1_3 = dist_coeffs_2.radial_coeffs_0.x;
        float k2_3 = dist_coeffs_2.radial_coeffs_0.y;
        float k3_3 = dist_coeffs_2.radial_coeffs_0.z;
        float k4_4 = dist_coeffs_2.radial_coeffs_0.w;
        float p1_5 = dist_coeffs_2.tangential_coeffs_0.x;
        float p2_5 = dist_coeffs_2.tangential_coeffs_0.y;
        float sx1_4 = dist_coeffs_2.thin_prism_coeffs_0.x;
        float sy1_4 = dist_coeffs_2.thin_prism_coeffs_0.y;
        float u_11 = _S2575.x;
        float v_11 = _S2575.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float _S2576 = 2.0f * p1_5;
        float _S2577 = 2.0f * p2_5;
        float2  _S2578 = _S2575 * make_float2 (1.0f + r2_11 * (k1_3 + r2_11 * (k2_3 + r2_11 * (k3_3 + r2_11 * k4_4)))) + make_float2 (_S2576 * u_11 * v_11 + p2_5 * (r2_11 + 2.0f * u_11 * u_11) + sx1_4 * r2_11, _S2577 * u_11 * v_11 + p1_5 * (r2_11 + 2.0f * v_11 * v_11) + sy1_4 * r2_11);
        *uv0_1 = make_float2 (fx_20 * _S2578.x + cx_20, fy_20 * _S2578.y + cy_20);
        float2  _S2579 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_9 = length_0(_S2579);
        float _S2580 = vert1_c_1.z;
        float theta_5 = (F32_atan2((r_9), (_S2580)));
        if(theta_5 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_5 * theta_5 / 3.0f) / _S2580;
        }
        else
        {
            k_4 = theta_5 / r_9;
        }
        float2  _S2581 = _S2579 * make_float2 (k_4);
        float u_12 = _S2581.x;
        float v_12 = _S2581.y;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float2  _S2582 = _S2581 * make_float2 (1.0f + r2_12 * (k1_3 + r2_12 * (k2_3 + r2_12 * (k3_3 + r2_12 * k4_4)))) + make_float2 (_S2576 * u_12 * v_12 + p2_5 * (r2_12 + 2.0f * u_12 * u_12) + sx1_4 * r2_12, _S2577 * u_12 * v_12 + p1_5 * (r2_12 + 2.0f * v_12 * v_12) + sy1_4 * r2_12);
        *uv1_1 = make_float2 (fx_20 * _S2582.x + cx_20, fy_20 * _S2582.y + cy_20);
        float2  _S2583 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_10 = length_0(_S2583);
        float _S2584 = vert2_c_1.z;
        float theta_6 = (F32_atan2((r_10), (_S2584)));
        if(theta_6 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_6 * theta_6 / 3.0f) / _S2584;
        }
        else
        {
            k_4 = theta_6 / r_10;
        }
        float2  _S2585 = _S2583 * make_float2 (k_4);
        float u_13 = _S2585.x;
        float v_13 = _S2585.y;
        float r2_13 = u_13 * u_13 + v_13 * v_13;
        float2  _S2586 = _S2585 * make_float2 (1.0f + r2_13 * (k1_3 + r2_13 * (k2_3 + r2_13 * (k3_3 + r2_13 * k4_4)))) + make_float2 (_S2576 * u_13 * v_13 + p2_5 * (r2_13 + 2.0f * u_13 * u_13) + sx1_4 * r2_13, _S2577 * u_13 * v_13 + p1_5 * (r2_13 + 2.0f * v_13 * v_13) + sy1_4 * r2_13);
        float _S2587 = fx_20 * _S2586.x + cx_20;
        float _S2588 = fy_20 * _S2586.y + cy_20;
        float2  _S2589 = make_float2 (_S2587, _S2588);
        *uv2_1 = _S2589;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S2589 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S2589)));
        float xmax_6 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S2587))) + offset_1;
        float xmin_6 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S2587))) - offset_1;
        float ymax_6 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S2588))) + offset_1;
        float ymin_6 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S2588))) - offset_1;
        if(xmax_6 <= 0.0f)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = xmin_6 >= float(image_width_16);
        }
        if(_S2566)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = ymax_6 <= 0.0f;
        }
        if(_S2566)
        {
            _S2566 = true;
        }
        else
        {
            _S2566 = ymin_6 >= float(image_height_16);
        }
        if(_S2566)
        {
            _S2566 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S2566 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S2566 = false;
                }
                if(_S2566)
                {
                    _S2566 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S2566 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S2566 = false;
                    }
                }
            }
            else
            {
                _S2566 = false;
            }
        }
        if(_S2566)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S2570, _S2571, _S2572) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S2590 = mean_20 - - mul_0(transpose_0(R_20), t_19);
        float _S2591 = _S2590.x;
        float _S2592 = _S2590.y;
        float _S2593 = _S2590.z;
        float norm_11 = (F32_sqrt((_S2591 * _S2591 + _S2592 * _S2592 + _S2593 * _S2593)));
        float x_44 = _S2591 / norm_11;
        float y_19 = _S2592 / norm_11;
        float z_16 = _S2593 / norm_11;
        float z2_38 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_44 * x_44 - y_19 * y_19;
        float fS1_16 = 2.0f * x_44 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_38 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_44) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_38 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_44) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_44 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_38 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_44) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_44 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S2594 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S2594);
        float3  _S2595 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S2596 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S2595 + _S2596 + make_float3 (0.5f), _S2594);
        (*rgb_13)[int(2)] = max_0(_S2595 - _S2596 + make_float3 (0.5f), _S2594);
        float3  _S2597 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S2597 * make_float3 (float(- (F32_sign((dot_0(_S2597, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_21, float3  t_20, float fx_21, float fy_21, float cx_21, float cy_21, float4  radial_coeffs_26, float2  tangential_coeffs_26, float2  thin_prism_coeffs_26, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_21, mean_21) + t_20;
    float _S2598 = scale_21.x;
    float sx_2 = (F32_exp((_S2598)));
    float _S2599 = scale_21.y;
    float sy_2 = (F32_exp((_S2599)));
    float sz_2 = scale_21.z - 0.5f * (_S2598 + _S2599);
    float x_45 = quat_22.y;
    float x2_22 = x_45 * x_45;
    float y2_22 = quat_22.z * quat_22.z;
    float z2_39 = quat_22.w * quat_22.w;
    float xy_22 = quat_22.y * quat_22.z;
    float xz_22 = quat_22.y * quat_22.w;
    float yz_22 = quat_22.z * quat_22.w;
    float wx_22 = quat_22.x * quat_22.y;
    float wy_22 = quat_22.x * quat_22.z;
    float wz_22 = quat_22.x * quat_22.w;
    Matrix<float, 3, 3>  _S2600 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_39), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_39), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
    float3  vert0_c_2 = mul_0(R_21, mul_0(_S2600, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_20;
    float3  vert1_c_2 = mul_0(R_21, mul_0(_S2600, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_20;
    float3  vert2_c_2 = mul_0(R_21, mul_0(_S2600, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_20;
    *uv0_2 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    *uv1_2 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    *uv2_2 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float2  _S2601 = make_float2 (fx_21, fy_21);
    float2  _S2602 = make_float2 (cx_21, cy_21);
    *uv0_2 = _S2601 * *uv0_2 + _S2602;
    *uv1_2 = _S2601 * *uv1_2 + _S2602;
    float2  _S2603 = _S2601 * *uv2_2 + _S2602;
    *uv2_2 = _S2603;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S2603 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S2603)));
    float _S2604 = _S2603.x;
    float _S2605 = _S2603.y;
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S2604))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S2605))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S2604))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S2605))) + offset_2)))));
    *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_2 = hardness_2;
    float3  _S2606 = mean_21 - - mul_0(transpose_0(R_21), t_20);
    float _S2607 = _S2606.x;
    float _S2608 = _S2606.y;
    float _S2609 = _S2606.z;
    float norm_12 = (F32_sqrt((_S2607 * _S2607 + _S2608 * _S2608 + _S2609 * _S2609)));
    float x_46 = _S2607 / norm_12;
    float y_20 = _S2608 / norm_12;
    float z_17 = _S2609 / norm_12;
    float z2_40 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_46 * x_46 - y_20 * y_20;
    float fS1_17 = 2.0f * x_46 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_40 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_46) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_40 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_46) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_40 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_46) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S2610 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S2610);
    float3  _S2611 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S2612 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S2611 + _S2612 + make_float3 (0.5f), _S2610);
    (*rgb_14)[int(2)] = max_0(_S2611 - _S2612 + make_float3 (0.5f), _S2610);
    float3  _S2613 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S2613 * make_float3 (float(- (F32_sign((dot_0(_S2613, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_22, float3  t_21, float fx_22, float fy_22, float cx_22, float cy_22, float4  radial_coeffs_27, float2  tangential_coeffs_27, float2  thin_prism_coeffs_27, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_22, mean_22) + t_21;
    float _S2614 = scale_22.x;
    float sx_3 = (F32_exp((_S2614)));
    float _S2615 = scale_22.y;
    float sy_3 = (F32_exp((_S2615)));
    float sz_3 = scale_22.z - 0.5f * (_S2614 + _S2615);
    float x_47 = quat_23.y;
    float x2_23 = x_47 * x_47;
    float y2_23 = quat_23.z * quat_23.z;
    float z2_41 = quat_23.w * quat_23.w;
    float xy_23 = quat_23.y * quat_23.z;
    float xz_23 = quat_23.y * quat_23.w;
    float yz_23 = quat_23.z * quat_23.w;
    float wx_23 = quat_23.x * quat_23.y;
    float wy_23 = quat_23.x * quat_23.z;
    float wz_23 = quat_23.x * quat_23.w;
    Matrix<float, 3, 3>  _S2616 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_41), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_41), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_3 = mul_0(R_22, mul_0(_S2616, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_21;
    float3  vert1_c_3 = mul_0(R_22, mul_0(_S2616, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_21;
    float3  vert2_c_3 = mul_0(R_22, mul_0(_S2616, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_21;
    CameraDistortion_0 dist_coeffs_3 = CameraDistortion_x24init_0(radial_coeffs_27, tangential_coeffs_27, thin_prism_coeffs_27);
    float2  _S2617 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_11 = length_0(_S2617);
    float _S2618 = vert0_c_3.z;
    float theta_7 = (F32_atan2((r_11), (_S2618)));
    float k_5;
    if(theta_7 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_7 * theta_7 / 3.0f) / _S2618;
    }
    else
    {
        k_5 = theta_7 / r_11;
    }
    float2  _S2619 = _S2617 * make_float2 (k_5);
    float k1_4 = dist_coeffs_3.radial_coeffs_0.x;
    float k2_4 = dist_coeffs_3.radial_coeffs_0.y;
    float k3_4 = dist_coeffs_3.radial_coeffs_0.z;
    float k4_5 = dist_coeffs_3.radial_coeffs_0.w;
    float p1_6 = dist_coeffs_3.tangential_coeffs_0.x;
    float p2_6 = dist_coeffs_3.tangential_coeffs_0.y;
    float sx1_5 = dist_coeffs_3.thin_prism_coeffs_0.x;
    float sy1_5 = dist_coeffs_3.thin_prism_coeffs_0.y;
    float u_14 = _S2619.x;
    float v_14 = _S2619.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float _S2620 = 2.0f * p1_6;
    float _S2621 = 2.0f * p2_6;
    float2  _S2622 = _S2619 * make_float2 (1.0f + r2_14 * (k1_4 + r2_14 * (k2_4 + r2_14 * (k3_4 + r2_14 * k4_5)))) + make_float2 (_S2620 * u_14 * v_14 + p2_6 * (r2_14 + 2.0f * u_14 * u_14) + sx1_5 * r2_14, _S2621 * u_14 * v_14 + p1_6 * (r2_14 + 2.0f * v_14 * v_14) + sy1_5 * r2_14);
    *uv0_3 = make_float2 (fx_22 * _S2622.x + cx_22, fy_22 * _S2622.y + cy_22);
    float2  _S2623 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_12 = length_0(_S2623);
    float _S2624 = vert1_c_3.z;
    float theta_8 = (F32_atan2((r_12), (_S2624)));
    if(theta_8 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_8 * theta_8 / 3.0f) / _S2624;
    }
    else
    {
        k_5 = theta_8 / r_12;
    }
    float2  _S2625 = _S2623 * make_float2 (k_5);
    float u_15 = _S2625.x;
    float v_15 = _S2625.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S2626 = _S2625 * make_float2 (1.0f + r2_15 * (k1_4 + r2_15 * (k2_4 + r2_15 * (k3_4 + r2_15 * k4_5)))) + make_float2 (_S2620 * u_15 * v_15 + p2_6 * (r2_15 + 2.0f * u_15 * u_15) + sx1_5 * r2_15, _S2621 * u_15 * v_15 + p1_6 * (r2_15 + 2.0f * v_15 * v_15) + sy1_5 * r2_15);
    *uv1_3 = make_float2 (fx_22 * _S2626.x + cx_22, fy_22 * _S2626.y + cy_22);
    float2  _S2627 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_13 = length_0(_S2627);
    float _S2628 = vert2_c_3.z;
    float theta_9 = (F32_atan2((r_13), (_S2628)));
    if(theta_9 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_9 * theta_9 / 3.0f) / _S2628;
    }
    else
    {
        k_5 = theta_9 / r_13;
    }
    float2  _S2629 = _S2627 * make_float2 (k_5);
    float u_16 = _S2629.x;
    float v_16 = _S2629.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float2  _S2630 = _S2629 * make_float2 (1.0f + r2_16 * (k1_4 + r2_16 * (k2_4 + r2_16 * (k3_4 + r2_16 * k4_5)))) + make_float2 (_S2620 * u_16 * v_16 + p2_6 * (r2_16 + 2.0f * u_16 * u_16) + sx1_5 * r2_16, _S2621 * u_16 * v_16 + p1_6 * (r2_16 + 2.0f * v_16 * v_16) + sy1_5 * r2_16);
    float _S2631 = fx_22 * _S2630.x + cx_22;
    float _S2632 = fy_22 * _S2630.y + cy_22;
    float2  _S2633 = make_float2 (_S2631, _S2632);
    *uv2_3 = _S2633;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S2633 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S2633)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S2631))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S2632))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S2631))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S2632))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S2634 = mean_22 - - mul_0(transpose_0(R_22), t_21);
    float _S2635 = _S2634.x;
    float _S2636 = _S2634.y;
    float _S2637 = _S2634.z;
    float norm_13 = (F32_sqrt((_S2635 * _S2635 + _S2636 * _S2636 + _S2637 * _S2637)));
    float x_48 = _S2635 / norm_13;
    float y_21 = _S2636 / norm_13;
    float z_18 = _S2637 / norm_13;
    float z2_42 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_48 * x_48 - y_21 * y_21;
    float fS1_18 = 2.0f * x_48 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_42 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_48) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_42 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_48) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_48 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_42 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_48) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_48 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S2638 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S2638);
    float3  _S2639 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S2640 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S2639 + _S2640 + make_float3 (0.5f), _S2638);
    (*rgb_15)[int(2)] = max_0(_S2639 - _S2640 + make_float3 (0.5f), _S2638);
    float3  _S2641 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S2641 * make_float3 (float(- (F32_sign((dot_0(_S2641, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2642, float3  _S2643)
{
    _d_log_vector_0(_S2642, _S2643);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S2644, float _S2645)
{
    _d_exp2_0(_S2644, _S2645);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S2646, float _S2647)
{
    _d_abs_0(_S2646, _S2647);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_23, float3  t_22, float fx_23, float fy_23, float cx_23, float cy_23, float4  radial_coeffs_28, float2  tangential_coeffs_28, float2  thin_prism_coeffs_28, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_1(R_23, mean_23) + t_22;
    float _S2648 = scale_23.x;
    float _S2649 = s_primal_ctx_exp_1(_S2648);
    float _S2650 = scale_23.y;
    float _S2651 = s_primal_ctx_exp_1(_S2650);
    float sz_4 = scale_23.z - 0.5f * (_S2648 + _S2650);
    float _S2652 = quat_24.y;
    float x2_24 = _S2652 * _S2652;
    float y2_24 = quat_24.z * quat_24.z;
    float z2_43 = quat_24.w * quat_24.w;
    float xy_24 = quat_24.y * quat_24.z;
    float xz_24 = quat_24.y * quat_24.w;
    float yz_24 = quat_24.z * quat_24.w;
    float wx_24 = quat_24.x * quat_24.y;
    float wy_24 = quat_24.x * quat_24.z;
    float wz_24 = quat_24.x * quat_24.w;
    Matrix<float, 3, 3>  _S2653 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_43), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_43), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  _S2654 = make_float3 (_S2649, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_1(_S2653, _S2654) + mean_23;
    float _S2655 = -0.5f + sz_4;
    float3  _S2656 = make_float3 (_S2649 * _S2655, _S2651, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_1(_S2653, _S2656) + mean_23;
    float _S2657 = -0.5f - sz_4;
    float3  _S2658 = make_float3 (_S2649 * _S2657, - _S2651, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_1(_S2653, _S2658) + mean_23;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_23, vert0_0) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_23, vert1_0) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_23, vert2_0) + t_22;
    float2  _S2659 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S2660 = vert0_c_4.z;
    float2  _S2661 = make_float2 (_S2660);
    float2  _S2662 = make_float2 (_S2660 * _S2660);
    float2  _S2663 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S2664 = vert1_c_4.z;
    float2  _S2665 = make_float2 (_S2664);
    float2  _S2666 = make_float2 (_S2664 * _S2664);
    float2  _S2667 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S2668 = vert2_c_4.z;
    float2  _S2669 = make_float2 (_S2668);
    float2  _S2670 = make_float2 (_S2668 * _S2668);
    float2  _S2671 = make_float2 (fx_23, fy_23);
    float2  _S2672 = make_float2 (cx_23, cy_23);
    float2  _S2673 = _S2671 * (_S2659 / make_float2 (_S2660)) + _S2672;
    float2  _S2674 = _S2671 * (_S2663 / make_float2 (_S2664)) + _S2672;
    float2  _S2675 = _S2671 * (_S2667 / make_float2 (_S2668)) + _S2672;
    float2  e0_4 = _S2674 - _S2673;
    float2  e1_4 = _S2675 - _S2674;
    float2  e2_0 = _S2673 - _S2675;
    float _S2676 = e0_4.x;
    float _S2677 = e1_4.y;
    float _S2678 = e0_4.y;
    float _S2679 = e1_4.x;
    float _S2680 = _S2676 * _S2677 - _S2678 * _S2679;
    float _S2681 = 1.0f - hardness_4.y;
    float _S2682 = -1.0f / _S2681;
    float _S2683 = _S2681 * _S2681;
    float _S2684 = _S2673.x;
    float _S2685 = _S2674.x;
    float _S2686 = s_primal_ctx_max_0(_S2684, _S2685);
    float _S2687 = _S2675.x;
    float _S2688 = s_primal_ctx_min_0(_S2684, _S2685);
    float _S2689 = _S2673.y;
    float _S2690 = _S2674.y;
    float _S2691 = s_primal_ctx_max_0(_S2689, _S2690);
    float _S2692 = _S2675.y;
    float _S2693 = s_primal_ctx_min_0(_S2689, _S2690);
    float3  _S2694 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S2695 = transpose_0(R_23);
    float3  _S2696 = mean_23 - - s_primal_ctx_mul_1(_S2695, t_22);
    float _S2697 = _S2696.x;
    float _S2698 = _S2696.y;
    float _S2699 = _S2696.z;
    float _S2700 = _S2697 * _S2697 + _S2698 * _S2698 + _S2699 * _S2699;
    float _S2701 = s_primal_ctx_sqrt_0(_S2700);
    float x_49 = _S2697 / _S2701;
    float3  _S2702 = make_float3 (x_49);
    float _S2703 = _S2701 * _S2701;
    float y_22 = _S2698 / _S2701;
    float z_19 = _S2699 / _S2701;
    float3  _S2704 = make_float3 (z_19);
    float _S2705 = - y_22;
    float3  _S2706 = make_float3 (_S2705);
    float z2_44 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_49 * x_49 - y_22 * y_22;
    float _S2707 = 2.0f * x_49;
    float fS1_19 = _S2707 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_44 - 0.31539157032966614f;
    float3  _S2708 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_49;
    float3  _S2709 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S2710 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S2711 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S2712 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_44 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S2713 = 1.86588168144226074f * z2_44 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S2713;
    float3  _S2714 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_49;
    float3  _S2715 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S2716 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S2717 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S2718 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_49 * fC1_19 - y_22 * fS1_19);
    float3  _S2719 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_49 * fS1_19 + y_22 * fC1_19);
    float3  _S2720 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2705) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_49) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S2721 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S2722 = make_float3 (0.0f);
    float3  _S2723 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S2724 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S2725 = make_float3 (_S2724);
    float3  _S2726 = (*ch_coeffs_4)[int(1)] * make_float3 (_S2724);
    float3  _S2727 = _S2723 + _S2726 + make_float3 (0.5f);
    float3  _S2728 = _S2723 - _S2726 + make_float3 (0.5f);
    float3  _S2729 = vert1_c_4 - vert0_c_4;
    float3  _S2730 = vert2_c_4 - vert0_c_4;
    float3  _S2731 = s_primal_ctx_cross_0(_S2729, _S2730);
    float3  _S2732 = normalize_0(_S2731);
    float3  _S2733 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S2732, mean_c_19)))))) * v_normal_0;
    float3  _S2734 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2735;
    (&_S2735)->primal_0 = _S2732;
    (&_S2735)->differential_0 = _S2734;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2736;
    (&_S2736)->primal_0 = mean_c_19;
    (&_S2736)->differential_0 = _S2734;
    s_bwd_prop_dot_0(&_S2735, &_S2736, 0.0f);
    float3  _S2737 = _S2733 + _S2735.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2738;
    (&_S2738)->primal_0 = _S2731;
    (&_S2738)->differential_0 = _S2734;
    s_bwd_normalize_impl_0(&_S2738, _S2737);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2739;
    (&_S2739)->primal_0 = _S2729;
    (&_S2739)->differential_0 = _S2734;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2740;
    (&_S2740)->primal_0 = _S2730;
    (&_S2740)->differential_0 = _S2734;
    s_bwd_prop_cross_0(&_S2739, &_S2740, _S2738.differential_0);
    float3  _S2741 = - _S2740.differential_0;
    float3  _S2742 = - _S2739.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2743;
    (&_S2743)->primal_0 = _S2728;
    (&_S2743)->differential_0 = _S2734;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2744;
    (&_S2744)->primal_0 = _S2722;
    (&_S2744)->differential_0 = _S2734;
    s_bwd_prop_max_0(&_S2743, &_S2744, (*v_rgb_6)[int(2)]);
    float3  _S2745 = - _S2743.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2746;
    (&_S2746)->primal_0 = _S2727;
    (&_S2746)->differential_0 = _S2734;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2747;
    (&_S2747)->primal_0 = _S2722;
    (&_S2747)->differential_0 = _S2734;
    s_bwd_prop_max_0(&_S2746, &_S2747, (*v_rgb_6)[int(1)]);
    float3  _S2748 = _S2725 * (_S2745 + _S2746.differential_0);
    float3  _S2749 = _S2743.differential_0 + _S2746.differential_0;
    float3  _S2750 = make_float3 (0.5f) * - _S2749;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2751;
    (&_S2751)->primal_0 = _S2721;
    (&_S2751)->differential_0 = _S2734;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2752;
    (&_S2752)->primal_0 = _S2722;
    (&_S2752)->differential_0 = _S2734;
    s_bwd_prop_max_0(&_S2751, &_S2752, (*v_rgb_6)[int(0)]);
    float3  _S2753 = _S2750 + _S2751.differential_0;
    float3  _S2754 = _S2749 + _S2751.differential_0;
    float3  _S2755 = _S2719 * _S2754;
    float3  _S2756 = (*sh_coeffs_19)[int(15)] * _S2754;
    float3  _S2757 = _S2717 * _S2754;
    float3  _S2758 = (*sh_coeffs_19)[int(14)] * _S2754;
    float3  _S2759 = _S2715 * _S2754;
    float3  _S2760 = (*sh_coeffs_19)[int(13)] * _S2754;
    float3  _S2761 = _S2714 * _S2754;
    float3  _S2762 = (*sh_coeffs_19)[int(12)] * _S2754;
    float3  _S2763 = _S2716 * _S2754;
    float3  _S2764 = (*sh_coeffs_19)[int(11)] * _S2754;
    float3  _S2765 = _S2718 * _S2754;
    float3  _S2766 = (*sh_coeffs_19)[int(10)] * _S2754;
    float3  _S2767 = _S2720 * _S2754;
    float3  _S2768 = (*sh_coeffs_19)[int(9)] * _S2754;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S2768.x + _S2768.y + _S2768.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S2756.x + _S2756.y + _S2756.z);
    float _S2769 = _S2766.x + _S2766.y + _S2766.z;
    float _S2770 = _S2758.x + _S2758.y + _S2758.z;
    float _S2771 = _S2764.x + _S2764.y + _S2764.z;
    float _S2772 = _S2760.x + _S2760.y + _S2760.z;
    float _S2773 = _S2762.x + _S2762.y + _S2762.z;
    float _S2774 = - s_diff_fC2_T_5;
    float3  _S2775 = _S2711 * _S2754;
    float3  _S2776 = (*sh_coeffs_19)[int(8)] * _S2754;
    float3  _S2777 = _S2709 * _S2754;
    float3  _S2778 = (*sh_coeffs_19)[int(7)] * _S2754;
    float3  _S2779 = _S2708 * _S2754;
    float3  _S2780 = (*sh_coeffs_19)[int(6)] * _S2754;
    float3  _S2781 = _S2710 * _S2754;
    float3  _S2782 = (*sh_coeffs_19)[int(5)] * _S2754;
    float3  _S2783 = _S2712 * _S2754;
    float3  _S2784 = (*sh_coeffs_19)[int(4)] * _S2754;
    float _S2785 = _S2782.x + _S2782.y + _S2782.z;
    float _S2786 = _S2778.x + _S2778.y + _S2778.z;
    float _S2787 = fTmp1B_19 * _S2769 + x_49 * s_diff_fS2_T_5 + y_22 * _S2774 + 0.54627424478530884f * (_S2784.x + _S2784.y + _S2784.z);
    float _S2788 = fTmp1B_19 * _S2770 + y_22 * s_diff_fS2_T_5 + x_49 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S2776.x + _S2776.y + _S2776.z);
    float _S2789 = y_22 * - _S2788;
    float _S2790 = x_49 * _S2788;
    float _S2791 = z_19 * (1.86588168144226074f * (z_19 * _S2773) + -2.28522896766662598f * (y_22 * _S2771 + x_49 * _S2772) + 0.94617468118667603f * (_S2780.x + _S2780.y + _S2780.z));
    float3  _S2792 = make_float3 (0.48860251903533936f) * _S2754;
    float3  _S2793 = - _S2792;
    float3  _S2794 = _S2702 * _S2793;
    float3  _S2795 = (*sh_coeffs_19)[int(3)] * _S2793;
    float3  _S2796 = _S2704 * _S2792;
    float3  _S2797 = (*sh_coeffs_19)[int(2)] * _S2792;
    float3  _S2798 = _S2706 * _S2792;
    float3  _S2799 = (*sh_coeffs_19)[int(1)] * _S2792;
    float _S2800 = (_S2713 * _S2773 + 1.44530570507049561f * (fS1_19 * _S2769 + fC1_19 * _S2770) + -1.09254848957061768f * (y_22 * _S2785 + x_49 * _S2786) + _S2791 + _S2791 + _S2797.x + _S2797.y + _S2797.z) / _S2703;
    float _S2801 = _S2701 * _S2800;
    float _S2802 = (fTmp0C_19 * _S2771 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S2774 + fTmp0B_19 * _S2785 + _S2707 * _S2787 + _S2789 + _S2789 + - (_S2799.x + _S2799.y + _S2799.z)) / _S2703;
    float _S2803 = _S2701 * _S2802;
    float _S2804 = (fTmp0C_19 * _S2772 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S2786 + 2.0f * (y_22 * _S2787) + _S2790 + _S2790 + _S2795.x + _S2795.y + _S2795.z) / _S2703;
    float _S2805 = _S2701 * _S2804;
    float _S2806 = _S2699 * - _S2800 + _S2698 * - _S2802 + _S2697 * - _S2804;
    DiffPair_float_0 _S2807;
    (&_S2807)->primal_0 = _S2700;
    (&_S2807)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2807, _S2806);
    float _S2808 = _S2699 * _S2807.differential_0;
    float _S2809 = _S2698 * _S2807.differential_0;
    float _S2810 = _S2697 * _S2807.differential_0;
    float3  _S2811 = make_float3 (0.282094806432724f) * _S2754;
    float3  _S2812 = make_float3 (_S2805 + _S2810 + _S2810, _S2803 + _S2809 + _S2809, _S2801 + _S2808 + _S2808);
    float3  _S2813 = - - _S2812;
    Matrix<float, 3, 3>  _S2814 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2815;
    (&_S2815)->primal_0 = _S2695;
    (&_S2815)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2816;
    (&_S2816)->primal_0 = t_22;
    (&_S2816)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2815, &_S2816, _S2813);
    Matrix<float, 3, 3>  _S2817 = transpose_0(_S2815.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2818;
    (&_S2818)->primal_0 = _S2694;
    (&_S2818)->differential_0 = _S2734;
    s_bwd_prop_log_1(&_S2818, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2819;
    (&_S2819)->primal_0 = vert2_c_4;
    (&_S2819)->differential_0 = _S2734;
    s_bwd_length_impl_0(&_S2819, _S2818.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2820;
    (&_S2820)->primal_0 = vert1_c_4;
    (&_S2820)->differential_0 = _S2734;
    s_bwd_length_impl_0(&_S2820, _S2818.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2821;
    (&_S2821)->primal_0 = vert0_c_4;
    (&_S2821)->differential_0 = _S2734;
    s_bwd_length_impl_0(&_S2821, _S2818.differential_0.x);
    DiffPair_float_0 _S2822;
    (&_S2822)->primal_0 = _S2693;
    (&_S2822)->differential_0 = 0.0f;
    DiffPair_float_0 _S2823;
    (&_S2823)->primal_0 = _S2692;
    (&_S2823)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2822, &_S2823, 0.0f);
    DiffPair_float_0 _S2824;
    (&_S2824)->primal_0 = _S2689;
    (&_S2824)->differential_0 = 0.0f;
    DiffPair_float_0 _S2825;
    (&_S2825)->primal_0 = _S2690;
    (&_S2825)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2824, &_S2825, _S2822.differential_0);
    DiffPair_float_0 _S2826;
    (&_S2826)->primal_0 = _S2691;
    (&_S2826)->differential_0 = 0.0f;
    DiffPair_float_0 _S2827;
    (&_S2827)->primal_0 = _S2692;
    (&_S2827)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2826, &_S2827, 0.0f);
    float _S2828 = _S2823.differential_0 + _S2827.differential_0;
    DiffPair_float_0 _S2829;
    (&_S2829)->primal_0 = _S2689;
    (&_S2829)->differential_0 = 0.0f;
    DiffPair_float_0 _S2830;
    (&_S2830)->primal_0 = _S2690;
    (&_S2830)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2829, &_S2830, _S2826.differential_0);
    float _S2831 = _S2825.differential_0 + _S2830.differential_0;
    float _S2832 = _S2824.differential_0 + _S2829.differential_0;
    DiffPair_float_0 _S2833;
    (&_S2833)->primal_0 = _S2688;
    (&_S2833)->differential_0 = 0.0f;
    DiffPair_float_0 _S2834;
    (&_S2834)->primal_0 = _S2687;
    (&_S2834)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2833, &_S2834, 0.0f);
    DiffPair_float_0 _S2835;
    (&_S2835)->primal_0 = _S2684;
    (&_S2835)->differential_0 = 0.0f;
    DiffPair_float_0 _S2836;
    (&_S2836)->primal_0 = _S2685;
    (&_S2836)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2835, &_S2836, _S2833.differential_0);
    DiffPair_float_0 _S2837;
    (&_S2837)->primal_0 = _S2686;
    (&_S2837)->differential_0 = 0.0f;
    DiffPair_float_0 _S2838;
    (&_S2838)->primal_0 = _S2687;
    (&_S2838)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2837, &_S2838, 0.0f);
    float _S2839 = _S2834.differential_0 + _S2838.differential_0;
    DiffPair_float_0 _S2840;
    (&_S2840)->primal_0 = _S2684;
    (&_S2840)->differential_0 = 0.0f;
    DiffPair_float_0 _S2841;
    (&_S2841)->primal_0 = _S2685;
    (&_S2841)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2840, &_S2841, _S2837.differential_0);
    float _S2842 = _S2836.differential_0 + _S2841.differential_0;
    float _S2843 = _S2835.differential_0 + _S2840.differential_0;
    DiffPair_float_0 _S2844;
    (&_S2844)->primal_0 = _S2682;
    (&_S2844)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2844, 0.0f);
    float _S2845 = - (-1.0f * - (_S2844.differential_0 / _S2683));
    float2  _S2846 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2847;
    (&_S2847)->primal_0 = e2_0;
    (&_S2847)->differential_0 = _S2846;
    s_bwd_length_impl_1(&_S2847, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2848;
    (&_S2848)->primal_0 = e1_4;
    (&_S2848)->differential_0 = _S2846;
    s_bwd_length_impl_1(&_S2848, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2849;
    (&_S2849)->primal_0 = e0_4;
    (&_S2849)->differential_0 = _S2846;
    s_bwd_length_impl_1(&_S2849, -0.0f);
    DiffPair_float_0 _S2850;
    (&_S2850)->primal_0 = _S2680;
    (&_S2850)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2850, 0.0f);
    float _S2851 = - _S2850.differential_0;
    float2  _S2852 = _S2848.differential_0 + make_float2 (_S2678 * _S2851, _S2676 * _S2850.differential_0);
    float2  _S2853 = _S2849.differential_0 + make_float2 (_S2677 * _S2850.differential_0, _S2679 * _S2851);
    float2  _S2854 = _S2671 * (v_uv2_0 + - _S2847.differential_0 + _S2852 + make_float2 (_S2839, _S2828)) / _S2670;
    float2  _S2855 = _S2667 * - _S2854;
    float2  _S2856 = _S2669 * _S2854;
    float2  _S2857 = _S2671 * (v_uv1_0 + - _S2852 + _S2853 + make_float2 (_S2842, _S2831)) / _S2666;
    float2  _S2858 = _S2663 * - _S2857;
    float2  _S2859 = _S2665 * _S2857;
    float _S2860 = _S2858.x + _S2858.y;
    float2  _S2861 = _S2671 * (v_uv0_0 + _S2847.differential_0 + - _S2853 + make_float2 (_S2843, _S2832)) / _S2662;
    float2  _S2862 = _S2659 * - _S2861;
    float2  _S2863 = _S2661 * _S2861;
    float _S2864 = _S2862.x + _S2862.y;
    float3  _S2865 = _S2740.differential_0 + _S2819.differential_0 + make_float3 (_S2856.x, _S2856.y, _S2855.x + _S2855.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2866;
    (&_S2866)->primal_0 = R_23;
    (&_S2866)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2867;
    (&_S2867)->primal_0 = vert2_0;
    (&_S2867)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2866, &_S2867, _S2865);
    float3  _S2868 = _S2739.differential_0 + _S2820.differential_0 + make_float3 (_S2859.x, _S2859.y, _S2860);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2869;
    (&_S2869)->primal_0 = R_23;
    (&_S2869)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2870;
    (&_S2870)->primal_0 = vert1_0;
    (&_S2870)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2869, &_S2870, _S2868);
    float3  _S2871 = _S2741 + _S2742 + _S2821.differential_0 + make_float3 (_S2863.x, _S2863.y, _S2864);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2872;
    (&_S2872)->primal_0 = R_23;
    (&_S2872)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2873;
    (&_S2873)->primal_0 = vert0_0;
    (&_S2873)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2872, &_S2873, _S2871);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2874;
    (&_S2874)->primal_0 = _S2653;
    (&_S2874)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2875;
    (&_S2875)->primal_0 = _S2658;
    (&_S2875)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2874, &_S2875, _S2867.differential_0);
    float _S2876 = - _S2875.differential_0.y;
    float _S2877 = _S2657 * _S2875.differential_0.x;
    float _S2878 = - (_S2649 * _S2875.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2879;
    (&_S2879)->primal_0 = _S2653;
    (&_S2879)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2880;
    (&_S2880)->primal_0 = _S2656;
    (&_S2880)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2879, &_S2880, _S2870.differential_0);
    float _S2881 = _S2649 * _S2880.differential_0.x;
    float _S2882 = _S2655 * _S2880.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2883;
    (&_S2883)->primal_0 = _S2653;
    (&_S2883)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2884;
    (&_S2884)->primal_0 = _S2654;
    (&_S2884)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2883, &_S2884, _S2873.differential_0);
    Matrix<float, 3, 3>  _S2885 = transpose_0(_S2874.differential_0 + _S2879.differential_0 + _S2883.differential_0);
    float _S2886 = 2.0f * - _S2885.rows[int(2)].z;
    float _S2887 = 2.0f * _S2885.rows[int(2)].y;
    float _S2888 = 2.0f * _S2885.rows[int(2)].x;
    float _S2889 = 2.0f * _S2885.rows[int(1)].z;
    float _S2890 = 2.0f * - _S2885.rows[int(1)].y;
    float _S2891 = 2.0f * _S2885.rows[int(1)].x;
    float _S2892 = 2.0f * _S2885.rows[int(0)].z;
    float _S2893 = 2.0f * _S2885.rows[int(0)].y;
    float _S2894 = 2.0f * - _S2885.rows[int(0)].x;
    float _S2895 = - _S2891 + _S2893;
    float _S2896 = _S2888 + - _S2892;
    float _S2897 = - _S2887 + _S2889;
    float _S2898 = _S2887 + _S2889;
    float _S2899 = _S2888 + _S2892;
    float _S2900 = _S2891 + _S2893;
    float _S2901 = quat_24.w * (_S2890 + _S2894);
    float _S2902 = quat_24.z * (_S2886 + _S2894);
    float _S2903 = quat_24.y * (_S2886 + _S2890);
    float _S2904 = quat_24.x * _S2895 + quat_24.z * _S2898 + quat_24.y * _S2899 + _S2901 + _S2901;
    float _S2905 = quat_24.x * _S2896 + quat_24.w * _S2898 + quat_24.y * _S2900 + _S2902 + _S2902;
    float _S2906 = quat_24.x * _S2897 + quat_24.w * _S2899 + quat_24.z * _S2900 + _S2903 + _S2903;
    float _S2907 = quat_24.w * _S2895 + quat_24.z * _S2896 + quat_24.y * _S2897;
    float _S2908 = _S2878 + _S2881;
    float _S2909 = 0.5f * - _S2908;
    float _S2910 = _S2876 + _S2880.differential_0.y;
    DiffPair_float_0 _S2911;
    (&_S2911)->primal_0 = _S2650;
    (&_S2911)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2911, _S2910);
    float _S2912 = _S2909 + _S2911.differential_0;
    float _S2913 = _S2877 + _S2882 + _S2884.differential_0.x;
    DiffPair_float_0 _S2914;
    (&_S2914)->primal_0 = _S2648;
    (&_S2914)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2914, _S2913);
    float _S2915 = _S2909 + _S2914.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2916;
    (&_S2916)->primal_0 = R_23;
    (&_S2916)->differential_0 = _S2814;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2917;
    (&_S2917)->primal_0 = mean_23;
    (&_S2917)->differential_0 = _S2734;
    s_bwd_prop_mul_1(&_S2916, &_S2917, _S2736.differential_0);
    float3  _S2918 = _S2816.differential_0 + _S2865 + _S2868 + _S2871 + _S2736.differential_0;
    Matrix<float, 3, 3>  _S2919 = _S2817 + _S2866.differential_0 + _S2869.differential_0 + _S2872.differential_0 + _S2916.differential_0;
    FixedArray<float3 , 2>  _S2920;
    _S2920[int(0)] = _S2734;
    _S2920[int(1)] = _S2734;
    _S2920[int(1)] = _S2748;
    _S2920[int(0)] = _S2753;
    FixedArray<float3 , 16>  _S2921;
    _S2921[int(0)] = _S2734;
    _S2921[int(1)] = _S2734;
    _S2921[int(2)] = _S2734;
    _S2921[int(3)] = _S2734;
    _S2921[int(4)] = _S2734;
    _S2921[int(5)] = _S2734;
    _S2921[int(6)] = _S2734;
    _S2921[int(7)] = _S2734;
    _S2921[int(8)] = _S2734;
    _S2921[int(9)] = _S2734;
    _S2921[int(10)] = _S2734;
    _S2921[int(11)] = _S2734;
    _S2921[int(12)] = _S2734;
    _S2921[int(13)] = _S2734;
    _S2921[int(14)] = _S2734;
    _S2921[int(15)] = _S2734;
    _S2921[int(15)] = _S2755;
    _S2921[int(14)] = _S2757;
    _S2921[int(13)] = _S2759;
    _S2921[int(12)] = _S2761;
    _S2921[int(11)] = _S2763;
    _S2921[int(10)] = _S2765;
    _S2921[int(9)] = _S2767;
    _S2921[int(8)] = _S2775;
    _S2921[int(7)] = _S2777;
    _S2921[int(6)] = _S2779;
    _S2921[int(5)] = _S2781;
    _S2921[int(4)] = _S2783;
    _S2921[int(3)] = _S2794;
    _S2921[int(2)] = _S2796;
    _S2921[int(1)] = _S2798;
    _S2921[int(0)] = _S2811;
    float2  _S2922 = v_out_hardness_0 + make_float2 (0.0f, _S2845);
    float3  _S2923 = make_float3 (_S2915, _S2912, _S2908);
    float4  _S2924 = make_float4 (0.0f);
    *&((&_S2924)->w) = _S2904;
    *&((&_S2924)->z) = _S2905;
    *&((&_S2924)->y) = _S2906;
    *&((&_S2924)->x) = _S2907;
    *v_mean_7 = _S2812 + _S2867.differential_0 + _S2870.differential_0 + _S2873.differential_0 + _S2917.differential_0;
    *v_quat_6 = _S2924;
    *v_scale_6 = _S2923;
    *v_hardness_0 = _S2922;
    *v_sh_coeffs_5 = _S2921;
    *v_ch_coeffs_0 = _S2920;
    *v_R_6 = _S2919;
    *v_t_6 = _S2918;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_24, float3  t_23, float fx_24, float fy_24, float cx_24, float cy_24, float4  radial_coeffs_29, float2  tangential_coeffs_29, float2  thin_prism_coeffs_29, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_1(R_24, mean_24) + t_23;
    float _S2925 = scale_24.x;
    float _S2926 = s_primal_ctx_exp_1(_S2925);
    float _S2927 = scale_24.y;
    float _S2928 = s_primal_ctx_exp_1(_S2927);
    float sz_5 = scale_24.z - 0.5f * (_S2925 + _S2927);
    float _S2929 = quat_25.y;
    float x2_25 = _S2929 * _S2929;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_45 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S2930 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_45), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_45), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S2931 = make_float3 (_S2926, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S2930, _S2931) + mean_24;
    float _S2932 = -0.5f + sz_5;
    float3  _S2933 = make_float3 (_S2926 * _S2932, _S2928, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S2930, _S2933) + mean_24;
    float _S2934 = -0.5f - sz_5;
    float3  _S2935 = make_float3 (_S2926 * _S2934, - _S2928, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S2930, _S2935) + mean_24;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_24, vert0_1) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_24, vert1_1) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_24, vert2_1) + t_23;
    float _S2936 = length_1(vert0_c_5);
    float _S2937 = length_1(vert1_c_5);
    float _S2938 = length_1(vert2_c_5);
    CameraDistortion_0 _S2939 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_29, tangential_coeffs_29, thin_prism_coeffs_29);
    float2  _S2940 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S2941 = length_0(_S2940);
    float _S2942 = vert0_c_5.z;
    float _S2943 = s_primal_ctx_atan2_0(_S2941, _S2942);
    bool _S2944 = _S2943 < 0.00100000004749745f;
    float k_6;
    float _S2945;
    float _S2946;
    float _S2947;
    if(_S2944)
    {
        float _S2948 = 1.0f - _S2943 * _S2943 / 3.0f;
        float _S2949 = _S2942 * _S2942;
        k_6 = _S2948 / _S2942;
        _S2945 = 0.0f;
        _S2946 = _S2949;
        _S2947 = _S2948;
    }
    else
    {
        float _S2950 = _S2941 * _S2941;
        k_6 = _S2943 / _S2941;
        _S2945 = _S2950;
        _S2946 = 0.0f;
        _S2947 = 0.0f;
    }
    float2  _S2951 = make_float2 (k_6);
    float2  _S2952 = _S2940 * make_float2 (k_6);
    float k1_5 = _S2939.radial_coeffs_0.x;
    float k2_5 = _S2939.radial_coeffs_0.y;
    float k3_5 = _S2939.radial_coeffs_0.z;
    float k4_6 = _S2939.radial_coeffs_0.w;
    float p1_7 = _S2939.tangential_coeffs_0.x;
    float p2_7 = _S2939.tangential_coeffs_0.y;
    float sx1_6 = _S2939.thin_prism_coeffs_0.x;
    float sy1_6 = _S2939.thin_prism_coeffs_0.y;
    float u_17 = _S2952.x;
    float v_17 = _S2952.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S2953 = k3_5 + r2_17 * k4_6;
    float _S2954 = k2_5 + r2_17 * _S2953;
    float _S2955 = k1_5 + r2_17 * _S2954;
    float radial_2 = 1.0f + r2_17 * _S2955;
    float2  _S2956 = make_float2 (radial_2);
    float _S2957 = 2.0f * p1_7;
    float _S2958 = _S2957 * u_17;
    float _S2959 = 2.0f * u_17;
    float _S2960 = r2_17 + _S2959 * u_17;
    float _S2961 = 2.0f * p2_7;
    float _S2962 = _S2961 * u_17;
    float _S2963 = 2.0f * v_17;
    float _S2964 = r2_17 + _S2963 * v_17;
    float2  _S2965 = _S2952 * make_float2 (radial_2) + make_float2 (_S2958 * v_17 + p2_7 * _S2960 + sx1_6 * r2_17, _S2962 * v_17 + p1_7 * _S2964 + sy1_6 * r2_17);
    float _S2966 = fx_24 * _S2965.x + cx_24;
    float _S2967 = fy_24 * _S2965.y + cy_24;
    float2  _S2968 = make_float2 (_S2966, _S2967);
    float2  _S2969 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S2970 = length_0(_S2969);
    float _S2971 = vert1_c_5.z;
    float _S2972 = s_primal_ctx_atan2_0(_S2970, _S2971);
    bool _S2973 = _S2972 < 0.00100000004749745f;
    float _S2974;
    float _S2975;
    float _S2976;
    if(_S2973)
    {
        float _S2977 = 1.0f - _S2972 * _S2972 / 3.0f;
        float _S2978 = _S2971 * _S2971;
        k_6 = _S2977 / _S2971;
        _S2974 = 0.0f;
        _S2975 = _S2978;
        _S2976 = _S2977;
    }
    else
    {
        float _S2979 = _S2970 * _S2970;
        k_6 = _S2972 / _S2970;
        _S2974 = _S2979;
        _S2975 = 0.0f;
        _S2976 = 0.0f;
    }
    float2  _S2980 = make_float2 (k_6);
    float2  _S2981 = _S2969 * make_float2 (k_6);
    float u_18 = _S2981.x;
    float v_18 = _S2981.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S2982 = k3_5 + r2_18 * k4_6;
    float _S2983 = k2_5 + r2_18 * _S2982;
    float _S2984 = k1_5 + r2_18 * _S2983;
    float radial_3 = 1.0f + r2_18 * _S2984;
    float2  _S2985 = make_float2 (radial_3);
    float _S2986 = _S2957 * u_18;
    float _S2987 = 2.0f * u_18;
    float _S2988 = r2_18 + _S2987 * u_18;
    float _S2989 = _S2961 * u_18;
    float _S2990 = 2.0f * v_18;
    float _S2991 = r2_18 + _S2990 * v_18;
    float2  _S2992 = _S2981 * make_float2 (radial_3) + make_float2 (_S2986 * v_18 + p2_7 * _S2988 + sx1_6 * r2_18, _S2989 * v_18 + p1_7 * _S2991 + sy1_6 * r2_18);
    float _S2993 = fx_24 * _S2992.x + cx_24;
    float _S2994 = fy_24 * _S2992.y + cy_24;
    float2  _S2995 = make_float2 (_S2993, _S2994);
    float2  _S2996 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S2997 = length_0(_S2996);
    float _S2998 = vert2_c_5.z;
    float _S2999 = s_primal_ctx_atan2_0(_S2997, _S2998);
    bool _S3000 = _S2999 < 0.00100000004749745f;
    float _S3001;
    float _S3002;
    float _S3003;
    if(_S3000)
    {
        float _S3004 = 1.0f - _S2999 * _S2999 / 3.0f;
        float _S3005 = _S2998 * _S2998;
        k_6 = _S3004 / _S2998;
        _S3001 = 0.0f;
        _S3002 = _S3005;
        _S3003 = _S3004;
    }
    else
    {
        float _S3006 = _S2997 * _S2997;
        k_6 = _S2999 / _S2997;
        _S3001 = _S3006;
        _S3002 = 0.0f;
        _S3003 = 0.0f;
    }
    float2  _S3007 = make_float2 (k_6);
    float2  _S3008 = _S2996 * make_float2 (k_6);
    float u_19 = _S3008.x;
    float v_19 = _S3008.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S3009 = k3_5 + r2_19 * k4_6;
    float _S3010 = k2_5 + r2_19 * _S3009;
    float _S3011 = k1_5 + r2_19 * _S3010;
    float radial_4 = 1.0f + r2_19 * _S3011;
    float2  _S3012 = make_float2 (radial_4);
    float _S3013 = _S2957 * u_19;
    float _S3014 = 2.0f * u_19;
    float _S3015 = r2_19 + _S3014 * u_19;
    float _S3016 = _S2961 * u_19;
    float _S3017 = 2.0f * v_19;
    float _S3018 = r2_19 + _S3017 * v_19;
    float2  _S3019 = _S3008 * make_float2 (radial_4) + make_float2 (_S3013 * v_19 + p2_7 * _S3015 + sx1_6 * r2_19, _S3016 * v_19 + p1_7 * _S3018 + sy1_6 * r2_19);
    float _S3020 = fx_24 * _S3019.x + cx_24;
    float _S3021 = fy_24 * _S3019.y + cy_24;
    float2  _S3022 = make_float2 (_S3020, _S3021);
    float2  e0_5 = _S2995 - _S2968;
    float2  e1_5 = _S3022 - _S2995;
    float2  e2_1 = _S2968 - _S3022;
    float _S3023 = e0_5.x;
    float _S3024 = e1_5.y;
    float _S3025 = e0_5.y;
    float _S3026 = e1_5.x;
    float _S3027 = _S3023 * _S3024 - _S3025 * _S3026;
    float _S3028 = 1.0f - hardness_5.y;
    float _S3029 = -1.0f / _S3028;
    float _S3030 = _S3028 * _S3028;
    float _S3031 = s_primal_ctx_max_0(_S2966, _S2993);
    float _S3032 = s_primal_ctx_min_0(_S2966, _S2993);
    float _S3033 = s_primal_ctx_max_0(_S2967, _S2994);
    float _S3034 = s_primal_ctx_min_0(_S2967, _S2994);
    float3  _S3035 = make_float3 (_S2936, _S2937, _S2938) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3036 = transpose_0(R_24);
    float3  _S3037 = mean_24 - - s_primal_ctx_mul_1(_S3036, t_23);
    float _S3038 = _S3037.x;
    float _S3039 = _S3037.y;
    float _S3040 = _S3037.z;
    float _S3041 = _S3038 * _S3038 + _S3039 * _S3039 + _S3040 * _S3040;
    float _S3042 = s_primal_ctx_sqrt_0(_S3041);
    float x_50 = _S3038 / _S3042;
    float3  _S3043 = make_float3 (x_50);
    float _S3044 = _S3042 * _S3042;
    float y_23 = _S3039 / _S3042;
    float z_20 = _S3040 / _S3042;
    float3  _S3045 = make_float3 (z_20);
    float _S3046 = - y_23;
    float3  _S3047 = make_float3 (_S3046);
    float z2_46 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_50 * x_50 - y_23 * y_23;
    float _S3048 = 2.0f * x_50;
    float fS1_20 = _S3048 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_46 - 0.31539157032966614f;
    float3  _S3049 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_50;
    float3  _S3050 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3051 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3052 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3053 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_46 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3054 = 1.86588168144226074f * z2_46 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3054;
    float3  _S3055 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_50;
    float3  _S3056 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3057 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3058 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3059 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_50 * fC1_20 - y_23 * fS1_20);
    float3  _S3060 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_50 * fS1_20 + y_23 * fC1_20);
    float3  _S3061 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3046) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_50) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3062 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3063 = make_float3 (0.0f);
    float3  _S3064 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3065 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3066 = make_float3 (_S3065);
    float3  _S3067 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3065);
    float3  _S3068 = _S3064 + _S3067 + make_float3 (0.5f);
    float3  _S3069 = _S3064 - _S3067 + make_float3 (0.5f);
    float3  _S3070 = vert1_c_5 - vert0_c_5;
    float3  _S3071 = vert2_c_5 - vert0_c_5;
    float3  _S3072 = s_primal_ctx_cross_0(_S3070, _S3071);
    float3  _S3073 = normalize_0(_S3072);
    float3  _S3074 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3073, mean_c_20)))))) * v_normal_1;
    float3  _S3075 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3076;
    (&_S3076)->primal_0 = _S3073;
    (&_S3076)->differential_0 = _S3075;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3077;
    (&_S3077)->primal_0 = mean_c_20;
    (&_S3077)->differential_0 = _S3075;
    s_bwd_prop_dot_0(&_S3076, &_S3077, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3078 = _S3077;
    float3  _S3079 = _S3074 + _S3076.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3080;
    (&_S3080)->primal_0 = _S3072;
    (&_S3080)->differential_0 = _S3075;
    s_bwd_normalize_impl_0(&_S3080, _S3079);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3081;
    (&_S3081)->primal_0 = _S3070;
    (&_S3081)->differential_0 = _S3075;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3082;
    (&_S3082)->primal_0 = _S3071;
    (&_S3082)->differential_0 = _S3075;
    s_bwd_prop_cross_0(&_S3081, &_S3082, _S3080.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3083 = _S3081;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3084 = _S3082;
    float3  _S3085 = - _S3082.differential_0;
    float3  _S3086 = - _S3081.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3087;
    (&_S3087)->primal_0 = _S3069;
    (&_S3087)->differential_0 = _S3075;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3088;
    (&_S3088)->primal_0 = _S3063;
    (&_S3088)->differential_0 = _S3075;
    s_bwd_prop_max_0(&_S3087, &_S3088, (*v_rgb_7)[int(2)]);
    float3  _S3089 = - _S3087.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3090;
    (&_S3090)->primal_0 = _S3068;
    (&_S3090)->differential_0 = _S3075;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3091;
    (&_S3091)->primal_0 = _S3063;
    (&_S3091)->differential_0 = _S3075;
    s_bwd_prop_max_0(&_S3090, &_S3091, (*v_rgb_7)[int(1)]);
    float3  _S3092 = _S3066 * (_S3089 + _S3090.differential_0);
    float3  _S3093 = _S3087.differential_0 + _S3090.differential_0;
    float3  _S3094 = make_float3 (0.5f) * - _S3093;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3095;
    (&_S3095)->primal_0 = _S3062;
    (&_S3095)->differential_0 = _S3075;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3096;
    (&_S3096)->primal_0 = _S3063;
    (&_S3096)->differential_0 = _S3075;
    s_bwd_prop_max_0(&_S3095, &_S3096, (*v_rgb_7)[int(0)]);
    float3  _S3097 = _S3094 + _S3095.differential_0;
    float3  _S3098 = _S3093 + _S3095.differential_0;
    float3  _S3099 = _S3060 * _S3098;
    float3  _S3100 = (*sh_coeffs_20)[int(15)] * _S3098;
    float3  _S3101 = _S3058 * _S3098;
    float3  _S3102 = (*sh_coeffs_20)[int(14)] * _S3098;
    float3  _S3103 = _S3056 * _S3098;
    float3  _S3104 = (*sh_coeffs_20)[int(13)] * _S3098;
    float3  _S3105 = _S3055 * _S3098;
    float3  _S3106 = (*sh_coeffs_20)[int(12)] * _S3098;
    float3  _S3107 = _S3057 * _S3098;
    float3  _S3108 = (*sh_coeffs_20)[int(11)] * _S3098;
    float3  _S3109 = _S3059 * _S3098;
    float3  _S3110 = (*sh_coeffs_20)[int(10)] * _S3098;
    float3  _S3111 = _S3061 * _S3098;
    float3  _S3112 = (*sh_coeffs_20)[int(9)] * _S3098;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3112.x + _S3112.y + _S3112.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3100.x + _S3100.y + _S3100.z);
    float _S3113 = _S3110.x + _S3110.y + _S3110.z;
    float _S3114 = _S3102.x + _S3102.y + _S3102.z;
    float _S3115 = _S3108.x + _S3108.y + _S3108.z;
    float _S3116 = _S3104.x + _S3104.y + _S3104.z;
    float _S3117 = _S3106.x + _S3106.y + _S3106.z;
    float _S3118 = - s_diff_fC2_T_6;
    float3  _S3119 = _S3052 * _S3098;
    float3  _S3120 = (*sh_coeffs_20)[int(8)] * _S3098;
    float3  _S3121 = _S3050 * _S3098;
    float3  _S3122 = (*sh_coeffs_20)[int(7)] * _S3098;
    float3  _S3123 = _S3049 * _S3098;
    float3  _S3124 = (*sh_coeffs_20)[int(6)] * _S3098;
    float3  _S3125 = _S3051 * _S3098;
    float3  _S3126 = (*sh_coeffs_20)[int(5)] * _S3098;
    float3  _S3127 = _S3053 * _S3098;
    float3  _S3128 = (*sh_coeffs_20)[int(4)] * _S3098;
    float _S3129 = _S3126.x + _S3126.y + _S3126.z;
    float _S3130 = _S3122.x + _S3122.y + _S3122.z;
    float _S3131 = fTmp1B_20 * _S3113 + x_50 * s_diff_fS2_T_6 + y_23 * _S3118 + 0.54627424478530884f * (_S3128.x + _S3128.y + _S3128.z);
    float _S3132 = fTmp1B_20 * _S3114 + y_23 * s_diff_fS2_T_6 + x_50 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3120.x + _S3120.y + _S3120.z);
    float _S3133 = y_23 * - _S3132;
    float _S3134 = x_50 * _S3132;
    float _S3135 = z_20 * (1.86588168144226074f * (z_20 * _S3117) + -2.28522896766662598f * (y_23 * _S3115 + x_50 * _S3116) + 0.94617468118667603f * (_S3124.x + _S3124.y + _S3124.z));
    float3  _S3136 = make_float3 (0.48860251903533936f) * _S3098;
    float3  _S3137 = - _S3136;
    float3  _S3138 = _S3043 * _S3137;
    float3  _S3139 = (*sh_coeffs_20)[int(3)] * _S3137;
    float3  _S3140 = _S3045 * _S3136;
    float3  _S3141 = (*sh_coeffs_20)[int(2)] * _S3136;
    float3  _S3142 = _S3047 * _S3136;
    float3  _S3143 = (*sh_coeffs_20)[int(1)] * _S3136;
    float _S3144 = (_S3054 * _S3117 + 1.44530570507049561f * (fS1_20 * _S3113 + fC1_20 * _S3114) + -1.09254848957061768f * (y_23 * _S3129 + x_50 * _S3130) + _S3135 + _S3135 + _S3141.x + _S3141.y + _S3141.z) / _S3044;
    float _S3145 = _S3042 * _S3144;
    float _S3146 = (fTmp0C_20 * _S3115 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3118 + fTmp0B_20 * _S3129 + _S3048 * _S3131 + _S3133 + _S3133 + - (_S3143.x + _S3143.y + _S3143.z)) / _S3044;
    float _S3147 = _S3042 * _S3146;
    float _S3148 = (fTmp0C_20 * _S3116 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3130 + 2.0f * (y_23 * _S3131) + _S3134 + _S3134 + _S3139.x + _S3139.y + _S3139.z) / _S3044;
    float _S3149 = _S3042 * _S3148;
    float _S3150 = _S3040 * - _S3144 + _S3039 * - _S3146 + _S3038 * - _S3148;
    DiffPair_float_0 _S3151;
    (&_S3151)->primal_0 = _S3041;
    (&_S3151)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3151, _S3150);
    float _S3152 = _S3040 * _S3151.differential_0;
    float _S3153 = _S3039 * _S3151.differential_0;
    float _S3154 = _S3038 * _S3151.differential_0;
    float3  _S3155 = make_float3 (0.282094806432724f) * _S3098;
    float3  _S3156 = make_float3 (_S3149 + _S3154 + _S3154, _S3147 + _S3153 + _S3153, _S3145 + _S3152 + _S3152);
    float3  _S3157 = - - _S3156;
    Matrix<float, 3, 3>  _S3158 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3159;
    (&_S3159)->primal_0 = _S3036;
    (&_S3159)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3160;
    (&_S3160)->primal_0 = t_23;
    (&_S3160)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3159, &_S3160, _S3157);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3161 = _S3160;
    Matrix<float, 3, 3>  _S3162 = transpose_0(_S3159.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3163;
    (&_S3163)->primal_0 = _S3035;
    (&_S3163)->differential_0 = _S3075;
    s_bwd_prop_log_1(&_S3163, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3164 = _S3163;
    DiffPair_float_0 _S3165;
    (&_S3165)->primal_0 = _S3034;
    (&_S3165)->differential_0 = 0.0f;
    DiffPair_float_0 _S3166;
    (&_S3166)->primal_0 = _S3021;
    (&_S3166)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3165, &_S3166, 0.0f);
    DiffPair_float_0 _S3167;
    (&_S3167)->primal_0 = _S2967;
    (&_S3167)->differential_0 = 0.0f;
    DiffPair_float_0 _S3168;
    (&_S3168)->primal_0 = _S2994;
    (&_S3168)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3167, &_S3168, _S3165.differential_0);
    DiffPair_float_0 _S3169;
    (&_S3169)->primal_0 = _S3033;
    (&_S3169)->differential_0 = 0.0f;
    DiffPair_float_0 _S3170;
    (&_S3170)->primal_0 = _S3021;
    (&_S3170)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3169, &_S3170, 0.0f);
    DiffPair_float_0 _S3171;
    (&_S3171)->primal_0 = _S2967;
    (&_S3171)->differential_0 = 0.0f;
    DiffPair_float_0 _S3172;
    (&_S3172)->primal_0 = _S2994;
    (&_S3172)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3171, &_S3172, _S3169.differential_0);
    DiffPair_float_0 _S3173;
    (&_S3173)->primal_0 = _S3032;
    (&_S3173)->differential_0 = 0.0f;
    DiffPair_float_0 _S3174;
    (&_S3174)->primal_0 = _S3020;
    (&_S3174)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3173, &_S3174, 0.0f);
    DiffPair_float_0 _S3175;
    (&_S3175)->primal_0 = _S2966;
    (&_S3175)->differential_0 = 0.0f;
    DiffPair_float_0 _S3176;
    (&_S3176)->primal_0 = _S2993;
    (&_S3176)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3175, &_S3176, _S3173.differential_0);
    DiffPair_float_0 _S3177;
    (&_S3177)->primal_0 = _S3031;
    (&_S3177)->differential_0 = 0.0f;
    DiffPair_float_0 _S3178;
    (&_S3178)->primal_0 = _S3020;
    (&_S3178)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3177, &_S3178, 0.0f);
    DiffPair_float_0 _S3179;
    (&_S3179)->primal_0 = _S2966;
    (&_S3179)->differential_0 = 0.0f;
    DiffPair_float_0 _S3180;
    (&_S3180)->primal_0 = _S2993;
    (&_S3180)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3179, &_S3180, _S3177.differential_0);
    DiffPair_float_0 _S3181;
    (&_S3181)->primal_0 = _S3029;
    (&_S3181)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3181, 0.0f);
    float _S3182 = - (-1.0f * - (_S3181.differential_0 / _S3030));
    float2  _S3183 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3184;
    (&_S3184)->primal_0 = e2_1;
    (&_S3184)->differential_0 = _S3183;
    s_bwd_length_impl_1(&_S3184, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3185;
    (&_S3185)->primal_0 = e1_5;
    (&_S3185)->differential_0 = _S3183;
    s_bwd_length_impl_1(&_S3185, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3186;
    (&_S3186)->primal_0 = e0_5;
    (&_S3186)->differential_0 = _S3183;
    s_bwd_length_impl_1(&_S3186, -0.0f);
    DiffPair_float_0 _S3187;
    (&_S3187)->primal_0 = _S3027;
    (&_S3187)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3187, 0.0f);
    float _S3188 = - _S3187.differential_0;
    float2  _S3189 = _S3185.differential_0 + make_float2 (_S3025 * _S3188, _S3023 * _S3187.differential_0);
    float2  _S3190 = _S3186.differential_0 + make_float2 (_S3024 * _S3187.differential_0, _S3026 * _S3188);
    float2  _S3191 = v_uv2_1 + - _S3184.differential_0 + _S3189;
    float _S3192 = fy_24 * (_S3166.differential_0 + _S3170.differential_0 + _S3191.y);
    float _S3193 = fx_24 * (_S3174.differential_0 + _S3178.differential_0 + _S3191.x);
    float2  _S3194 = make_float2 (_S3193, _S3192);
    float2  _S3195 = _S3008 * _S3194;
    float2  _S3196 = _S3012 * _S3194;
    float _S3197 = r2_19 * _S3192;
    float _S3198 = p1_7 * _S3192;
    float _S3199 = _S3018 * _S3192;
    float _S3200 = v_19 * _S3192;
    float _S3201 = u_19 * _S3200;
    float _S3202 = r2_19 * _S3193;
    float _S3203 = p2_7 * _S3193;
    float _S3204 = _S3015 * _S3193;
    float _S3205 = v_19 * _S3193;
    float _S3206 = u_19 * _S3205;
    float _S3207 = _S3195.x + _S3195.y;
    float _S3208 = r2_19 * _S3207;
    float _S3209 = r2_19 * _S3208;
    float _S3210 = r2_19 * _S3209;
    float _S3211 = r2_19 * _S3210;
    float _S3212 = sy1_6 * _S3192 + _S3198 + sx1_6 * _S3193 + _S3203 + _S3011 * _S3207 + _S3010 * _S3208 + _S3009 * _S3209 + k4_6 * _S3210;
    float _S3213 = v_19 * _S3212;
    float _S3214 = u_19 * _S3212;
    float _S3215 = _S3017 * _S3198 + 2.0f * (v_19 * _S3198) + _S3016 * _S3192 + _S3013 * _S3193 + _S3213 + _S3213;
    float _S3216 = _S2961 * _S3200 + _S3014 * _S3203 + 2.0f * (u_19 * _S3203) + _S2957 * _S3205 + _S3214 + _S3214;
    float2  _S3217 = v_uv0_1 + _S3184.differential_0 + - _S3190;
    float2  _S3218 = v_out_hardness_1 + make_float2 (0.0f, _S3182);
    float _S3219 = _S3176.differential_0 + _S3180.differential_0;
    float2  _S3220 = v_uv1_1 + - _S3189 + _S3190;
    float3  _S3221 = _S3085 + _S3086;
    FixedArray<float3 , 2>  _S3222;
    _S3222[int(0)] = _S3075;
    _S3222[int(1)] = _S3075;
    _S3222[int(1)] = _S3092;
    _S3222[int(0)] = _S3097;
    float3  _S3223 = _S3222[int(0)];
    float3  _S3224 = _S3222[int(1)];
    FixedArray<float3 , 16>  _S3225;
    _S3225[int(0)] = _S3075;
    _S3225[int(1)] = _S3075;
    _S3225[int(2)] = _S3075;
    _S3225[int(3)] = _S3075;
    _S3225[int(4)] = _S3075;
    _S3225[int(5)] = _S3075;
    _S3225[int(6)] = _S3075;
    _S3225[int(7)] = _S3075;
    _S3225[int(8)] = _S3075;
    _S3225[int(9)] = _S3075;
    _S3225[int(10)] = _S3075;
    _S3225[int(11)] = _S3075;
    _S3225[int(12)] = _S3075;
    _S3225[int(13)] = _S3075;
    _S3225[int(14)] = _S3075;
    _S3225[int(15)] = _S3075;
    _S3225[int(7)] = _S3121;
    _S3225[int(0)] = _S3155;
    _S3225[int(1)] = _S3142;
    _S3225[int(2)] = _S3140;
    _S3225[int(3)] = _S3138;
    _S3225[int(4)] = _S3127;
    _S3225[int(5)] = _S3125;
    _S3225[int(6)] = _S3123;
    _S3225[int(15)] = _S3099;
    _S3225[int(8)] = _S3119;
    _S3225[int(9)] = _S3111;
    _S3225[int(10)] = _S3109;
    _S3225[int(11)] = _S3107;
    _S3225[int(12)] = _S3105;
    _S3225[int(13)] = _S3103;
    _S3225[int(14)] = _S3101;
    float3  _S3226 = _S3225[int(0)];
    float3  _S3227 = _S3225[int(1)];
    float3  _S3228 = _S3225[int(2)];
    float3  _S3229 = _S3225[int(3)];
    float3  _S3230 = _S3225[int(4)];
    float3  _S3231 = _S3225[int(5)];
    float3  _S3232 = _S3225[int(6)];
    float3  _S3233 = _S3225[int(7)];
    float3  _S3234 = _S3225[int(8)];
    float3  _S3235 = _S3225[int(9)];
    float3  _S3236 = _S3225[int(10)];
    float3  _S3237 = _S3225[int(11)];
    float3  _S3238 = _S3225[int(12)];
    float3  _S3239 = _S3225[int(13)];
    float3  _S3240 = _S3225[int(14)];
    float3  _S3241 = _S3225[int(15)];
    float _S3242 = _S3175.differential_0 + _S3179.differential_0;
    float _S3243 = _S3167.differential_0 + _S3171.differential_0;
    float _S3244 = _S3168.differential_0 + _S3172.differential_0;
    float2  _S3245 = _S3196 + make_float2 (_S3216, _S3215);
    float2  _S3246 = _S2996 * _S3245;
    float2  _S3247 = _S3007 * _S3245;
    float _S3248 = _S3246.x + _S3246.y;
    if(_S3000)
    {
        float _S3249 = _S3248 / _S3002;
        float _S3250 = _S3003 * - _S3249;
        float _S3251 = _S2999 * (0.3333333432674408f * - (_S2998 * _S3249));
        k_6 = _S3251 + _S3251;
        _S3001 = _S3250;
        _S3002 = 0.0f;
    }
    else
    {
        float _S3252 = _S3248 / _S3001;
        float _S3253 = _S2999 * - _S3252;
        k_6 = _S2997 * _S3252;
        _S3001 = 0.0f;
        _S3002 = _S3253;
    }
    DiffPair_float_0 _S3254;
    (&_S3254)->primal_0 = _S2997;
    (&_S3254)->differential_0 = 0.0f;
    DiffPair_float_0 _S3255;
    (&_S3255)->primal_0 = _S2998;
    (&_S3255)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3254, &_S3255, k_6);
    float _S3256 = _S3255.differential_0 + _S3001;
    float _S3257 = _S3254.differential_0 + _S3002;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3258;
    (&_S3258)->primal_0 = _S2996;
    (&_S3258)->differential_0 = _S3183;
    s_bwd_length_impl_1(&_S3258, _S3257);
    float2  _S3259 = _S3258.differential_0 + _S3247;
    float _S3260 = fy_24 * (_S3220.y + _S3244);
    float _S3261 = fx_24 * (_S3220.x + _S3219);
    float2  _S3262 = make_float2 (_S3261, _S3260);
    float2  _S3263 = _S2981 * _S3262;
    float _S3264 = p1_7 * _S3260;
    float _S3265 = v_18 * _S3260;
    float _S3266 = p2_7 * _S3261;
    float _S3267 = v_18 * _S3261;
    float _S3268 = _S3263.x + _S3263.y;
    float _S3269 = r2_18 * _S3268;
    float _S3270 = r2_18 * _S3269;
    float _S3271 = r2_18 * _S3270;
    float _S3272 = sy1_6 * _S3260 + _S3264 + sx1_6 * _S3261 + _S3266 + _S2984 * _S3268 + _S2983 * _S3269 + _S2982 * _S3270 + k4_6 * _S3271;
    float _S3273 = v_18 * _S3272;
    float _S3274 = u_18 * _S3272;
    float3  _S3275 = _S3084.differential_0 + make_float3 (_S3259.x, _S3259.y, _S3256);
    float2  _S3276 = _S2985 * _S3262 + make_float2 (_S2961 * _S3265 + _S2987 * _S3266 + 2.0f * (u_18 * _S3266) + _S2957 * _S3267 + _S3274 + _S3274, _S2990 * _S3264 + 2.0f * (v_18 * _S3264) + _S2989 * _S3260 + _S2986 * _S3261 + _S3273 + _S3273);
    float _S3277 = u_18 * _S3265 + _S3201;
    float _S3278 = u_18 * _S3267 + _S3206;
    float _S3279 = r2_18 * _S3260 + _S3197;
    float _S3280 = r2_18 * _S3271 + _S3211;
    float _S3281 = _S3271 + _S3210;
    float _S3282 = _S2991 * _S3260 + _S3199;
    float _S3283 = _S3270 + _S3209;
    float _S3284 = r2_18 * _S3261 + _S3202;
    float _S3285 = _S2988 * _S3261 + _S3204;
    float _S3286 = _S3269 + _S3208;
    float2  _S3287 = _S2969 * _S3276;
    float2  _S3288 = _S2980 * _S3276;
    float _S3289 = _S3287.x + _S3287.y;
    if(_S2973)
    {
        float _S3290 = _S3289 / _S2975;
        float _S3291 = _S2976 * - _S3290;
        float _S3292 = _S2972 * (0.3333333432674408f * - (_S2971 * _S3290));
        k_6 = _S3292 + _S3292;
        _S2974 = _S3291;
        _S2975 = 0.0f;
    }
    else
    {
        float _S3293 = _S3289 / _S2974;
        float _S3294 = _S2972 * - _S3293;
        k_6 = _S2970 * _S3293;
        _S2974 = 0.0f;
        _S2975 = _S3294;
    }
    DiffPair_float_0 _S3295;
    (&_S3295)->primal_0 = _S2970;
    (&_S3295)->differential_0 = 0.0f;
    DiffPair_float_0 _S3296;
    (&_S3296)->primal_0 = _S2971;
    (&_S3296)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3295, &_S3296, k_6);
    float _S3297 = _S3296.differential_0 + _S2974;
    float _S3298 = _S3295.differential_0 + _S2975;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3299;
    (&_S3299)->primal_0 = _S2969;
    (&_S3299)->differential_0 = _S3183;
    s_bwd_length_impl_1(&_S3299, _S3298);
    float2  _S3300 = _S3299.differential_0 + _S3288;
    float _S3301 = fy_24 * (_S3217.y + _S3243);
    float _S3302 = fx_24 * (_S3217.x + _S3242);
    float2  _S3303 = make_float2 (_S3302, _S3301);
    float2  _S3304 = _S2952 * _S3303;
    float _S3305 = p1_7 * _S3301;
    float _S3306 = v_17 * _S3301;
    float _S3307 = p2_7 * _S3302;
    float _S3308 = v_17 * _S3302;
    float _S3309 = _S3304.x + _S3304.y;
    float _S3310 = r2_17 * _S3309;
    float _S3311 = r2_17 * _S3310;
    float _S3312 = r2_17 * _S3311;
    float _S3313 = sy1_6 * _S3301 + _S3305 + sx1_6 * _S3302 + _S3307 + _S2955 * _S3309 + _S2954 * _S3310 + _S2953 * _S3311 + k4_6 * _S3312;
    float _S3314 = v_17 * _S3313;
    float _S3315 = u_17 * _S3313;
    float2  _S3316 = make_float2 (r2_17 * _S3302 + _S3284, r2_17 * _S3301 + _S3279);
    float2  _S3317 = make_float2 (_S2964 * _S3301 + 2.0f * (u_17 * _S3308 + _S3278) + _S3282, 2.0f * (u_17 * _S3306 + _S3277) + _S2960 * _S3302 + _S3285);
    float4  _S3318 = make_float4 (_S3310 + _S3286, _S3311 + _S3283, _S3312 + _S3281, r2_17 * _S3312 + _S3280);
    float3  _S3319 = _S3083.differential_0 + make_float3 (_S3300.x, _S3300.y, _S3297);
    float2  _S3320 = _S2956 * _S3303 + make_float2 (_S2961 * _S3306 + _S2959 * _S3307 + 2.0f * (u_17 * _S3307) + _S2957 * _S3308 + _S3315 + _S3315, _S2963 * _S3305 + 2.0f * (v_17 * _S3305) + _S2962 * _S3301 + _S2958 * _S3302 + _S3314 + _S3314);
    CameraDistortion_0 _S3321 = CameraDistortion_x24_syn_dzero_0();
    (&_S3321)->thin_prism_coeffs_0 = _S3316;
    (&_S3321)->tangential_coeffs_0 = _S3317;
    (&_S3321)->radial_coeffs_0 = _S3318;
    float2  _S3322 = _S2940 * _S3320;
    float2  _S3323 = _S2951 * _S3320;
    float _S3324 = _S3322.x + _S3322.y;
    if(_S2944)
    {
        float _S3325 = _S3324 / _S2946;
        float _S3326 = _S2947 * - _S3325;
        float _S3327 = _S2943 * (0.3333333432674408f * - (_S2942 * _S3325));
        k_6 = _S3327 + _S3327;
        _S2945 = _S3326;
        _S2946 = 0.0f;
    }
    else
    {
        float _S3328 = _S3324 / _S2945;
        float _S3329 = _S2943 * - _S3328;
        k_6 = _S2941 * _S3328;
        _S2945 = 0.0f;
        _S2946 = _S3329;
    }
    DiffPair_float_0 _S3330;
    (&_S3330)->primal_0 = _S2941;
    (&_S3330)->differential_0 = 0.0f;
    DiffPair_float_0 _S3331;
    (&_S3331)->primal_0 = _S2942;
    (&_S3331)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3330, &_S3331, k_6);
    float _S3332 = _S3331.differential_0 + _S2945;
    float _S3333 = _S3330.differential_0 + _S2946;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3334;
    (&_S3334)->primal_0 = _S2940;
    (&_S3334)->differential_0 = _S3183;
    s_bwd_length_impl_1(&_S3334, _S3333);
    float2  _S3335 = _S3334.differential_0 + _S3323;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3336;
    (&_S3336)->primal_0 = vert2_c_5;
    (&_S3336)->differential_0 = _S3075;
    s_bwd_length_impl_0(&_S3336, _S3164.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3337;
    (&_S3337)->primal_0 = vert1_c_5;
    (&_S3337)->differential_0 = _S3075;
    s_bwd_length_impl_0(&_S3337, _S3164.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3338;
    (&_S3338)->primal_0 = vert0_c_5;
    (&_S3338)->differential_0 = _S3075;
    s_bwd_length_impl_0(&_S3338, _S3164.differential_0.x);
    float3  _S3339 = _S3336.differential_0 + _S3275;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3340;
    (&_S3340)->primal_0 = R_24;
    (&_S3340)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3341;
    (&_S3341)->primal_0 = vert2_1;
    (&_S3341)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3340, &_S3341, _S3339);
    float3  _S3342 = _S3337.differential_0 + _S3319;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3343;
    (&_S3343)->primal_0 = R_24;
    (&_S3343)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3344;
    (&_S3344)->primal_0 = vert1_1;
    (&_S3344)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3343, &_S3344, _S3342);
    float3  _S3345 = _S3338.differential_0 + _S3221 + make_float3 (_S3335.x, _S3335.y, _S3332);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3346;
    (&_S3346)->primal_0 = R_24;
    (&_S3346)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3347;
    (&_S3347)->primal_0 = vert0_1;
    (&_S3347)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3346, &_S3347, _S3345);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3348;
    (&_S3348)->primal_0 = _S2930;
    (&_S3348)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3349;
    (&_S3349)->primal_0 = _S2935;
    (&_S3349)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3348, &_S3349, _S3341.differential_0);
    float _S3350 = - _S3349.differential_0.y;
    float _S3351 = _S2934 * _S3349.differential_0.x;
    float _S3352 = - (_S2926 * _S3349.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3353;
    (&_S3353)->primal_0 = _S2930;
    (&_S3353)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3354;
    (&_S3354)->primal_0 = _S2933;
    (&_S3354)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3353, &_S3354, _S3344.differential_0);
    float _S3355 = _S2926 * _S3354.differential_0.x;
    float _S3356 = _S2932 * _S3354.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3357;
    (&_S3357)->primal_0 = _S2930;
    (&_S3357)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3358;
    (&_S3358)->primal_0 = _S2931;
    (&_S3358)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3357, &_S3358, _S3347.differential_0);
    Matrix<float, 3, 3>  _S3359 = transpose_0(_S3348.differential_0 + _S3353.differential_0 + _S3357.differential_0);
    float _S3360 = 2.0f * - _S3359.rows[int(2)].z;
    float _S3361 = 2.0f * _S3359.rows[int(2)].y;
    float _S3362 = 2.0f * _S3359.rows[int(2)].x;
    float _S3363 = 2.0f * _S3359.rows[int(1)].z;
    float _S3364 = 2.0f * - _S3359.rows[int(1)].y;
    float _S3365 = 2.0f * _S3359.rows[int(1)].x;
    float _S3366 = 2.0f * _S3359.rows[int(0)].z;
    float _S3367 = 2.0f * _S3359.rows[int(0)].y;
    float _S3368 = 2.0f * - _S3359.rows[int(0)].x;
    float _S3369 = - _S3365 + _S3367;
    float _S3370 = _S3362 + - _S3366;
    float _S3371 = - _S3361 + _S3363;
    float _S3372 = _S3361 + _S3363;
    float _S3373 = _S3362 + _S3366;
    float _S3374 = _S3365 + _S3367;
    float _S3375 = quat_25.w * (_S3364 + _S3368);
    float _S3376 = quat_25.z * (_S3360 + _S3368);
    float _S3377 = quat_25.y * (_S3360 + _S3364);
    float _S3378 = quat_25.x * _S3369 + quat_25.z * _S3372 + quat_25.y * _S3373 + _S3375 + _S3375;
    float _S3379 = quat_25.x * _S3370 + quat_25.w * _S3372 + quat_25.y * _S3374 + _S3376 + _S3376;
    float _S3380 = quat_25.x * _S3371 + quat_25.w * _S3373 + quat_25.z * _S3374 + _S3377 + _S3377;
    float _S3381 = quat_25.w * _S3369 + quat_25.z * _S3370 + quat_25.y * _S3371;
    float _S3382 = _S3352 + _S3355;
    float _S3383 = 0.5f * - _S3382;
    float _S3384 = _S3350 + _S3354.differential_0.y;
    DiffPair_float_0 _S3385;
    (&_S3385)->primal_0 = _S2927;
    (&_S3385)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3385, _S3384);
    float _S3386 = _S3383 + _S3385.differential_0;
    float _S3387 = _S3351 + _S3356 + _S3358.differential_0.x;
    DiffPair_float_0 _S3388;
    (&_S3388)->primal_0 = _S2925;
    (&_S3388)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3388, _S3387);
    float _S3389 = _S3383 + _S3388.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3390;
    (&_S3390)->primal_0 = mean_c_20;
    (&_S3390)->differential_0 = _S3075;
    s_bwd_length_impl_0(&_S3390, 0.0f);
    float3  _S3391 = _S3390.differential_0 + _S3078.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3392;
    (&_S3392)->primal_0 = R_24;
    (&_S3392)->differential_0 = _S3158;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3393;
    (&_S3393)->primal_0 = mean_24;
    (&_S3393)->differential_0 = _S3075;
    s_bwd_prop_mul_1(&_S3392, &_S3393, _S3391);
    float3  _S3394 = _S3339 + _S3342 + _S3345 + _S3391 + _S3161.differential_0;
    Matrix<float, 3, 3>  _S3395 = _S3340.differential_0 + _S3343.differential_0 + _S3346.differential_0 + _S3392.differential_0 + _S3162;
    float3  _S3396 = make_float3 (_S3389, _S3386, _S3382);
    float4  _S3397 = make_float4 (0.0f);
    *&((&_S3397)->w) = _S3378;
    *&((&_S3397)->z) = _S3379;
    *&((&_S3397)->y) = _S3380;
    *&((&_S3397)->x) = _S3381;
    float4  _S3398 = _S3397;
    float3  _S3399 = _S3341.differential_0 + _S3344.differential_0 + _S3347.differential_0 + _S3393.differential_0 + _S3156;
    *v_mean_8 = _S3399;
    *v_quat_7 = _S3398;
    *v_scale_7 = _S3396;
    *v_hardness_1 = _S3218;
    (*v_sh_coeffs_6)[int(0)] = _S3226;
    (*v_sh_coeffs_6)[int(1)] = _S3227;
    (*v_sh_coeffs_6)[int(2)] = _S3228;
    (*v_sh_coeffs_6)[int(3)] = _S3229;
    (*v_sh_coeffs_6)[int(4)] = _S3230;
    (*v_sh_coeffs_6)[int(5)] = _S3231;
    (*v_sh_coeffs_6)[int(6)] = _S3232;
    (*v_sh_coeffs_6)[int(7)] = _S3233;
    (*v_sh_coeffs_6)[int(8)] = _S3234;
    (*v_sh_coeffs_6)[int(9)] = _S3235;
    (*v_sh_coeffs_6)[int(10)] = _S3236;
    (*v_sh_coeffs_6)[int(11)] = _S3237;
    (*v_sh_coeffs_6)[int(12)] = _S3238;
    (*v_sh_coeffs_6)[int(13)] = _S3239;
    (*v_sh_coeffs_6)[int(14)] = _S3240;
    (*v_sh_coeffs_6)[int(15)] = _S3241;
    (*v_ch_coeffs_1)[int(0)] = _S3223;
    (*v_ch_coeffs_1)[int(1)] = _S3224;
    *v_R_7 = _S3395;
    *v_t_7 = _S3394;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_19)
{
    DiffPair_float_0 _S3400 = *dpx_16;
    bool _S3401;
    if(((*dpx_16).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S3401 = ((*dpx_16).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S3401 = false;
    }
    float _S3402;
    if(_S3401)
    {
        _S3402 = dOut_19;
    }
    else
    {
        _S3402 = 0.0f;
    }
    dpx_16->primal_0 = _S3400.primal_0;
    dpx_16->differential_0 = _S3402;
    DiffPair_float_0 _S3403 = *dpMin_0;
    if((_S3400.primal_0) < ((*dpMin_0).primal_0))
    {
        _S3402 = dOut_19;
    }
    else
    {
        _S3402 = 0.0f;
    }
    dpMin_0->primal_0 = _S3403.primal_0;
    dpMin_0->differential_0 = _S3402;
    DiffPair_float_0 _S3404 = *dpMax_0;
    if(((*dpx_16).primal_0) > ((*dpMax_0).primal_0))
    {
        _S3402 = dOut_19;
    }
    else
    {
        _S3402 = 0.0f;
    }
    dpMax_0->primal_0 = _S3404.primal_0;
    dpMax_0->differential_0 = _S3402;
    return;
}

inline __device__ float clamp_0(float x_51, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_51), (minBound_0)))), (maxBound_0)));
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
        DiffPair_float_0 _S3405 = *dpx_17;
        float _S3406 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S3406;
        float _S3407 = val_0 * (F32_log((_S3405.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3407;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3408 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3408))));
    float2  _S3409 = p_0 - v0_0;
    float2  _S3410 = normalize_1(e0_6);
    float2  _S3411 = p_0 - v1_0;
    float2  _S3412 = normalize_1(e1_6);
    float2  _S3413 = p_0 - v2_0;
    float2  _S3414 = normalize_1(e2_2);
    float _S3415 = hardness_6.x;
    float _S3416 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3409.x * _S3410.y - _S3409.y * _S3410.x)), (se_0 * (_S3411.x * _S3412.y - _S3411.y * _S3412.x))))), (se_0 * (_S3413.x * _S3414.y - _S3413.y * _S3414.x)))) / ((F32_abs((_S3408))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S3416))));
    float _S3417;
    if(a_1 <= 0.0f)
    {
        _S3417 = 0.0f;
    }
    else
    {
        _S3417 = (F32_min(((F32_pow((a_1), (_S3416)))), (0.99900001287460327f)));
    }
    return _S3415 * _S3417;
}

inline __device__ float s_primal_ctx_abs_0(float _S3418)
{
    return (F32_abs((_S3418)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S3419, float _S3420, float _S3421)
{
    return clamp_0(_S3419, _S3420, _S3421);
}

inline __device__ float s_primal_ctx_exp2_0(float _S3422)
{
    return (F32_exp2((_S3422)));
}

inline __device__ float s_primal_ctx_pow_0(float _S3423, float _S3424)
{
    return (F32_pow((_S3423), (_S3424)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S3425, DiffPair_float_0 * _S3426, float _S3427)
{
    _d_pow_0(_S3425, _S3426, _S3427);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S3428, DiffPair_float_0 * _S3429, DiffPair_float_0 * _S3430, float _S3431)
{
    _d_clamp_0(_S3428, _S3429, _S3430, _S3431);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_9)
{
    float _S3432 = length_0((*dpx_18).primal_0);
    float2  _S3433 = (*dpx_18).primal_0 * _s_dOut_9;
    float2  _S3434 = make_float2 (1.0f / _S3432) * _s_dOut_9;
    float _S3435 = - ((_S3433.x + _S3433.y) / (_S3432 * _S3432));
    float2  _S3436 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3437;
    (&_S3437)->primal_0 = (*dpx_18).primal_0;
    (&_S3437)->differential_0 = _S3436;
    s_bwd_length_impl_1(&_S3437, _S3435);
    float2  _S3438 = _S3434 + _S3437.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S3438;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3439, float2  _S3440)
{
    s_bwd_prop_normalize_impl_1(_S3439, _S3440);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_10)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S3441 = e0_7.x;
    float _S3442 = e1_7.y;
    float _S3443 = e0_7.y;
    float _S3444 = e1_7.x;
    float _S3445 = _S3441 * _S3442 - _S3443 * _S3444;
    float se_1 = float((F32_sign((_S3445))));
    float2  _S3446 = p_1 - (*dpv0_0).primal_0;
    float2  _S3447 = normalize_1(e0_7);
    float _S3448 = _S3446.x;
    float _S3449 = _S3447.y;
    float _S3450 = _S3446.y;
    float _S3451 = _S3447.x;
    float de0_0 = se_1 * (_S3448 * _S3449 - _S3450 * _S3451);
    float2  _S3452 = p_1 - (*dpv1_0).primal_0;
    float2  _S3453 = normalize_1(e1_7);
    float _S3454 = _S3452.x;
    float _S3455 = _S3453.y;
    float _S3456 = _S3452.y;
    float _S3457 = _S3453.x;
    float de1_0 = se_1 * (_S3454 * _S3455 - _S3456 * _S3457);
    float2  _S3458 = p_1 - (*dpv2_0).primal_0;
    float2  _S3459 = normalize_1(e2_3);
    float _S3460 = _S3458.x;
    float _S3461 = _S3459.y;
    float _S3462 = _S3458.y;
    float _S3463 = _S3459.x;
    float de2_0 = se_1 * (_S3460 * _S3461 - _S3462 * _S3463);
    float _S3464 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S3465 = s_primal_ctx_max_0(_S3464, de2_0);
    float _S3466 = s_primal_ctx_abs_0(_S3445);
    float _S3467 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S3466 / _S3467;
    float _S3468 = _S3467 * _S3467;
    float _S3469 = (*dphardness_0).primal_0.x;
    float _S3470 = (*dphardness_0).primal_0.y;
    float _S3471 = dmax_0 * dmax_0;
    float _S3472 = 1.0f + _S3465 / dmax_0;
    float _S3473 = 1.0f - s_primal_ctx_clamp_0(_S3470, 0.00499999988824129f, 0.98000001907348633f);
    float _S3474 = -1.0f / _S3473;
    float _S3475 = _S3473 * _S3473;
    float _S3476 = 1.0f - s_primal_ctx_exp2_0(_S3474);
    float a_2 = 1.0f - _S3472 * _S3476;
    bool _S3477 = a_2 <= 0.0f;
    float _S3478;
    float _S3479;
    if(_S3477)
    {
        _S3478 = 0.0f;
        _S3479 = 0.0f;
    }
    else
    {
        float _S3480 = s_primal_ctx_pow_0(a_2, _S3473);
        _S3478 = s_primal_ctx_min_0(_S3480, 0.99900001287460327f);
        _S3479 = _S3480;
    }
    float _S3481 = _S3469 * _s_dOut_10;
    float _S3482 = _S3478 * _s_dOut_10;
    if(_S3477)
    {
        _S3478 = 0.0f;
        _S3479 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S3483;
        (&_S3483)->primal_0 = _S3479;
        (&_S3483)->differential_0 = 0.0f;
        DiffPair_float_0 _S3484;
        (&_S3484)->primal_0 = 0.99900001287460327f;
        (&_S3484)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S3483, &_S3484, _S3481);
        DiffPair_float_0 _S3485;
        (&_S3485)->primal_0 = a_2;
        (&_S3485)->differential_0 = 0.0f;
        DiffPair_float_0 _S3486;
        (&_S3486)->primal_0 = _S3473;
        (&_S3486)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S3485, &_S3486, _S3483.differential_0);
        _S3478 = _S3485.differential_0;
        _S3479 = _S3486.differential_0;
    }
    float _S3487 = - _S3478;
    float _S3488 = _S3476 * _S3487;
    float _S3489 = - (_S3472 * _S3487);
    DiffPair_float_0 _S3490;
    (&_S3490)->primal_0 = _S3474;
    (&_S3490)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3490, _S3489);
    float _S3491 = - (-1.0f * - (_S3490.differential_0 / _S3475) + _S3479);
    float _S3492 = _S3488 / _S3471;
    float s_diff_dmax_T_0 = _S3465 * - _S3492;
    float _S3493 = dmax_0 * _S3492;
    DiffPair_float_0 _S3494;
    (&_S3494)->primal_0 = _S3470;
    (&_S3494)->differential_0 = 0.0f;
    DiffPair_float_0 _S3495;
    (&_S3495)->primal_0 = 0.00499999988824129f;
    (&_S3495)->differential_0 = 0.0f;
    DiffPair_float_0 _S3496;
    (&_S3496)->primal_0 = 0.98000001907348633f;
    (&_S3496)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3494, &_S3495, &_S3496, _S3491);
    float _S3497 = s_diff_dmax_T_0 / _S3468;
    float _S3498 = _S3466 * - _S3497;
    float _S3499 = _S3467 * _S3497;
    float2  _S3500 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3501;
    (&_S3501)->primal_0 = e2_3;
    (&_S3501)->differential_0 = _S3500;
    s_bwd_length_impl_1(&_S3501, _S3498);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3502;
    (&_S3502)->primal_0 = e1_7;
    (&_S3502)->differential_0 = _S3500;
    s_bwd_length_impl_1(&_S3502, _S3498);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3503;
    (&_S3503)->primal_0 = e0_7;
    (&_S3503)->differential_0 = _S3500;
    s_bwd_length_impl_1(&_S3503, _S3498);
    DiffPair_float_0 _S3504;
    (&_S3504)->primal_0 = _S3445;
    (&_S3504)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3504, _S3499);
    DiffPair_float_0 _S3505;
    (&_S3505)->primal_0 = _S3464;
    (&_S3505)->differential_0 = 0.0f;
    DiffPair_float_0 _S3506;
    (&_S3506)->primal_0 = de2_0;
    (&_S3506)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3505, &_S3506, _S3493);
    DiffPair_float_0 _S3507;
    (&_S3507)->primal_0 = de0_0;
    (&_S3507)->differential_0 = 0.0f;
    DiffPair_float_0 _S3508;
    (&_S3508)->primal_0 = de1_0;
    (&_S3508)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3507, &_S3508, _S3505.differential_0);
    float _S3509 = se_1 * _S3506.differential_0;
    float _S3510 = - _S3509;
    float _S3511 = _S3463 * _S3510;
    float _S3512 = _S3461 * _S3509;
    float2  _S3513 = make_float2 (_S3462 * _S3510, _S3460 * _S3509);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3514;
    (&_S3514)->primal_0 = e2_3;
    (&_S3514)->differential_0 = _S3500;
    s_bwd_normalize_impl_1(&_S3514, _S3513);
    float2  _S3515 = - make_float2 (_S3512, _S3511);
    float _S3516 = se_1 * _S3508.differential_0;
    float _S3517 = - _S3516;
    float _S3518 = _S3457 * _S3517;
    float _S3519 = _S3455 * _S3516;
    float2  _S3520 = make_float2 (_S3456 * _S3517, _S3454 * _S3516);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3521;
    (&_S3521)->primal_0 = e1_7;
    (&_S3521)->differential_0 = _S3500;
    s_bwd_normalize_impl_1(&_S3521, _S3520);
    float2  _S3522 = - make_float2 (_S3519, _S3518);
    float _S3523 = se_1 * _S3507.differential_0;
    float _S3524 = - _S3523;
    float _S3525 = _S3451 * _S3524;
    float _S3526 = _S3449 * _S3523;
    float2  _S3527 = make_float2 (_S3450 * _S3524, _S3448 * _S3523);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3528;
    (&_S3528)->primal_0 = e0_7;
    (&_S3528)->differential_0 = _S3500;
    s_bwd_normalize_impl_1(&_S3528, _S3527);
    float2  _S3529 = - make_float2 (_S3526, _S3525);
    float _S3530 = - _S3504.differential_0;
    float2  _S3531 = _S3501.differential_0 + _S3514.differential_0;
    float2  _S3532 = - _S3531;
    float2  _S3533 = _S3502.differential_0 + _S3521.differential_0 + make_float2 (_S3443 * _S3530, _S3441 * _S3504.differential_0);
    float2  _S3534 = - _S3533;
    float2  _S3535 = _S3503.differential_0 + _S3528.differential_0 + make_float2 (_S3442 * _S3504.differential_0, _S3444 * _S3530);
    float2  _S3536 = - _S3535;
    float2  _S3537 = make_float2 (_S3482, _S3494.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S3537;
    float2  _S3538 = _S3515 + _S3532 + _S3533;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S3538;
    float2  _S3539 = _S3522 + _S3534 + _S3535;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S3539;
    float2  _S3540 = _S3529 + _S3531 + _S3536;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S3540;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3541, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3542, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3543, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3544, float2  _S3545, float _S3546)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S3541, _S3542, _S3543, _S3544, _S3545, _S3546);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S3547 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S3547;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S3547;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S3547;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S3547;
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
    float2  _S3548 = p_3 - v0_2;
    float2  _S3549 = p_3 - v1_2;
    float2  _S3550 = p_3 - v2_2;
    float _S3551 = e0_8.x;
    float _S3552 = e1_8.y;
    float _S3553 = e0_8.y;
    float _S3554 = e1_8.x;
    float _S3555 = _S3551 * _S3552 - _S3553 * _S3554;
    float se_2 = float((F32_sign((_S3555))));
    float _S3556 = hardness_8.x;
    float _S3557 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S3548.x * _S3553 - _S3548.y * _S3551)), (se_2 * (_S3549.x * _S3552 - _S3549.y * _S3554))))), (se_2 * (_S3550.x * e2_4.y - _S3550.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S3548 - e0_8 * make_float2 (clamp_0(dot_1(_S3548, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S3549 - e1_8 * make_float2 (clamp_0(dot_1(_S3549, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S3550 - e2_4 * make_float2 (clamp_0(dot_1(_S3550, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S3555))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S3557))));
    float _S3558;
    if(a_3 <= 0.0f)
    {
        _S3558 = 0.0f;
    }
    else
    {
        _S3558 = (F32_min(((F32_pow((a_3), (_S3557)))), (0.99900001287460327f)));
    }
    return _S3556 * _S3558;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S3559, float2  _S3560)
{
    return dot_1(_S3559, _S3560);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3561, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3562, float _S3563)
{
    _d_dot_1(_S3561, _S3562, _S3563);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_11)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S3564 = p_4 - (*dpv0_1).primal_0;
    float _S3565 = s_primal_ctx_dot_1(_S3564, e0_9);
    float _S3566 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S3567 = _S3565 / _S3566;
    float _S3568 = _S3566 * _S3566;
    float _S3569 = s_primal_ctx_clamp_0(_S3567, 0.0f, 1.0f);
    float2  _S3570 = make_float2 (_S3569);
    float2  _S3571 = _S3564 - e0_9 * make_float2 (_S3569);
    float _S3572 = length_0(_S3571);
    float2  _S3573 = p_4 - (*dpv1_1).primal_0;
    float _S3574 = s_primal_ctx_dot_1(_S3573, e1_9);
    float _S3575 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S3576 = _S3574 / _S3575;
    float _S3577 = _S3575 * _S3575;
    float _S3578 = s_primal_ctx_clamp_0(_S3576, 0.0f, 1.0f);
    float2  _S3579 = make_float2 (_S3578);
    float2  _S3580 = _S3573 - e1_9 * make_float2 (_S3578);
    float _S3581 = length_0(_S3580);
    float2  _S3582 = p_4 - (*dpv2_1).primal_0;
    float _S3583 = s_primal_ctx_dot_1(_S3582, e2_5);
    float _S3584 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S3585 = _S3583 / _S3584;
    float _S3586 = _S3584 * _S3584;
    float _S3587 = s_primal_ctx_clamp_0(_S3585, 0.0f, 1.0f);
    float2  _S3588 = make_float2 (_S3587);
    float2  _S3589 = _S3582 - e2_5 * make_float2 (_S3587);
    float _S3590 = length_0(_S3589);
    float _S3591 = e0_9.x;
    float _S3592 = e1_9.y;
    float _S3593 = e0_9.y;
    float _S3594 = e1_9.x;
    float _S3595 = _S3591 * _S3592 - _S3593 * _S3594;
    float se_3 = float((F32_sign((_S3595))));
    float _S3596 = _S3564.x;
    float _S3597 = _S3564.y;
    float s0_0 = se_3 * (_S3596 * _S3593 - _S3597 * _S3591);
    float _S3598 = _S3573.x;
    float _S3599 = _S3573.y;
    float s1_0 = se_3 * (_S3598 * _S3592 - _S3599 * _S3594);
    float _S3600 = _S3582.x;
    float _S3601 = e2_5.y;
    float _S3602 = _S3582.y;
    float _S3603 = e2_5.x;
    float s2_0 = se_3 * (_S3600 * _S3601 - _S3602 * _S3603);
    float _S3604 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S3604, s2_0)))));
    float _S3605 = s_primal_ctx_min_0(_S3572, _S3581);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S3605, _S3590);
    float _S3606 = s_primal_ctx_abs_0(_S3595);
    float _S3607 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S3606 / _S3607;
    float _S3608 = _S3607 * _S3607;
    float _S3609 = (*dphardness_1).primal_0.x;
    float _S3610 = (*dphardness_1).primal_0.y;
    float _S3611 = dmax_1 * dmax_1;
    float _S3612 = 1.0f + dv_0 / dmax_1;
    float _S3613 = 1.0f - s_primal_ctx_clamp_0(_S3610, 0.00499999988824129f, 0.98000001907348633f);
    float _S3614 = -1.0f / _S3613;
    float _S3615 = _S3613 * _S3613;
    float _S3616 = 1.0f - s_primal_ctx_exp2_0(_S3614);
    float a_4 = 1.0f - _S3612 * _S3616;
    bool _S3617 = a_4 <= 0.0f;
    float _S3618;
    float _S3619;
    if(_S3617)
    {
        _S3618 = 0.0f;
        _S3619 = 0.0f;
    }
    else
    {
        float _S3620 = s_primal_ctx_pow_0(a_4, _S3613);
        _S3618 = s_primal_ctx_min_0(_S3620, 0.99900001287460327f);
        _S3619 = _S3620;
    }
    float _S3621 = _S3609 * _s_dOut_11;
    float _S3622 = _S3618 * _s_dOut_11;
    if(_S3617)
    {
        _S3618 = 0.0f;
        _S3619 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S3623;
        (&_S3623)->primal_0 = _S3619;
        (&_S3623)->differential_0 = 0.0f;
        DiffPair_float_0 _S3624;
        (&_S3624)->primal_0 = 0.99900001287460327f;
        (&_S3624)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S3623, &_S3624, _S3621);
        DiffPair_float_0 _S3625;
        (&_S3625)->primal_0 = a_4;
        (&_S3625)->differential_0 = 0.0f;
        DiffPair_float_0 _S3626;
        (&_S3626)->primal_0 = _S3613;
        (&_S3626)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S3625, &_S3626, _S3623.differential_0);
        _S3618 = _S3625.differential_0;
        _S3619 = _S3626.differential_0;
    }
    float _S3627 = - _S3618;
    float _S3628 = _S3616 * _S3627;
    float _S3629 = - (_S3612 * _S3627);
    DiffPair_float_0 _S3630;
    (&_S3630)->primal_0 = _S3614;
    (&_S3630)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3630, _S3629);
    float _S3631 = - (-1.0f * - (_S3630.differential_0 / _S3615) + _S3619);
    float _S3632 = _S3628 / _S3611;
    float s_diff_dmax_T_1 = dv_0 * - _S3632;
    float s_diff_dv_T_0 = dmax_1 * _S3632;
    DiffPair_float_0 _S3633;
    (&_S3633)->primal_0 = _S3610;
    (&_S3633)->differential_0 = 0.0f;
    DiffPair_float_0 _S3634;
    (&_S3634)->primal_0 = 0.00499999988824129f;
    (&_S3634)->differential_0 = 0.0f;
    DiffPair_float_0 _S3635;
    (&_S3635)->primal_0 = 0.98000001907348633f;
    (&_S3635)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3633, &_S3634, &_S3635, _S3631);
    float _S3636 = s_diff_dmax_T_1 / _S3608;
    float _S3637 = _S3606 * - _S3636;
    float _S3638 = _S3607 * _S3636;
    float2  _S3639 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3640;
    (&_S3640)->primal_0 = e2_5;
    (&_S3640)->differential_0 = _S3639;
    s_bwd_length_impl_1(&_S3640, _S3637);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3641;
    (&_S3641)->primal_0 = e1_9;
    (&_S3641)->differential_0 = _S3639;
    s_bwd_length_impl_1(&_S3641, _S3637);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3642;
    (&_S3642)->primal_0 = e0_9;
    (&_S3642)->differential_0 = _S3639;
    s_bwd_length_impl_1(&_S3642, _S3637);
    DiffPair_float_0 _S3643;
    (&_S3643)->primal_0 = _S3595;
    (&_S3643)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3643, _S3638);
    float _S3644 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S3645;
    (&_S3645)->primal_0 = _S3605;
    (&_S3645)->differential_0 = 0.0f;
    DiffPair_float_0 _S3646;
    (&_S3646)->primal_0 = _S3590;
    (&_S3646)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3645, &_S3646, _S3644);
    DiffPair_float_0 _S3647;
    (&_S3647)->primal_0 = _S3572;
    (&_S3647)->differential_0 = 0.0f;
    DiffPair_float_0 _S3648;
    (&_S3648)->primal_0 = _S3581;
    (&_S3648)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3647, &_S3648, _S3645.differential_0);
    DiffPair_float_0 _S3649;
    (&_S3649)->primal_0 = _S3604;
    (&_S3649)->differential_0 = 0.0f;
    DiffPair_float_0 _S3650;
    (&_S3650)->primal_0 = s2_0;
    (&_S3650)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3649, &_S3650, 0.0f);
    DiffPair_float_0 _S3651;
    (&_S3651)->primal_0 = s0_0;
    (&_S3651)->differential_0 = 0.0f;
    DiffPair_float_0 _S3652;
    (&_S3652)->primal_0 = s1_0;
    (&_S3652)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3651, &_S3652, _S3649.differential_0);
    float _S3653 = se_3 * _S3650.differential_0;
    float _S3654 = - _S3653;
    float _S3655 = _S3602 * _S3654;
    float _S3656 = _S3603 * _S3654;
    float _S3657 = _S3600 * _S3653;
    float _S3658 = _S3601 * _S3653;
    float _S3659 = se_3 * _S3652.differential_0;
    float _S3660 = - _S3659;
    float _S3661 = _S3594 * _S3660;
    float _S3662 = _S3592 * _S3659;
    float _S3663 = se_3 * _S3651.differential_0;
    float _S3664 = - _S3663;
    float _S3665 = _S3591 * _S3664;
    float _S3666 = _S3593 * _S3663;
    float _S3667 = - _S3643.differential_0;
    float _S3668 = _S3599 * _S3660 + _S3593 * _S3667;
    float _S3669 = _S3596 * _S3663 + _S3594 * _S3667;
    float _S3670 = _S3598 * _S3659 + _S3591 * _S3643.differential_0;
    float _S3671 = _S3597 * _S3664 + _S3592 * _S3643.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3672;
    (&_S3672)->primal_0 = _S3589;
    (&_S3672)->differential_0 = _S3639;
    s_bwd_length_impl_1(&_S3672, _S3646.differential_0);
    float2  _S3673 = - _S3672.differential_0;
    float2  _S3674 = e2_5 * _S3673;
    float2  _S3675 = _S3588 * _S3673;
    float _S3676 = _S3674.x + _S3674.y;
    DiffPair_float_0 _S3677;
    (&_S3677)->primal_0 = _S3585;
    (&_S3677)->differential_0 = 0.0f;
    DiffPair_float_0 _S3678;
    (&_S3678)->primal_0 = 0.0f;
    (&_S3678)->differential_0 = 0.0f;
    DiffPair_float_0 _S3679;
    (&_S3679)->primal_0 = 1.0f;
    (&_S3679)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3677, &_S3678, &_S3679, _S3676);
    float _S3680 = _S3677.differential_0 / _S3586;
    float _S3681 = _S3583 * - _S3680;
    float _S3682 = _S3584 * _S3680;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3683;
    (&_S3683)->primal_0 = e2_5;
    (&_S3683)->differential_0 = _S3639;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3684;
    (&_S3684)->primal_0 = e2_5;
    (&_S3684)->differential_0 = _S3639;
    s_bwd_prop_dot_1(&_S3683, &_S3684, _S3681);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3685;
    (&_S3685)->primal_0 = _S3582;
    (&_S3685)->differential_0 = _S3639;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3686;
    (&_S3686)->primal_0 = e2_5;
    (&_S3686)->differential_0 = _S3639;
    s_bwd_prop_dot_1(&_S3685, &_S3686, _S3682);
    float2  _S3687 = - (_S3672.differential_0 + _S3685.differential_0 + make_float2 (_S3658, _S3656));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3688;
    (&_S3688)->primal_0 = _S3580;
    (&_S3688)->differential_0 = _S3639;
    s_bwd_length_impl_1(&_S3688, _S3648.differential_0);
    float2  _S3689 = - _S3688.differential_0;
    float2  _S3690 = e1_9 * _S3689;
    float2  _S3691 = _S3579 * _S3689;
    float _S3692 = _S3690.x + _S3690.y;
    DiffPair_float_0 _S3693;
    (&_S3693)->primal_0 = _S3576;
    (&_S3693)->differential_0 = 0.0f;
    DiffPair_float_0 _S3694;
    (&_S3694)->primal_0 = 0.0f;
    (&_S3694)->differential_0 = 0.0f;
    DiffPair_float_0 _S3695;
    (&_S3695)->primal_0 = 1.0f;
    (&_S3695)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3693, &_S3694, &_S3695, _S3692);
    float _S3696 = _S3693.differential_0 / _S3577;
    float _S3697 = _S3574 * - _S3696;
    float _S3698 = _S3575 * _S3696;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3699;
    (&_S3699)->primal_0 = e1_9;
    (&_S3699)->differential_0 = _S3639;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3700;
    (&_S3700)->primal_0 = e1_9;
    (&_S3700)->differential_0 = _S3639;
    s_bwd_prop_dot_1(&_S3699, &_S3700, _S3697);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3701;
    (&_S3701)->primal_0 = _S3573;
    (&_S3701)->differential_0 = _S3639;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3702;
    (&_S3702)->primal_0 = e1_9;
    (&_S3702)->differential_0 = _S3639;
    s_bwd_prop_dot_1(&_S3701, &_S3702, _S3698);
    float2  _S3703 = - (_S3688.differential_0 + _S3701.differential_0 + make_float2 (_S3662, _S3661));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3704;
    (&_S3704)->primal_0 = _S3571;
    (&_S3704)->differential_0 = _S3639;
    s_bwd_length_impl_1(&_S3704, _S3647.differential_0);
    float2  _S3705 = - _S3704.differential_0;
    float2  _S3706 = e0_9 * _S3705;
    float2  _S3707 = _S3570 * _S3705;
    float _S3708 = _S3706.x + _S3706.y;
    DiffPair_float_0 _S3709;
    (&_S3709)->primal_0 = _S3567;
    (&_S3709)->differential_0 = 0.0f;
    DiffPair_float_0 _S3710;
    (&_S3710)->primal_0 = 0.0f;
    (&_S3710)->differential_0 = 0.0f;
    DiffPair_float_0 _S3711;
    (&_S3711)->primal_0 = 1.0f;
    (&_S3711)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3709, &_S3710, &_S3711, _S3708);
    float _S3712 = _S3709.differential_0 / _S3568;
    float _S3713 = _S3565 * - _S3712;
    float _S3714 = _S3566 * _S3712;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3715;
    (&_S3715)->primal_0 = e0_9;
    (&_S3715)->differential_0 = _S3639;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3716;
    (&_S3716)->primal_0 = e0_9;
    (&_S3716)->differential_0 = _S3639;
    s_bwd_prop_dot_1(&_S3715, &_S3716, _S3713);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3717;
    (&_S3717)->primal_0 = _S3564;
    (&_S3717)->differential_0 = _S3639;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3718;
    (&_S3718)->primal_0 = e0_9;
    (&_S3718)->differential_0 = _S3639;
    s_bwd_prop_dot_1(&_S3717, &_S3718, _S3714);
    float2  _S3719 = - (_S3704.differential_0 + _S3717.differential_0 + make_float2 (_S3666, _S3665));
    float2  _S3720 = _S3640.differential_0 + _S3675 + _S3684.differential_0 + _S3683.differential_0 + _S3686.differential_0 + make_float2 (_S3655, _S3657);
    float2  _S3721 = - _S3720;
    float2  _S3722 = _S3641.differential_0 + _S3691 + _S3700.differential_0 + _S3699.differential_0 + _S3702.differential_0 + make_float2 (_S3668, _S3670);
    float2  _S3723 = - _S3722;
    float2  _S3724 = _S3642.differential_0 + _S3707 + _S3716.differential_0 + _S3715.differential_0 + _S3718.differential_0 + make_float2 (_S3671, _S3669);
    float2  _S3725 = - _S3724;
    float2  _S3726 = make_float2 (_S3622, _S3633.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S3726;
    float2  _S3727 = _S3687 + _S3721 + _S3722;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S3727;
    float2  _S3728 = _S3703 + _S3723 + _S3724;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S3728;
    float2  _S3729 = _S3719 + _S3720 + _S3725;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S3729;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3730, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3731, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3732, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3733, float2  _S3734, float _S3735)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S3730, _S3731, _S3732, _S3733, _S3734, _S3735);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S3736 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S3736;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S3736;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S3736;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S3736;
    s_bwd_evaluate_alpha_opaque_triangle_precise_0(&dp_v0_1, &dp_v1_1, &dp_v2_1, &dp_hardness_1, p_5, v_alpha_2);
    *v_v0_1 = dp_v0_1.differential_0;
    *v_v1_1 = dp_v2_1.differential_0;
    *v_v2_1 = dp_v1_1.differential_0;
    *v_hardness_3 = dp_hardness_1.differential_0;
    return;
}

inline __device__ void evaluate_color_opaque_triangle(float2  v0_4, float2  v1_4, float2  v2_4, FixedArray<float3 , 3>  * colors_0, float3  depths_0, float2  p_6, float3  * color_6, float * depth_15)
{
    *color_6 = ((*colors_0)[int(0)] + (*colors_0)[int(1)] + (*colors_0)[int(2)]) / make_float3 (3.0f);
    *depth_15 = (depths_0.x + depths_0.y + depths_0.z) / 3.0f;
    return;
}

struct DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0
{
    FixedArray<float3 , 3>  primal_0;
    FixedArray<float3 , 3>  differential_0;
};

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_2, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpcolors_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepths_0, float2  p_7, float3  dpcolor_0, float dpdepth_1)
{
    float _S3737 = 0.3333333432674408f * dpdepth_1;
    float3  _S3738 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S3739 = make_float3 (0.0f);
    float3  _S3740 = _S3739;
    *&((&_S3740)->z) = _S3737;
    *&((&_S3740)->y) = _S3737;
    *&((&_S3740)->x) = _S3737;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S3740;
    FixedArray<float3 , 3>  _S3741;
    _S3741[int(0)] = _S3739;
    _S3741[int(1)] = _S3739;
    _S3741[int(2)] = _S3739;
    _S3741[int(2)] = _S3738;
    _S3741[int(1)] = _S3738;
    _S3741[int(0)] = _S3738;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S3741;
    float2  _S3742 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S3742;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S3742;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S3742;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3743, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3744, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3745, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S3746, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3747, float2  _S3748, float3  _S3749, float _S3750)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S3743, _S3744, _S3745, _S3746, _S3747, _S3748, _S3749, _S3750);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_8, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S3751 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S3751;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S3751;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S3751;
    float3  _S3752 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S3753 = { _S3752, _S3752, _S3752 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S3753;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S3752;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_8, v_color_0, v_depth_8);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_25, float4  quat_26, float3  scale_25, float2  hardness_10, FixedArray<float3 , 16>  * sh_coeffs_21, FixedArray<float3 , 2>  * ch_coeffs_6, Matrix<float, 3, 3>  R_25, float3  t_24, float fx_25, float fy_25, float cx_25, float cy_25, float4  radial_coeffs_30, float2  tangential_coeffs_30, float2  thin_prism_coeffs_30, uint image_width_21, uint image_height_21, float near_plane_14, float far_plane_14, int4  * aabb_xyxy_14, float * depth_16, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_4)
{
    for(;;)
    {
        float3  mean_c_21 = mul_0(R_25, mean_25) + t_24;
        float _S3754 = mean_c_21.z;
        bool _S3755;
        if(_S3754 < near_plane_14)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = _S3754 > far_plane_14;
        }
        if(_S3755)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3756 = scale_25.x;
        float sx_4 = (F32_exp((_S3756)));
        float _S3757 = scale_25.y;
        float sy_4 = (F32_exp((_S3757)));
        float sz_6 = scale_25.z - 0.5f * (_S3756 + _S3757);
        float x_52 = quat_26.y;
        float x2_26 = x_52 * x_52;
        float y2_26 = quat_26.z * quat_26.z;
        float z2_47 = quat_26.w * quat_26.w;
        float xy_26 = quat_26.y * quat_26.z;
        float xz_26 = quat_26.y * quat_26.w;
        float yz_26 = quat_26.z * quat_26.w;
        float wx_26 = quat_26.x * quat_26.y;
        float wy_26 = quat_26.x * quat_26.z;
        float wz_26 = quat_26.x * quat_26.w;
        Matrix<float, 3, 3>  _S3758 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_47), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_47), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
        float3  vert0_2 = mul_0(_S3758, make_float3 (sx_4, 0.0f, 0.0f)) + mean_25;
        float3  vert1_2 = mul_0(_S3758, make_float3 (sx_4 * (-0.5f + sz_6), sy_4, 0.0f)) + mean_25;
        float3  vert2_2 = mul_0(_S3758, make_float3 (sx_4 * (-0.5f - sz_6), - sy_4, 0.0f)) + mean_25;
        float3  vert0_c_6 = mul_0(R_25, vert0_2) + t_24;
        float3  vert1_c_6 = mul_0(R_25, vert1_2) + t_24;
        float3  vert2_c_6 = mul_0(R_25, vert2_2) + t_24;
        float _S3759 = vert0_c_6.z;
        float _S3760 = vert1_c_6.z;
        float _S3761 = vert2_c_6.z;
        if(_S3759 < near_plane_14)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = _S3759 > far_plane_14;
        }
        if(_S3755)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = _S3760 < near_plane_14;
        }
        if(_S3755)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = _S3760 > far_plane_14;
        }
        if(_S3755)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = _S3761 < near_plane_14;
        }
        if(_S3755)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = _S3761 > far_plane_14;
        }
        if(_S3755)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  _S3762 = make_float2 (fx_25, fy_25);
        float2  _S3763 = make_float2 (cx_25, cy_25);
        float2  _S3764 = _S3762 * (float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S3759)) + _S3763;
        float2  _S3765 = _S3762 * (float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S3760)) + _S3763;
        float2  _S3766 = _S3762 * (float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S3761)) + _S3763;
        float2  e0_10 = _S3765 - _S3764;
        float2  e1_10 = _S3766 - _S3765;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(_S3764 - _S3766)));
        float _S3767 = _S3764.x;
        float _S3768 = _S3765.x;
        float _S3769 = _S3766.x;
        float xmax_7 = (F32_max(((F32_max((_S3767), (_S3768)))), (_S3769))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S3767), (_S3768)))), (_S3769))) - offset_4;
        float _S3770 = _S3764.y;
        float _S3771 = _S3765.y;
        float _S3772 = _S3766.y;
        float ymax_7 = (F32_max(((F32_max((_S3770), (_S3771)))), (_S3772))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S3770), (_S3771)))), (_S3772))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = xmin_7 >= float(image_width_21);
        }
        if(_S3755)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = ymax_7 <= 0.0f;
        }
        if(_S3755)
        {
            _S3755 = true;
        }
        else
        {
            _S3755 = ymin_7 >= float(image_height_21);
        }
        if(_S3755)
        {
            _S3755 = true;
        }
        else
        {
            if(_S3754 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S3755 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S3755 = false;
                }
                if(_S3755)
                {
                    _S3755 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S3755 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S3755 = false;
                    }
                }
            }
            else
            {
                _S3755 = false;
            }
        }
        if(_S3755)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S3773 = mean_25 - - mul_0(transpose_0(R_25), t_24);
        float _S3774 = _S3773.x;
        float _S3775 = _S3773.y;
        float _S3776 = _S3773.z;
        float norm_14 = (F32_sqrt((_S3774 * _S3774 + _S3775 * _S3775 + _S3776 * _S3776)));
        float x_53 = _S3774 / norm_14;
        float y_24 = _S3775 / norm_14;
        float z_21 = _S3776 / norm_14;
        float z2_48 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_53 * x_53 - y_24 * y_24;
        float fS1_21 = 2.0f * x_53 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_48 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_53) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_48 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_53) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_53 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_48 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_53) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_53 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S3777 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S3777);
        float3  _S3778 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S3779 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S3778 + _S3779 + make_float3 (0.5f), _S3777);
        (*rgbs_0)[int(2)] = max_0(_S3778 - _S3779 + make_float3 (0.5f), _S3777);
        (*verts_0)[int(0)] = vert0_2;
        (*verts_0)[int(1)] = vert1_2;
        (*verts_0)[int(2)] = vert2_2;
        float3  _S3780 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S3780 * make_float3 (float(- (F32_sign((dot_0(_S3780, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_26, float3  t_25, float fx_26, float fy_26, float cx_26, float cy_26, float4  radial_coeffs_31, float2  tangential_coeffs_31, float2  thin_prism_coeffs_31, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_26, mean_26) + t_25;
        float _S3781 = length_1(mean_c_22);
        bool _S3782;
        if(_S3781 < near_plane_15)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = _S3781 > far_plane_15;
        }
        if(_S3782)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3783 = scale_26.x;
        float sx_5 = (F32_exp((_S3783)));
        float _S3784 = scale_26.y;
        float sy_5 = (F32_exp((_S3784)));
        float sz_7 = scale_26.z - 0.5f * (_S3783 + _S3784);
        float x_54 = quat_27.y;
        float x2_27 = x_54 * x_54;
        float y2_27 = quat_27.z * quat_27.z;
        float z2_49 = quat_27.w * quat_27.w;
        float xy_27 = quat_27.y * quat_27.z;
        float xz_27 = quat_27.y * quat_27.w;
        float yz_27 = quat_27.z * quat_27.w;
        float wx_27 = quat_27.x * quat_27.y;
        float wy_27 = quat_27.x * quat_27.z;
        float wz_27 = quat_27.x * quat_27.w;
        Matrix<float, 3, 3>  _S3785 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_49), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_49), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S3785, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S3785, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S3785, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_7 = mul_0(R_26, vert0_3) + t_25;
        float3  vert1_c_7 = mul_0(R_26, vert1_3) + t_25;
        float3  vert2_c_7 = mul_0(R_26, vert2_3) + t_25;
        float _S3786 = length_1(vert0_c_7);
        float _S3787 = length_1(vert1_c_7);
        float _S3788 = length_1(vert2_c_7);
        if(_S3786 < near_plane_15)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = _S3786 > far_plane_15;
        }
        if(_S3782)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = _S3787 < near_plane_15;
        }
        if(_S3782)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = _S3787 > far_plane_15;
        }
        if(_S3782)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = _S3788 < near_plane_15;
        }
        if(_S3782)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = _S3788 > far_plane_15;
        }
        if(_S3782)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_4 = CameraDistortion_x24init_0(radial_coeffs_31, tangential_coeffs_31, thin_prism_coeffs_31);
        float2  _S3789 = float2 {vert0_c_7.x, vert0_c_7.y};
        float r_14 = length_0(_S3789);
        float _S3790 = vert0_c_7.z;
        float theta_10 = (F32_atan2((r_14), (_S3790)));
        float k_7;
        if(theta_10 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_10 * theta_10 / 3.0f) / _S3790;
        }
        else
        {
            k_7 = theta_10 / r_14;
        }
        float2  _S3791 = _S3789 * make_float2 (k_7);
        float k1_6 = dist_coeffs_4.radial_coeffs_0.x;
        float k2_6 = dist_coeffs_4.radial_coeffs_0.y;
        float k3_6 = dist_coeffs_4.radial_coeffs_0.z;
        float k4_7 = dist_coeffs_4.radial_coeffs_0.w;
        float p1_8 = dist_coeffs_4.tangential_coeffs_0.x;
        float p2_8 = dist_coeffs_4.tangential_coeffs_0.y;
        float sx1_7 = dist_coeffs_4.thin_prism_coeffs_0.x;
        float sy1_7 = dist_coeffs_4.thin_prism_coeffs_0.y;
        float u_20 = _S3791.x;
        float v_20 = _S3791.y;
        float r2_20 = u_20 * u_20 + v_20 * v_20;
        float _S3792 = 2.0f * p1_8;
        float _S3793 = 2.0f * p2_8;
        float2  _S3794 = _S3791 * make_float2 (1.0f + r2_20 * (k1_6 + r2_20 * (k2_6 + r2_20 * (k3_6 + r2_20 * k4_7)))) + make_float2 (_S3792 * u_20 * v_20 + p2_8 * (r2_20 + 2.0f * u_20 * u_20) + sx1_7 * r2_20, _S3793 * u_20 * v_20 + p1_8 * (r2_20 + 2.0f * v_20 * v_20) + sy1_7 * r2_20);
        float _S3795 = fx_26 * _S3794.x + cx_26;
        float _S3796 = fy_26 * _S3794.y + cy_26;
        float2  _S3797 = make_float2 (_S3795, _S3796);
        float2  _S3798 = float2 {vert1_c_7.x, vert1_c_7.y};
        float r_15 = length_0(_S3798);
        float _S3799 = vert1_c_7.z;
        float theta_11 = (F32_atan2((r_15), (_S3799)));
        if(theta_11 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_11 * theta_11 / 3.0f) / _S3799;
        }
        else
        {
            k_7 = theta_11 / r_15;
        }
        float2  _S3800 = _S3798 * make_float2 (k_7);
        float u_21 = _S3800.x;
        float v_21 = _S3800.y;
        float r2_21 = u_21 * u_21 + v_21 * v_21;
        float2  _S3801 = _S3800 * make_float2 (1.0f + r2_21 * (k1_6 + r2_21 * (k2_6 + r2_21 * (k3_6 + r2_21 * k4_7)))) + make_float2 (_S3792 * u_21 * v_21 + p2_8 * (r2_21 + 2.0f * u_21 * u_21) + sx1_7 * r2_21, _S3793 * u_21 * v_21 + p1_8 * (r2_21 + 2.0f * v_21 * v_21) + sy1_7 * r2_21);
        float _S3802 = fx_26 * _S3801.x + cx_26;
        float _S3803 = fy_26 * _S3801.y + cy_26;
        float2  _S3804 = make_float2 (_S3802, _S3803);
        float2  _S3805 = float2 {vert2_c_7.x, vert2_c_7.y};
        float r_16 = length_0(_S3805);
        float _S3806 = vert2_c_7.z;
        float theta_12 = (F32_atan2((r_16), (_S3806)));
        if(theta_12 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_12 * theta_12 / 3.0f) / _S3806;
        }
        else
        {
            k_7 = theta_12 / r_16;
        }
        float2  _S3807 = _S3805 * make_float2 (k_7);
        float u_22 = _S3807.x;
        float v_22 = _S3807.y;
        float r2_22 = u_22 * u_22 + v_22 * v_22;
        float2  _S3808 = _S3807 * make_float2 (1.0f + r2_22 * (k1_6 + r2_22 * (k2_6 + r2_22 * (k3_6 + r2_22 * k4_7)))) + make_float2 (_S3792 * u_22 * v_22 + p2_8 * (r2_22 + 2.0f * u_22 * u_22) + sx1_7 * r2_22, _S3793 * u_22 * v_22 + p1_8 * (r2_22 + 2.0f * v_22 * v_22) + sy1_7 * r2_22);
        float _S3809 = fx_26 * _S3808.x + cx_26;
        float _S3810 = fy_26 * _S3808.y + cy_26;
        float2  _S3811 = make_float2 (_S3809, _S3810);
        float2  e0_11 = _S3804 - _S3797;
        float2  e1_11 = _S3811 - _S3804;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(_S3797 - _S3811)));
        float xmax_8 = (F32_max(((F32_max((_S3795), (_S3802)))), (_S3809))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S3795), (_S3802)))), (_S3809))) - offset_5;
        float ymax_8 = (F32_max(((F32_max((_S3796), (_S3803)))), (_S3810))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S3796), (_S3803)))), (_S3810))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = xmin_8 >= float(image_width_22);
        }
        if(_S3782)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = ymax_8 <= 0.0f;
        }
        if(_S3782)
        {
            _S3782 = true;
        }
        else
        {
            _S3782 = ymin_8 >= float(image_height_22);
        }
        if(_S3782)
        {
            _S3782 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S3782 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S3782 = false;
                }
                if(_S3782)
                {
                    _S3782 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S3782 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S3782 = false;
                    }
                }
            }
            else
            {
                _S3782 = false;
            }
        }
        if(_S3782)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S3812 = mean_26 - - mul_0(transpose_0(R_26), t_25);
        float _S3813 = _S3812.x;
        float _S3814 = _S3812.y;
        float _S3815 = _S3812.z;
        float norm_15 = (F32_sqrt((_S3813 * _S3813 + _S3814 * _S3814 + _S3815 * _S3815)));
        float x_55 = _S3813 / norm_15;
        float y_25 = _S3814 / norm_15;
        float z_22 = _S3815 / norm_15;
        float z2_50 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_55 * x_55 - y_25 * y_25;
        float fS1_22 = 2.0f * x_55 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_50 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_55) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_50 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_55) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_55 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_50 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_55) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_55 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S3816 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S3816);
        float3  _S3817 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S3818 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S3817 + _S3818 + make_float3 (0.5f), _S3816);
        (*rgbs_1)[int(2)] = max_0(_S3817 - _S3818 + make_float3 (0.5f), _S3816);
        (*verts_1)[int(0)] = vert0_3;
        (*verts_1)[int(1)] = vert1_3;
        (*verts_1)[int(2)] = vert2_3;
        float3  _S3819 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S3819 * make_float3 (float(- (F32_sign((dot_0(_S3819, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_27, float3  t_26, float fx_27, float fy_27, float cx_27, float cy_27, float4  radial_coeffs_32, float2  tangential_coeffs_32, float2  thin_prism_coeffs_32, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_27, mean_27) + t_26;
    float _S3820 = scale_27.x;
    float sx_6 = (F32_exp((_S3820)));
    float _S3821 = scale_27.y;
    float sy_6 = (F32_exp((_S3821)));
    float sz_8 = scale_27.z - 0.5f * (_S3820 + _S3821);
    float x_56 = quat_28.y;
    float x2_28 = x_56 * x_56;
    float y2_28 = quat_28.z * quat_28.z;
    float z2_51 = quat_28.w * quat_28.w;
    float xy_28 = quat_28.y * quat_28.z;
    float xz_28 = quat_28.y * quat_28.w;
    float yz_28 = quat_28.z * quat_28.w;
    float wx_28 = quat_28.x * quat_28.y;
    float wy_28 = quat_28.x * quat_28.z;
    float wz_28 = quat_28.x * quat_28.w;
    Matrix<float, 3, 3>  _S3822 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_51), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_51), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
    float3  vert0_4 = mul_0(_S3822, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
    float3  vert1_4 = mul_0(_S3822, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
    float3  vert2_4 = mul_0(_S3822, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
    float3  vert0_c_8 = mul_0(R_27, vert0_4) + t_26;
    float3  vert1_c_8 = mul_0(R_27, vert1_4) + t_26;
    float3  vert2_c_8 = mul_0(R_27, vert2_4) + t_26;
    float2  _S3823 = make_float2 (fx_27, fy_27);
    float2  _S3824 = make_float2 (cx_27, cy_27);
    float2  _S3825 = _S3823 * (float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z)) + _S3824;
    float2  _S3826 = _S3823 * (float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z)) + _S3824;
    float2  _S3827 = _S3823 * (float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z)) + _S3824;
    float2  e0_12 = _S3826 - _S3825;
    float2  e1_12 = _S3827 - _S3826;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(_S3825 - _S3827)));
    float _S3828 = _S3825.x;
    float _S3829 = _S3826.x;
    float _S3830 = _S3827.x;
    float _S3831 = _S3825.y;
    float _S3832 = _S3826.y;
    float _S3833 = _S3827.y;
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S3828), (_S3829)))), (_S3830))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S3831), (_S3832)))), (_S3833))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S3828), (_S3829)))), (_S3830))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S3831), (_S3832)))), (_S3833))) + offset_6)))));
    *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
    float3  _S3834 = mean_27 - - mul_0(transpose_0(R_27), t_26);
    float _S3835 = _S3834.x;
    float _S3836 = _S3834.y;
    float _S3837 = _S3834.z;
    float norm_16 = (F32_sqrt((_S3835 * _S3835 + _S3836 * _S3836 + _S3837 * _S3837)));
    float x_57 = _S3835 / norm_16;
    float y_26 = _S3836 / norm_16;
    float z_23 = _S3837 / norm_16;
    float z2_52 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_57 * x_57 - y_26 * y_26;
    float fS1_23 = 2.0f * x_57 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_52 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_57) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_52 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_57) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_57 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_52 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_57) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_57 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S3838 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S3838);
    float3  _S3839 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S3840 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S3839 + _S3840 + make_float3 (0.5f), _S3838);
    (*rgbs_2)[int(2)] = max_0(_S3839 - _S3840 + make_float3 (0.5f), _S3838);
    (*verts_2)[int(0)] = vert0_4;
    (*verts_2)[int(1)] = vert1_4;
    (*verts_2)[int(2)] = vert2_4;
    float3  _S3841 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S3841 * make_float3 (float(- (F32_sign((dot_0(_S3841, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_28, float3  t_27, float fx_28, float fy_28, float cx_28, float cy_28, float4  radial_coeffs_33, float2  tangential_coeffs_33, float2  thin_prism_coeffs_33, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_28, mean_28) + t_27;
    float _S3842 = scale_28.x;
    float sx_7 = (F32_exp((_S3842)));
    float _S3843 = scale_28.y;
    float sy_7 = (F32_exp((_S3843)));
    float sz_9 = scale_28.z - 0.5f * (_S3842 + _S3843);
    float x_58 = quat_29.y;
    float x2_29 = x_58 * x_58;
    float y2_29 = quat_29.z * quat_29.z;
    float z2_53 = quat_29.w * quat_29.w;
    float xy_29 = quat_29.y * quat_29.z;
    float xz_29 = quat_29.y * quat_29.w;
    float yz_29 = quat_29.z * quat_29.w;
    float wx_29 = quat_29.x * quat_29.y;
    float wy_29 = quat_29.x * quat_29.z;
    float wz_29 = quat_29.x * quat_29.w;
    Matrix<float, 3, 3>  _S3844 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_53), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_53), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S3844, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S3844, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S3844, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_9 = mul_0(R_28, vert0_5) + t_27;
    float3  vert1_c_9 = mul_0(R_28, vert1_5) + t_27;
    float3  vert2_c_9 = mul_0(R_28, vert2_5) + t_27;
    CameraDistortion_0 dist_coeffs_5 = CameraDistortion_x24init_0(radial_coeffs_33, tangential_coeffs_33, thin_prism_coeffs_33);
    float2  _S3845 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_17 = length_0(_S3845);
    float _S3846 = vert0_c_9.z;
    float theta_13 = (F32_atan2((r_17), (_S3846)));
    float k_8;
    if(theta_13 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_13 * theta_13 / 3.0f) / _S3846;
    }
    else
    {
        k_8 = theta_13 / r_17;
    }
    float2  _S3847 = _S3845 * make_float2 (k_8);
    float k1_7 = dist_coeffs_5.radial_coeffs_0.x;
    float k2_7 = dist_coeffs_5.radial_coeffs_0.y;
    float k3_7 = dist_coeffs_5.radial_coeffs_0.z;
    float k4_8 = dist_coeffs_5.radial_coeffs_0.w;
    float p1_9 = dist_coeffs_5.tangential_coeffs_0.x;
    float p2_9 = dist_coeffs_5.tangential_coeffs_0.y;
    float sx1_8 = dist_coeffs_5.thin_prism_coeffs_0.x;
    float sy1_8 = dist_coeffs_5.thin_prism_coeffs_0.y;
    float u_23 = _S3847.x;
    float v_23 = _S3847.y;
    float r2_23 = u_23 * u_23 + v_23 * v_23;
    float _S3848 = 2.0f * p1_9;
    float _S3849 = 2.0f * p2_9;
    float2  _S3850 = _S3847 * make_float2 (1.0f + r2_23 * (k1_7 + r2_23 * (k2_7 + r2_23 * (k3_7 + r2_23 * k4_8)))) + make_float2 (_S3848 * u_23 * v_23 + p2_9 * (r2_23 + 2.0f * u_23 * u_23) + sx1_8 * r2_23, _S3849 * u_23 * v_23 + p1_9 * (r2_23 + 2.0f * v_23 * v_23) + sy1_8 * r2_23);
    float _S3851 = fx_28 * _S3850.x + cx_28;
    float _S3852 = fy_28 * _S3850.y + cy_28;
    float2  _S3853 = make_float2 (_S3851, _S3852);
    float2  _S3854 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_18 = length_0(_S3854);
    float _S3855 = vert1_c_9.z;
    float theta_14 = (F32_atan2((r_18), (_S3855)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_14 * theta_14 / 3.0f) / _S3855;
    }
    else
    {
        k_8 = theta_14 / r_18;
    }
    float2  _S3856 = _S3854 * make_float2 (k_8);
    float u_24 = _S3856.x;
    float v_24 = _S3856.y;
    float r2_24 = u_24 * u_24 + v_24 * v_24;
    float2  _S3857 = _S3856 * make_float2 (1.0f + r2_24 * (k1_7 + r2_24 * (k2_7 + r2_24 * (k3_7 + r2_24 * k4_8)))) + make_float2 (_S3848 * u_24 * v_24 + p2_9 * (r2_24 + 2.0f * u_24 * u_24) + sx1_8 * r2_24, _S3849 * u_24 * v_24 + p1_9 * (r2_24 + 2.0f * v_24 * v_24) + sy1_8 * r2_24);
    float _S3858 = fx_28 * _S3857.x + cx_28;
    float _S3859 = fy_28 * _S3857.y + cy_28;
    float2  _S3860 = make_float2 (_S3858, _S3859);
    float2  _S3861 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_19 = length_0(_S3861);
    float _S3862 = vert2_c_9.z;
    float theta_15 = (F32_atan2((r_19), (_S3862)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_15 * theta_15 / 3.0f) / _S3862;
    }
    else
    {
        k_8 = theta_15 / r_19;
    }
    float2  _S3863 = _S3861 * make_float2 (k_8);
    float u_25 = _S3863.x;
    float v_25 = _S3863.y;
    float r2_25 = u_25 * u_25 + v_25 * v_25;
    float2  _S3864 = _S3863 * make_float2 (1.0f + r2_25 * (k1_7 + r2_25 * (k2_7 + r2_25 * (k3_7 + r2_25 * k4_8)))) + make_float2 (_S3848 * u_25 * v_25 + p2_9 * (r2_25 + 2.0f * u_25 * u_25) + sx1_8 * r2_25, _S3849 * u_25 * v_25 + p1_9 * (r2_25 + 2.0f * v_25 * v_25) + sy1_8 * r2_25);
    float _S3865 = fx_28 * _S3864.x + cx_28;
    float _S3866 = fy_28 * _S3864.y + cy_28;
    float2  _S3867 = make_float2 (_S3865, _S3866);
    float2  e0_13 = _S3860 - _S3853;
    float2  e1_13 = _S3867 - _S3860;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(_S3853 - _S3867)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S3851), (_S3858)))), (_S3865))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S3852), (_S3859)))), (_S3866))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S3851), (_S3858)))), (_S3865))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S3852), (_S3859)))), (_S3866))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S3868 = mean_28 - - mul_0(transpose_0(R_28), t_27);
    float _S3869 = _S3868.x;
    float _S3870 = _S3868.y;
    float _S3871 = _S3868.z;
    float norm_17 = (F32_sqrt((_S3869 * _S3869 + _S3870 * _S3870 + _S3871 * _S3871)));
    float x_59 = _S3869 / norm_17;
    float y_27 = _S3870 / norm_17;
    float z_24 = _S3871 / norm_17;
    float z2_54 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_59 * x_59 - y_27 * y_27;
    float fS1_24 = 2.0f * x_59 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_54 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_59) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_54 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_59) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_59 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_54 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_59) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_59 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S3872 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S3872);
    float3  _S3873 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S3874 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S3873 + _S3874 + make_float3 (0.5f), _S3872);
    (*rgbs_3)[int(2)] = max_0(_S3873 - _S3874 + make_float3 (0.5f), _S3872);
    (*verts_3)[int(0)] = vert0_5;
    (*verts_3)[int(1)] = vert1_5;
    (*verts_3)[int(2)] = vert2_5;
    float3  _S3875 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S3875 * make_float3 (float(- (F32_sign((dot_0(_S3875, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_29, float3  t_28, float fx_29, float fy_29, float cx_29, float cy_29, float4  radial_coeffs_34, float2  tangential_coeffs_34, float2  thin_prism_coeffs_34, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_1(R_29, mean_29) + t_28;
    float _S3876 = scale_29.x;
    float _S3877 = s_primal_ctx_exp_1(_S3876);
    float _S3878 = scale_29.y;
    float _S3879 = s_primal_ctx_exp_1(_S3878);
    float sz_10 = scale_29.z - 0.5f * (_S3876 + _S3878);
    float _S3880 = quat_30.y;
    float x2_30 = _S3880 * _S3880;
    float y2_30 = quat_30.z * quat_30.z;
    float z2_55 = quat_30.w * quat_30.w;
    float xy_30 = quat_30.y * quat_30.z;
    float xz_30 = quat_30.y * quat_30.w;
    float yz_30 = quat_30.z * quat_30.w;
    float wx_30 = quat_30.x * quat_30.y;
    float wy_30 = quat_30.x * quat_30.z;
    float wz_30 = quat_30.x * quat_30.w;
    Matrix<float, 3, 3>  _S3881 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_55), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_55), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  _S3882 = make_float3 (_S3877, 0.0f, 0.0f);
    float3  vert0_6 = s_primal_ctx_mul_1(_S3881, _S3882) + mean_29;
    float _S3883 = -0.5f + sz_10;
    float3  _S3884 = make_float3 (_S3877 * _S3883, _S3879, 0.0f);
    float3  vert1_6 = s_primal_ctx_mul_1(_S3881, _S3884) + mean_29;
    float _S3885 = -0.5f - sz_10;
    float3  _S3886 = make_float3 (_S3877 * _S3885, - _S3879, 0.0f);
    float3  vert2_6 = s_primal_ctx_mul_1(_S3881, _S3886) + mean_29;
    float3  vert0_c_10 = s_primal_ctx_mul_1(R_29, vert0_6) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_1(R_29, vert1_6) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_1(R_29, vert2_6) + t_28;
    float2  _S3887 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S3888 = vert0_c_10.z;
    float2  _S3889 = make_float2 (_S3888);
    float2  _S3890 = make_float2 (_S3888 * _S3888);
    float2  _S3891 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S3892 = vert1_c_10.z;
    float2  _S3893 = make_float2 (_S3892);
    float2  _S3894 = make_float2 (_S3892 * _S3892);
    float2  _S3895 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S3896 = vert2_c_10.z;
    float2  _S3897 = make_float2 (_S3896);
    float2  _S3898 = make_float2 (_S3896 * _S3896);
    float2  _S3899 = make_float2 (fx_29, fy_29);
    float2  _S3900 = make_float2 (cx_29, cy_29);
    float2  _S3901 = _S3899 * (_S3887 / make_float2 (_S3888)) + _S3900;
    float2  _S3902 = _S3899 * (_S3891 / make_float2 (_S3892)) + _S3900;
    float2  _S3903 = _S3899 * (_S3895 / make_float2 (_S3896)) + _S3900;
    float2  e0_14 = _S3902 - _S3901;
    float2  e1_14 = _S3903 - _S3902;
    float2  e2_6 = _S3901 - _S3903;
    float _S3904 = e0_14.x;
    float _S3905 = e1_14.y;
    float _S3906 = e0_14.y;
    float _S3907 = e1_14.x;
    float _S3908 = _S3904 * _S3905 - _S3906 * _S3907;
    float _S3909 = 1.0f - hardness_14.y;
    float _S3910 = -1.0f / _S3909;
    float _S3911 = _S3909 * _S3909;
    float _S3912 = _S3901.x;
    float _S3913 = _S3902.x;
    float _S3914 = s_primal_ctx_max_0(_S3912, _S3913);
    float _S3915 = _S3903.x;
    float _S3916 = s_primal_ctx_min_0(_S3912, _S3913);
    float _S3917 = _S3901.y;
    float _S3918 = _S3902.y;
    float _S3919 = s_primal_ctx_max_0(_S3917, _S3918);
    float _S3920 = _S3903.y;
    float _S3921 = s_primal_ctx_min_0(_S3917, _S3918);
    float3  _S3922 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S3923 = length_1(_S3922) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S3924 = transpose_0(R_29);
    float3  _S3925 = mean_29 - - s_primal_ctx_mul_1(_S3924, t_28);
    float _S3926 = _S3925.x;
    float _S3927 = _S3925.y;
    float _S3928 = _S3925.z;
    float _S3929 = _S3926 * _S3926 + _S3927 * _S3927 + _S3928 * _S3928;
    float _S3930 = s_primal_ctx_sqrt_0(_S3929);
    float x_60 = _S3926 / _S3930;
    float3  _S3931 = make_float3 (x_60);
    float _S3932 = _S3930 * _S3930;
    float y_28 = _S3927 / _S3930;
    float z_25 = _S3928 / _S3930;
    float3  _S3933 = make_float3 (z_25);
    float _S3934 = - y_28;
    float3  _S3935 = make_float3 (_S3934);
    float z2_56 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_60 * x_60 - y_28 * y_28;
    float _S3936 = 2.0f * x_60;
    float fS1_25 = _S3936 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_56 - 0.31539157032966614f;
    float3  _S3937 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_60;
    float3  _S3938 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S3939 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S3940 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S3941 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_56 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S3942 = 1.86588168144226074f * z2_56 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S3942;
    float3  _S3943 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_60;
    float3  _S3944 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S3945 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S3946 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S3947 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_60 * fC1_25 - y_28 * fS1_25);
    float3  _S3948 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_60 * fS1_25 + y_28 * fC1_25);
    float3  _S3949 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3934) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_60) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S3950 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S3951 = make_float3 (0.0f);
    float3  _S3952 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S3953 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3954 = make_float3 (_S3953);
    float3  _S3955 = (*ch_coeffs_10)[int(1)] * make_float3 (_S3953);
    float3  _S3956 = _S3952 + _S3955 + make_float3 (0.5f);
    float3  _S3957 = _S3952 - _S3955 + make_float3 (0.5f);
    float3  _S3958 = vert1_c_10 - vert0_c_10;
    float3  _S3959 = vert2_c_10 - vert0_c_10;
    float3  _S3960 = s_primal_ctx_cross_0(_S3958, _S3959);
    float3  _S3961 = normalize_0(_S3960);
    float3  _S3962 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3961, mean_c_25)))))) * v_normal_2;
    float3  _S3963 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3964;
    (&_S3964)->primal_0 = _S3961;
    (&_S3964)->differential_0 = _S3963;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3965;
    (&_S3965)->primal_0 = mean_c_25;
    (&_S3965)->differential_0 = _S3963;
    s_bwd_prop_dot_0(&_S3964, &_S3965, 0.0f);
    float3  _S3966 = _S3962 + _S3964.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3967;
    (&_S3967)->primal_0 = _S3960;
    (&_S3967)->differential_0 = _S3963;
    s_bwd_normalize_impl_0(&_S3967, _S3966);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3968;
    (&_S3968)->primal_0 = _S3958;
    (&_S3968)->differential_0 = _S3963;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3969;
    (&_S3969)->primal_0 = _S3959;
    (&_S3969)->differential_0 = _S3963;
    s_bwd_prop_cross_0(&_S3968, &_S3969, _S3967.differential_0);
    float3  _S3970 = - _S3969.differential_0;
    float3  _S3971 = - _S3968.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3972;
    (&_S3972)->primal_0 = _S3957;
    (&_S3972)->differential_0 = _S3963;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3973;
    (&_S3973)->primal_0 = _S3951;
    (&_S3973)->differential_0 = _S3963;
    s_bwd_prop_max_0(&_S3972, &_S3973, (*v_rgbs_0)[int(2)]);
    float3  _S3974 = - _S3972.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3975;
    (&_S3975)->primal_0 = _S3956;
    (&_S3975)->differential_0 = _S3963;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3976;
    (&_S3976)->primal_0 = _S3951;
    (&_S3976)->differential_0 = _S3963;
    s_bwd_prop_max_0(&_S3975, &_S3976, (*v_rgbs_0)[int(1)]);
    float3  _S3977 = _S3954 * (_S3974 + _S3975.differential_0);
    float3  _S3978 = _S3972.differential_0 + _S3975.differential_0;
    float3  _S3979 = make_float3 (0.5f) * - _S3978;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3980;
    (&_S3980)->primal_0 = _S3950;
    (&_S3980)->differential_0 = _S3963;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3981;
    (&_S3981)->primal_0 = _S3951;
    (&_S3981)->differential_0 = _S3963;
    s_bwd_prop_max_0(&_S3980, &_S3981, (*v_rgbs_0)[int(0)]);
    float3  _S3982 = _S3979 + _S3980.differential_0;
    float3  _S3983 = _S3978 + _S3980.differential_0;
    float3  _S3984 = _S3948 * _S3983;
    float3  _S3985 = (*sh_coeffs_25)[int(15)] * _S3983;
    float3  _S3986 = _S3946 * _S3983;
    float3  _S3987 = (*sh_coeffs_25)[int(14)] * _S3983;
    float3  _S3988 = _S3944 * _S3983;
    float3  _S3989 = (*sh_coeffs_25)[int(13)] * _S3983;
    float3  _S3990 = _S3943 * _S3983;
    float3  _S3991 = (*sh_coeffs_25)[int(12)] * _S3983;
    float3  _S3992 = _S3945 * _S3983;
    float3  _S3993 = (*sh_coeffs_25)[int(11)] * _S3983;
    float3  _S3994 = _S3947 * _S3983;
    float3  _S3995 = (*sh_coeffs_25)[int(10)] * _S3983;
    float3  _S3996 = _S3949 * _S3983;
    float3  _S3997 = (*sh_coeffs_25)[int(9)] * _S3983;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S3997.x + _S3997.y + _S3997.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S3985.x + _S3985.y + _S3985.z);
    float _S3998 = _S3995.x + _S3995.y + _S3995.z;
    float _S3999 = _S3987.x + _S3987.y + _S3987.z;
    float _S4000 = _S3993.x + _S3993.y + _S3993.z;
    float _S4001 = _S3989.x + _S3989.y + _S3989.z;
    float _S4002 = _S3991.x + _S3991.y + _S3991.z;
    float _S4003 = - s_diff_fC2_T_7;
    float3  _S4004 = _S3940 * _S3983;
    float3  _S4005 = (*sh_coeffs_25)[int(8)] * _S3983;
    float3  _S4006 = _S3938 * _S3983;
    float3  _S4007 = (*sh_coeffs_25)[int(7)] * _S3983;
    float3  _S4008 = _S3937 * _S3983;
    float3  _S4009 = (*sh_coeffs_25)[int(6)] * _S3983;
    float3  _S4010 = _S3939 * _S3983;
    float3  _S4011 = (*sh_coeffs_25)[int(5)] * _S3983;
    float3  _S4012 = _S3941 * _S3983;
    float3  _S4013 = (*sh_coeffs_25)[int(4)] * _S3983;
    float _S4014 = _S4011.x + _S4011.y + _S4011.z;
    float _S4015 = _S4007.x + _S4007.y + _S4007.z;
    float _S4016 = fTmp1B_25 * _S3998 + x_60 * s_diff_fS2_T_7 + y_28 * _S4003 + 0.54627424478530884f * (_S4013.x + _S4013.y + _S4013.z);
    float _S4017 = fTmp1B_25 * _S3999 + y_28 * s_diff_fS2_T_7 + x_60 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S4005.x + _S4005.y + _S4005.z);
    float _S4018 = y_28 * - _S4017;
    float _S4019 = x_60 * _S4017;
    float _S4020 = z_25 * (1.86588168144226074f * (z_25 * _S4002) + -2.28522896766662598f * (y_28 * _S4000 + x_60 * _S4001) + 0.94617468118667603f * (_S4009.x + _S4009.y + _S4009.z));
    float3  _S4021 = make_float3 (0.48860251903533936f) * _S3983;
    float3  _S4022 = - _S4021;
    float3  _S4023 = _S3931 * _S4022;
    float3  _S4024 = (*sh_coeffs_25)[int(3)] * _S4022;
    float3  _S4025 = _S3933 * _S4021;
    float3  _S4026 = (*sh_coeffs_25)[int(2)] * _S4021;
    float3  _S4027 = _S3935 * _S4021;
    float3  _S4028 = (*sh_coeffs_25)[int(1)] * _S4021;
    float _S4029 = (_S3942 * _S4002 + 1.44530570507049561f * (fS1_25 * _S3998 + fC1_25 * _S3999) + -1.09254848957061768f * (y_28 * _S4014 + x_60 * _S4015) + _S4020 + _S4020 + _S4026.x + _S4026.y + _S4026.z) / _S3932;
    float _S4030 = _S3930 * _S4029;
    float _S4031 = (fTmp0C_25 * _S4000 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S4003 + fTmp0B_25 * _S4014 + _S3936 * _S4016 + _S4018 + _S4018 + - (_S4028.x + _S4028.y + _S4028.z)) / _S3932;
    float _S4032 = _S3930 * _S4031;
    float _S4033 = (fTmp0C_25 * _S4001 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4015 + 2.0f * (y_28 * _S4016) + _S4019 + _S4019 + _S4024.x + _S4024.y + _S4024.z) / _S3932;
    float _S4034 = _S3930 * _S4033;
    float _S4035 = _S3928 * - _S4029 + _S3927 * - _S4031 + _S3926 * - _S4033;
    DiffPair_float_0 _S4036;
    (&_S4036)->primal_0 = _S3929;
    (&_S4036)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4036, _S4035);
    float _S4037 = _S3928 * _S4036.differential_0;
    float _S4038 = _S3927 * _S4036.differential_0;
    float _S4039 = _S3926 * _S4036.differential_0;
    float3  _S4040 = make_float3 (0.282094806432724f) * _S3983;
    float3  _S4041 = make_float3 (_S4034 + _S4039 + _S4039, _S4032 + _S4038 + _S4038, _S4030 + _S4037 + _S4037);
    float3  _S4042 = - - _S4041;
    Matrix<float, 3, 3>  _S4043 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4044;
    (&_S4044)->primal_0 = _S3924;
    (&_S4044)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4045;
    (&_S4045)->primal_0 = t_28;
    (&_S4045)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4044, &_S4045, _S4042);
    Matrix<float, 3, 3>  _S4046 = transpose_0(_S4044.differential_0);
    DiffPair_float_0 _S4047;
    (&_S4047)->primal_0 = _S3923;
    (&_S4047)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4047, v_depth_9);
    float _S4048 = 0.3333333432674408f * _S4047.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4049;
    (&_S4049)->primal_0 = _S3922;
    (&_S4049)->differential_0 = _S3963;
    s_bwd_length_impl_0(&_S4049, _S4048);
    DiffPair_float_0 _S4050;
    (&_S4050)->primal_0 = _S3921;
    (&_S4050)->differential_0 = 0.0f;
    DiffPair_float_0 _S4051;
    (&_S4051)->primal_0 = _S3920;
    (&_S4051)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4050, &_S4051, 0.0f);
    DiffPair_float_0 _S4052;
    (&_S4052)->primal_0 = _S3917;
    (&_S4052)->differential_0 = 0.0f;
    DiffPair_float_0 _S4053;
    (&_S4053)->primal_0 = _S3918;
    (&_S4053)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4052, &_S4053, _S4050.differential_0);
    DiffPair_float_0 _S4054;
    (&_S4054)->primal_0 = _S3919;
    (&_S4054)->differential_0 = 0.0f;
    DiffPair_float_0 _S4055;
    (&_S4055)->primal_0 = _S3920;
    (&_S4055)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4054, &_S4055, 0.0f);
    float _S4056 = _S4051.differential_0 + _S4055.differential_0;
    DiffPair_float_0 _S4057;
    (&_S4057)->primal_0 = _S3917;
    (&_S4057)->differential_0 = 0.0f;
    DiffPair_float_0 _S4058;
    (&_S4058)->primal_0 = _S3918;
    (&_S4058)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4057, &_S4058, _S4054.differential_0);
    float _S4059 = _S4053.differential_0 + _S4058.differential_0;
    float _S4060 = _S4052.differential_0 + _S4057.differential_0;
    DiffPair_float_0 _S4061;
    (&_S4061)->primal_0 = _S3916;
    (&_S4061)->differential_0 = 0.0f;
    DiffPair_float_0 _S4062;
    (&_S4062)->primal_0 = _S3915;
    (&_S4062)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4061, &_S4062, 0.0f);
    DiffPair_float_0 _S4063;
    (&_S4063)->primal_0 = _S3912;
    (&_S4063)->differential_0 = 0.0f;
    DiffPair_float_0 _S4064;
    (&_S4064)->primal_0 = _S3913;
    (&_S4064)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4063, &_S4064, _S4061.differential_0);
    DiffPair_float_0 _S4065;
    (&_S4065)->primal_0 = _S3914;
    (&_S4065)->differential_0 = 0.0f;
    DiffPair_float_0 _S4066;
    (&_S4066)->primal_0 = _S3915;
    (&_S4066)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4065, &_S4066, 0.0f);
    float _S4067 = _S4062.differential_0 + _S4066.differential_0;
    DiffPair_float_0 _S4068;
    (&_S4068)->primal_0 = _S3912;
    (&_S4068)->differential_0 = 0.0f;
    DiffPair_float_0 _S4069;
    (&_S4069)->primal_0 = _S3913;
    (&_S4069)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4068, &_S4069, _S4065.differential_0);
    float _S4070 = _S4064.differential_0 + _S4069.differential_0;
    float _S4071 = _S4063.differential_0 + _S4068.differential_0;
    DiffPair_float_0 _S4072;
    (&_S4072)->primal_0 = _S3910;
    (&_S4072)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4072, 0.0f);
    float _S4073 = - (-1.0f * - (_S4072.differential_0 / _S3911));
    float2  _S4074 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4075;
    (&_S4075)->primal_0 = e2_6;
    (&_S4075)->differential_0 = _S4074;
    s_bwd_length_impl_1(&_S4075, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4076;
    (&_S4076)->primal_0 = e1_14;
    (&_S4076)->differential_0 = _S4074;
    s_bwd_length_impl_1(&_S4076, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4077;
    (&_S4077)->primal_0 = e0_14;
    (&_S4077)->differential_0 = _S4074;
    s_bwd_length_impl_1(&_S4077, -0.0f);
    DiffPair_float_0 _S4078;
    (&_S4078)->primal_0 = _S3908;
    (&_S4078)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4078, 0.0f);
    float _S4079 = - _S4078.differential_0;
    float2  _S4080 = _S4076.differential_0 + make_float2 (_S3906 * _S4079, _S3904 * _S4078.differential_0);
    float2  _S4081 = _S4077.differential_0 + make_float2 (_S3905 * _S4078.differential_0, _S3907 * _S4079);
    float2  _S4082 = _S3899 * (- _S4075.differential_0 + _S4080 + make_float2 (_S4067, _S4056)) / _S3898;
    float2  _S4083 = _S3895 * - _S4082;
    float2  _S4084 = _S3897 * _S4082;
    float2  _S4085 = _S3899 * (- _S4080 + _S4081 + make_float2 (_S4070, _S4059)) / _S3894;
    float2  _S4086 = _S3891 * - _S4085;
    float2  _S4087 = _S3893 * _S4085;
    float _S4088 = _S4086.x + _S4086.y;
    float2  _S4089 = _S3899 * (_S4075.differential_0 + - _S4081 + make_float2 (_S4071, _S4060)) / _S3890;
    float2  _S4090 = _S3887 * - _S4089;
    float2  _S4091 = _S3889 * _S4089;
    float _S4092 = _S4090.x + _S4090.y;
    float3  _S4093 = _S3969.differential_0 + _S4049.differential_0 + make_float3 (_S4084.x, _S4084.y, _S4083.x + _S4083.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4094;
    (&_S4094)->primal_0 = R_29;
    (&_S4094)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4095;
    (&_S4095)->primal_0 = vert2_6;
    (&_S4095)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4094, &_S4095, _S4093);
    float3  _S4096 = _S3968.differential_0 + _S4049.differential_0 + make_float3 (_S4087.x, _S4087.y, _S4088);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4097;
    (&_S4097)->primal_0 = R_29;
    (&_S4097)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4098;
    (&_S4098)->primal_0 = vert1_6;
    (&_S4098)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4097, &_S4098, _S4096);
    float3  _S4099 = _S3970 + _S3971 + _S4049.differential_0 + make_float3 (_S4091.x, _S4091.y, _S4092);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4100;
    (&_S4100)->primal_0 = R_29;
    (&_S4100)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4101;
    (&_S4101)->primal_0 = vert0_6;
    (&_S4101)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4100, &_S4101, _S4099);
    float3  _S4102 = (*v_verts_0)[int(2)] + _S4095.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4103;
    (&_S4103)->primal_0 = _S3881;
    (&_S4103)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4104;
    (&_S4104)->primal_0 = _S3886;
    (&_S4104)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4103, &_S4104, _S4102);
    float _S4105 = - _S4104.differential_0.y;
    float _S4106 = _S3885 * _S4104.differential_0.x;
    float _S4107 = - (_S3877 * _S4104.differential_0.x);
    float3  _S4108 = (*v_verts_0)[int(1)] + _S4098.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4109;
    (&_S4109)->primal_0 = _S3881;
    (&_S4109)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4110;
    (&_S4110)->primal_0 = _S3884;
    (&_S4110)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4109, &_S4110, _S4108);
    float _S4111 = _S3877 * _S4110.differential_0.x;
    float _S4112 = _S3883 * _S4110.differential_0.x;
    float3  _S4113 = (*v_verts_0)[int(0)] + _S4101.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4114;
    (&_S4114)->primal_0 = _S3881;
    (&_S4114)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4115;
    (&_S4115)->primal_0 = _S3882;
    (&_S4115)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4114, &_S4115, _S4113);
    Matrix<float, 3, 3>  _S4116 = transpose_0(_S4103.differential_0 + _S4109.differential_0 + _S4114.differential_0);
    float _S4117 = 2.0f * - _S4116.rows[int(2)].z;
    float _S4118 = 2.0f * _S4116.rows[int(2)].y;
    float _S4119 = 2.0f * _S4116.rows[int(2)].x;
    float _S4120 = 2.0f * _S4116.rows[int(1)].z;
    float _S4121 = 2.0f * - _S4116.rows[int(1)].y;
    float _S4122 = 2.0f * _S4116.rows[int(1)].x;
    float _S4123 = 2.0f * _S4116.rows[int(0)].z;
    float _S4124 = 2.0f * _S4116.rows[int(0)].y;
    float _S4125 = 2.0f * - _S4116.rows[int(0)].x;
    float _S4126 = - _S4122 + _S4124;
    float _S4127 = _S4119 + - _S4123;
    float _S4128 = - _S4118 + _S4120;
    float _S4129 = _S4118 + _S4120;
    float _S4130 = _S4119 + _S4123;
    float _S4131 = _S4122 + _S4124;
    float _S4132 = quat_30.w * (_S4121 + _S4125);
    float _S4133 = quat_30.z * (_S4117 + _S4125);
    float _S4134 = quat_30.y * (_S4117 + _S4121);
    float _S4135 = quat_30.x * _S4126 + quat_30.z * _S4129 + quat_30.y * _S4130 + _S4132 + _S4132;
    float _S4136 = quat_30.x * _S4127 + quat_30.w * _S4129 + quat_30.y * _S4131 + _S4133 + _S4133;
    float _S4137 = quat_30.x * _S4128 + quat_30.w * _S4130 + quat_30.z * _S4131 + _S4134 + _S4134;
    float _S4138 = quat_30.w * _S4126 + quat_30.z * _S4127 + quat_30.y * _S4128;
    float _S4139 = _S4107 + _S4111;
    float _S4140 = 0.5f * - _S4139;
    float _S4141 = _S4105 + _S4110.differential_0.y;
    DiffPair_float_0 _S4142;
    (&_S4142)->primal_0 = _S3878;
    (&_S4142)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4142, _S4141);
    float _S4143 = _S4140 + _S4142.differential_0;
    float _S4144 = _S4106 + _S4112 + _S4115.differential_0.x;
    DiffPair_float_0 _S4145;
    (&_S4145)->primal_0 = _S3876;
    (&_S4145)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4145, _S4144);
    float _S4146 = _S4140 + _S4145.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4147;
    (&_S4147)->primal_0 = R_29;
    (&_S4147)->differential_0 = _S4043;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4148;
    (&_S4148)->primal_0 = mean_29;
    (&_S4148)->differential_0 = _S3963;
    s_bwd_prop_mul_1(&_S4147, &_S4148, _S3965.differential_0);
    float3  _S4149 = _S4045.differential_0 + _S4093 + _S4096 + _S4099 + _S3965.differential_0;
    Matrix<float, 3, 3>  _S4150 = _S4046 + _S4094.differential_0 + _S4097.differential_0 + _S4100.differential_0 + _S4147.differential_0;
    FixedArray<float3 , 2>  _S4151;
    _S4151[int(0)] = _S3963;
    _S4151[int(1)] = _S3963;
    _S4151[int(1)] = _S3977;
    _S4151[int(0)] = _S3982;
    FixedArray<float3 , 16>  _S4152;
    _S4152[int(0)] = _S3963;
    _S4152[int(1)] = _S3963;
    _S4152[int(2)] = _S3963;
    _S4152[int(3)] = _S3963;
    _S4152[int(4)] = _S3963;
    _S4152[int(5)] = _S3963;
    _S4152[int(6)] = _S3963;
    _S4152[int(7)] = _S3963;
    _S4152[int(8)] = _S3963;
    _S4152[int(9)] = _S3963;
    _S4152[int(10)] = _S3963;
    _S4152[int(11)] = _S3963;
    _S4152[int(12)] = _S3963;
    _S4152[int(13)] = _S3963;
    _S4152[int(14)] = _S3963;
    _S4152[int(15)] = _S3963;
    _S4152[int(15)] = _S3984;
    _S4152[int(14)] = _S3986;
    _S4152[int(13)] = _S3988;
    _S4152[int(12)] = _S3990;
    _S4152[int(11)] = _S3992;
    _S4152[int(10)] = _S3994;
    _S4152[int(9)] = _S3996;
    _S4152[int(8)] = _S4004;
    _S4152[int(7)] = _S4006;
    _S4152[int(6)] = _S4008;
    _S4152[int(5)] = _S4010;
    _S4152[int(4)] = _S4012;
    _S4152[int(3)] = _S4023;
    _S4152[int(2)] = _S4025;
    _S4152[int(1)] = _S4027;
    _S4152[int(0)] = _S4040;
    float2  _S4153 = make_float2 (0.0f, _S4073);
    float3  _S4154 = make_float3 (_S4146, _S4143, _S4139);
    float4  _S4155 = make_float4 (0.0f);
    *&((&_S4155)->w) = _S4135;
    *&((&_S4155)->z) = _S4136;
    *&((&_S4155)->y) = _S4137;
    *&((&_S4155)->x) = _S4138;
    *v_mean_9 = _S4041 + _S4102 + _S4108 + _S4113 + _S4148.differential_0;
    *v_quat_8 = _S4155;
    *v_scale_8 = _S4154;
    *v_hardness_4 = _S4153;
    *v_sh_coeffs_7 = _S4152;
    *v_ch_coeffs_2 = _S4151;
    *v_R_8 = _S4150;
    *v_t_8 = _S4149;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_30, float3  t_29, float fx_30, float fy_30, float cx_30, float cy_30, float4  radial_coeffs_35, float2  tangential_coeffs_35, float2  thin_prism_coeffs_35, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_1(R_30, mean_30) + t_29;
    float _S4156 = scale_30.x;
    float _S4157 = s_primal_ctx_exp_1(_S4156);
    float _S4158 = scale_30.y;
    float _S4159 = s_primal_ctx_exp_1(_S4158);
    float sz_11 = scale_30.z - 0.5f * (_S4156 + _S4158);
    float _S4160 = quat_31.y;
    float x2_31 = _S4160 * _S4160;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_57 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4161 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_57), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_57), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4162 = make_float3 (_S4157, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_1(_S4161, _S4162) + mean_30;
    float _S4163 = -0.5f + sz_11;
    float3  _S4164 = make_float3 (_S4157 * _S4163, _S4159, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_1(_S4161, _S4164) + mean_30;
    float _S4165 = -0.5f - sz_11;
    float3  _S4166 = make_float3 (_S4157 * _S4165, - _S4159, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_1(_S4161, _S4166) + mean_30;
    float3  vert0_c_11 = s_primal_ctx_mul_1(R_30, vert0_7) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_1(R_30, vert1_7) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_1(R_30, vert2_7) + t_29;
    CameraDistortion_0 _S4167 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_35, tangential_coeffs_35, thin_prism_coeffs_35);
    float2  _S4168 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4169 = length_0(_S4168);
    float _S4170 = vert0_c_11.z;
    float _S4171 = s_primal_ctx_atan2_0(_S4169, _S4170);
    bool _S4172 = _S4171 < 0.00100000004749745f;
    float k_9;
    float _S4173;
    float _S4174;
    float _S4175;
    if(_S4172)
    {
        float _S4176 = 1.0f - _S4171 * _S4171 / 3.0f;
        float _S4177 = _S4170 * _S4170;
        k_9 = _S4176 / _S4170;
        _S4173 = 0.0f;
        _S4174 = _S4177;
        _S4175 = _S4176;
    }
    else
    {
        float _S4178 = _S4169 * _S4169;
        k_9 = _S4171 / _S4169;
        _S4173 = _S4178;
        _S4174 = 0.0f;
        _S4175 = 0.0f;
    }
    float2  _S4179 = make_float2 (k_9);
    float2  _S4180 = _S4168 * make_float2 (k_9);
    float k1_8 = _S4167.radial_coeffs_0.x;
    float k2_8 = _S4167.radial_coeffs_0.y;
    float k3_8 = _S4167.radial_coeffs_0.z;
    float k4_9 = _S4167.radial_coeffs_0.w;
    float p1_10 = _S4167.tangential_coeffs_0.x;
    float p2_10 = _S4167.tangential_coeffs_0.y;
    float sx1_9 = _S4167.thin_prism_coeffs_0.x;
    float sy1_9 = _S4167.thin_prism_coeffs_0.y;
    float u_26 = _S4180.x;
    float v_26 = _S4180.y;
    float r2_26 = u_26 * u_26 + v_26 * v_26;
    float _S4181 = k3_8 + r2_26 * k4_9;
    float _S4182 = k2_8 + r2_26 * _S4181;
    float _S4183 = k1_8 + r2_26 * _S4182;
    float radial_5 = 1.0f + r2_26 * _S4183;
    float2  _S4184 = make_float2 (radial_5);
    float _S4185 = 2.0f * p1_10;
    float _S4186 = _S4185 * u_26;
    float _S4187 = 2.0f * u_26;
    float _S4188 = r2_26 + _S4187 * u_26;
    float _S4189 = 2.0f * p2_10;
    float _S4190 = _S4189 * u_26;
    float _S4191 = 2.0f * v_26;
    float _S4192 = r2_26 + _S4191 * v_26;
    float2  _S4193 = _S4180 * make_float2 (radial_5) + make_float2 (_S4186 * v_26 + p2_10 * _S4188 + sx1_9 * r2_26, _S4190 * v_26 + p1_10 * _S4192 + sy1_9 * r2_26);
    float _S4194 = fx_30 * _S4193.x + cx_30;
    float _S4195 = fy_30 * _S4193.y + cy_30;
    float2  _S4196 = make_float2 (_S4194, _S4195);
    float2  _S4197 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4198 = length_0(_S4197);
    float _S4199 = vert1_c_11.z;
    float _S4200 = s_primal_ctx_atan2_0(_S4198, _S4199);
    bool _S4201 = _S4200 < 0.00100000004749745f;
    float _S4202;
    float _S4203;
    float _S4204;
    if(_S4201)
    {
        float _S4205 = 1.0f - _S4200 * _S4200 / 3.0f;
        float _S4206 = _S4199 * _S4199;
        k_9 = _S4205 / _S4199;
        _S4202 = 0.0f;
        _S4203 = _S4206;
        _S4204 = _S4205;
    }
    else
    {
        float _S4207 = _S4198 * _S4198;
        k_9 = _S4200 / _S4198;
        _S4202 = _S4207;
        _S4203 = 0.0f;
        _S4204 = 0.0f;
    }
    float2  _S4208 = make_float2 (k_9);
    float2  _S4209 = _S4197 * make_float2 (k_9);
    float u_27 = _S4209.x;
    float v_27 = _S4209.y;
    float r2_27 = u_27 * u_27 + v_27 * v_27;
    float _S4210 = k3_8 + r2_27 * k4_9;
    float _S4211 = k2_8 + r2_27 * _S4210;
    float _S4212 = k1_8 + r2_27 * _S4211;
    float radial_6 = 1.0f + r2_27 * _S4212;
    float2  _S4213 = make_float2 (radial_6);
    float _S4214 = _S4185 * u_27;
    float _S4215 = 2.0f * u_27;
    float _S4216 = r2_27 + _S4215 * u_27;
    float _S4217 = _S4189 * u_27;
    float _S4218 = 2.0f * v_27;
    float _S4219 = r2_27 + _S4218 * v_27;
    float2  _S4220 = _S4209 * make_float2 (radial_6) + make_float2 (_S4214 * v_27 + p2_10 * _S4216 + sx1_9 * r2_27, _S4217 * v_27 + p1_10 * _S4219 + sy1_9 * r2_27);
    float _S4221 = fx_30 * _S4220.x + cx_30;
    float _S4222 = fy_30 * _S4220.y + cy_30;
    float2  _S4223 = make_float2 (_S4221, _S4222);
    float2  _S4224 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4225 = length_0(_S4224);
    float _S4226 = vert2_c_11.z;
    float _S4227 = s_primal_ctx_atan2_0(_S4225, _S4226);
    bool _S4228 = _S4227 < 0.00100000004749745f;
    float _S4229;
    float _S4230;
    float _S4231;
    if(_S4228)
    {
        float _S4232 = 1.0f - _S4227 * _S4227 / 3.0f;
        float _S4233 = _S4226 * _S4226;
        k_9 = _S4232 / _S4226;
        _S4229 = 0.0f;
        _S4230 = _S4233;
        _S4231 = _S4232;
    }
    else
    {
        float _S4234 = _S4225 * _S4225;
        k_9 = _S4227 / _S4225;
        _S4229 = _S4234;
        _S4230 = 0.0f;
        _S4231 = 0.0f;
    }
    float2  _S4235 = make_float2 (k_9);
    float2  _S4236 = _S4224 * make_float2 (k_9);
    float u_28 = _S4236.x;
    float v_28 = _S4236.y;
    float r2_28 = u_28 * u_28 + v_28 * v_28;
    float _S4237 = k3_8 + r2_28 * k4_9;
    float _S4238 = k2_8 + r2_28 * _S4237;
    float _S4239 = k1_8 + r2_28 * _S4238;
    float radial_7 = 1.0f + r2_28 * _S4239;
    float2  _S4240 = make_float2 (radial_7);
    float _S4241 = _S4185 * u_28;
    float _S4242 = 2.0f * u_28;
    float _S4243 = r2_28 + _S4242 * u_28;
    float _S4244 = _S4189 * u_28;
    float _S4245 = 2.0f * v_28;
    float _S4246 = r2_28 + _S4245 * v_28;
    float2  _S4247 = _S4236 * make_float2 (radial_7) + make_float2 (_S4241 * v_28 + p2_10 * _S4243 + sx1_9 * r2_28, _S4244 * v_28 + p1_10 * _S4246 + sy1_9 * r2_28);
    float _S4248 = fx_30 * _S4247.x + cx_30;
    float _S4249 = fy_30 * _S4247.y + cy_30;
    float2  _S4250 = make_float2 (_S4248, _S4249);
    float2  e0_15 = _S4223 - _S4196;
    float2  e1_15 = _S4250 - _S4223;
    float2  e2_7 = _S4196 - _S4250;
    float _S4251 = e0_15.x;
    float _S4252 = e1_15.y;
    float _S4253 = e0_15.y;
    float _S4254 = e1_15.x;
    float _S4255 = _S4251 * _S4252 - _S4253 * _S4254;
    float _S4256 = 1.0f - hardness_15.y;
    float _S4257 = -1.0f / _S4256;
    float _S4258 = _S4256 * _S4256;
    float _S4259 = s_primal_ctx_max_0(_S4194, _S4221);
    float _S4260 = s_primal_ctx_min_0(_S4194, _S4221);
    float _S4261 = s_primal_ctx_max_0(_S4195, _S4222);
    float _S4262 = s_primal_ctx_min_0(_S4195, _S4222);
    float3  _S4263 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S4264 = length_1(_S4263) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4265 = transpose_0(R_30);
    float3  _S4266 = mean_30 - - s_primal_ctx_mul_1(_S4265, t_29);
    float _S4267 = _S4266.x;
    float _S4268 = _S4266.y;
    float _S4269 = _S4266.z;
    float _S4270 = _S4267 * _S4267 + _S4268 * _S4268 + _S4269 * _S4269;
    float _S4271 = s_primal_ctx_sqrt_0(_S4270);
    float x_61 = _S4267 / _S4271;
    float3  _S4272 = make_float3 (x_61);
    float _S4273 = _S4271 * _S4271;
    float y_29 = _S4268 / _S4271;
    float z_26 = _S4269 / _S4271;
    float3  _S4274 = make_float3 (z_26);
    float _S4275 = - y_29;
    float3  _S4276 = make_float3 (_S4275);
    float z2_58 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_61 * x_61 - y_29 * y_29;
    float _S4277 = 2.0f * x_61;
    float fS1_26 = _S4277 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_58 - 0.31539157032966614f;
    float3  _S4278 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_61;
    float3  _S4279 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4280 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4281 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4282 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_58 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4283 = 1.86588168144226074f * z2_58 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4283;
    float3  _S4284 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_61;
    float3  _S4285 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4286 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4287 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4288 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_61 * fC1_26 - y_29 * fS1_26);
    float3  _S4289 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_61 * fS1_26 + y_29 * fC1_26);
    float3  _S4290 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4275) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_61) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4291 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4292 = make_float3 (0.0f);
    float3  _S4293 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4294 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4295 = make_float3 (_S4294);
    float3  _S4296 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4294);
    float3  _S4297 = _S4293 + _S4296 + make_float3 (0.5f);
    float3  _S4298 = _S4293 - _S4296 + make_float3 (0.5f);
    float3  _S4299 = vert1_c_11 - vert0_c_11;
    float3  _S4300 = vert2_c_11 - vert0_c_11;
    float3  _S4301 = s_primal_ctx_cross_0(_S4299, _S4300);
    float3  _S4302 = normalize_0(_S4301);
    float3  _S4303 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4302, mean_c_26)))))) * v_normal_3;
    float3  _S4304 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4305;
    (&_S4305)->primal_0 = _S4302;
    (&_S4305)->differential_0 = _S4304;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4306;
    (&_S4306)->primal_0 = mean_c_26;
    (&_S4306)->differential_0 = _S4304;
    s_bwd_prop_dot_0(&_S4305, &_S4306, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4307 = _S4306;
    float3  _S4308 = _S4303 + _S4305.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4309;
    (&_S4309)->primal_0 = _S4301;
    (&_S4309)->differential_0 = _S4304;
    s_bwd_normalize_impl_0(&_S4309, _S4308);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4310;
    (&_S4310)->primal_0 = _S4299;
    (&_S4310)->differential_0 = _S4304;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4311;
    (&_S4311)->primal_0 = _S4300;
    (&_S4311)->differential_0 = _S4304;
    s_bwd_prop_cross_0(&_S4310, &_S4311, _S4309.differential_0);
    float3  _S4312 = - _S4311.differential_0;
    float3  _S4313 = - _S4310.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4314;
    (&_S4314)->primal_0 = _S4298;
    (&_S4314)->differential_0 = _S4304;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4315;
    (&_S4315)->primal_0 = _S4292;
    (&_S4315)->differential_0 = _S4304;
    s_bwd_prop_max_0(&_S4314, &_S4315, (*v_rgbs_1)[int(2)]);
    float3  _S4316 = - _S4314.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4317;
    (&_S4317)->primal_0 = _S4297;
    (&_S4317)->differential_0 = _S4304;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4318;
    (&_S4318)->primal_0 = _S4292;
    (&_S4318)->differential_0 = _S4304;
    s_bwd_prop_max_0(&_S4317, &_S4318, (*v_rgbs_1)[int(1)]);
    float3  _S4319 = _S4295 * (_S4316 + _S4317.differential_0);
    float3  _S4320 = _S4314.differential_0 + _S4317.differential_0;
    float3  _S4321 = make_float3 (0.5f) * - _S4320;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4322;
    (&_S4322)->primal_0 = _S4291;
    (&_S4322)->differential_0 = _S4304;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4323;
    (&_S4323)->primal_0 = _S4292;
    (&_S4323)->differential_0 = _S4304;
    s_bwd_prop_max_0(&_S4322, &_S4323, (*v_rgbs_1)[int(0)]);
    float3  _S4324 = _S4321 + _S4322.differential_0;
    float3  _S4325 = _S4320 + _S4322.differential_0;
    float3  _S4326 = _S4289 * _S4325;
    float3  _S4327 = (*sh_coeffs_26)[int(15)] * _S4325;
    float3  _S4328 = _S4287 * _S4325;
    float3  _S4329 = (*sh_coeffs_26)[int(14)] * _S4325;
    float3  _S4330 = _S4285 * _S4325;
    float3  _S4331 = (*sh_coeffs_26)[int(13)] * _S4325;
    float3  _S4332 = _S4284 * _S4325;
    float3  _S4333 = (*sh_coeffs_26)[int(12)] * _S4325;
    float3  _S4334 = _S4286 * _S4325;
    float3  _S4335 = (*sh_coeffs_26)[int(11)] * _S4325;
    float3  _S4336 = _S4288 * _S4325;
    float3  _S4337 = (*sh_coeffs_26)[int(10)] * _S4325;
    float3  _S4338 = _S4290 * _S4325;
    float3  _S4339 = (*sh_coeffs_26)[int(9)] * _S4325;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4339.x + _S4339.y + _S4339.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4327.x + _S4327.y + _S4327.z);
    float _S4340 = _S4337.x + _S4337.y + _S4337.z;
    float _S4341 = _S4329.x + _S4329.y + _S4329.z;
    float _S4342 = _S4335.x + _S4335.y + _S4335.z;
    float _S4343 = _S4331.x + _S4331.y + _S4331.z;
    float _S4344 = _S4333.x + _S4333.y + _S4333.z;
    float _S4345 = - s_diff_fC2_T_8;
    float3  _S4346 = _S4281 * _S4325;
    float3  _S4347 = (*sh_coeffs_26)[int(8)] * _S4325;
    float3  _S4348 = _S4279 * _S4325;
    float3  _S4349 = (*sh_coeffs_26)[int(7)] * _S4325;
    float3  _S4350 = _S4278 * _S4325;
    float3  _S4351 = (*sh_coeffs_26)[int(6)] * _S4325;
    float3  _S4352 = _S4280 * _S4325;
    float3  _S4353 = (*sh_coeffs_26)[int(5)] * _S4325;
    float3  _S4354 = _S4282 * _S4325;
    float3  _S4355 = (*sh_coeffs_26)[int(4)] * _S4325;
    float _S4356 = _S4353.x + _S4353.y + _S4353.z;
    float _S4357 = _S4349.x + _S4349.y + _S4349.z;
    float _S4358 = fTmp1B_26 * _S4340 + x_61 * s_diff_fS2_T_8 + y_29 * _S4345 + 0.54627424478530884f * (_S4355.x + _S4355.y + _S4355.z);
    float _S4359 = fTmp1B_26 * _S4341 + y_29 * s_diff_fS2_T_8 + x_61 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S4347.x + _S4347.y + _S4347.z);
    float _S4360 = y_29 * - _S4359;
    float _S4361 = x_61 * _S4359;
    float _S4362 = z_26 * (1.86588168144226074f * (z_26 * _S4344) + -2.28522896766662598f * (y_29 * _S4342 + x_61 * _S4343) + 0.94617468118667603f * (_S4351.x + _S4351.y + _S4351.z));
    float3  _S4363 = make_float3 (0.48860251903533936f) * _S4325;
    float3  _S4364 = - _S4363;
    float3  _S4365 = _S4272 * _S4364;
    float3  _S4366 = (*sh_coeffs_26)[int(3)] * _S4364;
    float3  _S4367 = _S4274 * _S4363;
    float3  _S4368 = (*sh_coeffs_26)[int(2)] * _S4363;
    float3  _S4369 = _S4276 * _S4363;
    float3  _S4370 = (*sh_coeffs_26)[int(1)] * _S4363;
    float _S4371 = (_S4283 * _S4344 + 1.44530570507049561f * (fS1_26 * _S4340 + fC1_26 * _S4341) + -1.09254848957061768f * (y_29 * _S4356 + x_61 * _S4357) + _S4362 + _S4362 + _S4368.x + _S4368.y + _S4368.z) / _S4273;
    float _S4372 = _S4271 * _S4371;
    float _S4373 = (fTmp0C_26 * _S4342 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S4345 + fTmp0B_26 * _S4356 + _S4277 * _S4358 + _S4360 + _S4360 + - (_S4370.x + _S4370.y + _S4370.z)) / _S4273;
    float _S4374 = _S4271 * _S4373;
    float _S4375 = (fTmp0C_26 * _S4343 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S4357 + 2.0f * (y_29 * _S4358) + _S4361 + _S4361 + _S4366.x + _S4366.y + _S4366.z) / _S4273;
    float _S4376 = _S4271 * _S4375;
    float _S4377 = _S4269 * - _S4371 + _S4268 * - _S4373 + _S4267 * - _S4375;
    DiffPair_float_0 _S4378;
    (&_S4378)->primal_0 = _S4270;
    (&_S4378)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4378, _S4377);
    float _S4379 = _S4269 * _S4378.differential_0;
    float _S4380 = _S4268 * _S4378.differential_0;
    float _S4381 = _S4267 * _S4378.differential_0;
    float3  _S4382 = make_float3 (0.282094806432724f) * _S4325;
    float3  _S4383 = make_float3 (_S4376 + _S4381 + _S4381, _S4374 + _S4380 + _S4380, _S4372 + _S4379 + _S4379);
    float3  _S4384 = - - _S4383;
    Matrix<float, 3, 3>  _S4385 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4386;
    (&_S4386)->primal_0 = _S4265;
    (&_S4386)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4387;
    (&_S4387)->primal_0 = t_29;
    (&_S4387)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4386, &_S4387, _S4384);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4388 = _S4387;
    Matrix<float, 3, 3>  _S4389 = transpose_0(_S4386.differential_0);
    DiffPair_float_0 _S4390;
    (&_S4390)->primal_0 = _S4264;
    (&_S4390)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4390, v_depth_10);
    float _S4391 = 0.3333333432674408f * _S4390.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4392;
    (&_S4392)->primal_0 = _S4263;
    (&_S4392)->differential_0 = _S4304;
    s_bwd_length_impl_0(&_S4392, _S4391);
    DiffPair_float_0 _S4393;
    (&_S4393)->primal_0 = _S4262;
    (&_S4393)->differential_0 = 0.0f;
    DiffPair_float_0 _S4394;
    (&_S4394)->primal_0 = _S4249;
    (&_S4394)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4393, &_S4394, 0.0f);
    DiffPair_float_0 _S4395;
    (&_S4395)->primal_0 = _S4195;
    (&_S4395)->differential_0 = 0.0f;
    DiffPair_float_0 _S4396;
    (&_S4396)->primal_0 = _S4222;
    (&_S4396)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4395, &_S4396, _S4393.differential_0);
    DiffPair_float_0 _S4397;
    (&_S4397)->primal_0 = _S4261;
    (&_S4397)->differential_0 = 0.0f;
    DiffPair_float_0 _S4398;
    (&_S4398)->primal_0 = _S4249;
    (&_S4398)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4397, &_S4398, 0.0f);
    DiffPair_float_0 _S4399;
    (&_S4399)->primal_0 = _S4195;
    (&_S4399)->differential_0 = 0.0f;
    DiffPair_float_0 _S4400;
    (&_S4400)->primal_0 = _S4222;
    (&_S4400)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4399, &_S4400, _S4397.differential_0);
    DiffPair_float_0 _S4401;
    (&_S4401)->primal_0 = _S4260;
    (&_S4401)->differential_0 = 0.0f;
    DiffPair_float_0 _S4402;
    (&_S4402)->primal_0 = _S4248;
    (&_S4402)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4401, &_S4402, 0.0f);
    DiffPair_float_0 _S4403;
    (&_S4403)->primal_0 = _S4194;
    (&_S4403)->differential_0 = 0.0f;
    DiffPair_float_0 _S4404;
    (&_S4404)->primal_0 = _S4221;
    (&_S4404)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4403, &_S4404, _S4401.differential_0);
    DiffPair_float_0 _S4405;
    (&_S4405)->primal_0 = _S4259;
    (&_S4405)->differential_0 = 0.0f;
    DiffPair_float_0 _S4406;
    (&_S4406)->primal_0 = _S4248;
    (&_S4406)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4405, &_S4406, 0.0f);
    DiffPair_float_0 _S4407;
    (&_S4407)->primal_0 = _S4194;
    (&_S4407)->differential_0 = 0.0f;
    DiffPair_float_0 _S4408;
    (&_S4408)->primal_0 = _S4221;
    (&_S4408)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4407, &_S4408, _S4405.differential_0);
    DiffPair_float_0 _S4409;
    (&_S4409)->primal_0 = _S4257;
    (&_S4409)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4409, 0.0f);
    float _S4410 = - (-1.0f * - (_S4409.differential_0 / _S4258));
    float2  _S4411 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4412;
    (&_S4412)->primal_0 = e2_7;
    (&_S4412)->differential_0 = _S4411;
    s_bwd_length_impl_1(&_S4412, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4413;
    (&_S4413)->primal_0 = e1_15;
    (&_S4413)->differential_0 = _S4411;
    s_bwd_length_impl_1(&_S4413, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4414;
    (&_S4414)->primal_0 = e0_15;
    (&_S4414)->differential_0 = _S4411;
    s_bwd_length_impl_1(&_S4414, -0.0f);
    DiffPair_float_0 _S4415;
    (&_S4415)->primal_0 = _S4255;
    (&_S4415)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4415, 0.0f);
    float _S4416 = - _S4415.differential_0;
    float2  _S4417 = _S4413.differential_0 + make_float2 (_S4253 * _S4416, _S4251 * _S4415.differential_0);
    float2  _S4418 = - _S4417;
    float2  _S4419 = _S4414.differential_0 + make_float2 (_S4252 * _S4415.differential_0, _S4254 * _S4416);
    float2  _S4420 = - _S4419;
    float2  _S4421 = - _S4412.differential_0 + _S4417;
    float _S4422 = fy_30 * (_S4394.differential_0 + _S4398.differential_0 + _S4421.y);
    float _S4423 = fx_30 * (_S4402.differential_0 + _S4406.differential_0 + _S4421.x);
    float2  _S4424 = make_float2 (_S4423, _S4422);
    float2  _S4425 = _S4236 * _S4424;
    float2  _S4426 = _S4240 * _S4424;
    float _S4427 = r2_28 * _S4422;
    float _S4428 = p1_10 * _S4422;
    float _S4429 = _S4246 * _S4422;
    float _S4430 = v_28 * _S4422;
    float _S4431 = u_28 * _S4430;
    float _S4432 = r2_28 * _S4423;
    float _S4433 = p2_10 * _S4423;
    float _S4434 = _S4243 * _S4423;
    float _S4435 = v_28 * _S4423;
    float _S4436 = u_28 * _S4435;
    float _S4437 = _S4425.x + _S4425.y;
    float _S4438 = r2_28 * _S4437;
    float _S4439 = r2_28 * _S4438;
    float _S4440 = r2_28 * _S4439;
    float _S4441 = r2_28 * _S4440;
    float _S4442 = sy1_9 * _S4422 + _S4428 + sx1_9 * _S4423 + _S4433 + _S4239 * _S4437 + _S4238 * _S4438 + _S4237 * _S4439 + k4_9 * _S4440;
    float _S4443 = v_28 * _S4442;
    float _S4444 = u_28 * _S4442;
    float _S4445 = _S4245 * _S4428 + 2.0f * (v_28 * _S4428) + _S4244 * _S4422 + _S4241 * _S4423 + _S4443 + _S4443;
    float _S4446 = _S4189 * _S4430 + _S4242 * _S4433 + 2.0f * (u_28 * _S4433) + _S4185 * _S4435 + _S4444 + _S4444;
    float3  _S4447 = _S4310.differential_0 + _S4392.differential_0;
    float3  _S4448 = _S4312 + _S4313 + _S4392.differential_0;
    float3  _S4449 = _S4311.differential_0 + _S4392.differential_0;
    FixedArray<float3 , 2>  _S4450;
    _S4450[int(0)] = _S4304;
    _S4450[int(1)] = _S4304;
    _S4450[int(1)] = _S4319;
    _S4450[int(0)] = _S4324;
    float3  _S4451 = _S4450[int(0)];
    float3  _S4452 = _S4450[int(1)];
    FixedArray<float3 , 16>  _S4453;
    _S4453[int(0)] = _S4304;
    _S4453[int(1)] = _S4304;
    _S4453[int(2)] = _S4304;
    _S4453[int(3)] = _S4304;
    _S4453[int(4)] = _S4304;
    _S4453[int(5)] = _S4304;
    _S4453[int(6)] = _S4304;
    _S4453[int(7)] = _S4304;
    _S4453[int(8)] = _S4304;
    _S4453[int(9)] = _S4304;
    _S4453[int(10)] = _S4304;
    _S4453[int(11)] = _S4304;
    _S4453[int(12)] = _S4304;
    _S4453[int(13)] = _S4304;
    _S4453[int(14)] = _S4304;
    _S4453[int(15)] = _S4304;
    _S4453[int(7)] = _S4348;
    _S4453[int(0)] = _S4382;
    _S4453[int(1)] = _S4369;
    _S4453[int(2)] = _S4367;
    _S4453[int(3)] = _S4365;
    _S4453[int(4)] = _S4354;
    _S4453[int(5)] = _S4352;
    _S4453[int(6)] = _S4350;
    _S4453[int(15)] = _S4326;
    _S4453[int(8)] = _S4346;
    _S4453[int(9)] = _S4338;
    _S4453[int(10)] = _S4336;
    _S4453[int(11)] = _S4334;
    _S4453[int(12)] = _S4332;
    _S4453[int(13)] = _S4330;
    _S4453[int(14)] = _S4328;
    float3  _S4454 = _S4453[int(0)];
    float3  _S4455 = _S4453[int(1)];
    float3  _S4456 = _S4453[int(2)];
    float3  _S4457 = _S4453[int(3)];
    float3  _S4458 = _S4453[int(4)];
    float3  _S4459 = _S4453[int(5)];
    float3  _S4460 = _S4453[int(6)];
    float3  _S4461 = _S4453[int(7)];
    float3  _S4462 = _S4453[int(8)];
    float3  _S4463 = _S4453[int(9)];
    float3  _S4464 = _S4453[int(10)];
    float3  _S4465 = _S4453[int(11)];
    float3  _S4466 = _S4453[int(12)];
    float3  _S4467 = _S4453[int(13)];
    float3  _S4468 = _S4453[int(14)];
    float3  _S4469 = _S4453[int(15)];
    float _S4470 = _S4403.differential_0 + _S4407.differential_0;
    float2  _S4471 = _S4412.differential_0 + _S4420;
    float _S4472 = _S4396.differential_0 + _S4400.differential_0;
    float _S4473 = _S4395.differential_0 + _S4399.differential_0;
    float2  _S4474 = _S4418 + _S4419;
    float _S4475 = _S4404.differential_0 + _S4408.differential_0;
    float2  _S4476 = make_float2 (0.0f, _S4410);
    float2  _S4477 = _S4426 + make_float2 (_S4446, _S4445);
    float2  _S4478 = _S4224 * _S4477;
    float2  _S4479 = _S4235 * _S4477;
    float _S4480 = _S4478.x + _S4478.y;
    if(_S4228)
    {
        float _S4481 = _S4480 / _S4230;
        float _S4482 = _S4231 * - _S4481;
        float _S4483 = _S4227 * (0.3333333432674408f * - (_S4226 * _S4481));
        k_9 = _S4483 + _S4483;
        _S4229 = _S4482;
        _S4230 = 0.0f;
    }
    else
    {
        float _S4484 = _S4480 / _S4229;
        float _S4485 = _S4227 * - _S4484;
        k_9 = _S4225 * _S4484;
        _S4229 = 0.0f;
        _S4230 = _S4485;
    }
    DiffPair_float_0 _S4486;
    (&_S4486)->primal_0 = _S4225;
    (&_S4486)->differential_0 = 0.0f;
    DiffPair_float_0 _S4487;
    (&_S4487)->primal_0 = _S4226;
    (&_S4487)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4486, &_S4487, k_9);
    float _S4488 = _S4487.differential_0 + _S4229;
    float _S4489 = _S4486.differential_0 + _S4230;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4490;
    (&_S4490)->primal_0 = _S4224;
    (&_S4490)->differential_0 = _S4411;
    s_bwd_length_impl_1(&_S4490, _S4489);
    float2  _S4491 = _S4490.differential_0 + _S4479;
    float _S4492 = fy_30 * (_S4474.y + _S4472);
    float _S4493 = fx_30 * (_S4474.x + _S4475);
    float2  _S4494 = make_float2 (_S4493, _S4492);
    float2  _S4495 = _S4209 * _S4494;
    float _S4496 = p1_10 * _S4492;
    float _S4497 = v_27 * _S4492;
    float _S4498 = p2_10 * _S4493;
    float _S4499 = v_27 * _S4493;
    float _S4500 = _S4495.x + _S4495.y;
    float _S4501 = r2_27 * _S4500;
    float _S4502 = r2_27 * _S4501;
    float _S4503 = r2_27 * _S4502;
    float _S4504 = sy1_9 * _S4492 + _S4496 + sx1_9 * _S4493 + _S4498 + _S4212 * _S4500 + _S4211 * _S4501 + _S4210 * _S4502 + k4_9 * _S4503;
    float _S4505 = v_27 * _S4504;
    float _S4506 = u_27 * _S4504;
    float3  _S4507 = _S4449 + make_float3 (_S4491.x, _S4491.y, _S4488);
    float2  _S4508 = _S4213 * _S4494 + make_float2 (_S4189 * _S4497 + _S4215 * _S4498 + 2.0f * (u_27 * _S4498) + _S4185 * _S4499 + _S4506 + _S4506, _S4218 * _S4496 + 2.0f * (v_27 * _S4496) + _S4217 * _S4492 + _S4214 * _S4493 + _S4505 + _S4505);
    float _S4509 = u_27 * _S4497 + _S4431;
    float _S4510 = u_27 * _S4499 + _S4436;
    float _S4511 = r2_27 * _S4492 + _S4427;
    float _S4512 = r2_27 * _S4503 + _S4441;
    float _S4513 = _S4503 + _S4440;
    float _S4514 = _S4219 * _S4492 + _S4429;
    float _S4515 = _S4502 + _S4439;
    float _S4516 = r2_27 * _S4493 + _S4432;
    float _S4517 = _S4216 * _S4493 + _S4434;
    float _S4518 = _S4501 + _S4438;
    float2  _S4519 = _S4197 * _S4508;
    float2  _S4520 = _S4208 * _S4508;
    float _S4521 = _S4519.x + _S4519.y;
    if(_S4201)
    {
        float _S4522 = _S4521 / _S4203;
        float _S4523 = _S4204 * - _S4522;
        float _S4524 = _S4200 * (0.3333333432674408f * - (_S4199 * _S4522));
        k_9 = _S4524 + _S4524;
        _S4202 = _S4523;
        _S4203 = 0.0f;
    }
    else
    {
        float _S4525 = _S4521 / _S4202;
        float _S4526 = _S4200 * - _S4525;
        k_9 = _S4198 * _S4525;
        _S4202 = 0.0f;
        _S4203 = _S4526;
    }
    DiffPair_float_0 _S4527;
    (&_S4527)->primal_0 = _S4198;
    (&_S4527)->differential_0 = 0.0f;
    DiffPair_float_0 _S4528;
    (&_S4528)->primal_0 = _S4199;
    (&_S4528)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4527, &_S4528, k_9);
    float _S4529 = _S4528.differential_0 + _S4202;
    float _S4530 = _S4527.differential_0 + _S4203;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4531;
    (&_S4531)->primal_0 = _S4197;
    (&_S4531)->differential_0 = _S4411;
    s_bwd_length_impl_1(&_S4531, _S4530);
    float2  _S4532 = _S4531.differential_0 + _S4520;
    float _S4533 = fy_30 * (_S4471.y + _S4473);
    float _S4534 = fx_30 * (_S4471.x + _S4470);
    float2  _S4535 = make_float2 (_S4534, _S4533);
    float2  _S4536 = _S4180 * _S4535;
    float _S4537 = p1_10 * _S4533;
    float _S4538 = v_26 * _S4533;
    float _S4539 = p2_10 * _S4534;
    float _S4540 = v_26 * _S4534;
    float _S4541 = _S4536.x + _S4536.y;
    float _S4542 = r2_26 * _S4541;
    float _S4543 = r2_26 * _S4542;
    float _S4544 = r2_26 * _S4543;
    float _S4545 = sy1_9 * _S4533 + _S4537 + sx1_9 * _S4534 + _S4539 + _S4183 * _S4541 + _S4182 * _S4542 + _S4181 * _S4543 + k4_9 * _S4544;
    float _S4546 = v_26 * _S4545;
    float _S4547 = u_26 * _S4545;
    float2  _S4548 = make_float2 (r2_26 * _S4534 + _S4516, r2_26 * _S4533 + _S4511);
    float2  _S4549 = make_float2 (_S4192 * _S4533 + 2.0f * (u_26 * _S4540 + _S4510) + _S4514, 2.0f * (u_26 * _S4538 + _S4509) + _S4188 * _S4534 + _S4517);
    float4  _S4550 = make_float4 (_S4542 + _S4518, _S4543 + _S4515, _S4544 + _S4513, r2_26 * _S4544 + _S4512);
    float3  _S4551 = _S4447 + make_float3 (_S4532.x, _S4532.y, _S4529);
    float2  _S4552 = _S4184 * _S4535 + make_float2 (_S4189 * _S4538 + _S4187 * _S4539 + 2.0f * (u_26 * _S4539) + _S4185 * _S4540 + _S4547 + _S4547, _S4191 * _S4537 + 2.0f * (v_26 * _S4537) + _S4190 * _S4533 + _S4186 * _S4534 + _S4546 + _S4546);
    CameraDistortion_0 _S4553 = CameraDistortion_x24_syn_dzero_0();
    (&_S4553)->thin_prism_coeffs_0 = _S4548;
    (&_S4553)->tangential_coeffs_0 = _S4549;
    (&_S4553)->radial_coeffs_0 = _S4550;
    float2  _S4554 = _S4168 * _S4552;
    float2  _S4555 = _S4179 * _S4552;
    float _S4556 = _S4554.x + _S4554.y;
    if(_S4172)
    {
        float _S4557 = _S4556 / _S4174;
        float _S4558 = _S4175 * - _S4557;
        float _S4559 = _S4171 * (0.3333333432674408f * - (_S4170 * _S4557));
        k_9 = _S4559 + _S4559;
        _S4173 = _S4558;
        _S4174 = 0.0f;
    }
    else
    {
        float _S4560 = _S4556 / _S4173;
        float _S4561 = _S4171 * - _S4560;
        k_9 = _S4169 * _S4560;
        _S4173 = 0.0f;
        _S4174 = _S4561;
    }
    DiffPair_float_0 _S4562;
    (&_S4562)->primal_0 = _S4169;
    (&_S4562)->differential_0 = 0.0f;
    DiffPair_float_0 _S4563;
    (&_S4563)->primal_0 = _S4170;
    (&_S4563)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4562, &_S4563, k_9);
    float _S4564 = _S4563.differential_0 + _S4173;
    float _S4565 = _S4562.differential_0 + _S4174;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4566;
    (&_S4566)->primal_0 = _S4168;
    (&_S4566)->differential_0 = _S4411;
    s_bwd_length_impl_1(&_S4566, _S4565);
    float2  _S4567 = _S4566.differential_0 + _S4555;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4568;
    (&_S4568)->primal_0 = vert2_c_11;
    (&_S4568)->differential_0 = _S4304;
    s_bwd_length_impl_0(&_S4568, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4569;
    (&_S4569)->primal_0 = vert1_c_11;
    (&_S4569)->differential_0 = _S4304;
    s_bwd_length_impl_0(&_S4569, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4570;
    (&_S4570)->primal_0 = vert0_c_11;
    (&_S4570)->differential_0 = _S4304;
    s_bwd_length_impl_0(&_S4570, 0.0f);
    float3  _S4571 = _S4568.differential_0 + _S4507;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4572;
    (&_S4572)->primal_0 = R_30;
    (&_S4572)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4573;
    (&_S4573)->primal_0 = vert2_7;
    (&_S4573)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4572, &_S4573, _S4571);
    float3  _S4574 = _S4569.differential_0 + _S4551;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4575;
    (&_S4575)->primal_0 = R_30;
    (&_S4575)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4576;
    (&_S4576)->primal_0 = vert1_7;
    (&_S4576)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4575, &_S4576, _S4574);
    float3  _S4577 = _S4570.differential_0 + _S4448 + make_float3 (_S4567.x, _S4567.y, _S4564);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4578;
    (&_S4578)->primal_0 = R_30;
    (&_S4578)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4579;
    (&_S4579)->primal_0 = vert0_7;
    (&_S4579)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4578, &_S4579, _S4577);
    float3  _S4580 = _S4573.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4581;
    (&_S4581)->primal_0 = _S4161;
    (&_S4581)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4582;
    (&_S4582)->primal_0 = _S4166;
    (&_S4582)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4581, &_S4582, _S4580);
    float _S4583 = - _S4582.differential_0.y;
    float _S4584 = _S4165 * _S4582.differential_0.x;
    float _S4585 = - (_S4157 * _S4582.differential_0.x);
    float3  _S4586 = _S4576.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4587;
    (&_S4587)->primal_0 = _S4161;
    (&_S4587)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4588;
    (&_S4588)->primal_0 = _S4164;
    (&_S4588)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4587, &_S4588, _S4586);
    float _S4589 = _S4157 * _S4588.differential_0.x;
    float _S4590 = _S4163 * _S4588.differential_0.x;
    float3  _S4591 = _S4579.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4592;
    (&_S4592)->primal_0 = _S4161;
    (&_S4592)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4593;
    (&_S4593)->primal_0 = _S4162;
    (&_S4593)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4592, &_S4593, _S4591);
    Matrix<float, 3, 3>  _S4594 = transpose_0(_S4581.differential_0 + _S4587.differential_0 + _S4592.differential_0);
    float _S4595 = 2.0f * - _S4594.rows[int(2)].z;
    float _S4596 = 2.0f * _S4594.rows[int(2)].y;
    float _S4597 = 2.0f * _S4594.rows[int(2)].x;
    float _S4598 = 2.0f * _S4594.rows[int(1)].z;
    float _S4599 = 2.0f * - _S4594.rows[int(1)].y;
    float _S4600 = 2.0f * _S4594.rows[int(1)].x;
    float _S4601 = 2.0f * _S4594.rows[int(0)].z;
    float _S4602 = 2.0f * _S4594.rows[int(0)].y;
    float _S4603 = 2.0f * - _S4594.rows[int(0)].x;
    float _S4604 = - _S4600 + _S4602;
    float _S4605 = _S4597 + - _S4601;
    float _S4606 = - _S4596 + _S4598;
    float _S4607 = _S4596 + _S4598;
    float _S4608 = _S4597 + _S4601;
    float _S4609 = _S4600 + _S4602;
    float _S4610 = quat_31.w * (_S4599 + _S4603);
    float _S4611 = quat_31.z * (_S4595 + _S4603);
    float _S4612 = quat_31.y * (_S4595 + _S4599);
    float _S4613 = quat_31.x * _S4604 + quat_31.z * _S4607 + quat_31.y * _S4608 + _S4610 + _S4610;
    float _S4614 = quat_31.x * _S4605 + quat_31.w * _S4607 + quat_31.y * _S4609 + _S4611 + _S4611;
    float _S4615 = quat_31.x * _S4606 + quat_31.w * _S4608 + quat_31.z * _S4609 + _S4612 + _S4612;
    float _S4616 = quat_31.w * _S4604 + quat_31.z * _S4605 + quat_31.y * _S4606;
    float _S4617 = _S4585 + _S4589;
    float _S4618 = 0.5f * - _S4617;
    float _S4619 = _S4583 + _S4588.differential_0.y;
    DiffPair_float_0 _S4620;
    (&_S4620)->primal_0 = _S4158;
    (&_S4620)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4620, _S4619);
    float _S4621 = _S4618 + _S4620.differential_0;
    float _S4622 = _S4584 + _S4590 + _S4593.differential_0.x;
    DiffPair_float_0 _S4623;
    (&_S4623)->primal_0 = _S4156;
    (&_S4623)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4623, _S4622);
    float _S4624 = _S4618 + _S4623.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4625;
    (&_S4625)->primal_0 = mean_c_26;
    (&_S4625)->differential_0 = _S4304;
    s_bwd_length_impl_0(&_S4625, 0.0f);
    float3  _S4626 = _S4625.differential_0 + _S4307.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4627;
    (&_S4627)->primal_0 = R_30;
    (&_S4627)->differential_0 = _S4385;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4628;
    (&_S4628)->primal_0 = mean_30;
    (&_S4628)->differential_0 = _S4304;
    s_bwd_prop_mul_1(&_S4627, &_S4628, _S4626);
    float3  _S4629 = _S4571 + _S4574 + _S4577 + _S4626 + _S4388.differential_0;
    Matrix<float, 3, 3>  _S4630 = _S4572.differential_0 + _S4575.differential_0 + _S4578.differential_0 + _S4627.differential_0 + _S4389;
    float3  _S4631 = make_float3 (_S4624, _S4621, _S4617);
    float4  _S4632 = make_float4 (0.0f);
    *&((&_S4632)->w) = _S4613;
    *&((&_S4632)->z) = _S4614;
    *&((&_S4632)->y) = _S4615;
    *&((&_S4632)->x) = _S4616;
    float4  _S4633 = _S4632;
    float3  _S4634 = _S4580 + _S4586 + _S4591 + _S4628.differential_0 + _S4383;
    *v_mean_10 = _S4634;
    *v_quat_9 = _S4633;
    *v_scale_9 = _S4631;
    *v_hardness_5 = _S4476;
    (*v_sh_coeffs_8)[int(0)] = _S4454;
    (*v_sh_coeffs_8)[int(1)] = _S4455;
    (*v_sh_coeffs_8)[int(2)] = _S4456;
    (*v_sh_coeffs_8)[int(3)] = _S4457;
    (*v_sh_coeffs_8)[int(4)] = _S4458;
    (*v_sh_coeffs_8)[int(5)] = _S4459;
    (*v_sh_coeffs_8)[int(6)] = _S4460;
    (*v_sh_coeffs_8)[int(7)] = _S4461;
    (*v_sh_coeffs_8)[int(8)] = _S4462;
    (*v_sh_coeffs_8)[int(9)] = _S4463;
    (*v_sh_coeffs_8)[int(10)] = _S4464;
    (*v_sh_coeffs_8)[int(11)] = _S4465;
    (*v_sh_coeffs_8)[int(12)] = _S4466;
    (*v_sh_coeffs_8)[int(13)] = _S4467;
    (*v_sh_coeffs_8)[int(14)] = _S4468;
    (*v_sh_coeffs_8)[int(15)] = _S4469;
    (*v_ch_coeffs_3)[int(0)] = _S4451;
    (*v_ch_coeffs_3)[int(1)] = _S4452;
    *v_R_9 = _S4630;
    *v_t_9 = _S4629;
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle(FixedArray<float3 , 3>  * verts_4, float2  hardness_16, float3  ray_o_5, float3  ray_d_5)
{
    float3  v1v0_0 = (*verts_4)[int(1)] - (*verts_4)[int(0)];
    float3  v2v0_0 = (*verts_4)[int(2)] - (*verts_4)[int(0)];
    float3  rov0_0 = ray_o_5 - (*verts_4)[int(0)];
    float3  n_0 = cross_0(v1v0_0, v2v0_0);
    float3  q_2 = cross_0(rov0_0, ray_d_5);
    float d_0 = 1.0f / dot_0(ray_d_5, n_0);
    float u_29 = d_0 * dot_0(- q_2, v2v0_0);
    float v_29 = d_0 * dot_0(q_2, v1v0_0);
    float t_30 = d_0 * dot_0(- n_0, rov0_0);
    bool _S4635;
    if(u_29 >= 0.0f)
    {
        _S4635 = v_29 >= 0.0f;
    }
    else
    {
        _S4635 = false;
    }
    if(_S4635)
    {
        _S4635 = (u_29 + v_29) <= 1.0f;
    }
    else
    {
        _S4635 = false;
    }
    if(_S4635)
    {
        _S4635 = t_30 >= 0.0f;
    }
    else
    {
        _S4635 = false;
    }
    if(!_S4635)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_29), (v_29)))), ((F32_sqrt((0.5f))) * (1.0f - u_29 - v_29)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_16.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_16.x;
    float _S4636;
    if(opac_0 < 0.0f)
    {
        _S4636 = 0.0f;
    }
    else
    {
        _S4636 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S4636;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_4, float _s_dOut_12)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4637 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4638 = *dpray_d_4;
    float3  v1v0_1 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_1 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_1 = (*dpray_o_4).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S4639 = s_primal_ctx_cross_0(v1v0_1, v2v0_1);
    float3  _S4640 = s_primal_ctx_cross_0(rov0_1, (*dpray_d_4).primal_0);
    float _S4641 = s_primal_ctx_dot_0((*dpray_d_4).primal_0, _S4639);
    float d_1 = 1.0f / _S4641;
    float _S4642 = _S4641 * _S4641;
    float3  _S4643 = - _S4640;
    float _S4644 = s_primal_ctx_dot_0(_S4643, v2v0_1);
    float u_30 = d_1 * _S4644;
    float _S4645 = s_primal_ctx_dot_0(_S4640, v1v0_1);
    float v_30 = d_1 * _S4645;
    float3  _S4646 = - _S4639;
    float t_31 = d_1 * s_primal_ctx_dot_0(_S4646, rov0_1);
    bool _S4647;
    if(u_30 >= 0.0f)
    {
        _S4647 = v_30 >= 0.0f;
    }
    else
    {
        _S4647 = false;
    }
    if(_S4647)
    {
        _S4647 = (u_30 + v_30) <= 1.0f;
    }
    else
    {
        _S4647 = false;
    }
    if(_S4647)
    {
        _S4647 = t_31 >= 0.0f;
    }
    else
    {
        _S4647 = false;
    }
    bool _S4648 = !!_S4647;
    float _S4649;
    float _S4650;
    float _S4651;
    float _S4652;
    float _S4653;
    float _S4654;
    float _S4655;
    float _S4656;
    float _S4657;
    float _S4658;
    float _S4659;
    if(_S4648)
    {
        float _S4660 = s_primal_ctx_min_0(u_30, v_30);
        float _S4661 = s_primal_ctx_sqrt_0(0.5f);
        float _S4662 = _S4661 * (1.0f - u_30 - v_30);
        float _S4663 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S4660, _S4662) * _S4663;
        float _S4664 = _S4637.primal_0.y;
        float _S4665 = 1.0f - opac_1;
        float _S4666 = 1.0f - s_primal_ctx_clamp_0(_S4664, 0.0f, 0.99989998340606689f);
        float _S4667 = 1.0f / _S4666;
        float _S4668 = _S4666 * _S4666;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S4665, _S4667);
        float o_1 = _S4637.primal_0.x;
        bool _S4669 = opac_1 < 0.0f;
        if(_S4669)
        {
            _S4649 = 0.0f;
        }
        else
        {
            _S4649 = o_1 * w_1;
        }
        _S4647 = _S4669;
        _S4650 = o_1;
        _S4651 = w_1;
        _S4652 = _S4665;
        _S4653 = _S4667;
        _S4654 = _S4668;
        _S4655 = _S4664;
        _S4656 = _S4663;
        _S4657 = _S4660;
        _S4658 = _S4662;
        _S4659 = _S4661;
    }
    else
    {
        _S4647 = false;
        _S4649 = 0.0f;
        _S4650 = 0.0f;
        _S4651 = 0.0f;
        _S4652 = 0.0f;
        _S4653 = 0.0f;
        _S4654 = 0.0f;
        _S4655 = 0.0f;
        _S4656 = 0.0f;
        _S4657 = 0.0f;
        _S4658 = 0.0f;
        _S4659 = 0.0f;
    }
    float2  _S4670 = make_float2 (0.0f);
    float2  _S4671;
    if(_S4648)
    {
        if(_S4647)
        {
            _S4649 = 0.0f;
            _S4650 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S4672;
            (&_S4672)->primal_0 = _S4649;
            (&_S4672)->differential_0 = 0.0f;
            DiffPair_float_0 _S4673;
            (&_S4673)->primal_0 = 0.99500000476837158f;
            (&_S4673)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S4672, &_S4673, _s_dOut_12);
            float _S4674 = _S4650 * _S4672.differential_0;
            _S4649 = _S4651 * _S4672.differential_0;
            _S4650 = _S4674;
        }
        float _S4675 = - _S4650;
        DiffPair_float_0 _S4676;
        (&_S4676)->primal_0 = _S4652;
        (&_S4676)->differential_0 = 0.0f;
        DiffPair_float_0 _S4677;
        (&_S4677)->primal_0 = _S4653;
        (&_S4677)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4676, &_S4677, _S4675);
        float _S4678 = - - (_S4677.differential_0 / _S4654);
        float s_diff_opac_T_0 = - _S4676.differential_0;
        DiffPair_float_0 _S4679;
        (&_S4679)->primal_0 = _S4655;
        (&_S4679)->differential_0 = 0.0f;
        DiffPair_float_0 _S4680;
        (&_S4680)->primal_0 = 0.0f;
        (&_S4680)->differential_0 = 0.0f;
        DiffPair_float_0 _S4681;
        (&_S4681)->primal_0 = 0.99989998340606689f;
        (&_S4681)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S4679, &_S4680, &_S4681, _S4678);
        float _S4682 = _S4656 * s_diff_opac_T_0;
        DiffPair_float_0 _S4683;
        (&_S4683)->primal_0 = _S4657;
        (&_S4683)->differential_0 = 0.0f;
        DiffPair_float_0 _S4684;
        (&_S4684)->primal_0 = _S4658;
        (&_S4684)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4683, &_S4684, _S4682);
        float _S4685 = - (_S4659 * _S4684.differential_0);
        DiffPair_float_0 _S4686;
        (&_S4686)->primal_0 = u_30;
        (&_S4686)->differential_0 = 0.0f;
        DiffPair_float_0 _S4687;
        (&_S4687)->primal_0 = v_30;
        (&_S4687)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4686, &_S4687, _S4683.differential_0);
        float2  _S4688 = make_float2 (_S4649, _S4679.differential_0);
        float _S4689 = _S4685 + _S4687.differential_0;
        _S4649 = _S4685 + _S4686.differential_0;
        _S4650 = _S4689;
        _S4671 = _S4688;
    }
    else
    {
        _S4649 = 0.0f;
        _S4650 = 0.0f;
        _S4671 = _S4670;
    }
    float3  _S4690 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4691;
    (&_S4691)->primal_0 = _S4646;
    (&_S4691)->differential_0 = _S4690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4692;
    (&_S4692)->primal_0 = rov0_1;
    (&_S4692)->differential_0 = _S4690;
    s_bwd_prop_dot_0(&_S4691, &_S4692, 0.0f);
    float3  _S4693 = - _S4691.differential_0;
    float _S4694 = d_1 * _S4650;
    float _S4695 = _S4645 * _S4650;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4696;
    (&_S4696)->primal_0 = _S4640;
    (&_S4696)->differential_0 = _S4690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4697;
    (&_S4697)->primal_0 = v1v0_1;
    (&_S4697)->differential_0 = _S4690;
    s_bwd_prop_dot_0(&_S4696, &_S4697, _S4694);
    float _S4698 = d_1 * _S4649;
    float _S4699 = _S4644 * _S4649;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4700;
    (&_S4700)->primal_0 = _S4643;
    (&_S4700)->differential_0 = _S4690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4701;
    (&_S4701)->primal_0 = v2v0_1;
    (&_S4701)->differential_0 = _S4690;
    s_bwd_prop_dot_0(&_S4700, &_S4701, _S4698);
    float3  _S4702 = - _S4700.differential_0;
    float _S4703 = - ((_S4695 + _S4699) / _S4642);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4704;
    (&_S4704)->primal_0 = _S4638.primal_0;
    (&_S4704)->differential_0 = _S4690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4705;
    (&_S4705)->primal_0 = _S4639;
    (&_S4705)->differential_0 = _S4690;
    s_bwd_prop_dot_0(&_S4704, &_S4705, _S4703);
    float3  _S4706 = _S4696.differential_0 + _S4702;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4707;
    (&_S4707)->primal_0 = rov0_1;
    (&_S4707)->differential_0 = _S4690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4708;
    (&_S4708)->primal_0 = _S4638.primal_0;
    (&_S4708)->differential_0 = _S4690;
    s_bwd_prop_cross_0(&_S4707, &_S4708, _S4706);
    float3  _S4709 = _S4693 + _S4705.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4710;
    (&_S4710)->primal_0 = v1v0_1;
    (&_S4710)->differential_0 = _S4690;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4711;
    (&_S4711)->primal_0 = v2v0_1;
    (&_S4711)->differential_0 = _S4690;
    s_bwd_prop_cross_0(&_S4710, &_S4711, _S4709);
    float3  _S4712 = _S4692.differential_0 + _S4707.differential_0;
    float3  _S4713 = _S4701.differential_0 + _S4711.differential_0;
    float3  _S4714 = _S4697.differential_0 + _S4710.differential_0;
    float3  _S4715 = - _S4712 + - _S4713 + - _S4714;
    float3  _S4716 = _S4704.differential_0 + _S4708.differential_0;
    dpray_d_4->primal_0 = (*dpray_d_4).primal_0;
    dpray_d_4->differential_0 = _S4716;
    dpray_o_4->primal_0 = (*dpray_o_4).primal_0;
    dpray_o_4->differential_0 = _S4712;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S4671;
    FixedArray<float3 , 3>  _S4717;
    _S4717[int(0)] = _S4690;
    _S4717[int(1)] = _S4690;
    _S4717[int(2)] = _S4690;
    _S4717[int(2)] = _S4713;
    _S4717[int(0)] = _S4715;
    _S4717[int(1)] = _S4714;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S4717;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4718, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4719, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4720, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4721, float _S4722)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S4718, _S4719, _S4720, _S4721, _S4722);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_5, float2  hardness_17, float3  ray_o_6, float3  ray_d_6, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S4723 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4724 = { _S4723, _S4723, _S4723 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_5;
    (&dp_verts_0)->differential_0 = _S4724;
    float2  _S4725 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S4725;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_6;
    (&dp_ray_o_2)->differential_0 = _S4723;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_6;
    (&dp_ray_d_2)->differential_0 = _S4723;
    s_bwd_evaluate_alpha_opaque_triangle_0(&dp_verts_0, &dp_hardness_2, &dp_ray_o_2, &dp_ray_d_2, v_alpha_3);
    *v_verts_2 = (&dp_verts_0)->differential_0;
    *v_hardness_6 = dp_hardness_2.differential_0;
    *v_ray_o_3 = dp_ray_o_2.differential_0;
    *v_ray_d_3 = dp_ray_d_2.differential_0;
    return;
}

inline __device__ float evaluate_sorting_depth_opaque_triangle(FixedArray<float3 , 3>  * verts_6, FixedArray<float3 , 3>  * rgbs_4, float3  ray_o_7, float3  ray_d_7)
{
    float3  n_1 = cross_0((*verts_6)[int(1)] - (*verts_6)[int(0)], (*verts_6)[int(2)] - (*verts_6)[int(0)]);
    return 1.0f / dot_0(ray_d_7, n_1) * dot_0(- n_1, ray_o_7 - (*verts_6)[int(0)]);
}

inline __device__ void evaluate_color_opaque_triangle(FixedArray<float3 , 3>  * verts_7, FixedArray<float3 , 3>  * rgbs_5, float3  ray_o_8, float3  ray_d_8, float3  * color_13, float * depth_20)
{
    float3  v1v0_2 = (*verts_7)[int(1)] - (*verts_7)[int(0)];
    float3  v2v0_2 = (*verts_7)[int(2)] - (*verts_7)[int(0)];
    float3  rov0_2 = ray_o_8 - (*verts_7)[int(0)];
    float3  n_2 = cross_0(v1v0_2, v2v0_2);
    float3  q_3 = cross_0(rov0_2, ray_d_8);
    float d_2 = 1.0f / dot_0(ray_d_8, n_2);
    float u_31 = d_2 * dot_0(- q_3, v2v0_2);
    float v_31 = d_2 * dot_0(q_3, v1v0_2);
    *depth_20 = d_2 * dot_0(- n_2, rov0_2);
    *color_13 = (*rgbs_5)[int(0)] * make_float3 (1.0f - u_31 - v_31) + (*rgbs_5)[int(1)] * make_float3 (u_31) + (*rgbs_5)[int(2)] * make_float3 (v_31);
    *depth_20 = (F32_log(((F32_max((*depth_20), (9.999999960041972e-13f))))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float3  dpcolor_1, float dpdepth_2)
{
    float3  v1v0_3 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_3 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_3 = (*dpray_o_5).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S4726 = s_primal_ctx_cross_0(v1v0_3, v2v0_3);
    float3  _S4727 = s_primal_ctx_cross_0(rov0_3, (*dpray_d_5).primal_0);
    float _S4728 = s_primal_ctx_dot_0((*dpray_d_5).primal_0, _S4726);
    float d_3 = 1.0f / _S4728;
    float _S4729 = _S4728 * _S4728;
    float3  _S4730 = - _S4727;
    float _S4731 = s_primal_ctx_dot_0(_S4730, v2v0_3);
    float u_32 = d_3 * _S4731;
    float3  _S4732 = make_float3 (u_32);
    float _S4733 = s_primal_ctx_dot_0(_S4727, v1v0_3);
    float v_32 = d_3 * _S4733;
    float3  _S4734 = make_float3 (v_32);
    float3  _S4735 = - _S4726;
    float _S4736 = s_primal_ctx_dot_0(_S4735, rov0_3);
    float _S4737 = d_3 * _S4736;
    float3  _S4738 = make_float3 (1.0f - u_32 - v_32);
    DiffPair_float_0 _S4739;
    (&_S4739)->primal_0 = s_primal_ctx_max_0(_S4737, 9.999999960041972e-13f);
    (&_S4739)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4739, dpdepth_2);
    DiffPair_float_0 _S4740;
    (&_S4740)->primal_0 = _S4737;
    (&_S4740)->differential_0 = 0.0f;
    DiffPair_float_0 _S4741;
    (&_S4741)->primal_0 = 9.999999960041972e-13f;
    (&_S4741)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4740, &_S4741, _S4739.differential_0);
    float3  _S4742 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S4743 = _S4734 * dpcolor_1;
    float3  _S4744 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S4745 = _S4732 * dpcolor_1;
    float3  _S4746 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S4747 = _S4738 * dpcolor_1;
    float _S4748 = - (_S4746.x + _S4746.y + _S4746.z);
    float _S4749 = d_3 * _S4740.differential_0;
    float _S4750 = _S4736 * _S4740.differential_0;
    float3  _S4751 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4752;
    (&_S4752)->primal_0 = _S4735;
    (&_S4752)->differential_0 = _S4751;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4753;
    (&_S4753)->primal_0 = rov0_3;
    (&_S4753)->differential_0 = _S4751;
    s_bwd_prop_dot_0(&_S4752, &_S4753, _S4749);
    float3  _S4754 = - _S4752.differential_0;
    float _S4755 = _S4748 + _S4742.x + _S4742.y + _S4742.z;
    float _S4756 = d_3 * _S4755;
    float _S4757 = _S4733 * _S4755;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4758;
    (&_S4758)->primal_0 = _S4727;
    (&_S4758)->differential_0 = _S4751;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4759;
    (&_S4759)->primal_0 = v1v0_3;
    (&_S4759)->differential_0 = _S4751;
    s_bwd_prop_dot_0(&_S4758, &_S4759, _S4756);
    float _S4760 = _S4748 + _S4744.x + _S4744.y + _S4744.z;
    float _S4761 = d_3 * _S4760;
    float _S4762 = _S4731 * _S4760;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4763;
    (&_S4763)->primal_0 = _S4730;
    (&_S4763)->differential_0 = _S4751;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4764;
    (&_S4764)->primal_0 = v2v0_3;
    (&_S4764)->differential_0 = _S4751;
    s_bwd_prop_dot_0(&_S4763, &_S4764, _S4761);
    float3  _S4765 = - _S4763.differential_0;
    float _S4766 = - ((_S4750 + _S4757 + _S4762) / _S4729);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4767;
    (&_S4767)->primal_0 = (*dpray_d_5).primal_0;
    (&_S4767)->differential_0 = _S4751;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4768;
    (&_S4768)->primal_0 = _S4726;
    (&_S4768)->differential_0 = _S4751;
    s_bwd_prop_dot_0(&_S4767, &_S4768, _S4766);
    float3  _S4769 = _S4758.differential_0 + _S4765;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4770;
    (&_S4770)->primal_0 = rov0_3;
    (&_S4770)->differential_0 = _S4751;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4771;
    (&_S4771)->primal_0 = (*dpray_d_5).primal_0;
    (&_S4771)->differential_0 = _S4751;
    s_bwd_prop_cross_0(&_S4770, &_S4771, _S4769);
    float3  _S4772 = _S4754 + _S4768.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4773;
    (&_S4773)->primal_0 = v1v0_3;
    (&_S4773)->differential_0 = _S4751;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4774;
    (&_S4774)->primal_0 = v2v0_3;
    (&_S4774)->differential_0 = _S4751;
    s_bwd_prop_cross_0(&_S4773, &_S4774, _S4772);
    float3  _S4775 = _S4753.differential_0 + _S4770.differential_0;
    float3  _S4776 = _S4764.differential_0 + _S4774.differential_0;
    float3  _S4777 = _S4759.differential_0 + _S4773.differential_0;
    float3  _S4778 = - _S4775 + - _S4776 + - _S4777;
    float3  _S4779 = _S4767.differential_0 + _S4771.differential_0;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S4779;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S4775;
    FixedArray<float3 , 3>  _S4780;
    _S4780[int(0)] = _S4751;
    _S4780[int(1)] = _S4751;
    _S4780[int(2)] = _S4751;
    _S4780[int(2)] = _S4743;
    _S4780[int(1)] = _S4745;
    _S4780[int(0)] = _S4747;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S4780;
    FixedArray<float3 , 3>  _S4781;
    _S4781[int(0)] = _S4751;
    _S4781[int(1)] = _S4751;
    _S4781[int(2)] = _S4751;
    _S4781[int(2)] = _S4776;
    _S4781[int(0)] = _S4778;
    _S4781[int(1)] = _S4777;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S4781;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4782, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4783, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4784, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4785, float3  _S4786, float _S4787)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S4782, _S4783, _S4784, _S4785, _S4786, _S4787);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_8, FixedArray<float3 , 3>  * rgbs_6, float3  ray_o_9, float3  ray_d_9, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S4788 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4789 = { _S4788, _S4788, _S4788 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_8;
    (&dp_verts_1)->differential_0 = _S4789;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S4789;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_9;
    (&dp_ray_o_3)->differential_0 = _S4788;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_9;
    (&dp_ray_d_3)->differential_0 = _S4788;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

