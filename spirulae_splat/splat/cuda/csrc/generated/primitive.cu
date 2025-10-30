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

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_8, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_8, float3  dOut_4)
{
    float _S5 = (*right_8).primal_0.rows[int(0)].x * dOut_4.x;
    Matrix<float, 3, 3>  right_d_result_4;
    *&(((&right_d_result_4)->rows + (int(0)))->x) = (*left_8).primal_0.x * dOut_4.x;
    float sum_14 = _S5 + (*right_8).primal_0.rows[int(0)].y * dOut_4.y;
    *&(((&right_d_result_4)->rows + (int(0)))->y) = (*left_8).primal_0.x * dOut_4.y;
    float sum_15 = sum_14 + (*right_8).primal_0.rows[int(0)].z * dOut_4.z;
    *&(((&right_d_result_4)->rows + (int(0)))->z) = (*left_8).primal_0.x * dOut_4.z;
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = sum_15;
    float _S6 = (*right_8).primal_0.rows[int(1)].x * dOut_4.x;
    *&(((&right_d_result_4)->rows + (int(1)))->x) = (*left_8).primal_0.y * dOut_4.x;
    float sum_16 = _S6 + (*right_8).primal_0.rows[int(1)].y * dOut_4.y;
    *&(((&right_d_result_4)->rows + (int(1)))->y) = (*left_8).primal_0.y * dOut_4.y;
    float sum_17 = sum_16 + (*right_8).primal_0.rows[int(1)].z * dOut_4.z;
    *&(((&right_d_result_4)->rows + (int(1)))->z) = (*left_8).primal_0.y * dOut_4.z;
    *&((&left_d_result_4)->y) = sum_17;
    float _S7 = (*right_8).primal_0.rows[int(2)].x * dOut_4.x;
    *&(((&right_d_result_4)->rows + (int(2)))->x) = (*left_8).primal_0.z * dOut_4.x;
    float sum_18 = _S7 + (*right_8).primal_0.rows[int(2)].y * dOut_4.y;
    *&(((&right_d_result_4)->rows + (int(2)))->y) = (*left_8).primal_0.z * dOut_4.y;
    float sum_19 = sum_18 + (*right_8).primal_0.rows[int(2)].z * dOut_4.z;
    *&(((&right_d_result_4)->rows + (int(2)))->z) = (*left_8).primal_0.z * dOut_4.z;
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
        float _S8 = dist_coeffs_0->radial_coeffs_0.z + r2_0 * k4_0;
        float _S9 = dist_coeffs_0->radial_coeffs_0.y + r2_0 * _S8;
        float _S10 = dist_coeffs_0->radial_coeffs_0.x + r2_0 * _S9;
        float radial_0 = 1.0f + r2_0 * _S10;
        float _S11 = 2.0f * p1_0;
        float _S12 = _S11 * u_0;
        float _S13 = 2.0f * u_0;
        float _S14 = 2.0f * p2_0;
        float _S15 = _S14 * u_0;
        float _S16 = 2.0f * v_0;
        float2  _S17 = q_0 * make_float2 (radial_0) + make_float2 (_S12 * v_0 + p2_0 * (r2_0 + _S13 * u_0) + sx1_0 * r2_0, _S15 * v_0 + p1_0 * (r2_0 + _S16 * v_0) + sy1_0 * r2_0);
        float2  _S18 = make_float2 (0.0f);
        float2  seed_0 = _S18;
        *&((&seed_0)->x) = 1.0f;
        float2  _S19 = make_float2 (radial_0);
        float2  _S20 = q_0 * seed_0;
        float _S21 = p1_0 * seed_0.y;
        float _S22 = p2_0 * seed_0.x;
        float _S23 = _S20.x + _S20.y;
        float _S24 = r2_0 * _S23;
        float _S25 = r2_0 * _S24;
        float _S26 = sy1_0 * seed_0.y + _S21 + sx1_0 * seed_0.x + _S22 + _S10 * _S23 + _S9 * _S24 + _S8 * _S25 + k4_0 * (r2_0 * _S25);
        float _S27 = v_0 * _S26;
        float _S28 = u_0 * _S26;
        Matrix<float, 2, 2>  J_0;
        J_0[int(0)] = _S19 * seed_0 + make_float2 (_S14 * (v_0 * seed_0.y) + _S13 * _S22 + 2.0f * (u_0 * _S22) + _S11 * (v_0 * seed_0.x) + _S28 + _S28, _S16 * _S21 + 2.0f * (v_0 * _S21) + _S15 * seed_0.y + _S12 * seed_0.x + _S27 + _S27);
        float2  seed_1 = _S18;
        *&((&seed_1)->y) = 1.0f;
        float2  _S29 = q_0 * seed_1;
        float _S30 = p1_0 * seed_1.y;
        float _S31 = p2_0 * seed_1.x;
        float _S32 = _S29.x + _S29.y;
        float _S33 = r2_0 * _S32;
        float _S34 = r2_0 * _S33;
        float _S35 = sy1_0 * seed_1.y + _S30 + sx1_0 * seed_1.x + _S31 + _S10 * _S32 + _S9 * _S33 + _S8 * _S34 + k4_0 * (r2_0 * _S34);
        float _S36 = v_0 * _S35;
        float _S37 = u_0 * _S35;
        J_0[int(1)] = _S19 * seed_1 + make_float2 (_S14 * (v_0 * seed_1.y) + _S13 * _S31 + 2.0f * (u_0 * _S31) + _S11 * (v_0 * seed_1.x) + _S37 + _S37, _S16 * _S30 + 2.0f * (v_0 * _S30) + _S15 * seed_1.y + _S12 * seed_1.x + _S36 + _S36);
        float2  _S38 = _S17 - uv_0;
        float inv_det_0 = 1.0f / (J_0.rows[int(0)].x * J_0.rows[int(1)].y - J_0.rows[int(0)].y * J_0.rows[int(1)].x);
        float _S39 = _S38.x;
        float _S40 = _S38.y;
        float2  q_1 = q_0 - make_float2 ((_S39 * J_0.rows[int(1)].y - _S40 * J_0.rows[int(0)].y) * inv_det_0, (- _S39 * J_0.rows[int(1)].x + _S40 * J_0.rows[int(0)].x) * inv_det_0);
        i_5 = i_5 + int(1);
        q_0 = q_1;
    }
    return q_0;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpy_0, float dOut_5)
{
    DiffPair_float_0 _S41 = *dpx_0;
    float _S42;
    if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
    {
        _S42 = dOut_5;
    }
    else
    {
        if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
        {
            _S42 = 0.0f;
        }
        else
        {
            _S42 = 0.5f * dOut_5;
        }
    }
    dpx_0->primal_0 = _S41.primal_0;
    dpx_0->differential_0 = _S42;
    DiffPair_float_0 _S43 = *dpy_0;
    if(((*dpy_0).primal_0) > (_S41.primal_0))
    {
        _S42 = dOut_5;
    }
    else
    {
        if(((*dpy_0).primal_0) < ((*dpx_0).primal_0))
        {
            _S42 = 0.0f;
        }
        else
        {
            _S42 = 0.5f * dOut_5;
        }
    }
    dpy_0->primal_0 = _S43.primal_0;
    dpy_0->differential_0 = _S42;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_1, float dOut_6)
{
    float _S44 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_1).primal_0)))))) * dOut_6;
    dpx_1->primal_0 = (*dpx_1).primal_0;
    dpx_1->differential_0 = _S44;
    return;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_1, float dOut_7)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_1).primal_0.x * dOut_7;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_2).primal_0.x * dOut_7;
    *&((&x_d_result_0)->y) = (*dpy_1).primal_0.y * dOut_7;
    *&((&y_d_result_0)->y) = (*dpx_2).primal_0.y * dOut_7;
    *&((&x_d_result_0)->z) = (*dpy_1).primal_0.z * dOut_7;
    *&((&y_d_result_0)->z) = (*dpx_2).primal_0.z * dOut_7;
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

inline __device__ void _d_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpy_2, float dOut_8)
{
    float2  x_d_result_1;
    *&((&x_d_result_1)->x) = (*dpy_2).primal_0.x * dOut_8;
    float2  y_d_result_1;
    *&((&y_d_result_1)->x) = (*dpx_3).primal_0.x * dOut_8;
    *&((&x_d_result_1)->y) = (*dpy_2).primal_0.y * dOut_8;
    *&((&y_d_result_1)->y) = (*dpx_3).primal_0.y * dOut_8;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = x_d_result_1;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_1;
    return;
}

inline __device__ float dot_0(float3  x_6, float3  y_0)
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
        float result_11 = result_10 + _slang_vector_get_element(x_6, i_6) * _slang_vector_get_element(y_0, i_6);
        i_6 = i_6 + int(1);
        result_10 = result_11;
    }
    return result_10;
}

inline __device__ float dot_1(float2  x_7, float2  y_1)
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
        float result_13 = result_12 + _slang_vector_get_element(x_7, i_7) * _slang_vector_get_element(y_1, i_7);
        i_7 = i_7 + int(1);
        result_12 = result_13;
    }
    return result_12;
}

inline __device__ float length_0(float2  x_8)
{
    return (F32_sqrt((dot_1(x_8, x_8))));
}

inline __device__ float length_1(float3  x_9)
{
    return (F32_sqrt((dot_0(x_9, x_9))));
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float2  normalize_1(float2  x_11)
{
    return x_11 / make_float2 (length_0(x_11));
}

inline __device__ void generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_1, bool is_fisheye_0, float4  radial_coeffs_1, float2  tangential_coeffs_1, float2  thin_prism_coeffs_1, float3  * ray_o_0, float3  * ray_d_0)
{
    *ray_o_0 = - mul_7(t_1, R_2);
    CameraDistortion_0 _S45;
    (&_S45)->radial_coeffs_0 = radial_coeffs_1;
    (&_S45)->tangential_coeffs_0 = tangential_coeffs_1;
    (&_S45)->thin_prism_coeffs_0 = thin_prism_coeffs_1;
    float2  _S46 = undistort_point_0(uv_1, &_S45, int(8));
    float3  raydir_0;
    if(is_fisheye_0)
    {
        float theta_0 = length_0(_S46);
        float _S47;
        if(theta_0 < 0.00100000004749745f)
        {
            _S47 = 1.0f - theta_0 * theta_0 / 6.0f;
        }
        else
        {
            _S47 = (F32_sin((theta_0))) / theta_0;
        }
        float3  _S48 = make_float3 ((_S46 * make_float2 (_S47)).x, (_S46 * make_float2 (_S47)).y, (F32_cos((theta_0))));
        raydir_0 = _S48;
    }
    else
    {
        raydir_0 = make_float3 (_S46.x, _S46.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_7(raydir_0, R_2));
    return;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S49;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S50, Matrix<float, 3, 3>  _S51)
{
    return mul_7(_S50, _S51);
}

inline __device__ float s_primal_ctx_sin_0(float _S52)
{
    return (F32_sin((_S52)));
}

inline __device__ float s_primal_ctx_cos_0(float _S53)
{
    return (F32_cos((_S53)));
}

inline __device__ void s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_2, bool is_fisheye_1, float4  radial_coeffs_2, float2  tangential_coeffs_2, float2  thin_prism_coeffs_2, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S49 = make_float2 (0.0f);
    float3  _S54 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    CameraDistortion_0 _S55;
    (&_S55)->radial_coeffs_0 = radial_coeffs_2;
    (&_S55)->tangential_coeffs_0 = tangential_coeffs_2;
    (&_S55)->thin_prism_coeffs_0 = thin_prism_coeffs_2;
    float2  _S56 = undistort_point_0(uv_2, &_S55, int(8));
    _s_diff_ctx_0->_S49 = _S56;
    float3  raydir_1;
    if(is_fisheye_1)
    {
        float _S57 = length_0(_S56);
        float _S58;
        if(_S57 < 0.00100000004749745f)
        {
            _S58 = 1.0f - _S57 * _S57 / 6.0f;
        }
        else
        {
            _S58 = s_primal_ctx_sin_0(_S57) / _S57;
        }
        float3  _S59 = make_float3 ((_S56 * make_float2 (_S58)).x, (_S56 * make_float2 (_S58)).y, s_primal_ctx_cos_0(_S57));
        raydir_1 = _S59;
    }
    else
    {
        raydir_1 = make_float3 (_S56.x, _S56.y, 1.0f);
    }
    float3  _S60 = normalize_0(s_primal_ctx_mul_0(raydir_1, dpR_0));
    *dpray_o_0 = _S54;
    *dpray_d_0 = _S60;
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S61, float _S62)
{
    _d_sqrt_0(_S61, _S62);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_4, float _s_dOut_0)
{
    float _S63 = (*dpx_4).primal_0.x;
    float _S64 = (*dpx_4).primal_0.y;
    float _S65 = (*dpx_4).primal_0.z;
    DiffPair_float_0 _S66;
    (&_S66)->primal_0 = _S63 * _S63 + _S64 * _S64 + _S65 * _S65;
    (&_S66)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S66, _s_dOut_0);
    float _S67 = (*dpx_4).primal_0.z * _S66.differential_0;
    float _S68 = _S67 + _S67;
    float _S69 = (*dpx_4).primal_0.y * _S66.differential_0;
    float _S70 = _S69 + _S69;
    float _S71 = (*dpx_4).primal_0.x * _S66.differential_0;
    float _S72 = _S71 + _S71;
    float3  _S73 = make_float3 (0.0f);
    *&((&_S73)->z) = _S68;
    *&((&_S73)->y) = _S70;
    *&((&_S73)->x) = _S72;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = _S73;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S74, float _S75)
{
    s_bwd_prop_length_impl_0(_S74, _S75);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float3  _s_dOut_1)
{
    float _S76 = length_1((*dpx_5).primal_0);
    float3  _S77 = (*dpx_5).primal_0 * _s_dOut_1;
    float3  _S78 = make_float3 (1.0f / _S76) * _s_dOut_1;
    float _S79 = - ((_S77.x + _S77.y + _S77.z) / (_S76 * _S76));
    float3  _S80 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S81;
    (&_S81)->primal_0 = (*dpx_5).primal_0;
    (&_S81)->differential_0 = _S80;
    s_bwd_length_impl_0(&_S81, _S79);
    float3  _S82 = _S78 + _S81.differential_0;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S82;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S83, float3  _S84)
{
    s_bwd_prop_normalize_impl_0(_S83, _S84);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S85, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S86, float3  _S87)
{
    _d_mul_1(_S85, _S86, _S87);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_3, bool is_fisheye_2, float4  radial_coeffs_3, float2  tangential_coeffs_3, float2  thin_prism_coeffs_3, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S88 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S89 = *dpt_1;
    float3  raydir_2;
    if(is_fisheye_2)
    {
        float _S90 = length_0(_s_diff_ctx_1->_S49);
        float _S91;
        if(_S90 < 0.00100000004749745f)
        {
            _S91 = 1.0f - _S90 * _S90 / 6.0f;
        }
        else
        {
            _S91 = s_primal_ctx_sin_0(_S90) / _S90;
        }
        float3  _S92 = make_float3 ((_s_diff_ctx_1->_S49 * make_float2 (_S91)).x, (_s_diff_ctx_1->_S49 * make_float2 (_S91)).y, s_primal_ctx_cos_0(_S90));
        raydir_2 = _S92;
    }
    else
    {
        raydir_2 = make_float3 (_s_diff_ctx_1->_S49.x, _s_diff_ctx_1->_S49.y, 1.0f);
    }
    float3  _S93 = s_primal_ctx_mul_0(raydir_2, _S88.primal_0);
    float3  _S94 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S95;
    (&_S95)->primal_0 = _S93;
    (&_S95)->differential_0 = _S94;
    s_bwd_normalize_impl_0(&_S95, dpray_d_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S96;
    (&_S96)->primal_0 = raydir_2;
    (&_S96)->differential_0 = _S94;
    Matrix<float, 3, 3>  _S97 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S98;
    (&_S98)->primal_0 = _S88.primal_0;
    (&_S98)->differential_0 = _S97;
    s_bwd_prop_mul_0(&_S96, &_S98, _S95.differential_0);
    float3  _S99 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S100;
    (&_S100)->primal_0 = _S89.primal_0;
    (&_S100)->differential_0 = _S94;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S101;
    (&_S101)->primal_0 = _S88.primal_0;
    (&_S101)->differential_0 = _S97;
    s_bwd_prop_mul_0(&_S100, &_S101, _S99);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S100.differential_0;
    Matrix<float, 3, 3>  _S102 = _S101.differential_0 + _S98.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S102;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S103, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S104, float2  _S105, bool _S106, float4  _S107, float2  _S108, float2  _S109, float3  _S110, float3  _S111)
{
    float3  _S112;
    float3  _S113;
    s_bwd_prop_generate_ray_Intermediates_0 _S114;
    s_primal_ctx_generate_ray_0((*_S103).primal_0, (*_S104).primal_0, _S105, _S106, _S107, _S108, _S109, &_S112, &_S113, &_S114);
    s_bwd_prop_generate_ray_Intermediates_0 _S115 = _S114;
    s_bwd_prop_generate_ray_0(_S103, _S104, _S105, _S106, _S107, _S108, _S109, _S110, _S111, &_S115);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_4, bool is_fisheye_3, float4  radial_coeffs_4, float2  tangential_coeffs_4, float2  thin_prism_coeffs_4, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S116 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S116;
    float3  _S117 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S117;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_4, is_fisheye_3, radial_coeffs_4, tangential_coeffs_4, thin_prism_coeffs_4, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_6, DiffPair_float_0 * dpy_3, float dOut_9)
{
    DiffPair_float_0 _S118 = *dpx_6;
    float _S119;
    if(((*dpx_6).primal_0) < ((*dpy_3).primal_0))
    {
        _S119 = dOut_9;
    }
    else
    {
        if(((*dpx_6).primal_0) > ((*dpy_3).primal_0))
        {
            _S119 = 0.0f;
        }
        else
        {
            _S119 = 0.5f * dOut_9;
        }
    }
    dpx_6->primal_0 = _S118.primal_0;
    dpx_6->differential_0 = _S119;
    DiffPair_float_0 _S120 = *dpy_3;
    if(((*dpy_3).primal_0) < (_S118.primal_0))
    {
        _S119 = dOut_9;
    }
    else
    {
        if(((*dpy_3).primal_0) > ((*dpx_6).primal_0))
        {
            _S119 = 0.0f;
        }
        else
        {
            _S119 = 0.5f * dOut_9;
        }
    }
    dpy_3->primal_0 = _S120.primal_0;
    dpy_3->differential_0 = _S119;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S121 = float(width_0);
    float _S122 = float(height_0);
    float _S123 = 0.30000001192092896f * (0.5f * _S121 / fx_0);
    float _S124 = 0.30000001192092896f * (0.5f * _S122 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_1 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S121 - cx_0) / fx_0 + _S123), ((F32_max((- (cx_0 / fx_0 + _S123)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S122 - cy_0) / fy_0 + _S124), ((F32_max((- (cy_0 / fy_0 + _S124)), (mean3d_0.y * rz_0))))))) * rz2_0);
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
    CameraDistortion_0 _S125;
    (&_S125)->radial_coeffs_0 = radial_coeffs_5;
    (&_S125)->tangential_coeffs_0 = tangential_coeffs_5;
    (&_S125)->thin_prism_coeffs_0 = thin_prism_coeffs_5;
    return _S125;
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_7, float dOut_10)
{
    DiffPair_float_0 _S126 = *dpx_7;
    float _S127 = - (*dpy_4).primal_0 / ((*dpx_7).primal_0 * (*dpx_7).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S127;
    float _S128 = _S126.primal_0 / (_S126.primal_0 * _S126.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_10;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S128;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S129, float _S130)
{
    return (F32_atan2((_S129), (_S130)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S131, DiffPair_float_0 * _S132, float _S133)
{
    _d_atan2_0(_S131, _S132, _S133);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_8, float _s_dOut_2)
{
    float _S134 = (*dpx_8).primal_0.x;
    float _S135 = (*dpx_8).primal_0.y;
    DiffPair_float_0 _S136;
    (&_S136)->primal_0 = _S134 * _S134 + _S135 * _S135;
    (&_S136)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S136, _s_dOut_2);
    float _S137 = (*dpx_8).primal_0.y * _S136.differential_0;
    float _S138 = _S137 + _S137;
    float _S139 = (*dpx_8).primal_0.x * _S136.differential_0;
    float _S140 = _S139 + _S139;
    float2  _S141 = make_float2 (0.0f);
    *&((&_S141)->y) = _S138;
    *&((&_S141)->x) = _S140;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S141;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S142, float _S143)
{
    s_bwd_prop_length_impl_1(_S142, _S143);
    return;
}

inline __device__ void fisheye_proj_3dgs(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float fx_2, float fy_2, float cx_2, float cy_2, float4  radial_coeffs_6, float2  tangential_coeffs_6, float2  thin_prism_coeffs_6, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    CameraDistortion_0 dist_coeffs_1 = CameraDistortion_x24init_0(radial_coeffs_6, tangential_coeffs_6, thin_prism_coeffs_6);
    float2  _S144 = float2 {mean3d_2.x, mean3d_2.y};
    float r_6 = length_0(_S144);
    float _S145 = mean3d_2.z;
    float theta_1 = (F32_atan2((r_6), (_S145)));
    float k_0;
    if(theta_1 < 0.00100000004749745f)
    {
        k_0 = (1.0f - theta_1 * theta_1 / 3.0f) / _S145;
    }
    else
    {
        k_0 = theta_1 / r_6;
    }
    float2  _S146 = _S144 * make_float2 (k_0);
    float k1_0 = dist_coeffs_1.radial_coeffs_0.x;
    float k2_0 = dist_coeffs_1.radial_coeffs_0.y;
    float k3_0 = dist_coeffs_1.radial_coeffs_0.z;
    float k4_1 = dist_coeffs_1.radial_coeffs_0.w;
    float p1_1 = dist_coeffs_1.tangential_coeffs_0.x;
    float p2_1 = dist_coeffs_1.tangential_coeffs_0.y;
    float sx1_1 = dist_coeffs_1.thin_prism_coeffs_0.x;
    float sy1_1 = dist_coeffs_1.thin_prism_coeffs_0.y;
    float u_1 = _S146.x;
    float v_1 = _S146.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float _S147 = 2.0f * p1_1;
    float _S148 = 2.0f * p2_1;
    float2  _S149 = _S146 * make_float2 (1.0f + r2_1 * (k1_0 + r2_1 * (k2_0 + r2_1 * (k3_0 + r2_1 * k4_1)))) + make_float2 (_S147 * u_1 * v_1 + p2_1 * (r2_1 + 2.0f * u_1 * u_1) + sx1_1 * r2_1, _S148 * u_1 * v_1 + p1_1 * (r2_1 + 2.0f * v_1 * v_1) + sy1_1 * r2_1);
    *mean2d_2 = make_float2 (fx_2 * _S149.x + cx_2, fy_2 * _S149.y + cy_2);
    Matrix<float, 2, 3>  J_3;
    float2  _S150 = make_float2 (0.0f);
    float2  seed_2 = _S150;
    *&((&seed_2)->x) = 1.0f;
    float2  _S151 = seed_2;
    float _S152 = s_primal_ctx_atan2_0(r_6, _S145);
    bool _S153 = _S152 < 0.00100000004749745f;
    float _S154;
    float _S155;
    float _S156;
    if(_S153)
    {
        float _S157 = 1.0f - _S152 * _S152 / 3.0f;
        float _S158 = _S145 * _S145;
        k_0 = _S157 / _S145;
        _S154 = 0.0f;
        _S155 = _S158;
        _S156 = _S157;
    }
    else
    {
        float _S159 = r_6 * r_6;
        k_0 = _S152 / r_6;
        _S154 = _S159;
        _S155 = 0.0f;
        _S156 = 0.0f;
    }
    float2  _S160 = make_float2 (k_0);
    float2  _S161 = _S144 * make_float2 (k_0);
    float u_2 = _S161.x;
    float v_2 = _S161.y;
    float r2_2 = u_2 * u_2 + v_2 * v_2;
    float _S162 = k3_0 + r2_2 * k4_1;
    float _S163 = k2_0 + r2_2 * _S162;
    float _S164 = k1_0 + r2_2 * _S163;
    float _S165 = fy_2 * _S151.y;
    float _S166 = fx_2 * _S151.x;
    float2  _S167 = make_float2 (_S166, _S165);
    float2  _S168 = _S161 * _S167;
    float _S169 = p1_1 * _S165;
    float _S170 = p2_1 * _S166;
    float _S171 = _S168.x + _S168.y;
    float _S172 = r2_2 * _S171;
    float _S173 = r2_2 * _S172;
    float _S174 = sy1_1 * _S165 + _S169 + sx1_1 * _S166 + _S170 + _S164 * _S171 + _S163 * _S172 + _S162 * _S173 + k4_1 * (r2_2 * _S173);
    float _S175 = v_2 * _S174;
    float _S176 = u_2 * _S174;
    float2  _S177 = make_float2 (1.0f + r2_2 * _S164) * _S167 + make_float2 (_S148 * (v_2 * _S165) + 2.0f * u_2 * _S170 + 2.0f * (u_2 * _S170) + _S147 * (v_2 * _S166) + _S176 + _S176, 2.0f * v_2 * _S169 + 2.0f * (v_2 * _S169) + _S148 * u_2 * _S165 + _S147 * u_2 * _S166 + _S175 + _S175);
    float2  _S178 = _S144 * _S177;
    float2  _S179 = _S160 * _S177;
    float _S180 = _S178.x + _S178.y;
    if(_S153)
    {
        float _S181 = _S180 / _S155;
        float _S182 = _S156 * - _S181;
        float _S183 = _S152 * (0.3333333432674408f * - (_S145 * _S181));
        k_0 = _S183 + _S183;
        _S154 = _S182;
        _S155 = 0.0f;
    }
    else
    {
        float _S184 = _S180 / _S154;
        float _S185 = _S152 * - _S184;
        k_0 = r_6 * _S184;
        _S154 = 0.0f;
        _S155 = _S185;
    }
    DiffPair_float_0 _S186;
    (&_S186)->primal_0 = r_6;
    (&_S186)->differential_0 = 0.0f;
    DiffPair_float_0 _S187;
    (&_S187)->primal_0 = _S145;
    (&_S187)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S186, &_S187, k_0);
    float _S188 = _S187.differential_0 + _S154;
    float _S189 = _S186.differential_0 + _S155;
    float2  _S190 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S191;
    (&_S191)->primal_0 = _S144;
    (&_S191)->differential_0 = _S190;
    s_bwd_length_impl_1(&_S191, _S189);
    float2  _S192 = _S191.differential_0 + _S179;
    float3  _S193 = make_float3 (_S192.x, _S192.y, _S188);
    J_3[int(0)] = _S193;
    float2  seed_3 = _S150;
    *&((&seed_3)->y) = 1.0f;
    float2  _S194 = seed_3;
    if(_S153)
    {
        float _S195 = 1.0f - _S152 * _S152 / 3.0f;
        float _S196 = _S145 * _S145;
        k_0 = _S195 / _S145;
        _S154 = 0.0f;
        _S155 = _S196;
        _S156 = _S195;
    }
    else
    {
        float _S197 = r_6 * r_6;
        k_0 = _S152 / r_6;
        _S154 = _S197;
        _S155 = 0.0f;
        _S156 = 0.0f;
    }
    float2  _S198 = make_float2 (k_0);
    float2  _S199 = _S144 * make_float2 (k_0);
    float u_3 = _S199.x;
    float v_3 = _S199.y;
    float r2_3 = u_3 * u_3 + v_3 * v_3;
    float _S200 = k3_0 + r2_3 * k4_1;
    float _S201 = k2_0 + r2_3 * _S200;
    float _S202 = k1_0 + r2_3 * _S201;
    float _S203 = fy_2 * _S194.y;
    float _S204 = fx_2 * _S194.x;
    float2  _S205 = make_float2 (_S204, _S203);
    float2  _S206 = _S199 * _S205;
    float _S207 = p1_1 * _S203;
    float _S208 = p2_1 * _S204;
    float _S209 = _S206.x + _S206.y;
    float _S210 = r2_3 * _S209;
    float _S211 = r2_3 * _S210;
    float _S212 = sy1_1 * _S203 + _S207 + sx1_1 * _S204 + _S208 + _S202 * _S209 + _S201 * _S210 + _S200 * _S211 + k4_1 * (r2_3 * _S211);
    float _S213 = v_3 * _S212;
    float _S214 = u_3 * _S212;
    float2  _S215 = make_float2 (1.0f + r2_3 * _S202) * _S205 + make_float2 (_S148 * (v_3 * _S203) + 2.0f * u_3 * _S208 + 2.0f * (u_3 * _S208) + _S147 * (v_3 * _S204) + _S214 + _S214, 2.0f * v_3 * _S207 + 2.0f * (v_3 * _S207) + _S148 * u_3 * _S203 + _S147 * u_3 * _S204 + _S213 + _S213);
    float2  _S216 = _S144 * _S215;
    float2  _S217 = _S198 * _S215;
    float _S218 = _S216.x + _S216.y;
    if(_S153)
    {
        float _S219 = _S218 / _S155;
        float _S220 = _S156 * - _S219;
        float _S221 = _S152 * (0.3333333432674408f * - (_S145 * _S219));
        k_0 = _S221 + _S221;
        _S154 = _S220;
        _S155 = 0.0f;
    }
    else
    {
        float _S222 = _S218 / _S154;
        float _S223 = _S152 * - _S222;
        k_0 = r_6 * _S222;
        _S154 = 0.0f;
        _S155 = _S223;
    }
    DiffPair_float_0 _S224;
    (&_S224)->primal_0 = r_6;
    (&_S224)->differential_0 = 0.0f;
    DiffPair_float_0 _S225;
    (&_S225)->primal_0 = _S145;
    (&_S225)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S224, &_S225, k_0);
    float _S226 = _S225.differential_0 + _S154;
    float _S227 = _S224.differential_0 + _S155;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S228;
    (&_S228)->primal_0 = _S144;
    (&_S228)->differential_0 = _S190;
    s_bwd_length_impl_1(&_S228, _S227);
    float2  _S229 = _S228.differential_0 + _S217;
    float3  _S230 = make_float3 (_S229.x, _S229.y, _S226);
    J_3[int(1)] = _S230;
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
    float _S231 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S231;
    float _S232 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S232;
    float det_blur_0 = _S231 * _S232 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_9, float dOut_11)
{
    float _S233 = (F32_exp(((*dpx_9).primal_0))) * dOut_11;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S233;
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
    float3  _S234 = exp_0((*dpx_10).primal_0) * dOut_12;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S234;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_11, float dOut_13)
{
    float _S235 = 1.0f / (*dpx_11).primal_0 * dOut_13;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S235;
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

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_4, float fy_4, float cx_4, float cy_4, float4  radial_coeffs_7, float2  tangential_coeffs_7, float2  thin_prism_coeffs_7, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_4, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_4, mean_0) + t_3;
        float _S236 = mean_c_0.z;
        bool _S237;
        if(_S236 < near_plane_0)
        {
            _S237 = true;
        }
        else
        {
            _S237 = _S236 > far_plane_0;
        }
        if(_S237)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S238 = exp_0(scale_2);
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
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S238.x, 0.0f, 0.0f, 0.0f, _S238.y, 0.0f, 0.0f, 0.0f, _S238.z));
        Matrix<float, 3, 3>  _S239 = transpose_0(R_4);
        Matrix<float, 3, 3>  covar_c_0 = mul_4(mul_4(R_4, mul_4(M_2, transpose_0(M_2))), _S239);
        Matrix<float, 2, 2>  covar2d_0;
        float _S240 = float(image_width_0);
        float _S241 = float(image_height_0);
        float _S242 = 0.30000001192092896f * (0.5f * _S240 / fx_4);
        float _S243 = 0.30000001192092896f * (0.5f * _S241 / fy_4);
        float rz_2 = 1.0f / mean_c_0.z;
        float rz2_2 = rz_2 * rz_2;
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_4 * rz_2, 0.0f, - fx_4 * (mean_c_0.z * (F32_min(((_S240 - cx_4) / fx_4 + _S242), ((F32_max((- (cx_4 / fx_4 + _S242)), (mean_c_0.x * rz_2))))))) * rz2_2, 0.0f, fy_4 * rz_2, - fy_4 * (mean_c_0.z * (F32_min(((_S241 - cy_4) / fy_4 + _S243), ((F32_max((- (cy_4 / fy_4 + _S243)), (mean_c_0.y * rz_2))))))) * rz2_2);
        covar2d_0 = mul_6(mul_5(J_5, covar_c_0), transpose_1(J_5));
        *mean2d_4 = make_float2 (fx_4 * mean_c_0.x * rz_2 + cx_4, fy_4 * mean_c_0.y * rz_2 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S244 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S244;
        float _S245 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S245;
        float det_blur_1 = _S244 * _S245 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S246 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S237 = true;
        }
        else
        {
            _S237 = xmin_0 >= _S240;
        }
        if(_S237)
        {
            _S237 = true;
        }
        else
        {
            _S237 = ymax_0 <= 0.0f;
        }
        if(_S237)
        {
            _S237 = true;
        }
        else
        {
            _S237 = ymin_0 >= _S241;
        }
        if(_S237)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = 0.5f * (F32_log((dot_0(mean_c_0, mean_c_0) + 9.99999997475242708e-07f)));
        *conic_0 = make_float3 (_S246.rows[int(0)].x, _S246.rows[int(0)].y, _S246.rows[int(1)].y);
        float3  _S247 = mean_0 - - mul_0(_S239, t_3);
        float3  _S248 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S248;
        float _S249 = _S247.x;
        float _S250 = _S247.y;
        float _S251 = _S247.z;
        float norm_0 = (F32_sqrt((_S249 * _S249 + _S250 * _S250 + _S251 * _S251)));
        float x_15 = _S249 / norm_0;
        float y_3 = _S250 / norm_0;
        float z_0 = _S251 / norm_0;
        float3  _S252 = _S248 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_15) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S252;
        float z2_4 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_15 * x_15 - y_3 * y_3;
        float fS1_0 = 2.0f * x_15 * y_3;
        float3  _S253 = _S252 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_4 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_15) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S253;
        float fTmp0C_0 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgb_0 = max_0(_S253 + (make_float3 (-0.59004360437393188f * (x_15 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_4 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_15) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_15 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_5, float fy_5, float cx_5, float cy_5, float4  radial_coeffs_8, float2  tangential_coeffs_8, float2  thin_prism_coeffs_8, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_5, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_5, mean_1) + t_4;
        float _S254 = length_1(mean_c_1);
        bool _S255;
        if(_S254 < near_plane_1)
        {
            _S255 = true;
        }
        else
        {
            _S255 = _S254 > far_plane_1;
        }
        if(_S255)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S256 = exp_0(scale_3);
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
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_5), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_5), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S256.x, 0.0f, 0.0f, 0.0f, _S256.y, 0.0f, 0.0f, 0.0f, _S256.z));
        Matrix<float, 3, 3>  _S257 = transpose_0(R_5);
        Matrix<float, 3, 3>  covar_c_1 = mul_4(mul_4(R_5, mul_4(M_3, transpose_0(M_3))), _S257);
        Matrix<float, 2, 2>  covar2d_1;
        fisheye_proj_3dgs(mean_c_1, covar_c_1, fx_5, fy_5, cx_5, cy_5, radial_coeffs_8, tangential_coeffs_8, thin_prism_coeffs_8, &covar2d_1, mean2d_5);
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S258 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S258;
        float _S259 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S259;
        float det_blur_2 = _S258 * _S259 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S260 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S255 = true;
        }
        else
        {
            _S255 = xmin_1 >= float(image_width_1);
        }
        if(_S255)
        {
            _S255 = true;
        }
        else
        {
            _S255 = ymax_1 <= 0.0f;
        }
        if(_S255)
        {
            _S255 = true;
        }
        else
        {
            _S255 = ymin_1 >= float(image_height_1);
        }
        if(_S255)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = 0.5f * (F32_log((dot_0(mean_c_1, mean_c_1) + 9.99999997475242708e-07f)));
        *conic_1 = make_float3 (_S260.rows[int(0)].x, _S260.rows[int(0)].y, _S260.rows[int(1)].y);
        float3  _S261 = mean_1 - - mul_0(_S257, t_4);
        float3  _S262 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S262;
        float _S263 = _S261.x;
        float _S264 = _S261.y;
        float _S265 = _S261.z;
        float norm_1 = (F32_sqrt((_S263 * _S263 + _S264 * _S264 + _S265 * _S265)));
        float x_17 = _S263 / norm_1;
        float y_4 = _S264 / norm_1;
        float z_1 = _S265 / norm_1;
        float3  _S266 = _S262 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_1)[int(1)] + make_float3 (z_1) * (*sh_coeffs_1)[int(2)] - make_float3 (x_17) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S266;
        float z2_6 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_17 * x_17 - y_4 * y_4;
        float fS1_1 = 2.0f * x_17 * y_4;
        float3  _S267 = _S266 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_4) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_6 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_17) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S267;
        float fTmp0C_1 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgb_1 = max_0(_S267 + (make_float3 (-0.59004360437393188f * (x_17 * fS1_1 + y_4 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_4) * (*sh_coeffs_1)[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_6 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_17) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_1 - y_4 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_6, float fy_6, float cx_6, float cy_6, float4  radial_coeffs_9, float2  tangential_coeffs_9, float2  thin_prism_coeffs_9, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_6, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_6, mean_2) + t_5;
        float _S268 = mean_c_2.z;
        bool _S269;
        if(_S268 < near_plane_2)
        {
            _S269 = true;
        }
        else
        {
            _S269 = _S268 > far_plane_2;
        }
        if(_S269)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S270 = exp_0(scale_4);
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
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_7), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_7), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S270.x, 0.0f, 0.0f, 0.0f, _S270.y, 0.0f, 0.0f, 0.0f, _S270.z));
        Matrix<float, 3, 3>  _S271 = transpose_0(R_6);
        Matrix<float, 3, 3>  covar_c_2 = mul_4(mul_4(R_6, mul_4(M_4, transpose_0(M_4))), _S271);
        Matrix<float, 2, 2>  covar2d_2;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
        covar2d_2 = mul_6(mul_5(J_6, covar_c_2), transpose_1(J_6));
        *mean2d_6 = make_float2 (fx_6 * mean_c_2.x + cx_6, fy_6 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S272 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S272;
        float _S273 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S273;
        float det_blur_3 = _S272 * _S273 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S274 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S269 = true;
        }
        else
        {
            _S269 = xmin_2 >= float(image_width_2);
        }
        if(_S269)
        {
            _S269 = true;
        }
        else
        {
            _S269 = ymax_2 <= 0.0f;
        }
        if(_S269)
        {
            _S269 = true;
        }
        else
        {
            _S269 = ymin_2 >= float(image_height_2);
        }
        if(_S269)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = 0.5f * (F32_log((dot_0(mean_c_2, mean_c_2) + 9.99999997475242708e-07f)));
        *conic_2 = make_float3 (_S274.rows[int(0)].x, _S274.rows[int(0)].y, _S274.rows[int(1)].y);
        float3  _S275 = mean_2 - - mul_0(_S271, t_5);
        float3  _S276 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S276;
        float _S277 = _S275.x;
        float _S278 = _S275.y;
        float _S279 = _S275.z;
        float norm_2 = (F32_sqrt((_S277 * _S277 + _S278 * _S278 + _S279 * _S279)));
        float x_19 = _S277 / norm_2;
        float y_5 = _S278 / norm_2;
        float z_2 = _S279 / norm_2;
        float3  _S280 = _S276 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_2)[int(1)] + make_float3 (z_2) * (*sh_coeffs_2)[int(2)] - make_float3 (x_19) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S280;
        float z2_8 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_19 * x_19 - y_5 * y_5;
        float fS1_2 = 2.0f * x_19 * y_5;
        float3  _S281 = _S280 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_5) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_8 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_19) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S281;
        float fTmp0C_2 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        *rgb_2 = max_0(_S281 + (make_float3 (-0.59004360437393188f * (x_19 * fS1_2 + y_5 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_5) * (*sh_coeffs_2)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_8 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_19) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_19 * fC1_2 - y_5 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_persp(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_7, float fy_7, float cx_7, float cy_7, float4  radial_coeffs_10, float2  tangential_coeffs_10, float2  thin_prism_coeffs_10, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_7, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_0(R_7, mean_3) + t_6;
        float _S282 = mean_c_3.z;
        bool _S283;
        if(_S282 < near_plane_3)
        {
            _S283 = true;
        }
        else
        {
            _S283 = _S282 > far_plane_3;
        }
        if(_S283)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S284 = exp_0(scale_5);
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
        Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_9), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_9), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S284.x, 0.0f, 0.0f, 0.0f, _S284.y, 0.0f, 0.0f, 0.0f, _S284.z));
        Matrix<float, 3, 3>  _S285 = transpose_0(R_7);
        Matrix<float, 3, 3>  covar_c_3 = mul_4(mul_4(R_7, mul_4(M_5, transpose_0(M_5))), _S285);
        Matrix<float, 2, 2>  covar2d_3;
        float _S286 = float(image_width_3);
        float _S287 = float(image_height_3);
        float _S288 = 0.30000001192092896f * (0.5f * _S286 / fx_7);
        float _S289 = 0.30000001192092896f * (0.5f * _S287 / fy_7);
        float rz_3 = 1.0f / mean_c_3.z;
        float rz2_3 = rz_3 * rz_3;
        Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_7 * rz_3, 0.0f, - fx_7 * (mean_c_3.z * (F32_min(((_S286 - cx_7) / fx_7 + _S288), ((F32_max((- (cx_7 / fx_7 + _S288)), (mean_c_3.x * rz_3))))))) * rz2_3, 0.0f, fy_7 * rz_3, - fy_7 * (mean_c_3.z * (F32_min(((_S287 - cy_7) / fy_7 + _S289), ((F32_max((- (cy_7 / fy_7 + _S289)), (mean_c_3.y * rz_3))))))) * rz2_3);
        covar2d_3 = mul_6(mul_5(J_7, covar_c_3), transpose_1(J_7));
        *mean2d_7 = make_float2 (fx_7 * mean_c_3.x * rz_3 + cx_7, fy_7 * mean_c_3.y * rz_3 + cy_7);
        float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S290 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S290;
        float _S291 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S291;
        float det_blur_4 = _S290 * _S291 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
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
            _S283 = true;
        }
        else
        {
            _S283 = xmin_3 >= _S286;
        }
        if(_S283)
        {
            _S283 = true;
        }
        else
        {
            _S283 = ymax_3 <= 0.0f;
        }
        if(_S283)
        {
            _S283 = true;
        }
        else
        {
            _S283 = ymin_3 >= _S287;
        }
        if(_S283)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_3 = make_int4 (int(xmin_3), int(ymin_3), int(xmax_3), int(ymax_3));
        *depth_3 = 0.5f * (F32_log((dot_0(mean_c_3, mean_c_3) + 9.99999997475242708e-07f)));
        *conic_3 = exp_0(- scale_5);
        float3  _S292 = mean_3 - - mul_0(_S285, t_6);
        float3  _S293 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
        *rgb_3 = _S293;
        float _S294 = _S292.x;
        float _S295 = _S292.y;
        float _S296 = _S292.z;
        float norm_3 = (F32_sqrt((_S294 * _S294 + _S295 * _S295 + _S296 * _S296)));
        float x_21 = _S294 / norm_3;
        float y_6 = _S295 / norm_3;
        float z_3 = _S296 / norm_3;
        float3  _S297 = _S293 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_3)[int(1)] + make_float3 (z_3) * (*sh_coeffs_3)[int(2)] - make_float3 (x_21) * (*sh_coeffs_3)[int(3)]);
        *rgb_3 = _S297;
        float z2_10 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_21 * x_21 - y_6 * y_6;
        float fS1_3 = 2.0f * x_21 * y_6;
        float3  _S298 = _S297 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_6) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_10 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_21) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
        *rgb_3 = _S298;
        float fTmp0C_3 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        *rgb_3 = max_0(_S298 + (make_float3 (-0.59004360437393188f * (x_21 * fS1_3 + y_6 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_6) * (*sh_coeffs_3)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_10 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_21) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_21 * fC1_3 - y_6 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_fisheye(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_8, float fy_8, float cx_8, float cy_8, float4  radial_coeffs_11, float2  tangential_coeffs_11, float2  thin_prism_coeffs_11, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_8, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    for(;;)
    {
        float3  mean_c_4 = mul_0(R_8, mean_4) + t_7;
        float _S299 = length_1(mean_c_4);
        bool _S300;
        if(_S299 < near_plane_4)
        {
            _S300 = true;
        }
        else
        {
            _S300 = _S299 > far_plane_4;
        }
        if(_S300)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  _S301 = exp_0(scale_6);
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
        Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_11), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_11), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))), makeMatrix<float, 3, 3> (_S301.x, 0.0f, 0.0f, 0.0f, _S301.y, 0.0f, 0.0f, 0.0f, _S301.z));
        Matrix<float, 3, 3>  _S302 = transpose_0(R_8);
        Matrix<float, 3, 3>  covar_c_4 = mul_4(mul_4(R_8, mul_4(M_6, transpose_0(M_6))), _S302);
        Matrix<float, 2, 2>  covar2d_4;
        fisheye_proj_3dgs(mean_c_4, covar_c_4, fx_8, fy_8, cx_8, cy_8, radial_coeffs_11, tangential_coeffs_11, thin_prism_coeffs_11, &covar2d_4, mean2d_8);
        float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S303 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S303;
        float _S304 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S304;
        float det_blur_5 = _S303 * _S304 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
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
            _S300 = true;
        }
        else
        {
            _S300 = xmin_4 >= float(image_width_4);
        }
        if(_S300)
        {
            _S300 = true;
        }
        else
        {
            _S300 = ymax_4 <= 0.0f;
        }
        if(_S300)
        {
            _S300 = true;
        }
        else
        {
            _S300 = ymin_4 >= float(image_height_4);
        }
        if(_S300)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_4 = make_int4 (int(xmin_4), int(ymin_4), int(xmax_4), int(ymax_4));
        *depth_4 = 0.5f * (F32_log((dot_0(mean_c_4, mean_c_4) + 9.99999997475242708e-07f)));
        *conic_4 = exp_0(- scale_6);
        float3  _S305 = mean_4 - - mul_0(_S302, t_7);
        float3  _S306 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
        *rgb_4 = _S306;
        float _S307 = _S305.x;
        float _S308 = _S305.y;
        float _S309 = _S305.z;
        float norm_4 = (F32_sqrt((_S307 * _S307 + _S308 * _S308 + _S309 * _S309)));
        float x_23 = _S307 / norm_4;
        float y_7 = _S308 / norm_4;
        float z_4 = _S309 / norm_4;
        float3  _S310 = _S306 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_4)[int(1)] + make_float3 (z_4) * (*sh_coeffs_4)[int(2)] - make_float3 (x_23) * (*sh_coeffs_4)[int(3)]);
        *rgb_4 = _S310;
        float z2_12 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_23 * x_23 - y_7 * y_7;
        float fS1_4 = 2.0f * x_23 * y_7;
        float3  _S311 = _S310 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_7) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_23) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
        *rgb_4 = _S311;
        float fTmp0C_4 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        *rgb_4 = max_0(_S311 + (make_float3 (-0.59004360437393188f * (x_23 * fS1_4 + y_7 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_7) * (*sh_coeffs_4)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_23) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_23 * fC1_4 - y_7 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_9, float fy_9, float cx_9, float cy_9, float4  radial_coeffs_12, float2  tangential_coeffs_12, float2  thin_prism_coeffs_12, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_9, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_9, mean_5) + t_8;
    float3  _S312 = exp_0(scale_7);
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
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_13), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_13), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S312.x, 0.0f, 0.0f, 0.0f, _S312.y, 0.0f, 0.0f, 0.0f, _S312.z));
    Matrix<float, 3, 3>  _S313 = transpose_0(R_9);
    float _S314 = float(image_width_5);
    float _S315 = float(image_height_5);
    float _S316 = 0.30000001192092896f * (0.5f * _S314 / fx_9);
    float _S317 = 0.30000001192092896f * (0.5f * _S315 / fy_9);
    float rz_4 = 1.0f / mean_c_5.z;
    float rz2_4 = rz_4 * rz_4;
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_9 * rz_4, 0.0f, - fx_9 * (mean_c_5.z * (F32_min(((_S314 - cx_9) / fx_9 + _S316), ((F32_max((- (cx_9 / fx_9 + _S316)), (mean_c_5.x * rz_4))))))) * rz2_4, 0.0f, fy_9 * rz_4, - fy_9 * (mean_c_5.z * (F32_min(((_S315 - cy_9) / fy_9 + _S317), ((F32_max((- (cy_9 / fy_9 + _S317)), (mean_c_5.y * rz_4))))))) * rz2_4);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_8, mul_4(mul_4(R_9, mul_4(M_7, transpose_0(M_7))), _S313)), transpose_1(J_8));
    *mean2d_9 = make_float2 (fx_9 * mean_c_5.x * rz_4 + cx_9, fy_9 * mean_c_5.y * rz_4 + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S318 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S318;
    float _S319 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S319;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S318 * _S319 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S320 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_4, - covar2d_5.rows[int(0)].y * invdet_4, - covar2d_5.rows[int(1)].x * invdet_4, covar2d_5.rows[int(0)].x * invdet_4);
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
    *conic_5 = make_float3 (_S320.rows[int(0)].x, _S320.rows[int(0)].y, _S320.rows[int(1)].y);
    float3  _S321 = mean_5 - - mul_0(_S313, t_8);
    float3  _S322 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S322;
    float _S323 = _S321.x;
    float _S324 = _S321.y;
    float _S325 = _S321.z;
    float norm_5 = (F32_sqrt((_S323 * _S323 + _S324 * _S324 + _S325 * _S325)));
    float x_25 = _S323 / norm_5;
    float y_8 = _S324 / norm_5;
    float z_5 = _S325 / norm_5;
    float3  _S326 = _S322 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_5)[int(1)] + make_float3 (z_5) * (*sh_coeffs_5)[int(2)] - make_float3 (x_25) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S326;
    float z2_14 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_25 * x_25 - y_8 * y_8;
    float fS1_5 = 2.0f * x_25 * y_8;
    float3  _S327 = _S326 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_8) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_14 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_25) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S327;
    float fTmp0C_5 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    *rgb_5 = max_0(_S327 + (make_float3 (-0.59004360437393188f * (x_25 * fS1_5 + y_8 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_8) * (*sh_coeffs_5)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_14 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_25) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_25 * fC1_5 - y_8 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_10, float fy_10, float cx_10, float cy_10, float4  radial_coeffs_13, float2  tangential_coeffs_13, float2  thin_prism_coeffs_13, uint image_width_6, uint image_height_6, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * mean2d_10, float * depth_6, float3  * conic_6, float * opacity_6, float3  * rgb_6)
{
    float3  mean_c_6 = mul_0(R_10, mean_6) + t_9;
    float3  _S328 = exp_0(scale_8);
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
    Matrix<float, 3, 3>  M_8 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S328.x, 0.0f, 0.0f, 0.0f, _S328.y, 0.0f, 0.0f, 0.0f, _S328.z));
    Matrix<float, 3, 3>  _S329 = transpose_0(R_10);
    Matrix<float, 2, 2>  covar2d_6;
    fisheye_proj_3dgs(mean_c_6, mul_4(mul_4(R_10, mul_4(M_8, transpose_0(M_8))), _S329), fx_10, fy_10, cx_10, cy_10, radial_coeffs_13, tangential_coeffs_13, thin_prism_coeffs_13, &covar2d_6, mean2d_10);
    float det_orig_7 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
    float _S330 = *&(((&covar2d_6)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(0)))->x) = _S330;
    float _S331 = *&(((&covar2d_6)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(1)))->y) = _S331;
    float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / (_S330 * _S331 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
    Matrix<float, 2, 2>  _S332 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_5, - covar2d_6.rows[int(0)].y * invdet_5, - covar2d_6.rows[int(1)].x * invdet_5, covar2d_6.rows[int(0)].x * invdet_5);
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
    *conic_6 = make_float3 (_S332.rows[int(0)].x, _S332.rows[int(0)].y, _S332.rows[int(1)].y);
    float3  _S333 = mean_6 - - mul_0(_S329, t_9);
    float3  _S334 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)];
    *rgb_6 = _S334;
    float _S335 = _S333.x;
    float _S336 = _S333.y;
    float _S337 = _S333.z;
    float norm_6 = (F32_sqrt((_S335 * _S335 + _S336 * _S336 + _S337 * _S337)));
    float x_27 = _S335 / norm_6;
    float y_9 = _S336 / norm_6;
    float z_6 = _S337 / norm_6;
    float3  _S338 = _S334 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_6)[int(1)] + make_float3 (z_6) * (*sh_coeffs_6)[int(2)] - make_float3 (x_27) * (*sh_coeffs_6)[int(3)]);
    *rgb_6 = _S338;
    float z2_16 = z_6 * z_6;
    float fTmp0B_6 = -1.09254848957061768f * z_6;
    float fC1_6 = x_27 * x_27 - y_9 * y_9;
    float fS1_6 = 2.0f * x_27 * y_9;
    float3  _S339 = _S338 + (make_float3 (0.54627424478530884f * fS1_6) * (*sh_coeffs_6)[int(4)] + make_float3 (fTmp0B_6 * y_9) * (*sh_coeffs_6)[int(5)] + make_float3 (0.94617468118667603f * z2_16 - 0.31539157032966614f) * (*sh_coeffs_6)[int(6)] + make_float3 (fTmp0B_6 * x_27) * (*sh_coeffs_6)[int(7)] + make_float3 (0.54627424478530884f * fC1_6) * (*sh_coeffs_6)[int(8)]);
    *rgb_6 = _S339;
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_6;
    *rgb_6 = max_0(_S339 + (make_float3 (-0.59004360437393188f * (x_27 * fS1_6 + y_9 * fC1_6)) * (*sh_coeffs_6)[int(9)] + make_float3 (fTmp1B_6 * fS1_6) * (*sh_coeffs_6)[int(10)] + make_float3 (fTmp0C_6 * y_9) * (*sh_coeffs_6)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_16 - 1.11952900886535645f)) * (*sh_coeffs_6)[int(12)] + make_float3 (fTmp0C_6 * x_27) * (*sh_coeffs_6)[int(13)] + make_float3 (fTmp1B_6 * fC1_6) * (*sh_coeffs_6)[int(14)] + make_float3 (-0.59004360437393188f * (x_27 * fC1_6 - y_9 * fS1_6)) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_11, float fy_11, float cx_11, float cy_11, float4  radial_coeffs_14, float2  tangential_coeffs_14, float2  thin_prism_coeffs_14, uint image_width_7, uint image_height_7, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * mean2d_11, float * depth_7, float3  * conic_7, float * opacity_7, float3  * rgb_7)
{
    float3  mean_c_7 = mul_0(R_11, mean_7) + t_10;
    float3  _S340 = exp_0(scale_9);
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
    Matrix<float, 3, 3>  M_9 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S340.x, 0.0f, 0.0f, 0.0f, _S340.y, 0.0f, 0.0f, 0.0f, _S340.z));
    Matrix<float, 3, 3>  _S341 = transpose_0(R_11);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_11, 0.0f, 0.0f, 0.0f, fy_11, 0.0f);
    Matrix<float, 2, 2>  covar2d_7 = mul_6(mul_5(J_9, mul_4(mul_4(R_11, mul_4(M_9, transpose_0(M_9))), _S341)), transpose_1(J_9));
    *mean2d_11 = make_float2 (fx_11 * mean_c_7.x + cx_11, fy_11 * mean_c_7.y + cy_11);
    float det_orig_8 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
    float _S342 = *&(((&covar2d_7)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(0)))->x) = _S342;
    float _S343 = *&(((&covar2d_7)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(1)))->y) = _S343;
    float compensation_8 = (F32_sqrt(((F32_max((0.0f), (det_orig_8 / (_S342 * _S343 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
    Matrix<float, 2, 2>  _S344 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_6, - covar2d_7.rows[int(0)].y * invdet_6, - covar2d_7.rows[int(1)].x * invdet_6, covar2d_7.rows[int(0)].x * invdet_6);
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
    *conic_7 = make_float3 (_S344.rows[int(0)].x, _S344.rows[int(0)].y, _S344.rows[int(1)].y);
    float3  _S345 = mean_7 - - mul_0(_S341, t_10);
    float3  _S346 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)];
    *rgb_7 = _S346;
    float _S347 = _S345.x;
    float _S348 = _S345.y;
    float _S349 = _S345.z;
    float norm_7 = (F32_sqrt((_S347 * _S347 + _S348 * _S348 + _S349 * _S349)));
    float x_29 = _S347 / norm_7;
    float y_10 = _S348 / norm_7;
    float z_7 = _S349 / norm_7;
    float3  _S350 = _S346 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * (*sh_coeffs_7)[int(1)] + make_float3 (z_7) * (*sh_coeffs_7)[int(2)] - make_float3 (x_29) * (*sh_coeffs_7)[int(3)]);
    *rgb_7 = _S350;
    float z2_18 = z_7 * z_7;
    float fTmp0B_7 = -1.09254848957061768f * z_7;
    float fC1_7 = x_29 * x_29 - y_10 * y_10;
    float fS1_7 = 2.0f * x_29 * y_10;
    float3  _S351 = _S350 + (make_float3 (0.54627424478530884f * fS1_7) * (*sh_coeffs_7)[int(4)] + make_float3 (fTmp0B_7 * y_10) * (*sh_coeffs_7)[int(5)] + make_float3 (0.94617468118667603f * z2_18 - 0.31539157032966614f) * (*sh_coeffs_7)[int(6)] + make_float3 (fTmp0B_7 * x_29) * (*sh_coeffs_7)[int(7)] + make_float3 (0.54627424478530884f * fC1_7) * (*sh_coeffs_7)[int(8)]);
    *rgb_7 = _S351;
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_7;
    *rgb_7 = max_0(_S351 + (make_float3 (-0.59004360437393188f * (x_29 * fS1_7 + y_10 * fC1_7)) * (*sh_coeffs_7)[int(9)] + make_float3 (fTmp1B_7 * fS1_7) * (*sh_coeffs_7)[int(10)] + make_float3 (fTmp0C_7 * y_10) * (*sh_coeffs_7)[int(11)] + make_float3 (z_7 * (1.86588168144226074f * z2_18 - 1.11952900886535645f)) * (*sh_coeffs_7)[int(12)] + make_float3 (fTmp0C_7 * x_29) * (*sh_coeffs_7)[int(13)] + make_float3 (fTmp1B_7 * fC1_7) * (*sh_coeffs_7)[int(14)] + make_float3 (-0.59004360437393188f * (x_29 * fC1_7 - y_10 * fS1_7)) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_persp_differentiable(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_12, float fy_12, float cx_12, float cy_12, float4  radial_coeffs_15, float2  tangential_coeffs_15, float2  thin_prism_coeffs_15, uint image_width_8, uint image_height_8, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * mean2d_12, float * depth_8, float3  * conic_8, float * opacity_8, float3  * rgb_8)
{
    float3  mean_c_8 = mul_0(R_12, mean_8) + t_11;
    float3  _S352 = exp_0(scale_10);
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
    Matrix<float, 3, 3>  M_10 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))), makeMatrix<float, 3, 3> (_S352.x, 0.0f, 0.0f, 0.0f, _S352.y, 0.0f, 0.0f, 0.0f, _S352.z));
    Matrix<float, 3, 3>  _S353 = transpose_0(R_12);
    float _S354 = float(image_width_8);
    float _S355 = float(image_height_8);
    float _S356 = 0.30000001192092896f * (0.5f * _S354 / fx_12);
    float _S357 = 0.30000001192092896f * (0.5f * _S355 / fy_12);
    float rz_5 = 1.0f / mean_c_8.z;
    float rz2_5 = rz_5 * rz_5;
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (fx_12 * rz_5, 0.0f, - fx_12 * (mean_c_8.z * (F32_min(((_S354 - cx_12) / fx_12 + _S356), ((F32_max((- (cx_12 / fx_12 + _S356)), (mean_c_8.x * rz_5))))))) * rz2_5, 0.0f, fy_12 * rz_5, - fy_12 * (mean_c_8.z * (F32_min(((_S355 - cy_12) / fy_12 + _S357), ((F32_max((- (cy_12 / fy_12 + _S357)), (mean_c_8.y * rz_5))))))) * rz2_5);
    Matrix<float, 2, 2>  covar2d_8 = mul_6(mul_5(J_10, mul_4(mul_4(R_12, mul_4(M_10, transpose_0(M_10))), _S353)), transpose_1(J_10));
    *mean2d_12 = make_float2 (fx_12 * mean_c_8.x * rz_5 + cx_12, fy_12 * mean_c_8.y * rz_5 + cy_12);
    float det_orig_9 = *&(((&covar2d_8)->rows + (int(0)))->x) * *&(((&covar2d_8)->rows + (int(1)))->y) - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x);
    float _S358 = *&(((&covar2d_8)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(0)))->x) = _S358;
    float _S359 = *&(((&covar2d_8)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(1)))->y) = _S359;
    float compensation_9 = (F32_sqrt(((F32_max((0.0f), (det_orig_9 / (_S358 * _S359 - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x))))))));
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
    float3  _S360 = mean_8 - - mul_0(_S353, t_11);
    float3  _S361 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)];
    *rgb_8 = _S361;
    float _S362 = _S360.x;
    float _S363 = _S360.y;
    float _S364 = _S360.z;
    float norm_8 = (F32_sqrt((_S362 * _S362 + _S363 * _S363 + _S364 * _S364)));
    float x_31 = _S362 / norm_8;
    float y_11 = _S363 / norm_8;
    float z_8 = _S364 / norm_8;
    float3  _S365 = _S361 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_8)[int(1)] + make_float3 (z_8) * (*sh_coeffs_8)[int(2)] - make_float3 (x_31) * (*sh_coeffs_8)[int(3)]);
    *rgb_8 = _S365;
    float z2_20 = z_8 * z_8;
    float fTmp0B_8 = -1.09254848957061768f * z_8;
    float fC1_8 = x_31 * x_31 - y_11 * y_11;
    float fS1_8 = 2.0f * x_31 * y_11;
    float3  _S366 = _S365 + (make_float3 (0.54627424478530884f * fS1_8) * (*sh_coeffs_8)[int(4)] + make_float3 (fTmp0B_8 * y_11) * (*sh_coeffs_8)[int(5)] + make_float3 (0.94617468118667603f * z2_20 - 0.31539157032966614f) * (*sh_coeffs_8)[int(6)] + make_float3 (fTmp0B_8 * x_31) * (*sh_coeffs_8)[int(7)] + make_float3 (0.54627424478530884f * fC1_8) * (*sh_coeffs_8)[int(8)]);
    *rgb_8 = _S366;
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_8;
    *rgb_8 = max_0(_S366 + (make_float3 (-0.59004360437393188f * (x_31 * fS1_8 + y_11 * fC1_8)) * (*sh_coeffs_8)[int(9)] + make_float3 (fTmp1B_8 * fS1_8) * (*sh_coeffs_8)[int(10)] + make_float3 (fTmp0C_8 * y_11) * (*sh_coeffs_8)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_20 - 1.11952900886535645f)) * (*sh_coeffs_8)[int(12)] + make_float3 (fTmp0C_8 * x_31) * (*sh_coeffs_8)[int(13)] + make_float3 (fTmp1B_8 * fC1_8) * (*sh_coeffs_8)[int(14)] + make_float3 (-0.59004360437393188f * (x_31 * fC1_8 - y_11 * fS1_8)) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_fisheye_differentiable(bool antialiased_9, float3  mean_9, float4  quat_12, float3  scale_11, float in_opacity_9, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_13, float fy_13, float cx_13, float cy_13, float4  radial_coeffs_16, float2  tangential_coeffs_16, float2  thin_prism_coeffs_16, uint image_width_9, uint image_height_9, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * mean2d_13, float * depth_9, float3  * conic_9, float * opacity_9, float3  * rgb_9)
{
    float3  mean_c_9 = mul_0(R_13, mean_9) + t_12;
    float3  _S367 = exp_0(scale_11);
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
    Matrix<float, 3, 3>  M_11 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))), makeMatrix<float, 3, 3> (_S367.x, 0.0f, 0.0f, 0.0f, _S367.y, 0.0f, 0.0f, 0.0f, _S367.z));
    Matrix<float, 3, 3>  _S368 = transpose_0(R_13);
    Matrix<float, 2, 2>  covar2d_9;
    fisheye_proj_3dgs(mean_c_9, mul_4(mul_4(R_13, mul_4(M_11, transpose_0(M_11))), _S368), fx_13, fy_13, cx_13, cy_13, radial_coeffs_16, tangential_coeffs_16, thin_prism_coeffs_16, &covar2d_9, mean2d_13);
    float det_orig_10 = *&(((&covar2d_9)->rows + (int(0)))->x) * *&(((&covar2d_9)->rows + (int(1)))->y) - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x);
    float _S369 = *&(((&covar2d_9)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(0)))->x) = _S369;
    float _S370 = *&(((&covar2d_9)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(1)))->y) = _S370;
    float compensation_10 = (F32_sqrt(((F32_max((0.0f), (det_orig_10 / (_S369 * _S370 - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x))))))));
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
    float3  _S371 = mean_9 - - mul_0(_S368, t_12);
    float3  _S372 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
    *rgb_9 = _S372;
    float _S373 = _S371.x;
    float _S374 = _S371.y;
    float _S375 = _S371.z;
    float norm_9 = (F32_sqrt((_S373 * _S373 + _S374 * _S374 + _S375 * _S375)));
    float x_33 = _S373 / norm_9;
    float y_12 = _S374 / norm_9;
    float z_9 = _S375 / norm_9;
    float3  _S376 = _S372 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * (*sh_coeffs_9)[int(1)] + make_float3 (z_9) * (*sh_coeffs_9)[int(2)] - make_float3 (x_33) * (*sh_coeffs_9)[int(3)]);
    *rgb_9 = _S376;
    float z2_22 = z_9 * z_9;
    float fTmp0B_9 = -1.09254848957061768f * z_9;
    float fC1_9 = x_33 * x_33 - y_12 * y_12;
    float fS1_9 = 2.0f * x_33 * y_12;
    float3  _S377 = _S376 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_12) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_22 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_33) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
    *rgb_9 = _S377;
    float fTmp0C_9 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_9;
    *rgb_9 = max_0(_S377 + (make_float3 (-0.59004360437393188f * (x_33 * fS1_9 + y_12 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_12) * (*sh_coeffs_9)[int(11)] + make_float3 (z_9 * (1.86588168144226074f * z2_22 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_33) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_33 * fC1_9 - y_12 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S378, float3  _S379)
{
    return mul_0(_S378, _S379);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S380)
{
    return exp_0(_S380);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S381, Matrix<float, 3, 3>  _S382)
{
    return mul_4(_S381, _S382);
}

inline __device__ float s_primal_ctx_max_0(float _S383, float _S384)
{
    return (F32_max((_S383), (_S384)));
}

inline __device__ float s_primal_ctx_min_0(float _S385, float _S386)
{
    return (F32_min((_S385), (_S386)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S387, Matrix<float, 3, 3>  _S388)
{
    return mul_5(_S387, _S388);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_4(Matrix<float, 2, 3>  _S389, Matrix<float, 3, 2>  _S390)
{
    return mul_6(_S389, _S390);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S391)
{
    return (F32_sqrt((_S391)));
}

inline __device__ float s_primal_ctx_exp_1(float _S392)
{
    return (F32_exp((_S392)));
}

inline __device__ float s_primal_ctx_log_0(float _S393)
{
    return (F32_log((_S393)));
}

inline __device__ float s_primal_ctx_dot_0(float3  _S394, float3  _S395)
{
    return dot_0(_S394, _S395);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S396, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S397, float3  _S398)
{
    _d_max_vector_0(_S396, _S397, _S398);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S399, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S400, float3  _S401)
{
    _d_mul_0(_S399, _S400, _S401);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S402, float _S403)
{
    _d_log_0(_S402, _S403);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S404, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S405, float _S406)
{
    _d_dot_0(_S404, _S405, _S406);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S407, DiffPair_float_0 * _S408, float _S409)
{
    _d_min_0(_S407, _S408, _S409);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S410, float _S411)
{
    _d_exp_0(_S410, _S411);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S412, DiffPair_float_0 * _S413, float _S414)
{
    _d_max_0(_S412, _S413, _S414);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S415, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S416, Matrix<float, 2, 2>  _S417)
{
    mul_3(_S415, _S416, _S417);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S418, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S419, Matrix<float, 2, 3>  _S420)
{
    mul_2(_S418, _S419, _S420);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S421, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S422, Matrix<float, 3, 3>  _S423)
{
    mul_1(_S421, _S422, _S423);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S424, float3  _S425)
{
    _d_exp_vector_0(_S424, _S425);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_14, float fy_14, float cx_14, float cy_14, float4  radial_coeffs_17, float2  tangential_coeffs_17, float2  thin_prism_coeffs_17, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_10 = s_primal_ctx_mul_1(R_14, mean_10) + t_13;
    float3  _S426 = s_primal_ctx_exp_0(scale_12);
    float _S427 = quat_13.y;
    float x2_13 = _S427 * _S427;
    float y2_13 = quat_13.z * quat_13.z;
    float z2_23 = quat_13.w * quat_13.w;
    float xy_13 = quat_13.y * quat_13.z;
    float xz_13 = quat_13.y * quat_13.w;
    float yz_13 = quat_13.z * quat_13.w;
    float wx_13 = quat_13.x * quat_13.y;
    float wy_13 = quat_13.x * quat_13.z;
    float wz_13 = quat_13.x * quat_13.w;
    Matrix<float, 3, 3>  _S428 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S426.x, 0.0f, 0.0f, 0.0f, _S426.y, 0.0f, 0.0f, 0.0f, _S426.z);
    Matrix<float, 3, 3>  _S429 = s_primal_ctx_mul_2(_S428, S_0);
    Matrix<float, 3, 3>  _S430 = transpose_0(_S429);
    Matrix<float, 3, 3>  _S431 = s_primal_ctx_mul_2(_S429, _S430);
    Matrix<float, 3, 3>  _S432 = s_primal_ctx_mul_2(R_14, _S431);
    Matrix<float, 3, 3>  _S433 = transpose_0(R_14);
    Matrix<float, 3, 3>  _S434 = s_primal_ctx_mul_2(_S432, _S433);
    float _S435 = float(image_width_10);
    float _S436 = float(image_height_10);
    float _S437 = 0.30000001192092896f * (0.5f * _S435 / fx_14);
    float lim_x_pos_0 = (_S435 - cx_14) / fx_14 + _S437;
    float _S438 = 0.30000001192092896f * (0.5f * _S436 / fy_14);
    float lim_y_pos_0 = (_S436 - cy_14) / fy_14 + _S438;
    float rz_6 = 1.0f / mean_c_10.z;
    float _S439 = mean_c_10.z * mean_c_10.z;
    float rz2_6 = rz_6 * rz_6;
    float _S440 = - (cx_14 / fx_14 + _S437);
    float _S441 = mean_c_10.x * rz_6;
    float _S442 = s_primal_ctx_max_0(_S440, _S441);
    float _S443 = s_primal_ctx_min_0(lim_x_pos_0, _S442);
    float _S444 = - (cy_14 / fy_14 + _S438);
    float _S445 = mean_c_10.y * rz_6;
    float _S446 = s_primal_ctx_max_0(_S444, _S445);
    float _S447 = s_primal_ctx_min_0(lim_y_pos_0, _S446);
    float _S448 = - fx_14;
    float _S449 = _S448 * (mean_c_10.z * _S443);
    float _S450 = - fy_14;
    float _S451 = _S450 * (mean_c_10.z * _S447);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_14 * rz_6, 0.0f, _S449 * rz2_6, 0.0f, fy_14 * rz_6, _S451 * rz2_6);
    Matrix<float, 2, 3>  _S452 = s_primal_ctx_mul_3(J_11, _S434);
    Matrix<float, 3, 2>  _S453 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S454 = s_primal_ctx_mul_4(_S452, _S453);
    float _S455 = fx_14 * mean_c_10.x;
    float _S456 = fy_14 * mean_c_10.y;
    float _S457 = _S454.rows[int(0)].y * _S454.rows[int(1)].x;
    float det_orig_11 = _S454.rows[int(0)].x * _S454.rows[int(1)].y - _S457;
    float _S458 = _S454.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S459 = _S454;
    *&(((&_S459)->rows + (int(0)))->x) = _S458;
    float _S460 = _S454.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S459)->rows + (int(1)))->y) = _S460;
    Matrix<float, 2, 2>  _S461 = _S459;
    Matrix<float, 2, 2>  _S462 = _S459;
    float det_blur_6 = _S458 * _S460 - _S457;
    float _S463 = det_orig_11 / det_blur_6;
    float _S464 = det_blur_6 * det_blur_6;
    float _S465 = s_primal_ctx_max_0(0.0f, _S463);
    float _S466 = s_primal_ctx_sqrt_0(_S465);
    float invdet_7 = 1.0f / det_blur_6;
    float _S467 = - _S454.rows[int(0)].y;
    float _S468 = - _S454.rows[int(1)].x;
    float _S469 = - in_opacity_10;
    float _S470 = 1.0f + s_primal_ctx_exp_1(_S469);
    float _S471 = 1.0f / _S470;
    float _S472 = _S470 * _S470;
    float _S473;
    if(antialiased_10)
    {
        _S473 = _S471 * _S466;
    }
    else
    {
        _S473 = _S471;
    }
    float _S474 = _S473 / 0.00392156885936856f;
    float _S475 = 2.0f * s_primal_ctx_log_0(_S474);
    float _S476 = s_primal_ctx_sqrt_0(_S475);
    float _S477 = _S461.rows[int(0)].x;
    float _S478 = _S462.rows[int(1)].y;
    float _S479 = s_primal_ctx_dot_0(mean_c_10, mean_c_10) + 9.99999997475242708e-07f;
    float3  _S480 = mean_10 - - s_primal_ctx_mul_1(_S433, t_13);
    float _S481 = _S480.x;
    float _S482 = _S480.y;
    float _S483 = _S480.z;
    float _S484 = _S481 * _S481 + _S482 * _S482 + _S483 * _S483;
    float _S485 = s_primal_ctx_sqrt_0(_S484);
    float x_34 = _S481 / _S485;
    float3  _S486 = make_float3 (x_34);
    float _S487 = _S485 * _S485;
    float y_13 = _S482 / _S485;
    float z_10 = _S483 / _S485;
    float3  _S488 = make_float3 (z_10);
    float _S489 = - y_13;
    float3  _S490 = make_float3 (_S489);
    float z2_24 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_34 * x_34 - y_13 * y_13;
    float _S491 = 2.0f * x_34;
    float fS1_10 = _S491 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_24 - 0.31539157032966614f;
    float3  _S492 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_34;
    float3  _S493 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S494 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S495 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S496 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S497 = 1.86588168144226074f * z2_24 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S497;
    float3  _S498 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_34;
    float3  _S499 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S500 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S501 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S502 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_34 * fC1_10 - y_13 * fS1_10);
    float3  _S503 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_34 * fS1_10 + y_13 * fC1_10);
    float3  _S504 = make_float3 (pSH9_0);
    float3  _S505 = make_float3 (0.0f);
    float3  _S506 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S507;
    (&_S507)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S489) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_34) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S507)->differential_0 = _S506;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S508;
    (&_S508)->primal_0 = _S505;
    (&_S508)->differential_0 = _S506;
    s_bwd_prop_max_0(&_S507, &_S508, v_rgb_0);
    float3  _S509 = _S503 * _S507.differential_0;
    float3  _S510 = (*sh_coeffs_10)[int(15)] * _S507.differential_0;
    float3  _S511 = _S501 * _S507.differential_0;
    float3  _S512 = (*sh_coeffs_10)[int(14)] * _S507.differential_0;
    float3  _S513 = _S499 * _S507.differential_0;
    float3  _S514 = (*sh_coeffs_10)[int(13)] * _S507.differential_0;
    float3  _S515 = _S498 * _S507.differential_0;
    float3  _S516 = (*sh_coeffs_10)[int(12)] * _S507.differential_0;
    float3  _S517 = _S500 * _S507.differential_0;
    float3  _S518 = (*sh_coeffs_10)[int(11)] * _S507.differential_0;
    float3  _S519 = _S502 * _S507.differential_0;
    float3  _S520 = (*sh_coeffs_10)[int(10)] * _S507.differential_0;
    float3  _S521 = _S504 * _S507.differential_0;
    float3  _S522 = (*sh_coeffs_10)[int(9)] * _S507.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S522.x + _S522.y + _S522.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S510.x + _S510.y + _S510.z);
    float _S523 = _S520.x + _S520.y + _S520.z;
    float _S524 = _S512.x + _S512.y + _S512.z;
    float _S525 = _S518.x + _S518.y + _S518.z;
    float _S526 = _S514.x + _S514.y + _S514.z;
    float _S527 = _S516.x + _S516.y + _S516.z;
    float _S528 = - s_diff_fC2_T_0;
    float3  _S529 = _S495 * _S507.differential_0;
    float3  _S530 = (*sh_coeffs_10)[int(8)] * _S507.differential_0;
    float3  _S531 = _S493 * _S507.differential_0;
    float3  _S532 = (*sh_coeffs_10)[int(7)] * _S507.differential_0;
    float3  _S533 = _S492 * _S507.differential_0;
    float3  _S534 = (*sh_coeffs_10)[int(6)] * _S507.differential_0;
    float3  _S535 = _S494 * _S507.differential_0;
    float3  _S536 = (*sh_coeffs_10)[int(5)] * _S507.differential_0;
    float3  _S537 = _S496 * _S507.differential_0;
    float3  _S538 = (*sh_coeffs_10)[int(4)] * _S507.differential_0;
    float _S539 = _S536.x + _S536.y + _S536.z;
    float _S540 = _S532.x + _S532.y + _S532.z;
    float _S541 = fTmp1B_10 * _S523 + x_34 * s_diff_fS2_T_0 + y_13 * _S528 + 0.54627424478530884f * (_S538.x + _S538.y + _S538.z);
    float _S542 = fTmp1B_10 * _S524 + y_13 * s_diff_fS2_T_0 + x_34 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S530.x + _S530.y + _S530.z);
    float _S543 = y_13 * - _S542;
    float _S544 = x_34 * _S542;
    float _S545 = z_10 * (1.86588168144226074f * (z_10 * _S527) + -2.28522896766662598f * (y_13 * _S525 + x_34 * _S526) + 0.94617468118667603f * (_S534.x + _S534.y + _S534.z));
    float3  _S546 = make_float3 (0.48860251903533936f) * _S507.differential_0;
    float3  _S547 = - _S546;
    float3  _S548 = _S486 * _S547;
    float3  _S549 = (*sh_coeffs_10)[int(3)] * _S547;
    float3  _S550 = _S488 * _S546;
    float3  _S551 = (*sh_coeffs_10)[int(2)] * _S546;
    float3  _S552 = _S490 * _S546;
    float3  _S553 = (*sh_coeffs_10)[int(1)] * _S546;
    float _S554 = (_S497 * _S527 + 1.44530570507049561f * (fS1_10 * _S523 + fC1_10 * _S524) + -1.09254848957061768f * (y_13 * _S539 + x_34 * _S540) + _S545 + _S545 + _S551.x + _S551.y + _S551.z) / _S487;
    float _S555 = _S485 * _S554;
    float _S556 = (fTmp0C_10 * _S525 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S528 + fTmp0B_10 * _S539 + _S491 * _S541 + _S543 + _S543 + - (_S553.x + _S553.y + _S553.z)) / _S487;
    float _S557 = _S485 * _S556;
    float _S558 = (fTmp0C_10 * _S526 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S540 + 2.0f * (y_13 * _S541) + _S544 + _S544 + _S549.x + _S549.y + _S549.z) / _S487;
    float _S559 = _S485 * _S558;
    float _S560 = _S483 * - _S554 + _S482 * - _S556 + _S481 * - _S558;
    DiffPair_float_0 _S561;
    (&_S561)->primal_0 = _S484;
    (&_S561)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S561, _S560);
    float _S562 = _S483 * _S561.differential_0;
    float _S563 = _S482 * _S561.differential_0;
    float _S564 = _S481 * _S561.differential_0;
    float3  _S565 = make_float3 (0.282094806432724f) * _S507.differential_0;
    float3  _S566 = make_float3 (_S559 + _S564 + _S564, _S557 + _S563 + _S563, _S555 + _S562 + _S562);
    float3  _S567 = - - _S566;
    Matrix<float, 3, 3>  _S568 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S569;
    (&_S569)->primal_0 = _S433;
    (&_S569)->differential_0 = _S568;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S570;
    (&_S570)->primal_0 = t_13;
    (&_S570)->differential_0 = _S506;
    s_bwd_prop_mul_1(&_S569, &_S570, _S567);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S571 = _S569;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S572 = _S570;
    float2  _S573 = make_float2 (0.0f);
    float2  _S574 = _S573;
    *&((&_S574)->y) = v_conic_0.z;
    float2  _S575 = _S573;
    *&((&_S575)->y) = v_conic_0.y;
    *&((&_S575)->x) = v_conic_0.x;
    float _S576 = 0.5f * v_depth_0;
    DiffPair_float_0 _S577;
    (&_S577)->primal_0 = _S479;
    (&_S577)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S577, _S576);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S578;
    (&_S578)->primal_0 = mean_c_10;
    (&_S578)->differential_0 = _S506;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S579;
    (&_S579)->primal_0 = mean_c_10;
    (&_S579)->differential_0 = _S506;
    s_bwd_prop_dot_0(&_S578, &_S579, _S577.differential_0);
    DiffPair_float_0 _S580;
    (&_S580)->primal_0 = _S478;
    (&_S580)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S580, 0.0f);
    DiffPair_float_0 _S581;
    (&_S581)->primal_0 = _S477;
    (&_S581)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S581, 0.0f);
    DiffPair_float_0 _S582;
    (&_S582)->primal_0 = 3.32999992370605469f;
    (&_S582)->differential_0 = 0.0f;
    DiffPair_float_0 _S583;
    (&_S583)->primal_0 = _S476;
    (&_S583)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S582, &_S583, 0.0f);
    DiffPair_float_0 _S584;
    (&_S584)->primal_0 = _S475;
    (&_S584)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S584, _S583.differential_0);
    float _S585 = 2.0f * _S584.differential_0;
    DiffPair_float_0 _S586;
    (&_S586)->primal_0 = _S474;
    (&_S586)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S586, _S585);
    float _S587 = v_opacity_0 + 254.9999847412109375f * _S586.differential_0;
    Matrix<float, 2, 2>  _S588 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S589 = _S588;
    _S589[int(1)] = _S574;
    _S589[int(0)] = _S575;
    Matrix<float, 2, 2>  _S590 = _S589;
    FixedArray<float3 , 16>  _S591;
    _S591[int(0)] = _S506;
    _S591[int(1)] = _S506;
    _S591[int(2)] = _S506;
    _S591[int(3)] = _S506;
    _S591[int(4)] = _S506;
    _S591[int(5)] = _S506;
    _S591[int(6)] = _S506;
    _S591[int(7)] = _S506;
    _S591[int(8)] = _S506;
    _S591[int(9)] = _S506;
    _S591[int(10)] = _S506;
    _S591[int(11)] = _S506;
    _S591[int(12)] = _S506;
    _S591[int(13)] = _S506;
    _S591[int(14)] = _S506;
    _S591[int(15)] = _S506;
    _S591[int(7)] = _S531;
    _S591[int(0)] = _S565;
    _S591[int(1)] = _S552;
    _S591[int(2)] = _S550;
    _S591[int(3)] = _S548;
    _S591[int(4)] = _S537;
    _S591[int(5)] = _S535;
    _S591[int(6)] = _S533;
    _S591[int(15)] = _S509;
    _S591[int(8)] = _S529;
    _S591[int(9)] = _S521;
    _S591[int(10)] = _S519;
    _S591[int(11)] = _S517;
    _S591[int(12)] = _S515;
    _S591[int(13)] = _S513;
    _S591[int(14)] = _S511;
    float3  _S592 = _S591[int(0)];
    float3  _S593 = _S591[int(1)];
    float3  _S594 = _S591[int(2)];
    float3  _S595 = _S591[int(3)];
    float3  _S596 = _S591[int(4)];
    float3  _S597 = _S591[int(5)];
    float3  _S598 = _S591[int(6)];
    float3  _S599 = _S591[int(7)];
    float3  _S600 = _S591[int(8)];
    float3  _S601 = _S591[int(9)];
    float3  _S602 = _S591[int(10)];
    float3  _S603 = _S591[int(11)];
    float3  _S604 = _S591[int(12)];
    float3  _S605 = _S591[int(13)];
    float3  _S606 = _S591[int(14)];
    float3  _S607 = _S591[int(15)];
    float3  _S608 = _S579.differential_0 + _S578.differential_0;
    float2  _S609 = make_float2 (0.0f, _S580.differential_0);
    float2  _S610 = make_float2 (_S581.differential_0, 0.0f);
    float _S611;
    if(antialiased_10)
    {
        float _S612 = _S471 * _S587;
        _S473 = _S466 * _S587;
        _S611 = _S612;
    }
    else
    {
        _S473 = _S587;
        _S611 = 0.0f;
    }
    float _S613 = - (_S473 / _S472);
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = _S469;
    (&_S614)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S614, _S613);
    float _S615 = - _S614.differential_0;
    float _S616 = invdet_7 * _S590.rows[int(1)].y;
    float _S617 = - (invdet_7 * _S590.rows[int(1)].x);
    float _S618 = - (invdet_7 * _S590.rows[int(0)].y);
    float _S619 = invdet_7 * _S590.rows[int(0)].x;
    float _S620 = - ((_S458 * _S590.rows[int(1)].y + _S468 * _S590.rows[int(1)].x + _S467 * _S590.rows[int(0)].y + _S460 * _S590.rows[int(0)].x) / _S464);
    DiffPair_float_0 _S621;
    (&_S621)->primal_0 = _S465;
    (&_S621)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S621, _S611);
    DiffPair_float_0 _S622;
    (&_S622)->primal_0 = 0.0f;
    (&_S622)->differential_0 = 0.0f;
    DiffPair_float_0 _S623;
    (&_S623)->primal_0 = _S463;
    (&_S623)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S622, &_S623, _S621.differential_0);
    float _S624 = _S623.differential_0 / _S464;
    float s_diff_det_orig_T_0 = det_blur_6 * _S624;
    float _S625 = _S620 + det_orig_11 * - _S624;
    float _S626 = - _S625;
    float _S627 = _S458 * _S625;
    float _S628 = _S460 * _S625;
    Matrix<float, 2, 2>  _S629 = _S588;
    _S629[int(1)] = _S609;
    _S629[int(0)] = _S610;
    _S459 = _S629;
    *&(((&_S459)->rows + (int(1)))->y) = 0.0f;
    float _S630 = _S619 + _S627 + _S629.rows[int(1)].y;
    *&(((&_S459)->rows + (int(0)))->x) = 0.0f;
    float _S631 = _S616 + _S628 + _S629.rows[int(0)].x;
    float _S632 = _S626 + - s_diff_det_orig_T_0;
    float _S633 = _S617 + _S454.rows[int(0)].y * _S632;
    float _S634 = _S618 + _S454.rows[int(1)].x * _S632;
    float _S635 = _S454.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S636 = _S630 + _S454.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S637 = _S573;
    *&((&_S637)->x) = _S633;
    *&((&_S637)->y) = _S636;
    float _S638 = _S631 + _S635;
    float2  _S639 = _S573;
    *&((&_S639)->y) = _S634;
    *&((&_S639)->x) = _S638;
    float _S640 = _S456 * v_mean2d_0.y;
    float _S641 = fy_14 * (rz_6 * v_mean2d_0.y);
    float _S642 = _S455 * v_mean2d_0.x;
    float _S643 = fx_14 * (rz_6 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S644 = _S588;
    _S644[int(1)] = _S637;
    _S644[int(0)] = _S639;
    Matrix<float, 2, 2>  _S645 = _S459 + _S644;
    Matrix<float, 2, 3>  _S646 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S647;
    (&_S647)->primal_0 = _S452;
    (&_S647)->differential_0 = _S646;
    Matrix<float, 3, 2>  _S648 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S649;
    (&_S649)->primal_0 = _S453;
    (&_S649)->differential_0 = _S648;
    s_bwd_prop_mul_2(&_S647, &_S649, _S645);
    Matrix<float, 2, 3>  _S650 = transpose_2(_S649.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S651;
    (&_S651)->primal_0 = J_11;
    (&_S651)->differential_0 = _S646;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S652;
    (&_S652)->primal_0 = _S434;
    (&_S652)->differential_0 = _S568;
    s_bwd_prop_mul_3(&_S651, &_S652, _S647.differential_0);
    Matrix<float, 2, 3>  _S653 = _S650 + _S651.differential_0;
    float _S654 = _S451 * _S653.rows[int(1)].z;
    float s_diff_ty_T_0 = _S450 * (rz2_6 * _S653.rows[int(1)].z);
    float _S655 = fy_14 * _S653.rows[int(1)].y;
    float _S656 = _S449 * _S653.rows[int(0)].z;
    float s_diff_tx_T_0 = _S448 * (rz2_6 * _S653.rows[int(0)].z);
    float _S657 = fx_14 * _S653.rows[int(0)].x;
    float _S658 = mean_c_10.z * s_diff_ty_T_0;
    float _S659 = _S447 * s_diff_ty_T_0;
    DiffPair_float_0 _S660;
    (&_S660)->primal_0 = lim_y_pos_0;
    (&_S660)->differential_0 = 0.0f;
    DiffPair_float_0 _S661;
    (&_S661)->primal_0 = _S446;
    (&_S661)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S660, &_S661, _S658);
    DiffPair_float_0 _S662;
    (&_S662)->primal_0 = _S444;
    (&_S662)->differential_0 = 0.0f;
    DiffPair_float_0 _S663;
    (&_S663)->primal_0 = _S445;
    (&_S663)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S662, &_S663, _S661.differential_0);
    float _S664 = mean_c_10.y * _S663.differential_0;
    float _S665 = rz_6 * _S663.differential_0;
    float _S666 = mean_c_10.z * s_diff_tx_T_0;
    float _S667 = _S443 * s_diff_tx_T_0;
    DiffPair_float_0 _S668;
    (&_S668)->primal_0 = lim_x_pos_0;
    (&_S668)->differential_0 = 0.0f;
    DiffPair_float_0 _S669;
    (&_S669)->primal_0 = _S442;
    (&_S669)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S668, &_S669, _S666);
    DiffPair_float_0 _S670;
    (&_S670)->primal_0 = _S440;
    (&_S670)->differential_0 = 0.0f;
    DiffPair_float_0 _S671;
    (&_S671)->primal_0 = _S441;
    (&_S671)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S670, &_S671, _S669.differential_0);
    float _S672 = rz_6 * (_S654 + _S656);
    float _S673 = _S659 + _S667 + - ((_S640 + _S642 + _S655 + _S657 + _S664 + mean_c_10.x * _S671.differential_0 + _S672 + _S672) / _S439);
    float _S674 = _S641 + _S665;
    float _S675 = _S643 + rz_6 * _S671.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S676;
    (&_S676)->primal_0 = _S432;
    (&_S676)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S677;
    (&_S677)->primal_0 = _S433;
    (&_S677)->differential_0 = _S568;
    s_bwd_prop_mul_4(&_S676, &_S677, _S652.differential_0);
    Matrix<float, 3, 3>  _S678 = transpose_0(_S677.differential_0 + _S571.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S679;
    (&_S679)->primal_0 = R_14;
    (&_S679)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S680;
    (&_S680)->primal_0 = _S431;
    (&_S680)->differential_0 = _S568;
    s_bwd_prop_mul_4(&_S679, &_S680, _S676.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S681;
    (&_S681)->primal_0 = _S429;
    (&_S681)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S682;
    (&_S682)->primal_0 = _S430;
    (&_S682)->differential_0 = _S568;
    s_bwd_prop_mul_4(&_S681, &_S682, _S680.differential_0);
    Matrix<float, 3, 3>  _S683 = _S681.differential_0 + transpose_0(_S682.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S684;
    (&_S684)->primal_0 = _S428;
    (&_S684)->differential_0 = _S568;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S685;
    (&_S685)->primal_0 = S_0;
    (&_S685)->differential_0 = _S568;
    s_bwd_prop_mul_4(&_S684, &_S685, _S683);
    Matrix<float, 3, 3>  _S686 = transpose_0(_S684.differential_0);
    float _S687 = 2.0f * - _S686.rows[int(2)].z;
    float _S688 = 2.0f * _S686.rows[int(2)].y;
    float _S689 = 2.0f * _S686.rows[int(2)].x;
    float _S690 = 2.0f * _S686.rows[int(1)].z;
    float _S691 = 2.0f * - _S686.rows[int(1)].y;
    float _S692 = 2.0f * _S686.rows[int(1)].x;
    float _S693 = 2.0f * _S686.rows[int(0)].z;
    float _S694 = 2.0f * _S686.rows[int(0)].y;
    float _S695 = 2.0f * - _S686.rows[int(0)].x;
    float _S696 = - _S692 + _S694;
    float _S697 = _S689 + - _S693;
    float _S698 = - _S688 + _S690;
    float _S699 = _S688 + _S690;
    float _S700 = _S689 + _S693;
    float _S701 = _S692 + _S694;
    float _S702 = quat_13.w * (_S691 + _S695);
    float _S703 = quat_13.z * (_S687 + _S695);
    float _S704 = quat_13.y * (_S687 + _S691);
    float _S705 = quat_13.x * _S696 + quat_13.z * _S699 + quat_13.y * _S700 + _S702 + _S702;
    float _S706 = quat_13.x * _S697 + quat_13.w * _S699 + quat_13.y * _S701 + _S703 + _S703;
    float _S707 = quat_13.x * _S698 + quat_13.w * _S700 + quat_13.z * _S701 + _S704 + _S704;
    float _S708 = quat_13.w * _S696 + quat_13.z * _S697 + quat_13.y * _S698;
    float3  _S709 = _S506;
    *&((&_S709)->z) = _S685.differential_0.rows[int(2)].z;
    *&((&_S709)->y) = _S685.differential_0.rows[int(1)].y;
    *&((&_S709)->x) = _S685.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S710;
    (&_S710)->primal_0 = scale_12;
    (&_S710)->differential_0 = _S506;
    s_bwd_prop_exp_1(&_S710, _S709);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S711 = _S710;
    float3  _S712 = _S506;
    *&((&_S712)->z) = _S673;
    *&((&_S712)->y) = _S674;
    *&((&_S712)->x) = _S675;
    float3  _S713 = _S608 + _S712;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S714;
    (&_S714)->primal_0 = R_14;
    (&_S714)->differential_0 = _S568;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S715;
    (&_S715)->primal_0 = mean_10;
    (&_S715)->differential_0 = _S506;
    s_bwd_prop_mul_1(&_S714, &_S715, _S713);
    float3  _S716 = _S713 + _S572.differential_0;
    Matrix<float, 3, 3>  _S717 = _S678 + _S679.differential_0 + _S714.differential_0;
    float4  _S718 = make_float4 (0.0f);
    *&((&_S718)->w) = _S705;
    *&((&_S718)->z) = _S706;
    *&((&_S718)->y) = _S707;
    *&((&_S718)->x) = _S708;
    float4  _S719 = _S718;
    float3  _S720 = _S715.differential_0 + _S566;
    *v_mean_0 = _S720;
    *v_quat_0 = _S719;
    *v_scale_0 = _S711.differential_0;
    *v_in_opacity_0 = _S615;
    (*v_sh_coeffs_0)[int(0)] = _S592;
    (*v_sh_coeffs_0)[int(1)] = _S593;
    (*v_sh_coeffs_0)[int(2)] = _S594;
    (*v_sh_coeffs_0)[int(3)] = _S595;
    (*v_sh_coeffs_0)[int(4)] = _S596;
    (*v_sh_coeffs_0)[int(5)] = _S597;
    (*v_sh_coeffs_0)[int(6)] = _S598;
    (*v_sh_coeffs_0)[int(7)] = _S599;
    (*v_sh_coeffs_0)[int(8)] = _S600;
    (*v_sh_coeffs_0)[int(9)] = _S601;
    (*v_sh_coeffs_0)[int(10)] = _S602;
    (*v_sh_coeffs_0)[int(11)] = _S603;
    (*v_sh_coeffs_0)[int(12)] = _S604;
    (*v_sh_coeffs_0)[int(13)] = _S605;
    (*v_sh_coeffs_0)[int(14)] = _S606;
    (*v_sh_coeffs_0)[int(15)] = _S607;
    *v_R_1 = _S717;
    *v_t_1 = _S716;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S721;
    DiffPair_float_0 _S722;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S723;
    DiffPair_float_0 _S724;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S725;
    DiffPair_float_0 _S726;
    DiffPair_float_0 _S727;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S728;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S729;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S730;
};

inline __device__ CameraDistortion_0 s_primal_ctx_CameraDistortion_x24init_0(float4  dpradial_coeffs_0, float2  dptangential_coeffs_0, float2  dpthin_prism_coeffs_0)
{
    CameraDistortion_0 _S731 = { dpradial_coeffs_0, dptangential_coeffs_0, dpthin_prism_coeffs_0 };
    return _S731;
}

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S732, float _S733)
{
    return s_primal_ctx_atan2_0(_S732, _S733);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S734;
    DiffPair_float_0 _S735;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S736 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S734 = _S736;
    _s_diff_ctx_2->_S735 = _S736;
    (&_s_diff_ctx_2->_S734)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S734)->differential_0 = 0.0f;
    (&_s_diff_ctx_2->_S735)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S735)->differential_0 = 0.0f;
    DiffPair_float_0 _S737 = *dpdpy_0;
    _s_diff_ctx_2->_S734 = *dpdpy_0;
    DiffPair_float_0 _S738 = *dpdpx_0;
    _s_diff_ctx_2->_S735 = *dpdpx_0;
    float _S739 = _S738.primal_0 * _S738.primal_0 + _S737.primal_0 * _S737.primal_0;
    float _S740 = - _S737.primal_0 / _S739 * dpdOut_0;
    float _S741 = _S738.primal_0 / _S739 * dpdOut_0;
    dpdpy_0->primal_0 = _S737.primal_0;
    dpdpy_0->differential_0 = _S741;
    dpdpx_0->primal_0 = _S738.primal_0;
    dpdpx_0->differential_0 = _S740;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S742, DiffPair_float_0 * _S743, float _S744, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S745 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S721 = _S745;
    _s_diff_ctx_3->_S722 = _S745;
    (&_s_diff_ctx_3->_S721)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S721)->differential_0 = 0.0f;
    (&_s_diff_ctx_3->_S722)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S722)->differential_0 = 0.0f;
    DiffPair_float_0 _S746 = *_S742;
    _s_diff_ctx_3->_S721 = *_S742;
    DiffPair_float_0 _S747 = *_S743;
    _s_diff_ctx_3->_S722 = *_S743;
    DiffPair_float_0 _S748 = _S746;
    DiffPair_float_0 _S749 = _S747;
    s_bwd_prop_d_atan2_Intermediates_0 _S750;
    (&_S750)->_S734 = _S745;
    (&_S750)->_S735 = _S745;
    s_primal_ctx_d_atan2_0(&_S748, &_S749, _S744, &_S750);
    *_S742 = _S748;
    *_S743 = _S749;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S751;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S752;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S753;
    DiffPair_float_0 _S754;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S755;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S756;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S757 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S756 = _S757;
    (&_s_diff_ctx_4->_S756)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S756)->differential_0 = 0.0f;
    DiffPair_float_0 _S758 = *dpdpx_1;
    _s_diff_ctx_4->_S756 = *dpdpx_1;
    float _S759 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S758.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S758.primal_0;
    dpdpx_1->differential_0 = _S759;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S760, float _S761, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S762 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S752 = _S762;
    (&_s_diff_ctx_5->_S752)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S752)->differential_0 = 0.0f;
    DiffPair_float_0 _S763 = *_S760;
    _s_diff_ctx_5->_S752 = *_S760;
    DiffPair_float_0 _S764 = _S763;
    s_bwd_prop_d_sqrt_Intermediates_0 _S765;
    (&_S765)->_S756 = _S762;
    s_primal_ctx_d_sqrt_0(&_S764, _S761, &_S765);
    *_S760 = _S764;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_6)
{
    float2  _S766 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S767 = { _S766, _S766 };
    DiffPair_float_0 _S768 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S769 = { _S768 };
    _s_diff_ctx_6->_S753 = _S767;
    _s_diff_ctx_6->_S754 = _S768;
    _s_diff_ctx_6->_S755 = _S769;
    (&_s_diff_ctx_6->_S753)->primal_0 = _S766;
    (&_s_diff_ctx_6->_S753)->differential_0 = _S766;
    (&_s_diff_ctx_6->_S754)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S754)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S770 = *dpdpx_2;
    _s_diff_ctx_6->_S753 = *dpdpx_2;
    float _S771 = _S770.primal_0.x;
    float _S772 = _S770.primal_0.y;
    DiffPair_float_0 _S773;
    (&_S773)->primal_0 = _S771 * _S771 + _S772 * _S772;
    (&_S773)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S773, dp_s_dOut_0, &_s_diff_ctx_6->_S755);
    _s_diff_ctx_6->_S754 = _S773;
    float _S774 = _S770.primal_0.y * _S773.differential_0;
    float _S775 = _S774 + _S774;
    float _S776 = _S770.primal_0.x * _S773.differential_0;
    float _S777 = _S776 + _S776;
    float2  _S778 = _S766;
    *&((&_S778)->y) = _S775;
    *&((&_S778)->x) = _S777;
    dpdpx_2->primal_0 = _S770.primal_0;
    dpdpx_2->differential_0 = _S778;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S779, float _S780, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_7)
{
    float2  _S781 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S782 = { _S781, _S781 };
    _s_diff_ctx_7->_S751 = _S782;
    (&_s_diff_ctx_7->_S751)->primal_0 = _S781;
    (&_s_diff_ctx_7->_S751)->differential_0 = _S781;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S783 = *_S779;
    _s_diff_ctx_7->_S751 = *_S779;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S784 = _S783;
    DiffPair_float_0 _S785 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S786 = { _S785 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S787;
    (&_S787)->_S753 = _S782;
    (&_S787)->_S754 = _S785;
    (&_S787)->_S755 = _S786;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S784, _S780, &_S787);
    *_S779 = _S784;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, float4  dpradial_coeffs_1, float2  dptangential_coeffs_1, float2  dpthin_prism_coeffs_1, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_float_0 _S788 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S789 = { _S788, _S788 };
    _s_diff_ctx_8->_S723 = _S788;
    _s_diff_ctx_8->_S724 = _S788;
    _s_diff_ctx_8->_S725 = _S789;
    _s_diff_ctx_8->_S726 = _S788;
    _s_diff_ctx_8->_S727 = _S788;
    _s_diff_ctx_8->_S728 = _S789;
    (&_s_diff_ctx_8->_S723)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S723)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S724)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S724)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S726)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S726)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S727)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S727)->differential_0 = 0.0f;
    float2  _S790 = make_float2 (0.0f);
    CameraDistortion_0 _S791 = s_primal_ctx_CameraDistortion_x24init_0(dpradial_coeffs_1, dptangential_coeffs_1, dpthin_prism_coeffs_1);
    float2  _S792 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S793 = length_0(_S792);
    float _S794 = dpmean3d_0.z;
    float _S795 = s_primal_ctx_atan2_0(_S793, _S794);
    float k_1;
    if(_S795 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S795 * _S795 / 3.0f) / _S794;
    }
    else
    {
        k_1 = _S795 / _S793;
    }
    float2  _S796 = _S792 * make_float2 (k_1);
    float k1_1 = _S791.radial_coeffs_0.x;
    float k2_1 = _S791.radial_coeffs_0.y;
    float k3_1 = _S791.radial_coeffs_0.z;
    float k4_2 = _S791.radial_coeffs_0.w;
    float p1_2 = _S791.tangential_coeffs_0.x;
    float p2_2 = _S791.tangential_coeffs_0.y;
    float sx1_2 = _S791.thin_prism_coeffs_0.x;
    float sy1_2 = _S791.thin_prism_coeffs_0.y;
    float u_4 = _S796.x;
    float v_4 = _S796.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float _S797 = 2.0f * p1_2;
    float _S798 = 2.0f * p2_2;
    float2  _S799 = _S796 * make_float2 (1.0f + r2_4 * (k1_1 + r2_4 * (k2_1 + r2_4 * (k3_1 + r2_4 * k4_2)))) + make_float2 (_S797 * u_4 * v_4 + p2_2 * (r2_4 + 2.0f * u_4 * u_4) + sx1_2 * r2_4, _S798 * u_4 * v_4 + p1_2 * (r2_4 + 2.0f * v_4 * v_4) + sy1_2 * r2_4);
    float2  _S800 = make_float2 (dpfx_0 * _S799.x + dpcx_0, dpfy_0 * _S799.y + dpcy_0);
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (0.0f);
    float _S801 = s_primal_ctx_s_primal_ctx_atan2_0(_S793, _S794);
    bool _S802 = _S801 < 0.00100000004749745f;
    float _S803;
    float _S804;
    float _S805;
    if(_S802)
    {
        float _S806 = 1.0f - _S801 * _S801 / 3.0f;
        float _S807 = _S794 * _S794;
        k_1 = _S806 / _S794;
        _S803 = 0.0f;
        _S804 = _S807;
        _S805 = _S806;
    }
    else
    {
        float _S808 = _S793 * _S793;
        k_1 = _S801 / _S793;
        _S803 = _S808;
        _S804 = 0.0f;
        _S805 = 0.0f;
    }
    float2  _S809 = make_float2 (k_1);
    float2  _S810 = _S792 * make_float2 (k_1);
    float u_5 = _S810.x;
    float v_5 = _S810.y;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float _S811 = k3_1 + r2_5 * k4_2;
    float _S812 = k2_1 + r2_5 * _S811;
    float _S813 = k1_1 + r2_5 * _S812;
    float2  _S814 = make_float2 (dpfx_0, 0.0f);
    float2  _S815 = _S810 * _S814;
    float _S816 = p2_2 * dpfx_0;
    float _S817 = _S815.x + _S815.y;
    float _S818 = r2_5 * _S817;
    float _S819 = r2_5 * _S818;
    float _S820 = sx1_2 * dpfx_0 + _S816 + _S813 * _S817 + _S812 * _S818 + _S811 * _S819 + k4_2 * (r2_5 * _S819);
    float _S821 = v_5 * _S820;
    float _S822 = u_5 * _S820;
    float2  _S823 = make_float2 (1.0f + r2_5 * _S813) * _S814 + make_float2 (2.0f * u_5 * _S816 + 2.0f * (u_5 * _S816) + _S797 * (v_5 * dpfx_0) + _S822 + _S822, _S797 * u_5 * dpfx_0 + _S821 + _S821);
    float2  _S824 = _S792 * _S823;
    float2  _S825 = _S809 * _S823;
    float _S826 = _S824.x + _S824.y;
    if(_S802)
    {
        float _S827 = _S826 / _S804;
        float _S828 = _S805 * - _S827;
        float _S829 = _S801 * (0.3333333432674408f * - (_S794 * _S827));
        k_1 = _S829 + _S829;
        _S803 = _S828;
        _S804 = 0.0f;
    }
    else
    {
        float _S830 = _S826 / _S803;
        float _S831 = _S801 * - _S830;
        k_1 = _S793 * _S830;
        _S803 = 0.0f;
        _S804 = _S831;
    }
    DiffPair_float_0 _S832;
    (&_S832)->primal_0 = _S793;
    (&_S832)->differential_0 = 0.0f;
    DiffPair_float_0 _S833;
    (&_S833)->primal_0 = _S794;
    (&_S833)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S832, &_S833, k_1, &_s_diff_ctx_8->_S725);
    _s_diff_ctx_8->_S723 = _S832;
    _s_diff_ctx_8->_S724 = _S833;
    float _S834 = _S833.differential_0 + _S803;
    float _S835 = _S832.differential_0 + _S804;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S836;
    (&_S836)->primal_0 = _S792;
    (&_S836)->differential_0 = _S790;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S837 = { _S790, _S790 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S838;
    (&_S838)->_S751 = _S837;
    s_primal_ctx_s_bwd_length_impl_0(&_S836, _S835, &_S838);
    float2  _S839 = _S836.differential_0 + _S825;
    float3  _S840 = make_float3 (_S839.x, _S839.y, _S834);
    Matrix<float, 2, 3>  _S841 = J_12;
    _S841[int(0)] = _S840;
    if(_S802)
    {
        float _S842 = 1.0f - _S801 * _S801 / 3.0f;
        float _S843 = _S794 * _S794;
        k_1 = _S842 / _S794;
        _S803 = 0.0f;
        _S804 = _S843;
        _S805 = _S842;
    }
    else
    {
        float _S844 = _S793 * _S793;
        k_1 = _S801 / _S793;
        _S803 = _S844;
        _S804 = 0.0f;
        _S805 = 0.0f;
    }
    float2  _S845 = make_float2 (k_1);
    float2  _S846 = _S792 * make_float2 (k_1);
    float u_6 = _S846.x;
    float v_6 = _S846.y;
    float r2_6 = u_6 * u_6 + v_6 * v_6;
    float _S847 = k3_1 + r2_6 * k4_2;
    float _S848 = k2_1 + r2_6 * _S847;
    float _S849 = k1_1 + r2_6 * _S848;
    float2  _S850 = make_float2 (0.0f, dpfy_0);
    float2  _S851 = _S846 * _S850;
    float _S852 = p1_2 * dpfy_0;
    float _S853 = _S851.x + _S851.y;
    float _S854 = r2_6 * _S853;
    float _S855 = r2_6 * _S854;
    float _S856 = sy1_2 * dpfy_0 + _S852 + _S849 * _S853 + _S848 * _S854 + _S847 * _S855 + k4_2 * (r2_6 * _S855);
    float _S857 = v_6 * _S856;
    float _S858 = u_6 * _S856;
    float2  _S859 = make_float2 (1.0f + r2_6 * _S849) * _S850 + make_float2 (_S798 * (v_6 * dpfy_0) + _S858 + _S858, 2.0f * v_6 * _S852 + 2.0f * (v_6 * _S852) + _S798 * u_6 * dpfy_0 + _S857 + _S857);
    float2  _S860 = _S792 * _S859;
    float2  _S861 = _S845 * _S859;
    float _S862 = _S860.x + _S860.y;
    if(_S802)
    {
        float _S863 = _S862 / _S804;
        float _S864 = _S805 * - _S863;
        float _S865 = _S801 * (0.3333333432674408f * - (_S794 * _S863));
        k_1 = _S865 + _S865;
        _S803 = _S864;
        _S804 = 0.0f;
    }
    else
    {
        float _S866 = _S862 / _S803;
        float _S867 = _S801 * - _S866;
        k_1 = _S793 * _S866;
        _S803 = 0.0f;
        _S804 = _S867;
    }
    DiffPair_float_0 _S868;
    (&_S868)->primal_0 = _S793;
    (&_S868)->differential_0 = 0.0f;
    DiffPair_float_0 _S869;
    (&_S869)->primal_0 = _S794;
    (&_S869)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S868, &_S869, k_1, &_s_diff_ctx_8->_S728);
    _s_diff_ctx_8->_S726 = _S868;
    _s_diff_ctx_8->_S727 = _S869;
    float _S870 = _S869.differential_0 + _S803;
    float _S871 = _S868.differential_0 + _S804;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S872;
    (&_S872)->primal_0 = _S792;
    (&_S872)->differential_0 = _S790;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S873;
    (&_S873)->_S751 = _S837;
    s_primal_ctx_s_bwd_length_impl_0(&_S872, _S871, &_S873);
    float2  _S874 = _S872.differential_0 + _S861;
    float3  _S875 = make_float3 (_S874.x, _S874.y, _S870);
    _S841[int(1)] = _S875;
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S841, dpcov3d_0), transpose_1(_S841));
    *dpmean2d_0 = _S800;
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
    DiffPair_1 _S876 = *dpdpx_3;
    float _S877 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_9->_S756)->primal_0);
    float _S878 = s_primal_ctx_sqrt_0(_S877);
    float _S879 = 0.5f / _S878 * (*dpdpx_3).differential_0.differential_0;
    float _S880 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S878 * _S878));
    DiffPair_float_0 _S881;
    (&_S881)->primal_0 = _S877;
    (&_S881)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S881, _S880);
    DiffPair_float_0 _S882;
    (&_S882)->primal_0 = 1.00000001168609742e-07f;
    (&_S882)->differential_0 = 0.0f;
    DiffPair_float_0 _S883;
    (&_S883)->primal_0 = (&_s_diff_ctx_9->_S756)->primal_0;
    (&_S883)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S882, &_S883, _S881.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S883.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S879;
    dpdpx_3->primal_0 = _S876.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S884, DiffPair_float_0 * _S885, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_1 _S886 = *_S884;
    DiffPair_float_0 _S887 = _s_diff_ctx_10->_S752;
    DiffPair_float_0 _S888 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S889;
    (&_S889)->_S756 = _S888;
    s_primal_ctx_d_sqrt_0(&_S887, (*_S885).primal_0, &_S889);
    DiffPair_float_0 _S890 = { (*_S884).differential_0.primal_0, (*_S884).differential_0.differential_0 };
    DiffPair_1 _S891;
    (&_S891)->primal_0 = _s_diff_ctx_10->_S752;
    (&_S891)->differential_0 = _S890;
    DiffPair_float_0 _S892;
    (&_S892)->primal_0 = (*_S885).primal_0;
    (&_S892)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S893 = _S889;
    s_bwd_prop_d_sqrt_0(&_S891, &_S892, &_S893);
    DiffPair_float_0 _S894 = { _S891.differential_0.primal_0, _S891.differential_0.differential_0 };
    _S885->primal_0 = (*_S885).primal_0;
    _S885->differential_0 = _S892.differential_0;
    _S884->primal_0 = _S886.primal_0;
    _S884->differential_0 = _S894;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S895, float _s_dOut_3)
{
    DiffPair_float_0 _S896;
    (&_S896)->primal_0 = (*_S895).primal_0;
    (&_S896)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S896, _s_dOut_3);
    _S895->primal_0 = (*_S895).primal_0;
    _S895->differential_0 = _S896.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_0 _S897 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_11->_S753)->primal_0)->x) * *&((&(&_s_diff_ctx_11->_S753)->primal_0)->x) + *&((&(&_s_diff_ctx_11->_S753)->primal_0)->y) * *&((&(&_s_diff_ctx_11->_S753)->primal_0)->y);
    DiffPair_float_0 _S898 = { len_0, 0.0f };
    float2  _S899 = make_float2 (0.0f);
    float _S900 = (*dpdpx_5).differential_0.differential_0.x;
    float _S901 = _S900 + _S900;
    float _S902 = (&_s_diff_ctx_11->_S754)->differential_0 * _S901;
    float _S903 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S904 = (&_s_diff_ctx_11->_S754)->differential_0 * _S903;
    DiffPair_float_0 _S905 = { 0.0f, *&((&(&_s_diff_ctx_11->_S753)->primal_0)->x) * _S901 + *&((&(&_s_diff_ctx_11->_S753)->primal_0)->y) * _S903 };
    DiffPair_1 _S906;
    (&_S906)->primal_0 = _S898;
    (&_S906)->differential_0 = _S905;
    DiffPair_float_0 _S907;
    (&_S907)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S907)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S906, &_S907, &_s_diff_ctx_11->_S755);
    DiffPair_float_0 _S908;
    (&_S908)->primal_0 = len_0;
    (&_S908)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S908, 0.0f);
    float _S909 = _S906.differential_0.primal_0 + _S908.differential_0;
    float _S910 = *&((&(&_s_diff_ctx_11->_S753)->primal_0)->y) * _S909;
    float _S911 = _S904 + _S910 + _S910;
    float _S912 = *&((&(&_s_diff_ctx_11->_S753)->primal_0)->x) * _S909;
    float _S913 = _S902 + _S912 + _S912;
    float2  _S914 = _S899;
    *&((&_S914)->y) = _S911;
    *&((&_S914)->x) = _S913;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S897.differential_0.primal_0 + _S914, _S899 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S907.differential_0;
    dpdpx_5->primal_0 = _S897.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S915 = (*dpdpx_7).primal_0.x;
    float _S916 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S917;
    (&_S917)->primal_0 = _S915 * _S915 + _S916 * _S916;
    (&_S917)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S917, _s_dOut_4);
    float _S918 = (*dpdpx_7).primal_0.y * _S917.differential_0;
    float _S919 = _S918 + _S918;
    float _S920 = (*dpdpx_7).primal_0.x * _S917.differential_0;
    float _S921 = _S920 + _S920;
    float2  _S922 = make_float2 (0.0f);
    *&((&_S922)->y) = _S919;
    *&((&_S922)->x) = _S921;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S922;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S923, DiffPair_float_0 * _S924, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_0 _S925 = *_S923;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S926 = _s_diff_ctx_12->_S751;
    float2  _S927 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S928 = { _S927, _S927 };
    DiffPair_float_0 _S929 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S930 = { _S929 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S931;
    (&_S931)->_S753 = _S928;
    (&_S931)->_S754 = _S929;
    (&_S931)->_S755 = _S930;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S926, (*_S924).primal_0, &_S931);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S932 = { (*_S923).differential_0.primal_0, (*_S923).differential_0.differential_0 };
    DiffPair_0 _S933;
    (&_S933)->primal_0 = _s_diff_ctx_12->_S751;
    (&_S933)->differential_0 = _S932;
    DiffPair_float_0 _S934;
    (&_S934)->primal_0 = (*_S924).primal_0;
    (&_S934)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S935 = _S931;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S933, &_S934, &_S935);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S936;
    (&_S936)->primal_0 = (&_s_diff_ctx_12->_S751)->primal_0;
    (&_S936)->differential_0 = _S927;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S936, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S937 = { _S933.differential_0.primal_0 + _S936.differential_0, _S933.differential_0.differential_0 };
    _S924->primal_0 = (*_S924).primal_0;
    _S924->differential_0 = _S934.differential_0;
    _S923->primal_0 = _S925.primal_0;
    _S923->differential_0 = _S937;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_1 _S938 = *dpdpy_1;
    DiffPair_1 _S939 = *dpdpx_8;
    float _S940 = - (&_s_diff_ctx_13->_S734)->primal_0;
    float _S941 = (&_s_diff_ctx_13->_S735)->primal_0 * (&_s_diff_ctx_13->_S735)->primal_0 + (&_s_diff_ctx_13->_S734)->primal_0 * (&_s_diff_ctx_13->_S734)->primal_0;
    float _S942 = _S941 * _S941;
    float _S943 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S942;
    float _S944 = (&_s_diff_ctx_13->_S735)->primal_0 * - _S943;
    float _S945 = (&_s_diff_ctx_13->_S734)->primal_0 * _S944;
    float _S946 = (&_s_diff_ctx_13->_S735)->primal_0 * _S944;
    float _S947 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S942;
    float _S948 = _S940 * - _S947;
    float _S949 = (&_s_diff_ctx_13->_S734)->primal_0 * _S948;
    float _S950 = (&_s_diff_ctx_13->_S735)->primal_0 * _S948;
    DiffPair_float_0 dpdpx_9 = { _S950 + _S950 + ((*dpdpx_8).differential_0.primal_0 + (_S946 + _S946 + _S941 * _S943)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S945 + _S945 + (*dpdpy_1).differential_0.primal_0 + _S949 + _S949 + - (_S941 * _S947), 0.0f };
    float _S951 = (&_s_diff_ctx_13->_S735)->primal_0 / _S941 * (*dpdpy_1).differential_0.differential_0 + _S940 / _S941 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S951;
    dpdpy_1->primal_0 = _S938.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S939.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S952, DiffPair_1 * _S953, DiffPair_float_0 * _S954, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_1 _S955 = *_S952;
    DiffPair_1 _S956 = *_S953;
    DiffPair_float_0 _S957 = _s_diff_ctx_14->_S721;
    DiffPair_float_0 _S958 = _s_diff_ctx_14->_S722;
    DiffPair_float_0 _S959 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S960;
    (&_S960)->_S734 = _S959;
    (&_S960)->_S735 = _S959;
    s_primal_ctx_d_atan2_0(&_S957, &_S958, (*_S954).primal_0, &_S960);
    DiffPair_float_0 _S961 = { (*_S953).differential_0.primal_0, (*_S953).differential_0.differential_0 };
    DiffPair_float_0 _S962 = { (*_S952).differential_0.primal_0, (*_S952).differential_0.differential_0 };
    DiffPair_1 _S963;
    (&_S963)->primal_0 = _s_diff_ctx_14->_S721;
    (&_S963)->differential_0 = _S962;
    DiffPair_1 _S964;
    (&_S964)->primal_0 = _s_diff_ctx_14->_S722;
    (&_S964)->differential_0 = _S961;
    DiffPair_float_0 _S965;
    (&_S965)->primal_0 = (*_S954).primal_0;
    (&_S965)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S966 = _S960;
    s_bwd_prop_d_atan2_0(&_S963, &_S964, &_S965, &_S966);
    DiffPair_float_0 _S967 = { _S964.differential_0.primal_0, _S964.differential_0.differential_0 };
    DiffPair_float_0 _S968 = { _S963.differential_0.primal_0, _S963.differential_0.differential_0 };
    _S954->primal_0 = (*_S954).primal_0;
    _S954->differential_0 = _S965.differential_0;
    _S952->primal_0 = _S955.primal_0;
    _S952->differential_0 = _S968;
    _S953->primal_0 = _S956.primal_0;
    _S953->differential_0 = _S967;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S969, DiffPair_float_0 * _S970, float _s_dOut_5)
{
    DiffPair_float_0 _S971;
    (&_S971)->primal_0 = (*_S969).primal_0;
    (&_S971)->differential_0 = 0.0f;
    DiffPair_float_0 _S972;
    (&_S972)->primal_0 = (*_S970).primal_0;
    (&_S972)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S971, &_S972, _s_dOut_5);
    _S970->primal_0 = (*_S970).primal_0;
    _S970->differential_0 = _S972.differential_0;
    _S969->primal_0 = (*_S969).primal_0;
    _S969->differential_0 = _S971.differential_0;
    return;
}

inline __device__ void s_bwd_prop_CameraDistortion_x24init_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_2, CameraDistortion_0 * _s_dOut_6)
{
    float2  _S973 = _s_dOut_6->thin_prism_coeffs_0;
    dpthin_prism_coeffs_2->primal_0 = (*dpthin_prism_coeffs_2).primal_0;
    dpthin_prism_coeffs_2->differential_0 = _S973;
    float2  _S974 = _s_dOut_6->tangential_coeffs_0;
    dptangential_coeffs_2->primal_0 = (*dptangential_coeffs_2).primal_0;
    dptangential_coeffs_2->differential_0 = _S974;
    float4  _S975 = _s_dOut_6->radial_coeffs_0;
    dpradial_coeffs_2->primal_0 = (*dpradial_coeffs_2).primal_0;
    dpradial_coeffs_2->differential_0 = _S975;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_3, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S976 = *dpcov3d_1;
    DiffPair_float_0 _S977 = *dpfx_1;
    DiffPair_float_0 _S978 = *dpfy_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S979 = *dpradial_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S980 = *dptangential_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S981 = *dpthin_prism_coeffs_3;
    float2  _S982 = make_float2 (0.0f);
    CameraDistortion_0 _S983 = s_primal_ctx_CameraDistortion_x24init_0((*dpradial_coeffs_3).primal_0, (*dptangential_coeffs_3).primal_0, (*dpthin_prism_coeffs_3).primal_0);
    float2  _S984 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S985 = length_0(_S984);
    float _S986 = (*dpmean3d_1).primal_0.z;
    float _S987 = s_primal_ctx_atan2_0(_S985, _S986);
    bool _S988 = _S987 < 0.00100000004749745f;
    float k_2;
    float _S989;
    float _S990;
    float _S991;
    if(_S988)
    {
        float _S992 = 1.0f - _S987 * _S987 / 3.0f;
        float _S993 = _S986 * _S986;
        k_2 = _S992 / _S986;
        _S989 = 0.0f;
        _S990 = _S993;
        _S991 = _S992;
    }
    else
    {
        float _S994 = _S985 * _S985;
        k_2 = _S987 / _S985;
        _S989 = _S994;
        _S990 = 0.0f;
        _S991 = 0.0f;
    }
    float2  _S995 = make_float2 (k_2);
    float2  _S996 = _S984 * make_float2 (k_2);
    float k1_2 = _S983.radial_coeffs_0.x;
    float k2_2 = _S983.radial_coeffs_0.y;
    float k3_2 = _S983.radial_coeffs_0.z;
    float k4_3 = _S983.radial_coeffs_0.w;
    float p1_3 = _S983.tangential_coeffs_0.x;
    float p2_3 = _S983.tangential_coeffs_0.y;
    float sx1_3 = _S983.thin_prism_coeffs_0.x;
    float sy1_3 = _S983.thin_prism_coeffs_0.y;
    float u_7 = _S996.x;
    float v_7 = _S996.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S997 = k3_2 + r2_7 * k4_3;
    float _S998 = k2_2 + r2_7 * _S997;
    float _S999 = k1_2 + r2_7 * _S998;
    float radial_1 = 1.0f + r2_7 * _S999;
    float2  _S1000 = make_float2 (radial_1);
    float _S1001 = 2.0f * p1_3;
    float _S1002 = _S1001 * u_7;
    float _S1003 = 2.0f * u_7;
    float _S1004 = r2_7 + _S1003 * u_7;
    float _S1005 = 2.0f * p2_3;
    float _S1006 = _S1005 * u_7;
    float _S1007 = 2.0f * v_7;
    float _S1008 = r2_7 + _S1007 * v_7;
    float2  _S1009 = _S996 * make_float2 (radial_1) + make_float2 (_S1002 * v_7 + p2_3 * _S1004 + sx1_3 * r2_7, _S1006 * v_7 + p1_3 * _S1008 + sy1_3 * r2_7);
    float _S1010 = _S1009.x;
    float _S1011 = _S1009.y;
    Matrix<float, 2, 3>  J_13 = makeMatrix<float, 2, 3> (0.0f);
    float _S1012 = s_primal_ctx_s_primal_ctx_atan2_0(_S985, _S986);
    bool _S1013 = _S1012 < 0.00100000004749745f;
    float _S1014;
    float _S1015;
    float _S1016;
    if(_S1013)
    {
        float _S1017 = 1.0f - _S1012 * _S1012 / 3.0f;
        float _S1018 = _S986 * _S986;
        k_2 = _S1017 / _S986;
        _S1014 = 0.0f;
        _S1015 = _S1018;
        _S1016 = _S1017;
    }
    else
    {
        float _S1019 = _S985 * _S985;
        k_2 = _S1012 / _S985;
        _S1014 = _S1019;
        _S1015 = 0.0f;
        _S1016 = 0.0f;
    }
    float2  _S1020 = make_float2 (k_2);
    float2  _S1021 = _S984 * make_float2 (k_2);
    float u_8 = _S1021.x;
    float v_8 = _S1021.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S1022 = k3_2 + r2_8 * k4_3;
    float _S1023 = k2_2 + r2_8 * _S1022;
    float _S1024 = k1_2 + r2_8 * _S1023;
    float2  _S1025 = make_float2 (1.0f + r2_8 * _S1024);
    float _S1026 = _S1001 * u_8;
    float _S1027 = 2.0f * u_8;
    float2  _S1028 = make_float2 (_S977.primal_0, 0.0f);
    float2  _S1029 = _S1021 * _S1028;
    float _S1030 = p2_3 * _S977.primal_0;
    float _S1031 = v_8 * _S977.primal_0;
    float _S1032 = _S1029.x + _S1029.y;
    float _S1033 = r2_8 * _S1032;
    float _S1034 = r2_8 * _S1033;
    float _S1035 = r2_8 * _S1034;
    float _S1036 = sx1_3 * _S977.primal_0 + _S1030 + _S1024 * _S1032 + _S1023 * _S1033 + _S1022 * _S1034 + k4_3 * _S1035;
    float _S1037 = v_8 * _S1036;
    float _S1038 = u_8 * _S1036;
    float2  _S1039 = _S1025 * _S1028 + make_float2 (_S1027 * _S1030 + 2.0f * (u_8 * _S1030) + _S1001 * _S1031 + _S1038 + _S1038, _S1026 * _S977.primal_0 + _S1037 + _S1037);
    float2  _S1040 = _S984 * _S1039;
    float2  _S1041 = _S1020 * _S1039;
    float _S1042 = _S1040.x + _S1040.y;
    float k_3;
    float _S1043;
    float _S1044;
    float _S1045;
    float _S1046;
    float _S1047;
    float _S1048;
    float _S1049;
    float _S1050;
    if(_S1013)
    {
        float _S1051 = _S1042 / _S1015;
        float _S1052 = _S1015 * _S1015;
        float _S1053 = - _S1051;
        float _S1054 = _S1016 * _S1053;
        float _S1055 = 0.3333333432674408f * - (_S986 * _S1051);
        float _S1056 = _S1012 * _S1055;
        k_2 = _S1056 + _S1056;
        k_3 = _S1054;
        _S1043 = 0.0f;
        _S1044 = 0.0f;
        _S1045 = 0.0f;
        _S1046 = 0.0f;
        _S1047 = _S1055;
        _S1048 = _S1051;
        _S1049 = _S1053;
        _S1050 = _S1052;
    }
    else
    {
        float _S1057 = _S1042 / _S1014;
        float _S1058 = _S1014 * _S1014;
        float _S1059 = - _S1057;
        float _S1060 = _S1012 * _S1059;
        k_2 = _S985 * _S1057;
        k_3 = 0.0f;
        _S1043 = _S1060;
        _S1044 = _S1057;
        _S1045 = _S1059;
        _S1046 = _S1058;
        _S1047 = 0.0f;
        _S1048 = 0.0f;
        _S1049 = 0.0f;
        _S1050 = 0.0f;
    }
    DiffPair_float_0 _S1061 = { _S985, 0.0f };
    DiffPair_float_0 _S1062 = { _S986, 0.0f };
    float _S1063 = (&_s_diff_ctx_15->_S724)->differential_0 + k_3;
    float _S1064 = (&_s_diff_ctx_15->_S723)->differential_0 + _S1043;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1065 = { _S984, _S982 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1066;
    (&_S1066)->primal_0 = _S984;
    (&_S1066)->differential_0 = _S982;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1067 = { _S982, _S982 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1068;
    (&_S1068)->_S751 = _S1067;
    s_primal_ctx_s_bwd_length_impl_0(&_S1066, _S1064, &_S1068);
    float2  _S1069 = _S1066.differential_0 + _S1041;
    float3  _S1070 = make_float3 (_S1069.x, _S1069.y, _S1063);
    Matrix<float, 2, 3>  _S1071 = J_13;
    _S1071[int(0)] = _S1070;
    float _S1072;
    float _S1073;
    if(_S1013)
    {
        float _S1074 = 1.0f - _S1012 * _S1012 / 3.0f;
        float _S1075 = _S986 * _S986;
        k_3 = _S1074 / _S986;
        _S1043 = 0.0f;
        _S1072 = _S1075;
        _S1073 = _S1074;
    }
    else
    {
        float _S1076 = _S985 * _S985;
        k_3 = _S1012 / _S985;
        _S1043 = _S1076;
        _S1072 = 0.0f;
        _S1073 = 0.0f;
    }
    float2  _S1077 = make_float2 (k_3);
    float2  _S1078 = _S984 * make_float2 (k_3);
    float u_9 = _S1078.x;
    float v_9 = _S1078.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S1079 = k3_2 + r2_9 * k4_3;
    float _S1080 = k2_2 + r2_9 * _S1079;
    float _S1081 = k1_2 + r2_9 * _S1080;
    float2  _S1082 = make_float2 (1.0f + r2_9 * _S1081);
    float _S1083 = _S1005 * u_9;
    float _S1084 = 2.0f * v_9;
    float2  _S1085 = make_float2 (0.0f, _S978.primal_0);
    float2  _S1086 = _S1078 * _S1085;
    float _S1087 = p1_3 * _S978.primal_0;
    float _S1088 = v_9 * _S978.primal_0;
    float _S1089 = _S1086.x + _S1086.y;
    float _S1090 = r2_9 * _S1089;
    float _S1091 = r2_9 * _S1090;
    float _S1092 = r2_9 * _S1091;
    float _S1093 = sy1_3 * _S978.primal_0 + _S1087 + _S1081 * _S1089 + _S1080 * _S1090 + _S1079 * _S1091 + k4_3 * _S1092;
    float _S1094 = v_9 * _S1093;
    float _S1095 = u_9 * _S1093;
    float2  _S1096 = _S1082 * _S1085 + make_float2 (_S1005 * _S1088 + _S1095 + _S1095, _S1084 * _S1087 + 2.0f * (v_9 * _S1087) + _S1083 * _S978.primal_0 + _S1094 + _S1094);
    float2  _S1097 = _S984 * _S1096;
    float2  _S1098 = _S1077 * _S1096;
    float _S1099 = _S1097.x + _S1097.y;
    float _S1100;
    float _S1101;
    float _S1102;
    float _S1103;
    float _S1104;
    float _S1105;
    float _S1106;
    float _S1107;
    float _S1108;
    if(_S1013)
    {
        float _S1109 = _S1099 / _S1072;
        float _S1110 = _S1072 * _S1072;
        float _S1111 = - _S1109;
        float _S1112 = _S1073 * _S1111;
        float _S1113 = 0.3333333432674408f * - (_S986 * _S1109);
        float _S1114 = _S1012 * _S1113;
        k_3 = _S1114 + _S1114;
        _S1100 = _S1112;
        _S1101 = 0.0f;
        _S1102 = 0.0f;
        _S1103 = 0.0f;
        _S1104 = 0.0f;
        _S1105 = _S1113;
        _S1106 = _S1109;
        _S1107 = _S1111;
        _S1108 = _S1110;
    }
    else
    {
        float _S1115 = _S1099 / _S1043;
        float _S1116 = _S1043 * _S1043;
        float _S1117 = - _S1115;
        float _S1118 = _S1012 * _S1117;
        k_3 = _S985 * _S1115;
        _S1100 = 0.0f;
        _S1101 = _S1118;
        _S1102 = _S1115;
        _S1103 = _S1117;
        _S1104 = _S1116;
        _S1105 = 0.0f;
        _S1106 = 0.0f;
        _S1107 = 0.0f;
        _S1108 = 0.0f;
    }
    float _S1119 = (&_s_diff_ctx_15->_S727)->differential_0 + _S1100;
    float _S1120 = (&_s_diff_ctx_15->_S726)->differential_0 + _S1101;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1121;
    (&_S1121)->primal_0 = _S984;
    (&_S1121)->differential_0 = _S982;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1122;
    (&_S1122)->_S751 = _S1067;
    s_primal_ctx_s_bwd_length_impl_0(&_S1121, _S1120, &_S1122);
    float2  _S1123 = _S1121.differential_0 + _S1098;
    float3  _S1124 = make_float3 (_S1123.x, _S1123.y, _S1119);
    _S1071[int(1)] = _S1124;
    Matrix<float, 3, 2>  _S1125 = transpose_1(_S1071);
    CameraDistortion_0 _S1126 = CameraDistortion_x24_syn_dzero_0();
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1127;
    (&_S1127)->primal_0 = s_primal_ctx_mul_3(_S1071, _S976.primal_0);
    (&_S1127)->differential_0 = J_13;
    Matrix<float, 3, 2>  _S1128 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1129;
    (&_S1129)->primal_0 = _S1125;
    (&_S1129)->differential_0 = _S1128;
    s_bwd_prop_mul_2(&_S1127, &_S1129, dpcov2d_1);
    Matrix<float, 2, 3>  _S1130 = transpose_2(_S1129.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1131;
    (&_S1131)->primal_0 = _S1071;
    (&_S1131)->differential_0 = J_13;
    Matrix<float, 3, 3>  _S1132 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1133;
    (&_S1133)->primal_0 = _S976.primal_0;
    (&_S1133)->differential_0 = _S1132;
    s_bwd_prop_mul_3(&_S1131, &_S1133, _S1127.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1134 = _S1133;
    Matrix<float, 2, 3>  _S1135 = _S1130 + _S1131.differential_0;
    float2  _S1136 = _S982;
    *&((&_S1136)->y) = _S1135.rows[int(1)].y;
    *&((&_S1136)->x) = _S1135.rows[int(1)].x;
    float2  _S1137 = _S1136;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1138 = { _S982, _S1136 };
    DiffPair_0 _S1139;
    (&_S1139)->primal_0 = _S1065;
    (&_S1139)->differential_0 = _S1138;
    DiffPair_float_0 _S1140;
    (&_S1140)->primal_0 = _S1120;
    (&_S1140)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1141 = _S1122;
    s_bwd_prop_s_bwd_length_impl_0(&_S1139, &_S1140, &_S1141);
    DiffPair_0 _S1142 = _S1139;
    DiffPair_float_0 _S1143 = _S1140;
    DiffPair_float_0 _S1144 = { 0.0f, _S1135.rows[int(1)].z };
    DiffPair_float_0 _S1145 = { 0.0f, _S1140.differential_0 };
    DiffPair_1 _S1146;
    (&_S1146)->primal_0 = _S1061;
    (&_S1146)->differential_0 = _S1145;
    DiffPair_1 _S1147;
    (&_S1147)->primal_0 = _S1062;
    (&_S1147)->differential_0 = _S1144;
    DiffPair_float_0 _S1148;
    (&_S1148)->primal_0 = k_3;
    (&_S1148)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1146, &_S1147, &_S1148, &_s_diff_ctx_15->_S728);
    DiffPair_1 _S1149 = _S1146;
    DiffPair_1 _S1150 = _S1147;
    DiffPair_float_0 _S1151 = _S1148;
    if(_S1013)
    {
        float _S1152 = _S1151.differential_0 + _S1151.differential_0;
        float _S1153 = _S1105 * _S1152;
        float _S1154 = - (0.3333333432674408f * (_S1012 * _S1152));
        float _S1155 = _S1107 * _S1135.rows[int(1)].z;
        float _S1156 = (_S986 * _S1154 + - (_S1073 * _S1135.rows[int(1)].z)) / _S1108;
        float _S1157 = _S1099 * - _S1156;
        float _S1158 = _S1106 * _S1154 + _S1150.differential_0.primal_0;
        k_3 = _S1072 * _S1156;
        _S1100 = 0.0f;
        _S1101 = _S1157;
        _S1102 = _S1155;
        _S1103 = _S1149.differential_0.primal_0;
        _S1104 = _S1153;
        _S1105 = _S1158;
    }
    else
    {
        float _S1159 = _S1103 * _S1143.differential_0;
        float _S1160 = (_S985 * _S1151.differential_0 + - (_S1012 * _S1143.differential_0)) / _S1104;
        float _S1161 = _S1099 * - _S1160;
        float _S1162 = _S1102 * _S1151.differential_0 + _S1149.differential_0.primal_0;
        k_3 = _S1043 * _S1160;
        _S1100 = _S1161;
        _S1101 = 0.0f;
        _S1102 = 0.0f;
        _S1103 = _S1162;
        _S1104 = _S1159;
        _S1105 = _S1150.differential_0.primal_0;
    }
    float2  _S1163 = _S1077 * _S1137;
    float2  _S1164 = _S1096 * _S1137;
    float2  _S1165 = _S982;
    *&((&_S1165)->y) = k_3;
    *&((&_S1165)->x) = k_3;
    float2  _S1166 = _S1096 * _S1165;
    float2  _S1167 = _S1163 + _S984 * _S1165;
    float _S1168 = _S1167.x;
    float _S1169 = _S1168 + _S1168;
    float _S1170 = _S1093 * _S1169;
    float _S1171 = _S1167.y + _S1167.y;
    float _S1172 = _S1093 * _S1171;
    float _S1173 = u_9 * _S1169 + v_9 * _S1171;
    float _S1174 = k4_3 * _S1173;
    float _S1175 = _S1092 * _S1173;
    float _S1176 = _S1091 * _S1173;
    float _S1177 = _S1091 * _S1174;
    float _S1178 = _S1090 * _S1173;
    float _S1179 = _S1079 * _S1173 + r2_9 * _S1174;
    float _S1180 = _S1090 * _S1179;
    float _S1181 = _S1089 * _S1173;
    float _S1182 = _S1080 * _S1173 + r2_9 * _S1179;
    float _S1183 = _S1089 * _S1182;
    float _S1184 = _S1081 * _S1173 + r2_9 * _S1182;
    float _S1185 = _S1005 * _S1167.x;
    float _S1186 = _S1088 * _S1167.x;
    float _S1187 = v_9 * _S1185;
    float _S1188 = _S978.primal_0 * _S1185;
    float _S1189 = _S1083 * _S1167.y;
    float _S1190 = _S978.primal_0 * _S1167.y;
    float _S1191 = 2.0f * _S1167.y;
    float _S1192 = _S1087 * _S1191;
    float _S1193 = _S1087 * _S1167.y;
    float _S1194 = _S1173 + v_9 * _S1191 + _S1084 * _S1167.y;
    float _S1195 = p1_3 * _S1194;
    float _S1196 = _S978.primal_0 * _S1194;
    float _S1197 = sy1_3 * _S1173;
    float _S1198 = _S978.primal_0 * _S1173;
    float2  _S1199 = _S1082 * _S1167;
    float2  _S1200 = _S1085 * _S1167;
    float2  _S1201 = _S982;
    *&((&_S1201)->y) = _S1184;
    *&((&_S1201)->x) = _S1184;
    float _S1202 = _S1200.x + _S1200.y;
    float _S1203 = _S1181 + r2_9 * _S1202;
    float _S1204 = _S1178 + r2_9 * _S1203;
    float _S1205 = _S1176 + r2_9 * _S1204;
    float _S1206 = _S1177 + _S1180 + _S1183 + _S1081 * _S1202 + _S1080 * _S1203 + _S1079 * _S1204 + k4_3 * _S1205;
    float _S1207 = v_9 * _S1206;
    float _S1208 = u_9 * _S1206;
    float2  _S1209 = _S1166 + _S1142.differential_0.primal_0;
    float2  _S1210 = _S1085 * _S1201 + make_float2 (_S1170 + _S1005 * _S1190 + _S1208 + _S1208, _S1172 + _S1188 + _S1192 + 2.0f * _S1193 + _S1207 + _S1207);
    float _S1211 = _S1187 + _S1189 + _S1195 + _S1197 + (_S1199 + _S1078 * _S1201).y;
    float _S1212 = _S1186 + u_9 * _S1190;
    float _S1213 = _S1175 + r2_9 * _S1205;
    float2  _S1214 = _S984 * _S1210;
    float _S1215 = _S1164.x + _S1164.y + _S1214.x + _S1214.y;
    float2  _S1216 = _S1077 * _S1210 + _S1209;
    if(_S1013)
    {
        float _S1217 = _S986 * _S1101;
        float _S1218 = _S1215 / _S1072;
        float _S1219 = _S1012 * (0.3333333432674408f * - (_S1102 + _S986 * _S1218));
        float _S1220 = _S1217 + _S1217 + _S1073 * - _S1218 + _S1105;
        k_3 = _S1219 + _S1219 + _S1104;
        _S1043 = _S1220;
        _S1072 = _S1103;
    }
    else
    {
        float _S1221 = _S985 * _S1100;
        float _S1222 = _S1215 / _S1043;
        float _S1223 = _S1221 + _S1221 + _S1012 * - _S1222 + _S1103;
        k_3 = _S985 * _S1222 + _S1104;
        _S1043 = _S1105;
        _S1072 = _S1223;
    }
    DiffPair_float_0 _S1224;
    (&_S1224)->primal_0 = _S985;
    (&_S1224)->differential_0 = 0.0f;
    DiffPair_float_0 _S1225;
    (&_S1225)->primal_0 = _S986;
    (&_S1225)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1224, &_S1225, k_3);
    float _S1226 = _S1225.differential_0 + _S1043;
    float _S1227 = _S1224.differential_0 + _S1072;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1228;
    (&_S1228)->primal_0 = _S984;
    (&_S1228)->differential_0 = _S982;
    s_bwd_length_impl_1(&_S1228, _S1227);
    float2  _S1229 = _S1228.differential_0 + _S1216;
    float3  _S1230 = make_float3 (_S1229.x, _S1229.y, _S1226);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1231;
    (&_S1231)->primal_0 = _S984;
    (&_S1231)->differential_0 = _S982;
    s_bwd_length_impl_1(&_S1231, 0.0f);
    float3  _S1232 = _S1230 + make_float3 (_S1231.differential_0.x, _S1231.differential_0.y, 0.0f);
    float2  _S1233 = _S982;
    *&((&_S1233)->y) = _S1135.rows[int(0)].y;
    *&((&_S1233)->x) = _S1135.rows[int(0)].x;
    float2  _S1234 = _S1233;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1235 = { _S982, _S1233 };
    DiffPair_0 _S1236;
    (&_S1236)->primal_0 = _S1065;
    (&_S1236)->differential_0 = _S1235;
    DiffPair_float_0 _S1237;
    (&_S1237)->primal_0 = _S1064;
    (&_S1237)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1238 = _S1068;
    s_bwd_prop_s_bwd_length_impl_0(&_S1236, &_S1237, &_S1238);
    DiffPair_0 _S1239 = _S1236;
    DiffPair_float_0 _S1240 = _S1237;
    DiffPair_float_0 _S1241 = { 0.0f, _S1135.rows[int(0)].z };
    DiffPair_float_0 _S1242 = { 0.0f, _S1237.differential_0 };
    DiffPair_1 _S1243;
    (&_S1243)->primal_0 = _S1061;
    (&_S1243)->differential_0 = _S1242;
    DiffPair_1 _S1244;
    (&_S1244)->primal_0 = _S1062;
    (&_S1244)->differential_0 = _S1241;
    DiffPair_float_0 _S1245;
    (&_S1245)->primal_0 = k_2;
    (&_S1245)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1243, &_S1244, &_S1245, &_s_diff_ctx_15->_S725);
    DiffPair_1 _S1246 = _S1243;
    DiffPair_1 _S1247 = _S1244;
    DiffPair_float_0 _S1248 = _S1245;
    if(_S1013)
    {
        float _S1249 = _S1248.differential_0 + _S1248.differential_0;
        float _S1250 = _S1047 * _S1249;
        float _S1251 = - (0.3333333432674408f * (_S1012 * _S1249));
        float _S1252 = _S1049 * _S1135.rows[int(0)].z;
        float _S1253 = (_S986 * _S1251 + - (_S1016 * _S1135.rows[int(0)].z)) / _S1050;
        float _S1254 = _S1042 * - _S1253;
        float _S1255 = _S1048 * _S1251 + _S1247.differential_0.primal_0;
        k_2 = _S1015 * _S1253;
        k_3 = 0.0f;
        _S1043 = _S1254;
        _S1044 = _S1252;
        _S1045 = _S1246.differential_0.primal_0;
        _S1046 = _S1250;
        _S1047 = _S1255;
    }
    else
    {
        float _S1256 = _S1045 * _S1240.differential_0;
        float _S1257 = (_S985 * _S1248.differential_0 + - (_S1012 * _S1240.differential_0)) / _S1046;
        float _S1258 = _S1042 * - _S1257;
        float _S1259 = _S1044 * _S1248.differential_0 + _S1246.differential_0.primal_0;
        k_2 = _S1014 * _S1257;
        k_3 = _S1258;
        _S1043 = 0.0f;
        _S1044 = 0.0f;
        _S1045 = _S1259;
        _S1046 = _S1256;
        _S1047 = _S1247.differential_0.primal_0;
    }
    float2  _S1260 = _S1020 * _S1234;
    float2  _S1261 = _S1039 * _S1234;
    float2  _S1262 = _S982;
    *&((&_S1262)->y) = k_2;
    *&((&_S1262)->x) = k_2;
    float2  _S1263 = _S1039 * _S1262;
    float2  _S1264 = _S1260 + _S984 * _S1262;
    float _S1265 = _S1264.x;
    float _S1266 = _S1265 + _S1265;
    float _S1267 = _S1036 * _S1266;
    float _S1268 = _S1264.y + _S1264.y;
    float _S1269 = _S1036 * _S1268;
    float _S1270 = u_8 * _S1266 + v_8 * _S1268;
    float _S1271 = k4_3 * _S1270;
    float _S1272 = _S1035 * _S1270;
    float _S1273 = _S1034 * _S1270;
    float _S1274 = _S1034 * _S1271;
    float _S1275 = _S1033 * _S1270;
    float _S1276 = _S1022 * _S1270 + r2_8 * _S1271;
    float _S1277 = _S1033 * _S1276;
    float _S1278 = _S1032 * _S1270;
    float _S1279 = _S1023 * _S1270 + r2_8 * _S1276;
    float _S1280 = _S1032 * _S1279;
    float _S1281 = _S1024 * _S1270 + r2_8 * _S1279;
    float _S1282 = _S1001 * _S1264.x;
    float _S1283 = _S1031 * _S1264.x;
    float _S1284 = v_8 * _S1282;
    float _S1285 = _S977.primal_0 * _S1282;
    float _S1286 = _S1026 * _S1264.y;
    float _S1287 = _S977.primal_0 * _S1264.y;
    float _S1288 = 2.0f * _S1264.x;
    float _S1289 = _S1030 * _S1288;
    float _S1290 = _S1030 * _S1264.x;
    float _S1291 = _S1270 + u_8 * _S1288 + _S1027 * _S1264.x;
    float _S1292 = p2_3 * _S1291;
    float _S1293 = _S977.primal_0 * _S1291;
    float _S1294 = sx1_3 * _S1270;
    float _S1295 = _S977.primal_0 * _S1270;
    float2  _S1296 = _S1025 * _S1264;
    float2  _S1297 = _S1028 * _S1264;
    float2  _S1298 = _S982;
    *&((&_S1298)->y) = _S1281;
    *&((&_S1298)->x) = _S1281;
    float _S1299 = _S1297.x + _S1297.y;
    float _S1300 = _S1278 + r2_8 * _S1299;
    float _S1301 = _S1275 + r2_8 * _S1300;
    float _S1302 = _S1273 + r2_8 * _S1301;
    float _S1303 = _S1274 + _S1277 + _S1280 + _S1024 * _S1299 + _S1023 * _S1300 + _S1022 * _S1301 + k4_3 * _S1302;
    float _S1304 = v_8 * _S1303;
    float _S1305 = u_8 * _S1303;
    float2  _S1306 = _S1263 + _S1239.differential_0.primal_0;
    float _S1307 = _S1283 + u_8 * _S1287;
    float _S1308 = _S1300 + _S1203;
    float _S1309 = _S1302 + _S1205;
    float _S1310 = _S1284 + _S1286 + _S1292 + _S1294 + (_S1296 + _S1021 * _S1298).x;
    float2  _S1311 = _S1028 * _S1298 + make_float2 (_S1267 + _S1289 + 2.0f * _S1290 + _S1001 * _S1287 + _S1305 + _S1305, _S1269 + _S1285 + _S1304 + _S1304);
    float _S1312 = _S1301 + _S1204;
    float _S1313 = _S1272 + r2_8 * _S1302 + _S1213;
    float2  _S1314 = _S984 * _S1311;
    float _S1315 = _S1261.x + _S1261.y + _S1314.x + _S1314.y;
    float2  _S1316 = _S1020 * _S1311 + _S1306;
    if(_S1013)
    {
        float _S1317 = _S986 * _S1043;
        float _S1318 = _S1315 / _S1015;
        float _S1319 = _S1012 * (0.3333333432674408f * - (_S1044 + _S986 * _S1318));
        float _S1320 = _S1317 + _S1317 + _S1016 * - _S1318 + _S1047;
        k_2 = _S1319 + _S1319 + _S1046;
        _S1014 = _S1320;
        _S1015 = _S1045;
    }
    else
    {
        float _S1321 = _S985 * k_3;
        float _S1322 = _S1315 / _S1014;
        float _S1323 = _S1321 + _S1321 + _S1012 * - _S1322 + _S1045;
        k_2 = _S985 * _S1322 + _S1046;
        _S1014 = _S1047;
        _S1015 = _S1323;
    }
    DiffPair_float_0 _S1324;
    (&_S1324)->primal_0 = _S985;
    (&_S1324)->differential_0 = 0.0f;
    DiffPair_float_0 _S1325;
    (&_S1325)->primal_0 = _S986;
    (&_S1325)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1324, &_S1325, k_2);
    float _S1326 = _S1325.differential_0 + _S1014;
    float _S1327 = _S1324.differential_0 + _S1015;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1328;
    (&_S1328)->primal_0 = _S984;
    (&_S1328)->differential_0 = _S982;
    s_bwd_length_impl_1(&_S1328, _S1327);
    float2  _S1329 = _S1328.differential_0 + _S1316;
    float3  _S1330 = make_float3 (_S1329.x, _S1329.y, _S1326);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1331;
    (&_S1331)->primal_0 = _S984;
    (&_S1331)->differential_0 = _S982;
    s_bwd_length_impl_1(&_S1331, 0.0f);
    float _S1332 = _S978.primal_0 * dpmean2d_1.y;
    float _S1333 = _S977.primal_0 * dpmean2d_1.x;
    float2  _S1334 = make_float2 (_S1333, _S1332);
    float2  _S1335 = _S996 * _S1334;
    float _S1336 = p1_3 * _S1332;
    float _S1337 = v_7 * _S1332;
    float _S1338 = p2_3 * _S1333;
    float _S1339 = v_7 * _S1333;
    float _S1340 = _S1335.x + _S1335.y;
    float _S1341 = r2_7 * _S1340;
    float _S1342 = r2_7 * _S1341;
    float _S1343 = r2_7 * _S1342;
    float _S1344 = sy1_3 * _S1332 + _S1336 + sx1_3 * _S1333 + _S1338 + _S999 * _S1340 + _S998 * _S1341 + _S997 * _S1342 + k4_3 * _S1343;
    float _S1345 = v_7 * _S1344;
    float _S1346 = u_7 * _S1344;
    float2  _S1347 = make_float2 (r2_7 * _S1333 + _S1295, r2_7 * _S1332 + _S1198);
    float2  _S1348 = make_float2 (_S1008 * _S1332 + 2.0f * (u_7 * _S1339 + _S1307) + _S1196, 2.0f * (u_7 * _S1337 + _S1212) + _S1004 * _S1333 + _S1293);
    float4  _S1349 = make_float4 (_S1341 + _S1308, _S1342 + _S1312, _S1343 + _S1309, r2_7 * _S1343 + _S1313);
    float3  _S1350 = _S1330 + make_float3 (_S1331.differential_0.x, _S1331.differential_0.y, 0.0f) + _S1232;
    float _S1351 = _S1010 * dpmean2d_1.x + _S1310;
    float _S1352 = _S1011 * dpmean2d_1.y + _S1211;
    float2  _S1353 = _S1000 * _S1334 + make_float2 (_S1005 * _S1337 + _S1003 * _S1338 + 2.0f * (u_7 * _S1338) + _S1001 * _S1339 + _S1346 + _S1346, _S1007 * _S1336 + 2.0f * (v_7 * _S1336) + _S1006 * _S1332 + _S1002 * _S1333 + _S1345 + _S1345);
    CameraDistortion_0 _S1354 = _S1126;
    (&_S1354)->thin_prism_coeffs_0 = _S1347;
    (&_S1354)->tangential_coeffs_0 = _S1348;
    (&_S1354)->radial_coeffs_0 = _S1349;
    CameraDistortion_0 _S1355 = _S1126;
    CameraDistortion_0 _S1356 = _S1354;
    CameraDistortion_0 _S1357 = CameraDistortion_x24_syn_dadd_0(&_S1355, &_S1356);
    float2  _S1358 = _S984 * _S1353;
    float2  _S1359 = _S995 * _S1353;
    float _S1360 = _S1358.x + _S1358.y;
    if(_S988)
    {
        float _S1361 = _S1360 / _S990;
        float _S1362 = _S991 * - _S1361;
        float _S1363 = _S987 * (0.3333333432674408f * - (_S986 * _S1361));
        k_2 = _S1363 + _S1363;
        _S989 = _S1362;
        _S990 = 0.0f;
    }
    else
    {
        float _S1364 = _S1360 / _S989;
        float _S1365 = _S987 * - _S1364;
        k_2 = _S985 * _S1364;
        _S989 = 0.0f;
        _S990 = _S1365;
    }
    DiffPair_float_0 _S1366;
    (&_S1366)->primal_0 = _S985;
    (&_S1366)->differential_0 = 0.0f;
    DiffPair_float_0 _S1367;
    (&_S1367)->primal_0 = _S986;
    (&_S1367)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1366, &_S1367, k_2);
    float _S1368 = _S1367.differential_0 + _S989;
    float _S1369 = _S1366.differential_0 + _S990;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1370;
    (&_S1370)->primal_0 = _S984;
    (&_S1370)->differential_0 = _S982;
    s_bwd_length_impl_1(&_S1370, _S1369);
    float2  _S1371 = _S1370.differential_0 + _S1359;
    float4  _S1372 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1373;
    (&_S1373)->primal_0 = _S979.primal_0;
    (&_S1373)->differential_0 = _S1372;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1374;
    (&_S1374)->primal_0 = _S980.primal_0;
    (&_S1374)->differential_0 = _S982;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1375;
    (&_S1375)->primal_0 = _S981.primal_0;
    (&_S1375)->differential_0 = _S982;
    CameraDistortion_0 _S1376 = _S1357;
    s_bwd_prop_CameraDistortion_x24init_0(&_S1373, &_S1374, &_S1375, &_S1376);
    dpthin_prism_coeffs_3->primal_0 = (*dpthin_prism_coeffs_3).primal_0;
    dpthin_prism_coeffs_3->differential_0 = _S1375.differential_0;
    dptangential_coeffs_3->primal_0 = (*dptangential_coeffs_3).primal_0;
    dptangential_coeffs_3->differential_0 = _S1374.differential_0;
    dpradial_coeffs_3->primal_0 = (*dpradial_coeffs_3).primal_0;
    dpradial_coeffs_3->differential_0 = _S1373.differential_0;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = dpmean2d_1.y;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = dpmean2d_1.x;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S1352;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S1351;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1134.differential_0;
    float3  _S1377 = _S1350 + make_float3 (_S1371.x, _S1371.y, _S1368);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1377;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_15, float fy_15, float cx_15, float cy_15, float4  radial_coeffs_18, float2  tangential_coeffs_18, float2  thin_prism_coeffs_18, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1378 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1379 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1380 = { _S1379, _S1379 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1381 = { _S1379, _S1379, _S1380, _S1379, _S1379, _S1380 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1382;
    (&_S1382)->_S729 = _S1378;
    (&_S1382)->_S730 = _S1381;
    float3  mean_c_11 = s_primal_ctx_mul_1(R_15, mean_11) + t_14;
    float3  _S1383 = s_primal_ctx_exp_0(scale_13);
    float _S1384 = quat_14.y;
    float x2_14 = _S1384 * _S1384;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_25 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S1385 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1383.x, 0.0f, 0.0f, 0.0f, _S1383.y, 0.0f, 0.0f, 0.0f, _S1383.z);
    Matrix<float, 3, 3>  _S1386 = s_primal_ctx_mul_2(_S1385, S_1);
    Matrix<float, 3, 3>  _S1387 = transpose_0(_S1386);
    Matrix<float, 3, 3>  _S1388 = s_primal_ctx_mul_2(_S1386, _S1387);
    Matrix<float, 3, 3>  _S1389 = s_primal_ctx_mul_2(R_15, _S1388);
    Matrix<float, 3, 3>  _S1390 = transpose_0(R_15);
    Matrix<float, 3, 3>  _S1391 = s_primal_ctx_mul_2(_S1389, _S1390);
    Matrix<float, 2, 2>  _S1392 = _S1378;
    float2  _S1393 = make_float2 (0.0f);
    float2  _S1394 = _S1393;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1391, fx_15, fy_15, cx_15, cy_15, radial_coeffs_18, tangential_coeffs_18, thin_prism_coeffs_18, &_S1392, &_S1394, &(&_S1382)->_S730);
    (&_S1382)->_S729 = _S1392;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1395 = _S1382;
    float _S1396 = _S1382._S729.rows[int(0)].y * _S1382._S729.rows[int(1)].x;
    float det_orig_12 = _S1382._S729.rows[int(0)].x * _S1382._S729.rows[int(1)].y - _S1396;
    float _S1397 = _S1382._S729.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1398 = _S1382._S729;
    *&(((&_S1398)->rows + (int(0)))->x) = _S1397;
    float _S1399 = _S1382._S729.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1398)->rows + (int(1)))->y) = _S1399;
    Matrix<float, 2, 2>  _S1400 = _S1398;
    Matrix<float, 2, 2>  _S1401 = _S1398;
    float det_blur_7 = _S1397 * _S1399 - _S1396;
    float _S1402 = det_orig_12 / det_blur_7;
    float _S1403 = det_blur_7 * det_blur_7;
    float _S1404 = s_primal_ctx_max_0(0.0f, _S1402);
    float _S1405 = s_primal_ctx_sqrt_0(_S1404);
    float invdet_8 = 1.0f / det_blur_7;
    float _S1406 = - _S1382._S729.rows[int(0)].y;
    float _S1407 = - _S1382._S729.rows[int(1)].x;
    float _S1408 = - in_opacity_11;
    float _S1409 = 1.0f + s_primal_ctx_exp_1(_S1408);
    float _S1410 = 1.0f / _S1409;
    float _S1411 = _S1409 * _S1409;
    float _S1412;
    if(antialiased_11)
    {
        _S1412 = _S1410 * _S1405;
    }
    else
    {
        _S1412 = _S1410;
    }
    float _S1413 = _S1412 / 0.00392156885936856f;
    float _S1414 = 2.0f * s_primal_ctx_log_0(_S1413);
    float _S1415 = s_primal_ctx_sqrt_0(_S1414);
    float _S1416 = _S1400.rows[int(0)].x;
    float _S1417 = _S1401.rows[int(1)].y;
    float _S1418 = s_primal_ctx_dot_0(mean_c_11, mean_c_11) + 9.99999997475242708e-07f;
    float3  _S1419 = mean_11 - - s_primal_ctx_mul_1(_S1390, t_14);
    float _S1420 = _S1419.x;
    float _S1421 = _S1419.y;
    float _S1422 = _S1419.z;
    float _S1423 = _S1420 * _S1420 + _S1421 * _S1421 + _S1422 * _S1422;
    float _S1424 = s_primal_ctx_sqrt_0(_S1423);
    float x_35 = _S1420 / _S1424;
    float3  _S1425 = make_float3 (x_35);
    float _S1426 = _S1424 * _S1424;
    float y_14 = _S1421 / _S1424;
    float z_11 = _S1422 / _S1424;
    float3  _S1427 = make_float3 (z_11);
    float _S1428 = - y_14;
    float3  _S1429 = make_float3 (_S1428);
    float z2_26 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_35 * x_35 - y_14 * y_14;
    float _S1430 = 2.0f * x_35;
    float fS1_11 = _S1430 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_26 - 0.31539157032966614f;
    float3  _S1431 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_35;
    float3  _S1432 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S1433 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S1434 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S1435 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S1436 = 1.86588168144226074f * z2_26 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S1436;
    float3  _S1437 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_35;
    float3  _S1438 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S1439 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S1440 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S1441 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_35 * fC1_11 - y_14 * fS1_11);
    float3  _S1442 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_35 * fS1_11 + y_14 * fC1_11);
    float3  _S1443 = make_float3 (pSH9_1);
    float3  _S1444 = make_float3 (0.0f);
    float3  _S1445 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1446;
    (&_S1446)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1428) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_35) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S1446)->differential_0 = _S1445;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1447;
    (&_S1447)->primal_0 = _S1444;
    (&_S1447)->differential_0 = _S1445;
    s_bwd_prop_max_0(&_S1446, &_S1447, v_rgb_1);
    float3  _S1448 = _S1442 * _S1446.differential_0;
    float3  _S1449 = (*sh_coeffs_11)[int(15)] * _S1446.differential_0;
    float3  _S1450 = _S1440 * _S1446.differential_0;
    float3  _S1451 = (*sh_coeffs_11)[int(14)] * _S1446.differential_0;
    float3  _S1452 = _S1438 * _S1446.differential_0;
    float3  _S1453 = (*sh_coeffs_11)[int(13)] * _S1446.differential_0;
    float3  _S1454 = _S1437 * _S1446.differential_0;
    float3  _S1455 = (*sh_coeffs_11)[int(12)] * _S1446.differential_0;
    float3  _S1456 = _S1439 * _S1446.differential_0;
    float3  _S1457 = (*sh_coeffs_11)[int(11)] * _S1446.differential_0;
    float3  _S1458 = _S1441 * _S1446.differential_0;
    float3  _S1459 = (*sh_coeffs_11)[int(10)] * _S1446.differential_0;
    float3  _S1460 = _S1443 * _S1446.differential_0;
    float3  _S1461 = (*sh_coeffs_11)[int(9)] * _S1446.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1461.x + _S1461.y + _S1461.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1449.x + _S1449.y + _S1449.z);
    float _S1462 = _S1459.x + _S1459.y + _S1459.z;
    float _S1463 = _S1451.x + _S1451.y + _S1451.z;
    float _S1464 = _S1457.x + _S1457.y + _S1457.z;
    float _S1465 = _S1453.x + _S1453.y + _S1453.z;
    float _S1466 = _S1455.x + _S1455.y + _S1455.z;
    float _S1467 = - s_diff_fC2_T_1;
    float3  _S1468 = _S1434 * _S1446.differential_0;
    float3  _S1469 = (*sh_coeffs_11)[int(8)] * _S1446.differential_0;
    float3  _S1470 = _S1432 * _S1446.differential_0;
    float3  _S1471 = (*sh_coeffs_11)[int(7)] * _S1446.differential_0;
    float3  _S1472 = _S1431 * _S1446.differential_0;
    float3  _S1473 = (*sh_coeffs_11)[int(6)] * _S1446.differential_0;
    float3  _S1474 = _S1433 * _S1446.differential_0;
    float3  _S1475 = (*sh_coeffs_11)[int(5)] * _S1446.differential_0;
    float3  _S1476 = _S1435 * _S1446.differential_0;
    float3  _S1477 = (*sh_coeffs_11)[int(4)] * _S1446.differential_0;
    float _S1478 = _S1475.x + _S1475.y + _S1475.z;
    float _S1479 = _S1471.x + _S1471.y + _S1471.z;
    float _S1480 = fTmp1B_11 * _S1462 + x_35 * s_diff_fS2_T_1 + y_14 * _S1467 + 0.54627424478530884f * (_S1477.x + _S1477.y + _S1477.z);
    float _S1481 = fTmp1B_11 * _S1463 + y_14 * s_diff_fS2_T_1 + x_35 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1469.x + _S1469.y + _S1469.z);
    float _S1482 = y_14 * - _S1481;
    float _S1483 = x_35 * _S1481;
    float _S1484 = z_11 * (1.86588168144226074f * (z_11 * _S1466) + -2.28522896766662598f * (y_14 * _S1464 + x_35 * _S1465) + 0.94617468118667603f * (_S1473.x + _S1473.y + _S1473.z));
    float3  _S1485 = make_float3 (0.48860251903533936f) * _S1446.differential_0;
    float3  _S1486 = - _S1485;
    float3  _S1487 = _S1425 * _S1486;
    float3  _S1488 = (*sh_coeffs_11)[int(3)] * _S1486;
    float3  _S1489 = _S1427 * _S1485;
    float3  _S1490 = (*sh_coeffs_11)[int(2)] * _S1485;
    float3  _S1491 = _S1429 * _S1485;
    float3  _S1492 = (*sh_coeffs_11)[int(1)] * _S1485;
    float _S1493 = (_S1436 * _S1466 + 1.44530570507049561f * (fS1_11 * _S1462 + fC1_11 * _S1463) + -1.09254848957061768f * (y_14 * _S1478 + x_35 * _S1479) + _S1484 + _S1484 + _S1490.x + _S1490.y + _S1490.z) / _S1426;
    float _S1494 = _S1424 * _S1493;
    float _S1495 = (fTmp0C_11 * _S1464 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S1467 + fTmp0B_11 * _S1478 + _S1430 * _S1480 + _S1482 + _S1482 + - (_S1492.x + _S1492.y + _S1492.z)) / _S1426;
    float _S1496 = _S1424 * _S1495;
    float _S1497 = (fTmp0C_11 * _S1465 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S1479 + 2.0f * (y_14 * _S1480) + _S1483 + _S1483 + _S1488.x + _S1488.y + _S1488.z) / _S1426;
    float _S1498 = _S1424 * _S1497;
    float _S1499 = _S1422 * - _S1493 + _S1421 * - _S1495 + _S1420 * - _S1497;
    DiffPair_float_0 _S1500;
    (&_S1500)->primal_0 = _S1423;
    (&_S1500)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1500, _S1499);
    float _S1501 = _S1422 * _S1500.differential_0;
    float _S1502 = _S1421 * _S1500.differential_0;
    float _S1503 = _S1420 * _S1500.differential_0;
    float3  _S1504 = make_float3 (0.282094806432724f) * _S1446.differential_0;
    float3  _S1505 = make_float3 (_S1498 + _S1503 + _S1503, _S1496 + _S1502 + _S1502, _S1494 + _S1501 + _S1501);
    float3  _S1506 = - - _S1505;
    Matrix<float, 3, 3>  _S1507 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1508;
    (&_S1508)->primal_0 = _S1390;
    (&_S1508)->differential_0 = _S1507;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1509;
    (&_S1509)->primal_0 = t_14;
    (&_S1509)->differential_0 = _S1445;
    s_bwd_prop_mul_1(&_S1508, &_S1509, _S1506);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1510 = _S1508;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1511 = _S1509;
    float2  _S1512 = _S1393;
    *&((&_S1512)->y) = v_conic_1.z;
    float2  _S1513 = _S1393;
    *&((&_S1513)->y) = v_conic_1.y;
    *&((&_S1513)->x) = v_conic_1.x;
    float _S1514 = 0.5f * v_depth_1;
    DiffPair_float_0 _S1515;
    (&_S1515)->primal_0 = _S1418;
    (&_S1515)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1515, _S1514);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1516;
    (&_S1516)->primal_0 = mean_c_11;
    (&_S1516)->differential_0 = _S1445;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1517;
    (&_S1517)->primal_0 = mean_c_11;
    (&_S1517)->differential_0 = _S1445;
    s_bwd_prop_dot_0(&_S1516, &_S1517, _S1515.differential_0);
    DiffPair_float_0 _S1518;
    (&_S1518)->primal_0 = _S1417;
    (&_S1518)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1518, 0.0f);
    DiffPair_float_0 _S1519;
    (&_S1519)->primal_0 = _S1416;
    (&_S1519)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1519, 0.0f);
    DiffPair_float_0 _S1520;
    (&_S1520)->primal_0 = 3.32999992370605469f;
    (&_S1520)->differential_0 = 0.0f;
    DiffPair_float_0 _S1521;
    (&_S1521)->primal_0 = _S1415;
    (&_S1521)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1520, &_S1521, 0.0f);
    DiffPair_float_0 _S1522;
    (&_S1522)->primal_0 = _S1414;
    (&_S1522)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1522, _S1521.differential_0);
    float _S1523 = 2.0f * _S1522.differential_0;
    DiffPair_float_0 _S1524;
    (&_S1524)->primal_0 = _S1413;
    (&_S1524)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1524, _S1523);
    float2  _S1525 = make_float2 (_S1519.differential_0, 0.0f);
    float _S1526 = v_opacity_1 + 254.9999847412109375f * _S1524.differential_0;
    Matrix<float, 2, 2>  _S1527 = _S1378;
    _S1527[int(1)] = _S1512;
    _S1527[int(0)] = _S1513;
    Matrix<float, 2, 2>  _S1528 = _S1527;
    FixedArray<float3 , 16>  _S1529;
    _S1529[int(0)] = _S1445;
    _S1529[int(1)] = _S1445;
    _S1529[int(2)] = _S1445;
    _S1529[int(3)] = _S1445;
    _S1529[int(4)] = _S1445;
    _S1529[int(5)] = _S1445;
    _S1529[int(6)] = _S1445;
    _S1529[int(7)] = _S1445;
    _S1529[int(8)] = _S1445;
    _S1529[int(9)] = _S1445;
    _S1529[int(10)] = _S1445;
    _S1529[int(11)] = _S1445;
    _S1529[int(12)] = _S1445;
    _S1529[int(13)] = _S1445;
    _S1529[int(14)] = _S1445;
    _S1529[int(15)] = _S1445;
    _S1529[int(7)] = _S1470;
    _S1529[int(0)] = _S1504;
    _S1529[int(1)] = _S1491;
    _S1529[int(2)] = _S1489;
    _S1529[int(3)] = _S1487;
    _S1529[int(4)] = _S1476;
    _S1529[int(5)] = _S1474;
    _S1529[int(6)] = _S1472;
    _S1529[int(15)] = _S1448;
    _S1529[int(8)] = _S1468;
    _S1529[int(9)] = _S1460;
    _S1529[int(10)] = _S1458;
    _S1529[int(11)] = _S1456;
    _S1529[int(12)] = _S1454;
    _S1529[int(13)] = _S1452;
    _S1529[int(14)] = _S1450;
    float3  _S1530 = _S1529[int(0)];
    float3  _S1531 = _S1529[int(1)];
    float3  _S1532 = _S1529[int(2)];
    float3  _S1533 = _S1529[int(3)];
    float3  _S1534 = _S1529[int(4)];
    float3  _S1535 = _S1529[int(5)];
    float3  _S1536 = _S1529[int(6)];
    float3  _S1537 = _S1529[int(7)];
    float3  _S1538 = _S1529[int(8)];
    float3  _S1539 = _S1529[int(9)];
    float3  _S1540 = _S1529[int(10)];
    float3  _S1541 = _S1529[int(11)];
    float3  _S1542 = _S1529[int(12)];
    float3  _S1543 = _S1529[int(13)];
    float3  _S1544 = _S1529[int(14)];
    float3  _S1545 = _S1529[int(15)];
    float3  _S1546 = _S1517.differential_0 + _S1516.differential_0;
    float2  _S1547 = make_float2 (0.0f, _S1518.differential_0);
    float _S1548;
    if(antialiased_11)
    {
        float _S1549 = _S1410 * _S1526;
        _S1412 = _S1405 * _S1526;
        _S1548 = _S1549;
    }
    else
    {
        _S1412 = _S1526;
        _S1548 = 0.0f;
    }
    float _S1550 = - (_S1412 / _S1411);
    DiffPair_float_0 _S1551;
    (&_S1551)->primal_0 = _S1408;
    (&_S1551)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1551, _S1550);
    float _S1552 = - _S1551.differential_0;
    float _S1553 = invdet_8 * _S1528.rows[int(1)].y;
    float _S1554 = - (invdet_8 * _S1528.rows[int(1)].x);
    float _S1555 = - (invdet_8 * _S1528.rows[int(0)].y);
    float _S1556 = invdet_8 * _S1528.rows[int(0)].x;
    float _S1557 = - ((_S1397 * _S1528.rows[int(1)].y + _S1407 * _S1528.rows[int(1)].x + _S1406 * _S1528.rows[int(0)].y + _S1399 * _S1528.rows[int(0)].x) / _S1403);
    DiffPair_float_0 _S1558;
    (&_S1558)->primal_0 = _S1404;
    (&_S1558)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1558, _S1548);
    DiffPair_float_0 _S1559;
    (&_S1559)->primal_0 = 0.0f;
    (&_S1559)->differential_0 = 0.0f;
    DiffPair_float_0 _S1560;
    (&_S1560)->primal_0 = _S1402;
    (&_S1560)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1559, &_S1560, _S1558.differential_0);
    float _S1561 = _S1560.differential_0 / _S1403;
    float s_diff_det_orig_T_1 = det_blur_7 * _S1561;
    float _S1562 = _S1557 + det_orig_12 * - _S1561;
    float _S1563 = - _S1562;
    float _S1564 = _S1397 * _S1562;
    float _S1565 = _S1399 * _S1562;
    Matrix<float, 2, 2>  _S1566 = _S1378;
    _S1566[int(1)] = _S1547;
    _S1566[int(0)] = _S1525;
    _S1398 = _S1566;
    *&(((&_S1398)->rows + (int(1)))->y) = 0.0f;
    float _S1567 = _S1556 + _S1564 + _S1566.rows[int(1)].y;
    *&(((&_S1398)->rows + (int(0)))->x) = 0.0f;
    float _S1568 = _S1553 + _S1565 + _S1566.rows[int(0)].x;
    float _S1569 = _S1563 + - s_diff_det_orig_T_1;
    float _S1570 = _S1554 + _S1395._S729.rows[int(0)].y * _S1569;
    float _S1571 = _S1555 + _S1395._S729.rows[int(1)].x * _S1569;
    float _S1572 = _S1395._S729.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S1573 = _S1567 + _S1395._S729.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S1574 = _S1393;
    *&((&_S1574)->x) = _S1570;
    *&((&_S1574)->y) = _S1573;
    float _S1575 = _S1568 + _S1572;
    float2  _S1576 = _S1393;
    *&((&_S1576)->y) = _S1571;
    *&((&_S1576)->x) = _S1575;
    Matrix<float, 2, 2>  _S1577 = _S1378;
    _S1577[int(1)] = _S1574;
    _S1577[int(0)] = _S1576;
    Matrix<float, 2, 2>  _S1578 = _S1398 + _S1577;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1579;
    (&_S1579)->primal_0 = mean_c_11;
    (&_S1579)->differential_0 = _S1445;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1580;
    (&_S1580)->primal_0 = _S1391;
    (&_S1580)->differential_0 = _S1507;
    DiffPair_float_0 _S1581;
    (&_S1581)->primal_0 = fx_15;
    (&_S1581)->differential_0 = 0.0f;
    DiffPair_float_0 _S1582;
    (&_S1582)->primal_0 = fy_15;
    (&_S1582)->differential_0 = 0.0f;
    DiffPair_float_0 _S1583;
    (&_S1583)->primal_0 = cx_15;
    (&_S1583)->differential_0 = 0.0f;
    DiffPair_float_0 _S1584;
    (&_S1584)->primal_0 = cy_15;
    (&_S1584)->differential_0 = 0.0f;
    float4  _S1585 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1586;
    (&_S1586)->primal_0 = radial_coeffs_18;
    (&_S1586)->differential_0 = _S1585;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1587;
    (&_S1587)->primal_0 = tangential_coeffs_18;
    (&_S1587)->differential_0 = _S1393;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1588;
    (&_S1588)->primal_0 = thin_prism_coeffs_18;
    (&_S1588)->differential_0 = _S1393;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1589 = _S1395._S730;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S1579, &_S1580, &_S1581, &_S1582, &_S1583, &_S1584, &_S1586, &_S1587, &_S1588, _S1578, v_mean2d_1, &_S1589);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1590;
    (&_S1590)->primal_0 = _S1389;
    (&_S1590)->differential_0 = _S1507;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1591;
    (&_S1591)->primal_0 = _S1390;
    (&_S1591)->differential_0 = _S1507;
    s_bwd_prop_mul_4(&_S1590, &_S1591, _S1580.differential_0);
    Matrix<float, 3, 3>  _S1592 = transpose_0(_S1591.differential_0 + _S1510.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1593;
    (&_S1593)->primal_0 = R_15;
    (&_S1593)->differential_0 = _S1507;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1594;
    (&_S1594)->primal_0 = _S1388;
    (&_S1594)->differential_0 = _S1507;
    s_bwd_prop_mul_4(&_S1593, &_S1594, _S1590.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1595;
    (&_S1595)->primal_0 = _S1386;
    (&_S1595)->differential_0 = _S1507;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1596;
    (&_S1596)->primal_0 = _S1387;
    (&_S1596)->differential_0 = _S1507;
    s_bwd_prop_mul_4(&_S1595, &_S1596, _S1594.differential_0);
    Matrix<float, 3, 3>  _S1597 = _S1595.differential_0 + transpose_0(_S1596.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1598;
    (&_S1598)->primal_0 = _S1385;
    (&_S1598)->differential_0 = _S1507;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1599;
    (&_S1599)->primal_0 = S_1;
    (&_S1599)->differential_0 = _S1507;
    s_bwd_prop_mul_4(&_S1598, &_S1599, _S1597);
    Matrix<float, 3, 3>  _S1600 = transpose_0(_S1598.differential_0);
    float _S1601 = 2.0f * - _S1600.rows[int(2)].z;
    float _S1602 = 2.0f * _S1600.rows[int(2)].y;
    float _S1603 = 2.0f * _S1600.rows[int(2)].x;
    float _S1604 = 2.0f * _S1600.rows[int(1)].z;
    float _S1605 = 2.0f * - _S1600.rows[int(1)].y;
    float _S1606 = 2.0f * _S1600.rows[int(1)].x;
    float _S1607 = 2.0f * _S1600.rows[int(0)].z;
    float _S1608 = 2.0f * _S1600.rows[int(0)].y;
    float _S1609 = 2.0f * - _S1600.rows[int(0)].x;
    float _S1610 = - _S1606 + _S1608;
    float _S1611 = _S1603 + - _S1607;
    float _S1612 = - _S1602 + _S1604;
    float _S1613 = _S1602 + _S1604;
    float _S1614 = _S1603 + _S1607;
    float _S1615 = _S1606 + _S1608;
    float _S1616 = quat_14.w * (_S1605 + _S1609);
    float _S1617 = quat_14.z * (_S1601 + _S1609);
    float _S1618 = quat_14.y * (_S1601 + _S1605);
    float _S1619 = quat_14.x * _S1610 + quat_14.z * _S1613 + quat_14.y * _S1614 + _S1616 + _S1616;
    float _S1620 = quat_14.x * _S1611 + quat_14.w * _S1613 + quat_14.y * _S1615 + _S1617 + _S1617;
    float _S1621 = quat_14.x * _S1612 + quat_14.w * _S1614 + quat_14.z * _S1615 + _S1618 + _S1618;
    float _S1622 = quat_14.w * _S1610 + quat_14.z * _S1611 + quat_14.y * _S1612;
    float3  _S1623 = _S1445;
    *&((&_S1623)->z) = _S1599.differential_0.rows[int(2)].z;
    *&((&_S1623)->y) = _S1599.differential_0.rows[int(1)].y;
    *&((&_S1623)->x) = _S1599.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1624;
    (&_S1624)->primal_0 = scale_13;
    (&_S1624)->differential_0 = _S1445;
    s_bwd_prop_exp_1(&_S1624, _S1623);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1625 = _S1624;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1626;
    (&_S1626)->primal_0 = mean_c_11;
    (&_S1626)->differential_0 = _S1445;
    s_bwd_length_impl_0(&_S1626, 0.0f);
    float3  _S1627 = _S1579.differential_0 + _S1626.differential_0 + _S1546;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1628;
    (&_S1628)->primal_0 = R_15;
    (&_S1628)->differential_0 = _S1507;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1629;
    (&_S1629)->primal_0 = mean_11;
    (&_S1629)->differential_0 = _S1445;
    s_bwd_prop_mul_1(&_S1628, &_S1629, _S1627);
    float3  _S1630 = _S1627 + _S1511.differential_0;
    Matrix<float, 3, 3>  _S1631 = _S1592 + _S1593.differential_0 + _S1628.differential_0;
    float4  _S1632 = _S1585;
    *&((&_S1632)->w) = _S1619;
    *&((&_S1632)->z) = _S1620;
    *&((&_S1632)->y) = _S1621;
    *&((&_S1632)->x) = _S1622;
    float4  _S1633 = _S1632;
    float3  _S1634 = _S1629.differential_0 + _S1505;
    *v_mean_1 = _S1634;
    *v_quat_1 = _S1633;
    *v_scale_1 = _S1625.differential_0;
    *v_in_opacity_1 = _S1552;
    (*v_sh_coeffs_1)[int(0)] = _S1530;
    (*v_sh_coeffs_1)[int(1)] = _S1531;
    (*v_sh_coeffs_1)[int(2)] = _S1532;
    (*v_sh_coeffs_1)[int(3)] = _S1533;
    (*v_sh_coeffs_1)[int(4)] = _S1534;
    (*v_sh_coeffs_1)[int(5)] = _S1535;
    (*v_sh_coeffs_1)[int(6)] = _S1536;
    (*v_sh_coeffs_1)[int(7)] = _S1537;
    (*v_sh_coeffs_1)[int(8)] = _S1538;
    (*v_sh_coeffs_1)[int(9)] = _S1539;
    (*v_sh_coeffs_1)[int(10)] = _S1540;
    (*v_sh_coeffs_1)[int(11)] = _S1541;
    (*v_sh_coeffs_1)[int(12)] = _S1542;
    (*v_sh_coeffs_1)[int(13)] = _S1543;
    (*v_sh_coeffs_1)[int(14)] = _S1544;
    (*v_sh_coeffs_1)[int(15)] = _S1545;
    *v_R_2 = _S1631;
    *v_t_2 = _S1630;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_16, float fy_16, float cx_16, float cy_16, float4  radial_coeffs_19, float2  tangential_coeffs_19, float2  thin_prism_coeffs_19, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_1(R_16, mean_12) + t_15;
    float3  _S1635 = s_primal_ctx_exp_0(scale_14);
    float _S1636 = quat_15.y;
    float x2_15 = _S1636 * _S1636;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_27 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S1637 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S1635.x, 0.0f, 0.0f, 0.0f, _S1635.y, 0.0f, 0.0f, 0.0f, _S1635.z);
    Matrix<float, 3, 3>  _S1638 = s_primal_ctx_mul_2(_S1637, S_2);
    Matrix<float, 3, 3>  _S1639 = transpose_0(_S1638);
    Matrix<float, 3, 3>  _S1640 = s_primal_ctx_mul_2(_S1638, _S1639);
    Matrix<float, 3, 3>  _S1641 = s_primal_ctx_mul_2(R_16, _S1640);
    Matrix<float, 3, 3>  _S1642 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S1643 = s_primal_ctx_mul_2(_S1641, _S1642);
    Matrix<float, 2, 3>  J_14 = makeMatrix<float, 2, 3> (fx_16, 0.0f, 0.0f, 0.0f, fy_16, 0.0f);
    Matrix<float, 2, 3>  _S1644 = s_primal_ctx_mul_3(J_14, _S1643);
    Matrix<float, 3, 2>  _S1645 = transpose_1(J_14);
    Matrix<float, 2, 2>  _S1646 = s_primal_ctx_mul_4(_S1644, _S1645);
    float _S1647 = _S1646.rows[int(0)].y * _S1646.rows[int(1)].x;
    float det_orig_13 = _S1646.rows[int(0)].x * _S1646.rows[int(1)].y - _S1647;
    float _S1648 = _S1646.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1649 = _S1646;
    *&(((&_S1649)->rows + (int(0)))->x) = _S1648;
    float _S1650 = _S1646.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1649)->rows + (int(1)))->y) = _S1650;
    Matrix<float, 2, 2>  _S1651 = _S1649;
    Matrix<float, 2, 2>  _S1652 = _S1649;
    float det_blur_8 = _S1648 * _S1650 - _S1647;
    float _S1653 = det_orig_13 / det_blur_8;
    float _S1654 = det_blur_8 * det_blur_8;
    float _S1655 = s_primal_ctx_max_0(0.0f, _S1653);
    float _S1656 = s_primal_ctx_sqrt_0(_S1655);
    float invdet_9 = 1.0f / det_blur_8;
    float _S1657 = - _S1646.rows[int(0)].y;
    float _S1658 = - _S1646.rows[int(1)].x;
    float _S1659 = - in_opacity_12;
    float _S1660 = 1.0f + s_primal_ctx_exp_1(_S1659);
    float _S1661 = 1.0f / _S1660;
    float _S1662 = _S1660 * _S1660;
    float _S1663;
    if(antialiased_12)
    {
        _S1663 = _S1661 * _S1656;
    }
    else
    {
        _S1663 = _S1661;
    }
    float _S1664 = _S1663 / 0.00392156885936856f;
    float _S1665 = 2.0f * s_primal_ctx_log_0(_S1664);
    float _S1666 = s_primal_ctx_sqrt_0(_S1665);
    float _S1667 = _S1651.rows[int(0)].x;
    float _S1668 = _S1652.rows[int(1)].y;
    float _S1669 = s_primal_ctx_dot_0(mean_c_12, mean_c_12) + 9.99999997475242708e-07f;
    float3  _S1670 = mean_12 - - s_primal_ctx_mul_1(_S1642, t_15);
    float _S1671 = _S1670.x;
    float _S1672 = _S1670.y;
    float _S1673 = _S1670.z;
    float _S1674 = _S1671 * _S1671 + _S1672 * _S1672 + _S1673 * _S1673;
    float _S1675 = s_primal_ctx_sqrt_0(_S1674);
    float x_36 = _S1671 / _S1675;
    float3  _S1676 = make_float3 (x_36);
    float _S1677 = _S1675 * _S1675;
    float y_15 = _S1672 / _S1675;
    float z_12 = _S1673 / _S1675;
    float3  _S1678 = make_float3 (z_12);
    float _S1679 = - y_15;
    float3  _S1680 = make_float3 (_S1679);
    float z2_28 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_36 * x_36 - y_15 * y_15;
    float _S1681 = 2.0f * x_36;
    float fS1_12 = _S1681 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_28 - 0.31539157032966614f;
    float3  _S1682 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_36;
    float3  _S1683 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S1684 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S1685 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S1686 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S1687 = 1.86588168144226074f * z2_28 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S1687;
    float3  _S1688 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_36;
    float3  _S1689 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S1690 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S1691 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S1692 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_36 * fC1_12 - y_15 * fS1_12);
    float3  _S1693 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_36 * fS1_12 + y_15 * fC1_12);
    float3  _S1694 = make_float3 (pSH9_2);
    float3  _S1695 = make_float3 (0.0f);
    float3  _S1696 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1697;
    (&_S1697)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1679) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_36) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S1697)->differential_0 = _S1696;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1698;
    (&_S1698)->primal_0 = _S1695;
    (&_S1698)->differential_0 = _S1696;
    s_bwd_prop_max_0(&_S1697, &_S1698, v_rgb_2);
    float3  _S1699 = _S1693 * _S1697.differential_0;
    float3  _S1700 = (*sh_coeffs_12)[int(15)] * _S1697.differential_0;
    float3  _S1701 = _S1691 * _S1697.differential_0;
    float3  _S1702 = (*sh_coeffs_12)[int(14)] * _S1697.differential_0;
    float3  _S1703 = _S1689 * _S1697.differential_0;
    float3  _S1704 = (*sh_coeffs_12)[int(13)] * _S1697.differential_0;
    float3  _S1705 = _S1688 * _S1697.differential_0;
    float3  _S1706 = (*sh_coeffs_12)[int(12)] * _S1697.differential_0;
    float3  _S1707 = _S1690 * _S1697.differential_0;
    float3  _S1708 = (*sh_coeffs_12)[int(11)] * _S1697.differential_0;
    float3  _S1709 = _S1692 * _S1697.differential_0;
    float3  _S1710 = (*sh_coeffs_12)[int(10)] * _S1697.differential_0;
    float3  _S1711 = _S1694 * _S1697.differential_0;
    float3  _S1712 = (*sh_coeffs_12)[int(9)] * _S1697.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S1712.x + _S1712.y + _S1712.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S1700.x + _S1700.y + _S1700.z);
    float _S1713 = _S1710.x + _S1710.y + _S1710.z;
    float _S1714 = _S1702.x + _S1702.y + _S1702.z;
    float _S1715 = _S1708.x + _S1708.y + _S1708.z;
    float _S1716 = _S1704.x + _S1704.y + _S1704.z;
    float _S1717 = _S1706.x + _S1706.y + _S1706.z;
    float _S1718 = - s_diff_fC2_T_2;
    float3  _S1719 = _S1685 * _S1697.differential_0;
    float3  _S1720 = (*sh_coeffs_12)[int(8)] * _S1697.differential_0;
    float3  _S1721 = _S1683 * _S1697.differential_0;
    float3  _S1722 = (*sh_coeffs_12)[int(7)] * _S1697.differential_0;
    float3  _S1723 = _S1682 * _S1697.differential_0;
    float3  _S1724 = (*sh_coeffs_12)[int(6)] * _S1697.differential_0;
    float3  _S1725 = _S1684 * _S1697.differential_0;
    float3  _S1726 = (*sh_coeffs_12)[int(5)] * _S1697.differential_0;
    float3  _S1727 = _S1686 * _S1697.differential_0;
    float3  _S1728 = (*sh_coeffs_12)[int(4)] * _S1697.differential_0;
    float _S1729 = _S1726.x + _S1726.y + _S1726.z;
    float _S1730 = _S1722.x + _S1722.y + _S1722.z;
    float _S1731 = fTmp1B_12 * _S1713 + x_36 * s_diff_fS2_T_2 + y_15 * _S1718 + 0.54627424478530884f * (_S1728.x + _S1728.y + _S1728.z);
    float _S1732 = fTmp1B_12 * _S1714 + y_15 * s_diff_fS2_T_2 + x_36 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S1720.x + _S1720.y + _S1720.z);
    float _S1733 = y_15 * - _S1732;
    float _S1734 = x_36 * _S1732;
    float _S1735 = z_12 * (1.86588168144226074f * (z_12 * _S1717) + -2.28522896766662598f * (y_15 * _S1715 + x_36 * _S1716) + 0.94617468118667603f * (_S1724.x + _S1724.y + _S1724.z));
    float3  _S1736 = make_float3 (0.48860251903533936f) * _S1697.differential_0;
    float3  _S1737 = - _S1736;
    float3  _S1738 = _S1676 * _S1737;
    float3  _S1739 = (*sh_coeffs_12)[int(3)] * _S1737;
    float3  _S1740 = _S1678 * _S1736;
    float3  _S1741 = (*sh_coeffs_12)[int(2)] * _S1736;
    float3  _S1742 = _S1680 * _S1736;
    float3  _S1743 = (*sh_coeffs_12)[int(1)] * _S1736;
    float _S1744 = (_S1687 * _S1717 + 1.44530570507049561f * (fS1_12 * _S1713 + fC1_12 * _S1714) + -1.09254848957061768f * (y_15 * _S1729 + x_36 * _S1730) + _S1735 + _S1735 + _S1741.x + _S1741.y + _S1741.z) / _S1677;
    float _S1745 = _S1675 * _S1744;
    float _S1746 = (fTmp0C_12 * _S1715 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S1718 + fTmp0B_12 * _S1729 + _S1681 * _S1731 + _S1733 + _S1733 + - (_S1743.x + _S1743.y + _S1743.z)) / _S1677;
    float _S1747 = _S1675 * _S1746;
    float _S1748 = (fTmp0C_12 * _S1716 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S1730 + 2.0f * (y_15 * _S1731) + _S1734 + _S1734 + _S1739.x + _S1739.y + _S1739.z) / _S1677;
    float _S1749 = _S1675 * _S1748;
    float _S1750 = _S1673 * - _S1744 + _S1672 * - _S1746 + _S1671 * - _S1748;
    DiffPair_float_0 _S1751;
    (&_S1751)->primal_0 = _S1674;
    (&_S1751)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1751, _S1750);
    float _S1752 = _S1673 * _S1751.differential_0;
    float _S1753 = _S1672 * _S1751.differential_0;
    float _S1754 = _S1671 * _S1751.differential_0;
    float3  _S1755 = make_float3 (0.282094806432724f) * _S1697.differential_0;
    float3  _S1756 = make_float3 (_S1749 + _S1754 + _S1754, _S1747 + _S1753 + _S1753, _S1745 + _S1752 + _S1752);
    float3  _S1757 = - - _S1756;
    Matrix<float, 3, 3>  _S1758 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1759;
    (&_S1759)->primal_0 = _S1642;
    (&_S1759)->differential_0 = _S1758;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1760;
    (&_S1760)->primal_0 = t_15;
    (&_S1760)->differential_0 = _S1696;
    s_bwd_prop_mul_1(&_S1759, &_S1760, _S1757);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1761 = _S1759;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1762 = _S1760;
    float2  _S1763 = make_float2 (0.0f);
    float2  _S1764 = _S1763;
    *&((&_S1764)->y) = v_conic_2.z;
    float2  _S1765 = _S1763;
    *&((&_S1765)->y) = v_conic_2.y;
    *&((&_S1765)->x) = v_conic_2.x;
    float _S1766 = 0.5f * v_depth_2;
    DiffPair_float_0 _S1767;
    (&_S1767)->primal_0 = _S1669;
    (&_S1767)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1767, _S1766);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1768;
    (&_S1768)->primal_0 = mean_c_12;
    (&_S1768)->differential_0 = _S1696;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1769;
    (&_S1769)->primal_0 = mean_c_12;
    (&_S1769)->differential_0 = _S1696;
    s_bwd_prop_dot_0(&_S1768, &_S1769, _S1767.differential_0);
    DiffPair_float_0 _S1770;
    (&_S1770)->primal_0 = _S1668;
    (&_S1770)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1770, 0.0f);
    DiffPair_float_0 _S1771;
    (&_S1771)->primal_0 = _S1667;
    (&_S1771)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1771, 0.0f);
    DiffPair_float_0 _S1772;
    (&_S1772)->primal_0 = 3.32999992370605469f;
    (&_S1772)->differential_0 = 0.0f;
    DiffPair_float_0 _S1773;
    (&_S1773)->primal_0 = _S1666;
    (&_S1773)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1772, &_S1773, 0.0f);
    DiffPair_float_0 _S1774;
    (&_S1774)->primal_0 = _S1665;
    (&_S1774)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1774, _S1773.differential_0);
    float _S1775 = 2.0f * _S1774.differential_0;
    DiffPair_float_0 _S1776;
    (&_S1776)->primal_0 = _S1664;
    (&_S1776)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1776, _S1775);
    float _S1777 = v_opacity_2 + 254.9999847412109375f * _S1776.differential_0;
    Matrix<float, 2, 2>  _S1778 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1779 = _S1778;
    _S1779[int(1)] = _S1764;
    _S1779[int(0)] = _S1765;
    Matrix<float, 2, 2>  _S1780 = _S1779;
    FixedArray<float3 , 16>  _S1781;
    _S1781[int(0)] = _S1696;
    _S1781[int(1)] = _S1696;
    _S1781[int(2)] = _S1696;
    _S1781[int(3)] = _S1696;
    _S1781[int(4)] = _S1696;
    _S1781[int(5)] = _S1696;
    _S1781[int(6)] = _S1696;
    _S1781[int(7)] = _S1696;
    _S1781[int(8)] = _S1696;
    _S1781[int(9)] = _S1696;
    _S1781[int(10)] = _S1696;
    _S1781[int(11)] = _S1696;
    _S1781[int(12)] = _S1696;
    _S1781[int(13)] = _S1696;
    _S1781[int(14)] = _S1696;
    _S1781[int(15)] = _S1696;
    _S1781[int(7)] = _S1721;
    _S1781[int(0)] = _S1755;
    _S1781[int(1)] = _S1742;
    _S1781[int(2)] = _S1740;
    _S1781[int(3)] = _S1738;
    _S1781[int(4)] = _S1727;
    _S1781[int(5)] = _S1725;
    _S1781[int(6)] = _S1723;
    _S1781[int(15)] = _S1699;
    _S1781[int(8)] = _S1719;
    _S1781[int(9)] = _S1711;
    _S1781[int(10)] = _S1709;
    _S1781[int(11)] = _S1707;
    _S1781[int(12)] = _S1705;
    _S1781[int(13)] = _S1703;
    _S1781[int(14)] = _S1701;
    float3  _S1782 = _S1781[int(0)];
    float3  _S1783 = _S1781[int(1)];
    float3  _S1784 = _S1781[int(2)];
    float3  _S1785 = _S1781[int(3)];
    float3  _S1786 = _S1781[int(4)];
    float3  _S1787 = _S1781[int(5)];
    float3  _S1788 = _S1781[int(6)];
    float3  _S1789 = _S1781[int(7)];
    float3  _S1790 = _S1781[int(8)];
    float3  _S1791 = _S1781[int(9)];
    float3  _S1792 = _S1781[int(10)];
    float3  _S1793 = _S1781[int(11)];
    float3  _S1794 = _S1781[int(12)];
    float3  _S1795 = _S1781[int(13)];
    float3  _S1796 = _S1781[int(14)];
    float3  _S1797 = _S1781[int(15)];
    float3  _S1798 = _S1769.differential_0 + _S1768.differential_0;
    float2  _S1799 = make_float2 (0.0f, _S1770.differential_0);
    float2  _S1800 = make_float2 (_S1771.differential_0, 0.0f);
    float _S1801;
    if(antialiased_12)
    {
        float _S1802 = _S1661 * _S1777;
        _S1663 = _S1656 * _S1777;
        _S1801 = _S1802;
    }
    else
    {
        _S1663 = _S1777;
        _S1801 = 0.0f;
    }
    float _S1803 = - (_S1663 / _S1662);
    DiffPair_float_0 _S1804;
    (&_S1804)->primal_0 = _S1659;
    (&_S1804)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1804, _S1803);
    float _S1805 = - _S1804.differential_0;
    float _S1806 = invdet_9 * _S1780.rows[int(1)].y;
    float _S1807 = - (invdet_9 * _S1780.rows[int(1)].x);
    float _S1808 = - (invdet_9 * _S1780.rows[int(0)].y);
    float _S1809 = invdet_9 * _S1780.rows[int(0)].x;
    float _S1810 = - ((_S1648 * _S1780.rows[int(1)].y + _S1658 * _S1780.rows[int(1)].x + _S1657 * _S1780.rows[int(0)].y + _S1650 * _S1780.rows[int(0)].x) / _S1654);
    DiffPair_float_0 _S1811;
    (&_S1811)->primal_0 = _S1655;
    (&_S1811)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1811, _S1801);
    DiffPair_float_0 _S1812;
    (&_S1812)->primal_0 = 0.0f;
    (&_S1812)->differential_0 = 0.0f;
    DiffPair_float_0 _S1813;
    (&_S1813)->primal_0 = _S1653;
    (&_S1813)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1812, &_S1813, _S1811.differential_0);
    float _S1814 = _S1813.differential_0 / _S1654;
    float s_diff_det_orig_T_2 = det_blur_8 * _S1814;
    float _S1815 = _S1810 + det_orig_13 * - _S1814;
    float _S1816 = - _S1815;
    float _S1817 = _S1648 * _S1815;
    float _S1818 = _S1650 * _S1815;
    Matrix<float, 2, 2>  _S1819 = _S1778;
    _S1819[int(1)] = _S1799;
    _S1819[int(0)] = _S1800;
    _S1649 = _S1819;
    *&(((&_S1649)->rows + (int(1)))->y) = 0.0f;
    float _S1820 = _S1809 + _S1817 + _S1819.rows[int(1)].y;
    *&(((&_S1649)->rows + (int(0)))->x) = 0.0f;
    float _S1821 = _S1806 + _S1818 + _S1819.rows[int(0)].x;
    float _S1822 = _S1816 + - s_diff_det_orig_T_2;
    float _S1823 = _S1807 + _S1646.rows[int(0)].y * _S1822;
    float _S1824 = _S1808 + _S1646.rows[int(1)].x * _S1822;
    float _S1825 = _S1646.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S1826 = _S1820 + _S1646.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S1827 = _S1763;
    *&((&_S1827)->x) = _S1823;
    *&((&_S1827)->y) = _S1826;
    float _S1828 = _S1821 + _S1825;
    float2  _S1829 = _S1763;
    *&((&_S1829)->y) = _S1824;
    *&((&_S1829)->x) = _S1828;
    float _S1830 = fy_16 * v_mean2d_2.y;
    float _S1831 = fx_16 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S1832 = _S1778;
    _S1832[int(1)] = _S1827;
    _S1832[int(0)] = _S1829;
    Matrix<float, 2, 2>  _S1833 = _S1649 + _S1832;
    Matrix<float, 2, 3>  _S1834 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1835;
    (&_S1835)->primal_0 = _S1644;
    (&_S1835)->differential_0 = _S1834;
    Matrix<float, 3, 2>  _S1836 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1837;
    (&_S1837)->primal_0 = _S1645;
    (&_S1837)->differential_0 = _S1836;
    s_bwd_prop_mul_2(&_S1835, &_S1837, _S1833);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1838;
    (&_S1838)->primal_0 = J_14;
    (&_S1838)->differential_0 = _S1834;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1839;
    (&_S1839)->primal_0 = _S1643;
    (&_S1839)->differential_0 = _S1758;
    s_bwd_prop_mul_3(&_S1838, &_S1839, _S1835.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1840;
    (&_S1840)->primal_0 = _S1641;
    (&_S1840)->differential_0 = _S1758;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1841;
    (&_S1841)->primal_0 = _S1642;
    (&_S1841)->differential_0 = _S1758;
    s_bwd_prop_mul_4(&_S1840, &_S1841, _S1839.differential_0);
    Matrix<float, 3, 3>  _S1842 = transpose_0(_S1841.differential_0 + _S1761.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1843;
    (&_S1843)->primal_0 = R_16;
    (&_S1843)->differential_0 = _S1758;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1844;
    (&_S1844)->primal_0 = _S1640;
    (&_S1844)->differential_0 = _S1758;
    s_bwd_prop_mul_4(&_S1843, &_S1844, _S1840.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1845;
    (&_S1845)->primal_0 = _S1638;
    (&_S1845)->differential_0 = _S1758;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1846;
    (&_S1846)->primal_0 = _S1639;
    (&_S1846)->differential_0 = _S1758;
    s_bwd_prop_mul_4(&_S1845, &_S1846, _S1844.differential_0);
    Matrix<float, 3, 3>  _S1847 = _S1845.differential_0 + transpose_0(_S1846.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1848;
    (&_S1848)->primal_0 = _S1637;
    (&_S1848)->differential_0 = _S1758;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1849;
    (&_S1849)->primal_0 = S_2;
    (&_S1849)->differential_0 = _S1758;
    s_bwd_prop_mul_4(&_S1848, &_S1849, _S1847);
    Matrix<float, 3, 3>  _S1850 = transpose_0(_S1848.differential_0);
    float _S1851 = 2.0f * - _S1850.rows[int(2)].z;
    float _S1852 = 2.0f * _S1850.rows[int(2)].y;
    float _S1853 = 2.0f * _S1850.rows[int(2)].x;
    float _S1854 = 2.0f * _S1850.rows[int(1)].z;
    float _S1855 = 2.0f * - _S1850.rows[int(1)].y;
    float _S1856 = 2.0f * _S1850.rows[int(1)].x;
    float _S1857 = 2.0f * _S1850.rows[int(0)].z;
    float _S1858 = 2.0f * _S1850.rows[int(0)].y;
    float _S1859 = 2.0f * - _S1850.rows[int(0)].x;
    float _S1860 = - _S1856 + _S1858;
    float _S1861 = _S1853 + - _S1857;
    float _S1862 = - _S1852 + _S1854;
    float _S1863 = _S1852 + _S1854;
    float _S1864 = _S1853 + _S1857;
    float _S1865 = _S1856 + _S1858;
    float _S1866 = quat_15.w * (_S1855 + _S1859);
    float _S1867 = quat_15.z * (_S1851 + _S1859);
    float _S1868 = quat_15.y * (_S1851 + _S1855);
    float _S1869 = quat_15.x * _S1860 + quat_15.z * _S1863 + quat_15.y * _S1864 + _S1866 + _S1866;
    float _S1870 = quat_15.x * _S1861 + quat_15.w * _S1863 + quat_15.y * _S1865 + _S1867 + _S1867;
    float _S1871 = quat_15.x * _S1862 + quat_15.w * _S1864 + quat_15.z * _S1865 + _S1868 + _S1868;
    float _S1872 = quat_15.w * _S1860 + quat_15.z * _S1861 + quat_15.y * _S1862;
    float3  _S1873 = _S1696;
    *&((&_S1873)->z) = _S1849.differential_0.rows[int(2)].z;
    *&((&_S1873)->y) = _S1849.differential_0.rows[int(1)].y;
    *&((&_S1873)->x) = _S1849.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1874;
    (&_S1874)->primal_0 = scale_14;
    (&_S1874)->differential_0 = _S1696;
    s_bwd_prop_exp_1(&_S1874, _S1873);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1875 = _S1874;
    float3  _S1876 = _S1696;
    *&((&_S1876)->y) = _S1830;
    *&((&_S1876)->x) = _S1831;
    float3  _S1877 = _S1798 + _S1876;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1878;
    (&_S1878)->primal_0 = R_16;
    (&_S1878)->differential_0 = _S1758;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1879;
    (&_S1879)->primal_0 = mean_12;
    (&_S1879)->differential_0 = _S1696;
    s_bwd_prop_mul_1(&_S1878, &_S1879, _S1877);
    float3  _S1880 = _S1877 + _S1762.differential_0;
    Matrix<float, 3, 3>  _S1881 = _S1842 + _S1843.differential_0 + _S1878.differential_0;
    float4  _S1882 = make_float4 (0.0f);
    *&((&_S1882)->w) = _S1869;
    *&((&_S1882)->z) = _S1870;
    *&((&_S1882)->y) = _S1871;
    *&((&_S1882)->x) = _S1872;
    float4  _S1883 = _S1882;
    float3  _S1884 = _S1879.differential_0 + _S1756;
    *v_mean_2 = _S1884;
    *v_quat_2 = _S1883;
    *v_scale_2 = _S1875.differential_0;
    *v_in_opacity_2 = _S1805;
    (*v_sh_coeffs_2)[int(0)] = _S1782;
    (*v_sh_coeffs_2)[int(1)] = _S1783;
    (*v_sh_coeffs_2)[int(2)] = _S1784;
    (*v_sh_coeffs_2)[int(3)] = _S1785;
    (*v_sh_coeffs_2)[int(4)] = _S1786;
    (*v_sh_coeffs_2)[int(5)] = _S1787;
    (*v_sh_coeffs_2)[int(6)] = _S1788;
    (*v_sh_coeffs_2)[int(7)] = _S1789;
    (*v_sh_coeffs_2)[int(8)] = _S1790;
    (*v_sh_coeffs_2)[int(9)] = _S1791;
    (*v_sh_coeffs_2)[int(10)] = _S1792;
    (*v_sh_coeffs_2)[int(11)] = _S1793;
    (*v_sh_coeffs_2)[int(12)] = _S1794;
    (*v_sh_coeffs_2)[int(13)] = _S1795;
    (*v_sh_coeffs_2)[int(14)] = _S1796;
    (*v_sh_coeffs_2)[int(15)] = _S1797;
    *v_R_3 = _S1881;
    *v_t_3 = _S1880;
    return;
}

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_17, float fy_17, float cx_17, float cy_17, float4  radial_coeffs_20, float2  tangential_coeffs_20, float2  thin_prism_coeffs_20, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_13) + t_16;
    float3  _S1885 = s_primal_ctx_exp_0(scale_15);
    float _S1886 = quat_16.y;
    float x2_16 = _S1886 * _S1886;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_29 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S1887 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    Matrix<float, 3, 3>  S_3 = makeMatrix<float, 3, 3> (_S1885.x, 0.0f, 0.0f, 0.0f, _S1885.y, 0.0f, 0.0f, 0.0f, _S1885.z);
    Matrix<float, 3, 3>  _S1888 = s_primal_ctx_mul_2(_S1887, S_3);
    Matrix<float, 3, 3>  _S1889 = transpose_0(_S1888);
    Matrix<float, 3, 3>  _S1890 = s_primal_ctx_mul_2(_S1888, _S1889);
    Matrix<float, 3, 3>  _S1891 = s_primal_ctx_mul_2(R_17, _S1890);
    Matrix<float, 3, 3>  _S1892 = transpose_0(R_17);
    Matrix<float, 3, 3>  _S1893 = s_primal_ctx_mul_2(_S1891, _S1892);
    float _S1894 = float(image_width_13);
    float _S1895 = float(image_height_13);
    float _S1896 = 0.30000001192092896f * (0.5f * _S1894 / fx_17);
    float lim_x_pos_1 = (_S1894 - cx_17) / fx_17 + _S1896;
    float _S1897 = 0.30000001192092896f * (0.5f * _S1895 / fy_17);
    float lim_y_pos_1 = (_S1895 - cy_17) / fy_17 + _S1897;
    float rz_7 = 1.0f / mean_c_13.z;
    float _S1898 = mean_c_13.z * mean_c_13.z;
    float rz2_7 = rz_7 * rz_7;
    float _S1899 = - (cx_17 / fx_17 + _S1896);
    float _S1900 = mean_c_13.x * rz_7;
    float _S1901 = s_primal_ctx_max_0(_S1899, _S1900);
    float _S1902 = s_primal_ctx_min_0(lim_x_pos_1, _S1901);
    float _S1903 = - (cy_17 / fy_17 + _S1897);
    float _S1904 = mean_c_13.y * rz_7;
    float _S1905 = s_primal_ctx_max_0(_S1903, _S1904);
    float _S1906 = s_primal_ctx_min_0(lim_y_pos_1, _S1905);
    float _S1907 = - fx_17;
    float _S1908 = _S1907 * (mean_c_13.z * _S1902);
    float _S1909 = - fy_17;
    float _S1910 = _S1909 * (mean_c_13.z * _S1906);
    Matrix<float, 2, 3>  J_15 = makeMatrix<float, 2, 3> (fx_17 * rz_7, 0.0f, _S1908 * rz2_7, 0.0f, fy_17 * rz_7, _S1910 * rz2_7);
    Matrix<float, 2, 3>  _S1911 = s_primal_ctx_mul_3(J_15, _S1893);
    Matrix<float, 3, 2>  _S1912 = transpose_1(J_15);
    Matrix<float, 2, 2>  _S1913 = s_primal_ctx_mul_4(_S1911, _S1912);
    float _S1914 = fx_17 * mean_c_13.x;
    float _S1915 = fy_17 * mean_c_13.y;
    float _S1916 = _S1913.rows[int(0)].y * _S1913.rows[int(1)].x;
    float det_orig_14 = _S1913.rows[int(0)].x * _S1913.rows[int(1)].y - _S1916;
    float _S1917 = _S1913.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1918 = _S1913;
    *&(((&_S1918)->rows + (int(0)))->x) = _S1917;
    float _S1919 = _S1913.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1918)->rows + (int(1)))->y) = _S1919;
    Matrix<float, 2, 2>  _S1920 = _S1918;
    Matrix<float, 2, 2>  _S1921 = _S1918;
    float det_blur_9 = _S1917 * _S1919 - _S1916;
    float _S1922 = det_orig_14 / det_blur_9;
    float _S1923 = det_blur_9 * det_blur_9;
    float _S1924 = s_primal_ctx_max_0(0.0f, _S1922);
    float _S1925 = s_primal_ctx_sqrt_0(_S1924);
    float _S1926 = - in_opacity_13;
    float _S1927 = 1.0f + s_primal_ctx_exp_1(_S1926);
    float _S1928 = 1.0f / _S1927;
    float _S1929 = _S1927 * _S1927;
    float _S1930;
    if(antialiased_13)
    {
        _S1930 = _S1928 * _S1925;
    }
    else
    {
        _S1930 = _S1928;
    }
    float _S1931 = _S1930 / 0.00392156885936856f;
    float _S1932 = 2.0f * s_primal_ctx_log_0(_S1931);
    float _S1933 = s_primal_ctx_sqrt_0(_S1932);
    float _S1934 = _S1920.rows[int(0)].x;
    float _S1935 = _S1921.rows[int(1)].y;
    float _S1936 = s_primal_ctx_dot_0(mean_c_13, mean_c_13) + 9.99999997475242708e-07f;
    float3  _S1937 = - scale_15;
    float3  _S1938 = mean_13 - - s_primal_ctx_mul_1(_S1892, t_16);
    float _S1939 = _S1938.x;
    float _S1940 = _S1938.y;
    float _S1941 = _S1938.z;
    float _S1942 = _S1939 * _S1939 + _S1940 * _S1940 + _S1941 * _S1941;
    float _S1943 = s_primal_ctx_sqrt_0(_S1942);
    float x_37 = _S1939 / _S1943;
    float3  _S1944 = make_float3 (x_37);
    float _S1945 = _S1943 * _S1943;
    float y_16 = _S1940 / _S1943;
    float z_13 = _S1941 / _S1943;
    float3  _S1946 = make_float3 (z_13);
    float _S1947 = - y_16;
    float3  _S1948 = make_float3 (_S1947);
    float z2_30 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_37 * x_37 - y_16 * y_16;
    float _S1949 = 2.0f * x_37;
    float fS1_13 = _S1949 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S1950 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_37;
    float3  _S1951 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S1952 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S1953 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S1954 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S1955 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S1955;
    float3  _S1956 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_37;
    float3  _S1957 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S1958 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S1959 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S1960 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_37 * fC1_13 - y_16 * fS1_13);
    float3  _S1961 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_37 * fS1_13 + y_16 * fC1_13);
    float3  _S1962 = make_float3 (pSH9_3);
    float3  _S1963 = make_float3 (0.0f);
    float3  _S1964 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1965;
    (&_S1965)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1947) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_37) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S1965)->differential_0 = _S1964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1966;
    (&_S1966)->primal_0 = _S1963;
    (&_S1966)->differential_0 = _S1964;
    s_bwd_prop_max_0(&_S1965, &_S1966, v_rgb_3);
    float3  _S1967 = _S1961 * _S1965.differential_0;
    float3  _S1968 = (*sh_coeffs_13)[int(15)] * _S1965.differential_0;
    float3  _S1969 = _S1959 * _S1965.differential_0;
    float3  _S1970 = (*sh_coeffs_13)[int(14)] * _S1965.differential_0;
    float3  _S1971 = _S1957 * _S1965.differential_0;
    float3  _S1972 = (*sh_coeffs_13)[int(13)] * _S1965.differential_0;
    float3  _S1973 = _S1956 * _S1965.differential_0;
    float3  _S1974 = (*sh_coeffs_13)[int(12)] * _S1965.differential_0;
    float3  _S1975 = _S1958 * _S1965.differential_0;
    float3  _S1976 = (*sh_coeffs_13)[int(11)] * _S1965.differential_0;
    float3  _S1977 = _S1960 * _S1965.differential_0;
    float3  _S1978 = (*sh_coeffs_13)[int(10)] * _S1965.differential_0;
    float3  _S1979 = _S1962 * _S1965.differential_0;
    float3  _S1980 = (*sh_coeffs_13)[int(9)] * _S1965.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S1980.x + _S1980.y + _S1980.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S1968.x + _S1968.y + _S1968.z);
    float _S1981 = _S1978.x + _S1978.y + _S1978.z;
    float _S1982 = _S1970.x + _S1970.y + _S1970.z;
    float _S1983 = _S1976.x + _S1976.y + _S1976.z;
    float _S1984 = _S1972.x + _S1972.y + _S1972.z;
    float _S1985 = _S1974.x + _S1974.y + _S1974.z;
    float _S1986 = - s_diff_fC2_T_3;
    float3  _S1987 = _S1953 * _S1965.differential_0;
    float3  _S1988 = (*sh_coeffs_13)[int(8)] * _S1965.differential_0;
    float3  _S1989 = _S1951 * _S1965.differential_0;
    float3  _S1990 = (*sh_coeffs_13)[int(7)] * _S1965.differential_0;
    float3  _S1991 = _S1950 * _S1965.differential_0;
    float3  _S1992 = (*sh_coeffs_13)[int(6)] * _S1965.differential_0;
    float3  _S1993 = _S1952 * _S1965.differential_0;
    float3  _S1994 = (*sh_coeffs_13)[int(5)] * _S1965.differential_0;
    float3  _S1995 = _S1954 * _S1965.differential_0;
    float3  _S1996 = (*sh_coeffs_13)[int(4)] * _S1965.differential_0;
    float _S1997 = _S1994.x + _S1994.y + _S1994.z;
    float _S1998 = _S1990.x + _S1990.y + _S1990.z;
    float _S1999 = fTmp1B_13 * _S1981 + x_37 * s_diff_fS2_T_3 + y_16 * _S1986 + 0.54627424478530884f * (_S1996.x + _S1996.y + _S1996.z);
    float _S2000 = fTmp1B_13 * _S1982 + y_16 * s_diff_fS2_T_3 + x_37 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S1988.x + _S1988.y + _S1988.z);
    float _S2001 = y_16 * - _S2000;
    float _S2002 = x_37 * _S2000;
    float _S2003 = z_13 * (1.86588168144226074f * (z_13 * _S1985) + -2.28522896766662598f * (y_16 * _S1983 + x_37 * _S1984) + 0.94617468118667603f * (_S1992.x + _S1992.y + _S1992.z));
    float3  _S2004 = make_float3 (0.48860251903533936f) * _S1965.differential_0;
    float3  _S2005 = - _S2004;
    float3  _S2006 = _S1944 * _S2005;
    float3  _S2007 = (*sh_coeffs_13)[int(3)] * _S2005;
    float3  _S2008 = _S1946 * _S2004;
    float3  _S2009 = (*sh_coeffs_13)[int(2)] * _S2004;
    float3  _S2010 = _S1948 * _S2004;
    float3  _S2011 = (*sh_coeffs_13)[int(1)] * _S2004;
    float _S2012 = (_S1955 * _S1985 + 1.44530570507049561f * (fS1_13 * _S1981 + fC1_13 * _S1982) + -1.09254848957061768f * (y_16 * _S1997 + x_37 * _S1998) + _S2003 + _S2003 + _S2009.x + _S2009.y + _S2009.z) / _S1945;
    float _S2013 = _S1943 * _S2012;
    float _S2014 = (fTmp0C_13 * _S1983 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S1986 + fTmp0B_13 * _S1997 + _S1949 * _S1999 + _S2001 + _S2001 + - (_S2011.x + _S2011.y + _S2011.z)) / _S1945;
    float _S2015 = _S1943 * _S2014;
    float _S2016 = (fTmp0C_13 * _S1984 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S1998 + 2.0f * (y_16 * _S1999) + _S2002 + _S2002 + _S2007.x + _S2007.y + _S2007.z) / _S1945;
    float _S2017 = _S1943 * _S2016;
    float _S2018 = _S1941 * - _S2012 + _S1940 * - _S2014 + _S1939 * - _S2016;
    DiffPair_float_0 _S2019;
    (&_S2019)->primal_0 = _S1942;
    (&_S2019)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2019, _S2018);
    float _S2020 = _S1941 * _S2019.differential_0;
    float _S2021 = _S1940 * _S2019.differential_0;
    float _S2022 = _S1939 * _S2019.differential_0;
    float3  _S2023 = make_float3 (0.282094806432724f) * _S1965.differential_0;
    float3  _S2024 = make_float3 (_S2017 + _S2022 + _S2022, _S2015 + _S2021 + _S2021, _S2013 + _S2020 + _S2020);
    float3  _S2025 = - - _S2024;
    Matrix<float, 3, 3>  _S2026 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2027;
    (&_S2027)->primal_0 = _S1892;
    (&_S2027)->differential_0 = _S2026;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2028;
    (&_S2028)->primal_0 = t_16;
    (&_S2028)->differential_0 = _S1964;
    s_bwd_prop_mul_1(&_S2027, &_S2028, _S2025);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2029 = _S2027;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2030 = _S2028;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2031;
    (&_S2031)->primal_0 = _S1937;
    (&_S2031)->differential_0 = _S1964;
    s_bwd_prop_exp_1(&_S2031, v_conic_3);
    float3  _S2032 = - _S2031.differential_0;
    float _S2033 = 0.5f * v_depth_3;
    DiffPair_float_0 _S2034;
    (&_S2034)->primal_0 = _S1936;
    (&_S2034)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2034, _S2033);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2035;
    (&_S2035)->primal_0 = mean_c_13;
    (&_S2035)->differential_0 = _S1964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2036;
    (&_S2036)->primal_0 = mean_c_13;
    (&_S2036)->differential_0 = _S1964;
    s_bwd_prop_dot_0(&_S2035, &_S2036, _S2034.differential_0);
    DiffPair_float_0 _S2037;
    (&_S2037)->primal_0 = _S1935;
    (&_S2037)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2037, 0.0f);
    DiffPair_float_0 _S2038;
    (&_S2038)->primal_0 = _S1934;
    (&_S2038)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2038, 0.0f);
    DiffPair_float_0 _S2039;
    (&_S2039)->primal_0 = 3.32999992370605469f;
    (&_S2039)->differential_0 = 0.0f;
    DiffPair_float_0 _S2040;
    (&_S2040)->primal_0 = _S1933;
    (&_S2040)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2039, &_S2040, 0.0f);
    DiffPair_float_0 _S2041;
    (&_S2041)->primal_0 = _S1932;
    (&_S2041)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2041, _S2040.differential_0);
    float _S2042 = 2.0f * _S2041.differential_0;
    DiffPair_float_0 _S2043;
    (&_S2043)->primal_0 = _S1931;
    (&_S2043)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2043, _S2042);
    float _S2044 = v_opacity_3 + 254.9999847412109375f * _S2043.differential_0;
    FixedArray<float3 , 16>  _S2045;
    _S2045[int(0)] = _S1964;
    _S2045[int(1)] = _S1964;
    _S2045[int(2)] = _S1964;
    _S2045[int(3)] = _S1964;
    _S2045[int(4)] = _S1964;
    _S2045[int(5)] = _S1964;
    _S2045[int(6)] = _S1964;
    _S2045[int(7)] = _S1964;
    _S2045[int(8)] = _S1964;
    _S2045[int(9)] = _S1964;
    _S2045[int(10)] = _S1964;
    _S2045[int(11)] = _S1964;
    _S2045[int(12)] = _S1964;
    _S2045[int(13)] = _S1964;
    _S2045[int(14)] = _S1964;
    _S2045[int(15)] = _S1964;
    _S2045[int(7)] = _S1989;
    _S2045[int(0)] = _S2023;
    _S2045[int(1)] = _S2010;
    _S2045[int(2)] = _S2008;
    _S2045[int(3)] = _S2006;
    _S2045[int(4)] = _S1995;
    _S2045[int(5)] = _S1993;
    _S2045[int(6)] = _S1991;
    _S2045[int(15)] = _S1967;
    _S2045[int(8)] = _S1987;
    _S2045[int(9)] = _S1979;
    _S2045[int(10)] = _S1977;
    _S2045[int(11)] = _S1975;
    _S2045[int(12)] = _S1973;
    _S2045[int(13)] = _S1971;
    _S2045[int(14)] = _S1969;
    float3  _S2046 = _S2045[int(0)];
    float3  _S2047 = _S2045[int(1)];
    float3  _S2048 = _S2045[int(2)];
    float3  _S2049 = _S2045[int(3)];
    float3  _S2050 = _S2045[int(4)];
    float3  _S2051 = _S2045[int(5)];
    float3  _S2052 = _S2045[int(6)];
    float3  _S2053 = _S2045[int(7)];
    float3  _S2054 = _S2045[int(8)];
    float3  _S2055 = _S2045[int(9)];
    float3  _S2056 = _S2045[int(10)];
    float3  _S2057 = _S2045[int(11)];
    float3  _S2058 = _S2045[int(12)];
    float3  _S2059 = _S2045[int(13)];
    float3  _S2060 = _S2045[int(14)];
    float3  _S2061 = _S2045[int(15)];
    float3  _S2062 = _S2036.differential_0 + _S2035.differential_0;
    float2  _S2063 = make_float2 (0.0f, _S2037.differential_0);
    float2  _S2064 = make_float2 (_S2038.differential_0, 0.0f);
    float _S2065;
    if(antialiased_13)
    {
        float _S2066 = _S1928 * _S2044;
        _S1930 = _S1925 * _S2044;
        _S2065 = _S2066;
    }
    else
    {
        _S1930 = _S2044;
        _S2065 = 0.0f;
    }
    float _S2067 = - (_S1930 / _S1929);
    DiffPair_float_0 _S2068;
    (&_S2068)->primal_0 = _S1926;
    (&_S2068)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2068, _S2067);
    float _S2069 = - _S2068.differential_0;
    DiffPair_float_0 _S2070;
    (&_S2070)->primal_0 = _S1924;
    (&_S2070)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2070, _S2065);
    DiffPair_float_0 _S2071;
    (&_S2071)->primal_0 = 0.0f;
    (&_S2071)->differential_0 = 0.0f;
    DiffPair_float_0 _S2072;
    (&_S2072)->primal_0 = _S1922;
    (&_S2072)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2071, &_S2072, _S2070.differential_0);
    float _S2073 = _S2072.differential_0 / _S1923;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2073;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2073;
    float _S2074 = - s_diff_det_blur_T_0;
    float _S2075 = _S1917 * s_diff_det_blur_T_0;
    float _S2076 = _S1919 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2077 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2078 = _S2077;
    _S2078[int(1)] = _S2063;
    _S2078[int(0)] = _S2064;
    _S1918 = _S2078;
    *&(((&_S1918)->rows + (int(1)))->y) = 0.0f;
    float _S2079 = _S2075 + _S2078.rows[int(1)].y;
    *&(((&_S1918)->rows + (int(0)))->x) = 0.0f;
    float _S2080 = _S2076 + _S2078.rows[int(0)].x;
    float _S2081 = _S2074 + - s_diff_det_orig_T_3;
    float _S2082 = _S1913.rows[int(0)].y * _S2081;
    float _S2083 = _S1913.rows[int(1)].x * _S2081;
    float _S2084 = _S1913.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2085 = _S2079 + _S1913.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2086 = make_float2 (0.0f);
    float2  _S2087 = _S2086;
    *&((&_S2087)->x) = _S2082;
    *&((&_S2087)->y) = _S2085;
    float _S2088 = _S2080 + _S2084;
    float2  _S2089 = _S2086;
    *&((&_S2089)->y) = _S2083;
    *&((&_S2089)->x) = _S2088;
    float _S2090 = _S1915 * v_mean2d_3.y;
    float _S2091 = fy_17 * (rz_7 * v_mean2d_3.y);
    float _S2092 = _S1914 * v_mean2d_3.x;
    float _S2093 = fx_17 * (rz_7 * v_mean2d_3.x);
    Matrix<float, 2, 2>  _S2094 = _S2077;
    _S2094[int(1)] = _S2087;
    _S2094[int(0)] = _S2089;
    Matrix<float, 2, 2>  _S2095 = _S1918 + _S2094;
    Matrix<float, 2, 3>  _S2096 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2097;
    (&_S2097)->primal_0 = _S1911;
    (&_S2097)->differential_0 = _S2096;
    Matrix<float, 3, 2>  _S2098 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2099;
    (&_S2099)->primal_0 = _S1912;
    (&_S2099)->differential_0 = _S2098;
    s_bwd_prop_mul_2(&_S2097, &_S2099, _S2095);
    Matrix<float, 2, 3>  _S2100 = transpose_2(_S2099.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2101;
    (&_S2101)->primal_0 = J_15;
    (&_S2101)->differential_0 = _S2096;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2102;
    (&_S2102)->primal_0 = _S1893;
    (&_S2102)->differential_0 = _S2026;
    s_bwd_prop_mul_3(&_S2101, &_S2102, _S2097.differential_0);
    Matrix<float, 2, 3>  _S2103 = _S2100 + _S2101.differential_0;
    float _S2104 = _S1910 * _S2103.rows[int(1)].z;
    float s_diff_ty_T_1 = _S1909 * (rz2_7 * _S2103.rows[int(1)].z);
    float _S2105 = fy_17 * _S2103.rows[int(1)].y;
    float _S2106 = _S1908 * _S2103.rows[int(0)].z;
    float s_diff_tx_T_1 = _S1907 * (rz2_7 * _S2103.rows[int(0)].z);
    float _S2107 = fx_17 * _S2103.rows[int(0)].x;
    float _S2108 = mean_c_13.z * s_diff_ty_T_1;
    float _S2109 = _S1906 * s_diff_ty_T_1;
    DiffPair_float_0 _S2110;
    (&_S2110)->primal_0 = lim_y_pos_1;
    (&_S2110)->differential_0 = 0.0f;
    DiffPair_float_0 _S2111;
    (&_S2111)->primal_0 = _S1905;
    (&_S2111)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2110, &_S2111, _S2108);
    DiffPair_float_0 _S2112;
    (&_S2112)->primal_0 = _S1903;
    (&_S2112)->differential_0 = 0.0f;
    DiffPair_float_0 _S2113;
    (&_S2113)->primal_0 = _S1904;
    (&_S2113)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2112, &_S2113, _S2111.differential_0);
    float _S2114 = mean_c_13.y * _S2113.differential_0;
    float _S2115 = rz_7 * _S2113.differential_0;
    float _S2116 = mean_c_13.z * s_diff_tx_T_1;
    float _S2117 = _S1902 * s_diff_tx_T_1;
    DiffPair_float_0 _S2118;
    (&_S2118)->primal_0 = lim_x_pos_1;
    (&_S2118)->differential_0 = 0.0f;
    DiffPair_float_0 _S2119;
    (&_S2119)->primal_0 = _S1901;
    (&_S2119)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2118, &_S2119, _S2116);
    DiffPair_float_0 _S2120;
    (&_S2120)->primal_0 = _S1899;
    (&_S2120)->differential_0 = 0.0f;
    DiffPair_float_0 _S2121;
    (&_S2121)->primal_0 = _S1900;
    (&_S2121)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2120, &_S2121, _S2119.differential_0);
    float _S2122 = rz_7 * (_S2104 + _S2106);
    float _S2123 = _S2109 + _S2117 + - ((_S2090 + _S2092 + _S2105 + _S2107 + _S2114 + mean_c_13.x * _S2121.differential_0 + _S2122 + _S2122) / _S1898);
    float _S2124 = _S2091 + _S2115;
    float _S2125 = _S2093 + rz_7 * _S2121.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2126;
    (&_S2126)->primal_0 = _S1891;
    (&_S2126)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2127;
    (&_S2127)->primal_0 = _S1892;
    (&_S2127)->differential_0 = _S2026;
    s_bwd_prop_mul_4(&_S2126, &_S2127, _S2102.differential_0);
    Matrix<float, 3, 3>  _S2128 = transpose_0(_S2127.differential_0 + _S2029.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2129;
    (&_S2129)->primal_0 = R_17;
    (&_S2129)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2130;
    (&_S2130)->primal_0 = _S1890;
    (&_S2130)->differential_0 = _S2026;
    s_bwd_prop_mul_4(&_S2129, &_S2130, _S2126.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2131;
    (&_S2131)->primal_0 = _S1888;
    (&_S2131)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2132;
    (&_S2132)->primal_0 = _S1889;
    (&_S2132)->differential_0 = _S2026;
    s_bwd_prop_mul_4(&_S2131, &_S2132, _S2130.differential_0);
    Matrix<float, 3, 3>  _S2133 = _S2131.differential_0 + transpose_0(_S2132.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2134;
    (&_S2134)->primal_0 = _S1887;
    (&_S2134)->differential_0 = _S2026;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2135;
    (&_S2135)->primal_0 = S_3;
    (&_S2135)->differential_0 = _S2026;
    s_bwd_prop_mul_4(&_S2134, &_S2135, _S2133);
    Matrix<float, 3, 3>  _S2136 = transpose_0(_S2134.differential_0);
    float _S2137 = 2.0f * - _S2136.rows[int(2)].z;
    float _S2138 = 2.0f * _S2136.rows[int(2)].y;
    float _S2139 = 2.0f * _S2136.rows[int(2)].x;
    float _S2140 = 2.0f * _S2136.rows[int(1)].z;
    float _S2141 = 2.0f * - _S2136.rows[int(1)].y;
    float _S2142 = 2.0f * _S2136.rows[int(1)].x;
    float _S2143 = 2.0f * _S2136.rows[int(0)].z;
    float _S2144 = 2.0f * _S2136.rows[int(0)].y;
    float _S2145 = 2.0f * - _S2136.rows[int(0)].x;
    float _S2146 = - _S2142 + _S2144;
    float _S2147 = _S2139 + - _S2143;
    float _S2148 = - _S2138 + _S2140;
    float _S2149 = _S2138 + _S2140;
    float _S2150 = _S2139 + _S2143;
    float _S2151 = _S2142 + _S2144;
    float _S2152 = quat_16.w * (_S2141 + _S2145);
    float _S2153 = quat_16.z * (_S2137 + _S2145);
    float _S2154 = quat_16.y * (_S2137 + _S2141);
    float _S2155 = quat_16.x * _S2146 + quat_16.z * _S2149 + quat_16.y * _S2150 + _S2152 + _S2152;
    float _S2156 = quat_16.x * _S2147 + quat_16.w * _S2149 + quat_16.y * _S2151 + _S2153 + _S2153;
    float _S2157 = quat_16.x * _S2148 + quat_16.w * _S2150 + quat_16.z * _S2151 + _S2154 + _S2154;
    float _S2158 = quat_16.w * _S2146 + quat_16.z * _S2147 + quat_16.y * _S2148;
    float3  _S2159 = _S1964;
    *&((&_S2159)->z) = _S2135.differential_0.rows[int(2)].z;
    *&((&_S2159)->y) = _S2135.differential_0.rows[int(1)].y;
    *&((&_S2159)->x) = _S2135.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2160;
    (&_S2160)->primal_0 = scale_15;
    (&_S2160)->differential_0 = _S1964;
    s_bwd_prop_exp_1(&_S2160, _S2159);
    float3  _S2161 = _S1964;
    *&((&_S2161)->z) = _S2123;
    *&((&_S2161)->y) = _S2124;
    *&((&_S2161)->x) = _S2125;
    float3  _S2162 = _S2062 + _S2161;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2163;
    (&_S2163)->primal_0 = R_17;
    (&_S2163)->differential_0 = _S2026;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2164;
    (&_S2164)->primal_0 = mean_13;
    (&_S2164)->differential_0 = _S1964;
    s_bwd_prop_mul_1(&_S2163, &_S2164, _S2162);
    float3  _S2165 = _S2162 + _S2030.differential_0;
    Matrix<float, 3, 3>  _S2166 = _S2128 + _S2129.differential_0 + _S2163.differential_0;
    float3  _S2167 = _S2160.differential_0 + _S2032;
    float4  _S2168 = make_float4 (0.0f);
    *&((&_S2168)->w) = _S2155;
    *&((&_S2168)->z) = _S2156;
    *&((&_S2168)->y) = _S2157;
    *&((&_S2168)->x) = _S2158;
    float4  _S2169 = _S2168;
    float3  _S2170 = _S2164.differential_0 + _S2024;
    *v_mean_3 = _S2170;
    *v_quat_3 = _S2169;
    *v_scale_3 = _S2167;
    *v_in_opacity_3 = _S2069;
    (*v_sh_coeffs_3)[int(0)] = _S2046;
    (*v_sh_coeffs_3)[int(1)] = _S2047;
    (*v_sh_coeffs_3)[int(2)] = _S2048;
    (*v_sh_coeffs_3)[int(3)] = _S2049;
    (*v_sh_coeffs_3)[int(4)] = _S2050;
    (*v_sh_coeffs_3)[int(5)] = _S2051;
    (*v_sh_coeffs_3)[int(6)] = _S2052;
    (*v_sh_coeffs_3)[int(7)] = _S2053;
    (*v_sh_coeffs_3)[int(8)] = _S2054;
    (*v_sh_coeffs_3)[int(9)] = _S2055;
    (*v_sh_coeffs_3)[int(10)] = _S2056;
    (*v_sh_coeffs_3)[int(11)] = _S2057;
    (*v_sh_coeffs_3)[int(12)] = _S2058;
    (*v_sh_coeffs_3)[int(13)] = _S2059;
    (*v_sh_coeffs_3)[int(14)] = _S2060;
    (*v_sh_coeffs_3)[int(15)] = _S2061;
    *v_R_4 = _S2166;
    *v_t_4 = _S2165;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2171;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2172;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_18, float fy_18, float cx_18, float cy_18, float4  radial_coeffs_21, float2  tangential_coeffs_21, float2  thin_prism_coeffs_21, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2173 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S2174 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S2175 = { _S2174, _S2174 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2176 = { _S2174, _S2174, _S2175, _S2174, _S2174, _S2175 };
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2177;
    (&_S2177)->_S2171 = _S2173;
    (&_S2177)->_S2172 = _S2176;
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_14) + t_17;
    float3  _S2178 = s_primal_ctx_exp_0(scale_16);
    float _S2179 = quat_17.y;
    float x2_17 = _S2179 * _S2179;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_31 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2180 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    Matrix<float, 3, 3>  S_4 = makeMatrix<float, 3, 3> (_S2178.x, 0.0f, 0.0f, 0.0f, _S2178.y, 0.0f, 0.0f, 0.0f, _S2178.z);
    Matrix<float, 3, 3>  _S2181 = s_primal_ctx_mul_2(_S2180, S_4);
    Matrix<float, 3, 3>  _S2182 = transpose_0(_S2181);
    Matrix<float, 3, 3>  _S2183 = s_primal_ctx_mul_2(_S2181, _S2182);
    Matrix<float, 3, 3>  _S2184 = s_primal_ctx_mul_2(R_18, _S2183);
    Matrix<float, 3, 3>  _S2185 = transpose_0(R_18);
    Matrix<float, 3, 3>  _S2186 = s_primal_ctx_mul_2(_S2184, _S2185);
    Matrix<float, 2, 2>  _S2187 = _S2173;
    float2  _S2188 = make_float2 (0.0f);
    float2  _S2189 = _S2188;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_14, _S2186, fx_18, fy_18, cx_18, cy_18, radial_coeffs_21, tangential_coeffs_21, thin_prism_coeffs_21, &_S2187, &_S2189, &(&_S2177)->_S2172);
    (&_S2177)->_S2171 = _S2187;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2190 = _S2177;
    float _S2191 = _S2177._S2171.rows[int(0)].y * _S2177._S2171.rows[int(1)].x;
    float det_orig_15 = _S2177._S2171.rows[int(0)].x * _S2177._S2171.rows[int(1)].y - _S2191;
    float _S2192 = _S2177._S2171.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2193 = _S2177._S2171;
    *&(((&_S2193)->rows + (int(0)))->x) = _S2192;
    float _S2194 = _S2177._S2171.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2193)->rows + (int(1)))->y) = _S2194;
    Matrix<float, 2, 2>  _S2195 = _S2193;
    Matrix<float, 2, 2>  _S2196 = _S2193;
    float det_blur_10 = _S2192 * _S2194 - _S2191;
    float _S2197 = det_orig_15 / det_blur_10;
    float _S2198 = det_blur_10 * det_blur_10;
    float _S2199 = s_primal_ctx_max_0(0.0f, _S2197);
    float _S2200 = s_primal_ctx_sqrt_0(_S2199);
    float _S2201 = - in_opacity_14;
    float _S2202 = 1.0f + s_primal_ctx_exp_1(_S2201);
    float _S2203 = 1.0f / _S2202;
    float _S2204 = _S2202 * _S2202;
    float _S2205;
    if(antialiased_14)
    {
        _S2205 = _S2203 * _S2200;
    }
    else
    {
        _S2205 = _S2203;
    }
    float _S2206 = _S2205 / 0.00392156885936856f;
    float _S2207 = 2.0f * s_primal_ctx_log_0(_S2206);
    float _S2208 = s_primal_ctx_sqrt_0(_S2207);
    float _S2209 = _S2195.rows[int(0)].x;
    float _S2210 = _S2196.rows[int(1)].y;
    float _S2211 = s_primal_ctx_dot_0(mean_c_14, mean_c_14) + 9.99999997475242708e-07f;
    float3  _S2212 = - scale_16;
    float3  _S2213 = mean_14 - - s_primal_ctx_mul_1(_S2185, t_17);
    float _S2214 = _S2213.x;
    float _S2215 = _S2213.y;
    float _S2216 = _S2213.z;
    float _S2217 = _S2214 * _S2214 + _S2215 * _S2215 + _S2216 * _S2216;
    float _S2218 = s_primal_ctx_sqrt_0(_S2217);
    float x_38 = _S2214 / _S2218;
    float3  _S2219 = make_float3 (x_38);
    float _S2220 = _S2218 * _S2218;
    float y_17 = _S2215 / _S2218;
    float z_14 = _S2216 / _S2218;
    float3  _S2221 = make_float3 (z_14);
    float _S2222 = - y_17;
    float3  _S2223 = make_float3 (_S2222);
    float z2_32 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_38 * x_38 - y_17 * y_17;
    float _S2224 = 2.0f * x_38;
    float fS1_14 = _S2224 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2225 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_38;
    float3  _S2226 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2227 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2228 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2229 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2230 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2230;
    float3  _S2231 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_38;
    float3  _S2232 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2233 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2234 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2235 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_38 * fC1_14 - y_17 * fS1_14);
    float3  _S2236 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_38 * fS1_14 + y_17 * fC1_14);
    float3  _S2237 = make_float3 (pSH9_4);
    float3  _S2238 = make_float3 (0.0f);
    float3  _S2239 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2240;
    (&_S2240)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2222) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_38) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2240)->differential_0 = _S2239;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2241;
    (&_S2241)->primal_0 = _S2238;
    (&_S2241)->differential_0 = _S2239;
    s_bwd_prop_max_0(&_S2240, &_S2241, v_rgb_4);
    float3  _S2242 = _S2236 * _S2240.differential_0;
    float3  _S2243 = (*sh_coeffs_14)[int(15)] * _S2240.differential_0;
    float3  _S2244 = _S2234 * _S2240.differential_0;
    float3  _S2245 = (*sh_coeffs_14)[int(14)] * _S2240.differential_0;
    float3  _S2246 = _S2232 * _S2240.differential_0;
    float3  _S2247 = (*sh_coeffs_14)[int(13)] * _S2240.differential_0;
    float3  _S2248 = _S2231 * _S2240.differential_0;
    float3  _S2249 = (*sh_coeffs_14)[int(12)] * _S2240.differential_0;
    float3  _S2250 = _S2233 * _S2240.differential_0;
    float3  _S2251 = (*sh_coeffs_14)[int(11)] * _S2240.differential_0;
    float3  _S2252 = _S2235 * _S2240.differential_0;
    float3  _S2253 = (*sh_coeffs_14)[int(10)] * _S2240.differential_0;
    float3  _S2254 = _S2237 * _S2240.differential_0;
    float3  _S2255 = (*sh_coeffs_14)[int(9)] * _S2240.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2255.x + _S2255.y + _S2255.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2243.x + _S2243.y + _S2243.z);
    float _S2256 = _S2253.x + _S2253.y + _S2253.z;
    float _S2257 = _S2245.x + _S2245.y + _S2245.z;
    float _S2258 = _S2251.x + _S2251.y + _S2251.z;
    float _S2259 = _S2247.x + _S2247.y + _S2247.z;
    float _S2260 = _S2249.x + _S2249.y + _S2249.z;
    float _S2261 = - s_diff_fC2_T_4;
    float3  _S2262 = _S2228 * _S2240.differential_0;
    float3  _S2263 = (*sh_coeffs_14)[int(8)] * _S2240.differential_0;
    float3  _S2264 = _S2226 * _S2240.differential_0;
    float3  _S2265 = (*sh_coeffs_14)[int(7)] * _S2240.differential_0;
    float3  _S2266 = _S2225 * _S2240.differential_0;
    float3  _S2267 = (*sh_coeffs_14)[int(6)] * _S2240.differential_0;
    float3  _S2268 = _S2227 * _S2240.differential_0;
    float3  _S2269 = (*sh_coeffs_14)[int(5)] * _S2240.differential_0;
    float3  _S2270 = _S2229 * _S2240.differential_0;
    float3  _S2271 = (*sh_coeffs_14)[int(4)] * _S2240.differential_0;
    float _S2272 = _S2269.x + _S2269.y + _S2269.z;
    float _S2273 = _S2265.x + _S2265.y + _S2265.z;
    float _S2274 = fTmp1B_14 * _S2256 + x_38 * s_diff_fS2_T_4 + y_17 * _S2261 + 0.54627424478530884f * (_S2271.x + _S2271.y + _S2271.z);
    float _S2275 = fTmp1B_14 * _S2257 + y_17 * s_diff_fS2_T_4 + x_38 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2263.x + _S2263.y + _S2263.z);
    float _S2276 = y_17 * - _S2275;
    float _S2277 = x_38 * _S2275;
    float _S2278 = z_14 * (1.86588168144226074f * (z_14 * _S2260) + -2.28522896766662598f * (y_17 * _S2258 + x_38 * _S2259) + 0.94617468118667603f * (_S2267.x + _S2267.y + _S2267.z));
    float3  _S2279 = make_float3 (0.48860251903533936f) * _S2240.differential_0;
    float3  _S2280 = - _S2279;
    float3  _S2281 = _S2219 * _S2280;
    float3  _S2282 = (*sh_coeffs_14)[int(3)] * _S2280;
    float3  _S2283 = _S2221 * _S2279;
    float3  _S2284 = (*sh_coeffs_14)[int(2)] * _S2279;
    float3  _S2285 = _S2223 * _S2279;
    float3  _S2286 = (*sh_coeffs_14)[int(1)] * _S2279;
    float _S2287 = (_S2230 * _S2260 + 1.44530570507049561f * (fS1_14 * _S2256 + fC1_14 * _S2257) + -1.09254848957061768f * (y_17 * _S2272 + x_38 * _S2273) + _S2278 + _S2278 + _S2284.x + _S2284.y + _S2284.z) / _S2220;
    float _S2288 = _S2218 * _S2287;
    float _S2289 = (fTmp0C_14 * _S2258 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2261 + fTmp0B_14 * _S2272 + _S2224 * _S2274 + _S2276 + _S2276 + - (_S2286.x + _S2286.y + _S2286.z)) / _S2220;
    float _S2290 = _S2218 * _S2289;
    float _S2291 = (fTmp0C_14 * _S2259 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2273 + 2.0f * (y_17 * _S2274) + _S2277 + _S2277 + _S2282.x + _S2282.y + _S2282.z) / _S2220;
    float _S2292 = _S2218 * _S2291;
    float _S2293 = _S2216 * - _S2287 + _S2215 * - _S2289 + _S2214 * - _S2291;
    DiffPair_float_0 _S2294;
    (&_S2294)->primal_0 = _S2217;
    (&_S2294)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2294, _S2293);
    float _S2295 = _S2216 * _S2294.differential_0;
    float _S2296 = _S2215 * _S2294.differential_0;
    float _S2297 = _S2214 * _S2294.differential_0;
    float3  _S2298 = make_float3 (0.282094806432724f) * _S2240.differential_0;
    float3  _S2299 = make_float3 (_S2292 + _S2297 + _S2297, _S2290 + _S2296 + _S2296, _S2288 + _S2295 + _S2295);
    float3  _S2300 = - - _S2299;
    Matrix<float, 3, 3>  _S2301 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2302;
    (&_S2302)->primal_0 = _S2185;
    (&_S2302)->differential_0 = _S2301;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2303;
    (&_S2303)->primal_0 = t_17;
    (&_S2303)->differential_0 = _S2239;
    s_bwd_prop_mul_1(&_S2302, &_S2303, _S2300);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2304 = _S2302;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2305 = _S2303;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2306;
    (&_S2306)->primal_0 = _S2212;
    (&_S2306)->differential_0 = _S2239;
    s_bwd_prop_exp_1(&_S2306, v_conic_4);
    float3  _S2307 = - _S2306.differential_0;
    float _S2308 = 0.5f * v_depth_4;
    DiffPair_float_0 _S2309;
    (&_S2309)->primal_0 = _S2211;
    (&_S2309)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2309, _S2308);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2310;
    (&_S2310)->primal_0 = mean_c_14;
    (&_S2310)->differential_0 = _S2239;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2311;
    (&_S2311)->primal_0 = mean_c_14;
    (&_S2311)->differential_0 = _S2239;
    s_bwd_prop_dot_0(&_S2310, &_S2311, _S2309.differential_0);
    DiffPair_float_0 _S2312;
    (&_S2312)->primal_0 = _S2210;
    (&_S2312)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2312, 0.0f);
    DiffPair_float_0 _S2313;
    (&_S2313)->primal_0 = _S2209;
    (&_S2313)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2313, 0.0f);
    DiffPair_float_0 _S2314;
    (&_S2314)->primal_0 = 3.32999992370605469f;
    (&_S2314)->differential_0 = 0.0f;
    DiffPair_float_0 _S2315;
    (&_S2315)->primal_0 = _S2208;
    (&_S2315)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2314, &_S2315, 0.0f);
    DiffPair_float_0 _S2316;
    (&_S2316)->primal_0 = _S2207;
    (&_S2316)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2316, _S2315.differential_0);
    float _S2317 = 2.0f * _S2316.differential_0;
    DiffPair_float_0 _S2318;
    (&_S2318)->primal_0 = _S2206;
    (&_S2318)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2318, _S2317);
    float2  _S2319 = make_float2 (_S2313.differential_0, 0.0f);
    float _S2320 = v_opacity_4 + 254.9999847412109375f * _S2318.differential_0;
    FixedArray<float3 , 16>  _S2321;
    _S2321[int(0)] = _S2239;
    _S2321[int(1)] = _S2239;
    _S2321[int(2)] = _S2239;
    _S2321[int(3)] = _S2239;
    _S2321[int(4)] = _S2239;
    _S2321[int(5)] = _S2239;
    _S2321[int(6)] = _S2239;
    _S2321[int(7)] = _S2239;
    _S2321[int(8)] = _S2239;
    _S2321[int(9)] = _S2239;
    _S2321[int(10)] = _S2239;
    _S2321[int(11)] = _S2239;
    _S2321[int(12)] = _S2239;
    _S2321[int(13)] = _S2239;
    _S2321[int(14)] = _S2239;
    _S2321[int(15)] = _S2239;
    _S2321[int(7)] = _S2264;
    _S2321[int(0)] = _S2298;
    _S2321[int(1)] = _S2285;
    _S2321[int(2)] = _S2283;
    _S2321[int(3)] = _S2281;
    _S2321[int(4)] = _S2270;
    _S2321[int(5)] = _S2268;
    _S2321[int(6)] = _S2266;
    _S2321[int(15)] = _S2242;
    _S2321[int(8)] = _S2262;
    _S2321[int(9)] = _S2254;
    _S2321[int(10)] = _S2252;
    _S2321[int(11)] = _S2250;
    _S2321[int(12)] = _S2248;
    _S2321[int(13)] = _S2246;
    _S2321[int(14)] = _S2244;
    float3  _S2322 = _S2321[int(0)];
    float3  _S2323 = _S2321[int(1)];
    float3  _S2324 = _S2321[int(2)];
    float3  _S2325 = _S2321[int(3)];
    float3  _S2326 = _S2321[int(4)];
    float3  _S2327 = _S2321[int(5)];
    float3  _S2328 = _S2321[int(6)];
    float3  _S2329 = _S2321[int(7)];
    float3  _S2330 = _S2321[int(8)];
    float3  _S2331 = _S2321[int(9)];
    float3  _S2332 = _S2321[int(10)];
    float3  _S2333 = _S2321[int(11)];
    float3  _S2334 = _S2321[int(12)];
    float3  _S2335 = _S2321[int(13)];
    float3  _S2336 = _S2321[int(14)];
    float3  _S2337 = _S2321[int(15)];
    float3  _S2338 = _S2311.differential_0 + _S2310.differential_0;
    float2  _S2339 = make_float2 (0.0f, _S2312.differential_0);
    float _S2340;
    if(antialiased_14)
    {
        float _S2341 = _S2203 * _S2320;
        _S2205 = _S2200 * _S2320;
        _S2340 = _S2341;
    }
    else
    {
        _S2205 = _S2320;
        _S2340 = 0.0f;
    }
    float _S2342 = - (_S2205 / _S2204);
    DiffPair_float_0 _S2343;
    (&_S2343)->primal_0 = _S2201;
    (&_S2343)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2343, _S2342);
    float _S2344 = - _S2343.differential_0;
    DiffPair_float_0 _S2345;
    (&_S2345)->primal_0 = _S2199;
    (&_S2345)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2345, _S2340);
    DiffPair_float_0 _S2346;
    (&_S2346)->primal_0 = 0.0f;
    (&_S2346)->differential_0 = 0.0f;
    DiffPair_float_0 _S2347;
    (&_S2347)->primal_0 = _S2197;
    (&_S2347)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2346, &_S2347, _S2345.differential_0);
    float _S2348 = _S2347.differential_0 / _S2198;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2348;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2348;
    float _S2349 = - s_diff_det_blur_T_1;
    float _S2350 = _S2192 * s_diff_det_blur_T_1;
    float _S2351 = _S2194 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2352 = _S2173;
    _S2352[int(1)] = _S2339;
    _S2352[int(0)] = _S2319;
    _S2193 = _S2352;
    *&(((&_S2193)->rows + (int(1)))->y) = 0.0f;
    float _S2353 = _S2350 + _S2352.rows[int(1)].y;
    *&(((&_S2193)->rows + (int(0)))->x) = 0.0f;
    float _S2354 = _S2351 + _S2352.rows[int(0)].x;
    float _S2355 = _S2349 + - s_diff_det_orig_T_4;
    float _S2356 = _S2190._S2171.rows[int(0)].y * _S2355;
    float _S2357 = _S2190._S2171.rows[int(1)].x * _S2355;
    float _S2358 = _S2190._S2171.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2359 = _S2353 + _S2190._S2171.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2360 = _S2188;
    *&((&_S2360)->x) = _S2356;
    *&((&_S2360)->y) = _S2359;
    float _S2361 = _S2354 + _S2358;
    float2  _S2362 = _S2188;
    *&((&_S2362)->y) = _S2357;
    *&((&_S2362)->x) = _S2361;
    Matrix<float, 2, 2>  _S2363 = _S2173;
    _S2363[int(1)] = _S2360;
    _S2363[int(0)] = _S2362;
    Matrix<float, 2, 2>  _S2364 = _S2193 + _S2363;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2365;
    (&_S2365)->primal_0 = mean_c_14;
    (&_S2365)->differential_0 = _S2239;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2366;
    (&_S2366)->primal_0 = _S2186;
    (&_S2366)->differential_0 = _S2301;
    DiffPair_float_0 _S2367;
    (&_S2367)->primal_0 = fx_18;
    (&_S2367)->differential_0 = 0.0f;
    DiffPair_float_0 _S2368;
    (&_S2368)->primal_0 = fy_18;
    (&_S2368)->differential_0 = 0.0f;
    DiffPair_float_0 _S2369;
    (&_S2369)->primal_0 = cx_18;
    (&_S2369)->differential_0 = 0.0f;
    DiffPair_float_0 _S2370;
    (&_S2370)->primal_0 = cy_18;
    (&_S2370)->differential_0 = 0.0f;
    float4  _S2371 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2372;
    (&_S2372)->primal_0 = radial_coeffs_21;
    (&_S2372)->differential_0 = _S2371;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2373;
    (&_S2373)->primal_0 = tangential_coeffs_21;
    (&_S2373)->differential_0 = _S2188;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2374;
    (&_S2374)->primal_0 = thin_prism_coeffs_21;
    (&_S2374)->differential_0 = _S2188;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2375 = _S2190._S2172;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2365, &_S2366, &_S2367, &_S2368, &_S2369, &_S2370, &_S2372, &_S2373, &_S2374, _S2364, v_mean2d_4, &_S2375);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2376;
    (&_S2376)->primal_0 = _S2184;
    (&_S2376)->differential_0 = _S2301;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2377;
    (&_S2377)->primal_0 = _S2185;
    (&_S2377)->differential_0 = _S2301;
    s_bwd_prop_mul_4(&_S2376, &_S2377, _S2366.differential_0);
    Matrix<float, 3, 3>  _S2378 = transpose_0(_S2377.differential_0 + _S2304.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2379;
    (&_S2379)->primal_0 = R_18;
    (&_S2379)->differential_0 = _S2301;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2380;
    (&_S2380)->primal_0 = _S2183;
    (&_S2380)->differential_0 = _S2301;
    s_bwd_prop_mul_4(&_S2379, &_S2380, _S2376.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2381;
    (&_S2381)->primal_0 = _S2181;
    (&_S2381)->differential_0 = _S2301;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2382;
    (&_S2382)->primal_0 = _S2182;
    (&_S2382)->differential_0 = _S2301;
    s_bwd_prop_mul_4(&_S2381, &_S2382, _S2380.differential_0);
    Matrix<float, 3, 3>  _S2383 = _S2381.differential_0 + transpose_0(_S2382.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2384;
    (&_S2384)->primal_0 = _S2180;
    (&_S2384)->differential_0 = _S2301;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2385;
    (&_S2385)->primal_0 = S_4;
    (&_S2385)->differential_0 = _S2301;
    s_bwd_prop_mul_4(&_S2384, &_S2385, _S2383);
    Matrix<float, 3, 3>  _S2386 = transpose_0(_S2384.differential_0);
    float _S2387 = 2.0f * - _S2386.rows[int(2)].z;
    float _S2388 = 2.0f * _S2386.rows[int(2)].y;
    float _S2389 = 2.0f * _S2386.rows[int(2)].x;
    float _S2390 = 2.0f * _S2386.rows[int(1)].z;
    float _S2391 = 2.0f * - _S2386.rows[int(1)].y;
    float _S2392 = 2.0f * _S2386.rows[int(1)].x;
    float _S2393 = 2.0f * _S2386.rows[int(0)].z;
    float _S2394 = 2.0f * _S2386.rows[int(0)].y;
    float _S2395 = 2.0f * - _S2386.rows[int(0)].x;
    float _S2396 = - _S2392 + _S2394;
    float _S2397 = _S2389 + - _S2393;
    float _S2398 = - _S2388 + _S2390;
    float _S2399 = _S2388 + _S2390;
    float _S2400 = _S2389 + _S2393;
    float _S2401 = _S2392 + _S2394;
    float _S2402 = quat_17.w * (_S2391 + _S2395);
    float _S2403 = quat_17.z * (_S2387 + _S2395);
    float _S2404 = quat_17.y * (_S2387 + _S2391);
    float _S2405 = quat_17.x * _S2396 + quat_17.z * _S2399 + quat_17.y * _S2400 + _S2402 + _S2402;
    float _S2406 = quat_17.x * _S2397 + quat_17.w * _S2399 + quat_17.y * _S2401 + _S2403 + _S2403;
    float _S2407 = quat_17.x * _S2398 + quat_17.w * _S2400 + quat_17.z * _S2401 + _S2404 + _S2404;
    float _S2408 = quat_17.w * _S2396 + quat_17.z * _S2397 + quat_17.y * _S2398;
    float3  _S2409 = _S2239;
    *&((&_S2409)->z) = _S2385.differential_0.rows[int(2)].z;
    *&((&_S2409)->y) = _S2385.differential_0.rows[int(1)].y;
    *&((&_S2409)->x) = _S2385.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2410;
    (&_S2410)->primal_0 = scale_16;
    (&_S2410)->differential_0 = _S2239;
    s_bwd_prop_exp_1(&_S2410, _S2409);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2411;
    (&_S2411)->primal_0 = mean_c_14;
    (&_S2411)->differential_0 = _S2239;
    s_bwd_length_impl_0(&_S2411, 0.0f);
    float3  _S2412 = _S2365.differential_0 + _S2411.differential_0 + _S2338;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2413;
    (&_S2413)->primal_0 = R_18;
    (&_S2413)->differential_0 = _S2301;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2414;
    (&_S2414)->primal_0 = mean_14;
    (&_S2414)->differential_0 = _S2239;
    s_bwd_prop_mul_1(&_S2413, &_S2414, _S2412);
    float3  _S2415 = _S2412 + _S2305.differential_0;
    Matrix<float, 3, 3>  _S2416 = _S2378 + _S2379.differential_0 + _S2413.differential_0;
    float3  _S2417 = _S2410.differential_0 + _S2307;
    float4  _S2418 = _S2371;
    *&((&_S2418)->w) = _S2405;
    *&((&_S2418)->z) = _S2406;
    *&((&_S2418)->y) = _S2407;
    *&((&_S2418)->x) = _S2408;
    float4  _S2419 = _S2418;
    float3  _S2420 = _S2414.differential_0 + _S2299;
    *v_mean_4 = _S2420;
    *v_quat_4 = _S2419;
    *v_scale_4 = _S2417;
    *v_in_opacity_4 = _S2344;
    (*v_sh_coeffs_4)[int(0)] = _S2322;
    (*v_sh_coeffs_4)[int(1)] = _S2323;
    (*v_sh_coeffs_4)[int(2)] = _S2324;
    (*v_sh_coeffs_4)[int(3)] = _S2325;
    (*v_sh_coeffs_4)[int(4)] = _S2326;
    (*v_sh_coeffs_4)[int(5)] = _S2327;
    (*v_sh_coeffs_4)[int(6)] = _S2328;
    (*v_sh_coeffs_4)[int(7)] = _S2329;
    (*v_sh_coeffs_4)[int(8)] = _S2330;
    (*v_sh_coeffs_4)[int(9)] = _S2331;
    (*v_sh_coeffs_4)[int(10)] = _S2332;
    (*v_sh_coeffs_4)[int(11)] = _S2333;
    (*v_sh_coeffs_4)[int(12)] = _S2334;
    (*v_sh_coeffs_4)[int(13)] = _S2335;
    (*v_sh_coeffs_4)[int(14)] = _S2336;
    (*v_sh_coeffs_4)[int(15)] = _S2337;
    *v_R_5 = _S2416;
    *v_t_5 = _S2415;
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
    float _S2421 = (*dpquat_0).primal_0.y;
    float x2_19 = _S2421 * _S2421;
    float y2_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_34 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2422 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_34), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_34), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19))));
    Matrix<float, 3, 3>  _S2423 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2424;
    (&_S2424)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2424)->differential_0 = _S2423;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2425;
    (&_S2425)->primal_0 = _S2422;
    (&_S2425)->differential_0 = _S2423;
    s_bwd_prop_mul_4(&_S2424, &_S2425, _s_dOut_7);
    Matrix<float, 3, 3>  _S2426 = transpose_0(transpose_0(_S2425.differential_0));
    float _S2427 = 2.0f * - _S2426.rows[int(2)].z;
    float _S2428 = 2.0f * _S2426.rows[int(2)].y;
    float _S2429 = 2.0f * _S2426.rows[int(2)].x;
    float _S2430 = 2.0f * _S2426.rows[int(1)].z;
    float _S2431 = 2.0f * - _S2426.rows[int(1)].y;
    float _S2432 = 2.0f * _S2426.rows[int(1)].x;
    float _S2433 = 2.0f * _S2426.rows[int(0)].z;
    float _S2434 = 2.0f * _S2426.rows[int(0)].y;
    float _S2435 = 2.0f * - _S2426.rows[int(0)].x;
    float _S2436 = - _S2432 + _S2434;
    float _S2437 = _S2429 + - _S2433;
    float _S2438 = - _S2428 + _S2430;
    float _S2439 = _S2428 + _S2430;
    float _S2440 = _S2429 + _S2433;
    float _S2441 = _S2432 + _S2434;
    float _S2442 = (*dpquat_0).primal_0.w * (_S2431 + _S2435);
    float _S2443 = (*dpquat_0).primal_0.z * (_S2427 + _S2435);
    float _S2444 = (*dpquat_0).primal_0.y * (_S2427 + _S2431);
    float _S2445 = (*dpquat_0).primal_0.x * _S2436 + (*dpquat_0).primal_0.z * _S2439 + (*dpquat_0).primal_0.y * _S2440 + _S2442 + _S2442;
    float _S2446 = (*dpquat_0).primal_0.x * _S2437 + (*dpquat_0).primal_0.w * _S2439 + (*dpquat_0).primal_0.y * _S2441 + _S2443 + _S2443;
    float _S2447 = (*dpquat_0).primal_0.x * _S2438 + (*dpquat_0).primal_0.w * _S2440 + (*dpquat_0).primal_0.z * _S2441 + _S2444 + _S2444;
    float _S2448 = (*dpquat_0).primal_0.w * _S2436 + (*dpquat_0).primal_0.z * _S2437 + (*dpquat_0).primal_0.y * _S2438;
    float3  _S2449 = make_float3 (_S2424.differential_0.rows[int(0)].x, _S2424.differential_0.rows[int(1)].y, _S2424.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S2449;
    float4  _S2450 = make_float4 (0.0f);
    *&((&_S2450)->w) = _S2445;
    *&((&_S2450)->z) = _S2446;
    *&((&_S2450)->y) = _S2447;
    *&((&_S2450)->x) = _S2448;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S2450;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S2451, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2452, Matrix<float, 3, 3>  _S2453)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S2451, _S2452, _S2453);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_19, float3  scale_18, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S2454 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_19;
    (&dp_quat_0)->differential_0 = _S2454;
    float3  _S2455 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_18;
    (&dp_scale_0)->differential_0 = _S2455;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_15)
{
    float _S2456 = dOut_15.y;
    float _S2457 = dOut_15.z;
    float _S2458 = dOut_15.x;
    float _S2459 = (*a_0).primal_0.z * _S2456 + - (*a_0).primal_0.y * _S2457;
    float _S2460 = - (*a_0).primal_0.z * _S2458 + (*a_0).primal_0.x * _S2457;
    float _S2461 = (*a_0).primal_0.y * _S2458 + - (*a_0).primal_0.x * _S2456;
    float3  _S2462 = make_float3 (- (*b_0).primal_0.z * _S2456 + (*b_0).primal_0.y * _S2457, (*b_0).primal_0.z * _S2458 + - (*b_0).primal_0.x * _S2457, - (*b_0).primal_0.y * _S2458 + (*b_0).primal_0.x * _S2456);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S2462;
    float3  _S2463 = make_float3 (_S2459, _S2460, _S2461);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S2463;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S2464 = left_10.y;
    float _S2465 = right_10.z;
    float _S2466 = left_10.z;
    float _S2467 = right_10.y;
    float _S2468 = right_10.x;
    float _S2469 = left_10.x;
    return make_float3 (_S2464 * _S2465 - _S2466 * _S2467, _S2466 * _S2468 - _S2469 * _S2465, _S2469 * _S2467 - _S2464 * _S2468);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_15, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_1, float3  ray_d_1)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_1);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_1 - mean_15));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S2470, float3  _S2471)
{
    return cross_0(_S2470, _S2471);
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2472, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2473, float3  _S2474)
{
    _d_cross_0(_S2472, _S2473, _S2474);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_8)
{
    float3  _S2475 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S2476 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, _S2475);
    float3  _S2477 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, (*dpray_d_2).primal_0);
    float3  _S2478 = s_primal_ctx_cross_0(_S2477, _S2476);
    float _S2479 = -0.5f * s_primal_ctx_dot_0(_S2478, _S2478);
    float _S2480 = s_primal_ctx_dot_0(_S2477, _S2477);
    float _S2481 = _S2479 / _S2480;
    float _S2482 = _S2480 * _S2480;
    float _S2483 = (*dpopacity_0).primal_0 * _s_dOut_8;
    float _S2484 = s_primal_ctx_exp_1(_S2481) * _s_dOut_8;
    DiffPair_float_0 _S2485;
    (&_S2485)->primal_0 = _S2481;
    (&_S2485)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2485, _S2483);
    float _S2486 = _S2485.differential_0 / _S2482;
    float _S2487 = _S2479 * - _S2486;
    float _S2488 = _S2480 * _S2486;
    float3  _S2489 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2490;
    (&_S2490)->primal_0 = _S2477;
    (&_S2490)->differential_0 = _S2489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2491;
    (&_S2491)->primal_0 = _S2477;
    (&_S2491)->differential_0 = _S2489;
    s_bwd_prop_dot_0(&_S2490, &_S2491, _S2487);
    float _S2492 = -0.5f * _S2488;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2493;
    (&_S2493)->primal_0 = _S2478;
    (&_S2493)->differential_0 = _S2489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2494;
    (&_S2494)->primal_0 = _S2478;
    (&_S2494)->differential_0 = _S2489;
    s_bwd_prop_dot_0(&_S2493, &_S2494, _S2492);
    float3  _S2495 = _S2494.differential_0 + _S2493.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2496;
    (&_S2496)->primal_0 = _S2477;
    (&_S2496)->differential_0 = _S2489;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2497;
    (&_S2497)->primal_0 = _S2476;
    (&_S2497)->differential_0 = _S2489;
    s_bwd_prop_cross_0(&_S2496, &_S2497, _S2495);
    float3  _S2498 = _S2491.differential_0 + _S2490.differential_0 + _S2496.differential_0;
    Matrix<float, 3, 3>  _S2499 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2500;
    (&_S2500)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2500)->differential_0 = _S2499;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2501;
    (&_S2501)->primal_0 = (*dpray_d_2).primal_0;
    (&_S2501)->differential_0 = _S2489;
    s_bwd_prop_mul_1(&_S2500, &_S2501, _S2498);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2502;
    (&_S2502)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2502)->differential_0 = _S2499;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2503;
    (&_S2503)->primal_0 = _S2475;
    (&_S2503)->differential_0 = _S2489;
    s_bwd_prop_mul_1(&_S2502, &_S2503, _S2497.differential_0);
    float3  _S2504 = - _S2503.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S2501.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S2503.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S2484;
    Matrix<float, 3, 3>  _S2505 = _S2500.differential_0 + _S2502.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S2505;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S2504;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2506, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S2507, DiffPair_float_0 * _S2508, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2509, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2510, float _S2511)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S2506, _S2507, _S2508, _S2509, _S2510, _S2511);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S2512 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_16;
    (&dp_mean_0)->differential_0 = _S2512;
    Matrix<float, 3, 3>  _S2513 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S2513;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S2512;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S2512;
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
    float3  _S2514 = mean_17 - ray_o_3;
    *depth_10 = 0.5f * (F32_log((dot_0(_S2514, _S2514) + 9.99999997475242708e-07f)));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S2515 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    float _S2516 = 0.5f * dpdepth_0;
    DiffPair_float_0 _S2517;
    (&_S2517)->primal_0 = s_primal_ctx_dot_0(_S2515, _S2515) + 9.99999997475242708e-07f;
    (&_S2517)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2517, _S2516);
    float3  _S2518 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2519;
    (&_S2519)->primal_0 = _S2515;
    (&_S2519)->differential_0 = _S2518;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2520;
    (&_S2520)->primal_0 = _S2515;
    (&_S2520)->differential_0 = _S2518;
    s_bwd_prop_dot_0(&_S2519, &_S2520, _S2517.differential_0);
    float3  _S2521 = _S2520.differential_0 + _S2519.differential_0;
    float3  _S2522 = - _S2521;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S2518;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S2522;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S2523 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S2523;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S2521;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2524, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S2525, DiffPair_float_0 * _S2526, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2527, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2528, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2529, float3  _S2530, float _S2531)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S2524, _S2525, _S2526, _S2527, _S2528, _S2529, _S2530, _S2531);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S2532 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_18;
    (&dp_mean_1)->differential_0 = _S2532;
    Matrix<float, 3, 3>  _S2533 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S2533;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S2532;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S2532;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S2532;
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
    float _S2534 = _slang_select(((*dpx_13).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_13).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_16;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S2534;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_14, float dOut_17)
{
    float _S2535 = (F32_exp2(((*dpx_14).primal_0))) * 50.693145751953125f * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S2535;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float3  dOut_18)
{
    float3  _S2536 = make_float3 (1.0f) / (*dpx_15).primal_0 * dOut_18;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S2536;
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

inline __device__ void projection_opaque_triangle_persp(float3  mean_19, float4  quat_20, float3  scale_19, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_19, float3  t_18, float fx_19, float fy_19, float cx_19, float cy_19, float4  radial_coeffs_22, float2  tangential_coeffs_22, float2  thin_prism_coeffs_22, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_19, mean_19) + t_18;
        float _S2537 = mean_c_15.z;
        bool _S2538;
        if(_S2537 < near_plane_10)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = _S2537 > far_plane_10;
        }
        if(_S2538)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S2539 = scale_19.x;
        float sx_0 = (F32_exp((_S2539)));
        float _S2540 = scale_19.y;
        float sy_0 = (F32_exp((_S2540)));
        float sz_0 = scale_19.z - 0.5f * (_S2539 + _S2540);
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
        Matrix<float, 3, 3>  _S2541 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20)));
        float3  vert0_c_0 = mul_0(R_19, mul_0(_S2541, make_float3 (sx_0, 0.0f, 0.0f)) + mean_19) + t_18;
        float3  vert1_c_0 = mul_0(R_19, mul_0(_S2541, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_19) + t_18;
        float3  vert2_c_0 = mul_0(R_19, mul_0(_S2541, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_19) + t_18;
        float _S2542 = vert0_c_0.z;
        float _S2543 = vert1_c_0.z;
        float _S2544 = vert2_c_0.z;
        if(_S2542 < near_plane_10)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = _S2542 > far_plane_10;
        }
        if(_S2538)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = _S2543 < near_plane_10;
        }
        if(_S2538)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = _S2543 > far_plane_10;
        }
        if(_S2538)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = _S2544 < near_plane_10;
        }
        if(_S2538)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = _S2544 > far_plane_10;
        }
        if(_S2538)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S2542);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S2543);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S2544);
        float2  _S2545 = make_float2 (fx_19, fy_19);
        float2  _S2546 = make_float2 (cx_19, cy_19);
        *uv0_0 = _S2545 * *uv0_0 + _S2546;
        *uv1_0 = _S2545 * *uv1_0 + _S2546;
        float2  _S2547 = _S2545 * *uv2_0 + _S2546;
        *uv2_0 = _S2547;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S2547 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S2547)));
        float _S2548 = _S2547.x;
        float xmax_5 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S2548))) + offset_0;
        float xmin_5 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S2548))) - offset_0;
        float _S2549 = _S2547.y;
        float ymax_5 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S2549))) + offset_0;
        float ymin_5 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S2549))) - offset_0;
        if(xmax_5 <= 0.0f)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = xmin_5 >= float(image_width_15);
        }
        if(_S2538)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = ymax_5 <= 0.0f;
        }
        if(_S2538)
        {
            _S2538 = true;
        }
        else
        {
            _S2538 = ymin_5 >= float(image_height_15);
        }
        if(_S2538)
        {
            _S2538 = true;
        }
        else
        {
            if(_S2537 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S2538 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S2538 = false;
                }
                if(_S2538)
                {
                    _S2538 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S2538 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S2538 = false;
                    }
                }
            }
            else
            {
                _S2538 = false;
            }
        }
        if(_S2538)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S2550 = mean_19 - - mul_0(transpose_0(R_19), t_18);
        float _S2551 = _S2550.x;
        float _S2552 = _S2550.y;
        float _S2553 = _S2550.z;
        float norm_10 = (F32_sqrt((_S2551 * _S2551 + _S2552 * _S2552 + _S2553 * _S2553)));
        float x_42 = _S2551 / norm_10;
        float y_18 = _S2552 / norm_10;
        float z_15 = _S2553 / norm_10;
        float z2_36 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_42 * x_42 - y_18 * y_18;
        float fS1_15 = 2.0f * x_42 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_36 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_42) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_36 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_42) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_42 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_36 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_42) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_42 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S2554 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S2554);
        float3  _S2555 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S2556 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S2555 + _S2556 + make_float3 (0.5f), _S2554);
        (*rgb_12)[int(2)] = max_0(_S2555 - _S2556 + make_float3 (0.5f), _S2554);
        float3  _S2557 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S2557 * make_float3 (float(- (F32_sign((dot_0(_S2557, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_20, float3  t_19, float fx_20, float fy_20, float cx_20, float cy_20, float4  radial_coeffs_23, float2  tangential_coeffs_23, float2  thin_prism_coeffs_23, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_20, mean_20) + t_19;
        float _S2558 = length_1(mean_c_16);
        bool _S2559;
        if(_S2558 < near_plane_11)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = _S2558 > far_plane_11;
        }
        if(_S2559)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S2560 = scale_20.x;
        float sx_1 = (F32_exp((_S2560)));
        float _S2561 = scale_20.y;
        float sy_1 = (F32_exp((_S2561)));
        float sz_1 = scale_20.z - 0.5f * (_S2560 + _S2561);
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
        Matrix<float, 3, 3>  _S2562 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_37), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_37), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_1 = mul_0(R_20, mul_0(_S2562, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_19;
        float3  vert1_c_1 = mul_0(R_20, mul_0(_S2562, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_19;
        float3  vert2_c_1 = mul_0(R_20, mul_0(_S2562, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_19;
        float _S2563 = length_1(vert0_c_1);
        float _S2564 = length_1(vert1_c_1);
        float _S2565 = length_1(vert2_c_1);
        if(_S2563 < near_plane_11)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = _S2563 > far_plane_11;
        }
        if(_S2559)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = _S2564 < near_plane_11;
        }
        if(_S2559)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = _S2564 > far_plane_11;
        }
        if(_S2559)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = _S2565 < near_plane_11;
        }
        if(_S2559)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = _S2565 > far_plane_11;
        }
        if(_S2559)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_2 = CameraDistortion_x24init_0(radial_coeffs_23, tangential_coeffs_23, thin_prism_coeffs_23);
        float2  _S2566 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_7 = length_0(_S2566);
        float _S2567 = vert0_c_1.z;
        float theta_2 = (F32_atan2((r_7), (_S2567)));
        float k_4;
        if(theta_2 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_2 * theta_2 / 3.0f) / _S2567;
        }
        else
        {
            k_4 = theta_2 / r_7;
        }
        float2  _S2568 = _S2566 * make_float2 (k_4);
        float k1_3 = dist_coeffs_2.radial_coeffs_0.x;
        float k2_3 = dist_coeffs_2.radial_coeffs_0.y;
        float k3_3 = dist_coeffs_2.radial_coeffs_0.z;
        float k4_4 = dist_coeffs_2.radial_coeffs_0.w;
        float p1_4 = dist_coeffs_2.tangential_coeffs_0.x;
        float p2_4 = dist_coeffs_2.tangential_coeffs_0.y;
        float sx1_4 = dist_coeffs_2.thin_prism_coeffs_0.x;
        float sy1_4 = dist_coeffs_2.thin_prism_coeffs_0.y;
        float u_10 = _S2568.x;
        float v_10 = _S2568.y;
        float r2_10 = u_10 * u_10 + v_10 * v_10;
        float _S2569 = 2.0f * p1_4;
        float _S2570 = 2.0f * p2_4;
        float2  _S2571 = _S2568 * make_float2 (1.0f + r2_10 * (k1_3 + r2_10 * (k2_3 + r2_10 * (k3_3 + r2_10 * k4_4)))) + make_float2 (_S2569 * u_10 * v_10 + p2_4 * (r2_10 + 2.0f * u_10 * u_10) + sx1_4 * r2_10, _S2570 * u_10 * v_10 + p1_4 * (r2_10 + 2.0f * v_10 * v_10) + sy1_4 * r2_10);
        *uv0_1 = make_float2 (fx_20 * _S2571.x + cx_20, fy_20 * _S2571.y + cy_20);
        float2  _S2572 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_8 = length_0(_S2572);
        float _S2573 = vert1_c_1.z;
        float theta_3 = (F32_atan2((r_8), (_S2573)));
        if(theta_3 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_3 * theta_3 / 3.0f) / _S2573;
        }
        else
        {
            k_4 = theta_3 / r_8;
        }
        float2  _S2574 = _S2572 * make_float2 (k_4);
        float u_11 = _S2574.x;
        float v_11 = _S2574.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float2  _S2575 = _S2574 * make_float2 (1.0f + r2_11 * (k1_3 + r2_11 * (k2_3 + r2_11 * (k3_3 + r2_11 * k4_4)))) + make_float2 (_S2569 * u_11 * v_11 + p2_4 * (r2_11 + 2.0f * u_11 * u_11) + sx1_4 * r2_11, _S2570 * u_11 * v_11 + p1_4 * (r2_11 + 2.0f * v_11 * v_11) + sy1_4 * r2_11);
        *uv1_1 = make_float2 (fx_20 * _S2575.x + cx_20, fy_20 * _S2575.y + cy_20);
        float2  _S2576 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_9 = length_0(_S2576);
        float _S2577 = vert2_c_1.z;
        float theta_4 = (F32_atan2((r_9), (_S2577)));
        if(theta_4 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_4 * theta_4 / 3.0f) / _S2577;
        }
        else
        {
            k_4 = theta_4 / r_9;
        }
        float2  _S2578 = _S2576 * make_float2 (k_4);
        float u_12 = _S2578.x;
        float v_12 = _S2578.y;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float2  _S2579 = _S2578 * make_float2 (1.0f + r2_12 * (k1_3 + r2_12 * (k2_3 + r2_12 * (k3_3 + r2_12 * k4_4)))) + make_float2 (_S2569 * u_12 * v_12 + p2_4 * (r2_12 + 2.0f * u_12 * u_12) + sx1_4 * r2_12, _S2570 * u_12 * v_12 + p1_4 * (r2_12 + 2.0f * v_12 * v_12) + sy1_4 * r2_12);
        float _S2580 = fx_20 * _S2579.x + cx_20;
        float _S2581 = fy_20 * _S2579.y + cy_20;
        float2  _S2582 = make_float2 (_S2580, _S2581);
        *uv2_1 = _S2582;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S2582 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S2582)));
        float xmax_6 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S2580))) + offset_1;
        float xmin_6 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S2580))) - offset_1;
        float ymax_6 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S2581))) + offset_1;
        float ymin_6 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S2581))) - offset_1;
        if(xmax_6 <= 0.0f)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = xmin_6 >= float(image_width_16);
        }
        if(_S2559)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = ymax_6 <= 0.0f;
        }
        if(_S2559)
        {
            _S2559 = true;
        }
        else
        {
            _S2559 = ymin_6 >= float(image_height_16);
        }
        if(_S2559)
        {
            _S2559 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S2559 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S2559 = false;
                }
                if(_S2559)
                {
                    _S2559 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S2559 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S2559 = false;
                    }
                }
            }
            else
            {
                _S2559 = false;
            }
        }
        if(_S2559)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S2563, _S2564, _S2565) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S2583 = mean_20 - - mul_0(transpose_0(R_20), t_19);
        float _S2584 = _S2583.x;
        float _S2585 = _S2583.y;
        float _S2586 = _S2583.z;
        float norm_11 = (F32_sqrt((_S2584 * _S2584 + _S2585 * _S2585 + _S2586 * _S2586)));
        float x_44 = _S2584 / norm_11;
        float y_19 = _S2585 / norm_11;
        float z_16 = _S2586 / norm_11;
        float z2_38 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_44 * x_44 - y_19 * y_19;
        float fS1_16 = 2.0f * x_44 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_38 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_44) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_38 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_44) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_44 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_38 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_44) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_44 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S2587 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S2587);
        float3  _S2588 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S2589 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S2588 + _S2589 + make_float3 (0.5f), _S2587);
        (*rgb_13)[int(2)] = max_0(_S2588 - _S2589 + make_float3 (0.5f), _S2587);
        float3  _S2590 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S2590 * make_float3 (float(- (F32_sign((dot_0(_S2590, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_21, float3  t_20, float fx_21, float fy_21, float cx_21, float cy_21, float4  radial_coeffs_24, float2  tangential_coeffs_24, float2  thin_prism_coeffs_24, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_21, mean_21) + t_20;
    float _S2591 = scale_21.x;
    float sx_2 = (F32_exp((_S2591)));
    float _S2592 = scale_21.y;
    float sy_2 = (F32_exp((_S2592)));
    float sz_2 = scale_21.z - 0.5f * (_S2591 + _S2592);
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
    Matrix<float, 3, 3>  _S2593 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_39), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_39), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
    float3  vert0_c_2 = mul_0(R_21, mul_0(_S2593, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_20;
    float3  vert1_c_2 = mul_0(R_21, mul_0(_S2593, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_20;
    float3  vert2_c_2 = mul_0(R_21, mul_0(_S2593, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_20;
    *uv0_2 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    *uv1_2 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    *uv2_2 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float2  _S2594 = make_float2 (fx_21, fy_21);
    float2  _S2595 = make_float2 (cx_21, cy_21);
    *uv0_2 = _S2594 * *uv0_2 + _S2595;
    *uv1_2 = _S2594 * *uv1_2 + _S2595;
    float2  _S2596 = _S2594 * *uv2_2 + _S2595;
    *uv2_2 = _S2596;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S2596 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S2596)));
    float _S2597 = _S2596.x;
    float _S2598 = _S2596.y;
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S2597))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S2598))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S2597))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S2598))) + offset_2)))));
    *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_2 = hardness_2;
    float3  _S2599 = mean_21 - - mul_0(transpose_0(R_21), t_20);
    float _S2600 = _S2599.x;
    float _S2601 = _S2599.y;
    float _S2602 = _S2599.z;
    float norm_12 = (F32_sqrt((_S2600 * _S2600 + _S2601 * _S2601 + _S2602 * _S2602)));
    float x_46 = _S2600 / norm_12;
    float y_20 = _S2601 / norm_12;
    float z_17 = _S2602 / norm_12;
    float z2_40 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_46 * x_46 - y_20 * y_20;
    float fS1_17 = 2.0f * x_46 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_40 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_46) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_40 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_46) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_40 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_46) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S2603 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S2603);
    float3  _S2604 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S2605 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S2604 + _S2605 + make_float3 (0.5f), _S2603);
    (*rgb_14)[int(2)] = max_0(_S2604 - _S2605 + make_float3 (0.5f), _S2603);
    float3  _S2606 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S2606 * make_float3 (float(- (F32_sign((dot_0(_S2606, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_22, float3  t_21, float fx_22, float fy_22, float cx_22, float cy_22, float4  radial_coeffs_25, float2  tangential_coeffs_25, float2  thin_prism_coeffs_25, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_22, mean_22) + t_21;
    float _S2607 = scale_22.x;
    float sx_3 = (F32_exp((_S2607)));
    float _S2608 = scale_22.y;
    float sy_3 = (F32_exp((_S2608)));
    float sz_3 = scale_22.z - 0.5f * (_S2607 + _S2608);
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
    Matrix<float, 3, 3>  _S2609 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_41), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_41), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_3 = mul_0(R_22, mul_0(_S2609, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_21;
    float3  vert1_c_3 = mul_0(R_22, mul_0(_S2609, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_21;
    float3  vert2_c_3 = mul_0(R_22, mul_0(_S2609, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_21;
    CameraDistortion_0 dist_coeffs_3 = CameraDistortion_x24init_0(radial_coeffs_25, tangential_coeffs_25, thin_prism_coeffs_25);
    float2  _S2610 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_10 = length_0(_S2610);
    float _S2611 = vert0_c_3.z;
    float theta_5 = (F32_atan2((r_10), (_S2611)));
    float k_5;
    if(theta_5 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_5 * theta_5 / 3.0f) / _S2611;
    }
    else
    {
        k_5 = theta_5 / r_10;
    }
    float2  _S2612 = _S2610 * make_float2 (k_5);
    float k1_4 = dist_coeffs_3.radial_coeffs_0.x;
    float k2_4 = dist_coeffs_3.radial_coeffs_0.y;
    float k3_4 = dist_coeffs_3.radial_coeffs_0.z;
    float k4_5 = dist_coeffs_3.radial_coeffs_0.w;
    float p1_5 = dist_coeffs_3.tangential_coeffs_0.x;
    float p2_5 = dist_coeffs_3.tangential_coeffs_0.y;
    float sx1_5 = dist_coeffs_3.thin_prism_coeffs_0.x;
    float sy1_5 = dist_coeffs_3.thin_prism_coeffs_0.y;
    float u_13 = _S2612.x;
    float v_13 = _S2612.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float _S2613 = 2.0f * p1_5;
    float _S2614 = 2.0f * p2_5;
    float2  _S2615 = _S2612 * make_float2 (1.0f + r2_13 * (k1_4 + r2_13 * (k2_4 + r2_13 * (k3_4 + r2_13 * k4_5)))) + make_float2 (_S2613 * u_13 * v_13 + p2_5 * (r2_13 + 2.0f * u_13 * u_13) + sx1_5 * r2_13, _S2614 * u_13 * v_13 + p1_5 * (r2_13 + 2.0f * v_13 * v_13) + sy1_5 * r2_13);
    *uv0_3 = make_float2 (fx_22 * _S2615.x + cx_22, fy_22 * _S2615.y + cy_22);
    float2  _S2616 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_11 = length_0(_S2616);
    float _S2617 = vert1_c_3.z;
    float theta_6 = (F32_atan2((r_11), (_S2617)));
    if(theta_6 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_6 * theta_6 / 3.0f) / _S2617;
    }
    else
    {
        k_5 = theta_6 / r_11;
    }
    float2  _S2618 = _S2616 * make_float2 (k_5);
    float u_14 = _S2618.x;
    float v_14 = _S2618.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S2619 = _S2618 * make_float2 (1.0f + r2_14 * (k1_4 + r2_14 * (k2_4 + r2_14 * (k3_4 + r2_14 * k4_5)))) + make_float2 (_S2613 * u_14 * v_14 + p2_5 * (r2_14 + 2.0f * u_14 * u_14) + sx1_5 * r2_14, _S2614 * u_14 * v_14 + p1_5 * (r2_14 + 2.0f * v_14 * v_14) + sy1_5 * r2_14);
    *uv1_3 = make_float2 (fx_22 * _S2619.x + cx_22, fy_22 * _S2619.y + cy_22);
    float2  _S2620 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_12 = length_0(_S2620);
    float _S2621 = vert2_c_3.z;
    float theta_7 = (F32_atan2((r_12), (_S2621)));
    if(theta_7 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_7 * theta_7 / 3.0f) / _S2621;
    }
    else
    {
        k_5 = theta_7 / r_12;
    }
    float2  _S2622 = _S2620 * make_float2 (k_5);
    float u_15 = _S2622.x;
    float v_15 = _S2622.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S2623 = _S2622 * make_float2 (1.0f + r2_15 * (k1_4 + r2_15 * (k2_4 + r2_15 * (k3_4 + r2_15 * k4_5)))) + make_float2 (_S2613 * u_15 * v_15 + p2_5 * (r2_15 + 2.0f * u_15 * u_15) + sx1_5 * r2_15, _S2614 * u_15 * v_15 + p1_5 * (r2_15 + 2.0f * v_15 * v_15) + sy1_5 * r2_15);
    float _S2624 = fx_22 * _S2623.x + cx_22;
    float _S2625 = fy_22 * _S2623.y + cy_22;
    float2  _S2626 = make_float2 (_S2624, _S2625);
    *uv2_3 = _S2626;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S2626 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S2626)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S2624))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S2625))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S2624))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S2625))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S2627 = mean_22 - - mul_0(transpose_0(R_22), t_21);
    float _S2628 = _S2627.x;
    float _S2629 = _S2627.y;
    float _S2630 = _S2627.z;
    float norm_13 = (F32_sqrt((_S2628 * _S2628 + _S2629 * _S2629 + _S2630 * _S2630)));
    float x_48 = _S2628 / norm_13;
    float y_21 = _S2629 / norm_13;
    float z_18 = _S2630 / norm_13;
    float z2_42 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_48 * x_48 - y_21 * y_21;
    float fS1_18 = 2.0f * x_48 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_42 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_48) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_42 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_48) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_48 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_42 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_48) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_48 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S2631 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S2631);
    float3  _S2632 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S2633 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S2632 + _S2633 + make_float3 (0.5f), _S2631);
    (*rgb_15)[int(2)] = max_0(_S2632 - _S2633 + make_float3 (0.5f), _S2631);
    float3  _S2634 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S2634 * make_float3 (float(- (F32_sign((dot_0(_S2634, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2635, float3  _S2636)
{
    _d_log_vector_0(_S2635, _S2636);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S2637, float _S2638)
{
    _d_exp2_0(_S2637, _S2638);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S2639, float _S2640)
{
    _d_abs_0(_S2639, _S2640);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_23, float3  t_22, float fx_23, float fy_23, float cx_23, float cy_23, float4  radial_coeffs_26, float2  tangential_coeffs_26, float2  thin_prism_coeffs_26, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_1(R_23, mean_23) + t_22;
    float _S2641 = scale_23.x;
    float _S2642 = s_primal_ctx_exp_1(_S2641);
    float _S2643 = scale_23.y;
    float _S2644 = s_primal_ctx_exp_1(_S2643);
    float sz_4 = scale_23.z - 0.5f * (_S2641 + _S2643);
    float _S2645 = quat_24.y;
    float x2_24 = _S2645 * _S2645;
    float y2_24 = quat_24.z * quat_24.z;
    float z2_43 = quat_24.w * quat_24.w;
    float xy_24 = quat_24.y * quat_24.z;
    float xz_24 = quat_24.y * quat_24.w;
    float yz_24 = quat_24.z * quat_24.w;
    float wx_24 = quat_24.x * quat_24.y;
    float wy_24 = quat_24.x * quat_24.z;
    float wz_24 = quat_24.x * quat_24.w;
    Matrix<float, 3, 3>  _S2646 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_43), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_43), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  _S2647 = make_float3 (_S2642, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_1(_S2646, _S2647) + mean_23;
    float _S2648 = -0.5f + sz_4;
    float3  _S2649 = make_float3 (_S2642 * _S2648, _S2644, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_1(_S2646, _S2649) + mean_23;
    float _S2650 = -0.5f - sz_4;
    float3  _S2651 = make_float3 (_S2642 * _S2650, - _S2644, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_1(_S2646, _S2651) + mean_23;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_23, vert0_0) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_23, vert1_0) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_23, vert2_0) + t_22;
    float2  _S2652 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S2653 = vert0_c_4.z;
    float2  _S2654 = make_float2 (_S2653);
    float2  _S2655 = make_float2 (_S2653 * _S2653);
    float2  _S2656 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S2657 = vert1_c_4.z;
    float2  _S2658 = make_float2 (_S2657);
    float2  _S2659 = make_float2 (_S2657 * _S2657);
    float2  _S2660 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S2661 = vert2_c_4.z;
    float2  _S2662 = make_float2 (_S2661);
    float2  _S2663 = make_float2 (_S2661 * _S2661);
    float2  _S2664 = make_float2 (fx_23, fy_23);
    float2  _S2665 = make_float2 (cx_23, cy_23);
    float2  _S2666 = _S2664 * (_S2652 / make_float2 (_S2653)) + _S2665;
    float2  _S2667 = _S2664 * (_S2656 / make_float2 (_S2657)) + _S2665;
    float2  _S2668 = _S2664 * (_S2660 / make_float2 (_S2661)) + _S2665;
    float2  e0_4 = _S2667 - _S2666;
    float2  e1_4 = _S2668 - _S2667;
    float2  e2_0 = _S2666 - _S2668;
    float _S2669 = e0_4.x;
    float _S2670 = e1_4.y;
    float _S2671 = e0_4.y;
    float _S2672 = e1_4.x;
    float _S2673 = _S2669 * _S2670 - _S2671 * _S2672;
    float _S2674 = 1.0f - hardness_4.y;
    float _S2675 = -1.0f / _S2674;
    float _S2676 = _S2674 * _S2674;
    float _S2677 = _S2666.x;
    float _S2678 = _S2667.x;
    float _S2679 = s_primal_ctx_max_0(_S2677, _S2678);
    float _S2680 = _S2668.x;
    float _S2681 = s_primal_ctx_min_0(_S2677, _S2678);
    float _S2682 = _S2666.y;
    float _S2683 = _S2667.y;
    float _S2684 = s_primal_ctx_max_0(_S2682, _S2683);
    float _S2685 = _S2668.y;
    float _S2686 = s_primal_ctx_min_0(_S2682, _S2683);
    float3  _S2687 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S2688 = transpose_0(R_23);
    float3  _S2689 = mean_23 - - s_primal_ctx_mul_1(_S2688, t_22);
    float _S2690 = _S2689.x;
    float _S2691 = _S2689.y;
    float _S2692 = _S2689.z;
    float _S2693 = _S2690 * _S2690 + _S2691 * _S2691 + _S2692 * _S2692;
    float _S2694 = s_primal_ctx_sqrt_0(_S2693);
    float x_49 = _S2690 / _S2694;
    float3  _S2695 = make_float3 (x_49);
    float _S2696 = _S2694 * _S2694;
    float y_22 = _S2691 / _S2694;
    float z_19 = _S2692 / _S2694;
    float3  _S2697 = make_float3 (z_19);
    float _S2698 = - y_22;
    float3  _S2699 = make_float3 (_S2698);
    float z2_44 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_49 * x_49 - y_22 * y_22;
    float _S2700 = 2.0f * x_49;
    float fS1_19 = _S2700 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_44 - 0.31539157032966614f;
    float3  _S2701 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_49;
    float3  _S2702 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S2703 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S2704 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S2705 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_44 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S2706 = 1.86588168144226074f * z2_44 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S2706;
    float3  _S2707 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_49;
    float3  _S2708 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S2709 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S2710 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S2711 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_49 * fC1_19 - y_22 * fS1_19);
    float3  _S2712 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_49 * fS1_19 + y_22 * fC1_19);
    float3  _S2713 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2698) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_49) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S2714 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S2715 = make_float3 (0.0f);
    float3  _S2716 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S2717 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S2718 = make_float3 (_S2717);
    float3  _S2719 = (*ch_coeffs_4)[int(1)] * make_float3 (_S2717);
    float3  _S2720 = _S2716 + _S2719 + make_float3 (0.5f);
    float3  _S2721 = _S2716 - _S2719 + make_float3 (0.5f);
    float3  _S2722 = vert1_c_4 - vert0_c_4;
    float3  _S2723 = vert2_c_4 - vert0_c_4;
    float3  _S2724 = s_primal_ctx_cross_0(_S2722, _S2723);
    float3  _S2725 = normalize_0(_S2724);
    float3  _S2726 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S2725, mean_c_19)))))) * v_normal_0;
    float3  _S2727 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2728;
    (&_S2728)->primal_0 = _S2725;
    (&_S2728)->differential_0 = _S2727;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2729;
    (&_S2729)->primal_0 = mean_c_19;
    (&_S2729)->differential_0 = _S2727;
    s_bwd_prop_dot_0(&_S2728, &_S2729, 0.0f);
    float3  _S2730 = _S2726 + _S2728.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2731;
    (&_S2731)->primal_0 = _S2724;
    (&_S2731)->differential_0 = _S2727;
    s_bwd_normalize_impl_0(&_S2731, _S2730);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2732;
    (&_S2732)->primal_0 = _S2722;
    (&_S2732)->differential_0 = _S2727;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2733;
    (&_S2733)->primal_0 = _S2723;
    (&_S2733)->differential_0 = _S2727;
    s_bwd_prop_cross_0(&_S2732, &_S2733, _S2731.differential_0);
    float3  _S2734 = - _S2733.differential_0;
    float3  _S2735 = - _S2732.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2736;
    (&_S2736)->primal_0 = _S2721;
    (&_S2736)->differential_0 = _S2727;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2737;
    (&_S2737)->primal_0 = _S2715;
    (&_S2737)->differential_0 = _S2727;
    s_bwd_prop_max_0(&_S2736, &_S2737, (*v_rgb_6)[int(2)]);
    float3  _S2738 = - _S2736.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2739;
    (&_S2739)->primal_0 = _S2720;
    (&_S2739)->differential_0 = _S2727;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2740;
    (&_S2740)->primal_0 = _S2715;
    (&_S2740)->differential_0 = _S2727;
    s_bwd_prop_max_0(&_S2739, &_S2740, (*v_rgb_6)[int(1)]);
    float3  _S2741 = _S2718 * (_S2738 + _S2739.differential_0);
    float3  _S2742 = _S2736.differential_0 + _S2739.differential_0;
    float3  _S2743 = make_float3 (0.5f) * - _S2742;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2744;
    (&_S2744)->primal_0 = _S2714;
    (&_S2744)->differential_0 = _S2727;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2745;
    (&_S2745)->primal_0 = _S2715;
    (&_S2745)->differential_0 = _S2727;
    s_bwd_prop_max_0(&_S2744, &_S2745, (*v_rgb_6)[int(0)]);
    float3  _S2746 = _S2743 + _S2744.differential_0;
    float3  _S2747 = _S2742 + _S2744.differential_0;
    float3  _S2748 = _S2712 * _S2747;
    float3  _S2749 = (*sh_coeffs_19)[int(15)] * _S2747;
    float3  _S2750 = _S2710 * _S2747;
    float3  _S2751 = (*sh_coeffs_19)[int(14)] * _S2747;
    float3  _S2752 = _S2708 * _S2747;
    float3  _S2753 = (*sh_coeffs_19)[int(13)] * _S2747;
    float3  _S2754 = _S2707 * _S2747;
    float3  _S2755 = (*sh_coeffs_19)[int(12)] * _S2747;
    float3  _S2756 = _S2709 * _S2747;
    float3  _S2757 = (*sh_coeffs_19)[int(11)] * _S2747;
    float3  _S2758 = _S2711 * _S2747;
    float3  _S2759 = (*sh_coeffs_19)[int(10)] * _S2747;
    float3  _S2760 = _S2713 * _S2747;
    float3  _S2761 = (*sh_coeffs_19)[int(9)] * _S2747;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S2761.x + _S2761.y + _S2761.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S2749.x + _S2749.y + _S2749.z);
    float _S2762 = _S2759.x + _S2759.y + _S2759.z;
    float _S2763 = _S2751.x + _S2751.y + _S2751.z;
    float _S2764 = _S2757.x + _S2757.y + _S2757.z;
    float _S2765 = _S2753.x + _S2753.y + _S2753.z;
    float _S2766 = _S2755.x + _S2755.y + _S2755.z;
    float _S2767 = - s_diff_fC2_T_5;
    float3  _S2768 = _S2704 * _S2747;
    float3  _S2769 = (*sh_coeffs_19)[int(8)] * _S2747;
    float3  _S2770 = _S2702 * _S2747;
    float3  _S2771 = (*sh_coeffs_19)[int(7)] * _S2747;
    float3  _S2772 = _S2701 * _S2747;
    float3  _S2773 = (*sh_coeffs_19)[int(6)] * _S2747;
    float3  _S2774 = _S2703 * _S2747;
    float3  _S2775 = (*sh_coeffs_19)[int(5)] * _S2747;
    float3  _S2776 = _S2705 * _S2747;
    float3  _S2777 = (*sh_coeffs_19)[int(4)] * _S2747;
    float _S2778 = _S2775.x + _S2775.y + _S2775.z;
    float _S2779 = _S2771.x + _S2771.y + _S2771.z;
    float _S2780 = fTmp1B_19 * _S2762 + x_49 * s_diff_fS2_T_5 + y_22 * _S2767 + 0.54627424478530884f * (_S2777.x + _S2777.y + _S2777.z);
    float _S2781 = fTmp1B_19 * _S2763 + y_22 * s_diff_fS2_T_5 + x_49 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S2769.x + _S2769.y + _S2769.z);
    float _S2782 = y_22 * - _S2781;
    float _S2783 = x_49 * _S2781;
    float _S2784 = z_19 * (1.86588168144226074f * (z_19 * _S2766) + -2.28522896766662598f * (y_22 * _S2764 + x_49 * _S2765) + 0.94617468118667603f * (_S2773.x + _S2773.y + _S2773.z));
    float3  _S2785 = make_float3 (0.48860251903533936f) * _S2747;
    float3  _S2786 = - _S2785;
    float3  _S2787 = _S2695 * _S2786;
    float3  _S2788 = (*sh_coeffs_19)[int(3)] * _S2786;
    float3  _S2789 = _S2697 * _S2785;
    float3  _S2790 = (*sh_coeffs_19)[int(2)] * _S2785;
    float3  _S2791 = _S2699 * _S2785;
    float3  _S2792 = (*sh_coeffs_19)[int(1)] * _S2785;
    float _S2793 = (_S2706 * _S2766 + 1.44530570507049561f * (fS1_19 * _S2762 + fC1_19 * _S2763) + -1.09254848957061768f * (y_22 * _S2778 + x_49 * _S2779) + _S2784 + _S2784 + _S2790.x + _S2790.y + _S2790.z) / _S2696;
    float _S2794 = _S2694 * _S2793;
    float _S2795 = (fTmp0C_19 * _S2764 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S2767 + fTmp0B_19 * _S2778 + _S2700 * _S2780 + _S2782 + _S2782 + - (_S2792.x + _S2792.y + _S2792.z)) / _S2696;
    float _S2796 = _S2694 * _S2795;
    float _S2797 = (fTmp0C_19 * _S2765 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S2779 + 2.0f * (y_22 * _S2780) + _S2783 + _S2783 + _S2788.x + _S2788.y + _S2788.z) / _S2696;
    float _S2798 = _S2694 * _S2797;
    float _S2799 = _S2692 * - _S2793 + _S2691 * - _S2795 + _S2690 * - _S2797;
    DiffPair_float_0 _S2800;
    (&_S2800)->primal_0 = _S2693;
    (&_S2800)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2800, _S2799);
    float _S2801 = _S2692 * _S2800.differential_0;
    float _S2802 = _S2691 * _S2800.differential_0;
    float _S2803 = _S2690 * _S2800.differential_0;
    float3  _S2804 = make_float3 (0.282094806432724f) * _S2747;
    float3  _S2805 = make_float3 (_S2798 + _S2803 + _S2803, _S2796 + _S2802 + _S2802, _S2794 + _S2801 + _S2801);
    float3  _S2806 = - - _S2805;
    Matrix<float, 3, 3>  _S2807 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2808;
    (&_S2808)->primal_0 = _S2688;
    (&_S2808)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2809;
    (&_S2809)->primal_0 = t_22;
    (&_S2809)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2808, &_S2809, _S2806);
    Matrix<float, 3, 3>  _S2810 = transpose_0(_S2808.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2811;
    (&_S2811)->primal_0 = _S2687;
    (&_S2811)->differential_0 = _S2727;
    s_bwd_prop_log_1(&_S2811, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2812;
    (&_S2812)->primal_0 = vert2_c_4;
    (&_S2812)->differential_0 = _S2727;
    s_bwd_length_impl_0(&_S2812, _S2811.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2813;
    (&_S2813)->primal_0 = vert1_c_4;
    (&_S2813)->differential_0 = _S2727;
    s_bwd_length_impl_0(&_S2813, _S2811.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2814;
    (&_S2814)->primal_0 = vert0_c_4;
    (&_S2814)->differential_0 = _S2727;
    s_bwd_length_impl_0(&_S2814, _S2811.differential_0.x);
    DiffPair_float_0 _S2815;
    (&_S2815)->primal_0 = _S2686;
    (&_S2815)->differential_0 = 0.0f;
    DiffPair_float_0 _S2816;
    (&_S2816)->primal_0 = _S2685;
    (&_S2816)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2815, &_S2816, 0.0f);
    DiffPair_float_0 _S2817;
    (&_S2817)->primal_0 = _S2682;
    (&_S2817)->differential_0 = 0.0f;
    DiffPair_float_0 _S2818;
    (&_S2818)->primal_0 = _S2683;
    (&_S2818)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2817, &_S2818, _S2815.differential_0);
    DiffPair_float_0 _S2819;
    (&_S2819)->primal_0 = _S2684;
    (&_S2819)->differential_0 = 0.0f;
    DiffPair_float_0 _S2820;
    (&_S2820)->primal_0 = _S2685;
    (&_S2820)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2819, &_S2820, 0.0f);
    float _S2821 = _S2816.differential_0 + _S2820.differential_0;
    DiffPair_float_0 _S2822;
    (&_S2822)->primal_0 = _S2682;
    (&_S2822)->differential_0 = 0.0f;
    DiffPair_float_0 _S2823;
    (&_S2823)->primal_0 = _S2683;
    (&_S2823)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2822, &_S2823, _S2819.differential_0);
    float _S2824 = _S2818.differential_0 + _S2823.differential_0;
    float _S2825 = _S2817.differential_0 + _S2822.differential_0;
    DiffPair_float_0 _S2826;
    (&_S2826)->primal_0 = _S2681;
    (&_S2826)->differential_0 = 0.0f;
    DiffPair_float_0 _S2827;
    (&_S2827)->primal_0 = _S2680;
    (&_S2827)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2826, &_S2827, 0.0f);
    DiffPair_float_0 _S2828;
    (&_S2828)->primal_0 = _S2677;
    (&_S2828)->differential_0 = 0.0f;
    DiffPair_float_0 _S2829;
    (&_S2829)->primal_0 = _S2678;
    (&_S2829)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2828, &_S2829, _S2826.differential_0);
    DiffPair_float_0 _S2830;
    (&_S2830)->primal_0 = _S2679;
    (&_S2830)->differential_0 = 0.0f;
    DiffPair_float_0 _S2831;
    (&_S2831)->primal_0 = _S2680;
    (&_S2831)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2830, &_S2831, 0.0f);
    float _S2832 = _S2827.differential_0 + _S2831.differential_0;
    DiffPair_float_0 _S2833;
    (&_S2833)->primal_0 = _S2677;
    (&_S2833)->differential_0 = 0.0f;
    DiffPair_float_0 _S2834;
    (&_S2834)->primal_0 = _S2678;
    (&_S2834)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2833, &_S2834, _S2830.differential_0);
    float _S2835 = _S2829.differential_0 + _S2834.differential_0;
    float _S2836 = _S2828.differential_0 + _S2833.differential_0;
    DiffPair_float_0 _S2837;
    (&_S2837)->primal_0 = _S2675;
    (&_S2837)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2837, 0.0f);
    float _S2838 = - (-1.0f * - (_S2837.differential_0 / _S2676));
    float2  _S2839 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2840;
    (&_S2840)->primal_0 = e2_0;
    (&_S2840)->differential_0 = _S2839;
    s_bwd_length_impl_1(&_S2840, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2841;
    (&_S2841)->primal_0 = e1_4;
    (&_S2841)->differential_0 = _S2839;
    s_bwd_length_impl_1(&_S2841, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2842;
    (&_S2842)->primal_0 = e0_4;
    (&_S2842)->differential_0 = _S2839;
    s_bwd_length_impl_1(&_S2842, -0.0f);
    DiffPair_float_0 _S2843;
    (&_S2843)->primal_0 = _S2673;
    (&_S2843)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2843, 0.0f);
    float _S2844 = - _S2843.differential_0;
    float2  _S2845 = _S2841.differential_0 + make_float2 (_S2671 * _S2844, _S2669 * _S2843.differential_0);
    float2  _S2846 = _S2842.differential_0 + make_float2 (_S2670 * _S2843.differential_0, _S2672 * _S2844);
    float2  _S2847 = _S2664 * (v_uv2_0 + - _S2840.differential_0 + _S2845 + make_float2 (_S2832, _S2821)) / _S2663;
    float2  _S2848 = _S2660 * - _S2847;
    float2  _S2849 = _S2662 * _S2847;
    float2  _S2850 = _S2664 * (v_uv1_0 + - _S2845 + _S2846 + make_float2 (_S2835, _S2824)) / _S2659;
    float2  _S2851 = _S2656 * - _S2850;
    float2  _S2852 = _S2658 * _S2850;
    float _S2853 = _S2851.x + _S2851.y;
    float2  _S2854 = _S2664 * (v_uv0_0 + _S2840.differential_0 + - _S2846 + make_float2 (_S2836, _S2825)) / _S2655;
    float2  _S2855 = _S2652 * - _S2854;
    float2  _S2856 = _S2654 * _S2854;
    float _S2857 = _S2855.x + _S2855.y;
    float3  _S2858 = _S2733.differential_0 + _S2812.differential_0 + make_float3 (_S2849.x, _S2849.y, _S2848.x + _S2848.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2859;
    (&_S2859)->primal_0 = R_23;
    (&_S2859)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2860;
    (&_S2860)->primal_0 = vert2_0;
    (&_S2860)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2859, &_S2860, _S2858);
    float3  _S2861 = _S2732.differential_0 + _S2813.differential_0 + make_float3 (_S2852.x, _S2852.y, _S2853);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2862;
    (&_S2862)->primal_0 = R_23;
    (&_S2862)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2863;
    (&_S2863)->primal_0 = vert1_0;
    (&_S2863)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2862, &_S2863, _S2861);
    float3  _S2864 = _S2734 + _S2735 + _S2814.differential_0 + make_float3 (_S2856.x, _S2856.y, _S2857);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2865;
    (&_S2865)->primal_0 = R_23;
    (&_S2865)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2866;
    (&_S2866)->primal_0 = vert0_0;
    (&_S2866)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2865, &_S2866, _S2864);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2867;
    (&_S2867)->primal_0 = _S2646;
    (&_S2867)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2868;
    (&_S2868)->primal_0 = _S2651;
    (&_S2868)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2867, &_S2868, _S2860.differential_0);
    float _S2869 = - _S2868.differential_0.y;
    float _S2870 = _S2650 * _S2868.differential_0.x;
    float _S2871 = - (_S2642 * _S2868.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2872;
    (&_S2872)->primal_0 = _S2646;
    (&_S2872)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2873;
    (&_S2873)->primal_0 = _S2649;
    (&_S2873)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2872, &_S2873, _S2863.differential_0);
    float _S2874 = _S2642 * _S2873.differential_0.x;
    float _S2875 = _S2648 * _S2873.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2876;
    (&_S2876)->primal_0 = _S2646;
    (&_S2876)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2877;
    (&_S2877)->primal_0 = _S2647;
    (&_S2877)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2876, &_S2877, _S2866.differential_0);
    Matrix<float, 3, 3>  _S2878 = transpose_0(_S2867.differential_0 + _S2872.differential_0 + _S2876.differential_0);
    float _S2879 = 2.0f * - _S2878.rows[int(2)].z;
    float _S2880 = 2.0f * _S2878.rows[int(2)].y;
    float _S2881 = 2.0f * _S2878.rows[int(2)].x;
    float _S2882 = 2.0f * _S2878.rows[int(1)].z;
    float _S2883 = 2.0f * - _S2878.rows[int(1)].y;
    float _S2884 = 2.0f * _S2878.rows[int(1)].x;
    float _S2885 = 2.0f * _S2878.rows[int(0)].z;
    float _S2886 = 2.0f * _S2878.rows[int(0)].y;
    float _S2887 = 2.0f * - _S2878.rows[int(0)].x;
    float _S2888 = - _S2884 + _S2886;
    float _S2889 = _S2881 + - _S2885;
    float _S2890 = - _S2880 + _S2882;
    float _S2891 = _S2880 + _S2882;
    float _S2892 = _S2881 + _S2885;
    float _S2893 = _S2884 + _S2886;
    float _S2894 = quat_24.w * (_S2883 + _S2887);
    float _S2895 = quat_24.z * (_S2879 + _S2887);
    float _S2896 = quat_24.y * (_S2879 + _S2883);
    float _S2897 = quat_24.x * _S2888 + quat_24.z * _S2891 + quat_24.y * _S2892 + _S2894 + _S2894;
    float _S2898 = quat_24.x * _S2889 + quat_24.w * _S2891 + quat_24.y * _S2893 + _S2895 + _S2895;
    float _S2899 = quat_24.x * _S2890 + quat_24.w * _S2892 + quat_24.z * _S2893 + _S2896 + _S2896;
    float _S2900 = quat_24.w * _S2888 + quat_24.z * _S2889 + quat_24.y * _S2890;
    float _S2901 = _S2871 + _S2874;
    float _S2902 = 0.5f * - _S2901;
    float _S2903 = _S2869 + _S2873.differential_0.y;
    DiffPair_float_0 _S2904;
    (&_S2904)->primal_0 = _S2643;
    (&_S2904)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2904, _S2903);
    float _S2905 = _S2902 + _S2904.differential_0;
    float _S2906 = _S2870 + _S2875 + _S2877.differential_0.x;
    DiffPair_float_0 _S2907;
    (&_S2907)->primal_0 = _S2641;
    (&_S2907)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2907, _S2906);
    float _S2908 = _S2902 + _S2907.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2909;
    (&_S2909)->primal_0 = R_23;
    (&_S2909)->differential_0 = _S2807;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2910;
    (&_S2910)->primal_0 = mean_23;
    (&_S2910)->differential_0 = _S2727;
    s_bwd_prop_mul_1(&_S2909, &_S2910, _S2729.differential_0);
    float3  _S2911 = _S2809.differential_0 + _S2858 + _S2861 + _S2864 + _S2729.differential_0;
    Matrix<float, 3, 3>  _S2912 = _S2810 + _S2859.differential_0 + _S2862.differential_0 + _S2865.differential_0 + _S2909.differential_0;
    FixedArray<float3 , 2>  _S2913;
    _S2913[int(0)] = _S2727;
    _S2913[int(1)] = _S2727;
    _S2913[int(1)] = _S2741;
    _S2913[int(0)] = _S2746;
    FixedArray<float3 , 16>  _S2914;
    _S2914[int(0)] = _S2727;
    _S2914[int(1)] = _S2727;
    _S2914[int(2)] = _S2727;
    _S2914[int(3)] = _S2727;
    _S2914[int(4)] = _S2727;
    _S2914[int(5)] = _S2727;
    _S2914[int(6)] = _S2727;
    _S2914[int(7)] = _S2727;
    _S2914[int(8)] = _S2727;
    _S2914[int(9)] = _S2727;
    _S2914[int(10)] = _S2727;
    _S2914[int(11)] = _S2727;
    _S2914[int(12)] = _S2727;
    _S2914[int(13)] = _S2727;
    _S2914[int(14)] = _S2727;
    _S2914[int(15)] = _S2727;
    _S2914[int(15)] = _S2748;
    _S2914[int(14)] = _S2750;
    _S2914[int(13)] = _S2752;
    _S2914[int(12)] = _S2754;
    _S2914[int(11)] = _S2756;
    _S2914[int(10)] = _S2758;
    _S2914[int(9)] = _S2760;
    _S2914[int(8)] = _S2768;
    _S2914[int(7)] = _S2770;
    _S2914[int(6)] = _S2772;
    _S2914[int(5)] = _S2774;
    _S2914[int(4)] = _S2776;
    _S2914[int(3)] = _S2787;
    _S2914[int(2)] = _S2789;
    _S2914[int(1)] = _S2791;
    _S2914[int(0)] = _S2804;
    float2  _S2915 = v_out_hardness_0 + make_float2 (0.0f, _S2838);
    float3  _S2916 = make_float3 (_S2908, _S2905, _S2901);
    float4  _S2917 = make_float4 (0.0f);
    *&((&_S2917)->w) = _S2897;
    *&((&_S2917)->z) = _S2898;
    *&((&_S2917)->y) = _S2899;
    *&((&_S2917)->x) = _S2900;
    *v_mean_7 = _S2805 + _S2860.differential_0 + _S2863.differential_0 + _S2866.differential_0 + _S2910.differential_0;
    *v_quat_6 = _S2917;
    *v_scale_6 = _S2916;
    *v_hardness_0 = _S2915;
    *v_sh_coeffs_5 = _S2914;
    *v_ch_coeffs_0 = _S2913;
    *v_R_6 = _S2912;
    *v_t_6 = _S2911;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_24, float3  t_23, float fx_24, float fy_24, float cx_24, float cy_24, float4  radial_coeffs_27, float2  tangential_coeffs_27, float2  thin_prism_coeffs_27, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_1(R_24, mean_24) + t_23;
    float _S2918 = scale_24.x;
    float _S2919 = s_primal_ctx_exp_1(_S2918);
    float _S2920 = scale_24.y;
    float _S2921 = s_primal_ctx_exp_1(_S2920);
    float sz_5 = scale_24.z - 0.5f * (_S2918 + _S2920);
    float _S2922 = quat_25.y;
    float x2_25 = _S2922 * _S2922;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_45 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S2923 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_45), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_45), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S2924 = make_float3 (_S2919, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S2923, _S2924) + mean_24;
    float _S2925 = -0.5f + sz_5;
    float3  _S2926 = make_float3 (_S2919 * _S2925, _S2921, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S2923, _S2926) + mean_24;
    float _S2927 = -0.5f - sz_5;
    float3  _S2928 = make_float3 (_S2919 * _S2927, - _S2921, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S2923, _S2928) + mean_24;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_24, vert0_1) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_24, vert1_1) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_24, vert2_1) + t_23;
    float _S2929 = length_1(vert0_c_5);
    float _S2930 = length_1(vert1_c_5);
    float _S2931 = length_1(vert2_c_5);
    CameraDistortion_0 _S2932 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_27, tangential_coeffs_27, thin_prism_coeffs_27);
    float2  _S2933 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S2934 = length_0(_S2933);
    float _S2935 = vert0_c_5.z;
    float _S2936 = s_primal_ctx_atan2_0(_S2934, _S2935);
    bool _S2937 = _S2936 < 0.00100000004749745f;
    float k_6;
    float _S2938;
    float _S2939;
    float _S2940;
    if(_S2937)
    {
        float _S2941 = 1.0f - _S2936 * _S2936 / 3.0f;
        float _S2942 = _S2935 * _S2935;
        k_6 = _S2941 / _S2935;
        _S2938 = 0.0f;
        _S2939 = _S2942;
        _S2940 = _S2941;
    }
    else
    {
        float _S2943 = _S2934 * _S2934;
        k_6 = _S2936 / _S2934;
        _S2938 = _S2943;
        _S2939 = 0.0f;
        _S2940 = 0.0f;
    }
    float2  _S2944 = make_float2 (k_6);
    float2  _S2945 = _S2933 * make_float2 (k_6);
    float k1_5 = _S2932.radial_coeffs_0.x;
    float k2_5 = _S2932.radial_coeffs_0.y;
    float k3_5 = _S2932.radial_coeffs_0.z;
    float k4_6 = _S2932.radial_coeffs_0.w;
    float p1_6 = _S2932.tangential_coeffs_0.x;
    float p2_6 = _S2932.tangential_coeffs_0.y;
    float sx1_6 = _S2932.thin_prism_coeffs_0.x;
    float sy1_6 = _S2932.thin_prism_coeffs_0.y;
    float u_16 = _S2945.x;
    float v_16 = _S2945.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S2946 = k3_5 + r2_16 * k4_6;
    float _S2947 = k2_5 + r2_16 * _S2946;
    float _S2948 = k1_5 + r2_16 * _S2947;
    float radial_2 = 1.0f + r2_16 * _S2948;
    float2  _S2949 = make_float2 (radial_2);
    float _S2950 = 2.0f * p1_6;
    float _S2951 = _S2950 * u_16;
    float _S2952 = 2.0f * u_16;
    float _S2953 = r2_16 + _S2952 * u_16;
    float _S2954 = 2.0f * p2_6;
    float _S2955 = _S2954 * u_16;
    float _S2956 = 2.0f * v_16;
    float _S2957 = r2_16 + _S2956 * v_16;
    float2  _S2958 = _S2945 * make_float2 (radial_2) + make_float2 (_S2951 * v_16 + p2_6 * _S2953 + sx1_6 * r2_16, _S2955 * v_16 + p1_6 * _S2957 + sy1_6 * r2_16);
    float _S2959 = fx_24 * _S2958.x + cx_24;
    float _S2960 = fy_24 * _S2958.y + cy_24;
    float2  _S2961 = make_float2 (_S2959, _S2960);
    float2  _S2962 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S2963 = length_0(_S2962);
    float _S2964 = vert1_c_5.z;
    float _S2965 = s_primal_ctx_atan2_0(_S2963, _S2964);
    bool _S2966 = _S2965 < 0.00100000004749745f;
    float _S2967;
    float _S2968;
    float _S2969;
    if(_S2966)
    {
        float _S2970 = 1.0f - _S2965 * _S2965 / 3.0f;
        float _S2971 = _S2964 * _S2964;
        k_6 = _S2970 / _S2964;
        _S2967 = 0.0f;
        _S2968 = _S2971;
        _S2969 = _S2970;
    }
    else
    {
        float _S2972 = _S2963 * _S2963;
        k_6 = _S2965 / _S2963;
        _S2967 = _S2972;
        _S2968 = 0.0f;
        _S2969 = 0.0f;
    }
    float2  _S2973 = make_float2 (k_6);
    float2  _S2974 = _S2962 * make_float2 (k_6);
    float u_17 = _S2974.x;
    float v_17 = _S2974.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S2975 = k3_5 + r2_17 * k4_6;
    float _S2976 = k2_5 + r2_17 * _S2975;
    float _S2977 = k1_5 + r2_17 * _S2976;
    float radial_3 = 1.0f + r2_17 * _S2977;
    float2  _S2978 = make_float2 (radial_3);
    float _S2979 = _S2950 * u_17;
    float _S2980 = 2.0f * u_17;
    float _S2981 = r2_17 + _S2980 * u_17;
    float _S2982 = _S2954 * u_17;
    float _S2983 = 2.0f * v_17;
    float _S2984 = r2_17 + _S2983 * v_17;
    float2  _S2985 = _S2974 * make_float2 (radial_3) + make_float2 (_S2979 * v_17 + p2_6 * _S2981 + sx1_6 * r2_17, _S2982 * v_17 + p1_6 * _S2984 + sy1_6 * r2_17);
    float _S2986 = fx_24 * _S2985.x + cx_24;
    float _S2987 = fy_24 * _S2985.y + cy_24;
    float2  _S2988 = make_float2 (_S2986, _S2987);
    float2  _S2989 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S2990 = length_0(_S2989);
    float _S2991 = vert2_c_5.z;
    float _S2992 = s_primal_ctx_atan2_0(_S2990, _S2991);
    bool _S2993 = _S2992 < 0.00100000004749745f;
    float _S2994;
    float _S2995;
    float _S2996;
    if(_S2993)
    {
        float _S2997 = 1.0f - _S2992 * _S2992 / 3.0f;
        float _S2998 = _S2991 * _S2991;
        k_6 = _S2997 / _S2991;
        _S2994 = 0.0f;
        _S2995 = _S2998;
        _S2996 = _S2997;
    }
    else
    {
        float _S2999 = _S2990 * _S2990;
        k_6 = _S2992 / _S2990;
        _S2994 = _S2999;
        _S2995 = 0.0f;
        _S2996 = 0.0f;
    }
    float2  _S3000 = make_float2 (k_6);
    float2  _S3001 = _S2989 * make_float2 (k_6);
    float u_18 = _S3001.x;
    float v_18 = _S3001.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S3002 = k3_5 + r2_18 * k4_6;
    float _S3003 = k2_5 + r2_18 * _S3002;
    float _S3004 = k1_5 + r2_18 * _S3003;
    float radial_4 = 1.0f + r2_18 * _S3004;
    float2  _S3005 = make_float2 (radial_4);
    float _S3006 = _S2950 * u_18;
    float _S3007 = 2.0f * u_18;
    float _S3008 = r2_18 + _S3007 * u_18;
    float _S3009 = _S2954 * u_18;
    float _S3010 = 2.0f * v_18;
    float _S3011 = r2_18 + _S3010 * v_18;
    float2  _S3012 = _S3001 * make_float2 (radial_4) + make_float2 (_S3006 * v_18 + p2_6 * _S3008 + sx1_6 * r2_18, _S3009 * v_18 + p1_6 * _S3011 + sy1_6 * r2_18);
    float _S3013 = fx_24 * _S3012.x + cx_24;
    float _S3014 = fy_24 * _S3012.y + cy_24;
    float2  _S3015 = make_float2 (_S3013, _S3014);
    float2  e0_5 = _S2988 - _S2961;
    float2  e1_5 = _S3015 - _S2988;
    float2  e2_1 = _S2961 - _S3015;
    float _S3016 = e0_5.x;
    float _S3017 = e1_5.y;
    float _S3018 = e0_5.y;
    float _S3019 = e1_5.x;
    float _S3020 = _S3016 * _S3017 - _S3018 * _S3019;
    float _S3021 = 1.0f - hardness_5.y;
    float _S3022 = -1.0f / _S3021;
    float _S3023 = _S3021 * _S3021;
    float _S3024 = s_primal_ctx_max_0(_S2959, _S2986);
    float _S3025 = s_primal_ctx_min_0(_S2959, _S2986);
    float _S3026 = s_primal_ctx_max_0(_S2960, _S2987);
    float _S3027 = s_primal_ctx_min_0(_S2960, _S2987);
    float3  _S3028 = make_float3 (_S2929, _S2930, _S2931) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3029 = transpose_0(R_24);
    float3  _S3030 = mean_24 - - s_primal_ctx_mul_1(_S3029, t_23);
    float _S3031 = _S3030.x;
    float _S3032 = _S3030.y;
    float _S3033 = _S3030.z;
    float _S3034 = _S3031 * _S3031 + _S3032 * _S3032 + _S3033 * _S3033;
    float _S3035 = s_primal_ctx_sqrt_0(_S3034);
    float x_50 = _S3031 / _S3035;
    float3  _S3036 = make_float3 (x_50);
    float _S3037 = _S3035 * _S3035;
    float y_23 = _S3032 / _S3035;
    float z_20 = _S3033 / _S3035;
    float3  _S3038 = make_float3 (z_20);
    float _S3039 = - y_23;
    float3  _S3040 = make_float3 (_S3039);
    float z2_46 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_50 * x_50 - y_23 * y_23;
    float _S3041 = 2.0f * x_50;
    float fS1_20 = _S3041 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_46 - 0.31539157032966614f;
    float3  _S3042 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_50;
    float3  _S3043 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3044 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3045 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3046 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_46 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3047 = 1.86588168144226074f * z2_46 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3047;
    float3  _S3048 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_50;
    float3  _S3049 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3050 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3051 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3052 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_50 * fC1_20 - y_23 * fS1_20);
    float3  _S3053 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_50 * fS1_20 + y_23 * fC1_20);
    float3  _S3054 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3039) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_50) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3055 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3056 = make_float3 (0.0f);
    float3  _S3057 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3058 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3059 = make_float3 (_S3058);
    float3  _S3060 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3058);
    float3  _S3061 = _S3057 + _S3060 + make_float3 (0.5f);
    float3  _S3062 = _S3057 - _S3060 + make_float3 (0.5f);
    float3  _S3063 = vert1_c_5 - vert0_c_5;
    float3  _S3064 = vert2_c_5 - vert0_c_5;
    float3  _S3065 = s_primal_ctx_cross_0(_S3063, _S3064);
    float3  _S3066 = normalize_0(_S3065);
    float3  _S3067 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3066, mean_c_20)))))) * v_normal_1;
    float3  _S3068 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3069;
    (&_S3069)->primal_0 = _S3066;
    (&_S3069)->differential_0 = _S3068;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3070;
    (&_S3070)->primal_0 = mean_c_20;
    (&_S3070)->differential_0 = _S3068;
    s_bwd_prop_dot_0(&_S3069, &_S3070, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3071 = _S3070;
    float3  _S3072 = _S3067 + _S3069.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3073;
    (&_S3073)->primal_0 = _S3065;
    (&_S3073)->differential_0 = _S3068;
    s_bwd_normalize_impl_0(&_S3073, _S3072);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3074;
    (&_S3074)->primal_0 = _S3063;
    (&_S3074)->differential_0 = _S3068;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3075;
    (&_S3075)->primal_0 = _S3064;
    (&_S3075)->differential_0 = _S3068;
    s_bwd_prop_cross_0(&_S3074, &_S3075, _S3073.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3076 = _S3074;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3077 = _S3075;
    float3  _S3078 = - _S3075.differential_0;
    float3  _S3079 = - _S3074.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3080;
    (&_S3080)->primal_0 = _S3062;
    (&_S3080)->differential_0 = _S3068;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3081;
    (&_S3081)->primal_0 = _S3056;
    (&_S3081)->differential_0 = _S3068;
    s_bwd_prop_max_0(&_S3080, &_S3081, (*v_rgb_7)[int(2)]);
    float3  _S3082 = - _S3080.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3083;
    (&_S3083)->primal_0 = _S3061;
    (&_S3083)->differential_0 = _S3068;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3084;
    (&_S3084)->primal_0 = _S3056;
    (&_S3084)->differential_0 = _S3068;
    s_bwd_prop_max_0(&_S3083, &_S3084, (*v_rgb_7)[int(1)]);
    float3  _S3085 = _S3059 * (_S3082 + _S3083.differential_0);
    float3  _S3086 = _S3080.differential_0 + _S3083.differential_0;
    float3  _S3087 = make_float3 (0.5f) * - _S3086;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3088;
    (&_S3088)->primal_0 = _S3055;
    (&_S3088)->differential_0 = _S3068;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3089;
    (&_S3089)->primal_0 = _S3056;
    (&_S3089)->differential_0 = _S3068;
    s_bwd_prop_max_0(&_S3088, &_S3089, (*v_rgb_7)[int(0)]);
    float3  _S3090 = _S3087 + _S3088.differential_0;
    float3  _S3091 = _S3086 + _S3088.differential_0;
    float3  _S3092 = _S3053 * _S3091;
    float3  _S3093 = (*sh_coeffs_20)[int(15)] * _S3091;
    float3  _S3094 = _S3051 * _S3091;
    float3  _S3095 = (*sh_coeffs_20)[int(14)] * _S3091;
    float3  _S3096 = _S3049 * _S3091;
    float3  _S3097 = (*sh_coeffs_20)[int(13)] * _S3091;
    float3  _S3098 = _S3048 * _S3091;
    float3  _S3099 = (*sh_coeffs_20)[int(12)] * _S3091;
    float3  _S3100 = _S3050 * _S3091;
    float3  _S3101 = (*sh_coeffs_20)[int(11)] * _S3091;
    float3  _S3102 = _S3052 * _S3091;
    float3  _S3103 = (*sh_coeffs_20)[int(10)] * _S3091;
    float3  _S3104 = _S3054 * _S3091;
    float3  _S3105 = (*sh_coeffs_20)[int(9)] * _S3091;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3105.x + _S3105.y + _S3105.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3093.x + _S3093.y + _S3093.z);
    float _S3106 = _S3103.x + _S3103.y + _S3103.z;
    float _S3107 = _S3095.x + _S3095.y + _S3095.z;
    float _S3108 = _S3101.x + _S3101.y + _S3101.z;
    float _S3109 = _S3097.x + _S3097.y + _S3097.z;
    float _S3110 = _S3099.x + _S3099.y + _S3099.z;
    float _S3111 = - s_diff_fC2_T_6;
    float3  _S3112 = _S3045 * _S3091;
    float3  _S3113 = (*sh_coeffs_20)[int(8)] * _S3091;
    float3  _S3114 = _S3043 * _S3091;
    float3  _S3115 = (*sh_coeffs_20)[int(7)] * _S3091;
    float3  _S3116 = _S3042 * _S3091;
    float3  _S3117 = (*sh_coeffs_20)[int(6)] * _S3091;
    float3  _S3118 = _S3044 * _S3091;
    float3  _S3119 = (*sh_coeffs_20)[int(5)] * _S3091;
    float3  _S3120 = _S3046 * _S3091;
    float3  _S3121 = (*sh_coeffs_20)[int(4)] * _S3091;
    float _S3122 = _S3119.x + _S3119.y + _S3119.z;
    float _S3123 = _S3115.x + _S3115.y + _S3115.z;
    float _S3124 = fTmp1B_20 * _S3106 + x_50 * s_diff_fS2_T_6 + y_23 * _S3111 + 0.54627424478530884f * (_S3121.x + _S3121.y + _S3121.z);
    float _S3125 = fTmp1B_20 * _S3107 + y_23 * s_diff_fS2_T_6 + x_50 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3113.x + _S3113.y + _S3113.z);
    float _S3126 = y_23 * - _S3125;
    float _S3127 = x_50 * _S3125;
    float _S3128 = z_20 * (1.86588168144226074f * (z_20 * _S3110) + -2.28522896766662598f * (y_23 * _S3108 + x_50 * _S3109) + 0.94617468118667603f * (_S3117.x + _S3117.y + _S3117.z));
    float3  _S3129 = make_float3 (0.48860251903533936f) * _S3091;
    float3  _S3130 = - _S3129;
    float3  _S3131 = _S3036 * _S3130;
    float3  _S3132 = (*sh_coeffs_20)[int(3)] * _S3130;
    float3  _S3133 = _S3038 * _S3129;
    float3  _S3134 = (*sh_coeffs_20)[int(2)] * _S3129;
    float3  _S3135 = _S3040 * _S3129;
    float3  _S3136 = (*sh_coeffs_20)[int(1)] * _S3129;
    float _S3137 = (_S3047 * _S3110 + 1.44530570507049561f * (fS1_20 * _S3106 + fC1_20 * _S3107) + -1.09254848957061768f * (y_23 * _S3122 + x_50 * _S3123) + _S3128 + _S3128 + _S3134.x + _S3134.y + _S3134.z) / _S3037;
    float _S3138 = _S3035 * _S3137;
    float _S3139 = (fTmp0C_20 * _S3108 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3111 + fTmp0B_20 * _S3122 + _S3041 * _S3124 + _S3126 + _S3126 + - (_S3136.x + _S3136.y + _S3136.z)) / _S3037;
    float _S3140 = _S3035 * _S3139;
    float _S3141 = (fTmp0C_20 * _S3109 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3123 + 2.0f * (y_23 * _S3124) + _S3127 + _S3127 + _S3132.x + _S3132.y + _S3132.z) / _S3037;
    float _S3142 = _S3035 * _S3141;
    float _S3143 = _S3033 * - _S3137 + _S3032 * - _S3139 + _S3031 * - _S3141;
    DiffPair_float_0 _S3144;
    (&_S3144)->primal_0 = _S3034;
    (&_S3144)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3144, _S3143);
    float _S3145 = _S3033 * _S3144.differential_0;
    float _S3146 = _S3032 * _S3144.differential_0;
    float _S3147 = _S3031 * _S3144.differential_0;
    float3  _S3148 = make_float3 (0.282094806432724f) * _S3091;
    float3  _S3149 = make_float3 (_S3142 + _S3147 + _S3147, _S3140 + _S3146 + _S3146, _S3138 + _S3145 + _S3145);
    float3  _S3150 = - - _S3149;
    Matrix<float, 3, 3>  _S3151 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3152;
    (&_S3152)->primal_0 = _S3029;
    (&_S3152)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3153;
    (&_S3153)->primal_0 = t_23;
    (&_S3153)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3152, &_S3153, _S3150);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3154 = _S3153;
    Matrix<float, 3, 3>  _S3155 = transpose_0(_S3152.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3156;
    (&_S3156)->primal_0 = _S3028;
    (&_S3156)->differential_0 = _S3068;
    s_bwd_prop_log_1(&_S3156, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3157 = _S3156;
    DiffPair_float_0 _S3158;
    (&_S3158)->primal_0 = _S3027;
    (&_S3158)->differential_0 = 0.0f;
    DiffPair_float_0 _S3159;
    (&_S3159)->primal_0 = _S3014;
    (&_S3159)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3158, &_S3159, 0.0f);
    DiffPair_float_0 _S3160;
    (&_S3160)->primal_0 = _S2960;
    (&_S3160)->differential_0 = 0.0f;
    DiffPair_float_0 _S3161;
    (&_S3161)->primal_0 = _S2987;
    (&_S3161)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3160, &_S3161, _S3158.differential_0);
    DiffPair_float_0 _S3162;
    (&_S3162)->primal_0 = _S3026;
    (&_S3162)->differential_0 = 0.0f;
    DiffPair_float_0 _S3163;
    (&_S3163)->primal_0 = _S3014;
    (&_S3163)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3162, &_S3163, 0.0f);
    DiffPair_float_0 _S3164;
    (&_S3164)->primal_0 = _S2960;
    (&_S3164)->differential_0 = 0.0f;
    DiffPair_float_0 _S3165;
    (&_S3165)->primal_0 = _S2987;
    (&_S3165)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3164, &_S3165, _S3162.differential_0);
    DiffPair_float_0 _S3166;
    (&_S3166)->primal_0 = _S3025;
    (&_S3166)->differential_0 = 0.0f;
    DiffPair_float_0 _S3167;
    (&_S3167)->primal_0 = _S3013;
    (&_S3167)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3166, &_S3167, 0.0f);
    DiffPair_float_0 _S3168;
    (&_S3168)->primal_0 = _S2959;
    (&_S3168)->differential_0 = 0.0f;
    DiffPair_float_0 _S3169;
    (&_S3169)->primal_0 = _S2986;
    (&_S3169)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3168, &_S3169, _S3166.differential_0);
    DiffPair_float_0 _S3170;
    (&_S3170)->primal_0 = _S3024;
    (&_S3170)->differential_0 = 0.0f;
    DiffPair_float_0 _S3171;
    (&_S3171)->primal_0 = _S3013;
    (&_S3171)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3170, &_S3171, 0.0f);
    DiffPair_float_0 _S3172;
    (&_S3172)->primal_0 = _S2959;
    (&_S3172)->differential_0 = 0.0f;
    DiffPair_float_0 _S3173;
    (&_S3173)->primal_0 = _S2986;
    (&_S3173)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3172, &_S3173, _S3170.differential_0);
    DiffPair_float_0 _S3174;
    (&_S3174)->primal_0 = _S3022;
    (&_S3174)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3174, 0.0f);
    float _S3175 = - (-1.0f * - (_S3174.differential_0 / _S3023));
    float2  _S3176 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3177;
    (&_S3177)->primal_0 = e2_1;
    (&_S3177)->differential_0 = _S3176;
    s_bwd_length_impl_1(&_S3177, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3178;
    (&_S3178)->primal_0 = e1_5;
    (&_S3178)->differential_0 = _S3176;
    s_bwd_length_impl_1(&_S3178, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3179;
    (&_S3179)->primal_0 = e0_5;
    (&_S3179)->differential_0 = _S3176;
    s_bwd_length_impl_1(&_S3179, -0.0f);
    DiffPair_float_0 _S3180;
    (&_S3180)->primal_0 = _S3020;
    (&_S3180)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3180, 0.0f);
    float _S3181 = - _S3180.differential_0;
    float2  _S3182 = _S3178.differential_0 + make_float2 (_S3018 * _S3181, _S3016 * _S3180.differential_0);
    float2  _S3183 = _S3179.differential_0 + make_float2 (_S3017 * _S3180.differential_0, _S3019 * _S3181);
    float2  _S3184 = v_uv2_1 + - _S3177.differential_0 + _S3182;
    float _S3185 = fy_24 * (_S3159.differential_0 + _S3163.differential_0 + _S3184.y);
    float _S3186 = fx_24 * (_S3167.differential_0 + _S3171.differential_0 + _S3184.x);
    float2  _S3187 = make_float2 (_S3186, _S3185);
    float2  _S3188 = _S3001 * _S3187;
    float2  _S3189 = _S3005 * _S3187;
    float _S3190 = r2_18 * _S3185;
    float _S3191 = p1_6 * _S3185;
    float _S3192 = _S3011 * _S3185;
    float _S3193 = v_18 * _S3185;
    float _S3194 = u_18 * _S3193;
    float _S3195 = r2_18 * _S3186;
    float _S3196 = p2_6 * _S3186;
    float _S3197 = _S3008 * _S3186;
    float _S3198 = v_18 * _S3186;
    float _S3199 = u_18 * _S3198;
    float _S3200 = _S3188.x + _S3188.y;
    float _S3201 = r2_18 * _S3200;
    float _S3202 = r2_18 * _S3201;
    float _S3203 = r2_18 * _S3202;
    float _S3204 = r2_18 * _S3203;
    float _S3205 = sy1_6 * _S3185 + _S3191 + sx1_6 * _S3186 + _S3196 + _S3004 * _S3200 + _S3003 * _S3201 + _S3002 * _S3202 + k4_6 * _S3203;
    float _S3206 = v_18 * _S3205;
    float _S3207 = u_18 * _S3205;
    float _S3208 = _S3010 * _S3191 + 2.0f * (v_18 * _S3191) + _S3009 * _S3185 + _S3006 * _S3186 + _S3206 + _S3206;
    float _S3209 = _S2954 * _S3193 + _S3007 * _S3196 + 2.0f * (u_18 * _S3196) + _S2950 * _S3198 + _S3207 + _S3207;
    float2  _S3210 = v_uv0_1 + _S3177.differential_0 + - _S3183;
    float2  _S3211 = v_out_hardness_1 + make_float2 (0.0f, _S3175);
    float _S3212 = _S3169.differential_0 + _S3173.differential_0;
    float2  _S3213 = v_uv1_1 + - _S3182 + _S3183;
    float3  _S3214 = _S3078 + _S3079;
    FixedArray<float3 , 2>  _S3215;
    _S3215[int(0)] = _S3068;
    _S3215[int(1)] = _S3068;
    _S3215[int(1)] = _S3085;
    _S3215[int(0)] = _S3090;
    float3  _S3216 = _S3215[int(0)];
    float3  _S3217 = _S3215[int(1)];
    FixedArray<float3 , 16>  _S3218;
    _S3218[int(0)] = _S3068;
    _S3218[int(1)] = _S3068;
    _S3218[int(2)] = _S3068;
    _S3218[int(3)] = _S3068;
    _S3218[int(4)] = _S3068;
    _S3218[int(5)] = _S3068;
    _S3218[int(6)] = _S3068;
    _S3218[int(7)] = _S3068;
    _S3218[int(8)] = _S3068;
    _S3218[int(9)] = _S3068;
    _S3218[int(10)] = _S3068;
    _S3218[int(11)] = _S3068;
    _S3218[int(12)] = _S3068;
    _S3218[int(13)] = _S3068;
    _S3218[int(14)] = _S3068;
    _S3218[int(15)] = _S3068;
    _S3218[int(7)] = _S3114;
    _S3218[int(0)] = _S3148;
    _S3218[int(1)] = _S3135;
    _S3218[int(2)] = _S3133;
    _S3218[int(3)] = _S3131;
    _S3218[int(4)] = _S3120;
    _S3218[int(5)] = _S3118;
    _S3218[int(6)] = _S3116;
    _S3218[int(15)] = _S3092;
    _S3218[int(8)] = _S3112;
    _S3218[int(9)] = _S3104;
    _S3218[int(10)] = _S3102;
    _S3218[int(11)] = _S3100;
    _S3218[int(12)] = _S3098;
    _S3218[int(13)] = _S3096;
    _S3218[int(14)] = _S3094;
    float3  _S3219 = _S3218[int(0)];
    float3  _S3220 = _S3218[int(1)];
    float3  _S3221 = _S3218[int(2)];
    float3  _S3222 = _S3218[int(3)];
    float3  _S3223 = _S3218[int(4)];
    float3  _S3224 = _S3218[int(5)];
    float3  _S3225 = _S3218[int(6)];
    float3  _S3226 = _S3218[int(7)];
    float3  _S3227 = _S3218[int(8)];
    float3  _S3228 = _S3218[int(9)];
    float3  _S3229 = _S3218[int(10)];
    float3  _S3230 = _S3218[int(11)];
    float3  _S3231 = _S3218[int(12)];
    float3  _S3232 = _S3218[int(13)];
    float3  _S3233 = _S3218[int(14)];
    float3  _S3234 = _S3218[int(15)];
    float _S3235 = _S3168.differential_0 + _S3172.differential_0;
    float _S3236 = _S3160.differential_0 + _S3164.differential_0;
    float _S3237 = _S3161.differential_0 + _S3165.differential_0;
    float2  _S3238 = _S3189 + make_float2 (_S3209, _S3208);
    float2  _S3239 = _S2989 * _S3238;
    float2  _S3240 = _S3000 * _S3238;
    float _S3241 = _S3239.x + _S3239.y;
    if(_S2993)
    {
        float _S3242 = _S3241 / _S2995;
        float _S3243 = _S2996 * - _S3242;
        float _S3244 = _S2992 * (0.3333333432674408f * - (_S2991 * _S3242));
        k_6 = _S3244 + _S3244;
        _S2994 = _S3243;
        _S2995 = 0.0f;
    }
    else
    {
        float _S3245 = _S3241 / _S2994;
        float _S3246 = _S2992 * - _S3245;
        k_6 = _S2990 * _S3245;
        _S2994 = 0.0f;
        _S2995 = _S3246;
    }
    DiffPair_float_0 _S3247;
    (&_S3247)->primal_0 = _S2990;
    (&_S3247)->differential_0 = 0.0f;
    DiffPair_float_0 _S3248;
    (&_S3248)->primal_0 = _S2991;
    (&_S3248)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3247, &_S3248, k_6);
    float _S3249 = _S3248.differential_0 + _S2994;
    float _S3250 = _S3247.differential_0 + _S2995;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3251;
    (&_S3251)->primal_0 = _S2989;
    (&_S3251)->differential_0 = _S3176;
    s_bwd_length_impl_1(&_S3251, _S3250);
    float2  _S3252 = _S3251.differential_0 + _S3240;
    float _S3253 = fy_24 * (_S3213.y + _S3237);
    float _S3254 = fx_24 * (_S3213.x + _S3212);
    float2  _S3255 = make_float2 (_S3254, _S3253);
    float2  _S3256 = _S2974 * _S3255;
    float _S3257 = p1_6 * _S3253;
    float _S3258 = v_17 * _S3253;
    float _S3259 = p2_6 * _S3254;
    float _S3260 = v_17 * _S3254;
    float _S3261 = _S3256.x + _S3256.y;
    float _S3262 = r2_17 * _S3261;
    float _S3263 = r2_17 * _S3262;
    float _S3264 = r2_17 * _S3263;
    float _S3265 = sy1_6 * _S3253 + _S3257 + sx1_6 * _S3254 + _S3259 + _S2977 * _S3261 + _S2976 * _S3262 + _S2975 * _S3263 + k4_6 * _S3264;
    float _S3266 = v_17 * _S3265;
    float _S3267 = u_17 * _S3265;
    float3  _S3268 = _S3077.differential_0 + make_float3 (_S3252.x, _S3252.y, _S3249);
    float2  _S3269 = _S2978 * _S3255 + make_float2 (_S2954 * _S3258 + _S2980 * _S3259 + 2.0f * (u_17 * _S3259) + _S2950 * _S3260 + _S3267 + _S3267, _S2983 * _S3257 + 2.0f * (v_17 * _S3257) + _S2982 * _S3253 + _S2979 * _S3254 + _S3266 + _S3266);
    float _S3270 = u_17 * _S3258 + _S3194;
    float _S3271 = u_17 * _S3260 + _S3199;
    float _S3272 = r2_17 * _S3253 + _S3190;
    float _S3273 = r2_17 * _S3264 + _S3204;
    float _S3274 = _S3264 + _S3203;
    float _S3275 = _S2984 * _S3253 + _S3192;
    float _S3276 = _S3263 + _S3202;
    float _S3277 = r2_17 * _S3254 + _S3195;
    float _S3278 = _S2981 * _S3254 + _S3197;
    float _S3279 = _S3262 + _S3201;
    float2  _S3280 = _S2962 * _S3269;
    float2  _S3281 = _S2973 * _S3269;
    float _S3282 = _S3280.x + _S3280.y;
    if(_S2966)
    {
        float _S3283 = _S3282 / _S2968;
        float _S3284 = _S2969 * - _S3283;
        float _S3285 = _S2965 * (0.3333333432674408f * - (_S2964 * _S3283));
        k_6 = _S3285 + _S3285;
        _S2967 = _S3284;
        _S2968 = 0.0f;
    }
    else
    {
        float _S3286 = _S3282 / _S2967;
        float _S3287 = _S2965 * - _S3286;
        k_6 = _S2963 * _S3286;
        _S2967 = 0.0f;
        _S2968 = _S3287;
    }
    DiffPair_float_0 _S3288;
    (&_S3288)->primal_0 = _S2963;
    (&_S3288)->differential_0 = 0.0f;
    DiffPair_float_0 _S3289;
    (&_S3289)->primal_0 = _S2964;
    (&_S3289)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3288, &_S3289, k_6);
    float _S3290 = _S3289.differential_0 + _S2967;
    float _S3291 = _S3288.differential_0 + _S2968;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3292;
    (&_S3292)->primal_0 = _S2962;
    (&_S3292)->differential_0 = _S3176;
    s_bwd_length_impl_1(&_S3292, _S3291);
    float2  _S3293 = _S3292.differential_0 + _S3281;
    float _S3294 = fy_24 * (_S3210.y + _S3236);
    float _S3295 = fx_24 * (_S3210.x + _S3235);
    float2  _S3296 = make_float2 (_S3295, _S3294);
    float2  _S3297 = _S2945 * _S3296;
    float _S3298 = p1_6 * _S3294;
    float _S3299 = v_16 * _S3294;
    float _S3300 = p2_6 * _S3295;
    float _S3301 = v_16 * _S3295;
    float _S3302 = _S3297.x + _S3297.y;
    float _S3303 = r2_16 * _S3302;
    float _S3304 = r2_16 * _S3303;
    float _S3305 = r2_16 * _S3304;
    float _S3306 = sy1_6 * _S3294 + _S3298 + sx1_6 * _S3295 + _S3300 + _S2948 * _S3302 + _S2947 * _S3303 + _S2946 * _S3304 + k4_6 * _S3305;
    float _S3307 = v_16 * _S3306;
    float _S3308 = u_16 * _S3306;
    float2  _S3309 = make_float2 (r2_16 * _S3295 + _S3277, r2_16 * _S3294 + _S3272);
    float2  _S3310 = make_float2 (_S2957 * _S3294 + 2.0f * (u_16 * _S3301 + _S3271) + _S3275, 2.0f * (u_16 * _S3299 + _S3270) + _S2953 * _S3295 + _S3278);
    float4  _S3311 = make_float4 (_S3303 + _S3279, _S3304 + _S3276, _S3305 + _S3274, r2_16 * _S3305 + _S3273);
    float3  _S3312 = _S3076.differential_0 + make_float3 (_S3293.x, _S3293.y, _S3290);
    float2  _S3313 = _S2949 * _S3296 + make_float2 (_S2954 * _S3299 + _S2952 * _S3300 + 2.0f * (u_16 * _S3300) + _S2950 * _S3301 + _S3308 + _S3308, _S2956 * _S3298 + 2.0f * (v_16 * _S3298) + _S2955 * _S3294 + _S2951 * _S3295 + _S3307 + _S3307);
    CameraDistortion_0 _S3314 = CameraDistortion_x24_syn_dzero_0();
    (&_S3314)->thin_prism_coeffs_0 = _S3309;
    (&_S3314)->tangential_coeffs_0 = _S3310;
    (&_S3314)->radial_coeffs_0 = _S3311;
    float2  _S3315 = _S2933 * _S3313;
    float2  _S3316 = _S2944 * _S3313;
    float _S3317 = _S3315.x + _S3315.y;
    if(_S2937)
    {
        float _S3318 = _S3317 / _S2939;
        float _S3319 = _S2940 * - _S3318;
        float _S3320 = _S2936 * (0.3333333432674408f * - (_S2935 * _S3318));
        k_6 = _S3320 + _S3320;
        _S2938 = _S3319;
        _S2939 = 0.0f;
    }
    else
    {
        float _S3321 = _S3317 / _S2938;
        float _S3322 = _S2936 * - _S3321;
        k_6 = _S2934 * _S3321;
        _S2938 = 0.0f;
        _S2939 = _S3322;
    }
    DiffPair_float_0 _S3323;
    (&_S3323)->primal_0 = _S2934;
    (&_S3323)->differential_0 = 0.0f;
    DiffPair_float_0 _S3324;
    (&_S3324)->primal_0 = _S2935;
    (&_S3324)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3323, &_S3324, k_6);
    float _S3325 = _S3324.differential_0 + _S2938;
    float _S3326 = _S3323.differential_0 + _S2939;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3327;
    (&_S3327)->primal_0 = _S2933;
    (&_S3327)->differential_0 = _S3176;
    s_bwd_length_impl_1(&_S3327, _S3326);
    float2  _S3328 = _S3327.differential_0 + _S3316;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3329;
    (&_S3329)->primal_0 = vert2_c_5;
    (&_S3329)->differential_0 = _S3068;
    s_bwd_length_impl_0(&_S3329, _S3157.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3330;
    (&_S3330)->primal_0 = vert1_c_5;
    (&_S3330)->differential_0 = _S3068;
    s_bwd_length_impl_0(&_S3330, _S3157.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3331;
    (&_S3331)->primal_0 = vert0_c_5;
    (&_S3331)->differential_0 = _S3068;
    s_bwd_length_impl_0(&_S3331, _S3157.differential_0.x);
    float3  _S3332 = _S3329.differential_0 + _S3268;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3333;
    (&_S3333)->primal_0 = R_24;
    (&_S3333)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3334;
    (&_S3334)->primal_0 = vert2_1;
    (&_S3334)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3333, &_S3334, _S3332);
    float3  _S3335 = _S3330.differential_0 + _S3312;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3336;
    (&_S3336)->primal_0 = R_24;
    (&_S3336)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3337;
    (&_S3337)->primal_0 = vert1_1;
    (&_S3337)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3336, &_S3337, _S3335);
    float3  _S3338 = _S3331.differential_0 + _S3214 + make_float3 (_S3328.x, _S3328.y, _S3325);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3339;
    (&_S3339)->primal_0 = R_24;
    (&_S3339)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3340;
    (&_S3340)->primal_0 = vert0_1;
    (&_S3340)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3339, &_S3340, _S3338);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3341;
    (&_S3341)->primal_0 = _S2923;
    (&_S3341)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3342;
    (&_S3342)->primal_0 = _S2928;
    (&_S3342)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3341, &_S3342, _S3334.differential_0);
    float _S3343 = - _S3342.differential_0.y;
    float _S3344 = _S2927 * _S3342.differential_0.x;
    float _S3345 = - (_S2919 * _S3342.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3346;
    (&_S3346)->primal_0 = _S2923;
    (&_S3346)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3347;
    (&_S3347)->primal_0 = _S2926;
    (&_S3347)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3346, &_S3347, _S3337.differential_0);
    float _S3348 = _S2919 * _S3347.differential_0.x;
    float _S3349 = _S2925 * _S3347.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3350;
    (&_S3350)->primal_0 = _S2923;
    (&_S3350)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3351;
    (&_S3351)->primal_0 = _S2924;
    (&_S3351)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3350, &_S3351, _S3340.differential_0);
    Matrix<float, 3, 3>  _S3352 = transpose_0(_S3341.differential_0 + _S3346.differential_0 + _S3350.differential_0);
    float _S3353 = 2.0f * - _S3352.rows[int(2)].z;
    float _S3354 = 2.0f * _S3352.rows[int(2)].y;
    float _S3355 = 2.0f * _S3352.rows[int(2)].x;
    float _S3356 = 2.0f * _S3352.rows[int(1)].z;
    float _S3357 = 2.0f * - _S3352.rows[int(1)].y;
    float _S3358 = 2.0f * _S3352.rows[int(1)].x;
    float _S3359 = 2.0f * _S3352.rows[int(0)].z;
    float _S3360 = 2.0f * _S3352.rows[int(0)].y;
    float _S3361 = 2.0f * - _S3352.rows[int(0)].x;
    float _S3362 = - _S3358 + _S3360;
    float _S3363 = _S3355 + - _S3359;
    float _S3364 = - _S3354 + _S3356;
    float _S3365 = _S3354 + _S3356;
    float _S3366 = _S3355 + _S3359;
    float _S3367 = _S3358 + _S3360;
    float _S3368 = quat_25.w * (_S3357 + _S3361);
    float _S3369 = quat_25.z * (_S3353 + _S3361);
    float _S3370 = quat_25.y * (_S3353 + _S3357);
    float _S3371 = quat_25.x * _S3362 + quat_25.z * _S3365 + quat_25.y * _S3366 + _S3368 + _S3368;
    float _S3372 = quat_25.x * _S3363 + quat_25.w * _S3365 + quat_25.y * _S3367 + _S3369 + _S3369;
    float _S3373 = quat_25.x * _S3364 + quat_25.w * _S3366 + quat_25.z * _S3367 + _S3370 + _S3370;
    float _S3374 = quat_25.w * _S3362 + quat_25.z * _S3363 + quat_25.y * _S3364;
    float _S3375 = _S3345 + _S3348;
    float _S3376 = 0.5f * - _S3375;
    float _S3377 = _S3343 + _S3347.differential_0.y;
    DiffPair_float_0 _S3378;
    (&_S3378)->primal_0 = _S2920;
    (&_S3378)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3378, _S3377);
    float _S3379 = _S3376 + _S3378.differential_0;
    float _S3380 = _S3344 + _S3349 + _S3351.differential_0.x;
    DiffPair_float_0 _S3381;
    (&_S3381)->primal_0 = _S2918;
    (&_S3381)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3381, _S3380);
    float _S3382 = _S3376 + _S3381.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3383;
    (&_S3383)->primal_0 = mean_c_20;
    (&_S3383)->differential_0 = _S3068;
    s_bwd_length_impl_0(&_S3383, 0.0f);
    float3  _S3384 = _S3383.differential_0 + _S3071.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3385;
    (&_S3385)->primal_0 = R_24;
    (&_S3385)->differential_0 = _S3151;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3386;
    (&_S3386)->primal_0 = mean_24;
    (&_S3386)->differential_0 = _S3068;
    s_bwd_prop_mul_1(&_S3385, &_S3386, _S3384);
    float3  _S3387 = _S3332 + _S3335 + _S3338 + _S3384 + _S3154.differential_0;
    Matrix<float, 3, 3>  _S3388 = _S3333.differential_0 + _S3336.differential_0 + _S3339.differential_0 + _S3385.differential_0 + _S3155;
    float3  _S3389 = make_float3 (_S3382, _S3379, _S3375);
    float4  _S3390 = make_float4 (0.0f);
    *&((&_S3390)->w) = _S3371;
    *&((&_S3390)->z) = _S3372;
    *&((&_S3390)->y) = _S3373;
    *&((&_S3390)->x) = _S3374;
    float4  _S3391 = _S3390;
    float3  _S3392 = _S3334.differential_0 + _S3337.differential_0 + _S3340.differential_0 + _S3386.differential_0 + _S3149;
    *v_mean_8 = _S3392;
    *v_quat_7 = _S3391;
    *v_scale_7 = _S3389;
    *v_hardness_1 = _S3211;
    (*v_sh_coeffs_6)[int(0)] = _S3219;
    (*v_sh_coeffs_6)[int(1)] = _S3220;
    (*v_sh_coeffs_6)[int(2)] = _S3221;
    (*v_sh_coeffs_6)[int(3)] = _S3222;
    (*v_sh_coeffs_6)[int(4)] = _S3223;
    (*v_sh_coeffs_6)[int(5)] = _S3224;
    (*v_sh_coeffs_6)[int(6)] = _S3225;
    (*v_sh_coeffs_6)[int(7)] = _S3226;
    (*v_sh_coeffs_6)[int(8)] = _S3227;
    (*v_sh_coeffs_6)[int(9)] = _S3228;
    (*v_sh_coeffs_6)[int(10)] = _S3229;
    (*v_sh_coeffs_6)[int(11)] = _S3230;
    (*v_sh_coeffs_6)[int(12)] = _S3231;
    (*v_sh_coeffs_6)[int(13)] = _S3232;
    (*v_sh_coeffs_6)[int(14)] = _S3233;
    (*v_sh_coeffs_6)[int(15)] = _S3234;
    (*v_ch_coeffs_1)[int(0)] = _S3216;
    (*v_ch_coeffs_1)[int(1)] = _S3217;
    *v_R_7 = _S3388;
    *v_t_7 = _S3387;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_19)
{
    DiffPair_float_0 _S3393 = *dpx_16;
    bool _S3394;
    if(((*dpx_16).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S3394 = ((*dpx_16).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S3394 = false;
    }
    float _S3395;
    if(_S3394)
    {
        _S3395 = dOut_19;
    }
    else
    {
        _S3395 = 0.0f;
    }
    dpx_16->primal_0 = _S3393.primal_0;
    dpx_16->differential_0 = _S3395;
    DiffPair_float_0 _S3396 = *dpMin_0;
    if((_S3393.primal_0) < ((*dpMin_0).primal_0))
    {
        _S3395 = dOut_19;
    }
    else
    {
        _S3395 = 0.0f;
    }
    dpMin_0->primal_0 = _S3396.primal_0;
    dpMin_0->differential_0 = _S3395;
    DiffPair_float_0 _S3397 = *dpMax_0;
    if(((*dpx_16).primal_0) > ((*dpMax_0).primal_0))
    {
        _S3395 = dOut_19;
    }
    else
    {
        _S3395 = 0.0f;
    }
    dpMax_0->primal_0 = _S3397.primal_0;
    dpMax_0->differential_0 = _S3395;
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
        DiffPair_float_0 _S3398 = *dpx_17;
        float _S3399 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S3399;
        float _S3400 = val_0 * (F32_log((_S3398.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3400;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3401 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3401))));
    float2  _S3402 = p_0 - v0_0;
    float2  _S3403 = normalize_1(e0_6);
    float2  _S3404 = p_0 - v1_0;
    float2  _S3405 = normalize_1(e1_6);
    float2  _S3406 = p_0 - v2_0;
    float2  _S3407 = normalize_1(e2_2);
    float _S3408 = hardness_6.x;
    float _S3409 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3402.x * _S3403.y - _S3402.y * _S3403.x)), (se_0 * (_S3404.x * _S3405.y - _S3404.y * _S3405.x))))), (se_0 * (_S3406.x * _S3407.y - _S3406.y * _S3407.x)))) / ((F32_abs((_S3401))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S3409))));
    float _S3410;
    if(a_1 <= 0.0f)
    {
        _S3410 = 0.0f;
    }
    else
    {
        _S3410 = (F32_min(((F32_pow((a_1), (_S3409)))), (0.99900001287460327f)));
    }
    return _S3408 * _S3410;
}

inline __device__ float s_primal_ctx_abs_0(float _S3411)
{
    return (F32_abs((_S3411)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S3412, float _S3413, float _S3414)
{
    return clamp_0(_S3412, _S3413, _S3414);
}

inline __device__ float s_primal_ctx_exp2_0(float _S3415)
{
    return (F32_exp2((_S3415)));
}

inline __device__ float s_primal_ctx_pow_0(float _S3416, float _S3417)
{
    return (F32_pow((_S3416), (_S3417)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S3418, DiffPair_float_0 * _S3419, float _S3420)
{
    _d_pow_0(_S3418, _S3419, _S3420);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S3421, DiffPair_float_0 * _S3422, DiffPair_float_0 * _S3423, float _S3424)
{
    _d_clamp_0(_S3421, _S3422, _S3423, _S3424);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_9)
{
    float _S3425 = length_0((*dpx_18).primal_0);
    float2  _S3426 = (*dpx_18).primal_0 * _s_dOut_9;
    float2  _S3427 = make_float2 (1.0f / _S3425) * _s_dOut_9;
    float _S3428 = - ((_S3426.x + _S3426.y) / (_S3425 * _S3425));
    float2  _S3429 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3430;
    (&_S3430)->primal_0 = (*dpx_18).primal_0;
    (&_S3430)->differential_0 = _S3429;
    s_bwd_length_impl_1(&_S3430, _S3428);
    float2  _S3431 = _S3427 + _S3430.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S3431;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3432, float2  _S3433)
{
    s_bwd_prop_normalize_impl_1(_S3432, _S3433);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_10)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S3434 = e0_7.x;
    float _S3435 = e1_7.y;
    float _S3436 = e0_7.y;
    float _S3437 = e1_7.x;
    float _S3438 = _S3434 * _S3435 - _S3436 * _S3437;
    float se_1 = float((F32_sign((_S3438))));
    float2  _S3439 = p_1 - (*dpv0_0).primal_0;
    float2  _S3440 = normalize_1(e0_7);
    float _S3441 = _S3439.x;
    float _S3442 = _S3440.y;
    float _S3443 = _S3439.y;
    float _S3444 = _S3440.x;
    float de0_0 = se_1 * (_S3441 * _S3442 - _S3443 * _S3444);
    float2  _S3445 = p_1 - (*dpv1_0).primal_0;
    float2  _S3446 = normalize_1(e1_7);
    float _S3447 = _S3445.x;
    float _S3448 = _S3446.y;
    float _S3449 = _S3445.y;
    float _S3450 = _S3446.x;
    float de1_0 = se_1 * (_S3447 * _S3448 - _S3449 * _S3450);
    float2  _S3451 = p_1 - (*dpv2_0).primal_0;
    float2  _S3452 = normalize_1(e2_3);
    float _S3453 = _S3451.x;
    float _S3454 = _S3452.y;
    float _S3455 = _S3451.y;
    float _S3456 = _S3452.x;
    float de2_0 = se_1 * (_S3453 * _S3454 - _S3455 * _S3456);
    float _S3457 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S3458 = s_primal_ctx_max_0(_S3457, de2_0);
    float _S3459 = s_primal_ctx_abs_0(_S3438);
    float _S3460 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S3459 / _S3460;
    float _S3461 = _S3460 * _S3460;
    float _S3462 = (*dphardness_0).primal_0.x;
    float _S3463 = (*dphardness_0).primal_0.y;
    float _S3464 = dmax_0 * dmax_0;
    float _S3465 = 1.0f + _S3458 / dmax_0;
    float _S3466 = 1.0f - s_primal_ctx_clamp_0(_S3463, 0.00499999988824129f, 0.98000001907348633f);
    float _S3467 = -1.0f / _S3466;
    float _S3468 = _S3466 * _S3466;
    float _S3469 = 1.0f - s_primal_ctx_exp2_0(_S3467);
    float a_2 = 1.0f - _S3465 * _S3469;
    bool _S3470 = a_2 <= 0.0f;
    float _S3471;
    float _S3472;
    if(_S3470)
    {
        _S3471 = 0.0f;
        _S3472 = 0.0f;
    }
    else
    {
        float _S3473 = s_primal_ctx_pow_0(a_2, _S3466);
        _S3471 = s_primal_ctx_min_0(_S3473, 0.99900001287460327f);
        _S3472 = _S3473;
    }
    float _S3474 = _S3462 * _s_dOut_10;
    float _S3475 = _S3471 * _s_dOut_10;
    if(_S3470)
    {
        _S3471 = 0.0f;
        _S3472 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S3476;
        (&_S3476)->primal_0 = _S3472;
        (&_S3476)->differential_0 = 0.0f;
        DiffPair_float_0 _S3477;
        (&_S3477)->primal_0 = 0.99900001287460327f;
        (&_S3477)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S3476, &_S3477, _S3474);
        DiffPair_float_0 _S3478;
        (&_S3478)->primal_0 = a_2;
        (&_S3478)->differential_0 = 0.0f;
        DiffPair_float_0 _S3479;
        (&_S3479)->primal_0 = _S3466;
        (&_S3479)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S3478, &_S3479, _S3476.differential_0);
        _S3471 = _S3478.differential_0;
        _S3472 = _S3479.differential_0;
    }
    float _S3480 = - _S3471;
    float _S3481 = _S3469 * _S3480;
    float _S3482 = - (_S3465 * _S3480);
    DiffPair_float_0 _S3483;
    (&_S3483)->primal_0 = _S3467;
    (&_S3483)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3483, _S3482);
    float _S3484 = - (-1.0f * - (_S3483.differential_0 / _S3468) + _S3472);
    float _S3485 = _S3481 / _S3464;
    float s_diff_dmax_T_0 = _S3458 * - _S3485;
    float _S3486 = dmax_0 * _S3485;
    DiffPair_float_0 _S3487;
    (&_S3487)->primal_0 = _S3463;
    (&_S3487)->differential_0 = 0.0f;
    DiffPair_float_0 _S3488;
    (&_S3488)->primal_0 = 0.00499999988824129f;
    (&_S3488)->differential_0 = 0.0f;
    DiffPair_float_0 _S3489;
    (&_S3489)->primal_0 = 0.98000001907348633f;
    (&_S3489)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3487, &_S3488, &_S3489, _S3484);
    float _S3490 = s_diff_dmax_T_0 / _S3461;
    float _S3491 = _S3459 * - _S3490;
    float _S3492 = _S3460 * _S3490;
    float2  _S3493 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3494;
    (&_S3494)->primal_0 = e2_3;
    (&_S3494)->differential_0 = _S3493;
    s_bwd_length_impl_1(&_S3494, _S3491);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3495;
    (&_S3495)->primal_0 = e1_7;
    (&_S3495)->differential_0 = _S3493;
    s_bwd_length_impl_1(&_S3495, _S3491);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3496;
    (&_S3496)->primal_0 = e0_7;
    (&_S3496)->differential_0 = _S3493;
    s_bwd_length_impl_1(&_S3496, _S3491);
    DiffPair_float_0 _S3497;
    (&_S3497)->primal_0 = _S3438;
    (&_S3497)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3497, _S3492);
    DiffPair_float_0 _S3498;
    (&_S3498)->primal_0 = _S3457;
    (&_S3498)->differential_0 = 0.0f;
    DiffPair_float_0 _S3499;
    (&_S3499)->primal_0 = de2_0;
    (&_S3499)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3498, &_S3499, _S3486);
    DiffPair_float_0 _S3500;
    (&_S3500)->primal_0 = de0_0;
    (&_S3500)->differential_0 = 0.0f;
    DiffPair_float_0 _S3501;
    (&_S3501)->primal_0 = de1_0;
    (&_S3501)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3500, &_S3501, _S3498.differential_0);
    float _S3502 = se_1 * _S3499.differential_0;
    float _S3503 = - _S3502;
    float _S3504 = _S3456 * _S3503;
    float _S3505 = _S3454 * _S3502;
    float2  _S3506 = make_float2 (_S3455 * _S3503, _S3453 * _S3502);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3507;
    (&_S3507)->primal_0 = e2_3;
    (&_S3507)->differential_0 = _S3493;
    s_bwd_normalize_impl_1(&_S3507, _S3506);
    float2  _S3508 = - make_float2 (_S3505, _S3504);
    float _S3509 = se_1 * _S3501.differential_0;
    float _S3510 = - _S3509;
    float _S3511 = _S3450 * _S3510;
    float _S3512 = _S3448 * _S3509;
    float2  _S3513 = make_float2 (_S3449 * _S3510, _S3447 * _S3509);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3514;
    (&_S3514)->primal_0 = e1_7;
    (&_S3514)->differential_0 = _S3493;
    s_bwd_normalize_impl_1(&_S3514, _S3513);
    float2  _S3515 = - make_float2 (_S3512, _S3511);
    float _S3516 = se_1 * _S3500.differential_0;
    float _S3517 = - _S3516;
    float _S3518 = _S3444 * _S3517;
    float _S3519 = _S3442 * _S3516;
    float2  _S3520 = make_float2 (_S3443 * _S3517, _S3441 * _S3516);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3521;
    (&_S3521)->primal_0 = e0_7;
    (&_S3521)->differential_0 = _S3493;
    s_bwd_normalize_impl_1(&_S3521, _S3520);
    float2  _S3522 = - make_float2 (_S3519, _S3518);
    float _S3523 = - _S3497.differential_0;
    float2  _S3524 = _S3494.differential_0 + _S3507.differential_0;
    float2  _S3525 = - _S3524;
    float2  _S3526 = _S3495.differential_0 + _S3514.differential_0 + make_float2 (_S3436 * _S3523, _S3434 * _S3497.differential_0);
    float2  _S3527 = - _S3526;
    float2  _S3528 = _S3496.differential_0 + _S3521.differential_0 + make_float2 (_S3435 * _S3497.differential_0, _S3437 * _S3523);
    float2  _S3529 = - _S3528;
    float2  _S3530 = make_float2 (_S3475, _S3487.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S3530;
    float2  _S3531 = _S3508 + _S3525 + _S3526;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S3531;
    float2  _S3532 = _S3515 + _S3527 + _S3528;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S3532;
    float2  _S3533 = _S3522 + _S3524 + _S3529;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S3533;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3534, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3535, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3536, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3537, float2  _S3538, float _S3539)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S3534, _S3535, _S3536, _S3537, _S3538, _S3539);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S3540 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S3540;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S3540;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S3540;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S3540;
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
    float2  _S3541 = p_3 - v0_2;
    float2  _S3542 = p_3 - v1_2;
    float2  _S3543 = p_3 - v2_2;
    float _S3544 = e0_8.x;
    float _S3545 = e1_8.y;
    float _S3546 = e0_8.y;
    float _S3547 = e1_8.x;
    float _S3548 = _S3544 * _S3545 - _S3546 * _S3547;
    float se_2 = float((F32_sign((_S3548))));
    float _S3549 = hardness_8.x;
    float _S3550 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S3541.x * _S3546 - _S3541.y * _S3544)), (se_2 * (_S3542.x * _S3545 - _S3542.y * _S3547))))), (se_2 * (_S3543.x * e2_4.y - _S3543.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S3541 - e0_8 * make_float2 (clamp_0(dot_1(_S3541, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S3542 - e1_8 * make_float2 (clamp_0(dot_1(_S3542, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S3543 - e2_4 * make_float2 (clamp_0(dot_1(_S3543, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S3548))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S3550))));
    float _S3551;
    if(a_3 <= 0.0f)
    {
        _S3551 = 0.0f;
    }
    else
    {
        _S3551 = (F32_min(((F32_pow((a_3), (_S3550)))), (0.99900001287460327f)));
    }
    return _S3549 * _S3551;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S3552, float2  _S3553)
{
    return dot_1(_S3552, _S3553);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3554, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3555, float _S3556)
{
    _d_dot_1(_S3554, _S3555, _S3556);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_11)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S3557 = p_4 - (*dpv0_1).primal_0;
    float _S3558 = s_primal_ctx_dot_1(_S3557, e0_9);
    float _S3559 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S3560 = _S3558 / _S3559;
    float _S3561 = _S3559 * _S3559;
    float _S3562 = s_primal_ctx_clamp_0(_S3560, 0.0f, 1.0f);
    float2  _S3563 = make_float2 (_S3562);
    float2  _S3564 = _S3557 - e0_9 * make_float2 (_S3562);
    float _S3565 = length_0(_S3564);
    float2  _S3566 = p_4 - (*dpv1_1).primal_0;
    float _S3567 = s_primal_ctx_dot_1(_S3566, e1_9);
    float _S3568 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S3569 = _S3567 / _S3568;
    float _S3570 = _S3568 * _S3568;
    float _S3571 = s_primal_ctx_clamp_0(_S3569, 0.0f, 1.0f);
    float2  _S3572 = make_float2 (_S3571);
    float2  _S3573 = _S3566 - e1_9 * make_float2 (_S3571);
    float _S3574 = length_0(_S3573);
    float2  _S3575 = p_4 - (*dpv2_1).primal_0;
    float _S3576 = s_primal_ctx_dot_1(_S3575, e2_5);
    float _S3577 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S3578 = _S3576 / _S3577;
    float _S3579 = _S3577 * _S3577;
    float _S3580 = s_primal_ctx_clamp_0(_S3578, 0.0f, 1.0f);
    float2  _S3581 = make_float2 (_S3580);
    float2  _S3582 = _S3575 - e2_5 * make_float2 (_S3580);
    float _S3583 = length_0(_S3582);
    float _S3584 = e0_9.x;
    float _S3585 = e1_9.y;
    float _S3586 = e0_9.y;
    float _S3587 = e1_9.x;
    float _S3588 = _S3584 * _S3585 - _S3586 * _S3587;
    float se_3 = float((F32_sign((_S3588))));
    float _S3589 = _S3557.x;
    float _S3590 = _S3557.y;
    float s0_0 = se_3 * (_S3589 * _S3586 - _S3590 * _S3584);
    float _S3591 = _S3566.x;
    float _S3592 = _S3566.y;
    float s1_0 = se_3 * (_S3591 * _S3585 - _S3592 * _S3587);
    float _S3593 = _S3575.x;
    float _S3594 = e2_5.y;
    float _S3595 = _S3575.y;
    float _S3596 = e2_5.x;
    float s2_0 = se_3 * (_S3593 * _S3594 - _S3595 * _S3596);
    float _S3597 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S3597, s2_0)))));
    float _S3598 = s_primal_ctx_min_0(_S3565, _S3574);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S3598, _S3583);
    float _S3599 = s_primal_ctx_abs_0(_S3588);
    float _S3600 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S3599 / _S3600;
    float _S3601 = _S3600 * _S3600;
    float _S3602 = (*dphardness_1).primal_0.x;
    float _S3603 = (*dphardness_1).primal_0.y;
    float _S3604 = dmax_1 * dmax_1;
    float _S3605 = 1.0f + dv_0 / dmax_1;
    float _S3606 = 1.0f - s_primal_ctx_clamp_0(_S3603, 0.00499999988824129f, 0.98000001907348633f);
    float _S3607 = -1.0f / _S3606;
    float _S3608 = _S3606 * _S3606;
    float _S3609 = 1.0f - s_primal_ctx_exp2_0(_S3607);
    float a_4 = 1.0f - _S3605 * _S3609;
    bool _S3610 = a_4 <= 0.0f;
    float _S3611;
    float _S3612;
    if(_S3610)
    {
        _S3611 = 0.0f;
        _S3612 = 0.0f;
    }
    else
    {
        float _S3613 = s_primal_ctx_pow_0(a_4, _S3606);
        _S3611 = s_primal_ctx_min_0(_S3613, 0.99900001287460327f);
        _S3612 = _S3613;
    }
    float _S3614 = _S3602 * _s_dOut_11;
    float _S3615 = _S3611 * _s_dOut_11;
    if(_S3610)
    {
        _S3611 = 0.0f;
        _S3612 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S3616;
        (&_S3616)->primal_0 = _S3612;
        (&_S3616)->differential_0 = 0.0f;
        DiffPair_float_0 _S3617;
        (&_S3617)->primal_0 = 0.99900001287460327f;
        (&_S3617)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S3616, &_S3617, _S3614);
        DiffPair_float_0 _S3618;
        (&_S3618)->primal_0 = a_4;
        (&_S3618)->differential_0 = 0.0f;
        DiffPair_float_0 _S3619;
        (&_S3619)->primal_0 = _S3606;
        (&_S3619)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S3618, &_S3619, _S3616.differential_0);
        _S3611 = _S3618.differential_0;
        _S3612 = _S3619.differential_0;
    }
    float _S3620 = - _S3611;
    float _S3621 = _S3609 * _S3620;
    float _S3622 = - (_S3605 * _S3620);
    DiffPair_float_0 _S3623;
    (&_S3623)->primal_0 = _S3607;
    (&_S3623)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3623, _S3622);
    float _S3624 = - (-1.0f * - (_S3623.differential_0 / _S3608) + _S3612);
    float _S3625 = _S3621 / _S3604;
    float s_diff_dmax_T_1 = dv_0 * - _S3625;
    float s_diff_dv_T_0 = dmax_1 * _S3625;
    DiffPair_float_0 _S3626;
    (&_S3626)->primal_0 = _S3603;
    (&_S3626)->differential_0 = 0.0f;
    DiffPair_float_0 _S3627;
    (&_S3627)->primal_0 = 0.00499999988824129f;
    (&_S3627)->differential_0 = 0.0f;
    DiffPair_float_0 _S3628;
    (&_S3628)->primal_0 = 0.98000001907348633f;
    (&_S3628)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3626, &_S3627, &_S3628, _S3624);
    float _S3629 = s_diff_dmax_T_1 / _S3601;
    float _S3630 = _S3599 * - _S3629;
    float _S3631 = _S3600 * _S3629;
    float2  _S3632 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3633;
    (&_S3633)->primal_0 = e2_5;
    (&_S3633)->differential_0 = _S3632;
    s_bwd_length_impl_1(&_S3633, _S3630);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3634;
    (&_S3634)->primal_0 = e1_9;
    (&_S3634)->differential_0 = _S3632;
    s_bwd_length_impl_1(&_S3634, _S3630);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3635;
    (&_S3635)->primal_0 = e0_9;
    (&_S3635)->differential_0 = _S3632;
    s_bwd_length_impl_1(&_S3635, _S3630);
    DiffPair_float_0 _S3636;
    (&_S3636)->primal_0 = _S3588;
    (&_S3636)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3636, _S3631);
    float _S3637 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S3638;
    (&_S3638)->primal_0 = _S3598;
    (&_S3638)->differential_0 = 0.0f;
    DiffPair_float_0 _S3639;
    (&_S3639)->primal_0 = _S3583;
    (&_S3639)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3638, &_S3639, _S3637);
    DiffPair_float_0 _S3640;
    (&_S3640)->primal_0 = _S3565;
    (&_S3640)->differential_0 = 0.0f;
    DiffPair_float_0 _S3641;
    (&_S3641)->primal_0 = _S3574;
    (&_S3641)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3640, &_S3641, _S3638.differential_0);
    DiffPair_float_0 _S3642;
    (&_S3642)->primal_0 = _S3597;
    (&_S3642)->differential_0 = 0.0f;
    DiffPair_float_0 _S3643;
    (&_S3643)->primal_0 = s2_0;
    (&_S3643)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3642, &_S3643, 0.0f);
    DiffPair_float_0 _S3644;
    (&_S3644)->primal_0 = s0_0;
    (&_S3644)->differential_0 = 0.0f;
    DiffPair_float_0 _S3645;
    (&_S3645)->primal_0 = s1_0;
    (&_S3645)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3644, &_S3645, _S3642.differential_0);
    float _S3646 = se_3 * _S3643.differential_0;
    float _S3647 = - _S3646;
    float _S3648 = _S3595 * _S3647;
    float _S3649 = _S3596 * _S3647;
    float _S3650 = _S3593 * _S3646;
    float _S3651 = _S3594 * _S3646;
    float _S3652 = se_3 * _S3645.differential_0;
    float _S3653 = - _S3652;
    float _S3654 = _S3587 * _S3653;
    float _S3655 = _S3585 * _S3652;
    float _S3656 = se_3 * _S3644.differential_0;
    float _S3657 = - _S3656;
    float _S3658 = _S3584 * _S3657;
    float _S3659 = _S3586 * _S3656;
    float _S3660 = - _S3636.differential_0;
    float _S3661 = _S3592 * _S3653 + _S3586 * _S3660;
    float _S3662 = _S3589 * _S3656 + _S3587 * _S3660;
    float _S3663 = _S3591 * _S3652 + _S3584 * _S3636.differential_0;
    float _S3664 = _S3590 * _S3657 + _S3585 * _S3636.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3665;
    (&_S3665)->primal_0 = _S3582;
    (&_S3665)->differential_0 = _S3632;
    s_bwd_length_impl_1(&_S3665, _S3639.differential_0);
    float2  _S3666 = - _S3665.differential_0;
    float2  _S3667 = e2_5 * _S3666;
    float2  _S3668 = _S3581 * _S3666;
    float _S3669 = _S3667.x + _S3667.y;
    DiffPair_float_0 _S3670;
    (&_S3670)->primal_0 = _S3578;
    (&_S3670)->differential_0 = 0.0f;
    DiffPair_float_0 _S3671;
    (&_S3671)->primal_0 = 0.0f;
    (&_S3671)->differential_0 = 0.0f;
    DiffPair_float_0 _S3672;
    (&_S3672)->primal_0 = 1.0f;
    (&_S3672)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3670, &_S3671, &_S3672, _S3669);
    float _S3673 = _S3670.differential_0 / _S3579;
    float _S3674 = _S3576 * - _S3673;
    float _S3675 = _S3577 * _S3673;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3676;
    (&_S3676)->primal_0 = e2_5;
    (&_S3676)->differential_0 = _S3632;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3677;
    (&_S3677)->primal_0 = e2_5;
    (&_S3677)->differential_0 = _S3632;
    s_bwd_prop_dot_1(&_S3676, &_S3677, _S3674);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3678;
    (&_S3678)->primal_0 = _S3575;
    (&_S3678)->differential_0 = _S3632;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3679;
    (&_S3679)->primal_0 = e2_5;
    (&_S3679)->differential_0 = _S3632;
    s_bwd_prop_dot_1(&_S3678, &_S3679, _S3675);
    float2  _S3680 = - (_S3665.differential_0 + _S3678.differential_0 + make_float2 (_S3651, _S3649));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3681;
    (&_S3681)->primal_0 = _S3573;
    (&_S3681)->differential_0 = _S3632;
    s_bwd_length_impl_1(&_S3681, _S3641.differential_0);
    float2  _S3682 = - _S3681.differential_0;
    float2  _S3683 = e1_9 * _S3682;
    float2  _S3684 = _S3572 * _S3682;
    float _S3685 = _S3683.x + _S3683.y;
    DiffPair_float_0 _S3686;
    (&_S3686)->primal_0 = _S3569;
    (&_S3686)->differential_0 = 0.0f;
    DiffPair_float_0 _S3687;
    (&_S3687)->primal_0 = 0.0f;
    (&_S3687)->differential_0 = 0.0f;
    DiffPair_float_0 _S3688;
    (&_S3688)->primal_0 = 1.0f;
    (&_S3688)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3686, &_S3687, &_S3688, _S3685);
    float _S3689 = _S3686.differential_0 / _S3570;
    float _S3690 = _S3567 * - _S3689;
    float _S3691 = _S3568 * _S3689;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3692;
    (&_S3692)->primal_0 = e1_9;
    (&_S3692)->differential_0 = _S3632;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3693;
    (&_S3693)->primal_0 = e1_9;
    (&_S3693)->differential_0 = _S3632;
    s_bwd_prop_dot_1(&_S3692, &_S3693, _S3690);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3694;
    (&_S3694)->primal_0 = _S3566;
    (&_S3694)->differential_0 = _S3632;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3695;
    (&_S3695)->primal_0 = e1_9;
    (&_S3695)->differential_0 = _S3632;
    s_bwd_prop_dot_1(&_S3694, &_S3695, _S3691);
    float2  _S3696 = - (_S3681.differential_0 + _S3694.differential_0 + make_float2 (_S3655, _S3654));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3697;
    (&_S3697)->primal_0 = _S3564;
    (&_S3697)->differential_0 = _S3632;
    s_bwd_length_impl_1(&_S3697, _S3640.differential_0);
    float2  _S3698 = - _S3697.differential_0;
    float2  _S3699 = e0_9 * _S3698;
    float2  _S3700 = _S3563 * _S3698;
    float _S3701 = _S3699.x + _S3699.y;
    DiffPair_float_0 _S3702;
    (&_S3702)->primal_0 = _S3560;
    (&_S3702)->differential_0 = 0.0f;
    DiffPair_float_0 _S3703;
    (&_S3703)->primal_0 = 0.0f;
    (&_S3703)->differential_0 = 0.0f;
    DiffPair_float_0 _S3704;
    (&_S3704)->primal_0 = 1.0f;
    (&_S3704)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3702, &_S3703, &_S3704, _S3701);
    float _S3705 = _S3702.differential_0 / _S3561;
    float _S3706 = _S3558 * - _S3705;
    float _S3707 = _S3559 * _S3705;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3708;
    (&_S3708)->primal_0 = e0_9;
    (&_S3708)->differential_0 = _S3632;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3709;
    (&_S3709)->primal_0 = e0_9;
    (&_S3709)->differential_0 = _S3632;
    s_bwd_prop_dot_1(&_S3708, &_S3709, _S3706);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3710;
    (&_S3710)->primal_0 = _S3557;
    (&_S3710)->differential_0 = _S3632;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3711;
    (&_S3711)->primal_0 = e0_9;
    (&_S3711)->differential_0 = _S3632;
    s_bwd_prop_dot_1(&_S3710, &_S3711, _S3707);
    float2  _S3712 = - (_S3697.differential_0 + _S3710.differential_0 + make_float2 (_S3659, _S3658));
    float2  _S3713 = _S3633.differential_0 + _S3668 + _S3677.differential_0 + _S3676.differential_0 + _S3679.differential_0 + make_float2 (_S3648, _S3650);
    float2  _S3714 = - _S3713;
    float2  _S3715 = _S3634.differential_0 + _S3684 + _S3693.differential_0 + _S3692.differential_0 + _S3695.differential_0 + make_float2 (_S3661, _S3663);
    float2  _S3716 = - _S3715;
    float2  _S3717 = _S3635.differential_0 + _S3700 + _S3709.differential_0 + _S3708.differential_0 + _S3711.differential_0 + make_float2 (_S3664, _S3662);
    float2  _S3718 = - _S3717;
    float2  _S3719 = make_float2 (_S3615, _S3626.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S3719;
    float2  _S3720 = _S3680 + _S3714 + _S3715;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S3720;
    float2  _S3721 = _S3696 + _S3716 + _S3717;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S3721;
    float2  _S3722 = _S3712 + _S3713 + _S3718;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S3722;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3723, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3724, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3725, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3726, float2  _S3727, float _S3728)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S3723, _S3724, _S3725, _S3726, _S3727, _S3728);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S3729 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S3729;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S3729;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S3729;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S3729;
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
    float _S3730 = 0.3333333432674408f * dpdepth_1;
    float3  _S3731 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S3732 = make_float3 (0.0f);
    float3  _S3733 = _S3732;
    *&((&_S3733)->z) = _S3730;
    *&((&_S3733)->y) = _S3730;
    *&((&_S3733)->x) = _S3730;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S3733;
    FixedArray<float3 , 3>  _S3734;
    _S3734[int(0)] = _S3732;
    _S3734[int(1)] = _S3732;
    _S3734[int(2)] = _S3732;
    _S3734[int(2)] = _S3731;
    _S3734[int(1)] = _S3731;
    _S3734[int(0)] = _S3731;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S3734;
    float2  _S3735 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S3735;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S3735;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S3735;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3736, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3737, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3738, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S3739, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3740, float2  _S3741, float3  _S3742, float _S3743)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S3736, _S3737, _S3738, _S3739, _S3740, _S3741, _S3742, _S3743);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_8, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S3744 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S3744;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S3744;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S3744;
    float3  _S3745 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S3746 = { _S3745, _S3745, _S3745 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S3746;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S3745;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_8, v_color_0, v_depth_8);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_25, float4  quat_26, float3  scale_25, float2  hardness_10, FixedArray<float3 , 16>  * sh_coeffs_21, FixedArray<float3 , 2>  * ch_coeffs_6, Matrix<float, 3, 3>  R_25, float3  t_24, float fx_25, float fy_25, float cx_25, float cy_25, float4  radial_coeffs_28, float2  tangential_coeffs_28, float2  thin_prism_coeffs_28, uint image_width_21, uint image_height_21, float near_plane_14, float far_plane_14, int4  * aabb_xyxy_14, float * depth_16, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_4)
{
    for(;;)
    {
        float3  mean_c_21 = mul_0(R_25, mean_25) + t_24;
        float _S3747 = mean_c_21.z;
        bool _S3748;
        if(_S3747 < near_plane_14)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = _S3747 > far_plane_14;
        }
        if(_S3748)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3749 = scale_25.x;
        float sx_4 = (F32_exp((_S3749)));
        float _S3750 = scale_25.y;
        float sy_4 = (F32_exp((_S3750)));
        float sz_6 = scale_25.z - 0.5f * (_S3749 + _S3750);
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
        Matrix<float, 3, 3>  _S3751 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_47), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_47), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
        float3  vert0_2 = mul_0(_S3751, make_float3 (sx_4, 0.0f, 0.0f)) + mean_25;
        float3  vert1_2 = mul_0(_S3751, make_float3 (sx_4 * (-0.5f + sz_6), sy_4, 0.0f)) + mean_25;
        float3  vert2_2 = mul_0(_S3751, make_float3 (sx_4 * (-0.5f - sz_6), - sy_4, 0.0f)) + mean_25;
        float3  vert0_c_6 = mul_0(R_25, vert0_2) + t_24;
        float3  vert1_c_6 = mul_0(R_25, vert1_2) + t_24;
        float3  vert2_c_6 = mul_0(R_25, vert2_2) + t_24;
        float _S3752 = vert0_c_6.z;
        float _S3753 = vert1_c_6.z;
        float _S3754 = vert2_c_6.z;
        if(_S3752 < near_plane_14)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = _S3752 > far_plane_14;
        }
        if(_S3748)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = _S3753 < near_plane_14;
        }
        if(_S3748)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = _S3753 > far_plane_14;
        }
        if(_S3748)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = _S3754 < near_plane_14;
        }
        if(_S3748)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = _S3754 > far_plane_14;
        }
        if(_S3748)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  _S3755 = make_float2 (fx_25, fy_25);
        float2  _S3756 = make_float2 (cx_25, cy_25);
        float2  _S3757 = _S3755 * (float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S3752)) + _S3756;
        float2  _S3758 = _S3755 * (float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S3753)) + _S3756;
        float2  _S3759 = _S3755 * (float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S3754)) + _S3756;
        float2  e0_10 = _S3758 - _S3757;
        float2  e1_10 = _S3759 - _S3758;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(_S3757 - _S3759)));
        float _S3760 = _S3757.x;
        float _S3761 = _S3758.x;
        float _S3762 = _S3759.x;
        float xmax_7 = (F32_max(((F32_max((_S3760), (_S3761)))), (_S3762))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S3760), (_S3761)))), (_S3762))) - offset_4;
        float _S3763 = _S3757.y;
        float _S3764 = _S3758.y;
        float _S3765 = _S3759.y;
        float ymax_7 = (F32_max(((F32_max((_S3763), (_S3764)))), (_S3765))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S3763), (_S3764)))), (_S3765))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = xmin_7 >= float(image_width_21);
        }
        if(_S3748)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = ymax_7 <= 0.0f;
        }
        if(_S3748)
        {
            _S3748 = true;
        }
        else
        {
            _S3748 = ymin_7 >= float(image_height_21);
        }
        if(_S3748)
        {
            _S3748 = true;
        }
        else
        {
            if(_S3747 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S3748 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S3748 = false;
                }
                if(_S3748)
                {
                    _S3748 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S3748 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S3748 = false;
                    }
                }
            }
            else
            {
                _S3748 = false;
            }
        }
        if(_S3748)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S3766 = mean_25 - - mul_0(transpose_0(R_25), t_24);
        float _S3767 = _S3766.x;
        float _S3768 = _S3766.y;
        float _S3769 = _S3766.z;
        float norm_14 = (F32_sqrt((_S3767 * _S3767 + _S3768 * _S3768 + _S3769 * _S3769)));
        float x_53 = _S3767 / norm_14;
        float y_24 = _S3768 / norm_14;
        float z_21 = _S3769 / norm_14;
        float z2_48 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_53 * x_53 - y_24 * y_24;
        float fS1_21 = 2.0f * x_53 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_48 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_53) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_48 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_53) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_53 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_48 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_53) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_53 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S3770 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S3770);
        float3  _S3771 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S3772 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S3771 + _S3772 + make_float3 (0.5f), _S3770);
        (*rgbs_0)[int(2)] = max_0(_S3771 - _S3772 + make_float3 (0.5f), _S3770);
        (*verts_0)[int(0)] = vert0_2;
        (*verts_0)[int(1)] = vert1_2;
        (*verts_0)[int(2)] = vert2_2;
        float3  _S3773 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S3773 * make_float3 (float(- (F32_sign((dot_0(_S3773, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_26, float3  t_25, float fx_26, float fy_26, float cx_26, float cy_26, float4  radial_coeffs_29, float2  tangential_coeffs_29, float2  thin_prism_coeffs_29, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_26, mean_26) + t_25;
        float _S3774 = length_1(mean_c_22);
        bool _S3775;
        if(_S3774 < near_plane_15)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = _S3774 > far_plane_15;
        }
        if(_S3775)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3776 = scale_26.x;
        float sx_5 = (F32_exp((_S3776)));
        float _S3777 = scale_26.y;
        float sy_5 = (F32_exp((_S3777)));
        float sz_7 = scale_26.z - 0.5f * (_S3776 + _S3777);
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
        Matrix<float, 3, 3>  _S3778 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_49), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_49), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S3778, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S3778, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S3778, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_7 = mul_0(R_26, vert0_3) + t_25;
        float3  vert1_c_7 = mul_0(R_26, vert1_3) + t_25;
        float3  vert2_c_7 = mul_0(R_26, vert2_3) + t_25;
        float _S3779 = length_1(vert0_c_7);
        float _S3780 = length_1(vert1_c_7);
        float _S3781 = length_1(vert2_c_7);
        if(_S3779 < near_plane_15)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = _S3779 > far_plane_15;
        }
        if(_S3775)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = _S3780 < near_plane_15;
        }
        if(_S3775)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = _S3780 > far_plane_15;
        }
        if(_S3775)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = _S3781 < near_plane_15;
        }
        if(_S3775)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = _S3781 > far_plane_15;
        }
        if(_S3775)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_4 = CameraDistortion_x24init_0(radial_coeffs_29, tangential_coeffs_29, thin_prism_coeffs_29);
        float2  _S3782 = float2 {vert0_c_7.x, vert0_c_7.y};
        float r_13 = length_0(_S3782);
        float _S3783 = vert0_c_7.z;
        float theta_8 = (F32_atan2((r_13), (_S3783)));
        float k_7;
        if(theta_8 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_8 * theta_8 / 3.0f) / _S3783;
        }
        else
        {
            k_7 = theta_8 / r_13;
        }
        float2  _S3784 = _S3782 * make_float2 (k_7);
        float k1_6 = dist_coeffs_4.radial_coeffs_0.x;
        float k2_6 = dist_coeffs_4.radial_coeffs_0.y;
        float k3_6 = dist_coeffs_4.radial_coeffs_0.z;
        float k4_7 = dist_coeffs_4.radial_coeffs_0.w;
        float p1_7 = dist_coeffs_4.tangential_coeffs_0.x;
        float p2_7 = dist_coeffs_4.tangential_coeffs_0.y;
        float sx1_7 = dist_coeffs_4.thin_prism_coeffs_0.x;
        float sy1_7 = dist_coeffs_4.thin_prism_coeffs_0.y;
        float u_19 = _S3784.x;
        float v_19 = _S3784.y;
        float r2_19 = u_19 * u_19 + v_19 * v_19;
        float _S3785 = 2.0f * p1_7;
        float _S3786 = 2.0f * p2_7;
        float2  _S3787 = _S3784 * make_float2 (1.0f + r2_19 * (k1_6 + r2_19 * (k2_6 + r2_19 * (k3_6 + r2_19 * k4_7)))) + make_float2 (_S3785 * u_19 * v_19 + p2_7 * (r2_19 + 2.0f * u_19 * u_19) + sx1_7 * r2_19, _S3786 * u_19 * v_19 + p1_7 * (r2_19 + 2.0f * v_19 * v_19) + sy1_7 * r2_19);
        float _S3788 = fx_26 * _S3787.x + cx_26;
        float _S3789 = fy_26 * _S3787.y + cy_26;
        float2  _S3790 = make_float2 (_S3788, _S3789);
        float2  _S3791 = float2 {vert1_c_7.x, vert1_c_7.y};
        float r_14 = length_0(_S3791);
        float _S3792 = vert1_c_7.z;
        float theta_9 = (F32_atan2((r_14), (_S3792)));
        if(theta_9 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_9 * theta_9 / 3.0f) / _S3792;
        }
        else
        {
            k_7 = theta_9 / r_14;
        }
        float2  _S3793 = _S3791 * make_float2 (k_7);
        float u_20 = _S3793.x;
        float v_20 = _S3793.y;
        float r2_20 = u_20 * u_20 + v_20 * v_20;
        float2  _S3794 = _S3793 * make_float2 (1.0f + r2_20 * (k1_6 + r2_20 * (k2_6 + r2_20 * (k3_6 + r2_20 * k4_7)))) + make_float2 (_S3785 * u_20 * v_20 + p2_7 * (r2_20 + 2.0f * u_20 * u_20) + sx1_7 * r2_20, _S3786 * u_20 * v_20 + p1_7 * (r2_20 + 2.0f * v_20 * v_20) + sy1_7 * r2_20);
        float _S3795 = fx_26 * _S3794.x + cx_26;
        float _S3796 = fy_26 * _S3794.y + cy_26;
        float2  _S3797 = make_float2 (_S3795, _S3796);
        float2  _S3798 = float2 {vert2_c_7.x, vert2_c_7.y};
        float r_15 = length_0(_S3798);
        float _S3799 = vert2_c_7.z;
        float theta_10 = (F32_atan2((r_15), (_S3799)));
        if(theta_10 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_10 * theta_10 / 3.0f) / _S3799;
        }
        else
        {
            k_7 = theta_10 / r_15;
        }
        float2  _S3800 = _S3798 * make_float2 (k_7);
        float u_21 = _S3800.x;
        float v_21 = _S3800.y;
        float r2_21 = u_21 * u_21 + v_21 * v_21;
        float2  _S3801 = _S3800 * make_float2 (1.0f + r2_21 * (k1_6 + r2_21 * (k2_6 + r2_21 * (k3_6 + r2_21 * k4_7)))) + make_float2 (_S3785 * u_21 * v_21 + p2_7 * (r2_21 + 2.0f * u_21 * u_21) + sx1_7 * r2_21, _S3786 * u_21 * v_21 + p1_7 * (r2_21 + 2.0f * v_21 * v_21) + sy1_7 * r2_21);
        float _S3802 = fx_26 * _S3801.x + cx_26;
        float _S3803 = fy_26 * _S3801.y + cy_26;
        float2  _S3804 = make_float2 (_S3802, _S3803);
        float2  e0_11 = _S3797 - _S3790;
        float2  e1_11 = _S3804 - _S3797;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(_S3790 - _S3804)));
        float xmax_8 = (F32_max(((F32_max((_S3788), (_S3795)))), (_S3802))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S3788), (_S3795)))), (_S3802))) - offset_5;
        float ymax_8 = (F32_max(((F32_max((_S3789), (_S3796)))), (_S3803))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S3789), (_S3796)))), (_S3803))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = xmin_8 >= float(image_width_22);
        }
        if(_S3775)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = ymax_8 <= 0.0f;
        }
        if(_S3775)
        {
            _S3775 = true;
        }
        else
        {
            _S3775 = ymin_8 >= float(image_height_22);
        }
        if(_S3775)
        {
            _S3775 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S3775 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S3775 = false;
                }
                if(_S3775)
                {
                    _S3775 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S3775 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S3775 = false;
                    }
                }
            }
            else
            {
                _S3775 = false;
            }
        }
        if(_S3775)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S3805 = mean_26 - - mul_0(transpose_0(R_26), t_25);
        float _S3806 = _S3805.x;
        float _S3807 = _S3805.y;
        float _S3808 = _S3805.z;
        float norm_15 = (F32_sqrt((_S3806 * _S3806 + _S3807 * _S3807 + _S3808 * _S3808)));
        float x_55 = _S3806 / norm_15;
        float y_25 = _S3807 / norm_15;
        float z_22 = _S3808 / norm_15;
        float z2_50 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_55 * x_55 - y_25 * y_25;
        float fS1_22 = 2.0f * x_55 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_50 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_55) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_50 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_55) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_55 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_50 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_55) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_55 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S3809 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S3809);
        float3  _S3810 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S3811 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S3810 + _S3811 + make_float3 (0.5f), _S3809);
        (*rgbs_1)[int(2)] = max_0(_S3810 - _S3811 + make_float3 (0.5f), _S3809);
        (*verts_1)[int(0)] = vert0_3;
        (*verts_1)[int(1)] = vert1_3;
        (*verts_1)[int(2)] = vert2_3;
        float3  _S3812 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S3812 * make_float3 (float(- (F32_sign((dot_0(_S3812, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_27, float3  t_26, float fx_27, float fy_27, float cx_27, float cy_27, float4  radial_coeffs_30, float2  tangential_coeffs_30, float2  thin_prism_coeffs_30, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_27, mean_27) + t_26;
    float _S3813 = scale_27.x;
    float sx_6 = (F32_exp((_S3813)));
    float _S3814 = scale_27.y;
    float sy_6 = (F32_exp((_S3814)));
    float sz_8 = scale_27.z - 0.5f * (_S3813 + _S3814);
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
    Matrix<float, 3, 3>  _S3815 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_51), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_51), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
    float3  vert0_4 = mul_0(_S3815, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
    float3  vert1_4 = mul_0(_S3815, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
    float3  vert2_4 = mul_0(_S3815, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
    float3  vert0_c_8 = mul_0(R_27, vert0_4) + t_26;
    float3  vert1_c_8 = mul_0(R_27, vert1_4) + t_26;
    float3  vert2_c_8 = mul_0(R_27, vert2_4) + t_26;
    float2  _S3816 = make_float2 (fx_27, fy_27);
    float2  _S3817 = make_float2 (cx_27, cy_27);
    float2  _S3818 = _S3816 * (float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z)) + _S3817;
    float2  _S3819 = _S3816 * (float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z)) + _S3817;
    float2  _S3820 = _S3816 * (float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z)) + _S3817;
    float2  e0_12 = _S3819 - _S3818;
    float2  e1_12 = _S3820 - _S3819;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(_S3818 - _S3820)));
    float _S3821 = _S3818.x;
    float _S3822 = _S3819.x;
    float _S3823 = _S3820.x;
    float _S3824 = _S3818.y;
    float _S3825 = _S3819.y;
    float _S3826 = _S3820.y;
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S3821), (_S3822)))), (_S3823))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S3824), (_S3825)))), (_S3826))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S3821), (_S3822)))), (_S3823))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S3824), (_S3825)))), (_S3826))) + offset_6)))));
    *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
    float3  _S3827 = mean_27 - - mul_0(transpose_0(R_27), t_26);
    float _S3828 = _S3827.x;
    float _S3829 = _S3827.y;
    float _S3830 = _S3827.z;
    float norm_16 = (F32_sqrt((_S3828 * _S3828 + _S3829 * _S3829 + _S3830 * _S3830)));
    float x_57 = _S3828 / norm_16;
    float y_26 = _S3829 / norm_16;
    float z_23 = _S3830 / norm_16;
    float z2_52 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_57 * x_57 - y_26 * y_26;
    float fS1_23 = 2.0f * x_57 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_52 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_57) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_52 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_57) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_57 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_52 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_57) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_57 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S3831 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S3831);
    float3  _S3832 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S3833 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S3832 + _S3833 + make_float3 (0.5f), _S3831);
    (*rgbs_2)[int(2)] = max_0(_S3832 - _S3833 + make_float3 (0.5f), _S3831);
    (*verts_2)[int(0)] = vert0_4;
    (*verts_2)[int(1)] = vert1_4;
    (*verts_2)[int(2)] = vert2_4;
    float3  _S3834 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S3834 * make_float3 (float(- (F32_sign((dot_0(_S3834, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_28, float3  t_27, float fx_28, float fy_28, float cx_28, float cy_28, float4  radial_coeffs_31, float2  tangential_coeffs_31, float2  thin_prism_coeffs_31, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_28, mean_28) + t_27;
    float _S3835 = scale_28.x;
    float sx_7 = (F32_exp((_S3835)));
    float _S3836 = scale_28.y;
    float sy_7 = (F32_exp((_S3836)));
    float sz_9 = scale_28.z - 0.5f * (_S3835 + _S3836);
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
    Matrix<float, 3, 3>  _S3837 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_53), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_53), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S3837, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S3837, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S3837, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_9 = mul_0(R_28, vert0_5) + t_27;
    float3  vert1_c_9 = mul_0(R_28, vert1_5) + t_27;
    float3  vert2_c_9 = mul_0(R_28, vert2_5) + t_27;
    CameraDistortion_0 dist_coeffs_5 = CameraDistortion_x24init_0(radial_coeffs_31, tangential_coeffs_31, thin_prism_coeffs_31);
    float2  _S3838 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_16 = length_0(_S3838);
    float _S3839 = vert0_c_9.z;
    float theta_11 = (F32_atan2((r_16), (_S3839)));
    float k_8;
    if(theta_11 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_11 * theta_11 / 3.0f) / _S3839;
    }
    else
    {
        k_8 = theta_11 / r_16;
    }
    float2  _S3840 = _S3838 * make_float2 (k_8);
    float k1_7 = dist_coeffs_5.radial_coeffs_0.x;
    float k2_7 = dist_coeffs_5.radial_coeffs_0.y;
    float k3_7 = dist_coeffs_5.radial_coeffs_0.z;
    float k4_8 = dist_coeffs_5.radial_coeffs_0.w;
    float p1_8 = dist_coeffs_5.tangential_coeffs_0.x;
    float p2_8 = dist_coeffs_5.tangential_coeffs_0.y;
    float sx1_8 = dist_coeffs_5.thin_prism_coeffs_0.x;
    float sy1_8 = dist_coeffs_5.thin_prism_coeffs_0.y;
    float u_22 = _S3840.x;
    float v_22 = _S3840.y;
    float r2_22 = u_22 * u_22 + v_22 * v_22;
    float _S3841 = 2.0f * p1_8;
    float _S3842 = 2.0f * p2_8;
    float2  _S3843 = _S3840 * make_float2 (1.0f + r2_22 * (k1_7 + r2_22 * (k2_7 + r2_22 * (k3_7 + r2_22 * k4_8)))) + make_float2 (_S3841 * u_22 * v_22 + p2_8 * (r2_22 + 2.0f * u_22 * u_22) + sx1_8 * r2_22, _S3842 * u_22 * v_22 + p1_8 * (r2_22 + 2.0f * v_22 * v_22) + sy1_8 * r2_22);
    float _S3844 = fx_28 * _S3843.x + cx_28;
    float _S3845 = fy_28 * _S3843.y + cy_28;
    float2  _S3846 = make_float2 (_S3844, _S3845);
    float2  _S3847 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_17 = length_0(_S3847);
    float _S3848 = vert1_c_9.z;
    float theta_12 = (F32_atan2((r_17), (_S3848)));
    if(theta_12 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_12 * theta_12 / 3.0f) / _S3848;
    }
    else
    {
        k_8 = theta_12 / r_17;
    }
    float2  _S3849 = _S3847 * make_float2 (k_8);
    float u_23 = _S3849.x;
    float v_23 = _S3849.y;
    float r2_23 = u_23 * u_23 + v_23 * v_23;
    float2  _S3850 = _S3849 * make_float2 (1.0f + r2_23 * (k1_7 + r2_23 * (k2_7 + r2_23 * (k3_7 + r2_23 * k4_8)))) + make_float2 (_S3841 * u_23 * v_23 + p2_8 * (r2_23 + 2.0f * u_23 * u_23) + sx1_8 * r2_23, _S3842 * u_23 * v_23 + p1_8 * (r2_23 + 2.0f * v_23 * v_23) + sy1_8 * r2_23);
    float _S3851 = fx_28 * _S3850.x + cx_28;
    float _S3852 = fy_28 * _S3850.y + cy_28;
    float2  _S3853 = make_float2 (_S3851, _S3852);
    float2  _S3854 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_18 = length_0(_S3854);
    float _S3855 = vert2_c_9.z;
    float theta_13 = (F32_atan2((r_18), (_S3855)));
    if(theta_13 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_13 * theta_13 / 3.0f) / _S3855;
    }
    else
    {
        k_8 = theta_13 / r_18;
    }
    float2  _S3856 = _S3854 * make_float2 (k_8);
    float u_24 = _S3856.x;
    float v_24 = _S3856.y;
    float r2_24 = u_24 * u_24 + v_24 * v_24;
    float2  _S3857 = _S3856 * make_float2 (1.0f + r2_24 * (k1_7 + r2_24 * (k2_7 + r2_24 * (k3_7 + r2_24 * k4_8)))) + make_float2 (_S3841 * u_24 * v_24 + p2_8 * (r2_24 + 2.0f * u_24 * u_24) + sx1_8 * r2_24, _S3842 * u_24 * v_24 + p1_8 * (r2_24 + 2.0f * v_24 * v_24) + sy1_8 * r2_24);
    float _S3858 = fx_28 * _S3857.x + cx_28;
    float _S3859 = fy_28 * _S3857.y + cy_28;
    float2  _S3860 = make_float2 (_S3858, _S3859);
    float2  e0_13 = _S3853 - _S3846;
    float2  e1_13 = _S3860 - _S3853;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(_S3846 - _S3860)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S3844), (_S3851)))), (_S3858))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S3845), (_S3852)))), (_S3859))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S3844), (_S3851)))), (_S3858))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S3845), (_S3852)))), (_S3859))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S3861 = mean_28 - - mul_0(transpose_0(R_28), t_27);
    float _S3862 = _S3861.x;
    float _S3863 = _S3861.y;
    float _S3864 = _S3861.z;
    float norm_17 = (F32_sqrt((_S3862 * _S3862 + _S3863 * _S3863 + _S3864 * _S3864)));
    float x_59 = _S3862 / norm_17;
    float y_27 = _S3863 / norm_17;
    float z_24 = _S3864 / norm_17;
    float z2_54 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_59 * x_59 - y_27 * y_27;
    float fS1_24 = 2.0f * x_59 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_54 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_59) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_54 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_59) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_59 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_54 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_59) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_59 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S3865 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S3865);
    float3  _S3866 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S3867 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S3866 + _S3867 + make_float3 (0.5f), _S3865);
    (*rgbs_3)[int(2)] = max_0(_S3866 - _S3867 + make_float3 (0.5f), _S3865);
    (*verts_3)[int(0)] = vert0_5;
    (*verts_3)[int(1)] = vert1_5;
    (*verts_3)[int(2)] = vert2_5;
    float3  _S3868 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S3868 * make_float3 (float(- (F32_sign((dot_0(_S3868, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_29, float3  t_28, float fx_29, float fy_29, float cx_29, float cy_29, float4  radial_coeffs_32, float2  tangential_coeffs_32, float2  thin_prism_coeffs_32, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_1(R_29, mean_29) + t_28;
    float _S3869 = scale_29.x;
    float _S3870 = s_primal_ctx_exp_1(_S3869);
    float _S3871 = scale_29.y;
    float _S3872 = s_primal_ctx_exp_1(_S3871);
    float sz_10 = scale_29.z - 0.5f * (_S3869 + _S3871);
    float _S3873 = quat_30.y;
    float x2_30 = _S3873 * _S3873;
    float y2_30 = quat_30.z * quat_30.z;
    float z2_55 = quat_30.w * quat_30.w;
    float xy_30 = quat_30.y * quat_30.z;
    float xz_30 = quat_30.y * quat_30.w;
    float yz_30 = quat_30.z * quat_30.w;
    float wx_30 = quat_30.x * quat_30.y;
    float wy_30 = quat_30.x * quat_30.z;
    float wz_30 = quat_30.x * quat_30.w;
    Matrix<float, 3, 3>  _S3874 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_55), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_55), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  _S3875 = make_float3 (_S3870, 0.0f, 0.0f);
    float3  vert0_6 = s_primal_ctx_mul_1(_S3874, _S3875) + mean_29;
    float _S3876 = -0.5f + sz_10;
    float3  _S3877 = make_float3 (_S3870 * _S3876, _S3872, 0.0f);
    float3  vert1_6 = s_primal_ctx_mul_1(_S3874, _S3877) + mean_29;
    float _S3878 = -0.5f - sz_10;
    float3  _S3879 = make_float3 (_S3870 * _S3878, - _S3872, 0.0f);
    float3  vert2_6 = s_primal_ctx_mul_1(_S3874, _S3879) + mean_29;
    float3  vert0_c_10 = s_primal_ctx_mul_1(R_29, vert0_6) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_1(R_29, vert1_6) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_1(R_29, vert2_6) + t_28;
    float2  _S3880 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S3881 = vert0_c_10.z;
    float2  _S3882 = make_float2 (_S3881);
    float2  _S3883 = make_float2 (_S3881 * _S3881);
    float2  _S3884 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S3885 = vert1_c_10.z;
    float2  _S3886 = make_float2 (_S3885);
    float2  _S3887 = make_float2 (_S3885 * _S3885);
    float2  _S3888 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S3889 = vert2_c_10.z;
    float2  _S3890 = make_float2 (_S3889);
    float2  _S3891 = make_float2 (_S3889 * _S3889);
    float2  _S3892 = make_float2 (fx_29, fy_29);
    float2  _S3893 = make_float2 (cx_29, cy_29);
    float2  _S3894 = _S3892 * (_S3880 / make_float2 (_S3881)) + _S3893;
    float2  _S3895 = _S3892 * (_S3884 / make_float2 (_S3885)) + _S3893;
    float2  _S3896 = _S3892 * (_S3888 / make_float2 (_S3889)) + _S3893;
    float2  e0_14 = _S3895 - _S3894;
    float2  e1_14 = _S3896 - _S3895;
    float2  e2_6 = _S3894 - _S3896;
    float _S3897 = e0_14.x;
    float _S3898 = e1_14.y;
    float _S3899 = e0_14.y;
    float _S3900 = e1_14.x;
    float _S3901 = _S3897 * _S3898 - _S3899 * _S3900;
    float _S3902 = 1.0f - hardness_14.y;
    float _S3903 = -1.0f / _S3902;
    float _S3904 = _S3902 * _S3902;
    float _S3905 = _S3894.x;
    float _S3906 = _S3895.x;
    float _S3907 = s_primal_ctx_max_0(_S3905, _S3906);
    float _S3908 = _S3896.x;
    float _S3909 = s_primal_ctx_min_0(_S3905, _S3906);
    float _S3910 = _S3894.y;
    float _S3911 = _S3895.y;
    float _S3912 = s_primal_ctx_max_0(_S3910, _S3911);
    float _S3913 = _S3896.y;
    float _S3914 = s_primal_ctx_min_0(_S3910, _S3911);
    float3  _S3915 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S3916 = length_1(_S3915) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S3917 = transpose_0(R_29);
    float3  _S3918 = mean_29 - - s_primal_ctx_mul_1(_S3917, t_28);
    float _S3919 = _S3918.x;
    float _S3920 = _S3918.y;
    float _S3921 = _S3918.z;
    float _S3922 = _S3919 * _S3919 + _S3920 * _S3920 + _S3921 * _S3921;
    float _S3923 = s_primal_ctx_sqrt_0(_S3922);
    float x_60 = _S3919 / _S3923;
    float3  _S3924 = make_float3 (x_60);
    float _S3925 = _S3923 * _S3923;
    float y_28 = _S3920 / _S3923;
    float z_25 = _S3921 / _S3923;
    float3  _S3926 = make_float3 (z_25);
    float _S3927 = - y_28;
    float3  _S3928 = make_float3 (_S3927);
    float z2_56 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_60 * x_60 - y_28 * y_28;
    float _S3929 = 2.0f * x_60;
    float fS1_25 = _S3929 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_56 - 0.31539157032966614f;
    float3  _S3930 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_60;
    float3  _S3931 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S3932 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S3933 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S3934 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_56 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S3935 = 1.86588168144226074f * z2_56 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S3935;
    float3  _S3936 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_60;
    float3  _S3937 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S3938 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S3939 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S3940 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_60 * fC1_25 - y_28 * fS1_25);
    float3  _S3941 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_60 * fS1_25 + y_28 * fC1_25);
    float3  _S3942 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3927) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_60) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S3943 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S3944 = make_float3 (0.0f);
    float3  _S3945 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S3946 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3947 = make_float3 (_S3946);
    float3  _S3948 = (*ch_coeffs_10)[int(1)] * make_float3 (_S3946);
    float3  _S3949 = _S3945 + _S3948 + make_float3 (0.5f);
    float3  _S3950 = _S3945 - _S3948 + make_float3 (0.5f);
    float3  _S3951 = vert1_c_10 - vert0_c_10;
    float3  _S3952 = vert2_c_10 - vert0_c_10;
    float3  _S3953 = s_primal_ctx_cross_0(_S3951, _S3952);
    float3  _S3954 = normalize_0(_S3953);
    float3  _S3955 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3954, mean_c_25)))))) * v_normal_2;
    float3  _S3956 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3957;
    (&_S3957)->primal_0 = _S3954;
    (&_S3957)->differential_0 = _S3956;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3958;
    (&_S3958)->primal_0 = mean_c_25;
    (&_S3958)->differential_0 = _S3956;
    s_bwd_prop_dot_0(&_S3957, &_S3958, 0.0f);
    float3  _S3959 = _S3955 + _S3957.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3960;
    (&_S3960)->primal_0 = _S3953;
    (&_S3960)->differential_0 = _S3956;
    s_bwd_normalize_impl_0(&_S3960, _S3959);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3961;
    (&_S3961)->primal_0 = _S3951;
    (&_S3961)->differential_0 = _S3956;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3962;
    (&_S3962)->primal_0 = _S3952;
    (&_S3962)->differential_0 = _S3956;
    s_bwd_prop_cross_0(&_S3961, &_S3962, _S3960.differential_0);
    float3  _S3963 = - _S3962.differential_0;
    float3  _S3964 = - _S3961.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3965;
    (&_S3965)->primal_0 = _S3950;
    (&_S3965)->differential_0 = _S3956;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3966;
    (&_S3966)->primal_0 = _S3944;
    (&_S3966)->differential_0 = _S3956;
    s_bwd_prop_max_0(&_S3965, &_S3966, (*v_rgbs_0)[int(2)]);
    float3  _S3967 = - _S3965.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3968;
    (&_S3968)->primal_0 = _S3949;
    (&_S3968)->differential_0 = _S3956;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3969;
    (&_S3969)->primal_0 = _S3944;
    (&_S3969)->differential_0 = _S3956;
    s_bwd_prop_max_0(&_S3968, &_S3969, (*v_rgbs_0)[int(1)]);
    float3  _S3970 = _S3947 * (_S3967 + _S3968.differential_0);
    float3  _S3971 = _S3965.differential_0 + _S3968.differential_0;
    float3  _S3972 = make_float3 (0.5f) * - _S3971;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3973;
    (&_S3973)->primal_0 = _S3943;
    (&_S3973)->differential_0 = _S3956;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3974;
    (&_S3974)->primal_0 = _S3944;
    (&_S3974)->differential_0 = _S3956;
    s_bwd_prop_max_0(&_S3973, &_S3974, (*v_rgbs_0)[int(0)]);
    float3  _S3975 = _S3972 + _S3973.differential_0;
    float3  _S3976 = _S3971 + _S3973.differential_0;
    float3  _S3977 = _S3941 * _S3976;
    float3  _S3978 = (*sh_coeffs_25)[int(15)] * _S3976;
    float3  _S3979 = _S3939 * _S3976;
    float3  _S3980 = (*sh_coeffs_25)[int(14)] * _S3976;
    float3  _S3981 = _S3937 * _S3976;
    float3  _S3982 = (*sh_coeffs_25)[int(13)] * _S3976;
    float3  _S3983 = _S3936 * _S3976;
    float3  _S3984 = (*sh_coeffs_25)[int(12)] * _S3976;
    float3  _S3985 = _S3938 * _S3976;
    float3  _S3986 = (*sh_coeffs_25)[int(11)] * _S3976;
    float3  _S3987 = _S3940 * _S3976;
    float3  _S3988 = (*sh_coeffs_25)[int(10)] * _S3976;
    float3  _S3989 = _S3942 * _S3976;
    float3  _S3990 = (*sh_coeffs_25)[int(9)] * _S3976;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S3990.x + _S3990.y + _S3990.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S3978.x + _S3978.y + _S3978.z);
    float _S3991 = _S3988.x + _S3988.y + _S3988.z;
    float _S3992 = _S3980.x + _S3980.y + _S3980.z;
    float _S3993 = _S3986.x + _S3986.y + _S3986.z;
    float _S3994 = _S3982.x + _S3982.y + _S3982.z;
    float _S3995 = _S3984.x + _S3984.y + _S3984.z;
    float _S3996 = - s_diff_fC2_T_7;
    float3  _S3997 = _S3933 * _S3976;
    float3  _S3998 = (*sh_coeffs_25)[int(8)] * _S3976;
    float3  _S3999 = _S3931 * _S3976;
    float3  _S4000 = (*sh_coeffs_25)[int(7)] * _S3976;
    float3  _S4001 = _S3930 * _S3976;
    float3  _S4002 = (*sh_coeffs_25)[int(6)] * _S3976;
    float3  _S4003 = _S3932 * _S3976;
    float3  _S4004 = (*sh_coeffs_25)[int(5)] * _S3976;
    float3  _S4005 = _S3934 * _S3976;
    float3  _S4006 = (*sh_coeffs_25)[int(4)] * _S3976;
    float _S4007 = _S4004.x + _S4004.y + _S4004.z;
    float _S4008 = _S4000.x + _S4000.y + _S4000.z;
    float _S4009 = fTmp1B_25 * _S3991 + x_60 * s_diff_fS2_T_7 + y_28 * _S3996 + 0.54627424478530884f * (_S4006.x + _S4006.y + _S4006.z);
    float _S4010 = fTmp1B_25 * _S3992 + y_28 * s_diff_fS2_T_7 + x_60 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S3998.x + _S3998.y + _S3998.z);
    float _S4011 = y_28 * - _S4010;
    float _S4012 = x_60 * _S4010;
    float _S4013 = z_25 * (1.86588168144226074f * (z_25 * _S3995) + -2.28522896766662598f * (y_28 * _S3993 + x_60 * _S3994) + 0.94617468118667603f * (_S4002.x + _S4002.y + _S4002.z));
    float3  _S4014 = make_float3 (0.48860251903533936f) * _S3976;
    float3  _S4015 = - _S4014;
    float3  _S4016 = _S3924 * _S4015;
    float3  _S4017 = (*sh_coeffs_25)[int(3)] * _S4015;
    float3  _S4018 = _S3926 * _S4014;
    float3  _S4019 = (*sh_coeffs_25)[int(2)] * _S4014;
    float3  _S4020 = _S3928 * _S4014;
    float3  _S4021 = (*sh_coeffs_25)[int(1)] * _S4014;
    float _S4022 = (_S3935 * _S3995 + 1.44530570507049561f * (fS1_25 * _S3991 + fC1_25 * _S3992) + -1.09254848957061768f * (y_28 * _S4007 + x_60 * _S4008) + _S4013 + _S4013 + _S4019.x + _S4019.y + _S4019.z) / _S3925;
    float _S4023 = _S3923 * _S4022;
    float _S4024 = (fTmp0C_25 * _S3993 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S3996 + fTmp0B_25 * _S4007 + _S3929 * _S4009 + _S4011 + _S4011 + - (_S4021.x + _S4021.y + _S4021.z)) / _S3925;
    float _S4025 = _S3923 * _S4024;
    float _S4026 = (fTmp0C_25 * _S3994 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4008 + 2.0f * (y_28 * _S4009) + _S4012 + _S4012 + _S4017.x + _S4017.y + _S4017.z) / _S3925;
    float _S4027 = _S3923 * _S4026;
    float _S4028 = _S3921 * - _S4022 + _S3920 * - _S4024 + _S3919 * - _S4026;
    DiffPair_float_0 _S4029;
    (&_S4029)->primal_0 = _S3922;
    (&_S4029)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4029, _S4028);
    float _S4030 = _S3921 * _S4029.differential_0;
    float _S4031 = _S3920 * _S4029.differential_0;
    float _S4032 = _S3919 * _S4029.differential_0;
    float3  _S4033 = make_float3 (0.282094806432724f) * _S3976;
    float3  _S4034 = make_float3 (_S4027 + _S4032 + _S4032, _S4025 + _S4031 + _S4031, _S4023 + _S4030 + _S4030);
    float3  _S4035 = - - _S4034;
    Matrix<float, 3, 3>  _S4036 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4037;
    (&_S4037)->primal_0 = _S3917;
    (&_S4037)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4038;
    (&_S4038)->primal_0 = t_28;
    (&_S4038)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4037, &_S4038, _S4035);
    Matrix<float, 3, 3>  _S4039 = transpose_0(_S4037.differential_0);
    DiffPair_float_0 _S4040;
    (&_S4040)->primal_0 = _S3916;
    (&_S4040)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4040, v_depth_9);
    float _S4041 = 0.3333333432674408f * _S4040.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4042;
    (&_S4042)->primal_0 = _S3915;
    (&_S4042)->differential_0 = _S3956;
    s_bwd_length_impl_0(&_S4042, _S4041);
    DiffPair_float_0 _S4043;
    (&_S4043)->primal_0 = _S3914;
    (&_S4043)->differential_0 = 0.0f;
    DiffPair_float_0 _S4044;
    (&_S4044)->primal_0 = _S3913;
    (&_S4044)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4043, &_S4044, 0.0f);
    DiffPair_float_0 _S4045;
    (&_S4045)->primal_0 = _S3910;
    (&_S4045)->differential_0 = 0.0f;
    DiffPair_float_0 _S4046;
    (&_S4046)->primal_0 = _S3911;
    (&_S4046)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4045, &_S4046, _S4043.differential_0);
    DiffPair_float_0 _S4047;
    (&_S4047)->primal_0 = _S3912;
    (&_S4047)->differential_0 = 0.0f;
    DiffPair_float_0 _S4048;
    (&_S4048)->primal_0 = _S3913;
    (&_S4048)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4047, &_S4048, 0.0f);
    float _S4049 = _S4044.differential_0 + _S4048.differential_0;
    DiffPair_float_0 _S4050;
    (&_S4050)->primal_0 = _S3910;
    (&_S4050)->differential_0 = 0.0f;
    DiffPair_float_0 _S4051;
    (&_S4051)->primal_0 = _S3911;
    (&_S4051)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4050, &_S4051, _S4047.differential_0);
    float _S4052 = _S4046.differential_0 + _S4051.differential_0;
    float _S4053 = _S4045.differential_0 + _S4050.differential_0;
    DiffPair_float_0 _S4054;
    (&_S4054)->primal_0 = _S3909;
    (&_S4054)->differential_0 = 0.0f;
    DiffPair_float_0 _S4055;
    (&_S4055)->primal_0 = _S3908;
    (&_S4055)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4054, &_S4055, 0.0f);
    DiffPair_float_0 _S4056;
    (&_S4056)->primal_0 = _S3905;
    (&_S4056)->differential_0 = 0.0f;
    DiffPair_float_0 _S4057;
    (&_S4057)->primal_0 = _S3906;
    (&_S4057)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4056, &_S4057, _S4054.differential_0);
    DiffPair_float_0 _S4058;
    (&_S4058)->primal_0 = _S3907;
    (&_S4058)->differential_0 = 0.0f;
    DiffPair_float_0 _S4059;
    (&_S4059)->primal_0 = _S3908;
    (&_S4059)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4058, &_S4059, 0.0f);
    float _S4060 = _S4055.differential_0 + _S4059.differential_0;
    DiffPair_float_0 _S4061;
    (&_S4061)->primal_0 = _S3905;
    (&_S4061)->differential_0 = 0.0f;
    DiffPair_float_0 _S4062;
    (&_S4062)->primal_0 = _S3906;
    (&_S4062)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4061, &_S4062, _S4058.differential_0);
    float _S4063 = _S4057.differential_0 + _S4062.differential_0;
    float _S4064 = _S4056.differential_0 + _S4061.differential_0;
    DiffPair_float_0 _S4065;
    (&_S4065)->primal_0 = _S3903;
    (&_S4065)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4065, 0.0f);
    float _S4066 = - (-1.0f * - (_S4065.differential_0 / _S3904));
    float2  _S4067 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4068;
    (&_S4068)->primal_0 = e2_6;
    (&_S4068)->differential_0 = _S4067;
    s_bwd_length_impl_1(&_S4068, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4069;
    (&_S4069)->primal_0 = e1_14;
    (&_S4069)->differential_0 = _S4067;
    s_bwd_length_impl_1(&_S4069, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4070;
    (&_S4070)->primal_0 = e0_14;
    (&_S4070)->differential_0 = _S4067;
    s_bwd_length_impl_1(&_S4070, -0.0f);
    DiffPair_float_0 _S4071;
    (&_S4071)->primal_0 = _S3901;
    (&_S4071)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4071, 0.0f);
    float _S4072 = - _S4071.differential_0;
    float2  _S4073 = _S4069.differential_0 + make_float2 (_S3899 * _S4072, _S3897 * _S4071.differential_0);
    float2  _S4074 = _S4070.differential_0 + make_float2 (_S3898 * _S4071.differential_0, _S3900 * _S4072);
    float2  _S4075 = _S3892 * (- _S4068.differential_0 + _S4073 + make_float2 (_S4060, _S4049)) / _S3891;
    float2  _S4076 = _S3888 * - _S4075;
    float2  _S4077 = _S3890 * _S4075;
    float2  _S4078 = _S3892 * (- _S4073 + _S4074 + make_float2 (_S4063, _S4052)) / _S3887;
    float2  _S4079 = _S3884 * - _S4078;
    float2  _S4080 = _S3886 * _S4078;
    float _S4081 = _S4079.x + _S4079.y;
    float2  _S4082 = _S3892 * (_S4068.differential_0 + - _S4074 + make_float2 (_S4064, _S4053)) / _S3883;
    float2  _S4083 = _S3880 * - _S4082;
    float2  _S4084 = _S3882 * _S4082;
    float _S4085 = _S4083.x + _S4083.y;
    float3  _S4086 = _S3962.differential_0 + _S4042.differential_0 + make_float3 (_S4077.x, _S4077.y, _S4076.x + _S4076.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4087;
    (&_S4087)->primal_0 = R_29;
    (&_S4087)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4088;
    (&_S4088)->primal_0 = vert2_6;
    (&_S4088)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4087, &_S4088, _S4086);
    float3  _S4089 = _S3961.differential_0 + _S4042.differential_0 + make_float3 (_S4080.x, _S4080.y, _S4081);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4090;
    (&_S4090)->primal_0 = R_29;
    (&_S4090)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4091;
    (&_S4091)->primal_0 = vert1_6;
    (&_S4091)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4090, &_S4091, _S4089);
    float3  _S4092 = _S3963 + _S3964 + _S4042.differential_0 + make_float3 (_S4084.x, _S4084.y, _S4085);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4093;
    (&_S4093)->primal_0 = R_29;
    (&_S4093)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4094;
    (&_S4094)->primal_0 = vert0_6;
    (&_S4094)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4093, &_S4094, _S4092);
    float3  _S4095 = (*v_verts_0)[int(2)] + _S4088.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4096;
    (&_S4096)->primal_0 = _S3874;
    (&_S4096)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4097;
    (&_S4097)->primal_0 = _S3879;
    (&_S4097)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4096, &_S4097, _S4095);
    float _S4098 = - _S4097.differential_0.y;
    float _S4099 = _S3878 * _S4097.differential_0.x;
    float _S4100 = - (_S3870 * _S4097.differential_0.x);
    float3  _S4101 = (*v_verts_0)[int(1)] + _S4091.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4102;
    (&_S4102)->primal_0 = _S3874;
    (&_S4102)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4103;
    (&_S4103)->primal_0 = _S3877;
    (&_S4103)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4102, &_S4103, _S4101);
    float _S4104 = _S3870 * _S4103.differential_0.x;
    float _S4105 = _S3876 * _S4103.differential_0.x;
    float3  _S4106 = (*v_verts_0)[int(0)] + _S4094.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4107;
    (&_S4107)->primal_0 = _S3874;
    (&_S4107)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4108;
    (&_S4108)->primal_0 = _S3875;
    (&_S4108)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4107, &_S4108, _S4106);
    Matrix<float, 3, 3>  _S4109 = transpose_0(_S4096.differential_0 + _S4102.differential_0 + _S4107.differential_0);
    float _S4110 = 2.0f * - _S4109.rows[int(2)].z;
    float _S4111 = 2.0f * _S4109.rows[int(2)].y;
    float _S4112 = 2.0f * _S4109.rows[int(2)].x;
    float _S4113 = 2.0f * _S4109.rows[int(1)].z;
    float _S4114 = 2.0f * - _S4109.rows[int(1)].y;
    float _S4115 = 2.0f * _S4109.rows[int(1)].x;
    float _S4116 = 2.0f * _S4109.rows[int(0)].z;
    float _S4117 = 2.0f * _S4109.rows[int(0)].y;
    float _S4118 = 2.0f * - _S4109.rows[int(0)].x;
    float _S4119 = - _S4115 + _S4117;
    float _S4120 = _S4112 + - _S4116;
    float _S4121 = - _S4111 + _S4113;
    float _S4122 = _S4111 + _S4113;
    float _S4123 = _S4112 + _S4116;
    float _S4124 = _S4115 + _S4117;
    float _S4125 = quat_30.w * (_S4114 + _S4118);
    float _S4126 = quat_30.z * (_S4110 + _S4118);
    float _S4127 = quat_30.y * (_S4110 + _S4114);
    float _S4128 = quat_30.x * _S4119 + quat_30.z * _S4122 + quat_30.y * _S4123 + _S4125 + _S4125;
    float _S4129 = quat_30.x * _S4120 + quat_30.w * _S4122 + quat_30.y * _S4124 + _S4126 + _S4126;
    float _S4130 = quat_30.x * _S4121 + quat_30.w * _S4123 + quat_30.z * _S4124 + _S4127 + _S4127;
    float _S4131 = quat_30.w * _S4119 + quat_30.z * _S4120 + quat_30.y * _S4121;
    float _S4132 = _S4100 + _S4104;
    float _S4133 = 0.5f * - _S4132;
    float _S4134 = _S4098 + _S4103.differential_0.y;
    DiffPair_float_0 _S4135;
    (&_S4135)->primal_0 = _S3871;
    (&_S4135)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4135, _S4134);
    float _S4136 = _S4133 + _S4135.differential_0;
    float _S4137 = _S4099 + _S4105 + _S4108.differential_0.x;
    DiffPair_float_0 _S4138;
    (&_S4138)->primal_0 = _S3869;
    (&_S4138)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4138, _S4137);
    float _S4139 = _S4133 + _S4138.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4140;
    (&_S4140)->primal_0 = R_29;
    (&_S4140)->differential_0 = _S4036;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4141;
    (&_S4141)->primal_0 = mean_29;
    (&_S4141)->differential_0 = _S3956;
    s_bwd_prop_mul_1(&_S4140, &_S4141, _S3958.differential_0);
    float3  _S4142 = _S4038.differential_0 + _S4086 + _S4089 + _S4092 + _S3958.differential_0;
    Matrix<float, 3, 3>  _S4143 = _S4039 + _S4087.differential_0 + _S4090.differential_0 + _S4093.differential_0 + _S4140.differential_0;
    FixedArray<float3 , 2>  _S4144;
    _S4144[int(0)] = _S3956;
    _S4144[int(1)] = _S3956;
    _S4144[int(1)] = _S3970;
    _S4144[int(0)] = _S3975;
    FixedArray<float3 , 16>  _S4145;
    _S4145[int(0)] = _S3956;
    _S4145[int(1)] = _S3956;
    _S4145[int(2)] = _S3956;
    _S4145[int(3)] = _S3956;
    _S4145[int(4)] = _S3956;
    _S4145[int(5)] = _S3956;
    _S4145[int(6)] = _S3956;
    _S4145[int(7)] = _S3956;
    _S4145[int(8)] = _S3956;
    _S4145[int(9)] = _S3956;
    _S4145[int(10)] = _S3956;
    _S4145[int(11)] = _S3956;
    _S4145[int(12)] = _S3956;
    _S4145[int(13)] = _S3956;
    _S4145[int(14)] = _S3956;
    _S4145[int(15)] = _S3956;
    _S4145[int(15)] = _S3977;
    _S4145[int(14)] = _S3979;
    _S4145[int(13)] = _S3981;
    _S4145[int(12)] = _S3983;
    _S4145[int(11)] = _S3985;
    _S4145[int(10)] = _S3987;
    _S4145[int(9)] = _S3989;
    _S4145[int(8)] = _S3997;
    _S4145[int(7)] = _S3999;
    _S4145[int(6)] = _S4001;
    _S4145[int(5)] = _S4003;
    _S4145[int(4)] = _S4005;
    _S4145[int(3)] = _S4016;
    _S4145[int(2)] = _S4018;
    _S4145[int(1)] = _S4020;
    _S4145[int(0)] = _S4033;
    float2  _S4146 = make_float2 (0.0f, _S4066);
    float3  _S4147 = make_float3 (_S4139, _S4136, _S4132);
    float4  _S4148 = make_float4 (0.0f);
    *&((&_S4148)->w) = _S4128;
    *&((&_S4148)->z) = _S4129;
    *&((&_S4148)->y) = _S4130;
    *&((&_S4148)->x) = _S4131;
    *v_mean_9 = _S4034 + _S4095 + _S4101 + _S4106 + _S4141.differential_0;
    *v_quat_8 = _S4148;
    *v_scale_8 = _S4147;
    *v_hardness_4 = _S4146;
    *v_sh_coeffs_7 = _S4145;
    *v_ch_coeffs_2 = _S4144;
    *v_R_8 = _S4143;
    *v_t_8 = _S4142;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_30, float3  t_29, float fx_30, float fy_30, float cx_30, float cy_30, float4  radial_coeffs_33, float2  tangential_coeffs_33, float2  thin_prism_coeffs_33, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_1(R_30, mean_30) + t_29;
    float _S4149 = scale_30.x;
    float _S4150 = s_primal_ctx_exp_1(_S4149);
    float _S4151 = scale_30.y;
    float _S4152 = s_primal_ctx_exp_1(_S4151);
    float sz_11 = scale_30.z - 0.5f * (_S4149 + _S4151);
    float _S4153 = quat_31.y;
    float x2_31 = _S4153 * _S4153;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_57 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4154 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_57), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_57), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4155 = make_float3 (_S4150, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_1(_S4154, _S4155) + mean_30;
    float _S4156 = -0.5f + sz_11;
    float3  _S4157 = make_float3 (_S4150 * _S4156, _S4152, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_1(_S4154, _S4157) + mean_30;
    float _S4158 = -0.5f - sz_11;
    float3  _S4159 = make_float3 (_S4150 * _S4158, - _S4152, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_1(_S4154, _S4159) + mean_30;
    float3  vert0_c_11 = s_primal_ctx_mul_1(R_30, vert0_7) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_1(R_30, vert1_7) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_1(R_30, vert2_7) + t_29;
    CameraDistortion_0 _S4160 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_33, tangential_coeffs_33, thin_prism_coeffs_33);
    float2  _S4161 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4162 = length_0(_S4161);
    float _S4163 = vert0_c_11.z;
    float _S4164 = s_primal_ctx_atan2_0(_S4162, _S4163);
    bool _S4165 = _S4164 < 0.00100000004749745f;
    float k_9;
    float _S4166;
    float _S4167;
    float _S4168;
    if(_S4165)
    {
        float _S4169 = 1.0f - _S4164 * _S4164 / 3.0f;
        float _S4170 = _S4163 * _S4163;
        k_9 = _S4169 / _S4163;
        _S4166 = 0.0f;
        _S4167 = _S4170;
        _S4168 = _S4169;
    }
    else
    {
        float _S4171 = _S4162 * _S4162;
        k_9 = _S4164 / _S4162;
        _S4166 = _S4171;
        _S4167 = 0.0f;
        _S4168 = 0.0f;
    }
    float2  _S4172 = make_float2 (k_9);
    float2  _S4173 = _S4161 * make_float2 (k_9);
    float k1_8 = _S4160.radial_coeffs_0.x;
    float k2_8 = _S4160.radial_coeffs_0.y;
    float k3_8 = _S4160.radial_coeffs_0.z;
    float k4_9 = _S4160.radial_coeffs_0.w;
    float p1_9 = _S4160.tangential_coeffs_0.x;
    float p2_9 = _S4160.tangential_coeffs_0.y;
    float sx1_9 = _S4160.thin_prism_coeffs_0.x;
    float sy1_9 = _S4160.thin_prism_coeffs_0.y;
    float u_25 = _S4173.x;
    float v_25 = _S4173.y;
    float r2_25 = u_25 * u_25 + v_25 * v_25;
    float _S4174 = k3_8 + r2_25 * k4_9;
    float _S4175 = k2_8 + r2_25 * _S4174;
    float _S4176 = k1_8 + r2_25 * _S4175;
    float radial_5 = 1.0f + r2_25 * _S4176;
    float2  _S4177 = make_float2 (radial_5);
    float _S4178 = 2.0f * p1_9;
    float _S4179 = _S4178 * u_25;
    float _S4180 = 2.0f * u_25;
    float _S4181 = r2_25 + _S4180 * u_25;
    float _S4182 = 2.0f * p2_9;
    float _S4183 = _S4182 * u_25;
    float _S4184 = 2.0f * v_25;
    float _S4185 = r2_25 + _S4184 * v_25;
    float2  _S4186 = _S4173 * make_float2 (radial_5) + make_float2 (_S4179 * v_25 + p2_9 * _S4181 + sx1_9 * r2_25, _S4183 * v_25 + p1_9 * _S4185 + sy1_9 * r2_25);
    float _S4187 = fx_30 * _S4186.x + cx_30;
    float _S4188 = fy_30 * _S4186.y + cy_30;
    float2  _S4189 = make_float2 (_S4187, _S4188);
    float2  _S4190 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4191 = length_0(_S4190);
    float _S4192 = vert1_c_11.z;
    float _S4193 = s_primal_ctx_atan2_0(_S4191, _S4192);
    bool _S4194 = _S4193 < 0.00100000004749745f;
    float _S4195;
    float _S4196;
    float _S4197;
    if(_S4194)
    {
        float _S4198 = 1.0f - _S4193 * _S4193 / 3.0f;
        float _S4199 = _S4192 * _S4192;
        k_9 = _S4198 / _S4192;
        _S4195 = 0.0f;
        _S4196 = _S4199;
        _S4197 = _S4198;
    }
    else
    {
        float _S4200 = _S4191 * _S4191;
        k_9 = _S4193 / _S4191;
        _S4195 = _S4200;
        _S4196 = 0.0f;
        _S4197 = 0.0f;
    }
    float2  _S4201 = make_float2 (k_9);
    float2  _S4202 = _S4190 * make_float2 (k_9);
    float u_26 = _S4202.x;
    float v_26 = _S4202.y;
    float r2_26 = u_26 * u_26 + v_26 * v_26;
    float _S4203 = k3_8 + r2_26 * k4_9;
    float _S4204 = k2_8 + r2_26 * _S4203;
    float _S4205 = k1_8 + r2_26 * _S4204;
    float radial_6 = 1.0f + r2_26 * _S4205;
    float2  _S4206 = make_float2 (radial_6);
    float _S4207 = _S4178 * u_26;
    float _S4208 = 2.0f * u_26;
    float _S4209 = r2_26 + _S4208 * u_26;
    float _S4210 = _S4182 * u_26;
    float _S4211 = 2.0f * v_26;
    float _S4212 = r2_26 + _S4211 * v_26;
    float2  _S4213 = _S4202 * make_float2 (radial_6) + make_float2 (_S4207 * v_26 + p2_9 * _S4209 + sx1_9 * r2_26, _S4210 * v_26 + p1_9 * _S4212 + sy1_9 * r2_26);
    float _S4214 = fx_30 * _S4213.x + cx_30;
    float _S4215 = fy_30 * _S4213.y + cy_30;
    float2  _S4216 = make_float2 (_S4214, _S4215);
    float2  _S4217 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4218 = length_0(_S4217);
    float _S4219 = vert2_c_11.z;
    float _S4220 = s_primal_ctx_atan2_0(_S4218, _S4219);
    bool _S4221 = _S4220 < 0.00100000004749745f;
    float _S4222;
    float _S4223;
    float _S4224;
    if(_S4221)
    {
        float _S4225 = 1.0f - _S4220 * _S4220 / 3.0f;
        float _S4226 = _S4219 * _S4219;
        k_9 = _S4225 / _S4219;
        _S4222 = 0.0f;
        _S4223 = _S4226;
        _S4224 = _S4225;
    }
    else
    {
        float _S4227 = _S4218 * _S4218;
        k_9 = _S4220 / _S4218;
        _S4222 = _S4227;
        _S4223 = 0.0f;
        _S4224 = 0.0f;
    }
    float2  _S4228 = make_float2 (k_9);
    float2  _S4229 = _S4217 * make_float2 (k_9);
    float u_27 = _S4229.x;
    float v_27 = _S4229.y;
    float r2_27 = u_27 * u_27 + v_27 * v_27;
    float _S4230 = k3_8 + r2_27 * k4_9;
    float _S4231 = k2_8 + r2_27 * _S4230;
    float _S4232 = k1_8 + r2_27 * _S4231;
    float radial_7 = 1.0f + r2_27 * _S4232;
    float2  _S4233 = make_float2 (radial_7);
    float _S4234 = _S4178 * u_27;
    float _S4235 = 2.0f * u_27;
    float _S4236 = r2_27 + _S4235 * u_27;
    float _S4237 = _S4182 * u_27;
    float _S4238 = 2.0f * v_27;
    float _S4239 = r2_27 + _S4238 * v_27;
    float2  _S4240 = _S4229 * make_float2 (radial_7) + make_float2 (_S4234 * v_27 + p2_9 * _S4236 + sx1_9 * r2_27, _S4237 * v_27 + p1_9 * _S4239 + sy1_9 * r2_27);
    float _S4241 = fx_30 * _S4240.x + cx_30;
    float _S4242 = fy_30 * _S4240.y + cy_30;
    float2  _S4243 = make_float2 (_S4241, _S4242);
    float2  e0_15 = _S4216 - _S4189;
    float2  e1_15 = _S4243 - _S4216;
    float2  e2_7 = _S4189 - _S4243;
    float _S4244 = e0_15.x;
    float _S4245 = e1_15.y;
    float _S4246 = e0_15.y;
    float _S4247 = e1_15.x;
    float _S4248 = _S4244 * _S4245 - _S4246 * _S4247;
    float _S4249 = 1.0f - hardness_15.y;
    float _S4250 = -1.0f / _S4249;
    float _S4251 = _S4249 * _S4249;
    float _S4252 = s_primal_ctx_max_0(_S4187, _S4214);
    float _S4253 = s_primal_ctx_min_0(_S4187, _S4214);
    float _S4254 = s_primal_ctx_max_0(_S4188, _S4215);
    float _S4255 = s_primal_ctx_min_0(_S4188, _S4215);
    float3  _S4256 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S4257 = length_1(_S4256) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4258 = transpose_0(R_30);
    float3  _S4259 = mean_30 - - s_primal_ctx_mul_1(_S4258, t_29);
    float _S4260 = _S4259.x;
    float _S4261 = _S4259.y;
    float _S4262 = _S4259.z;
    float _S4263 = _S4260 * _S4260 + _S4261 * _S4261 + _S4262 * _S4262;
    float _S4264 = s_primal_ctx_sqrt_0(_S4263);
    float x_61 = _S4260 / _S4264;
    float3  _S4265 = make_float3 (x_61);
    float _S4266 = _S4264 * _S4264;
    float y_29 = _S4261 / _S4264;
    float z_26 = _S4262 / _S4264;
    float3  _S4267 = make_float3 (z_26);
    float _S4268 = - y_29;
    float3  _S4269 = make_float3 (_S4268);
    float z2_58 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_61 * x_61 - y_29 * y_29;
    float _S4270 = 2.0f * x_61;
    float fS1_26 = _S4270 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_58 - 0.31539157032966614f;
    float3  _S4271 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_61;
    float3  _S4272 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4273 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4274 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4275 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_58 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4276 = 1.86588168144226074f * z2_58 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4276;
    float3  _S4277 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_61;
    float3  _S4278 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4279 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4280 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4281 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_61 * fC1_26 - y_29 * fS1_26);
    float3  _S4282 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_61 * fS1_26 + y_29 * fC1_26);
    float3  _S4283 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4268) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_61) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4284 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4285 = make_float3 (0.0f);
    float3  _S4286 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4287 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4288 = make_float3 (_S4287);
    float3  _S4289 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4287);
    float3  _S4290 = _S4286 + _S4289 + make_float3 (0.5f);
    float3  _S4291 = _S4286 - _S4289 + make_float3 (0.5f);
    float3  _S4292 = vert1_c_11 - vert0_c_11;
    float3  _S4293 = vert2_c_11 - vert0_c_11;
    float3  _S4294 = s_primal_ctx_cross_0(_S4292, _S4293);
    float3  _S4295 = normalize_0(_S4294);
    float3  _S4296 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4295, mean_c_26)))))) * v_normal_3;
    float3  _S4297 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4298;
    (&_S4298)->primal_0 = _S4295;
    (&_S4298)->differential_0 = _S4297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4299;
    (&_S4299)->primal_0 = mean_c_26;
    (&_S4299)->differential_0 = _S4297;
    s_bwd_prop_dot_0(&_S4298, &_S4299, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4300 = _S4299;
    float3  _S4301 = _S4296 + _S4298.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4302;
    (&_S4302)->primal_0 = _S4294;
    (&_S4302)->differential_0 = _S4297;
    s_bwd_normalize_impl_0(&_S4302, _S4301);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4303;
    (&_S4303)->primal_0 = _S4292;
    (&_S4303)->differential_0 = _S4297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4304;
    (&_S4304)->primal_0 = _S4293;
    (&_S4304)->differential_0 = _S4297;
    s_bwd_prop_cross_0(&_S4303, &_S4304, _S4302.differential_0);
    float3  _S4305 = - _S4304.differential_0;
    float3  _S4306 = - _S4303.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4307;
    (&_S4307)->primal_0 = _S4291;
    (&_S4307)->differential_0 = _S4297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4308;
    (&_S4308)->primal_0 = _S4285;
    (&_S4308)->differential_0 = _S4297;
    s_bwd_prop_max_0(&_S4307, &_S4308, (*v_rgbs_1)[int(2)]);
    float3  _S4309 = - _S4307.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4310;
    (&_S4310)->primal_0 = _S4290;
    (&_S4310)->differential_0 = _S4297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4311;
    (&_S4311)->primal_0 = _S4285;
    (&_S4311)->differential_0 = _S4297;
    s_bwd_prop_max_0(&_S4310, &_S4311, (*v_rgbs_1)[int(1)]);
    float3  _S4312 = _S4288 * (_S4309 + _S4310.differential_0);
    float3  _S4313 = _S4307.differential_0 + _S4310.differential_0;
    float3  _S4314 = make_float3 (0.5f) * - _S4313;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4315;
    (&_S4315)->primal_0 = _S4284;
    (&_S4315)->differential_0 = _S4297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4316;
    (&_S4316)->primal_0 = _S4285;
    (&_S4316)->differential_0 = _S4297;
    s_bwd_prop_max_0(&_S4315, &_S4316, (*v_rgbs_1)[int(0)]);
    float3  _S4317 = _S4314 + _S4315.differential_0;
    float3  _S4318 = _S4313 + _S4315.differential_0;
    float3  _S4319 = _S4282 * _S4318;
    float3  _S4320 = (*sh_coeffs_26)[int(15)] * _S4318;
    float3  _S4321 = _S4280 * _S4318;
    float3  _S4322 = (*sh_coeffs_26)[int(14)] * _S4318;
    float3  _S4323 = _S4278 * _S4318;
    float3  _S4324 = (*sh_coeffs_26)[int(13)] * _S4318;
    float3  _S4325 = _S4277 * _S4318;
    float3  _S4326 = (*sh_coeffs_26)[int(12)] * _S4318;
    float3  _S4327 = _S4279 * _S4318;
    float3  _S4328 = (*sh_coeffs_26)[int(11)] * _S4318;
    float3  _S4329 = _S4281 * _S4318;
    float3  _S4330 = (*sh_coeffs_26)[int(10)] * _S4318;
    float3  _S4331 = _S4283 * _S4318;
    float3  _S4332 = (*sh_coeffs_26)[int(9)] * _S4318;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4332.x + _S4332.y + _S4332.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4320.x + _S4320.y + _S4320.z);
    float _S4333 = _S4330.x + _S4330.y + _S4330.z;
    float _S4334 = _S4322.x + _S4322.y + _S4322.z;
    float _S4335 = _S4328.x + _S4328.y + _S4328.z;
    float _S4336 = _S4324.x + _S4324.y + _S4324.z;
    float _S4337 = _S4326.x + _S4326.y + _S4326.z;
    float _S4338 = - s_diff_fC2_T_8;
    float3  _S4339 = _S4274 * _S4318;
    float3  _S4340 = (*sh_coeffs_26)[int(8)] * _S4318;
    float3  _S4341 = _S4272 * _S4318;
    float3  _S4342 = (*sh_coeffs_26)[int(7)] * _S4318;
    float3  _S4343 = _S4271 * _S4318;
    float3  _S4344 = (*sh_coeffs_26)[int(6)] * _S4318;
    float3  _S4345 = _S4273 * _S4318;
    float3  _S4346 = (*sh_coeffs_26)[int(5)] * _S4318;
    float3  _S4347 = _S4275 * _S4318;
    float3  _S4348 = (*sh_coeffs_26)[int(4)] * _S4318;
    float _S4349 = _S4346.x + _S4346.y + _S4346.z;
    float _S4350 = _S4342.x + _S4342.y + _S4342.z;
    float _S4351 = fTmp1B_26 * _S4333 + x_61 * s_diff_fS2_T_8 + y_29 * _S4338 + 0.54627424478530884f * (_S4348.x + _S4348.y + _S4348.z);
    float _S4352 = fTmp1B_26 * _S4334 + y_29 * s_diff_fS2_T_8 + x_61 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S4340.x + _S4340.y + _S4340.z);
    float _S4353 = y_29 * - _S4352;
    float _S4354 = x_61 * _S4352;
    float _S4355 = z_26 * (1.86588168144226074f * (z_26 * _S4337) + -2.28522896766662598f * (y_29 * _S4335 + x_61 * _S4336) + 0.94617468118667603f * (_S4344.x + _S4344.y + _S4344.z));
    float3  _S4356 = make_float3 (0.48860251903533936f) * _S4318;
    float3  _S4357 = - _S4356;
    float3  _S4358 = _S4265 * _S4357;
    float3  _S4359 = (*sh_coeffs_26)[int(3)] * _S4357;
    float3  _S4360 = _S4267 * _S4356;
    float3  _S4361 = (*sh_coeffs_26)[int(2)] * _S4356;
    float3  _S4362 = _S4269 * _S4356;
    float3  _S4363 = (*sh_coeffs_26)[int(1)] * _S4356;
    float _S4364 = (_S4276 * _S4337 + 1.44530570507049561f * (fS1_26 * _S4333 + fC1_26 * _S4334) + -1.09254848957061768f * (y_29 * _S4349 + x_61 * _S4350) + _S4355 + _S4355 + _S4361.x + _S4361.y + _S4361.z) / _S4266;
    float _S4365 = _S4264 * _S4364;
    float _S4366 = (fTmp0C_26 * _S4335 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S4338 + fTmp0B_26 * _S4349 + _S4270 * _S4351 + _S4353 + _S4353 + - (_S4363.x + _S4363.y + _S4363.z)) / _S4266;
    float _S4367 = _S4264 * _S4366;
    float _S4368 = (fTmp0C_26 * _S4336 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S4350 + 2.0f * (y_29 * _S4351) + _S4354 + _S4354 + _S4359.x + _S4359.y + _S4359.z) / _S4266;
    float _S4369 = _S4264 * _S4368;
    float _S4370 = _S4262 * - _S4364 + _S4261 * - _S4366 + _S4260 * - _S4368;
    DiffPair_float_0 _S4371;
    (&_S4371)->primal_0 = _S4263;
    (&_S4371)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4371, _S4370);
    float _S4372 = _S4262 * _S4371.differential_0;
    float _S4373 = _S4261 * _S4371.differential_0;
    float _S4374 = _S4260 * _S4371.differential_0;
    float3  _S4375 = make_float3 (0.282094806432724f) * _S4318;
    float3  _S4376 = make_float3 (_S4369 + _S4374 + _S4374, _S4367 + _S4373 + _S4373, _S4365 + _S4372 + _S4372);
    float3  _S4377 = - - _S4376;
    Matrix<float, 3, 3>  _S4378 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4379;
    (&_S4379)->primal_0 = _S4258;
    (&_S4379)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4380;
    (&_S4380)->primal_0 = t_29;
    (&_S4380)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4379, &_S4380, _S4377);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4381 = _S4380;
    Matrix<float, 3, 3>  _S4382 = transpose_0(_S4379.differential_0);
    DiffPair_float_0 _S4383;
    (&_S4383)->primal_0 = _S4257;
    (&_S4383)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4383, v_depth_10);
    float _S4384 = 0.3333333432674408f * _S4383.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4385;
    (&_S4385)->primal_0 = _S4256;
    (&_S4385)->differential_0 = _S4297;
    s_bwd_length_impl_0(&_S4385, _S4384);
    DiffPair_float_0 _S4386;
    (&_S4386)->primal_0 = _S4255;
    (&_S4386)->differential_0 = 0.0f;
    DiffPair_float_0 _S4387;
    (&_S4387)->primal_0 = _S4242;
    (&_S4387)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4386, &_S4387, 0.0f);
    DiffPair_float_0 _S4388;
    (&_S4388)->primal_0 = _S4188;
    (&_S4388)->differential_0 = 0.0f;
    DiffPair_float_0 _S4389;
    (&_S4389)->primal_0 = _S4215;
    (&_S4389)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4388, &_S4389, _S4386.differential_0);
    DiffPair_float_0 _S4390;
    (&_S4390)->primal_0 = _S4254;
    (&_S4390)->differential_0 = 0.0f;
    DiffPair_float_0 _S4391;
    (&_S4391)->primal_0 = _S4242;
    (&_S4391)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4390, &_S4391, 0.0f);
    DiffPair_float_0 _S4392;
    (&_S4392)->primal_0 = _S4188;
    (&_S4392)->differential_0 = 0.0f;
    DiffPair_float_0 _S4393;
    (&_S4393)->primal_0 = _S4215;
    (&_S4393)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4392, &_S4393, _S4390.differential_0);
    DiffPair_float_0 _S4394;
    (&_S4394)->primal_0 = _S4253;
    (&_S4394)->differential_0 = 0.0f;
    DiffPair_float_0 _S4395;
    (&_S4395)->primal_0 = _S4241;
    (&_S4395)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4394, &_S4395, 0.0f);
    DiffPair_float_0 _S4396;
    (&_S4396)->primal_0 = _S4187;
    (&_S4396)->differential_0 = 0.0f;
    DiffPair_float_0 _S4397;
    (&_S4397)->primal_0 = _S4214;
    (&_S4397)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4396, &_S4397, _S4394.differential_0);
    DiffPair_float_0 _S4398;
    (&_S4398)->primal_0 = _S4252;
    (&_S4398)->differential_0 = 0.0f;
    DiffPair_float_0 _S4399;
    (&_S4399)->primal_0 = _S4241;
    (&_S4399)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4398, &_S4399, 0.0f);
    DiffPair_float_0 _S4400;
    (&_S4400)->primal_0 = _S4187;
    (&_S4400)->differential_0 = 0.0f;
    DiffPair_float_0 _S4401;
    (&_S4401)->primal_0 = _S4214;
    (&_S4401)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4400, &_S4401, _S4398.differential_0);
    DiffPair_float_0 _S4402;
    (&_S4402)->primal_0 = _S4250;
    (&_S4402)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4402, 0.0f);
    float _S4403 = - (-1.0f * - (_S4402.differential_0 / _S4251));
    float2  _S4404 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4405;
    (&_S4405)->primal_0 = e2_7;
    (&_S4405)->differential_0 = _S4404;
    s_bwd_length_impl_1(&_S4405, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4406;
    (&_S4406)->primal_0 = e1_15;
    (&_S4406)->differential_0 = _S4404;
    s_bwd_length_impl_1(&_S4406, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4407;
    (&_S4407)->primal_0 = e0_15;
    (&_S4407)->differential_0 = _S4404;
    s_bwd_length_impl_1(&_S4407, -0.0f);
    DiffPair_float_0 _S4408;
    (&_S4408)->primal_0 = _S4248;
    (&_S4408)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4408, 0.0f);
    float _S4409 = - _S4408.differential_0;
    float2  _S4410 = _S4406.differential_0 + make_float2 (_S4246 * _S4409, _S4244 * _S4408.differential_0);
    float2  _S4411 = - _S4410;
    float2  _S4412 = _S4407.differential_0 + make_float2 (_S4245 * _S4408.differential_0, _S4247 * _S4409);
    float2  _S4413 = - _S4412;
    float2  _S4414 = - _S4405.differential_0 + _S4410;
    float _S4415 = fy_30 * (_S4387.differential_0 + _S4391.differential_0 + _S4414.y);
    float _S4416 = fx_30 * (_S4395.differential_0 + _S4399.differential_0 + _S4414.x);
    float2  _S4417 = make_float2 (_S4416, _S4415);
    float2  _S4418 = _S4229 * _S4417;
    float2  _S4419 = _S4233 * _S4417;
    float _S4420 = r2_27 * _S4415;
    float _S4421 = p1_9 * _S4415;
    float _S4422 = _S4239 * _S4415;
    float _S4423 = v_27 * _S4415;
    float _S4424 = u_27 * _S4423;
    float _S4425 = r2_27 * _S4416;
    float _S4426 = p2_9 * _S4416;
    float _S4427 = _S4236 * _S4416;
    float _S4428 = v_27 * _S4416;
    float _S4429 = u_27 * _S4428;
    float _S4430 = _S4418.x + _S4418.y;
    float _S4431 = r2_27 * _S4430;
    float _S4432 = r2_27 * _S4431;
    float _S4433 = r2_27 * _S4432;
    float _S4434 = r2_27 * _S4433;
    float _S4435 = sy1_9 * _S4415 + _S4421 + sx1_9 * _S4416 + _S4426 + _S4232 * _S4430 + _S4231 * _S4431 + _S4230 * _S4432 + k4_9 * _S4433;
    float _S4436 = v_27 * _S4435;
    float _S4437 = u_27 * _S4435;
    float _S4438 = _S4238 * _S4421 + 2.0f * (v_27 * _S4421) + _S4237 * _S4415 + _S4234 * _S4416 + _S4436 + _S4436;
    float _S4439 = _S4182 * _S4423 + _S4235 * _S4426 + 2.0f * (u_27 * _S4426) + _S4178 * _S4428 + _S4437 + _S4437;
    float3  _S4440 = _S4303.differential_0 + _S4385.differential_0;
    float3  _S4441 = _S4305 + _S4306 + _S4385.differential_0;
    float3  _S4442 = _S4304.differential_0 + _S4385.differential_0;
    FixedArray<float3 , 2>  _S4443;
    _S4443[int(0)] = _S4297;
    _S4443[int(1)] = _S4297;
    _S4443[int(1)] = _S4312;
    _S4443[int(0)] = _S4317;
    float3  _S4444 = _S4443[int(0)];
    float3  _S4445 = _S4443[int(1)];
    FixedArray<float3 , 16>  _S4446;
    _S4446[int(0)] = _S4297;
    _S4446[int(1)] = _S4297;
    _S4446[int(2)] = _S4297;
    _S4446[int(3)] = _S4297;
    _S4446[int(4)] = _S4297;
    _S4446[int(5)] = _S4297;
    _S4446[int(6)] = _S4297;
    _S4446[int(7)] = _S4297;
    _S4446[int(8)] = _S4297;
    _S4446[int(9)] = _S4297;
    _S4446[int(10)] = _S4297;
    _S4446[int(11)] = _S4297;
    _S4446[int(12)] = _S4297;
    _S4446[int(13)] = _S4297;
    _S4446[int(14)] = _S4297;
    _S4446[int(15)] = _S4297;
    _S4446[int(7)] = _S4341;
    _S4446[int(0)] = _S4375;
    _S4446[int(1)] = _S4362;
    _S4446[int(2)] = _S4360;
    _S4446[int(3)] = _S4358;
    _S4446[int(4)] = _S4347;
    _S4446[int(5)] = _S4345;
    _S4446[int(6)] = _S4343;
    _S4446[int(15)] = _S4319;
    _S4446[int(8)] = _S4339;
    _S4446[int(9)] = _S4331;
    _S4446[int(10)] = _S4329;
    _S4446[int(11)] = _S4327;
    _S4446[int(12)] = _S4325;
    _S4446[int(13)] = _S4323;
    _S4446[int(14)] = _S4321;
    float3  _S4447 = _S4446[int(0)];
    float3  _S4448 = _S4446[int(1)];
    float3  _S4449 = _S4446[int(2)];
    float3  _S4450 = _S4446[int(3)];
    float3  _S4451 = _S4446[int(4)];
    float3  _S4452 = _S4446[int(5)];
    float3  _S4453 = _S4446[int(6)];
    float3  _S4454 = _S4446[int(7)];
    float3  _S4455 = _S4446[int(8)];
    float3  _S4456 = _S4446[int(9)];
    float3  _S4457 = _S4446[int(10)];
    float3  _S4458 = _S4446[int(11)];
    float3  _S4459 = _S4446[int(12)];
    float3  _S4460 = _S4446[int(13)];
    float3  _S4461 = _S4446[int(14)];
    float3  _S4462 = _S4446[int(15)];
    float _S4463 = _S4396.differential_0 + _S4400.differential_0;
    float2  _S4464 = _S4405.differential_0 + _S4413;
    float _S4465 = _S4389.differential_0 + _S4393.differential_0;
    float _S4466 = _S4388.differential_0 + _S4392.differential_0;
    float2  _S4467 = _S4411 + _S4412;
    float _S4468 = _S4397.differential_0 + _S4401.differential_0;
    float2  _S4469 = make_float2 (0.0f, _S4403);
    float2  _S4470 = _S4419 + make_float2 (_S4439, _S4438);
    float2  _S4471 = _S4217 * _S4470;
    float2  _S4472 = _S4228 * _S4470;
    float _S4473 = _S4471.x + _S4471.y;
    if(_S4221)
    {
        float _S4474 = _S4473 / _S4223;
        float _S4475 = _S4224 * - _S4474;
        float _S4476 = _S4220 * (0.3333333432674408f * - (_S4219 * _S4474));
        k_9 = _S4476 + _S4476;
        _S4222 = _S4475;
        _S4223 = 0.0f;
    }
    else
    {
        float _S4477 = _S4473 / _S4222;
        float _S4478 = _S4220 * - _S4477;
        k_9 = _S4218 * _S4477;
        _S4222 = 0.0f;
        _S4223 = _S4478;
    }
    DiffPair_float_0 _S4479;
    (&_S4479)->primal_0 = _S4218;
    (&_S4479)->differential_0 = 0.0f;
    DiffPair_float_0 _S4480;
    (&_S4480)->primal_0 = _S4219;
    (&_S4480)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4479, &_S4480, k_9);
    float _S4481 = _S4480.differential_0 + _S4222;
    float _S4482 = _S4479.differential_0 + _S4223;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4483;
    (&_S4483)->primal_0 = _S4217;
    (&_S4483)->differential_0 = _S4404;
    s_bwd_length_impl_1(&_S4483, _S4482);
    float2  _S4484 = _S4483.differential_0 + _S4472;
    float _S4485 = fy_30 * (_S4467.y + _S4465);
    float _S4486 = fx_30 * (_S4467.x + _S4468);
    float2  _S4487 = make_float2 (_S4486, _S4485);
    float2  _S4488 = _S4202 * _S4487;
    float _S4489 = p1_9 * _S4485;
    float _S4490 = v_26 * _S4485;
    float _S4491 = p2_9 * _S4486;
    float _S4492 = v_26 * _S4486;
    float _S4493 = _S4488.x + _S4488.y;
    float _S4494 = r2_26 * _S4493;
    float _S4495 = r2_26 * _S4494;
    float _S4496 = r2_26 * _S4495;
    float _S4497 = sy1_9 * _S4485 + _S4489 + sx1_9 * _S4486 + _S4491 + _S4205 * _S4493 + _S4204 * _S4494 + _S4203 * _S4495 + k4_9 * _S4496;
    float _S4498 = v_26 * _S4497;
    float _S4499 = u_26 * _S4497;
    float3  _S4500 = _S4442 + make_float3 (_S4484.x, _S4484.y, _S4481);
    float2  _S4501 = _S4206 * _S4487 + make_float2 (_S4182 * _S4490 + _S4208 * _S4491 + 2.0f * (u_26 * _S4491) + _S4178 * _S4492 + _S4499 + _S4499, _S4211 * _S4489 + 2.0f * (v_26 * _S4489) + _S4210 * _S4485 + _S4207 * _S4486 + _S4498 + _S4498);
    float _S4502 = u_26 * _S4490 + _S4424;
    float _S4503 = u_26 * _S4492 + _S4429;
    float _S4504 = r2_26 * _S4485 + _S4420;
    float _S4505 = r2_26 * _S4496 + _S4434;
    float _S4506 = _S4496 + _S4433;
    float _S4507 = _S4212 * _S4485 + _S4422;
    float _S4508 = _S4495 + _S4432;
    float _S4509 = r2_26 * _S4486 + _S4425;
    float _S4510 = _S4209 * _S4486 + _S4427;
    float _S4511 = _S4494 + _S4431;
    float2  _S4512 = _S4190 * _S4501;
    float2  _S4513 = _S4201 * _S4501;
    float _S4514 = _S4512.x + _S4512.y;
    if(_S4194)
    {
        float _S4515 = _S4514 / _S4196;
        float _S4516 = _S4197 * - _S4515;
        float _S4517 = _S4193 * (0.3333333432674408f * - (_S4192 * _S4515));
        k_9 = _S4517 + _S4517;
        _S4195 = _S4516;
        _S4196 = 0.0f;
    }
    else
    {
        float _S4518 = _S4514 / _S4195;
        float _S4519 = _S4193 * - _S4518;
        k_9 = _S4191 * _S4518;
        _S4195 = 0.0f;
        _S4196 = _S4519;
    }
    DiffPair_float_0 _S4520;
    (&_S4520)->primal_0 = _S4191;
    (&_S4520)->differential_0 = 0.0f;
    DiffPair_float_0 _S4521;
    (&_S4521)->primal_0 = _S4192;
    (&_S4521)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4520, &_S4521, k_9);
    float _S4522 = _S4521.differential_0 + _S4195;
    float _S4523 = _S4520.differential_0 + _S4196;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4524;
    (&_S4524)->primal_0 = _S4190;
    (&_S4524)->differential_0 = _S4404;
    s_bwd_length_impl_1(&_S4524, _S4523);
    float2  _S4525 = _S4524.differential_0 + _S4513;
    float _S4526 = fy_30 * (_S4464.y + _S4466);
    float _S4527 = fx_30 * (_S4464.x + _S4463);
    float2  _S4528 = make_float2 (_S4527, _S4526);
    float2  _S4529 = _S4173 * _S4528;
    float _S4530 = p1_9 * _S4526;
    float _S4531 = v_25 * _S4526;
    float _S4532 = p2_9 * _S4527;
    float _S4533 = v_25 * _S4527;
    float _S4534 = _S4529.x + _S4529.y;
    float _S4535 = r2_25 * _S4534;
    float _S4536 = r2_25 * _S4535;
    float _S4537 = r2_25 * _S4536;
    float _S4538 = sy1_9 * _S4526 + _S4530 + sx1_9 * _S4527 + _S4532 + _S4176 * _S4534 + _S4175 * _S4535 + _S4174 * _S4536 + k4_9 * _S4537;
    float _S4539 = v_25 * _S4538;
    float _S4540 = u_25 * _S4538;
    float2  _S4541 = make_float2 (r2_25 * _S4527 + _S4509, r2_25 * _S4526 + _S4504);
    float2  _S4542 = make_float2 (_S4185 * _S4526 + 2.0f * (u_25 * _S4533 + _S4503) + _S4507, 2.0f * (u_25 * _S4531 + _S4502) + _S4181 * _S4527 + _S4510);
    float4  _S4543 = make_float4 (_S4535 + _S4511, _S4536 + _S4508, _S4537 + _S4506, r2_25 * _S4537 + _S4505);
    float3  _S4544 = _S4440 + make_float3 (_S4525.x, _S4525.y, _S4522);
    float2  _S4545 = _S4177 * _S4528 + make_float2 (_S4182 * _S4531 + _S4180 * _S4532 + 2.0f * (u_25 * _S4532) + _S4178 * _S4533 + _S4540 + _S4540, _S4184 * _S4530 + 2.0f * (v_25 * _S4530) + _S4183 * _S4526 + _S4179 * _S4527 + _S4539 + _S4539);
    CameraDistortion_0 _S4546 = CameraDistortion_x24_syn_dzero_0();
    (&_S4546)->thin_prism_coeffs_0 = _S4541;
    (&_S4546)->tangential_coeffs_0 = _S4542;
    (&_S4546)->radial_coeffs_0 = _S4543;
    float2  _S4547 = _S4161 * _S4545;
    float2  _S4548 = _S4172 * _S4545;
    float _S4549 = _S4547.x + _S4547.y;
    if(_S4165)
    {
        float _S4550 = _S4549 / _S4167;
        float _S4551 = _S4168 * - _S4550;
        float _S4552 = _S4164 * (0.3333333432674408f * - (_S4163 * _S4550));
        k_9 = _S4552 + _S4552;
        _S4166 = _S4551;
        _S4167 = 0.0f;
    }
    else
    {
        float _S4553 = _S4549 / _S4166;
        float _S4554 = _S4164 * - _S4553;
        k_9 = _S4162 * _S4553;
        _S4166 = 0.0f;
        _S4167 = _S4554;
    }
    DiffPair_float_0 _S4555;
    (&_S4555)->primal_0 = _S4162;
    (&_S4555)->differential_0 = 0.0f;
    DiffPair_float_0 _S4556;
    (&_S4556)->primal_0 = _S4163;
    (&_S4556)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4555, &_S4556, k_9);
    float _S4557 = _S4556.differential_0 + _S4166;
    float _S4558 = _S4555.differential_0 + _S4167;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4559;
    (&_S4559)->primal_0 = _S4161;
    (&_S4559)->differential_0 = _S4404;
    s_bwd_length_impl_1(&_S4559, _S4558);
    float2  _S4560 = _S4559.differential_0 + _S4548;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4561;
    (&_S4561)->primal_0 = vert2_c_11;
    (&_S4561)->differential_0 = _S4297;
    s_bwd_length_impl_0(&_S4561, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4562;
    (&_S4562)->primal_0 = vert1_c_11;
    (&_S4562)->differential_0 = _S4297;
    s_bwd_length_impl_0(&_S4562, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4563;
    (&_S4563)->primal_0 = vert0_c_11;
    (&_S4563)->differential_0 = _S4297;
    s_bwd_length_impl_0(&_S4563, 0.0f);
    float3  _S4564 = _S4561.differential_0 + _S4500;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4565;
    (&_S4565)->primal_0 = R_30;
    (&_S4565)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4566;
    (&_S4566)->primal_0 = vert2_7;
    (&_S4566)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4565, &_S4566, _S4564);
    float3  _S4567 = _S4562.differential_0 + _S4544;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4568;
    (&_S4568)->primal_0 = R_30;
    (&_S4568)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4569;
    (&_S4569)->primal_0 = vert1_7;
    (&_S4569)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4568, &_S4569, _S4567);
    float3  _S4570 = _S4563.differential_0 + _S4441 + make_float3 (_S4560.x, _S4560.y, _S4557);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4571;
    (&_S4571)->primal_0 = R_30;
    (&_S4571)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4572;
    (&_S4572)->primal_0 = vert0_7;
    (&_S4572)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4571, &_S4572, _S4570);
    float3  _S4573 = _S4566.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4574;
    (&_S4574)->primal_0 = _S4154;
    (&_S4574)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4575;
    (&_S4575)->primal_0 = _S4159;
    (&_S4575)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4574, &_S4575, _S4573);
    float _S4576 = - _S4575.differential_0.y;
    float _S4577 = _S4158 * _S4575.differential_0.x;
    float _S4578 = - (_S4150 * _S4575.differential_0.x);
    float3  _S4579 = _S4569.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4580;
    (&_S4580)->primal_0 = _S4154;
    (&_S4580)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4581;
    (&_S4581)->primal_0 = _S4157;
    (&_S4581)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4580, &_S4581, _S4579);
    float _S4582 = _S4150 * _S4581.differential_0.x;
    float _S4583 = _S4156 * _S4581.differential_0.x;
    float3  _S4584 = _S4572.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4585;
    (&_S4585)->primal_0 = _S4154;
    (&_S4585)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4586;
    (&_S4586)->primal_0 = _S4155;
    (&_S4586)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4585, &_S4586, _S4584);
    Matrix<float, 3, 3>  _S4587 = transpose_0(_S4574.differential_0 + _S4580.differential_0 + _S4585.differential_0);
    float _S4588 = 2.0f * - _S4587.rows[int(2)].z;
    float _S4589 = 2.0f * _S4587.rows[int(2)].y;
    float _S4590 = 2.0f * _S4587.rows[int(2)].x;
    float _S4591 = 2.0f * _S4587.rows[int(1)].z;
    float _S4592 = 2.0f * - _S4587.rows[int(1)].y;
    float _S4593 = 2.0f * _S4587.rows[int(1)].x;
    float _S4594 = 2.0f * _S4587.rows[int(0)].z;
    float _S4595 = 2.0f * _S4587.rows[int(0)].y;
    float _S4596 = 2.0f * - _S4587.rows[int(0)].x;
    float _S4597 = - _S4593 + _S4595;
    float _S4598 = _S4590 + - _S4594;
    float _S4599 = - _S4589 + _S4591;
    float _S4600 = _S4589 + _S4591;
    float _S4601 = _S4590 + _S4594;
    float _S4602 = _S4593 + _S4595;
    float _S4603 = quat_31.w * (_S4592 + _S4596);
    float _S4604 = quat_31.z * (_S4588 + _S4596);
    float _S4605 = quat_31.y * (_S4588 + _S4592);
    float _S4606 = quat_31.x * _S4597 + quat_31.z * _S4600 + quat_31.y * _S4601 + _S4603 + _S4603;
    float _S4607 = quat_31.x * _S4598 + quat_31.w * _S4600 + quat_31.y * _S4602 + _S4604 + _S4604;
    float _S4608 = quat_31.x * _S4599 + quat_31.w * _S4601 + quat_31.z * _S4602 + _S4605 + _S4605;
    float _S4609 = quat_31.w * _S4597 + quat_31.z * _S4598 + quat_31.y * _S4599;
    float _S4610 = _S4578 + _S4582;
    float _S4611 = 0.5f * - _S4610;
    float _S4612 = _S4576 + _S4581.differential_0.y;
    DiffPair_float_0 _S4613;
    (&_S4613)->primal_0 = _S4151;
    (&_S4613)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4613, _S4612);
    float _S4614 = _S4611 + _S4613.differential_0;
    float _S4615 = _S4577 + _S4583 + _S4586.differential_0.x;
    DiffPair_float_0 _S4616;
    (&_S4616)->primal_0 = _S4149;
    (&_S4616)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4616, _S4615);
    float _S4617 = _S4611 + _S4616.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4618;
    (&_S4618)->primal_0 = mean_c_26;
    (&_S4618)->differential_0 = _S4297;
    s_bwd_length_impl_0(&_S4618, 0.0f);
    float3  _S4619 = _S4618.differential_0 + _S4300.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4620;
    (&_S4620)->primal_0 = R_30;
    (&_S4620)->differential_0 = _S4378;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4621;
    (&_S4621)->primal_0 = mean_30;
    (&_S4621)->differential_0 = _S4297;
    s_bwd_prop_mul_1(&_S4620, &_S4621, _S4619);
    float3  _S4622 = _S4564 + _S4567 + _S4570 + _S4619 + _S4381.differential_0;
    Matrix<float, 3, 3>  _S4623 = _S4565.differential_0 + _S4568.differential_0 + _S4571.differential_0 + _S4620.differential_0 + _S4382;
    float3  _S4624 = make_float3 (_S4617, _S4614, _S4610);
    float4  _S4625 = make_float4 (0.0f);
    *&((&_S4625)->w) = _S4606;
    *&((&_S4625)->z) = _S4607;
    *&((&_S4625)->y) = _S4608;
    *&((&_S4625)->x) = _S4609;
    float4  _S4626 = _S4625;
    float3  _S4627 = _S4573 + _S4579 + _S4584 + _S4621.differential_0 + _S4376;
    *v_mean_10 = _S4627;
    *v_quat_9 = _S4626;
    *v_scale_9 = _S4624;
    *v_hardness_5 = _S4469;
    (*v_sh_coeffs_8)[int(0)] = _S4447;
    (*v_sh_coeffs_8)[int(1)] = _S4448;
    (*v_sh_coeffs_8)[int(2)] = _S4449;
    (*v_sh_coeffs_8)[int(3)] = _S4450;
    (*v_sh_coeffs_8)[int(4)] = _S4451;
    (*v_sh_coeffs_8)[int(5)] = _S4452;
    (*v_sh_coeffs_8)[int(6)] = _S4453;
    (*v_sh_coeffs_8)[int(7)] = _S4454;
    (*v_sh_coeffs_8)[int(8)] = _S4455;
    (*v_sh_coeffs_8)[int(9)] = _S4456;
    (*v_sh_coeffs_8)[int(10)] = _S4457;
    (*v_sh_coeffs_8)[int(11)] = _S4458;
    (*v_sh_coeffs_8)[int(12)] = _S4459;
    (*v_sh_coeffs_8)[int(13)] = _S4460;
    (*v_sh_coeffs_8)[int(14)] = _S4461;
    (*v_sh_coeffs_8)[int(15)] = _S4462;
    (*v_ch_coeffs_3)[int(0)] = _S4444;
    (*v_ch_coeffs_3)[int(1)] = _S4445;
    *v_R_9 = _S4623;
    *v_t_9 = _S4622;
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
    float u_28 = d_0 * dot_0(- q_2, v2v0_0);
    float v_28 = d_0 * dot_0(q_2, v1v0_0);
    float t_30 = d_0 * dot_0(- n_0, rov0_0);
    bool _S4628;
    if(u_28 >= 0.0f)
    {
        _S4628 = v_28 >= 0.0f;
    }
    else
    {
        _S4628 = false;
    }
    if(_S4628)
    {
        _S4628 = (u_28 + v_28) <= 1.0f;
    }
    else
    {
        _S4628 = false;
    }
    if(_S4628)
    {
        _S4628 = t_30 >= 0.0f;
    }
    else
    {
        _S4628 = false;
    }
    if(!_S4628)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_28), (v_28)))), ((F32_sqrt((0.5f))) * (1.0f - u_28 - v_28)))) * (2.0f + (F32_sqrt((2.0f))));
    float h_0 = clamp_0(hardness_16.y, 0.0f, 0.95999997854232788f);
    float o_0 = hardness_16.x;
    float _S4629;
    if(opac_0 < 0.0f)
    {
        _S4629 = 0.0f;
    }
    else
    {
        _S4629 = (F32_min((o_0 * (F32_pow((opac_0), (1.0f - h_0)))), (0.99500000476837158f)));
    }
    return _S4629;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_4, float _s_dOut_12)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4630 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4631 = *dpray_d_4;
    float3  v1v0_1 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_1 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_1 = (*dpray_o_4).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S4632 = s_primal_ctx_cross_0(v1v0_1, v2v0_1);
    float3  _S4633 = s_primal_ctx_cross_0(rov0_1, (*dpray_d_4).primal_0);
    float _S4634 = s_primal_ctx_dot_0((*dpray_d_4).primal_0, _S4632);
    float d_1 = 1.0f / _S4634;
    float _S4635 = _S4634 * _S4634;
    float3  _S4636 = - _S4633;
    float _S4637 = s_primal_ctx_dot_0(_S4636, v2v0_1);
    float u_29 = d_1 * _S4637;
    float _S4638 = s_primal_ctx_dot_0(_S4633, v1v0_1);
    float v_29 = d_1 * _S4638;
    float3  _S4639 = - _S4632;
    float t_31 = d_1 * s_primal_ctx_dot_0(_S4639, rov0_1);
    bool _S4640;
    if(u_29 >= 0.0f)
    {
        _S4640 = v_29 >= 0.0f;
    }
    else
    {
        _S4640 = false;
    }
    if(_S4640)
    {
        _S4640 = (u_29 + v_29) <= 1.0f;
    }
    else
    {
        _S4640 = false;
    }
    if(_S4640)
    {
        _S4640 = t_31 >= 0.0f;
    }
    else
    {
        _S4640 = false;
    }
    bool _S4641 = !!_S4640;
    float _S4642;
    float _S4643;
    float _S4644;
    float _S4645;
    float _S4646;
    float _S4647;
    float _S4648;
    float _S4649;
    float _S4650;
    float _S4651;
    if(_S4641)
    {
        float _S4652 = s_primal_ctx_min_0(u_29, v_29);
        float _S4653 = s_primal_ctx_sqrt_0(0.5f);
        float _S4654 = _S4653 * (1.0f - u_29 - v_29);
        float _S4655 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S4652, _S4654) * _S4655;
        float _S4656 = _S4630.primal_0.y;
        float _S4657 = s_primal_ctx_clamp_0(_S4656, 0.0f, 0.95999997854232788f);
        float o_1 = _S4630.primal_0.x;
        bool _S4658 = opac_1 < 0.0f;
        if(_S4658)
        {
            _S4642 = 0.0f;
            _S4643 = 0.0f;
            _S4644 = 0.0f;
        }
        else
        {
            float _S4659 = 1.0f - _S4657;
            float _S4660 = s_primal_ctx_pow_0(opac_1, _S4659);
            _S4642 = o_1 * _S4660;
            _S4643 = _S4660;
            _S4644 = _S4659;
        }
        float _S4661 = _S4643;
        float _S4662 = _S4644;
        _S4640 = _S4658;
        _S4643 = o_1;
        _S4644 = _S4661;
        _S4645 = opac_1;
        _S4646 = _S4662;
        _S4647 = _S4656;
        _S4648 = _S4655;
        _S4649 = _S4652;
        _S4650 = _S4654;
        _S4651 = _S4653;
    }
    else
    {
        _S4640 = false;
        _S4642 = 0.0f;
        _S4643 = 0.0f;
        _S4644 = 0.0f;
        _S4645 = 0.0f;
        _S4646 = 0.0f;
        _S4647 = 0.0f;
        _S4648 = 0.0f;
        _S4649 = 0.0f;
        _S4650 = 0.0f;
        _S4651 = 0.0f;
    }
    float2  _S4663 = make_float2 (0.0f);
    float2  _S4664;
    if(_S4641)
    {
        if(_S4640)
        {
            _S4642 = 0.0f;
            _S4643 = 0.0f;
            _S4644 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S4665;
            (&_S4665)->primal_0 = _S4642;
            (&_S4665)->differential_0 = 0.0f;
            DiffPair_float_0 _S4666;
            (&_S4666)->primal_0 = 0.99500000476837158f;
            (&_S4666)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S4665, &_S4666, _s_dOut_12);
            float _S4667 = _S4643 * _S4665.differential_0;
            float _S4668 = _S4644 * _S4665.differential_0;
            DiffPair_float_0 _S4669;
            (&_S4669)->primal_0 = _S4645;
            (&_S4669)->differential_0 = 0.0f;
            DiffPair_float_0 _S4670;
            (&_S4670)->primal_0 = _S4646;
            (&_S4670)->differential_0 = 0.0f;
            s_bwd_prop_pow_0(&_S4669, &_S4670, _S4667);
            float _S4671 = - _S4670.differential_0;
            _S4642 = _S4668;
            _S4643 = _S4671;
            _S4644 = _S4669.differential_0;
        }
        DiffPair_float_0 _S4672;
        (&_S4672)->primal_0 = _S4647;
        (&_S4672)->differential_0 = 0.0f;
        DiffPair_float_0 _S4673;
        (&_S4673)->primal_0 = 0.0f;
        (&_S4673)->differential_0 = 0.0f;
        DiffPair_float_0 _S4674;
        (&_S4674)->primal_0 = 0.95999997854232788f;
        (&_S4674)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S4672, &_S4673, &_S4674, _S4643);
        float _S4675 = _S4648 * _S4644;
        DiffPair_float_0 _S4676;
        (&_S4676)->primal_0 = _S4649;
        (&_S4676)->differential_0 = 0.0f;
        DiffPair_float_0 _S4677;
        (&_S4677)->primal_0 = _S4650;
        (&_S4677)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4676, &_S4677, _S4675);
        float _S4678 = - (_S4651 * _S4677.differential_0);
        DiffPair_float_0 _S4679;
        (&_S4679)->primal_0 = u_29;
        (&_S4679)->differential_0 = 0.0f;
        DiffPair_float_0 _S4680;
        (&_S4680)->primal_0 = v_29;
        (&_S4680)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4679, &_S4680, _S4676.differential_0);
        float2  _S4681 = make_float2 (_S4642, _S4672.differential_0);
        float _S4682 = _S4678 + _S4680.differential_0;
        _S4642 = _S4678 + _S4679.differential_0;
        _S4643 = _S4682;
        _S4664 = _S4681;
    }
    else
    {
        _S4642 = 0.0f;
        _S4643 = 0.0f;
        _S4664 = _S4663;
    }
    float3  _S4683 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4684;
    (&_S4684)->primal_0 = _S4639;
    (&_S4684)->differential_0 = _S4683;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4685;
    (&_S4685)->primal_0 = rov0_1;
    (&_S4685)->differential_0 = _S4683;
    s_bwd_prop_dot_0(&_S4684, &_S4685, 0.0f);
    float3  _S4686 = - _S4684.differential_0;
    float _S4687 = d_1 * _S4643;
    float _S4688 = _S4638 * _S4643;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4689;
    (&_S4689)->primal_0 = _S4633;
    (&_S4689)->differential_0 = _S4683;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4690;
    (&_S4690)->primal_0 = v1v0_1;
    (&_S4690)->differential_0 = _S4683;
    s_bwd_prop_dot_0(&_S4689, &_S4690, _S4687);
    float _S4691 = d_1 * _S4642;
    float _S4692 = _S4637 * _S4642;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4693;
    (&_S4693)->primal_0 = _S4636;
    (&_S4693)->differential_0 = _S4683;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4694;
    (&_S4694)->primal_0 = v2v0_1;
    (&_S4694)->differential_0 = _S4683;
    s_bwd_prop_dot_0(&_S4693, &_S4694, _S4691);
    float3  _S4695 = - _S4693.differential_0;
    float _S4696 = - ((_S4688 + _S4692) / _S4635);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4697;
    (&_S4697)->primal_0 = _S4631.primal_0;
    (&_S4697)->differential_0 = _S4683;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4698;
    (&_S4698)->primal_0 = _S4632;
    (&_S4698)->differential_0 = _S4683;
    s_bwd_prop_dot_0(&_S4697, &_S4698, _S4696);
    float3  _S4699 = _S4689.differential_0 + _S4695;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4700;
    (&_S4700)->primal_0 = rov0_1;
    (&_S4700)->differential_0 = _S4683;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4701;
    (&_S4701)->primal_0 = _S4631.primal_0;
    (&_S4701)->differential_0 = _S4683;
    s_bwd_prop_cross_0(&_S4700, &_S4701, _S4699);
    float3  _S4702 = _S4686 + _S4698.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4703;
    (&_S4703)->primal_0 = v1v0_1;
    (&_S4703)->differential_0 = _S4683;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4704;
    (&_S4704)->primal_0 = v2v0_1;
    (&_S4704)->differential_0 = _S4683;
    s_bwd_prop_cross_0(&_S4703, &_S4704, _S4702);
    float3  _S4705 = _S4685.differential_0 + _S4700.differential_0;
    float3  _S4706 = _S4694.differential_0 + _S4704.differential_0;
    float3  _S4707 = _S4690.differential_0 + _S4703.differential_0;
    float3  _S4708 = - _S4705 + - _S4706 + - _S4707;
    float3  _S4709 = _S4697.differential_0 + _S4701.differential_0;
    dpray_d_4->primal_0 = (*dpray_d_4).primal_0;
    dpray_d_4->differential_0 = _S4709;
    dpray_o_4->primal_0 = (*dpray_o_4).primal_0;
    dpray_o_4->differential_0 = _S4705;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S4664;
    FixedArray<float3 , 3>  _S4710;
    _S4710[int(0)] = _S4683;
    _S4710[int(1)] = _S4683;
    _S4710[int(2)] = _S4683;
    _S4710[int(2)] = _S4706;
    _S4710[int(0)] = _S4708;
    _S4710[int(1)] = _S4707;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S4710;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4711, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4712, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4713, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4714, float _S4715)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S4711, _S4712, _S4713, _S4714, _S4715);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_5, float2  hardness_17, float3  ray_o_6, float3  ray_d_6, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S4716 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4717 = { _S4716, _S4716, _S4716 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_5;
    (&dp_verts_0)->differential_0 = _S4717;
    float2  _S4718 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S4718;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_6;
    (&dp_ray_o_2)->differential_0 = _S4716;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_6;
    (&dp_ray_d_2)->differential_0 = _S4716;
    s_bwd_evaluate_alpha_opaque_triangle_0(&dp_verts_0, &dp_hardness_2, &dp_ray_o_2, &dp_ray_d_2, v_alpha_3);
    *v_verts_2 = (&dp_verts_0)->differential_0;
    *v_hardness_6 = dp_hardness_2.differential_0;
    *v_ray_o_3 = dp_ray_o_2.differential_0;
    *v_ray_d_3 = dp_ray_d_2.differential_0;
    return;
}

inline __device__ void evaluate_color_opaque_triangle(FixedArray<float3 , 3>  * verts_6, FixedArray<float3 , 3>  * rgbs_4, float3  ray_o_7, float3  ray_d_7, float3  * color_13, float * depth_20)
{
    float3  v1v0_2 = (*verts_6)[int(1)] - (*verts_6)[int(0)];
    float3  v2v0_2 = (*verts_6)[int(2)] - (*verts_6)[int(0)];
    float3  rov0_2 = ray_o_7 - (*verts_6)[int(0)];
    float3  n_1 = cross_0(v1v0_2, v2v0_2);
    float3  q_3 = cross_0(rov0_2, ray_d_7);
    float d_2 = 1.0f / dot_0(ray_d_7, n_1);
    float u_30 = d_2 * dot_0(- q_3, v2v0_2);
    float v_30 = d_2 * dot_0(q_3, v1v0_2);
    *depth_20 = d_2 * dot_0(- n_1, rov0_2);
    *color_13 = (*rgbs_4)[int(0)] * make_float3 (1.0f - u_30 - v_30) + (*rgbs_4)[int(1)] * make_float3 (u_30) + (*rgbs_4)[int(2)] * make_float3 (v_30);
    *depth_20 = (F32_log(((F32_max((*depth_20), (9.999999960041972e-13f))))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float3  dpcolor_1, float dpdepth_2)
{
    float3  v1v0_3 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_3 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_3 = (*dpray_o_5).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S4719 = s_primal_ctx_cross_0(v1v0_3, v2v0_3);
    float3  _S4720 = s_primal_ctx_cross_0(rov0_3, (*dpray_d_5).primal_0);
    float _S4721 = s_primal_ctx_dot_0((*dpray_d_5).primal_0, _S4719);
    float d_3 = 1.0f / _S4721;
    float _S4722 = _S4721 * _S4721;
    float3  _S4723 = - _S4720;
    float _S4724 = s_primal_ctx_dot_0(_S4723, v2v0_3);
    float u_31 = d_3 * _S4724;
    float3  _S4725 = make_float3 (u_31);
    float _S4726 = s_primal_ctx_dot_0(_S4720, v1v0_3);
    float v_31 = d_3 * _S4726;
    float3  _S4727 = make_float3 (v_31);
    float3  _S4728 = - _S4719;
    float _S4729 = s_primal_ctx_dot_0(_S4728, rov0_3);
    float _S4730 = d_3 * _S4729;
    float3  _S4731 = make_float3 (1.0f - u_31 - v_31);
    DiffPair_float_0 _S4732;
    (&_S4732)->primal_0 = s_primal_ctx_max_0(_S4730, 9.999999960041972e-13f);
    (&_S4732)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4732, dpdepth_2);
    DiffPair_float_0 _S4733;
    (&_S4733)->primal_0 = _S4730;
    (&_S4733)->differential_0 = 0.0f;
    DiffPair_float_0 _S4734;
    (&_S4734)->primal_0 = 9.999999960041972e-13f;
    (&_S4734)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4733, &_S4734, _S4732.differential_0);
    float3  _S4735 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S4736 = _S4727 * dpcolor_1;
    float3  _S4737 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S4738 = _S4725 * dpcolor_1;
    float3  _S4739 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S4740 = _S4731 * dpcolor_1;
    float _S4741 = - (_S4739.x + _S4739.y + _S4739.z);
    float _S4742 = d_3 * _S4733.differential_0;
    float _S4743 = _S4729 * _S4733.differential_0;
    float3  _S4744 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4745;
    (&_S4745)->primal_0 = _S4728;
    (&_S4745)->differential_0 = _S4744;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4746;
    (&_S4746)->primal_0 = rov0_3;
    (&_S4746)->differential_0 = _S4744;
    s_bwd_prop_dot_0(&_S4745, &_S4746, _S4742);
    float3  _S4747 = - _S4745.differential_0;
    float _S4748 = _S4741 + _S4735.x + _S4735.y + _S4735.z;
    float _S4749 = d_3 * _S4748;
    float _S4750 = _S4726 * _S4748;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4751;
    (&_S4751)->primal_0 = _S4720;
    (&_S4751)->differential_0 = _S4744;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4752;
    (&_S4752)->primal_0 = v1v0_3;
    (&_S4752)->differential_0 = _S4744;
    s_bwd_prop_dot_0(&_S4751, &_S4752, _S4749);
    float _S4753 = _S4741 + _S4737.x + _S4737.y + _S4737.z;
    float _S4754 = d_3 * _S4753;
    float _S4755 = _S4724 * _S4753;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4756;
    (&_S4756)->primal_0 = _S4723;
    (&_S4756)->differential_0 = _S4744;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4757;
    (&_S4757)->primal_0 = v2v0_3;
    (&_S4757)->differential_0 = _S4744;
    s_bwd_prop_dot_0(&_S4756, &_S4757, _S4754);
    float3  _S4758 = - _S4756.differential_0;
    float _S4759 = - ((_S4743 + _S4750 + _S4755) / _S4722);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4760;
    (&_S4760)->primal_0 = (*dpray_d_5).primal_0;
    (&_S4760)->differential_0 = _S4744;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4761;
    (&_S4761)->primal_0 = _S4719;
    (&_S4761)->differential_0 = _S4744;
    s_bwd_prop_dot_0(&_S4760, &_S4761, _S4759);
    float3  _S4762 = _S4751.differential_0 + _S4758;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4763;
    (&_S4763)->primal_0 = rov0_3;
    (&_S4763)->differential_0 = _S4744;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4764;
    (&_S4764)->primal_0 = (*dpray_d_5).primal_0;
    (&_S4764)->differential_0 = _S4744;
    s_bwd_prop_cross_0(&_S4763, &_S4764, _S4762);
    float3  _S4765 = _S4747 + _S4761.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4766;
    (&_S4766)->primal_0 = v1v0_3;
    (&_S4766)->differential_0 = _S4744;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4767;
    (&_S4767)->primal_0 = v2v0_3;
    (&_S4767)->differential_0 = _S4744;
    s_bwd_prop_cross_0(&_S4766, &_S4767, _S4765);
    float3  _S4768 = _S4746.differential_0 + _S4763.differential_0;
    float3  _S4769 = _S4757.differential_0 + _S4767.differential_0;
    float3  _S4770 = _S4752.differential_0 + _S4766.differential_0;
    float3  _S4771 = - _S4768 + - _S4769 + - _S4770;
    float3  _S4772 = _S4760.differential_0 + _S4764.differential_0;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S4772;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S4768;
    FixedArray<float3 , 3>  _S4773;
    _S4773[int(0)] = _S4744;
    _S4773[int(1)] = _S4744;
    _S4773[int(2)] = _S4744;
    _S4773[int(2)] = _S4736;
    _S4773[int(1)] = _S4738;
    _S4773[int(0)] = _S4740;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S4773;
    FixedArray<float3 , 3>  _S4774;
    _S4774[int(0)] = _S4744;
    _S4774[int(1)] = _S4744;
    _S4774[int(2)] = _S4744;
    _S4774[int(2)] = _S4769;
    _S4774[int(0)] = _S4771;
    _S4774[int(1)] = _S4770;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S4774;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4775, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4776, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4777, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4778, float3  _S4779, float _S4780)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S4775, _S4776, _S4777, _S4778, _S4779, _S4780);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_7, FixedArray<float3 , 3>  * rgbs_5, float3  ray_o_8, float3  ray_d_8, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S4781 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4782 = { _S4781, _S4781, _S4781 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_7;
    (&dp_verts_1)->differential_0 = _S4782;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_5;
    (&dp_rgbs_0)->differential_0 = _S4782;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_8;
    (&dp_ray_o_3)->differential_0 = _S4781;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_8;
    (&dp_ray_d_3)->differential_0 = _S4781;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

