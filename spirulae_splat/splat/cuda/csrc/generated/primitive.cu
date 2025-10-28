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
        *depth_0 = length_1(mean_c_0);
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
        *depth_1 = _S254;
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
        *depth_2 = length_1(mean_c_2);
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
        *depth_3 = length_1(mean_c_3);
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
        *depth_4 = _S299;
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
    *depth_5 = length_1(mean_c_5);
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
    *depth_6 = length_1(mean_c_6);
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
    *depth_7 = length_1(mean_c_7);
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
    *depth_8 = length_1(mean_c_8);
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
    *depth_9 = length_1(mean_c_9);
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

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S394, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S395, float3  _S396)
{
    _d_max_vector_0(_S394, _S395, _S396);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S397, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S398, float3  _S399)
{
    _d_mul_0(_S397, _S398, _S399);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S400, DiffPair_float_0 * _S401, float _S402)
{
    _d_min_0(_S400, _S401, _S402);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S403, float _S404)
{
    _d_log_0(_S403, _S404);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S405, float _S406)
{
    _d_exp_0(_S405, _S406);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S407, DiffPair_float_0 * _S408, float _S409)
{
    _d_max_0(_S407, _S408, _S409);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S410, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S411, Matrix<float, 2, 2>  _S412)
{
    mul_3(_S410, _S411, _S412);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S413, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S414, Matrix<float, 2, 3>  _S415)
{
    mul_2(_S413, _S414, _S415);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S416, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S417, Matrix<float, 3, 3>  _S418)
{
    mul_1(_S416, _S417, _S418);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S419, float3  _S420)
{
    _d_exp_vector_0(_S419, _S420);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_14, float fy_14, float cx_14, float cy_14, float4  radial_coeffs_17, float2  tangential_coeffs_17, float2  thin_prism_coeffs_17, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_10 = s_primal_ctx_mul_1(R_14, mean_10) + t_13;
    float3  _S421 = s_primal_ctx_exp_0(scale_12);
    float _S422 = quat_13.y;
    float x2_13 = _S422 * _S422;
    float y2_13 = quat_13.z * quat_13.z;
    float z2_23 = quat_13.w * quat_13.w;
    float xy_13 = quat_13.y * quat_13.z;
    float xz_13 = quat_13.y * quat_13.w;
    float yz_13 = quat_13.z * quat_13.w;
    float wx_13 = quat_13.x * quat_13.y;
    float wy_13 = quat_13.x * quat_13.z;
    float wz_13 = quat_13.x * quat_13.w;
    Matrix<float, 3, 3>  _S423 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S421.x, 0.0f, 0.0f, 0.0f, _S421.y, 0.0f, 0.0f, 0.0f, _S421.z);
    Matrix<float, 3, 3>  _S424 = s_primal_ctx_mul_2(_S423, S_0);
    Matrix<float, 3, 3>  _S425 = transpose_0(_S424);
    Matrix<float, 3, 3>  _S426 = s_primal_ctx_mul_2(_S424, _S425);
    Matrix<float, 3, 3>  _S427 = s_primal_ctx_mul_2(R_14, _S426);
    Matrix<float, 3, 3>  _S428 = transpose_0(R_14);
    Matrix<float, 3, 3>  _S429 = s_primal_ctx_mul_2(_S427, _S428);
    float _S430 = float(image_width_10);
    float _S431 = float(image_height_10);
    float _S432 = 0.30000001192092896f * (0.5f * _S430 / fx_14);
    float lim_x_pos_0 = (_S430 - cx_14) / fx_14 + _S432;
    float _S433 = 0.30000001192092896f * (0.5f * _S431 / fy_14);
    float lim_y_pos_0 = (_S431 - cy_14) / fy_14 + _S433;
    float rz_6 = 1.0f / mean_c_10.z;
    float _S434 = mean_c_10.z * mean_c_10.z;
    float rz2_6 = rz_6 * rz_6;
    float _S435 = - (cx_14 / fx_14 + _S432);
    float _S436 = mean_c_10.x * rz_6;
    float _S437 = s_primal_ctx_max_0(_S435, _S436);
    float _S438 = s_primal_ctx_min_0(lim_x_pos_0, _S437);
    float _S439 = - (cy_14 / fy_14 + _S433);
    float _S440 = mean_c_10.y * rz_6;
    float _S441 = s_primal_ctx_max_0(_S439, _S440);
    float _S442 = s_primal_ctx_min_0(lim_y_pos_0, _S441);
    float _S443 = - fx_14;
    float _S444 = _S443 * (mean_c_10.z * _S438);
    float _S445 = - fy_14;
    float _S446 = _S445 * (mean_c_10.z * _S442);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_14 * rz_6, 0.0f, _S444 * rz2_6, 0.0f, fy_14 * rz_6, _S446 * rz2_6);
    Matrix<float, 2, 3>  _S447 = s_primal_ctx_mul_3(J_11, _S429);
    Matrix<float, 3, 2>  _S448 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S449 = s_primal_ctx_mul_4(_S447, _S448);
    float _S450 = fx_14 * mean_c_10.x;
    float _S451 = fy_14 * mean_c_10.y;
    float _S452 = _S449.rows[int(0)].y * _S449.rows[int(1)].x;
    float det_orig_11 = _S449.rows[int(0)].x * _S449.rows[int(1)].y - _S452;
    float _S453 = _S449.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S454 = _S449;
    *&(((&_S454)->rows + (int(0)))->x) = _S453;
    float _S455 = _S449.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S454)->rows + (int(1)))->y) = _S455;
    Matrix<float, 2, 2>  _S456 = _S454;
    Matrix<float, 2, 2>  _S457 = _S454;
    float det_blur_6 = _S453 * _S455 - _S452;
    float _S458 = det_orig_11 / det_blur_6;
    float _S459 = det_blur_6 * det_blur_6;
    float _S460 = s_primal_ctx_max_0(0.0f, _S458);
    float _S461 = s_primal_ctx_sqrt_0(_S460);
    float invdet_7 = 1.0f / det_blur_6;
    float _S462 = - _S449.rows[int(0)].y;
    float _S463 = - _S449.rows[int(1)].x;
    float _S464 = - in_opacity_10;
    float _S465 = 1.0f + s_primal_ctx_exp_1(_S464);
    float _S466 = 1.0f / _S465;
    float _S467 = _S465 * _S465;
    float _S468;
    if(antialiased_10)
    {
        _S468 = _S466 * _S461;
    }
    else
    {
        _S468 = _S466;
    }
    float _S469 = _S468 / 0.00392156885936856f;
    float _S470 = 2.0f * s_primal_ctx_log_0(_S469);
    float _S471 = s_primal_ctx_sqrt_0(_S470);
    float _S472 = _S456.rows[int(0)].x;
    float _S473 = _S457.rows[int(1)].y;
    float3  _S474 = mean_10 - - s_primal_ctx_mul_1(_S428, t_13);
    float _S475 = _S474.x;
    float _S476 = _S474.y;
    float _S477 = _S474.z;
    float _S478 = _S475 * _S475 + _S476 * _S476 + _S477 * _S477;
    float _S479 = s_primal_ctx_sqrt_0(_S478);
    float x_34 = _S475 / _S479;
    float3  _S480 = make_float3 (x_34);
    float _S481 = _S479 * _S479;
    float y_13 = _S476 / _S479;
    float z_10 = _S477 / _S479;
    float3  _S482 = make_float3 (z_10);
    float _S483 = - y_13;
    float3  _S484 = make_float3 (_S483);
    float z2_24 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_34 * x_34 - y_13 * y_13;
    float _S485 = 2.0f * x_34;
    float fS1_10 = _S485 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_24 - 0.31539157032966614f;
    float3  _S486 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_34;
    float3  _S487 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S488 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S489 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S490 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S491 = 1.86588168144226074f * z2_24 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S491;
    float3  _S492 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_34;
    float3  _S493 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S494 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S495 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S496 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_34 * fC1_10 - y_13 * fS1_10);
    float3  _S497 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_34 * fS1_10 + y_13 * fC1_10);
    float3  _S498 = make_float3 (pSH9_0);
    float3  _S499 = make_float3 (0.0f);
    float3  _S500 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S501;
    (&_S501)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S483) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_34) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S501)->differential_0 = _S500;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S502;
    (&_S502)->primal_0 = _S499;
    (&_S502)->differential_0 = _S500;
    s_bwd_prop_max_0(&_S501, &_S502, v_rgb_0);
    float3  _S503 = _S497 * _S501.differential_0;
    float3  _S504 = (*sh_coeffs_10)[int(15)] * _S501.differential_0;
    float3  _S505 = _S495 * _S501.differential_0;
    float3  _S506 = (*sh_coeffs_10)[int(14)] * _S501.differential_0;
    float3  _S507 = _S493 * _S501.differential_0;
    float3  _S508 = (*sh_coeffs_10)[int(13)] * _S501.differential_0;
    float3  _S509 = _S492 * _S501.differential_0;
    float3  _S510 = (*sh_coeffs_10)[int(12)] * _S501.differential_0;
    float3  _S511 = _S494 * _S501.differential_0;
    float3  _S512 = (*sh_coeffs_10)[int(11)] * _S501.differential_0;
    float3  _S513 = _S496 * _S501.differential_0;
    float3  _S514 = (*sh_coeffs_10)[int(10)] * _S501.differential_0;
    float3  _S515 = _S498 * _S501.differential_0;
    float3  _S516 = (*sh_coeffs_10)[int(9)] * _S501.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S516.x + _S516.y + _S516.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S504.x + _S504.y + _S504.z);
    float _S517 = _S514.x + _S514.y + _S514.z;
    float _S518 = _S506.x + _S506.y + _S506.z;
    float _S519 = _S512.x + _S512.y + _S512.z;
    float _S520 = _S508.x + _S508.y + _S508.z;
    float _S521 = _S510.x + _S510.y + _S510.z;
    float _S522 = - s_diff_fC2_T_0;
    float3  _S523 = _S489 * _S501.differential_0;
    float3  _S524 = (*sh_coeffs_10)[int(8)] * _S501.differential_0;
    float3  _S525 = _S487 * _S501.differential_0;
    float3  _S526 = (*sh_coeffs_10)[int(7)] * _S501.differential_0;
    float3  _S527 = _S486 * _S501.differential_0;
    float3  _S528 = (*sh_coeffs_10)[int(6)] * _S501.differential_0;
    float3  _S529 = _S488 * _S501.differential_0;
    float3  _S530 = (*sh_coeffs_10)[int(5)] * _S501.differential_0;
    float3  _S531 = _S490 * _S501.differential_0;
    float3  _S532 = (*sh_coeffs_10)[int(4)] * _S501.differential_0;
    float _S533 = _S530.x + _S530.y + _S530.z;
    float _S534 = _S526.x + _S526.y + _S526.z;
    float _S535 = fTmp1B_10 * _S517 + x_34 * s_diff_fS2_T_0 + y_13 * _S522 + 0.54627424478530884f * (_S532.x + _S532.y + _S532.z);
    float _S536 = fTmp1B_10 * _S518 + y_13 * s_diff_fS2_T_0 + x_34 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S524.x + _S524.y + _S524.z);
    float _S537 = y_13 * - _S536;
    float _S538 = x_34 * _S536;
    float _S539 = z_10 * (1.86588168144226074f * (z_10 * _S521) + -2.28522896766662598f * (y_13 * _S519 + x_34 * _S520) + 0.94617468118667603f * (_S528.x + _S528.y + _S528.z));
    float3  _S540 = make_float3 (0.48860251903533936f) * _S501.differential_0;
    float3  _S541 = - _S540;
    float3  _S542 = _S480 * _S541;
    float3  _S543 = (*sh_coeffs_10)[int(3)] * _S541;
    float3  _S544 = _S482 * _S540;
    float3  _S545 = (*sh_coeffs_10)[int(2)] * _S540;
    float3  _S546 = _S484 * _S540;
    float3  _S547 = (*sh_coeffs_10)[int(1)] * _S540;
    float _S548 = (_S491 * _S521 + 1.44530570507049561f * (fS1_10 * _S517 + fC1_10 * _S518) + -1.09254848957061768f * (y_13 * _S533 + x_34 * _S534) + _S539 + _S539 + _S545.x + _S545.y + _S545.z) / _S481;
    float _S549 = _S479 * _S548;
    float _S550 = (fTmp0C_10 * _S519 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S522 + fTmp0B_10 * _S533 + _S485 * _S535 + _S537 + _S537 + - (_S547.x + _S547.y + _S547.z)) / _S481;
    float _S551 = _S479 * _S550;
    float _S552 = (fTmp0C_10 * _S520 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S534 + 2.0f * (y_13 * _S535) + _S538 + _S538 + _S543.x + _S543.y + _S543.z) / _S481;
    float _S553 = _S479 * _S552;
    float _S554 = _S477 * - _S548 + _S476 * - _S550 + _S475 * - _S552;
    DiffPair_float_0 _S555;
    (&_S555)->primal_0 = _S478;
    (&_S555)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S555, _S554);
    float _S556 = _S477 * _S555.differential_0;
    float _S557 = _S476 * _S555.differential_0;
    float _S558 = _S475 * _S555.differential_0;
    float3  _S559 = make_float3 (0.282094806432724f) * _S501.differential_0;
    float3  _S560 = make_float3 (_S553 + _S558 + _S558, _S551 + _S557 + _S557, _S549 + _S556 + _S556);
    float3  _S561 = - - _S560;
    Matrix<float, 3, 3>  _S562 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S563;
    (&_S563)->primal_0 = _S428;
    (&_S563)->differential_0 = _S562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S564;
    (&_S564)->primal_0 = t_13;
    (&_S564)->differential_0 = _S500;
    s_bwd_prop_mul_1(&_S563, &_S564, _S561);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S565 = _S563;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S566 = _S564;
    float2  _S567 = make_float2 (0.0f);
    float2  _S568 = _S567;
    *&((&_S568)->y) = v_conic_0.z;
    float2  _S569 = _S567;
    *&((&_S569)->y) = v_conic_0.y;
    *&((&_S569)->x) = v_conic_0.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S570;
    (&_S570)->primal_0 = mean_c_10;
    (&_S570)->differential_0 = _S500;
    s_bwd_length_impl_0(&_S570, v_depth_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S571 = _S570;
    DiffPair_float_0 _S572;
    (&_S572)->primal_0 = _S473;
    (&_S572)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S572, 0.0f);
    DiffPair_float_0 _S573;
    (&_S573)->primal_0 = _S472;
    (&_S573)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S573, 0.0f);
    DiffPair_float_0 _S574;
    (&_S574)->primal_0 = 3.32999992370605469f;
    (&_S574)->differential_0 = 0.0f;
    DiffPair_float_0 _S575;
    (&_S575)->primal_0 = _S471;
    (&_S575)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S574, &_S575, 0.0f);
    DiffPair_float_0 _S576;
    (&_S576)->primal_0 = _S470;
    (&_S576)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S576, _S575.differential_0);
    float _S577 = 2.0f * _S576.differential_0;
    DiffPair_float_0 _S578;
    (&_S578)->primal_0 = _S469;
    (&_S578)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S578, _S577);
    float _S579 = v_opacity_0 + 254.9999847412109375f * _S578.differential_0;
    Matrix<float, 2, 2>  _S580 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S581 = _S580;
    _S581[int(1)] = _S568;
    _S581[int(0)] = _S569;
    Matrix<float, 2, 2>  _S582 = _S581;
    FixedArray<float3 , 16>  _S583;
    _S583[int(0)] = _S500;
    _S583[int(1)] = _S500;
    _S583[int(2)] = _S500;
    _S583[int(3)] = _S500;
    _S583[int(4)] = _S500;
    _S583[int(5)] = _S500;
    _S583[int(6)] = _S500;
    _S583[int(7)] = _S500;
    _S583[int(8)] = _S500;
    _S583[int(9)] = _S500;
    _S583[int(10)] = _S500;
    _S583[int(11)] = _S500;
    _S583[int(12)] = _S500;
    _S583[int(13)] = _S500;
    _S583[int(14)] = _S500;
    _S583[int(15)] = _S500;
    _S583[int(7)] = _S525;
    _S583[int(0)] = _S559;
    _S583[int(1)] = _S546;
    _S583[int(2)] = _S544;
    _S583[int(3)] = _S542;
    _S583[int(4)] = _S531;
    _S583[int(5)] = _S529;
    _S583[int(6)] = _S527;
    _S583[int(15)] = _S503;
    _S583[int(8)] = _S523;
    _S583[int(9)] = _S515;
    _S583[int(10)] = _S513;
    _S583[int(11)] = _S511;
    _S583[int(12)] = _S509;
    _S583[int(13)] = _S507;
    _S583[int(14)] = _S505;
    float3  _S584 = _S583[int(0)];
    float3  _S585 = _S583[int(1)];
    float3  _S586 = _S583[int(2)];
    float3  _S587 = _S583[int(3)];
    float3  _S588 = _S583[int(4)];
    float3  _S589 = _S583[int(5)];
    float3  _S590 = _S583[int(6)];
    float3  _S591 = _S583[int(7)];
    float3  _S592 = _S583[int(8)];
    float3  _S593 = _S583[int(9)];
    float3  _S594 = _S583[int(10)];
    float3  _S595 = _S583[int(11)];
    float3  _S596 = _S583[int(12)];
    float3  _S597 = _S583[int(13)];
    float3  _S598 = _S583[int(14)];
    float3  _S599 = _S583[int(15)];
    float2  _S600 = make_float2 (0.0f, _S572.differential_0);
    float2  _S601 = make_float2 (_S573.differential_0, 0.0f);
    float _S602;
    if(antialiased_10)
    {
        float _S603 = _S466 * _S579;
        _S468 = _S461 * _S579;
        _S602 = _S603;
    }
    else
    {
        _S468 = _S579;
        _S602 = 0.0f;
    }
    float _S604 = - (_S468 / _S467);
    DiffPair_float_0 _S605;
    (&_S605)->primal_0 = _S464;
    (&_S605)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S605, _S604);
    float _S606 = - _S605.differential_0;
    float _S607 = invdet_7 * _S582.rows[int(1)].y;
    float _S608 = - (invdet_7 * _S582.rows[int(1)].x);
    float _S609 = - (invdet_7 * _S582.rows[int(0)].y);
    float _S610 = invdet_7 * _S582.rows[int(0)].x;
    float _S611 = - ((_S453 * _S582.rows[int(1)].y + _S463 * _S582.rows[int(1)].x + _S462 * _S582.rows[int(0)].y + _S455 * _S582.rows[int(0)].x) / _S459);
    DiffPair_float_0 _S612;
    (&_S612)->primal_0 = _S460;
    (&_S612)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S612, _S602);
    DiffPair_float_0 _S613;
    (&_S613)->primal_0 = 0.0f;
    (&_S613)->differential_0 = 0.0f;
    DiffPair_float_0 _S614;
    (&_S614)->primal_0 = _S458;
    (&_S614)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S613, &_S614, _S612.differential_0);
    float _S615 = _S614.differential_0 / _S459;
    float s_diff_det_orig_T_0 = det_blur_6 * _S615;
    float _S616 = _S611 + det_orig_11 * - _S615;
    float _S617 = - _S616;
    float _S618 = _S453 * _S616;
    float _S619 = _S455 * _S616;
    Matrix<float, 2, 2>  _S620 = _S580;
    _S620[int(1)] = _S600;
    _S620[int(0)] = _S601;
    _S454 = _S620;
    *&(((&_S454)->rows + (int(1)))->y) = 0.0f;
    float _S621 = _S610 + _S618 + _S620.rows[int(1)].y;
    *&(((&_S454)->rows + (int(0)))->x) = 0.0f;
    float _S622 = _S607 + _S619 + _S620.rows[int(0)].x;
    float _S623 = _S617 + - s_diff_det_orig_T_0;
    float _S624 = _S608 + _S449.rows[int(0)].y * _S623;
    float _S625 = _S609 + _S449.rows[int(1)].x * _S623;
    float _S626 = _S449.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S627 = _S621 + _S449.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S628 = _S567;
    *&((&_S628)->x) = _S624;
    *&((&_S628)->y) = _S627;
    float _S629 = _S622 + _S626;
    float2  _S630 = _S567;
    *&((&_S630)->y) = _S625;
    *&((&_S630)->x) = _S629;
    float _S631 = _S451 * v_mean2d_0.y;
    float _S632 = fy_14 * (rz_6 * v_mean2d_0.y);
    float _S633 = _S450 * v_mean2d_0.x;
    float _S634 = fx_14 * (rz_6 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S635 = _S580;
    _S635[int(1)] = _S628;
    _S635[int(0)] = _S630;
    Matrix<float, 2, 2>  _S636 = _S454 + _S635;
    Matrix<float, 2, 3>  _S637 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S638;
    (&_S638)->primal_0 = _S447;
    (&_S638)->differential_0 = _S637;
    Matrix<float, 3, 2>  _S639 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S640;
    (&_S640)->primal_0 = _S448;
    (&_S640)->differential_0 = _S639;
    s_bwd_prop_mul_2(&_S638, &_S640, _S636);
    Matrix<float, 2, 3>  _S641 = transpose_2(_S640.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S642;
    (&_S642)->primal_0 = J_11;
    (&_S642)->differential_0 = _S637;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S643;
    (&_S643)->primal_0 = _S429;
    (&_S643)->differential_0 = _S562;
    s_bwd_prop_mul_3(&_S642, &_S643, _S638.differential_0);
    Matrix<float, 2, 3>  _S644 = _S641 + _S642.differential_0;
    float _S645 = _S446 * _S644.rows[int(1)].z;
    float s_diff_ty_T_0 = _S445 * (rz2_6 * _S644.rows[int(1)].z);
    float _S646 = fy_14 * _S644.rows[int(1)].y;
    float _S647 = _S444 * _S644.rows[int(0)].z;
    float s_diff_tx_T_0 = _S443 * (rz2_6 * _S644.rows[int(0)].z);
    float _S648 = fx_14 * _S644.rows[int(0)].x;
    float _S649 = mean_c_10.z * s_diff_ty_T_0;
    float _S650 = _S442 * s_diff_ty_T_0;
    DiffPair_float_0 _S651;
    (&_S651)->primal_0 = lim_y_pos_0;
    (&_S651)->differential_0 = 0.0f;
    DiffPair_float_0 _S652;
    (&_S652)->primal_0 = _S441;
    (&_S652)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S651, &_S652, _S649);
    DiffPair_float_0 _S653;
    (&_S653)->primal_0 = _S439;
    (&_S653)->differential_0 = 0.0f;
    DiffPair_float_0 _S654;
    (&_S654)->primal_0 = _S440;
    (&_S654)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S653, &_S654, _S652.differential_0);
    float _S655 = mean_c_10.y * _S654.differential_0;
    float _S656 = rz_6 * _S654.differential_0;
    float _S657 = mean_c_10.z * s_diff_tx_T_0;
    float _S658 = _S438 * s_diff_tx_T_0;
    DiffPair_float_0 _S659;
    (&_S659)->primal_0 = lim_x_pos_0;
    (&_S659)->differential_0 = 0.0f;
    DiffPair_float_0 _S660;
    (&_S660)->primal_0 = _S437;
    (&_S660)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S659, &_S660, _S657);
    DiffPair_float_0 _S661;
    (&_S661)->primal_0 = _S435;
    (&_S661)->differential_0 = 0.0f;
    DiffPair_float_0 _S662;
    (&_S662)->primal_0 = _S436;
    (&_S662)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S661, &_S662, _S660.differential_0);
    float _S663 = rz_6 * (_S645 + _S647);
    float _S664 = _S650 + _S658 + - ((_S631 + _S633 + _S646 + _S648 + _S655 + mean_c_10.x * _S662.differential_0 + _S663 + _S663) / _S434);
    float _S665 = _S632 + _S656;
    float _S666 = _S634 + rz_6 * _S662.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S667;
    (&_S667)->primal_0 = _S427;
    (&_S667)->differential_0 = _S562;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S668;
    (&_S668)->primal_0 = _S428;
    (&_S668)->differential_0 = _S562;
    s_bwd_prop_mul_4(&_S667, &_S668, _S643.differential_0);
    Matrix<float, 3, 3>  _S669 = transpose_0(_S668.differential_0 + _S565.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S670;
    (&_S670)->primal_0 = R_14;
    (&_S670)->differential_0 = _S562;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S671;
    (&_S671)->primal_0 = _S426;
    (&_S671)->differential_0 = _S562;
    s_bwd_prop_mul_4(&_S670, &_S671, _S667.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S672;
    (&_S672)->primal_0 = _S424;
    (&_S672)->differential_0 = _S562;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S673;
    (&_S673)->primal_0 = _S425;
    (&_S673)->differential_0 = _S562;
    s_bwd_prop_mul_4(&_S672, &_S673, _S671.differential_0);
    Matrix<float, 3, 3>  _S674 = _S672.differential_0 + transpose_0(_S673.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S675;
    (&_S675)->primal_0 = _S423;
    (&_S675)->differential_0 = _S562;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S676;
    (&_S676)->primal_0 = S_0;
    (&_S676)->differential_0 = _S562;
    s_bwd_prop_mul_4(&_S675, &_S676, _S674);
    Matrix<float, 3, 3>  _S677 = transpose_0(_S675.differential_0);
    float _S678 = 2.0f * - _S677.rows[int(2)].z;
    float _S679 = 2.0f * _S677.rows[int(2)].y;
    float _S680 = 2.0f * _S677.rows[int(2)].x;
    float _S681 = 2.0f * _S677.rows[int(1)].z;
    float _S682 = 2.0f * - _S677.rows[int(1)].y;
    float _S683 = 2.0f * _S677.rows[int(1)].x;
    float _S684 = 2.0f * _S677.rows[int(0)].z;
    float _S685 = 2.0f * _S677.rows[int(0)].y;
    float _S686 = 2.0f * - _S677.rows[int(0)].x;
    float _S687 = - _S683 + _S685;
    float _S688 = _S680 + - _S684;
    float _S689 = - _S679 + _S681;
    float _S690 = _S679 + _S681;
    float _S691 = _S680 + _S684;
    float _S692 = _S683 + _S685;
    float _S693 = quat_13.w * (_S682 + _S686);
    float _S694 = quat_13.z * (_S678 + _S686);
    float _S695 = quat_13.y * (_S678 + _S682);
    float _S696 = quat_13.x * _S687 + quat_13.z * _S690 + quat_13.y * _S691 + _S693 + _S693;
    float _S697 = quat_13.x * _S688 + quat_13.w * _S690 + quat_13.y * _S692 + _S694 + _S694;
    float _S698 = quat_13.x * _S689 + quat_13.w * _S691 + quat_13.z * _S692 + _S695 + _S695;
    float _S699 = quat_13.w * _S687 + quat_13.z * _S688 + quat_13.y * _S689;
    float3  _S700 = _S500;
    *&((&_S700)->z) = _S676.differential_0.rows[int(2)].z;
    *&((&_S700)->y) = _S676.differential_0.rows[int(1)].y;
    *&((&_S700)->x) = _S676.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S701;
    (&_S701)->primal_0 = scale_12;
    (&_S701)->differential_0 = _S500;
    s_bwd_prop_exp_1(&_S701, _S700);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S702 = _S701;
    float3  _S703 = _S500;
    *&((&_S703)->z) = _S664;
    *&((&_S703)->y) = _S665;
    *&((&_S703)->x) = _S666;
    float3  _S704 = _S571.differential_0 + _S703;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S705;
    (&_S705)->primal_0 = R_14;
    (&_S705)->differential_0 = _S562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S706;
    (&_S706)->primal_0 = mean_10;
    (&_S706)->differential_0 = _S500;
    s_bwd_prop_mul_1(&_S705, &_S706, _S704);
    float3  _S707 = _S704 + _S566.differential_0;
    Matrix<float, 3, 3>  _S708 = _S669 + _S670.differential_0 + _S705.differential_0;
    float4  _S709 = make_float4 (0.0f);
    *&((&_S709)->w) = _S696;
    *&((&_S709)->z) = _S697;
    *&((&_S709)->y) = _S698;
    *&((&_S709)->x) = _S699;
    float4  _S710 = _S709;
    float3  _S711 = _S706.differential_0 + _S560;
    *v_mean_0 = _S711;
    *v_quat_0 = _S710;
    *v_scale_0 = _S702.differential_0;
    *v_in_opacity_0 = _S606;
    (*v_sh_coeffs_0)[int(0)] = _S584;
    (*v_sh_coeffs_0)[int(1)] = _S585;
    (*v_sh_coeffs_0)[int(2)] = _S586;
    (*v_sh_coeffs_0)[int(3)] = _S587;
    (*v_sh_coeffs_0)[int(4)] = _S588;
    (*v_sh_coeffs_0)[int(5)] = _S589;
    (*v_sh_coeffs_0)[int(6)] = _S590;
    (*v_sh_coeffs_0)[int(7)] = _S591;
    (*v_sh_coeffs_0)[int(8)] = _S592;
    (*v_sh_coeffs_0)[int(9)] = _S593;
    (*v_sh_coeffs_0)[int(10)] = _S594;
    (*v_sh_coeffs_0)[int(11)] = _S595;
    (*v_sh_coeffs_0)[int(12)] = _S596;
    (*v_sh_coeffs_0)[int(13)] = _S597;
    (*v_sh_coeffs_0)[int(14)] = _S598;
    (*v_sh_coeffs_0)[int(15)] = _S599;
    *v_R_1 = _S708;
    *v_t_1 = _S707;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S712;
    DiffPair_float_0 _S713;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S714;
    DiffPair_float_0 _S715;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S716;
    DiffPair_float_0 _S717;
    DiffPair_float_0 _S718;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S719;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S720;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S721;
};

inline __device__ CameraDistortion_0 s_primal_ctx_CameraDistortion_x24init_0(float4  dpradial_coeffs_0, float2  dptangential_coeffs_0, float2  dpthin_prism_coeffs_0)
{
    CameraDistortion_0 _S722 = { dpradial_coeffs_0, dptangential_coeffs_0, dpthin_prism_coeffs_0 };
    return _S722;
}

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S723, float _S724)
{
    return s_primal_ctx_atan2_0(_S723, _S724);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S725;
    DiffPair_float_0 _S726;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S727 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S725 = _S727;
    _s_diff_ctx_2->_S726 = _S727;
    (&_s_diff_ctx_2->_S725)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S725)->differential_0 = 0.0f;
    (&_s_diff_ctx_2->_S726)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S726)->differential_0 = 0.0f;
    DiffPair_float_0 _S728 = *dpdpy_0;
    _s_diff_ctx_2->_S725 = *dpdpy_0;
    DiffPair_float_0 _S729 = *dpdpx_0;
    _s_diff_ctx_2->_S726 = *dpdpx_0;
    float _S730 = _S729.primal_0 * _S729.primal_0 + _S728.primal_0 * _S728.primal_0;
    float _S731 = - _S728.primal_0 / _S730 * dpdOut_0;
    float _S732 = _S729.primal_0 / _S730 * dpdOut_0;
    dpdpy_0->primal_0 = _S728.primal_0;
    dpdpy_0->differential_0 = _S732;
    dpdpx_0->primal_0 = _S729.primal_0;
    dpdpx_0->differential_0 = _S731;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S733, DiffPair_float_0 * _S734, float _S735, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S736 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S712 = _S736;
    _s_diff_ctx_3->_S713 = _S736;
    (&_s_diff_ctx_3->_S712)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S712)->differential_0 = 0.0f;
    (&_s_diff_ctx_3->_S713)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S713)->differential_0 = 0.0f;
    DiffPair_float_0 _S737 = *_S733;
    _s_diff_ctx_3->_S712 = *_S733;
    DiffPair_float_0 _S738 = *_S734;
    _s_diff_ctx_3->_S713 = *_S734;
    DiffPair_float_0 _S739 = _S737;
    DiffPair_float_0 _S740 = _S738;
    s_bwd_prop_d_atan2_Intermediates_0 _S741;
    (&_S741)->_S725 = _S736;
    (&_S741)->_S726 = _S736;
    s_primal_ctx_d_atan2_0(&_S739, &_S740, _S735, &_S741);
    *_S733 = _S739;
    *_S734 = _S740;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S742;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S743;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S744;
    DiffPair_float_0 _S745;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S746;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S747;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S748 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S747 = _S748;
    (&_s_diff_ctx_4->_S747)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S747)->differential_0 = 0.0f;
    DiffPair_float_0 _S749 = *dpdpx_1;
    _s_diff_ctx_4->_S747 = *dpdpx_1;
    float _S750 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S749.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S749.primal_0;
    dpdpx_1->differential_0 = _S750;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S751, float _S752, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S753 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S743 = _S753;
    (&_s_diff_ctx_5->_S743)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S743)->differential_0 = 0.0f;
    DiffPair_float_0 _S754 = *_S751;
    _s_diff_ctx_5->_S743 = *_S751;
    DiffPair_float_0 _S755 = _S754;
    s_bwd_prop_d_sqrt_Intermediates_0 _S756;
    (&_S756)->_S747 = _S753;
    s_primal_ctx_d_sqrt_0(&_S755, _S752, &_S756);
    *_S751 = _S755;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_6)
{
    float2  _S757 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S758 = { _S757, _S757 };
    DiffPair_float_0 _S759 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S760 = { _S759 };
    _s_diff_ctx_6->_S744 = _S758;
    _s_diff_ctx_6->_S745 = _S759;
    _s_diff_ctx_6->_S746 = _S760;
    (&_s_diff_ctx_6->_S744)->primal_0 = _S757;
    (&_s_diff_ctx_6->_S744)->differential_0 = _S757;
    (&_s_diff_ctx_6->_S745)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S745)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S761 = *dpdpx_2;
    _s_diff_ctx_6->_S744 = *dpdpx_2;
    float _S762 = _S761.primal_0.x;
    float _S763 = _S761.primal_0.y;
    DiffPair_float_0 _S764;
    (&_S764)->primal_0 = _S762 * _S762 + _S763 * _S763;
    (&_S764)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S764, dp_s_dOut_0, &_s_diff_ctx_6->_S746);
    _s_diff_ctx_6->_S745 = _S764;
    float _S765 = _S761.primal_0.y * _S764.differential_0;
    float _S766 = _S765 + _S765;
    float _S767 = _S761.primal_0.x * _S764.differential_0;
    float _S768 = _S767 + _S767;
    float2  _S769 = _S757;
    *&((&_S769)->y) = _S766;
    *&((&_S769)->x) = _S768;
    dpdpx_2->primal_0 = _S761.primal_0;
    dpdpx_2->differential_0 = _S769;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S770, float _S771, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_7)
{
    float2  _S772 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S773 = { _S772, _S772 };
    _s_diff_ctx_7->_S742 = _S773;
    (&_s_diff_ctx_7->_S742)->primal_0 = _S772;
    (&_s_diff_ctx_7->_S742)->differential_0 = _S772;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S774 = *_S770;
    _s_diff_ctx_7->_S742 = *_S770;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S775 = _S774;
    DiffPair_float_0 _S776 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S777 = { _S776 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S778;
    (&_S778)->_S744 = _S773;
    (&_S778)->_S745 = _S776;
    (&_S778)->_S746 = _S777;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S775, _S771, &_S778);
    *_S770 = _S775;
    return;
}

inline __device__ void s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float dpfx_0, float dpfy_0, float dpcx_0, float dpcy_0, float4  dpradial_coeffs_1, float2  dptangential_coeffs_1, float2  dpthin_prism_coeffs_1, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_float_0 _S779 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S780 = { _S779, _S779 };
    _s_diff_ctx_8->_S714 = _S779;
    _s_diff_ctx_8->_S715 = _S779;
    _s_diff_ctx_8->_S716 = _S780;
    _s_diff_ctx_8->_S717 = _S779;
    _s_diff_ctx_8->_S718 = _S779;
    _s_diff_ctx_8->_S719 = _S780;
    (&_s_diff_ctx_8->_S714)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S714)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S715)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S715)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S717)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S717)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S718)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S718)->differential_0 = 0.0f;
    float2  _S781 = make_float2 (0.0f);
    CameraDistortion_0 _S782 = s_primal_ctx_CameraDistortion_x24init_0(dpradial_coeffs_1, dptangential_coeffs_1, dpthin_prism_coeffs_1);
    float2  _S783 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S784 = length_0(_S783);
    float _S785 = dpmean3d_0.z;
    float _S786 = s_primal_ctx_atan2_0(_S784, _S785);
    float k_1;
    if(_S786 < 0.00100000004749745f)
    {
        k_1 = (1.0f - _S786 * _S786 / 3.0f) / _S785;
    }
    else
    {
        k_1 = _S786 / _S784;
    }
    float2  _S787 = _S783 * make_float2 (k_1);
    float k1_1 = _S782.radial_coeffs_0.x;
    float k2_1 = _S782.radial_coeffs_0.y;
    float k3_1 = _S782.radial_coeffs_0.z;
    float k4_2 = _S782.radial_coeffs_0.w;
    float p1_2 = _S782.tangential_coeffs_0.x;
    float p2_2 = _S782.tangential_coeffs_0.y;
    float sx1_2 = _S782.thin_prism_coeffs_0.x;
    float sy1_2 = _S782.thin_prism_coeffs_0.y;
    float u_4 = _S787.x;
    float v_4 = _S787.y;
    float r2_4 = u_4 * u_4 + v_4 * v_4;
    float _S788 = 2.0f * p1_2;
    float _S789 = 2.0f * p2_2;
    float2  _S790 = _S787 * make_float2 (1.0f + r2_4 * (k1_1 + r2_4 * (k2_1 + r2_4 * (k3_1 + r2_4 * k4_2)))) + make_float2 (_S788 * u_4 * v_4 + p2_2 * (r2_4 + 2.0f * u_4 * u_4) + sx1_2 * r2_4, _S789 * u_4 * v_4 + p1_2 * (r2_4 + 2.0f * v_4 * v_4) + sy1_2 * r2_4);
    float2  _S791 = make_float2 (dpfx_0 * _S790.x + dpcx_0, dpfy_0 * _S790.y + dpcy_0);
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (0.0f);
    float _S792 = s_primal_ctx_s_primal_ctx_atan2_0(_S784, _S785);
    bool _S793 = _S792 < 0.00100000004749745f;
    float _S794;
    float _S795;
    float _S796;
    if(_S793)
    {
        float _S797 = 1.0f - _S792 * _S792 / 3.0f;
        float _S798 = _S785 * _S785;
        k_1 = _S797 / _S785;
        _S794 = 0.0f;
        _S795 = _S798;
        _S796 = _S797;
    }
    else
    {
        float _S799 = _S784 * _S784;
        k_1 = _S792 / _S784;
        _S794 = _S799;
        _S795 = 0.0f;
        _S796 = 0.0f;
    }
    float2  _S800 = make_float2 (k_1);
    float2  _S801 = _S783 * make_float2 (k_1);
    float u_5 = _S801.x;
    float v_5 = _S801.y;
    float r2_5 = u_5 * u_5 + v_5 * v_5;
    float _S802 = k3_1 + r2_5 * k4_2;
    float _S803 = k2_1 + r2_5 * _S802;
    float _S804 = k1_1 + r2_5 * _S803;
    float2  _S805 = make_float2 (dpfx_0, 0.0f);
    float2  _S806 = _S801 * _S805;
    float _S807 = p2_2 * dpfx_0;
    float _S808 = _S806.x + _S806.y;
    float _S809 = r2_5 * _S808;
    float _S810 = r2_5 * _S809;
    float _S811 = sx1_2 * dpfx_0 + _S807 + _S804 * _S808 + _S803 * _S809 + _S802 * _S810 + k4_2 * (r2_5 * _S810);
    float _S812 = v_5 * _S811;
    float _S813 = u_5 * _S811;
    float2  _S814 = make_float2 (1.0f + r2_5 * _S804) * _S805 + make_float2 (2.0f * u_5 * _S807 + 2.0f * (u_5 * _S807) + _S788 * (v_5 * dpfx_0) + _S813 + _S813, _S788 * u_5 * dpfx_0 + _S812 + _S812);
    float2  _S815 = _S783 * _S814;
    float2  _S816 = _S800 * _S814;
    float _S817 = _S815.x + _S815.y;
    if(_S793)
    {
        float _S818 = _S817 / _S795;
        float _S819 = _S796 * - _S818;
        float _S820 = _S792 * (0.3333333432674408f * - (_S785 * _S818));
        k_1 = _S820 + _S820;
        _S794 = _S819;
        _S795 = 0.0f;
    }
    else
    {
        float _S821 = _S817 / _S794;
        float _S822 = _S792 * - _S821;
        k_1 = _S784 * _S821;
        _S794 = 0.0f;
        _S795 = _S822;
    }
    DiffPair_float_0 _S823;
    (&_S823)->primal_0 = _S784;
    (&_S823)->differential_0 = 0.0f;
    DiffPair_float_0 _S824;
    (&_S824)->primal_0 = _S785;
    (&_S824)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S823, &_S824, k_1, &_s_diff_ctx_8->_S716);
    _s_diff_ctx_8->_S714 = _S823;
    _s_diff_ctx_8->_S715 = _S824;
    float _S825 = _S824.differential_0 + _S794;
    float _S826 = _S823.differential_0 + _S795;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S827;
    (&_S827)->primal_0 = _S783;
    (&_S827)->differential_0 = _S781;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S828 = { _S781, _S781 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S829;
    (&_S829)->_S742 = _S828;
    s_primal_ctx_s_bwd_length_impl_0(&_S827, _S826, &_S829);
    float2  _S830 = _S827.differential_0 + _S816;
    float3  _S831 = make_float3 (_S830.x, _S830.y, _S825);
    Matrix<float, 2, 3>  _S832 = J_12;
    _S832[int(0)] = _S831;
    if(_S793)
    {
        float _S833 = 1.0f - _S792 * _S792 / 3.0f;
        float _S834 = _S785 * _S785;
        k_1 = _S833 / _S785;
        _S794 = 0.0f;
        _S795 = _S834;
        _S796 = _S833;
    }
    else
    {
        float _S835 = _S784 * _S784;
        k_1 = _S792 / _S784;
        _S794 = _S835;
        _S795 = 0.0f;
        _S796 = 0.0f;
    }
    float2  _S836 = make_float2 (k_1);
    float2  _S837 = _S783 * make_float2 (k_1);
    float u_6 = _S837.x;
    float v_6 = _S837.y;
    float r2_6 = u_6 * u_6 + v_6 * v_6;
    float _S838 = k3_1 + r2_6 * k4_2;
    float _S839 = k2_1 + r2_6 * _S838;
    float _S840 = k1_1 + r2_6 * _S839;
    float2  _S841 = make_float2 (0.0f, dpfy_0);
    float2  _S842 = _S837 * _S841;
    float _S843 = p1_2 * dpfy_0;
    float _S844 = _S842.x + _S842.y;
    float _S845 = r2_6 * _S844;
    float _S846 = r2_6 * _S845;
    float _S847 = sy1_2 * dpfy_0 + _S843 + _S840 * _S844 + _S839 * _S845 + _S838 * _S846 + k4_2 * (r2_6 * _S846);
    float _S848 = v_6 * _S847;
    float _S849 = u_6 * _S847;
    float2  _S850 = make_float2 (1.0f + r2_6 * _S840) * _S841 + make_float2 (_S789 * (v_6 * dpfy_0) + _S849 + _S849, 2.0f * v_6 * _S843 + 2.0f * (v_6 * _S843) + _S789 * u_6 * dpfy_0 + _S848 + _S848);
    float2  _S851 = _S783 * _S850;
    float2  _S852 = _S836 * _S850;
    float _S853 = _S851.x + _S851.y;
    if(_S793)
    {
        float _S854 = _S853 / _S795;
        float _S855 = _S796 * - _S854;
        float _S856 = _S792 * (0.3333333432674408f * - (_S785 * _S854));
        k_1 = _S856 + _S856;
        _S794 = _S855;
        _S795 = 0.0f;
    }
    else
    {
        float _S857 = _S853 / _S794;
        float _S858 = _S792 * - _S857;
        k_1 = _S784 * _S857;
        _S794 = 0.0f;
        _S795 = _S858;
    }
    DiffPair_float_0 _S859;
    (&_S859)->primal_0 = _S784;
    (&_S859)->differential_0 = 0.0f;
    DiffPair_float_0 _S860;
    (&_S860)->primal_0 = _S785;
    (&_S860)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S859, &_S860, k_1, &_s_diff_ctx_8->_S719);
    _s_diff_ctx_8->_S717 = _S859;
    _s_diff_ctx_8->_S718 = _S860;
    float _S861 = _S860.differential_0 + _S794;
    float _S862 = _S859.differential_0 + _S795;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S863;
    (&_S863)->primal_0 = _S783;
    (&_S863)->differential_0 = _S781;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S864;
    (&_S864)->_S742 = _S828;
    s_primal_ctx_s_bwd_length_impl_0(&_S863, _S862, &_S864);
    float2  _S865 = _S863.differential_0 + _S852;
    float3  _S866 = make_float3 (_S865.x, _S865.y, _S861);
    _S832[int(1)] = _S866;
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S832, dpcov3d_0), transpose_1(_S832));
    *dpmean2d_0 = _S791;
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
    DiffPair_1 _S867 = *dpdpx_3;
    float _S868 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_9->_S747)->primal_0);
    float _S869 = s_primal_ctx_sqrt_0(_S868);
    float _S870 = 0.5f / _S869 * (*dpdpx_3).differential_0.differential_0;
    float _S871 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S869 * _S869));
    DiffPair_float_0 _S872;
    (&_S872)->primal_0 = _S868;
    (&_S872)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S872, _S871);
    DiffPair_float_0 _S873;
    (&_S873)->primal_0 = 1.00000001168609742e-07f;
    (&_S873)->differential_0 = 0.0f;
    DiffPair_float_0 _S874;
    (&_S874)->primal_0 = (&_s_diff_ctx_9->_S747)->primal_0;
    (&_S874)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S873, &_S874, _S872.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S874.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S870;
    dpdpx_3->primal_0 = _S867.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S875, DiffPair_float_0 * _S876, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_1 _S877 = *_S875;
    DiffPair_float_0 _S878 = _s_diff_ctx_10->_S743;
    DiffPair_float_0 _S879 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S880;
    (&_S880)->_S747 = _S879;
    s_primal_ctx_d_sqrt_0(&_S878, (*_S876).primal_0, &_S880);
    DiffPair_float_0 _S881 = { (*_S875).differential_0.primal_0, (*_S875).differential_0.differential_0 };
    DiffPair_1 _S882;
    (&_S882)->primal_0 = _s_diff_ctx_10->_S743;
    (&_S882)->differential_0 = _S881;
    DiffPair_float_0 _S883;
    (&_S883)->primal_0 = (*_S876).primal_0;
    (&_S883)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S884 = _S880;
    s_bwd_prop_d_sqrt_0(&_S882, &_S883, &_S884);
    DiffPair_float_0 _S885 = { _S882.differential_0.primal_0, _S882.differential_0.differential_0 };
    _S876->primal_0 = (*_S876).primal_0;
    _S876->differential_0 = _S883.differential_0;
    _S875->primal_0 = _S877.primal_0;
    _S875->differential_0 = _S885;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S886, float _s_dOut_3)
{
    DiffPair_float_0 _S887;
    (&_S887)->primal_0 = (*_S886).primal_0;
    (&_S887)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S887, _s_dOut_3);
    _S886->primal_0 = (*_S886).primal_0;
    _S886->differential_0 = _S887.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_0 _S888 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_11->_S744)->primal_0)->x) * *&((&(&_s_diff_ctx_11->_S744)->primal_0)->x) + *&((&(&_s_diff_ctx_11->_S744)->primal_0)->y) * *&((&(&_s_diff_ctx_11->_S744)->primal_0)->y);
    DiffPair_float_0 _S889 = { len_0, 0.0f };
    float2  _S890 = make_float2 (0.0f);
    float _S891 = (*dpdpx_5).differential_0.differential_0.x;
    float _S892 = _S891 + _S891;
    float _S893 = (&_s_diff_ctx_11->_S745)->differential_0 * _S892;
    float _S894 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S895 = (&_s_diff_ctx_11->_S745)->differential_0 * _S894;
    DiffPair_float_0 _S896 = { 0.0f, *&((&(&_s_diff_ctx_11->_S744)->primal_0)->x) * _S892 + *&((&(&_s_diff_ctx_11->_S744)->primal_0)->y) * _S894 };
    DiffPair_1 _S897;
    (&_S897)->primal_0 = _S889;
    (&_S897)->differential_0 = _S896;
    DiffPair_float_0 _S898;
    (&_S898)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S898)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S897, &_S898, &_s_diff_ctx_11->_S746);
    DiffPair_float_0 _S899;
    (&_S899)->primal_0 = len_0;
    (&_S899)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S899, 0.0f);
    float _S900 = _S897.differential_0.primal_0 + _S899.differential_0;
    float _S901 = *&((&(&_s_diff_ctx_11->_S744)->primal_0)->y) * _S900;
    float _S902 = _S895 + _S901 + _S901;
    float _S903 = *&((&(&_s_diff_ctx_11->_S744)->primal_0)->x) * _S900;
    float _S904 = _S893 + _S903 + _S903;
    float2  _S905 = _S890;
    *&((&_S905)->y) = _S902;
    *&((&_S905)->x) = _S904;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S888.differential_0.primal_0 + _S905, _S890 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S898.differential_0;
    dpdpx_5->primal_0 = _S888.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S906 = (*dpdpx_7).primal_0.x;
    float _S907 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S908;
    (&_S908)->primal_0 = _S906 * _S906 + _S907 * _S907;
    (&_S908)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S908, _s_dOut_4);
    float _S909 = (*dpdpx_7).primal_0.y * _S908.differential_0;
    float _S910 = _S909 + _S909;
    float _S911 = (*dpdpx_7).primal_0.x * _S908.differential_0;
    float _S912 = _S911 + _S911;
    float2  _S913 = make_float2 (0.0f);
    *&((&_S913)->y) = _S910;
    *&((&_S913)->x) = _S912;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S913;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S914, DiffPair_float_0 * _S915, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_0 _S916 = *_S914;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S917 = _s_diff_ctx_12->_S742;
    float2  _S918 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S919 = { _S918, _S918 };
    DiffPair_float_0 _S920 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S921 = { _S920 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S922;
    (&_S922)->_S744 = _S919;
    (&_S922)->_S745 = _S920;
    (&_S922)->_S746 = _S921;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S917, (*_S915).primal_0, &_S922);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S923 = { (*_S914).differential_0.primal_0, (*_S914).differential_0.differential_0 };
    DiffPair_0 _S924;
    (&_S924)->primal_0 = _s_diff_ctx_12->_S742;
    (&_S924)->differential_0 = _S923;
    DiffPair_float_0 _S925;
    (&_S925)->primal_0 = (*_S915).primal_0;
    (&_S925)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S926 = _S922;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S924, &_S925, &_S926);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S927;
    (&_S927)->primal_0 = (&_s_diff_ctx_12->_S742)->primal_0;
    (&_S927)->differential_0 = _S918;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S927, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S928 = { _S924.differential_0.primal_0 + _S927.differential_0, _S924.differential_0.differential_0 };
    _S915->primal_0 = (*_S915).primal_0;
    _S915->differential_0 = _S925.differential_0;
    _S914->primal_0 = _S916.primal_0;
    _S914->differential_0 = _S928;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_1 _S929 = *dpdpy_1;
    DiffPair_1 _S930 = *dpdpx_8;
    float _S931 = - (&_s_diff_ctx_13->_S725)->primal_0;
    float _S932 = (&_s_diff_ctx_13->_S726)->primal_0 * (&_s_diff_ctx_13->_S726)->primal_0 + (&_s_diff_ctx_13->_S725)->primal_0 * (&_s_diff_ctx_13->_S725)->primal_0;
    float _S933 = _S932 * _S932;
    float _S934 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S933;
    float _S935 = (&_s_diff_ctx_13->_S726)->primal_0 * - _S934;
    float _S936 = (&_s_diff_ctx_13->_S725)->primal_0 * _S935;
    float _S937 = (&_s_diff_ctx_13->_S726)->primal_0 * _S935;
    float _S938 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S933;
    float _S939 = _S931 * - _S938;
    float _S940 = (&_s_diff_ctx_13->_S725)->primal_0 * _S939;
    float _S941 = (&_s_diff_ctx_13->_S726)->primal_0 * _S939;
    DiffPair_float_0 dpdpx_9 = { _S941 + _S941 + ((*dpdpx_8).differential_0.primal_0 + (_S937 + _S937 + _S932 * _S934)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S936 + _S936 + (*dpdpy_1).differential_0.primal_0 + _S940 + _S940 + - (_S932 * _S938), 0.0f };
    float _S942 = (&_s_diff_ctx_13->_S726)->primal_0 / _S932 * (*dpdpy_1).differential_0.differential_0 + _S931 / _S932 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S942;
    dpdpy_1->primal_0 = _S929.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S930.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S943, DiffPair_1 * _S944, DiffPair_float_0 * _S945, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_1 _S946 = *_S943;
    DiffPair_1 _S947 = *_S944;
    DiffPair_float_0 _S948 = _s_diff_ctx_14->_S712;
    DiffPair_float_0 _S949 = _s_diff_ctx_14->_S713;
    DiffPair_float_0 _S950 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S951;
    (&_S951)->_S725 = _S950;
    (&_S951)->_S726 = _S950;
    s_primal_ctx_d_atan2_0(&_S948, &_S949, (*_S945).primal_0, &_S951);
    DiffPair_float_0 _S952 = { (*_S944).differential_0.primal_0, (*_S944).differential_0.differential_0 };
    DiffPair_float_0 _S953 = { (*_S943).differential_0.primal_0, (*_S943).differential_0.differential_0 };
    DiffPair_1 _S954;
    (&_S954)->primal_0 = _s_diff_ctx_14->_S712;
    (&_S954)->differential_0 = _S953;
    DiffPair_1 _S955;
    (&_S955)->primal_0 = _s_diff_ctx_14->_S713;
    (&_S955)->differential_0 = _S952;
    DiffPair_float_0 _S956;
    (&_S956)->primal_0 = (*_S945).primal_0;
    (&_S956)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S957 = _S951;
    s_bwd_prop_d_atan2_0(&_S954, &_S955, &_S956, &_S957);
    DiffPair_float_0 _S958 = { _S955.differential_0.primal_0, _S955.differential_0.differential_0 };
    DiffPair_float_0 _S959 = { _S954.differential_0.primal_0, _S954.differential_0.differential_0 };
    _S945->primal_0 = (*_S945).primal_0;
    _S945->differential_0 = _S956.differential_0;
    _S943->primal_0 = _S946.primal_0;
    _S943->differential_0 = _S959;
    _S944->primal_0 = _S947.primal_0;
    _S944->differential_0 = _S958;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S960, DiffPair_float_0 * _S961, float _s_dOut_5)
{
    DiffPair_float_0 _S962;
    (&_S962)->primal_0 = (*_S960).primal_0;
    (&_S962)->differential_0 = 0.0f;
    DiffPair_float_0 _S963;
    (&_S963)->primal_0 = (*_S961).primal_0;
    (&_S963)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S962, &_S963, _s_dOut_5);
    _S961->primal_0 = (*_S961).primal_0;
    _S961->differential_0 = _S963.differential_0;
    _S960->primal_0 = (*_S960).primal_0;
    _S960->differential_0 = _S962.differential_0;
    return;
}

inline __device__ void s_bwd_prop_CameraDistortion_x24init_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_2, CameraDistortion_0 * _s_dOut_6)
{
    float2  _S964 = _s_dOut_6->thin_prism_coeffs_0;
    dpthin_prism_coeffs_2->primal_0 = (*dpthin_prism_coeffs_2).primal_0;
    dpthin_prism_coeffs_2->differential_0 = _S964;
    float2  _S965 = _s_dOut_6->tangential_coeffs_0;
    dptangential_coeffs_2->primal_0 = (*dptangential_coeffs_2).primal_0;
    dptangential_coeffs_2->differential_0 = _S965;
    float4  _S966 = _s_dOut_6->radial_coeffs_0;
    dpradial_coeffs_2->primal_0 = (*dpradial_coeffs_2).primal_0;
    dpradial_coeffs_2->differential_0 = _S966;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_float_0 * dpfx_1, DiffPair_float_0 * dpfy_1, DiffPair_float_0 * dpcx_1, DiffPair_float_0 * dpcy_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpradial_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dptangential_coeffs_3, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpthin_prism_coeffs_3, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S967 = *dpcov3d_1;
    DiffPair_float_0 _S968 = *dpfx_1;
    DiffPair_float_0 _S969 = *dpfy_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S970 = *dpradial_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S971 = *dptangential_coeffs_3;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S972 = *dpthin_prism_coeffs_3;
    float2  _S973 = make_float2 (0.0f);
    CameraDistortion_0 _S974 = s_primal_ctx_CameraDistortion_x24init_0((*dpradial_coeffs_3).primal_0, (*dptangential_coeffs_3).primal_0, (*dpthin_prism_coeffs_3).primal_0);
    float2  _S975 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S976 = length_0(_S975);
    float _S977 = (*dpmean3d_1).primal_0.z;
    float _S978 = s_primal_ctx_atan2_0(_S976, _S977);
    bool _S979 = _S978 < 0.00100000004749745f;
    float k_2;
    float _S980;
    float _S981;
    float _S982;
    if(_S979)
    {
        float _S983 = 1.0f - _S978 * _S978 / 3.0f;
        float _S984 = _S977 * _S977;
        k_2 = _S983 / _S977;
        _S980 = 0.0f;
        _S981 = _S984;
        _S982 = _S983;
    }
    else
    {
        float _S985 = _S976 * _S976;
        k_2 = _S978 / _S976;
        _S980 = _S985;
        _S981 = 0.0f;
        _S982 = 0.0f;
    }
    float2  _S986 = make_float2 (k_2);
    float2  _S987 = _S975 * make_float2 (k_2);
    float k1_2 = _S974.radial_coeffs_0.x;
    float k2_2 = _S974.radial_coeffs_0.y;
    float k3_2 = _S974.radial_coeffs_0.z;
    float k4_3 = _S974.radial_coeffs_0.w;
    float p1_3 = _S974.tangential_coeffs_0.x;
    float p2_3 = _S974.tangential_coeffs_0.y;
    float sx1_3 = _S974.thin_prism_coeffs_0.x;
    float sy1_3 = _S974.thin_prism_coeffs_0.y;
    float u_7 = _S987.x;
    float v_7 = _S987.y;
    float r2_7 = u_7 * u_7 + v_7 * v_7;
    float _S988 = k3_2 + r2_7 * k4_3;
    float _S989 = k2_2 + r2_7 * _S988;
    float _S990 = k1_2 + r2_7 * _S989;
    float radial_1 = 1.0f + r2_7 * _S990;
    float2  _S991 = make_float2 (radial_1);
    float _S992 = 2.0f * p1_3;
    float _S993 = _S992 * u_7;
    float _S994 = 2.0f * u_7;
    float _S995 = r2_7 + _S994 * u_7;
    float _S996 = 2.0f * p2_3;
    float _S997 = _S996 * u_7;
    float _S998 = 2.0f * v_7;
    float _S999 = r2_7 + _S998 * v_7;
    float2  _S1000 = _S987 * make_float2 (radial_1) + make_float2 (_S993 * v_7 + p2_3 * _S995 + sx1_3 * r2_7, _S997 * v_7 + p1_3 * _S999 + sy1_3 * r2_7);
    float _S1001 = _S1000.x;
    float _S1002 = _S1000.y;
    Matrix<float, 2, 3>  J_13 = makeMatrix<float, 2, 3> (0.0f);
    float _S1003 = s_primal_ctx_s_primal_ctx_atan2_0(_S976, _S977);
    bool _S1004 = _S1003 < 0.00100000004749745f;
    float _S1005;
    float _S1006;
    float _S1007;
    if(_S1004)
    {
        float _S1008 = 1.0f - _S1003 * _S1003 / 3.0f;
        float _S1009 = _S977 * _S977;
        k_2 = _S1008 / _S977;
        _S1005 = 0.0f;
        _S1006 = _S1009;
        _S1007 = _S1008;
    }
    else
    {
        float _S1010 = _S976 * _S976;
        k_2 = _S1003 / _S976;
        _S1005 = _S1010;
        _S1006 = 0.0f;
        _S1007 = 0.0f;
    }
    float2  _S1011 = make_float2 (k_2);
    float2  _S1012 = _S975 * make_float2 (k_2);
    float u_8 = _S1012.x;
    float v_8 = _S1012.y;
    float r2_8 = u_8 * u_8 + v_8 * v_8;
    float _S1013 = k3_2 + r2_8 * k4_3;
    float _S1014 = k2_2 + r2_8 * _S1013;
    float _S1015 = k1_2 + r2_8 * _S1014;
    float2  _S1016 = make_float2 (1.0f + r2_8 * _S1015);
    float _S1017 = _S992 * u_8;
    float _S1018 = 2.0f * u_8;
    float2  _S1019 = make_float2 (_S968.primal_0, 0.0f);
    float2  _S1020 = _S1012 * _S1019;
    float _S1021 = p2_3 * _S968.primal_0;
    float _S1022 = v_8 * _S968.primal_0;
    float _S1023 = _S1020.x + _S1020.y;
    float _S1024 = r2_8 * _S1023;
    float _S1025 = r2_8 * _S1024;
    float _S1026 = r2_8 * _S1025;
    float _S1027 = sx1_3 * _S968.primal_0 + _S1021 + _S1015 * _S1023 + _S1014 * _S1024 + _S1013 * _S1025 + k4_3 * _S1026;
    float _S1028 = v_8 * _S1027;
    float _S1029 = u_8 * _S1027;
    float2  _S1030 = _S1016 * _S1019 + make_float2 (_S1018 * _S1021 + 2.0f * (u_8 * _S1021) + _S992 * _S1022 + _S1029 + _S1029, _S1017 * _S968.primal_0 + _S1028 + _S1028);
    float2  _S1031 = _S975 * _S1030;
    float2  _S1032 = _S1011 * _S1030;
    float _S1033 = _S1031.x + _S1031.y;
    float k_3;
    float _S1034;
    float _S1035;
    float _S1036;
    float _S1037;
    float _S1038;
    float _S1039;
    float _S1040;
    float _S1041;
    if(_S1004)
    {
        float _S1042 = _S1033 / _S1006;
        float _S1043 = _S1006 * _S1006;
        float _S1044 = - _S1042;
        float _S1045 = _S1007 * _S1044;
        float _S1046 = 0.3333333432674408f * - (_S977 * _S1042);
        float _S1047 = _S1003 * _S1046;
        k_2 = _S1047 + _S1047;
        k_3 = _S1045;
        _S1034 = 0.0f;
        _S1035 = 0.0f;
        _S1036 = 0.0f;
        _S1037 = 0.0f;
        _S1038 = _S1046;
        _S1039 = _S1042;
        _S1040 = _S1044;
        _S1041 = _S1043;
    }
    else
    {
        float _S1048 = _S1033 / _S1005;
        float _S1049 = _S1005 * _S1005;
        float _S1050 = - _S1048;
        float _S1051 = _S1003 * _S1050;
        k_2 = _S976 * _S1048;
        k_3 = 0.0f;
        _S1034 = _S1051;
        _S1035 = _S1048;
        _S1036 = _S1050;
        _S1037 = _S1049;
        _S1038 = 0.0f;
        _S1039 = 0.0f;
        _S1040 = 0.0f;
        _S1041 = 0.0f;
    }
    DiffPair_float_0 _S1052 = { _S976, 0.0f };
    DiffPair_float_0 _S1053 = { _S977, 0.0f };
    float _S1054 = (&_s_diff_ctx_15->_S715)->differential_0 + k_3;
    float _S1055 = (&_s_diff_ctx_15->_S714)->differential_0 + _S1034;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1056 = { _S975, _S973 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1057;
    (&_S1057)->primal_0 = _S975;
    (&_S1057)->differential_0 = _S973;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1058 = { _S973, _S973 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1059;
    (&_S1059)->_S742 = _S1058;
    s_primal_ctx_s_bwd_length_impl_0(&_S1057, _S1055, &_S1059);
    float2  _S1060 = _S1057.differential_0 + _S1032;
    float3  _S1061 = make_float3 (_S1060.x, _S1060.y, _S1054);
    Matrix<float, 2, 3>  _S1062 = J_13;
    _S1062[int(0)] = _S1061;
    float _S1063;
    float _S1064;
    if(_S1004)
    {
        float _S1065 = 1.0f - _S1003 * _S1003 / 3.0f;
        float _S1066 = _S977 * _S977;
        k_3 = _S1065 / _S977;
        _S1034 = 0.0f;
        _S1063 = _S1066;
        _S1064 = _S1065;
    }
    else
    {
        float _S1067 = _S976 * _S976;
        k_3 = _S1003 / _S976;
        _S1034 = _S1067;
        _S1063 = 0.0f;
        _S1064 = 0.0f;
    }
    float2  _S1068 = make_float2 (k_3);
    float2  _S1069 = _S975 * make_float2 (k_3);
    float u_9 = _S1069.x;
    float v_9 = _S1069.y;
    float r2_9 = u_9 * u_9 + v_9 * v_9;
    float _S1070 = k3_2 + r2_9 * k4_3;
    float _S1071 = k2_2 + r2_9 * _S1070;
    float _S1072 = k1_2 + r2_9 * _S1071;
    float2  _S1073 = make_float2 (1.0f + r2_9 * _S1072);
    float _S1074 = _S996 * u_9;
    float _S1075 = 2.0f * v_9;
    float2  _S1076 = make_float2 (0.0f, _S969.primal_0);
    float2  _S1077 = _S1069 * _S1076;
    float _S1078 = p1_3 * _S969.primal_0;
    float _S1079 = v_9 * _S969.primal_0;
    float _S1080 = _S1077.x + _S1077.y;
    float _S1081 = r2_9 * _S1080;
    float _S1082 = r2_9 * _S1081;
    float _S1083 = r2_9 * _S1082;
    float _S1084 = sy1_3 * _S969.primal_0 + _S1078 + _S1072 * _S1080 + _S1071 * _S1081 + _S1070 * _S1082 + k4_3 * _S1083;
    float _S1085 = v_9 * _S1084;
    float _S1086 = u_9 * _S1084;
    float2  _S1087 = _S1073 * _S1076 + make_float2 (_S996 * _S1079 + _S1086 + _S1086, _S1075 * _S1078 + 2.0f * (v_9 * _S1078) + _S1074 * _S969.primal_0 + _S1085 + _S1085);
    float2  _S1088 = _S975 * _S1087;
    float2  _S1089 = _S1068 * _S1087;
    float _S1090 = _S1088.x + _S1088.y;
    float _S1091;
    float _S1092;
    float _S1093;
    float _S1094;
    float _S1095;
    float _S1096;
    float _S1097;
    float _S1098;
    float _S1099;
    if(_S1004)
    {
        float _S1100 = _S1090 / _S1063;
        float _S1101 = _S1063 * _S1063;
        float _S1102 = - _S1100;
        float _S1103 = _S1064 * _S1102;
        float _S1104 = 0.3333333432674408f * - (_S977 * _S1100);
        float _S1105 = _S1003 * _S1104;
        k_3 = _S1105 + _S1105;
        _S1091 = _S1103;
        _S1092 = 0.0f;
        _S1093 = 0.0f;
        _S1094 = 0.0f;
        _S1095 = 0.0f;
        _S1096 = _S1104;
        _S1097 = _S1100;
        _S1098 = _S1102;
        _S1099 = _S1101;
    }
    else
    {
        float _S1106 = _S1090 / _S1034;
        float _S1107 = _S1034 * _S1034;
        float _S1108 = - _S1106;
        float _S1109 = _S1003 * _S1108;
        k_3 = _S976 * _S1106;
        _S1091 = 0.0f;
        _S1092 = _S1109;
        _S1093 = _S1106;
        _S1094 = _S1108;
        _S1095 = _S1107;
        _S1096 = 0.0f;
        _S1097 = 0.0f;
        _S1098 = 0.0f;
        _S1099 = 0.0f;
    }
    float _S1110 = (&_s_diff_ctx_15->_S718)->differential_0 + _S1091;
    float _S1111 = (&_s_diff_ctx_15->_S717)->differential_0 + _S1092;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1112;
    (&_S1112)->primal_0 = _S975;
    (&_S1112)->differential_0 = _S973;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1113;
    (&_S1113)->_S742 = _S1058;
    s_primal_ctx_s_bwd_length_impl_0(&_S1112, _S1111, &_S1113);
    float2  _S1114 = _S1112.differential_0 + _S1089;
    float3  _S1115 = make_float3 (_S1114.x, _S1114.y, _S1110);
    _S1062[int(1)] = _S1115;
    Matrix<float, 3, 2>  _S1116 = transpose_1(_S1062);
    CameraDistortion_0 _S1117 = CameraDistortion_x24_syn_dzero_0();
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1118;
    (&_S1118)->primal_0 = s_primal_ctx_mul_3(_S1062, _S967.primal_0);
    (&_S1118)->differential_0 = J_13;
    Matrix<float, 3, 2>  _S1119 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1120;
    (&_S1120)->primal_0 = _S1116;
    (&_S1120)->differential_0 = _S1119;
    s_bwd_prop_mul_2(&_S1118, &_S1120, dpcov2d_1);
    Matrix<float, 2, 3>  _S1121 = transpose_2(_S1120.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1122;
    (&_S1122)->primal_0 = _S1062;
    (&_S1122)->differential_0 = J_13;
    Matrix<float, 3, 3>  _S1123 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1124;
    (&_S1124)->primal_0 = _S967.primal_0;
    (&_S1124)->differential_0 = _S1123;
    s_bwd_prop_mul_3(&_S1122, &_S1124, _S1118.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1125 = _S1124;
    Matrix<float, 2, 3>  _S1126 = _S1121 + _S1122.differential_0;
    float2  _S1127 = _S973;
    *&((&_S1127)->y) = _S1126.rows[int(1)].y;
    *&((&_S1127)->x) = _S1126.rows[int(1)].x;
    float2  _S1128 = _S1127;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1129 = { _S973, _S1127 };
    DiffPair_0 _S1130;
    (&_S1130)->primal_0 = _S1056;
    (&_S1130)->differential_0 = _S1129;
    DiffPair_float_0 _S1131;
    (&_S1131)->primal_0 = _S1111;
    (&_S1131)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1132 = _S1113;
    s_bwd_prop_s_bwd_length_impl_0(&_S1130, &_S1131, &_S1132);
    DiffPair_0 _S1133 = _S1130;
    DiffPair_float_0 _S1134 = _S1131;
    DiffPair_float_0 _S1135 = { 0.0f, _S1126.rows[int(1)].z };
    DiffPair_float_0 _S1136 = { 0.0f, _S1131.differential_0 };
    DiffPair_1 _S1137;
    (&_S1137)->primal_0 = _S1052;
    (&_S1137)->differential_0 = _S1136;
    DiffPair_1 _S1138;
    (&_S1138)->primal_0 = _S1053;
    (&_S1138)->differential_0 = _S1135;
    DiffPair_float_0 _S1139;
    (&_S1139)->primal_0 = k_3;
    (&_S1139)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1137, &_S1138, &_S1139, &_s_diff_ctx_15->_S719);
    DiffPair_1 _S1140 = _S1137;
    DiffPair_1 _S1141 = _S1138;
    DiffPair_float_0 _S1142 = _S1139;
    if(_S1004)
    {
        float _S1143 = _S1142.differential_0 + _S1142.differential_0;
        float _S1144 = _S1096 * _S1143;
        float _S1145 = - (0.3333333432674408f * (_S1003 * _S1143));
        float _S1146 = _S1098 * _S1126.rows[int(1)].z;
        float _S1147 = (_S977 * _S1145 + - (_S1064 * _S1126.rows[int(1)].z)) / _S1099;
        float _S1148 = _S1090 * - _S1147;
        float _S1149 = _S1097 * _S1145 + _S1141.differential_0.primal_0;
        k_3 = _S1063 * _S1147;
        _S1091 = 0.0f;
        _S1092 = _S1148;
        _S1093 = _S1146;
        _S1094 = _S1140.differential_0.primal_0;
        _S1095 = _S1144;
        _S1096 = _S1149;
    }
    else
    {
        float _S1150 = _S1094 * _S1134.differential_0;
        float _S1151 = (_S976 * _S1142.differential_0 + - (_S1003 * _S1134.differential_0)) / _S1095;
        float _S1152 = _S1090 * - _S1151;
        float _S1153 = _S1093 * _S1142.differential_0 + _S1140.differential_0.primal_0;
        k_3 = _S1034 * _S1151;
        _S1091 = _S1152;
        _S1092 = 0.0f;
        _S1093 = 0.0f;
        _S1094 = _S1153;
        _S1095 = _S1150;
        _S1096 = _S1141.differential_0.primal_0;
    }
    float2  _S1154 = _S1068 * _S1128;
    float2  _S1155 = _S1087 * _S1128;
    float2  _S1156 = _S973;
    *&((&_S1156)->y) = k_3;
    *&((&_S1156)->x) = k_3;
    float2  _S1157 = _S1087 * _S1156;
    float2  _S1158 = _S1154 + _S975 * _S1156;
    float _S1159 = _S1158.x;
    float _S1160 = _S1159 + _S1159;
    float _S1161 = _S1084 * _S1160;
    float _S1162 = _S1158.y + _S1158.y;
    float _S1163 = _S1084 * _S1162;
    float _S1164 = u_9 * _S1160 + v_9 * _S1162;
    float _S1165 = k4_3 * _S1164;
    float _S1166 = _S1083 * _S1164;
    float _S1167 = _S1082 * _S1164;
    float _S1168 = _S1082 * _S1165;
    float _S1169 = _S1081 * _S1164;
    float _S1170 = _S1070 * _S1164 + r2_9 * _S1165;
    float _S1171 = _S1081 * _S1170;
    float _S1172 = _S1080 * _S1164;
    float _S1173 = _S1071 * _S1164 + r2_9 * _S1170;
    float _S1174 = _S1080 * _S1173;
    float _S1175 = _S1072 * _S1164 + r2_9 * _S1173;
    float _S1176 = _S996 * _S1158.x;
    float _S1177 = _S1079 * _S1158.x;
    float _S1178 = v_9 * _S1176;
    float _S1179 = _S969.primal_0 * _S1176;
    float _S1180 = _S1074 * _S1158.y;
    float _S1181 = _S969.primal_0 * _S1158.y;
    float _S1182 = 2.0f * _S1158.y;
    float _S1183 = _S1078 * _S1182;
    float _S1184 = _S1078 * _S1158.y;
    float _S1185 = _S1164 + v_9 * _S1182 + _S1075 * _S1158.y;
    float _S1186 = p1_3 * _S1185;
    float _S1187 = _S969.primal_0 * _S1185;
    float _S1188 = sy1_3 * _S1164;
    float _S1189 = _S969.primal_0 * _S1164;
    float2  _S1190 = _S1073 * _S1158;
    float2  _S1191 = _S1076 * _S1158;
    float2  _S1192 = _S973;
    *&((&_S1192)->y) = _S1175;
    *&((&_S1192)->x) = _S1175;
    float _S1193 = _S1191.x + _S1191.y;
    float _S1194 = _S1172 + r2_9 * _S1193;
    float _S1195 = _S1169 + r2_9 * _S1194;
    float _S1196 = _S1167 + r2_9 * _S1195;
    float _S1197 = _S1168 + _S1171 + _S1174 + _S1072 * _S1193 + _S1071 * _S1194 + _S1070 * _S1195 + k4_3 * _S1196;
    float _S1198 = v_9 * _S1197;
    float _S1199 = u_9 * _S1197;
    float2  _S1200 = _S1157 + _S1133.differential_0.primal_0;
    float2  _S1201 = _S1076 * _S1192 + make_float2 (_S1161 + _S996 * _S1181 + _S1199 + _S1199, _S1163 + _S1179 + _S1183 + 2.0f * _S1184 + _S1198 + _S1198);
    float _S1202 = _S1178 + _S1180 + _S1186 + _S1188 + (_S1190 + _S1069 * _S1192).y;
    float _S1203 = _S1177 + u_9 * _S1181;
    float _S1204 = _S1166 + r2_9 * _S1196;
    float2  _S1205 = _S975 * _S1201;
    float _S1206 = _S1155.x + _S1155.y + _S1205.x + _S1205.y;
    float2  _S1207 = _S1068 * _S1201 + _S1200;
    if(_S1004)
    {
        float _S1208 = _S977 * _S1092;
        float _S1209 = _S1206 / _S1063;
        float _S1210 = _S1003 * (0.3333333432674408f * - (_S1093 + _S977 * _S1209));
        float _S1211 = _S1208 + _S1208 + _S1064 * - _S1209 + _S1096;
        k_3 = _S1210 + _S1210 + _S1095;
        _S1034 = _S1211;
        _S1063 = _S1094;
    }
    else
    {
        float _S1212 = _S976 * _S1091;
        float _S1213 = _S1206 / _S1034;
        float _S1214 = _S1212 + _S1212 + _S1003 * - _S1213 + _S1094;
        k_3 = _S976 * _S1213 + _S1095;
        _S1034 = _S1096;
        _S1063 = _S1214;
    }
    DiffPair_float_0 _S1215;
    (&_S1215)->primal_0 = _S976;
    (&_S1215)->differential_0 = 0.0f;
    DiffPair_float_0 _S1216;
    (&_S1216)->primal_0 = _S977;
    (&_S1216)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1215, &_S1216, k_3);
    float _S1217 = _S1216.differential_0 + _S1034;
    float _S1218 = _S1215.differential_0 + _S1063;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1219;
    (&_S1219)->primal_0 = _S975;
    (&_S1219)->differential_0 = _S973;
    s_bwd_length_impl_1(&_S1219, _S1218);
    float2  _S1220 = _S1219.differential_0 + _S1207;
    float3  _S1221 = make_float3 (_S1220.x, _S1220.y, _S1217);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1222;
    (&_S1222)->primal_0 = _S975;
    (&_S1222)->differential_0 = _S973;
    s_bwd_length_impl_1(&_S1222, 0.0f);
    float3  _S1223 = _S1221 + make_float3 (_S1222.differential_0.x, _S1222.differential_0.y, 0.0f);
    float2  _S1224 = _S973;
    *&((&_S1224)->y) = _S1126.rows[int(0)].y;
    *&((&_S1224)->x) = _S1126.rows[int(0)].x;
    float2  _S1225 = _S1224;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1226 = { _S973, _S1224 };
    DiffPair_0 _S1227;
    (&_S1227)->primal_0 = _S1056;
    (&_S1227)->differential_0 = _S1226;
    DiffPair_float_0 _S1228;
    (&_S1228)->primal_0 = _S1055;
    (&_S1228)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1229 = _S1059;
    s_bwd_prop_s_bwd_length_impl_0(&_S1227, &_S1228, &_S1229);
    DiffPair_0 _S1230 = _S1227;
    DiffPair_float_0 _S1231 = _S1228;
    DiffPair_float_0 _S1232 = { 0.0f, _S1126.rows[int(0)].z };
    DiffPair_float_0 _S1233 = { 0.0f, _S1228.differential_0 };
    DiffPair_1 _S1234;
    (&_S1234)->primal_0 = _S1052;
    (&_S1234)->differential_0 = _S1233;
    DiffPair_1 _S1235;
    (&_S1235)->primal_0 = _S1053;
    (&_S1235)->differential_0 = _S1232;
    DiffPair_float_0 _S1236;
    (&_S1236)->primal_0 = k_2;
    (&_S1236)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1234, &_S1235, &_S1236, &_s_diff_ctx_15->_S716);
    DiffPair_1 _S1237 = _S1234;
    DiffPair_1 _S1238 = _S1235;
    DiffPair_float_0 _S1239 = _S1236;
    if(_S1004)
    {
        float _S1240 = _S1239.differential_0 + _S1239.differential_0;
        float _S1241 = _S1038 * _S1240;
        float _S1242 = - (0.3333333432674408f * (_S1003 * _S1240));
        float _S1243 = _S1040 * _S1126.rows[int(0)].z;
        float _S1244 = (_S977 * _S1242 + - (_S1007 * _S1126.rows[int(0)].z)) / _S1041;
        float _S1245 = _S1033 * - _S1244;
        float _S1246 = _S1039 * _S1242 + _S1238.differential_0.primal_0;
        k_2 = _S1006 * _S1244;
        k_3 = 0.0f;
        _S1034 = _S1245;
        _S1035 = _S1243;
        _S1036 = _S1237.differential_0.primal_0;
        _S1037 = _S1241;
        _S1038 = _S1246;
    }
    else
    {
        float _S1247 = _S1036 * _S1231.differential_0;
        float _S1248 = (_S976 * _S1239.differential_0 + - (_S1003 * _S1231.differential_0)) / _S1037;
        float _S1249 = _S1033 * - _S1248;
        float _S1250 = _S1035 * _S1239.differential_0 + _S1237.differential_0.primal_0;
        k_2 = _S1005 * _S1248;
        k_3 = _S1249;
        _S1034 = 0.0f;
        _S1035 = 0.0f;
        _S1036 = _S1250;
        _S1037 = _S1247;
        _S1038 = _S1238.differential_0.primal_0;
    }
    float2  _S1251 = _S1011 * _S1225;
    float2  _S1252 = _S1030 * _S1225;
    float2  _S1253 = _S973;
    *&((&_S1253)->y) = k_2;
    *&((&_S1253)->x) = k_2;
    float2  _S1254 = _S1030 * _S1253;
    float2  _S1255 = _S1251 + _S975 * _S1253;
    float _S1256 = _S1255.x;
    float _S1257 = _S1256 + _S1256;
    float _S1258 = _S1027 * _S1257;
    float _S1259 = _S1255.y + _S1255.y;
    float _S1260 = _S1027 * _S1259;
    float _S1261 = u_8 * _S1257 + v_8 * _S1259;
    float _S1262 = k4_3 * _S1261;
    float _S1263 = _S1026 * _S1261;
    float _S1264 = _S1025 * _S1261;
    float _S1265 = _S1025 * _S1262;
    float _S1266 = _S1024 * _S1261;
    float _S1267 = _S1013 * _S1261 + r2_8 * _S1262;
    float _S1268 = _S1024 * _S1267;
    float _S1269 = _S1023 * _S1261;
    float _S1270 = _S1014 * _S1261 + r2_8 * _S1267;
    float _S1271 = _S1023 * _S1270;
    float _S1272 = _S1015 * _S1261 + r2_8 * _S1270;
    float _S1273 = _S992 * _S1255.x;
    float _S1274 = _S1022 * _S1255.x;
    float _S1275 = v_8 * _S1273;
    float _S1276 = _S968.primal_0 * _S1273;
    float _S1277 = _S1017 * _S1255.y;
    float _S1278 = _S968.primal_0 * _S1255.y;
    float _S1279 = 2.0f * _S1255.x;
    float _S1280 = _S1021 * _S1279;
    float _S1281 = _S1021 * _S1255.x;
    float _S1282 = _S1261 + u_8 * _S1279 + _S1018 * _S1255.x;
    float _S1283 = p2_3 * _S1282;
    float _S1284 = _S968.primal_0 * _S1282;
    float _S1285 = sx1_3 * _S1261;
    float _S1286 = _S968.primal_0 * _S1261;
    float2  _S1287 = _S1016 * _S1255;
    float2  _S1288 = _S1019 * _S1255;
    float2  _S1289 = _S973;
    *&((&_S1289)->y) = _S1272;
    *&((&_S1289)->x) = _S1272;
    float _S1290 = _S1288.x + _S1288.y;
    float _S1291 = _S1269 + r2_8 * _S1290;
    float _S1292 = _S1266 + r2_8 * _S1291;
    float _S1293 = _S1264 + r2_8 * _S1292;
    float _S1294 = _S1265 + _S1268 + _S1271 + _S1015 * _S1290 + _S1014 * _S1291 + _S1013 * _S1292 + k4_3 * _S1293;
    float _S1295 = v_8 * _S1294;
    float _S1296 = u_8 * _S1294;
    float2  _S1297 = _S1254 + _S1230.differential_0.primal_0;
    float _S1298 = _S1274 + u_8 * _S1278;
    float _S1299 = _S1291 + _S1194;
    float _S1300 = _S1293 + _S1196;
    float _S1301 = _S1275 + _S1277 + _S1283 + _S1285 + (_S1287 + _S1012 * _S1289).x;
    float2  _S1302 = _S1019 * _S1289 + make_float2 (_S1258 + _S1280 + 2.0f * _S1281 + _S992 * _S1278 + _S1296 + _S1296, _S1260 + _S1276 + _S1295 + _S1295);
    float _S1303 = _S1292 + _S1195;
    float _S1304 = _S1263 + r2_8 * _S1293 + _S1204;
    float2  _S1305 = _S975 * _S1302;
    float _S1306 = _S1252.x + _S1252.y + _S1305.x + _S1305.y;
    float2  _S1307 = _S1011 * _S1302 + _S1297;
    if(_S1004)
    {
        float _S1308 = _S977 * _S1034;
        float _S1309 = _S1306 / _S1006;
        float _S1310 = _S1003 * (0.3333333432674408f * - (_S1035 + _S977 * _S1309));
        float _S1311 = _S1308 + _S1308 + _S1007 * - _S1309 + _S1038;
        k_2 = _S1310 + _S1310 + _S1037;
        _S1005 = _S1311;
        _S1006 = _S1036;
    }
    else
    {
        float _S1312 = _S976 * k_3;
        float _S1313 = _S1306 / _S1005;
        float _S1314 = _S1312 + _S1312 + _S1003 * - _S1313 + _S1036;
        k_2 = _S976 * _S1313 + _S1037;
        _S1005 = _S1038;
        _S1006 = _S1314;
    }
    DiffPair_float_0 _S1315;
    (&_S1315)->primal_0 = _S976;
    (&_S1315)->differential_0 = 0.0f;
    DiffPair_float_0 _S1316;
    (&_S1316)->primal_0 = _S977;
    (&_S1316)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1315, &_S1316, k_2);
    float _S1317 = _S1316.differential_0 + _S1005;
    float _S1318 = _S1315.differential_0 + _S1006;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1319;
    (&_S1319)->primal_0 = _S975;
    (&_S1319)->differential_0 = _S973;
    s_bwd_length_impl_1(&_S1319, _S1318);
    float2  _S1320 = _S1319.differential_0 + _S1307;
    float3  _S1321 = make_float3 (_S1320.x, _S1320.y, _S1317);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1322;
    (&_S1322)->primal_0 = _S975;
    (&_S1322)->differential_0 = _S973;
    s_bwd_length_impl_1(&_S1322, 0.0f);
    float _S1323 = _S969.primal_0 * dpmean2d_1.y;
    float _S1324 = _S968.primal_0 * dpmean2d_1.x;
    float2  _S1325 = make_float2 (_S1324, _S1323);
    float2  _S1326 = _S987 * _S1325;
    float _S1327 = p1_3 * _S1323;
    float _S1328 = v_7 * _S1323;
    float _S1329 = p2_3 * _S1324;
    float _S1330 = v_7 * _S1324;
    float _S1331 = _S1326.x + _S1326.y;
    float _S1332 = r2_7 * _S1331;
    float _S1333 = r2_7 * _S1332;
    float _S1334 = r2_7 * _S1333;
    float _S1335 = sy1_3 * _S1323 + _S1327 + sx1_3 * _S1324 + _S1329 + _S990 * _S1331 + _S989 * _S1332 + _S988 * _S1333 + k4_3 * _S1334;
    float _S1336 = v_7 * _S1335;
    float _S1337 = u_7 * _S1335;
    float2  _S1338 = make_float2 (r2_7 * _S1324 + _S1286, r2_7 * _S1323 + _S1189);
    float2  _S1339 = make_float2 (_S999 * _S1323 + 2.0f * (u_7 * _S1330 + _S1298) + _S1187, 2.0f * (u_7 * _S1328 + _S1203) + _S995 * _S1324 + _S1284);
    float4  _S1340 = make_float4 (_S1332 + _S1299, _S1333 + _S1303, _S1334 + _S1300, r2_7 * _S1334 + _S1304);
    float3  _S1341 = _S1321 + make_float3 (_S1322.differential_0.x, _S1322.differential_0.y, 0.0f) + _S1223;
    float _S1342 = _S1001 * dpmean2d_1.x + _S1301;
    float _S1343 = _S1002 * dpmean2d_1.y + _S1202;
    float2  _S1344 = _S991 * _S1325 + make_float2 (_S996 * _S1328 + _S994 * _S1329 + 2.0f * (u_7 * _S1329) + _S992 * _S1330 + _S1337 + _S1337, _S998 * _S1327 + 2.0f * (v_7 * _S1327) + _S997 * _S1323 + _S993 * _S1324 + _S1336 + _S1336);
    CameraDistortion_0 _S1345 = _S1117;
    (&_S1345)->thin_prism_coeffs_0 = _S1338;
    (&_S1345)->tangential_coeffs_0 = _S1339;
    (&_S1345)->radial_coeffs_0 = _S1340;
    CameraDistortion_0 _S1346 = _S1117;
    CameraDistortion_0 _S1347 = _S1345;
    CameraDistortion_0 _S1348 = CameraDistortion_x24_syn_dadd_0(&_S1346, &_S1347);
    float2  _S1349 = _S975 * _S1344;
    float2  _S1350 = _S986 * _S1344;
    float _S1351 = _S1349.x + _S1349.y;
    if(_S979)
    {
        float _S1352 = _S1351 / _S981;
        float _S1353 = _S982 * - _S1352;
        float _S1354 = _S978 * (0.3333333432674408f * - (_S977 * _S1352));
        k_2 = _S1354 + _S1354;
        _S980 = _S1353;
        _S981 = 0.0f;
    }
    else
    {
        float _S1355 = _S1351 / _S980;
        float _S1356 = _S978 * - _S1355;
        k_2 = _S976 * _S1355;
        _S980 = 0.0f;
        _S981 = _S1356;
    }
    DiffPair_float_0 _S1357;
    (&_S1357)->primal_0 = _S976;
    (&_S1357)->differential_0 = 0.0f;
    DiffPair_float_0 _S1358;
    (&_S1358)->primal_0 = _S977;
    (&_S1358)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1357, &_S1358, k_2);
    float _S1359 = _S1358.differential_0 + _S980;
    float _S1360 = _S1357.differential_0 + _S981;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1361;
    (&_S1361)->primal_0 = _S975;
    (&_S1361)->differential_0 = _S973;
    s_bwd_length_impl_1(&_S1361, _S1360);
    float2  _S1362 = _S1361.differential_0 + _S1350;
    float4  _S1363 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1364;
    (&_S1364)->primal_0 = _S970.primal_0;
    (&_S1364)->differential_0 = _S1363;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1365;
    (&_S1365)->primal_0 = _S971.primal_0;
    (&_S1365)->differential_0 = _S973;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1366;
    (&_S1366)->primal_0 = _S972.primal_0;
    (&_S1366)->differential_0 = _S973;
    CameraDistortion_0 _S1367 = _S1348;
    s_bwd_prop_CameraDistortion_x24init_0(&_S1364, &_S1365, &_S1366, &_S1367);
    dpthin_prism_coeffs_3->primal_0 = (*dpthin_prism_coeffs_3).primal_0;
    dpthin_prism_coeffs_3->differential_0 = _S1366.differential_0;
    dptangential_coeffs_3->primal_0 = (*dptangential_coeffs_3).primal_0;
    dptangential_coeffs_3->differential_0 = _S1365.differential_0;
    dpradial_coeffs_3->primal_0 = (*dpradial_coeffs_3).primal_0;
    dpradial_coeffs_3->differential_0 = _S1364.differential_0;
    dpcy_1->primal_0 = (*dpcy_1).primal_0;
    dpcy_1->differential_0 = dpmean2d_1.y;
    dpcx_1->primal_0 = (*dpcx_1).primal_0;
    dpcx_1->differential_0 = dpmean2d_1.x;
    dpfy_1->primal_0 = (*dpfy_1).primal_0;
    dpfy_1->differential_0 = _S1343;
    dpfx_1->primal_0 = (*dpfx_1).primal_0;
    dpfx_1->differential_0 = _S1342;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1125.differential_0;
    float3  _S1368 = _S1341 + make_float3 (_S1362.x, _S1362.y, _S1359);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1368;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_15, float fy_15, float cx_15, float cy_15, float4  radial_coeffs_18, float2  tangential_coeffs_18, float2  thin_prism_coeffs_18, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1369 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1370 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1371 = { _S1370, _S1370 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1372 = { _S1370, _S1370, _S1371, _S1370, _S1370, _S1371 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1373;
    (&_S1373)->_S720 = _S1369;
    (&_S1373)->_S721 = _S1372;
    float3  mean_c_11 = s_primal_ctx_mul_1(R_15, mean_11) + t_14;
    float3  _S1374 = s_primal_ctx_exp_0(scale_13);
    float _S1375 = quat_14.y;
    float x2_14 = _S1375 * _S1375;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_25 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S1376 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1374.x, 0.0f, 0.0f, 0.0f, _S1374.y, 0.0f, 0.0f, 0.0f, _S1374.z);
    Matrix<float, 3, 3>  _S1377 = s_primal_ctx_mul_2(_S1376, S_1);
    Matrix<float, 3, 3>  _S1378 = transpose_0(_S1377);
    Matrix<float, 3, 3>  _S1379 = s_primal_ctx_mul_2(_S1377, _S1378);
    Matrix<float, 3, 3>  _S1380 = s_primal_ctx_mul_2(R_15, _S1379);
    Matrix<float, 3, 3>  _S1381 = transpose_0(R_15);
    Matrix<float, 3, 3>  _S1382 = s_primal_ctx_mul_2(_S1380, _S1381);
    Matrix<float, 2, 2>  _S1383 = _S1369;
    float2  _S1384 = make_float2 (0.0f);
    float2  _S1385 = _S1384;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1382, fx_15, fy_15, cx_15, cy_15, radial_coeffs_18, tangential_coeffs_18, thin_prism_coeffs_18, &_S1383, &_S1385, &(&_S1373)->_S721);
    (&_S1373)->_S720 = _S1383;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1386 = _S1373;
    float _S1387 = _S1373._S720.rows[int(0)].y * _S1373._S720.rows[int(1)].x;
    float det_orig_12 = _S1373._S720.rows[int(0)].x * _S1373._S720.rows[int(1)].y - _S1387;
    float _S1388 = _S1373._S720.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1389 = _S1373._S720;
    *&(((&_S1389)->rows + (int(0)))->x) = _S1388;
    float _S1390 = _S1373._S720.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1389)->rows + (int(1)))->y) = _S1390;
    Matrix<float, 2, 2>  _S1391 = _S1389;
    Matrix<float, 2, 2>  _S1392 = _S1389;
    float det_blur_7 = _S1388 * _S1390 - _S1387;
    float _S1393 = det_orig_12 / det_blur_7;
    float _S1394 = det_blur_7 * det_blur_7;
    float _S1395 = s_primal_ctx_max_0(0.0f, _S1393);
    float _S1396 = s_primal_ctx_sqrt_0(_S1395);
    float invdet_8 = 1.0f / det_blur_7;
    float _S1397 = - _S1373._S720.rows[int(0)].y;
    float _S1398 = - _S1373._S720.rows[int(1)].x;
    float _S1399 = - in_opacity_11;
    float _S1400 = 1.0f + s_primal_ctx_exp_1(_S1399);
    float _S1401 = 1.0f / _S1400;
    float _S1402 = _S1400 * _S1400;
    float _S1403;
    if(antialiased_11)
    {
        _S1403 = _S1401 * _S1396;
    }
    else
    {
        _S1403 = _S1401;
    }
    float _S1404 = _S1403 / 0.00392156885936856f;
    float _S1405 = 2.0f * s_primal_ctx_log_0(_S1404);
    float _S1406 = s_primal_ctx_sqrt_0(_S1405);
    float _S1407 = _S1391.rows[int(0)].x;
    float _S1408 = _S1392.rows[int(1)].y;
    float3  _S1409 = mean_11 - - s_primal_ctx_mul_1(_S1381, t_14);
    float _S1410 = _S1409.x;
    float _S1411 = _S1409.y;
    float _S1412 = _S1409.z;
    float _S1413 = _S1410 * _S1410 + _S1411 * _S1411 + _S1412 * _S1412;
    float _S1414 = s_primal_ctx_sqrt_0(_S1413);
    float x_35 = _S1410 / _S1414;
    float3  _S1415 = make_float3 (x_35);
    float _S1416 = _S1414 * _S1414;
    float y_14 = _S1411 / _S1414;
    float z_11 = _S1412 / _S1414;
    float3  _S1417 = make_float3 (z_11);
    float _S1418 = - y_14;
    float3  _S1419 = make_float3 (_S1418);
    float z2_26 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_35 * x_35 - y_14 * y_14;
    float _S1420 = 2.0f * x_35;
    float fS1_11 = _S1420 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_26 - 0.31539157032966614f;
    float3  _S1421 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_35;
    float3  _S1422 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S1423 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S1424 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S1425 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S1426 = 1.86588168144226074f * z2_26 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S1426;
    float3  _S1427 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_35;
    float3  _S1428 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S1429 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S1430 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S1431 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_35 * fC1_11 - y_14 * fS1_11);
    float3  _S1432 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_35 * fS1_11 + y_14 * fC1_11);
    float3  _S1433 = make_float3 (pSH9_1);
    float3  _S1434 = make_float3 (0.0f);
    float3  _S1435 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1436;
    (&_S1436)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1418) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_35) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S1436)->differential_0 = _S1435;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1437;
    (&_S1437)->primal_0 = _S1434;
    (&_S1437)->differential_0 = _S1435;
    s_bwd_prop_max_0(&_S1436, &_S1437, v_rgb_1);
    float3  _S1438 = _S1432 * _S1436.differential_0;
    float3  _S1439 = (*sh_coeffs_11)[int(15)] * _S1436.differential_0;
    float3  _S1440 = _S1430 * _S1436.differential_0;
    float3  _S1441 = (*sh_coeffs_11)[int(14)] * _S1436.differential_0;
    float3  _S1442 = _S1428 * _S1436.differential_0;
    float3  _S1443 = (*sh_coeffs_11)[int(13)] * _S1436.differential_0;
    float3  _S1444 = _S1427 * _S1436.differential_0;
    float3  _S1445 = (*sh_coeffs_11)[int(12)] * _S1436.differential_0;
    float3  _S1446 = _S1429 * _S1436.differential_0;
    float3  _S1447 = (*sh_coeffs_11)[int(11)] * _S1436.differential_0;
    float3  _S1448 = _S1431 * _S1436.differential_0;
    float3  _S1449 = (*sh_coeffs_11)[int(10)] * _S1436.differential_0;
    float3  _S1450 = _S1433 * _S1436.differential_0;
    float3  _S1451 = (*sh_coeffs_11)[int(9)] * _S1436.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1451.x + _S1451.y + _S1451.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1439.x + _S1439.y + _S1439.z);
    float _S1452 = _S1449.x + _S1449.y + _S1449.z;
    float _S1453 = _S1441.x + _S1441.y + _S1441.z;
    float _S1454 = _S1447.x + _S1447.y + _S1447.z;
    float _S1455 = _S1443.x + _S1443.y + _S1443.z;
    float _S1456 = _S1445.x + _S1445.y + _S1445.z;
    float _S1457 = - s_diff_fC2_T_1;
    float3  _S1458 = _S1424 * _S1436.differential_0;
    float3  _S1459 = (*sh_coeffs_11)[int(8)] * _S1436.differential_0;
    float3  _S1460 = _S1422 * _S1436.differential_0;
    float3  _S1461 = (*sh_coeffs_11)[int(7)] * _S1436.differential_0;
    float3  _S1462 = _S1421 * _S1436.differential_0;
    float3  _S1463 = (*sh_coeffs_11)[int(6)] * _S1436.differential_0;
    float3  _S1464 = _S1423 * _S1436.differential_0;
    float3  _S1465 = (*sh_coeffs_11)[int(5)] * _S1436.differential_0;
    float3  _S1466 = _S1425 * _S1436.differential_0;
    float3  _S1467 = (*sh_coeffs_11)[int(4)] * _S1436.differential_0;
    float _S1468 = _S1465.x + _S1465.y + _S1465.z;
    float _S1469 = _S1461.x + _S1461.y + _S1461.z;
    float _S1470 = fTmp1B_11 * _S1452 + x_35 * s_diff_fS2_T_1 + y_14 * _S1457 + 0.54627424478530884f * (_S1467.x + _S1467.y + _S1467.z);
    float _S1471 = fTmp1B_11 * _S1453 + y_14 * s_diff_fS2_T_1 + x_35 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1459.x + _S1459.y + _S1459.z);
    float _S1472 = y_14 * - _S1471;
    float _S1473 = x_35 * _S1471;
    float _S1474 = z_11 * (1.86588168144226074f * (z_11 * _S1456) + -2.28522896766662598f * (y_14 * _S1454 + x_35 * _S1455) + 0.94617468118667603f * (_S1463.x + _S1463.y + _S1463.z));
    float3  _S1475 = make_float3 (0.48860251903533936f) * _S1436.differential_0;
    float3  _S1476 = - _S1475;
    float3  _S1477 = _S1415 * _S1476;
    float3  _S1478 = (*sh_coeffs_11)[int(3)] * _S1476;
    float3  _S1479 = _S1417 * _S1475;
    float3  _S1480 = (*sh_coeffs_11)[int(2)] * _S1475;
    float3  _S1481 = _S1419 * _S1475;
    float3  _S1482 = (*sh_coeffs_11)[int(1)] * _S1475;
    float _S1483 = (_S1426 * _S1456 + 1.44530570507049561f * (fS1_11 * _S1452 + fC1_11 * _S1453) + -1.09254848957061768f * (y_14 * _S1468 + x_35 * _S1469) + _S1474 + _S1474 + _S1480.x + _S1480.y + _S1480.z) / _S1416;
    float _S1484 = _S1414 * _S1483;
    float _S1485 = (fTmp0C_11 * _S1454 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S1457 + fTmp0B_11 * _S1468 + _S1420 * _S1470 + _S1472 + _S1472 + - (_S1482.x + _S1482.y + _S1482.z)) / _S1416;
    float _S1486 = _S1414 * _S1485;
    float _S1487 = (fTmp0C_11 * _S1455 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S1469 + 2.0f * (y_14 * _S1470) + _S1473 + _S1473 + _S1478.x + _S1478.y + _S1478.z) / _S1416;
    float _S1488 = _S1414 * _S1487;
    float _S1489 = _S1412 * - _S1483 + _S1411 * - _S1485 + _S1410 * - _S1487;
    DiffPair_float_0 _S1490;
    (&_S1490)->primal_0 = _S1413;
    (&_S1490)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1490, _S1489);
    float _S1491 = _S1412 * _S1490.differential_0;
    float _S1492 = _S1411 * _S1490.differential_0;
    float _S1493 = _S1410 * _S1490.differential_0;
    float3  _S1494 = make_float3 (0.282094806432724f) * _S1436.differential_0;
    float3  _S1495 = make_float3 (_S1488 + _S1493 + _S1493, _S1486 + _S1492 + _S1492, _S1484 + _S1491 + _S1491);
    float3  _S1496 = - - _S1495;
    Matrix<float, 3, 3>  _S1497 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1498;
    (&_S1498)->primal_0 = _S1381;
    (&_S1498)->differential_0 = _S1497;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1499;
    (&_S1499)->primal_0 = t_14;
    (&_S1499)->differential_0 = _S1435;
    s_bwd_prop_mul_1(&_S1498, &_S1499, _S1496);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1500 = _S1498;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1501 = _S1499;
    float2  _S1502 = _S1384;
    *&((&_S1502)->y) = v_conic_1.z;
    float2  _S1503 = _S1384;
    *&((&_S1503)->y) = v_conic_1.y;
    *&((&_S1503)->x) = v_conic_1.x;
    DiffPair_float_0 _S1504;
    (&_S1504)->primal_0 = _S1408;
    (&_S1504)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1504, 0.0f);
    DiffPair_float_0 _S1505;
    (&_S1505)->primal_0 = _S1407;
    (&_S1505)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1505, 0.0f);
    DiffPair_float_0 _S1506;
    (&_S1506)->primal_0 = 3.32999992370605469f;
    (&_S1506)->differential_0 = 0.0f;
    DiffPair_float_0 _S1507;
    (&_S1507)->primal_0 = _S1406;
    (&_S1507)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1506, &_S1507, 0.0f);
    DiffPair_float_0 _S1508;
    (&_S1508)->primal_0 = _S1405;
    (&_S1508)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1508, _S1507.differential_0);
    float _S1509 = 2.0f * _S1508.differential_0;
    DiffPair_float_0 _S1510;
    (&_S1510)->primal_0 = _S1404;
    (&_S1510)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1510, _S1509);
    float2  _S1511 = make_float2 (_S1505.differential_0, 0.0f);
    float _S1512 = v_opacity_1 + 254.9999847412109375f * _S1510.differential_0;
    FixedArray<float3 , 16>  _S1513;
    _S1513[int(0)] = _S1435;
    _S1513[int(1)] = _S1435;
    _S1513[int(2)] = _S1435;
    _S1513[int(3)] = _S1435;
    _S1513[int(4)] = _S1435;
    _S1513[int(5)] = _S1435;
    _S1513[int(6)] = _S1435;
    _S1513[int(7)] = _S1435;
    _S1513[int(8)] = _S1435;
    _S1513[int(9)] = _S1435;
    _S1513[int(10)] = _S1435;
    _S1513[int(11)] = _S1435;
    _S1513[int(12)] = _S1435;
    _S1513[int(13)] = _S1435;
    _S1513[int(14)] = _S1435;
    _S1513[int(15)] = _S1435;
    _S1513[int(7)] = _S1460;
    _S1513[int(0)] = _S1494;
    _S1513[int(1)] = _S1481;
    _S1513[int(2)] = _S1479;
    _S1513[int(3)] = _S1477;
    _S1513[int(4)] = _S1466;
    _S1513[int(5)] = _S1464;
    _S1513[int(6)] = _S1462;
    _S1513[int(15)] = _S1438;
    _S1513[int(8)] = _S1458;
    _S1513[int(9)] = _S1450;
    _S1513[int(10)] = _S1448;
    _S1513[int(11)] = _S1446;
    _S1513[int(12)] = _S1444;
    _S1513[int(13)] = _S1442;
    _S1513[int(14)] = _S1440;
    float3  _S1514 = _S1513[int(0)];
    float3  _S1515 = _S1513[int(1)];
    float3  _S1516 = _S1513[int(2)];
    float3  _S1517 = _S1513[int(3)];
    float3  _S1518 = _S1513[int(4)];
    float3  _S1519 = _S1513[int(5)];
    float3  _S1520 = _S1513[int(6)];
    float3  _S1521 = _S1513[int(7)];
    float3  _S1522 = _S1513[int(8)];
    float3  _S1523 = _S1513[int(9)];
    float3  _S1524 = _S1513[int(10)];
    float3  _S1525 = _S1513[int(11)];
    float3  _S1526 = _S1513[int(12)];
    float3  _S1527 = _S1513[int(13)];
    float3  _S1528 = _S1513[int(14)];
    float3  _S1529 = _S1513[int(15)];
    Matrix<float, 2, 2>  _S1530 = _S1369;
    _S1530[int(1)] = _S1502;
    _S1530[int(0)] = _S1503;
    Matrix<float, 2, 2>  _S1531 = _S1530;
    float2  _S1532 = make_float2 (0.0f, _S1504.differential_0);
    float _S1533;
    if(antialiased_11)
    {
        float _S1534 = _S1401 * _S1512;
        _S1403 = _S1396 * _S1512;
        _S1533 = _S1534;
    }
    else
    {
        _S1403 = _S1512;
        _S1533 = 0.0f;
    }
    float _S1535 = - (_S1403 / _S1402);
    DiffPair_float_0 _S1536;
    (&_S1536)->primal_0 = _S1399;
    (&_S1536)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1536, _S1535);
    float _S1537 = - _S1536.differential_0;
    float _S1538 = invdet_8 * _S1531.rows[int(1)].y;
    float _S1539 = - (invdet_8 * _S1531.rows[int(1)].x);
    float _S1540 = - (invdet_8 * _S1531.rows[int(0)].y);
    float _S1541 = invdet_8 * _S1531.rows[int(0)].x;
    float _S1542 = - ((_S1388 * _S1531.rows[int(1)].y + _S1398 * _S1531.rows[int(1)].x + _S1397 * _S1531.rows[int(0)].y + _S1390 * _S1531.rows[int(0)].x) / _S1394);
    DiffPair_float_0 _S1543;
    (&_S1543)->primal_0 = _S1395;
    (&_S1543)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1543, _S1533);
    DiffPair_float_0 _S1544;
    (&_S1544)->primal_0 = 0.0f;
    (&_S1544)->differential_0 = 0.0f;
    DiffPair_float_0 _S1545;
    (&_S1545)->primal_0 = _S1393;
    (&_S1545)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1544, &_S1545, _S1543.differential_0);
    float _S1546 = _S1545.differential_0 / _S1394;
    float s_diff_det_orig_T_1 = det_blur_7 * _S1546;
    float _S1547 = _S1542 + det_orig_12 * - _S1546;
    float _S1548 = - _S1547;
    float _S1549 = _S1388 * _S1547;
    float _S1550 = _S1390 * _S1547;
    Matrix<float, 2, 2>  _S1551 = _S1369;
    _S1551[int(1)] = _S1532;
    _S1551[int(0)] = _S1511;
    _S1389 = _S1551;
    *&(((&_S1389)->rows + (int(1)))->y) = 0.0f;
    float _S1552 = _S1541 + _S1549 + _S1551.rows[int(1)].y;
    *&(((&_S1389)->rows + (int(0)))->x) = 0.0f;
    float _S1553 = _S1538 + _S1550 + _S1551.rows[int(0)].x;
    float _S1554 = _S1548 + - s_diff_det_orig_T_1;
    float _S1555 = _S1539 + _S1386._S720.rows[int(0)].y * _S1554;
    float _S1556 = _S1540 + _S1386._S720.rows[int(1)].x * _S1554;
    float _S1557 = _S1386._S720.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S1558 = _S1552 + _S1386._S720.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S1559 = _S1384;
    *&((&_S1559)->x) = _S1555;
    *&((&_S1559)->y) = _S1558;
    float _S1560 = _S1553 + _S1557;
    float2  _S1561 = _S1384;
    *&((&_S1561)->y) = _S1556;
    *&((&_S1561)->x) = _S1560;
    Matrix<float, 2, 2>  _S1562 = _S1369;
    _S1562[int(1)] = _S1559;
    _S1562[int(0)] = _S1561;
    Matrix<float, 2, 2>  _S1563 = _S1389 + _S1562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1564;
    (&_S1564)->primal_0 = mean_c_11;
    (&_S1564)->differential_0 = _S1435;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1565;
    (&_S1565)->primal_0 = _S1382;
    (&_S1565)->differential_0 = _S1497;
    DiffPair_float_0 _S1566;
    (&_S1566)->primal_0 = fx_15;
    (&_S1566)->differential_0 = 0.0f;
    DiffPair_float_0 _S1567;
    (&_S1567)->primal_0 = fy_15;
    (&_S1567)->differential_0 = 0.0f;
    DiffPair_float_0 _S1568;
    (&_S1568)->primal_0 = cx_15;
    (&_S1568)->differential_0 = 0.0f;
    DiffPair_float_0 _S1569;
    (&_S1569)->primal_0 = cy_15;
    (&_S1569)->differential_0 = 0.0f;
    float4  _S1570 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1571;
    (&_S1571)->primal_0 = radial_coeffs_18;
    (&_S1571)->differential_0 = _S1570;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1572;
    (&_S1572)->primal_0 = tangential_coeffs_18;
    (&_S1572)->differential_0 = _S1384;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1573;
    (&_S1573)->primal_0 = thin_prism_coeffs_18;
    (&_S1573)->differential_0 = _S1384;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1574 = _S1386._S721;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S1564, &_S1565, &_S1566, &_S1567, &_S1568, &_S1569, &_S1571, &_S1572, &_S1573, _S1563, v_mean2d_1, &_S1574);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1575;
    (&_S1575)->primal_0 = _S1380;
    (&_S1575)->differential_0 = _S1497;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1576;
    (&_S1576)->primal_0 = _S1381;
    (&_S1576)->differential_0 = _S1497;
    s_bwd_prop_mul_4(&_S1575, &_S1576, _S1565.differential_0);
    Matrix<float, 3, 3>  _S1577 = transpose_0(_S1576.differential_0 + _S1500.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1578;
    (&_S1578)->primal_0 = R_15;
    (&_S1578)->differential_0 = _S1497;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1579;
    (&_S1579)->primal_0 = _S1379;
    (&_S1579)->differential_0 = _S1497;
    s_bwd_prop_mul_4(&_S1578, &_S1579, _S1575.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1580;
    (&_S1580)->primal_0 = _S1377;
    (&_S1580)->differential_0 = _S1497;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1581;
    (&_S1581)->primal_0 = _S1378;
    (&_S1581)->differential_0 = _S1497;
    s_bwd_prop_mul_4(&_S1580, &_S1581, _S1579.differential_0);
    Matrix<float, 3, 3>  _S1582 = _S1580.differential_0 + transpose_0(_S1581.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1583;
    (&_S1583)->primal_0 = _S1376;
    (&_S1583)->differential_0 = _S1497;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1584;
    (&_S1584)->primal_0 = S_1;
    (&_S1584)->differential_0 = _S1497;
    s_bwd_prop_mul_4(&_S1583, &_S1584, _S1582);
    Matrix<float, 3, 3>  _S1585 = transpose_0(_S1583.differential_0);
    float _S1586 = 2.0f * - _S1585.rows[int(2)].z;
    float _S1587 = 2.0f * _S1585.rows[int(2)].y;
    float _S1588 = 2.0f * _S1585.rows[int(2)].x;
    float _S1589 = 2.0f * _S1585.rows[int(1)].z;
    float _S1590 = 2.0f * - _S1585.rows[int(1)].y;
    float _S1591 = 2.0f * _S1585.rows[int(1)].x;
    float _S1592 = 2.0f * _S1585.rows[int(0)].z;
    float _S1593 = 2.0f * _S1585.rows[int(0)].y;
    float _S1594 = 2.0f * - _S1585.rows[int(0)].x;
    float _S1595 = - _S1591 + _S1593;
    float _S1596 = _S1588 + - _S1592;
    float _S1597 = - _S1587 + _S1589;
    float _S1598 = _S1587 + _S1589;
    float _S1599 = _S1588 + _S1592;
    float _S1600 = _S1591 + _S1593;
    float _S1601 = quat_14.w * (_S1590 + _S1594);
    float _S1602 = quat_14.z * (_S1586 + _S1594);
    float _S1603 = quat_14.y * (_S1586 + _S1590);
    float _S1604 = quat_14.x * _S1595 + quat_14.z * _S1598 + quat_14.y * _S1599 + _S1601 + _S1601;
    float _S1605 = quat_14.x * _S1596 + quat_14.w * _S1598 + quat_14.y * _S1600 + _S1602 + _S1602;
    float _S1606 = quat_14.x * _S1597 + quat_14.w * _S1599 + quat_14.z * _S1600 + _S1603 + _S1603;
    float _S1607 = quat_14.w * _S1595 + quat_14.z * _S1596 + quat_14.y * _S1597;
    float3  _S1608 = _S1435;
    *&((&_S1608)->z) = _S1584.differential_0.rows[int(2)].z;
    *&((&_S1608)->y) = _S1584.differential_0.rows[int(1)].y;
    *&((&_S1608)->x) = _S1584.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1609;
    (&_S1609)->primal_0 = scale_13;
    (&_S1609)->differential_0 = _S1435;
    s_bwd_prop_exp_1(&_S1609, _S1608);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1610 = _S1609;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1611;
    (&_S1611)->primal_0 = mean_c_11;
    (&_S1611)->differential_0 = _S1435;
    s_bwd_length_impl_0(&_S1611, v_depth_1);
    float3  _S1612 = _S1564.differential_0 + _S1611.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1613;
    (&_S1613)->primal_0 = R_15;
    (&_S1613)->differential_0 = _S1497;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1614;
    (&_S1614)->primal_0 = mean_11;
    (&_S1614)->differential_0 = _S1435;
    s_bwd_prop_mul_1(&_S1613, &_S1614, _S1612);
    float3  _S1615 = _S1612 + _S1501.differential_0;
    Matrix<float, 3, 3>  _S1616 = _S1577 + _S1578.differential_0 + _S1613.differential_0;
    float4  _S1617 = _S1570;
    *&((&_S1617)->w) = _S1604;
    *&((&_S1617)->z) = _S1605;
    *&((&_S1617)->y) = _S1606;
    *&((&_S1617)->x) = _S1607;
    float4  _S1618 = _S1617;
    float3  _S1619 = _S1614.differential_0 + _S1495;
    *v_mean_1 = _S1619;
    *v_quat_1 = _S1618;
    *v_scale_1 = _S1610.differential_0;
    *v_in_opacity_1 = _S1537;
    (*v_sh_coeffs_1)[int(0)] = _S1514;
    (*v_sh_coeffs_1)[int(1)] = _S1515;
    (*v_sh_coeffs_1)[int(2)] = _S1516;
    (*v_sh_coeffs_1)[int(3)] = _S1517;
    (*v_sh_coeffs_1)[int(4)] = _S1518;
    (*v_sh_coeffs_1)[int(5)] = _S1519;
    (*v_sh_coeffs_1)[int(6)] = _S1520;
    (*v_sh_coeffs_1)[int(7)] = _S1521;
    (*v_sh_coeffs_1)[int(8)] = _S1522;
    (*v_sh_coeffs_1)[int(9)] = _S1523;
    (*v_sh_coeffs_1)[int(10)] = _S1524;
    (*v_sh_coeffs_1)[int(11)] = _S1525;
    (*v_sh_coeffs_1)[int(12)] = _S1526;
    (*v_sh_coeffs_1)[int(13)] = _S1527;
    (*v_sh_coeffs_1)[int(14)] = _S1528;
    (*v_sh_coeffs_1)[int(15)] = _S1529;
    *v_R_2 = _S1616;
    *v_t_2 = _S1615;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_16, float fy_16, float cx_16, float cy_16, float4  radial_coeffs_19, float2  tangential_coeffs_19, float2  thin_prism_coeffs_19, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_1(R_16, mean_12) + t_15;
    float3  _S1620 = s_primal_ctx_exp_0(scale_14);
    float _S1621 = quat_15.y;
    float x2_15 = _S1621 * _S1621;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_27 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S1622 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S1620.x, 0.0f, 0.0f, 0.0f, _S1620.y, 0.0f, 0.0f, 0.0f, _S1620.z);
    Matrix<float, 3, 3>  _S1623 = s_primal_ctx_mul_2(_S1622, S_2);
    Matrix<float, 3, 3>  _S1624 = transpose_0(_S1623);
    Matrix<float, 3, 3>  _S1625 = s_primal_ctx_mul_2(_S1623, _S1624);
    Matrix<float, 3, 3>  _S1626 = s_primal_ctx_mul_2(R_16, _S1625);
    Matrix<float, 3, 3>  _S1627 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S1628 = s_primal_ctx_mul_2(_S1626, _S1627);
    Matrix<float, 2, 3>  J_14 = makeMatrix<float, 2, 3> (fx_16, 0.0f, 0.0f, 0.0f, fy_16, 0.0f);
    Matrix<float, 2, 3>  _S1629 = s_primal_ctx_mul_3(J_14, _S1628);
    Matrix<float, 3, 2>  _S1630 = transpose_1(J_14);
    Matrix<float, 2, 2>  _S1631 = s_primal_ctx_mul_4(_S1629, _S1630);
    float _S1632 = _S1631.rows[int(0)].y * _S1631.rows[int(1)].x;
    float det_orig_13 = _S1631.rows[int(0)].x * _S1631.rows[int(1)].y - _S1632;
    float _S1633 = _S1631.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1634 = _S1631;
    *&(((&_S1634)->rows + (int(0)))->x) = _S1633;
    float _S1635 = _S1631.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1634)->rows + (int(1)))->y) = _S1635;
    Matrix<float, 2, 2>  _S1636 = _S1634;
    Matrix<float, 2, 2>  _S1637 = _S1634;
    float det_blur_8 = _S1633 * _S1635 - _S1632;
    float _S1638 = det_orig_13 / det_blur_8;
    float _S1639 = det_blur_8 * det_blur_8;
    float _S1640 = s_primal_ctx_max_0(0.0f, _S1638);
    float _S1641 = s_primal_ctx_sqrt_0(_S1640);
    float invdet_9 = 1.0f / det_blur_8;
    float _S1642 = - _S1631.rows[int(0)].y;
    float _S1643 = - _S1631.rows[int(1)].x;
    float _S1644 = - in_opacity_12;
    float _S1645 = 1.0f + s_primal_ctx_exp_1(_S1644);
    float _S1646 = 1.0f / _S1645;
    float _S1647 = _S1645 * _S1645;
    float _S1648;
    if(antialiased_12)
    {
        _S1648 = _S1646 * _S1641;
    }
    else
    {
        _S1648 = _S1646;
    }
    float _S1649 = _S1648 / 0.00392156885936856f;
    float _S1650 = 2.0f * s_primal_ctx_log_0(_S1649);
    float _S1651 = s_primal_ctx_sqrt_0(_S1650);
    float _S1652 = _S1636.rows[int(0)].x;
    float _S1653 = _S1637.rows[int(1)].y;
    float3  _S1654 = mean_12 - - s_primal_ctx_mul_1(_S1627, t_15);
    float _S1655 = _S1654.x;
    float _S1656 = _S1654.y;
    float _S1657 = _S1654.z;
    float _S1658 = _S1655 * _S1655 + _S1656 * _S1656 + _S1657 * _S1657;
    float _S1659 = s_primal_ctx_sqrt_0(_S1658);
    float x_36 = _S1655 / _S1659;
    float3  _S1660 = make_float3 (x_36);
    float _S1661 = _S1659 * _S1659;
    float y_15 = _S1656 / _S1659;
    float z_12 = _S1657 / _S1659;
    float3  _S1662 = make_float3 (z_12);
    float _S1663 = - y_15;
    float3  _S1664 = make_float3 (_S1663);
    float z2_28 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_36 * x_36 - y_15 * y_15;
    float _S1665 = 2.0f * x_36;
    float fS1_12 = _S1665 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_28 - 0.31539157032966614f;
    float3  _S1666 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_36;
    float3  _S1667 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S1668 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S1669 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S1670 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S1671 = 1.86588168144226074f * z2_28 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S1671;
    float3  _S1672 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_36;
    float3  _S1673 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S1674 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S1675 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S1676 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_36 * fC1_12 - y_15 * fS1_12);
    float3  _S1677 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_36 * fS1_12 + y_15 * fC1_12);
    float3  _S1678 = make_float3 (pSH9_2);
    float3  _S1679 = make_float3 (0.0f);
    float3  _S1680 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1681;
    (&_S1681)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1663) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_36) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S1681)->differential_0 = _S1680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1682;
    (&_S1682)->primal_0 = _S1679;
    (&_S1682)->differential_0 = _S1680;
    s_bwd_prop_max_0(&_S1681, &_S1682, v_rgb_2);
    float3  _S1683 = _S1677 * _S1681.differential_0;
    float3  _S1684 = (*sh_coeffs_12)[int(15)] * _S1681.differential_0;
    float3  _S1685 = _S1675 * _S1681.differential_0;
    float3  _S1686 = (*sh_coeffs_12)[int(14)] * _S1681.differential_0;
    float3  _S1687 = _S1673 * _S1681.differential_0;
    float3  _S1688 = (*sh_coeffs_12)[int(13)] * _S1681.differential_0;
    float3  _S1689 = _S1672 * _S1681.differential_0;
    float3  _S1690 = (*sh_coeffs_12)[int(12)] * _S1681.differential_0;
    float3  _S1691 = _S1674 * _S1681.differential_0;
    float3  _S1692 = (*sh_coeffs_12)[int(11)] * _S1681.differential_0;
    float3  _S1693 = _S1676 * _S1681.differential_0;
    float3  _S1694 = (*sh_coeffs_12)[int(10)] * _S1681.differential_0;
    float3  _S1695 = _S1678 * _S1681.differential_0;
    float3  _S1696 = (*sh_coeffs_12)[int(9)] * _S1681.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S1696.x + _S1696.y + _S1696.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S1684.x + _S1684.y + _S1684.z);
    float _S1697 = _S1694.x + _S1694.y + _S1694.z;
    float _S1698 = _S1686.x + _S1686.y + _S1686.z;
    float _S1699 = _S1692.x + _S1692.y + _S1692.z;
    float _S1700 = _S1688.x + _S1688.y + _S1688.z;
    float _S1701 = _S1690.x + _S1690.y + _S1690.z;
    float _S1702 = - s_diff_fC2_T_2;
    float3  _S1703 = _S1669 * _S1681.differential_0;
    float3  _S1704 = (*sh_coeffs_12)[int(8)] * _S1681.differential_0;
    float3  _S1705 = _S1667 * _S1681.differential_0;
    float3  _S1706 = (*sh_coeffs_12)[int(7)] * _S1681.differential_0;
    float3  _S1707 = _S1666 * _S1681.differential_0;
    float3  _S1708 = (*sh_coeffs_12)[int(6)] * _S1681.differential_0;
    float3  _S1709 = _S1668 * _S1681.differential_0;
    float3  _S1710 = (*sh_coeffs_12)[int(5)] * _S1681.differential_0;
    float3  _S1711 = _S1670 * _S1681.differential_0;
    float3  _S1712 = (*sh_coeffs_12)[int(4)] * _S1681.differential_0;
    float _S1713 = _S1710.x + _S1710.y + _S1710.z;
    float _S1714 = _S1706.x + _S1706.y + _S1706.z;
    float _S1715 = fTmp1B_12 * _S1697 + x_36 * s_diff_fS2_T_2 + y_15 * _S1702 + 0.54627424478530884f * (_S1712.x + _S1712.y + _S1712.z);
    float _S1716 = fTmp1B_12 * _S1698 + y_15 * s_diff_fS2_T_2 + x_36 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S1704.x + _S1704.y + _S1704.z);
    float _S1717 = y_15 * - _S1716;
    float _S1718 = x_36 * _S1716;
    float _S1719 = z_12 * (1.86588168144226074f * (z_12 * _S1701) + -2.28522896766662598f * (y_15 * _S1699 + x_36 * _S1700) + 0.94617468118667603f * (_S1708.x + _S1708.y + _S1708.z));
    float3  _S1720 = make_float3 (0.48860251903533936f) * _S1681.differential_0;
    float3  _S1721 = - _S1720;
    float3  _S1722 = _S1660 * _S1721;
    float3  _S1723 = (*sh_coeffs_12)[int(3)] * _S1721;
    float3  _S1724 = _S1662 * _S1720;
    float3  _S1725 = (*sh_coeffs_12)[int(2)] * _S1720;
    float3  _S1726 = _S1664 * _S1720;
    float3  _S1727 = (*sh_coeffs_12)[int(1)] * _S1720;
    float _S1728 = (_S1671 * _S1701 + 1.44530570507049561f * (fS1_12 * _S1697 + fC1_12 * _S1698) + -1.09254848957061768f * (y_15 * _S1713 + x_36 * _S1714) + _S1719 + _S1719 + _S1725.x + _S1725.y + _S1725.z) / _S1661;
    float _S1729 = _S1659 * _S1728;
    float _S1730 = (fTmp0C_12 * _S1699 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S1702 + fTmp0B_12 * _S1713 + _S1665 * _S1715 + _S1717 + _S1717 + - (_S1727.x + _S1727.y + _S1727.z)) / _S1661;
    float _S1731 = _S1659 * _S1730;
    float _S1732 = (fTmp0C_12 * _S1700 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S1714 + 2.0f * (y_15 * _S1715) + _S1718 + _S1718 + _S1723.x + _S1723.y + _S1723.z) / _S1661;
    float _S1733 = _S1659 * _S1732;
    float _S1734 = _S1657 * - _S1728 + _S1656 * - _S1730 + _S1655 * - _S1732;
    DiffPair_float_0 _S1735;
    (&_S1735)->primal_0 = _S1658;
    (&_S1735)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1735, _S1734);
    float _S1736 = _S1657 * _S1735.differential_0;
    float _S1737 = _S1656 * _S1735.differential_0;
    float _S1738 = _S1655 * _S1735.differential_0;
    float3  _S1739 = make_float3 (0.282094806432724f) * _S1681.differential_0;
    float3  _S1740 = make_float3 (_S1733 + _S1738 + _S1738, _S1731 + _S1737 + _S1737, _S1729 + _S1736 + _S1736);
    float3  _S1741 = - - _S1740;
    Matrix<float, 3, 3>  _S1742 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1743;
    (&_S1743)->primal_0 = _S1627;
    (&_S1743)->differential_0 = _S1742;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1744;
    (&_S1744)->primal_0 = t_15;
    (&_S1744)->differential_0 = _S1680;
    s_bwd_prop_mul_1(&_S1743, &_S1744, _S1741);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1745 = _S1743;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1746 = _S1744;
    float2  _S1747 = make_float2 (0.0f);
    float2  _S1748 = _S1747;
    *&((&_S1748)->y) = v_conic_2.z;
    float2  _S1749 = _S1747;
    *&((&_S1749)->y) = v_conic_2.y;
    *&((&_S1749)->x) = v_conic_2.x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1750;
    (&_S1750)->primal_0 = mean_c_12;
    (&_S1750)->differential_0 = _S1680;
    s_bwd_length_impl_0(&_S1750, v_depth_2);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1751 = _S1750;
    DiffPair_float_0 _S1752;
    (&_S1752)->primal_0 = _S1653;
    (&_S1752)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1752, 0.0f);
    DiffPair_float_0 _S1753;
    (&_S1753)->primal_0 = _S1652;
    (&_S1753)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1753, 0.0f);
    DiffPair_float_0 _S1754;
    (&_S1754)->primal_0 = 3.32999992370605469f;
    (&_S1754)->differential_0 = 0.0f;
    DiffPair_float_0 _S1755;
    (&_S1755)->primal_0 = _S1651;
    (&_S1755)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1754, &_S1755, 0.0f);
    DiffPair_float_0 _S1756;
    (&_S1756)->primal_0 = _S1650;
    (&_S1756)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1756, _S1755.differential_0);
    float _S1757 = 2.0f * _S1756.differential_0;
    DiffPair_float_0 _S1758;
    (&_S1758)->primal_0 = _S1649;
    (&_S1758)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1758, _S1757);
    float _S1759 = v_opacity_2 + 254.9999847412109375f * _S1758.differential_0;
    Matrix<float, 2, 2>  _S1760 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1761 = _S1760;
    _S1761[int(1)] = _S1748;
    _S1761[int(0)] = _S1749;
    Matrix<float, 2, 2>  _S1762 = _S1761;
    FixedArray<float3 , 16>  _S1763;
    _S1763[int(0)] = _S1680;
    _S1763[int(1)] = _S1680;
    _S1763[int(2)] = _S1680;
    _S1763[int(3)] = _S1680;
    _S1763[int(4)] = _S1680;
    _S1763[int(5)] = _S1680;
    _S1763[int(6)] = _S1680;
    _S1763[int(7)] = _S1680;
    _S1763[int(8)] = _S1680;
    _S1763[int(9)] = _S1680;
    _S1763[int(10)] = _S1680;
    _S1763[int(11)] = _S1680;
    _S1763[int(12)] = _S1680;
    _S1763[int(13)] = _S1680;
    _S1763[int(14)] = _S1680;
    _S1763[int(15)] = _S1680;
    _S1763[int(7)] = _S1705;
    _S1763[int(0)] = _S1739;
    _S1763[int(1)] = _S1726;
    _S1763[int(2)] = _S1724;
    _S1763[int(3)] = _S1722;
    _S1763[int(4)] = _S1711;
    _S1763[int(5)] = _S1709;
    _S1763[int(6)] = _S1707;
    _S1763[int(15)] = _S1683;
    _S1763[int(8)] = _S1703;
    _S1763[int(9)] = _S1695;
    _S1763[int(10)] = _S1693;
    _S1763[int(11)] = _S1691;
    _S1763[int(12)] = _S1689;
    _S1763[int(13)] = _S1687;
    _S1763[int(14)] = _S1685;
    float3  _S1764 = _S1763[int(0)];
    float3  _S1765 = _S1763[int(1)];
    float3  _S1766 = _S1763[int(2)];
    float3  _S1767 = _S1763[int(3)];
    float3  _S1768 = _S1763[int(4)];
    float3  _S1769 = _S1763[int(5)];
    float3  _S1770 = _S1763[int(6)];
    float3  _S1771 = _S1763[int(7)];
    float3  _S1772 = _S1763[int(8)];
    float3  _S1773 = _S1763[int(9)];
    float3  _S1774 = _S1763[int(10)];
    float3  _S1775 = _S1763[int(11)];
    float3  _S1776 = _S1763[int(12)];
    float3  _S1777 = _S1763[int(13)];
    float3  _S1778 = _S1763[int(14)];
    float3  _S1779 = _S1763[int(15)];
    float2  _S1780 = make_float2 (0.0f, _S1752.differential_0);
    float2  _S1781 = make_float2 (_S1753.differential_0, 0.0f);
    float _S1782;
    if(antialiased_12)
    {
        float _S1783 = _S1646 * _S1759;
        _S1648 = _S1641 * _S1759;
        _S1782 = _S1783;
    }
    else
    {
        _S1648 = _S1759;
        _S1782 = 0.0f;
    }
    float _S1784 = - (_S1648 / _S1647);
    DiffPair_float_0 _S1785;
    (&_S1785)->primal_0 = _S1644;
    (&_S1785)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1785, _S1784);
    float _S1786 = - _S1785.differential_0;
    float _S1787 = invdet_9 * _S1762.rows[int(1)].y;
    float _S1788 = - (invdet_9 * _S1762.rows[int(1)].x);
    float _S1789 = - (invdet_9 * _S1762.rows[int(0)].y);
    float _S1790 = invdet_9 * _S1762.rows[int(0)].x;
    float _S1791 = - ((_S1633 * _S1762.rows[int(1)].y + _S1643 * _S1762.rows[int(1)].x + _S1642 * _S1762.rows[int(0)].y + _S1635 * _S1762.rows[int(0)].x) / _S1639);
    DiffPair_float_0 _S1792;
    (&_S1792)->primal_0 = _S1640;
    (&_S1792)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1792, _S1782);
    DiffPair_float_0 _S1793;
    (&_S1793)->primal_0 = 0.0f;
    (&_S1793)->differential_0 = 0.0f;
    DiffPair_float_0 _S1794;
    (&_S1794)->primal_0 = _S1638;
    (&_S1794)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1793, &_S1794, _S1792.differential_0);
    float _S1795 = _S1794.differential_0 / _S1639;
    float s_diff_det_orig_T_2 = det_blur_8 * _S1795;
    float _S1796 = _S1791 + det_orig_13 * - _S1795;
    float _S1797 = - _S1796;
    float _S1798 = _S1633 * _S1796;
    float _S1799 = _S1635 * _S1796;
    Matrix<float, 2, 2>  _S1800 = _S1760;
    _S1800[int(1)] = _S1780;
    _S1800[int(0)] = _S1781;
    _S1634 = _S1800;
    *&(((&_S1634)->rows + (int(1)))->y) = 0.0f;
    float _S1801 = _S1790 + _S1798 + _S1800.rows[int(1)].y;
    *&(((&_S1634)->rows + (int(0)))->x) = 0.0f;
    float _S1802 = _S1787 + _S1799 + _S1800.rows[int(0)].x;
    float _S1803 = _S1797 + - s_diff_det_orig_T_2;
    float _S1804 = _S1788 + _S1631.rows[int(0)].y * _S1803;
    float _S1805 = _S1789 + _S1631.rows[int(1)].x * _S1803;
    float _S1806 = _S1631.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S1807 = _S1801 + _S1631.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S1808 = _S1747;
    *&((&_S1808)->x) = _S1804;
    *&((&_S1808)->y) = _S1807;
    float _S1809 = _S1802 + _S1806;
    float2  _S1810 = _S1747;
    *&((&_S1810)->y) = _S1805;
    *&((&_S1810)->x) = _S1809;
    float _S1811 = fy_16 * v_mean2d_2.y;
    float _S1812 = fx_16 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S1813 = _S1760;
    _S1813[int(1)] = _S1808;
    _S1813[int(0)] = _S1810;
    Matrix<float, 2, 2>  _S1814 = _S1634 + _S1813;
    Matrix<float, 2, 3>  _S1815 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1816;
    (&_S1816)->primal_0 = _S1629;
    (&_S1816)->differential_0 = _S1815;
    Matrix<float, 3, 2>  _S1817 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1818;
    (&_S1818)->primal_0 = _S1630;
    (&_S1818)->differential_0 = _S1817;
    s_bwd_prop_mul_2(&_S1816, &_S1818, _S1814);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1819;
    (&_S1819)->primal_0 = J_14;
    (&_S1819)->differential_0 = _S1815;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1820;
    (&_S1820)->primal_0 = _S1628;
    (&_S1820)->differential_0 = _S1742;
    s_bwd_prop_mul_3(&_S1819, &_S1820, _S1816.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1821;
    (&_S1821)->primal_0 = _S1626;
    (&_S1821)->differential_0 = _S1742;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1822;
    (&_S1822)->primal_0 = _S1627;
    (&_S1822)->differential_0 = _S1742;
    s_bwd_prop_mul_4(&_S1821, &_S1822, _S1820.differential_0);
    Matrix<float, 3, 3>  _S1823 = transpose_0(_S1822.differential_0 + _S1745.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1824;
    (&_S1824)->primal_0 = R_16;
    (&_S1824)->differential_0 = _S1742;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1825;
    (&_S1825)->primal_0 = _S1625;
    (&_S1825)->differential_0 = _S1742;
    s_bwd_prop_mul_4(&_S1824, &_S1825, _S1821.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1826;
    (&_S1826)->primal_0 = _S1623;
    (&_S1826)->differential_0 = _S1742;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1827;
    (&_S1827)->primal_0 = _S1624;
    (&_S1827)->differential_0 = _S1742;
    s_bwd_prop_mul_4(&_S1826, &_S1827, _S1825.differential_0);
    Matrix<float, 3, 3>  _S1828 = _S1826.differential_0 + transpose_0(_S1827.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1829;
    (&_S1829)->primal_0 = _S1622;
    (&_S1829)->differential_0 = _S1742;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1830;
    (&_S1830)->primal_0 = S_2;
    (&_S1830)->differential_0 = _S1742;
    s_bwd_prop_mul_4(&_S1829, &_S1830, _S1828);
    Matrix<float, 3, 3>  _S1831 = transpose_0(_S1829.differential_0);
    float _S1832 = 2.0f * - _S1831.rows[int(2)].z;
    float _S1833 = 2.0f * _S1831.rows[int(2)].y;
    float _S1834 = 2.0f * _S1831.rows[int(2)].x;
    float _S1835 = 2.0f * _S1831.rows[int(1)].z;
    float _S1836 = 2.0f * - _S1831.rows[int(1)].y;
    float _S1837 = 2.0f * _S1831.rows[int(1)].x;
    float _S1838 = 2.0f * _S1831.rows[int(0)].z;
    float _S1839 = 2.0f * _S1831.rows[int(0)].y;
    float _S1840 = 2.0f * - _S1831.rows[int(0)].x;
    float _S1841 = - _S1837 + _S1839;
    float _S1842 = _S1834 + - _S1838;
    float _S1843 = - _S1833 + _S1835;
    float _S1844 = _S1833 + _S1835;
    float _S1845 = _S1834 + _S1838;
    float _S1846 = _S1837 + _S1839;
    float _S1847 = quat_15.w * (_S1836 + _S1840);
    float _S1848 = quat_15.z * (_S1832 + _S1840);
    float _S1849 = quat_15.y * (_S1832 + _S1836);
    float _S1850 = quat_15.x * _S1841 + quat_15.z * _S1844 + quat_15.y * _S1845 + _S1847 + _S1847;
    float _S1851 = quat_15.x * _S1842 + quat_15.w * _S1844 + quat_15.y * _S1846 + _S1848 + _S1848;
    float _S1852 = quat_15.x * _S1843 + quat_15.w * _S1845 + quat_15.z * _S1846 + _S1849 + _S1849;
    float _S1853 = quat_15.w * _S1841 + quat_15.z * _S1842 + quat_15.y * _S1843;
    float3  _S1854 = _S1680;
    *&((&_S1854)->z) = _S1830.differential_0.rows[int(2)].z;
    *&((&_S1854)->y) = _S1830.differential_0.rows[int(1)].y;
    *&((&_S1854)->x) = _S1830.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1855;
    (&_S1855)->primal_0 = scale_14;
    (&_S1855)->differential_0 = _S1680;
    s_bwd_prop_exp_1(&_S1855, _S1854);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1856 = _S1855;
    float3  _S1857 = _S1680;
    *&((&_S1857)->y) = _S1811;
    *&((&_S1857)->x) = _S1812;
    float3  _S1858 = _S1751.differential_0 + _S1857;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1859;
    (&_S1859)->primal_0 = R_16;
    (&_S1859)->differential_0 = _S1742;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1860;
    (&_S1860)->primal_0 = mean_12;
    (&_S1860)->differential_0 = _S1680;
    s_bwd_prop_mul_1(&_S1859, &_S1860, _S1858);
    float3  _S1861 = _S1858 + _S1746.differential_0;
    Matrix<float, 3, 3>  _S1862 = _S1823 + _S1824.differential_0 + _S1859.differential_0;
    float4  _S1863 = make_float4 (0.0f);
    *&((&_S1863)->w) = _S1850;
    *&((&_S1863)->z) = _S1851;
    *&((&_S1863)->y) = _S1852;
    *&((&_S1863)->x) = _S1853;
    float4  _S1864 = _S1863;
    float3  _S1865 = _S1860.differential_0 + _S1740;
    *v_mean_2 = _S1865;
    *v_quat_2 = _S1864;
    *v_scale_2 = _S1856.differential_0;
    *v_in_opacity_2 = _S1786;
    (*v_sh_coeffs_2)[int(0)] = _S1764;
    (*v_sh_coeffs_2)[int(1)] = _S1765;
    (*v_sh_coeffs_2)[int(2)] = _S1766;
    (*v_sh_coeffs_2)[int(3)] = _S1767;
    (*v_sh_coeffs_2)[int(4)] = _S1768;
    (*v_sh_coeffs_2)[int(5)] = _S1769;
    (*v_sh_coeffs_2)[int(6)] = _S1770;
    (*v_sh_coeffs_2)[int(7)] = _S1771;
    (*v_sh_coeffs_2)[int(8)] = _S1772;
    (*v_sh_coeffs_2)[int(9)] = _S1773;
    (*v_sh_coeffs_2)[int(10)] = _S1774;
    (*v_sh_coeffs_2)[int(11)] = _S1775;
    (*v_sh_coeffs_2)[int(12)] = _S1776;
    (*v_sh_coeffs_2)[int(13)] = _S1777;
    (*v_sh_coeffs_2)[int(14)] = _S1778;
    (*v_sh_coeffs_2)[int(15)] = _S1779;
    *v_R_3 = _S1862;
    *v_t_3 = _S1861;
    return;
}

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_17, float fy_17, float cx_17, float cy_17, float4  radial_coeffs_20, float2  tangential_coeffs_20, float2  thin_prism_coeffs_20, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_13) + t_16;
    float3  _S1866 = s_primal_ctx_exp_0(scale_15);
    float _S1867 = quat_16.y;
    float x2_16 = _S1867 * _S1867;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_29 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S1868 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    Matrix<float, 3, 3>  S_3 = makeMatrix<float, 3, 3> (_S1866.x, 0.0f, 0.0f, 0.0f, _S1866.y, 0.0f, 0.0f, 0.0f, _S1866.z);
    Matrix<float, 3, 3>  _S1869 = s_primal_ctx_mul_2(_S1868, S_3);
    Matrix<float, 3, 3>  _S1870 = transpose_0(_S1869);
    Matrix<float, 3, 3>  _S1871 = s_primal_ctx_mul_2(_S1869, _S1870);
    Matrix<float, 3, 3>  _S1872 = s_primal_ctx_mul_2(R_17, _S1871);
    Matrix<float, 3, 3>  _S1873 = transpose_0(R_17);
    Matrix<float, 3, 3>  _S1874 = s_primal_ctx_mul_2(_S1872, _S1873);
    float _S1875 = float(image_width_13);
    float _S1876 = float(image_height_13);
    float _S1877 = 0.30000001192092896f * (0.5f * _S1875 / fx_17);
    float lim_x_pos_1 = (_S1875 - cx_17) / fx_17 + _S1877;
    float _S1878 = 0.30000001192092896f * (0.5f * _S1876 / fy_17);
    float lim_y_pos_1 = (_S1876 - cy_17) / fy_17 + _S1878;
    float rz_7 = 1.0f / mean_c_13.z;
    float _S1879 = mean_c_13.z * mean_c_13.z;
    float rz2_7 = rz_7 * rz_7;
    float _S1880 = - (cx_17 / fx_17 + _S1877);
    float _S1881 = mean_c_13.x * rz_7;
    float _S1882 = s_primal_ctx_max_0(_S1880, _S1881);
    float _S1883 = s_primal_ctx_min_0(lim_x_pos_1, _S1882);
    float _S1884 = - (cy_17 / fy_17 + _S1878);
    float _S1885 = mean_c_13.y * rz_7;
    float _S1886 = s_primal_ctx_max_0(_S1884, _S1885);
    float _S1887 = s_primal_ctx_min_0(lim_y_pos_1, _S1886);
    float _S1888 = - fx_17;
    float _S1889 = _S1888 * (mean_c_13.z * _S1883);
    float _S1890 = - fy_17;
    float _S1891 = _S1890 * (mean_c_13.z * _S1887);
    Matrix<float, 2, 3>  J_15 = makeMatrix<float, 2, 3> (fx_17 * rz_7, 0.0f, _S1889 * rz2_7, 0.0f, fy_17 * rz_7, _S1891 * rz2_7);
    Matrix<float, 2, 3>  _S1892 = s_primal_ctx_mul_3(J_15, _S1874);
    Matrix<float, 3, 2>  _S1893 = transpose_1(J_15);
    Matrix<float, 2, 2>  _S1894 = s_primal_ctx_mul_4(_S1892, _S1893);
    float _S1895 = fx_17 * mean_c_13.x;
    float _S1896 = fy_17 * mean_c_13.y;
    float _S1897 = _S1894.rows[int(0)].y * _S1894.rows[int(1)].x;
    float det_orig_14 = _S1894.rows[int(0)].x * _S1894.rows[int(1)].y - _S1897;
    float _S1898 = _S1894.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1899 = _S1894;
    *&(((&_S1899)->rows + (int(0)))->x) = _S1898;
    float _S1900 = _S1894.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1899)->rows + (int(1)))->y) = _S1900;
    Matrix<float, 2, 2>  _S1901 = _S1899;
    Matrix<float, 2, 2>  _S1902 = _S1899;
    float det_blur_9 = _S1898 * _S1900 - _S1897;
    float _S1903 = det_orig_14 / det_blur_9;
    float _S1904 = det_blur_9 * det_blur_9;
    float _S1905 = s_primal_ctx_max_0(0.0f, _S1903);
    float _S1906 = s_primal_ctx_sqrt_0(_S1905);
    float _S1907 = - in_opacity_13;
    float _S1908 = 1.0f + s_primal_ctx_exp_1(_S1907);
    float _S1909 = 1.0f / _S1908;
    float _S1910 = _S1908 * _S1908;
    float _S1911;
    if(antialiased_13)
    {
        _S1911 = _S1909 * _S1906;
    }
    else
    {
        _S1911 = _S1909;
    }
    float _S1912 = _S1911 / 0.00392156885936856f;
    float _S1913 = 2.0f * s_primal_ctx_log_0(_S1912);
    float _S1914 = s_primal_ctx_sqrt_0(_S1913);
    float _S1915 = _S1901.rows[int(0)].x;
    float _S1916 = _S1902.rows[int(1)].y;
    float3  _S1917 = - scale_15;
    float3  _S1918 = mean_13 - - s_primal_ctx_mul_1(_S1873, t_16);
    float _S1919 = _S1918.x;
    float _S1920 = _S1918.y;
    float _S1921 = _S1918.z;
    float _S1922 = _S1919 * _S1919 + _S1920 * _S1920 + _S1921 * _S1921;
    float _S1923 = s_primal_ctx_sqrt_0(_S1922);
    float x_37 = _S1919 / _S1923;
    float3  _S1924 = make_float3 (x_37);
    float _S1925 = _S1923 * _S1923;
    float y_16 = _S1920 / _S1923;
    float z_13 = _S1921 / _S1923;
    float3  _S1926 = make_float3 (z_13);
    float _S1927 = - y_16;
    float3  _S1928 = make_float3 (_S1927);
    float z2_30 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_37 * x_37 - y_16 * y_16;
    float _S1929 = 2.0f * x_37;
    float fS1_13 = _S1929 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S1930 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_37;
    float3  _S1931 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S1932 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S1933 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S1934 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S1935 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S1935;
    float3  _S1936 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_37;
    float3  _S1937 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S1938 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S1939 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S1940 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_37 * fC1_13 - y_16 * fS1_13);
    float3  _S1941 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_37 * fS1_13 + y_16 * fC1_13);
    float3  _S1942 = make_float3 (pSH9_3);
    float3  _S1943 = make_float3 (0.0f);
    float3  _S1944 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1945;
    (&_S1945)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1927) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_37) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S1945)->differential_0 = _S1944;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1946;
    (&_S1946)->primal_0 = _S1943;
    (&_S1946)->differential_0 = _S1944;
    s_bwd_prop_max_0(&_S1945, &_S1946, v_rgb_3);
    float3  _S1947 = _S1941 * _S1945.differential_0;
    float3  _S1948 = (*sh_coeffs_13)[int(15)] * _S1945.differential_0;
    float3  _S1949 = _S1939 * _S1945.differential_0;
    float3  _S1950 = (*sh_coeffs_13)[int(14)] * _S1945.differential_0;
    float3  _S1951 = _S1937 * _S1945.differential_0;
    float3  _S1952 = (*sh_coeffs_13)[int(13)] * _S1945.differential_0;
    float3  _S1953 = _S1936 * _S1945.differential_0;
    float3  _S1954 = (*sh_coeffs_13)[int(12)] * _S1945.differential_0;
    float3  _S1955 = _S1938 * _S1945.differential_0;
    float3  _S1956 = (*sh_coeffs_13)[int(11)] * _S1945.differential_0;
    float3  _S1957 = _S1940 * _S1945.differential_0;
    float3  _S1958 = (*sh_coeffs_13)[int(10)] * _S1945.differential_0;
    float3  _S1959 = _S1942 * _S1945.differential_0;
    float3  _S1960 = (*sh_coeffs_13)[int(9)] * _S1945.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S1960.x + _S1960.y + _S1960.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S1948.x + _S1948.y + _S1948.z);
    float _S1961 = _S1958.x + _S1958.y + _S1958.z;
    float _S1962 = _S1950.x + _S1950.y + _S1950.z;
    float _S1963 = _S1956.x + _S1956.y + _S1956.z;
    float _S1964 = _S1952.x + _S1952.y + _S1952.z;
    float _S1965 = _S1954.x + _S1954.y + _S1954.z;
    float _S1966 = - s_diff_fC2_T_3;
    float3  _S1967 = _S1933 * _S1945.differential_0;
    float3  _S1968 = (*sh_coeffs_13)[int(8)] * _S1945.differential_0;
    float3  _S1969 = _S1931 * _S1945.differential_0;
    float3  _S1970 = (*sh_coeffs_13)[int(7)] * _S1945.differential_0;
    float3  _S1971 = _S1930 * _S1945.differential_0;
    float3  _S1972 = (*sh_coeffs_13)[int(6)] * _S1945.differential_0;
    float3  _S1973 = _S1932 * _S1945.differential_0;
    float3  _S1974 = (*sh_coeffs_13)[int(5)] * _S1945.differential_0;
    float3  _S1975 = _S1934 * _S1945.differential_0;
    float3  _S1976 = (*sh_coeffs_13)[int(4)] * _S1945.differential_0;
    float _S1977 = _S1974.x + _S1974.y + _S1974.z;
    float _S1978 = _S1970.x + _S1970.y + _S1970.z;
    float _S1979 = fTmp1B_13 * _S1961 + x_37 * s_diff_fS2_T_3 + y_16 * _S1966 + 0.54627424478530884f * (_S1976.x + _S1976.y + _S1976.z);
    float _S1980 = fTmp1B_13 * _S1962 + y_16 * s_diff_fS2_T_3 + x_37 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S1968.x + _S1968.y + _S1968.z);
    float _S1981 = y_16 * - _S1980;
    float _S1982 = x_37 * _S1980;
    float _S1983 = z_13 * (1.86588168144226074f * (z_13 * _S1965) + -2.28522896766662598f * (y_16 * _S1963 + x_37 * _S1964) + 0.94617468118667603f * (_S1972.x + _S1972.y + _S1972.z));
    float3  _S1984 = make_float3 (0.48860251903533936f) * _S1945.differential_0;
    float3  _S1985 = - _S1984;
    float3  _S1986 = _S1924 * _S1985;
    float3  _S1987 = (*sh_coeffs_13)[int(3)] * _S1985;
    float3  _S1988 = _S1926 * _S1984;
    float3  _S1989 = (*sh_coeffs_13)[int(2)] * _S1984;
    float3  _S1990 = _S1928 * _S1984;
    float3  _S1991 = (*sh_coeffs_13)[int(1)] * _S1984;
    float _S1992 = (_S1935 * _S1965 + 1.44530570507049561f * (fS1_13 * _S1961 + fC1_13 * _S1962) + -1.09254848957061768f * (y_16 * _S1977 + x_37 * _S1978) + _S1983 + _S1983 + _S1989.x + _S1989.y + _S1989.z) / _S1925;
    float _S1993 = _S1923 * _S1992;
    float _S1994 = (fTmp0C_13 * _S1963 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S1966 + fTmp0B_13 * _S1977 + _S1929 * _S1979 + _S1981 + _S1981 + - (_S1991.x + _S1991.y + _S1991.z)) / _S1925;
    float _S1995 = _S1923 * _S1994;
    float _S1996 = (fTmp0C_13 * _S1964 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S1978 + 2.0f * (y_16 * _S1979) + _S1982 + _S1982 + _S1987.x + _S1987.y + _S1987.z) / _S1925;
    float _S1997 = _S1923 * _S1996;
    float _S1998 = _S1921 * - _S1992 + _S1920 * - _S1994 + _S1919 * - _S1996;
    DiffPair_float_0 _S1999;
    (&_S1999)->primal_0 = _S1922;
    (&_S1999)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1999, _S1998);
    float _S2000 = _S1921 * _S1999.differential_0;
    float _S2001 = _S1920 * _S1999.differential_0;
    float _S2002 = _S1919 * _S1999.differential_0;
    float3  _S2003 = make_float3 (0.282094806432724f) * _S1945.differential_0;
    float3  _S2004 = make_float3 (_S1997 + _S2002 + _S2002, _S1995 + _S2001 + _S2001, _S1993 + _S2000 + _S2000);
    float3  _S2005 = - - _S2004;
    Matrix<float, 3, 3>  _S2006 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2007;
    (&_S2007)->primal_0 = _S1873;
    (&_S2007)->differential_0 = _S2006;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2008;
    (&_S2008)->primal_0 = t_16;
    (&_S2008)->differential_0 = _S1944;
    s_bwd_prop_mul_1(&_S2007, &_S2008, _S2005);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2009 = _S2007;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2010 = _S2008;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2011;
    (&_S2011)->primal_0 = _S1917;
    (&_S2011)->differential_0 = _S1944;
    s_bwd_prop_exp_1(&_S2011, v_conic_3);
    float3  _S2012 = - _S2011.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2013;
    (&_S2013)->primal_0 = mean_c_13;
    (&_S2013)->differential_0 = _S1944;
    s_bwd_length_impl_0(&_S2013, v_depth_3);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2014 = _S2013;
    DiffPair_float_0 _S2015;
    (&_S2015)->primal_0 = _S1916;
    (&_S2015)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2015, 0.0f);
    DiffPair_float_0 _S2016;
    (&_S2016)->primal_0 = _S1915;
    (&_S2016)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2016, 0.0f);
    DiffPair_float_0 _S2017;
    (&_S2017)->primal_0 = 3.32999992370605469f;
    (&_S2017)->differential_0 = 0.0f;
    DiffPair_float_0 _S2018;
    (&_S2018)->primal_0 = _S1914;
    (&_S2018)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2017, &_S2018, 0.0f);
    DiffPair_float_0 _S2019;
    (&_S2019)->primal_0 = _S1913;
    (&_S2019)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2019, _S2018.differential_0);
    float _S2020 = 2.0f * _S2019.differential_0;
    DiffPair_float_0 _S2021;
    (&_S2021)->primal_0 = _S1912;
    (&_S2021)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2021, _S2020);
    float _S2022 = v_opacity_3 + 254.9999847412109375f * _S2021.differential_0;
    FixedArray<float3 , 16>  _S2023;
    _S2023[int(0)] = _S1944;
    _S2023[int(1)] = _S1944;
    _S2023[int(2)] = _S1944;
    _S2023[int(3)] = _S1944;
    _S2023[int(4)] = _S1944;
    _S2023[int(5)] = _S1944;
    _S2023[int(6)] = _S1944;
    _S2023[int(7)] = _S1944;
    _S2023[int(8)] = _S1944;
    _S2023[int(9)] = _S1944;
    _S2023[int(10)] = _S1944;
    _S2023[int(11)] = _S1944;
    _S2023[int(12)] = _S1944;
    _S2023[int(13)] = _S1944;
    _S2023[int(14)] = _S1944;
    _S2023[int(15)] = _S1944;
    _S2023[int(7)] = _S1969;
    _S2023[int(0)] = _S2003;
    _S2023[int(1)] = _S1990;
    _S2023[int(2)] = _S1988;
    _S2023[int(3)] = _S1986;
    _S2023[int(4)] = _S1975;
    _S2023[int(5)] = _S1973;
    _S2023[int(6)] = _S1971;
    _S2023[int(15)] = _S1947;
    _S2023[int(8)] = _S1967;
    _S2023[int(9)] = _S1959;
    _S2023[int(10)] = _S1957;
    _S2023[int(11)] = _S1955;
    _S2023[int(12)] = _S1953;
    _S2023[int(13)] = _S1951;
    _S2023[int(14)] = _S1949;
    float3  _S2024 = _S2023[int(0)];
    float3  _S2025 = _S2023[int(1)];
    float3  _S2026 = _S2023[int(2)];
    float3  _S2027 = _S2023[int(3)];
    float3  _S2028 = _S2023[int(4)];
    float3  _S2029 = _S2023[int(5)];
    float3  _S2030 = _S2023[int(6)];
    float3  _S2031 = _S2023[int(7)];
    float3  _S2032 = _S2023[int(8)];
    float3  _S2033 = _S2023[int(9)];
    float3  _S2034 = _S2023[int(10)];
    float3  _S2035 = _S2023[int(11)];
    float3  _S2036 = _S2023[int(12)];
    float3  _S2037 = _S2023[int(13)];
    float3  _S2038 = _S2023[int(14)];
    float3  _S2039 = _S2023[int(15)];
    float2  _S2040 = make_float2 (0.0f, _S2015.differential_0);
    float2  _S2041 = make_float2 (_S2016.differential_0, 0.0f);
    float _S2042;
    if(antialiased_13)
    {
        float _S2043 = _S1909 * _S2022;
        _S1911 = _S1906 * _S2022;
        _S2042 = _S2043;
    }
    else
    {
        _S1911 = _S2022;
        _S2042 = 0.0f;
    }
    float _S2044 = - (_S1911 / _S1910);
    DiffPair_float_0 _S2045;
    (&_S2045)->primal_0 = _S1907;
    (&_S2045)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2045, _S2044);
    float _S2046 = - _S2045.differential_0;
    DiffPair_float_0 _S2047;
    (&_S2047)->primal_0 = _S1905;
    (&_S2047)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2047, _S2042);
    DiffPair_float_0 _S2048;
    (&_S2048)->primal_0 = 0.0f;
    (&_S2048)->differential_0 = 0.0f;
    DiffPair_float_0 _S2049;
    (&_S2049)->primal_0 = _S1903;
    (&_S2049)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2048, &_S2049, _S2047.differential_0);
    float _S2050 = _S2049.differential_0 / _S1904;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2050;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2050;
    float _S2051 = - s_diff_det_blur_T_0;
    float _S2052 = _S1898 * s_diff_det_blur_T_0;
    float _S2053 = _S1900 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2054 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2055 = _S2054;
    _S2055[int(1)] = _S2040;
    _S2055[int(0)] = _S2041;
    _S1899 = _S2055;
    *&(((&_S1899)->rows + (int(1)))->y) = 0.0f;
    float _S2056 = _S2052 + _S2055.rows[int(1)].y;
    *&(((&_S1899)->rows + (int(0)))->x) = 0.0f;
    float _S2057 = _S2053 + _S2055.rows[int(0)].x;
    float _S2058 = _S2051 + - s_diff_det_orig_T_3;
    float _S2059 = _S1894.rows[int(0)].y * _S2058;
    float _S2060 = _S1894.rows[int(1)].x * _S2058;
    float _S2061 = _S1894.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2062 = _S2056 + _S1894.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2063 = make_float2 (0.0f);
    float2  _S2064 = _S2063;
    *&((&_S2064)->x) = _S2059;
    *&((&_S2064)->y) = _S2062;
    float _S2065 = _S2057 + _S2061;
    float2  _S2066 = _S2063;
    *&((&_S2066)->y) = _S2060;
    *&((&_S2066)->x) = _S2065;
    float _S2067 = _S1896 * v_mean2d_3.y;
    float _S2068 = fy_17 * (rz_7 * v_mean2d_3.y);
    float _S2069 = _S1895 * v_mean2d_3.x;
    float _S2070 = fx_17 * (rz_7 * v_mean2d_3.x);
    Matrix<float, 2, 2>  _S2071 = _S2054;
    _S2071[int(1)] = _S2064;
    _S2071[int(0)] = _S2066;
    Matrix<float, 2, 2>  _S2072 = _S1899 + _S2071;
    Matrix<float, 2, 3>  _S2073 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2074;
    (&_S2074)->primal_0 = _S1892;
    (&_S2074)->differential_0 = _S2073;
    Matrix<float, 3, 2>  _S2075 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2076;
    (&_S2076)->primal_0 = _S1893;
    (&_S2076)->differential_0 = _S2075;
    s_bwd_prop_mul_2(&_S2074, &_S2076, _S2072);
    Matrix<float, 2, 3>  _S2077 = transpose_2(_S2076.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2078;
    (&_S2078)->primal_0 = J_15;
    (&_S2078)->differential_0 = _S2073;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2079;
    (&_S2079)->primal_0 = _S1874;
    (&_S2079)->differential_0 = _S2006;
    s_bwd_prop_mul_3(&_S2078, &_S2079, _S2074.differential_0);
    Matrix<float, 2, 3>  _S2080 = _S2077 + _S2078.differential_0;
    float _S2081 = _S1891 * _S2080.rows[int(1)].z;
    float s_diff_ty_T_1 = _S1890 * (rz2_7 * _S2080.rows[int(1)].z);
    float _S2082 = fy_17 * _S2080.rows[int(1)].y;
    float _S2083 = _S1889 * _S2080.rows[int(0)].z;
    float s_diff_tx_T_1 = _S1888 * (rz2_7 * _S2080.rows[int(0)].z);
    float _S2084 = fx_17 * _S2080.rows[int(0)].x;
    float _S2085 = mean_c_13.z * s_diff_ty_T_1;
    float _S2086 = _S1887 * s_diff_ty_T_1;
    DiffPair_float_0 _S2087;
    (&_S2087)->primal_0 = lim_y_pos_1;
    (&_S2087)->differential_0 = 0.0f;
    DiffPair_float_0 _S2088;
    (&_S2088)->primal_0 = _S1886;
    (&_S2088)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2087, &_S2088, _S2085);
    DiffPair_float_0 _S2089;
    (&_S2089)->primal_0 = _S1884;
    (&_S2089)->differential_0 = 0.0f;
    DiffPair_float_0 _S2090;
    (&_S2090)->primal_0 = _S1885;
    (&_S2090)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2089, &_S2090, _S2088.differential_0);
    float _S2091 = mean_c_13.y * _S2090.differential_0;
    float _S2092 = rz_7 * _S2090.differential_0;
    float _S2093 = mean_c_13.z * s_diff_tx_T_1;
    float _S2094 = _S1883 * s_diff_tx_T_1;
    DiffPair_float_0 _S2095;
    (&_S2095)->primal_0 = lim_x_pos_1;
    (&_S2095)->differential_0 = 0.0f;
    DiffPair_float_0 _S2096;
    (&_S2096)->primal_0 = _S1882;
    (&_S2096)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2095, &_S2096, _S2093);
    DiffPair_float_0 _S2097;
    (&_S2097)->primal_0 = _S1880;
    (&_S2097)->differential_0 = 0.0f;
    DiffPair_float_0 _S2098;
    (&_S2098)->primal_0 = _S1881;
    (&_S2098)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2097, &_S2098, _S2096.differential_0);
    float _S2099 = rz_7 * (_S2081 + _S2083);
    float _S2100 = _S2086 + _S2094 + - ((_S2067 + _S2069 + _S2082 + _S2084 + _S2091 + mean_c_13.x * _S2098.differential_0 + _S2099 + _S2099) / _S1879);
    float _S2101 = _S2068 + _S2092;
    float _S2102 = _S2070 + rz_7 * _S2098.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2103;
    (&_S2103)->primal_0 = _S1872;
    (&_S2103)->differential_0 = _S2006;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2104;
    (&_S2104)->primal_0 = _S1873;
    (&_S2104)->differential_0 = _S2006;
    s_bwd_prop_mul_4(&_S2103, &_S2104, _S2079.differential_0);
    Matrix<float, 3, 3>  _S2105 = transpose_0(_S2104.differential_0 + _S2009.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2106;
    (&_S2106)->primal_0 = R_17;
    (&_S2106)->differential_0 = _S2006;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2107;
    (&_S2107)->primal_0 = _S1871;
    (&_S2107)->differential_0 = _S2006;
    s_bwd_prop_mul_4(&_S2106, &_S2107, _S2103.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2108;
    (&_S2108)->primal_0 = _S1869;
    (&_S2108)->differential_0 = _S2006;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2109;
    (&_S2109)->primal_0 = _S1870;
    (&_S2109)->differential_0 = _S2006;
    s_bwd_prop_mul_4(&_S2108, &_S2109, _S2107.differential_0);
    Matrix<float, 3, 3>  _S2110 = _S2108.differential_0 + transpose_0(_S2109.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2111;
    (&_S2111)->primal_0 = _S1868;
    (&_S2111)->differential_0 = _S2006;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2112;
    (&_S2112)->primal_0 = S_3;
    (&_S2112)->differential_0 = _S2006;
    s_bwd_prop_mul_4(&_S2111, &_S2112, _S2110);
    Matrix<float, 3, 3>  _S2113 = transpose_0(_S2111.differential_0);
    float _S2114 = 2.0f * - _S2113.rows[int(2)].z;
    float _S2115 = 2.0f * _S2113.rows[int(2)].y;
    float _S2116 = 2.0f * _S2113.rows[int(2)].x;
    float _S2117 = 2.0f * _S2113.rows[int(1)].z;
    float _S2118 = 2.0f * - _S2113.rows[int(1)].y;
    float _S2119 = 2.0f * _S2113.rows[int(1)].x;
    float _S2120 = 2.0f * _S2113.rows[int(0)].z;
    float _S2121 = 2.0f * _S2113.rows[int(0)].y;
    float _S2122 = 2.0f * - _S2113.rows[int(0)].x;
    float _S2123 = - _S2119 + _S2121;
    float _S2124 = _S2116 + - _S2120;
    float _S2125 = - _S2115 + _S2117;
    float _S2126 = _S2115 + _S2117;
    float _S2127 = _S2116 + _S2120;
    float _S2128 = _S2119 + _S2121;
    float _S2129 = quat_16.w * (_S2118 + _S2122);
    float _S2130 = quat_16.z * (_S2114 + _S2122);
    float _S2131 = quat_16.y * (_S2114 + _S2118);
    float _S2132 = quat_16.x * _S2123 + quat_16.z * _S2126 + quat_16.y * _S2127 + _S2129 + _S2129;
    float _S2133 = quat_16.x * _S2124 + quat_16.w * _S2126 + quat_16.y * _S2128 + _S2130 + _S2130;
    float _S2134 = quat_16.x * _S2125 + quat_16.w * _S2127 + quat_16.z * _S2128 + _S2131 + _S2131;
    float _S2135 = quat_16.w * _S2123 + quat_16.z * _S2124 + quat_16.y * _S2125;
    float3  _S2136 = _S1944;
    *&((&_S2136)->z) = _S2112.differential_0.rows[int(2)].z;
    *&((&_S2136)->y) = _S2112.differential_0.rows[int(1)].y;
    *&((&_S2136)->x) = _S2112.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2137;
    (&_S2137)->primal_0 = scale_15;
    (&_S2137)->differential_0 = _S1944;
    s_bwd_prop_exp_1(&_S2137, _S2136);
    float3  _S2138 = _S1944;
    *&((&_S2138)->z) = _S2100;
    *&((&_S2138)->y) = _S2101;
    *&((&_S2138)->x) = _S2102;
    float3  _S2139 = _S2014.differential_0 + _S2138;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2140;
    (&_S2140)->primal_0 = R_17;
    (&_S2140)->differential_0 = _S2006;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2141;
    (&_S2141)->primal_0 = mean_13;
    (&_S2141)->differential_0 = _S1944;
    s_bwd_prop_mul_1(&_S2140, &_S2141, _S2139);
    float3  _S2142 = _S2139 + _S2010.differential_0;
    Matrix<float, 3, 3>  _S2143 = _S2105 + _S2106.differential_0 + _S2140.differential_0;
    float3  _S2144 = _S2137.differential_0 + _S2012;
    float4  _S2145 = make_float4 (0.0f);
    *&((&_S2145)->w) = _S2132;
    *&((&_S2145)->z) = _S2133;
    *&((&_S2145)->y) = _S2134;
    *&((&_S2145)->x) = _S2135;
    float4  _S2146 = _S2145;
    float3  _S2147 = _S2141.differential_0 + _S2004;
    *v_mean_3 = _S2147;
    *v_quat_3 = _S2146;
    *v_scale_3 = _S2144;
    *v_in_opacity_3 = _S2046;
    (*v_sh_coeffs_3)[int(0)] = _S2024;
    (*v_sh_coeffs_3)[int(1)] = _S2025;
    (*v_sh_coeffs_3)[int(2)] = _S2026;
    (*v_sh_coeffs_3)[int(3)] = _S2027;
    (*v_sh_coeffs_3)[int(4)] = _S2028;
    (*v_sh_coeffs_3)[int(5)] = _S2029;
    (*v_sh_coeffs_3)[int(6)] = _S2030;
    (*v_sh_coeffs_3)[int(7)] = _S2031;
    (*v_sh_coeffs_3)[int(8)] = _S2032;
    (*v_sh_coeffs_3)[int(9)] = _S2033;
    (*v_sh_coeffs_3)[int(10)] = _S2034;
    (*v_sh_coeffs_3)[int(11)] = _S2035;
    (*v_sh_coeffs_3)[int(12)] = _S2036;
    (*v_sh_coeffs_3)[int(13)] = _S2037;
    (*v_sh_coeffs_3)[int(14)] = _S2038;
    (*v_sh_coeffs_3)[int(15)] = _S2039;
    *v_R_4 = _S2143;
    *v_t_4 = _S2142;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2148;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2149;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_18, float fy_18, float cx_18, float cy_18, float4  radial_coeffs_21, float2  tangential_coeffs_21, float2  thin_prism_coeffs_21, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2150 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S2151 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S2152 = { _S2151, _S2151 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2153 = { _S2151, _S2151, _S2152, _S2151, _S2151, _S2152 };
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2154;
    (&_S2154)->_S2148 = _S2150;
    (&_S2154)->_S2149 = _S2153;
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_14) + t_17;
    float3  _S2155 = s_primal_ctx_exp_0(scale_16);
    float _S2156 = quat_17.y;
    float x2_17 = _S2156 * _S2156;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_31 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2157 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17)));
    Matrix<float, 3, 3>  S_4 = makeMatrix<float, 3, 3> (_S2155.x, 0.0f, 0.0f, 0.0f, _S2155.y, 0.0f, 0.0f, 0.0f, _S2155.z);
    Matrix<float, 3, 3>  _S2158 = s_primal_ctx_mul_2(_S2157, S_4);
    Matrix<float, 3, 3>  _S2159 = transpose_0(_S2158);
    Matrix<float, 3, 3>  _S2160 = s_primal_ctx_mul_2(_S2158, _S2159);
    Matrix<float, 3, 3>  _S2161 = s_primal_ctx_mul_2(R_18, _S2160);
    Matrix<float, 3, 3>  _S2162 = transpose_0(R_18);
    Matrix<float, 3, 3>  _S2163 = s_primal_ctx_mul_2(_S2161, _S2162);
    Matrix<float, 2, 2>  _S2164 = _S2150;
    float2  _S2165 = make_float2 (0.0f);
    float2  _S2166 = _S2165;
    s_primal_ctx_fisheye_proj_3dgs_0(mean_c_14, _S2163, fx_18, fy_18, cx_18, cy_18, radial_coeffs_21, tangential_coeffs_21, thin_prism_coeffs_21, &_S2164, &_S2166, &(&_S2154)->_S2149);
    (&_S2154)->_S2148 = _S2164;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2167 = _S2154;
    float _S2168 = _S2154._S2148.rows[int(0)].y * _S2154._S2148.rows[int(1)].x;
    float det_orig_15 = _S2154._S2148.rows[int(0)].x * _S2154._S2148.rows[int(1)].y - _S2168;
    float _S2169 = _S2154._S2148.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2170 = _S2154._S2148;
    *&(((&_S2170)->rows + (int(0)))->x) = _S2169;
    float _S2171 = _S2154._S2148.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2170)->rows + (int(1)))->y) = _S2171;
    Matrix<float, 2, 2>  _S2172 = _S2170;
    Matrix<float, 2, 2>  _S2173 = _S2170;
    float det_blur_10 = _S2169 * _S2171 - _S2168;
    float _S2174 = det_orig_15 / det_blur_10;
    float _S2175 = det_blur_10 * det_blur_10;
    float _S2176 = s_primal_ctx_max_0(0.0f, _S2174);
    float _S2177 = s_primal_ctx_sqrt_0(_S2176);
    float _S2178 = - in_opacity_14;
    float _S2179 = 1.0f + s_primal_ctx_exp_1(_S2178);
    float _S2180 = 1.0f / _S2179;
    float _S2181 = _S2179 * _S2179;
    float _S2182;
    if(antialiased_14)
    {
        _S2182 = _S2180 * _S2177;
    }
    else
    {
        _S2182 = _S2180;
    }
    float _S2183 = _S2182 / 0.00392156885936856f;
    float _S2184 = 2.0f * s_primal_ctx_log_0(_S2183);
    float _S2185 = s_primal_ctx_sqrt_0(_S2184);
    float _S2186 = _S2172.rows[int(0)].x;
    float _S2187 = _S2173.rows[int(1)].y;
    float3  _S2188 = - scale_16;
    float3  _S2189 = mean_14 - - s_primal_ctx_mul_1(_S2162, t_17);
    float _S2190 = _S2189.x;
    float _S2191 = _S2189.y;
    float _S2192 = _S2189.z;
    float _S2193 = _S2190 * _S2190 + _S2191 * _S2191 + _S2192 * _S2192;
    float _S2194 = s_primal_ctx_sqrt_0(_S2193);
    float x_38 = _S2190 / _S2194;
    float3  _S2195 = make_float3 (x_38);
    float _S2196 = _S2194 * _S2194;
    float y_17 = _S2191 / _S2194;
    float z_14 = _S2192 / _S2194;
    float3  _S2197 = make_float3 (z_14);
    float _S2198 = - y_17;
    float3  _S2199 = make_float3 (_S2198);
    float z2_32 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_38 * x_38 - y_17 * y_17;
    float _S2200 = 2.0f * x_38;
    float fS1_14 = _S2200 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2201 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_38;
    float3  _S2202 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2203 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2204 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2205 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2206 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2206;
    float3  _S2207 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_38;
    float3  _S2208 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2209 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2210 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2211 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_38 * fC1_14 - y_17 * fS1_14);
    float3  _S2212 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_38 * fS1_14 + y_17 * fC1_14);
    float3  _S2213 = make_float3 (pSH9_4);
    float3  _S2214 = make_float3 (0.0f);
    float3  _S2215 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2216;
    (&_S2216)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2198) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_38) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2216)->differential_0 = _S2215;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2217;
    (&_S2217)->primal_0 = _S2214;
    (&_S2217)->differential_0 = _S2215;
    s_bwd_prop_max_0(&_S2216, &_S2217, v_rgb_4);
    float3  _S2218 = _S2212 * _S2216.differential_0;
    float3  _S2219 = (*sh_coeffs_14)[int(15)] * _S2216.differential_0;
    float3  _S2220 = _S2210 * _S2216.differential_0;
    float3  _S2221 = (*sh_coeffs_14)[int(14)] * _S2216.differential_0;
    float3  _S2222 = _S2208 * _S2216.differential_0;
    float3  _S2223 = (*sh_coeffs_14)[int(13)] * _S2216.differential_0;
    float3  _S2224 = _S2207 * _S2216.differential_0;
    float3  _S2225 = (*sh_coeffs_14)[int(12)] * _S2216.differential_0;
    float3  _S2226 = _S2209 * _S2216.differential_0;
    float3  _S2227 = (*sh_coeffs_14)[int(11)] * _S2216.differential_0;
    float3  _S2228 = _S2211 * _S2216.differential_0;
    float3  _S2229 = (*sh_coeffs_14)[int(10)] * _S2216.differential_0;
    float3  _S2230 = _S2213 * _S2216.differential_0;
    float3  _S2231 = (*sh_coeffs_14)[int(9)] * _S2216.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2231.x + _S2231.y + _S2231.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2219.x + _S2219.y + _S2219.z);
    float _S2232 = _S2229.x + _S2229.y + _S2229.z;
    float _S2233 = _S2221.x + _S2221.y + _S2221.z;
    float _S2234 = _S2227.x + _S2227.y + _S2227.z;
    float _S2235 = _S2223.x + _S2223.y + _S2223.z;
    float _S2236 = _S2225.x + _S2225.y + _S2225.z;
    float _S2237 = - s_diff_fC2_T_4;
    float3  _S2238 = _S2204 * _S2216.differential_0;
    float3  _S2239 = (*sh_coeffs_14)[int(8)] * _S2216.differential_0;
    float3  _S2240 = _S2202 * _S2216.differential_0;
    float3  _S2241 = (*sh_coeffs_14)[int(7)] * _S2216.differential_0;
    float3  _S2242 = _S2201 * _S2216.differential_0;
    float3  _S2243 = (*sh_coeffs_14)[int(6)] * _S2216.differential_0;
    float3  _S2244 = _S2203 * _S2216.differential_0;
    float3  _S2245 = (*sh_coeffs_14)[int(5)] * _S2216.differential_0;
    float3  _S2246 = _S2205 * _S2216.differential_0;
    float3  _S2247 = (*sh_coeffs_14)[int(4)] * _S2216.differential_0;
    float _S2248 = _S2245.x + _S2245.y + _S2245.z;
    float _S2249 = _S2241.x + _S2241.y + _S2241.z;
    float _S2250 = fTmp1B_14 * _S2232 + x_38 * s_diff_fS2_T_4 + y_17 * _S2237 + 0.54627424478530884f * (_S2247.x + _S2247.y + _S2247.z);
    float _S2251 = fTmp1B_14 * _S2233 + y_17 * s_diff_fS2_T_4 + x_38 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2239.x + _S2239.y + _S2239.z);
    float _S2252 = y_17 * - _S2251;
    float _S2253 = x_38 * _S2251;
    float _S2254 = z_14 * (1.86588168144226074f * (z_14 * _S2236) + -2.28522896766662598f * (y_17 * _S2234 + x_38 * _S2235) + 0.94617468118667603f * (_S2243.x + _S2243.y + _S2243.z));
    float3  _S2255 = make_float3 (0.48860251903533936f) * _S2216.differential_0;
    float3  _S2256 = - _S2255;
    float3  _S2257 = _S2195 * _S2256;
    float3  _S2258 = (*sh_coeffs_14)[int(3)] * _S2256;
    float3  _S2259 = _S2197 * _S2255;
    float3  _S2260 = (*sh_coeffs_14)[int(2)] * _S2255;
    float3  _S2261 = _S2199 * _S2255;
    float3  _S2262 = (*sh_coeffs_14)[int(1)] * _S2255;
    float _S2263 = (_S2206 * _S2236 + 1.44530570507049561f * (fS1_14 * _S2232 + fC1_14 * _S2233) + -1.09254848957061768f * (y_17 * _S2248 + x_38 * _S2249) + _S2254 + _S2254 + _S2260.x + _S2260.y + _S2260.z) / _S2196;
    float _S2264 = _S2194 * _S2263;
    float _S2265 = (fTmp0C_14 * _S2234 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2237 + fTmp0B_14 * _S2248 + _S2200 * _S2250 + _S2252 + _S2252 + - (_S2262.x + _S2262.y + _S2262.z)) / _S2196;
    float _S2266 = _S2194 * _S2265;
    float _S2267 = (fTmp0C_14 * _S2235 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2249 + 2.0f * (y_17 * _S2250) + _S2253 + _S2253 + _S2258.x + _S2258.y + _S2258.z) / _S2196;
    float _S2268 = _S2194 * _S2267;
    float _S2269 = _S2192 * - _S2263 + _S2191 * - _S2265 + _S2190 * - _S2267;
    DiffPair_float_0 _S2270;
    (&_S2270)->primal_0 = _S2193;
    (&_S2270)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2270, _S2269);
    float _S2271 = _S2192 * _S2270.differential_0;
    float _S2272 = _S2191 * _S2270.differential_0;
    float _S2273 = _S2190 * _S2270.differential_0;
    float3  _S2274 = make_float3 (0.282094806432724f) * _S2216.differential_0;
    float3  _S2275 = make_float3 (_S2268 + _S2273 + _S2273, _S2266 + _S2272 + _S2272, _S2264 + _S2271 + _S2271);
    float3  _S2276 = - - _S2275;
    Matrix<float, 3, 3>  _S2277 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2278;
    (&_S2278)->primal_0 = _S2162;
    (&_S2278)->differential_0 = _S2277;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2279;
    (&_S2279)->primal_0 = t_17;
    (&_S2279)->differential_0 = _S2215;
    s_bwd_prop_mul_1(&_S2278, &_S2279, _S2276);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2280 = _S2278;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2281 = _S2279;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2282;
    (&_S2282)->primal_0 = _S2188;
    (&_S2282)->differential_0 = _S2215;
    s_bwd_prop_exp_1(&_S2282, v_conic_4);
    float3  _S2283 = - _S2282.differential_0;
    DiffPair_float_0 _S2284;
    (&_S2284)->primal_0 = _S2187;
    (&_S2284)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2284, 0.0f);
    DiffPair_float_0 _S2285;
    (&_S2285)->primal_0 = _S2186;
    (&_S2285)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2285, 0.0f);
    DiffPair_float_0 _S2286;
    (&_S2286)->primal_0 = 3.32999992370605469f;
    (&_S2286)->differential_0 = 0.0f;
    DiffPair_float_0 _S2287;
    (&_S2287)->primal_0 = _S2185;
    (&_S2287)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2286, &_S2287, 0.0f);
    DiffPair_float_0 _S2288;
    (&_S2288)->primal_0 = _S2184;
    (&_S2288)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2288, _S2287.differential_0);
    float _S2289 = 2.0f * _S2288.differential_0;
    DiffPair_float_0 _S2290;
    (&_S2290)->primal_0 = _S2183;
    (&_S2290)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2290, _S2289);
    float2  _S2291 = make_float2 (_S2285.differential_0, 0.0f);
    float _S2292 = v_opacity_4 + 254.9999847412109375f * _S2290.differential_0;
    FixedArray<float3 , 16>  _S2293;
    _S2293[int(0)] = _S2215;
    _S2293[int(1)] = _S2215;
    _S2293[int(2)] = _S2215;
    _S2293[int(3)] = _S2215;
    _S2293[int(4)] = _S2215;
    _S2293[int(5)] = _S2215;
    _S2293[int(6)] = _S2215;
    _S2293[int(7)] = _S2215;
    _S2293[int(8)] = _S2215;
    _S2293[int(9)] = _S2215;
    _S2293[int(10)] = _S2215;
    _S2293[int(11)] = _S2215;
    _S2293[int(12)] = _S2215;
    _S2293[int(13)] = _S2215;
    _S2293[int(14)] = _S2215;
    _S2293[int(15)] = _S2215;
    _S2293[int(7)] = _S2240;
    _S2293[int(0)] = _S2274;
    _S2293[int(1)] = _S2261;
    _S2293[int(2)] = _S2259;
    _S2293[int(3)] = _S2257;
    _S2293[int(4)] = _S2246;
    _S2293[int(5)] = _S2244;
    _S2293[int(6)] = _S2242;
    _S2293[int(15)] = _S2218;
    _S2293[int(8)] = _S2238;
    _S2293[int(9)] = _S2230;
    _S2293[int(10)] = _S2228;
    _S2293[int(11)] = _S2226;
    _S2293[int(12)] = _S2224;
    _S2293[int(13)] = _S2222;
    _S2293[int(14)] = _S2220;
    float3  _S2294 = _S2293[int(0)];
    float3  _S2295 = _S2293[int(1)];
    float3  _S2296 = _S2293[int(2)];
    float3  _S2297 = _S2293[int(3)];
    float3  _S2298 = _S2293[int(4)];
    float3  _S2299 = _S2293[int(5)];
    float3  _S2300 = _S2293[int(6)];
    float3  _S2301 = _S2293[int(7)];
    float3  _S2302 = _S2293[int(8)];
    float3  _S2303 = _S2293[int(9)];
    float3  _S2304 = _S2293[int(10)];
    float3  _S2305 = _S2293[int(11)];
    float3  _S2306 = _S2293[int(12)];
    float3  _S2307 = _S2293[int(13)];
    float3  _S2308 = _S2293[int(14)];
    float3  _S2309 = _S2293[int(15)];
    float2  _S2310 = make_float2 (0.0f, _S2284.differential_0);
    float _S2311;
    if(antialiased_14)
    {
        float _S2312 = _S2180 * _S2292;
        _S2182 = _S2177 * _S2292;
        _S2311 = _S2312;
    }
    else
    {
        _S2182 = _S2292;
        _S2311 = 0.0f;
    }
    float _S2313 = - (_S2182 / _S2181);
    DiffPair_float_0 _S2314;
    (&_S2314)->primal_0 = _S2178;
    (&_S2314)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2314, _S2313);
    float _S2315 = - _S2314.differential_0;
    DiffPair_float_0 _S2316;
    (&_S2316)->primal_0 = _S2176;
    (&_S2316)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2316, _S2311);
    DiffPair_float_0 _S2317;
    (&_S2317)->primal_0 = 0.0f;
    (&_S2317)->differential_0 = 0.0f;
    DiffPair_float_0 _S2318;
    (&_S2318)->primal_0 = _S2174;
    (&_S2318)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2317, &_S2318, _S2316.differential_0);
    float _S2319 = _S2318.differential_0 / _S2175;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2319;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2319;
    float _S2320 = - s_diff_det_blur_T_1;
    float _S2321 = _S2169 * s_diff_det_blur_T_1;
    float _S2322 = _S2171 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2323 = _S2150;
    _S2323[int(1)] = _S2310;
    _S2323[int(0)] = _S2291;
    _S2170 = _S2323;
    *&(((&_S2170)->rows + (int(1)))->y) = 0.0f;
    float _S2324 = _S2321 + _S2323.rows[int(1)].y;
    *&(((&_S2170)->rows + (int(0)))->x) = 0.0f;
    float _S2325 = _S2322 + _S2323.rows[int(0)].x;
    float _S2326 = _S2320 + - s_diff_det_orig_T_4;
    float _S2327 = _S2167._S2148.rows[int(0)].y * _S2326;
    float _S2328 = _S2167._S2148.rows[int(1)].x * _S2326;
    float _S2329 = _S2167._S2148.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2330 = _S2324 + _S2167._S2148.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2331 = _S2165;
    *&((&_S2331)->x) = _S2327;
    *&((&_S2331)->y) = _S2330;
    float _S2332 = _S2325 + _S2329;
    float2  _S2333 = _S2165;
    *&((&_S2333)->y) = _S2328;
    *&((&_S2333)->x) = _S2332;
    Matrix<float, 2, 2>  _S2334 = _S2150;
    _S2334[int(1)] = _S2331;
    _S2334[int(0)] = _S2333;
    Matrix<float, 2, 2>  _S2335 = _S2170 + _S2334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2336;
    (&_S2336)->primal_0 = mean_c_14;
    (&_S2336)->differential_0 = _S2215;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2337;
    (&_S2337)->primal_0 = _S2163;
    (&_S2337)->differential_0 = _S2277;
    DiffPair_float_0 _S2338;
    (&_S2338)->primal_0 = fx_18;
    (&_S2338)->differential_0 = 0.0f;
    DiffPair_float_0 _S2339;
    (&_S2339)->primal_0 = fy_18;
    (&_S2339)->differential_0 = 0.0f;
    DiffPair_float_0 _S2340;
    (&_S2340)->primal_0 = cx_18;
    (&_S2340)->differential_0 = 0.0f;
    DiffPair_float_0 _S2341;
    (&_S2341)->primal_0 = cy_18;
    (&_S2341)->differential_0 = 0.0f;
    float4  _S2342 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2343;
    (&_S2343)->primal_0 = radial_coeffs_21;
    (&_S2343)->differential_0 = _S2342;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2344;
    (&_S2344)->primal_0 = tangential_coeffs_21;
    (&_S2344)->differential_0 = _S2165;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2345;
    (&_S2345)->primal_0 = thin_prism_coeffs_21;
    (&_S2345)->differential_0 = _S2165;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2346 = _S2167._S2149;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2336, &_S2337, &_S2338, &_S2339, &_S2340, &_S2341, &_S2343, &_S2344, &_S2345, _S2335, v_mean2d_4, &_S2346);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2347;
    (&_S2347)->primal_0 = _S2161;
    (&_S2347)->differential_0 = _S2277;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2348;
    (&_S2348)->primal_0 = _S2162;
    (&_S2348)->differential_0 = _S2277;
    s_bwd_prop_mul_4(&_S2347, &_S2348, _S2337.differential_0);
    Matrix<float, 3, 3>  _S2349 = transpose_0(_S2348.differential_0 + _S2280.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2350;
    (&_S2350)->primal_0 = R_18;
    (&_S2350)->differential_0 = _S2277;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2351;
    (&_S2351)->primal_0 = _S2160;
    (&_S2351)->differential_0 = _S2277;
    s_bwd_prop_mul_4(&_S2350, &_S2351, _S2347.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2352;
    (&_S2352)->primal_0 = _S2158;
    (&_S2352)->differential_0 = _S2277;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2353;
    (&_S2353)->primal_0 = _S2159;
    (&_S2353)->differential_0 = _S2277;
    s_bwd_prop_mul_4(&_S2352, &_S2353, _S2351.differential_0);
    Matrix<float, 3, 3>  _S2354 = _S2352.differential_0 + transpose_0(_S2353.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2355;
    (&_S2355)->primal_0 = _S2157;
    (&_S2355)->differential_0 = _S2277;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2356;
    (&_S2356)->primal_0 = S_4;
    (&_S2356)->differential_0 = _S2277;
    s_bwd_prop_mul_4(&_S2355, &_S2356, _S2354);
    Matrix<float, 3, 3>  _S2357 = transpose_0(_S2355.differential_0);
    float _S2358 = 2.0f * - _S2357.rows[int(2)].z;
    float _S2359 = 2.0f * _S2357.rows[int(2)].y;
    float _S2360 = 2.0f * _S2357.rows[int(2)].x;
    float _S2361 = 2.0f * _S2357.rows[int(1)].z;
    float _S2362 = 2.0f * - _S2357.rows[int(1)].y;
    float _S2363 = 2.0f * _S2357.rows[int(1)].x;
    float _S2364 = 2.0f * _S2357.rows[int(0)].z;
    float _S2365 = 2.0f * _S2357.rows[int(0)].y;
    float _S2366 = 2.0f * - _S2357.rows[int(0)].x;
    float _S2367 = - _S2363 + _S2365;
    float _S2368 = _S2360 + - _S2364;
    float _S2369 = - _S2359 + _S2361;
    float _S2370 = _S2359 + _S2361;
    float _S2371 = _S2360 + _S2364;
    float _S2372 = _S2363 + _S2365;
    float _S2373 = quat_17.w * (_S2362 + _S2366);
    float _S2374 = quat_17.z * (_S2358 + _S2366);
    float _S2375 = quat_17.y * (_S2358 + _S2362);
    float _S2376 = quat_17.x * _S2367 + quat_17.z * _S2370 + quat_17.y * _S2371 + _S2373 + _S2373;
    float _S2377 = quat_17.x * _S2368 + quat_17.w * _S2370 + quat_17.y * _S2372 + _S2374 + _S2374;
    float _S2378 = quat_17.x * _S2369 + quat_17.w * _S2371 + quat_17.z * _S2372 + _S2375 + _S2375;
    float _S2379 = quat_17.w * _S2367 + quat_17.z * _S2368 + quat_17.y * _S2369;
    float3  _S2380 = _S2215;
    *&((&_S2380)->z) = _S2356.differential_0.rows[int(2)].z;
    *&((&_S2380)->y) = _S2356.differential_0.rows[int(1)].y;
    *&((&_S2380)->x) = _S2356.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2381;
    (&_S2381)->primal_0 = scale_16;
    (&_S2381)->differential_0 = _S2215;
    s_bwd_prop_exp_1(&_S2381, _S2380);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2382;
    (&_S2382)->primal_0 = mean_c_14;
    (&_S2382)->differential_0 = _S2215;
    s_bwd_length_impl_0(&_S2382, v_depth_4);
    float3  _S2383 = _S2336.differential_0 + _S2382.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2384;
    (&_S2384)->primal_0 = R_18;
    (&_S2384)->differential_0 = _S2277;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2385;
    (&_S2385)->primal_0 = mean_14;
    (&_S2385)->differential_0 = _S2215;
    s_bwd_prop_mul_1(&_S2384, &_S2385, _S2383);
    float3  _S2386 = _S2383 + _S2281.differential_0;
    Matrix<float, 3, 3>  _S2387 = _S2349 + _S2350.differential_0 + _S2384.differential_0;
    float3  _S2388 = _S2381.differential_0 + _S2283;
    float4  _S2389 = _S2342;
    *&((&_S2389)->w) = _S2376;
    *&((&_S2389)->z) = _S2377;
    *&((&_S2389)->y) = _S2378;
    *&((&_S2389)->x) = _S2379;
    float4  _S2390 = _S2389;
    float3  _S2391 = _S2385.differential_0 + _S2275;
    *v_mean_4 = _S2391;
    *v_quat_4 = _S2390;
    *v_scale_4 = _S2388;
    *v_in_opacity_4 = _S2315;
    (*v_sh_coeffs_4)[int(0)] = _S2294;
    (*v_sh_coeffs_4)[int(1)] = _S2295;
    (*v_sh_coeffs_4)[int(2)] = _S2296;
    (*v_sh_coeffs_4)[int(3)] = _S2297;
    (*v_sh_coeffs_4)[int(4)] = _S2298;
    (*v_sh_coeffs_4)[int(5)] = _S2299;
    (*v_sh_coeffs_4)[int(6)] = _S2300;
    (*v_sh_coeffs_4)[int(7)] = _S2301;
    (*v_sh_coeffs_4)[int(8)] = _S2302;
    (*v_sh_coeffs_4)[int(9)] = _S2303;
    (*v_sh_coeffs_4)[int(10)] = _S2304;
    (*v_sh_coeffs_4)[int(11)] = _S2305;
    (*v_sh_coeffs_4)[int(12)] = _S2306;
    (*v_sh_coeffs_4)[int(13)] = _S2307;
    (*v_sh_coeffs_4)[int(14)] = _S2308;
    (*v_sh_coeffs_4)[int(15)] = _S2309;
    *v_R_5 = _S2387;
    *v_t_5 = _S2386;
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
    float _S2392 = (*dpquat_0).primal_0.y;
    float x2_19 = _S2392 * _S2392;
    float y2_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_34 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2393 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_34), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_34), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19))));
    Matrix<float, 3, 3>  _S2394 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2395;
    (&_S2395)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2395)->differential_0 = _S2394;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2396;
    (&_S2396)->primal_0 = _S2393;
    (&_S2396)->differential_0 = _S2394;
    s_bwd_prop_mul_4(&_S2395, &_S2396, _s_dOut_7);
    Matrix<float, 3, 3>  _S2397 = transpose_0(transpose_0(_S2396.differential_0));
    float _S2398 = 2.0f * - _S2397.rows[int(2)].z;
    float _S2399 = 2.0f * _S2397.rows[int(2)].y;
    float _S2400 = 2.0f * _S2397.rows[int(2)].x;
    float _S2401 = 2.0f * _S2397.rows[int(1)].z;
    float _S2402 = 2.0f * - _S2397.rows[int(1)].y;
    float _S2403 = 2.0f * _S2397.rows[int(1)].x;
    float _S2404 = 2.0f * _S2397.rows[int(0)].z;
    float _S2405 = 2.0f * _S2397.rows[int(0)].y;
    float _S2406 = 2.0f * - _S2397.rows[int(0)].x;
    float _S2407 = - _S2403 + _S2405;
    float _S2408 = _S2400 + - _S2404;
    float _S2409 = - _S2399 + _S2401;
    float _S2410 = _S2399 + _S2401;
    float _S2411 = _S2400 + _S2404;
    float _S2412 = _S2403 + _S2405;
    float _S2413 = (*dpquat_0).primal_0.w * (_S2402 + _S2406);
    float _S2414 = (*dpquat_0).primal_0.z * (_S2398 + _S2406);
    float _S2415 = (*dpquat_0).primal_0.y * (_S2398 + _S2402);
    float _S2416 = (*dpquat_0).primal_0.x * _S2407 + (*dpquat_0).primal_0.z * _S2410 + (*dpquat_0).primal_0.y * _S2411 + _S2413 + _S2413;
    float _S2417 = (*dpquat_0).primal_0.x * _S2408 + (*dpquat_0).primal_0.w * _S2410 + (*dpquat_0).primal_0.y * _S2412 + _S2414 + _S2414;
    float _S2418 = (*dpquat_0).primal_0.x * _S2409 + (*dpquat_0).primal_0.w * _S2411 + (*dpquat_0).primal_0.z * _S2412 + _S2415 + _S2415;
    float _S2419 = (*dpquat_0).primal_0.w * _S2407 + (*dpquat_0).primal_0.z * _S2408 + (*dpquat_0).primal_0.y * _S2409;
    float3  _S2420 = make_float3 (_S2395.differential_0.rows[int(0)].x, _S2395.differential_0.rows[int(1)].y, _S2395.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S2420;
    float4  _S2421 = make_float4 (0.0f);
    *&((&_S2421)->w) = _S2416;
    *&((&_S2421)->z) = _S2417;
    *&((&_S2421)->y) = _S2418;
    *&((&_S2421)->x) = _S2419;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S2421;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S2422, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2423, Matrix<float, 3, 3>  _S2424)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S2422, _S2423, _S2424);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_19, float3  scale_18, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S2425 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_19;
    (&dp_quat_0)->differential_0 = _S2425;
    float3  _S2426 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_18;
    (&dp_scale_0)->differential_0 = _S2426;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_15)
{
    float _S2427 = dOut_15.y;
    float _S2428 = dOut_15.z;
    float _S2429 = dOut_15.x;
    float _S2430 = (*a_0).primal_0.z * _S2427 + - (*a_0).primal_0.y * _S2428;
    float _S2431 = - (*a_0).primal_0.z * _S2429 + (*a_0).primal_0.x * _S2428;
    float _S2432 = (*a_0).primal_0.y * _S2429 + - (*a_0).primal_0.x * _S2427;
    float3  _S2433 = make_float3 (- (*b_0).primal_0.z * _S2427 + (*b_0).primal_0.y * _S2428, (*b_0).primal_0.z * _S2429 + - (*b_0).primal_0.x * _S2428, - (*b_0).primal_0.y * _S2429 + (*b_0).primal_0.x * _S2427);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S2433;
    float3  _S2434 = make_float3 (_S2430, _S2431, _S2432);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S2434;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S2435 = left_10.y;
    float _S2436 = right_10.z;
    float _S2437 = left_10.z;
    float _S2438 = right_10.y;
    float _S2439 = right_10.x;
    float _S2440 = left_10.x;
    return make_float3 (_S2435 * _S2436 - _S2437 * _S2438, _S2437 * _S2439 - _S2440 * _S2436, _S2440 * _S2438 - _S2435 * _S2439);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_15, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_1, float3  ray_d_1)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_1);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_1 - mean_15));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S2441, float3  _S2442)
{
    return cross_0(_S2441, _S2442);
}

inline __device__ float s_primal_ctx_dot_0(float3  _S2443, float3  _S2444)
{
    return dot_0(_S2443, _S2444);
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2445, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2446, float _S2447)
{
    _d_dot_0(_S2445, _S2446, _S2447);
    return;
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2448, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2449, float3  _S2450)
{
    _d_cross_0(_S2448, _S2449, _S2450);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_8)
{
    float3  _S2451 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S2452 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, _S2451);
    float3  _S2453 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, (*dpray_d_2).primal_0);
    float3  _S2454 = s_primal_ctx_cross_0(_S2453, _S2452);
    float _S2455 = -0.5f * s_primal_ctx_dot_0(_S2454, _S2454);
    float _S2456 = s_primal_ctx_dot_0(_S2453, _S2453);
    float _S2457 = _S2455 / _S2456;
    float _S2458 = _S2456 * _S2456;
    float _S2459 = (*dpopacity_0).primal_0 * _s_dOut_8;
    float _S2460 = s_primal_ctx_exp_1(_S2457) * _s_dOut_8;
    DiffPair_float_0 _S2461;
    (&_S2461)->primal_0 = _S2457;
    (&_S2461)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2461, _S2459);
    float _S2462 = _S2461.differential_0 / _S2458;
    float _S2463 = _S2455 * - _S2462;
    float _S2464 = _S2456 * _S2462;
    float3  _S2465 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2466;
    (&_S2466)->primal_0 = _S2453;
    (&_S2466)->differential_0 = _S2465;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2467;
    (&_S2467)->primal_0 = _S2453;
    (&_S2467)->differential_0 = _S2465;
    s_bwd_prop_dot_0(&_S2466, &_S2467, _S2463);
    float _S2468 = -0.5f * _S2464;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2469;
    (&_S2469)->primal_0 = _S2454;
    (&_S2469)->differential_0 = _S2465;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2470;
    (&_S2470)->primal_0 = _S2454;
    (&_S2470)->differential_0 = _S2465;
    s_bwd_prop_dot_0(&_S2469, &_S2470, _S2468);
    float3  _S2471 = _S2470.differential_0 + _S2469.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2472;
    (&_S2472)->primal_0 = _S2453;
    (&_S2472)->differential_0 = _S2465;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2473;
    (&_S2473)->primal_0 = _S2452;
    (&_S2473)->differential_0 = _S2465;
    s_bwd_prop_cross_0(&_S2472, &_S2473, _S2471);
    float3  _S2474 = _S2467.differential_0 + _S2466.differential_0 + _S2472.differential_0;
    Matrix<float, 3, 3>  _S2475 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2476;
    (&_S2476)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2476)->differential_0 = _S2475;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2477;
    (&_S2477)->primal_0 = (*dpray_d_2).primal_0;
    (&_S2477)->differential_0 = _S2465;
    s_bwd_prop_mul_1(&_S2476, &_S2477, _S2474);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2478;
    (&_S2478)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2478)->differential_0 = _S2475;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2479;
    (&_S2479)->primal_0 = _S2451;
    (&_S2479)->differential_0 = _S2465;
    s_bwd_prop_mul_1(&_S2478, &_S2479, _S2473.differential_0);
    float3  _S2480 = - _S2479.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S2477.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S2479.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S2460;
    Matrix<float, 3, 3>  _S2481 = _S2476.differential_0 + _S2478.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S2481;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S2480;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2482, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S2483, DiffPair_float_0 * _S2484, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2485, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2486, float _S2487)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S2482, _S2483, _S2484, _S2485, _S2486, _S2487);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S2488 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_16;
    (&dp_mean_0)->differential_0 = _S2488;
    Matrix<float, 3, 3>  _S2489 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S2489;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S2488;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S2488;
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
    *depth_10 = length_1(mean_17 - ray_o_3);
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S2490 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2491;
    (&_S2491)->primal_0 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    (&_S2491)->differential_0 = _S2490;
    s_bwd_length_impl_0(&_S2491, dpdepth_0);
    float3  _S2492 = - _S2491.differential_0;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S2490;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S2492;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S2493 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S2493;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S2491.differential_0;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2494, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S2495, DiffPair_float_0 * _S2496, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2497, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2498, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2499, float3  _S2500, float _S2501)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S2494, _S2495, _S2496, _S2497, _S2498, _S2499, _S2500, _S2501);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S2502 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_18;
    (&dp_mean_1)->differential_0 = _S2502;
    Matrix<float, 3, 3>  _S2503 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S2503;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S2502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S2502;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S2502;
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
    float _S2504 = _slang_select(((*dpx_13).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_13).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_16;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S2504;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_14, float dOut_17)
{
    float _S2505 = (F32_exp2(((*dpx_14).primal_0))) * 50.693145751953125f * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S2505;
    return;
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_19, float4  quat_20, float3  scale_19, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_19, float3  t_18, float fx_19, float fy_19, float cx_19, float cy_19, float4  radial_coeffs_22, float2  tangential_coeffs_22, float2  thin_prism_coeffs_22, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_19, mean_19) + t_18;
        float _S2506 = mean_c_15.z;
        bool _S2507;
        if(_S2506 < near_plane_10)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = _S2506 > far_plane_10;
        }
        if(_S2507)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S2508 = scale_19.x;
        float sx_0 = (F32_exp((_S2508)));
        float _S2509 = scale_19.y;
        float sy_0 = (F32_exp((_S2509)));
        float sz_0 = scale_19.z - 0.5f * (_S2508 + _S2509);
        float x_40 = quat_20.y;
        float x2_20 = x_40 * x_40;
        float y2_20 = quat_20.z * quat_20.z;
        float z2_35 = quat_20.w * quat_20.w;
        float xy_20 = quat_20.y * quat_20.z;
        float xz_20 = quat_20.y * quat_20.w;
        float yz_20 = quat_20.z * quat_20.w;
        float wx_20 = quat_20.x * quat_20.y;
        float wy_20 = quat_20.x * quat_20.z;
        float wz_20 = quat_20.x * quat_20.w;
        Matrix<float, 3, 3>  _S2510 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20)));
        float3  vert0_c_0 = mul_0(R_19, mul_0(_S2510, make_float3 (sx_0, 0.0f, 0.0f)) + mean_19) + t_18;
        float3  vert1_c_0 = mul_0(R_19, mul_0(_S2510, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_19) + t_18;
        float3  vert2_c_0 = mul_0(R_19, mul_0(_S2510, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_19) + t_18;
        float _S2511 = vert0_c_0.z;
        float _S2512 = vert1_c_0.z;
        float _S2513 = vert2_c_0.z;
        if(_S2511 < near_plane_10)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = _S2511 > far_plane_10;
        }
        if(_S2507)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = _S2512 < near_plane_10;
        }
        if(_S2507)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = _S2512 > far_plane_10;
        }
        if(_S2507)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = _S2513 < near_plane_10;
        }
        if(_S2507)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = _S2513 > far_plane_10;
        }
        if(_S2507)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S2511);
        *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S2512);
        *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S2513);
        float2  _S2514 = make_float2 (fx_19, fy_19);
        float2  _S2515 = make_float2 (cx_19, cy_19);
        *uv0_0 = _S2514 * *uv0_0 + _S2515;
        *uv1_0 = _S2514 * *uv1_0 + _S2515;
        float2  _S2516 = _S2514 * *uv2_0 + _S2515;
        *uv2_0 = _S2516;
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = _S2516 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - _S2516)));
        float _S2517 = _S2516.x;
        float xmax_5 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), (_S2517))) + offset_0;
        float xmin_5 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), (_S2517))) - offset_0;
        float _S2518 = _S2516.y;
        float ymax_5 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), (_S2518))) + offset_0;
        float ymin_5 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), (_S2518))) - offset_0;
        if(xmax_5 <= 0.0f)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = xmin_5 >= float(image_width_15);
        }
        if(_S2507)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = ymax_5 <= 0.0f;
        }
        if(_S2507)
        {
            _S2507 = true;
        }
        else
        {
            _S2507 = ymin_5 >= float(image_height_15);
        }
        if(_S2507)
        {
            _S2507 = true;
        }
        else
        {
            if(_S2506 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S2507 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S2507 = false;
                }
                if(_S2507)
                {
                    _S2507 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S2507 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S2507 = false;
                    }
                }
            }
            else
            {
                _S2507 = false;
            }
        }
        if(_S2507)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0));
        *out_hardness_0 = hardness_0;
        float3  _S2519 = mean_19 - - mul_0(transpose_0(R_19), t_18);
        float _S2520 = _S2519.x;
        float _S2521 = _S2519.y;
        float _S2522 = _S2519.z;
        float norm_10 = (F32_sqrt((_S2520 * _S2520 + _S2521 * _S2521 + _S2522 * _S2522)));
        float x_41 = _S2520 / norm_10;
        float y_18 = _S2521 / norm_10;
        float z_15 = _S2522 / norm_10;
        float z2_36 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_41 * x_41 - y_18 * y_18;
        float fS1_15 = 2.0f * x_41 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_36 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_41) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_36 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_41) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_41 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_36 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_41) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_41 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S2523 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S2523);
        float3  _S2524 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S2525 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S2524 + _S2525 + make_float3 (0.5f), _S2523);
        (*rgb_12)[int(2)] = max_0(_S2524 - _S2525 + make_float3 (0.5f), _S2523);
        float3  _S2526 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S2526 * make_float3 (float(- (F32_sign((dot_0(_S2526, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_20, float3  t_19, float fx_20, float fy_20, float cx_20, float cy_20, float4  radial_coeffs_23, float2  tangential_coeffs_23, float2  thin_prism_coeffs_23, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_20, mean_20) + t_19;
        float _S2527 = length_1(mean_c_16);
        bool _S2528;
        if(_S2527 < near_plane_11)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = _S2527 > far_plane_11;
        }
        if(_S2528)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S2529 = scale_20.x;
        float sx_1 = (F32_exp((_S2529)));
        float _S2530 = scale_20.y;
        float sy_1 = (F32_exp((_S2530)));
        float sz_1 = scale_20.z - 0.5f * (_S2529 + _S2530);
        float x_42 = quat_21.y;
        float x2_21 = x_42 * x_42;
        float y2_21 = quat_21.z * quat_21.z;
        float z2_37 = quat_21.w * quat_21.w;
        float xy_21 = quat_21.y * quat_21.z;
        float xz_21 = quat_21.y * quat_21.w;
        float yz_21 = quat_21.z * quat_21.w;
        float wx_21 = quat_21.x * quat_21.y;
        float wy_21 = quat_21.x * quat_21.z;
        float wz_21 = quat_21.x * quat_21.w;
        Matrix<float, 3, 3>  _S2531 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_37), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_37), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_1 = mul_0(R_20, mul_0(_S2531, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_19;
        float3  vert1_c_1 = mul_0(R_20, mul_0(_S2531, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_19;
        float3  vert2_c_1 = mul_0(R_20, mul_0(_S2531, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_19;
        float _S2532 = length_1(vert0_c_1);
        float _S2533 = length_1(vert1_c_1);
        float _S2534 = length_1(vert2_c_1);
        if(_S2532 < near_plane_11)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = _S2532 > far_plane_11;
        }
        if(_S2528)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = _S2533 < near_plane_11;
        }
        if(_S2528)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = _S2533 > far_plane_11;
        }
        if(_S2528)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = _S2534 < near_plane_11;
        }
        if(_S2528)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = _S2534 > far_plane_11;
        }
        if(_S2528)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_2 = CameraDistortion_x24init_0(radial_coeffs_23, tangential_coeffs_23, thin_prism_coeffs_23);
        float2  _S2535 = float2 {vert0_c_1.x, vert0_c_1.y};
        float r_7 = length_0(_S2535);
        float _S2536 = vert0_c_1.z;
        float theta_2 = (F32_atan2((r_7), (_S2536)));
        float k_4;
        if(theta_2 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_2 * theta_2 / 3.0f) / _S2536;
        }
        else
        {
            k_4 = theta_2 / r_7;
        }
        float2  _S2537 = _S2535 * make_float2 (k_4);
        float k1_3 = dist_coeffs_2.radial_coeffs_0.x;
        float k2_3 = dist_coeffs_2.radial_coeffs_0.y;
        float k3_3 = dist_coeffs_2.radial_coeffs_0.z;
        float k4_4 = dist_coeffs_2.radial_coeffs_0.w;
        float p1_4 = dist_coeffs_2.tangential_coeffs_0.x;
        float p2_4 = dist_coeffs_2.tangential_coeffs_0.y;
        float sx1_4 = dist_coeffs_2.thin_prism_coeffs_0.x;
        float sy1_4 = dist_coeffs_2.thin_prism_coeffs_0.y;
        float u_10 = _S2537.x;
        float v_10 = _S2537.y;
        float r2_10 = u_10 * u_10 + v_10 * v_10;
        float _S2538 = 2.0f * p1_4;
        float _S2539 = 2.0f * p2_4;
        float2  _S2540 = _S2537 * make_float2 (1.0f + r2_10 * (k1_3 + r2_10 * (k2_3 + r2_10 * (k3_3 + r2_10 * k4_4)))) + make_float2 (_S2538 * u_10 * v_10 + p2_4 * (r2_10 + 2.0f * u_10 * u_10) + sx1_4 * r2_10, _S2539 * u_10 * v_10 + p1_4 * (r2_10 + 2.0f * v_10 * v_10) + sy1_4 * r2_10);
        *uv0_1 = make_float2 (fx_20 * _S2540.x + cx_20, fy_20 * _S2540.y + cy_20);
        float2  _S2541 = float2 {vert1_c_1.x, vert1_c_1.y};
        float r_8 = length_0(_S2541);
        float _S2542 = vert1_c_1.z;
        float theta_3 = (F32_atan2((r_8), (_S2542)));
        if(theta_3 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_3 * theta_3 / 3.0f) / _S2542;
        }
        else
        {
            k_4 = theta_3 / r_8;
        }
        float2  _S2543 = _S2541 * make_float2 (k_4);
        float u_11 = _S2543.x;
        float v_11 = _S2543.y;
        float r2_11 = u_11 * u_11 + v_11 * v_11;
        float2  _S2544 = _S2543 * make_float2 (1.0f + r2_11 * (k1_3 + r2_11 * (k2_3 + r2_11 * (k3_3 + r2_11 * k4_4)))) + make_float2 (_S2538 * u_11 * v_11 + p2_4 * (r2_11 + 2.0f * u_11 * u_11) + sx1_4 * r2_11, _S2539 * u_11 * v_11 + p1_4 * (r2_11 + 2.0f * v_11 * v_11) + sy1_4 * r2_11);
        *uv1_1 = make_float2 (fx_20 * _S2544.x + cx_20, fy_20 * _S2544.y + cy_20);
        float2  _S2545 = float2 {vert2_c_1.x, vert2_c_1.y};
        float r_9 = length_0(_S2545);
        float _S2546 = vert2_c_1.z;
        float theta_4 = (F32_atan2((r_9), (_S2546)));
        if(theta_4 < 0.00100000004749745f)
        {
            k_4 = (1.0f - theta_4 * theta_4 / 3.0f) / _S2546;
        }
        else
        {
            k_4 = theta_4 / r_9;
        }
        float2  _S2547 = _S2545 * make_float2 (k_4);
        float u_12 = _S2547.x;
        float v_12 = _S2547.y;
        float r2_12 = u_12 * u_12 + v_12 * v_12;
        float2  _S2548 = _S2547 * make_float2 (1.0f + r2_12 * (k1_3 + r2_12 * (k2_3 + r2_12 * (k3_3 + r2_12 * k4_4)))) + make_float2 (_S2538 * u_12 * v_12 + p2_4 * (r2_12 + 2.0f * u_12 * u_12) + sx1_4 * r2_12, _S2539 * u_12 * v_12 + p1_4 * (r2_12 + 2.0f * v_12 * v_12) + sy1_4 * r2_12);
        float _S2549 = fx_20 * _S2548.x + cx_20;
        float _S2550 = fy_20 * _S2548.y + cy_20;
        float2  _S2551 = make_float2 (_S2549, _S2550);
        *uv2_1 = _S2551;
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = _S2551 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - _S2551)));
        float xmax_6 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), (_S2549))) + offset_1;
        float xmin_6 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), (_S2549))) - offset_1;
        float ymax_6 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), (_S2550))) + offset_1;
        float ymin_6 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), (_S2550))) - offset_1;
        if(xmax_6 <= 0.0f)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = xmin_6 >= float(image_width_16);
        }
        if(_S2528)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = ymax_6 <= 0.0f;
        }
        if(_S2528)
        {
            _S2528 = true;
        }
        else
        {
            _S2528 = ymin_6 >= float(image_height_16);
        }
        if(_S2528)
        {
            _S2528 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S2528 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S2528 = false;
                }
                if(_S2528)
                {
                    _S2528 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S2528 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S2528 = false;
                    }
                }
            }
            else
            {
                _S2528 = false;
            }
        }
        if(_S2528)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = make_float3 (_S2532, _S2533, _S2534);
        *out_hardness_1 = hardness_1;
        float3  _S2552 = mean_20 - - mul_0(transpose_0(R_20), t_19);
        float _S2553 = _S2552.x;
        float _S2554 = _S2552.y;
        float _S2555 = _S2552.z;
        float norm_11 = (F32_sqrt((_S2553 * _S2553 + _S2554 * _S2554 + _S2555 * _S2555)));
        float x_43 = _S2553 / norm_11;
        float y_19 = _S2554 / norm_11;
        float z_16 = _S2555 / norm_11;
        float z2_38 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_43 * x_43 - y_19 * y_19;
        float fS1_16 = 2.0f * x_43 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_38 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_43) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_38 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_43) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_43 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_38 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_43) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_43 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S2556 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S2556);
        float3  _S2557 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S2558 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S2557 + _S2558 + make_float3 (0.5f), _S2556);
        (*rgb_13)[int(2)] = max_0(_S2557 - _S2558 + make_float3 (0.5f), _S2556);
        float3  _S2559 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S2559 * make_float3 (float(- (F32_sign((dot_0(_S2559, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_21, float3  t_20, float fx_21, float fy_21, float cx_21, float cy_21, float4  radial_coeffs_24, float2  tangential_coeffs_24, float2  thin_prism_coeffs_24, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_21, mean_21) + t_20;
    float _S2560 = scale_21.x;
    float sx_2 = (F32_exp((_S2560)));
    float _S2561 = scale_21.y;
    float sy_2 = (F32_exp((_S2561)));
    float sz_2 = scale_21.z - 0.5f * (_S2560 + _S2561);
    float x_44 = quat_22.y;
    float x2_22 = x_44 * x_44;
    float y2_22 = quat_22.z * quat_22.z;
    float z2_39 = quat_22.w * quat_22.w;
    float xy_22 = quat_22.y * quat_22.z;
    float xz_22 = quat_22.y * quat_22.w;
    float yz_22 = quat_22.z * quat_22.w;
    float wx_22 = quat_22.x * quat_22.y;
    float wy_22 = quat_22.x * quat_22.z;
    float wz_22 = quat_22.x * quat_22.w;
    Matrix<float, 3, 3>  _S2562 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_39), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_39), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
    float3  vert0_c_2 = mul_0(R_21, mul_0(_S2562, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_20;
    float3  vert1_c_2 = mul_0(R_21, mul_0(_S2562, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_20;
    float3  vert2_c_2 = mul_0(R_21, mul_0(_S2562, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_20;
    *uv0_2 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    *uv1_2 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    *uv2_2 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float2  _S2563 = make_float2 (fx_21, fy_21);
    float2  _S2564 = make_float2 (cx_21, cy_21);
    *uv0_2 = _S2563 * *uv0_2 + _S2564;
    *uv1_2 = _S2563 * *uv1_2 + _S2564;
    float2  _S2565 = _S2563 * *uv2_2 + _S2564;
    *uv2_2 = _S2565;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S2565 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S2565)));
    float _S2566 = _S2565.x;
    float _S2567 = _S2565.y;
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S2566))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S2567))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S2566))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S2567))) + offset_2)))));
    *depth_13 = make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2));
    *out_hardness_2 = hardness_2;
    float3  _S2568 = mean_21 - - mul_0(transpose_0(R_21), t_20);
    float _S2569 = _S2568.x;
    float _S2570 = _S2568.y;
    float _S2571 = _S2568.z;
    float norm_12 = (F32_sqrt((_S2569 * _S2569 + _S2570 * _S2570 + _S2571 * _S2571)));
    float x_45 = _S2569 / norm_12;
    float y_20 = _S2570 / norm_12;
    float z_17 = _S2571 / norm_12;
    float z2_40 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_45 * x_45 - y_20 * y_20;
    float fS1_17 = 2.0f * x_45 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_40 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_45) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_40 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_45) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_45 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_40 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_45) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_45 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S2572 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S2572);
    float3  _S2573 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S2574 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S2573 + _S2574 + make_float3 (0.5f), _S2572);
    (*rgb_14)[int(2)] = max_0(_S2573 - _S2574 + make_float3 (0.5f), _S2572);
    float3  _S2575 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S2575 * make_float3 (float(- (F32_sign((dot_0(_S2575, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_22, float3  t_21, float fx_22, float fy_22, float cx_22, float cy_22, float4  radial_coeffs_25, float2  tangential_coeffs_25, float2  thin_prism_coeffs_25, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_22, mean_22) + t_21;
    float _S2576 = scale_22.x;
    float sx_3 = (F32_exp((_S2576)));
    float _S2577 = scale_22.y;
    float sy_3 = (F32_exp((_S2577)));
    float sz_3 = scale_22.z - 0.5f * (_S2576 + _S2577);
    float x_46 = quat_23.y;
    float x2_23 = x_46 * x_46;
    float y2_23 = quat_23.z * quat_23.z;
    float z2_41 = quat_23.w * quat_23.w;
    float xy_23 = quat_23.y * quat_23.z;
    float xz_23 = quat_23.y * quat_23.w;
    float yz_23 = quat_23.z * quat_23.w;
    float wx_23 = quat_23.x * quat_23.y;
    float wy_23 = quat_23.x * quat_23.z;
    float wz_23 = quat_23.x * quat_23.w;
    Matrix<float, 3, 3>  _S2578 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_41), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_41), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_3 = mul_0(R_22, mul_0(_S2578, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_21;
    float3  vert1_c_3 = mul_0(R_22, mul_0(_S2578, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_21;
    float3  vert2_c_3 = mul_0(R_22, mul_0(_S2578, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_21;
    CameraDistortion_0 dist_coeffs_3 = CameraDistortion_x24init_0(radial_coeffs_25, tangential_coeffs_25, thin_prism_coeffs_25);
    float2  _S2579 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_10 = length_0(_S2579);
    float _S2580 = vert0_c_3.z;
    float theta_5 = (F32_atan2((r_10), (_S2580)));
    float k_5;
    if(theta_5 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_5 * theta_5 / 3.0f) / _S2580;
    }
    else
    {
        k_5 = theta_5 / r_10;
    }
    float2  _S2581 = _S2579 * make_float2 (k_5);
    float k1_4 = dist_coeffs_3.radial_coeffs_0.x;
    float k2_4 = dist_coeffs_3.radial_coeffs_0.y;
    float k3_4 = dist_coeffs_3.radial_coeffs_0.z;
    float k4_5 = dist_coeffs_3.radial_coeffs_0.w;
    float p1_5 = dist_coeffs_3.tangential_coeffs_0.x;
    float p2_5 = dist_coeffs_3.tangential_coeffs_0.y;
    float sx1_5 = dist_coeffs_3.thin_prism_coeffs_0.x;
    float sy1_5 = dist_coeffs_3.thin_prism_coeffs_0.y;
    float u_13 = _S2581.x;
    float v_13 = _S2581.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float _S2582 = 2.0f * p1_5;
    float _S2583 = 2.0f * p2_5;
    float2  _S2584 = _S2581 * make_float2 (1.0f + r2_13 * (k1_4 + r2_13 * (k2_4 + r2_13 * (k3_4 + r2_13 * k4_5)))) + make_float2 (_S2582 * u_13 * v_13 + p2_5 * (r2_13 + 2.0f * u_13 * u_13) + sx1_5 * r2_13, _S2583 * u_13 * v_13 + p1_5 * (r2_13 + 2.0f * v_13 * v_13) + sy1_5 * r2_13);
    *uv0_3 = make_float2 (fx_22 * _S2584.x + cx_22, fy_22 * _S2584.y + cy_22);
    float2  _S2585 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_11 = length_0(_S2585);
    float _S2586 = vert1_c_3.z;
    float theta_6 = (F32_atan2((r_11), (_S2586)));
    if(theta_6 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_6 * theta_6 / 3.0f) / _S2586;
    }
    else
    {
        k_5 = theta_6 / r_11;
    }
    float2  _S2587 = _S2585 * make_float2 (k_5);
    float u_14 = _S2587.x;
    float v_14 = _S2587.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S2588 = _S2587 * make_float2 (1.0f + r2_14 * (k1_4 + r2_14 * (k2_4 + r2_14 * (k3_4 + r2_14 * k4_5)))) + make_float2 (_S2582 * u_14 * v_14 + p2_5 * (r2_14 + 2.0f * u_14 * u_14) + sx1_5 * r2_14, _S2583 * u_14 * v_14 + p1_5 * (r2_14 + 2.0f * v_14 * v_14) + sy1_5 * r2_14);
    *uv1_3 = make_float2 (fx_22 * _S2588.x + cx_22, fy_22 * _S2588.y + cy_22);
    float2  _S2589 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_12 = length_0(_S2589);
    float _S2590 = vert2_c_3.z;
    float theta_7 = (F32_atan2((r_12), (_S2590)));
    if(theta_7 < 0.00100000004749745f)
    {
        k_5 = (1.0f - theta_7 * theta_7 / 3.0f) / _S2590;
    }
    else
    {
        k_5 = theta_7 / r_12;
    }
    float2  _S2591 = _S2589 * make_float2 (k_5);
    float u_15 = _S2591.x;
    float v_15 = _S2591.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S2592 = _S2591 * make_float2 (1.0f + r2_15 * (k1_4 + r2_15 * (k2_4 + r2_15 * (k3_4 + r2_15 * k4_5)))) + make_float2 (_S2582 * u_15 * v_15 + p2_5 * (r2_15 + 2.0f * u_15 * u_15) + sx1_5 * r2_15, _S2583 * u_15 * v_15 + p1_5 * (r2_15 + 2.0f * v_15 * v_15) + sy1_5 * r2_15);
    float _S2593 = fx_22 * _S2592.x + cx_22;
    float _S2594 = fy_22 * _S2592.y + cy_22;
    float2  _S2595 = make_float2 (_S2593, _S2594);
    *uv2_3 = _S2595;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S2595 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S2595)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S2593))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S2594))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S2593))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S2594))) + offset_3)))));
    *depth_14 = make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3));
    *out_hardness_3 = hardness_3;
    float3  _S2596 = mean_22 - - mul_0(transpose_0(R_22), t_21);
    float _S2597 = _S2596.x;
    float _S2598 = _S2596.y;
    float _S2599 = _S2596.z;
    float norm_13 = (F32_sqrt((_S2597 * _S2597 + _S2598 * _S2598 + _S2599 * _S2599)));
    float x_47 = _S2597 / norm_13;
    float y_21 = _S2598 / norm_13;
    float z_18 = _S2599 / norm_13;
    float z2_42 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_47 * x_47 - y_21 * y_21;
    float fS1_18 = 2.0f * x_47 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_42 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_47) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_42 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_47) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_47 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_42 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_47) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_47 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S2600 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S2600);
    float3  _S2601 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S2602 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S2601 + _S2602 + make_float3 (0.5f), _S2600);
    (*rgb_15)[int(2)] = max_0(_S2601 - _S2602 + make_float3 (0.5f), _S2600);
    float3  _S2603 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S2603 * make_float3 (float(- (F32_sign((dot_0(_S2603, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S2604, float _S2605)
{
    _d_exp2_0(_S2604, _S2605);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S2606, float _S2607)
{
    _d_abs_0(_S2606, _S2607);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_23, float3  t_22, float fx_23, float fy_23, float cx_23, float cy_23, float4  radial_coeffs_26, float2  tangential_coeffs_26, float2  thin_prism_coeffs_26, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_1(R_23, mean_23) + t_22;
    float _S2608 = scale_23.x;
    float _S2609 = s_primal_ctx_exp_1(_S2608);
    float _S2610 = scale_23.y;
    float _S2611 = s_primal_ctx_exp_1(_S2610);
    float sz_4 = scale_23.z - 0.5f * (_S2608 + _S2610);
    float _S2612 = quat_24.y;
    float x2_24 = _S2612 * _S2612;
    float y2_24 = quat_24.z * quat_24.z;
    float z2_43 = quat_24.w * quat_24.w;
    float xy_24 = quat_24.y * quat_24.z;
    float xz_24 = quat_24.y * quat_24.w;
    float yz_24 = quat_24.z * quat_24.w;
    float wx_24 = quat_24.x * quat_24.y;
    float wy_24 = quat_24.x * quat_24.z;
    float wz_24 = quat_24.x * quat_24.w;
    Matrix<float, 3, 3>  _S2613 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_43), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_43), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  _S2614 = make_float3 (_S2609, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_1(_S2613, _S2614) + mean_23;
    float _S2615 = -0.5f + sz_4;
    float3  _S2616 = make_float3 (_S2609 * _S2615, _S2611, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_1(_S2613, _S2616) + mean_23;
    float _S2617 = -0.5f - sz_4;
    float3  _S2618 = make_float3 (_S2609 * _S2617, - _S2611, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_1(_S2613, _S2618) + mean_23;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_23, vert0_0) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_23, vert1_0) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_23, vert2_0) + t_22;
    float2  _S2619 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S2620 = vert0_c_4.z;
    float2  _S2621 = make_float2 (_S2620);
    float2  _S2622 = make_float2 (_S2620 * _S2620);
    float2  _S2623 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S2624 = vert1_c_4.z;
    float2  _S2625 = make_float2 (_S2624);
    float2  _S2626 = make_float2 (_S2624 * _S2624);
    float2  _S2627 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S2628 = vert2_c_4.z;
    float2  _S2629 = make_float2 (_S2628);
    float2  _S2630 = make_float2 (_S2628 * _S2628);
    float2  _S2631 = make_float2 (fx_23, fy_23);
    float2  _S2632 = make_float2 (cx_23, cy_23);
    float2  _S2633 = _S2631 * (_S2619 / make_float2 (_S2620)) + _S2632;
    float2  _S2634 = _S2631 * (_S2623 / make_float2 (_S2624)) + _S2632;
    float2  _S2635 = _S2631 * (_S2627 / make_float2 (_S2628)) + _S2632;
    float2  e0_4 = _S2634 - _S2633;
    float2  e1_4 = _S2635 - _S2634;
    float2  e2_0 = _S2633 - _S2635;
    float _S2636 = e0_4.x;
    float _S2637 = e1_4.y;
    float _S2638 = e0_4.y;
    float _S2639 = e1_4.x;
    float _S2640 = _S2636 * _S2637 - _S2638 * _S2639;
    float _S2641 = 1.0f - hardness_4.y;
    float _S2642 = -1.0f / _S2641;
    float _S2643 = _S2641 * _S2641;
    float _S2644 = _S2633.x;
    float _S2645 = _S2634.x;
    float _S2646 = s_primal_ctx_max_0(_S2644, _S2645);
    float _S2647 = _S2635.x;
    float _S2648 = s_primal_ctx_min_0(_S2644, _S2645);
    float _S2649 = _S2633.y;
    float _S2650 = _S2634.y;
    float _S2651 = s_primal_ctx_max_0(_S2649, _S2650);
    float _S2652 = _S2635.y;
    float _S2653 = s_primal_ctx_min_0(_S2649, _S2650);
    Matrix<float, 3, 3>  _S2654 = transpose_0(R_23);
    float3  _S2655 = mean_23 - - s_primal_ctx_mul_1(_S2654, t_22);
    float _S2656 = _S2655.x;
    float _S2657 = _S2655.y;
    float _S2658 = _S2655.z;
    float _S2659 = _S2656 * _S2656 + _S2657 * _S2657 + _S2658 * _S2658;
    float _S2660 = s_primal_ctx_sqrt_0(_S2659);
    float x_48 = _S2656 / _S2660;
    float3  _S2661 = make_float3 (x_48);
    float _S2662 = _S2660 * _S2660;
    float y_22 = _S2657 / _S2660;
    float z_19 = _S2658 / _S2660;
    float3  _S2663 = make_float3 (z_19);
    float _S2664 = - y_22;
    float3  _S2665 = make_float3 (_S2664);
    float z2_44 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_48 * x_48 - y_22 * y_22;
    float _S2666 = 2.0f * x_48;
    float fS1_19 = _S2666 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_44 - 0.31539157032966614f;
    float3  _S2667 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_48;
    float3  _S2668 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S2669 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S2670 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S2671 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_44 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S2672 = 1.86588168144226074f * z2_44 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S2672;
    float3  _S2673 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_48;
    float3  _S2674 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S2675 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S2676 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S2677 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_48 * fC1_19 - y_22 * fS1_19);
    float3  _S2678 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_48 * fS1_19 + y_22 * fC1_19);
    float3  _S2679 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2664) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_48) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S2680 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S2681 = make_float3 (0.0f);
    float3  _S2682 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S2683 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S2684 = make_float3 (_S2683);
    float3  _S2685 = (*ch_coeffs_4)[int(1)] * make_float3 (_S2683);
    float3  _S2686 = _S2682 + _S2685 + make_float3 (0.5f);
    float3  _S2687 = _S2682 - _S2685 + make_float3 (0.5f);
    float3  _S2688 = vert1_c_4 - vert0_c_4;
    float3  _S2689 = vert2_c_4 - vert0_c_4;
    float3  _S2690 = s_primal_ctx_cross_0(_S2688, _S2689);
    float3  _S2691 = normalize_0(_S2690);
    float3  _S2692 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S2691, mean_c_19)))))) * v_normal_0;
    float3  _S2693 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2694;
    (&_S2694)->primal_0 = _S2691;
    (&_S2694)->differential_0 = _S2693;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2695;
    (&_S2695)->primal_0 = mean_c_19;
    (&_S2695)->differential_0 = _S2693;
    s_bwd_prop_dot_0(&_S2694, &_S2695, 0.0f);
    float3  _S2696 = _S2692 + _S2694.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2697;
    (&_S2697)->primal_0 = _S2690;
    (&_S2697)->differential_0 = _S2693;
    s_bwd_normalize_impl_0(&_S2697, _S2696);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2698;
    (&_S2698)->primal_0 = _S2688;
    (&_S2698)->differential_0 = _S2693;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2699;
    (&_S2699)->primal_0 = _S2689;
    (&_S2699)->differential_0 = _S2693;
    s_bwd_prop_cross_0(&_S2698, &_S2699, _S2697.differential_0);
    float3  _S2700 = - _S2699.differential_0;
    float3  _S2701 = - _S2698.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2702;
    (&_S2702)->primal_0 = _S2687;
    (&_S2702)->differential_0 = _S2693;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2703;
    (&_S2703)->primal_0 = _S2681;
    (&_S2703)->differential_0 = _S2693;
    s_bwd_prop_max_0(&_S2702, &_S2703, (*v_rgb_6)[int(2)]);
    float3  _S2704 = - _S2702.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2705;
    (&_S2705)->primal_0 = _S2686;
    (&_S2705)->differential_0 = _S2693;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2706;
    (&_S2706)->primal_0 = _S2681;
    (&_S2706)->differential_0 = _S2693;
    s_bwd_prop_max_0(&_S2705, &_S2706, (*v_rgb_6)[int(1)]);
    float3  _S2707 = _S2684 * (_S2704 + _S2705.differential_0);
    float3  _S2708 = _S2702.differential_0 + _S2705.differential_0;
    float3  _S2709 = make_float3 (0.5f) * - _S2708;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2710;
    (&_S2710)->primal_0 = _S2680;
    (&_S2710)->differential_0 = _S2693;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2711;
    (&_S2711)->primal_0 = _S2681;
    (&_S2711)->differential_0 = _S2693;
    s_bwd_prop_max_0(&_S2710, &_S2711, (*v_rgb_6)[int(0)]);
    float3  _S2712 = _S2709 + _S2710.differential_0;
    float3  _S2713 = _S2708 + _S2710.differential_0;
    float3  _S2714 = _S2678 * _S2713;
    float3  _S2715 = (*sh_coeffs_19)[int(15)] * _S2713;
    float3  _S2716 = _S2676 * _S2713;
    float3  _S2717 = (*sh_coeffs_19)[int(14)] * _S2713;
    float3  _S2718 = _S2674 * _S2713;
    float3  _S2719 = (*sh_coeffs_19)[int(13)] * _S2713;
    float3  _S2720 = _S2673 * _S2713;
    float3  _S2721 = (*sh_coeffs_19)[int(12)] * _S2713;
    float3  _S2722 = _S2675 * _S2713;
    float3  _S2723 = (*sh_coeffs_19)[int(11)] * _S2713;
    float3  _S2724 = _S2677 * _S2713;
    float3  _S2725 = (*sh_coeffs_19)[int(10)] * _S2713;
    float3  _S2726 = _S2679 * _S2713;
    float3  _S2727 = (*sh_coeffs_19)[int(9)] * _S2713;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S2727.x + _S2727.y + _S2727.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S2715.x + _S2715.y + _S2715.z);
    float _S2728 = _S2725.x + _S2725.y + _S2725.z;
    float _S2729 = _S2717.x + _S2717.y + _S2717.z;
    float _S2730 = _S2723.x + _S2723.y + _S2723.z;
    float _S2731 = _S2719.x + _S2719.y + _S2719.z;
    float _S2732 = _S2721.x + _S2721.y + _S2721.z;
    float _S2733 = - s_diff_fC2_T_5;
    float3  _S2734 = _S2670 * _S2713;
    float3  _S2735 = (*sh_coeffs_19)[int(8)] * _S2713;
    float3  _S2736 = _S2668 * _S2713;
    float3  _S2737 = (*sh_coeffs_19)[int(7)] * _S2713;
    float3  _S2738 = _S2667 * _S2713;
    float3  _S2739 = (*sh_coeffs_19)[int(6)] * _S2713;
    float3  _S2740 = _S2669 * _S2713;
    float3  _S2741 = (*sh_coeffs_19)[int(5)] * _S2713;
    float3  _S2742 = _S2671 * _S2713;
    float3  _S2743 = (*sh_coeffs_19)[int(4)] * _S2713;
    float _S2744 = _S2741.x + _S2741.y + _S2741.z;
    float _S2745 = _S2737.x + _S2737.y + _S2737.z;
    float _S2746 = fTmp1B_19 * _S2728 + x_48 * s_diff_fS2_T_5 + y_22 * _S2733 + 0.54627424478530884f * (_S2743.x + _S2743.y + _S2743.z);
    float _S2747 = fTmp1B_19 * _S2729 + y_22 * s_diff_fS2_T_5 + x_48 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S2735.x + _S2735.y + _S2735.z);
    float _S2748 = y_22 * - _S2747;
    float _S2749 = x_48 * _S2747;
    float _S2750 = z_19 * (1.86588168144226074f * (z_19 * _S2732) + -2.28522896766662598f * (y_22 * _S2730 + x_48 * _S2731) + 0.94617468118667603f * (_S2739.x + _S2739.y + _S2739.z));
    float3  _S2751 = make_float3 (0.48860251903533936f) * _S2713;
    float3  _S2752 = - _S2751;
    float3  _S2753 = _S2661 * _S2752;
    float3  _S2754 = (*sh_coeffs_19)[int(3)] * _S2752;
    float3  _S2755 = _S2663 * _S2751;
    float3  _S2756 = (*sh_coeffs_19)[int(2)] * _S2751;
    float3  _S2757 = _S2665 * _S2751;
    float3  _S2758 = (*sh_coeffs_19)[int(1)] * _S2751;
    float _S2759 = (_S2672 * _S2732 + 1.44530570507049561f * (fS1_19 * _S2728 + fC1_19 * _S2729) + -1.09254848957061768f * (y_22 * _S2744 + x_48 * _S2745) + _S2750 + _S2750 + _S2756.x + _S2756.y + _S2756.z) / _S2662;
    float _S2760 = _S2660 * _S2759;
    float _S2761 = (fTmp0C_19 * _S2730 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S2733 + fTmp0B_19 * _S2744 + _S2666 * _S2746 + _S2748 + _S2748 + - (_S2758.x + _S2758.y + _S2758.z)) / _S2662;
    float _S2762 = _S2660 * _S2761;
    float _S2763 = (fTmp0C_19 * _S2731 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S2745 + 2.0f * (y_22 * _S2746) + _S2749 + _S2749 + _S2754.x + _S2754.y + _S2754.z) / _S2662;
    float _S2764 = _S2660 * _S2763;
    float _S2765 = _S2658 * - _S2759 + _S2657 * - _S2761 + _S2656 * - _S2763;
    DiffPair_float_0 _S2766;
    (&_S2766)->primal_0 = _S2659;
    (&_S2766)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2766, _S2765);
    float _S2767 = _S2658 * _S2766.differential_0;
    float _S2768 = _S2657 * _S2766.differential_0;
    float _S2769 = _S2656 * _S2766.differential_0;
    float3  _S2770 = make_float3 (0.282094806432724f) * _S2713;
    float3  _S2771 = make_float3 (_S2764 + _S2769 + _S2769, _S2762 + _S2768 + _S2768, _S2760 + _S2767 + _S2767);
    float3  _S2772 = - - _S2771;
    Matrix<float, 3, 3>  _S2773 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2774;
    (&_S2774)->primal_0 = _S2654;
    (&_S2774)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2775;
    (&_S2775)->primal_0 = t_22;
    (&_S2775)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2774, &_S2775, _S2772);
    Matrix<float, 3, 3>  _S2776 = transpose_0(_S2774.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2777;
    (&_S2777)->primal_0 = vert2_c_4;
    (&_S2777)->differential_0 = _S2693;
    s_bwd_length_impl_0(&_S2777, v_depth_6.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2778;
    (&_S2778)->primal_0 = vert1_c_4;
    (&_S2778)->differential_0 = _S2693;
    s_bwd_length_impl_0(&_S2778, v_depth_6.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2779;
    (&_S2779)->primal_0 = vert0_c_4;
    (&_S2779)->differential_0 = _S2693;
    s_bwd_length_impl_0(&_S2779, v_depth_6.x);
    DiffPair_float_0 _S2780;
    (&_S2780)->primal_0 = _S2653;
    (&_S2780)->differential_0 = 0.0f;
    DiffPair_float_0 _S2781;
    (&_S2781)->primal_0 = _S2652;
    (&_S2781)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2780, &_S2781, 0.0f);
    DiffPair_float_0 _S2782;
    (&_S2782)->primal_0 = _S2649;
    (&_S2782)->differential_0 = 0.0f;
    DiffPair_float_0 _S2783;
    (&_S2783)->primal_0 = _S2650;
    (&_S2783)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2782, &_S2783, _S2780.differential_0);
    DiffPair_float_0 _S2784;
    (&_S2784)->primal_0 = _S2651;
    (&_S2784)->differential_0 = 0.0f;
    DiffPair_float_0 _S2785;
    (&_S2785)->primal_0 = _S2652;
    (&_S2785)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2784, &_S2785, 0.0f);
    float _S2786 = _S2781.differential_0 + _S2785.differential_0;
    DiffPair_float_0 _S2787;
    (&_S2787)->primal_0 = _S2649;
    (&_S2787)->differential_0 = 0.0f;
    DiffPair_float_0 _S2788;
    (&_S2788)->primal_0 = _S2650;
    (&_S2788)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2787, &_S2788, _S2784.differential_0);
    float _S2789 = _S2783.differential_0 + _S2788.differential_0;
    float _S2790 = _S2782.differential_0 + _S2787.differential_0;
    DiffPair_float_0 _S2791;
    (&_S2791)->primal_0 = _S2648;
    (&_S2791)->differential_0 = 0.0f;
    DiffPair_float_0 _S2792;
    (&_S2792)->primal_0 = _S2647;
    (&_S2792)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2791, &_S2792, 0.0f);
    DiffPair_float_0 _S2793;
    (&_S2793)->primal_0 = _S2644;
    (&_S2793)->differential_0 = 0.0f;
    DiffPair_float_0 _S2794;
    (&_S2794)->primal_0 = _S2645;
    (&_S2794)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2793, &_S2794, _S2791.differential_0);
    DiffPair_float_0 _S2795;
    (&_S2795)->primal_0 = _S2646;
    (&_S2795)->differential_0 = 0.0f;
    DiffPair_float_0 _S2796;
    (&_S2796)->primal_0 = _S2647;
    (&_S2796)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2795, &_S2796, 0.0f);
    float _S2797 = _S2792.differential_0 + _S2796.differential_0;
    DiffPair_float_0 _S2798;
    (&_S2798)->primal_0 = _S2644;
    (&_S2798)->differential_0 = 0.0f;
    DiffPair_float_0 _S2799;
    (&_S2799)->primal_0 = _S2645;
    (&_S2799)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2798, &_S2799, _S2795.differential_0);
    float _S2800 = _S2794.differential_0 + _S2799.differential_0;
    float _S2801 = _S2793.differential_0 + _S2798.differential_0;
    DiffPair_float_0 _S2802;
    (&_S2802)->primal_0 = _S2642;
    (&_S2802)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S2802, 0.0f);
    float _S2803 = - (-1.0f * - (_S2802.differential_0 / _S2643));
    float2  _S2804 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2805;
    (&_S2805)->primal_0 = e2_0;
    (&_S2805)->differential_0 = _S2804;
    s_bwd_length_impl_1(&_S2805, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2806;
    (&_S2806)->primal_0 = e1_4;
    (&_S2806)->differential_0 = _S2804;
    s_bwd_length_impl_1(&_S2806, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S2807;
    (&_S2807)->primal_0 = e0_4;
    (&_S2807)->differential_0 = _S2804;
    s_bwd_length_impl_1(&_S2807, -0.0f);
    DiffPair_float_0 _S2808;
    (&_S2808)->primal_0 = _S2640;
    (&_S2808)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S2808, 0.0f);
    float _S2809 = - _S2808.differential_0;
    float2  _S2810 = _S2806.differential_0 + make_float2 (_S2638 * _S2809, _S2636 * _S2808.differential_0);
    float2  _S2811 = _S2807.differential_0 + make_float2 (_S2637 * _S2808.differential_0, _S2639 * _S2809);
    float2  _S2812 = _S2631 * (v_uv2_0 + - _S2805.differential_0 + _S2810 + make_float2 (_S2797, _S2786)) / _S2630;
    float2  _S2813 = _S2627 * - _S2812;
    float2  _S2814 = _S2629 * _S2812;
    float2  _S2815 = _S2631 * (v_uv1_0 + - _S2810 + _S2811 + make_float2 (_S2800, _S2789)) / _S2626;
    float2  _S2816 = _S2623 * - _S2815;
    float2  _S2817 = _S2625 * _S2815;
    float _S2818 = _S2816.x + _S2816.y;
    float2  _S2819 = _S2631 * (v_uv0_0 + _S2805.differential_0 + - _S2811 + make_float2 (_S2801, _S2790)) / _S2622;
    float2  _S2820 = _S2619 * - _S2819;
    float2  _S2821 = _S2621 * _S2819;
    float _S2822 = _S2820.x + _S2820.y;
    float3  _S2823 = _S2699.differential_0 + _S2777.differential_0 + make_float3 (_S2814.x, _S2814.y, _S2813.x + _S2813.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2824;
    (&_S2824)->primal_0 = R_23;
    (&_S2824)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2825;
    (&_S2825)->primal_0 = vert2_0;
    (&_S2825)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2824, &_S2825, _S2823);
    float3  _S2826 = _S2698.differential_0 + _S2778.differential_0 + make_float3 (_S2817.x, _S2817.y, _S2818);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2827;
    (&_S2827)->primal_0 = R_23;
    (&_S2827)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2828;
    (&_S2828)->primal_0 = vert1_0;
    (&_S2828)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2827, &_S2828, _S2826);
    float3  _S2829 = _S2700 + _S2701 + _S2779.differential_0 + make_float3 (_S2821.x, _S2821.y, _S2822);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2830;
    (&_S2830)->primal_0 = R_23;
    (&_S2830)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2831;
    (&_S2831)->primal_0 = vert0_0;
    (&_S2831)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2830, &_S2831, _S2829);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2832;
    (&_S2832)->primal_0 = _S2613;
    (&_S2832)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2833;
    (&_S2833)->primal_0 = _S2618;
    (&_S2833)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2832, &_S2833, _S2825.differential_0);
    float _S2834 = - _S2833.differential_0.y;
    float _S2835 = _S2617 * _S2833.differential_0.x;
    float _S2836 = - (_S2609 * _S2833.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2837;
    (&_S2837)->primal_0 = _S2613;
    (&_S2837)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2838;
    (&_S2838)->primal_0 = _S2616;
    (&_S2838)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2837, &_S2838, _S2828.differential_0);
    float _S2839 = _S2609 * _S2838.differential_0.x;
    float _S2840 = _S2615 * _S2838.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2841;
    (&_S2841)->primal_0 = _S2613;
    (&_S2841)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2842;
    (&_S2842)->primal_0 = _S2614;
    (&_S2842)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2841, &_S2842, _S2831.differential_0);
    Matrix<float, 3, 3>  _S2843 = transpose_0(_S2832.differential_0 + _S2837.differential_0 + _S2841.differential_0);
    float _S2844 = 2.0f * - _S2843.rows[int(2)].z;
    float _S2845 = 2.0f * _S2843.rows[int(2)].y;
    float _S2846 = 2.0f * _S2843.rows[int(2)].x;
    float _S2847 = 2.0f * _S2843.rows[int(1)].z;
    float _S2848 = 2.0f * - _S2843.rows[int(1)].y;
    float _S2849 = 2.0f * _S2843.rows[int(1)].x;
    float _S2850 = 2.0f * _S2843.rows[int(0)].z;
    float _S2851 = 2.0f * _S2843.rows[int(0)].y;
    float _S2852 = 2.0f * - _S2843.rows[int(0)].x;
    float _S2853 = - _S2849 + _S2851;
    float _S2854 = _S2846 + - _S2850;
    float _S2855 = - _S2845 + _S2847;
    float _S2856 = _S2845 + _S2847;
    float _S2857 = _S2846 + _S2850;
    float _S2858 = _S2849 + _S2851;
    float _S2859 = quat_24.w * (_S2848 + _S2852);
    float _S2860 = quat_24.z * (_S2844 + _S2852);
    float _S2861 = quat_24.y * (_S2844 + _S2848);
    float _S2862 = quat_24.x * _S2853 + quat_24.z * _S2856 + quat_24.y * _S2857 + _S2859 + _S2859;
    float _S2863 = quat_24.x * _S2854 + quat_24.w * _S2856 + quat_24.y * _S2858 + _S2860 + _S2860;
    float _S2864 = quat_24.x * _S2855 + quat_24.w * _S2857 + quat_24.z * _S2858 + _S2861 + _S2861;
    float _S2865 = quat_24.w * _S2853 + quat_24.z * _S2854 + quat_24.y * _S2855;
    float _S2866 = _S2836 + _S2839;
    float _S2867 = 0.5f * - _S2866;
    float _S2868 = _S2834 + _S2838.differential_0.y;
    DiffPair_float_0 _S2869;
    (&_S2869)->primal_0 = _S2610;
    (&_S2869)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2869, _S2868);
    float _S2870 = _S2867 + _S2869.differential_0;
    float _S2871 = _S2835 + _S2840 + _S2842.differential_0.x;
    DiffPair_float_0 _S2872;
    (&_S2872)->primal_0 = _S2608;
    (&_S2872)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2872, _S2871);
    float _S2873 = _S2867 + _S2872.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2874;
    (&_S2874)->primal_0 = R_23;
    (&_S2874)->differential_0 = _S2773;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2875;
    (&_S2875)->primal_0 = mean_23;
    (&_S2875)->differential_0 = _S2693;
    s_bwd_prop_mul_1(&_S2874, &_S2875, _S2695.differential_0);
    float3  _S2876 = _S2775.differential_0 + _S2823 + _S2826 + _S2829 + _S2695.differential_0;
    Matrix<float, 3, 3>  _S2877 = _S2776 + _S2824.differential_0 + _S2827.differential_0 + _S2830.differential_0 + _S2874.differential_0;
    FixedArray<float3 , 2>  _S2878;
    _S2878[int(0)] = _S2693;
    _S2878[int(1)] = _S2693;
    _S2878[int(1)] = _S2707;
    _S2878[int(0)] = _S2712;
    FixedArray<float3 , 16>  _S2879;
    _S2879[int(0)] = _S2693;
    _S2879[int(1)] = _S2693;
    _S2879[int(2)] = _S2693;
    _S2879[int(3)] = _S2693;
    _S2879[int(4)] = _S2693;
    _S2879[int(5)] = _S2693;
    _S2879[int(6)] = _S2693;
    _S2879[int(7)] = _S2693;
    _S2879[int(8)] = _S2693;
    _S2879[int(9)] = _S2693;
    _S2879[int(10)] = _S2693;
    _S2879[int(11)] = _S2693;
    _S2879[int(12)] = _S2693;
    _S2879[int(13)] = _S2693;
    _S2879[int(14)] = _S2693;
    _S2879[int(15)] = _S2693;
    _S2879[int(15)] = _S2714;
    _S2879[int(14)] = _S2716;
    _S2879[int(13)] = _S2718;
    _S2879[int(12)] = _S2720;
    _S2879[int(11)] = _S2722;
    _S2879[int(10)] = _S2724;
    _S2879[int(9)] = _S2726;
    _S2879[int(8)] = _S2734;
    _S2879[int(7)] = _S2736;
    _S2879[int(6)] = _S2738;
    _S2879[int(5)] = _S2740;
    _S2879[int(4)] = _S2742;
    _S2879[int(3)] = _S2753;
    _S2879[int(2)] = _S2755;
    _S2879[int(1)] = _S2757;
    _S2879[int(0)] = _S2770;
    float2  _S2880 = v_out_hardness_0 + make_float2 (0.0f, _S2803);
    float3  _S2881 = make_float3 (_S2873, _S2870, _S2866);
    float4  _S2882 = make_float4 (0.0f);
    *&((&_S2882)->w) = _S2862;
    *&((&_S2882)->z) = _S2863;
    *&((&_S2882)->y) = _S2864;
    *&((&_S2882)->x) = _S2865;
    *v_mean_7 = _S2771 + _S2825.differential_0 + _S2828.differential_0 + _S2831.differential_0 + _S2875.differential_0;
    *v_quat_6 = _S2882;
    *v_scale_6 = _S2881;
    *v_hardness_0 = _S2880;
    *v_sh_coeffs_5 = _S2879;
    *v_ch_coeffs_0 = _S2878;
    *v_R_6 = _S2877;
    *v_t_6 = _S2876;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_24, float3  t_23, float fx_24, float fy_24, float cx_24, float cy_24, float4  radial_coeffs_27, float2  tangential_coeffs_27, float2  thin_prism_coeffs_27, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_1(R_24, mean_24) + t_23;
    float _S2883 = scale_24.x;
    float _S2884 = s_primal_ctx_exp_1(_S2883);
    float _S2885 = scale_24.y;
    float _S2886 = s_primal_ctx_exp_1(_S2885);
    float sz_5 = scale_24.z - 0.5f * (_S2883 + _S2885);
    float _S2887 = quat_25.y;
    float x2_25 = _S2887 * _S2887;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_45 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S2888 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_45), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_45), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S2889 = make_float3 (_S2884, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S2888, _S2889) + mean_24;
    float _S2890 = -0.5f + sz_5;
    float3  _S2891 = make_float3 (_S2884 * _S2890, _S2886, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S2888, _S2891) + mean_24;
    float _S2892 = -0.5f - sz_5;
    float3  _S2893 = make_float3 (_S2884 * _S2892, - _S2886, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S2888, _S2893) + mean_24;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_24, vert0_1) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_24, vert1_1) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_24, vert2_1) + t_23;
    CameraDistortion_0 _S2894 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_27, tangential_coeffs_27, thin_prism_coeffs_27);
    float2  _S2895 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S2896 = length_0(_S2895);
    float _S2897 = vert0_c_5.z;
    float _S2898 = s_primal_ctx_atan2_0(_S2896, _S2897);
    bool _S2899 = _S2898 < 0.00100000004749745f;
    float k_6;
    float _S2900;
    float _S2901;
    float _S2902;
    if(_S2899)
    {
        float _S2903 = 1.0f - _S2898 * _S2898 / 3.0f;
        float _S2904 = _S2897 * _S2897;
        k_6 = _S2903 / _S2897;
        _S2900 = 0.0f;
        _S2901 = _S2904;
        _S2902 = _S2903;
    }
    else
    {
        float _S2905 = _S2896 * _S2896;
        k_6 = _S2898 / _S2896;
        _S2900 = _S2905;
        _S2901 = 0.0f;
        _S2902 = 0.0f;
    }
    float2  _S2906 = make_float2 (k_6);
    float2  _S2907 = _S2895 * make_float2 (k_6);
    float k1_5 = _S2894.radial_coeffs_0.x;
    float k2_5 = _S2894.radial_coeffs_0.y;
    float k3_5 = _S2894.radial_coeffs_0.z;
    float k4_6 = _S2894.radial_coeffs_0.w;
    float p1_6 = _S2894.tangential_coeffs_0.x;
    float p2_6 = _S2894.tangential_coeffs_0.y;
    float sx1_6 = _S2894.thin_prism_coeffs_0.x;
    float sy1_6 = _S2894.thin_prism_coeffs_0.y;
    float u_16 = _S2907.x;
    float v_16 = _S2907.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float _S2908 = k3_5 + r2_16 * k4_6;
    float _S2909 = k2_5 + r2_16 * _S2908;
    float _S2910 = k1_5 + r2_16 * _S2909;
    float radial_2 = 1.0f + r2_16 * _S2910;
    float2  _S2911 = make_float2 (radial_2);
    float _S2912 = 2.0f * p1_6;
    float _S2913 = _S2912 * u_16;
    float _S2914 = 2.0f * u_16;
    float _S2915 = r2_16 + _S2914 * u_16;
    float _S2916 = 2.0f * p2_6;
    float _S2917 = _S2916 * u_16;
    float _S2918 = 2.0f * v_16;
    float _S2919 = r2_16 + _S2918 * v_16;
    float2  _S2920 = _S2907 * make_float2 (radial_2) + make_float2 (_S2913 * v_16 + p2_6 * _S2915 + sx1_6 * r2_16, _S2917 * v_16 + p1_6 * _S2919 + sy1_6 * r2_16);
    float _S2921 = fx_24 * _S2920.x + cx_24;
    float _S2922 = fy_24 * _S2920.y + cy_24;
    float2  _S2923 = make_float2 (_S2921, _S2922);
    float2  _S2924 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S2925 = length_0(_S2924);
    float _S2926 = vert1_c_5.z;
    float _S2927 = s_primal_ctx_atan2_0(_S2925, _S2926);
    bool _S2928 = _S2927 < 0.00100000004749745f;
    float _S2929;
    float _S2930;
    float _S2931;
    if(_S2928)
    {
        float _S2932 = 1.0f - _S2927 * _S2927 / 3.0f;
        float _S2933 = _S2926 * _S2926;
        k_6 = _S2932 / _S2926;
        _S2929 = 0.0f;
        _S2930 = _S2933;
        _S2931 = _S2932;
    }
    else
    {
        float _S2934 = _S2925 * _S2925;
        k_6 = _S2927 / _S2925;
        _S2929 = _S2934;
        _S2930 = 0.0f;
        _S2931 = 0.0f;
    }
    float2  _S2935 = make_float2 (k_6);
    float2  _S2936 = _S2924 * make_float2 (k_6);
    float u_17 = _S2936.x;
    float v_17 = _S2936.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float _S2937 = k3_5 + r2_17 * k4_6;
    float _S2938 = k2_5 + r2_17 * _S2937;
    float _S2939 = k1_5 + r2_17 * _S2938;
    float radial_3 = 1.0f + r2_17 * _S2939;
    float2  _S2940 = make_float2 (radial_3);
    float _S2941 = _S2912 * u_17;
    float _S2942 = 2.0f * u_17;
    float _S2943 = r2_17 + _S2942 * u_17;
    float _S2944 = _S2916 * u_17;
    float _S2945 = 2.0f * v_17;
    float _S2946 = r2_17 + _S2945 * v_17;
    float2  _S2947 = _S2936 * make_float2 (radial_3) + make_float2 (_S2941 * v_17 + p2_6 * _S2943 + sx1_6 * r2_17, _S2944 * v_17 + p1_6 * _S2946 + sy1_6 * r2_17);
    float _S2948 = fx_24 * _S2947.x + cx_24;
    float _S2949 = fy_24 * _S2947.y + cy_24;
    float2  _S2950 = make_float2 (_S2948, _S2949);
    float2  _S2951 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S2952 = length_0(_S2951);
    float _S2953 = vert2_c_5.z;
    float _S2954 = s_primal_ctx_atan2_0(_S2952, _S2953);
    bool _S2955 = _S2954 < 0.00100000004749745f;
    float _S2956;
    float _S2957;
    float _S2958;
    if(_S2955)
    {
        float _S2959 = 1.0f - _S2954 * _S2954 / 3.0f;
        float _S2960 = _S2953 * _S2953;
        k_6 = _S2959 / _S2953;
        _S2956 = 0.0f;
        _S2957 = _S2960;
        _S2958 = _S2959;
    }
    else
    {
        float _S2961 = _S2952 * _S2952;
        k_6 = _S2954 / _S2952;
        _S2956 = _S2961;
        _S2957 = 0.0f;
        _S2958 = 0.0f;
    }
    float2  _S2962 = make_float2 (k_6);
    float2  _S2963 = _S2951 * make_float2 (k_6);
    float u_18 = _S2963.x;
    float v_18 = _S2963.y;
    float r2_18 = u_18 * u_18 + v_18 * v_18;
    float _S2964 = k3_5 + r2_18 * k4_6;
    float _S2965 = k2_5 + r2_18 * _S2964;
    float _S2966 = k1_5 + r2_18 * _S2965;
    float radial_4 = 1.0f + r2_18 * _S2966;
    float2  _S2967 = make_float2 (radial_4);
    float _S2968 = _S2912 * u_18;
    float _S2969 = 2.0f * u_18;
    float _S2970 = r2_18 + _S2969 * u_18;
    float _S2971 = _S2916 * u_18;
    float _S2972 = 2.0f * v_18;
    float _S2973 = r2_18 + _S2972 * v_18;
    float2  _S2974 = _S2963 * make_float2 (radial_4) + make_float2 (_S2968 * v_18 + p2_6 * _S2970 + sx1_6 * r2_18, _S2971 * v_18 + p1_6 * _S2973 + sy1_6 * r2_18);
    float _S2975 = fx_24 * _S2974.x + cx_24;
    float _S2976 = fy_24 * _S2974.y + cy_24;
    float2  _S2977 = make_float2 (_S2975, _S2976);
    float2  e0_5 = _S2950 - _S2923;
    float2  e1_5 = _S2977 - _S2950;
    float2  e2_1 = _S2923 - _S2977;
    float _S2978 = e0_5.x;
    float _S2979 = e1_5.y;
    float _S2980 = e0_5.y;
    float _S2981 = e1_5.x;
    float _S2982 = _S2978 * _S2979 - _S2980 * _S2981;
    float _S2983 = 1.0f - hardness_5.y;
    float _S2984 = -1.0f / _S2983;
    float _S2985 = _S2983 * _S2983;
    float _S2986 = s_primal_ctx_max_0(_S2921, _S2948);
    float _S2987 = s_primal_ctx_min_0(_S2921, _S2948);
    float _S2988 = s_primal_ctx_max_0(_S2922, _S2949);
    float _S2989 = s_primal_ctx_min_0(_S2922, _S2949);
    Matrix<float, 3, 3>  _S2990 = transpose_0(R_24);
    float3  _S2991 = mean_24 - - s_primal_ctx_mul_1(_S2990, t_23);
    float _S2992 = _S2991.x;
    float _S2993 = _S2991.y;
    float _S2994 = _S2991.z;
    float _S2995 = _S2992 * _S2992 + _S2993 * _S2993 + _S2994 * _S2994;
    float _S2996 = s_primal_ctx_sqrt_0(_S2995);
    float x_49 = _S2992 / _S2996;
    float3  _S2997 = make_float3 (x_49);
    float _S2998 = _S2996 * _S2996;
    float y_23 = _S2993 / _S2996;
    float z_20 = _S2994 / _S2996;
    float3  _S2999 = make_float3 (z_20);
    float _S3000 = - y_23;
    float3  _S3001 = make_float3 (_S3000);
    float z2_46 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_49 * x_49 - y_23 * y_23;
    float _S3002 = 2.0f * x_49;
    float fS1_20 = _S3002 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_46 - 0.31539157032966614f;
    float3  _S3003 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_49;
    float3  _S3004 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3005 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3006 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3007 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_46 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3008 = 1.86588168144226074f * z2_46 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3008;
    float3  _S3009 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_49;
    float3  _S3010 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3011 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3012 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3013 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_49 * fC1_20 - y_23 * fS1_20);
    float3  _S3014 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_49 * fS1_20 + y_23 * fC1_20);
    float3  _S3015 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3000) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_49) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3016 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3017 = make_float3 (0.0f);
    float3  _S3018 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3019 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3020 = make_float3 (_S3019);
    float3  _S3021 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3019);
    float3  _S3022 = _S3018 + _S3021 + make_float3 (0.5f);
    float3  _S3023 = _S3018 - _S3021 + make_float3 (0.5f);
    float3  _S3024 = vert1_c_5 - vert0_c_5;
    float3  _S3025 = vert2_c_5 - vert0_c_5;
    float3  _S3026 = s_primal_ctx_cross_0(_S3024, _S3025);
    float3  _S3027 = normalize_0(_S3026);
    float3  _S3028 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3027, mean_c_20)))))) * v_normal_1;
    float3  _S3029 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3030;
    (&_S3030)->primal_0 = _S3027;
    (&_S3030)->differential_0 = _S3029;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3031;
    (&_S3031)->primal_0 = mean_c_20;
    (&_S3031)->differential_0 = _S3029;
    s_bwd_prop_dot_0(&_S3030, &_S3031, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3032 = _S3031;
    float3  _S3033 = _S3028 + _S3030.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3034;
    (&_S3034)->primal_0 = _S3026;
    (&_S3034)->differential_0 = _S3029;
    s_bwd_normalize_impl_0(&_S3034, _S3033);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3035;
    (&_S3035)->primal_0 = _S3024;
    (&_S3035)->differential_0 = _S3029;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3036;
    (&_S3036)->primal_0 = _S3025;
    (&_S3036)->differential_0 = _S3029;
    s_bwd_prop_cross_0(&_S3035, &_S3036, _S3034.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3037 = _S3035;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3038 = _S3036;
    float3  _S3039 = - _S3036.differential_0;
    float3  _S3040 = - _S3035.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3041;
    (&_S3041)->primal_0 = _S3023;
    (&_S3041)->differential_0 = _S3029;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3042;
    (&_S3042)->primal_0 = _S3017;
    (&_S3042)->differential_0 = _S3029;
    s_bwd_prop_max_0(&_S3041, &_S3042, (*v_rgb_7)[int(2)]);
    float3  _S3043 = - _S3041.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3044;
    (&_S3044)->primal_0 = _S3022;
    (&_S3044)->differential_0 = _S3029;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3045;
    (&_S3045)->primal_0 = _S3017;
    (&_S3045)->differential_0 = _S3029;
    s_bwd_prop_max_0(&_S3044, &_S3045, (*v_rgb_7)[int(1)]);
    float3  _S3046 = _S3020 * (_S3043 + _S3044.differential_0);
    float3  _S3047 = _S3041.differential_0 + _S3044.differential_0;
    float3  _S3048 = make_float3 (0.5f) * - _S3047;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3049;
    (&_S3049)->primal_0 = _S3016;
    (&_S3049)->differential_0 = _S3029;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3050;
    (&_S3050)->primal_0 = _S3017;
    (&_S3050)->differential_0 = _S3029;
    s_bwd_prop_max_0(&_S3049, &_S3050, (*v_rgb_7)[int(0)]);
    float3  _S3051 = _S3048 + _S3049.differential_0;
    float3  _S3052 = _S3047 + _S3049.differential_0;
    float3  _S3053 = _S3014 * _S3052;
    float3  _S3054 = (*sh_coeffs_20)[int(15)] * _S3052;
    float3  _S3055 = _S3012 * _S3052;
    float3  _S3056 = (*sh_coeffs_20)[int(14)] * _S3052;
    float3  _S3057 = _S3010 * _S3052;
    float3  _S3058 = (*sh_coeffs_20)[int(13)] * _S3052;
    float3  _S3059 = _S3009 * _S3052;
    float3  _S3060 = (*sh_coeffs_20)[int(12)] * _S3052;
    float3  _S3061 = _S3011 * _S3052;
    float3  _S3062 = (*sh_coeffs_20)[int(11)] * _S3052;
    float3  _S3063 = _S3013 * _S3052;
    float3  _S3064 = (*sh_coeffs_20)[int(10)] * _S3052;
    float3  _S3065 = _S3015 * _S3052;
    float3  _S3066 = (*sh_coeffs_20)[int(9)] * _S3052;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3066.x + _S3066.y + _S3066.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3054.x + _S3054.y + _S3054.z);
    float _S3067 = _S3064.x + _S3064.y + _S3064.z;
    float _S3068 = _S3056.x + _S3056.y + _S3056.z;
    float _S3069 = _S3062.x + _S3062.y + _S3062.z;
    float _S3070 = _S3058.x + _S3058.y + _S3058.z;
    float _S3071 = _S3060.x + _S3060.y + _S3060.z;
    float _S3072 = - s_diff_fC2_T_6;
    float3  _S3073 = _S3006 * _S3052;
    float3  _S3074 = (*sh_coeffs_20)[int(8)] * _S3052;
    float3  _S3075 = _S3004 * _S3052;
    float3  _S3076 = (*sh_coeffs_20)[int(7)] * _S3052;
    float3  _S3077 = _S3003 * _S3052;
    float3  _S3078 = (*sh_coeffs_20)[int(6)] * _S3052;
    float3  _S3079 = _S3005 * _S3052;
    float3  _S3080 = (*sh_coeffs_20)[int(5)] * _S3052;
    float3  _S3081 = _S3007 * _S3052;
    float3  _S3082 = (*sh_coeffs_20)[int(4)] * _S3052;
    float _S3083 = _S3080.x + _S3080.y + _S3080.z;
    float _S3084 = _S3076.x + _S3076.y + _S3076.z;
    float _S3085 = fTmp1B_20 * _S3067 + x_49 * s_diff_fS2_T_6 + y_23 * _S3072 + 0.54627424478530884f * (_S3082.x + _S3082.y + _S3082.z);
    float _S3086 = fTmp1B_20 * _S3068 + y_23 * s_diff_fS2_T_6 + x_49 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3074.x + _S3074.y + _S3074.z);
    float _S3087 = y_23 * - _S3086;
    float _S3088 = x_49 * _S3086;
    float _S3089 = z_20 * (1.86588168144226074f * (z_20 * _S3071) + -2.28522896766662598f * (y_23 * _S3069 + x_49 * _S3070) + 0.94617468118667603f * (_S3078.x + _S3078.y + _S3078.z));
    float3  _S3090 = make_float3 (0.48860251903533936f) * _S3052;
    float3  _S3091 = - _S3090;
    float3  _S3092 = _S2997 * _S3091;
    float3  _S3093 = (*sh_coeffs_20)[int(3)] * _S3091;
    float3  _S3094 = _S2999 * _S3090;
    float3  _S3095 = (*sh_coeffs_20)[int(2)] * _S3090;
    float3  _S3096 = _S3001 * _S3090;
    float3  _S3097 = (*sh_coeffs_20)[int(1)] * _S3090;
    float _S3098 = (_S3008 * _S3071 + 1.44530570507049561f * (fS1_20 * _S3067 + fC1_20 * _S3068) + -1.09254848957061768f * (y_23 * _S3083 + x_49 * _S3084) + _S3089 + _S3089 + _S3095.x + _S3095.y + _S3095.z) / _S2998;
    float _S3099 = _S2996 * _S3098;
    float _S3100 = (fTmp0C_20 * _S3069 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3072 + fTmp0B_20 * _S3083 + _S3002 * _S3085 + _S3087 + _S3087 + - (_S3097.x + _S3097.y + _S3097.z)) / _S2998;
    float _S3101 = _S2996 * _S3100;
    float _S3102 = (fTmp0C_20 * _S3070 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3084 + 2.0f * (y_23 * _S3085) + _S3088 + _S3088 + _S3093.x + _S3093.y + _S3093.z) / _S2998;
    float _S3103 = _S2996 * _S3102;
    float _S3104 = _S2994 * - _S3098 + _S2993 * - _S3100 + _S2992 * - _S3102;
    DiffPair_float_0 _S3105;
    (&_S3105)->primal_0 = _S2995;
    (&_S3105)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3105, _S3104);
    float _S3106 = _S2994 * _S3105.differential_0;
    float _S3107 = _S2993 * _S3105.differential_0;
    float _S3108 = _S2992 * _S3105.differential_0;
    float3  _S3109 = make_float3 (0.282094806432724f) * _S3052;
    float3  _S3110 = make_float3 (_S3103 + _S3108 + _S3108, _S3101 + _S3107 + _S3107, _S3099 + _S3106 + _S3106);
    float3  _S3111 = - - _S3110;
    Matrix<float, 3, 3>  _S3112 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3113;
    (&_S3113)->primal_0 = _S2990;
    (&_S3113)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3114;
    (&_S3114)->primal_0 = t_23;
    (&_S3114)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3113, &_S3114, _S3111);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3115 = _S3114;
    Matrix<float, 3, 3>  _S3116 = transpose_0(_S3113.differential_0);
    DiffPair_float_0 _S3117;
    (&_S3117)->primal_0 = _S2989;
    (&_S3117)->differential_0 = 0.0f;
    DiffPair_float_0 _S3118;
    (&_S3118)->primal_0 = _S2976;
    (&_S3118)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3117, &_S3118, 0.0f);
    DiffPair_float_0 _S3119;
    (&_S3119)->primal_0 = _S2922;
    (&_S3119)->differential_0 = 0.0f;
    DiffPair_float_0 _S3120;
    (&_S3120)->primal_0 = _S2949;
    (&_S3120)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3119, &_S3120, _S3117.differential_0);
    DiffPair_float_0 _S3121;
    (&_S3121)->primal_0 = _S2988;
    (&_S3121)->differential_0 = 0.0f;
    DiffPair_float_0 _S3122;
    (&_S3122)->primal_0 = _S2976;
    (&_S3122)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3121, &_S3122, 0.0f);
    DiffPair_float_0 _S3123;
    (&_S3123)->primal_0 = _S2922;
    (&_S3123)->differential_0 = 0.0f;
    DiffPair_float_0 _S3124;
    (&_S3124)->primal_0 = _S2949;
    (&_S3124)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3123, &_S3124, _S3121.differential_0);
    DiffPair_float_0 _S3125;
    (&_S3125)->primal_0 = _S2987;
    (&_S3125)->differential_0 = 0.0f;
    DiffPair_float_0 _S3126;
    (&_S3126)->primal_0 = _S2975;
    (&_S3126)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3125, &_S3126, 0.0f);
    DiffPair_float_0 _S3127;
    (&_S3127)->primal_0 = _S2921;
    (&_S3127)->differential_0 = 0.0f;
    DiffPair_float_0 _S3128;
    (&_S3128)->primal_0 = _S2948;
    (&_S3128)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3127, &_S3128, _S3125.differential_0);
    DiffPair_float_0 _S3129;
    (&_S3129)->primal_0 = _S2986;
    (&_S3129)->differential_0 = 0.0f;
    DiffPair_float_0 _S3130;
    (&_S3130)->primal_0 = _S2975;
    (&_S3130)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3129, &_S3130, 0.0f);
    DiffPair_float_0 _S3131;
    (&_S3131)->primal_0 = _S2921;
    (&_S3131)->differential_0 = 0.0f;
    DiffPair_float_0 _S3132;
    (&_S3132)->primal_0 = _S2948;
    (&_S3132)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3131, &_S3132, _S3129.differential_0);
    DiffPair_float_0 _S3133;
    (&_S3133)->primal_0 = _S2984;
    (&_S3133)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3133, 0.0f);
    float _S3134 = - (-1.0f * - (_S3133.differential_0 / _S2985));
    float2  _S3135 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3136;
    (&_S3136)->primal_0 = e2_1;
    (&_S3136)->differential_0 = _S3135;
    s_bwd_length_impl_1(&_S3136, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3137;
    (&_S3137)->primal_0 = e1_5;
    (&_S3137)->differential_0 = _S3135;
    s_bwd_length_impl_1(&_S3137, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3138;
    (&_S3138)->primal_0 = e0_5;
    (&_S3138)->differential_0 = _S3135;
    s_bwd_length_impl_1(&_S3138, -0.0f);
    DiffPair_float_0 _S3139;
    (&_S3139)->primal_0 = _S2982;
    (&_S3139)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3139, 0.0f);
    float _S3140 = - _S3139.differential_0;
    float2  _S3141 = _S3137.differential_0 + make_float2 (_S2980 * _S3140, _S2978 * _S3139.differential_0);
    float2  _S3142 = _S3138.differential_0 + make_float2 (_S2979 * _S3139.differential_0, _S2981 * _S3140);
    float2  _S3143 = v_uv2_1 + - _S3136.differential_0 + _S3141;
    float _S3144 = fy_24 * (_S3118.differential_0 + _S3122.differential_0 + _S3143.y);
    float _S3145 = fx_24 * (_S3126.differential_0 + _S3130.differential_0 + _S3143.x);
    float2  _S3146 = make_float2 (_S3145, _S3144);
    float2  _S3147 = _S2963 * _S3146;
    float2  _S3148 = _S2967 * _S3146;
    float _S3149 = r2_18 * _S3144;
    float _S3150 = p1_6 * _S3144;
    float _S3151 = _S2973 * _S3144;
    float _S3152 = v_18 * _S3144;
    float _S3153 = u_18 * _S3152;
    float _S3154 = r2_18 * _S3145;
    float _S3155 = p2_6 * _S3145;
    float _S3156 = _S2970 * _S3145;
    float _S3157 = v_18 * _S3145;
    float _S3158 = u_18 * _S3157;
    float _S3159 = _S3147.x + _S3147.y;
    float _S3160 = r2_18 * _S3159;
    float _S3161 = r2_18 * _S3160;
    float _S3162 = r2_18 * _S3161;
    float _S3163 = r2_18 * _S3162;
    float _S3164 = sy1_6 * _S3144 + _S3150 + sx1_6 * _S3145 + _S3155 + _S2966 * _S3159 + _S2965 * _S3160 + _S2964 * _S3161 + k4_6 * _S3162;
    float _S3165 = v_18 * _S3164;
    float _S3166 = u_18 * _S3164;
    float _S3167 = _S2972 * _S3150 + 2.0f * (v_18 * _S3150) + _S2971 * _S3144 + _S2968 * _S3145 + _S3165 + _S3165;
    float _S3168 = _S2916 * _S3152 + _S2969 * _S3155 + 2.0f * (u_18 * _S3155) + _S2912 * _S3157 + _S3166 + _S3166;
    float2  _S3169 = v_uv0_1 + _S3136.differential_0 + - _S3142;
    float2  _S3170 = v_out_hardness_1 + make_float2 (0.0f, _S3134);
    float _S3171 = _S3128.differential_0 + _S3132.differential_0;
    float2  _S3172 = v_uv1_1 + - _S3141 + _S3142;
    float3  _S3173 = _S3039 + _S3040;
    FixedArray<float3 , 2>  _S3174;
    _S3174[int(0)] = _S3029;
    _S3174[int(1)] = _S3029;
    _S3174[int(1)] = _S3046;
    _S3174[int(0)] = _S3051;
    float3  _S3175 = _S3174[int(0)];
    float3  _S3176 = _S3174[int(1)];
    FixedArray<float3 , 16>  _S3177;
    _S3177[int(0)] = _S3029;
    _S3177[int(1)] = _S3029;
    _S3177[int(2)] = _S3029;
    _S3177[int(3)] = _S3029;
    _S3177[int(4)] = _S3029;
    _S3177[int(5)] = _S3029;
    _S3177[int(6)] = _S3029;
    _S3177[int(7)] = _S3029;
    _S3177[int(8)] = _S3029;
    _S3177[int(9)] = _S3029;
    _S3177[int(10)] = _S3029;
    _S3177[int(11)] = _S3029;
    _S3177[int(12)] = _S3029;
    _S3177[int(13)] = _S3029;
    _S3177[int(14)] = _S3029;
    _S3177[int(15)] = _S3029;
    _S3177[int(7)] = _S3075;
    _S3177[int(0)] = _S3109;
    _S3177[int(1)] = _S3096;
    _S3177[int(2)] = _S3094;
    _S3177[int(3)] = _S3092;
    _S3177[int(4)] = _S3081;
    _S3177[int(5)] = _S3079;
    _S3177[int(6)] = _S3077;
    _S3177[int(15)] = _S3053;
    _S3177[int(8)] = _S3073;
    _S3177[int(9)] = _S3065;
    _S3177[int(10)] = _S3063;
    _S3177[int(11)] = _S3061;
    _S3177[int(12)] = _S3059;
    _S3177[int(13)] = _S3057;
    _S3177[int(14)] = _S3055;
    float3  _S3178 = _S3177[int(0)];
    float3  _S3179 = _S3177[int(1)];
    float3  _S3180 = _S3177[int(2)];
    float3  _S3181 = _S3177[int(3)];
    float3  _S3182 = _S3177[int(4)];
    float3  _S3183 = _S3177[int(5)];
    float3  _S3184 = _S3177[int(6)];
    float3  _S3185 = _S3177[int(7)];
    float3  _S3186 = _S3177[int(8)];
    float3  _S3187 = _S3177[int(9)];
    float3  _S3188 = _S3177[int(10)];
    float3  _S3189 = _S3177[int(11)];
    float3  _S3190 = _S3177[int(12)];
    float3  _S3191 = _S3177[int(13)];
    float3  _S3192 = _S3177[int(14)];
    float3  _S3193 = _S3177[int(15)];
    float _S3194 = _S3127.differential_0 + _S3131.differential_0;
    float _S3195 = _S3119.differential_0 + _S3123.differential_0;
    float _S3196 = _S3120.differential_0 + _S3124.differential_0;
    float2  _S3197 = _S3148 + make_float2 (_S3168, _S3167);
    float2  _S3198 = _S2951 * _S3197;
    float2  _S3199 = _S2962 * _S3197;
    float _S3200 = _S3198.x + _S3198.y;
    if(_S2955)
    {
        float _S3201 = _S3200 / _S2957;
        float _S3202 = _S2958 * - _S3201;
        float _S3203 = _S2954 * (0.3333333432674408f * - (_S2953 * _S3201));
        k_6 = _S3203 + _S3203;
        _S2956 = _S3202;
        _S2957 = 0.0f;
    }
    else
    {
        float _S3204 = _S3200 / _S2956;
        float _S3205 = _S2954 * - _S3204;
        k_6 = _S2952 * _S3204;
        _S2956 = 0.0f;
        _S2957 = _S3205;
    }
    DiffPair_float_0 _S3206;
    (&_S3206)->primal_0 = _S2952;
    (&_S3206)->differential_0 = 0.0f;
    DiffPair_float_0 _S3207;
    (&_S3207)->primal_0 = _S2953;
    (&_S3207)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3206, &_S3207, k_6);
    float _S3208 = _S3207.differential_0 + _S2956;
    float _S3209 = _S3206.differential_0 + _S2957;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3210;
    (&_S3210)->primal_0 = _S2951;
    (&_S3210)->differential_0 = _S3135;
    s_bwd_length_impl_1(&_S3210, _S3209);
    float2  _S3211 = _S3210.differential_0 + _S3199;
    float _S3212 = fy_24 * (_S3172.y + _S3196);
    float _S3213 = fx_24 * (_S3172.x + _S3171);
    float2  _S3214 = make_float2 (_S3213, _S3212);
    float2  _S3215 = _S2936 * _S3214;
    float _S3216 = p1_6 * _S3212;
    float _S3217 = v_17 * _S3212;
    float _S3218 = p2_6 * _S3213;
    float _S3219 = v_17 * _S3213;
    float _S3220 = _S3215.x + _S3215.y;
    float _S3221 = r2_17 * _S3220;
    float _S3222 = r2_17 * _S3221;
    float _S3223 = r2_17 * _S3222;
    float _S3224 = sy1_6 * _S3212 + _S3216 + sx1_6 * _S3213 + _S3218 + _S2939 * _S3220 + _S2938 * _S3221 + _S2937 * _S3222 + k4_6 * _S3223;
    float _S3225 = v_17 * _S3224;
    float _S3226 = u_17 * _S3224;
    float3  _S3227 = _S3038.differential_0 + make_float3 (_S3211.x, _S3211.y, _S3208);
    float2  _S3228 = _S2940 * _S3214 + make_float2 (_S2916 * _S3217 + _S2942 * _S3218 + 2.0f * (u_17 * _S3218) + _S2912 * _S3219 + _S3226 + _S3226, _S2945 * _S3216 + 2.0f * (v_17 * _S3216) + _S2944 * _S3212 + _S2941 * _S3213 + _S3225 + _S3225);
    float _S3229 = u_17 * _S3217 + _S3153;
    float _S3230 = u_17 * _S3219 + _S3158;
    float _S3231 = r2_17 * _S3212 + _S3149;
    float _S3232 = r2_17 * _S3223 + _S3163;
    float _S3233 = _S3223 + _S3162;
    float _S3234 = _S2946 * _S3212 + _S3151;
    float _S3235 = _S3222 + _S3161;
    float _S3236 = r2_17 * _S3213 + _S3154;
    float _S3237 = _S2943 * _S3213 + _S3156;
    float _S3238 = _S3221 + _S3160;
    float2  _S3239 = _S2924 * _S3228;
    float2  _S3240 = _S2935 * _S3228;
    float _S3241 = _S3239.x + _S3239.y;
    if(_S2928)
    {
        float _S3242 = _S3241 / _S2930;
        float _S3243 = _S2931 * - _S3242;
        float _S3244 = _S2927 * (0.3333333432674408f * - (_S2926 * _S3242));
        k_6 = _S3244 + _S3244;
        _S2929 = _S3243;
        _S2930 = 0.0f;
    }
    else
    {
        float _S3245 = _S3241 / _S2929;
        float _S3246 = _S2927 * - _S3245;
        k_6 = _S2925 * _S3245;
        _S2929 = 0.0f;
        _S2930 = _S3246;
    }
    DiffPair_float_0 _S3247;
    (&_S3247)->primal_0 = _S2925;
    (&_S3247)->differential_0 = 0.0f;
    DiffPair_float_0 _S3248;
    (&_S3248)->primal_0 = _S2926;
    (&_S3248)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3247, &_S3248, k_6);
    float _S3249 = _S3248.differential_0 + _S2929;
    float _S3250 = _S3247.differential_0 + _S2930;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3251;
    (&_S3251)->primal_0 = _S2924;
    (&_S3251)->differential_0 = _S3135;
    s_bwd_length_impl_1(&_S3251, _S3250);
    float2  _S3252 = _S3251.differential_0 + _S3240;
    float _S3253 = fy_24 * (_S3169.y + _S3195);
    float _S3254 = fx_24 * (_S3169.x + _S3194);
    float2  _S3255 = make_float2 (_S3254, _S3253);
    float2  _S3256 = _S2907 * _S3255;
    float _S3257 = p1_6 * _S3253;
    float _S3258 = v_16 * _S3253;
    float _S3259 = p2_6 * _S3254;
    float _S3260 = v_16 * _S3254;
    float _S3261 = _S3256.x + _S3256.y;
    float _S3262 = r2_16 * _S3261;
    float _S3263 = r2_16 * _S3262;
    float _S3264 = r2_16 * _S3263;
    float _S3265 = sy1_6 * _S3253 + _S3257 + sx1_6 * _S3254 + _S3259 + _S2910 * _S3261 + _S2909 * _S3262 + _S2908 * _S3263 + k4_6 * _S3264;
    float _S3266 = v_16 * _S3265;
    float _S3267 = u_16 * _S3265;
    float2  _S3268 = make_float2 (r2_16 * _S3254 + _S3236, r2_16 * _S3253 + _S3231);
    float2  _S3269 = make_float2 (_S2919 * _S3253 + 2.0f * (u_16 * _S3260 + _S3230) + _S3234, 2.0f * (u_16 * _S3258 + _S3229) + _S2915 * _S3254 + _S3237);
    float4  _S3270 = make_float4 (_S3262 + _S3238, _S3263 + _S3235, _S3264 + _S3233, r2_16 * _S3264 + _S3232);
    float3  _S3271 = _S3037.differential_0 + make_float3 (_S3252.x, _S3252.y, _S3249);
    float2  _S3272 = _S2911 * _S3255 + make_float2 (_S2916 * _S3258 + _S2914 * _S3259 + 2.0f * (u_16 * _S3259) + _S2912 * _S3260 + _S3267 + _S3267, _S2918 * _S3257 + 2.0f * (v_16 * _S3257) + _S2917 * _S3253 + _S2913 * _S3254 + _S3266 + _S3266);
    CameraDistortion_0 _S3273 = CameraDistortion_x24_syn_dzero_0();
    (&_S3273)->thin_prism_coeffs_0 = _S3268;
    (&_S3273)->tangential_coeffs_0 = _S3269;
    (&_S3273)->radial_coeffs_0 = _S3270;
    float2  _S3274 = _S2895 * _S3272;
    float2  _S3275 = _S2906 * _S3272;
    float _S3276 = _S3274.x + _S3274.y;
    if(_S2899)
    {
        float _S3277 = _S3276 / _S2901;
        float _S3278 = _S2902 * - _S3277;
        float _S3279 = _S2898 * (0.3333333432674408f * - (_S2897 * _S3277));
        k_6 = _S3279 + _S3279;
        _S2900 = _S3278;
        _S2901 = 0.0f;
    }
    else
    {
        float _S3280 = _S3276 / _S2900;
        float _S3281 = _S2898 * - _S3280;
        k_6 = _S2896 * _S3280;
        _S2900 = 0.0f;
        _S2901 = _S3281;
    }
    DiffPair_float_0 _S3282;
    (&_S3282)->primal_0 = _S2896;
    (&_S3282)->differential_0 = 0.0f;
    DiffPair_float_0 _S3283;
    (&_S3283)->primal_0 = _S2897;
    (&_S3283)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3282, &_S3283, k_6);
    float _S3284 = _S3283.differential_0 + _S2900;
    float _S3285 = _S3282.differential_0 + _S2901;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3286;
    (&_S3286)->primal_0 = _S2895;
    (&_S3286)->differential_0 = _S3135;
    s_bwd_length_impl_1(&_S3286, _S3285);
    float2  _S3287 = _S3286.differential_0 + _S3275;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3288;
    (&_S3288)->primal_0 = vert2_c_5;
    (&_S3288)->differential_0 = _S3029;
    s_bwd_length_impl_0(&_S3288, v_depth_7.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3289;
    (&_S3289)->primal_0 = vert1_c_5;
    (&_S3289)->differential_0 = _S3029;
    s_bwd_length_impl_0(&_S3289, v_depth_7.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3290;
    (&_S3290)->primal_0 = vert0_c_5;
    (&_S3290)->differential_0 = _S3029;
    s_bwd_length_impl_0(&_S3290, v_depth_7.x);
    float3  _S3291 = _S3288.differential_0 + _S3227;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3292;
    (&_S3292)->primal_0 = R_24;
    (&_S3292)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3293;
    (&_S3293)->primal_0 = vert2_1;
    (&_S3293)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3292, &_S3293, _S3291);
    float3  _S3294 = _S3289.differential_0 + _S3271;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3295;
    (&_S3295)->primal_0 = R_24;
    (&_S3295)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3296;
    (&_S3296)->primal_0 = vert1_1;
    (&_S3296)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3295, &_S3296, _S3294);
    float3  _S3297 = _S3290.differential_0 + _S3173 + make_float3 (_S3287.x, _S3287.y, _S3284);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3298;
    (&_S3298)->primal_0 = R_24;
    (&_S3298)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3299;
    (&_S3299)->primal_0 = vert0_1;
    (&_S3299)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3298, &_S3299, _S3297);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3300;
    (&_S3300)->primal_0 = _S2888;
    (&_S3300)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3301;
    (&_S3301)->primal_0 = _S2893;
    (&_S3301)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3300, &_S3301, _S3293.differential_0);
    float _S3302 = - _S3301.differential_0.y;
    float _S3303 = _S2892 * _S3301.differential_0.x;
    float _S3304 = - (_S2884 * _S3301.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3305;
    (&_S3305)->primal_0 = _S2888;
    (&_S3305)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3306;
    (&_S3306)->primal_0 = _S2891;
    (&_S3306)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3305, &_S3306, _S3296.differential_0);
    float _S3307 = _S2884 * _S3306.differential_0.x;
    float _S3308 = _S2890 * _S3306.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3309;
    (&_S3309)->primal_0 = _S2888;
    (&_S3309)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3310;
    (&_S3310)->primal_0 = _S2889;
    (&_S3310)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3309, &_S3310, _S3299.differential_0);
    Matrix<float, 3, 3>  _S3311 = transpose_0(_S3300.differential_0 + _S3305.differential_0 + _S3309.differential_0);
    float _S3312 = 2.0f * - _S3311.rows[int(2)].z;
    float _S3313 = 2.0f * _S3311.rows[int(2)].y;
    float _S3314 = 2.0f * _S3311.rows[int(2)].x;
    float _S3315 = 2.0f * _S3311.rows[int(1)].z;
    float _S3316 = 2.0f * - _S3311.rows[int(1)].y;
    float _S3317 = 2.0f * _S3311.rows[int(1)].x;
    float _S3318 = 2.0f * _S3311.rows[int(0)].z;
    float _S3319 = 2.0f * _S3311.rows[int(0)].y;
    float _S3320 = 2.0f * - _S3311.rows[int(0)].x;
    float _S3321 = - _S3317 + _S3319;
    float _S3322 = _S3314 + - _S3318;
    float _S3323 = - _S3313 + _S3315;
    float _S3324 = _S3313 + _S3315;
    float _S3325 = _S3314 + _S3318;
    float _S3326 = _S3317 + _S3319;
    float _S3327 = quat_25.w * (_S3316 + _S3320);
    float _S3328 = quat_25.z * (_S3312 + _S3320);
    float _S3329 = quat_25.y * (_S3312 + _S3316);
    float _S3330 = quat_25.x * _S3321 + quat_25.z * _S3324 + quat_25.y * _S3325 + _S3327 + _S3327;
    float _S3331 = quat_25.x * _S3322 + quat_25.w * _S3324 + quat_25.y * _S3326 + _S3328 + _S3328;
    float _S3332 = quat_25.x * _S3323 + quat_25.w * _S3325 + quat_25.z * _S3326 + _S3329 + _S3329;
    float _S3333 = quat_25.w * _S3321 + quat_25.z * _S3322 + quat_25.y * _S3323;
    float _S3334 = _S3304 + _S3307;
    float _S3335 = 0.5f * - _S3334;
    float _S3336 = _S3302 + _S3306.differential_0.y;
    DiffPair_float_0 _S3337;
    (&_S3337)->primal_0 = _S2885;
    (&_S3337)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3337, _S3336);
    float _S3338 = _S3335 + _S3337.differential_0;
    float _S3339 = _S3303 + _S3308 + _S3310.differential_0.x;
    DiffPair_float_0 _S3340;
    (&_S3340)->primal_0 = _S2883;
    (&_S3340)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3340, _S3339);
    float _S3341 = _S3335 + _S3340.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3342;
    (&_S3342)->primal_0 = mean_c_20;
    (&_S3342)->differential_0 = _S3029;
    s_bwd_length_impl_0(&_S3342, 0.0f);
    float3  _S3343 = _S3342.differential_0 + _S3032.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3344;
    (&_S3344)->primal_0 = R_24;
    (&_S3344)->differential_0 = _S3112;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3345;
    (&_S3345)->primal_0 = mean_24;
    (&_S3345)->differential_0 = _S3029;
    s_bwd_prop_mul_1(&_S3344, &_S3345, _S3343);
    float3  _S3346 = _S3291 + _S3294 + _S3297 + _S3343 + _S3115.differential_0;
    Matrix<float, 3, 3>  _S3347 = _S3292.differential_0 + _S3295.differential_0 + _S3298.differential_0 + _S3344.differential_0 + _S3116;
    float3  _S3348 = make_float3 (_S3341, _S3338, _S3334);
    float4  _S3349 = make_float4 (0.0f);
    *&((&_S3349)->w) = _S3330;
    *&((&_S3349)->z) = _S3331;
    *&((&_S3349)->y) = _S3332;
    *&((&_S3349)->x) = _S3333;
    float4  _S3350 = _S3349;
    float3  _S3351 = _S3293.differential_0 + _S3296.differential_0 + _S3299.differential_0 + _S3345.differential_0 + _S3110;
    *v_mean_8 = _S3351;
    *v_quat_7 = _S3350;
    *v_scale_7 = _S3348;
    *v_hardness_1 = _S3170;
    (*v_sh_coeffs_6)[int(0)] = _S3178;
    (*v_sh_coeffs_6)[int(1)] = _S3179;
    (*v_sh_coeffs_6)[int(2)] = _S3180;
    (*v_sh_coeffs_6)[int(3)] = _S3181;
    (*v_sh_coeffs_6)[int(4)] = _S3182;
    (*v_sh_coeffs_6)[int(5)] = _S3183;
    (*v_sh_coeffs_6)[int(6)] = _S3184;
    (*v_sh_coeffs_6)[int(7)] = _S3185;
    (*v_sh_coeffs_6)[int(8)] = _S3186;
    (*v_sh_coeffs_6)[int(9)] = _S3187;
    (*v_sh_coeffs_6)[int(10)] = _S3188;
    (*v_sh_coeffs_6)[int(11)] = _S3189;
    (*v_sh_coeffs_6)[int(12)] = _S3190;
    (*v_sh_coeffs_6)[int(13)] = _S3191;
    (*v_sh_coeffs_6)[int(14)] = _S3192;
    (*v_sh_coeffs_6)[int(15)] = _S3193;
    (*v_ch_coeffs_1)[int(0)] = _S3175;
    (*v_ch_coeffs_1)[int(1)] = _S3176;
    *v_R_7 = _S3347;
    *v_t_7 = _S3346;
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_15, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_18)
{
    DiffPair_float_0 _S3352 = *dpx_15;
    bool _S3353;
    if(((*dpx_15).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S3353 = ((*dpx_15).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S3353 = false;
    }
    float _S3354;
    if(_S3353)
    {
        _S3354 = dOut_18;
    }
    else
    {
        _S3354 = 0.0f;
    }
    dpx_15->primal_0 = _S3352.primal_0;
    dpx_15->differential_0 = _S3354;
    DiffPair_float_0 _S3355 = *dpMin_0;
    if((_S3352.primal_0) < ((*dpMin_0).primal_0))
    {
        _S3354 = dOut_18;
    }
    else
    {
        _S3354 = 0.0f;
    }
    dpMin_0->primal_0 = _S3355.primal_0;
    dpMin_0->differential_0 = _S3354;
    DiffPair_float_0 _S3356 = *dpMax_0;
    if(((*dpx_15).primal_0) > ((*dpMax_0).primal_0))
    {
        _S3354 = dOut_18;
    }
    else
    {
        _S3354 = 0.0f;
    }
    dpMax_0->primal_0 = _S3356.primal_0;
    dpMax_0->differential_0 = _S3354;
    return;
}

inline __device__ float clamp_0(float x_50, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_50), (minBound_0)))), (maxBound_0)));
}

inline __device__ void _d_pow_0(DiffPair_float_0 * dpx_16, DiffPair_float_0 * dpy_6, float dOut_19)
{
    if(((*dpx_16).primal_0) < 9.99999997475242708e-07f)
    {
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = 0.0f;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = 0.0f;
    }
    else
    {
        float val_0 = (F32_pow(((*dpx_16).primal_0), ((*dpy_6).primal_0)));
        DiffPair_float_0 _S3357 = *dpx_16;
        float _S3358 = val_0 * (*dpy_6).primal_0 / (*dpx_16).primal_0 * dOut_19;
        dpx_16->primal_0 = (*dpx_16).primal_0;
        dpx_16->differential_0 = _S3358;
        float _S3359 = val_0 * (F32_log((_S3357.primal_0))) * dOut_19;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3359;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_0)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3360 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3360))));
    float2  _S3361 = p_0 - v0_0;
    float2  _S3362 = normalize_1(e0_6);
    float2  _S3363 = p_0 - v1_0;
    float2  _S3364 = normalize_1(e1_6);
    float2  _S3365 = p_0 - v2_0;
    float2  _S3366 = normalize_1(e2_2);
    float _S3367 = hardness_6.x;
    float _S3368 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3361.x * _S3362.y - _S3361.y * _S3362.x)), (se_0 * (_S3363.x * _S3364.y - _S3363.y * _S3364.x))))), (se_0 * (_S3365.x * _S3366.y - _S3365.y * _S3366.x)))) / ((F32_abs((_S3360))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S3368))));
    float _S3369;
    if(a_1 <= 0.0f)
    {
        _S3369 = 0.0f;
    }
    else
    {
        _S3369 = (F32_min(((F32_pow((a_1), (_S3368)))), (0.99900001287460327f)));
    }
    return _S3367 * _S3369;
}

inline __device__ float s_primal_ctx_abs_0(float _S3370)
{
    return (F32_abs((_S3370)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S3371, float _S3372, float _S3373)
{
    return clamp_0(_S3371, _S3372, _S3373);
}

inline __device__ float s_primal_ctx_exp2_0(float _S3374)
{
    return (F32_exp2((_S3374)));
}

inline __device__ float s_primal_ctx_pow_0(float _S3375, float _S3376)
{
    return (F32_pow((_S3375), (_S3376)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S3377, DiffPair_float_0 * _S3378, float _S3379)
{
    _d_pow_0(_S3377, _S3378, _S3379);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S3380, DiffPair_float_0 * _S3381, DiffPair_float_0 * _S3382, float _S3383)
{
    _d_clamp_0(_S3380, _S3381, _S3382, _S3383);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_17, float2  _s_dOut_9)
{
    float _S3384 = length_0((*dpx_17).primal_0);
    float2  _S3385 = (*dpx_17).primal_0 * _s_dOut_9;
    float2  _S3386 = make_float2 (1.0f / _S3384) * _s_dOut_9;
    float _S3387 = - ((_S3385.x + _S3385.y) / (_S3384 * _S3384));
    float2  _S3388 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3389;
    (&_S3389)->primal_0 = (*dpx_17).primal_0;
    (&_S3389)->differential_0 = _S3388;
    s_bwd_length_impl_1(&_S3389, _S3387);
    float2  _S3390 = _S3386 + _S3389.differential_0;
    dpx_17->primal_0 = (*dpx_17).primal_0;
    dpx_17->differential_0 = _S3390;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3391, float2  _S3392)
{
    s_bwd_prop_normalize_impl_1(_S3391, _S3392);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_1, float _s_dOut_10)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S3393 = e0_7.x;
    float _S3394 = e1_7.y;
    float _S3395 = e0_7.y;
    float _S3396 = e1_7.x;
    float _S3397 = _S3393 * _S3394 - _S3395 * _S3396;
    float se_1 = float((F32_sign((_S3397))));
    float2  _S3398 = p_1 - (*dpv0_0).primal_0;
    float2  _S3399 = normalize_1(e0_7);
    float _S3400 = _S3398.x;
    float _S3401 = _S3399.y;
    float _S3402 = _S3398.y;
    float _S3403 = _S3399.x;
    float de0_0 = se_1 * (_S3400 * _S3401 - _S3402 * _S3403);
    float2  _S3404 = p_1 - (*dpv1_0).primal_0;
    float2  _S3405 = normalize_1(e1_7);
    float _S3406 = _S3404.x;
    float _S3407 = _S3405.y;
    float _S3408 = _S3404.y;
    float _S3409 = _S3405.x;
    float de1_0 = se_1 * (_S3406 * _S3407 - _S3408 * _S3409);
    float2  _S3410 = p_1 - (*dpv2_0).primal_0;
    float2  _S3411 = normalize_1(e2_3);
    float _S3412 = _S3410.x;
    float _S3413 = _S3411.y;
    float _S3414 = _S3410.y;
    float _S3415 = _S3411.x;
    float de2_0 = se_1 * (_S3412 * _S3413 - _S3414 * _S3415);
    float _S3416 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S3417 = s_primal_ctx_max_0(_S3416, de2_0);
    float _S3418 = s_primal_ctx_abs_0(_S3397);
    float _S3419 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S3418 / _S3419;
    float _S3420 = _S3419 * _S3419;
    float _S3421 = (*dphardness_0).primal_0.x;
    float _S3422 = (*dphardness_0).primal_0.y;
    float _S3423 = dmax_0 * dmax_0;
    float _S3424 = 1.0f + _S3417 / dmax_0;
    float _S3425 = 1.0f - s_primal_ctx_clamp_0(_S3422, 0.00499999988824129f, 0.98000001907348633f);
    float _S3426 = -1.0f / _S3425;
    float _S3427 = _S3425 * _S3425;
    float _S3428 = 1.0f - s_primal_ctx_exp2_0(_S3426);
    float a_2 = 1.0f - _S3424 * _S3428;
    bool _S3429 = a_2 <= 0.0f;
    float _S3430;
    float _S3431;
    if(_S3429)
    {
        _S3430 = 0.0f;
        _S3431 = 0.0f;
    }
    else
    {
        float _S3432 = s_primal_ctx_pow_0(a_2, _S3425);
        _S3430 = s_primal_ctx_min_0(_S3432, 0.99900001287460327f);
        _S3431 = _S3432;
    }
    float _S3433 = _S3421 * _s_dOut_10;
    float _S3434 = _S3430 * _s_dOut_10;
    if(_S3429)
    {
        _S3430 = 0.0f;
        _S3431 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S3435;
        (&_S3435)->primal_0 = _S3431;
        (&_S3435)->differential_0 = 0.0f;
        DiffPair_float_0 _S3436;
        (&_S3436)->primal_0 = 0.99900001287460327f;
        (&_S3436)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S3435, &_S3436, _S3433);
        DiffPair_float_0 _S3437;
        (&_S3437)->primal_0 = a_2;
        (&_S3437)->differential_0 = 0.0f;
        DiffPair_float_0 _S3438;
        (&_S3438)->primal_0 = _S3425;
        (&_S3438)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S3437, &_S3438, _S3435.differential_0);
        _S3430 = _S3437.differential_0;
        _S3431 = _S3438.differential_0;
    }
    float _S3439 = - _S3430;
    float _S3440 = _S3428 * _S3439;
    float _S3441 = - (_S3424 * _S3439);
    DiffPair_float_0 _S3442;
    (&_S3442)->primal_0 = _S3426;
    (&_S3442)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3442, _S3441);
    float _S3443 = - (-1.0f * - (_S3442.differential_0 / _S3427) + _S3431);
    float _S3444 = _S3440 / _S3423;
    float s_diff_dmax_T_0 = _S3417 * - _S3444;
    float _S3445 = dmax_0 * _S3444;
    DiffPair_float_0 _S3446;
    (&_S3446)->primal_0 = _S3422;
    (&_S3446)->differential_0 = 0.0f;
    DiffPair_float_0 _S3447;
    (&_S3447)->primal_0 = 0.00499999988824129f;
    (&_S3447)->differential_0 = 0.0f;
    DiffPair_float_0 _S3448;
    (&_S3448)->primal_0 = 0.98000001907348633f;
    (&_S3448)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3446, &_S3447, &_S3448, _S3443);
    float _S3449 = s_diff_dmax_T_0 / _S3420;
    float _S3450 = _S3418 * - _S3449;
    float _S3451 = _S3419 * _S3449;
    float2  _S3452 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3453;
    (&_S3453)->primal_0 = e2_3;
    (&_S3453)->differential_0 = _S3452;
    s_bwd_length_impl_1(&_S3453, _S3450);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3454;
    (&_S3454)->primal_0 = e1_7;
    (&_S3454)->differential_0 = _S3452;
    s_bwd_length_impl_1(&_S3454, _S3450);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3455;
    (&_S3455)->primal_0 = e0_7;
    (&_S3455)->differential_0 = _S3452;
    s_bwd_length_impl_1(&_S3455, _S3450);
    DiffPair_float_0 _S3456;
    (&_S3456)->primal_0 = _S3397;
    (&_S3456)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3456, _S3451);
    DiffPair_float_0 _S3457;
    (&_S3457)->primal_0 = _S3416;
    (&_S3457)->differential_0 = 0.0f;
    DiffPair_float_0 _S3458;
    (&_S3458)->primal_0 = de2_0;
    (&_S3458)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3457, &_S3458, _S3445);
    DiffPair_float_0 _S3459;
    (&_S3459)->primal_0 = de0_0;
    (&_S3459)->differential_0 = 0.0f;
    DiffPair_float_0 _S3460;
    (&_S3460)->primal_0 = de1_0;
    (&_S3460)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3459, &_S3460, _S3457.differential_0);
    float _S3461 = se_1 * _S3458.differential_0;
    float _S3462 = - _S3461;
    float _S3463 = _S3415 * _S3462;
    float _S3464 = _S3413 * _S3461;
    float2  _S3465 = make_float2 (_S3414 * _S3462, _S3412 * _S3461);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3466;
    (&_S3466)->primal_0 = e2_3;
    (&_S3466)->differential_0 = _S3452;
    s_bwd_normalize_impl_1(&_S3466, _S3465);
    float2  _S3467 = - make_float2 (_S3464, _S3463);
    float _S3468 = se_1 * _S3460.differential_0;
    float _S3469 = - _S3468;
    float _S3470 = _S3409 * _S3469;
    float _S3471 = _S3407 * _S3468;
    float2  _S3472 = make_float2 (_S3408 * _S3469, _S3406 * _S3468);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3473;
    (&_S3473)->primal_0 = e1_7;
    (&_S3473)->differential_0 = _S3452;
    s_bwd_normalize_impl_1(&_S3473, _S3472);
    float2  _S3474 = - make_float2 (_S3471, _S3470);
    float _S3475 = se_1 * _S3459.differential_0;
    float _S3476 = - _S3475;
    float _S3477 = _S3403 * _S3476;
    float _S3478 = _S3401 * _S3475;
    float2  _S3479 = make_float2 (_S3402 * _S3476, _S3400 * _S3475);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3480;
    (&_S3480)->primal_0 = e0_7;
    (&_S3480)->differential_0 = _S3452;
    s_bwd_normalize_impl_1(&_S3480, _S3479);
    float2  _S3481 = - make_float2 (_S3478, _S3477);
    float _S3482 = - _S3456.differential_0;
    float2  _S3483 = _S3453.differential_0 + _S3466.differential_0;
    float2  _S3484 = - _S3483;
    float2  _S3485 = _S3454.differential_0 + _S3473.differential_0 + make_float2 (_S3395 * _S3482, _S3393 * _S3456.differential_0);
    float2  _S3486 = - _S3485;
    float2  _S3487 = _S3455.differential_0 + _S3480.differential_0 + make_float2 (_S3394 * _S3456.differential_0, _S3396 * _S3482);
    float2  _S3488 = - _S3487;
    float2  _S3489 = make_float2 (_S3434, _S3446.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S3489;
    float2  _S3490 = _S3467 + _S3484 + _S3485;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S3490;
    float2  _S3491 = _S3474 + _S3486 + _S3487;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S3491;
    float2  _S3492 = _S3481 + _S3483 + _S3488;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S3492;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3493, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3494, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3495, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3496, float2  _S3497, float _S3498)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S3493, _S3494, _S3495, _S3496, _S3497, _S3498);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_2, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S3499 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S3499;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S3499;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S3499;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S3499;
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
    float2  _S3500 = p_3 - v0_2;
    float2  _S3501 = p_3 - v1_2;
    float2  _S3502 = p_3 - v2_2;
    float _S3503 = e0_8.x;
    float _S3504 = e1_8.y;
    float _S3505 = e0_8.y;
    float _S3506 = e1_8.x;
    float _S3507 = _S3503 * _S3504 - _S3505 * _S3506;
    float se_2 = float((F32_sign((_S3507))));
    float _S3508 = hardness_8.x;
    float _S3509 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S3500.x * _S3505 - _S3500.y * _S3503)), (se_2 * (_S3501.x * _S3504 - _S3501.y * _S3506))))), (se_2 * (_S3502.x * e2_4.y - _S3502.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S3500 - e0_8 * make_float2 (clamp_0(dot_1(_S3500, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S3501 - e1_8 * make_float2 (clamp_0(dot_1(_S3501, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S3502 - e2_4 * make_float2 (clamp_0(dot_1(_S3502, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S3507))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S3509))));
    float _S3510;
    if(a_3 <= 0.0f)
    {
        _S3510 = 0.0f;
    }
    else
    {
        _S3510 = (F32_min(((F32_pow((a_3), (_S3509)))), (0.99900001287460327f)));
    }
    return _S3508 * _S3510;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S3511, float2  _S3512)
{
    return dot_1(_S3511, _S3512);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3513, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3514, float _S3515)
{
    _d_dot_1(_S3513, _S3514, _S3515);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_4, float _s_dOut_11)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S3516 = p_4 - (*dpv0_1).primal_0;
    float _S3517 = s_primal_ctx_dot_1(_S3516, e0_9);
    float _S3518 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S3519 = _S3517 / _S3518;
    float _S3520 = _S3518 * _S3518;
    float _S3521 = s_primal_ctx_clamp_0(_S3519, 0.0f, 1.0f);
    float2  _S3522 = make_float2 (_S3521);
    float2  _S3523 = _S3516 - e0_9 * make_float2 (_S3521);
    float _S3524 = length_0(_S3523);
    float2  _S3525 = p_4 - (*dpv1_1).primal_0;
    float _S3526 = s_primal_ctx_dot_1(_S3525, e1_9);
    float _S3527 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S3528 = _S3526 / _S3527;
    float _S3529 = _S3527 * _S3527;
    float _S3530 = s_primal_ctx_clamp_0(_S3528, 0.0f, 1.0f);
    float2  _S3531 = make_float2 (_S3530);
    float2  _S3532 = _S3525 - e1_9 * make_float2 (_S3530);
    float _S3533 = length_0(_S3532);
    float2  _S3534 = p_4 - (*dpv2_1).primal_0;
    float _S3535 = s_primal_ctx_dot_1(_S3534, e2_5);
    float _S3536 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S3537 = _S3535 / _S3536;
    float _S3538 = _S3536 * _S3536;
    float _S3539 = s_primal_ctx_clamp_0(_S3537, 0.0f, 1.0f);
    float2  _S3540 = make_float2 (_S3539);
    float2  _S3541 = _S3534 - e2_5 * make_float2 (_S3539);
    float _S3542 = length_0(_S3541);
    float _S3543 = e0_9.x;
    float _S3544 = e1_9.y;
    float _S3545 = e0_9.y;
    float _S3546 = e1_9.x;
    float _S3547 = _S3543 * _S3544 - _S3545 * _S3546;
    float se_3 = float((F32_sign((_S3547))));
    float _S3548 = _S3516.x;
    float _S3549 = _S3516.y;
    float s0_0 = se_3 * (_S3548 * _S3545 - _S3549 * _S3543);
    float _S3550 = _S3525.x;
    float _S3551 = _S3525.y;
    float s1_0 = se_3 * (_S3550 * _S3544 - _S3551 * _S3546);
    float _S3552 = _S3534.x;
    float _S3553 = e2_5.y;
    float _S3554 = _S3534.y;
    float _S3555 = e2_5.x;
    float s2_0 = se_3 * (_S3552 * _S3553 - _S3554 * _S3555);
    float _S3556 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S3556, s2_0)))));
    float _S3557 = s_primal_ctx_min_0(_S3524, _S3533);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S3557, _S3542);
    float _S3558 = s_primal_ctx_abs_0(_S3547);
    float _S3559 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S3558 / _S3559;
    float _S3560 = _S3559 * _S3559;
    float _S3561 = (*dphardness_1).primal_0.x;
    float _S3562 = (*dphardness_1).primal_0.y;
    float _S3563 = dmax_1 * dmax_1;
    float _S3564 = 1.0f + dv_0 / dmax_1;
    float _S3565 = 1.0f - s_primal_ctx_clamp_0(_S3562, 0.00499999988824129f, 0.98000001907348633f);
    float _S3566 = -1.0f / _S3565;
    float _S3567 = _S3565 * _S3565;
    float _S3568 = 1.0f - s_primal_ctx_exp2_0(_S3566);
    float a_4 = 1.0f - _S3564 * _S3568;
    bool _S3569 = a_4 <= 0.0f;
    float _S3570;
    float _S3571;
    if(_S3569)
    {
        _S3570 = 0.0f;
        _S3571 = 0.0f;
    }
    else
    {
        float _S3572 = s_primal_ctx_pow_0(a_4, _S3565);
        _S3570 = s_primal_ctx_min_0(_S3572, 0.99900001287460327f);
        _S3571 = _S3572;
    }
    float _S3573 = _S3561 * _s_dOut_11;
    float _S3574 = _S3570 * _s_dOut_11;
    if(_S3569)
    {
        _S3570 = 0.0f;
        _S3571 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S3575;
        (&_S3575)->primal_0 = _S3571;
        (&_S3575)->differential_0 = 0.0f;
        DiffPair_float_0 _S3576;
        (&_S3576)->primal_0 = 0.99900001287460327f;
        (&_S3576)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S3575, &_S3576, _S3573);
        DiffPair_float_0 _S3577;
        (&_S3577)->primal_0 = a_4;
        (&_S3577)->differential_0 = 0.0f;
        DiffPair_float_0 _S3578;
        (&_S3578)->primal_0 = _S3565;
        (&_S3578)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S3577, &_S3578, _S3575.differential_0);
        _S3570 = _S3577.differential_0;
        _S3571 = _S3578.differential_0;
    }
    float _S3579 = - _S3570;
    float _S3580 = _S3568 * _S3579;
    float _S3581 = - (_S3564 * _S3579);
    DiffPair_float_0 _S3582;
    (&_S3582)->primal_0 = _S3566;
    (&_S3582)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3582, _S3581);
    float _S3583 = - (-1.0f * - (_S3582.differential_0 / _S3567) + _S3571);
    float _S3584 = _S3580 / _S3563;
    float s_diff_dmax_T_1 = dv_0 * - _S3584;
    float s_diff_dv_T_0 = dmax_1 * _S3584;
    DiffPair_float_0 _S3585;
    (&_S3585)->primal_0 = _S3562;
    (&_S3585)->differential_0 = 0.0f;
    DiffPair_float_0 _S3586;
    (&_S3586)->primal_0 = 0.00499999988824129f;
    (&_S3586)->differential_0 = 0.0f;
    DiffPair_float_0 _S3587;
    (&_S3587)->primal_0 = 0.98000001907348633f;
    (&_S3587)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3585, &_S3586, &_S3587, _S3583);
    float _S3588 = s_diff_dmax_T_1 / _S3560;
    float _S3589 = _S3558 * - _S3588;
    float _S3590 = _S3559 * _S3588;
    float2  _S3591 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3592;
    (&_S3592)->primal_0 = e2_5;
    (&_S3592)->differential_0 = _S3591;
    s_bwd_length_impl_1(&_S3592, _S3589);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3593;
    (&_S3593)->primal_0 = e1_9;
    (&_S3593)->differential_0 = _S3591;
    s_bwd_length_impl_1(&_S3593, _S3589);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3594;
    (&_S3594)->primal_0 = e0_9;
    (&_S3594)->differential_0 = _S3591;
    s_bwd_length_impl_1(&_S3594, _S3589);
    DiffPair_float_0 _S3595;
    (&_S3595)->primal_0 = _S3547;
    (&_S3595)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3595, _S3590);
    float _S3596 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S3597;
    (&_S3597)->primal_0 = _S3557;
    (&_S3597)->differential_0 = 0.0f;
    DiffPair_float_0 _S3598;
    (&_S3598)->primal_0 = _S3542;
    (&_S3598)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3597, &_S3598, _S3596);
    DiffPair_float_0 _S3599;
    (&_S3599)->primal_0 = _S3524;
    (&_S3599)->differential_0 = 0.0f;
    DiffPair_float_0 _S3600;
    (&_S3600)->primal_0 = _S3533;
    (&_S3600)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3599, &_S3600, _S3597.differential_0);
    DiffPair_float_0 _S3601;
    (&_S3601)->primal_0 = _S3556;
    (&_S3601)->differential_0 = 0.0f;
    DiffPair_float_0 _S3602;
    (&_S3602)->primal_0 = s2_0;
    (&_S3602)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3601, &_S3602, 0.0f);
    DiffPair_float_0 _S3603;
    (&_S3603)->primal_0 = s0_0;
    (&_S3603)->differential_0 = 0.0f;
    DiffPair_float_0 _S3604;
    (&_S3604)->primal_0 = s1_0;
    (&_S3604)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3603, &_S3604, _S3601.differential_0);
    float _S3605 = se_3 * _S3602.differential_0;
    float _S3606 = - _S3605;
    float _S3607 = _S3554 * _S3606;
    float _S3608 = _S3555 * _S3606;
    float _S3609 = _S3552 * _S3605;
    float _S3610 = _S3553 * _S3605;
    float _S3611 = se_3 * _S3604.differential_0;
    float _S3612 = - _S3611;
    float _S3613 = _S3546 * _S3612;
    float _S3614 = _S3544 * _S3611;
    float _S3615 = se_3 * _S3603.differential_0;
    float _S3616 = - _S3615;
    float _S3617 = _S3543 * _S3616;
    float _S3618 = _S3545 * _S3615;
    float _S3619 = - _S3595.differential_0;
    float _S3620 = _S3551 * _S3612 + _S3545 * _S3619;
    float _S3621 = _S3548 * _S3615 + _S3546 * _S3619;
    float _S3622 = _S3550 * _S3611 + _S3543 * _S3595.differential_0;
    float _S3623 = _S3549 * _S3616 + _S3544 * _S3595.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3624;
    (&_S3624)->primal_0 = _S3541;
    (&_S3624)->differential_0 = _S3591;
    s_bwd_length_impl_1(&_S3624, _S3598.differential_0);
    float2  _S3625 = - _S3624.differential_0;
    float2  _S3626 = e2_5 * _S3625;
    float2  _S3627 = _S3540 * _S3625;
    float _S3628 = _S3626.x + _S3626.y;
    DiffPair_float_0 _S3629;
    (&_S3629)->primal_0 = _S3537;
    (&_S3629)->differential_0 = 0.0f;
    DiffPair_float_0 _S3630;
    (&_S3630)->primal_0 = 0.0f;
    (&_S3630)->differential_0 = 0.0f;
    DiffPair_float_0 _S3631;
    (&_S3631)->primal_0 = 1.0f;
    (&_S3631)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3629, &_S3630, &_S3631, _S3628);
    float _S3632 = _S3629.differential_0 / _S3538;
    float _S3633 = _S3535 * - _S3632;
    float _S3634 = _S3536 * _S3632;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3635;
    (&_S3635)->primal_0 = e2_5;
    (&_S3635)->differential_0 = _S3591;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3636;
    (&_S3636)->primal_0 = e2_5;
    (&_S3636)->differential_0 = _S3591;
    s_bwd_prop_dot_1(&_S3635, &_S3636, _S3633);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3637;
    (&_S3637)->primal_0 = _S3534;
    (&_S3637)->differential_0 = _S3591;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3638;
    (&_S3638)->primal_0 = e2_5;
    (&_S3638)->differential_0 = _S3591;
    s_bwd_prop_dot_1(&_S3637, &_S3638, _S3634);
    float2  _S3639 = - (_S3624.differential_0 + _S3637.differential_0 + make_float2 (_S3610, _S3608));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3640;
    (&_S3640)->primal_0 = _S3532;
    (&_S3640)->differential_0 = _S3591;
    s_bwd_length_impl_1(&_S3640, _S3600.differential_0);
    float2  _S3641 = - _S3640.differential_0;
    float2  _S3642 = e1_9 * _S3641;
    float2  _S3643 = _S3531 * _S3641;
    float _S3644 = _S3642.x + _S3642.y;
    DiffPair_float_0 _S3645;
    (&_S3645)->primal_0 = _S3528;
    (&_S3645)->differential_0 = 0.0f;
    DiffPair_float_0 _S3646;
    (&_S3646)->primal_0 = 0.0f;
    (&_S3646)->differential_0 = 0.0f;
    DiffPair_float_0 _S3647;
    (&_S3647)->primal_0 = 1.0f;
    (&_S3647)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3645, &_S3646, &_S3647, _S3644);
    float _S3648 = _S3645.differential_0 / _S3529;
    float _S3649 = _S3526 * - _S3648;
    float _S3650 = _S3527 * _S3648;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3651;
    (&_S3651)->primal_0 = e1_9;
    (&_S3651)->differential_0 = _S3591;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3652;
    (&_S3652)->primal_0 = e1_9;
    (&_S3652)->differential_0 = _S3591;
    s_bwd_prop_dot_1(&_S3651, &_S3652, _S3649);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3653;
    (&_S3653)->primal_0 = _S3525;
    (&_S3653)->differential_0 = _S3591;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3654;
    (&_S3654)->primal_0 = e1_9;
    (&_S3654)->differential_0 = _S3591;
    s_bwd_prop_dot_1(&_S3653, &_S3654, _S3650);
    float2  _S3655 = - (_S3640.differential_0 + _S3653.differential_0 + make_float2 (_S3614, _S3613));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3656;
    (&_S3656)->primal_0 = _S3523;
    (&_S3656)->differential_0 = _S3591;
    s_bwd_length_impl_1(&_S3656, _S3599.differential_0);
    float2  _S3657 = - _S3656.differential_0;
    float2  _S3658 = e0_9 * _S3657;
    float2  _S3659 = _S3522 * _S3657;
    float _S3660 = _S3658.x + _S3658.y;
    DiffPair_float_0 _S3661;
    (&_S3661)->primal_0 = _S3519;
    (&_S3661)->differential_0 = 0.0f;
    DiffPair_float_0 _S3662;
    (&_S3662)->primal_0 = 0.0f;
    (&_S3662)->differential_0 = 0.0f;
    DiffPair_float_0 _S3663;
    (&_S3663)->primal_0 = 1.0f;
    (&_S3663)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S3661, &_S3662, &_S3663, _S3660);
    float _S3664 = _S3661.differential_0 / _S3520;
    float _S3665 = _S3517 * - _S3664;
    float _S3666 = _S3518 * _S3664;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3667;
    (&_S3667)->primal_0 = e0_9;
    (&_S3667)->differential_0 = _S3591;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3668;
    (&_S3668)->primal_0 = e0_9;
    (&_S3668)->differential_0 = _S3591;
    s_bwd_prop_dot_1(&_S3667, &_S3668, _S3665);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3669;
    (&_S3669)->primal_0 = _S3516;
    (&_S3669)->differential_0 = _S3591;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3670;
    (&_S3670)->primal_0 = e0_9;
    (&_S3670)->differential_0 = _S3591;
    s_bwd_prop_dot_1(&_S3669, &_S3670, _S3666);
    float2  _S3671 = - (_S3656.differential_0 + _S3669.differential_0 + make_float2 (_S3618, _S3617));
    float2  _S3672 = _S3592.differential_0 + _S3627 + _S3636.differential_0 + _S3635.differential_0 + _S3638.differential_0 + make_float2 (_S3607, _S3609);
    float2  _S3673 = - _S3672;
    float2  _S3674 = _S3593.differential_0 + _S3643 + _S3652.differential_0 + _S3651.differential_0 + _S3654.differential_0 + make_float2 (_S3620, _S3622);
    float2  _S3675 = - _S3674;
    float2  _S3676 = _S3594.differential_0 + _S3659 + _S3668.differential_0 + _S3667.differential_0 + _S3670.differential_0 + make_float2 (_S3623, _S3621);
    float2  _S3677 = - _S3676;
    float2  _S3678 = make_float2 (_S3574, _S3585.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S3678;
    float2  _S3679 = _S3639 + _S3673 + _S3674;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S3679;
    float2  _S3680 = _S3655 + _S3675 + _S3676;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S3680;
    float2  _S3681 = _S3671 + _S3672 + _S3677;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S3681;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3682, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3683, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3684, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3685, float2  _S3686, float _S3687)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S3682, _S3683, _S3684, _S3685, _S3686, _S3687);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_5, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S3688 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S3688;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S3688;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S3688;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S3688;
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
    float _S3689 = 0.3333333432674408f * dpdepth_1;
    float3  _S3690 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S3691 = make_float3 (0.0f);
    float3  _S3692 = _S3691;
    *&((&_S3692)->z) = _S3689;
    *&((&_S3692)->y) = _S3689;
    *&((&_S3692)->x) = _S3689;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S3692;
    FixedArray<float3 , 3>  _S3693;
    _S3693[int(0)] = _S3691;
    _S3693[int(1)] = _S3691;
    _S3693[int(2)] = _S3691;
    _S3693[int(2)] = _S3690;
    _S3693[int(1)] = _S3690;
    _S3693[int(0)] = _S3690;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S3693;
    float2  _S3694 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S3694;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S3694;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S3694;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3695, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3696, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3697, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S3698, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3699, float2  _S3700, float3  _S3701, float _S3702)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S3695, _S3696, _S3697, _S3698, _S3699, _S3700, _S3701, _S3702);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_8, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S3703 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S3703;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S3703;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S3703;
    float3  _S3704 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S3705 = { _S3704, _S3704, _S3704 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S3705;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S3704;
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
        float _S3706 = mean_c_21.z;
        bool _S3707;
        if(_S3706 < near_plane_14)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = _S3706 > far_plane_14;
        }
        if(_S3707)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3708 = scale_25.x;
        float sx_4 = (F32_exp((_S3708)));
        float _S3709 = scale_25.y;
        float sy_4 = (F32_exp((_S3709)));
        float sz_6 = scale_25.z - 0.5f * (_S3708 + _S3709);
        float x_51 = quat_26.y;
        float x2_26 = x_51 * x_51;
        float y2_26 = quat_26.z * quat_26.z;
        float z2_47 = quat_26.w * quat_26.w;
        float xy_26 = quat_26.y * quat_26.z;
        float xz_26 = quat_26.y * quat_26.w;
        float yz_26 = quat_26.z * quat_26.w;
        float wx_26 = quat_26.x * quat_26.y;
        float wy_26 = quat_26.x * quat_26.z;
        float wz_26 = quat_26.x * quat_26.w;
        Matrix<float, 3, 3>  _S3710 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_47), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_47), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
        float3  vert0_2 = mul_0(_S3710, make_float3 (sx_4, 0.0f, 0.0f)) + mean_25;
        float3  vert1_2 = mul_0(_S3710, make_float3 (sx_4 * (-0.5f + sz_6), sy_4, 0.0f)) + mean_25;
        float3  vert2_2 = mul_0(_S3710, make_float3 (sx_4 * (-0.5f - sz_6), - sy_4, 0.0f)) + mean_25;
        float3  vert0_c_6 = mul_0(R_25, vert0_2) + t_24;
        float3  vert1_c_6 = mul_0(R_25, vert1_2) + t_24;
        float3  vert2_c_6 = mul_0(R_25, vert2_2) + t_24;
        float _S3711 = vert0_c_6.z;
        float _S3712 = vert1_c_6.z;
        float _S3713 = vert2_c_6.z;
        if(_S3711 < near_plane_14)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = _S3711 > far_plane_14;
        }
        if(_S3707)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = _S3712 < near_plane_14;
        }
        if(_S3707)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = _S3712 > far_plane_14;
        }
        if(_S3707)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = _S3713 < near_plane_14;
        }
        if(_S3707)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = _S3713 > far_plane_14;
        }
        if(_S3707)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  _S3714 = make_float2 (fx_25, fy_25);
        float2  _S3715 = make_float2 (cx_25, cy_25);
        float2  _S3716 = _S3714 * (float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S3711)) + _S3715;
        float2  _S3717 = _S3714 * (float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S3712)) + _S3715;
        float2  _S3718 = _S3714 * (float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S3713)) + _S3715;
        float2  e0_10 = _S3717 - _S3716;
        float2  e1_10 = _S3718 - _S3717;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(_S3716 - _S3718)));
        float _S3719 = _S3716.x;
        float _S3720 = _S3717.x;
        float _S3721 = _S3718.x;
        float xmax_7 = (F32_max(((F32_max((_S3719), (_S3720)))), (_S3721))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S3719), (_S3720)))), (_S3721))) - offset_4;
        float _S3722 = _S3716.y;
        float _S3723 = _S3717.y;
        float _S3724 = _S3718.y;
        float ymax_7 = (F32_max(((F32_max((_S3722), (_S3723)))), (_S3724))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S3722), (_S3723)))), (_S3724))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = xmin_7 >= float(image_width_21);
        }
        if(_S3707)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = ymax_7 <= 0.0f;
        }
        if(_S3707)
        {
            _S3707 = true;
        }
        else
        {
            _S3707 = ymin_7 >= float(image_height_21);
        }
        if(_S3707)
        {
            _S3707 = true;
        }
        else
        {
            if(_S3706 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S3707 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S3707 = false;
                }
                if(_S3707)
                {
                    _S3707 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S3707 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S3707 = false;
                    }
                }
            }
            else
            {
                _S3707 = false;
            }
        }
        if(_S3707)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f;
        float3  _S3725 = mean_25 - - mul_0(transpose_0(R_25), t_24);
        float _S3726 = _S3725.x;
        float _S3727 = _S3725.y;
        float _S3728 = _S3725.z;
        float norm_14 = (F32_sqrt((_S3726 * _S3726 + _S3727 * _S3727 + _S3728 * _S3728)));
        float x_52 = _S3726 / norm_14;
        float y_24 = _S3727 / norm_14;
        float z_21 = _S3728 / norm_14;
        float z2_48 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_52 * x_52 - y_24 * y_24;
        float fS1_21 = 2.0f * x_52 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_48 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_52) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_48 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_52) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_52 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_48 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_52) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_52 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S3729 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S3729);
        float3  _S3730 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S3731 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S3730 + _S3731 + make_float3 (0.5f), _S3729);
        (*rgbs_0)[int(2)] = max_0(_S3730 - _S3731 + make_float3 (0.5f), _S3729);
        (*verts_0)[int(0)] = vert0_2;
        (*verts_0)[int(1)] = vert1_2;
        (*verts_0)[int(2)] = vert2_2;
        float3  _S3732 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S3732 * make_float3 (float(- (F32_sign((dot_0(_S3732, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_26, float3  t_25, float fx_26, float fy_26, float cx_26, float cy_26, float4  radial_coeffs_29, float2  tangential_coeffs_29, float2  thin_prism_coeffs_29, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_26, mean_26) + t_25;
        float _S3733 = length_1(mean_c_22);
        bool _S3734;
        if(_S3733 < near_plane_15)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = _S3733 > far_plane_15;
        }
        if(_S3734)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3735 = scale_26.x;
        float sx_5 = (F32_exp((_S3735)));
        float _S3736 = scale_26.y;
        float sy_5 = (F32_exp((_S3736)));
        float sz_7 = scale_26.z - 0.5f * (_S3735 + _S3736);
        float x_53 = quat_27.y;
        float x2_27 = x_53 * x_53;
        float y2_27 = quat_27.z * quat_27.z;
        float z2_49 = quat_27.w * quat_27.w;
        float xy_27 = quat_27.y * quat_27.z;
        float xz_27 = quat_27.y * quat_27.w;
        float yz_27 = quat_27.z * quat_27.w;
        float wx_27 = quat_27.x * quat_27.y;
        float wy_27 = quat_27.x * quat_27.z;
        float wz_27 = quat_27.x * quat_27.w;
        Matrix<float, 3, 3>  _S3737 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_49), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_49), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S3737, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S3737, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S3737, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_7 = mul_0(R_26, vert0_3) + t_25;
        float3  vert1_c_7 = mul_0(R_26, vert1_3) + t_25;
        float3  vert2_c_7 = mul_0(R_26, vert2_3) + t_25;
        float _S3738 = length_1(vert0_c_7);
        float _S3739 = length_1(vert1_c_7);
        float _S3740 = length_1(vert2_c_7);
        if(_S3738 < near_plane_15)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = _S3738 > far_plane_15;
        }
        if(_S3734)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = _S3739 < near_plane_15;
        }
        if(_S3734)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = _S3739 > far_plane_15;
        }
        if(_S3734)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = _S3740 < near_plane_15;
        }
        if(_S3734)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = _S3740 > far_plane_15;
        }
        if(_S3734)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        CameraDistortion_0 dist_coeffs_4 = CameraDistortion_x24init_0(radial_coeffs_29, tangential_coeffs_29, thin_prism_coeffs_29);
        float2  _S3741 = float2 {vert0_c_7.x, vert0_c_7.y};
        float r_13 = length_0(_S3741);
        float _S3742 = vert0_c_7.z;
        float theta_8 = (F32_atan2((r_13), (_S3742)));
        float k_7;
        if(theta_8 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_8 * theta_8 / 3.0f) / _S3742;
        }
        else
        {
            k_7 = theta_8 / r_13;
        }
        float2  _S3743 = _S3741 * make_float2 (k_7);
        float k1_6 = dist_coeffs_4.radial_coeffs_0.x;
        float k2_6 = dist_coeffs_4.radial_coeffs_0.y;
        float k3_6 = dist_coeffs_4.radial_coeffs_0.z;
        float k4_7 = dist_coeffs_4.radial_coeffs_0.w;
        float p1_7 = dist_coeffs_4.tangential_coeffs_0.x;
        float p2_7 = dist_coeffs_4.tangential_coeffs_0.y;
        float sx1_7 = dist_coeffs_4.thin_prism_coeffs_0.x;
        float sy1_7 = dist_coeffs_4.thin_prism_coeffs_0.y;
        float u_19 = _S3743.x;
        float v_19 = _S3743.y;
        float r2_19 = u_19 * u_19 + v_19 * v_19;
        float _S3744 = 2.0f * p1_7;
        float _S3745 = 2.0f * p2_7;
        float2  _S3746 = _S3743 * make_float2 (1.0f + r2_19 * (k1_6 + r2_19 * (k2_6 + r2_19 * (k3_6 + r2_19 * k4_7)))) + make_float2 (_S3744 * u_19 * v_19 + p2_7 * (r2_19 + 2.0f * u_19 * u_19) + sx1_7 * r2_19, _S3745 * u_19 * v_19 + p1_7 * (r2_19 + 2.0f * v_19 * v_19) + sy1_7 * r2_19);
        float _S3747 = fx_26 * _S3746.x + cx_26;
        float _S3748 = fy_26 * _S3746.y + cy_26;
        float2  _S3749 = make_float2 (_S3747, _S3748);
        float2  _S3750 = float2 {vert1_c_7.x, vert1_c_7.y};
        float r_14 = length_0(_S3750);
        float _S3751 = vert1_c_7.z;
        float theta_9 = (F32_atan2((r_14), (_S3751)));
        if(theta_9 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_9 * theta_9 / 3.0f) / _S3751;
        }
        else
        {
            k_7 = theta_9 / r_14;
        }
        float2  _S3752 = _S3750 * make_float2 (k_7);
        float u_20 = _S3752.x;
        float v_20 = _S3752.y;
        float r2_20 = u_20 * u_20 + v_20 * v_20;
        float2  _S3753 = _S3752 * make_float2 (1.0f + r2_20 * (k1_6 + r2_20 * (k2_6 + r2_20 * (k3_6 + r2_20 * k4_7)))) + make_float2 (_S3744 * u_20 * v_20 + p2_7 * (r2_20 + 2.0f * u_20 * u_20) + sx1_7 * r2_20, _S3745 * u_20 * v_20 + p1_7 * (r2_20 + 2.0f * v_20 * v_20) + sy1_7 * r2_20);
        float _S3754 = fx_26 * _S3753.x + cx_26;
        float _S3755 = fy_26 * _S3753.y + cy_26;
        float2  _S3756 = make_float2 (_S3754, _S3755);
        float2  _S3757 = float2 {vert2_c_7.x, vert2_c_7.y};
        float r_15 = length_0(_S3757);
        float _S3758 = vert2_c_7.z;
        float theta_10 = (F32_atan2((r_15), (_S3758)));
        if(theta_10 < 0.00100000004749745f)
        {
            k_7 = (1.0f - theta_10 * theta_10 / 3.0f) / _S3758;
        }
        else
        {
            k_7 = theta_10 / r_15;
        }
        float2  _S3759 = _S3757 * make_float2 (k_7);
        float u_21 = _S3759.x;
        float v_21 = _S3759.y;
        float r2_21 = u_21 * u_21 + v_21 * v_21;
        float2  _S3760 = _S3759 * make_float2 (1.0f + r2_21 * (k1_6 + r2_21 * (k2_6 + r2_21 * (k3_6 + r2_21 * k4_7)))) + make_float2 (_S3744 * u_21 * v_21 + p2_7 * (r2_21 + 2.0f * u_21 * u_21) + sx1_7 * r2_21, _S3745 * u_21 * v_21 + p1_7 * (r2_21 + 2.0f * v_21 * v_21) + sy1_7 * r2_21);
        float _S3761 = fx_26 * _S3760.x + cx_26;
        float _S3762 = fy_26 * _S3760.y + cy_26;
        float2  _S3763 = make_float2 (_S3761, _S3762);
        float2  e0_11 = _S3756 - _S3749;
        float2  e1_11 = _S3763 - _S3756;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(_S3749 - _S3763)));
        float xmax_8 = (F32_max(((F32_max((_S3747), (_S3754)))), (_S3761))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S3747), (_S3754)))), (_S3761))) - offset_5;
        float ymax_8 = (F32_max(((F32_max((_S3748), (_S3755)))), (_S3762))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S3748), (_S3755)))), (_S3762))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = xmin_8 >= float(image_width_22);
        }
        if(_S3734)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = ymax_8 <= 0.0f;
        }
        if(_S3734)
        {
            _S3734 = true;
        }
        else
        {
            _S3734 = ymin_8 >= float(image_height_22);
        }
        if(_S3734)
        {
            _S3734 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S3734 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S3734 = false;
                }
                if(_S3734)
                {
                    _S3734 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S3734 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S3734 = false;
                    }
                }
            }
            else
            {
                _S3734 = false;
            }
        }
        if(_S3734)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f;
        float3  _S3764 = mean_26 - - mul_0(transpose_0(R_26), t_25);
        float _S3765 = _S3764.x;
        float _S3766 = _S3764.y;
        float _S3767 = _S3764.z;
        float norm_15 = (F32_sqrt((_S3765 * _S3765 + _S3766 * _S3766 + _S3767 * _S3767)));
        float x_54 = _S3765 / norm_15;
        float y_25 = _S3766 / norm_15;
        float z_22 = _S3767 / norm_15;
        float z2_50 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_54 * x_54 - y_25 * y_25;
        float fS1_22 = 2.0f * x_54 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_50 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_54) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_50 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_54) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_54 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_50 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_54) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_54 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S3768 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S3768);
        float3  _S3769 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S3770 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S3769 + _S3770 + make_float3 (0.5f), _S3768);
        (*rgbs_1)[int(2)] = max_0(_S3769 - _S3770 + make_float3 (0.5f), _S3768);
        (*verts_1)[int(0)] = vert0_3;
        (*verts_1)[int(1)] = vert1_3;
        (*verts_1)[int(2)] = vert2_3;
        float3  _S3771 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S3771 * make_float3 (float(- (F32_sign((dot_0(_S3771, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_27, float3  t_26, float fx_27, float fy_27, float cx_27, float cy_27, float4  radial_coeffs_30, float2  tangential_coeffs_30, float2  thin_prism_coeffs_30, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_27, mean_27) + t_26;
    float _S3772 = scale_27.x;
    float sx_6 = (F32_exp((_S3772)));
    float _S3773 = scale_27.y;
    float sy_6 = (F32_exp((_S3773)));
    float sz_8 = scale_27.z - 0.5f * (_S3772 + _S3773);
    float x_55 = quat_28.y;
    float x2_28 = x_55 * x_55;
    float y2_28 = quat_28.z * quat_28.z;
    float z2_51 = quat_28.w * quat_28.w;
    float xy_28 = quat_28.y * quat_28.z;
    float xz_28 = quat_28.y * quat_28.w;
    float yz_28 = quat_28.z * quat_28.w;
    float wx_28 = quat_28.x * quat_28.y;
    float wy_28 = quat_28.x * quat_28.z;
    float wz_28 = quat_28.x * quat_28.w;
    Matrix<float, 3, 3>  _S3774 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_51), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_51), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
    float3  vert0_4 = mul_0(_S3774, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
    float3  vert1_4 = mul_0(_S3774, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
    float3  vert2_4 = mul_0(_S3774, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
    float3  vert0_c_8 = mul_0(R_27, vert0_4) + t_26;
    float3  vert1_c_8 = mul_0(R_27, vert1_4) + t_26;
    float3  vert2_c_8 = mul_0(R_27, vert2_4) + t_26;
    float2  _S3775 = make_float2 (fx_27, fy_27);
    float2  _S3776 = make_float2 (cx_27, cy_27);
    float2  _S3777 = _S3775 * (float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z)) + _S3776;
    float2  _S3778 = _S3775 * (float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z)) + _S3776;
    float2  _S3779 = _S3775 * (float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z)) + _S3776;
    float2  e0_12 = _S3778 - _S3777;
    float2  e1_12 = _S3779 - _S3778;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(_S3777 - _S3779)));
    float _S3780 = _S3777.x;
    float _S3781 = _S3778.x;
    float _S3782 = _S3779.x;
    float _S3783 = _S3777.y;
    float _S3784 = _S3778.y;
    float _S3785 = _S3779.y;
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S3780), (_S3781)))), (_S3782))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S3783), (_S3784)))), (_S3785))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S3780), (_S3781)))), (_S3782))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S3783), (_S3784)))), (_S3785))) + offset_6)))));
    *depth_18 = length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f;
    float3  _S3786 = mean_27 - - mul_0(transpose_0(R_27), t_26);
    float _S3787 = _S3786.x;
    float _S3788 = _S3786.y;
    float _S3789 = _S3786.z;
    float norm_16 = (F32_sqrt((_S3787 * _S3787 + _S3788 * _S3788 + _S3789 * _S3789)));
    float x_56 = _S3787 / norm_16;
    float y_26 = _S3788 / norm_16;
    float z_23 = _S3789 / norm_16;
    float z2_52 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_56 * x_56 - y_26 * y_26;
    float fS1_23 = 2.0f * x_56 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_52 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_56) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_52 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_56) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_56 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_52 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_56) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_56 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S3790 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S3790);
    float3  _S3791 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S3792 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S3791 + _S3792 + make_float3 (0.5f), _S3790);
    (*rgbs_2)[int(2)] = max_0(_S3791 - _S3792 + make_float3 (0.5f), _S3790);
    (*verts_2)[int(0)] = vert0_4;
    (*verts_2)[int(1)] = vert1_4;
    (*verts_2)[int(2)] = vert2_4;
    float3  _S3793 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S3793 * make_float3 (float(- (F32_sign((dot_0(_S3793, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_28, float3  t_27, float fx_28, float fy_28, float cx_28, float cy_28, float4  radial_coeffs_31, float2  tangential_coeffs_31, float2  thin_prism_coeffs_31, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_28, mean_28) + t_27;
    float _S3794 = scale_28.x;
    float sx_7 = (F32_exp((_S3794)));
    float _S3795 = scale_28.y;
    float sy_7 = (F32_exp((_S3795)));
    float sz_9 = scale_28.z - 0.5f * (_S3794 + _S3795);
    float x_57 = quat_29.y;
    float x2_29 = x_57 * x_57;
    float y2_29 = quat_29.z * quat_29.z;
    float z2_53 = quat_29.w * quat_29.w;
    float xy_29 = quat_29.y * quat_29.z;
    float xz_29 = quat_29.y * quat_29.w;
    float yz_29 = quat_29.z * quat_29.w;
    float wx_29 = quat_29.x * quat_29.y;
    float wy_29 = quat_29.x * quat_29.z;
    float wz_29 = quat_29.x * quat_29.w;
    Matrix<float, 3, 3>  _S3796 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_53), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_53), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S3796, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S3796, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S3796, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_9 = mul_0(R_28, vert0_5) + t_27;
    float3  vert1_c_9 = mul_0(R_28, vert1_5) + t_27;
    float3  vert2_c_9 = mul_0(R_28, vert2_5) + t_27;
    CameraDistortion_0 dist_coeffs_5 = CameraDistortion_x24init_0(radial_coeffs_31, tangential_coeffs_31, thin_prism_coeffs_31);
    float2  _S3797 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_16 = length_0(_S3797);
    float _S3798 = vert0_c_9.z;
    float theta_11 = (F32_atan2((r_16), (_S3798)));
    float k_8;
    if(theta_11 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_11 * theta_11 / 3.0f) / _S3798;
    }
    else
    {
        k_8 = theta_11 / r_16;
    }
    float2  _S3799 = _S3797 * make_float2 (k_8);
    float k1_7 = dist_coeffs_5.radial_coeffs_0.x;
    float k2_7 = dist_coeffs_5.radial_coeffs_0.y;
    float k3_7 = dist_coeffs_5.radial_coeffs_0.z;
    float k4_8 = dist_coeffs_5.radial_coeffs_0.w;
    float p1_8 = dist_coeffs_5.tangential_coeffs_0.x;
    float p2_8 = dist_coeffs_5.tangential_coeffs_0.y;
    float sx1_8 = dist_coeffs_5.thin_prism_coeffs_0.x;
    float sy1_8 = dist_coeffs_5.thin_prism_coeffs_0.y;
    float u_22 = _S3799.x;
    float v_22 = _S3799.y;
    float r2_22 = u_22 * u_22 + v_22 * v_22;
    float _S3800 = 2.0f * p1_8;
    float _S3801 = 2.0f * p2_8;
    float2  _S3802 = _S3799 * make_float2 (1.0f + r2_22 * (k1_7 + r2_22 * (k2_7 + r2_22 * (k3_7 + r2_22 * k4_8)))) + make_float2 (_S3800 * u_22 * v_22 + p2_8 * (r2_22 + 2.0f * u_22 * u_22) + sx1_8 * r2_22, _S3801 * u_22 * v_22 + p1_8 * (r2_22 + 2.0f * v_22 * v_22) + sy1_8 * r2_22);
    float _S3803 = fx_28 * _S3802.x + cx_28;
    float _S3804 = fy_28 * _S3802.y + cy_28;
    float2  _S3805 = make_float2 (_S3803, _S3804);
    float2  _S3806 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_17 = length_0(_S3806);
    float _S3807 = vert1_c_9.z;
    float theta_12 = (F32_atan2((r_17), (_S3807)));
    if(theta_12 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_12 * theta_12 / 3.0f) / _S3807;
    }
    else
    {
        k_8 = theta_12 / r_17;
    }
    float2  _S3808 = _S3806 * make_float2 (k_8);
    float u_23 = _S3808.x;
    float v_23 = _S3808.y;
    float r2_23 = u_23 * u_23 + v_23 * v_23;
    float2  _S3809 = _S3808 * make_float2 (1.0f + r2_23 * (k1_7 + r2_23 * (k2_7 + r2_23 * (k3_7 + r2_23 * k4_8)))) + make_float2 (_S3800 * u_23 * v_23 + p2_8 * (r2_23 + 2.0f * u_23 * u_23) + sx1_8 * r2_23, _S3801 * u_23 * v_23 + p1_8 * (r2_23 + 2.0f * v_23 * v_23) + sy1_8 * r2_23);
    float _S3810 = fx_28 * _S3809.x + cx_28;
    float _S3811 = fy_28 * _S3809.y + cy_28;
    float2  _S3812 = make_float2 (_S3810, _S3811);
    float2  _S3813 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_18 = length_0(_S3813);
    float _S3814 = vert2_c_9.z;
    float theta_13 = (F32_atan2((r_18), (_S3814)));
    if(theta_13 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_13 * theta_13 / 3.0f) / _S3814;
    }
    else
    {
        k_8 = theta_13 / r_18;
    }
    float2  _S3815 = _S3813 * make_float2 (k_8);
    float u_24 = _S3815.x;
    float v_24 = _S3815.y;
    float r2_24 = u_24 * u_24 + v_24 * v_24;
    float2  _S3816 = _S3815 * make_float2 (1.0f + r2_24 * (k1_7 + r2_24 * (k2_7 + r2_24 * (k3_7 + r2_24 * k4_8)))) + make_float2 (_S3800 * u_24 * v_24 + p2_8 * (r2_24 + 2.0f * u_24 * u_24) + sx1_8 * r2_24, _S3801 * u_24 * v_24 + p1_8 * (r2_24 + 2.0f * v_24 * v_24) + sy1_8 * r2_24);
    float _S3817 = fx_28 * _S3816.x + cx_28;
    float _S3818 = fy_28 * _S3816.y + cy_28;
    float2  _S3819 = make_float2 (_S3817, _S3818);
    float2  e0_13 = _S3812 - _S3805;
    float2  e1_13 = _S3819 - _S3812;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(_S3805 - _S3819)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S3803), (_S3810)))), (_S3817))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S3804), (_S3811)))), (_S3818))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S3803), (_S3810)))), (_S3817))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S3804), (_S3811)))), (_S3818))) + offset_7)))));
    *depth_19 = length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f;
    float3  _S3820 = mean_28 - - mul_0(transpose_0(R_28), t_27);
    float _S3821 = _S3820.x;
    float _S3822 = _S3820.y;
    float _S3823 = _S3820.z;
    float norm_17 = (F32_sqrt((_S3821 * _S3821 + _S3822 * _S3822 + _S3823 * _S3823)));
    float x_58 = _S3821 / norm_17;
    float y_27 = _S3822 / norm_17;
    float z_24 = _S3823 / norm_17;
    float z2_54 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_58 * x_58 - y_27 * y_27;
    float fS1_24 = 2.0f * x_58 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_54 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_58) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_54 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_58) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_58 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_54 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_58) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_58 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S3824 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S3824);
    float3  _S3825 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S3826 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S3825 + _S3826 + make_float3 (0.5f), _S3824);
    (*rgbs_3)[int(2)] = max_0(_S3825 - _S3826 + make_float3 (0.5f), _S3824);
    (*verts_3)[int(0)] = vert0_5;
    (*verts_3)[int(1)] = vert1_5;
    (*verts_3)[int(2)] = vert2_5;
    float3  _S3827 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S3827 * make_float3 (float(- (F32_sign((dot_0(_S3827, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_29, float3  t_28, float fx_29, float fy_29, float cx_29, float cy_29, float4  radial_coeffs_32, float2  tangential_coeffs_32, float2  thin_prism_coeffs_32, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_1(R_29, mean_29) + t_28;
    float _S3828 = scale_29.x;
    float _S3829 = s_primal_ctx_exp_1(_S3828);
    float _S3830 = scale_29.y;
    float _S3831 = s_primal_ctx_exp_1(_S3830);
    float sz_10 = scale_29.z - 0.5f * (_S3828 + _S3830);
    float _S3832 = quat_30.y;
    float x2_30 = _S3832 * _S3832;
    float y2_30 = quat_30.z * quat_30.z;
    float z2_55 = quat_30.w * quat_30.w;
    float xy_30 = quat_30.y * quat_30.z;
    float xz_30 = quat_30.y * quat_30.w;
    float yz_30 = quat_30.z * quat_30.w;
    float wx_30 = quat_30.x * quat_30.y;
    float wy_30 = quat_30.x * quat_30.z;
    float wz_30 = quat_30.x * quat_30.w;
    Matrix<float, 3, 3>  _S3833 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_55), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_55), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  _S3834 = make_float3 (_S3829, 0.0f, 0.0f);
    float3  vert0_6 = s_primal_ctx_mul_1(_S3833, _S3834) + mean_29;
    float _S3835 = -0.5f + sz_10;
    float3  _S3836 = make_float3 (_S3829 * _S3835, _S3831, 0.0f);
    float3  vert1_6 = s_primal_ctx_mul_1(_S3833, _S3836) + mean_29;
    float _S3837 = -0.5f - sz_10;
    float3  _S3838 = make_float3 (_S3829 * _S3837, - _S3831, 0.0f);
    float3  vert2_6 = s_primal_ctx_mul_1(_S3833, _S3838) + mean_29;
    float3  vert0_c_10 = s_primal_ctx_mul_1(R_29, vert0_6) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_1(R_29, vert1_6) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_1(R_29, vert2_6) + t_28;
    float2  _S3839 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S3840 = vert0_c_10.z;
    float2  _S3841 = make_float2 (_S3840);
    float2  _S3842 = make_float2 (_S3840 * _S3840);
    float2  _S3843 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S3844 = vert1_c_10.z;
    float2  _S3845 = make_float2 (_S3844);
    float2  _S3846 = make_float2 (_S3844 * _S3844);
    float2  _S3847 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S3848 = vert2_c_10.z;
    float2  _S3849 = make_float2 (_S3848);
    float2  _S3850 = make_float2 (_S3848 * _S3848);
    float2  _S3851 = make_float2 (fx_29, fy_29);
    float2  _S3852 = make_float2 (cx_29, cy_29);
    float2  _S3853 = _S3851 * (_S3839 / make_float2 (_S3840)) + _S3852;
    float2  _S3854 = _S3851 * (_S3843 / make_float2 (_S3844)) + _S3852;
    float2  _S3855 = _S3851 * (_S3847 / make_float2 (_S3848)) + _S3852;
    float2  e0_14 = _S3854 - _S3853;
    float2  e1_14 = _S3855 - _S3854;
    float2  e2_6 = _S3853 - _S3855;
    float _S3856 = e0_14.x;
    float _S3857 = e1_14.y;
    float _S3858 = e0_14.y;
    float _S3859 = e1_14.x;
    float _S3860 = _S3856 * _S3857 - _S3858 * _S3859;
    float _S3861 = 1.0f - hardness_14.y;
    float _S3862 = -1.0f / _S3861;
    float _S3863 = _S3861 * _S3861;
    float _S3864 = _S3853.x;
    float _S3865 = _S3854.x;
    float _S3866 = s_primal_ctx_max_0(_S3864, _S3865);
    float _S3867 = _S3855.x;
    float _S3868 = s_primal_ctx_min_0(_S3864, _S3865);
    float _S3869 = _S3853.y;
    float _S3870 = _S3854.y;
    float _S3871 = s_primal_ctx_max_0(_S3869, _S3870);
    float _S3872 = _S3855.y;
    float _S3873 = s_primal_ctx_min_0(_S3869, _S3870);
    float3  _S3874 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    Matrix<float, 3, 3>  _S3875 = transpose_0(R_29);
    float3  _S3876 = mean_29 - - s_primal_ctx_mul_1(_S3875, t_28);
    float _S3877 = _S3876.x;
    float _S3878 = _S3876.y;
    float _S3879 = _S3876.z;
    float _S3880 = _S3877 * _S3877 + _S3878 * _S3878 + _S3879 * _S3879;
    float _S3881 = s_primal_ctx_sqrt_0(_S3880);
    float x_59 = _S3877 / _S3881;
    float3  _S3882 = make_float3 (x_59);
    float _S3883 = _S3881 * _S3881;
    float y_28 = _S3878 / _S3881;
    float z_25 = _S3879 / _S3881;
    float3  _S3884 = make_float3 (z_25);
    float _S3885 = - y_28;
    float3  _S3886 = make_float3 (_S3885);
    float z2_56 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_59 * x_59 - y_28 * y_28;
    float _S3887 = 2.0f * x_59;
    float fS1_25 = _S3887 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_56 - 0.31539157032966614f;
    float3  _S3888 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_59;
    float3  _S3889 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S3890 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S3891 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S3892 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_56 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S3893 = 1.86588168144226074f * z2_56 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S3893;
    float3  _S3894 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_59;
    float3  _S3895 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S3896 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S3897 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S3898 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_59 * fC1_25 - y_28 * fS1_25);
    float3  _S3899 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_59 * fS1_25 + y_28 * fC1_25);
    float3  _S3900 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3885) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_59) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S3901 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S3902 = make_float3 (0.0f);
    float3  _S3903 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S3904 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3905 = make_float3 (_S3904);
    float3  _S3906 = (*ch_coeffs_10)[int(1)] * make_float3 (_S3904);
    float3  _S3907 = _S3903 + _S3906 + make_float3 (0.5f);
    float3  _S3908 = _S3903 - _S3906 + make_float3 (0.5f);
    float3  _S3909 = vert1_c_10 - vert0_c_10;
    float3  _S3910 = vert2_c_10 - vert0_c_10;
    float3  _S3911 = s_primal_ctx_cross_0(_S3909, _S3910);
    float3  _S3912 = normalize_0(_S3911);
    float3  _S3913 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3912, mean_c_25)))))) * v_normal_2;
    float3  _S3914 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3915;
    (&_S3915)->primal_0 = _S3912;
    (&_S3915)->differential_0 = _S3914;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3916;
    (&_S3916)->primal_0 = mean_c_25;
    (&_S3916)->differential_0 = _S3914;
    s_bwd_prop_dot_0(&_S3915, &_S3916, 0.0f);
    float3  _S3917 = _S3913 + _S3915.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3918;
    (&_S3918)->primal_0 = _S3911;
    (&_S3918)->differential_0 = _S3914;
    s_bwd_normalize_impl_0(&_S3918, _S3917);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3919;
    (&_S3919)->primal_0 = _S3909;
    (&_S3919)->differential_0 = _S3914;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3920;
    (&_S3920)->primal_0 = _S3910;
    (&_S3920)->differential_0 = _S3914;
    s_bwd_prop_cross_0(&_S3919, &_S3920, _S3918.differential_0);
    float3  _S3921 = - _S3920.differential_0;
    float3  _S3922 = - _S3919.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3923;
    (&_S3923)->primal_0 = _S3908;
    (&_S3923)->differential_0 = _S3914;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3924;
    (&_S3924)->primal_0 = _S3902;
    (&_S3924)->differential_0 = _S3914;
    s_bwd_prop_max_0(&_S3923, &_S3924, (*v_rgbs_0)[int(2)]);
    float3  _S3925 = - _S3923.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3926;
    (&_S3926)->primal_0 = _S3907;
    (&_S3926)->differential_0 = _S3914;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3927;
    (&_S3927)->primal_0 = _S3902;
    (&_S3927)->differential_0 = _S3914;
    s_bwd_prop_max_0(&_S3926, &_S3927, (*v_rgbs_0)[int(1)]);
    float3  _S3928 = _S3905 * (_S3925 + _S3926.differential_0);
    float3  _S3929 = _S3923.differential_0 + _S3926.differential_0;
    float3  _S3930 = make_float3 (0.5f) * - _S3929;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3931;
    (&_S3931)->primal_0 = _S3901;
    (&_S3931)->differential_0 = _S3914;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3932;
    (&_S3932)->primal_0 = _S3902;
    (&_S3932)->differential_0 = _S3914;
    s_bwd_prop_max_0(&_S3931, &_S3932, (*v_rgbs_0)[int(0)]);
    float3  _S3933 = _S3930 + _S3931.differential_0;
    float3  _S3934 = _S3929 + _S3931.differential_0;
    float3  _S3935 = _S3899 * _S3934;
    float3  _S3936 = (*sh_coeffs_25)[int(15)] * _S3934;
    float3  _S3937 = _S3897 * _S3934;
    float3  _S3938 = (*sh_coeffs_25)[int(14)] * _S3934;
    float3  _S3939 = _S3895 * _S3934;
    float3  _S3940 = (*sh_coeffs_25)[int(13)] * _S3934;
    float3  _S3941 = _S3894 * _S3934;
    float3  _S3942 = (*sh_coeffs_25)[int(12)] * _S3934;
    float3  _S3943 = _S3896 * _S3934;
    float3  _S3944 = (*sh_coeffs_25)[int(11)] * _S3934;
    float3  _S3945 = _S3898 * _S3934;
    float3  _S3946 = (*sh_coeffs_25)[int(10)] * _S3934;
    float3  _S3947 = _S3900 * _S3934;
    float3  _S3948 = (*sh_coeffs_25)[int(9)] * _S3934;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S3948.x + _S3948.y + _S3948.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S3936.x + _S3936.y + _S3936.z);
    float _S3949 = _S3946.x + _S3946.y + _S3946.z;
    float _S3950 = _S3938.x + _S3938.y + _S3938.z;
    float _S3951 = _S3944.x + _S3944.y + _S3944.z;
    float _S3952 = _S3940.x + _S3940.y + _S3940.z;
    float _S3953 = _S3942.x + _S3942.y + _S3942.z;
    float _S3954 = - s_diff_fC2_T_7;
    float3  _S3955 = _S3891 * _S3934;
    float3  _S3956 = (*sh_coeffs_25)[int(8)] * _S3934;
    float3  _S3957 = _S3889 * _S3934;
    float3  _S3958 = (*sh_coeffs_25)[int(7)] * _S3934;
    float3  _S3959 = _S3888 * _S3934;
    float3  _S3960 = (*sh_coeffs_25)[int(6)] * _S3934;
    float3  _S3961 = _S3890 * _S3934;
    float3  _S3962 = (*sh_coeffs_25)[int(5)] * _S3934;
    float3  _S3963 = _S3892 * _S3934;
    float3  _S3964 = (*sh_coeffs_25)[int(4)] * _S3934;
    float _S3965 = _S3962.x + _S3962.y + _S3962.z;
    float _S3966 = _S3958.x + _S3958.y + _S3958.z;
    float _S3967 = fTmp1B_25 * _S3949 + x_59 * s_diff_fS2_T_7 + y_28 * _S3954 + 0.54627424478530884f * (_S3964.x + _S3964.y + _S3964.z);
    float _S3968 = fTmp1B_25 * _S3950 + y_28 * s_diff_fS2_T_7 + x_59 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S3956.x + _S3956.y + _S3956.z);
    float _S3969 = y_28 * - _S3968;
    float _S3970 = x_59 * _S3968;
    float _S3971 = z_25 * (1.86588168144226074f * (z_25 * _S3953) + -2.28522896766662598f * (y_28 * _S3951 + x_59 * _S3952) + 0.94617468118667603f * (_S3960.x + _S3960.y + _S3960.z));
    float3  _S3972 = make_float3 (0.48860251903533936f) * _S3934;
    float3  _S3973 = - _S3972;
    float3  _S3974 = _S3882 * _S3973;
    float3  _S3975 = (*sh_coeffs_25)[int(3)] * _S3973;
    float3  _S3976 = _S3884 * _S3972;
    float3  _S3977 = (*sh_coeffs_25)[int(2)] * _S3972;
    float3  _S3978 = _S3886 * _S3972;
    float3  _S3979 = (*sh_coeffs_25)[int(1)] * _S3972;
    float _S3980 = (_S3893 * _S3953 + 1.44530570507049561f * (fS1_25 * _S3949 + fC1_25 * _S3950) + -1.09254848957061768f * (y_28 * _S3965 + x_59 * _S3966) + _S3971 + _S3971 + _S3977.x + _S3977.y + _S3977.z) / _S3883;
    float _S3981 = _S3881 * _S3980;
    float _S3982 = (fTmp0C_25 * _S3951 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S3954 + fTmp0B_25 * _S3965 + _S3887 * _S3967 + _S3969 + _S3969 + - (_S3979.x + _S3979.y + _S3979.z)) / _S3883;
    float _S3983 = _S3881 * _S3982;
    float _S3984 = (fTmp0C_25 * _S3952 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S3966 + 2.0f * (y_28 * _S3967) + _S3970 + _S3970 + _S3975.x + _S3975.y + _S3975.z) / _S3883;
    float _S3985 = _S3881 * _S3984;
    float _S3986 = _S3879 * - _S3980 + _S3878 * - _S3982 + _S3877 * - _S3984;
    DiffPair_float_0 _S3987;
    (&_S3987)->primal_0 = _S3880;
    (&_S3987)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3987, _S3986);
    float _S3988 = _S3879 * _S3987.differential_0;
    float _S3989 = _S3878 * _S3987.differential_0;
    float _S3990 = _S3877 * _S3987.differential_0;
    float3  _S3991 = make_float3 (0.282094806432724f) * _S3934;
    float3  _S3992 = make_float3 (_S3985 + _S3990 + _S3990, _S3983 + _S3989 + _S3989, _S3981 + _S3988 + _S3988);
    float3  _S3993 = - - _S3992;
    Matrix<float, 3, 3>  _S3994 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3995;
    (&_S3995)->primal_0 = _S3875;
    (&_S3995)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3996;
    (&_S3996)->primal_0 = t_28;
    (&_S3996)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S3995, &_S3996, _S3993);
    Matrix<float, 3, 3>  _S3997 = transpose_0(_S3995.differential_0);
    float _S3998 = 0.3333333432674408f * v_depth_9;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3999;
    (&_S3999)->primal_0 = _S3874;
    (&_S3999)->differential_0 = _S3914;
    s_bwd_length_impl_0(&_S3999, _S3998);
    DiffPair_float_0 _S4000;
    (&_S4000)->primal_0 = _S3873;
    (&_S4000)->differential_0 = 0.0f;
    DiffPair_float_0 _S4001;
    (&_S4001)->primal_0 = _S3872;
    (&_S4001)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4000, &_S4001, 0.0f);
    DiffPair_float_0 _S4002;
    (&_S4002)->primal_0 = _S3869;
    (&_S4002)->differential_0 = 0.0f;
    DiffPair_float_0 _S4003;
    (&_S4003)->primal_0 = _S3870;
    (&_S4003)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4002, &_S4003, _S4000.differential_0);
    DiffPair_float_0 _S4004;
    (&_S4004)->primal_0 = _S3871;
    (&_S4004)->differential_0 = 0.0f;
    DiffPair_float_0 _S4005;
    (&_S4005)->primal_0 = _S3872;
    (&_S4005)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4004, &_S4005, 0.0f);
    float _S4006 = _S4001.differential_0 + _S4005.differential_0;
    DiffPair_float_0 _S4007;
    (&_S4007)->primal_0 = _S3869;
    (&_S4007)->differential_0 = 0.0f;
    DiffPair_float_0 _S4008;
    (&_S4008)->primal_0 = _S3870;
    (&_S4008)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4007, &_S4008, _S4004.differential_0);
    float _S4009 = _S4003.differential_0 + _S4008.differential_0;
    float _S4010 = _S4002.differential_0 + _S4007.differential_0;
    DiffPair_float_0 _S4011;
    (&_S4011)->primal_0 = _S3868;
    (&_S4011)->differential_0 = 0.0f;
    DiffPair_float_0 _S4012;
    (&_S4012)->primal_0 = _S3867;
    (&_S4012)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4011, &_S4012, 0.0f);
    DiffPair_float_0 _S4013;
    (&_S4013)->primal_0 = _S3864;
    (&_S4013)->differential_0 = 0.0f;
    DiffPair_float_0 _S4014;
    (&_S4014)->primal_0 = _S3865;
    (&_S4014)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4013, &_S4014, _S4011.differential_0);
    DiffPair_float_0 _S4015;
    (&_S4015)->primal_0 = _S3866;
    (&_S4015)->differential_0 = 0.0f;
    DiffPair_float_0 _S4016;
    (&_S4016)->primal_0 = _S3867;
    (&_S4016)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4015, &_S4016, 0.0f);
    float _S4017 = _S4012.differential_0 + _S4016.differential_0;
    DiffPair_float_0 _S4018;
    (&_S4018)->primal_0 = _S3864;
    (&_S4018)->differential_0 = 0.0f;
    DiffPair_float_0 _S4019;
    (&_S4019)->primal_0 = _S3865;
    (&_S4019)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4018, &_S4019, _S4015.differential_0);
    float _S4020 = _S4014.differential_0 + _S4019.differential_0;
    float _S4021 = _S4013.differential_0 + _S4018.differential_0;
    DiffPair_float_0 _S4022;
    (&_S4022)->primal_0 = _S3862;
    (&_S4022)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4022, 0.0f);
    float _S4023 = - (-1.0f * - (_S4022.differential_0 / _S3863));
    float2  _S4024 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4025;
    (&_S4025)->primal_0 = e2_6;
    (&_S4025)->differential_0 = _S4024;
    s_bwd_length_impl_1(&_S4025, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4026;
    (&_S4026)->primal_0 = e1_14;
    (&_S4026)->differential_0 = _S4024;
    s_bwd_length_impl_1(&_S4026, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4027;
    (&_S4027)->primal_0 = e0_14;
    (&_S4027)->differential_0 = _S4024;
    s_bwd_length_impl_1(&_S4027, -0.0f);
    DiffPair_float_0 _S4028;
    (&_S4028)->primal_0 = _S3860;
    (&_S4028)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4028, 0.0f);
    float _S4029 = - _S4028.differential_0;
    float2  _S4030 = _S4026.differential_0 + make_float2 (_S3858 * _S4029, _S3856 * _S4028.differential_0);
    float2  _S4031 = _S4027.differential_0 + make_float2 (_S3857 * _S4028.differential_0, _S3859 * _S4029);
    float2  _S4032 = _S3851 * (- _S4025.differential_0 + _S4030 + make_float2 (_S4017, _S4006)) / _S3850;
    float2  _S4033 = _S3847 * - _S4032;
    float2  _S4034 = _S3849 * _S4032;
    float2  _S4035 = _S3851 * (- _S4030 + _S4031 + make_float2 (_S4020, _S4009)) / _S3846;
    float2  _S4036 = _S3843 * - _S4035;
    float2  _S4037 = _S3845 * _S4035;
    float _S4038 = _S4036.x + _S4036.y;
    float2  _S4039 = _S3851 * (_S4025.differential_0 + - _S4031 + make_float2 (_S4021, _S4010)) / _S3842;
    float2  _S4040 = _S3839 * - _S4039;
    float2  _S4041 = _S3841 * _S4039;
    float _S4042 = _S4040.x + _S4040.y;
    float3  _S4043 = _S3920.differential_0 + _S3999.differential_0 + make_float3 (_S4034.x, _S4034.y, _S4033.x + _S4033.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4044;
    (&_S4044)->primal_0 = R_29;
    (&_S4044)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4045;
    (&_S4045)->primal_0 = vert2_6;
    (&_S4045)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S4044, &_S4045, _S4043);
    float3  _S4046 = _S3919.differential_0 + _S3999.differential_0 + make_float3 (_S4037.x, _S4037.y, _S4038);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4047;
    (&_S4047)->primal_0 = R_29;
    (&_S4047)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4048;
    (&_S4048)->primal_0 = vert1_6;
    (&_S4048)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S4047, &_S4048, _S4046);
    float3  _S4049 = _S3921 + _S3922 + _S3999.differential_0 + make_float3 (_S4041.x, _S4041.y, _S4042);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4050;
    (&_S4050)->primal_0 = R_29;
    (&_S4050)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4051;
    (&_S4051)->primal_0 = vert0_6;
    (&_S4051)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S4050, &_S4051, _S4049);
    float3  _S4052 = (*v_verts_0)[int(2)] + _S4045.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4053;
    (&_S4053)->primal_0 = _S3833;
    (&_S4053)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4054;
    (&_S4054)->primal_0 = _S3838;
    (&_S4054)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S4053, &_S4054, _S4052);
    float _S4055 = - _S4054.differential_0.y;
    float _S4056 = _S3837 * _S4054.differential_0.x;
    float _S4057 = - (_S3829 * _S4054.differential_0.x);
    float3  _S4058 = (*v_verts_0)[int(1)] + _S4048.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4059;
    (&_S4059)->primal_0 = _S3833;
    (&_S4059)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4060;
    (&_S4060)->primal_0 = _S3836;
    (&_S4060)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S4059, &_S4060, _S4058);
    float _S4061 = _S3829 * _S4060.differential_0.x;
    float _S4062 = _S3835 * _S4060.differential_0.x;
    float3  _S4063 = (*v_verts_0)[int(0)] + _S4051.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4064;
    (&_S4064)->primal_0 = _S3833;
    (&_S4064)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4065;
    (&_S4065)->primal_0 = _S3834;
    (&_S4065)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S4064, &_S4065, _S4063);
    Matrix<float, 3, 3>  _S4066 = transpose_0(_S4053.differential_0 + _S4059.differential_0 + _S4064.differential_0);
    float _S4067 = 2.0f * - _S4066.rows[int(2)].z;
    float _S4068 = 2.0f * _S4066.rows[int(2)].y;
    float _S4069 = 2.0f * _S4066.rows[int(2)].x;
    float _S4070 = 2.0f * _S4066.rows[int(1)].z;
    float _S4071 = 2.0f * - _S4066.rows[int(1)].y;
    float _S4072 = 2.0f * _S4066.rows[int(1)].x;
    float _S4073 = 2.0f * _S4066.rows[int(0)].z;
    float _S4074 = 2.0f * _S4066.rows[int(0)].y;
    float _S4075 = 2.0f * - _S4066.rows[int(0)].x;
    float _S4076 = - _S4072 + _S4074;
    float _S4077 = _S4069 + - _S4073;
    float _S4078 = - _S4068 + _S4070;
    float _S4079 = _S4068 + _S4070;
    float _S4080 = _S4069 + _S4073;
    float _S4081 = _S4072 + _S4074;
    float _S4082 = quat_30.w * (_S4071 + _S4075);
    float _S4083 = quat_30.z * (_S4067 + _S4075);
    float _S4084 = quat_30.y * (_S4067 + _S4071);
    float _S4085 = quat_30.x * _S4076 + quat_30.z * _S4079 + quat_30.y * _S4080 + _S4082 + _S4082;
    float _S4086 = quat_30.x * _S4077 + quat_30.w * _S4079 + quat_30.y * _S4081 + _S4083 + _S4083;
    float _S4087 = quat_30.x * _S4078 + quat_30.w * _S4080 + quat_30.z * _S4081 + _S4084 + _S4084;
    float _S4088 = quat_30.w * _S4076 + quat_30.z * _S4077 + quat_30.y * _S4078;
    float _S4089 = _S4057 + _S4061;
    float _S4090 = 0.5f * - _S4089;
    float _S4091 = _S4055 + _S4060.differential_0.y;
    DiffPair_float_0 _S4092;
    (&_S4092)->primal_0 = _S3830;
    (&_S4092)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4092, _S4091);
    float _S4093 = _S4090 + _S4092.differential_0;
    float _S4094 = _S4056 + _S4062 + _S4065.differential_0.x;
    DiffPair_float_0 _S4095;
    (&_S4095)->primal_0 = _S3828;
    (&_S4095)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4095, _S4094);
    float _S4096 = _S4090 + _S4095.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4097;
    (&_S4097)->primal_0 = R_29;
    (&_S4097)->differential_0 = _S3994;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4098;
    (&_S4098)->primal_0 = mean_29;
    (&_S4098)->differential_0 = _S3914;
    s_bwd_prop_mul_1(&_S4097, &_S4098, _S3916.differential_0);
    float3  _S4099 = _S3996.differential_0 + _S4043 + _S4046 + _S4049 + _S3916.differential_0;
    Matrix<float, 3, 3>  _S4100 = _S3997 + _S4044.differential_0 + _S4047.differential_0 + _S4050.differential_0 + _S4097.differential_0;
    FixedArray<float3 , 2>  _S4101;
    _S4101[int(0)] = _S3914;
    _S4101[int(1)] = _S3914;
    _S4101[int(1)] = _S3928;
    _S4101[int(0)] = _S3933;
    FixedArray<float3 , 16>  _S4102;
    _S4102[int(0)] = _S3914;
    _S4102[int(1)] = _S3914;
    _S4102[int(2)] = _S3914;
    _S4102[int(3)] = _S3914;
    _S4102[int(4)] = _S3914;
    _S4102[int(5)] = _S3914;
    _S4102[int(6)] = _S3914;
    _S4102[int(7)] = _S3914;
    _S4102[int(8)] = _S3914;
    _S4102[int(9)] = _S3914;
    _S4102[int(10)] = _S3914;
    _S4102[int(11)] = _S3914;
    _S4102[int(12)] = _S3914;
    _S4102[int(13)] = _S3914;
    _S4102[int(14)] = _S3914;
    _S4102[int(15)] = _S3914;
    _S4102[int(15)] = _S3935;
    _S4102[int(14)] = _S3937;
    _S4102[int(13)] = _S3939;
    _S4102[int(12)] = _S3941;
    _S4102[int(11)] = _S3943;
    _S4102[int(10)] = _S3945;
    _S4102[int(9)] = _S3947;
    _S4102[int(8)] = _S3955;
    _S4102[int(7)] = _S3957;
    _S4102[int(6)] = _S3959;
    _S4102[int(5)] = _S3961;
    _S4102[int(4)] = _S3963;
    _S4102[int(3)] = _S3974;
    _S4102[int(2)] = _S3976;
    _S4102[int(1)] = _S3978;
    _S4102[int(0)] = _S3991;
    float2  _S4103 = make_float2 (0.0f, _S4023);
    float3  _S4104 = make_float3 (_S4096, _S4093, _S4089);
    float4  _S4105 = make_float4 (0.0f);
    *&((&_S4105)->w) = _S4085;
    *&((&_S4105)->z) = _S4086;
    *&((&_S4105)->y) = _S4087;
    *&((&_S4105)->x) = _S4088;
    *v_mean_9 = _S3992 + _S4052 + _S4058 + _S4063 + _S4098.differential_0;
    *v_quat_8 = _S4105;
    *v_scale_8 = _S4104;
    *v_hardness_4 = _S4103;
    *v_sh_coeffs_7 = _S4102;
    *v_ch_coeffs_2 = _S4101;
    *v_R_8 = _S4100;
    *v_t_8 = _S4099;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_30, float3  t_29, float fx_30, float fy_30, float cx_30, float cy_30, float4  radial_coeffs_33, float2  tangential_coeffs_33, float2  thin_prism_coeffs_33, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_1(R_30, mean_30) + t_29;
    float _S4106 = scale_30.x;
    float _S4107 = s_primal_ctx_exp_1(_S4106);
    float _S4108 = scale_30.y;
    float _S4109 = s_primal_ctx_exp_1(_S4108);
    float sz_11 = scale_30.z - 0.5f * (_S4106 + _S4108);
    float _S4110 = quat_31.y;
    float x2_31 = _S4110 * _S4110;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_57 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4111 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_57), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_57), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4112 = make_float3 (_S4107, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_1(_S4111, _S4112) + mean_30;
    float _S4113 = -0.5f + sz_11;
    float3  _S4114 = make_float3 (_S4107 * _S4113, _S4109, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_1(_S4111, _S4114) + mean_30;
    float _S4115 = -0.5f - sz_11;
    float3  _S4116 = make_float3 (_S4107 * _S4115, - _S4109, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_1(_S4111, _S4116) + mean_30;
    float3  vert0_c_11 = s_primal_ctx_mul_1(R_30, vert0_7) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_1(R_30, vert1_7) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_1(R_30, vert2_7) + t_29;
    CameraDistortion_0 _S4117 = s_primal_ctx_CameraDistortion_x24init_0(radial_coeffs_33, tangential_coeffs_33, thin_prism_coeffs_33);
    float2  _S4118 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4119 = length_0(_S4118);
    float _S4120 = vert0_c_11.z;
    float _S4121 = s_primal_ctx_atan2_0(_S4119, _S4120);
    bool _S4122 = _S4121 < 0.00100000004749745f;
    float k_9;
    float _S4123;
    float _S4124;
    float _S4125;
    if(_S4122)
    {
        float _S4126 = 1.0f - _S4121 * _S4121 / 3.0f;
        float _S4127 = _S4120 * _S4120;
        k_9 = _S4126 / _S4120;
        _S4123 = 0.0f;
        _S4124 = _S4127;
        _S4125 = _S4126;
    }
    else
    {
        float _S4128 = _S4119 * _S4119;
        k_9 = _S4121 / _S4119;
        _S4123 = _S4128;
        _S4124 = 0.0f;
        _S4125 = 0.0f;
    }
    float2  _S4129 = make_float2 (k_9);
    float2  _S4130 = _S4118 * make_float2 (k_9);
    float k1_8 = _S4117.radial_coeffs_0.x;
    float k2_8 = _S4117.radial_coeffs_0.y;
    float k3_8 = _S4117.radial_coeffs_0.z;
    float k4_9 = _S4117.radial_coeffs_0.w;
    float p1_9 = _S4117.tangential_coeffs_0.x;
    float p2_9 = _S4117.tangential_coeffs_0.y;
    float sx1_9 = _S4117.thin_prism_coeffs_0.x;
    float sy1_9 = _S4117.thin_prism_coeffs_0.y;
    float u_25 = _S4130.x;
    float v_25 = _S4130.y;
    float r2_25 = u_25 * u_25 + v_25 * v_25;
    float _S4131 = k3_8 + r2_25 * k4_9;
    float _S4132 = k2_8 + r2_25 * _S4131;
    float _S4133 = k1_8 + r2_25 * _S4132;
    float radial_5 = 1.0f + r2_25 * _S4133;
    float2  _S4134 = make_float2 (radial_5);
    float _S4135 = 2.0f * p1_9;
    float _S4136 = _S4135 * u_25;
    float _S4137 = 2.0f * u_25;
    float _S4138 = r2_25 + _S4137 * u_25;
    float _S4139 = 2.0f * p2_9;
    float _S4140 = _S4139 * u_25;
    float _S4141 = 2.0f * v_25;
    float _S4142 = r2_25 + _S4141 * v_25;
    float2  _S4143 = _S4130 * make_float2 (radial_5) + make_float2 (_S4136 * v_25 + p2_9 * _S4138 + sx1_9 * r2_25, _S4140 * v_25 + p1_9 * _S4142 + sy1_9 * r2_25);
    float _S4144 = fx_30 * _S4143.x + cx_30;
    float _S4145 = fy_30 * _S4143.y + cy_30;
    float2  _S4146 = make_float2 (_S4144, _S4145);
    float2  _S4147 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4148 = length_0(_S4147);
    float _S4149 = vert1_c_11.z;
    float _S4150 = s_primal_ctx_atan2_0(_S4148, _S4149);
    bool _S4151 = _S4150 < 0.00100000004749745f;
    float _S4152;
    float _S4153;
    float _S4154;
    if(_S4151)
    {
        float _S4155 = 1.0f - _S4150 * _S4150 / 3.0f;
        float _S4156 = _S4149 * _S4149;
        k_9 = _S4155 / _S4149;
        _S4152 = 0.0f;
        _S4153 = _S4156;
        _S4154 = _S4155;
    }
    else
    {
        float _S4157 = _S4148 * _S4148;
        k_9 = _S4150 / _S4148;
        _S4152 = _S4157;
        _S4153 = 0.0f;
        _S4154 = 0.0f;
    }
    float2  _S4158 = make_float2 (k_9);
    float2  _S4159 = _S4147 * make_float2 (k_9);
    float u_26 = _S4159.x;
    float v_26 = _S4159.y;
    float r2_26 = u_26 * u_26 + v_26 * v_26;
    float _S4160 = k3_8 + r2_26 * k4_9;
    float _S4161 = k2_8 + r2_26 * _S4160;
    float _S4162 = k1_8 + r2_26 * _S4161;
    float radial_6 = 1.0f + r2_26 * _S4162;
    float2  _S4163 = make_float2 (radial_6);
    float _S4164 = _S4135 * u_26;
    float _S4165 = 2.0f * u_26;
    float _S4166 = r2_26 + _S4165 * u_26;
    float _S4167 = _S4139 * u_26;
    float _S4168 = 2.0f * v_26;
    float _S4169 = r2_26 + _S4168 * v_26;
    float2  _S4170 = _S4159 * make_float2 (radial_6) + make_float2 (_S4164 * v_26 + p2_9 * _S4166 + sx1_9 * r2_26, _S4167 * v_26 + p1_9 * _S4169 + sy1_9 * r2_26);
    float _S4171 = fx_30 * _S4170.x + cx_30;
    float _S4172 = fy_30 * _S4170.y + cy_30;
    float2  _S4173 = make_float2 (_S4171, _S4172);
    float2  _S4174 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4175 = length_0(_S4174);
    float _S4176 = vert2_c_11.z;
    float _S4177 = s_primal_ctx_atan2_0(_S4175, _S4176);
    bool _S4178 = _S4177 < 0.00100000004749745f;
    float _S4179;
    float _S4180;
    float _S4181;
    if(_S4178)
    {
        float _S4182 = 1.0f - _S4177 * _S4177 / 3.0f;
        float _S4183 = _S4176 * _S4176;
        k_9 = _S4182 / _S4176;
        _S4179 = 0.0f;
        _S4180 = _S4183;
        _S4181 = _S4182;
    }
    else
    {
        float _S4184 = _S4175 * _S4175;
        k_9 = _S4177 / _S4175;
        _S4179 = _S4184;
        _S4180 = 0.0f;
        _S4181 = 0.0f;
    }
    float2  _S4185 = make_float2 (k_9);
    float2  _S4186 = _S4174 * make_float2 (k_9);
    float u_27 = _S4186.x;
    float v_27 = _S4186.y;
    float r2_27 = u_27 * u_27 + v_27 * v_27;
    float _S4187 = k3_8 + r2_27 * k4_9;
    float _S4188 = k2_8 + r2_27 * _S4187;
    float _S4189 = k1_8 + r2_27 * _S4188;
    float radial_7 = 1.0f + r2_27 * _S4189;
    float2  _S4190 = make_float2 (radial_7);
    float _S4191 = _S4135 * u_27;
    float _S4192 = 2.0f * u_27;
    float _S4193 = r2_27 + _S4192 * u_27;
    float _S4194 = _S4139 * u_27;
    float _S4195 = 2.0f * v_27;
    float _S4196 = r2_27 + _S4195 * v_27;
    float2  _S4197 = _S4186 * make_float2 (radial_7) + make_float2 (_S4191 * v_27 + p2_9 * _S4193 + sx1_9 * r2_27, _S4194 * v_27 + p1_9 * _S4196 + sy1_9 * r2_27);
    float _S4198 = fx_30 * _S4197.x + cx_30;
    float _S4199 = fy_30 * _S4197.y + cy_30;
    float2  _S4200 = make_float2 (_S4198, _S4199);
    float2  e0_15 = _S4173 - _S4146;
    float2  e1_15 = _S4200 - _S4173;
    float2  e2_7 = _S4146 - _S4200;
    float _S4201 = e0_15.x;
    float _S4202 = e1_15.y;
    float _S4203 = e0_15.y;
    float _S4204 = e1_15.x;
    float _S4205 = _S4201 * _S4202 - _S4203 * _S4204;
    float _S4206 = 1.0f - hardness_15.y;
    float _S4207 = -1.0f / _S4206;
    float _S4208 = _S4206 * _S4206;
    float _S4209 = s_primal_ctx_max_0(_S4144, _S4171);
    float _S4210 = s_primal_ctx_min_0(_S4144, _S4171);
    float _S4211 = s_primal_ctx_max_0(_S4145, _S4172);
    float _S4212 = s_primal_ctx_min_0(_S4145, _S4172);
    float3  _S4213 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    Matrix<float, 3, 3>  _S4214 = transpose_0(R_30);
    float3  _S4215 = mean_30 - - s_primal_ctx_mul_1(_S4214, t_29);
    float _S4216 = _S4215.x;
    float _S4217 = _S4215.y;
    float _S4218 = _S4215.z;
    float _S4219 = _S4216 * _S4216 + _S4217 * _S4217 + _S4218 * _S4218;
    float _S4220 = s_primal_ctx_sqrt_0(_S4219);
    float x_60 = _S4216 / _S4220;
    float3  _S4221 = make_float3 (x_60);
    float _S4222 = _S4220 * _S4220;
    float y_29 = _S4217 / _S4220;
    float z_26 = _S4218 / _S4220;
    float3  _S4223 = make_float3 (z_26);
    float _S4224 = - y_29;
    float3  _S4225 = make_float3 (_S4224);
    float z2_58 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_60 * x_60 - y_29 * y_29;
    float _S4226 = 2.0f * x_60;
    float fS1_26 = _S4226 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_58 - 0.31539157032966614f;
    float3  _S4227 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_60;
    float3  _S4228 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4229 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4230 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4231 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_58 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4232 = 1.86588168144226074f * z2_58 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4232;
    float3  _S4233 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_60;
    float3  _S4234 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4235 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4236 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4237 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_60 * fC1_26 - y_29 * fS1_26);
    float3  _S4238 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_60 * fS1_26 + y_29 * fC1_26);
    float3  _S4239 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4224) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_60) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4240 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4241 = make_float3 (0.0f);
    float3  _S4242 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4243 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4244 = make_float3 (_S4243);
    float3  _S4245 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4243);
    float3  _S4246 = _S4242 + _S4245 + make_float3 (0.5f);
    float3  _S4247 = _S4242 - _S4245 + make_float3 (0.5f);
    float3  _S4248 = vert1_c_11 - vert0_c_11;
    float3  _S4249 = vert2_c_11 - vert0_c_11;
    float3  _S4250 = s_primal_ctx_cross_0(_S4248, _S4249);
    float3  _S4251 = normalize_0(_S4250);
    float3  _S4252 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4251, mean_c_26)))))) * v_normal_3;
    float3  _S4253 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4254;
    (&_S4254)->primal_0 = _S4251;
    (&_S4254)->differential_0 = _S4253;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4255;
    (&_S4255)->primal_0 = mean_c_26;
    (&_S4255)->differential_0 = _S4253;
    s_bwd_prop_dot_0(&_S4254, &_S4255, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4256 = _S4255;
    float3  _S4257 = _S4252 + _S4254.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4258;
    (&_S4258)->primal_0 = _S4250;
    (&_S4258)->differential_0 = _S4253;
    s_bwd_normalize_impl_0(&_S4258, _S4257);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4259;
    (&_S4259)->primal_0 = _S4248;
    (&_S4259)->differential_0 = _S4253;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4260;
    (&_S4260)->primal_0 = _S4249;
    (&_S4260)->differential_0 = _S4253;
    s_bwd_prop_cross_0(&_S4259, &_S4260, _S4258.differential_0);
    float3  _S4261 = - _S4260.differential_0;
    float3  _S4262 = - _S4259.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4263;
    (&_S4263)->primal_0 = _S4247;
    (&_S4263)->differential_0 = _S4253;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4264;
    (&_S4264)->primal_0 = _S4241;
    (&_S4264)->differential_0 = _S4253;
    s_bwd_prop_max_0(&_S4263, &_S4264, (*v_rgbs_1)[int(2)]);
    float3  _S4265 = - _S4263.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4266;
    (&_S4266)->primal_0 = _S4246;
    (&_S4266)->differential_0 = _S4253;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4267;
    (&_S4267)->primal_0 = _S4241;
    (&_S4267)->differential_0 = _S4253;
    s_bwd_prop_max_0(&_S4266, &_S4267, (*v_rgbs_1)[int(1)]);
    float3  _S4268 = _S4244 * (_S4265 + _S4266.differential_0);
    float3  _S4269 = _S4263.differential_0 + _S4266.differential_0;
    float3  _S4270 = make_float3 (0.5f) * - _S4269;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4271;
    (&_S4271)->primal_0 = _S4240;
    (&_S4271)->differential_0 = _S4253;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4272;
    (&_S4272)->primal_0 = _S4241;
    (&_S4272)->differential_0 = _S4253;
    s_bwd_prop_max_0(&_S4271, &_S4272, (*v_rgbs_1)[int(0)]);
    float3  _S4273 = _S4270 + _S4271.differential_0;
    float3  _S4274 = _S4269 + _S4271.differential_0;
    float3  _S4275 = _S4238 * _S4274;
    float3  _S4276 = (*sh_coeffs_26)[int(15)] * _S4274;
    float3  _S4277 = _S4236 * _S4274;
    float3  _S4278 = (*sh_coeffs_26)[int(14)] * _S4274;
    float3  _S4279 = _S4234 * _S4274;
    float3  _S4280 = (*sh_coeffs_26)[int(13)] * _S4274;
    float3  _S4281 = _S4233 * _S4274;
    float3  _S4282 = (*sh_coeffs_26)[int(12)] * _S4274;
    float3  _S4283 = _S4235 * _S4274;
    float3  _S4284 = (*sh_coeffs_26)[int(11)] * _S4274;
    float3  _S4285 = _S4237 * _S4274;
    float3  _S4286 = (*sh_coeffs_26)[int(10)] * _S4274;
    float3  _S4287 = _S4239 * _S4274;
    float3  _S4288 = (*sh_coeffs_26)[int(9)] * _S4274;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4288.x + _S4288.y + _S4288.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4276.x + _S4276.y + _S4276.z);
    float _S4289 = _S4286.x + _S4286.y + _S4286.z;
    float _S4290 = _S4278.x + _S4278.y + _S4278.z;
    float _S4291 = _S4284.x + _S4284.y + _S4284.z;
    float _S4292 = _S4280.x + _S4280.y + _S4280.z;
    float _S4293 = _S4282.x + _S4282.y + _S4282.z;
    float _S4294 = - s_diff_fC2_T_8;
    float3  _S4295 = _S4230 * _S4274;
    float3  _S4296 = (*sh_coeffs_26)[int(8)] * _S4274;
    float3  _S4297 = _S4228 * _S4274;
    float3  _S4298 = (*sh_coeffs_26)[int(7)] * _S4274;
    float3  _S4299 = _S4227 * _S4274;
    float3  _S4300 = (*sh_coeffs_26)[int(6)] * _S4274;
    float3  _S4301 = _S4229 * _S4274;
    float3  _S4302 = (*sh_coeffs_26)[int(5)] * _S4274;
    float3  _S4303 = _S4231 * _S4274;
    float3  _S4304 = (*sh_coeffs_26)[int(4)] * _S4274;
    float _S4305 = _S4302.x + _S4302.y + _S4302.z;
    float _S4306 = _S4298.x + _S4298.y + _S4298.z;
    float _S4307 = fTmp1B_26 * _S4289 + x_60 * s_diff_fS2_T_8 + y_29 * _S4294 + 0.54627424478530884f * (_S4304.x + _S4304.y + _S4304.z);
    float _S4308 = fTmp1B_26 * _S4290 + y_29 * s_diff_fS2_T_8 + x_60 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S4296.x + _S4296.y + _S4296.z);
    float _S4309 = y_29 * - _S4308;
    float _S4310 = x_60 * _S4308;
    float _S4311 = z_26 * (1.86588168144226074f * (z_26 * _S4293) + -2.28522896766662598f * (y_29 * _S4291 + x_60 * _S4292) + 0.94617468118667603f * (_S4300.x + _S4300.y + _S4300.z));
    float3  _S4312 = make_float3 (0.48860251903533936f) * _S4274;
    float3  _S4313 = - _S4312;
    float3  _S4314 = _S4221 * _S4313;
    float3  _S4315 = (*sh_coeffs_26)[int(3)] * _S4313;
    float3  _S4316 = _S4223 * _S4312;
    float3  _S4317 = (*sh_coeffs_26)[int(2)] * _S4312;
    float3  _S4318 = _S4225 * _S4312;
    float3  _S4319 = (*sh_coeffs_26)[int(1)] * _S4312;
    float _S4320 = (_S4232 * _S4293 + 1.44530570507049561f * (fS1_26 * _S4289 + fC1_26 * _S4290) + -1.09254848957061768f * (y_29 * _S4305 + x_60 * _S4306) + _S4311 + _S4311 + _S4317.x + _S4317.y + _S4317.z) / _S4222;
    float _S4321 = _S4220 * _S4320;
    float _S4322 = (fTmp0C_26 * _S4291 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S4294 + fTmp0B_26 * _S4305 + _S4226 * _S4307 + _S4309 + _S4309 + - (_S4319.x + _S4319.y + _S4319.z)) / _S4222;
    float _S4323 = _S4220 * _S4322;
    float _S4324 = (fTmp0C_26 * _S4292 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S4306 + 2.0f * (y_29 * _S4307) + _S4310 + _S4310 + _S4315.x + _S4315.y + _S4315.z) / _S4222;
    float _S4325 = _S4220 * _S4324;
    float _S4326 = _S4218 * - _S4320 + _S4217 * - _S4322 + _S4216 * - _S4324;
    DiffPair_float_0 _S4327;
    (&_S4327)->primal_0 = _S4219;
    (&_S4327)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4327, _S4326);
    float _S4328 = _S4218 * _S4327.differential_0;
    float _S4329 = _S4217 * _S4327.differential_0;
    float _S4330 = _S4216 * _S4327.differential_0;
    float3  _S4331 = make_float3 (0.282094806432724f) * _S4274;
    float3  _S4332 = make_float3 (_S4325 + _S4330 + _S4330, _S4323 + _S4329 + _S4329, _S4321 + _S4328 + _S4328);
    float3  _S4333 = - - _S4332;
    Matrix<float, 3, 3>  _S4334 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4335;
    (&_S4335)->primal_0 = _S4214;
    (&_S4335)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4336;
    (&_S4336)->primal_0 = t_29;
    (&_S4336)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4335, &_S4336, _S4333);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4337 = _S4336;
    Matrix<float, 3, 3>  _S4338 = transpose_0(_S4335.differential_0);
    float _S4339 = 0.3333333432674408f * v_depth_10;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4340;
    (&_S4340)->primal_0 = _S4213;
    (&_S4340)->differential_0 = _S4253;
    s_bwd_length_impl_0(&_S4340, _S4339);
    DiffPair_float_0 _S4341;
    (&_S4341)->primal_0 = _S4212;
    (&_S4341)->differential_0 = 0.0f;
    DiffPair_float_0 _S4342;
    (&_S4342)->primal_0 = _S4199;
    (&_S4342)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4341, &_S4342, 0.0f);
    DiffPair_float_0 _S4343;
    (&_S4343)->primal_0 = _S4145;
    (&_S4343)->differential_0 = 0.0f;
    DiffPair_float_0 _S4344;
    (&_S4344)->primal_0 = _S4172;
    (&_S4344)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4343, &_S4344, _S4341.differential_0);
    DiffPair_float_0 _S4345;
    (&_S4345)->primal_0 = _S4211;
    (&_S4345)->differential_0 = 0.0f;
    DiffPair_float_0 _S4346;
    (&_S4346)->primal_0 = _S4199;
    (&_S4346)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4345, &_S4346, 0.0f);
    DiffPair_float_0 _S4347;
    (&_S4347)->primal_0 = _S4145;
    (&_S4347)->differential_0 = 0.0f;
    DiffPair_float_0 _S4348;
    (&_S4348)->primal_0 = _S4172;
    (&_S4348)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4347, &_S4348, _S4345.differential_0);
    DiffPair_float_0 _S4349;
    (&_S4349)->primal_0 = _S4210;
    (&_S4349)->differential_0 = 0.0f;
    DiffPair_float_0 _S4350;
    (&_S4350)->primal_0 = _S4198;
    (&_S4350)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4349, &_S4350, 0.0f);
    DiffPair_float_0 _S4351;
    (&_S4351)->primal_0 = _S4144;
    (&_S4351)->differential_0 = 0.0f;
    DiffPair_float_0 _S4352;
    (&_S4352)->primal_0 = _S4171;
    (&_S4352)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4351, &_S4352, _S4349.differential_0);
    DiffPair_float_0 _S4353;
    (&_S4353)->primal_0 = _S4209;
    (&_S4353)->differential_0 = 0.0f;
    DiffPair_float_0 _S4354;
    (&_S4354)->primal_0 = _S4198;
    (&_S4354)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4353, &_S4354, 0.0f);
    DiffPair_float_0 _S4355;
    (&_S4355)->primal_0 = _S4144;
    (&_S4355)->differential_0 = 0.0f;
    DiffPair_float_0 _S4356;
    (&_S4356)->primal_0 = _S4171;
    (&_S4356)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4355, &_S4356, _S4353.differential_0);
    DiffPair_float_0 _S4357;
    (&_S4357)->primal_0 = _S4207;
    (&_S4357)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4357, 0.0f);
    float _S4358 = - (-1.0f * - (_S4357.differential_0 / _S4208));
    float2  _S4359 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4360;
    (&_S4360)->primal_0 = e2_7;
    (&_S4360)->differential_0 = _S4359;
    s_bwd_length_impl_1(&_S4360, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4361;
    (&_S4361)->primal_0 = e1_15;
    (&_S4361)->differential_0 = _S4359;
    s_bwd_length_impl_1(&_S4361, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4362;
    (&_S4362)->primal_0 = e0_15;
    (&_S4362)->differential_0 = _S4359;
    s_bwd_length_impl_1(&_S4362, -0.0f);
    DiffPair_float_0 _S4363;
    (&_S4363)->primal_0 = _S4205;
    (&_S4363)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4363, 0.0f);
    float _S4364 = - _S4363.differential_0;
    float2  _S4365 = _S4361.differential_0 + make_float2 (_S4203 * _S4364, _S4201 * _S4363.differential_0);
    float2  _S4366 = - _S4365;
    float2  _S4367 = _S4362.differential_0 + make_float2 (_S4202 * _S4363.differential_0, _S4204 * _S4364);
    float2  _S4368 = - _S4367;
    float2  _S4369 = - _S4360.differential_0 + _S4365;
    float _S4370 = fy_30 * (_S4342.differential_0 + _S4346.differential_0 + _S4369.y);
    float _S4371 = fx_30 * (_S4350.differential_0 + _S4354.differential_0 + _S4369.x);
    float2  _S4372 = make_float2 (_S4371, _S4370);
    float2  _S4373 = _S4186 * _S4372;
    float2  _S4374 = _S4190 * _S4372;
    float _S4375 = r2_27 * _S4370;
    float _S4376 = p1_9 * _S4370;
    float _S4377 = _S4196 * _S4370;
    float _S4378 = v_27 * _S4370;
    float _S4379 = u_27 * _S4378;
    float _S4380 = r2_27 * _S4371;
    float _S4381 = p2_9 * _S4371;
    float _S4382 = _S4193 * _S4371;
    float _S4383 = v_27 * _S4371;
    float _S4384 = u_27 * _S4383;
    float _S4385 = _S4373.x + _S4373.y;
    float _S4386 = r2_27 * _S4385;
    float _S4387 = r2_27 * _S4386;
    float _S4388 = r2_27 * _S4387;
    float _S4389 = r2_27 * _S4388;
    float _S4390 = sy1_9 * _S4370 + _S4376 + sx1_9 * _S4371 + _S4381 + _S4189 * _S4385 + _S4188 * _S4386 + _S4187 * _S4387 + k4_9 * _S4388;
    float _S4391 = v_27 * _S4390;
    float _S4392 = u_27 * _S4390;
    float _S4393 = _S4195 * _S4376 + 2.0f * (v_27 * _S4376) + _S4194 * _S4370 + _S4191 * _S4371 + _S4391 + _S4391;
    float _S4394 = _S4139 * _S4378 + _S4192 * _S4381 + 2.0f * (u_27 * _S4381) + _S4135 * _S4383 + _S4392 + _S4392;
    float3  _S4395 = _S4259.differential_0 + _S4340.differential_0;
    float3  _S4396 = _S4261 + _S4262 + _S4340.differential_0;
    float3  _S4397 = _S4260.differential_0 + _S4340.differential_0;
    FixedArray<float3 , 2>  _S4398;
    _S4398[int(0)] = _S4253;
    _S4398[int(1)] = _S4253;
    _S4398[int(1)] = _S4268;
    _S4398[int(0)] = _S4273;
    float3  _S4399 = _S4398[int(0)];
    float3  _S4400 = _S4398[int(1)];
    FixedArray<float3 , 16>  _S4401;
    _S4401[int(0)] = _S4253;
    _S4401[int(1)] = _S4253;
    _S4401[int(2)] = _S4253;
    _S4401[int(3)] = _S4253;
    _S4401[int(4)] = _S4253;
    _S4401[int(5)] = _S4253;
    _S4401[int(6)] = _S4253;
    _S4401[int(7)] = _S4253;
    _S4401[int(8)] = _S4253;
    _S4401[int(9)] = _S4253;
    _S4401[int(10)] = _S4253;
    _S4401[int(11)] = _S4253;
    _S4401[int(12)] = _S4253;
    _S4401[int(13)] = _S4253;
    _S4401[int(14)] = _S4253;
    _S4401[int(15)] = _S4253;
    _S4401[int(7)] = _S4297;
    _S4401[int(0)] = _S4331;
    _S4401[int(1)] = _S4318;
    _S4401[int(2)] = _S4316;
    _S4401[int(3)] = _S4314;
    _S4401[int(4)] = _S4303;
    _S4401[int(5)] = _S4301;
    _S4401[int(6)] = _S4299;
    _S4401[int(15)] = _S4275;
    _S4401[int(8)] = _S4295;
    _S4401[int(9)] = _S4287;
    _S4401[int(10)] = _S4285;
    _S4401[int(11)] = _S4283;
    _S4401[int(12)] = _S4281;
    _S4401[int(13)] = _S4279;
    _S4401[int(14)] = _S4277;
    float3  _S4402 = _S4401[int(0)];
    float3  _S4403 = _S4401[int(1)];
    float3  _S4404 = _S4401[int(2)];
    float3  _S4405 = _S4401[int(3)];
    float3  _S4406 = _S4401[int(4)];
    float3  _S4407 = _S4401[int(5)];
    float3  _S4408 = _S4401[int(6)];
    float3  _S4409 = _S4401[int(7)];
    float3  _S4410 = _S4401[int(8)];
    float3  _S4411 = _S4401[int(9)];
    float3  _S4412 = _S4401[int(10)];
    float3  _S4413 = _S4401[int(11)];
    float3  _S4414 = _S4401[int(12)];
    float3  _S4415 = _S4401[int(13)];
    float3  _S4416 = _S4401[int(14)];
    float3  _S4417 = _S4401[int(15)];
    float _S4418 = _S4351.differential_0 + _S4355.differential_0;
    float2  _S4419 = _S4360.differential_0 + _S4368;
    float _S4420 = _S4344.differential_0 + _S4348.differential_0;
    float _S4421 = _S4343.differential_0 + _S4347.differential_0;
    float2  _S4422 = _S4366 + _S4367;
    float _S4423 = _S4352.differential_0 + _S4356.differential_0;
    float2  _S4424 = make_float2 (0.0f, _S4358);
    float2  _S4425 = _S4374 + make_float2 (_S4394, _S4393);
    float2  _S4426 = _S4174 * _S4425;
    float2  _S4427 = _S4185 * _S4425;
    float _S4428 = _S4426.x + _S4426.y;
    if(_S4178)
    {
        float _S4429 = _S4428 / _S4180;
        float _S4430 = _S4181 * - _S4429;
        float _S4431 = _S4177 * (0.3333333432674408f * - (_S4176 * _S4429));
        k_9 = _S4431 + _S4431;
        _S4179 = _S4430;
        _S4180 = 0.0f;
    }
    else
    {
        float _S4432 = _S4428 / _S4179;
        float _S4433 = _S4177 * - _S4432;
        k_9 = _S4175 * _S4432;
        _S4179 = 0.0f;
        _S4180 = _S4433;
    }
    DiffPair_float_0 _S4434;
    (&_S4434)->primal_0 = _S4175;
    (&_S4434)->differential_0 = 0.0f;
    DiffPair_float_0 _S4435;
    (&_S4435)->primal_0 = _S4176;
    (&_S4435)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4434, &_S4435, k_9);
    float _S4436 = _S4435.differential_0 + _S4179;
    float _S4437 = _S4434.differential_0 + _S4180;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4438;
    (&_S4438)->primal_0 = _S4174;
    (&_S4438)->differential_0 = _S4359;
    s_bwd_length_impl_1(&_S4438, _S4437);
    float2  _S4439 = _S4438.differential_0 + _S4427;
    float _S4440 = fy_30 * (_S4422.y + _S4420);
    float _S4441 = fx_30 * (_S4422.x + _S4423);
    float2  _S4442 = make_float2 (_S4441, _S4440);
    float2  _S4443 = _S4159 * _S4442;
    float _S4444 = p1_9 * _S4440;
    float _S4445 = v_26 * _S4440;
    float _S4446 = p2_9 * _S4441;
    float _S4447 = v_26 * _S4441;
    float _S4448 = _S4443.x + _S4443.y;
    float _S4449 = r2_26 * _S4448;
    float _S4450 = r2_26 * _S4449;
    float _S4451 = r2_26 * _S4450;
    float _S4452 = sy1_9 * _S4440 + _S4444 + sx1_9 * _S4441 + _S4446 + _S4162 * _S4448 + _S4161 * _S4449 + _S4160 * _S4450 + k4_9 * _S4451;
    float _S4453 = v_26 * _S4452;
    float _S4454 = u_26 * _S4452;
    float3  _S4455 = _S4397 + make_float3 (_S4439.x, _S4439.y, _S4436);
    float2  _S4456 = _S4163 * _S4442 + make_float2 (_S4139 * _S4445 + _S4165 * _S4446 + 2.0f * (u_26 * _S4446) + _S4135 * _S4447 + _S4454 + _S4454, _S4168 * _S4444 + 2.0f * (v_26 * _S4444) + _S4167 * _S4440 + _S4164 * _S4441 + _S4453 + _S4453);
    float _S4457 = u_26 * _S4445 + _S4379;
    float _S4458 = u_26 * _S4447 + _S4384;
    float _S4459 = r2_26 * _S4440 + _S4375;
    float _S4460 = r2_26 * _S4451 + _S4389;
    float _S4461 = _S4451 + _S4388;
    float _S4462 = _S4169 * _S4440 + _S4377;
    float _S4463 = _S4450 + _S4387;
    float _S4464 = r2_26 * _S4441 + _S4380;
    float _S4465 = _S4166 * _S4441 + _S4382;
    float _S4466 = _S4449 + _S4386;
    float2  _S4467 = _S4147 * _S4456;
    float2  _S4468 = _S4158 * _S4456;
    float _S4469 = _S4467.x + _S4467.y;
    if(_S4151)
    {
        float _S4470 = _S4469 / _S4153;
        float _S4471 = _S4154 * - _S4470;
        float _S4472 = _S4150 * (0.3333333432674408f * - (_S4149 * _S4470));
        k_9 = _S4472 + _S4472;
        _S4152 = _S4471;
        _S4153 = 0.0f;
    }
    else
    {
        float _S4473 = _S4469 / _S4152;
        float _S4474 = _S4150 * - _S4473;
        k_9 = _S4148 * _S4473;
        _S4152 = 0.0f;
        _S4153 = _S4474;
    }
    DiffPair_float_0 _S4475;
    (&_S4475)->primal_0 = _S4148;
    (&_S4475)->differential_0 = 0.0f;
    DiffPair_float_0 _S4476;
    (&_S4476)->primal_0 = _S4149;
    (&_S4476)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4475, &_S4476, k_9);
    float _S4477 = _S4476.differential_0 + _S4152;
    float _S4478 = _S4475.differential_0 + _S4153;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4479;
    (&_S4479)->primal_0 = _S4147;
    (&_S4479)->differential_0 = _S4359;
    s_bwd_length_impl_1(&_S4479, _S4478);
    float2  _S4480 = _S4479.differential_0 + _S4468;
    float _S4481 = fy_30 * (_S4419.y + _S4421);
    float _S4482 = fx_30 * (_S4419.x + _S4418);
    float2  _S4483 = make_float2 (_S4482, _S4481);
    float2  _S4484 = _S4130 * _S4483;
    float _S4485 = p1_9 * _S4481;
    float _S4486 = v_25 * _S4481;
    float _S4487 = p2_9 * _S4482;
    float _S4488 = v_25 * _S4482;
    float _S4489 = _S4484.x + _S4484.y;
    float _S4490 = r2_25 * _S4489;
    float _S4491 = r2_25 * _S4490;
    float _S4492 = r2_25 * _S4491;
    float _S4493 = sy1_9 * _S4481 + _S4485 + sx1_9 * _S4482 + _S4487 + _S4133 * _S4489 + _S4132 * _S4490 + _S4131 * _S4491 + k4_9 * _S4492;
    float _S4494 = v_25 * _S4493;
    float _S4495 = u_25 * _S4493;
    float2  _S4496 = make_float2 (r2_25 * _S4482 + _S4464, r2_25 * _S4481 + _S4459);
    float2  _S4497 = make_float2 (_S4142 * _S4481 + 2.0f * (u_25 * _S4488 + _S4458) + _S4462, 2.0f * (u_25 * _S4486 + _S4457) + _S4138 * _S4482 + _S4465);
    float4  _S4498 = make_float4 (_S4490 + _S4466, _S4491 + _S4463, _S4492 + _S4461, r2_25 * _S4492 + _S4460);
    float3  _S4499 = _S4395 + make_float3 (_S4480.x, _S4480.y, _S4477);
    float2  _S4500 = _S4134 * _S4483 + make_float2 (_S4139 * _S4486 + _S4137 * _S4487 + 2.0f * (u_25 * _S4487) + _S4135 * _S4488 + _S4495 + _S4495, _S4141 * _S4485 + 2.0f * (v_25 * _S4485) + _S4140 * _S4481 + _S4136 * _S4482 + _S4494 + _S4494);
    CameraDistortion_0 _S4501 = CameraDistortion_x24_syn_dzero_0();
    (&_S4501)->thin_prism_coeffs_0 = _S4496;
    (&_S4501)->tangential_coeffs_0 = _S4497;
    (&_S4501)->radial_coeffs_0 = _S4498;
    float2  _S4502 = _S4118 * _S4500;
    float2  _S4503 = _S4129 * _S4500;
    float _S4504 = _S4502.x + _S4502.y;
    if(_S4122)
    {
        float _S4505 = _S4504 / _S4124;
        float _S4506 = _S4125 * - _S4505;
        float _S4507 = _S4121 * (0.3333333432674408f * - (_S4120 * _S4505));
        k_9 = _S4507 + _S4507;
        _S4123 = _S4506;
        _S4124 = 0.0f;
    }
    else
    {
        float _S4508 = _S4504 / _S4123;
        float _S4509 = _S4121 * - _S4508;
        k_9 = _S4119 * _S4508;
        _S4123 = 0.0f;
        _S4124 = _S4509;
    }
    DiffPair_float_0 _S4510;
    (&_S4510)->primal_0 = _S4119;
    (&_S4510)->differential_0 = 0.0f;
    DiffPair_float_0 _S4511;
    (&_S4511)->primal_0 = _S4120;
    (&_S4511)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S4510, &_S4511, k_9);
    float _S4512 = _S4511.differential_0 + _S4123;
    float _S4513 = _S4510.differential_0 + _S4124;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4514;
    (&_S4514)->primal_0 = _S4118;
    (&_S4514)->differential_0 = _S4359;
    s_bwd_length_impl_1(&_S4514, _S4513);
    float2  _S4515 = _S4514.differential_0 + _S4503;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4516;
    (&_S4516)->primal_0 = vert2_c_11;
    (&_S4516)->differential_0 = _S4253;
    s_bwd_length_impl_0(&_S4516, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4517;
    (&_S4517)->primal_0 = vert1_c_11;
    (&_S4517)->differential_0 = _S4253;
    s_bwd_length_impl_0(&_S4517, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4518;
    (&_S4518)->primal_0 = vert0_c_11;
    (&_S4518)->differential_0 = _S4253;
    s_bwd_length_impl_0(&_S4518, 0.0f);
    float3  _S4519 = _S4516.differential_0 + _S4455;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4520;
    (&_S4520)->primal_0 = R_30;
    (&_S4520)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4521;
    (&_S4521)->primal_0 = vert2_7;
    (&_S4521)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4520, &_S4521, _S4519);
    float3  _S4522 = _S4517.differential_0 + _S4499;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4523;
    (&_S4523)->primal_0 = R_30;
    (&_S4523)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4524;
    (&_S4524)->primal_0 = vert1_7;
    (&_S4524)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4523, &_S4524, _S4522);
    float3  _S4525 = _S4518.differential_0 + _S4396 + make_float3 (_S4515.x, _S4515.y, _S4512);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4526;
    (&_S4526)->primal_0 = R_30;
    (&_S4526)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4527;
    (&_S4527)->primal_0 = vert0_7;
    (&_S4527)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4526, &_S4527, _S4525);
    float3  _S4528 = _S4521.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4529;
    (&_S4529)->primal_0 = _S4111;
    (&_S4529)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4530;
    (&_S4530)->primal_0 = _S4116;
    (&_S4530)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4529, &_S4530, _S4528);
    float _S4531 = - _S4530.differential_0.y;
    float _S4532 = _S4115 * _S4530.differential_0.x;
    float _S4533 = - (_S4107 * _S4530.differential_0.x);
    float3  _S4534 = _S4524.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4535;
    (&_S4535)->primal_0 = _S4111;
    (&_S4535)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4536;
    (&_S4536)->primal_0 = _S4114;
    (&_S4536)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4535, &_S4536, _S4534);
    float _S4537 = _S4107 * _S4536.differential_0.x;
    float _S4538 = _S4113 * _S4536.differential_0.x;
    float3  _S4539 = _S4527.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4540;
    (&_S4540)->primal_0 = _S4111;
    (&_S4540)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4541;
    (&_S4541)->primal_0 = _S4112;
    (&_S4541)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4540, &_S4541, _S4539);
    Matrix<float, 3, 3>  _S4542 = transpose_0(_S4529.differential_0 + _S4535.differential_0 + _S4540.differential_0);
    float _S4543 = 2.0f * - _S4542.rows[int(2)].z;
    float _S4544 = 2.0f * _S4542.rows[int(2)].y;
    float _S4545 = 2.0f * _S4542.rows[int(2)].x;
    float _S4546 = 2.0f * _S4542.rows[int(1)].z;
    float _S4547 = 2.0f * - _S4542.rows[int(1)].y;
    float _S4548 = 2.0f * _S4542.rows[int(1)].x;
    float _S4549 = 2.0f * _S4542.rows[int(0)].z;
    float _S4550 = 2.0f * _S4542.rows[int(0)].y;
    float _S4551 = 2.0f * - _S4542.rows[int(0)].x;
    float _S4552 = - _S4548 + _S4550;
    float _S4553 = _S4545 + - _S4549;
    float _S4554 = - _S4544 + _S4546;
    float _S4555 = _S4544 + _S4546;
    float _S4556 = _S4545 + _S4549;
    float _S4557 = _S4548 + _S4550;
    float _S4558 = quat_31.w * (_S4547 + _S4551);
    float _S4559 = quat_31.z * (_S4543 + _S4551);
    float _S4560 = quat_31.y * (_S4543 + _S4547);
    float _S4561 = quat_31.x * _S4552 + quat_31.z * _S4555 + quat_31.y * _S4556 + _S4558 + _S4558;
    float _S4562 = quat_31.x * _S4553 + quat_31.w * _S4555 + quat_31.y * _S4557 + _S4559 + _S4559;
    float _S4563 = quat_31.x * _S4554 + quat_31.w * _S4556 + quat_31.z * _S4557 + _S4560 + _S4560;
    float _S4564 = quat_31.w * _S4552 + quat_31.z * _S4553 + quat_31.y * _S4554;
    float _S4565 = _S4533 + _S4537;
    float _S4566 = 0.5f * - _S4565;
    float _S4567 = _S4531 + _S4536.differential_0.y;
    DiffPair_float_0 _S4568;
    (&_S4568)->primal_0 = _S4108;
    (&_S4568)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4568, _S4567);
    float _S4569 = _S4566 + _S4568.differential_0;
    float _S4570 = _S4532 + _S4538 + _S4541.differential_0.x;
    DiffPair_float_0 _S4571;
    (&_S4571)->primal_0 = _S4106;
    (&_S4571)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4571, _S4570);
    float _S4572 = _S4566 + _S4571.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4573;
    (&_S4573)->primal_0 = mean_c_26;
    (&_S4573)->differential_0 = _S4253;
    s_bwd_length_impl_0(&_S4573, 0.0f);
    float3  _S4574 = _S4573.differential_0 + _S4256.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4575;
    (&_S4575)->primal_0 = R_30;
    (&_S4575)->differential_0 = _S4334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4576;
    (&_S4576)->primal_0 = mean_30;
    (&_S4576)->differential_0 = _S4253;
    s_bwd_prop_mul_1(&_S4575, &_S4576, _S4574);
    float3  _S4577 = _S4519 + _S4522 + _S4525 + _S4574 + _S4337.differential_0;
    Matrix<float, 3, 3>  _S4578 = _S4520.differential_0 + _S4523.differential_0 + _S4526.differential_0 + _S4575.differential_0 + _S4338;
    float3  _S4579 = make_float3 (_S4572, _S4569, _S4565);
    float4  _S4580 = make_float4 (0.0f);
    *&((&_S4580)->w) = _S4561;
    *&((&_S4580)->z) = _S4562;
    *&((&_S4580)->y) = _S4563;
    *&((&_S4580)->x) = _S4564;
    float4  _S4581 = _S4580;
    float3  _S4582 = _S4528 + _S4534 + _S4539 + _S4576.differential_0 + _S4332;
    *v_mean_10 = _S4582;
    *v_quat_9 = _S4581;
    *v_scale_9 = _S4579;
    *v_hardness_5 = _S4424;
    (*v_sh_coeffs_8)[int(0)] = _S4402;
    (*v_sh_coeffs_8)[int(1)] = _S4403;
    (*v_sh_coeffs_8)[int(2)] = _S4404;
    (*v_sh_coeffs_8)[int(3)] = _S4405;
    (*v_sh_coeffs_8)[int(4)] = _S4406;
    (*v_sh_coeffs_8)[int(5)] = _S4407;
    (*v_sh_coeffs_8)[int(6)] = _S4408;
    (*v_sh_coeffs_8)[int(7)] = _S4409;
    (*v_sh_coeffs_8)[int(8)] = _S4410;
    (*v_sh_coeffs_8)[int(9)] = _S4411;
    (*v_sh_coeffs_8)[int(10)] = _S4412;
    (*v_sh_coeffs_8)[int(11)] = _S4413;
    (*v_sh_coeffs_8)[int(12)] = _S4414;
    (*v_sh_coeffs_8)[int(13)] = _S4415;
    (*v_sh_coeffs_8)[int(14)] = _S4416;
    (*v_sh_coeffs_8)[int(15)] = _S4417;
    (*v_ch_coeffs_3)[int(0)] = _S4399;
    (*v_ch_coeffs_3)[int(1)] = _S4400;
    *v_R_9 = _S4578;
    *v_t_9 = _S4577;
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
    bool _S4583;
    if(u_28 >= 0.0f)
    {
        _S4583 = v_28 >= 0.0f;
    }
    else
    {
        _S4583 = false;
    }
    if(_S4583)
    {
        _S4583 = (u_28 + v_28) <= 1.0f;
    }
    else
    {
        _S4583 = false;
    }
    if(_S4583)
    {
        _S4583 = t_30 >= 0.0f;
    }
    else
    {
        _S4583 = false;
    }
    if(!_S4583)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_28), (v_28)))), ((F32_sqrt((0.5f))) * (1.0f - u_28 - v_28)))) * (2.0f + (F32_sqrt((2.0f))));
    float h_0 = clamp_0(hardness_16.y, 0.0f, 0.95999997854232788f);
    float o_0 = hardness_16.x;
    float _S4584;
    if(opac_0 < 0.0f)
    {
        _S4584 = 0.0f;
    }
    else
    {
        _S4584 = (F32_min((o_0 * (F32_pow((opac_0), (1.0f - h_0)))), (0.99500000476837158f)));
    }
    return _S4584;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_4, float _s_dOut_12)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4585 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4586 = *dpray_d_4;
    float3  v1v0_1 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_1 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_1 = (*dpray_o_4).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S4587 = s_primal_ctx_cross_0(v1v0_1, v2v0_1);
    float3  _S4588 = s_primal_ctx_cross_0(rov0_1, (*dpray_d_4).primal_0);
    float _S4589 = s_primal_ctx_dot_0((*dpray_d_4).primal_0, _S4587);
    float d_1 = 1.0f / _S4589;
    float _S4590 = _S4589 * _S4589;
    float3  _S4591 = - _S4588;
    float _S4592 = s_primal_ctx_dot_0(_S4591, v2v0_1);
    float u_29 = d_1 * _S4592;
    float _S4593 = s_primal_ctx_dot_0(_S4588, v1v0_1);
    float v_29 = d_1 * _S4593;
    float3  _S4594 = - _S4587;
    float t_31 = d_1 * s_primal_ctx_dot_0(_S4594, rov0_1);
    bool _S4595;
    if(u_29 >= 0.0f)
    {
        _S4595 = v_29 >= 0.0f;
    }
    else
    {
        _S4595 = false;
    }
    if(_S4595)
    {
        _S4595 = (u_29 + v_29) <= 1.0f;
    }
    else
    {
        _S4595 = false;
    }
    if(_S4595)
    {
        _S4595 = t_31 >= 0.0f;
    }
    else
    {
        _S4595 = false;
    }
    bool _S4596 = !!_S4595;
    float _S4597;
    float _S4598;
    float _S4599;
    float _S4600;
    float _S4601;
    float _S4602;
    float _S4603;
    float _S4604;
    float _S4605;
    float _S4606;
    if(_S4596)
    {
        float _S4607 = s_primal_ctx_min_0(u_29, v_29);
        float _S4608 = s_primal_ctx_sqrt_0(0.5f);
        float _S4609 = _S4608 * (1.0f - u_29 - v_29);
        float _S4610 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S4607, _S4609) * _S4610;
        float _S4611 = _S4585.primal_0.y;
        float _S4612 = s_primal_ctx_clamp_0(_S4611, 0.0f, 0.95999997854232788f);
        float o_1 = _S4585.primal_0.x;
        bool _S4613 = opac_1 < 0.0f;
        if(_S4613)
        {
            _S4597 = 0.0f;
            _S4598 = 0.0f;
            _S4599 = 0.0f;
        }
        else
        {
            float _S4614 = 1.0f - _S4612;
            float _S4615 = s_primal_ctx_pow_0(opac_1, _S4614);
            _S4597 = o_1 * _S4615;
            _S4598 = _S4615;
            _S4599 = _S4614;
        }
        float _S4616 = _S4598;
        float _S4617 = _S4599;
        _S4595 = _S4613;
        _S4598 = o_1;
        _S4599 = _S4616;
        _S4600 = opac_1;
        _S4601 = _S4617;
        _S4602 = _S4611;
        _S4603 = _S4610;
        _S4604 = _S4607;
        _S4605 = _S4609;
        _S4606 = _S4608;
    }
    else
    {
        _S4595 = false;
        _S4597 = 0.0f;
        _S4598 = 0.0f;
        _S4599 = 0.0f;
        _S4600 = 0.0f;
        _S4601 = 0.0f;
        _S4602 = 0.0f;
        _S4603 = 0.0f;
        _S4604 = 0.0f;
        _S4605 = 0.0f;
        _S4606 = 0.0f;
    }
    float2  _S4618 = make_float2 (0.0f);
    float2  _S4619;
    if(_S4596)
    {
        if(_S4595)
        {
            _S4597 = 0.0f;
            _S4598 = 0.0f;
            _S4599 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S4620;
            (&_S4620)->primal_0 = _S4597;
            (&_S4620)->differential_0 = 0.0f;
            DiffPair_float_0 _S4621;
            (&_S4621)->primal_0 = 0.99500000476837158f;
            (&_S4621)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S4620, &_S4621, _s_dOut_12);
            float _S4622 = _S4598 * _S4620.differential_0;
            float _S4623 = _S4599 * _S4620.differential_0;
            DiffPair_float_0 _S4624;
            (&_S4624)->primal_0 = _S4600;
            (&_S4624)->differential_0 = 0.0f;
            DiffPair_float_0 _S4625;
            (&_S4625)->primal_0 = _S4601;
            (&_S4625)->differential_0 = 0.0f;
            s_bwd_prop_pow_0(&_S4624, &_S4625, _S4622);
            float _S4626 = - _S4625.differential_0;
            _S4597 = _S4623;
            _S4598 = _S4626;
            _S4599 = _S4624.differential_0;
        }
        DiffPair_float_0 _S4627;
        (&_S4627)->primal_0 = _S4602;
        (&_S4627)->differential_0 = 0.0f;
        DiffPair_float_0 _S4628;
        (&_S4628)->primal_0 = 0.0f;
        (&_S4628)->differential_0 = 0.0f;
        DiffPair_float_0 _S4629;
        (&_S4629)->primal_0 = 0.95999997854232788f;
        (&_S4629)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S4627, &_S4628, &_S4629, _S4598);
        float _S4630 = _S4603 * _S4599;
        DiffPair_float_0 _S4631;
        (&_S4631)->primal_0 = _S4604;
        (&_S4631)->differential_0 = 0.0f;
        DiffPair_float_0 _S4632;
        (&_S4632)->primal_0 = _S4605;
        (&_S4632)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4631, &_S4632, _S4630);
        float _S4633 = - (_S4606 * _S4632.differential_0);
        DiffPair_float_0 _S4634;
        (&_S4634)->primal_0 = u_29;
        (&_S4634)->differential_0 = 0.0f;
        DiffPair_float_0 _S4635;
        (&_S4635)->primal_0 = v_29;
        (&_S4635)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4634, &_S4635, _S4631.differential_0);
        float2  _S4636 = make_float2 (_S4597, _S4627.differential_0);
        float _S4637 = _S4633 + _S4635.differential_0;
        _S4597 = _S4633 + _S4634.differential_0;
        _S4598 = _S4637;
        _S4619 = _S4636;
    }
    else
    {
        _S4597 = 0.0f;
        _S4598 = 0.0f;
        _S4619 = _S4618;
    }
    float3  _S4638 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4639;
    (&_S4639)->primal_0 = _S4594;
    (&_S4639)->differential_0 = _S4638;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4640;
    (&_S4640)->primal_0 = rov0_1;
    (&_S4640)->differential_0 = _S4638;
    s_bwd_prop_dot_0(&_S4639, &_S4640, 0.0f);
    float3  _S4641 = - _S4639.differential_0;
    float _S4642 = d_1 * _S4598;
    float _S4643 = _S4593 * _S4598;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4644;
    (&_S4644)->primal_0 = _S4588;
    (&_S4644)->differential_0 = _S4638;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4645;
    (&_S4645)->primal_0 = v1v0_1;
    (&_S4645)->differential_0 = _S4638;
    s_bwd_prop_dot_0(&_S4644, &_S4645, _S4642);
    float _S4646 = d_1 * _S4597;
    float _S4647 = _S4592 * _S4597;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4648;
    (&_S4648)->primal_0 = _S4591;
    (&_S4648)->differential_0 = _S4638;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4649;
    (&_S4649)->primal_0 = v2v0_1;
    (&_S4649)->differential_0 = _S4638;
    s_bwd_prop_dot_0(&_S4648, &_S4649, _S4646);
    float3  _S4650 = - _S4648.differential_0;
    float _S4651 = - ((_S4643 + _S4647) / _S4590);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4652;
    (&_S4652)->primal_0 = _S4586.primal_0;
    (&_S4652)->differential_0 = _S4638;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4653;
    (&_S4653)->primal_0 = _S4587;
    (&_S4653)->differential_0 = _S4638;
    s_bwd_prop_dot_0(&_S4652, &_S4653, _S4651);
    float3  _S4654 = _S4644.differential_0 + _S4650;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4655;
    (&_S4655)->primal_0 = rov0_1;
    (&_S4655)->differential_0 = _S4638;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4656;
    (&_S4656)->primal_0 = _S4586.primal_0;
    (&_S4656)->differential_0 = _S4638;
    s_bwd_prop_cross_0(&_S4655, &_S4656, _S4654);
    float3  _S4657 = _S4641 + _S4653.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4658;
    (&_S4658)->primal_0 = v1v0_1;
    (&_S4658)->differential_0 = _S4638;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4659;
    (&_S4659)->primal_0 = v2v0_1;
    (&_S4659)->differential_0 = _S4638;
    s_bwd_prop_cross_0(&_S4658, &_S4659, _S4657);
    float3  _S4660 = _S4640.differential_0 + _S4655.differential_0;
    float3  _S4661 = _S4649.differential_0 + _S4659.differential_0;
    float3  _S4662 = _S4645.differential_0 + _S4658.differential_0;
    float3  _S4663 = - _S4660 + - _S4661 + - _S4662;
    float3  _S4664 = _S4652.differential_0 + _S4656.differential_0;
    dpray_d_4->primal_0 = (*dpray_d_4).primal_0;
    dpray_d_4->differential_0 = _S4664;
    dpray_o_4->primal_0 = (*dpray_o_4).primal_0;
    dpray_o_4->differential_0 = _S4660;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S4619;
    FixedArray<float3 , 3>  _S4665;
    _S4665[int(0)] = _S4638;
    _S4665[int(1)] = _S4638;
    _S4665[int(2)] = _S4638;
    _S4665[int(2)] = _S4661;
    _S4665[int(0)] = _S4663;
    _S4665[int(1)] = _S4662;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S4665;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4666, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4667, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4668, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4669, float _S4670)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S4666, _S4667, _S4668, _S4669, _S4670);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_5, float2  hardness_17, float3  ray_o_6, float3  ray_d_6, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S4671 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4672 = { _S4671, _S4671, _S4671 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_5;
    (&dp_verts_0)->differential_0 = _S4672;
    float2  _S4673 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S4673;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_6;
    (&dp_ray_o_2)->differential_0 = _S4671;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_6;
    (&dp_ray_d_2)->differential_0 = _S4671;
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
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float3  dpcolor_1, float dpdepth_2)
{
    float3  v1v0_3 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_3 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_3 = (*dpray_o_5).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S4674 = s_primal_ctx_cross_0(v1v0_3, v2v0_3);
    float3  _S4675 = s_primal_ctx_cross_0(rov0_3, (*dpray_d_5).primal_0);
    float _S4676 = s_primal_ctx_dot_0((*dpray_d_5).primal_0, _S4674);
    float d_3 = 1.0f / _S4676;
    float _S4677 = _S4676 * _S4676;
    float3  _S4678 = - _S4675;
    float _S4679 = s_primal_ctx_dot_0(_S4678, v2v0_3);
    float u_31 = d_3 * _S4679;
    float _S4680 = s_primal_ctx_dot_0(_S4675, v1v0_3);
    float v_31 = d_3 * _S4680;
    float3  _S4681 = - _S4674;
    float3  _S4682 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S4683 = make_float3 (v_31) * dpcolor_1;
    float3  _S4684 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S4685 = make_float3 (u_31) * dpcolor_1;
    float3  _S4686 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S4687 = make_float3 (1.0f - u_31 - v_31) * dpcolor_1;
    float _S4688 = - (_S4686.x + _S4686.y + _S4686.z);
    float _S4689 = d_3 * dpdepth_2;
    float _S4690 = s_primal_ctx_dot_0(_S4681, rov0_3) * dpdepth_2;
    float3  _S4691 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4692;
    (&_S4692)->primal_0 = _S4681;
    (&_S4692)->differential_0 = _S4691;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4693;
    (&_S4693)->primal_0 = rov0_3;
    (&_S4693)->differential_0 = _S4691;
    s_bwd_prop_dot_0(&_S4692, &_S4693, _S4689);
    float3  _S4694 = - _S4692.differential_0;
    float _S4695 = _S4688 + _S4682.x + _S4682.y + _S4682.z;
    float _S4696 = d_3 * _S4695;
    float _S4697 = _S4680 * _S4695;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4698;
    (&_S4698)->primal_0 = _S4675;
    (&_S4698)->differential_0 = _S4691;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4699;
    (&_S4699)->primal_0 = v1v0_3;
    (&_S4699)->differential_0 = _S4691;
    s_bwd_prop_dot_0(&_S4698, &_S4699, _S4696);
    float _S4700 = _S4688 + _S4684.x + _S4684.y + _S4684.z;
    float _S4701 = d_3 * _S4700;
    float _S4702 = _S4679 * _S4700;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4703;
    (&_S4703)->primal_0 = _S4678;
    (&_S4703)->differential_0 = _S4691;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4704;
    (&_S4704)->primal_0 = v2v0_3;
    (&_S4704)->differential_0 = _S4691;
    s_bwd_prop_dot_0(&_S4703, &_S4704, _S4701);
    float3  _S4705 = - _S4703.differential_0;
    float _S4706 = - ((_S4690 + _S4697 + _S4702) / _S4677);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4707;
    (&_S4707)->primal_0 = (*dpray_d_5).primal_0;
    (&_S4707)->differential_0 = _S4691;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4708;
    (&_S4708)->primal_0 = _S4674;
    (&_S4708)->differential_0 = _S4691;
    s_bwd_prop_dot_0(&_S4707, &_S4708, _S4706);
    float3  _S4709 = _S4698.differential_0 + _S4705;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4710;
    (&_S4710)->primal_0 = rov0_3;
    (&_S4710)->differential_0 = _S4691;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4711;
    (&_S4711)->primal_0 = (*dpray_d_5).primal_0;
    (&_S4711)->differential_0 = _S4691;
    s_bwd_prop_cross_0(&_S4710, &_S4711, _S4709);
    float3  _S4712 = _S4694 + _S4708.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4713;
    (&_S4713)->primal_0 = v1v0_3;
    (&_S4713)->differential_0 = _S4691;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4714;
    (&_S4714)->primal_0 = v2v0_3;
    (&_S4714)->differential_0 = _S4691;
    s_bwd_prop_cross_0(&_S4713, &_S4714, _S4712);
    float3  _S4715 = _S4693.differential_0 + _S4710.differential_0;
    float3  _S4716 = _S4704.differential_0 + _S4714.differential_0;
    float3  _S4717 = _S4699.differential_0 + _S4713.differential_0;
    float3  _S4718 = - _S4715 + - _S4716 + - _S4717;
    float3  _S4719 = _S4707.differential_0 + _S4711.differential_0;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S4719;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S4715;
    FixedArray<float3 , 3>  _S4720;
    _S4720[int(0)] = _S4691;
    _S4720[int(1)] = _S4691;
    _S4720[int(2)] = _S4691;
    _S4720[int(2)] = _S4683;
    _S4720[int(1)] = _S4685;
    _S4720[int(0)] = _S4687;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S4720;
    FixedArray<float3 , 3>  _S4721;
    _S4721[int(0)] = _S4691;
    _S4721[int(1)] = _S4691;
    _S4721[int(2)] = _S4691;
    _S4721[int(2)] = _S4716;
    _S4721[int(0)] = _S4718;
    _S4721[int(1)] = _S4717;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S4721;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4722, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4723, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4724, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4725, float3  _S4726, float _S4727)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S4722, _S4723, _S4724, _S4725, _S4726, _S4727);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_7, FixedArray<float3 , 3>  * rgbs_5, float3  ray_o_8, float3  ray_d_8, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S4728 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4729 = { _S4728, _S4728, _S4728 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_7;
    (&dp_verts_1)->differential_0 = _S4729;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_5;
    (&dp_rgbs_0)->differential_0 = _S4729;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_8;
    (&dp_ray_o_3)->differential_0 = _S4728;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_8;
    (&dp_ray_d_3)->differential_0 = _S4728;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

