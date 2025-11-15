#pragma once

#include "slang.cuh"

inline __device__ Matrix<float, 3, 3>  transpose_0(Matrix<float, 3, 3>  x_0)
{
    Matrix<float, 3, 3>  result_0;
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
            *_slang_vector_get_element_ptr(((&result_0)->rows + (r_0)), c_0) = _slang_vector_get_element(x_0.rows[c_0], r_0);
            c_0 = c_0 + int(1);
        }
        r_0 = r_0 + int(1);
    }
    return result_0;
}

inline __device__ Matrix<float, 3, 2>  transpose_1(Matrix<float, 2, 3>  x_1)
{
    Matrix<float, 3, 2>  result_1;
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
            *_slang_vector_get_element_ptr(((&result_1)->rows + (r_1)), c_1) = _slang_vector_get_element(x_1.rows[c_1], r_1);
            c_1 = c_1 + int(1);
        }
        r_1 = r_1 + int(1);
    }
    return result_1;
}

inline __device__ Matrix<float, 2, 3>  transpose_2(Matrix<float, 3, 2>  x_2)
{
    Matrix<float, 2, 3>  result_2;
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
            *_slang_vector_get_element_ptr(((&result_2)->rows + (r_2)), c_2) = _slang_vector_get_element(x_2.rows[c_2], r_2);
            c_2 = c_2 + int(1);
        }
        r_2 = r_2 + int(1);
    }
    return result_2;
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
    float _S1 = (*left_0).primal_0.rows[int(0)].x * dOut_0.x;
    Matrix<float, 3, 3>  left_d_result_0;
    *&(((&left_d_result_0)->rows + (int(0)))->x) = (*right_0).primal_0.x * dOut_0.x;
    float sum_0 = _S1 + (*left_0).primal_0.rows[int(1)].x * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->x) = (*right_0).primal_0.x * dOut_0.y;
    float sum_1 = sum_0 + (*left_0).primal_0.rows[int(2)].x * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->x) = (*right_0).primal_0.x * dOut_0.z;
    float3  right_d_result_0;
    *&((&right_d_result_0)->x) = sum_1;
    float _S2 = (*left_0).primal_0.rows[int(0)].y * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->y) = (*right_0).primal_0.y * dOut_0.x;
    float sum_2 = _S2 + (*left_0).primal_0.rows[int(1)].y * dOut_0.y;
    *&(((&left_d_result_0)->rows + (int(1)))->y) = (*right_0).primal_0.y * dOut_0.y;
    float sum_3 = sum_2 + (*left_0).primal_0.rows[int(2)].y * dOut_0.z;
    *&(((&left_d_result_0)->rows + (int(2)))->y) = (*right_0).primal_0.y * dOut_0.z;
    *&((&right_d_result_0)->y) = sum_3;
    float _S3 = (*left_0).primal_0.rows[int(0)].z * dOut_0.x;
    *&(((&left_d_result_0)->rows + (int(0)))->z) = (*right_0).primal_0.z * dOut_0.x;
    float sum_4 = _S3 + (*left_0).primal_0.rows[int(1)].z * dOut_0.y;
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
    float3  result_3;
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
        *_slang_vector_get_element_ptr(&result_3, i_0) = sum_6;
        i_0 = i_0 + int(1);
    }
    return result_3;
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
    Matrix<float, 3, 3>  result_4;
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
            *_slang_vector_get_element_ptr(((&result_4)->rows + (r_3)), c_3) = sum_8;
            c_3 = c_3 + int(1);
        }
        r_3 = r_3 + int(1);
    }
    return result_4;
}

inline __device__ Matrix<float, 2, 3>  mul_5(Matrix<float, 2, 3>  left_6, Matrix<float, 3, 3>  right_6)
{
    Matrix<float, 2, 3>  result_5;
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
            *_slang_vector_get_element_ptr(((&result_5)->rows + (r_4)), c_4) = sum_10;
            c_4 = c_4 + int(1);
        }
        r_4 = r_4 + int(1);
    }
    return result_5;
}

inline __device__ Matrix<float, 2, 2>  mul_6(Matrix<float, 2, 3>  left_7, Matrix<float, 3, 2>  right_7)
{
    Matrix<float, 2, 2>  result_6;
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
            *_slang_vector_get_element_ptr(((&result_6)->rows + (r_5)), c_5) = sum_12;
            c_5 = c_5 + int(1);
        }
        r_5 = r_5 + int(1);
    }
    return result_6;
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

inline __device__ float determinant_0(Matrix<float, 2, 2>  m_1)
{
    return m_1.rows[int(0)].x * m_1.rows[int(1)].y - m_1.rows[int(0)].y * m_1.rows[int(1)].x;
}

struct DiffPair_float_0
{
    float primal_0;
    float differential_0;
};

inline __device__ void _d_min_0(DiffPair_float_0 * dpx_0, DiffPair_float_0 * dpy_0, float dOut_4)
{
    DiffPair_float_0 _S4 = *dpx_0;
    float _S5;
    if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
    {
        _S5 = dOut_4;
    }
    else
    {
        if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
        {
            _S5 = 0.0f;
        }
        else
        {
            _S5 = 0.5f * dOut_4;
        }
    }
    dpx_0->primal_0 = _S4.primal_0;
    dpx_0->differential_0 = _S5;
    DiffPair_float_0 _S6 = *dpy_0;
    if(((*dpy_0).primal_0) < (_S4.primal_0))
    {
        _S5 = dOut_4;
    }
    else
    {
        if(((*dpy_0).primal_0) > ((*dpx_0).primal_0))
        {
            _S5 = 0.0f;
        }
        else
        {
            _S5 = 0.5f * dOut_4;
        }
    }
    dpy_0->primal_0 = _S6.primal_0;
    dpy_0->differential_0 = _S5;
    return;
}

inline __device__ bool is_valid_distortion(float2  uv_0, FixedArray<float, 10>  * dist_coeffs_0)
{
    float2  _S7 = make_float2 (0.0f);
    float2  seed_0 = _S7;
    *&((&seed_0)->x) = 1.0f;
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float _S8 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
    float _S9 = (*dist_coeffs_0)[int(1)] + r2_0 * _S8;
    float _S10 = (*dist_coeffs_0)[int(0)] + r2_0 * _S9;
    float2  _S11 = make_float2 (1.0f + r2_0 * _S10);
    float _S12 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S13 = _S12 * u_0;
    float _S14 = 2.0f * u_0;
    float _S15 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S16 = _S15 * u_0;
    float _S17 = 2.0f * v_0;
    float2  _S18 = seed_0 + make_float2 ((*dist_coeffs_0)[int(8)] * seed_0.x, (*dist_coeffs_0)[int(9)] * seed_0.x);
    float2  _S19 = uv_0 * _S18;
    float _S20 = (*dist_coeffs_0)[int(4)] * _S18.y;
    float _S21 = (*dist_coeffs_0)[int(5)] * _S18.x;
    float _S22 = _S19.x + _S19.y;
    float _S23 = r2_0 * _S22;
    float _S24 = r2_0 * _S23;
    float _S25 = (*dist_coeffs_0)[int(7)] * _S18.y + _S20 + (*dist_coeffs_0)[int(6)] * _S18.x + _S21 + _S10 * _S22 + _S9 * _S23 + _S8 * _S24 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S24);
    float _S26 = v_0 * _S25;
    float _S27 = u_0 * _S25;
    Matrix<float, 2, 2>  J_0;
    J_0[int(0)] = _S11 * _S18 + make_float2 (_S15 * (v_0 * _S18.y) + _S14 * _S21 + 2.0f * (u_0 * _S21) + _S12 * (v_0 * _S18.x) + _S27 + _S27, _S17 * _S20 + 2.0f * (v_0 * _S20) + _S16 * _S18.y + _S13 * _S18.x + _S26 + _S26);
    float2  seed_1 = _S7;
    *&((&seed_1)->y) = 1.0f;
    float2  _S28 = seed_1 + make_float2 ((*dist_coeffs_0)[int(8)] * seed_1.x, (*dist_coeffs_0)[int(9)] * seed_1.x);
    float2  _S29 = uv_0 * _S28;
    float _S30 = (*dist_coeffs_0)[int(4)] * _S28.y;
    float _S31 = (*dist_coeffs_0)[int(5)] * _S28.x;
    float _S32 = _S29.x + _S29.y;
    float _S33 = r2_0 * _S32;
    float _S34 = r2_0 * _S33;
    float _S35 = (*dist_coeffs_0)[int(7)] * _S28.y + _S30 + (*dist_coeffs_0)[int(6)] * _S28.x + _S31 + _S10 * _S32 + _S9 * _S33 + _S8 * _S34 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S34);
    float _S36 = v_0 * _S35;
    float _S37 = u_0 * _S35;
    J_0[int(1)] = _S11 * _S28 + make_float2 (_S15 * (v_0 * _S28.y) + _S14 * _S31 + 2.0f * (u_0 * _S31) + _S12 * (v_0 * _S28.x) + _S37 + _S37, _S17 * _S30 + 2.0f * (v_0 * _S30) + _S16 * _S28.y + _S13 * _S28.x + _S36 + _S36);
    return (F32_min((determinant_0(J_0)), ((F32_min((J_0.rows[int(0)].x), (J_0.rows[int(1)].y)))))) > 0.0f;
}

inline __device__ void _d_max_0(DiffPair_float_0 * dpx_1, DiffPair_float_0 * dpy_1, float dOut_5)
{
    DiffPair_float_0 _S38 = *dpx_1;
    float _S39;
    if(((*dpx_1).primal_0) > ((*dpy_1).primal_0))
    {
        _S39 = dOut_5;
    }
    else
    {
        if(((*dpx_1).primal_0) < ((*dpy_1).primal_0))
        {
            _S39 = 0.0f;
        }
        else
        {
            _S39 = 0.5f * dOut_5;
        }
    }
    dpx_1->primal_0 = _S38.primal_0;
    dpx_1->differential_0 = _S39;
    DiffPair_float_0 _S40 = *dpy_1;
    if(((*dpy_1).primal_0) > (_S38.primal_0))
    {
        _S39 = dOut_5;
    }
    else
    {
        if(((*dpy_1).primal_0) < ((*dpx_1).primal_0))
        {
            _S39 = 0.0f;
        }
        else
        {
            _S39 = 0.5f * dOut_5;
        }
    }
    dpy_1->primal_0 = _S40.primal_0;
    dpy_1->differential_0 = _S39;
    return;
}

inline __device__ void _d_sqrt_0(DiffPair_float_0 * dpx_2, float dOut_6)
{
    float _S41 = 0.5f / (F32_sqrt(((F32_max((1.00000001168609742e-07f), ((*dpx_2).primal_0)))))) * dOut_6;
    dpx_2->primal_0 = (*dpx_2).primal_0;
    dpx_2->differential_0 = _S41;
    return;
}

inline __device__ void _d_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_2, float dOut_7)
{
    float3  x_d_result_0;
    *&((&x_d_result_0)->x) = (*dpy_2).primal_0.x * dOut_7;
    float3  y_d_result_0;
    *&((&y_d_result_0)->x) = (*dpx_3).primal_0.x * dOut_7;
    *&((&x_d_result_0)->y) = (*dpy_2).primal_0.y * dOut_7;
    *&((&y_d_result_0)->y) = (*dpx_3).primal_0.y * dOut_7;
    *&((&x_d_result_0)->z) = (*dpy_2).primal_0.z * dOut_7;
    *&((&y_d_result_0)->z) = (*dpx_3).primal_0.z * dOut_7;
    dpx_3->primal_0 = (*dpx_3).primal_0;
    dpx_3->differential_0 = x_d_result_0;
    dpy_2->primal_0 = (*dpy_2).primal_0;
    dpy_2->differential_0 = y_d_result_0;
    return;
}

struct DiffPair_vectorx3Cfloatx2C2x3E_0
{
    float2  primal_0;
    float2  differential_0;
};

inline __device__ void _d_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_4, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpy_3, float dOut_8)
{
    float2  x_d_result_1;
    *&((&x_d_result_1)->x) = (*dpy_3).primal_0.x * dOut_8;
    float2  y_d_result_1;
    *&((&y_d_result_1)->x) = (*dpx_4).primal_0.x * dOut_8;
    *&((&x_d_result_1)->y) = (*dpy_3).primal_0.y * dOut_8;
    *&((&y_d_result_1)->y) = (*dpx_4).primal_0.y * dOut_8;
    dpx_4->primal_0 = (*dpx_4).primal_0;
    dpx_4->differential_0 = x_d_result_1;
    dpy_3->primal_0 = (*dpy_3).primal_0;
    dpy_3->differential_0 = y_d_result_1;
    return;
}

inline __device__ float dot_0(float3  x_6, float3  y_0)
{
    int i_4 = int(0);
    float result_7 = 0.0f;
    for(;;)
    {
        if(i_4 < int(3))
        {
        }
        else
        {
            break;
        }
        float result_8 = result_7 + _slang_vector_get_element(x_6, i_4) * _slang_vector_get_element(y_0, i_4);
        i_4 = i_4 + int(1);
        result_7 = result_8;
    }
    return result_7;
}

inline __device__ float dot_1(float2  x_7, float2  y_1)
{
    int i_5 = int(0);
    float result_9 = 0.0f;
    for(;;)
    {
        if(i_5 < int(2))
        {
        }
        else
        {
            break;
        }
        float result_10 = result_9 + _slang_vector_get_element(x_7, i_5) * _slang_vector_get_element(y_1, i_5);
        i_5 = i_5 + int(1);
        result_9 = result_10;
    }
    return result_9;
}

inline __device__ float length_0(float2  x_8)
{
    return (F32_sqrt((dot_1(x_8, x_8))));
}

inline __device__ float length_1(float3  x_9)
{
    return (F32_sqrt((dot_0(x_9, x_9))));
}

inline __device__ float2  distort_point(float2  uv_1, bool is_fisheye_0, FixedArray<float, 10>  * dist_coeffs_1)
{
    float2  _S42;
    if(is_fisheye_0)
    {
        float r_6 = length_0(uv_1);
        float theta_0 = (F32_atan((r_6)));
        float _S43;
        if(r_6 < 0.00100000004749745f)
        {
            _S43 = 1.0f - r_6 * r_6 / 3.0f;
        }
        else
        {
            _S43 = theta_0 / r_6;
        }
        _S42 = uv_1 * make_float2 (_S43);
    }
    else
    {
        _S42 = uv_1;
    }
    float u_1 = _S42.x;
    float v_1 = _S42.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S44 = _S42 * make_float2 (1.0f + r2_1 * ((*dist_coeffs_1)[int(0)] + r2_1 * ((*dist_coeffs_1)[int(1)] + r2_1 * ((*dist_coeffs_1)[int(2)] + r2_1 * (*dist_coeffs_1)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_1)[int(4)] * u_1 * v_1 + (*dist_coeffs_1)[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + (*dist_coeffs_1)[int(6)] * r2_1, 2.0f * (*dist_coeffs_1)[int(5)] * u_1 * v_1 + (*dist_coeffs_1)[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + (*dist_coeffs_1)[int(7)] * r2_1);
    return _S44 + make_float2 ((*dist_coeffs_1)[int(8)] * _S44.x + (*dist_coeffs_1)[int(9)] * _S44.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_2, FixedArray<float, 10>  * dist_coeffs_2, int maxiter_0, float2  * uv_undist_0)
{
    int i_6 = int(0);
    float2  q_0 = uv_2;
    for(;;)
    {
        if(i_6 < maxiter_0)
        {
        }
        else
        {
            break;
        }
        float u_2 = q_0.x;
        float v_2 = q_0.y;
        float r2_2 = u_2 * u_2 + v_2 * v_2;
        float _S45 = (*dist_coeffs_2)[int(2)] + r2_2 * (*dist_coeffs_2)[int(3)];
        float _S46 = (*dist_coeffs_2)[int(1)] + r2_2 * _S45;
        float _S47 = (*dist_coeffs_2)[int(0)] + r2_2 * _S46;
        float radial_0 = 1.0f + r2_2 * _S47;
        float _S48 = 2.0f * (*dist_coeffs_2)[int(4)];
        float _S49 = _S48 * u_2;
        float _S50 = 2.0f * u_2;
        float _S51 = 2.0f * (*dist_coeffs_2)[int(5)];
        float _S52 = _S51 * u_2;
        float _S53 = 2.0f * v_2;
        float2  _S54 = q_0 * make_float2 (radial_0) + make_float2 (_S49 * v_2 + (*dist_coeffs_2)[int(5)] * (r2_2 + _S50 * u_2) + (*dist_coeffs_2)[int(6)] * r2_2, _S52 * v_2 + (*dist_coeffs_2)[int(4)] * (r2_2 + _S53 * v_2) + (*dist_coeffs_2)[int(7)] * r2_2);
        float2  r_7 = _S54 + make_float2 ((*dist_coeffs_2)[int(8)] * _S54.x + (*dist_coeffs_2)[int(9)] * _S54.y, 0.0f) - uv_2;
        float2  _S55 = make_float2 (0.0f);
        float2  seed_2 = _S55;
        *&((&seed_2)->x) = 1.0f;
        float2  _S56 = make_float2 (radial_0);
        float2  _S57 = seed_2 + make_float2 ((*dist_coeffs_2)[int(8)] * seed_2.x, (*dist_coeffs_2)[int(9)] * seed_2.x);
        float2  _S58 = q_0 * _S57;
        float _S59 = (*dist_coeffs_2)[int(4)] * _S57.y;
        float _S60 = (*dist_coeffs_2)[int(5)] * _S57.x;
        float _S61 = _S58.x + _S58.y;
        float _S62 = r2_2 * _S61;
        float _S63 = r2_2 * _S62;
        float _S64 = (*dist_coeffs_2)[int(7)] * _S57.y + _S59 + (*dist_coeffs_2)[int(6)] * _S57.x + _S60 + _S47 * _S61 + _S46 * _S62 + _S45 * _S63 + (*dist_coeffs_2)[int(3)] * (r2_2 * _S63);
        float _S65 = v_2 * _S64;
        float _S66 = u_2 * _S64;
        Matrix<float, 2, 2>  J_1;
        J_1[int(0)] = _S56 * _S57 + make_float2 (_S51 * (v_2 * _S57.y) + _S50 * _S60 + 2.0f * (u_2 * _S60) + _S48 * (v_2 * _S57.x) + _S66 + _S66, _S53 * _S59 + 2.0f * (v_2 * _S59) + _S52 * _S57.y + _S49 * _S57.x + _S65 + _S65);
        float2  seed_3 = _S55;
        *&((&seed_3)->y) = 1.0f;
        float2  _S67 = seed_3 + make_float2 ((*dist_coeffs_2)[int(8)] * seed_3.x, (*dist_coeffs_2)[int(9)] * seed_3.x);
        float2  _S68 = q_0 * _S67;
        float _S69 = (*dist_coeffs_2)[int(4)] * _S67.y;
        float _S70 = (*dist_coeffs_2)[int(5)] * _S67.x;
        float _S71 = _S68.x + _S68.y;
        float _S72 = r2_2 * _S71;
        float _S73 = r2_2 * _S72;
        float _S74 = (*dist_coeffs_2)[int(7)] * _S67.y + _S69 + (*dist_coeffs_2)[int(6)] * _S67.x + _S70 + _S47 * _S71 + _S46 * _S72 + _S45 * _S73 + (*dist_coeffs_2)[int(3)] * (r2_2 * _S73);
        float _S75 = v_2 * _S74;
        float _S76 = u_2 * _S74;
        J_1[int(1)] = _S56 * _S67 + make_float2 (_S51 * (v_2 * _S67.y) + _S50 * _S70 + 2.0f * (u_2 * _S70) + _S48 * (v_2 * _S67.x) + _S76 + _S76, _S53 * _S69 + 2.0f * (v_2 * _S69) + _S52 * _S67.y + _S49 * _S67.x + _S75 + _S75);
        float inv_det_0 = 1.0f / (J_1.rows[int(0)].x * J_1.rows[int(1)].y - J_1.rows[int(0)].y * J_1.rows[int(1)].x);
        float _S77 = r_7.x;
        float _S78 = r_7.y;
        float2  q_1 = q_0 - make_float2 ((_S77 * J_1.rows[int(1)].y - _S78 * J_1.rows[int(0)].y) * inv_det_0, (- _S77 * J_1.rows[int(1)].x + _S78 * J_1.rows[int(0)].x) * inv_det_0);
        i_6 = i_6 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    bool _S79 = is_valid_distortion(q_0, dist_coeffs_2);
    bool _S80;
    if(_S79)
    {
        float u_3 = (*uv_undist_0).x;
        float v_3 = (*uv_undist_0).y;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float2  _S81 = *uv_undist_0 * make_float2 (1.0f + r2_3 * ((*dist_coeffs_2)[int(0)] + r2_3 * ((*dist_coeffs_2)[int(1)] + r2_3 * ((*dist_coeffs_2)[int(2)] + r2_3 * (*dist_coeffs_2)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_2)[int(4)] * u_3 * v_3 + (*dist_coeffs_2)[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + (*dist_coeffs_2)[int(6)] * r2_3, 2.0f * (*dist_coeffs_2)[int(5)] * u_3 * v_3 + (*dist_coeffs_2)[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + (*dist_coeffs_2)[int(7)] * r2_3);
        _S80 = (length_0(_S81 + make_float2 ((*dist_coeffs_2)[int(8)] * _S81.x + (*dist_coeffs_2)[int(9)] * _S81.y, 0.0f) - uv_2)) < 0.00999999977648258f;
    }
    else
    {
        _S80 = false;
    }
    return _S80;
}

inline __device__ bool undistort_point(float2  uv_3, bool is_fisheye_1, FixedArray<float, 10>  * dist_coeffs_3, float2  * uv_undist_1)
{
    float2  _S82 = uv_3;
    bool _S83 = undistort_point_0(uv_3, dist_coeffs_3, int(8), &_S82);
    if(!_S83)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S84 = _S82;
        float theta_1 = length_0(_S82);
        float _S85;
        if(theta_1 < 0.00100000004749745f)
        {
            _S85 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S85 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S86 = make_float3 ((_S84 * make_float2 (_S85)).x, (_S84 * make_float2 (_S85)).y, (F32_cos((theta_1))));
        raydir_0 = _S86;
    }
    else
    {
        raydir_0 = make_float3 (_S82.x, _S82.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_4, bool is_fisheye_2, FixedArray<float, 10>  * dist_coeffs_4, float3  * raydir_1)
{
    float2  _S87 = uv_4;
    bool _S88 = undistort_point_0(uv_4, dist_coeffs_4, int(8), &_S87);
    if(!_S88)
    {
        return false;
    }
    float3  _S89;
    if(is_fisheye_2)
    {
        float2  _S90 = _S87;
        float theta_2 = length_0(_S87);
        float _S91;
        if(theta_2 < 0.00100000004749745f)
        {
            _S91 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S91 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S92 = make_float3 ((_S90 * make_float2 (_S91)).x, (_S90 * make_float2 (_S91)).y, (F32_cos((theta_2))));
        _S89 = _S92;
    }
    else
    {
        _S89 = make_float3 (_S87.x, _S87.y, 1.0f);
    }
    *raydir_1 = _S89;
    return true;
}

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_8, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_8, float3  dOut_9)
{
    float _S93 = (*right_8).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 3, 3>  right_d_result_4;
    *&(((&right_d_result_4)->rows + (int(0)))->x) = (*left_8).primal_0.x * dOut_9.x;
    float sum_14 = _S93 + (*right_8).primal_0.rows[int(0)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(0)))->y) = (*left_8).primal_0.x * dOut_9.y;
    float sum_15 = sum_14 + (*right_8).primal_0.rows[int(0)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(0)))->z) = (*left_8).primal_0.x * dOut_9.z;
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = sum_15;
    float _S94 = (*right_8).primal_0.rows[int(1)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(1)))->x) = (*left_8).primal_0.y * dOut_9.x;
    float sum_16 = _S94 + (*right_8).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(1)))->y) = (*left_8).primal_0.y * dOut_9.y;
    float sum_17 = sum_16 + (*right_8).primal_0.rows[int(1)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(1)))->z) = (*left_8).primal_0.y * dOut_9.z;
    *&((&left_d_result_4)->y) = sum_17;
    float _S95 = (*right_8).primal_0.rows[int(2)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(2)))->x) = (*left_8).primal_0.z * dOut_9.x;
    float sum_18 = _S95 + (*right_8).primal_0.rows[int(2)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(2)))->y) = (*left_8).primal_0.z * dOut_9.y;
    float sum_19 = sum_18 + (*right_8).primal_0.rows[int(2)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(2)))->z) = (*left_8).primal_0.z * dOut_9.z;
    *&((&left_d_result_4)->z) = sum_19;
    left_8->primal_0 = (*left_8).primal_0;
    left_8->differential_0 = left_d_result_4;
    right_8->primal_0 = (*right_8).primal_0;
    right_8->differential_0 = right_d_result_4;
    return;
}

inline __device__ float3  mul_7(float3  left_9, Matrix<float, 3, 3>  right_9)
{
    float3  result_11;
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
        *_slang_vector_get_element_ptr(&result_11, j_1) = sum_20;
        j_1 = j_1 + int(1);
    }
    return result_11;
}

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float2  normalize_1(float2  x_11)
{
    return x_11 / make_float2 (length_0(x_11));
}

inline __device__ bool generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_5, bool is_fisheye_3, FixedArray<float, 10>  * dist_coeffs_5, float3  * ray_o_0, float3  * ray_d_0)
{
    float2  _S96 = uv_5;
    *ray_o_0 = - mul_7(t_1, R_2);
    bool _S97 = undistort_point_0(uv_5, dist_coeffs_5, int(8), &_S96);
    if(!_S97)
    {
        return false;
    }
    float3  raydir_2;
    if(is_fisheye_3)
    {
        float2  _S98 = _S96;
        float theta_3 = length_0(_S96);
        float _S99;
        if(theta_3 < 0.00100000004749745f)
        {
            _S99 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S99 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S100 = make_float3 ((_S98 * make_float2 (_S99)).x, (_S98 * make_float2 (_S99)).y, (F32_cos((theta_3))));
        raydir_2 = _S100;
    }
    else
    {
        raydir_2 = make_float3 (_S96.x, _S96.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_7(raydir_2, R_2));
    return true;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S101;
    bool _S102;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S103, Matrix<float, 3, 3>  _S104)
{
    return mul_7(_S103, _S104);
}

inline __device__ float s_primal_ctx_sin_0(float _S105)
{
    return (F32_sin((_S105)));
}

inline __device__ float s_primal_ctx_cos_0(float _S106)
{
    return (F32_cos((_S106)));
}

inline __device__ bool s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_6, bool is_fisheye_4, FixedArray<float, 10>  * dist_coeffs_6, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S101 = make_float2 (0.0f);
    _s_diff_ctx_0->_S102 = false;
    float3  _S107 = make_float3 (0.0f);
    float3  _S108 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    float2  _S109 = uv_6;
    bool _S110 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S109);
    _s_diff_ctx_0->_S101 = _S109;
    _s_diff_ctx_0->_S102 = _S110;
    float2  _S111 = _S109;
    float3  raydir_3;
    bool _S112;
    if(!!_S110)
    {
        if(is_fisheye_4)
        {
            float _S113 = length_0(_S111);
            float _S114;
            if(_S113 < 0.00100000004749745f)
            {
                _S114 = 1.0f - _S113 * _S113 / 6.0f;
            }
            else
            {
                _S114 = s_primal_ctx_sin_0(_S113) / _S113;
            }
            float3  _S115 = make_float3 ((_S111 * make_float2 (_S114)).x, (_S111 * make_float2 (_S114)).y, s_primal_ctx_cos_0(_S113));
            raydir_3 = _S115;
        }
        else
        {
            raydir_3 = make_float3 (_S111.x, _S111.y, 1.0f);
        }
        float3  _S116 = normalize_0(s_primal_ctx_mul_0(raydir_3, dpR_0));
        _S112 = true;
        raydir_3 = _S116;
    }
    else
    {
        _S112 = false;
        raydir_3 = _S107;
    }
    *dpray_o_0 = _S108;
    *dpray_d_0 = raydir_3;
    return _S112;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S117, float _S118)
{
    _d_sqrt_0(_S117, _S118);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float _s_dOut_0)
{
    float _S119 = (*dpx_5).primal_0.x;
    float _S120 = (*dpx_5).primal_0.y;
    float _S121 = (*dpx_5).primal_0.z;
    DiffPair_float_0 _S122;
    (&_S122)->primal_0 = _S119 * _S119 + _S120 * _S120 + _S121 * _S121;
    (&_S122)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S122, _s_dOut_0);
    float _S123 = (*dpx_5).primal_0.z * _S122.differential_0;
    float _S124 = _S123 + _S123;
    float _S125 = (*dpx_5).primal_0.y * _S122.differential_0;
    float _S126 = _S125 + _S125;
    float _S127 = (*dpx_5).primal_0.x * _S122.differential_0;
    float _S128 = _S127 + _S127;
    float3  _S129 = make_float3 (0.0f);
    *&((&_S129)->z) = _S124;
    *&((&_S129)->y) = _S126;
    *&((&_S129)->x) = _S128;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S129;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S130, float _S131)
{
    s_bwd_prop_length_impl_0(_S130, _S131);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float3  _s_dOut_1)
{
    float _S132 = length_1((*dpx_6).primal_0);
    float3  _S133 = (*dpx_6).primal_0 * _s_dOut_1;
    float3  _S134 = make_float3 (1.0f / _S132) * _s_dOut_1;
    float _S135 = - ((_S133.x + _S133.y + _S133.z) / (_S132 * _S132));
    float3  _S136 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S137;
    (&_S137)->primal_0 = (*dpx_6).primal_0;
    (&_S137)->differential_0 = _S136;
    s_bwd_length_impl_0(&_S137, _S135);
    float3  _S138 = _S134 + _S137.differential_0;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S138;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S139, float3  _S140)
{
    s_bwd_prop_normalize_impl_0(_S139, _S140);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S141, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S142, float3  _S143)
{
    _d_mul_1(_S141, _S142, _S143);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_7, bool is_fisheye_5, FixedArray<float, 10>  * dist_coeffs_7, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S144 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S145 = *dpt_1;
    float3  _S146 = make_float3 (0.0f);
    bool _S147 = !!_s_diff_ctx_1->_S102;
    float3  raydir_4;
    float3  _S148;
    if(_S147)
    {
        if(is_fisheye_5)
        {
            float _S149 = length_0(_s_diff_ctx_1->_S101);
            float _S150;
            if(_S149 < 0.00100000004749745f)
            {
                _S150 = 1.0f - _S149 * _S149 / 6.0f;
            }
            else
            {
                _S150 = s_primal_ctx_sin_0(_S149) / _S149;
            }
            float3  _S151 = make_float3 ((_s_diff_ctx_1->_S101 * make_float2 (_S150)).x, (_s_diff_ctx_1->_S101 * make_float2 (_S150)).y, s_primal_ctx_cos_0(_S149));
            raydir_4 = _S151;
        }
        else
        {
            raydir_4 = make_float3 (_s_diff_ctx_1->_S101.x, _s_diff_ctx_1->_S101.y, 1.0f);
        }
        float3  _S152 = raydir_4;
        raydir_4 = s_primal_ctx_mul_0(raydir_4, _S144.primal_0);
        _S148 = _S152;
    }
    else
    {
        raydir_4 = _S146;
        _S148 = _S146;
    }
    Matrix<float, 3, 3>  _S153 = makeMatrix<float, 3, 3> (0.0f);
    Matrix<float, 3, 3>  _S154;
    if(_S147)
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S155;
        (&_S155)->primal_0 = raydir_4;
        (&_S155)->differential_0 = _S146;
        s_bwd_normalize_impl_0(&_S155, dpray_d_1);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S156;
        (&_S156)->primal_0 = _S148;
        (&_S156)->differential_0 = _S146;
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S157;
        (&_S157)->primal_0 = _S144.primal_0;
        (&_S157)->differential_0 = _S153;
        s_bwd_prop_mul_0(&_S156, &_S157, _S155.differential_0);
        _S154 = _S157.differential_0;
    }
    else
    {
        _S154 = _S153;
    }
    float3  _S158 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S159;
    (&_S159)->primal_0 = _S145.primal_0;
    (&_S159)->differential_0 = _S146;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S160;
    (&_S160)->primal_0 = _S144.primal_0;
    (&_S160)->differential_0 = _S153;
    s_bwd_prop_mul_0(&_S159, &_S160, _S158);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S159.differential_0;
    Matrix<float, 3, 3>  _S161 = _S160.differential_0 + _S154;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S161;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S162, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S163, float2  _S164, bool _S165, FixedArray<float, 10>  * _S166, float3  _S167, float3  _S168)
{
    float3  _S169;
    float3  _S170;
    s_bwd_prop_generate_ray_Intermediates_0 _S171;
    bool _S172 = s_primal_ctx_generate_ray_0((*_S162).primal_0, (*_S163).primal_0, _S164, _S165, _S166, &_S169, &_S170, &_S171);
    s_bwd_prop_generate_ray_Intermediates_0 _S173 = _S171;
    s_bwd_prop_generate_ray_0(_S162, _S163, _S164, _S165, _S166, _S167, _S168, &_S173);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_8, bool is_fisheye_6, FixedArray<float, 10>  * dist_coeffs_8, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S174 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S174;
    float3  _S175 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S175;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_8, is_fisheye_6, dist_coeffs_8, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S176 = float(width_0);
    float _S177 = float(height_0);
    float _S178 = 0.30000001192092896f * (0.5f * _S176 / fx_0);
    float _S179 = 0.30000001192092896f * (0.5f * _S177 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_2 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S176 - cx_0) / fx_0 + _S178), ((F32_max((- (cx_0 / fx_0 + _S178)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S177 - cy_0) / fy_0 + _S179), ((F32_max((- (cy_0 / fy_0 + _S179)), (mean3d_0.y * rz_0))))))) * rz2_0);
    *cov2d_0 = mul_6(mul_5(J_2, cov3d_0), transpose_1(J_2));
    *mean2d_0 = make_float2 (fx_0 * mean3d_0.x * rz_0 + cx_0, fy_0 * mean3d_0.y * rz_0 + cy_0);
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_7, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_10)
{
    DiffPair_float_0 _S180 = *dpx_7;
    bool _S181;
    if(((*dpx_7).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S181 = ((*dpx_7).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S181 = false;
    }
    float _S182;
    if(_S181)
    {
        _S182 = dOut_10;
    }
    else
    {
        _S182 = 0.0f;
    }
    dpx_7->primal_0 = _S180.primal_0;
    dpx_7->differential_0 = _S182;
    DiffPair_float_0 _S183 = *dpMin_0;
    if((_S180.primal_0) < ((*dpMin_0).primal_0))
    {
        _S182 = dOut_10;
    }
    else
    {
        _S182 = 0.0f;
    }
    dpMin_0->primal_0 = _S183.primal_0;
    dpMin_0->differential_0 = _S182;
    DiffPair_float_0 _S184 = *dpMax_0;
    if(((*dpx_7).primal_0) > ((*dpMax_0).primal_0))
    {
        _S182 = dOut_10;
    }
    else
    {
        _S182 = 0.0f;
    }
    dpMax_0->primal_0 = _S184.primal_0;
    dpMax_0->differential_0 = _S182;
    return;
}

inline __device__ float clamp_0(float x_12, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_12), (minBound_0)))), (maxBound_0)));
}

struct SigmaPoints_0
{
    FixedArray<float3 , 7>  p_0;
    FixedArray<float, 7>  w_mean_0;
    FixedArray<float, 7>  w_cov_0;
};

inline __device__ bool persp_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_9, uint width_1, uint height_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float2  * _S185;
    float2  * _S186;
    float2  * _S187;
    bool _S188;
    float2  * _S189;
    float2  * _S190;
    float2  * _S191;
    bool _S192;
    float2  * _S193;
    bool _S194;
    int2  _S195 = make_int2 (int(0));
    float2  _S196 = make_float2 ((float)_S195.x, (float)_S195.y);
    *mean2d_1 = _S196;
    *cov2d_1 = makeMatrix<float, 2, 2> (0.0f);
    float fx_1 = intrins_0.x;
    float fy_1 = intrins_0.y;
    float _S197 = float(width_1);
    float _S198 = float(height_1);
    float _S199 = 0.30000001192092896f * (0.5f * _S197 / fx_1) * fx_1;
    float lim_x_pos_0 = _S197 + _S199;
    float _S200 = 0.30000001192092896f * (0.5f * _S198 / fy_1) * fy_1;
    float lim_y_pos_0 = _S198 + _S200;
    FixedArray<float2 , 7>  proj_points_0;
    for(;;)
    {
        bool _S201;
        _S185 = &proj_points_0[int(0)];
        for(;;)
        {
            float _S202 = sigmas_0->p_0[int(0)].z;
            proj_points_0[int(0)] = float2 {sigmas_0->p_0[int(0)].x, sigmas_0->p_0[int(0)].y} / make_float2 (_S202);
            if(_S202 < 0.0f)
            {
                _S201 = true;
            }
            else
            {
                bool _S203 = is_valid_distortion(proj_points_0[int(0)], dist_coeffs_9);
                _S201 = !_S203;
            }
            if(_S201)
            {
                break;
            }
            float u_4 = proj_points_0[int(0)].x;
            float v_4 = proj_points_0[int(0)].y;
            float r2_4 = u_4 * u_4 + v_4 * v_4;
            float2  _S204 = proj_points_0[int(0)] * make_float2 (1.0f + r2_4 * ((*dist_coeffs_9)[int(0)] + r2_4 * ((*dist_coeffs_9)[int(1)] + r2_4 * ((*dist_coeffs_9)[int(2)] + r2_4 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_4 * v_4 + (*dist_coeffs_9)[int(5)] * (r2_4 + 2.0f * u_4 * u_4) + (*dist_coeffs_9)[int(6)] * r2_4, 2.0f * (*dist_coeffs_9)[int(5)] * u_4 * v_4 + (*dist_coeffs_9)[int(4)] * (r2_4 + 2.0f * v_4 * v_4) + (*dist_coeffs_9)[int(7)] * r2_4);
            float2  _S205 = _S204 + make_float2 ((*dist_coeffs_9)[int(8)] * _S204.x + (*dist_coeffs_9)[int(9)] * _S204.y, 0.0f);
            proj_points_0[int(0)] = make_float2 (fx_1 * _S205.x + intrins_0.z, fy_1 * _S205.y + intrins_0.w);
            break;
        }
        bool all_valid_0 = true & (!_S201);
        _S186 = &proj_points_0[int(1)];
        for(;;)
        {
            float _S206 = sigmas_0->p_0[int(1)].z;
            proj_points_0[int(1)] = float2 {sigmas_0->p_0[int(1)].x, sigmas_0->p_0[int(1)].y} / make_float2 (_S206);
            if(_S206 < 0.0f)
            {
                _S201 = true;
            }
            else
            {
                bool _S207 = is_valid_distortion(proj_points_0[int(1)], dist_coeffs_9);
                _S201 = !_S207;
            }
            if(_S201)
            {
                break;
            }
            float u_5 = proj_points_0[int(1)].x;
            float v_5 = proj_points_0[int(1)].y;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float2  _S208 = proj_points_0[int(1)] * make_float2 (1.0f + r2_5 * ((*dist_coeffs_9)[int(0)] + r2_5 * ((*dist_coeffs_9)[int(1)] + r2_5 * ((*dist_coeffs_9)[int(2)] + r2_5 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_5 * v_5 + (*dist_coeffs_9)[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + (*dist_coeffs_9)[int(6)] * r2_5, 2.0f * (*dist_coeffs_9)[int(5)] * u_5 * v_5 + (*dist_coeffs_9)[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + (*dist_coeffs_9)[int(7)] * r2_5);
            float2  _S209 = _S208 + make_float2 ((*dist_coeffs_9)[int(8)] * _S208.x + (*dist_coeffs_9)[int(9)] * _S208.y, 0.0f);
            proj_points_0[int(1)] = make_float2 (fx_1 * _S209.x + intrins_0.z, fy_1 * _S209.y + intrins_0.w);
            break;
        }
        bool all_valid_1 = all_valid_0 & (!_S201);
        for(;;)
        {
            _S187 = &proj_points_0[int(2)];
            for(;;)
            {
                float _S210 = sigmas_0->p_0[int(2)].z;
                proj_points_0[int(2)] = float2 {sigmas_0->p_0[int(2)].x, sigmas_0->p_0[int(2)].y} / make_float2 (_S210);
                if(_S210 < 0.0f)
                {
                    _S201 = true;
                }
                else
                {
                    bool _S211 = is_valid_distortion(proj_points_0[int(2)], dist_coeffs_9);
                    _S201 = !_S211;
                }
                if(_S201)
                {
                    break;
                }
                float u_6 = proj_points_0[int(2)].x;
                float v_6 = proj_points_0[int(2)].y;
                float r2_6 = u_6 * u_6 + v_6 * v_6;
                float2  _S212 = proj_points_0[int(2)] * make_float2 (1.0f + r2_6 * ((*dist_coeffs_9)[int(0)] + r2_6 * ((*dist_coeffs_9)[int(1)] + r2_6 * ((*dist_coeffs_9)[int(2)] + r2_6 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_6 * v_6 + (*dist_coeffs_9)[int(5)] * (r2_6 + 2.0f * u_6 * u_6) + (*dist_coeffs_9)[int(6)] * r2_6, 2.0f * (*dist_coeffs_9)[int(5)] * u_6 * v_6 + (*dist_coeffs_9)[int(4)] * (r2_6 + 2.0f * v_6 * v_6) + (*dist_coeffs_9)[int(7)] * r2_6);
                float2  _S213 = _S212 + make_float2 ((*dist_coeffs_9)[int(8)] * _S212.x + (*dist_coeffs_9)[int(9)] * _S212.y, 0.0f);
                proj_points_0[int(2)] = make_float2 (fx_1 * _S213.x + intrins_0.z, fy_1 * _S213.y + intrins_0.w);
                break;
            }
            _S188 = all_valid_1 & (!_S201);
            break;
        }
        _S189 = &proj_points_0[int(3)];
        for(;;)
        {
            float _S214 = sigmas_0->p_0[int(3)].z;
            proj_points_0[int(3)] = float2 {sigmas_0->p_0[int(3)].x, sigmas_0->p_0[int(3)].y} / make_float2 (_S214);
            if(_S214 < 0.0f)
            {
                _S201 = true;
            }
            else
            {
                bool _S215 = is_valid_distortion(proj_points_0[int(3)], dist_coeffs_9);
                _S201 = !_S215;
            }
            if(_S201)
            {
                break;
            }
            float u_7 = proj_points_0[int(3)].x;
            float v_7 = proj_points_0[int(3)].y;
            float r2_7 = u_7 * u_7 + v_7 * v_7;
            float2  _S216 = proj_points_0[int(3)] * make_float2 (1.0f + r2_7 * ((*dist_coeffs_9)[int(0)] + r2_7 * ((*dist_coeffs_9)[int(1)] + r2_7 * ((*dist_coeffs_9)[int(2)] + r2_7 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_7 * v_7 + (*dist_coeffs_9)[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + (*dist_coeffs_9)[int(6)] * r2_7, 2.0f * (*dist_coeffs_9)[int(5)] * u_7 * v_7 + (*dist_coeffs_9)[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + (*dist_coeffs_9)[int(7)] * r2_7);
            float2  _S217 = _S216 + make_float2 ((*dist_coeffs_9)[int(8)] * _S216.x + (*dist_coeffs_9)[int(9)] * _S216.y, 0.0f);
            proj_points_0[int(3)] = make_float2 (fx_1 * _S217.x + intrins_0.z, fy_1 * _S217.y + intrins_0.w);
            break;
        }
        bool all_valid_2 = _S188 & (!_S201);
        _S190 = &proj_points_0[int(4)];
        for(;;)
        {
            float _S218 = sigmas_0->p_0[int(4)].z;
            proj_points_0[int(4)] = float2 {sigmas_0->p_0[int(4)].x, sigmas_0->p_0[int(4)].y} / make_float2 (_S218);
            if(_S218 < 0.0f)
            {
                _S201 = true;
            }
            else
            {
                bool _S219 = is_valid_distortion(proj_points_0[int(4)], dist_coeffs_9);
                _S201 = !_S219;
            }
            if(_S201)
            {
                break;
            }
            float u_8 = proj_points_0[int(4)].x;
            float v_8 = proj_points_0[int(4)].y;
            float r2_8 = u_8 * u_8 + v_8 * v_8;
            float2  _S220 = proj_points_0[int(4)] * make_float2 (1.0f + r2_8 * ((*dist_coeffs_9)[int(0)] + r2_8 * ((*dist_coeffs_9)[int(1)] + r2_8 * ((*dist_coeffs_9)[int(2)] + r2_8 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_8 * v_8 + (*dist_coeffs_9)[int(5)] * (r2_8 + 2.0f * u_8 * u_8) + (*dist_coeffs_9)[int(6)] * r2_8, 2.0f * (*dist_coeffs_9)[int(5)] * u_8 * v_8 + (*dist_coeffs_9)[int(4)] * (r2_8 + 2.0f * v_8 * v_8) + (*dist_coeffs_9)[int(7)] * r2_8);
            float2  _S221 = _S220 + make_float2 ((*dist_coeffs_9)[int(8)] * _S220.x + (*dist_coeffs_9)[int(9)] * _S220.y, 0.0f);
            proj_points_0[int(4)] = make_float2 (fx_1 * _S221.x + intrins_0.z, fy_1 * _S221.y + intrins_0.w);
            break;
        }
        bool all_valid_3 = all_valid_2 & (!_S201);
        for(;;)
        {
            _S191 = &proj_points_0[int(5)];
            for(;;)
            {
                float _S222 = sigmas_0->p_0[int(5)].z;
                proj_points_0[int(5)] = float2 {sigmas_0->p_0[int(5)].x, sigmas_0->p_0[int(5)].y} / make_float2 (_S222);
                if(_S222 < 0.0f)
                {
                    _S201 = true;
                }
                else
                {
                    bool _S223 = is_valid_distortion(proj_points_0[int(5)], dist_coeffs_9);
                    _S201 = !_S223;
                }
                if(_S201)
                {
                    break;
                }
                float u_9 = proj_points_0[int(5)].x;
                float v_9 = proj_points_0[int(5)].y;
                float r2_9 = u_9 * u_9 + v_9 * v_9;
                float2  _S224 = proj_points_0[int(5)] * make_float2 (1.0f + r2_9 * ((*dist_coeffs_9)[int(0)] + r2_9 * ((*dist_coeffs_9)[int(1)] + r2_9 * ((*dist_coeffs_9)[int(2)] + r2_9 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_9 * v_9 + (*dist_coeffs_9)[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + (*dist_coeffs_9)[int(6)] * r2_9, 2.0f * (*dist_coeffs_9)[int(5)] * u_9 * v_9 + (*dist_coeffs_9)[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + (*dist_coeffs_9)[int(7)] * r2_9);
                float2  _S225 = _S224 + make_float2 ((*dist_coeffs_9)[int(8)] * _S224.x + (*dist_coeffs_9)[int(9)] * _S224.y, 0.0f);
                proj_points_0[int(5)] = make_float2 (fx_1 * _S225.x + intrins_0.z, fy_1 * _S225.y + intrins_0.w);
                break;
            }
            _S192 = all_valid_3 & (!_S201);
            break;
        }
        _S193 = &proj_points_0[int(6)];
        for(;;)
        {
            float _S226 = sigmas_0->p_0[int(6)].z;
            proj_points_0[int(6)] = float2 {sigmas_0->p_0[int(6)].x, sigmas_0->p_0[int(6)].y} / make_float2 (_S226);
            if(_S226 < 0.0f)
            {
                _S201 = true;
            }
            else
            {
                bool _S227 = is_valid_distortion(proj_points_0[int(6)], dist_coeffs_9);
                _S201 = !_S227;
            }
            if(_S201)
            {
                break;
            }
            float u_10 = proj_points_0[int(6)].x;
            float v_10 = proj_points_0[int(6)].y;
            float r2_10 = u_10 * u_10 + v_10 * v_10;
            float2  _S228 = proj_points_0[int(6)] * make_float2 (1.0f + r2_10 * ((*dist_coeffs_9)[int(0)] + r2_10 * ((*dist_coeffs_9)[int(1)] + r2_10 * ((*dist_coeffs_9)[int(2)] + r2_10 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_10 * v_10 + (*dist_coeffs_9)[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + (*dist_coeffs_9)[int(6)] * r2_10, 2.0f * (*dist_coeffs_9)[int(5)] * u_10 * v_10 + (*dist_coeffs_9)[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + (*dist_coeffs_9)[int(7)] * r2_10);
            float2  _S229 = _S228 + make_float2 ((*dist_coeffs_9)[int(8)] * _S228.x + (*dist_coeffs_9)[int(9)] * _S228.y, 0.0f);
            proj_points_0[int(6)] = make_float2 (fx_1 * _S229.x + intrins_0.z, fy_1 * _S229.y + intrins_0.w);
            break;
        }
        _S194 = _S192 & (!_S201);
        break;
    }
    if(!_S194)
    {
        return false;
    }
    float2  _S230 = *mean2d_1 + make_float2 (sigmas_0->w_mean_0[int(0)]) * *_S185;
    *mean2d_1 = _S230;
    float2  _S231 = _S230 + make_float2 (sigmas_0->w_mean_0[int(1)]) * *_S186;
    *mean2d_1 = _S231;
    float2  _S232 = _S231 + make_float2 (sigmas_0->w_mean_0[int(2)]) * *_S187;
    *mean2d_1 = _S232;
    float2  _S233 = _S232 + make_float2 (sigmas_0->w_mean_0[int(3)]) * *_S189;
    *mean2d_1 = _S233;
    float2  _S234 = _S233 + make_float2 (sigmas_0->w_mean_0[int(4)]) * *_S190;
    *mean2d_1 = _S234;
    float2  _S235 = _S234 + make_float2 (sigmas_0->w_mean_0[int(5)]) * *_S191;
    *mean2d_1 = _S235;
    float2  _S236 = _S235 + make_float2 (sigmas_0->w_mean_0[int(6)]) * *_S193;
    *mean2d_1 = _S236;
    float _S237 = - _S199;
    float _S238 = - _S200;
    float2  _S239 = make_float2 (clamp_0(_S236.x, _S237, lim_x_pos_0), clamp_0(_S236.y, _S238, lim_y_pos_0));
    float2  d_0 = make_float2 (clamp_0((*_S185).x, _S237, lim_x_pos_0), clamp_0((*_S185).y, _S238, lim_y_pos_0)) - _S239;
    float _S240 = d_0.x;
    float _S241 = d_0.y;
    float _S242 = _S240 * _S241;
    Matrix<float, 2, 2>  _S243 = *cov2d_1 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S240 * _S240, _S242, _S242, _S241 * _S241);
    *cov2d_1 = _S243;
    float2  d_1 = make_float2 (clamp_0((*_S186).x, _S237, lim_x_pos_0), clamp_0((*_S186).y, _S238, lim_y_pos_0)) - _S239;
    float _S244 = d_1.x;
    float _S245 = d_1.y;
    float _S246 = _S244 * _S245;
    Matrix<float, 2, 2>  _S247 = _S243 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S244 * _S244, _S246, _S246, _S245 * _S245);
    *cov2d_1 = _S247;
    float2  d_2 = make_float2 (clamp_0((*_S187).x, _S237, lim_x_pos_0), clamp_0((*_S187).y, _S238, lim_y_pos_0)) - _S239;
    float _S248 = d_2.x;
    float _S249 = d_2.y;
    float _S250 = _S248 * _S249;
    Matrix<float, 2, 2>  _S251 = _S247 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S248 * _S248, _S250, _S250, _S249 * _S249);
    *cov2d_1 = _S251;
    float2  d_3 = make_float2 (clamp_0((*_S189).x, _S237, lim_x_pos_0), clamp_0((*_S189).y, _S238, lim_y_pos_0)) - _S239;
    float _S252 = d_3.x;
    float _S253 = d_3.y;
    float _S254 = _S252 * _S253;
    Matrix<float, 2, 2>  _S255 = _S251 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S252 * _S252, _S254, _S254, _S253 * _S253);
    *cov2d_1 = _S255;
    float2  d_4 = make_float2 (clamp_0((*_S190).x, _S237, lim_x_pos_0), clamp_0((*_S190).y, _S238, lim_y_pos_0)) - _S239;
    float _S256 = d_4.x;
    float _S257 = d_4.y;
    float _S258 = _S256 * _S257;
    Matrix<float, 2, 2>  _S259 = _S255 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S256 * _S256, _S258, _S258, _S257 * _S257);
    *cov2d_1 = _S259;
    float2  d_5 = make_float2 (clamp_0((*_S191).x, _S237, lim_x_pos_0), clamp_0((*_S191).y, _S238, lim_y_pos_0)) - _S239;
    float _S260 = d_5.x;
    float _S261 = d_5.y;
    float _S262 = _S260 * _S261;
    Matrix<float, 2, 2>  _S263 = _S259 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S260 * _S260, _S262, _S262, _S261 * _S261);
    *cov2d_1 = _S263;
    float2  d_6 = make_float2 (clamp_0((*_S193).x, _S237, lim_x_pos_0), clamp_0((*_S193).y, _S238, lim_y_pos_0)) - _S239;
    float _S264 = d_6.x;
    float _S265 = d_6.y;
    float _S266 = _S264 * _S265;
    *cov2d_1 = _S263 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S264 * _S264, _S266, _S266, _S265 * _S265);
    return true;
}

inline __device__ bool persp_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_1, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_10, uint width_2, uint height_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    int2  _S267 = make_int2 (int(0));
    float2  _S268 = make_float2 ((float)_S267.x, (float)_S267.y);
    *mean2d_2 = _S268;
    *cov2d_2 = makeMatrix<float, 2, 2> (0.0f);
    float fx_2 = intrins_1.x;
    float fy_2 = intrins_1.y;
    float _S269 = float(width_2);
    float _S270 = float(height_2);
    float _S271 = 0.30000001192092896f * (0.5f * _S269 / fx_2) * fx_2;
    float lim_x_pos_1 = _S269 + _S271;
    float _S272 = 0.30000001192092896f * (0.5f * _S270 / fy_2) * fy_2;
    float lim_y_pos_1 = _S270 + _S272;
    FixedArray<float2 , 7>  proj_points_1;
    float2  _S273 = float2 {sigmas_1->p_0[int(0)].x, sigmas_1->p_0[int(0)].y} / make_float2 (sigmas_1->p_0[int(0)].z);
    float u_11 = _S273.x;
    float v_11 = _S273.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float _S274 = 2.0f * (*dist_coeffs_10)[int(4)];
    float _S275 = 2.0f * (*dist_coeffs_10)[int(5)];
    float2  _S276 = _S273 * make_float2 (1.0f + r2_11 * ((*dist_coeffs_10)[int(0)] + r2_11 * ((*dist_coeffs_10)[int(1)] + r2_11 * ((*dist_coeffs_10)[int(2)] + r2_11 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S274 * u_11 * v_11 + (*dist_coeffs_10)[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + (*dist_coeffs_10)[int(6)] * r2_11, _S275 * u_11 * v_11 + (*dist_coeffs_10)[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + (*dist_coeffs_10)[int(7)] * r2_11);
    float2  _S277 = _S276 + make_float2 ((*dist_coeffs_10)[int(8)] * _S276.x + (*dist_coeffs_10)[int(9)] * _S276.y, 0.0f);
    float cx_1 = intrins_1.z;
    float cy_1 = intrins_1.w;
    float _S278 = fx_2 * _S277.x + cx_1;
    float _S279 = fy_2 * _S277.y + cy_1;
    float2  _S280 = make_float2 (_S278, _S279);
    proj_points_1[int(0)] = _S280;
    float2  _S281 = float2 {sigmas_1->p_0[int(1)].x, sigmas_1->p_0[int(1)].y} / make_float2 (sigmas_1->p_0[int(1)].z);
    float u_12 = _S281.x;
    float v_12 = _S281.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float2  _S282 = _S281 * make_float2 (1.0f + r2_12 * ((*dist_coeffs_10)[int(0)] + r2_12 * ((*dist_coeffs_10)[int(1)] + r2_12 * ((*dist_coeffs_10)[int(2)] + r2_12 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S274 * u_12 * v_12 + (*dist_coeffs_10)[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + (*dist_coeffs_10)[int(6)] * r2_12, _S275 * u_12 * v_12 + (*dist_coeffs_10)[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + (*dist_coeffs_10)[int(7)] * r2_12);
    float2  _S283 = _S282 + make_float2 ((*dist_coeffs_10)[int(8)] * _S282.x + (*dist_coeffs_10)[int(9)] * _S282.y, 0.0f);
    float _S284 = fx_2 * _S283.x + cx_1;
    float _S285 = fy_2 * _S283.y + cy_1;
    float2  _S286 = make_float2 (_S284, _S285);
    proj_points_1[int(1)] = _S286;
    float2  _S287 = float2 {sigmas_1->p_0[int(2)].x, sigmas_1->p_0[int(2)].y} / make_float2 (sigmas_1->p_0[int(2)].z);
    float u_13 = _S287.x;
    float v_13 = _S287.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S288 = _S287 * make_float2 (1.0f + r2_13 * ((*dist_coeffs_10)[int(0)] + r2_13 * ((*dist_coeffs_10)[int(1)] + r2_13 * ((*dist_coeffs_10)[int(2)] + r2_13 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S274 * u_13 * v_13 + (*dist_coeffs_10)[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + (*dist_coeffs_10)[int(6)] * r2_13, _S275 * u_13 * v_13 + (*dist_coeffs_10)[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + (*dist_coeffs_10)[int(7)] * r2_13);
    float2  _S289 = _S288 + make_float2 ((*dist_coeffs_10)[int(8)] * _S288.x + (*dist_coeffs_10)[int(9)] * _S288.y, 0.0f);
    float _S290 = fx_2 * _S289.x + cx_1;
    float _S291 = fy_2 * _S289.y + cy_1;
    float2  _S292 = make_float2 (_S290, _S291);
    proj_points_1[int(2)] = _S292;
    float2  _S293 = float2 {sigmas_1->p_0[int(3)].x, sigmas_1->p_0[int(3)].y} / make_float2 (sigmas_1->p_0[int(3)].z);
    float u_14 = _S293.x;
    float v_14 = _S293.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S294 = _S293 * make_float2 (1.0f + r2_14 * ((*dist_coeffs_10)[int(0)] + r2_14 * ((*dist_coeffs_10)[int(1)] + r2_14 * ((*dist_coeffs_10)[int(2)] + r2_14 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S274 * u_14 * v_14 + (*dist_coeffs_10)[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + (*dist_coeffs_10)[int(6)] * r2_14, _S275 * u_14 * v_14 + (*dist_coeffs_10)[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + (*dist_coeffs_10)[int(7)] * r2_14);
    float2  _S295 = _S294 + make_float2 ((*dist_coeffs_10)[int(8)] * _S294.x + (*dist_coeffs_10)[int(9)] * _S294.y, 0.0f);
    float _S296 = fx_2 * _S295.x + cx_1;
    float _S297 = fy_2 * _S295.y + cy_1;
    float2  _S298 = make_float2 (_S296, _S297);
    proj_points_1[int(3)] = _S298;
    float2  _S299 = float2 {sigmas_1->p_0[int(4)].x, sigmas_1->p_0[int(4)].y} / make_float2 (sigmas_1->p_0[int(4)].z);
    float u_15 = _S299.x;
    float v_15 = _S299.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S300 = _S299 * make_float2 (1.0f + r2_15 * ((*dist_coeffs_10)[int(0)] + r2_15 * ((*dist_coeffs_10)[int(1)] + r2_15 * ((*dist_coeffs_10)[int(2)] + r2_15 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S274 * u_15 * v_15 + (*dist_coeffs_10)[int(5)] * (r2_15 + 2.0f * u_15 * u_15) + (*dist_coeffs_10)[int(6)] * r2_15, _S275 * u_15 * v_15 + (*dist_coeffs_10)[int(4)] * (r2_15 + 2.0f * v_15 * v_15) + (*dist_coeffs_10)[int(7)] * r2_15);
    float2  _S301 = _S300 + make_float2 ((*dist_coeffs_10)[int(8)] * _S300.x + (*dist_coeffs_10)[int(9)] * _S300.y, 0.0f);
    float _S302 = fx_2 * _S301.x + cx_1;
    float _S303 = fy_2 * _S301.y + cy_1;
    float2  _S304 = make_float2 (_S302, _S303);
    proj_points_1[int(4)] = _S304;
    float2  _S305 = float2 {sigmas_1->p_0[int(5)].x, sigmas_1->p_0[int(5)].y} / make_float2 (sigmas_1->p_0[int(5)].z);
    float u_16 = _S305.x;
    float v_16 = _S305.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float2  _S306 = _S305 * make_float2 (1.0f + r2_16 * ((*dist_coeffs_10)[int(0)] + r2_16 * ((*dist_coeffs_10)[int(1)] + r2_16 * ((*dist_coeffs_10)[int(2)] + r2_16 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S274 * u_16 * v_16 + (*dist_coeffs_10)[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + (*dist_coeffs_10)[int(6)] * r2_16, _S275 * u_16 * v_16 + (*dist_coeffs_10)[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + (*dist_coeffs_10)[int(7)] * r2_16);
    float2  _S307 = _S306 + make_float2 ((*dist_coeffs_10)[int(8)] * _S306.x + (*dist_coeffs_10)[int(9)] * _S306.y, 0.0f);
    float _S308 = fx_2 * _S307.x + cx_1;
    float _S309 = fy_2 * _S307.y + cy_1;
    float2  _S310 = make_float2 (_S308, _S309);
    proj_points_1[int(5)] = _S310;
    float2  _S311 = float2 {sigmas_1->p_0[int(6)].x, sigmas_1->p_0[int(6)].y} / make_float2 (sigmas_1->p_0[int(6)].z);
    float u_17 = _S311.x;
    float v_17 = _S311.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float2  _S312 = _S311 * make_float2 (1.0f + r2_17 * ((*dist_coeffs_10)[int(0)] + r2_17 * ((*dist_coeffs_10)[int(1)] + r2_17 * ((*dist_coeffs_10)[int(2)] + r2_17 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S274 * u_17 * v_17 + (*dist_coeffs_10)[int(5)] * (r2_17 + 2.0f * u_17 * u_17) + (*dist_coeffs_10)[int(6)] * r2_17, _S275 * u_17 * v_17 + (*dist_coeffs_10)[int(4)] * (r2_17 + 2.0f * v_17 * v_17) + (*dist_coeffs_10)[int(7)] * r2_17);
    float2  _S313 = _S312 + make_float2 ((*dist_coeffs_10)[int(8)] * _S312.x + (*dist_coeffs_10)[int(9)] * _S312.y, 0.0f);
    float _S314 = fx_2 * _S313.x + cx_1;
    float _S315 = fy_2 * _S313.y + cy_1;
    float2  _S316 = make_float2 (_S314, _S315);
    proj_points_1[int(6)] = _S316;
    float2  _S317 = *mean2d_2 + make_float2 (sigmas_1->w_mean_0[int(0)]) * _S280;
    *mean2d_2 = _S317;
    float2  _S318 = _S317 + make_float2 (sigmas_1->w_mean_0[int(1)]) * _S286;
    *mean2d_2 = _S318;
    float2  _S319 = _S318 + make_float2 (sigmas_1->w_mean_0[int(2)]) * _S292;
    *mean2d_2 = _S319;
    float2  _S320 = _S319 + make_float2 (sigmas_1->w_mean_0[int(3)]) * _S298;
    *mean2d_2 = _S320;
    float2  _S321 = _S320 + make_float2 (sigmas_1->w_mean_0[int(4)]) * _S304;
    *mean2d_2 = _S321;
    float2  _S322 = _S321 + make_float2 (sigmas_1->w_mean_0[int(5)]) * _S310;
    *mean2d_2 = _S322;
    float2  _S323 = _S322 + make_float2 (sigmas_1->w_mean_0[int(6)]) * _S316;
    *mean2d_2 = _S323;
    float _S324 = - _S271;
    float _S325 = - _S272;
    float2  _S326 = make_float2 (clamp_0(_S323.x, _S324, lim_x_pos_1), clamp_0(_S323.y, _S325, lim_y_pos_1));
    float2  d_7 = make_float2 (clamp_0(_S278, _S324, lim_x_pos_1), clamp_0(_S279, _S325, lim_y_pos_1)) - _S326;
    float _S327 = d_7.x;
    float _S328 = d_7.y;
    float _S329 = _S327 * _S328;
    Matrix<float, 2, 2>  _S330 = *cov2d_2 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S327 * _S327, _S329, _S329, _S328 * _S328);
    *cov2d_2 = _S330;
    float2  d_8 = make_float2 (clamp_0(_S284, _S324, lim_x_pos_1), clamp_0(_S285, _S325, lim_y_pos_1)) - _S326;
    float _S331 = d_8.x;
    float _S332 = d_8.y;
    float _S333 = _S331 * _S332;
    Matrix<float, 2, 2>  _S334 = _S330 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S331 * _S331, _S333, _S333, _S332 * _S332);
    *cov2d_2 = _S334;
    float2  d_9 = make_float2 (clamp_0(_S290, _S324, lim_x_pos_1), clamp_0(_S291, _S325, lim_y_pos_1)) - _S326;
    float _S335 = d_9.x;
    float _S336 = d_9.y;
    float _S337 = _S335 * _S336;
    Matrix<float, 2, 2>  _S338 = _S334 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S335 * _S335, _S337, _S337, _S336 * _S336);
    *cov2d_2 = _S338;
    float2  d_10 = make_float2 (clamp_0(_S296, _S324, lim_x_pos_1), clamp_0(_S297, _S325, lim_y_pos_1)) - _S326;
    float _S339 = d_10.x;
    float _S340 = d_10.y;
    float _S341 = _S339 * _S340;
    Matrix<float, 2, 2>  _S342 = _S338 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S339 * _S339, _S341, _S341, _S340 * _S340);
    *cov2d_2 = _S342;
    float2  d_11 = make_float2 (clamp_0(_S302, _S324, lim_x_pos_1), clamp_0(_S303, _S325, lim_y_pos_1)) - _S326;
    float _S343 = d_11.x;
    float _S344 = d_11.y;
    float _S345 = _S343 * _S344;
    Matrix<float, 2, 2>  _S346 = _S342 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S343 * _S343, _S345, _S345, _S344 * _S344);
    *cov2d_2 = _S346;
    float2  d_12 = make_float2 (clamp_0(_S308, _S324, lim_x_pos_1), clamp_0(_S309, _S325, lim_y_pos_1)) - _S326;
    float _S347 = d_12.x;
    float _S348 = d_12.y;
    float _S349 = _S347 * _S348;
    Matrix<float, 2, 2>  _S350 = _S346 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S347 * _S347, _S349, _S349, _S348 * _S348);
    *cov2d_2 = _S350;
    float2  d_13 = make_float2 (clamp_0(_S314, _S324, lim_x_pos_1), clamp_0(_S315, _S325, lim_y_pos_1)) - _S326;
    float _S351 = d_13.x;
    float _S352 = d_13.y;
    float _S353 = _S351 * _S352;
    *cov2d_2 = _S350 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S351 * _S351, _S353, _S353, _S352 * _S352);
    return true;
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_8, float dOut_11)
{
    DiffPair_float_0 _S354 = *dpx_8;
    float _S355 = - (*dpy_4).primal_0 / ((*dpx_8).primal_0 * (*dpx_8).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_11;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S355;
    float _S356 = _S354.primal_0 / (_S354.primal_0 * _S354.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_11;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S356;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S357, float _S358)
{
    return (F32_atan2((_S357), (_S358)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S359, DiffPair_float_0 * _S360, float _S361)
{
    _d_atan2_0(_S359, _S360, _S361);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_9, float _s_dOut_2)
{
    float _S362 = (*dpx_9).primal_0.x;
    float _S363 = (*dpx_9).primal_0.y;
    DiffPair_float_0 _S364;
    (&_S364)->primal_0 = _S362 * _S362 + _S363 * _S363;
    (&_S364)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S364, _s_dOut_2);
    float _S365 = (*dpx_9).primal_0.y * _S364.differential_0;
    float _S366 = _S365 + _S365;
    float _S367 = (*dpx_9).primal_0.x * _S364.differential_0;
    float _S368 = _S367 + _S367;
    float2  _S369 = make_float2 (0.0f);
    *&((&_S369)->y) = _S366;
    *&((&_S369)->x) = _S368;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S369;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S370, float _S371)
{
    s_bwd_prop_length_impl_1(_S370, _S371);
    return;
}

inline __device__ bool fisheye_proj_3dgs_0(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_11, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    float k_0;
    float2  _S372;
    float _S373;
    float _S374;
    bool _S375;
    for(;;)
    {
        float2  _S376 = float2 {mean3d_1.x, mean3d_1.y};
        _S372 = _S376;
        float r_8 = length_0(_S376);
        _S373 = r_8;
        float _S377 = mean3d_1.z;
        _S374 = _S377;
        float theta_4 = (F32_atan2((r_8), (_S377)));
        if(theta_4 < 0.00100000004749745f)
        {
            k_0 = (1.0f - theta_4 * theta_4 / 3.0f) / _S377;
        }
        else
        {
            k_0 = theta_4 / r_8;
        }
        float2  _S378 = _S376 * make_float2 (k_0);
        *mean2d_3 = _S378;
        bool _S379 = is_valid_distortion(_S378, dist_coeffs_11);
        bool _S380 = !_S379;
        _S375 = _S380;
        if(_S380)
        {
            break;
        }
        float u_18 = (*mean2d_3).x;
        float v_18 = (*mean2d_3).y;
        float r2_18 = u_18 * u_18 + v_18 * v_18;
        float2  _S381 = *mean2d_3 * make_float2 (1.0f + r2_18 * ((*dist_coeffs_11)[int(0)] + r2_18 * ((*dist_coeffs_11)[int(1)] + r2_18 * ((*dist_coeffs_11)[int(2)] + r2_18 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_18 * v_18 + (*dist_coeffs_11)[int(5)] * (r2_18 + 2.0f * u_18 * u_18) + (*dist_coeffs_11)[int(6)] * r2_18, 2.0f * (*dist_coeffs_11)[int(5)] * u_18 * v_18 + (*dist_coeffs_11)[int(4)] * (r2_18 + 2.0f * v_18 * v_18) + (*dist_coeffs_11)[int(7)] * r2_18);
        float2  _S382 = _S381 + make_float2 ((*dist_coeffs_11)[int(8)] * _S381.x + (*dist_coeffs_11)[int(9)] * _S381.y, 0.0f);
        *mean2d_3 = make_float2 (intrins_2.x * _S382.x + intrins_2.z, intrins_2.y * _S382.y + intrins_2.w);
        break;
    }
    if(!!_S375)
    {
        return false;
    }
    Matrix<float, 2, 3>  J_3;
    float2  _S383 = make_float2 (0.0f);
    float2  seed_4 = _S383;
    *&((&seed_4)->x) = 1.0f;
    float2  _S384 = seed_4;
    float _S385 = s_primal_ctx_atan2_0(_S373, _S374);
    bool _S386 = _S385 < 0.00100000004749745f;
    float _S387;
    float _S388;
    float _S389;
    if(_S386)
    {
        float _S390 = 1.0f - _S385 * _S385 / 3.0f;
        float _S391 = _S374 * _S374;
        k_0 = _S390 / _S374;
        _S387 = 0.0f;
        _S388 = _S391;
        _S389 = _S390;
    }
    else
    {
        float _S392 = _S373 * _S373;
        k_0 = _S385 / _S373;
        _S387 = _S392;
        _S388 = 0.0f;
        _S389 = 0.0f;
    }
    float2  _S393 = make_float2 (k_0);
    float2  _S394 = _S372 * make_float2 (k_0);
    float u_19 = _S394.x;
    float v_19 = _S394.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S395 = (*dist_coeffs_11)[int(2)] + r2_19 * (*dist_coeffs_11)[int(3)];
    float _S396 = (*dist_coeffs_11)[int(1)] + r2_19 * _S395;
    float _S397 = (*dist_coeffs_11)[int(0)] + r2_19 * _S396;
    float _S398 = 2.0f * (*dist_coeffs_11)[int(4)];
    float _S399 = 2.0f * (*dist_coeffs_11)[int(5)];
    float fx_3 = intrins_2.x;
    float fy_3 = intrins_2.y;
    float _S400 = fx_3 * _S384.x;
    float2  _S401 = make_float2 (_S400, fy_3 * _S384.y) + make_float2 ((*dist_coeffs_11)[int(8)] * _S400, (*dist_coeffs_11)[int(9)] * _S400);
    float2  _S402 = _S394 * _S401;
    float _S403 = (*dist_coeffs_11)[int(4)] * _S401.y;
    float _S404 = (*dist_coeffs_11)[int(5)] * _S401.x;
    float _S405 = _S402.x + _S402.y;
    float _S406 = r2_19 * _S405;
    float _S407 = r2_19 * _S406;
    float _S408 = (*dist_coeffs_11)[int(7)] * _S401.y + _S403 + (*dist_coeffs_11)[int(6)] * _S401.x + _S404 + _S397 * _S405 + _S396 * _S406 + _S395 * _S407 + (*dist_coeffs_11)[int(3)] * (r2_19 * _S407);
    float _S409 = v_19 * _S408;
    float _S410 = u_19 * _S408;
    float2  _S411 = make_float2 (1.0f + r2_19 * _S397) * _S401 + make_float2 (_S399 * (v_19 * _S401.y) + 2.0f * u_19 * _S404 + 2.0f * (u_19 * _S404) + _S398 * (v_19 * _S401.x) + _S410 + _S410, 2.0f * v_19 * _S403 + 2.0f * (v_19 * _S403) + _S399 * u_19 * _S401.y + _S398 * u_19 * _S401.x + _S409 + _S409);
    float2  _S412 = _S372 * _S411;
    float2  _S413 = _S393 * _S411;
    float _S414 = _S412.x + _S412.y;
    if(_S386)
    {
        float _S415 = _S414 / _S388;
        float _S416 = _S389 * - _S415;
        float _S417 = _S385 * (0.3333333432674408f * - (_S374 * _S415));
        k_0 = _S417 + _S417;
        _S387 = _S416;
        _S388 = 0.0f;
    }
    else
    {
        float _S418 = _S414 / _S387;
        float _S419 = _S385 * - _S418;
        k_0 = _S373 * _S418;
        _S387 = 0.0f;
        _S388 = _S419;
    }
    float _S420 = _S373;
    DiffPair_float_0 _S421;
    (&_S421)->primal_0 = _S373;
    (&_S421)->differential_0 = 0.0f;
    float _S422 = _S374;
    DiffPair_float_0 _S423;
    (&_S423)->primal_0 = _S374;
    (&_S423)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S421, &_S423, k_0);
    float _S424 = _S423.differential_0 + _S387;
    float _S425 = _S421.differential_0 + _S388;
    float2  _S426 = make_float2 (0.0f);
    float2  _S427 = _S372;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S428;
    (&_S428)->primal_0 = _S372;
    (&_S428)->differential_0 = _S426;
    s_bwd_length_impl_1(&_S428, _S425);
    float2  _S429 = _S428.differential_0 + _S413;
    float3  _S430 = make_float3 (_S429.x, _S429.y, _S424);
    J_3[int(0)] = _S430;
    float2  seed_5 = _S383;
    *&((&seed_5)->y) = 1.0f;
    float2  _S431 = seed_5;
    if(_S386)
    {
        float _S432 = 1.0f - _S385 * _S385 / 3.0f;
        float _S433 = _S374 * _S374;
        k_0 = _S432 / _S374;
        _S387 = 0.0f;
        _S388 = _S433;
        _S389 = _S432;
    }
    else
    {
        float _S434 = _S373 * _S373;
        k_0 = _S385 / _S373;
        _S387 = _S434;
        _S388 = 0.0f;
        _S389 = 0.0f;
    }
    float2  _S435 = make_float2 (k_0);
    float2  _S436 = _S372 * make_float2 (k_0);
    float u_20 = _S436.x;
    float v_20 = _S436.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S437 = (*dist_coeffs_11)[int(2)] + r2_20 * (*dist_coeffs_11)[int(3)];
    float _S438 = (*dist_coeffs_11)[int(1)] + r2_20 * _S437;
    float _S439 = (*dist_coeffs_11)[int(0)] + r2_20 * _S438;
    float _S440 = fx_3 * _S431.x;
    float2  _S441 = make_float2 (_S440, fy_3 * _S431.y) + make_float2 ((*dist_coeffs_11)[int(8)] * _S440, (*dist_coeffs_11)[int(9)] * _S440);
    float2  _S442 = _S436 * _S441;
    float _S443 = (*dist_coeffs_11)[int(4)] * _S441.y;
    float _S444 = (*dist_coeffs_11)[int(5)] * _S441.x;
    float _S445 = _S442.x + _S442.y;
    float _S446 = r2_20 * _S445;
    float _S447 = r2_20 * _S446;
    float _S448 = (*dist_coeffs_11)[int(7)] * _S441.y + _S443 + (*dist_coeffs_11)[int(6)] * _S441.x + _S444 + _S439 * _S445 + _S438 * _S446 + _S437 * _S447 + (*dist_coeffs_11)[int(3)] * (r2_20 * _S447);
    float _S449 = v_20 * _S448;
    float _S450 = u_20 * _S448;
    float2  _S451 = make_float2 (1.0f + r2_20 * _S439) * _S441 + make_float2 (_S399 * (v_20 * _S441.y) + 2.0f * u_20 * _S444 + 2.0f * (u_20 * _S444) + _S398 * (v_20 * _S441.x) + _S450 + _S450, 2.0f * v_20 * _S443 + 2.0f * (v_20 * _S443) + _S399 * u_20 * _S441.y + _S398 * u_20 * _S441.x + _S449 + _S449);
    float2  _S452 = _S372 * _S451;
    float2  _S453 = _S435 * _S451;
    float _S454 = _S452.x + _S452.y;
    if(_S386)
    {
        float _S455 = _S454 / _S388;
        float _S456 = _S389 * - _S455;
        float _S457 = _S385 * (0.3333333432674408f * - (_S374 * _S455));
        k_0 = _S457 + _S457;
        _S387 = _S456;
        _S388 = 0.0f;
    }
    else
    {
        float _S458 = _S454 / _S387;
        float _S459 = _S385 * - _S458;
        k_0 = _S373 * _S458;
        _S387 = 0.0f;
        _S388 = _S459;
    }
    DiffPair_float_0 _S460;
    (&_S460)->primal_0 = _S420;
    (&_S460)->differential_0 = 0.0f;
    DiffPair_float_0 _S461;
    (&_S461)->primal_0 = _S422;
    (&_S461)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S460, &_S461, k_0);
    float _S462 = _S461.differential_0 + _S387;
    float _S463 = _S460.differential_0 + _S388;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S464;
    (&_S464)->primal_0 = _S427;
    (&_S464)->differential_0 = _S426;
    s_bwd_length_impl_1(&_S464, _S463);
    float2  _S465 = _S464.differential_0 + _S453;
    float3  _S466 = make_float3 (_S465.x, _S465.y, _S462);
    J_3[int(1)] = _S466;
    *cov2d_3 = mul_6(mul_5(J_3, cov3d_1), transpose_1(J_3));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_1(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_12, Matrix<float, 2, 2>  * cov2d_4, float2  * mean2d_4)
{
    float2  _S467 = float2 {mean3d_2.x, mean3d_2.y};
    float r_9 = length_0(_S467);
    float _S468 = mean3d_2.z;
    float theta_5 = (F32_atan2((r_9), (_S468)));
    float k_1;
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S468;
    }
    else
    {
        k_1 = theta_5 / r_9;
    }
    float2  _S469 = _S467 * make_float2 (k_1);
    float u_21 = _S469.x;
    float v_21 = _S469.y;
    float r2_21 = u_21 * u_21 + v_21 * v_21;
    float _S470 = 2.0f * (*dist_coeffs_12)[int(4)];
    float _S471 = 2.0f * (*dist_coeffs_12)[int(5)];
    float2  _S472 = _S469 * make_float2 (1.0f + r2_21 * ((*dist_coeffs_12)[int(0)] + r2_21 * ((*dist_coeffs_12)[int(1)] + r2_21 * ((*dist_coeffs_12)[int(2)] + r2_21 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S470 * u_21 * v_21 + (*dist_coeffs_12)[int(5)] * (r2_21 + 2.0f * u_21 * u_21) + (*dist_coeffs_12)[int(6)] * r2_21, _S471 * u_21 * v_21 + (*dist_coeffs_12)[int(4)] * (r2_21 + 2.0f * v_21 * v_21) + (*dist_coeffs_12)[int(7)] * r2_21);
    float2  _S473 = _S472 + make_float2 ((*dist_coeffs_12)[int(8)] * _S472.x + (*dist_coeffs_12)[int(9)] * _S472.y, 0.0f);
    float fx_4 = intrins_3.x;
    float fy_4 = intrins_3.y;
    *mean2d_4 = make_float2 (fx_4 * _S473.x + intrins_3.z, fy_4 * _S473.y + intrins_3.w);
    Matrix<float, 2, 3>  J_4;
    float2  _S474 = make_float2 (0.0f);
    float2  seed_6 = _S474;
    *&((&seed_6)->x) = 1.0f;
    float2  _S475 = seed_6;
    float _S476 = s_primal_ctx_atan2_0(r_9, _S468);
    bool _S477 = _S476 < 0.00100000004749745f;
    float _S478;
    float _S479;
    float _S480;
    if(_S477)
    {
        float _S481 = 1.0f - _S476 * _S476 / 3.0f;
        float _S482 = _S468 * _S468;
        k_1 = _S481 / _S468;
        _S478 = 0.0f;
        _S479 = _S482;
        _S480 = _S481;
    }
    else
    {
        float _S483 = r_9 * r_9;
        k_1 = _S476 / r_9;
        _S478 = _S483;
        _S479 = 0.0f;
        _S480 = 0.0f;
    }
    float2  _S484 = make_float2 (k_1);
    float2  _S485 = _S467 * make_float2 (k_1);
    float u_22 = _S485.x;
    float v_22 = _S485.y;
    float r2_22 = u_22 * u_22 + v_22 * v_22;
    float _S486 = (*dist_coeffs_12)[int(2)] + r2_22 * (*dist_coeffs_12)[int(3)];
    float _S487 = (*dist_coeffs_12)[int(1)] + r2_22 * _S486;
    float _S488 = (*dist_coeffs_12)[int(0)] + r2_22 * _S487;
    float _S489 = fx_4 * _S475.x;
    float2  _S490 = make_float2 (_S489, fy_4 * _S475.y) + make_float2 ((*dist_coeffs_12)[int(8)] * _S489, (*dist_coeffs_12)[int(9)] * _S489);
    float2  _S491 = _S485 * _S490;
    float _S492 = (*dist_coeffs_12)[int(4)] * _S490.y;
    float _S493 = (*dist_coeffs_12)[int(5)] * _S490.x;
    float _S494 = _S491.x + _S491.y;
    float _S495 = r2_22 * _S494;
    float _S496 = r2_22 * _S495;
    float _S497 = (*dist_coeffs_12)[int(7)] * _S490.y + _S492 + (*dist_coeffs_12)[int(6)] * _S490.x + _S493 + _S488 * _S494 + _S487 * _S495 + _S486 * _S496 + (*dist_coeffs_12)[int(3)] * (r2_22 * _S496);
    float _S498 = v_22 * _S497;
    float _S499 = u_22 * _S497;
    float2  _S500 = make_float2 (1.0f + r2_22 * _S488) * _S490 + make_float2 (_S471 * (v_22 * _S490.y) + 2.0f * u_22 * _S493 + 2.0f * (u_22 * _S493) + _S470 * (v_22 * _S490.x) + _S499 + _S499, 2.0f * v_22 * _S492 + 2.0f * (v_22 * _S492) + _S471 * u_22 * _S490.y + _S470 * u_22 * _S490.x + _S498 + _S498);
    float2  _S501 = _S467 * _S500;
    float2  _S502 = _S484 * _S500;
    float _S503 = _S501.x + _S501.y;
    if(_S477)
    {
        float _S504 = _S503 / _S479;
        float _S505 = _S480 * - _S504;
        float _S506 = _S476 * (0.3333333432674408f * - (_S468 * _S504));
        k_1 = _S506 + _S506;
        _S478 = _S505;
        _S479 = 0.0f;
    }
    else
    {
        float _S507 = _S503 / _S478;
        float _S508 = _S476 * - _S507;
        k_1 = r_9 * _S507;
        _S478 = 0.0f;
        _S479 = _S508;
    }
    DiffPair_float_0 _S509;
    (&_S509)->primal_0 = r_9;
    (&_S509)->differential_0 = 0.0f;
    DiffPair_float_0 _S510;
    (&_S510)->primal_0 = _S468;
    (&_S510)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S509, &_S510, k_1);
    float _S511 = _S510.differential_0 + _S478;
    float _S512 = _S509.differential_0 + _S479;
    float2  _S513 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S514;
    (&_S514)->primal_0 = _S467;
    (&_S514)->differential_0 = _S513;
    s_bwd_length_impl_1(&_S514, _S512);
    float2  _S515 = _S514.differential_0 + _S502;
    float3  _S516 = make_float3 (_S515.x, _S515.y, _S511);
    J_4[int(0)] = _S516;
    float2  seed_7 = _S474;
    *&((&seed_7)->y) = 1.0f;
    float2  _S517 = seed_7;
    if(_S477)
    {
        float _S518 = 1.0f - _S476 * _S476 / 3.0f;
        float _S519 = _S468 * _S468;
        k_1 = _S518 / _S468;
        _S478 = 0.0f;
        _S479 = _S519;
        _S480 = _S518;
    }
    else
    {
        float _S520 = r_9 * r_9;
        k_1 = _S476 / r_9;
        _S478 = _S520;
        _S479 = 0.0f;
        _S480 = 0.0f;
    }
    float2  _S521 = make_float2 (k_1);
    float2  _S522 = _S467 * make_float2 (k_1);
    float u_23 = _S522.x;
    float v_23 = _S522.y;
    float r2_23 = u_23 * u_23 + v_23 * v_23;
    float _S523 = (*dist_coeffs_12)[int(2)] + r2_23 * (*dist_coeffs_12)[int(3)];
    float _S524 = (*dist_coeffs_12)[int(1)] + r2_23 * _S523;
    float _S525 = (*dist_coeffs_12)[int(0)] + r2_23 * _S524;
    float _S526 = fx_4 * _S517.x;
    float2  _S527 = make_float2 (_S526, fy_4 * _S517.y) + make_float2 ((*dist_coeffs_12)[int(8)] * _S526, (*dist_coeffs_12)[int(9)] * _S526);
    float2  _S528 = _S522 * _S527;
    float _S529 = (*dist_coeffs_12)[int(4)] * _S527.y;
    float _S530 = (*dist_coeffs_12)[int(5)] * _S527.x;
    float _S531 = _S528.x + _S528.y;
    float _S532 = r2_23 * _S531;
    float _S533 = r2_23 * _S532;
    float _S534 = (*dist_coeffs_12)[int(7)] * _S527.y + _S529 + (*dist_coeffs_12)[int(6)] * _S527.x + _S530 + _S525 * _S531 + _S524 * _S532 + _S523 * _S533 + (*dist_coeffs_12)[int(3)] * (r2_23 * _S533);
    float _S535 = v_23 * _S534;
    float _S536 = u_23 * _S534;
    float2  _S537 = make_float2 (1.0f + r2_23 * _S525) * _S527 + make_float2 (_S471 * (v_23 * _S527.y) + 2.0f * u_23 * _S530 + 2.0f * (u_23 * _S530) + _S470 * (v_23 * _S527.x) + _S536 + _S536, 2.0f * v_23 * _S529 + 2.0f * (v_23 * _S529) + _S471 * u_23 * _S527.y + _S470 * u_23 * _S527.x + _S535 + _S535);
    float2  _S538 = _S467 * _S537;
    float2  _S539 = _S521 * _S537;
    float _S540 = _S538.x + _S538.y;
    if(_S477)
    {
        float _S541 = _S540 / _S479;
        float _S542 = _S480 * - _S541;
        float _S543 = _S476 * (0.3333333432674408f * - (_S468 * _S541));
        k_1 = _S543 + _S543;
        _S478 = _S542;
        _S479 = 0.0f;
    }
    else
    {
        float _S544 = _S540 / _S478;
        float _S545 = _S476 * - _S544;
        k_1 = r_9 * _S544;
        _S478 = 0.0f;
        _S479 = _S545;
    }
    DiffPair_float_0 _S546;
    (&_S546)->primal_0 = r_9;
    (&_S546)->differential_0 = 0.0f;
    DiffPair_float_0 _S547;
    (&_S547)->primal_0 = _S468;
    (&_S547)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S546, &_S547, k_1);
    float _S548 = _S547.differential_0 + _S478;
    float _S549 = _S546.differential_0 + _S479;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S550;
    (&_S550)->primal_0 = _S467;
    (&_S550)->differential_0 = _S513;
    s_bwd_length_impl_1(&_S550, _S549);
    float2  _S551 = _S550.differential_0 + _S539;
    float3  _S552 = make_float3 (_S551.x, _S551.y, _S548);
    J_4[int(1)] = _S552;
    *cov2d_4 = mul_6(mul_5(J_4, cov3d_2), transpose_1(J_4));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_2, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_13, Matrix<float, 2, 2>  * cov2d_5, float2  * mean2d_5)
{
    float2  * _S553;
    bool _S554;
    float2  * _S555;
    bool _S556;
    float2  * _S557;
    bool _S558;
    bool _S559;
    float2  * _S560;
    bool _S561;
    float2  * _S562;
    bool _S563;
    float2  * _S564;
    bool _S565;
    bool _S566;
    float2  * _S567;
    bool _S568;
    bool _S569;
    int2  _S570 = make_int2 (int(0));
    float2  _S571 = make_float2 ((float)_S570.x, (float)_S570.y);
    *mean2d_5 = _S571;
    *cov2d_5 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_2;
    for(;;)
    {
        float k_2;
        _S553 = &proj_points_2[int(0)];
        for(;;)
        {
            float2  _S572 = float2 {sigmas_2->p_0[int(0)].x, sigmas_2->p_0[int(0)].y};
            float r_10 = length_0(_S572);
            float _S573 = sigmas_2->p_0[int(0)].z;
            float theta_6 = (F32_atan2((r_10), (_S573)));
            if(theta_6 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_6 * theta_6 / 3.0f) / _S573;
            }
            else
            {
                k_2 = theta_6 / r_10;
            }
            float2  _S574 = _S572 * make_float2 (k_2);
            proj_points_2[int(0)] = _S574;
            bool _S575 = is_valid_distortion(_S574, dist_coeffs_13);
            bool _S576 = !_S575;
            _S554 = _S576;
            if(_S576)
            {
                break;
            }
            float u_24 = proj_points_2[int(0)].x;
            float v_24 = proj_points_2[int(0)].y;
            float r2_24 = u_24 * u_24 + v_24 * v_24;
            float2  _S577 = proj_points_2[int(0)] * make_float2 (1.0f + r2_24 * ((*dist_coeffs_13)[int(0)] + r2_24 * ((*dist_coeffs_13)[int(1)] + r2_24 * ((*dist_coeffs_13)[int(2)] + r2_24 * (*dist_coeffs_13)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_13)[int(4)] * u_24 * v_24 + (*dist_coeffs_13)[int(5)] * (r2_24 + 2.0f * u_24 * u_24) + (*dist_coeffs_13)[int(6)] * r2_24, 2.0f * (*dist_coeffs_13)[int(5)] * u_24 * v_24 + (*dist_coeffs_13)[int(4)] * (r2_24 + 2.0f * v_24 * v_24) + (*dist_coeffs_13)[int(7)] * r2_24);
            float2  _S578 = _S577 + make_float2 ((*dist_coeffs_13)[int(8)] * _S577.x + (*dist_coeffs_13)[int(9)] * _S577.y, 0.0f);
            proj_points_2[int(0)] = make_float2 (intrins_4.x * _S578.x + intrins_4.z, intrins_4.y * _S578.y + intrins_4.w);
            break;
        }
        bool all_valid_4 = true & (!_S554);
        _S555 = &proj_points_2[int(1)];
        for(;;)
        {
            float2  _S579 = float2 {sigmas_2->p_0[int(1)].x, sigmas_2->p_0[int(1)].y};
            float r_11 = length_0(_S579);
            float _S580 = sigmas_2->p_0[int(1)].z;
            float theta_7 = (F32_atan2((r_11), (_S580)));
            if(theta_7 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_7 * theta_7 / 3.0f) / _S580;
            }
            else
            {
                k_2 = theta_7 / r_11;
            }
            float2  _S581 = _S579 * make_float2 (k_2);
            proj_points_2[int(1)] = _S581;
            bool _S582 = is_valid_distortion(_S581, dist_coeffs_13);
            bool _S583 = !_S582;
            _S556 = _S583;
            if(_S583)
            {
                break;
            }
            float u_25 = proj_points_2[int(1)].x;
            float v_25 = proj_points_2[int(1)].y;
            float r2_25 = u_25 * u_25 + v_25 * v_25;
            float2  _S584 = proj_points_2[int(1)] * make_float2 (1.0f + r2_25 * ((*dist_coeffs_13)[int(0)] + r2_25 * ((*dist_coeffs_13)[int(1)] + r2_25 * ((*dist_coeffs_13)[int(2)] + r2_25 * (*dist_coeffs_13)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_13)[int(4)] * u_25 * v_25 + (*dist_coeffs_13)[int(5)] * (r2_25 + 2.0f * u_25 * u_25) + (*dist_coeffs_13)[int(6)] * r2_25, 2.0f * (*dist_coeffs_13)[int(5)] * u_25 * v_25 + (*dist_coeffs_13)[int(4)] * (r2_25 + 2.0f * v_25 * v_25) + (*dist_coeffs_13)[int(7)] * r2_25);
            float2  _S585 = _S584 + make_float2 ((*dist_coeffs_13)[int(8)] * _S584.x + (*dist_coeffs_13)[int(9)] * _S584.y, 0.0f);
            proj_points_2[int(1)] = make_float2 (intrins_4.x * _S585.x + intrins_4.z, intrins_4.y * _S585.y + intrins_4.w);
            break;
        }
        bool all_valid_5 = all_valid_4 & (!_S556);
        for(;;)
        {
            _S557 = &proj_points_2[int(2)];
            for(;;)
            {
                float2  _S586 = float2 {sigmas_2->p_0[int(2)].x, sigmas_2->p_0[int(2)].y};
                float r_12 = length_0(_S586);
                float _S587 = sigmas_2->p_0[int(2)].z;
                float theta_8 = (F32_atan2((r_12), (_S587)));
                if(theta_8 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_8 * theta_8 / 3.0f) / _S587;
                }
                else
                {
                    k_2 = theta_8 / r_12;
                }
                float2  _S588 = _S586 * make_float2 (k_2);
                proj_points_2[int(2)] = _S588;
                bool _S589 = is_valid_distortion(_S588, dist_coeffs_13);
                bool _S590 = !_S589;
                _S558 = _S590;
                if(_S590)
                {
                    break;
                }
                float u_26 = proj_points_2[int(2)].x;
                float v_26 = proj_points_2[int(2)].y;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float2  _S591 = proj_points_2[int(2)] * make_float2 (1.0f + r2_26 * ((*dist_coeffs_13)[int(0)] + r2_26 * ((*dist_coeffs_13)[int(1)] + r2_26 * ((*dist_coeffs_13)[int(2)] + r2_26 * (*dist_coeffs_13)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_13)[int(4)] * u_26 * v_26 + (*dist_coeffs_13)[int(5)] * (r2_26 + 2.0f * u_26 * u_26) + (*dist_coeffs_13)[int(6)] * r2_26, 2.0f * (*dist_coeffs_13)[int(5)] * u_26 * v_26 + (*dist_coeffs_13)[int(4)] * (r2_26 + 2.0f * v_26 * v_26) + (*dist_coeffs_13)[int(7)] * r2_26);
                float2  _S592 = _S591 + make_float2 ((*dist_coeffs_13)[int(8)] * _S591.x + (*dist_coeffs_13)[int(9)] * _S591.y, 0.0f);
                proj_points_2[int(2)] = make_float2 (intrins_4.x * _S592.x + intrins_4.z, intrins_4.y * _S592.y + intrins_4.w);
                break;
            }
            _S559 = all_valid_5 & (!_S558);
            break;
        }
        _S560 = &proj_points_2[int(3)];
        for(;;)
        {
            float2  _S593 = float2 {sigmas_2->p_0[int(3)].x, sigmas_2->p_0[int(3)].y};
            float r_13 = length_0(_S593);
            float _S594 = sigmas_2->p_0[int(3)].z;
            float theta_9 = (F32_atan2((r_13), (_S594)));
            if(theta_9 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_9 * theta_9 / 3.0f) / _S594;
            }
            else
            {
                k_2 = theta_9 / r_13;
            }
            float2  _S595 = _S593 * make_float2 (k_2);
            proj_points_2[int(3)] = _S595;
            bool _S596 = is_valid_distortion(_S595, dist_coeffs_13);
            bool _S597 = !_S596;
            _S561 = _S597;
            if(_S597)
            {
                break;
            }
            float u_27 = proj_points_2[int(3)].x;
            float v_27 = proj_points_2[int(3)].y;
            float r2_27 = u_27 * u_27 + v_27 * v_27;
            float2  _S598 = proj_points_2[int(3)] * make_float2 (1.0f + r2_27 * ((*dist_coeffs_13)[int(0)] + r2_27 * ((*dist_coeffs_13)[int(1)] + r2_27 * ((*dist_coeffs_13)[int(2)] + r2_27 * (*dist_coeffs_13)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_13)[int(4)] * u_27 * v_27 + (*dist_coeffs_13)[int(5)] * (r2_27 + 2.0f * u_27 * u_27) + (*dist_coeffs_13)[int(6)] * r2_27, 2.0f * (*dist_coeffs_13)[int(5)] * u_27 * v_27 + (*dist_coeffs_13)[int(4)] * (r2_27 + 2.0f * v_27 * v_27) + (*dist_coeffs_13)[int(7)] * r2_27);
            float2  _S599 = _S598 + make_float2 ((*dist_coeffs_13)[int(8)] * _S598.x + (*dist_coeffs_13)[int(9)] * _S598.y, 0.0f);
            proj_points_2[int(3)] = make_float2 (intrins_4.x * _S599.x + intrins_4.z, intrins_4.y * _S599.y + intrins_4.w);
            break;
        }
        bool all_valid_6 = _S559 & (!_S561);
        _S562 = &proj_points_2[int(4)];
        for(;;)
        {
            float2  _S600 = float2 {sigmas_2->p_0[int(4)].x, sigmas_2->p_0[int(4)].y};
            float r_14 = length_0(_S600);
            float _S601 = sigmas_2->p_0[int(4)].z;
            float theta_10 = (F32_atan2((r_14), (_S601)));
            if(theta_10 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_10 * theta_10 / 3.0f) / _S601;
            }
            else
            {
                k_2 = theta_10 / r_14;
            }
            float2  _S602 = _S600 * make_float2 (k_2);
            proj_points_2[int(4)] = _S602;
            bool _S603 = is_valid_distortion(_S602, dist_coeffs_13);
            bool _S604 = !_S603;
            _S563 = _S604;
            if(_S604)
            {
                break;
            }
            float u_28 = proj_points_2[int(4)].x;
            float v_28 = proj_points_2[int(4)].y;
            float r2_28 = u_28 * u_28 + v_28 * v_28;
            float2  _S605 = proj_points_2[int(4)] * make_float2 (1.0f + r2_28 * ((*dist_coeffs_13)[int(0)] + r2_28 * ((*dist_coeffs_13)[int(1)] + r2_28 * ((*dist_coeffs_13)[int(2)] + r2_28 * (*dist_coeffs_13)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_13)[int(4)] * u_28 * v_28 + (*dist_coeffs_13)[int(5)] * (r2_28 + 2.0f * u_28 * u_28) + (*dist_coeffs_13)[int(6)] * r2_28, 2.0f * (*dist_coeffs_13)[int(5)] * u_28 * v_28 + (*dist_coeffs_13)[int(4)] * (r2_28 + 2.0f * v_28 * v_28) + (*dist_coeffs_13)[int(7)] * r2_28);
            float2  _S606 = _S605 + make_float2 ((*dist_coeffs_13)[int(8)] * _S605.x + (*dist_coeffs_13)[int(9)] * _S605.y, 0.0f);
            proj_points_2[int(4)] = make_float2 (intrins_4.x * _S606.x + intrins_4.z, intrins_4.y * _S606.y + intrins_4.w);
            break;
        }
        bool all_valid_7 = all_valid_6 & (!_S563);
        for(;;)
        {
            _S564 = &proj_points_2[int(5)];
            for(;;)
            {
                float2  _S607 = float2 {sigmas_2->p_0[int(5)].x, sigmas_2->p_0[int(5)].y};
                float r_15 = length_0(_S607);
                float _S608 = sigmas_2->p_0[int(5)].z;
                float theta_11 = (F32_atan2((r_15), (_S608)));
                if(theta_11 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_11 * theta_11 / 3.0f) / _S608;
                }
                else
                {
                    k_2 = theta_11 / r_15;
                }
                float2  _S609 = _S607 * make_float2 (k_2);
                proj_points_2[int(5)] = _S609;
                bool _S610 = is_valid_distortion(_S609, dist_coeffs_13);
                bool _S611 = !_S610;
                _S565 = _S611;
                if(_S611)
                {
                    break;
                }
                float u_29 = proj_points_2[int(5)].x;
                float v_29 = proj_points_2[int(5)].y;
                float r2_29 = u_29 * u_29 + v_29 * v_29;
                float2  _S612 = proj_points_2[int(5)] * make_float2 (1.0f + r2_29 * ((*dist_coeffs_13)[int(0)] + r2_29 * ((*dist_coeffs_13)[int(1)] + r2_29 * ((*dist_coeffs_13)[int(2)] + r2_29 * (*dist_coeffs_13)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_13)[int(4)] * u_29 * v_29 + (*dist_coeffs_13)[int(5)] * (r2_29 + 2.0f * u_29 * u_29) + (*dist_coeffs_13)[int(6)] * r2_29, 2.0f * (*dist_coeffs_13)[int(5)] * u_29 * v_29 + (*dist_coeffs_13)[int(4)] * (r2_29 + 2.0f * v_29 * v_29) + (*dist_coeffs_13)[int(7)] * r2_29);
                float2  _S613 = _S612 + make_float2 ((*dist_coeffs_13)[int(8)] * _S612.x + (*dist_coeffs_13)[int(9)] * _S612.y, 0.0f);
                proj_points_2[int(5)] = make_float2 (intrins_4.x * _S613.x + intrins_4.z, intrins_4.y * _S613.y + intrins_4.w);
                break;
            }
            _S566 = all_valid_7 & (!_S565);
            break;
        }
        _S567 = &proj_points_2[int(6)];
        for(;;)
        {
            float2  _S614 = float2 {sigmas_2->p_0[int(6)].x, sigmas_2->p_0[int(6)].y};
            float r_16 = length_0(_S614);
            float _S615 = sigmas_2->p_0[int(6)].z;
            float theta_12 = (F32_atan2((r_16), (_S615)));
            if(theta_12 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_12 * theta_12 / 3.0f) / _S615;
            }
            else
            {
                k_2 = theta_12 / r_16;
            }
            float2  _S616 = _S614 * make_float2 (k_2);
            proj_points_2[int(6)] = _S616;
            bool _S617 = is_valid_distortion(_S616, dist_coeffs_13);
            bool _S618 = !_S617;
            _S568 = _S618;
            if(_S618)
            {
                break;
            }
            float u_30 = proj_points_2[int(6)].x;
            float v_30 = proj_points_2[int(6)].y;
            float r2_30 = u_30 * u_30 + v_30 * v_30;
            float2  _S619 = proj_points_2[int(6)] * make_float2 (1.0f + r2_30 * ((*dist_coeffs_13)[int(0)] + r2_30 * ((*dist_coeffs_13)[int(1)] + r2_30 * ((*dist_coeffs_13)[int(2)] + r2_30 * (*dist_coeffs_13)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_13)[int(4)] * u_30 * v_30 + (*dist_coeffs_13)[int(5)] * (r2_30 + 2.0f * u_30 * u_30) + (*dist_coeffs_13)[int(6)] * r2_30, 2.0f * (*dist_coeffs_13)[int(5)] * u_30 * v_30 + (*dist_coeffs_13)[int(4)] * (r2_30 + 2.0f * v_30 * v_30) + (*dist_coeffs_13)[int(7)] * r2_30);
            float2  _S620 = _S619 + make_float2 ((*dist_coeffs_13)[int(8)] * _S619.x + (*dist_coeffs_13)[int(9)] * _S619.y, 0.0f);
            proj_points_2[int(6)] = make_float2 (intrins_4.x * _S620.x + intrins_4.z, intrins_4.y * _S620.y + intrins_4.w);
            break;
        }
        _S569 = _S566 & (!_S568);
        break;
    }
    if(!_S569)
    {
        return false;
    }
    float2  _S621 = *mean2d_5 + make_float2 (sigmas_2->w_mean_0[int(0)]) * *_S553;
    *mean2d_5 = _S621;
    float2  _S622 = _S621 + make_float2 (sigmas_2->w_mean_0[int(1)]) * *_S555;
    *mean2d_5 = _S622;
    float2  _S623 = _S622 + make_float2 (sigmas_2->w_mean_0[int(2)]) * *_S557;
    *mean2d_5 = _S623;
    float2  _S624 = _S623 + make_float2 (sigmas_2->w_mean_0[int(3)]) * *_S560;
    *mean2d_5 = _S624;
    float2  _S625 = _S624 + make_float2 (sigmas_2->w_mean_0[int(4)]) * *_S562;
    *mean2d_5 = _S625;
    float2  _S626 = _S625 + make_float2 (sigmas_2->w_mean_0[int(5)]) * *_S564;
    *mean2d_5 = _S626;
    float2  _S627 = _S626 + make_float2 (sigmas_2->w_mean_0[int(6)]) * *_S567;
    *mean2d_5 = _S627;
    float2  d_14 = *_S553 - _S627;
    float _S628 = d_14.x;
    float _S629 = d_14.y;
    float _S630 = _S628 * _S629;
    Matrix<float, 2, 2>  _S631 = *cov2d_5 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S628 * _S628, _S630, _S630, _S629 * _S629);
    *cov2d_5 = _S631;
    float2  d_15 = *_S555 - *mean2d_5;
    float _S632 = d_15.x;
    float _S633 = d_15.y;
    float _S634 = _S632 * _S633;
    Matrix<float, 2, 2>  _S635 = _S631 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S632 * _S632, _S634, _S634, _S633 * _S633);
    *cov2d_5 = _S635;
    float2  d_16 = *_S557 - *mean2d_5;
    float _S636 = d_16.x;
    float _S637 = d_16.y;
    float _S638 = _S636 * _S637;
    Matrix<float, 2, 2>  _S639 = _S635 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S636 * _S636, _S638, _S638, _S637 * _S637);
    *cov2d_5 = _S639;
    float2  d_17 = *_S560 - *mean2d_5;
    float _S640 = d_17.x;
    float _S641 = d_17.y;
    float _S642 = _S640 * _S641;
    Matrix<float, 2, 2>  _S643 = _S639 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S640 * _S640, _S642, _S642, _S641 * _S641);
    *cov2d_5 = _S643;
    float2  d_18 = *_S562 - *mean2d_5;
    float _S644 = d_18.x;
    float _S645 = d_18.y;
    float _S646 = _S644 * _S645;
    Matrix<float, 2, 2>  _S647 = _S643 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S644 * _S644, _S646, _S646, _S645 * _S645);
    *cov2d_5 = _S647;
    float2  d_19 = *_S564 - *mean2d_5;
    float _S648 = d_19.x;
    float _S649 = d_19.y;
    float _S650 = _S648 * _S649;
    Matrix<float, 2, 2>  _S651 = _S647 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S648 * _S648, _S650, _S650, _S649 * _S649);
    *cov2d_5 = _S651;
    float2  d_20 = *_S567 - *mean2d_5;
    float _S652 = d_20.x;
    float _S653 = d_20.y;
    float _S654 = _S652 * _S653;
    *cov2d_5 = _S651 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S652 * _S652, _S654, _S654, _S653 * _S653);
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_3, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_14, Matrix<float, 2, 2>  * cov2d_6, float2  * mean2d_6)
{
    int2  _S655 = make_int2 (int(0));
    float2  _S656 = make_float2 ((float)_S655.x, (float)_S655.y);
    *mean2d_6 = _S656;
    *cov2d_6 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_3;
    float2  _S657 = float2 {sigmas_3->p_0[int(0)].x, sigmas_3->p_0[int(0)].y};
    float r_17 = length_0(_S657);
    float _S658 = sigmas_3->p_0[int(0)].z;
    float theta_13 = (F32_atan2((r_17), (_S658)));
    float k_3;
    if(theta_13 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_13 * theta_13 / 3.0f) / _S658;
    }
    else
    {
        k_3 = theta_13 / r_17;
    }
    float2  _S659 = _S657 * make_float2 (k_3);
    float u_31 = _S659.x;
    float v_31 = _S659.y;
    float r2_31 = u_31 * u_31 + v_31 * v_31;
    float _S660 = 2.0f * (*dist_coeffs_14)[int(4)];
    float _S661 = 2.0f * (*dist_coeffs_14)[int(5)];
    float2  _S662 = _S659 * make_float2 (1.0f + r2_31 * ((*dist_coeffs_14)[int(0)] + r2_31 * ((*dist_coeffs_14)[int(1)] + r2_31 * ((*dist_coeffs_14)[int(2)] + r2_31 * (*dist_coeffs_14)[int(3)])))) + make_float2 (_S660 * u_31 * v_31 + (*dist_coeffs_14)[int(5)] * (r2_31 + 2.0f * u_31 * u_31) + (*dist_coeffs_14)[int(6)] * r2_31, _S661 * u_31 * v_31 + (*dist_coeffs_14)[int(4)] * (r2_31 + 2.0f * v_31 * v_31) + (*dist_coeffs_14)[int(7)] * r2_31);
    float2  _S663 = _S662 + make_float2 ((*dist_coeffs_14)[int(8)] * _S662.x + (*dist_coeffs_14)[int(9)] * _S662.y, 0.0f);
    float fx_5 = intrins_5.x;
    float fy_5 = intrins_5.y;
    float cx_2 = intrins_5.z;
    float cy_2 = intrins_5.w;
    proj_points_3[int(0)] = make_float2 (fx_5 * _S663.x + cx_2, fy_5 * _S663.y + cy_2);
    float2  _S664 = float2 {sigmas_3->p_0[int(1)].x, sigmas_3->p_0[int(1)].y};
    float r_18 = length_0(_S664);
    float _S665 = sigmas_3->p_0[int(1)].z;
    float theta_14 = (F32_atan2((r_18), (_S665)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_14 * theta_14 / 3.0f) / _S665;
    }
    else
    {
        k_3 = theta_14 / r_18;
    }
    float2  _S666 = _S664 * make_float2 (k_3);
    float u_32 = _S666.x;
    float v_32 = _S666.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float2  _S667 = _S666 * make_float2 (1.0f + r2_32 * ((*dist_coeffs_14)[int(0)] + r2_32 * ((*dist_coeffs_14)[int(1)] + r2_32 * ((*dist_coeffs_14)[int(2)] + r2_32 * (*dist_coeffs_14)[int(3)])))) + make_float2 (_S660 * u_32 * v_32 + (*dist_coeffs_14)[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + (*dist_coeffs_14)[int(6)] * r2_32, _S661 * u_32 * v_32 + (*dist_coeffs_14)[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + (*dist_coeffs_14)[int(7)] * r2_32);
    float2  _S668 = _S667 + make_float2 ((*dist_coeffs_14)[int(8)] * _S667.x + (*dist_coeffs_14)[int(9)] * _S667.y, 0.0f);
    proj_points_3[int(1)] = make_float2 (fx_5 * _S668.x + cx_2, fy_5 * _S668.y + cy_2);
    float2  _S669 = float2 {sigmas_3->p_0[int(2)].x, sigmas_3->p_0[int(2)].y};
    float r_19 = length_0(_S669);
    float _S670 = sigmas_3->p_0[int(2)].z;
    float theta_15 = (F32_atan2((r_19), (_S670)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_15 * theta_15 / 3.0f) / _S670;
    }
    else
    {
        k_3 = theta_15 / r_19;
    }
    float2  _S671 = _S669 * make_float2 (k_3);
    float u_33 = _S671.x;
    float v_33 = _S671.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S672 = _S671 * make_float2 (1.0f + r2_33 * ((*dist_coeffs_14)[int(0)] + r2_33 * ((*dist_coeffs_14)[int(1)] + r2_33 * ((*dist_coeffs_14)[int(2)] + r2_33 * (*dist_coeffs_14)[int(3)])))) + make_float2 (_S660 * u_33 * v_33 + (*dist_coeffs_14)[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + (*dist_coeffs_14)[int(6)] * r2_33, _S661 * u_33 * v_33 + (*dist_coeffs_14)[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + (*dist_coeffs_14)[int(7)] * r2_33);
    float2  _S673 = _S672 + make_float2 ((*dist_coeffs_14)[int(8)] * _S672.x + (*dist_coeffs_14)[int(9)] * _S672.y, 0.0f);
    proj_points_3[int(2)] = make_float2 (fx_5 * _S673.x + cx_2, fy_5 * _S673.y + cy_2);
    float2  _S674 = float2 {sigmas_3->p_0[int(3)].x, sigmas_3->p_0[int(3)].y};
    float r_20 = length_0(_S674);
    float _S675 = sigmas_3->p_0[int(3)].z;
    float theta_16 = (F32_atan2((r_20), (_S675)));
    if(theta_16 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_16 * theta_16 / 3.0f) / _S675;
    }
    else
    {
        k_3 = theta_16 / r_20;
    }
    float2  _S676 = _S674 * make_float2 (k_3);
    float u_34 = _S676.x;
    float v_34 = _S676.y;
    float r2_34 = u_34 * u_34 + v_34 * v_34;
    float2  _S677 = _S676 * make_float2 (1.0f + r2_34 * ((*dist_coeffs_14)[int(0)] + r2_34 * ((*dist_coeffs_14)[int(1)] + r2_34 * ((*dist_coeffs_14)[int(2)] + r2_34 * (*dist_coeffs_14)[int(3)])))) + make_float2 (_S660 * u_34 * v_34 + (*dist_coeffs_14)[int(5)] * (r2_34 + 2.0f * u_34 * u_34) + (*dist_coeffs_14)[int(6)] * r2_34, _S661 * u_34 * v_34 + (*dist_coeffs_14)[int(4)] * (r2_34 + 2.0f * v_34 * v_34) + (*dist_coeffs_14)[int(7)] * r2_34);
    float2  _S678 = _S677 + make_float2 ((*dist_coeffs_14)[int(8)] * _S677.x + (*dist_coeffs_14)[int(9)] * _S677.y, 0.0f);
    proj_points_3[int(3)] = make_float2 (fx_5 * _S678.x + cx_2, fy_5 * _S678.y + cy_2);
    float2  _S679 = float2 {sigmas_3->p_0[int(4)].x, sigmas_3->p_0[int(4)].y};
    float r_21 = length_0(_S679);
    float _S680 = sigmas_3->p_0[int(4)].z;
    float theta_17 = (F32_atan2((r_21), (_S680)));
    if(theta_17 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_17 * theta_17 / 3.0f) / _S680;
    }
    else
    {
        k_3 = theta_17 / r_21;
    }
    float2  _S681 = _S679 * make_float2 (k_3);
    float u_35 = _S681.x;
    float v_35 = _S681.y;
    float r2_35 = u_35 * u_35 + v_35 * v_35;
    float2  _S682 = _S681 * make_float2 (1.0f + r2_35 * ((*dist_coeffs_14)[int(0)] + r2_35 * ((*dist_coeffs_14)[int(1)] + r2_35 * ((*dist_coeffs_14)[int(2)] + r2_35 * (*dist_coeffs_14)[int(3)])))) + make_float2 (_S660 * u_35 * v_35 + (*dist_coeffs_14)[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + (*dist_coeffs_14)[int(6)] * r2_35, _S661 * u_35 * v_35 + (*dist_coeffs_14)[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + (*dist_coeffs_14)[int(7)] * r2_35);
    float2  _S683 = _S682 + make_float2 ((*dist_coeffs_14)[int(8)] * _S682.x + (*dist_coeffs_14)[int(9)] * _S682.y, 0.0f);
    proj_points_3[int(4)] = make_float2 (fx_5 * _S683.x + cx_2, fy_5 * _S683.y + cy_2);
    float2  _S684 = float2 {sigmas_3->p_0[int(5)].x, sigmas_3->p_0[int(5)].y};
    float r_22 = length_0(_S684);
    float _S685 = sigmas_3->p_0[int(5)].z;
    float theta_18 = (F32_atan2((r_22), (_S685)));
    if(theta_18 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_18 * theta_18 / 3.0f) / _S685;
    }
    else
    {
        k_3 = theta_18 / r_22;
    }
    float2  _S686 = _S684 * make_float2 (k_3);
    float u_36 = _S686.x;
    float v_36 = _S686.y;
    float r2_36 = u_36 * u_36 + v_36 * v_36;
    float2  _S687 = _S686 * make_float2 (1.0f + r2_36 * ((*dist_coeffs_14)[int(0)] + r2_36 * ((*dist_coeffs_14)[int(1)] + r2_36 * ((*dist_coeffs_14)[int(2)] + r2_36 * (*dist_coeffs_14)[int(3)])))) + make_float2 (_S660 * u_36 * v_36 + (*dist_coeffs_14)[int(5)] * (r2_36 + 2.0f * u_36 * u_36) + (*dist_coeffs_14)[int(6)] * r2_36, _S661 * u_36 * v_36 + (*dist_coeffs_14)[int(4)] * (r2_36 + 2.0f * v_36 * v_36) + (*dist_coeffs_14)[int(7)] * r2_36);
    float2  _S688 = _S687 + make_float2 ((*dist_coeffs_14)[int(8)] * _S687.x + (*dist_coeffs_14)[int(9)] * _S687.y, 0.0f);
    proj_points_3[int(5)] = make_float2 (fx_5 * _S688.x + cx_2, fy_5 * _S688.y + cy_2);
    float2  _S689 = float2 {sigmas_3->p_0[int(6)].x, sigmas_3->p_0[int(6)].y};
    float r_23 = length_0(_S689);
    float _S690 = sigmas_3->p_0[int(6)].z;
    float theta_19 = (F32_atan2((r_23), (_S690)));
    if(theta_19 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_19 * theta_19 / 3.0f) / _S690;
    }
    else
    {
        k_3 = theta_19 / r_23;
    }
    float2  _S691 = _S689 * make_float2 (k_3);
    float u_37 = _S691.x;
    float v_37 = _S691.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float2  _S692 = _S691 * make_float2 (1.0f + r2_37 * ((*dist_coeffs_14)[int(0)] + r2_37 * ((*dist_coeffs_14)[int(1)] + r2_37 * ((*dist_coeffs_14)[int(2)] + r2_37 * (*dist_coeffs_14)[int(3)])))) + make_float2 (_S660 * u_37 * v_37 + (*dist_coeffs_14)[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + (*dist_coeffs_14)[int(6)] * r2_37, _S661 * u_37 * v_37 + (*dist_coeffs_14)[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + (*dist_coeffs_14)[int(7)] * r2_37);
    float2  _S693 = _S692 + make_float2 ((*dist_coeffs_14)[int(8)] * _S692.x + (*dist_coeffs_14)[int(9)] * _S692.y, 0.0f);
    float2  _S694 = make_float2 (fx_5 * _S693.x + cx_2, fy_5 * _S693.y + cy_2);
    proj_points_3[int(6)] = _S694;
    float2  _S695 = *mean2d_6 + make_float2 (sigmas_3->w_mean_0[int(0)]) * proj_points_3[int(0)];
    *mean2d_6 = _S695;
    float2  _S696 = _S695 + make_float2 (sigmas_3->w_mean_0[int(1)]) * proj_points_3[int(1)];
    *mean2d_6 = _S696;
    float2  _S697 = _S696 + make_float2 (sigmas_3->w_mean_0[int(2)]) * proj_points_3[int(2)];
    *mean2d_6 = _S697;
    float2  _S698 = _S697 + make_float2 (sigmas_3->w_mean_0[int(3)]) * proj_points_3[int(3)];
    *mean2d_6 = _S698;
    float2  _S699 = _S698 + make_float2 (sigmas_3->w_mean_0[int(4)]) * proj_points_3[int(4)];
    *mean2d_6 = _S699;
    float2  _S700 = _S699 + make_float2 (sigmas_3->w_mean_0[int(5)]) * proj_points_3[int(5)];
    *mean2d_6 = _S700;
    float2  _S701 = _S700 + make_float2 (sigmas_3->w_mean_0[int(6)]) * _S694;
    *mean2d_6 = _S701;
    float2  d_21 = proj_points_3[int(0)] - _S701;
    float _S702 = d_21.x;
    float _S703 = d_21.y;
    float _S704 = _S702 * _S703;
    Matrix<float, 2, 2>  _S705 = *cov2d_6 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S702 * _S702, _S704, _S704, _S703 * _S703);
    *cov2d_6 = _S705;
    float2  d_22 = proj_points_3[int(1)] - *mean2d_6;
    float _S706 = d_22.x;
    float _S707 = d_22.y;
    float _S708 = _S706 * _S707;
    Matrix<float, 2, 2>  _S709 = _S705 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S706 * _S706, _S708, _S708, _S707 * _S707);
    *cov2d_6 = _S709;
    float2  d_23 = proj_points_3[int(2)] - *mean2d_6;
    float _S710 = d_23.x;
    float _S711 = d_23.y;
    float _S712 = _S710 * _S711;
    Matrix<float, 2, 2>  _S713 = _S709 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S710 * _S710, _S712, _S712, _S711 * _S711);
    *cov2d_6 = _S713;
    float2  d_24 = proj_points_3[int(3)] - *mean2d_6;
    float _S714 = d_24.x;
    float _S715 = d_24.y;
    float _S716 = _S714 * _S715;
    Matrix<float, 2, 2>  _S717 = _S713 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S714 * _S714, _S716, _S716, _S715 * _S715);
    *cov2d_6 = _S717;
    float2  d_25 = proj_points_3[int(4)] - *mean2d_6;
    float _S718 = d_25.x;
    float _S719 = d_25.y;
    float _S720 = _S718 * _S719;
    Matrix<float, 2, 2>  _S721 = _S717 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S718 * _S718, _S720, _S720, _S719 * _S719);
    *cov2d_6 = _S721;
    float2  d_26 = proj_points_3[int(5)] - *mean2d_6;
    float _S722 = d_26.x;
    float _S723 = d_26.y;
    float _S724 = _S722 * _S723;
    Matrix<float, 2, 2>  _S725 = _S721 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S722 * _S722, _S724, _S724, _S723 * _S723);
    *cov2d_6 = _S725;
    float2  d_27 = _S694 - *mean2d_6;
    float _S726 = d_27.x;
    float _S727 = d_27.y;
    float _S728 = _S726 * _S727;
    *cov2d_6 = _S725 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S726 * _S726, _S728, _S728, _S727 * _S727);
    return true;
}

inline __device__ void ortho_proj_3dgs(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float fx_6, float fy_6, float cx_3, float cy_3, Matrix<float, 2, 2>  * cov2d_7, float2  * mean2d_7)
{
    Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
    *cov2d_7 = mul_6(mul_5(J_5, cov3d_3), transpose_1(J_5));
    *mean2d_7 = make_float2 (fx_6 * mean3d_3.x + cx_3, fy_6 * mean3d_3.y + cy_3);
    return;
}

inline __device__ float add_blur(float eps2d_0, Matrix<float, 2, 2>  * covar_1, float * compensation_0)
{
    float det_orig_0 = *&((covar_1->rows + (int(0)))->x) * *&((covar_1->rows + (int(1)))->y) - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    float _S729 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S729;
    float _S730 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S730;
    float det_blur_0 = _S729 * _S730 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_10, float dOut_12)
{
    float _S731 = (F32_exp(((*dpx_10).primal_0))) * dOut_12;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S731;
    return;
}

inline __device__ float3  exp_0(float3  x_13)
{
    float3  result_12;
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
        *_slang_vector_get_element_ptr(&result_12, i_8) = (F32_exp((_slang_vector_get_element(x_13, i_8))));
        i_8 = i_8 + int(1);
    }
    return result_12;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, float3  dOut_13)
{
    float3  _S732 = exp_0((*dpx_11).primal_0) * dOut_13;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S732;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_12, float dOut_14)
{
    float _S733 = 1.0f / (*dpx_12).primal_0 * dOut_14;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S733;
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

inline __device__ float3  max_0(float3  x_14, float3  y_2)
{
    float3  result_13;
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
        *_slang_vector_get_element_ptr(&result_13, i_9) = (F32_max((_slang_vector_get_element(x_14, i_9)), (_slang_vector_get_element(y_2, i_9))));
        i_9 = i_9 + int(1);
    }
    return result_13;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_0, float4  quat_3, float3  scale_2, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_7, float fy_7, float cx_4, float cy_4, FixedArray<float, 10>  * dist_coeffs_15, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_8, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_4, mean_0) + t_3;
        float _S734 = mean_c_0.z;
        bool _S735;
        if(_S734 < near_plane_0)
        {
            _S735 = true;
        }
        else
        {
            _S735 = _S734 > far_plane_0;
        }
        if(_S735)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_0;
        float3  _S736 = exp_0(scale_2);
        float x_15 = quat_3.y;
        float x2_3 = x_15 * x_15;
        float y2_3 = quat_3.z * quat_3.z;
        float z2_3 = quat_3.w * quat_3.w;
        float xy_3 = quat_3.y * quat_3.z;
        float xz_3 = quat_3.y * quat_3.w;
        float yz_3 = quat_3.z * quat_3.w;
        float wx_3 = quat_3.x * quat_3.y;
        float wy_3 = quat_3.x * quat_3.z;
        float wz_3 = quat_3.x * quat_3.w;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3))), makeMatrix<float, 3, 3> (_S736.x, 0.0f, 0.0f, 0.0f, _S736.y, 0.0f, 0.0f, 0.0f, _S736.z));
        Matrix<float, 3, 3>  _S737 = transpose_0(R_4);
        float _S738 = float(image_width_0);
        float _S739 = float(image_height_0);
        float _S740 = 0.30000001192092896f * (0.5f * _S738 / fx_7);
        float _S741 = 0.30000001192092896f * (0.5f * _S739 / fy_7);
        float rz_1 = 1.0f / mean_c_0.z;
        float rz2_1 = rz_1 * rz_1;
        Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_7 * rz_1, 0.0f, - fx_7 * (mean_c_0.z * (F32_min(((_S738 - cx_4) / fx_7 + _S740), ((F32_max((- (cx_4 / fx_7 + _S740)), (mean_c_0.x * rz_1))))))) * rz2_1, 0.0f, fy_7 * rz_1, - fy_7 * (mean_c_0.z * (F32_min(((_S739 - cy_4) / fy_7 + _S741), ((F32_max((- (cy_4 / fy_7 + _S741)), (mean_c_0.y * rz_1))))))) * rz2_1);
        covar2d_0 = mul_6(mul_5(J_6, mul_4(mul_4(R_4, mul_4(M_2, transpose_0(M_2))), _S737)), transpose_1(J_6));
        *mean2d_8 = make_float2 (fx_7 * mean_c_0.x * rz_1 + cx_4, fy_7 * mean_c_0.y * rz_1 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S742 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S742;
        float _S743 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S743;
        float det_blur_1 = _S742 * _S743 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S744 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
        float xmin_0 = (F32_floor(((*mean2d_8).x - radius_x_0)));
        float xmax_0 = (F32_ceil(((*mean2d_8).x + radius_x_0)));
        float ymin_0 = (F32_floor(((*mean2d_8).y - radius_y_0)));
        float ymax_0 = (F32_ceil(((*mean2d_8).y + radius_y_0)));
        if(xmax_0 <= 0.0f)
        {
            _S735 = true;
        }
        else
        {
            _S735 = xmin_0 >= _S738;
        }
        if(_S735)
        {
            _S735 = true;
        }
        else
        {
            _S735 = ymax_0 <= 0.0f;
        }
        if(_S735)
        {
            _S735 = true;
        }
        else
        {
            _S735 = ymin_0 >= _S739;
        }
        if(_S735)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = 0.5f * (F32_log((dot_0(mean_c_0, mean_c_0) + 9.99999997475242708e-07f)));
        *conic_0 = make_float3 (_S744.rows[int(0)].x, _S744.rows[int(0)].y, _S744.rows[int(1)].y);
        float3  _S745 = mean_0 - - mul_0(_S737, t_3);
        float3  _S746 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S746;
        float _S747 = _S745.x;
        float _S748 = _S745.y;
        float _S749 = _S745.z;
        float norm_0 = (F32_sqrt((_S747 * _S747 + _S748 * _S748 + _S749 * _S749)));
        float x_16 = _S747 / norm_0;
        float y_3 = _S748 / norm_0;
        float z_0 = _S749 / norm_0;
        float3  _S750 = _S746 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_16) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S750;
        float z2_4 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_16 * x_16 - y_3 * y_3;
        float fS1_0 = 2.0f * x_16 * y_3;
        float3  _S751 = _S750 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_4 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_16) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S751;
        float fTmp0C_0 = -2.28522896766662598f * z2_4 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgb_0 = max_0(_S751 + (make_float3 (-0.59004360437393188f * (x_16 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_4 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_16) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_16 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_8, float fy_8, float cx_5, float cy_5, FixedArray<float, 10>  * dist_coeffs_16, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_9, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_5, mean_1) + t_4;
        float _S752 = length_1(mean_c_1);
        bool _S753;
        if(_S752 < near_plane_1)
        {
            _S753 = true;
        }
        else
        {
            _S753 = _S752 > far_plane_1;
        }
        if(_S753)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_6 = make_float4 (fx_8, fy_8, cx_5, cy_5);
        Matrix<float, 2, 2>  covar2d_1;
        float3  _S754 = exp_0(scale_3);
        float x_17 = quat_4.y;
        float x2_4 = x_17 * x_17;
        float y2_4 = quat_4.z * quat_4.z;
        float z2_5 = quat_4.w * quat_4.w;
        float xy_4 = quat_4.y * quat_4.z;
        float xz_4 = quat_4.y * quat_4.w;
        float yz_4 = quat_4.z * quat_4.w;
        float wx_4 = quat_4.x * quat_4.y;
        float wy_4 = quat_4.x * quat_4.z;
        float wz_4 = quat_4.x * quat_4.w;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_5), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_5), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S754.x, 0.0f, 0.0f, 0.0f, _S754.y, 0.0f, 0.0f, 0.0f, _S754.z));
        Matrix<float, 3, 3>  _S755 = transpose_0(R_5);
        bool _S756 = fisheye_proj_3dgs_0(mean_c_1, mul_4(mul_4(R_5, mul_4(M_3, transpose_0(M_3))), _S755), intrins_6, dist_coeffs_16, &covar2d_1, mean2d_9);
        if(!(true & _S756))
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S757 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S757;
        float _S758 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S758;
        float det_blur_2 = _S757 * _S758 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S759 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
        float xmin_1 = (F32_floor(((*mean2d_9).x - radius_x_1)));
        float xmax_1 = (F32_ceil(((*mean2d_9).x + radius_x_1)));
        float ymin_1 = (F32_floor(((*mean2d_9).y - radius_y_1)));
        float ymax_1 = (F32_ceil(((*mean2d_9).y + radius_y_1)));
        if(xmax_1 <= 0.0f)
        {
            _S753 = true;
        }
        else
        {
            _S753 = xmin_1 >= float(image_width_1);
        }
        if(_S753)
        {
            _S753 = true;
        }
        else
        {
            _S753 = ymax_1 <= 0.0f;
        }
        if(_S753)
        {
            _S753 = true;
        }
        else
        {
            _S753 = ymin_1 >= float(image_height_1);
        }
        if(_S753)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = 0.5f * (F32_log((dot_0(mean_c_1, mean_c_1) + 9.99999997475242708e-07f)));
        *conic_1 = make_float3 (_S759.rows[int(0)].x, _S759.rows[int(0)].y, _S759.rows[int(1)].y);
        float3  _S760 = mean_1 - - mul_0(_S755, t_4);
        float3  _S761 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S761;
        float _S762 = _S760.x;
        float _S763 = _S760.y;
        float _S764 = _S760.z;
        float norm_1 = (F32_sqrt((_S762 * _S762 + _S763 * _S763 + _S764 * _S764)));
        float x_18 = _S762 / norm_1;
        float y_4 = _S763 / norm_1;
        float z_1 = _S764 / norm_1;
        float3  _S765 = _S761 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_1)[int(1)] + make_float3 (z_1) * (*sh_coeffs_1)[int(2)] - make_float3 (x_18) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S765;
        float z2_6 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_18 * x_18 - y_4 * y_4;
        float fS1_1 = 2.0f * x_18 * y_4;
        float3  _S766 = _S765 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_4) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_6 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_18) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S766;
        float fTmp0C_1 = -2.28522896766662598f * z2_6 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgb_1 = max_0(_S766 + (make_float3 (-0.59004360437393188f * (x_18 * fS1_1 + y_4 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_4) * (*sh_coeffs_1)[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_6 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_18) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_18 * fC1_1 - y_4 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_9, float fy_9, float cx_6, float cy_6, FixedArray<float, 10>  * dist_coeffs_17, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_10, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_6, mean_2) + t_5;
        float _S767 = mean_c_2.z;
        bool _S768;
        if(_S767 < near_plane_2)
        {
            _S768 = true;
        }
        else
        {
            _S768 = _S767 > far_plane_2;
        }
        if(_S768)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_2;
        float3  _S769 = exp_0(scale_4);
        float x_19 = quat_5.y;
        float x2_5 = x_19 * x_19;
        float y2_5 = quat_5.z * quat_5.z;
        float z2_7 = quat_5.w * quat_5.w;
        float xy_5 = quat_5.y * quat_5.z;
        float xz_5 = quat_5.y * quat_5.w;
        float yz_5 = quat_5.z * quat_5.w;
        float wx_5 = quat_5.x * quat_5.y;
        float wy_5 = quat_5.x * quat_5.z;
        float wz_5 = quat_5.x * quat_5.w;
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_7), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_7), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S769.x, 0.0f, 0.0f, 0.0f, _S769.y, 0.0f, 0.0f, 0.0f, _S769.z));
        Matrix<float, 3, 3>  _S770 = transpose_0(R_6);
        Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
        covar2d_2 = mul_6(mul_5(J_7, mul_4(mul_4(R_6, mul_4(M_4, transpose_0(M_4))), _S770)), transpose_1(J_7));
        *mean2d_10 = make_float2 (fx_9 * mean_c_2.x + cx_6, fy_9 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S771 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S771;
        float _S772 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S772;
        float det_blur_3 = _S771 * _S772 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S773 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
        float xmin_2 = (F32_floor(((*mean2d_10).x - radius_x_2)));
        float xmax_2 = (F32_ceil(((*mean2d_10).x + radius_x_2)));
        float ymin_2 = (F32_floor(((*mean2d_10).y - radius_y_2)));
        float ymax_2 = (F32_ceil(((*mean2d_10).y + radius_y_2)));
        if(xmax_2 <= 0.0f)
        {
            _S768 = true;
        }
        else
        {
            _S768 = xmin_2 >= float(image_width_2);
        }
        if(_S768)
        {
            _S768 = true;
        }
        else
        {
            _S768 = ymax_2 <= 0.0f;
        }
        if(_S768)
        {
            _S768 = true;
        }
        else
        {
            _S768 = ymin_2 >= float(image_height_2);
        }
        if(_S768)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = 0.5f * (F32_log((dot_0(mean_c_2, mean_c_2) + 9.99999997475242708e-07f)));
        *conic_2 = make_float3 (_S773.rows[int(0)].x, _S773.rows[int(0)].y, _S773.rows[int(1)].y);
        float3  _S774 = mean_2 - - mul_0(_S770, t_5);
        float3  _S775 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S775;
        float _S776 = _S774.x;
        float _S777 = _S774.y;
        float _S778 = _S774.z;
        float norm_2 = (F32_sqrt((_S776 * _S776 + _S777 * _S777 + _S778 * _S778)));
        float x_20 = _S776 / norm_2;
        float y_5 = _S777 / norm_2;
        float z_2 = _S778 / norm_2;
        float3  _S779 = _S775 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_2)[int(1)] + make_float3 (z_2) * (*sh_coeffs_2)[int(2)] - make_float3 (x_20) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S779;
        float z2_8 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_20 * x_20 - y_5 * y_5;
        float fS1_2 = 2.0f * x_20 * y_5;
        float3  _S780 = _S779 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_5) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_8 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_20) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S780;
        float fTmp0C_2 = -2.28522896766662598f * z2_8 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        *rgb_2 = max_0(_S780 + (make_float3 (-0.59004360437393188f * (x_20 * fS1_2 + y_5 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_5) * (*sh_coeffs_2)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_8 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_20) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_20 * fC1_2 - y_5 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_persp(bool antialiased_3, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_10, float fy_10, float cx_7, float cy_7, FixedArray<float, 10>  * dist_coeffs_18, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_11, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_0(R_7, mean_3) + t_6;
        float _S781 = mean_c_3.z;
        bool _S782;
        if(_S781 < near_plane_3)
        {
            _S782 = true;
        }
        else
        {
            _S782 = _S781 > far_plane_3;
        }
        if(_S782)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_7 = make_float4 (fx_10, fy_10, cx_7, cy_7);
        Matrix<float, 2, 2>  covar2d_3;
        float3  _S783 = exp_0(scale_5);
        float x_21 = quat_6.y;
        float x2_6 = x_21 * x_21;
        float y2_6 = quat_6.z * quat_6.z;
        float z2_9 = quat_6.w * quat_6.w;
        float xy_6 = quat_6.y * quat_6.z;
        float xz_6 = quat_6.y * quat_6.w;
        float yz_6 = quat_6.z * quat_6.w;
        float wx_6 = quat_6.x * quat_6.y;
        float wy_6 = quat_6.x * quat_6.z;
        float wz_6 = quat_6.x * quat_6.w;
        Matrix<float, 3, 3>  _S784 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_9), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_9), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))));
        SigmaPoints_0 ret_0;
        (&ret_0)->p_0[int(0)] = mean_3;
        (&ret_0)->w_mean_0[int(0)] = 0.0f;
        (&ret_0)->w_cov_0[int(0)] = 2.0f;
        float _S785 = (F32_sqrt((3.0f)));
        float3  delta_0 = make_float3 (_S785 * _S783.x) * _S784.rows[0U];
        float3  _S786 = mean_3 + delta_0;
        float3  _S787 = mean_3 - delta_0;
        float3  delta_1 = make_float3 (_S785 * _S783.y) * _S784.rows[1U];
        float3  _S788 = mean_3 + delta_1;
        float3  _S789 = mean_3 - delta_1;
        float3  delta_2 = make_float3 (_S785 * _S783.z) * _S784.rows[2U];
        float3  _S790 = mean_3 + delta_2;
        float3  _S791 = mean_3 - delta_2;
        (&ret_0)->w_mean_0[1U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[1U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[2U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[2U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[3U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[3U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[4U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[4U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[5U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[5U] = 0.1666666716337204f;
        (&ret_0)->w_mean_0[6U] = 0.1666666716337204f;
        (&ret_0)->w_cov_0[6U] = 0.1666666716337204f;
        (&ret_0)->p_0[0U] = mul_0(R_7, (&ret_0)->p_0[0U]) + t_6;
        (&ret_0)->p_0[1U] = mul_0(R_7, _S786) + t_6;
        (&ret_0)->p_0[2U] = mul_0(R_7, _S788) + t_6;
        (&ret_0)->p_0[3U] = mul_0(R_7, _S790) + t_6;
        (&ret_0)->p_0[4U] = mul_0(R_7, _S787) + t_6;
        (&ret_0)->p_0[5U] = mul_0(R_7, _S789) + t_6;
        (&ret_0)->p_0[6U] = mul_0(R_7, _S791) + t_6;
        SigmaPoints_0 _S792 = ret_0;
        bool _S793 = persp_proj_3dgs_ut_0(&_S792, intrins_7, dist_coeffs_18, image_width_3, image_height_3, &covar2d_3, mean2d_11);
        if(!(true & _S793))
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S794 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S794;
        float _S795 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S795;
        float det_blur_4 = _S794 * _S795 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
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
        float xmin_3 = (F32_floor(((*mean2d_11).x - radius_x_3)));
        float xmax_3 = (F32_ceil(((*mean2d_11).x + radius_x_3)));
        float ymin_3 = (F32_floor(((*mean2d_11).y - radius_y_3)));
        float ymax_3 = (F32_ceil(((*mean2d_11).y + radius_y_3)));
        if(xmax_3 <= 0.0f)
        {
            _S782 = true;
        }
        else
        {
            _S782 = xmin_3 >= float(image_width_3);
        }
        if(_S782)
        {
            _S782 = true;
        }
        else
        {
            _S782 = ymax_3 <= 0.0f;
        }
        if(_S782)
        {
            _S782 = true;
        }
        else
        {
            _S782 = ymin_3 >= float(image_height_3);
        }
        if(_S782)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_3 = make_int4 (int(xmin_3), int(ymin_3), int(xmax_3), int(ymax_3));
        *depth_3 = 0.5f * (F32_log((dot_0(mean_c_3, mean_c_3) + 9.99999997475242708e-07f)));
        *conic_3 = exp_0(- scale_5);
        float3  _S796 = mean_3 - - mul_0(transpose_0(R_7), t_6);
        float3  _S797 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
        *rgb_3 = _S797;
        float _S798 = _S796.x;
        float _S799 = _S796.y;
        float _S800 = _S796.z;
        float norm_3 = (F32_sqrt((_S798 * _S798 + _S799 * _S799 + _S800 * _S800)));
        float x_22 = _S798 / norm_3;
        float y_6 = _S799 / norm_3;
        float z_3 = _S800 / norm_3;
        float3  _S801 = _S797 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_3)[int(1)] + make_float3 (z_3) * (*sh_coeffs_3)[int(2)] - make_float3 (x_22) * (*sh_coeffs_3)[int(3)]);
        *rgb_3 = _S801;
        float z2_10 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_22 * x_22 - y_6 * y_6;
        float fS1_3 = 2.0f * x_22 * y_6;
        float3  _S802 = _S801 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_6) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_10 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_22) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
        *rgb_3 = _S802;
        float fTmp0C_3 = -2.28522896766662598f * z2_10 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        *rgb_3 = max_0(_S802 + (make_float3 (-0.59004360437393188f * (x_22 * fS1_3 + y_6 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_6) * (*sh_coeffs_3)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_10 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_22) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_22 * fC1_3 - y_6 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_fisheye(bool antialiased_4, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_11, float fy_11, float cx_8, float cy_8, FixedArray<float, 10>  * dist_coeffs_19, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_12, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    for(;;)
    {
        float3  mean_c_4 = mul_0(R_8, mean_4) + t_7;
        float _S803 = length_1(mean_c_4);
        bool _S804;
        if(_S803 < near_plane_4)
        {
            _S804 = true;
        }
        else
        {
            _S804 = _S803 > far_plane_4;
        }
        if(_S804)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_8 = make_float4 (fx_11, fy_11, cx_8, cy_8);
        Matrix<float, 2, 2>  covar2d_4;
        float3  _S805 = exp_0(scale_6);
        float x_23 = quat_7.y;
        float x2_7 = x_23 * x_23;
        float y2_7 = quat_7.z * quat_7.z;
        float z2_11 = quat_7.w * quat_7.w;
        float xy_7 = quat_7.y * quat_7.z;
        float xz_7 = quat_7.y * quat_7.w;
        float yz_7 = quat_7.z * quat_7.w;
        float wx_7 = quat_7.x * quat_7.y;
        float wy_7 = quat_7.x * quat_7.z;
        float wz_7 = quat_7.x * quat_7.w;
        Matrix<float, 3, 3>  _S806 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_11), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_11), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))));
        SigmaPoints_0 ret_1;
        (&ret_1)->p_0[int(0)] = mean_4;
        (&ret_1)->w_mean_0[int(0)] = 0.0f;
        (&ret_1)->w_cov_0[int(0)] = 2.0f;
        float _S807 = (F32_sqrt((3.0f)));
        float3  delta_3 = make_float3 (_S807 * _S805.x) * _S806.rows[0U];
        float3  _S808 = mean_4 + delta_3;
        float3  _S809 = mean_4 - delta_3;
        float3  delta_4 = make_float3 (_S807 * _S805.y) * _S806.rows[1U];
        float3  _S810 = mean_4 + delta_4;
        float3  _S811 = mean_4 - delta_4;
        float3  delta_5 = make_float3 (_S807 * _S805.z) * _S806.rows[2U];
        float3  _S812 = mean_4 + delta_5;
        float3  _S813 = mean_4 - delta_5;
        (&ret_1)->w_mean_0[1U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[1U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[2U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[2U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[3U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[3U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[4U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[4U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[5U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[5U] = 0.1666666716337204f;
        (&ret_1)->w_mean_0[6U] = 0.1666666716337204f;
        (&ret_1)->w_cov_0[6U] = 0.1666666716337204f;
        (&ret_1)->p_0[0U] = mul_0(R_8, (&ret_1)->p_0[0U]) + t_7;
        (&ret_1)->p_0[1U] = mul_0(R_8, _S808) + t_7;
        (&ret_1)->p_0[2U] = mul_0(R_8, _S810) + t_7;
        (&ret_1)->p_0[3U] = mul_0(R_8, _S812) + t_7;
        (&ret_1)->p_0[4U] = mul_0(R_8, _S809) + t_7;
        (&ret_1)->p_0[5U] = mul_0(R_8, _S811) + t_7;
        (&ret_1)->p_0[6U] = mul_0(R_8, _S813) + t_7;
        SigmaPoints_0 _S814 = ret_1;
        bool _S815 = fisheye_proj_3dgs_ut_0(&_S814, intrins_8, dist_coeffs_19, &covar2d_4, mean2d_12);
        if(!(true & _S815))
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S816 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S816;
        float _S817 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S817;
        float det_blur_5 = _S816 * _S817 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
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
        float xmin_4 = (F32_floor(((*mean2d_12).x - radius_x_4)));
        float xmax_4 = (F32_ceil(((*mean2d_12).x + radius_x_4)));
        float ymin_4 = (F32_floor(((*mean2d_12).y - radius_y_4)));
        float ymax_4 = (F32_ceil(((*mean2d_12).y + radius_y_4)));
        if(xmax_4 <= 0.0f)
        {
            _S804 = true;
        }
        else
        {
            _S804 = xmin_4 >= float(image_width_4);
        }
        if(_S804)
        {
            _S804 = true;
        }
        else
        {
            _S804 = ymax_4 <= 0.0f;
        }
        if(_S804)
        {
            _S804 = true;
        }
        else
        {
            _S804 = ymin_4 >= float(image_height_4);
        }
        if(_S804)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_4 = make_int4 (int(xmin_4), int(ymin_4), int(xmax_4), int(ymax_4));
        *depth_4 = 0.5f * (F32_log((dot_0(mean_c_4, mean_c_4) + 9.99999997475242708e-07f)));
        *conic_4 = exp_0(- scale_6);
        float3  _S818 = mean_4 - - mul_0(transpose_0(R_8), t_7);
        float3  _S819 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
        *rgb_4 = _S819;
        float _S820 = _S818.x;
        float _S821 = _S818.y;
        float _S822 = _S818.z;
        float norm_4 = (F32_sqrt((_S820 * _S820 + _S821 * _S821 + _S822 * _S822)));
        float x_24 = _S820 / norm_4;
        float y_7 = _S821 / norm_4;
        float z_4 = _S822 / norm_4;
        float3  _S823 = _S819 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_4)[int(1)] + make_float3 (z_4) * (*sh_coeffs_4)[int(2)] - make_float3 (x_24) * (*sh_coeffs_4)[int(3)]);
        *rgb_4 = _S823;
        float z2_12 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_24 * x_24 - y_7 * y_7;
        float fS1_4 = 2.0f * x_24 * y_7;
        float3  _S824 = _S823 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_7) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_12 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_24) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
        *rgb_4 = _S824;
        float fTmp0C_4 = -2.28522896766662598f * z2_12 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        *rgb_4 = max_0(_S824 + (make_float3 (-0.59004360437393188f * (x_24 * fS1_4 + y_7 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_7) * (*sh_coeffs_4)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_12 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_24) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_24 * fC1_4 - y_7 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_5, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_12, float fy_12, float cx_9, float cy_9, FixedArray<float, 10>  * dist_coeffs_20, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_13, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_9, mean_5) + t_8;
    float3  _S825 = exp_0(scale_7);
    float x_25 = quat_8.y;
    float x2_8 = x_25 * x_25;
    float y2_8 = quat_8.z * quat_8.z;
    float z2_13 = quat_8.w * quat_8.w;
    float xy_8 = quat_8.y * quat_8.z;
    float xz_8 = quat_8.y * quat_8.w;
    float yz_8 = quat_8.z * quat_8.w;
    float wx_8 = quat_8.x * quat_8.y;
    float wy_8 = quat_8.x * quat_8.z;
    float wz_8 = quat_8.x * quat_8.w;
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_13), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_13), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))), makeMatrix<float, 3, 3> (_S825.x, 0.0f, 0.0f, 0.0f, _S825.y, 0.0f, 0.0f, 0.0f, _S825.z));
    Matrix<float, 3, 3>  _S826 = transpose_0(R_9);
    float _S827 = float(image_width_5);
    float _S828 = float(image_height_5);
    float _S829 = 0.30000001192092896f * (0.5f * _S827 / fx_12);
    float _S830 = 0.30000001192092896f * (0.5f * _S828 / fy_12);
    float rz_2 = 1.0f / mean_c_5.z;
    float rz2_2 = rz_2 * rz_2;
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_12 * rz_2, 0.0f, - fx_12 * (mean_c_5.z * (F32_min(((_S827 - cx_9) / fx_12 + _S829), ((F32_max((- (cx_9 / fx_12 + _S829)), (mean_c_5.x * rz_2))))))) * rz2_2, 0.0f, fy_12 * rz_2, - fy_12 * (mean_c_5.z * (F32_min(((_S828 - cy_9) / fy_12 + _S830), ((F32_max((- (cy_9 / fy_12 + _S830)), (mean_c_5.y * rz_2))))))) * rz2_2);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_8, mul_4(mul_4(R_9, mul_4(M_5, transpose_0(M_5))), _S826)), transpose_1(J_8));
    *mean2d_13 = make_float2 (fx_12 * mean_c_5.x * rz_2 + cx_9, fy_12 * mean_c_5.y * rz_2 + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S831 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S831;
    float _S832 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S832;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S831 * _S832 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S833 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_4, - covar2d_5.rows[int(0)].y * invdet_4, - covar2d_5.rows[int(1)].x * invdet_4, covar2d_5.rows[int(0)].x * invdet_4);
    *opacity_5 = 1.0f / (1.0f + (F32_exp((- in_opacity_5))));
    if(antialiased_5)
    {
        *opacity_5 = *opacity_5 * compensation_6;
    }
    float extend_5 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_5 / 0.00392156885936856f)))))))));
    float radius_x_5 = extend_5 * (F32_sqrt((covar2d_5[int(0)].x)));
    float radius_y_5 = extend_5 * (F32_sqrt((covar2d_5[int(1)].y)));
    *aabb_xyxy_5 = make_int4 (int((F32_floor(((*mean2d_13).x - radius_x_5)))), int((F32_floor(((*mean2d_13).y - radius_y_5)))), int((F32_ceil(((*mean2d_13).x + radius_x_5)))), int((F32_ceil(((*mean2d_13).y + radius_y_5)))));
    *depth_5 = 0.5f * (F32_log((dot_0(mean_c_5, mean_c_5) + 9.99999997475242708e-07f)));
    *conic_5 = make_float3 (_S833.rows[int(0)].x, _S833.rows[int(0)].y, _S833.rows[int(1)].y);
    float3  _S834 = mean_5 - - mul_0(_S826, t_8);
    float3  _S835 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S835;
    float _S836 = _S834.x;
    float _S837 = _S834.y;
    float _S838 = _S834.z;
    float norm_5 = (F32_sqrt((_S836 * _S836 + _S837 * _S837 + _S838 * _S838)));
    float x_26 = _S836 / norm_5;
    float y_8 = _S837 / norm_5;
    float z_5 = _S838 / norm_5;
    float3  _S839 = _S835 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_5)[int(1)] + make_float3 (z_5) * (*sh_coeffs_5)[int(2)] - make_float3 (x_26) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S839;
    float z2_14 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_26 * x_26 - y_8 * y_8;
    float fS1_5 = 2.0f * x_26 * y_8;
    float3  _S840 = _S839 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_8) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_14 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_26) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S840;
    float fTmp0C_5 = -2.28522896766662598f * z2_14 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    *rgb_5 = max_0(_S840 + (make_float3 (-0.59004360437393188f * (x_26 * fS1_5 + y_8 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_8) * (*sh_coeffs_5)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_14 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_26) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_26 * fC1_5 - y_8 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_6, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_13, float fy_13, float cx_10, float cy_10, FixedArray<float, 10>  * dist_coeffs_21, uint image_width_6, uint image_height_6, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * mean2d_14, float * depth_6, float3  * conic_6, float * opacity_6, float3  * rgb_6)
{
    float3  mean_c_6 = mul_0(R_10, mean_6) + t_9;
    float3  _S841 = exp_0(scale_8);
    float x_27 = quat_9.y;
    float x2_9 = x_27 * x_27;
    float y2_9 = quat_9.z * quat_9.z;
    float z2_15 = quat_9.w * quat_9.w;
    float xy_9 = quat_9.y * quat_9.z;
    float xz_9 = quat_9.y * quat_9.w;
    float yz_9 = quat_9.z * quat_9.w;
    float wx_9 = quat_9.x * quat_9.y;
    float wy_9 = quat_9.x * quat_9.z;
    float wz_9 = quat_9.x * quat_9.w;
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_15), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_15), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S841.x, 0.0f, 0.0f, 0.0f, _S841.y, 0.0f, 0.0f, 0.0f, _S841.z));
    Matrix<float, 3, 3>  _S842 = transpose_0(R_10);
    Matrix<float, 2, 2>  covar2d_6;
    bool _S843 = fisheye_proj_3dgs_1(mean_c_6, mul_4(mul_4(R_10, mul_4(M_6, transpose_0(M_6))), _S842), make_float4 (fx_13, fy_13, cx_10, cy_10), dist_coeffs_21, &covar2d_6, mean2d_14);
    float det_orig_7 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
    float _S844 = *&(((&covar2d_6)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(0)))->x) = _S844;
    float _S845 = *&(((&covar2d_6)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(1)))->y) = _S845;
    float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / (_S844 * _S845 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
    Matrix<float, 2, 2>  _S846 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_5, - covar2d_6.rows[int(0)].y * invdet_5, - covar2d_6.rows[int(1)].x * invdet_5, covar2d_6.rows[int(0)].x * invdet_5);
    *opacity_6 = 1.0f / (1.0f + (F32_exp((- in_opacity_6))));
    if(antialiased_6)
    {
        *opacity_6 = *opacity_6 * compensation_7;
    }
    float extend_6 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_6 / 0.00392156885936856f)))))))));
    float radius_x_6 = extend_6 * (F32_sqrt((covar2d_6[int(0)].x)));
    float radius_y_6 = extend_6 * (F32_sqrt((covar2d_6[int(1)].y)));
    *aabb_xyxy_6 = make_int4 (int((F32_floor(((*mean2d_14).x - radius_x_6)))), int((F32_floor(((*mean2d_14).y - radius_y_6)))), int((F32_ceil(((*mean2d_14).x + radius_x_6)))), int((F32_ceil(((*mean2d_14).y + radius_y_6)))));
    *depth_6 = 0.5f * (F32_log((dot_0(mean_c_6, mean_c_6) + 9.99999997475242708e-07f)));
    *conic_6 = make_float3 (_S846.rows[int(0)].x, _S846.rows[int(0)].y, _S846.rows[int(1)].y);
    float3  _S847 = mean_6 - - mul_0(_S842, t_9);
    float3  _S848 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)];
    *rgb_6 = _S848;
    float _S849 = _S847.x;
    float _S850 = _S847.y;
    float _S851 = _S847.z;
    float norm_6 = (F32_sqrt((_S849 * _S849 + _S850 * _S850 + _S851 * _S851)));
    float x_28 = _S849 / norm_6;
    float y_9 = _S850 / norm_6;
    float z_6 = _S851 / norm_6;
    float3  _S852 = _S848 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_6)[int(1)] + make_float3 (z_6) * (*sh_coeffs_6)[int(2)] - make_float3 (x_28) * (*sh_coeffs_6)[int(3)]);
    *rgb_6 = _S852;
    float z2_16 = z_6 * z_6;
    float fTmp0B_6 = -1.09254848957061768f * z_6;
    float fC1_6 = x_28 * x_28 - y_9 * y_9;
    float fS1_6 = 2.0f * x_28 * y_9;
    float3  _S853 = _S852 + (make_float3 (0.54627424478530884f * fS1_6) * (*sh_coeffs_6)[int(4)] + make_float3 (fTmp0B_6 * y_9) * (*sh_coeffs_6)[int(5)] + make_float3 (0.94617468118667603f * z2_16 - 0.31539157032966614f) * (*sh_coeffs_6)[int(6)] + make_float3 (fTmp0B_6 * x_28) * (*sh_coeffs_6)[int(7)] + make_float3 (0.54627424478530884f * fC1_6) * (*sh_coeffs_6)[int(8)]);
    *rgb_6 = _S853;
    float fTmp0C_6 = -2.28522896766662598f * z2_16 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_6;
    *rgb_6 = max_0(_S853 + (make_float3 (-0.59004360437393188f * (x_28 * fS1_6 + y_9 * fC1_6)) * (*sh_coeffs_6)[int(9)] + make_float3 (fTmp1B_6 * fS1_6) * (*sh_coeffs_6)[int(10)] + make_float3 (fTmp0C_6 * y_9) * (*sh_coeffs_6)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_16 - 1.11952900886535645f)) * (*sh_coeffs_6)[int(12)] + make_float3 (fTmp0C_6 * x_28) * (*sh_coeffs_6)[int(13)] + make_float3 (fTmp1B_6 * fC1_6) * (*sh_coeffs_6)[int(14)] + make_float3 (-0.59004360437393188f * (x_28 * fC1_6 - y_9 * fS1_6)) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_7, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_14, float fy_14, float cx_11, float cy_11, FixedArray<float, 10>  * dist_coeffs_22, uint image_width_7, uint image_height_7, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * mean2d_15, float * depth_7, float3  * conic_7, float * opacity_7, float3  * rgb_7)
{
    float3  mean_c_7 = mul_0(R_11, mean_7) + t_10;
    float3  _S854 = exp_0(scale_9);
    float x_29 = quat_10.y;
    float x2_10 = x_29 * x_29;
    float y2_10 = quat_10.z * quat_10.z;
    float z2_17 = quat_10.w * quat_10.w;
    float xy_10 = quat_10.y * quat_10.z;
    float xz_10 = quat_10.y * quat_10.w;
    float yz_10 = quat_10.z * quat_10.w;
    float wx_10 = quat_10.x * quat_10.y;
    float wy_10 = quat_10.x * quat_10.z;
    float wz_10 = quat_10.x * quat_10.w;
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_17), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_17), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S854.x, 0.0f, 0.0f, 0.0f, _S854.y, 0.0f, 0.0f, 0.0f, _S854.z));
    Matrix<float, 3, 3>  _S855 = transpose_0(R_11);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (fx_14, 0.0f, 0.0f, 0.0f, fy_14, 0.0f);
    Matrix<float, 2, 2>  covar2d_7 = mul_6(mul_5(J_9, mul_4(mul_4(R_11, mul_4(M_7, transpose_0(M_7))), _S855)), transpose_1(J_9));
    *mean2d_15 = make_float2 (fx_14 * mean_c_7.x + cx_11, fy_14 * mean_c_7.y + cy_11);
    float det_orig_8 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
    float _S856 = *&(((&covar2d_7)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(0)))->x) = _S856;
    float _S857 = *&(((&covar2d_7)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(1)))->y) = _S857;
    float compensation_8 = (F32_sqrt(((F32_max((0.0f), (det_orig_8 / (_S856 * _S857 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
    Matrix<float, 2, 2>  _S858 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_6, - covar2d_7.rows[int(0)].y * invdet_6, - covar2d_7.rows[int(1)].x * invdet_6, covar2d_7.rows[int(0)].x * invdet_6);
    *opacity_7 = 1.0f / (1.0f + (F32_exp((- in_opacity_7))));
    if(antialiased_7)
    {
        *opacity_7 = *opacity_7 * compensation_8;
    }
    float extend_7 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_7 / 0.00392156885936856f)))))))));
    float radius_x_7 = extend_7 * (F32_sqrt((covar2d_7[int(0)].x)));
    float radius_y_7 = extend_7 * (F32_sqrt((covar2d_7[int(1)].y)));
    *aabb_xyxy_7 = make_int4 (int((F32_floor(((*mean2d_15).x - radius_x_7)))), int((F32_floor(((*mean2d_15).y - radius_y_7)))), int((F32_ceil(((*mean2d_15).x + radius_x_7)))), int((F32_ceil(((*mean2d_15).y + radius_y_7)))));
    *depth_7 = 0.5f * (F32_log((dot_0(mean_c_7, mean_c_7) + 9.99999997475242708e-07f)));
    *conic_7 = make_float3 (_S858.rows[int(0)].x, _S858.rows[int(0)].y, _S858.rows[int(1)].y);
    float3  _S859 = mean_7 - - mul_0(_S855, t_10);
    float3  _S860 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)];
    *rgb_7 = _S860;
    float _S861 = _S859.x;
    float _S862 = _S859.y;
    float _S863 = _S859.z;
    float norm_7 = (F32_sqrt((_S861 * _S861 + _S862 * _S862 + _S863 * _S863)));
    float x_30 = _S861 / norm_7;
    float y_10 = _S862 / norm_7;
    float z_7 = _S863 / norm_7;
    float3  _S864 = _S860 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * (*sh_coeffs_7)[int(1)] + make_float3 (z_7) * (*sh_coeffs_7)[int(2)] - make_float3 (x_30) * (*sh_coeffs_7)[int(3)]);
    *rgb_7 = _S864;
    float z2_18 = z_7 * z_7;
    float fTmp0B_7 = -1.09254848957061768f * z_7;
    float fC1_7 = x_30 * x_30 - y_10 * y_10;
    float fS1_7 = 2.0f * x_30 * y_10;
    float3  _S865 = _S864 + (make_float3 (0.54627424478530884f * fS1_7) * (*sh_coeffs_7)[int(4)] + make_float3 (fTmp0B_7 * y_10) * (*sh_coeffs_7)[int(5)] + make_float3 (0.94617468118667603f * z2_18 - 0.31539157032966614f) * (*sh_coeffs_7)[int(6)] + make_float3 (fTmp0B_7 * x_30) * (*sh_coeffs_7)[int(7)] + make_float3 (0.54627424478530884f * fC1_7) * (*sh_coeffs_7)[int(8)]);
    *rgb_7 = _S865;
    float fTmp0C_7 = -2.28522896766662598f * z2_18 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_7;
    *rgb_7 = max_0(_S865 + (make_float3 (-0.59004360437393188f * (x_30 * fS1_7 + y_10 * fC1_7)) * (*sh_coeffs_7)[int(9)] + make_float3 (fTmp1B_7 * fS1_7) * (*sh_coeffs_7)[int(10)] + make_float3 (fTmp0C_7 * y_10) * (*sh_coeffs_7)[int(11)] + make_float3 (z_7 * (1.86588168144226074f * z2_18 - 1.11952900886535645f)) * (*sh_coeffs_7)[int(12)] + make_float3 (fTmp0C_7 * x_30) * (*sh_coeffs_7)[int(13)] + make_float3 (fTmp1B_7 * fC1_7) * (*sh_coeffs_7)[int(14)] + make_float3 (-0.59004360437393188f * (x_30 * fC1_7 - y_10 * fS1_7)) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_persp_differentiable(bool antialiased_8, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_15, float fy_15, float cx_12, float cy_12, FixedArray<float, 10>  * dist_coeffs_23, uint image_width_8, uint image_height_8, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * mean2d_16, float * depth_8, float3  * conic_8, float * opacity_8, float3  * rgb_8)
{
    float3  mean_c_8 = mul_0(R_12, mean_8) + t_11;
    float4  intrins_9 = make_float4 (fx_15, fy_15, cx_12, cy_12);
    float3  _S866 = exp_0(scale_10);
    float x_31 = quat_11.y;
    float x2_11 = x_31 * x_31;
    float y2_11 = quat_11.z * quat_11.z;
    float z2_19 = quat_11.w * quat_11.w;
    float xy_11 = quat_11.y * quat_11.z;
    float xz_11 = quat_11.y * quat_11.w;
    float yz_11 = quat_11.z * quat_11.w;
    float wx_11 = quat_11.x * quat_11.y;
    float wy_11 = quat_11.x * quat_11.z;
    float wz_11 = quat_11.x * quat_11.w;
    Matrix<float, 3, 3>  _S867 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_19), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_19), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))));
    SigmaPoints_0 ret_2;
    (&ret_2)->p_0[int(0)] = mean_8;
    (&ret_2)->w_mean_0[int(0)] = 0.0f;
    (&ret_2)->w_cov_0[int(0)] = 2.0f;
    float _S868 = (F32_sqrt((3.0f)));
    float3  delta_6 = make_float3 (_S868 * _S866.x) * _S867.rows[0U];
    float3  _S869 = mean_8 + delta_6;
    float3  _S870 = mean_8 - delta_6;
    float3  delta_7 = make_float3 (_S868 * _S866.y) * _S867.rows[1U];
    float3  _S871 = mean_8 + delta_7;
    float3  _S872 = mean_8 - delta_7;
    float3  delta_8 = make_float3 (_S868 * _S866.z) * _S867.rows[2U];
    float3  _S873 = mean_8 + delta_8;
    float3  _S874 = mean_8 - delta_8;
    (&ret_2)->w_mean_0[1U] = 0.1666666716337204f;
    (&ret_2)->w_cov_0[1U] = 0.1666666716337204f;
    (&ret_2)->w_mean_0[2U] = 0.1666666716337204f;
    (&ret_2)->w_cov_0[2U] = 0.1666666716337204f;
    (&ret_2)->w_mean_0[3U] = 0.1666666716337204f;
    (&ret_2)->w_cov_0[3U] = 0.1666666716337204f;
    (&ret_2)->w_mean_0[4U] = 0.1666666716337204f;
    (&ret_2)->w_cov_0[4U] = 0.1666666716337204f;
    (&ret_2)->w_mean_0[5U] = 0.1666666716337204f;
    (&ret_2)->w_cov_0[5U] = 0.1666666716337204f;
    (&ret_2)->w_mean_0[6U] = 0.1666666716337204f;
    (&ret_2)->w_cov_0[6U] = 0.1666666716337204f;
    (&ret_2)->p_0[0U] = mul_0(R_12, (&ret_2)->p_0[0U]) + t_11;
    (&ret_2)->p_0[1U] = mul_0(R_12, _S869) + t_11;
    (&ret_2)->p_0[2U] = mul_0(R_12, _S871) + t_11;
    (&ret_2)->p_0[3U] = mul_0(R_12, _S873) + t_11;
    (&ret_2)->p_0[4U] = mul_0(R_12, _S870) + t_11;
    (&ret_2)->p_0[5U] = mul_0(R_12, _S872) + t_11;
    (&ret_2)->p_0[6U] = mul_0(R_12, _S874) + t_11;
    SigmaPoints_0 _S875 = ret_2;
    Matrix<float, 2, 2>  covar2d_8;
    bool _S876 = persp_proj_3dgs_ut_1(&_S875, intrins_9, dist_coeffs_23, image_width_8, image_height_8, &covar2d_8, mean2d_16);
    float det_orig_9 = *&(((&covar2d_8)->rows + (int(0)))->x) * *&(((&covar2d_8)->rows + (int(1)))->y) - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x);
    float _S877 = *&(((&covar2d_8)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(0)))->x) = _S877;
    float _S878 = *&(((&covar2d_8)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(1)))->y) = _S878;
    float compensation_9 = (F32_sqrt(((F32_max((0.0f), (det_orig_9 / (_S877 * _S878 - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x))))))));
    *opacity_8 = 1.0f / (1.0f + (F32_exp((- in_opacity_8))));
    if(antialiased_8)
    {
        *opacity_8 = *opacity_8 * compensation_9;
    }
    float extend_8 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_8 / 0.00392156885936856f)))))))));
    float radius_x_8 = extend_8 * (F32_sqrt((covar2d_8[int(0)].x)));
    float radius_y_8 = extend_8 * (F32_sqrt((covar2d_8[int(1)].y)));
    *aabb_xyxy_8 = make_int4 (int((F32_floor(((*mean2d_16).x - radius_x_8)))), int((F32_floor(((*mean2d_16).y - radius_y_8)))), int((F32_ceil(((*mean2d_16).x + radius_x_8)))), int((F32_ceil(((*mean2d_16).y + radius_y_8)))));
    *depth_8 = 0.5f * (F32_log((dot_0(mean_c_8, mean_c_8) + 9.99999997475242708e-07f)));
    *conic_8 = exp_0(- scale_10);
    float3  _S879 = mean_8 - - mul_0(transpose_0(R_12), t_11);
    float3  _S880 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)];
    *rgb_8 = _S880;
    float _S881 = _S879.x;
    float _S882 = _S879.y;
    float _S883 = _S879.z;
    float norm_8 = (F32_sqrt((_S881 * _S881 + _S882 * _S882 + _S883 * _S883)));
    float x_32 = _S881 / norm_8;
    float y_11 = _S882 / norm_8;
    float z_8 = _S883 / norm_8;
    float3  _S884 = _S880 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_8)[int(1)] + make_float3 (z_8) * (*sh_coeffs_8)[int(2)] - make_float3 (x_32) * (*sh_coeffs_8)[int(3)]);
    *rgb_8 = _S884;
    float z2_20 = z_8 * z_8;
    float fTmp0B_8 = -1.09254848957061768f * z_8;
    float fC1_8 = x_32 * x_32 - y_11 * y_11;
    float fS1_8 = 2.0f * x_32 * y_11;
    float3  _S885 = _S884 + (make_float3 (0.54627424478530884f * fS1_8) * (*sh_coeffs_8)[int(4)] + make_float3 (fTmp0B_8 * y_11) * (*sh_coeffs_8)[int(5)] + make_float3 (0.94617468118667603f * z2_20 - 0.31539157032966614f) * (*sh_coeffs_8)[int(6)] + make_float3 (fTmp0B_8 * x_32) * (*sh_coeffs_8)[int(7)] + make_float3 (0.54627424478530884f * fC1_8) * (*sh_coeffs_8)[int(8)]);
    *rgb_8 = _S885;
    float fTmp0C_8 = -2.28522896766662598f * z2_20 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_8;
    *rgb_8 = max_0(_S885 + (make_float3 (-0.59004360437393188f * (x_32 * fS1_8 + y_11 * fC1_8)) * (*sh_coeffs_8)[int(9)] + make_float3 (fTmp1B_8 * fS1_8) * (*sh_coeffs_8)[int(10)] + make_float3 (fTmp0C_8 * y_11) * (*sh_coeffs_8)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_20 - 1.11952900886535645f)) * (*sh_coeffs_8)[int(12)] + make_float3 (fTmp0C_8 * x_32) * (*sh_coeffs_8)[int(13)] + make_float3 (fTmp1B_8 * fC1_8) * (*sh_coeffs_8)[int(14)] + make_float3 (-0.59004360437393188f * (x_32 * fC1_8 - y_11 * fS1_8)) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_fisheye_differentiable(bool antialiased_9, float3  mean_9, float4  quat_12, float3  scale_11, float in_opacity_9, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_16, float fy_16, float cx_13, float cy_13, FixedArray<float, 10>  * dist_coeffs_24, uint image_width_9, uint image_height_9, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * mean2d_17, float * depth_9, float3  * conic_9, float * opacity_9, float3  * rgb_9)
{
    float3  mean_c_9 = mul_0(R_13, mean_9) + t_12;
    float4  intrins_10 = make_float4 (fx_16, fy_16, cx_13, cy_13);
    float3  _S886 = exp_0(scale_11);
    float x_33 = quat_12.y;
    float x2_12 = x_33 * x_33;
    float y2_12 = quat_12.z * quat_12.z;
    float z2_21 = quat_12.w * quat_12.w;
    float xy_12 = quat_12.y * quat_12.z;
    float xz_12 = quat_12.y * quat_12.w;
    float yz_12 = quat_12.z * quat_12.w;
    float wx_12 = quat_12.x * quat_12.y;
    float wy_12 = quat_12.x * quat_12.z;
    float wz_12 = quat_12.x * quat_12.w;
    Matrix<float, 3, 3>  _S887 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_21), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_21), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))));
    SigmaPoints_0 ret_3;
    (&ret_3)->p_0[int(0)] = mean_9;
    (&ret_3)->w_mean_0[int(0)] = 0.0f;
    (&ret_3)->w_cov_0[int(0)] = 2.0f;
    float _S888 = (F32_sqrt((3.0f)));
    float3  delta_9 = make_float3 (_S888 * _S886.x) * _S887.rows[0U];
    float3  _S889 = mean_9 + delta_9;
    float3  _S890 = mean_9 - delta_9;
    float3  delta_10 = make_float3 (_S888 * _S886.y) * _S887.rows[1U];
    float3  _S891 = mean_9 + delta_10;
    float3  _S892 = mean_9 - delta_10;
    float3  delta_11 = make_float3 (_S888 * _S886.z) * _S887.rows[2U];
    float3  _S893 = mean_9 + delta_11;
    float3  _S894 = mean_9 - delta_11;
    (&ret_3)->w_mean_0[1U] = 0.1666666716337204f;
    (&ret_3)->w_cov_0[1U] = 0.1666666716337204f;
    (&ret_3)->w_mean_0[2U] = 0.1666666716337204f;
    (&ret_3)->w_cov_0[2U] = 0.1666666716337204f;
    (&ret_3)->w_mean_0[3U] = 0.1666666716337204f;
    (&ret_3)->w_cov_0[3U] = 0.1666666716337204f;
    (&ret_3)->w_mean_0[4U] = 0.1666666716337204f;
    (&ret_3)->w_cov_0[4U] = 0.1666666716337204f;
    (&ret_3)->w_mean_0[5U] = 0.1666666716337204f;
    (&ret_3)->w_cov_0[5U] = 0.1666666716337204f;
    (&ret_3)->w_mean_0[6U] = 0.1666666716337204f;
    (&ret_3)->w_cov_0[6U] = 0.1666666716337204f;
    (&ret_3)->p_0[0U] = mul_0(R_13, (&ret_3)->p_0[0U]) + t_12;
    (&ret_3)->p_0[1U] = mul_0(R_13, _S889) + t_12;
    (&ret_3)->p_0[2U] = mul_0(R_13, _S891) + t_12;
    (&ret_3)->p_0[3U] = mul_0(R_13, _S893) + t_12;
    (&ret_3)->p_0[4U] = mul_0(R_13, _S890) + t_12;
    (&ret_3)->p_0[5U] = mul_0(R_13, _S892) + t_12;
    (&ret_3)->p_0[6U] = mul_0(R_13, _S894) + t_12;
    SigmaPoints_0 _S895 = ret_3;
    Matrix<float, 2, 2>  covar2d_9;
    bool _S896 = fisheye_proj_3dgs_ut_1(&_S895, intrins_10, dist_coeffs_24, &covar2d_9, mean2d_17);
    float det_orig_10 = *&(((&covar2d_9)->rows + (int(0)))->x) * *&(((&covar2d_9)->rows + (int(1)))->y) - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x);
    float _S897 = *&(((&covar2d_9)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(0)))->x) = _S897;
    float _S898 = *&(((&covar2d_9)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(1)))->y) = _S898;
    float compensation_10 = (F32_sqrt(((F32_max((0.0f), (det_orig_10 / (_S897 * _S898 - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x))))))));
    *opacity_9 = 1.0f / (1.0f + (F32_exp((- in_opacity_9))));
    if(antialiased_9)
    {
        *opacity_9 = *opacity_9 * compensation_10;
    }
    float extend_9 = (F32_min((3.32999992370605469f), ((F32_sqrt((2.0f * (F32_log((*opacity_9 / 0.00392156885936856f)))))))));
    float radius_x_9 = extend_9 * (F32_sqrt((covar2d_9[int(0)].x)));
    float radius_y_9 = extend_9 * (F32_sqrt((covar2d_9[int(1)].y)));
    *aabb_xyxy_9 = make_int4 (int((F32_floor(((*mean2d_17).x - radius_x_9)))), int((F32_floor(((*mean2d_17).y - radius_y_9)))), int((F32_ceil(((*mean2d_17).x + radius_x_9)))), int((F32_ceil(((*mean2d_17).y + radius_y_9)))));
    *depth_9 = 0.5f * (F32_log((dot_0(mean_c_9, mean_c_9) + 9.99999997475242708e-07f)));
    *conic_9 = exp_0(- scale_11);
    float3  _S899 = mean_9 - - mul_0(transpose_0(R_13), t_12);
    float3  _S900 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
    *rgb_9 = _S900;
    float _S901 = _S899.x;
    float _S902 = _S899.y;
    float _S903 = _S899.z;
    float norm_9 = (F32_sqrt((_S901 * _S901 + _S902 * _S902 + _S903 * _S903)));
    float x_34 = _S901 / norm_9;
    float y_12 = _S902 / norm_9;
    float z_9 = _S903 / norm_9;
    float3  _S904 = _S900 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * (*sh_coeffs_9)[int(1)] + make_float3 (z_9) * (*sh_coeffs_9)[int(2)] - make_float3 (x_34) * (*sh_coeffs_9)[int(3)]);
    *rgb_9 = _S904;
    float z2_22 = z_9 * z_9;
    float fTmp0B_9 = -1.09254848957061768f * z_9;
    float fC1_9 = x_34 * x_34 - y_12 * y_12;
    float fS1_9 = 2.0f * x_34 * y_12;
    float3  _S905 = _S904 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_12) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_22 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_34) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
    *rgb_9 = _S905;
    float fTmp0C_9 = -2.28522896766662598f * z2_22 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_9;
    *rgb_9 = max_0(_S905 + (make_float3 (-0.59004360437393188f * (x_34 * fS1_9 + y_12 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_12) * (*sh_coeffs_9)[int(11)] + make_float3 (z_9 * (1.86588168144226074f * z2_22 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_34) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_34 * fC1_9 - y_12 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S906, float3  _S907)
{
    return mul_0(_S906, _S907);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S908)
{
    return exp_0(_S908);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S909, Matrix<float, 3, 3>  _S910)
{
    return mul_4(_S909, _S910);
}

inline __device__ float s_primal_ctx_max_0(float _S911, float _S912)
{
    return (F32_max((_S911), (_S912)));
}

inline __device__ float s_primal_ctx_min_0(float _S913, float _S914)
{
    return (F32_min((_S913), (_S914)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S915, Matrix<float, 3, 3>  _S916)
{
    return mul_5(_S915, _S916);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_4(Matrix<float, 2, 3>  _S917, Matrix<float, 3, 2>  _S918)
{
    return mul_6(_S917, _S918);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S919)
{
    return (F32_sqrt((_S919)));
}

inline __device__ float s_primal_ctx_exp_1(float _S920)
{
    return (F32_exp((_S920)));
}

inline __device__ float s_primal_ctx_log_0(float _S921)
{
    return (F32_log((_S921)));
}

inline __device__ float s_primal_ctx_dot_0(float3  _S922, float3  _S923)
{
    return dot_0(_S922, _S923);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S924, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S925, float3  _S926)
{
    _d_max_vector_0(_S924, _S925, _S926);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S927, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S928, float3  _S929)
{
    _d_mul_0(_S927, _S928, _S929);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S930, float _S931)
{
    _d_log_0(_S930, _S931);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S932, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S933, float _S934)
{
    _d_dot_0(_S932, _S933, _S934);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S935, DiffPair_float_0 * _S936, float _S937)
{
    _d_min_0(_S935, _S936, _S937);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S938, float _S939)
{
    _d_exp_0(_S938, _S939);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S940, DiffPair_float_0 * _S941, float _S942)
{
    _d_max_0(_S940, _S941, _S942);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S943, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S944, Matrix<float, 2, 2>  _S945)
{
    mul_3(_S943, _S944, _S945);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S946, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S947, Matrix<float, 2, 3>  _S948)
{
    mul_2(_S946, _S947, _S948);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S949, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S950, Matrix<float, 3, 3>  _S951)
{
    mul_1(_S949, _S950, _S951);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S952, float3  _S953)
{
    _d_exp_vector_0(_S952, _S953);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_17, float fy_17, float cx_14, float cy_14, FixedArray<float, 10>  * dist_coeffs_25, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_10 = s_primal_ctx_mul_1(R_14, mean_10) + t_13;
    float3  _S954 = s_primal_ctx_exp_0(scale_12);
    float _S955 = quat_13.y;
    float x2_13 = _S955 * _S955;
    float y2_13 = quat_13.z * quat_13.z;
    float z2_23 = quat_13.w * quat_13.w;
    float xy_13 = quat_13.y * quat_13.z;
    float xz_13 = quat_13.y * quat_13.w;
    float yz_13 = quat_13.z * quat_13.w;
    float wx_13 = quat_13.x * quat_13.y;
    float wy_13 = quat_13.x * quat_13.z;
    float wz_13 = quat_13.x * quat_13.w;
    Matrix<float, 3, 3>  _S956 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_23), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_23), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S954.x, 0.0f, 0.0f, 0.0f, _S954.y, 0.0f, 0.0f, 0.0f, _S954.z);
    Matrix<float, 3, 3>  _S957 = s_primal_ctx_mul_2(_S956, S_0);
    Matrix<float, 3, 3>  _S958 = transpose_0(_S957);
    Matrix<float, 3, 3>  _S959 = s_primal_ctx_mul_2(_S957, _S958);
    Matrix<float, 3, 3>  _S960 = s_primal_ctx_mul_2(R_14, _S959);
    Matrix<float, 3, 3>  _S961 = transpose_0(R_14);
    Matrix<float, 3, 3>  _S962 = s_primal_ctx_mul_2(_S960, _S961);
    float _S963 = float(image_width_10);
    float _S964 = float(image_height_10);
    float _S965 = 0.30000001192092896f * (0.5f * _S963 / fx_17);
    float lim_x_pos_2 = (_S963 - cx_14) / fx_17 + _S965;
    float _S966 = 0.30000001192092896f * (0.5f * _S964 / fy_17);
    float lim_y_pos_2 = (_S964 - cy_14) / fy_17 + _S966;
    float rz_3 = 1.0f / mean_c_10.z;
    float _S967 = mean_c_10.z * mean_c_10.z;
    float rz2_3 = rz_3 * rz_3;
    float _S968 = - (cx_14 / fx_17 + _S965);
    float _S969 = mean_c_10.x * rz_3;
    float _S970 = s_primal_ctx_max_0(_S968, _S969);
    float _S971 = s_primal_ctx_min_0(lim_x_pos_2, _S970);
    float _S972 = - (cy_14 / fy_17 + _S966);
    float _S973 = mean_c_10.y * rz_3;
    float _S974 = s_primal_ctx_max_0(_S972, _S973);
    float _S975 = s_primal_ctx_min_0(lim_y_pos_2, _S974);
    float _S976 = - fx_17;
    float _S977 = _S976 * (mean_c_10.z * _S971);
    float _S978 = - fy_17;
    float _S979 = _S978 * (mean_c_10.z * _S975);
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (fx_17 * rz_3, 0.0f, _S977 * rz2_3, 0.0f, fy_17 * rz_3, _S979 * rz2_3);
    Matrix<float, 2, 3>  _S980 = s_primal_ctx_mul_3(J_10, _S962);
    Matrix<float, 3, 2>  _S981 = transpose_1(J_10);
    Matrix<float, 2, 2>  _S982 = s_primal_ctx_mul_4(_S980, _S981);
    float _S983 = fx_17 * mean_c_10.x;
    float _S984 = fy_17 * mean_c_10.y;
    float _S985 = _S982.rows[int(0)].y * _S982.rows[int(1)].x;
    float det_orig_11 = _S982.rows[int(0)].x * _S982.rows[int(1)].y - _S985;
    float _S986 = _S982.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S987 = _S982;
    *&(((&_S987)->rows + (int(0)))->x) = _S986;
    float _S988 = _S982.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S987)->rows + (int(1)))->y) = _S988;
    Matrix<float, 2, 2>  _S989 = _S987;
    Matrix<float, 2, 2>  _S990 = _S987;
    float det_blur_6 = _S986 * _S988 - _S985;
    float _S991 = det_orig_11 / det_blur_6;
    float _S992 = det_blur_6 * det_blur_6;
    float _S993 = s_primal_ctx_max_0(0.0f, _S991);
    float _S994 = s_primal_ctx_sqrt_0(_S993);
    float invdet_7 = 1.0f / det_blur_6;
    float _S995 = - _S982.rows[int(0)].y;
    float _S996 = - _S982.rows[int(1)].x;
    float _S997 = - in_opacity_10;
    float _S998 = 1.0f + s_primal_ctx_exp_1(_S997);
    float _S999 = 1.0f / _S998;
    float _S1000 = _S998 * _S998;
    float _S1001;
    if(antialiased_10)
    {
        _S1001 = _S999 * _S994;
    }
    else
    {
        _S1001 = _S999;
    }
    float _S1002 = _S1001 / 0.00392156885936856f;
    float _S1003 = 2.0f * s_primal_ctx_log_0(_S1002);
    float _S1004 = s_primal_ctx_sqrt_0(_S1003);
    float _S1005 = _S989.rows[int(0)].x;
    float _S1006 = _S990.rows[int(1)].y;
    float _S1007 = s_primal_ctx_dot_0(mean_c_10, mean_c_10) + 9.99999997475242708e-07f;
    float3  _S1008 = mean_10 - - s_primal_ctx_mul_1(_S961, t_13);
    float _S1009 = _S1008.x;
    float _S1010 = _S1008.y;
    float _S1011 = _S1008.z;
    float _S1012 = _S1009 * _S1009 + _S1010 * _S1010 + _S1011 * _S1011;
    float _S1013 = s_primal_ctx_sqrt_0(_S1012);
    float x_35 = _S1009 / _S1013;
    float3  _S1014 = make_float3 (x_35);
    float _S1015 = _S1013 * _S1013;
    float y_13 = _S1010 / _S1013;
    float z_10 = _S1011 / _S1013;
    float3  _S1016 = make_float3 (z_10);
    float _S1017 = - y_13;
    float3  _S1018 = make_float3 (_S1017);
    float z2_24 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_35 * x_35 - y_13 * y_13;
    float _S1019 = 2.0f * x_35;
    float fS1_10 = _S1019 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_24 - 0.31539157032966614f;
    float3  _S1020 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_35;
    float3  _S1021 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S1022 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S1023 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S1024 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_24 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S1025 = 1.86588168144226074f * z2_24 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S1025;
    float3  _S1026 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_35;
    float3  _S1027 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S1028 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S1029 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S1030 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_35 * fC1_10 - y_13 * fS1_10);
    float3  _S1031 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_35 * fS1_10 + y_13 * fC1_10);
    float3  _S1032 = make_float3 (pSH9_0);
    float3  _S1033 = make_float3 (0.0f);
    float3  _S1034 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1035;
    (&_S1035)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1017) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_35) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S1035)->differential_0 = _S1034;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1036;
    (&_S1036)->primal_0 = _S1033;
    (&_S1036)->differential_0 = _S1034;
    s_bwd_prop_max_0(&_S1035, &_S1036, v_rgb_0);
    float3  _S1037 = _S1031 * _S1035.differential_0;
    float3  _S1038 = (*sh_coeffs_10)[int(15)] * _S1035.differential_0;
    float3  _S1039 = _S1029 * _S1035.differential_0;
    float3  _S1040 = (*sh_coeffs_10)[int(14)] * _S1035.differential_0;
    float3  _S1041 = _S1027 * _S1035.differential_0;
    float3  _S1042 = (*sh_coeffs_10)[int(13)] * _S1035.differential_0;
    float3  _S1043 = _S1026 * _S1035.differential_0;
    float3  _S1044 = (*sh_coeffs_10)[int(12)] * _S1035.differential_0;
    float3  _S1045 = _S1028 * _S1035.differential_0;
    float3  _S1046 = (*sh_coeffs_10)[int(11)] * _S1035.differential_0;
    float3  _S1047 = _S1030 * _S1035.differential_0;
    float3  _S1048 = (*sh_coeffs_10)[int(10)] * _S1035.differential_0;
    float3  _S1049 = _S1032 * _S1035.differential_0;
    float3  _S1050 = (*sh_coeffs_10)[int(9)] * _S1035.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S1050.x + _S1050.y + _S1050.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S1038.x + _S1038.y + _S1038.z);
    float _S1051 = _S1048.x + _S1048.y + _S1048.z;
    float _S1052 = _S1040.x + _S1040.y + _S1040.z;
    float _S1053 = _S1046.x + _S1046.y + _S1046.z;
    float _S1054 = _S1042.x + _S1042.y + _S1042.z;
    float _S1055 = _S1044.x + _S1044.y + _S1044.z;
    float _S1056 = - s_diff_fC2_T_0;
    float3  _S1057 = _S1023 * _S1035.differential_0;
    float3  _S1058 = (*sh_coeffs_10)[int(8)] * _S1035.differential_0;
    float3  _S1059 = _S1021 * _S1035.differential_0;
    float3  _S1060 = (*sh_coeffs_10)[int(7)] * _S1035.differential_0;
    float3  _S1061 = _S1020 * _S1035.differential_0;
    float3  _S1062 = (*sh_coeffs_10)[int(6)] * _S1035.differential_0;
    float3  _S1063 = _S1022 * _S1035.differential_0;
    float3  _S1064 = (*sh_coeffs_10)[int(5)] * _S1035.differential_0;
    float3  _S1065 = _S1024 * _S1035.differential_0;
    float3  _S1066 = (*sh_coeffs_10)[int(4)] * _S1035.differential_0;
    float _S1067 = _S1064.x + _S1064.y + _S1064.z;
    float _S1068 = _S1060.x + _S1060.y + _S1060.z;
    float _S1069 = fTmp1B_10 * _S1051 + x_35 * s_diff_fS2_T_0 + y_13 * _S1056 + 0.54627424478530884f * (_S1066.x + _S1066.y + _S1066.z);
    float _S1070 = fTmp1B_10 * _S1052 + y_13 * s_diff_fS2_T_0 + x_35 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S1058.x + _S1058.y + _S1058.z);
    float _S1071 = y_13 * - _S1070;
    float _S1072 = x_35 * _S1070;
    float _S1073 = z_10 * (1.86588168144226074f * (z_10 * _S1055) + -2.28522896766662598f * (y_13 * _S1053 + x_35 * _S1054) + 0.94617468118667603f * (_S1062.x + _S1062.y + _S1062.z));
    float3  _S1074 = make_float3 (0.48860251903533936f) * _S1035.differential_0;
    float3  _S1075 = - _S1074;
    float3  _S1076 = _S1014 * _S1075;
    float3  _S1077 = (*sh_coeffs_10)[int(3)] * _S1075;
    float3  _S1078 = _S1016 * _S1074;
    float3  _S1079 = (*sh_coeffs_10)[int(2)] * _S1074;
    float3  _S1080 = _S1018 * _S1074;
    float3  _S1081 = (*sh_coeffs_10)[int(1)] * _S1074;
    float _S1082 = (_S1025 * _S1055 + 1.44530570507049561f * (fS1_10 * _S1051 + fC1_10 * _S1052) + -1.09254848957061768f * (y_13 * _S1067 + x_35 * _S1068) + _S1073 + _S1073 + _S1079.x + _S1079.y + _S1079.z) / _S1015;
    float _S1083 = _S1013 * _S1082;
    float _S1084 = (fTmp0C_10 * _S1053 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S1056 + fTmp0B_10 * _S1067 + _S1019 * _S1069 + _S1071 + _S1071 + - (_S1081.x + _S1081.y + _S1081.z)) / _S1015;
    float _S1085 = _S1013 * _S1084;
    float _S1086 = (fTmp0C_10 * _S1054 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S1068 + 2.0f * (y_13 * _S1069) + _S1072 + _S1072 + _S1077.x + _S1077.y + _S1077.z) / _S1015;
    float _S1087 = _S1013 * _S1086;
    float _S1088 = _S1011 * - _S1082 + _S1010 * - _S1084 + _S1009 * - _S1086;
    DiffPair_float_0 _S1089;
    (&_S1089)->primal_0 = _S1012;
    (&_S1089)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1089, _S1088);
    float _S1090 = _S1011 * _S1089.differential_0;
    float _S1091 = _S1010 * _S1089.differential_0;
    float _S1092 = _S1009 * _S1089.differential_0;
    float3  _S1093 = make_float3 (0.282094806432724f) * _S1035.differential_0;
    float3  _S1094 = make_float3 (_S1087 + _S1092 + _S1092, _S1085 + _S1091 + _S1091, _S1083 + _S1090 + _S1090);
    float3  _S1095 = - - _S1094;
    Matrix<float, 3, 3>  _S1096 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1097;
    (&_S1097)->primal_0 = _S961;
    (&_S1097)->differential_0 = _S1096;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1098;
    (&_S1098)->primal_0 = t_13;
    (&_S1098)->differential_0 = _S1034;
    s_bwd_prop_mul_1(&_S1097, &_S1098, _S1095);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1099 = _S1097;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1100 = _S1098;
    float2  _S1101 = make_float2 (0.0f);
    float2  _S1102 = _S1101;
    *&((&_S1102)->y) = v_conic_0.z;
    float2  _S1103 = _S1101;
    *&((&_S1103)->y) = v_conic_0.y;
    *&((&_S1103)->x) = v_conic_0.x;
    float _S1104 = 0.5f * v_depth_0;
    DiffPair_float_0 _S1105;
    (&_S1105)->primal_0 = _S1007;
    (&_S1105)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1105, _S1104);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1106;
    (&_S1106)->primal_0 = mean_c_10;
    (&_S1106)->differential_0 = _S1034;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1107;
    (&_S1107)->primal_0 = mean_c_10;
    (&_S1107)->differential_0 = _S1034;
    s_bwd_prop_dot_0(&_S1106, &_S1107, _S1105.differential_0);
    DiffPair_float_0 _S1108;
    (&_S1108)->primal_0 = _S1006;
    (&_S1108)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1108, 0.0f);
    DiffPair_float_0 _S1109;
    (&_S1109)->primal_0 = _S1005;
    (&_S1109)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1109, 0.0f);
    DiffPair_float_0 _S1110;
    (&_S1110)->primal_0 = 3.32999992370605469f;
    (&_S1110)->differential_0 = 0.0f;
    DiffPair_float_0 _S1111;
    (&_S1111)->primal_0 = _S1004;
    (&_S1111)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1110, &_S1111, 0.0f);
    DiffPair_float_0 _S1112;
    (&_S1112)->primal_0 = _S1003;
    (&_S1112)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1112, _S1111.differential_0);
    float _S1113 = 2.0f * _S1112.differential_0;
    DiffPair_float_0 _S1114;
    (&_S1114)->primal_0 = _S1002;
    (&_S1114)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1114, _S1113);
    float _S1115 = v_opacity_0 + 254.9999847412109375f * _S1114.differential_0;
    Matrix<float, 2, 2>  _S1116 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1117 = _S1116;
    _S1117[int(1)] = _S1102;
    _S1117[int(0)] = _S1103;
    Matrix<float, 2, 2>  _S1118 = _S1117;
    FixedArray<float3 , 16>  _S1119;
    _S1119[int(0)] = _S1034;
    _S1119[int(1)] = _S1034;
    _S1119[int(2)] = _S1034;
    _S1119[int(3)] = _S1034;
    _S1119[int(4)] = _S1034;
    _S1119[int(5)] = _S1034;
    _S1119[int(6)] = _S1034;
    _S1119[int(7)] = _S1034;
    _S1119[int(8)] = _S1034;
    _S1119[int(9)] = _S1034;
    _S1119[int(10)] = _S1034;
    _S1119[int(11)] = _S1034;
    _S1119[int(12)] = _S1034;
    _S1119[int(13)] = _S1034;
    _S1119[int(14)] = _S1034;
    _S1119[int(15)] = _S1034;
    _S1119[int(7)] = _S1059;
    _S1119[int(0)] = _S1093;
    _S1119[int(1)] = _S1080;
    _S1119[int(2)] = _S1078;
    _S1119[int(3)] = _S1076;
    _S1119[int(4)] = _S1065;
    _S1119[int(5)] = _S1063;
    _S1119[int(6)] = _S1061;
    _S1119[int(15)] = _S1037;
    _S1119[int(8)] = _S1057;
    _S1119[int(9)] = _S1049;
    _S1119[int(10)] = _S1047;
    _S1119[int(11)] = _S1045;
    _S1119[int(12)] = _S1043;
    _S1119[int(13)] = _S1041;
    _S1119[int(14)] = _S1039;
    float3  _S1120 = _S1119[int(0)];
    float3  _S1121 = _S1119[int(1)];
    float3  _S1122 = _S1119[int(2)];
    float3  _S1123 = _S1119[int(3)];
    float3  _S1124 = _S1119[int(4)];
    float3  _S1125 = _S1119[int(5)];
    float3  _S1126 = _S1119[int(6)];
    float3  _S1127 = _S1119[int(7)];
    float3  _S1128 = _S1119[int(8)];
    float3  _S1129 = _S1119[int(9)];
    float3  _S1130 = _S1119[int(10)];
    float3  _S1131 = _S1119[int(11)];
    float3  _S1132 = _S1119[int(12)];
    float3  _S1133 = _S1119[int(13)];
    float3  _S1134 = _S1119[int(14)];
    float3  _S1135 = _S1119[int(15)];
    float3  _S1136 = _S1107.differential_0 + _S1106.differential_0;
    float2  _S1137 = make_float2 (0.0f, _S1108.differential_0);
    float2  _S1138 = make_float2 (_S1109.differential_0, 0.0f);
    float _S1139;
    if(antialiased_10)
    {
        float _S1140 = _S999 * _S1115;
        _S1001 = _S994 * _S1115;
        _S1139 = _S1140;
    }
    else
    {
        _S1001 = _S1115;
        _S1139 = 0.0f;
    }
    float _S1141 = - (_S1001 / _S1000);
    DiffPair_float_0 _S1142;
    (&_S1142)->primal_0 = _S997;
    (&_S1142)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1142, _S1141);
    float _S1143 = - _S1142.differential_0;
    float _S1144 = invdet_7 * _S1118.rows[int(1)].y;
    float _S1145 = - (invdet_7 * _S1118.rows[int(1)].x);
    float _S1146 = - (invdet_7 * _S1118.rows[int(0)].y);
    float _S1147 = invdet_7 * _S1118.rows[int(0)].x;
    float _S1148 = - ((_S986 * _S1118.rows[int(1)].y + _S996 * _S1118.rows[int(1)].x + _S995 * _S1118.rows[int(0)].y + _S988 * _S1118.rows[int(0)].x) / _S992);
    DiffPair_float_0 _S1149;
    (&_S1149)->primal_0 = _S993;
    (&_S1149)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1149, _S1139);
    DiffPair_float_0 _S1150;
    (&_S1150)->primal_0 = 0.0f;
    (&_S1150)->differential_0 = 0.0f;
    DiffPair_float_0 _S1151;
    (&_S1151)->primal_0 = _S991;
    (&_S1151)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1150, &_S1151, _S1149.differential_0);
    float _S1152 = _S1151.differential_0 / _S992;
    float s_diff_det_orig_T_0 = det_blur_6 * _S1152;
    float _S1153 = _S1148 + det_orig_11 * - _S1152;
    float _S1154 = - _S1153;
    float _S1155 = _S986 * _S1153;
    float _S1156 = _S988 * _S1153;
    Matrix<float, 2, 2>  _S1157 = _S1116;
    _S1157[int(1)] = _S1137;
    _S1157[int(0)] = _S1138;
    _S987 = _S1157;
    *&(((&_S987)->rows + (int(1)))->y) = 0.0f;
    float _S1158 = _S1147 + _S1155 + _S1157.rows[int(1)].y;
    *&(((&_S987)->rows + (int(0)))->x) = 0.0f;
    float _S1159 = _S1144 + _S1156 + _S1157.rows[int(0)].x;
    float _S1160 = _S1154 + - s_diff_det_orig_T_0;
    float _S1161 = _S1145 + _S982.rows[int(0)].y * _S1160;
    float _S1162 = _S1146 + _S982.rows[int(1)].x * _S1160;
    float _S1163 = _S982.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S1164 = _S1158 + _S982.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S1165 = _S1101;
    *&((&_S1165)->x) = _S1161;
    *&((&_S1165)->y) = _S1164;
    float _S1166 = _S1159 + _S1163;
    float2  _S1167 = _S1101;
    *&((&_S1167)->y) = _S1162;
    *&((&_S1167)->x) = _S1166;
    float _S1168 = _S984 * v_mean2d_0.y;
    float _S1169 = fy_17 * (rz_3 * v_mean2d_0.y);
    float _S1170 = _S983 * v_mean2d_0.x;
    float _S1171 = fx_17 * (rz_3 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S1172 = _S1116;
    _S1172[int(1)] = _S1165;
    _S1172[int(0)] = _S1167;
    Matrix<float, 2, 2>  _S1173 = _S987 + _S1172;
    Matrix<float, 2, 3>  _S1174 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1175;
    (&_S1175)->primal_0 = _S980;
    (&_S1175)->differential_0 = _S1174;
    Matrix<float, 3, 2>  _S1176 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1177;
    (&_S1177)->primal_0 = _S981;
    (&_S1177)->differential_0 = _S1176;
    s_bwd_prop_mul_2(&_S1175, &_S1177, _S1173);
    Matrix<float, 2, 3>  _S1178 = transpose_2(_S1177.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1179;
    (&_S1179)->primal_0 = J_10;
    (&_S1179)->differential_0 = _S1174;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1180;
    (&_S1180)->primal_0 = _S962;
    (&_S1180)->differential_0 = _S1096;
    s_bwd_prop_mul_3(&_S1179, &_S1180, _S1175.differential_0);
    Matrix<float, 2, 3>  _S1181 = _S1178 + _S1179.differential_0;
    float _S1182 = _S979 * _S1181.rows[int(1)].z;
    float s_diff_ty_T_0 = _S978 * (rz2_3 * _S1181.rows[int(1)].z);
    float _S1183 = fy_17 * _S1181.rows[int(1)].y;
    float _S1184 = _S977 * _S1181.rows[int(0)].z;
    float s_diff_tx_T_0 = _S976 * (rz2_3 * _S1181.rows[int(0)].z);
    float _S1185 = fx_17 * _S1181.rows[int(0)].x;
    float _S1186 = mean_c_10.z * s_diff_ty_T_0;
    float _S1187 = _S975 * s_diff_ty_T_0;
    DiffPair_float_0 _S1188;
    (&_S1188)->primal_0 = lim_y_pos_2;
    (&_S1188)->differential_0 = 0.0f;
    DiffPair_float_0 _S1189;
    (&_S1189)->primal_0 = _S974;
    (&_S1189)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1188, &_S1189, _S1186);
    DiffPair_float_0 _S1190;
    (&_S1190)->primal_0 = _S972;
    (&_S1190)->differential_0 = 0.0f;
    DiffPair_float_0 _S1191;
    (&_S1191)->primal_0 = _S973;
    (&_S1191)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1190, &_S1191, _S1189.differential_0);
    float _S1192 = mean_c_10.y * _S1191.differential_0;
    float _S1193 = rz_3 * _S1191.differential_0;
    float _S1194 = mean_c_10.z * s_diff_tx_T_0;
    float _S1195 = _S971 * s_diff_tx_T_0;
    DiffPair_float_0 _S1196;
    (&_S1196)->primal_0 = lim_x_pos_2;
    (&_S1196)->differential_0 = 0.0f;
    DiffPair_float_0 _S1197;
    (&_S1197)->primal_0 = _S970;
    (&_S1197)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1196, &_S1197, _S1194);
    DiffPair_float_0 _S1198;
    (&_S1198)->primal_0 = _S968;
    (&_S1198)->differential_0 = 0.0f;
    DiffPair_float_0 _S1199;
    (&_S1199)->primal_0 = _S969;
    (&_S1199)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1198, &_S1199, _S1197.differential_0);
    float _S1200 = rz_3 * (_S1182 + _S1184);
    float _S1201 = _S1187 + _S1195 + - ((_S1168 + _S1170 + _S1183 + _S1185 + _S1192 + mean_c_10.x * _S1199.differential_0 + _S1200 + _S1200) / _S967);
    float _S1202 = _S1169 + _S1193;
    float _S1203 = _S1171 + rz_3 * _S1199.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1204;
    (&_S1204)->primal_0 = _S960;
    (&_S1204)->differential_0 = _S1096;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1205;
    (&_S1205)->primal_0 = _S961;
    (&_S1205)->differential_0 = _S1096;
    s_bwd_prop_mul_4(&_S1204, &_S1205, _S1180.differential_0);
    Matrix<float, 3, 3>  _S1206 = transpose_0(_S1205.differential_0 + _S1099.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1207;
    (&_S1207)->primal_0 = R_14;
    (&_S1207)->differential_0 = _S1096;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1208;
    (&_S1208)->primal_0 = _S959;
    (&_S1208)->differential_0 = _S1096;
    s_bwd_prop_mul_4(&_S1207, &_S1208, _S1204.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1209;
    (&_S1209)->primal_0 = _S957;
    (&_S1209)->differential_0 = _S1096;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1210;
    (&_S1210)->primal_0 = _S958;
    (&_S1210)->differential_0 = _S1096;
    s_bwd_prop_mul_4(&_S1209, &_S1210, _S1208.differential_0);
    Matrix<float, 3, 3>  _S1211 = _S1209.differential_0 + transpose_0(_S1210.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1212;
    (&_S1212)->primal_0 = _S956;
    (&_S1212)->differential_0 = _S1096;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1213;
    (&_S1213)->primal_0 = S_0;
    (&_S1213)->differential_0 = _S1096;
    s_bwd_prop_mul_4(&_S1212, &_S1213, _S1211);
    Matrix<float, 3, 3>  _S1214 = transpose_0(_S1212.differential_0);
    float _S1215 = 2.0f * - _S1214.rows[int(2)].z;
    float _S1216 = 2.0f * _S1214.rows[int(2)].y;
    float _S1217 = 2.0f * _S1214.rows[int(2)].x;
    float _S1218 = 2.0f * _S1214.rows[int(1)].z;
    float _S1219 = 2.0f * - _S1214.rows[int(1)].y;
    float _S1220 = 2.0f * _S1214.rows[int(1)].x;
    float _S1221 = 2.0f * _S1214.rows[int(0)].z;
    float _S1222 = 2.0f * _S1214.rows[int(0)].y;
    float _S1223 = 2.0f * - _S1214.rows[int(0)].x;
    float _S1224 = - _S1220 + _S1222;
    float _S1225 = _S1217 + - _S1221;
    float _S1226 = - _S1216 + _S1218;
    float _S1227 = _S1216 + _S1218;
    float _S1228 = _S1217 + _S1221;
    float _S1229 = _S1220 + _S1222;
    float _S1230 = quat_13.w * (_S1219 + _S1223);
    float _S1231 = quat_13.z * (_S1215 + _S1223);
    float _S1232 = quat_13.y * (_S1215 + _S1219);
    float _S1233 = quat_13.x * _S1224 + quat_13.z * _S1227 + quat_13.y * _S1228 + _S1230 + _S1230;
    float _S1234 = quat_13.x * _S1225 + quat_13.w * _S1227 + quat_13.y * _S1229 + _S1231 + _S1231;
    float _S1235 = quat_13.x * _S1226 + quat_13.w * _S1228 + quat_13.z * _S1229 + _S1232 + _S1232;
    float _S1236 = quat_13.w * _S1224 + quat_13.z * _S1225 + quat_13.y * _S1226;
    float3  _S1237 = _S1034;
    *&((&_S1237)->z) = _S1213.differential_0.rows[int(2)].z;
    *&((&_S1237)->y) = _S1213.differential_0.rows[int(1)].y;
    *&((&_S1237)->x) = _S1213.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1238;
    (&_S1238)->primal_0 = scale_12;
    (&_S1238)->differential_0 = _S1034;
    s_bwd_prop_exp_1(&_S1238, _S1237);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1239 = _S1238;
    float3  _S1240 = _S1034;
    *&((&_S1240)->z) = _S1201;
    *&((&_S1240)->y) = _S1202;
    *&((&_S1240)->x) = _S1203;
    float3  _S1241 = _S1136 + _S1240;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1242;
    (&_S1242)->primal_0 = R_14;
    (&_S1242)->differential_0 = _S1096;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1243;
    (&_S1243)->primal_0 = mean_10;
    (&_S1243)->differential_0 = _S1034;
    s_bwd_prop_mul_1(&_S1242, &_S1243, _S1241);
    float3  _S1244 = _S1241 + _S1100.differential_0;
    Matrix<float, 3, 3>  _S1245 = _S1206 + _S1207.differential_0 + _S1242.differential_0;
    float4  _S1246 = make_float4 (0.0f);
    *&((&_S1246)->w) = _S1233;
    *&((&_S1246)->z) = _S1234;
    *&((&_S1246)->y) = _S1235;
    *&((&_S1246)->x) = _S1236;
    float4  _S1247 = _S1246;
    float3  _S1248 = _S1243.differential_0 + _S1094;
    *v_mean_0 = _S1248;
    *v_quat_0 = _S1247;
    *v_scale_0 = _S1239.differential_0;
    *v_in_opacity_0 = _S1143;
    (*v_sh_coeffs_0)[int(0)] = _S1120;
    (*v_sh_coeffs_0)[int(1)] = _S1121;
    (*v_sh_coeffs_0)[int(2)] = _S1122;
    (*v_sh_coeffs_0)[int(3)] = _S1123;
    (*v_sh_coeffs_0)[int(4)] = _S1124;
    (*v_sh_coeffs_0)[int(5)] = _S1125;
    (*v_sh_coeffs_0)[int(6)] = _S1126;
    (*v_sh_coeffs_0)[int(7)] = _S1127;
    (*v_sh_coeffs_0)[int(8)] = _S1128;
    (*v_sh_coeffs_0)[int(9)] = _S1129;
    (*v_sh_coeffs_0)[int(10)] = _S1130;
    (*v_sh_coeffs_0)[int(11)] = _S1131;
    (*v_sh_coeffs_0)[int(12)] = _S1132;
    (*v_sh_coeffs_0)[int(13)] = _S1133;
    (*v_sh_coeffs_0)[int(14)] = _S1134;
    (*v_sh_coeffs_0)[int(15)] = _S1135;
    *v_R_1 = _S1245;
    *v_t_1 = _S1244;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S1249;
    DiffPair_float_0 _S1250;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S1251;
    DiffPair_float_0 _S1252;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1253;
    DiffPair_float_0 _S1254;
    DiffPair_float_0 _S1255;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1256;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S1257;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1258;
};

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S1259, float _S1260)
{
    return s_primal_ctx_atan2_0(_S1259, _S1260);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S1261;
    DiffPair_float_0 _S1262;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S1263 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S1261 = _S1263;
    _s_diff_ctx_2->_S1262 = _S1263;
    (&_s_diff_ctx_2->_S1261)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S1261)->differential_0 = 0.0f;
    (&_s_diff_ctx_2->_S1262)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S1262)->differential_0 = 0.0f;
    DiffPair_float_0 _S1264 = *dpdpy_0;
    _s_diff_ctx_2->_S1261 = *dpdpy_0;
    DiffPair_float_0 _S1265 = *dpdpx_0;
    _s_diff_ctx_2->_S1262 = *dpdpx_0;
    float _S1266 = _S1265.primal_0 * _S1265.primal_0 + _S1264.primal_0 * _S1264.primal_0;
    float _S1267 = - _S1264.primal_0 / _S1266 * dpdOut_0;
    float _S1268 = _S1265.primal_0 / _S1266 * dpdOut_0;
    dpdpy_0->primal_0 = _S1264.primal_0;
    dpdpy_0->differential_0 = _S1268;
    dpdpx_0->primal_0 = _S1265.primal_0;
    dpdpx_0->differential_0 = _S1267;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S1269, DiffPair_float_0 * _S1270, float _S1271, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S1272 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S1249 = _S1272;
    _s_diff_ctx_3->_S1250 = _S1272;
    (&_s_diff_ctx_3->_S1249)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S1249)->differential_0 = 0.0f;
    (&_s_diff_ctx_3->_S1250)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S1250)->differential_0 = 0.0f;
    DiffPair_float_0 _S1273 = *_S1269;
    _s_diff_ctx_3->_S1249 = *_S1269;
    DiffPair_float_0 _S1274 = *_S1270;
    _s_diff_ctx_3->_S1250 = *_S1270;
    DiffPair_float_0 _S1275 = _S1273;
    DiffPair_float_0 _S1276 = _S1274;
    s_bwd_prop_d_atan2_Intermediates_0 _S1277;
    (&_S1277)->_S1261 = _S1272;
    (&_S1277)->_S1262 = _S1272;
    s_primal_ctx_d_atan2_0(&_S1275, &_S1276, _S1271, &_S1277);
    *_S1269 = _S1275;
    *_S1270 = _S1276;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1278;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1279;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1280;
    DiffPair_float_0 _S1281;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1282;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1283;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S1284 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S1283 = _S1284;
    (&_s_diff_ctx_4->_S1283)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S1283)->differential_0 = 0.0f;
    DiffPair_float_0 _S1285 = *dpdpx_1;
    _s_diff_ctx_4->_S1283 = *dpdpx_1;
    float _S1286 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S1285.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S1285.primal_0;
    dpdpx_1->differential_0 = _S1286;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S1287, float _S1288, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S1289 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S1279 = _S1289;
    (&_s_diff_ctx_5->_S1279)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S1279)->differential_0 = 0.0f;
    DiffPair_float_0 _S1290 = *_S1287;
    _s_diff_ctx_5->_S1279 = *_S1287;
    DiffPair_float_0 _S1291 = _S1290;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1292;
    (&_S1292)->_S1283 = _S1289;
    s_primal_ctx_d_sqrt_0(&_S1291, _S1288, &_S1292);
    *_S1287 = _S1291;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_6)
{
    float2  _S1293 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1294 = { _S1293, _S1293 };
    DiffPair_float_0 _S1295 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1296 = { _S1295 };
    _s_diff_ctx_6->_S1280 = _S1294;
    _s_diff_ctx_6->_S1281 = _S1295;
    _s_diff_ctx_6->_S1282 = _S1296;
    (&_s_diff_ctx_6->_S1280)->primal_0 = _S1293;
    (&_s_diff_ctx_6->_S1280)->differential_0 = _S1293;
    (&_s_diff_ctx_6->_S1281)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1281)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1297 = *dpdpx_2;
    _s_diff_ctx_6->_S1280 = *dpdpx_2;
    float _S1298 = _S1297.primal_0.x;
    float _S1299 = _S1297.primal_0.y;
    DiffPair_float_0 _S1300;
    (&_S1300)->primal_0 = _S1298 * _S1298 + _S1299 * _S1299;
    (&_S1300)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S1300, dp_s_dOut_0, &_s_diff_ctx_6->_S1282);
    _s_diff_ctx_6->_S1281 = _S1300;
    float _S1301 = _S1297.primal_0.y * _S1300.differential_0;
    float _S1302 = _S1301 + _S1301;
    float _S1303 = _S1297.primal_0.x * _S1300.differential_0;
    float _S1304 = _S1303 + _S1303;
    float2  _S1305 = _S1293;
    *&((&_S1305)->y) = _S1302;
    *&((&_S1305)->x) = _S1304;
    dpdpx_2->primal_0 = _S1297.primal_0;
    dpdpx_2->differential_0 = _S1305;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1306, float _S1307, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_7)
{
    float2  _S1308 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1309 = { _S1308, _S1308 };
    _s_diff_ctx_7->_S1278 = _S1309;
    (&_s_diff_ctx_7->_S1278)->primal_0 = _S1308;
    (&_s_diff_ctx_7->_S1278)->differential_0 = _S1308;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1310 = *_S1306;
    _s_diff_ctx_7->_S1278 = *_S1306;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1311 = _S1310;
    DiffPair_float_0 _S1312 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1313 = { _S1312 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1314;
    (&_S1314)->_S1280 = _S1309;
    (&_S1314)->_S1281 = _S1312;
    (&_S1314)->_S1282 = _S1313;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1311, _S1307, &_S1314);
    *_S1306 = _S1311;
    return;
}

inline __device__ bool s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float4  dpintrins_0, FixedArray<float, 10>  * dpdist_coeffs_0, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_float_0 _S1315 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1316 = { _S1315, _S1315 };
    _s_diff_ctx_8->_S1251 = _S1315;
    _s_diff_ctx_8->_S1252 = _S1315;
    _s_diff_ctx_8->_S1253 = _S1316;
    _s_diff_ctx_8->_S1254 = _S1315;
    _s_diff_ctx_8->_S1255 = _S1315;
    _s_diff_ctx_8->_S1256 = _S1316;
    (&_s_diff_ctx_8->_S1251)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1251)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1252)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1252)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1254)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1254)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1255)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1255)->differential_0 = 0.0f;
    float2  _S1317 = make_float2 (0.0f);
    float2  _S1318 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S1319 = length_0(_S1318);
    float _S1320 = dpmean3d_0.z;
    float _S1321 = s_primal_ctx_atan2_0(_S1319, _S1320);
    float k_4;
    if(_S1321 < 0.00100000004749745f)
    {
        k_4 = (1.0f - _S1321 * _S1321 / 3.0f) / _S1320;
    }
    else
    {
        k_4 = _S1321 / _S1319;
    }
    float2  _S1322 = _S1318 * make_float2 (k_4);
    float u_38 = _S1322.x;
    float v_38 = _S1322.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float _S1323 = 2.0f * (*dpdist_coeffs_0)[int(4)];
    float _S1324 = 2.0f * (*dpdist_coeffs_0)[int(5)];
    float2  _S1325 = _S1322 * make_float2 (1.0f + r2_38 * ((*dpdist_coeffs_0)[int(0)] + r2_38 * ((*dpdist_coeffs_0)[int(1)] + r2_38 * ((*dpdist_coeffs_0)[int(2)] + r2_38 * (*dpdist_coeffs_0)[int(3)])))) + make_float2 (_S1323 * u_38 * v_38 + (*dpdist_coeffs_0)[int(5)] * (r2_38 + 2.0f * u_38 * u_38) + (*dpdist_coeffs_0)[int(6)] * r2_38, _S1324 * u_38 * v_38 + (*dpdist_coeffs_0)[int(4)] * (r2_38 + 2.0f * v_38 * v_38) + (*dpdist_coeffs_0)[int(7)] * r2_38);
    float2  _S1326 = _S1325 + make_float2 ((*dpdist_coeffs_0)[int(8)] * _S1325.x + (*dpdist_coeffs_0)[int(9)] * _S1325.y, 0.0f);
    float fx_18 = dpintrins_0.x;
    float fy_18 = dpintrins_0.y;
    float2  _S1327 = make_float2 (fx_18 * _S1326.x + dpintrins_0.z, fy_18 * _S1326.y + dpintrins_0.w);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (0.0f);
    float _S1328 = s_primal_ctx_s_primal_ctx_atan2_0(_S1319, _S1320);
    bool _S1329 = _S1328 < 0.00100000004749745f;
    float _S1330;
    float _S1331;
    float _S1332;
    if(_S1329)
    {
        float _S1333 = 1.0f - _S1328 * _S1328 / 3.0f;
        float _S1334 = _S1320 * _S1320;
        k_4 = _S1333 / _S1320;
        _S1330 = 0.0f;
        _S1331 = _S1334;
        _S1332 = _S1333;
    }
    else
    {
        float _S1335 = _S1319 * _S1319;
        k_4 = _S1328 / _S1319;
        _S1330 = _S1335;
        _S1331 = 0.0f;
        _S1332 = 0.0f;
    }
    float2  _S1336 = make_float2 (k_4);
    float2  _S1337 = _S1318 * make_float2 (k_4);
    float u_39 = _S1337.x;
    float v_39 = _S1337.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float _S1338 = (*dpdist_coeffs_0)[int(2)] + r2_39 * (*dpdist_coeffs_0)[int(3)];
    float _S1339 = (*dpdist_coeffs_0)[int(1)] + r2_39 * _S1338;
    float _S1340 = (*dpdist_coeffs_0)[int(0)] + r2_39 * _S1339;
    float2  _S1341 = make_float2 (fx_18, 0.0f) + make_float2 ((*dpdist_coeffs_0)[int(8)] * fx_18, (*dpdist_coeffs_0)[int(9)] * fx_18);
    float2  _S1342 = _S1337 * _S1341;
    float _S1343 = (*dpdist_coeffs_0)[int(4)] * _S1341.y;
    float _S1344 = (*dpdist_coeffs_0)[int(5)] * _S1341.x;
    float _S1345 = _S1342.x + _S1342.y;
    float _S1346 = r2_39 * _S1345;
    float _S1347 = r2_39 * _S1346;
    float _S1348 = (*dpdist_coeffs_0)[int(7)] * _S1341.y + _S1343 + (*dpdist_coeffs_0)[int(6)] * _S1341.x + _S1344 + _S1340 * _S1345 + _S1339 * _S1346 + _S1338 * _S1347 + (*dpdist_coeffs_0)[int(3)] * (r2_39 * _S1347);
    float _S1349 = v_39 * _S1348;
    float _S1350 = u_39 * _S1348;
    float2  _S1351 = make_float2 (1.0f + r2_39 * _S1340) * _S1341 + make_float2 (_S1324 * (v_39 * _S1341.y) + 2.0f * u_39 * _S1344 + 2.0f * (u_39 * _S1344) + _S1323 * (v_39 * _S1341.x) + _S1350 + _S1350, 2.0f * v_39 * _S1343 + 2.0f * (v_39 * _S1343) + _S1324 * u_39 * _S1341.y + _S1323 * u_39 * _S1341.x + _S1349 + _S1349);
    float2  _S1352 = _S1318 * _S1351;
    float2  _S1353 = _S1336 * _S1351;
    float _S1354 = _S1352.x + _S1352.y;
    if(_S1329)
    {
        float _S1355 = _S1354 / _S1331;
        float _S1356 = _S1332 * - _S1355;
        float _S1357 = _S1328 * (0.3333333432674408f * - (_S1320 * _S1355));
        k_4 = _S1357 + _S1357;
        _S1330 = _S1356;
        _S1331 = 0.0f;
    }
    else
    {
        float _S1358 = _S1354 / _S1330;
        float _S1359 = _S1328 * - _S1358;
        k_4 = _S1319 * _S1358;
        _S1330 = 0.0f;
        _S1331 = _S1359;
    }
    DiffPair_float_0 _S1360;
    (&_S1360)->primal_0 = _S1319;
    (&_S1360)->differential_0 = 0.0f;
    DiffPair_float_0 _S1361;
    (&_S1361)->primal_0 = _S1320;
    (&_S1361)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1360, &_S1361, k_4, &_s_diff_ctx_8->_S1253);
    _s_diff_ctx_8->_S1251 = _S1360;
    _s_diff_ctx_8->_S1252 = _S1361;
    float _S1362 = _S1361.differential_0 + _S1330;
    float _S1363 = _S1360.differential_0 + _S1331;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1364;
    (&_S1364)->primal_0 = _S1318;
    (&_S1364)->differential_0 = _S1317;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1365 = { _S1317, _S1317 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1366;
    (&_S1366)->_S1278 = _S1365;
    s_primal_ctx_s_bwd_length_impl_0(&_S1364, _S1363, &_S1366);
    float2  _S1367 = _S1364.differential_0 + _S1353;
    float3  _S1368 = make_float3 (_S1367.x, _S1367.y, _S1362);
    Matrix<float, 2, 3>  _S1369 = J_11;
    _S1369[int(0)] = _S1368;
    if(_S1329)
    {
        float _S1370 = 1.0f - _S1328 * _S1328 / 3.0f;
        float _S1371 = _S1320 * _S1320;
        k_4 = _S1370 / _S1320;
        _S1330 = 0.0f;
        _S1331 = _S1371;
        _S1332 = _S1370;
    }
    else
    {
        float _S1372 = _S1319 * _S1319;
        k_4 = _S1328 / _S1319;
        _S1330 = _S1372;
        _S1331 = 0.0f;
        _S1332 = 0.0f;
    }
    float2  _S1373 = make_float2 (k_4);
    float2  _S1374 = _S1318 * make_float2 (k_4);
    float u_40 = _S1374.x;
    float v_40 = _S1374.y;
    float r2_40 = u_40 * u_40 + v_40 * v_40;
    float _S1375 = (*dpdist_coeffs_0)[int(2)] + r2_40 * (*dpdist_coeffs_0)[int(3)];
    float _S1376 = (*dpdist_coeffs_0)[int(1)] + r2_40 * _S1375;
    float _S1377 = (*dpdist_coeffs_0)[int(0)] + r2_40 * _S1376;
    float2  _S1378 = make_float2 (0.0f, fy_18);
    float2  _S1379 = _S1374 * _S1378;
    float _S1380 = (*dpdist_coeffs_0)[int(4)] * fy_18;
    float _S1381 = _S1379.x + _S1379.y;
    float _S1382 = r2_40 * _S1381;
    float _S1383 = r2_40 * _S1382;
    float _S1384 = (*dpdist_coeffs_0)[int(7)] * fy_18 + _S1380 + _S1377 * _S1381 + _S1376 * _S1382 + _S1375 * _S1383 + (*dpdist_coeffs_0)[int(3)] * (r2_40 * _S1383);
    float _S1385 = v_40 * _S1384;
    float _S1386 = u_40 * _S1384;
    float2  _S1387 = make_float2 (1.0f + r2_40 * _S1377) * _S1378 + make_float2 (_S1324 * (v_40 * fy_18) + _S1386 + _S1386, 2.0f * v_40 * _S1380 + 2.0f * (v_40 * _S1380) + _S1324 * u_40 * fy_18 + _S1385 + _S1385);
    float2  _S1388 = _S1318 * _S1387;
    float2  _S1389 = _S1373 * _S1387;
    float _S1390 = _S1388.x + _S1388.y;
    if(_S1329)
    {
        float _S1391 = _S1390 / _S1331;
        float _S1392 = _S1332 * - _S1391;
        float _S1393 = _S1328 * (0.3333333432674408f * - (_S1320 * _S1391));
        k_4 = _S1393 + _S1393;
        _S1330 = _S1392;
        _S1331 = 0.0f;
    }
    else
    {
        float _S1394 = _S1390 / _S1330;
        float _S1395 = _S1328 * - _S1394;
        k_4 = _S1319 * _S1394;
        _S1330 = 0.0f;
        _S1331 = _S1395;
    }
    DiffPair_float_0 _S1396;
    (&_S1396)->primal_0 = _S1319;
    (&_S1396)->differential_0 = 0.0f;
    DiffPair_float_0 _S1397;
    (&_S1397)->primal_0 = _S1320;
    (&_S1397)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1396, &_S1397, k_4, &_s_diff_ctx_8->_S1256);
    _s_diff_ctx_8->_S1254 = _S1396;
    _s_diff_ctx_8->_S1255 = _S1397;
    float _S1398 = _S1397.differential_0 + _S1330;
    float _S1399 = _S1396.differential_0 + _S1331;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1400;
    (&_S1400)->primal_0 = _S1318;
    (&_S1400)->differential_0 = _S1317;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1401;
    (&_S1401)->_S1278 = _S1365;
    s_primal_ctx_s_bwd_length_impl_0(&_S1400, _S1399, &_S1401);
    float2  _S1402 = _S1400.differential_0 + _S1389;
    float3  _S1403 = make_float3 (_S1402.x, _S1402.y, _S1398);
    _S1369[int(1)] = _S1403;
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S1369, dpcov3d_0), transpose_1(_S1369));
    *dpmean2d_0 = _S1327;
    return true;
}

struct DiffPair_vectorx3Cfloatx2C4x3E_0
{
    float4  primal_0;
    float4  differential_0;
};

struct DiffPair_arrayx3Cfloatx2C10x3E_0
{
    FixedArray<float, 10>  primal_0;
    FixedArray<float, 10>  differential_0;
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
    DiffPair_1 _S1404 = *dpdpx_3;
    float _S1405 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_9->_S1283)->primal_0);
    float _S1406 = s_primal_ctx_sqrt_0(_S1405);
    float _S1407 = 0.5f / _S1406 * (*dpdpx_3).differential_0.differential_0;
    float _S1408 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S1406 * _S1406));
    DiffPair_float_0 _S1409;
    (&_S1409)->primal_0 = _S1405;
    (&_S1409)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1409, _S1408);
    DiffPair_float_0 _S1410;
    (&_S1410)->primal_0 = 1.00000001168609742e-07f;
    (&_S1410)->differential_0 = 0.0f;
    DiffPair_float_0 _S1411;
    (&_S1411)->primal_0 = (&_s_diff_ctx_9->_S1283)->primal_0;
    (&_S1411)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1410, &_S1411, _S1409.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S1411.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S1407;
    dpdpx_3->primal_0 = _S1404.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S1412, DiffPair_float_0 * _S1413, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_1 _S1414 = *_S1412;
    DiffPair_float_0 _S1415 = _s_diff_ctx_10->_S1279;
    DiffPair_float_0 _S1416 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S1417;
    (&_S1417)->_S1283 = _S1416;
    s_primal_ctx_d_sqrt_0(&_S1415, (*_S1413).primal_0, &_S1417);
    DiffPair_float_0 _S1418 = { (*_S1412).differential_0.primal_0, (*_S1412).differential_0.differential_0 };
    DiffPair_1 _S1419;
    (&_S1419)->primal_0 = _s_diff_ctx_10->_S1279;
    (&_S1419)->differential_0 = _S1418;
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = (*_S1413).primal_0;
    (&_S1420)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1421 = _S1417;
    s_bwd_prop_d_sqrt_0(&_S1419, &_S1420, &_S1421);
    DiffPair_float_0 _S1422 = { _S1419.differential_0.primal_0, _S1419.differential_0.differential_0 };
    _S1413->primal_0 = (*_S1413).primal_0;
    _S1413->differential_0 = _S1420.differential_0;
    _S1412->primal_0 = _S1414.primal_0;
    _S1412->differential_0 = _S1422;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S1423, float _s_dOut_3)
{
    DiffPair_float_0 _S1424;
    (&_S1424)->primal_0 = (*_S1423).primal_0;
    (&_S1424)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1424, _s_dOut_3);
    _S1423->primal_0 = (*_S1423).primal_0;
    _S1423->differential_0 = _S1424.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_0 _S1425 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->x) * *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->x) + *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->y) * *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->y);
    DiffPair_float_0 _S1426 = { len_0, 0.0f };
    float2  _S1427 = make_float2 (0.0f);
    float _S1428 = (*dpdpx_5).differential_0.differential_0.x;
    float _S1429 = _S1428 + _S1428;
    float _S1430 = (&_s_diff_ctx_11->_S1281)->differential_0 * _S1429;
    float _S1431 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S1432 = (&_s_diff_ctx_11->_S1281)->differential_0 * _S1431;
    DiffPair_float_0 _S1433 = { 0.0f, *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->x) * _S1429 + *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->y) * _S1431 };
    DiffPair_1 _S1434;
    (&_S1434)->primal_0 = _S1426;
    (&_S1434)->differential_0 = _S1433;
    DiffPair_float_0 _S1435;
    (&_S1435)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S1435)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S1434, &_S1435, &_s_diff_ctx_11->_S1282);
    DiffPair_float_0 _S1436;
    (&_S1436)->primal_0 = len_0;
    (&_S1436)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1436, 0.0f);
    float _S1437 = _S1434.differential_0.primal_0 + _S1436.differential_0;
    float _S1438 = *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->y) * _S1437;
    float _S1439 = _S1432 + _S1438 + _S1438;
    float _S1440 = *&((&(&_s_diff_ctx_11->_S1280)->primal_0)->x) * _S1437;
    float _S1441 = _S1430 + _S1440 + _S1440;
    float2  _S1442 = _S1427;
    *&((&_S1442)->y) = _S1439;
    *&((&_S1442)->x) = _S1441;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S1425.differential_0.primal_0 + _S1442, _S1427 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S1435.differential_0;
    dpdpx_5->primal_0 = _S1425.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S1443 = (*dpdpx_7).primal_0.x;
    float _S1444 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S1445;
    (&_S1445)->primal_0 = _S1443 * _S1443 + _S1444 * _S1444;
    (&_S1445)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1445, _s_dOut_4);
    float _S1446 = (*dpdpx_7).primal_0.y * _S1445.differential_0;
    float _S1447 = _S1446 + _S1446;
    float _S1448 = (*dpdpx_7).primal_0.x * _S1445.differential_0;
    float _S1449 = _S1448 + _S1448;
    float2  _S1450 = make_float2 (0.0f);
    *&((&_S1450)->y) = _S1447;
    *&((&_S1450)->x) = _S1449;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S1450;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S1451, DiffPair_float_0 * _S1452, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_0 _S1453 = *_S1451;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1454 = _s_diff_ctx_12->_S1278;
    float2  _S1455 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1456 = { _S1455, _S1455 };
    DiffPair_float_0 _S1457 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1458 = { _S1457 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1459;
    (&_S1459)->_S1280 = _S1456;
    (&_S1459)->_S1281 = _S1457;
    (&_S1459)->_S1282 = _S1458;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1454, (*_S1452).primal_0, &_S1459);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1460 = { (*_S1451).differential_0.primal_0, (*_S1451).differential_0.differential_0 };
    DiffPair_0 _S1461;
    (&_S1461)->primal_0 = _s_diff_ctx_12->_S1278;
    (&_S1461)->differential_0 = _S1460;
    DiffPair_float_0 _S1462;
    (&_S1462)->primal_0 = (*_S1452).primal_0;
    (&_S1462)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1463 = _S1459;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S1461, &_S1462, &_S1463);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1464;
    (&_S1464)->primal_0 = (&_s_diff_ctx_12->_S1278)->primal_0;
    (&_S1464)->differential_0 = _S1455;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S1464, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1465 = { _S1461.differential_0.primal_0 + _S1464.differential_0, _S1461.differential_0.differential_0 };
    _S1452->primal_0 = (*_S1452).primal_0;
    _S1452->differential_0 = _S1462.differential_0;
    _S1451->primal_0 = _S1453.primal_0;
    _S1451->differential_0 = _S1465;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_1 _S1466 = *dpdpy_1;
    DiffPair_1 _S1467 = *dpdpx_8;
    float _S1468 = - (&_s_diff_ctx_13->_S1261)->primal_0;
    float _S1469 = (&_s_diff_ctx_13->_S1262)->primal_0 * (&_s_diff_ctx_13->_S1262)->primal_0 + (&_s_diff_ctx_13->_S1261)->primal_0 * (&_s_diff_ctx_13->_S1261)->primal_0;
    float _S1470 = _S1469 * _S1469;
    float _S1471 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S1470;
    float _S1472 = (&_s_diff_ctx_13->_S1262)->primal_0 * - _S1471;
    float _S1473 = (&_s_diff_ctx_13->_S1261)->primal_0 * _S1472;
    float _S1474 = (&_s_diff_ctx_13->_S1262)->primal_0 * _S1472;
    float _S1475 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S1470;
    float _S1476 = _S1468 * - _S1475;
    float _S1477 = (&_s_diff_ctx_13->_S1261)->primal_0 * _S1476;
    float _S1478 = (&_s_diff_ctx_13->_S1262)->primal_0 * _S1476;
    DiffPair_float_0 dpdpx_9 = { _S1478 + _S1478 + ((*dpdpx_8).differential_0.primal_0 + (_S1474 + _S1474 + _S1469 * _S1471)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S1473 + _S1473 + (*dpdpy_1).differential_0.primal_0 + _S1477 + _S1477 + - (_S1469 * _S1475), 0.0f };
    float _S1479 = (&_s_diff_ctx_13->_S1262)->primal_0 / _S1469 * (*dpdpy_1).differential_0.differential_0 + _S1468 / _S1469 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S1479;
    dpdpy_1->primal_0 = _S1466.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S1467.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S1480, DiffPair_1 * _S1481, DiffPair_float_0 * _S1482, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_1 _S1483 = *_S1480;
    DiffPair_1 _S1484 = *_S1481;
    DiffPair_float_0 _S1485 = _s_diff_ctx_14->_S1249;
    DiffPair_float_0 _S1486 = _s_diff_ctx_14->_S1250;
    DiffPair_float_0 _S1487 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S1488;
    (&_S1488)->_S1261 = _S1487;
    (&_S1488)->_S1262 = _S1487;
    s_primal_ctx_d_atan2_0(&_S1485, &_S1486, (*_S1482).primal_0, &_S1488);
    DiffPair_float_0 _S1489 = { (*_S1481).differential_0.primal_0, (*_S1481).differential_0.differential_0 };
    DiffPair_float_0 _S1490 = { (*_S1480).differential_0.primal_0, (*_S1480).differential_0.differential_0 };
    DiffPair_1 _S1491;
    (&_S1491)->primal_0 = _s_diff_ctx_14->_S1249;
    (&_S1491)->differential_0 = _S1490;
    DiffPair_1 _S1492;
    (&_S1492)->primal_0 = _s_diff_ctx_14->_S1250;
    (&_S1492)->differential_0 = _S1489;
    DiffPair_float_0 _S1493;
    (&_S1493)->primal_0 = (*_S1482).primal_0;
    (&_S1493)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S1494 = _S1488;
    s_bwd_prop_d_atan2_0(&_S1491, &_S1492, &_S1493, &_S1494);
    DiffPair_float_0 _S1495 = { _S1492.differential_0.primal_0, _S1492.differential_0.differential_0 };
    DiffPair_float_0 _S1496 = { _S1491.differential_0.primal_0, _S1491.differential_0.differential_0 };
    _S1482->primal_0 = (*_S1482).primal_0;
    _S1482->differential_0 = _S1493.differential_0;
    _S1480->primal_0 = _S1483.primal_0;
    _S1480->differential_0 = _S1496;
    _S1481->primal_0 = _S1484.primal_0;
    _S1481->differential_0 = _S1495;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S1497, DiffPair_float_0 * _S1498, float _s_dOut_5)
{
    DiffPair_float_0 _S1499;
    (&_S1499)->primal_0 = (*_S1497).primal_0;
    (&_S1499)->differential_0 = 0.0f;
    DiffPair_float_0 _S1500;
    (&_S1500)->primal_0 = (*_S1498).primal_0;
    (&_S1500)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1499, &_S1500, _s_dOut_5);
    _S1498->primal_0 = (*_S1498).primal_0;
    _S1498->differential_0 = _S1500.differential_0;
    _S1497->primal_0 = (*_S1497).primal_0;
    _S1497->differential_0 = _S1499.differential_0;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpintrins_1, DiffPair_arrayx3Cfloatx2C10x3E_0 * dpdist_coeffs_1, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1501 = *dpcov3d_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1502 = *dpintrins_1;
    FixedArray<float, 10>  _S1503 = dpdist_coeffs_1->primal_0;
    float2  _S1504 = make_float2 (0.0f);
    float2  _S1505 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S1506 = length_0(_S1505);
    float _S1507 = (*dpmean3d_1).primal_0.z;
    float _S1508 = s_primal_ctx_atan2_0(_S1506, _S1507);
    bool _S1509 = _S1508 < 0.00100000004749745f;
    float k_5;
    float _S1510;
    float _S1511;
    float _S1512;
    if(_S1509)
    {
        float _S1513 = 1.0f - _S1508 * _S1508 / 3.0f;
        float _S1514 = _S1507 * _S1507;
        k_5 = _S1513 / _S1507;
        _S1510 = 0.0f;
        _S1511 = _S1514;
        _S1512 = _S1513;
    }
    else
    {
        float _S1515 = _S1506 * _S1506;
        k_5 = _S1508 / _S1506;
        _S1510 = _S1515;
        _S1511 = 0.0f;
        _S1512 = 0.0f;
    }
    float2  _S1516 = make_float2 (k_5);
    float2  _S1517 = _S1505 * make_float2 (k_5);
    float u_41 = _S1517.x;
    float v_41 = _S1517.y;
    float r2_41 = u_41 * u_41 + v_41 * v_41;
    float _S1518 = _S1503[int(2)] + r2_41 * _S1503[int(3)];
    float _S1519 = _S1503[int(1)] + r2_41 * _S1518;
    float _S1520 = _S1503[int(0)] + r2_41 * _S1519;
    float radial_1 = 1.0f + r2_41 * _S1520;
    float2  _S1521 = make_float2 (radial_1);
    float _S1522 = 2.0f * _S1503[int(4)];
    float _S1523 = _S1522 * u_41;
    float _S1524 = 2.0f * u_41;
    float _S1525 = r2_41 + _S1524 * u_41;
    float _S1526 = 2.0f * _S1503[int(5)];
    float _S1527 = _S1526 * u_41;
    float _S1528 = 2.0f * v_41;
    float _S1529 = r2_41 + _S1528 * v_41;
    float2  _S1530 = _S1517 * make_float2 (radial_1) + make_float2 (_S1523 * v_41 + _S1503[int(5)] * _S1525 + _S1503[int(6)] * r2_41, _S1527 * v_41 + _S1503[int(4)] * _S1529 + _S1503[int(7)] * r2_41);
    float _S1531 = _S1530.x;
    float _S1532 = _S1530.y;
    float2  _S1533 = _S1530 + make_float2 (_S1503[int(8)] * _S1531 + _S1503[int(9)] * _S1532, 0.0f);
    float fx_19 = _S1502.primal_0.x;
    float fy_19 = _S1502.primal_0.y;
    float _S1534 = _S1533.x;
    float _S1535 = _S1533.y;
    Matrix<float, 2, 3>  J_12 = makeMatrix<float, 2, 3> (0.0f);
    float _S1536 = s_primal_ctx_s_primal_ctx_atan2_0(_S1506, _S1507);
    bool _S1537 = _S1536 < 0.00100000004749745f;
    float _S1538;
    float _S1539;
    float _S1540;
    if(_S1537)
    {
        float _S1541 = 1.0f - _S1536 * _S1536 / 3.0f;
        float _S1542 = _S1507 * _S1507;
        k_5 = _S1541 / _S1507;
        _S1538 = 0.0f;
        _S1539 = _S1542;
        _S1540 = _S1541;
    }
    else
    {
        float _S1543 = _S1506 * _S1506;
        k_5 = _S1536 / _S1506;
        _S1538 = _S1543;
        _S1539 = 0.0f;
        _S1540 = 0.0f;
    }
    float2  _S1544 = make_float2 (k_5);
    float2  _S1545 = _S1505 * make_float2 (k_5);
    float u_42 = _S1545.x;
    float v_42 = _S1545.y;
    float r2_42 = u_42 * u_42 + v_42 * v_42;
    float _S1546 = _S1503[int(2)] + r2_42 * _S1503[int(3)];
    float _S1547 = _S1503[int(1)] + r2_42 * _S1546;
    float _S1548 = _S1503[int(0)] + r2_42 * _S1547;
    float2  _S1549 = make_float2 (1.0f + r2_42 * _S1548);
    float _S1550 = _S1522 * u_42;
    float _S1551 = 2.0f * u_42;
    float _S1552 = _S1526 * u_42;
    float _S1553 = 2.0f * v_42;
    float2  _S1554 = make_float2 (fx_19, 0.0f) + make_float2 (_S1503[int(8)] * fx_19, _S1503[int(9)] * fx_19);
    float2  _S1555 = _S1545 * _S1554;
    float _S1556 = _S1503[int(4)] * _S1554.y;
    float _S1557 = v_42 * _S1554.y;
    float _S1558 = _S1503[int(5)] * _S1554.x;
    float _S1559 = v_42 * _S1554.x;
    float _S1560 = _S1555.x + _S1555.y;
    float _S1561 = r2_42 * _S1560;
    float _S1562 = r2_42 * _S1561;
    float _S1563 = r2_42 * _S1562;
    float _S1564 = _S1503[int(7)] * _S1554.y + _S1556 + _S1503[int(6)] * _S1554.x + _S1558 + _S1548 * _S1560 + _S1547 * _S1561 + _S1546 * _S1562 + _S1503[int(3)] * _S1563;
    float _S1565 = v_42 * _S1564;
    float _S1566 = u_42 * _S1564;
    float2  _S1567 = _S1549 * _S1554 + make_float2 (_S1526 * _S1557 + _S1551 * _S1558 + 2.0f * (u_42 * _S1558) + _S1522 * _S1559 + _S1566 + _S1566, _S1553 * _S1556 + 2.0f * (v_42 * _S1556) + _S1552 * _S1554.y + _S1550 * _S1554.x + _S1565 + _S1565);
    float2  _S1568 = _S1505 * _S1567;
    float2  _S1569 = _S1544 * _S1567;
    float _S1570 = _S1568.x + _S1568.y;
    float k_6;
    float _S1571;
    float _S1572;
    float _S1573;
    float _S1574;
    float _S1575;
    float _S1576;
    float _S1577;
    float _S1578;
    if(_S1537)
    {
        float _S1579 = _S1570 / _S1539;
        float _S1580 = _S1539 * _S1539;
        float _S1581 = - _S1579;
        float _S1582 = _S1540 * _S1581;
        float _S1583 = 0.3333333432674408f * - (_S1507 * _S1579);
        float _S1584 = _S1536 * _S1583;
        k_5 = _S1584 + _S1584;
        k_6 = _S1582;
        _S1571 = 0.0f;
        _S1572 = 0.0f;
        _S1573 = 0.0f;
        _S1574 = 0.0f;
        _S1575 = _S1583;
        _S1576 = _S1579;
        _S1577 = _S1581;
        _S1578 = _S1580;
    }
    else
    {
        float _S1585 = _S1570 / _S1538;
        float _S1586 = _S1538 * _S1538;
        float _S1587 = - _S1585;
        float _S1588 = _S1536 * _S1587;
        k_5 = _S1506 * _S1585;
        k_6 = 0.0f;
        _S1571 = _S1588;
        _S1572 = _S1585;
        _S1573 = _S1587;
        _S1574 = _S1586;
        _S1575 = 0.0f;
        _S1576 = 0.0f;
        _S1577 = 0.0f;
        _S1578 = 0.0f;
    }
    DiffPair_float_0 _S1589 = { _S1506, 0.0f };
    DiffPair_float_0 _S1590 = { _S1507, 0.0f };
    float _S1591 = (&_s_diff_ctx_15->_S1252)->differential_0 + k_6;
    float _S1592 = (&_s_diff_ctx_15->_S1251)->differential_0 + _S1571;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1593 = { _S1505, _S1504 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1594;
    (&_S1594)->primal_0 = _S1505;
    (&_S1594)->differential_0 = _S1504;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1595 = { _S1504, _S1504 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1596;
    (&_S1596)->_S1278 = _S1595;
    s_primal_ctx_s_bwd_length_impl_0(&_S1594, _S1592, &_S1596);
    float2  _S1597 = _S1594.differential_0 + _S1569;
    float3  _S1598 = make_float3 (_S1597.x, _S1597.y, _S1591);
    Matrix<float, 2, 3>  _S1599 = J_12;
    _S1599[int(0)] = _S1598;
    float _S1600;
    float _S1601;
    if(_S1537)
    {
        float _S1602 = 1.0f - _S1536 * _S1536 / 3.0f;
        float _S1603 = _S1507 * _S1507;
        k_6 = _S1602 / _S1507;
        _S1571 = 0.0f;
        _S1600 = _S1603;
        _S1601 = _S1602;
    }
    else
    {
        float _S1604 = _S1506 * _S1506;
        k_6 = _S1536 / _S1506;
        _S1571 = _S1604;
        _S1600 = 0.0f;
        _S1601 = 0.0f;
    }
    float2  _S1605 = make_float2 (k_6);
    float2  _S1606 = _S1505 * make_float2 (k_6);
    float u_43 = _S1606.x;
    float v_43 = _S1606.y;
    float r2_43 = u_43 * u_43 + v_43 * v_43;
    float _S1607 = _S1503[int(2)] + r2_43 * _S1503[int(3)];
    float _S1608 = _S1503[int(1)] + r2_43 * _S1607;
    float _S1609 = _S1503[int(0)] + r2_43 * _S1608;
    float2  _S1610 = make_float2 (1.0f + r2_43 * _S1609);
    float _S1611 = _S1522 * u_43;
    float _S1612 = 2.0f * u_43;
    float _S1613 = _S1526 * u_43;
    float _S1614 = 2.0f * v_43;
    float2  _S1615 = make_float2 (0.0f, fy_19);
    float2  _S1616 = _S1606 * _S1615;
    float _S1617 = _S1503[int(4)] * fy_19;
    float _S1618 = v_43 * fy_19;
    float _S1619 = _S1616.x + _S1616.y;
    float _S1620 = r2_43 * _S1619;
    float _S1621 = r2_43 * _S1620;
    float _S1622 = r2_43 * _S1621;
    float _S1623 = _S1503[int(7)] * fy_19 + _S1617 + _S1609 * _S1619 + _S1608 * _S1620 + _S1607 * _S1621 + _S1503[int(3)] * _S1622;
    float _S1624 = v_43 * _S1623;
    float _S1625 = u_43 * _S1623;
    float2  _S1626 = _S1610 * _S1615 + make_float2 (_S1526 * _S1618 + _S1625 + _S1625, _S1614 * _S1617 + 2.0f * (v_43 * _S1617) + _S1613 * fy_19 + _S1624 + _S1624);
    float2  _S1627 = _S1505 * _S1626;
    float2  _S1628 = _S1605 * _S1626;
    float _S1629 = _S1627.x + _S1627.y;
    float _S1630;
    float _S1631;
    float _S1632;
    float _S1633;
    float _S1634;
    float _S1635;
    float _S1636;
    float _S1637;
    float _S1638;
    if(_S1537)
    {
        float _S1639 = _S1629 / _S1600;
        float _S1640 = _S1600 * _S1600;
        float _S1641 = - _S1639;
        float _S1642 = _S1601 * _S1641;
        float _S1643 = 0.3333333432674408f * - (_S1507 * _S1639);
        float _S1644 = _S1536 * _S1643;
        k_6 = _S1644 + _S1644;
        _S1630 = _S1642;
        _S1631 = 0.0f;
        _S1632 = 0.0f;
        _S1633 = 0.0f;
        _S1634 = 0.0f;
        _S1635 = _S1643;
        _S1636 = _S1639;
        _S1637 = _S1641;
        _S1638 = _S1640;
    }
    else
    {
        float _S1645 = _S1629 / _S1571;
        float _S1646 = _S1571 * _S1571;
        float _S1647 = - _S1645;
        float _S1648 = _S1536 * _S1647;
        k_6 = _S1506 * _S1645;
        _S1630 = 0.0f;
        _S1631 = _S1648;
        _S1632 = _S1645;
        _S1633 = _S1647;
        _S1634 = _S1646;
        _S1635 = 0.0f;
        _S1636 = 0.0f;
        _S1637 = 0.0f;
        _S1638 = 0.0f;
    }
    float _S1649 = (&_s_diff_ctx_15->_S1255)->differential_0 + _S1630;
    float _S1650 = (&_s_diff_ctx_15->_S1254)->differential_0 + _S1631;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1651;
    (&_S1651)->primal_0 = _S1505;
    (&_S1651)->differential_0 = _S1504;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1652;
    (&_S1652)->_S1278 = _S1595;
    s_primal_ctx_s_bwd_length_impl_0(&_S1651, _S1650, &_S1652);
    float2  _S1653 = _S1651.differential_0 + _S1628;
    float3  _S1654 = make_float3 (_S1653.x, _S1653.y, _S1649);
    _S1599[int(1)] = _S1654;
    Matrix<float, 3, 2>  _S1655 = transpose_1(_S1599);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1656;
    (&_S1656)->primal_0 = s_primal_ctx_mul_3(_S1599, _S1501.primal_0);
    (&_S1656)->differential_0 = J_12;
    Matrix<float, 3, 2>  _S1657 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1658;
    (&_S1658)->primal_0 = _S1655;
    (&_S1658)->differential_0 = _S1657;
    s_bwd_prop_mul_2(&_S1656, &_S1658, dpcov2d_1);
    Matrix<float, 2, 3>  _S1659 = transpose_2(_S1658.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1660;
    (&_S1660)->primal_0 = _S1599;
    (&_S1660)->differential_0 = J_12;
    Matrix<float, 3, 3>  _S1661 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1662;
    (&_S1662)->primal_0 = _S1501.primal_0;
    (&_S1662)->differential_0 = _S1661;
    s_bwd_prop_mul_3(&_S1660, &_S1662, _S1656.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1663 = _S1662;
    Matrix<float, 2, 3>  _S1664 = _S1659 + _S1660.differential_0;
    float2  _S1665 = _S1504;
    *&((&_S1665)->y) = _S1664.rows[int(1)].y;
    *&((&_S1665)->x) = _S1664.rows[int(1)].x;
    float2  _S1666 = _S1665;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1667 = { _S1504, _S1665 };
    DiffPair_0 _S1668;
    (&_S1668)->primal_0 = _S1593;
    (&_S1668)->differential_0 = _S1667;
    DiffPair_float_0 _S1669;
    (&_S1669)->primal_0 = _S1650;
    (&_S1669)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1670 = _S1652;
    s_bwd_prop_s_bwd_length_impl_0(&_S1668, &_S1669, &_S1670);
    DiffPair_0 _S1671 = _S1668;
    DiffPair_float_0 _S1672 = _S1669;
    DiffPair_float_0 _S1673 = { 0.0f, _S1664.rows[int(1)].z };
    DiffPair_float_0 _S1674 = { 0.0f, _S1669.differential_0 };
    DiffPair_1 _S1675;
    (&_S1675)->primal_0 = _S1589;
    (&_S1675)->differential_0 = _S1674;
    DiffPair_1 _S1676;
    (&_S1676)->primal_0 = _S1590;
    (&_S1676)->differential_0 = _S1673;
    DiffPair_float_0 _S1677;
    (&_S1677)->primal_0 = k_6;
    (&_S1677)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1675, &_S1676, &_S1677, &_s_diff_ctx_15->_S1256);
    DiffPair_1 _S1678 = _S1675;
    DiffPair_1 _S1679 = _S1676;
    DiffPair_float_0 _S1680 = _S1677;
    if(_S1537)
    {
        float _S1681 = _S1680.differential_0 + _S1680.differential_0;
        float _S1682 = _S1635 * _S1681;
        float _S1683 = - (0.3333333432674408f * (_S1536 * _S1681));
        float _S1684 = _S1637 * _S1664.rows[int(1)].z;
        float _S1685 = (_S1507 * _S1683 + - (_S1601 * _S1664.rows[int(1)].z)) / _S1638;
        float _S1686 = _S1629 * - _S1685;
        float _S1687 = _S1636 * _S1683 + _S1679.differential_0.primal_0;
        k_6 = _S1600 * _S1685;
        _S1630 = 0.0f;
        _S1631 = _S1686;
        _S1632 = _S1684;
        _S1633 = _S1678.differential_0.primal_0;
        _S1634 = _S1682;
        _S1635 = _S1687;
    }
    else
    {
        float _S1688 = _S1633 * _S1672.differential_0;
        float _S1689 = (_S1506 * _S1680.differential_0 + - (_S1536 * _S1672.differential_0)) / _S1634;
        float _S1690 = _S1629 * - _S1689;
        float _S1691 = _S1632 * _S1680.differential_0 + _S1678.differential_0.primal_0;
        k_6 = _S1571 * _S1689;
        _S1630 = _S1690;
        _S1631 = 0.0f;
        _S1632 = 0.0f;
        _S1633 = _S1691;
        _S1634 = _S1688;
        _S1635 = _S1679.differential_0.primal_0;
    }
    float2  _S1692 = _S1605 * _S1666;
    float2  _S1693 = _S1626 * _S1666;
    float2  _S1694 = _S1504;
    *&((&_S1694)->y) = k_6;
    *&((&_S1694)->x) = k_6;
    float2  _S1695 = _S1626 * _S1694;
    float2  _S1696 = _S1692 + _S1505 * _S1694;
    float _S1697 = _S1696.x;
    float _S1698 = _S1697 + _S1697;
    float _S1699 = _S1623 * _S1698;
    float _S1700 = _S1696.y + _S1696.y;
    float _S1701 = _S1623 * _S1700;
    float _S1702 = u_43 * _S1698 + v_43 * _S1700;
    float _S1703 = _S1503[int(3)] * _S1702;
    float _S1704 = _S1622 * _S1702;
    float _S1705 = _S1621 * _S1702;
    float _S1706 = _S1621 * _S1703;
    float _S1707 = _S1620 * _S1702;
    float _S1708 = _S1607 * _S1702 + r2_43 * _S1703;
    float _S1709 = _S1620 * _S1708;
    float _S1710 = _S1619 * _S1702;
    float _S1711 = _S1608 * _S1702 + r2_43 * _S1708;
    float _S1712 = _S1619 * _S1711;
    float _S1713 = _S1609 * _S1702 + r2_43 * _S1711;
    float _S1714 = v_43 * (_S1522 * _S1696.x);
    float _S1715 = _S1611 * _S1696.y;
    float _S1716 = _S1503[int(5)] * (_S1702 + u_43 * (2.0f * _S1696.x) + _S1612 * _S1696.x);
    float _S1717 = _S1503[int(6)] * _S1702;
    float _S1718 = _S1526 * _S1696.x;
    float _S1719 = _S1618 * _S1696.x;
    float _S1720 = v_43 * _S1718;
    float _S1721 = fy_19 * _S1718;
    float _S1722 = _S1613 * _S1696.y;
    float _S1723 = fy_19 * _S1696.y;
    float _S1724 = 2.0f * _S1696.y;
    float _S1725 = _S1617 * _S1724;
    float _S1726 = _S1617 * _S1696.y;
    float _S1727 = _S1702 + v_43 * _S1724 + _S1614 * _S1696.y;
    float _S1728 = _S1503[int(4)] * _S1727;
    float _S1729 = fy_19 * _S1727;
    float _S1730 = _S1503[int(7)] * _S1702;
    float _S1731 = fy_19 * _S1702;
    float2  _S1732 = _S1610 * _S1696;
    float2  _S1733 = _S1615 * _S1696;
    float2  _S1734 = _S1504;
    *&((&_S1734)->y) = _S1713;
    *&((&_S1734)->x) = _S1713;
    float2  _S1735 = _S1615 * _S1734;
    float _S1736 = _S1720 + _S1722 + _S1728 + _S1730;
    float _S1737 = _S1714 + _S1715 + _S1716 + _S1717;
    float2  _S1738 = _S1732 + _S1606 * _S1734;
    float2  _S1739 = _S1504;
    *&((&_S1739)->y) = _S1736;
    *&((&_S1739)->x) = _S1737;
    float2  _S1740 = _S1738 + _S1739;
    float _S1741 = _S1733.x + _S1733.y;
    float _S1742 = _S1710 + r2_43 * _S1741;
    float _S1743 = _S1707 + r2_43 * _S1742;
    float _S1744 = _S1705 + r2_43 * _S1743;
    float _S1745 = _S1706 + _S1709 + _S1712 + _S1609 * _S1741 + _S1608 * _S1742 + _S1607 * _S1743 + _S1503[int(3)] * _S1744;
    float _S1746 = v_43 * _S1745;
    float _S1747 = u_43 * _S1745;
    float2  _S1748 = _S1695 + _S1671.differential_0.primal_0;
    float2  _S1749 = _S1735 + make_float2 (_S1699 + _S1526 * _S1723 + _S1747 + _S1747, _S1701 + _S1721 + _S1725 + 2.0f * _S1726 + _S1746 + _S1746);
    float _S1750 = _S1719 + u_43 * _S1723;
    float _S1751 = _S1704 + r2_43 * _S1744;
    float2  _S1752 = _S1505 * _S1749;
    float _S1753 = _S1693.x + _S1693.y + _S1752.x + _S1752.y;
    float2  _S1754 = _S1605 * _S1749 + _S1748;
    if(_S1537)
    {
        float _S1755 = _S1507 * _S1631;
        float _S1756 = _S1753 / _S1600;
        float _S1757 = _S1536 * (0.3333333432674408f * - (_S1632 + _S1507 * _S1756));
        float _S1758 = _S1755 + _S1755 + _S1601 * - _S1756 + _S1635;
        k_6 = _S1757 + _S1757 + _S1634;
        _S1571 = _S1758;
        _S1600 = _S1633;
    }
    else
    {
        float _S1759 = _S1506 * _S1630;
        float _S1760 = _S1753 / _S1571;
        float _S1761 = _S1759 + _S1759 + _S1536 * - _S1760 + _S1633;
        k_6 = _S1506 * _S1760 + _S1634;
        _S1571 = _S1635;
        _S1600 = _S1761;
    }
    DiffPair_float_0 _S1762;
    (&_S1762)->primal_0 = _S1506;
    (&_S1762)->differential_0 = 0.0f;
    DiffPair_float_0 _S1763;
    (&_S1763)->primal_0 = _S1507;
    (&_S1763)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1762, &_S1763, k_6);
    float _S1764 = _S1763.differential_0 + _S1571;
    float _S1765 = _S1762.differential_0 + _S1600;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1766;
    (&_S1766)->primal_0 = _S1505;
    (&_S1766)->differential_0 = _S1504;
    s_bwd_length_impl_1(&_S1766, _S1765);
    float2  _S1767 = _S1766.differential_0 + _S1754;
    float3  _S1768 = make_float3 (_S1767.x, _S1767.y, _S1764);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1769;
    (&_S1769)->primal_0 = _S1505;
    (&_S1769)->differential_0 = _S1504;
    s_bwd_length_impl_1(&_S1769, 0.0f);
    float3  _S1770 = _S1768 + make_float3 (_S1769.differential_0.x, _S1769.differential_0.y, 0.0f);
    float2  _S1771 = _S1504;
    *&((&_S1771)->y) = _S1664.rows[int(0)].y;
    *&((&_S1771)->x) = _S1664.rows[int(0)].x;
    float2  _S1772 = _S1771;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1773 = { _S1504, _S1771 };
    DiffPair_0 _S1774;
    (&_S1774)->primal_0 = _S1593;
    (&_S1774)->differential_0 = _S1773;
    DiffPair_float_0 _S1775;
    (&_S1775)->primal_0 = _S1592;
    (&_S1775)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1776 = _S1596;
    s_bwd_prop_s_bwd_length_impl_0(&_S1774, &_S1775, &_S1776);
    DiffPair_0 _S1777 = _S1774;
    DiffPair_float_0 _S1778 = _S1775;
    DiffPair_float_0 _S1779 = { 0.0f, _S1664.rows[int(0)].z };
    DiffPair_float_0 _S1780 = { 0.0f, _S1775.differential_0 };
    DiffPair_1 _S1781;
    (&_S1781)->primal_0 = _S1589;
    (&_S1781)->differential_0 = _S1780;
    DiffPair_1 _S1782;
    (&_S1782)->primal_0 = _S1590;
    (&_S1782)->differential_0 = _S1779;
    DiffPair_float_0 _S1783;
    (&_S1783)->primal_0 = k_5;
    (&_S1783)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1781, &_S1782, &_S1783, &_s_diff_ctx_15->_S1253);
    DiffPair_1 _S1784 = _S1781;
    DiffPair_1 _S1785 = _S1782;
    DiffPair_float_0 _S1786 = _S1783;
    if(_S1537)
    {
        float _S1787 = _S1786.differential_0 + _S1786.differential_0;
        float _S1788 = _S1575 * _S1787;
        float _S1789 = - (0.3333333432674408f * (_S1536 * _S1787));
        float _S1790 = _S1577 * _S1664.rows[int(0)].z;
        float _S1791 = (_S1507 * _S1789 + - (_S1540 * _S1664.rows[int(0)].z)) / _S1578;
        float _S1792 = _S1570 * - _S1791;
        float _S1793 = _S1576 * _S1789 + _S1785.differential_0.primal_0;
        k_5 = _S1539 * _S1791;
        k_6 = 0.0f;
        _S1571 = _S1792;
        _S1572 = _S1790;
        _S1573 = _S1784.differential_0.primal_0;
        _S1574 = _S1788;
        _S1575 = _S1793;
    }
    else
    {
        float _S1794 = _S1573 * _S1778.differential_0;
        float _S1795 = (_S1506 * _S1786.differential_0 + - (_S1536 * _S1778.differential_0)) / _S1574;
        float _S1796 = _S1570 * - _S1795;
        float _S1797 = _S1572 * _S1786.differential_0 + _S1784.differential_0.primal_0;
        k_5 = _S1538 * _S1795;
        k_6 = _S1796;
        _S1571 = 0.0f;
        _S1572 = 0.0f;
        _S1573 = _S1797;
        _S1574 = _S1794;
        _S1575 = _S1785.differential_0.primal_0;
    }
    float2  _S1798 = _S1544 * _S1772;
    float2  _S1799 = _S1567 * _S1772;
    float2  _S1800 = _S1504;
    *&((&_S1800)->y) = k_5;
    *&((&_S1800)->x) = k_5;
    float2  _S1801 = _S1567 * _S1800;
    float2  _S1802 = _S1798 + _S1505 * _S1800;
    float _S1803 = _S1802.x;
    float _S1804 = _S1803 + _S1803;
    float _S1805 = _S1564 * _S1804;
    float _S1806 = _S1802.y + _S1802.y;
    float _S1807 = _S1564 * _S1806;
    float _S1808 = u_42 * _S1804 + v_42 * _S1806;
    float _S1809 = _S1503[int(3)] * _S1808;
    float _S1810 = _S1563 * _S1808;
    float _S1811 = _S1562 * _S1808;
    float _S1812 = _S1562 * _S1809;
    float _S1813 = _S1561 * _S1808;
    float _S1814 = _S1546 * _S1808 + r2_42 * _S1809;
    float _S1815 = _S1561 * _S1814;
    float _S1816 = _S1560 * _S1808;
    float _S1817 = _S1547 * _S1808 + r2_42 * _S1814;
    float _S1818 = _S1560 * _S1817;
    float _S1819 = _S1548 * _S1808 + r2_42 * _S1817;
    float _S1820 = _S1522 * _S1802.x;
    float _S1821 = _S1559 * _S1802.x;
    float _S1822 = v_42 * _S1820;
    float _S1823 = _S1554.x * _S1820;
    float _S1824 = _S1550 * _S1802.y;
    float _S1825 = _S1554.x * _S1802.y;
    float _S1826 = 2.0f * _S1802.x;
    float _S1827 = _S1558 * _S1826;
    float _S1828 = _S1558 * _S1802.x;
    float _S1829 = _S1808 + u_42 * _S1826 + _S1551 * _S1802.x;
    float _S1830 = _S1503[int(5)] * _S1829;
    float _S1831 = _S1554.x * _S1829;
    float _S1832 = _S1503[int(6)] * _S1808;
    float _S1833 = _S1554.x * _S1808;
    float _S1834 = _S1526 * _S1802.x;
    float _S1835 = _S1557 * _S1802.x;
    float _S1836 = v_42 * _S1834;
    float _S1837 = _S1554.y * _S1834;
    float _S1838 = _S1552 * _S1802.y;
    float _S1839 = _S1554.y * _S1802.y;
    float _S1840 = 2.0f * _S1802.y;
    float _S1841 = _S1556 * _S1840;
    float _S1842 = _S1556 * _S1802.y;
    float _S1843 = _S1808 + v_42 * _S1840 + _S1553 * _S1802.y;
    float _S1844 = _S1503[int(4)] * _S1843;
    float _S1845 = _S1554.y * _S1843;
    float _S1846 = _S1503[int(7)] * _S1808;
    float _S1847 = _S1554.y * _S1808;
    float2  _S1848 = _S1549 * _S1802;
    float2  _S1849 = _S1554 * _S1802;
    float2  _S1850 = _S1504;
    *&((&_S1850)->y) = _S1819;
    *&((&_S1850)->x) = _S1819;
    float2  _S1851 = _S1554 * _S1850;
    float _S1852 = _S1836 + _S1838 + _S1844 + _S1846;
    float _S1853 = _S1822 + _S1824 + _S1830 + _S1832;
    float2  _S1854 = _S1848 + _S1545 * _S1850;
    float2  _S1855 = _S1504;
    *&((&_S1855)->y) = _S1852;
    *&((&_S1855)->x) = _S1853;
    float2  _S1856 = _S1854 + _S1855;
    float _S1857 = fx_19 * _S1856.x;
    float _S1858 = fx_19 * _S1856.y;
    float _S1859 = _S1849.x + _S1849.y;
    float _S1860 = _S1816 + r2_42 * _S1859;
    float _S1861 = _S1813 + r2_42 * _S1860;
    float _S1862 = _S1811 + r2_42 * _S1861;
    float _S1863 = _S1812 + _S1815 + _S1818 + _S1548 * _S1859 + _S1547 * _S1860 + _S1546 * _S1861 + _S1503[int(3)] * _S1862;
    float _S1864 = v_42 * _S1863;
    float _S1865 = u_42 * _S1863;
    float2  _S1866 = _S1801 + _S1777.differential_0.primal_0;
    float _S1867 = _S1847 + _S1731;
    float _S1868 = _S1821 + u_42 * _S1825;
    float _S1869 = _S1503[int(8)] * _S1856.x + _S1503[int(9)] * _S1856.y + _S1856.x;
    float _S1870 = _S1861 + _S1743;
    float _S1871 = _S1860 + _S1742;
    float2  _S1872 = _S1851 + make_float2 (_S1805 + _S1827 + _S1526 * _S1839 + 2.0f * _S1828 + _S1522 * _S1825 + _S1865 + _S1865, _S1807 + _S1823 + _S1837 + _S1841 + 2.0f * _S1842 + _S1864 + _S1864);
    float _S1873 = _S1835 + u_42 * _S1839 + _S1750;
    float _S1874 = _S1810 + r2_42 * _S1862 + _S1751;
    float _S1875 = _S1862 + _S1744;
    float _S1876 = _S1845 + _S1729;
    float2  _S1877 = _S1505 * _S1872;
    float _S1878 = _S1799.x + _S1799.y + _S1877.x + _S1877.y;
    float2  _S1879 = _S1544 * _S1872 + _S1866;
    if(_S1537)
    {
        float _S1880 = _S1507 * _S1571;
        float _S1881 = _S1878 / _S1539;
        float _S1882 = _S1536 * (0.3333333432674408f * - (_S1572 + _S1507 * _S1881));
        float _S1883 = _S1880 + _S1880 + _S1540 * - _S1881 + _S1575;
        k_5 = _S1882 + _S1882 + _S1574;
        _S1538 = _S1883;
        _S1539 = _S1573;
    }
    else
    {
        float _S1884 = _S1506 * k_6;
        float _S1885 = _S1878 / _S1538;
        float _S1886 = _S1884 + _S1884 + _S1536 * - _S1885 + _S1573;
        k_5 = _S1506 * _S1885 + _S1574;
        _S1538 = _S1575;
        _S1539 = _S1886;
    }
    DiffPair_float_0 _S1887;
    (&_S1887)->primal_0 = _S1506;
    (&_S1887)->differential_0 = 0.0f;
    DiffPair_float_0 _S1888;
    (&_S1888)->primal_0 = _S1507;
    (&_S1888)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1887, &_S1888, k_5);
    float _S1889 = _S1888.differential_0 + _S1538;
    float _S1890 = _S1887.differential_0 + _S1539;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1891;
    (&_S1891)->primal_0 = _S1505;
    (&_S1891)->differential_0 = _S1504;
    s_bwd_length_impl_1(&_S1891, _S1890);
    float2  _S1892 = _S1891.differential_0 + _S1879;
    float3  _S1893 = make_float3 (_S1892.x, _S1892.y, _S1889);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1894;
    (&_S1894)->primal_0 = _S1505;
    (&_S1894)->differential_0 = _S1504;
    s_bwd_length_impl_1(&_S1894, 0.0f);
    float _S1895 = fx_19 * dpmean2d_1.x;
    float2  _S1896 = make_float2 (_S1895, fy_19 * dpmean2d_1.y) + make_float2 (_S1503[int(8)] * _S1895, _S1503[int(9)] * _S1895);
    float2  _S1897 = _S1517 * _S1896;
    float2  _S1898 = _S1521 * _S1896;
    float _S1899 = _S1503[int(4)] * _S1896.y;
    float _S1900 = v_41 * _S1896.y;
    float _S1901 = _S1503[int(5)] * _S1896.x;
    float _S1902 = v_41 * _S1896.x;
    float _S1903 = _S1897.x + _S1897.y;
    float _S1904 = r2_41 * _S1903;
    float _S1905 = r2_41 * _S1904;
    float _S1906 = r2_41 * _S1905;
    float _S1907 = _S1503[int(7)] * _S1896.y + _S1899 + _S1503[int(6)] * _S1896.x + _S1901 + _S1520 * _S1903 + _S1519 * _S1904 + _S1518 * _S1905 + _S1503[int(3)] * _S1906;
    float _S1908 = v_41 * _S1907;
    float _S1909 = u_41 * _S1907;
    float _S1910 = _S1528 * _S1899 + 2.0f * (v_41 * _S1899) + _S1527 * _S1896.y + _S1523 * _S1896.x + _S1908 + _S1908;
    float _S1911 = _S1526 * _S1900 + _S1524 * _S1901 + 2.0f * (u_41 * _S1901) + _S1522 * _S1902 + _S1909 + _S1909;
    float _S1912 = _S1532 * _S1895 + _S1858;
    float _S1913 = _S1531 * _S1895 + _S1857;
    float _S1914 = r2_41 * _S1896.y + _S1867;
    float _S1915 = r2_41 * _S1896.x + _S1833;
    float _S1916 = 2.0f * (u_41 * _S1900 + _S1873) + _S1525 * _S1896.x + _S1831;
    float _S1917 = _S1529 * _S1896.y + 2.0f * (u_41 * _S1902 + _S1868) + _S1876;
    float _S1918 = r2_41 * _S1906 + _S1874;
    float _S1919 = _S1906 + _S1875;
    float _S1920 = _S1905 + _S1870;
    float _S1921 = _S1904 + _S1871;
    float3  _S1922 = _S1893 + make_float3 (_S1894.differential_0.x, _S1894.differential_0.y, 0.0f) + _S1770;
    float4  _S1923 = make_float4 (_S1534 * dpmean2d_1.x + _S1869, _S1535 * dpmean2d_1.y + _S1740.y, dpmean2d_1.x, dpmean2d_1.y);
    FixedArray<float, 10>  _S1924;
    _S1924[int(0)] = 0.0f;
    _S1924[int(1)] = 0.0f;
    _S1924[int(2)] = 0.0f;
    _S1924[int(3)] = 0.0f;
    _S1924[int(4)] = 0.0f;
    _S1924[int(5)] = 0.0f;
    _S1924[int(6)] = 0.0f;
    _S1924[int(7)] = 0.0f;
    _S1924[int(8)] = 0.0f;
    _S1924[int(9)] = 0.0f;
    _S1924[int(9)] = _S1912;
    _S1924[int(8)] = _S1913;
    _S1924[int(7)] = _S1914;
    _S1924[int(6)] = _S1915;
    _S1924[int(5)] = _S1916;
    _S1924[int(4)] = _S1917;
    _S1924[int(3)] = _S1918;
    _S1924[int(2)] = _S1919;
    _S1924[int(1)] = _S1920;
    _S1924[int(0)] = _S1921;
    FixedArray<float, 10>  _S1925 = {
        _S1924[int(0)], _S1924[int(1)], _S1924[int(2)], _S1924[int(3)], _S1924[int(4)], _S1924[int(5)], _S1924[int(6)], _S1924[int(7)], _S1924[int(8)], _S1924[int(9)]
    };
    float2  _S1926 = _S1898 + make_float2 (_S1911, _S1910);
    float2  _S1927 = _S1505 * _S1926;
    float2  _S1928 = _S1516 * _S1926;
    float _S1929 = _S1927.x + _S1927.y;
    if(_S1509)
    {
        float _S1930 = _S1929 / _S1511;
        float _S1931 = _S1512 * - _S1930;
        float _S1932 = _S1508 * (0.3333333432674408f * - (_S1507 * _S1930));
        k_5 = _S1932 + _S1932;
        _S1510 = _S1931;
        _S1511 = 0.0f;
    }
    else
    {
        float _S1933 = _S1929 / _S1510;
        float _S1934 = _S1508 * - _S1933;
        k_5 = _S1506 * _S1933;
        _S1510 = 0.0f;
        _S1511 = _S1934;
    }
    DiffPair_float_0 _S1935;
    (&_S1935)->primal_0 = _S1506;
    (&_S1935)->differential_0 = 0.0f;
    DiffPair_float_0 _S1936;
    (&_S1936)->primal_0 = _S1507;
    (&_S1936)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1935, &_S1936, k_5);
    float _S1937 = _S1936.differential_0 + _S1510;
    float _S1938 = _S1935.differential_0 + _S1511;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1939;
    (&_S1939)->primal_0 = _S1505;
    (&_S1939)->differential_0 = _S1504;
    s_bwd_length_impl_1(&_S1939, _S1938);
    float2  _S1940 = _S1939.differential_0 + _S1928;
    dpdist_coeffs_1->primal_0 = dpdist_coeffs_1->primal_0;
    dpdist_coeffs_1->differential_0 = _S1925;
    dpintrins_1->primal_0 = (*dpintrins_1).primal_0;
    dpintrins_1->differential_0 = _S1923;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1663.differential_0;
    float3  _S1941 = _S1922 + make_float3 (_S1940.x, _S1940.y, _S1937);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1941;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_20, float fy_20, float cx_15, float cy_15, FixedArray<float, 10>  * dist_coeffs_26, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1942 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1943 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1944 = { _S1943, _S1943 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1945 = { _S1943, _S1943, _S1944, _S1943, _S1943, _S1944 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1946;
    (&_S1946)->_S1257 = _S1942;
    (&_S1946)->_S1258 = _S1945;
    float3  mean_c_11 = s_primal_ctx_mul_1(R_15, mean_11) + t_14;
    float4  intrins_11 = make_float4 (fx_20, fy_20, cx_15, cy_15);
    float3  _S1947 = s_primal_ctx_exp_0(scale_13);
    float _S1948 = quat_14.y;
    float x2_14 = _S1948 * _S1948;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_25 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S1949 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_25), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_25), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1947.x, 0.0f, 0.0f, 0.0f, _S1947.y, 0.0f, 0.0f, 0.0f, _S1947.z);
    Matrix<float, 3, 3>  _S1950 = s_primal_ctx_mul_2(_S1949, S_1);
    Matrix<float, 3, 3>  _S1951 = transpose_0(_S1950);
    Matrix<float, 3, 3>  _S1952 = s_primal_ctx_mul_2(_S1950, _S1951);
    Matrix<float, 3, 3>  _S1953 = s_primal_ctx_mul_2(R_15, _S1952);
    Matrix<float, 3, 3>  _S1954 = transpose_0(R_15);
    Matrix<float, 3, 3>  _S1955 = s_primal_ctx_mul_2(_S1953, _S1954);
    Matrix<float, 2, 2>  _S1956 = _S1942;
    float2  _S1957 = make_float2 (0.0f);
    float2  _S1958 = _S1957;
    bool _S1959 = s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1955, intrins_11, dist_coeffs_26, &_S1956, &_S1958, &(&_S1946)->_S1258);
    (&_S1946)->_S1257 = _S1956;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1960 = _S1946;
    float _S1961 = _S1946._S1257.rows[int(0)].y * _S1946._S1257.rows[int(1)].x;
    float det_orig_12 = _S1946._S1257.rows[int(0)].x * _S1946._S1257.rows[int(1)].y - _S1961;
    float _S1962 = _S1946._S1257.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1963 = _S1946._S1257;
    *&(((&_S1963)->rows + (int(0)))->x) = _S1962;
    float _S1964 = _S1946._S1257.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1963)->rows + (int(1)))->y) = _S1964;
    Matrix<float, 2, 2>  _S1965 = _S1963;
    Matrix<float, 2, 2>  _S1966 = _S1963;
    float det_blur_7 = _S1962 * _S1964 - _S1961;
    float _S1967 = det_orig_12 / det_blur_7;
    float _S1968 = det_blur_7 * det_blur_7;
    float _S1969 = s_primal_ctx_max_0(0.0f, _S1967);
    float _S1970 = s_primal_ctx_sqrt_0(_S1969);
    float invdet_8 = 1.0f / det_blur_7;
    float _S1971 = - _S1946._S1257.rows[int(0)].y;
    float _S1972 = - _S1946._S1257.rows[int(1)].x;
    float _S1973 = - in_opacity_11;
    float _S1974 = 1.0f + s_primal_ctx_exp_1(_S1973);
    float _S1975 = 1.0f / _S1974;
    float _S1976 = _S1974 * _S1974;
    float _S1977;
    if(antialiased_11)
    {
        _S1977 = _S1975 * _S1970;
    }
    else
    {
        _S1977 = _S1975;
    }
    float _S1978 = _S1977 / 0.00392156885936856f;
    float _S1979 = 2.0f * s_primal_ctx_log_0(_S1978);
    float _S1980 = s_primal_ctx_sqrt_0(_S1979);
    float _S1981 = _S1965.rows[int(0)].x;
    float _S1982 = _S1966.rows[int(1)].y;
    float _S1983 = s_primal_ctx_dot_0(mean_c_11, mean_c_11) + 9.99999997475242708e-07f;
    float3  _S1984 = mean_11 - - s_primal_ctx_mul_1(_S1954, t_14);
    float _S1985 = _S1984.x;
    float _S1986 = _S1984.y;
    float _S1987 = _S1984.z;
    float _S1988 = _S1985 * _S1985 + _S1986 * _S1986 + _S1987 * _S1987;
    float _S1989 = s_primal_ctx_sqrt_0(_S1988);
    float x_36 = _S1985 / _S1989;
    float3  _S1990 = make_float3 (x_36);
    float _S1991 = _S1989 * _S1989;
    float y_14 = _S1986 / _S1989;
    float z_11 = _S1987 / _S1989;
    float3  _S1992 = make_float3 (z_11);
    float _S1993 = - y_14;
    float3  _S1994 = make_float3 (_S1993);
    float z2_26 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_36 * x_36 - y_14 * y_14;
    float _S1995 = 2.0f * x_36;
    float fS1_11 = _S1995 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_26 - 0.31539157032966614f;
    float3  _S1996 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_36;
    float3  _S1997 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S1998 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S1999 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S2000 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_26 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S2001 = 1.86588168144226074f * z2_26 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S2001;
    float3  _S2002 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_36;
    float3  _S2003 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S2004 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S2005 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S2006 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_36 * fC1_11 - y_14 * fS1_11);
    float3  _S2007 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_36 * fS1_11 + y_14 * fC1_11);
    float3  _S2008 = make_float3 (pSH9_1);
    float3  _S2009 = make_float3 (0.0f);
    float3  _S2010 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2011;
    (&_S2011)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1993) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_36) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S2011)->differential_0 = _S2010;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2012;
    (&_S2012)->primal_0 = _S2009;
    (&_S2012)->differential_0 = _S2010;
    s_bwd_prop_max_0(&_S2011, &_S2012, v_rgb_1);
    float3  _S2013 = _S2007 * _S2011.differential_0;
    float3  _S2014 = (*sh_coeffs_11)[int(15)] * _S2011.differential_0;
    float3  _S2015 = _S2005 * _S2011.differential_0;
    float3  _S2016 = (*sh_coeffs_11)[int(14)] * _S2011.differential_0;
    float3  _S2017 = _S2003 * _S2011.differential_0;
    float3  _S2018 = (*sh_coeffs_11)[int(13)] * _S2011.differential_0;
    float3  _S2019 = _S2002 * _S2011.differential_0;
    float3  _S2020 = (*sh_coeffs_11)[int(12)] * _S2011.differential_0;
    float3  _S2021 = _S2004 * _S2011.differential_0;
    float3  _S2022 = (*sh_coeffs_11)[int(11)] * _S2011.differential_0;
    float3  _S2023 = _S2006 * _S2011.differential_0;
    float3  _S2024 = (*sh_coeffs_11)[int(10)] * _S2011.differential_0;
    float3  _S2025 = _S2008 * _S2011.differential_0;
    float3  _S2026 = (*sh_coeffs_11)[int(9)] * _S2011.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S2026.x + _S2026.y + _S2026.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S2014.x + _S2014.y + _S2014.z);
    float _S2027 = _S2024.x + _S2024.y + _S2024.z;
    float _S2028 = _S2016.x + _S2016.y + _S2016.z;
    float _S2029 = _S2022.x + _S2022.y + _S2022.z;
    float _S2030 = _S2018.x + _S2018.y + _S2018.z;
    float _S2031 = _S2020.x + _S2020.y + _S2020.z;
    float _S2032 = - s_diff_fC2_T_1;
    float3  _S2033 = _S1999 * _S2011.differential_0;
    float3  _S2034 = (*sh_coeffs_11)[int(8)] * _S2011.differential_0;
    float3  _S2035 = _S1997 * _S2011.differential_0;
    float3  _S2036 = (*sh_coeffs_11)[int(7)] * _S2011.differential_0;
    float3  _S2037 = _S1996 * _S2011.differential_0;
    float3  _S2038 = (*sh_coeffs_11)[int(6)] * _S2011.differential_0;
    float3  _S2039 = _S1998 * _S2011.differential_0;
    float3  _S2040 = (*sh_coeffs_11)[int(5)] * _S2011.differential_0;
    float3  _S2041 = _S2000 * _S2011.differential_0;
    float3  _S2042 = (*sh_coeffs_11)[int(4)] * _S2011.differential_0;
    float _S2043 = _S2040.x + _S2040.y + _S2040.z;
    float _S2044 = _S2036.x + _S2036.y + _S2036.z;
    float _S2045 = fTmp1B_11 * _S2027 + x_36 * s_diff_fS2_T_1 + y_14 * _S2032 + 0.54627424478530884f * (_S2042.x + _S2042.y + _S2042.z);
    float _S2046 = fTmp1B_11 * _S2028 + y_14 * s_diff_fS2_T_1 + x_36 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S2034.x + _S2034.y + _S2034.z);
    float _S2047 = y_14 * - _S2046;
    float _S2048 = x_36 * _S2046;
    float _S2049 = z_11 * (1.86588168144226074f * (z_11 * _S2031) + -2.28522896766662598f * (y_14 * _S2029 + x_36 * _S2030) + 0.94617468118667603f * (_S2038.x + _S2038.y + _S2038.z));
    float3  _S2050 = make_float3 (0.48860251903533936f) * _S2011.differential_0;
    float3  _S2051 = - _S2050;
    float3  _S2052 = _S1990 * _S2051;
    float3  _S2053 = (*sh_coeffs_11)[int(3)] * _S2051;
    float3  _S2054 = _S1992 * _S2050;
    float3  _S2055 = (*sh_coeffs_11)[int(2)] * _S2050;
    float3  _S2056 = _S1994 * _S2050;
    float3  _S2057 = (*sh_coeffs_11)[int(1)] * _S2050;
    float _S2058 = (_S2001 * _S2031 + 1.44530570507049561f * (fS1_11 * _S2027 + fC1_11 * _S2028) + -1.09254848957061768f * (y_14 * _S2043 + x_36 * _S2044) + _S2049 + _S2049 + _S2055.x + _S2055.y + _S2055.z) / _S1991;
    float _S2059 = _S1989 * _S2058;
    float _S2060 = (fTmp0C_11 * _S2029 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S2032 + fTmp0B_11 * _S2043 + _S1995 * _S2045 + _S2047 + _S2047 + - (_S2057.x + _S2057.y + _S2057.z)) / _S1991;
    float _S2061 = _S1989 * _S2060;
    float _S2062 = (fTmp0C_11 * _S2030 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S2044 + 2.0f * (y_14 * _S2045) + _S2048 + _S2048 + _S2053.x + _S2053.y + _S2053.z) / _S1991;
    float _S2063 = _S1989 * _S2062;
    float _S2064 = _S1987 * - _S2058 + _S1986 * - _S2060 + _S1985 * - _S2062;
    DiffPair_float_0 _S2065;
    (&_S2065)->primal_0 = _S1988;
    (&_S2065)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2065, _S2064);
    float _S2066 = _S1987 * _S2065.differential_0;
    float _S2067 = _S1986 * _S2065.differential_0;
    float _S2068 = _S1985 * _S2065.differential_0;
    float3  _S2069 = make_float3 (0.282094806432724f) * _S2011.differential_0;
    float3  _S2070 = make_float3 (_S2063 + _S2068 + _S2068, _S2061 + _S2067 + _S2067, _S2059 + _S2066 + _S2066);
    float3  _S2071 = - - _S2070;
    Matrix<float, 3, 3>  _S2072 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2073;
    (&_S2073)->primal_0 = _S1954;
    (&_S2073)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2074;
    (&_S2074)->primal_0 = t_14;
    (&_S2074)->differential_0 = _S2010;
    s_bwd_prop_mul_1(&_S2073, &_S2074, _S2071);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2075 = _S2073;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2076 = _S2074;
    float2  _S2077 = _S1957;
    *&((&_S2077)->y) = v_conic_1.z;
    float2  _S2078 = _S1957;
    *&((&_S2078)->y) = v_conic_1.y;
    *&((&_S2078)->x) = v_conic_1.x;
    float _S2079 = 0.5f * v_depth_1;
    DiffPair_float_0 _S2080;
    (&_S2080)->primal_0 = _S1983;
    (&_S2080)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2080, _S2079);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2081;
    (&_S2081)->primal_0 = mean_c_11;
    (&_S2081)->differential_0 = _S2010;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2082;
    (&_S2082)->primal_0 = mean_c_11;
    (&_S2082)->differential_0 = _S2010;
    s_bwd_prop_dot_0(&_S2081, &_S2082, _S2080.differential_0);
    DiffPair_float_0 _S2083;
    (&_S2083)->primal_0 = _S1982;
    (&_S2083)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2083, 0.0f);
    DiffPair_float_0 _S2084;
    (&_S2084)->primal_0 = _S1981;
    (&_S2084)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2084, 0.0f);
    DiffPair_float_0 _S2085;
    (&_S2085)->primal_0 = 3.32999992370605469f;
    (&_S2085)->differential_0 = 0.0f;
    DiffPair_float_0 _S2086;
    (&_S2086)->primal_0 = _S1980;
    (&_S2086)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2085, &_S2086, 0.0f);
    DiffPair_float_0 _S2087;
    (&_S2087)->primal_0 = _S1979;
    (&_S2087)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2087, _S2086.differential_0);
    float _S2088 = 2.0f * _S2087.differential_0;
    DiffPair_float_0 _S2089;
    (&_S2089)->primal_0 = _S1978;
    (&_S2089)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2089, _S2088);
    float2  _S2090 = make_float2 (_S2084.differential_0, 0.0f);
    float _S2091 = v_opacity_1 + 254.9999847412109375f * _S2089.differential_0;
    Matrix<float, 2, 2>  _S2092 = _S1942;
    _S2092[int(1)] = _S2077;
    _S2092[int(0)] = _S2078;
    Matrix<float, 2, 2>  _S2093 = _S2092;
    FixedArray<float3 , 16>  _S2094;
    _S2094[int(0)] = _S2010;
    _S2094[int(1)] = _S2010;
    _S2094[int(2)] = _S2010;
    _S2094[int(3)] = _S2010;
    _S2094[int(4)] = _S2010;
    _S2094[int(5)] = _S2010;
    _S2094[int(6)] = _S2010;
    _S2094[int(7)] = _S2010;
    _S2094[int(8)] = _S2010;
    _S2094[int(9)] = _S2010;
    _S2094[int(10)] = _S2010;
    _S2094[int(11)] = _S2010;
    _S2094[int(12)] = _S2010;
    _S2094[int(13)] = _S2010;
    _S2094[int(14)] = _S2010;
    _S2094[int(15)] = _S2010;
    _S2094[int(7)] = _S2035;
    _S2094[int(0)] = _S2069;
    _S2094[int(1)] = _S2056;
    _S2094[int(2)] = _S2054;
    _S2094[int(3)] = _S2052;
    _S2094[int(4)] = _S2041;
    _S2094[int(5)] = _S2039;
    _S2094[int(6)] = _S2037;
    _S2094[int(15)] = _S2013;
    _S2094[int(8)] = _S2033;
    _S2094[int(9)] = _S2025;
    _S2094[int(10)] = _S2023;
    _S2094[int(11)] = _S2021;
    _S2094[int(12)] = _S2019;
    _S2094[int(13)] = _S2017;
    _S2094[int(14)] = _S2015;
    float3  _S2095 = _S2094[int(0)];
    float3  _S2096 = _S2094[int(1)];
    float3  _S2097 = _S2094[int(2)];
    float3  _S2098 = _S2094[int(3)];
    float3  _S2099 = _S2094[int(4)];
    float3  _S2100 = _S2094[int(5)];
    float3  _S2101 = _S2094[int(6)];
    float3  _S2102 = _S2094[int(7)];
    float3  _S2103 = _S2094[int(8)];
    float3  _S2104 = _S2094[int(9)];
    float3  _S2105 = _S2094[int(10)];
    float3  _S2106 = _S2094[int(11)];
    float3  _S2107 = _S2094[int(12)];
    float3  _S2108 = _S2094[int(13)];
    float3  _S2109 = _S2094[int(14)];
    float3  _S2110 = _S2094[int(15)];
    float3  _S2111 = _S2082.differential_0 + _S2081.differential_0;
    float2  _S2112 = make_float2 (0.0f, _S2083.differential_0);
    float _S2113;
    if(antialiased_11)
    {
        float _S2114 = _S1975 * _S2091;
        _S1977 = _S1970 * _S2091;
        _S2113 = _S2114;
    }
    else
    {
        _S1977 = _S2091;
        _S2113 = 0.0f;
    }
    float _S2115 = - (_S1977 / _S1976);
    DiffPair_float_0 _S2116;
    (&_S2116)->primal_0 = _S1973;
    (&_S2116)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2116, _S2115);
    float _S2117 = - _S2116.differential_0;
    float _S2118 = invdet_8 * _S2093.rows[int(1)].y;
    float _S2119 = - (invdet_8 * _S2093.rows[int(1)].x);
    float _S2120 = - (invdet_8 * _S2093.rows[int(0)].y);
    float _S2121 = invdet_8 * _S2093.rows[int(0)].x;
    float _S2122 = - ((_S1962 * _S2093.rows[int(1)].y + _S1972 * _S2093.rows[int(1)].x + _S1971 * _S2093.rows[int(0)].y + _S1964 * _S2093.rows[int(0)].x) / _S1968);
    DiffPair_float_0 _S2123;
    (&_S2123)->primal_0 = _S1969;
    (&_S2123)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2123, _S2113);
    DiffPair_float_0 _S2124;
    (&_S2124)->primal_0 = 0.0f;
    (&_S2124)->differential_0 = 0.0f;
    DiffPair_float_0 _S2125;
    (&_S2125)->primal_0 = _S1967;
    (&_S2125)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2124, &_S2125, _S2123.differential_0);
    float _S2126 = _S2125.differential_0 / _S1968;
    float s_diff_det_orig_T_1 = det_blur_7 * _S2126;
    float _S2127 = _S2122 + det_orig_12 * - _S2126;
    float _S2128 = - _S2127;
    float _S2129 = _S1962 * _S2127;
    float _S2130 = _S1964 * _S2127;
    Matrix<float, 2, 2>  _S2131 = _S1942;
    _S2131[int(1)] = _S2112;
    _S2131[int(0)] = _S2090;
    _S1963 = _S2131;
    *&(((&_S1963)->rows + (int(1)))->y) = 0.0f;
    float _S2132 = _S2121 + _S2129 + _S2131.rows[int(1)].y;
    *&(((&_S1963)->rows + (int(0)))->x) = 0.0f;
    float _S2133 = _S2118 + _S2130 + _S2131.rows[int(0)].x;
    float _S2134 = _S2128 + - s_diff_det_orig_T_1;
    float _S2135 = _S2119 + _S1960._S1257.rows[int(0)].y * _S2134;
    float _S2136 = _S2120 + _S1960._S1257.rows[int(1)].x * _S2134;
    float _S2137 = _S1960._S1257.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S2138 = _S2132 + _S1960._S1257.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S2139 = _S1957;
    *&((&_S2139)->x) = _S2135;
    *&((&_S2139)->y) = _S2138;
    float _S2140 = _S2133 + _S2137;
    float2  _S2141 = _S1957;
    *&((&_S2141)->y) = _S2136;
    *&((&_S2141)->x) = _S2140;
    Matrix<float, 2, 2>  _S2142 = _S1942;
    _S2142[int(1)] = _S2139;
    _S2142[int(0)] = _S2141;
    Matrix<float, 2, 2>  _S2143 = _S1963 + _S2142;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2144;
    (&_S2144)->primal_0 = mean_c_11;
    (&_S2144)->differential_0 = _S2010;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2145;
    (&_S2145)->primal_0 = _S1955;
    (&_S2145)->differential_0 = _S2072;
    float4  _S2146 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2147;
    (&_S2147)->primal_0 = intrins_11;
    (&_S2147)->differential_0 = _S2146;
    FixedArray<float, 10>  _S2148 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S2149;
    (&_S2149)->primal_0 = *dist_coeffs_26;
    (&_S2149)->differential_0 = _S2148;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2150 = _S1960._S1258;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2144, &_S2145, &_S2147, &_S2149, _S2143, v_mean2d_1, &_S2150);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2151;
    (&_S2151)->primal_0 = _S1953;
    (&_S2151)->differential_0 = _S2072;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2152;
    (&_S2152)->primal_0 = _S1954;
    (&_S2152)->differential_0 = _S2072;
    s_bwd_prop_mul_4(&_S2151, &_S2152, _S2145.differential_0);
    Matrix<float, 3, 3>  _S2153 = transpose_0(_S2152.differential_0 + _S2075.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2154;
    (&_S2154)->primal_0 = R_15;
    (&_S2154)->differential_0 = _S2072;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2155;
    (&_S2155)->primal_0 = _S1952;
    (&_S2155)->differential_0 = _S2072;
    s_bwd_prop_mul_4(&_S2154, &_S2155, _S2151.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2156;
    (&_S2156)->primal_0 = _S1950;
    (&_S2156)->differential_0 = _S2072;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2157;
    (&_S2157)->primal_0 = _S1951;
    (&_S2157)->differential_0 = _S2072;
    s_bwd_prop_mul_4(&_S2156, &_S2157, _S2155.differential_0);
    Matrix<float, 3, 3>  _S2158 = _S2156.differential_0 + transpose_0(_S2157.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2159;
    (&_S2159)->primal_0 = _S1949;
    (&_S2159)->differential_0 = _S2072;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2160;
    (&_S2160)->primal_0 = S_1;
    (&_S2160)->differential_0 = _S2072;
    s_bwd_prop_mul_4(&_S2159, &_S2160, _S2158);
    Matrix<float, 3, 3>  _S2161 = transpose_0(_S2159.differential_0);
    float _S2162 = 2.0f * - _S2161.rows[int(2)].z;
    float _S2163 = 2.0f * _S2161.rows[int(2)].y;
    float _S2164 = 2.0f * _S2161.rows[int(2)].x;
    float _S2165 = 2.0f * _S2161.rows[int(1)].z;
    float _S2166 = 2.0f * - _S2161.rows[int(1)].y;
    float _S2167 = 2.0f * _S2161.rows[int(1)].x;
    float _S2168 = 2.0f * _S2161.rows[int(0)].z;
    float _S2169 = 2.0f * _S2161.rows[int(0)].y;
    float _S2170 = 2.0f * - _S2161.rows[int(0)].x;
    float _S2171 = - _S2167 + _S2169;
    float _S2172 = _S2164 + - _S2168;
    float _S2173 = - _S2163 + _S2165;
    float _S2174 = _S2163 + _S2165;
    float _S2175 = _S2164 + _S2168;
    float _S2176 = _S2167 + _S2169;
    float _S2177 = quat_14.w * (_S2166 + _S2170);
    float _S2178 = quat_14.z * (_S2162 + _S2170);
    float _S2179 = quat_14.y * (_S2162 + _S2166);
    float _S2180 = quat_14.x * _S2171 + quat_14.z * _S2174 + quat_14.y * _S2175 + _S2177 + _S2177;
    float _S2181 = quat_14.x * _S2172 + quat_14.w * _S2174 + quat_14.y * _S2176 + _S2178 + _S2178;
    float _S2182 = quat_14.x * _S2173 + quat_14.w * _S2175 + quat_14.z * _S2176 + _S2179 + _S2179;
    float _S2183 = quat_14.w * _S2171 + quat_14.z * _S2172 + quat_14.y * _S2173;
    float3  _S2184 = _S2010;
    *&((&_S2184)->z) = _S2160.differential_0.rows[int(2)].z;
    *&((&_S2184)->y) = _S2160.differential_0.rows[int(1)].y;
    *&((&_S2184)->x) = _S2160.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2185;
    (&_S2185)->primal_0 = scale_13;
    (&_S2185)->differential_0 = _S2010;
    s_bwd_prop_exp_1(&_S2185, _S2184);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2186 = _S2185;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2187;
    (&_S2187)->primal_0 = mean_c_11;
    (&_S2187)->differential_0 = _S2010;
    s_bwd_length_impl_0(&_S2187, 0.0f);
    float3  _S2188 = _S2144.differential_0 + _S2187.differential_0 + _S2111;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2189;
    (&_S2189)->primal_0 = R_15;
    (&_S2189)->differential_0 = _S2072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2190;
    (&_S2190)->primal_0 = mean_11;
    (&_S2190)->differential_0 = _S2010;
    s_bwd_prop_mul_1(&_S2189, &_S2190, _S2188);
    float3  _S2191 = _S2188 + _S2076.differential_0;
    Matrix<float, 3, 3>  _S2192 = _S2153 + _S2154.differential_0 + _S2189.differential_0;
    float4  _S2193 = _S2146;
    *&((&_S2193)->w) = _S2180;
    *&((&_S2193)->z) = _S2181;
    *&((&_S2193)->y) = _S2182;
    *&((&_S2193)->x) = _S2183;
    float4  _S2194 = _S2193;
    float3  _S2195 = _S2190.differential_0 + _S2070;
    *v_mean_1 = _S2195;
    *v_quat_1 = _S2194;
    *v_scale_1 = _S2186.differential_0;
    *v_in_opacity_1 = _S2117;
    (*v_sh_coeffs_1)[int(0)] = _S2095;
    (*v_sh_coeffs_1)[int(1)] = _S2096;
    (*v_sh_coeffs_1)[int(2)] = _S2097;
    (*v_sh_coeffs_1)[int(3)] = _S2098;
    (*v_sh_coeffs_1)[int(4)] = _S2099;
    (*v_sh_coeffs_1)[int(5)] = _S2100;
    (*v_sh_coeffs_1)[int(6)] = _S2101;
    (*v_sh_coeffs_1)[int(7)] = _S2102;
    (*v_sh_coeffs_1)[int(8)] = _S2103;
    (*v_sh_coeffs_1)[int(9)] = _S2104;
    (*v_sh_coeffs_1)[int(10)] = _S2105;
    (*v_sh_coeffs_1)[int(11)] = _S2106;
    (*v_sh_coeffs_1)[int(12)] = _S2107;
    (*v_sh_coeffs_1)[int(13)] = _S2108;
    (*v_sh_coeffs_1)[int(14)] = _S2109;
    (*v_sh_coeffs_1)[int(15)] = _S2110;
    *v_R_2 = _S2192;
    *v_t_2 = _S2191;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_21, float fy_21, float cx_16, float cy_16, FixedArray<float, 10>  * dist_coeffs_27, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_1(R_16, mean_12) + t_15;
    float3  _S2196 = s_primal_ctx_exp_0(scale_14);
    float _S2197 = quat_15.y;
    float x2_15 = _S2197 * _S2197;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_27 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S2198 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_27), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_27), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S2196.x, 0.0f, 0.0f, 0.0f, _S2196.y, 0.0f, 0.0f, 0.0f, _S2196.z);
    Matrix<float, 3, 3>  _S2199 = s_primal_ctx_mul_2(_S2198, S_2);
    Matrix<float, 3, 3>  _S2200 = transpose_0(_S2199);
    Matrix<float, 3, 3>  _S2201 = s_primal_ctx_mul_2(_S2199, _S2200);
    Matrix<float, 3, 3>  _S2202 = s_primal_ctx_mul_2(R_16, _S2201);
    Matrix<float, 3, 3>  _S2203 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S2204 = s_primal_ctx_mul_2(_S2202, _S2203);
    Matrix<float, 2, 3>  J_13 = makeMatrix<float, 2, 3> (fx_21, 0.0f, 0.0f, 0.0f, fy_21, 0.0f);
    Matrix<float, 2, 3>  _S2205 = s_primal_ctx_mul_3(J_13, _S2204);
    Matrix<float, 3, 2>  _S2206 = transpose_1(J_13);
    Matrix<float, 2, 2>  _S2207 = s_primal_ctx_mul_4(_S2205, _S2206);
    float _S2208 = _S2207.rows[int(0)].y * _S2207.rows[int(1)].x;
    float det_orig_13 = _S2207.rows[int(0)].x * _S2207.rows[int(1)].y - _S2208;
    float _S2209 = _S2207.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2210 = _S2207;
    *&(((&_S2210)->rows + (int(0)))->x) = _S2209;
    float _S2211 = _S2207.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2210)->rows + (int(1)))->y) = _S2211;
    Matrix<float, 2, 2>  _S2212 = _S2210;
    Matrix<float, 2, 2>  _S2213 = _S2210;
    float det_blur_8 = _S2209 * _S2211 - _S2208;
    float _S2214 = det_orig_13 / det_blur_8;
    float _S2215 = det_blur_8 * det_blur_8;
    float _S2216 = s_primal_ctx_max_0(0.0f, _S2214);
    float _S2217 = s_primal_ctx_sqrt_0(_S2216);
    float invdet_9 = 1.0f / det_blur_8;
    float _S2218 = - _S2207.rows[int(0)].y;
    float _S2219 = - _S2207.rows[int(1)].x;
    float _S2220 = - in_opacity_12;
    float _S2221 = 1.0f + s_primal_ctx_exp_1(_S2220);
    float _S2222 = 1.0f / _S2221;
    float _S2223 = _S2221 * _S2221;
    float _S2224;
    if(antialiased_12)
    {
        _S2224 = _S2222 * _S2217;
    }
    else
    {
        _S2224 = _S2222;
    }
    float _S2225 = _S2224 / 0.00392156885936856f;
    float _S2226 = 2.0f * s_primal_ctx_log_0(_S2225);
    float _S2227 = s_primal_ctx_sqrt_0(_S2226);
    float _S2228 = _S2212.rows[int(0)].x;
    float _S2229 = _S2213.rows[int(1)].y;
    float _S2230 = s_primal_ctx_dot_0(mean_c_12, mean_c_12) + 9.99999997475242708e-07f;
    float3  _S2231 = mean_12 - - s_primal_ctx_mul_1(_S2203, t_15);
    float _S2232 = _S2231.x;
    float _S2233 = _S2231.y;
    float _S2234 = _S2231.z;
    float _S2235 = _S2232 * _S2232 + _S2233 * _S2233 + _S2234 * _S2234;
    float _S2236 = s_primal_ctx_sqrt_0(_S2235);
    float x_37 = _S2232 / _S2236;
    float3  _S2237 = make_float3 (x_37);
    float _S2238 = _S2236 * _S2236;
    float y_15 = _S2233 / _S2236;
    float z_12 = _S2234 / _S2236;
    float3  _S2239 = make_float3 (z_12);
    float _S2240 = - y_15;
    float3  _S2241 = make_float3 (_S2240);
    float z2_28 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_37 * x_37 - y_15 * y_15;
    float _S2242 = 2.0f * x_37;
    float fS1_12 = _S2242 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_28 - 0.31539157032966614f;
    float3  _S2243 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_37;
    float3  _S2244 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S2245 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S2246 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S2247 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_28 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S2248 = 1.86588168144226074f * z2_28 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S2248;
    float3  _S2249 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_37;
    float3  _S2250 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S2251 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S2252 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S2253 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_37 * fC1_12 - y_15 * fS1_12);
    float3  _S2254 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_37 * fS1_12 + y_15 * fC1_12);
    float3  _S2255 = make_float3 (pSH9_2);
    float3  _S2256 = make_float3 (0.0f);
    float3  _S2257 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2258;
    (&_S2258)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2240) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_37) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S2258)->differential_0 = _S2257;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2259;
    (&_S2259)->primal_0 = _S2256;
    (&_S2259)->differential_0 = _S2257;
    s_bwd_prop_max_0(&_S2258, &_S2259, v_rgb_2);
    float3  _S2260 = _S2254 * _S2258.differential_0;
    float3  _S2261 = (*sh_coeffs_12)[int(15)] * _S2258.differential_0;
    float3  _S2262 = _S2252 * _S2258.differential_0;
    float3  _S2263 = (*sh_coeffs_12)[int(14)] * _S2258.differential_0;
    float3  _S2264 = _S2250 * _S2258.differential_0;
    float3  _S2265 = (*sh_coeffs_12)[int(13)] * _S2258.differential_0;
    float3  _S2266 = _S2249 * _S2258.differential_0;
    float3  _S2267 = (*sh_coeffs_12)[int(12)] * _S2258.differential_0;
    float3  _S2268 = _S2251 * _S2258.differential_0;
    float3  _S2269 = (*sh_coeffs_12)[int(11)] * _S2258.differential_0;
    float3  _S2270 = _S2253 * _S2258.differential_0;
    float3  _S2271 = (*sh_coeffs_12)[int(10)] * _S2258.differential_0;
    float3  _S2272 = _S2255 * _S2258.differential_0;
    float3  _S2273 = (*sh_coeffs_12)[int(9)] * _S2258.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S2273.x + _S2273.y + _S2273.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S2261.x + _S2261.y + _S2261.z);
    float _S2274 = _S2271.x + _S2271.y + _S2271.z;
    float _S2275 = _S2263.x + _S2263.y + _S2263.z;
    float _S2276 = _S2269.x + _S2269.y + _S2269.z;
    float _S2277 = _S2265.x + _S2265.y + _S2265.z;
    float _S2278 = _S2267.x + _S2267.y + _S2267.z;
    float _S2279 = - s_diff_fC2_T_2;
    float3  _S2280 = _S2246 * _S2258.differential_0;
    float3  _S2281 = (*sh_coeffs_12)[int(8)] * _S2258.differential_0;
    float3  _S2282 = _S2244 * _S2258.differential_0;
    float3  _S2283 = (*sh_coeffs_12)[int(7)] * _S2258.differential_0;
    float3  _S2284 = _S2243 * _S2258.differential_0;
    float3  _S2285 = (*sh_coeffs_12)[int(6)] * _S2258.differential_0;
    float3  _S2286 = _S2245 * _S2258.differential_0;
    float3  _S2287 = (*sh_coeffs_12)[int(5)] * _S2258.differential_0;
    float3  _S2288 = _S2247 * _S2258.differential_0;
    float3  _S2289 = (*sh_coeffs_12)[int(4)] * _S2258.differential_0;
    float _S2290 = _S2287.x + _S2287.y + _S2287.z;
    float _S2291 = _S2283.x + _S2283.y + _S2283.z;
    float _S2292 = fTmp1B_12 * _S2274 + x_37 * s_diff_fS2_T_2 + y_15 * _S2279 + 0.54627424478530884f * (_S2289.x + _S2289.y + _S2289.z);
    float _S2293 = fTmp1B_12 * _S2275 + y_15 * s_diff_fS2_T_2 + x_37 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S2281.x + _S2281.y + _S2281.z);
    float _S2294 = y_15 * - _S2293;
    float _S2295 = x_37 * _S2293;
    float _S2296 = z_12 * (1.86588168144226074f * (z_12 * _S2278) + -2.28522896766662598f * (y_15 * _S2276 + x_37 * _S2277) + 0.94617468118667603f * (_S2285.x + _S2285.y + _S2285.z));
    float3  _S2297 = make_float3 (0.48860251903533936f) * _S2258.differential_0;
    float3  _S2298 = - _S2297;
    float3  _S2299 = _S2237 * _S2298;
    float3  _S2300 = (*sh_coeffs_12)[int(3)] * _S2298;
    float3  _S2301 = _S2239 * _S2297;
    float3  _S2302 = (*sh_coeffs_12)[int(2)] * _S2297;
    float3  _S2303 = _S2241 * _S2297;
    float3  _S2304 = (*sh_coeffs_12)[int(1)] * _S2297;
    float _S2305 = (_S2248 * _S2278 + 1.44530570507049561f * (fS1_12 * _S2274 + fC1_12 * _S2275) + -1.09254848957061768f * (y_15 * _S2290 + x_37 * _S2291) + _S2296 + _S2296 + _S2302.x + _S2302.y + _S2302.z) / _S2238;
    float _S2306 = _S2236 * _S2305;
    float _S2307 = (fTmp0C_12 * _S2276 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S2279 + fTmp0B_12 * _S2290 + _S2242 * _S2292 + _S2294 + _S2294 + - (_S2304.x + _S2304.y + _S2304.z)) / _S2238;
    float _S2308 = _S2236 * _S2307;
    float _S2309 = (fTmp0C_12 * _S2277 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S2291 + 2.0f * (y_15 * _S2292) + _S2295 + _S2295 + _S2300.x + _S2300.y + _S2300.z) / _S2238;
    float _S2310 = _S2236 * _S2309;
    float _S2311 = _S2234 * - _S2305 + _S2233 * - _S2307 + _S2232 * - _S2309;
    DiffPair_float_0 _S2312;
    (&_S2312)->primal_0 = _S2235;
    (&_S2312)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2312, _S2311);
    float _S2313 = _S2234 * _S2312.differential_0;
    float _S2314 = _S2233 * _S2312.differential_0;
    float _S2315 = _S2232 * _S2312.differential_0;
    float3  _S2316 = make_float3 (0.282094806432724f) * _S2258.differential_0;
    float3  _S2317 = make_float3 (_S2310 + _S2315 + _S2315, _S2308 + _S2314 + _S2314, _S2306 + _S2313 + _S2313);
    float3  _S2318 = - - _S2317;
    Matrix<float, 3, 3>  _S2319 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2320;
    (&_S2320)->primal_0 = _S2203;
    (&_S2320)->differential_0 = _S2319;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2321;
    (&_S2321)->primal_0 = t_15;
    (&_S2321)->differential_0 = _S2257;
    s_bwd_prop_mul_1(&_S2320, &_S2321, _S2318);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2322 = _S2320;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2323 = _S2321;
    float2  _S2324 = make_float2 (0.0f);
    float2  _S2325 = _S2324;
    *&((&_S2325)->y) = v_conic_2.z;
    float2  _S2326 = _S2324;
    *&((&_S2326)->y) = v_conic_2.y;
    *&((&_S2326)->x) = v_conic_2.x;
    float _S2327 = 0.5f * v_depth_2;
    DiffPair_float_0 _S2328;
    (&_S2328)->primal_0 = _S2230;
    (&_S2328)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2328, _S2327);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2329;
    (&_S2329)->primal_0 = mean_c_12;
    (&_S2329)->differential_0 = _S2257;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2330;
    (&_S2330)->primal_0 = mean_c_12;
    (&_S2330)->differential_0 = _S2257;
    s_bwd_prop_dot_0(&_S2329, &_S2330, _S2328.differential_0);
    DiffPair_float_0 _S2331;
    (&_S2331)->primal_0 = _S2229;
    (&_S2331)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2331, 0.0f);
    DiffPair_float_0 _S2332;
    (&_S2332)->primal_0 = _S2228;
    (&_S2332)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2332, 0.0f);
    DiffPair_float_0 _S2333;
    (&_S2333)->primal_0 = 3.32999992370605469f;
    (&_S2333)->differential_0 = 0.0f;
    DiffPair_float_0 _S2334;
    (&_S2334)->primal_0 = _S2227;
    (&_S2334)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2333, &_S2334, 0.0f);
    DiffPair_float_0 _S2335;
    (&_S2335)->primal_0 = _S2226;
    (&_S2335)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2335, _S2334.differential_0);
    float _S2336 = 2.0f * _S2335.differential_0;
    DiffPair_float_0 _S2337;
    (&_S2337)->primal_0 = _S2225;
    (&_S2337)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2337, _S2336);
    float _S2338 = v_opacity_2 + 254.9999847412109375f * _S2337.differential_0;
    Matrix<float, 2, 2>  _S2339 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2340 = _S2339;
    _S2340[int(1)] = _S2325;
    _S2340[int(0)] = _S2326;
    Matrix<float, 2, 2>  _S2341 = _S2340;
    FixedArray<float3 , 16>  _S2342;
    _S2342[int(0)] = _S2257;
    _S2342[int(1)] = _S2257;
    _S2342[int(2)] = _S2257;
    _S2342[int(3)] = _S2257;
    _S2342[int(4)] = _S2257;
    _S2342[int(5)] = _S2257;
    _S2342[int(6)] = _S2257;
    _S2342[int(7)] = _S2257;
    _S2342[int(8)] = _S2257;
    _S2342[int(9)] = _S2257;
    _S2342[int(10)] = _S2257;
    _S2342[int(11)] = _S2257;
    _S2342[int(12)] = _S2257;
    _S2342[int(13)] = _S2257;
    _S2342[int(14)] = _S2257;
    _S2342[int(15)] = _S2257;
    _S2342[int(7)] = _S2282;
    _S2342[int(0)] = _S2316;
    _S2342[int(1)] = _S2303;
    _S2342[int(2)] = _S2301;
    _S2342[int(3)] = _S2299;
    _S2342[int(4)] = _S2288;
    _S2342[int(5)] = _S2286;
    _S2342[int(6)] = _S2284;
    _S2342[int(15)] = _S2260;
    _S2342[int(8)] = _S2280;
    _S2342[int(9)] = _S2272;
    _S2342[int(10)] = _S2270;
    _S2342[int(11)] = _S2268;
    _S2342[int(12)] = _S2266;
    _S2342[int(13)] = _S2264;
    _S2342[int(14)] = _S2262;
    float3  _S2343 = _S2342[int(0)];
    float3  _S2344 = _S2342[int(1)];
    float3  _S2345 = _S2342[int(2)];
    float3  _S2346 = _S2342[int(3)];
    float3  _S2347 = _S2342[int(4)];
    float3  _S2348 = _S2342[int(5)];
    float3  _S2349 = _S2342[int(6)];
    float3  _S2350 = _S2342[int(7)];
    float3  _S2351 = _S2342[int(8)];
    float3  _S2352 = _S2342[int(9)];
    float3  _S2353 = _S2342[int(10)];
    float3  _S2354 = _S2342[int(11)];
    float3  _S2355 = _S2342[int(12)];
    float3  _S2356 = _S2342[int(13)];
    float3  _S2357 = _S2342[int(14)];
    float3  _S2358 = _S2342[int(15)];
    float3  _S2359 = _S2330.differential_0 + _S2329.differential_0;
    float2  _S2360 = make_float2 (0.0f, _S2331.differential_0);
    float2  _S2361 = make_float2 (_S2332.differential_0, 0.0f);
    float _S2362;
    if(antialiased_12)
    {
        float _S2363 = _S2222 * _S2338;
        _S2224 = _S2217 * _S2338;
        _S2362 = _S2363;
    }
    else
    {
        _S2224 = _S2338;
        _S2362 = 0.0f;
    }
    float _S2364 = - (_S2224 / _S2223);
    DiffPair_float_0 _S2365;
    (&_S2365)->primal_0 = _S2220;
    (&_S2365)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2365, _S2364);
    float _S2366 = - _S2365.differential_0;
    float _S2367 = invdet_9 * _S2341.rows[int(1)].y;
    float _S2368 = - (invdet_9 * _S2341.rows[int(1)].x);
    float _S2369 = - (invdet_9 * _S2341.rows[int(0)].y);
    float _S2370 = invdet_9 * _S2341.rows[int(0)].x;
    float _S2371 = - ((_S2209 * _S2341.rows[int(1)].y + _S2219 * _S2341.rows[int(1)].x + _S2218 * _S2341.rows[int(0)].y + _S2211 * _S2341.rows[int(0)].x) / _S2215);
    DiffPair_float_0 _S2372;
    (&_S2372)->primal_0 = _S2216;
    (&_S2372)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2372, _S2362);
    DiffPair_float_0 _S2373;
    (&_S2373)->primal_0 = 0.0f;
    (&_S2373)->differential_0 = 0.0f;
    DiffPair_float_0 _S2374;
    (&_S2374)->primal_0 = _S2214;
    (&_S2374)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2373, &_S2374, _S2372.differential_0);
    float _S2375 = _S2374.differential_0 / _S2215;
    float s_diff_det_orig_T_2 = det_blur_8 * _S2375;
    float _S2376 = _S2371 + det_orig_13 * - _S2375;
    float _S2377 = - _S2376;
    float _S2378 = _S2209 * _S2376;
    float _S2379 = _S2211 * _S2376;
    Matrix<float, 2, 2>  _S2380 = _S2339;
    _S2380[int(1)] = _S2360;
    _S2380[int(0)] = _S2361;
    _S2210 = _S2380;
    *&(((&_S2210)->rows + (int(1)))->y) = 0.0f;
    float _S2381 = _S2370 + _S2378 + _S2380.rows[int(1)].y;
    *&(((&_S2210)->rows + (int(0)))->x) = 0.0f;
    float _S2382 = _S2367 + _S2379 + _S2380.rows[int(0)].x;
    float _S2383 = _S2377 + - s_diff_det_orig_T_2;
    float _S2384 = _S2368 + _S2207.rows[int(0)].y * _S2383;
    float _S2385 = _S2369 + _S2207.rows[int(1)].x * _S2383;
    float _S2386 = _S2207.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S2387 = _S2381 + _S2207.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S2388 = _S2324;
    *&((&_S2388)->x) = _S2384;
    *&((&_S2388)->y) = _S2387;
    float _S2389 = _S2382 + _S2386;
    float2  _S2390 = _S2324;
    *&((&_S2390)->y) = _S2385;
    *&((&_S2390)->x) = _S2389;
    float _S2391 = fy_21 * v_mean2d_2.y;
    float _S2392 = fx_21 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S2393 = _S2339;
    _S2393[int(1)] = _S2388;
    _S2393[int(0)] = _S2390;
    Matrix<float, 2, 2>  _S2394 = _S2210 + _S2393;
    Matrix<float, 2, 3>  _S2395 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2396;
    (&_S2396)->primal_0 = _S2205;
    (&_S2396)->differential_0 = _S2395;
    Matrix<float, 3, 2>  _S2397 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2398;
    (&_S2398)->primal_0 = _S2206;
    (&_S2398)->differential_0 = _S2397;
    s_bwd_prop_mul_2(&_S2396, &_S2398, _S2394);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2399;
    (&_S2399)->primal_0 = J_13;
    (&_S2399)->differential_0 = _S2395;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2400;
    (&_S2400)->primal_0 = _S2204;
    (&_S2400)->differential_0 = _S2319;
    s_bwd_prop_mul_3(&_S2399, &_S2400, _S2396.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2401;
    (&_S2401)->primal_0 = _S2202;
    (&_S2401)->differential_0 = _S2319;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2402;
    (&_S2402)->primal_0 = _S2203;
    (&_S2402)->differential_0 = _S2319;
    s_bwd_prop_mul_4(&_S2401, &_S2402, _S2400.differential_0);
    Matrix<float, 3, 3>  _S2403 = transpose_0(_S2402.differential_0 + _S2322.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2404;
    (&_S2404)->primal_0 = R_16;
    (&_S2404)->differential_0 = _S2319;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2405;
    (&_S2405)->primal_0 = _S2201;
    (&_S2405)->differential_0 = _S2319;
    s_bwd_prop_mul_4(&_S2404, &_S2405, _S2401.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2406;
    (&_S2406)->primal_0 = _S2199;
    (&_S2406)->differential_0 = _S2319;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2407;
    (&_S2407)->primal_0 = _S2200;
    (&_S2407)->differential_0 = _S2319;
    s_bwd_prop_mul_4(&_S2406, &_S2407, _S2405.differential_0);
    Matrix<float, 3, 3>  _S2408 = _S2406.differential_0 + transpose_0(_S2407.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2409;
    (&_S2409)->primal_0 = _S2198;
    (&_S2409)->differential_0 = _S2319;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2410;
    (&_S2410)->primal_0 = S_2;
    (&_S2410)->differential_0 = _S2319;
    s_bwd_prop_mul_4(&_S2409, &_S2410, _S2408);
    Matrix<float, 3, 3>  _S2411 = transpose_0(_S2409.differential_0);
    float _S2412 = 2.0f * - _S2411.rows[int(2)].z;
    float _S2413 = 2.0f * _S2411.rows[int(2)].y;
    float _S2414 = 2.0f * _S2411.rows[int(2)].x;
    float _S2415 = 2.0f * _S2411.rows[int(1)].z;
    float _S2416 = 2.0f * - _S2411.rows[int(1)].y;
    float _S2417 = 2.0f * _S2411.rows[int(1)].x;
    float _S2418 = 2.0f * _S2411.rows[int(0)].z;
    float _S2419 = 2.0f * _S2411.rows[int(0)].y;
    float _S2420 = 2.0f * - _S2411.rows[int(0)].x;
    float _S2421 = - _S2417 + _S2419;
    float _S2422 = _S2414 + - _S2418;
    float _S2423 = - _S2413 + _S2415;
    float _S2424 = _S2413 + _S2415;
    float _S2425 = _S2414 + _S2418;
    float _S2426 = _S2417 + _S2419;
    float _S2427 = quat_15.w * (_S2416 + _S2420);
    float _S2428 = quat_15.z * (_S2412 + _S2420);
    float _S2429 = quat_15.y * (_S2412 + _S2416);
    float _S2430 = quat_15.x * _S2421 + quat_15.z * _S2424 + quat_15.y * _S2425 + _S2427 + _S2427;
    float _S2431 = quat_15.x * _S2422 + quat_15.w * _S2424 + quat_15.y * _S2426 + _S2428 + _S2428;
    float _S2432 = quat_15.x * _S2423 + quat_15.w * _S2425 + quat_15.z * _S2426 + _S2429 + _S2429;
    float _S2433 = quat_15.w * _S2421 + quat_15.z * _S2422 + quat_15.y * _S2423;
    float3  _S2434 = _S2257;
    *&((&_S2434)->z) = _S2410.differential_0.rows[int(2)].z;
    *&((&_S2434)->y) = _S2410.differential_0.rows[int(1)].y;
    *&((&_S2434)->x) = _S2410.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2435;
    (&_S2435)->primal_0 = scale_14;
    (&_S2435)->differential_0 = _S2257;
    s_bwd_prop_exp_1(&_S2435, _S2434);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2436 = _S2435;
    float3  _S2437 = _S2257;
    *&((&_S2437)->y) = _S2391;
    *&((&_S2437)->x) = _S2392;
    float3  _S2438 = _S2359 + _S2437;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2439;
    (&_S2439)->primal_0 = R_16;
    (&_S2439)->differential_0 = _S2319;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2440;
    (&_S2440)->primal_0 = mean_12;
    (&_S2440)->differential_0 = _S2257;
    s_bwd_prop_mul_1(&_S2439, &_S2440, _S2438);
    float3  _S2441 = _S2438 + _S2323.differential_0;
    Matrix<float, 3, 3>  _S2442 = _S2403 + _S2404.differential_0 + _S2439.differential_0;
    float4  _S2443 = make_float4 (0.0f);
    *&((&_S2443)->w) = _S2430;
    *&((&_S2443)->z) = _S2431;
    *&((&_S2443)->y) = _S2432;
    *&((&_S2443)->x) = _S2433;
    float4  _S2444 = _S2443;
    float3  _S2445 = _S2440.differential_0 + _S2317;
    *v_mean_2 = _S2445;
    *v_quat_2 = _S2444;
    *v_scale_2 = _S2436.differential_0;
    *v_in_opacity_2 = _S2366;
    (*v_sh_coeffs_2)[int(0)] = _S2343;
    (*v_sh_coeffs_2)[int(1)] = _S2344;
    (*v_sh_coeffs_2)[int(2)] = _S2345;
    (*v_sh_coeffs_2)[int(3)] = _S2346;
    (*v_sh_coeffs_2)[int(4)] = _S2347;
    (*v_sh_coeffs_2)[int(5)] = _S2348;
    (*v_sh_coeffs_2)[int(6)] = _S2349;
    (*v_sh_coeffs_2)[int(7)] = _S2350;
    (*v_sh_coeffs_2)[int(8)] = _S2351;
    (*v_sh_coeffs_2)[int(9)] = _S2352;
    (*v_sh_coeffs_2)[int(10)] = _S2353;
    (*v_sh_coeffs_2)[int(11)] = _S2354;
    (*v_sh_coeffs_2)[int(12)] = _S2355;
    (*v_sh_coeffs_2)[int(13)] = _S2356;
    (*v_sh_coeffs_2)[int(14)] = _S2357;
    (*v_sh_coeffs_2)[int(15)] = _S2358;
    *v_R_3 = _S2442;
    *v_t_3 = _S2441;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2446;
};

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_22, float fy_22, float cx_17, float cy_17, FixedArray<float, 10>  * dist_coeffs_28, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    Matrix<float, 2, 2>  _S2447 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2448;
    (&_S2448)->_S2446 = _S2447;
    float2  _S2449 = make_float2 (0.0f);
    float3  _S2450 = make_float3 (0.0f);
    float4  intrins_12 = make_float4 (fx_22, fy_22, cx_17, cy_17);
    float3  _S2451 = s_primal_ctx_exp_0(scale_15);
    float _S2452 = quat_16.y;
    float x2_16 = _S2452 * _S2452;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_29 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S2453 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_29), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_29), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16))));
    FixedArray<float3 , 7>  _S2454 = {
        _S2450, _S2450, _S2450, _S2450, _S2450, _S2450, _S2450
    };
    FixedArray<float, 7>  _S2455 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2456;
    (&_S2456)->p_0 = _S2454;
    (&_S2456)->w_mean_0 = _S2455;
    (&_S2456)->w_cov_0 = _S2455;
    (&_S2456)->p_0[int(0)] = mean_13;
    SigmaPoints_0 _S2457 = _S2456;
    (&_S2457)->w_mean_0[int(0)] = 0.0f;
    (&_S2457)->w_cov_0[int(0)] = 2.0f;
    float _S2458 = s_primal_ctx_sqrt_0(3.0f);
    float _S2459 = _S2458 * _S2451.x;
    float3  delta_12 = make_float3 (_S2459) * _S2453.rows[0U];
    float3  _S2460 = mean_13 + delta_12;
    float3  _S2461 = mean_13 - delta_12;
    float _S2462 = _S2458 * _S2451.y;
    float3  delta_13 = make_float3 (_S2462) * _S2453.rows[1U];
    float3  _S2463 = mean_13 + delta_13;
    float3  _S2464 = mean_13 - delta_13;
    float _S2465 = _S2458 * _S2451.z;
    float3  delta_14 = make_float3 (_S2465) * _S2453.rows[2U];
    float3  _S2466 = mean_13 + delta_14;
    float3  _S2467 = mean_13 - delta_14;
    (&_S2457)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2457)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2457)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2457)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2457)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2457)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2457)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2457)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2457)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2457)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2457)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2457)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2468 = _S2456;
    (&_S2457)->p_0[0U] = s_primal_ctx_mul_1(R_17, _S2456.p_0[0U]) + t_16;
    (&_S2457)->p_0[1U] = s_primal_ctx_mul_1(R_17, _S2460) + t_16;
    (&_S2457)->p_0[2U] = s_primal_ctx_mul_1(R_17, _S2463) + t_16;
    (&_S2457)->p_0[3U] = s_primal_ctx_mul_1(R_17, _S2466) + t_16;
    (&_S2457)->p_0[4U] = s_primal_ctx_mul_1(R_17, _S2461) + t_16;
    (&_S2457)->p_0[5U] = s_primal_ctx_mul_1(R_17, _S2464) + t_16;
    (&_S2457)->p_0[6U] = s_primal_ctx_mul_1(R_17, _S2467) + t_16;
    float2  _S2469 = _S2449;
    Matrix<float, 2, 2>  _S2470 = _S2447;
    SigmaPoints_0 _S2471 = _S2457;
    bool _S2472 = persp_proj_3dgs_ut_1(&_S2471, intrins_12, dist_coeffs_28, image_width_13, image_height_13, &_S2470, &_S2469);
    (&_S2448)->_S2446 = _S2470;
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2473 = _S2448;
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_13) + t_16;
    float3  _S2474 = make_float3 (_S2459);
    float3  _S2475 = make_float3 (_S2462);
    float3  _S2476 = make_float3 (_S2465);
    float _S2477 = _S2448._S2446.rows[int(0)].y * _S2448._S2446.rows[int(1)].x;
    float det_orig_14 = _S2448._S2446.rows[int(0)].x * _S2448._S2446.rows[int(1)].y - _S2477;
    float _S2478 = _S2448._S2446.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2479 = _S2448._S2446;
    *&(((&_S2479)->rows + (int(0)))->x) = _S2478;
    float _S2480 = _S2448._S2446.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2479)->rows + (int(1)))->y) = _S2480;
    Matrix<float, 2, 2>  _S2481 = _S2479;
    Matrix<float, 2, 2>  _S2482 = _S2479;
    float det_blur_9 = _S2478 * _S2480 - _S2477;
    float _S2483 = det_orig_14 / det_blur_9;
    float _S2484 = det_blur_9 * det_blur_9;
    float _S2485 = s_primal_ctx_max_0(0.0f, _S2483);
    float _S2486 = s_primal_ctx_sqrt_0(_S2485);
    float _S2487 = - in_opacity_13;
    float _S2488 = 1.0f + s_primal_ctx_exp_1(_S2487);
    float _S2489 = 1.0f / _S2488;
    float _S2490 = _S2488 * _S2488;
    float _S2491;
    if(antialiased_13)
    {
        _S2491 = _S2489 * _S2486;
    }
    else
    {
        _S2491 = _S2489;
    }
    float _S2492 = _S2491 / 0.00392156885936856f;
    float _S2493 = 2.0f * s_primal_ctx_log_0(_S2492);
    float _S2494 = s_primal_ctx_sqrt_0(_S2493);
    float _S2495 = _S2481.rows[int(0)].x;
    float _S2496 = _S2482.rows[int(1)].y;
    float _S2497 = s_primal_ctx_dot_0(mean_c_13, mean_c_13) + 9.99999997475242708e-07f;
    float3  _S2498 = - scale_15;
    Matrix<float, 3, 3>  _S2499 = transpose_0(R_17);
    float3  _S2500 = mean_13 - - s_primal_ctx_mul_1(_S2499, t_16);
    float _S2501 = _S2500.x;
    float _S2502 = _S2500.y;
    float _S2503 = _S2500.z;
    float _S2504 = _S2501 * _S2501 + _S2502 * _S2502 + _S2503 * _S2503;
    float _S2505 = s_primal_ctx_sqrt_0(_S2504);
    float x_38 = _S2501 / _S2505;
    float3  _S2506 = make_float3 (x_38);
    float _S2507 = _S2505 * _S2505;
    float y_16 = _S2502 / _S2505;
    float z_13 = _S2503 / _S2505;
    float3  _S2508 = make_float3 (z_13);
    float _S2509 = - y_16;
    float3  _S2510 = make_float3 (_S2509);
    float z2_30 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_38 * x_38 - y_16 * y_16;
    float _S2511 = 2.0f * x_38;
    float fS1_13 = _S2511 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_30 - 0.31539157032966614f;
    float3  _S2512 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_38;
    float3  _S2513 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S2514 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S2515 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S2516 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_30 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S2517 = 1.86588168144226074f * z2_30 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S2517;
    float3  _S2518 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_38;
    float3  _S2519 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S2520 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S2521 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S2522 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_38 * fC1_13 - y_16 * fS1_13);
    float3  _S2523 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_38 * fS1_13 + y_16 * fC1_13);
    float3  _S2524 = make_float3 (pSH9_3);
    float3  _S2525 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2526;
    (&_S2526)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2509) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_38) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S2526)->differential_0 = _S2450;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2527;
    (&_S2527)->primal_0 = _S2525;
    (&_S2527)->differential_0 = _S2450;
    s_bwd_prop_max_0(&_S2526, &_S2527, v_rgb_3);
    float3  _S2528 = _S2523 * _S2526.differential_0;
    float3  _S2529 = (*sh_coeffs_13)[int(15)] * _S2526.differential_0;
    float3  _S2530 = _S2521 * _S2526.differential_0;
    float3  _S2531 = (*sh_coeffs_13)[int(14)] * _S2526.differential_0;
    float3  _S2532 = _S2519 * _S2526.differential_0;
    float3  _S2533 = (*sh_coeffs_13)[int(13)] * _S2526.differential_0;
    float3  _S2534 = _S2518 * _S2526.differential_0;
    float3  _S2535 = (*sh_coeffs_13)[int(12)] * _S2526.differential_0;
    float3  _S2536 = _S2520 * _S2526.differential_0;
    float3  _S2537 = (*sh_coeffs_13)[int(11)] * _S2526.differential_0;
    float3  _S2538 = _S2522 * _S2526.differential_0;
    float3  _S2539 = (*sh_coeffs_13)[int(10)] * _S2526.differential_0;
    float3  _S2540 = _S2524 * _S2526.differential_0;
    float3  _S2541 = (*sh_coeffs_13)[int(9)] * _S2526.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2541.x + _S2541.y + _S2541.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S2529.x + _S2529.y + _S2529.z);
    float _S2542 = _S2539.x + _S2539.y + _S2539.z;
    float _S2543 = _S2531.x + _S2531.y + _S2531.z;
    float _S2544 = _S2537.x + _S2537.y + _S2537.z;
    float _S2545 = _S2533.x + _S2533.y + _S2533.z;
    float _S2546 = _S2535.x + _S2535.y + _S2535.z;
    float _S2547 = - s_diff_fC2_T_3;
    float3  _S2548 = _S2515 * _S2526.differential_0;
    float3  _S2549 = (*sh_coeffs_13)[int(8)] * _S2526.differential_0;
    float3  _S2550 = _S2513 * _S2526.differential_0;
    float3  _S2551 = (*sh_coeffs_13)[int(7)] * _S2526.differential_0;
    float3  _S2552 = _S2512 * _S2526.differential_0;
    float3  _S2553 = (*sh_coeffs_13)[int(6)] * _S2526.differential_0;
    float3  _S2554 = _S2514 * _S2526.differential_0;
    float3  _S2555 = (*sh_coeffs_13)[int(5)] * _S2526.differential_0;
    float3  _S2556 = _S2516 * _S2526.differential_0;
    float3  _S2557 = (*sh_coeffs_13)[int(4)] * _S2526.differential_0;
    float _S2558 = _S2555.x + _S2555.y + _S2555.z;
    float _S2559 = _S2551.x + _S2551.y + _S2551.z;
    float _S2560 = fTmp1B_13 * _S2542 + x_38 * s_diff_fS2_T_3 + y_16 * _S2547 + 0.54627424478530884f * (_S2557.x + _S2557.y + _S2557.z);
    float _S2561 = fTmp1B_13 * _S2543 + y_16 * s_diff_fS2_T_3 + x_38 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2549.x + _S2549.y + _S2549.z);
    float _S2562 = y_16 * - _S2561;
    float _S2563 = x_38 * _S2561;
    float _S2564 = z_13 * (1.86588168144226074f * (z_13 * _S2546) + -2.28522896766662598f * (y_16 * _S2544 + x_38 * _S2545) + 0.94617468118667603f * (_S2553.x + _S2553.y + _S2553.z));
    float3  _S2565 = make_float3 (0.48860251903533936f) * _S2526.differential_0;
    float3  _S2566 = - _S2565;
    float3  _S2567 = _S2506 * _S2566;
    float3  _S2568 = (*sh_coeffs_13)[int(3)] * _S2566;
    float3  _S2569 = _S2508 * _S2565;
    float3  _S2570 = (*sh_coeffs_13)[int(2)] * _S2565;
    float3  _S2571 = _S2510 * _S2565;
    float3  _S2572 = (*sh_coeffs_13)[int(1)] * _S2565;
    float _S2573 = (_S2517 * _S2546 + 1.44530570507049561f * (fS1_13 * _S2542 + fC1_13 * _S2543) + -1.09254848957061768f * (y_16 * _S2558 + x_38 * _S2559) + _S2564 + _S2564 + _S2570.x + _S2570.y + _S2570.z) / _S2507;
    float _S2574 = _S2505 * _S2573;
    float _S2575 = (fTmp0C_13 * _S2544 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2547 + fTmp0B_13 * _S2558 + _S2511 * _S2560 + _S2562 + _S2562 + - (_S2572.x + _S2572.y + _S2572.z)) / _S2507;
    float _S2576 = _S2505 * _S2575;
    float _S2577 = (fTmp0C_13 * _S2545 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2559 + 2.0f * (y_16 * _S2560) + _S2563 + _S2563 + _S2568.x + _S2568.y + _S2568.z) / _S2507;
    float _S2578 = _S2505 * _S2577;
    float _S2579 = _S2503 * - _S2573 + _S2502 * - _S2575 + _S2501 * - _S2577;
    DiffPair_float_0 _S2580;
    (&_S2580)->primal_0 = _S2504;
    (&_S2580)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2580, _S2579);
    float _S2581 = _S2503 * _S2580.differential_0;
    float _S2582 = _S2502 * _S2580.differential_0;
    float _S2583 = _S2501 * _S2580.differential_0;
    float3  _S2584 = make_float3 (0.282094806432724f) * _S2526.differential_0;
    float3  _S2585 = make_float3 (_S2578 + _S2583 + _S2583, _S2576 + _S2582 + _S2582, _S2574 + _S2581 + _S2581);
    float3  _S2586 = - - _S2585;
    Matrix<float, 3, 3>  _S2587 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2588;
    (&_S2588)->primal_0 = _S2499;
    (&_S2588)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2589;
    (&_S2589)->primal_0 = t_16;
    (&_S2589)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2588, &_S2589, _S2586);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2590 = _S2589;
    Matrix<float, 3, 3>  _S2591 = transpose_0(_S2588.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2592;
    (&_S2592)->primal_0 = _S2498;
    (&_S2592)->differential_0 = _S2450;
    s_bwd_prop_exp_1(&_S2592, v_conic_3);
    float3  _S2593 = - _S2592.differential_0;
    float _S2594 = 0.5f * v_depth_3;
    DiffPair_float_0 _S2595;
    (&_S2595)->primal_0 = _S2497;
    (&_S2595)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2595, _S2594);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2596;
    (&_S2596)->primal_0 = mean_c_13;
    (&_S2596)->differential_0 = _S2450;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2597;
    (&_S2597)->primal_0 = mean_c_13;
    (&_S2597)->differential_0 = _S2450;
    s_bwd_prop_dot_0(&_S2596, &_S2597, _S2595.differential_0);
    DiffPair_float_0 _S2598;
    (&_S2598)->primal_0 = _S2496;
    (&_S2598)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2598, 0.0f);
    DiffPair_float_0 _S2599;
    (&_S2599)->primal_0 = _S2495;
    (&_S2599)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2599, 0.0f);
    DiffPair_float_0 _S2600;
    (&_S2600)->primal_0 = 3.32999992370605469f;
    (&_S2600)->differential_0 = 0.0f;
    DiffPair_float_0 _S2601;
    (&_S2601)->primal_0 = _S2494;
    (&_S2601)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2600, &_S2601, 0.0f);
    DiffPair_float_0 _S2602;
    (&_S2602)->primal_0 = _S2493;
    (&_S2602)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2602, _S2601.differential_0);
    float _S2603 = 2.0f * _S2602.differential_0;
    DiffPair_float_0 _S2604;
    (&_S2604)->primal_0 = _S2492;
    (&_S2604)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2604, _S2603);
    float _S2605 = v_opacity_3 + 254.9999847412109375f * _S2604.differential_0;
    FixedArray<float3 , 16>  _S2606;
    _S2606[int(0)] = _S2450;
    _S2606[int(1)] = _S2450;
    _S2606[int(2)] = _S2450;
    _S2606[int(3)] = _S2450;
    _S2606[int(4)] = _S2450;
    _S2606[int(5)] = _S2450;
    _S2606[int(6)] = _S2450;
    _S2606[int(7)] = _S2450;
    _S2606[int(8)] = _S2450;
    _S2606[int(9)] = _S2450;
    _S2606[int(10)] = _S2450;
    _S2606[int(11)] = _S2450;
    _S2606[int(12)] = _S2450;
    _S2606[int(13)] = _S2450;
    _S2606[int(14)] = _S2450;
    _S2606[int(15)] = _S2450;
    _S2606[int(7)] = _S2550;
    _S2606[int(0)] = _S2584;
    _S2606[int(1)] = _S2571;
    _S2606[int(2)] = _S2569;
    _S2606[int(3)] = _S2567;
    _S2606[int(4)] = _S2556;
    _S2606[int(5)] = _S2554;
    _S2606[int(6)] = _S2552;
    _S2606[int(15)] = _S2528;
    _S2606[int(8)] = _S2548;
    _S2606[int(9)] = _S2540;
    _S2606[int(10)] = _S2538;
    _S2606[int(11)] = _S2536;
    _S2606[int(12)] = _S2534;
    _S2606[int(13)] = _S2532;
    _S2606[int(14)] = _S2530;
    float3  _S2607 = _S2606[int(0)];
    float3  _S2608 = _S2606[int(1)];
    float3  _S2609 = _S2606[int(2)];
    float3  _S2610 = _S2606[int(3)];
    float3  _S2611 = _S2606[int(4)];
    float3  _S2612 = _S2606[int(5)];
    float3  _S2613 = _S2606[int(6)];
    float3  _S2614 = _S2606[int(7)];
    float3  _S2615 = _S2606[int(8)];
    float3  _S2616 = _S2606[int(9)];
    float3  _S2617 = _S2606[int(10)];
    float3  _S2618 = _S2606[int(11)];
    float3  _S2619 = _S2606[int(12)];
    float3  _S2620 = _S2606[int(13)];
    float3  _S2621 = _S2606[int(14)];
    float3  _S2622 = _S2606[int(15)];
    float3  _S2623 = _S2597.differential_0 + _S2596.differential_0;
    float2  _S2624 = make_float2 (0.0f, _S2598.differential_0);
    float2  _S2625 = make_float2 (_S2599.differential_0, 0.0f);
    float _S2626;
    if(antialiased_13)
    {
        float _S2627 = _S2489 * _S2605;
        _S2491 = _S2486 * _S2605;
        _S2626 = _S2627;
    }
    else
    {
        _S2491 = _S2605;
        _S2626 = 0.0f;
    }
    float _S2628 = - (_S2491 / _S2490);
    DiffPair_float_0 _S2629;
    (&_S2629)->primal_0 = _S2487;
    (&_S2629)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2629, _S2628);
    float _S2630 = - _S2629.differential_0;
    DiffPair_float_0 _S2631;
    (&_S2631)->primal_0 = _S2485;
    (&_S2631)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2631, _S2626);
    DiffPair_float_0 _S2632;
    (&_S2632)->primal_0 = 0.0f;
    (&_S2632)->differential_0 = 0.0f;
    DiffPair_float_0 _S2633;
    (&_S2633)->primal_0 = _S2483;
    (&_S2633)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2632, &_S2633, _S2631.differential_0);
    float _S2634 = _S2633.differential_0 / _S2484;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2634;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2634;
    float _S2635 = - s_diff_det_blur_T_0;
    float _S2636 = _S2478 * s_diff_det_blur_T_0;
    float _S2637 = _S2480 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2638 = _S2447;
    _S2638[int(1)] = _S2624;
    _S2638[int(0)] = _S2625;
    float _S2639 = _S2637 + _S2638.rows[int(0)].x;
    float _S2640 = _S2635 + - s_diff_det_orig_T_3;
    float _S2641 = _S2473._S2446.rows[int(0)].y * _S2640;
    float _S2642 = _S2473._S2446.rows[int(1)].x * _S2640;
    float _S2643 = _S2473._S2446.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2644 = _S2636 + _S2638.rows[int(1)].y + _S2473._S2446.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2645 = _S2449;
    *&((&_S2645)->x) = _S2641;
    *&((&_S2645)->y) = _S2644;
    float _S2646 = _S2639 + _S2643;
    float2  _S2647 = _S2449;
    *&((&_S2647)->y) = _S2642;
    *&((&_S2647)->x) = _S2646;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2648;
    (&_S2648)->primal_0 = R_17;
    (&_S2648)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2649;
    (&_S2649)->primal_0 = _S2467;
    (&_S2649)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2648, &_S2649, _S2450);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2650;
    (&_S2650)->primal_0 = R_17;
    (&_S2650)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2651;
    (&_S2651)->primal_0 = _S2464;
    (&_S2651)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2650, &_S2651, _S2450);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2652;
    (&_S2652)->primal_0 = R_17;
    (&_S2652)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2653;
    (&_S2653)->primal_0 = _S2461;
    (&_S2653)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2652, &_S2653, _S2450);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2654;
    (&_S2654)->primal_0 = R_17;
    (&_S2654)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2655;
    (&_S2655)->primal_0 = _S2466;
    (&_S2655)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2654, &_S2655, _S2450);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2656;
    (&_S2656)->primal_0 = R_17;
    (&_S2656)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2657;
    (&_S2657)->primal_0 = _S2463;
    (&_S2657)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2656, &_S2657, _S2450);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2658;
    (&_S2658)->primal_0 = R_17;
    (&_S2658)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2659;
    (&_S2659)->primal_0 = _S2460;
    (&_S2659)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2658, &_S2659, _S2450);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2660;
    (&_S2660)->primal_0 = R_17;
    (&_S2660)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2661;
    (&_S2661)->primal_0 = _S2468.p_0[0U];
    (&_S2661)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2660, &_S2661, _S2450);
    float3  _S2662 = - _S2649.differential_0 + _S2655.differential_0;
    float3  _S2663 = _S2476 * _S2662;
    float3  _S2664 = _S2453.rows[2U] * _S2662;
    float _S2665 = _S2458 * (_S2664.x + _S2664.y + _S2664.z);
    float3  _S2666 = - _S2651.differential_0 + _S2657.differential_0;
    float3  _S2667 = _S2475 * _S2666;
    float3  _S2668 = _S2453.rows[1U] * _S2666;
    float _S2669 = _S2458 * (_S2668.x + _S2668.y + _S2668.z);
    float3  _S2670 = - _S2653.differential_0 + _S2659.differential_0;
    float3  _S2671 = _S2474 * _S2670;
    float3  _S2672 = _S2453.rows[0U] * _S2670;
    float _S2673 = _S2458 * (_S2672.x + _S2672.y + _S2672.z);
    Matrix<float, 3, 3>  _S2674 = _S2587;
    _S2674[2U] = _S2663;
    _S2674[1U] = _S2667;
    _S2674[0U] = _S2671;
    Matrix<float, 3, 3>  _S2675 = transpose_0(transpose_0(_S2674));
    float _S2676 = 2.0f * - _S2675.rows[int(2)].z;
    float _S2677 = 2.0f * _S2675.rows[int(2)].y;
    float _S2678 = 2.0f * _S2675.rows[int(2)].x;
    float _S2679 = 2.0f * _S2675.rows[int(1)].z;
    float _S2680 = 2.0f * - _S2675.rows[int(1)].y;
    float _S2681 = 2.0f * _S2675.rows[int(1)].x;
    float _S2682 = 2.0f * _S2675.rows[int(0)].z;
    float _S2683 = 2.0f * _S2675.rows[int(0)].y;
    float _S2684 = 2.0f * - _S2675.rows[int(0)].x;
    float _S2685 = - _S2681 + _S2683;
    float _S2686 = _S2678 + - _S2682;
    float _S2687 = - _S2677 + _S2679;
    float _S2688 = _S2677 + _S2679;
    float _S2689 = _S2678 + _S2682;
    float _S2690 = _S2681 + _S2683;
    float _S2691 = quat_16.w * (_S2680 + _S2684);
    float _S2692 = quat_16.z * (_S2676 + _S2684);
    float _S2693 = quat_16.y * (_S2676 + _S2680);
    float _S2694 = quat_16.x * _S2685 + quat_16.z * _S2688 + quat_16.y * _S2689 + _S2691 + _S2691;
    float _S2695 = quat_16.x * _S2686 + quat_16.w * _S2688 + quat_16.y * _S2690 + _S2692 + _S2692;
    float _S2696 = quat_16.x * _S2687 + quat_16.w * _S2689 + quat_16.z * _S2690 + _S2693 + _S2693;
    float _S2697 = quat_16.w * _S2685 + quat_16.z * _S2686 + quat_16.y * _S2687;
    float3  _S2698 = _S2450;
    *&((&_S2698)->z) = _S2665;
    *&((&_S2698)->y) = _S2669;
    *&((&_S2698)->x) = _S2673;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2699;
    (&_S2699)->primal_0 = scale_15;
    (&_S2699)->differential_0 = _S2450;
    s_bwd_prop_exp_1(&_S2699, _S2698);
    Matrix<float, 2, 2>  _S2700 = _S2447;
    _S2700[int(1)] = _S2645;
    _S2700[int(0)] = _S2647;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2701;
    (&_S2701)->primal_0 = R_17;
    (&_S2701)->differential_0 = _S2587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2702;
    (&_S2702)->primal_0 = mean_13;
    (&_S2702)->differential_0 = _S2450;
    s_bwd_prop_mul_1(&_S2701, &_S2702, _S2623);
    float3  _S2703 = _S2623 + _S2590.differential_0;
    Matrix<float, 3, 3>  _S2704 = _S2648.differential_0 + _S2650.differential_0 + _S2652.differential_0 + _S2654.differential_0 + _S2656.differential_0 + _S2658.differential_0 + _S2660.differential_0 + _S2701.differential_0 + _S2591;
    float3  _S2705 = _S2699.differential_0 + _S2593;
    float4  _S2706 = make_float4 (0.0f);
    *&((&_S2706)->w) = _S2694;
    *&((&_S2706)->z) = _S2695;
    *&((&_S2706)->y) = _S2696;
    *&((&_S2706)->x) = _S2697;
    float4  _S2707 = _S2706;
    float3  _S2708 = _S2649.differential_0 + _S2655.differential_0 + _S2651.differential_0 + _S2657.differential_0 + _S2653.differential_0 + _S2659.differential_0 + _S2702.differential_0 + _S2585;
    *v_mean_3 = _S2708;
    *v_quat_3 = _S2707;
    *v_scale_3 = _S2705;
    *v_in_opacity_3 = _S2630;
    (*v_sh_coeffs_3)[int(0)] = _S2607;
    (*v_sh_coeffs_3)[int(1)] = _S2608;
    (*v_sh_coeffs_3)[int(2)] = _S2609;
    (*v_sh_coeffs_3)[int(3)] = _S2610;
    (*v_sh_coeffs_3)[int(4)] = _S2611;
    (*v_sh_coeffs_3)[int(5)] = _S2612;
    (*v_sh_coeffs_3)[int(6)] = _S2613;
    (*v_sh_coeffs_3)[int(7)] = _S2614;
    (*v_sh_coeffs_3)[int(8)] = _S2615;
    (*v_sh_coeffs_3)[int(9)] = _S2616;
    (*v_sh_coeffs_3)[int(10)] = _S2617;
    (*v_sh_coeffs_3)[int(11)] = _S2618;
    (*v_sh_coeffs_3)[int(12)] = _S2619;
    (*v_sh_coeffs_3)[int(13)] = _S2620;
    (*v_sh_coeffs_3)[int(14)] = _S2621;
    (*v_sh_coeffs_3)[int(15)] = _S2622;
    *v_R_4 = _S2704;
    *v_t_4 = _S2703;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2709;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_23, float fy_23, float cx_18, float cy_18, FixedArray<float, 10>  * dist_coeffs_29, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2710 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2711;
    (&_S2711)->_S2709 = _S2710;
    float2  _S2712 = make_float2 (0.0f);
    float3  _S2713 = make_float3 (0.0f);
    float4  intrins_13 = make_float4 (fx_23, fy_23, cx_18, cy_18);
    float3  _S2714 = s_primal_ctx_exp_0(scale_16);
    float _S2715 = quat_17.y;
    float x2_17 = _S2715 * _S2715;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_31 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2716 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_31), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_31), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17))));
    FixedArray<float3 , 7>  _S2717 = {
        _S2713, _S2713, _S2713, _S2713, _S2713, _S2713, _S2713
    };
    FixedArray<float, 7>  _S2718 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2719;
    (&_S2719)->p_0 = _S2717;
    (&_S2719)->w_mean_0 = _S2718;
    (&_S2719)->w_cov_0 = _S2718;
    (&_S2719)->p_0[int(0)] = mean_14;
    SigmaPoints_0 _S2720 = _S2719;
    (&_S2720)->w_mean_0[int(0)] = 0.0f;
    (&_S2720)->w_cov_0[int(0)] = 2.0f;
    float _S2721 = s_primal_ctx_sqrt_0(3.0f);
    float _S2722 = _S2721 * _S2714.x;
    float3  delta_15 = make_float3 (_S2722) * _S2716.rows[0U];
    float3  _S2723 = mean_14 + delta_15;
    float3  _S2724 = mean_14 - delta_15;
    float _S2725 = _S2721 * _S2714.y;
    float3  delta_16 = make_float3 (_S2725) * _S2716.rows[1U];
    float3  _S2726 = mean_14 + delta_16;
    float3  _S2727 = mean_14 - delta_16;
    float _S2728 = _S2721 * _S2714.z;
    float3  delta_17 = make_float3 (_S2728) * _S2716.rows[2U];
    float3  _S2729 = mean_14 + delta_17;
    float3  _S2730 = mean_14 - delta_17;
    (&_S2720)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2720)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2720)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2720)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2720)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2720)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2720)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2720)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2720)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2720)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2720)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2720)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2731 = _S2719;
    (&_S2720)->p_0[0U] = s_primal_ctx_mul_1(R_18, _S2719.p_0[0U]) + t_17;
    (&_S2720)->p_0[1U] = s_primal_ctx_mul_1(R_18, _S2723) + t_17;
    (&_S2720)->p_0[2U] = s_primal_ctx_mul_1(R_18, _S2726) + t_17;
    (&_S2720)->p_0[3U] = s_primal_ctx_mul_1(R_18, _S2729) + t_17;
    (&_S2720)->p_0[4U] = s_primal_ctx_mul_1(R_18, _S2724) + t_17;
    (&_S2720)->p_0[5U] = s_primal_ctx_mul_1(R_18, _S2727) + t_17;
    (&_S2720)->p_0[6U] = s_primal_ctx_mul_1(R_18, _S2730) + t_17;
    float2  _S2732 = _S2712;
    Matrix<float, 2, 2>  _S2733 = _S2710;
    SigmaPoints_0 _S2734 = _S2720;
    bool _S2735 = fisheye_proj_3dgs_ut_1(&_S2734, intrins_13, dist_coeffs_29, &_S2733, &_S2732);
    (&_S2711)->_S2709 = _S2733;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2736 = _S2711;
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_14) + t_17;
    float3  _S2737 = make_float3 (_S2722);
    float3  _S2738 = make_float3 (_S2725);
    float3  _S2739 = make_float3 (_S2728);
    float _S2740 = _S2711._S2709.rows[int(0)].y * _S2711._S2709.rows[int(1)].x;
    float det_orig_15 = _S2711._S2709.rows[int(0)].x * _S2711._S2709.rows[int(1)].y - _S2740;
    float _S2741 = _S2711._S2709.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2742 = _S2711._S2709;
    *&(((&_S2742)->rows + (int(0)))->x) = _S2741;
    float _S2743 = _S2711._S2709.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2742)->rows + (int(1)))->y) = _S2743;
    Matrix<float, 2, 2>  _S2744 = _S2742;
    Matrix<float, 2, 2>  _S2745 = _S2742;
    float det_blur_10 = _S2741 * _S2743 - _S2740;
    float _S2746 = det_orig_15 / det_blur_10;
    float _S2747 = det_blur_10 * det_blur_10;
    float _S2748 = s_primal_ctx_max_0(0.0f, _S2746);
    float _S2749 = s_primal_ctx_sqrt_0(_S2748);
    float _S2750 = - in_opacity_14;
    float _S2751 = 1.0f + s_primal_ctx_exp_1(_S2750);
    float _S2752 = 1.0f / _S2751;
    float _S2753 = _S2751 * _S2751;
    float _S2754;
    if(antialiased_14)
    {
        _S2754 = _S2752 * _S2749;
    }
    else
    {
        _S2754 = _S2752;
    }
    float _S2755 = _S2754 / 0.00392156885936856f;
    float _S2756 = 2.0f * s_primal_ctx_log_0(_S2755);
    float _S2757 = s_primal_ctx_sqrt_0(_S2756);
    float _S2758 = _S2744.rows[int(0)].x;
    float _S2759 = _S2745.rows[int(1)].y;
    float _S2760 = s_primal_ctx_dot_0(mean_c_14, mean_c_14) + 9.99999997475242708e-07f;
    float3  _S2761 = - scale_16;
    Matrix<float, 3, 3>  _S2762 = transpose_0(R_18);
    float3  _S2763 = mean_14 - - s_primal_ctx_mul_1(_S2762, t_17);
    float _S2764 = _S2763.x;
    float _S2765 = _S2763.y;
    float _S2766 = _S2763.z;
    float _S2767 = _S2764 * _S2764 + _S2765 * _S2765 + _S2766 * _S2766;
    float _S2768 = s_primal_ctx_sqrt_0(_S2767);
    float x_39 = _S2764 / _S2768;
    float3  _S2769 = make_float3 (x_39);
    float _S2770 = _S2768 * _S2768;
    float y_17 = _S2765 / _S2768;
    float z_14 = _S2766 / _S2768;
    float3  _S2771 = make_float3 (z_14);
    float _S2772 = - y_17;
    float3  _S2773 = make_float3 (_S2772);
    float z2_32 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_39 * x_39 - y_17 * y_17;
    float _S2774 = 2.0f * x_39;
    float fS1_14 = _S2774 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_32 - 0.31539157032966614f;
    float3  _S2775 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_39;
    float3  _S2776 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2777 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2778 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2779 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_32 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2780 = 1.86588168144226074f * z2_32 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2780;
    float3  _S2781 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_39;
    float3  _S2782 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2783 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2784 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2785 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_39 * fC1_14 - y_17 * fS1_14);
    float3  _S2786 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_39 * fS1_14 + y_17 * fC1_14);
    float3  _S2787 = make_float3 (pSH9_4);
    float3  _S2788 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2789;
    (&_S2789)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2772) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_39) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2789)->differential_0 = _S2713;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2790;
    (&_S2790)->primal_0 = _S2788;
    (&_S2790)->differential_0 = _S2713;
    s_bwd_prop_max_0(&_S2789, &_S2790, v_rgb_4);
    float3  _S2791 = _S2786 * _S2789.differential_0;
    float3  _S2792 = (*sh_coeffs_14)[int(15)] * _S2789.differential_0;
    float3  _S2793 = _S2784 * _S2789.differential_0;
    float3  _S2794 = (*sh_coeffs_14)[int(14)] * _S2789.differential_0;
    float3  _S2795 = _S2782 * _S2789.differential_0;
    float3  _S2796 = (*sh_coeffs_14)[int(13)] * _S2789.differential_0;
    float3  _S2797 = _S2781 * _S2789.differential_0;
    float3  _S2798 = (*sh_coeffs_14)[int(12)] * _S2789.differential_0;
    float3  _S2799 = _S2783 * _S2789.differential_0;
    float3  _S2800 = (*sh_coeffs_14)[int(11)] * _S2789.differential_0;
    float3  _S2801 = _S2785 * _S2789.differential_0;
    float3  _S2802 = (*sh_coeffs_14)[int(10)] * _S2789.differential_0;
    float3  _S2803 = _S2787 * _S2789.differential_0;
    float3  _S2804 = (*sh_coeffs_14)[int(9)] * _S2789.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2804.x + _S2804.y + _S2804.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2792.x + _S2792.y + _S2792.z);
    float _S2805 = _S2802.x + _S2802.y + _S2802.z;
    float _S2806 = _S2794.x + _S2794.y + _S2794.z;
    float _S2807 = _S2800.x + _S2800.y + _S2800.z;
    float _S2808 = _S2796.x + _S2796.y + _S2796.z;
    float _S2809 = _S2798.x + _S2798.y + _S2798.z;
    float _S2810 = - s_diff_fC2_T_4;
    float3  _S2811 = _S2778 * _S2789.differential_0;
    float3  _S2812 = (*sh_coeffs_14)[int(8)] * _S2789.differential_0;
    float3  _S2813 = _S2776 * _S2789.differential_0;
    float3  _S2814 = (*sh_coeffs_14)[int(7)] * _S2789.differential_0;
    float3  _S2815 = _S2775 * _S2789.differential_0;
    float3  _S2816 = (*sh_coeffs_14)[int(6)] * _S2789.differential_0;
    float3  _S2817 = _S2777 * _S2789.differential_0;
    float3  _S2818 = (*sh_coeffs_14)[int(5)] * _S2789.differential_0;
    float3  _S2819 = _S2779 * _S2789.differential_0;
    float3  _S2820 = (*sh_coeffs_14)[int(4)] * _S2789.differential_0;
    float _S2821 = _S2818.x + _S2818.y + _S2818.z;
    float _S2822 = _S2814.x + _S2814.y + _S2814.z;
    float _S2823 = fTmp1B_14 * _S2805 + x_39 * s_diff_fS2_T_4 + y_17 * _S2810 + 0.54627424478530884f * (_S2820.x + _S2820.y + _S2820.z);
    float _S2824 = fTmp1B_14 * _S2806 + y_17 * s_diff_fS2_T_4 + x_39 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2812.x + _S2812.y + _S2812.z);
    float _S2825 = y_17 * - _S2824;
    float _S2826 = x_39 * _S2824;
    float _S2827 = z_14 * (1.86588168144226074f * (z_14 * _S2809) + -2.28522896766662598f * (y_17 * _S2807 + x_39 * _S2808) + 0.94617468118667603f * (_S2816.x + _S2816.y + _S2816.z));
    float3  _S2828 = make_float3 (0.48860251903533936f) * _S2789.differential_0;
    float3  _S2829 = - _S2828;
    float3  _S2830 = _S2769 * _S2829;
    float3  _S2831 = (*sh_coeffs_14)[int(3)] * _S2829;
    float3  _S2832 = _S2771 * _S2828;
    float3  _S2833 = (*sh_coeffs_14)[int(2)] * _S2828;
    float3  _S2834 = _S2773 * _S2828;
    float3  _S2835 = (*sh_coeffs_14)[int(1)] * _S2828;
    float _S2836 = (_S2780 * _S2809 + 1.44530570507049561f * (fS1_14 * _S2805 + fC1_14 * _S2806) + -1.09254848957061768f * (y_17 * _S2821 + x_39 * _S2822) + _S2827 + _S2827 + _S2833.x + _S2833.y + _S2833.z) / _S2770;
    float _S2837 = _S2768 * _S2836;
    float _S2838 = (fTmp0C_14 * _S2807 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2810 + fTmp0B_14 * _S2821 + _S2774 * _S2823 + _S2825 + _S2825 + - (_S2835.x + _S2835.y + _S2835.z)) / _S2770;
    float _S2839 = _S2768 * _S2838;
    float _S2840 = (fTmp0C_14 * _S2808 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2822 + 2.0f * (y_17 * _S2823) + _S2826 + _S2826 + _S2831.x + _S2831.y + _S2831.z) / _S2770;
    float _S2841 = _S2768 * _S2840;
    float _S2842 = _S2766 * - _S2836 + _S2765 * - _S2838 + _S2764 * - _S2840;
    DiffPair_float_0 _S2843;
    (&_S2843)->primal_0 = _S2767;
    (&_S2843)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2843, _S2842);
    float _S2844 = _S2766 * _S2843.differential_0;
    float _S2845 = _S2765 * _S2843.differential_0;
    float _S2846 = _S2764 * _S2843.differential_0;
    float3  _S2847 = make_float3 (0.282094806432724f) * _S2789.differential_0;
    float3  _S2848 = make_float3 (_S2841 + _S2846 + _S2846, _S2839 + _S2845 + _S2845, _S2837 + _S2844 + _S2844);
    float3  _S2849 = - - _S2848;
    Matrix<float, 3, 3>  _S2850 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2851;
    (&_S2851)->primal_0 = _S2762;
    (&_S2851)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2852;
    (&_S2852)->primal_0 = t_17;
    (&_S2852)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2851, &_S2852, _S2849);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2853 = _S2852;
    Matrix<float, 3, 3>  _S2854 = transpose_0(_S2851.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2855;
    (&_S2855)->primal_0 = _S2761;
    (&_S2855)->differential_0 = _S2713;
    s_bwd_prop_exp_1(&_S2855, v_conic_4);
    float3  _S2856 = - _S2855.differential_0;
    float _S2857 = 0.5f * v_depth_4;
    DiffPair_float_0 _S2858;
    (&_S2858)->primal_0 = _S2760;
    (&_S2858)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2858, _S2857);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2859;
    (&_S2859)->primal_0 = mean_c_14;
    (&_S2859)->differential_0 = _S2713;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2860;
    (&_S2860)->primal_0 = mean_c_14;
    (&_S2860)->differential_0 = _S2713;
    s_bwd_prop_dot_0(&_S2859, &_S2860, _S2858.differential_0);
    DiffPair_float_0 _S2861;
    (&_S2861)->primal_0 = _S2759;
    (&_S2861)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2861, 0.0f);
    DiffPair_float_0 _S2862;
    (&_S2862)->primal_0 = _S2758;
    (&_S2862)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2862, 0.0f);
    DiffPair_float_0 _S2863;
    (&_S2863)->primal_0 = 3.32999992370605469f;
    (&_S2863)->differential_0 = 0.0f;
    DiffPair_float_0 _S2864;
    (&_S2864)->primal_0 = _S2757;
    (&_S2864)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2863, &_S2864, 0.0f);
    DiffPair_float_0 _S2865;
    (&_S2865)->primal_0 = _S2756;
    (&_S2865)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2865, _S2864.differential_0);
    float _S2866 = 2.0f * _S2865.differential_0;
    DiffPair_float_0 _S2867;
    (&_S2867)->primal_0 = _S2755;
    (&_S2867)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2867, _S2866);
    float _S2868 = v_opacity_4 + 254.9999847412109375f * _S2867.differential_0;
    FixedArray<float3 , 16>  _S2869;
    _S2869[int(0)] = _S2713;
    _S2869[int(1)] = _S2713;
    _S2869[int(2)] = _S2713;
    _S2869[int(3)] = _S2713;
    _S2869[int(4)] = _S2713;
    _S2869[int(5)] = _S2713;
    _S2869[int(6)] = _S2713;
    _S2869[int(7)] = _S2713;
    _S2869[int(8)] = _S2713;
    _S2869[int(9)] = _S2713;
    _S2869[int(10)] = _S2713;
    _S2869[int(11)] = _S2713;
    _S2869[int(12)] = _S2713;
    _S2869[int(13)] = _S2713;
    _S2869[int(14)] = _S2713;
    _S2869[int(15)] = _S2713;
    _S2869[int(7)] = _S2813;
    _S2869[int(0)] = _S2847;
    _S2869[int(1)] = _S2834;
    _S2869[int(2)] = _S2832;
    _S2869[int(3)] = _S2830;
    _S2869[int(4)] = _S2819;
    _S2869[int(5)] = _S2817;
    _S2869[int(6)] = _S2815;
    _S2869[int(15)] = _S2791;
    _S2869[int(8)] = _S2811;
    _S2869[int(9)] = _S2803;
    _S2869[int(10)] = _S2801;
    _S2869[int(11)] = _S2799;
    _S2869[int(12)] = _S2797;
    _S2869[int(13)] = _S2795;
    _S2869[int(14)] = _S2793;
    float3  _S2870 = _S2869[int(0)];
    float3  _S2871 = _S2869[int(1)];
    float3  _S2872 = _S2869[int(2)];
    float3  _S2873 = _S2869[int(3)];
    float3  _S2874 = _S2869[int(4)];
    float3  _S2875 = _S2869[int(5)];
    float3  _S2876 = _S2869[int(6)];
    float3  _S2877 = _S2869[int(7)];
    float3  _S2878 = _S2869[int(8)];
    float3  _S2879 = _S2869[int(9)];
    float3  _S2880 = _S2869[int(10)];
    float3  _S2881 = _S2869[int(11)];
    float3  _S2882 = _S2869[int(12)];
    float3  _S2883 = _S2869[int(13)];
    float3  _S2884 = _S2869[int(14)];
    float3  _S2885 = _S2869[int(15)];
    float3  _S2886 = _S2860.differential_0 + _S2859.differential_0;
    float2  _S2887 = make_float2 (0.0f, _S2861.differential_0);
    float2  _S2888 = make_float2 (_S2862.differential_0, 0.0f);
    float _S2889;
    if(antialiased_14)
    {
        float _S2890 = _S2752 * _S2868;
        _S2754 = _S2749 * _S2868;
        _S2889 = _S2890;
    }
    else
    {
        _S2754 = _S2868;
        _S2889 = 0.0f;
    }
    float _S2891 = - (_S2754 / _S2753);
    DiffPair_float_0 _S2892;
    (&_S2892)->primal_0 = _S2750;
    (&_S2892)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2892, _S2891);
    float _S2893 = - _S2892.differential_0;
    DiffPair_float_0 _S2894;
    (&_S2894)->primal_0 = _S2748;
    (&_S2894)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2894, _S2889);
    DiffPair_float_0 _S2895;
    (&_S2895)->primal_0 = 0.0f;
    (&_S2895)->differential_0 = 0.0f;
    DiffPair_float_0 _S2896;
    (&_S2896)->primal_0 = _S2746;
    (&_S2896)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2895, &_S2896, _S2894.differential_0);
    float _S2897 = _S2896.differential_0 / _S2747;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2897;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2897;
    float _S2898 = - s_diff_det_blur_T_1;
    float _S2899 = _S2741 * s_diff_det_blur_T_1;
    float _S2900 = _S2743 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2901 = _S2710;
    _S2901[int(1)] = _S2887;
    _S2901[int(0)] = _S2888;
    float _S2902 = _S2900 + _S2901.rows[int(0)].x;
    float _S2903 = _S2898 + - s_diff_det_orig_T_4;
    float _S2904 = _S2736._S2709.rows[int(0)].y * _S2903;
    float _S2905 = _S2736._S2709.rows[int(1)].x * _S2903;
    float _S2906 = _S2736._S2709.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2907 = _S2899 + _S2901.rows[int(1)].y + _S2736._S2709.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2908 = _S2712;
    *&((&_S2908)->x) = _S2904;
    *&((&_S2908)->y) = _S2907;
    float _S2909 = _S2902 + _S2906;
    float2  _S2910 = _S2712;
    *&((&_S2910)->y) = _S2905;
    *&((&_S2910)->x) = _S2909;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2911;
    (&_S2911)->primal_0 = R_18;
    (&_S2911)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2912;
    (&_S2912)->primal_0 = _S2730;
    (&_S2912)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2911, &_S2912, _S2713);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2913;
    (&_S2913)->primal_0 = R_18;
    (&_S2913)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2914;
    (&_S2914)->primal_0 = _S2727;
    (&_S2914)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2913, &_S2914, _S2713);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2915;
    (&_S2915)->primal_0 = R_18;
    (&_S2915)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2916;
    (&_S2916)->primal_0 = _S2724;
    (&_S2916)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2915, &_S2916, _S2713);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2917;
    (&_S2917)->primal_0 = R_18;
    (&_S2917)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2918;
    (&_S2918)->primal_0 = _S2729;
    (&_S2918)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2917, &_S2918, _S2713);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2919;
    (&_S2919)->primal_0 = R_18;
    (&_S2919)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2920;
    (&_S2920)->primal_0 = _S2726;
    (&_S2920)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2919, &_S2920, _S2713);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2921;
    (&_S2921)->primal_0 = R_18;
    (&_S2921)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2922;
    (&_S2922)->primal_0 = _S2723;
    (&_S2922)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2921, &_S2922, _S2713);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2923;
    (&_S2923)->primal_0 = R_18;
    (&_S2923)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2924;
    (&_S2924)->primal_0 = _S2731.p_0[0U];
    (&_S2924)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2923, &_S2924, _S2713);
    float3  _S2925 = - _S2912.differential_0 + _S2918.differential_0;
    float3  _S2926 = _S2739 * _S2925;
    float3  _S2927 = _S2716.rows[2U] * _S2925;
    float _S2928 = _S2721 * (_S2927.x + _S2927.y + _S2927.z);
    float3  _S2929 = - _S2914.differential_0 + _S2920.differential_0;
    float3  _S2930 = _S2738 * _S2929;
    float3  _S2931 = _S2716.rows[1U] * _S2929;
    float _S2932 = _S2721 * (_S2931.x + _S2931.y + _S2931.z);
    float3  _S2933 = - _S2916.differential_0 + _S2922.differential_0;
    float3  _S2934 = _S2737 * _S2933;
    float3  _S2935 = _S2716.rows[0U] * _S2933;
    float _S2936 = _S2721 * (_S2935.x + _S2935.y + _S2935.z);
    Matrix<float, 3, 3>  _S2937 = _S2850;
    _S2937[2U] = _S2926;
    _S2937[1U] = _S2930;
    _S2937[0U] = _S2934;
    Matrix<float, 3, 3>  _S2938 = transpose_0(transpose_0(_S2937));
    float _S2939 = 2.0f * - _S2938.rows[int(2)].z;
    float _S2940 = 2.0f * _S2938.rows[int(2)].y;
    float _S2941 = 2.0f * _S2938.rows[int(2)].x;
    float _S2942 = 2.0f * _S2938.rows[int(1)].z;
    float _S2943 = 2.0f * - _S2938.rows[int(1)].y;
    float _S2944 = 2.0f * _S2938.rows[int(1)].x;
    float _S2945 = 2.0f * _S2938.rows[int(0)].z;
    float _S2946 = 2.0f * _S2938.rows[int(0)].y;
    float _S2947 = 2.0f * - _S2938.rows[int(0)].x;
    float _S2948 = - _S2944 + _S2946;
    float _S2949 = _S2941 + - _S2945;
    float _S2950 = - _S2940 + _S2942;
    float _S2951 = _S2940 + _S2942;
    float _S2952 = _S2941 + _S2945;
    float _S2953 = _S2944 + _S2946;
    float _S2954 = quat_17.w * (_S2943 + _S2947);
    float _S2955 = quat_17.z * (_S2939 + _S2947);
    float _S2956 = quat_17.y * (_S2939 + _S2943);
    float _S2957 = quat_17.x * _S2948 + quat_17.z * _S2951 + quat_17.y * _S2952 + _S2954 + _S2954;
    float _S2958 = quat_17.x * _S2949 + quat_17.w * _S2951 + quat_17.y * _S2953 + _S2955 + _S2955;
    float _S2959 = quat_17.x * _S2950 + quat_17.w * _S2952 + quat_17.z * _S2953 + _S2956 + _S2956;
    float _S2960 = quat_17.w * _S2948 + quat_17.z * _S2949 + quat_17.y * _S2950;
    float3  _S2961 = _S2713;
    *&((&_S2961)->z) = _S2928;
    *&((&_S2961)->y) = _S2932;
    *&((&_S2961)->x) = _S2936;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2962;
    (&_S2962)->primal_0 = scale_16;
    (&_S2962)->differential_0 = _S2713;
    s_bwd_prop_exp_1(&_S2962, _S2961);
    Matrix<float, 2, 2>  _S2963 = _S2710;
    _S2963[int(1)] = _S2908;
    _S2963[int(0)] = _S2910;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2964;
    (&_S2964)->primal_0 = mean_c_14;
    (&_S2964)->differential_0 = _S2713;
    s_bwd_length_impl_0(&_S2964, 0.0f);
    float3  _S2965 = _S2964.differential_0 + _S2886;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2966;
    (&_S2966)->primal_0 = R_18;
    (&_S2966)->differential_0 = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2967;
    (&_S2967)->primal_0 = mean_14;
    (&_S2967)->differential_0 = _S2713;
    s_bwd_prop_mul_1(&_S2966, &_S2967, _S2965);
    float3  _S2968 = _S2965 + _S2853.differential_0;
    Matrix<float, 3, 3>  _S2969 = _S2911.differential_0 + _S2913.differential_0 + _S2915.differential_0 + _S2917.differential_0 + _S2919.differential_0 + _S2921.differential_0 + _S2923.differential_0 + _S2966.differential_0 + _S2854;
    float3  _S2970 = _S2962.differential_0 + _S2856;
    float4  _S2971 = make_float4 (0.0f);
    *&((&_S2971)->w) = _S2957;
    *&((&_S2971)->z) = _S2958;
    *&((&_S2971)->y) = _S2959;
    *&((&_S2971)->x) = _S2960;
    float4  _S2972 = _S2971;
    float3  _S2973 = _S2912.differential_0 + _S2918.differential_0 + _S2914.differential_0 + _S2920.differential_0 + _S2916.differential_0 + _S2922.differential_0 + _S2967.differential_0 + _S2848;
    *v_mean_4 = _S2973;
    *v_quat_4 = _S2972;
    *v_scale_4 = _S2970;
    *v_in_opacity_4 = _S2893;
    (*v_sh_coeffs_4)[int(0)] = _S2870;
    (*v_sh_coeffs_4)[int(1)] = _S2871;
    (*v_sh_coeffs_4)[int(2)] = _S2872;
    (*v_sh_coeffs_4)[int(3)] = _S2873;
    (*v_sh_coeffs_4)[int(4)] = _S2874;
    (*v_sh_coeffs_4)[int(5)] = _S2875;
    (*v_sh_coeffs_4)[int(6)] = _S2876;
    (*v_sh_coeffs_4)[int(7)] = _S2877;
    (*v_sh_coeffs_4)[int(8)] = _S2878;
    (*v_sh_coeffs_4)[int(9)] = _S2879;
    (*v_sh_coeffs_4)[int(10)] = _S2880;
    (*v_sh_coeffs_4)[int(11)] = _S2881;
    (*v_sh_coeffs_4)[int(12)] = _S2882;
    (*v_sh_coeffs_4)[int(13)] = _S2883;
    (*v_sh_coeffs_4)[int(14)] = _S2884;
    (*v_sh_coeffs_4)[int(15)] = _S2885;
    *v_R_5 = _S2969;
    *v_t_5 = _S2968;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_18, float3  scale_17)
{
    float x_40 = quat_18.y;
    float x2_18 = x_40 * x_40;
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

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_6)
{
    float _S2974 = (*dpquat_0).primal_0.y;
    float x2_19 = _S2974 * _S2974;
    float y2_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_34 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_19 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_19 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_19 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2975 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_34), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_34), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19))));
    Matrix<float, 3, 3>  _S2976 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2977;
    (&_S2977)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2977)->differential_0 = _S2976;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2978;
    (&_S2978)->primal_0 = _S2975;
    (&_S2978)->differential_0 = _S2976;
    s_bwd_prop_mul_4(&_S2977, &_S2978, _s_dOut_6);
    Matrix<float, 3, 3>  _S2979 = transpose_0(transpose_0(_S2978.differential_0));
    float _S2980 = 2.0f * - _S2979.rows[int(2)].z;
    float _S2981 = 2.0f * _S2979.rows[int(2)].y;
    float _S2982 = 2.0f * _S2979.rows[int(2)].x;
    float _S2983 = 2.0f * _S2979.rows[int(1)].z;
    float _S2984 = 2.0f * - _S2979.rows[int(1)].y;
    float _S2985 = 2.0f * _S2979.rows[int(1)].x;
    float _S2986 = 2.0f * _S2979.rows[int(0)].z;
    float _S2987 = 2.0f * _S2979.rows[int(0)].y;
    float _S2988 = 2.0f * - _S2979.rows[int(0)].x;
    float _S2989 = - _S2985 + _S2987;
    float _S2990 = _S2982 + - _S2986;
    float _S2991 = - _S2981 + _S2983;
    float _S2992 = _S2981 + _S2983;
    float _S2993 = _S2982 + _S2986;
    float _S2994 = _S2985 + _S2987;
    float _S2995 = (*dpquat_0).primal_0.w * (_S2984 + _S2988);
    float _S2996 = (*dpquat_0).primal_0.z * (_S2980 + _S2988);
    float _S2997 = (*dpquat_0).primal_0.y * (_S2980 + _S2984);
    float _S2998 = (*dpquat_0).primal_0.x * _S2989 + (*dpquat_0).primal_0.z * _S2992 + (*dpquat_0).primal_0.y * _S2993 + _S2995 + _S2995;
    float _S2999 = (*dpquat_0).primal_0.x * _S2990 + (*dpquat_0).primal_0.w * _S2992 + (*dpquat_0).primal_0.y * _S2994 + _S2996 + _S2996;
    float _S3000 = (*dpquat_0).primal_0.x * _S2991 + (*dpquat_0).primal_0.w * _S2993 + (*dpquat_0).primal_0.z * _S2994 + _S2997 + _S2997;
    float _S3001 = (*dpquat_0).primal_0.w * _S2989 + (*dpquat_0).primal_0.z * _S2990 + (*dpquat_0).primal_0.y * _S2991;
    float3  _S3002 = make_float3 (_S2977.differential_0.rows[int(0)].x, _S2977.differential_0.rows[int(1)].y, _S2977.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S3002;
    float4  _S3003 = make_float4 (0.0f);
    *&((&_S3003)->w) = _S2998;
    *&((&_S3003)->z) = _S2999;
    *&((&_S3003)->y) = _S3000;
    *&((&_S3003)->x) = _S3001;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S3003;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S3004, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3005, Matrix<float, 3, 3>  _S3006)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S3004, _S3005, _S3006);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_19, float3  scale_18, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S3007 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_19;
    (&dp_quat_0)->differential_0 = _S3007;
    float3  _S3008 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_18;
    (&dp_scale_0)->differential_0 = _S3008;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_16)
{
    float _S3009 = dOut_16.y;
    float _S3010 = dOut_16.z;
    float _S3011 = dOut_16.x;
    float _S3012 = (*a_0).primal_0.z * _S3009 + - (*a_0).primal_0.y * _S3010;
    float _S3013 = - (*a_0).primal_0.z * _S3011 + (*a_0).primal_0.x * _S3010;
    float _S3014 = (*a_0).primal_0.y * _S3011 + - (*a_0).primal_0.x * _S3009;
    float3  _S3015 = make_float3 (- (*b_0).primal_0.z * _S3009 + (*b_0).primal_0.y * _S3010, (*b_0).primal_0.z * _S3011 + - (*b_0).primal_0.x * _S3010, - (*b_0).primal_0.y * _S3011 + (*b_0).primal_0.x * _S3009);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S3015;
    float3  _S3016 = make_float3 (_S3012, _S3013, _S3014);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S3016;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S3017 = left_10.y;
    float _S3018 = right_10.z;
    float _S3019 = left_10.z;
    float _S3020 = right_10.y;
    float _S3021 = right_10.x;
    float _S3022 = left_10.x;
    return make_float3 (_S3017 * _S3018 - _S3019 * _S3020, _S3019 * _S3021 - _S3022 * _S3018, _S3022 * _S3020 - _S3017 * _S3021);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_15, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_1, float3  ray_d_1)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_1);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_1 - mean_15));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S3023, float3  _S3024)
{
    return cross_0(_S3023, _S3024);
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3025, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3026, float3  _S3027)
{
    _d_cross_0(_S3025, _S3026, _S3027);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_7)
{
    float3  _S3028 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S3029 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, _S3028);
    float3  _S3030 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, (*dpray_d_2).primal_0);
    float3  _S3031 = s_primal_ctx_cross_0(_S3030, _S3029);
    float _S3032 = -0.5f * s_primal_ctx_dot_0(_S3031, _S3031);
    float _S3033 = s_primal_ctx_dot_0(_S3030, _S3030);
    float _S3034 = _S3032 / _S3033;
    float _S3035 = _S3033 * _S3033;
    float _S3036 = (*dpopacity_0).primal_0 * _s_dOut_7;
    float _S3037 = s_primal_ctx_exp_1(_S3034) * _s_dOut_7;
    DiffPair_float_0 _S3038;
    (&_S3038)->primal_0 = _S3034;
    (&_S3038)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3038, _S3036);
    float _S3039 = _S3038.differential_0 / _S3035;
    float _S3040 = _S3032 * - _S3039;
    float _S3041 = _S3033 * _S3039;
    float3  _S3042 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3043;
    (&_S3043)->primal_0 = _S3030;
    (&_S3043)->differential_0 = _S3042;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3044;
    (&_S3044)->primal_0 = _S3030;
    (&_S3044)->differential_0 = _S3042;
    s_bwd_prop_dot_0(&_S3043, &_S3044, _S3040);
    float _S3045 = -0.5f * _S3041;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3046;
    (&_S3046)->primal_0 = _S3031;
    (&_S3046)->differential_0 = _S3042;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3047;
    (&_S3047)->primal_0 = _S3031;
    (&_S3047)->differential_0 = _S3042;
    s_bwd_prop_dot_0(&_S3046, &_S3047, _S3045);
    float3  _S3048 = _S3047.differential_0 + _S3046.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3049;
    (&_S3049)->primal_0 = _S3030;
    (&_S3049)->differential_0 = _S3042;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3050;
    (&_S3050)->primal_0 = _S3029;
    (&_S3050)->differential_0 = _S3042;
    s_bwd_prop_cross_0(&_S3049, &_S3050, _S3048);
    float3  _S3051 = _S3044.differential_0 + _S3043.differential_0 + _S3049.differential_0;
    Matrix<float, 3, 3>  _S3052 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3053;
    (&_S3053)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3053)->differential_0 = _S3052;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3054;
    (&_S3054)->primal_0 = (*dpray_d_2).primal_0;
    (&_S3054)->differential_0 = _S3042;
    s_bwd_prop_mul_1(&_S3053, &_S3054, _S3051);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3055;
    (&_S3055)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3055)->differential_0 = _S3052;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3056;
    (&_S3056)->primal_0 = _S3028;
    (&_S3056)->differential_0 = _S3042;
    s_bwd_prop_mul_1(&_S3055, &_S3056, _S3050.differential_0);
    float3  _S3057 = - _S3056.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S3054.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S3056.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S3037;
    Matrix<float, 3, 3>  _S3058 = _S3053.differential_0 + _S3055.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S3058;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S3057;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3059, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3060, DiffPair_float_0 * _S3061, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3062, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3063, float _S3064)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S3059, _S3060, _S3061, _S3062, _S3063, _S3064);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S3065 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_16;
    (&dp_mean_0)->differential_0 = _S3065;
    Matrix<float, 3, 3>  _S3066 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S3066;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S3065;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S3065;
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
    float3  _S3067 = mean_17 - ray_o_3;
    *depth_10 = 0.5f * (F32_log((dot_0(_S3067, _S3067) + 9.99999997475242708e-07f)));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S3068 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    float _S3069 = 0.5f * dpdepth_0;
    DiffPair_float_0 _S3070;
    (&_S3070)->primal_0 = s_primal_ctx_dot_0(_S3068, _S3068) + 9.99999997475242708e-07f;
    (&_S3070)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3070, _S3069);
    float3  _S3071 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3072;
    (&_S3072)->primal_0 = _S3068;
    (&_S3072)->differential_0 = _S3071;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3073;
    (&_S3073)->primal_0 = _S3068;
    (&_S3073)->differential_0 = _S3071;
    s_bwd_prop_dot_0(&_S3072, &_S3073, _S3070.differential_0);
    float3  _S3074 = _S3073.differential_0 + _S3072.differential_0;
    float3  _S3075 = - _S3074;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S3071;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S3075;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S3076 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S3076;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S3074;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3077, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3078, DiffPair_float_0 * _S3079, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3080, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3081, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3082, float3  _S3083, float _S3084)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S3077, _S3078, _S3079, _S3080, _S3081, _S3082, _S3083, _S3084);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S3085 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_18;
    (&dp_mean_1)->differential_0 = _S3085;
    Matrix<float, 3, 3>  _S3086 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S3086;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S3085;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S3085;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S3085;
    s_bwd_evaluate_color_3dgs_0(&dp_mean_1, &dp_iscl_rot_1, &dp_opacity_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_5);
    *v_mean_6 = dp_mean_1.differential_0;
    *v_iscl_rot_2 = dp_iscl_rot_1.differential_0;
    *v_opacity_6 = dp_opacity_1.differential_0;
    *v_rgb_5 = dp_rgb_0.differential_0;
    *v_ray_o_2 = dp_ray_o_1.differential_0;
    *v_ray_d_2 = dp_ray_d_1.differential_0;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_14, float dOut_17)
{
    float _S3087 = _slang_select(((*dpx_14).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_14).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S3087;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_15, float dOut_18)
{
    float _S3088 = (F32_exp2(((*dpx_15).primal_0))) * 50.693145751953125f * dOut_18;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S3088;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  dOut_19)
{
    float3  _S3089 = make_float3 (1.0f) / (*dpx_16).primal_0 * dOut_19;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S3089;
    return;
}

inline __device__ float3  log_0(float3  x_41)
{
    float3  result_14;
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
        *_slang_vector_get_element_ptr(&result_14, i_10) = (F32_log((_slang_vector_get_element(x_41, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_14;
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_19, float4  quat_20, float3  scale_19, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_19, float3  t_18, float fx_24, float fy_24, float cx_19, float cy_19, FixedArray<float, 10>  * dist_coeffs_30, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_19, mean_19) + t_18;
        float _S3090 = mean_c_15.z;
        bool _S3091;
        if(_S3090 < near_plane_10)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = _S3090 > far_plane_10;
        }
        if(_S3091)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3092 = scale_19.x;
        float sx_0 = (F32_exp((_S3092)));
        float _S3093 = scale_19.y;
        float sy_0 = (F32_exp((_S3093)));
        float sz_0 = scale_19.z - 0.5f * (_S3092 + _S3093);
        float x_42 = quat_20.y;
        float x2_20 = x_42 * x_42;
        float y2_20 = quat_20.z * quat_20.z;
        float z2_35 = quat_20.w * quat_20.w;
        float xy_20 = quat_20.y * quat_20.z;
        float xz_20 = quat_20.y * quat_20.w;
        float yz_20 = quat_20.z * quat_20.w;
        float wx_20 = quat_20.x * quat_20.y;
        float wy_20 = quat_20.x * quat_20.z;
        float wz_20 = quat_20.x * quat_20.w;
        Matrix<float, 3, 3>  _S3094 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20)));
        float3  vert0_c_0 = mul_0(R_19, mul_0(_S3094, make_float3 (sx_0, 0.0f, 0.0f)) + mean_19) + t_18;
        float3  vert1_c_0 = mul_0(R_19, mul_0(_S3094, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_19) + t_18;
        float3  vert2_c_0 = mul_0(R_19, mul_0(_S3094, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_19) + t_18;
        float _S3095 = vert0_c_0.z;
        float _S3096 = vert1_c_0.z;
        float _S3097 = vert2_c_0.z;
        if(_S3095 < near_plane_10)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = _S3095 > far_plane_10;
        }
        if(_S3091)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = _S3096 < near_plane_10;
        }
        if(_S3091)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = _S3096 > far_plane_10;
        }
        if(_S3091)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = _S3097 < near_plane_10;
        }
        if(_S3091)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = _S3097 > far_plane_10;
        }
        if(_S3091)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        for(;;)
        {
            *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S3095);
            if(_S3095 < 0.0f)
            {
                _S3091 = true;
            }
            else
            {
                bool _S3098 = is_valid_distortion(*uv0_0, dist_coeffs_30);
                _S3091 = !_S3098;
            }
            if(_S3091)
            {
                break;
            }
            float u_44 = (*uv0_0).x;
            float v_44 = (*uv0_0).y;
            float r2_44 = u_44 * u_44 + v_44 * v_44;
            float2  _S3099 = *uv0_0 * make_float2 (1.0f + r2_44 * ((*dist_coeffs_30)[int(0)] + r2_44 * ((*dist_coeffs_30)[int(1)] + r2_44 * ((*dist_coeffs_30)[int(2)] + r2_44 * (*dist_coeffs_30)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_30)[int(4)] * u_44 * v_44 + (*dist_coeffs_30)[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + (*dist_coeffs_30)[int(6)] * r2_44, 2.0f * (*dist_coeffs_30)[int(5)] * u_44 * v_44 + (*dist_coeffs_30)[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + (*dist_coeffs_30)[int(7)] * r2_44);
            float2  _S3100 = _S3099 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3099.x + (*dist_coeffs_30)[int(9)] * _S3099.y, 0.0f);
            *uv0_0 = make_float2 (fx_24 * _S3100.x + cx_19, fy_24 * _S3100.y + cy_19);
            break;
        }
        bool all_valid_8 = true & (!_S3091);
        for(;;)
        {
            *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S3096);
            if(_S3096 < 0.0f)
            {
                _S3091 = true;
            }
            else
            {
                bool _S3101 = is_valid_distortion(*uv1_0, dist_coeffs_30);
                _S3091 = !_S3101;
            }
            if(_S3091)
            {
                break;
            }
            float u_45 = (*uv1_0).x;
            float v_45 = (*uv1_0).y;
            float r2_45 = u_45 * u_45 + v_45 * v_45;
            float2  _S3102 = *uv1_0 * make_float2 (1.0f + r2_45 * ((*dist_coeffs_30)[int(0)] + r2_45 * ((*dist_coeffs_30)[int(1)] + r2_45 * ((*dist_coeffs_30)[int(2)] + r2_45 * (*dist_coeffs_30)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_30)[int(4)] * u_45 * v_45 + (*dist_coeffs_30)[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + (*dist_coeffs_30)[int(6)] * r2_45, 2.0f * (*dist_coeffs_30)[int(5)] * u_45 * v_45 + (*dist_coeffs_30)[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + (*dist_coeffs_30)[int(7)] * r2_45);
            float2  _S3103 = _S3102 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3102.x + (*dist_coeffs_30)[int(9)] * _S3102.y, 0.0f);
            *uv1_0 = make_float2 (fx_24 * _S3103.x + cx_19, fy_24 * _S3103.y + cy_19);
            break;
        }
        bool all_valid_9 = all_valid_8 & (!_S3091);
        for(;;)
        {
            *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S3097);
            if(_S3097 < 0.0f)
            {
                _S3091 = true;
            }
            else
            {
                bool _S3104 = is_valid_distortion(*uv2_0, dist_coeffs_30);
                _S3091 = !_S3104;
            }
            if(_S3091)
            {
                break;
            }
            float u_46 = (*uv2_0).x;
            float v_46 = (*uv2_0).y;
            float r2_46 = u_46 * u_46 + v_46 * v_46;
            float2  _S3105 = *uv2_0 * make_float2 (1.0f + r2_46 * ((*dist_coeffs_30)[int(0)] + r2_46 * ((*dist_coeffs_30)[int(1)] + r2_46 * ((*dist_coeffs_30)[int(2)] + r2_46 * (*dist_coeffs_30)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_30)[int(4)] * u_46 * v_46 + (*dist_coeffs_30)[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + (*dist_coeffs_30)[int(6)] * r2_46, 2.0f * (*dist_coeffs_30)[int(5)] * u_46 * v_46 + (*dist_coeffs_30)[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + (*dist_coeffs_30)[int(7)] * r2_46);
            float2  _S3106 = _S3105 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3105.x + (*dist_coeffs_30)[int(9)] * _S3105.y, 0.0f);
            *uv2_0 = make_float2 (fx_24 * _S3106.x + cx_19, fy_24 * _S3106.y + cy_19);
            break;
        }
        if(!(all_valid_9 & (!_S3091)))
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_0 = *uv1_0 - *uv0_0;
        float2  e1_0 = *uv2_0 - *uv1_0;
        float offset_0 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_0.y))))) - 1.0f) * ((F32_abs((e0_0.x * e1_0.y - e0_0.y * e1_0.x))) / (length_0(e0_0) + length_0(e1_0) + length_0(*uv0_0 - *uv2_0)));
        float xmax_5 = (F32_max(((F32_max(((*uv0_0).x), ((*uv1_0).x)))), ((*uv2_0).x))) + offset_0;
        float xmin_5 = (F32_min(((F32_min(((*uv0_0).x), ((*uv1_0).x)))), ((*uv2_0).x))) - offset_0;
        float ymax_5 = (F32_max(((F32_max(((*uv0_0).y), ((*uv1_0).y)))), ((*uv2_0).y))) + offset_0;
        float ymin_5 = (F32_min(((F32_min(((*uv0_0).y), ((*uv1_0).y)))), ((*uv2_0).y))) - offset_0;
        if(xmax_5 <= 0.0f)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = xmin_5 >= float(image_width_15);
        }
        if(_S3091)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = ymax_5 <= 0.0f;
        }
        if(_S3091)
        {
            _S3091 = true;
        }
        else
        {
            _S3091 = ymin_5 >= float(image_height_15);
        }
        if(_S3091)
        {
            _S3091 = true;
        }
        else
        {
            if(_S3090 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S3091 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S3091 = false;
                }
                if(_S3091)
                {
                    _S3091 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S3091 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S3091 = false;
                    }
                }
            }
            else
            {
                _S3091 = false;
            }
        }
        if(_S3091)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S3107 = mean_19 - - mul_0(transpose_0(R_19), t_18);
        float _S3108 = _S3107.x;
        float _S3109 = _S3107.y;
        float _S3110 = _S3107.z;
        float norm_10 = (F32_sqrt((_S3108 * _S3108 + _S3109 * _S3109 + _S3110 * _S3110)));
        float x_43 = _S3108 / norm_10;
        float y_18 = _S3109 / norm_10;
        float z_15 = _S3110 / norm_10;
        float z2_36 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_43 * x_43 - y_18 * y_18;
        float fS1_15 = 2.0f * x_43 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_36 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_43) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_36 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_43) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_43 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_36 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_43) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_43 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S3111 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S3111);
        float3  _S3112 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S3113 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S3112 + _S3113 + make_float3 (0.5f), _S3111);
        (*rgb_12)[int(2)] = max_0(_S3112 - _S3113 + make_float3 (0.5f), _S3111);
        float3  _S3114 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S3114 * make_float3 (float(- (F32_sign((dot_0(_S3114, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_20, float3  t_19, float fx_25, float fy_25, float cx_20, float cy_20, FixedArray<float, 10>  * dist_coeffs_31, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    bool _S3115;
    bool _S3116;
    bool _S3117;
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_20, mean_20) + t_19;
        float _S3118 = length_1(mean_c_16);
        bool _S3119;
        if(_S3118 < near_plane_11)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = _S3118 > far_plane_11;
        }
        if(_S3119)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3120 = scale_20.x;
        float sx_1 = (F32_exp((_S3120)));
        float _S3121 = scale_20.y;
        float sy_1 = (F32_exp((_S3121)));
        float sz_1 = scale_20.z - 0.5f * (_S3120 + _S3121);
        float x_44 = quat_21.y;
        float x2_21 = x_44 * x_44;
        float y2_21 = quat_21.z * quat_21.z;
        float z2_37 = quat_21.w * quat_21.w;
        float xy_21 = quat_21.y * quat_21.z;
        float xz_21 = quat_21.y * quat_21.w;
        float yz_21 = quat_21.z * quat_21.w;
        float wx_21 = quat_21.x * quat_21.y;
        float wy_21 = quat_21.x * quat_21.z;
        float wz_21 = quat_21.x * quat_21.w;
        Matrix<float, 3, 3>  _S3122 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_37), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_37), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_1 = mul_0(R_20, mul_0(_S3122, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_19;
        float3  vert1_c_1 = mul_0(R_20, mul_0(_S3122, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_19;
        float3  vert2_c_1 = mul_0(R_20, mul_0(_S3122, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_19;
        float _S3123 = length_1(vert0_c_1);
        float _S3124 = length_1(vert1_c_1);
        float _S3125 = length_1(vert2_c_1);
        if(_S3123 < near_plane_11)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = _S3123 > far_plane_11;
        }
        if(_S3119)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = _S3124 < near_plane_11;
        }
        if(_S3119)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = _S3124 > far_plane_11;
        }
        if(_S3119)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = _S3125 < near_plane_11;
        }
        if(_S3119)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = _S3125 > far_plane_11;
        }
        if(_S3119)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float k_7;
        for(;;)
        {
            float2  _S3126 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_24 = length_0(_S3126);
            float _S3127 = vert0_c_1.z;
            float theta_20 = (F32_atan2((r_24), (_S3127)));
            if(theta_20 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_20 * theta_20 / 3.0f) / _S3127;
            }
            else
            {
                k_7 = theta_20 / r_24;
            }
            float2  _S3128 = _S3126 * make_float2 (k_7);
            *uv0_1 = _S3128;
            bool _S3129 = is_valid_distortion(_S3128, dist_coeffs_31);
            bool _S3130 = !_S3129;
            _S3115 = _S3130;
            if(_S3130)
            {
                break;
            }
            float u_47 = (*uv0_1).x;
            float v_47 = (*uv0_1).y;
            float r2_47 = u_47 * u_47 + v_47 * v_47;
            float2  _S3131 = *uv0_1 * make_float2 (1.0f + r2_47 * ((*dist_coeffs_31)[int(0)] + r2_47 * ((*dist_coeffs_31)[int(1)] + r2_47 * ((*dist_coeffs_31)[int(2)] + r2_47 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_47 * v_47 + (*dist_coeffs_31)[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + (*dist_coeffs_31)[int(6)] * r2_47, 2.0f * (*dist_coeffs_31)[int(5)] * u_47 * v_47 + (*dist_coeffs_31)[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + (*dist_coeffs_31)[int(7)] * r2_47);
            float2  _S3132 = _S3131 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3131.x + (*dist_coeffs_31)[int(9)] * _S3131.y, 0.0f);
            *uv0_1 = make_float2 (fx_25 * _S3132.x + cx_20, fy_25 * _S3132.y + cy_20);
            break;
        }
        bool all_valid_10 = true & (!_S3115);
        for(;;)
        {
            float2  _S3133 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_25 = length_0(_S3133);
            float _S3134 = vert1_c_1.z;
            float theta_21 = (F32_atan2((r_25), (_S3134)));
            if(theta_21 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_21 * theta_21 / 3.0f) / _S3134;
            }
            else
            {
                k_7 = theta_21 / r_25;
            }
            float2  _S3135 = _S3133 * make_float2 (k_7);
            *uv1_1 = _S3135;
            bool _S3136 = is_valid_distortion(_S3135, dist_coeffs_31);
            bool _S3137 = !_S3136;
            _S3116 = _S3137;
            if(_S3137)
            {
                break;
            }
            float u_48 = (*uv1_1).x;
            float v_48 = (*uv1_1).y;
            float r2_48 = u_48 * u_48 + v_48 * v_48;
            float2  _S3138 = *uv1_1 * make_float2 (1.0f + r2_48 * ((*dist_coeffs_31)[int(0)] + r2_48 * ((*dist_coeffs_31)[int(1)] + r2_48 * ((*dist_coeffs_31)[int(2)] + r2_48 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_48 * v_48 + (*dist_coeffs_31)[int(5)] * (r2_48 + 2.0f * u_48 * u_48) + (*dist_coeffs_31)[int(6)] * r2_48, 2.0f * (*dist_coeffs_31)[int(5)] * u_48 * v_48 + (*dist_coeffs_31)[int(4)] * (r2_48 + 2.0f * v_48 * v_48) + (*dist_coeffs_31)[int(7)] * r2_48);
            float2  _S3139 = _S3138 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3138.x + (*dist_coeffs_31)[int(9)] * _S3138.y, 0.0f);
            *uv1_1 = make_float2 (fx_25 * _S3139.x + cx_20, fy_25 * _S3139.y + cy_20);
            break;
        }
        bool all_valid_11 = all_valid_10 & (!_S3116);
        for(;;)
        {
            float2  _S3140 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_26 = length_0(_S3140);
            float _S3141 = vert2_c_1.z;
            float theta_22 = (F32_atan2((r_26), (_S3141)));
            if(theta_22 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_22 * theta_22 / 3.0f) / _S3141;
            }
            else
            {
                k_7 = theta_22 / r_26;
            }
            float2  _S3142 = _S3140 * make_float2 (k_7);
            *uv2_1 = _S3142;
            bool _S3143 = is_valid_distortion(_S3142, dist_coeffs_31);
            bool _S3144 = !_S3143;
            _S3117 = _S3144;
            if(_S3144)
            {
                break;
            }
            float u_49 = (*uv2_1).x;
            float v_49 = (*uv2_1).y;
            float r2_49 = u_49 * u_49 + v_49 * v_49;
            float2  _S3145 = *uv2_1 * make_float2 (1.0f + r2_49 * ((*dist_coeffs_31)[int(0)] + r2_49 * ((*dist_coeffs_31)[int(1)] + r2_49 * ((*dist_coeffs_31)[int(2)] + r2_49 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_49 * v_49 + (*dist_coeffs_31)[int(5)] * (r2_49 + 2.0f * u_49 * u_49) + (*dist_coeffs_31)[int(6)] * r2_49, 2.0f * (*dist_coeffs_31)[int(5)] * u_49 * v_49 + (*dist_coeffs_31)[int(4)] * (r2_49 + 2.0f * v_49 * v_49) + (*dist_coeffs_31)[int(7)] * r2_49);
            float2  _S3146 = _S3145 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3145.x + (*dist_coeffs_31)[int(9)] * _S3145.y, 0.0f);
            *uv2_1 = make_float2 (fx_25 * _S3146.x + cx_20, fy_25 * _S3146.y + cy_20);
            break;
        }
        if(!(all_valid_11 & (!_S3117)))
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_1 = *uv1_1 - *uv0_1;
        float2  e1_1 = *uv2_1 - *uv1_1;
        float offset_1 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_1.y))))) - 1.0f) * ((F32_abs((e0_1.x * e1_1.y - e0_1.y * e1_1.x))) / (length_0(e0_1) + length_0(e1_1) + length_0(*uv0_1 - *uv2_1)));
        float xmax_6 = (F32_max(((F32_max(((*uv0_1).x), ((*uv1_1).x)))), ((*uv2_1).x))) + offset_1;
        float xmin_6 = (F32_min(((F32_min(((*uv0_1).x), ((*uv1_1).x)))), ((*uv2_1).x))) - offset_1;
        float ymax_6 = (F32_max(((F32_max(((*uv0_1).y), ((*uv1_1).y)))), ((*uv2_1).y))) + offset_1;
        float ymin_6 = (F32_min(((F32_min(((*uv0_1).y), ((*uv1_1).y)))), ((*uv2_1).y))) - offset_1;
        if(xmax_6 <= 0.0f)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = xmin_6 >= float(image_width_16);
        }
        if(_S3119)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = ymax_6 <= 0.0f;
        }
        if(_S3119)
        {
            _S3119 = true;
        }
        else
        {
            _S3119 = ymin_6 >= float(image_height_16);
        }
        if(_S3119)
        {
            _S3119 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S3119 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S3119 = false;
                }
                if(_S3119)
                {
                    _S3119 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S3119 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S3119 = false;
                    }
                }
            }
            else
            {
                _S3119 = false;
            }
        }
        if(_S3119)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S3123, _S3124, _S3125) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S3147 = mean_20 - - mul_0(transpose_0(R_20), t_19);
        float _S3148 = _S3147.x;
        float _S3149 = _S3147.y;
        float _S3150 = _S3147.z;
        float norm_11 = (F32_sqrt((_S3148 * _S3148 + _S3149 * _S3149 + _S3150 * _S3150)));
        float x_45 = _S3148 / norm_11;
        float y_19 = _S3149 / norm_11;
        float z_16 = _S3150 / norm_11;
        float z2_38 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_45 * x_45 - y_19 * y_19;
        float fS1_16 = 2.0f * x_45 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_38 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_45) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_38 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_45) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_45 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_38 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_45) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_45 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S3151 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S3151);
        float3  _S3152 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S3153 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S3152 + _S3153 + make_float3 (0.5f), _S3151);
        (*rgb_13)[int(2)] = max_0(_S3152 - _S3153 + make_float3 (0.5f), _S3151);
        float3  _S3154 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S3154 * make_float3 (float(- (F32_sign((dot_0(_S3154, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_21, float3  t_20, float fx_26, float fy_26, float cx_21, float cy_21, FixedArray<float, 10>  * dist_coeffs_32, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_21, mean_21) + t_20;
    float _S3155 = scale_21.x;
    float sx_2 = (F32_exp((_S3155)));
    float _S3156 = scale_21.y;
    float sy_2 = (F32_exp((_S3156)));
    float sz_2 = scale_21.z - 0.5f * (_S3155 + _S3156);
    float x_46 = quat_22.y;
    float x2_22 = x_46 * x_46;
    float y2_22 = quat_22.z * quat_22.z;
    float z2_39 = quat_22.w * quat_22.w;
    float xy_22 = quat_22.y * quat_22.z;
    float xz_22 = quat_22.y * quat_22.w;
    float yz_22 = quat_22.z * quat_22.w;
    float wx_22 = quat_22.x * quat_22.y;
    float wy_22 = quat_22.x * quat_22.z;
    float wz_22 = quat_22.x * quat_22.w;
    Matrix<float, 3, 3>  _S3157 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_39), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_39), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
    float3  vert0_c_2 = mul_0(R_21, mul_0(_S3157, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_20;
    float3  vert1_c_2 = mul_0(R_21, mul_0(_S3157, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_20;
    float3  vert2_c_2 = mul_0(R_21, mul_0(_S3157, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_20;
    float2  _S3158 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_50 = _S3158.x;
    float v_50 = _S3158.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S3159 = 2.0f * (*dist_coeffs_32)[int(4)];
    float _S3160 = 2.0f * (*dist_coeffs_32)[int(5)];
    float2  _S3161 = _S3158 * make_float2 (1.0f + r2_50 * ((*dist_coeffs_32)[int(0)] + r2_50 * ((*dist_coeffs_32)[int(1)] + r2_50 * ((*dist_coeffs_32)[int(2)] + r2_50 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S3159 * u_50 * v_50 + (*dist_coeffs_32)[int(5)] * (r2_50 + 2.0f * u_50 * u_50) + (*dist_coeffs_32)[int(6)] * r2_50, _S3160 * u_50 * v_50 + (*dist_coeffs_32)[int(4)] * (r2_50 + 2.0f * v_50 * v_50) + (*dist_coeffs_32)[int(7)] * r2_50);
    float2  _S3162 = _S3161 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3161.x + (*dist_coeffs_32)[int(9)] * _S3161.y, 0.0f);
    *uv0_2 = make_float2 (fx_26 * _S3162.x + cx_21, fy_26 * _S3162.y + cy_21);
    float2  _S3163 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_51 = _S3163.x;
    float v_51 = _S3163.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float2  _S3164 = _S3163 * make_float2 (1.0f + r2_51 * ((*dist_coeffs_32)[int(0)] + r2_51 * ((*dist_coeffs_32)[int(1)] + r2_51 * ((*dist_coeffs_32)[int(2)] + r2_51 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S3159 * u_51 * v_51 + (*dist_coeffs_32)[int(5)] * (r2_51 + 2.0f * u_51 * u_51) + (*dist_coeffs_32)[int(6)] * r2_51, _S3160 * u_51 * v_51 + (*dist_coeffs_32)[int(4)] * (r2_51 + 2.0f * v_51 * v_51) + (*dist_coeffs_32)[int(7)] * r2_51);
    float2  _S3165 = _S3164 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3164.x + (*dist_coeffs_32)[int(9)] * _S3164.y, 0.0f);
    *uv1_2 = make_float2 (fx_26 * _S3165.x + cx_21, fy_26 * _S3165.y + cy_21);
    float2  _S3166 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_52 = _S3166.x;
    float v_52 = _S3166.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float2  _S3167 = _S3166 * make_float2 (1.0f + r2_52 * ((*dist_coeffs_32)[int(0)] + r2_52 * ((*dist_coeffs_32)[int(1)] + r2_52 * ((*dist_coeffs_32)[int(2)] + r2_52 * (*dist_coeffs_32)[int(3)])))) + make_float2 (_S3159 * u_52 * v_52 + (*dist_coeffs_32)[int(5)] * (r2_52 + 2.0f * u_52 * u_52) + (*dist_coeffs_32)[int(6)] * r2_52, _S3160 * u_52 * v_52 + (*dist_coeffs_32)[int(4)] * (r2_52 + 2.0f * v_52 * v_52) + (*dist_coeffs_32)[int(7)] * r2_52);
    float2  _S3168 = _S3167 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3167.x + (*dist_coeffs_32)[int(9)] * _S3167.y, 0.0f);
    float _S3169 = fx_26 * _S3168.x + cx_21;
    float _S3170 = fy_26 * _S3168.y + cy_21;
    float2  _S3171 = make_float2 (_S3169, _S3170);
    *uv2_2 = _S3171;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S3171 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S3171)));
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S3169))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S3170))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S3169))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S3170))) + offset_2)))));
    *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_2 = hardness_2;
    float3  _S3172 = mean_21 - - mul_0(transpose_0(R_21), t_20);
    float _S3173 = _S3172.x;
    float _S3174 = _S3172.y;
    float _S3175 = _S3172.z;
    float norm_12 = (F32_sqrt((_S3173 * _S3173 + _S3174 * _S3174 + _S3175 * _S3175)));
    float x_47 = _S3173 / norm_12;
    float y_20 = _S3174 / norm_12;
    float z_17 = _S3175 / norm_12;
    float z2_40 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_47 * x_47 - y_20 * y_20;
    float fS1_17 = 2.0f * x_47 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_40 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_47) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_40 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_47) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_47 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_40 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_47) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_47 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S3176 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S3176);
    float3  _S3177 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S3178 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S3177 + _S3178 + make_float3 (0.5f), _S3176);
    (*rgb_14)[int(2)] = max_0(_S3177 - _S3178 + make_float3 (0.5f), _S3176);
    float3  _S3179 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S3179 * make_float3 (float(- (F32_sign((dot_0(_S3179, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_22, float3  t_21, float fx_27, float fy_27, float cx_22, float cy_22, FixedArray<float, 10>  * dist_coeffs_33, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_22, mean_22) + t_21;
    float _S3180 = scale_22.x;
    float sx_3 = (F32_exp((_S3180)));
    float _S3181 = scale_22.y;
    float sy_3 = (F32_exp((_S3181)));
    float sz_3 = scale_22.z - 0.5f * (_S3180 + _S3181);
    float x_48 = quat_23.y;
    float x2_23 = x_48 * x_48;
    float y2_23 = quat_23.z * quat_23.z;
    float z2_41 = quat_23.w * quat_23.w;
    float xy_23 = quat_23.y * quat_23.z;
    float xz_23 = quat_23.y * quat_23.w;
    float yz_23 = quat_23.z * quat_23.w;
    float wx_23 = quat_23.x * quat_23.y;
    float wy_23 = quat_23.x * quat_23.z;
    float wz_23 = quat_23.x * quat_23.w;
    Matrix<float, 3, 3>  _S3182 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_41), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_41), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_3 = mul_0(R_22, mul_0(_S3182, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_21;
    float3  vert1_c_3 = mul_0(R_22, mul_0(_S3182, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_21;
    float3  vert2_c_3 = mul_0(R_22, mul_0(_S3182, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_21;
    float2  _S3183 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_27 = length_0(_S3183);
    float _S3184 = vert0_c_3.z;
    float theta_23 = (F32_atan2((r_27), (_S3184)));
    float k_8;
    if(theta_23 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_23 * theta_23 / 3.0f) / _S3184;
    }
    else
    {
        k_8 = theta_23 / r_27;
    }
    float2  _S3185 = _S3183 * make_float2 (k_8);
    float u_53 = _S3185.x;
    float v_53 = _S3185.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S3186 = 2.0f * (*dist_coeffs_33)[int(4)];
    float _S3187 = 2.0f * (*dist_coeffs_33)[int(5)];
    float2  _S3188 = _S3185 * make_float2 (1.0f + r2_53 * ((*dist_coeffs_33)[int(0)] + r2_53 * ((*dist_coeffs_33)[int(1)] + r2_53 * ((*dist_coeffs_33)[int(2)] + r2_53 * (*dist_coeffs_33)[int(3)])))) + make_float2 (_S3186 * u_53 * v_53 + (*dist_coeffs_33)[int(5)] * (r2_53 + 2.0f * u_53 * u_53) + (*dist_coeffs_33)[int(6)] * r2_53, _S3187 * u_53 * v_53 + (*dist_coeffs_33)[int(4)] * (r2_53 + 2.0f * v_53 * v_53) + (*dist_coeffs_33)[int(7)] * r2_53);
    float2  _S3189 = _S3188 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3188.x + (*dist_coeffs_33)[int(9)] * _S3188.y, 0.0f);
    *uv0_3 = make_float2 (fx_27 * _S3189.x + cx_22, fy_27 * _S3189.y + cy_22);
    float2  _S3190 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_28 = length_0(_S3190);
    float _S3191 = vert1_c_3.z;
    float theta_24 = (F32_atan2((r_28), (_S3191)));
    if(theta_24 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_24 * theta_24 / 3.0f) / _S3191;
    }
    else
    {
        k_8 = theta_24 / r_28;
    }
    float2  _S3192 = _S3190 * make_float2 (k_8);
    float u_54 = _S3192.x;
    float v_54 = _S3192.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float2  _S3193 = _S3192 * make_float2 (1.0f + r2_54 * ((*dist_coeffs_33)[int(0)] + r2_54 * ((*dist_coeffs_33)[int(1)] + r2_54 * ((*dist_coeffs_33)[int(2)] + r2_54 * (*dist_coeffs_33)[int(3)])))) + make_float2 (_S3186 * u_54 * v_54 + (*dist_coeffs_33)[int(5)] * (r2_54 + 2.0f * u_54 * u_54) + (*dist_coeffs_33)[int(6)] * r2_54, _S3187 * u_54 * v_54 + (*dist_coeffs_33)[int(4)] * (r2_54 + 2.0f * v_54 * v_54) + (*dist_coeffs_33)[int(7)] * r2_54);
    float2  _S3194 = _S3193 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3193.x + (*dist_coeffs_33)[int(9)] * _S3193.y, 0.0f);
    *uv1_3 = make_float2 (fx_27 * _S3194.x + cx_22, fy_27 * _S3194.y + cy_22);
    float2  _S3195 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_29 = length_0(_S3195);
    float _S3196 = vert2_c_3.z;
    float theta_25 = (F32_atan2((r_29), (_S3196)));
    if(theta_25 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_25 * theta_25 / 3.0f) / _S3196;
    }
    else
    {
        k_8 = theta_25 / r_29;
    }
    float2  _S3197 = _S3195 * make_float2 (k_8);
    float u_55 = _S3197.x;
    float v_55 = _S3197.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float2  _S3198 = _S3197 * make_float2 (1.0f + r2_55 * ((*dist_coeffs_33)[int(0)] + r2_55 * ((*dist_coeffs_33)[int(1)] + r2_55 * ((*dist_coeffs_33)[int(2)] + r2_55 * (*dist_coeffs_33)[int(3)])))) + make_float2 (_S3186 * u_55 * v_55 + (*dist_coeffs_33)[int(5)] * (r2_55 + 2.0f * u_55 * u_55) + (*dist_coeffs_33)[int(6)] * r2_55, _S3187 * u_55 * v_55 + (*dist_coeffs_33)[int(4)] * (r2_55 + 2.0f * v_55 * v_55) + (*dist_coeffs_33)[int(7)] * r2_55);
    float2  _S3199 = _S3198 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3198.x + (*dist_coeffs_33)[int(9)] * _S3198.y, 0.0f);
    float _S3200 = fx_27 * _S3199.x + cx_22;
    float _S3201 = fy_27 * _S3199.y + cy_22;
    float2  _S3202 = make_float2 (_S3200, _S3201);
    *uv2_3 = _S3202;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S3202 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S3202)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S3200))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S3201))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S3200))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S3201))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S3203 = mean_22 - - mul_0(transpose_0(R_22), t_21);
    float _S3204 = _S3203.x;
    float _S3205 = _S3203.y;
    float _S3206 = _S3203.z;
    float norm_13 = (F32_sqrt((_S3204 * _S3204 + _S3205 * _S3205 + _S3206 * _S3206)));
    float x_49 = _S3204 / norm_13;
    float y_21 = _S3205 / norm_13;
    float z_18 = _S3206 / norm_13;
    float z2_42 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_49 * x_49 - y_21 * y_21;
    float fS1_18 = 2.0f * x_49 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_42 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_49) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_42 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_49) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_49 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_42 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_49) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_49 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S3207 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S3207);
    float3  _S3208 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S3209 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S3208 + _S3209 + make_float3 (0.5f), _S3207);
    (*rgb_15)[int(2)] = max_0(_S3208 - _S3209 + make_float3 (0.5f), _S3207);
    float3  _S3210 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S3210 * make_float3 (float(- (F32_sign((dot_0(_S3210, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3211, float3  _S3212)
{
    _d_log_vector_0(_S3211, _S3212);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S3213, float _S3214)
{
    _d_exp2_0(_S3213, _S3214);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S3215, float _S3216)
{
    _d_abs_0(_S3215, _S3216);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_23, float3  t_22, float fx_28, float fy_28, float cx_23, float cy_23, FixedArray<float, 10>  * dist_coeffs_34, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_1(R_23, mean_23) + t_22;
    float _S3217 = scale_23.x;
    float _S3218 = s_primal_ctx_exp_1(_S3217);
    float _S3219 = scale_23.y;
    float _S3220 = s_primal_ctx_exp_1(_S3219);
    float sz_4 = scale_23.z - 0.5f * (_S3217 + _S3219);
    float _S3221 = quat_24.y;
    float x2_24 = _S3221 * _S3221;
    float y2_24 = quat_24.z * quat_24.z;
    float z2_43 = quat_24.w * quat_24.w;
    float xy_24 = quat_24.y * quat_24.z;
    float xz_24 = quat_24.y * quat_24.w;
    float yz_24 = quat_24.z * quat_24.w;
    float wx_24 = quat_24.x * quat_24.y;
    float wy_24 = quat_24.x * quat_24.z;
    float wz_24 = quat_24.x * quat_24.w;
    Matrix<float, 3, 3>  _S3222 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_43), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_43), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  _S3223 = make_float3 (_S3218, 0.0f, 0.0f);
    float3  vert0_0 = s_primal_ctx_mul_1(_S3222, _S3223) + mean_23;
    float _S3224 = -0.5f + sz_4;
    float3  _S3225 = make_float3 (_S3218 * _S3224, _S3220, 0.0f);
    float3  vert1_0 = s_primal_ctx_mul_1(_S3222, _S3225) + mean_23;
    float _S3226 = -0.5f - sz_4;
    float3  _S3227 = make_float3 (_S3218 * _S3226, - _S3220, 0.0f);
    float3  vert2_0 = s_primal_ctx_mul_1(_S3222, _S3227) + mean_23;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_23, vert0_0) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_23, vert1_0) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_23, vert2_0) + t_22;
    float2  _S3228 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S3229 = vert0_c_4.z;
    float2  _S3230 = make_float2 (_S3229);
    float2  _S3231 = _S3228 / make_float2 (_S3229);
    float2  _S3232 = make_float2 (_S3229 * _S3229);
    float u_56 = _S3231.x;
    float v_56 = _S3231.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S3233 = (*dist_coeffs_34)[int(2)] + r2_56 * (*dist_coeffs_34)[int(3)];
    float _S3234 = (*dist_coeffs_34)[int(1)] + r2_56 * _S3233;
    float _S3235 = (*dist_coeffs_34)[int(0)] + r2_56 * _S3234;
    float radial_2 = 1.0f + r2_56 * _S3235;
    float2  _S3236 = make_float2 (radial_2);
    float _S3237 = 2.0f * (*dist_coeffs_34)[int(4)];
    float _S3238 = _S3237 * u_56;
    float _S3239 = 2.0f * u_56;
    float _S3240 = 2.0f * (*dist_coeffs_34)[int(5)];
    float _S3241 = _S3240 * u_56;
    float _S3242 = 2.0f * v_56;
    float2  _S3243 = _S3231 * make_float2 (radial_2) + make_float2 (_S3238 * v_56 + (*dist_coeffs_34)[int(5)] * (r2_56 + _S3239 * u_56) + (*dist_coeffs_34)[int(6)] * r2_56, _S3241 * v_56 + (*dist_coeffs_34)[int(4)] * (r2_56 + _S3242 * v_56) + (*dist_coeffs_34)[int(7)] * r2_56);
    float2  _S3244 = _S3243 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3243.x + (*dist_coeffs_34)[int(9)] * _S3243.y, 0.0f);
    float _S3245 = fx_28 * _S3244.x + cx_23;
    float _S3246 = fy_28 * _S3244.y + cy_23;
    float2  _S3247 = make_float2 (_S3245, _S3246);
    float2  _S3248 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S3249 = vert1_c_4.z;
    float2  _S3250 = make_float2 (_S3249);
    float2  _S3251 = _S3248 / make_float2 (_S3249);
    float2  _S3252 = make_float2 (_S3249 * _S3249);
    float u_57 = _S3251.x;
    float v_57 = _S3251.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S3253 = (*dist_coeffs_34)[int(2)] + r2_57 * (*dist_coeffs_34)[int(3)];
    float _S3254 = (*dist_coeffs_34)[int(1)] + r2_57 * _S3253;
    float _S3255 = (*dist_coeffs_34)[int(0)] + r2_57 * _S3254;
    float radial_3 = 1.0f + r2_57 * _S3255;
    float2  _S3256 = make_float2 (radial_3);
    float _S3257 = _S3237 * u_57;
    float _S3258 = 2.0f * u_57;
    float _S3259 = _S3240 * u_57;
    float _S3260 = 2.0f * v_57;
    float2  _S3261 = _S3251 * make_float2 (radial_3) + make_float2 (_S3257 * v_57 + (*dist_coeffs_34)[int(5)] * (r2_57 + _S3258 * u_57) + (*dist_coeffs_34)[int(6)] * r2_57, _S3259 * v_57 + (*dist_coeffs_34)[int(4)] * (r2_57 + _S3260 * v_57) + (*dist_coeffs_34)[int(7)] * r2_57);
    float2  _S3262 = _S3261 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3261.x + (*dist_coeffs_34)[int(9)] * _S3261.y, 0.0f);
    float _S3263 = fx_28 * _S3262.x + cx_23;
    float _S3264 = fy_28 * _S3262.y + cy_23;
    float2  _S3265 = make_float2 (_S3263, _S3264);
    float2  _S3266 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S3267 = vert2_c_4.z;
    float2  _S3268 = make_float2 (_S3267);
    float2  _S3269 = _S3266 / make_float2 (_S3267);
    float2  _S3270 = make_float2 (_S3267 * _S3267);
    float u_58 = _S3269.x;
    float v_58 = _S3269.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S3271 = (*dist_coeffs_34)[int(2)] + r2_58 * (*dist_coeffs_34)[int(3)];
    float _S3272 = (*dist_coeffs_34)[int(1)] + r2_58 * _S3271;
    float _S3273 = (*dist_coeffs_34)[int(0)] + r2_58 * _S3272;
    float radial_4 = 1.0f + r2_58 * _S3273;
    float2  _S3274 = make_float2 (radial_4);
    float _S3275 = _S3237 * u_58;
    float _S3276 = 2.0f * u_58;
    float _S3277 = _S3240 * u_58;
    float _S3278 = 2.0f * v_58;
    float2  _S3279 = _S3269 * make_float2 (radial_4) + make_float2 (_S3275 * v_58 + (*dist_coeffs_34)[int(5)] * (r2_58 + _S3276 * u_58) + (*dist_coeffs_34)[int(6)] * r2_58, _S3277 * v_58 + (*dist_coeffs_34)[int(4)] * (r2_58 + _S3278 * v_58) + (*dist_coeffs_34)[int(7)] * r2_58);
    float2  _S3280 = _S3279 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3279.x + (*dist_coeffs_34)[int(9)] * _S3279.y, 0.0f);
    float _S3281 = fx_28 * _S3280.x + cx_23;
    float _S3282 = fy_28 * _S3280.y + cy_23;
    float2  _S3283 = make_float2 (_S3281, _S3282);
    float2  e0_4 = _S3265 - _S3247;
    float2  e1_4 = _S3283 - _S3265;
    float2  e2_0 = _S3247 - _S3283;
    float _S3284 = e0_4.x;
    float _S3285 = e1_4.y;
    float _S3286 = e0_4.y;
    float _S3287 = e1_4.x;
    float _S3288 = _S3284 * _S3285 - _S3286 * _S3287;
    float _S3289 = 1.0f - hardness_4.y;
    float _S3290 = -1.0f / _S3289;
    float _S3291 = _S3289 * _S3289;
    float _S3292 = s_primal_ctx_max_0(_S3245, _S3263);
    float _S3293 = s_primal_ctx_min_0(_S3245, _S3263);
    float _S3294 = s_primal_ctx_max_0(_S3246, _S3264);
    float _S3295 = s_primal_ctx_min_0(_S3246, _S3264);
    float3  _S3296 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3297 = transpose_0(R_23);
    float3  _S3298 = mean_23 - - s_primal_ctx_mul_1(_S3297, t_22);
    float _S3299 = _S3298.x;
    float _S3300 = _S3298.y;
    float _S3301 = _S3298.z;
    float _S3302 = _S3299 * _S3299 + _S3300 * _S3300 + _S3301 * _S3301;
    float _S3303 = s_primal_ctx_sqrt_0(_S3302);
    float x_50 = _S3299 / _S3303;
    float3  _S3304 = make_float3 (x_50);
    float _S3305 = _S3303 * _S3303;
    float y_22 = _S3300 / _S3303;
    float z_19 = _S3301 / _S3303;
    float3  _S3306 = make_float3 (z_19);
    float _S3307 = - y_22;
    float3  _S3308 = make_float3 (_S3307);
    float z2_44 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_50 * x_50 - y_22 * y_22;
    float _S3309 = 2.0f * x_50;
    float fS1_19 = _S3309 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_44 - 0.31539157032966614f;
    float3  _S3310 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_50;
    float3  _S3311 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S3312 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S3313 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S3314 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_44 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S3315 = 1.86588168144226074f * z2_44 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S3315;
    float3  _S3316 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_50;
    float3  _S3317 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S3318 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S3319 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S3320 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_50 * fC1_19 - y_22 * fS1_19);
    float3  _S3321 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_50 * fS1_19 + y_22 * fC1_19);
    float3  _S3322 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3307) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_50) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S3323 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S3324 = make_float3 (0.0f);
    float3  _S3325 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S3326 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3327 = make_float3 (_S3326);
    float3  _S3328 = (*ch_coeffs_4)[int(1)] * make_float3 (_S3326);
    float3  _S3329 = _S3325 + _S3328 + make_float3 (0.5f);
    float3  _S3330 = _S3325 - _S3328 + make_float3 (0.5f);
    float3  _S3331 = vert1_c_4 - vert0_c_4;
    float3  _S3332 = vert2_c_4 - vert0_c_4;
    float3  _S3333 = s_primal_ctx_cross_0(_S3331, _S3332);
    float3  _S3334 = normalize_0(_S3333);
    float3  _S3335 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3334, mean_c_19)))))) * v_normal_0;
    float3  _S3336 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3337;
    (&_S3337)->primal_0 = _S3334;
    (&_S3337)->differential_0 = _S3336;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3338;
    (&_S3338)->primal_0 = mean_c_19;
    (&_S3338)->differential_0 = _S3336;
    s_bwd_prop_dot_0(&_S3337, &_S3338, 0.0f);
    float3  _S3339 = _S3335 + _S3337.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3340;
    (&_S3340)->primal_0 = _S3333;
    (&_S3340)->differential_0 = _S3336;
    s_bwd_normalize_impl_0(&_S3340, _S3339);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3341;
    (&_S3341)->primal_0 = _S3331;
    (&_S3341)->differential_0 = _S3336;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3342;
    (&_S3342)->primal_0 = _S3332;
    (&_S3342)->differential_0 = _S3336;
    s_bwd_prop_cross_0(&_S3341, &_S3342, _S3340.differential_0);
    float3  _S3343 = - _S3342.differential_0;
    float3  _S3344 = - _S3341.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3345;
    (&_S3345)->primal_0 = _S3330;
    (&_S3345)->differential_0 = _S3336;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3346;
    (&_S3346)->primal_0 = _S3324;
    (&_S3346)->differential_0 = _S3336;
    s_bwd_prop_max_0(&_S3345, &_S3346, (*v_rgb_6)[int(2)]);
    float3  _S3347 = - _S3345.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3348;
    (&_S3348)->primal_0 = _S3329;
    (&_S3348)->differential_0 = _S3336;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3349;
    (&_S3349)->primal_0 = _S3324;
    (&_S3349)->differential_0 = _S3336;
    s_bwd_prop_max_0(&_S3348, &_S3349, (*v_rgb_6)[int(1)]);
    float3  _S3350 = _S3327 * (_S3347 + _S3348.differential_0);
    float3  _S3351 = _S3345.differential_0 + _S3348.differential_0;
    float3  _S3352 = make_float3 (0.5f) * - _S3351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3353;
    (&_S3353)->primal_0 = _S3323;
    (&_S3353)->differential_0 = _S3336;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3354;
    (&_S3354)->primal_0 = _S3324;
    (&_S3354)->differential_0 = _S3336;
    s_bwd_prop_max_0(&_S3353, &_S3354, (*v_rgb_6)[int(0)]);
    float3  _S3355 = _S3352 + _S3353.differential_0;
    float3  _S3356 = _S3351 + _S3353.differential_0;
    float3  _S3357 = _S3321 * _S3356;
    float3  _S3358 = (*sh_coeffs_19)[int(15)] * _S3356;
    float3  _S3359 = _S3319 * _S3356;
    float3  _S3360 = (*sh_coeffs_19)[int(14)] * _S3356;
    float3  _S3361 = _S3317 * _S3356;
    float3  _S3362 = (*sh_coeffs_19)[int(13)] * _S3356;
    float3  _S3363 = _S3316 * _S3356;
    float3  _S3364 = (*sh_coeffs_19)[int(12)] * _S3356;
    float3  _S3365 = _S3318 * _S3356;
    float3  _S3366 = (*sh_coeffs_19)[int(11)] * _S3356;
    float3  _S3367 = _S3320 * _S3356;
    float3  _S3368 = (*sh_coeffs_19)[int(10)] * _S3356;
    float3  _S3369 = _S3322 * _S3356;
    float3  _S3370 = (*sh_coeffs_19)[int(9)] * _S3356;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S3370.x + _S3370.y + _S3370.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S3358.x + _S3358.y + _S3358.z);
    float _S3371 = _S3368.x + _S3368.y + _S3368.z;
    float _S3372 = _S3360.x + _S3360.y + _S3360.z;
    float _S3373 = _S3366.x + _S3366.y + _S3366.z;
    float _S3374 = _S3362.x + _S3362.y + _S3362.z;
    float _S3375 = _S3364.x + _S3364.y + _S3364.z;
    float _S3376 = - s_diff_fC2_T_5;
    float3  _S3377 = _S3313 * _S3356;
    float3  _S3378 = (*sh_coeffs_19)[int(8)] * _S3356;
    float3  _S3379 = _S3311 * _S3356;
    float3  _S3380 = (*sh_coeffs_19)[int(7)] * _S3356;
    float3  _S3381 = _S3310 * _S3356;
    float3  _S3382 = (*sh_coeffs_19)[int(6)] * _S3356;
    float3  _S3383 = _S3312 * _S3356;
    float3  _S3384 = (*sh_coeffs_19)[int(5)] * _S3356;
    float3  _S3385 = _S3314 * _S3356;
    float3  _S3386 = (*sh_coeffs_19)[int(4)] * _S3356;
    float _S3387 = _S3384.x + _S3384.y + _S3384.z;
    float _S3388 = _S3380.x + _S3380.y + _S3380.z;
    float _S3389 = fTmp1B_19 * _S3371 + x_50 * s_diff_fS2_T_5 + y_22 * _S3376 + 0.54627424478530884f * (_S3386.x + _S3386.y + _S3386.z);
    float _S3390 = fTmp1B_19 * _S3372 + y_22 * s_diff_fS2_T_5 + x_50 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S3378.x + _S3378.y + _S3378.z);
    float _S3391 = y_22 * - _S3390;
    float _S3392 = x_50 * _S3390;
    float _S3393 = z_19 * (1.86588168144226074f * (z_19 * _S3375) + -2.28522896766662598f * (y_22 * _S3373 + x_50 * _S3374) + 0.94617468118667603f * (_S3382.x + _S3382.y + _S3382.z));
    float3  _S3394 = make_float3 (0.48860251903533936f) * _S3356;
    float3  _S3395 = - _S3394;
    float3  _S3396 = _S3304 * _S3395;
    float3  _S3397 = (*sh_coeffs_19)[int(3)] * _S3395;
    float3  _S3398 = _S3306 * _S3394;
    float3  _S3399 = (*sh_coeffs_19)[int(2)] * _S3394;
    float3  _S3400 = _S3308 * _S3394;
    float3  _S3401 = (*sh_coeffs_19)[int(1)] * _S3394;
    float _S3402 = (_S3315 * _S3375 + 1.44530570507049561f * (fS1_19 * _S3371 + fC1_19 * _S3372) + -1.09254848957061768f * (y_22 * _S3387 + x_50 * _S3388) + _S3393 + _S3393 + _S3399.x + _S3399.y + _S3399.z) / _S3305;
    float _S3403 = _S3303 * _S3402;
    float _S3404 = (fTmp0C_19 * _S3373 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S3376 + fTmp0B_19 * _S3387 + _S3309 * _S3389 + _S3391 + _S3391 + - (_S3401.x + _S3401.y + _S3401.z)) / _S3305;
    float _S3405 = _S3303 * _S3404;
    float _S3406 = (fTmp0C_19 * _S3374 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S3388 + 2.0f * (y_22 * _S3389) + _S3392 + _S3392 + _S3397.x + _S3397.y + _S3397.z) / _S3305;
    float _S3407 = _S3303 * _S3406;
    float _S3408 = _S3301 * - _S3402 + _S3300 * - _S3404 + _S3299 * - _S3406;
    DiffPair_float_0 _S3409;
    (&_S3409)->primal_0 = _S3302;
    (&_S3409)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3409, _S3408);
    float _S3410 = _S3301 * _S3409.differential_0;
    float _S3411 = _S3300 * _S3409.differential_0;
    float _S3412 = _S3299 * _S3409.differential_0;
    float3  _S3413 = make_float3 (0.282094806432724f) * _S3356;
    float3  _S3414 = make_float3 (_S3407 + _S3412 + _S3412, _S3405 + _S3411 + _S3411, _S3403 + _S3410 + _S3410);
    float3  _S3415 = - - _S3414;
    Matrix<float, 3, 3>  _S3416 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3417;
    (&_S3417)->primal_0 = _S3297;
    (&_S3417)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3418;
    (&_S3418)->primal_0 = t_22;
    (&_S3418)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3417, &_S3418, _S3415);
    Matrix<float, 3, 3>  _S3419 = transpose_0(_S3417.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3420;
    (&_S3420)->primal_0 = _S3296;
    (&_S3420)->differential_0 = _S3336;
    s_bwd_prop_log_1(&_S3420, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3421;
    (&_S3421)->primal_0 = vert2_c_4;
    (&_S3421)->differential_0 = _S3336;
    s_bwd_length_impl_0(&_S3421, _S3420.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3422;
    (&_S3422)->primal_0 = vert1_c_4;
    (&_S3422)->differential_0 = _S3336;
    s_bwd_length_impl_0(&_S3422, _S3420.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3423;
    (&_S3423)->primal_0 = vert0_c_4;
    (&_S3423)->differential_0 = _S3336;
    s_bwd_length_impl_0(&_S3423, _S3420.differential_0.x);
    DiffPair_float_0 _S3424;
    (&_S3424)->primal_0 = _S3295;
    (&_S3424)->differential_0 = 0.0f;
    DiffPair_float_0 _S3425;
    (&_S3425)->primal_0 = _S3282;
    (&_S3425)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3424, &_S3425, 0.0f);
    DiffPair_float_0 _S3426;
    (&_S3426)->primal_0 = _S3246;
    (&_S3426)->differential_0 = 0.0f;
    DiffPair_float_0 _S3427;
    (&_S3427)->primal_0 = _S3264;
    (&_S3427)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3426, &_S3427, _S3424.differential_0);
    DiffPair_float_0 _S3428;
    (&_S3428)->primal_0 = _S3294;
    (&_S3428)->differential_0 = 0.0f;
    DiffPair_float_0 _S3429;
    (&_S3429)->primal_0 = _S3282;
    (&_S3429)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3428, &_S3429, 0.0f);
    DiffPair_float_0 _S3430;
    (&_S3430)->primal_0 = _S3246;
    (&_S3430)->differential_0 = 0.0f;
    DiffPair_float_0 _S3431;
    (&_S3431)->primal_0 = _S3264;
    (&_S3431)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3430, &_S3431, _S3428.differential_0);
    DiffPair_float_0 _S3432;
    (&_S3432)->primal_0 = _S3293;
    (&_S3432)->differential_0 = 0.0f;
    DiffPair_float_0 _S3433;
    (&_S3433)->primal_0 = _S3281;
    (&_S3433)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3432, &_S3433, 0.0f);
    DiffPair_float_0 _S3434;
    (&_S3434)->primal_0 = _S3245;
    (&_S3434)->differential_0 = 0.0f;
    DiffPair_float_0 _S3435;
    (&_S3435)->primal_0 = _S3263;
    (&_S3435)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3434, &_S3435, _S3432.differential_0);
    DiffPair_float_0 _S3436;
    (&_S3436)->primal_0 = _S3292;
    (&_S3436)->differential_0 = 0.0f;
    DiffPair_float_0 _S3437;
    (&_S3437)->primal_0 = _S3281;
    (&_S3437)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3436, &_S3437, 0.0f);
    DiffPair_float_0 _S3438;
    (&_S3438)->primal_0 = _S3245;
    (&_S3438)->differential_0 = 0.0f;
    DiffPair_float_0 _S3439;
    (&_S3439)->primal_0 = _S3263;
    (&_S3439)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3438, &_S3439, _S3436.differential_0);
    DiffPair_float_0 _S3440;
    (&_S3440)->primal_0 = _S3290;
    (&_S3440)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3440, 0.0f);
    float _S3441 = - (-1.0f * - (_S3440.differential_0 / _S3291));
    float2  _S3442 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3443;
    (&_S3443)->primal_0 = e2_0;
    (&_S3443)->differential_0 = _S3442;
    s_bwd_length_impl_1(&_S3443, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3444;
    (&_S3444)->primal_0 = e1_4;
    (&_S3444)->differential_0 = _S3442;
    s_bwd_length_impl_1(&_S3444, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3445;
    (&_S3445)->primal_0 = e0_4;
    (&_S3445)->differential_0 = _S3442;
    s_bwd_length_impl_1(&_S3445, -0.0f);
    DiffPair_float_0 _S3446;
    (&_S3446)->primal_0 = _S3288;
    (&_S3446)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3446, 0.0f);
    float _S3447 = - _S3446.differential_0;
    float2  _S3448 = _S3444.differential_0 + make_float2 (_S3286 * _S3447, _S3284 * _S3446.differential_0);
    float2  _S3449 = _S3445.differential_0 + make_float2 (_S3285 * _S3446.differential_0, _S3287 * _S3447);
    float2  _S3450 = v_uv2_0 + - _S3443.differential_0 + _S3448;
    float _S3451 = fx_28 * (_S3433.differential_0 + _S3437.differential_0 + _S3450.x);
    float2  _S3452 = make_float2 (_S3451, fy_28 * (_S3425.differential_0 + _S3429.differential_0 + _S3450.y)) + make_float2 ((*dist_coeffs_34)[int(8)] * _S3451, (*dist_coeffs_34)[int(9)] * _S3451);
    float2  _S3453 = _S3269 * _S3452;
    float _S3454 = (*dist_coeffs_34)[int(4)] * _S3452.y;
    float _S3455 = (*dist_coeffs_34)[int(5)] * _S3452.x;
    float _S3456 = _S3453.x + _S3453.y;
    float _S3457 = r2_58 * _S3456;
    float _S3458 = r2_58 * _S3457;
    float _S3459 = (*dist_coeffs_34)[int(7)] * _S3452.y + _S3454 + (*dist_coeffs_34)[int(6)] * _S3452.x + _S3455 + _S3273 * _S3456 + _S3272 * _S3457 + _S3271 * _S3458 + (*dist_coeffs_34)[int(3)] * (r2_58 * _S3458);
    float _S3460 = v_58 * _S3459;
    float _S3461 = u_58 * _S3459;
    float2  _S3462 = (_S3274 * _S3452 + make_float2 (_S3240 * (v_58 * _S3452.y) + _S3276 * _S3455 + 2.0f * (u_58 * _S3455) + _S3237 * (v_58 * _S3452.x) + _S3461 + _S3461, _S3278 * _S3454 + 2.0f * (v_58 * _S3454) + _S3277 * _S3452.y + _S3275 * _S3452.x + _S3460 + _S3460)) / _S3270;
    float2  _S3463 = _S3266 * - _S3462;
    float2  _S3464 = _S3268 * _S3462;
    float2  _S3465 = v_uv1_0 + - _S3448 + _S3449;
    float _S3466 = fx_28 * (_S3435.differential_0 + _S3439.differential_0 + _S3465.x);
    float2  _S3467 = make_float2 (_S3466, fy_28 * (_S3427.differential_0 + _S3431.differential_0 + _S3465.y)) + make_float2 ((*dist_coeffs_34)[int(8)] * _S3466, (*dist_coeffs_34)[int(9)] * _S3466);
    float2  _S3468 = _S3251 * _S3467;
    float _S3469 = (*dist_coeffs_34)[int(4)] * _S3467.y;
    float _S3470 = (*dist_coeffs_34)[int(5)] * _S3467.x;
    float _S3471 = _S3468.x + _S3468.y;
    float _S3472 = r2_57 * _S3471;
    float _S3473 = r2_57 * _S3472;
    float _S3474 = (*dist_coeffs_34)[int(7)] * _S3467.y + _S3469 + (*dist_coeffs_34)[int(6)] * _S3467.x + _S3470 + _S3255 * _S3471 + _S3254 * _S3472 + _S3253 * _S3473 + (*dist_coeffs_34)[int(3)] * (r2_57 * _S3473);
    float _S3475 = v_57 * _S3474;
    float _S3476 = u_57 * _S3474;
    float2  _S3477 = (_S3256 * _S3467 + make_float2 (_S3240 * (v_57 * _S3467.y) + _S3258 * _S3470 + 2.0f * (u_57 * _S3470) + _S3237 * (v_57 * _S3467.x) + _S3476 + _S3476, _S3260 * _S3469 + 2.0f * (v_57 * _S3469) + _S3259 * _S3467.y + _S3257 * _S3467.x + _S3475 + _S3475)) / _S3252;
    float2  _S3478 = _S3248 * - _S3477;
    float2  _S3479 = _S3250 * _S3477;
    float _S3480 = _S3478.x + _S3478.y;
    float2  _S3481 = v_uv0_0 + _S3443.differential_0 + - _S3449;
    float _S3482 = fx_28 * (_S3434.differential_0 + _S3438.differential_0 + _S3481.x);
    float2  _S3483 = make_float2 (_S3482, fy_28 * (_S3426.differential_0 + _S3430.differential_0 + _S3481.y)) + make_float2 ((*dist_coeffs_34)[int(8)] * _S3482, (*dist_coeffs_34)[int(9)] * _S3482);
    float2  _S3484 = _S3231 * _S3483;
    float _S3485 = (*dist_coeffs_34)[int(4)] * _S3483.y;
    float _S3486 = (*dist_coeffs_34)[int(5)] * _S3483.x;
    float _S3487 = _S3484.x + _S3484.y;
    float _S3488 = r2_56 * _S3487;
    float _S3489 = r2_56 * _S3488;
    float _S3490 = (*dist_coeffs_34)[int(7)] * _S3483.y + _S3485 + (*dist_coeffs_34)[int(6)] * _S3483.x + _S3486 + _S3235 * _S3487 + _S3234 * _S3488 + _S3233 * _S3489 + (*dist_coeffs_34)[int(3)] * (r2_56 * _S3489);
    float _S3491 = v_56 * _S3490;
    float _S3492 = u_56 * _S3490;
    float2  _S3493 = (_S3236 * _S3483 + make_float2 (_S3240 * (v_56 * _S3483.y) + _S3239 * _S3486 + 2.0f * (u_56 * _S3486) + _S3237 * (v_56 * _S3483.x) + _S3492 + _S3492, _S3242 * _S3485 + 2.0f * (v_56 * _S3485) + _S3241 * _S3483.y + _S3238 * _S3483.x + _S3491 + _S3491)) / _S3232;
    float2  _S3494 = _S3228 * - _S3493;
    float2  _S3495 = _S3230 * _S3493;
    float _S3496 = _S3494.x + _S3494.y;
    float3  _S3497 = _S3342.differential_0 + _S3421.differential_0 + make_float3 (_S3464.x, _S3464.y, _S3463.x + _S3463.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3498;
    (&_S3498)->primal_0 = R_23;
    (&_S3498)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3499;
    (&_S3499)->primal_0 = vert2_0;
    (&_S3499)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3498, &_S3499, _S3497);
    float3  _S3500 = _S3341.differential_0 + _S3422.differential_0 + make_float3 (_S3479.x, _S3479.y, _S3480);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3501;
    (&_S3501)->primal_0 = R_23;
    (&_S3501)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3502;
    (&_S3502)->primal_0 = vert1_0;
    (&_S3502)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3501, &_S3502, _S3500);
    float3  _S3503 = _S3343 + _S3344 + _S3423.differential_0 + make_float3 (_S3495.x, _S3495.y, _S3496);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3504;
    (&_S3504)->primal_0 = R_23;
    (&_S3504)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3505;
    (&_S3505)->primal_0 = vert0_0;
    (&_S3505)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3504, &_S3505, _S3503);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3506;
    (&_S3506)->primal_0 = _S3222;
    (&_S3506)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3507;
    (&_S3507)->primal_0 = _S3227;
    (&_S3507)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3506, &_S3507, _S3499.differential_0);
    float _S3508 = - _S3507.differential_0.y;
    float _S3509 = _S3226 * _S3507.differential_0.x;
    float _S3510 = - (_S3218 * _S3507.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3511;
    (&_S3511)->primal_0 = _S3222;
    (&_S3511)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3512;
    (&_S3512)->primal_0 = _S3225;
    (&_S3512)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3511, &_S3512, _S3502.differential_0);
    float _S3513 = _S3218 * _S3512.differential_0.x;
    float _S3514 = _S3224 * _S3512.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3515;
    (&_S3515)->primal_0 = _S3222;
    (&_S3515)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3516;
    (&_S3516)->primal_0 = _S3223;
    (&_S3516)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3515, &_S3516, _S3505.differential_0);
    Matrix<float, 3, 3>  _S3517 = transpose_0(_S3506.differential_0 + _S3511.differential_0 + _S3515.differential_0);
    float _S3518 = 2.0f * - _S3517.rows[int(2)].z;
    float _S3519 = 2.0f * _S3517.rows[int(2)].y;
    float _S3520 = 2.0f * _S3517.rows[int(2)].x;
    float _S3521 = 2.0f * _S3517.rows[int(1)].z;
    float _S3522 = 2.0f * - _S3517.rows[int(1)].y;
    float _S3523 = 2.0f * _S3517.rows[int(1)].x;
    float _S3524 = 2.0f * _S3517.rows[int(0)].z;
    float _S3525 = 2.0f * _S3517.rows[int(0)].y;
    float _S3526 = 2.0f * - _S3517.rows[int(0)].x;
    float _S3527 = - _S3523 + _S3525;
    float _S3528 = _S3520 + - _S3524;
    float _S3529 = - _S3519 + _S3521;
    float _S3530 = _S3519 + _S3521;
    float _S3531 = _S3520 + _S3524;
    float _S3532 = _S3523 + _S3525;
    float _S3533 = quat_24.w * (_S3522 + _S3526);
    float _S3534 = quat_24.z * (_S3518 + _S3526);
    float _S3535 = quat_24.y * (_S3518 + _S3522);
    float _S3536 = quat_24.x * _S3527 + quat_24.z * _S3530 + quat_24.y * _S3531 + _S3533 + _S3533;
    float _S3537 = quat_24.x * _S3528 + quat_24.w * _S3530 + quat_24.y * _S3532 + _S3534 + _S3534;
    float _S3538 = quat_24.x * _S3529 + quat_24.w * _S3531 + quat_24.z * _S3532 + _S3535 + _S3535;
    float _S3539 = quat_24.w * _S3527 + quat_24.z * _S3528 + quat_24.y * _S3529;
    float _S3540 = _S3510 + _S3513;
    float _S3541 = 0.5f * - _S3540;
    float _S3542 = _S3508 + _S3512.differential_0.y;
    DiffPair_float_0 _S3543;
    (&_S3543)->primal_0 = _S3219;
    (&_S3543)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3543, _S3542);
    float _S3544 = _S3541 + _S3543.differential_0;
    float _S3545 = _S3509 + _S3514 + _S3516.differential_0.x;
    DiffPair_float_0 _S3546;
    (&_S3546)->primal_0 = _S3217;
    (&_S3546)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3546, _S3545);
    float _S3547 = _S3541 + _S3546.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3548;
    (&_S3548)->primal_0 = R_23;
    (&_S3548)->differential_0 = _S3416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3549;
    (&_S3549)->primal_0 = mean_23;
    (&_S3549)->differential_0 = _S3336;
    s_bwd_prop_mul_1(&_S3548, &_S3549, _S3338.differential_0);
    float3  _S3550 = _S3418.differential_0 + _S3497 + _S3500 + _S3503 + _S3338.differential_0;
    Matrix<float, 3, 3>  _S3551 = _S3419 + _S3498.differential_0 + _S3501.differential_0 + _S3504.differential_0 + _S3548.differential_0;
    FixedArray<float3 , 2>  _S3552;
    _S3552[int(0)] = _S3336;
    _S3552[int(1)] = _S3336;
    _S3552[int(1)] = _S3350;
    _S3552[int(0)] = _S3355;
    FixedArray<float3 , 16>  _S3553;
    _S3553[int(0)] = _S3336;
    _S3553[int(1)] = _S3336;
    _S3553[int(2)] = _S3336;
    _S3553[int(3)] = _S3336;
    _S3553[int(4)] = _S3336;
    _S3553[int(5)] = _S3336;
    _S3553[int(6)] = _S3336;
    _S3553[int(7)] = _S3336;
    _S3553[int(8)] = _S3336;
    _S3553[int(9)] = _S3336;
    _S3553[int(10)] = _S3336;
    _S3553[int(11)] = _S3336;
    _S3553[int(12)] = _S3336;
    _S3553[int(13)] = _S3336;
    _S3553[int(14)] = _S3336;
    _S3553[int(15)] = _S3336;
    _S3553[int(15)] = _S3357;
    _S3553[int(14)] = _S3359;
    _S3553[int(13)] = _S3361;
    _S3553[int(12)] = _S3363;
    _S3553[int(11)] = _S3365;
    _S3553[int(10)] = _S3367;
    _S3553[int(9)] = _S3369;
    _S3553[int(8)] = _S3377;
    _S3553[int(7)] = _S3379;
    _S3553[int(6)] = _S3381;
    _S3553[int(5)] = _S3383;
    _S3553[int(4)] = _S3385;
    _S3553[int(3)] = _S3396;
    _S3553[int(2)] = _S3398;
    _S3553[int(1)] = _S3400;
    _S3553[int(0)] = _S3413;
    float2  _S3554 = v_out_hardness_0 + make_float2 (0.0f, _S3441);
    float3  _S3555 = make_float3 (_S3547, _S3544, _S3540);
    float4  _S3556 = make_float4 (0.0f);
    *&((&_S3556)->w) = _S3536;
    *&((&_S3556)->z) = _S3537;
    *&((&_S3556)->y) = _S3538;
    *&((&_S3556)->x) = _S3539;
    *v_mean_7 = _S3414 + _S3499.differential_0 + _S3502.differential_0 + _S3505.differential_0 + _S3549.differential_0;
    *v_quat_6 = _S3556;
    *v_scale_6 = _S3555;
    *v_hardness_0 = _S3554;
    *v_sh_coeffs_5 = _S3553;
    *v_ch_coeffs_0 = _S3552;
    *v_R_6 = _S3551;
    *v_t_6 = _S3550;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_24, float3  t_23, float fx_29, float fy_29, float cx_24, float cy_24, FixedArray<float, 10>  * dist_coeffs_35, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_1(R_24, mean_24) + t_23;
    float _S3557 = scale_24.x;
    float _S3558 = s_primal_ctx_exp_1(_S3557);
    float _S3559 = scale_24.y;
    float _S3560 = s_primal_ctx_exp_1(_S3559);
    float sz_5 = scale_24.z - 0.5f * (_S3557 + _S3559);
    float _S3561 = quat_25.y;
    float x2_25 = _S3561 * _S3561;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_45 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S3562 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_45), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_45), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S3563 = make_float3 (_S3558, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S3562, _S3563) + mean_24;
    float _S3564 = -0.5f + sz_5;
    float3  _S3565 = make_float3 (_S3558 * _S3564, _S3560, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S3562, _S3565) + mean_24;
    float _S3566 = -0.5f - sz_5;
    float3  _S3567 = make_float3 (_S3558 * _S3566, - _S3560, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S3562, _S3567) + mean_24;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_24, vert0_1) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_24, vert1_1) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_24, vert2_1) + t_23;
    float _S3568 = length_1(vert0_c_5);
    float _S3569 = length_1(vert1_c_5);
    float _S3570 = length_1(vert2_c_5);
    float2  _S3571 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S3572 = length_0(_S3571);
    float _S3573 = vert0_c_5.z;
    float _S3574 = s_primal_ctx_atan2_0(_S3572, _S3573);
    bool _S3575 = _S3574 < 0.00100000004749745f;
    float k_9;
    float _S3576;
    float _S3577;
    float _S3578;
    if(_S3575)
    {
        float _S3579 = 1.0f - _S3574 * _S3574 / 3.0f;
        float _S3580 = _S3573 * _S3573;
        k_9 = _S3579 / _S3573;
        _S3576 = 0.0f;
        _S3577 = _S3580;
        _S3578 = _S3579;
    }
    else
    {
        float _S3581 = _S3572 * _S3572;
        k_9 = _S3574 / _S3572;
        _S3576 = _S3581;
        _S3577 = 0.0f;
        _S3578 = 0.0f;
    }
    float2  _S3582 = make_float2 (k_9);
    float2  _S3583 = _S3571 * make_float2 (k_9);
    float u_59 = _S3583.x;
    float v_59 = _S3583.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S3584 = (*dist_coeffs_35)[int(2)] + r2_59 * (*dist_coeffs_35)[int(3)];
    float _S3585 = (*dist_coeffs_35)[int(1)] + r2_59 * _S3584;
    float _S3586 = (*dist_coeffs_35)[int(0)] + r2_59 * _S3585;
    float radial_5 = 1.0f + r2_59 * _S3586;
    float2  _S3587 = make_float2 (radial_5);
    float _S3588 = 2.0f * (*dist_coeffs_35)[int(4)];
    float _S3589 = _S3588 * u_59;
    float _S3590 = 2.0f * u_59;
    float _S3591 = 2.0f * (*dist_coeffs_35)[int(5)];
    float _S3592 = _S3591 * u_59;
    float _S3593 = 2.0f * v_59;
    float2  _S3594 = _S3583 * make_float2 (radial_5) + make_float2 (_S3589 * v_59 + (*dist_coeffs_35)[int(5)] * (r2_59 + _S3590 * u_59) + (*dist_coeffs_35)[int(6)] * r2_59, _S3592 * v_59 + (*dist_coeffs_35)[int(4)] * (r2_59 + _S3593 * v_59) + (*dist_coeffs_35)[int(7)] * r2_59);
    float2  _S3595 = _S3594 + make_float2 ((*dist_coeffs_35)[int(8)] * _S3594.x + (*dist_coeffs_35)[int(9)] * _S3594.y, 0.0f);
    float _S3596 = fx_29 * _S3595.x + cx_24;
    float _S3597 = fy_29 * _S3595.y + cy_24;
    float2  _S3598 = make_float2 (_S3596, _S3597);
    float2  _S3599 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S3600 = length_0(_S3599);
    float _S3601 = vert1_c_5.z;
    float _S3602 = s_primal_ctx_atan2_0(_S3600, _S3601);
    bool _S3603 = _S3602 < 0.00100000004749745f;
    float _S3604;
    float _S3605;
    float _S3606;
    if(_S3603)
    {
        float _S3607 = 1.0f - _S3602 * _S3602 / 3.0f;
        float _S3608 = _S3601 * _S3601;
        k_9 = _S3607 / _S3601;
        _S3604 = 0.0f;
        _S3605 = _S3608;
        _S3606 = _S3607;
    }
    else
    {
        float _S3609 = _S3600 * _S3600;
        k_9 = _S3602 / _S3600;
        _S3604 = _S3609;
        _S3605 = 0.0f;
        _S3606 = 0.0f;
    }
    float2  _S3610 = make_float2 (k_9);
    float2  _S3611 = _S3599 * make_float2 (k_9);
    float u_60 = _S3611.x;
    float v_60 = _S3611.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S3612 = (*dist_coeffs_35)[int(2)] + r2_60 * (*dist_coeffs_35)[int(3)];
    float _S3613 = (*dist_coeffs_35)[int(1)] + r2_60 * _S3612;
    float _S3614 = (*dist_coeffs_35)[int(0)] + r2_60 * _S3613;
    float radial_6 = 1.0f + r2_60 * _S3614;
    float2  _S3615 = make_float2 (radial_6);
    float _S3616 = _S3588 * u_60;
    float _S3617 = 2.0f * u_60;
    float _S3618 = _S3591 * u_60;
    float _S3619 = 2.0f * v_60;
    float2  _S3620 = _S3611 * make_float2 (radial_6) + make_float2 (_S3616 * v_60 + (*dist_coeffs_35)[int(5)] * (r2_60 + _S3617 * u_60) + (*dist_coeffs_35)[int(6)] * r2_60, _S3618 * v_60 + (*dist_coeffs_35)[int(4)] * (r2_60 + _S3619 * v_60) + (*dist_coeffs_35)[int(7)] * r2_60);
    float2  _S3621 = _S3620 + make_float2 ((*dist_coeffs_35)[int(8)] * _S3620.x + (*dist_coeffs_35)[int(9)] * _S3620.y, 0.0f);
    float _S3622 = fx_29 * _S3621.x + cx_24;
    float _S3623 = fy_29 * _S3621.y + cy_24;
    float2  _S3624 = make_float2 (_S3622, _S3623);
    float2  _S3625 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S3626 = length_0(_S3625);
    float _S3627 = vert2_c_5.z;
    float _S3628 = s_primal_ctx_atan2_0(_S3626, _S3627);
    bool _S3629 = _S3628 < 0.00100000004749745f;
    float _S3630;
    float _S3631;
    float _S3632;
    if(_S3629)
    {
        float _S3633 = 1.0f - _S3628 * _S3628 / 3.0f;
        float _S3634 = _S3627 * _S3627;
        k_9 = _S3633 / _S3627;
        _S3630 = 0.0f;
        _S3631 = _S3634;
        _S3632 = _S3633;
    }
    else
    {
        float _S3635 = _S3626 * _S3626;
        k_9 = _S3628 / _S3626;
        _S3630 = _S3635;
        _S3631 = 0.0f;
        _S3632 = 0.0f;
    }
    float2  _S3636 = make_float2 (k_9);
    float2  _S3637 = _S3625 * make_float2 (k_9);
    float u_61 = _S3637.x;
    float v_61 = _S3637.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S3638 = (*dist_coeffs_35)[int(2)] + r2_61 * (*dist_coeffs_35)[int(3)];
    float _S3639 = (*dist_coeffs_35)[int(1)] + r2_61 * _S3638;
    float _S3640 = (*dist_coeffs_35)[int(0)] + r2_61 * _S3639;
    float radial_7 = 1.0f + r2_61 * _S3640;
    float2  _S3641 = make_float2 (radial_7);
    float _S3642 = _S3588 * u_61;
    float _S3643 = 2.0f * u_61;
    float _S3644 = _S3591 * u_61;
    float _S3645 = 2.0f * v_61;
    float2  _S3646 = _S3637 * make_float2 (radial_7) + make_float2 (_S3642 * v_61 + (*dist_coeffs_35)[int(5)] * (r2_61 + _S3643 * u_61) + (*dist_coeffs_35)[int(6)] * r2_61, _S3644 * v_61 + (*dist_coeffs_35)[int(4)] * (r2_61 + _S3645 * v_61) + (*dist_coeffs_35)[int(7)] * r2_61);
    float2  _S3647 = _S3646 + make_float2 ((*dist_coeffs_35)[int(8)] * _S3646.x + (*dist_coeffs_35)[int(9)] * _S3646.y, 0.0f);
    float _S3648 = fx_29 * _S3647.x + cx_24;
    float _S3649 = fy_29 * _S3647.y + cy_24;
    float2  _S3650 = make_float2 (_S3648, _S3649);
    float2  e0_5 = _S3624 - _S3598;
    float2  e1_5 = _S3650 - _S3624;
    float2  e2_1 = _S3598 - _S3650;
    float _S3651 = e0_5.x;
    float _S3652 = e1_5.y;
    float _S3653 = e0_5.y;
    float _S3654 = e1_5.x;
    float _S3655 = _S3651 * _S3652 - _S3653 * _S3654;
    float _S3656 = 1.0f - hardness_5.y;
    float _S3657 = -1.0f / _S3656;
    float _S3658 = _S3656 * _S3656;
    float _S3659 = s_primal_ctx_max_0(_S3596, _S3622);
    float _S3660 = s_primal_ctx_min_0(_S3596, _S3622);
    float _S3661 = s_primal_ctx_max_0(_S3597, _S3623);
    float _S3662 = s_primal_ctx_min_0(_S3597, _S3623);
    float3  _S3663 = make_float3 (_S3568, _S3569, _S3570) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3664 = transpose_0(R_24);
    float3  _S3665 = mean_24 - - s_primal_ctx_mul_1(_S3664, t_23);
    float _S3666 = _S3665.x;
    float _S3667 = _S3665.y;
    float _S3668 = _S3665.z;
    float _S3669 = _S3666 * _S3666 + _S3667 * _S3667 + _S3668 * _S3668;
    float _S3670 = s_primal_ctx_sqrt_0(_S3669);
    float x_51 = _S3666 / _S3670;
    float3  _S3671 = make_float3 (x_51);
    float _S3672 = _S3670 * _S3670;
    float y_23 = _S3667 / _S3670;
    float z_20 = _S3668 / _S3670;
    float3  _S3673 = make_float3 (z_20);
    float _S3674 = - y_23;
    float3  _S3675 = make_float3 (_S3674);
    float z2_46 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_51 * x_51 - y_23 * y_23;
    float _S3676 = 2.0f * x_51;
    float fS1_20 = _S3676 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_46 - 0.31539157032966614f;
    float3  _S3677 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_51;
    float3  _S3678 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3679 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3680 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3681 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_46 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3682 = 1.86588168144226074f * z2_46 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3682;
    float3  _S3683 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_51;
    float3  _S3684 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3685 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3686 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3687 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_51 * fC1_20 - y_23 * fS1_20);
    float3  _S3688 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_51 * fS1_20 + y_23 * fC1_20);
    float3  _S3689 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3674) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_51) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3690 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3691 = make_float3 (0.0f);
    float3  _S3692 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3693 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3694 = make_float3 (_S3693);
    float3  _S3695 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3693);
    float3  _S3696 = _S3692 + _S3695 + make_float3 (0.5f);
    float3  _S3697 = _S3692 - _S3695 + make_float3 (0.5f);
    float3  _S3698 = vert1_c_5 - vert0_c_5;
    float3  _S3699 = vert2_c_5 - vert0_c_5;
    float3  _S3700 = s_primal_ctx_cross_0(_S3698, _S3699);
    float3  _S3701 = normalize_0(_S3700);
    float3  _S3702 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3701, mean_c_20)))))) * v_normal_1;
    float3  _S3703 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3704;
    (&_S3704)->primal_0 = _S3701;
    (&_S3704)->differential_0 = _S3703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3705;
    (&_S3705)->primal_0 = mean_c_20;
    (&_S3705)->differential_0 = _S3703;
    s_bwd_prop_dot_0(&_S3704, &_S3705, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3706 = _S3705;
    float3  _S3707 = _S3702 + _S3704.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3708;
    (&_S3708)->primal_0 = _S3700;
    (&_S3708)->differential_0 = _S3703;
    s_bwd_normalize_impl_0(&_S3708, _S3707);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3709;
    (&_S3709)->primal_0 = _S3698;
    (&_S3709)->differential_0 = _S3703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3710;
    (&_S3710)->primal_0 = _S3699;
    (&_S3710)->differential_0 = _S3703;
    s_bwd_prop_cross_0(&_S3709, &_S3710, _S3708.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3711 = _S3709;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3712 = _S3710;
    float3  _S3713 = - _S3710.differential_0;
    float3  _S3714 = - _S3709.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3715;
    (&_S3715)->primal_0 = _S3697;
    (&_S3715)->differential_0 = _S3703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3716;
    (&_S3716)->primal_0 = _S3691;
    (&_S3716)->differential_0 = _S3703;
    s_bwd_prop_max_0(&_S3715, &_S3716, (*v_rgb_7)[int(2)]);
    float3  _S3717 = - _S3715.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3718;
    (&_S3718)->primal_0 = _S3696;
    (&_S3718)->differential_0 = _S3703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3719;
    (&_S3719)->primal_0 = _S3691;
    (&_S3719)->differential_0 = _S3703;
    s_bwd_prop_max_0(&_S3718, &_S3719, (*v_rgb_7)[int(1)]);
    float3  _S3720 = _S3694 * (_S3717 + _S3718.differential_0);
    float3  _S3721 = _S3715.differential_0 + _S3718.differential_0;
    float3  _S3722 = make_float3 (0.5f) * - _S3721;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3723;
    (&_S3723)->primal_0 = _S3690;
    (&_S3723)->differential_0 = _S3703;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3724;
    (&_S3724)->primal_0 = _S3691;
    (&_S3724)->differential_0 = _S3703;
    s_bwd_prop_max_0(&_S3723, &_S3724, (*v_rgb_7)[int(0)]);
    float3  _S3725 = _S3722 + _S3723.differential_0;
    float3  _S3726 = _S3721 + _S3723.differential_0;
    float3  _S3727 = _S3688 * _S3726;
    float3  _S3728 = (*sh_coeffs_20)[int(15)] * _S3726;
    float3  _S3729 = _S3686 * _S3726;
    float3  _S3730 = (*sh_coeffs_20)[int(14)] * _S3726;
    float3  _S3731 = _S3684 * _S3726;
    float3  _S3732 = (*sh_coeffs_20)[int(13)] * _S3726;
    float3  _S3733 = _S3683 * _S3726;
    float3  _S3734 = (*sh_coeffs_20)[int(12)] * _S3726;
    float3  _S3735 = _S3685 * _S3726;
    float3  _S3736 = (*sh_coeffs_20)[int(11)] * _S3726;
    float3  _S3737 = _S3687 * _S3726;
    float3  _S3738 = (*sh_coeffs_20)[int(10)] * _S3726;
    float3  _S3739 = _S3689 * _S3726;
    float3  _S3740 = (*sh_coeffs_20)[int(9)] * _S3726;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3740.x + _S3740.y + _S3740.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3728.x + _S3728.y + _S3728.z);
    float _S3741 = _S3738.x + _S3738.y + _S3738.z;
    float _S3742 = _S3730.x + _S3730.y + _S3730.z;
    float _S3743 = _S3736.x + _S3736.y + _S3736.z;
    float _S3744 = _S3732.x + _S3732.y + _S3732.z;
    float _S3745 = _S3734.x + _S3734.y + _S3734.z;
    float _S3746 = - s_diff_fC2_T_6;
    float3  _S3747 = _S3680 * _S3726;
    float3  _S3748 = (*sh_coeffs_20)[int(8)] * _S3726;
    float3  _S3749 = _S3678 * _S3726;
    float3  _S3750 = (*sh_coeffs_20)[int(7)] * _S3726;
    float3  _S3751 = _S3677 * _S3726;
    float3  _S3752 = (*sh_coeffs_20)[int(6)] * _S3726;
    float3  _S3753 = _S3679 * _S3726;
    float3  _S3754 = (*sh_coeffs_20)[int(5)] * _S3726;
    float3  _S3755 = _S3681 * _S3726;
    float3  _S3756 = (*sh_coeffs_20)[int(4)] * _S3726;
    float _S3757 = _S3754.x + _S3754.y + _S3754.z;
    float _S3758 = _S3750.x + _S3750.y + _S3750.z;
    float _S3759 = fTmp1B_20 * _S3741 + x_51 * s_diff_fS2_T_6 + y_23 * _S3746 + 0.54627424478530884f * (_S3756.x + _S3756.y + _S3756.z);
    float _S3760 = fTmp1B_20 * _S3742 + y_23 * s_diff_fS2_T_6 + x_51 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3748.x + _S3748.y + _S3748.z);
    float _S3761 = y_23 * - _S3760;
    float _S3762 = x_51 * _S3760;
    float _S3763 = z_20 * (1.86588168144226074f * (z_20 * _S3745) + -2.28522896766662598f * (y_23 * _S3743 + x_51 * _S3744) + 0.94617468118667603f * (_S3752.x + _S3752.y + _S3752.z));
    float3  _S3764 = make_float3 (0.48860251903533936f) * _S3726;
    float3  _S3765 = - _S3764;
    float3  _S3766 = _S3671 * _S3765;
    float3  _S3767 = (*sh_coeffs_20)[int(3)] * _S3765;
    float3  _S3768 = _S3673 * _S3764;
    float3  _S3769 = (*sh_coeffs_20)[int(2)] * _S3764;
    float3  _S3770 = _S3675 * _S3764;
    float3  _S3771 = (*sh_coeffs_20)[int(1)] * _S3764;
    float _S3772 = (_S3682 * _S3745 + 1.44530570507049561f * (fS1_20 * _S3741 + fC1_20 * _S3742) + -1.09254848957061768f * (y_23 * _S3757 + x_51 * _S3758) + _S3763 + _S3763 + _S3769.x + _S3769.y + _S3769.z) / _S3672;
    float _S3773 = _S3670 * _S3772;
    float _S3774 = (fTmp0C_20 * _S3743 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3746 + fTmp0B_20 * _S3757 + _S3676 * _S3759 + _S3761 + _S3761 + - (_S3771.x + _S3771.y + _S3771.z)) / _S3672;
    float _S3775 = _S3670 * _S3774;
    float _S3776 = (fTmp0C_20 * _S3744 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3758 + 2.0f * (y_23 * _S3759) + _S3762 + _S3762 + _S3767.x + _S3767.y + _S3767.z) / _S3672;
    float _S3777 = _S3670 * _S3776;
    float _S3778 = _S3668 * - _S3772 + _S3667 * - _S3774 + _S3666 * - _S3776;
    DiffPair_float_0 _S3779;
    (&_S3779)->primal_0 = _S3669;
    (&_S3779)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3779, _S3778);
    float _S3780 = _S3668 * _S3779.differential_0;
    float _S3781 = _S3667 * _S3779.differential_0;
    float _S3782 = _S3666 * _S3779.differential_0;
    float3  _S3783 = make_float3 (0.282094806432724f) * _S3726;
    float3  _S3784 = make_float3 (_S3777 + _S3782 + _S3782, _S3775 + _S3781 + _S3781, _S3773 + _S3780 + _S3780);
    float3  _S3785 = - - _S3784;
    Matrix<float, 3, 3>  _S3786 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3787;
    (&_S3787)->primal_0 = _S3664;
    (&_S3787)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3788;
    (&_S3788)->primal_0 = t_23;
    (&_S3788)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3787, &_S3788, _S3785);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3789 = _S3788;
    Matrix<float, 3, 3>  _S3790 = transpose_0(_S3787.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3791;
    (&_S3791)->primal_0 = _S3663;
    (&_S3791)->differential_0 = _S3703;
    s_bwd_prop_log_1(&_S3791, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3792 = _S3791;
    DiffPair_float_0 _S3793;
    (&_S3793)->primal_0 = _S3662;
    (&_S3793)->differential_0 = 0.0f;
    DiffPair_float_0 _S3794;
    (&_S3794)->primal_0 = _S3649;
    (&_S3794)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3793, &_S3794, 0.0f);
    DiffPair_float_0 _S3795;
    (&_S3795)->primal_0 = _S3597;
    (&_S3795)->differential_0 = 0.0f;
    DiffPair_float_0 _S3796;
    (&_S3796)->primal_0 = _S3623;
    (&_S3796)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3795, &_S3796, _S3793.differential_0);
    DiffPair_float_0 _S3797;
    (&_S3797)->primal_0 = _S3661;
    (&_S3797)->differential_0 = 0.0f;
    DiffPair_float_0 _S3798;
    (&_S3798)->primal_0 = _S3649;
    (&_S3798)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3797, &_S3798, 0.0f);
    DiffPair_float_0 _S3799;
    (&_S3799)->primal_0 = _S3597;
    (&_S3799)->differential_0 = 0.0f;
    DiffPair_float_0 _S3800;
    (&_S3800)->primal_0 = _S3623;
    (&_S3800)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3799, &_S3800, _S3797.differential_0);
    DiffPair_float_0 _S3801;
    (&_S3801)->primal_0 = _S3660;
    (&_S3801)->differential_0 = 0.0f;
    DiffPair_float_0 _S3802;
    (&_S3802)->primal_0 = _S3648;
    (&_S3802)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3801, &_S3802, 0.0f);
    DiffPair_float_0 _S3803;
    (&_S3803)->primal_0 = _S3596;
    (&_S3803)->differential_0 = 0.0f;
    DiffPair_float_0 _S3804;
    (&_S3804)->primal_0 = _S3622;
    (&_S3804)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3803, &_S3804, _S3801.differential_0);
    DiffPair_float_0 _S3805;
    (&_S3805)->primal_0 = _S3659;
    (&_S3805)->differential_0 = 0.0f;
    DiffPair_float_0 _S3806;
    (&_S3806)->primal_0 = _S3648;
    (&_S3806)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3805, &_S3806, 0.0f);
    DiffPair_float_0 _S3807;
    (&_S3807)->primal_0 = _S3596;
    (&_S3807)->differential_0 = 0.0f;
    DiffPair_float_0 _S3808;
    (&_S3808)->primal_0 = _S3622;
    (&_S3808)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3807, &_S3808, _S3805.differential_0);
    DiffPair_float_0 _S3809;
    (&_S3809)->primal_0 = _S3657;
    (&_S3809)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3809, 0.0f);
    float _S3810 = - (-1.0f * - (_S3809.differential_0 / _S3658));
    float2  _S3811 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3812;
    (&_S3812)->primal_0 = e2_1;
    (&_S3812)->differential_0 = _S3811;
    s_bwd_length_impl_1(&_S3812, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3813;
    (&_S3813)->primal_0 = e1_5;
    (&_S3813)->differential_0 = _S3811;
    s_bwd_length_impl_1(&_S3813, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3814;
    (&_S3814)->primal_0 = e0_5;
    (&_S3814)->differential_0 = _S3811;
    s_bwd_length_impl_1(&_S3814, -0.0f);
    DiffPair_float_0 _S3815;
    (&_S3815)->primal_0 = _S3655;
    (&_S3815)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3815, 0.0f);
    float _S3816 = - _S3815.differential_0;
    float2  _S3817 = _S3813.differential_0 + make_float2 (_S3653 * _S3816, _S3651 * _S3815.differential_0);
    float2  _S3818 = _S3814.differential_0 + make_float2 (_S3652 * _S3815.differential_0, _S3654 * _S3816);
    float2  _S3819 = v_uv2_1 + - _S3812.differential_0 + _S3817;
    float _S3820 = fx_29 * (_S3802.differential_0 + _S3806.differential_0 + _S3819.x);
    float2  _S3821 = make_float2 (_S3820, fy_29 * (_S3794.differential_0 + _S3798.differential_0 + _S3819.y)) + make_float2 ((*dist_coeffs_35)[int(8)] * _S3820, (*dist_coeffs_35)[int(9)] * _S3820);
    float2  _S3822 = _S3637 * _S3821;
    float2  _S3823 = _S3641 * _S3821;
    float _S3824 = (*dist_coeffs_35)[int(4)] * _S3821.y;
    float _S3825 = (*dist_coeffs_35)[int(5)] * _S3821.x;
    float _S3826 = _S3822.x + _S3822.y;
    float _S3827 = r2_61 * _S3826;
    float _S3828 = r2_61 * _S3827;
    float _S3829 = (*dist_coeffs_35)[int(7)] * _S3821.y + _S3824 + (*dist_coeffs_35)[int(6)] * _S3821.x + _S3825 + _S3640 * _S3826 + _S3639 * _S3827 + _S3638 * _S3828 + (*dist_coeffs_35)[int(3)] * (r2_61 * _S3828);
    float _S3830 = v_61 * _S3829;
    float _S3831 = u_61 * _S3829;
    float _S3832 = _S3645 * _S3824 + 2.0f * (v_61 * _S3824) + _S3644 * _S3821.y + _S3642 * _S3821.x + _S3830 + _S3830;
    float _S3833 = _S3591 * (v_61 * _S3821.y) + _S3643 * _S3825 + 2.0f * (u_61 * _S3825) + _S3588 * (v_61 * _S3821.x) + _S3831 + _S3831;
    float2  _S3834 = v_uv0_1 + _S3812.differential_0 + - _S3818;
    float2  _S3835 = v_out_hardness_1 + make_float2 (0.0f, _S3810);
    float _S3836 = _S3804.differential_0 + _S3808.differential_0;
    float2  _S3837 = v_uv1_1 + - _S3817 + _S3818;
    float3  _S3838 = _S3713 + _S3714;
    FixedArray<float3 , 2>  _S3839;
    _S3839[int(0)] = _S3703;
    _S3839[int(1)] = _S3703;
    _S3839[int(1)] = _S3720;
    _S3839[int(0)] = _S3725;
    float3  _S3840 = _S3839[int(0)];
    float3  _S3841 = _S3839[int(1)];
    FixedArray<float3 , 16>  _S3842;
    _S3842[int(0)] = _S3703;
    _S3842[int(1)] = _S3703;
    _S3842[int(2)] = _S3703;
    _S3842[int(3)] = _S3703;
    _S3842[int(4)] = _S3703;
    _S3842[int(5)] = _S3703;
    _S3842[int(6)] = _S3703;
    _S3842[int(7)] = _S3703;
    _S3842[int(8)] = _S3703;
    _S3842[int(9)] = _S3703;
    _S3842[int(10)] = _S3703;
    _S3842[int(11)] = _S3703;
    _S3842[int(12)] = _S3703;
    _S3842[int(13)] = _S3703;
    _S3842[int(14)] = _S3703;
    _S3842[int(15)] = _S3703;
    _S3842[int(7)] = _S3749;
    _S3842[int(0)] = _S3783;
    _S3842[int(1)] = _S3770;
    _S3842[int(2)] = _S3768;
    _S3842[int(3)] = _S3766;
    _S3842[int(4)] = _S3755;
    _S3842[int(5)] = _S3753;
    _S3842[int(6)] = _S3751;
    _S3842[int(15)] = _S3727;
    _S3842[int(8)] = _S3747;
    _S3842[int(9)] = _S3739;
    _S3842[int(10)] = _S3737;
    _S3842[int(11)] = _S3735;
    _S3842[int(12)] = _S3733;
    _S3842[int(13)] = _S3731;
    _S3842[int(14)] = _S3729;
    float3  _S3843 = _S3842[int(0)];
    float3  _S3844 = _S3842[int(1)];
    float3  _S3845 = _S3842[int(2)];
    float3  _S3846 = _S3842[int(3)];
    float3  _S3847 = _S3842[int(4)];
    float3  _S3848 = _S3842[int(5)];
    float3  _S3849 = _S3842[int(6)];
    float3  _S3850 = _S3842[int(7)];
    float3  _S3851 = _S3842[int(8)];
    float3  _S3852 = _S3842[int(9)];
    float3  _S3853 = _S3842[int(10)];
    float3  _S3854 = _S3842[int(11)];
    float3  _S3855 = _S3842[int(12)];
    float3  _S3856 = _S3842[int(13)];
    float3  _S3857 = _S3842[int(14)];
    float3  _S3858 = _S3842[int(15)];
    float _S3859 = _S3803.differential_0 + _S3807.differential_0;
    float _S3860 = _S3795.differential_0 + _S3799.differential_0;
    float _S3861 = _S3796.differential_0 + _S3800.differential_0;
    float2  _S3862 = _S3823 + make_float2 (_S3833, _S3832);
    float2  _S3863 = _S3625 * _S3862;
    float2  _S3864 = _S3636 * _S3862;
    float _S3865 = _S3863.x + _S3863.y;
    if(_S3629)
    {
        float _S3866 = _S3865 / _S3631;
        float _S3867 = _S3632 * - _S3866;
        float _S3868 = _S3628 * (0.3333333432674408f * - (_S3627 * _S3866));
        k_9 = _S3868 + _S3868;
        _S3630 = _S3867;
        _S3631 = 0.0f;
    }
    else
    {
        float _S3869 = _S3865 / _S3630;
        float _S3870 = _S3628 * - _S3869;
        k_9 = _S3626 * _S3869;
        _S3630 = 0.0f;
        _S3631 = _S3870;
    }
    DiffPair_float_0 _S3871;
    (&_S3871)->primal_0 = _S3626;
    (&_S3871)->differential_0 = 0.0f;
    DiffPair_float_0 _S3872;
    (&_S3872)->primal_0 = _S3627;
    (&_S3872)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3871, &_S3872, k_9);
    float _S3873 = _S3872.differential_0 + _S3630;
    float _S3874 = _S3871.differential_0 + _S3631;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3875;
    (&_S3875)->primal_0 = _S3625;
    (&_S3875)->differential_0 = _S3811;
    s_bwd_length_impl_1(&_S3875, _S3874);
    float2  _S3876 = _S3875.differential_0 + _S3864;
    float _S3877 = fx_29 * (_S3837.x + _S3836);
    float2  _S3878 = make_float2 (_S3877, fy_29 * (_S3837.y + _S3861)) + make_float2 ((*dist_coeffs_35)[int(8)] * _S3877, (*dist_coeffs_35)[int(9)] * _S3877);
    float2  _S3879 = _S3611 * _S3878;
    float _S3880 = (*dist_coeffs_35)[int(4)] * _S3878.y;
    float _S3881 = (*dist_coeffs_35)[int(5)] * _S3878.x;
    float _S3882 = _S3879.x + _S3879.y;
    float _S3883 = r2_60 * _S3882;
    float _S3884 = r2_60 * _S3883;
    float _S3885 = (*dist_coeffs_35)[int(7)] * _S3878.y + _S3880 + (*dist_coeffs_35)[int(6)] * _S3878.x + _S3881 + _S3614 * _S3882 + _S3613 * _S3883 + _S3612 * _S3884 + (*dist_coeffs_35)[int(3)] * (r2_60 * _S3884);
    float _S3886 = v_60 * _S3885;
    float _S3887 = u_60 * _S3885;
    float3  _S3888 = _S3712.differential_0 + make_float3 (_S3876.x, _S3876.y, _S3873);
    float2  _S3889 = _S3615 * _S3878 + make_float2 (_S3591 * (v_60 * _S3878.y) + _S3617 * _S3881 + 2.0f * (u_60 * _S3881) + _S3588 * (v_60 * _S3878.x) + _S3887 + _S3887, _S3619 * _S3880 + 2.0f * (v_60 * _S3880) + _S3618 * _S3878.y + _S3616 * _S3878.x + _S3886 + _S3886);
    float2  _S3890 = _S3599 * _S3889;
    float2  _S3891 = _S3610 * _S3889;
    float _S3892 = _S3890.x + _S3890.y;
    if(_S3603)
    {
        float _S3893 = _S3892 / _S3605;
        float _S3894 = _S3606 * - _S3893;
        float _S3895 = _S3602 * (0.3333333432674408f * - (_S3601 * _S3893));
        k_9 = _S3895 + _S3895;
        _S3604 = _S3894;
        _S3605 = 0.0f;
    }
    else
    {
        float _S3896 = _S3892 / _S3604;
        float _S3897 = _S3602 * - _S3896;
        k_9 = _S3600 * _S3896;
        _S3604 = 0.0f;
        _S3605 = _S3897;
    }
    DiffPair_float_0 _S3898;
    (&_S3898)->primal_0 = _S3600;
    (&_S3898)->differential_0 = 0.0f;
    DiffPair_float_0 _S3899;
    (&_S3899)->primal_0 = _S3601;
    (&_S3899)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3898, &_S3899, k_9);
    float _S3900 = _S3899.differential_0 + _S3604;
    float _S3901 = _S3898.differential_0 + _S3605;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3902;
    (&_S3902)->primal_0 = _S3599;
    (&_S3902)->differential_0 = _S3811;
    s_bwd_length_impl_1(&_S3902, _S3901);
    float2  _S3903 = _S3902.differential_0 + _S3891;
    float _S3904 = fx_29 * (_S3834.x + _S3859);
    float2  _S3905 = make_float2 (_S3904, fy_29 * (_S3834.y + _S3860)) + make_float2 ((*dist_coeffs_35)[int(8)] * _S3904, (*dist_coeffs_35)[int(9)] * _S3904);
    float2  _S3906 = _S3583 * _S3905;
    float _S3907 = (*dist_coeffs_35)[int(4)] * _S3905.y;
    float _S3908 = (*dist_coeffs_35)[int(5)] * _S3905.x;
    float _S3909 = _S3906.x + _S3906.y;
    float _S3910 = r2_59 * _S3909;
    float _S3911 = r2_59 * _S3910;
    float _S3912 = (*dist_coeffs_35)[int(7)] * _S3905.y + _S3907 + (*dist_coeffs_35)[int(6)] * _S3905.x + _S3908 + _S3586 * _S3909 + _S3585 * _S3910 + _S3584 * _S3911 + (*dist_coeffs_35)[int(3)] * (r2_59 * _S3911);
    float _S3913 = v_59 * _S3912;
    float _S3914 = u_59 * _S3912;
    float3  _S3915 = _S3711.differential_0 + make_float3 (_S3903.x, _S3903.y, _S3900);
    float2  _S3916 = _S3587 * _S3905 + make_float2 (_S3591 * (v_59 * _S3905.y) + _S3590 * _S3908 + 2.0f * (u_59 * _S3908) + _S3588 * (v_59 * _S3905.x) + _S3914 + _S3914, _S3593 * _S3907 + 2.0f * (v_59 * _S3907) + _S3592 * _S3905.y + _S3589 * _S3905.x + _S3913 + _S3913);
    float2  _S3917 = _S3571 * _S3916;
    float2  _S3918 = _S3582 * _S3916;
    float _S3919 = _S3917.x + _S3917.y;
    if(_S3575)
    {
        float _S3920 = _S3919 / _S3577;
        float _S3921 = _S3578 * - _S3920;
        float _S3922 = _S3574 * (0.3333333432674408f * - (_S3573 * _S3920));
        k_9 = _S3922 + _S3922;
        _S3576 = _S3921;
        _S3577 = 0.0f;
    }
    else
    {
        float _S3923 = _S3919 / _S3576;
        float _S3924 = _S3574 * - _S3923;
        k_9 = _S3572 * _S3923;
        _S3576 = 0.0f;
        _S3577 = _S3924;
    }
    DiffPair_float_0 _S3925;
    (&_S3925)->primal_0 = _S3572;
    (&_S3925)->differential_0 = 0.0f;
    DiffPair_float_0 _S3926;
    (&_S3926)->primal_0 = _S3573;
    (&_S3926)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3925, &_S3926, k_9);
    float _S3927 = _S3926.differential_0 + _S3576;
    float _S3928 = _S3925.differential_0 + _S3577;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3929;
    (&_S3929)->primal_0 = _S3571;
    (&_S3929)->differential_0 = _S3811;
    s_bwd_length_impl_1(&_S3929, _S3928);
    float2  _S3930 = _S3929.differential_0 + _S3918;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3931;
    (&_S3931)->primal_0 = vert2_c_5;
    (&_S3931)->differential_0 = _S3703;
    s_bwd_length_impl_0(&_S3931, _S3792.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3932;
    (&_S3932)->primal_0 = vert1_c_5;
    (&_S3932)->differential_0 = _S3703;
    s_bwd_length_impl_0(&_S3932, _S3792.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3933;
    (&_S3933)->primal_0 = vert0_c_5;
    (&_S3933)->differential_0 = _S3703;
    s_bwd_length_impl_0(&_S3933, _S3792.differential_0.x);
    float3  _S3934 = _S3931.differential_0 + _S3888;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3935;
    (&_S3935)->primal_0 = R_24;
    (&_S3935)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3936;
    (&_S3936)->primal_0 = vert2_1;
    (&_S3936)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3935, &_S3936, _S3934);
    float3  _S3937 = _S3932.differential_0 + _S3915;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3938;
    (&_S3938)->primal_0 = R_24;
    (&_S3938)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3939;
    (&_S3939)->primal_0 = vert1_1;
    (&_S3939)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3938, &_S3939, _S3937);
    float3  _S3940 = _S3933.differential_0 + _S3838 + make_float3 (_S3930.x, _S3930.y, _S3927);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3941;
    (&_S3941)->primal_0 = R_24;
    (&_S3941)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3942;
    (&_S3942)->primal_0 = vert0_1;
    (&_S3942)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3941, &_S3942, _S3940);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3943;
    (&_S3943)->primal_0 = _S3562;
    (&_S3943)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3944;
    (&_S3944)->primal_0 = _S3567;
    (&_S3944)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3943, &_S3944, _S3936.differential_0);
    float _S3945 = - _S3944.differential_0.y;
    float _S3946 = _S3566 * _S3944.differential_0.x;
    float _S3947 = - (_S3558 * _S3944.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3948;
    (&_S3948)->primal_0 = _S3562;
    (&_S3948)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3949;
    (&_S3949)->primal_0 = _S3565;
    (&_S3949)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3948, &_S3949, _S3939.differential_0);
    float _S3950 = _S3558 * _S3949.differential_0.x;
    float _S3951 = _S3564 * _S3949.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3952;
    (&_S3952)->primal_0 = _S3562;
    (&_S3952)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3953;
    (&_S3953)->primal_0 = _S3563;
    (&_S3953)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3952, &_S3953, _S3942.differential_0);
    Matrix<float, 3, 3>  _S3954 = transpose_0(_S3943.differential_0 + _S3948.differential_0 + _S3952.differential_0);
    float _S3955 = 2.0f * - _S3954.rows[int(2)].z;
    float _S3956 = 2.0f * _S3954.rows[int(2)].y;
    float _S3957 = 2.0f * _S3954.rows[int(2)].x;
    float _S3958 = 2.0f * _S3954.rows[int(1)].z;
    float _S3959 = 2.0f * - _S3954.rows[int(1)].y;
    float _S3960 = 2.0f * _S3954.rows[int(1)].x;
    float _S3961 = 2.0f * _S3954.rows[int(0)].z;
    float _S3962 = 2.0f * _S3954.rows[int(0)].y;
    float _S3963 = 2.0f * - _S3954.rows[int(0)].x;
    float _S3964 = - _S3960 + _S3962;
    float _S3965 = _S3957 + - _S3961;
    float _S3966 = - _S3956 + _S3958;
    float _S3967 = _S3956 + _S3958;
    float _S3968 = _S3957 + _S3961;
    float _S3969 = _S3960 + _S3962;
    float _S3970 = quat_25.w * (_S3959 + _S3963);
    float _S3971 = quat_25.z * (_S3955 + _S3963);
    float _S3972 = quat_25.y * (_S3955 + _S3959);
    float _S3973 = quat_25.x * _S3964 + quat_25.z * _S3967 + quat_25.y * _S3968 + _S3970 + _S3970;
    float _S3974 = quat_25.x * _S3965 + quat_25.w * _S3967 + quat_25.y * _S3969 + _S3971 + _S3971;
    float _S3975 = quat_25.x * _S3966 + quat_25.w * _S3968 + quat_25.z * _S3969 + _S3972 + _S3972;
    float _S3976 = quat_25.w * _S3964 + quat_25.z * _S3965 + quat_25.y * _S3966;
    float _S3977 = _S3947 + _S3950;
    float _S3978 = 0.5f * - _S3977;
    float _S3979 = _S3945 + _S3949.differential_0.y;
    DiffPair_float_0 _S3980;
    (&_S3980)->primal_0 = _S3559;
    (&_S3980)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3980, _S3979);
    float _S3981 = _S3978 + _S3980.differential_0;
    float _S3982 = _S3946 + _S3951 + _S3953.differential_0.x;
    DiffPair_float_0 _S3983;
    (&_S3983)->primal_0 = _S3557;
    (&_S3983)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3983, _S3982);
    float _S3984 = _S3978 + _S3983.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3985;
    (&_S3985)->primal_0 = mean_c_20;
    (&_S3985)->differential_0 = _S3703;
    s_bwd_length_impl_0(&_S3985, 0.0f);
    float3  _S3986 = _S3985.differential_0 + _S3706.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3987;
    (&_S3987)->primal_0 = R_24;
    (&_S3987)->differential_0 = _S3786;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3988;
    (&_S3988)->primal_0 = mean_24;
    (&_S3988)->differential_0 = _S3703;
    s_bwd_prop_mul_1(&_S3987, &_S3988, _S3986);
    float3  _S3989 = _S3934 + _S3937 + _S3940 + _S3986 + _S3789.differential_0;
    Matrix<float, 3, 3>  _S3990 = _S3935.differential_0 + _S3938.differential_0 + _S3941.differential_0 + _S3987.differential_0 + _S3790;
    float3  _S3991 = make_float3 (_S3984, _S3981, _S3977);
    float4  _S3992 = make_float4 (0.0f);
    *&((&_S3992)->w) = _S3973;
    *&((&_S3992)->z) = _S3974;
    *&((&_S3992)->y) = _S3975;
    *&((&_S3992)->x) = _S3976;
    float4  _S3993 = _S3992;
    float3  _S3994 = _S3936.differential_0 + _S3939.differential_0 + _S3942.differential_0 + _S3988.differential_0 + _S3784;
    *v_mean_8 = _S3994;
    *v_quat_7 = _S3993;
    *v_scale_7 = _S3991;
    *v_hardness_1 = _S3835;
    (*v_sh_coeffs_6)[int(0)] = _S3843;
    (*v_sh_coeffs_6)[int(1)] = _S3844;
    (*v_sh_coeffs_6)[int(2)] = _S3845;
    (*v_sh_coeffs_6)[int(3)] = _S3846;
    (*v_sh_coeffs_6)[int(4)] = _S3847;
    (*v_sh_coeffs_6)[int(5)] = _S3848;
    (*v_sh_coeffs_6)[int(6)] = _S3849;
    (*v_sh_coeffs_6)[int(7)] = _S3850;
    (*v_sh_coeffs_6)[int(8)] = _S3851;
    (*v_sh_coeffs_6)[int(9)] = _S3852;
    (*v_sh_coeffs_6)[int(10)] = _S3853;
    (*v_sh_coeffs_6)[int(11)] = _S3854;
    (*v_sh_coeffs_6)[int(12)] = _S3855;
    (*v_sh_coeffs_6)[int(13)] = _S3856;
    (*v_sh_coeffs_6)[int(14)] = _S3857;
    (*v_sh_coeffs_6)[int(15)] = _S3858;
    (*v_ch_coeffs_1)[int(0)] = _S3840;
    (*v_ch_coeffs_1)[int(1)] = _S3841;
    *v_R_7 = _S3990;
    *v_t_7 = _S3989;
    return;
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
        DiffPair_float_0 _S3995 = *dpx_17;
        float _S3996 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S3996;
        float _S3997 = val_0 * (F32_log((_S3995.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3997;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_1)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3998 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3998))));
    float2  _S3999 = p_1 - v0_0;
    float2  _S4000 = normalize_1(e0_6);
    float2  _S4001 = p_1 - v1_0;
    float2  _S4002 = normalize_1(e1_6);
    float2  _S4003 = p_1 - v2_0;
    float2  _S4004 = normalize_1(e2_2);
    float _S4005 = hardness_6.x;
    float _S4006 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3999.x * _S4000.y - _S3999.y * _S4000.x)), (se_0 * (_S4001.x * _S4002.y - _S4001.y * _S4002.x))))), (se_0 * (_S4003.x * _S4004.y - _S4003.y * _S4004.x)))) / ((F32_abs((_S3998))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S4006))));
    float _S4007;
    if(a_1 <= 0.0f)
    {
        _S4007 = 0.0f;
    }
    else
    {
        _S4007 = (F32_min(((F32_pow((a_1), (_S4006)))), (0.99900001287460327f)));
    }
    return _S4005 * _S4007;
}

inline __device__ float s_primal_ctx_abs_0(float _S4008)
{
    return (F32_abs((_S4008)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S4009, float _S4010, float _S4011)
{
    return clamp_0(_S4009, _S4010, _S4011);
}

inline __device__ float s_primal_ctx_exp2_0(float _S4012)
{
    return (F32_exp2((_S4012)));
}

inline __device__ float s_primal_ctx_pow_0(float _S4013, float _S4014)
{
    return (F32_pow((_S4013), (_S4014)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S4015, DiffPair_float_0 * _S4016, float _S4017)
{
    _d_pow_0(_S4015, _S4016, _S4017);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S4018, DiffPair_float_0 * _S4019, DiffPair_float_0 * _S4020, float _S4021)
{
    _d_clamp_0(_S4018, _S4019, _S4020, _S4021);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_8)
{
    float _S4022 = length_0((*dpx_18).primal_0);
    float2  _S4023 = (*dpx_18).primal_0 * _s_dOut_8;
    float2  _S4024 = make_float2 (1.0f / _S4022) * _s_dOut_8;
    float _S4025 = - ((_S4023.x + _S4023.y) / (_S4022 * _S4022));
    float2  _S4026 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4027;
    (&_S4027)->primal_0 = (*dpx_18).primal_0;
    (&_S4027)->differential_0 = _S4026;
    s_bwd_length_impl_1(&_S4027, _S4025);
    float2  _S4028 = _S4024 + _S4027.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S4028;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4029, float2  _S4030)
{
    s_bwd_prop_normalize_impl_1(_S4029, _S4030);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_2, float _s_dOut_9)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S4031 = e0_7.x;
    float _S4032 = e1_7.y;
    float _S4033 = e0_7.y;
    float _S4034 = e1_7.x;
    float _S4035 = _S4031 * _S4032 - _S4033 * _S4034;
    float se_1 = float((F32_sign((_S4035))));
    float2  _S4036 = p_2 - (*dpv0_0).primal_0;
    float2  _S4037 = normalize_1(e0_7);
    float _S4038 = _S4036.x;
    float _S4039 = _S4037.y;
    float _S4040 = _S4036.y;
    float _S4041 = _S4037.x;
    float de0_0 = se_1 * (_S4038 * _S4039 - _S4040 * _S4041);
    float2  _S4042 = p_2 - (*dpv1_0).primal_0;
    float2  _S4043 = normalize_1(e1_7);
    float _S4044 = _S4042.x;
    float _S4045 = _S4043.y;
    float _S4046 = _S4042.y;
    float _S4047 = _S4043.x;
    float de1_0 = se_1 * (_S4044 * _S4045 - _S4046 * _S4047);
    float2  _S4048 = p_2 - (*dpv2_0).primal_0;
    float2  _S4049 = normalize_1(e2_3);
    float _S4050 = _S4048.x;
    float _S4051 = _S4049.y;
    float _S4052 = _S4048.y;
    float _S4053 = _S4049.x;
    float de2_0 = se_1 * (_S4050 * _S4051 - _S4052 * _S4053);
    float _S4054 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S4055 = s_primal_ctx_max_0(_S4054, de2_0);
    float _S4056 = s_primal_ctx_abs_0(_S4035);
    float _S4057 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S4056 / _S4057;
    float _S4058 = _S4057 * _S4057;
    float _S4059 = (*dphardness_0).primal_0.x;
    float _S4060 = (*dphardness_0).primal_0.y;
    float _S4061 = dmax_0 * dmax_0;
    float _S4062 = 1.0f + _S4055 / dmax_0;
    float _S4063 = 1.0f - s_primal_ctx_clamp_0(_S4060, 0.00499999988824129f, 0.98000001907348633f);
    float _S4064 = -1.0f / _S4063;
    float _S4065 = _S4063 * _S4063;
    float _S4066 = 1.0f - s_primal_ctx_exp2_0(_S4064);
    float a_2 = 1.0f - _S4062 * _S4066;
    bool _S4067 = a_2 <= 0.0f;
    float _S4068;
    float _S4069;
    if(_S4067)
    {
        _S4068 = 0.0f;
        _S4069 = 0.0f;
    }
    else
    {
        float _S4070 = s_primal_ctx_pow_0(a_2, _S4063);
        _S4068 = s_primal_ctx_min_0(_S4070, 0.99900001287460327f);
        _S4069 = _S4070;
    }
    float _S4071 = _S4059 * _s_dOut_9;
    float _S4072 = _S4068 * _s_dOut_9;
    if(_S4067)
    {
        _S4068 = 0.0f;
        _S4069 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4073;
        (&_S4073)->primal_0 = _S4069;
        (&_S4073)->differential_0 = 0.0f;
        DiffPair_float_0 _S4074;
        (&_S4074)->primal_0 = 0.99900001287460327f;
        (&_S4074)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4073, &_S4074, _S4071);
        DiffPair_float_0 _S4075;
        (&_S4075)->primal_0 = a_2;
        (&_S4075)->differential_0 = 0.0f;
        DiffPair_float_0 _S4076;
        (&_S4076)->primal_0 = _S4063;
        (&_S4076)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4075, &_S4076, _S4073.differential_0);
        _S4068 = _S4075.differential_0;
        _S4069 = _S4076.differential_0;
    }
    float _S4077 = - _S4068;
    float _S4078 = _S4066 * _S4077;
    float _S4079 = - (_S4062 * _S4077);
    DiffPair_float_0 _S4080;
    (&_S4080)->primal_0 = _S4064;
    (&_S4080)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4080, _S4079);
    float _S4081 = - (-1.0f * - (_S4080.differential_0 / _S4065) + _S4069);
    float _S4082 = _S4078 / _S4061;
    float s_diff_dmax_T_0 = _S4055 * - _S4082;
    float _S4083 = dmax_0 * _S4082;
    DiffPair_float_0 _S4084;
    (&_S4084)->primal_0 = _S4060;
    (&_S4084)->differential_0 = 0.0f;
    DiffPair_float_0 _S4085;
    (&_S4085)->primal_0 = 0.00499999988824129f;
    (&_S4085)->differential_0 = 0.0f;
    DiffPair_float_0 _S4086;
    (&_S4086)->primal_0 = 0.98000001907348633f;
    (&_S4086)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4084, &_S4085, &_S4086, _S4081);
    float _S4087 = s_diff_dmax_T_0 / _S4058;
    float _S4088 = _S4056 * - _S4087;
    float _S4089 = _S4057 * _S4087;
    float2  _S4090 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4091;
    (&_S4091)->primal_0 = e2_3;
    (&_S4091)->differential_0 = _S4090;
    s_bwd_length_impl_1(&_S4091, _S4088);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4092;
    (&_S4092)->primal_0 = e1_7;
    (&_S4092)->differential_0 = _S4090;
    s_bwd_length_impl_1(&_S4092, _S4088);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4093;
    (&_S4093)->primal_0 = e0_7;
    (&_S4093)->differential_0 = _S4090;
    s_bwd_length_impl_1(&_S4093, _S4088);
    DiffPair_float_0 _S4094;
    (&_S4094)->primal_0 = _S4035;
    (&_S4094)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4094, _S4089);
    DiffPair_float_0 _S4095;
    (&_S4095)->primal_0 = _S4054;
    (&_S4095)->differential_0 = 0.0f;
    DiffPair_float_0 _S4096;
    (&_S4096)->primal_0 = de2_0;
    (&_S4096)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4095, &_S4096, _S4083);
    DiffPair_float_0 _S4097;
    (&_S4097)->primal_0 = de0_0;
    (&_S4097)->differential_0 = 0.0f;
    DiffPair_float_0 _S4098;
    (&_S4098)->primal_0 = de1_0;
    (&_S4098)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4097, &_S4098, _S4095.differential_0);
    float _S4099 = se_1 * _S4096.differential_0;
    float _S4100 = - _S4099;
    float _S4101 = _S4053 * _S4100;
    float _S4102 = _S4051 * _S4099;
    float2  _S4103 = make_float2 (_S4052 * _S4100, _S4050 * _S4099);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4104;
    (&_S4104)->primal_0 = e2_3;
    (&_S4104)->differential_0 = _S4090;
    s_bwd_normalize_impl_1(&_S4104, _S4103);
    float2  _S4105 = - make_float2 (_S4102, _S4101);
    float _S4106 = se_1 * _S4098.differential_0;
    float _S4107 = - _S4106;
    float _S4108 = _S4047 * _S4107;
    float _S4109 = _S4045 * _S4106;
    float2  _S4110 = make_float2 (_S4046 * _S4107, _S4044 * _S4106);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4111;
    (&_S4111)->primal_0 = e1_7;
    (&_S4111)->differential_0 = _S4090;
    s_bwd_normalize_impl_1(&_S4111, _S4110);
    float2  _S4112 = - make_float2 (_S4109, _S4108);
    float _S4113 = se_1 * _S4097.differential_0;
    float _S4114 = - _S4113;
    float _S4115 = _S4041 * _S4114;
    float _S4116 = _S4039 * _S4113;
    float2  _S4117 = make_float2 (_S4040 * _S4114, _S4038 * _S4113);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4118;
    (&_S4118)->primal_0 = e0_7;
    (&_S4118)->differential_0 = _S4090;
    s_bwd_normalize_impl_1(&_S4118, _S4117);
    float2  _S4119 = - make_float2 (_S4116, _S4115);
    float _S4120 = - _S4094.differential_0;
    float2  _S4121 = _S4091.differential_0 + _S4104.differential_0;
    float2  _S4122 = - _S4121;
    float2  _S4123 = _S4092.differential_0 + _S4111.differential_0 + make_float2 (_S4033 * _S4120, _S4031 * _S4094.differential_0);
    float2  _S4124 = - _S4123;
    float2  _S4125 = _S4093.differential_0 + _S4118.differential_0 + make_float2 (_S4032 * _S4094.differential_0, _S4034 * _S4120);
    float2  _S4126 = - _S4125;
    float2  _S4127 = make_float2 (_S4072, _S4084.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S4127;
    float2  _S4128 = _S4105 + _S4122 + _S4123;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S4128;
    float2  _S4129 = _S4112 + _S4124 + _S4125;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S4129;
    float2  _S4130 = _S4119 + _S4121 + _S4126;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S4130;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4131, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4132, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4133, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4134, float2  _S4135, float _S4136)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S4131, _S4132, _S4133, _S4134, _S4135, _S4136);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_3, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S4137 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S4137;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S4137;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S4137;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S4137;
    s_bwd_evaluate_alpha_opaque_triangle_fast_0(&dp_v0_0, &dp_v1_0, &dp_v2_0, &dp_hardness_0, p_3, v_alpha_1);
    *v_v0_0 = dp_v0_0.differential_0;
    *v_v1_0 = dp_v2_0.differential_0;
    *v_v2_0 = dp_v1_0.differential_0;
    *v_hardness_2 = dp_hardness_0.differential_0;
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_precise(float2  v0_2, float2  v1_2, float2  v2_2, float2  hardness_8, float2  p_4)
{
    float2  e0_8 = v1_2 - v0_2;
    float2  e1_8 = v2_2 - v1_2;
    float2  e2_4 = v0_2 - v2_2;
    float2  _S4138 = p_4 - v0_2;
    float2  _S4139 = p_4 - v1_2;
    float2  _S4140 = p_4 - v2_2;
    float _S4141 = e0_8.x;
    float _S4142 = e1_8.y;
    float _S4143 = e0_8.y;
    float _S4144 = e1_8.x;
    float _S4145 = _S4141 * _S4142 - _S4143 * _S4144;
    float se_2 = float((F32_sign((_S4145))));
    float _S4146 = hardness_8.x;
    float _S4147 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S4138.x * _S4143 - _S4138.y * _S4141)), (se_2 * (_S4139.x * _S4142 - _S4139.y * _S4144))))), (se_2 * (_S4140.x * e2_4.y - _S4140.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S4138 - e0_8 * make_float2 (clamp_0(dot_1(_S4138, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S4139 - e1_8 * make_float2 (clamp_0(dot_1(_S4139, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S4140 - e2_4 * make_float2 (clamp_0(dot_1(_S4140, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S4145))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S4147))));
    float _S4148;
    if(a_3 <= 0.0f)
    {
        _S4148 = 0.0f;
    }
    else
    {
        _S4148 = (F32_min(((F32_pow((a_3), (_S4147)))), (0.99900001287460327f)));
    }
    return _S4146 * _S4148;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S4149, float2  _S4150)
{
    return dot_1(_S4149, _S4150);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4151, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4152, float _S4153)
{
    _d_dot_1(_S4151, _S4152, _S4153);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_5, float _s_dOut_10)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S4154 = p_5 - (*dpv0_1).primal_0;
    float _S4155 = s_primal_ctx_dot_1(_S4154, e0_9);
    float _S4156 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S4157 = _S4155 / _S4156;
    float _S4158 = _S4156 * _S4156;
    float _S4159 = s_primal_ctx_clamp_0(_S4157, 0.0f, 1.0f);
    float2  _S4160 = make_float2 (_S4159);
    float2  _S4161 = _S4154 - e0_9 * make_float2 (_S4159);
    float _S4162 = length_0(_S4161);
    float2  _S4163 = p_5 - (*dpv1_1).primal_0;
    float _S4164 = s_primal_ctx_dot_1(_S4163, e1_9);
    float _S4165 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S4166 = _S4164 / _S4165;
    float _S4167 = _S4165 * _S4165;
    float _S4168 = s_primal_ctx_clamp_0(_S4166, 0.0f, 1.0f);
    float2  _S4169 = make_float2 (_S4168);
    float2  _S4170 = _S4163 - e1_9 * make_float2 (_S4168);
    float _S4171 = length_0(_S4170);
    float2  _S4172 = p_5 - (*dpv2_1).primal_0;
    float _S4173 = s_primal_ctx_dot_1(_S4172, e2_5);
    float _S4174 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S4175 = _S4173 / _S4174;
    float _S4176 = _S4174 * _S4174;
    float _S4177 = s_primal_ctx_clamp_0(_S4175, 0.0f, 1.0f);
    float2  _S4178 = make_float2 (_S4177);
    float2  _S4179 = _S4172 - e2_5 * make_float2 (_S4177);
    float _S4180 = length_0(_S4179);
    float _S4181 = e0_9.x;
    float _S4182 = e1_9.y;
    float _S4183 = e0_9.y;
    float _S4184 = e1_9.x;
    float _S4185 = _S4181 * _S4182 - _S4183 * _S4184;
    float se_3 = float((F32_sign((_S4185))));
    float _S4186 = _S4154.x;
    float _S4187 = _S4154.y;
    float s0_0 = se_3 * (_S4186 * _S4183 - _S4187 * _S4181);
    float _S4188 = _S4163.x;
    float _S4189 = _S4163.y;
    float s1_0 = se_3 * (_S4188 * _S4182 - _S4189 * _S4184);
    float _S4190 = _S4172.x;
    float _S4191 = e2_5.y;
    float _S4192 = _S4172.y;
    float _S4193 = e2_5.x;
    float s2_0 = se_3 * (_S4190 * _S4191 - _S4192 * _S4193);
    float _S4194 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S4194, s2_0)))));
    float _S4195 = s_primal_ctx_min_0(_S4162, _S4171);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S4195, _S4180);
    float _S4196 = s_primal_ctx_abs_0(_S4185);
    float _S4197 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S4196 / _S4197;
    float _S4198 = _S4197 * _S4197;
    float _S4199 = (*dphardness_1).primal_0.x;
    float _S4200 = (*dphardness_1).primal_0.y;
    float _S4201 = dmax_1 * dmax_1;
    float _S4202 = 1.0f + dv_0 / dmax_1;
    float _S4203 = 1.0f - s_primal_ctx_clamp_0(_S4200, 0.00499999988824129f, 0.98000001907348633f);
    float _S4204 = -1.0f / _S4203;
    float _S4205 = _S4203 * _S4203;
    float _S4206 = 1.0f - s_primal_ctx_exp2_0(_S4204);
    float a_4 = 1.0f - _S4202 * _S4206;
    bool _S4207 = a_4 <= 0.0f;
    float _S4208;
    float _S4209;
    if(_S4207)
    {
        _S4208 = 0.0f;
        _S4209 = 0.0f;
    }
    else
    {
        float _S4210 = s_primal_ctx_pow_0(a_4, _S4203);
        _S4208 = s_primal_ctx_min_0(_S4210, 0.99900001287460327f);
        _S4209 = _S4210;
    }
    float _S4211 = _S4199 * _s_dOut_10;
    float _S4212 = _S4208 * _s_dOut_10;
    if(_S4207)
    {
        _S4208 = 0.0f;
        _S4209 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4213;
        (&_S4213)->primal_0 = _S4209;
        (&_S4213)->differential_0 = 0.0f;
        DiffPair_float_0 _S4214;
        (&_S4214)->primal_0 = 0.99900001287460327f;
        (&_S4214)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4213, &_S4214, _S4211);
        DiffPair_float_0 _S4215;
        (&_S4215)->primal_0 = a_4;
        (&_S4215)->differential_0 = 0.0f;
        DiffPair_float_0 _S4216;
        (&_S4216)->primal_0 = _S4203;
        (&_S4216)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4215, &_S4216, _S4213.differential_0);
        _S4208 = _S4215.differential_0;
        _S4209 = _S4216.differential_0;
    }
    float _S4217 = - _S4208;
    float _S4218 = _S4206 * _S4217;
    float _S4219 = - (_S4202 * _S4217);
    DiffPair_float_0 _S4220;
    (&_S4220)->primal_0 = _S4204;
    (&_S4220)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4220, _S4219);
    float _S4221 = - (-1.0f * - (_S4220.differential_0 / _S4205) + _S4209);
    float _S4222 = _S4218 / _S4201;
    float s_diff_dmax_T_1 = dv_0 * - _S4222;
    float s_diff_dv_T_0 = dmax_1 * _S4222;
    DiffPair_float_0 _S4223;
    (&_S4223)->primal_0 = _S4200;
    (&_S4223)->differential_0 = 0.0f;
    DiffPair_float_0 _S4224;
    (&_S4224)->primal_0 = 0.00499999988824129f;
    (&_S4224)->differential_0 = 0.0f;
    DiffPair_float_0 _S4225;
    (&_S4225)->primal_0 = 0.98000001907348633f;
    (&_S4225)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4223, &_S4224, &_S4225, _S4221);
    float _S4226 = s_diff_dmax_T_1 / _S4198;
    float _S4227 = _S4196 * - _S4226;
    float _S4228 = _S4197 * _S4226;
    float2  _S4229 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4230;
    (&_S4230)->primal_0 = e2_5;
    (&_S4230)->differential_0 = _S4229;
    s_bwd_length_impl_1(&_S4230, _S4227);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4231;
    (&_S4231)->primal_0 = e1_9;
    (&_S4231)->differential_0 = _S4229;
    s_bwd_length_impl_1(&_S4231, _S4227);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4232;
    (&_S4232)->primal_0 = e0_9;
    (&_S4232)->differential_0 = _S4229;
    s_bwd_length_impl_1(&_S4232, _S4227);
    DiffPair_float_0 _S4233;
    (&_S4233)->primal_0 = _S4185;
    (&_S4233)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4233, _S4228);
    float _S4234 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S4235;
    (&_S4235)->primal_0 = _S4195;
    (&_S4235)->differential_0 = 0.0f;
    DiffPair_float_0 _S4236;
    (&_S4236)->primal_0 = _S4180;
    (&_S4236)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4235, &_S4236, _S4234);
    DiffPair_float_0 _S4237;
    (&_S4237)->primal_0 = _S4162;
    (&_S4237)->differential_0 = 0.0f;
    DiffPair_float_0 _S4238;
    (&_S4238)->primal_0 = _S4171;
    (&_S4238)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4237, &_S4238, _S4235.differential_0);
    DiffPair_float_0 _S4239;
    (&_S4239)->primal_0 = _S4194;
    (&_S4239)->differential_0 = 0.0f;
    DiffPair_float_0 _S4240;
    (&_S4240)->primal_0 = s2_0;
    (&_S4240)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4239, &_S4240, 0.0f);
    DiffPair_float_0 _S4241;
    (&_S4241)->primal_0 = s0_0;
    (&_S4241)->differential_0 = 0.0f;
    DiffPair_float_0 _S4242;
    (&_S4242)->primal_0 = s1_0;
    (&_S4242)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4241, &_S4242, _S4239.differential_0);
    float _S4243 = se_3 * _S4240.differential_0;
    float _S4244 = - _S4243;
    float _S4245 = _S4192 * _S4244;
    float _S4246 = _S4193 * _S4244;
    float _S4247 = _S4190 * _S4243;
    float _S4248 = _S4191 * _S4243;
    float _S4249 = se_3 * _S4242.differential_0;
    float _S4250 = - _S4249;
    float _S4251 = _S4184 * _S4250;
    float _S4252 = _S4182 * _S4249;
    float _S4253 = se_3 * _S4241.differential_0;
    float _S4254 = - _S4253;
    float _S4255 = _S4181 * _S4254;
    float _S4256 = _S4183 * _S4253;
    float _S4257 = - _S4233.differential_0;
    float _S4258 = _S4189 * _S4250 + _S4183 * _S4257;
    float _S4259 = _S4186 * _S4253 + _S4184 * _S4257;
    float _S4260 = _S4188 * _S4249 + _S4181 * _S4233.differential_0;
    float _S4261 = _S4187 * _S4254 + _S4182 * _S4233.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4262;
    (&_S4262)->primal_0 = _S4179;
    (&_S4262)->differential_0 = _S4229;
    s_bwd_length_impl_1(&_S4262, _S4236.differential_0);
    float2  _S4263 = - _S4262.differential_0;
    float2  _S4264 = e2_5 * _S4263;
    float2  _S4265 = _S4178 * _S4263;
    float _S4266 = _S4264.x + _S4264.y;
    DiffPair_float_0 _S4267;
    (&_S4267)->primal_0 = _S4175;
    (&_S4267)->differential_0 = 0.0f;
    DiffPair_float_0 _S4268;
    (&_S4268)->primal_0 = 0.0f;
    (&_S4268)->differential_0 = 0.0f;
    DiffPair_float_0 _S4269;
    (&_S4269)->primal_0 = 1.0f;
    (&_S4269)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4267, &_S4268, &_S4269, _S4266);
    float _S4270 = _S4267.differential_0 / _S4176;
    float _S4271 = _S4173 * - _S4270;
    float _S4272 = _S4174 * _S4270;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4273;
    (&_S4273)->primal_0 = e2_5;
    (&_S4273)->differential_0 = _S4229;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4274;
    (&_S4274)->primal_0 = e2_5;
    (&_S4274)->differential_0 = _S4229;
    s_bwd_prop_dot_1(&_S4273, &_S4274, _S4271);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4275;
    (&_S4275)->primal_0 = _S4172;
    (&_S4275)->differential_0 = _S4229;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4276;
    (&_S4276)->primal_0 = e2_5;
    (&_S4276)->differential_0 = _S4229;
    s_bwd_prop_dot_1(&_S4275, &_S4276, _S4272);
    float2  _S4277 = - (_S4262.differential_0 + _S4275.differential_0 + make_float2 (_S4248, _S4246));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4278;
    (&_S4278)->primal_0 = _S4170;
    (&_S4278)->differential_0 = _S4229;
    s_bwd_length_impl_1(&_S4278, _S4238.differential_0);
    float2  _S4279 = - _S4278.differential_0;
    float2  _S4280 = e1_9 * _S4279;
    float2  _S4281 = _S4169 * _S4279;
    float _S4282 = _S4280.x + _S4280.y;
    DiffPair_float_0 _S4283;
    (&_S4283)->primal_0 = _S4166;
    (&_S4283)->differential_0 = 0.0f;
    DiffPair_float_0 _S4284;
    (&_S4284)->primal_0 = 0.0f;
    (&_S4284)->differential_0 = 0.0f;
    DiffPair_float_0 _S4285;
    (&_S4285)->primal_0 = 1.0f;
    (&_S4285)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4283, &_S4284, &_S4285, _S4282);
    float _S4286 = _S4283.differential_0 / _S4167;
    float _S4287 = _S4164 * - _S4286;
    float _S4288 = _S4165 * _S4286;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4289;
    (&_S4289)->primal_0 = e1_9;
    (&_S4289)->differential_0 = _S4229;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4290;
    (&_S4290)->primal_0 = e1_9;
    (&_S4290)->differential_0 = _S4229;
    s_bwd_prop_dot_1(&_S4289, &_S4290, _S4287);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4291;
    (&_S4291)->primal_0 = _S4163;
    (&_S4291)->differential_0 = _S4229;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4292;
    (&_S4292)->primal_0 = e1_9;
    (&_S4292)->differential_0 = _S4229;
    s_bwd_prop_dot_1(&_S4291, &_S4292, _S4288);
    float2  _S4293 = - (_S4278.differential_0 + _S4291.differential_0 + make_float2 (_S4252, _S4251));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4294;
    (&_S4294)->primal_0 = _S4161;
    (&_S4294)->differential_0 = _S4229;
    s_bwd_length_impl_1(&_S4294, _S4237.differential_0);
    float2  _S4295 = - _S4294.differential_0;
    float2  _S4296 = e0_9 * _S4295;
    float2  _S4297 = _S4160 * _S4295;
    float _S4298 = _S4296.x + _S4296.y;
    DiffPair_float_0 _S4299;
    (&_S4299)->primal_0 = _S4157;
    (&_S4299)->differential_0 = 0.0f;
    DiffPair_float_0 _S4300;
    (&_S4300)->primal_0 = 0.0f;
    (&_S4300)->differential_0 = 0.0f;
    DiffPair_float_0 _S4301;
    (&_S4301)->primal_0 = 1.0f;
    (&_S4301)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4299, &_S4300, &_S4301, _S4298);
    float _S4302 = _S4299.differential_0 / _S4158;
    float _S4303 = _S4155 * - _S4302;
    float _S4304 = _S4156 * _S4302;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4305;
    (&_S4305)->primal_0 = e0_9;
    (&_S4305)->differential_0 = _S4229;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4306;
    (&_S4306)->primal_0 = e0_9;
    (&_S4306)->differential_0 = _S4229;
    s_bwd_prop_dot_1(&_S4305, &_S4306, _S4303);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4307;
    (&_S4307)->primal_0 = _S4154;
    (&_S4307)->differential_0 = _S4229;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4308;
    (&_S4308)->primal_0 = e0_9;
    (&_S4308)->differential_0 = _S4229;
    s_bwd_prop_dot_1(&_S4307, &_S4308, _S4304);
    float2  _S4309 = - (_S4294.differential_0 + _S4307.differential_0 + make_float2 (_S4256, _S4255));
    float2  _S4310 = _S4230.differential_0 + _S4265 + _S4274.differential_0 + _S4273.differential_0 + _S4276.differential_0 + make_float2 (_S4245, _S4247);
    float2  _S4311 = - _S4310;
    float2  _S4312 = _S4231.differential_0 + _S4281 + _S4290.differential_0 + _S4289.differential_0 + _S4292.differential_0 + make_float2 (_S4258, _S4260);
    float2  _S4313 = - _S4312;
    float2  _S4314 = _S4232.differential_0 + _S4297 + _S4306.differential_0 + _S4305.differential_0 + _S4308.differential_0 + make_float2 (_S4261, _S4259);
    float2  _S4315 = - _S4314;
    float2  _S4316 = make_float2 (_S4212, _S4223.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S4316;
    float2  _S4317 = _S4277 + _S4311 + _S4312;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S4317;
    float2  _S4318 = _S4293 + _S4313 + _S4314;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S4318;
    float2  _S4319 = _S4309 + _S4310 + _S4315;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S4319;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4320, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4321, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4322, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4323, float2  _S4324, float _S4325)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S4320, _S4321, _S4322, _S4323, _S4324, _S4325);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_6, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S4326 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S4326;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S4326;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S4326;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S4326;
    s_bwd_evaluate_alpha_opaque_triangle_precise_0(&dp_v0_1, &dp_v1_1, &dp_v2_1, &dp_hardness_1, p_6, v_alpha_2);
    *v_v0_1 = dp_v0_1.differential_0;
    *v_v1_1 = dp_v2_1.differential_0;
    *v_v2_1 = dp_v1_1.differential_0;
    *v_hardness_3 = dp_hardness_1.differential_0;
    return;
}

inline __device__ void evaluate_color_opaque_triangle(float2  v0_4, float2  v1_4, float2  v2_4, FixedArray<float3 , 3>  * colors_0, float3  depths_0, float2  p_7, float3  * color_6, float * depth_15)
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

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_2, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_2, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpcolors_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpdepths_0, float2  p_8, float3  dpcolor_0, float dpdepth_1)
{
    float _S4327 = 0.3333333432674408f * dpdepth_1;
    float3  _S4328 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S4329 = make_float3 (0.0f);
    float3  _S4330 = _S4329;
    *&((&_S4330)->z) = _S4327;
    *&((&_S4330)->y) = _S4327;
    *&((&_S4330)->x) = _S4327;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S4330;
    FixedArray<float3 , 3>  _S4331;
    _S4331[int(0)] = _S4329;
    _S4331[int(1)] = _S4329;
    _S4331[int(2)] = _S4329;
    _S4331[int(2)] = _S4328;
    _S4331[int(1)] = _S4328;
    _S4331[int(0)] = _S4328;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S4331;
    float2  _S4332 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S4332;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S4332;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S4332;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4333, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4334, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4335, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4336, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4337, float2  _S4338, float3  _S4339, float _S4340)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S4333, _S4334, _S4335, _S4336, _S4337, _S4338, _S4339, _S4340);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_9, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S4341 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S4341;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S4341;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S4341;
    float3  _S4342 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4343 = { _S4342, _S4342, _S4342 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S4343;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S4342;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_9, v_color_0, v_depth_8);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_25, float4  quat_26, float3  scale_25, float3  * vert0_2, float3  * vert1_2, float3  * vert2_2)
{
    float _S4344 = scale_25.x;
    float sx_4 = (F32_exp((_S4344)));
    float _S4345 = scale_25.y;
    float sy_4 = (F32_exp((_S4345)));
    float sz_6 = scale_25.z - 0.5f * (_S4344 + _S4345);
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
    Matrix<float, 3, 3>  _S4346 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_47), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_47), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
    *vert0_2 = mul_0(_S4346, make_float3 (sx_4, 0.0f, 0.0f)) + mean_25;
    *vert1_2 = mul_0(_S4346, make_float3 (sx_4 * (-0.5f + sz_6), sy_4, 0.0f)) + mean_25;
    *vert2_2 = mul_0(_S4346, make_float3 (sx_4 * (-0.5f - sz_6), - sy_4, 0.0f)) + mean_25;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_10, FixedArray<float3 , 16>  * sh_coeffs_21, FixedArray<float3 , 2>  * ch_coeffs_6, Matrix<float, 3, 3>  R_25, float3  t_24, float fx_30, float fy_30, float cx_25, float cy_25, FixedArray<float, 10>  * dist_coeffs_36, uint image_width_21, uint image_height_21, float near_plane_14, float far_plane_14, int4  * aabb_xyxy_14, float * depth_16, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_4)
{
    for(;;)
    {
        float3  mean_c_21 = mul_0(R_25, mean_26) + t_24;
        float _S4347 = mean_c_21.z;
        bool _S4348;
        if(_S4347 < near_plane_14)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = _S4347 > far_plane_14;
        }
        if(_S4348)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4349 = scale_26.x;
        float sx_5 = (F32_exp((_S4349)));
        float _S4350 = scale_26.y;
        float sy_5 = (F32_exp((_S4350)));
        float sz_7 = scale_26.z - 0.5f * (_S4349 + _S4350);
        float x_53 = quat_27.y;
        float x2_27 = x_53 * x_53;
        float y2_27 = quat_27.z * quat_27.z;
        float z2_48 = quat_27.w * quat_27.w;
        float xy_27 = quat_27.y * quat_27.z;
        float xz_27 = quat_27.y * quat_27.w;
        float yz_27 = quat_27.z * quat_27.w;
        float wx_27 = quat_27.x * quat_27.y;
        float wy_27 = quat_27.x * quat_27.z;
        float wz_27 = quat_27.x * quat_27.w;
        Matrix<float, 3, 3>  _S4351 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_48), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_48), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S4351, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S4351, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S4351, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_6 = mul_0(R_25, vert0_3) + t_24;
        float3  vert1_c_6 = mul_0(R_25, vert1_3) + t_24;
        float3  vert2_c_6 = mul_0(R_25, vert2_3) + t_24;
        float _S4352 = vert0_c_6.z;
        float _S4353 = vert1_c_6.z;
        float _S4354 = vert2_c_6.z;
        if(_S4352 < near_plane_14)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = _S4352 > far_plane_14;
        }
        if(_S4348)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = _S4353 < near_plane_14;
        }
        if(_S4348)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = _S4353 > far_plane_14;
        }
        if(_S4348)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = _S4354 < near_plane_14;
        }
        if(_S4348)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = _S4354 > far_plane_14;
        }
        if(_S4348)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_4;
        for(;;)
        {
            float2  uv0_5 = float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S4352);
            if(_S4352 < 0.0f)
            {
                _S4348 = true;
            }
            else
            {
                bool _S4355 = is_valid_distortion(uv0_5, dist_coeffs_36);
                _S4348 = !_S4355;
            }
            if(_S4348)
            {
                uv0_4 = uv0_5;
                break;
            }
            float u_62 = uv0_5.x;
            float v_62 = uv0_5.y;
            float r2_62 = u_62 * u_62 + v_62 * v_62;
            float2  _S4356 = uv0_5 * make_float2 (1.0f + r2_62 * ((*dist_coeffs_36)[int(0)] + r2_62 * ((*dist_coeffs_36)[int(1)] + r2_62 * ((*dist_coeffs_36)[int(2)] + r2_62 * (*dist_coeffs_36)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_36)[int(4)] * u_62 * v_62 + (*dist_coeffs_36)[int(5)] * (r2_62 + 2.0f * u_62 * u_62) + (*dist_coeffs_36)[int(6)] * r2_62, 2.0f * (*dist_coeffs_36)[int(5)] * u_62 * v_62 + (*dist_coeffs_36)[int(4)] * (r2_62 + 2.0f * v_62 * v_62) + (*dist_coeffs_36)[int(7)] * r2_62);
            float2  _S4357 = _S4356 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4356.x + (*dist_coeffs_36)[int(9)] * _S4356.y, 0.0f);
            uv0_4 = make_float2 (fx_30 * _S4357.x + cx_25, fy_30 * _S4357.y + cy_25);
            break;
        }
        float2  uv1_4;
        bool all_valid_12 = true & (!_S4348);
        for(;;)
        {
            float2  uv1_5 = float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S4353);
            if(_S4353 < 0.0f)
            {
                _S4348 = true;
            }
            else
            {
                bool _S4358 = is_valid_distortion(uv1_5, dist_coeffs_36);
                _S4348 = !_S4358;
            }
            if(_S4348)
            {
                uv1_4 = uv1_5;
                break;
            }
            float u_63 = uv1_5.x;
            float v_63 = uv1_5.y;
            float r2_63 = u_63 * u_63 + v_63 * v_63;
            float2  _S4359 = uv1_5 * make_float2 (1.0f + r2_63 * ((*dist_coeffs_36)[int(0)] + r2_63 * ((*dist_coeffs_36)[int(1)] + r2_63 * ((*dist_coeffs_36)[int(2)] + r2_63 * (*dist_coeffs_36)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_36)[int(4)] * u_63 * v_63 + (*dist_coeffs_36)[int(5)] * (r2_63 + 2.0f * u_63 * u_63) + (*dist_coeffs_36)[int(6)] * r2_63, 2.0f * (*dist_coeffs_36)[int(5)] * u_63 * v_63 + (*dist_coeffs_36)[int(4)] * (r2_63 + 2.0f * v_63 * v_63) + (*dist_coeffs_36)[int(7)] * r2_63);
            float2  _S4360 = _S4359 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4359.x + (*dist_coeffs_36)[int(9)] * _S4359.y, 0.0f);
            uv1_4 = make_float2 (fx_30 * _S4360.x + cx_25, fy_30 * _S4360.y + cy_25);
            break;
        }
        float2  uv2_4;
        bool all_valid_13 = all_valid_12 & (!_S4348);
        for(;;)
        {
            float2  uv2_5 = float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S4354);
            if(_S4354 < 0.0f)
            {
                _S4348 = true;
            }
            else
            {
                bool _S4361 = is_valid_distortion(uv2_5, dist_coeffs_36);
                _S4348 = !_S4361;
            }
            if(_S4348)
            {
                uv2_4 = uv2_5;
                break;
            }
            float u_64 = uv2_5.x;
            float v_64 = uv2_5.y;
            float r2_64 = u_64 * u_64 + v_64 * v_64;
            float2  _S4362 = uv2_5 * make_float2 (1.0f + r2_64 * ((*dist_coeffs_36)[int(0)] + r2_64 * ((*dist_coeffs_36)[int(1)] + r2_64 * ((*dist_coeffs_36)[int(2)] + r2_64 * (*dist_coeffs_36)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_36)[int(4)] * u_64 * v_64 + (*dist_coeffs_36)[int(5)] * (r2_64 + 2.0f * u_64 * u_64) + (*dist_coeffs_36)[int(6)] * r2_64, 2.0f * (*dist_coeffs_36)[int(5)] * u_64 * v_64 + (*dist_coeffs_36)[int(4)] * (r2_64 + 2.0f * v_64 * v_64) + (*dist_coeffs_36)[int(7)] * r2_64);
            float2  _S4363 = _S4362 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4362.x + (*dist_coeffs_36)[int(9)] * _S4362.y, 0.0f);
            uv2_4 = make_float2 (fx_30 * _S4363.x + cx_25, fy_30 * _S4363.y + cy_25);
            break;
        }
        if(!(all_valid_13 & (!_S4348)))
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_10 = uv1_4 - uv0_4;
        float2  e1_10 = uv2_4 - uv1_4;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(uv0_4 - uv2_4)));
        float _S4364 = uv0_4.x;
        float _S4365 = uv1_4.x;
        float _S4366 = uv2_4.x;
        float xmax_7 = (F32_max(((F32_max((_S4364), (_S4365)))), (_S4366))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S4364), (_S4365)))), (_S4366))) - offset_4;
        float _S4367 = uv0_4.y;
        float _S4368 = uv1_4.y;
        float _S4369 = uv2_4.y;
        float ymax_7 = (F32_max(((F32_max((_S4367), (_S4368)))), (_S4369))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S4367), (_S4368)))), (_S4369))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = xmin_7 >= float(image_width_21);
        }
        if(_S4348)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = ymax_7 <= 0.0f;
        }
        if(_S4348)
        {
            _S4348 = true;
        }
        else
        {
            _S4348 = ymin_7 >= float(image_height_21);
        }
        if(_S4348)
        {
            _S4348 = true;
        }
        else
        {
            if(_S4347 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S4348 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S4348 = false;
                }
                if(_S4348)
                {
                    _S4348 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S4348 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S4348 = false;
                    }
                }
            }
            else
            {
                _S4348 = false;
            }
        }
        if(_S4348)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4370 = mean_26 - - mul_0(transpose_0(R_25), t_24);
        float _S4371 = _S4370.x;
        float _S4372 = _S4370.y;
        float _S4373 = _S4370.z;
        float norm_14 = (F32_sqrt((_S4371 * _S4371 + _S4372 * _S4372 + _S4373 * _S4373)));
        float x_54 = _S4371 / norm_14;
        float y_24 = _S4372 / norm_14;
        float z_21 = _S4373 / norm_14;
        float z2_49 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_54 * x_54 - y_24 * y_24;
        float fS1_21 = 2.0f * x_54 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_49 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_54) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_49 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_54) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_54 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_49 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_54) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_54 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S4374 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S4374);
        float3  _S4375 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S4376 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S4375 + _S4376 + make_float3 (0.5f), _S4374);
        (*rgbs_0)[int(2)] = max_0(_S4375 - _S4376 + make_float3 (0.5f), _S4374);
        (*verts_0)[int(0)] = vert0_3;
        (*verts_0)[int(1)] = vert1_3;
        (*verts_0)[int(2)] = vert2_3;
        float3  _S4377 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S4377 * make_float3 (float(- (F32_sign((dot_0(_S4377, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_26, float3  t_25, float fx_31, float fy_31, float cx_26, float cy_26, FixedArray<float, 10>  * dist_coeffs_37, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    bool _S4378;
    bool _S4379;
    bool _S4380;
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_26, mean_27) + t_25;
        float _S4381 = length_1(mean_c_22);
        bool _S4382;
        if(_S4381 < near_plane_15)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = _S4381 > far_plane_15;
        }
        if(_S4382)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4383 = scale_27.x;
        float sx_6 = (F32_exp((_S4383)));
        float _S4384 = scale_27.y;
        float sy_6 = (F32_exp((_S4384)));
        float sz_8 = scale_27.z - 0.5f * (_S4383 + _S4384);
        float x_55 = quat_28.y;
        float x2_28 = x_55 * x_55;
        float y2_28 = quat_28.z * quat_28.z;
        float z2_50 = quat_28.w * quat_28.w;
        float xy_28 = quat_28.y * quat_28.z;
        float xz_28 = quat_28.y * quat_28.w;
        float yz_28 = quat_28.z * quat_28.w;
        float wx_28 = quat_28.x * quat_28.y;
        float wy_28 = quat_28.x * quat_28.z;
        float wz_28 = quat_28.x * quat_28.w;
        Matrix<float, 3, 3>  _S4385 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_50), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_50), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
        float3  vert0_4 = mul_0(_S4385, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
        float3  vert1_4 = mul_0(_S4385, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
        float3  vert2_4 = mul_0(_S4385, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
        float3  vert0_c_7 = mul_0(R_26, vert0_4) + t_25;
        float3  vert1_c_7 = mul_0(R_26, vert1_4) + t_25;
        float3  vert2_c_7 = mul_0(R_26, vert2_4) + t_25;
        float _S4386 = length_1(vert0_c_7);
        float _S4387 = length_1(vert1_c_7);
        float _S4388 = length_1(vert2_c_7);
        if(_S4386 < near_plane_15)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = _S4386 > far_plane_15;
        }
        if(_S4382)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = _S4387 < near_plane_15;
        }
        if(_S4382)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = _S4387 > far_plane_15;
        }
        if(_S4382)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = _S4388 < near_plane_15;
        }
        if(_S4382)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = _S4388 > far_plane_15;
        }
        if(_S4382)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_6;
        float k_10;
        for(;;)
        {
            float2  _S4389 = float2 {vert0_c_7.x, vert0_c_7.y};
            float r_30 = length_0(_S4389);
            float _S4390 = vert0_c_7.z;
            float theta_26 = (F32_atan2((r_30), (_S4390)));
            if(theta_26 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_26 * theta_26 / 3.0f) / _S4390;
            }
            else
            {
                k_10 = theta_26 / r_30;
            }
            float2  uv0_7 = _S4389 * make_float2 (k_10);
            bool _S4391 = is_valid_distortion(uv0_7, dist_coeffs_37);
            bool _S4392 = !_S4391;
            _S4378 = _S4392;
            if(_S4392)
            {
                uv0_6 = uv0_7;
                break;
            }
            float u_65 = uv0_7.x;
            float v_65 = uv0_7.y;
            float r2_65 = u_65 * u_65 + v_65 * v_65;
            float2  _S4393 = uv0_7 * make_float2 (1.0f + r2_65 * ((*dist_coeffs_37)[int(0)] + r2_65 * ((*dist_coeffs_37)[int(1)] + r2_65 * ((*dist_coeffs_37)[int(2)] + r2_65 * (*dist_coeffs_37)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_37)[int(4)] * u_65 * v_65 + (*dist_coeffs_37)[int(5)] * (r2_65 + 2.0f * u_65 * u_65) + (*dist_coeffs_37)[int(6)] * r2_65, 2.0f * (*dist_coeffs_37)[int(5)] * u_65 * v_65 + (*dist_coeffs_37)[int(4)] * (r2_65 + 2.0f * v_65 * v_65) + (*dist_coeffs_37)[int(7)] * r2_65);
            float2  _S4394 = _S4393 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4393.x + (*dist_coeffs_37)[int(9)] * _S4393.y, 0.0f);
            uv0_6 = make_float2 (fx_31 * _S4394.x + cx_26, fy_31 * _S4394.y + cy_26);
            break;
        }
        float2  uv1_6;
        bool all_valid_14 = true & (!_S4378);
        for(;;)
        {
            float2  _S4395 = float2 {vert1_c_7.x, vert1_c_7.y};
            float r_31 = length_0(_S4395);
            float _S4396 = vert1_c_7.z;
            float theta_27 = (F32_atan2((r_31), (_S4396)));
            if(theta_27 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_27 * theta_27 / 3.0f) / _S4396;
            }
            else
            {
                k_10 = theta_27 / r_31;
            }
            float2  uv1_7 = _S4395 * make_float2 (k_10);
            bool _S4397 = is_valid_distortion(uv1_7, dist_coeffs_37);
            bool _S4398 = !_S4397;
            _S4379 = _S4398;
            if(_S4398)
            {
                uv1_6 = uv1_7;
                break;
            }
            float u_66 = uv1_7.x;
            float v_66 = uv1_7.y;
            float r2_66 = u_66 * u_66 + v_66 * v_66;
            float2  _S4399 = uv1_7 * make_float2 (1.0f + r2_66 * ((*dist_coeffs_37)[int(0)] + r2_66 * ((*dist_coeffs_37)[int(1)] + r2_66 * ((*dist_coeffs_37)[int(2)] + r2_66 * (*dist_coeffs_37)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_37)[int(4)] * u_66 * v_66 + (*dist_coeffs_37)[int(5)] * (r2_66 + 2.0f * u_66 * u_66) + (*dist_coeffs_37)[int(6)] * r2_66, 2.0f * (*dist_coeffs_37)[int(5)] * u_66 * v_66 + (*dist_coeffs_37)[int(4)] * (r2_66 + 2.0f * v_66 * v_66) + (*dist_coeffs_37)[int(7)] * r2_66);
            float2  _S4400 = _S4399 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4399.x + (*dist_coeffs_37)[int(9)] * _S4399.y, 0.0f);
            uv1_6 = make_float2 (fx_31 * _S4400.x + cx_26, fy_31 * _S4400.y + cy_26);
            break;
        }
        float2  uv2_6;
        bool all_valid_15 = all_valid_14 & (!_S4379);
        for(;;)
        {
            float2  _S4401 = float2 {vert2_c_7.x, vert2_c_7.y};
            float r_32 = length_0(_S4401);
            float _S4402 = vert2_c_7.z;
            float theta_28 = (F32_atan2((r_32), (_S4402)));
            if(theta_28 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_28 * theta_28 / 3.0f) / _S4402;
            }
            else
            {
                k_10 = theta_28 / r_32;
            }
            float2  uv2_7 = _S4401 * make_float2 (k_10);
            bool _S4403 = is_valid_distortion(uv2_7, dist_coeffs_37);
            bool _S4404 = !_S4403;
            _S4380 = _S4404;
            if(_S4404)
            {
                uv2_6 = uv2_7;
                break;
            }
            float u_67 = uv2_7.x;
            float v_67 = uv2_7.y;
            float r2_67 = u_67 * u_67 + v_67 * v_67;
            float2  _S4405 = uv2_7 * make_float2 (1.0f + r2_67 * ((*dist_coeffs_37)[int(0)] + r2_67 * ((*dist_coeffs_37)[int(1)] + r2_67 * ((*dist_coeffs_37)[int(2)] + r2_67 * (*dist_coeffs_37)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_37)[int(4)] * u_67 * v_67 + (*dist_coeffs_37)[int(5)] * (r2_67 + 2.0f * u_67 * u_67) + (*dist_coeffs_37)[int(6)] * r2_67, 2.0f * (*dist_coeffs_37)[int(5)] * u_67 * v_67 + (*dist_coeffs_37)[int(4)] * (r2_67 + 2.0f * v_67 * v_67) + (*dist_coeffs_37)[int(7)] * r2_67);
            float2  _S4406 = _S4405 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4405.x + (*dist_coeffs_37)[int(9)] * _S4405.y, 0.0f);
            uv2_6 = make_float2 (fx_31 * _S4406.x + cx_26, fy_31 * _S4406.y + cy_26);
            break;
        }
        if(!(all_valid_15 & (!_S4380)))
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_11 = uv1_6 - uv0_6;
        float2  e1_11 = uv2_6 - uv1_6;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(uv0_6 - uv2_6)));
        float _S4407 = uv0_6.x;
        float _S4408 = uv1_6.x;
        float _S4409 = uv2_6.x;
        float xmax_8 = (F32_max(((F32_max((_S4407), (_S4408)))), (_S4409))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S4407), (_S4408)))), (_S4409))) - offset_5;
        float _S4410 = uv0_6.y;
        float _S4411 = uv1_6.y;
        float _S4412 = uv2_6.y;
        float ymax_8 = (F32_max(((F32_max((_S4410), (_S4411)))), (_S4412))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S4410), (_S4411)))), (_S4412))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = xmin_8 >= float(image_width_22);
        }
        if(_S4382)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = ymax_8 <= 0.0f;
        }
        if(_S4382)
        {
            _S4382 = true;
        }
        else
        {
            _S4382 = ymin_8 >= float(image_height_22);
        }
        if(_S4382)
        {
            _S4382 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S4382 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S4382 = false;
                }
                if(_S4382)
                {
                    _S4382 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S4382 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S4382 = false;
                    }
                }
            }
            else
            {
                _S4382 = false;
            }
        }
        if(_S4382)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4413 = mean_27 - - mul_0(transpose_0(R_26), t_25);
        float _S4414 = _S4413.x;
        float _S4415 = _S4413.y;
        float _S4416 = _S4413.z;
        float norm_15 = (F32_sqrt((_S4414 * _S4414 + _S4415 * _S4415 + _S4416 * _S4416)));
        float x_56 = _S4414 / norm_15;
        float y_25 = _S4415 / norm_15;
        float z_22 = _S4416 / norm_15;
        float z2_51 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_56 * x_56 - y_25 * y_25;
        float fS1_22 = 2.0f * x_56 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_51 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_56) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_51 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_56) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_56 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_51 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_56) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_56 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S4417 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S4417);
        float3  _S4418 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S4419 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S4418 + _S4419 + make_float3 (0.5f), _S4417);
        (*rgbs_1)[int(2)] = max_0(_S4418 - _S4419 + make_float3 (0.5f), _S4417);
        (*verts_1)[int(0)] = vert0_4;
        (*verts_1)[int(1)] = vert1_4;
        (*verts_1)[int(2)] = vert2_4;
        float3  _S4420 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S4420 * make_float3 (float(- (F32_sign((dot_0(_S4420, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_27, float3  t_26, float fx_32, float fy_32, float cx_27, float cy_27, FixedArray<float, 10>  * dist_coeffs_38, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_27, mean_28) + t_26;
    float _S4421 = scale_28.x;
    float sx_7 = (F32_exp((_S4421)));
    float _S4422 = scale_28.y;
    float sy_7 = (F32_exp((_S4422)));
    float sz_9 = scale_28.z - 0.5f * (_S4421 + _S4422);
    float x_57 = quat_29.y;
    float x2_29 = x_57 * x_57;
    float y2_29 = quat_29.z * quat_29.z;
    float z2_52 = quat_29.w * quat_29.w;
    float xy_29 = quat_29.y * quat_29.z;
    float xz_29 = quat_29.y * quat_29.w;
    float yz_29 = quat_29.z * quat_29.w;
    float wx_29 = quat_29.x * quat_29.y;
    float wy_29 = quat_29.x * quat_29.z;
    float wz_29 = quat_29.x * quat_29.w;
    Matrix<float, 3, 3>  _S4423 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_52), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_52), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S4423, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S4423, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S4423, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_8 = mul_0(R_27, vert0_5) + t_26;
    float3  vert1_c_8 = mul_0(R_27, vert1_5) + t_26;
    float3  vert2_c_8 = mul_0(R_27, vert2_5) + t_26;
    float2  _S4424 = float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z);
    float u_68 = _S4424.x;
    float v_68 = _S4424.y;
    float r2_68 = u_68 * u_68 + v_68 * v_68;
    float _S4425 = 2.0f * (*dist_coeffs_38)[int(4)];
    float _S4426 = 2.0f * (*dist_coeffs_38)[int(5)];
    float2  _S4427 = _S4424 * make_float2 (1.0f + r2_68 * ((*dist_coeffs_38)[int(0)] + r2_68 * ((*dist_coeffs_38)[int(1)] + r2_68 * ((*dist_coeffs_38)[int(2)] + r2_68 * (*dist_coeffs_38)[int(3)])))) + make_float2 (_S4425 * u_68 * v_68 + (*dist_coeffs_38)[int(5)] * (r2_68 + 2.0f * u_68 * u_68) + (*dist_coeffs_38)[int(6)] * r2_68, _S4426 * u_68 * v_68 + (*dist_coeffs_38)[int(4)] * (r2_68 + 2.0f * v_68 * v_68) + (*dist_coeffs_38)[int(7)] * r2_68);
    float2  _S4428 = _S4427 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4427.x + (*dist_coeffs_38)[int(9)] * _S4427.y, 0.0f);
    float _S4429 = fx_32 * _S4428.x + cx_27;
    float _S4430 = fy_32 * _S4428.y + cy_27;
    float2  uv0_8 = make_float2 (_S4429, _S4430);
    float2  _S4431 = float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z);
    float u_69 = _S4431.x;
    float v_69 = _S4431.y;
    float r2_69 = u_69 * u_69 + v_69 * v_69;
    float2  _S4432 = _S4431 * make_float2 (1.0f + r2_69 * ((*dist_coeffs_38)[int(0)] + r2_69 * ((*dist_coeffs_38)[int(1)] + r2_69 * ((*dist_coeffs_38)[int(2)] + r2_69 * (*dist_coeffs_38)[int(3)])))) + make_float2 (_S4425 * u_69 * v_69 + (*dist_coeffs_38)[int(5)] * (r2_69 + 2.0f * u_69 * u_69) + (*dist_coeffs_38)[int(6)] * r2_69, _S4426 * u_69 * v_69 + (*dist_coeffs_38)[int(4)] * (r2_69 + 2.0f * v_69 * v_69) + (*dist_coeffs_38)[int(7)] * r2_69);
    float2  _S4433 = _S4432 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4432.x + (*dist_coeffs_38)[int(9)] * _S4432.y, 0.0f);
    float _S4434 = fx_32 * _S4433.x + cx_27;
    float _S4435 = fy_32 * _S4433.y + cy_27;
    float2  uv1_8 = make_float2 (_S4434, _S4435);
    float2  _S4436 = float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z);
    float u_70 = _S4436.x;
    float v_70 = _S4436.y;
    float r2_70 = u_70 * u_70 + v_70 * v_70;
    float2  _S4437 = _S4436 * make_float2 (1.0f + r2_70 * ((*dist_coeffs_38)[int(0)] + r2_70 * ((*dist_coeffs_38)[int(1)] + r2_70 * ((*dist_coeffs_38)[int(2)] + r2_70 * (*dist_coeffs_38)[int(3)])))) + make_float2 (_S4425 * u_70 * v_70 + (*dist_coeffs_38)[int(5)] * (r2_70 + 2.0f * u_70 * u_70) + (*dist_coeffs_38)[int(6)] * r2_70, _S4426 * u_70 * v_70 + (*dist_coeffs_38)[int(4)] * (r2_70 + 2.0f * v_70 * v_70) + (*dist_coeffs_38)[int(7)] * r2_70);
    float2  _S4438 = _S4437 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4437.x + (*dist_coeffs_38)[int(9)] * _S4437.y, 0.0f);
    float _S4439 = fx_32 * _S4438.x + cx_27;
    float _S4440 = fy_32 * _S4438.y + cy_27;
    float2  uv2_8 = make_float2 (_S4439, _S4440);
    float2  e0_12 = uv1_8 - uv0_8;
    float2  e1_12 = uv2_8 - uv1_8;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(uv0_8 - uv2_8)));
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4429), (_S4434)))), (_S4439))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S4430), (_S4435)))), (_S4440))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4429), (_S4434)))), (_S4439))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4430), (_S4435)))), (_S4440))) + offset_6)))));
    *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4441 = mean_28 - - mul_0(transpose_0(R_27), t_26);
    float _S4442 = _S4441.x;
    float _S4443 = _S4441.y;
    float _S4444 = _S4441.z;
    float norm_16 = (F32_sqrt((_S4442 * _S4442 + _S4443 * _S4443 + _S4444 * _S4444)));
    float x_58 = _S4442 / norm_16;
    float y_26 = _S4443 / norm_16;
    float z_23 = _S4444 / norm_16;
    float z2_53 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_58 * x_58 - y_26 * y_26;
    float fS1_23 = 2.0f * x_58 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_53 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_58) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_53 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_58) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_58 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_53 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_58) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_58 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S4445 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S4445);
    float3  _S4446 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S4447 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S4446 + _S4447 + make_float3 (0.5f), _S4445);
    (*rgbs_2)[int(2)] = max_0(_S4446 - _S4447 + make_float3 (0.5f), _S4445);
    (*verts_2)[int(0)] = vert0_5;
    (*verts_2)[int(1)] = vert1_5;
    (*verts_2)[int(2)] = vert2_5;
    float3  _S4448 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S4448 * make_float3 (float(- (F32_sign((dot_0(_S4448, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_28, float3  t_27, float fx_33, float fy_33, float cx_28, float cy_28, FixedArray<float, 10>  * dist_coeffs_39, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_28, mean_29) + t_27;
    float _S4449 = scale_29.x;
    float sx_8 = (F32_exp((_S4449)));
    float _S4450 = scale_29.y;
    float sy_8 = (F32_exp((_S4450)));
    float sz_10 = scale_29.z - 0.5f * (_S4449 + _S4450);
    float x_59 = quat_30.y;
    float x2_30 = x_59 * x_59;
    float y2_30 = quat_30.z * quat_30.z;
    float z2_54 = quat_30.w * quat_30.w;
    float xy_30 = quat_30.y * quat_30.z;
    float xz_30 = quat_30.y * quat_30.w;
    float yz_30 = quat_30.z * quat_30.w;
    float wx_30 = quat_30.x * quat_30.y;
    float wy_30 = quat_30.x * quat_30.z;
    float wz_30 = quat_30.x * quat_30.w;
    Matrix<float, 3, 3>  _S4451 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_54), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_54), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  vert0_6 = mul_0(_S4451, make_float3 (sx_8, 0.0f, 0.0f)) + mean_29;
    float3  vert1_6 = mul_0(_S4451, make_float3 (sx_8 * (-0.5f + sz_10), sy_8, 0.0f)) + mean_29;
    float3  vert2_6 = mul_0(_S4451, make_float3 (sx_8 * (-0.5f - sz_10), - sy_8, 0.0f)) + mean_29;
    float3  vert0_c_9 = mul_0(R_28, vert0_6) + t_27;
    float3  vert1_c_9 = mul_0(R_28, vert1_6) + t_27;
    float3  vert2_c_9 = mul_0(R_28, vert2_6) + t_27;
    float2  _S4452 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_33 = length_0(_S4452);
    float _S4453 = vert0_c_9.z;
    float theta_29 = (F32_atan2((r_33), (_S4453)));
    float k_11;
    if(theta_29 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_29 * theta_29 / 3.0f) / _S4453;
    }
    else
    {
        k_11 = theta_29 / r_33;
    }
    float2  _S4454 = _S4452 * make_float2 (k_11);
    float u_71 = _S4454.x;
    float v_71 = _S4454.y;
    float r2_71 = u_71 * u_71 + v_71 * v_71;
    float _S4455 = 2.0f * (*dist_coeffs_39)[int(4)];
    float _S4456 = 2.0f * (*dist_coeffs_39)[int(5)];
    float2  _S4457 = _S4454 * make_float2 (1.0f + r2_71 * ((*dist_coeffs_39)[int(0)] + r2_71 * ((*dist_coeffs_39)[int(1)] + r2_71 * ((*dist_coeffs_39)[int(2)] + r2_71 * (*dist_coeffs_39)[int(3)])))) + make_float2 (_S4455 * u_71 * v_71 + (*dist_coeffs_39)[int(5)] * (r2_71 + 2.0f * u_71 * u_71) + (*dist_coeffs_39)[int(6)] * r2_71, _S4456 * u_71 * v_71 + (*dist_coeffs_39)[int(4)] * (r2_71 + 2.0f * v_71 * v_71) + (*dist_coeffs_39)[int(7)] * r2_71);
    float2  _S4458 = _S4457 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4457.x + (*dist_coeffs_39)[int(9)] * _S4457.y, 0.0f);
    float _S4459 = fx_33 * _S4458.x + cx_28;
    float _S4460 = fy_33 * _S4458.y + cy_28;
    float2  uv0_9 = make_float2 (_S4459, _S4460);
    float2  _S4461 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_34 = length_0(_S4461);
    float _S4462 = vert1_c_9.z;
    float theta_30 = (F32_atan2((r_34), (_S4462)));
    if(theta_30 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_30 * theta_30 / 3.0f) / _S4462;
    }
    else
    {
        k_11 = theta_30 / r_34;
    }
    float2  _S4463 = _S4461 * make_float2 (k_11);
    float u_72 = _S4463.x;
    float v_72 = _S4463.y;
    float r2_72 = u_72 * u_72 + v_72 * v_72;
    float2  _S4464 = _S4463 * make_float2 (1.0f + r2_72 * ((*dist_coeffs_39)[int(0)] + r2_72 * ((*dist_coeffs_39)[int(1)] + r2_72 * ((*dist_coeffs_39)[int(2)] + r2_72 * (*dist_coeffs_39)[int(3)])))) + make_float2 (_S4455 * u_72 * v_72 + (*dist_coeffs_39)[int(5)] * (r2_72 + 2.0f * u_72 * u_72) + (*dist_coeffs_39)[int(6)] * r2_72, _S4456 * u_72 * v_72 + (*dist_coeffs_39)[int(4)] * (r2_72 + 2.0f * v_72 * v_72) + (*dist_coeffs_39)[int(7)] * r2_72);
    float2  _S4465 = _S4464 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4464.x + (*dist_coeffs_39)[int(9)] * _S4464.y, 0.0f);
    float _S4466 = fx_33 * _S4465.x + cx_28;
    float _S4467 = fy_33 * _S4465.y + cy_28;
    float2  uv1_9 = make_float2 (_S4466, _S4467);
    float2  _S4468 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_35 = length_0(_S4468);
    float _S4469 = vert2_c_9.z;
    float theta_31 = (F32_atan2((r_35), (_S4469)));
    if(theta_31 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_31 * theta_31 / 3.0f) / _S4469;
    }
    else
    {
        k_11 = theta_31 / r_35;
    }
    float2  _S4470 = _S4468 * make_float2 (k_11);
    float u_73 = _S4470.x;
    float v_73 = _S4470.y;
    float r2_73 = u_73 * u_73 + v_73 * v_73;
    float2  _S4471 = _S4470 * make_float2 (1.0f + r2_73 * ((*dist_coeffs_39)[int(0)] + r2_73 * ((*dist_coeffs_39)[int(1)] + r2_73 * ((*dist_coeffs_39)[int(2)] + r2_73 * (*dist_coeffs_39)[int(3)])))) + make_float2 (_S4455 * u_73 * v_73 + (*dist_coeffs_39)[int(5)] * (r2_73 + 2.0f * u_73 * u_73) + (*dist_coeffs_39)[int(6)] * r2_73, _S4456 * u_73 * v_73 + (*dist_coeffs_39)[int(4)] * (r2_73 + 2.0f * v_73 * v_73) + (*dist_coeffs_39)[int(7)] * r2_73);
    float2  _S4472 = _S4471 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4471.x + (*dist_coeffs_39)[int(9)] * _S4471.y, 0.0f);
    float _S4473 = fx_33 * _S4472.x + cx_28;
    float _S4474 = fy_33 * _S4472.y + cy_28;
    float2  uv2_9 = make_float2 (_S4473, _S4474);
    float2  e0_13 = uv1_9 - uv0_9;
    float2  e1_13 = uv2_9 - uv1_9;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(uv0_9 - uv2_9)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4459), (_S4466)))), (_S4473))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S4460), (_S4467)))), (_S4474))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4459), (_S4466)))), (_S4473))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4460), (_S4467)))), (_S4474))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4475 = mean_29 - - mul_0(transpose_0(R_28), t_27);
    float _S4476 = _S4475.x;
    float _S4477 = _S4475.y;
    float _S4478 = _S4475.z;
    float norm_17 = (F32_sqrt((_S4476 * _S4476 + _S4477 * _S4477 + _S4478 * _S4478)));
    float x_60 = _S4476 / norm_17;
    float y_27 = _S4477 / norm_17;
    float z_24 = _S4478 / norm_17;
    float z2_55 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_60 * x_60 - y_27 * y_27;
    float fS1_24 = 2.0f * x_60 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_55 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_60) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_55 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_60) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_60 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_55 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_60) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_60 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S4479 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S4479);
    float3  _S4480 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S4481 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S4480 + _S4481 + make_float3 (0.5f), _S4479);
    (*rgbs_3)[int(2)] = max_0(_S4480 - _S4481 + make_float3 (0.5f), _S4479);
    (*verts_3)[int(0)] = vert0_6;
    (*verts_3)[int(1)] = vert1_6;
    (*verts_3)[int(2)] = vert2_6;
    float3  _S4482 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S4482 * make_float3 (float(- (F32_sign((dot_0(_S4482, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_29, float3  t_28, float fx_34, float fy_34, float cx_29, float cy_29, FixedArray<float, 10>  * dist_coeffs_40, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_1(R_29, mean_30) + t_28;
    float _S4483 = scale_30.x;
    float _S4484 = s_primal_ctx_exp_1(_S4483);
    float _S4485 = scale_30.y;
    float _S4486 = s_primal_ctx_exp_1(_S4485);
    float sz_11 = scale_30.z - 0.5f * (_S4483 + _S4485);
    float _S4487 = quat_31.y;
    float x2_31 = _S4487 * _S4487;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_56 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4488 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_56), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_56), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4489 = make_float3 (_S4484, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_1(_S4488, _S4489) + mean_30;
    float _S4490 = -0.5f + sz_11;
    float3  _S4491 = make_float3 (_S4484 * _S4490, _S4486, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_1(_S4488, _S4491) + mean_30;
    float _S4492 = -0.5f - sz_11;
    float3  _S4493 = make_float3 (_S4484 * _S4492, - _S4486, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_1(_S4488, _S4493) + mean_30;
    float3  vert0_c_10 = s_primal_ctx_mul_1(R_29, vert0_7) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_1(R_29, vert1_7) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_1(R_29, vert2_7) + t_28;
    float2  _S4494 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S4495 = vert0_c_10.z;
    float2  _S4496 = make_float2 (_S4495);
    float2  _S4497 = _S4494 / make_float2 (_S4495);
    float2  _S4498 = make_float2 (_S4495 * _S4495);
    float u_74 = _S4497.x;
    float v_74 = _S4497.y;
    float r2_74 = u_74 * u_74 + v_74 * v_74;
    float _S4499 = (*dist_coeffs_40)[int(2)] + r2_74 * (*dist_coeffs_40)[int(3)];
    float _S4500 = (*dist_coeffs_40)[int(1)] + r2_74 * _S4499;
    float _S4501 = (*dist_coeffs_40)[int(0)] + r2_74 * _S4500;
    float radial_8 = 1.0f + r2_74 * _S4501;
    float2  _S4502 = make_float2 (radial_8);
    float _S4503 = 2.0f * (*dist_coeffs_40)[int(4)];
    float _S4504 = _S4503 * u_74;
    float _S4505 = 2.0f * u_74;
    float _S4506 = 2.0f * (*dist_coeffs_40)[int(5)];
    float _S4507 = _S4506 * u_74;
    float _S4508 = 2.0f * v_74;
    float2  _S4509 = _S4497 * make_float2 (radial_8) + make_float2 (_S4504 * v_74 + (*dist_coeffs_40)[int(5)] * (r2_74 + _S4505 * u_74) + (*dist_coeffs_40)[int(6)] * r2_74, _S4507 * v_74 + (*dist_coeffs_40)[int(4)] * (r2_74 + _S4508 * v_74) + (*dist_coeffs_40)[int(7)] * r2_74);
    float2  _S4510 = _S4509 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4509.x + (*dist_coeffs_40)[int(9)] * _S4509.y, 0.0f);
    float _S4511 = fx_34 * _S4510.x + cx_29;
    float _S4512 = fy_34 * _S4510.y + cy_29;
    float2  uv0_10 = make_float2 (_S4511, _S4512);
    float2  _S4513 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S4514 = vert1_c_10.z;
    float2  _S4515 = make_float2 (_S4514);
    float2  _S4516 = _S4513 / make_float2 (_S4514);
    float2  _S4517 = make_float2 (_S4514 * _S4514);
    float u_75 = _S4516.x;
    float v_75 = _S4516.y;
    float r2_75 = u_75 * u_75 + v_75 * v_75;
    float _S4518 = (*dist_coeffs_40)[int(2)] + r2_75 * (*dist_coeffs_40)[int(3)];
    float _S4519 = (*dist_coeffs_40)[int(1)] + r2_75 * _S4518;
    float _S4520 = (*dist_coeffs_40)[int(0)] + r2_75 * _S4519;
    float radial_9 = 1.0f + r2_75 * _S4520;
    float2  _S4521 = make_float2 (radial_9);
    float _S4522 = _S4503 * u_75;
    float _S4523 = 2.0f * u_75;
    float _S4524 = _S4506 * u_75;
    float _S4525 = 2.0f * v_75;
    float2  _S4526 = _S4516 * make_float2 (radial_9) + make_float2 (_S4522 * v_75 + (*dist_coeffs_40)[int(5)] * (r2_75 + _S4523 * u_75) + (*dist_coeffs_40)[int(6)] * r2_75, _S4524 * v_75 + (*dist_coeffs_40)[int(4)] * (r2_75 + _S4525 * v_75) + (*dist_coeffs_40)[int(7)] * r2_75);
    float2  _S4527 = _S4526 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4526.x + (*dist_coeffs_40)[int(9)] * _S4526.y, 0.0f);
    float _S4528 = fx_34 * _S4527.x + cx_29;
    float _S4529 = fy_34 * _S4527.y + cy_29;
    float2  uv1_10 = make_float2 (_S4528, _S4529);
    float2  _S4530 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S4531 = vert2_c_10.z;
    float2  _S4532 = make_float2 (_S4531);
    float2  _S4533 = _S4530 / make_float2 (_S4531);
    float2  _S4534 = make_float2 (_S4531 * _S4531);
    float u_76 = _S4533.x;
    float v_76 = _S4533.y;
    float r2_76 = u_76 * u_76 + v_76 * v_76;
    float _S4535 = (*dist_coeffs_40)[int(2)] + r2_76 * (*dist_coeffs_40)[int(3)];
    float _S4536 = (*dist_coeffs_40)[int(1)] + r2_76 * _S4535;
    float _S4537 = (*dist_coeffs_40)[int(0)] + r2_76 * _S4536;
    float radial_10 = 1.0f + r2_76 * _S4537;
    float2  _S4538 = make_float2 (radial_10);
    float _S4539 = _S4503 * u_76;
    float _S4540 = 2.0f * u_76;
    float _S4541 = _S4506 * u_76;
    float _S4542 = 2.0f * v_76;
    float2  _S4543 = _S4533 * make_float2 (radial_10) + make_float2 (_S4539 * v_76 + (*dist_coeffs_40)[int(5)] * (r2_76 + _S4540 * u_76) + (*dist_coeffs_40)[int(6)] * r2_76, _S4541 * v_76 + (*dist_coeffs_40)[int(4)] * (r2_76 + _S4542 * v_76) + (*dist_coeffs_40)[int(7)] * r2_76);
    float2  _S4544 = _S4543 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4543.x + (*dist_coeffs_40)[int(9)] * _S4543.y, 0.0f);
    float _S4545 = fx_34 * _S4544.x + cx_29;
    float _S4546 = fy_34 * _S4544.y + cy_29;
    float2  uv2_10 = make_float2 (_S4545, _S4546);
    float2  e0_14 = uv1_10 - uv0_10;
    float2  e1_14 = uv2_10 - uv1_10;
    float2  e2_6 = uv0_10 - uv2_10;
    float _S4547 = e0_14.x;
    float _S4548 = e1_14.y;
    float _S4549 = e0_14.y;
    float _S4550 = e1_14.x;
    float _S4551 = _S4547 * _S4548 - _S4549 * _S4550;
    float _S4552 = 1.0f - hardness_14.y;
    float _S4553 = -1.0f / _S4552;
    float _S4554 = _S4552 * _S4552;
    float _S4555 = s_primal_ctx_max_0(_S4511, _S4528);
    float _S4556 = s_primal_ctx_min_0(_S4511, _S4528);
    float _S4557 = s_primal_ctx_max_0(_S4512, _S4529);
    float _S4558 = s_primal_ctx_min_0(_S4512, _S4529);
    float3  _S4559 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S4560 = length_1(_S4559) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4561 = transpose_0(R_29);
    float3  _S4562 = mean_30 - - s_primal_ctx_mul_1(_S4561, t_28);
    float _S4563 = _S4562.x;
    float _S4564 = _S4562.y;
    float _S4565 = _S4562.z;
    float _S4566 = _S4563 * _S4563 + _S4564 * _S4564 + _S4565 * _S4565;
    float _S4567 = s_primal_ctx_sqrt_0(_S4566);
    float x_61 = _S4563 / _S4567;
    float3  _S4568 = make_float3 (x_61);
    float _S4569 = _S4567 * _S4567;
    float y_28 = _S4564 / _S4567;
    float z_25 = _S4565 / _S4567;
    float3  _S4570 = make_float3 (z_25);
    float _S4571 = - y_28;
    float3  _S4572 = make_float3 (_S4571);
    float z2_57 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_61 * x_61 - y_28 * y_28;
    float _S4573 = 2.0f * x_61;
    float fS1_25 = _S4573 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_57 - 0.31539157032966614f;
    float3  _S4574 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_61;
    float3  _S4575 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S4576 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S4577 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S4578 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_57 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S4579 = 1.86588168144226074f * z2_57 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S4579;
    float3  _S4580 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_61;
    float3  _S4581 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S4582 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S4583 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S4584 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_61 * fC1_25 - y_28 * fS1_25);
    float3  _S4585 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_61 * fS1_25 + y_28 * fC1_25);
    float3  _S4586 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4571) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_61) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S4587 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S4588 = make_float3 (0.0f);
    float3  _S4589 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S4590 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4591 = make_float3 (_S4590);
    float3  _S4592 = (*ch_coeffs_10)[int(1)] * make_float3 (_S4590);
    float3  _S4593 = _S4589 + _S4592 + make_float3 (0.5f);
    float3  _S4594 = _S4589 - _S4592 + make_float3 (0.5f);
    float3  _S4595 = vert1_c_10 - vert0_c_10;
    float3  _S4596 = vert2_c_10 - vert0_c_10;
    float3  _S4597 = s_primal_ctx_cross_0(_S4595, _S4596);
    float3  _S4598 = normalize_0(_S4597);
    float3  _S4599 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4598, mean_c_25)))))) * v_normal_2;
    float3  _S4600 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4601;
    (&_S4601)->primal_0 = _S4598;
    (&_S4601)->differential_0 = _S4600;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4602;
    (&_S4602)->primal_0 = mean_c_25;
    (&_S4602)->differential_0 = _S4600;
    s_bwd_prop_dot_0(&_S4601, &_S4602, 0.0f);
    float3  _S4603 = _S4599 + _S4601.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4604;
    (&_S4604)->primal_0 = _S4597;
    (&_S4604)->differential_0 = _S4600;
    s_bwd_normalize_impl_0(&_S4604, _S4603);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4605;
    (&_S4605)->primal_0 = _S4595;
    (&_S4605)->differential_0 = _S4600;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4606;
    (&_S4606)->primal_0 = _S4596;
    (&_S4606)->differential_0 = _S4600;
    s_bwd_prop_cross_0(&_S4605, &_S4606, _S4604.differential_0);
    float3  _S4607 = - _S4606.differential_0;
    float3  _S4608 = - _S4605.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4609;
    (&_S4609)->primal_0 = _S4594;
    (&_S4609)->differential_0 = _S4600;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4610;
    (&_S4610)->primal_0 = _S4588;
    (&_S4610)->differential_0 = _S4600;
    s_bwd_prop_max_0(&_S4609, &_S4610, (*v_rgbs_0)[int(2)]);
    float3  _S4611 = - _S4609.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4612;
    (&_S4612)->primal_0 = _S4593;
    (&_S4612)->differential_0 = _S4600;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4613;
    (&_S4613)->primal_0 = _S4588;
    (&_S4613)->differential_0 = _S4600;
    s_bwd_prop_max_0(&_S4612, &_S4613, (*v_rgbs_0)[int(1)]);
    float3  _S4614 = _S4591 * (_S4611 + _S4612.differential_0);
    float3  _S4615 = _S4609.differential_0 + _S4612.differential_0;
    float3  _S4616 = make_float3 (0.5f) * - _S4615;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4617;
    (&_S4617)->primal_0 = _S4587;
    (&_S4617)->differential_0 = _S4600;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4618;
    (&_S4618)->primal_0 = _S4588;
    (&_S4618)->differential_0 = _S4600;
    s_bwd_prop_max_0(&_S4617, &_S4618, (*v_rgbs_0)[int(0)]);
    float3  _S4619 = _S4616 + _S4617.differential_0;
    float3  _S4620 = _S4615 + _S4617.differential_0;
    float3  _S4621 = _S4585 * _S4620;
    float3  _S4622 = (*sh_coeffs_25)[int(15)] * _S4620;
    float3  _S4623 = _S4583 * _S4620;
    float3  _S4624 = (*sh_coeffs_25)[int(14)] * _S4620;
    float3  _S4625 = _S4581 * _S4620;
    float3  _S4626 = (*sh_coeffs_25)[int(13)] * _S4620;
    float3  _S4627 = _S4580 * _S4620;
    float3  _S4628 = (*sh_coeffs_25)[int(12)] * _S4620;
    float3  _S4629 = _S4582 * _S4620;
    float3  _S4630 = (*sh_coeffs_25)[int(11)] * _S4620;
    float3  _S4631 = _S4584 * _S4620;
    float3  _S4632 = (*sh_coeffs_25)[int(10)] * _S4620;
    float3  _S4633 = _S4586 * _S4620;
    float3  _S4634 = (*sh_coeffs_25)[int(9)] * _S4620;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S4634.x + _S4634.y + _S4634.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S4622.x + _S4622.y + _S4622.z);
    float _S4635 = _S4632.x + _S4632.y + _S4632.z;
    float _S4636 = _S4624.x + _S4624.y + _S4624.z;
    float _S4637 = _S4630.x + _S4630.y + _S4630.z;
    float _S4638 = _S4626.x + _S4626.y + _S4626.z;
    float _S4639 = _S4628.x + _S4628.y + _S4628.z;
    float _S4640 = - s_diff_fC2_T_7;
    float3  _S4641 = _S4577 * _S4620;
    float3  _S4642 = (*sh_coeffs_25)[int(8)] * _S4620;
    float3  _S4643 = _S4575 * _S4620;
    float3  _S4644 = (*sh_coeffs_25)[int(7)] * _S4620;
    float3  _S4645 = _S4574 * _S4620;
    float3  _S4646 = (*sh_coeffs_25)[int(6)] * _S4620;
    float3  _S4647 = _S4576 * _S4620;
    float3  _S4648 = (*sh_coeffs_25)[int(5)] * _S4620;
    float3  _S4649 = _S4578 * _S4620;
    float3  _S4650 = (*sh_coeffs_25)[int(4)] * _S4620;
    float _S4651 = _S4648.x + _S4648.y + _S4648.z;
    float _S4652 = _S4644.x + _S4644.y + _S4644.z;
    float _S4653 = fTmp1B_25 * _S4635 + x_61 * s_diff_fS2_T_7 + y_28 * _S4640 + 0.54627424478530884f * (_S4650.x + _S4650.y + _S4650.z);
    float _S4654 = fTmp1B_25 * _S4636 + y_28 * s_diff_fS2_T_7 + x_61 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S4642.x + _S4642.y + _S4642.z);
    float _S4655 = y_28 * - _S4654;
    float _S4656 = x_61 * _S4654;
    float _S4657 = z_25 * (1.86588168144226074f * (z_25 * _S4639) + -2.28522896766662598f * (y_28 * _S4637 + x_61 * _S4638) + 0.94617468118667603f * (_S4646.x + _S4646.y + _S4646.z));
    float3  _S4658 = make_float3 (0.48860251903533936f) * _S4620;
    float3  _S4659 = - _S4658;
    float3  _S4660 = _S4568 * _S4659;
    float3  _S4661 = (*sh_coeffs_25)[int(3)] * _S4659;
    float3  _S4662 = _S4570 * _S4658;
    float3  _S4663 = (*sh_coeffs_25)[int(2)] * _S4658;
    float3  _S4664 = _S4572 * _S4658;
    float3  _S4665 = (*sh_coeffs_25)[int(1)] * _S4658;
    float _S4666 = (_S4579 * _S4639 + 1.44530570507049561f * (fS1_25 * _S4635 + fC1_25 * _S4636) + -1.09254848957061768f * (y_28 * _S4651 + x_61 * _S4652) + _S4657 + _S4657 + _S4663.x + _S4663.y + _S4663.z) / _S4569;
    float _S4667 = _S4567 * _S4666;
    float _S4668 = (fTmp0C_25 * _S4637 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S4640 + fTmp0B_25 * _S4651 + _S4573 * _S4653 + _S4655 + _S4655 + - (_S4665.x + _S4665.y + _S4665.z)) / _S4569;
    float _S4669 = _S4567 * _S4668;
    float _S4670 = (fTmp0C_25 * _S4638 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4652 + 2.0f * (y_28 * _S4653) + _S4656 + _S4656 + _S4661.x + _S4661.y + _S4661.z) / _S4569;
    float _S4671 = _S4567 * _S4670;
    float _S4672 = _S4565 * - _S4666 + _S4564 * - _S4668 + _S4563 * - _S4670;
    DiffPair_float_0 _S4673;
    (&_S4673)->primal_0 = _S4566;
    (&_S4673)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4673, _S4672);
    float _S4674 = _S4565 * _S4673.differential_0;
    float _S4675 = _S4564 * _S4673.differential_0;
    float _S4676 = _S4563 * _S4673.differential_0;
    float3  _S4677 = make_float3 (0.282094806432724f) * _S4620;
    float3  _S4678 = make_float3 (_S4671 + _S4676 + _S4676, _S4669 + _S4675 + _S4675, _S4667 + _S4674 + _S4674);
    float3  _S4679 = - - _S4678;
    Matrix<float, 3, 3>  _S4680 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4681;
    (&_S4681)->primal_0 = _S4561;
    (&_S4681)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4682;
    (&_S4682)->primal_0 = t_28;
    (&_S4682)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4681, &_S4682, _S4679);
    Matrix<float, 3, 3>  _S4683 = transpose_0(_S4681.differential_0);
    DiffPair_float_0 _S4684;
    (&_S4684)->primal_0 = _S4560;
    (&_S4684)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4684, v_depth_9);
    float _S4685 = 0.3333333432674408f * _S4684.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4686;
    (&_S4686)->primal_0 = _S4559;
    (&_S4686)->differential_0 = _S4600;
    s_bwd_length_impl_0(&_S4686, _S4685);
    DiffPair_float_0 _S4687;
    (&_S4687)->primal_0 = _S4558;
    (&_S4687)->differential_0 = 0.0f;
    DiffPair_float_0 _S4688;
    (&_S4688)->primal_0 = _S4546;
    (&_S4688)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4687, &_S4688, 0.0f);
    DiffPair_float_0 _S4689;
    (&_S4689)->primal_0 = _S4512;
    (&_S4689)->differential_0 = 0.0f;
    DiffPair_float_0 _S4690;
    (&_S4690)->primal_0 = _S4529;
    (&_S4690)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4689, &_S4690, _S4687.differential_0);
    DiffPair_float_0 _S4691;
    (&_S4691)->primal_0 = _S4557;
    (&_S4691)->differential_0 = 0.0f;
    DiffPair_float_0 _S4692;
    (&_S4692)->primal_0 = _S4546;
    (&_S4692)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4691, &_S4692, 0.0f);
    DiffPair_float_0 _S4693;
    (&_S4693)->primal_0 = _S4512;
    (&_S4693)->differential_0 = 0.0f;
    DiffPair_float_0 _S4694;
    (&_S4694)->primal_0 = _S4529;
    (&_S4694)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4693, &_S4694, _S4691.differential_0);
    DiffPair_float_0 _S4695;
    (&_S4695)->primal_0 = _S4556;
    (&_S4695)->differential_0 = 0.0f;
    DiffPair_float_0 _S4696;
    (&_S4696)->primal_0 = _S4545;
    (&_S4696)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4695, &_S4696, 0.0f);
    DiffPair_float_0 _S4697;
    (&_S4697)->primal_0 = _S4511;
    (&_S4697)->differential_0 = 0.0f;
    DiffPair_float_0 _S4698;
    (&_S4698)->primal_0 = _S4528;
    (&_S4698)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4697, &_S4698, _S4695.differential_0);
    DiffPair_float_0 _S4699;
    (&_S4699)->primal_0 = _S4555;
    (&_S4699)->differential_0 = 0.0f;
    DiffPair_float_0 _S4700;
    (&_S4700)->primal_0 = _S4545;
    (&_S4700)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4699, &_S4700, 0.0f);
    DiffPair_float_0 _S4701;
    (&_S4701)->primal_0 = _S4511;
    (&_S4701)->differential_0 = 0.0f;
    DiffPair_float_0 _S4702;
    (&_S4702)->primal_0 = _S4528;
    (&_S4702)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4701, &_S4702, _S4699.differential_0);
    DiffPair_float_0 _S4703;
    (&_S4703)->primal_0 = _S4553;
    (&_S4703)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4703, 0.0f);
    float _S4704 = - (-1.0f * - (_S4703.differential_0 / _S4554));
    float2  _S4705 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4706;
    (&_S4706)->primal_0 = e2_6;
    (&_S4706)->differential_0 = _S4705;
    s_bwd_length_impl_1(&_S4706, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4707;
    (&_S4707)->primal_0 = e1_14;
    (&_S4707)->differential_0 = _S4705;
    s_bwd_length_impl_1(&_S4707, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4708;
    (&_S4708)->primal_0 = e0_14;
    (&_S4708)->differential_0 = _S4705;
    s_bwd_length_impl_1(&_S4708, -0.0f);
    DiffPair_float_0 _S4709;
    (&_S4709)->primal_0 = _S4551;
    (&_S4709)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4709, 0.0f);
    float _S4710 = - _S4709.differential_0;
    float2  _S4711 = _S4707.differential_0 + make_float2 (_S4549 * _S4710, _S4547 * _S4709.differential_0);
    float2  _S4712 = _S4708.differential_0 + make_float2 (_S4548 * _S4709.differential_0, _S4550 * _S4710);
    float2  _S4713 = - _S4706.differential_0 + _S4711;
    float _S4714 = fx_34 * (_S4696.differential_0 + _S4700.differential_0 + _S4713.x);
    float2  _S4715 = make_float2 (_S4714, fy_34 * (_S4688.differential_0 + _S4692.differential_0 + _S4713.y)) + make_float2 ((*dist_coeffs_40)[int(8)] * _S4714, (*dist_coeffs_40)[int(9)] * _S4714);
    float2  _S4716 = _S4533 * _S4715;
    float _S4717 = (*dist_coeffs_40)[int(4)] * _S4715.y;
    float _S4718 = (*dist_coeffs_40)[int(5)] * _S4715.x;
    float _S4719 = _S4716.x + _S4716.y;
    float _S4720 = r2_76 * _S4719;
    float _S4721 = r2_76 * _S4720;
    float _S4722 = (*dist_coeffs_40)[int(7)] * _S4715.y + _S4717 + (*dist_coeffs_40)[int(6)] * _S4715.x + _S4718 + _S4537 * _S4719 + _S4536 * _S4720 + _S4535 * _S4721 + (*dist_coeffs_40)[int(3)] * (r2_76 * _S4721);
    float _S4723 = v_76 * _S4722;
    float _S4724 = u_76 * _S4722;
    float2  _S4725 = (_S4538 * _S4715 + make_float2 (_S4506 * (v_76 * _S4715.y) + _S4540 * _S4718 + 2.0f * (u_76 * _S4718) + _S4503 * (v_76 * _S4715.x) + _S4724 + _S4724, _S4542 * _S4717 + 2.0f * (v_76 * _S4717) + _S4541 * _S4715.y + _S4539 * _S4715.x + _S4723 + _S4723)) / _S4534;
    float2  _S4726 = _S4530 * - _S4725;
    float2  _S4727 = _S4532 * _S4725;
    float2  _S4728 = - _S4711 + _S4712;
    float _S4729 = fx_34 * (_S4698.differential_0 + _S4702.differential_0 + _S4728.x);
    float2  _S4730 = make_float2 (_S4729, fy_34 * (_S4690.differential_0 + _S4694.differential_0 + _S4728.y)) + make_float2 ((*dist_coeffs_40)[int(8)] * _S4729, (*dist_coeffs_40)[int(9)] * _S4729);
    float2  _S4731 = _S4516 * _S4730;
    float _S4732 = (*dist_coeffs_40)[int(4)] * _S4730.y;
    float _S4733 = (*dist_coeffs_40)[int(5)] * _S4730.x;
    float _S4734 = _S4731.x + _S4731.y;
    float _S4735 = r2_75 * _S4734;
    float _S4736 = r2_75 * _S4735;
    float _S4737 = (*dist_coeffs_40)[int(7)] * _S4730.y + _S4732 + (*dist_coeffs_40)[int(6)] * _S4730.x + _S4733 + _S4520 * _S4734 + _S4519 * _S4735 + _S4518 * _S4736 + (*dist_coeffs_40)[int(3)] * (r2_75 * _S4736);
    float _S4738 = v_75 * _S4737;
    float _S4739 = u_75 * _S4737;
    float2  _S4740 = (_S4521 * _S4730 + make_float2 (_S4506 * (v_75 * _S4730.y) + _S4523 * _S4733 + 2.0f * (u_75 * _S4733) + _S4503 * (v_75 * _S4730.x) + _S4739 + _S4739, _S4525 * _S4732 + 2.0f * (v_75 * _S4732) + _S4524 * _S4730.y + _S4522 * _S4730.x + _S4738 + _S4738)) / _S4517;
    float2  _S4741 = _S4513 * - _S4740;
    float2  _S4742 = _S4515 * _S4740;
    float _S4743 = _S4741.x + _S4741.y;
    float2  _S4744 = _S4706.differential_0 + - _S4712;
    float _S4745 = fx_34 * (_S4697.differential_0 + _S4701.differential_0 + _S4744.x);
    float2  _S4746 = make_float2 (_S4745, fy_34 * (_S4689.differential_0 + _S4693.differential_0 + _S4744.y)) + make_float2 ((*dist_coeffs_40)[int(8)] * _S4745, (*dist_coeffs_40)[int(9)] * _S4745);
    float2  _S4747 = _S4497 * _S4746;
    float _S4748 = (*dist_coeffs_40)[int(4)] * _S4746.y;
    float _S4749 = (*dist_coeffs_40)[int(5)] * _S4746.x;
    float _S4750 = _S4747.x + _S4747.y;
    float _S4751 = r2_74 * _S4750;
    float _S4752 = r2_74 * _S4751;
    float _S4753 = (*dist_coeffs_40)[int(7)] * _S4746.y + _S4748 + (*dist_coeffs_40)[int(6)] * _S4746.x + _S4749 + _S4501 * _S4750 + _S4500 * _S4751 + _S4499 * _S4752 + (*dist_coeffs_40)[int(3)] * (r2_74 * _S4752);
    float _S4754 = v_74 * _S4753;
    float _S4755 = u_74 * _S4753;
    float2  _S4756 = (_S4502 * _S4746 + make_float2 (_S4506 * (v_74 * _S4746.y) + _S4505 * _S4749 + 2.0f * (u_74 * _S4749) + _S4503 * (v_74 * _S4746.x) + _S4755 + _S4755, _S4508 * _S4748 + 2.0f * (v_74 * _S4748) + _S4507 * _S4746.y + _S4504 * _S4746.x + _S4754 + _S4754)) / _S4498;
    float2  _S4757 = _S4494 * - _S4756;
    float2  _S4758 = _S4496 * _S4756;
    float _S4759 = _S4757.x + _S4757.y;
    float3  _S4760 = _S4606.differential_0 + _S4686.differential_0 + make_float3 (_S4727.x, _S4727.y, _S4726.x + _S4726.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4761;
    (&_S4761)->primal_0 = R_29;
    (&_S4761)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4762;
    (&_S4762)->primal_0 = vert2_7;
    (&_S4762)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4761, &_S4762, _S4760);
    float3  _S4763 = _S4605.differential_0 + _S4686.differential_0 + make_float3 (_S4742.x, _S4742.y, _S4743);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4764;
    (&_S4764)->primal_0 = R_29;
    (&_S4764)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4765;
    (&_S4765)->primal_0 = vert1_7;
    (&_S4765)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4764, &_S4765, _S4763);
    float3  _S4766 = _S4607 + _S4608 + _S4686.differential_0 + make_float3 (_S4758.x, _S4758.y, _S4759);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4767;
    (&_S4767)->primal_0 = R_29;
    (&_S4767)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4768;
    (&_S4768)->primal_0 = vert0_7;
    (&_S4768)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4767, &_S4768, _S4766);
    float3  _S4769 = (*v_verts_0)[int(2)] + _S4762.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4770;
    (&_S4770)->primal_0 = _S4488;
    (&_S4770)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4771;
    (&_S4771)->primal_0 = _S4493;
    (&_S4771)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4770, &_S4771, _S4769);
    float _S4772 = - _S4771.differential_0.y;
    float _S4773 = _S4492 * _S4771.differential_0.x;
    float _S4774 = - (_S4484 * _S4771.differential_0.x);
    float3  _S4775 = (*v_verts_0)[int(1)] + _S4765.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4776;
    (&_S4776)->primal_0 = _S4488;
    (&_S4776)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4777;
    (&_S4777)->primal_0 = _S4491;
    (&_S4777)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4776, &_S4777, _S4775);
    float _S4778 = _S4484 * _S4777.differential_0.x;
    float _S4779 = _S4490 * _S4777.differential_0.x;
    float3  _S4780 = (*v_verts_0)[int(0)] + _S4768.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4781;
    (&_S4781)->primal_0 = _S4488;
    (&_S4781)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4782;
    (&_S4782)->primal_0 = _S4489;
    (&_S4782)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4781, &_S4782, _S4780);
    Matrix<float, 3, 3>  _S4783 = transpose_0(_S4770.differential_0 + _S4776.differential_0 + _S4781.differential_0);
    float _S4784 = 2.0f * - _S4783.rows[int(2)].z;
    float _S4785 = 2.0f * _S4783.rows[int(2)].y;
    float _S4786 = 2.0f * _S4783.rows[int(2)].x;
    float _S4787 = 2.0f * _S4783.rows[int(1)].z;
    float _S4788 = 2.0f * - _S4783.rows[int(1)].y;
    float _S4789 = 2.0f * _S4783.rows[int(1)].x;
    float _S4790 = 2.0f * _S4783.rows[int(0)].z;
    float _S4791 = 2.0f * _S4783.rows[int(0)].y;
    float _S4792 = 2.0f * - _S4783.rows[int(0)].x;
    float _S4793 = - _S4789 + _S4791;
    float _S4794 = _S4786 + - _S4790;
    float _S4795 = - _S4785 + _S4787;
    float _S4796 = _S4785 + _S4787;
    float _S4797 = _S4786 + _S4790;
    float _S4798 = _S4789 + _S4791;
    float _S4799 = quat_31.w * (_S4788 + _S4792);
    float _S4800 = quat_31.z * (_S4784 + _S4792);
    float _S4801 = quat_31.y * (_S4784 + _S4788);
    float _S4802 = quat_31.x * _S4793 + quat_31.z * _S4796 + quat_31.y * _S4797 + _S4799 + _S4799;
    float _S4803 = quat_31.x * _S4794 + quat_31.w * _S4796 + quat_31.y * _S4798 + _S4800 + _S4800;
    float _S4804 = quat_31.x * _S4795 + quat_31.w * _S4797 + quat_31.z * _S4798 + _S4801 + _S4801;
    float _S4805 = quat_31.w * _S4793 + quat_31.z * _S4794 + quat_31.y * _S4795;
    float _S4806 = _S4774 + _S4778;
    float _S4807 = 0.5f * - _S4806;
    float _S4808 = _S4772 + _S4777.differential_0.y;
    DiffPair_float_0 _S4809;
    (&_S4809)->primal_0 = _S4485;
    (&_S4809)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4809, _S4808);
    float _S4810 = _S4807 + _S4809.differential_0;
    float _S4811 = _S4773 + _S4779 + _S4782.differential_0.x;
    DiffPair_float_0 _S4812;
    (&_S4812)->primal_0 = _S4483;
    (&_S4812)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4812, _S4811);
    float _S4813 = _S4807 + _S4812.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4814;
    (&_S4814)->primal_0 = R_29;
    (&_S4814)->differential_0 = _S4680;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4815;
    (&_S4815)->primal_0 = mean_30;
    (&_S4815)->differential_0 = _S4600;
    s_bwd_prop_mul_1(&_S4814, &_S4815, _S4602.differential_0);
    float3  _S4816 = _S4682.differential_0 + _S4760 + _S4763 + _S4766 + _S4602.differential_0;
    Matrix<float, 3, 3>  _S4817 = _S4683 + _S4761.differential_0 + _S4764.differential_0 + _S4767.differential_0 + _S4814.differential_0;
    FixedArray<float3 , 2>  _S4818;
    _S4818[int(0)] = _S4600;
    _S4818[int(1)] = _S4600;
    _S4818[int(1)] = _S4614;
    _S4818[int(0)] = _S4619;
    FixedArray<float3 , 16>  _S4819;
    _S4819[int(0)] = _S4600;
    _S4819[int(1)] = _S4600;
    _S4819[int(2)] = _S4600;
    _S4819[int(3)] = _S4600;
    _S4819[int(4)] = _S4600;
    _S4819[int(5)] = _S4600;
    _S4819[int(6)] = _S4600;
    _S4819[int(7)] = _S4600;
    _S4819[int(8)] = _S4600;
    _S4819[int(9)] = _S4600;
    _S4819[int(10)] = _S4600;
    _S4819[int(11)] = _S4600;
    _S4819[int(12)] = _S4600;
    _S4819[int(13)] = _S4600;
    _S4819[int(14)] = _S4600;
    _S4819[int(15)] = _S4600;
    _S4819[int(15)] = _S4621;
    _S4819[int(14)] = _S4623;
    _S4819[int(13)] = _S4625;
    _S4819[int(12)] = _S4627;
    _S4819[int(11)] = _S4629;
    _S4819[int(10)] = _S4631;
    _S4819[int(9)] = _S4633;
    _S4819[int(8)] = _S4641;
    _S4819[int(7)] = _S4643;
    _S4819[int(6)] = _S4645;
    _S4819[int(5)] = _S4647;
    _S4819[int(4)] = _S4649;
    _S4819[int(3)] = _S4660;
    _S4819[int(2)] = _S4662;
    _S4819[int(1)] = _S4664;
    _S4819[int(0)] = _S4677;
    float2  _S4820 = make_float2 (0.0f, _S4704);
    float3  _S4821 = make_float3 (_S4813, _S4810, _S4806);
    float4  _S4822 = make_float4 (0.0f);
    *&((&_S4822)->w) = _S4802;
    *&((&_S4822)->z) = _S4803;
    *&((&_S4822)->y) = _S4804;
    *&((&_S4822)->x) = _S4805;
    *v_mean_9 = _S4678 + _S4769 + _S4775 + _S4780 + _S4815.differential_0;
    *v_quat_8 = _S4822;
    *v_scale_8 = _S4821;
    *v_hardness_4 = _S4820;
    *v_sh_coeffs_7 = _S4819;
    *v_ch_coeffs_2 = _S4818;
    *v_R_8 = _S4817;
    *v_t_8 = _S4816;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_31, float4  quat_32, float3  scale_31, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_30, float3  t_29, float fx_35, float fy_35, float cx_30, float cy_30, FixedArray<float, 10>  * dist_coeffs_41, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_1(R_30, mean_31) + t_29;
    float _S4823 = scale_31.x;
    float _S4824 = s_primal_ctx_exp_1(_S4823);
    float _S4825 = scale_31.y;
    float _S4826 = s_primal_ctx_exp_1(_S4825);
    float sz_12 = scale_31.z - 0.5f * (_S4823 + _S4825);
    float _S4827 = quat_32.y;
    float x2_32 = _S4827 * _S4827;
    float y2_32 = quat_32.z * quat_32.z;
    float z2_58 = quat_32.w * quat_32.w;
    float xy_32 = quat_32.y * quat_32.z;
    float xz_32 = quat_32.y * quat_32.w;
    float yz_32 = quat_32.z * quat_32.w;
    float wx_32 = quat_32.x * quat_32.y;
    float wy_32 = quat_32.x * quat_32.z;
    float wz_32 = quat_32.x * quat_32.w;
    Matrix<float, 3, 3>  _S4828 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_32 + z2_58), 2.0f * (xy_32 + wz_32), 2.0f * (xz_32 - wy_32), 2.0f * (xy_32 - wz_32), 1.0f - 2.0f * (x2_32 + z2_58), 2.0f * (yz_32 + wx_32), 2.0f * (xz_32 + wy_32), 2.0f * (yz_32 - wx_32), 1.0f - 2.0f * (x2_32 + y2_32)));
    float3  _S4829 = make_float3 (_S4824, 0.0f, 0.0f);
    float3  vert0_8 = s_primal_ctx_mul_1(_S4828, _S4829) + mean_31;
    float _S4830 = -0.5f + sz_12;
    float3  _S4831 = make_float3 (_S4824 * _S4830, _S4826, 0.0f);
    float3  vert1_8 = s_primal_ctx_mul_1(_S4828, _S4831) + mean_31;
    float _S4832 = -0.5f - sz_12;
    float3  _S4833 = make_float3 (_S4824 * _S4832, - _S4826, 0.0f);
    float3  vert2_8 = s_primal_ctx_mul_1(_S4828, _S4833) + mean_31;
    float3  vert0_c_11 = s_primal_ctx_mul_1(R_30, vert0_8) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_1(R_30, vert1_8) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_1(R_30, vert2_8) + t_29;
    float2  _S4834 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4835 = length_0(_S4834);
    float _S4836 = vert0_c_11.z;
    float _S4837 = s_primal_ctx_atan2_0(_S4835, _S4836);
    bool _S4838 = _S4837 < 0.00100000004749745f;
    float k_12;
    float _S4839;
    float _S4840;
    float _S4841;
    if(_S4838)
    {
        float _S4842 = 1.0f - _S4837 * _S4837 / 3.0f;
        float _S4843 = _S4836 * _S4836;
        k_12 = _S4842 / _S4836;
        _S4839 = 0.0f;
        _S4840 = _S4843;
        _S4841 = _S4842;
    }
    else
    {
        float _S4844 = _S4835 * _S4835;
        k_12 = _S4837 / _S4835;
        _S4839 = _S4844;
        _S4840 = 0.0f;
        _S4841 = 0.0f;
    }
    float2  _S4845 = make_float2 (k_12);
    float2  _S4846 = _S4834 * make_float2 (k_12);
    float u_77 = _S4846.x;
    float v_77 = _S4846.y;
    float r2_77 = u_77 * u_77 + v_77 * v_77;
    float _S4847 = (*dist_coeffs_41)[int(2)] + r2_77 * (*dist_coeffs_41)[int(3)];
    float _S4848 = (*dist_coeffs_41)[int(1)] + r2_77 * _S4847;
    float _S4849 = (*dist_coeffs_41)[int(0)] + r2_77 * _S4848;
    float radial_11 = 1.0f + r2_77 * _S4849;
    float2  _S4850 = make_float2 (radial_11);
    float _S4851 = 2.0f * (*dist_coeffs_41)[int(4)];
    float _S4852 = _S4851 * u_77;
    float _S4853 = 2.0f * u_77;
    float _S4854 = 2.0f * (*dist_coeffs_41)[int(5)];
    float _S4855 = _S4854 * u_77;
    float _S4856 = 2.0f * v_77;
    float2  _S4857 = _S4846 * make_float2 (radial_11) + make_float2 (_S4852 * v_77 + (*dist_coeffs_41)[int(5)] * (r2_77 + _S4853 * u_77) + (*dist_coeffs_41)[int(6)] * r2_77, _S4855 * v_77 + (*dist_coeffs_41)[int(4)] * (r2_77 + _S4856 * v_77) + (*dist_coeffs_41)[int(7)] * r2_77);
    float2  _S4858 = _S4857 + make_float2 ((*dist_coeffs_41)[int(8)] * _S4857.x + (*dist_coeffs_41)[int(9)] * _S4857.y, 0.0f);
    float _S4859 = fx_35 * _S4858.x + cx_30;
    float _S4860 = fy_35 * _S4858.y + cy_30;
    float2  uv0_11 = make_float2 (_S4859, _S4860);
    float2  _S4861 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4862 = length_0(_S4861);
    float _S4863 = vert1_c_11.z;
    float _S4864 = s_primal_ctx_atan2_0(_S4862, _S4863);
    bool _S4865 = _S4864 < 0.00100000004749745f;
    float _S4866;
    float _S4867;
    float _S4868;
    if(_S4865)
    {
        float _S4869 = 1.0f - _S4864 * _S4864 / 3.0f;
        float _S4870 = _S4863 * _S4863;
        k_12 = _S4869 / _S4863;
        _S4866 = 0.0f;
        _S4867 = _S4870;
        _S4868 = _S4869;
    }
    else
    {
        float _S4871 = _S4862 * _S4862;
        k_12 = _S4864 / _S4862;
        _S4866 = _S4871;
        _S4867 = 0.0f;
        _S4868 = 0.0f;
    }
    float2  _S4872 = make_float2 (k_12);
    float2  _S4873 = _S4861 * make_float2 (k_12);
    float u_78 = _S4873.x;
    float v_78 = _S4873.y;
    float r2_78 = u_78 * u_78 + v_78 * v_78;
    float _S4874 = (*dist_coeffs_41)[int(2)] + r2_78 * (*dist_coeffs_41)[int(3)];
    float _S4875 = (*dist_coeffs_41)[int(1)] + r2_78 * _S4874;
    float _S4876 = (*dist_coeffs_41)[int(0)] + r2_78 * _S4875;
    float radial_12 = 1.0f + r2_78 * _S4876;
    float2  _S4877 = make_float2 (radial_12);
    float _S4878 = _S4851 * u_78;
    float _S4879 = 2.0f * u_78;
    float _S4880 = _S4854 * u_78;
    float _S4881 = 2.0f * v_78;
    float2  _S4882 = _S4873 * make_float2 (radial_12) + make_float2 (_S4878 * v_78 + (*dist_coeffs_41)[int(5)] * (r2_78 + _S4879 * u_78) + (*dist_coeffs_41)[int(6)] * r2_78, _S4880 * v_78 + (*dist_coeffs_41)[int(4)] * (r2_78 + _S4881 * v_78) + (*dist_coeffs_41)[int(7)] * r2_78);
    float2  _S4883 = _S4882 + make_float2 ((*dist_coeffs_41)[int(8)] * _S4882.x + (*dist_coeffs_41)[int(9)] * _S4882.y, 0.0f);
    float _S4884 = fx_35 * _S4883.x + cx_30;
    float _S4885 = fy_35 * _S4883.y + cy_30;
    float2  uv1_11 = make_float2 (_S4884, _S4885);
    float2  _S4886 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4887 = length_0(_S4886);
    float _S4888 = vert2_c_11.z;
    float _S4889 = s_primal_ctx_atan2_0(_S4887, _S4888);
    bool _S4890 = _S4889 < 0.00100000004749745f;
    float _S4891;
    float _S4892;
    float _S4893;
    if(_S4890)
    {
        float _S4894 = 1.0f - _S4889 * _S4889 / 3.0f;
        float _S4895 = _S4888 * _S4888;
        k_12 = _S4894 / _S4888;
        _S4891 = 0.0f;
        _S4892 = _S4895;
        _S4893 = _S4894;
    }
    else
    {
        float _S4896 = _S4887 * _S4887;
        k_12 = _S4889 / _S4887;
        _S4891 = _S4896;
        _S4892 = 0.0f;
        _S4893 = 0.0f;
    }
    float2  _S4897 = make_float2 (k_12);
    float2  _S4898 = _S4886 * make_float2 (k_12);
    float u_79 = _S4898.x;
    float v_79 = _S4898.y;
    float r2_79 = u_79 * u_79 + v_79 * v_79;
    float _S4899 = (*dist_coeffs_41)[int(2)] + r2_79 * (*dist_coeffs_41)[int(3)];
    float _S4900 = (*dist_coeffs_41)[int(1)] + r2_79 * _S4899;
    float _S4901 = (*dist_coeffs_41)[int(0)] + r2_79 * _S4900;
    float radial_13 = 1.0f + r2_79 * _S4901;
    float2  _S4902 = make_float2 (radial_13);
    float _S4903 = _S4851 * u_79;
    float _S4904 = 2.0f * u_79;
    float _S4905 = _S4854 * u_79;
    float _S4906 = 2.0f * v_79;
    float2  _S4907 = _S4898 * make_float2 (radial_13) + make_float2 (_S4903 * v_79 + (*dist_coeffs_41)[int(5)] * (r2_79 + _S4904 * u_79) + (*dist_coeffs_41)[int(6)] * r2_79, _S4905 * v_79 + (*dist_coeffs_41)[int(4)] * (r2_79 + _S4906 * v_79) + (*dist_coeffs_41)[int(7)] * r2_79);
    float2  _S4908 = _S4907 + make_float2 ((*dist_coeffs_41)[int(8)] * _S4907.x + (*dist_coeffs_41)[int(9)] * _S4907.y, 0.0f);
    float _S4909 = fx_35 * _S4908.x + cx_30;
    float _S4910 = fy_35 * _S4908.y + cy_30;
    float2  uv2_11 = make_float2 (_S4909, _S4910);
    float2  e0_15 = uv1_11 - uv0_11;
    float2  e1_15 = uv2_11 - uv1_11;
    float2  e2_7 = uv0_11 - uv2_11;
    float _S4911 = e0_15.x;
    float _S4912 = e1_15.y;
    float _S4913 = e0_15.y;
    float _S4914 = e1_15.x;
    float _S4915 = _S4911 * _S4912 - _S4913 * _S4914;
    float _S4916 = 1.0f - hardness_15.y;
    float _S4917 = -1.0f / _S4916;
    float _S4918 = _S4916 * _S4916;
    float _S4919 = s_primal_ctx_max_0(_S4859, _S4884);
    float _S4920 = s_primal_ctx_min_0(_S4859, _S4884);
    float _S4921 = s_primal_ctx_max_0(_S4860, _S4885);
    float _S4922 = s_primal_ctx_min_0(_S4860, _S4885);
    float3  _S4923 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S4924 = length_1(_S4923) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4925 = transpose_0(R_30);
    float3  _S4926 = mean_31 - - s_primal_ctx_mul_1(_S4925, t_29);
    float _S4927 = _S4926.x;
    float _S4928 = _S4926.y;
    float _S4929 = _S4926.z;
    float _S4930 = _S4927 * _S4927 + _S4928 * _S4928 + _S4929 * _S4929;
    float _S4931 = s_primal_ctx_sqrt_0(_S4930);
    float x_62 = _S4927 / _S4931;
    float3  _S4932 = make_float3 (x_62);
    float _S4933 = _S4931 * _S4931;
    float y_29 = _S4928 / _S4931;
    float z_26 = _S4929 / _S4931;
    float3  _S4934 = make_float3 (z_26);
    float _S4935 = - y_29;
    float3  _S4936 = make_float3 (_S4935);
    float z2_59 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_62 * x_62 - y_29 * y_29;
    float _S4937 = 2.0f * x_62;
    float fS1_26 = _S4937 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_59 - 0.31539157032966614f;
    float3  _S4938 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_62;
    float3  _S4939 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4940 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4941 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4942 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_59 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4943 = 1.86588168144226074f * z2_59 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4943;
    float3  _S4944 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_62;
    float3  _S4945 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4946 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4947 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4948 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_62 * fC1_26 - y_29 * fS1_26);
    float3  _S4949 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_62 * fS1_26 + y_29 * fC1_26);
    float3  _S4950 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4935) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_62) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4951 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4952 = make_float3 (0.0f);
    float3  _S4953 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4954 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4955 = make_float3 (_S4954);
    float3  _S4956 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4954);
    float3  _S4957 = _S4953 + _S4956 + make_float3 (0.5f);
    float3  _S4958 = _S4953 - _S4956 + make_float3 (0.5f);
    float3  _S4959 = vert1_c_11 - vert0_c_11;
    float3  _S4960 = vert2_c_11 - vert0_c_11;
    float3  _S4961 = s_primal_ctx_cross_0(_S4959, _S4960);
    float3  _S4962 = normalize_0(_S4961);
    float3  _S4963 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4962, mean_c_26)))))) * v_normal_3;
    float3  _S4964 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4965;
    (&_S4965)->primal_0 = _S4962;
    (&_S4965)->differential_0 = _S4964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4966;
    (&_S4966)->primal_0 = mean_c_26;
    (&_S4966)->differential_0 = _S4964;
    s_bwd_prop_dot_0(&_S4965, &_S4966, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4967 = _S4966;
    float3  _S4968 = _S4963 + _S4965.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4969;
    (&_S4969)->primal_0 = _S4961;
    (&_S4969)->differential_0 = _S4964;
    s_bwd_normalize_impl_0(&_S4969, _S4968);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4970;
    (&_S4970)->primal_0 = _S4959;
    (&_S4970)->differential_0 = _S4964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4971;
    (&_S4971)->primal_0 = _S4960;
    (&_S4971)->differential_0 = _S4964;
    s_bwd_prop_cross_0(&_S4970, &_S4971, _S4969.differential_0);
    float3  _S4972 = - _S4971.differential_0;
    float3  _S4973 = - _S4970.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4974;
    (&_S4974)->primal_0 = _S4958;
    (&_S4974)->differential_0 = _S4964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4975;
    (&_S4975)->primal_0 = _S4952;
    (&_S4975)->differential_0 = _S4964;
    s_bwd_prop_max_0(&_S4974, &_S4975, (*v_rgbs_1)[int(2)]);
    float3  _S4976 = - _S4974.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4977;
    (&_S4977)->primal_0 = _S4957;
    (&_S4977)->differential_0 = _S4964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4978;
    (&_S4978)->primal_0 = _S4952;
    (&_S4978)->differential_0 = _S4964;
    s_bwd_prop_max_0(&_S4977, &_S4978, (*v_rgbs_1)[int(1)]);
    float3  _S4979 = _S4955 * (_S4976 + _S4977.differential_0);
    float3  _S4980 = _S4974.differential_0 + _S4977.differential_0;
    float3  _S4981 = make_float3 (0.5f) * - _S4980;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4982;
    (&_S4982)->primal_0 = _S4951;
    (&_S4982)->differential_0 = _S4964;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4983;
    (&_S4983)->primal_0 = _S4952;
    (&_S4983)->differential_0 = _S4964;
    s_bwd_prop_max_0(&_S4982, &_S4983, (*v_rgbs_1)[int(0)]);
    float3  _S4984 = _S4981 + _S4982.differential_0;
    float3  _S4985 = _S4980 + _S4982.differential_0;
    float3  _S4986 = _S4949 * _S4985;
    float3  _S4987 = (*sh_coeffs_26)[int(15)] * _S4985;
    float3  _S4988 = _S4947 * _S4985;
    float3  _S4989 = (*sh_coeffs_26)[int(14)] * _S4985;
    float3  _S4990 = _S4945 * _S4985;
    float3  _S4991 = (*sh_coeffs_26)[int(13)] * _S4985;
    float3  _S4992 = _S4944 * _S4985;
    float3  _S4993 = (*sh_coeffs_26)[int(12)] * _S4985;
    float3  _S4994 = _S4946 * _S4985;
    float3  _S4995 = (*sh_coeffs_26)[int(11)] * _S4985;
    float3  _S4996 = _S4948 * _S4985;
    float3  _S4997 = (*sh_coeffs_26)[int(10)] * _S4985;
    float3  _S4998 = _S4950 * _S4985;
    float3  _S4999 = (*sh_coeffs_26)[int(9)] * _S4985;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4999.x + _S4999.y + _S4999.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4987.x + _S4987.y + _S4987.z);
    float _S5000 = _S4997.x + _S4997.y + _S4997.z;
    float _S5001 = _S4989.x + _S4989.y + _S4989.z;
    float _S5002 = _S4995.x + _S4995.y + _S4995.z;
    float _S5003 = _S4991.x + _S4991.y + _S4991.z;
    float _S5004 = _S4993.x + _S4993.y + _S4993.z;
    float _S5005 = - s_diff_fC2_T_8;
    float3  _S5006 = _S4941 * _S4985;
    float3  _S5007 = (*sh_coeffs_26)[int(8)] * _S4985;
    float3  _S5008 = _S4939 * _S4985;
    float3  _S5009 = (*sh_coeffs_26)[int(7)] * _S4985;
    float3  _S5010 = _S4938 * _S4985;
    float3  _S5011 = (*sh_coeffs_26)[int(6)] * _S4985;
    float3  _S5012 = _S4940 * _S4985;
    float3  _S5013 = (*sh_coeffs_26)[int(5)] * _S4985;
    float3  _S5014 = _S4942 * _S4985;
    float3  _S5015 = (*sh_coeffs_26)[int(4)] * _S4985;
    float _S5016 = _S5013.x + _S5013.y + _S5013.z;
    float _S5017 = _S5009.x + _S5009.y + _S5009.z;
    float _S5018 = fTmp1B_26 * _S5000 + x_62 * s_diff_fS2_T_8 + y_29 * _S5005 + 0.54627424478530884f * (_S5015.x + _S5015.y + _S5015.z);
    float _S5019 = fTmp1B_26 * _S5001 + y_29 * s_diff_fS2_T_8 + x_62 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S5007.x + _S5007.y + _S5007.z);
    float _S5020 = y_29 * - _S5019;
    float _S5021 = x_62 * _S5019;
    float _S5022 = z_26 * (1.86588168144226074f * (z_26 * _S5004) + -2.28522896766662598f * (y_29 * _S5002 + x_62 * _S5003) + 0.94617468118667603f * (_S5011.x + _S5011.y + _S5011.z));
    float3  _S5023 = make_float3 (0.48860251903533936f) * _S4985;
    float3  _S5024 = - _S5023;
    float3  _S5025 = _S4932 * _S5024;
    float3  _S5026 = (*sh_coeffs_26)[int(3)] * _S5024;
    float3  _S5027 = _S4934 * _S5023;
    float3  _S5028 = (*sh_coeffs_26)[int(2)] * _S5023;
    float3  _S5029 = _S4936 * _S5023;
    float3  _S5030 = (*sh_coeffs_26)[int(1)] * _S5023;
    float _S5031 = (_S4943 * _S5004 + 1.44530570507049561f * (fS1_26 * _S5000 + fC1_26 * _S5001) + -1.09254848957061768f * (y_29 * _S5016 + x_62 * _S5017) + _S5022 + _S5022 + _S5028.x + _S5028.y + _S5028.z) / _S4933;
    float _S5032 = _S4931 * _S5031;
    float _S5033 = (fTmp0C_26 * _S5002 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S5005 + fTmp0B_26 * _S5016 + _S4937 * _S5018 + _S5020 + _S5020 + - (_S5030.x + _S5030.y + _S5030.z)) / _S4933;
    float _S5034 = _S4931 * _S5033;
    float _S5035 = (fTmp0C_26 * _S5003 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S5017 + 2.0f * (y_29 * _S5018) + _S5021 + _S5021 + _S5026.x + _S5026.y + _S5026.z) / _S4933;
    float _S5036 = _S4931 * _S5035;
    float _S5037 = _S4929 * - _S5031 + _S4928 * - _S5033 + _S4927 * - _S5035;
    DiffPair_float_0 _S5038;
    (&_S5038)->primal_0 = _S4930;
    (&_S5038)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5038, _S5037);
    float _S5039 = _S4929 * _S5038.differential_0;
    float _S5040 = _S4928 * _S5038.differential_0;
    float _S5041 = _S4927 * _S5038.differential_0;
    float3  _S5042 = make_float3 (0.282094806432724f) * _S4985;
    float3  _S5043 = make_float3 (_S5036 + _S5041 + _S5041, _S5034 + _S5040 + _S5040, _S5032 + _S5039 + _S5039);
    float3  _S5044 = - - _S5043;
    Matrix<float, 3, 3>  _S5045 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5046;
    (&_S5046)->primal_0 = _S4925;
    (&_S5046)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5047;
    (&_S5047)->primal_0 = t_29;
    (&_S5047)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5046, &_S5047, _S5044);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5048 = _S5047;
    Matrix<float, 3, 3>  _S5049 = transpose_0(_S5046.differential_0);
    DiffPair_float_0 _S5050;
    (&_S5050)->primal_0 = _S4924;
    (&_S5050)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5050, v_depth_10);
    float _S5051 = 0.3333333432674408f * _S5050.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5052;
    (&_S5052)->primal_0 = _S4923;
    (&_S5052)->differential_0 = _S4964;
    s_bwd_length_impl_0(&_S5052, _S5051);
    DiffPair_float_0 _S5053;
    (&_S5053)->primal_0 = _S4922;
    (&_S5053)->differential_0 = 0.0f;
    DiffPair_float_0 _S5054;
    (&_S5054)->primal_0 = _S4910;
    (&_S5054)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5053, &_S5054, 0.0f);
    DiffPair_float_0 _S5055;
    (&_S5055)->primal_0 = _S4860;
    (&_S5055)->differential_0 = 0.0f;
    DiffPair_float_0 _S5056;
    (&_S5056)->primal_0 = _S4885;
    (&_S5056)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5055, &_S5056, _S5053.differential_0);
    DiffPair_float_0 _S5057;
    (&_S5057)->primal_0 = _S4921;
    (&_S5057)->differential_0 = 0.0f;
    DiffPair_float_0 _S5058;
    (&_S5058)->primal_0 = _S4910;
    (&_S5058)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5057, &_S5058, 0.0f);
    DiffPair_float_0 _S5059;
    (&_S5059)->primal_0 = _S4860;
    (&_S5059)->differential_0 = 0.0f;
    DiffPair_float_0 _S5060;
    (&_S5060)->primal_0 = _S4885;
    (&_S5060)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5059, &_S5060, _S5057.differential_0);
    DiffPair_float_0 _S5061;
    (&_S5061)->primal_0 = _S4920;
    (&_S5061)->differential_0 = 0.0f;
    DiffPair_float_0 _S5062;
    (&_S5062)->primal_0 = _S4909;
    (&_S5062)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5061, &_S5062, 0.0f);
    DiffPair_float_0 _S5063;
    (&_S5063)->primal_0 = _S4859;
    (&_S5063)->differential_0 = 0.0f;
    DiffPair_float_0 _S5064;
    (&_S5064)->primal_0 = _S4884;
    (&_S5064)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5063, &_S5064, _S5061.differential_0);
    DiffPair_float_0 _S5065;
    (&_S5065)->primal_0 = _S4919;
    (&_S5065)->differential_0 = 0.0f;
    DiffPair_float_0 _S5066;
    (&_S5066)->primal_0 = _S4909;
    (&_S5066)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5065, &_S5066, 0.0f);
    DiffPair_float_0 _S5067;
    (&_S5067)->primal_0 = _S4859;
    (&_S5067)->differential_0 = 0.0f;
    DiffPair_float_0 _S5068;
    (&_S5068)->primal_0 = _S4884;
    (&_S5068)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5067, &_S5068, _S5065.differential_0);
    DiffPair_float_0 _S5069;
    (&_S5069)->primal_0 = _S4917;
    (&_S5069)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S5069, 0.0f);
    float _S5070 = - (-1.0f * - (_S5069.differential_0 / _S4918));
    float2  _S5071 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5072;
    (&_S5072)->primal_0 = e2_7;
    (&_S5072)->differential_0 = _S5071;
    s_bwd_length_impl_1(&_S5072, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5073;
    (&_S5073)->primal_0 = e1_15;
    (&_S5073)->differential_0 = _S5071;
    s_bwd_length_impl_1(&_S5073, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5074;
    (&_S5074)->primal_0 = e0_15;
    (&_S5074)->differential_0 = _S5071;
    s_bwd_length_impl_1(&_S5074, -0.0f);
    DiffPair_float_0 _S5075;
    (&_S5075)->primal_0 = _S4915;
    (&_S5075)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S5075, 0.0f);
    float _S5076 = - _S5075.differential_0;
    float2  _S5077 = _S5073.differential_0 + make_float2 (_S4913 * _S5076, _S4911 * _S5075.differential_0);
    float2  _S5078 = - _S5077;
    float2  _S5079 = _S5074.differential_0 + make_float2 (_S4912 * _S5075.differential_0, _S4914 * _S5076);
    float2  _S5080 = - _S5079;
    float2  _S5081 = - _S5072.differential_0 + _S5077;
    float _S5082 = fx_35 * (_S5062.differential_0 + _S5066.differential_0 + _S5081.x);
    float2  _S5083 = make_float2 (_S5082, fy_35 * (_S5054.differential_0 + _S5058.differential_0 + _S5081.y)) + make_float2 ((*dist_coeffs_41)[int(8)] * _S5082, (*dist_coeffs_41)[int(9)] * _S5082);
    float2  _S5084 = _S4898 * _S5083;
    float2  _S5085 = _S4902 * _S5083;
    float _S5086 = (*dist_coeffs_41)[int(4)] * _S5083.y;
    float _S5087 = (*dist_coeffs_41)[int(5)] * _S5083.x;
    float _S5088 = _S5084.x + _S5084.y;
    float _S5089 = r2_79 * _S5088;
    float _S5090 = r2_79 * _S5089;
    float _S5091 = (*dist_coeffs_41)[int(7)] * _S5083.y + _S5086 + (*dist_coeffs_41)[int(6)] * _S5083.x + _S5087 + _S4901 * _S5088 + _S4900 * _S5089 + _S4899 * _S5090 + (*dist_coeffs_41)[int(3)] * (r2_79 * _S5090);
    float _S5092 = v_79 * _S5091;
    float _S5093 = u_79 * _S5091;
    float _S5094 = _S4906 * _S5086 + 2.0f * (v_79 * _S5086) + _S4905 * _S5083.y + _S4903 * _S5083.x + _S5092 + _S5092;
    float _S5095 = _S4854 * (v_79 * _S5083.y) + _S4904 * _S5087 + 2.0f * (u_79 * _S5087) + _S4851 * (v_79 * _S5083.x) + _S5093 + _S5093;
    float3  _S5096 = _S4970.differential_0 + _S5052.differential_0;
    float3  _S5097 = _S4972 + _S4973 + _S5052.differential_0;
    float3  _S5098 = _S4971.differential_0 + _S5052.differential_0;
    FixedArray<float3 , 2>  _S5099;
    _S5099[int(0)] = _S4964;
    _S5099[int(1)] = _S4964;
    _S5099[int(1)] = _S4979;
    _S5099[int(0)] = _S4984;
    float3  _S5100 = _S5099[int(0)];
    float3  _S5101 = _S5099[int(1)];
    FixedArray<float3 , 16>  _S5102;
    _S5102[int(0)] = _S4964;
    _S5102[int(1)] = _S4964;
    _S5102[int(2)] = _S4964;
    _S5102[int(3)] = _S4964;
    _S5102[int(4)] = _S4964;
    _S5102[int(5)] = _S4964;
    _S5102[int(6)] = _S4964;
    _S5102[int(7)] = _S4964;
    _S5102[int(8)] = _S4964;
    _S5102[int(9)] = _S4964;
    _S5102[int(10)] = _S4964;
    _S5102[int(11)] = _S4964;
    _S5102[int(12)] = _S4964;
    _S5102[int(13)] = _S4964;
    _S5102[int(14)] = _S4964;
    _S5102[int(15)] = _S4964;
    _S5102[int(7)] = _S5008;
    _S5102[int(0)] = _S5042;
    _S5102[int(1)] = _S5029;
    _S5102[int(2)] = _S5027;
    _S5102[int(3)] = _S5025;
    _S5102[int(4)] = _S5014;
    _S5102[int(5)] = _S5012;
    _S5102[int(6)] = _S5010;
    _S5102[int(15)] = _S4986;
    _S5102[int(8)] = _S5006;
    _S5102[int(9)] = _S4998;
    _S5102[int(10)] = _S4996;
    _S5102[int(11)] = _S4994;
    _S5102[int(12)] = _S4992;
    _S5102[int(13)] = _S4990;
    _S5102[int(14)] = _S4988;
    float3  _S5103 = _S5102[int(0)];
    float3  _S5104 = _S5102[int(1)];
    float3  _S5105 = _S5102[int(2)];
    float3  _S5106 = _S5102[int(3)];
    float3  _S5107 = _S5102[int(4)];
    float3  _S5108 = _S5102[int(5)];
    float3  _S5109 = _S5102[int(6)];
    float3  _S5110 = _S5102[int(7)];
    float3  _S5111 = _S5102[int(8)];
    float3  _S5112 = _S5102[int(9)];
    float3  _S5113 = _S5102[int(10)];
    float3  _S5114 = _S5102[int(11)];
    float3  _S5115 = _S5102[int(12)];
    float3  _S5116 = _S5102[int(13)];
    float3  _S5117 = _S5102[int(14)];
    float3  _S5118 = _S5102[int(15)];
    float _S5119 = _S5063.differential_0 + _S5067.differential_0;
    float2  _S5120 = _S5072.differential_0 + _S5080;
    float _S5121 = _S5056.differential_0 + _S5060.differential_0;
    float _S5122 = _S5055.differential_0 + _S5059.differential_0;
    float2  _S5123 = _S5078 + _S5079;
    float _S5124 = _S5064.differential_0 + _S5068.differential_0;
    float2  _S5125 = make_float2 (0.0f, _S5070);
    float2  _S5126 = _S5085 + make_float2 (_S5095, _S5094);
    float2  _S5127 = _S4886 * _S5126;
    float2  _S5128 = _S4897 * _S5126;
    float _S5129 = _S5127.x + _S5127.y;
    if(_S4890)
    {
        float _S5130 = _S5129 / _S4892;
        float _S5131 = _S4893 * - _S5130;
        float _S5132 = _S4889 * (0.3333333432674408f * - (_S4888 * _S5130));
        k_12 = _S5132 + _S5132;
        _S4891 = _S5131;
        _S4892 = 0.0f;
    }
    else
    {
        float _S5133 = _S5129 / _S4891;
        float _S5134 = _S4889 * - _S5133;
        k_12 = _S4887 * _S5133;
        _S4891 = 0.0f;
        _S4892 = _S5134;
    }
    DiffPair_float_0 _S5135;
    (&_S5135)->primal_0 = _S4887;
    (&_S5135)->differential_0 = 0.0f;
    DiffPair_float_0 _S5136;
    (&_S5136)->primal_0 = _S4888;
    (&_S5136)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5135, &_S5136, k_12);
    float _S5137 = _S5136.differential_0 + _S4891;
    float _S5138 = _S5135.differential_0 + _S4892;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5139;
    (&_S5139)->primal_0 = _S4886;
    (&_S5139)->differential_0 = _S5071;
    s_bwd_length_impl_1(&_S5139, _S5138);
    float2  _S5140 = _S5139.differential_0 + _S5128;
    float _S5141 = fx_35 * (_S5123.x + _S5124);
    float2  _S5142 = make_float2 (_S5141, fy_35 * (_S5123.y + _S5121)) + make_float2 ((*dist_coeffs_41)[int(8)] * _S5141, (*dist_coeffs_41)[int(9)] * _S5141);
    float2  _S5143 = _S4873 * _S5142;
    float _S5144 = (*dist_coeffs_41)[int(4)] * _S5142.y;
    float _S5145 = (*dist_coeffs_41)[int(5)] * _S5142.x;
    float _S5146 = _S5143.x + _S5143.y;
    float _S5147 = r2_78 * _S5146;
    float _S5148 = r2_78 * _S5147;
    float _S5149 = (*dist_coeffs_41)[int(7)] * _S5142.y + _S5144 + (*dist_coeffs_41)[int(6)] * _S5142.x + _S5145 + _S4876 * _S5146 + _S4875 * _S5147 + _S4874 * _S5148 + (*dist_coeffs_41)[int(3)] * (r2_78 * _S5148);
    float _S5150 = v_78 * _S5149;
    float _S5151 = u_78 * _S5149;
    float3  _S5152 = _S5098 + make_float3 (_S5140.x, _S5140.y, _S5137);
    float2  _S5153 = _S4877 * _S5142 + make_float2 (_S4854 * (v_78 * _S5142.y) + _S4879 * _S5145 + 2.0f * (u_78 * _S5145) + _S4851 * (v_78 * _S5142.x) + _S5151 + _S5151, _S4881 * _S5144 + 2.0f * (v_78 * _S5144) + _S4880 * _S5142.y + _S4878 * _S5142.x + _S5150 + _S5150);
    float2  _S5154 = _S4861 * _S5153;
    float2  _S5155 = _S4872 * _S5153;
    float _S5156 = _S5154.x + _S5154.y;
    if(_S4865)
    {
        float _S5157 = _S5156 / _S4867;
        float _S5158 = _S4868 * - _S5157;
        float _S5159 = _S4864 * (0.3333333432674408f * - (_S4863 * _S5157));
        k_12 = _S5159 + _S5159;
        _S4866 = _S5158;
        _S4867 = 0.0f;
    }
    else
    {
        float _S5160 = _S5156 / _S4866;
        float _S5161 = _S4864 * - _S5160;
        k_12 = _S4862 * _S5160;
        _S4866 = 0.0f;
        _S4867 = _S5161;
    }
    DiffPair_float_0 _S5162;
    (&_S5162)->primal_0 = _S4862;
    (&_S5162)->differential_0 = 0.0f;
    DiffPair_float_0 _S5163;
    (&_S5163)->primal_0 = _S4863;
    (&_S5163)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5162, &_S5163, k_12);
    float _S5164 = _S5163.differential_0 + _S4866;
    float _S5165 = _S5162.differential_0 + _S4867;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5166;
    (&_S5166)->primal_0 = _S4861;
    (&_S5166)->differential_0 = _S5071;
    s_bwd_length_impl_1(&_S5166, _S5165);
    float2  _S5167 = _S5166.differential_0 + _S5155;
    float _S5168 = fx_35 * (_S5120.x + _S5119);
    float2  _S5169 = make_float2 (_S5168, fy_35 * (_S5120.y + _S5122)) + make_float2 ((*dist_coeffs_41)[int(8)] * _S5168, (*dist_coeffs_41)[int(9)] * _S5168);
    float2  _S5170 = _S4846 * _S5169;
    float _S5171 = (*dist_coeffs_41)[int(4)] * _S5169.y;
    float _S5172 = (*dist_coeffs_41)[int(5)] * _S5169.x;
    float _S5173 = _S5170.x + _S5170.y;
    float _S5174 = r2_77 * _S5173;
    float _S5175 = r2_77 * _S5174;
    float _S5176 = (*dist_coeffs_41)[int(7)] * _S5169.y + _S5171 + (*dist_coeffs_41)[int(6)] * _S5169.x + _S5172 + _S4849 * _S5173 + _S4848 * _S5174 + _S4847 * _S5175 + (*dist_coeffs_41)[int(3)] * (r2_77 * _S5175);
    float _S5177 = v_77 * _S5176;
    float _S5178 = u_77 * _S5176;
    float3  _S5179 = _S5096 + make_float3 (_S5167.x, _S5167.y, _S5164);
    float2  _S5180 = _S4850 * _S5169 + make_float2 (_S4854 * (v_77 * _S5169.y) + _S4853 * _S5172 + 2.0f * (u_77 * _S5172) + _S4851 * (v_77 * _S5169.x) + _S5178 + _S5178, _S4856 * _S5171 + 2.0f * (v_77 * _S5171) + _S4855 * _S5169.y + _S4852 * _S5169.x + _S5177 + _S5177);
    float2  _S5181 = _S4834 * _S5180;
    float2  _S5182 = _S4845 * _S5180;
    float _S5183 = _S5181.x + _S5181.y;
    if(_S4838)
    {
        float _S5184 = _S5183 / _S4840;
        float _S5185 = _S4841 * - _S5184;
        float _S5186 = _S4837 * (0.3333333432674408f * - (_S4836 * _S5184));
        k_12 = _S5186 + _S5186;
        _S4839 = _S5185;
        _S4840 = 0.0f;
    }
    else
    {
        float _S5187 = _S5183 / _S4839;
        float _S5188 = _S4837 * - _S5187;
        k_12 = _S4835 * _S5187;
        _S4839 = 0.0f;
        _S4840 = _S5188;
    }
    DiffPair_float_0 _S5189;
    (&_S5189)->primal_0 = _S4835;
    (&_S5189)->differential_0 = 0.0f;
    DiffPair_float_0 _S5190;
    (&_S5190)->primal_0 = _S4836;
    (&_S5190)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5189, &_S5190, k_12);
    float _S5191 = _S5190.differential_0 + _S4839;
    float _S5192 = _S5189.differential_0 + _S4840;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5193;
    (&_S5193)->primal_0 = _S4834;
    (&_S5193)->differential_0 = _S5071;
    s_bwd_length_impl_1(&_S5193, _S5192);
    float2  _S5194 = _S5193.differential_0 + _S5182;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5195;
    (&_S5195)->primal_0 = vert2_c_11;
    (&_S5195)->differential_0 = _S4964;
    s_bwd_length_impl_0(&_S5195, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5196;
    (&_S5196)->primal_0 = vert1_c_11;
    (&_S5196)->differential_0 = _S4964;
    s_bwd_length_impl_0(&_S5196, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5197;
    (&_S5197)->primal_0 = vert0_c_11;
    (&_S5197)->differential_0 = _S4964;
    s_bwd_length_impl_0(&_S5197, 0.0f);
    float3  _S5198 = _S5195.differential_0 + _S5152;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5199;
    (&_S5199)->primal_0 = R_30;
    (&_S5199)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5200;
    (&_S5200)->primal_0 = vert2_8;
    (&_S5200)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5199, &_S5200, _S5198);
    float3  _S5201 = _S5196.differential_0 + _S5179;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5202;
    (&_S5202)->primal_0 = R_30;
    (&_S5202)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5203;
    (&_S5203)->primal_0 = vert1_8;
    (&_S5203)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5202, &_S5203, _S5201);
    float3  _S5204 = _S5197.differential_0 + _S5097 + make_float3 (_S5194.x, _S5194.y, _S5191);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5205;
    (&_S5205)->primal_0 = R_30;
    (&_S5205)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5206;
    (&_S5206)->primal_0 = vert0_8;
    (&_S5206)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5205, &_S5206, _S5204);
    float3  _S5207 = _S5200.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5208;
    (&_S5208)->primal_0 = _S4828;
    (&_S5208)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5209;
    (&_S5209)->primal_0 = _S4833;
    (&_S5209)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5208, &_S5209, _S5207);
    float _S5210 = - _S5209.differential_0.y;
    float _S5211 = _S4832 * _S5209.differential_0.x;
    float _S5212 = - (_S4824 * _S5209.differential_0.x);
    float3  _S5213 = _S5203.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5214;
    (&_S5214)->primal_0 = _S4828;
    (&_S5214)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5215;
    (&_S5215)->primal_0 = _S4831;
    (&_S5215)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5214, &_S5215, _S5213);
    float _S5216 = _S4824 * _S5215.differential_0.x;
    float _S5217 = _S4830 * _S5215.differential_0.x;
    float3  _S5218 = _S5206.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5219;
    (&_S5219)->primal_0 = _S4828;
    (&_S5219)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5220;
    (&_S5220)->primal_0 = _S4829;
    (&_S5220)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5219, &_S5220, _S5218);
    Matrix<float, 3, 3>  _S5221 = transpose_0(_S5208.differential_0 + _S5214.differential_0 + _S5219.differential_0);
    float _S5222 = 2.0f * - _S5221.rows[int(2)].z;
    float _S5223 = 2.0f * _S5221.rows[int(2)].y;
    float _S5224 = 2.0f * _S5221.rows[int(2)].x;
    float _S5225 = 2.0f * _S5221.rows[int(1)].z;
    float _S5226 = 2.0f * - _S5221.rows[int(1)].y;
    float _S5227 = 2.0f * _S5221.rows[int(1)].x;
    float _S5228 = 2.0f * _S5221.rows[int(0)].z;
    float _S5229 = 2.0f * _S5221.rows[int(0)].y;
    float _S5230 = 2.0f * - _S5221.rows[int(0)].x;
    float _S5231 = - _S5227 + _S5229;
    float _S5232 = _S5224 + - _S5228;
    float _S5233 = - _S5223 + _S5225;
    float _S5234 = _S5223 + _S5225;
    float _S5235 = _S5224 + _S5228;
    float _S5236 = _S5227 + _S5229;
    float _S5237 = quat_32.w * (_S5226 + _S5230);
    float _S5238 = quat_32.z * (_S5222 + _S5230);
    float _S5239 = quat_32.y * (_S5222 + _S5226);
    float _S5240 = quat_32.x * _S5231 + quat_32.z * _S5234 + quat_32.y * _S5235 + _S5237 + _S5237;
    float _S5241 = quat_32.x * _S5232 + quat_32.w * _S5234 + quat_32.y * _S5236 + _S5238 + _S5238;
    float _S5242 = quat_32.x * _S5233 + quat_32.w * _S5235 + quat_32.z * _S5236 + _S5239 + _S5239;
    float _S5243 = quat_32.w * _S5231 + quat_32.z * _S5232 + quat_32.y * _S5233;
    float _S5244 = _S5212 + _S5216;
    float _S5245 = 0.5f * - _S5244;
    float _S5246 = _S5210 + _S5215.differential_0.y;
    DiffPair_float_0 _S5247;
    (&_S5247)->primal_0 = _S4825;
    (&_S5247)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5247, _S5246);
    float _S5248 = _S5245 + _S5247.differential_0;
    float _S5249 = _S5211 + _S5217 + _S5220.differential_0.x;
    DiffPair_float_0 _S5250;
    (&_S5250)->primal_0 = _S4823;
    (&_S5250)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5250, _S5249);
    float _S5251 = _S5245 + _S5250.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5252;
    (&_S5252)->primal_0 = mean_c_26;
    (&_S5252)->differential_0 = _S4964;
    s_bwd_length_impl_0(&_S5252, 0.0f);
    float3  _S5253 = _S5252.differential_0 + _S4967.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5254;
    (&_S5254)->primal_0 = R_30;
    (&_S5254)->differential_0 = _S5045;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5255;
    (&_S5255)->primal_0 = mean_31;
    (&_S5255)->differential_0 = _S4964;
    s_bwd_prop_mul_1(&_S5254, &_S5255, _S5253);
    float3  _S5256 = _S5198 + _S5201 + _S5204 + _S5253 + _S5048.differential_0;
    Matrix<float, 3, 3>  _S5257 = _S5199.differential_0 + _S5202.differential_0 + _S5205.differential_0 + _S5254.differential_0 + _S5049;
    float3  _S5258 = make_float3 (_S5251, _S5248, _S5244);
    float4  _S5259 = make_float4 (0.0f);
    *&((&_S5259)->w) = _S5240;
    *&((&_S5259)->z) = _S5241;
    *&((&_S5259)->y) = _S5242;
    *&((&_S5259)->x) = _S5243;
    float4  _S5260 = _S5259;
    float3  _S5261 = _S5207 + _S5213 + _S5218 + _S5255.differential_0 + _S5043;
    *v_mean_10 = _S5261;
    *v_quat_9 = _S5260;
    *v_scale_9 = _S5258;
    *v_hardness_5 = _S5125;
    (*v_sh_coeffs_8)[int(0)] = _S5103;
    (*v_sh_coeffs_8)[int(1)] = _S5104;
    (*v_sh_coeffs_8)[int(2)] = _S5105;
    (*v_sh_coeffs_8)[int(3)] = _S5106;
    (*v_sh_coeffs_8)[int(4)] = _S5107;
    (*v_sh_coeffs_8)[int(5)] = _S5108;
    (*v_sh_coeffs_8)[int(6)] = _S5109;
    (*v_sh_coeffs_8)[int(7)] = _S5110;
    (*v_sh_coeffs_8)[int(8)] = _S5111;
    (*v_sh_coeffs_8)[int(9)] = _S5112;
    (*v_sh_coeffs_8)[int(10)] = _S5113;
    (*v_sh_coeffs_8)[int(11)] = _S5114;
    (*v_sh_coeffs_8)[int(12)] = _S5115;
    (*v_sh_coeffs_8)[int(13)] = _S5116;
    (*v_sh_coeffs_8)[int(14)] = _S5117;
    (*v_sh_coeffs_8)[int(15)] = _S5118;
    (*v_ch_coeffs_3)[int(0)] = _S5100;
    (*v_ch_coeffs_3)[int(1)] = _S5101;
    *v_R_9 = _S5257;
    *v_t_9 = _S5256;
    return;
}

inline __device__ bool ray_triangle_intersection_uvt(float3  ray_o_5, float3  ray_d_5, FixedArray<float3 , 3>  * verts_4, float * u_80, float * v_80, float * t_30)
{
    float3  v1v0_0 = (*verts_4)[int(1)] - (*verts_4)[int(0)];
    float3  v2v0_0 = (*verts_4)[int(2)] - (*verts_4)[int(0)];
    float3  rov0_0 = ray_o_5 - (*verts_4)[int(0)];
    float3  n_0 = cross_0(v1v0_0, v2v0_0);
    float3  q_2 = cross_0(rov0_0, ray_d_5);
    float d_28 = 1.0f / dot_0(ray_d_5, n_0);
    *u_80 = d_28 * dot_0(- q_2, v2v0_0);
    *v_80 = d_28 * dot_0(q_2, v1v0_0);
    *t_30 = d_28 * dot_0(- n_0, rov0_0);
    bool _S5262;
    if((*u_80) >= 0.0f)
    {
        _S5262 = (*v_80) >= 0.0f;
    }
    else
    {
        _S5262 = false;
    }
    if(_S5262)
    {
        _S5262 = (*u_80 + *v_80) <= 1.0f;
    }
    else
    {
        _S5262 = false;
    }
    if(_S5262)
    {
        _S5262 = (*t_30) >= 0.0f;
    }
    else
    {
        _S5262 = false;
    }
    return _S5262;
}

inline __device__ float evaluate_alpha_opaque_triangle(FixedArray<float3 , 3>  * verts_5, float2  hardness_16, float3  ray_o_6, float3  ray_d_6)
{
    float3  v1v0_1 = (*verts_5)[int(1)] - (*verts_5)[int(0)];
    float3  v2v0_1 = (*verts_5)[int(2)] - (*verts_5)[int(0)];
    float3  rov0_1 = ray_o_6 - (*verts_5)[int(0)];
    float3  n_1 = cross_0(v1v0_1, v2v0_1);
    float3  q_3 = cross_0(rov0_1, ray_d_6);
    float d_29 = 1.0f / dot_0(ray_d_6, n_1);
    float u_81 = d_29 * dot_0(- q_3, v2v0_1);
    float v_81 = d_29 * dot_0(q_3, v1v0_1);
    float t_31 = d_29 * dot_0(- n_1, rov0_1);
    bool _S5263;
    if(u_81 >= 0.0f)
    {
        _S5263 = v_81 >= 0.0f;
    }
    else
    {
        _S5263 = false;
    }
    if(_S5263)
    {
        _S5263 = (u_81 + v_81) <= 1.0f;
    }
    else
    {
        _S5263 = false;
    }
    if(_S5263)
    {
        _S5263 = t_31 >= 0.0f;
    }
    else
    {
        _S5263 = false;
    }
    if(!_S5263)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_81), (v_81)))), ((F32_sqrt((0.5f))) * (1.0f - u_81 - v_81)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_16.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_16.x;
    float _S5264;
    if(opac_0 < 0.0f)
    {
        _S5264 = 0.0f;
    }
    else
    {
        _S5264 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S5264;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_4, float _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5265 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5266 = *dpray_d_4;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_4).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S5267 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S5268 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_4).primal_0);
    float _S5269 = s_primal_ctx_dot_0((*dpray_d_4).primal_0, _S5267);
    float d_30 = 1.0f / _S5269;
    float _S5270 = _S5269 * _S5269;
    float3  _S5271 = - _S5268;
    float _S5272 = s_primal_ctx_dot_0(_S5271, v2v0_2);
    float u_82 = d_30 * _S5272;
    float _S5273 = s_primal_ctx_dot_0(_S5268, v1v0_2);
    float v_82 = d_30 * _S5273;
    float3  _S5274 = - _S5267;
    float t_32 = d_30 * s_primal_ctx_dot_0(_S5274, rov0_2);
    bool _S5275;
    if(u_82 >= 0.0f)
    {
        _S5275 = v_82 >= 0.0f;
    }
    else
    {
        _S5275 = false;
    }
    if(_S5275)
    {
        _S5275 = (u_82 + v_82) <= 1.0f;
    }
    else
    {
        _S5275 = false;
    }
    if(_S5275)
    {
        _S5275 = t_32 >= 0.0f;
    }
    else
    {
        _S5275 = false;
    }
    bool _S5276 = !!_S5275;
    float _S5277;
    float _S5278;
    float _S5279;
    float _S5280;
    float _S5281;
    float _S5282;
    float _S5283;
    float _S5284;
    float _S5285;
    float _S5286;
    float _S5287;
    if(_S5276)
    {
        float _S5288 = s_primal_ctx_min_0(u_82, v_82);
        float _S5289 = s_primal_ctx_sqrt_0(0.5f);
        float _S5290 = _S5289 * (1.0f - u_82 - v_82);
        float _S5291 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S5288, _S5290) * _S5291;
        float _S5292 = _S5265.primal_0.y;
        float _S5293 = 1.0f - opac_1;
        float _S5294 = 1.0f - s_primal_ctx_clamp_0(_S5292, 0.0f, 0.99989998340606689f);
        float _S5295 = 1.0f / _S5294;
        float _S5296 = _S5294 * _S5294;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S5293, _S5295);
        float o_1 = _S5265.primal_0.x;
        bool _S5297 = opac_1 < 0.0f;
        if(_S5297)
        {
            _S5277 = 0.0f;
        }
        else
        {
            _S5277 = o_1 * w_1;
        }
        _S5275 = _S5297;
        _S5278 = o_1;
        _S5279 = w_1;
        _S5280 = _S5293;
        _S5281 = _S5295;
        _S5282 = _S5296;
        _S5283 = _S5292;
        _S5284 = _S5291;
        _S5285 = _S5288;
        _S5286 = _S5290;
        _S5287 = _S5289;
    }
    else
    {
        _S5275 = false;
        _S5277 = 0.0f;
        _S5278 = 0.0f;
        _S5279 = 0.0f;
        _S5280 = 0.0f;
        _S5281 = 0.0f;
        _S5282 = 0.0f;
        _S5283 = 0.0f;
        _S5284 = 0.0f;
        _S5285 = 0.0f;
        _S5286 = 0.0f;
        _S5287 = 0.0f;
    }
    float2  _S5298 = make_float2 (0.0f);
    float2  _S5299;
    if(_S5276)
    {
        if(_S5275)
        {
            _S5277 = 0.0f;
            _S5278 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S5300;
            (&_S5300)->primal_0 = _S5277;
            (&_S5300)->differential_0 = 0.0f;
            DiffPair_float_0 _S5301;
            (&_S5301)->primal_0 = 0.99500000476837158f;
            (&_S5301)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S5300, &_S5301, _s_dOut_11);
            float _S5302 = _S5278 * _S5300.differential_0;
            _S5277 = _S5279 * _S5300.differential_0;
            _S5278 = _S5302;
        }
        float _S5303 = - _S5278;
        DiffPair_float_0 _S5304;
        (&_S5304)->primal_0 = _S5280;
        (&_S5304)->differential_0 = 0.0f;
        DiffPair_float_0 _S5305;
        (&_S5305)->primal_0 = _S5281;
        (&_S5305)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S5304, &_S5305, _S5303);
        float _S5306 = - - (_S5305.differential_0 / _S5282);
        float s_diff_opac_T_0 = - _S5304.differential_0;
        DiffPair_float_0 _S5307;
        (&_S5307)->primal_0 = _S5283;
        (&_S5307)->differential_0 = 0.0f;
        DiffPair_float_0 _S5308;
        (&_S5308)->primal_0 = 0.0f;
        (&_S5308)->differential_0 = 0.0f;
        DiffPair_float_0 _S5309;
        (&_S5309)->primal_0 = 0.99989998340606689f;
        (&_S5309)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S5307, &_S5308, &_S5309, _S5306);
        float _S5310 = _S5284 * s_diff_opac_T_0;
        DiffPair_float_0 _S5311;
        (&_S5311)->primal_0 = _S5285;
        (&_S5311)->differential_0 = 0.0f;
        DiffPair_float_0 _S5312;
        (&_S5312)->primal_0 = _S5286;
        (&_S5312)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5311, &_S5312, _S5310);
        float _S5313 = - (_S5287 * _S5312.differential_0);
        DiffPair_float_0 _S5314;
        (&_S5314)->primal_0 = u_82;
        (&_S5314)->differential_0 = 0.0f;
        DiffPair_float_0 _S5315;
        (&_S5315)->primal_0 = v_82;
        (&_S5315)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5314, &_S5315, _S5311.differential_0);
        float2  _S5316 = make_float2 (_S5277, _S5307.differential_0);
        float _S5317 = _S5313 + _S5315.differential_0;
        _S5277 = _S5313 + _S5314.differential_0;
        _S5278 = _S5317;
        _S5299 = _S5316;
    }
    else
    {
        _S5277 = 0.0f;
        _S5278 = 0.0f;
        _S5299 = _S5298;
    }
    float3  _S5318 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5319;
    (&_S5319)->primal_0 = _S5274;
    (&_S5319)->differential_0 = _S5318;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5320;
    (&_S5320)->primal_0 = rov0_2;
    (&_S5320)->differential_0 = _S5318;
    s_bwd_prop_dot_0(&_S5319, &_S5320, 0.0f);
    float3  _S5321 = - _S5319.differential_0;
    float _S5322 = d_30 * _S5278;
    float _S5323 = _S5273 * _S5278;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5324;
    (&_S5324)->primal_0 = _S5268;
    (&_S5324)->differential_0 = _S5318;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5325;
    (&_S5325)->primal_0 = v1v0_2;
    (&_S5325)->differential_0 = _S5318;
    s_bwd_prop_dot_0(&_S5324, &_S5325, _S5322);
    float _S5326 = d_30 * _S5277;
    float _S5327 = _S5272 * _S5277;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5328;
    (&_S5328)->primal_0 = _S5271;
    (&_S5328)->differential_0 = _S5318;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5329;
    (&_S5329)->primal_0 = v2v0_2;
    (&_S5329)->differential_0 = _S5318;
    s_bwd_prop_dot_0(&_S5328, &_S5329, _S5326);
    float3  _S5330 = - _S5328.differential_0;
    float _S5331 = - ((_S5323 + _S5327) / _S5270);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5332;
    (&_S5332)->primal_0 = _S5266.primal_0;
    (&_S5332)->differential_0 = _S5318;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5333;
    (&_S5333)->primal_0 = _S5267;
    (&_S5333)->differential_0 = _S5318;
    s_bwd_prop_dot_0(&_S5332, &_S5333, _S5331);
    float3  _S5334 = _S5324.differential_0 + _S5330;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5335;
    (&_S5335)->primal_0 = rov0_2;
    (&_S5335)->differential_0 = _S5318;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5336;
    (&_S5336)->primal_0 = _S5266.primal_0;
    (&_S5336)->differential_0 = _S5318;
    s_bwd_prop_cross_0(&_S5335, &_S5336, _S5334);
    float3  _S5337 = _S5321 + _S5333.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5338;
    (&_S5338)->primal_0 = v1v0_2;
    (&_S5338)->differential_0 = _S5318;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5339;
    (&_S5339)->primal_0 = v2v0_2;
    (&_S5339)->differential_0 = _S5318;
    s_bwd_prop_cross_0(&_S5338, &_S5339, _S5337);
    float3  _S5340 = _S5320.differential_0 + _S5335.differential_0;
    float3  _S5341 = _S5329.differential_0 + _S5339.differential_0;
    float3  _S5342 = _S5325.differential_0 + _S5338.differential_0;
    float3  _S5343 = - _S5340 + - _S5341 + - _S5342;
    float3  _S5344 = _S5332.differential_0 + _S5336.differential_0;
    dpray_d_4->primal_0 = (*dpray_d_4).primal_0;
    dpray_d_4->differential_0 = _S5344;
    dpray_o_4->primal_0 = (*dpray_o_4).primal_0;
    dpray_o_4->differential_0 = _S5340;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S5299;
    FixedArray<float3 , 3>  _S5345;
    _S5345[int(0)] = _S5318;
    _S5345[int(1)] = _S5318;
    _S5345[int(2)] = _S5318;
    _S5345[int(2)] = _S5341;
    _S5345[int(0)] = _S5343;
    _S5345[int(1)] = _S5342;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S5345;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5346, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S5347, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5348, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5349, float _S5350)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S5346, _S5347, _S5348, _S5349, _S5350);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_6, float2  hardness_17, float3  ray_o_7, float3  ray_d_7, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S5351 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5352 = { _S5351, _S5351, _S5351 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_6;
    (&dp_verts_0)->differential_0 = _S5352;
    float2  _S5353 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S5353;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_7;
    (&dp_ray_o_2)->differential_0 = _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_7;
    (&dp_ray_d_2)->differential_0 = _S5351;
    s_bwd_evaluate_alpha_opaque_triangle_0(&dp_verts_0, &dp_hardness_2, &dp_ray_o_2, &dp_ray_d_2, v_alpha_3);
    *v_verts_2 = (&dp_verts_0)->differential_0;
    *v_hardness_6 = dp_hardness_2.differential_0;
    *v_ray_o_3 = dp_ray_o_2.differential_0;
    *v_ray_d_3 = dp_ray_d_2.differential_0;
    return;
}

inline __device__ float evaluate_sorting_depth_opaque_triangle(FixedArray<float3 , 3>  * verts_7, FixedArray<float3 , 3>  * rgbs_4, float3  ray_o_8, float3  ray_d_8)
{
    float3  n_2 = cross_0((*verts_7)[int(1)] - (*verts_7)[int(0)], (*verts_7)[int(2)] - (*verts_7)[int(0)]);
    return 1.0f / dot_0(ray_d_8, n_2) * dot_0(- n_2, ray_o_8 - (*verts_7)[int(0)]);
}

inline __device__ void evaluate_color_opaque_triangle(FixedArray<float3 , 3>  * verts_8, FixedArray<float3 , 3>  * rgbs_5, float3  ray_o_9, float3  ray_d_9, float3  * color_13, float * depth_20)
{
    float3  v1v0_3 = (*verts_8)[int(1)] - (*verts_8)[int(0)];
    float3  v2v0_3 = (*verts_8)[int(2)] - (*verts_8)[int(0)];
    float3  rov0_3 = ray_o_9 - (*verts_8)[int(0)];
    float3  n_3 = cross_0(v1v0_3, v2v0_3);
    float3  q_4 = cross_0(rov0_3, ray_d_9);
    float d_31 = 1.0f / dot_0(ray_d_9, n_3);
    float u_83 = d_31 * dot_0(- q_4, v2v0_3);
    float v_83 = d_31 * dot_0(q_4, v1v0_3);
    *depth_20 = d_31 * dot_0(- n_3, rov0_3);
    *color_13 = (*rgbs_5)[int(0)] * make_float3 (1.0f - u_83 - v_83) + (*rgbs_5)[int(1)] * make_float3 (u_83) + (*rgbs_5)[int(2)] * make_float3 (v_83);
    *depth_20 = (F32_log(((F32_max((*depth_20), (9.999999960041972e-13f))))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float3  dpcolor_1, float dpdepth_2)
{
    float3  v1v0_4 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_4 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_4 = (*dpray_o_5).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S5354 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S5355 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_5).primal_0);
    float _S5356 = s_primal_ctx_dot_0((*dpray_d_5).primal_0, _S5354);
    float d_32 = 1.0f / _S5356;
    float _S5357 = _S5356 * _S5356;
    float3  _S5358 = - _S5355;
    float _S5359 = s_primal_ctx_dot_0(_S5358, v2v0_4);
    float u_84 = d_32 * _S5359;
    float3  _S5360 = make_float3 (u_84);
    float _S5361 = s_primal_ctx_dot_0(_S5355, v1v0_4);
    float v_84 = d_32 * _S5361;
    float3  _S5362 = make_float3 (v_84);
    float3  _S5363 = - _S5354;
    float _S5364 = s_primal_ctx_dot_0(_S5363, rov0_4);
    float _S5365 = d_32 * _S5364;
    float3  _S5366 = make_float3 (1.0f - u_84 - v_84);
    DiffPair_float_0 _S5367;
    (&_S5367)->primal_0 = s_primal_ctx_max_0(_S5365, 9.999999960041972e-13f);
    (&_S5367)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5367, dpdepth_2);
    DiffPair_float_0 _S5368;
    (&_S5368)->primal_0 = _S5365;
    (&_S5368)->differential_0 = 0.0f;
    DiffPair_float_0 _S5369;
    (&_S5369)->primal_0 = 9.999999960041972e-13f;
    (&_S5369)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5368, &_S5369, _S5367.differential_0);
    float3  _S5370 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S5371 = _S5362 * dpcolor_1;
    float3  _S5372 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S5373 = _S5360 * dpcolor_1;
    float3  _S5374 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S5375 = _S5366 * dpcolor_1;
    float _S5376 = - (_S5374.x + _S5374.y + _S5374.z);
    float _S5377 = d_32 * _S5368.differential_0;
    float _S5378 = _S5364 * _S5368.differential_0;
    float3  _S5379 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5380;
    (&_S5380)->primal_0 = _S5363;
    (&_S5380)->differential_0 = _S5379;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5381;
    (&_S5381)->primal_0 = rov0_4;
    (&_S5381)->differential_0 = _S5379;
    s_bwd_prop_dot_0(&_S5380, &_S5381, _S5377);
    float3  _S5382 = - _S5380.differential_0;
    float _S5383 = _S5376 + _S5370.x + _S5370.y + _S5370.z;
    float _S5384 = d_32 * _S5383;
    float _S5385 = _S5361 * _S5383;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5386;
    (&_S5386)->primal_0 = _S5355;
    (&_S5386)->differential_0 = _S5379;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5387;
    (&_S5387)->primal_0 = v1v0_4;
    (&_S5387)->differential_0 = _S5379;
    s_bwd_prop_dot_0(&_S5386, &_S5387, _S5384);
    float _S5388 = _S5376 + _S5372.x + _S5372.y + _S5372.z;
    float _S5389 = d_32 * _S5388;
    float _S5390 = _S5359 * _S5388;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5391;
    (&_S5391)->primal_0 = _S5358;
    (&_S5391)->differential_0 = _S5379;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5392;
    (&_S5392)->primal_0 = v2v0_4;
    (&_S5392)->differential_0 = _S5379;
    s_bwd_prop_dot_0(&_S5391, &_S5392, _S5389);
    float3  _S5393 = - _S5391.differential_0;
    float _S5394 = - ((_S5378 + _S5385 + _S5390) / _S5357);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5395;
    (&_S5395)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5395)->differential_0 = _S5379;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5396;
    (&_S5396)->primal_0 = _S5354;
    (&_S5396)->differential_0 = _S5379;
    s_bwd_prop_dot_0(&_S5395, &_S5396, _S5394);
    float3  _S5397 = _S5386.differential_0 + _S5393;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5398;
    (&_S5398)->primal_0 = rov0_4;
    (&_S5398)->differential_0 = _S5379;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5399;
    (&_S5399)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5399)->differential_0 = _S5379;
    s_bwd_prop_cross_0(&_S5398, &_S5399, _S5397);
    float3  _S5400 = _S5382 + _S5396.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5401;
    (&_S5401)->primal_0 = v1v0_4;
    (&_S5401)->differential_0 = _S5379;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5402;
    (&_S5402)->primal_0 = v2v0_4;
    (&_S5402)->differential_0 = _S5379;
    s_bwd_prop_cross_0(&_S5401, &_S5402, _S5400);
    float3  _S5403 = _S5381.differential_0 + _S5398.differential_0;
    float3  _S5404 = _S5392.differential_0 + _S5402.differential_0;
    float3  _S5405 = _S5387.differential_0 + _S5401.differential_0;
    float3  _S5406 = - _S5403 + - _S5404 + - _S5405;
    float3  _S5407 = _S5395.differential_0 + _S5399.differential_0;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S5407;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S5403;
    FixedArray<float3 , 3>  _S5408;
    _S5408[int(0)] = _S5379;
    _S5408[int(1)] = _S5379;
    _S5408[int(2)] = _S5379;
    _S5408[int(2)] = _S5371;
    _S5408[int(1)] = _S5373;
    _S5408[int(0)] = _S5375;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S5408;
    FixedArray<float3 , 3>  _S5409;
    _S5409[int(0)] = _S5379;
    _S5409[int(1)] = _S5379;
    _S5409[int(2)] = _S5379;
    _S5409[int(2)] = _S5404;
    _S5409[int(0)] = _S5406;
    _S5409[int(1)] = _S5405;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S5409;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5410, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5411, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5412, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5413, float3  _S5414, float _S5415)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S5410, _S5411, _S5412, _S5413, _S5414, _S5415);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_9, FixedArray<float3 , 3>  * rgbs_6, float3  ray_o_10, float3  ray_d_10, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S5416 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5417 = { _S5416, _S5416, _S5416 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_9;
    (&dp_verts_1)->differential_0 = _S5417;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S5417;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_10;
    (&dp_ray_o_3)->differential_0 = _S5416;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_10;
    (&dp_ray_d_3)->differential_0 = _S5416;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

