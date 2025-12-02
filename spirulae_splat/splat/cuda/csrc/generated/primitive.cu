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

inline __device__ Matrix<float, 2, 2>  camera_distortion_jac_0(float2  uv_0, FixedArray<float, 10>  * dist_coeffs_0)
{
    float u_0 = uv_0.x;
    float v_0 = uv_0.y;
    float r2_0 = u_0 * u_0 + v_0 * v_0;
    float _S4 = (*dist_coeffs_0)[int(2)] + r2_0 * (*dist_coeffs_0)[int(3)];
    float _S5 = (*dist_coeffs_0)[int(1)] + r2_0 * _S4;
    float _S6 = (*dist_coeffs_0)[int(0)] + r2_0 * _S5;
    float2  _S7 = make_float2 (1.0f + r2_0 * _S6);
    float _S8 = 2.0f * (*dist_coeffs_0)[int(4)];
    float _S9 = _S8 * u_0;
    float _S10 = 2.0f * u_0;
    float _S11 = 2.0f * (*dist_coeffs_0)[int(5)];
    float _S12 = _S11 * u_0;
    float _S13 = 2.0f * v_0;
    float2  _S14 = make_float2 (1.0f, 0.0f) + make_float2 ((*dist_coeffs_0)[int(8)], (*dist_coeffs_0)[int(9)]);
    float2  _S15 = uv_0 * _S14;
    float _S16 = (*dist_coeffs_0)[int(4)] * _S14.y;
    float _S17 = (*dist_coeffs_0)[int(5)] * _S14.x;
    float _S18 = _S15.x + _S15.y;
    float _S19 = r2_0 * _S18;
    float _S20 = r2_0 * _S19;
    float _S21 = (*dist_coeffs_0)[int(7)] * _S14.y + _S16 + (*dist_coeffs_0)[int(6)] * _S14.x + _S17 + _S6 * _S18 + _S5 * _S19 + _S4 * _S20 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S20);
    float _S22 = v_0 * _S21;
    float _S23 = u_0 * _S21;
    float2  _S24 = make_float2 (0.0f, 1.0f) + make_float2 (0.0f, 0.0f);
    float2  _S25 = uv_0 * _S24;
    float _S26 = (*dist_coeffs_0)[int(4)] * _S24.y;
    float _S27 = (*dist_coeffs_0)[int(5)] * _S24.x;
    float _S28 = _S25.x + _S25.y;
    float _S29 = r2_0 * _S28;
    float _S30 = r2_0 * _S29;
    float _S31 = (*dist_coeffs_0)[int(7)] * _S24.y + _S26 + (*dist_coeffs_0)[int(6)] * _S24.x + _S27 + _S6 * _S28 + _S5 * _S29 + _S4 * _S30 + (*dist_coeffs_0)[int(3)] * (r2_0 * _S30);
    float _S32 = v_0 * _S31;
    float _S33 = u_0 * _S31;
    return makeMatrix<float, 2, 2> (_S7 * _S14 + make_float2 (_S11 * (v_0 * _S14.y) + _S10 * _S17 + 2.0f * (u_0 * _S17) + _S8 * (v_0 * _S14.x) + _S23 + _S23, _S13 * _S16 + 2.0f * (v_0 * _S16) + _S12 * _S14.y + _S9 * _S14.x + _S22 + _S22), _S7 * _S24 + make_float2 (_S11 * (v_0 * _S24.y) + _S10 * _S27 + 2.0f * (u_0 * _S27) + _S8 * (v_0 * _S24.x) + _S33 + _S33, _S13 * _S26 + 2.0f * (v_0 * _S26) + _S12 * _S24.y + _S9 * _S24.x + _S32 + _S32));
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
    DiffPair_float_0 _S34 = *dpx_0;
    float _S35;
    if(((*dpx_0).primal_0) < ((*dpy_0).primal_0))
    {
        _S35 = dOut_4;
    }
    else
    {
        if(((*dpx_0).primal_0) > ((*dpy_0).primal_0))
        {
            _S35 = 0.0f;
        }
        else
        {
            _S35 = 0.5f * dOut_4;
        }
    }
    dpx_0->primal_0 = _S34.primal_0;
    dpx_0->differential_0 = _S35;
    DiffPair_float_0 _S36 = *dpy_0;
    if(((*dpy_0).primal_0) < (_S34.primal_0))
    {
        _S35 = dOut_4;
    }
    else
    {
        if(((*dpy_0).primal_0) > ((*dpx_0).primal_0))
        {
            _S35 = 0.0f;
        }
        else
        {
            _S35 = 0.5f * dOut_4;
        }
    }
    dpy_0->primal_0 = _S36.primal_0;
    dpy_0->differential_0 = _S35;
    return;
}

inline __device__ bool is_valid_distortion(float2  uv_1, FixedArray<float, 10>  * dist_coeffs_1)
{
    Matrix<float, 2, 2>  _S37 = camera_distortion_jac_0(uv_1, dist_coeffs_1);
    return (F32_min((determinant_0(_S37)), ((F32_min((_S37.rows[int(0)].x), (_S37.rows[int(1)].y)))))) > 0.0f;
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

inline __device__ float2  distort_point(float2  uv_2, bool is_fisheye_0, FixedArray<float, 10>  * dist_coeffs_2)
{
    float2  _S42;
    if(is_fisheye_0)
    {
        float r_6 = length_0(uv_2);
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
        _S42 = uv_2 * make_float2 (_S43);
    }
    else
    {
        _S42 = uv_2;
    }
    float u_1 = _S42.x;
    float v_1 = _S42.y;
    float r2_1 = u_1 * u_1 + v_1 * v_1;
    float2  _S44 = _S42 * make_float2 (1.0f + r2_1 * ((*dist_coeffs_2)[int(0)] + r2_1 * ((*dist_coeffs_2)[int(1)] + r2_1 * ((*dist_coeffs_2)[int(2)] + r2_1 * (*dist_coeffs_2)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_2)[int(4)] * u_1 * v_1 + (*dist_coeffs_2)[int(5)] * (r2_1 + 2.0f * u_1 * u_1) + (*dist_coeffs_2)[int(6)] * r2_1, 2.0f * (*dist_coeffs_2)[int(5)] * u_1 * v_1 + (*dist_coeffs_2)[int(4)] * (r2_1 + 2.0f * v_1 * v_1) + (*dist_coeffs_2)[int(7)] * r2_1);
    return _S44 + make_float2 ((*dist_coeffs_2)[int(8)] * _S44.x + (*dist_coeffs_2)[int(9)] * _S44.y, 0.0f);
}

inline __device__ bool undistort_point_0(float2  uv_3, FixedArray<float, 10>  * dist_coeffs_3, int maxiter_0, float2  * uv_undist_0)
{
    int i_6 = int(0);
    float2  q_0 = uv_3;
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
        float2  _S45 = q_0 * make_float2 (1.0f + r2_2 * ((*dist_coeffs_3)[int(0)] + r2_2 * ((*dist_coeffs_3)[int(1)] + r2_2 * ((*dist_coeffs_3)[int(2)] + r2_2 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_2 * v_2 + (*dist_coeffs_3)[int(5)] * (r2_2 + 2.0f * u_2 * u_2) + (*dist_coeffs_3)[int(6)] * r2_2, 2.0f * (*dist_coeffs_3)[int(5)] * u_2 * v_2 + (*dist_coeffs_3)[int(4)] * (r2_2 + 2.0f * v_2 * v_2) + (*dist_coeffs_3)[int(7)] * r2_2);
        float2  r_7 = _S45 + make_float2 ((*dist_coeffs_3)[int(8)] * _S45.x + (*dist_coeffs_3)[int(9)] * _S45.y, 0.0f) - uv_3;
        Matrix<float, 2, 2>  _S46 = camera_distortion_jac_0(q_0, dist_coeffs_3);
        float inv_det_0 = 1.0f / (_S46.rows[int(0)].x * _S46.rows[int(1)].y - _S46.rows[int(0)].y * _S46.rows[int(1)].x);
        float _S47 = r_7.x;
        float _S48 = r_7.y;
        float2  q_1 = q_0 - make_float2 ((_S47 * _S46.rows[int(1)].y - _S48 * _S46.rows[int(0)].y) * inv_det_0, (- _S47 * _S46.rows[int(1)].x + _S48 * _S46.rows[int(0)].x) * inv_det_0);
        i_6 = i_6 + int(1);
        q_0 = q_1;
    }
    *uv_undist_0 = q_0;
    Matrix<float, 2, 2>  _S49 = camera_distortion_jac_0(q_0, dist_coeffs_3);
    bool _S50;
    if((F32_min((determinant_0(_S49)), ((F32_min((_S49.rows[int(0)].x), (_S49.rows[int(1)].y)))))) > 0.0f)
    {
        float u_3 = (*uv_undist_0).x;
        float v_3 = (*uv_undist_0).y;
        float r2_3 = u_3 * u_3 + v_3 * v_3;
        float2  _S51 = *uv_undist_0 * make_float2 (1.0f + r2_3 * ((*dist_coeffs_3)[int(0)] + r2_3 * ((*dist_coeffs_3)[int(1)] + r2_3 * ((*dist_coeffs_3)[int(2)] + r2_3 * (*dist_coeffs_3)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_3)[int(4)] * u_3 * v_3 + (*dist_coeffs_3)[int(5)] * (r2_3 + 2.0f * u_3 * u_3) + (*dist_coeffs_3)[int(6)] * r2_3, 2.0f * (*dist_coeffs_3)[int(5)] * u_3 * v_3 + (*dist_coeffs_3)[int(4)] * (r2_3 + 2.0f * v_3 * v_3) + (*dist_coeffs_3)[int(7)] * r2_3);
        _S50 = (length_0(_S51 + make_float2 ((*dist_coeffs_3)[int(8)] * _S51.x + (*dist_coeffs_3)[int(9)] * _S51.y, 0.0f) - uv_3)) < 0.00999999977648258f;
    }
    else
    {
        _S50 = false;
    }
    return _S50;
}

inline __device__ bool undistort_point(float2  uv_4, bool is_fisheye_1, FixedArray<float, 10>  * dist_coeffs_4, float2  * uv_undist_1)
{
    float2  _S52 = uv_4;
    bool _S53 = undistort_point_0(uv_4, dist_coeffs_4, int(8), &_S52);
    if(!_S53)
    {
        return false;
    }
    float3  raydir_0;
    if(is_fisheye_1)
    {
        float2  _S54 = _S52;
        float theta_1 = length_0(_S52);
        float _S55;
        if(theta_1 < 0.00100000004749745f)
        {
            _S55 = 1.0f - theta_1 * theta_1 / 6.0f;
        }
        else
        {
            _S55 = (F32_sin((theta_1))) / theta_1;
        }
        float3  _S56 = make_float3 ((_S54 * make_float2 (_S55)).x, (_S54 * make_float2 (_S55)).y, (F32_cos((theta_1))));
        raydir_0 = _S56;
    }
    else
    {
        raydir_0 = make_float3 (_S52.x, _S52.y, 1.0f);
    }
    *uv_undist_1 = float2 {raydir_0.x, raydir_0.y} / make_float2 ((F32_max((raydir_0.z), (9.999999960041972e-13f))));
    return true;
}

inline __device__ bool unproject_point(float2  uv_5, bool is_fisheye_2, FixedArray<float, 10>  * dist_coeffs_5, float3  * raydir_1)
{
    float2  _S57 = uv_5;
    int3  _S58 = make_int3 (int(0));
    float3  _S59 = make_float3 ((float)_S58.x, (float)_S58.y, (float)_S58.z);
    *raydir_1 = _S59;
    bool _S60 = undistort_point_0(uv_5, dist_coeffs_5, int(8), &_S57);
    if(!_S60)
    {
        return false;
    }
    float3  _S61;
    if(is_fisheye_2)
    {
        float2  _S62 = _S57;
        float theta_2 = length_0(_S57);
        float _S63;
        if(theta_2 < 0.00100000004749745f)
        {
            _S63 = 1.0f - theta_2 * theta_2 / 6.0f;
        }
        else
        {
            _S63 = (F32_sin((theta_2))) / theta_2;
        }
        float3  _S64 = make_float3 ((_S62 * make_float2 (_S63)).x, (_S62 * make_float2 (_S63)).y, (F32_cos((theta_2))));
        _S61 = _S64;
    }
    else
    {
        _S61 = make_float3 (_S57.x, _S57.y, 1.0f);
    }
    *raydir_1 = _S61;
    return true;
}

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_8, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_8, float3  dOut_9)
{
    float _S65 = (*right_8).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 3, 3>  right_d_result_4;
    *&(((&right_d_result_4)->rows + (int(0)))->x) = (*left_8).primal_0.x * dOut_9.x;
    float sum_14 = _S65 + (*right_8).primal_0.rows[int(0)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(0)))->y) = (*left_8).primal_0.x * dOut_9.y;
    float sum_15 = sum_14 + (*right_8).primal_0.rows[int(0)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(0)))->z) = (*left_8).primal_0.x * dOut_9.z;
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = sum_15;
    float _S66 = (*right_8).primal_0.rows[int(1)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(1)))->x) = (*left_8).primal_0.y * dOut_9.x;
    float sum_16 = _S66 + (*right_8).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(1)))->y) = (*left_8).primal_0.y * dOut_9.y;
    float sum_17 = sum_16 + (*right_8).primal_0.rows[int(1)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(1)))->z) = (*left_8).primal_0.y * dOut_9.z;
    *&((&left_d_result_4)->y) = sum_17;
    float _S67 = (*right_8).primal_0.rows[int(2)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(2)))->x) = (*left_8).primal_0.z * dOut_9.x;
    float sum_18 = _S67 + (*right_8).primal_0.rows[int(2)].y * dOut_9.y;
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

inline __device__ bool generate_ray(Matrix<float, 3, 3>  R_2, float3  t_1, float2  uv_6, bool is_fisheye_3, FixedArray<float, 10>  * dist_coeffs_6, float3  * ray_o_0, float3  * ray_d_0)
{
    float2  _S68 = uv_6;
    *ray_o_0 = - mul_7(t_1, R_2);
    bool _S69 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S68);
    if(!_S69)
    {
        return false;
    }
    float3  raydir_2;
    if(is_fisheye_3)
    {
        float2  _S70 = _S68;
        float theta_3 = length_0(_S68);
        float _S71;
        if(theta_3 < 0.00100000004749745f)
        {
            _S71 = 1.0f - theta_3 * theta_3 / 6.0f;
        }
        else
        {
            _S71 = (F32_sin((theta_3))) / theta_3;
        }
        float3  _S72 = make_float3 ((_S70 * make_float2 (_S71)).x, (_S70 * make_float2 (_S71)).y, (F32_cos((theta_3))));
        raydir_2 = _S72;
    }
    else
    {
        raydir_2 = make_float3 (_S68.x, _S68.y, 1.0f);
    }
    *ray_d_0 = normalize_0(mul_7(raydir_2, R_2));
    return true;
}

struct s_bwd_prop_generate_ray_Intermediates_0
{
    float2  _S73;
    bool _S74;
};

inline __device__ float3  s_primal_ctx_mul_0(float3  _S75, Matrix<float, 3, 3>  _S76)
{
    return mul_7(_S75, _S76);
}

inline __device__ float s_primal_ctx_sin_0(float _S77)
{
    return (F32_sin((_S77)));
}

inline __device__ float s_primal_ctx_cos_0(float _S78)
{
    return (F32_cos((_S78)));
}

inline __device__ bool s_primal_ctx_generate_ray_0(Matrix<float, 3, 3>  dpR_0, float3  dpt_0, float2  uv_7, bool is_fisheye_4, FixedArray<float, 10>  * dist_coeffs_7, float3  * dpray_o_0, float3  * dpray_d_0, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_0)
{
    _s_diff_ctx_0->_S73 = make_float2 (0.0f);
    _s_diff_ctx_0->_S74 = false;
    float3  _S79 = make_float3 (0.0f);
    float3  _S80 = - s_primal_ctx_mul_0(dpt_0, dpR_0);
    float2  _S81 = uv_7;
    bool _S82 = undistort_point_0(uv_7, dist_coeffs_7, int(8), &_S81);
    _s_diff_ctx_0->_S73 = _S81;
    _s_diff_ctx_0->_S74 = _S82;
    float2  _S83 = _S81;
    float3  raydir_3;
    bool _S84;
    if(!!_S82)
    {
        if(is_fisheye_4)
        {
            float _S85 = length_0(_S83);
            float _S86;
            if(_S85 < 0.00100000004749745f)
            {
                _S86 = 1.0f - _S85 * _S85 / 6.0f;
            }
            else
            {
                _S86 = s_primal_ctx_sin_0(_S85) / _S85;
            }
            float3  _S87 = make_float3 ((_S83 * make_float2 (_S86)).x, (_S83 * make_float2 (_S86)).y, s_primal_ctx_cos_0(_S85));
            raydir_3 = _S87;
        }
        else
        {
            raydir_3 = make_float3 (_S83.x, _S83.y, 1.0f);
        }
        float3  _S88 = normalize_0(s_primal_ctx_mul_0(raydir_3, dpR_0));
        _S84 = true;
        raydir_3 = _S88;
    }
    else
    {
        _S84 = false;
        raydir_3 = _S79;
    }
    *dpray_o_0 = _S80;
    *dpray_d_0 = raydir_3;
    return _S84;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S89, float _S90)
{
    _d_sqrt_0(_S89, _S90);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_5, float _s_dOut_0)
{
    float _S91 = (*dpx_5).primal_0.x;
    float _S92 = (*dpx_5).primal_0.y;
    float _S93 = (*dpx_5).primal_0.z;
    DiffPair_float_0 _S94;
    (&_S94)->primal_0 = _S91 * _S91 + _S92 * _S92 + _S93 * _S93;
    (&_S94)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S94, _s_dOut_0);
    float _S95 = (*dpx_5).primal_0.z * _S94.differential_0;
    float _S96 = _S95 + _S95;
    float _S97 = (*dpx_5).primal_0.y * _S94.differential_0;
    float _S98 = _S97 + _S97;
    float _S99 = (*dpx_5).primal_0.x * _S94.differential_0;
    float _S100 = _S99 + _S99;
    float3  _S101 = make_float3 (0.0f);
    *&((&_S101)->z) = _S96;
    *&((&_S101)->y) = _S98;
    *&((&_S101)->x) = _S100;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S101;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S102, float _S103)
{
    s_bwd_prop_length_impl_0(_S102, _S103);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_6, float3  _s_dOut_1)
{
    float _S104 = length_1((*dpx_6).primal_0);
    float3  _S105 = (*dpx_6).primal_0 * _s_dOut_1;
    float3  _S106 = make_float3 (1.0f / _S104) * _s_dOut_1;
    float _S107 = - ((_S105.x + _S105.y + _S105.z) / (_S104 * _S104));
    float3  _S108 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S109;
    (&_S109)->primal_0 = (*dpx_6).primal_0;
    (&_S109)->differential_0 = _S108;
    s_bwd_length_impl_0(&_S109, _S107);
    float3  _S110 = _S106 + _S109.differential_0;
    dpx_6->primal_0 = (*dpx_6).primal_0;
    dpx_6->differential_0 = _S110;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S111, float3  _S112)
{
    s_bwd_prop_normalize_impl_0(_S111, _S112);
    return;
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S113, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S114, float3  _S115)
{
    _d_mul_1(_S113, _S114, _S115);
    return;
}

inline __device__ void s_bwd_prop_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_1, float2  uv_8, bool is_fisheye_5, FixedArray<float, 10>  * dist_coeffs_8, float3  dpray_o_1, float3  dpray_d_1, s_bwd_prop_generate_ray_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S116 = *dpR_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S117 = *dpt_1;
    float3  _S118 = make_float3 (0.0f);
    bool _S119 = !!_s_diff_ctx_1->_S74;
    float3  raydir_4;
    float3  _S120;
    if(_S119)
    {
        if(is_fisheye_5)
        {
            float _S121 = length_0(_s_diff_ctx_1->_S73);
            float _S122;
            if(_S121 < 0.00100000004749745f)
            {
                _S122 = 1.0f - _S121 * _S121 / 6.0f;
            }
            else
            {
                _S122 = s_primal_ctx_sin_0(_S121) / _S121;
            }
            float3  _S123 = make_float3 ((_s_diff_ctx_1->_S73 * make_float2 (_S122)).x, (_s_diff_ctx_1->_S73 * make_float2 (_S122)).y, s_primal_ctx_cos_0(_S121));
            raydir_4 = _S123;
        }
        else
        {
            raydir_4 = make_float3 (_s_diff_ctx_1->_S73.x, _s_diff_ctx_1->_S73.y, 1.0f);
        }
        float3  _S124 = raydir_4;
        raydir_4 = s_primal_ctx_mul_0(raydir_4, _S116.primal_0);
        _S120 = _S124;
    }
    else
    {
        raydir_4 = _S118;
        _S120 = _S118;
    }
    Matrix<float, 3, 3>  _S125 = makeMatrix<float, 3, 3> (0.0f);
    Matrix<float, 3, 3>  _S126;
    if(_S119)
    {
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S127;
        (&_S127)->primal_0 = raydir_4;
        (&_S127)->differential_0 = _S118;
        s_bwd_normalize_impl_0(&_S127, dpray_d_1);
        DiffPair_vectorx3Cfloatx2C3x3E_0 _S128;
        (&_S128)->primal_0 = _S120;
        (&_S128)->differential_0 = _S118;
        DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S129;
        (&_S129)->primal_0 = _S116.primal_0;
        (&_S129)->differential_0 = _S125;
        s_bwd_prop_mul_0(&_S128, &_S129, _S127.differential_0);
        _S126 = _S129.differential_0;
    }
    else
    {
        _S126 = _S125;
    }
    float3  _S130 = - dpray_o_1;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S131;
    (&_S131)->primal_0 = _S117.primal_0;
    (&_S131)->differential_0 = _S118;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S132;
    (&_S132)->primal_0 = _S116.primal_0;
    (&_S132)->differential_0 = _S125;
    s_bwd_prop_mul_0(&_S131, &_S132, _S130);
    dpt_1->primal_0 = (*dpt_1).primal_0;
    dpt_1->differential_0 = _S131.differential_0;
    Matrix<float, 3, 3>  _S133 = _S132.differential_0 + _S126;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S133;
    return;
}

inline __device__ void s_bwd_generate_ray_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S134, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S135, float2  _S136, bool _S137, FixedArray<float, 10>  * _S138, float3  _S139, float3  _S140)
{
    float3  _S141;
    float3  _S142;
    s_bwd_prop_generate_ray_Intermediates_0 _S143;
    bool _S144 = s_primal_ctx_generate_ray_0((*_S134).primal_0, (*_S135).primal_0, _S136, _S137, _S138, &_S141, &_S142, &_S143);
    s_bwd_prop_generate_ray_Intermediates_0 _S145 = _S143;
    s_bwd_prop_generate_ray_0(_S134, _S135, _S136, _S137, _S138, _S139, _S140, &_S145);
    return;
}

inline __device__ void generate_ray_vjp(Matrix<float, 3, 3>  R_3, float3  t_2, float2  uv_9, bool is_fisheye_6, FixedArray<float, 10>  * dist_coeffs_9, float3  v_ray_o_0, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S146 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_3;
    (&dp_R_0)->differential_0 = _S146;
    float3  _S147 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S147;
    s_bwd_generate_ray_0(&dp_R_0, &dp_t_0, uv_9, is_fisheye_6, dist_coeffs_9, v_ray_o_0, v_ray_d_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_7, float dOut_10)
{
    float _S148 = (F32_exp(((*dpx_7).primal_0))) * dOut_10;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S148;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_3, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S149 = scale_2.x;
    float sx_0 = (F32_exp((_S149)));
    float _S150 = scale_2.y;
    float sy_0 = (F32_exp((_S150)));
    float sz_0 = scale_2.z - 0.5f * (_S149 + _S150);
    float x_12 = quat_3.y;
    float x2_3 = x_12 * x_12;
    float y2_3 = quat_3.z * quat_3.z;
    float z2_3 = quat_3.w * quat_3.w;
    float xy_3 = quat_3.y * quat_3.z;
    float xz_3 = quat_3.y * quat_3.w;
    float yz_3 = quat_3.z * quat_3.w;
    float wx_3 = quat_3.x * quat_3.y;
    float wy_3 = quat_3.x * quat_3.z;
    float wz_3 = quat_3.x * quat_3.w;
    Matrix<float, 3, 3>  _S151 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S151, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S151, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S151, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S152 = float(width_0);
    float _S153 = float(height_0);
    float _S154 = 0.30000001192092896f * (0.5f * _S152 / fx_0);
    float _S155 = 0.30000001192092896f * (0.5f * _S153 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_0 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S152 - cx_0) / fx_0 + _S154), ((F32_max((- (cx_0 / fx_0 + _S154)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S153 - cy_0) / fy_0 + _S155), ((F32_max((- (cy_0 / fy_0 + _S155)), (mean3d_0.y * rz_0))))))) * rz2_0);
    *cov2d_0 = mul_6(mul_5(J_0, cov3d_0), transpose_1(J_0));
    *mean2d_0 = make_float2 (fx_0 * mean3d_0.x * rz_0 + cx_0, fy_0 * mean3d_0.y * rz_0 + cy_0);
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_8, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S156 = *dpx_8;
    bool _S157;
    if(((*dpx_8).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S157 = ((*dpx_8).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S157 = false;
    }
    float _S158;
    if(_S157)
    {
        _S158 = dOut_11;
    }
    else
    {
        _S158 = 0.0f;
    }
    dpx_8->primal_0 = _S156.primal_0;
    dpx_8->differential_0 = _S158;
    DiffPair_float_0 _S159 = *dpMin_0;
    if((_S156.primal_0) < ((*dpMin_0).primal_0))
    {
        _S158 = dOut_11;
    }
    else
    {
        _S158 = 0.0f;
    }
    dpMin_0->primal_0 = _S159.primal_0;
    dpMin_0->differential_0 = _S158;
    DiffPair_float_0 _S160 = *dpMax_0;
    if(((*dpx_8).primal_0) > ((*dpMax_0).primal_0))
    {
        _S158 = dOut_11;
    }
    else
    {
        _S158 = 0.0f;
    }
    dpMax_0->primal_0 = _S160.primal_0;
    dpMax_0->differential_0 = _S158;
    return;
}

inline __device__ float clamp_0(float x_13, float minBound_0, float maxBound_0)
{
    return (F32_min(((F32_max((x_13), (minBound_0)))), (maxBound_0)));
}

struct SigmaPoints_0
{
    FixedArray<float3 , 7>  p_0;
    FixedArray<float, 7>  w_mean_0;
    FixedArray<float, 7>  w_cov_0;
};

inline __device__ bool persp_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_10, uint width_1, uint height_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float2  * _S161;
    float2  * _S162;
    float2  * _S163;
    bool _S164;
    float2  * _S165;
    float2  * _S166;
    float2  * _S167;
    bool _S168;
    float2  * _S169;
    bool _S170;
    int2  _S171 = make_int2 (int(0));
    float2  _S172 = make_float2 ((float)_S171.x, (float)_S171.y);
    *mean2d_1 = _S172;
    *cov2d_1 = makeMatrix<float, 2, 2> (0.0f);
    float fx_1 = intrins_0.x;
    float fy_1 = intrins_0.y;
    float _S173 = float(width_1);
    float _S174 = float(height_1);
    float _S175 = 0.30000001192092896f * (0.5f * _S173 / fx_1) * fx_1;
    float lim_x_pos_0 = _S173 + _S175;
    float _S176 = 0.30000001192092896f * (0.5f * _S174 / fy_1) * fy_1;
    float lim_y_pos_0 = _S174 + _S176;
    FixedArray<float2 , 7>  proj_points_0;
    for(;;)
    {
        bool _S177;
        _S161 = &proj_points_0[int(0)];
        for(;;)
        {
            float _S178 = sigmas_0->p_0[int(0)].z;
            proj_points_0[int(0)] = float2 {sigmas_0->p_0[int(0)].x, sigmas_0->p_0[int(0)].y} / make_float2 (_S178);
            if(_S178 < 0.0f)
            {
                _S177 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S179 = camera_distortion_jac_0(proj_points_0[int(0)], dist_coeffs_10);
                _S177 = !((F32_min((determinant_0(_S179)), ((F32_min((_S179.rows[int(0)].x), (_S179.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S177)
            {
                break;
            }
            float u_4 = proj_points_0[int(0)].x;
            float v_4 = proj_points_0[int(0)].y;
            float r2_4 = u_4 * u_4 + v_4 * v_4;
            float2  _S180 = proj_points_0[int(0)] * make_float2 (1.0f + r2_4 * ((*dist_coeffs_10)[int(0)] + r2_4 * ((*dist_coeffs_10)[int(1)] + r2_4 * ((*dist_coeffs_10)[int(2)] + r2_4 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_4 * v_4 + (*dist_coeffs_10)[int(5)] * (r2_4 + 2.0f * u_4 * u_4) + (*dist_coeffs_10)[int(6)] * r2_4, 2.0f * (*dist_coeffs_10)[int(5)] * u_4 * v_4 + (*dist_coeffs_10)[int(4)] * (r2_4 + 2.0f * v_4 * v_4) + (*dist_coeffs_10)[int(7)] * r2_4);
            float2  _S181 = _S180 + make_float2 ((*dist_coeffs_10)[int(8)] * _S180.x + (*dist_coeffs_10)[int(9)] * _S180.y, 0.0f);
            proj_points_0[int(0)] = make_float2 (fx_1 * _S181.x + intrins_0.z, fy_1 * _S181.y + intrins_0.w);
            break;
        }
        bool all_valid_0 = true & (!_S177);
        _S162 = &proj_points_0[int(1)];
        for(;;)
        {
            float _S182 = sigmas_0->p_0[int(1)].z;
            proj_points_0[int(1)] = float2 {sigmas_0->p_0[int(1)].x, sigmas_0->p_0[int(1)].y} / make_float2 (_S182);
            if(_S182 < 0.0f)
            {
                _S177 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S183 = camera_distortion_jac_0(proj_points_0[int(1)], dist_coeffs_10);
                _S177 = !((F32_min((determinant_0(_S183)), ((F32_min((_S183.rows[int(0)].x), (_S183.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S177)
            {
                break;
            }
            float u_5 = proj_points_0[int(1)].x;
            float v_5 = proj_points_0[int(1)].y;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float2  _S184 = proj_points_0[int(1)] * make_float2 (1.0f + r2_5 * ((*dist_coeffs_10)[int(0)] + r2_5 * ((*dist_coeffs_10)[int(1)] + r2_5 * ((*dist_coeffs_10)[int(2)] + r2_5 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_5 * v_5 + (*dist_coeffs_10)[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + (*dist_coeffs_10)[int(6)] * r2_5, 2.0f * (*dist_coeffs_10)[int(5)] * u_5 * v_5 + (*dist_coeffs_10)[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + (*dist_coeffs_10)[int(7)] * r2_5);
            float2  _S185 = _S184 + make_float2 ((*dist_coeffs_10)[int(8)] * _S184.x + (*dist_coeffs_10)[int(9)] * _S184.y, 0.0f);
            proj_points_0[int(1)] = make_float2 (fx_1 * _S185.x + intrins_0.z, fy_1 * _S185.y + intrins_0.w);
            break;
        }
        bool all_valid_1 = all_valid_0 & (!_S177);
        for(;;)
        {
            _S163 = &proj_points_0[int(2)];
            for(;;)
            {
                float _S186 = sigmas_0->p_0[int(2)].z;
                proj_points_0[int(2)] = float2 {sigmas_0->p_0[int(2)].x, sigmas_0->p_0[int(2)].y} / make_float2 (_S186);
                if(_S186 < 0.0f)
                {
                    _S177 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S187 = camera_distortion_jac_0(proj_points_0[int(2)], dist_coeffs_10);
                    _S177 = !((F32_min((determinant_0(_S187)), ((F32_min((_S187.rows[int(0)].x), (_S187.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S177)
                {
                    break;
                }
                float u_6 = proj_points_0[int(2)].x;
                float v_6 = proj_points_0[int(2)].y;
                float r2_6 = u_6 * u_6 + v_6 * v_6;
                float2  _S188 = proj_points_0[int(2)] * make_float2 (1.0f + r2_6 * ((*dist_coeffs_10)[int(0)] + r2_6 * ((*dist_coeffs_10)[int(1)] + r2_6 * ((*dist_coeffs_10)[int(2)] + r2_6 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_6 * v_6 + (*dist_coeffs_10)[int(5)] * (r2_6 + 2.0f * u_6 * u_6) + (*dist_coeffs_10)[int(6)] * r2_6, 2.0f * (*dist_coeffs_10)[int(5)] * u_6 * v_6 + (*dist_coeffs_10)[int(4)] * (r2_6 + 2.0f * v_6 * v_6) + (*dist_coeffs_10)[int(7)] * r2_6);
                float2  _S189 = _S188 + make_float2 ((*dist_coeffs_10)[int(8)] * _S188.x + (*dist_coeffs_10)[int(9)] * _S188.y, 0.0f);
                proj_points_0[int(2)] = make_float2 (fx_1 * _S189.x + intrins_0.z, fy_1 * _S189.y + intrins_0.w);
                break;
            }
            _S164 = all_valid_1 & (!_S177);
            break;
        }
        _S165 = &proj_points_0[int(3)];
        for(;;)
        {
            float _S190 = sigmas_0->p_0[int(3)].z;
            proj_points_0[int(3)] = float2 {sigmas_0->p_0[int(3)].x, sigmas_0->p_0[int(3)].y} / make_float2 (_S190);
            if(_S190 < 0.0f)
            {
                _S177 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S191 = camera_distortion_jac_0(proj_points_0[int(3)], dist_coeffs_10);
                _S177 = !((F32_min((determinant_0(_S191)), ((F32_min((_S191.rows[int(0)].x), (_S191.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S177)
            {
                break;
            }
            float u_7 = proj_points_0[int(3)].x;
            float v_7 = proj_points_0[int(3)].y;
            float r2_7 = u_7 * u_7 + v_7 * v_7;
            float2  _S192 = proj_points_0[int(3)] * make_float2 (1.0f + r2_7 * ((*dist_coeffs_10)[int(0)] + r2_7 * ((*dist_coeffs_10)[int(1)] + r2_7 * ((*dist_coeffs_10)[int(2)] + r2_7 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_7 * v_7 + (*dist_coeffs_10)[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + (*dist_coeffs_10)[int(6)] * r2_7, 2.0f * (*dist_coeffs_10)[int(5)] * u_7 * v_7 + (*dist_coeffs_10)[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + (*dist_coeffs_10)[int(7)] * r2_7);
            float2  _S193 = _S192 + make_float2 ((*dist_coeffs_10)[int(8)] * _S192.x + (*dist_coeffs_10)[int(9)] * _S192.y, 0.0f);
            proj_points_0[int(3)] = make_float2 (fx_1 * _S193.x + intrins_0.z, fy_1 * _S193.y + intrins_0.w);
            break;
        }
        bool all_valid_2 = _S164 & (!_S177);
        _S166 = &proj_points_0[int(4)];
        for(;;)
        {
            float _S194 = sigmas_0->p_0[int(4)].z;
            proj_points_0[int(4)] = float2 {sigmas_0->p_0[int(4)].x, sigmas_0->p_0[int(4)].y} / make_float2 (_S194);
            if(_S194 < 0.0f)
            {
                _S177 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S195 = camera_distortion_jac_0(proj_points_0[int(4)], dist_coeffs_10);
                _S177 = !((F32_min((determinant_0(_S195)), ((F32_min((_S195.rows[int(0)].x), (_S195.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S177)
            {
                break;
            }
            float u_8 = proj_points_0[int(4)].x;
            float v_8 = proj_points_0[int(4)].y;
            float r2_8 = u_8 * u_8 + v_8 * v_8;
            float2  _S196 = proj_points_0[int(4)] * make_float2 (1.0f + r2_8 * ((*dist_coeffs_10)[int(0)] + r2_8 * ((*dist_coeffs_10)[int(1)] + r2_8 * ((*dist_coeffs_10)[int(2)] + r2_8 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_8 * v_8 + (*dist_coeffs_10)[int(5)] * (r2_8 + 2.0f * u_8 * u_8) + (*dist_coeffs_10)[int(6)] * r2_8, 2.0f * (*dist_coeffs_10)[int(5)] * u_8 * v_8 + (*dist_coeffs_10)[int(4)] * (r2_8 + 2.0f * v_8 * v_8) + (*dist_coeffs_10)[int(7)] * r2_8);
            float2  _S197 = _S196 + make_float2 ((*dist_coeffs_10)[int(8)] * _S196.x + (*dist_coeffs_10)[int(9)] * _S196.y, 0.0f);
            proj_points_0[int(4)] = make_float2 (fx_1 * _S197.x + intrins_0.z, fy_1 * _S197.y + intrins_0.w);
            break;
        }
        bool all_valid_3 = all_valid_2 & (!_S177);
        for(;;)
        {
            _S167 = &proj_points_0[int(5)];
            for(;;)
            {
                float _S198 = sigmas_0->p_0[int(5)].z;
                proj_points_0[int(5)] = float2 {sigmas_0->p_0[int(5)].x, sigmas_0->p_0[int(5)].y} / make_float2 (_S198);
                if(_S198 < 0.0f)
                {
                    _S177 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S199 = camera_distortion_jac_0(proj_points_0[int(5)], dist_coeffs_10);
                    _S177 = !((F32_min((determinant_0(_S199)), ((F32_min((_S199.rows[int(0)].x), (_S199.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S177)
                {
                    break;
                }
                float u_9 = proj_points_0[int(5)].x;
                float v_9 = proj_points_0[int(5)].y;
                float r2_9 = u_9 * u_9 + v_9 * v_9;
                float2  _S200 = proj_points_0[int(5)] * make_float2 (1.0f + r2_9 * ((*dist_coeffs_10)[int(0)] + r2_9 * ((*dist_coeffs_10)[int(1)] + r2_9 * ((*dist_coeffs_10)[int(2)] + r2_9 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_9 * v_9 + (*dist_coeffs_10)[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + (*dist_coeffs_10)[int(6)] * r2_9, 2.0f * (*dist_coeffs_10)[int(5)] * u_9 * v_9 + (*dist_coeffs_10)[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + (*dist_coeffs_10)[int(7)] * r2_9);
                float2  _S201 = _S200 + make_float2 ((*dist_coeffs_10)[int(8)] * _S200.x + (*dist_coeffs_10)[int(9)] * _S200.y, 0.0f);
                proj_points_0[int(5)] = make_float2 (fx_1 * _S201.x + intrins_0.z, fy_1 * _S201.y + intrins_0.w);
                break;
            }
            _S168 = all_valid_3 & (!_S177);
            break;
        }
        _S169 = &proj_points_0[int(6)];
        for(;;)
        {
            float _S202 = sigmas_0->p_0[int(6)].z;
            proj_points_0[int(6)] = float2 {sigmas_0->p_0[int(6)].x, sigmas_0->p_0[int(6)].y} / make_float2 (_S202);
            if(_S202 < 0.0f)
            {
                _S177 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S203 = camera_distortion_jac_0(proj_points_0[int(6)], dist_coeffs_10);
                _S177 = !((F32_min((determinant_0(_S203)), ((F32_min((_S203.rows[int(0)].x), (_S203.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S177)
            {
                break;
            }
            float u_10 = proj_points_0[int(6)].x;
            float v_10 = proj_points_0[int(6)].y;
            float r2_10 = u_10 * u_10 + v_10 * v_10;
            float2  _S204 = proj_points_0[int(6)] * make_float2 (1.0f + r2_10 * ((*dist_coeffs_10)[int(0)] + r2_10 * ((*dist_coeffs_10)[int(1)] + r2_10 * ((*dist_coeffs_10)[int(2)] + r2_10 * (*dist_coeffs_10)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_10)[int(4)] * u_10 * v_10 + (*dist_coeffs_10)[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + (*dist_coeffs_10)[int(6)] * r2_10, 2.0f * (*dist_coeffs_10)[int(5)] * u_10 * v_10 + (*dist_coeffs_10)[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + (*dist_coeffs_10)[int(7)] * r2_10);
            float2  _S205 = _S204 + make_float2 ((*dist_coeffs_10)[int(8)] * _S204.x + (*dist_coeffs_10)[int(9)] * _S204.y, 0.0f);
            proj_points_0[int(6)] = make_float2 (fx_1 * _S205.x + intrins_0.z, fy_1 * _S205.y + intrins_0.w);
            break;
        }
        _S170 = _S168 & (!_S177);
        break;
    }
    if(!_S170)
    {
        return false;
    }
    float2  _S206 = *mean2d_1 + make_float2 (sigmas_0->w_mean_0[int(0)]) * *_S161;
    *mean2d_1 = _S206;
    float2  _S207 = _S206 + make_float2 (sigmas_0->w_mean_0[int(1)]) * *_S162;
    *mean2d_1 = _S207;
    float2  _S208 = _S207 + make_float2 (sigmas_0->w_mean_0[int(2)]) * *_S163;
    *mean2d_1 = _S208;
    float2  _S209 = _S208 + make_float2 (sigmas_0->w_mean_0[int(3)]) * *_S165;
    *mean2d_1 = _S209;
    float2  _S210 = _S209 + make_float2 (sigmas_0->w_mean_0[int(4)]) * *_S166;
    *mean2d_1 = _S210;
    float2  _S211 = _S210 + make_float2 (sigmas_0->w_mean_0[int(5)]) * *_S167;
    *mean2d_1 = _S211;
    float2  _S212 = _S211 + make_float2 (sigmas_0->w_mean_0[int(6)]) * *_S169;
    *mean2d_1 = _S212;
    float _S213 = - _S175;
    float _S214 = - _S176;
    float2  _S215 = make_float2 (clamp_0(_S212.x, _S213, lim_x_pos_0), clamp_0(_S212.y, _S214, lim_y_pos_0));
    float2  d_0 = make_float2 (clamp_0((*_S161).x, _S213, lim_x_pos_0), clamp_0((*_S161).y, _S214, lim_y_pos_0)) - _S215;
    float _S216 = d_0.x;
    float _S217 = d_0.y;
    float _S218 = _S216 * _S217;
    Matrix<float, 2, 2>  _S219 = *cov2d_1 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S216 * _S216, _S218, _S218, _S217 * _S217);
    *cov2d_1 = _S219;
    float2  d_1 = make_float2 (clamp_0((*_S162).x, _S213, lim_x_pos_0), clamp_0((*_S162).y, _S214, lim_y_pos_0)) - _S215;
    float _S220 = d_1.x;
    float _S221 = d_1.y;
    float _S222 = _S220 * _S221;
    Matrix<float, 2, 2>  _S223 = _S219 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S220 * _S220, _S222, _S222, _S221 * _S221);
    *cov2d_1 = _S223;
    float2  d_2 = make_float2 (clamp_0((*_S163).x, _S213, lim_x_pos_0), clamp_0((*_S163).y, _S214, lim_y_pos_0)) - _S215;
    float _S224 = d_2.x;
    float _S225 = d_2.y;
    float _S226 = _S224 * _S225;
    Matrix<float, 2, 2>  _S227 = _S223 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S224 * _S224, _S226, _S226, _S225 * _S225);
    *cov2d_1 = _S227;
    float2  d_3 = make_float2 (clamp_0((*_S165).x, _S213, lim_x_pos_0), clamp_0((*_S165).y, _S214, lim_y_pos_0)) - _S215;
    float _S228 = d_3.x;
    float _S229 = d_3.y;
    float _S230 = _S228 * _S229;
    Matrix<float, 2, 2>  _S231 = _S227 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S228 * _S228, _S230, _S230, _S229 * _S229);
    *cov2d_1 = _S231;
    float2  d_4 = make_float2 (clamp_0((*_S166).x, _S213, lim_x_pos_0), clamp_0((*_S166).y, _S214, lim_y_pos_0)) - _S215;
    float _S232 = d_4.x;
    float _S233 = d_4.y;
    float _S234 = _S232 * _S233;
    Matrix<float, 2, 2>  _S235 = _S231 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S232 * _S232, _S234, _S234, _S233 * _S233);
    *cov2d_1 = _S235;
    float2  d_5 = make_float2 (clamp_0((*_S167).x, _S213, lim_x_pos_0), clamp_0((*_S167).y, _S214, lim_y_pos_0)) - _S215;
    float _S236 = d_5.x;
    float _S237 = d_5.y;
    float _S238 = _S236 * _S237;
    Matrix<float, 2, 2>  _S239 = _S235 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S236 * _S236, _S238, _S238, _S237 * _S237);
    *cov2d_1 = _S239;
    float2  d_6 = make_float2 (clamp_0((*_S169).x, _S213, lim_x_pos_0), clamp_0((*_S169).y, _S214, lim_y_pos_0)) - _S215;
    float _S240 = d_6.x;
    float _S241 = d_6.y;
    float _S242 = _S240 * _S241;
    *cov2d_1 = _S239 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S240 * _S240, _S242, _S242, _S241 * _S241);
    return true;
}

inline __device__ bool persp_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_1, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_11, uint width_2, uint height_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    int2  _S243 = make_int2 (int(0));
    float2  _S244 = make_float2 ((float)_S243.x, (float)_S243.y);
    *mean2d_2 = _S244;
    *cov2d_2 = makeMatrix<float, 2, 2> (0.0f);
    float fx_2 = intrins_1.x;
    float fy_2 = intrins_1.y;
    float _S245 = float(width_2);
    float _S246 = float(height_2);
    float _S247 = 0.30000001192092896f * (0.5f * _S245 / fx_2) * fx_2;
    float lim_x_pos_1 = _S245 + _S247;
    float _S248 = 0.30000001192092896f * (0.5f * _S246 / fy_2) * fy_2;
    float lim_y_pos_1 = _S246 + _S248;
    FixedArray<float2 , 7>  proj_points_1;
    float2  _S249 = float2 {sigmas_1->p_0[int(0)].x, sigmas_1->p_0[int(0)].y} / make_float2 (sigmas_1->p_0[int(0)].z);
    float u_11 = _S249.x;
    float v_11 = _S249.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float _S250 = 2.0f * (*dist_coeffs_11)[int(4)];
    float _S251 = 2.0f * (*dist_coeffs_11)[int(5)];
    float2  _S252 = _S249 * make_float2 (1.0f + r2_11 * ((*dist_coeffs_11)[int(0)] + r2_11 * ((*dist_coeffs_11)[int(1)] + r2_11 * ((*dist_coeffs_11)[int(2)] + r2_11 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S250 * u_11 * v_11 + (*dist_coeffs_11)[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + (*dist_coeffs_11)[int(6)] * r2_11, _S251 * u_11 * v_11 + (*dist_coeffs_11)[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + (*dist_coeffs_11)[int(7)] * r2_11);
    float2  _S253 = _S252 + make_float2 ((*dist_coeffs_11)[int(8)] * _S252.x + (*dist_coeffs_11)[int(9)] * _S252.y, 0.0f);
    float cx_1 = intrins_1.z;
    float cy_1 = intrins_1.w;
    float _S254 = fx_2 * _S253.x + cx_1;
    float _S255 = fy_2 * _S253.y + cy_1;
    float2  _S256 = make_float2 (_S254, _S255);
    proj_points_1[int(0)] = _S256;
    float2  _S257 = float2 {sigmas_1->p_0[int(1)].x, sigmas_1->p_0[int(1)].y} / make_float2 (sigmas_1->p_0[int(1)].z);
    float u_12 = _S257.x;
    float v_12 = _S257.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float2  _S258 = _S257 * make_float2 (1.0f + r2_12 * ((*dist_coeffs_11)[int(0)] + r2_12 * ((*dist_coeffs_11)[int(1)] + r2_12 * ((*dist_coeffs_11)[int(2)] + r2_12 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S250 * u_12 * v_12 + (*dist_coeffs_11)[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + (*dist_coeffs_11)[int(6)] * r2_12, _S251 * u_12 * v_12 + (*dist_coeffs_11)[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + (*dist_coeffs_11)[int(7)] * r2_12);
    float2  _S259 = _S258 + make_float2 ((*dist_coeffs_11)[int(8)] * _S258.x + (*dist_coeffs_11)[int(9)] * _S258.y, 0.0f);
    float _S260 = fx_2 * _S259.x + cx_1;
    float _S261 = fy_2 * _S259.y + cy_1;
    float2  _S262 = make_float2 (_S260, _S261);
    proj_points_1[int(1)] = _S262;
    float2  _S263 = float2 {sigmas_1->p_0[int(2)].x, sigmas_1->p_0[int(2)].y} / make_float2 (sigmas_1->p_0[int(2)].z);
    float u_13 = _S263.x;
    float v_13 = _S263.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S264 = _S263 * make_float2 (1.0f + r2_13 * ((*dist_coeffs_11)[int(0)] + r2_13 * ((*dist_coeffs_11)[int(1)] + r2_13 * ((*dist_coeffs_11)[int(2)] + r2_13 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S250 * u_13 * v_13 + (*dist_coeffs_11)[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + (*dist_coeffs_11)[int(6)] * r2_13, _S251 * u_13 * v_13 + (*dist_coeffs_11)[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + (*dist_coeffs_11)[int(7)] * r2_13);
    float2  _S265 = _S264 + make_float2 ((*dist_coeffs_11)[int(8)] * _S264.x + (*dist_coeffs_11)[int(9)] * _S264.y, 0.0f);
    float _S266 = fx_2 * _S265.x + cx_1;
    float _S267 = fy_2 * _S265.y + cy_1;
    float2  _S268 = make_float2 (_S266, _S267);
    proj_points_1[int(2)] = _S268;
    float2  _S269 = float2 {sigmas_1->p_0[int(3)].x, sigmas_1->p_0[int(3)].y} / make_float2 (sigmas_1->p_0[int(3)].z);
    float u_14 = _S269.x;
    float v_14 = _S269.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S270 = _S269 * make_float2 (1.0f + r2_14 * ((*dist_coeffs_11)[int(0)] + r2_14 * ((*dist_coeffs_11)[int(1)] + r2_14 * ((*dist_coeffs_11)[int(2)] + r2_14 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S250 * u_14 * v_14 + (*dist_coeffs_11)[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + (*dist_coeffs_11)[int(6)] * r2_14, _S251 * u_14 * v_14 + (*dist_coeffs_11)[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + (*dist_coeffs_11)[int(7)] * r2_14);
    float2  _S271 = _S270 + make_float2 ((*dist_coeffs_11)[int(8)] * _S270.x + (*dist_coeffs_11)[int(9)] * _S270.y, 0.0f);
    float _S272 = fx_2 * _S271.x + cx_1;
    float _S273 = fy_2 * _S271.y + cy_1;
    float2  _S274 = make_float2 (_S272, _S273);
    proj_points_1[int(3)] = _S274;
    float2  _S275 = float2 {sigmas_1->p_0[int(4)].x, sigmas_1->p_0[int(4)].y} / make_float2 (sigmas_1->p_0[int(4)].z);
    float u_15 = _S275.x;
    float v_15 = _S275.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S276 = _S275 * make_float2 (1.0f + r2_15 * ((*dist_coeffs_11)[int(0)] + r2_15 * ((*dist_coeffs_11)[int(1)] + r2_15 * ((*dist_coeffs_11)[int(2)] + r2_15 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S250 * u_15 * v_15 + (*dist_coeffs_11)[int(5)] * (r2_15 + 2.0f * u_15 * u_15) + (*dist_coeffs_11)[int(6)] * r2_15, _S251 * u_15 * v_15 + (*dist_coeffs_11)[int(4)] * (r2_15 + 2.0f * v_15 * v_15) + (*dist_coeffs_11)[int(7)] * r2_15);
    float2  _S277 = _S276 + make_float2 ((*dist_coeffs_11)[int(8)] * _S276.x + (*dist_coeffs_11)[int(9)] * _S276.y, 0.0f);
    float _S278 = fx_2 * _S277.x + cx_1;
    float _S279 = fy_2 * _S277.y + cy_1;
    float2  _S280 = make_float2 (_S278, _S279);
    proj_points_1[int(4)] = _S280;
    float2  _S281 = float2 {sigmas_1->p_0[int(5)].x, sigmas_1->p_0[int(5)].y} / make_float2 (sigmas_1->p_0[int(5)].z);
    float u_16 = _S281.x;
    float v_16 = _S281.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float2  _S282 = _S281 * make_float2 (1.0f + r2_16 * ((*dist_coeffs_11)[int(0)] + r2_16 * ((*dist_coeffs_11)[int(1)] + r2_16 * ((*dist_coeffs_11)[int(2)] + r2_16 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S250 * u_16 * v_16 + (*dist_coeffs_11)[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + (*dist_coeffs_11)[int(6)] * r2_16, _S251 * u_16 * v_16 + (*dist_coeffs_11)[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + (*dist_coeffs_11)[int(7)] * r2_16);
    float2  _S283 = _S282 + make_float2 ((*dist_coeffs_11)[int(8)] * _S282.x + (*dist_coeffs_11)[int(9)] * _S282.y, 0.0f);
    float _S284 = fx_2 * _S283.x + cx_1;
    float _S285 = fy_2 * _S283.y + cy_1;
    float2  _S286 = make_float2 (_S284, _S285);
    proj_points_1[int(5)] = _S286;
    float2  _S287 = float2 {sigmas_1->p_0[int(6)].x, sigmas_1->p_0[int(6)].y} / make_float2 (sigmas_1->p_0[int(6)].z);
    float u_17 = _S287.x;
    float v_17 = _S287.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float2  _S288 = _S287 * make_float2 (1.0f + r2_17 * ((*dist_coeffs_11)[int(0)] + r2_17 * ((*dist_coeffs_11)[int(1)] + r2_17 * ((*dist_coeffs_11)[int(2)] + r2_17 * (*dist_coeffs_11)[int(3)])))) + make_float2 (_S250 * u_17 * v_17 + (*dist_coeffs_11)[int(5)] * (r2_17 + 2.0f * u_17 * u_17) + (*dist_coeffs_11)[int(6)] * r2_17, _S251 * u_17 * v_17 + (*dist_coeffs_11)[int(4)] * (r2_17 + 2.0f * v_17 * v_17) + (*dist_coeffs_11)[int(7)] * r2_17);
    float2  _S289 = _S288 + make_float2 ((*dist_coeffs_11)[int(8)] * _S288.x + (*dist_coeffs_11)[int(9)] * _S288.y, 0.0f);
    float _S290 = fx_2 * _S289.x + cx_1;
    float _S291 = fy_2 * _S289.y + cy_1;
    float2  _S292 = make_float2 (_S290, _S291);
    proj_points_1[int(6)] = _S292;
    float2  _S293 = *mean2d_2 + make_float2 (sigmas_1->w_mean_0[int(0)]) * _S256;
    *mean2d_2 = _S293;
    float2  _S294 = _S293 + make_float2 (sigmas_1->w_mean_0[int(1)]) * _S262;
    *mean2d_2 = _S294;
    float2  _S295 = _S294 + make_float2 (sigmas_1->w_mean_0[int(2)]) * _S268;
    *mean2d_2 = _S295;
    float2  _S296 = _S295 + make_float2 (sigmas_1->w_mean_0[int(3)]) * _S274;
    *mean2d_2 = _S296;
    float2  _S297 = _S296 + make_float2 (sigmas_1->w_mean_0[int(4)]) * _S280;
    *mean2d_2 = _S297;
    float2  _S298 = _S297 + make_float2 (sigmas_1->w_mean_0[int(5)]) * _S286;
    *mean2d_2 = _S298;
    float2  _S299 = _S298 + make_float2 (sigmas_1->w_mean_0[int(6)]) * _S292;
    *mean2d_2 = _S299;
    float _S300 = - _S247;
    float _S301 = - _S248;
    float2  _S302 = make_float2 (clamp_0(_S299.x, _S300, lim_x_pos_1), clamp_0(_S299.y, _S301, lim_y_pos_1));
    float2  d_7 = make_float2 (clamp_0(_S254, _S300, lim_x_pos_1), clamp_0(_S255, _S301, lim_y_pos_1)) - _S302;
    float _S303 = d_7.x;
    float _S304 = d_7.y;
    float _S305 = _S303 * _S304;
    Matrix<float, 2, 2>  _S306 = *cov2d_2 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S303 * _S303, _S305, _S305, _S304 * _S304);
    *cov2d_2 = _S306;
    float2  d_8 = make_float2 (clamp_0(_S260, _S300, lim_x_pos_1), clamp_0(_S261, _S301, lim_y_pos_1)) - _S302;
    float _S307 = d_8.x;
    float _S308 = d_8.y;
    float _S309 = _S307 * _S308;
    Matrix<float, 2, 2>  _S310 = _S306 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S307 * _S307, _S309, _S309, _S308 * _S308);
    *cov2d_2 = _S310;
    float2  d_9 = make_float2 (clamp_0(_S266, _S300, lim_x_pos_1), clamp_0(_S267, _S301, lim_y_pos_1)) - _S302;
    float _S311 = d_9.x;
    float _S312 = d_9.y;
    float _S313 = _S311 * _S312;
    Matrix<float, 2, 2>  _S314 = _S310 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S311 * _S311, _S313, _S313, _S312 * _S312);
    *cov2d_2 = _S314;
    float2  d_10 = make_float2 (clamp_0(_S272, _S300, lim_x_pos_1), clamp_0(_S273, _S301, lim_y_pos_1)) - _S302;
    float _S315 = d_10.x;
    float _S316 = d_10.y;
    float _S317 = _S315 * _S316;
    Matrix<float, 2, 2>  _S318 = _S314 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S315 * _S315, _S317, _S317, _S316 * _S316);
    *cov2d_2 = _S318;
    float2  d_11 = make_float2 (clamp_0(_S278, _S300, lim_x_pos_1), clamp_0(_S279, _S301, lim_y_pos_1)) - _S302;
    float _S319 = d_11.x;
    float _S320 = d_11.y;
    float _S321 = _S319 * _S320;
    Matrix<float, 2, 2>  _S322 = _S318 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S319 * _S319, _S321, _S321, _S320 * _S320);
    *cov2d_2 = _S322;
    float2  d_12 = make_float2 (clamp_0(_S284, _S300, lim_x_pos_1), clamp_0(_S285, _S301, lim_y_pos_1)) - _S302;
    float _S323 = d_12.x;
    float _S324 = d_12.y;
    float _S325 = _S323 * _S324;
    Matrix<float, 2, 2>  _S326 = _S322 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S323 * _S323, _S325, _S325, _S324 * _S324);
    *cov2d_2 = _S326;
    float2  d_13 = make_float2 (clamp_0(_S290, _S300, lim_x_pos_1), clamp_0(_S291, _S301, lim_y_pos_1)) - _S302;
    float _S327 = d_13.x;
    float _S328 = d_13.y;
    float _S329 = _S327 * _S328;
    *cov2d_2 = _S326 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S327 * _S327, _S329, _S329, _S328 * _S328);
    return true;
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_9, float dOut_12)
{
    DiffPair_float_0 _S330 = *dpx_9;
    float _S331 = - (*dpy_4).primal_0 / ((*dpx_9).primal_0 * (*dpx_9).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_12;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S331;
    float _S332 = _S330.primal_0 / (_S330.primal_0 * _S330.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_12;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S332;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S333, float _S334)
{
    return (F32_atan2((_S333), (_S334)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S335, DiffPair_float_0 * _S336, float _S337)
{
    _d_atan2_0(_S335, _S336, _S337);
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_10, float _s_dOut_2)
{
    float _S338 = (*dpx_10).primal_0.x;
    float _S339 = (*dpx_10).primal_0.y;
    DiffPair_float_0 _S340;
    (&_S340)->primal_0 = _S338 * _S338 + _S339 * _S339;
    (&_S340)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S340, _s_dOut_2);
    float _S341 = (*dpx_10).primal_0.y * _S340.differential_0;
    float _S342 = _S341 + _S341;
    float _S343 = (*dpx_10).primal_0.x * _S340.differential_0;
    float _S344 = _S343 + _S343;
    float2  _S345 = make_float2 (0.0f);
    *&((&_S345)->y) = _S342;
    *&((&_S345)->x) = _S344;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S345;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S346, float _S347)
{
    s_bwd_prop_length_impl_1(_S346, _S347);
    return;
}

inline __device__ bool fisheye_proj_3dgs_0(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_12, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    float k_0;
    float2  _S348;
    float _S349;
    float _S350;
    bool _S351;
    for(;;)
    {
        float2  _S352 = float2 {mean3d_1.x, mean3d_1.y};
        _S348 = _S352;
        float r_8 = length_0(_S352);
        _S349 = r_8;
        float _S353 = mean3d_1.z;
        _S350 = _S353;
        float theta_4 = (F32_atan2((r_8), (_S353)));
        if(theta_4 < 0.00100000004749745f)
        {
            k_0 = (1.0f - theta_4 * theta_4 / 3.0f) / _S353;
        }
        else
        {
            k_0 = theta_4 / r_8;
        }
        float2  _S354 = _S352 * make_float2 (k_0);
        *mean2d_3 = _S354;
        Matrix<float, 2, 2>  _S355 = camera_distortion_jac_0(_S354, dist_coeffs_12);
        bool _S356 = !((F32_min((determinant_0(_S355)), ((F32_min((_S355.rows[int(0)].x), (_S355.rows[int(1)].y)))))) > 0.0f);
        _S351 = _S356;
        if(_S356)
        {
            break;
        }
        float u_18 = (*mean2d_3).x;
        float v_18 = (*mean2d_3).y;
        float r2_18 = u_18 * u_18 + v_18 * v_18;
        float2  _S357 = *mean2d_3 * make_float2 (1.0f + r2_18 * ((*dist_coeffs_12)[int(0)] + r2_18 * ((*dist_coeffs_12)[int(1)] + r2_18 * ((*dist_coeffs_12)[int(2)] + r2_18 * (*dist_coeffs_12)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_12)[int(4)] * u_18 * v_18 + (*dist_coeffs_12)[int(5)] * (r2_18 + 2.0f * u_18 * u_18) + (*dist_coeffs_12)[int(6)] * r2_18, 2.0f * (*dist_coeffs_12)[int(5)] * u_18 * v_18 + (*dist_coeffs_12)[int(4)] * (r2_18 + 2.0f * v_18 * v_18) + (*dist_coeffs_12)[int(7)] * r2_18);
        float2  _S358 = _S357 + make_float2 ((*dist_coeffs_12)[int(8)] * _S357.x + (*dist_coeffs_12)[int(9)] * _S357.y, 0.0f);
        *mean2d_3 = make_float2 (intrins_2.x * _S358.x + intrins_2.z, intrins_2.y * _S358.y + intrins_2.w);
        break;
    }
    if(!!_S351)
    {
        return false;
    }
    Matrix<float, 2, 3>  J_1;
    float2  _S359 = make_float2 (0.0f);
    float2  seed_0 = _S359;
    *&((&seed_0)->x) = 1.0f;
    float2  _S360 = seed_0;
    float _S361 = s_primal_ctx_atan2_0(_S349, _S350);
    bool _S362 = _S361 < 0.00100000004749745f;
    float _S363;
    float _S364;
    float _S365;
    if(_S362)
    {
        float _S366 = 1.0f - _S361 * _S361 / 3.0f;
        float _S367 = _S350 * _S350;
        k_0 = _S366 / _S350;
        _S363 = 0.0f;
        _S364 = _S367;
        _S365 = _S366;
    }
    else
    {
        float _S368 = _S349 * _S349;
        k_0 = _S361 / _S349;
        _S363 = _S368;
        _S364 = 0.0f;
        _S365 = 0.0f;
    }
    float2  _S369 = make_float2 (k_0);
    float2  _S370 = _S348 * make_float2 (k_0);
    float u_19 = _S370.x;
    float v_19 = _S370.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S371 = (*dist_coeffs_12)[int(2)] + r2_19 * (*dist_coeffs_12)[int(3)];
    float _S372 = (*dist_coeffs_12)[int(1)] + r2_19 * _S371;
    float _S373 = (*dist_coeffs_12)[int(0)] + r2_19 * _S372;
    float _S374 = 2.0f * (*dist_coeffs_12)[int(4)];
    float _S375 = 2.0f * (*dist_coeffs_12)[int(5)];
    float fx_3 = intrins_2.x;
    float fy_3 = intrins_2.y;
    float _S376 = fx_3 * _S360.x;
    float2  _S377 = make_float2 (_S376, fy_3 * _S360.y) + make_float2 ((*dist_coeffs_12)[int(8)] * _S376, (*dist_coeffs_12)[int(9)] * _S376);
    float2  _S378 = _S370 * _S377;
    float _S379 = (*dist_coeffs_12)[int(4)] * _S377.y;
    float _S380 = (*dist_coeffs_12)[int(5)] * _S377.x;
    float _S381 = _S378.x + _S378.y;
    float _S382 = r2_19 * _S381;
    float _S383 = r2_19 * _S382;
    float _S384 = (*dist_coeffs_12)[int(7)] * _S377.y + _S379 + (*dist_coeffs_12)[int(6)] * _S377.x + _S380 + _S373 * _S381 + _S372 * _S382 + _S371 * _S383 + (*dist_coeffs_12)[int(3)] * (r2_19 * _S383);
    float _S385 = v_19 * _S384;
    float _S386 = u_19 * _S384;
    float2  _S387 = make_float2 (1.0f + r2_19 * _S373) * _S377 + make_float2 (_S375 * (v_19 * _S377.y) + 2.0f * u_19 * _S380 + 2.0f * (u_19 * _S380) + _S374 * (v_19 * _S377.x) + _S386 + _S386, 2.0f * v_19 * _S379 + 2.0f * (v_19 * _S379) + _S375 * u_19 * _S377.y + _S374 * u_19 * _S377.x + _S385 + _S385);
    float2  _S388 = _S348 * _S387;
    float2  _S389 = _S369 * _S387;
    float _S390 = _S388.x + _S388.y;
    if(_S362)
    {
        float _S391 = _S390 / _S364;
        float _S392 = _S365 * - _S391;
        float _S393 = _S361 * (0.3333333432674408f * - (_S350 * _S391));
        k_0 = _S393 + _S393;
        _S363 = _S392;
        _S364 = 0.0f;
    }
    else
    {
        float _S394 = _S390 / _S363;
        float _S395 = _S361 * - _S394;
        k_0 = _S349 * _S394;
        _S363 = 0.0f;
        _S364 = _S395;
    }
    float _S396 = _S349;
    DiffPair_float_0 _S397;
    (&_S397)->primal_0 = _S349;
    (&_S397)->differential_0 = 0.0f;
    float _S398 = _S350;
    DiffPair_float_0 _S399;
    (&_S399)->primal_0 = _S350;
    (&_S399)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S397, &_S399, k_0);
    float _S400 = _S399.differential_0 + _S363;
    float _S401 = _S397.differential_0 + _S364;
    float2  _S402 = make_float2 (0.0f);
    float2  _S403 = _S348;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S404;
    (&_S404)->primal_0 = _S348;
    (&_S404)->differential_0 = _S402;
    s_bwd_length_impl_1(&_S404, _S401);
    float2  _S405 = _S404.differential_0 + _S389;
    float3  _S406 = make_float3 (_S405.x, _S405.y, _S400);
    J_1[int(0)] = _S406;
    float2  seed_1 = _S359;
    *&((&seed_1)->y) = 1.0f;
    float2  _S407 = seed_1;
    if(_S362)
    {
        float _S408 = 1.0f - _S361 * _S361 / 3.0f;
        float _S409 = _S350 * _S350;
        k_0 = _S408 / _S350;
        _S363 = 0.0f;
        _S364 = _S409;
        _S365 = _S408;
    }
    else
    {
        float _S410 = _S349 * _S349;
        k_0 = _S361 / _S349;
        _S363 = _S410;
        _S364 = 0.0f;
        _S365 = 0.0f;
    }
    float2  _S411 = make_float2 (k_0);
    float2  _S412 = _S348 * make_float2 (k_0);
    float u_20 = _S412.x;
    float v_20 = _S412.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S413 = (*dist_coeffs_12)[int(2)] + r2_20 * (*dist_coeffs_12)[int(3)];
    float _S414 = (*dist_coeffs_12)[int(1)] + r2_20 * _S413;
    float _S415 = (*dist_coeffs_12)[int(0)] + r2_20 * _S414;
    float _S416 = fx_3 * _S407.x;
    float2  _S417 = make_float2 (_S416, fy_3 * _S407.y) + make_float2 ((*dist_coeffs_12)[int(8)] * _S416, (*dist_coeffs_12)[int(9)] * _S416);
    float2  _S418 = _S412 * _S417;
    float _S419 = (*dist_coeffs_12)[int(4)] * _S417.y;
    float _S420 = (*dist_coeffs_12)[int(5)] * _S417.x;
    float _S421 = _S418.x + _S418.y;
    float _S422 = r2_20 * _S421;
    float _S423 = r2_20 * _S422;
    float _S424 = (*dist_coeffs_12)[int(7)] * _S417.y + _S419 + (*dist_coeffs_12)[int(6)] * _S417.x + _S420 + _S415 * _S421 + _S414 * _S422 + _S413 * _S423 + (*dist_coeffs_12)[int(3)] * (r2_20 * _S423);
    float _S425 = v_20 * _S424;
    float _S426 = u_20 * _S424;
    float2  _S427 = make_float2 (1.0f + r2_20 * _S415) * _S417 + make_float2 (_S375 * (v_20 * _S417.y) + 2.0f * u_20 * _S420 + 2.0f * (u_20 * _S420) + _S374 * (v_20 * _S417.x) + _S426 + _S426, 2.0f * v_20 * _S419 + 2.0f * (v_20 * _S419) + _S375 * u_20 * _S417.y + _S374 * u_20 * _S417.x + _S425 + _S425);
    float2  _S428 = _S348 * _S427;
    float2  _S429 = _S411 * _S427;
    float _S430 = _S428.x + _S428.y;
    if(_S362)
    {
        float _S431 = _S430 / _S364;
        float _S432 = _S365 * - _S431;
        float _S433 = _S361 * (0.3333333432674408f * - (_S350 * _S431));
        k_0 = _S433 + _S433;
        _S363 = _S432;
        _S364 = 0.0f;
    }
    else
    {
        float _S434 = _S430 / _S363;
        float _S435 = _S361 * - _S434;
        k_0 = _S349 * _S434;
        _S363 = 0.0f;
        _S364 = _S435;
    }
    DiffPair_float_0 _S436;
    (&_S436)->primal_0 = _S396;
    (&_S436)->differential_0 = 0.0f;
    DiffPair_float_0 _S437;
    (&_S437)->primal_0 = _S398;
    (&_S437)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S436, &_S437, k_0);
    float _S438 = _S437.differential_0 + _S363;
    float _S439 = _S436.differential_0 + _S364;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S440;
    (&_S440)->primal_0 = _S403;
    (&_S440)->differential_0 = _S402;
    s_bwd_length_impl_1(&_S440, _S439);
    float2  _S441 = _S440.differential_0 + _S429;
    float3  _S442 = make_float3 (_S441.x, _S441.y, _S438);
    J_1[int(1)] = _S442;
    *cov2d_3 = mul_6(mul_5(J_1, cov3d_1), transpose_1(J_1));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_1(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_13, Matrix<float, 2, 2>  * cov2d_4, float2  * mean2d_4)
{
    float2  _S443 = float2 {mean3d_2.x, mean3d_2.y};
    float r_9 = length_0(_S443);
    float _S444 = mean3d_2.z;
    float theta_5 = (F32_atan2((r_9), (_S444)));
    float k_1;
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S444;
    }
    else
    {
        k_1 = theta_5 / r_9;
    }
    float2  _S445 = _S443 * make_float2 (k_1);
    float u_21 = _S445.x;
    float v_21 = _S445.y;
    float r2_21 = u_21 * u_21 + v_21 * v_21;
    float _S446 = 2.0f * (*dist_coeffs_13)[int(4)];
    float _S447 = 2.0f * (*dist_coeffs_13)[int(5)];
    float2  _S448 = _S445 * make_float2 (1.0f + r2_21 * ((*dist_coeffs_13)[int(0)] + r2_21 * ((*dist_coeffs_13)[int(1)] + r2_21 * ((*dist_coeffs_13)[int(2)] + r2_21 * (*dist_coeffs_13)[int(3)])))) + make_float2 (_S446 * u_21 * v_21 + (*dist_coeffs_13)[int(5)] * (r2_21 + 2.0f * u_21 * u_21) + (*dist_coeffs_13)[int(6)] * r2_21, _S447 * u_21 * v_21 + (*dist_coeffs_13)[int(4)] * (r2_21 + 2.0f * v_21 * v_21) + (*dist_coeffs_13)[int(7)] * r2_21);
    float2  _S449 = _S448 + make_float2 ((*dist_coeffs_13)[int(8)] * _S448.x + (*dist_coeffs_13)[int(9)] * _S448.y, 0.0f);
    float fx_4 = intrins_3.x;
    float fy_4 = intrins_3.y;
    *mean2d_4 = make_float2 (fx_4 * _S449.x + intrins_3.z, fy_4 * _S449.y + intrins_3.w);
    Matrix<float, 2, 3>  J_2;
    float2  _S450 = make_float2 (0.0f);
    float2  seed_2 = _S450;
    *&((&seed_2)->x) = 1.0f;
    float2  _S451 = seed_2;
    float _S452 = s_primal_ctx_atan2_0(r_9, _S444);
    bool _S453 = _S452 < 0.00100000004749745f;
    float _S454;
    float _S455;
    float _S456;
    if(_S453)
    {
        float _S457 = 1.0f - _S452 * _S452 / 3.0f;
        float _S458 = _S444 * _S444;
        k_1 = _S457 / _S444;
        _S454 = 0.0f;
        _S455 = _S458;
        _S456 = _S457;
    }
    else
    {
        float _S459 = r_9 * r_9;
        k_1 = _S452 / r_9;
        _S454 = _S459;
        _S455 = 0.0f;
        _S456 = 0.0f;
    }
    float2  _S460 = make_float2 (k_1);
    float2  _S461 = _S443 * make_float2 (k_1);
    float u_22 = _S461.x;
    float v_22 = _S461.y;
    float r2_22 = u_22 * u_22 + v_22 * v_22;
    float _S462 = (*dist_coeffs_13)[int(2)] + r2_22 * (*dist_coeffs_13)[int(3)];
    float _S463 = (*dist_coeffs_13)[int(1)] + r2_22 * _S462;
    float _S464 = (*dist_coeffs_13)[int(0)] + r2_22 * _S463;
    float _S465 = fx_4 * _S451.x;
    float2  _S466 = make_float2 (_S465, fy_4 * _S451.y) + make_float2 ((*dist_coeffs_13)[int(8)] * _S465, (*dist_coeffs_13)[int(9)] * _S465);
    float2  _S467 = _S461 * _S466;
    float _S468 = (*dist_coeffs_13)[int(4)] * _S466.y;
    float _S469 = (*dist_coeffs_13)[int(5)] * _S466.x;
    float _S470 = _S467.x + _S467.y;
    float _S471 = r2_22 * _S470;
    float _S472 = r2_22 * _S471;
    float _S473 = (*dist_coeffs_13)[int(7)] * _S466.y + _S468 + (*dist_coeffs_13)[int(6)] * _S466.x + _S469 + _S464 * _S470 + _S463 * _S471 + _S462 * _S472 + (*dist_coeffs_13)[int(3)] * (r2_22 * _S472);
    float _S474 = v_22 * _S473;
    float _S475 = u_22 * _S473;
    float2  _S476 = make_float2 (1.0f + r2_22 * _S464) * _S466 + make_float2 (_S447 * (v_22 * _S466.y) + 2.0f * u_22 * _S469 + 2.0f * (u_22 * _S469) + _S446 * (v_22 * _S466.x) + _S475 + _S475, 2.0f * v_22 * _S468 + 2.0f * (v_22 * _S468) + _S447 * u_22 * _S466.y + _S446 * u_22 * _S466.x + _S474 + _S474);
    float2  _S477 = _S443 * _S476;
    float2  _S478 = _S460 * _S476;
    float _S479 = _S477.x + _S477.y;
    if(_S453)
    {
        float _S480 = _S479 / _S455;
        float _S481 = _S456 * - _S480;
        float _S482 = _S452 * (0.3333333432674408f * - (_S444 * _S480));
        k_1 = _S482 + _S482;
        _S454 = _S481;
        _S455 = 0.0f;
    }
    else
    {
        float _S483 = _S479 / _S454;
        float _S484 = _S452 * - _S483;
        k_1 = r_9 * _S483;
        _S454 = 0.0f;
        _S455 = _S484;
    }
    DiffPair_float_0 _S485;
    (&_S485)->primal_0 = r_9;
    (&_S485)->differential_0 = 0.0f;
    DiffPair_float_0 _S486;
    (&_S486)->primal_0 = _S444;
    (&_S486)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S485, &_S486, k_1);
    float _S487 = _S486.differential_0 + _S454;
    float _S488 = _S485.differential_0 + _S455;
    float2  _S489 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S490;
    (&_S490)->primal_0 = _S443;
    (&_S490)->differential_0 = _S489;
    s_bwd_length_impl_1(&_S490, _S488);
    float2  _S491 = _S490.differential_0 + _S478;
    float3  _S492 = make_float3 (_S491.x, _S491.y, _S487);
    J_2[int(0)] = _S492;
    float2  seed_3 = _S450;
    *&((&seed_3)->y) = 1.0f;
    float2  _S493 = seed_3;
    if(_S453)
    {
        float _S494 = 1.0f - _S452 * _S452 / 3.0f;
        float _S495 = _S444 * _S444;
        k_1 = _S494 / _S444;
        _S454 = 0.0f;
        _S455 = _S495;
        _S456 = _S494;
    }
    else
    {
        float _S496 = r_9 * r_9;
        k_1 = _S452 / r_9;
        _S454 = _S496;
        _S455 = 0.0f;
        _S456 = 0.0f;
    }
    float2  _S497 = make_float2 (k_1);
    float2  _S498 = _S443 * make_float2 (k_1);
    float u_23 = _S498.x;
    float v_23 = _S498.y;
    float r2_23 = u_23 * u_23 + v_23 * v_23;
    float _S499 = (*dist_coeffs_13)[int(2)] + r2_23 * (*dist_coeffs_13)[int(3)];
    float _S500 = (*dist_coeffs_13)[int(1)] + r2_23 * _S499;
    float _S501 = (*dist_coeffs_13)[int(0)] + r2_23 * _S500;
    float _S502 = fx_4 * _S493.x;
    float2  _S503 = make_float2 (_S502, fy_4 * _S493.y) + make_float2 ((*dist_coeffs_13)[int(8)] * _S502, (*dist_coeffs_13)[int(9)] * _S502);
    float2  _S504 = _S498 * _S503;
    float _S505 = (*dist_coeffs_13)[int(4)] * _S503.y;
    float _S506 = (*dist_coeffs_13)[int(5)] * _S503.x;
    float _S507 = _S504.x + _S504.y;
    float _S508 = r2_23 * _S507;
    float _S509 = r2_23 * _S508;
    float _S510 = (*dist_coeffs_13)[int(7)] * _S503.y + _S505 + (*dist_coeffs_13)[int(6)] * _S503.x + _S506 + _S501 * _S507 + _S500 * _S508 + _S499 * _S509 + (*dist_coeffs_13)[int(3)] * (r2_23 * _S509);
    float _S511 = v_23 * _S510;
    float _S512 = u_23 * _S510;
    float2  _S513 = make_float2 (1.0f + r2_23 * _S501) * _S503 + make_float2 (_S447 * (v_23 * _S503.y) + 2.0f * u_23 * _S506 + 2.0f * (u_23 * _S506) + _S446 * (v_23 * _S503.x) + _S512 + _S512, 2.0f * v_23 * _S505 + 2.0f * (v_23 * _S505) + _S447 * u_23 * _S503.y + _S446 * u_23 * _S503.x + _S511 + _S511);
    float2  _S514 = _S443 * _S513;
    float2  _S515 = _S497 * _S513;
    float _S516 = _S514.x + _S514.y;
    if(_S453)
    {
        float _S517 = _S516 / _S455;
        float _S518 = _S456 * - _S517;
        float _S519 = _S452 * (0.3333333432674408f * - (_S444 * _S517));
        k_1 = _S519 + _S519;
        _S454 = _S518;
        _S455 = 0.0f;
    }
    else
    {
        float _S520 = _S516 / _S454;
        float _S521 = _S452 * - _S520;
        k_1 = r_9 * _S520;
        _S454 = 0.0f;
        _S455 = _S521;
    }
    DiffPair_float_0 _S522;
    (&_S522)->primal_0 = r_9;
    (&_S522)->differential_0 = 0.0f;
    DiffPair_float_0 _S523;
    (&_S523)->primal_0 = _S444;
    (&_S523)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S522, &_S523, k_1);
    float _S524 = _S523.differential_0 + _S454;
    float _S525 = _S522.differential_0 + _S455;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S526;
    (&_S526)->primal_0 = _S443;
    (&_S526)->differential_0 = _S489;
    s_bwd_length_impl_1(&_S526, _S525);
    float2  _S527 = _S526.differential_0 + _S515;
    float3  _S528 = make_float3 (_S527.x, _S527.y, _S524);
    J_2[int(1)] = _S528;
    *cov2d_4 = mul_6(mul_5(J_2, cov3d_2), transpose_1(J_2));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_2, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_14, Matrix<float, 2, 2>  * cov2d_5, float2  * mean2d_5)
{
    float2  * _S529;
    bool _S530;
    float2  * _S531;
    bool _S532;
    float2  * _S533;
    bool _S534;
    bool _S535;
    float2  * _S536;
    bool _S537;
    float2  * _S538;
    bool _S539;
    float2  * _S540;
    bool _S541;
    bool _S542;
    float2  * _S543;
    bool _S544;
    bool _S545;
    int2  _S546 = make_int2 (int(0));
    float2  _S547 = make_float2 ((float)_S546.x, (float)_S546.y);
    *mean2d_5 = _S547;
    *cov2d_5 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_2;
    for(;;)
    {
        float k_2;
        _S529 = &proj_points_2[int(0)];
        for(;;)
        {
            float2  _S548 = float2 {sigmas_2->p_0[int(0)].x, sigmas_2->p_0[int(0)].y};
            float r_10 = length_0(_S548);
            float _S549 = sigmas_2->p_0[int(0)].z;
            float theta_6 = (F32_atan2((r_10), (_S549)));
            if(theta_6 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_6 * theta_6 / 3.0f) / _S549;
            }
            else
            {
                k_2 = theta_6 / r_10;
            }
            float2  _S550 = _S548 * make_float2 (k_2);
            proj_points_2[int(0)] = _S550;
            Matrix<float, 2, 2>  _S551 = camera_distortion_jac_0(_S550, dist_coeffs_14);
            bool _S552 = !((F32_min((determinant_0(_S551)), ((F32_min((_S551.rows[int(0)].x), (_S551.rows[int(1)].y)))))) > 0.0f);
            _S530 = _S552;
            if(_S552)
            {
                break;
            }
            float u_24 = proj_points_2[int(0)].x;
            float v_24 = proj_points_2[int(0)].y;
            float r2_24 = u_24 * u_24 + v_24 * v_24;
            float2  _S553 = proj_points_2[int(0)] * make_float2 (1.0f + r2_24 * ((*dist_coeffs_14)[int(0)] + r2_24 * ((*dist_coeffs_14)[int(1)] + r2_24 * ((*dist_coeffs_14)[int(2)] + r2_24 * (*dist_coeffs_14)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_14)[int(4)] * u_24 * v_24 + (*dist_coeffs_14)[int(5)] * (r2_24 + 2.0f * u_24 * u_24) + (*dist_coeffs_14)[int(6)] * r2_24, 2.0f * (*dist_coeffs_14)[int(5)] * u_24 * v_24 + (*dist_coeffs_14)[int(4)] * (r2_24 + 2.0f * v_24 * v_24) + (*dist_coeffs_14)[int(7)] * r2_24);
            float2  _S554 = _S553 + make_float2 ((*dist_coeffs_14)[int(8)] * _S553.x + (*dist_coeffs_14)[int(9)] * _S553.y, 0.0f);
            proj_points_2[int(0)] = make_float2 (intrins_4.x * _S554.x + intrins_4.z, intrins_4.y * _S554.y + intrins_4.w);
            break;
        }
        bool all_valid_4 = true & (!_S530);
        _S531 = &proj_points_2[int(1)];
        for(;;)
        {
            float2  _S555 = float2 {sigmas_2->p_0[int(1)].x, sigmas_2->p_0[int(1)].y};
            float r_11 = length_0(_S555);
            float _S556 = sigmas_2->p_0[int(1)].z;
            float theta_7 = (F32_atan2((r_11), (_S556)));
            if(theta_7 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_7 * theta_7 / 3.0f) / _S556;
            }
            else
            {
                k_2 = theta_7 / r_11;
            }
            float2  _S557 = _S555 * make_float2 (k_2);
            proj_points_2[int(1)] = _S557;
            Matrix<float, 2, 2>  _S558 = camera_distortion_jac_0(_S557, dist_coeffs_14);
            bool _S559 = !((F32_min((determinant_0(_S558)), ((F32_min((_S558.rows[int(0)].x), (_S558.rows[int(1)].y)))))) > 0.0f);
            _S532 = _S559;
            if(_S559)
            {
                break;
            }
            float u_25 = proj_points_2[int(1)].x;
            float v_25 = proj_points_2[int(1)].y;
            float r2_25 = u_25 * u_25 + v_25 * v_25;
            float2  _S560 = proj_points_2[int(1)] * make_float2 (1.0f + r2_25 * ((*dist_coeffs_14)[int(0)] + r2_25 * ((*dist_coeffs_14)[int(1)] + r2_25 * ((*dist_coeffs_14)[int(2)] + r2_25 * (*dist_coeffs_14)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_14)[int(4)] * u_25 * v_25 + (*dist_coeffs_14)[int(5)] * (r2_25 + 2.0f * u_25 * u_25) + (*dist_coeffs_14)[int(6)] * r2_25, 2.0f * (*dist_coeffs_14)[int(5)] * u_25 * v_25 + (*dist_coeffs_14)[int(4)] * (r2_25 + 2.0f * v_25 * v_25) + (*dist_coeffs_14)[int(7)] * r2_25);
            float2  _S561 = _S560 + make_float2 ((*dist_coeffs_14)[int(8)] * _S560.x + (*dist_coeffs_14)[int(9)] * _S560.y, 0.0f);
            proj_points_2[int(1)] = make_float2 (intrins_4.x * _S561.x + intrins_4.z, intrins_4.y * _S561.y + intrins_4.w);
            break;
        }
        bool all_valid_5 = all_valid_4 & (!_S532);
        for(;;)
        {
            _S533 = &proj_points_2[int(2)];
            for(;;)
            {
                float2  _S562 = float2 {sigmas_2->p_0[int(2)].x, sigmas_2->p_0[int(2)].y};
                float r_12 = length_0(_S562);
                float _S563 = sigmas_2->p_0[int(2)].z;
                float theta_8 = (F32_atan2((r_12), (_S563)));
                if(theta_8 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_8 * theta_8 / 3.0f) / _S563;
                }
                else
                {
                    k_2 = theta_8 / r_12;
                }
                float2  _S564 = _S562 * make_float2 (k_2);
                proj_points_2[int(2)] = _S564;
                Matrix<float, 2, 2>  _S565 = camera_distortion_jac_0(_S564, dist_coeffs_14);
                bool _S566 = !((F32_min((determinant_0(_S565)), ((F32_min((_S565.rows[int(0)].x), (_S565.rows[int(1)].y)))))) > 0.0f);
                _S534 = _S566;
                if(_S566)
                {
                    break;
                }
                float u_26 = proj_points_2[int(2)].x;
                float v_26 = proj_points_2[int(2)].y;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float2  _S567 = proj_points_2[int(2)] * make_float2 (1.0f + r2_26 * ((*dist_coeffs_14)[int(0)] + r2_26 * ((*dist_coeffs_14)[int(1)] + r2_26 * ((*dist_coeffs_14)[int(2)] + r2_26 * (*dist_coeffs_14)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_14)[int(4)] * u_26 * v_26 + (*dist_coeffs_14)[int(5)] * (r2_26 + 2.0f * u_26 * u_26) + (*dist_coeffs_14)[int(6)] * r2_26, 2.0f * (*dist_coeffs_14)[int(5)] * u_26 * v_26 + (*dist_coeffs_14)[int(4)] * (r2_26 + 2.0f * v_26 * v_26) + (*dist_coeffs_14)[int(7)] * r2_26);
                float2  _S568 = _S567 + make_float2 ((*dist_coeffs_14)[int(8)] * _S567.x + (*dist_coeffs_14)[int(9)] * _S567.y, 0.0f);
                proj_points_2[int(2)] = make_float2 (intrins_4.x * _S568.x + intrins_4.z, intrins_4.y * _S568.y + intrins_4.w);
                break;
            }
            _S535 = all_valid_5 & (!_S534);
            break;
        }
        _S536 = &proj_points_2[int(3)];
        for(;;)
        {
            float2  _S569 = float2 {sigmas_2->p_0[int(3)].x, sigmas_2->p_0[int(3)].y};
            float r_13 = length_0(_S569);
            float _S570 = sigmas_2->p_0[int(3)].z;
            float theta_9 = (F32_atan2((r_13), (_S570)));
            if(theta_9 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_9 * theta_9 / 3.0f) / _S570;
            }
            else
            {
                k_2 = theta_9 / r_13;
            }
            float2  _S571 = _S569 * make_float2 (k_2);
            proj_points_2[int(3)] = _S571;
            Matrix<float, 2, 2>  _S572 = camera_distortion_jac_0(_S571, dist_coeffs_14);
            bool _S573 = !((F32_min((determinant_0(_S572)), ((F32_min((_S572.rows[int(0)].x), (_S572.rows[int(1)].y)))))) > 0.0f);
            _S537 = _S573;
            if(_S573)
            {
                break;
            }
            float u_27 = proj_points_2[int(3)].x;
            float v_27 = proj_points_2[int(3)].y;
            float r2_27 = u_27 * u_27 + v_27 * v_27;
            float2  _S574 = proj_points_2[int(3)] * make_float2 (1.0f + r2_27 * ((*dist_coeffs_14)[int(0)] + r2_27 * ((*dist_coeffs_14)[int(1)] + r2_27 * ((*dist_coeffs_14)[int(2)] + r2_27 * (*dist_coeffs_14)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_14)[int(4)] * u_27 * v_27 + (*dist_coeffs_14)[int(5)] * (r2_27 + 2.0f * u_27 * u_27) + (*dist_coeffs_14)[int(6)] * r2_27, 2.0f * (*dist_coeffs_14)[int(5)] * u_27 * v_27 + (*dist_coeffs_14)[int(4)] * (r2_27 + 2.0f * v_27 * v_27) + (*dist_coeffs_14)[int(7)] * r2_27);
            float2  _S575 = _S574 + make_float2 ((*dist_coeffs_14)[int(8)] * _S574.x + (*dist_coeffs_14)[int(9)] * _S574.y, 0.0f);
            proj_points_2[int(3)] = make_float2 (intrins_4.x * _S575.x + intrins_4.z, intrins_4.y * _S575.y + intrins_4.w);
            break;
        }
        bool all_valid_6 = _S535 & (!_S537);
        _S538 = &proj_points_2[int(4)];
        for(;;)
        {
            float2  _S576 = float2 {sigmas_2->p_0[int(4)].x, sigmas_2->p_0[int(4)].y};
            float r_14 = length_0(_S576);
            float _S577 = sigmas_2->p_0[int(4)].z;
            float theta_10 = (F32_atan2((r_14), (_S577)));
            if(theta_10 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_10 * theta_10 / 3.0f) / _S577;
            }
            else
            {
                k_2 = theta_10 / r_14;
            }
            float2  _S578 = _S576 * make_float2 (k_2);
            proj_points_2[int(4)] = _S578;
            Matrix<float, 2, 2>  _S579 = camera_distortion_jac_0(_S578, dist_coeffs_14);
            bool _S580 = !((F32_min((determinant_0(_S579)), ((F32_min((_S579.rows[int(0)].x), (_S579.rows[int(1)].y)))))) > 0.0f);
            _S539 = _S580;
            if(_S580)
            {
                break;
            }
            float u_28 = proj_points_2[int(4)].x;
            float v_28 = proj_points_2[int(4)].y;
            float r2_28 = u_28 * u_28 + v_28 * v_28;
            float2  _S581 = proj_points_2[int(4)] * make_float2 (1.0f + r2_28 * ((*dist_coeffs_14)[int(0)] + r2_28 * ((*dist_coeffs_14)[int(1)] + r2_28 * ((*dist_coeffs_14)[int(2)] + r2_28 * (*dist_coeffs_14)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_14)[int(4)] * u_28 * v_28 + (*dist_coeffs_14)[int(5)] * (r2_28 + 2.0f * u_28 * u_28) + (*dist_coeffs_14)[int(6)] * r2_28, 2.0f * (*dist_coeffs_14)[int(5)] * u_28 * v_28 + (*dist_coeffs_14)[int(4)] * (r2_28 + 2.0f * v_28 * v_28) + (*dist_coeffs_14)[int(7)] * r2_28);
            float2  _S582 = _S581 + make_float2 ((*dist_coeffs_14)[int(8)] * _S581.x + (*dist_coeffs_14)[int(9)] * _S581.y, 0.0f);
            proj_points_2[int(4)] = make_float2 (intrins_4.x * _S582.x + intrins_4.z, intrins_4.y * _S582.y + intrins_4.w);
            break;
        }
        bool all_valid_7 = all_valid_6 & (!_S539);
        for(;;)
        {
            _S540 = &proj_points_2[int(5)];
            for(;;)
            {
                float2  _S583 = float2 {sigmas_2->p_0[int(5)].x, sigmas_2->p_0[int(5)].y};
                float r_15 = length_0(_S583);
                float _S584 = sigmas_2->p_0[int(5)].z;
                float theta_11 = (F32_atan2((r_15), (_S584)));
                if(theta_11 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_11 * theta_11 / 3.0f) / _S584;
                }
                else
                {
                    k_2 = theta_11 / r_15;
                }
                float2  _S585 = _S583 * make_float2 (k_2);
                proj_points_2[int(5)] = _S585;
                Matrix<float, 2, 2>  _S586 = camera_distortion_jac_0(_S585, dist_coeffs_14);
                bool _S587 = !((F32_min((determinant_0(_S586)), ((F32_min((_S586.rows[int(0)].x), (_S586.rows[int(1)].y)))))) > 0.0f);
                _S541 = _S587;
                if(_S587)
                {
                    break;
                }
                float u_29 = proj_points_2[int(5)].x;
                float v_29 = proj_points_2[int(5)].y;
                float r2_29 = u_29 * u_29 + v_29 * v_29;
                float2  _S588 = proj_points_2[int(5)] * make_float2 (1.0f + r2_29 * ((*dist_coeffs_14)[int(0)] + r2_29 * ((*dist_coeffs_14)[int(1)] + r2_29 * ((*dist_coeffs_14)[int(2)] + r2_29 * (*dist_coeffs_14)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_14)[int(4)] * u_29 * v_29 + (*dist_coeffs_14)[int(5)] * (r2_29 + 2.0f * u_29 * u_29) + (*dist_coeffs_14)[int(6)] * r2_29, 2.0f * (*dist_coeffs_14)[int(5)] * u_29 * v_29 + (*dist_coeffs_14)[int(4)] * (r2_29 + 2.0f * v_29 * v_29) + (*dist_coeffs_14)[int(7)] * r2_29);
                float2  _S589 = _S588 + make_float2 ((*dist_coeffs_14)[int(8)] * _S588.x + (*dist_coeffs_14)[int(9)] * _S588.y, 0.0f);
                proj_points_2[int(5)] = make_float2 (intrins_4.x * _S589.x + intrins_4.z, intrins_4.y * _S589.y + intrins_4.w);
                break;
            }
            _S542 = all_valid_7 & (!_S541);
            break;
        }
        _S543 = &proj_points_2[int(6)];
        for(;;)
        {
            float2  _S590 = float2 {sigmas_2->p_0[int(6)].x, sigmas_2->p_0[int(6)].y};
            float r_16 = length_0(_S590);
            float _S591 = sigmas_2->p_0[int(6)].z;
            float theta_12 = (F32_atan2((r_16), (_S591)));
            if(theta_12 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_12 * theta_12 / 3.0f) / _S591;
            }
            else
            {
                k_2 = theta_12 / r_16;
            }
            float2  _S592 = _S590 * make_float2 (k_2);
            proj_points_2[int(6)] = _S592;
            Matrix<float, 2, 2>  _S593 = camera_distortion_jac_0(_S592, dist_coeffs_14);
            bool _S594 = !((F32_min((determinant_0(_S593)), ((F32_min((_S593.rows[int(0)].x), (_S593.rows[int(1)].y)))))) > 0.0f);
            _S544 = _S594;
            if(_S594)
            {
                break;
            }
            float u_30 = proj_points_2[int(6)].x;
            float v_30 = proj_points_2[int(6)].y;
            float r2_30 = u_30 * u_30 + v_30 * v_30;
            float2  _S595 = proj_points_2[int(6)] * make_float2 (1.0f + r2_30 * ((*dist_coeffs_14)[int(0)] + r2_30 * ((*dist_coeffs_14)[int(1)] + r2_30 * ((*dist_coeffs_14)[int(2)] + r2_30 * (*dist_coeffs_14)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_14)[int(4)] * u_30 * v_30 + (*dist_coeffs_14)[int(5)] * (r2_30 + 2.0f * u_30 * u_30) + (*dist_coeffs_14)[int(6)] * r2_30, 2.0f * (*dist_coeffs_14)[int(5)] * u_30 * v_30 + (*dist_coeffs_14)[int(4)] * (r2_30 + 2.0f * v_30 * v_30) + (*dist_coeffs_14)[int(7)] * r2_30);
            float2  _S596 = _S595 + make_float2 ((*dist_coeffs_14)[int(8)] * _S595.x + (*dist_coeffs_14)[int(9)] * _S595.y, 0.0f);
            proj_points_2[int(6)] = make_float2 (intrins_4.x * _S596.x + intrins_4.z, intrins_4.y * _S596.y + intrins_4.w);
            break;
        }
        _S545 = _S542 & (!_S544);
        break;
    }
    if(!_S545)
    {
        return false;
    }
    float2  _S597 = *mean2d_5 + make_float2 (sigmas_2->w_mean_0[int(0)]) * *_S529;
    *mean2d_5 = _S597;
    float2  _S598 = _S597 + make_float2 (sigmas_2->w_mean_0[int(1)]) * *_S531;
    *mean2d_5 = _S598;
    float2  _S599 = _S598 + make_float2 (sigmas_2->w_mean_0[int(2)]) * *_S533;
    *mean2d_5 = _S599;
    float2  _S600 = _S599 + make_float2 (sigmas_2->w_mean_0[int(3)]) * *_S536;
    *mean2d_5 = _S600;
    float2  _S601 = _S600 + make_float2 (sigmas_2->w_mean_0[int(4)]) * *_S538;
    *mean2d_5 = _S601;
    float2  _S602 = _S601 + make_float2 (sigmas_2->w_mean_0[int(5)]) * *_S540;
    *mean2d_5 = _S602;
    float2  _S603 = _S602 + make_float2 (sigmas_2->w_mean_0[int(6)]) * *_S543;
    *mean2d_5 = _S603;
    float2  d_14 = *_S529 - _S603;
    float _S604 = d_14.x;
    float _S605 = d_14.y;
    float _S606 = _S604 * _S605;
    Matrix<float, 2, 2>  _S607 = *cov2d_5 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S604 * _S604, _S606, _S606, _S605 * _S605);
    *cov2d_5 = _S607;
    float2  d_15 = *_S531 - *mean2d_5;
    float _S608 = d_15.x;
    float _S609 = d_15.y;
    float _S610 = _S608 * _S609;
    Matrix<float, 2, 2>  _S611 = _S607 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S608 * _S608, _S610, _S610, _S609 * _S609);
    *cov2d_5 = _S611;
    float2  d_16 = *_S533 - *mean2d_5;
    float _S612 = d_16.x;
    float _S613 = d_16.y;
    float _S614 = _S612 * _S613;
    Matrix<float, 2, 2>  _S615 = _S611 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S612 * _S612, _S614, _S614, _S613 * _S613);
    *cov2d_5 = _S615;
    float2  d_17 = *_S536 - *mean2d_5;
    float _S616 = d_17.x;
    float _S617 = d_17.y;
    float _S618 = _S616 * _S617;
    Matrix<float, 2, 2>  _S619 = _S615 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S616 * _S616, _S618, _S618, _S617 * _S617);
    *cov2d_5 = _S619;
    float2  d_18 = *_S538 - *mean2d_5;
    float _S620 = d_18.x;
    float _S621 = d_18.y;
    float _S622 = _S620 * _S621;
    Matrix<float, 2, 2>  _S623 = _S619 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S620 * _S620, _S622, _S622, _S621 * _S621);
    *cov2d_5 = _S623;
    float2  d_19 = *_S540 - *mean2d_5;
    float _S624 = d_19.x;
    float _S625 = d_19.y;
    float _S626 = _S624 * _S625;
    Matrix<float, 2, 2>  _S627 = _S623 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S624 * _S624, _S626, _S626, _S625 * _S625);
    *cov2d_5 = _S627;
    float2  d_20 = *_S543 - *mean2d_5;
    float _S628 = d_20.x;
    float _S629 = d_20.y;
    float _S630 = _S628 * _S629;
    *cov2d_5 = _S627 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S628 * _S628, _S630, _S630, _S629 * _S629);
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_3, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_15, Matrix<float, 2, 2>  * cov2d_6, float2  * mean2d_6)
{
    int2  _S631 = make_int2 (int(0));
    float2  _S632 = make_float2 ((float)_S631.x, (float)_S631.y);
    *mean2d_6 = _S632;
    *cov2d_6 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_3;
    float2  _S633 = float2 {sigmas_3->p_0[int(0)].x, sigmas_3->p_0[int(0)].y};
    float r_17 = length_0(_S633);
    float _S634 = sigmas_3->p_0[int(0)].z;
    float theta_13 = (F32_atan2((r_17), (_S634)));
    float k_3;
    if(theta_13 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_13 * theta_13 / 3.0f) / _S634;
    }
    else
    {
        k_3 = theta_13 / r_17;
    }
    float2  _S635 = _S633 * make_float2 (k_3);
    float u_31 = _S635.x;
    float v_31 = _S635.y;
    float r2_31 = u_31 * u_31 + v_31 * v_31;
    float _S636 = 2.0f * (*dist_coeffs_15)[int(4)];
    float _S637 = 2.0f * (*dist_coeffs_15)[int(5)];
    float2  _S638 = _S635 * make_float2 (1.0f + r2_31 * ((*dist_coeffs_15)[int(0)] + r2_31 * ((*dist_coeffs_15)[int(1)] + r2_31 * ((*dist_coeffs_15)[int(2)] + r2_31 * (*dist_coeffs_15)[int(3)])))) + make_float2 (_S636 * u_31 * v_31 + (*dist_coeffs_15)[int(5)] * (r2_31 + 2.0f * u_31 * u_31) + (*dist_coeffs_15)[int(6)] * r2_31, _S637 * u_31 * v_31 + (*dist_coeffs_15)[int(4)] * (r2_31 + 2.0f * v_31 * v_31) + (*dist_coeffs_15)[int(7)] * r2_31);
    float2  _S639 = _S638 + make_float2 ((*dist_coeffs_15)[int(8)] * _S638.x + (*dist_coeffs_15)[int(9)] * _S638.y, 0.0f);
    float fx_5 = intrins_5.x;
    float fy_5 = intrins_5.y;
    float cx_2 = intrins_5.z;
    float cy_2 = intrins_5.w;
    proj_points_3[int(0)] = make_float2 (fx_5 * _S639.x + cx_2, fy_5 * _S639.y + cy_2);
    float2  _S640 = float2 {sigmas_3->p_0[int(1)].x, sigmas_3->p_0[int(1)].y};
    float r_18 = length_0(_S640);
    float _S641 = sigmas_3->p_0[int(1)].z;
    float theta_14 = (F32_atan2((r_18), (_S641)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_14 * theta_14 / 3.0f) / _S641;
    }
    else
    {
        k_3 = theta_14 / r_18;
    }
    float2  _S642 = _S640 * make_float2 (k_3);
    float u_32 = _S642.x;
    float v_32 = _S642.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float2  _S643 = _S642 * make_float2 (1.0f + r2_32 * ((*dist_coeffs_15)[int(0)] + r2_32 * ((*dist_coeffs_15)[int(1)] + r2_32 * ((*dist_coeffs_15)[int(2)] + r2_32 * (*dist_coeffs_15)[int(3)])))) + make_float2 (_S636 * u_32 * v_32 + (*dist_coeffs_15)[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + (*dist_coeffs_15)[int(6)] * r2_32, _S637 * u_32 * v_32 + (*dist_coeffs_15)[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + (*dist_coeffs_15)[int(7)] * r2_32);
    float2  _S644 = _S643 + make_float2 ((*dist_coeffs_15)[int(8)] * _S643.x + (*dist_coeffs_15)[int(9)] * _S643.y, 0.0f);
    proj_points_3[int(1)] = make_float2 (fx_5 * _S644.x + cx_2, fy_5 * _S644.y + cy_2);
    float2  _S645 = float2 {sigmas_3->p_0[int(2)].x, sigmas_3->p_0[int(2)].y};
    float r_19 = length_0(_S645);
    float _S646 = sigmas_3->p_0[int(2)].z;
    float theta_15 = (F32_atan2((r_19), (_S646)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_15 * theta_15 / 3.0f) / _S646;
    }
    else
    {
        k_3 = theta_15 / r_19;
    }
    float2  _S647 = _S645 * make_float2 (k_3);
    float u_33 = _S647.x;
    float v_33 = _S647.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S648 = _S647 * make_float2 (1.0f + r2_33 * ((*dist_coeffs_15)[int(0)] + r2_33 * ((*dist_coeffs_15)[int(1)] + r2_33 * ((*dist_coeffs_15)[int(2)] + r2_33 * (*dist_coeffs_15)[int(3)])))) + make_float2 (_S636 * u_33 * v_33 + (*dist_coeffs_15)[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + (*dist_coeffs_15)[int(6)] * r2_33, _S637 * u_33 * v_33 + (*dist_coeffs_15)[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + (*dist_coeffs_15)[int(7)] * r2_33);
    float2  _S649 = _S648 + make_float2 ((*dist_coeffs_15)[int(8)] * _S648.x + (*dist_coeffs_15)[int(9)] * _S648.y, 0.0f);
    proj_points_3[int(2)] = make_float2 (fx_5 * _S649.x + cx_2, fy_5 * _S649.y + cy_2);
    float2  _S650 = float2 {sigmas_3->p_0[int(3)].x, sigmas_3->p_0[int(3)].y};
    float r_20 = length_0(_S650);
    float _S651 = sigmas_3->p_0[int(3)].z;
    float theta_16 = (F32_atan2((r_20), (_S651)));
    if(theta_16 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_16 * theta_16 / 3.0f) / _S651;
    }
    else
    {
        k_3 = theta_16 / r_20;
    }
    float2  _S652 = _S650 * make_float2 (k_3);
    float u_34 = _S652.x;
    float v_34 = _S652.y;
    float r2_34 = u_34 * u_34 + v_34 * v_34;
    float2  _S653 = _S652 * make_float2 (1.0f + r2_34 * ((*dist_coeffs_15)[int(0)] + r2_34 * ((*dist_coeffs_15)[int(1)] + r2_34 * ((*dist_coeffs_15)[int(2)] + r2_34 * (*dist_coeffs_15)[int(3)])))) + make_float2 (_S636 * u_34 * v_34 + (*dist_coeffs_15)[int(5)] * (r2_34 + 2.0f * u_34 * u_34) + (*dist_coeffs_15)[int(6)] * r2_34, _S637 * u_34 * v_34 + (*dist_coeffs_15)[int(4)] * (r2_34 + 2.0f * v_34 * v_34) + (*dist_coeffs_15)[int(7)] * r2_34);
    float2  _S654 = _S653 + make_float2 ((*dist_coeffs_15)[int(8)] * _S653.x + (*dist_coeffs_15)[int(9)] * _S653.y, 0.0f);
    proj_points_3[int(3)] = make_float2 (fx_5 * _S654.x + cx_2, fy_5 * _S654.y + cy_2);
    float2  _S655 = float2 {sigmas_3->p_0[int(4)].x, sigmas_3->p_0[int(4)].y};
    float r_21 = length_0(_S655);
    float _S656 = sigmas_3->p_0[int(4)].z;
    float theta_17 = (F32_atan2((r_21), (_S656)));
    if(theta_17 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_17 * theta_17 / 3.0f) / _S656;
    }
    else
    {
        k_3 = theta_17 / r_21;
    }
    float2  _S657 = _S655 * make_float2 (k_3);
    float u_35 = _S657.x;
    float v_35 = _S657.y;
    float r2_35 = u_35 * u_35 + v_35 * v_35;
    float2  _S658 = _S657 * make_float2 (1.0f + r2_35 * ((*dist_coeffs_15)[int(0)] + r2_35 * ((*dist_coeffs_15)[int(1)] + r2_35 * ((*dist_coeffs_15)[int(2)] + r2_35 * (*dist_coeffs_15)[int(3)])))) + make_float2 (_S636 * u_35 * v_35 + (*dist_coeffs_15)[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + (*dist_coeffs_15)[int(6)] * r2_35, _S637 * u_35 * v_35 + (*dist_coeffs_15)[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + (*dist_coeffs_15)[int(7)] * r2_35);
    float2  _S659 = _S658 + make_float2 ((*dist_coeffs_15)[int(8)] * _S658.x + (*dist_coeffs_15)[int(9)] * _S658.y, 0.0f);
    proj_points_3[int(4)] = make_float2 (fx_5 * _S659.x + cx_2, fy_5 * _S659.y + cy_2);
    float2  _S660 = float2 {sigmas_3->p_0[int(5)].x, sigmas_3->p_0[int(5)].y};
    float r_22 = length_0(_S660);
    float _S661 = sigmas_3->p_0[int(5)].z;
    float theta_18 = (F32_atan2((r_22), (_S661)));
    if(theta_18 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_18 * theta_18 / 3.0f) / _S661;
    }
    else
    {
        k_3 = theta_18 / r_22;
    }
    float2  _S662 = _S660 * make_float2 (k_3);
    float u_36 = _S662.x;
    float v_36 = _S662.y;
    float r2_36 = u_36 * u_36 + v_36 * v_36;
    float2  _S663 = _S662 * make_float2 (1.0f + r2_36 * ((*dist_coeffs_15)[int(0)] + r2_36 * ((*dist_coeffs_15)[int(1)] + r2_36 * ((*dist_coeffs_15)[int(2)] + r2_36 * (*dist_coeffs_15)[int(3)])))) + make_float2 (_S636 * u_36 * v_36 + (*dist_coeffs_15)[int(5)] * (r2_36 + 2.0f * u_36 * u_36) + (*dist_coeffs_15)[int(6)] * r2_36, _S637 * u_36 * v_36 + (*dist_coeffs_15)[int(4)] * (r2_36 + 2.0f * v_36 * v_36) + (*dist_coeffs_15)[int(7)] * r2_36);
    float2  _S664 = _S663 + make_float2 ((*dist_coeffs_15)[int(8)] * _S663.x + (*dist_coeffs_15)[int(9)] * _S663.y, 0.0f);
    proj_points_3[int(5)] = make_float2 (fx_5 * _S664.x + cx_2, fy_5 * _S664.y + cy_2);
    float2  _S665 = float2 {sigmas_3->p_0[int(6)].x, sigmas_3->p_0[int(6)].y};
    float r_23 = length_0(_S665);
    float _S666 = sigmas_3->p_0[int(6)].z;
    float theta_19 = (F32_atan2((r_23), (_S666)));
    if(theta_19 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_19 * theta_19 / 3.0f) / _S666;
    }
    else
    {
        k_3 = theta_19 / r_23;
    }
    float2  _S667 = _S665 * make_float2 (k_3);
    float u_37 = _S667.x;
    float v_37 = _S667.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float2  _S668 = _S667 * make_float2 (1.0f + r2_37 * ((*dist_coeffs_15)[int(0)] + r2_37 * ((*dist_coeffs_15)[int(1)] + r2_37 * ((*dist_coeffs_15)[int(2)] + r2_37 * (*dist_coeffs_15)[int(3)])))) + make_float2 (_S636 * u_37 * v_37 + (*dist_coeffs_15)[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + (*dist_coeffs_15)[int(6)] * r2_37, _S637 * u_37 * v_37 + (*dist_coeffs_15)[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + (*dist_coeffs_15)[int(7)] * r2_37);
    float2  _S669 = _S668 + make_float2 ((*dist_coeffs_15)[int(8)] * _S668.x + (*dist_coeffs_15)[int(9)] * _S668.y, 0.0f);
    float2  _S670 = make_float2 (fx_5 * _S669.x + cx_2, fy_5 * _S669.y + cy_2);
    proj_points_3[int(6)] = _S670;
    float2  _S671 = *mean2d_6 + make_float2 (sigmas_3->w_mean_0[int(0)]) * proj_points_3[int(0)];
    *mean2d_6 = _S671;
    float2  _S672 = _S671 + make_float2 (sigmas_3->w_mean_0[int(1)]) * proj_points_3[int(1)];
    *mean2d_6 = _S672;
    float2  _S673 = _S672 + make_float2 (sigmas_3->w_mean_0[int(2)]) * proj_points_3[int(2)];
    *mean2d_6 = _S673;
    float2  _S674 = _S673 + make_float2 (sigmas_3->w_mean_0[int(3)]) * proj_points_3[int(3)];
    *mean2d_6 = _S674;
    float2  _S675 = _S674 + make_float2 (sigmas_3->w_mean_0[int(4)]) * proj_points_3[int(4)];
    *mean2d_6 = _S675;
    float2  _S676 = _S675 + make_float2 (sigmas_3->w_mean_0[int(5)]) * proj_points_3[int(5)];
    *mean2d_6 = _S676;
    float2  _S677 = _S676 + make_float2 (sigmas_3->w_mean_0[int(6)]) * _S670;
    *mean2d_6 = _S677;
    float2  d_21 = proj_points_3[int(0)] - _S677;
    float _S678 = d_21.x;
    float _S679 = d_21.y;
    float _S680 = _S678 * _S679;
    Matrix<float, 2, 2>  _S681 = *cov2d_6 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S678 * _S678, _S680, _S680, _S679 * _S679);
    *cov2d_6 = _S681;
    float2  d_22 = proj_points_3[int(1)] - *mean2d_6;
    float _S682 = d_22.x;
    float _S683 = d_22.y;
    float _S684 = _S682 * _S683;
    Matrix<float, 2, 2>  _S685 = _S681 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S682 * _S682, _S684, _S684, _S683 * _S683);
    *cov2d_6 = _S685;
    float2  d_23 = proj_points_3[int(2)] - *mean2d_6;
    float _S686 = d_23.x;
    float _S687 = d_23.y;
    float _S688 = _S686 * _S687;
    Matrix<float, 2, 2>  _S689 = _S685 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S686 * _S686, _S688, _S688, _S687 * _S687);
    *cov2d_6 = _S689;
    float2  d_24 = proj_points_3[int(3)] - *mean2d_6;
    float _S690 = d_24.x;
    float _S691 = d_24.y;
    float _S692 = _S690 * _S691;
    Matrix<float, 2, 2>  _S693 = _S689 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S690 * _S690, _S692, _S692, _S691 * _S691);
    *cov2d_6 = _S693;
    float2  d_25 = proj_points_3[int(4)] - *mean2d_6;
    float _S694 = d_25.x;
    float _S695 = d_25.y;
    float _S696 = _S694 * _S695;
    Matrix<float, 2, 2>  _S697 = _S693 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S694 * _S694, _S696, _S696, _S695 * _S695);
    *cov2d_6 = _S697;
    float2  d_26 = proj_points_3[int(5)] - *mean2d_6;
    float _S698 = d_26.x;
    float _S699 = d_26.y;
    float _S700 = _S698 * _S699;
    Matrix<float, 2, 2>  _S701 = _S697 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S698 * _S698, _S700, _S700, _S699 * _S699);
    *cov2d_6 = _S701;
    float2  d_27 = _S670 - *mean2d_6;
    float _S702 = d_27.x;
    float _S703 = d_27.y;
    float _S704 = _S702 * _S703;
    *cov2d_6 = _S701 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S702 * _S702, _S704, _S704, _S703 * _S703);
    return true;
}

inline __device__ void ortho_proj_3dgs(float3  mean3d_3, Matrix<float, 3, 3>  cov3d_3, float fx_6, float fy_6, float cx_3, float cy_3, Matrix<float, 2, 2>  * cov2d_7, float2  * mean2d_7)
{
    Matrix<float, 2, 3>  J_3 = makeMatrix<float, 2, 3> (fx_6, 0.0f, 0.0f, 0.0f, fy_6, 0.0f);
    *cov2d_7 = mul_6(mul_5(J_3, cov3d_3), transpose_1(J_3));
    *mean2d_7 = make_float2 (fx_6 * mean3d_3.x + cx_3, fy_6 * mean3d_3.y + cy_3);
    return;
}

inline __device__ float add_blur(float eps2d_0, Matrix<float, 2, 2>  * covar_1, float * compensation_0)
{
    float det_orig_0 = *&((covar_1->rows + (int(0)))->x) * *&((covar_1->rows + (int(1)))->y) - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    float _S705 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S705;
    float _S706 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S706;
    float det_blur_0 = _S705 * _S706 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
    *compensation_0 = (F32_sqrt(((F32_max((0.0f), (det_orig_0 / det_blur_0))))));
    return det_blur_0;
}

inline __device__ float3  exp_0(float3  x_14)
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
        *_slang_vector_get_element_ptr(&result_12, i_8) = (F32_exp((_slang_vector_get_element(x_14, i_8))));
        i_8 = i_8 + int(1);
    }
    return result_12;
}

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, float3  dOut_13)
{
    float3  _S707 = exp_0((*dpx_11).primal_0) * dOut_13;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = _S707;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_12, float dOut_14)
{
    float _S708 = 1.0f / (*dpx_12).primal_0 * dOut_14;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S708;
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

inline __device__ float3  max_0(float3  x_15, float3  y_2)
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
        *_slang_vector_get_element_ptr(&result_13, i_9) = (F32_max((_slang_vector_get_element(x_15, i_9)), (_slang_vector_get_element(y_2, i_9))));
        i_9 = i_9 + int(1);
    }
    return result_13;
}

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_4, float3  t_3, float fx_7, float fy_7, float cx_4, float cy_4, FixedArray<float, 10>  * dist_coeffs_16, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_8, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_4, mean_1) + t_3;
        float _S709 = mean_c_0.z;
        bool _S710;
        if(_S709 < near_plane_0)
        {
            _S710 = true;
        }
        else
        {
            _S710 = _S709 > far_plane_0;
        }
        if(_S710)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_0;
        float3  _S711 = exp_0(scale_3);
        float x_16 = quat_4.y;
        float x2_4 = x_16 * x_16;
        float y2_4 = quat_4.z * quat_4.z;
        float z2_4 = quat_4.w * quat_4.w;
        float xy_4 = quat_4.y * quat_4.z;
        float xz_4 = quat_4.y * quat_4.w;
        float yz_4 = quat_4.z * quat_4.w;
        float wx_4 = quat_4.x * quat_4.y;
        float wy_4 = quat_4.x * quat_4.z;
        float wz_4 = quat_4.x * quat_4.w;
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S711.x, 0.0f, 0.0f, 0.0f, _S711.y, 0.0f, 0.0f, 0.0f, _S711.z));
        Matrix<float, 3, 3>  _S712 = transpose_0(R_4);
        float _S713 = float(image_width_0);
        float _S714 = float(image_height_0);
        float _S715 = 0.30000001192092896f * (0.5f * _S713 / fx_7);
        float _S716 = 0.30000001192092896f * (0.5f * _S714 / fy_7);
        float rz_1 = 1.0f / mean_c_0.z;
        float rz2_1 = rz_1 * rz_1;
        Matrix<float, 2, 3>  J_4 = makeMatrix<float, 2, 3> (fx_7 * rz_1, 0.0f, - fx_7 * (mean_c_0.z * (F32_min(((_S713 - cx_4) / fx_7 + _S715), ((F32_max((- (cx_4 / fx_7 + _S715)), (mean_c_0.x * rz_1))))))) * rz2_1, 0.0f, fy_7 * rz_1, - fy_7 * (mean_c_0.z * (F32_min(((_S714 - cy_4) / fy_7 + _S716), ((F32_max((- (cy_4 / fy_7 + _S716)), (mean_c_0.y * rz_1))))))) * rz2_1);
        covar2d_0 = mul_6(mul_5(J_4, mul_4(mul_4(R_4, mul_4(M_2, transpose_0(M_2))), _S712)), transpose_1(J_4));
        *mean2d_8 = make_float2 (fx_7 * mean_c_0.x * rz_1 + cx_4, fy_7 * mean_c_0.y * rz_1 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S717 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S717;
        float _S718 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S718;
        float det_blur_1 = _S717 * _S718 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S719 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S710 = true;
        }
        else
        {
            _S710 = xmin_0 >= _S713;
        }
        if(_S710)
        {
            _S710 = true;
        }
        else
        {
            _S710 = ymax_0 <= 0.0f;
        }
        if(_S710)
        {
            _S710 = true;
        }
        else
        {
            _S710 = ymin_0 >= _S714;
        }
        if(_S710)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = 0.5f * (F32_log((dot_0(mean_c_0, mean_c_0) + 9.99999997475242708e-07f)));
        *conic_0 = make_float3 (_S719.rows[int(0)].x, _S719.rows[int(0)].y, _S719.rows[int(1)].y);
        float3  _S720 = mean_1 - - mul_0(_S712, t_3);
        float3  _S721 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S721;
        float _S722 = _S720.x;
        float _S723 = _S720.y;
        float _S724 = _S720.z;
        float norm_0 = (F32_sqrt((_S722 * _S722 + _S723 * _S723 + _S724 * _S724)));
        float x_17 = _S722 / norm_0;
        float y_3 = _S723 / norm_0;
        float z_0 = _S724 / norm_0;
        float3  _S725 = _S721 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_17) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S725;
        float z2_5 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_17 * x_17 - y_3 * y_3;
        float fS1_0 = 2.0f * x_17 * y_3;
        float3  _S726 = _S725 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_17) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S726;
        float fTmp0C_0 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgb_0 = max_0(_S726 + (make_float3 (-0.59004360437393188f * (x_17 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_17) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_5, float3  t_4, float fx_8, float fy_8, float cx_5, float cy_5, FixedArray<float, 10>  * dist_coeffs_17, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_9, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_5, mean_2) + t_4;
        float _S727 = length_1(mean_c_1);
        bool _S728;
        if(_S727 < near_plane_1)
        {
            _S728 = true;
        }
        else
        {
            _S728 = _S727 > far_plane_1;
        }
        if(_S728)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_6 = make_float4 (fx_8, fy_8, cx_5, cy_5);
        Matrix<float, 2, 2>  covar2d_1;
        float3  _S729 = exp_0(scale_4);
        float x_18 = quat_5.y;
        float x2_5 = x_18 * x_18;
        float y2_5 = quat_5.z * quat_5.z;
        float z2_6 = quat_5.w * quat_5.w;
        float xy_5 = quat_5.y * quat_5.z;
        float xz_5 = quat_5.y * quat_5.w;
        float yz_5 = quat_5.z * quat_5.w;
        float wx_5 = quat_5.x * quat_5.y;
        float wy_5 = quat_5.x * quat_5.z;
        float wz_5 = quat_5.x * quat_5.w;
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_6), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_6), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S729.x, 0.0f, 0.0f, 0.0f, _S729.y, 0.0f, 0.0f, 0.0f, _S729.z));
        Matrix<float, 3, 3>  _S730 = transpose_0(R_5);
        bool _S731 = fisheye_proj_3dgs_0(mean_c_1, mul_4(mul_4(R_5, mul_4(M_3, transpose_0(M_3))), _S730), intrins_6, dist_coeffs_17, &covar2d_1, mean2d_9);
        if(!(true & _S731))
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S732 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S732;
        float _S733 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S733;
        float det_blur_2 = _S732 * _S733 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S734 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S728 = true;
        }
        else
        {
            _S728 = xmin_1 >= float(image_width_1);
        }
        if(_S728)
        {
            _S728 = true;
        }
        else
        {
            _S728 = ymax_1 <= 0.0f;
        }
        if(_S728)
        {
            _S728 = true;
        }
        else
        {
            _S728 = ymin_1 >= float(image_height_1);
        }
        if(_S728)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = 0.5f * (F32_log((dot_0(mean_c_1, mean_c_1) + 9.99999997475242708e-07f)));
        *conic_1 = make_float3 (_S734.rows[int(0)].x, _S734.rows[int(0)].y, _S734.rows[int(1)].y);
        float3  _S735 = mean_2 - - mul_0(_S730, t_4);
        float3  _S736 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S736;
        float _S737 = _S735.x;
        float _S738 = _S735.y;
        float _S739 = _S735.z;
        float norm_1 = (F32_sqrt((_S737 * _S737 + _S738 * _S738 + _S739 * _S739)));
        float x_19 = _S737 / norm_1;
        float y_4 = _S738 / norm_1;
        float z_1 = _S739 / norm_1;
        float3  _S740 = _S736 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_1)[int(1)] + make_float3 (z_1) * (*sh_coeffs_1)[int(2)] - make_float3 (x_19) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S740;
        float z2_7 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_19 * x_19 - y_4 * y_4;
        float fS1_1 = 2.0f * x_19 * y_4;
        float3  _S741 = _S740 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_4) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_19) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S741;
        float fTmp0C_1 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgb_1 = max_0(_S741 + (make_float3 (-0.59004360437393188f * (x_19 * fS1_1 + y_4 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_4) * (*sh_coeffs_1)[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_19) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_19 * fC1_1 - y_4 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_6, float3  t_5, float fx_9, float fy_9, float cx_6, float cy_6, FixedArray<float, 10>  * dist_coeffs_18, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_10, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_6, mean_3) + t_5;
        float _S742 = mean_c_2.z;
        bool _S743;
        if(_S742 < near_plane_2)
        {
            _S743 = true;
        }
        else
        {
            _S743 = _S742 > far_plane_2;
        }
        if(_S743)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_2;
        float3  _S744 = exp_0(scale_5);
        float x_20 = quat_6.y;
        float x2_6 = x_20 * x_20;
        float y2_6 = quat_6.z * quat_6.z;
        float z2_8 = quat_6.w * quat_6.w;
        float xy_6 = quat_6.y * quat_6.z;
        float xz_6 = quat_6.y * quat_6.w;
        float yz_6 = quat_6.z * quat_6.w;
        float wx_6 = quat_6.x * quat_6.y;
        float wy_6 = quat_6.x * quat_6.z;
        float wz_6 = quat_6.x * quat_6.w;
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_8), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_8), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S744.x, 0.0f, 0.0f, 0.0f, _S744.y, 0.0f, 0.0f, 0.0f, _S744.z));
        Matrix<float, 3, 3>  _S745 = transpose_0(R_6);
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
        covar2d_2 = mul_6(mul_5(J_5, mul_4(mul_4(R_6, mul_4(M_4, transpose_0(M_4))), _S745)), transpose_1(J_5));
        *mean2d_10 = make_float2 (fx_9 * mean_c_2.x + cx_6, fy_9 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S746 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S746;
        float _S747 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S747;
        float det_blur_3 = _S746 * _S747 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S748 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S743 = true;
        }
        else
        {
            _S743 = xmin_2 >= float(image_width_2);
        }
        if(_S743)
        {
            _S743 = true;
        }
        else
        {
            _S743 = ymax_2 <= 0.0f;
        }
        if(_S743)
        {
            _S743 = true;
        }
        else
        {
            _S743 = ymin_2 >= float(image_height_2);
        }
        if(_S743)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = 0.5f * (F32_log((dot_0(mean_c_2, mean_c_2) + 9.99999997475242708e-07f)));
        *conic_2 = make_float3 (_S748.rows[int(0)].x, _S748.rows[int(0)].y, _S748.rows[int(1)].y);
        float3  _S749 = mean_3 - - mul_0(_S745, t_5);
        float3  _S750 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S750;
        float _S751 = _S749.x;
        float _S752 = _S749.y;
        float _S753 = _S749.z;
        float norm_2 = (F32_sqrt((_S751 * _S751 + _S752 * _S752 + _S753 * _S753)));
        float x_21 = _S751 / norm_2;
        float y_5 = _S752 / norm_2;
        float z_2 = _S753 / norm_2;
        float3  _S754 = _S750 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_2)[int(1)] + make_float3 (z_2) * (*sh_coeffs_2)[int(2)] - make_float3 (x_21) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S754;
        float z2_9 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_21 * x_21 - y_5 * y_5;
        float fS1_2 = 2.0f * x_21 * y_5;
        float3  _S755 = _S754 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_5) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_9 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_21) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S755;
        float fTmp0C_2 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        *rgb_2 = max_0(_S755 + (make_float3 (-0.59004360437393188f * (x_21 * fS1_2 + y_5 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_5) * (*sh_coeffs_2)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_9 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_21) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_21 * fC1_2 - y_5 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_persp(bool antialiased_3, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_7, float3  t_6, float fx_10, float fy_10, float cx_7, float cy_7, FixedArray<float, 10>  * dist_coeffs_19, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_11, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_0(R_7, mean_4) + t_6;
        float _S756 = mean_c_3.z;
        bool _S757;
        if(_S756 < near_plane_3)
        {
            _S757 = true;
        }
        else
        {
            _S757 = _S756 > far_plane_3;
        }
        if(_S757)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_7 = make_float4 (fx_10, fy_10, cx_7, cy_7);
        Matrix<float, 2, 2>  covar2d_3;
        float3  _S758 = exp_0(scale_6);
        float x_22 = quat_7.y;
        float x2_7 = x_22 * x_22;
        float y2_7 = quat_7.z * quat_7.z;
        float z2_10 = quat_7.w * quat_7.w;
        float xy_7 = quat_7.y * quat_7.z;
        float xz_7 = quat_7.y * quat_7.w;
        float yz_7 = quat_7.z * quat_7.w;
        float wx_7 = quat_7.x * quat_7.y;
        float wy_7 = quat_7.x * quat_7.z;
        float wz_7 = quat_7.x * quat_7.w;
        Matrix<float, 3, 3>  _S759 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_10), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_10), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))));
        SigmaPoints_0 ret_0;
        (&ret_0)->p_0[int(0)] = mean_4;
        (&ret_0)->w_mean_0[int(0)] = 0.0f;
        (&ret_0)->w_cov_0[int(0)] = 2.0f;
        float _S760 = (F32_sqrt((3.0f)));
        float3  delta_0 = make_float3 (_S760 * _S758.x) * _S759.rows[0U];
        float3  _S761 = mean_4 + delta_0;
        float3  _S762 = mean_4 - delta_0;
        float3  delta_1 = make_float3 (_S760 * _S758.y) * _S759.rows[1U];
        float3  _S763 = mean_4 + delta_1;
        float3  _S764 = mean_4 - delta_1;
        float3  delta_2 = make_float3 (_S760 * _S758.z) * _S759.rows[2U];
        float3  _S765 = mean_4 + delta_2;
        float3  _S766 = mean_4 - delta_2;
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
        (&ret_0)->p_0[1U] = mul_0(R_7, _S761) + t_6;
        (&ret_0)->p_0[2U] = mul_0(R_7, _S763) + t_6;
        (&ret_0)->p_0[3U] = mul_0(R_7, _S765) + t_6;
        (&ret_0)->p_0[4U] = mul_0(R_7, _S762) + t_6;
        (&ret_0)->p_0[5U] = mul_0(R_7, _S764) + t_6;
        (&ret_0)->p_0[6U] = mul_0(R_7, _S766) + t_6;
        SigmaPoints_0 _S767 = ret_0;
        bool _S768 = persp_proj_3dgs_ut_0(&_S767, intrins_7, dist_coeffs_19, image_width_3, image_height_3, &covar2d_3, mean2d_11);
        if(!(true & _S768))
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S769 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S769;
        float _S770 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S770;
        float det_blur_4 = _S769 * _S770 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
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
            _S757 = true;
        }
        else
        {
            _S757 = xmin_3 >= float(image_width_3);
        }
        if(_S757)
        {
            _S757 = true;
        }
        else
        {
            _S757 = ymax_3 <= 0.0f;
        }
        if(_S757)
        {
            _S757 = true;
        }
        else
        {
            _S757 = ymin_3 >= float(image_height_3);
        }
        if(_S757)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_3 = make_int4 (int(xmin_3), int(ymin_3), int(xmax_3), int(ymax_3));
        *depth_3 = 0.5f * (F32_log((dot_0(mean_c_3, mean_c_3) + 9.99999997475242708e-07f)));
        *conic_3 = exp_0(- scale_6);
        float3  _S771 = mean_4 - - mul_0(transpose_0(R_7), t_6);
        float3  _S772 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
        *rgb_3 = _S772;
        float _S773 = _S771.x;
        float _S774 = _S771.y;
        float _S775 = _S771.z;
        float norm_3 = (F32_sqrt((_S773 * _S773 + _S774 * _S774 + _S775 * _S775)));
        float x_23 = _S773 / norm_3;
        float y_6 = _S774 / norm_3;
        float z_3 = _S775 / norm_3;
        float3  _S776 = _S772 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_3)[int(1)] + make_float3 (z_3) * (*sh_coeffs_3)[int(2)] - make_float3 (x_23) * (*sh_coeffs_3)[int(3)]);
        *rgb_3 = _S776;
        float z2_11 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_23 * x_23 - y_6 * y_6;
        float fS1_3 = 2.0f * x_23 * y_6;
        float3  _S777 = _S776 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_6) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_11 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_23) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
        *rgb_3 = _S777;
        float fTmp0C_3 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        *rgb_3 = max_0(_S777 + (make_float3 (-0.59004360437393188f * (x_23 * fS1_3 + y_6 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_6) * (*sh_coeffs_3)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_11 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_23) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_23 * fC1_3 - y_6 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_fisheye(bool antialiased_4, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_8, float3  t_7, float fx_11, float fy_11, float cx_8, float cy_8, FixedArray<float, 10>  * dist_coeffs_20, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_12, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    for(;;)
    {
        float3  mean_c_4 = mul_0(R_8, mean_5) + t_7;
        float _S778 = length_1(mean_c_4);
        bool _S779;
        if(_S778 < near_plane_4)
        {
            _S779 = true;
        }
        else
        {
            _S779 = _S778 > far_plane_4;
        }
        if(_S779)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_8 = make_float4 (fx_11, fy_11, cx_8, cy_8);
        Matrix<float, 2, 2>  covar2d_4;
        float3  _S780 = exp_0(scale_7);
        float x_24 = quat_8.y;
        float x2_8 = x_24 * x_24;
        float y2_8 = quat_8.z * quat_8.z;
        float z2_12 = quat_8.w * quat_8.w;
        float xy_8 = quat_8.y * quat_8.z;
        float xz_8 = quat_8.y * quat_8.w;
        float yz_8 = quat_8.z * quat_8.w;
        float wx_8 = quat_8.x * quat_8.y;
        float wy_8 = quat_8.x * quat_8.z;
        float wz_8 = quat_8.x * quat_8.w;
        Matrix<float, 3, 3>  _S781 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_12), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_12), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))));
        SigmaPoints_0 ret_1;
        (&ret_1)->p_0[int(0)] = mean_5;
        (&ret_1)->w_mean_0[int(0)] = 0.0f;
        (&ret_1)->w_cov_0[int(0)] = 2.0f;
        float _S782 = (F32_sqrt((3.0f)));
        float3  delta_3 = make_float3 (_S782 * _S780.x) * _S781.rows[0U];
        float3  _S783 = mean_5 + delta_3;
        float3  _S784 = mean_5 - delta_3;
        float3  delta_4 = make_float3 (_S782 * _S780.y) * _S781.rows[1U];
        float3  _S785 = mean_5 + delta_4;
        float3  _S786 = mean_5 - delta_4;
        float3  delta_5 = make_float3 (_S782 * _S780.z) * _S781.rows[2U];
        float3  _S787 = mean_5 + delta_5;
        float3  _S788 = mean_5 - delta_5;
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
        (&ret_1)->p_0[1U] = mul_0(R_8, _S783) + t_7;
        (&ret_1)->p_0[2U] = mul_0(R_8, _S785) + t_7;
        (&ret_1)->p_0[3U] = mul_0(R_8, _S787) + t_7;
        (&ret_1)->p_0[4U] = mul_0(R_8, _S784) + t_7;
        (&ret_1)->p_0[5U] = mul_0(R_8, _S786) + t_7;
        (&ret_1)->p_0[6U] = mul_0(R_8, _S788) + t_7;
        SigmaPoints_0 _S789 = ret_1;
        bool _S790 = fisheye_proj_3dgs_ut_0(&_S789, intrins_8, dist_coeffs_20, &covar2d_4, mean2d_12);
        if(!(true & _S790))
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S791 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S791;
        float _S792 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S792;
        float det_blur_5 = _S791 * _S792 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
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
            _S779 = true;
        }
        else
        {
            _S779 = xmin_4 >= float(image_width_4);
        }
        if(_S779)
        {
            _S779 = true;
        }
        else
        {
            _S779 = ymax_4 <= 0.0f;
        }
        if(_S779)
        {
            _S779 = true;
        }
        else
        {
            _S779 = ymin_4 >= float(image_height_4);
        }
        if(_S779)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_4 = make_int4 (int(xmin_4), int(ymin_4), int(xmax_4), int(ymax_4));
        *depth_4 = 0.5f * (F32_log((dot_0(mean_c_4, mean_c_4) + 9.99999997475242708e-07f)));
        *conic_4 = exp_0(- scale_7);
        float3  _S793 = mean_5 - - mul_0(transpose_0(R_8), t_7);
        float3  _S794 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
        *rgb_4 = _S794;
        float _S795 = _S793.x;
        float _S796 = _S793.y;
        float _S797 = _S793.z;
        float norm_4 = (F32_sqrt((_S795 * _S795 + _S796 * _S796 + _S797 * _S797)));
        float x_25 = _S795 / norm_4;
        float y_7 = _S796 / norm_4;
        float z_4 = _S797 / norm_4;
        float3  _S798 = _S794 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_4)[int(1)] + make_float3 (z_4) * (*sh_coeffs_4)[int(2)] - make_float3 (x_25) * (*sh_coeffs_4)[int(3)]);
        *rgb_4 = _S798;
        float z2_13 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_25 * x_25 - y_7 * y_7;
        float fS1_4 = 2.0f * x_25 * y_7;
        float3  _S799 = _S798 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_7) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_13 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_25) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
        *rgb_4 = _S799;
        float fTmp0C_4 = -2.28522896766662598f * z2_13 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        *rgb_4 = max_0(_S799 + (make_float3 (-0.59004360437393188f * (x_25 * fS1_4 + y_7 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_7) * (*sh_coeffs_4)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_13 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_25) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_25 * fC1_4 - y_7 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_5, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_9, float3  t_8, float fx_12, float fy_12, float cx_9, float cy_9, FixedArray<float, 10>  * dist_coeffs_21, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_13, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_9, mean_6) + t_8;
    float3  _S800 = exp_0(scale_8);
    float x_26 = quat_9.y;
    float x2_9 = x_26 * x_26;
    float y2_9 = quat_9.z * quat_9.z;
    float z2_14 = quat_9.w * quat_9.w;
    float xy_9 = quat_9.y * quat_9.z;
    float xz_9 = quat_9.y * quat_9.w;
    float yz_9 = quat_9.z * quat_9.w;
    float wx_9 = quat_9.x * quat_9.y;
    float wy_9 = quat_9.x * quat_9.z;
    float wz_9 = quat_9.x * quat_9.w;
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_14), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_14), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S800.x, 0.0f, 0.0f, 0.0f, _S800.y, 0.0f, 0.0f, 0.0f, _S800.z));
    Matrix<float, 3, 3>  _S801 = transpose_0(R_9);
    float _S802 = float(image_width_5);
    float _S803 = float(image_height_5);
    float _S804 = 0.30000001192092896f * (0.5f * _S802 / fx_12);
    float _S805 = 0.30000001192092896f * (0.5f * _S803 / fy_12);
    float rz_2 = 1.0f / mean_c_5.z;
    float rz2_2 = rz_2 * rz_2;
    Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_12 * rz_2, 0.0f, - fx_12 * (mean_c_5.z * (F32_min(((_S802 - cx_9) / fx_12 + _S804), ((F32_max((- (cx_9 / fx_12 + _S804)), (mean_c_5.x * rz_2))))))) * rz2_2, 0.0f, fy_12 * rz_2, - fy_12 * (mean_c_5.z * (F32_min(((_S803 - cy_9) / fy_12 + _S805), ((F32_max((- (cy_9 / fy_12 + _S805)), (mean_c_5.y * rz_2))))))) * rz2_2);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_6, mul_4(mul_4(R_9, mul_4(M_5, transpose_0(M_5))), _S801)), transpose_1(J_6));
    *mean2d_13 = make_float2 (fx_12 * mean_c_5.x * rz_2 + cx_9, fy_12 * mean_c_5.y * rz_2 + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S806 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S806;
    float _S807 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S807;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S806 * _S807 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S808 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_4, - covar2d_5.rows[int(0)].y * invdet_4, - covar2d_5.rows[int(1)].x * invdet_4, covar2d_5.rows[int(0)].x * invdet_4);
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
    *conic_5 = make_float3 (_S808.rows[int(0)].x, _S808.rows[int(0)].y, _S808.rows[int(1)].y);
    float3  _S809 = mean_6 - - mul_0(_S801, t_8);
    float3  _S810 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S810;
    float _S811 = _S809.x;
    float _S812 = _S809.y;
    float _S813 = _S809.z;
    float norm_5 = (F32_sqrt((_S811 * _S811 + _S812 * _S812 + _S813 * _S813)));
    float x_27 = _S811 / norm_5;
    float y_8 = _S812 / norm_5;
    float z_5 = _S813 / norm_5;
    float3  _S814 = _S810 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_5)[int(1)] + make_float3 (z_5) * (*sh_coeffs_5)[int(2)] - make_float3 (x_27) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S814;
    float z2_15 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_27 * x_27 - y_8 * y_8;
    float fS1_5 = 2.0f * x_27 * y_8;
    float3  _S815 = _S814 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_8) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_15 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_27) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S815;
    float fTmp0C_5 = -2.28522896766662598f * z2_15 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    *rgb_5 = max_0(_S815 + (make_float3 (-0.59004360437393188f * (x_27 * fS1_5 + y_8 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_8) * (*sh_coeffs_5)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_15 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_27) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_27 * fC1_5 - y_8 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_6, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_10, float3  t_9, float fx_13, float fy_13, float cx_10, float cy_10, FixedArray<float, 10>  * dist_coeffs_22, uint image_width_6, uint image_height_6, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * mean2d_14, float * depth_6, float3  * conic_6, float * opacity_6, float3  * rgb_6)
{
    float3  mean_c_6 = mul_0(R_10, mean_7) + t_9;
    float3  _S816 = exp_0(scale_9);
    float x_28 = quat_10.y;
    float x2_10 = x_28 * x_28;
    float y2_10 = quat_10.z * quat_10.z;
    float z2_16 = quat_10.w * quat_10.w;
    float xy_10 = quat_10.y * quat_10.z;
    float xz_10 = quat_10.y * quat_10.w;
    float yz_10 = quat_10.z * quat_10.w;
    float wx_10 = quat_10.x * quat_10.y;
    float wy_10 = quat_10.x * quat_10.z;
    float wz_10 = quat_10.x * quat_10.w;
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_16), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_16), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S816.x, 0.0f, 0.0f, 0.0f, _S816.y, 0.0f, 0.0f, 0.0f, _S816.z));
    Matrix<float, 3, 3>  _S817 = transpose_0(R_10);
    Matrix<float, 2, 2>  covar2d_6;
    bool _S818 = fisheye_proj_3dgs_1(mean_c_6, mul_4(mul_4(R_10, mul_4(M_6, transpose_0(M_6))), _S817), make_float4 (fx_13, fy_13, cx_10, cy_10), dist_coeffs_22, &covar2d_6, mean2d_14);
    float det_orig_7 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
    float _S819 = *&(((&covar2d_6)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(0)))->x) = _S819;
    float _S820 = *&(((&covar2d_6)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(1)))->y) = _S820;
    float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / (_S819 * _S820 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
    Matrix<float, 2, 2>  _S821 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_5, - covar2d_6.rows[int(0)].y * invdet_5, - covar2d_6.rows[int(1)].x * invdet_5, covar2d_6.rows[int(0)].x * invdet_5);
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
    *conic_6 = make_float3 (_S821.rows[int(0)].x, _S821.rows[int(0)].y, _S821.rows[int(1)].y);
    float3  _S822 = mean_7 - - mul_0(_S817, t_9);
    float3  _S823 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)];
    *rgb_6 = _S823;
    float _S824 = _S822.x;
    float _S825 = _S822.y;
    float _S826 = _S822.z;
    float norm_6 = (F32_sqrt((_S824 * _S824 + _S825 * _S825 + _S826 * _S826)));
    float x_29 = _S824 / norm_6;
    float y_9 = _S825 / norm_6;
    float z_6 = _S826 / norm_6;
    float3  _S827 = _S823 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_6)[int(1)] + make_float3 (z_6) * (*sh_coeffs_6)[int(2)] - make_float3 (x_29) * (*sh_coeffs_6)[int(3)]);
    *rgb_6 = _S827;
    float z2_17 = z_6 * z_6;
    float fTmp0B_6 = -1.09254848957061768f * z_6;
    float fC1_6 = x_29 * x_29 - y_9 * y_9;
    float fS1_6 = 2.0f * x_29 * y_9;
    float3  _S828 = _S827 + (make_float3 (0.54627424478530884f * fS1_6) * (*sh_coeffs_6)[int(4)] + make_float3 (fTmp0B_6 * y_9) * (*sh_coeffs_6)[int(5)] + make_float3 (0.94617468118667603f * z2_17 - 0.31539157032966614f) * (*sh_coeffs_6)[int(6)] + make_float3 (fTmp0B_6 * x_29) * (*sh_coeffs_6)[int(7)] + make_float3 (0.54627424478530884f * fC1_6) * (*sh_coeffs_6)[int(8)]);
    *rgb_6 = _S828;
    float fTmp0C_6 = -2.28522896766662598f * z2_17 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_6;
    *rgb_6 = max_0(_S828 + (make_float3 (-0.59004360437393188f * (x_29 * fS1_6 + y_9 * fC1_6)) * (*sh_coeffs_6)[int(9)] + make_float3 (fTmp1B_6 * fS1_6) * (*sh_coeffs_6)[int(10)] + make_float3 (fTmp0C_6 * y_9) * (*sh_coeffs_6)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_17 - 1.11952900886535645f)) * (*sh_coeffs_6)[int(12)] + make_float3 (fTmp0C_6 * x_29) * (*sh_coeffs_6)[int(13)] + make_float3 (fTmp1B_6 * fC1_6) * (*sh_coeffs_6)[int(14)] + make_float3 (-0.59004360437393188f * (x_29 * fC1_6 - y_9 * fS1_6)) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_7, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_11, float3  t_10, float fx_14, float fy_14, float cx_11, float cy_11, FixedArray<float, 10>  * dist_coeffs_23, uint image_width_7, uint image_height_7, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * mean2d_15, float * depth_7, float3  * conic_7, float * opacity_7, float3  * rgb_7)
{
    float3  mean_c_7 = mul_0(R_11, mean_8) + t_10;
    float3  _S829 = exp_0(scale_10);
    float x_30 = quat_11.y;
    float x2_11 = x_30 * x_30;
    float y2_11 = quat_11.z * quat_11.z;
    float z2_18 = quat_11.w * quat_11.w;
    float xy_11 = quat_11.y * quat_11.z;
    float xz_11 = quat_11.y * quat_11.w;
    float yz_11 = quat_11.z * quat_11.w;
    float wx_11 = quat_11.x * quat_11.y;
    float wy_11 = quat_11.x * quat_11.z;
    float wz_11 = quat_11.x * quat_11.w;
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_18), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_18), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))), makeMatrix<float, 3, 3> (_S829.x, 0.0f, 0.0f, 0.0f, _S829.y, 0.0f, 0.0f, 0.0f, _S829.z));
    Matrix<float, 3, 3>  _S830 = transpose_0(R_11);
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_14, 0.0f, 0.0f, 0.0f, fy_14, 0.0f);
    Matrix<float, 2, 2>  covar2d_7 = mul_6(mul_5(J_7, mul_4(mul_4(R_11, mul_4(M_7, transpose_0(M_7))), _S830)), transpose_1(J_7));
    *mean2d_15 = make_float2 (fx_14 * mean_c_7.x + cx_11, fy_14 * mean_c_7.y + cy_11);
    float det_orig_8 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
    float _S831 = *&(((&covar2d_7)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(0)))->x) = _S831;
    float _S832 = *&(((&covar2d_7)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(1)))->y) = _S832;
    float compensation_8 = (F32_sqrt(((F32_max((0.0f), (det_orig_8 / (_S831 * _S832 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
    Matrix<float, 2, 2>  _S833 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_6, - covar2d_7.rows[int(0)].y * invdet_6, - covar2d_7.rows[int(1)].x * invdet_6, covar2d_7.rows[int(0)].x * invdet_6);
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
    *conic_7 = make_float3 (_S833.rows[int(0)].x, _S833.rows[int(0)].y, _S833.rows[int(1)].y);
    float3  _S834 = mean_8 - - mul_0(_S830, t_10);
    float3  _S835 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)];
    *rgb_7 = _S835;
    float _S836 = _S834.x;
    float _S837 = _S834.y;
    float _S838 = _S834.z;
    float norm_7 = (F32_sqrt((_S836 * _S836 + _S837 * _S837 + _S838 * _S838)));
    float x_31 = _S836 / norm_7;
    float y_10 = _S837 / norm_7;
    float z_7 = _S838 / norm_7;
    float3  _S839 = _S835 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * (*sh_coeffs_7)[int(1)] + make_float3 (z_7) * (*sh_coeffs_7)[int(2)] - make_float3 (x_31) * (*sh_coeffs_7)[int(3)]);
    *rgb_7 = _S839;
    float z2_19 = z_7 * z_7;
    float fTmp0B_7 = -1.09254848957061768f * z_7;
    float fC1_7 = x_31 * x_31 - y_10 * y_10;
    float fS1_7 = 2.0f * x_31 * y_10;
    float3  _S840 = _S839 + (make_float3 (0.54627424478530884f * fS1_7) * (*sh_coeffs_7)[int(4)] + make_float3 (fTmp0B_7 * y_10) * (*sh_coeffs_7)[int(5)] + make_float3 (0.94617468118667603f * z2_19 - 0.31539157032966614f) * (*sh_coeffs_7)[int(6)] + make_float3 (fTmp0B_7 * x_31) * (*sh_coeffs_7)[int(7)] + make_float3 (0.54627424478530884f * fC1_7) * (*sh_coeffs_7)[int(8)]);
    *rgb_7 = _S840;
    float fTmp0C_7 = -2.28522896766662598f * z2_19 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_7;
    *rgb_7 = max_0(_S840 + (make_float3 (-0.59004360437393188f * (x_31 * fS1_7 + y_10 * fC1_7)) * (*sh_coeffs_7)[int(9)] + make_float3 (fTmp1B_7 * fS1_7) * (*sh_coeffs_7)[int(10)] + make_float3 (fTmp0C_7 * y_10) * (*sh_coeffs_7)[int(11)] + make_float3 (z_7 * (1.86588168144226074f * z2_19 - 1.11952900886535645f)) * (*sh_coeffs_7)[int(12)] + make_float3 (fTmp0C_7 * x_31) * (*sh_coeffs_7)[int(13)] + make_float3 (fTmp1B_7 * fC1_7) * (*sh_coeffs_7)[int(14)] + make_float3 (-0.59004360437393188f * (x_31 * fC1_7 - y_10 * fS1_7)) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_persp_differentiable(bool antialiased_8, float3  mean_9, float4  quat_12, float3  scale_11, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_12, float3  t_11, float fx_15, float fy_15, float cx_12, float cy_12, FixedArray<float, 10>  * dist_coeffs_24, uint image_width_8, uint image_height_8, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * mean2d_16, float * depth_8, float3  * conic_8, float * opacity_8, float3  * rgb_8)
{
    float3  mean_c_8 = mul_0(R_12, mean_9) + t_11;
    float4  intrins_9 = make_float4 (fx_15, fy_15, cx_12, cy_12);
    float3  _S841 = exp_0(scale_11);
    float x_32 = quat_12.y;
    float x2_12 = x_32 * x_32;
    float y2_12 = quat_12.z * quat_12.z;
    float z2_20 = quat_12.w * quat_12.w;
    float xy_12 = quat_12.y * quat_12.z;
    float xz_12 = quat_12.y * quat_12.w;
    float yz_12 = quat_12.z * quat_12.w;
    float wx_12 = quat_12.x * quat_12.y;
    float wy_12 = quat_12.x * quat_12.z;
    float wz_12 = quat_12.x * quat_12.w;
    Matrix<float, 3, 3>  _S842 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_20), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_20), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))));
    SigmaPoints_0 ret_2;
    (&ret_2)->p_0[int(0)] = mean_9;
    (&ret_2)->w_mean_0[int(0)] = 0.0f;
    (&ret_2)->w_cov_0[int(0)] = 2.0f;
    float _S843 = (F32_sqrt((3.0f)));
    float3  delta_6 = make_float3 (_S843 * _S841.x) * _S842.rows[0U];
    float3  _S844 = mean_9 + delta_6;
    float3  _S845 = mean_9 - delta_6;
    float3  delta_7 = make_float3 (_S843 * _S841.y) * _S842.rows[1U];
    float3  _S846 = mean_9 + delta_7;
    float3  _S847 = mean_9 - delta_7;
    float3  delta_8 = make_float3 (_S843 * _S841.z) * _S842.rows[2U];
    float3  _S848 = mean_9 + delta_8;
    float3  _S849 = mean_9 - delta_8;
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
    (&ret_2)->p_0[1U] = mul_0(R_12, _S844) + t_11;
    (&ret_2)->p_0[2U] = mul_0(R_12, _S846) + t_11;
    (&ret_2)->p_0[3U] = mul_0(R_12, _S848) + t_11;
    (&ret_2)->p_0[4U] = mul_0(R_12, _S845) + t_11;
    (&ret_2)->p_0[5U] = mul_0(R_12, _S847) + t_11;
    (&ret_2)->p_0[6U] = mul_0(R_12, _S849) + t_11;
    SigmaPoints_0 _S850 = ret_2;
    Matrix<float, 2, 2>  covar2d_8;
    bool _S851 = persp_proj_3dgs_ut_1(&_S850, intrins_9, dist_coeffs_24, image_width_8, image_height_8, &covar2d_8, mean2d_16);
    float det_orig_9 = *&(((&covar2d_8)->rows + (int(0)))->x) * *&(((&covar2d_8)->rows + (int(1)))->y) - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x);
    float _S852 = *&(((&covar2d_8)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(0)))->x) = _S852;
    float _S853 = *&(((&covar2d_8)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(1)))->y) = _S853;
    float compensation_9 = (F32_sqrt(((F32_max((0.0f), (det_orig_9 / (_S852 * _S853 - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x))))))));
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
    *conic_8 = exp_0(- scale_11);
    float3  _S854 = mean_9 - - mul_0(transpose_0(R_12), t_11);
    float3  _S855 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)];
    *rgb_8 = _S855;
    float _S856 = _S854.x;
    float _S857 = _S854.y;
    float _S858 = _S854.z;
    float norm_8 = (F32_sqrt((_S856 * _S856 + _S857 * _S857 + _S858 * _S858)));
    float x_33 = _S856 / norm_8;
    float y_11 = _S857 / norm_8;
    float z_8 = _S858 / norm_8;
    float3  _S859 = _S855 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_8)[int(1)] + make_float3 (z_8) * (*sh_coeffs_8)[int(2)] - make_float3 (x_33) * (*sh_coeffs_8)[int(3)]);
    *rgb_8 = _S859;
    float z2_21 = z_8 * z_8;
    float fTmp0B_8 = -1.09254848957061768f * z_8;
    float fC1_8 = x_33 * x_33 - y_11 * y_11;
    float fS1_8 = 2.0f * x_33 * y_11;
    float3  _S860 = _S859 + (make_float3 (0.54627424478530884f * fS1_8) * (*sh_coeffs_8)[int(4)] + make_float3 (fTmp0B_8 * y_11) * (*sh_coeffs_8)[int(5)] + make_float3 (0.94617468118667603f * z2_21 - 0.31539157032966614f) * (*sh_coeffs_8)[int(6)] + make_float3 (fTmp0B_8 * x_33) * (*sh_coeffs_8)[int(7)] + make_float3 (0.54627424478530884f * fC1_8) * (*sh_coeffs_8)[int(8)]);
    *rgb_8 = _S860;
    float fTmp0C_8 = -2.28522896766662598f * z2_21 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_8;
    *rgb_8 = max_0(_S860 + (make_float3 (-0.59004360437393188f * (x_33 * fS1_8 + y_11 * fC1_8)) * (*sh_coeffs_8)[int(9)] + make_float3 (fTmp1B_8 * fS1_8) * (*sh_coeffs_8)[int(10)] + make_float3 (fTmp0C_8 * y_11) * (*sh_coeffs_8)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_21 - 1.11952900886535645f)) * (*sh_coeffs_8)[int(12)] + make_float3 (fTmp0C_8 * x_33) * (*sh_coeffs_8)[int(13)] + make_float3 (fTmp1B_8 * fC1_8) * (*sh_coeffs_8)[int(14)] + make_float3 (-0.59004360437393188f * (x_33 * fC1_8 - y_11 * fS1_8)) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_fisheye_differentiable(bool antialiased_9, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_9, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_13, float3  t_12, float fx_16, float fy_16, float cx_13, float cy_13, FixedArray<float, 10>  * dist_coeffs_25, uint image_width_9, uint image_height_9, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * mean2d_17, float * depth_9, float3  * conic_9, float * opacity_9, float3  * rgb_9)
{
    float3  mean_c_9 = mul_0(R_13, mean_10) + t_12;
    float4  intrins_10 = make_float4 (fx_16, fy_16, cx_13, cy_13);
    float3  _S861 = exp_0(scale_12);
    float x_34 = quat_13.y;
    float x2_13 = x_34 * x_34;
    float y2_13 = quat_13.z * quat_13.z;
    float z2_22 = quat_13.w * quat_13.w;
    float xy_13 = quat_13.y * quat_13.z;
    float xz_13 = quat_13.y * quat_13.w;
    float yz_13 = quat_13.z * quat_13.w;
    float wx_13 = quat_13.x * quat_13.y;
    float wy_13 = quat_13.x * quat_13.z;
    float wz_13 = quat_13.x * quat_13.w;
    Matrix<float, 3, 3>  _S862 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_22), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_22), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13))));
    SigmaPoints_0 ret_3;
    (&ret_3)->p_0[int(0)] = mean_10;
    (&ret_3)->w_mean_0[int(0)] = 0.0f;
    (&ret_3)->w_cov_0[int(0)] = 2.0f;
    float _S863 = (F32_sqrt((3.0f)));
    float3  delta_9 = make_float3 (_S863 * _S861.x) * _S862.rows[0U];
    float3  _S864 = mean_10 + delta_9;
    float3  _S865 = mean_10 - delta_9;
    float3  delta_10 = make_float3 (_S863 * _S861.y) * _S862.rows[1U];
    float3  _S866 = mean_10 + delta_10;
    float3  _S867 = mean_10 - delta_10;
    float3  delta_11 = make_float3 (_S863 * _S861.z) * _S862.rows[2U];
    float3  _S868 = mean_10 + delta_11;
    float3  _S869 = mean_10 - delta_11;
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
    (&ret_3)->p_0[1U] = mul_0(R_13, _S864) + t_12;
    (&ret_3)->p_0[2U] = mul_0(R_13, _S866) + t_12;
    (&ret_3)->p_0[3U] = mul_0(R_13, _S868) + t_12;
    (&ret_3)->p_0[4U] = mul_0(R_13, _S865) + t_12;
    (&ret_3)->p_0[5U] = mul_0(R_13, _S867) + t_12;
    (&ret_3)->p_0[6U] = mul_0(R_13, _S869) + t_12;
    SigmaPoints_0 _S870 = ret_3;
    Matrix<float, 2, 2>  covar2d_9;
    bool _S871 = fisheye_proj_3dgs_ut_1(&_S870, intrins_10, dist_coeffs_25, &covar2d_9, mean2d_17);
    float det_orig_10 = *&(((&covar2d_9)->rows + (int(0)))->x) * *&(((&covar2d_9)->rows + (int(1)))->y) - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x);
    float _S872 = *&(((&covar2d_9)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(0)))->x) = _S872;
    float _S873 = *&(((&covar2d_9)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(1)))->y) = _S873;
    float compensation_10 = (F32_sqrt(((F32_max((0.0f), (det_orig_10 / (_S872 * _S873 - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x))))))));
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
    *conic_9 = exp_0(- scale_12);
    float3  _S874 = mean_10 - - mul_0(transpose_0(R_13), t_12);
    float3  _S875 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
    *rgb_9 = _S875;
    float _S876 = _S874.x;
    float _S877 = _S874.y;
    float _S878 = _S874.z;
    float norm_9 = (F32_sqrt((_S876 * _S876 + _S877 * _S877 + _S878 * _S878)));
    float x_35 = _S876 / norm_9;
    float y_12 = _S877 / norm_9;
    float z_9 = _S878 / norm_9;
    float3  _S879 = _S875 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * (*sh_coeffs_9)[int(1)] + make_float3 (z_9) * (*sh_coeffs_9)[int(2)] - make_float3 (x_35) * (*sh_coeffs_9)[int(3)]);
    *rgb_9 = _S879;
    float z2_23 = z_9 * z_9;
    float fTmp0B_9 = -1.09254848957061768f * z_9;
    float fC1_9 = x_35 * x_35 - y_12 * y_12;
    float fS1_9 = 2.0f * x_35 * y_12;
    float3  _S880 = _S879 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_12) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_23 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_35) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
    *rgb_9 = _S880;
    float fTmp0C_9 = -2.28522896766662598f * z2_23 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_9;
    *rgb_9 = max_0(_S880 + (make_float3 (-0.59004360437393188f * (x_35 * fS1_9 + y_12 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_12) * (*sh_coeffs_9)[int(11)] + make_float3 (z_9 * (1.86588168144226074f * z2_23 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_35) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_35 * fC1_9 - y_12 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S881, float3  _S882)
{
    return mul_0(_S881, _S882);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S883)
{
    return exp_0(_S883);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_2(Matrix<float, 3, 3>  _S884, Matrix<float, 3, 3>  _S885)
{
    return mul_4(_S884, _S885);
}

inline __device__ float s_primal_ctx_max_0(float _S886, float _S887)
{
    return (F32_max((_S886), (_S887)));
}

inline __device__ float s_primal_ctx_min_0(float _S888, float _S889)
{
    return (F32_min((_S888), (_S889)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S890, Matrix<float, 3, 3>  _S891)
{
    return mul_5(_S890, _S891);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_4(Matrix<float, 2, 3>  _S892, Matrix<float, 3, 2>  _S893)
{
    return mul_6(_S892, _S893);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S894)
{
    return (F32_sqrt((_S894)));
}

inline __device__ float s_primal_ctx_exp_1(float _S895)
{
    return (F32_exp((_S895)));
}

inline __device__ float s_primal_ctx_log_0(float _S896)
{
    return (F32_log((_S896)));
}

inline __device__ float s_primal_ctx_dot_0(float3  _S897, float3  _S898)
{
    return dot_0(_S897, _S898);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S899, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S900, float3  _S901)
{
    _d_max_vector_0(_S899, _S900, _S901);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S902, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S903, float3  _S904)
{
    _d_mul_0(_S902, _S903, _S904);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S905, float _S906)
{
    _d_log_0(_S905, _S906);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S907, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S908, float _S909)
{
    _d_dot_0(_S907, _S908, _S909);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S910, DiffPair_float_0 * _S911, float _S912)
{
    _d_min_0(_S910, _S911, _S912);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S913, float _S914)
{
    _d_exp_0(_S913, _S914);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S915, DiffPair_float_0 * _S916, float _S917)
{
    _d_max_0(_S915, _S916, _S917);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S918, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S919, Matrix<float, 2, 2>  _S920)
{
    mul_3(_S918, _S919, _S920);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S921, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S922, Matrix<float, 2, 3>  _S923)
{
    mul_2(_S921, _S922, _S923);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S924, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S925, Matrix<float, 3, 3>  _S926)
{
    mul_1(_S924, _S925, _S926);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S927, float3  _S928)
{
    _d_exp_vector_0(_S927, _S928);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_14, float3  t_13, float fx_17, float fy_17, float cx_14, float cy_14, FixedArray<float, 10>  * dist_coeffs_26, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_t_1)
{
    float3  mean_c_10 = s_primal_ctx_mul_1(R_14, mean_11) + t_13;
    float3  _S929 = s_primal_ctx_exp_0(scale_13);
    float _S930 = quat_14.y;
    float x2_14 = _S930 * _S930;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_24 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S931 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_24), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_24), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S929.x, 0.0f, 0.0f, 0.0f, _S929.y, 0.0f, 0.0f, 0.0f, _S929.z);
    Matrix<float, 3, 3>  _S932 = s_primal_ctx_mul_2(_S931, S_0);
    Matrix<float, 3, 3>  _S933 = transpose_0(_S932);
    Matrix<float, 3, 3>  _S934 = s_primal_ctx_mul_2(_S932, _S933);
    Matrix<float, 3, 3>  _S935 = s_primal_ctx_mul_2(R_14, _S934);
    Matrix<float, 3, 3>  _S936 = transpose_0(R_14);
    Matrix<float, 3, 3>  _S937 = s_primal_ctx_mul_2(_S935, _S936);
    float _S938 = float(image_width_10);
    float _S939 = float(image_height_10);
    float _S940 = 0.30000001192092896f * (0.5f * _S938 / fx_17);
    float lim_x_pos_2 = (_S938 - cx_14) / fx_17 + _S940;
    float _S941 = 0.30000001192092896f * (0.5f * _S939 / fy_17);
    float lim_y_pos_2 = (_S939 - cy_14) / fy_17 + _S941;
    float rz_3 = 1.0f / mean_c_10.z;
    float _S942 = mean_c_10.z * mean_c_10.z;
    float rz2_3 = rz_3 * rz_3;
    float _S943 = - (cx_14 / fx_17 + _S940);
    float _S944 = mean_c_10.x * rz_3;
    float _S945 = s_primal_ctx_max_0(_S943, _S944);
    float _S946 = s_primal_ctx_min_0(lim_x_pos_2, _S945);
    float _S947 = - (cy_14 / fy_17 + _S941);
    float _S948 = mean_c_10.y * rz_3;
    float _S949 = s_primal_ctx_max_0(_S947, _S948);
    float _S950 = s_primal_ctx_min_0(lim_y_pos_2, _S949);
    float _S951 = - fx_17;
    float _S952 = _S951 * (mean_c_10.z * _S946);
    float _S953 = - fy_17;
    float _S954 = _S953 * (mean_c_10.z * _S950);
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_17 * rz_3, 0.0f, _S952 * rz2_3, 0.0f, fy_17 * rz_3, _S954 * rz2_3);
    Matrix<float, 2, 3>  _S955 = s_primal_ctx_mul_3(J_8, _S937);
    Matrix<float, 3, 2>  _S956 = transpose_1(J_8);
    Matrix<float, 2, 2>  _S957 = s_primal_ctx_mul_4(_S955, _S956);
    float _S958 = fx_17 * mean_c_10.x;
    float _S959 = fy_17 * mean_c_10.y;
    float _S960 = _S957.rows[int(0)].y * _S957.rows[int(1)].x;
    float det_orig_11 = _S957.rows[int(0)].x * _S957.rows[int(1)].y - _S960;
    float _S961 = _S957.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S962 = _S957;
    *&(((&_S962)->rows + (int(0)))->x) = _S961;
    float _S963 = _S957.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S962)->rows + (int(1)))->y) = _S963;
    Matrix<float, 2, 2>  _S964 = _S962;
    Matrix<float, 2, 2>  _S965 = _S962;
    float det_blur_6 = _S961 * _S963 - _S960;
    float _S966 = det_orig_11 / det_blur_6;
    float _S967 = det_blur_6 * det_blur_6;
    float _S968 = s_primal_ctx_max_0(0.0f, _S966);
    float _S969 = s_primal_ctx_sqrt_0(_S968);
    float invdet_7 = 1.0f / det_blur_6;
    float _S970 = - _S957.rows[int(0)].y;
    float _S971 = - _S957.rows[int(1)].x;
    float _S972 = - in_opacity_10;
    float _S973 = 1.0f + s_primal_ctx_exp_1(_S972);
    float _S974 = 1.0f / _S973;
    float _S975 = _S973 * _S973;
    float _S976;
    if(antialiased_10)
    {
        _S976 = _S974 * _S969;
    }
    else
    {
        _S976 = _S974;
    }
    float _S977 = _S976 / 0.00392156885936856f;
    float _S978 = 2.0f * s_primal_ctx_log_0(_S977);
    float _S979 = s_primal_ctx_sqrt_0(_S978);
    float _S980 = _S964.rows[int(0)].x;
    float _S981 = _S965.rows[int(1)].y;
    float _S982 = s_primal_ctx_dot_0(mean_c_10, mean_c_10) + 9.99999997475242708e-07f;
    float3  _S983 = mean_11 - - s_primal_ctx_mul_1(_S936, t_13);
    float _S984 = _S983.x;
    float _S985 = _S983.y;
    float _S986 = _S983.z;
    float _S987 = _S984 * _S984 + _S985 * _S985 + _S986 * _S986;
    float _S988 = s_primal_ctx_sqrt_0(_S987);
    float x_36 = _S984 / _S988;
    float3  _S989 = make_float3 (x_36);
    float _S990 = _S988 * _S988;
    float y_13 = _S985 / _S988;
    float z_10 = _S986 / _S988;
    float3  _S991 = make_float3 (z_10);
    float _S992 = - y_13;
    float3  _S993 = make_float3 (_S992);
    float z2_25 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_36 * x_36 - y_13 * y_13;
    float _S994 = 2.0f * x_36;
    float fS1_10 = _S994 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_25 - 0.31539157032966614f;
    float3  _S995 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_36;
    float3  _S996 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S997 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S998 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S999 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_25 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S1000 = 1.86588168144226074f * z2_25 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S1000;
    float3  _S1001 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_36;
    float3  _S1002 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S1003 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S1004 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S1005 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_36 * fC1_10 - y_13 * fS1_10);
    float3  _S1006 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_36 * fS1_10 + y_13 * fC1_10);
    float3  _S1007 = make_float3 (pSH9_0);
    float3  _S1008 = make_float3 (0.0f);
    float3  _S1009 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1010;
    (&_S1010)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S992) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_36) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S1010)->differential_0 = _S1009;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1011;
    (&_S1011)->primal_0 = _S1008;
    (&_S1011)->differential_0 = _S1009;
    s_bwd_prop_max_0(&_S1010, &_S1011, v_rgb_0);
    float3  _S1012 = _S1006 * _S1010.differential_0;
    float3  _S1013 = (*sh_coeffs_10)[int(15)] * _S1010.differential_0;
    float3  _S1014 = _S1004 * _S1010.differential_0;
    float3  _S1015 = (*sh_coeffs_10)[int(14)] * _S1010.differential_0;
    float3  _S1016 = _S1002 * _S1010.differential_0;
    float3  _S1017 = (*sh_coeffs_10)[int(13)] * _S1010.differential_0;
    float3  _S1018 = _S1001 * _S1010.differential_0;
    float3  _S1019 = (*sh_coeffs_10)[int(12)] * _S1010.differential_0;
    float3  _S1020 = _S1003 * _S1010.differential_0;
    float3  _S1021 = (*sh_coeffs_10)[int(11)] * _S1010.differential_0;
    float3  _S1022 = _S1005 * _S1010.differential_0;
    float3  _S1023 = (*sh_coeffs_10)[int(10)] * _S1010.differential_0;
    float3  _S1024 = _S1007 * _S1010.differential_0;
    float3  _S1025 = (*sh_coeffs_10)[int(9)] * _S1010.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S1025.x + _S1025.y + _S1025.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S1013.x + _S1013.y + _S1013.z);
    float _S1026 = _S1023.x + _S1023.y + _S1023.z;
    float _S1027 = _S1015.x + _S1015.y + _S1015.z;
    float _S1028 = _S1021.x + _S1021.y + _S1021.z;
    float _S1029 = _S1017.x + _S1017.y + _S1017.z;
    float _S1030 = _S1019.x + _S1019.y + _S1019.z;
    float _S1031 = - s_diff_fC2_T_0;
    float3  _S1032 = _S998 * _S1010.differential_0;
    float3  _S1033 = (*sh_coeffs_10)[int(8)] * _S1010.differential_0;
    float3  _S1034 = _S996 * _S1010.differential_0;
    float3  _S1035 = (*sh_coeffs_10)[int(7)] * _S1010.differential_0;
    float3  _S1036 = _S995 * _S1010.differential_0;
    float3  _S1037 = (*sh_coeffs_10)[int(6)] * _S1010.differential_0;
    float3  _S1038 = _S997 * _S1010.differential_0;
    float3  _S1039 = (*sh_coeffs_10)[int(5)] * _S1010.differential_0;
    float3  _S1040 = _S999 * _S1010.differential_0;
    float3  _S1041 = (*sh_coeffs_10)[int(4)] * _S1010.differential_0;
    float _S1042 = _S1039.x + _S1039.y + _S1039.z;
    float _S1043 = _S1035.x + _S1035.y + _S1035.z;
    float _S1044 = fTmp1B_10 * _S1026 + x_36 * s_diff_fS2_T_0 + y_13 * _S1031 + 0.54627424478530884f * (_S1041.x + _S1041.y + _S1041.z);
    float _S1045 = fTmp1B_10 * _S1027 + y_13 * s_diff_fS2_T_0 + x_36 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S1033.x + _S1033.y + _S1033.z);
    float _S1046 = y_13 * - _S1045;
    float _S1047 = x_36 * _S1045;
    float _S1048 = z_10 * (1.86588168144226074f * (z_10 * _S1030) + -2.28522896766662598f * (y_13 * _S1028 + x_36 * _S1029) + 0.94617468118667603f * (_S1037.x + _S1037.y + _S1037.z));
    float3  _S1049 = make_float3 (0.48860251903533936f) * _S1010.differential_0;
    float3  _S1050 = - _S1049;
    float3  _S1051 = _S989 * _S1050;
    float3  _S1052 = (*sh_coeffs_10)[int(3)] * _S1050;
    float3  _S1053 = _S991 * _S1049;
    float3  _S1054 = (*sh_coeffs_10)[int(2)] * _S1049;
    float3  _S1055 = _S993 * _S1049;
    float3  _S1056 = (*sh_coeffs_10)[int(1)] * _S1049;
    float _S1057 = (_S1000 * _S1030 + 1.44530570507049561f * (fS1_10 * _S1026 + fC1_10 * _S1027) + -1.09254848957061768f * (y_13 * _S1042 + x_36 * _S1043) + _S1048 + _S1048 + _S1054.x + _S1054.y + _S1054.z) / _S990;
    float _S1058 = _S988 * _S1057;
    float _S1059 = (fTmp0C_10 * _S1028 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S1031 + fTmp0B_10 * _S1042 + _S994 * _S1044 + _S1046 + _S1046 + - (_S1056.x + _S1056.y + _S1056.z)) / _S990;
    float _S1060 = _S988 * _S1059;
    float _S1061 = (fTmp0C_10 * _S1029 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S1043 + 2.0f * (y_13 * _S1044) + _S1047 + _S1047 + _S1052.x + _S1052.y + _S1052.z) / _S990;
    float _S1062 = _S988 * _S1061;
    float _S1063 = _S986 * - _S1057 + _S985 * - _S1059 + _S984 * - _S1061;
    DiffPair_float_0 _S1064;
    (&_S1064)->primal_0 = _S987;
    (&_S1064)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1064, _S1063);
    float _S1065 = _S986 * _S1064.differential_0;
    float _S1066 = _S985 * _S1064.differential_0;
    float _S1067 = _S984 * _S1064.differential_0;
    float3  _S1068 = make_float3 (0.282094806432724f) * _S1010.differential_0;
    float3  _S1069 = make_float3 (_S1062 + _S1067 + _S1067, _S1060 + _S1066 + _S1066, _S1058 + _S1065 + _S1065);
    float3  _S1070 = - - _S1069;
    Matrix<float, 3, 3>  _S1071 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1072;
    (&_S1072)->primal_0 = _S936;
    (&_S1072)->differential_0 = _S1071;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1073;
    (&_S1073)->primal_0 = t_13;
    (&_S1073)->differential_0 = _S1009;
    s_bwd_prop_mul_1(&_S1072, &_S1073, _S1070);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1074 = _S1072;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1075 = _S1073;
    float2  _S1076 = make_float2 (0.0f);
    float2  _S1077 = _S1076;
    *&((&_S1077)->y) = v_conic_0.z;
    float2  _S1078 = _S1076;
    *&((&_S1078)->y) = v_conic_0.y;
    *&((&_S1078)->x) = v_conic_0.x;
    float _S1079 = 0.5f * v_depth_0;
    DiffPair_float_0 _S1080;
    (&_S1080)->primal_0 = _S982;
    (&_S1080)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1080, _S1079);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1081;
    (&_S1081)->primal_0 = mean_c_10;
    (&_S1081)->differential_0 = _S1009;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1082;
    (&_S1082)->primal_0 = mean_c_10;
    (&_S1082)->differential_0 = _S1009;
    s_bwd_prop_dot_0(&_S1081, &_S1082, _S1080.differential_0);
    DiffPair_float_0 _S1083;
    (&_S1083)->primal_0 = _S981;
    (&_S1083)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1083, 0.0f);
    DiffPair_float_0 _S1084;
    (&_S1084)->primal_0 = _S980;
    (&_S1084)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1084, 0.0f);
    DiffPair_float_0 _S1085;
    (&_S1085)->primal_0 = 3.32999992370605469f;
    (&_S1085)->differential_0 = 0.0f;
    DiffPair_float_0 _S1086;
    (&_S1086)->primal_0 = _S979;
    (&_S1086)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1085, &_S1086, 0.0f);
    DiffPair_float_0 _S1087;
    (&_S1087)->primal_0 = _S978;
    (&_S1087)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1087, _S1086.differential_0);
    float _S1088 = 2.0f * _S1087.differential_0;
    DiffPair_float_0 _S1089;
    (&_S1089)->primal_0 = _S977;
    (&_S1089)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1089, _S1088);
    float _S1090 = v_opacity_0 + 254.9999847412109375f * _S1089.differential_0;
    Matrix<float, 2, 2>  _S1091 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1092 = _S1091;
    _S1092[int(1)] = _S1077;
    _S1092[int(0)] = _S1078;
    Matrix<float, 2, 2>  _S1093 = _S1092;
    FixedArray<float3 , 16>  _S1094;
    _S1094[int(0)] = _S1009;
    _S1094[int(1)] = _S1009;
    _S1094[int(2)] = _S1009;
    _S1094[int(3)] = _S1009;
    _S1094[int(4)] = _S1009;
    _S1094[int(5)] = _S1009;
    _S1094[int(6)] = _S1009;
    _S1094[int(7)] = _S1009;
    _S1094[int(8)] = _S1009;
    _S1094[int(9)] = _S1009;
    _S1094[int(10)] = _S1009;
    _S1094[int(11)] = _S1009;
    _S1094[int(12)] = _S1009;
    _S1094[int(13)] = _S1009;
    _S1094[int(14)] = _S1009;
    _S1094[int(15)] = _S1009;
    _S1094[int(7)] = _S1034;
    _S1094[int(0)] = _S1068;
    _S1094[int(1)] = _S1055;
    _S1094[int(2)] = _S1053;
    _S1094[int(3)] = _S1051;
    _S1094[int(4)] = _S1040;
    _S1094[int(5)] = _S1038;
    _S1094[int(6)] = _S1036;
    _S1094[int(15)] = _S1012;
    _S1094[int(8)] = _S1032;
    _S1094[int(9)] = _S1024;
    _S1094[int(10)] = _S1022;
    _S1094[int(11)] = _S1020;
    _S1094[int(12)] = _S1018;
    _S1094[int(13)] = _S1016;
    _S1094[int(14)] = _S1014;
    float3  _S1095 = _S1094[int(0)];
    float3  _S1096 = _S1094[int(1)];
    float3  _S1097 = _S1094[int(2)];
    float3  _S1098 = _S1094[int(3)];
    float3  _S1099 = _S1094[int(4)];
    float3  _S1100 = _S1094[int(5)];
    float3  _S1101 = _S1094[int(6)];
    float3  _S1102 = _S1094[int(7)];
    float3  _S1103 = _S1094[int(8)];
    float3  _S1104 = _S1094[int(9)];
    float3  _S1105 = _S1094[int(10)];
    float3  _S1106 = _S1094[int(11)];
    float3  _S1107 = _S1094[int(12)];
    float3  _S1108 = _S1094[int(13)];
    float3  _S1109 = _S1094[int(14)];
    float3  _S1110 = _S1094[int(15)];
    float3  _S1111 = _S1082.differential_0 + _S1081.differential_0;
    float2  _S1112 = make_float2 (0.0f, _S1083.differential_0);
    float2  _S1113 = make_float2 (_S1084.differential_0, 0.0f);
    float _S1114;
    if(antialiased_10)
    {
        float _S1115 = _S974 * _S1090;
        _S976 = _S969 * _S1090;
        _S1114 = _S1115;
    }
    else
    {
        _S976 = _S1090;
        _S1114 = 0.0f;
    }
    float _S1116 = - (_S976 / _S975);
    DiffPair_float_0 _S1117;
    (&_S1117)->primal_0 = _S972;
    (&_S1117)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1117, _S1116);
    float _S1118 = - _S1117.differential_0;
    float _S1119 = invdet_7 * _S1093.rows[int(1)].y;
    float _S1120 = - (invdet_7 * _S1093.rows[int(1)].x);
    float _S1121 = - (invdet_7 * _S1093.rows[int(0)].y);
    float _S1122 = invdet_7 * _S1093.rows[int(0)].x;
    float _S1123 = - ((_S961 * _S1093.rows[int(1)].y + _S971 * _S1093.rows[int(1)].x + _S970 * _S1093.rows[int(0)].y + _S963 * _S1093.rows[int(0)].x) / _S967);
    DiffPair_float_0 _S1124;
    (&_S1124)->primal_0 = _S968;
    (&_S1124)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1124, _S1114);
    DiffPair_float_0 _S1125;
    (&_S1125)->primal_0 = 0.0f;
    (&_S1125)->differential_0 = 0.0f;
    DiffPair_float_0 _S1126;
    (&_S1126)->primal_0 = _S966;
    (&_S1126)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1125, &_S1126, _S1124.differential_0);
    float _S1127 = _S1126.differential_0 / _S967;
    float s_diff_det_orig_T_0 = det_blur_6 * _S1127;
    float _S1128 = _S1123 + det_orig_11 * - _S1127;
    float _S1129 = - _S1128;
    float _S1130 = _S961 * _S1128;
    float _S1131 = _S963 * _S1128;
    Matrix<float, 2, 2>  _S1132 = _S1091;
    _S1132[int(1)] = _S1112;
    _S1132[int(0)] = _S1113;
    _S962 = _S1132;
    *&(((&_S962)->rows + (int(1)))->y) = 0.0f;
    float _S1133 = _S1122 + _S1130 + _S1132.rows[int(1)].y;
    *&(((&_S962)->rows + (int(0)))->x) = 0.0f;
    float _S1134 = _S1119 + _S1131 + _S1132.rows[int(0)].x;
    float _S1135 = _S1129 + - s_diff_det_orig_T_0;
    float _S1136 = _S1120 + _S957.rows[int(0)].y * _S1135;
    float _S1137 = _S1121 + _S957.rows[int(1)].x * _S1135;
    float _S1138 = _S957.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S1139 = _S1133 + _S957.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S1140 = _S1076;
    *&((&_S1140)->x) = _S1136;
    *&((&_S1140)->y) = _S1139;
    float _S1141 = _S1134 + _S1138;
    float2  _S1142 = _S1076;
    *&((&_S1142)->y) = _S1137;
    *&((&_S1142)->x) = _S1141;
    float _S1143 = _S959 * v_mean2d_0.y;
    float _S1144 = fy_17 * (rz_3 * v_mean2d_0.y);
    float _S1145 = _S958 * v_mean2d_0.x;
    float _S1146 = fx_17 * (rz_3 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S1147 = _S1091;
    _S1147[int(1)] = _S1140;
    _S1147[int(0)] = _S1142;
    Matrix<float, 2, 2>  _S1148 = _S962 + _S1147;
    Matrix<float, 2, 3>  _S1149 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1150;
    (&_S1150)->primal_0 = _S955;
    (&_S1150)->differential_0 = _S1149;
    Matrix<float, 3, 2>  _S1151 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1152;
    (&_S1152)->primal_0 = _S956;
    (&_S1152)->differential_0 = _S1151;
    s_bwd_prop_mul_2(&_S1150, &_S1152, _S1148);
    Matrix<float, 2, 3>  _S1153 = transpose_2(_S1152.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1154;
    (&_S1154)->primal_0 = J_8;
    (&_S1154)->differential_0 = _S1149;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1155;
    (&_S1155)->primal_0 = _S937;
    (&_S1155)->differential_0 = _S1071;
    s_bwd_prop_mul_3(&_S1154, &_S1155, _S1150.differential_0);
    Matrix<float, 2, 3>  _S1156 = _S1153 + _S1154.differential_0;
    float _S1157 = _S954 * _S1156.rows[int(1)].z;
    float s_diff_ty_T_0 = _S953 * (rz2_3 * _S1156.rows[int(1)].z);
    float _S1158 = fy_17 * _S1156.rows[int(1)].y;
    float _S1159 = _S952 * _S1156.rows[int(0)].z;
    float s_diff_tx_T_0 = _S951 * (rz2_3 * _S1156.rows[int(0)].z);
    float _S1160 = fx_17 * _S1156.rows[int(0)].x;
    float _S1161 = mean_c_10.z * s_diff_ty_T_0;
    float _S1162 = _S950 * s_diff_ty_T_0;
    DiffPair_float_0 _S1163;
    (&_S1163)->primal_0 = lim_y_pos_2;
    (&_S1163)->differential_0 = 0.0f;
    DiffPair_float_0 _S1164;
    (&_S1164)->primal_0 = _S949;
    (&_S1164)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1163, &_S1164, _S1161);
    DiffPair_float_0 _S1165;
    (&_S1165)->primal_0 = _S947;
    (&_S1165)->differential_0 = 0.0f;
    DiffPair_float_0 _S1166;
    (&_S1166)->primal_0 = _S948;
    (&_S1166)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1165, &_S1166, _S1164.differential_0);
    float _S1167 = mean_c_10.y * _S1166.differential_0;
    float _S1168 = rz_3 * _S1166.differential_0;
    float _S1169 = mean_c_10.z * s_diff_tx_T_0;
    float _S1170 = _S946 * s_diff_tx_T_0;
    DiffPair_float_0 _S1171;
    (&_S1171)->primal_0 = lim_x_pos_2;
    (&_S1171)->differential_0 = 0.0f;
    DiffPair_float_0 _S1172;
    (&_S1172)->primal_0 = _S945;
    (&_S1172)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1171, &_S1172, _S1169);
    DiffPair_float_0 _S1173;
    (&_S1173)->primal_0 = _S943;
    (&_S1173)->differential_0 = 0.0f;
    DiffPair_float_0 _S1174;
    (&_S1174)->primal_0 = _S944;
    (&_S1174)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1173, &_S1174, _S1172.differential_0);
    float _S1175 = rz_3 * (_S1157 + _S1159);
    float _S1176 = _S1162 + _S1170 + - ((_S1143 + _S1145 + _S1158 + _S1160 + _S1167 + mean_c_10.x * _S1174.differential_0 + _S1175 + _S1175) / _S942);
    float _S1177 = _S1144 + _S1168;
    float _S1178 = _S1146 + rz_3 * _S1174.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1179;
    (&_S1179)->primal_0 = _S935;
    (&_S1179)->differential_0 = _S1071;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1180;
    (&_S1180)->primal_0 = _S936;
    (&_S1180)->differential_0 = _S1071;
    s_bwd_prop_mul_4(&_S1179, &_S1180, _S1155.differential_0);
    Matrix<float, 3, 3>  _S1181 = transpose_0(_S1180.differential_0 + _S1074.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1182;
    (&_S1182)->primal_0 = R_14;
    (&_S1182)->differential_0 = _S1071;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1183;
    (&_S1183)->primal_0 = _S934;
    (&_S1183)->differential_0 = _S1071;
    s_bwd_prop_mul_4(&_S1182, &_S1183, _S1179.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1184;
    (&_S1184)->primal_0 = _S932;
    (&_S1184)->differential_0 = _S1071;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1185;
    (&_S1185)->primal_0 = _S933;
    (&_S1185)->differential_0 = _S1071;
    s_bwd_prop_mul_4(&_S1184, &_S1185, _S1183.differential_0);
    Matrix<float, 3, 3>  _S1186 = _S1184.differential_0 + transpose_0(_S1185.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1187;
    (&_S1187)->primal_0 = _S931;
    (&_S1187)->differential_0 = _S1071;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1188;
    (&_S1188)->primal_0 = S_0;
    (&_S1188)->differential_0 = _S1071;
    s_bwd_prop_mul_4(&_S1187, &_S1188, _S1186);
    Matrix<float, 3, 3>  _S1189 = transpose_0(_S1187.differential_0);
    float _S1190 = 2.0f * - _S1189.rows[int(2)].z;
    float _S1191 = 2.0f * _S1189.rows[int(2)].y;
    float _S1192 = 2.0f * _S1189.rows[int(2)].x;
    float _S1193 = 2.0f * _S1189.rows[int(1)].z;
    float _S1194 = 2.0f * - _S1189.rows[int(1)].y;
    float _S1195 = 2.0f * _S1189.rows[int(1)].x;
    float _S1196 = 2.0f * _S1189.rows[int(0)].z;
    float _S1197 = 2.0f * _S1189.rows[int(0)].y;
    float _S1198 = 2.0f * - _S1189.rows[int(0)].x;
    float _S1199 = - _S1195 + _S1197;
    float _S1200 = _S1192 + - _S1196;
    float _S1201 = - _S1191 + _S1193;
    float _S1202 = _S1191 + _S1193;
    float _S1203 = _S1192 + _S1196;
    float _S1204 = _S1195 + _S1197;
    float _S1205 = quat_14.w * (_S1194 + _S1198);
    float _S1206 = quat_14.z * (_S1190 + _S1198);
    float _S1207 = quat_14.y * (_S1190 + _S1194);
    float _S1208 = quat_14.x * _S1199 + quat_14.z * _S1202 + quat_14.y * _S1203 + _S1205 + _S1205;
    float _S1209 = quat_14.x * _S1200 + quat_14.w * _S1202 + quat_14.y * _S1204 + _S1206 + _S1206;
    float _S1210 = quat_14.x * _S1201 + quat_14.w * _S1203 + quat_14.z * _S1204 + _S1207 + _S1207;
    float _S1211 = quat_14.w * _S1199 + quat_14.z * _S1200 + quat_14.y * _S1201;
    float3  _S1212 = _S1009;
    *&((&_S1212)->z) = _S1188.differential_0.rows[int(2)].z;
    *&((&_S1212)->y) = _S1188.differential_0.rows[int(1)].y;
    *&((&_S1212)->x) = _S1188.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1213;
    (&_S1213)->primal_0 = scale_13;
    (&_S1213)->differential_0 = _S1009;
    s_bwd_prop_exp_1(&_S1213, _S1212);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1214 = _S1213;
    float3  _S1215 = _S1009;
    *&((&_S1215)->z) = _S1176;
    *&((&_S1215)->y) = _S1177;
    *&((&_S1215)->x) = _S1178;
    float3  _S1216 = _S1111 + _S1215;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1217;
    (&_S1217)->primal_0 = R_14;
    (&_S1217)->differential_0 = _S1071;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1218;
    (&_S1218)->primal_0 = mean_11;
    (&_S1218)->differential_0 = _S1009;
    s_bwd_prop_mul_1(&_S1217, &_S1218, _S1216);
    float3  _S1219 = _S1216 + _S1075.differential_0;
    Matrix<float, 3, 3>  _S1220 = _S1181 + _S1182.differential_0 + _S1217.differential_0;
    float4  _S1221 = make_float4 (0.0f);
    *&((&_S1221)->w) = _S1208;
    *&((&_S1221)->z) = _S1209;
    *&((&_S1221)->y) = _S1210;
    *&((&_S1221)->x) = _S1211;
    float4  _S1222 = _S1221;
    float3  _S1223 = _S1218.differential_0 + _S1069;
    *v_mean_0 = _S1223;
    *v_quat_0 = _S1222;
    *v_scale_0 = _S1214.differential_0;
    *v_in_opacity_0 = _S1118;
    (*v_sh_coeffs_0)[int(0)] = _S1095;
    (*v_sh_coeffs_0)[int(1)] = _S1096;
    (*v_sh_coeffs_0)[int(2)] = _S1097;
    (*v_sh_coeffs_0)[int(3)] = _S1098;
    (*v_sh_coeffs_0)[int(4)] = _S1099;
    (*v_sh_coeffs_0)[int(5)] = _S1100;
    (*v_sh_coeffs_0)[int(6)] = _S1101;
    (*v_sh_coeffs_0)[int(7)] = _S1102;
    (*v_sh_coeffs_0)[int(8)] = _S1103;
    (*v_sh_coeffs_0)[int(9)] = _S1104;
    (*v_sh_coeffs_0)[int(10)] = _S1105;
    (*v_sh_coeffs_0)[int(11)] = _S1106;
    (*v_sh_coeffs_0)[int(12)] = _S1107;
    (*v_sh_coeffs_0)[int(13)] = _S1108;
    (*v_sh_coeffs_0)[int(14)] = _S1109;
    (*v_sh_coeffs_0)[int(15)] = _S1110;
    *v_R_1 = _S1220;
    *v_t_1 = _S1219;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S1224;
    DiffPair_float_0 _S1225;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S1226;
    DiffPair_float_0 _S1227;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1228;
    DiffPair_float_0 _S1229;
    DiffPair_float_0 _S1230;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1231;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S1232;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1233;
};

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S1234, float _S1235)
{
    return s_primal_ctx_atan2_0(_S1234, _S1235);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S1236;
    DiffPair_float_0 _S1237;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S1238 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S1236 = _S1238;
    _s_diff_ctx_2->_S1237 = _S1238;
    (&_s_diff_ctx_2->_S1236)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S1236)->differential_0 = 0.0f;
    (&_s_diff_ctx_2->_S1237)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S1237)->differential_0 = 0.0f;
    DiffPair_float_0 _S1239 = *dpdpy_0;
    _s_diff_ctx_2->_S1236 = *dpdpy_0;
    DiffPair_float_0 _S1240 = *dpdpx_0;
    _s_diff_ctx_2->_S1237 = *dpdpx_0;
    float _S1241 = _S1240.primal_0 * _S1240.primal_0 + _S1239.primal_0 * _S1239.primal_0;
    float _S1242 = - _S1239.primal_0 / _S1241 * dpdOut_0;
    float _S1243 = _S1240.primal_0 / _S1241 * dpdOut_0;
    dpdpy_0->primal_0 = _S1239.primal_0;
    dpdpy_0->differential_0 = _S1243;
    dpdpx_0->primal_0 = _S1240.primal_0;
    dpdpx_0->differential_0 = _S1242;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S1244, DiffPair_float_0 * _S1245, float _S1246, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S1247 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S1224 = _S1247;
    _s_diff_ctx_3->_S1225 = _S1247;
    (&_s_diff_ctx_3->_S1224)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S1224)->differential_0 = 0.0f;
    (&_s_diff_ctx_3->_S1225)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S1225)->differential_0 = 0.0f;
    DiffPair_float_0 _S1248 = *_S1244;
    _s_diff_ctx_3->_S1224 = *_S1244;
    DiffPair_float_0 _S1249 = *_S1245;
    _s_diff_ctx_3->_S1225 = *_S1245;
    DiffPair_float_0 _S1250 = _S1248;
    DiffPair_float_0 _S1251 = _S1249;
    s_bwd_prop_d_atan2_Intermediates_0 _S1252;
    (&_S1252)->_S1236 = _S1247;
    (&_S1252)->_S1237 = _S1247;
    s_primal_ctx_d_atan2_0(&_S1250, &_S1251, _S1246, &_S1252);
    *_S1244 = _S1250;
    *_S1245 = _S1251;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1253;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1254;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1255;
    DiffPair_float_0 _S1256;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1257;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1258;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_4)
{
    DiffPair_float_0 _S1259 = { 0.0f, 0.0f };
    _s_diff_ctx_4->_S1258 = _S1259;
    (&_s_diff_ctx_4->_S1258)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S1258)->differential_0 = 0.0f;
    DiffPair_float_0 _S1260 = *dpdpx_1;
    _s_diff_ctx_4->_S1258 = *dpdpx_1;
    float _S1261 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S1260.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S1260.primal_0;
    dpdpx_1->differential_0 = _S1261;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S1262, float _S1263, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_5)
{
    DiffPair_float_0 _S1264 = { 0.0f, 0.0f };
    _s_diff_ctx_5->_S1254 = _S1264;
    (&_s_diff_ctx_5->_S1254)->primal_0 = 0.0f;
    (&_s_diff_ctx_5->_S1254)->differential_0 = 0.0f;
    DiffPair_float_0 _S1265 = *_S1262;
    _s_diff_ctx_5->_S1254 = *_S1262;
    DiffPair_float_0 _S1266 = _S1265;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1267;
    (&_S1267)->_S1258 = _S1264;
    s_primal_ctx_d_sqrt_0(&_S1266, _S1263, &_S1267);
    *_S1262 = _S1266;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_6)
{
    float2  _S1268 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1269 = { _S1268, _S1268 };
    DiffPair_float_0 _S1270 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1271 = { _S1270 };
    _s_diff_ctx_6->_S1255 = _S1269;
    _s_diff_ctx_6->_S1256 = _S1270;
    _s_diff_ctx_6->_S1257 = _S1271;
    (&_s_diff_ctx_6->_S1255)->primal_0 = _S1268;
    (&_s_diff_ctx_6->_S1255)->differential_0 = _S1268;
    (&_s_diff_ctx_6->_S1256)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1256)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1272 = *dpdpx_2;
    _s_diff_ctx_6->_S1255 = *dpdpx_2;
    float _S1273 = _S1272.primal_0.x;
    float _S1274 = _S1272.primal_0.y;
    DiffPair_float_0 _S1275;
    (&_S1275)->primal_0 = _S1273 * _S1273 + _S1274 * _S1274;
    (&_S1275)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S1275, dp_s_dOut_0, &_s_diff_ctx_6->_S1257);
    _s_diff_ctx_6->_S1256 = _S1275;
    float _S1276 = _S1272.primal_0.y * _S1275.differential_0;
    float _S1277 = _S1276 + _S1276;
    float _S1278 = _S1272.primal_0.x * _S1275.differential_0;
    float _S1279 = _S1278 + _S1278;
    float2  _S1280 = _S1268;
    *&((&_S1280)->y) = _S1277;
    *&((&_S1280)->x) = _S1279;
    dpdpx_2->primal_0 = _S1272.primal_0;
    dpdpx_2->differential_0 = _S1280;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1281, float _S1282, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_7)
{
    float2  _S1283 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1284 = { _S1283, _S1283 };
    _s_diff_ctx_7->_S1253 = _S1284;
    (&_s_diff_ctx_7->_S1253)->primal_0 = _S1283;
    (&_s_diff_ctx_7->_S1253)->differential_0 = _S1283;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1285 = *_S1281;
    _s_diff_ctx_7->_S1253 = *_S1281;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1286 = _S1285;
    DiffPair_float_0 _S1287 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1288 = { _S1287 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1289;
    (&_S1289)->_S1255 = _S1284;
    (&_S1289)->_S1256 = _S1287;
    (&_S1289)->_S1257 = _S1288;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1286, _S1282, &_S1289);
    *_S1281 = _S1286;
    return;
}

inline __device__ bool s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float4  dpintrins_0, FixedArray<float, 10>  * dpdist_coeffs_0, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_float_0 _S1290 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1291 = { _S1290, _S1290 };
    _s_diff_ctx_8->_S1226 = _S1290;
    _s_diff_ctx_8->_S1227 = _S1290;
    _s_diff_ctx_8->_S1228 = _S1291;
    _s_diff_ctx_8->_S1229 = _S1290;
    _s_diff_ctx_8->_S1230 = _S1290;
    _s_diff_ctx_8->_S1231 = _S1291;
    (&_s_diff_ctx_8->_S1226)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1226)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1227)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1227)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1229)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1229)->differential_0 = 0.0f;
    (&_s_diff_ctx_8->_S1230)->primal_0 = 0.0f;
    (&_s_diff_ctx_8->_S1230)->differential_0 = 0.0f;
    float2  _S1292 = make_float2 (0.0f);
    float2  _S1293 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S1294 = length_0(_S1293);
    float _S1295 = dpmean3d_0.z;
    float _S1296 = s_primal_ctx_atan2_0(_S1294, _S1295);
    float k_4;
    if(_S1296 < 0.00100000004749745f)
    {
        k_4 = (1.0f - _S1296 * _S1296 / 3.0f) / _S1295;
    }
    else
    {
        k_4 = _S1296 / _S1294;
    }
    float2  _S1297 = _S1293 * make_float2 (k_4);
    float u_38 = _S1297.x;
    float v_38 = _S1297.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float _S1298 = 2.0f * (*dpdist_coeffs_0)[int(4)];
    float _S1299 = 2.0f * (*dpdist_coeffs_0)[int(5)];
    float2  _S1300 = _S1297 * make_float2 (1.0f + r2_38 * ((*dpdist_coeffs_0)[int(0)] + r2_38 * ((*dpdist_coeffs_0)[int(1)] + r2_38 * ((*dpdist_coeffs_0)[int(2)] + r2_38 * (*dpdist_coeffs_0)[int(3)])))) + make_float2 (_S1298 * u_38 * v_38 + (*dpdist_coeffs_0)[int(5)] * (r2_38 + 2.0f * u_38 * u_38) + (*dpdist_coeffs_0)[int(6)] * r2_38, _S1299 * u_38 * v_38 + (*dpdist_coeffs_0)[int(4)] * (r2_38 + 2.0f * v_38 * v_38) + (*dpdist_coeffs_0)[int(7)] * r2_38);
    float2  _S1301 = _S1300 + make_float2 ((*dpdist_coeffs_0)[int(8)] * _S1300.x + (*dpdist_coeffs_0)[int(9)] * _S1300.y, 0.0f);
    float fx_18 = dpintrins_0.x;
    float fy_18 = dpintrins_0.y;
    float2  _S1302 = make_float2 (fx_18 * _S1301.x + dpintrins_0.z, fy_18 * _S1301.y + dpintrins_0.w);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (0.0f);
    float _S1303 = s_primal_ctx_s_primal_ctx_atan2_0(_S1294, _S1295);
    bool _S1304 = _S1303 < 0.00100000004749745f;
    float _S1305;
    float _S1306;
    float _S1307;
    if(_S1304)
    {
        float _S1308 = 1.0f - _S1303 * _S1303 / 3.0f;
        float _S1309 = _S1295 * _S1295;
        k_4 = _S1308 / _S1295;
        _S1305 = 0.0f;
        _S1306 = _S1309;
        _S1307 = _S1308;
    }
    else
    {
        float _S1310 = _S1294 * _S1294;
        k_4 = _S1303 / _S1294;
        _S1305 = _S1310;
        _S1306 = 0.0f;
        _S1307 = 0.0f;
    }
    float2  _S1311 = make_float2 (k_4);
    float2  _S1312 = _S1293 * make_float2 (k_4);
    float u_39 = _S1312.x;
    float v_39 = _S1312.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float _S1313 = (*dpdist_coeffs_0)[int(2)] + r2_39 * (*dpdist_coeffs_0)[int(3)];
    float _S1314 = (*dpdist_coeffs_0)[int(1)] + r2_39 * _S1313;
    float _S1315 = (*dpdist_coeffs_0)[int(0)] + r2_39 * _S1314;
    float2  _S1316 = make_float2 (fx_18, 0.0f) + make_float2 ((*dpdist_coeffs_0)[int(8)] * fx_18, (*dpdist_coeffs_0)[int(9)] * fx_18);
    float2  _S1317 = _S1312 * _S1316;
    float _S1318 = (*dpdist_coeffs_0)[int(4)] * _S1316.y;
    float _S1319 = (*dpdist_coeffs_0)[int(5)] * _S1316.x;
    float _S1320 = _S1317.x + _S1317.y;
    float _S1321 = r2_39 * _S1320;
    float _S1322 = r2_39 * _S1321;
    float _S1323 = (*dpdist_coeffs_0)[int(7)] * _S1316.y + _S1318 + (*dpdist_coeffs_0)[int(6)] * _S1316.x + _S1319 + _S1315 * _S1320 + _S1314 * _S1321 + _S1313 * _S1322 + (*dpdist_coeffs_0)[int(3)] * (r2_39 * _S1322);
    float _S1324 = v_39 * _S1323;
    float _S1325 = u_39 * _S1323;
    float2  _S1326 = make_float2 (1.0f + r2_39 * _S1315) * _S1316 + make_float2 (_S1299 * (v_39 * _S1316.y) + 2.0f * u_39 * _S1319 + 2.0f * (u_39 * _S1319) + _S1298 * (v_39 * _S1316.x) + _S1325 + _S1325, 2.0f * v_39 * _S1318 + 2.0f * (v_39 * _S1318) + _S1299 * u_39 * _S1316.y + _S1298 * u_39 * _S1316.x + _S1324 + _S1324);
    float2  _S1327 = _S1293 * _S1326;
    float2  _S1328 = _S1311 * _S1326;
    float _S1329 = _S1327.x + _S1327.y;
    if(_S1304)
    {
        float _S1330 = _S1329 / _S1306;
        float _S1331 = _S1307 * - _S1330;
        float _S1332 = _S1303 * (0.3333333432674408f * - (_S1295 * _S1330));
        k_4 = _S1332 + _S1332;
        _S1305 = _S1331;
        _S1306 = 0.0f;
    }
    else
    {
        float _S1333 = _S1329 / _S1305;
        float _S1334 = _S1303 * - _S1333;
        k_4 = _S1294 * _S1333;
        _S1305 = 0.0f;
        _S1306 = _S1334;
    }
    DiffPair_float_0 _S1335;
    (&_S1335)->primal_0 = _S1294;
    (&_S1335)->differential_0 = 0.0f;
    DiffPair_float_0 _S1336;
    (&_S1336)->primal_0 = _S1295;
    (&_S1336)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1335, &_S1336, k_4, &_s_diff_ctx_8->_S1228);
    _s_diff_ctx_8->_S1226 = _S1335;
    _s_diff_ctx_8->_S1227 = _S1336;
    float _S1337 = _S1336.differential_0 + _S1305;
    float _S1338 = _S1335.differential_0 + _S1306;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1339;
    (&_S1339)->primal_0 = _S1293;
    (&_S1339)->differential_0 = _S1292;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1340 = { _S1292, _S1292 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1341;
    (&_S1341)->_S1253 = _S1340;
    s_primal_ctx_s_bwd_length_impl_0(&_S1339, _S1338, &_S1341);
    float2  _S1342 = _S1339.differential_0 + _S1328;
    float3  _S1343 = make_float3 (_S1342.x, _S1342.y, _S1337);
    Matrix<float, 2, 3>  _S1344 = J_9;
    _S1344[int(0)] = _S1343;
    if(_S1304)
    {
        float _S1345 = 1.0f - _S1303 * _S1303 / 3.0f;
        float _S1346 = _S1295 * _S1295;
        k_4 = _S1345 / _S1295;
        _S1305 = 0.0f;
        _S1306 = _S1346;
        _S1307 = _S1345;
    }
    else
    {
        float _S1347 = _S1294 * _S1294;
        k_4 = _S1303 / _S1294;
        _S1305 = _S1347;
        _S1306 = 0.0f;
        _S1307 = 0.0f;
    }
    float2  _S1348 = make_float2 (k_4);
    float2  _S1349 = _S1293 * make_float2 (k_4);
    float u_40 = _S1349.x;
    float v_40 = _S1349.y;
    float r2_40 = u_40 * u_40 + v_40 * v_40;
    float _S1350 = (*dpdist_coeffs_0)[int(2)] + r2_40 * (*dpdist_coeffs_0)[int(3)];
    float _S1351 = (*dpdist_coeffs_0)[int(1)] + r2_40 * _S1350;
    float _S1352 = (*dpdist_coeffs_0)[int(0)] + r2_40 * _S1351;
    float2  _S1353 = make_float2 (0.0f, fy_18);
    float2  _S1354 = _S1349 * _S1353;
    float _S1355 = (*dpdist_coeffs_0)[int(4)] * fy_18;
    float _S1356 = _S1354.x + _S1354.y;
    float _S1357 = r2_40 * _S1356;
    float _S1358 = r2_40 * _S1357;
    float _S1359 = (*dpdist_coeffs_0)[int(7)] * fy_18 + _S1355 + _S1352 * _S1356 + _S1351 * _S1357 + _S1350 * _S1358 + (*dpdist_coeffs_0)[int(3)] * (r2_40 * _S1358);
    float _S1360 = v_40 * _S1359;
    float _S1361 = u_40 * _S1359;
    float2  _S1362 = make_float2 (1.0f + r2_40 * _S1352) * _S1353 + make_float2 (_S1299 * (v_40 * fy_18) + _S1361 + _S1361, 2.0f * v_40 * _S1355 + 2.0f * (v_40 * _S1355) + _S1299 * u_40 * fy_18 + _S1360 + _S1360);
    float2  _S1363 = _S1293 * _S1362;
    float2  _S1364 = _S1348 * _S1362;
    float _S1365 = _S1363.x + _S1363.y;
    if(_S1304)
    {
        float _S1366 = _S1365 / _S1306;
        float _S1367 = _S1307 * - _S1366;
        float _S1368 = _S1303 * (0.3333333432674408f * - (_S1295 * _S1366));
        k_4 = _S1368 + _S1368;
        _S1305 = _S1367;
        _S1306 = 0.0f;
    }
    else
    {
        float _S1369 = _S1365 / _S1305;
        float _S1370 = _S1303 * - _S1369;
        k_4 = _S1294 * _S1369;
        _S1305 = 0.0f;
        _S1306 = _S1370;
    }
    DiffPair_float_0 _S1371;
    (&_S1371)->primal_0 = _S1294;
    (&_S1371)->differential_0 = 0.0f;
    DiffPair_float_0 _S1372;
    (&_S1372)->primal_0 = _S1295;
    (&_S1372)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1371, &_S1372, k_4, &_s_diff_ctx_8->_S1231);
    _s_diff_ctx_8->_S1229 = _S1371;
    _s_diff_ctx_8->_S1230 = _S1372;
    float _S1373 = _S1372.differential_0 + _S1305;
    float _S1374 = _S1371.differential_0 + _S1306;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1375;
    (&_S1375)->primal_0 = _S1293;
    (&_S1375)->differential_0 = _S1292;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1376;
    (&_S1376)->_S1253 = _S1340;
    s_primal_ctx_s_bwd_length_impl_0(&_S1375, _S1374, &_S1376);
    float2  _S1377 = _S1375.differential_0 + _S1364;
    float3  _S1378 = make_float3 (_S1377.x, _S1377.y, _S1373);
    _S1344[int(1)] = _S1378;
    *dpcov2d_0 = s_primal_ctx_mul_4(s_primal_ctx_mul_3(_S1344, dpcov3d_0), transpose_1(_S1344));
    *dpmean2d_0 = _S1302;
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
    DiffPair_1 _S1379 = *dpdpx_3;
    float _S1380 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_9->_S1258)->primal_0);
    float _S1381 = s_primal_ctx_sqrt_0(_S1380);
    float _S1382 = 0.5f / _S1381 * (*dpdpx_3).differential_0.differential_0;
    float _S1383 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S1381 * _S1381));
    DiffPair_float_0 _S1384;
    (&_S1384)->primal_0 = _S1380;
    (&_S1384)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1384, _S1383);
    DiffPair_float_0 _S1385;
    (&_S1385)->primal_0 = 1.00000001168609742e-07f;
    (&_S1385)->differential_0 = 0.0f;
    DiffPair_float_0 _S1386;
    (&_S1386)->primal_0 = (&_s_diff_ctx_9->_S1258)->primal_0;
    (&_S1386)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1385, &_S1386, _S1384.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S1386.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S1382;
    dpdpx_3->primal_0 = _S1379.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S1387, DiffPair_float_0 * _S1388, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_1 _S1389 = *_S1387;
    DiffPair_float_0 _S1390 = _s_diff_ctx_10->_S1254;
    DiffPair_float_0 _S1391 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S1392;
    (&_S1392)->_S1258 = _S1391;
    s_primal_ctx_d_sqrt_0(&_S1390, (*_S1388).primal_0, &_S1392);
    DiffPair_float_0 _S1393 = { (*_S1387).differential_0.primal_0, (*_S1387).differential_0.differential_0 };
    DiffPair_1 _S1394;
    (&_S1394)->primal_0 = _s_diff_ctx_10->_S1254;
    (&_S1394)->differential_0 = _S1393;
    DiffPair_float_0 _S1395;
    (&_S1395)->primal_0 = (*_S1388).primal_0;
    (&_S1395)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1396 = _S1392;
    s_bwd_prop_d_sqrt_0(&_S1394, &_S1395, &_S1396);
    DiffPair_float_0 _S1397 = { _S1394.differential_0.primal_0, _S1394.differential_0.differential_0 };
    _S1388->primal_0 = (*_S1388).primal_0;
    _S1388->differential_0 = _S1395.differential_0;
    _S1387->primal_0 = _S1389.primal_0;
    _S1387->differential_0 = _S1397;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S1398, float _s_dOut_3)
{
    DiffPair_float_0 _S1399;
    (&_S1399)->primal_0 = (*_S1398).primal_0;
    (&_S1399)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1399, _s_dOut_3);
    _S1398->primal_0 = (*_S1398).primal_0;
    _S1398->differential_0 = _S1399.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_0 _S1400 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->x) * *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->x) + *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->y) * *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->y);
    DiffPair_float_0 _S1401 = { len_0, 0.0f };
    float2  _S1402 = make_float2 (0.0f);
    float _S1403 = (*dpdpx_5).differential_0.differential_0.x;
    float _S1404 = _S1403 + _S1403;
    float _S1405 = (&_s_diff_ctx_11->_S1256)->differential_0 * _S1404;
    float _S1406 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S1407 = (&_s_diff_ctx_11->_S1256)->differential_0 * _S1406;
    DiffPair_float_0 _S1408 = { 0.0f, *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->x) * _S1404 + *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->y) * _S1406 };
    DiffPair_1 _S1409;
    (&_S1409)->primal_0 = _S1401;
    (&_S1409)->differential_0 = _S1408;
    DiffPair_float_0 _S1410;
    (&_S1410)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S1410)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S1409, &_S1410, &_s_diff_ctx_11->_S1257);
    DiffPair_float_0 _S1411;
    (&_S1411)->primal_0 = len_0;
    (&_S1411)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1411, 0.0f);
    float _S1412 = _S1409.differential_0.primal_0 + _S1411.differential_0;
    float _S1413 = *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->y) * _S1412;
    float _S1414 = _S1407 + _S1413 + _S1413;
    float _S1415 = *&((&(&_s_diff_ctx_11->_S1255)->primal_0)->x) * _S1412;
    float _S1416 = _S1405 + _S1415 + _S1415;
    float2  _S1417 = _S1402;
    *&((&_S1417)->y) = _S1414;
    *&((&_S1417)->x) = _S1416;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S1400.differential_0.primal_0 + _S1417, _S1402 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S1410.differential_0;
    dpdpx_5->primal_0 = _S1400.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S1418 = (*dpdpx_7).primal_0.x;
    float _S1419 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = _S1418 * _S1418 + _S1419 * _S1419;
    (&_S1420)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1420, _s_dOut_4);
    float _S1421 = (*dpdpx_7).primal_0.y * _S1420.differential_0;
    float _S1422 = _S1421 + _S1421;
    float _S1423 = (*dpdpx_7).primal_0.x * _S1420.differential_0;
    float _S1424 = _S1423 + _S1423;
    float2  _S1425 = make_float2 (0.0f);
    *&((&_S1425)->y) = _S1422;
    *&((&_S1425)->x) = _S1424;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S1425;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S1426, DiffPair_float_0 * _S1427, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_0 _S1428 = *_S1426;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1429 = _s_diff_ctx_12->_S1253;
    float2  _S1430 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1431 = { _S1430, _S1430 };
    DiffPair_float_0 _S1432 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1433 = { _S1432 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1434;
    (&_S1434)->_S1255 = _S1431;
    (&_S1434)->_S1256 = _S1432;
    (&_S1434)->_S1257 = _S1433;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1429, (*_S1427).primal_0, &_S1434);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1435 = { (*_S1426).differential_0.primal_0, (*_S1426).differential_0.differential_0 };
    DiffPair_0 _S1436;
    (&_S1436)->primal_0 = _s_diff_ctx_12->_S1253;
    (&_S1436)->differential_0 = _S1435;
    DiffPair_float_0 _S1437;
    (&_S1437)->primal_0 = (*_S1427).primal_0;
    (&_S1437)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1438 = _S1434;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S1436, &_S1437, &_S1438);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1439;
    (&_S1439)->primal_0 = (&_s_diff_ctx_12->_S1253)->primal_0;
    (&_S1439)->differential_0 = _S1430;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S1439, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1440 = { _S1436.differential_0.primal_0 + _S1439.differential_0, _S1436.differential_0.differential_0 };
    _S1427->primal_0 = (*_S1427).primal_0;
    _S1427->differential_0 = _S1437.differential_0;
    _S1426->primal_0 = _S1428.primal_0;
    _S1426->differential_0 = _S1440;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_1 _S1441 = *dpdpy_1;
    DiffPair_1 _S1442 = *dpdpx_8;
    float _S1443 = - (&_s_diff_ctx_13->_S1236)->primal_0;
    float _S1444 = (&_s_diff_ctx_13->_S1237)->primal_0 * (&_s_diff_ctx_13->_S1237)->primal_0 + (&_s_diff_ctx_13->_S1236)->primal_0 * (&_s_diff_ctx_13->_S1236)->primal_0;
    float _S1445 = _S1444 * _S1444;
    float _S1446 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S1445;
    float _S1447 = (&_s_diff_ctx_13->_S1237)->primal_0 * - _S1446;
    float _S1448 = (&_s_diff_ctx_13->_S1236)->primal_0 * _S1447;
    float _S1449 = (&_s_diff_ctx_13->_S1237)->primal_0 * _S1447;
    float _S1450 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S1445;
    float _S1451 = _S1443 * - _S1450;
    float _S1452 = (&_s_diff_ctx_13->_S1236)->primal_0 * _S1451;
    float _S1453 = (&_s_diff_ctx_13->_S1237)->primal_0 * _S1451;
    DiffPair_float_0 dpdpx_9 = { _S1453 + _S1453 + ((*dpdpx_8).differential_0.primal_0 + (_S1449 + _S1449 + _S1444 * _S1446)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S1448 + _S1448 + (*dpdpy_1).differential_0.primal_0 + _S1452 + _S1452 + - (_S1444 * _S1450), 0.0f };
    float _S1454 = (&_s_diff_ctx_13->_S1237)->primal_0 / _S1444 * (*dpdpy_1).differential_0.differential_0 + _S1443 / _S1444 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S1454;
    dpdpy_1->primal_0 = _S1441.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S1442.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S1455, DiffPair_1 * _S1456, DiffPair_float_0 * _S1457, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_14)
{
    DiffPair_1 _S1458 = *_S1455;
    DiffPair_1 _S1459 = *_S1456;
    DiffPair_float_0 _S1460 = _s_diff_ctx_14->_S1224;
    DiffPair_float_0 _S1461 = _s_diff_ctx_14->_S1225;
    DiffPair_float_0 _S1462 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S1463;
    (&_S1463)->_S1236 = _S1462;
    (&_S1463)->_S1237 = _S1462;
    s_primal_ctx_d_atan2_0(&_S1460, &_S1461, (*_S1457).primal_0, &_S1463);
    DiffPair_float_0 _S1464 = { (*_S1456).differential_0.primal_0, (*_S1456).differential_0.differential_0 };
    DiffPair_float_0 _S1465 = { (*_S1455).differential_0.primal_0, (*_S1455).differential_0.differential_0 };
    DiffPair_1 _S1466;
    (&_S1466)->primal_0 = _s_diff_ctx_14->_S1224;
    (&_S1466)->differential_0 = _S1465;
    DiffPair_1 _S1467;
    (&_S1467)->primal_0 = _s_diff_ctx_14->_S1225;
    (&_S1467)->differential_0 = _S1464;
    DiffPair_float_0 _S1468;
    (&_S1468)->primal_0 = (*_S1457).primal_0;
    (&_S1468)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S1469 = _S1463;
    s_bwd_prop_d_atan2_0(&_S1466, &_S1467, &_S1468, &_S1469);
    DiffPair_float_0 _S1470 = { _S1467.differential_0.primal_0, _S1467.differential_0.differential_0 };
    DiffPair_float_0 _S1471 = { _S1466.differential_0.primal_0, _S1466.differential_0.differential_0 };
    _S1457->primal_0 = (*_S1457).primal_0;
    _S1457->differential_0 = _S1468.differential_0;
    _S1455->primal_0 = _S1458.primal_0;
    _S1455->differential_0 = _S1471;
    _S1456->primal_0 = _S1459.primal_0;
    _S1456->differential_0 = _S1470;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S1472, DiffPair_float_0 * _S1473, float _s_dOut_5)
{
    DiffPair_float_0 _S1474;
    (&_S1474)->primal_0 = (*_S1472).primal_0;
    (&_S1474)->differential_0 = 0.0f;
    DiffPair_float_0 _S1475;
    (&_S1475)->primal_0 = (*_S1473).primal_0;
    (&_S1475)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1474, &_S1475, _s_dOut_5);
    _S1473->primal_0 = (*_S1473).primal_0;
    _S1473->differential_0 = _S1475.differential_0;
    _S1472->primal_0 = (*_S1472).primal_0;
    _S1472->differential_0 = _S1474.differential_0;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpintrins_1, DiffPair_arrayx3Cfloatx2C10x3E_0 * dpdist_coeffs_1, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_15)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1476 = *dpcov3d_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1477 = *dpintrins_1;
    FixedArray<float, 10>  _S1478 = dpdist_coeffs_1->primal_0;
    float2  _S1479 = make_float2 (0.0f);
    float2  _S1480 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S1481 = length_0(_S1480);
    float _S1482 = (*dpmean3d_1).primal_0.z;
    float _S1483 = s_primal_ctx_atan2_0(_S1481, _S1482);
    bool _S1484 = _S1483 < 0.00100000004749745f;
    float k_5;
    float _S1485;
    float _S1486;
    float _S1487;
    if(_S1484)
    {
        float _S1488 = 1.0f - _S1483 * _S1483 / 3.0f;
        float _S1489 = _S1482 * _S1482;
        k_5 = _S1488 / _S1482;
        _S1485 = 0.0f;
        _S1486 = _S1489;
        _S1487 = _S1488;
    }
    else
    {
        float _S1490 = _S1481 * _S1481;
        k_5 = _S1483 / _S1481;
        _S1485 = _S1490;
        _S1486 = 0.0f;
        _S1487 = 0.0f;
    }
    float2  _S1491 = make_float2 (k_5);
    float2  _S1492 = _S1480 * make_float2 (k_5);
    float u_41 = _S1492.x;
    float v_41 = _S1492.y;
    float r2_41 = u_41 * u_41 + v_41 * v_41;
    float _S1493 = _S1478[int(2)] + r2_41 * _S1478[int(3)];
    float _S1494 = _S1478[int(1)] + r2_41 * _S1493;
    float _S1495 = _S1478[int(0)] + r2_41 * _S1494;
    float radial_0 = 1.0f + r2_41 * _S1495;
    float2  _S1496 = make_float2 (radial_0);
    float _S1497 = 2.0f * _S1478[int(4)];
    float _S1498 = _S1497 * u_41;
    float _S1499 = 2.0f * u_41;
    float _S1500 = r2_41 + _S1499 * u_41;
    float _S1501 = 2.0f * _S1478[int(5)];
    float _S1502 = _S1501 * u_41;
    float _S1503 = 2.0f * v_41;
    float _S1504 = r2_41 + _S1503 * v_41;
    float2  _S1505 = _S1492 * make_float2 (radial_0) + make_float2 (_S1498 * v_41 + _S1478[int(5)] * _S1500 + _S1478[int(6)] * r2_41, _S1502 * v_41 + _S1478[int(4)] * _S1504 + _S1478[int(7)] * r2_41);
    float _S1506 = _S1505.x;
    float _S1507 = _S1505.y;
    float2  _S1508 = _S1505 + make_float2 (_S1478[int(8)] * _S1506 + _S1478[int(9)] * _S1507, 0.0f);
    float fx_19 = _S1477.primal_0.x;
    float fy_19 = _S1477.primal_0.y;
    float _S1509 = _S1508.x;
    float _S1510 = _S1508.y;
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float _S1511 = s_primal_ctx_s_primal_ctx_atan2_0(_S1481, _S1482);
    bool _S1512 = _S1511 < 0.00100000004749745f;
    float _S1513;
    float _S1514;
    float _S1515;
    if(_S1512)
    {
        float _S1516 = 1.0f - _S1511 * _S1511 / 3.0f;
        float _S1517 = _S1482 * _S1482;
        k_5 = _S1516 / _S1482;
        _S1513 = 0.0f;
        _S1514 = _S1517;
        _S1515 = _S1516;
    }
    else
    {
        float _S1518 = _S1481 * _S1481;
        k_5 = _S1511 / _S1481;
        _S1513 = _S1518;
        _S1514 = 0.0f;
        _S1515 = 0.0f;
    }
    float2  _S1519 = make_float2 (k_5);
    float2  _S1520 = _S1480 * make_float2 (k_5);
    float u_42 = _S1520.x;
    float v_42 = _S1520.y;
    float r2_42 = u_42 * u_42 + v_42 * v_42;
    float _S1521 = _S1478[int(2)] + r2_42 * _S1478[int(3)];
    float _S1522 = _S1478[int(1)] + r2_42 * _S1521;
    float _S1523 = _S1478[int(0)] + r2_42 * _S1522;
    float2  _S1524 = make_float2 (1.0f + r2_42 * _S1523);
    float _S1525 = _S1497 * u_42;
    float _S1526 = 2.0f * u_42;
    float _S1527 = _S1501 * u_42;
    float _S1528 = 2.0f * v_42;
    float2  _S1529 = make_float2 (fx_19, 0.0f) + make_float2 (_S1478[int(8)] * fx_19, _S1478[int(9)] * fx_19);
    float2  _S1530 = _S1520 * _S1529;
    float _S1531 = _S1478[int(4)] * _S1529.y;
    float _S1532 = v_42 * _S1529.y;
    float _S1533 = _S1478[int(5)] * _S1529.x;
    float _S1534 = v_42 * _S1529.x;
    float _S1535 = _S1530.x + _S1530.y;
    float _S1536 = r2_42 * _S1535;
    float _S1537 = r2_42 * _S1536;
    float _S1538 = r2_42 * _S1537;
    float _S1539 = _S1478[int(7)] * _S1529.y + _S1531 + _S1478[int(6)] * _S1529.x + _S1533 + _S1523 * _S1535 + _S1522 * _S1536 + _S1521 * _S1537 + _S1478[int(3)] * _S1538;
    float _S1540 = v_42 * _S1539;
    float _S1541 = u_42 * _S1539;
    float2  _S1542 = _S1524 * _S1529 + make_float2 (_S1501 * _S1532 + _S1526 * _S1533 + 2.0f * (u_42 * _S1533) + _S1497 * _S1534 + _S1541 + _S1541, _S1528 * _S1531 + 2.0f * (v_42 * _S1531) + _S1527 * _S1529.y + _S1525 * _S1529.x + _S1540 + _S1540);
    float2  _S1543 = _S1480 * _S1542;
    float2  _S1544 = _S1519 * _S1542;
    float _S1545 = _S1543.x + _S1543.y;
    float k_6;
    float _S1546;
    float _S1547;
    float _S1548;
    float _S1549;
    float _S1550;
    float _S1551;
    float _S1552;
    float _S1553;
    if(_S1512)
    {
        float _S1554 = _S1545 / _S1514;
        float _S1555 = _S1514 * _S1514;
        float _S1556 = - _S1554;
        float _S1557 = _S1515 * _S1556;
        float _S1558 = 0.3333333432674408f * - (_S1482 * _S1554);
        float _S1559 = _S1511 * _S1558;
        k_5 = _S1559 + _S1559;
        k_6 = _S1557;
        _S1546 = 0.0f;
        _S1547 = 0.0f;
        _S1548 = 0.0f;
        _S1549 = 0.0f;
        _S1550 = _S1558;
        _S1551 = _S1554;
        _S1552 = _S1556;
        _S1553 = _S1555;
    }
    else
    {
        float _S1560 = _S1545 / _S1513;
        float _S1561 = _S1513 * _S1513;
        float _S1562 = - _S1560;
        float _S1563 = _S1511 * _S1562;
        k_5 = _S1481 * _S1560;
        k_6 = 0.0f;
        _S1546 = _S1563;
        _S1547 = _S1560;
        _S1548 = _S1562;
        _S1549 = _S1561;
        _S1550 = 0.0f;
        _S1551 = 0.0f;
        _S1552 = 0.0f;
        _S1553 = 0.0f;
    }
    DiffPair_float_0 _S1564 = { _S1481, 0.0f };
    DiffPair_float_0 _S1565 = { _S1482, 0.0f };
    float _S1566 = (&_s_diff_ctx_15->_S1227)->differential_0 + k_6;
    float _S1567 = (&_s_diff_ctx_15->_S1226)->differential_0 + _S1546;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1568 = { _S1480, _S1479 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1569;
    (&_S1569)->primal_0 = _S1480;
    (&_S1569)->differential_0 = _S1479;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1570 = { _S1479, _S1479 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1571;
    (&_S1571)->_S1253 = _S1570;
    s_primal_ctx_s_bwd_length_impl_0(&_S1569, _S1567, &_S1571);
    float2  _S1572 = _S1569.differential_0 + _S1544;
    float3  _S1573 = make_float3 (_S1572.x, _S1572.y, _S1566);
    Matrix<float, 2, 3>  _S1574 = J_10;
    _S1574[int(0)] = _S1573;
    float _S1575;
    float _S1576;
    if(_S1512)
    {
        float _S1577 = 1.0f - _S1511 * _S1511 / 3.0f;
        float _S1578 = _S1482 * _S1482;
        k_6 = _S1577 / _S1482;
        _S1546 = 0.0f;
        _S1575 = _S1578;
        _S1576 = _S1577;
    }
    else
    {
        float _S1579 = _S1481 * _S1481;
        k_6 = _S1511 / _S1481;
        _S1546 = _S1579;
        _S1575 = 0.0f;
        _S1576 = 0.0f;
    }
    float2  _S1580 = make_float2 (k_6);
    float2  _S1581 = _S1480 * make_float2 (k_6);
    float u_43 = _S1581.x;
    float v_43 = _S1581.y;
    float r2_43 = u_43 * u_43 + v_43 * v_43;
    float _S1582 = _S1478[int(2)] + r2_43 * _S1478[int(3)];
    float _S1583 = _S1478[int(1)] + r2_43 * _S1582;
    float _S1584 = _S1478[int(0)] + r2_43 * _S1583;
    float2  _S1585 = make_float2 (1.0f + r2_43 * _S1584);
    float _S1586 = _S1497 * u_43;
    float _S1587 = 2.0f * u_43;
    float _S1588 = _S1501 * u_43;
    float _S1589 = 2.0f * v_43;
    float2  _S1590 = make_float2 (0.0f, fy_19);
    float2  _S1591 = _S1581 * _S1590;
    float _S1592 = _S1478[int(4)] * fy_19;
    float _S1593 = v_43 * fy_19;
    float _S1594 = _S1591.x + _S1591.y;
    float _S1595 = r2_43 * _S1594;
    float _S1596 = r2_43 * _S1595;
    float _S1597 = r2_43 * _S1596;
    float _S1598 = _S1478[int(7)] * fy_19 + _S1592 + _S1584 * _S1594 + _S1583 * _S1595 + _S1582 * _S1596 + _S1478[int(3)] * _S1597;
    float _S1599 = v_43 * _S1598;
    float _S1600 = u_43 * _S1598;
    float2  _S1601 = _S1585 * _S1590 + make_float2 (_S1501 * _S1593 + _S1600 + _S1600, _S1589 * _S1592 + 2.0f * (v_43 * _S1592) + _S1588 * fy_19 + _S1599 + _S1599);
    float2  _S1602 = _S1480 * _S1601;
    float2  _S1603 = _S1580 * _S1601;
    float _S1604 = _S1602.x + _S1602.y;
    float _S1605;
    float _S1606;
    float _S1607;
    float _S1608;
    float _S1609;
    float _S1610;
    float _S1611;
    float _S1612;
    float _S1613;
    if(_S1512)
    {
        float _S1614 = _S1604 / _S1575;
        float _S1615 = _S1575 * _S1575;
        float _S1616 = - _S1614;
        float _S1617 = _S1576 * _S1616;
        float _S1618 = 0.3333333432674408f * - (_S1482 * _S1614);
        float _S1619 = _S1511 * _S1618;
        k_6 = _S1619 + _S1619;
        _S1605 = _S1617;
        _S1606 = 0.0f;
        _S1607 = 0.0f;
        _S1608 = 0.0f;
        _S1609 = 0.0f;
        _S1610 = _S1618;
        _S1611 = _S1614;
        _S1612 = _S1616;
        _S1613 = _S1615;
    }
    else
    {
        float _S1620 = _S1604 / _S1546;
        float _S1621 = _S1546 * _S1546;
        float _S1622 = - _S1620;
        float _S1623 = _S1511 * _S1622;
        k_6 = _S1481 * _S1620;
        _S1605 = 0.0f;
        _S1606 = _S1623;
        _S1607 = _S1620;
        _S1608 = _S1622;
        _S1609 = _S1621;
        _S1610 = 0.0f;
        _S1611 = 0.0f;
        _S1612 = 0.0f;
        _S1613 = 0.0f;
    }
    float _S1624 = (&_s_diff_ctx_15->_S1230)->differential_0 + _S1605;
    float _S1625 = (&_s_diff_ctx_15->_S1229)->differential_0 + _S1606;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1626;
    (&_S1626)->primal_0 = _S1480;
    (&_S1626)->differential_0 = _S1479;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1627;
    (&_S1627)->_S1253 = _S1570;
    s_primal_ctx_s_bwd_length_impl_0(&_S1626, _S1625, &_S1627);
    float2  _S1628 = _S1626.differential_0 + _S1603;
    float3  _S1629 = make_float3 (_S1628.x, _S1628.y, _S1624);
    _S1574[int(1)] = _S1629;
    Matrix<float, 3, 2>  _S1630 = transpose_1(_S1574);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1631;
    (&_S1631)->primal_0 = s_primal_ctx_mul_3(_S1574, _S1476.primal_0);
    (&_S1631)->differential_0 = J_10;
    Matrix<float, 3, 2>  _S1632 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1633;
    (&_S1633)->primal_0 = _S1630;
    (&_S1633)->differential_0 = _S1632;
    s_bwd_prop_mul_2(&_S1631, &_S1633, dpcov2d_1);
    Matrix<float, 2, 3>  _S1634 = transpose_2(_S1633.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1635;
    (&_S1635)->primal_0 = _S1574;
    (&_S1635)->differential_0 = J_10;
    Matrix<float, 3, 3>  _S1636 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1637;
    (&_S1637)->primal_0 = _S1476.primal_0;
    (&_S1637)->differential_0 = _S1636;
    s_bwd_prop_mul_3(&_S1635, &_S1637, _S1631.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1638 = _S1637;
    Matrix<float, 2, 3>  _S1639 = _S1634 + _S1635.differential_0;
    float2  _S1640 = _S1479;
    *&((&_S1640)->y) = _S1639.rows[int(1)].y;
    *&((&_S1640)->x) = _S1639.rows[int(1)].x;
    float2  _S1641 = _S1640;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1642 = { _S1479, _S1640 };
    DiffPair_0 _S1643;
    (&_S1643)->primal_0 = _S1568;
    (&_S1643)->differential_0 = _S1642;
    DiffPair_float_0 _S1644;
    (&_S1644)->primal_0 = _S1625;
    (&_S1644)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1645 = _S1627;
    s_bwd_prop_s_bwd_length_impl_0(&_S1643, &_S1644, &_S1645);
    DiffPair_0 _S1646 = _S1643;
    DiffPair_float_0 _S1647 = _S1644;
    DiffPair_float_0 _S1648 = { 0.0f, _S1639.rows[int(1)].z };
    DiffPair_float_0 _S1649 = { 0.0f, _S1644.differential_0 };
    DiffPair_1 _S1650;
    (&_S1650)->primal_0 = _S1564;
    (&_S1650)->differential_0 = _S1649;
    DiffPair_1 _S1651;
    (&_S1651)->primal_0 = _S1565;
    (&_S1651)->differential_0 = _S1648;
    DiffPair_float_0 _S1652;
    (&_S1652)->primal_0 = k_6;
    (&_S1652)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1650, &_S1651, &_S1652, &_s_diff_ctx_15->_S1231);
    DiffPair_1 _S1653 = _S1650;
    DiffPair_1 _S1654 = _S1651;
    DiffPair_float_0 _S1655 = _S1652;
    if(_S1512)
    {
        float _S1656 = _S1655.differential_0 + _S1655.differential_0;
        float _S1657 = _S1610 * _S1656;
        float _S1658 = - (0.3333333432674408f * (_S1511 * _S1656));
        float _S1659 = _S1612 * _S1639.rows[int(1)].z;
        float _S1660 = (_S1482 * _S1658 + - (_S1576 * _S1639.rows[int(1)].z)) / _S1613;
        float _S1661 = _S1604 * - _S1660;
        float _S1662 = _S1611 * _S1658 + _S1654.differential_0.primal_0;
        k_6 = _S1575 * _S1660;
        _S1605 = 0.0f;
        _S1606 = _S1661;
        _S1607 = _S1659;
        _S1608 = _S1653.differential_0.primal_0;
        _S1609 = _S1657;
        _S1610 = _S1662;
    }
    else
    {
        float _S1663 = _S1608 * _S1647.differential_0;
        float _S1664 = (_S1481 * _S1655.differential_0 + - (_S1511 * _S1647.differential_0)) / _S1609;
        float _S1665 = _S1604 * - _S1664;
        float _S1666 = _S1607 * _S1655.differential_0 + _S1653.differential_0.primal_0;
        k_6 = _S1546 * _S1664;
        _S1605 = _S1665;
        _S1606 = 0.0f;
        _S1607 = 0.0f;
        _S1608 = _S1666;
        _S1609 = _S1663;
        _S1610 = _S1654.differential_0.primal_0;
    }
    float2  _S1667 = _S1580 * _S1641;
    float2  _S1668 = _S1601 * _S1641;
    float2  _S1669 = _S1479;
    *&((&_S1669)->y) = k_6;
    *&((&_S1669)->x) = k_6;
    float2  _S1670 = _S1601 * _S1669;
    float2  _S1671 = _S1667 + _S1480 * _S1669;
    float _S1672 = _S1671.x;
    float _S1673 = _S1672 + _S1672;
    float _S1674 = _S1598 * _S1673;
    float _S1675 = _S1671.y + _S1671.y;
    float _S1676 = _S1598 * _S1675;
    float _S1677 = u_43 * _S1673 + v_43 * _S1675;
    float _S1678 = _S1478[int(3)] * _S1677;
    float _S1679 = _S1597 * _S1677;
    float _S1680 = _S1596 * _S1677;
    float _S1681 = _S1596 * _S1678;
    float _S1682 = _S1595 * _S1677;
    float _S1683 = _S1582 * _S1677 + r2_43 * _S1678;
    float _S1684 = _S1595 * _S1683;
    float _S1685 = _S1594 * _S1677;
    float _S1686 = _S1583 * _S1677 + r2_43 * _S1683;
    float _S1687 = _S1594 * _S1686;
    float _S1688 = _S1584 * _S1677 + r2_43 * _S1686;
    float _S1689 = v_43 * (_S1497 * _S1671.x);
    float _S1690 = _S1586 * _S1671.y;
    float _S1691 = _S1478[int(5)] * (_S1677 + u_43 * (2.0f * _S1671.x) + _S1587 * _S1671.x);
    float _S1692 = _S1478[int(6)] * _S1677;
    float _S1693 = _S1501 * _S1671.x;
    float _S1694 = _S1593 * _S1671.x;
    float _S1695 = v_43 * _S1693;
    float _S1696 = fy_19 * _S1693;
    float _S1697 = _S1588 * _S1671.y;
    float _S1698 = fy_19 * _S1671.y;
    float _S1699 = 2.0f * _S1671.y;
    float _S1700 = _S1592 * _S1699;
    float _S1701 = _S1592 * _S1671.y;
    float _S1702 = _S1677 + v_43 * _S1699 + _S1589 * _S1671.y;
    float _S1703 = _S1478[int(4)] * _S1702;
    float _S1704 = fy_19 * _S1702;
    float _S1705 = _S1478[int(7)] * _S1677;
    float _S1706 = fy_19 * _S1677;
    float2  _S1707 = _S1585 * _S1671;
    float2  _S1708 = _S1590 * _S1671;
    float2  _S1709 = _S1479;
    *&((&_S1709)->y) = _S1688;
    *&((&_S1709)->x) = _S1688;
    float2  _S1710 = _S1590 * _S1709;
    float _S1711 = _S1695 + _S1697 + _S1703 + _S1705;
    float _S1712 = _S1689 + _S1690 + _S1691 + _S1692;
    float2  _S1713 = _S1707 + _S1581 * _S1709;
    float2  _S1714 = _S1479;
    *&((&_S1714)->y) = _S1711;
    *&((&_S1714)->x) = _S1712;
    float2  _S1715 = _S1713 + _S1714;
    float _S1716 = _S1708.x + _S1708.y;
    float _S1717 = _S1685 + r2_43 * _S1716;
    float _S1718 = _S1682 + r2_43 * _S1717;
    float _S1719 = _S1680 + r2_43 * _S1718;
    float _S1720 = _S1681 + _S1684 + _S1687 + _S1584 * _S1716 + _S1583 * _S1717 + _S1582 * _S1718 + _S1478[int(3)] * _S1719;
    float _S1721 = v_43 * _S1720;
    float _S1722 = u_43 * _S1720;
    float2  _S1723 = _S1670 + _S1646.differential_0.primal_0;
    float2  _S1724 = _S1710 + make_float2 (_S1674 + _S1501 * _S1698 + _S1722 + _S1722, _S1676 + _S1696 + _S1700 + 2.0f * _S1701 + _S1721 + _S1721);
    float _S1725 = _S1694 + u_43 * _S1698;
    float _S1726 = _S1679 + r2_43 * _S1719;
    float2  _S1727 = _S1480 * _S1724;
    float _S1728 = _S1668.x + _S1668.y + _S1727.x + _S1727.y;
    float2  _S1729 = _S1580 * _S1724 + _S1723;
    if(_S1512)
    {
        float _S1730 = _S1482 * _S1606;
        float _S1731 = _S1728 / _S1575;
        float _S1732 = _S1511 * (0.3333333432674408f * - (_S1607 + _S1482 * _S1731));
        float _S1733 = _S1730 + _S1730 + _S1576 * - _S1731 + _S1610;
        k_6 = _S1732 + _S1732 + _S1609;
        _S1546 = _S1733;
        _S1575 = _S1608;
    }
    else
    {
        float _S1734 = _S1481 * _S1605;
        float _S1735 = _S1728 / _S1546;
        float _S1736 = _S1734 + _S1734 + _S1511 * - _S1735 + _S1608;
        k_6 = _S1481 * _S1735 + _S1609;
        _S1546 = _S1610;
        _S1575 = _S1736;
    }
    DiffPair_float_0 _S1737;
    (&_S1737)->primal_0 = _S1481;
    (&_S1737)->differential_0 = 0.0f;
    DiffPair_float_0 _S1738;
    (&_S1738)->primal_0 = _S1482;
    (&_S1738)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1737, &_S1738, k_6);
    float _S1739 = _S1738.differential_0 + _S1546;
    float _S1740 = _S1737.differential_0 + _S1575;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1741;
    (&_S1741)->primal_0 = _S1480;
    (&_S1741)->differential_0 = _S1479;
    s_bwd_length_impl_1(&_S1741, _S1740);
    float2  _S1742 = _S1741.differential_0 + _S1729;
    float3  _S1743 = make_float3 (_S1742.x, _S1742.y, _S1739);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1744;
    (&_S1744)->primal_0 = _S1480;
    (&_S1744)->differential_0 = _S1479;
    s_bwd_length_impl_1(&_S1744, 0.0f);
    float3  _S1745 = _S1743 + make_float3 (_S1744.differential_0.x, _S1744.differential_0.y, 0.0f);
    float2  _S1746 = _S1479;
    *&((&_S1746)->y) = _S1639.rows[int(0)].y;
    *&((&_S1746)->x) = _S1639.rows[int(0)].x;
    float2  _S1747 = _S1746;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1748 = { _S1479, _S1746 };
    DiffPair_0 _S1749;
    (&_S1749)->primal_0 = _S1568;
    (&_S1749)->differential_0 = _S1748;
    DiffPair_float_0 _S1750;
    (&_S1750)->primal_0 = _S1567;
    (&_S1750)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1751 = _S1571;
    s_bwd_prop_s_bwd_length_impl_0(&_S1749, &_S1750, &_S1751);
    DiffPair_0 _S1752 = _S1749;
    DiffPair_float_0 _S1753 = _S1750;
    DiffPair_float_0 _S1754 = { 0.0f, _S1639.rows[int(0)].z };
    DiffPair_float_0 _S1755 = { 0.0f, _S1750.differential_0 };
    DiffPair_1 _S1756;
    (&_S1756)->primal_0 = _S1564;
    (&_S1756)->differential_0 = _S1755;
    DiffPair_1 _S1757;
    (&_S1757)->primal_0 = _S1565;
    (&_S1757)->differential_0 = _S1754;
    DiffPair_float_0 _S1758;
    (&_S1758)->primal_0 = k_5;
    (&_S1758)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1756, &_S1757, &_S1758, &_s_diff_ctx_15->_S1228);
    DiffPair_1 _S1759 = _S1756;
    DiffPair_1 _S1760 = _S1757;
    DiffPair_float_0 _S1761 = _S1758;
    if(_S1512)
    {
        float _S1762 = _S1761.differential_0 + _S1761.differential_0;
        float _S1763 = _S1550 * _S1762;
        float _S1764 = - (0.3333333432674408f * (_S1511 * _S1762));
        float _S1765 = _S1552 * _S1639.rows[int(0)].z;
        float _S1766 = (_S1482 * _S1764 + - (_S1515 * _S1639.rows[int(0)].z)) / _S1553;
        float _S1767 = _S1545 * - _S1766;
        float _S1768 = _S1551 * _S1764 + _S1760.differential_0.primal_0;
        k_5 = _S1514 * _S1766;
        k_6 = 0.0f;
        _S1546 = _S1767;
        _S1547 = _S1765;
        _S1548 = _S1759.differential_0.primal_0;
        _S1549 = _S1763;
        _S1550 = _S1768;
    }
    else
    {
        float _S1769 = _S1548 * _S1753.differential_0;
        float _S1770 = (_S1481 * _S1761.differential_0 + - (_S1511 * _S1753.differential_0)) / _S1549;
        float _S1771 = _S1545 * - _S1770;
        float _S1772 = _S1547 * _S1761.differential_0 + _S1759.differential_0.primal_0;
        k_5 = _S1513 * _S1770;
        k_6 = _S1771;
        _S1546 = 0.0f;
        _S1547 = 0.0f;
        _S1548 = _S1772;
        _S1549 = _S1769;
        _S1550 = _S1760.differential_0.primal_0;
    }
    float2  _S1773 = _S1519 * _S1747;
    float2  _S1774 = _S1542 * _S1747;
    float2  _S1775 = _S1479;
    *&((&_S1775)->y) = k_5;
    *&((&_S1775)->x) = k_5;
    float2  _S1776 = _S1542 * _S1775;
    float2  _S1777 = _S1773 + _S1480 * _S1775;
    float _S1778 = _S1777.x;
    float _S1779 = _S1778 + _S1778;
    float _S1780 = _S1539 * _S1779;
    float _S1781 = _S1777.y + _S1777.y;
    float _S1782 = _S1539 * _S1781;
    float _S1783 = u_42 * _S1779 + v_42 * _S1781;
    float _S1784 = _S1478[int(3)] * _S1783;
    float _S1785 = _S1538 * _S1783;
    float _S1786 = _S1537 * _S1783;
    float _S1787 = _S1537 * _S1784;
    float _S1788 = _S1536 * _S1783;
    float _S1789 = _S1521 * _S1783 + r2_42 * _S1784;
    float _S1790 = _S1536 * _S1789;
    float _S1791 = _S1535 * _S1783;
    float _S1792 = _S1522 * _S1783 + r2_42 * _S1789;
    float _S1793 = _S1535 * _S1792;
    float _S1794 = _S1523 * _S1783 + r2_42 * _S1792;
    float _S1795 = _S1497 * _S1777.x;
    float _S1796 = _S1534 * _S1777.x;
    float _S1797 = v_42 * _S1795;
    float _S1798 = _S1529.x * _S1795;
    float _S1799 = _S1525 * _S1777.y;
    float _S1800 = _S1529.x * _S1777.y;
    float _S1801 = 2.0f * _S1777.x;
    float _S1802 = _S1533 * _S1801;
    float _S1803 = _S1533 * _S1777.x;
    float _S1804 = _S1783 + u_42 * _S1801 + _S1526 * _S1777.x;
    float _S1805 = _S1478[int(5)] * _S1804;
    float _S1806 = _S1529.x * _S1804;
    float _S1807 = _S1478[int(6)] * _S1783;
    float _S1808 = _S1529.x * _S1783;
    float _S1809 = _S1501 * _S1777.x;
    float _S1810 = _S1532 * _S1777.x;
    float _S1811 = v_42 * _S1809;
    float _S1812 = _S1529.y * _S1809;
    float _S1813 = _S1527 * _S1777.y;
    float _S1814 = _S1529.y * _S1777.y;
    float _S1815 = 2.0f * _S1777.y;
    float _S1816 = _S1531 * _S1815;
    float _S1817 = _S1531 * _S1777.y;
    float _S1818 = _S1783 + v_42 * _S1815 + _S1528 * _S1777.y;
    float _S1819 = _S1478[int(4)] * _S1818;
    float _S1820 = _S1529.y * _S1818;
    float _S1821 = _S1478[int(7)] * _S1783;
    float _S1822 = _S1529.y * _S1783;
    float2  _S1823 = _S1524 * _S1777;
    float2  _S1824 = _S1529 * _S1777;
    float2  _S1825 = _S1479;
    *&((&_S1825)->y) = _S1794;
    *&((&_S1825)->x) = _S1794;
    float2  _S1826 = _S1529 * _S1825;
    float _S1827 = _S1811 + _S1813 + _S1819 + _S1821;
    float _S1828 = _S1797 + _S1799 + _S1805 + _S1807;
    float2  _S1829 = _S1823 + _S1520 * _S1825;
    float2  _S1830 = _S1479;
    *&((&_S1830)->y) = _S1827;
    *&((&_S1830)->x) = _S1828;
    float2  _S1831 = _S1829 + _S1830;
    float _S1832 = fx_19 * _S1831.x;
    float _S1833 = fx_19 * _S1831.y;
    float _S1834 = _S1824.x + _S1824.y;
    float _S1835 = _S1791 + r2_42 * _S1834;
    float _S1836 = _S1788 + r2_42 * _S1835;
    float _S1837 = _S1786 + r2_42 * _S1836;
    float _S1838 = _S1787 + _S1790 + _S1793 + _S1523 * _S1834 + _S1522 * _S1835 + _S1521 * _S1836 + _S1478[int(3)] * _S1837;
    float _S1839 = v_42 * _S1838;
    float _S1840 = u_42 * _S1838;
    float2  _S1841 = _S1776 + _S1752.differential_0.primal_0;
    float _S1842 = _S1822 + _S1706;
    float _S1843 = _S1796 + u_42 * _S1800;
    float _S1844 = _S1478[int(8)] * _S1831.x + _S1478[int(9)] * _S1831.y + _S1831.x;
    float _S1845 = _S1836 + _S1718;
    float _S1846 = _S1835 + _S1717;
    float2  _S1847 = _S1826 + make_float2 (_S1780 + _S1802 + _S1501 * _S1814 + 2.0f * _S1803 + _S1497 * _S1800 + _S1840 + _S1840, _S1782 + _S1798 + _S1812 + _S1816 + 2.0f * _S1817 + _S1839 + _S1839);
    float _S1848 = _S1810 + u_42 * _S1814 + _S1725;
    float _S1849 = _S1785 + r2_42 * _S1837 + _S1726;
    float _S1850 = _S1837 + _S1719;
    float _S1851 = _S1820 + _S1704;
    float2  _S1852 = _S1480 * _S1847;
    float _S1853 = _S1774.x + _S1774.y + _S1852.x + _S1852.y;
    float2  _S1854 = _S1519 * _S1847 + _S1841;
    if(_S1512)
    {
        float _S1855 = _S1482 * _S1546;
        float _S1856 = _S1853 / _S1514;
        float _S1857 = _S1511 * (0.3333333432674408f * - (_S1547 + _S1482 * _S1856));
        float _S1858 = _S1855 + _S1855 + _S1515 * - _S1856 + _S1550;
        k_5 = _S1857 + _S1857 + _S1549;
        _S1513 = _S1858;
        _S1514 = _S1548;
    }
    else
    {
        float _S1859 = _S1481 * k_6;
        float _S1860 = _S1853 / _S1513;
        float _S1861 = _S1859 + _S1859 + _S1511 * - _S1860 + _S1548;
        k_5 = _S1481 * _S1860 + _S1549;
        _S1513 = _S1550;
        _S1514 = _S1861;
    }
    DiffPair_float_0 _S1862;
    (&_S1862)->primal_0 = _S1481;
    (&_S1862)->differential_0 = 0.0f;
    DiffPair_float_0 _S1863;
    (&_S1863)->primal_0 = _S1482;
    (&_S1863)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1862, &_S1863, k_5);
    float _S1864 = _S1863.differential_0 + _S1513;
    float _S1865 = _S1862.differential_0 + _S1514;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1866;
    (&_S1866)->primal_0 = _S1480;
    (&_S1866)->differential_0 = _S1479;
    s_bwd_length_impl_1(&_S1866, _S1865);
    float2  _S1867 = _S1866.differential_0 + _S1854;
    float3  _S1868 = make_float3 (_S1867.x, _S1867.y, _S1864);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1869;
    (&_S1869)->primal_0 = _S1480;
    (&_S1869)->differential_0 = _S1479;
    s_bwd_length_impl_1(&_S1869, 0.0f);
    float _S1870 = fx_19 * dpmean2d_1.x;
    float2  _S1871 = make_float2 (_S1870, fy_19 * dpmean2d_1.y) + make_float2 (_S1478[int(8)] * _S1870, _S1478[int(9)] * _S1870);
    float2  _S1872 = _S1492 * _S1871;
    float2  _S1873 = _S1496 * _S1871;
    float _S1874 = _S1478[int(4)] * _S1871.y;
    float _S1875 = v_41 * _S1871.y;
    float _S1876 = _S1478[int(5)] * _S1871.x;
    float _S1877 = v_41 * _S1871.x;
    float _S1878 = _S1872.x + _S1872.y;
    float _S1879 = r2_41 * _S1878;
    float _S1880 = r2_41 * _S1879;
    float _S1881 = r2_41 * _S1880;
    float _S1882 = _S1478[int(7)] * _S1871.y + _S1874 + _S1478[int(6)] * _S1871.x + _S1876 + _S1495 * _S1878 + _S1494 * _S1879 + _S1493 * _S1880 + _S1478[int(3)] * _S1881;
    float _S1883 = v_41 * _S1882;
    float _S1884 = u_41 * _S1882;
    float _S1885 = _S1503 * _S1874 + 2.0f * (v_41 * _S1874) + _S1502 * _S1871.y + _S1498 * _S1871.x + _S1883 + _S1883;
    float _S1886 = _S1501 * _S1875 + _S1499 * _S1876 + 2.0f * (u_41 * _S1876) + _S1497 * _S1877 + _S1884 + _S1884;
    float _S1887 = _S1507 * _S1870 + _S1833;
    float _S1888 = _S1506 * _S1870 + _S1832;
    float _S1889 = r2_41 * _S1871.y + _S1842;
    float _S1890 = r2_41 * _S1871.x + _S1808;
    float _S1891 = 2.0f * (u_41 * _S1875 + _S1848) + _S1500 * _S1871.x + _S1806;
    float _S1892 = _S1504 * _S1871.y + 2.0f * (u_41 * _S1877 + _S1843) + _S1851;
    float _S1893 = r2_41 * _S1881 + _S1849;
    float _S1894 = _S1881 + _S1850;
    float _S1895 = _S1880 + _S1845;
    float _S1896 = _S1879 + _S1846;
    float3  _S1897 = _S1868 + make_float3 (_S1869.differential_0.x, _S1869.differential_0.y, 0.0f) + _S1745;
    float4  _S1898 = make_float4 (_S1509 * dpmean2d_1.x + _S1844, _S1510 * dpmean2d_1.y + _S1715.y, dpmean2d_1.x, dpmean2d_1.y);
    FixedArray<float, 10>  _S1899;
    _S1899[int(0)] = 0.0f;
    _S1899[int(1)] = 0.0f;
    _S1899[int(2)] = 0.0f;
    _S1899[int(3)] = 0.0f;
    _S1899[int(4)] = 0.0f;
    _S1899[int(5)] = 0.0f;
    _S1899[int(6)] = 0.0f;
    _S1899[int(7)] = 0.0f;
    _S1899[int(8)] = 0.0f;
    _S1899[int(9)] = 0.0f;
    _S1899[int(9)] = _S1887;
    _S1899[int(8)] = _S1888;
    _S1899[int(7)] = _S1889;
    _S1899[int(6)] = _S1890;
    _S1899[int(5)] = _S1891;
    _S1899[int(4)] = _S1892;
    _S1899[int(3)] = _S1893;
    _S1899[int(2)] = _S1894;
    _S1899[int(1)] = _S1895;
    _S1899[int(0)] = _S1896;
    FixedArray<float, 10>  _S1900 = {
        _S1899[int(0)], _S1899[int(1)], _S1899[int(2)], _S1899[int(3)], _S1899[int(4)], _S1899[int(5)], _S1899[int(6)], _S1899[int(7)], _S1899[int(8)], _S1899[int(9)]
    };
    float2  _S1901 = _S1873 + make_float2 (_S1886, _S1885);
    float2  _S1902 = _S1480 * _S1901;
    float2  _S1903 = _S1491 * _S1901;
    float _S1904 = _S1902.x + _S1902.y;
    if(_S1484)
    {
        float _S1905 = _S1904 / _S1486;
        float _S1906 = _S1487 * - _S1905;
        float _S1907 = _S1483 * (0.3333333432674408f * - (_S1482 * _S1905));
        k_5 = _S1907 + _S1907;
        _S1485 = _S1906;
        _S1486 = 0.0f;
    }
    else
    {
        float _S1908 = _S1904 / _S1485;
        float _S1909 = _S1483 * - _S1908;
        k_5 = _S1481 * _S1908;
        _S1485 = 0.0f;
        _S1486 = _S1909;
    }
    DiffPair_float_0 _S1910;
    (&_S1910)->primal_0 = _S1481;
    (&_S1910)->differential_0 = 0.0f;
    DiffPair_float_0 _S1911;
    (&_S1911)->primal_0 = _S1482;
    (&_S1911)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1910, &_S1911, k_5);
    float _S1912 = _S1911.differential_0 + _S1485;
    float _S1913 = _S1910.differential_0 + _S1486;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1914;
    (&_S1914)->primal_0 = _S1480;
    (&_S1914)->differential_0 = _S1479;
    s_bwd_length_impl_1(&_S1914, _S1913);
    float2  _S1915 = _S1914.differential_0 + _S1903;
    dpdist_coeffs_1->primal_0 = dpdist_coeffs_1->primal_0;
    dpdist_coeffs_1->differential_0 = _S1900;
    dpintrins_1->primal_0 = (*dpintrins_1).primal_0;
    dpintrins_1->differential_0 = _S1898;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1638.differential_0;
    float3  _S1916 = _S1897 + make_float3 (_S1915.x, _S1915.y, _S1912);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1916;
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_15, float3  t_14, float fx_20, float fy_20, float cx_15, float cy_15, FixedArray<float, 10>  * dist_coeffs_27, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1917 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1918 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1919 = { _S1918, _S1918 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1920 = { _S1918, _S1918, _S1919, _S1918, _S1918, _S1919 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1921;
    (&_S1921)->_S1232 = _S1917;
    (&_S1921)->_S1233 = _S1920;
    float3  mean_c_11 = s_primal_ctx_mul_1(R_15, mean_12) + t_14;
    float4  intrins_11 = make_float4 (fx_20, fy_20, cx_15, cy_15);
    float3  _S1922 = s_primal_ctx_exp_0(scale_14);
    float _S1923 = quat_15.y;
    float x2_15 = _S1923 * _S1923;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_26 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S1924 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_26), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_26), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1922.x, 0.0f, 0.0f, 0.0f, _S1922.y, 0.0f, 0.0f, 0.0f, _S1922.z);
    Matrix<float, 3, 3>  _S1925 = s_primal_ctx_mul_2(_S1924, S_1);
    Matrix<float, 3, 3>  _S1926 = transpose_0(_S1925);
    Matrix<float, 3, 3>  _S1927 = s_primal_ctx_mul_2(_S1925, _S1926);
    Matrix<float, 3, 3>  _S1928 = s_primal_ctx_mul_2(R_15, _S1927);
    Matrix<float, 3, 3>  _S1929 = transpose_0(R_15);
    Matrix<float, 3, 3>  _S1930 = s_primal_ctx_mul_2(_S1928, _S1929);
    Matrix<float, 2, 2>  _S1931 = _S1917;
    float2  _S1932 = make_float2 (0.0f);
    float2  _S1933 = _S1932;
    bool _S1934 = s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1930, intrins_11, dist_coeffs_27, &_S1931, &_S1933, &(&_S1921)->_S1233);
    (&_S1921)->_S1232 = _S1931;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1935 = _S1921;
    float _S1936 = _S1921._S1232.rows[int(0)].y * _S1921._S1232.rows[int(1)].x;
    float det_orig_12 = _S1921._S1232.rows[int(0)].x * _S1921._S1232.rows[int(1)].y - _S1936;
    float _S1937 = _S1921._S1232.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1938 = _S1921._S1232;
    *&(((&_S1938)->rows + (int(0)))->x) = _S1937;
    float _S1939 = _S1921._S1232.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1938)->rows + (int(1)))->y) = _S1939;
    Matrix<float, 2, 2>  _S1940 = _S1938;
    Matrix<float, 2, 2>  _S1941 = _S1938;
    float det_blur_7 = _S1937 * _S1939 - _S1936;
    float _S1942 = det_orig_12 / det_blur_7;
    float _S1943 = det_blur_7 * det_blur_7;
    float _S1944 = s_primal_ctx_max_0(0.0f, _S1942);
    float _S1945 = s_primal_ctx_sqrt_0(_S1944);
    float invdet_8 = 1.0f / det_blur_7;
    float _S1946 = - _S1921._S1232.rows[int(0)].y;
    float _S1947 = - _S1921._S1232.rows[int(1)].x;
    float _S1948 = - in_opacity_11;
    float _S1949 = 1.0f + s_primal_ctx_exp_1(_S1948);
    float _S1950 = 1.0f / _S1949;
    float _S1951 = _S1949 * _S1949;
    float _S1952;
    if(antialiased_11)
    {
        _S1952 = _S1950 * _S1945;
    }
    else
    {
        _S1952 = _S1950;
    }
    float _S1953 = _S1952 / 0.00392156885936856f;
    float _S1954 = 2.0f * s_primal_ctx_log_0(_S1953);
    float _S1955 = s_primal_ctx_sqrt_0(_S1954);
    float _S1956 = _S1940.rows[int(0)].x;
    float _S1957 = _S1941.rows[int(1)].y;
    float _S1958 = s_primal_ctx_dot_0(mean_c_11, mean_c_11) + 9.99999997475242708e-07f;
    float3  _S1959 = mean_12 - - s_primal_ctx_mul_1(_S1929, t_14);
    float _S1960 = _S1959.x;
    float _S1961 = _S1959.y;
    float _S1962 = _S1959.z;
    float _S1963 = _S1960 * _S1960 + _S1961 * _S1961 + _S1962 * _S1962;
    float _S1964 = s_primal_ctx_sqrt_0(_S1963);
    float x_37 = _S1960 / _S1964;
    float3  _S1965 = make_float3 (x_37);
    float _S1966 = _S1964 * _S1964;
    float y_14 = _S1961 / _S1964;
    float z_11 = _S1962 / _S1964;
    float3  _S1967 = make_float3 (z_11);
    float _S1968 = - y_14;
    float3  _S1969 = make_float3 (_S1968);
    float z2_27 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_37 * x_37 - y_14 * y_14;
    float _S1970 = 2.0f * x_37;
    float fS1_11 = _S1970 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_27 - 0.31539157032966614f;
    float3  _S1971 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_37;
    float3  _S1972 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S1973 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S1974 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S1975 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_27 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S1976 = 1.86588168144226074f * z2_27 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S1976;
    float3  _S1977 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_37;
    float3  _S1978 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S1979 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S1980 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S1981 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_37 * fC1_11 - y_14 * fS1_11);
    float3  _S1982 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_37 * fS1_11 + y_14 * fC1_11);
    float3  _S1983 = make_float3 (pSH9_1);
    float3  _S1984 = make_float3 (0.0f);
    float3  _S1985 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1986;
    (&_S1986)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1968) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_37) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S1986)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1987;
    (&_S1987)->primal_0 = _S1984;
    (&_S1987)->differential_0 = _S1985;
    s_bwd_prop_max_0(&_S1986, &_S1987, v_rgb_1);
    float3  _S1988 = _S1982 * _S1986.differential_0;
    float3  _S1989 = (*sh_coeffs_11)[int(15)] * _S1986.differential_0;
    float3  _S1990 = _S1980 * _S1986.differential_0;
    float3  _S1991 = (*sh_coeffs_11)[int(14)] * _S1986.differential_0;
    float3  _S1992 = _S1978 * _S1986.differential_0;
    float3  _S1993 = (*sh_coeffs_11)[int(13)] * _S1986.differential_0;
    float3  _S1994 = _S1977 * _S1986.differential_0;
    float3  _S1995 = (*sh_coeffs_11)[int(12)] * _S1986.differential_0;
    float3  _S1996 = _S1979 * _S1986.differential_0;
    float3  _S1997 = (*sh_coeffs_11)[int(11)] * _S1986.differential_0;
    float3  _S1998 = _S1981 * _S1986.differential_0;
    float3  _S1999 = (*sh_coeffs_11)[int(10)] * _S1986.differential_0;
    float3  _S2000 = _S1983 * _S1986.differential_0;
    float3  _S2001 = (*sh_coeffs_11)[int(9)] * _S1986.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S2001.x + _S2001.y + _S2001.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1989.x + _S1989.y + _S1989.z);
    float _S2002 = _S1999.x + _S1999.y + _S1999.z;
    float _S2003 = _S1991.x + _S1991.y + _S1991.z;
    float _S2004 = _S1997.x + _S1997.y + _S1997.z;
    float _S2005 = _S1993.x + _S1993.y + _S1993.z;
    float _S2006 = _S1995.x + _S1995.y + _S1995.z;
    float _S2007 = - s_diff_fC2_T_1;
    float3  _S2008 = _S1974 * _S1986.differential_0;
    float3  _S2009 = (*sh_coeffs_11)[int(8)] * _S1986.differential_0;
    float3  _S2010 = _S1972 * _S1986.differential_0;
    float3  _S2011 = (*sh_coeffs_11)[int(7)] * _S1986.differential_0;
    float3  _S2012 = _S1971 * _S1986.differential_0;
    float3  _S2013 = (*sh_coeffs_11)[int(6)] * _S1986.differential_0;
    float3  _S2014 = _S1973 * _S1986.differential_0;
    float3  _S2015 = (*sh_coeffs_11)[int(5)] * _S1986.differential_0;
    float3  _S2016 = _S1975 * _S1986.differential_0;
    float3  _S2017 = (*sh_coeffs_11)[int(4)] * _S1986.differential_0;
    float _S2018 = _S2015.x + _S2015.y + _S2015.z;
    float _S2019 = _S2011.x + _S2011.y + _S2011.z;
    float _S2020 = fTmp1B_11 * _S2002 + x_37 * s_diff_fS2_T_1 + y_14 * _S2007 + 0.54627424478530884f * (_S2017.x + _S2017.y + _S2017.z);
    float _S2021 = fTmp1B_11 * _S2003 + y_14 * s_diff_fS2_T_1 + x_37 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S2009.x + _S2009.y + _S2009.z);
    float _S2022 = y_14 * - _S2021;
    float _S2023 = x_37 * _S2021;
    float _S2024 = z_11 * (1.86588168144226074f * (z_11 * _S2006) + -2.28522896766662598f * (y_14 * _S2004 + x_37 * _S2005) + 0.94617468118667603f * (_S2013.x + _S2013.y + _S2013.z));
    float3  _S2025 = make_float3 (0.48860251903533936f) * _S1986.differential_0;
    float3  _S2026 = - _S2025;
    float3  _S2027 = _S1965 * _S2026;
    float3  _S2028 = (*sh_coeffs_11)[int(3)] * _S2026;
    float3  _S2029 = _S1967 * _S2025;
    float3  _S2030 = (*sh_coeffs_11)[int(2)] * _S2025;
    float3  _S2031 = _S1969 * _S2025;
    float3  _S2032 = (*sh_coeffs_11)[int(1)] * _S2025;
    float _S2033 = (_S1976 * _S2006 + 1.44530570507049561f * (fS1_11 * _S2002 + fC1_11 * _S2003) + -1.09254848957061768f * (y_14 * _S2018 + x_37 * _S2019) + _S2024 + _S2024 + _S2030.x + _S2030.y + _S2030.z) / _S1966;
    float _S2034 = _S1964 * _S2033;
    float _S2035 = (fTmp0C_11 * _S2004 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S2007 + fTmp0B_11 * _S2018 + _S1970 * _S2020 + _S2022 + _S2022 + - (_S2032.x + _S2032.y + _S2032.z)) / _S1966;
    float _S2036 = _S1964 * _S2035;
    float _S2037 = (fTmp0C_11 * _S2005 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S2019 + 2.0f * (y_14 * _S2020) + _S2023 + _S2023 + _S2028.x + _S2028.y + _S2028.z) / _S1966;
    float _S2038 = _S1964 * _S2037;
    float _S2039 = _S1962 * - _S2033 + _S1961 * - _S2035 + _S1960 * - _S2037;
    DiffPair_float_0 _S2040;
    (&_S2040)->primal_0 = _S1963;
    (&_S2040)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2040, _S2039);
    float _S2041 = _S1962 * _S2040.differential_0;
    float _S2042 = _S1961 * _S2040.differential_0;
    float _S2043 = _S1960 * _S2040.differential_0;
    float3  _S2044 = make_float3 (0.282094806432724f) * _S1986.differential_0;
    float3  _S2045 = make_float3 (_S2038 + _S2043 + _S2043, _S2036 + _S2042 + _S2042, _S2034 + _S2041 + _S2041);
    float3  _S2046 = - - _S2045;
    Matrix<float, 3, 3>  _S2047 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2048;
    (&_S2048)->primal_0 = _S1929;
    (&_S2048)->differential_0 = _S2047;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2049;
    (&_S2049)->primal_0 = t_14;
    (&_S2049)->differential_0 = _S1985;
    s_bwd_prop_mul_1(&_S2048, &_S2049, _S2046);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2050 = _S2048;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2051 = _S2049;
    float2  _S2052 = _S1932;
    *&((&_S2052)->y) = v_conic_1.z;
    float2  _S2053 = _S1932;
    *&((&_S2053)->y) = v_conic_1.y;
    *&((&_S2053)->x) = v_conic_1.x;
    float _S2054 = 0.5f * v_depth_1;
    DiffPair_float_0 _S2055;
    (&_S2055)->primal_0 = _S1958;
    (&_S2055)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2055, _S2054);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2056;
    (&_S2056)->primal_0 = mean_c_11;
    (&_S2056)->differential_0 = _S1985;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2057;
    (&_S2057)->primal_0 = mean_c_11;
    (&_S2057)->differential_0 = _S1985;
    s_bwd_prop_dot_0(&_S2056, &_S2057, _S2055.differential_0);
    DiffPair_float_0 _S2058;
    (&_S2058)->primal_0 = _S1957;
    (&_S2058)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2058, 0.0f);
    DiffPair_float_0 _S2059;
    (&_S2059)->primal_0 = _S1956;
    (&_S2059)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2059, 0.0f);
    DiffPair_float_0 _S2060;
    (&_S2060)->primal_0 = 3.32999992370605469f;
    (&_S2060)->differential_0 = 0.0f;
    DiffPair_float_0 _S2061;
    (&_S2061)->primal_0 = _S1955;
    (&_S2061)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2060, &_S2061, 0.0f);
    DiffPair_float_0 _S2062;
    (&_S2062)->primal_0 = _S1954;
    (&_S2062)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2062, _S2061.differential_0);
    float _S2063 = 2.0f * _S2062.differential_0;
    DiffPair_float_0 _S2064;
    (&_S2064)->primal_0 = _S1953;
    (&_S2064)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2064, _S2063);
    float2  _S2065 = make_float2 (_S2059.differential_0, 0.0f);
    float _S2066 = v_opacity_1 + 254.9999847412109375f * _S2064.differential_0;
    Matrix<float, 2, 2>  _S2067 = _S1917;
    _S2067[int(1)] = _S2052;
    _S2067[int(0)] = _S2053;
    Matrix<float, 2, 2>  _S2068 = _S2067;
    FixedArray<float3 , 16>  _S2069;
    _S2069[int(0)] = _S1985;
    _S2069[int(1)] = _S1985;
    _S2069[int(2)] = _S1985;
    _S2069[int(3)] = _S1985;
    _S2069[int(4)] = _S1985;
    _S2069[int(5)] = _S1985;
    _S2069[int(6)] = _S1985;
    _S2069[int(7)] = _S1985;
    _S2069[int(8)] = _S1985;
    _S2069[int(9)] = _S1985;
    _S2069[int(10)] = _S1985;
    _S2069[int(11)] = _S1985;
    _S2069[int(12)] = _S1985;
    _S2069[int(13)] = _S1985;
    _S2069[int(14)] = _S1985;
    _S2069[int(15)] = _S1985;
    _S2069[int(7)] = _S2010;
    _S2069[int(0)] = _S2044;
    _S2069[int(1)] = _S2031;
    _S2069[int(2)] = _S2029;
    _S2069[int(3)] = _S2027;
    _S2069[int(4)] = _S2016;
    _S2069[int(5)] = _S2014;
    _S2069[int(6)] = _S2012;
    _S2069[int(15)] = _S1988;
    _S2069[int(8)] = _S2008;
    _S2069[int(9)] = _S2000;
    _S2069[int(10)] = _S1998;
    _S2069[int(11)] = _S1996;
    _S2069[int(12)] = _S1994;
    _S2069[int(13)] = _S1992;
    _S2069[int(14)] = _S1990;
    float3  _S2070 = _S2069[int(0)];
    float3  _S2071 = _S2069[int(1)];
    float3  _S2072 = _S2069[int(2)];
    float3  _S2073 = _S2069[int(3)];
    float3  _S2074 = _S2069[int(4)];
    float3  _S2075 = _S2069[int(5)];
    float3  _S2076 = _S2069[int(6)];
    float3  _S2077 = _S2069[int(7)];
    float3  _S2078 = _S2069[int(8)];
    float3  _S2079 = _S2069[int(9)];
    float3  _S2080 = _S2069[int(10)];
    float3  _S2081 = _S2069[int(11)];
    float3  _S2082 = _S2069[int(12)];
    float3  _S2083 = _S2069[int(13)];
    float3  _S2084 = _S2069[int(14)];
    float3  _S2085 = _S2069[int(15)];
    float3  _S2086 = _S2057.differential_0 + _S2056.differential_0;
    float2  _S2087 = make_float2 (0.0f, _S2058.differential_0);
    float _S2088;
    if(antialiased_11)
    {
        float _S2089 = _S1950 * _S2066;
        _S1952 = _S1945 * _S2066;
        _S2088 = _S2089;
    }
    else
    {
        _S1952 = _S2066;
        _S2088 = 0.0f;
    }
    float _S2090 = - (_S1952 / _S1951);
    DiffPair_float_0 _S2091;
    (&_S2091)->primal_0 = _S1948;
    (&_S2091)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2091, _S2090);
    float _S2092 = - _S2091.differential_0;
    float _S2093 = invdet_8 * _S2068.rows[int(1)].y;
    float _S2094 = - (invdet_8 * _S2068.rows[int(1)].x);
    float _S2095 = - (invdet_8 * _S2068.rows[int(0)].y);
    float _S2096 = invdet_8 * _S2068.rows[int(0)].x;
    float _S2097 = - ((_S1937 * _S2068.rows[int(1)].y + _S1947 * _S2068.rows[int(1)].x + _S1946 * _S2068.rows[int(0)].y + _S1939 * _S2068.rows[int(0)].x) / _S1943);
    DiffPair_float_0 _S2098;
    (&_S2098)->primal_0 = _S1944;
    (&_S2098)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2098, _S2088);
    DiffPair_float_0 _S2099;
    (&_S2099)->primal_0 = 0.0f;
    (&_S2099)->differential_0 = 0.0f;
    DiffPair_float_0 _S2100;
    (&_S2100)->primal_0 = _S1942;
    (&_S2100)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2099, &_S2100, _S2098.differential_0);
    float _S2101 = _S2100.differential_0 / _S1943;
    float s_diff_det_orig_T_1 = det_blur_7 * _S2101;
    float _S2102 = _S2097 + det_orig_12 * - _S2101;
    float _S2103 = - _S2102;
    float _S2104 = _S1937 * _S2102;
    float _S2105 = _S1939 * _S2102;
    Matrix<float, 2, 2>  _S2106 = _S1917;
    _S2106[int(1)] = _S2087;
    _S2106[int(0)] = _S2065;
    _S1938 = _S2106;
    *&(((&_S1938)->rows + (int(1)))->y) = 0.0f;
    float _S2107 = _S2096 + _S2104 + _S2106.rows[int(1)].y;
    *&(((&_S1938)->rows + (int(0)))->x) = 0.0f;
    float _S2108 = _S2093 + _S2105 + _S2106.rows[int(0)].x;
    float _S2109 = _S2103 + - s_diff_det_orig_T_1;
    float _S2110 = _S2094 + _S1935._S1232.rows[int(0)].y * _S2109;
    float _S2111 = _S2095 + _S1935._S1232.rows[int(1)].x * _S2109;
    float _S2112 = _S1935._S1232.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S2113 = _S2107 + _S1935._S1232.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S2114 = _S1932;
    *&((&_S2114)->x) = _S2110;
    *&((&_S2114)->y) = _S2113;
    float _S2115 = _S2108 + _S2112;
    float2  _S2116 = _S1932;
    *&((&_S2116)->y) = _S2111;
    *&((&_S2116)->x) = _S2115;
    Matrix<float, 2, 2>  _S2117 = _S1917;
    _S2117[int(1)] = _S2114;
    _S2117[int(0)] = _S2116;
    Matrix<float, 2, 2>  _S2118 = _S1938 + _S2117;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2119;
    (&_S2119)->primal_0 = mean_c_11;
    (&_S2119)->differential_0 = _S1985;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2120;
    (&_S2120)->primal_0 = _S1930;
    (&_S2120)->differential_0 = _S2047;
    float4  _S2121 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2122;
    (&_S2122)->primal_0 = intrins_11;
    (&_S2122)->differential_0 = _S2121;
    FixedArray<float, 10>  _S2123 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S2124;
    (&_S2124)->primal_0 = *dist_coeffs_27;
    (&_S2124)->differential_0 = _S2123;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2125 = _S1935._S1233;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2119, &_S2120, &_S2122, &_S2124, _S2118, v_mean2d_1, &_S2125);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2126;
    (&_S2126)->primal_0 = _S1928;
    (&_S2126)->differential_0 = _S2047;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2127;
    (&_S2127)->primal_0 = _S1929;
    (&_S2127)->differential_0 = _S2047;
    s_bwd_prop_mul_4(&_S2126, &_S2127, _S2120.differential_0);
    Matrix<float, 3, 3>  _S2128 = transpose_0(_S2127.differential_0 + _S2050.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2129;
    (&_S2129)->primal_0 = R_15;
    (&_S2129)->differential_0 = _S2047;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2130;
    (&_S2130)->primal_0 = _S1927;
    (&_S2130)->differential_0 = _S2047;
    s_bwd_prop_mul_4(&_S2129, &_S2130, _S2126.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2131;
    (&_S2131)->primal_0 = _S1925;
    (&_S2131)->differential_0 = _S2047;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2132;
    (&_S2132)->primal_0 = _S1926;
    (&_S2132)->differential_0 = _S2047;
    s_bwd_prop_mul_4(&_S2131, &_S2132, _S2130.differential_0);
    Matrix<float, 3, 3>  _S2133 = _S2131.differential_0 + transpose_0(_S2132.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2134;
    (&_S2134)->primal_0 = _S1924;
    (&_S2134)->differential_0 = _S2047;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2135;
    (&_S2135)->primal_0 = S_1;
    (&_S2135)->differential_0 = _S2047;
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
    float _S2152 = quat_15.w * (_S2141 + _S2145);
    float _S2153 = quat_15.z * (_S2137 + _S2145);
    float _S2154 = quat_15.y * (_S2137 + _S2141);
    float _S2155 = quat_15.x * _S2146 + quat_15.z * _S2149 + quat_15.y * _S2150 + _S2152 + _S2152;
    float _S2156 = quat_15.x * _S2147 + quat_15.w * _S2149 + quat_15.y * _S2151 + _S2153 + _S2153;
    float _S2157 = quat_15.x * _S2148 + quat_15.w * _S2150 + quat_15.z * _S2151 + _S2154 + _S2154;
    float _S2158 = quat_15.w * _S2146 + quat_15.z * _S2147 + quat_15.y * _S2148;
    float3  _S2159 = _S1985;
    *&((&_S2159)->z) = _S2135.differential_0.rows[int(2)].z;
    *&((&_S2159)->y) = _S2135.differential_0.rows[int(1)].y;
    *&((&_S2159)->x) = _S2135.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2160;
    (&_S2160)->primal_0 = scale_14;
    (&_S2160)->differential_0 = _S1985;
    s_bwd_prop_exp_1(&_S2160, _S2159);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2161 = _S2160;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2162;
    (&_S2162)->primal_0 = mean_c_11;
    (&_S2162)->differential_0 = _S1985;
    s_bwd_length_impl_0(&_S2162, 0.0f);
    float3  _S2163 = _S2119.differential_0 + _S2162.differential_0 + _S2086;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2164;
    (&_S2164)->primal_0 = R_15;
    (&_S2164)->differential_0 = _S2047;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2165;
    (&_S2165)->primal_0 = mean_12;
    (&_S2165)->differential_0 = _S1985;
    s_bwd_prop_mul_1(&_S2164, &_S2165, _S2163);
    float3  _S2166 = _S2163 + _S2051.differential_0;
    Matrix<float, 3, 3>  _S2167 = _S2128 + _S2129.differential_0 + _S2164.differential_0;
    float4  _S2168 = _S2121;
    *&((&_S2168)->w) = _S2155;
    *&((&_S2168)->z) = _S2156;
    *&((&_S2168)->y) = _S2157;
    *&((&_S2168)->x) = _S2158;
    float4  _S2169 = _S2168;
    float3  _S2170 = _S2165.differential_0 + _S2045;
    *v_mean_1 = _S2170;
    *v_quat_1 = _S2169;
    *v_scale_1 = _S2161.differential_0;
    *v_in_opacity_1 = _S2092;
    (*v_sh_coeffs_1)[int(0)] = _S2070;
    (*v_sh_coeffs_1)[int(1)] = _S2071;
    (*v_sh_coeffs_1)[int(2)] = _S2072;
    (*v_sh_coeffs_1)[int(3)] = _S2073;
    (*v_sh_coeffs_1)[int(4)] = _S2074;
    (*v_sh_coeffs_1)[int(5)] = _S2075;
    (*v_sh_coeffs_1)[int(6)] = _S2076;
    (*v_sh_coeffs_1)[int(7)] = _S2077;
    (*v_sh_coeffs_1)[int(8)] = _S2078;
    (*v_sh_coeffs_1)[int(9)] = _S2079;
    (*v_sh_coeffs_1)[int(10)] = _S2080;
    (*v_sh_coeffs_1)[int(11)] = _S2081;
    (*v_sh_coeffs_1)[int(12)] = _S2082;
    (*v_sh_coeffs_1)[int(13)] = _S2083;
    (*v_sh_coeffs_1)[int(14)] = _S2084;
    (*v_sh_coeffs_1)[int(15)] = _S2085;
    *v_R_2 = _S2167;
    *v_t_2 = _S2166;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_16, float3  t_15, float fx_21, float fy_21, float cx_16, float cy_16, FixedArray<float, 10>  * dist_coeffs_28, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_1(R_16, mean_13) + t_15;
    float3  _S2171 = s_primal_ctx_exp_0(scale_15);
    float _S2172 = quat_16.y;
    float x2_16 = _S2172 * _S2172;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_28 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S2173 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_28), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_28), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S2171.x, 0.0f, 0.0f, 0.0f, _S2171.y, 0.0f, 0.0f, 0.0f, _S2171.z);
    Matrix<float, 3, 3>  _S2174 = s_primal_ctx_mul_2(_S2173, S_2);
    Matrix<float, 3, 3>  _S2175 = transpose_0(_S2174);
    Matrix<float, 3, 3>  _S2176 = s_primal_ctx_mul_2(_S2174, _S2175);
    Matrix<float, 3, 3>  _S2177 = s_primal_ctx_mul_2(R_16, _S2176);
    Matrix<float, 3, 3>  _S2178 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S2179 = s_primal_ctx_mul_2(_S2177, _S2178);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_21, 0.0f, 0.0f, 0.0f, fy_21, 0.0f);
    Matrix<float, 2, 3>  _S2180 = s_primal_ctx_mul_3(J_11, _S2179);
    Matrix<float, 3, 2>  _S2181 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S2182 = s_primal_ctx_mul_4(_S2180, _S2181);
    float _S2183 = _S2182.rows[int(0)].y * _S2182.rows[int(1)].x;
    float det_orig_13 = _S2182.rows[int(0)].x * _S2182.rows[int(1)].y - _S2183;
    float _S2184 = _S2182.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2185 = _S2182;
    *&(((&_S2185)->rows + (int(0)))->x) = _S2184;
    float _S2186 = _S2182.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2185)->rows + (int(1)))->y) = _S2186;
    Matrix<float, 2, 2>  _S2187 = _S2185;
    Matrix<float, 2, 2>  _S2188 = _S2185;
    float det_blur_8 = _S2184 * _S2186 - _S2183;
    float _S2189 = det_orig_13 / det_blur_8;
    float _S2190 = det_blur_8 * det_blur_8;
    float _S2191 = s_primal_ctx_max_0(0.0f, _S2189);
    float _S2192 = s_primal_ctx_sqrt_0(_S2191);
    float invdet_9 = 1.0f / det_blur_8;
    float _S2193 = - _S2182.rows[int(0)].y;
    float _S2194 = - _S2182.rows[int(1)].x;
    float _S2195 = - in_opacity_12;
    float _S2196 = 1.0f + s_primal_ctx_exp_1(_S2195);
    float _S2197 = 1.0f / _S2196;
    float _S2198 = _S2196 * _S2196;
    float _S2199;
    if(antialiased_12)
    {
        _S2199 = _S2197 * _S2192;
    }
    else
    {
        _S2199 = _S2197;
    }
    float _S2200 = _S2199 / 0.00392156885936856f;
    float _S2201 = 2.0f * s_primal_ctx_log_0(_S2200);
    float _S2202 = s_primal_ctx_sqrt_0(_S2201);
    float _S2203 = _S2187.rows[int(0)].x;
    float _S2204 = _S2188.rows[int(1)].y;
    float _S2205 = s_primal_ctx_dot_0(mean_c_12, mean_c_12) + 9.99999997475242708e-07f;
    float3  _S2206 = mean_13 - - s_primal_ctx_mul_1(_S2178, t_15);
    float _S2207 = _S2206.x;
    float _S2208 = _S2206.y;
    float _S2209 = _S2206.z;
    float _S2210 = _S2207 * _S2207 + _S2208 * _S2208 + _S2209 * _S2209;
    float _S2211 = s_primal_ctx_sqrt_0(_S2210);
    float x_38 = _S2207 / _S2211;
    float3  _S2212 = make_float3 (x_38);
    float _S2213 = _S2211 * _S2211;
    float y_15 = _S2208 / _S2211;
    float z_12 = _S2209 / _S2211;
    float3  _S2214 = make_float3 (z_12);
    float _S2215 = - y_15;
    float3  _S2216 = make_float3 (_S2215);
    float z2_29 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_38 * x_38 - y_15 * y_15;
    float _S2217 = 2.0f * x_38;
    float fS1_12 = _S2217 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_29 - 0.31539157032966614f;
    float3  _S2218 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_38;
    float3  _S2219 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S2220 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S2221 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S2222 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_29 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S2223 = 1.86588168144226074f * z2_29 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S2223;
    float3  _S2224 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_38;
    float3  _S2225 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S2226 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S2227 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S2228 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_38 * fC1_12 - y_15 * fS1_12);
    float3  _S2229 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_38 * fS1_12 + y_15 * fC1_12);
    float3  _S2230 = make_float3 (pSH9_2);
    float3  _S2231 = make_float3 (0.0f);
    float3  _S2232 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2233;
    (&_S2233)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2215) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_38) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S2233)->differential_0 = _S2232;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2234;
    (&_S2234)->primal_0 = _S2231;
    (&_S2234)->differential_0 = _S2232;
    s_bwd_prop_max_0(&_S2233, &_S2234, v_rgb_2);
    float3  _S2235 = _S2229 * _S2233.differential_0;
    float3  _S2236 = (*sh_coeffs_12)[int(15)] * _S2233.differential_0;
    float3  _S2237 = _S2227 * _S2233.differential_0;
    float3  _S2238 = (*sh_coeffs_12)[int(14)] * _S2233.differential_0;
    float3  _S2239 = _S2225 * _S2233.differential_0;
    float3  _S2240 = (*sh_coeffs_12)[int(13)] * _S2233.differential_0;
    float3  _S2241 = _S2224 * _S2233.differential_0;
    float3  _S2242 = (*sh_coeffs_12)[int(12)] * _S2233.differential_0;
    float3  _S2243 = _S2226 * _S2233.differential_0;
    float3  _S2244 = (*sh_coeffs_12)[int(11)] * _S2233.differential_0;
    float3  _S2245 = _S2228 * _S2233.differential_0;
    float3  _S2246 = (*sh_coeffs_12)[int(10)] * _S2233.differential_0;
    float3  _S2247 = _S2230 * _S2233.differential_0;
    float3  _S2248 = (*sh_coeffs_12)[int(9)] * _S2233.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S2248.x + _S2248.y + _S2248.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S2236.x + _S2236.y + _S2236.z);
    float _S2249 = _S2246.x + _S2246.y + _S2246.z;
    float _S2250 = _S2238.x + _S2238.y + _S2238.z;
    float _S2251 = _S2244.x + _S2244.y + _S2244.z;
    float _S2252 = _S2240.x + _S2240.y + _S2240.z;
    float _S2253 = _S2242.x + _S2242.y + _S2242.z;
    float _S2254 = - s_diff_fC2_T_2;
    float3  _S2255 = _S2221 * _S2233.differential_0;
    float3  _S2256 = (*sh_coeffs_12)[int(8)] * _S2233.differential_0;
    float3  _S2257 = _S2219 * _S2233.differential_0;
    float3  _S2258 = (*sh_coeffs_12)[int(7)] * _S2233.differential_0;
    float3  _S2259 = _S2218 * _S2233.differential_0;
    float3  _S2260 = (*sh_coeffs_12)[int(6)] * _S2233.differential_0;
    float3  _S2261 = _S2220 * _S2233.differential_0;
    float3  _S2262 = (*sh_coeffs_12)[int(5)] * _S2233.differential_0;
    float3  _S2263 = _S2222 * _S2233.differential_0;
    float3  _S2264 = (*sh_coeffs_12)[int(4)] * _S2233.differential_0;
    float _S2265 = _S2262.x + _S2262.y + _S2262.z;
    float _S2266 = _S2258.x + _S2258.y + _S2258.z;
    float _S2267 = fTmp1B_12 * _S2249 + x_38 * s_diff_fS2_T_2 + y_15 * _S2254 + 0.54627424478530884f * (_S2264.x + _S2264.y + _S2264.z);
    float _S2268 = fTmp1B_12 * _S2250 + y_15 * s_diff_fS2_T_2 + x_38 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S2256.x + _S2256.y + _S2256.z);
    float _S2269 = y_15 * - _S2268;
    float _S2270 = x_38 * _S2268;
    float _S2271 = z_12 * (1.86588168144226074f * (z_12 * _S2253) + -2.28522896766662598f * (y_15 * _S2251 + x_38 * _S2252) + 0.94617468118667603f * (_S2260.x + _S2260.y + _S2260.z));
    float3  _S2272 = make_float3 (0.48860251903533936f) * _S2233.differential_0;
    float3  _S2273 = - _S2272;
    float3  _S2274 = _S2212 * _S2273;
    float3  _S2275 = (*sh_coeffs_12)[int(3)] * _S2273;
    float3  _S2276 = _S2214 * _S2272;
    float3  _S2277 = (*sh_coeffs_12)[int(2)] * _S2272;
    float3  _S2278 = _S2216 * _S2272;
    float3  _S2279 = (*sh_coeffs_12)[int(1)] * _S2272;
    float _S2280 = (_S2223 * _S2253 + 1.44530570507049561f * (fS1_12 * _S2249 + fC1_12 * _S2250) + -1.09254848957061768f * (y_15 * _S2265 + x_38 * _S2266) + _S2271 + _S2271 + _S2277.x + _S2277.y + _S2277.z) / _S2213;
    float _S2281 = _S2211 * _S2280;
    float _S2282 = (fTmp0C_12 * _S2251 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S2254 + fTmp0B_12 * _S2265 + _S2217 * _S2267 + _S2269 + _S2269 + - (_S2279.x + _S2279.y + _S2279.z)) / _S2213;
    float _S2283 = _S2211 * _S2282;
    float _S2284 = (fTmp0C_12 * _S2252 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S2266 + 2.0f * (y_15 * _S2267) + _S2270 + _S2270 + _S2275.x + _S2275.y + _S2275.z) / _S2213;
    float _S2285 = _S2211 * _S2284;
    float _S2286 = _S2209 * - _S2280 + _S2208 * - _S2282 + _S2207 * - _S2284;
    DiffPair_float_0 _S2287;
    (&_S2287)->primal_0 = _S2210;
    (&_S2287)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2287, _S2286);
    float _S2288 = _S2209 * _S2287.differential_0;
    float _S2289 = _S2208 * _S2287.differential_0;
    float _S2290 = _S2207 * _S2287.differential_0;
    float3  _S2291 = make_float3 (0.282094806432724f) * _S2233.differential_0;
    float3  _S2292 = make_float3 (_S2285 + _S2290 + _S2290, _S2283 + _S2289 + _S2289, _S2281 + _S2288 + _S2288);
    float3  _S2293 = - - _S2292;
    Matrix<float, 3, 3>  _S2294 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2295;
    (&_S2295)->primal_0 = _S2178;
    (&_S2295)->differential_0 = _S2294;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2296;
    (&_S2296)->primal_0 = t_15;
    (&_S2296)->differential_0 = _S2232;
    s_bwd_prop_mul_1(&_S2295, &_S2296, _S2293);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2297 = _S2295;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2298 = _S2296;
    float2  _S2299 = make_float2 (0.0f);
    float2  _S2300 = _S2299;
    *&((&_S2300)->y) = v_conic_2.z;
    float2  _S2301 = _S2299;
    *&((&_S2301)->y) = v_conic_2.y;
    *&((&_S2301)->x) = v_conic_2.x;
    float _S2302 = 0.5f * v_depth_2;
    DiffPair_float_0 _S2303;
    (&_S2303)->primal_0 = _S2205;
    (&_S2303)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2303, _S2302);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2304;
    (&_S2304)->primal_0 = mean_c_12;
    (&_S2304)->differential_0 = _S2232;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2305;
    (&_S2305)->primal_0 = mean_c_12;
    (&_S2305)->differential_0 = _S2232;
    s_bwd_prop_dot_0(&_S2304, &_S2305, _S2303.differential_0);
    DiffPair_float_0 _S2306;
    (&_S2306)->primal_0 = _S2204;
    (&_S2306)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2306, 0.0f);
    DiffPair_float_0 _S2307;
    (&_S2307)->primal_0 = _S2203;
    (&_S2307)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2307, 0.0f);
    DiffPair_float_0 _S2308;
    (&_S2308)->primal_0 = 3.32999992370605469f;
    (&_S2308)->differential_0 = 0.0f;
    DiffPair_float_0 _S2309;
    (&_S2309)->primal_0 = _S2202;
    (&_S2309)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2308, &_S2309, 0.0f);
    DiffPair_float_0 _S2310;
    (&_S2310)->primal_0 = _S2201;
    (&_S2310)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2310, _S2309.differential_0);
    float _S2311 = 2.0f * _S2310.differential_0;
    DiffPair_float_0 _S2312;
    (&_S2312)->primal_0 = _S2200;
    (&_S2312)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2312, _S2311);
    float _S2313 = v_opacity_2 + 254.9999847412109375f * _S2312.differential_0;
    Matrix<float, 2, 2>  _S2314 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2315 = _S2314;
    _S2315[int(1)] = _S2300;
    _S2315[int(0)] = _S2301;
    Matrix<float, 2, 2>  _S2316 = _S2315;
    FixedArray<float3 , 16>  _S2317;
    _S2317[int(0)] = _S2232;
    _S2317[int(1)] = _S2232;
    _S2317[int(2)] = _S2232;
    _S2317[int(3)] = _S2232;
    _S2317[int(4)] = _S2232;
    _S2317[int(5)] = _S2232;
    _S2317[int(6)] = _S2232;
    _S2317[int(7)] = _S2232;
    _S2317[int(8)] = _S2232;
    _S2317[int(9)] = _S2232;
    _S2317[int(10)] = _S2232;
    _S2317[int(11)] = _S2232;
    _S2317[int(12)] = _S2232;
    _S2317[int(13)] = _S2232;
    _S2317[int(14)] = _S2232;
    _S2317[int(15)] = _S2232;
    _S2317[int(7)] = _S2257;
    _S2317[int(0)] = _S2291;
    _S2317[int(1)] = _S2278;
    _S2317[int(2)] = _S2276;
    _S2317[int(3)] = _S2274;
    _S2317[int(4)] = _S2263;
    _S2317[int(5)] = _S2261;
    _S2317[int(6)] = _S2259;
    _S2317[int(15)] = _S2235;
    _S2317[int(8)] = _S2255;
    _S2317[int(9)] = _S2247;
    _S2317[int(10)] = _S2245;
    _S2317[int(11)] = _S2243;
    _S2317[int(12)] = _S2241;
    _S2317[int(13)] = _S2239;
    _S2317[int(14)] = _S2237;
    float3  _S2318 = _S2317[int(0)];
    float3  _S2319 = _S2317[int(1)];
    float3  _S2320 = _S2317[int(2)];
    float3  _S2321 = _S2317[int(3)];
    float3  _S2322 = _S2317[int(4)];
    float3  _S2323 = _S2317[int(5)];
    float3  _S2324 = _S2317[int(6)];
    float3  _S2325 = _S2317[int(7)];
    float3  _S2326 = _S2317[int(8)];
    float3  _S2327 = _S2317[int(9)];
    float3  _S2328 = _S2317[int(10)];
    float3  _S2329 = _S2317[int(11)];
    float3  _S2330 = _S2317[int(12)];
    float3  _S2331 = _S2317[int(13)];
    float3  _S2332 = _S2317[int(14)];
    float3  _S2333 = _S2317[int(15)];
    float3  _S2334 = _S2305.differential_0 + _S2304.differential_0;
    float2  _S2335 = make_float2 (0.0f, _S2306.differential_0);
    float2  _S2336 = make_float2 (_S2307.differential_0, 0.0f);
    float _S2337;
    if(antialiased_12)
    {
        float _S2338 = _S2197 * _S2313;
        _S2199 = _S2192 * _S2313;
        _S2337 = _S2338;
    }
    else
    {
        _S2199 = _S2313;
        _S2337 = 0.0f;
    }
    float _S2339 = - (_S2199 / _S2198);
    DiffPair_float_0 _S2340;
    (&_S2340)->primal_0 = _S2195;
    (&_S2340)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2340, _S2339);
    float _S2341 = - _S2340.differential_0;
    float _S2342 = invdet_9 * _S2316.rows[int(1)].y;
    float _S2343 = - (invdet_9 * _S2316.rows[int(1)].x);
    float _S2344 = - (invdet_9 * _S2316.rows[int(0)].y);
    float _S2345 = invdet_9 * _S2316.rows[int(0)].x;
    float _S2346 = - ((_S2184 * _S2316.rows[int(1)].y + _S2194 * _S2316.rows[int(1)].x + _S2193 * _S2316.rows[int(0)].y + _S2186 * _S2316.rows[int(0)].x) / _S2190);
    DiffPair_float_0 _S2347;
    (&_S2347)->primal_0 = _S2191;
    (&_S2347)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2347, _S2337);
    DiffPair_float_0 _S2348;
    (&_S2348)->primal_0 = 0.0f;
    (&_S2348)->differential_0 = 0.0f;
    DiffPair_float_0 _S2349;
    (&_S2349)->primal_0 = _S2189;
    (&_S2349)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2348, &_S2349, _S2347.differential_0);
    float _S2350 = _S2349.differential_0 / _S2190;
    float s_diff_det_orig_T_2 = det_blur_8 * _S2350;
    float _S2351 = _S2346 + det_orig_13 * - _S2350;
    float _S2352 = - _S2351;
    float _S2353 = _S2184 * _S2351;
    float _S2354 = _S2186 * _S2351;
    Matrix<float, 2, 2>  _S2355 = _S2314;
    _S2355[int(1)] = _S2335;
    _S2355[int(0)] = _S2336;
    _S2185 = _S2355;
    *&(((&_S2185)->rows + (int(1)))->y) = 0.0f;
    float _S2356 = _S2345 + _S2353 + _S2355.rows[int(1)].y;
    *&(((&_S2185)->rows + (int(0)))->x) = 0.0f;
    float _S2357 = _S2342 + _S2354 + _S2355.rows[int(0)].x;
    float _S2358 = _S2352 + - s_diff_det_orig_T_2;
    float _S2359 = _S2343 + _S2182.rows[int(0)].y * _S2358;
    float _S2360 = _S2344 + _S2182.rows[int(1)].x * _S2358;
    float _S2361 = _S2182.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S2362 = _S2356 + _S2182.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S2363 = _S2299;
    *&((&_S2363)->x) = _S2359;
    *&((&_S2363)->y) = _S2362;
    float _S2364 = _S2357 + _S2361;
    float2  _S2365 = _S2299;
    *&((&_S2365)->y) = _S2360;
    *&((&_S2365)->x) = _S2364;
    float _S2366 = fy_21 * v_mean2d_2.y;
    float _S2367 = fx_21 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S2368 = _S2314;
    _S2368[int(1)] = _S2363;
    _S2368[int(0)] = _S2365;
    Matrix<float, 2, 2>  _S2369 = _S2185 + _S2368;
    Matrix<float, 2, 3>  _S2370 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2371;
    (&_S2371)->primal_0 = _S2180;
    (&_S2371)->differential_0 = _S2370;
    Matrix<float, 3, 2>  _S2372 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2373;
    (&_S2373)->primal_0 = _S2181;
    (&_S2373)->differential_0 = _S2372;
    s_bwd_prop_mul_2(&_S2371, &_S2373, _S2369);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2374;
    (&_S2374)->primal_0 = J_11;
    (&_S2374)->differential_0 = _S2370;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2375;
    (&_S2375)->primal_0 = _S2179;
    (&_S2375)->differential_0 = _S2294;
    s_bwd_prop_mul_3(&_S2374, &_S2375, _S2371.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2376;
    (&_S2376)->primal_0 = _S2177;
    (&_S2376)->differential_0 = _S2294;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2377;
    (&_S2377)->primal_0 = _S2178;
    (&_S2377)->differential_0 = _S2294;
    s_bwd_prop_mul_4(&_S2376, &_S2377, _S2375.differential_0);
    Matrix<float, 3, 3>  _S2378 = transpose_0(_S2377.differential_0 + _S2297.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2379;
    (&_S2379)->primal_0 = R_16;
    (&_S2379)->differential_0 = _S2294;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2380;
    (&_S2380)->primal_0 = _S2176;
    (&_S2380)->differential_0 = _S2294;
    s_bwd_prop_mul_4(&_S2379, &_S2380, _S2376.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2381;
    (&_S2381)->primal_0 = _S2174;
    (&_S2381)->differential_0 = _S2294;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2382;
    (&_S2382)->primal_0 = _S2175;
    (&_S2382)->differential_0 = _S2294;
    s_bwd_prop_mul_4(&_S2381, &_S2382, _S2380.differential_0);
    Matrix<float, 3, 3>  _S2383 = _S2381.differential_0 + transpose_0(_S2382.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2384;
    (&_S2384)->primal_0 = _S2173;
    (&_S2384)->differential_0 = _S2294;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2385;
    (&_S2385)->primal_0 = S_2;
    (&_S2385)->differential_0 = _S2294;
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
    float _S2402 = quat_16.w * (_S2391 + _S2395);
    float _S2403 = quat_16.z * (_S2387 + _S2395);
    float _S2404 = quat_16.y * (_S2387 + _S2391);
    float _S2405 = quat_16.x * _S2396 + quat_16.z * _S2399 + quat_16.y * _S2400 + _S2402 + _S2402;
    float _S2406 = quat_16.x * _S2397 + quat_16.w * _S2399 + quat_16.y * _S2401 + _S2403 + _S2403;
    float _S2407 = quat_16.x * _S2398 + quat_16.w * _S2400 + quat_16.z * _S2401 + _S2404 + _S2404;
    float _S2408 = quat_16.w * _S2396 + quat_16.z * _S2397 + quat_16.y * _S2398;
    float3  _S2409 = _S2232;
    *&((&_S2409)->z) = _S2385.differential_0.rows[int(2)].z;
    *&((&_S2409)->y) = _S2385.differential_0.rows[int(1)].y;
    *&((&_S2409)->x) = _S2385.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2410;
    (&_S2410)->primal_0 = scale_15;
    (&_S2410)->differential_0 = _S2232;
    s_bwd_prop_exp_1(&_S2410, _S2409);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2411 = _S2410;
    float3  _S2412 = _S2232;
    *&((&_S2412)->y) = _S2366;
    *&((&_S2412)->x) = _S2367;
    float3  _S2413 = _S2334 + _S2412;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2414;
    (&_S2414)->primal_0 = R_16;
    (&_S2414)->differential_0 = _S2294;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2415;
    (&_S2415)->primal_0 = mean_13;
    (&_S2415)->differential_0 = _S2232;
    s_bwd_prop_mul_1(&_S2414, &_S2415, _S2413);
    float3  _S2416 = _S2413 + _S2298.differential_0;
    Matrix<float, 3, 3>  _S2417 = _S2378 + _S2379.differential_0 + _S2414.differential_0;
    float4  _S2418 = make_float4 (0.0f);
    *&((&_S2418)->w) = _S2405;
    *&((&_S2418)->z) = _S2406;
    *&((&_S2418)->y) = _S2407;
    *&((&_S2418)->x) = _S2408;
    float4  _S2419 = _S2418;
    float3  _S2420 = _S2415.differential_0 + _S2292;
    *v_mean_2 = _S2420;
    *v_quat_2 = _S2419;
    *v_scale_2 = _S2411.differential_0;
    *v_in_opacity_2 = _S2341;
    (*v_sh_coeffs_2)[int(0)] = _S2318;
    (*v_sh_coeffs_2)[int(1)] = _S2319;
    (*v_sh_coeffs_2)[int(2)] = _S2320;
    (*v_sh_coeffs_2)[int(3)] = _S2321;
    (*v_sh_coeffs_2)[int(4)] = _S2322;
    (*v_sh_coeffs_2)[int(5)] = _S2323;
    (*v_sh_coeffs_2)[int(6)] = _S2324;
    (*v_sh_coeffs_2)[int(7)] = _S2325;
    (*v_sh_coeffs_2)[int(8)] = _S2326;
    (*v_sh_coeffs_2)[int(9)] = _S2327;
    (*v_sh_coeffs_2)[int(10)] = _S2328;
    (*v_sh_coeffs_2)[int(11)] = _S2329;
    (*v_sh_coeffs_2)[int(12)] = _S2330;
    (*v_sh_coeffs_2)[int(13)] = _S2331;
    (*v_sh_coeffs_2)[int(14)] = _S2332;
    (*v_sh_coeffs_2)[int(15)] = _S2333;
    *v_R_3 = _S2417;
    *v_t_3 = _S2416;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2421;
};

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_17, float3  t_16, float fx_22, float fy_22, float cx_17, float cy_17, FixedArray<float, 10>  * dist_coeffs_29, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_4)
{
    Matrix<float, 2, 2>  _S2422 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2423;
    (&_S2423)->_S2421 = _S2422;
    float2  _S2424 = make_float2 (0.0f);
    float3  _S2425 = make_float3 (0.0f);
    float4  intrins_12 = make_float4 (fx_22, fy_22, cx_17, cy_17);
    float3  _S2426 = s_primal_ctx_exp_0(scale_16);
    float _S2427 = quat_17.y;
    float x2_17 = _S2427 * _S2427;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_30 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2428 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_30), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_30), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17))));
    FixedArray<float3 , 7>  _S2429 = {
        _S2425, _S2425, _S2425, _S2425, _S2425, _S2425, _S2425
    };
    FixedArray<float, 7>  _S2430 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2431;
    (&_S2431)->p_0 = _S2429;
    (&_S2431)->w_mean_0 = _S2430;
    (&_S2431)->w_cov_0 = _S2430;
    (&_S2431)->p_0[int(0)] = mean_14;
    SigmaPoints_0 _S2432 = _S2431;
    (&_S2432)->w_mean_0[int(0)] = 0.0f;
    (&_S2432)->w_cov_0[int(0)] = 2.0f;
    float _S2433 = s_primal_ctx_sqrt_0(3.0f);
    float _S2434 = _S2433 * _S2426.x;
    float3  delta_12 = make_float3 (_S2434) * _S2428.rows[0U];
    float3  _S2435 = mean_14 + delta_12;
    float3  _S2436 = mean_14 - delta_12;
    float _S2437 = _S2433 * _S2426.y;
    float3  delta_13 = make_float3 (_S2437) * _S2428.rows[1U];
    float3  _S2438 = mean_14 + delta_13;
    float3  _S2439 = mean_14 - delta_13;
    float _S2440 = _S2433 * _S2426.z;
    float3  delta_14 = make_float3 (_S2440) * _S2428.rows[2U];
    float3  _S2441 = mean_14 + delta_14;
    float3  _S2442 = mean_14 - delta_14;
    (&_S2432)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2432)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2432)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2432)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2432)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2432)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2432)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2432)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2432)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2432)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2432)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2432)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2443 = _S2431;
    (&_S2432)->p_0[0U] = s_primal_ctx_mul_1(R_17, _S2431.p_0[0U]) + t_16;
    (&_S2432)->p_0[1U] = s_primal_ctx_mul_1(R_17, _S2435) + t_16;
    (&_S2432)->p_0[2U] = s_primal_ctx_mul_1(R_17, _S2438) + t_16;
    (&_S2432)->p_0[3U] = s_primal_ctx_mul_1(R_17, _S2441) + t_16;
    (&_S2432)->p_0[4U] = s_primal_ctx_mul_1(R_17, _S2436) + t_16;
    (&_S2432)->p_0[5U] = s_primal_ctx_mul_1(R_17, _S2439) + t_16;
    (&_S2432)->p_0[6U] = s_primal_ctx_mul_1(R_17, _S2442) + t_16;
    float2  _S2444 = _S2424;
    Matrix<float, 2, 2>  _S2445 = _S2422;
    SigmaPoints_0 _S2446 = _S2432;
    bool _S2447 = persp_proj_3dgs_ut_1(&_S2446, intrins_12, dist_coeffs_29, image_width_13, image_height_13, &_S2445, &_S2444);
    (&_S2423)->_S2421 = _S2445;
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2448 = _S2423;
    float3  mean_c_13 = s_primal_ctx_mul_1(R_17, mean_14) + t_16;
    float3  _S2449 = make_float3 (_S2434);
    float3  _S2450 = make_float3 (_S2437);
    float3  _S2451 = make_float3 (_S2440);
    float _S2452 = _S2423._S2421.rows[int(0)].y * _S2423._S2421.rows[int(1)].x;
    float det_orig_14 = _S2423._S2421.rows[int(0)].x * _S2423._S2421.rows[int(1)].y - _S2452;
    float _S2453 = _S2423._S2421.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2454 = _S2423._S2421;
    *&(((&_S2454)->rows + (int(0)))->x) = _S2453;
    float _S2455 = _S2423._S2421.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2454)->rows + (int(1)))->y) = _S2455;
    Matrix<float, 2, 2>  _S2456 = _S2454;
    Matrix<float, 2, 2>  _S2457 = _S2454;
    float det_blur_9 = _S2453 * _S2455 - _S2452;
    float _S2458 = det_orig_14 / det_blur_9;
    float _S2459 = det_blur_9 * det_blur_9;
    float _S2460 = s_primal_ctx_max_0(0.0f, _S2458);
    float _S2461 = s_primal_ctx_sqrt_0(_S2460);
    float _S2462 = - in_opacity_13;
    float _S2463 = 1.0f + s_primal_ctx_exp_1(_S2462);
    float _S2464 = 1.0f / _S2463;
    float _S2465 = _S2463 * _S2463;
    float _S2466;
    if(antialiased_13)
    {
        _S2466 = _S2464 * _S2461;
    }
    else
    {
        _S2466 = _S2464;
    }
    float _S2467 = _S2466 / 0.00392156885936856f;
    float _S2468 = 2.0f * s_primal_ctx_log_0(_S2467);
    float _S2469 = s_primal_ctx_sqrt_0(_S2468);
    float _S2470 = _S2456.rows[int(0)].x;
    float _S2471 = _S2457.rows[int(1)].y;
    float _S2472 = s_primal_ctx_dot_0(mean_c_13, mean_c_13) + 9.99999997475242708e-07f;
    float3  _S2473 = - scale_16;
    Matrix<float, 3, 3>  _S2474 = transpose_0(R_17);
    float3  _S2475 = mean_14 - - s_primal_ctx_mul_1(_S2474, t_16);
    float _S2476 = _S2475.x;
    float _S2477 = _S2475.y;
    float _S2478 = _S2475.z;
    float _S2479 = _S2476 * _S2476 + _S2477 * _S2477 + _S2478 * _S2478;
    float _S2480 = s_primal_ctx_sqrt_0(_S2479);
    float x_39 = _S2476 / _S2480;
    float3  _S2481 = make_float3 (x_39);
    float _S2482 = _S2480 * _S2480;
    float y_16 = _S2477 / _S2480;
    float z_13 = _S2478 / _S2480;
    float3  _S2483 = make_float3 (z_13);
    float _S2484 = - y_16;
    float3  _S2485 = make_float3 (_S2484);
    float z2_31 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_39 * x_39 - y_16 * y_16;
    float _S2486 = 2.0f * x_39;
    float fS1_13 = _S2486 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_31 - 0.31539157032966614f;
    float3  _S2487 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_39;
    float3  _S2488 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S2489 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S2490 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S2491 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_31 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S2492 = 1.86588168144226074f * z2_31 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S2492;
    float3  _S2493 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_39;
    float3  _S2494 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S2495 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S2496 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S2497 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_39 * fC1_13 - y_16 * fS1_13);
    float3  _S2498 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_39 * fS1_13 + y_16 * fC1_13);
    float3  _S2499 = make_float3 (pSH9_3);
    float3  _S2500 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2501;
    (&_S2501)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2484) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_39) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S2501)->differential_0 = _S2425;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2502;
    (&_S2502)->primal_0 = _S2500;
    (&_S2502)->differential_0 = _S2425;
    s_bwd_prop_max_0(&_S2501, &_S2502, v_rgb_3);
    float3  _S2503 = _S2498 * _S2501.differential_0;
    float3  _S2504 = (*sh_coeffs_13)[int(15)] * _S2501.differential_0;
    float3  _S2505 = _S2496 * _S2501.differential_0;
    float3  _S2506 = (*sh_coeffs_13)[int(14)] * _S2501.differential_0;
    float3  _S2507 = _S2494 * _S2501.differential_0;
    float3  _S2508 = (*sh_coeffs_13)[int(13)] * _S2501.differential_0;
    float3  _S2509 = _S2493 * _S2501.differential_0;
    float3  _S2510 = (*sh_coeffs_13)[int(12)] * _S2501.differential_0;
    float3  _S2511 = _S2495 * _S2501.differential_0;
    float3  _S2512 = (*sh_coeffs_13)[int(11)] * _S2501.differential_0;
    float3  _S2513 = _S2497 * _S2501.differential_0;
    float3  _S2514 = (*sh_coeffs_13)[int(10)] * _S2501.differential_0;
    float3  _S2515 = _S2499 * _S2501.differential_0;
    float3  _S2516 = (*sh_coeffs_13)[int(9)] * _S2501.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2516.x + _S2516.y + _S2516.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S2504.x + _S2504.y + _S2504.z);
    float _S2517 = _S2514.x + _S2514.y + _S2514.z;
    float _S2518 = _S2506.x + _S2506.y + _S2506.z;
    float _S2519 = _S2512.x + _S2512.y + _S2512.z;
    float _S2520 = _S2508.x + _S2508.y + _S2508.z;
    float _S2521 = _S2510.x + _S2510.y + _S2510.z;
    float _S2522 = - s_diff_fC2_T_3;
    float3  _S2523 = _S2490 * _S2501.differential_0;
    float3  _S2524 = (*sh_coeffs_13)[int(8)] * _S2501.differential_0;
    float3  _S2525 = _S2488 * _S2501.differential_0;
    float3  _S2526 = (*sh_coeffs_13)[int(7)] * _S2501.differential_0;
    float3  _S2527 = _S2487 * _S2501.differential_0;
    float3  _S2528 = (*sh_coeffs_13)[int(6)] * _S2501.differential_0;
    float3  _S2529 = _S2489 * _S2501.differential_0;
    float3  _S2530 = (*sh_coeffs_13)[int(5)] * _S2501.differential_0;
    float3  _S2531 = _S2491 * _S2501.differential_0;
    float3  _S2532 = (*sh_coeffs_13)[int(4)] * _S2501.differential_0;
    float _S2533 = _S2530.x + _S2530.y + _S2530.z;
    float _S2534 = _S2526.x + _S2526.y + _S2526.z;
    float _S2535 = fTmp1B_13 * _S2517 + x_39 * s_diff_fS2_T_3 + y_16 * _S2522 + 0.54627424478530884f * (_S2532.x + _S2532.y + _S2532.z);
    float _S2536 = fTmp1B_13 * _S2518 + y_16 * s_diff_fS2_T_3 + x_39 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2524.x + _S2524.y + _S2524.z);
    float _S2537 = y_16 * - _S2536;
    float _S2538 = x_39 * _S2536;
    float _S2539 = z_13 * (1.86588168144226074f * (z_13 * _S2521) + -2.28522896766662598f * (y_16 * _S2519 + x_39 * _S2520) + 0.94617468118667603f * (_S2528.x + _S2528.y + _S2528.z));
    float3  _S2540 = make_float3 (0.48860251903533936f) * _S2501.differential_0;
    float3  _S2541 = - _S2540;
    float3  _S2542 = _S2481 * _S2541;
    float3  _S2543 = (*sh_coeffs_13)[int(3)] * _S2541;
    float3  _S2544 = _S2483 * _S2540;
    float3  _S2545 = (*sh_coeffs_13)[int(2)] * _S2540;
    float3  _S2546 = _S2485 * _S2540;
    float3  _S2547 = (*sh_coeffs_13)[int(1)] * _S2540;
    float _S2548 = (_S2492 * _S2521 + 1.44530570507049561f * (fS1_13 * _S2517 + fC1_13 * _S2518) + -1.09254848957061768f * (y_16 * _S2533 + x_39 * _S2534) + _S2539 + _S2539 + _S2545.x + _S2545.y + _S2545.z) / _S2482;
    float _S2549 = _S2480 * _S2548;
    float _S2550 = (fTmp0C_13 * _S2519 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2522 + fTmp0B_13 * _S2533 + _S2486 * _S2535 + _S2537 + _S2537 + - (_S2547.x + _S2547.y + _S2547.z)) / _S2482;
    float _S2551 = _S2480 * _S2550;
    float _S2552 = (fTmp0C_13 * _S2520 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2534 + 2.0f * (y_16 * _S2535) + _S2538 + _S2538 + _S2543.x + _S2543.y + _S2543.z) / _S2482;
    float _S2553 = _S2480 * _S2552;
    float _S2554 = _S2478 * - _S2548 + _S2477 * - _S2550 + _S2476 * - _S2552;
    DiffPair_float_0 _S2555;
    (&_S2555)->primal_0 = _S2479;
    (&_S2555)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2555, _S2554);
    float _S2556 = _S2478 * _S2555.differential_0;
    float _S2557 = _S2477 * _S2555.differential_0;
    float _S2558 = _S2476 * _S2555.differential_0;
    float3  _S2559 = make_float3 (0.282094806432724f) * _S2501.differential_0;
    float3  _S2560 = make_float3 (_S2553 + _S2558 + _S2558, _S2551 + _S2557 + _S2557, _S2549 + _S2556 + _S2556);
    float3  _S2561 = - - _S2560;
    Matrix<float, 3, 3>  _S2562 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2563;
    (&_S2563)->primal_0 = _S2474;
    (&_S2563)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2564;
    (&_S2564)->primal_0 = t_16;
    (&_S2564)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2563, &_S2564, _S2561);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2565 = _S2564;
    Matrix<float, 3, 3>  _S2566 = transpose_0(_S2563.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2567;
    (&_S2567)->primal_0 = _S2473;
    (&_S2567)->differential_0 = _S2425;
    s_bwd_prop_exp_1(&_S2567, v_conic_3);
    float3  _S2568 = - _S2567.differential_0;
    float _S2569 = 0.5f * v_depth_3;
    DiffPair_float_0 _S2570;
    (&_S2570)->primal_0 = _S2472;
    (&_S2570)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2570, _S2569);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2571;
    (&_S2571)->primal_0 = mean_c_13;
    (&_S2571)->differential_0 = _S2425;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2572;
    (&_S2572)->primal_0 = mean_c_13;
    (&_S2572)->differential_0 = _S2425;
    s_bwd_prop_dot_0(&_S2571, &_S2572, _S2570.differential_0);
    DiffPair_float_0 _S2573;
    (&_S2573)->primal_0 = _S2471;
    (&_S2573)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2573, 0.0f);
    DiffPair_float_0 _S2574;
    (&_S2574)->primal_0 = _S2470;
    (&_S2574)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2574, 0.0f);
    DiffPair_float_0 _S2575;
    (&_S2575)->primal_0 = 3.32999992370605469f;
    (&_S2575)->differential_0 = 0.0f;
    DiffPair_float_0 _S2576;
    (&_S2576)->primal_0 = _S2469;
    (&_S2576)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2575, &_S2576, 0.0f);
    DiffPair_float_0 _S2577;
    (&_S2577)->primal_0 = _S2468;
    (&_S2577)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2577, _S2576.differential_0);
    float _S2578 = 2.0f * _S2577.differential_0;
    DiffPair_float_0 _S2579;
    (&_S2579)->primal_0 = _S2467;
    (&_S2579)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2579, _S2578);
    float _S2580 = v_opacity_3 + 254.9999847412109375f * _S2579.differential_0;
    FixedArray<float3 , 16>  _S2581;
    _S2581[int(0)] = _S2425;
    _S2581[int(1)] = _S2425;
    _S2581[int(2)] = _S2425;
    _S2581[int(3)] = _S2425;
    _S2581[int(4)] = _S2425;
    _S2581[int(5)] = _S2425;
    _S2581[int(6)] = _S2425;
    _S2581[int(7)] = _S2425;
    _S2581[int(8)] = _S2425;
    _S2581[int(9)] = _S2425;
    _S2581[int(10)] = _S2425;
    _S2581[int(11)] = _S2425;
    _S2581[int(12)] = _S2425;
    _S2581[int(13)] = _S2425;
    _S2581[int(14)] = _S2425;
    _S2581[int(15)] = _S2425;
    _S2581[int(7)] = _S2525;
    _S2581[int(0)] = _S2559;
    _S2581[int(1)] = _S2546;
    _S2581[int(2)] = _S2544;
    _S2581[int(3)] = _S2542;
    _S2581[int(4)] = _S2531;
    _S2581[int(5)] = _S2529;
    _S2581[int(6)] = _S2527;
    _S2581[int(15)] = _S2503;
    _S2581[int(8)] = _S2523;
    _S2581[int(9)] = _S2515;
    _S2581[int(10)] = _S2513;
    _S2581[int(11)] = _S2511;
    _S2581[int(12)] = _S2509;
    _S2581[int(13)] = _S2507;
    _S2581[int(14)] = _S2505;
    float3  _S2582 = _S2581[int(0)];
    float3  _S2583 = _S2581[int(1)];
    float3  _S2584 = _S2581[int(2)];
    float3  _S2585 = _S2581[int(3)];
    float3  _S2586 = _S2581[int(4)];
    float3  _S2587 = _S2581[int(5)];
    float3  _S2588 = _S2581[int(6)];
    float3  _S2589 = _S2581[int(7)];
    float3  _S2590 = _S2581[int(8)];
    float3  _S2591 = _S2581[int(9)];
    float3  _S2592 = _S2581[int(10)];
    float3  _S2593 = _S2581[int(11)];
    float3  _S2594 = _S2581[int(12)];
    float3  _S2595 = _S2581[int(13)];
    float3  _S2596 = _S2581[int(14)];
    float3  _S2597 = _S2581[int(15)];
    float3  _S2598 = _S2572.differential_0 + _S2571.differential_0;
    float2  _S2599 = make_float2 (0.0f, _S2573.differential_0);
    float2  _S2600 = make_float2 (_S2574.differential_0, 0.0f);
    float _S2601;
    if(antialiased_13)
    {
        float _S2602 = _S2464 * _S2580;
        _S2466 = _S2461 * _S2580;
        _S2601 = _S2602;
    }
    else
    {
        _S2466 = _S2580;
        _S2601 = 0.0f;
    }
    float _S2603 = - (_S2466 / _S2465);
    DiffPair_float_0 _S2604;
    (&_S2604)->primal_0 = _S2462;
    (&_S2604)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2604, _S2603);
    float _S2605 = - _S2604.differential_0;
    DiffPair_float_0 _S2606;
    (&_S2606)->primal_0 = _S2460;
    (&_S2606)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2606, _S2601);
    DiffPair_float_0 _S2607;
    (&_S2607)->primal_0 = 0.0f;
    (&_S2607)->differential_0 = 0.0f;
    DiffPair_float_0 _S2608;
    (&_S2608)->primal_0 = _S2458;
    (&_S2608)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2607, &_S2608, _S2606.differential_0);
    float _S2609 = _S2608.differential_0 / _S2459;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2609;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2609;
    float _S2610 = - s_diff_det_blur_T_0;
    float _S2611 = _S2453 * s_diff_det_blur_T_0;
    float _S2612 = _S2455 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2613 = _S2422;
    _S2613[int(1)] = _S2599;
    _S2613[int(0)] = _S2600;
    float _S2614 = _S2612 + _S2613.rows[int(0)].x;
    float _S2615 = _S2610 + - s_diff_det_orig_T_3;
    float _S2616 = _S2448._S2421.rows[int(0)].y * _S2615;
    float _S2617 = _S2448._S2421.rows[int(1)].x * _S2615;
    float _S2618 = _S2448._S2421.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2619 = _S2611 + _S2613.rows[int(1)].y + _S2448._S2421.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2620 = _S2424;
    *&((&_S2620)->x) = _S2616;
    *&((&_S2620)->y) = _S2619;
    float _S2621 = _S2614 + _S2618;
    float2  _S2622 = _S2424;
    *&((&_S2622)->y) = _S2617;
    *&((&_S2622)->x) = _S2621;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2623;
    (&_S2623)->primal_0 = R_17;
    (&_S2623)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2624;
    (&_S2624)->primal_0 = _S2442;
    (&_S2624)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2623, &_S2624, _S2425);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2625;
    (&_S2625)->primal_0 = R_17;
    (&_S2625)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2626;
    (&_S2626)->primal_0 = _S2439;
    (&_S2626)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2625, &_S2626, _S2425);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2627;
    (&_S2627)->primal_0 = R_17;
    (&_S2627)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2628;
    (&_S2628)->primal_0 = _S2436;
    (&_S2628)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2627, &_S2628, _S2425);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2629;
    (&_S2629)->primal_0 = R_17;
    (&_S2629)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2630;
    (&_S2630)->primal_0 = _S2441;
    (&_S2630)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2629, &_S2630, _S2425);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2631;
    (&_S2631)->primal_0 = R_17;
    (&_S2631)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2632;
    (&_S2632)->primal_0 = _S2438;
    (&_S2632)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2631, &_S2632, _S2425);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2633;
    (&_S2633)->primal_0 = R_17;
    (&_S2633)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2634;
    (&_S2634)->primal_0 = _S2435;
    (&_S2634)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2633, &_S2634, _S2425);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2635;
    (&_S2635)->primal_0 = R_17;
    (&_S2635)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2636;
    (&_S2636)->primal_0 = _S2443.p_0[0U];
    (&_S2636)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2635, &_S2636, _S2425);
    float3  _S2637 = - _S2624.differential_0 + _S2630.differential_0;
    float3  _S2638 = _S2451 * _S2637;
    float3  _S2639 = _S2428.rows[2U] * _S2637;
    float _S2640 = _S2433 * (_S2639.x + _S2639.y + _S2639.z);
    float3  _S2641 = - _S2626.differential_0 + _S2632.differential_0;
    float3  _S2642 = _S2450 * _S2641;
    float3  _S2643 = _S2428.rows[1U] * _S2641;
    float _S2644 = _S2433 * (_S2643.x + _S2643.y + _S2643.z);
    float3  _S2645 = - _S2628.differential_0 + _S2634.differential_0;
    float3  _S2646 = _S2449 * _S2645;
    float3  _S2647 = _S2428.rows[0U] * _S2645;
    float _S2648 = _S2433 * (_S2647.x + _S2647.y + _S2647.z);
    Matrix<float, 3, 3>  _S2649 = _S2562;
    _S2649[2U] = _S2638;
    _S2649[1U] = _S2642;
    _S2649[0U] = _S2646;
    Matrix<float, 3, 3>  _S2650 = transpose_0(transpose_0(_S2649));
    float _S2651 = 2.0f * - _S2650.rows[int(2)].z;
    float _S2652 = 2.0f * _S2650.rows[int(2)].y;
    float _S2653 = 2.0f * _S2650.rows[int(2)].x;
    float _S2654 = 2.0f * _S2650.rows[int(1)].z;
    float _S2655 = 2.0f * - _S2650.rows[int(1)].y;
    float _S2656 = 2.0f * _S2650.rows[int(1)].x;
    float _S2657 = 2.0f * _S2650.rows[int(0)].z;
    float _S2658 = 2.0f * _S2650.rows[int(0)].y;
    float _S2659 = 2.0f * - _S2650.rows[int(0)].x;
    float _S2660 = - _S2656 + _S2658;
    float _S2661 = _S2653 + - _S2657;
    float _S2662 = - _S2652 + _S2654;
    float _S2663 = _S2652 + _S2654;
    float _S2664 = _S2653 + _S2657;
    float _S2665 = _S2656 + _S2658;
    float _S2666 = quat_17.w * (_S2655 + _S2659);
    float _S2667 = quat_17.z * (_S2651 + _S2659);
    float _S2668 = quat_17.y * (_S2651 + _S2655);
    float _S2669 = quat_17.x * _S2660 + quat_17.z * _S2663 + quat_17.y * _S2664 + _S2666 + _S2666;
    float _S2670 = quat_17.x * _S2661 + quat_17.w * _S2663 + quat_17.y * _S2665 + _S2667 + _S2667;
    float _S2671 = quat_17.x * _S2662 + quat_17.w * _S2664 + quat_17.z * _S2665 + _S2668 + _S2668;
    float _S2672 = quat_17.w * _S2660 + quat_17.z * _S2661 + quat_17.y * _S2662;
    float3  _S2673 = _S2425;
    *&((&_S2673)->z) = _S2640;
    *&((&_S2673)->y) = _S2644;
    *&((&_S2673)->x) = _S2648;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2674;
    (&_S2674)->primal_0 = scale_16;
    (&_S2674)->differential_0 = _S2425;
    s_bwd_prop_exp_1(&_S2674, _S2673);
    Matrix<float, 2, 2>  _S2675 = _S2422;
    _S2675[int(1)] = _S2620;
    _S2675[int(0)] = _S2622;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2676;
    (&_S2676)->primal_0 = R_17;
    (&_S2676)->differential_0 = _S2562;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2677;
    (&_S2677)->primal_0 = mean_14;
    (&_S2677)->differential_0 = _S2425;
    s_bwd_prop_mul_1(&_S2676, &_S2677, _S2598);
    float3  _S2678 = _S2598 + _S2565.differential_0;
    Matrix<float, 3, 3>  _S2679 = _S2623.differential_0 + _S2625.differential_0 + _S2627.differential_0 + _S2629.differential_0 + _S2631.differential_0 + _S2633.differential_0 + _S2635.differential_0 + _S2676.differential_0 + _S2566;
    float3  _S2680 = _S2674.differential_0 + _S2568;
    float4  _S2681 = make_float4 (0.0f);
    *&((&_S2681)->w) = _S2669;
    *&((&_S2681)->z) = _S2670;
    *&((&_S2681)->y) = _S2671;
    *&((&_S2681)->x) = _S2672;
    float4  _S2682 = _S2681;
    float3  _S2683 = _S2624.differential_0 + _S2630.differential_0 + _S2626.differential_0 + _S2632.differential_0 + _S2628.differential_0 + _S2634.differential_0 + _S2677.differential_0 + _S2560;
    *v_mean_3 = _S2683;
    *v_quat_3 = _S2682;
    *v_scale_3 = _S2680;
    *v_in_opacity_3 = _S2605;
    (*v_sh_coeffs_3)[int(0)] = _S2582;
    (*v_sh_coeffs_3)[int(1)] = _S2583;
    (*v_sh_coeffs_3)[int(2)] = _S2584;
    (*v_sh_coeffs_3)[int(3)] = _S2585;
    (*v_sh_coeffs_3)[int(4)] = _S2586;
    (*v_sh_coeffs_3)[int(5)] = _S2587;
    (*v_sh_coeffs_3)[int(6)] = _S2588;
    (*v_sh_coeffs_3)[int(7)] = _S2589;
    (*v_sh_coeffs_3)[int(8)] = _S2590;
    (*v_sh_coeffs_3)[int(9)] = _S2591;
    (*v_sh_coeffs_3)[int(10)] = _S2592;
    (*v_sh_coeffs_3)[int(11)] = _S2593;
    (*v_sh_coeffs_3)[int(12)] = _S2594;
    (*v_sh_coeffs_3)[int(13)] = _S2595;
    (*v_sh_coeffs_3)[int(14)] = _S2596;
    (*v_sh_coeffs_3)[int(15)] = _S2597;
    *v_R_4 = _S2679;
    *v_t_4 = _S2678;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2684;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_15, float4  quat_18, float3  scale_17, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_18, float3  t_17, float fx_23, float fy_23, float cx_18, float cy_18, FixedArray<float, 10>  * dist_coeffs_30, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2685 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2686;
    (&_S2686)->_S2684 = _S2685;
    float2  _S2687 = make_float2 (0.0f);
    float3  _S2688 = make_float3 (0.0f);
    float4  intrins_13 = make_float4 (fx_23, fy_23, cx_18, cy_18);
    float3  _S2689 = s_primal_ctx_exp_0(scale_17);
    float _S2690 = quat_18.y;
    float x2_18 = _S2690 * _S2690;
    float y2_18 = quat_18.z * quat_18.z;
    float z2_32 = quat_18.w * quat_18.w;
    float xy_18 = quat_18.y * quat_18.z;
    float xz_18 = quat_18.y * quat_18.w;
    float yz_18 = quat_18.z * quat_18.w;
    float wx_18 = quat_18.x * quat_18.y;
    float wy_18 = quat_18.x * quat_18.z;
    float wz_18 = quat_18.x * quat_18.w;
    Matrix<float, 3, 3>  _S2691 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_18 + z2_32), 2.0f * (xy_18 + wz_18), 2.0f * (xz_18 - wy_18), 2.0f * (xy_18 - wz_18), 1.0f - 2.0f * (x2_18 + z2_32), 2.0f * (yz_18 + wx_18), 2.0f * (xz_18 + wy_18), 2.0f * (yz_18 - wx_18), 1.0f - 2.0f * (x2_18 + y2_18))));
    FixedArray<float3 , 7>  _S2692 = {
        _S2688, _S2688, _S2688, _S2688, _S2688, _S2688, _S2688
    };
    FixedArray<float, 7>  _S2693 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2694;
    (&_S2694)->p_0 = _S2692;
    (&_S2694)->w_mean_0 = _S2693;
    (&_S2694)->w_cov_0 = _S2693;
    (&_S2694)->p_0[int(0)] = mean_15;
    SigmaPoints_0 _S2695 = _S2694;
    (&_S2695)->w_mean_0[int(0)] = 0.0f;
    (&_S2695)->w_cov_0[int(0)] = 2.0f;
    float _S2696 = s_primal_ctx_sqrt_0(3.0f);
    float _S2697 = _S2696 * _S2689.x;
    float3  delta_15 = make_float3 (_S2697) * _S2691.rows[0U];
    float3  _S2698 = mean_15 + delta_15;
    float3  _S2699 = mean_15 - delta_15;
    float _S2700 = _S2696 * _S2689.y;
    float3  delta_16 = make_float3 (_S2700) * _S2691.rows[1U];
    float3  _S2701 = mean_15 + delta_16;
    float3  _S2702 = mean_15 - delta_16;
    float _S2703 = _S2696 * _S2689.z;
    float3  delta_17 = make_float3 (_S2703) * _S2691.rows[2U];
    float3  _S2704 = mean_15 + delta_17;
    float3  _S2705 = mean_15 - delta_17;
    (&_S2695)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2695)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2695)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2695)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2695)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2695)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2695)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2695)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2695)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2695)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2695)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2695)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2706 = _S2694;
    (&_S2695)->p_0[0U] = s_primal_ctx_mul_1(R_18, _S2694.p_0[0U]) + t_17;
    (&_S2695)->p_0[1U] = s_primal_ctx_mul_1(R_18, _S2698) + t_17;
    (&_S2695)->p_0[2U] = s_primal_ctx_mul_1(R_18, _S2701) + t_17;
    (&_S2695)->p_0[3U] = s_primal_ctx_mul_1(R_18, _S2704) + t_17;
    (&_S2695)->p_0[4U] = s_primal_ctx_mul_1(R_18, _S2699) + t_17;
    (&_S2695)->p_0[5U] = s_primal_ctx_mul_1(R_18, _S2702) + t_17;
    (&_S2695)->p_0[6U] = s_primal_ctx_mul_1(R_18, _S2705) + t_17;
    float2  _S2707 = _S2687;
    Matrix<float, 2, 2>  _S2708 = _S2685;
    SigmaPoints_0 _S2709 = _S2695;
    bool _S2710 = fisheye_proj_3dgs_ut_1(&_S2709, intrins_13, dist_coeffs_30, &_S2708, &_S2707);
    (&_S2686)->_S2684 = _S2708;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2711 = _S2686;
    float3  mean_c_14 = s_primal_ctx_mul_1(R_18, mean_15) + t_17;
    float3  _S2712 = make_float3 (_S2697);
    float3  _S2713 = make_float3 (_S2700);
    float3  _S2714 = make_float3 (_S2703);
    float _S2715 = _S2686._S2684.rows[int(0)].y * _S2686._S2684.rows[int(1)].x;
    float det_orig_15 = _S2686._S2684.rows[int(0)].x * _S2686._S2684.rows[int(1)].y - _S2715;
    float _S2716 = _S2686._S2684.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2717 = _S2686._S2684;
    *&(((&_S2717)->rows + (int(0)))->x) = _S2716;
    float _S2718 = _S2686._S2684.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2717)->rows + (int(1)))->y) = _S2718;
    Matrix<float, 2, 2>  _S2719 = _S2717;
    Matrix<float, 2, 2>  _S2720 = _S2717;
    float det_blur_10 = _S2716 * _S2718 - _S2715;
    float _S2721 = det_orig_15 / det_blur_10;
    float _S2722 = det_blur_10 * det_blur_10;
    float _S2723 = s_primal_ctx_max_0(0.0f, _S2721);
    float _S2724 = s_primal_ctx_sqrt_0(_S2723);
    float _S2725 = - in_opacity_14;
    float _S2726 = 1.0f + s_primal_ctx_exp_1(_S2725);
    float _S2727 = 1.0f / _S2726;
    float _S2728 = _S2726 * _S2726;
    float _S2729;
    if(antialiased_14)
    {
        _S2729 = _S2727 * _S2724;
    }
    else
    {
        _S2729 = _S2727;
    }
    float _S2730 = _S2729 / 0.00392156885936856f;
    float _S2731 = 2.0f * s_primal_ctx_log_0(_S2730);
    float _S2732 = s_primal_ctx_sqrt_0(_S2731);
    float _S2733 = _S2719.rows[int(0)].x;
    float _S2734 = _S2720.rows[int(1)].y;
    float _S2735 = s_primal_ctx_dot_0(mean_c_14, mean_c_14) + 9.99999997475242708e-07f;
    float3  _S2736 = - scale_17;
    Matrix<float, 3, 3>  _S2737 = transpose_0(R_18);
    float3  _S2738 = mean_15 - - s_primal_ctx_mul_1(_S2737, t_17);
    float _S2739 = _S2738.x;
    float _S2740 = _S2738.y;
    float _S2741 = _S2738.z;
    float _S2742 = _S2739 * _S2739 + _S2740 * _S2740 + _S2741 * _S2741;
    float _S2743 = s_primal_ctx_sqrt_0(_S2742);
    float x_40 = _S2739 / _S2743;
    float3  _S2744 = make_float3 (x_40);
    float _S2745 = _S2743 * _S2743;
    float y_17 = _S2740 / _S2743;
    float z_14 = _S2741 / _S2743;
    float3  _S2746 = make_float3 (z_14);
    float _S2747 = - y_17;
    float3  _S2748 = make_float3 (_S2747);
    float z2_33 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_40 * x_40 - y_17 * y_17;
    float _S2749 = 2.0f * x_40;
    float fS1_14 = _S2749 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_33 - 0.31539157032966614f;
    float3  _S2750 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_40;
    float3  _S2751 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2752 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2753 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2754 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_33 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2755 = 1.86588168144226074f * z2_33 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2755;
    float3  _S2756 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_40;
    float3  _S2757 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2758 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2759 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2760 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_40 * fC1_14 - y_17 * fS1_14);
    float3  _S2761 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_40 * fS1_14 + y_17 * fC1_14);
    float3  _S2762 = make_float3 (pSH9_4);
    float3  _S2763 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2764;
    (&_S2764)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2747) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_40) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2764)->differential_0 = _S2688;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2765;
    (&_S2765)->primal_0 = _S2763;
    (&_S2765)->differential_0 = _S2688;
    s_bwd_prop_max_0(&_S2764, &_S2765, v_rgb_4);
    float3  _S2766 = _S2761 * _S2764.differential_0;
    float3  _S2767 = (*sh_coeffs_14)[int(15)] * _S2764.differential_0;
    float3  _S2768 = _S2759 * _S2764.differential_0;
    float3  _S2769 = (*sh_coeffs_14)[int(14)] * _S2764.differential_0;
    float3  _S2770 = _S2757 * _S2764.differential_0;
    float3  _S2771 = (*sh_coeffs_14)[int(13)] * _S2764.differential_0;
    float3  _S2772 = _S2756 * _S2764.differential_0;
    float3  _S2773 = (*sh_coeffs_14)[int(12)] * _S2764.differential_0;
    float3  _S2774 = _S2758 * _S2764.differential_0;
    float3  _S2775 = (*sh_coeffs_14)[int(11)] * _S2764.differential_0;
    float3  _S2776 = _S2760 * _S2764.differential_0;
    float3  _S2777 = (*sh_coeffs_14)[int(10)] * _S2764.differential_0;
    float3  _S2778 = _S2762 * _S2764.differential_0;
    float3  _S2779 = (*sh_coeffs_14)[int(9)] * _S2764.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2779.x + _S2779.y + _S2779.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2767.x + _S2767.y + _S2767.z);
    float _S2780 = _S2777.x + _S2777.y + _S2777.z;
    float _S2781 = _S2769.x + _S2769.y + _S2769.z;
    float _S2782 = _S2775.x + _S2775.y + _S2775.z;
    float _S2783 = _S2771.x + _S2771.y + _S2771.z;
    float _S2784 = _S2773.x + _S2773.y + _S2773.z;
    float _S2785 = - s_diff_fC2_T_4;
    float3  _S2786 = _S2753 * _S2764.differential_0;
    float3  _S2787 = (*sh_coeffs_14)[int(8)] * _S2764.differential_0;
    float3  _S2788 = _S2751 * _S2764.differential_0;
    float3  _S2789 = (*sh_coeffs_14)[int(7)] * _S2764.differential_0;
    float3  _S2790 = _S2750 * _S2764.differential_0;
    float3  _S2791 = (*sh_coeffs_14)[int(6)] * _S2764.differential_0;
    float3  _S2792 = _S2752 * _S2764.differential_0;
    float3  _S2793 = (*sh_coeffs_14)[int(5)] * _S2764.differential_0;
    float3  _S2794 = _S2754 * _S2764.differential_0;
    float3  _S2795 = (*sh_coeffs_14)[int(4)] * _S2764.differential_0;
    float _S2796 = _S2793.x + _S2793.y + _S2793.z;
    float _S2797 = _S2789.x + _S2789.y + _S2789.z;
    float _S2798 = fTmp1B_14 * _S2780 + x_40 * s_diff_fS2_T_4 + y_17 * _S2785 + 0.54627424478530884f * (_S2795.x + _S2795.y + _S2795.z);
    float _S2799 = fTmp1B_14 * _S2781 + y_17 * s_diff_fS2_T_4 + x_40 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2787.x + _S2787.y + _S2787.z);
    float _S2800 = y_17 * - _S2799;
    float _S2801 = x_40 * _S2799;
    float _S2802 = z_14 * (1.86588168144226074f * (z_14 * _S2784) + -2.28522896766662598f * (y_17 * _S2782 + x_40 * _S2783) + 0.94617468118667603f * (_S2791.x + _S2791.y + _S2791.z));
    float3  _S2803 = make_float3 (0.48860251903533936f) * _S2764.differential_0;
    float3  _S2804 = - _S2803;
    float3  _S2805 = _S2744 * _S2804;
    float3  _S2806 = (*sh_coeffs_14)[int(3)] * _S2804;
    float3  _S2807 = _S2746 * _S2803;
    float3  _S2808 = (*sh_coeffs_14)[int(2)] * _S2803;
    float3  _S2809 = _S2748 * _S2803;
    float3  _S2810 = (*sh_coeffs_14)[int(1)] * _S2803;
    float _S2811 = (_S2755 * _S2784 + 1.44530570507049561f * (fS1_14 * _S2780 + fC1_14 * _S2781) + -1.09254848957061768f * (y_17 * _S2796 + x_40 * _S2797) + _S2802 + _S2802 + _S2808.x + _S2808.y + _S2808.z) / _S2745;
    float _S2812 = _S2743 * _S2811;
    float _S2813 = (fTmp0C_14 * _S2782 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2785 + fTmp0B_14 * _S2796 + _S2749 * _S2798 + _S2800 + _S2800 + - (_S2810.x + _S2810.y + _S2810.z)) / _S2745;
    float _S2814 = _S2743 * _S2813;
    float _S2815 = (fTmp0C_14 * _S2783 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2797 + 2.0f * (y_17 * _S2798) + _S2801 + _S2801 + _S2806.x + _S2806.y + _S2806.z) / _S2745;
    float _S2816 = _S2743 * _S2815;
    float _S2817 = _S2741 * - _S2811 + _S2740 * - _S2813 + _S2739 * - _S2815;
    DiffPair_float_0 _S2818;
    (&_S2818)->primal_0 = _S2742;
    (&_S2818)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2818, _S2817);
    float _S2819 = _S2741 * _S2818.differential_0;
    float _S2820 = _S2740 * _S2818.differential_0;
    float _S2821 = _S2739 * _S2818.differential_0;
    float3  _S2822 = make_float3 (0.282094806432724f) * _S2764.differential_0;
    float3  _S2823 = make_float3 (_S2816 + _S2821 + _S2821, _S2814 + _S2820 + _S2820, _S2812 + _S2819 + _S2819);
    float3  _S2824 = - - _S2823;
    Matrix<float, 3, 3>  _S2825 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2826;
    (&_S2826)->primal_0 = _S2737;
    (&_S2826)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2827;
    (&_S2827)->primal_0 = t_17;
    (&_S2827)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2826, &_S2827, _S2824);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2828 = _S2827;
    Matrix<float, 3, 3>  _S2829 = transpose_0(_S2826.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2830;
    (&_S2830)->primal_0 = _S2736;
    (&_S2830)->differential_0 = _S2688;
    s_bwd_prop_exp_1(&_S2830, v_conic_4);
    float3  _S2831 = - _S2830.differential_0;
    float _S2832 = 0.5f * v_depth_4;
    DiffPair_float_0 _S2833;
    (&_S2833)->primal_0 = _S2735;
    (&_S2833)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2833, _S2832);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2834;
    (&_S2834)->primal_0 = mean_c_14;
    (&_S2834)->differential_0 = _S2688;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2835;
    (&_S2835)->primal_0 = mean_c_14;
    (&_S2835)->differential_0 = _S2688;
    s_bwd_prop_dot_0(&_S2834, &_S2835, _S2833.differential_0);
    DiffPair_float_0 _S2836;
    (&_S2836)->primal_0 = _S2734;
    (&_S2836)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2836, 0.0f);
    DiffPair_float_0 _S2837;
    (&_S2837)->primal_0 = _S2733;
    (&_S2837)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2837, 0.0f);
    DiffPair_float_0 _S2838;
    (&_S2838)->primal_0 = 3.32999992370605469f;
    (&_S2838)->differential_0 = 0.0f;
    DiffPair_float_0 _S2839;
    (&_S2839)->primal_0 = _S2732;
    (&_S2839)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2838, &_S2839, 0.0f);
    DiffPair_float_0 _S2840;
    (&_S2840)->primal_0 = _S2731;
    (&_S2840)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2840, _S2839.differential_0);
    float _S2841 = 2.0f * _S2840.differential_0;
    DiffPair_float_0 _S2842;
    (&_S2842)->primal_0 = _S2730;
    (&_S2842)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2842, _S2841);
    float _S2843 = v_opacity_4 + 254.9999847412109375f * _S2842.differential_0;
    FixedArray<float3 , 16>  _S2844;
    _S2844[int(0)] = _S2688;
    _S2844[int(1)] = _S2688;
    _S2844[int(2)] = _S2688;
    _S2844[int(3)] = _S2688;
    _S2844[int(4)] = _S2688;
    _S2844[int(5)] = _S2688;
    _S2844[int(6)] = _S2688;
    _S2844[int(7)] = _S2688;
    _S2844[int(8)] = _S2688;
    _S2844[int(9)] = _S2688;
    _S2844[int(10)] = _S2688;
    _S2844[int(11)] = _S2688;
    _S2844[int(12)] = _S2688;
    _S2844[int(13)] = _S2688;
    _S2844[int(14)] = _S2688;
    _S2844[int(15)] = _S2688;
    _S2844[int(7)] = _S2788;
    _S2844[int(0)] = _S2822;
    _S2844[int(1)] = _S2809;
    _S2844[int(2)] = _S2807;
    _S2844[int(3)] = _S2805;
    _S2844[int(4)] = _S2794;
    _S2844[int(5)] = _S2792;
    _S2844[int(6)] = _S2790;
    _S2844[int(15)] = _S2766;
    _S2844[int(8)] = _S2786;
    _S2844[int(9)] = _S2778;
    _S2844[int(10)] = _S2776;
    _S2844[int(11)] = _S2774;
    _S2844[int(12)] = _S2772;
    _S2844[int(13)] = _S2770;
    _S2844[int(14)] = _S2768;
    float3  _S2845 = _S2844[int(0)];
    float3  _S2846 = _S2844[int(1)];
    float3  _S2847 = _S2844[int(2)];
    float3  _S2848 = _S2844[int(3)];
    float3  _S2849 = _S2844[int(4)];
    float3  _S2850 = _S2844[int(5)];
    float3  _S2851 = _S2844[int(6)];
    float3  _S2852 = _S2844[int(7)];
    float3  _S2853 = _S2844[int(8)];
    float3  _S2854 = _S2844[int(9)];
    float3  _S2855 = _S2844[int(10)];
    float3  _S2856 = _S2844[int(11)];
    float3  _S2857 = _S2844[int(12)];
    float3  _S2858 = _S2844[int(13)];
    float3  _S2859 = _S2844[int(14)];
    float3  _S2860 = _S2844[int(15)];
    float3  _S2861 = _S2835.differential_0 + _S2834.differential_0;
    float2  _S2862 = make_float2 (0.0f, _S2836.differential_0);
    float2  _S2863 = make_float2 (_S2837.differential_0, 0.0f);
    float _S2864;
    if(antialiased_14)
    {
        float _S2865 = _S2727 * _S2843;
        _S2729 = _S2724 * _S2843;
        _S2864 = _S2865;
    }
    else
    {
        _S2729 = _S2843;
        _S2864 = 0.0f;
    }
    float _S2866 = - (_S2729 / _S2728);
    DiffPair_float_0 _S2867;
    (&_S2867)->primal_0 = _S2725;
    (&_S2867)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2867, _S2866);
    float _S2868 = - _S2867.differential_0;
    DiffPair_float_0 _S2869;
    (&_S2869)->primal_0 = _S2723;
    (&_S2869)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2869, _S2864);
    DiffPair_float_0 _S2870;
    (&_S2870)->primal_0 = 0.0f;
    (&_S2870)->differential_0 = 0.0f;
    DiffPair_float_0 _S2871;
    (&_S2871)->primal_0 = _S2721;
    (&_S2871)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2870, &_S2871, _S2869.differential_0);
    float _S2872 = _S2871.differential_0 / _S2722;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2872;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2872;
    float _S2873 = - s_diff_det_blur_T_1;
    float _S2874 = _S2716 * s_diff_det_blur_T_1;
    float _S2875 = _S2718 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2876 = _S2685;
    _S2876[int(1)] = _S2862;
    _S2876[int(0)] = _S2863;
    float _S2877 = _S2875 + _S2876.rows[int(0)].x;
    float _S2878 = _S2873 + - s_diff_det_orig_T_4;
    float _S2879 = _S2711._S2684.rows[int(0)].y * _S2878;
    float _S2880 = _S2711._S2684.rows[int(1)].x * _S2878;
    float _S2881 = _S2711._S2684.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2882 = _S2874 + _S2876.rows[int(1)].y + _S2711._S2684.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2883 = _S2687;
    *&((&_S2883)->x) = _S2879;
    *&((&_S2883)->y) = _S2882;
    float _S2884 = _S2877 + _S2881;
    float2  _S2885 = _S2687;
    *&((&_S2885)->y) = _S2880;
    *&((&_S2885)->x) = _S2884;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2886;
    (&_S2886)->primal_0 = R_18;
    (&_S2886)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2887;
    (&_S2887)->primal_0 = _S2705;
    (&_S2887)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2886, &_S2887, _S2688);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2888;
    (&_S2888)->primal_0 = R_18;
    (&_S2888)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2889;
    (&_S2889)->primal_0 = _S2702;
    (&_S2889)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2888, &_S2889, _S2688);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2890;
    (&_S2890)->primal_0 = R_18;
    (&_S2890)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2891;
    (&_S2891)->primal_0 = _S2699;
    (&_S2891)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2890, &_S2891, _S2688);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2892;
    (&_S2892)->primal_0 = R_18;
    (&_S2892)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2893;
    (&_S2893)->primal_0 = _S2704;
    (&_S2893)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2892, &_S2893, _S2688);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2894;
    (&_S2894)->primal_0 = R_18;
    (&_S2894)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2895;
    (&_S2895)->primal_0 = _S2701;
    (&_S2895)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2894, &_S2895, _S2688);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2896;
    (&_S2896)->primal_0 = R_18;
    (&_S2896)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2897;
    (&_S2897)->primal_0 = _S2698;
    (&_S2897)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2896, &_S2897, _S2688);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2898;
    (&_S2898)->primal_0 = R_18;
    (&_S2898)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2899;
    (&_S2899)->primal_0 = _S2706.p_0[0U];
    (&_S2899)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2898, &_S2899, _S2688);
    float3  _S2900 = - _S2887.differential_0 + _S2893.differential_0;
    float3  _S2901 = _S2714 * _S2900;
    float3  _S2902 = _S2691.rows[2U] * _S2900;
    float _S2903 = _S2696 * (_S2902.x + _S2902.y + _S2902.z);
    float3  _S2904 = - _S2889.differential_0 + _S2895.differential_0;
    float3  _S2905 = _S2713 * _S2904;
    float3  _S2906 = _S2691.rows[1U] * _S2904;
    float _S2907 = _S2696 * (_S2906.x + _S2906.y + _S2906.z);
    float3  _S2908 = - _S2891.differential_0 + _S2897.differential_0;
    float3  _S2909 = _S2712 * _S2908;
    float3  _S2910 = _S2691.rows[0U] * _S2908;
    float _S2911 = _S2696 * (_S2910.x + _S2910.y + _S2910.z);
    Matrix<float, 3, 3>  _S2912 = _S2825;
    _S2912[2U] = _S2901;
    _S2912[1U] = _S2905;
    _S2912[0U] = _S2909;
    Matrix<float, 3, 3>  _S2913 = transpose_0(transpose_0(_S2912));
    float _S2914 = 2.0f * - _S2913.rows[int(2)].z;
    float _S2915 = 2.0f * _S2913.rows[int(2)].y;
    float _S2916 = 2.0f * _S2913.rows[int(2)].x;
    float _S2917 = 2.0f * _S2913.rows[int(1)].z;
    float _S2918 = 2.0f * - _S2913.rows[int(1)].y;
    float _S2919 = 2.0f * _S2913.rows[int(1)].x;
    float _S2920 = 2.0f * _S2913.rows[int(0)].z;
    float _S2921 = 2.0f * _S2913.rows[int(0)].y;
    float _S2922 = 2.0f * - _S2913.rows[int(0)].x;
    float _S2923 = - _S2919 + _S2921;
    float _S2924 = _S2916 + - _S2920;
    float _S2925 = - _S2915 + _S2917;
    float _S2926 = _S2915 + _S2917;
    float _S2927 = _S2916 + _S2920;
    float _S2928 = _S2919 + _S2921;
    float _S2929 = quat_18.w * (_S2918 + _S2922);
    float _S2930 = quat_18.z * (_S2914 + _S2922);
    float _S2931 = quat_18.y * (_S2914 + _S2918);
    float _S2932 = quat_18.x * _S2923 + quat_18.z * _S2926 + quat_18.y * _S2927 + _S2929 + _S2929;
    float _S2933 = quat_18.x * _S2924 + quat_18.w * _S2926 + quat_18.y * _S2928 + _S2930 + _S2930;
    float _S2934 = quat_18.x * _S2925 + quat_18.w * _S2927 + quat_18.z * _S2928 + _S2931 + _S2931;
    float _S2935 = quat_18.w * _S2923 + quat_18.z * _S2924 + quat_18.y * _S2925;
    float3  _S2936 = _S2688;
    *&((&_S2936)->z) = _S2903;
    *&((&_S2936)->y) = _S2907;
    *&((&_S2936)->x) = _S2911;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2937;
    (&_S2937)->primal_0 = scale_17;
    (&_S2937)->differential_0 = _S2688;
    s_bwd_prop_exp_1(&_S2937, _S2936);
    Matrix<float, 2, 2>  _S2938 = _S2685;
    _S2938[int(1)] = _S2883;
    _S2938[int(0)] = _S2885;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2939;
    (&_S2939)->primal_0 = mean_c_14;
    (&_S2939)->differential_0 = _S2688;
    s_bwd_length_impl_0(&_S2939, 0.0f);
    float3  _S2940 = _S2939.differential_0 + _S2861;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2941;
    (&_S2941)->primal_0 = R_18;
    (&_S2941)->differential_0 = _S2825;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2942;
    (&_S2942)->primal_0 = mean_15;
    (&_S2942)->differential_0 = _S2688;
    s_bwd_prop_mul_1(&_S2941, &_S2942, _S2940);
    float3  _S2943 = _S2940 + _S2828.differential_0;
    Matrix<float, 3, 3>  _S2944 = _S2886.differential_0 + _S2888.differential_0 + _S2890.differential_0 + _S2892.differential_0 + _S2894.differential_0 + _S2896.differential_0 + _S2898.differential_0 + _S2941.differential_0 + _S2829;
    float3  _S2945 = _S2937.differential_0 + _S2831;
    float4  _S2946 = make_float4 (0.0f);
    *&((&_S2946)->w) = _S2932;
    *&((&_S2946)->z) = _S2933;
    *&((&_S2946)->y) = _S2934;
    *&((&_S2946)->x) = _S2935;
    float4  _S2947 = _S2946;
    float3  _S2948 = _S2887.differential_0 + _S2893.differential_0 + _S2889.differential_0 + _S2895.differential_0 + _S2891.differential_0 + _S2897.differential_0 + _S2942.differential_0 + _S2823;
    *v_mean_4 = _S2948;
    *v_quat_4 = _S2947;
    *v_scale_4 = _S2945;
    *v_in_opacity_4 = _S2868;
    (*v_sh_coeffs_4)[int(0)] = _S2845;
    (*v_sh_coeffs_4)[int(1)] = _S2846;
    (*v_sh_coeffs_4)[int(2)] = _S2847;
    (*v_sh_coeffs_4)[int(3)] = _S2848;
    (*v_sh_coeffs_4)[int(4)] = _S2849;
    (*v_sh_coeffs_4)[int(5)] = _S2850;
    (*v_sh_coeffs_4)[int(6)] = _S2851;
    (*v_sh_coeffs_4)[int(7)] = _S2852;
    (*v_sh_coeffs_4)[int(8)] = _S2853;
    (*v_sh_coeffs_4)[int(9)] = _S2854;
    (*v_sh_coeffs_4)[int(10)] = _S2855;
    (*v_sh_coeffs_4)[int(11)] = _S2856;
    (*v_sh_coeffs_4)[int(12)] = _S2857;
    (*v_sh_coeffs_4)[int(13)] = _S2858;
    (*v_sh_coeffs_4)[int(14)] = _S2859;
    (*v_sh_coeffs_4)[int(15)] = _S2860;
    *v_R_5 = _S2944;
    *v_t_5 = _S2943;
    return;
}

inline __device__ Matrix<float, 3, 3>  compute_3dgut_iscl_rot(float4  quat_19, float3  scale_18)
{
    float x_41 = quat_19.y;
    float x2_19 = x_41 * x_41;
    float y2_19 = quat_19.z * quat_19.z;
    float z2_34 = quat_19.w * quat_19.w;
    float xy_19 = quat_19.y * quat_19.z;
    float xz_19 = quat_19.y * quat_19.w;
    float yz_19 = quat_19.z * quat_19.w;
    float wx_19 = quat_19.x * quat_19.y;
    float wy_19 = quat_19.x * quat_19.z;
    float wz_19 = quat_19.x * quat_19.w;
    return mul_4(makeMatrix<float, 3, 3> (scale_18.x, 0.0f, 0.0f, 0.0f, scale_18.y, 0.0f, 0.0f, 0.0f, scale_18.z), transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_19 + z2_34), 2.0f * (xy_19 + wz_19), 2.0f * (xz_19 - wy_19), 2.0f * (xy_19 - wz_19), 1.0f - 2.0f * (x2_19 + z2_34), 2.0f * (yz_19 + wx_19), 2.0f * (xz_19 + wy_19), 2.0f * (yz_19 - wx_19), 1.0f - 2.0f * (x2_19 + y2_19)))));
}

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_6)
{
    float _S2949 = (*dpquat_0).primal_0.y;
    float x2_20 = _S2949 * _S2949;
    float y2_20 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_35 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_20 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_20 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_20 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_20 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_20 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_20 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2950 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20))));
    Matrix<float, 3, 3>  _S2951 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2952;
    (&_S2952)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2952)->differential_0 = _S2951;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2953;
    (&_S2953)->primal_0 = _S2950;
    (&_S2953)->differential_0 = _S2951;
    s_bwd_prop_mul_4(&_S2952, &_S2953, _s_dOut_6);
    Matrix<float, 3, 3>  _S2954 = transpose_0(transpose_0(_S2953.differential_0));
    float _S2955 = 2.0f * - _S2954.rows[int(2)].z;
    float _S2956 = 2.0f * _S2954.rows[int(2)].y;
    float _S2957 = 2.0f * _S2954.rows[int(2)].x;
    float _S2958 = 2.0f * _S2954.rows[int(1)].z;
    float _S2959 = 2.0f * - _S2954.rows[int(1)].y;
    float _S2960 = 2.0f * _S2954.rows[int(1)].x;
    float _S2961 = 2.0f * _S2954.rows[int(0)].z;
    float _S2962 = 2.0f * _S2954.rows[int(0)].y;
    float _S2963 = 2.0f * - _S2954.rows[int(0)].x;
    float _S2964 = - _S2960 + _S2962;
    float _S2965 = _S2957 + - _S2961;
    float _S2966 = - _S2956 + _S2958;
    float _S2967 = _S2956 + _S2958;
    float _S2968 = _S2957 + _S2961;
    float _S2969 = _S2960 + _S2962;
    float _S2970 = (*dpquat_0).primal_0.w * (_S2959 + _S2963);
    float _S2971 = (*dpquat_0).primal_0.z * (_S2955 + _S2963);
    float _S2972 = (*dpquat_0).primal_0.y * (_S2955 + _S2959);
    float _S2973 = (*dpquat_0).primal_0.x * _S2964 + (*dpquat_0).primal_0.z * _S2967 + (*dpquat_0).primal_0.y * _S2968 + _S2970 + _S2970;
    float _S2974 = (*dpquat_0).primal_0.x * _S2965 + (*dpquat_0).primal_0.w * _S2967 + (*dpquat_0).primal_0.y * _S2969 + _S2971 + _S2971;
    float _S2975 = (*dpquat_0).primal_0.x * _S2966 + (*dpquat_0).primal_0.w * _S2968 + (*dpquat_0).primal_0.z * _S2969 + _S2972 + _S2972;
    float _S2976 = (*dpquat_0).primal_0.w * _S2964 + (*dpquat_0).primal_0.z * _S2965 + (*dpquat_0).primal_0.y * _S2966;
    float3  _S2977 = make_float3 (_S2952.differential_0.rows[int(0)].x, _S2952.differential_0.rows[int(1)].y, _S2952.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S2977;
    float4  _S2978 = make_float4 (0.0f);
    *&((&_S2978)->w) = _S2973;
    *&((&_S2978)->z) = _S2974;
    *&((&_S2978)->y) = _S2975;
    *&((&_S2978)->x) = _S2976;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S2978;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S2979, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2980, Matrix<float, 3, 3>  _S2981)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S2979, _S2980, _S2981);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_20, float3  scale_19, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S2982 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_20;
    (&dp_quat_0)->differential_0 = _S2982;
    float3  _S2983 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_19;
    (&dp_scale_0)->differential_0 = _S2983;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_16)
{
    float _S2984 = dOut_16.y;
    float _S2985 = dOut_16.z;
    float _S2986 = dOut_16.x;
    float _S2987 = (*a_0).primal_0.z * _S2984 + - (*a_0).primal_0.y * _S2985;
    float _S2988 = - (*a_0).primal_0.z * _S2986 + (*a_0).primal_0.x * _S2985;
    float _S2989 = (*a_0).primal_0.y * _S2986 + - (*a_0).primal_0.x * _S2984;
    float3  _S2990 = make_float3 (- (*b_0).primal_0.z * _S2984 + (*b_0).primal_0.y * _S2985, (*b_0).primal_0.z * _S2986 + - (*b_0).primal_0.x * _S2985, - (*b_0).primal_0.y * _S2986 + (*b_0).primal_0.x * _S2984);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S2990;
    float3  _S2991 = make_float3 (_S2987, _S2988, _S2989);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S2991;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S2992 = left_10.y;
    float _S2993 = right_10.z;
    float _S2994 = left_10.z;
    float _S2995 = right_10.y;
    float _S2996 = right_10.x;
    float _S2997 = left_10.x;
    return make_float3 (_S2992 * _S2993 - _S2994 * _S2995, _S2994 * _S2996 - _S2997 * _S2993, _S2997 * _S2995 - _S2992 * _S2996);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_1, float3  ray_d_1)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_1);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_1 - mean_16));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S2998, float3  _S2999)
{
    return cross_0(_S2998, _S2999);
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3000, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3001, float3  _S3002)
{
    _d_cross_0(_S3000, _S3001, _S3002);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_7)
{
    float3  _S3003 = (*dpray_o_2).primal_0 - (*dpmean_0).primal_0;
    float3  _S3004 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, _S3003);
    float3  _S3005 = s_primal_ctx_mul_1((*dpiscl_rot_0).primal_0, (*dpray_d_2).primal_0);
    float3  _S3006 = s_primal_ctx_cross_0(_S3005, _S3004);
    float _S3007 = -0.5f * s_primal_ctx_dot_0(_S3006, _S3006);
    float _S3008 = s_primal_ctx_dot_0(_S3005, _S3005);
    float _S3009 = _S3007 / _S3008;
    float _S3010 = _S3008 * _S3008;
    float _S3011 = (*dpopacity_0).primal_0 * _s_dOut_7;
    float _S3012 = s_primal_ctx_exp_1(_S3009) * _s_dOut_7;
    DiffPair_float_0 _S3013;
    (&_S3013)->primal_0 = _S3009;
    (&_S3013)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3013, _S3011);
    float _S3014 = _S3013.differential_0 / _S3010;
    float _S3015 = _S3007 * - _S3014;
    float _S3016 = _S3008 * _S3014;
    float3  _S3017 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3018;
    (&_S3018)->primal_0 = _S3005;
    (&_S3018)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3019;
    (&_S3019)->primal_0 = _S3005;
    (&_S3019)->differential_0 = _S3017;
    s_bwd_prop_dot_0(&_S3018, &_S3019, _S3015);
    float _S3020 = -0.5f * _S3016;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3021;
    (&_S3021)->primal_0 = _S3006;
    (&_S3021)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3022;
    (&_S3022)->primal_0 = _S3006;
    (&_S3022)->differential_0 = _S3017;
    s_bwd_prop_dot_0(&_S3021, &_S3022, _S3020);
    float3  _S3023 = _S3022.differential_0 + _S3021.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3024;
    (&_S3024)->primal_0 = _S3005;
    (&_S3024)->differential_0 = _S3017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3025;
    (&_S3025)->primal_0 = _S3004;
    (&_S3025)->differential_0 = _S3017;
    s_bwd_prop_cross_0(&_S3024, &_S3025, _S3023);
    float3  _S3026 = _S3019.differential_0 + _S3018.differential_0 + _S3024.differential_0;
    Matrix<float, 3, 3>  _S3027 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3028;
    (&_S3028)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3028)->differential_0 = _S3027;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3029;
    (&_S3029)->primal_0 = (*dpray_d_2).primal_0;
    (&_S3029)->differential_0 = _S3017;
    s_bwd_prop_mul_1(&_S3028, &_S3029, _S3026);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3030;
    (&_S3030)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S3030)->differential_0 = _S3027;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3031;
    (&_S3031)->primal_0 = _S3003;
    (&_S3031)->differential_0 = _S3017;
    s_bwd_prop_mul_1(&_S3030, &_S3031, _S3025.differential_0);
    float3  _S3032 = - _S3031.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S3029.differential_0;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S3031.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S3012;
    Matrix<float, 3, 3>  _S3033 = _S3028.differential_0 + _S3030.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S3033;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S3032;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3034, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3035, DiffPair_float_0 * _S3036, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3037, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3038, float _S3039)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S3034, _S3035, _S3036, _S3037, _S3038, _S3039);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_17, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_2, float3  ray_d_2, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S3040 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_17;
    (&dp_mean_0)->differential_0 = _S3040;
    Matrix<float, 3, 3>  _S3041 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S3041;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_2;
    (&dp_ray_o_0)->differential_0 = _S3040;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_2;
    (&dp_ray_d_0)->differential_0 = _S3040;
    s_bwd_evaluate_alpha_3dgs_0(&dp_mean_0, &dp_iscl_rot_0, &dp_opacity_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_mean_5 = dp_mean_0.differential_0;
    *v_iscl_rot_1 = dp_iscl_rot_0.differential_0;
    *v_opacity_5 = dp_opacity_0.differential_0;
    *v_ray_o_1 = dp_ray_o_0.differential_0;
    *v_ray_d_1 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_3dgs(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_2, float opacity_12, float3  rgb_10, float3  ray_o_3, float3  ray_d_3, float3  * out_rgb_0, float * depth_10)
{
    *out_rgb_0 = rgb_10;
    float3  _S3042 = mean_18 - ray_o_3;
    *depth_10 = 0.5f * (F32_log((dot_0(_S3042, _S3042) + 9.99999997475242708e-07f)));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S3043 = (*dpmean_1).primal_0 - (*dpray_o_3).primal_0;
    float _S3044 = 0.5f * dpdepth_0;
    DiffPair_float_0 _S3045;
    (&_S3045)->primal_0 = s_primal_ctx_dot_0(_S3043, _S3043) + 9.99999997475242708e-07f;
    (&_S3045)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3045, _S3044);
    float3  _S3046 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3047;
    (&_S3047)->primal_0 = _S3043;
    (&_S3047)->differential_0 = _S3046;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3048;
    (&_S3048)->primal_0 = _S3043;
    (&_S3048)->differential_0 = _S3046;
    s_bwd_prop_dot_0(&_S3047, &_S3048, _S3045.differential_0);
    float3  _S3049 = _S3048.differential_0 + _S3047.differential_0;
    float3  _S3050 = - _S3049;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S3046;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S3050;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S3051 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S3051;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S3049;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3052, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3053, DiffPair_float_0 * _S3054, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3055, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3056, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3057, float3  _S3058, float _S3059)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S3052, _S3053, _S3054, _S3055, _S3056, _S3057, _S3058, _S3059);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_19, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_4, float3  ray_d_4, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S3060 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_19;
    (&dp_mean_1)->differential_0 = _S3060;
    Matrix<float, 3, 3>  _S3061 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S3061;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S3060;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_4;
    (&dp_ray_o_1)->differential_0 = _S3060;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_4;
    (&dp_ray_d_1)->differential_0 = _S3060;
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
    float _S3062 = _slang_select(((*dpx_14).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_14).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_17;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S3062;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_15, float dOut_18)
{
    float _S3063 = (F32_exp2(((*dpx_15).primal_0))) * 50.693145751953125f * dOut_18;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S3063;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  dOut_19)
{
    float3  _S3064 = make_float3 (1.0f) / (*dpx_16).primal_0 * dOut_19;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S3064;
    return;
}

inline __device__ float3  log_0(float3  x_42)
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
        *_slang_vector_get_element_ptr(&result_14, i_10) = (F32_log((_slang_vector_get_element(x_42, i_10))));
        i_10 = i_10 + int(1);
    }
    return result_14;
}

inline __device__ void projection_opaque_triangle_persp(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_19, float3  t_18, float fx_24, float fy_24, float cx_19, float cy_19, FixedArray<float, 10>  * dist_coeffs_31, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_19, mean_20) + t_18;
        float _S3065 = mean_c_15.z;
        bool _S3066;
        if(_S3065 < near_plane_10)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = _S3065 > far_plane_10;
        }
        if(_S3066)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3067 = scale_20.x;
        float sx_1 = (F32_exp((_S3067)));
        float _S3068 = scale_20.y;
        float sy_1 = (F32_exp((_S3068)));
        float sz_1 = scale_20.z - 0.5f * (_S3067 + _S3068);
        float x_43 = quat_21.y;
        float x2_21 = x_43 * x_43;
        float y2_21 = quat_21.z * quat_21.z;
        float z2_36 = quat_21.w * quat_21.w;
        float xy_21 = quat_21.y * quat_21.z;
        float xz_21 = quat_21.y * quat_21.w;
        float yz_21 = quat_21.z * quat_21.w;
        float wx_21 = quat_21.x * quat_21.y;
        float wy_21 = quat_21.x * quat_21.z;
        float wz_21 = quat_21.x * quat_21.w;
        Matrix<float, 3, 3>  _S3069 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_36), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_36), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_0 = mul_0(R_19, mul_0(_S3069, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_18;
        float3  vert1_c_0 = mul_0(R_19, mul_0(_S3069, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_18;
        float3  vert2_c_0 = mul_0(R_19, mul_0(_S3069, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_18;
        float _S3070 = vert0_c_0.z;
        float _S3071 = vert1_c_0.z;
        float _S3072 = vert2_c_0.z;
        if(_S3070 < near_plane_10)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = _S3070 > far_plane_10;
        }
        if(_S3066)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = _S3071 < near_plane_10;
        }
        if(_S3066)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = _S3071 > far_plane_10;
        }
        if(_S3066)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = _S3072 < near_plane_10;
        }
        if(_S3066)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = _S3072 > far_plane_10;
        }
        if(_S3066)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        for(;;)
        {
            *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S3070);
            if(_S3070 < 0.0f)
            {
                _S3066 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3073 = camera_distortion_jac_0(*uv0_0, dist_coeffs_31);
                _S3066 = !((F32_min((determinant_0(_S3073)), ((F32_min((_S3073.rows[int(0)].x), (_S3073.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3066)
            {
                break;
            }
            float u_44 = (*uv0_0).x;
            float v_44 = (*uv0_0).y;
            float r2_44 = u_44 * u_44 + v_44 * v_44;
            float2  _S3074 = *uv0_0 * make_float2 (1.0f + r2_44 * ((*dist_coeffs_31)[int(0)] + r2_44 * ((*dist_coeffs_31)[int(1)] + r2_44 * ((*dist_coeffs_31)[int(2)] + r2_44 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_44 * v_44 + (*dist_coeffs_31)[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + (*dist_coeffs_31)[int(6)] * r2_44, 2.0f * (*dist_coeffs_31)[int(5)] * u_44 * v_44 + (*dist_coeffs_31)[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + (*dist_coeffs_31)[int(7)] * r2_44);
            float2  _S3075 = _S3074 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3074.x + (*dist_coeffs_31)[int(9)] * _S3074.y, 0.0f);
            *uv0_0 = make_float2 (fx_24 * _S3075.x + cx_19, fy_24 * _S3075.y + cy_19);
            break;
        }
        bool all_valid_8 = true & (!_S3066);
        for(;;)
        {
            *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S3071);
            if(_S3071 < 0.0f)
            {
                _S3066 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3076 = camera_distortion_jac_0(*uv1_0, dist_coeffs_31);
                _S3066 = !((F32_min((determinant_0(_S3076)), ((F32_min((_S3076.rows[int(0)].x), (_S3076.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3066)
            {
                break;
            }
            float u_45 = (*uv1_0).x;
            float v_45 = (*uv1_0).y;
            float r2_45 = u_45 * u_45 + v_45 * v_45;
            float2  _S3077 = *uv1_0 * make_float2 (1.0f + r2_45 * ((*dist_coeffs_31)[int(0)] + r2_45 * ((*dist_coeffs_31)[int(1)] + r2_45 * ((*dist_coeffs_31)[int(2)] + r2_45 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_45 * v_45 + (*dist_coeffs_31)[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + (*dist_coeffs_31)[int(6)] * r2_45, 2.0f * (*dist_coeffs_31)[int(5)] * u_45 * v_45 + (*dist_coeffs_31)[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + (*dist_coeffs_31)[int(7)] * r2_45);
            float2  _S3078 = _S3077 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3077.x + (*dist_coeffs_31)[int(9)] * _S3077.y, 0.0f);
            *uv1_0 = make_float2 (fx_24 * _S3078.x + cx_19, fy_24 * _S3078.y + cy_19);
            break;
        }
        bool all_valid_9 = all_valid_8 & (!_S3066);
        for(;;)
        {
            *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S3072);
            if(_S3072 < 0.0f)
            {
                _S3066 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3079 = camera_distortion_jac_0(*uv2_0, dist_coeffs_31);
                _S3066 = !((F32_min((determinant_0(_S3079)), ((F32_min((_S3079.rows[int(0)].x), (_S3079.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3066)
            {
                break;
            }
            float u_46 = (*uv2_0).x;
            float v_46 = (*uv2_0).y;
            float r2_46 = u_46 * u_46 + v_46 * v_46;
            float2  _S3080 = *uv2_0 * make_float2 (1.0f + r2_46 * ((*dist_coeffs_31)[int(0)] + r2_46 * ((*dist_coeffs_31)[int(1)] + r2_46 * ((*dist_coeffs_31)[int(2)] + r2_46 * (*dist_coeffs_31)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_31)[int(4)] * u_46 * v_46 + (*dist_coeffs_31)[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + (*dist_coeffs_31)[int(6)] * r2_46, 2.0f * (*dist_coeffs_31)[int(5)] * u_46 * v_46 + (*dist_coeffs_31)[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + (*dist_coeffs_31)[int(7)] * r2_46);
            float2  _S3081 = _S3080 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3080.x + (*dist_coeffs_31)[int(9)] * _S3080.y, 0.0f);
            *uv2_0 = make_float2 (fx_24 * _S3081.x + cx_19, fy_24 * _S3081.y + cy_19);
            break;
        }
        if(!(all_valid_9 & (!_S3066)))
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
            _S3066 = true;
        }
        else
        {
            _S3066 = xmin_5 >= float(image_width_15);
        }
        if(_S3066)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = ymax_5 <= 0.0f;
        }
        if(_S3066)
        {
            _S3066 = true;
        }
        else
        {
            _S3066 = ymin_5 >= float(image_height_15);
        }
        if(_S3066)
        {
            _S3066 = true;
        }
        else
        {
            if(_S3065 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S3066 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S3066 = false;
                }
                if(_S3066)
                {
                    _S3066 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S3066 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S3066 = false;
                    }
                }
            }
            else
            {
                _S3066 = false;
            }
        }
        if(_S3066)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S3082 = mean_20 - - mul_0(transpose_0(R_19), t_18);
        float _S3083 = _S3082.x;
        float _S3084 = _S3082.y;
        float _S3085 = _S3082.z;
        float norm_10 = (F32_sqrt((_S3083 * _S3083 + _S3084 * _S3084 + _S3085 * _S3085)));
        float x_44 = _S3083 / norm_10;
        float y_18 = _S3084 / norm_10;
        float z_15 = _S3085 / norm_10;
        float z2_37 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_44 * x_44 - y_18 * y_18;
        float fS1_15 = 2.0f * x_44 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_37 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_44) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_37 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_44) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_44 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_37 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_44) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_44 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S3086 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S3086);
        float3  _S3087 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S3088 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S3087 + _S3088 + make_float3 (0.5f), _S3086);
        (*rgb_12)[int(2)] = max_0(_S3087 - _S3088 + make_float3 (0.5f), _S3086);
        float3  _S3089 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S3089 * make_float3 (float(- (F32_sign((dot_0(_S3089, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_20, float3  t_19, float fx_25, float fy_25, float cx_20, float cy_20, FixedArray<float, 10>  * dist_coeffs_32, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    bool _S3090;
    bool _S3091;
    bool _S3092;
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_20, mean_21) + t_19;
        float _S3093 = length_1(mean_c_16);
        bool _S3094;
        if(_S3093 < near_plane_11)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = _S3093 > far_plane_11;
        }
        if(_S3094)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3095 = scale_21.x;
        float sx_2 = (F32_exp((_S3095)));
        float _S3096 = scale_21.y;
        float sy_2 = (F32_exp((_S3096)));
        float sz_2 = scale_21.z - 0.5f * (_S3095 + _S3096);
        float x_45 = quat_22.y;
        float x2_22 = x_45 * x_45;
        float y2_22 = quat_22.z * quat_22.z;
        float z2_38 = quat_22.w * quat_22.w;
        float xy_22 = quat_22.y * quat_22.z;
        float xz_22 = quat_22.y * quat_22.w;
        float yz_22 = quat_22.z * quat_22.w;
        float wx_22 = quat_22.x * quat_22.y;
        float wy_22 = quat_22.x * quat_22.z;
        float wz_22 = quat_22.x * quat_22.w;
        Matrix<float, 3, 3>  _S3097 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_38), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_38), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
        float3  vert0_c_1 = mul_0(R_20, mul_0(_S3097, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_19;
        float3  vert1_c_1 = mul_0(R_20, mul_0(_S3097, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_19;
        float3  vert2_c_1 = mul_0(R_20, mul_0(_S3097, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_19;
        float _S3098 = length_1(vert0_c_1);
        float _S3099 = length_1(vert1_c_1);
        float _S3100 = length_1(vert2_c_1);
        if(_S3098 < near_plane_11)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = _S3098 > far_plane_11;
        }
        if(_S3094)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = _S3099 < near_plane_11;
        }
        if(_S3094)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = _S3099 > far_plane_11;
        }
        if(_S3094)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = _S3100 < near_plane_11;
        }
        if(_S3094)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = _S3100 > far_plane_11;
        }
        if(_S3094)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float k_7;
        for(;;)
        {
            float2  _S3101 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_24 = length_0(_S3101);
            float _S3102 = vert0_c_1.z;
            float theta_20 = (F32_atan2((r_24), (_S3102)));
            if(theta_20 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_20 * theta_20 / 3.0f) / _S3102;
            }
            else
            {
                k_7 = theta_20 / r_24;
            }
            float2  _S3103 = _S3101 * make_float2 (k_7);
            *uv0_1 = _S3103;
            Matrix<float, 2, 2>  _S3104 = camera_distortion_jac_0(_S3103, dist_coeffs_32);
            bool _S3105 = !((F32_min((determinant_0(_S3104)), ((F32_min((_S3104.rows[int(0)].x), (_S3104.rows[int(1)].y)))))) > 0.0f);
            _S3090 = _S3105;
            if(_S3105)
            {
                break;
            }
            float u_47 = (*uv0_1).x;
            float v_47 = (*uv0_1).y;
            float r2_47 = u_47 * u_47 + v_47 * v_47;
            float2  _S3106 = *uv0_1 * make_float2 (1.0f + r2_47 * ((*dist_coeffs_32)[int(0)] + r2_47 * ((*dist_coeffs_32)[int(1)] + r2_47 * ((*dist_coeffs_32)[int(2)] + r2_47 * (*dist_coeffs_32)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_32)[int(4)] * u_47 * v_47 + (*dist_coeffs_32)[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + (*dist_coeffs_32)[int(6)] * r2_47, 2.0f * (*dist_coeffs_32)[int(5)] * u_47 * v_47 + (*dist_coeffs_32)[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + (*dist_coeffs_32)[int(7)] * r2_47);
            float2  _S3107 = _S3106 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3106.x + (*dist_coeffs_32)[int(9)] * _S3106.y, 0.0f);
            *uv0_1 = make_float2 (fx_25 * _S3107.x + cx_20, fy_25 * _S3107.y + cy_20);
            break;
        }
        bool all_valid_10 = true & (!_S3090);
        for(;;)
        {
            float2  _S3108 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_25 = length_0(_S3108);
            float _S3109 = vert1_c_1.z;
            float theta_21 = (F32_atan2((r_25), (_S3109)));
            if(theta_21 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_21 * theta_21 / 3.0f) / _S3109;
            }
            else
            {
                k_7 = theta_21 / r_25;
            }
            float2  _S3110 = _S3108 * make_float2 (k_7);
            *uv1_1 = _S3110;
            Matrix<float, 2, 2>  _S3111 = camera_distortion_jac_0(_S3110, dist_coeffs_32);
            bool _S3112 = !((F32_min((determinant_0(_S3111)), ((F32_min((_S3111.rows[int(0)].x), (_S3111.rows[int(1)].y)))))) > 0.0f);
            _S3091 = _S3112;
            if(_S3112)
            {
                break;
            }
            float u_48 = (*uv1_1).x;
            float v_48 = (*uv1_1).y;
            float r2_48 = u_48 * u_48 + v_48 * v_48;
            float2  _S3113 = *uv1_1 * make_float2 (1.0f + r2_48 * ((*dist_coeffs_32)[int(0)] + r2_48 * ((*dist_coeffs_32)[int(1)] + r2_48 * ((*dist_coeffs_32)[int(2)] + r2_48 * (*dist_coeffs_32)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_32)[int(4)] * u_48 * v_48 + (*dist_coeffs_32)[int(5)] * (r2_48 + 2.0f * u_48 * u_48) + (*dist_coeffs_32)[int(6)] * r2_48, 2.0f * (*dist_coeffs_32)[int(5)] * u_48 * v_48 + (*dist_coeffs_32)[int(4)] * (r2_48 + 2.0f * v_48 * v_48) + (*dist_coeffs_32)[int(7)] * r2_48);
            float2  _S3114 = _S3113 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3113.x + (*dist_coeffs_32)[int(9)] * _S3113.y, 0.0f);
            *uv1_1 = make_float2 (fx_25 * _S3114.x + cx_20, fy_25 * _S3114.y + cy_20);
            break;
        }
        bool all_valid_11 = all_valid_10 & (!_S3091);
        for(;;)
        {
            float2  _S3115 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_26 = length_0(_S3115);
            float _S3116 = vert2_c_1.z;
            float theta_22 = (F32_atan2((r_26), (_S3116)));
            if(theta_22 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_22 * theta_22 / 3.0f) / _S3116;
            }
            else
            {
                k_7 = theta_22 / r_26;
            }
            float2  _S3117 = _S3115 * make_float2 (k_7);
            *uv2_1 = _S3117;
            Matrix<float, 2, 2>  _S3118 = camera_distortion_jac_0(_S3117, dist_coeffs_32);
            bool _S3119 = !((F32_min((determinant_0(_S3118)), ((F32_min((_S3118.rows[int(0)].x), (_S3118.rows[int(1)].y)))))) > 0.0f);
            _S3092 = _S3119;
            if(_S3119)
            {
                break;
            }
            float u_49 = (*uv2_1).x;
            float v_49 = (*uv2_1).y;
            float r2_49 = u_49 * u_49 + v_49 * v_49;
            float2  _S3120 = *uv2_1 * make_float2 (1.0f + r2_49 * ((*dist_coeffs_32)[int(0)] + r2_49 * ((*dist_coeffs_32)[int(1)] + r2_49 * ((*dist_coeffs_32)[int(2)] + r2_49 * (*dist_coeffs_32)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_32)[int(4)] * u_49 * v_49 + (*dist_coeffs_32)[int(5)] * (r2_49 + 2.0f * u_49 * u_49) + (*dist_coeffs_32)[int(6)] * r2_49, 2.0f * (*dist_coeffs_32)[int(5)] * u_49 * v_49 + (*dist_coeffs_32)[int(4)] * (r2_49 + 2.0f * v_49 * v_49) + (*dist_coeffs_32)[int(7)] * r2_49);
            float2  _S3121 = _S3120 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3120.x + (*dist_coeffs_32)[int(9)] * _S3120.y, 0.0f);
            *uv2_1 = make_float2 (fx_25 * _S3121.x + cx_20, fy_25 * _S3121.y + cy_20);
            break;
        }
        if(!(all_valid_11 & (!_S3092)))
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
            _S3094 = true;
        }
        else
        {
            _S3094 = xmin_6 >= float(image_width_16);
        }
        if(_S3094)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = ymax_6 <= 0.0f;
        }
        if(_S3094)
        {
            _S3094 = true;
        }
        else
        {
            _S3094 = ymin_6 >= float(image_height_16);
        }
        if(_S3094)
        {
            _S3094 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S3094 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S3094 = false;
                }
                if(_S3094)
                {
                    _S3094 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S3094 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S3094 = false;
                    }
                }
            }
            else
            {
                _S3094 = false;
            }
        }
        if(_S3094)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S3098, _S3099, _S3100) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S3122 = mean_21 - - mul_0(transpose_0(R_20), t_19);
        float _S3123 = _S3122.x;
        float _S3124 = _S3122.y;
        float _S3125 = _S3122.z;
        float norm_11 = (F32_sqrt((_S3123 * _S3123 + _S3124 * _S3124 + _S3125 * _S3125)));
        float x_46 = _S3123 / norm_11;
        float y_19 = _S3124 / norm_11;
        float z_16 = _S3125 / norm_11;
        float z2_39 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_46 * x_46 - y_19 * y_19;
        float fS1_16 = 2.0f * x_46 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_39 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_46) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_39 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_46) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_39 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_46) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S3126 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S3126);
        float3  _S3127 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S3128 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S3127 + _S3128 + make_float3 (0.5f), _S3126);
        (*rgb_13)[int(2)] = max_0(_S3127 - _S3128 + make_float3 (0.5f), _S3126);
        float3  _S3129 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S3129 * make_float3 (float(- (F32_sign((dot_0(_S3129, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_21, float3  t_20, float fx_26, float fy_26, float cx_21, float cy_21, FixedArray<float, 10>  * dist_coeffs_33, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_21, mean_22) + t_20;
    float _S3130 = scale_22.x;
    float sx_3 = (F32_exp((_S3130)));
    float _S3131 = scale_22.y;
    float sy_3 = (F32_exp((_S3131)));
    float sz_3 = scale_22.z - 0.5f * (_S3130 + _S3131);
    float x_47 = quat_23.y;
    float x2_23 = x_47 * x_47;
    float y2_23 = quat_23.z * quat_23.z;
    float z2_40 = quat_23.w * quat_23.w;
    float xy_23 = quat_23.y * quat_23.z;
    float xz_23 = quat_23.y * quat_23.w;
    float yz_23 = quat_23.z * quat_23.w;
    float wx_23 = quat_23.x * quat_23.y;
    float wy_23 = quat_23.x * quat_23.z;
    float wz_23 = quat_23.x * quat_23.w;
    Matrix<float, 3, 3>  _S3132 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_40), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_40), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_2 = mul_0(R_21, mul_0(_S3132, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_20;
    float3  vert1_c_2 = mul_0(R_21, mul_0(_S3132, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_20;
    float3  vert2_c_2 = mul_0(R_21, mul_0(_S3132, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_20;
    float2  _S3133 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_50 = _S3133.x;
    float v_50 = _S3133.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S3134 = 2.0f * (*dist_coeffs_33)[int(4)];
    float _S3135 = 2.0f * (*dist_coeffs_33)[int(5)];
    float2  _S3136 = _S3133 * make_float2 (1.0f + r2_50 * ((*dist_coeffs_33)[int(0)] + r2_50 * ((*dist_coeffs_33)[int(1)] + r2_50 * ((*dist_coeffs_33)[int(2)] + r2_50 * (*dist_coeffs_33)[int(3)])))) + make_float2 (_S3134 * u_50 * v_50 + (*dist_coeffs_33)[int(5)] * (r2_50 + 2.0f * u_50 * u_50) + (*dist_coeffs_33)[int(6)] * r2_50, _S3135 * u_50 * v_50 + (*dist_coeffs_33)[int(4)] * (r2_50 + 2.0f * v_50 * v_50) + (*dist_coeffs_33)[int(7)] * r2_50);
    float2  _S3137 = _S3136 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3136.x + (*dist_coeffs_33)[int(9)] * _S3136.y, 0.0f);
    *uv0_2 = make_float2 (fx_26 * _S3137.x + cx_21, fy_26 * _S3137.y + cy_21);
    float2  _S3138 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_51 = _S3138.x;
    float v_51 = _S3138.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float2  _S3139 = _S3138 * make_float2 (1.0f + r2_51 * ((*dist_coeffs_33)[int(0)] + r2_51 * ((*dist_coeffs_33)[int(1)] + r2_51 * ((*dist_coeffs_33)[int(2)] + r2_51 * (*dist_coeffs_33)[int(3)])))) + make_float2 (_S3134 * u_51 * v_51 + (*dist_coeffs_33)[int(5)] * (r2_51 + 2.0f * u_51 * u_51) + (*dist_coeffs_33)[int(6)] * r2_51, _S3135 * u_51 * v_51 + (*dist_coeffs_33)[int(4)] * (r2_51 + 2.0f * v_51 * v_51) + (*dist_coeffs_33)[int(7)] * r2_51);
    float2  _S3140 = _S3139 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3139.x + (*dist_coeffs_33)[int(9)] * _S3139.y, 0.0f);
    *uv1_2 = make_float2 (fx_26 * _S3140.x + cx_21, fy_26 * _S3140.y + cy_21);
    float2  _S3141 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_52 = _S3141.x;
    float v_52 = _S3141.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float2  _S3142 = _S3141 * make_float2 (1.0f + r2_52 * ((*dist_coeffs_33)[int(0)] + r2_52 * ((*dist_coeffs_33)[int(1)] + r2_52 * ((*dist_coeffs_33)[int(2)] + r2_52 * (*dist_coeffs_33)[int(3)])))) + make_float2 (_S3134 * u_52 * v_52 + (*dist_coeffs_33)[int(5)] * (r2_52 + 2.0f * u_52 * u_52) + (*dist_coeffs_33)[int(6)] * r2_52, _S3135 * u_52 * v_52 + (*dist_coeffs_33)[int(4)] * (r2_52 + 2.0f * v_52 * v_52) + (*dist_coeffs_33)[int(7)] * r2_52);
    float2  _S3143 = _S3142 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3142.x + (*dist_coeffs_33)[int(9)] * _S3142.y, 0.0f);
    float _S3144 = fx_26 * _S3143.x + cx_21;
    float _S3145 = fy_26 * _S3143.y + cy_21;
    float2  _S3146 = make_float2 (_S3144, _S3145);
    *uv2_2 = _S3146;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S3146 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S3146)));
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S3144))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S3145))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S3144))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S3145))) + offset_2)))));
    *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_2 = hardness_2;
    float3  _S3147 = mean_22 - - mul_0(transpose_0(R_21), t_20);
    float _S3148 = _S3147.x;
    float _S3149 = _S3147.y;
    float _S3150 = _S3147.z;
    float norm_12 = (F32_sqrt((_S3148 * _S3148 + _S3149 * _S3149 + _S3150 * _S3150)));
    float x_48 = _S3148 / norm_12;
    float y_20 = _S3149 / norm_12;
    float z_17 = _S3150 / norm_12;
    float z2_41 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_48 * x_48 - y_20 * y_20;
    float fS1_17 = 2.0f * x_48 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_41 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_48) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_41 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_48) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_48 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_41 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_48) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_48 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S3151 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S3151);
    float3  _S3152 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S3153 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S3152 + _S3153 + make_float3 (0.5f), _S3151);
    (*rgb_14)[int(2)] = max_0(_S3152 - _S3153 + make_float3 (0.5f), _S3151);
    float3  _S3154 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S3154 * make_float3 (float(- (F32_sign((dot_0(_S3154, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_22, float3  t_21, float fx_27, float fy_27, float cx_22, float cy_22, FixedArray<float, 10>  * dist_coeffs_34, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_22, mean_23) + t_21;
    float _S3155 = scale_23.x;
    float sx_4 = (F32_exp((_S3155)));
    float _S3156 = scale_23.y;
    float sy_4 = (F32_exp((_S3156)));
    float sz_4 = scale_23.z - 0.5f * (_S3155 + _S3156);
    float x_49 = quat_24.y;
    float x2_24 = x_49 * x_49;
    float y2_24 = quat_24.z * quat_24.z;
    float z2_42 = quat_24.w * quat_24.w;
    float xy_24 = quat_24.y * quat_24.z;
    float xz_24 = quat_24.y * quat_24.w;
    float yz_24 = quat_24.z * quat_24.w;
    float wx_24 = quat_24.x * quat_24.y;
    float wy_24 = quat_24.x * quat_24.z;
    float wz_24 = quat_24.x * quat_24.w;
    Matrix<float, 3, 3>  _S3157 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_42), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_42), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  vert0_c_3 = mul_0(R_22, mul_0(_S3157, make_float3 (sx_4, 0.0f, 0.0f)) + mean_23) + t_21;
    float3  vert1_c_3 = mul_0(R_22, mul_0(_S3157, make_float3 (sx_4 * (-0.5f + sz_4), sy_4, 0.0f)) + mean_23) + t_21;
    float3  vert2_c_3 = mul_0(R_22, mul_0(_S3157, make_float3 (sx_4 * (-0.5f - sz_4), - sy_4, 0.0f)) + mean_23) + t_21;
    float2  _S3158 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_27 = length_0(_S3158);
    float _S3159 = vert0_c_3.z;
    float theta_23 = (F32_atan2((r_27), (_S3159)));
    float k_8;
    if(theta_23 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_23 * theta_23 / 3.0f) / _S3159;
    }
    else
    {
        k_8 = theta_23 / r_27;
    }
    float2  _S3160 = _S3158 * make_float2 (k_8);
    float u_53 = _S3160.x;
    float v_53 = _S3160.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S3161 = 2.0f * (*dist_coeffs_34)[int(4)];
    float _S3162 = 2.0f * (*dist_coeffs_34)[int(5)];
    float2  _S3163 = _S3160 * make_float2 (1.0f + r2_53 * ((*dist_coeffs_34)[int(0)] + r2_53 * ((*dist_coeffs_34)[int(1)] + r2_53 * ((*dist_coeffs_34)[int(2)] + r2_53 * (*dist_coeffs_34)[int(3)])))) + make_float2 (_S3161 * u_53 * v_53 + (*dist_coeffs_34)[int(5)] * (r2_53 + 2.0f * u_53 * u_53) + (*dist_coeffs_34)[int(6)] * r2_53, _S3162 * u_53 * v_53 + (*dist_coeffs_34)[int(4)] * (r2_53 + 2.0f * v_53 * v_53) + (*dist_coeffs_34)[int(7)] * r2_53);
    float2  _S3164 = _S3163 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3163.x + (*dist_coeffs_34)[int(9)] * _S3163.y, 0.0f);
    *uv0_3 = make_float2 (fx_27 * _S3164.x + cx_22, fy_27 * _S3164.y + cy_22);
    float2  _S3165 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_28 = length_0(_S3165);
    float _S3166 = vert1_c_3.z;
    float theta_24 = (F32_atan2((r_28), (_S3166)));
    if(theta_24 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_24 * theta_24 / 3.0f) / _S3166;
    }
    else
    {
        k_8 = theta_24 / r_28;
    }
    float2  _S3167 = _S3165 * make_float2 (k_8);
    float u_54 = _S3167.x;
    float v_54 = _S3167.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float2  _S3168 = _S3167 * make_float2 (1.0f + r2_54 * ((*dist_coeffs_34)[int(0)] + r2_54 * ((*dist_coeffs_34)[int(1)] + r2_54 * ((*dist_coeffs_34)[int(2)] + r2_54 * (*dist_coeffs_34)[int(3)])))) + make_float2 (_S3161 * u_54 * v_54 + (*dist_coeffs_34)[int(5)] * (r2_54 + 2.0f * u_54 * u_54) + (*dist_coeffs_34)[int(6)] * r2_54, _S3162 * u_54 * v_54 + (*dist_coeffs_34)[int(4)] * (r2_54 + 2.0f * v_54 * v_54) + (*dist_coeffs_34)[int(7)] * r2_54);
    float2  _S3169 = _S3168 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3168.x + (*dist_coeffs_34)[int(9)] * _S3168.y, 0.0f);
    *uv1_3 = make_float2 (fx_27 * _S3169.x + cx_22, fy_27 * _S3169.y + cy_22);
    float2  _S3170 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_29 = length_0(_S3170);
    float _S3171 = vert2_c_3.z;
    float theta_25 = (F32_atan2((r_29), (_S3171)));
    if(theta_25 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_25 * theta_25 / 3.0f) / _S3171;
    }
    else
    {
        k_8 = theta_25 / r_29;
    }
    float2  _S3172 = _S3170 * make_float2 (k_8);
    float u_55 = _S3172.x;
    float v_55 = _S3172.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float2  _S3173 = _S3172 * make_float2 (1.0f + r2_55 * ((*dist_coeffs_34)[int(0)] + r2_55 * ((*dist_coeffs_34)[int(1)] + r2_55 * ((*dist_coeffs_34)[int(2)] + r2_55 * (*dist_coeffs_34)[int(3)])))) + make_float2 (_S3161 * u_55 * v_55 + (*dist_coeffs_34)[int(5)] * (r2_55 + 2.0f * u_55 * u_55) + (*dist_coeffs_34)[int(6)] * r2_55, _S3162 * u_55 * v_55 + (*dist_coeffs_34)[int(4)] * (r2_55 + 2.0f * v_55 * v_55) + (*dist_coeffs_34)[int(7)] * r2_55);
    float2  _S3174 = _S3173 + make_float2 ((*dist_coeffs_34)[int(8)] * _S3173.x + (*dist_coeffs_34)[int(9)] * _S3173.y, 0.0f);
    float _S3175 = fx_27 * _S3174.x + cx_22;
    float _S3176 = fy_27 * _S3174.y + cy_22;
    float2  _S3177 = make_float2 (_S3175, _S3176);
    *uv2_3 = _S3177;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S3177 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S3177)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S3175))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S3176))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S3175))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S3176))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S3178 = mean_23 - - mul_0(transpose_0(R_22), t_21);
    float _S3179 = _S3178.x;
    float _S3180 = _S3178.y;
    float _S3181 = _S3178.z;
    float norm_13 = (F32_sqrt((_S3179 * _S3179 + _S3180 * _S3180 + _S3181 * _S3181)));
    float x_50 = _S3179 / norm_13;
    float y_21 = _S3180 / norm_13;
    float z_18 = _S3181 / norm_13;
    float z2_43 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_50 * x_50 - y_21 * y_21;
    float fS1_18 = 2.0f * x_50 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_43 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_50) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_43 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_50) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_50 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_43 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_50) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_50 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S3182 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S3182);
    float3  _S3183 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S3184 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S3183 + _S3184 + make_float3 (0.5f), _S3182);
    (*rgb_15)[int(2)] = max_0(_S3183 - _S3184 + make_float3 (0.5f), _S3182);
    float3  _S3185 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S3185 * make_float3 (float(- (F32_sign((dot_0(_S3185, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3186, float3  _S3187)
{
    _d_log_vector_0(_S3186, _S3187);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S3188, float _S3189)
{
    _d_exp2_0(_S3188, _S3189);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S3190, float _S3191)
{
    _d_abs_0(_S3190, _S3191);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_23, float3  t_22, float fx_28, float fy_28, float cx_23, float cy_23, FixedArray<float, 10>  * dist_coeffs_35, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_1(R_23, mean_24) + t_22;
    float _S3192 = scale_24.x;
    float _S3193 = s_primal_ctx_exp_1(_S3192);
    float _S3194 = scale_24.y;
    float _S3195 = s_primal_ctx_exp_1(_S3194);
    float sz_5 = scale_24.z - 0.5f * (_S3192 + _S3194);
    float _S3196 = quat_25.y;
    float x2_25 = _S3196 * _S3196;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_44 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S3197 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_44), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_44), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S3198 = make_float3 (_S3193, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_1(_S3197, _S3198) + mean_24;
    float _S3199 = -0.5f + sz_5;
    float3  _S3200 = make_float3 (_S3193 * _S3199, _S3195, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_1(_S3197, _S3200) + mean_24;
    float _S3201 = -0.5f - sz_5;
    float3  _S3202 = make_float3 (_S3193 * _S3201, - _S3195, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_1(_S3197, _S3202) + mean_24;
    float3  vert0_c_4 = s_primal_ctx_mul_1(R_23, vert0_1) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_1(R_23, vert1_1) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_1(R_23, vert2_1) + t_22;
    float2  _S3203 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S3204 = vert0_c_4.z;
    float2  _S3205 = make_float2 (_S3204);
    float2  _S3206 = _S3203 / make_float2 (_S3204);
    float2  _S3207 = make_float2 (_S3204 * _S3204);
    float u_56 = _S3206.x;
    float v_56 = _S3206.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S3208 = (*dist_coeffs_35)[int(2)] + r2_56 * (*dist_coeffs_35)[int(3)];
    float _S3209 = (*dist_coeffs_35)[int(1)] + r2_56 * _S3208;
    float _S3210 = (*dist_coeffs_35)[int(0)] + r2_56 * _S3209;
    float radial_1 = 1.0f + r2_56 * _S3210;
    float2  _S3211 = make_float2 (radial_1);
    float _S3212 = 2.0f * (*dist_coeffs_35)[int(4)];
    float _S3213 = _S3212 * u_56;
    float _S3214 = 2.0f * u_56;
    float _S3215 = 2.0f * (*dist_coeffs_35)[int(5)];
    float _S3216 = _S3215 * u_56;
    float _S3217 = 2.0f * v_56;
    float2  _S3218 = _S3206 * make_float2 (radial_1) + make_float2 (_S3213 * v_56 + (*dist_coeffs_35)[int(5)] * (r2_56 + _S3214 * u_56) + (*dist_coeffs_35)[int(6)] * r2_56, _S3216 * v_56 + (*dist_coeffs_35)[int(4)] * (r2_56 + _S3217 * v_56) + (*dist_coeffs_35)[int(7)] * r2_56);
    float2  _S3219 = _S3218 + make_float2 ((*dist_coeffs_35)[int(8)] * _S3218.x + (*dist_coeffs_35)[int(9)] * _S3218.y, 0.0f);
    float _S3220 = fx_28 * _S3219.x + cx_23;
    float _S3221 = fy_28 * _S3219.y + cy_23;
    float2  _S3222 = make_float2 (_S3220, _S3221);
    float2  _S3223 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S3224 = vert1_c_4.z;
    float2  _S3225 = make_float2 (_S3224);
    float2  _S3226 = _S3223 / make_float2 (_S3224);
    float2  _S3227 = make_float2 (_S3224 * _S3224);
    float u_57 = _S3226.x;
    float v_57 = _S3226.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S3228 = (*dist_coeffs_35)[int(2)] + r2_57 * (*dist_coeffs_35)[int(3)];
    float _S3229 = (*dist_coeffs_35)[int(1)] + r2_57 * _S3228;
    float _S3230 = (*dist_coeffs_35)[int(0)] + r2_57 * _S3229;
    float radial_2 = 1.0f + r2_57 * _S3230;
    float2  _S3231 = make_float2 (radial_2);
    float _S3232 = _S3212 * u_57;
    float _S3233 = 2.0f * u_57;
    float _S3234 = _S3215 * u_57;
    float _S3235 = 2.0f * v_57;
    float2  _S3236 = _S3226 * make_float2 (radial_2) + make_float2 (_S3232 * v_57 + (*dist_coeffs_35)[int(5)] * (r2_57 + _S3233 * u_57) + (*dist_coeffs_35)[int(6)] * r2_57, _S3234 * v_57 + (*dist_coeffs_35)[int(4)] * (r2_57 + _S3235 * v_57) + (*dist_coeffs_35)[int(7)] * r2_57);
    float2  _S3237 = _S3236 + make_float2 ((*dist_coeffs_35)[int(8)] * _S3236.x + (*dist_coeffs_35)[int(9)] * _S3236.y, 0.0f);
    float _S3238 = fx_28 * _S3237.x + cx_23;
    float _S3239 = fy_28 * _S3237.y + cy_23;
    float2  _S3240 = make_float2 (_S3238, _S3239);
    float2  _S3241 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S3242 = vert2_c_4.z;
    float2  _S3243 = make_float2 (_S3242);
    float2  _S3244 = _S3241 / make_float2 (_S3242);
    float2  _S3245 = make_float2 (_S3242 * _S3242);
    float u_58 = _S3244.x;
    float v_58 = _S3244.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S3246 = (*dist_coeffs_35)[int(2)] + r2_58 * (*dist_coeffs_35)[int(3)];
    float _S3247 = (*dist_coeffs_35)[int(1)] + r2_58 * _S3246;
    float _S3248 = (*dist_coeffs_35)[int(0)] + r2_58 * _S3247;
    float radial_3 = 1.0f + r2_58 * _S3248;
    float2  _S3249 = make_float2 (radial_3);
    float _S3250 = _S3212 * u_58;
    float _S3251 = 2.0f * u_58;
    float _S3252 = _S3215 * u_58;
    float _S3253 = 2.0f * v_58;
    float2  _S3254 = _S3244 * make_float2 (radial_3) + make_float2 (_S3250 * v_58 + (*dist_coeffs_35)[int(5)] * (r2_58 + _S3251 * u_58) + (*dist_coeffs_35)[int(6)] * r2_58, _S3252 * v_58 + (*dist_coeffs_35)[int(4)] * (r2_58 + _S3253 * v_58) + (*dist_coeffs_35)[int(7)] * r2_58);
    float2  _S3255 = _S3254 + make_float2 ((*dist_coeffs_35)[int(8)] * _S3254.x + (*dist_coeffs_35)[int(9)] * _S3254.y, 0.0f);
    float _S3256 = fx_28 * _S3255.x + cx_23;
    float _S3257 = fy_28 * _S3255.y + cy_23;
    float2  _S3258 = make_float2 (_S3256, _S3257);
    float2  e0_4 = _S3240 - _S3222;
    float2  e1_4 = _S3258 - _S3240;
    float2  e2_0 = _S3222 - _S3258;
    float _S3259 = e0_4.x;
    float _S3260 = e1_4.y;
    float _S3261 = e0_4.y;
    float _S3262 = e1_4.x;
    float _S3263 = _S3259 * _S3260 - _S3261 * _S3262;
    float _S3264 = 1.0f - hardness_4.y;
    float _S3265 = -1.0f / _S3264;
    float _S3266 = _S3264 * _S3264;
    float _S3267 = s_primal_ctx_max_0(_S3220, _S3238);
    float _S3268 = s_primal_ctx_min_0(_S3220, _S3238);
    float _S3269 = s_primal_ctx_max_0(_S3221, _S3239);
    float _S3270 = s_primal_ctx_min_0(_S3221, _S3239);
    float3  _S3271 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3272 = transpose_0(R_23);
    float3  _S3273 = mean_24 - - s_primal_ctx_mul_1(_S3272, t_22);
    float _S3274 = _S3273.x;
    float _S3275 = _S3273.y;
    float _S3276 = _S3273.z;
    float _S3277 = _S3274 * _S3274 + _S3275 * _S3275 + _S3276 * _S3276;
    float _S3278 = s_primal_ctx_sqrt_0(_S3277);
    float x_51 = _S3274 / _S3278;
    float3  _S3279 = make_float3 (x_51);
    float _S3280 = _S3278 * _S3278;
    float y_22 = _S3275 / _S3278;
    float z_19 = _S3276 / _S3278;
    float3  _S3281 = make_float3 (z_19);
    float _S3282 = - y_22;
    float3  _S3283 = make_float3 (_S3282);
    float z2_45 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_51 * x_51 - y_22 * y_22;
    float _S3284 = 2.0f * x_51;
    float fS1_19 = _S3284 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_45 - 0.31539157032966614f;
    float3  _S3285 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_51;
    float3  _S3286 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S3287 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S3288 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S3289 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_45 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S3290 = 1.86588168144226074f * z2_45 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S3290;
    float3  _S3291 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_51;
    float3  _S3292 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S3293 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S3294 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S3295 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_51 * fC1_19 - y_22 * fS1_19);
    float3  _S3296 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_51 * fS1_19 + y_22 * fC1_19);
    float3  _S3297 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3282) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_51) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S3298 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S3299 = make_float3 (0.0f);
    float3  _S3300 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S3301 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3302 = make_float3 (_S3301);
    float3  _S3303 = (*ch_coeffs_4)[int(1)] * make_float3 (_S3301);
    float3  _S3304 = _S3300 + _S3303 + make_float3 (0.5f);
    float3  _S3305 = _S3300 - _S3303 + make_float3 (0.5f);
    float3  _S3306 = vert1_c_4 - vert0_c_4;
    float3  _S3307 = vert2_c_4 - vert0_c_4;
    float3  _S3308 = s_primal_ctx_cross_0(_S3306, _S3307);
    float3  _S3309 = normalize_0(_S3308);
    float3  _S3310 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3309, mean_c_19)))))) * v_normal_0;
    float3  _S3311 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3312;
    (&_S3312)->primal_0 = _S3309;
    (&_S3312)->differential_0 = _S3311;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3313;
    (&_S3313)->primal_0 = mean_c_19;
    (&_S3313)->differential_0 = _S3311;
    s_bwd_prop_dot_0(&_S3312, &_S3313, 0.0f);
    float3  _S3314 = _S3310 + _S3312.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3315;
    (&_S3315)->primal_0 = _S3308;
    (&_S3315)->differential_0 = _S3311;
    s_bwd_normalize_impl_0(&_S3315, _S3314);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3316;
    (&_S3316)->primal_0 = _S3306;
    (&_S3316)->differential_0 = _S3311;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3317;
    (&_S3317)->primal_0 = _S3307;
    (&_S3317)->differential_0 = _S3311;
    s_bwd_prop_cross_0(&_S3316, &_S3317, _S3315.differential_0);
    float3  _S3318 = - _S3317.differential_0;
    float3  _S3319 = - _S3316.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3320;
    (&_S3320)->primal_0 = _S3305;
    (&_S3320)->differential_0 = _S3311;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3321;
    (&_S3321)->primal_0 = _S3299;
    (&_S3321)->differential_0 = _S3311;
    s_bwd_prop_max_0(&_S3320, &_S3321, (*v_rgb_6)[int(2)]);
    float3  _S3322 = - _S3320.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3323;
    (&_S3323)->primal_0 = _S3304;
    (&_S3323)->differential_0 = _S3311;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3324;
    (&_S3324)->primal_0 = _S3299;
    (&_S3324)->differential_0 = _S3311;
    s_bwd_prop_max_0(&_S3323, &_S3324, (*v_rgb_6)[int(1)]);
    float3  _S3325 = _S3302 * (_S3322 + _S3323.differential_0);
    float3  _S3326 = _S3320.differential_0 + _S3323.differential_0;
    float3  _S3327 = make_float3 (0.5f) * - _S3326;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3328;
    (&_S3328)->primal_0 = _S3298;
    (&_S3328)->differential_0 = _S3311;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3329;
    (&_S3329)->primal_0 = _S3299;
    (&_S3329)->differential_0 = _S3311;
    s_bwd_prop_max_0(&_S3328, &_S3329, (*v_rgb_6)[int(0)]);
    float3  _S3330 = _S3327 + _S3328.differential_0;
    float3  _S3331 = _S3326 + _S3328.differential_0;
    float3  _S3332 = _S3296 * _S3331;
    float3  _S3333 = (*sh_coeffs_19)[int(15)] * _S3331;
    float3  _S3334 = _S3294 * _S3331;
    float3  _S3335 = (*sh_coeffs_19)[int(14)] * _S3331;
    float3  _S3336 = _S3292 * _S3331;
    float3  _S3337 = (*sh_coeffs_19)[int(13)] * _S3331;
    float3  _S3338 = _S3291 * _S3331;
    float3  _S3339 = (*sh_coeffs_19)[int(12)] * _S3331;
    float3  _S3340 = _S3293 * _S3331;
    float3  _S3341 = (*sh_coeffs_19)[int(11)] * _S3331;
    float3  _S3342 = _S3295 * _S3331;
    float3  _S3343 = (*sh_coeffs_19)[int(10)] * _S3331;
    float3  _S3344 = _S3297 * _S3331;
    float3  _S3345 = (*sh_coeffs_19)[int(9)] * _S3331;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S3345.x + _S3345.y + _S3345.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S3333.x + _S3333.y + _S3333.z);
    float _S3346 = _S3343.x + _S3343.y + _S3343.z;
    float _S3347 = _S3335.x + _S3335.y + _S3335.z;
    float _S3348 = _S3341.x + _S3341.y + _S3341.z;
    float _S3349 = _S3337.x + _S3337.y + _S3337.z;
    float _S3350 = _S3339.x + _S3339.y + _S3339.z;
    float _S3351 = - s_diff_fC2_T_5;
    float3  _S3352 = _S3288 * _S3331;
    float3  _S3353 = (*sh_coeffs_19)[int(8)] * _S3331;
    float3  _S3354 = _S3286 * _S3331;
    float3  _S3355 = (*sh_coeffs_19)[int(7)] * _S3331;
    float3  _S3356 = _S3285 * _S3331;
    float3  _S3357 = (*sh_coeffs_19)[int(6)] * _S3331;
    float3  _S3358 = _S3287 * _S3331;
    float3  _S3359 = (*sh_coeffs_19)[int(5)] * _S3331;
    float3  _S3360 = _S3289 * _S3331;
    float3  _S3361 = (*sh_coeffs_19)[int(4)] * _S3331;
    float _S3362 = _S3359.x + _S3359.y + _S3359.z;
    float _S3363 = _S3355.x + _S3355.y + _S3355.z;
    float _S3364 = fTmp1B_19 * _S3346 + x_51 * s_diff_fS2_T_5 + y_22 * _S3351 + 0.54627424478530884f * (_S3361.x + _S3361.y + _S3361.z);
    float _S3365 = fTmp1B_19 * _S3347 + y_22 * s_diff_fS2_T_5 + x_51 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S3353.x + _S3353.y + _S3353.z);
    float _S3366 = y_22 * - _S3365;
    float _S3367 = x_51 * _S3365;
    float _S3368 = z_19 * (1.86588168144226074f * (z_19 * _S3350) + -2.28522896766662598f * (y_22 * _S3348 + x_51 * _S3349) + 0.94617468118667603f * (_S3357.x + _S3357.y + _S3357.z));
    float3  _S3369 = make_float3 (0.48860251903533936f) * _S3331;
    float3  _S3370 = - _S3369;
    float3  _S3371 = _S3279 * _S3370;
    float3  _S3372 = (*sh_coeffs_19)[int(3)] * _S3370;
    float3  _S3373 = _S3281 * _S3369;
    float3  _S3374 = (*sh_coeffs_19)[int(2)] * _S3369;
    float3  _S3375 = _S3283 * _S3369;
    float3  _S3376 = (*sh_coeffs_19)[int(1)] * _S3369;
    float _S3377 = (_S3290 * _S3350 + 1.44530570507049561f * (fS1_19 * _S3346 + fC1_19 * _S3347) + -1.09254848957061768f * (y_22 * _S3362 + x_51 * _S3363) + _S3368 + _S3368 + _S3374.x + _S3374.y + _S3374.z) / _S3280;
    float _S3378 = _S3278 * _S3377;
    float _S3379 = (fTmp0C_19 * _S3348 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S3351 + fTmp0B_19 * _S3362 + _S3284 * _S3364 + _S3366 + _S3366 + - (_S3376.x + _S3376.y + _S3376.z)) / _S3280;
    float _S3380 = _S3278 * _S3379;
    float _S3381 = (fTmp0C_19 * _S3349 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S3363 + 2.0f * (y_22 * _S3364) + _S3367 + _S3367 + _S3372.x + _S3372.y + _S3372.z) / _S3280;
    float _S3382 = _S3278 * _S3381;
    float _S3383 = _S3276 * - _S3377 + _S3275 * - _S3379 + _S3274 * - _S3381;
    DiffPair_float_0 _S3384;
    (&_S3384)->primal_0 = _S3277;
    (&_S3384)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3384, _S3383);
    float _S3385 = _S3276 * _S3384.differential_0;
    float _S3386 = _S3275 * _S3384.differential_0;
    float _S3387 = _S3274 * _S3384.differential_0;
    float3  _S3388 = make_float3 (0.282094806432724f) * _S3331;
    float3  _S3389 = make_float3 (_S3382 + _S3387 + _S3387, _S3380 + _S3386 + _S3386, _S3378 + _S3385 + _S3385);
    float3  _S3390 = - - _S3389;
    Matrix<float, 3, 3>  _S3391 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3392;
    (&_S3392)->primal_0 = _S3272;
    (&_S3392)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3393;
    (&_S3393)->primal_0 = t_22;
    (&_S3393)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3392, &_S3393, _S3390);
    Matrix<float, 3, 3>  _S3394 = transpose_0(_S3392.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3395;
    (&_S3395)->primal_0 = _S3271;
    (&_S3395)->differential_0 = _S3311;
    s_bwd_prop_log_1(&_S3395, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3396;
    (&_S3396)->primal_0 = vert2_c_4;
    (&_S3396)->differential_0 = _S3311;
    s_bwd_length_impl_0(&_S3396, _S3395.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3397;
    (&_S3397)->primal_0 = vert1_c_4;
    (&_S3397)->differential_0 = _S3311;
    s_bwd_length_impl_0(&_S3397, _S3395.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3398;
    (&_S3398)->primal_0 = vert0_c_4;
    (&_S3398)->differential_0 = _S3311;
    s_bwd_length_impl_0(&_S3398, _S3395.differential_0.x);
    DiffPair_float_0 _S3399;
    (&_S3399)->primal_0 = _S3270;
    (&_S3399)->differential_0 = 0.0f;
    DiffPair_float_0 _S3400;
    (&_S3400)->primal_0 = _S3257;
    (&_S3400)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3399, &_S3400, 0.0f);
    DiffPair_float_0 _S3401;
    (&_S3401)->primal_0 = _S3221;
    (&_S3401)->differential_0 = 0.0f;
    DiffPair_float_0 _S3402;
    (&_S3402)->primal_0 = _S3239;
    (&_S3402)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3401, &_S3402, _S3399.differential_0);
    DiffPair_float_0 _S3403;
    (&_S3403)->primal_0 = _S3269;
    (&_S3403)->differential_0 = 0.0f;
    DiffPair_float_0 _S3404;
    (&_S3404)->primal_0 = _S3257;
    (&_S3404)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3403, &_S3404, 0.0f);
    DiffPair_float_0 _S3405;
    (&_S3405)->primal_0 = _S3221;
    (&_S3405)->differential_0 = 0.0f;
    DiffPair_float_0 _S3406;
    (&_S3406)->primal_0 = _S3239;
    (&_S3406)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3405, &_S3406, _S3403.differential_0);
    DiffPair_float_0 _S3407;
    (&_S3407)->primal_0 = _S3268;
    (&_S3407)->differential_0 = 0.0f;
    DiffPair_float_0 _S3408;
    (&_S3408)->primal_0 = _S3256;
    (&_S3408)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3407, &_S3408, 0.0f);
    DiffPair_float_0 _S3409;
    (&_S3409)->primal_0 = _S3220;
    (&_S3409)->differential_0 = 0.0f;
    DiffPair_float_0 _S3410;
    (&_S3410)->primal_0 = _S3238;
    (&_S3410)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3409, &_S3410, _S3407.differential_0);
    DiffPair_float_0 _S3411;
    (&_S3411)->primal_0 = _S3267;
    (&_S3411)->differential_0 = 0.0f;
    DiffPair_float_0 _S3412;
    (&_S3412)->primal_0 = _S3256;
    (&_S3412)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3411, &_S3412, 0.0f);
    DiffPair_float_0 _S3413;
    (&_S3413)->primal_0 = _S3220;
    (&_S3413)->differential_0 = 0.0f;
    DiffPair_float_0 _S3414;
    (&_S3414)->primal_0 = _S3238;
    (&_S3414)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3413, &_S3414, _S3411.differential_0);
    DiffPair_float_0 _S3415;
    (&_S3415)->primal_0 = _S3265;
    (&_S3415)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3415, 0.0f);
    float _S3416 = - (-1.0f * - (_S3415.differential_0 / _S3266));
    float2  _S3417 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3418;
    (&_S3418)->primal_0 = e2_0;
    (&_S3418)->differential_0 = _S3417;
    s_bwd_length_impl_1(&_S3418, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3419;
    (&_S3419)->primal_0 = e1_4;
    (&_S3419)->differential_0 = _S3417;
    s_bwd_length_impl_1(&_S3419, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3420;
    (&_S3420)->primal_0 = e0_4;
    (&_S3420)->differential_0 = _S3417;
    s_bwd_length_impl_1(&_S3420, -0.0f);
    DiffPair_float_0 _S3421;
    (&_S3421)->primal_0 = _S3263;
    (&_S3421)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3421, 0.0f);
    float _S3422 = - _S3421.differential_0;
    float2  _S3423 = _S3419.differential_0 + make_float2 (_S3261 * _S3422, _S3259 * _S3421.differential_0);
    float2  _S3424 = _S3420.differential_0 + make_float2 (_S3260 * _S3421.differential_0, _S3262 * _S3422);
    float2  _S3425 = v_uv2_0 + - _S3418.differential_0 + _S3423;
    float _S3426 = fx_28 * (_S3408.differential_0 + _S3412.differential_0 + _S3425.x);
    float2  _S3427 = make_float2 (_S3426, fy_28 * (_S3400.differential_0 + _S3404.differential_0 + _S3425.y)) + make_float2 ((*dist_coeffs_35)[int(8)] * _S3426, (*dist_coeffs_35)[int(9)] * _S3426);
    float2  _S3428 = _S3244 * _S3427;
    float _S3429 = (*dist_coeffs_35)[int(4)] * _S3427.y;
    float _S3430 = (*dist_coeffs_35)[int(5)] * _S3427.x;
    float _S3431 = _S3428.x + _S3428.y;
    float _S3432 = r2_58 * _S3431;
    float _S3433 = r2_58 * _S3432;
    float _S3434 = (*dist_coeffs_35)[int(7)] * _S3427.y + _S3429 + (*dist_coeffs_35)[int(6)] * _S3427.x + _S3430 + _S3248 * _S3431 + _S3247 * _S3432 + _S3246 * _S3433 + (*dist_coeffs_35)[int(3)] * (r2_58 * _S3433);
    float _S3435 = v_58 * _S3434;
    float _S3436 = u_58 * _S3434;
    float2  _S3437 = (_S3249 * _S3427 + make_float2 (_S3215 * (v_58 * _S3427.y) + _S3251 * _S3430 + 2.0f * (u_58 * _S3430) + _S3212 * (v_58 * _S3427.x) + _S3436 + _S3436, _S3253 * _S3429 + 2.0f * (v_58 * _S3429) + _S3252 * _S3427.y + _S3250 * _S3427.x + _S3435 + _S3435)) / _S3245;
    float2  _S3438 = _S3241 * - _S3437;
    float2  _S3439 = _S3243 * _S3437;
    float2  _S3440 = v_uv1_0 + - _S3423 + _S3424;
    float _S3441 = fx_28 * (_S3410.differential_0 + _S3414.differential_0 + _S3440.x);
    float2  _S3442 = make_float2 (_S3441, fy_28 * (_S3402.differential_0 + _S3406.differential_0 + _S3440.y)) + make_float2 ((*dist_coeffs_35)[int(8)] * _S3441, (*dist_coeffs_35)[int(9)] * _S3441);
    float2  _S3443 = _S3226 * _S3442;
    float _S3444 = (*dist_coeffs_35)[int(4)] * _S3442.y;
    float _S3445 = (*dist_coeffs_35)[int(5)] * _S3442.x;
    float _S3446 = _S3443.x + _S3443.y;
    float _S3447 = r2_57 * _S3446;
    float _S3448 = r2_57 * _S3447;
    float _S3449 = (*dist_coeffs_35)[int(7)] * _S3442.y + _S3444 + (*dist_coeffs_35)[int(6)] * _S3442.x + _S3445 + _S3230 * _S3446 + _S3229 * _S3447 + _S3228 * _S3448 + (*dist_coeffs_35)[int(3)] * (r2_57 * _S3448);
    float _S3450 = v_57 * _S3449;
    float _S3451 = u_57 * _S3449;
    float2  _S3452 = (_S3231 * _S3442 + make_float2 (_S3215 * (v_57 * _S3442.y) + _S3233 * _S3445 + 2.0f * (u_57 * _S3445) + _S3212 * (v_57 * _S3442.x) + _S3451 + _S3451, _S3235 * _S3444 + 2.0f * (v_57 * _S3444) + _S3234 * _S3442.y + _S3232 * _S3442.x + _S3450 + _S3450)) / _S3227;
    float2  _S3453 = _S3223 * - _S3452;
    float2  _S3454 = _S3225 * _S3452;
    float _S3455 = _S3453.x + _S3453.y;
    float2  _S3456 = v_uv0_0 + _S3418.differential_0 + - _S3424;
    float _S3457 = fx_28 * (_S3409.differential_0 + _S3413.differential_0 + _S3456.x);
    float2  _S3458 = make_float2 (_S3457, fy_28 * (_S3401.differential_0 + _S3405.differential_0 + _S3456.y)) + make_float2 ((*dist_coeffs_35)[int(8)] * _S3457, (*dist_coeffs_35)[int(9)] * _S3457);
    float2  _S3459 = _S3206 * _S3458;
    float _S3460 = (*dist_coeffs_35)[int(4)] * _S3458.y;
    float _S3461 = (*dist_coeffs_35)[int(5)] * _S3458.x;
    float _S3462 = _S3459.x + _S3459.y;
    float _S3463 = r2_56 * _S3462;
    float _S3464 = r2_56 * _S3463;
    float _S3465 = (*dist_coeffs_35)[int(7)] * _S3458.y + _S3460 + (*dist_coeffs_35)[int(6)] * _S3458.x + _S3461 + _S3210 * _S3462 + _S3209 * _S3463 + _S3208 * _S3464 + (*dist_coeffs_35)[int(3)] * (r2_56 * _S3464);
    float _S3466 = v_56 * _S3465;
    float _S3467 = u_56 * _S3465;
    float2  _S3468 = (_S3211 * _S3458 + make_float2 (_S3215 * (v_56 * _S3458.y) + _S3214 * _S3461 + 2.0f * (u_56 * _S3461) + _S3212 * (v_56 * _S3458.x) + _S3467 + _S3467, _S3217 * _S3460 + 2.0f * (v_56 * _S3460) + _S3216 * _S3458.y + _S3213 * _S3458.x + _S3466 + _S3466)) / _S3207;
    float2  _S3469 = _S3203 * - _S3468;
    float2  _S3470 = _S3205 * _S3468;
    float _S3471 = _S3469.x + _S3469.y;
    float3  _S3472 = _S3317.differential_0 + _S3396.differential_0 + make_float3 (_S3439.x, _S3439.y, _S3438.x + _S3438.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3473;
    (&_S3473)->primal_0 = R_23;
    (&_S3473)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3474;
    (&_S3474)->primal_0 = vert2_1;
    (&_S3474)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3473, &_S3474, _S3472);
    float3  _S3475 = _S3316.differential_0 + _S3397.differential_0 + make_float3 (_S3454.x, _S3454.y, _S3455);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3476;
    (&_S3476)->primal_0 = R_23;
    (&_S3476)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3477;
    (&_S3477)->primal_0 = vert1_1;
    (&_S3477)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3476, &_S3477, _S3475);
    float3  _S3478 = _S3318 + _S3319 + _S3398.differential_0 + make_float3 (_S3470.x, _S3470.y, _S3471);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3479;
    (&_S3479)->primal_0 = R_23;
    (&_S3479)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3480;
    (&_S3480)->primal_0 = vert0_1;
    (&_S3480)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3479, &_S3480, _S3478);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3481;
    (&_S3481)->primal_0 = _S3197;
    (&_S3481)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3482;
    (&_S3482)->primal_0 = _S3202;
    (&_S3482)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3481, &_S3482, _S3474.differential_0);
    float _S3483 = - _S3482.differential_0.y;
    float _S3484 = _S3201 * _S3482.differential_0.x;
    float _S3485 = - (_S3193 * _S3482.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3486;
    (&_S3486)->primal_0 = _S3197;
    (&_S3486)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3487;
    (&_S3487)->primal_0 = _S3200;
    (&_S3487)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3486, &_S3487, _S3477.differential_0);
    float _S3488 = _S3193 * _S3487.differential_0.x;
    float _S3489 = _S3199 * _S3487.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3490;
    (&_S3490)->primal_0 = _S3197;
    (&_S3490)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3491;
    (&_S3491)->primal_0 = _S3198;
    (&_S3491)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3490, &_S3491, _S3480.differential_0);
    Matrix<float, 3, 3>  _S3492 = transpose_0(_S3481.differential_0 + _S3486.differential_0 + _S3490.differential_0);
    float _S3493 = 2.0f * - _S3492.rows[int(2)].z;
    float _S3494 = 2.0f * _S3492.rows[int(2)].y;
    float _S3495 = 2.0f * _S3492.rows[int(2)].x;
    float _S3496 = 2.0f * _S3492.rows[int(1)].z;
    float _S3497 = 2.0f * - _S3492.rows[int(1)].y;
    float _S3498 = 2.0f * _S3492.rows[int(1)].x;
    float _S3499 = 2.0f * _S3492.rows[int(0)].z;
    float _S3500 = 2.0f * _S3492.rows[int(0)].y;
    float _S3501 = 2.0f * - _S3492.rows[int(0)].x;
    float _S3502 = - _S3498 + _S3500;
    float _S3503 = _S3495 + - _S3499;
    float _S3504 = - _S3494 + _S3496;
    float _S3505 = _S3494 + _S3496;
    float _S3506 = _S3495 + _S3499;
    float _S3507 = _S3498 + _S3500;
    float _S3508 = quat_25.w * (_S3497 + _S3501);
    float _S3509 = quat_25.z * (_S3493 + _S3501);
    float _S3510 = quat_25.y * (_S3493 + _S3497);
    float _S3511 = quat_25.x * _S3502 + quat_25.z * _S3505 + quat_25.y * _S3506 + _S3508 + _S3508;
    float _S3512 = quat_25.x * _S3503 + quat_25.w * _S3505 + quat_25.y * _S3507 + _S3509 + _S3509;
    float _S3513 = quat_25.x * _S3504 + quat_25.w * _S3506 + quat_25.z * _S3507 + _S3510 + _S3510;
    float _S3514 = quat_25.w * _S3502 + quat_25.z * _S3503 + quat_25.y * _S3504;
    float _S3515 = _S3485 + _S3488;
    float _S3516 = 0.5f * - _S3515;
    float _S3517 = _S3483 + _S3487.differential_0.y;
    DiffPair_float_0 _S3518;
    (&_S3518)->primal_0 = _S3194;
    (&_S3518)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3518, _S3517);
    float _S3519 = _S3516 + _S3518.differential_0;
    float _S3520 = _S3484 + _S3489 + _S3491.differential_0.x;
    DiffPair_float_0 _S3521;
    (&_S3521)->primal_0 = _S3192;
    (&_S3521)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3521, _S3520);
    float _S3522 = _S3516 + _S3521.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3523;
    (&_S3523)->primal_0 = R_23;
    (&_S3523)->differential_0 = _S3391;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3524;
    (&_S3524)->primal_0 = mean_24;
    (&_S3524)->differential_0 = _S3311;
    s_bwd_prop_mul_1(&_S3523, &_S3524, _S3313.differential_0);
    float3  _S3525 = _S3393.differential_0 + _S3472 + _S3475 + _S3478 + _S3313.differential_0;
    Matrix<float, 3, 3>  _S3526 = _S3394 + _S3473.differential_0 + _S3476.differential_0 + _S3479.differential_0 + _S3523.differential_0;
    FixedArray<float3 , 2>  _S3527;
    _S3527[int(0)] = _S3311;
    _S3527[int(1)] = _S3311;
    _S3527[int(1)] = _S3325;
    _S3527[int(0)] = _S3330;
    FixedArray<float3 , 16>  _S3528;
    _S3528[int(0)] = _S3311;
    _S3528[int(1)] = _S3311;
    _S3528[int(2)] = _S3311;
    _S3528[int(3)] = _S3311;
    _S3528[int(4)] = _S3311;
    _S3528[int(5)] = _S3311;
    _S3528[int(6)] = _S3311;
    _S3528[int(7)] = _S3311;
    _S3528[int(8)] = _S3311;
    _S3528[int(9)] = _S3311;
    _S3528[int(10)] = _S3311;
    _S3528[int(11)] = _S3311;
    _S3528[int(12)] = _S3311;
    _S3528[int(13)] = _S3311;
    _S3528[int(14)] = _S3311;
    _S3528[int(15)] = _S3311;
    _S3528[int(15)] = _S3332;
    _S3528[int(14)] = _S3334;
    _S3528[int(13)] = _S3336;
    _S3528[int(12)] = _S3338;
    _S3528[int(11)] = _S3340;
    _S3528[int(10)] = _S3342;
    _S3528[int(9)] = _S3344;
    _S3528[int(8)] = _S3352;
    _S3528[int(7)] = _S3354;
    _S3528[int(6)] = _S3356;
    _S3528[int(5)] = _S3358;
    _S3528[int(4)] = _S3360;
    _S3528[int(3)] = _S3371;
    _S3528[int(2)] = _S3373;
    _S3528[int(1)] = _S3375;
    _S3528[int(0)] = _S3388;
    float2  _S3529 = v_out_hardness_0 + make_float2 (0.0f, _S3416);
    float3  _S3530 = make_float3 (_S3522, _S3519, _S3515);
    float4  _S3531 = make_float4 (0.0f);
    *&((&_S3531)->w) = _S3511;
    *&((&_S3531)->z) = _S3512;
    *&((&_S3531)->y) = _S3513;
    *&((&_S3531)->x) = _S3514;
    *v_mean_7 = _S3389 + _S3474.differential_0 + _S3477.differential_0 + _S3480.differential_0 + _S3524.differential_0;
    *v_quat_6 = _S3531;
    *v_scale_6 = _S3530;
    *v_hardness_0 = _S3529;
    *v_sh_coeffs_5 = _S3528;
    *v_ch_coeffs_0 = _S3527;
    *v_R_6 = _S3526;
    *v_t_6 = _S3525;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_25, float4  quat_26, float3  scale_25, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_24, float3  t_23, float fx_29, float fy_29, float cx_24, float cy_24, FixedArray<float, 10>  * dist_coeffs_36, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_1(R_24, mean_25) + t_23;
    float _S3532 = scale_25.x;
    float _S3533 = s_primal_ctx_exp_1(_S3532);
    float _S3534 = scale_25.y;
    float _S3535 = s_primal_ctx_exp_1(_S3534);
    float sz_6 = scale_25.z - 0.5f * (_S3532 + _S3534);
    float _S3536 = quat_26.y;
    float x2_26 = _S3536 * _S3536;
    float y2_26 = quat_26.z * quat_26.z;
    float z2_46 = quat_26.w * quat_26.w;
    float xy_26 = quat_26.y * quat_26.z;
    float xz_26 = quat_26.y * quat_26.w;
    float yz_26 = quat_26.z * quat_26.w;
    float wx_26 = quat_26.x * quat_26.y;
    float wy_26 = quat_26.x * quat_26.z;
    float wz_26 = quat_26.x * quat_26.w;
    Matrix<float, 3, 3>  _S3537 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_46), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_46), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
    float3  _S3538 = make_float3 (_S3533, 0.0f, 0.0f);
    float3  vert0_2 = s_primal_ctx_mul_1(_S3537, _S3538) + mean_25;
    float _S3539 = -0.5f + sz_6;
    float3  _S3540 = make_float3 (_S3533 * _S3539, _S3535, 0.0f);
    float3  vert1_2 = s_primal_ctx_mul_1(_S3537, _S3540) + mean_25;
    float _S3541 = -0.5f - sz_6;
    float3  _S3542 = make_float3 (_S3533 * _S3541, - _S3535, 0.0f);
    float3  vert2_2 = s_primal_ctx_mul_1(_S3537, _S3542) + mean_25;
    float3  vert0_c_5 = s_primal_ctx_mul_1(R_24, vert0_2) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_1(R_24, vert1_2) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_1(R_24, vert2_2) + t_23;
    float _S3543 = length_1(vert0_c_5);
    float _S3544 = length_1(vert1_c_5);
    float _S3545 = length_1(vert2_c_5);
    float2  _S3546 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S3547 = length_0(_S3546);
    float _S3548 = vert0_c_5.z;
    float _S3549 = s_primal_ctx_atan2_0(_S3547, _S3548);
    bool _S3550 = _S3549 < 0.00100000004749745f;
    float k_9;
    float _S3551;
    float _S3552;
    float _S3553;
    if(_S3550)
    {
        float _S3554 = 1.0f - _S3549 * _S3549 / 3.0f;
        float _S3555 = _S3548 * _S3548;
        k_9 = _S3554 / _S3548;
        _S3551 = 0.0f;
        _S3552 = _S3555;
        _S3553 = _S3554;
    }
    else
    {
        float _S3556 = _S3547 * _S3547;
        k_9 = _S3549 / _S3547;
        _S3551 = _S3556;
        _S3552 = 0.0f;
        _S3553 = 0.0f;
    }
    float2  _S3557 = make_float2 (k_9);
    float2  _S3558 = _S3546 * make_float2 (k_9);
    float u_59 = _S3558.x;
    float v_59 = _S3558.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S3559 = (*dist_coeffs_36)[int(2)] + r2_59 * (*dist_coeffs_36)[int(3)];
    float _S3560 = (*dist_coeffs_36)[int(1)] + r2_59 * _S3559;
    float _S3561 = (*dist_coeffs_36)[int(0)] + r2_59 * _S3560;
    float radial_4 = 1.0f + r2_59 * _S3561;
    float2  _S3562 = make_float2 (radial_4);
    float _S3563 = 2.0f * (*dist_coeffs_36)[int(4)];
    float _S3564 = _S3563 * u_59;
    float _S3565 = 2.0f * u_59;
    float _S3566 = 2.0f * (*dist_coeffs_36)[int(5)];
    float _S3567 = _S3566 * u_59;
    float _S3568 = 2.0f * v_59;
    float2  _S3569 = _S3558 * make_float2 (radial_4) + make_float2 (_S3564 * v_59 + (*dist_coeffs_36)[int(5)] * (r2_59 + _S3565 * u_59) + (*dist_coeffs_36)[int(6)] * r2_59, _S3567 * v_59 + (*dist_coeffs_36)[int(4)] * (r2_59 + _S3568 * v_59) + (*dist_coeffs_36)[int(7)] * r2_59);
    float2  _S3570 = _S3569 + make_float2 ((*dist_coeffs_36)[int(8)] * _S3569.x + (*dist_coeffs_36)[int(9)] * _S3569.y, 0.0f);
    float _S3571 = fx_29 * _S3570.x + cx_24;
    float _S3572 = fy_29 * _S3570.y + cy_24;
    float2  _S3573 = make_float2 (_S3571, _S3572);
    float2  _S3574 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S3575 = length_0(_S3574);
    float _S3576 = vert1_c_5.z;
    float _S3577 = s_primal_ctx_atan2_0(_S3575, _S3576);
    bool _S3578 = _S3577 < 0.00100000004749745f;
    float _S3579;
    float _S3580;
    float _S3581;
    if(_S3578)
    {
        float _S3582 = 1.0f - _S3577 * _S3577 / 3.0f;
        float _S3583 = _S3576 * _S3576;
        k_9 = _S3582 / _S3576;
        _S3579 = 0.0f;
        _S3580 = _S3583;
        _S3581 = _S3582;
    }
    else
    {
        float _S3584 = _S3575 * _S3575;
        k_9 = _S3577 / _S3575;
        _S3579 = _S3584;
        _S3580 = 0.0f;
        _S3581 = 0.0f;
    }
    float2  _S3585 = make_float2 (k_9);
    float2  _S3586 = _S3574 * make_float2 (k_9);
    float u_60 = _S3586.x;
    float v_60 = _S3586.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S3587 = (*dist_coeffs_36)[int(2)] + r2_60 * (*dist_coeffs_36)[int(3)];
    float _S3588 = (*dist_coeffs_36)[int(1)] + r2_60 * _S3587;
    float _S3589 = (*dist_coeffs_36)[int(0)] + r2_60 * _S3588;
    float radial_5 = 1.0f + r2_60 * _S3589;
    float2  _S3590 = make_float2 (radial_5);
    float _S3591 = _S3563 * u_60;
    float _S3592 = 2.0f * u_60;
    float _S3593 = _S3566 * u_60;
    float _S3594 = 2.0f * v_60;
    float2  _S3595 = _S3586 * make_float2 (radial_5) + make_float2 (_S3591 * v_60 + (*dist_coeffs_36)[int(5)] * (r2_60 + _S3592 * u_60) + (*dist_coeffs_36)[int(6)] * r2_60, _S3593 * v_60 + (*dist_coeffs_36)[int(4)] * (r2_60 + _S3594 * v_60) + (*dist_coeffs_36)[int(7)] * r2_60);
    float2  _S3596 = _S3595 + make_float2 ((*dist_coeffs_36)[int(8)] * _S3595.x + (*dist_coeffs_36)[int(9)] * _S3595.y, 0.0f);
    float _S3597 = fx_29 * _S3596.x + cx_24;
    float _S3598 = fy_29 * _S3596.y + cy_24;
    float2  _S3599 = make_float2 (_S3597, _S3598);
    float2  _S3600 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S3601 = length_0(_S3600);
    float _S3602 = vert2_c_5.z;
    float _S3603 = s_primal_ctx_atan2_0(_S3601, _S3602);
    bool _S3604 = _S3603 < 0.00100000004749745f;
    float _S3605;
    float _S3606;
    float _S3607;
    if(_S3604)
    {
        float _S3608 = 1.0f - _S3603 * _S3603 / 3.0f;
        float _S3609 = _S3602 * _S3602;
        k_9 = _S3608 / _S3602;
        _S3605 = 0.0f;
        _S3606 = _S3609;
        _S3607 = _S3608;
    }
    else
    {
        float _S3610 = _S3601 * _S3601;
        k_9 = _S3603 / _S3601;
        _S3605 = _S3610;
        _S3606 = 0.0f;
        _S3607 = 0.0f;
    }
    float2  _S3611 = make_float2 (k_9);
    float2  _S3612 = _S3600 * make_float2 (k_9);
    float u_61 = _S3612.x;
    float v_61 = _S3612.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S3613 = (*dist_coeffs_36)[int(2)] + r2_61 * (*dist_coeffs_36)[int(3)];
    float _S3614 = (*dist_coeffs_36)[int(1)] + r2_61 * _S3613;
    float _S3615 = (*dist_coeffs_36)[int(0)] + r2_61 * _S3614;
    float radial_6 = 1.0f + r2_61 * _S3615;
    float2  _S3616 = make_float2 (radial_6);
    float _S3617 = _S3563 * u_61;
    float _S3618 = 2.0f * u_61;
    float _S3619 = _S3566 * u_61;
    float _S3620 = 2.0f * v_61;
    float2  _S3621 = _S3612 * make_float2 (radial_6) + make_float2 (_S3617 * v_61 + (*dist_coeffs_36)[int(5)] * (r2_61 + _S3618 * u_61) + (*dist_coeffs_36)[int(6)] * r2_61, _S3619 * v_61 + (*dist_coeffs_36)[int(4)] * (r2_61 + _S3620 * v_61) + (*dist_coeffs_36)[int(7)] * r2_61);
    float2  _S3622 = _S3621 + make_float2 ((*dist_coeffs_36)[int(8)] * _S3621.x + (*dist_coeffs_36)[int(9)] * _S3621.y, 0.0f);
    float _S3623 = fx_29 * _S3622.x + cx_24;
    float _S3624 = fy_29 * _S3622.y + cy_24;
    float2  _S3625 = make_float2 (_S3623, _S3624);
    float2  e0_5 = _S3599 - _S3573;
    float2  e1_5 = _S3625 - _S3599;
    float2  e2_1 = _S3573 - _S3625;
    float _S3626 = e0_5.x;
    float _S3627 = e1_5.y;
    float _S3628 = e0_5.y;
    float _S3629 = e1_5.x;
    float _S3630 = _S3626 * _S3627 - _S3628 * _S3629;
    float _S3631 = 1.0f - hardness_5.y;
    float _S3632 = -1.0f / _S3631;
    float _S3633 = _S3631 * _S3631;
    float _S3634 = s_primal_ctx_max_0(_S3571, _S3597);
    float _S3635 = s_primal_ctx_min_0(_S3571, _S3597);
    float _S3636 = s_primal_ctx_max_0(_S3572, _S3598);
    float _S3637 = s_primal_ctx_min_0(_S3572, _S3598);
    float3  _S3638 = make_float3 (_S3543, _S3544, _S3545) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3639 = transpose_0(R_24);
    float3  _S3640 = mean_25 - - s_primal_ctx_mul_1(_S3639, t_23);
    float _S3641 = _S3640.x;
    float _S3642 = _S3640.y;
    float _S3643 = _S3640.z;
    float _S3644 = _S3641 * _S3641 + _S3642 * _S3642 + _S3643 * _S3643;
    float _S3645 = s_primal_ctx_sqrt_0(_S3644);
    float x_52 = _S3641 / _S3645;
    float3  _S3646 = make_float3 (x_52);
    float _S3647 = _S3645 * _S3645;
    float y_23 = _S3642 / _S3645;
    float z_20 = _S3643 / _S3645;
    float3  _S3648 = make_float3 (z_20);
    float _S3649 = - y_23;
    float3  _S3650 = make_float3 (_S3649);
    float z2_47 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_52 * x_52 - y_23 * y_23;
    float _S3651 = 2.0f * x_52;
    float fS1_20 = _S3651 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_47 - 0.31539157032966614f;
    float3  _S3652 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_52;
    float3  _S3653 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3654 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3655 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3656 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_47 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3657 = 1.86588168144226074f * z2_47 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3657;
    float3  _S3658 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_52;
    float3  _S3659 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3660 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3661 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3662 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_52 * fC1_20 - y_23 * fS1_20);
    float3  _S3663 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_52 * fS1_20 + y_23 * fC1_20);
    float3  _S3664 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3649) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_52) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3665 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3666 = make_float3 (0.0f);
    float3  _S3667 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3668 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3669 = make_float3 (_S3668);
    float3  _S3670 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3668);
    float3  _S3671 = _S3667 + _S3670 + make_float3 (0.5f);
    float3  _S3672 = _S3667 - _S3670 + make_float3 (0.5f);
    float3  _S3673 = vert1_c_5 - vert0_c_5;
    float3  _S3674 = vert2_c_5 - vert0_c_5;
    float3  _S3675 = s_primal_ctx_cross_0(_S3673, _S3674);
    float3  _S3676 = normalize_0(_S3675);
    float3  _S3677 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3676, mean_c_20)))))) * v_normal_1;
    float3  _S3678 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3679;
    (&_S3679)->primal_0 = _S3676;
    (&_S3679)->differential_0 = _S3678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3680;
    (&_S3680)->primal_0 = mean_c_20;
    (&_S3680)->differential_0 = _S3678;
    s_bwd_prop_dot_0(&_S3679, &_S3680, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3681 = _S3680;
    float3  _S3682 = _S3677 + _S3679.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3683;
    (&_S3683)->primal_0 = _S3675;
    (&_S3683)->differential_0 = _S3678;
    s_bwd_normalize_impl_0(&_S3683, _S3682);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3684;
    (&_S3684)->primal_0 = _S3673;
    (&_S3684)->differential_0 = _S3678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3685;
    (&_S3685)->primal_0 = _S3674;
    (&_S3685)->differential_0 = _S3678;
    s_bwd_prop_cross_0(&_S3684, &_S3685, _S3683.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3686 = _S3684;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3687 = _S3685;
    float3  _S3688 = - _S3685.differential_0;
    float3  _S3689 = - _S3684.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3690;
    (&_S3690)->primal_0 = _S3672;
    (&_S3690)->differential_0 = _S3678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3691;
    (&_S3691)->primal_0 = _S3666;
    (&_S3691)->differential_0 = _S3678;
    s_bwd_prop_max_0(&_S3690, &_S3691, (*v_rgb_7)[int(2)]);
    float3  _S3692 = - _S3690.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3693;
    (&_S3693)->primal_0 = _S3671;
    (&_S3693)->differential_0 = _S3678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3694;
    (&_S3694)->primal_0 = _S3666;
    (&_S3694)->differential_0 = _S3678;
    s_bwd_prop_max_0(&_S3693, &_S3694, (*v_rgb_7)[int(1)]);
    float3  _S3695 = _S3669 * (_S3692 + _S3693.differential_0);
    float3  _S3696 = _S3690.differential_0 + _S3693.differential_0;
    float3  _S3697 = make_float3 (0.5f) * - _S3696;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3698;
    (&_S3698)->primal_0 = _S3665;
    (&_S3698)->differential_0 = _S3678;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3699;
    (&_S3699)->primal_0 = _S3666;
    (&_S3699)->differential_0 = _S3678;
    s_bwd_prop_max_0(&_S3698, &_S3699, (*v_rgb_7)[int(0)]);
    float3  _S3700 = _S3697 + _S3698.differential_0;
    float3  _S3701 = _S3696 + _S3698.differential_0;
    float3  _S3702 = _S3663 * _S3701;
    float3  _S3703 = (*sh_coeffs_20)[int(15)] * _S3701;
    float3  _S3704 = _S3661 * _S3701;
    float3  _S3705 = (*sh_coeffs_20)[int(14)] * _S3701;
    float3  _S3706 = _S3659 * _S3701;
    float3  _S3707 = (*sh_coeffs_20)[int(13)] * _S3701;
    float3  _S3708 = _S3658 * _S3701;
    float3  _S3709 = (*sh_coeffs_20)[int(12)] * _S3701;
    float3  _S3710 = _S3660 * _S3701;
    float3  _S3711 = (*sh_coeffs_20)[int(11)] * _S3701;
    float3  _S3712 = _S3662 * _S3701;
    float3  _S3713 = (*sh_coeffs_20)[int(10)] * _S3701;
    float3  _S3714 = _S3664 * _S3701;
    float3  _S3715 = (*sh_coeffs_20)[int(9)] * _S3701;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3715.x + _S3715.y + _S3715.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3703.x + _S3703.y + _S3703.z);
    float _S3716 = _S3713.x + _S3713.y + _S3713.z;
    float _S3717 = _S3705.x + _S3705.y + _S3705.z;
    float _S3718 = _S3711.x + _S3711.y + _S3711.z;
    float _S3719 = _S3707.x + _S3707.y + _S3707.z;
    float _S3720 = _S3709.x + _S3709.y + _S3709.z;
    float _S3721 = - s_diff_fC2_T_6;
    float3  _S3722 = _S3655 * _S3701;
    float3  _S3723 = (*sh_coeffs_20)[int(8)] * _S3701;
    float3  _S3724 = _S3653 * _S3701;
    float3  _S3725 = (*sh_coeffs_20)[int(7)] * _S3701;
    float3  _S3726 = _S3652 * _S3701;
    float3  _S3727 = (*sh_coeffs_20)[int(6)] * _S3701;
    float3  _S3728 = _S3654 * _S3701;
    float3  _S3729 = (*sh_coeffs_20)[int(5)] * _S3701;
    float3  _S3730 = _S3656 * _S3701;
    float3  _S3731 = (*sh_coeffs_20)[int(4)] * _S3701;
    float _S3732 = _S3729.x + _S3729.y + _S3729.z;
    float _S3733 = _S3725.x + _S3725.y + _S3725.z;
    float _S3734 = fTmp1B_20 * _S3716 + x_52 * s_diff_fS2_T_6 + y_23 * _S3721 + 0.54627424478530884f * (_S3731.x + _S3731.y + _S3731.z);
    float _S3735 = fTmp1B_20 * _S3717 + y_23 * s_diff_fS2_T_6 + x_52 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3723.x + _S3723.y + _S3723.z);
    float _S3736 = y_23 * - _S3735;
    float _S3737 = x_52 * _S3735;
    float _S3738 = z_20 * (1.86588168144226074f * (z_20 * _S3720) + -2.28522896766662598f * (y_23 * _S3718 + x_52 * _S3719) + 0.94617468118667603f * (_S3727.x + _S3727.y + _S3727.z));
    float3  _S3739 = make_float3 (0.48860251903533936f) * _S3701;
    float3  _S3740 = - _S3739;
    float3  _S3741 = _S3646 * _S3740;
    float3  _S3742 = (*sh_coeffs_20)[int(3)] * _S3740;
    float3  _S3743 = _S3648 * _S3739;
    float3  _S3744 = (*sh_coeffs_20)[int(2)] * _S3739;
    float3  _S3745 = _S3650 * _S3739;
    float3  _S3746 = (*sh_coeffs_20)[int(1)] * _S3739;
    float _S3747 = (_S3657 * _S3720 + 1.44530570507049561f * (fS1_20 * _S3716 + fC1_20 * _S3717) + -1.09254848957061768f * (y_23 * _S3732 + x_52 * _S3733) + _S3738 + _S3738 + _S3744.x + _S3744.y + _S3744.z) / _S3647;
    float _S3748 = _S3645 * _S3747;
    float _S3749 = (fTmp0C_20 * _S3718 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3721 + fTmp0B_20 * _S3732 + _S3651 * _S3734 + _S3736 + _S3736 + - (_S3746.x + _S3746.y + _S3746.z)) / _S3647;
    float _S3750 = _S3645 * _S3749;
    float _S3751 = (fTmp0C_20 * _S3719 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3733 + 2.0f * (y_23 * _S3734) + _S3737 + _S3737 + _S3742.x + _S3742.y + _S3742.z) / _S3647;
    float _S3752 = _S3645 * _S3751;
    float _S3753 = _S3643 * - _S3747 + _S3642 * - _S3749 + _S3641 * - _S3751;
    DiffPair_float_0 _S3754;
    (&_S3754)->primal_0 = _S3644;
    (&_S3754)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3754, _S3753);
    float _S3755 = _S3643 * _S3754.differential_0;
    float _S3756 = _S3642 * _S3754.differential_0;
    float _S3757 = _S3641 * _S3754.differential_0;
    float3  _S3758 = make_float3 (0.282094806432724f) * _S3701;
    float3  _S3759 = make_float3 (_S3752 + _S3757 + _S3757, _S3750 + _S3756 + _S3756, _S3748 + _S3755 + _S3755);
    float3  _S3760 = - - _S3759;
    Matrix<float, 3, 3>  _S3761 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3762;
    (&_S3762)->primal_0 = _S3639;
    (&_S3762)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3763;
    (&_S3763)->primal_0 = t_23;
    (&_S3763)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3762, &_S3763, _S3760);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3764 = _S3763;
    Matrix<float, 3, 3>  _S3765 = transpose_0(_S3762.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3766;
    (&_S3766)->primal_0 = _S3638;
    (&_S3766)->differential_0 = _S3678;
    s_bwd_prop_log_1(&_S3766, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3767 = _S3766;
    DiffPair_float_0 _S3768;
    (&_S3768)->primal_0 = _S3637;
    (&_S3768)->differential_0 = 0.0f;
    DiffPair_float_0 _S3769;
    (&_S3769)->primal_0 = _S3624;
    (&_S3769)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3768, &_S3769, 0.0f);
    DiffPair_float_0 _S3770;
    (&_S3770)->primal_0 = _S3572;
    (&_S3770)->differential_0 = 0.0f;
    DiffPair_float_0 _S3771;
    (&_S3771)->primal_0 = _S3598;
    (&_S3771)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3770, &_S3771, _S3768.differential_0);
    DiffPair_float_0 _S3772;
    (&_S3772)->primal_0 = _S3636;
    (&_S3772)->differential_0 = 0.0f;
    DiffPair_float_0 _S3773;
    (&_S3773)->primal_0 = _S3624;
    (&_S3773)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3772, &_S3773, 0.0f);
    DiffPair_float_0 _S3774;
    (&_S3774)->primal_0 = _S3572;
    (&_S3774)->differential_0 = 0.0f;
    DiffPair_float_0 _S3775;
    (&_S3775)->primal_0 = _S3598;
    (&_S3775)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3774, &_S3775, _S3772.differential_0);
    DiffPair_float_0 _S3776;
    (&_S3776)->primal_0 = _S3635;
    (&_S3776)->differential_0 = 0.0f;
    DiffPair_float_0 _S3777;
    (&_S3777)->primal_0 = _S3623;
    (&_S3777)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3776, &_S3777, 0.0f);
    DiffPair_float_0 _S3778;
    (&_S3778)->primal_0 = _S3571;
    (&_S3778)->differential_0 = 0.0f;
    DiffPair_float_0 _S3779;
    (&_S3779)->primal_0 = _S3597;
    (&_S3779)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3778, &_S3779, _S3776.differential_0);
    DiffPair_float_0 _S3780;
    (&_S3780)->primal_0 = _S3634;
    (&_S3780)->differential_0 = 0.0f;
    DiffPair_float_0 _S3781;
    (&_S3781)->primal_0 = _S3623;
    (&_S3781)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3780, &_S3781, 0.0f);
    DiffPair_float_0 _S3782;
    (&_S3782)->primal_0 = _S3571;
    (&_S3782)->differential_0 = 0.0f;
    DiffPair_float_0 _S3783;
    (&_S3783)->primal_0 = _S3597;
    (&_S3783)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3782, &_S3783, _S3780.differential_0);
    DiffPair_float_0 _S3784;
    (&_S3784)->primal_0 = _S3632;
    (&_S3784)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3784, 0.0f);
    float _S3785 = - (-1.0f * - (_S3784.differential_0 / _S3633));
    float2  _S3786 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3787;
    (&_S3787)->primal_0 = e2_1;
    (&_S3787)->differential_0 = _S3786;
    s_bwd_length_impl_1(&_S3787, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3788;
    (&_S3788)->primal_0 = e1_5;
    (&_S3788)->differential_0 = _S3786;
    s_bwd_length_impl_1(&_S3788, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3789;
    (&_S3789)->primal_0 = e0_5;
    (&_S3789)->differential_0 = _S3786;
    s_bwd_length_impl_1(&_S3789, -0.0f);
    DiffPair_float_0 _S3790;
    (&_S3790)->primal_0 = _S3630;
    (&_S3790)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3790, 0.0f);
    float _S3791 = - _S3790.differential_0;
    float2  _S3792 = _S3788.differential_0 + make_float2 (_S3628 * _S3791, _S3626 * _S3790.differential_0);
    float2  _S3793 = _S3789.differential_0 + make_float2 (_S3627 * _S3790.differential_0, _S3629 * _S3791);
    float2  _S3794 = v_uv2_1 + - _S3787.differential_0 + _S3792;
    float _S3795 = fx_29 * (_S3777.differential_0 + _S3781.differential_0 + _S3794.x);
    float2  _S3796 = make_float2 (_S3795, fy_29 * (_S3769.differential_0 + _S3773.differential_0 + _S3794.y)) + make_float2 ((*dist_coeffs_36)[int(8)] * _S3795, (*dist_coeffs_36)[int(9)] * _S3795);
    float2  _S3797 = _S3612 * _S3796;
    float2  _S3798 = _S3616 * _S3796;
    float _S3799 = (*dist_coeffs_36)[int(4)] * _S3796.y;
    float _S3800 = (*dist_coeffs_36)[int(5)] * _S3796.x;
    float _S3801 = _S3797.x + _S3797.y;
    float _S3802 = r2_61 * _S3801;
    float _S3803 = r2_61 * _S3802;
    float _S3804 = (*dist_coeffs_36)[int(7)] * _S3796.y + _S3799 + (*dist_coeffs_36)[int(6)] * _S3796.x + _S3800 + _S3615 * _S3801 + _S3614 * _S3802 + _S3613 * _S3803 + (*dist_coeffs_36)[int(3)] * (r2_61 * _S3803);
    float _S3805 = v_61 * _S3804;
    float _S3806 = u_61 * _S3804;
    float _S3807 = _S3620 * _S3799 + 2.0f * (v_61 * _S3799) + _S3619 * _S3796.y + _S3617 * _S3796.x + _S3805 + _S3805;
    float _S3808 = _S3566 * (v_61 * _S3796.y) + _S3618 * _S3800 + 2.0f * (u_61 * _S3800) + _S3563 * (v_61 * _S3796.x) + _S3806 + _S3806;
    float2  _S3809 = v_uv0_1 + _S3787.differential_0 + - _S3793;
    float2  _S3810 = v_out_hardness_1 + make_float2 (0.0f, _S3785);
    float _S3811 = _S3779.differential_0 + _S3783.differential_0;
    float2  _S3812 = v_uv1_1 + - _S3792 + _S3793;
    float3  _S3813 = _S3688 + _S3689;
    FixedArray<float3 , 2>  _S3814;
    _S3814[int(0)] = _S3678;
    _S3814[int(1)] = _S3678;
    _S3814[int(1)] = _S3695;
    _S3814[int(0)] = _S3700;
    float3  _S3815 = _S3814[int(0)];
    float3  _S3816 = _S3814[int(1)];
    FixedArray<float3 , 16>  _S3817;
    _S3817[int(0)] = _S3678;
    _S3817[int(1)] = _S3678;
    _S3817[int(2)] = _S3678;
    _S3817[int(3)] = _S3678;
    _S3817[int(4)] = _S3678;
    _S3817[int(5)] = _S3678;
    _S3817[int(6)] = _S3678;
    _S3817[int(7)] = _S3678;
    _S3817[int(8)] = _S3678;
    _S3817[int(9)] = _S3678;
    _S3817[int(10)] = _S3678;
    _S3817[int(11)] = _S3678;
    _S3817[int(12)] = _S3678;
    _S3817[int(13)] = _S3678;
    _S3817[int(14)] = _S3678;
    _S3817[int(15)] = _S3678;
    _S3817[int(7)] = _S3724;
    _S3817[int(0)] = _S3758;
    _S3817[int(1)] = _S3745;
    _S3817[int(2)] = _S3743;
    _S3817[int(3)] = _S3741;
    _S3817[int(4)] = _S3730;
    _S3817[int(5)] = _S3728;
    _S3817[int(6)] = _S3726;
    _S3817[int(15)] = _S3702;
    _S3817[int(8)] = _S3722;
    _S3817[int(9)] = _S3714;
    _S3817[int(10)] = _S3712;
    _S3817[int(11)] = _S3710;
    _S3817[int(12)] = _S3708;
    _S3817[int(13)] = _S3706;
    _S3817[int(14)] = _S3704;
    float3  _S3818 = _S3817[int(0)];
    float3  _S3819 = _S3817[int(1)];
    float3  _S3820 = _S3817[int(2)];
    float3  _S3821 = _S3817[int(3)];
    float3  _S3822 = _S3817[int(4)];
    float3  _S3823 = _S3817[int(5)];
    float3  _S3824 = _S3817[int(6)];
    float3  _S3825 = _S3817[int(7)];
    float3  _S3826 = _S3817[int(8)];
    float3  _S3827 = _S3817[int(9)];
    float3  _S3828 = _S3817[int(10)];
    float3  _S3829 = _S3817[int(11)];
    float3  _S3830 = _S3817[int(12)];
    float3  _S3831 = _S3817[int(13)];
    float3  _S3832 = _S3817[int(14)];
    float3  _S3833 = _S3817[int(15)];
    float _S3834 = _S3778.differential_0 + _S3782.differential_0;
    float _S3835 = _S3770.differential_0 + _S3774.differential_0;
    float _S3836 = _S3771.differential_0 + _S3775.differential_0;
    float2  _S3837 = _S3798 + make_float2 (_S3808, _S3807);
    float2  _S3838 = _S3600 * _S3837;
    float2  _S3839 = _S3611 * _S3837;
    float _S3840 = _S3838.x + _S3838.y;
    if(_S3604)
    {
        float _S3841 = _S3840 / _S3606;
        float _S3842 = _S3607 * - _S3841;
        float _S3843 = _S3603 * (0.3333333432674408f * - (_S3602 * _S3841));
        k_9 = _S3843 + _S3843;
        _S3605 = _S3842;
        _S3606 = 0.0f;
    }
    else
    {
        float _S3844 = _S3840 / _S3605;
        float _S3845 = _S3603 * - _S3844;
        k_9 = _S3601 * _S3844;
        _S3605 = 0.0f;
        _S3606 = _S3845;
    }
    DiffPair_float_0 _S3846;
    (&_S3846)->primal_0 = _S3601;
    (&_S3846)->differential_0 = 0.0f;
    DiffPair_float_0 _S3847;
    (&_S3847)->primal_0 = _S3602;
    (&_S3847)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3846, &_S3847, k_9);
    float _S3848 = _S3847.differential_0 + _S3605;
    float _S3849 = _S3846.differential_0 + _S3606;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3850;
    (&_S3850)->primal_0 = _S3600;
    (&_S3850)->differential_0 = _S3786;
    s_bwd_length_impl_1(&_S3850, _S3849);
    float2  _S3851 = _S3850.differential_0 + _S3839;
    float _S3852 = fx_29 * (_S3812.x + _S3811);
    float2  _S3853 = make_float2 (_S3852, fy_29 * (_S3812.y + _S3836)) + make_float2 ((*dist_coeffs_36)[int(8)] * _S3852, (*dist_coeffs_36)[int(9)] * _S3852);
    float2  _S3854 = _S3586 * _S3853;
    float _S3855 = (*dist_coeffs_36)[int(4)] * _S3853.y;
    float _S3856 = (*dist_coeffs_36)[int(5)] * _S3853.x;
    float _S3857 = _S3854.x + _S3854.y;
    float _S3858 = r2_60 * _S3857;
    float _S3859 = r2_60 * _S3858;
    float _S3860 = (*dist_coeffs_36)[int(7)] * _S3853.y + _S3855 + (*dist_coeffs_36)[int(6)] * _S3853.x + _S3856 + _S3589 * _S3857 + _S3588 * _S3858 + _S3587 * _S3859 + (*dist_coeffs_36)[int(3)] * (r2_60 * _S3859);
    float _S3861 = v_60 * _S3860;
    float _S3862 = u_60 * _S3860;
    float3  _S3863 = _S3687.differential_0 + make_float3 (_S3851.x, _S3851.y, _S3848);
    float2  _S3864 = _S3590 * _S3853 + make_float2 (_S3566 * (v_60 * _S3853.y) + _S3592 * _S3856 + 2.0f * (u_60 * _S3856) + _S3563 * (v_60 * _S3853.x) + _S3862 + _S3862, _S3594 * _S3855 + 2.0f * (v_60 * _S3855) + _S3593 * _S3853.y + _S3591 * _S3853.x + _S3861 + _S3861);
    float2  _S3865 = _S3574 * _S3864;
    float2  _S3866 = _S3585 * _S3864;
    float _S3867 = _S3865.x + _S3865.y;
    if(_S3578)
    {
        float _S3868 = _S3867 / _S3580;
        float _S3869 = _S3581 * - _S3868;
        float _S3870 = _S3577 * (0.3333333432674408f * - (_S3576 * _S3868));
        k_9 = _S3870 + _S3870;
        _S3579 = _S3869;
        _S3580 = 0.0f;
    }
    else
    {
        float _S3871 = _S3867 / _S3579;
        float _S3872 = _S3577 * - _S3871;
        k_9 = _S3575 * _S3871;
        _S3579 = 0.0f;
        _S3580 = _S3872;
    }
    DiffPair_float_0 _S3873;
    (&_S3873)->primal_0 = _S3575;
    (&_S3873)->differential_0 = 0.0f;
    DiffPair_float_0 _S3874;
    (&_S3874)->primal_0 = _S3576;
    (&_S3874)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3873, &_S3874, k_9);
    float _S3875 = _S3874.differential_0 + _S3579;
    float _S3876 = _S3873.differential_0 + _S3580;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3877;
    (&_S3877)->primal_0 = _S3574;
    (&_S3877)->differential_0 = _S3786;
    s_bwd_length_impl_1(&_S3877, _S3876);
    float2  _S3878 = _S3877.differential_0 + _S3866;
    float _S3879 = fx_29 * (_S3809.x + _S3834);
    float2  _S3880 = make_float2 (_S3879, fy_29 * (_S3809.y + _S3835)) + make_float2 ((*dist_coeffs_36)[int(8)] * _S3879, (*dist_coeffs_36)[int(9)] * _S3879);
    float2  _S3881 = _S3558 * _S3880;
    float _S3882 = (*dist_coeffs_36)[int(4)] * _S3880.y;
    float _S3883 = (*dist_coeffs_36)[int(5)] * _S3880.x;
    float _S3884 = _S3881.x + _S3881.y;
    float _S3885 = r2_59 * _S3884;
    float _S3886 = r2_59 * _S3885;
    float _S3887 = (*dist_coeffs_36)[int(7)] * _S3880.y + _S3882 + (*dist_coeffs_36)[int(6)] * _S3880.x + _S3883 + _S3561 * _S3884 + _S3560 * _S3885 + _S3559 * _S3886 + (*dist_coeffs_36)[int(3)] * (r2_59 * _S3886);
    float _S3888 = v_59 * _S3887;
    float _S3889 = u_59 * _S3887;
    float3  _S3890 = _S3686.differential_0 + make_float3 (_S3878.x, _S3878.y, _S3875);
    float2  _S3891 = _S3562 * _S3880 + make_float2 (_S3566 * (v_59 * _S3880.y) + _S3565 * _S3883 + 2.0f * (u_59 * _S3883) + _S3563 * (v_59 * _S3880.x) + _S3889 + _S3889, _S3568 * _S3882 + 2.0f * (v_59 * _S3882) + _S3567 * _S3880.y + _S3564 * _S3880.x + _S3888 + _S3888);
    float2  _S3892 = _S3546 * _S3891;
    float2  _S3893 = _S3557 * _S3891;
    float _S3894 = _S3892.x + _S3892.y;
    if(_S3550)
    {
        float _S3895 = _S3894 / _S3552;
        float _S3896 = _S3553 * - _S3895;
        float _S3897 = _S3549 * (0.3333333432674408f * - (_S3548 * _S3895));
        k_9 = _S3897 + _S3897;
        _S3551 = _S3896;
        _S3552 = 0.0f;
    }
    else
    {
        float _S3898 = _S3894 / _S3551;
        float _S3899 = _S3549 * - _S3898;
        k_9 = _S3547 * _S3898;
        _S3551 = 0.0f;
        _S3552 = _S3899;
    }
    DiffPair_float_0 _S3900;
    (&_S3900)->primal_0 = _S3547;
    (&_S3900)->differential_0 = 0.0f;
    DiffPair_float_0 _S3901;
    (&_S3901)->primal_0 = _S3548;
    (&_S3901)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3900, &_S3901, k_9);
    float _S3902 = _S3901.differential_0 + _S3551;
    float _S3903 = _S3900.differential_0 + _S3552;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3904;
    (&_S3904)->primal_0 = _S3546;
    (&_S3904)->differential_0 = _S3786;
    s_bwd_length_impl_1(&_S3904, _S3903);
    float2  _S3905 = _S3904.differential_0 + _S3893;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3906;
    (&_S3906)->primal_0 = vert2_c_5;
    (&_S3906)->differential_0 = _S3678;
    s_bwd_length_impl_0(&_S3906, _S3767.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3907;
    (&_S3907)->primal_0 = vert1_c_5;
    (&_S3907)->differential_0 = _S3678;
    s_bwd_length_impl_0(&_S3907, _S3767.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3908;
    (&_S3908)->primal_0 = vert0_c_5;
    (&_S3908)->differential_0 = _S3678;
    s_bwd_length_impl_0(&_S3908, _S3767.differential_0.x);
    float3  _S3909 = _S3906.differential_0 + _S3863;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3910;
    (&_S3910)->primal_0 = R_24;
    (&_S3910)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3911;
    (&_S3911)->primal_0 = vert2_2;
    (&_S3911)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3910, &_S3911, _S3909);
    float3  _S3912 = _S3907.differential_0 + _S3890;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3913;
    (&_S3913)->primal_0 = R_24;
    (&_S3913)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3914;
    (&_S3914)->primal_0 = vert1_2;
    (&_S3914)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3913, &_S3914, _S3912);
    float3  _S3915 = _S3908.differential_0 + _S3813 + make_float3 (_S3905.x, _S3905.y, _S3902);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3916;
    (&_S3916)->primal_0 = R_24;
    (&_S3916)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3917;
    (&_S3917)->primal_0 = vert0_2;
    (&_S3917)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3916, &_S3917, _S3915);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3918;
    (&_S3918)->primal_0 = _S3537;
    (&_S3918)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3919;
    (&_S3919)->primal_0 = _S3542;
    (&_S3919)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3918, &_S3919, _S3911.differential_0);
    float _S3920 = - _S3919.differential_0.y;
    float _S3921 = _S3541 * _S3919.differential_0.x;
    float _S3922 = - (_S3533 * _S3919.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3923;
    (&_S3923)->primal_0 = _S3537;
    (&_S3923)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3924;
    (&_S3924)->primal_0 = _S3540;
    (&_S3924)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3923, &_S3924, _S3914.differential_0);
    float _S3925 = _S3533 * _S3924.differential_0.x;
    float _S3926 = _S3539 * _S3924.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3927;
    (&_S3927)->primal_0 = _S3537;
    (&_S3927)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3928;
    (&_S3928)->primal_0 = _S3538;
    (&_S3928)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3927, &_S3928, _S3917.differential_0);
    Matrix<float, 3, 3>  _S3929 = transpose_0(_S3918.differential_0 + _S3923.differential_0 + _S3927.differential_0);
    float _S3930 = 2.0f * - _S3929.rows[int(2)].z;
    float _S3931 = 2.0f * _S3929.rows[int(2)].y;
    float _S3932 = 2.0f * _S3929.rows[int(2)].x;
    float _S3933 = 2.0f * _S3929.rows[int(1)].z;
    float _S3934 = 2.0f * - _S3929.rows[int(1)].y;
    float _S3935 = 2.0f * _S3929.rows[int(1)].x;
    float _S3936 = 2.0f * _S3929.rows[int(0)].z;
    float _S3937 = 2.0f * _S3929.rows[int(0)].y;
    float _S3938 = 2.0f * - _S3929.rows[int(0)].x;
    float _S3939 = - _S3935 + _S3937;
    float _S3940 = _S3932 + - _S3936;
    float _S3941 = - _S3931 + _S3933;
    float _S3942 = _S3931 + _S3933;
    float _S3943 = _S3932 + _S3936;
    float _S3944 = _S3935 + _S3937;
    float _S3945 = quat_26.w * (_S3934 + _S3938);
    float _S3946 = quat_26.z * (_S3930 + _S3938);
    float _S3947 = quat_26.y * (_S3930 + _S3934);
    float _S3948 = quat_26.x * _S3939 + quat_26.z * _S3942 + quat_26.y * _S3943 + _S3945 + _S3945;
    float _S3949 = quat_26.x * _S3940 + quat_26.w * _S3942 + quat_26.y * _S3944 + _S3946 + _S3946;
    float _S3950 = quat_26.x * _S3941 + quat_26.w * _S3943 + quat_26.z * _S3944 + _S3947 + _S3947;
    float _S3951 = quat_26.w * _S3939 + quat_26.z * _S3940 + quat_26.y * _S3941;
    float _S3952 = _S3922 + _S3925;
    float _S3953 = 0.5f * - _S3952;
    float _S3954 = _S3920 + _S3924.differential_0.y;
    DiffPair_float_0 _S3955;
    (&_S3955)->primal_0 = _S3534;
    (&_S3955)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3955, _S3954);
    float _S3956 = _S3953 + _S3955.differential_0;
    float _S3957 = _S3921 + _S3926 + _S3928.differential_0.x;
    DiffPair_float_0 _S3958;
    (&_S3958)->primal_0 = _S3532;
    (&_S3958)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3958, _S3957);
    float _S3959 = _S3953 + _S3958.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3960;
    (&_S3960)->primal_0 = mean_c_20;
    (&_S3960)->differential_0 = _S3678;
    s_bwd_length_impl_0(&_S3960, 0.0f);
    float3  _S3961 = _S3960.differential_0 + _S3681.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3962;
    (&_S3962)->primal_0 = R_24;
    (&_S3962)->differential_0 = _S3761;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3963;
    (&_S3963)->primal_0 = mean_25;
    (&_S3963)->differential_0 = _S3678;
    s_bwd_prop_mul_1(&_S3962, &_S3963, _S3961);
    float3  _S3964 = _S3909 + _S3912 + _S3915 + _S3961 + _S3764.differential_0;
    Matrix<float, 3, 3>  _S3965 = _S3910.differential_0 + _S3913.differential_0 + _S3916.differential_0 + _S3962.differential_0 + _S3765;
    float3  _S3966 = make_float3 (_S3959, _S3956, _S3952);
    float4  _S3967 = make_float4 (0.0f);
    *&((&_S3967)->w) = _S3948;
    *&((&_S3967)->z) = _S3949;
    *&((&_S3967)->y) = _S3950;
    *&((&_S3967)->x) = _S3951;
    float4  _S3968 = _S3967;
    float3  _S3969 = _S3911.differential_0 + _S3914.differential_0 + _S3917.differential_0 + _S3963.differential_0 + _S3759;
    *v_mean_8 = _S3969;
    *v_quat_7 = _S3968;
    *v_scale_7 = _S3966;
    *v_hardness_1 = _S3810;
    (*v_sh_coeffs_6)[int(0)] = _S3818;
    (*v_sh_coeffs_6)[int(1)] = _S3819;
    (*v_sh_coeffs_6)[int(2)] = _S3820;
    (*v_sh_coeffs_6)[int(3)] = _S3821;
    (*v_sh_coeffs_6)[int(4)] = _S3822;
    (*v_sh_coeffs_6)[int(5)] = _S3823;
    (*v_sh_coeffs_6)[int(6)] = _S3824;
    (*v_sh_coeffs_6)[int(7)] = _S3825;
    (*v_sh_coeffs_6)[int(8)] = _S3826;
    (*v_sh_coeffs_6)[int(9)] = _S3827;
    (*v_sh_coeffs_6)[int(10)] = _S3828;
    (*v_sh_coeffs_6)[int(11)] = _S3829;
    (*v_sh_coeffs_6)[int(12)] = _S3830;
    (*v_sh_coeffs_6)[int(13)] = _S3831;
    (*v_sh_coeffs_6)[int(14)] = _S3832;
    (*v_sh_coeffs_6)[int(15)] = _S3833;
    (*v_ch_coeffs_1)[int(0)] = _S3815;
    (*v_ch_coeffs_1)[int(1)] = _S3816;
    *v_R_7 = _S3965;
    *v_t_7 = _S3964;
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
        DiffPair_float_0 _S3970 = *dpx_17;
        float _S3971 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S3971;
        float _S3972 = val_0 * (F32_log((_S3970.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3972;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_1)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3973 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3973))));
    float2  _S3974 = p_1 - v0_0;
    float2  _S3975 = normalize_1(e0_6);
    float2  _S3976 = p_1 - v1_0;
    float2  _S3977 = normalize_1(e1_6);
    float2  _S3978 = p_1 - v2_0;
    float2  _S3979 = normalize_1(e2_2);
    float _S3980 = hardness_6.x;
    float _S3981 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3974.x * _S3975.y - _S3974.y * _S3975.x)), (se_0 * (_S3976.x * _S3977.y - _S3976.y * _S3977.x))))), (se_0 * (_S3978.x * _S3979.y - _S3978.y * _S3979.x)))) / ((F32_abs((_S3973))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S3981))));
    float _S3982;
    if(a_1 <= 0.0f)
    {
        _S3982 = 0.0f;
    }
    else
    {
        _S3982 = (F32_min(((F32_pow((a_1), (_S3981)))), (0.99900001287460327f)));
    }
    return _S3980 * _S3982;
}

inline __device__ float s_primal_ctx_abs_0(float _S3983)
{
    return (F32_abs((_S3983)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S3984, float _S3985, float _S3986)
{
    return clamp_0(_S3984, _S3985, _S3986);
}

inline __device__ float s_primal_ctx_exp2_0(float _S3987)
{
    return (F32_exp2((_S3987)));
}

inline __device__ float s_primal_ctx_pow_0(float _S3988, float _S3989)
{
    return (F32_pow((_S3988), (_S3989)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S3990, DiffPair_float_0 * _S3991, float _S3992)
{
    _d_pow_0(_S3990, _S3991, _S3992);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S3993, DiffPair_float_0 * _S3994, DiffPair_float_0 * _S3995, float _S3996)
{
    _d_clamp_0(_S3993, _S3994, _S3995, _S3996);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_8)
{
    float _S3997 = length_0((*dpx_18).primal_0);
    float2  _S3998 = (*dpx_18).primal_0 * _s_dOut_8;
    float2  _S3999 = make_float2 (1.0f / _S3997) * _s_dOut_8;
    float _S4000 = - ((_S3998.x + _S3998.y) / (_S3997 * _S3997));
    float2  _S4001 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4002;
    (&_S4002)->primal_0 = (*dpx_18).primal_0;
    (&_S4002)->differential_0 = _S4001;
    s_bwd_length_impl_1(&_S4002, _S4000);
    float2  _S4003 = _S3999 + _S4002.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S4003;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4004, float2  _S4005)
{
    s_bwd_prop_normalize_impl_1(_S4004, _S4005);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_2, float _s_dOut_9)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S4006 = e0_7.x;
    float _S4007 = e1_7.y;
    float _S4008 = e0_7.y;
    float _S4009 = e1_7.x;
    float _S4010 = _S4006 * _S4007 - _S4008 * _S4009;
    float se_1 = float((F32_sign((_S4010))));
    float2  _S4011 = p_2 - (*dpv0_0).primal_0;
    float2  _S4012 = normalize_1(e0_7);
    float _S4013 = _S4011.x;
    float _S4014 = _S4012.y;
    float _S4015 = _S4011.y;
    float _S4016 = _S4012.x;
    float de0_0 = se_1 * (_S4013 * _S4014 - _S4015 * _S4016);
    float2  _S4017 = p_2 - (*dpv1_0).primal_0;
    float2  _S4018 = normalize_1(e1_7);
    float _S4019 = _S4017.x;
    float _S4020 = _S4018.y;
    float _S4021 = _S4017.y;
    float _S4022 = _S4018.x;
    float de1_0 = se_1 * (_S4019 * _S4020 - _S4021 * _S4022);
    float2  _S4023 = p_2 - (*dpv2_0).primal_0;
    float2  _S4024 = normalize_1(e2_3);
    float _S4025 = _S4023.x;
    float _S4026 = _S4024.y;
    float _S4027 = _S4023.y;
    float _S4028 = _S4024.x;
    float de2_0 = se_1 * (_S4025 * _S4026 - _S4027 * _S4028);
    float _S4029 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S4030 = s_primal_ctx_max_0(_S4029, de2_0);
    float _S4031 = s_primal_ctx_abs_0(_S4010);
    float _S4032 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S4031 / _S4032;
    float _S4033 = _S4032 * _S4032;
    float _S4034 = (*dphardness_0).primal_0.x;
    float _S4035 = (*dphardness_0).primal_0.y;
    float _S4036 = dmax_0 * dmax_0;
    float _S4037 = 1.0f + _S4030 / dmax_0;
    float _S4038 = 1.0f - s_primal_ctx_clamp_0(_S4035, 0.00499999988824129f, 0.98000001907348633f);
    float _S4039 = -1.0f / _S4038;
    float _S4040 = _S4038 * _S4038;
    float _S4041 = 1.0f - s_primal_ctx_exp2_0(_S4039);
    float a_2 = 1.0f - _S4037 * _S4041;
    bool _S4042 = a_2 <= 0.0f;
    float _S4043;
    float _S4044;
    if(_S4042)
    {
        _S4043 = 0.0f;
        _S4044 = 0.0f;
    }
    else
    {
        float _S4045 = s_primal_ctx_pow_0(a_2, _S4038);
        _S4043 = s_primal_ctx_min_0(_S4045, 0.99900001287460327f);
        _S4044 = _S4045;
    }
    float _S4046 = _S4034 * _s_dOut_9;
    float _S4047 = _S4043 * _s_dOut_9;
    if(_S4042)
    {
        _S4043 = 0.0f;
        _S4044 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4048;
        (&_S4048)->primal_0 = _S4044;
        (&_S4048)->differential_0 = 0.0f;
        DiffPair_float_0 _S4049;
        (&_S4049)->primal_0 = 0.99900001287460327f;
        (&_S4049)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4048, &_S4049, _S4046);
        DiffPair_float_0 _S4050;
        (&_S4050)->primal_0 = a_2;
        (&_S4050)->differential_0 = 0.0f;
        DiffPair_float_0 _S4051;
        (&_S4051)->primal_0 = _S4038;
        (&_S4051)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4050, &_S4051, _S4048.differential_0);
        _S4043 = _S4050.differential_0;
        _S4044 = _S4051.differential_0;
    }
    float _S4052 = - _S4043;
    float _S4053 = _S4041 * _S4052;
    float _S4054 = - (_S4037 * _S4052);
    DiffPair_float_0 _S4055;
    (&_S4055)->primal_0 = _S4039;
    (&_S4055)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4055, _S4054);
    float _S4056 = - (-1.0f * - (_S4055.differential_0 / _S4040) + _S4044);
    float _S4057 = _S4053 / _S4036;
    float s_diff_dmax_T_0 = _S4030 * - _S4057;
    float _S4058 = dmax_0 * _S4057;
    DiffPair_float_0 _S4059;
    (&_S4059)->primal_0 = _S4035;
    (&_S4059)->differential_0 = 0.0f;
    DiffPair_float_0 _S4060;
    (&_S4060)->primal_0 = 0.00499999988824129f;
    (&_S4060)->differential_0 = 0.0f;
    DiffPair_float_0 _S4061;
    (&_S4061)->primal_0 = 0.98000001907348633f;
    (&_S4061)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4059, &_S4060, &_S4061, _S4056);
    float _S4062 = s_diff_dmax_T_0 / _S4033;
    float _S4063 = _S4031 * - _S4062;
    float _S4064 = _S4032 * _S4062;
    float2  _S4065 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4066;
    (&_S4066)->primal_0 = e2_3;
    (&_S4066)->differential_0 = _S4065;
    s_bwd_length_impl_1(&_S4066, _S4063);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4067;
    (&_S4067)->primal_0 = e1_7;
    (&_S4067)->differential_0 = _S4065;
    s_bwd_length_impl_1(&_S4067, _S4063);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4068;
    (&_S4068)->primal_0 = e0_7;
    (&_S4068)->differential_0 = _S4065;
    s_bwd_length_impl_1(&_S4068, _S4063);
    DiffPair_float_0 _S4069;
    (&_S4069)->primal_0 = _S4010;
    (&_S4069)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4069, _S4064);
    DiffPair_float_0 _S4070;
    (&_S4070)->primal_0 = _S4029;
    (&_S4070)->differential_0 = 0.0f;
    DiffPair_float_0 _S4071;
    (&_S4071)->primal_0 = de2_0;
    (&_S4071)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4070, &_S4071, _S4058);
    DiffPair_float_0 _S4072;
    (&_S4072)->primal_0 = de0_0;
    (&_S4072)->differential_0 = 0.0f;
    DiffPair_float_0 _S4073;
    (&_S4073)->primal_0 = de1_0;
    (&_S4073)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4072, &_S4073, _S4070.differential_0);
    float _S4074 = se_1 * _S4071.differential_0;
    float _S4075 = - _S4074;
    float _S4076 = _S4028 * _S4075;
    float _S4077 = _S4026 * _S4074;
    float2  _S4078 = make_float2 (_S4027 * _S4075, _S4025 * _S4074);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4079;
    (&_S4079)->primal_0 = e2_3;
    (&_S4079)->differential_0 = _S4065;
    s_bwd_normalize_impl_1(&_S4079, _S4078);
    float2  _S4080 = - make_float2 (_S4077, _S4076);
    float _S4081 = se_1 * _S4073.differential_0;
    float _S4082 = - _S4081;
    float _S4083 = _S4022 * _S4082;
    float _S4084 = _S4020 * _S4081;
    float2  _S4085 = make_float2 (_S4021 * _S4082, _S4019 * _S4081);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4086;
    (&_S4086)->primal_0 = e1_7;
    (&_S4086)->differential_0 = _S4065;
    s_bwd_normalize_impl_1(&_S4086, _S4085);
    float2  _S4087 = - make_float2 (_S4084, _S4083);
    float _S4088 = se_1 * _S4072.differential_0;
    float _S4089 = - _S4088;
    float _S4090 = _S4016 * _S4089;
    float _S4091 = _S4014 * _S4088;
    float2  _S4092 = make_float2 (_S4015 * _S4089, _S4013 * _S4088);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4093;
    (&_S4093)->primal_0 = e0_7;
    (&_S4093)->differential_0 = _S4065;
    s_bwd_normalize_impl_1(&_S4093, _S4092);
    float2  _S4094 = - make_float2 (_S4091, _S4090);
    float _S4095 = - _S4069.differential_0;
    float2  _S4096 = _S4066.differential_0 + _S4079.differential_0;
    float2  _S4097 = - _S4096;
    float2  _S4098 = _S4067.differential_0 + _S4086.differential_0 + make_float2 (_S4008 * _S4095, _S4006 * _S4069.differential_0);
    float2  _S4099 = - _S4098;
    float2  _S4100 = _S4068.differential_0 + _S4093.differential_0 + make_float2 (_S4007 * _S4069.differential_0, _S4009 * _S4095);
    float2  _S4101 = - _S4100;
    float2  _S4102 = make_float2 (_S4047, _S4059.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S4102;
    float2  _S4103 = _S4080 + _S4097 + _S4098;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S4103;
    float2  _S4104 = _S4087 + _S4099 + _S4100;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S4104;
    float2  _S4105 = _S4094 + _S4096 + _S4101;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S4105;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4106, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4107, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4108, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4109, float2  _S4110, float _S4111)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S4106, _S4107, _S4108, _S4109, _S4110, _S4111);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_3, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S4112 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S4112;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S4112;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S4112;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S4112;
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
    float2  _S4113 = p_4 - v0_2;
    float2  _S4114 = p_4 - v1_2;
    float2  _S4115 = p_4 - v2_2;
    float _S4116 = e0_8.x;
    float _S4117 = e1_8.y;
    float _S4118 = e0_8.y;
    float _S4119 = e1_8.x;
    float _S4120 = _S4116 * _S4117 - _S4118 * _S4119;
    float se_2 = float((F32_sign((_S4120))));
    float _S4121 = hardness_8.x;
    float _S4122 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S4113.x * _S4118 - _S4113.y * _S4116)), (se_2 * (_S4114.x * _S4117 - _S4114.y * _S4119))))), (se_2 * (_S4115.x * e2_4.y - _S4115.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S4113 - e0_8 * make_float2 (clamp_0(dot_1(_S4113, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S4114 - e1_8 * make_float2 (clamp_0(dot_1(_S4114, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S4115 - e2_4 * make_float2 (clamp_0(dot_1(_S4115, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S4120))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S4122))));
    float _S4123;
    if(a_3 <= 0.0f)
    {
        _S4123 = 0.0f;
    }
    else
    {
        _S4123 = (F32_min(((F32_pow((a_3), (_S4122)))), (0.99900001287460327f)));
    }
    return _S4121 * _S4123;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S4124, float2  _S4125)
{
    return dot_1(_S4124, _S4125);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4126, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4127, float _S4128)
{
    _d_dot_1(_S4126, _S4127, _S4128);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_5, float _s_dOut_10)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S4129 = p_5 - (*dpv0_1).primal_0;
    float _S4130 = s_primal_ctx_dot_1(_S4129, e0_9);
    float _S4131 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S4132 = _S4130 / _S4131;
    float _S4133 = _S4131 * _S4131;
    float _S4134 = s_primal_ctx_clamp_0(_S4132, 0.0f, 1.0f);
    float2  _S4135 = make_float2 (_S4134);
    float2  _S4136 = _S4129 - e0_9 * make_float2 (_S4134);
    float _S4137 = length_0(_S4136);
    float2  _S4138 = p_5 - (*dpv1_1).primal_0;
    float _S4139 = s_primal_ctx_dot_1(_S4138, e1_9);
    float _S4140 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S4141 = _S4139 / _S4140;
    float _S4142 = _S4140 * _S4140;
    float _S4143 = s_primal_ctx_clamp_0(_S4141, 0.0f, 1.0f);
    float2  _S4144 = make_float2 (_S4143);
    float2  _S4145 = _S4138 - e1_9 * make_float2 (_S4143);
    float _S4146 = length_0(_S4145);
    float2  _S4147 = p_5 - (*dpv2_1).primal_0;
    float _S4148 = s_primal_ctx_dot_1(_S4147, e2_5);
    float _S4149 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S4150 = _S4148 / _S4149;
    float _S4151 = _S4149 * _S4149;
    float _S4152 = s_primal_ctx_clamp_0(_S4150, 0.0f, 1.0f);
    float2  _S4153 = make_float2 (_S4152);
    float2  _S4154 = _S4147 - e2_5 * make_float2 (_S4152);
    float _S4155 = length_0(_S4154);
    float _S4156 = e0_9.x;
    float _S4157 = e1_9.y;
    float _S4158 = e0_9.y;
    float _S4159 = e1_9.x;
    float _S4160 = _S4156 * _S4157 - _S4158 * _S4159;
    float se_3 = float((F32_sign((_S4160))));
    float _S4161 = _S4129.x;
    float _S4162 = _S4129.y;
    float s0_0 = se_3 * (_S4161 * _S4158 - _S4162 * _S4156);
    float _S4163 = _S4138.x;
    float _S4164 = _S4138.y;
    float s1_0 = se_3 * (_S4163 * _S4157 - _S4164 * _S4159);
    float _S4165 = _S4147.x;
    float _S4166 = e2_5.y;
    float _S4167 = _S4147.y;
    float _S4168 = e2_5.x;
    float s2_0 = se_3 * (_S4165 * _S4166 - _S4167 * _S4168);
    float _S4169 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S4169, s2_0)))));
    float _S4170 = s_primal_ctx_min_0(_S4137, _S4146);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S4170, _S4155);
    float _S4171 = s_primal_ctx_abs_0(_S4160);
    float _S4172 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S4171 / _S4172;
    float _S4173 = _S4172 * _S4172;
    float _S4174 = (*dphardness_1).primal_0.x;
    float _S4175 = (*dphardness_1).primal_0.y;
    float _S4176 = dmax_1 * dmax_1;
    float _S4177 = 1.0f + dv_0 / dmax_1;
    float _S4178 = 1.0f - s_primal_ctx_clamp_0(_S4175, 0.00499999988824129f, 0.98000001907348633f);
    float _S4179 = -1.0f / _S4178;
    float _S4180 = _S4178 * _S4178;
    float _S4181 = 1.0f - s_primal_ctx_exp2_0(_S4179);
    float a_4 = 1.0f - _S4177 * _S4181;
    bool _S4182 = a_4 <= 0.0f;
    float _S4183;
    float _S4184;
    if(_S4182)
    {
        _S4183 = 0.0f;
        _S4184 = 0.0f;
    }
    else
    {
        float _S4185 = s_primal_ctx_pow_0(a_4, _S4178);
        _S4183 = s_primal_ctx_min_0(_S4185, 0.99900001287460327f);
        _S4184 = _S4185;
    }
    float _S4186 = _S4174 * _s_dOut_10;
    float _S4187 = _S4183 * _s_dOut_10;
    if(_S4182)
    {
        _S4183 = 0.0f;
        _S4184 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4188;
        (&_S4188)->primal_0 = _S4184;
        (&_S4188)->differential_0 = 0.0f;
        DiffPair_float_0 _S4189;
        (&_S4189)->primal_0 = 0.99900001287460327f;
        (&_S4189)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4188, &_S4189, _S4186);
        DiffPair_float_0 _S4190;
        (&_S4190)->primal_0 = a_4;
        (&_S4190)->differential_0 = 0.0f;
        DiffPair_float_0 _S4191;
        (&_S4191)->primal_0 = _S4178;
        (&_S4191)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4190, &_S4191, _S4188.differential_0);
        _S4183 = _S4190.differential_0;
        _S4184 = _S4191.differential_0;
    }
    float _S4192 = - _S4183;
    float _S4193 = _S4181 * _S4192;
    float _S4194 = - (_S4177 * _S4192);
    DiffPair_float_0 _S4195;
    (&_S4195)->primal_0 = _S4179;
    (&_S4195)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4195, _S4194);
    float _S4196 = - (-1.0f * - (_S4195.differential_0 / _S4180) + _S4184);
    float _S4197 = _S4193 / _S4176;
    float s_diff_dmax_T_1 = dv_0 * - _S4197;
    float s_diff_dv_T_0 = dmax_1 * _S4197;
    DiffPair_float_0 _S4198;
    (&_S4198)->primal_0 = _S4175;
    (&_S4198)->differential_0 = 0.0f;
    DiffPair_float_0 _S4199;
    (&_S4199)->primal_0 = 0.00499999988824129f;
    (&_S4199)->differential_0 = 0.0f;
    DiffPair_float_0 _S4200;
    (&_S4200)->primal_0 = 0.98000001907348633f;
    (&_S4200)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4198, &_S4199, &_S4200, _S4196);
    float _S4201 = s_diff_dmax_T_1 / _S4173;
    float _S4202 = _S4171 * - _S4201;
    float _S4203 = _S4172 * _S4201;
    float2  _S4204 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4205;
    (&_S4205)->primal_0 = e2_5;
    (&_S4205)->differential_0 = _S4204;
    s_bwd_length_impl_1(&_S4205, _S4202);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4206;
    (&_S4206)->primal_0 = e1_9;
    (&_S4206)->differential_0 = _S4204;
    s_bwd_length_impl_1(&_S4206, _S4202);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4207;
    (&_S4207)->primal_0 = e0_9;
    (&_S4207)->differential_0 = _S4204;
    s_bwd_length_impl_1(&_S4207, _S4202);
    DiffPair_float_0 _S4208;
    (&_S4208)->primal_0 = _S4160;
    (&_S4208)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4208, _S4203);
    float _S4209 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S4210;
    (&_S4210)->primal_0 = _S4170;
    (&_S4210)->differential_0 = 0.0f;
    DiffPair_float_0 _S4211;
    (&_S4211)->primal_0 = _S4155;
    (&_S4211)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4210, &_S4211, _S4209);
    DiffPair_float_0 _S4212;
    (&_S4212)->primal_0 = _S4137;
    (&_S4212)->differential_0 = 0.0f;
    DiffPair_float_0 _S4213;
    (&_S4213)->primal_0 = _S4146;
    (&_S4213)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4212, &_S4213, _S4210.differential_0);
    DiffPair_float_0 _S4214;
    (&_S4214)->primal_0 = _S4169;
    (&_S4214)->differential_0 = 0.0f;
    DiffPair_float_0 _S4215;
    (&_S4215)->primal_0 = s2_0;
    (&_S4215)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4214, &_S4215, 0.0f);
    DiffPair_float_0 _S4216;
    (&_S4216)->primal_0 = s0_0;
    (&_S4216)->differential_0 = 0.0f;
    DiffPair_float_0 _S4217;
    (&_S4217)->primal_0 = s1_0;
    (&_S4217)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4216, &_S4217, _S4214.differential_0);
    float _S4218 = se_3 * _S4215.differential_0;
    float _S4219 = - _S4218;
    float _S4220 = _S4167 * _S4219;
    float _S4221 = _S4168 * _S4219;
    float _S4222 = _S4165 * _S4218;
    float _S4223 = _S4166 * _S4218;
    float _S4224 = se_3 * _S4217.differential_0;
    float _S4225 = - _S4224;
    float _S4226 = _S4159 * _S4225;
    float _S4227 = _S4157 * _S4224;
    float _S4228 = se_3 * _S4216.differential_0;
    float _S4229 = - _S4228;
    float _S4230 = _S4156 * _S4229;
    float _S4231 = _S4158 * _S4228;
    float _S4232 = - _S4208.differential_0;
    float _S4233 = _S4164 * _S4225 + _S4158 * _S4232;
    float _S4234 = _S4161 * _S4228 + _S4159 * _S4232;
    float _S4235 = _S4163 * _S4224 + _S4156 * _S4208.differential_0;
    float _S4236 = _S4162 * _S4229 + _S4157 * _S4208.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4237;
    (&_S4237)->primal_0 = _S4154;
    (&_S4237)->differential_0 = _S4204;
    s_bwd_length_impl_1(&_S4237, _S4211.differential_0);
    float2  _S4238 = - _S4237.differential_0;
    float2  _S4239 = e2_5 * _S4238;
    float2  _S4240 = _S4153 * _S4238;
    float _S4241 = _S4239.x + _S4239.y;
    DiffPair_float_0 _S4242;
    (&_S4242)->primal_0 = _S4150;
    (&_S4242)->differential_0 = 0.0f;
    DiffPair_float_0 _S4243;
    (&_S4243)->primal_0 = 0.0f;
    (&_S4243)->differential_0 = 0.0f;
    DiffPair_float_0 _S4244;
    (&_S4244)->primal_0 = 1.0f;
    (&_S4244)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4242, &_S4243, &_S4244, _S4241);
    float _S4245 = _S4242.differential_0 / _S4151;
    float _S4246 = _S4148 * - _S4245;
    float _S4247 = _S4149 * _S4245;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4248;
    (&_S4248)->primal_0 = e2_5;
    (&_S4248)->differential_0 = _S4204;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4249;
    (&_S4249)->primal_0 = e2_5;
    (&_S4249)->differential_0 = _S4204;
    s_bwd_prop_dot_1(&_S4248, &_S4249, _S4246);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4250;
    (&_S4250)->primal_0 = _S4147;
    (&_S4250)->differential_0 = _S4204;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4251;
    (&_S4251)->primal_0 = e2_5;
    (&_S4251)->differential_0 = _S4204;
    s_bwd_prop_dot_1(&_S4250, &_S4251, _S4247);
    float2  _S4252 = - (_S4237.differential_0 + _S4250.differential_0 + make_float2 (_S4223, _S4221));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4253;
    (&_S4253)->primal_0 = _S4145;
    (&_S4253)->differential_0 = _S4204;
    s_bwd_length_impl_1(&_S4253, _S4213.differential_0);
    float2  _S4254 = - _S4253.differential_0;
    float2  _S4255 = e1_9 * _S4254;
    float2  _S4256 = _S4144 * _S4254;
    float _S4257 = _S4255.x + _S4255.y;
    DiffPair_float_0 _S4258;
    (&_S4258)->primal_0 = _S4141;
    (&_S4258)->differential_0 = 0.0f;
    DiffPair_float_0 _S4259;
    (&_S4259)->primal_0 = 0.0f;
    (&_S4259)->differential_0 = 0.0f;
    DiffPair_float_0 _S4260;
    (&_S4260)->primal_0 = 1.0f;
    (&_S4260)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4258, &_S4259, &_S4260, _S4257);
    float _S4261 = _S4258.differential_0 / _S4142;
    float _S4262 = _S4139 * - _S4261;
    float _S4263 = _S4140 * _S4261;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4264;
    (&_S4264)->primal_0 = e1_9;
    (&_S4264)->differential_0 = _S4204;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4265;
    (&_S4265)->primal_0 = e1_9;
    (&_S4265)->differential_0 = _S4204;
    s_bwd_prop_dot_1(&_S4264, &_S4265, _S4262);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4266;
    (&_S4266)->primal_0 = _S4138;
    (&_S4266)->differential_0 = _S4204;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4267;
    (&_S4267)->primal_0 = e1_9;
    (&_S4267)->differential_0 = _S4204;
    s_bwd_prop_dot_1(&_S4266, &_S4267, _S4263);
    float2  _S4268 = - (_S4253.differential_0 + _S4266.differential_0 + make_float2 (_S4227, _S4226));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4269;
    (&_S4269)->primal_0 = _S4136;
    (&_S4269)->differential_0 = _S4204;
    s_bwd_length_impl_1(&_S4269, _S4212.differential_0);
    float2  _S4270 = - _S4269.differential_0;
    float2  _S4271 = e0_9 * _S4270;
    float2  _S4272 = _S4135 * _S4270;
    float _S4273 = _S4271.x + _S4271.y;
    DiffPair_float_0 _S4274;
    (&_S4274)->primal_0 = _S4132;
    (&_S4274)->differential_0 = 0.0f;
    DiffPair_float_0 _S4275;
    (&_S4275)->primal_0 = 0.0f;
    (&_S4275)->differential_0 = 0.0f;
    DiffPair_float_0 _S4276;
    (&_S4276)->primal_0 = 1.0f;
    (&_S4276)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4274, &_S4275, &_S4276, _S4273);
    float _S4277 = _S4274.differential_0 / _S4133;
    float _S4278 = _S4130 * - _S4277;
    float _S4279 = _S4131 * _S4277;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4280;
    (&_S4280)->primal_0 = e0_9;
    (&_S4280)->differential_0 = _S4204;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4281;
    (&_S4281)->primal_0 = e0_9;
    (&_S4281)->differential_0 = _S4204;
    s_bwd_prop_dot_1(&_S4280, &_S4281, _S4278);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4282;
    (&_S4282)->primal_0 = _S4129;
    (&_S4282)->differential_0 = _S4204;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4283;
    (&_S4283)->primal_0 = e0_9;
    (&_S4283)->differential_0 = _S4204;
    s_bwd_prop_dot_1(&_S4282, &_S4283, _S4279);
    float2  _S4284 = - (_S4269.differential_0 + _S4282.differential_0 + make_float2 (_S4231, _S4230));
    float2  _S4285 = _S4205.differential_0 + _S4240 + _S4249.differential_0 + _S4248.differential_0 + _S4251.differential_0 + make_float2 (_S4220, _S4222);
    float2  _S4286 = - _S4285;
    float2  _S4287 = _S4206.differential_0 + _S4256 + _S4265.differential_0 + _S4264.differential_0 + _S4267.differential_0 + make_float2 (_S4233, _S4235);
    float2  _S4288 = - _S4287;
    float2  _S4289 = _S4207.differential_0 + _S4272 + _S4281.differential_0 + _S4280.differential_0 + _S4283.differential_0 + make_float2 (_S4236, _S4234);
    float2  _S4290 = - _S4289;
    float2  _S4291 = make_float2 (_S4187, _S4198.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S4291;
    float2  _S4292 = _S4252 + _S4286 + _S4287;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S4292;
    float2  _S4293 = _S4268 + _S4288 + _S4289;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S4293;
    float2  _S4294 = _S4284 + _S4285 + _S4290;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S4294;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4295, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4296, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4297, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4298, float2  _S4299, float _S4300)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S4295, _S4296, _S4297, _S4298, _S4299, _S4300);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_6, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S4301 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S4301;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S4301;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S4301;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S4301;
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
    float _S4302 = 0.3333333432674408f * dpdepth_1;
    float3  _S4303 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S4304 = make_float3 (0.0f);
    float3  _S4305 = _S4304;
    *&((&_S4305)->z) = _S4302;
    *&((&_S4305)->y) = _S4302;
    *&((&_S4305)->x) = _S4302;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S4305;
    FixedArray<float3 , 3>  _S4306;
    _S4306[int(0)] = _S4304;
    _S4306[int(1)] = _S4304;
    _S4306[int(2)] = _S4304;
    _S4306[int(2)] = _S4303;
    _S4306[int(1)] = _S4303;
    _S4306[int(0)] = _S4303;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S4306;
    float2  _S4307 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S4307;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S4307;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S4307;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4308, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4309, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4310, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4311, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4312, float2  _S4313, float3  _S4314, float _S4315)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S4308, _S4309, _S4310, _S4311, _S4312, _S4313, _S4314, _S4315);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_9, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S4316 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S4316;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S4316;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S4316;
    float3  _S4317 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4318 = { _S4317, _S4317, _S4317 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S4318;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S4317;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_9, v_color_0, v_depth_8);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_10, FixedArray<float3 , 16>  * sh_coeffs_21, FixedArray<float3 , 2>  * ch_coeffs_6, Matrix<float, 3, 3>  R_25, float3  t_24, float fx_30, float fy_30, float cx_25, float cy_25, FixedArray<float, 10>  * dist_coeffs_37, uint image_width_21, uint image_height_21, float near_plane_14, float far_plane_14, int4  * aabb_xyxy_14, float * depth_16, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_4)
{
    for(;;)
    {
        float3  mean_c_21 = mul_0(R_25, mean_26) + t_24;
        float _S4319 = mean_c_21.z;
        bool _S4320;
        if(_S4319 < near_plane_14)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = _S4319 > far_plane_14;
        }
        if(_S4320)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4321 = scale_26.x;
        float sx_5 = (F32_exp((_S4321)));
        float _S4322 = scale_26.y;
        float sy_5 = (F32_exp((_S4322)));
        float sz_7 = scale_26.z - 0.5f * (_S4321 + _S4322);
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
        Matrix<float, 3, 3>  _S4323 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_48), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_48), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S4323, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S4323, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S4323, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_6 = mul_0(R_25, vert0_3) + t_24;
        float3  vert1_c_6 = mul_0(R_25, vert1_3) + t_24;
        float3  vert2_c_6 = mul_0(R_25, vert2_3) + t_24;
        float _S4324 = vert0_c_6.z;
        float _S4325 = vert1_c_6.z;
        float _S4326 = vert2_c_6.z;
        if(_S4324 < near_plane_14)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = _S4324 > far_plane_14;
        }
        if(_S4320)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = _S4325 < near_plane_14;
        }
        if(_S4320)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = _S4325 > far_plane_14;
        }
        if(_S4320)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = _S4326 < near_plane_14;
        }
        if(_S4320)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = _S4326 > far_plane_14;
        }
        if(_S4320)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_4;
        for(;;)
        {
            float2  uv0_5 = float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S4324);
            if(_S4324 < 0.0f)
            {
                _S4320 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4327 = camera_distortion_jac_0(uv0_5, dist_coeffs_37);
                _S4320 = !((F32_min((determinant_0(_S4327)), ((F32_min((_S4327.rows[int(0)].x), (_S4327.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4320)
            {
                uv0_4 = uv0_5;
                break;
            }
            float u_62 = uv0_5.x;
            float v_62 = uv0_5.y;
            float r2_62 = u_62 * u_62 + v_62 * v_62;
            float2  _S4328 = uv0_5 * make_float2 (1.0f + r2_62 * ((*dist_coeffs_37)[int(0)] + r2_62 * ((*dist_coeffs_37)[int(1)] + r2_62 * ((*dist_coeffs_37)[int(2)] + r2_62 * (*dist_coeffs_37)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_37)[int(4)] * u_62 * v_62 + (*dist_coeffs_37)[int(5)] * (r2_62 + 2.0f * u_62 * u_62) + (*dist_coeffs_37)[int(6)] * r2_62, 2.0f * (*dist_coeffs_37)[int(5)] * u_62 * v_62 + (*dist_coeffs_37)[int(4)] * (r2_62 + 2.0f * v_62 * v_62) + (*dist_coeffs_37)[int(7)] * r2_62);
            float2  _S4329 = _S4328 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4328.x + (*dist_coeffs_37)[int(9)] * _S4328.y, 0.0f);
            uv0_4 = make_float2 (fx_30 * _S4329.x + cx_25, fy_30 * _S4329.y + cy_25);
            break;
        }
        float2  uv1_4;
        bool all_valid_12 = true & (!_S4320);
        for(;;)
        {
            float2  uv1_5 = float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S4325);
            if(_S4325 < 0.0f)
            {
                _S4320 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4330 = camera_distortion_jac_0(uv1_5, dist_coeffs_37);
                _S4320 = !((F32_min((determinant_0(_S4330)), ((F32_min((_S4330.rows[int(0)].x), (_S4330.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4320)
            {
                uv1_4 = uv1_5;
                break;
            }
            float u_63 = uv1_5.x;
            float v_63 = uv1_5.y;
            float r2_63 = u_63 * u_63 + v_63 * v_63;
            float2  _S4331 = uv1_5 * make_float2 (1.0f + r2_63 * ((*dist_coeffs_37)[int(0)] + r2_63 * ((*dist_coeffs_37)[int(1)] + r2_63 * ((*dist_coeffs_37)[int(2)] + r2_63 * (*dist_coeffs_37)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_37)[int(4)] * u_63 * v_63 + (*dist_coeffs_37)[int(5)] * (r2_63 + 2.0f * u_63 * u_63) + (*dist_coeffs_37)[int(6)] * r2_63, 2.0f * (*dist_coeffs_37)[int(5)] * u_63 * v_63 + (*dist_coeffs_37)[int(4)] * (r2_63 + 2.0f * v_63 * v_63) + (*dist_coeffs_37)[int(7)] * r2_63);
            float2  _S4332 = _S4331 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4331.x + (*dist_coeffs_37)[int(9)] * _S4331.y, 0.0f);
            uv1_4 = make_float2 (fx_30 * _S4332.x + cx_25, fy_30 * _S4332.y + cy_25);
            break;
        }
        float2  uv2_4;
        bool all_valid_13 = all_valid_12 & (!_S4320);
        for(;;)
        {
            float2  uv2_5 = float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S4326);
            if(_S4326 < 0.0f)
            {
                _S4320 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4333 = camera_distortion_jac_0(uv2_5, dist_coeffs_37);
                _S4320 = !((F32_min((determinant_0(_S4333)), ((F32_min((_S4333.rows[int(0)].x), (_S4333.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4320)
            {
                uv2_4 = uv2_5;
                break;
            }
            float u_64 = uv2_5.x;
            float v_64 = uv2_5.y;
            float r2_64 = u_64 * u_64 + v_64 * v_64;
            float2  _S4334 = uv2_5 * make_float2 (1.0f + r2_64 * ((*dist_coeffs_37)[int(0)] + r2_64 * ((*dist_coeffs_37)[int(1)] + r2_64 * ((*dist_coeffs_37)[int(2)] + r2_64 * (*dist_coeffs_37)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_37)[int(4)] * u_64 * v_64 + (*dist_coeffs_37)[int(5)] * (r2_64 + 2.0f * u_64 * u_64) + (*dist_coeffs_37)[int(6)] * r2_64, 2.0f * (*dist_coeffs_37)[int(5)] * u_64 * v_64 + (*dist_coeffs_37)[int(4)] * (r2_64 + 2.0f * v_64 * v_64) + (*dist_coeffs_37)[int(7)] * r2_64);
            float2  _S4335 = _S4334 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4334.x + (*dist_coeffs_37)[int(9)] * _S4334.y, 0.0f);
            uv2_4 = make_float2 (fx_30 * _S4335.x + cx_25, fy_30 * _S4335.y + cy_25);
            break;
        }
        if(!(all_valid_13 & (!_S4320)))
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_10 = uv1_4 - uv0_4;
        float2  e1_10 = uv2_4 - uv1_4;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(uv0_4 - uv2_4)));
        float _S4336 = uv0_4.x;
        float _S4337 = uv1_4.x;
        float _S4338 = uv2_4.x;
        float xmax_7 = (F32_max(((F32_max((_S4336), (_S4337)))), (_S4338))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S4336), (_S4337)))), (_S4338))) - offset_4;
        float _S4339 = uv0_4.y;
        float _S4340 = uv1_4.y;
        float _S4341 = uv2_4.y;
        float ymax_7 = (F32_max(((F32_max((_S4339), (_S4340)))), (_S4341))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S4339), (_S4340)))), (_S4341))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = xmin_7 >= float(image_width_21);
        }
        if(_S4320)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = ymax_7 <= 0.0f;
        }
        if(_S4320)
        {
            _S4320 = true;
        }
        else
        {
            _S4320 = ymin_7 >= float(image_height_21);
        }
        if(_S4320)
        {
            _S4320 = true;
        }
        else
        {
            if(_S4319 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S4320 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S4320 = false;
                }
                if(_S4320)
                {
                    _S4320 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S4320 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S4320 = false;
                    }
                }
            }
            else
            {
                _S4320 = false;
            }
        }
        if(_S4320)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4342 = mean_26 - - mul_0(transpose_0(R_25), t_24);
        float _S4343 = _S4342.x;
        float _S4344 = _S4342.y;
        float _S4345 = _S4342.z;
        float norm_14 = (F32_sqrt((_S4343 * _S4343 + _S4344 * _S4344 + _S4345 * _S4345)));
        float x_54 = _S4343 / norm_14;
        float y_24 = _S4344 / norm_14;
        float z_21 = _S4345 / norm_14;
        float z2_49 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_54 * x_54 - y_24 * y_24;
        float fS1_21 = 2.0f * x_54 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_49 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_54) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_49 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_54) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_54 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_49 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_54) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_54 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S4346 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S4346);
        float3  _S4347 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S4348 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S4347 + _S4348 + make_float3 (0.5f), _S4346);
        (*rgbs_0)[int(2)] = max_0(_S4347 - _S4348 + make_float3 (0.5f), _S4346);
        (*verts_0)[int(0)] = vert0_3;
        (*verts_0)[int(1)] = vert1_3;
        (*verts_0)[int(2)] = vert2_3;
        float3  _S4349 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S4349 * make_float3 (float(- (F32_sign((dot_0(_S4349, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_26, float3  t_25, float fx_31, float fy_31, float cx_26, float cy_26, FixedArray<float, 10>  * dist_coeffs_38, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    bool _S4350;
    bool _S4351;
    bool _S4352;
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_26, mean_27) + t_25;
        float _S4353 = length_1(mean_c_22);
        bool _S4354;
        if(_S4353 < near_plane_15)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = _S4353 > far_plane_15;
        }
        if(_S4354)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4355 = scale_27.x;
        float sx_6 = (F32_exp((_S4355)));
        float _S4356 = scale_27.y;
        float sy_6 = (F32_exp((_S4356)));
        float sz_8 = scale_27.z - 0.5f * (_S4355 + _S4356);
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
        Matrix<float, 3, 3>  _S4357 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_50), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_50), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
        float3  vert0_4 = mul_0(_S4357, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
        float3  vert1_4 = mul_0(_S4357, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
        float3  vert2_4 = mul_0(_S4357, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
        float3  vert0_c_7 = mul_0(R_26, vert0_4) + t_25;
        float3  vert1_c_7 = mul_0(R_26, vert1_4) + t_25;
        float3  vert2_c_7 = mul_0(R_26, vert2_4) + t_25;
        float _S4358 = length_1(vert0_c_7);
        float _S4359 = length_1(vert1_c_7);
        float _S4360 = length_1(vert2_c_7);
        if(_S4358 < near_plane_15)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = _S4358 > far_plane_15;
        }
        if(_S4354)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = _S4359 < near_plane_15;
        }
        if(_S4354)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = _S4359 > far_plane_15;
        }
        if(_S4354)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = _S4360 < near_plane_15;
        }
        if(_S4354)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = _S4360 > far_plane_15;
        }
        if(_S4354)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_6;
        float k_10;
        for(;;)
        {
            float2  _S4361 = float2 {vert0_c_7.x, vert0_c_7.y};
            float r_30 = length_0(_S4361);
            float _S4362 = vert0_c_7.z;
            float theta_26 = (F32_atan2((r_30), (_S4362)));
            if(theta_26 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_26 * theta_26 / 3.0f) / _S4362;
            }
            else
            {
                k_10 = theta_26 / r_30;
            }
            float2  uv0_7 = _S4361 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4363 = camera_distortion_jac_0(uv0_7, dist_coeffs_38);
            bool _S4364 = !((F32_min((determinant_0(_S4363)), ((F32_min((_S4363.rows[int(0)].x), (_S4363.rows[int(1)].y)))))) > 0.0f);
            _S4350 = _S4364;
            if(_S4364)
            {
                uv0_6 = uv0_7;
                break;
            }
            float u_65 = uv0_7.x;
            float v_65 = uv0_7.y;
            float r2_65 = u_65 * u_65 + v_65 * v_65;
            float2  _S4365 = uv0_7 * make_float2 (1.0f + r2_65 * ((*dist_coeffs_38)[int(0)] + r2_65 * ((*dist_coeffs_38)[int(1)] + r2_65 * ((*dist_coeffs_38)[int(2)] + r2_65 * (*dist_coeffs_38)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_38)[int(4)] * u_65 * v_65 + (*dist_coeffs_38)[int(5)] * (r2_65 + 2.0f * u_65 * u_65) + (*dist_coeffs_38)[int(6)] * r2_65, 2.0f * (*dist_coeffs_38)[int(5)] * u_65 * v_65 + (*dist_coeffs_38)[int(4)] * (r2_65 + 2.0f * v_65 * v_65) + (*dist_coeffs_38)[int(7)] * r2_65);
            float2  _S4366 = _S4365 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4365.x + (*dist_coeffs_38)[int(9)] * _S4365.y, 0.0f);
            uv0_6 = make_float2 (fx_31 * _S4366.x + cx_26, fy_31 * _S4366.y + cy_26);
            break;
        }
        float2  uv1_6;
        bool all_valid_14 = true & (!_S4350);
        for(;;)
        {
            float2  _S4367 = float2 {vert1_c_7.x, vert1_c_7.y};
            float r_31 = length_0(_S4367);
            float _S4368 = vert1_c_7.z;
            float theta_27 = (F32_atan2((r_31), (_S4368)));
            if(theta_27 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_27 * theta_27 / 3.0f) / _S4368;
            }
            else
            {
                k_10 = theta_27 / r_31;
            }
            float2  uv1_7 = _S4367 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4369 = camera_distortion_jac_0(uv1_7, dist_coeffs_38);
            bool _S4370 = !((F32_min((determinant_0(_S4369)), ((F32_min((_S4369.rows[int(0)].x), (_S4369.rows[int(1)].y)))))) > 0.0f);
            _S4351 = _S4370;
            if(_S4370)
            {
                uv1_6 = uv1_7;
                break;
            }
            float u_66 = uv1_7.x;
            float v_66 = uv1_7.y;
            float r2_66 = u_66 * u_66 + v_66 * v_66;
            float2  _S4371 = uv1_7 * make_float2 (1.0f + r2_66 * ((*dist_coeffs_38)[int(0)] + r2_66 * ((*dist_coeffs_38)[int(1)] + r2_66 * ((*dist_coeffs_38)[int(2)] + r2_66 * (*dist_coeffs_38)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_38)[int(4)] * u_66 * v_66 + (*dist_coeffs_38)[int(5)] * (r2_66 + 2.0f * u_66 * u_66) + (*dist_coeffs_38)[int(6)] * r2_66, 2.0f * (*dist_coeffs_38)[int(5)] * u_66 * v_66 + (*dist_coeffs_38)[int(4)] * (r2_66 + 2.0f * v_66 * v_66) + (*dist_coeffs_38)[int(7)] * r2_66);
            float2  _S4372 = _S4371 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4371.x + (*dist_coeffs_38)[int(9)] * _S4371.y, 0.0f);
            uv1_6 = make_float2 (fx_31 * _S4372.x + cx_26, fy_31 * _S4372.y + cy_26);
            break;
        }
        float2  uv2_6;
        bool all_valid_15 = all_valid_14 & (!_S4351);
        for(;;)
        {
            float2  _S4373 = float2 {vert2_c_7.x, vert2_c_7.y};
            float r_32 = length_0(_S4373);
            float _S4374 = vert2_c_7.z;
            float theta_28 = (F32_atan2((r_32), (_S4374)));
            if(theta_28 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_28 * theta_28 / 3.0f) / _S4374;
            }
            else
            {
                k_10 = theta_28 / r_32;
            }
            float2  uv2_7 = _S4373 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4375 = camera_distortion_jac_0(uv2_7, dist_coeffs_38);
            bool _S4376 = !((F32_min((determinant_0(_S4375)), ((F32_min((_S4375.rows[int(0)].x), (_S4375.rows[int(1)].y)))))) > 0.0f);
            _S4352 = _S4376;
            if(_S4376)
            {
                uv2_6 = uv2_7;
                break;
            }
            float u_67 = uv2_7.x;
            float v_67 = uv2_7.y;
            float r2_67 = u_67 * u_67 + v_67 * v_67;
            float2  _S4377 = uv2_7 * make_float2 (1.0f + r2_67 * ((*dist_coeffs_38)[int(0)] + r2_67 * ((*dist_coeffs_38)[int(1)] + r2_67 * ((*dist_coeffs_38)[int(2)] + r2_67 * (*dist_coeffs_38)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_38)[int(4)] * u_67 * v_67 + (*dist_coeffs_38)[int(5)] * (r2_67 + 2.0f * u_67 * u_67) + (*dist_coeffs_38)[int(6)] * r2_67, 2.0f * (*dist_coeffs_38)[int(5)] * u_67 * v_67 + (*dist_coeffs_38)[int(4)] * (r2_67 + 2.0f * v_67 * v_67) + (*dist_coeffs_38)[int(7)] * r2_67);
            float2  _S4378 = _S4377 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4377.x + (*dist_coeffs_38)[int(9)] * _S4377.y, 0.0f);
            uv2_6 = make_float2 (fx_31 * _S4378.x + cx_26, fy_31 * _S4378.y + cy_26);
            break;
        }
        if(!(all_valid_15 & (!_S4352)))
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_11 = uv1_6 - uv0_6;
        float2  e1_11 = uv2_6 - uv1_6;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(uv0_6 - uv2_6)));
        float _S4379 = uv0_6.x;
        float _S4380 = uv1_6.x;
        float _S4381 = uv2_6.x;
        float xmax_8 = (F32_max(((F32_max((_S4379), (_S4380)))), (_S4381))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S4379), (_S4380)))), (_S4381))) - offset_5;
        float _S4382 = uv0_6.y;
        float _S4383 = uv1_6.y;
        float _S4384 = uv2_6.y;
        float ymax_8 = (F32_max(((F32_max((_S4382), (_S4383)))), (_S4384))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S4382), (_S4383)))), (_S4384))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = xmin_8 >= float(image_width_22);
        }
        if(_S4354)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = ymax_8 <= 0.0f;
        }
        if(_S4354)
        {
            _S4354 = true;
        }
        else
        {
            _S4354 = ymin_8 >= float(image_height_22);
        }
        if(_S4354)
        {
            _S4354 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S4354 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S4354 = false;
                }
                if(_S4354)
                {
                    _S4354 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S4354 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S4354 = false;
                    }
                }
            }
            else
            {
                _S4354 = false;
            }
        }
        if(_S4354)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4385 = mean_27 - - mul_0(transpose_0(R_26), t_25);
        float _S4386 = _S4385.x;
        float _S4387 = _S4385.y;
        float _S4388 = _S4385.z;
        float norm_15 = (F32_sqrt((_S4386 * _S4386 + _S4387 * _S4387 + _S4388 * _S4388)));
        float x_56 = _S4386 / norm_15;
        float y_25 = _S4387 / norm_15;
        float z_22 = _S4388 / norm_15;
        float z2_51 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_56 * x_56 - y_25 * y_25;
        float fS1_22 = 2.0f * x_56 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_51 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_56) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_51 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_56) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_56 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_51 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_56) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_56 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S4389 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S4389);
        float3  _S4390 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S4391 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S4390 + _S4391 + make_float3 (0.5f), _S4389);
        (*rgbs_1)[int(2)] = max_0(_S4390 - _S4391 + make_float3 (0.5f), _S4389);
        (*verts_1)[int(0)] = vert0_4;
        (*verts_1)[int(1)] = vert1_4;
        (*verts_1)[int(2)] = vert2_4;
        float3  _S4392 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S4392 * make_float3 (float(- (F32_sign((dot_0(_S4392, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_27, float3  t_26, float fx_32, float fy_32, float cx_27, float cy_27, FixedArray<float, 10>  * dist_coeffs_39, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_27, mean_28) + t_26;
    float _S4393 = scale_28.x;
    float sx_7 = (F32_exp((_S4393)));
    float _S4394 = scale_28.y;
    float sy_7 = (F32_exp((_S4394)));
    float sz_9 = scale_28.z - 0.5f * (_S4393 + _S4394);
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
    Matrix<float, 3, 3>  _S4395 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_52), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_52), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S4395, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S4395, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S4395, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_8 = mul_0(R_27, vert0_5) + t_26;
    float3  vert1_c_8 = mul_0(R_27, vert1_5) + t_26;
    float3  vert2_c_8 = mul_0(R_27, vert2_5) + t_26;
    float2  _S4396 = float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z);
    float u_68 = _S4396.x;
    float v_68 = _S4396.y;
    float r2_68 = u_68 * u_68 + v_68 * v_68;
    float _S4397 = 2.0f * (*dist_coeffs_39)[int(4)];
    float _S4398 = 2.0f * (*dist_coeffs_39)[int(5)];
    float2  _S4399 = _S4396 * make_float2 (1.0f + r2_68 * ((*dist_coeffs_39)[int(0)] + r2_68 * ((*dist_coeffs_39)[int(1)] + r2_68 * ((*dist_coeffs_39)[int(2)] + r2_68 * (*dist_coeffs_39)[int(3)])))) + make_float2 (_S4397 * u_68 * v_68 + (*dist_coeffs_39)[int(5)] * (r2_68 + 2.0f * u_68 * u_68) + (*dist_coeffs_39)[int(6)] * r2_68, _S4398 * u_68 * v_68 + (*dist_coeffs_39)[int(4)] * (r2_68 + 2.0f * v_68 * v_68) + (*dist_coeffs_39)[int(7)] * r2_68);
    float2  _S4400 = _S4399 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4399.x + (*dist_coeffs_39)[int(9)] * _S4399.y, 0.0f);
    float _S4401 = fx_32 * _S4400.x + cx_27;
    float _S4402 = fy_32 * _S4400.y + cy_27;
    float2  uv0_8 = make_float2 (_S4401, _S4402);
    float2  _S4403 = float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z);
    float u_69 = _S4403.x;
    float v_69 = _S4403.y;
    float r2_69 = u_69 * u_69 + v_69 * v_69;
    float2  _S4404 = _S4403 * make_float2 (1.0f + r2_69 * ((*dist_coeffs_39)[int(0)] + r2_69 * ((*dist_coeffs_39)[int(1)] + r2_69 * ((*dist_coeffs_39)[int(2)] + r2_69 * (*dist_coeffs_39)[int(3)])))) + make_float2 (_S4397 * u_69 * v_69 + (*dist_coeffs_39)[int(5)] * (r2_69 + 2.0f * u_69 * u_69) + (*dist_coeffs_39)[int(6)] * r2_69, _S4398 * u_69 * v_69 + (*dist_coeffs_39)[int(4)] * (r2_69 + 2.0f * v_69 * v_69) + (*dist_coeffs_39)[int(7)] * r2_69);
    float2  _S4405 = _S4404 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4404.x + (*dist_coeffs_39)[int(9)] * _S4404.y, 0.0f);
    float _S4406 = fx_32 * _S4405.x + cx_27;
    float _S4407 = fy_32 * _S4405.y + cy_27;
    float2  uv1_8 = make_float2 (_S4406, _S4407);
    float2  _S4408 = float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z);
    float u_70 = _S4408.x;
    float v_70 = _S4408.y;
    float r2_70 = u_70 * u_70 + v_70 * v_70;
    float2  _S4409 = _S4408 * make_float2 (1.0f + r2_70 * ((*dist_coeffs_39)[int(0)] + r2_70 * ((*dist_coeffs_39)[int(1)] + r2_70 * ((*dist_coeffs_39)[int(2)] + r2_70 * (*dist_coeffs_39)[int(3)])))) + make_float2 (_S4397 * u_70 * v_70 + (*dist_coeffs_39)[int(5)] * (r2_70 + 2.0f * u_70 * u_70) + (*dist_coeffs_39)[int(6)] * r2_70, _S4398 * u_70 * v_70 + (*dist_coeffs_39)[int(4)] * (r2_70 + 2.0f * v_70 * v_70) + (*dist_coeffs_39)[int(7)] * r2_70);
    float2  _S4410 = _S4409 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4409.x + (*dist_coeffs_39)[int(9)] * _S4409.y, 0.0f);
    float _S4411 = fx_32 * _S4410.x + cx_27;
    float _S4412 = fy_32 * _S4410.y + cy_27;
    float2  uv2_8 = make_float2 (_S4411, _S4412);
    float2  e0_12 = uv1_8 - uv0_8;
    float2  e1_12 = uv2_8 - uv1_8;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(uv0_8 - uv2_8)));
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4401), (_S4406)))), (_S4411))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S4402), (_S4407)))), (_S4412))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4401), (_S4406)))), (_S4411))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4402), (_S4407)))), (_S4412))) + offset_6)))));
    *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4413 = mean_28 - - mul_0(transpose_0(R_27), t_26);
    float _S4414 = _S4413.x;
    float _S4415 = _S4413.y;
    float _S4416 = _S4413.z;
    float norm_16 = (F32_sqrt((_S4414 * _S4414 + _S4415 * _S4415 + _S4416 * _S4416)));
    float x_58 = _S4414 / norm_16;
    float y_26 = _S4415 / norm_16;
    float z_23 = _S4416 / norm_16;
    float z2_53 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_58 * x_58 - y_26 * y_26;
    float fS1_23 = 2.0f * x_58 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_53 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_58) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_53 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_58) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_58 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_53 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_58) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_58 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S4417 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S4417);
    float3  _S4418 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S4419 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S4418 + _S4419 + make_float3 (0.5f), _S4417);
    (*rgbs_2)[int(2)] = max_0(_S4418 - _S4419 + make_float3 (0.5f), _S4417);
    (*verts_2)[int(0)] = vert0_5;
    (*verts_2)[int(1)] = vert1_5;
    (*verts_2)[int(2)] = vert2_5;
    float3  _S4420 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S4420 * make_float3 (float(- (F32_sign((dot_0(_S4420, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_28, float3  t_27, float fx_33, float fy_33, float cx_28, float cy_28, FixedArray<float, 10>  * dist_coeffs_40, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_28, mean_29) + t_27;
    float _S4421 = scale_29.x;
    float sx_8 = (F32_exp((_S4421)));
    float _S4422 = scale_29.y;
    float sy_8 = (F32_exp((_S4422)));
    float sz_10 = scale_29.z - 0.5f * (_S4421 + _S4422);
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
    Matrix<float, 3, 3>  _S4423 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_54), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_54), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  vert0_6 = mul_0(_S4423, make_float3 (sx_8, 0.0f, 0.0f)) + mean_29;
    float3  vert1_6 = mul_0(_S4423, make_float3 (sx_8 * (-0.5f + sz_10), sy_8, 0.0f)) + mean_29;
    float3  vert2_6 = mul_0(_S4423, make_float3 (sx_8 * (-0.5f - sz_10), - sy_8, 0.0f)) + mean_29;
    float3  vert0_c_9 = mul_0(R_28, vert0_6) + t_27;
    float3  vert1_c_9 = mul_0(R_28, vert1_6) + t_27;
    float3  vert2_c_9 = mul_0(R_28, vert2_6) + t_27;
    float2  _S4424 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_33 = length_0(_S4424);
    float _S4425 = vert0_c_9.z;
    float theta_29 = (F32_atan2((r_33), (_S4425)));
    float k_11;
    if(theta_29 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_29 * theta_29 / 3.0f) / _S4425;
    }
    else
    {
        k_11 = theta_29 / r_33;
    }
    float2  _S4426 = _S4424 * make_float2 (k_11);
    float u_71 = _S4426.x;
    float v_71 = _S4426.y;
    float r2_71 = u_71 * u_71 + v_71 * v_71;
    float _S4427 = 2.0f * (*dist_coeffs_40)[int(4)];
    float _S4428 = 2.0f * (*dist_coeffs_40)[int(5)];
    float2  _S4429 = _S4426 * make_float2 (1.0f + r2_71 * ((*dist_coeffs_40)[int(0)] + r2_71 * ((*dist_coeffs_40)[int(1)] + r2_71 * ((*dist_coeffs_40)[int(2)] + r2_71 * (*dist_coeffs_40)[int(3)])))) + make_float2 (_S4427 * u_71 * v_71 + (*dist_coeffs_40)[int(5)] * (r2_71 + 2.0f * u_71 * u_71) + (*dist_coeffs_40)[int(6)] * r2_71, _S4428 * u_71 * v_71 + (*dist_coeffs_40)[int(4)] * (r2_71 + 2.0f * v_71 * v_71) + (*dist_coeffs_40)[int(7)] * r2_71);
    float2  _S4430 = _S4429 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4429.x + (*dist_coeffs_40)[int(9)] * _S4429.y, 0.0f);
    float _S4431 = fx_33 * _S4430.x + cx_28;
    float _S4432 = fy_33 * _S4430.y + cy_28;
    float2  uv0_9 = make_float2 (_S4431, _S4432);
    float2  _S4433 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_34 = length_0(_S4433);
    float _S4434 = vert1_c_9.z;
    float theta_30 = (F32_atan2((r_34), (_S4434)));
    if(theta_30 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_30 * theta_30 / 3.0f) / _S4434;
    }
    else
    {
        k_11 = theta_30 / r_34;
    }
    float2  _S4435 = _S4433 * make_float2 (k_11);
    float u_72 = _S4435.x;
    float v_72 = _S4435.y;
    float r2_72 = u_72 * u_72 + v_72 * v_72;
    float2  _S4436 = _S4435 * make_float2 (1.0f + r2_72 * ((*dist_coeffs_40)[int(0)] + r2_72 * ((*dist_coeffs_40)[int(1)] + r2_72 * ((*dist_coeffs_40)[int(2)] + r2_72 * (*dist_coeffs_40)[int(3)])))) + make_float2 (_S4427 * u_72 * v_72 + (*dist_coeffs_40)[int(5)] * (r2_72 + 2.0f * u_72 * u_72) + (*dist_coeffs_40)[int(6)] * r2_72, _S4428 * u_72 * v_72 + (*dist_coeffs_40)[int(4)] * (r2_72 + 2.0f * v_72 * v_72) + (*dist_coeffs_40)[int(7)] * r2_72);
    float2  _S4437 = _S4436 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4436.x + (*dist_coeffs_40)[int(9)] * _S4436.y, 0.0f);
    float _S4438 = fx_33 * _S4437.x + cx_28;
    float _S4439 = fy_33 * _S4437.y + cy_28;
    float2  uv1_9 = make_float2 (_S4438, _S4439);
    float2  _S4440 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_35 = length_0(_S4440);
    float _S4441 = vert2_c_9.z;
    float theta_31 = (F32_atan2((r_35), (_S4441)));
    if(theta_31 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_31 * theta_31 / 3.0f) / _S4441;
    }
    else
    {
        k_11 = theta_31 / r_35;
    }
    float2  _S4442 = _S4440 * make_float2 (k_11);
    float u_73 = _S4442.x;
    float v_73 = _S4442.y;
    float r2_73 = u_73 * u_73 + v_73 * v_73;
    float2  _S4443 = _S4442 * make_float2 (1.0f + r2_73 * ((*dist_coeffs_40)[int(0)] + r2_73 * ((*dist_coeffs_40)[int(1)] + r2_73 * ((*dist_coeffs_40)[int(2)] + r2_73 * (*dist_coeffs_40)[int(3)])))) + make_float2 (_S4427 * u_73 * v_73 + (*dist_coeffs_40)[int(5)] * (r2_73 + 2.0f * u_73 * u_73) + (*dist_coeffs_40)[int(6)] * r2_73, _S4428 * u_73 * v_73 + (*dist_coeffs_40)[int(4)] * (r2_73 + 2.0f * v_73 * v_73) + (*dist_coeffs_40)[int(7)] * r2_73);
    float2  _S4444 = _S4443 + make_float2 ((*dist_coeffs_40)[int(8)] * _S4443.x + (*dist_coeffs_40)[int(9)] * _S4443.y, 0.0f);
    float _S4445 = fx_33 * _S4444.x + cx_28;
    float _S4446 = fy_33 * _S4444.y + cy_28;
    float2  uv2_9 = make_float2 (_S4445, _S4446);
    float2  e0_13 = uv1_9 - uv0_9;
    float2  e1_13 = uv2_9 - uv1_9;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(uv0_9 - uv2_9)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4431), (_S4438)))), (_S4445))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S4432), (_S4439)))), (_S4446))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4431), (_S4438)))), (_S4445))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4432), (_S4439)))), (_S4446))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4447 = mean_29 - - mul_0(transpose_0(R_28), t_27);
    float _S4448 = _S4447.x;
    float _S4449 = _S4447.y;
    float _S4450 = _S4447.z;
    float norm_17 = (F32_sqrt((_S4448 * _S4448 + _S4449 * _S4449 + _S4450 * _S4450)));
    float x_60 = _S4448 / norm_17;
    float y_27 = _S4449 / norm_17;
    float z_24 = _S4450 / norm_17;
    float z2_55 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_60 * x_60 - y_27 * y_27;
    float fS1_24 = 2.0f * x_60 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_55 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_60) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_55 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_60) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_60 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_55 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_60) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_60 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S4451 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S4451);
    float3  _S4452 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S4453 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S4452 + _S4453 + make_float3 (0.5f), _S4451);
    (*rgbs_3)[int(2)] = max_0(_S4452 - _S4453 + make_float3 (0.5f), _S4451);
    (*verts_3)[int(0)] = vert0_6;
    (*verts_3)[int(1)] = vert1_6;
    (*verts_3)[int(2)] = vert2_6;
    float3  _S4454 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S4454 * make_float3 (float(- (F32_sign((dot_0(_S4454, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_29, float3  t_28, float fx_34, float fy_34, float cx_29, float cy_29, FixedArray<float, 10>  * dist_coeffs_41, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_1(R_29, mean_30) + t_28;
    float _S4455 = scale_30.x;
    float _S4456 = s_primal_ctx_exp_1(_S4455);
    float _S4457 = scale_30.y;
    float _S4458 = s_primal_ctx_exp_1(_S4457);
    float sz_11 = scale_30.z - 0.5f * (_S4455 + _S4457);
    float _S4459 = quat_31.y;
    float x2_31 = _S4459 * _S4459;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_56 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4460 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_56), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_56), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4461 = make_float3 (_S4456, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_1(_S4460, _S4461) + mean_30;
    float _S4462 = -0.5f + sz_11;
    float3  _S4463 = make_float3 (_S4456 * _S4462, _S4458, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_1(_S4460, _S4463) + mean_30;
    float _S4464 = -0.5f - sz_11;
    float3  _S4465 = make_float3 (_S4456 * _S4464, - _S4458, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_1(_S4460, _S4465) + mean_30;
    float3  vert0_c_10 = s_primal_ctx_mul_1(R_29, vert0_7) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_1(R_29, vert1_7) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_1(R_29, vert2_7) + t_28;
    float2  _S4466 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S4467 = vert0_c_10.z;
    float2  _S4468 = make_float2 (_S4467);
    float2  _S4469 = _S4466 / make_float2 (_S4467);
    float2  _S4470 = make_float2 (_S4467 * _S4467);
    float u_74 = _S4469.x;
    float v_74 = _S4469.y;
    float r2_74 = u_74 * u_74 + v_74 * v_74;
    float _S4471 = (*dist_coeffs_41)[int(2)] + r2_74 * (*dist_coeffs_41)[int(3)];
    float _S4472 = (*dist_coeffs_41)[int(1)] + r2_74 * _S4471;
    float _S4473 = (*dist_coeffs_41)[int(0)] + r2_74 * _S4472;
    float radial_7 = 1.0f + r2_74 * _S4473;
    float2  _S4474 = make_float2 (radial_7);
    float _S4475 = 2.0f * (*dist_coeffs_41)[int(4)];
    float _S4476 = _S4475 * u_74;
    float _S4477 = 2.0f * u_74;
    float _S4478 = 2.0f * (*dist_coeffs_41)[int(5)];
    float _S4479 = _S4478 * u_74;
    float _S4480 = 2.0f * v_74;
    float2  _S4481 = _S4469 * make_float2 (radial_7) + make_float2 (_S4476 * v_74 + (*dist_coeffs_41)[int(5)] * (r2_74 + _S4477 * u_74) + (*dist_coeffs_41)[int(6)] * r2_74, _S4479 * v_74 + (*dist_coeffs_41)[int(4)] * (r2_74 + _S4480 * v_74) + (*dist_coeffs_41)[int(7)] * r2_74);
    float2  _S4482 = _S4481 + make_float2 ((*dist_coeffs_41)[int(8)] * _S4481.x + (*dist_coeffs_41)[int(9)] * _S4481.y, 0.0f);
    float _S4483 = fx_34 * _S4482.x + cx_29;
    float _S4484 = fy_34 * _S4482.y + cy_29;
    float2  uv0_10 = make_float2 (_S4483, _S4484);
    float2  _S4485 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S4486 = vert1_c_10.z;
    float2  _S4487 = make_float2 (_S4486);
    float2  _S4488 = _S4485 / make_float2 (_S4486);
    float2  _S4489 = make_float2 (_S4486 * _S4486);
    float u_75 = _S4488.x;
    float v_75 = _S4488.y;
    float r2_75 = u_75 * u_75 + v_75 * v_75;
    float _S4490 = (*dist_coeffs_41)[int(2)] + r2_75 * (*dist_coeffs_41)[int(3)];
    float _S4491 = (*dist_coeffs_41)[int(1)] + r2_75 * _S4490;
    float _S4492 = (*dist_coeffs_41)[int(0)] + r2_75 * _S4491;
    float radial_8 = 1.0f + r2_75 * _S4492;
    float2  _S4493 = make_float2 (radial_8);
    float _S4494 = _S4475 * u_75;
    float _S4495 = 2.0f * u_75;
    float _S4496 = _S4478 * u_75;
    float _S4497 = 2.0f * v_75;
    float2  _S4498 = _S4488 * make_float2 (radial_8) + make_float2 (_S4494 * v_75 + (*dist_coeffs_41)[int(5)] * (r2_75 + _S4495 * u_75) + (*dist_coeffs_41)[int(6)] * r2_75, _S4496 * v_75 + (*dist_coeffs_41)[int(4)] * (r2_75 + _S4497 * v_75) + (*dist_coeffs_41)[int(7)] * r2_75);
    float2  _S4499 = _S4498 + make_float2 ((*dist_coeffs_41)[int(8)] * _S4498.x + (*dist_coeffs_41)[int(9)] * _S4498.y, 0.0f);
    float _S4500 = fx_34 * _S4499.x + cx_29;
    float _S4501 = fy_34 * _S4499.y + cy_29;
    float2  uv1_10 = make_float2 (_S4500, _S4501);
    float2  _S4502 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S4503 = vert2_c_10.z;
    float2  _S4504 = make_float2 (_S4503);
    float2  _S4505 = _S4502 / make_float2 (_S4503);
    float2  _S4506 = make_float2 (_S4503 * _S4503);
    float u_76 = _S4505.x;
    float v_76 = _S4505.y;
    float r2_76 = u_76 * u_76 + v_76 * v_76;
    float _S4507 = (*dist_coeffs_41)[int(2)] + r2_76 * (*dist_coeffs_41)[int(3)];
    float _S4508 = (*dist_coeffs_41)[int(1)] + r2_76 * _S4507;
    float _S4509 = (*dist_coeffs_41)[int(0)] + r2_76 * _S4508;
    float radial_9 = 1.0f + r2_76 * _S4509;
    float2  _S4510 = make_float2 (radial_9);
    float _S4511 = _S4475 * u_76;
    float _S4512 = 2.0f * u_76;
    float _S4513 = _S4478 * u_76;
    float _S4514 = 2.0f * v_76;
    float2  _S4515 = _S4505 * make_float2 (radial_9) + make_float2 (_S4511 * v_76 + (*dist_coeffs_41)[int(5)] * (r2_76 + _S4512 * u_76) + (*dist_coeffs_41)[int(6)] * r2_76, _S4513 * v_76 + (*dist_coeffs_41)[int(4)] * (r2_76 + _S4514 * v_76) + (*dist_coeffs_41)[int(7)] * r2_76);
    float2  _S4516 = _S4515 + make_float2 ((*dist_coeffs_41)[int(8)] * _S4515.x + (*dist_coeffs_41)[int(9)] * _S4515.y, 0.0f);
    float _S4517 = fx_34 * _S4516.x + cx_29;
    float _S4518 = fy_34 * _S4516.y + cy_29;
    float2  uv2_10 = make_float2 (_S4517, _S4518);
    float2  e0_14 = uv1_10 - uv0_10;
    float2  e1_14 = uv2_10 - uv1_10;
    float2  e2_6 = uv0_10 - uv2_10;
    float _S4519 = e0_14.x;
    float _S4520 = e1_14.y;
    float _S4521 = e0_14.y;
    float _S4522 = e1_14.x;
    float _S4523 = _S4519 * _S4520 - _S4521 * _S4522;
    float _S4524 = 1.0f - hardness_14.y;
    float _S4525 = -1.0f / _S4524;
    float _S4526 = _S4524 * _S4524;
    float _S4527 = s_primal_ctx_max_0(_S4483, _S4500);
    float _S4528 = s_primal_ctx_min_0(_S4483, _S4500);
    float _S4529 = s_primal_ctx_max_0(_S4484, _S4501);
    float _S4530 = s_primal_ctx_min_0(_S4484, _S4501);
    float3  _S4531 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S4532 = length_1(_S4531) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4533 = transpose_0(R_29);
    float3  _S4534 = mean_30 - - s_primal_ctx_mul_1(_S4533, t_28);
    float _S4535 = _S4534.x;
    float _S4536 = _S4534.y;
    float _S4537 = _S4534.z;
    float _S4538 = _S4535 * _S4535 + _S4536 * _S4536 + _S4537 * _S4537;
    float _S4539 = s_primal_ctx_sqrt_0(_S4538);
    float x_61 = _S4535 / _S4539;
    float3  _S4540 = make_float3 (x_61);
    float _S4541 = _S4539 * _S4539;
    float y_28 = _S4536 / _S4539;
    float z_25 = _S4537 / _S4539;
    float3  _S4542 = make_float3 (z_25);
    float _S4543 = - y_28;
    float3  _S4544 = make_float3 (_S4543);
    float z2_57 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_61 * x_61 - y_28 * y_28;
    float _S4545 = 2.0f * x_61;
    float fS1_25 = _S4545 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_57 - 0.31539157032966614f;
    float3  _S4546 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_61;
    float3  _S4547 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S4548 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S4549 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S4550 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_57 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S4551 = 1.86588168144226074f * z2_57 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S4551;
    float3  _S4552 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_61;
    float3  _S4553 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S4554 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S4555 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S4556 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_61 * fC1_25 - y_28 * fS1_25);
    float3  _S4557 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_61 * fS1_25 + y_28 * fC1_25);
    float3  _S4558 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4543) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_61) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S4559 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S4560 = make_float3 (0.0f);
    float3  _S4561 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S4562 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4563 = make_float3 (_S4562);
    float3  _S4564 = (*ch_coeffs_10)[int(1)] * make_float3 (_S4562);
    float3  _S4565 = _S4561 + _S4564 + make_float3 (0.5f);
    float3  _S4566 = _S4561 - _S4564 + make_float3 (0.5f);
    float3  _S4567 = vert1_c_10 - vert0_c_10;
    float3  _S4568 = vert2_c_10 - vert0_c_10;
    float3  _S4569 = s_primal_ctx_cross_0(_S4567, _S4568);
    float3  _S4570 = normalize_0(_S4569);
    float3  _S4571 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4570, mean_c_25)))))) * v_normal_2;
    float3  _S4572 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4573;
    (&_S4573)->primal_0 = _S4570;
    (&_S4573)->differential_0 = _S4572;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4574;
    (&_S4574)->primal_0 = mean_c_25;
    (&_S4574)->differential_0 = _S4572;
    s_bwd_prop_dot_0(&_S4573, &_S4574, 0.0f);
    float3  _S4575 = _S4571 + _S4573.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4576;
    (&_S4576)->primal_0 = _S4569;
    (&_S4576)->differential_0 = _S4572;
    s_bwd_normalize_impl_0(&_S4576, _S4575);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4577;
    (&_S4577)->primal_0 = _S4567;
    (&_S4577)->differential_0 = _S4572;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4578;
    (&_S4578)->primal_0 = _S4568;
    (&_S4578)->differential_0 = _S4572;
    s_bwd_prop_cross_0(&_S4577, &_S4578, _S4576.differential_0);
    float3  _S4579 = - _S4578.differential_0;
    float3  _S4580 = - _S4577.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4581;
    (&_S4581)->primal_0 = _S4566;
    (&_S4581)->differential_0 = _S4572;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4582;
    (&_S4582)->primal_0 = _S4560;
    (&_S4582)->differential_0 = _S4572;
    s_bwd_prop_max_0(&_S4581, &_S4582, (*v_rgbs_0)[int(2)]);
    float3  _S4583 = - _S4581.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4584;
    (&_S4584)->primal_0 = _S4565;
    (&_S4584)->differential_0 = _S4572;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4585;
    (&_S4585)->primal_0 = _S4560;
    (&_S4585)->differential_0 = _S4572;
    s_bwd_prop_max_0(&_S4584, &_S4585, (*v_rgbs_0)[int(1)]);
    float3  _S4586 = _S4563 * (_S4583 + _S4584.differential_0);
    float3  _S4587 = _S4581.differential_0 + _S4584.differential_0;
    float3  _S4588 = make_float3 (0.5f) * - _S4587;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4589;
    (&_S4589)->primal_0 = _S4559;
    (&_S4589)->differential_0 = _S4572;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4590;
    (&_S4590)->primal_0 = _S4560;
    (&_S4590)->differential_0 = _S4572;
    s_bwd_prop_max_0(&_S4589, &_S4590, (*v_rgbs_0)[int(0)]);
    float3  _S4591 = _S4588 + _S4589.differential_0;
    float3  _S4592 = _S4587 + _S4589.differential_0;
    float3  _S4593 = _S4557 * _S4592;
    float3  _S4594 = (*sh_coeffs_25)[int(15)] * _S4592;
    float3  _S4595 = _S4555 * _S4592;
    float3  _S4596 = (*sh_coeffs_25)[int(14)] * _S4592;
    float3  _S4597 = _S4553 * _S4592;
    float3  _S4598 = (*sh_coeffs_25)[int(13)] * _S4592;
    float3  _S4599 = _S4552 * _S4592;
    float3  _S4600 = (*sh_coeffs_25)[int(12)] * _S4592;
    float3  _S4601 = _S4554 * _S4592;
    float3  _S4602 = (*sh_coeffs_25)[int(11)] * _S4592;
    float3  _S4603 = _S4556 * _S4592;
    float3  _S4604 = (*sh_coeffs_25)[int(10)] * _S4592;
    float3  _S4605 = _S4558 * _S4592;
    float3  _S4606 = (*sh_coeffs_25)[int(9)] * _S4592;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S4606.x + _S4606.y + _S4606.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S4594.x + _S4594.y + _S4594.z);
    float _S4607 = _S4604.x + _S4604.y + _S4604.z;
    float _S4608 = _S4596.x + _S4596.y + _S4596.z;
    float _S4609 = _S4602.x + _S4602.y + _S4602.z;
    float _S4610 = _S4598.x + _S4598.y + _S4598.z;
    float _S4611 = _S4600.x + _S4600.y + _S4600.z;
    float _S4612 = - s_diff_fC2_T_7;
    float3  _S4613 = _S4549 * _S4592;
    float3  _S4614 = (*sh_coeffs_25)[int(8)] * _S4592;
    float3  _S4615 = _S4547 * _S4592;
    float3  _S4616 = (*sh_coeffs_25)[int(7)] * _S4592;
    float3  _S4617 = _S4546 * _S4592;
    float3  _S4618 = (*sh_coeffs_25)[int(6)] * _S4592;
    float3  _S4619 = _S4548 * _S4592;
    float3  _S4620 = (*sh_coeffs_25)[int(5)] * _S4592;
    float3  _S4621 = _S4550 * _S4592;
    float3  _S4622 = (*sh_coeffs_25)[int(4)] * _S4592;
    float _S4623 = _S4620.x + _S4620.y + _S4620.z;
    float _S4624 = _S4616.x + _S4616.y + _S4616.z;
    float _S4625 = fTmp1B_25 * _S4607 + x_61 * s_diff_fS2_T_7 + y_28 * _S4612 + 0.54627424478530884f * (_S4622.x + _S4622.y + _S4622.z);
    float _S4626 = fTmp1B_25 * _S4608 + y_28 * s_diff_fS2_T_7 + x_61 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S4614.x + _S4614.y + _S4614.z);
    float _S4627 = y_28 * - _S4626;
    float _S4628 = x_61 * _S4626;
    float _S4629 = z_25 * (1.86588168144226074f * (z_25 * _S4611) + -2.28522896766662598f * (y_28 * _S4609 + x_61 * _S4610) + 0.94617468118667603f * (_S4618.x + _S4618.y + _S4618.z));
    float3  _S4630 = make_float3 (0.48860251903533936f) * _S4592;
    float3  _S4631 = - _S4630;
    float3  _S4632 = _S4540 * _S4631;
    float3  _S4633 = (*sh_coeffs_25)[int(3)] * _S4631;
    float3  _S4634 = _S4542 * _S4630;
    float3  _S4635 = (*sh_coeffs_25)[int(2)] * _S4630;
    float3  _S4636 = _S4544 * _S4630;
    float3  _S4637 = (*sh_coeffs_25)[int(1)] * _S4630;
    float _S4638 = (_S4551 * _S4611 + 1.44530570507049561f * (fS1_25 * _S4607 + fC1_25 * _S4608) + -1.09254848957061768f * (y_28 * _S4623 + x_61 * _S4624) + _S4629 + _S4629 + _S4635.x + _S4635.y + _S4635.z) / _S4541;
    float _S4639 = _S4539 * _S4638;
    float _S4640 = (fTmp0C_25 * _S4609 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S4612 + fTmp0B_25 * _S4623 + _S4545 * _S4625 + _S4627 + _S4627 + - (_S4637.x + _S4637.y + _S4637.z)) / _S4541;
    float _S4641 = _S4539 * _S4640;
    float _S4642 = (fTmp0C_25 * _S4610 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4624 + 2.0f * (y_28 * _S4625) + _S4628 + _S4628 + _S4633.x + _S4633.y + _S4633.z) / _S4541;
    float _S4643 = _S4539 * _S4642;
    float _S4644 = _S4537 * - _S4638 + _S4536 * - _S4640 + _S4535 * - _S4642;
    DiffPair_float_0 _S4645;
    (&_S4645)->primal_0 = _S4538;
    (&_S4645)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4645, _S4644);
    float _S4646 = _S4537 * _S4645.differential_0;
    float _S4647 = _S4536 * _S4645.differential_0;
    float _S4648 = _S4535 * _S4645.differential_0;
    float3  _S4649 = make_float3 (0.282094806432724f) * _S4592;
    float3  _S4650 = make_float3 (_S4643 + _S4648 + _S4648, _S4641 + _S4647 + _S4647, _S4639 + _S4646 + _S4646);
    float3  _S4651 = - - _S4650;
    Matrix<float, 3, 3>  _S4652 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4653;
    (&_S4653)->primal_0 = _S4533;
    (&_S4653)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4654;
    (&_S4654)->primal_0 = t_28;
    (&_S4654)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4653, &_S4654, _S4651);
    Matrix<float, 3, 3>  _S4655 = transpose_0(_S4653.differential_0);
    DiffPair_float_0 _S4656;
    (&_S4656)->primal_0 = _S4532;
    (&_S4656)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4656, v_depth_9);
    float _S4657 = 0.3333333432674408f * _S4656.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4658;
    (&_S4658)->primal_0 = _S4531;
    (&_S4658)->differential_0 = _S4572;
    s_bwd_length_impl_0(&_S4658, _S4657);
    DiffPair_float_0 _S4659;
    (&_S4659)->primal_0 = _S4530;
    (&_S4659)->differential_0 = 0.0f;
    DiffPair_float_0 _S4660;
    (&_S4660)->primal_0 = _S4518;
    (&_S4660)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4659, &_S4660, 0.0f);
    DiffPair_float_0 _S4661;
    (&_S4661)->primal_0 = _S4484;
    (&_S4661)->differential_0 = 0.0f;
    DiffPair_float_0 _S4662;
    (&_S4662)->primal_0 = _S4501;
    (&_S4662)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4661, &_S4662, _S4659.differential_0);
    DiffPair_float_0 _S4663;
    (&_S4663)->primal_0 = _S4529;
    (&_S4663)->differential_0 = 0.0f;
    DiffPair_float_0 _S4664;
    (&_S4664)->primal_0 = _S4518;
    (&_S4664)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4663, &_S4664, 0.0f);
    DiffPair_float_0 _S4665;
    (&_S4665)->primal_0 = _S4484;
    (&_S4665)->differential_0 = 0.0f;
    DiffPair_float_0 _S4666;
    (&_S4666)->primal_0 = _S4501;
    (&_S4666)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4665, &_S4666, _S4663.differential_0);
    DiffPair_float_0 _S4667;
    (&_S4667)->primal_0 = _S4528;
    (&_S4667)->differential_0 = 0.0f;
    DiffPair_float_0 _S4668;
    (&_S4668)->primal_0 = _S4517;
    (&_S4668)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4667, &_S4668, 0.0f);
    DiffPair_float_0 _S4669;
    (&_S4669)->primal_0 = _S4483;
    (&_S4669)->differential_0 = 0.0f;
    DiffPair_float_0 _S4670;
    (&_S4670)->primal_0 = _S4500;
    (&_S4670)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4669, &_S4670, _S4667.differential_0);
    DiffPair_float_0 _S4671;
    (&_S4671)->primal_0 = _S4527;
    (&_S4671)->differential_0 = 0.0f;
    DiffPair_float_0 _S4672;
    (&_S4672)->primal_0 = _S4517;
    (&_S4672)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4671, &_S4672, 0.0f);
    DiffPair_float_0 _S4673;
    (&_S4673)->primal_0 = _S4483;
    (&_S4673)->differential_0 = 0.0f;
    DiffPair_float_0 _S4674;
    (&_S4674)->primal_0 = _S4500;
    (&_S4674)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4673, &_S4674, _S4671.differential_0);
    DiffPair_float_0 _S4675;
    (&_S4675)->primal_0 = _S4525;
    (&_S4675)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4675, 0.0f);
    float _S4676 = - (-1.0f * - (_S4675.differential_0 / _S4526));
    float2  _S4677 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4678;
    (&_S4678)->primal_0 = e2_6;
    (&_S4678)->differential_0 = _S4677;
    s_bwd_length_impl_1(&_S4678, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4679;
    (&_S4679)->primal_0 = e1_14;
    (&_S4679)->differential_0 = _S4677;
    s_bwd_length_impl_1(&_S4679, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4680;
    (&_S4680)->primal_0 = e0_14;
    (&_S4680)->differential_0 = _S4677;
    s_bwd_length_impl_1(&_S4680, -0.0f);
    DiffPair_float_0 _S4681;
    (&_S4681)->primal_0 = _S4523;
    (&_S4681)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4681, 0.0f);
    float _S4682 = - _S4681.differential_0;
    float2  _S4683 = _S4679.differential_0 + make_float2 (_S4521 * _S4682, _S4519 * _S4681.differential_0);
    float2  _S4684 = _S4680.differential_0 + make_float2 (_S4520 * _S4681.differential_0, _S4522 * _S4682);
    float2  _S4685 = - _S4678.differential_0 + _S4683;
    float _S4686 = fx_34 * (_S4668.differential_0 + _S4672.differential_0 + _S4685.x);
    float2  _S4687 = make_float2 (_S4686, fy_34 * (_S4660.differential_0 + _S4664.differential_0 + _S4685.y)) + make_float2 ((*dist_coeffs_41)[int(8)] * _S4686, (*dist_coeffs_41)[int(9)] * _S4686);
    float2  _S4688 = _S4505 * _S4687;
    float _S4689 = (*dist_coeffs_41)[int(4)] * _S4687.y;
    float _S4690 = (*dist_coeffs_41)[int(5)] * _S4687.x;
    float _S4691 = _S4688.x + _S4688.y;
    float _S4692 = r2_76 * _S4691;
    float _S4693 = r2_76 * _S4692;
    float _S4694 = (*dist_coeffs_41)[int(7)] * _S4687.y + _S4689 + (*dist_coeffs_41)[int(6)] * _S4687.x + _S4690 + _S4509 * _S4691 + _S4508 * _S4692 + _S4507 * _S4693 + (*dist_coeffs_41)[int(3)] * (r2_76 * _S4693);
    float _S4695 = v_76 * _S4694;
    float _S4696 = u_76 * _S4694;
    float2  _S4697 = (_S4510 * _S4687 + make_float2 (_S4478 * (v_76 * _S4687.y) + _S4512 * _S4690 + 2.0f * (u_76 * _S4690) + _S4475 * (v_76 * _S4687.x) + _S4696 + _S4696, _S4514 * _S4689 + 2.0f * (v_76 * _S4689) + _S4513 * _S4687.y + _S4511 * _S4687.x + _S4695 + _S4695)) / _S4506;
    float2  _S4698 = _S4502 * - _S4697;
    float2  _S4699 = _S4504 * _S4697;
    float2  _S4700 = - _S4683 + _S4684;
    float _S4701 = fx_34 * (_S4670.differential_0 + _S4674.differential_0 + _S4700.x);
    float2  _S4702 = make_float2 (_S4701, fy_34 * (_S4662.differential_0 + _S4666.differential_0 + _S4700.y)) + make_float2 ((*dist_coeffs_41)[int(8)] * _S4701, (*dist_coeffs_41)[int(9)] * _S4701);
    float2  _S4703 = _S4488 * _S4702;
    float _S4704 = (*dist_coeffs_41)[int(4)] * _S4702.y;
    float _S4705 = (*dist_coeffs_41)[int(5)] * _S4702.x;
    float _S4706 = _S4703.x + _S4703.y;
    float _S4707 = r2_75 * _S4706;
    float _S4708 = r2_75 * _S4707;
    float _S4709 = (*dist_coeffs_41)[int(7)] * _S4702.y + _S4704 + (*dist_coeffs_41)[int(6)] * _S4702.x + _S4705 + _S4492 * _S4706 + _S4491 * _S4707 + _S4490 * _S4708 + (*dist_coeffs_41)[int(3)] * (r2_75 * _S4708);
    float _S4710 = v_75 * _S4709;
    float _S4711 = u_75 * _S4709;
    float2  _S4712 = (_S4493 * _S4702 + make_float2 (_S4478 * (v_75 * _S4702.y) + _S4495 * _S4705 + 2.0f * (u_75 * _S4705) + _S4475 * (v_75 * _S4702.x) + _S4711 + _S4711, _S4497 * _S4704 + 2.0f * (v_75 * _S4704) + _S4496 * _S4702.y + _S4494 * _S4702.x + _S4710 + _S4710)) / _S4489;
    float2  _S4713 = _S4485 * - _S4712;
    float2  _S4714 = _S4487 * _S4712;
    float _S4715 = _S4713.x + _S4713.y;
    float2  _S4716 = _S4678.differential_0 + - _S4684;
    float _S4717 = fx_34 * (_S4669.differential_0 + _S4673.differential_0 + _S4716.x);
    float2  _S4718 = make_float2 (_S4717, fy_34 * (_S4661.differential_0 + _S4665.differential_0 + _S4716.y)) + make_float2 ((*dist_coeffs_41)[int(8)] * _S4717, (*dist_coeffs_41)[int(9)] * _S4717);
    float2  _S4719 = _S4469 * _S4718;
    float _S4720 = (*dist_coeffs_41)[int(4)] * _S4718.y;
    float _S4721 = (*dist_coeffs_41)[int(5)] * _S4718.x;
    float _S4722 = _S4719.x + _S4719.y;
    float _S4723 = r2_74 * _S4722;
    float _S4724 = r2_74 * _S4723;
    float _S4725 = (*dist_coeffs_41)[int(7)] * _S4718.y + _S4720 + (*dist_coeffs_41)[int(6)] * _S4718.x + _S4721 + _S4473 * _S4722 + _S4472 * _S4723 + _S4471 * _S4724 + (*dist_coeffs_41)[int(3)] * (r2_74 * _S4724);
    float _S4726 = v_74 * _S4725;
    float _S4727 = u_74 * _S4725;
    float2  _S4728 = (_S4474 * _S4718 + make_float2 (_S4478 * (v_74 * _S4718.y) + _S4477 * _S4721 + 2.0f * (u_74 * _S4721) + _S4475 * (v_74 * _S4718.x) + _S4727 + _S4727, _S4480 * _S4720 + 2.0f * (v_74 * _S4720) + _S4479 * _S4718.y + _S4476 * _S4718.x + _S4726 + _S4726)) / _S4470;
    float2  _S4729 = _S4466 * - _S4728;
    float2  _S4730 = _S4468 * _S4728;
    float _S4731 = _S4729.x + _S4729.y;
    float3  _S4732 = _S4578.differential_0 + _S4658.differential_0 + make_float3 (_S4699.x, _S4699.y, _S4698.x + _S4698.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4733;
    (&_S4733)->primal_0 = R_29;
    (&_S4733)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4734;
    (&_S4734)->primal_0 = vert2_7;
    (&_S4734)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4733, &_S4734, _S4732);
    float3  _S4735 = _S4577.differential_0 + _S4658.differential_0 + make_float3 (_S4714.x, _S4714.y, _S4715);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4736;
    (&_S4736)->primal_0 = R_29;
    (&_S4736)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4737;
    (&_S4737)->primal_0 = vert1_7;
    (&_S4737)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4736, &_S4737, _S4735);
    float3  _S4738 = _S4579 + _S4580 + _S4658.differential_0 + make_float3 (_S4730.x, _S4730.y, _S4731);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4739;
    (&_S4739)->primal_0 = R_29;
    (&_S4739)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4740;
    (&_S4740)->primal_0 = vert0_7;
    (&_S4740)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4739, &_S4740, _S4738);
    float3  _S4741 = (*v_verts_0)[int(2)] + _S4734.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4742;
    (&_S4742)->primal_0 = _S4460;
    (&_S4742)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4743;
    (&_S4743)->primal_0 = _S4465;
    (&_S4743)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4742, &_S4743, _S4741);
    float _S4744 = - _S4743.differential_0.y;
    float _S4745 = _S4464 * _S4743.differential_0.x;
    float _S4746 = - (_S4456 * _S4743.differential_0.x);
    float3  _S4747 = (*v_verts_0)[int(1)] + _S4737.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4748;
    (&_S4748)->primal_0 = _S4460;
    (&_S4748)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4749;
    (&_S4749)->primal_0 = _S4463;
    (&_S4749)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4748, &_S4749, _S4747);
    float _S4750 = _S4456 * _S4749.differential_0.x;
    float _S4751 = _S4462 * _S4749.differential_0.x;
    float3  _S4752 = (*v_verts_0)[int(0)] + _S4740.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4753;
    (&_S4753)->primal_0 = _S4460;
    (&_S4753)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4754;
    (&_S4754)->primal_0 = _S4461;
    (&_S4754)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4753, &_S4754, _S4752);
    Matrix<float, 3, 3>  _S4755 = transpose_0(_S4742.differential_0 + _S4748.differential_0 + _S4753.differential_0);
    float _S4756 = 2.0f * - _S4755.rows[int(2)].z;
    float _S4757 = 2.0f * _S4755.rows[int(2)].y;
    float _S4758 = 2.0f * _S4755.rows[int(2)].x;
    float _S4759 = 2.0f * _S4755.rows[int(1)].z;
    float _S4760 = 2.0f * - _S4755.rows[int(1)].y;
    float _S4761 = 2.0f * _S4755.rows[int(1)].x;
    float _S4762 = 2.0f * _S4755.rows[int(0)].z;
    float _S4763 = 2.0f * _S4755.rows[int(0)].y;
    float _S4764 = 2.0f * - _S4755.rows[int(0)].x;
    float _S4765 = - _S4761 + _S4763;
    float _S4766 = _S4758 + - _S4762;
    float _S4767 = - _S4757 + _S4759;
    float _S4768 = _S4757 + _S4759;
    float _S4769 = _S4758 + _S4762;
    float _S4770 = _S4761 + _S4763;
    float _S4771 = quat_31.w * (_S4760 + _S4764);
    float _S4772 = quat_31.z * (_S4756 + _S4764);
    float _S4773 = quat_31.y * (_S4756 + _S4760);
    float _S4774 = quat_31.x * _S4765 + quat_31.z * _S4768 + quat_31.y * _S4769 + _S4771 + _S4771;
    float _S4775 = quat_31.x * _S4766 + quat_31.w * _S4768 + quat_31.y * _S4770 + _S4772 + _S4772;
    float _S4776 = quat_31.x * _S4767 + quat_31.w * _S4769 + quat_31.z * _S4770 + _S4773 + _S4773;
    float _S4777 = quat_31.w * _S4765 + quat_31.z * _S4766 + quat_31.y * _S4767;
    float _S4778 = _S4746 + _S4750;
    float _S4779 = 0.5f * - _S4778;
    float _S4780 = _S4744 + _S4749.differential_0.y;
    DiffPair_float_0 _S4781;
    (&_S4781)->primal_0 = _S4457;
    (&_S4781)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4781, _S4780);
    float _S4782 = _S4779 + _S4781.differential_0;
    float _S4783 = _S4745 + _S4751 + _S4754.differential_0.x;
    DiffPair_float_0 _S4784;
    (&_S4784)->primal_0 = _S4455;
    (&_S4784)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4784, _S4783);
    float _S4785 = _S4779 + _S4784.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4786;
    (&_S4786)->primal_0 = R_29;
    (&_S4786)->differential_0 = _S4652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4787;
    (&_S4787)->primal_0 = mean_30;
    (&_S4787)->differential_0 = _S4572;
    s_bwd_prop_mul_1(&_S4786, &_S4787, _S4574.differential_0);
    float3  _S4788 = _S4654.differential_0 + _S4732 + _S4735 + _S4738 + _S4574.differential_0;
    Matrix<float, 3, 3>  _S4789 = _S4655 + _S4733.differential_0 + _S4736.differential_0 + _S4739.differential_0 + _S4786.differential_0;
    FixedArray<float3 , 2>  _S4790;
    _S4790[int(0)] = _S4572;
    _S4790[int(1)] = _S4572;
    _S4790[int(1)] = _S4586;
    _S4790[int(0)] = _S4591;
    FixedArray<float3 , 16>  _S4791;
    _S4791[int(0)] = _S4572;
    _S4791[int(1)] = _S4572;
    _S4791[int(2)] = _S4572;
    _S4791[int(3)] = _S4572;
    _S4791[int(4)] = _S4572;
    _S4791[int(5)] = _S4572;
    _S4791[int(6)] = _S4572;
    _S4791[int(7)] = _S4572;
    _S4791[int(8)] = _S4572;
    _S4791[int(9)] = _S4572;
    _S4791[int(10)] = _S4572;
    _S4791[int(11)] = _S4572;
    _S4791[int(12)] = _S4572;
    _S4791[int(13)] = _S4572;
    _S4791[int(14)] = _S4572;
    _S4791[int(15)] = _S4572;
    _S4791[int(15)] = _S4593;
    _S4791[int(14)] = _S4595;
    _S4791[int(13)] = _S4597;
    _S4791[int(12)] = _S4599;
    _S4791[int(11)] = _S4601;
    _S4791[int(10)] = _S4603;
    _S4791[int(9)] = _S4605;
    _S4791[int(8)] = _S4613;
    _S4791[int(7)] = _S4615;
    _S4791[int(6)] = _S4617;
    _S4791[int(5)] = _S4619;
    _S4791[int(4)] = _S4621;
    _S4791[int(3)] = _S4632;
    _S4791[int(2)] = _S4634;
    _S4791[int(1)] = _S4636;
    _S4791[int(0)] = _S4649;
    float2  _S4792 = make_float2 (0.0f, _S4676);
    float3  _S4793 = make_float3 (_S4785, _S4782, _S4778);
    float4  _S4794 = make_float4 (0.0f);
    *&((&_S4794)->w) = _S4774;
    *&((&_S4794)->z) = _S4775;
    *&((&_S4794)->y) = _S4776;
    *&((&_S4794)->x) = _S4777;
    *v_mean_9 = _S4650 + _S4741 + _S4747 + _S4752 + _S4787.differential_0;
    *v_quat_8 = _S4794;
    *v_scale_8 = _S4793;
    *v_hardness_4 = _S4792;
    *v_sh_coeffs_7 = _S4791;
    *v_ch_coeffs_2 = _S4790;
    *v_R_8 = _S4789;
    *v_t_8 = _S4788;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_31, float4  quat_32, float3  scale_31, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_30, float3  t_29, float fx_35, float fy_35, float cx_30, float cy_30, FixedArray<float, 10>  * dist_coeffs_42, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_1(R_30, mean_31) + t_29;
    float _S4795 = scale_31.x;
    float _S4796 = s_primal_ctx_exp_1(_S4795);
    float _S4797 = scale_31.y;
    float _S4798 = s_primal_ctx_exp_1(_S4797);
    float sz_12 = scale_31.z - 0.5f * (_S4795 + _S4797);
    float _S4799 = quat_32.y;
    float x2_32 = _S4799 * _S4799;
    float y2_32 = quat_32.z * quat_32.z;
    float z2_58 = quat_32.w * quat_32.w;
    float xy_32 = quat_32.y * quat_32.z;
    float xz_32 = quat_32.y * quat_32.w;
    float yz_32 = quat_32.z * quat_32.w;
    float wx_32 = quat_32.x * quat_32.y;
    float wy_32 = quat_32.x * quat_32.z;
    float wz_32 = quat_32.x * quat_32.w;
    Matrix<float, 3, 3>  _S4800 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_32 + z2_58), 2.0f * (xy_32 + wz_32), 2.0f * (xz_32 - wy_32), 2.0f * (xy_32 - wz_32), 1.0f - 2.0f * (x2_32 + z2_58), 2.0f * (yz_32 + wx_32), 2.0f * (xz_32 + wy_32), 2.0f * (yz_32 - wx_32), 1.0f - 2.0f * (x2_32 + y2_32)));
    float3  _S4801 = make_float3 (_S4796, 0.0f, 0.0f);
    float3  vert0_8 = s_primal_ctx_mul_1(_S4800, _S4801) + mean_31;
    float _S4802 = -0.5f + sz_12;
    float3  _S4803 = make_float3 (_S4796 * _S4802, _S4798, 0.0f);
    float3  vert1_8 = s_primal_ctx_mul_1(_S4800, _S4803) + mean_31;
    float _S4804 = -0.5f - sz_12;
    float3  _S4805 = make_float3 (_S4796 * _S4804, - _S4798, 0.0f);
    float3  vert2_8 = s_primal_ctx_mul_1(_S4800, _S4805) + mean_31;
    float3  vert0_c_11 = s_primal_ctx_mul_1(R_30, vert0_8) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_1(R_30, vert1_8) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_1(R_30, vert2_8) + t_29;
    float2  _S4806 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4807 = length_0(_S4806);
    float _S4808 = vert0_c_11.z;
    float _S4809 = s_primal_ctx_atan2_0(_S4807, _S4808);
    bool _S4810 = _S4809 < 0.00100000004749745f;
    float k_12;
    float _S4811;
    float _S4812;
    float _S4813;
    if(_S4810)
    {
        float _S4814 = 1.0f - _S4809 * _S4809 / 3.0f;
        float _S4815 = _S4808 * _S4808;
        k_12 = _S4814 / _S4808;
        _S4811 = 0.0f;
        _S4812 = _S4815;
        _S4813 = _S4814;
    }
    else
    {
        float _S4816 = _S4807 * _S4807;
        k_12 = _S4809 / _S4807;
        _S4811 = _S4816;
        _S4812 = 0.0f;
        _S4813 = 0.0f;
    }
    float2  _S4817 = make_float2 (k_12);
    float2  _S4818 = _S4806 * make_float2 (k_12);
    float u_77 = _S4818.x;
    float v_77 = _S4818.y;
    float r2_77 = u_77 * u_77 + v_77 * v_77;
    float _S4819 = (*dist_coeffs_42)[int(2)] + r2_77 * (*dist_coeffs_42)[int(3)];
    float _S4820 = (*dist_coeffs_42)[int(1)] + r2_77 * _S4819;
    float _S4821 = (*dist_coeffs_42)[int(0)] + r2_77 * _S4820;
    float radial_10 = 1.0f + r2_77 * _S4821;
    float2  _S4822 = make_float2 (radial_10);
    float _S4823 = 2.0f * (*dist_coeffs_42)[int(4)];
    float _S4824 = _S4823 * u_77;
    float _S4825 = 2.0f * u_77;
    float _S4826 = 2.0f * (*dist_coeffs_42)[int(5)];
    float _S4827 = _S4826 * u_77;
    float _S4828 = 2.0f * v_77;
    float2  _S4829 = _S4818 * make_float2 (radial_10) + make_float2 (_S4824 * v_77 + (*dist_coeffs_42)[int(5)] * (r2_77 + _S4825 * u_77) + (*dist_coeffs_42)[int(6)] * r2_77, _S4827 * v_77 + (*dist_coeffs_42)[int(4)] * (r2_77 + _S4828 * v_77) + (*dist_coeffs_42)[int(7)] * r2_77);
    float2  _S4830 = _S4829 + make_float2 ((*dist_coeffs_42)[int(8)] * _S4829.x + (*dist_coeffs_42)[int(9)] * _S4829.y, 0.0f);
    float _S4831 = fx_35 * _S4830.x + cx_30;
    float _S4832 = fy_35 * _S4830.y + cy_30;
    float2  uv0_11 = make_float2 (_S4831, _S4832);
    float2  _S4833 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4834 = length_0(_S4833);
    float _S4835 = vert1_c_11.z;
    float _S4836 = s_primal_ctx_atan2_0(_S4834, _S4835);
    bool _S4837 = _S4836 < 0.00100000004749745f;
    float _S4838;
    float _S4839;
    float _S4840;
    if(_S4837)
    {
        float _S4841 = 1.0f - _S4836 * _S4836 / 3.0f;
        float _S4842 = _S4835 * _S4835;
        k_12 = _S4841 / _S4835;
        _S4838 = 0.0f;
        _S4839 = _S4842;
        _S4840 = _S4841;
    }
    else
    {
        float _S4843 = _S4834 * _S4834;
        k_12 = _S4836 / _S4834;
        _S4838 = _S4843;
        _S4839 = 0.0f;
        _S4840 = 0.0f;
    }
    float2  _S4844 = make_float2 (k_12);
    float2  _S4845 = _S4833 * make_float2 (k_12);
    float u_78 = _S4845.x;
    float v_78 = _S4845.y;
    float r2_78 = u_78 * u_78 + v_78 * v_78;
    float _S4846 = (*dist_coeffs_42)[int(2)] + r2_78 * (*dist_coeffs_42)[int(3)];
    float _S4847 = (*dist_coeffs_42)[int(1)] + r2_78 * _S4846;
    float _S4848 = (*dist_coeffs_42)[int(0)] + r2_78 * _S4847;
    float radial_11 = 1.0f + r2_78 * _S4848;
    float2  _S4849 = make_float2 (radial_11);
    float _S4850 = _S4823 * u_78;
    float _S4851 = 2.0f * u_78;
    float _S4852 = _S4826 * u_78;
    float _S4853 = 2.0f * v_78;
    float2  _S4854 = _S4845 * make_float2 (radial_11) + make_float2 (_S4850 * v_78 + (*dist_coeffs_42)[int(5)] * (r2_78 + _S4851 * u_78) + (*dist_coeffs_42)[int(6)] * r2_78, _S4852 * v_78 + (*dist_coeffs_42)[int(4)] * (r2_78 + _S4853 * v_78) + (*dist_coeffs_42)[int(7)] * r2_78);
    float2  _S4855 = _S4854 + make_float2 ((*dist_coeffs_42)[int(8)] * _S4854.x + (*dist_coeffs_42)[int(9)] * _S4854.y, 0.0f);
    float _S4856 = fx_35 * _S4855.x + cx_30;
    float _S4857 = fy_35 * _S4855.y + cy_30;
    float2  uv1_11 = make_float2 (_S4856, _S4857);
    float2  _S4858 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4859 = length_0(_S4858);
    float _S4860 = vert2_c_11.z;
    float _S4861 = s_primal_ctx_atan2_0(_S4859, _S4860);
    bool _S4862 = _S4861 < 0.00100000004749745f;
    float _S4863;
    float _S4864;
    float _S4865;
    if(_S4862)
    {
        float _S4866 = 1.0f - _S4861 * _S4861 / 3.0f;
        float _S4867 = _S4860 * _S4860;
        k_12 = _S4866 / _S4860;
        _S4863 = 0.0f;
        _S4864 = _S4867;
        _S4865 = _S4866;
    }
    else
    {
        float _S4868 = _S4859 * _S4859;
        k_12 = _S4861 / _S4859;
        _S4863 = _S4868;
        _S4864 = 0.0f;
        _S4865 = 0.0f;
    }
    float2  _S4869 = make_float2 (k_12);
    float2  _S4870 = _S4858 * make_float2 (k_12);
    float u_79 = _S4870.x;
    float v_79 = _S4870.y;
    float r2_79 = u_79 * u_79 + v_79 * v_79;
    float _S4871 = (*dist_coeffs_42)[int(2)] + r2_79 * (*dist_coeffs_42)[int(3)];
    float _S4872 = (*dist_coeffs_42)[int(1)] + r2_79 * _S4871;
    float _S4873 = (*dist_coeffs_42)[int(0)] + r2_79 * _S4872;
    float radial_12 = 1.0f + r2_79 * _S4873;
    float2  _S4874 = make_float2 (radial_12);
    float _S4875 = _S4823 * u_79;
    float _S4876 = 2.0f * u_79;
    float _S4877 = _S4826 * u_79;
    float _S4878 = 2.0f * v_79;
    float2  _S4879 = _S4870 * make_float2 (radial_12) + make_float2 (_S4875 * v_79 + (*dist_coeffs_42)[int(5)] * (r2_79 + _S4876 * u_79) + (*dist_coeffs_42)[int(6)] * r2_79, _S4877 * v_79 + (*dist_coeffs_42)[int(4)] * (r2_79 + _S4878 * v_79) + (*dist_coeffs_42)[int(7)] * r2_79);
    float2  _S4880 = _S4879 + make_float2 ((*dist_coeffs_42)[int(8)] * _S4879.x + (*dist_coeffs_42)[int(9)] * _S4879.y, 0.0f);
    float _S4881 = fx_35 * _S4880.x + cx_30;
    float _S4882 = fy_35 * _S4880.y + cy_30;
    float2  uv2_11 = make_float2 (_S4881, _S4882);
    float2  e0_15 = uv1_11 - uv0_11;
    float2  e1_15 = uv2_11 - uv1_11;
    float2  e2_7 = uv0_11 - uv2_11;
    float _S4883 = e0_15.x;
    float _S4884 = e1_15.y;
    float _S4885 = e0_15.y;
    float _S4886 = e1_15.x;
    float _S4887 = _S4883 * _S4884 - _S4885 * _S4886;
    float _S4888 = 1.0f - hardness_15.y;
    float _S4889 = -1.0f / _S4888;
    float _S4890 = _S4888 * _S4888;
    float _S4891 = s_primal_ctx_max_0(_S4831, _S4856);
    float _S4892 = s_primal_ctx_min_0(_S4831, _S4856);
    float _S4893 = s_primal_ctx_max_0(_S4832, _S4857);
    float _S4894 = s_primal_ctx_min_0(_S4832, _S4857);
    float3  _S4895 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S4896 = length_1(_S4895) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4897 = transpose_0(R_30);
    float3  _S4898 = mean_31 - - s_primal_ctx_mul_1(_S4897, t_29);
    float _S4899 = _S4898.x;
    float _S4900 = _S4898.y;
    float _S4901 = _S4898.z;
    float _S4902 = _S4899 * _S4899 + _S4900 * _S4900 + _S4901 * _S4901;
    float _S4903 = s_primal_ctx_sqrt_0(_S4902);
    float x_62 = _S4899 / _S4903;
    float3  _S4904 = make_float3 (x_62);
    float _S4905 = _S4903 * _S4903;
    float y_29 = _S4900 / _S4903;
    float z_26 = _S4901 / _S4903;
    float3  _S4906 = make_float3 (z_26);
    float _S4907 = - y_29;
    float3  _S4908 = make_float3 (_S4907);
    float z2_59 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_62 * x_62 - y_29 * y_29;
    float _S4909 = 2.0f * x_62;
    float fS1_26 = _S4909 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_59 - 0.31539157032966614f;
    float3  _S4910 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_62;
    float3  _S4911 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4912 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4913 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4914 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_59 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4915 = 1.86588168144226074f * z2_59 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4915;
    float3  _S4916 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_62;
    float3  _S4917 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4918 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4919 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4920 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_62 * fC1_26 - y_29 * fS1_26);
    float3  _S4921 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_62 * fS1_26 + y_29 * fC1_26);
    float3  _S4922 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4907) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_62) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4923 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4924 = make_float3 (0.0f);
    float3  _S4925 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4926 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4927 = make_float3 (_S4926);
    float3  _S4928 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4926);
    float3  _S4929 = _S4925 + _S4928 + make_float3 (0.5f);
    float3  _S4930 = _S4925 - _S4928 + make_float3 (0.5f);
    float3  _S4931 = vert1_c_11 - vert0_c_11;
    float3  _S4932 = vert2_c_11 - vert0_c_11;
    float3  _S4933 = s_primal_ctx_cross_0(_S4931, _S4932);
    float3  _S4934 = normalize_0(_S4933);
    float3  _S4935 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4934, mean_c_26)))))) * v_normal_3;
    float3  _S4936 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4937;
    (&_S4937)->primal_0 = _S4934;
    (&_S4937)->differential_0 = _S4936;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4938;
    (&_S4938)->primal_0 = mean_c_26;
    (&_S4938)->differential_0 = _S4936;
    s_bwd_prop_dot_0(&_S4937, &_S4938, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4939 = _S4938;
    float3  _S4940 = _S4935 + _S4937.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4941;
    (&_S4941)->primal_0 = _S4933;
    (&_S4941)->differential_0 = _S4936;
    s_bwd_normalize_impl_0(&_S4941, _S4940);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4942;
    (&_S4942)->primal_0 = _S4931;
    (&_S4942)->differential_0 = _S4936;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4943;
    (&_S4943)->primal_0 = _S4932;
    (&_S4943)->differential_0 = _S4936;
    s_bwd_prop_cross_0(&_S4942, &_S4943, _S4941.differential_0);
    float3  _S4944 = - _S4943.differential_0;
    float3  _S4945 = - _S4942.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4946;
    (&_S4946)->primal_0 = _S4930;
    (&_S4946)->differential_0 = _S4936;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4947;
    (&_S4947)->primal_0 = _S4924;
    (&_S4947)->differential_0 = _S4936;
    s_bwd_prop_max_0(&_S4946, &_S4947, (*v_rgbs_1)[int(2)]);
    float3  _S4948 = - _S4946.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4949;
    (&_S4949)->primal_0 = _S4929;
    (&_S4949)->differential_0 = _S4936;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4950;
    (&_S4950)->primal_0 = _S4924;
    (&_S4950)->differential_0 = _S4936;
    s_bwd_prop_max_0(&_S4949, &_S4950, (*v_rgbs_1)[int(1)]);
    float3  _S4951 = _S4927 * (_S4948 + _S4949.differential_0);
    float3  _S4952 = _S4946.differential_0 + _S4949.differential_0;
    float3  _S4953 = make_float3 (0.5f) * - _S4952;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4954;
    (&_S4954)->primal_0 = _S4923;
    (&_S4954)->differential_0 = _S4936;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4955;
    (&_S4955)->primal_0 = _S4924;
    (&_S4955)->differential_0 = _S4936;
    s_bwd_prop_max_0(&_S4954, &_S4955, (*v_rgbs_1)[int(0)]);
    float3  _S4956 = _S4953 + _S4954.differential_0;
    float3  _S4957 = _S4952 + _S4954.differential_0;
    float3  _S4958 = _S4921 * _S4957;
    float3  _S4959 = (*sh_coeffs_26)[int(15)] * _S4957;
    float3  _S4960 = _S4919 * _S4957;
    float3  _S4961 = (*sh_coeffs_26)[int(14)] * _S4957;
    float3  _S4962 = _S4917 * _S4957;
    float3  _S4963 = (*sh_coeffs_26)[int(13)] * _S4957;
    float3  _S4964 = _S4916 * _S4957;
    float3  _S4965 = (*sh_coeffs_26)[int(12)] * _S4957;
    float3  _S4966 = _S4918 * _S4957;
    float3  _S4967 = (*sh_coeffs_26)[int(11)] * _S4957;
    float3  _S4968 = _S4920 * _S4957;
    float3  _S4969 = (*sh_coeffs_26)[int(10)] * _S4957;
    float3  _S4970 = _S4922 * _S4957;
    float3  _S4971 = (*sh_coeffs_26)[int(9)] * _S4957;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4971.x + _S4971.y + _S4971.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4959.x + _S4959.y + _S4959.z);
    float _S4972 = _S4969.x + _S4969.y + _S4969.z;
    float _S4973 = _S4961.x + _S4961.y + _S4961.z;
    float _S4974 = _S4967.x + _S4967.y + _S4967.z;
    float _S4975 = _S4963.x + _S4963.y + _S4963.z;
    float _S4976 = _S4965.x + _S4965.y + _S4965.z;
    float _S4977 = - s_diff_fC2_T_8;
    float3  _S4978 = _S4913 * _S4957;
    float3  _S4979 = (*sh_coeffs_26)[int(8)] * _S4957;
    float3  _S4980 = _S4911 * _S4957;
    float3  _S4981 = (*sh_coeffs_26)[int(7)] * _S4957;
    float3  _S4982 = _S4910 * _S4957;
    float3  _S4983 = (*sh_coeffs_26)[int(6)] * _S4957;
    float3  _S4984 = _S4912 * _S4957;
    float3  _S4985 = (*sh_coeffs_26)[int(5)] * _S4957;
    float3  _S4986 = _S4914 * _S4957;
    float3  _S4987 = (*sh_coeffs_26)[int(4)] * _S4957;
    float _S4988 = _S4985.x + _S4985.y + _S4985.z;
    float _S4989 = _S4981.x + _S4981.y + _S4981.z;
    float _S4990 = fTmp1B_26 * _S4972 + x_62 * s_diff_fS2_T_8 + y_29 * _S4977 + 0.54627424478530884f * (_S4987.x + _S4987.y + _S4987.z);
    float _S4991 = fTmp1B_26 * _S4973 + y_29 * s_diff_fS2_T_8 + x_62 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S4979.x + _S4979.y + _S4979.z);
    float _S4992 = y_29 * - _S4991;
    float _S4993 = x_62 * _S4991;
    float _S4994 = z_26 * (1.86588168144226074f * (z_26 * _S4976) + -2.28522896766662598f * (y_29 * _S4974 + x_62 * _S4975) + 0.94617468118667603f * (_S4983.x + _S4983.y + _S4983.z));
    float3  _S4995 = make_float3 (0.48860251903533936f) * _S4957;
    float3  _S4996 = - _S4995;
    float3  _S4997 = _S4904 * _S4996;
    float3  _S4998 = (*sh_coeffs_26)[int(3)] * _S4996;
    float3  _S4999 = _S4906 * _S4995;
    float3  _S5000 = (*sh_coeffs_26)[int(2)] * _S4995;
    float3  _S5001 = _S4908 * _S4995;
    float3  _S5002 = (*sh_coeffs_26)[int(1)] * _S4995;
    float _S5003 = (_S4915 * _S4976 + 1.44530570507049561f * (fS1_26 * _S4972 + fC1_26 * _S4973) + -1.09254848957061768f * (y_29 * _S4988 + x_62 * _S4989) + _S4994 + _S4994 + _S5000.x + _S5000.y + _S5000.z) / _S4905;
    float _S5004 = _S4903 * _S5003;
    float _S5005 = (fTmp0C_26 * _S4974 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S4977 + fTmp0B_26 * _S4988 + _S4909 * _S4990 + _S4992 + _S4992 + - (_S5002.x + _S5002.y + _S5002.z)) / _S4905;
    float _S5006 = _S4903 * _S5005;
    float _S5007 = (fTmp0C_26 * _S4975 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S4989 + 2.0f * (y_29 * _S4990) + _S4993 + _S4993 + _S4998.x + _S4998.y + _S4998.z) / _S4905;
    float _S5008 = _S4903 * _S5007;
    float _S5009 = _S4901 * - _S5003 + _S4900 * - _S5005 + _S4899 * - _S5007;
    DiffPair_float_0 _S5010;
    (&_S5010)->primal_0 = _S4902;
    (&_S5010)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5010, _S5009);
    float _S5011 = _S4901 * _S5010.differential_0;
    float _S5012 = _S4900 * _S5010.differential_0;
    float _S5013 = _S4899 * _S5010.differential_0;
    float3  _S5014 = make_float3 (0.282094806432724f) * _S4957;
    float3  _S5015 = make_float3 (_S5008 + _S5013 + _S5013, _S5006 + _S5012 + _S5012, _S5004 + _S5011 + _S5011);
    float3  _S5016 = - - _S5015;
    Matrix<float, 3, 3>  _S5017 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5018;
    (&_S5018)->primal_0 = _S4897;
    (&_S5018)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5019;
    (&_S5019)->primal_0 = t_29;
    (&_S5019)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5018, &_S5019, _S5016);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5020 = _S5019;
    Matrix<float, 3, 3>  _S5021 = transpose_0(_S5018.differential_0);
    DiffPair_float_0 _S5022;
    (&_S5022)->primal_0 = _S4896;
    (&_S5022)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5022, v_depth_10);
    float _S5023 = 0.3333333432674408f * _S5022.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5024;
    (&_S5024)->primal_0 = _S4895;
    (&_S5024)->differential_0 = _S4936;
    s_bwd_length_impl_0(&_S5024, _S5023);
    DiffPair_float_0 _S5025;
    (&_S5025)->primal_0 = _S4894;
    (&_S5025)->differential_0 = 0.0f;
    DiffPair_float_0 _S5026;
    (&_S5026)->primal_0 = _S4882;
    (&_S5026)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5025, &_S5026, 0.0f);
    DiffPair_float_0 _S5027;
    (&_S5027)->primal_0 = _S4832;
    (&_S5027)->differential_0 = 0.0f;
    DiffPair_float_0 _S5028;
    (&_S5028)->primal_0 = _S4857;
    (&_S5028)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5027, &_S5028, _S5025.differential_0);
    DiffPair_float_0 _S5029;
    (&_S5029)->primal_0 = _S4893;
    (&_S5029)->differential_0 = 0.0f;
    DiffPair_float_0 _S5030;
    (&_S5030)->primal_0 = _S4882;
    (&_S5030)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5029, &_S5030, 0.0f);
    DiffPair_float_0 _S5031;
    (&_S5031)->primal_0 = _S4832;
    (&_S5031)->differential_0 = 0.0f;
    DiffPair_float_0 _S5032;
    (&_S5032)->primal_0 = _S4857;
    (&_S5032)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5031, &_S5032, _S5029.differential_0);
    DiffPair_float_0 _S5033;
    (&_S5033)->primal_0 = _S4892;
    (&_S5033)->differential_0 = 0.0f;
    DiffPair_float_0 _S5034;
    (&_S5034)->primal_0 = _S4881;
    (&_S5034)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5033, &_S5034, 0.0f);
    DiffPair_float_0 _S5035;
    (&_S5035)->primal_0 = _S4831;
    (&_S5035)->differential_0 = 0.0f;
    DiffPair_float_0 _S5036;
    (&_S5036)->primal_0 = _S4856;
    (&_S5036)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5035, &_S5036, _S5033.differential_0);
    DiffPair_float_0 _S5037;
    (&_S5037)->primal_0 = _S4891;
    (&_S5037)->differential_0 = 0.0f;
    DiffPair_float_0 _S5038;
    (&_S5038)->primal_0 = _S4881;
    (&_S5038)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5037, &_S5038, 0.0f);
    DiffPair_float_0 _S5039;
    (&_S5039)->primal_0 = _S4831;
    (&_S5039)->differential_0 = 0.0f;
    DiffPair_float_0 _S5040;
    (&_S5040)->primal_0 = _S4856;
    (&_S5040)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5039, &_S5040, _S5037.differential_0);
    DiffPair_float_0 _S5041;
    (&_S5041)->primal_0 = _S4889;
    (&_S5041)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S5041, 0.0f);
    float _S5042 = - (-1.0f * - (_S5041.differential_0 / _S4890));
    float2  _S5043 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5044;
    (&_S5044)->primal_0 = e2_7;
    (&_S5044)->differential_0 = _S5043;
    s_bwd_length_impl_1(&_S5044, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5045;
    (&_S5045)->primal_0 = e1_15;
    (&_S5045)->differential_0 = _S5043;
    s_bwd_length_impl_1(&_S5045, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5046;
    (&_S5046)->primal_0 = e0_15;
    (&_S5046)->differential_0 = _S5043;
    s_bwd_length_impl_1(&_S5046, -0.0f);
    DiffPair_float_0 _S5047;
    (&_S5047)->primal_0 = _S4887;
    (&_S5047)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S5047, 0.0f);
    float _S5048 = - _S5047.differential_0;
    float2  _S5049 = _S5045.differential_0 + make_float2 (_S4885 * _S5048, _S4883 * _S5047.differential_0);
    float2  _S5050 = - _S5049;
    float2  _S5051 = _S5046.differential_0 + make_float2 (_S4884 * _S5047.differential_0, _S4886 * _S5048);
    float2  _S5052 = - _S5051;
    float2  _S5053 = - _S5044.differential_0 + _S5049;
    float _S5054 = fx_35 * (_S5034.differential_0 + _S5038.differential_0 + _S5053.x);
    float2  _S5055 = make_float2 (_S5054, fy_35 * (_S5026.differential_0 + _S5030.differential_0 + _S5053.y)) + make_float2 ((*dist_coeffs_42)[int(8)] * _S5054, (*dist_coeffs_42)[int(9)] * _S5054);
    float2  _S5056 = _S4870 * _S5055;
    float2  _S5057 = _S4874 * _S5055;
    float _S5058 = (*dist_coeffs_42)[int(4)] * _S5055.y;
    float _S5059 = (*dist_coeffs_42)[int(5)] * _S5055.x;
    float _S5060 = _S5056.x + _S5056.y;
    float _S5061 = r2_79 * _S5060;
    float _S5062 = r2_79 * _S5061;
    float _S5063 = (*dist_coeffs_42)[int(7)] * _S5055.y + _S5058 + (*dist_coeffs_42)[int(6)] * _S5055.x + _S5059 + _S4873 * _S5060 + _S4872 * _S5061 + _S4871 * _S5062 + (*dist_coeffs_42)[int(3)] * (r2_79 * _S5062);
    float _S5064 = v_79 * _S5063;
    float _S5065 = u_79 * _S5063;
    float _S5066 = _S4878 * _S5058 + 2.0f * (v_79 * _S5058) + _S4877 * _S5055.y + _S4875 * _S5055.x + _S5064 + _S5064;
    float _S5067 = _S4826 * (v_79 * _S5055.y) + _S4876 * _S5059 + 2.0f * (u_79 * _S5059) + _S4823 * (v_79 * _S5055.x) + _S5065 + _S5065;
    float3  _S5068 = _S4942.differential_0 + _S5024.differential_0;
    float3  _S5069 = _S4944 + _S4945 + _S5024.differential_0;
    float3  _S5070 = _S4943.differential_0 + _S5024.differential_0;
    FixedArray<float3 , 2>  _S5071;
    _S5071[int(0)] = _S4936;
    _S5071[int(1)] = _S4936;
    _S5071[int(1)] = _S4951;
    _S5071[int(0)] = _S4956;
    float3  _S5072 = _S5071[int(0)];
    float3  _S5073 = _S5071[int(1)];
    FixedArray<float3 , 16>  _S5074;
    _S5074[int(0)] = _S4936;
    _S5074[int(1)] = _S4936;
    _S5074[int(2)] = _S4936;
    _S5074[int(3)] = _S4936;
    _S5074[int(4)] = _S4936;
    _S5074[int(5)] = _S4936;
    _S5074[int(6)] = _S4936;
    _S5074[int(7)] = _S4936;
    _S5074[int(8)] = _S4936;
    _S5074[int(9)] = _S4936;
    _S5074[int(10)] = _S4936;
    _S5074[int(11)] = _S4936;
    _S5074[int(12)] = _S4936;
    _S5074[int(13)] = _S4936;
    _S5074[int(14)] = _S4936;
    _S5074[int(15)] = _S4936;
    _S5074[int(7)] = _S4980;
    _S5074[int(0)] = _S5014;
    _S5074[int(1)] = _S5001;
    _S5074[int(2)] = _S4999;
    _S5074[int(3)] = _S4997;
    _S5074[int(4)] = _S4986;
    _S5074[int(5)] = _S4984;
    _S5074[int(6)] = _S4982;
    _S5074[int(15)] = _S4958;
    _S5074[int(8)] = _S4978;
    _S5074[int(9)] = _S4970;
    _S5074[int(10)] = _S4968;
    _S5074[int(11)] = _S4966;
    _S5074[int(12)] = _S4964;
    _S5074[int(13)] = _S4962;
    _S5074[int(14)] = _S4960;
    float3  _S5075 = _S5074[int(0)];
    float3  _S5076 = _S5074[int(1)];
    float3  _S5077 = _S5074[int(2)];
    float3  _S5078 = _S5074[int(3)];
    float3  _S5079 = _S5074[int(4)];
    float3  _S5080 = _S5074[int(5)];
    float3  _S5081 = _S5074[int(6)];
    float3  _S5082 = _S5074[int(7)];
    float3  _S5083 = _S5074[int(8)];
    float3  _S5084 = _S5074[int(9)];
    float3  _S5085 = _S5074[int(10)];
    float3  _S5086 = _S5074[int(11)];
    float3  _S5087 = _S5074[int(12)];
    float3  _S5088 = _S5074[int(13)];
    float3  _S5089 = _S5074[int(14)];
    float3  _S5090 = _S5074[int(15)];
    float _S5091 = _S5035.differential_0 + _S5039.differential_0;
    float2  _S5092 = _S5044.differential_0 + _S5052;
    float _S5093 = _S5028.differential_0 + _S5032.differential_0;
    float _S5094 = _S5027.differential_0 + _S5031.differential_0;
    float2  _S5095 = _S5050 + _S5051;
    float _S5096 = _S5036.differential_0 + _S5040.differential_0;
    float2  _S5097 = make_float2 (0.0f, _S5042);
    float2  _S5098 = _S5057 + make_float2 (_S5067, _S5066);
    float2  _S5099 = _S4858 * _S5098;
    float2  _S5100 = _S4869 * _S5098;
    float _S5101 = _S5099.x + _S5099.y;
    if(_S4862)
    {
        float _S5102 = _S5101 / _S4864;
        float _S5103 = _S4865 * - _S5102;
        float _S5104 = _S4861 * (0.3333333432674408f * - (_S4860 * _S5102));
        k_12 = _S5104 + _S5104;
        _S4863 = _S5103;
        _S4864 = 0.0f;
    }
    else
    {
        float _S5105 = _S5101 / _S4863;
        float _S5106 = _S4861 * - _S5105;
        k_12 = _S4859 * _S5105;
        _S4863 = 0.0f;
        _S4864 = _S5106;
    }
    DiffPair_float_0 _S5107;
    (&_S5107)->primal_0 = _S4859;
    (&_S5107)->differential_0 = 0.0f;
    DiffPair_float_0 _S5108;
    (&_S5108)->primal_0 = _S4860;
    (&_S5108)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5107, &_S5108, k_12);
    float _S5109 = _S5108.differential_0 + _S4863;
    float _S5110 = _S5107.differential_0 + _S4864;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5111;
    (&_S5111)->primal_0 = _S4858;
    (&_S5111)->differential_0 = _S5043;
    s_bwd_length_impl_1(&_S5111, _S5110);
    float2  _S5112 = _S5111.differential_0 + _S5100;
    float _S5113 = fx_35 * (_S5095.x + _S5096);
    float2  _S5114 = make_float2 (_S5113, fy_35 * (_S5095.y + _S5093)) + make_float2 ((*dist_coeffs_42)[int(8)] * _S5113, (*dist_coeffs_42)[int(9)] * _S5113);
    float2  _S5115 = _S4845 * _S5114;
    float _S5116 = (*dist_coeffs_42)[int(4)] * _S5114.y;
    float _S5117 = (*dist_coeffs_42)[int(5)] * _S5114.x;
    float _S5118 = _S5115.x + _S5115.y;
    float _S5119 = r2_78 * _S5118;
    float _S5120 = r2_78 * _S5119;
    float _S5121 = (*dist_coeffs_42)[int(7)] * _S5114.y + _S5116 + (*dist_coeffs_42)[int(6)] * _S5114.x + _S5117 + _S4848 * _S5118 + _S4847 * _S5119 + _S4846 * _S5120 + (*dist_coeffs_42)[int(3)] * (r2_78 * _S5120);
    float _S5122 = v_78 * _S5121;
    float _S5123 = u_78 * _S5121;
    float3  _S5124 = _S5070 + make_float3 (_S5112.x, _S5112.y, _S5109);
    float2  _S5125 = _S4849 * _S5114 + make_float2 (_S4826 * (v_78 * _S5114.y) + _S4851 * _S5117 + 2.0f * (u_78 * _S5117) + _S4823 * (v_78 * _S5114.x) + _S5123 + _S5123, _S4853 * _S5116 + 2.0f * (v_78 * _S5116) + _S4852 * _S5114.y + _S4850 * _S5114.x + _S5122 + _S5122);
    float2  _S5126 = _S4833 * _S5125;
    float2  _S5127 = _S4844 * _S5125;
    float _S5128 = _S5126.x + _S5126.y;
    if(_S4837)
    {
        float _S5129 = _S5128 / _S4839;
        float _S5130 = _S4840 * - _S5129;
        float _S5131 = _S4836 * (0.3333333432674408f * - (_S4835 * _S5129));
        k_12 = _S5131 + _S5131;
        _S4838 = _S5130;
        _S4839 = 0.0f;
    }
    else
    {
        float _S5132 = _S5128 / _S4838;
        float _S5133 = _S4836 * - _S5132;
        k_12 = _S4834 * _S5132;
        _S4838 = 0.0f;
        _S4839 = _S5133;
    }
    DiffPair_float_0 _S5134;
    (&_S5134)->primal_0 = _S4834;
    (&_S5134)->differential_0 = 0.0f;
    DiffPair_float_0 _S5135;
    (&_S5135)->primal_0 = _S4835;
    (&_S5135)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5134, &_S5135, k_12);
    float _S5136 = _S5135.differential_0 + _S4838;
    float _S5137 = _S5134.differential_0 + _S4839;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5138;
    (&_S5138)->primal_0 = _S4833;
    (&_S5138)->differential_0 = _S5043;
    s_bwd_length_impl_1(&_S5138, _S5137);
    float2  _S5139 = _S5138.differential_0 + _S5127;
    float _S5140 = fx_35 * (_S5092.x + _S5091);
    float2  _S5141 = make_float2 (_S5140, fy_35 * (_S5092.y + _S5094)) + make_float2 ((*dist_coeffs_42)[int(8)] * _S5140, (*dist_coeffs_42)[int(9)] * _S5140);
    float2  _S5142 = _S4818 * _S5141;
    float _S5143 = (*dist_coeffs_42)[int(4)] * _S5141.y;
    float _S5144 = (*dist_coeffs_42)[int(5)] * _S5141.x;
    float _S5145 = _S5142.x + _S5142.y;
    float _S5146 = r2_77 * _S5145;
    float _S5147 = r2_77 * _S5146;
    float _S5148 = (*dist_coeffs_42)[int(7)] * _S5141.y + _S5143 + (*dist_coeffs_42)[int(6)] * _S5141.x + _S5144 + _S4821 * _S5145 + _S4820 * _S5146 + _S4819 * _S5147 + (*dist_coeffs_42)[int(3)] * (r2_77 * _S5147);
    float _S5149 = v_77 * _S5148;
    float _S5150 = u_77 * _S5148;
    float3  _S5151 = _S5068 + make_float3 (_S5139.x, _S5139.y, _S5136);
    float2  _S5152 = _S4822 * _S5141 + make_float2 (_S4826 * (v_77 * _S5141.y) + _S4825 * _S5144 + 2.0f * (u_77 * _S5144) + _S4823 * (v_77 * _S5141.x) + _S5150 + _S5150, _S4828 * _S5143 + 2.0f * (v_77 * _S5143) + _S4827 * _S5141.y + _S4824 * _S5141.x + _S5149 + _S5149);
    float2  _S5153 = _S4806 * _S5152;
    float2  _S5154 = _S4817 * _S5152;
    float _S5155 = _S5153.x + _S5153.y;
    if(_S4810)
    {
        float _S5156 = _S5155 / _S4812;
        float _S5157 = _S4813 * - _S5156;
        float _S5158 = _S4809 * (0.3333333432674408f * - (_S4808 * _S5156));
        k_12 = _S5158 + _S5158;
        _S4811 = _S5157;
        _S4812 = 0.0f;
    }
    else
    {
        float _S5159 = _S5155 / _S4811;
        float _S5160 = _S4809 * - _S5159;
        k_12 = _S4807 * _S5159;
        _S4811 = 0.0f;
        _S4812 = _S5160;
    }
    DiffPair_float_0 _S5161;
    (&_S5161)->primal_0 = _S4807;
    (&_S5161)->differential_0 = 0.0f;
    DiffPair_float_0 _S5162;
    (&_S5162)->primal_0 = _S4808;
    (&_S5162)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5161, &_S5162, k_12);
    float _S5163 = _S5162.differential_0 + _S4811;
    float _S5164 = _S5161.differential_0 + _S4812;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5165;
    (&_S5165)->primal_0 = _S4806;
    (&_S5165)->differential_0 = _S5043;
    s_bwd_length_impl_1(&_S5165, _S5164);
    float2  _S5166 = _S5165.differential_0 + _S5154;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5167;
    (&_S5167)->primal_0 = vert2_c_11;
    (&_S5167)->differential_0 = _S4936;
    s_bwd_length_impl_0(&_S5167, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5168;
    (&_S5168)->primal_0 = vert1_c_11;
    (&_S5168)->differential_0 = _S4936;
    s_bwd_length_impl_0(&_S5168, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5169;
    (&_S5169)->primal_0 = vert0_c_11;
    (&_S5169)->differential_0 = _S4936;
    s_bwd_length_impl_0(&_S5169, 0.0f);
    float3  _S5170 = _S5167.differential_0 + _S5124;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5171;
    (&_S5171)->primal_0 = R_30;
    (&_S5171)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5172;
    (&_S5172)->primal_0 = vert2_8;
    (&_S5172)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5171, &_S5172, _S5170);
    float3  _S5173 = _S5168.differential_0 + _S5151;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5174;
    (&_S5174)->primal_0 = R_30;
    (&_S5174)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5175;
    (&_S5175)->primal_0 = vert1_8;
    (&_S5175)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5174, &_S5175, _S5173);
    float3  _S5176 = _S5169.differential_0 + _S5069 + make_float3 (_S5166.x, _S5166.y, _S5163);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5177;
    (&_S5177)->primal_0 = R_30;
    (&_S5177)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5178;
    (&_S5178)->primal_0 = vert0_8;
    (&_S5178)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5177, &_S5178, _S5176);
    float3  _S5179 = _S5172.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5180;
    (&_S5180)->primal_0 = _S4800;
    (&_S5180)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5181;
    (&_S5181)->primal_0 = _S4805;
    (&_S5181)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5180, &_S5181, _S5179);
    float _S5182 = - _S5181.differential_0.y;
    float _S5183 = _S4804 * _S5181.differential_0.x;
    float _S5184 = - (_S4796 * _S5181.differential_0.x);
    float3  _S5185 = _S5175.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5186;
    (&_S5186)->primal_0 = _S4800;
    (&_S5186)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5187;
    (&_S5187)->primal_0 = _S4803;
    (&_S5187)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5186, &_S5187, _S5185);
    float _S5188 = _S4796 * _S5187.differential_0.x;
    float _S5189 = _S4802 * _S5187.differential_0.x;
    float3  _S5190 = _S5178.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5191;
    (&_S5191)->primal_0 = _S4800;
    (&_S5191)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5192;
    (&_S5192)->primal_0 = _S4801;
    (&_S5192)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5191, &_S5192, _S5190);
    Matrix<float, 3, 3>  _S5193 = transpose_0(_S5180.differential_0 + _S5186.differential_0 + _S5191.differential_0);
    float _S5194 = 2.0f * - _S5193.rows[int(2)].z;
    float _S5195 = 2.0f * _S5193.rows[int(2)].y;
    float _S5196 = 2.0f * _S5193.rows[int(2)].x;
    float _S5197 = 2.0f * _S5193.rows[int(1)].z;
    float _S5198 = 2.0f * - _S5193.rows[int(1)].y;
    float _S5199 = 2.0f * _S5193.rows[int(1)].x;
    float _S5200 = 2.0f * _S5193.rows[int(0)].z;
    float _S5201 = 2.0f * _S5193.rows[int(0)].y;
    float _S5202 = 2.0f * - _S5193.rows[int(0)].x;
    float _S5203 = - _S5199 + _S5201;
    float _S5204 = _S5196 + - _S5200;
    float _S5205 = - _S5195 + _S5197;
    float _S5206 = _S5195 + _S5197;
    float _S5207 = _S5196 + _S5200;
    float _S5208 = _S5199 + _S5201;
    float _S5209 = quat_32.w * (_S5198 + _S5202);
    float _S5210 = quat_32.z * (_S5194 + _S5202);
    float _S5211 = quat_32.y * (_S5194 + _S5198);
    float _S5212 = quat_32.x * _S5203 + quat_32.z * _S5206 + quat_32.y * _S5207 + _S5209 + _S5209;
    float _S5213 = quat_32.x * _S5204 + quat_32.w * _S5206 + quat_32.y * _S5208 + _S5210 + _S5210;
    float _S5214 = quat_32.x * _S5205 + quat_32.w * _S5207 + quat_32.z * _S5208 + _S5211 + _S5211;
    float _S5215 = quat_32.w * _S5203 + quat_32.z * _S5204 + quat_32.y * _S5205;
    float _S5216 = _S5184 + _S5188;
    float _S5217 = 0.5f * - _S5216;
    float _S5218 = _S5182 + _S5187.differential_0.y;
    DiffPair_float_0 _S5219;
    (&_S5219)->primal_0 = _S4797;
    (&_S5219)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5219, _S5218);
    float _S5220 = _S5217 + _S5219.differential_0;
    float _S5221 = _S5183 + _S5189 + _S5192.differential_0.x;
    DiffPair_float_0 _S5222;
    (&_S5222)->primal_0 = _S4795;
    (&_S5222)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5222, _S5221);
    float _S5223 = _S5217 + _S5222.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5224;
    (&_S5224)->primal_0 = mean_c_26;
    (&_S5224)->differential_0 = _S4936;
    s_bwd_length_impl_0(&_S5224, 0.0f);
    float3  _S5225 = _S5224.differential_0 + _S4939.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5226;
    (&_S5226)->primal_0 = R_30;
    (&_S5226)->differential_0 = _S5017;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5227;
    (&_S5227)->primal_0 = mean_31;
    (&_S5227)->differential_0 = _S4936;
    s_bwd_prop_mul_1(&_S5226, &_S5227, _S5225);
    float3  _S5228 = _S5170 + _S5173 + _S5176 + _S5225 + _S5020.differential_0;
    Matrix<float, 3, 3>  _S5229 = _S5171.differential_0 + _S5174.differential_0 + _S5177.differential_0 + _S5226.differential_0 + _S5021;
    float3  _S5230 = make_float3 (_S5223, _S5220, _S5216);
    float4  _S5231 = make_float4 (0.0f);
    *&((&_S5231)->w) = _S5212;
    *&((&_S5231)->z) = _S5213;
    *&((&_S5231)->y) = _S5214;
    *&((&_S5231)->x) = _S5215;
    float4  _S5232 = _S5231;
    float3  _S5233 = _S5179 + _S5185 + _S5190 + _S5227.differential_0 + _S5015;
    *v_mean_10 = _S5233;
    *v_quat_9 = _S5232;
    *v_scale_9 = _S5230;
    *v_hardness_5 = _S5097;
    (*v_sh_coeffs_8)[int(0)] = _S5075;
    (*v_sh_coeffs_8)[int(1)] = _S5076;
    (*v_sh_coeffs_8)[int(2)] = _S5077;
    (*v_sh_coeffs_8)[int(3)] = _S5078;
    (*v_sh_coeffs_8)[int(4)] = _S5079;
    (*v_sh_coeffs_8)[int(5)] = _S5080;
    (*v_sh_coeffs_8)[int(6)] = _S5081;
    (*v_sh_coeffs_8)[int(7)] = _S5082;
    (*v_sh_coeffs_8)[int(8)] = _S5083;
    (*v_sh_coeffs_8)[int(9)] = _S5084;
    (*v_sh_coeffs_8)[int(10)] = _S5085;
    (*v_sh_coeffs_8)[int(11)] = _S5086;
    (*v_sh_coeffs_8)[int(12)] = _S5087;
    (*v_sh_coeffs_8)[int(13)] = _S5088;
    (*v_sh_coeffs_8)[int(14)] = _S5089;
    (*v_sh_coeffs_8)[int(15)] = _S5090;
    (*v_ch_coeffs_3)[int(0)] = _S5072;
    (*v_ch_coeffs_3)[int(1)] = _S5073;
    *v_R_9 = _S5229;
    *v_t_9 = _S5228;
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
    bool _S5234;
    if((*u_80) >= 0.0f)
    {
        _S5234 = (*v_80) >= 0.0f;
    }
    else
    {
        _S5234 = false;
    }
    if(_S5234)
    {
        _S5234 = (*u_80 + *v_80) <= 1.0f;
    }
    else
    {
        _S5234 = false;
    }
    if(_S5234)
    {
        _S5234 = (*t_30) >= 0.0f;
    }
    else
    {
        _S5234 = false;
    }
    return _S5234;
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
    bool _S5235;
    if(u_81 >= 0.0f)
    {
        _S5235 = v_81 >= 0.0f;
    }
    else
    {
        _S5235 = false;
    }
    if(_S5235)
    {
        _S5235 = (u_81 + v_81) <= 1.0f;
    }
    else
    {
        _S5235 = false;
    }
    if(_S5235)
    {
        _S5235 = t_31 >= 0.0f;
    }
    else
    {
        _S5235 = false;
    }
    if(!_S5235)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_81), (v_81)))), ((F32_sqrt((0.5f))) * (1.0f - u_81 - v_81)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_16.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_16.x;
    float _S5236;
    if(opac_0 < 0.0f)
    {
        _S5236 = 0.0f;
    }
    else
    {
        _S5236 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S5236;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_4, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_4, float _s_dOut_11)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5237 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5238 = *dpray_d_4;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_4).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S5239 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S5240 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_4).primal_0);
    float _S5241 = s_primal_ctx_dot_0((*dpray_d_4).primal_0, _S5239);
    float d_30 = 1.0f / _S5241;
    float _S5242 = _S5241 * _S5241;
    float3  _S5243 = - _S5240;
    float _S5244 = s_primal_ctx_dot_0(_S5243, v2v0_2);
    float u_82 = d_30 * _S5244;
    float _S5245 = s_primal_ctx_dot_0(_S5240, v1v0_2);
    float v_82 = d_30 * _S5245;
    float3  _S5246 = - _S5239;
    float t_32 = d_30 * s_primal_ctx_dot_0(_S5246, rov0_2);
    bool _S5247;
    if(u_82 >= 0.0f)
    {
        _S5247 = v_82 >= 0.0f;
    }
    else
    {
        _S5247 = false;
    }
    if(_S5247)
    {
        _S5247 = (u_82 + v_82) <= 1.0f;
    }
    else
    {
        _S5247 = false;
    }
    if(_S5247)
    {
        _S5247 = t_32 >= 0.0f;
    }
    else
    {
        _S5247 = false;
    }
    bool _S5248 = !!_S5247;
    float _S5249;
    float _S5250;
    float _S5251;
    float _S5252;
    float _S5253;
    float _S5254;
    float _S5255;
    float _S5256;
    float _S5257;
    float _S5258;
    float _S5259;
    if(_S5248)
    {
        float _S5260 = s_primal_ctx_min_0(u_82, v_82);
        float _S5261 = s_primal_ctx_sqrt_0(0.5f);
        float _S5262 = _S5261 * (1.0f - u_82 - v_82);
        float _S5263 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S5260, _S5262) * _S5263;
        float _S5264 = _S5237.primal_0.y;
        float _S5265 = 1.0f - opac_1;
        float _S5266 = 1.0f - s_primal_ctx_clamp_0(_S5264, 0.0f, 0.99989998340606689f);
        float _S5267 = 1.0f / _S5266;
        float _S5268 = _S5266 * _S5266;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S5265, _S5267);
        float o_1 = _S5237.primal_0.x;
        bool _S5269 = opac_1 < 0.0f;
        if(_S5269)
        {
            _S5249 = 0.0f;
        }
        else
        {
            _S5249 = o_1 * w_1;
        }
        _S5247 = _S5269;
        _S5250 = o_1;
        _S5251 = w_1;
        _S5252 = _S5265;
        _S5253 = _S5267;
        _S5254 = _S5268;
        _S5255 = _S5264;
        _S5256 = _S5263;
        _S5257 = _S5260;
        _S5258 = _S5262;
        _S5259 = _S5261;
    }
    else
    {
        _S5247 = false;
        _S5249 = 0.0f;
        _S5250 = 0.0f;
        _S5251 = 0.0f;
        _S5252 = 0.0f;
        _S5253 = 0.0f;
        _S5254 = 0.0f;
        _S5255 = 0.0f;
        _S5256 = 0.0f;
        _S5257 = 0.0f;
        _S5258 = 0.0f;
        _S5259 = 0.0f;
    }
    float2  _S5270 = make_float2 (0.0f);
    float2  _S5271;
    if(_S5248)
    {
        if(_S5247)
        {
            _S5249 = 0.0f;
            _S5250 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S5272;
            (&_S5272)->primal_0 = _S5249;
            (&_S5272)->differential_0 = 0.0f;
            DiffPair_float_0 _S5273;
            (&_S5273)->primal_0 = 0.99500000476837158f;
            (&_S5273)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S5272, &_S5273, _s_dOut_11);
            float _S5274 = _S5250 * _S5272.differential_0;
            _S5249 = _S5251 * _S5272.differential_0;
            _S5250 = _S5274;
        }
        float _S5275 = - _S5250;
        DiffPair_float_0 _S5276;
        (&_S5276)->primal_0 = _S5252;
        (&_S5276)->differential_0 = 0.0f;
        DiffPair_float_0 _S5277;
        (&_S5277)->primal_0 = _S5253;
        (&_S5277)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S5276, &_S5277, _S5275);
        float _S5278 = - - (_S5277.differential_0 / _S5254);
        float s_diff_opac_T_0 = - _S5276.differential_0;
        DiffPair_float_0 _S5279;
        (&_S5279)->primal_0 = _S5255;
        (&_S5279)->differential_0 = 0.0f;
        DiffPair_float_0 _S5280;
        (&_S5280)->primal_0 = 0.0f;
        (&_S5280)->differential_0 = 0.0f;
        DiffPair_float_0 _S5281;
        (&_S5281)->primal_0 = 0.99989998340606689f;
        (&_S5281)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S5279, &_S5280, &_S5281, _S5278);
        float _S5282 = _S5256 * s_diff_opac_T_0;
        DiffPair_float_0 _S5283;
        (&_S5283)->primal_0 = _S5257;
        (&_S5283)->differential_0 = 0.0f;
        DiffPair_float_0 _S5284;
        (&_S5284)->primal_0 = _S5258;
        (&_S5284)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5283, &_S5284, _S5282);
        float _S5285 = - (_S5259 * _S5284.differential_0);
        DiffPair_float_0 _S5286;
        (&_S5286)->primal_0 = u_82;
        (&_S5286)->differential_0 = 0.0f;
        DiffPair_float_0 _S5287;
        (&_S5287)->primal_0 = v_82;
        (&_S5287)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5286, &_S5287, _S5283.differential_0);
        float2  _S5288 = make_float2 (_S5249, _S5279.differential_0);
        float _S5289 = _S5285 + _S5287.differential_0;
        _S5249 = _S5285 + _S5286.differential_0;
        _S5250 = _S5289;
        _S5271 = _S5288;
    }
    else
    {
        _S5249 = 0.0f;
        _S5250 = 0.0f;
        _S5271 = _S5270;
    }
    float3  _S5290 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5291;
    (&_S5291)->primal_0 = _S5246;
    (&_S5291)->differential_0 = _S5290;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5292;
    (&_S5292)->primal_0 = rov0_2;
    (&_S5292)->differential_0 = _S5290;
    s_bwd_prop_dot_0(&_S5291, &_S5292, 0.0f);
    float3  _S5293 = - _S5291.differential_0;
    float _S5294 = d_30 * _S5250;
    float _S5295 = _S5245 * _S5250;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5296;
    (&_S5296)->primal_0 = _S5240;
    (&_S5296)->differential_0 = _S5290;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5297;
    (&_S5297)->primal_0 = v1v0_2;
    (&_S5297)->differential_0 = _S5290;
    s_bwd_prop_dot_0(&_S5296, &_S5297, _S5294);
    float _S5298 = d_30 * _S5249;
    float _S5299 = _S5244 * _S5249;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5300;
    (&_S5300)->primal_0 = _S5243;
    (&_S5300)->differential_0 = _S5290;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5301;
    (&_S5301)->primal_0 = v2v0_2;
    (&_S5301)->differential_0 = _S5290;
    s_bwd_prop_dot_0(&_S5300, &_S5301, _S5298);
    float3  _S5302 = - _S5300.differential_0;
    float _S5303 = - ((_S5295 + _S5299) / _S5242);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5304;
    (&_S5304)->primal_0 = _S5238.primal_0;
    (&_S5304)->differential_0 = _S5290;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5305;
    (&_S5305)->primal_0 = _S5239;
    (&_S5305)->differential_0 = _S5290;
    s_bwd_prop_dot_0(&_S5304, &_S5305, _S5303);
    float3  _S5306 = _S5296.differential_0 + _S5302;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5307;
    (&_S5307)->primal_0 = rov0_2;
    (&_S5307)->differential_0 = _S5290;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5308;
    (&_S5308)->primal_0 = _S5238.primal_0;
    (&_S5308)->differential_0 = _S5290;
    s_bwd_prop_cross_0(&_S5307, &_S5308, _S5306);
    float3  _S5309 = _S5293 + _S5305.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5310;
    (&_S5310)->primal_0 = v1v0_2;
    (&_S5310)->differential_0 = _S5290;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5311;
    (&_S5311)->primal_0 = v2v0_2;
    (&_S5311)->differential_0 = _S5290;
    s_bwd_prop_cross_0(&_S5310, &_S5311, _S5309);
    float3  _S5312 = _S5292.differential_0 + _S5307.differential_0;
    float3  _S5313 = _S5301.differential_0 + _S5311.differential_0;
    float3  _S5314 = _S5297.differential_0 + _S5310.differential_0;
    float3  _S5315 = - _S5312 + - _S5313 + - _S5314;
    float3  _S5316 = _S5304.differential_0 + _S5308.differential_0;
    dpray_d_4->primal_0 = (*dpray_d_4).primal_0;
    dpray_d_4->differential_0 = _S5316;
    dpray_o_4->primal_0 = (*dpray_o_4).primal_0;
    dpray_o_4->differential_0 = _S5312;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S5271;
    FixedArray<float3 , 3>  _S5317;
    _S5317[int(0)] = _S5290;
    _S5317[int(1)] = _S5290;
    _S5317[int(2)] = _S5290;
    _S5317[int(2)] = _S5313;
    _S5317[int(0)] = _S5315;
    _S5317[int(1)] = _S5314;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S5317;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5318, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S5319, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5320, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5321, float _S5322)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S5318, _S5319, _S5320, _S5321, _S5322);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_6, float2  hardness_17, float3  ray_o_7, float3  ray_d_7, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S5323 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5324 = { _S5323, _S5323, _S5323 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_6;
    (&dp_verts_0)->differential_0 = _S5324;
    float2  _S5325 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S5325;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_7;
    (&dp_ray_o_2)->differential_0 = _S5323;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_7;
    (&dp_ray_d_2)->differential_0 = _S5323;
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
    float3  _S5326 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S5327 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_5).primal_0);
    float _S5328 = s_primal_ctx_dot_0((*dpray_d_5).primal_0, _S5326);
    float d_32 = 1.0f / _S5328;
    float _S5329 = _S5328 * _S5328;
    float3  _S5330 = - _S5327;
    float _S5331 = s_primal_ctx_dot_0(_S5330, v2v0_4);
    float u_84 = d_32 * _S5331;
    float3  _S5332 = make_float3 (u_84);
    float _S5333 = s_primal_ctx_dot_0(_S5327, v1v0_4);
    float v_84 = d_32 * _S5333;
    float3  _S5334 = make_float3 (v_84);
    float3  _S5335 = - _S5326;
    float _S5336 = s_primal_ctx_dot_0(_S5335, rov0_4);
    float _S5337 = d_32 * _S5336;
    float3  _S5338 = make_float3 (1.0f - u_84 - v_84);
    DiffPair_float_0 _S5339;
    (&_S5339)->primal_0 = s_primal_ctx_max_0(_S5337, 9.999999960041972e-13f);
    (&_S5339)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5339, dpdepth_2);
    DiffPair_float_0 _S5340;
    (&_S5340)->primal_0 = _S5337;
    (&_S5340)->differential_0 = 0.0f;
    DiffPair_float_0 _S5341;
    (&_S5341)->primal_0 = 9.999999960041972e-13f;
    (&_S5341)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5340, &_S5341, _S5339.differential_0);
    float3  _S5342 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S5343 = _S5334 * dpcolor_1;
    float3  _S5344 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S5345 = _S5332 * dpcolor_1;
    float3  _S5346 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S5347 = _S5338 * dpcolor_1;
    float _S5348 = - (_S5346.x + _S5346.y + _S5346.z);
    float _S5349 = d_32 * _S5340.differential_0;
    float _S5350 = _S5336 * _S5340.differential_0;
    float3  _S5351 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5352;
    (&_S5352)->primal_0 = _S5335;
    (&_S5352)->differential_0 = _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5353;
    (&_S5353)->primal_0 = rov0_4;
    (&_S5353)->differential_0 = _S5351;
    s_bwd_prop_dot_0(&_S5352, &_S5353, _S5349);
    float3  _S5354 = - _S5352.differential_0;
    float _S5355 = _S5348 + _S5342.x + _S5342.y + _S5342.z;
    float _S5356 = d_32 * _S5355;
    float _S5357 = _S5333 * _S5355;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5358;
    (&_S5358)->primal_0 = _S5327;
    (&_S5358)->differential_0 = _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5359;
    (&_S5359)->primal_0 = v1v0_4;
    (&_S5359)->differential_0 = _S5351;
    s_bwd_prop_dot_0(&_S5358, &_S5359, _S5356);
    float _S5360 = _S5348 + _S5344.x + _S5344.y + _S5344.z;
    float _S5361 = d_32 * _S5360;
    float _S5362 = _S5331 * _S5360;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5363;
    (&_S5363)->primal_0 = _S5330;
    (&_S5363)->differential_0 = _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5364;
    (&_S5364)->primal_0 = v2v0_4;
    (&_S5364)->differential_0 = _S5351;
    s_bwd_prop_dot_0(&_S5363, &_S5364, _S5361);
    float3  _S5365 = - _S5363.differential_0;
    float _S5366 = - ((_S5350 + _S5357 + _S5362) / _S5329);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5367;
    (&_S5367)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5367)->differential_0 = _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5368;
    (&_S5368)->primal_0 = _S5326;
    (&_S5368)->differential_0 = _S5351;
    s_bwd_prop_dot_0(&_S5367, &_S5368, _S5366);
    float3  _S5369 = _S5358.differential_0 + _S5365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5370;
    (&_S5370)->primal_0 = rov0_4;
    (&_S5370)->differential_0 = _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5371;
    (&_S5371)->primal_0 = (*dpray_d_5).primal_0;
    (&_S5371)->differential_0 = _S5351;
    s_bwd_prop_cross_0(&_S5370, &_S5371, _S5369);
    float3  _S5372 = _S5354 + _S5368.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5373;
    (&_S5373)->primal_0 = v1v0_4;
    (&_S5373)->differential_0 = _S5351;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5374;
    (&_S5374)->primal_0 = v2v0_4;
    (&_S5374)->differential_0 = _S5351;
    s_bwd_prop_cross_0(&_S5373, &_S5374, _S5372);
    float3  _S5375 = _S5353.differential_0 + _S5370.differential_0;
    float3  _S5376 = _S5364.differential_0 + _S5374.differential_0;
    float3  _S5377 = _S5359.differential_0 + _S5373.differential_0;
    float3  _S5378 = - _S5375 + - _S5376 + - _S5377;
    float3  _S5379 = _S5367.differential_0 + _S5371.differential_0;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S5379;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S5375;
    FixedArray<float3 , 3>  _S5380;
    _S5380[int(0)] = _S5351;
    _S5380[int(1)] = _S5351;
    _S5380[int(2)] = _S5351;
    _S5380[int(2)] = _S5343;
    _S5380[int(1)] = _S5345;
    _S5380[int(0)] = _S5347;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S5380;
    FixedArray<float3 , 3>  _S5381;
    _S5381[int(0)] = _S5351;
    _S5381[int(1)] = _S5351;
    _S5381[int(2)] = _S5351;
    _S5381[int(2)] = _S5376;
    _S5381[int(0)] = _S5378;
    _S5381[int(1)] = _S5377;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S5381;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5382, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5383, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5384, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5385, float3  _S5386, float _S5387)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S5382, _S5383, _S5384, _S5385, _S5386, _S5387);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_9, FixedArray<float3 , 3>  * rgbs_6, float3  ray_o_10, float3  ray_d_10, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S5388 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5389 = { _S5388, _S5388, _S5388 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_9;
    (&dp_verts_1)->differential_0 = _S5389;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S5389;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_10;
    (&dp_ray_o_3)->differential_0 = _S5388;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_10;
    (&dp_ray_d_3)->differential_0 = _S5388;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

inline __device__ void projection_voxel_eval3d_persp(float3  pos_0, float size_0, FixedArray<float, 8>  * densities_0, FixedArray<float3 , 16>  * sh_coeffs_27, Matrix<float, 3, 3>  R_31, float3  t_33, float fx_36, float fy_36, float cx_31, float cy_31, FixedArray<float, 10>  * dist_coeffs_43, uint image_width_27, uint image_height_27, float near_plane_18, float far_plane_18, int4  * aabb_xyxy_18, float * depth_21, float3  * rgbs_7)
{
    float2  * _S5390;
    float2  * _S5391;
    float2  * _S5392;
    float2  * _S5393;
    float2  * _S5394;
    float2  * _S5395;
    float2  * _S5396;
    float2  * _S5397;
    bool _S5398;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_0;
        float3  _S5399 = mul_0(R_31, pos_0) + t_33;
        pos_c_0[int(0)] = _S5399;
        float _S5400 = _S5399.z;
        float _S5401 = (F32_min((far_plane_18), (_S5400)));
        float _S5402 = (F32_max((near_plane_18), (_S5400)));
        float3  _S5403 = mul_0(R_31, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 0.0f)) + t_33;
        pos_c_0[int(1)] = _S5403;
        float _S5404 = _S5403.z;
        float _S5405 = (F32_min((_S5401), (_S5404)));
        float _S5406 = (F32_max((_S5402), (_S5404)));
        float3  _S5407 = mul_0(R_31, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 0.0f)) + t_33;
        pos_c_0[int(2)] = _S5407;
        float _S5408 = _S5407.z;
        float _S5409 = (F32_min((_S5405), (_S5408)));
        float _S5410 = (F32_max((_S5406), (_S5408)));
        float3  _S5411 = mul_0(R_31, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 1.0f, 0.0f)) + t_33;
        pos_c_0[int(3)] = _S5411;
        float _S5412 = _S5411.z;
        float _S5413 = (F32_min((_S5409), (_S5412)));
        float _S5414 = (F32_max((_S5410), (_S5412)));
        float3  _S5415 = mul_0(R_31, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 0.0f, 1.0f)) + t_33;
        pos_c_0[int(4)] = _S5415;
        float _S5416 = _S5415.z;
        float _S5417 = (F32_min((_S5413), (_S5416)));
        float _S5418 = (F32_max((_S5414), (_S5416)));
        float3  _S5419 = mul_0(R_31, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 1.0f)) + t_33;
        pos_c_0[int(5)] = _S5419;
        float _S5420 = _S5419.z;
        float _S5421 = (F32_min((_S5417), (_S5420)));
        float _S5422 = (F32_max((_S5418), (_S5420)));
        float3  _S5423 = mul_0(R_31, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 1.0f)) + t_33;
        pos_c_0[int(6)] = _S5423;
        float _S5424 = _S5423.z;
        float _S5425 = (F32_min((_S5421), (_S5424)));
        float _S5426 = (F32_max((_S5422), (_S5424)));
        float3  _S5427 = mul_0(R_31, pos_0 + make_float3 (size_0)) + t_33;
        pos_c_0[int(7)] = _S5427;
        float _S5428 = _S5427.z;
        float _S5429 = (F32_min((_S5425), (_S5428)));
        float _S5430 = (F32_max((_S5426), (_S5428)));
        bool _S5431;
        if(_S5429 < near_plane_18)
        {
            _S5431 = true;
        }
        else
        {
            _S5431 = _S5430 > far_plane_18;
        }
        if(_S5431)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  mean_c_27 = mul_0(R_31, pos_0 + make_float3 (0.5f * size_0)) + t_33;
        FixedArray<float2 , 8>  uv_10;
        for(;;)
        {
            float3  _S5432 = pos_c_0[int(0)];
            _S5390 = &uv_10[int(0)];
            for(;;)
            {
                float _S5433 = _S5432.z;
                uv_10[int(0)] = float2 {_S5432.x, _S5432.y} / make_float2 (_S5433);
                if(_S5433 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5434 = camera_distortion_jac_0(uv_10[int(0)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5434)), ((F32_min((_S5434.rows[int(0)].x), (_S5434.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_85 = uv_10[int(0)].x;
                float v_85 = uv_10[int(0)].y;
                float r2_80 = u_85 * u_85 + v_85 * v_85;
                float2  _S5435 = uv_10[int(0)] * make_float2 (1.0f + r2_80 * ((*dist_coeffs_43)[int(0)] + r2_80 * ((*dist_coeffs_43)[int(1)] + r2_80 * ((*dist_coeffs_43)[int(2)] + r2_80 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_85 * v_85 + (*dist_coeffs_43)[int(5)] * (r2_80 + 2.0f * u_85 * u_85) + (*dist_coeffs_43)[int(6)] * r2_80, 2.0f * (*dist_coeffs_43)[int(5)] * u_85 * v_85 + (*dist_coeffs_43)[int(4)] * (r2_80 + 2.0f * v_85 * v_85) + (*dist_coeffs_43)[int(7)] * r2_80);
                float2  _S5436 = _S5435 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5435.x + (*dist_coeffs_43)[int(9)] * _S5435.y, 0.0f);
                uv_10[int(0)] = make_float2 (fx_36 * _S5436.x + cx_31, fy_36 * _S5436.y + cy_31);
                break;
            }
            bool all_valid_16 = true & (!_S5431);
            float3  _S5437 = pos_c_0[int(1)];
            _S5391 = &uv_10[int(1)];
            for(;;)
            {
                float _S5438 = _S5437.z;
                uv_10[int(1)] = float2 {_S5437.x, _S5437.y} / make_float2 (_S5438);
                if(_S5438 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5439 = camera_distortion_jac_0(uv_10[int(1)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5439)), ((F32_min((_S5439.rows[int(0)].x), (_S5439.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_86 = uv_10[int(1)].x;
                float v_86 = uv_10[int(1)].y;
                float r2_81 = u_86 * u_86 + v_86 * v_86;
                float2  _S5440 = uv_10[int(1)] * make_float2 (1.0f + r2_81 * ((*dist_coeffs_43)[int(0)] + r2_81 * ((*dist_coeffs_43)[int(1)] + r2_81 * ((*dist_coeffs_43)[int(2)] + r2_81 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_86 * v_86 + (*dist_coeffs_43)[int(5)] * (r2_81 + 2.0f * u_86 * u_86) + (*dist_coeffs_43)[int(6)] * r2_81, 2.0f * (*dist_coeffs_43)[int(5)] * u_86 * v_86 + (*dist_coeffs_43)[int(4)] * (r2_81 + 2.0f * v_86 * v_86) + (*dist_coeffs_43)[int(7)] * r2_81);
                float2  _S5441 = _S5440 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5440.x + (*dist_coeffs_43)[int(9)] * _S5440.y, 0.0f);
                uv_10[int(1)] = make_float2 (fx_36 * _S5441.x + cx_31, fy_36 * _S5441.y + cy_31);
                break;
            }
            bool all_valid_17 = all_valid_16 & (!_S5431);
            float3  _S5442 = pos_c_0[int(2)];
            _S5392 = &uv_10[int(2)];
            for(;;)
            {
                float _S5443 = _S5442.z;
                uv_10[int(2)] = float2 {_S5442.x, _S5442.y} / make_float2 (_S5443);
                if(_S5443 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5444 = camera_distortion_jac_0(uv_10[int(2)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5444)), ((F32_min((_S5444.rows[int(0)].x), (_S5444.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_87 = uv_10[int(2)].x;
                float v_87 = uv_10[int(2)].y;
                float r2_82 = u_87 * u_87 + v_87 * v_87;
                float2  _S5445 = uv_10[int(2)] * make_float2 (1.0f + r2_82 * ((*dist_coeffs_43)[int(0)] + r2_82 * ((*dist_coeffs_43)[int(1)] + r2_82 * ((*dist_coeffs_43)[int(2)] + r2_82 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_87 * v_87 + (*dist_coeffs_43)[int(5)] * (r2_82 + 2.0f * u_87 * u_87) + (*dist_coeffs_43)[int(6)] * r2_82, 2.0f * (*dist_coeffs_43)[int(5)] * u_87 * v_87 + (*dist_coeffs_43)[int(4)] * (r2_82 + 2.0f * v_87 * v_87) + (*dist_coeffs_43)[int(7)] * r2_82);
                float2  _S5446 = _S5445 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5445.x + (*dist_coeffs_43)[int(9)] * _S5445.y, 0.0f);
                uv_10[int(2)] = make_float2 (fx_36 * _S5446.x + cx_31, fy_36 * _S5446.y + cy_31);
                break;
            }
            bool all_valid_18 = all_valid_17 & (!_S5431);
            float3  _S5447 = pos_c_0[int(3)];
            _S5393 = &uv_10[int(3)];
            for(;;)
            {
                float _S5448 = _S5447.z;
                uv_10[int(3)] = float2 {_S5447.x, _S5447.y} / make_float2 (_S5448);
                if(_S5448 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5449 = camera_distortion_jac_0(uv_10[int(3)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5449)), ((F32_min((_S5449.rows[int(0)].x), (_S5449.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_88 = uv_10[int(3)].x;
                float v_88 = uv_10[int(3)].y;
                float r2_83 = u_88 * u_88 + v_88 * v_88;
                float2  _S5450 = uv_10[int(3)] * make_float2 (1.0f + r2_83 * ((*dist_coeffs_43)[int(0)] + r2_83 * ((*dist_coeffs_43)[int(1)] + r2_83 * ((*dist_coeffs_43)[int(2)] + r2_83 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_88 * v_88 + (*dist_coeffs_43)[int(5)] * (r2_83 + 2.0f * u_88 * u_88) + (*dist_coeffs_43)[int(6)] * r2_83, 2.0f * (*dist_coeffs_43)[int(5)] * u_88 * v_88 + (*dist_coeffs_43)[int(4)] * (r2_83 + 2.0f * v_88 * v_88) + (*dist_coeffs_43)[int(7)] * r2_83);
                float2  _S5451 = _S5450 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5450.x + (*dist_coeffs_43)[int(9)] * _S5450.y, 0.0f);
                uv_10[int(3)] = make_float2 (fx_36 * _S5451.x + cx_31, fy_36 * _S5451.y + cy_31);
                break;
            }
            bool all_valid_19 = all_valid_18 & (!_S5431);
            float3  _S5452 = pos_c_0[int(4)];
            _S5394 = &uv_10[int(4)];
            for(;;)
            {
                float _S5453 = _S5452.z;
                uv_10[int(4)] = float2 {_S5452.x, _S5452.y} / make_float2 (_S5453);
                if(_S5453 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5454 = camera_distortion_jac_0(uv_10[int(4)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5454)), ((F32_min((_S5454.rows[int(0)].x), (_S5454.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_89 = uv_10[int(4)].x;
                float v_89 = uv_10[int(4)].y;
                float r2_84 = u_89 * u_89 + v_89 * v_89;
                float2  _S5455 = uv_10[int(4)] * make_float2 (1.0f + r2_84 * ((*dist_coeffs_43)[int(0)] + r2_84 * ((*dist_coeffs_43)[int(1)] + r2_84 * ((*dist_coeffs_43)[int(2)] + r2_84 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_89 * v_89 + (*dist_coeffs_43)[int(5)] * (r2_84 + 2.0f * u_89 * u_89) + (*dist_coeffs_43)[int(6)] * r2_84, 2.0f * (*dist_coeffs_43)[int(5)] * u_89 * v_89 + (*dist_coeffs_43)[int(4)] * (r2_84 + 2.0f * v_89 * v_89) + (*dist_coeffs_43)[int(7)] * r2_84);
                float2  _S5456 = _S5455 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5455.x + (*dist_coeffs_43)[int(9)] * _S5455.y, 0.0f);
                uv_10[int(4)] = make_float2 (fx_36 * _S5456.x + cx_31, fy_36 * _S5456.y + cy_31);
                break;
            }
            bool all_valid_20 = all_valid_19 & (!_S5431);
            float3  _S5457 = pos_c_0[int(5)];
            _S5395 = &uv_10[int(5)];
            for(;;)
            {
                float _S5458 = _S5457.z;
                uv_10[int(5)] = float2 {_S5457.x, _S5457.y} / make_float2 (_S5458);
                if(_S5458 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5459 = camera_distortion_jac_0(uv_10[int(5)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5459)), ((F32_min((_S5459.rows[int(0)].x), (_S5459.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_90 = uv_10[int(5)].x;
                float v_90 = uv_10[int(5)].y;
                float r2_85 = u_90 * u_90 + v_90 * v_90;
                float2  _S5460 = uv_10[int(5)] * make_float2 (1.0f + r2_85 * ((*dist_coeffs_43)[int(0)] + r2_85 * ((*dist_coeffs_43)[int(1)] + r2_85 * ((*dist_coeffs_43)[int(2)] + r2_85 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_90 * v_90 + (*dist_coeffs_43)[int(5)] * (r2_85 + 2.0f * u_90 * u_90) + (*dist_coeffs_43)[int(6)] * r2_85, 2.0f * (*dist_coeffs_43)[int(5)] * u_90 * v_90 + (*dist_coeffs_43)[int(4)] * (r2_85 + 2.0f * v_90 * v_90) + (*dist_coeffs_43)[int(7)] * r2_85);
                float2  _S5461 = _S5460 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5460.x + (*dist_coeffs_43)[int(9)] * _S5460.y, 0.0f);
                uv_10[int(5)] = make_float2 (fx_36 * _S5461.x + cx_31, fy_36 * _S5461.y + cy_31);
                break;
            }
            bool all_valid_21 = all_valid_20 & (!_S5431);
            float3  _S5462 = pos_c_0[int(6)];
            _S5396 = &uv_10[int(6)];
            for(;;)
            {
                float _S5463 = _S5462.z;
                uv_10[int(6)] = float2 {_S5462.x, _S5462.y} / make_float2 (_S5463);
                if(_S5463 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5464 = camera_distortion_jac_0(uv_10[int(6)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5464)), ((F32_min((_S5464.rows[int(0)].x), (_S5464.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_91 = uv_10[int(6)].x;
                float v_91 = uv_10[int(6)].y;
                float r2_86 = u_91 * u_91 + v_91 * v_91;
                float2  _S5465 = uv_10[int(6)] * make_float2 (1.0f + r2_86 * ((*dist_coeffs_43)[int(0)] + r2_86 * ((*dist_coeffs_43)[int(1)] + r2_86 * ((*dist_coeffs_43)[int(2)] + r2_86 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_91 * v_91 + (*dist_coeffs_43)[int(5)] * (r2_86 + 2.0f * u_91 * u_91) + (*dist_coeffs_43)[int(6)] * r2_86, 2.0f * (*dist_coeffs_43)[int(5)] * u_91 * v_91 + (*dist_coeffs_43)[int(4)] * (r2_86 + 2.0f * v_91 * v_91) + (*dist_coeffs_43)[int(7)] * r2_86);
                float2  _S5466 = _S5465 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5465.x + (*dist_coeffs_43)[int(9)] * _S5465.y, 0.0f);
                uv_10[int(6)] = make_float2 (fx_36 * _S5466.x + cx_31, fy_36 * _S5466.y + cy_31);
                break;
            }
            bool all_valid_22 = all_valid_21 & (!_S5431);
            float3  _S5467 = pos_c_0[int(7)];
            _S5397 = &uv_10[int(7)];
            for(;;)
            {
                float _S5468 = _S5467.z;
                uv_10[int(7)] = float2 {_S5467.x, _S5467.y} / make_float2 (_S5468);
                if(_S5468 < 0.0f)
                {
                    _S5431 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5469 = camera_distortion_jac_0(uv_10[int(7)], dist_coeffs_43);
                    _S5431 = !((F32_min((determinant_0(_S5469)), ((F32_min((_S5469.rows[int(0)].x), (_S5469.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5431)
                {
                    break;
                }
                float u_92 = uv_10[int(7)].x;
                float v_92 = uv_10[int(7)].y;
                float r2_87 = u_92 * u_92 + v_92 * v_92;
                float2  _S5470 = uv_10[int(7)] * make_float2 (1.0f + r2_87 * ((*dist_coeffs_43)[int(0)] + r2_87 * ((*dist_coeffs_43)[int(1)] + r2_87 * ((*dist_coeffs_43)[int(2)] + r2_87 * (*dist_coeffs_43)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_43)[int(4)] * u_92 * v_92 + (*dist_coeffs_43)[int(5)] * (r2_87 + 2.0f * u_92 * u_92) + (*dist_coeffs_43)[int(6)] * r2_87, 2.0f * (*dist_coeffs_43)[int(5)] * u_92 * v_92 + (*dist_coeffs_43)[int(4)] * (r2_87 + 2.0f * v_92 * v_92) + (*dist_coeffs_43)[int(7)] * r2_87);
                float2  _S5471 = _S5470 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5470.x + (*dist_coeffs_43)[int(9)] * _S5470.y, 0.0f);
                uv_10[int(7)] = make_float2 (fx_36 * _S5471.x + cx_31, fy_36 * _S5471.y + cy_31);
                break;
            }
            _S5398 = all_valid_22 & (!_S5431);
            break;
        }
        if(!_S5398)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S5472 = (*_S5390).x;
        float _S5473 = (*_S5390).x;
        float _S5474 = (*_S5390).y;
        float _S5475 = (*_S5390).y;
        float _S5476 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5472), ((*_S5391).x)))), ((*_S5392).x)))), ((*_S5393).x)))), ((*_S5394).x)))), ((*_S5395).x)))), ((*_S5396).x)))), ((*_S5397).x)));
        float _S5477 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5473), ((*_S5391).x)))), ((*_S5392).x)))), ((*_S5393).x)))), ((*_S5394).x)))), ((*_S5395).x)))), ((*_S5396).x)))), ((*_S5397).x)));
        float _S5478 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5474), ((*_S5391).y)))), ((*_S5392).y)))), ((*_S5393).y)))), ((*_S5394).y)))), ((*_S5395).y)))), ((*_S5396).y)))), ((*_S5397).y)));
        float _S5479 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5475), ((*_S5391).y)))), ((*_S5392).y)))), ((*_S5393).y)))), ((*_S5394).y)))), ((*_S5395).y)))), ((*_S5396).y)))), ((*_S5397).y)));
        if(_S5476 <= 0.0f)
        {
            _S5431 = true;
        }
        else
        {
            _S5431 = _S5477 >= float(image_width_27);
        }
        if(_S5431)
        {
            _S5431 = true;
        }
        else
        {
            _S5431 = _S5478 <= 0.0f;
        }
        if(_S5431)
        {
            _S5431 = true;
        }
        else
        {
            _S5431 = _S5479 >= float(image_height_27);
        }
        if(_S5431)
        {
            _S5431 = true;
        }
        else
        {
            if(_S5429 <= 0.0f)
            {
                if(_S5477 <= 0.0f)
                {
                    _S5431 = _S5476 >= float(image_width_27);
                }
                else
                {
                    _S5431 = false;
                }
                if(_S5431)
                {
                    _S5431 = true;
                }
                else
                {
                    if(_S5479 <= 0.0f)
                    {
                        _S5431 = _S5478 >= float(image_width_27);
                    }
                    else
                    {
                        _S5431 = false;
                    }
                }
            }
            else
            {
                _S5431 = false;
            }
        }
        if(_S5431)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_18 = make_int4 (int((F32_floor((_S5477)))), int((F32_floor((_S5479)))), int((F32_ceil((_S5476)))), int((F32_ceil((_S5478)))));
        *depth_21 = mean_c_27.z;
        float3  _S5480 = mean_c_27 - - mul_0(transpose_0(R_31), t_33);
        float3  _S5481 = make_float3 (0.282094806432724f) * (*sh_coeffs_27)[int(0)];
        *rgbs_7 = _S5481;
        float _S5482 = _S5480.x;
        float _S5483 = _S5480.y;
        float _S5484 = _S5480.z;
        float norm_18 = (F32_sqrt((_S5482 * _S5482 + _S5483 * _S5483 + _S5484 * _S5484)));
        float x_63 = _S5482 / norm_18;
        float y_30 = _S5483 / norm_18;
        float z_27 = _S5484 / norm_18;
        float3  _S5485 = _S5481 + make_float3 (0.48860251903533936f) * (make_float3 (- y_30) * (*sh_coeffs_27)[int(1)] + make_float3 (z_27) * (*sh_coeffs_27)[int(2)] - make_float3 (x_63) * (*sh_coeffs_27)[int(3)]);
        *rgbs_7 = _S5485;
        float z2_60 = z_27 * z_27;
        float fTmp0B_27 = -1.09254848957061768f * z_27;
        float fC1_27 = x_63 * x_63 - y_30 * y_30;
        float fS1_27 = 2.0f * x_63 * y_30;
        float3  _S5486 = _S5485 + (make_float3 (0.54627424478530884f * fS1_27) * (*sh_coeffs_27)[int(4)] + make_float3 (fTmp0B_27 * y_30) * (*sh_coeffs_27)[int(5)] + make_float3 (0.94617468118667603f * z2_60 - 0.31539157032966614f) * (*sh_coeffs_27)[int(6)] + make_float3 (fTmp0B_27 * x_63) * (*sh_coeffs_27)[int(7)] + make_float3 (0.54627424478530884f * fC1_27) * (*sh_coeffs_27)[int(8)]);
        *rgbs_7 = _S5486;
        float fTmp0C_27 = -2.28522896766662598f * z2_60 + 0.4570457935333252f;
        float fTmp1B_27 = 1.44530570507049561f * z_27;
        *rgbs_7 = max_0(_S5486 + (make_float3 (-0.59004360437393188f * (x_63 * fS1_27 + y_30 * fC1_27)) * (*sh_coeffs_27)[int(9)] + make_float3 (fTmp1B_27 * fS1_27) * (*sh_coeffs_27)[int(10)] + make_float3 (fTmp0C_27 * y_30) * (*sh_coeffs_27)[int(11)] + make_float3 (z_27 * (1.86588168144226074f * z2_60 - 1.11952900886535645f)) * (*sh_coeffs_27)[int(12)] + make_float3 (fTmp0C_27 * x_63) * (*sh_coeffs_27)[int(13)] + make_float3 (fTmp1B_27 * fC1_27) * (*sh_coeffs_27)[int(14)] + make_float3 (-0.59004360437393188f * (x_63 * fC1_27 - y_30 * fS1_27)) * (*sh_coeffs_27)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye(float3  pos_1, float size_1, FixedArray<float, 8>  * densities_1, FixedArray<float3 , 16>  * sh_coeffs_28, Matrix<float, 3, 3>  R_32, float3  t_34, float fx_37, float fy_37, float cx_32, float cy_32, FixedArray<float, 10>  * dist_coeffs_44, uint image_width_28, uint image_height_28, float near_plane_19, float far_plane_19, int4  * aabb_xyxy_19, float * depth_22, float3  * rgbs_8)
{
    float2  * _S5487;
    bool _S5488;
    float2  * _S5489;
    bool _S5490;
    float2  * _S5491;
    bool _S5492;
    float2  * _S5493;
    bool _S5494;
    float2  * _S5495;
    bool _S5496;
    float2  * _S5497;
    bool _S5498;
    float2  * _S5499;
    bool _S5500;
    float2  * _S5501;
    bool _S5502;
    bool _S5503;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_1;
        float3  _S5504 = mul_0(R_32, pos_1) + t_34;
        pos_c_1[int(0)] = _S5504;
        float _S5505 = length_1(_S5504);
        float _S5506 = (F32_min((far_plane_19), (_S5505)));
        float _S5507 = (F32_max((near_plane_19), (_S5505)));
        float3  _S5508 = mul_0(R_32, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 0.0f)) + t_34;
        pos_c_1[int(1)] = _S5508;
        float _S5509 = length_1(_S5508);
        float _S5510 = (F32_min((_S5506), (_S5509)));
        float _S5511 = (F32_max((_S5507), (_S5509)));
        float3  _S5512 = mul_0(R_32, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 0.0f)) + t_34;
        pos_c_1[int(2)] = _S5512;
        float _S5513 = length_1(_S5512);
        float _S5514 = (F32_min((_S5510), (_S5513)));
        float _S5515 = (F32_max((_S5511), (_S5513)));
        float3  _S5516 = mul_0(R_32, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 1.0f, 0.0f)) + t_34;
        pos_c_1[int(3)] = _S5516;
        float _S5517 = length_1(_S5516);
        float _S5518 = (F32_min((_S5514), (_S5517)));
        float _S5519 = (F32_max((_S5515), (_S5517)));
        float3  _S5520 = mul_0(R_32, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 0.0f, 1.0f)) + t_34;
        pos_c_1[int(4)] = _S5520;
        float _S5521 = length_1(_S5520);
        float _S5522 = (F32_min((_S5518), (_S5521)));
        float _S5523 = (F32_max((_S5519), (_S5521)));
        float3  _S5524 = mul_0(R_32, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 1.0f)) + t_34;
        pos_c_1[int(5)] = _S5524;
        float _S5525 = length_1(_S5524);
        float _S5526 = (F32_min((_S5522), (_S5525)));
        float _S5527 = (F32_max((_S5523), (_S5525)));
        float3  _S5528 = mul_0(R_32, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 1.0f)) + t_34;
        pos_c_1[int(6)] = _S5528;
        float _S5529 = length_1(_S5528);
        float _S5530 = (F32_min((_S5526), (_S5529)));
        float _S5531 = (F32_max((_S5527), (_S5529)));
        float3  _S5532 = mul_0(R_32, pos_1 + make_float3 (size_1)) + t_34;
        pos_c_1[int(7)] = _S5532;
        float _S5533 = length_1(_S5532);
        float _S5534 = (F32_min((_S5530), (_S5533)));
        float _S5535 = (F32_max((_S5531), (_S5533)));
        bool _S5536;
        if(_S5534 < near_plane_19)
        {
            _S5536 = true;
        }
        else
        {
            _S5536 = _S5535 > far_plane_19;
        }
        if(_S5536)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  mean_c_28 = mul_0(R_32, pos_1 + make_float3 (0.5f * size_1)) + t_34;
        FixedArray<float2 , 8>  uv_11;
        for(;;)
        {
            float k_13;
            float3  _S5537 = pos_c_1[int(0)];
            _S5487 = &uv_11[int(0)];
            for(;;)
            {
                float2  _S5538 = float2 {_S5537.x, _S5537.y};
                float r_36 = length_0(_S5538);
                float _S5539 = _S5537.z;
                float theta_32 = (F32_atan2((r_36), (_S5539)));
                if(theta_32 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_32 * theta_32 / 3.0f) / _S5539;
                }
                else
                {
                    k_13 = theta_32 / r_36;
                }
                float2  _S5540 = _S5538 * make_float2 (k_13);
                uv_11[int(0)] = _S5540;
                Matrix<float, 2, 2>  _S5541 = camera_distortion_jac_0(_S5540, dist_coeffs_44);
                bool _S5542 = !((F32_min((determinant_0(_S5541)), ((F32_min((_S5541.rows[int(0)].x), (_S5541.rows[int(1)].y)))))) > 0.0f);
                _S5488 = _S5542;
                if(_S5542)
                {
                    break;
                }
                float u_93 = uv_11[int(0)].x;
                float v_93 = uv_11[int(0)].y;
                float r2_88 = u_93 * u_93 + v_93 * v_93;
                float2  _S5543 = uv_11[int(0)] * make_float2 (1.0f + r2_88 * ((*dist_coeffs_44)[int(0)] + r2_88 * ((*dist_coeffs_44)[int(1)] + r2_88 * ((*dist_coeffs_44)[int(2)] + r2_88 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_93 * v_93 + (*dist_coeffs_44)[int(5)] * (r2_88 + 2.0f * u_93 * u_93) + (*dist_coeffs_44)[int(6)] * r2_88, 2.0f * (*dist_coeffs_44)[int(5)] * u_93 * v_93 + (*dist_coeffs_44)[int(4)] * (r2_88 + 2.0f * v_93 * v_93) + (*dist_coeffs_44)[int(7)] * r2_88);
                float2  _S5544 = _S5543 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5543.x + (*dist_coeffs_44)[int(9)] * _S5543.y, 0.0f);
                uv_11[int(0)] = make_float2 (fx_37 * _S5544.x + cx_32, fy_37 * _S5544.y + cy_32);
                break;
            }
            bool all_valid_23 = true & (!_S5488);
            float3  _S5545 = pos_c_1[int(1)];
            _S5489 = &uv_11[int(1)];
            for(;;)
            {
                float2  _S5546 = float2 {_S5545.x, _S5545.y};
                float r_37 = length_0(_S5546);
                float _S5547 = _S5545.z;
                float theta_33 = (F32_atan2((r_37), (_S5547)));
                if(theta_33 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_33 * theta_33 / 3.0f) / _S5547;
                }
                else
                {
                    k_13 = theta_33 / r_37;
                }
                float2  _S5548 = _S5546 * make_float2 (k_13);
                uv_11[int(1)] = _S5548;
                Matrix<float, 2, 2>  _S5549 = camera_distortion_jac_0(_S5548, dist_coeffs_44);
                bool _S5550 = !((F32_min((determinant_0(_S5549)), ((F32_min((_S5549.rows[int(0)].x), (_S5549.rows[int(1)].y)))))) > 0.0f);
                _S5490 = _S5550;
                if(_S5550)
                {
                    break;
                }
                float u_94 = uv_11[int(1)].x;
                float v_94 = uv_11[int(1)].y;
                float r2_89 = u_94 * u_94 + v_94 * v_94;
                float2  _S5551 = uv_11[int(1)] * make_float2 (1.0f + r2_89 * ((*dist_coeffs_44)[int(0)] + r2_89 * ((*dist_coeffs_44)[int(1)] + r2_89 * ((*dist_coeffs_44)[int(2)] + r2_89 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_94 * v_94 + (*dist_coeffs_44)[int(5)] * (r2_89 + 2.0f * u_94 * u_94) + (*dist_coeffs_44)[int(6)] * r2_89, 2.0f * (*dist_coeffs_44)[int(5)] * u_94 * v_94 + (*dist_coeffs_44)[int(4)] * (r2_89 + 2.0f * v_94 * v_94) + (*dist_coeffs_44)[int(7)] * r2_89);
                float2  _S5552 = _S5551 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5551.x + (*dist_coeffs_44)[int(9)] * _S5551.y, 0.0f);
                uv_11[int(1)] = make_float2 (fx_37 * _S5552.x + cx_32, fy_37 * _S5552.y + cy_32);
                break;
            }
            bool all_valid_24 = all_valid_23 & (!_S5490);
            float3  _S5553 = pos_c_1[int(2)];
            _S5491 = &uv_11[int(2)];
            for(;;)
            {
                float2  _S5554 = float2 {_S5553.x, _S5553.y};
                float r_38 = length_0(_S5554);
                float _S5555 = _S5553.z;
                float theta_34 = (F32_atan2((r_38), (_S5555)));
                if(theta_34 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_34 * theta_34 / 3.0f) / _S5555;
                }
                else
                {
                    k_13 = theta_34 / r_38;
                }
                float2  _S5556 = _S5554 * make_float2 (k_13);
                uv_11[int(2)] = _S5556;
                Matrix<float, 2, 2>  _S5557 = camera_distortion_jac_0(_S5556, dist_coeffs_44);
                bool _S5558 = !((F32_min((determinant_0(_S5557)), ((F32_min((_S5557.rows[int(0)].x), (_S5557.rows[int(1)].y)))))) > 0.0f);
                _S5492 = _S5558;
                if(_S5558)
                {
                    break;
                }
                float u_95 = uv_11[int(2)].x;
                float v_95 = uv_11[int(2)].y;
                float r2_90 = u_95 * u_95 + v_95 * v_95;
                float2  _S5559 = uv_11[int(2)] * make_float2 (1.0f + r2_90 * ((*dist_coeffs_44)[int(0)] + r2_90 * ((*dist_coeffs_44)[int(1)] + r2_90 * ((*dist_coeffs_44)[int(2)] + r2_90 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_95 * v_95 + (*dist_coeffs_44)[int(5)] * (r2_90 + 2.0f * u_95 * u_95) + (*dist_coeffs_44)[int(6)] * r2_90, 2.0f * (*dist_coeffs_44)[int(5)] * u_95 * v_95 + (*dist_coeffs_44)[int(4)] * (r2_90 + 2.0f * v_95 * v_95) + (*dist_coeffs_44)[int(7)] * r2_90);
                float2  _S5560 = _S5559 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5559.x + (*dist_coeffs_44)[int(9)] * _S5559.y, 0.0f);
                uv_11[int(2)] = make_float2 (fx_37 * _S5560.x + cx_32, fy_37 * _S5560.y + cy_32);
                break;
            }
            bool all_valid_25 = all_valid_24 & (!_S5492);
            float3  _S5561 = pos_c_1[int(3)];
            _S5493 = &uv_11[int(3)];
            for(;;)
            {
                float2  _S5562 = float2 {_S5561.x, _S5561.y};
                float r_39 = length_0(_S5562);
                float _S5563 = _S5561.z;
                float theta_35 = (F32_atan2((r_39), (_S5563)));
                if(theta_35 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_35 * theta_35 / 3.0f) / _S5563;
                }
                else
                {
                    k_13 = theta_35 / r_39;
                }
                float2  _S5564 = _S5562 * make_float2 (k_13);
                uv_11[int(3)] = _S5564;
                Matrix<float, 2, 2>  _S5565 = camera_distortion_jac_0(_S5564, dist_coeffs_44);
                bool _S5566 = !((F32_min((determinant_0(_S5565)), ((F32_min((_S5565.rows[int(0)].x), (_S5565.rows[int(1)].y)))))) > 0.0f);
                _S5494 = _S5566;
                if(_S5566)
                {
                    break;
                }
                float u_96 = uv_11[int(3)].x;
                float v_96 = uv_11[int(3)].y;
                float r2_91 = u_96 * u_96 + v_96 * v_96;
                float2  _S5567 = uv_11[int(3)] * make_float2 (1.0f + r2_91 * ((*dist_coeffs_44)[int(0)] + r2_91 * ((*dist_coeffs_44)[int(1)] + r2_91 * ((*dist_coeffs_44)[int(2)] + r2_91 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_96 * v_96 + (*dist_coeffs_44)[int(5)] * (r2_91 + 2.0f * u_96 * u_96) + (*dist_coeffs_44)[int(6)] * r2_91, 2.0f * (*dist_coeffs_44)[int(5)] * u_96 * v_96 + (*dist_coeffs_44)[int(4)] * (r2_91 + 2.0f * v_96 * v_96) + (*dist_coeffs_44)[int(7)] * r2_91);
                float2  _S5568 = _S5567 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5567.x + (*dist_coeffs_44)[int(9)] * _S5567.y, 0.0f);
                uv_11[int(3)] = make_float2 (fx_37 * _S5568.x + cx_32, fy_37 * _S5568.y + cy_32);
                break;
            }
            bool all_valid_26 = all_valid_25 & (!_S5494);
            float3  _S5569 = pos_c_1[int(4)];
            _S5495 = &uv_11[int(4)];
            for(;;)
            {
                float2  _S5570 = float2 {_S5569.x, _S5569.y};
                float r_40 = length_0(_S5570);
                float _S5571 = _S5569.z;
                float theta_36 = (F32_atan2((r_40), (_S5571)));
                if(theta_36 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_36 * theta_36 / 3.0f) / _S5571;
                }
                else
                {
                    k_13 = theta_36 / r_40;
                }
                float2  _S5572 = _S5570 * make_float2 (k_13);
                uv_11[int(4)] = _S5572;
                Matrix<float, 2, 2>  _S5573 = camera_distortion_jac_0(_S5572, dist_coeffs_44);
                bool _S5574 = !((F32_min((determinant_0(_S5573)), ((F32_min((_S5573.rows[int(0)].x), (_S5573.rows[int(1)].y)))))) > 0.0f);
                _S5496 = _S5574;
                if(_S5574)
                {
                    break;
                }
                float u_97 = uv_11[int(4)].x;
                float v_97 = uv_11[int(4)].y;
                float r2_92 = u_97 * u_97 + v_97 * v_97;
                float2  _S5575 = uv_11[int(4)] * make_float2 (1.0f + r2_92 * ((*dist_coeffs_44)[int(0)] + r2_92 * ((*dist_coeffs_44)[int(1)] + r2_92 * ((*dist_coeffs_44)[int(2)] + r2_92 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_97 * v_97 + (*dist_coeffs_44)[int(5)] * (r2_92 + 2.0f * u_97 * u_97) + (*dist_coeffs_44)[int(6)] * r2_92, 2.0f * (*dist_coeffs_44)[int(5)] * u_97 * v_97 + (*dist_coeffs_44)[int(4)] * (r2_92 + 2.0f * v_97 * v_97) + (*dist_coeffs_44)[int(7)] * r2_92);
                float2  _S5576 = _S5575 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5575.x + (*dist_coeffs_44)[int(9)] * _S5575.y, 0.0f);
                uv_11[int(4)] = make_float2 (fx_37 * _S5576.x + cx_32, fy_37 * _S5576.y + cy_32);
                break;
            }
            bool all_valid_27 = all_valid_26 & (!_S5496);
            float3  _S5577 = pos_c_1[int(5)];
            _S5497 = &uv_11[int(5)];
            for(;;)
            {
                float2  _S5578 = float2 {_S5577.x, _S5577.y};
                float r_41 = length_0(_S5578);
                float _S5579 = _S5577.z;
                float theta_37 = (F32_atan2((r_41), (_S5579)));
                if(theta_37 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_37 * theta_37 / 3.0f) / _S5579;
                }
                else
                {
                    k_13 = theta_37 / r_41;
                }
                float2  _S5580 = _S5578 * make_float2 (k_13);
                uv_11[int(5)] = _S5580;
                Matrix<float, 2, 2>  _S5581 = camera_distortion_jac_0(_S5580, dist_coeffs_44);
                bool _S5582 = !((F32_min((determinant_0(_S5581)), ((F32_min((_S5581.rows[int(0)].x), (_S5581.rows[int(1)].y)))))) > 0.0f);
                _S5498 = _S5582;
                if(_S5582)
                {
                    break;
                }
                float u_98 = uv_11[int(5)].x;
                float v_98 = uv_11[int(5)].y;
                float r2_93 = u_98 * u_98 + v_98 * v_98;
                float2  _S5583 = uv_11[int(5)] * make_float2 (1.0f + r2_93 * ((*dist_coeffs_44)[int(0)] + r2_93 * ((*dist_coeffs_44)[int(1)] + r2_93 * ((*dist_coeffs_44)[int(2)] + r2_93 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_98 * v_98 + (*dist_coeffs_44)[int(5)] * (r2_93 + 2.0f * u_98 * u_98) + (*dist_coeffs_44)[int(6)] * r2_93, 2.0f * (*dist_coeffs_44)[int(5)] * u_98 * v_98 + (*dist_coeffs_44)[int(4)] * (r2_93 + 2.0f * v_98 * v_98) + (*dist_coeffs_44)[int(7)] * r2_93);
                float2  _S5584 = _S5583 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5583.x + (*dist_coeffs_44)[int(9)] * _S5583.y, 0.0f);
                uv_11[int(5)] = make_float2 (fx_37 * _S5584.x + cx_32, fy_37 * _S5584.y + cy_32);
                break;
            }
            bool all_valid_28 = all_valid_27 & (!_S5498);
            float3  _S5585 = pos_c_1[int(6)];
            _S5499 = &uv_11[int(6)];
            for(;;)
            {
                float2  _S5586 = float2 {_S5585.x, _S5585.y};
                float r_42 = length_0(_S5586);
                float _S5587 = _S5585.z;
                float theta_38 = (F32_atan2((r_42), (_S5587)));
                if(theta_38 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_38 * theta_38 / 3.0f) / _S5587;
                }
                else
                {
                    k_13 = theta_38 / r_42;
                }
                float2  _S5588 = _S5586 * make_float2 (k_13);
                uv_11[int(6)] = _S5588;
                Matrix<float, 2, 2>  _S5589 = camera_distortion_jac_0(_S5588, dist_coeffs_44);
                bool _S5590 = !((F32_min((determinant_0(_S5589)), ((F32_min((_S5589.rows[int(0)].x), (_S5589.rows[int(1)].y)))))) > 0.0f);
                _S5500 = _S5590;
                if(_S5590)
                {
                    break;
                }
                float u_99 = uv_11[int(6)].x;
                float v_99 = uv_11[int(6)].y;
                float r2_94 = u_99 * u_99 + v_99 * v_99;
                float2  _S5591 = uv_11[int(6)] * make_float2 (1.0f + r2_94 * ((*dist_coeffs_44)[int(0)] + r2_94 * ((*dist_coeffs_44)[int(1)] + r2_94 * ((*dist_coeffs_44)[int(2)] + r2_94 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_99 * v_99 + (*dist_coeffs_44)[int(5)] * (r2_94 + 2.0f * u_99 * u_99) + (*dist_coeffs_44)[int(6)] * r2_94, 2.0f * (*dist_coeffs_44)[int(5)] * u_99 * v_99 + (*dist_coeffs_44)[int(4)] * (r2_94 + 2.0f * v_99 * v_99) + (*dist_coeffs_44)[int(7)] * r2_94);
                float2  _S5592 = _S5591 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5591.x + (*dist_coeffs_44)[int(9)] * _S5591.y, 0.0f);
                uv_11[int(6)] = make_float2 (fx_37 * _S5592.x + cx_32, fy_37 * _S5592.y + cy_32);
                break;
            }
            bool all_valid_29 = all_valid_28 & (!_S5500);
            float3  _S5593 = pos_c_1[int(7)];
            _S5501 = &uv_11[int(7)];
            for(;;)
            {
                float2  _S5594 = float2 {_S5593.x, _S5593.y};
                float r_43 = length_0(_S5594);
                float _S5595 = _S5593.z;
                float theta_39 = (F32_atan2((r_43), (_S5595)));
                if(theta_39 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_39 * theta_39 / 3.0f) / _S5595;
                }
                else
                {
                    k_13 = theta_39 / r_43;
                }
                float2  _S5596 = _S5594 * make_float2 (k_13);
                uv_11[int(7)] = _S5596;
                Matrix<float, 2, 2>  _S5597 = camera_distortion_jac_0(_S5596, dist_coeffs_44);
                bool _S5598 = !((F32_min((determinant_0(_S5597)), ((F32_min((_S5597.rows[int(0)].x), (_S5597.rows[int(1)].y)))))) > 0.0f);
                _S5502 = _S5598;
                if(_S5598)
                {
                    break;
                }
                float u_100 = uv_11[int(7)].x;
                float v_100 = uv_11[int(7)].y;
                float r2_95 = u_100 * u_100 + v_100 * v_100;
                float2  _S5599 = uv_11[int(7)] * make_float2 (1.0f + r2_95 * ((*dist_coeffs_44)[int(0)] + r2_95 * ((*dist_coeffs_44)[int(1)] + r2_95 * ((*dist_coeffs_44)[int(2)] + r2_95 * (*dist_coeffs_44)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_44)[int(4)] * u_100 * v_100 + (*dist_coeffs_44)[int(5)] * (r2_95 + 2.0f * u_100 * u_100) + (*dist_coeffs_44)[int(6)] * r2_95, 2.0f * (*dist_coeffs_44)[int(5)] * u_100 * v_100 + (*dist_coeffs_44)[int(4)] * (r2_95 + 2.0f * v_100 * v_100) + (*dist_coeffs_44)[int(7)] * r2_95);
                float2  _S5600 = _S5599 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5599.x + (*dist_coeffs_44)[int(9)] * _S5599.y, 0.0f);
                uv_11[int(7)] = make_float2 (fx_37 * _S5600.x + cx_32, fy_37 * _S5600.y + cy_32);
                break;
            }
            _S5503 = all_valid_29 & (!_S5502);
            break;
        }
        if(!_S5503)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S5601 = (*_S5487).x;
        float _S5602 = (*_S5487).x;
        float _S5603 = (*_S5487).y;
        float _S5604 = (*_S5487).y;
        float _S5605 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5601), ((*_S5489).x)))), ((*_S5491).x)))), ((*_S5493).x)))), ((*_S5495).x)))), ((*_S5497).x)))), ((*_S5499).x)))), ((*_S5501).x)));
        float _S5606 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5602), ((*_S5489).x)))), ((*_S5491).x)))), ((*_S5493).x)))), ((*_S5495).x)))), ((*_S5497).x)))), ((*_S5499).x)))), ((*_S5501).x)));
        float _S5607 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5603), ((*_S5489).y)))), ((*_S5491).y)))), ((*_S5493).y)))), ((*_S5495).y)))), ((*_S5497).y)))), ((*_S5499).y)))), ((*_S5501).y)));
        float _S5608 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5604), ((*_S5489).y)))), ((*_S5491).y)))), ((*_S5493).y)))), ((*_S5495).y)))), ((*_S5497).y)))), ((*_S5499).y)))), ((*_S5501).y)));
        if(_S5605 <= 0.0f)
        {
            _S5536 = true;
        }
        else
        {
            _S5536 = _S5606 >= float(image_width_28);
        }
        if(_S5536)
        {
            _S5536 = true;
        }
        else
        {
            _S5536 = _S5607 <= 0.0f;
        }
        if(_S5536)
        {
            _S5536 = true;
        }
        else
        {
            _S5536 = _S5608 >= float(image_height_28);
        }
        if(_S5536)
        {
            _S5536 = true;
        }
        else
        {
            if(_S5534 <= 0.0f)
            {
                if(_S5606 <= 0.0f)
                {
                    _S5536 = _S5605 >= float(image_width_28);
                }
                else
                {
                    _S5536 = false;
                }
                if(_S5536)
                {
                    _S5536 = true;
                }
                else
                {
                    if(_S5608 <= 0.0f)
                    {
                        _S5536 = _S5607 >= float(image_width_28);
                    }
                    else
                    {
                        _S5536 = false;
                    }
                }
            }
            else
            {
                _S5536 = false;
            }
        }
        if(_S5536)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_19 = make_int4 (int((F32_floor((_S5606)))), int((F32_floor((_S5608)))), int((F32_ceil((_S5605)))), int((F32_ceil((_S5607)))));
        *depth_22 = mean_c_28.z;
        float3  _S5609 = mean_c_28 - - mul_0(transpose_0(R_32), t_34);
        float3  _S5610 = make_float3 (0.282094806432724f) * (*sh_coeffs_28)[int(0)];
        *rgbs_8 = _S5610;
        float _S5611 = _S5609.x;
        float _S5612 = _S5609.y;
        float _S5613 = _S5609.z;
        float norm_19 = (F32_sqrt((_S5611 * _S5611 + _S5612 * _S5612 + _S5613 * _S5613)));
        float x_64 = _S5611 / norm_19;
        float y_31 = _S5612 / norm_19;
        float z_28 = _S5613 / norm_19;
        float3  _S5614 = _S5610 + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * (*sh_coeffs_28)[int(1)] + make_float3 (z_28) * (*sh_coeffs_28)[int(2)] - make_float3 (x_64) * (*sh_coeffs_28)[int(3)]);
        *rgbs_8 = _S5614;
        float z2_61 = z_28 * z_28;
        float fTmp0B_28 = -1.09254848957061768f * z_28;
        float fC1_28 = x_64 * x_64 - y_31 * y_31;
        float fS1_28 = 2.0f * x_64 * y_31;
        float3  _S5615 = _S5614 + (make_float3 (0.54627424478530884f * fS1_28) * (*sh_coeffs_28)[int(4)] + make_float3 (fTmp0B_28 * y_31) * (*sh_coeffs_28)[int(5)] + make_float3 (0.94617468118667603f * z2_61 - 0.31539157032966614f) * (*sh_coeffs_28)[int(6)] + make_float3 (fTmp0B_28 * x_64) * (*sh_coeffs_28)[int(7)] + make_float3 (0.54627424478530884f * fC1_28) * (*sh_coeffs_28)[int(8)]);
        *rgbs_8 = _S5615;
        float fTmp0C_28 = -2.28522896766662598f * z2_61 + 0.4570457935333252f;
        float fTmp1B_28 = 1.44530570507049561f * z_28;
        *rgbs_8 = max_0(_S5615 + (make_float3 (-0.59004360437393188f * (x_64 * fS1_28 + y_31 * fC1_28)) * (*sh_coeffs_28)[int(9)] + make_float3 (fTmp1B_28 * fS1_28) * (*sh_coeffs_28)[int(10)] + make_float3 (fTmp0C_28 * y_31) * (*sh_coeffs_28)[int(11)] + make_float3 (z_28 * (1.86588168144226074f * z2_61 - 1.11952900886535645f)) * (*sh_coeffs_28)[int(12)] + make_float3 (fTmp0C_28 * x_64) * (*sh_coeffs_28)[int(13)] + make_float3 (fTmp1B_28 * fC1_28) * (*sh_coeffs_28)[int(14)] + make_float3 (-0.59004360437393188f * (x_64 * fC1_28 - y_31 * fS1_28)) * (*sh_coeffs_28)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_voxel_eval3d_persp_differentiable(float3  pos_2, float size_2, FixedArray<float, 8>  * densities_2, FixedArray<float3 , 16>  * sh_coeffs_29, Matrix<float, 3, 3>  R_33, float3  t_35, float fx_38, float fy_38, float cx_33, float cy_33, FixedArray<float, 10>  * dist_coeffs_45, uint image_width_29, uint image_height_29, float near_plane_20, float far_plane_20, int4  * aabb_xyxy_20, float * depth_23, float3  * rgbs_9)
{
    FixedArray<float3 , 8>  pos_c_2;
    float3  _S5616 = mul_0(R_33, pos_2) + t_35;
    pos_c_2[int(0)] = _S5616;
    float3  _S5617 = mul_0(R_33, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 0.0f)) + t_35;
    pos_c_2[int(1)] = _S5617;
    float3  _S5618 = mul_0(R_33, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 0.0f)) + t_35;
    pos_c_2[int(2)] = _S5618;
    float3  _S5619 = mul_0(R_33, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 1.0f, 0.0f)) + t_35;
    pos_c_2[int(3)] = _S5619;
    float3  _S5620 = mul_0(R_33, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 0.0f, 1.0f)) + t_35;
    pos_c_2[int(4)] = _S5620;
    float3  _S5621 = mul_0(R_33, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 1.0f)) + t_35;
    pos_c_2[int(5)] = _S5621;
    float3  _S5622 = mul_0(R_33, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 1.0f)) + t_35;
    pos_c_2[int(6)] = _S5622;
    float3  _S5623 = mul_0(R_33, pos_2 + make_float3 (size_2)) + t_35;
    pos_c_2[int(7)] = _S5623;
    float3  mean_c_29 = mul_0(R_33, pos_2 + make_float3 (0.5f * size_2)) + t_35;
    FixedArray<float2 , 8>  uv_12;
    float2  _S5624 = float2 {_S5616.x, _S5616.y} / make_float2 (_S5616.z);
    float u_101 = _S5624.x;
    float v_101 = _S5624.y;
    float r2_96 = u_101 * u_101 + v_101 * v_101;
    float _S5625 = 2.0f * (*dist_coeffs_45)[int(4)];
    float _S5626 = 2.0f * (*dist_coeffs_45)[int(5)];
    float2  _S5627 = _S5624 * make_float2 (1.0f + r2_96 * ((*dist_coeffs_45)[int(0)] + r2_96 * ((*dist_coeffs_45)[int(1)] + r2_96 * ((*dist_coeffs_45)[int(2)] + r2_96 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_101 * v_101 + (*dist_coeffs_45)[int(5)] * (r2_96 + 2.0f * u_101 * u_101) + (*dist_coeffs_45)[int(6)] * r2_96, _S5626 * u_101 * v_101 + (*dist_coeffs_45)[int(4)] * (r2_96 + 2.0f * v_101 * v_101) + (*dist_coeffs_45)[int(7)] * r2_96);
    float2  _S5628 = _S5627 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5627.x + (*dist_coeffs_45)[int(9)] * _S5627.y, 0.0f);
    float _S5629 = fx_38 * _S5628.x + cx_33;
    float _S5630 = fy_38 * _S5628.y + cy_33;
    uv_12[int(0)] = make_float2 (_S5629, _S5630);
    float2  _S5631 = float2 {_S5617.x, _S5617.y} / make_float2 (_S5617.z);
    float u_102 = _S5631.x;
    float v_102 = _S5631.y;
    float r2_97 = u_102 * u_102 + v_102 * v_102;
    float2  _S5632 = _S5631 * make_float2 (1.0f + r2_97 * ((*dist_coeffs_45)[int(0)] + r2_97 * ((*dist_coeffs_45)[int(1)] + r2_97 * ((*dist_coeffs_45)[int(2)] + r2_97 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_102 * v_102 + (*dist_coeffs_45)[int(5)] * (r2_97 + 2.0f * u_102 * u_102) + (*dist_coeffs_45)[int(6)] * r2_97, _S5626 * u_102 * v_102 + (*dist_coeffs_45)[int(4)] * (r2_97 + 2.0f * v_102 * v_102) + (*dist_coeffs_45)[int(7)] * r2_97);
    float2  _S5633 = _S5632 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5632.x + (*dist_coeffs_45)[int(9)] * _S5632.y, 0.0f);
    float _S5634 = fx_38 * _S5633.x + cx_33;
    float _S5635 = fy_38 * _S5633.y + cy_33;
    uv_12[int(1)] = make_float2 (_S5634, _S5635);
    float2  _S5636 = float2 {_S5618.x, _S5618.y} / make_float2 (_S5618.z);
    float u_103 = _S5636.x;
    float v_103 = _S5636.y;
    float r2_98 = u_103 * u_103 + v_103 * v_103;
    float2  _S5637 = _S5636 * make_float2 (1.0f + r2_98 * ((*dist_coeffs_45)[int(0)] + r2_98 * ((*dist_coeffs_45)[int(1)] + r2_98 * ((*dist_coeffs_45)[int(2)] + r2_98 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_103 * v_103 + (*dist_coeffs_45)[int(5)] * (r2_98 + 2.0f * u_103 * u_103) + (*dist_coeffs_45)[int(6)] * r2_98, _S5626 * u_103 * v_103 + (*dist_coeffs_45)[int(4)] * (r2_98 + 2.0f * v_103 * v_103) + (*dist_coeffs_45)[int(7)] * r2_98);
    float2  _S5638 = _S5637 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5637.x + (*dist_coeffs_45)[int(9)] * _S5637.y, 0.0f);
    float _S5639 = fx_38 * _S5638.x + cx_33;
    float _S5640 = fy_38 * _S5638.y + cy_33;
    uv_12[int(2)] = make_float2 (_S5639, _S5640);
    float2  _S5641 = float2 {_S5619.x, _S5619.y} / make_float2 (_S5619.z);
    float u_104 = _S5641.x;
    float v_104 = _S5641.y;
    float r2_99 = u_104 * u_104 + v_104 * v_104;
    float2  _S5642 = _S5641 * make_float2 (1.0f + r2_99 * ((*dist_coeffs_45)[int(0)] + r2_99 * ((*dist_coeffs_45)[int(1)] + r2_99 * ((*dist_coeffs_45)[int(2)] + r2_99 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_104 * v_104 + (*dist_coeffs_45)[int(5)] * (r2_99 + 2.0f * u_104 * u_104) + (*dist_coeffs_45)[int(6)] * r2_99, _S5626 * u_104 * v_104 + (*dist_coeffs_45)[int(4)] * (r2_99 + 2.0f * v_104 * v_104) + (*dist_coeffs_45)[int(7)] * r2_99);
    float2  _S5643 = _S5642 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5642.x + (*dist_coeffs_45)[int(9)] * _S5642.y, 0.0f);
    float _S5644 = fx_38 * _S5643.x + cx_33;
    float _S5645 = fy_38 * _S5643.y + cy_33;
    uv_12[int(3)] = make_float2 (_S5644, _S5645);
    float2  _S5646 = float2 {_S5620.x, _S5620.y} / make_float2 (_S5620.z);
    float u_105 = _S5646.x;
    float v_105 = _S5646.y;
    float r2_100 = u_105 * u_105 + v_105 * v_105;
    float2  _S5647 = _S5646 * make_float2 (1.0f + r2_100 * ((*dist_coeffs_45)[int(0)] + r2_100 * ((*dist_coeffs_45)[int(1)] + r2_100 * ((*dist_coeffs_45)[int(2)] + r2_100 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_105 * v_105 + (*dist_coeffs_45)[int(5)] * (r2_100 + 2.0f * u_105 * u_105) + (*dist_coeffs_45)[int(6)] * r2_100, _S5626 * u_105 * v_105 + (*dist_coeffs_45)[int(4)] * (r2_100 + 2.0f * v_105 * v_105) + (*dist_coeffs_45)[int(7)] * r2_100);
    float2  _S5648 = _S5647 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5647.x + (*dist_coeffs_45)[int(9)] * _S5647.y, 0.0f);
    float _S5649 = fx_38 * _S5648.x + cx_33;
    float _S5650 = fy_38 * _S5648.y + cy_33;
    uv_12[int(4)] = make_float2 (_S5649, _S5650);
    float2  _S5651 = float2 {_S5621.x, _S5621.y} / make_float2 (_S5621.z);
    float u_106 = _S5651.x;
    float v_106 = _S5651.y;
    float r2_101 = u_106 * u_106 + v_106 * v_106;
    float2  _S5652 = _S5651 * make_float2 (1.0f + r2_101 * ((*dist_coeffs_45)[int(0)] + r2_101 * ((*dist_coeffs_45)[int(1)] + r2_101 * ((*dist_coeffs_45)[int(2)] + r2_101 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_106 * v_106 + (*dist_coeffs_45)[int(5)] * (r2_101 + 2.0f * u_106 * u_106) + (*dist_coeffs_45)[int(6)] * r2_101, _S5626 * u_106 * v_106 + (*dist_coeffs_45)[int(4)] * (r2_101 + 2.0f * v_106 * v_106) + (*dist_coeffs_45)[int(7)] * r2_101);
    float2  _S5653 = _S5652 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5652.x + (*dist_coeffs_45)[int(9)] * _S5652.y, 0.0f);
    float _S5654 = fx_38 * _S5653.x + cx_33;
    float _S5655 = fy_38 * _S5653.y + cy_33;
    uv_12[int(5)] = make_float2 (_S5654, _S5655);
    float2  _S5656 = float2 {_S5622.x, _S5622.y} / make_float2 (_S5622.z);
    float u_107 = _S5656.x;
    float v_107 = _S5656.y;
    float r2_102 = u_107 * u_107 + v_107 * v_107;
    float2  _S5657 = _S5656 * make_float2 (1.0f + r2_102 * ((*dist_coeffs_45)[int(0)] + r2_102 * ((*dist_coeffs_45)[int(1)] + r2_102 * ((*dist_coeffs_45)[int(2)] + r2_102 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_107 * v_107 + (*dist_coeffs_45)[int(5)] * (r2_102 + 2.0f * u_107 * u_107) + (*dist_coeffs_45)[int(6)] * r2_102, _S5626 * u_107 * v_107 + (*dist_coeffs_45)[int(4)] * (r2_102 + 2.0f * v_107 * v_107) + (*dist_coeffs_45)[int(7)] * r2_102);
    float2  _S5658 = _S5657 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5657.x + (*dist_coeffs_45)[int(9)] * _S5657.y, 0.0f);
    float _S5659 = fx_38 * _S5658.x + cx_33;
    float _S5660 = fy_38 * _S5658.y + cy_33;
    uv_12[int(6)] = make_float2 (_S5659, _S5660);
    float2  _S5661 = float2 {_S5623.x, _S5623.y} / make_float2 (_S5623.z);
    float u_108 = _S5661.x;
    float v_108 = _S5661.y;
    float r2_103 = u_108 * u_108 + v_108 * v_108;
    float2  _S5662 = _S5661 * make_float2 (1.0f + r2_103 * ((*dist_coeffs_45)[int(0)] + r2_103 * ((*dist_coeffs_45)[int(1)] + r2_103 * ((*dist_coeffs_45)[int(2)] + r2_103 * (*dist_coeffs_45)[int(3)])))) + make_float2 (_S5625 * u_108 * v_108 + (*dist_coeffs_45)[int(5)] * (r2_103 + 2.0f * u_108 * u_108) + (*dist_coeffs_45)[int(6)] * r2_103, _S5626 * u_108 * v_108 + (*dist_coeffs_45)[int(4)] * (r2_103 + 2.0f * v_108 * v_108) + (*dist_coeffs_45)[int(7)] * r2_103);
    float2  _S5663 = _S5662 + make_float2 ((*dist_coeffs_45)[int(8)] * _S5662.x + (*dist_coeffs_45)[int(9)] * _S5662.y, 0.0f);
    float _S5664 = fx_38 * _S5663.x + cx_33;
    float _S5665 = fy_38 * _S5663.y + cy_33;
    uv_12[int(7)] = make_float2 (_S5664, _S5665);
    *aabb_xyxy_20 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5629), (_S5634)))), (_S5639)))), (_S5644)))), (_S5649)))), (_S5654)))), (_S5659)))), (_S5664))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5630), (_S5635)))), (_S5640)))), (_S5645)))), (_S5650)))), (_S5655)))), (_S5660)))), (_S5665))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5629), (_S5634)))), (_S5639)))), (_S5644)))), (_S5649)))), (_S5654)))), (_S5659)))), (_S5664))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5630), (_S5635)))), (_S5640)))), (_S5645)))), (_S5650)))), (_S5655)))), (_S5660)))), (_S5665))))))));
    *depth_23 = mean_c_29.z;
    float3  _S5666 = mean_c_29 - - mul_0(transpose_0(R_33), t_35);
    float3  _S5667 = make_float3 (0.282094806432724f) * (*sh_coeffs_29)[int(0)];
    *rgbs_9 = _S5667;
    float _S5668 = _S5666.x;
    float _S5669 = _S5666.y;
    float _S5670 = _S5666.z;
    float norm_20 = (F32_sqrt((_S5668 * _S5668 + _S5669 * _S5669 + _S5670 * _S5670)));
    float x_65 = _S5668 / norm_20;
    float y_32 = _S5669 / norm_20;
    float z_29 = _S5670 / norm_20;
    float3  _S5671 = _S5667 + make_float3 (0.48860251903533936f) * (make_float3 (- y_32) * (*sh_coeffs_29)[int(1)] + make_float3 (z_29) * (*sh_coeffs_29)[int(2)] - make_float3 (x_65) * (*sh_coeffs_29)[int(3)]);
    *rgbs_9 = _S5671;
    float z2_62 = z_29 * z_29;
    float fTmp0B_29 = -1.09254848957061768f * z_29;
    float fC1_29 = x_65 * x_65 - y_32 * y_32;
    float fS1_29 = 2.0f * x_65 * y_32;
    float3  _S5672 = _S5671 + (make_float3 (0.54627424478530884f * fS1_29) * (*sh_coeffs_29)[int(4)] + make_float3 (fTmp0B_29 * y_32) * (*sh_coeffs_29)[int(5)] + make_float3 (0.94617468118667603f * z2_62 - 0.31539157032966614f) * (*sh_coeffs_29)[int(6)] + make_float3 (fTmp0B_29 * x_65) * (*sh_coeffs_29)[int(7)] + make_float3 (0.54627424478530884f * fC1_29) * (*sh_coeffs_29)[int(8)]);
    *rgbs_9 = _S5672;
    float fTmp0C_29 = -2.28522896766662598f * z2_62 + 0.4570457935333252f;
    float fTmp1B_29 = 1.44530570507049561f * z_29;
    *rgbs_9 = max_0(_S5672 + (make_float3 (-0.59004360437393188f * (x_65 * fS1_29 + y_32 * fC1_29)) * (*sh_coeffs_29)[int(9)] + make_float3 (fTmp1B_29 * fS1_29) * (*sh_coeffs_29)[int(10)] + make_float3 (fTmp0C_29 * y_32) * (*sh_coeffs_29)[int(11)] + make_float3 (z_29 * (1.86588168144226074f * z2_62 - 1.11952900886535645f)) * (*sh_coeffs_29)[int(12)] + make_float3 (fTmp0C_29 * x_65) * (*sh_coeffs_29)[int(13)] + make_float3 (fTmp1B_29 * fC1_29) * (*sh_coeffs_29)[int(14)] + make_float3 (-0.59004360437393188f * (x_65 * fC1_29 - y_32 * fS1_29)) * (*sh_coeffs_29)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_voxel_eval3d_fisheye_differentiable(float3  pos_3, float size_3, FixedArray<float, 8>  * densities_3, FixedArray<float3 , 16>  * sh_coeffs_30, Matrix<float, 3, 3>  R_34, float3  t_36, float fx_39, float fy_39, float cx_34, float cy_34, FixedArray<float, 10>  * dist_coeffs_46, uint image_width_30, uint image_height_30, float near_plane_21, float far_plane_21, int4  * aabb_xyxy_21, float * depth_24, float3  * rgbs_10)
{
    FixedArray<float3 , 8>  pos_c_3;
    float3  _S5673 = mul_0(R_34, pos_3) + t_36;
    pos_c_3[int(0)] = _S5673;
    pos_c_3[int(1)] = mul_0(R_34, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 0.0f)) + t_36;
    pos_c_3[int(2)] = mul_0(R_34, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 0.0f)) + t_36;
    pos_c_3[int(3)] = mul_0(R_34, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 1.0f, 0.0f)) + t_36;
    pos_c_3[int(4)] = mul_0(R_34, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 0.0f, 1.0f)) + t_36;
    pos_c_3[int(5)] = mul_0(R_34, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 1.0f)) + t_36;
    pos_c_3[int(6)] = mul_0(R_34, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 1.0f)) + t_36;
    pos_c_3[int(7)] = mul_0(R_34, pos_3 + make_float3 (size_3)) + t_36;
    float3  mean_c_30 = mul_0(R_34, pos_3 + make_float3 (0.5f * size_3)) + t_36;
    FixedArray<float2 , 8>  uv_13;
    float2  _S5674 = float2 {_S5673.x, _S5673.y};
    float r_44 = length_0(_S5674);
    float _S5675 = _S5673.z;
    float theta_40 = (F32_atan2((r_44), (_S5675)));
    float k_14;
    if(theta_40 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_40 * theta_40 / 3.0f) / _S5675;
    }
    else
    {
        k_14 = theta_40 / r_44;
    }
    float2  _S5676 = _S5674 * make_float2 (k_14);
    float u_109 = _S5676.x;
    float v_109 = _S5676.y;
    float r2_104 = u_109 * u_109 + v_109 * v_109;
    float _S5677 = 2.0f * (*dist_coeffs_46)[int(4)];
    float _S5678 = 2.0f * (*dist_coeffs_46)[int(5)];
    float2  _S5679 = _S5676 * make_float2 (1.0f + r2_104 * ((*dist_coeffs_46)[int(0)] + r2_104 * ((*dist_coeffs_46)[int(1)] + r2_104 * ((*dist_coeffs_46)[int(2)] + r2_104 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_109 * v_109 + (*dist_coeffs_46)[int(5)] * (r2_104 + 2.0f * u_109 * u_109) + (*dist_coeffs_46)[int(6)] * r2_104, _S5678 * u_109 * v_109 + (*dist_coeffs_46)[int(4)] * (r2_104 + 2.0f * v_109 * v_109) + (*dist_coeffs_46)[int(7)] * r2_104);
    float2  _S5680 = _S5679 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5679.x + (*dist_coeffs_46)[int(9)] * _S5679.y, 0.0f);
    uv_13[int(0)] = make_float2 (fx_39 * _S5680.x + cx_34, fy_39 * _S5680.y + cy_34);
    float2  _S5681 = float2 {pos_c_3[int(1)].x, pos_c_3[int(1)].y};
    float r_45 = length_0(_S5681);
    float _S5682 = pos_c_3[int(1)].z;
    float theta_41 = (F32_atan2((r_45), (_S5682)));
    if(theta_41 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_41 * theta_41 / 3.0f) / _S5682;
    }
    else
    {
        k_14 = theta_41 / r_45;
    }
    float2  _S5683 = _S5681 * make_float2 (k_14);
    float u_110 = _S5683.x;
    float v_110 = _S5683.y;
    float r2_105 = u_110 * u_110 + v_110 * v_110;
    float2  _S5684 = _S5683 * make_float2 (1.0f + r2_105 * ((*dist_coeffs_46)[int(0)] + r2_105 * ((*dist_coeffs_46)[int(1)] + r2_105 * ((*dist_coeffs_46)[int(2)] + r2_105 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_110 * v_110 + (*dist_coeffs_46)[int(5)] * (r2_105 + 2.0f * u_110 * u_110) + (*dist_coeffs_46)[int(6)] * r2_105, _S5678 * u_110 * v_110 + (*dist_coeffs_46)[int(4)] * (r2_105 + 2.0f * v_110 * v_110) + (*dist_coeffs_46)[int(7)] * r2_105);
    float2  _S5685 = _S5684 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5684.x + (*dist_coeffs_46)[int(9)] * _S5684.y, 0.0f);
    uv_13[int(1)] = make_float2 (fx_39 * _S5685.x + cx_34, fy_39 * _S5685.y + cy_34);
    float2  _S5686 = float2 {pos_c_3[int(2)].x, pos_c_3[int(2)].y};
    float r_46 = length_0(_S5686);
    float _S5687 = pos_c_3[int(2)].z;
    float theta_42 = (F32_atan2((r_46), (_S5687)));
    if(theta_42 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_42 * theta_42 / 3.0f) / _S5687;
    }
    else
    {
        k_14 = theta_42 / r_46;
    }
    float2  _S5688 = _S5686 * make_float2 (k_14);
    float u_111 = _S5688.x;
    float v_111 = _S5688.y;
    float r2_106 = u_111 * u_111 + v_111 * v_111;
    float2  _S5689 = _S5688 * make_float2 (1.0f + r2_106 * ((*dist_coeffs_46)[int(0)] + r2_106 * ((*dist_coeffs_46)[int(1)] + r2_106 * ((*dist_coeffs_46)[int(2)] + r2_106 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_111 * v_111 + (*dist_coeffs_46)[int(5)] * (r2_106 + 2.0f * u_111 * u_111) + (*dist_coeffs_46)[int(6)] * r2_106, _S5678 * u_111 * v_111 + (*dist_coeffs_46)[int(4)] * (r2_106 + 2.0f * v_111 * v_111) + (*dist_coeffs_46)[int(7)] * r2_106);
    float2  _S5690 = _S5689 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5689.x + (*dist_coeffs_46)[int(9)] * _S5689.y, 0.0f);
    uv_13[int(2)] = make_float2 (fx_39 * _S5690.x + cx_34, fy_39 * _S5690.y + cy_34);
    float2  _S5691 = float2 {pos_c_3[int(3)].x, pos_c_3[int(3)].y};
    float r_47 = length_0(_S5691);
    float _S5692 = pos_c_3[int(3)].z;
    float theta_43 = (F32_atan2((r_47), (_S5692)));
    if(theta_43 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_43 * theta_43 / 3.0f) / _S5692;
    }
    else
    {
        k_14 = theta_43 / r_47;
    }
    float2  _S5693 = _S5691 * make_float2 (k_14);
    float u_112 = _S5693.x;
    float v_112 = _S5693.y;
    float r2_107 = u_112 * u_112 + v_112 * v_112;
    float2  _S5694 = _S5693 * make_float2 (1.0f + r2_107 * ((*dist_coeffs_46)[int(0)] + r2_107 * ((*dist_coeffs_46)[int(1)] + r2_107 * ((*dist_coeffs_46)[int(2)] + r2_107 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_112 * v_112 + (*dist_coeffs_46)[int(5)] * (r2_107 + 2.0f * u_112 * u_112) + (*dist_coeffs_46)[int(6)] * r2_107, _S5678 * u_112 * v_112 + (*dist_coeffs_46)[int(4)] * (r2_107 + 2.0f * v_112 * v_112) + (*dist_coeffs_46)[int(7)] * r2_107);
    float2  _S5695 = _S5694 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5694.x + (*dist_coeffs_46)[int(9)] * _S5694.y, 0.0f);
    uv_13[int(3)] = make_float2 (fx_39 * _S5695.x + cx_34, fy_39 * _S5695.y + cy_34);
    float2  _S5696 = float2 {pos_c_3[int(4)].x, pos_c_3[int(4)].y};
    float r_48 = length_0(_S5696);
    float _S5697 = pos_c_3[int(4)].z;
    float theta_44 = (F32_atan2((r_48), (_S5697)));
    if(theta_44 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_44 * theta_44 / 3.0f) / _S5697;
    }
    else
    {
        k_14 = theta_44 / r_48;
    }
    float2  _S5698 = _S5696 * make_float2 (k_14);
    float u_113 = _S5698.x;
    float v_113 = _S5698.y;
    float r2_108 = u_113 * u_113 + v_113 * v_113;
    float2  _S5699 = _S5698 * make_float2 (1.0f + r2_108 * ((*dist_coeffs_46)[int(0)] + r2_108 * ((*dist_coeffs_46)[int(1)] + r2_108 * ((*dist_coeffs_46)[int(2)] + r2_108 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_113 * v_113 + (*dist_coeffs_46)[int(5)] * (r2_108 + 2.0f * u_113 * u_113) + (*dist_coeffs_46)[int(6)] * r2_108, _S5678 * u_113 * v_113 + (*dist_coeffs_46)[int(4)] * (r2_108 + 2.0f * v_113 * v_113) + (*dist_coeffs_46)[int(7)] * r2_108);
    float2  _S5700 = _S5699 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5699.x + (*dist_coeffs_46)[int(9)] * _S5699.y, 0.0f);
    uv_13[int(4)] = make_float2 (fx_39 * _S5700.x + cx_34, fy_39 * _S5700.y + cy_34);
    float2  _S5701 = float2 {pos_c_3[int(5)].x, pos_c_3[int(5)].y};
    float r_49 = length_0(_S5701);
    float _S5702 = pos_c_3[int(5)].z;
    float theta_45 = (F32_atan2((r_49), (_S5702)));
    if(theta_45 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_45 * theta_45 / 3.0f) / _S5702;
    }
    else
    {
        k_14 = theta_45 / r_49;
    }
    float2  _S5703 = _S5701 * make_float2 (k_14);
    float u_114 = _S5703.x;
    float v_114 = _S5703.y;
    float r2_109 = u_114 * u_114 + v_114 * v_114;
    float2  _S5704 = _S5703 * make_float2 (1.0f + r2_109 * ((*dist_coeffs_46)[int(0)] + r2_109 * ((*dist_coeffs_46)[int(1)] + r2_109 * ((*dist_coeffs_46)[int(2)] + r2_109 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_114 * v_114 + (*dist_coeffs_46)[int(5)] * (r2_109 + 2.0f * u_114 * u_114) + (*dist_coeffs_46)[int(6)] * r2_109, _S5678 * u_114 * v_114 + (*dist_coeffs_46)[int(4)] * (r2_109 + 2.0f * v_114 * v_114) + (*dist_coeffs_46)[int(7)] * r2_109);
    float2  _S5705 = _S5704 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5704.x + (*dist_coeffs_46)[int(9)] * _S5704.y, 0.0f);
    uv_13[int(5)] = make_float2 (fx_39 * _S5705.x + cx_34, fy_39 * _S5705.y + cy_34);
    float2  _S5706 = float2 {pos_c_3[int(6)].x, pos_c_3[int(6)].y};
    float r_50 = length_0(_S5706);
    float _S5707 = pos_c_3[int(6)].z;
    float theta_46 = (F32_atan2((r_50), (_S5707)));
    if(theta_46 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_46 * theta_46 / 3.0f) / _S5707;
    }
    else
    {
        k_14 = theta_46 / r_50;
    }
    float2  _S5708 = _S5706 * make_float2 (k_14);
    float u_115 = _S5708.x;
    float v_115 = _S5708.y;
    float r2_110 = u_115 * u_115 + v_115 * v_115;
    float2  _S5709 = _S5708 * make_float2 (1.0f + r2_110 * ((*dist_coeffs_46)[int(0)] + r2_110 * ((*dist_coeffs_46)[int(1)] + r2_110 * ((*dist_coeffs_46)[int(2)] + r2_110 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_115 * v_115 + (*dist_coeffs_46)[int(5)] * (r2_110 + 2.0f * u_115 * u_115) + (*dist_coeffs_46)[int(6)] * r2_110, _S5678 * u_115 * v_115 + (*dist_coeffs_46)[int(4)] * (r2_110 + 2.0f * v_115 * v_115) + (*dist_coeffs_46)[int(7)] * r2_110);
    float2  _S5710 = _S5709 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5709.x + (*dist_coeffs_46)[int(9)] * _S5709.y, 0.0f);
    uv_13[int(6)] = make_float2 (fx_39 * _S5710.x + cx_34, fy_39 * _S5710.y + cy_34);
    float2  _S5711 = float2 {pos_c_3[int(7)].x, pos_c_3[int(7)].y};
    float r_51 = length_0(_S5711);
    float _S5712 = pos_c_3[int(7)].z;
    float theta_47 = (F32_atan2((r_51), (_S5712)));
    if(theta_47 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_47 * theta_47 / 3.0f) / _S5712;
    }
    else
    {
        k_14 = theta_47 / r_51;
    }
    float2  _S5713 = _S5711 * make_float2 (k_14);
    float u_116 = _S5713.x;
    float v_116 = _S5713.y;
    float r2_111 = u_116 * u_116 + v_116 * v_116;
    float2  _S5714 = _S5713 * make_float2 (1.0f + r2_111 * ((*dist_coeffs_46)[int(0)] + r2_111 * ((*dist_coeffs_46)[int(1)] + r2_111 * ((*dist_coeffs_46)[int(2)] + r2_111 * (*dist_coeffs_46)[int(3)])))) + make_float2 (_S5677 * u_116 * v_116 + (*dist_coeffs_46)[int(5)] * (r2_111 + 2.0f * u_116 * u_116) + (*dist_coeffs_46)[int(6)] * r2_111, _S5678 * u_116 * v_116 + (*dist_coeffs_46)[int(4)] * (r2_111 + 2.0f * v_116 * v_116) + (*dist_coeffs_46)[int(7)] * r2_111);
    float2  _S5715 = _S5714 + make_float2 ((*dist_coeffs_46)[int(8)] * _S5714.x + (*dist_coeffs_46)[int(9)] * _S5714.y, 0.0f);
    float _S5716 = fx_39 * _S5715.x + cx_34;
    float _S5717 = fy_39 * _S5715.y + cy_34;
    uv_13[int(7)] = make_float2 (_S5716, _S5717);
    *aabb_xyxy_21 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_13[int(0)].x), (uv_13[int(1)].x)))), (uv_13[int(2)].x)))), (uv_13[int(3)].x)))), (uv_13[int(4)].x)))), (uv_13[int(5)].x)))), (uv_13[int(6)].x)))), (_S5716))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_13[int(0)].y), (uv_13[int(1)].y)))), (uv_13[int(2)].y)))), (uv_13[int(3)].y)))), (uv_13[int(4)].y)))), (uv_13[int(5)].y)))), (uv_13[int(6)].y)))), (_S5717))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_13[int(0)].x), (uv_13[int(1)].x)))), (uv_13[int(2)].x)))), (uv_13[int(3)].x)))), (uv_13[int(4)].x)))), (uv_13[int(5)].x)))), (uv_13[int(6)].x)))), (_S5716))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_13[int(0)].y), (uv_13[int(1)].y)))), (uv_13[int(2)].y)))), (uv_13[int(3)].y)))), (uv_13[int(4)].y)))), (uv_13[int(5)].y)))), (uv_13[int(6)].y)))), (_S5717))))))));
    *depth_24 = mean_c_30.z;
    float3  _S5718 = mean_c_30 - - mul_0(transpose_0(R_34), t_36);
    float3  _S5719 = make_float3 (0.282094806432724f) * (*sh_coeffs_30)[int(0)];
    *rgbs_10 = _S5719;
    float _S5720 = _S5718.x;
    float _S5721 = _S5718.y;
    float _S5722 = _S5718.z;
    float norm_21 = (F32_sqrt((_S5720 * _S5720 + _S5721 * _S5721 + _S5722 * _S5722)));
    float x_66 = _S5720 / norm_21;
    float y_33 = _S5721 / norm_21;
    float z_30 = _S5722 / norm_21;
    float3  _S5723 = _S5719 + make_float3 (0.48860251903533936f) * (make_float3 (- y_33) * (*sh_coeffs_30)[int(1)] + make_float3 (z_30) * (*sh_coeffs_30)[int(2)] - make_float3 (x_66) * (*sh_coeffs_30)[int(3)]);
    *rgbs_10 = _S5723;
    float z2_63 = z_30 * z_30;
    float fTmp0B_30 = -1.09254848957061768f * z_30;
    float fC1_30 = x_66 * x_66 - y_33 * y_33;
    float fS1_30 = 2.0f * x_66 * y_33;
    float3  _S5724 = _S5723 + (make_float3 (0.54627424478530884f * fS1_30) * (*sh_coeffs_30)[int(4)] + make_float3 (fTmp0B_30 * y_33) * (*sh_coeffs_30)[int(5)] + make_float3 (0.94617468118667603f * z2_63 - 0.31539157032966614f) * (*sh_coeffs_30)[int(6)] + make_float3 (fTmp0B_30 * x_66) * (*sh_coeffs_30)[int(7)] + make_float3 (0.54627424478530884f * fC1_30) * (*sh_coeffs_30)[int(8)]);
    *rgbs_10 = _S5724;
    float fTmp0C_30 = -2.28522896766662598f * z2_63 + 0.4570457935333252f;
    float fTmp1B_30 = 1.44530570507049561f * z_30;
    *rgbs_10 = max_0(_S5724 + (make_float3 (-0.59004360437393188f * (x_66 * fS1_30 + y_33 * fC1_30)) * (*sh_coeffs_30)[int(9)] + make_float3 (fTmp1B_30 * fS1_30) * (*sh_coeffs_30)[int(10)] + make_float3 (fTmp0C_30 * y_33) * (*sh_coeffs_30)[int(11)] + make_float3 (z_30 * (1.86588168144226074f * z2_63 - 1.11952900886535645f)) * (*sh_coeffs_30)[int(12)] + make_float3 (fTmp0C_30 * x_66) * (*sh_coeffs_30)[int(13)] + make_float3 (fTmp1B_30 * fC1_30) * (*sh_coeffs_30)[int(14)] + make_float3 (-0.59004360437393188f * (x_66 * fC1_30 - y_33 * fS1_30)) * (*sh_coeffs_30)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void projection_voxel_eval3d_persp_vjp(float3  pos_4, float size_4, FixedArray<float, 8>  * densities_4, FixedArray<float3 , 16>  * sh_coeffs_31, Matrix<float, 3, 3>  R_35, float3  t_37, float fx_40, float fy_40, float cx_35, float cy_35, FixedArray<float, 10>  * dist_coeffs_47, uint image_width_31, uint image_height_31, float3  v_rgb_8, FixedArray<float, 8>  * v_densities_0, FixedArray<float3 , 16>  * v_sh_coeffs_9, Matrix<float, 3, 3>  * v_R_10, float3  * v_t_10)
{
    float3  _S5725 = s_primal_ctx_mul_1(R_35, pos_4) + t_37;
    float _S5726 = _S5725.z;
    float2  _S5727 = make_float2 (_S5726);
    float _S5728 = s_primal_ctx_min_0(1.00000001504746622e+30f, _S5726);
    float _S5729 = s_primal_ctx_max_0(0.0f, _S5726);
    float3  pos_i_0 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S5730 = s_primal_ctx_mul_1(R_35, pos_i_0) + t_37;
    float _S5731 = _S5730.z;
    float2  _S5732 = make_float2 (_S5731);
    float _S5733 = s_primal_ctx_min_0(_S5728, _S5731);
    float _S5734 = s_primal_ctx_max_0(_S5729, _S5731);
    float3  pos_i_1 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S5735 = s_primal_ctx_mul_1(R_35, pos_i_1) + t_37;
    float _S5736 = _S5735.z;
    float2  _S5737 = make_float2 (_S5736);
    float _S5738 = s_primal_ctx_min_0(_S5733, _S5736);
    float _S5739 = s_primal_ctx_max_0(_S5734, _S5736);
    float3  pos_i_2 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S5740 = s_primal_ctx_mul_1(R_35, pos_i_2) + t_37;
    float _S5741 = _S5740.z;
    float2  _S5742 = make_float2 (_S5741);
    float _S5743 = s_primal_ctx_min_0(_S5738, _S5741);
    float _S5744 = s_primal_ctx_max_0(_S5739, _S5741);
    float3  pos_i_3 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S5745 = s_primal_ctx_mul_1(R_35, pos_i_3) + t_37;
    float _S5746 = _S5745.z;
    float2  _S5747 = make_float2 (_S5746);
    float _S5748 = s_primal_ctx_min_0(_S5743, _S5746);
    float _S5749 = s_primal_ctx_max_0(_S5744, _S5746);
    float3  pos_i_4 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S5750 = s_primal_ctx_mul_1(R_35, pos_i_4) + t_37;
    float _S5751 = _S5750.z;
    float2  _S5752 = make_float2 (_S5751);
    float _S5753 = s_primal_ctx_min_0(_S5748, _S5751);
    float _S5754 = s_primal_ctx_max_0(_S5749, _S5751);
    float3  pos_i_5 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S5755 = s_primal_ctx_mul_1(R_35, pos_i_5) + t_37;
    float _S5756 = _S5755.z;
    float2  _S5757 = make_float2 (_S5756);
    float _S5758 = s_primal_ctx_min_0(_S5753, _S5756);
    float _S5759 = s_primal_ctx_max_0(_S5754, _S5756);
    float3  pos_i_6 = pos_4 + make_float3 (size_4);
    float3  _S5760 = s_primal_ctx_mul_1(R_35, pos_i_6) + t_37;
    float _S5761 = _S5760.z;
    float2  _S5762 = make_float2 (_S5761);
    float3  _S5763 = pos_4 + make_float3 (0.5f * size_4);
    float2  _S5764 = float2 {_S5725.x, _S5725.y};
    float2  _S5765 = _S5764 / make_float2 (_S5726);
    float2  _S5766 = make_float2 (_S5726 * _S5726);
    float u_117 = _S5765.x;
    float v_117 = _S5765.y;
    float r2_112 = u_117 * u_117 + v_117 * v_117;
    float _S5767 = (*dist_coeffs_47)[int(2)] + r2_112 * (*dist_coeffs_47)[int(3)];
    float _S5768 = (*dist_coeffs_47)[int(1)] + r2_112 * _S5767;
    float _S5769 = (*dist_coeffs_47)[int(0)] + r2_112 * _S5768;
    float radial_13 = 1.0f + r2_112 * _S5769;
    float2  _S5770 = make_float2 (radial_13);
    float _S5771 = 2.0f * (*dist_coeffs_47)[int(4)];
    float _S5772 = _S5771 * u_117;
    float _S5773 = 2.0f * u_117;
    float _S5774 = 2.0f * (*dist_coeffs_47)[int(5)];
    float _S5775 = _S5774 * u_117;
    float _S5776 = 2.0f * v_117;
    float2  _S5777 = _S5765 * make_float2 (radial_13) + make_float2 (_S5772 * v_117 + (*dist_coeffs_47)[int(5)] * (r2_112 + _S5773 * u_117) + (*dist_coeffs_47)[int(6)] * r2_112, _S5775 * v_117 + (*dist_coeffs_47)[int(4)] * (r2_112 + _S5776 * v_117) + (*dist_coeffs_47)[int(7)] * r2_112);
    float2  _S5778 = _S5777 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5777.x + (*dist_coeffs_47)[int(9)] * _S5777.y, 0.0f);
    float _S5779 = fx_40 * _S5778.x + cx_35;
    float _S5780 = fy_40 * _S5778.y + cy_35;
    float2  _S5781 = float2 {_S5730.x, _S5730.y};
    float2  _S5782 = _S5781 / make_float2 (_S5731);
    float2  _S5783 = make_float2 (_S5731 * _S5731);
    float u_118 = _S5782.x;
    float v_118 = _S5782.y;
    float r2_113 = u_118 * u_118 + v_118 * v_118;
    float _S5784 = (*dist_coeffs_47)[int(2)] + r2_113 * (*dist_coeffs_47)[int(3)];
    float _S5785 = (*dist_coeffs_47)[int(1)] + r2_113 * _S5784;
    float _S5786 = (*dist_coeffs_47)[int(0)] + r2_113 * _S5785;
    float radial_14 = 1.0f + r2_113 * _S5786;
    float2  _S5787 = make_float2 (radial_14);
    float _S5788 = _S5771 * u_118;
    float _S5789 = 2.0f * u_118;
    float _S5790 = _S5774 * u_118;
    float _S5791 = 2.0f * v_118;
    float2  _S5792 = _S5782 * make_float2 (radial_14) + make_float2 (_S5788 * v_118 + (*dist_coeffs_47)[int(5)] * (r2_113 + _S5789 * u_118) + (*dist_coeffs_47)[int(6)] * r2_113, _S5790 * v_118 + (*dist_coeffs_47)[int(4)] * (r2_113 + _S5791 * v_118) + (*dist_coeffs_47)[int(7)] * r2_113);
    float2  _S5793 = _S5792 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5792.x + (*dist_coeffs_47)[int(9)] * _S5792.y, 0.0f);
    float _S5794 = fx_40 * _S5793.x + cx_35;
    float _S5795 = fy_40 * _S5793.y + cy_35;
    float2  _S5796 = float2 {_S5735.x, _S5735.y};
    float2  _S5797 = _S5796 / make_float2 (_S5736);
    float2  _S5798 = make_float2 (_S5736 * _S5736);
    float u_119 = _S5797.x;
    float v_119 = _S5797.y;
    float r2_114 = u_119 * u_119 + v_119 * v_119;
    float _S5799 = (*dist_coeffs_47)[int(2)] + r2_114 * (*dist_coeffs_47)[int(3)];
    float _S5800 = (*dist_coeffs_47)[int(1)] + r2_114 * _S5799;
    float _S5801 = (*dist_coeffs_47)[int(0)] + r2_114 * _S5800;
    float radial_15 = 1.0f + r2_114 * _S5801;
    float2  _S5802 = make_float2 (radial_15);
    float _S5803 = _S5771 * u_119;
    float _S5804 = 2.0f * u_119;
    float _S5805 = _S5774 * u_119;
    float _S5806 = 2.0f * v_119;
    float2  _S5807 = _S5797 * make_float2 (radial_15) + make_float2 (_S5803 * v_119 + (*dist_coeffs_47)[int(5)] * (r2_114 + _S5804 * u_119) + (*dist_coeffs_47)[int(6)] * r2_114, _S5805 * v_119 + (*dist_coeffs_47)[int(4)] * (r2_114 + _S5806 * v_119) + (*dist_coeffs_47)[int(7)] * r2_114);
    float2  _S5808 = _S5807 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5807.x + (*dist_coeffs_47)[int(9)] * _S5807.y, 0.0f);
    float _S5809 = fx_40 * _S5808.x + cx_35;
    float _S5810 = fy_40 * _S5808.y + cy_35;
    float2  _S5811 = float2 {_S5740.x, _S5740.y};
    float2  _S5812 = _S5811 / make_float2 (_S5741);
    float2  _S5813 = make_float2 (_S5741 * _S5741);
    float u_120 = _S5812.x;
    float v_120 = _S5812.y;
    float r2_115 = u_120 * u_120 + v_120 * v_120;
    float _S5814 = (*dist_coeffs_47)[int(2)] + r2_115 * (*dist_coeffs_47)[int(3)];
    float _S5815 = (*dist_coeffs_47)[int(1)] + r2_115 * _S5814;
    float _S5816 = (*dist_coeffs_47)[int(0)] + r2_115 * _S5815;
    float radial_16 = 1.0f + r2_115 * _S5816;
    float2  _S5817 = make_float2 (radial_16);
    float _S5818 = _S5771 * u_120;
    float _S5819 = 2.0f * u_120;
    float _S5820 = _S5774 * u_120;
    float _S5821 = 2.0f * v_120;
    float2  _S5822 = _S5812 * make_float2 (radial_16) + make_float2 (_S5818 * v_120 + (*dist_coeffs_47)[int(5)] * (r2_115 + _S5819 * u_120) + (*dist_coeffs_47)[int(6)] * r2_115, _S5820 * v_120 + (*dist_coeffs_47)[int(4)] * (r2_115 + _S5821 * v_120) + (*dist_coeffs_47)[int(7)] * r2_115);
    float2  _S5823 = _S5822 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5822.x + (*dist_coeffs_47)[int(9)] * _S5822.y, 0.0f);
    float _S5824 = fx_40 * _S5823.x + cx_35;
    float _S5825 = fy_40 * _S5823.y + cy_35;
    float2  _S5826 = float2 {_S5745.x, _S5745.y};
    float2  _S5827 = _S5826 / make_float2 (_S5746);
    float2  _S5828 = make_float2 (_S5746 * _S5746);
    float u_121 = _S5827.x;
    float v_121 = _S5827.y;
    float r2_116 = u_121 * u_121 + v_121 * v_121;
    float _S5829 = (*dist_coeffs_47)[int(2)] + r2_116 * (*dist_coeffs_47)[int(3)];
    float _S5830 = (*dist_coeffs_47)[int(1)] + r2_116 * _S5829;
    float _S5831 = (*dist_coeffs_47)[int(0)] + r2_116 * _S5830;
    float radial_17 = 1.0f + r2_116 * _S5831;
    float2  _S5832 = make_float2 (radial_17);
    float _S5833 = _S5771 * u_121;
    float _S5834 = 2.0f * u_121;
    float _S5835 = _S5774 * u_121;
    float _S5836 = 2.0f * v_121;
    float2  _S5837 = _S5827 * make_float2 (radial_17) + make_float2 (_S5833 * v_121 + (*dist_coeffs_47)[int(5)] * (r2_116 + _S5834 * u_121) + (*dist_coeffs_47)[int(6)] * r2_116, _S5835 * v_121 + (*dist_coeffs_47)[int(4)] * (r2_116 + _S5836 * v_121) + (*dist_coeffs_47)[int(7)] * r2_116);
    float2  _S5838 = _S5837 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5837.x + (*dist_coeffs_47)[int(9)] * _S5837.y, 0.0f);
    float _S5839 = fx_40 * _S5838.x + cx_35;
    float _S5840 = fy_40 * _S5838.y + cy_35;
    float2  _S5841 = float2 {_S5750.x, _S5750.y};
    float2  _S5842 = _S5841 / make_float2 (_S5751);
    float2  _S5843 = make_float2 (_S5751 * _S5751);
    float u_122 = _S5842.x;
    float v_122 = _S5842.y;
    float r2_117 = u_122 * u_122 + v_122 * v_122;
    float _S5844 = (*dist_coeffs_47)[int(2)] + r2_117 * (*dist_coeffs_47)[int(3)];
    float _S5845 = (*dist_coeffs_47)[int(1)] + r2_117 * _S5844;
    float _S5846 = (*dist_coeffs_47)[int(0)] + r2_117 * _S5845;
    float radial_18 = 1.0f + r2_117 * _S5846;
    float2  _S5847 = make_float2 (radial_18);
    float _S5848 = _S5771 * u_122;
    float _S5849 = 2.0f * u_122;
    float _S5850 = _S5774 * u_122;
    float _S5851 = 2.0f * v_122;
    float2  _S5852 = _S5842 * make_float2 (radial_18) + make_float2 (_S5848 * v_122 + (*dist_coeffs_47)[int(5)] * (r2_117 + _S5849 * u_122) + (*dist_coeffs_47)[int(6)] * r2_117, _S5850 * v_122 + (*dist_coeffs_47)[int(4)] * (r2_117 + _S5851 * v_122) + (*dist_coeffs_47)[int(7)] * r2_117);
    float2  _S5853 = _S5852 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5852.x + (*dist_coeffs_47)[int(9)] * _S5852.y, 0.0f);
    float _S5854 = fx_40 * _S5853.x + cx_35;
    float _S5855 = fy_40 * _S5853.y + cy_35;
    float2  _S5856 = float2 {_S5755.x, _S5755.y};
    float2  _S5857 = _S5856 / make_float2 (_S5756);
    float2  _S5858 = make_float2 (_S5756 * _S5756);
    float u_123 = _S5857.x;
    float v_123 = _S5857.y;
    float r2_118 = u_123 * u_123 + v_123 * v_123;
    float _S5859 = (*dist_coeffs_47)[int(2)] + r2_118 * (*dist_coeffs_47)[int(3)];
    float _S5860 = (*dist_coeffs_47)[int(1)] + r2_118 * _S5859;
    float _S5861 = (*dist_coeffs_47)[int(0)] + r2_118 * _S5860;
    float radial_19 = 1.0f + r2_118 * _S5861;
    float2  _S5862 = make_float2 (radial_19);
    float _S5863 = _S5771 * u_123;
    float _S5864 = 2.0f * u_123;
    float _S5865 = _S5774 * u_123;
    float _S5866 = 2.0f * v_123;
    float2  _S5867 = _S5857 * make_float2 (radial_19) + make_float2 (_S5863 * v_123 + (*dist_coeffs_47)[int(5)] * (r2_118 + _S5864 * u_123) + (*dist_coeffs_47)[int(6)] * r2_118, _S5865 * v_123 + (*dist_coeffs_47)[int(4)] * (r2_118 + _S5866 * v_123) + (*dist_coeffs_47)[int(7)] * r2_118);
    float2  _S5868 = _S5867 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5867.x + (*dist_coeffs_47)[int(9)] * _S5867.y, 0.0f);
    float _S5869 = fx_40 * _S5868.x + cx_35;
    float _S5870 = fy_40 * _S5868.y + cy_35;
    float2  _S5871 = float2 {_S5760.x, _S5760.y};
    float2  _S5872 = _S5871 / make_float2 (_S5761);
    float2  _S5873 = make_float2 (_S5761 * _S5761);
    float u_124 = _S5872.x;
    float v_124 = _S5872.y;
    float r2_119 = u_124 * u_124 + v_124 * v_124;
    float _S5874 = (*dist_coeffs_47)[int(2)] + r2_119 * (*dist_coeffs_47)[int(3)];
    float _S5875 = (*dist_coeffs_47)[int(1)] + r2_119 * _S5874;
    float _S5876 = (*dist_coeffs_47)[int(0)] + r2_119 * _S5875;
    float radial_20 = 1.0f + r2_119 * _S5876;
    float2  _S5877 = make_float2 (radial_20);
    float _S5878 = _S5771 * u_124;
    float _S5879 = 2.0f * u_124;
    float _S5880 = _S5774 * u_124;
    float _S5881 = 2.0f * v_124;
    float2  _S5882 = _S5872 * make_float2 (radial_20) + make_float2 (_S5878 * v_124 + (*dist_coeffs_47)[int(5)] * (r2_119 + _S5879 * u_124) + (*dist_coeffs_47)[int(6)] * r2_119, _S5880 * v_124 + (*dist_coeffs_47)[int(4)] * (r2_119 + _S5881 * v_124) + (*dist_coeffs_47)[int(7)] * r2_119);
    float2  _S5883 = _S5882 + make_float2 ((*dist_coeffs_47)[int(8)] * _S5882.x + (*dist_coeffs_47)[int(9)] * _S5882.y, 0.0f);
    float _S5884 = fx_40 * _S5883.x + cx_35;
    float _S5885 = fy_40 * _S5883.y + cy_35;
    float _S5886 = s_primal_ctx_max_0(_S5779, _S5794);
    float _S5887 = s_primal_ctx_min_0(_S5779, _S5794);
    float _S5888 = s_primal_ctx_max_0(_S5780, _S5795);
    float _S5889 = s_primal_ctx_min_0(_S5780, _S5795);
    float _S5890 = s_primal_ctx_max_0(_S5886, _S5809);
    float _S5891 = s_primal_ctx_min_0(_S5887, _S5809);
    float _S5892 = s_primal_ctx_max_0(_S5888, _S5810);
    float _S5893 = s_primal_ctx_min_0(_S5889, _S5810);
    float _S5894 = s_primal_ctx_max_0(_S5890, _S5824);
    float _S5895 = s_primal_ctx_min_0(_S5891, _S5824);
    float _S5896 = s_primal_ctx_max_0(_S5892, _S5825);
    float _S5897 = s_primal_ctx_min_0(_S5893, _S5825);
    float _S5898 = s_primal_ctx_max_0(_S5894, _S5839);
    float _S5899 = s_primal_ctx_min_0(_S5895, _S5839);
    float _S5900 = s_primal_ctx_max_0(_S5896, _S5840);
    float _S5901 = s_primal_ctx_min_0(_S5897, _S5840);
    float _S5902 = s_primal_ctx_max_0(_S5898, _S5854);
    float _S5903 = s_primal_ctx_min_0(_S5899, _S5854);
    float _S5904 = s_primal_ctx_max_0(_S5900, _S5855);
    float _S5905 = s_primal_ctx_min_0(_S5901, _S5855);
    float _S5906 = s_primal_ctx_max_0(_S5902, _S5869);
    float _S5907 = s_primal_ctx_min_0(_S5903, _S5869);
    float _S5908 = s_primal_ctx_max_0(_S5904, _S5870);
    float _S5909 = s_primal_ctx_min_0(_S5905, _S5870);
    Matrix<float, 3, 3>  _S5910 = transpose_0(R_35);
    float3  _S5911 = s_primal_ctx_mul_1(R_35, _S5763) + t_37 - - s_primal_ctx_mul_1(_S5910, t_37);
    float _S5912 = _S5911.x;
    float _S5913 = _S5911.y;
    float _S5914 = _S5911.z;
    float _S5915 = _S5912 * _S5912 + _S5913 * _S5913 + _S5914 * _S5914;
    float _S5916 = s_primal_ctx_sqrt_0(_S5915);
    float x_67 = _S5912 / _S5916;
    float3  _S5917 = make_float3 (x_67);
    float _S5918 = _S5916 * _S5916;
    float y_34 = _S5913 / _S5916;
    float z_31 = _S5914 / _S5916;
    float3  _S5919 = make_float3 (z_31);
    float _S5920 = - y_34;
    float3  _S5921 = make_float3 (_S5920);
    float z2_64 = z_31 * z_31;
    float fTmp0B_31 = -1.09254848957061768f * z_31;
    float fC1_31 = x_67 * x_67 - y_34 * y_34;
    float _S5922 = 2.0f * x_67;
    float fS1_31 = _S5922 * y_34;
    float pSH6_9 = 0.94617468118667603f * z2_64 - 0.31539157032966614f;
    float3  _S5923 = make_float3 (pSH6_9);
    float pSH7_9 = fTmp0B_31 * x_67;
    float3  _S5924 = make_float3 (pSH7_9);
    float pSH5_9 = fTmp0B_31 * y_34;
    float3  _S5925 = make_float3 (pSH5_9);
    float pSH8_9 = 0.54627424478530884f * fC1_31;
    float3  _S5926 = make_float3 (pSH8_9);
    float pSH4_9 = 0.54627424478530884f * fS1_31;
    float3  _S5927 = make_float3 (pSH4_9);
    float fTmp0C_31 = -2.28522896766662598f * z2_64 + 0.4570457935333252f;
    float fTmp1B_31 = 1.44530570507049561f * z_31;
    float _S5928 = 1.86588168144226074f * z2_64 - 1.11952900886535645f;
    float pSH12_9 = z_31 * _S5928;
    float3  _S5929 = make_float3 (pSH12_9);
    float pSH13_9 = fTmp0C_31 * x_67;
    float3  _S5930 = make_float3 (pSH13_9);
    float pSH11_9 = fTmp0C_31 * y_34;
    float3  _S5931 = make_float3 (pSH11_9);
    float pSH14_9 = fTmp1B_31 * fC1_31;
    float3  _S5932 = make_float3 (pSH14_9);
    float pSH10_9 = fTmp1B_31 * fS1_31;
    float3  _S5933 = make_float3 (pSH10_9);
    float pSH15_9 = -0.59004360437393188f * (x_67 * fC1_31 - y_34 * fS1_31);
    float3  _S5934 = make_float3 (pSH15_9);
    float pSH9_9 = -0.59004360437393188f * (x_67 * fS1_31 + y_34 * fC1_31);
    float3  _S5935 = make_float3 (pSH9_9);
    float3  _S5936 = make_float3 (0.0f);
    float3  _S5937 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5938;
    (&_S5938)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_31)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S5920) * (*sh_coeffs_31)[int(1)] + make_float3 (z_31) * (*sh_coeffs_31)[int(2)] - make_float3 (x_67) * (*sh_coeffs_31)[int(3)]) + (make_float3 (pSH4_9) * (*sh_coeffs_31)[int(4)] + make_float3 (pSH5_9) * (*sh_coeffs_31)[int(5)] + make_float3 (pSH6_9) * (*sh_coeffs_31)[int(6)] + make_float3 (pSH7_9) * (*sh_coeffs_31)[int(7)] + make_float3 (pSH8_9) * (*sh_coeffs_31)[int(8)]) + (make_float3 (pSH9_9) * (*sh_coeffs_31)[int(9)] + make_float3 (pSH10_9) * (*sh_coeffs_31)[int(10)] + make_float3 (pSH11_9) * (*sh_coeffs_31)[int(11)] + make_float3 (pSH12_9) * (*sh_coeffs_31)[int(12)] + make_float3 (pSH13_9) * (*sh_coeffs_31)[int(13)] + make_float3 (pSH14_9) * (*sh_coeffs_31)[int(14)] + make_float3 (pSH15_9) * (*sh_coeffs_31)[int(15)]) + make_float3 (0.5f);
    (&_S5938)->differential_0 = _S5937;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5939;
    (&_S5939)->primal_0 = _S5936;
    (&_S5939)->differential_0 = _S5937;
    s_bwd_prop_max_0(&_S5938, &_S5939, v_rgb_8);
    float3  _S5940 = _S5934 * _S5938.differential_0;
    float3  _S5941 = (*sh_coeffs_31)[int(15)] * _S5938.differential_0;
    float3  _S5942 = _S5932 * _S5938.differential_0;
    float3  _S5943 = (*sh_coeffs_31)[int(14)] * _S5938.differential_0;
    float3  _S5944 = _S5930 * _S5938.differential_0;
    float3  _S5945 = (*sh_coeffs_31)[int(13)] * _S5938.differential_0;
    float3  _S5946 = _S5929 * _S5938.differential_0;
    float3  _S5947 = (*sh_coeffs_31)[int(12)] * _S5938.differential_0;
    float3  _S5948 = _S5931 * _S5938.differential_0;
    float3  _S5949 = (*sh_coeffs_31)[int(11)] * _S5938.differential_0;
    float3  _S5950 = _S5933 * _S5938.differential_0;
    float3  _S5951 = (*sh_coeffs_31)[int(10)] * _S5938.differential_0;
    float3  _S5952 = _S5935 * _S5938.differential_0;
    float3  _S5953 = (*sh_coeffs_31)[int(9)] * _S5938.differential_0;
    float s_diff_fS2_T_9 = -0.59004360437393188f * (_S5953.x + _S5953.y + _S5953.z);
    float s_diff_fC2_T_9 = -0.59004360437393188f * (_S5941.x + _S5941.y + _S5941.z);
    float _S5954 = _S5951.x + _S5951.y + _S5951.z;
    float _S5955 = _S5943.x + _S5943.y + _S5943.z;
    float _S5956 = _S5949.x + _S5949.y + _S5949.z;
    float _S5957 = _S5945.x + _S5945.y + _S5945.z;
    float _S5958 = _S5947.x + _S5947.y + _S5947.z;
    float _S5959 = - s_diff_fC2_T_9;
    float3  _S5960 = _S5926 * _S5938.differential_0;
    float3  _S5961 = (*sh_coeffs_31)[int(8)] * _S5938.differential_0;
    float3  _S5962 = _S5924 * _S5938.differential_0;
    float3  _S5963 = (*sh_coeffs_31)[int(7)] * _S5938.differential_0;
    float3  _S5964 = _S5923 * _S5938.differential_0;
    float3  _S5965 = (*sh_coeffs_31)[int(6)] * _S5938.differential_0;
    float3  _S5966 = _S5925 * _S5938.differential_0;
    float3  _S5967 = (*sh_coeffs_31)[int(5)] * _S5938.differential_0;
    float3  _S5968 = _S5927 * _S5938.differential_0;
    float3  _S5969 = (*sh_coeffs_31)[int(4)] * _S5938.differential_0;
    float _S5970 = _S5967.x + _S5967.y + _S5967.z;
    float _S5971 = _S5963.x + _S5963.y + _S5963.z;
    float _S5972 = fTmp1B_31 * _S5954 + x_67 * s_diff_fS2_T_9 + y_34 * _S5959 + 0.54627424478530884f * (_S5969.x + _S5969.y + _S5969.z);
    float _S5973 = fTmp1B_31 * _S5955 + y_34 * s_diff_fS2_T_9 + x_67 * s_diff_fC2_T_9 + 0.54627424478530884f * (_S5961.x + _S5961.y + _S5961.z);
    float _S5974 = y_34 * - _S5973;
    float _S5975 = x_67 * _S5973;
    float _S5976 = z_31 * (1.86588168144226074f * (z_31 * _S5958) + -2.28522896766662598f * (y_34 * _S5956 + x_67 * _S5957) + 0.94617468118667603f * (_S5965.x + _S5965.y + _S5965.z));
    float3  _S5977 = make_float3 (0.48860251903533936f) * _S5938.differential_0;
    float3  _S5978 = - _S5977;
    float3  _S5979 = _S5917 * _S5978;
    float3  _S5980 = (*sh_coeffs_31)[int(3)] * _S5978;
    float3  _S5981 = _S5919 * _S5977;
    float3  _S5982 = (*sh_coeffs_31)[int(2)] * _S5977;
    float3  _S5983 = _S5921 * _S5977;
    float3  _S5984 = (*sh_coeffs_31)[int(1)] * _S5977;
    float _S5985 = (_S5928 * _S5958 + 1.44530570507049561f * (fS1_31 * _S5954 + fC1_31 * _S5955) + -1.09254848957061768f * (y_34 * _S5970 + x_67 * _S5971) + _S5976 + _S5976 + _S5982.x + _S5982.y + _S5982.z) / _S5918;
    float _S5986 = _S5916 * _S5985;
    float _S5987 = (fTmp0C_31 * _S5956 + fC1_31 * s_diff_fS2_T_9 + fS1_31 * _S5959 + fTmp0B_31 * _S5970 + _S5922 * _S5972 + _S5974 + _S5974 + - (_S5984.x + _S5984.y + _S5984.z)) / _S5918;
    float _S5988 = _S5916 * _S5987;
    float _S5989 = (fTmp0C_31 * _S5957 + fS1_31 * s_diff_fS2_T_9 + fC1_31 * s_diff_fC2_T_9 + fTmp0B_31 * _S5971 + 2.0f * (y_34 * _S5972) + _S5975 + _S5975 + _S5980.x + _S5980.y + _S5980.z) / _S5918;
    float _S5990 = _S5916 * _S5989;
    float _S5991 = _S5914 * - _S5985 + _S5913 * - _S5987 + _S5912 * - _S5989;
    DiffPair_float_0 _S5992;
    (&_S5992)->primal_0 = _S5915;
    (&_S5992)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5992, _S5991);
    float _S5993 = _S5914 * _S5992.differential_0;
    float _S5994 = _S5913 * _S5992.differential_0;
    float _S5995 = _S5912 * _S5992.differential_0;
    float3  _S5996 = make_float3 (0.282094806432724f) * _S5938.differential_0;
    float3  _S5997 = make_float3 (_S5990 + _S5995 + _S5995, _S5988 + _S5994 + _S5994, _S5986 + _S5993 + _S5993);
    float3  _S5998 = - - _S5997;
    Matrix<float, 3, 3>  _S5999 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6000;
    (&_S6000)->primal_0 = _S5910;
    (&_S6000)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6001;
    (&_S6001)->primal_0 = t_37;
    (&_S6001)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6000, &_S6001, _S5998);
    Matrix<float, 3, 3>  _S6002 = transpose_0(_S6000.differential_0);
    DiffPair_float_0 _S6003;
    (&_S6003)->primal_0 = _S5909;
    (&_S6003)->differential_0 = 0.0f;
    DiffPair_float_0 _S6004;
    (&_S6004)->primal_0 = _S5885;
    (&_S6004)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6003, &_S6004, 0.0f);
    DiffPair_float_0 _S6005;
    (&_S6005)->primal_0 = _S5908;
    (&_S6005)->differential_0 = 0.0f;
    DiffPair_float_0 _S6006;
    (&_S6006)->primal_0 = _S5885;
    (&_S6006)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6005, &_S6006, 0.0f);
    DiffPair_float_0 _S6007;
    (&_S6007)->primal_0 = _S5907;
    (&_S6007)->differential_0 = 0.0f;
    DiffPair_float_0 _S6008;
    (&_S6008)->primal_0 = _S5884;
    (&_S6008)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6007, &_S6008, 0.0f);
    DiffPair_float_0 _S6009;
    (&_S6009)->primal_0 = _S5906;
    (&_S6009)->differential_0 = 0.0f;
    DiffPair_float_0 _S6010;
    (&_S6010)->primal_0 = _S5884;
    (&_S6010)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6009, &_S6010, 0.0f);
    DiffPair_float_0 _S6011;
    (&_S6011)->primal_0 = _S5905;
    (&_S6011)->differential_0 = 0.0f;
    DiffPair_float_0 _S6012;
    (&_S6012)->primal_0 = _S5870;
    (&_S6012)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6011, &_S6012, _S6003.differential_0);
    DiffPair_float_0 _S6013;
    (&_S6013)->primal_0 = _S5904;
    (&_S6013)->differential_0 = 0.0f;
    DiffPair_float_0 _S6014;
    (&_S6014)->primal_0 = _S5870;
    (&_S6014)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6013, &_S6014, _S6005.differential_0);
    DiffPair_float_0 _S6015;
    (&_S6015)->primal_0 = _S5903;
    (&_S6015)->differential_0 = 0.0f;
    DiffPair_float_0 _S6016;
    (&_S6016)->primal_0 = _S5869;
    (&_S6016)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6015, &_S6016, _S6007.differential_0);
    DiffPair_float_0 _S6017;
    (&_S6017)->primal_0 = _S5902;
    (&_S6017)->differential_0 = 0.0f;
    DiffPair_float_0 _S6018;
    (&_S6018)->primal_0 = _S5869;
    (&_S6018)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6017, &_S6018, _S6009.differential_0);
    DiffPair_float_0 _S6019;
    (&_S6019)->primal_0 = _S5901;
    (&_S6019)->differential_0 = 0.0f;
    DiffPair_float_0 _S6020;
    (&_S6020)->primal_0 = _S5855;
    (&_S6020)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6019, &_S6020, _S6011.differential_0);
    DiffPair_float_0 _S6021;
    (&_S6021)->primal_0 = _S5900;
    (&_S6021)->differential_0 = 0.0f;
    DiffPair_float_0 _S6022;
    (&_S6022)->primal_0 = _S5855;
    (&_S6022)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6021, &_S6022, _S6013.differential_0);
    DiffPair_float_0 _S6023;
    (&_S6023)->primal_0 = _S5899;
    (&_S6023)->differential_0 = 0.0f;
    DiffPair_float_0 _S6024;
    (&_S6024)->primal_0 = _S5854;
    (&_S6024)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6023, &_S6024, _S6015.differential_0);
    DiffPair_float_0 _S6025;
    (&_S6025)->primal_0 = _S5898;
    (&_S6025)->differential_0 = 0.0f;
    DiffPair_float_0 _S6026;
    (&_S6026)->primal_0 = _S5854;
    (&_S6026)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6025, &_S6026, _S6017.differential_0);
    DiffPair_float_0 _S6027;
    (&_S6027)->primal_0 = _S5897;
    (&_S6027)->differential_0 = 0.0f;
    DiffPair_float_0 _S6028;
    (&_S6028)->primal_0 = _S5840;
    (&_S6028)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6027, &_S6028, _S6019.differential_0);
    DiffPair_float_0 _S6029;
    (&_S6029)->primal_0 = _S5896;
    (&_S6029)->differential_0 = 0.0f;
    DiffPair_float_0 _S6030;
    (&_S6030)->primal_0 = _S5840;
    (&_S6030)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6029, &_S6030, _S6021.differential_0);
    DiffPair_float_0 _S6031;
    (&_S6031)->primal_0 = _S5895;
    (&_S6031)->differential_0 = 0.0f;
    DiffPair_float_0 _S6032;
    (&_S6032)->primal_0 = _S5839;
    (&_S6032)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6031, &_S6032, _S6023.differential_0);
    DiffPair_float_0 _S6033;
    (&_S6033)->primal_0 = _S5894;
    (&_S6033)->differential_0 = 0.0f;
    DiffPair_float_0 _S6034;
    (&_S6034)->primal_0 = _S5839;
    (&_S6034)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6033, &_S6034, _S6025.differential_0);
    DiffPair_float_0 _S6035;
    (&_S6035)->primal_0 = _S5893;
    (&_S6035)->differential_0 = 0.0f;
    DiffPair_float_0 _S6036;
    (&_S6036)->primal_0 = _S5825;
    (&_S6036)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6035, &_S6036, _S6027.differential_0);
    DiffPair_float_0 _S6037;
    (&_S6037)->primal_0 = _S5892;
    (&_S6037)->differential_0 = 0.0f;
    DiffPair_float_0 _S6038;
    (&_S6038)->primal_0 = _S5825;
    (&_S6038)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6037, &_S6038, _S6029.differential_0);
    DiffPair_float_0 _S6039;
    (&_S6039)->primal_0 = _S5891;
    (&_S6039)->differential_0 = 0.0f;
    DiffPair_float_0 _S6040;
    (&_S6040)->primal_0 = _S5824;
    (&_S6040)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6039, &_S6040, _S6031.differential_0);
    DiffPair_float_0 _S6041;
    (&_S6041)->primal_0 = _S5890;
    (&_S6041)->differential_0 = 0.0f;
    DiffPair_float_0 _S6042;
    (&_S6042)->primal_0 = _S5824;
    (&_S6042)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6041, &_S6042, _S6033.differential_0);
    DiffPair_float_0 _S6043;
    (&_S6043)->primal_0 = _S5889;
    (&_S6043)->differential_0 = 0.0f;
    DiffPair_float_0 _S6044;
    (&_S6044)->primal_0 = _S5810;
    (&_S6044)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6043, &_S6044, _S6035.differential_0);
    DiffPair_float_0 _S6045;
    (&_S6045)->primal_0 = _S5888;
    (&_S6045)->differential_0 = 0.0f;
    DiffPair_float_0 _S6046;
    (&_S6046)->primal_0 = _S5810;
    (&_S6046)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6045, &_S6046, _S6037.differential_0);
    DiffPair_float_0 _S6047;
    (&_S6047)->primal_0 = _S5887;
    (&_S6047)->differential_0 = 0.0f;
    DiffPair_float_0 _S6048;
    (&_S6048)->primal_0 = _S5809;
    (&_S6048)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6047, &_S6048, _S6039.differential_0);
    DiffPair_float_0 _S6049;
    (&_S6049)->primal_0 = _S5886;
    (&_S6049)->differential_0 = 0.0f;
    DiffPair_float_0 _S6050;
    (&_S6050)->primal_0 = _S5809;
    (&_S6050)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6049, &_S6050, _S6041.differential_0);
    DiffPair_float_0 _S6051;
    (&_S6051)->primal_0 = _S5780;
    (&_S6051)->differential_0 = 0.0f;
    DiffPair_float_0 _S6052;
    (&_S6052)->primal_0 = _S5795;
    (&_S6052)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6051, &_S6052, _S6043.differential_0);
    DiffPair_float_0 _S6053;
    (&_S6053)->primal_0 = _S5780;
    (&_S6053)->differential_0 = 0.0f;
    DiffPair_float_0 _S6054;
    (&_S6054)->primal_0 = _S5795;
    (&_S6054)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6053, &_S6054, _S6045.differential_0);
    DiffPair_float_0 _S6055;
    (&_S6055)->primal_0 = _S5779;
    (&_S6055)->differential_0 = 0.0f;
    DiffPair_float_0 _S6056;
    (&_S6056)->primal_0 = _S5794;
    (&_S6056)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6055, &_S6056, _S6047.differential_0);
    DiffPair_float_0 _S6057;
    (&_S6057)->primal_0 = _S5779;
    (&_S6057)->differential_0 = 0.0f;
    DiffPair_float_0 _S6058;
    (&_S6058)->primal_0 = _S5794;
    (&_S6058)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6057, &_S6058, _S6049.differential_0);
    float _S6059 = fx_40 * (_S6008.differential_0 + _S6010.differential_0);
    float2  _S6060 = make_float2 (_S6059, fy_40 * (_S6004.differential_0 + _S6006.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6059, (*dist_coeffs_47)[int(9)] * _S6059);
    float2  _S6061 = _S5872 * _S6060;
    float _S6062 = (*dist_coeffs_47)[int(4)] * _S6060.y;
    float _S6063 = (*dist_coeffs_47)[int(5)] * _S6060.x;
    float _S6064 = _S6061.x + _S6061.y;
    float _S6065 = r2_119 * _S6064;
    float _S6066 = r2_119 * _S6065;
    float _S6067 = (*dist_coeffs_47)[int(7)] * _S6060.y + _S6062 + (*dist_coeffs_47)[int(6)] * _S6060.x + _S6063 + _S5876 * _S6064 + _S5875 * _S6065 + _S5874 * _S6066 + (*dist_coeffs_47)[int(3)] * (r2_119 * _S6066);
    float _S6068 = v_124 * _S6067;
    float _S6069 = u_124 * _S6067;
    float2  _S6070 = (_S5877 * _S6060 + make_float2 (_S5774 * (v_124 * _S6060.y) + _S5879 * _S6063 + 2.0f * (u_124 * _S6063) + _S5771 * (v_124 * _S6060.x) + _S6069 + _S6069, _S5881 * _S6062 + 2.0f * (v_124 * _S6062) + _S5880 * _S6060.y + _S5878 * _S6060.x + _S6068 + _S6068)) / _S5873;
    float2  _S6071 = _S5871 * - _S6070;
    float2  _S6072 = _S5762 * _S6070;
    float _S6073 = fx_40 * (_S6016.differential_0 + _S6018.differential_0);
    float2  _S6074 = make_float2 (_S6073, fy_40 * (_S6012.differential_0 + _S6014.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6073, (*dist_coeffs_47)[int(9)] * _S6073);
    float2  _S6075 = _S5857 * _S6074;
    float _S6076 = (*dist_coeffs_47)[int(4)] * _S6074.y;
    float _S6077 = (*dist_coeffs_47)[int(5)] * _S6074.x;
    float _S6078 = _S6075.x + _S6075.y;
    float _S6079 = r2_118 * _S6078;
    float _S6080 = r2_118 * _S6079;
    float _S6081 = (*dist_coeffs_47)[int(7)] * _S6074.y + _S6076 + (*dist_coeffs_47)[int(6)] * _S6074.x + _S6077 + _S5861 * _S6078 + _S5860 * _S6079 + _S5859 * _S6080 + (*dist_coeffs_47)[int(3)] * (r2_118 * _S6080);
    float _S6082 = v_123 * _S6081;
    float _S6083 = u_123 * _S6081;
    float2  _S6084 = (_S5862 * _S6074 + make_float2 (_S5774 * (v_123 * _S6074.y) + _S5864 * _S6077 + 2.0f * (u_123 * _S6077) + _S5771 * (v_123 * _S6074.x) + _S6083 + _S6083, _S5866 * _S6076 + 2.0f * (v_123 * _S6076) + _S5865 * _S6074.y + _S5863 * _S6074.x + _S6082 + _S6082)) / _S5858;
    float2  _S6085 = _S5856 * - _S6084;
    float2  _S6086 = _S5757 * _S6084;
    float _S6087 = fx_40 * (_S6024.differential_0 + _S6026.differential_0);
    float2  _S6088 = make_float2 (_S6087, fy_40 * (_S6020.differential_0 + _S6022.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6087, (*dist_coeffs_47)[int(9)] * _S6087);
    float2  _S6089 = _S5842 * _S6088;
    float _S6090 = (*dist_coeffs_47)[int(4)] * _S6088.y;
    float _S6091 = (*dist_coeffs_47)[int(5)] * _S6088.x;
    float _S6092 = _S6089.x + _S6089.y;
    float _S6093 = r2_117 * _S6092;
    float _S6094 = r2_117 * _S6093;
    float _S6095 = (*dist_coeffs_47)[int(7)] * _S6088.y + _S6090 + (*dist_coeffs_47)[int(6)] * _S6088.x + _S6091 + _S5846 * _S6092 + _S5845 * _S6093 + _S5844 * _S6094 + (*dist_coeffs_47)[int(3)] * (r2_117 * _S6094);
    float _S6096 = v_122 * _S6095;
    float _S6097 = u_122 * _S6095;
    float2  _S6098 = (_S5847 * _S6088 + make_float2 (_S5774 * (v_122 * _S6088.y) + _S5849 * _S6091 + 2.0f * (u_122 * _S6091) + _S5771 * (v_122 * _S6088.x) + _S6097 + _S6097, _S5851 * _S6090 + 2.0f * (v_122 * _S6090) + _S5850 * _S6088.y + _S5848 * _S6088.x + _S6096 + _S6096)) / _S5843;
    float2  _S6099 = _S5841 * - _S6098;
    float2  _S6100 = _S5752 * _S6098;
    float _S6101 = fx_40 * (_S6032.differential_0 + _S6034.differential_0);
    float2  _S6102 = make_float2 (_S6101, fy_40 * (_S6028.differential_0 + _S6030.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6101, (*dist_coeffs_47)[int(9)] * _S6101);
    float2  _S6103 = _S5827 * _S6102;
    float _S6104 = (*dist_coeffs_47)[int(4)] * _S6102.y;
    float _S6105 = (*dist_coeffs_47)[int(5)] * _S6102.x;
    float _S6106 = _S6103.x + _S6103.y;
    float _S6107 = r2_116 * _S6106;
    float _S6108 = r2_116 * _S6107;
    float _S6109 = (*dist_coeffs_47)[int(7)] * _S6102.y + _S6104 + (*dist_coeffs_47)[int(6)] * _S6102.x + _S6105 + _S5831 * _S6106 + _S5830 * _S6107 + _S5829 * _S6108 + (*dist_coeffs_47)[int(3)] * (r2_116 * _S6108);
    float _S6110 = v_121 * _S6109;
    float _S6111 = u_121 * _S6109;
    float2  _S6112 = (_S5832 * _S6102 + make_float2 (_S5774 * (v_121 * _S6102.y) + _S5834 * _S6105 + 2.0f * (u_121 * _S6105) + _S5771 * (v_121 * _S6102.x) + _S6111 + _S6111, _S5836 * _S6104 + 2.0f * (v_121 * _S6104) + _S5835 * _S6102.y + _S5833 * _S6102.x + _S6110 + _S6110)) / _S5828;
    float2  _S6113 = _S5826 * - _S6112;
    float2  _S6114 = _S5747 * _S6112;
    float _S6115 = fx_40 * (_S6040.differential_0 + _S6042.differential_0);
    float2  _S6116 = make_float2 (_S6115, fy_40 * (_S6036.differential_0 + _S6038.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6115, (*dist_coeffs_47)[int(9)] * _S6115);
    float2  _S6117 = _S5812 * _S6116;
    float _S6118 = (*dist_coeffs_47)[int(4)] * _S6116.y;
    float _S6119 = (*dist_coeffs_47)[int(5)] * _S6116.x;
    float _S6120 = _S6117.x + _S6117.y;
    float _S6121 = r2_115 * _S6120;
    float _S6122 = r2_115 * _S6121;
    float _S6123 = (*dist_coeffs_47)[int(7)] * _S6116.y + _S6118 + (*dist_coeffs_47)[int(6)] * _S6116.x + _S6119 + _S5816 * _S6120 + _S5815 * _S6121 + _S5814 * _S6122 + (*dist_coeffs_47)[int(3)] * (r2_115 * _S6122);
    float _S6124 = v_120 * _S6123;
    float _S6125 = u_120 * _S6123;
    float2  _S6126 = (_S5817 * _S6116 + make_float2 (_S5774 * (v_120 * _S6116.y) + _S5819 * _S6119 + 2.0f * (u_120 * _S6119) + _S5771 * (v_120 * _S6116.x) + _S6125 + _S6125, _S5821 * _S6118 + 2.0f * (v_120 * _S6118) + _S5820 * _S6116.y + _S5818 * _S6116.x + _S6124 + _S6124)) / _S5813;
    float2  _S6127 = _S5811 * - _S6126;
    float2  _S6128 = _S5742 * _S6126;
    float _S6129 = fx_40 * (_S6048.differential_0 + _S6050.differential_0);
    float2  _S6130 = make_float2 (_S6129, fy_40 * (_S6044.differential_0 + _S6046.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6129, (*dist_coeffs_47)[int(9)] * _S6129);
    float2  _S6131 = _S5797 * _S6130;
    float _S6132 = (*dist_coeffs_47)[int(4)] * _S6130.y;
    float _S6133 = (*dist_coeffs_47)[int(5)] * _S6130.x;
    float _S6134 = _S6131.x + _S6131.y;
    float _S6135 = r2_114 * _S6134;
    float _S6136 = r2_114 * _S6135;
    float _S6137 = (*dist_coeffs_47)[int(7)] * _S6130.y + _S6132 + (*dist_coeffs_47)[int(6)] * _S6130.x + _S6133 + _S5801 * _S6134 + _S5800 * _S6135 + _S5799 * _S6136 + (*dist_coeffs_47)[int(3)] * (r2_114 * _S6136);
    float _S6138 = v_119 * _S6137;
    float _S6139 = u_119 * _S6137;
    float2  _S6140 = (_S5802 * _S6130 + make_float2 (_S5774 * (v_119 * _S6130.y) + _S5804 * _S6133 + 2.0f * (u_119 * _S6133) + _S5771 * (v_119 * _S6130.x) + _S6139 + _S6139, _S5806 * _S6132 + 2.0f * (v_119 * _S6132) + _S5805 * _S6130.y + _S5803 * _S6130.x + _S6138 + _S6138)) / _S5798;
    float2  _S6141 = _S5796 * - _S6140;
    float2  _S6142 = _S5737 * _S6140;
    float _S6143 = fx_40 * (_S6056.differential_0 + _S6058.differential_0);
    float2  _S6144 = make_float2 (_S6143, fy_40 * (_S6052.differential_0 + _S6054.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6143, (*dist_coeffs_47)[int(9)] * _S6143);
    float2  _S6145 = _S5782 * _S6144;
    float _S6146 = (*dist_coeffs_47)[int(4)] * _S6144.y;
    float _S6147 = (*dist_coeffs_47)[int(5)] * _S6144.x;
    float _S6148 = _S6145.x + _S6145.y;
    float _S6149 = r2_113 * _S6148;
    float _S6150 = r2_113 * _S6149;
    float _S6151 = (*dist_coeffs_47)[int(7)] * _S6144.y + _S6146 + (*dist_coeffs_47)[int(6)] * _S6144.x + _S6147 + _S5786 * _S6148 + _S5785 * _S6149 + _S5784 * _S6150 + (*dist_coeffs_47)[int(3)] * (r2_113 * _S6150);
    float _S6152 = v_118 * _S6151;
    float _S6153 = u_118 * _S6151;
    float2  _S6154 = (_S5787 * _S6144 + make_float2 (_S5774 * (v_118 * _S6144.y) + _S5789 * _S6147 + 2.0f * (u_118 * _S6147) + _S5771 * (v_118 * _S6144.x) + _S6153 + _S6153, _S5791 * _S6146 + 2.0f * (v_118 * _S6146) + _S5790 * _S6144.y + _S5788 * _S6144.x + _S6152 + _S6152)) / _S5783;
    float2  _S6155 = _S5781 * - _S6154;
    float2  _S6156 = _S5732 * _S6154;
    float _S6157 = fx_40 * (_S6055.differential_0 + _S6057.differential_0);
    float2  _S6158 = make_float2 (_S6157, fy_40 * (_S6051.differential_0 + _S6053.differential_0)) + make_float2 ((*dist_coeffs_47)[int(8)] * _S6157, (*dist_coeffs_47)[int(9)] * _S6157);
    float2  _S6159 = _S5765 * _S6158;
    float _S6160 = (*dist_coeffs_47)[int(4)] * _S6158.y;
    float _S6161 = (*dist_coeffs_47)[int(5)] * _S6158.x;
    float _S6162 = _S6159.x + _S6159.y;
    float _S6163 = r2_112 * _S6162;
    float _S6164 = r2_112 * _S6163;
    float _S6165 = (*dist_coeffs_47)[int(7)] * _S6158.y + _S6160 + (*dist_coeffs_47)[int(6)] * _S6158.x + _S6161 + _S5769 * _S6162 + _S5768 * _S6163 + _S5767 * _S6164 + (*dist_coeffs_47)[int(3)] * (r2_112 * _S6164);
    float _S6166 = v_117 * _S6165;
    float _S6167 = u_117 * _S6165;
    float2  _S6168 = (_S5770 * _S6158 + make_float2 (_S5774 * (v_117 * _S6158.y) + _S5773 * _S6161 + 2.0f * (u_117 * _S6161) + _S5771 * (v_117 * _S6158.x) + _S6167 + _S6167, _S5776 * _S6160 + 2.0f * (v_117 * _S6160) + _S5775 * _S6158.y + _S5772 * _S6158.x + _S6166 + _S6166)) / _S5766;
    float2  _S6169 = _S5764 * - _S6168;
    float2  _S6170 = _S5727 * _S6168;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6171;
    (&_S6171)->primal_0 = R_35;
    (&_S6171)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6172;
    (&_S6172)->primal_0 = _S5763;
    (&_S6172)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6171, &_S6172, _S5997);
    DiffPair_float_0 _S6173;
    (&_S6173)->primal_0 = _S5759;
    (&_S6173)->differential_0 = 0.0f;
    DiffPair_float_0 _S6174;
    (&_S6174)->primal_0 = _S5761;
    (&_S6174)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6173, &_S6174, 0.0f);
    DiffPair_float_0 _S6175;
    (&_S6175)->primal_0 = _S5758;
    (&_S6175)->differential_0 = 0.0f;
    DiffPair_float_0 _S6176;
    (&_S6176)->primal_0 = _S5761;
    (&_S6176)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6175, &_S6176, 0.0f);
    float3  _S6177 = make_float3 (_S6072.x, _S6072.y, _S6174.differential_0 + _S6176.differential_0 + _S6071.x + _S6071.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6178;
    (&_S6178)->primal_0 = R_35;
    (&_S6178)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6179;
    (&_S6179)->primal_0 = pos_i_6;
    (&_S6179)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6178, &_S6179, _S6177);
    DiffPair_float_0 _S6180;
    (&_S6180)->primal_0 = _S5754;
    (&_S6180)->differential_0 = 0.0f;
    DiffPair_float_0 _S6181;
    (&_S6181)->primal_0 = _S5756;
    (&_S6181)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6180, &_S6181, _S6173.differential_0);
    DiffPair_float_0 _S6182;
    (&_S6182)->primal_0 = _S5753;
    (&_S6182)->differential_0 = 0.0f;
    DiffPair_float_0 _S6183;
    (&_S6183)->primal_0 = _S5756;
    (&_S6183)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6182, &_S6183, _S6175.differential_0);
    float3  _S6184 = make_float3 (_S6086.x, _S6086.y, _S6181.differential_0 + _S6183.differential_0 + _S6085.x + _S6085.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6185;
    (&_S6185)->primal_0 = R_35;
    (&_S6185)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6186;
    (&_S6186)->primal_0 = pos_i_5;
    (&_S6186)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6185, &_S6186, _S6184);
    DiffPair_float_0 _S6187;
    (&_S6187)->primal_0 = _S5749;
    (&_S6187)->differential_0 = 0.0f;
    DiffPair_float_0 _S6188;
    (&_S6188)->primal_0 = _S5751;
    (&_S6188)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6187, &_S6188, _S6180.differential_0);
    DiffPair_float_0 _S6189;
    (&_S6189)->primal_0 = _S5748;
    (&_S6189)->differential_0 = 0.0f;
    DiffPair_float_0 _S6190;
    (&_S6190)->primal_0 = _S5751;
    (&_S6190)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6189, &_S6190, _S6182.differential_0);
    float3  _S6191 = make_float3 (_S6100.x, _S6100.y, _S6188.differential_0 + _S6190.differential_0 + _S6099.x + _S6099.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6192;
    (&_S6192)->primal_0 = R_35;
    (&_S6192)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6193;
    (&_S6193)->primal_0 = pos_i_4;
    (&_S6193)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6192, &_S6193, _S6191);
    DiffPair_float_0 _S6194;
    (&_S6194)->primal_0 = _S5744;
    (&_S6194)->differential_0 = 0.0f;
    DiffPair_float_0 _S6195;
    (&_S6195)->primal_0 = _S5746;
    (&_S6195)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6194, &_S6195, _S6187.differential_0);
    DiffPair_float_0 _S6196;
    (&_S6196)->primal_0 = _S5743;
    (&_S6196)->differential_0 = 0.0f;
    DiffPair_float_0 _S6197;
    (&_S6197)->primal_0 = _S5746;
    (&_S6197)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6196, &_S6197, _S6189.differential_0);
    float3  _S6198 = make_float3 (_S6114.x, _S6114.y, _S6195.differential_0 + _S6197.differential_0 + _S6113.x + _S6113.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6199;
    (&_S6199)->primal_0 = R_35;
    (&_S6199)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6200;
    (&_S6200)->primal_0 = pos_i_3;
    (&_S6200)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6199, &_S6200, _S6198);
    DiffPair_float_0 _S6201;
    (&_S6201)->primal_0 = _S5739;
    (&_S6201)->differential_0 = 0.0f;
    DiffPair_float_0 _S6202;
    (&_S6202)->primal_0 = _S5741;
    (&_S6202)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6201, &_S6202, _S6194.differential_0);
    DiffPair_float_0 _S6203;
    (&_S6203)->primal_0 = _S5738;
    (&_S6203)->differential_0 = 0.0f;
    DiffPair_float_0 _S6204;
    (&_S6204)->primal_0 = _S5741;
    (&_S6204)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6203, &_S6204, _S6196.differential_0);
    float3  _S6205 = make_float3 (_S6128.x, _S6128.y, _S6202.differential_0 + _S6204.differential_0 + _S6127.x + _S6127.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6206;
    (&_S6206)->primal_0 = R_35;
    (&_S6206)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6207;
    (&_S6207)->primal_0 = pos_i_2;
    (&_S6207)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6206, &_S6207, _S6205);
    DiffPair_float_0 _S6208;
    (&_S6208)->primal_0 = _S5734;
    (&_S6208)->differential_0 = 0.0f;
    DiffPair_float_0 _S6209;
    (&_S6209)->primal_0 = _S5736;
    (&_S6209)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6208, &_S6209, _S6201.differential_0);
    DiffPair_float_0 _S6210;
    (&_S6210)->primal_0 = _S5733;
    (&_S6210)->differential_0 = 0.0f;
    DiffPair_float_0 _S6211;
    (&_S6211)->primal_0 = _S5736;
    (&_S6211)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6210, &_S6211, _S6203.differential_0);
    float3  _S6212 = make_float3 (_S6142.x, _S6142.y, _S6209.differential_0 + _S6211.differential_0 + _S6141.x + _S6141.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6213;
    (&_S6213)->primal_0 = R_35;
    (&_S6213)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6214;
    (&_S6214)->primal_0 = pos_i_1;
    (&_S6214)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6213, &_S6214, _S6212);
    DiffPair_float_0 _S6215;
    (&_S6215)->primal_0 = _S5729;
    (&_S6215)->differential_0 = 0.0f;
    DiffPair_float_0 _S6216;
    (&_S6216)->primal_0 = _S5731;
    (&_S6216)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6215, &_S6216, _S6208.differential_0);
    DiffPair_float_0 _S6217;
    (&_S6217)->primal_0 = _S5728;
    (&_S6217)->differential_0 = 0.0f;
    DiffPair_float_0 _S6218;
    (&_S6218)->primal_0 = _S5731;
    (&_S6218)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6217, &_S6218, _S6210.differential_0);
    float3  _S6219 = make_float3 (_S6156.x, _S6156.y, _S6216.differential_0 + _S6218.differential_0 + _S6155.x + _S6155.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6220;
    (&_S6220)->primal_0 = R_35;
    (&_S6220)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6221;
    (&_S6221)->primal_0 = pos_i_0;
    (&_S6221)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6220, &_S6221, _S6219);
    DiffPair_float_0 _S6222;
    (&_S6222)->primal_0 = 0.0f;
    (&_S6222)->differential_0 = 0.0f;
    DiffPair_float_0 _S6223;
    (&_S6223)->primal_0 = _S5726;
    (&_S6223)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6222, &_S6223, _S6215.differential_0);
    DiffPair_float_0 _S6224;
    (&_S6224)->primal_0 = 1.00000001504746622e+30f;
    (&_S6224)->differential_0 = 0.0f;
    DiffPair_float_0 _S6225;
    (&_S6225)->primal_0 = _S5726;
    (&_S6225)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6224, &_S6225, _S6217.differential_0);
    float3  _S6226 = make_float3 (_S6170.x, _S6170.y, _S6223.differential_0 + _S6225.differential_0 + _S6169.x + _S6169.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6227;
    (&_S6227)->primal_0 = R_35;
    (&_S6227)->differential_0 = _S5999;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6228;
    (&_S6228)->primal_0 = pos_4;
    (&_S6228)->differential_0 = _S5937;
    s_bwd_prop_mul_1(&_S6227, &_S6228, _S6226);
    float3  _S6229 = _S6001.differential_0 + _S5997 + _S6177 + _S6184 + _S6191 + _S6198 + _S6205 + _S6212 + _S6219 + _S6226;
    Matrix<float, 3, 3>  _S6230 = _S6002 + _S6171.differential_0 + _S6178.differential_0 + _S6185.differential_0 + _S6192.differential_0 + _S6199.differential_0 + _S6206.differential_0 + _S6213.differential_0 + _S6220.differential_0 + _S6227.differential_0;
    FixedArray<float3 , 16>  _S6231;
    _S6231[int(0)] = _S5937;
    _S6231[int(1)] = _S5937;
    _S6231[int(2)] = _S5937;
    _S6231[int(3)] = _S5937;
    _S6231[int(4)] = _S5937;
    _S6231[int(5)] = _S5937;
    _S6231[int(6)] = _S5937;
    _S6231[int(7)] = _S5937;
    _S6231[int(8)] = _S5937;
    _S6231[int(9)] = _S5937;
    _S6231[int(10)] = _S5937;
    _S6231[int(11)] = _S5937;
    _S6231[int(12)] = _S5937;
    _S6231[int(13)] = _S5937;
    _S6231[int(14)] = _S5937;
    _S6231[int(15)] = _S5937;
    _S6231[int(15)] = _S5940;
    _S6231[int(14)] = _S5942;
    _S6231[int(13)] = _S5944;
    _S6231[int(12)] = _S5946;
    _S6231[int(11)] = _S5948;
    _S6231[int(10)] = _S5950;
    _S6231[int(9)] = _S5952;
    _S6231[int(8)] = _S5960;
    _S6231[int(7)] = _S5962;
    _S6231[int(6)] = _S5964;
    _S6231[int(5)] = _S5966;
    _S6231[int(4)] = _S5968;
    _S6231[int(3)] = _S5979;
    _S6231[int(2)] = _S5981;
    _S6231[int(1)] = _S5983;
    _S6231[int(0)] = _S5996;
    (*v_densities_0)[int(0)] = 0.0f;
    (*v_densities_0)[int(1)] = 0.0f;
    (*v_densities_0)[int(2)] = 0.0f;
    (*v_densities_0)[int(3)] = 0.0f;
    (*v_densities_0)[int(4)] = 0.0f;
    (*v_densities_0)[int(5)] = 0.0f;
    (*v_densities_0)[int(6)] = 0.0f;
    (*v_densities_0)[int(7)] = 0.0f;
    *v_sh_coeffs_9 = _S6231;
    *v_R_10 = _S6230;
    *v_t_10 = _S6229;
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye_vjp(float3  pos_5, float size_5, FixedArray<float, 8>  * densities_5, FixedArray<float3 , 16>  * sh_coeffs_32, Matrix<float, 3, 3>  R_36, float3  t_38, float fx_41, float fy_41, float cx_36, float cy_36, FixedArray<float, 10>  * dist_coeffs_48, uint image_width_32, uint image_height_32, float3  v_rgb_9, FixedArray<float, 8>  * v_densities_1, FixedArray<float3 , 16>  * v_sh_coeffs_10, Matrix<float, 3, 3>  * v_R_11, float3  * v_t_11)
{
    float3  _S6232 = s_primal_ctx_mul_1(R_36, pos_5) + t_38;
    float _S6233 = length_1(_S6232);
    float _S6234 = s_primal_ctx_min_0(1.00000001504746622e+30f, _S6233);
    float _S6235 = s_primal_ctx_max_0(0.0f, _S6233);
    float3  pos_i_7 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S6236 = s_primal_ctx_mul_1(R_36, pos_i_7) + t_38;
    float _S6237 = length_1(_S6236);
    float _S6238 = s_primal_ctx_min_0(_S6234, _S6237);
    float _S6239 = s_primal_ctx_max_0(_S6235, _S6237);
    float3  pos_i_8 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S6240 = s_primal_ctx_mul_1(R_36, pos_i_8) + t_38;
    float _S6241 = length_1(_S6240);
    float _S6242 = s_primal_ctx_min_0(_S6238, _S6241);
    float _S6243 = s_primal_ctx_max_0(_S6239, _S6241);
    float3  pos_i_9 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S6244 = s_primal_ctx_mul_1(R_36, pos_i_9) + t_38;
    float _S6245 = length_1(_S6244);
    float _S6246 = s_primal_ctx_min_0(_S6242, _S6245);
    float _S6247 = s_primal_ctx_max_0(_S6243, _S6245);
    float3  pos_i_10 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S6248 = s_primal_ctx_mul_1(R_36, pos_i_10) + t_38;
    float _S6249 = length_1(_S6248);
    float _S6250 = s_primal_ctx_min_0(_S6246, _S6249);
    float _S6251 = s_primal_ctx_max_0(_S6247, _S6249);
    float3  pos_i_11 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S6252 = s_primal_ctx_mul_1(R_36, pos_i_11) + t_38;
    float _S6253 = length_1(_S6252);
    float _S6254 = s_primal_ctx_min_0(_S6250, _S6253);
    float _S6255 = s_primal_ctx_max_0(_S6251, _S6253);
    float3  pos_i_12 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S6256 = s_primal_ctx_mul_1(R_36, pos_i_12) + t_38;
    float _S6257 = length_1(_S6256);
    float _S6258 = s_primal_ctx_min_0(_S6254, _S6257);
    float _S6259 = s_primal_ctx_max_0(_S6255, _S6257);
    float3  pos_i_13 = pos_5 + make_float3 (size_5);
    float3  _S6260 = s_primal_ctx_mul_1(R_36, pos_i_13) + t_38;
    float _S6261 = length_1(_S6260);
    float3  _S6262 = pos_5 + make_float3 (0.5f * size_5);
    float3  mean_c_31 = s_primal_ctx_mul_1(R_36, _S6262) + t_38;
    float2  _S6263 = float2 {_S6232.x, _S6232.y};
    float _S6264 = length_0(_S6263);
    float _S6265 = _S6232.z;
    float _S6266 = s_primal_ctx_atan2_0(_S6264, _S6265);
    bool _S6267 = _S6266 < 0.00100000004749745f;
    float k_15;
    float _S6268;
    float _S6269;
    float _S6270;
    if(_S6267)
    {
        float _S6271 = 1.0f - _S6266 * _S6266 / 3.0f;
        float _S6272 = _S6265 * _S6265;
        k_15 = _S6271 / _S6265;
        _S6268 = 0.0f;
        _S6269 = _S6272;
        _S6270 = _S6271;
    }
    else
    {
        float _S6273 = _S6264 * _S6264;
        k_15 = _S6266 / _S6264;
        _S6268 = _S6273;
        _S6269 = 0.0f;
        _S6270 = 0.0f;
    }
    float2  _S6274 = make_float2 (k_15);
    float2  _S6275 = _S6263 * make_float2 (k_15);
    float u_125 = _S6275.x;
    float v_125 = _S6275.y;
    float r2_120 = u_125 * u_125 + v_125 * v_125;
    float _S6276 = (*dist_coeffs_48)[int(2)] + r2_120 * (*dist_coeffs_48)[int(3)];
    float _S6277 = (*dist_coeffs_48)[int(1)] + r2_120 * _S6276;
    float _S6278 = (*dist_coeffs_48)[int(0)] + r2_120 * _S6277;
    float radial_21 = 1.0f + r2_120 * _S6278;
    float2  _S6279 = make_float2 (radial_21);
    float _S6280 = 2.0f * (*dist_coeffs_48)[int(4)];
    float _S6281 = _S6280 * u_125;
    float _S6282 = 2.0f * u_125;
    float _S6283 = 2.0f * (*dist_coeffs_48)[int(5)];
    float _S6284 = _S6283 * u_125;
    float _S6285 = 2.0f * v_125;
    float2  _S6286 = _S6275 * make_float2 (radial_21) + make_float2 (_S6281 * v_125 + (*dist_coeffs_48)[int(5)] * (r2_120 + _S6282 * u_125) + (*dist_coeffs_48)[int(6)] * r2_120, _S6284 * v_125 + (*dist_coeffs_48)[int(4)] * (r2_120 + _S6285 * v_125) + (*dist_coeffs_48)[int(7)] * r2_120);
    float2  _S6287 = _S6286 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6286.x + (*dist_coeffs_48)[int(9)] * _S6286.y, 0.0f);
    float _S6288 = fx_41 * _S6287.x + cx_36;
    float _S6289 = fy_41 * _S6287.y + cy_36;
    float2  _S6290 = float2 {_S6236.x, _S6236.y};
    float _S6291 = length_0(_S6290);
    float _S6292 = _S6236.z;
    float _S6293 = s_primal_ctx_atan2_0(_S6291, _S6292);
    bool _S6294 = _S6293 < 0.00100000004749745f;
    float _S6295;
    float _S6296;
    float _S6297;
    if(_S6294)
    {
        float _S6298 = 1.0f - _S6293 * _S6293 / 3.0f;
        float _S6299 = _S6292 * _S6292;
        k_15 = _S6298 / _S6292;
        _S6295 = 0.0f;
        _S6296 = _S6299;
        _S6297 = _S6298;
    }
    else
    {
        float _S6300 = _S6291 * _S6291;
        k_15 = _S6293 / _S6291;
        _S6295 = _S6300;
        _S6296 = 0.0f;
        _S6297 = 0.0f;
    }
    float2  _S6301 = make_float2 (k_15);
    float2  _S6302 = _S6290 * make_float2 (k_15);
    float u_126 = _S6302.x;
    float v_126 = _S6302.y;
    float r2_121 = u_126 * u_126 + v_126 * v_126;
    float _S6303 = (*dist_coeffs_48)[int(2)] + r2_121 * (*dist_coeffs_48)[int(3)];
    float _S6304 = (*dist_coeffs_48)[int(1)] + r2_121 * _S6303;
    float _S6305 = (*dist_coeffs_48)[int(0)] + r2_121 * _S6304;
    float radial_22 = 1.0f + r2_121 * _S6305;
    float2  _S6306 = make_float2 (radial_22);
    float _S6307 = _S6280 * u_126;
    float _S6308 = 2.0f * u_126;
    float _S6309 = _S6283 * u_126;
    float _S6310 = 2.0f * v_126;
    float2  _S6311 = _S6302 * make_float2 (radial_22) + make_float2 (_S6307 * v_126 + (*dist_coeffs_48)[int(5)] * (r2_121 + _S6308 * u_126) + (*dist_coeffs_48)[int(6)] * r2_121, _S6309 * v_126 + (*dist_coeffs_48)[int(4)] * (r2_121 + _S6310 * v_126) + (*dist_coeffs_48)[int(7)] * r2_121);
    float2  _S6312 = _S6311 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6311.x + (*dist_coeffs_48)[int(9)] * _S6311.y, 0.0f);
    float _S6313 = fx_41 * _S6312.x + cx_36;
    float _S6314 = fy_41 * _S6312.y + cy_36;
    float2  _S6315 = float2 {_S6240.x, _S6240.y};
    float _S6316 = length_0(_S6315);
    float _S6317 = _S6240.z;
    float _S6318 = s_primal_ctx_atan2_0(_S6316, _S6317);
    bool _S6319 = _S6318 < 0.00100000004749745f;
    float _S6320;
    float _S6321;
    float _S6322;
    if(_S6319)
    {
        float _S6323 = 1.0f - _S6318 * _S6318 / 3.0f;
        float _S6324 = _S6317 * _S6317;
        k_15 = _S6323 / _S6317;
        _S6320 = 0.0f;
        _S6321 = _S6324;
        _S6322 = _S6323;
    }
    else
    {
        float _S6325 = _S6316 * _S6316;
        k_15 = _S6318 / _S6316;
        _S6320 = _S6325;
        _S6321 = 0.0f;
        _S6322 = 0.0f;
    }
    float2  _S6326 = make_float2 (k_15);
    float2  _S6327 = _S6315 * make_float2 (k_15);
    float u_127 = _S6327.x;
    float v_127 = _S6327.y;
    float r2_122 = u_127 * u_127 + v_127 * v_127;
    float _S6328 = (*dist_coeffs_48)[int(2)] + r2_122 * (*dist_coeffs_48)[int(3)];
    float _S6329 = (*dist_coeffs_48)[int(1)] + r2_122 * _S6328;
    float _S6330 = (*dist_coeffs_48)[int(0)] + r2_122 * _S6329;
    float radial_23 = 1.0f + r2_122 * _S6330;
    float2  _S6331 = make_float2 (radial_23);
    float _S6332 = _S6280 * u_127;
    float _S6333 = 2.0f * u_127;
    float _S6334 = _S6283 * u_127;
    float _S6335 = 2.0f * v_127;
    float2  _S6336 = _S6327 * make_float2 (radial_23) + make_float2 (_S6332 * v_127 + (*dist_coeffs_48)[int(5)] * (r2_122 + _S6333 * u_127) + (*dist_coeffs_48)[int(6)] * r2_122, _S6334 * v_127 + (*dist_coeffs_48)[int(4)] * (r2_122 + _S6335 * v_127) + (*dist_coeffs_48)[int(7)] * r2_122);
    float2  _S6337 = _S6336 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6336.x + (*dist_coeffs_48)[int(9)] * _S6336.y, 0.0f);
    float _S6338 = fx_41 * _S6337.x + cx_36;
    float _S6339 = fy_41 * _S6337.y + cy_36;
    float2  _S6340 = float2 {_S6244.x, _S6244.y};
    float _S6341 = length_0(_S6340);
    float _S6342 = _S6244.z;
    float _S6343 = s_primal_ctx_atan2_0(_S6341, _S6342);
    bool _S6344 = _S6343 < 0.00100000004749745f;
    float _S6345;
    float _S6346;
    float _S6347;
    if(_S6344)
    {
        float _S6348 = 1.0f - _S6343 * _S6343 / 3.0f;
        float _S6349 = _S6342 * _S6342;
        k_15 = _S6348 / _S6342;
        _S6345 = 0.0f;
        _S6346 = _S6349;
        _S6347 = _S6348;
    }
    else
    {
        float _S6350 = _S6341 * _S6341;
        k_15 = _S6343 / _S6341;
        _S6345 = _S6350;
        _S6346 = 0.0f;
        _S6347 = 0.0f;
    }
    float2  _S6351 = make_float2 (k_15);
    float2  _S6352 = _S6340 * make_float2 (k_15);
    float u_128 = _S6352.x;
    float v_128 = _S6352.y;
    float r2_123 = u_128 * u_128 + v_128 * v_128;
    float _S6353 = (*dist_coeffs_48)[int(2)] + r2_123 * (*dist_coeffs_48)[int(3)];
    float _S6354 = (*dist_coeffs_48)[int(1)] + r2_123 * _S6353;
    float _S6355 = (*dist_coeffs_48)[int(0)] + r2_123 * _S6354;
    float radial_24 = 1.0f + r2_123 * _S6355;
    float2  _S6356 = make_float2 (radial_24);
    float _S6357 = _S6280 * u_128;
    float _S6358 = 2.0f * u_128;
    float _S6359 = _S6283 * u_128;
    float _S6360 = 2.0f * v_128;
    float2  _S6361 = _S6352 * make_float2 (radial_24) + make_float2 (_S6357 * v_128 + (*dist_coeffs_48)[int(5)] * (r2_123 + _S6358 * u_128) + (*dist_coeffs_48)[int(6)] * r2_123, _S6359 * v_128 + (*dist_coeffs_48)[int(4)] * (r2_123 + _S6360 * v_128) + (*dist_coeffs_48)[int(7)] * r2_123);
    float2  _S6362 = _S6361 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6361.x + (*dist_coeffs_48)[int(9)] * _S6361.y, 0.0f);
    float _S6363 = fx_41 * _S6362.x + cx_36;
    float _S6364 = fy_41 * _S6362.y + cy_36;
    float2  _S6365 = float2 {_S6248.x, _S6248.y};
    float _S6366 = length_0(_S6365);
    float _S6367 = _S6248.z;
    float _S6368 = s_primal_ctx_atan2_0(_S6366, _S6367);
    bool _S6369 = _S6368 < 0.00100000004749745f;
    float _S6370;
    float _S6371;
    float _S6372;
    if(_S6369)
    {
        float _S6373 = 1.0f - _S6368 * _S6368 / 3.0f;
        float _S6374 = _S6367 * _S6367;
        k_15 = _S6373 / _S6367;
        _S6370 = 0.0f;
        _S6371 = _S6374;
        _S6372 = _S6373;
    }
    else
    {
        float _S6375 = _S6366 * _S6366;
        k_15 = _S6368 / _S6366;
        _S6370 = _S6375;
        _S6371 = 0.0f;
        _S6372 = 0.0f;
    }
    float2  _S6376 = make_float2 (k_15);
    float2  _S6377 = _S6365 * make_float2 (k_15);
    float u_129 = _S6377.x;
    float v_129 = _S6377.y;
    float r2_124 = u_129 * u_129 + v_129 * v_129;
    float _S6378 = (*dist_coeffs_48)[int(2)] + r2_124 * (*dist_coeffs_48)[int(3)];
    float _S6379 = (*dist_coeffs_48)[int(1)] + r2_124 * _S6378;
    float _S6380 = (*dist_coeffs_48)[int(0)] + r2_124 * _S6379;
    float radial_25 = 1.0f + r2_124 * _S6380;
    float2  _S6381 = make_float2 (radial_25);
    float _S6382 = _S6280 * u_129;
    float _S6383 = 2.0f * u_129;
    float _S6384 = _S6283 * u_129;
    float _S6385 = 2.0f * v_129;
    float2  _S6386 = _S6377 * make_float2 (radial_25) + make_float2 (_S6382 * v_129 + (*dist_coeffs_48)[int(5)] * (r2_124 + _S6383 * u_129) + (*dist_coeffs_48)[int(6)] * r2_124, _S6384 * v_129 + (*dist_coeffs_48)[int(4)] * (r2_124 + _S6385 * v_129) + (*dist_coeffs_48)[int(7)] * r2_124);
    float2  _S6387 = _S6386 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6386.x + (*dist_coeffs_48)[int(9)] * _S6386.y, 0.0f);
    float _S6388 = fx_41 * _S6387.x + cx_36;
    float _S6389 = fy_41 * _S6387.y + cy_36;
    float2  _S6390 = float2 {_S6252.x, _S6252.y};
    float _S6391 = length_0(_S6390);
    float _S6392 = _S6252.z;
    float _S6393 = s_primal_ctx_atan2_0(_S6391, _S6392);
    bool _S6394 = _S6393 < 0.00100000004749745f;
    float _S6395;
    float _S6396;
    float _S6397;
    if(_S6394)
    {
        float _S6398 = 1.0f - _S6393 * _S6393 / 3.0f;
        float _S6399 = _S6392 * _S6392;
        k_15 = _S6398 / _S6392;
        _S6395 = 0.0f;
        _S6396 = _S6399;
        _S6397 = _S6398;
    }
    else
    {
        float _S6400 = _S6391 * _S6391;
        k_15 = _S6393 / _S6391;
        _S6395 = _S6400;
        _S6396 = 0.0f;
        _S6397 = 0.0f;
    }
    float2  _S6401 = make_float2 (k_15);
    float2  _S6402 = _S6390 * make_float2 (k_15);
    float u_130 = _S6402.x;
    float v_130 = _S6402.y;
    float r2_125 = u_130 * u_130 + v_130 * v_130;
    float _S6403 = (*dist_coeffs_48)[int(2)] + r2_125 * (*dist_coeffs_48)[int(3)];
    float _S6404 = (*dist_coeffs_48)[int(1)] + r2_125 * _S6403;
    float _S6405 = (*dist_coeffs_48)[int(0)] + r2_125 * _S6404;
    float radial_26 = 1.0f + r2_125 * _S6405;
    float2  _S6406 = make_float2 (radial_26);
    float _S6407 = _S6280 * u_130;
    float _S6408 = 2.0f * u_130;
    float _S6409 = _S6283 * u_130;
    float _S6410 = 2.0f * v_130;
    float2  _S6411 = _S6402 * make_float2 (radial_26) + make_float2 (_S6407 * v_130 + (*dist_coeffs_48)[int(5)] * (r2_125 + _S6408 * u_130) + (*dist_coeffs_48)[int(6)] * r2_125, _S6409 * v_130 + (*dist_coeffs_48)[int(4)] * (r2_125 + _S6410 * v_130) + (*dist_coeffs_48)[int(7)] * r2_125);
    float2  _S6412 = _S6411 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6411.x + (*dist_coeffs_48)[int(9)] * _S6411.y, 0.0f);
    float _S6413 = fx_41 * _S6412.x + cx_36;
    float _S6414 = fy_41 * _S6412.y + cy_36;
    float2  _S6415 = float2 {_S6256.x, _S6256.y};
    float _S6416 = length_0(_S6415);
    float _S6417 = _S6256.z;
    float _S6418 = s_primal_ctx_atan2_0(_S6416, _S6417);
    bool _S6419 = _S6418 < 0.00100000004749745f;
    float _S6420;
    float _S6421;
    float _S6422;
    if(_S6419)
    {
        float _S6423 = 1.0f - _S6418 * _S6418 / 3.0f;
        float _S6424 = _S6417 * _S6417;
        k_15 = _S6423 / _S6417;
        _S6420 = 0.0f;
        _S6421 = _S6424;
        _S6422 = _S6423;
    }
    else
    {
        float _S6425 = _S6416 * _S6416;
        k_15 = _S6418 / _S6416;
        _S6420 = _S6425;
        _S6421 = 0.0f;
        _S6422 = 0.0f;
    }
    float2  _S6426 = make_float2 (k_15);
    float2  _S6427 = _S6415 * make_float2 (k_15);
    float u_131 = _S6427.x;
    float v_131 = _S6427.y;
    float r2_126 = u_131 * u_131 + v_131 * v_131;
    float _S6428 = (*dist_coeffs_48)[int(2)] + r2_126 * (*dist_coeffs_48)[int(3)];
    float _S6429 = (*dist_coeffs_48)[int(1)] + r2_126 * _S6428;
    float _S6430 = (*dist_coeffs_48)[int(0)] + r2_126 * _S6429;
    float radial_27 = 1.0f + r2_126 * _S6430;
    float2  _S6431 = make_float2 (radial_27);
    float _S6432 = _S6280 * u_131;
    float _S6433 = 2.0f * u_131;
    float _S6434 = _S6283 * u_131;
    float _S6435 = 2.0f * v_131;
    float2  _S6436 = _S6427 * make_float2 (radial_27) + make_float2 (_S6432 * v_131 + (*dist_coeffs_48)[int(5)] * (r2_126 + _S6433 * u_131) + (*dist_coeffs_48)[int(6)] * r2_126, _S6434 * v_131 + (*dist_coeffs_48)[int(4)] * (r2_126 + _S6435 * v_131) + (*dist_coeffs_48)[int(7)] * r2_126);
    float2  _S6437 = _S6436 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6436.x + (*dist_coeffs_48)[int(9)] * _S6436.y, 0.0f);
    float _S6438 = fx_41 * _S6437.x + cx_36;
    float _S6439 = fy_41 * _S6437.y + cy_36;
    float2  _S6440 = float2 {_S6260.x, _S6260.y};
    float _S6441 = length_0(_S6440);
    float _S6442 = _S6260.z;
    float _S6443 = s_primal_ctx_atan2_0(_S6441, _S6442);
    bool _S6444 = _S6443 < 0.00100000004749745f;
    float _S6445;
    float _S6446;
    float _S6447;
    if(_S6444)
    {
        float _S6448 = 1.0f - _S6443 * _S6443 / 3.0f;
        float _S6449 = _S6442 * _S6442;
        k_15 = _S6448 / _S6442;
        _S6445 = 0.0f;
        _S6446 = _S6449;
        _S6447 = _S6448;
    }
    else
    {
        float _S6450 = _S6441 * _S6441;
        k_15 = _S6443 / _S6441;
        _S6445 = _S6450;
        _S6446 = 0.0f;
        _S6447 = 0.0f;
    }
    float2  _S6451 = make_float2 (k_15);
    float2  _S6452 = _S6440 * make_float2 (k_15);
    float u_132 = _S6452.x;
    float v_132 = _S6452.y;
    float r2_127 = u_132 * u_132 + v_132 * v_132;
    float _S6453 = (*dist_coeffs_48)[int(2)] + r2_127 * (*dist_coeffs_48)[int(3)];
    float _S6454 = (*dist_coeffs_48)[int(1)] + r2_127 * _S6453;
    float _S6455 = (*dist_coeffs_48)[int(0)] + r2_127 * _S6454;
    float radial_28 = 1.0f + r2_127 * _S6455;
    float2  _S6456 = make_float2 (radial_28);
    float _S6457 = _S6280 * u_132;
    float _S6458 = 2.0f * u_132;
    float _S6459 = _S6283 * u_132;
    float _S6460 = 2.0f * v_132;
    float2  _S6461 = _S6452 * make_float2 (radial_28) + make_float2 (_S6457 * v_132 + (*dist_coeffs_48)[int(5)] * (r2_127 + _S6458 * u_132) + (*dist_coeffs_48)[int(6)] * r2_127, _S6459 * v_132 + (*dist_coeffs_48)[int(4)] * (r2_127 + _S6460 * v_132) + (*dist_coeffs_48)[int(7)] * r2_127);
    float2  _S6462 = _S6461 + make_float2 ((*dist_coeffs_48)[int(8)] * _S6461.x + (*dist_coeffs_48)[int(9)] * _S6461.y, 0.0f);
    float _S6463 = fx_41 * _S6462.x + cx_36;
    float _S6464 = fy_41 * _S6462.y + cy_36;
    float _S6465 = s_primal_ctx_max_0(_S6288, _S6313);
    float _S6466 = s_primal_ctx_min_0(_S6288, _S6313);
    float _S6467 = s_primal_ctx_max_0(_S6289, _S6314);
    float _S6468 = s_primal_ctx_min_0(_S6289, _S6314);
    float _S6469 = s_primal_ctx_max_0(_S6465, _S6338);
    float _S6470 = s_primal_ctx_min_0(_S6466, _S6338);
    float _S6471 = s_primal_ctx_max_0(_S6467, _S6339);
    float _S6472 = s_primal_ctx_min_0(_S6468, _S6339);
    float _S6473 = s_primal_ctx_max_0(_S6469, _S6363);
    float _S6474 = s_primal_ctx_min_0(_S6470, _S6363);
    float _S6475 = s_primal_ctx_max_0(_S6471, _S6364);
    float _S6476 = s_primal_ctx_min_0(_S6472, _S6364);
    float _S6477 = s_primal_ctx_max_0(_S6473, _S6388);
    float _S6478 = s_primal_ctx_min_0(_S6474, _S6388);
    float _S6479 = s_primal_ctx_max_0(_S6475, _S6389);
    float _S6480 = s_primal_ctx_min_0(_S6476, _S6389);
    float _S6481 = s_primal_ctx_max_0(_S6477, _S6413);
    float _S6482 = s_primal_ctx_min_0(_S6478, _S6413);
    float _S6483 = s_primal_ctx_max_0(_S6479, _S6414);
    float _S6484 = s_primal_ctx_min_0(_S6480, _S6414);
    float _S6485 = s_primal_ctx_max_0(_S6481, _S6438);
    float _S6486 = s_primal_ctx_min_0(_S6482, _S6438);
    float _S6487 = s_primal_ctx_max_0(_S6483, _S6439);
    float _S6488 = s_primal_ctx_min_0(_S6484, _S6439);
    Matrix<float, 3, 3>  _S6489 = transpose_0(R_36);
    float3  _S6490 = mean_c_31 - - s_primal_ctx_mul_1(_S6489, t_38);
    float _S6491 = _S6490.x;
    float _S6492 = _S6490.y;
    float _S6493 = _S6490.z;
    float _S6494 = _S6491 * _S6491 + _S6492 * _S6492 + _S6493 * _S6493;
    float _S6495 = s_primal_ctx_sqrt_0(_S6494);
    float x_68 = _S6491 / _S6495;
    float3  _S6496 = make_float3 (x_68);
    float _S6497 = _S6495 * _S6495;
    float y_35 = _S6492 / _S6495;
    float z_32 = _S6493 / _S6495;
    float3  _S6498 = make_float3 (z_32);
    float _S6499 = - y_35;
    float3  _S6500 = make_float3 (_S6499);
    float z2_65 = z_32 * z_32;
    float fTmp0B_32 = -1.09254848957061768f * z_32;
    float fC1_32 = x_68 * x_68 - y_35 * y_35;
    float _S6501 = 2.0f * x_68;
    float fS1_32 = _S6501 * y_35;
    float pSH6_10 = 0.94617468118667603f * z2_65 - 0.31539157032966614f;
    float3  _S6502 = make_float3 (pSH6_10);
    float pSH7_10 = fTmp0B_32 * x_68;
    float3  _S6503 = make_float3 (pSH7_10);
    float pSH5_10 = fTmp0B_32 * y_35;
    float3  _S6504 = make_float3 (pSH5_10);
    float pSH8_10 = 0.54627424478530884f * fC1_32;
    float3  _S6505 = make_float3 (pSH8_10);
    float pSH4_10 = 0.54627424478530884f * fS1_32;
    float3  _S6506 = make_float3 (pSH4_10);
    float fTmp0C_32 = -2.28522896766662598f * z2_65 + 0.4570457935333252f;
    float fTmp1B_32 = 1.44530570507049561f * z_32;
    float _S6507 = 1.86588168144226074f * z2_65 - 1.11952900886535645f;
    float pSH12_10 = z_32 * _S6507;
    float3  _S6508 = make_float3 (pSH12_10);
    float pSH13_10 = fTmp0C_32 * x_68;
    float3  _S6509 = make_float3 (pSH13_10);
    float pSH11_10 = fTmp0C_32 * y_35;
    float3  _S6510 = make_float3 (pSH11_10);
    float pSH14_10 = fTmp1B_32 * fC1_32;
    float3  _S6511 = make_float3 (pSH14_10);
    float pSH10_10 = fTmp1B_32 * fS1_32;
    float3  _S6512 = make_float3 (pSH10_10);
    float pSH15_10 = -0.59004360437393188f * (x_68 * fC1_32 - y_35 * fS1_32);
    float3  _S6513 = make_float3 (pSH15_10);
    float pSH9_10 = -0.59004360437393188f * (x_68 * fS1_32 + y_35 * fC1_32);
    float3  _S6514 = make_float3 (pSH9_10);
    float3  _S6515 = make_float3 (0.0f);
    float3  _S6516 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6517;
    (&_S6517)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_32)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S6499) * (*sh_coeffs_32)[int(1)] + make_float3 (z_32) * (*sh_coeffs_32)[int(2)] - make_float3 (x_68) * (*sh_coeffs_32)[int(3)]) + (make_float3 (pSH4_10) * (*sh_coeffs_32)[int(4)] + make_float3 (pSH5_10) * (*sh_coeffs_32)[int(5)] + make_float3 (pSH6_10) * (*sh_coeffs_32)[int(6)] + make_float3 (pSH7_10) * (*sh_coeffs_32)[int(7)] + make_float3 (pSH8_10) * (*sh_coeffs_32)[int(8)]) + (make_float3 (pSH9_10) * (*sh_coeffs_32)[int(9)] + make_float3 (pSH10_10) * (*sh_coeffs_32)[int(10)] + make_float3 (pSH11_10) * (*sh_coeffs_32)[int(11)] + make_float3 (pSH12_10) * (*sh_coeffs_32)[int(12)] + make_float3 (pSH13_10) * (*sh_coeffs_32)[int(13)] + make_float3 (pSH14_10) * (*sh_coeffs_32)[int(14)] + make_float3 (pSH15_10) * (*sh_coeffs_32)[int(15)]) + make_float3 (0.5f);
    (&_S6517)->differential_0 = _S6516;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6518;
    (&_S6518)->primal_0 = _S6515;
    (&_S6518)->differential_0 = _S6516;
    s_bwd_prop_max_0(&_S6517, &_S6518, v_rgb_9);
    float3  _S6519 = _S6513 * _S6517.differential_0;
    float3  _S6520 = (*sh_coeffs_32)[int(15)] * _S6517.differential_0;
    float3  _S6521 = _S6511 * _S6517.differential_0;
    float3  _S6522 = (*sh_coeffs_32)[int(14)] * _S6517.differential_0;
    float3  _S6523 = _S6509 * _S6517.differential_0;
    float3  _S6524 = (*sh_coeffs_32)[int(13)] * _S6517.differential_0;
    float3  _S6525 = _S6508 * _S6517.differential_0;
    float3  _S6526 = (*sh_coeffs_32)[int(12)] * _S6517.differential_0;
    float3  _S6527 = _S6510 * _S6517.differential_0;
    float3  _S6528 = (*sh_coeffs_32)[int(11)] * _S6517.differential_0;
    float3  _S6529 = _S6512 * _S6517.differential_0;
    float3  _S6530 = (*sh_coeffs_32)[int(10)] * _S6517.differential_0;
    float3  _S6531 = _S6514 * _S6517.differential_0;
    float3  _S6532 = (*sh_coeffs_32)[int(9)] * _S6517.differential_0;
    float s_diff_fS2_T_10 = -0.59004360437393188f * (_S6532.x + _S6532.y + _S6532.z);
    float s_diff_fC2_T_10 = -0.59004360437393188f * (_S6520.x + _S6520.y + _S6520.z);
    float _S6533 = _S6530.x + _S6530.y + _S6530.z;
    float _S6534 = _S6522.x + _S6522.y + _S6522.z;
    float _S6535 = _S6528.x + _S6528.y + _S6528.z;
    float _S6536 = _S6524.x + _S6524.y + _S6524.z;
    float _S6537 = _S6526.x + _S6526.y + _S6526.z;
    float _S6538 = - s_diff_fC2_T_10;
    float3  _S6539 = _S6505 * _S6517.differential_0;
    float3  _S6540 = (*sh_coeffs_32)[int(8)] * _S6517.differential_0;
    float3  _S6541 = _S6503 * _S6517.differential_0;
    float3  _S6542 = (*sh_coeffs_32)[int(7)] * _S6517.differential_0;
    float3  _S6543 = _S6502 * _S6517.differential_0;
    float3  _S6544 = (*sh_coeffs_32)[int(6)] * _S6517.differential_0;
    float3  _S6545 = _S6504 * _S6517.differential_0;
    float3  _S6546 = (*sh_coeffs_32)[int(5)] * _S6517.differential_0;
    float3  _S6547 = _S6506 * _S6517.differential_0;
    float3  _S6548 = (*sh_coeffs_32)[int(4)] * _S6517.differential_0;
    float _S6549 = _S6546.x + _S6546.y + _S6546.z;
    float _S6550 = _S6542.x + _S6542.y + _S6542.z;
    float _S6551 = fTmp1B_32 * _S6533 + x_68 * s_diff_fS2_T_10 + y_35 * _S6538 + 0.54627424478530884f * (_S6548.x + _S6548.y + _S6548.z);
    float _S6552 = fTmp1B_32 * _S6534 + y_35 * s_diff_fS2_T_10 + x_68 * s_diff_fC2_T_10 + 0.54627424478530884f * (_S6540.x + _S6540.y + _S6540.z);
    float _S6553 = y_35 * - _S6552;
    float _S6554 = x_68 * _S6552;
    float _S6555 = z_32 * (1.86588168144226074f * (z_32 * _S6537) + -2.28522896766662598f * (y_35 * _S6535 + x_68 * _S6536) + 0.94617468118667603f * (_S6544.x + _S6544.y + _S6544.z));
    float3  _S6556 = make_float3 (0.48860251903533936f) * _S6517.differential_0;
    float3  _S6557 = - _S6556;
    float3  _S6558 = _S6496 * _S6557;
    float3  _S6559 = (*sh_coeffs_32)[int(3)] * _S6557;
    float3  _S6560 = _S6498 * _S6556;
    float3  _S6561 = (*sh_coeffs_32)[int(2)] * _S6556;
    float3  _S6562 = _S6500 * _S6556;
    float3  _S6563 = (*sh_coeffs_32)[int(1)] * _S6556;
    float _S6564 = (_S6507 * _S6537 + 1.44530570507049561f * (fS1_32 * _S6533 + fC1_32 * _S6534) + -1.09254848957061768f * (y_35 * _S6549 + x_68 * _S6550) + _S6555 + _S6555 + _S6561.x + _S6561.y + _S6561.z) / _S6497;
    float _S6565 = _S6495 * _S6564;
    float _S6566 = (fTmp0C_32 * _S6535 + fC1_32 * s_diff_fS2_T_10 + fS1_32 * _S6538 + fTmp0B_32 * _S6549 + _S6501 * _S6551 + _S6553 + _S6553 + - (_S6563.x + _S6563.y + _S6563.z)) / _S6497;
    float _S6567 = _S6495 * _S6566;
    float _S6568 = (fTmp0C_32 * _S6536 + fS1_32 * s_diff_fS2_T_10 + fC1_32 * s_diff_fC2_T_10 + fTmp0B_32 * _S6550 + 2.0f * (y_35 * _S6551) + _S6554 + _S6554 + _S6559.x + _S6559.y + _S6559.z) / _S6497;
    float _S6569 = _S6495 * _S6568;
    float _S6570 = _S6493 * - _S6564 + _S6492 * - _S6566 + _S6491 * - _S6568;
    DiffPair_float_0 _S6571;
    (&_S6571)->primal_0 = _S6494;
    (&_S6571)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S6571, _S6570);
    float _S6572 = _S6493 * _S6571.differential_0;
    float _S6573 = _S6492 * _S6571.differential_0;
    float _S6574 = _S6491 * _S6571.differential_0;
    float3  _S6575 = make_float3 (0.282094806432724f) * _S6517.differential_0;
    float3  _S6576 = make_float3 (_S6569 + _S6574 + _S6574, _S6567 + _S6573 + _S6573, _S6565 + _S6572 + _S6572);
    float3  _S6577 = - - _S6576;
    Matrix<float, 3, 3>  _S6578 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6579;
    (&_S6579)->primal_0 = _S6489;
    (&_S6579)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6580;
    (&_S6580)->primal_0 = t_38;
    (&_S6580)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6579, &_S6580, _S6577);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6581 = _S6580;
    Matrix<float, 3, 3>  _S6582 = transpose_0(_S6579.differential_0);
    DiffPair_float_0 _S6583;
    (&_S6583)->primal_0 = _S6488;
    (&_S6583)->differential_0 = 0.0f;
    DiffPair_float_0 _S6584;
    (&_S6584)->primal_0 = _S6464;
    (&_S6584)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6583, &_S6584, 0.0f);
    DiffPair_float_0 _S6585;
    (&_S6585)->primal_0 = _S6487;
    (&_S6585)->differential_0 = 0.0f;
    DiffPair_float_0 _S6586;
    (&_S6586)->primal_0 = _S6464;
    (&_S6586)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6585, &_S6586, 0.0f);
    DiffPair_float_0 _S6587;
    (&_S6587)->primal_0 = _S6486;
    (&_S6587)->differential_0 = 0.0f;
    DiffPair_float_0 _S6588;
    (&_S6588)->primal_0 = _S6463;
    (&_S6588)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6587, &_S6588, 0.0f);
    DiffPair_float_0 _S6589;
    (&_S6589)->primal_0 = _S6485;
    (&_S6589)->differential_0 = 0.0f;
    DiffPair_float_0 _S6590;
    (&_S6590)->primal_0 = _S6463;
    (&_S6590)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6589, &_S6590, 0.0f);
    DiffPair_float_0 _S6591;
    (&_S6591)->primal_0 = _S6484;
    (&_S6591)->differential_0 = 0.0f;
    DiffPair_float_0 _S6592;
    (&_S6592)->primal_0 = _S6439;
    (&_S6592)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6591, &_S6592, _S6583.differential_0);
    DiffPair_float_0 _S6593;
    (&_S6593)->primal_0 = _S6483;
    (&_S6593)->differential_0 = 0.0f;
    DiffPair_float_0 _S6594;
    (&_S6594)->primal_0 = _S6439;
    (&_S6594)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6593, &_S6594, _S6585.differential_0);
    DiffPair_float_0 _S6595;
    (&_S6595)->primal_0 = _S6482;
    (&_S6595)->differential_0 = 0.0f;
    DiffPair_float_0 _S6596;
    (&_S6596)->primal_0 = _S6438;
    (&_S6596)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6595, &_S6596, _S6587.differential_0);
    DiffPair_float_0 _S6597;
    (&_S6597)->primal_0 = _S6481;
    (&_S6597)->differential_0 = 0.0f;
    DiffPair_float_0 _S6598;
    (&_S6598)->primal_0 = _S6438;
    (&_S6598)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6597, &_S6598, _S6589.differential_0);
    DiffPair_float_0 _S6599;
    (&_S6599)->primal_0 = _S6480;
    (&_S6599)->differential_0 = 0.0f;
    DiffPair_float_0 _S6600;
    (&_S6600)->primal_0 = _S6414;
    (&_S6600)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6599, &_S6600, _S6591.differential_0);
    DiffPair_float_0 _S6601;
    (&_S6601)->primal_0 = _S6479;
    (&_S6601)->differential_0 = 0.0f;
    DiffPair_float_0 _S6602;
    (&_S6602)->primal_0 = _S6414;
    (&_S6602)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6601, &_S6602, _S6593.differential_0);
    DiffPair_float_0 _S6603;
    (&_S6603)->primal_0 = _S6478;
    (&_S6603)->differential_0 = 0.0f;
    DiffPair_float_0 _S6604;
    (&_S6604)->primal_0 = _S6413;
    (&_S6604)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6603, &_S6604, _S6595.differential_0);
    DiffPair_float_0 _S6605;
    (&_S6605)->primal_0 = _S6477;
    (&_S6605)->differential_0 = 0.0f;
    DiffPair_float_0 _S6606;
    (&_S6606)->primal_0 = _S6413;
    (&_S6606)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6605, &_S6606, _S6597.differential_0);
    DiffPair_float_0 _S6607;
    (&_S6607)->primal_0 = _S6476;
    (&_S6607)->differential_0 = 0.0f;
    DiffPair_float_0 _S6608;
    (&_S6608)->primal_0 = _S6389;
    (&_S6608)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6607, &_S6608, _S6599.differential_0);
    DiffPair_float_0 _S6609;
    (&_S6609)->primal_0 = _S6475;
    (&_S6609)->differential_0 = 0.0f;
    DiffPair_float_0 _S6610;
    (&_S6610)->primal_0 = _S6389;
    (&_S6610)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6609, &_S6610, _S6601.differential_0);
    DiffPair_float_0 _S6611;
    (&_S6611)->primal_0 = _S6474;
    (&_S6611)->differential_0 = 0.0f;
    DiffPair_float_0 _S6612;
    (&_S6612)->primal_0 = _S6388;
    (&_S6612)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6611, &_S6612, _S6603.differential_0);
    DiffPair_float_0 _S6613;
    (&_S6613)->primal_0 = _S6473;
    (&_S6613)->differential_0 = 0.0f;
    DiffPair_float_0 _S6614;
    (&_S6614)->primal_0 = _S6388;
    (&_S6614)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6613, &_S6614, _S6605.differential_0);
    DiffPair_float_0 _S6615;
    (&_S6615)->primal_0 = _S6472;
    (&_S6615)->differential_0 = 0.0f;
    DiffPair_float_0 _S6616;
    (&_S6616)->primal_0 = _S6364;
    (&_S6616)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6615, &_S6616, _S6607.differential_0);
    DiffPair_float_0 _S6617;
    (&_S6617)->primal_0 = _S6471;
    (&_S6617)->differential_0 = 0.0f;
    DiffPair_float_0 _S6618;
    (&_S6618)->primal_0 = _S6364;
    (&_S6618)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6617, &_S6618, _S6609.differential_0);
    DiffPair_float_0 _S6619;
    (&_S6619)->primal_0 = _S6470;
    (&_S6619)->differential_0 = 0.0f;
    DiffPair_float_0 _S6620;
    (&_S6620)->primal_0 = _S6363;
    (&_S6620)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6619, &_S6620, _S6611.differential_0);
    DiffPair_float_0 _S6621;
    (&_S6621)->primal_0 = _S6469;
    (&_S6621)->differential_0 = 0.0f;
    DiffPair_float_0 _S6622;
    (&_S6622)->primal_0 = _S6363;
    (&_S6622)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6621, &_S6622, _S6613.differential_0);
    DiffPair_float_0 _S6623;
    (&_S6623)->primal_0 = _S6468;
    (&_S6623)->differential_0 = 0.0f;
    DiffPair_float_0 _S6624;
    (&_S6624)->primal_0 = _S6339;
    (&_S6624)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6623, &_S6624, _S6615.differential_0);
    DiffPair_float_0 _S6625;
    (&_S6625)->primal_0 = _S6467;
    (&_S6625)->differential_0 = 0.0f;
    DiffPair_float_0 _S6626;
    (&_S6626)->primal_0 = _S6339;
    (&_S6626)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6625, &_S6626, _S6617.differential_0);
    DiffPair_float_0 _S6627;
    (&_S6627)->primal_0 = _S6466;
    (&_S6627)->differential_0 = 0.0f;
    DiffPair_float_0 _S6628;
    (&_S6628)->primal_0 = _S6338;
    (&_S6628)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6627, &_S6628, _S6619.differential_0);
    DiffPair_float_0 _S6629;
    (&_S6629)->primal_0 = _S6465;
    (&_S6629)->differential_0 = 0.0f;
    DiffPair_float_0 _S6630;
    (&_S6630)->primal_0 = _S6338;
    (&_S6630)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6629, &_S6630, _S6621.differential_0);
    DiffPair_float_0 _S6631;
    (&_S6631)->primal_0 = _S6289;
    (&_S6631)->differential_0 = 0.0f;
    DiffPair_float_0 _S6632;
    (&_S6632)->primal_0 = _S6314;
    (&_S6632)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6631, &_S6632, _S6623.differential_0);
    DiffPair_float_0 _S6633;
    (&_S6633)->primal_0 = _S6289;
    (&_S6633)->differential_0 = 0.0f;
    DiffPair_float_0 _S6634;
    (&_S6634)->primal_0 = _S6314;
    (&_S6634)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6633, &_S6634, _S6625.differential_0);
    DiffPair_float_0 _S6635;
    (&_S6635)->primal_0 = _S6288;
    (&_S6635)->differential_0 = 0.0f;
    DiffPair_float_0 _S6636;
    (&_S6636)->primal_0 = _S6313;
    (&_S6636)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6635, &_S6636, _S6627.differential_0);
    DiffPair_float_0 _S6637;
    (&_S6637)->primal_0 = _S6288;
    (&_S6637)->differential_0 = 0.0f;
    DiffPair_float_0 _S6638;
    (&_S6638)->primal_0 = _S6313;
    (&_S6638)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6637, &_S6638, _S6629.differential_0);
    float _S6639 = fx_41 * (_S6588.differential_0 + _S6590.differential_0);
    float2  _S6640 = make_float2 (_S6639, fy_41 * (_S6584.differential_0 + _S6586.differential_0)) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6639, (*dist_coeffs_48)[int(9)] * _S6639);
    float2  _S6641 = _S6452 * _S6640;
    float2  _S6642 = _S6456 * _S6640;
    float _S6643 = (*dist_coeffs_48)[int(4)] * _S6640.y;
    float _S6644 = (*dist_coeffs_48)[int(5)] * _S6640.x;
    float _S6645 = _S6641.x + _S6641.y;
    float _S6646 = r2_127 * _S6645;
    float _S6647 = r2_127 * _S6646;
    float _S6648 = (*dist_coeffs_48)[int(7)] * _S6640.y + _S6643 + (*dist_coeffs_48)[int(6)] * _S6640.x + _S6644 + _S6455 * _S6645 + _S6454 * _S6646 + _S6453 * _S6647 + (*dist_coeffs_48)[int(3)] * (r2_127 * _S6647);
    float _S6649 = v_132 * _S6648;
    float _S6650 = u_132 * _S6648;
    float _S6651 = _S6460 * _S6643 + 2.0f * (v_132 * _S6643) + _S6459 * _S6640.y + _S6457 * _S6640.x + _S6649 + _S6649;
    float _S6652 = _S6283 * (v_132 * _S6640.y) + _S6458 * _S6644 + 2.0f * (u_132 * _S6644) + _S6280 * (v_132 * _S6640.x) + _S6650 + _S6650;
    FixedArray<float3 , 16>  _S6653;
    _S6653[int(0)] = _S6516;
    _S6653[int(1)] = _S6516;
    _S6653[int(2)] = _S6516;
    _S6653[int(3)] = _S6516;
    _S6653[int(4)] = _S6516;
    _S6653[int(5)] = _S6516;
    _S6653[int(6)] = _S6516;
    _S6653[int(7)] = _S6516;
    _S6653[int(8)] = _S6516;
    _S6653[int(9)] = _S6516;
    _S6653[int(10)] = _S6516;
    _S6653[int(11)] = _S6516;
    _S6653[int(12)] = _S6516;
    _S6653[int(13)] = _S6516;
    _S6653[int(14)] = _S6516;
    _S6653[int(15)] = _S6516;
    _S6653[int(7)] = _S6541;
    _S6653[int(0)] = _S6575;
    _S6653[int(1)] = _S6562;
    _S6653[int(2)] = _S6560;
    _S6653[int(3)] = _S6558;
    _S6653[int(4)] = _S6547;
    _S6653[int(5)] = _S6545;
    _S6653[int(6)] = _S6543;
    _S6653[int(15)] = _S6519;
    _S6653[int(8)] = _S6539;
    _S6653[int(9)] = _S6531;
    _S6653[int(10)] = _S6529;
    _S6653[int(11)] = _S6527;
    _S6653[int(12)] = _S6525;
    _S6653[int(13)] = _S6523;
    _S6653[int(14)] = _S6521;
    float3  _S6654 = _S6653[int(0)];
    float3  _S6655 = _S6653[int(1)];
    float3  _S6656 = _S6653[int(2)];
    float3  _S6657 = _S6653[int(3)];
    float3  _S6658 = _S6653[int(4)];
    float3  _S6659 = _S6653[int(5)];
    float3  _S6660 = _S6653[int(6)];
    float3  _S6661 = _S6653[int(7)];
    float3  _S6662 = _S6653[int(8)];
    float3  _S6663 = _S6653[int(9)];
    float3  _S6664 = _S6653[int(10)];
    float3  _S6665 = _S6653[int(11)];
    float3  _S6666 = _S6653[int(12)];
    float3  _S6667 = _S6653[int(13)];
    float3  _S6668 = _S6653[int(14)];
    float3  _S6669 = _S6653[int(15)];
    float _S6670 = _S6636.differential_0 + _S6638.differential_0;
    float _S6671 = _S6631.differential_0 + _S6633.differential_0;
    float _S6672 = _S6632.differential_0 + _S6634.differential_0;
    float _S6673 = _S6628.differential_0 + _S6630.differential_0;
    float _S6674 = _S6635.differential_0 + _S6637.differential_0;
    float _S6675 = _S6592.differential_0 + _S6594.differential_0;
    float _S6676 = _S6596.differential_0 + _S6598.differential_0;
    float _S6677 = _S6600.differential_0 + _S6602.differential_0;
    float _S6678 = _S6604.differential_0 + _S6606.differential_0;
    float _S6679 = _S6608.differential_0 + _S6610.differential_0;
    float _S6680 = _S6612.differential_0 + _S6614.differential_0;
    float _S6681 = _S6616.differential_0 + _S6618.differential_0;
    float _S6682 = _S6620.differential_0 + _S6622.differential_0;
    float _S6683 = _S6624.differential_0 + _S6626.differential_0;
    float2  _S6684 = _S6642 + make_float2 (_S6652, _S6651);
    float2  _S6685 = _S6440 * _S6684;
    float2  _S6686 = _S6451 * _S6684;
    float _S6687 = _S6685.x + _S6685.y;
    if(_S6444)
    {
        float _S6688 = _S6687 / _S6446;
        float _S6689 = _S6447 * - _S6688;
        float _S6690 = _S6443 * (0.3333333432674408f * - (_S6442 * _S6688));
        k_15 = _S6690 + _S6690;
        _S6445 = _S6689;
        _S6446 = 0.0f;
    }
    else
    {
        float _S6691 = _S6687 / _S6445;
        float _S6692 = _S6443 * - _S6691;
        k_15 = _S6441 * _S6691;
        _S6445 = 0.0f;
        _S6446 = _S6692;
    }
    DiffPair_float_0 _S6693;
    (&_S6693)->primal_0 = _S6441;
    (&_S6693)->differential_0 = 0.0f;
    DiffPair_float_0 _S6694;
    (&_S6694)->primal_0 = _S6442;
    (&_S6694)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6693, &_S6694, k_15);
    float _S6695 = _S6694.differential_0 + _S6445;
    float _S6696 = _S6693.differential_0 + _S6446;
    float2  _S6697 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6698;
    (&_S6698)->primal_0 = _S6440;
    (&_S6698)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6698, _S6696);
    float2  _S6699 = _S6698.differential_0 + _S6686;
    float _S6700 = fx_41 * _S6676;
    float2  _S6701 = make_float2 (_S6700, fy_41 * _S6675) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6700, (*dist_coeffs_48)[int(9)] * _S6700);
    float2  _S6702 = _S6427 * _S6701;
    float _S6703 = (*dist_coeffs_48)[int(4)] * _S6701.y;
    float _S6704 = (*dist_coeffs_48)[int(5)] * _S6701.x;
    float _S6705 = _S6702.x + _S6702.y;
    float _S6706 = r2_126 * _S6705;
    float _S6707 = r2_126 * _S6706;
    float _S6708 = (*dist_coeffs_48)[int(7)] * _S6701.y + _S6703 + (*dist_coeffs_48)[int(6)] * _S6701.x + _S6704 + _S6430 * _S6705 + _S6429 * _S6706 + _S6428 * _S6707 + (*dist_coeffs_48)[int(3)] * (r2_126 * _S6707);
    float _S6709 = v_131 * _S6708;
    float _S6710 = u_131 * _S6708;
    float3  _S6711 = make_float3 (_S6699.x, _S6699.y, _S6695);
    float2  _S6712 = _S6431 * _S6701 + make_float2 (_S6283 * (v_131 * _S6701.y) + _S6433 * _S6704 + 2.0f * (u_131 * _S6704) + _S6280 * (v_131 * _S6701.x) + _S6710 + _S6710, _S6435 * _S6703 + 2.0f * (v_131 * _S6703) + _S6434 * _S6701.y + _S6432 * _S6701.x + _S6709 + _S6709);
    float2  _S6713 = _S6415 * _S6712;
    float2  _S6714 = _S6426 * _S6712;
    float _S6715 = _S6713.x + _S6713.y;
    if(_S6419)
    {
        float _S6716 = _S6715 / _S6421;
        float _S6717 = _S6422 * - _S6716;
        float _S6718 = _S6418 * (0.3333333432674408f * - (_S6417 * _S6716));
        k_15 = _S6718 + _S6718;
        _S6420 = _S6717;
        _S6421 = 0.0f;
    }
    else
    {
        float _S6719 = _S6715 / _S6420;
        float _S6720 = _S6418 * - _S6719;
        k_15 = _S6416 * _S6719;
        _S6420 = 0.0f;
        _S6421 = _S6720;
    }
    DiffPair_float_0 _S6721;
    (&_S6721)->primal_0 = _S6416;
    (&_S6721)->differential_0 = 0.0f;
    DiffPair_float_0 _S6722;
    (&_S6722)->primal_0 = _S6417;
    (&_S6722)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6721, &_S6722, k_15);
    float _S6723 = _S6722.differential_0 + _S6420;
    float _S6724 = _S6721.differential_0 + _S6421;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6725;
    (&_S6725)->primal_0 = _S6415;
    (&_S6725)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6725, _S6724);
    float2  _S6726 = _S6725.differential_0 + _S6714;
    float _S6727 = fx_41 * _S6678;
    float2  _S6728 = make_float2 (_S6727, fy_41 * _S6677) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6727, (*dist_coeffs_48)[int(9)] * _S6727);
    float2  _S6729 = _S6402 * _S6728;
    float _S6730 = (*dist_coeffs_48)[int(4)] * _S6728.y;
    float _S6731 = (*dist_coeffs_48)[int(5)] * _S6728.x;
    float _S6732 = _S6729.x + _S6729.y;
    float _S6733 = r2_125 * _S6732;
    float _S6734 = r2_125 * _S6733;
    float _S6735 = (*dist_coeffs_48)[int(7)] * _S6728.y + _S6730 + (*dist_coeffs_48)[int(6)] * _S6728.x + _S6731 + _S6405 * _S6732 + _S6404 * _S6733 + _S6403 * _S6734 + (*dist_coeffs_48)[int(3)] * (r2_125 * _S6734);
    float _S6736 = v_130 * _S6735;
    float _S6737 = u_130 * _S6735;
    float3  _S6738 = make_float3 (_S6726.x, _S6726.y, _S6723);
    float2  _S6739 = _S6406 * _S6728 + make_float2 (_S6283 * (v_130 * _S6728.y) + _S6408 * _S6731 + 2.0f * (u_130 * _S6731) + _S6280 * (v_130 * _S6728.x) + _S6737 + _S6737, _S6410 * _S6730 + 2.0f * (v_130 * _S6730) + _S6409 * _S6728.y + _S6407 * _S6728.x + _S6736 + _S6736);
    float2  _S6740 = _S6390 * _S6739;
    float2  _S6741 = _S6401 * _S6739;
    float _S6742 = _S6740.x + _S6740.y;
    if(_S6394)
    {
        float _S6743 = _S6742 / _S6396;
        float _S6744 = _S6397 * - _S6743;
        float _S6745 = _S6393 * (0.3333333432674408f * - (_S6392 * _S6743));
        k_15 = _S6745 + _S6745;
        _S6395 = _S6744;
        _S6396 = 0.0f;
    }
    else
    {
        float _S6746 = _S6742 / _S6395;
        float _S6747 = _S6393 * - _S6746;
        k_15 = _S6391 * _S6746;
        _S6395 = 0.0f;
        _S6396 = _S6747;
    }
    DiffPair_float_0 _S6748;
    (&_S6748)->primal_0 = _S6391;
    (&_S6748)->differential_0 = 0.0f;
    DiffPair_float_0 _S6749;
    (&_S6749)->primal_0 = _S6392;
    (&_S6749)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6748, &_S6749, k_15);
    float _S6750 = _S6749.differential_0 + _S6395;
    float _S6751 = _S6748.differential_0 + _S6396;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6752;
    (&_S6752)->primal_0 = _S6390;
    (&_S6752)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6752, _S6751);
    float2  _S6753 = _S6752.differential_0 + _S6741;
    float _S6754 = fx_41 * _S6680;
    float2  _S6755 = make_float2 (_S6754, fy_41 * _S6679) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6754, (*dist_coeffs_48)[int(9)] * _S6754);
    float2  _S6756 = _S6377 * _S6755;
    float _S6757 = (*dist_coeffs_48)[int(4)] * _S6755.y;
    float _S6758 = (*dist_coeffs_48)[int(5)] * _S6755.x;
    float _S6759 = _S6756.x + _S6756.y;
    float _S6760 = r2_124 * _S6759;
    float _S6761 = r2_124 * _S6760;
    float _S6762 = (*dist_coeffs_48)[int(7)] * _S6755.y + _S6757 + (*dist_coeffs_48)[int(6)] * _S6755.x + _S6758 + _S6380 * _S6759 + _S6379 * _S6760 + _S6378 * _S6761 + (*dist_coeffs_48)[int(3)] * (r2_124 * _S6761);
    float _S6763 = v_129 * _S6762;
    float _S6764 = u_129 * _S6762;
    float3  _S6765 = make_float3 (_S6753.x, _S6753.y, _S6750);
    float2  _S6766 = _S6381 * _S6755 + make_float2 (_S6283 * (v_129 * _S6755.y) + _S6383 * _S6758 + 2.0f * (u_129 * _S6758) + _S6280 * (v_129 * _S6755.x) + _S6764 + _S6764, _S6385 * _S6757 + 2.0f * (v_129 * _S6757) + _S6384 * _S6755.y + _S6382 * _S6755.x + _S6763 + _S6763);
    float2  _S6767 = _S6365 * _S6766;
    float2  _S6768 = _S6376 * _S6766;
    float _S6769 = _S6767.x + _S6767.y;
    if(_S6369)
    {
        float _S6770 = _S6769 / _S6371;
        float _S6771 = _S6372 * - _S6770;
        float _S6772 = _S6368 * (0.3333333432674408f * - (_S6367 * _S6770));
        k_15 = _S6772 + _S6772;
        _S6370 = _S6771;
        _S6371 = 0.0f;
    }
    else
    {
        float _S6773 = _S6769 / _S6370;
        float _S6774 = _S6368 * - _S6773;
        k_15 = _S6366 * _S6773;
        _S6370 = 0.0f;
        _S6371 = _S6774;
    }
    DiffPair_float_0 _S6775;
    (&_S6775)->primal_0 = _S6366;
    (&_S6775)->differential_0 = 0.0f;
    DiffPair_float_0 _S6776;
    (&_S6776)->primal_0 = _S6367;
    (&_S6776)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6775, &_S6776, k_15);
    float _S6777 = _S6776.differential_0 + _S6370;
    float _S6778 = _S6775.differential_0 + _S6371;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6779;
    (&_S6779)->primal_0 = _S6365;
    (&_S6779)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6779, _S6778);
    float2  _S6780 = _S6779.differential_0 + _S6768;
    float _S6781 = fx_41 * _S6682;
    float2  _S6782 = make_float2 (_S6781, fy_41 * _S6681) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6781, (*dist_coeffs_48)[int(9)] * _S6781);
    float2  _S6783 = _S6352 * _S6782;
    float _S6784 = (*dist_coeffs_48)[int(4)] * _S6782.y;
    float _S6785 = (*dist_coeffs_48)[int(5)] * _S6782.x;
    float _S6786 = _S6783.x + _S6783.y;
    float _S6787 = r2_123 * _S6786;
    float _S6788 = r2_123 * _S6787;
    float _S6789 = (*dist_coeffs_48)[int(7)] * _S6782.y + _S6784 + (*dist_coeffs_48)[int(6)] * _S6782.x + _S6785 + _S6355 * _S6786 + _S6354 * _S6787 + _S6353 * _S6788 + (*dist_coeffs_48)[int(3)] * (r2_123 * _S6788);
    float _S6790 = v_128 * _S6789;
    float _S6791 = u_128 * _S6789;
    float3  _S6792 = make_float3 (_S6780.x, _S6780.y, _S6777);
    float2  _S6793 = _S6356 * _S6782 + make_float2 (_S6283 * (v_128 * _S6782.y) + _S6358 * _S6785 + 2.0f * (u_128 * _S6785) + _S6280 * (v_128 * _S6782.x) + _S6791 + _S6791, _S6360 * _S6784 + 2.0f * (v_128 * _S6784) + _S6359 * _S6782.y + _S6357 * _S6782.x + _S6790 + _S6790);
    float2  _S6794 = _S6340 * _S6793;
    float2  _S6795 = _S6351 * _S6793;
    float _S6796 = _S6794.x + _S6794.y;
    if(_S6344)
    {
        float _S6797 = _S6796 / _S6346;
        float _S6798 = _S6347 * - _S6797;
        float _S6799 = _S6343 * (0.3333333432674408f * - (_S6342 * _S6797));
        k_15 = _S6799 + _S6799;
        _S6345 = _S6798;
        _S6346 = 0.0f;
    }
    else
    {
        float _S6800 = _S6796 / _S6345;
        float _S6801 = _S6343 * - _S6800;
        k_15 = _S6341 * _S6800;
        _S6345 = 0.0f;
        _S6346 = _S6801;
    }
    DiffPair_float_0 _S6802;
    (&_S6802)->primal_0 = _S6341;
    (&_S6802)->differential_0 = 0.0f;
    DiffPair_float_0 _S6803;
    (&_S6803)->primal_0 = _S6342;
    (&_S6803)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6802, &_S6803, k_15);
    float _S6804 = _S6803.differential_0 + _S6345;
    float _S6805 = _S6802.differential_0 + _S6346;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6806;
    (&_S6806)->primal_0 = _S6340;
    (&_S6806)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6806, _S6805);
    float2  _S6807 = _S6806.differential_0 + _S6795;
    float _S6808 = fx_41 * _S6673;
    float2  _S6809 = make_float2 (_S6808, fy_41 * _S6683) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6808, (*dist_coeffs_48)[int(9)] * _S6808);
    float2  _S6810 = _S6327 * _S6809;
    float _S6811 = (*dist_coeffs_48)[int(4)] * _S6809.y;
    float _S6812 = (*dist_coeffs_48)[int(5)] * _S6809.x;
    float _S6813 = _S6810.x + _S6810.y;
    float _S6814 = r2_122 * _S6813;
    float _S6815 = r2_122 * _S6814;
    float _S6816 = (*dist_coeffs_48)[int(7)] * _S6809.y + _S6811 + (*dist_coeffs_48)[int(6)] * _S6809.x + _S6812 + _S6330 * _S6813 + _S6329 * _S6814 + _S6328 * _S6815 + (*dist_coeffs_48)[int(3)] * (r2_122 * _S6815);
    float _S6817 = v_127 * _S6816;
    float _S6818 = u_127 * _S6816;
    float3  _S6819 = make_float3 (_S6807.x, _S6807.y, _S6804);
    float2  _S6820 = _S6331 * _S6809 + make_float2 (_S6283 * (v_127 * _S6809.y) + _S6333 * _S6812 + 2.0f * (u_127 * _S6812) + _S6280 * (v_127 * _S6809.x) + _S6818 + _S6818, _S6335 * _S6811 + 2.0f * (v_127 * _S6811) + _S6334 * _S6809.y + _S6332 * _S6809.x + _S6817 + _S6817);
    float2  _S6821 = _S6315 * _S6820;
    float2  _S6822 = _S6326 * _S6820;
    float _S6823 = _S6821.x + _S6821.y;
    if(_S6319)
    {
        float _S6824 = _S6823 / _S6321;
        float _S6825 = _S6322 * - _S6824;
        float _S6826 = _S6318 * (0.3333333432674408f * - (_S6317 * _S6824));
        k_15 = _S6826 + _S6826;
        _S6320 = _S6825;
        _S6321 = 0.0f;
    }
    else
    {
        float _S6827 = _S6823 / _S6320;
        float _S6828 = _S6318 * - _S6827;
        k_15 = _S6316 * _S6827;
        _S6320 = 0.0f;
        _S6321 = _S6828;
    }
    DiffPair_float_0 _S6829;
    (&_S6829)->primal_0 = _S6316;
    (&_S6829)->differential_0 = 0.0f;
    DiffPair_float_0 _S6830;
    (&_S6830)->primal_0 = _S6317;
    (&_S6830)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6829, &_S6830, k_15);
    float _S6831 = _S6830.differential_0 + _S6320;
    float _S6832 = _S6829.differential_0 + _S6321;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6833;
    (&_S6833)->primal_0 = _S6315;
    (&_S6833)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6833, _S6832);
    float2  _S6834 = _S6833.differential_0 + _S6822;
    float _S6835 = fx_41 * _S6670;
    float2  _S6836 = make_float2 (_S6835, fy_41 * _S6672) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6835, (*dist_coeffs_48)[int(9)] * _S6835);
    float2  _S6837 = _S6302 * _S6836;
    float _S6838 = (*dist_coeffs_48)[int(4)] * _S6836.y;
    float _S6839 = (*dist_coeffs_48)[int(5)] * _S6836.x;
    float _S6840 = _S6837.x + _S6837.y;
    float _S6841 = r2_121 * _S6840;
    float _S6842 = r2_121 * _S6841;
    float _S6843 = (*dist_coeffs_48)[int(7)] * _S6836.y + _S6838 + (*dist_coeffs_48)[int(6)] * _S6836.x + _S6839 + _S6305 * _S6840 + _S6304 * _S6841 + _S6303 * _S6842 + (*dist_coeffs_48)[int(3)] * (r2_121 * _S6842);
    float _S6844 = v_126 * _S6843;
    float _S6845 = u_126 * _S6843;
    float3  _S6846 = make_float3 (_S6834.x, _S6834.y, _S6831);
    float2  _S6847 = _S6306 * _S6836 + make_float2 (_S6283 * (v_126 * _S6836.y) + _S6308 * _S6839 + 2.0f * (u_126 * _S6839) + _S6280 * (v_126 * _S6836.x) + _S6845 + _S6845, _S6310 * _S6838 + 2.0f * (v_126 * _S6838) + _S6309 * _S6836.y + _S6307 * _S6836.x + _S6844 + _S6844);
    float2  _S6848 = _S6290 * _S6847;
    float2  _S6849 = _S6301 * _S6847;
    float _S6850 = _S6848.x + _S6848.y;
    if(_S6294)
    {
        float _S6851 = _S6850 / _S6296;
        float _S6852 = _S6297 * - _S6851;
        float _S6853 = _S6293 * (0.3333333432674408f * - (_S6292 * _S6851));
        k_15 = _S6853 + _S6853;
        _S6295 = _S6852;
        _S6296 = 0.0f;
    }
    else
    {
        float _S6854 = _S6850 / _S6295;
        float _S6855 = _S6293 * - _S6854;
        k_15 = _S6291 * _S6854;
        _S6295 = 0.0f;
        _S6296 = _S6855;
    }
    DiffPair_float_0 _S6856;
    (&_S6856)->primal_0 = _S6291;
    (&_S6856)->differential_0 = 0.0f;
    DiffPair_float_0 _S6857;
    (&_S6857)->primal_0 = _S6292;
    (&_S6857)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6856, &_S6857, k_15);
    float _S6858 = _S6857.differential_0 + _S6295;
    float _S6859 = _S6856.differential_0 + _S6296;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6860;
    (&_S6860)->primal_0 = _S6290;
    (&_S6860)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6860, _S6859);
    float2  _S6861 = _S6860.differential_0 + _S6849;
    float _S6862 = fx_41 * _S6674;
    float2  _S6863 = make_float2 (_S6862, fy_41 * _S6671) + make_float2 ((*dist_coeffs_48)[int(8)] * _S6862, (*dist_coeffs_48)[int(9)] * _S6862);
    float2  _S6864 = _S6275 * _S6863;
    float _S6865 = (*dist_coeffs_48)[int(4)] * _S6863.y;
    float _S6866 = (*dist_coeffs_48)[int(5)] * _S6863.x;
    float _S6867 = _S6864.x + _S6864.y;
    float _S6868 = r2_120 * _S6867;
    float _S6869 = r2_120 * _S6868;
    float _S6870 = (*dist_coeffs_48)[int(7)] * _S6863.y + _S6865 + (*dist_coeffs_48)[int(6)] * _S6863.x + _S6866 + _S6278 * _S6867 + _S6277 * _S6868 + _S6276 * _S6869 + (*dist_coeffs_48)[int(3)] * (r2_120 * _S6869);
    float _S6871 = v_125 * _S6870;
    float _S6872 = u_125 * _S6870;
    float3  _S6873 = make_float3 (_S6861.x, _S6861.y, _S6858);
    float2  _S6874 = _S6279 * _S6863 + make_float2 (_S6283 * (v_125 * _S6863.y) + _S6282 * _S6866 + 2.0f * (u_125 * _S6866) + _S6280 * (v_125 * _S6863.x) + _S6872 + _S6872, _S6285 * _S6865 + 2.0f * (v_125 * _S6865) + _S6284 * _S6863.y + _S6281 * _S6863.x + _S6871 + _S6871);
    float2  _S6875 = _S6263 * _S6874;
    float2  _S6876 = _S6274 * _S6874;
    float _S6877 = _S6875.x + _S6875.y;
    if(_S6267)
    {
        float _S6878 = _S6877 / _S6269;
        float _S6879 = _S6270 * - _S6878;
        float _S6880 = _S6266 * (0.3333333432674408f * - (_S6265 * _S6878));
        k_15 = _S6880 + _S6880;
        _S6268 = _S6879;
        _S6269 = 0.0f;
    }
    else
    {
        float _S6881 = _S6877 / _S6268;
        float _S6882 = _S6266 * - _S6881;
        k_15 = _S6264 * _S6881;
        _S6268 = 0.0f;
        _S6269 = _S6882;
    }
    DiffPair_float_0 _S6883;
    (&_S6883)->primal_0 = _S6264;
    (&_S6883)->differential_0 = 0.0f;
    DiffPair_float_0 _S6884;
    (&_S6884)->primal_0 = _S6265;
    (&_S6884)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6883, &_S6884, k_15);
    float _S6885 = _S6884.differential_0 + _S6268;
    float _S6886 = _S6883.differential_0 + _S6269;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6887;
    (&_S6887)->primal_0 = _S6263;
    (&_S6887)->differential_0 = _S6697;
    s_bwd_length_impl_1(&_S6887, _S6886);
    float2  _S6888 = _S6887.differential_0 + _S6876;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6889;
    (&_S6889)->primal_0 = R_36;
    (&_S6889)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6890;
    (&_S6890)->primal_0 = _S6262;
    (&_S6890)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6889, &_S6890, _S6576);
    DiffPair_float_0 _S6891;
    (&_S6891)->primal_0 = _S6259;
    (&_S6891)->differential_0 = 0.0f;
    DiffPair_float_0 _S6892;
    (&_S6892)->primal_0 = _S6261;
    (&_S6892)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6891, &_S6892, 0.0f);
    DiffPair_float_0 _S6893;
    (&_S6893)->primal_0 = _S6258;
    (&_S6893)->differential_0 = 0.0f;
    DiffPair_float_0 _S6894;
    (&_S6894)->primal_0 = _S6261;
    (&_S6894)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6893, &_S6894, 0.0f);
    float _S6895 = _S6892.differential_0 + _S6894.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6896;
    (&_S6896)->primal_0 = _S6260;
    (&_S6896)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6896, _S6895);
    float3  _S6897 = _S6896.differential_0 + _S6711;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6898;
    (&_S6898)->primal_0 = R_36;
    (&_S6898)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6899;
    (&_S6899)->primal_0 = pos_i_13;
    (&_S6899)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6898, &_S6899, _S6897);
    DiffPair_float_0 _S6900;
    (&_S6900)->primal_0 = _S6255;
    (&_S6900)->differential_0 = 0.0f;
    DiffPair_float_0 _S6901;
    (&_S6901)->primal_0 = _S6257;
    (&_S6901)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6900, &_S6901, _S6891.differential_0);
    DiffPair_float_0 _S6902;
    (&_S6902)->primal_0 = _S6254;
    (&_S6902)->differential_0 = 0.0f;
    DiffPair_float_0 _S6903;
    (&_S6903)->primal_0 = _S6257;
    (&_S6903)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6902, &_S6903, _S6893.differential_0);
    float _S6904 = _S6901.differential_0 + _S6903.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6905;
    (&_S6905)->primal_0 = _S6256;
    (&_S6905)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6905, _S6904);
    float3  _S6906 = _S6905.differential_0 + _S6738;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6907;
    (&_S6907)->primal_0 = R_36;
    (&_S6907)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6908;
    (&_S6908)->primal_0 = pos_i_12;
    (&_S6908)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6907, &_S6908, _S6906);
    DiffPair_float_0 _S6909;
    (&_S6909)->primal_0 = _S6251;
    (&_S6909)->differential_0 = 0.0f;
    DiffPair_float_0 _S6910;
    (&_S6910)->primal_0 = _S6253;
    (&_S6910)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6909, &_S6910, _S6900.differential_0);
    DiffPair_float_0 _S6911;
    (&_S6911)->primal_0 = _S6250;
    (&_S6911)->differential_0 = 0.0f;
    DiffPair_float_0 _S6912;
    (&_S6912)->primal_0 = _S6253;
    (&_S6912)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6911, &_S6912, _S6902.differential_0);
    float _S6913 = _S6910.differential_0 + _S6912.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6914;
    (&_S6914)->primal_0 = _S6252;
    (&_S6914)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6914, _S6913);
    float3  _S6915 = _S6914.differential_0 + _S6765;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6916;
    (&_S6916)->primal_0 = R_36;
    (&_S6916)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6917;
    (&_S6917)->primal_0 = pos_i_11;
    (&_S6917)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6916, &_S6917, _S6915);
    DiffPair_float_0 _S6918;
    (&_S6918)->primal_0 = _S6247;
    (&_S6918)->differential_0 = 0.0f;
    DiffPair_float_0 _S6919;
    (&_S6919)->primal_0 = _S6249;
    (&_S6919)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6918, &_S6919, _S6909.differential_0);
    DiffPair_float_0 _S6920;
    (&_S6920)->primal_0 = _S6246;
    (&_S6920)->differential_0 = 0.0f;
    DiffPair_float_0 _S6921;
    (&_S6921)->primal_0 = _S6249;
    (&_S6921)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6920, &_S6921, _S6911.differential_0);
    float _S6922 = _S6919.differential_0 + _S6921.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6923;
    (&_S6923)->primal_0 = _S6248;
    (&_S6923)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6923, _S6922);
    float3  _S6924 = _S6923.differential_0 + _S6792;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6925;
    (&_S6925)->primal_0 = R_36;
    (&_S6925)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6926;
    (&_S6926)->primal_0 = pos_i_10;
    (&_S6926)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6925, &_S6926, _S6924);
    DiffPair_float_0 _S6927;
    (&_S6927)->primal_0 = _S6243;
    (&_S6927)->differential_0 = 0.0f;
    DiffPair_float_0 _S6928;
    (&_S6928)->primal_0 = _S6245;
    (&_S6928)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6927, &_S6928, _S6918.differential_0);
    DiffPair_float_0 _S6929;
    (&_S6929)->primal_0 = _S6242;
    (&_S6929)->differential_0 = 0.0f;
    DiffPair_float_0 _S6930;
    (&_S6930)->primal_0 = _S6245;
    (&_S6930)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6929, &_S6930, _S6920.differential_0);
    float _S6931 = _S6928.differential_0 + _S6930.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6932;
    (&_S6932)->primal_0 = _S6244;
    (&_S6932)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6932, _S6931);
    float3  _S6933 = _S6932.differential_0 + _S6819;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6934;
    (&_S6934)->primal_0 = R_36;
    (&_S6934)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6935;
    (&_S6935)->primal_0 = pos_i_9;
    (&_S6935)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6934, &_S6935, _S6933);
    DiffPair_float_0 _S6936;
    (&_S6936)->primal_0 = _S6239;
    (&_S6936)->differential_0 = 0.0f;
    DiffPair_float_0 _S6937;
    (&_S6937)->primal_0 = _S6241;
    (&_S6937)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6936, &_S6937, _S6927.differential_0);
    DiffPair_float_0 _S6938;
    (&_S6938)->primal_0 = _S6238;
    (&_S6938)->differential_0 = 0.0f;
    DiffPair_float_0 _S6939;
    (&_S6939)->primal_0 = _S6241;
    (&_S6939)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6938, &_S6939, _S6929.differential_0);
    float _S6940 = _S6937.differential_0 + _S6939.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6941;
    (&_S6941)->primal_0 = _S6240;
    (&_S6941)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6941, _S6940);
    float3  _S6942 = _S6941.differential_0 + _S6846;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6943;
    (&_S6943)->primal_0 = R_36;
    (&_S6943)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6944;
    (&_S6944)->primal_0 = pos_i_8;
    (&_S6944)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6943, &_S6944, _S6942);
    DiffPair_float_0 _S6945;
    (&_S6945)->primal_0 = _S6235;
    (&_S6945)->differential_0 = 0.0f;
    DiffPair_float_0 _S6946;
    (&_S6946)->primal_0 = _S6237;
    (&_S6946)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6945, &_S6946, _S6936.differential_0);
    DiffPair_float_0 _S6947;
    (&_S6947)->primal_0 = _S6234;
    (&_S6947)->differential_0 = 0.0f;
    DiffPair_float_0 _S6948;
    (&_S6948)->primal_0 = _S6237;
    (&_S6948)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6947, &_S6948, _S6938.differential_0);
    float _S6949 = _S6946.differential_0 + _S6948.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6950;
    (&_S6950)->primal_0 = _S6236;
    (&_S6950)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6950, _S6949);
    float3  _S6951 = _S6950.differential_0 + _S6873;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6952;
    (&_S6952)->primal_0 = R_36;
    (&_S6952)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6953;
    (&_S6953)->primal_0 = pos_i_7;
    (&_S6953)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6952, &_S6953, _S6951);
    DiffPair_float_0 _S6954;
    (&_S6954)->primal_0 = 0.0f;
    (&_S6954)->differential_0 = 0.0f;
    DiffPair_float_0 _S6955;
    (&_S6955)->primal_0 = _S6233;
    (&_S6955)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6954, &_S6955, _S6945.differential_0);
    DiffPair_float_0 _S6956;
    (&_S6956)->primal_0 = 1.00000001504746622e+30f;
    (&_S6956)->differential_0 = 0.0f;
    DiffPair_float_0 _S6957;
    (&_S6957)->primal_0 = _S6233;
    (&_S6957)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6956, &_S6957, _S6947.differential_0);
    float _S6958 = _S6955.differential_0 + _S6957.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6959;
    (&_S6959)->primal_0 = _S6232;
    (&_S6959)->differential_0 = _S6516;
    s_bwd_length_impl_0(&_S6959, _S6958);
    float3  _S6960 = _S6959.differential_0 + make_float3 (_S6888.x, _S6888.y, _S6885);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6961;
    (&_S6961)->primal_0 = R_36;
    (&_S6961)->differential_0 = _S6578;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6962;
    (&_S6962)->primal_0 = pos_5;
    (&_S6962)->differential_0 = _S6516;
    s_bwd_prop_mul_1(&_S6961, &_S6962, _S6960);
    float3  _S6963 = _S6576 + _S6897 + _S6906 + _S6915 + _S6924 + _S6933 + _S6942 + _S6951 + _S6960 + _S6581.differential_0;
    Matrix<float, 3, 3>  _S6964 = _S6889.differential_0 + _S6898.differential_0 + _S6907.differential_0 + _S6916.differential_0 + _S6925.differential_0 + _S6934.differential_0 + _S6943.differential_0 + _S6952.differential_0 + _S6961.differential_0 + _S6582;
    (*v_densities_1)[int(0)] = 0.0f;
    (*v_densities_1)[int(1)] = 0.0f;
    (*v_densities_1)[int(2)] = 0.0f;
    (*v_densities_1)[int(3)] = 0.0f;
    (*v_densities_1)[int(4)] = 0.0f;
    (*v_densities_1)[int(5)] = 0.0f;
    (*v_densities_1)[int(6)] = 0.0f;
    (*v_densities_1)[int(7)] = 0.0f;
    (*v_sh_coeffs_10)[int(0)] = _S6654;
    (*v_sh_coeffs_10)[int(1)] = _S6655;
    (*v_sh_coeffs_10)[int(2)] = _S6656;
    (*v_sh_coeffs_10)[int(3)] = _S6657;
    (*v_sh_coeffs_10)[int(4)] = _S6658;
    (*v_sh_coeffs_10)[int(5)] = _S6659;
    (*v_sh_coeffs_10)[int(6)] = _S6660;
    (*v_sh_coeffs_10)[int(7)] = _S6661;
    (*v_sh_coeffs_10)[int(8)] = _S6662;
    (*v_sh_coeffs_10)[int(9)] = _S6663;
    (*v_sh_coeffs_10)[int(10)] = _S6664;
    (*v_sh_coeffs_10)[int(11)] = _S6665;
    (*v_sh_coeffs_10)[int(12)] = _S6666;
    (*v_sh_coeffs_10)[int(13)] = _S6667;
    (*v_sh_coeffs_10)[int(14)] = _S6668;
    (*v_sh_coeffs_10)[int(15)] = _S6669;
    *v_R_11 = _S6964;
    *v_t_11 = _S6963;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_19, float3  dOut_21)
{
    float3  _S6965 = _slang_select(((*dpx_19).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_19).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_21;
    dpx_19->primal_0 = (*dpx_19).primal_0;
    dpx_19->differential_0 = _S6965;
    return;
}

inline __device__ float3  abs_0(float3  x_69)
{
    float3  result_15;
    int i_11 = int(0);
    for(;;)
    {
        if(i_11 < int(3))
        {
        }
        else
        {
            break;
        }
        *_slang_vector_get_element_ptr(&result_15, i_11) = (F32_abs((_slang_vector_get_element(x_69, i_11))));
        i_11 = i_11 + int(1);
    }
    return result_15;
}

inline __device__ bool ray_aabb_intersection(float3  ray_o_11, float3  ray_d_11, float3  center_0, float radius_0, float * t0_0, float * t1_0)
{
    float3  m_2 = make_float3 (1.0f) / ray_d_11;
    float3  k_16 = abs_0(m_2) * make_float3 (radius_0);
    float3  _S6966 = - (m_2 * (ray_o_11 - center_0));
    float3  ta_0 = _S6966 - k_16;
    float3  tb_0 = _S6966 + k_16;
    *t0_0 = (F32_max(((F32_max((ta_0.x), (ta_0.y)))), ((F32_max((ta_0.z), (0.0f))))));
    float _S6967 = (F32_min(((F32_min((tb_0.x), (tb_0.y)))), (tb_0.z)));
    *t1_0 = _S6967;
    return (*t0_0) < _S6967;
}

inline __device__ float evaluate_alpha_voxel(float3  pos_6, float size_6, FixedArray<float, 8>  * densities_6, float3  ray_o_12, float3  ray_d_12)
{
    float _S6968 = 0.5f * size_6;
    float3  m_3 = make_float3 (1.0f) / ray_d_12;
    float3  k_17 = abs_0(m_3) * make_float3 (_S6968);
    float3  _S6969 = - (m_3 * (ray_o_12 - (pos_6 + make_float3 (_S6968))));
    float3  ta_1 = _S6969 - k_17;
    float3  tb_1 = _S6969 + k_17;
    float t0_1 = (F32_max(((F32_max((ta_1.x), (ta_1.y)))), ((F32_max((ta_1.z), (0.0f))))));
    float t1_1 = (F32_min(((F32_min((tb_1.x), (tb_1.y)))), (tb_1.z)));
    if(!(t0_1 < t1_1))
    {
        return 0.0f;
    }
    return (F32_min((1.0f - (F32_exp((- (t1_1 - t0_1) * (F32_max((((*densities_6)[int(0)] + (*densities_6)[int(1)] + (*densities_6)[int(2)] + (*densities_6)[int(3)] + (*densities_6)[int(4)] + (*densities_6)[int(5)] + (*densities_6)[int(6)] + (*densities_6)[int(7)]) / 8.0f), (0.0f))))))), (0.99900001287460327f)));
}

struct DiffPair_arrayx3Cfloatx2C8x3E_0
{
    FixedArray<float, 8>  primal_0;
    FixedArray<float, 8>  differential_0;
};

inline __device__ float3  s_primal_ctx_abs_1(float3  _S6970)
{
    return abs_0(_S6970);
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S6971, float3  _S6972)
{
    _d_abs_vector_0(_S6971, _S6972);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_voxel_0(float3  pos_7, float size_7, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_6, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_6, float _s_dOut_12)
{
    FixedArray<float, 8>  _S6973 = dpdensities_0->primal_0;
    float _S6974 = 0.5f * size_7;
    float3  _S6975 = make_float3 (_S6974);
    float3  m_4 = make_float3 (1.0f) / (*dpray_d_6).primal_0;
    float3  _S6976 = (*dpray_d_6).primal_0 * (*dpray_d_6).primal_0;
    float3  _S6977 = (*dpray_o_6).primal_0 - (pos_7 + make_float3 (_S6974));
    float3  k_18 = s_primal_ctx_abs_1(m_4) * make_float3 (_S6974);
    float3  _S6978 = - (m_4 * _S6977);
    float3  ta_2 = _S6978 - k_18;
    float3  tb_2 = _S6978 + k_18;
    float _S6979 = ta_2.x;
    float _S6980 = ta_2.y;
    float _S6981 = s_primal_ctx_max_0(_S6979, _S6980);
    float _S6982 = ta_2.z;
    float _S6983 = s_primal_ctx_max_0(_S6982, 0.0f);
    float _S6984 = s_primal_ctx_max_0(_S6981, _S6983);
    float _S6985 = tb_2.x;
    float _S6986 = tb_2.y;
    float _S6987 = s_primal_ctx_min_0(_S6985, _S6986);
    float _S6988 = tb_2.z;
    float _S6989 = s_primal_ctx_min_0(_S6987, _S6988);
    bool _S6990 = !!(_S6984 < _S6989);
    float _S6991;
    float _S6992;
    float _S6993;
    float _S6994;
    float _S6995;
    if(_S6990)
    {
        float density_0 = (_S6973[int(0)] + _S6973[int(1)] + _S6973[int(2)] + _S6973[int(3)] + _S6973[int(4)] + _S6973[int(5)] + _S6973[int(6)] + _S6973[int(7)]) / 8.0f;
        float _S6996 = - (_S6989 - _S6984);
        float _S6997 = s_primal_ctx_max_0(density_0, 0.0f);
        float _S6998 = _S6996 * _S6997;
        _S6991 = 1.0f - s_primal_ctx_exp_1(_S6998);
        _S6992 = _S6998;
        _S6993 = _S6996;
        _S6994 = _S6997;
        _S6995 = density_0;
    }
    else
    {
        _S6991 = 0.0f;
        _S6992 = 0.0f;
        _S6993 = 0.0f;
        _S6994 = 0.0f;
        _S6995 = 0.0f;
    }
    FixedArray<float, 8>  _S6999;
    if(_S6990)
    {
        DiffPair_float_0 _S7000;
        (&_S7000)->primal_0 = _S6991;
        (&_S7000)->differential_0 = 0.0f;
        DiffPair_float_0 _S7001;
        (&_S7001)->primal_0 = 0.99900001287460327f;
        (&_S7001)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S7000, &_S7001, _s_dOut_12);
        float _S7002 = - _S7000.differential_0;
        DiffPair_float_0 _S7003;
        (&_S7003)->primal_0 = _S6992;
        (&_S7003)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S7003, _S7002);
        float _S7004 = _S6993 * _S7003.differential_0;
        float _S7005 = _S6994 * _S7003.differential_0;
        DiffPair_float_0 _S7006;
        (&_S7006)->primal_0 = _S6995;
        (&_S7006)->differential_0 = 0.0f;
        DiffPair_float_0 _S7007;
        (&_S7007)->primal_0 = 0.0f;
        (&_S7007)->differential_0 = 0.0f;
        s_bwd_prop_max_1(&_S7006, &_S7007, _S7004);
        float _S7008 = - _S7005;
        float _S7009 = - _S7008;
        float _S7010 = 0.125f * _S7006.differential_0;
        FixedArray<float, 8>  _S7011;
        _S7011[int(0)] = 0.0f;
        _S7011[int(1)] = 0.0f;
        _S7011[int(2)] = 0.0f;
        _S7011[int(3)] = 0.0f;
        _S7011[int(4)] = 0.0f;
        _S7011[int(5)] = 0.0f;
        _S7011[int(6)] = 0.0f;
        _S7011[int(7)] = 0.0f;
        _S7011[int(7)] = _S7010;
        _S7011[int(6)] = _S7010;
        _S7011[int(5)] = _S7010;
        _S7011[int(4)] = _S7010;
        _S7011[int(3)] = _S7010;
        _S7011[int(2)] = _S7010;
        _S7011[int(1)] = _S7010;
        _S7011[int(0)] = _S7010;
        _S6991 = _S7008;
        _S6992 = _S7009;
        _S6999[int(0)] = _S7011[int(0)];
        _S6999[int(1)] = _S7011[int(1)];
        _S6999[int(2)] = _S7011[int(2)];
        _S6999[int(3)] = _S7011[int(3)];
        _S6999[int(4)] = _S7011[int(4)];
        _S6999[int(5)] = _S7011[int(5)];
        _S6999[int(6)] = _S7011[int(6)];
        _S6999[int(7)] = _S7011[int(7)];
    }
    else
    {
        _S6991 = 0.0f;
        _S6992 = 0.0f;
        _S6999[int(0)] = 0.0f;
        _S6999[int(1)] = 0.0f;
        _S6999[int(2)] = 0.0f;
        _S6999[int(3)] = 0.0f;
        _S6999[int(4)] = 0.0f;
        _S6999[int(5)] = 0.0f;
        _S6999[int(6)] = 0.0f;
        _S6999[int(7)] = 0.0f;
    }
    DiffPair_float_0 _S7012;
    (&_S7012)->primal_0 = _S6987;
    (&_S7012)->differential_0 = 0.0f;
    DiffPair_float_0 _S7013;
    (&_S7013)->primal_0 = _S6988;
    (&_S7013)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7012, &_S7013, _S6991);
    DiffPair_float_0 _S7014;
    (&_S7014)->primal_0 = _S6985;
    (&_S7014)->differential_0 = 0.0f;
    DiffPair_float_0 _S7015;
    (&_S7015)->primal_0 = _S6986;
    (&_S7015)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7014, &_S7015, _S7012.differential_0);
    DiffPair_float_0 _S7016;
    (&_S7016)->primal_0 = _S6981;
    (&_S7016)->differential_0 = 0.0f;
    DiffPair_float_0 _S7017;
    (&_S7017)->primal_0 = _S6983;
    (&_S7017)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7016, &_S7017, _S6992);
    DiffPair_float_0 _S7018;
    (&_S7018)->primal_0 = _S6982;
    (&_S7018)->differential_0 = 0.0f;
    DiffPair_float_0 _S7019;
    (&_S7019)->primal_0 = 0.0f;
    (&_S7019)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7018, &_S7019, _S7017.differential_0);
    DiffPair_float_0 _S7020;
    (&_S7020)->primal_0 = _S6979;
    (&_S7020)->differential_0 = 0.0f;
    DiffPair_float_0 _S7021;
    (&_S7021)->primal_0 = _S6980;
    (&_S7021)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7020, &_S7021, _S7016.differential_0);
    float3  s_diff_tb_T_0 = make_float3 (_S7014.differential_0, _S7015.differential_0, _S7013.differential_0);
    float3  s_diff_ta_T_0 = make_float3 (_S7020.differential_0, _S7021.differential_0, _S7018.differential_0);
    float3  s_diff_n_T_0 = - (s_diff_tb_T_0 + s_diff_ta_T_0);
    float3  _S7022 = _S6975 * (s_diff_tb_T_0 + - s_diff_ta_T_0);
    float3  _S7023 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7024;
    (&_S7024)->primal_0 = m_4;
    (&_S7024)->differential_0 = _S7023;
    s_bwd_prop_abs_1(&_S7024, _S7022);
    float3  _S7025 = m_4 * s_diff_n_T_0;
    float3  _S7026 = - ((_S7024.differential_0 + _S6977 * s_diff_n_T_0) / _S6976);
    dpray_d_6->primal_0 = (*dpray_d_6).primal_0;
    dpray_d_6->differential_0 = _S7026;
    dpray_o_6->primal_0 = (*dpray_o_6).primal_0;
    dpray_o_6->differential_0 = _S7025;
    dpdensities_0->primal_0 = dpdensities_0->primal_0;
    dpdensities_0->differential_0 = _S6999;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_voxel_0(float3  _S7027, float _S7028, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S7029, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7030, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7031, float _S7032)
{
    s_bwd_prop_evaluate_alpha_voxel_0(_S7027, _S7028, _S7029, _S7030, _S7031, _S7032);
    return;
}

inline __device__ void evaluate_alpha_voxel_vjp(float3  pos_8, float size_8, FixedArray<float, 8>  * densities_7, float3  ray_o_13, float3  ray_d_13, float v_alpha_4, FixedArray<float, 8>  * v_densities_2, float3  * v_ray_o_5, float3  * v_ray_d_5)
{
    FixedArray<float, 8>  _S7033 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_0;
    (&dp_densities_0)->primal_0 = *densities_7;
    (&dp_densities_0)->differential_0 = _S7033;
    float3  _S7034 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_4;
    (&dp_ray_o_4)->primal_0 = ray_o_13;
    (&dp_ray_o_4)->differential_0 = _S7034;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_4;
    (&dp_ray_d_4)->primal_0 = ray_d_13;
    (&dp_ray_d_4)->differential_0 = _S7034;
    s_bwd_evaluate_alpha_voxel_0(pos_8, size_8, &dp_densities_0, &dp_ray_o_4, &dp_ray_d_4, v_alpha_4);
    *v_densities_2 = (&dp_densities_0)->differential_0;
    *v_ray_o_5 = dp_ray_o_4.differential_0;
    *v_ray_d_5 = dp_ray_d_4.differential_0;
    return;
}

inline __device__ void evaluate_color_voxel(float3  pos_9, float size_9, float3  rgb_16, float3  ray_o_14, float3  ray_d_14, float3  * out_rgb_1, float * depth_25)
{
    *out_rgb_1 = rgb_16;
    float _S7035 = 0.5f * size_9;
    float3  m_5 = make_float3 (1.0f) / ray_d_14;
    float3  k_19 = abs_0(m_5) * make_float3 (_S7035);
    float3  _S7036 = - (m_5 * (ray_o_14 - (pos_9 + make_float3 (_S7035))));
    float3  ta_3 = _S7036 - k_19;
    float3  tb_3 = _S7036 + k_19;
    *depth_25 = 0.5f * ((F32_max(((F32_max((ta_3.x), (ta_3.y)))), ((F32_max((ta_3.z), (0.0f)))))) + (F32_min(((F32_min((tb_3.x), (tb_3.y)))), (tb_3.z))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_voxel_0(float3  pos_10, float size_10, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_7, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_7, float3  dpout_rgb_1, float dpdepth_3)
{
    float _S7037 = 0.5f * size_10;
    float3  _S7038 = make_float3 (_S7037);
    float3  m_6 = make_float3 (1.0f) / (*dpray_d_7).primal_0;
    float3  _S7039 = (*dpray_d_7).primal_0 * (*dpray_d_7).primal_0;
    float3  _S7040 = (*dpray_o_7).primal_0 - (pos_10 + make_float3 (_S7037));
    float3  k_20 = s_primal_ctx_abs_1(m_6) * make_float3 (_S7037);
    float3  _S7041 = - (m_6 * _S7040);
    float3  ta_4 = _S7041 - k_20;
    float3  tb_4 = _S7041 + k_20;
    float _S7042 = ta_4.x;
    float _S7043 = ta_4.y;
    float _S7044 = s_primal_ctx_max_0(_S7042, _S7043);
    float _S7045 = ta_4.z;
    float _S7046 = s_primal_ctx_max_0(_S7045, 0.0f);
    float _S7047 = tb_4.x;
    float _S7048 = tb_4.y;
    float _S7049 = tb_4.z;
    float _S7050 = 0.5f * dpdepth_3;
    DiffPair_float_0 _S7051;
    (&_S7051)->primal_0 = s_primal_ctx_min_0(_S7047, _S7048);
    (&_S7051)->differential_0 = 0.0f;
    DiffPair_float_0 _S7052;
    (&_S7052)->primal_0 = _S7049;
    (&_S7052)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7051, &_S7052, _S7050);
    DiffPair_float_0 _S7053;
    (&_S7053)->primal_0 = _S7047;
    (&_S7053)->differential_0 = 0.0f;
    DiffPair_float_0 _S7054;
    (&_S7054)->primal_0 = _S7048;
    (&_S7054)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7053, &_S7054, _S7051.differential_0);
    DiffPair_float_0 _S7055;
    (&_S7055)->primal_0 = _S7044;
    (&_S7055)->differential_0 = 0.0f;
    DiffPair_float_0 _S7056;
    (&_S7056)->primal_0 = _S7046;
    (&_S7056)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7055, &_S7056, _S7050);
    DiffPair_float_0 _S7057;
    (&_S7057)->primal_0 = _S7045;
    (&_S7057)->differential_0 = 0.0f;
    DiffPair_float_0 _S7058;
    (&_S7058)->primal_0 = 0.0f;
    (&_S7058)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7057, &_S7058, _S7056.differential_0);
    DiffPair_float_0 _S7059;
    (&_S7059)->primal_0 = _S7042;
    (&_S7059)->differential_0 = 0.0f;
    DiffPair_float_0 _S7060;
    (&_S7060)->primal_0 = _S7043;
    (&_S7060)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7059, &_S7060, _S7055.differential_0);
    float3  s_diff_tb_T_1 = make_float3 (_S7053.differential_0, _S7054.differential_0, _S7052.differential_0);
    float3  s_diff_ta_T_1 = make_float3 (_S7059.differential_0, _S7060.differential_0, _S7057.differential_0);
    float3  s_diff_n_T_1 = - (s_diff_tb_T_1 + s_diff_ta_T_1);
    float3  _S7061 = _S7038 * (s_diff_tb_T_1 + - s_diff_ta_T_1);
    float3  _S7062 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7063;
    (&_S7063)->primal_0 = m_6;
    (&_S7063)->differential_0 = _S7062;
    s_bwd_prop_abs_1(&_S7063, _S7061);
    float3  _S7064 = m_6 * s_diff_n_T_1;
    float3  _S7065 = - ((_S7063.differential_0 + _S7040 * s_diff_n_T_1) / _S7039);
    dpray_d_7->primal_0 = (*dpray_d_7).primal_0;
    dpray_d_7->differential_0 = _S7065;
    dpray_o_7->primal_0 = (*dpray_o_7).primal_0;
    dpray_o_7->differential_0 = _S7064;
    dprgb_1->primal_0 = (*dprgb_1).primal_0;
    dprgb_1->differential_0 = dpout_rgb_1;
    return;
}

inline __device__ void s_bwd_evaluate_color_voxel_0(float3  _S7066, float _S7067, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7068, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7069, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7070, float3  _S7071, float _S7072)
{
    s_bwd_prop_evaluate_color_voxel_0(_S7066, _S7067, _S7068, _S7069, _S7070, _S7071, _S7072);
    return;
}

inline __device__ void evaluate_color_voxel_vjp(float3  pos_11, float size_11, float3  rgb_17, float3  ray_o_15, float3  ray_d_15, float3  v_out_rgb_1, float v_depth_12, float3  * v_rgb_10, float3  * v_ray_o_6, float3  * v_ray_d_6)
{
    float3  _S7073 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_1;
    (&dp_rgb_1)->primal_0 = rgb_17;
    (&dp_rgb_1)->differential_0 = _S7073;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_5;
    (&dp_ray_o_5)->primal_0 = ray_o_15;
    (&dp_ray_o_5)->differential_0 = _S7073;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_5;
    (&dp_ray_d_5)->primal_0 = ray_d_15;
    (&dp_ray_d_5)->differential_0 = _S7073;
    s_bwd_evaluate_color_voxel_0(pos_11, size_11, &dp_rgb_1, &dp_ray_o_5, &dp_ray_d_5, v_out_rgb_1, v_depth_12);
    *v_rgb_10 = dp_rgb_1.differential_0;
    *v_ray_o_6 = dp_ray_o_5.differential_0;
    *v_ray_d_6 = dp_ray_d_5.differential_0;
    return;
}

