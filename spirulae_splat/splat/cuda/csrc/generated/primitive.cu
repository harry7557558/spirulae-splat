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

inline __device__ float3  normalize_0(float3  x_10)
{
    return x_10 / make_float3 (length_1(x_10));
}

inline __device__ float2  normalize_1(float2  x_11)
{
    return x_11 / make_float2 (length_0(x_11));
}

inline __device__ bool generate_ray(float2  uv_6, bool is_fisheye_3, FixedArray<float, 10>  * dist_coeffs_6, float3  * raydir_2)
{
    float2  _S65 = uv_6;
    bool _S66 = undistort_point_0(uv_6, dist_coeffs_6, int(8), &_S65);
    if(!_S66)
    {
        int3  _S67 = make_int3 (int(0));
        float3  _S68 = make_float3 ((float)_S67.x, (float)_S67.y, (float)_S67.z);
        *raydir_2 = _S68;
        return false;
    }
    float3  _S69;
    if(is_fisheye_3)
    {
        float2  _S70 = _S65;
        float theta_3 = length_0(_S65);
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
        _S69 = _S72;
    }
    else
    {
        _S69 = make_float3 (_S65.x, _S65.y, 1.0f);
    }
    *raydir_2 = normalize_0(_S69);
    return true;
}

inline __device__ void _d_mul_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * left_8, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * right_8, float3  dOut_9)
{
    float _S73 = (*right_8).primal_0.rows[int(0)].x * dOut_9.x;
    Matrix<float, 3, 3>  right_d_result_4;
    *&(((&right_d_result_4)->rows + (int(0)))->x) = (*left_8).primal_0.x * dOut_9.x;
    float sum_14 = _S73 + (*right_8).primal_0.rows[int(0)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(0)))->y) = (*left_8).primal_0.x * dOut_9.y;
    float sum_15 = sum_14 + (*right_8).primal_0.rows[int(0)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(0)))->z) = (*left_8).primal_0.x * dOut_9.z;
    float3  left_d_result_4;
    *&((&left_d_result_4)->x) = sum_15;
    float _S74 = (*right_8).primal_0.rows[int(1)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(1)))->x) = (*left_8).primal_0.y * dOut_9.x;
    float sum_16 = _S74 + (*right_8).primal_0.rows[int(1)].y * dOut_9.y;
    *&(((&right_d_result_4)->rows + (int(1)))->y) = (*left_8).primal_0.y * dOut_9.y;
    float sum_17 = sum_16 + (*right_8).primal_0.rows[int(1)].z * dOut_9.z;
    *&(((&right_d_result_4)->rows + (int(1)))->z) = (*left_8).primal_0.y * dOut_9.z;
    *&((&left_d_result_4)->y) = sum_17;
    float _S75 = (*right_8).primal_0.rows[int(2)].x * dOut_9.x;
    *&(((&right_d_result_4)->rows + (int(2)))->x) = (*left_8).primal_0.z * dOut_9.x;
    float sum_18 = _S75 + (*right_8).primal_0.rows[int(2)].y * dOut_9.y;
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

inline __device__ float3  transform_ray_o(Matrix<float, 3, 3>  R_2, float3  t_1)
{
    return - mul_7(t_1, R_2);
}

inline __device__ float3  transform_ray_d(Matrix<float, 3, 3>  R_3, float3  raydir_3)
{
    return mul_7(raydir_3, R_3);
}

inline __device__ void s_bwd_prop_mul_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S76, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S77, float3  _S78)
{
    _d_mul_1(_S76, _S77, _S78);
    return;
}

inline __device__ void s_bwd_prop_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpt_0, float3  _s_dOut_0)
{
    float3  _S79 = - _s_dOut_0;
    float3  _S80 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S81;
    (&_S81)->primal_0 = (*dpt_0).primal_0;
    (&_S81)->differential_0 = _S80;
    Matrix<float, 3, 3>  _S82 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S83;
    (&_S83)->primal_0 = (*dpR_0).primal_0;
    (&_S83)->differential_0 = _S82;
    s_bwd_prop_mul_0(&_S81, &_S83, _S79);
    dpt_0->primal_0 = (*dpt_0).primal_0;
    dpt_0->differential_0 = _S81.differential_0;
    dpR_0->primal_0 = (*dpR_0).primal_0;
    dpR_0->differential_0 = _S83.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_o_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S84, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S85, float3  _S86)
{
    s_bwd_prop_transform_ray_o_0(_S84, _S85, _S86);
    return;
}

inline __device__ void transform_ray_o_vjp(Matrix<float, 3, 3>  R_4, float3  t_2, float3  v_ray_o_0, Matrix<float, 3, 3>  * v_R_0, float3  * v_t_0)
{
    Matrix<float, 3, 3>  _S87 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_0;
    (&dp_R_0)->primal_0 = R_4;
    (&dp_R_0)->differential_0 = _S87;
    float3  _S88 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_t_0;
    (&dp_t_0)->primal_0 = t_2;
    (&dp_t_0)->differential_0 = _S88;
    s_bwd_transform_ray_o_0(&dp_R_0, &dp_t_0, v_ray_o_0);
    *v_R_0 = dp_R_0.differential_0;
    *v_t_0 = dp_t_0.differential_0;
    return;
}

inline __device__ void s_bwd_prop_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpR_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpraydir_0, float3  _s_dOut_1)
{
    float3  _S89 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S90;
    (&_S90)->primal_0 = (*dpraydir_0).primal_0;
    (&_S90)->differential_0 = _S89;
    Matrix<float, 3, 3>  _S91 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S92;
    (&_S92)->primal_0 = (*dpR_1).primal_0;
    (&_S92)->differential_0 = _S91;
    s_bwd_prop_mul_0(&_S90, &_S92, _s_dOut_1);
    dpraydir_0->primal_0 = (*dpraydir_0).primal_0;
    dpraydir_0->differential_0 = _S90.differential_0;
    dpR_1->primal_0 = (*dpR_1).primal_0;
    dpR_1->differential_0 = _S92.differential_0;
    return;
}

inline __device__ void s_bwd_transform_ray_d_0(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S93, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S94, float3  _S95)
{
    s_bwd_prop_transform_ray_d_0(_S93, _S94, _S95);
    return;
}

inline __device__ void transform_ray_d_vjp(Matrix<float, 3, 3>  R_5, float3  raydir_4, float3  v_ray_d_0, Matrix<float, 3, 3>  * v_R_1, float3  * v_raydir_0)
{
    Matrix<float, 3, 3>  _S96 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_R_1;
    (&dp_R_1)->primal_0 = R_5;
    (&dp_R_1)->differential_0 = _S96;
    float3  _S97 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_raydir_0;
    (&dp_raydir_0)->primal_0 = raydir_4;
    (&dp_raydir_0)->differential_0 = _S97;
    s_bwd_transform_ray_d_0(&dp_R_1, &dp_raydir_0, v_ray_d_0);
    *v_R_1 = dp_R_1.differential_0;
    *v_raydir_0 = dp_raydir_0.differential_0;
    return;
}

inline __device__ void _d_exp_0(DiffPair_float_0 * dpx_5, float dOut_10)
{
    float _S98 = (F32_exp(((*dpx_5).primal_0))) * dOut_10;
    dpx_5->primal_0 = (*dpx_5).primal_0;
    dpx_5->differential_0 = _S98;
    return;
}

inline __device__ void map_opaque_triangle(float3  mean_0, float4  quat_3, float3  scale_2, float3  * vert0_0, float3  * vert1_0, float3  * vert2_0)
{
    float _S99 = scale_2.x;
    float sx_0 = (F32_exp((_S99)));
    float _S100 = scale_2.y;
    float sy_0 = (F32_exp((_S100)));
    float sz_0 = scale_2.z - 0.5f * (_S99 + _S100);
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
    Matrix<float, 3, 3>  _S101 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_3 + z2_3), 2.0f * (xy_3 + wz_3), 2.0f * (xz_3 - wy_3), 2.0f * (xy_3 - wz_3), 1.0f - 2.0f * (x2_3 + z2_3), 2.0f * (yz_3 + wx_3), 2.0f * (xz_3 + wy_3), 2.0f * (yz_3 - wx_3), 1.0f - 2.0f * (x2_3 + y2_3)));
    *vert0_0 = mul_0(_S101, make_float3 (sx_0, 0.0f, 0.0f)) + mean_0;
    *vert1_0 = mul_0(_S101, make_float3 (sx_0 * (-0.5f + sz_0), sy_0, 0.0f)) + mean_0;
    *vert2_0 = mul_0(_S101, make_float3 (sx_0 * (-0.5f - sz_0), - sy_0, 0.0f)) + mean_0;
    return;
}

inline __device__ void persp_proj_3dgs(float3  mean3d_0, Matrix<float, 3, 3>  cov3d_0, float fx_0, float fy_0, float cx_0, float cy_0, uint width_0, uint height_0, Matrix<float, 2, 2>  * cov2d_0, float2  * mean2d_0)
{
    float _S102 = float(width_0);
    float _S103 = float(height_0);
    float _S104 = 0.30000001192092896f * (0.5f * _S102 / fx_0);
    float _S105 = 0.30000001192092896f * (0.5f * _S103 / fy_0);
    float rz_0 = 1.0f / mean3d_0.z;
    float rz2_0 = rz_0 * rz_0;
    Matrix<float, 2, 3>  J_0 = makeMatrix<float, 2, 3> (fx_0 * rz_0, 0.0f, - fx_0 * (mean3d_0.z * (F32_min(((_S102 - cx_0) / fx_0 + _S104), ((F32_max((- (cx_0 / fx_0 + _S104)), (mean3d_0.x * rz_0))))))) * rz2_0, 0.0f, fy_0 * rz_0, - fy_0 * (mean3d_0.z * (F32_min(((_S103 - cy_0) / fy_0 + _S105), ((F32_max((- (cy_0 / fy_0 + _S105)), (mean3d_0.y * rz_0))))))) * rz2_0);
    *cov2d_0 = mul_6(mul_5(J_0, cov3d_0), transpose_1(J_0));
    *mean2d_0 = make_float2 (fx_0 * mean3d_0.x * rz_0 + cx_0, fy_0 * mean3d_0.y * rz_0 + cy_0);
    return;
}

inline __device__ void _d_clamp_0(DiffPair_float_0 * dpx_6, DiffPair_float_0 * dpMin_0, DiffPair_float_0 * dpMax_0, float dOut_11)
{
    DiffPair_float_0 _S106 = *dpx_6;
    bool _S107;
    if(((*dpx_6).primal_0) >= ((*dpMin_0).primal_0))
    {
        _S107 = ((*dpx_6).primal_0) <= ((*dpMax_0).primal_0);
    }
    else
    {
        _S107 = false;
    }
    float _S108;
    if(_S107)
    {
        _S108 = dOut_11;
    }
    else
    {
        _S108 = 0.0f;
    }
    dpx_6->primal_0 = _S106.primal_0;
    dpx_6->differential_0 = _S108;
    DiffPair_float_0 _S109 = *dpMin_0;
    if((_S106.primal_0) < ((*dpMin_0).primal_0))
    {
        _S108 = dOut_11;
    }
    else
    {
        _S108 = 0.0f;
    }
    dpMin_0->primal_0 = _S109.primal_0;
    dpMin_0->differential_0 = _S108;
    DiffPair_float_0 _S110 = *dpMax_0;
    if(((*dpx_6).primal_0) > ((*dpMax_0).primal_0))
    {
        _S108 = dOut_11;
    }
    else
    {
        _S108 = 0.0f;
    }
    dpMax_0->primal_0 = _S110.primal_0;
    dpMax_0->differential_0 = _S108;
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

inline __device__ bool persp_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_0, float4  intrins_0, FixedArray<float, 10>  * dist_coeffs_7, uint width_1, uint height_1, Matrix<float, 2, 2>  * cov2d_1, float2  * mean2d_1)
{
    float2  * _S111;
    float2  * _S112;
    float2  * _S113;
    bool _S114;
    float2  * _S115;
    float2  * _S116;
    float2  * _S117;
    bool _S118;
    float2  * _S119;
    bool _S120;
    int2  _S121 = make_int2 (int(0));
    float2  _S122 = make_float2 ((float)_S121.x, (float)_S121.y);
    *mean2d_1 = _S122;
    *cov2d_1 = makeMatrix<float, 2, 2> (0.0f);
    float fx_1 = intrins_0.x;
    float fy_1 = intrins_0.y;
    float _S123 = float(width_1);
    float _S124 = float(height_1);
    float _S125 = 0.30000001192092896f * (0.5f * _S123 / fx_1) * fx_1;
    float lim_x_pos_0 = _S123 + _S125;
    float _S126 = 0.30000001192092896f * (0.5f * _S124 / fy_1) * fy_1;
    float lim_y_pos_0 = _S124 + _S126;
    FixedArray<float2 , 7>  proj_points_0;
    for(;;)
    {
        bool _S127;
        _S111 = &proj_points_0[int(0)];
        for(;;)
        {
            float _S128 = sigmas_0->p_0[int(0)].z;
            proj_points_0[int(0)] = float2 {sigmas_0->p_0[int(0)].x, sigmas_0->p_0[int(0)].y} / make_float2 (_S128);
            if(_S128 < 0.0f)
            {
                _S127 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S129 = camera_distortion_jac_0(proj_points_0[int(0)], dist_coeffs_7);
                _S127 = !((F32_min((determinant_0(_S129)), ((F32_min((_S129.rows[int(0)].x), (_S129.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S127)
            {
                break;
            }
            float u_4 = proj_points_0[int(0)].x;
            float v_4 = proj_points_0[int(0)].y;
            float r2_4 = u_4 * u_4 + v_4 * v_4;
            float2  _S130 = proj_points_0[int(0)] * make_float2 (1.0f + r2_4 * ((*dist_coeffs_7)[int(0)] + r2_4 * ((*dist_coeffs_7)[int(1)] + r2_4 * ((*dist_coeffs_7)[int(2)] + r2_4 * (*dist_coeffs_7)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_7)[int(4)] * u_4 * v_4 + (*dist_coeffs_7)[int(5)] * (r2_4 + 2.0f * u_4 * u_4) + (*dist_coeffs_7)[int(6)] * r2_4, 2.0f * (*dist_coeffs_7)[int(5)] * u_4 * v_4 + (*dist_coeffs_7)[int(4)] * (r2_4 + 2.0f * v_4 * v_4) + (*dist_coeffs_7)[int(7)] * r2_4);
            float2  _S131 = _S130 + make_float2 ((*dist_coeffs_7)[int(8)] * _S130.x + (*dist_coeffs_7)[int(9)] * _S130.y, 0.0f);
            proj_points_0[int(0)] = make_float2 (fx_1 * _S131.x + intrins_0.z, fy_1 * _S131.y + intrins_0.w);
            break;
        }
        bool all_valid_0 = true & (!_S127);
        _S112 = &proj_points_0[int(1)];
        for(;;)
        {
            float _S132 = sigmas_0->p_0[int(1)].z;
            proj_points_0[int(1)] = float2 {sigmas_0->p_0[int(1)].x, sigmas_0->p_0[int(1)].y} / make_float2 (_S132);
            if(_S132 < 0.0f)
            {
                _S127 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S133 = camera_distortion_jac_0(proj_points_0[int(1)], dist_coeffs_7);
                _S127 = !((F32_min((determinant_0(_S133)), ((F32_min((_S133.rows[int(0)].x), (_S133.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S127)
            {
                break;
            }
            float u_5 = proj_points_0[int(1)].x;
            float v_5 = proj_points_0[int(1)].y;
            float r2_5 = u_5 * u_5 + v_5 * v_5;
            float2  _S134 = proj_points_0[int(1)] * make_float2 (1.0f + r2_5 * ((*dist_coeffs_7)[int(0)] + r2_5 * ((*dist_coeffs_7)[int(1)] + r2_5 * ((*dist_coeffs_7)[int(2)] + r2_5 * (*dist_coeffs_7)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_7)[int(4)] * u_5 * v_5 + (*dist_coeffs_7)[int(5)] * (r2_5 + 2.0f * u_5 * u_5) + (*dist_coeffs_7)[int(6)] * r2_5, 2.0f * (*dist_coeffs_7)[int(5)] * u_5 * v_5 + (*dist_coeffs_7)[int(4)] * (r2_5 + 2.0f * v_5 * v_5) + (*dist_coeffs_7)[int(7)] * r2_5);
            float2  _S135 = _S134 + make_float2 ((*dist_coeffs_7)[int(8)] * _S134.x + (*dist_coeffs_7)[int(9)] * _S134.y, 0.0f);
            proj_points_0[int(1)] = make_float2 (fx_1 * _S135.x + intrins_0.z, fy_1 * _S135.y + intrins_0.w);
            break;
        }
        bool all_valid_1 = all_valid_0 & (!_S127);
        for(;;)
        {
            _S113 = &proj_points_0[int(2)];
            for(;;)
            {
                float _S136 = sigmas_0->p_0[int(2)].z;
                proj_points_0[int(2)] = float2 {sigmas_0->p_0[int(2)].x, sigmas_0->p_0[int(2)].y} / make_float2 (_S136);
                if(_S136 < 0.0f)
                {
                    _S127 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S137 = camera_distortion_jac_0(proj_points_0[int(2)], dist_coeffs_7);
                    _S127 = !((F32_min((determinant_0(_S137)), ((F32_min((_S137.rows[int(0)].x), (_S137.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S127)
                {
                    break;
                }
                float u_6 = proj_points_0[int(2)].x;
                float v_6 = proj_points_0[int(2)].y;
                float r2_6 = u_6 * u_6 + v_6 * v_6;
                float2  _S138 = proj_points_0[int(2)] * make_float2 (1.0f + r2_6 * ((*dist_coeffs_7)[int(0)] + r2_6 * ((*dist_coeffs_7)[int(1)] + r2_6 * ((*dist_coeffs_7)[int(2)] + r2_6 * (*dist_coeffs_7)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_7)[int(4)] * u_6 * v_6 + (*dist_coeffs_7)[int(5)] * (r2_6 + 2.0f * u_6 * u_6) + (*dist_coeffs_7)[int(6)] * r2_6, 2.0f * (*dist_coeffs_7)[int(5)] * u_6 * v_6 + (*dist_coeffs_7)[int(4)] * (r2_6 + 2.0f * v_6 * v_6) + (*dist_coeffs_7)[int(7)] * r2_6);
                float2  _S139 = _S138 + make_float2 ((*dist_coeffs_7)[int(8)] * _S138.x + (*dist_coeffs_7)[int(9)] * _S138.y, 0.0f);
                proj_points_0[int(2)] = make_float2 (fx_1 * _S139.x + intrins_0.z, fy_1 * _S139.y + intrins_0.w);
                break;
            }
            _S114 = all_valid_1 & (!_S127);
            break;
        }
        _S115 = &proj_points_0[int(3)];
        for(;;)
        {
            float _S140 = sigmas_0->p_0[int(3)].z;
            proj_points_0[int(3)] = float2 {sigmas_0->p_0[int(3)].x, sigmas_0->p_0[int(3)].y} / make_float2 (_S140);
            if(_S140 < 0.0f)
            {
                _S127 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S141 = camera_distortion_jac_0(proj_points_0[int(3)], dist_coeffs_7);
                _S127 = !((F32_min((determinant_0(_S141)), ((F32_min((_S141.rows[int(0)].x), (_S141.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S127)
            {
                break;
            }
            float u_7 = proj_points_0[int(3)].x;
            float v_7 = proj_points_0[int(3)].y;
            float r2_7 = u_7 * u_7 + v_7 * v_7;
            float2  _S142 = proj_points_0[int(3)] * make_float2 (1.0f + r2_7 * ((*dist_coeffs_7)[int(0)] + r2_7 * ((*dist_coeffs_7)[int(1)] + r2_7 * ((*dist_coeffs_7)[int(2)] + r2_7 * (*dist_coeffs_7)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_7)[int(4)] * u_7 * v_7 + (*dist_coeffs_7)[int(5)] * (r2_7 + 2.0f * u_7 * u_7) + (*dist_coeffs_7)[int(6)] * r2_7, 2.0f * (*dist_coeffs_7)[int(5)] * u_7 * v_7 + (*dist_coeffs_7)[int(4)] * (r2_7 + 2.0f * v_7 * v_7) + (*dist_coeffs_7)[int(7)] * r2_7);
            float2  _S143 = _S142 + make_float2 ((*dist_coeffs_7)[int(8)] * _S142.x + (*dist_coeffs_7)[int(9)] * _S142.y, 0.0f);
            proj_points_0[int(3)] = make_float2 (fx_1 * _S143.x + intrins_0.z, fy_1 * _S143.y + intrins_0.w);
            break;
        }
        bool all_valid_2 = _S114 & (!_S127);
        _S116 = &proj_points_0[int(4)];
        for(;;)
        {
            float _S144 = sigmas_0->p_0[int(4)].z;
            proj_points_0[int(4)] = float2 {sigmas_0->p_0[int(4)].x, sigmas_0->p_0[int(4)].y} / make_float2 (_S144);
            if(_S144 < 0.0f)
            {
                _S127 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S145 = camera_distortion_jac_0(proj_points_0[int(4)], dist_coeffs_7);
                _S127 = !((F32_min((determinant_0(_S145)), ((F32_min((_S145.rows[int(0)].x), (_S145.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S127)
            {
                break;
            }
            float u_8 = proj_points_0[int(4)].x;
            float v_8 = proj_points_0[int(4)].y;
            float r2_8 = u_8 * u_8 + v_8 * v_8;
            float2  _S146 = proj_points_0[int(4)] * make_float2 (1.0f + r2_8 * ((*dist_coeffs_7)[int(0)] + r2_8 * ((*dist_coeffs_7)[int(1)] + r2_8 * ((*dist_coeffs_7)[int(2)] + r2_8 * (*dist_coeffs_7)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_7)[int(4)] * u_8 * v_8 + (*dist_coeffs_7)[int(5)] * (r2_8 + 2.0f * u_8 * u_8) + (*dist_coeffs_7)[int(6)] * r2_8, 2.0f * (*dist_coeffs_7)[int(5)] * u_8 * v_8 + (*dist_coeffs_7)[int(4)] * (r2_8 + 2.0f * v_8 * v_8) + (*dist_coeffs_7)[int(7)] * r2_8);
            float2  _S147 = _S146 + make_float2 ((*dist_coeffs_7)[int(8)] * _S146.x + (*dist_coeffs_7)[int(9)] * _S146.y, 0.0f);
            proj_points_0[int(4)] = make_float2 (fx_1 * _S147.x + intrins_0.z, fy_1 * _S147.y + intrins_0.w);
            break;
        }
        bool all_valid_3 = all_valid_2 & (!_S127);
        for(;;)
        {
            _S117 = &proj_points_0[int(5)];
            for(;;)
            {
                float _S148 = sigmas_0->p_0[int(5)].z;
                proj_points_0[int(5)] = float2 {sigmas_0->p_0[int(5)].x, sigmas_0->p_0[int(5)].y} / make_float2 (_S148);
                if(_S148 < 0.0f)
                {
                    _S127 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S149 = camera_distortion_jac_0(proj_points_0[int(5)], dist_coeffs_7);
                    _S127 = !((F32_min((determinant_0(_S149)), ((F32_min((_S149.rows[int(0)].x), (_S149.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S127)
                {
                    break;
                }
                float u_9 = proj_points_0[int(5)].x;
                float v_9 = proj_points_0[int(5)].y;
                float r2_9 = u_9 * u_9 + v_9 * v_9;
                float2  _S150 = proj_points_0[int(5)] * make_float2 (1.0f + r2_9 * ((*dist_coeffs_7)[int(0)] + r2_9 * ((*dist_coeffs_7)[int(1)] + r2_9 * ((*dist_coeffs_7)[int(2)] + r2_9 * (*dist_coeffs_7)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_7)[int(4)] * u_9 * v_9 + (*dist_coeffs_7)[int(5)] * (r2_9 + 2.0f * u_9 * u_9) + (*dist_coeffs_7)[int(6)] * r2_9, 2.0f * (*dist_coeffs_7)[int(5)] * u_9 * v_9 + (*dist_coeffs_7)[int(4)] * (r2_9 + 2.0f * v_9 * v_9) + (*dist_coeffs_7)[int(7)] * r2_9);
                float2  _S151 = _S150 + make_float2 ((*dist_coeffs_7)[int(8)] * _S150.x + (*dist_coeffs_7)[int(9)] * _S150.y, 0.0f);
                proj_points_0[int(5)] = make_float2 (fx_1 * _S151.x + intrins_0.z, fy_1 * _S151.y + intrins_0.w);
                break;
            }
            _S118 = all_valid_3 & (!_S127);
            break;
        }
        _S119 = &proj_points_0[int(6)];
        for(;;)
        {
            float _S152 = sigmas_0->p_0[int(6)].z;
            proj_points_0[int(6)] = float2 {sigmas_0->p_0[int(6)].x, sigmas_0->p_0[int(6)].y} / make_float2 (_S152);
            if(_S152 < 0.0f)
            {
                _S127 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S153 = camera_distortion_jac_0(proj_points_0[int(6)], dist_coeffs_7);
                _S127 = !((F32_min((determinant_0(_S153)), ((F32_min((_S153.rows[int(0)].x), (_S153.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S127)
            {
                break;
            }
            float u_10 = proj_points_0[int(6)].x;
            float v_10 = proj_points_0[int(6)].y;
            float r2_10 = u_10 * u_10 + v_10 * v_10;
            float2  _S154 = proj_points_0[int(6)] * make_float2 (1.0f + r2_10 * ((*dist_coeffs_7)[int(0)] + r2_10 * ((*dist_coeffs_7)[int(1)] + r2_10 * ((*dist_coeffs_7)[int(2)] + r2_10 * (*dist_coeffs_7)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_7)[int(4)] * u_10 * v_10 + (*dist_coeffs_7)[int(5)] * (r2_10 + 2.0f * u_10 * u_10) + (*dist_coeffs_7)[int(6)] * r2_10, 2.0f * (*dist_coeffs_7)[int(5)] * u_10 * v_10 + (*dist_coeffs_7)[int(4)] * (r2_10 + 2.0f * v_10 * v_10) + (*dist_coeffs_7)[int(7)] * r2_10);
            float2  _S155 = _S154 + make_float2 ((*dist_coeffs_7)[int(8)] * _S154.x + (*dist_coeffs_7)[int(9)] * _S154.y, 0.0f);
            proj_points_0[int(6)] = make_float2 (fx_1 * _S155.x + intrins_0.z, fy_1 * _S155.y + intrins_0.w);
            break;
        }
        _S120 = _S118 & (!_S127);
        break;
    }
    if(!_S120)
    {
        return false;
    }
    float2  _S156 = *mean2d_1 + make_float2 (sigmas_0->w_mean_0[int(0)]) * *_S111;
    *mean2d_1 = _S156;
    float2  _S157 = _S156 + make_float2 (sigmas_0->w_mean_0[int(1)]) * *_S112;
    *mean2d_1 = _S157;
    float2  _S158 = _S157 + make_float2 (sigmas_0->w_mean_0[int(2)]) * *_S113;
    *mean2d_1 = _S158;
    float2  _S159 = _S158 + make_float2 (sigmas_0->w_mean_0[int(3)]) * *_S115;
    *mean2d_1 = _S159;
    float2  _S160 = _S159 + make_float2 (sigmas_0->w_mean_0[int(4)]) * *_S116;
    *mean2d_1 = _S160;
    float2  _S161 = _S160 + make_float2 (sigmas_0->w_mean_0[int(5)]) * *_S117;
    *mean2d_1 = _S161;
    float2  _S162 = _S161 + make_float2 (sigmas_0->w_mean_0[int(6)]) * *_S119;
    *mean2d_1 = _S162;
    float _S163 = - _S125;
    float _S164 = - _S126;
    float2  _S165 = make_float2 (clamp_0(_S162.x, _S163, lim_x_pos_0), clamp_0(_S162.y, _S164, lim_y_pos_0));
    float2  d_0 = make_float2 (clamp_0((*_S111).x, _S163, lim_x_pos_0), clamp_0((*_S111).y, _S164, lim_y_pos_0)) - _S165;
    float _S166 = d_0.x;
    float _S167 = d_0.y;
    float _S168 = _S166 * _S167;
    Matrix<float, 2, 2>  _S169 = *cov2d_1 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S166 * _S166, _S168, _S168, _S167 * _S167);
    *cov2d_1 = _S169;
    float2  d_1 = make_float2 (clamp_0((*_S112).x, _S163, lim_x_pos_0), clamp_0((*_S112).y, _S164, lim_y_pos_0)) - _S165;
    float _S170 = d_1.x;
    float _S171 = d_1.y;
    float _S172 = _S170 * _S171;
    Matrix<float, 2, 2>  _S173 = _S169 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S170 * _S170, _S172, _S172, _S171 * _S171);
    *cov2d_1 = _S173;
    float2  d_2 = make_float2 (clamp_0((*_S113).x, _S163, lim_x_pos_0), clamp_0((*_S113).y, _S164, lim_y_pos_0)) - _S165;
    float _S174 = d_2.x;
    float _S175 = d_2.y;
    float _S176 = _S174 * _S175;
    Matrix<float, 2, 2>  _S177 = _S173 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S174 * _S174, _S176, _S176, _S175 * _S175);
    *cov2d_1 = _S177;
    float2  d_3 = make_float2 (clamp_0((*_S115).x, _S163, lim_x_pos_0), clamp_0((*_S115).y, _S164, lim_y_pos_0)) - _S165;
    float _S178 = d_3.x;
    float _S179 = d_3.y;
    float _S180 = _S178 * _S179;
    Matrix<float, 2, 2>  _S181 = _S177 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S178 * _S178, _S180, _S180, _S179 * _S179);
    *cov2d_1 = _S181;
    float2  d_4 = make_float2 (clamp_0((*_S116).x, _S163, lim_x_pos_0), clamp_0((*_S116).y, _S164, lim_y_pos_0)) - _S165;
    float _S182 = d_4.x;
    float _S183 = d_4.y;
    float _S184 = _S182 * _S183;
    Matrix<float, 2, 2>  _S185 = _S181 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S182 * _S182, _S184, _S184, _S183 * _S183);
    *cov2d_1 = _S185;
    float2  d_5 = make_float2 (clamp_0((*_S117).x, _S163, lim_x_pos_0), clamp_0((*_S117).y, _S164, lim_y_pos_0)) - _S165;
    float _S186 = d_5.x;
    float _S187 = d_5.y;
    float _S188 = _S186 * _S187;
    Matrix<float, 2, 2>  _S189 = _S185 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S186 * _S186, _S188, _S188, _S187 * _S187);
    *cov2d_1 = _S189;
    float2  d_6 = make_float2 (clamp_0((*_S119).x, _S163, lim_x_pos_0), clamp_0((*_S119).y, _S164, lim_y_pos_0)) - _S165;
    float _S190 = d_6.x;
    float _S191 = d_6.y;
    float _S192 = _S190 * _S191;
    *cov2d_1 = _S189 + makeMatrix<float, 2, 2> (sigmas_0->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S190 * _S190, _S192, _S192, _S191 * _S191);
    return true;
}

inline __device__ bool persp_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_1, float4  intrins_1, FixedArray<float, 10>  * dist_coeffs_8, uint width_2, uint height_2, Matrix<float, 2, 2>  * cov2d_2, float2  * mean2d_2)
{
    int2  _S193 = make_int2 (int(0));
    float2  _S194 = make_float2 ((float)_S193.x, (float)_S193.y);
    *mean2d_2 = _S194;
    *cov2d_2 = makeMatrix<float, 2, 2> (0.0f);
    float fx_2 = intrins_1.x;
    float fy_2 = intrins_1.y;
    float _S195 = float(width_2);
    float _S196 = float(height_2);
    float _S197 = 0.30000001192092896f * (0.5f * _S195 / fx_2) * fx_2;
    float lim_x_pos_1 = _S195 + _S197;
    float _S198 = 0.30000001192092896f * (0.5f * _S196 / fy_2) * fy_2;
    float lim_y_pos_1 = _S196 + _S198;
    FixedArray<float2 , 7>  proj_points_1;
    float2  _S199 = float2 {sigmas_1->p_0[int(0)].x, sigmas_1->p_0[int(0)].y} / make_float2 (sigmas_1->p_0[int(0)].z);
    float u_11 = _S199.x;
    float v_11 = _S199.y;
    float r2_11 = u_11 * u_11 + v_11 * v_11;
    float _S200 = 2.0f * (*dist_coeffs_8)[int(4)];
    float _S201 = 2.0f * (*dist_coeffs_8)[int(5)];
    float2  _S202 = _S199 * make_float2 (1.0f + r2_11 * ((*dist_coeffs_8)[int(0)] + r2_11 * ((*dist_coeffs_8)[int(1)] + r2_11 * ((*dist_coeffs_8)[int(2)] + r2_11 * (*dist_coeffs_8)[int(3)])))) + make_float2 (_S200 * u_11 * v_11 + (*dist_coeffs_8)[int(5)] * (r2_11 + 2.0f * u_11 * u_11) + (*dist_coeffs_8)[int(6)] * r2_11, _S201 * u_11 * v_11 + (*dist_coeffs_8)[int(4)] * (r2_11 + 2.0f * v_11 * v_11) + (*dist_coeffs_8)[int(7)] * r2_11);
    float2  _S203 = _S202 + make_float2 ((*dist_coeffs_8)[int(8)] * _S202.x + (*dist_coeffs_8)[int(9)] * _S202.y, 0.0f);
    float cx_1 = intrins_1.z;
    float cy_1 = intrins_1.w;
    float _S204 = fx_2 * _S203.x + cx_1;
    float _S205 = fy_2 * _S203.y + cy_1;
    float2  _S206 = make_float2 (_S204, _S205);
    proj_points_1[int(0)] = _S206;
    float2  _S207 = float2 {sigmas_1->p_0[int(1)].x, sigmas_1->p_0[int(1)].y} / make_float2 (sigmas_1->p_0[int(1)].z);
    float u_12 = _S207.x;
    float v_12 = _S207.y;
    float r2_12 = u_12 * u_12 + v_12 * v_12;
    float2  _S208 = _S207 * make_float2 (1.0f + r2_12 * ((*dist_coeffs_8)[int(0)] + r2_12 * ((*dist_coeffs_8)[int(1)] + r2_12 * ((*dist_coeffs_8)[int(2)] + r2_12 * (*dist_coeffs_8)[int(3)])))) + make_float2 (_S200 * u_12 * v_12 + (*dist_coeffs_8)[int(5)] * (r2_12 + 2.0f * u_12 * u_12) + (*dist_coeffs_8)[int(6)] * r2_12, _S201 * u_12 * v_12 + (*dist_coeffs_8)[int(4)] * (r2_12 + 2.0f * v_12 * v_12) + (*dist_coeffs_8)[int(7)] * r2_12);
    float2  _S209 = _S208 + make_float2 ((*dist_coeffs_8)[int(8)] * _S208.x + (*dist_coeffs_8)[int(9)] * _S208.y, 0.0f);
    float _S210 = fx_2 * _S209.x + cx_1;
    float _S211 = fy_2 * _S209.y + cy_1;
    float2  _S212 = make_float2 (_S210, _S211);
    proj_points_1[int(1)] = _S212;
    float2  _S213 = float2 {sigmas_1->p_0[int(2)].x, sigmas_1->p_0[int(2)].y} / make_float2 (sigmas_1->p_0[int(2)].z);
    float u_13 = _S213.x;
    float v_13 = _S213.y;
    float r2_13 = u_13 * u_13 + v_13 * v_13;
    float2  _S214 = _S213 * make_float2 (1.0f + r2_13 * ((*dist_coeffs_8)[int(0)] + r2_13 * ((*dist_coeffs_8)[int(1)] + r2_13 * ((*dist_coeffs_8)[int(2)] + r2_13 * (*dist_coeffs_8)[int(3)])))) + make_float2 (_S200 * u_13 * v_13 + (*dist_coeffs_8)[int(5)] * (r2_13 + 2.0f * u_13 * u_13) + (*dist_coeffs_8)[int(6)] * r2_13, _S201 * u_13 * v_13 + (*dist_coeffs_8)[int(4)] * (r2_13 + 2.0f * v_13 * v_13) + (*dist_coeffs_8)[int(7)] * r2_13);
    float2  _S215 = _S214 + make_float2 ((*dist_coeffs_8)[int(8)] * _S214.x + (*dist_coeffs_8)[int(9)] * _S214.y, 0.0f);
    float _S216 = fx_2 * _S215.x + cx_1;
    float _S217 = fy_2 * _S215.y + cy_1;
    float2  _S218 = make_float2 (_S216, _S217);
    proj_points_1[int(2)] = _S218;
    float2  _S219 = float2 {sigmas_1->p_0[int(3)].x, sigmas_1->p_0[int(3)].y} / make_float2 (sigmas_1->p_0[int(3)].z);
    float u_14 = _S219.x;
    float v_14 = _S219.y;
    float r2_14 = u_14 * u_14 + v_14 * v_14;
    float2  _S220 = _S219 * make_float2 (1.0f + r2_14 * ((*dist_coeffs_8)[int(0)] + r2_14 * ((*dist_coeffs_8)[int(1)] + r2_14 * ((*dist_coeffs_8)[int(2)] + r2_14 * (*dist_coeffs_8)[int(3)])))) + make_float2 (_S200 * u_14 * v_14 + (*dist_coeffs_8)[int(5)] * (r2_14 + 2.0f * u_14 * u_14) + (*dist_coeffs_8)[int(6)] * r2_14, _S201 * u_14 * v_14 + (*dist_coeffs_8)[int(4)] * (r2_14 + 2.0f * v_14 * v_14) + (*dist_coeffs_8)[int(7)] * r2_14);
    float2  _S221 = _S220 + make_float2 ((*dist_coeffs_8)[int(8)] * _S220.x + (*dist_coeffs_8)[int(9)] * _S220.y, 0.0f);
    float _S222 = fx_2 * _S221.x + cx_1;
    float _S223 = fy_2 * _S221.y + cy_1;
    float2  _S224 = make_float2 (_S222, _S223);
    proj_points_1[int(3)] = _S224;
    float2  _S225 = float2 {sigmas_1->p_0[int(4)].x, sigmas_1->p_0[int(4)].y} / make_float2 (sigmas_1->p_0[int(4)].z);
    float u_15 = _S225.x;
    float v_15 = _S225.y;
    float r2_15 = u_15 * u_15 + v_15 * v_15;
    float2  _S226 = _S225 * make_float2 (1.0f + r2_15 * ((*dist_coeffs_8)[int(0)] + r2_15 * ((*dist_coeffs_8)[int(1)] + r2_15 * ((*dist_coeffs_8)[int(2)] + r2_15 * (*dist_coeffs_8)[int(3)])))) + make_float2 (_S200 * u_15 * v_15 + (*dist_coeffs_8)[int(5)] * (r2_15 + 2.0f * u_15 * u_15) + (*dist_coeffs_8)[int(6)] * r2_15, _S201 * u_15 * v_15 + (*dist_coeffs_8)[int(4)] * (r2_15 + 2.0f * v_15 * v_15) + (*dist_coeffs_8)[int(7)] * r2_15);
    float2  _S227 = _S226 + make_float2 ((*dist_coeffs_8)[int(8)] * _S226.x + (*dist_coeffs_8)[int(9)] * _S226.y, 0.0f);
    float _S228 = fx_2 * _S227.x + cx_1;
    float _S229 = fy_2 * _S227.y + cy_1;
    float2  _S230 = make_float2 (_S228, _S229);
    proj_points_1[int(4)] = _S230;
    float2  _S231 = float2 {sigmas_1->p_0[int(5)].x, sigmas_1->p_0[int(5)].y} / make_float2 (sigmas_1->p_0[int(5)].z);
    float u_16 = _S231.x;
    float v_16 = _S231.y;
    float r2_16 = u_16 * u_16 + v_16 * v_16;
    float2  _S232 = _S231 * make_float2 (1.0f + r2_16 * ((*dist_coeffs_8)[int(0)] + r2_16 * ((*dist_coeffs_8)[int(1)] + r2_16 * ((*dist_coeffs_8)[int(2)] + r2_16 * (*dist_coeffs_8)[int(3)])))) + make_float2 (_S200 * u_16 * v_16 + (*dist_coeffs_8)[int(5)] * (r2_16 + 2.0f * u_16 * u_16) + (*dist_coeffs_8)[int(6)] * r2_16, _S201 * u_16 * v_16 + (*dist_coeffs_8)[int(4)] * (r2_16 + 2.0f * v_16 * v_16) + (*dist_coeffs_8)[int(7)] * r2_16);
    float2  _S233 = _S232 + make_float2 ((*dist_coeffs_8)[int(8)] * _S232.x + (*dist_coeffs_8)[int(9)] * _S232.y, 0.0f);
    float _S234 = fx_2 * _S233.x + cx_1;
    float _S235 = fy_2 * _S233.y + cy_1;
    float2  _S236 = make_float2 (_S234, _S235);
    proj_points_1[int(5)] = _S236;
    float2  _S237 = float2 {sigmas_1->p_0[int(6)].x, sigmas_1->p_0[int(6)].y} / make_float2 (sigmas_1->p_0[int(6)].z);
    float u_17 = _S237.x;
    float v_17 = _S237.y;
    float r2_17 = u_17 * u_17 + v_17 * v_17;
    float2  _S238 = _S237 * make_float2 (1.0f + r2_17 * ((*dist_coeffs_8)[int(0)] + r2_17 * ((*dist_coeffs_8)[int(1)] + r2_17 * ((*dist_coeffs_8)[int(2)] + r2_17 * (*dist_coeffs_8)[int(3)])))) + make_float2 (_S200 * u_17 * v_17 + (*dist_coeffs_8)[int(5)] * (r2_17 + 2.0f * u_17 * u_17) + (*dist_coeffs_8)[int(6)] * r2_17, _S201 * u_17 * v_17 + (*dist_coeffs_8)[int(4)] * (r2_17 + 2.0f * v_17 * v_17) + (*dist_coeffs_8)[int(7)] * r2_17);
    float2  _S239 = _S238 + make_float2 ((*dist_coeffs_8)[int(8)] * _S238.x + (*dist_coeffs_8)[int(9)] * _S238.y, 0.0f);
    float _S240 = fx_2 * _S239.x + cx_1;
    float _S241 = fy_2 * _S239.y + cy_1;
    float2  _S242 = make_float2 (_S240, _S241);
    proj_points_1[int(6)] = _S242;
    float2  _S243 = *mean2d_2 + make_float2 (sigmas_1->w_mean_0[int(0)]) * _S206;
    *mean2d_2 = _S243;
    float2  _S244 = _S243 + make_float2 (sigmas_1->w_mean_0[int(1)]) * _S212;
    *mean2d_2 = _S244;
    float2  _S245 = _S244 + make_float2 (sigmas_1->w_mean_0[int(2)]) * _S218;
    *mean2d_2 = _S245;
    float2  _S246 = _S245 + make_float2 (sigmas_1->w_mean_0[int(3)]) * _S224;
    *mean2d_2 = _S246;
    float2  _S247 = _S246 + make_float2 (sigmas_1->w_mean_0[int(4)]) * _S230;
    *mean2d_2 = _S247;
    float2  _S248 = _S247 + make_float2 (sigmas_1->w_mean_0[int(5)]) * _S236;
    *mean2d_2 = _S248;
    float2  _S249 = _S248 + make_float2 (sigmas_1->w_mean_0[int(6)]) * _S242;
    *mean2d_2 = _S249;
    float _S250 = - _S197;
    float _S251 = - _S198;
    float2  _S252 = make_float2 (clamp_0(_S249.x, _S250, lim_x_pos_1), clamp_0(_S249.y, _S251, lim_y_pos_1));
    float2  d_7 = make_float2 (clamp_0(_S204, _S250, lim_x_pos_1), clamp_0(_S205, _S251, lim_y_pos_1)) - _S252;
    float _S253 = d_7.x;
    float _S254 = d_7.y;
    float _S255 = _S253 * _S254;
    Matrix<float, 2, 2>  _S256 = *cov2d_2 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S253 * _S253, _S255, _S255, _S254 * _S254);
    *cov2d_2 = _S256;
    float2  d_8 = make_float2 (clamp_0(_S210, _S250, lim_x_pos_1), clamp_0(_S211, _S251, lim_y_pos_1)) - _S252;
    float _S257 = d_8.x;
    float _S258 = d_8.y;
    float _S259 = _S257 * _S258;
    Matrix<float, 2, 2>  _S260 = _S256 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S257 * _S257, _S259, _S259, _S258 * _S258);
    *cov2d_2 = _S260;
    float2  d_9 = make_float2 (clamp_0(_S216, _S250, lim_x_pos_1), clamp_0(_S217, _S251, lim_y_pos_1)) - _S252;
    float _S261 = d_9.x;
    float _S262 = d_9.y;
    float _S263 = _S261 * _S262;
    Matrix<float, 2, 2>  _S264 = _S260 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S261 * _S261, _S263, _S263, _S262 * _S262);
    *cov2d_2 = _S264;
    float2  d_10 = make_float2 (clamp_0(_S222, _S250, lim_x_pos_1), clamp_0(_S223, _S251, lim_y_pos_1)) - _S252;
    float _S265 = d_10.x;
    float _S266 = d_10.y;
    float _S267 = _S265 * _S266;
    Matrix<float, 2, 2>  _S268 = _S264 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S265 * _S265, _S267, _S267, _S266 * _S266);
    *cov2d_2 = _S268;
    float2  d_11 = make_float2 (clamp_0(_S228, _S250, lim_x_pos_1), clamp_0(_S229, _S251, lim_y_pos_1)) - _S252;
    float _S269 = d_11.x;
    float _S270 = d_11.y;
    float _S271 = _S269 * _S270;
    Matrix<float, 2, 2>  _S272 = _S268 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S269 * _S269, _S271, _S271, _S270 * _S270);
    *cov2d_2 = _S272;
    float2  d_12 = make_float2 (clamp_0(_S234, _S250, lim_x_pos_1), clamp_0(_S235, _S251, lim_y_pos_1)) - _S252;
    float _S273 = d_12.x;
    float _S274 = d_12.y;
    float _S275 = _S273 * _S274;
    Matrix<float, 2, 2>  _S276 = _S272 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S273 * _S273, _S275, _S275, _S274 * _S274);
    *cov2d_2 = _S276;
    float2  d_13 = make_float2 (clamp_0(_S240, _S250, lim_x_pos_1), clamp_0(_S241, _S251, lim_y_pos_1)) - _S252;
    float _S277 = d_13.x;
    float _S278 = d_13.y;
    float _S279 = _S277 * _S278;
    *cov2d_2 = _S276 + makeMatrix<float, 2, 2> (sigmas_1->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S277 * _S277, _S279, _S279, _S278 * _S278);
    return true;
}

inline __device__ void _d_atan2_0(DiffPair_float_0 * dpy_4, DiffPair_float_0 * dpx_7, float dOut_12)
{
    DiffPair_float_0 _S280 = *dpx_7;
    float _S281 = - (*dpy_4).primal_0 / ((*dpx_7).primal_0 * (*dpx_7).primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_12;
    dpx_7->primal_0 = (*dpx_7).primal_0;
    dpx_7->differential_0 = _S281;
    float _S282 = _S280.primal_0 / (_S280.primal_0 * _S280.primal_0 + (*dpy_4).primal_0 * (*dpy_4).primal_0) * dOut_12;
    dpy_4->primal_0 = (*dpy_4).primal_0;
    dpy_4->differential_0 = _S282;
    return;
}

inline __device__ float s_primal_ctx_atan2_0(float _S283, float _S284)
{
    return (F32_atan2((_S283), (_S284)));
}

inline __device__ void s_bwd_prop_atan2_0(DiffPair_float_0 * _S285, DiffPair_float_0 * _S286, float _S287)
{
    _d_atan2_0(_S285, _S286, _S287);
    return;
}

inline __device__ void s_bwd_prop_sqrt_0(DiffPair_float_0 * _S288, float _S289)
{
    _d_sqrt_0(_S288, _S289);
    return;
}

inline __device__ void s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_8, float _s_dOut_2)
{
    float _S290 = (*dpx_8).primal_0.x;
    float _S291 = (*dpx_8).primal_0.y;
    DiffPair_float_0 _S292;
    (&_S292)->primal_0 = _S290 * _S290 + _S291 * _S291;
    (&_S292)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S292, _s_dOut_2);
    float _S293 = (*dpx_8).primal_0.y * _S292.differential_0;
    float _S294 = _S293 + _S293;
    float _S295 = (*dpx_8).primal_0.x * _S292.differential_0;
    float _S296 = _S295 + _S295;
    float2  _S297 = make_float2 (0.0f);
    *&((&_S297)->y) = _S294;
    *&((&_S297)->x) = _S296;
    dpx_8->primal_0 = (*dpx_8).primal_0;
    dpx_8->differential_0 = _S297;
    return;
}

inline __device__ void s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S298, float _S299)
{
    s_bwd_prop_length_impl_0(_S298, _S299);
    return;
}

inline __device__ bool fisheye_proj_3dgs_0(float3  mean3d_1, Matrix<float, 3, 3>  cov3d_1, float4  intrins_2, FixedArray<float, 10>  * dist_coeffs_9, Matrix<float, 2, 2>  * cov2d_3, float2  * mean2d_3)
{
    float k_0;
    float2  _S300;
    float _S301;
    float _S302;
    bool _S303;
    for(;;)
    {
        float2  _S304 = float2 {mean3d_1.x, mean3d_1.y};
        _S300 = _S304;
        float r_8 = length_0(_S304);
        _S301 = r_8;
        float _S305 = mean3d_1.z;
        _S302 = _S305;
        float theta_4 = (F32_atan2((r_8), (_S305)));
        if(theta_4 < 0.00100000004749745f)
        {
            k_0 = (1.0f - theta_4 * theta_4 / 3.0f) / _S305;
        }
        else
        {
            k_0 = theta_4 / r_8;
        }
        float2  _S306 = _S304 * make_float2 (k_0);
        *mean2d_3 = _S306;
        Matrix<float, 2, 2>  _S307 = camera_distortion_jac_0(_S306, dist_coeffs_9);
        bool _S308 = !((F32_min((determinant_0(_S307)), ((F32_min((_S307.rows[int(0)].x), (_S307.rows[int(1)].y)))))) > 0.0f);
        _S303 = _S308;
        if(_S308)
        {
            break;
        }
        float u_18 = (*mean2d_3).x;
        float v_18 = (*mean2d_3).y;
        float r2_18 = u_18 * u_18 + v_18 * v_18;
        float2  _S309 = *mean2d_3 * make_float2 (1.0f + r2_18 * ((*dist_coeffs_9)[int(0)] + r2_18 * ((*dist_coeffs_9)[int(1)] + r2_18 * ((*dist_coeffs_9)[int(2)] + r2_18 * (*dist_coeffs_9)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_9)[int(4)] * u_18 * v_18 + (*dist_coeffs_9)[int(5)] * (r2_18 + 2.0f * u_18 * u_18) + (*dist_coeffs_9)[int(6)] * r2_18, 2.0f * (*dist_coeffs_9)[int(5)] * u_18 * v_18 + (*dist_coeffs_9)[int(4)] * (r2_18 + 2.0f * v_18 * v_18) + (*dist_coeffs_9)[int(7)] * r2_18);
        float2  _S310 = _S309 + make_float2 ((*dist_coeffs_9)[int(8)] * _S309.x + (*dist_coeffs_9)[int(9)] * _S309.y, 0.0f);
        *mean2d_3 = make_float2 (intrins_2.x * _S310.x + intrins_2.z, intrins_2.y * _S310.y + intrins_2.w);
        break;
    }
    if(!!_S303)
    {
        return false;
    }
    Matrix<float, 2, 3>  J_1;
    float2  _S311 = make_float2 (0.0f);
    float2  seed_0 = _S311;
    *&((&seed_0)->x) = 1.0f;
    float2  _S312 = seed_0;
    float _S313 = s_primal_ctx_atan2_0(_S301, _S302);
    bool _S314 = _S313 < 0.00100000004749745f;
    float _S315;
    float _S316;
    float _S317;
    if(_S314)
    {
        float _S318 = 1.0f - _S313 * _S313 / 3.0f;
        float _S319 = _S302 * _S302;
        k_0 = _S318 / _S302;
        _S315 = 0.0f;
        _S316 = _S319;
        _S317 = _S318;
    }
    else
    {
        float _S320 = _S301 * _S301;
        k_0 = _S313 / _S301;
        _S315 = _S320;
        _S316 = 0.0f;
        _S317 = 0.0f;
    }
    float2  _S321 = make_float2 (k_0);
    float2  _S322 = _S300 * make_float2 (k_0);
    float u_19 = _S322.x;
    float v_19 = _S322.y;
    float r2_19 = u_19 * u_19 + v_19 * v_19;
    float _S323 = (*dist_coeffs_9)[int(2)] + r2_19 * (*dist_coeffs_9)[int(3)];
    float _S324 = (*dist_coeffs_9)[int(1)] + r2_19 * _S323;
    float _S325 = (*dist_coeffs_9)[int(0)] + r2_19 * _S324;
    float _S326 = 2.0f * (*dist_coeffs_9)[int(4)];
    float _S327 = 2.0f * (*dist_coeffs_9)[int(5)];
    float fx_3 = intrins_2.x;
    float fy_3 = intrins_2.y;
    float _S328 = fx_3 * _S312.x;
    float2  _S329 = make_float2 (_S328, fy_3 * _S312.y) + make_float2 ((*dist_coeffs_9)[int(8)] * _S328, (*dist_coeffs_9)[int(9)] * _S328);
    float2  _S330 = _S322 * _S329;
    float _S331 = (*dist_coeffs_9)[int(4)] * _S329.y;
    float _S332 = (*dist_coeffs_9)[int(5)] * _S329.x;
    float _S333 = _S330.x + _S330.y;
    float _S334 = r2_19 * _S333;
    float _S335 = r2_19 * _S334;
    float _S336 = (*dist_coeffs_9)[int(7)] * _S329.y + _S331 + (*dist_coeffs_9)[int(6)] * _S329.x + _S332 + _S325 * _S333 + _S324 * _S334 + _S323 * _S335 + (*dist_coeffs_9)[int(3)] * (r2_19 * _S335);
    float _S337 = v_19 * _S336;
    float _S338 = u_19 * _S336;
    float2  _S339 = make_float2 (1.0f + r2_19 * _S325) * _S329 + make_float2 (_S327 * (v_19 * _S329.y) + 2.0f * u_19 * _S332 + 2.0f * (u_19 * _S332) + _S326 * (v_19 * _S329.x) + _S338 + _S338, 2.0f * v_19 * _S331 + 2.0f * (v_19 * _S331) + _S327 * u_19 * _S329.y + _S326 * u_19 * _S329.x + _S337 + _S337);
    float2  _S340 = _S300 * _S339;
    float2  _S341 = _S321 * _S339;
    float _S342 = _S340.x + _S340.y;
    if(_S314)
    {
        float _S343 = _S342 / _S316;
        float _S344 = _S317 * - _S343;
        float _S345 = _S313 * (0.3333333432674408f * - (_S302 * _S343));
        k_0 = _S345 + _S345;
        _S315 = _S344;
        _S316 = 0.0f;
    }
    else
    {
        float _S346 = _S342 / _S315;
        float _S347 = _S313 * - _S346;
        k_0 = _S301 * _S346;
        _S315 = 0.0f;
        _S316 = _S347;
    }
    float _S348 = _S301;
    DiffPair_float_0 _S349;
    (&_S349)->primal_0 = _S301;
    (&_S349)->differential_0 = 0.0f;
    float _S350 = _S302;
    DiffPair_float_0 _S351;
    (&_S351)->primal_0 = _S302;
    (&_S351)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S349, &_S351, k_0);
    float _S352 = _S351.differential_0 + _S315;
    float _S353 = _S349.differential_0 + _S316;
    float2  _S354 = make_float2 (0.0f);
    float2  _S355 = _S300;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S356;
    (&_S356)->primal_0 = _S300;
    (&_S356)->differential_0 = _S354;
    s_bwd_length_impl_0(&_S356, _S353);
    float2  _S357 = _S356.differential_0 + _S341;
    float3  _S358 = make_float3 (_S357.x, _S357.y, _S352);
    J_1[int(0)] = _S358;
    float2  seed_1 = _S311;
    *&((&seed_1)->y) = 1.0f;
    float2  _S359 = seed_1;
    if(_S314)
    {
        float _S360 = 1.0f - _S313 * _S313 / 3.0f;
        float _S361 = _S302 * _S302;
        k_0 = _S360 / _S302;
        _S315 = 0.0f;
        _S316 = _S361;
        _S317 = _S360;
    }
    else
    {
        float _S362 = _S301 * _S301;
        k_0 = _S313 / _S301;
        _S315 = _S362;
        _S316 = 0.0f;
        _S317 = 0.0f;
    }
    float2  _S363 = make_float2 (k_0);
    float2  _S364 = _S300 * make_float2 (k_0);
    float u_20 = _S364.x;
    float v_20 = _S364.y;
    float r2_20 = u_20 * u_20 + v_20 * v_20;
    float _S365 = (*dist_coeffs_9)[int(2)] + r2_20 * (*dist_coeffs_9)[int(3)];
    float _S366 = (*dist_coeffs_9)[int(1)] + r2_20 * _S365;
    float _S367 = (*dist_coeffs_9)[int(0)] + r2_20 * _S366;
    float _S368 = fx_3 * _S359.x;
    float2  _S369 = make_float2 (_S368, fy_3 * _S359.y) + make_float2 ((*dist_coeffs_9)[int(8)] * _S368, (*dist_coeffs_9)[int(9)] * _S368);
    float2  _S370 = _S364 * _S369;
    float _S371 = (*dist_coeffs_9)[int(4)] * _S369.y;
    float _S372 = (*dist_coeffs_9)[int(5)] * _S369.x;
    float _S373 = _S370.x + _S370.y;
    float _S374 = r2_20 * _S373;
    float _S375 = r2_20 * _S374;
    float _S376 = (*dist_coeffs_9)[int(7)] * _S369.y + _S371 + (*dist_coeffs_9)[int(6)] * _S369.x + _S372 + _S367 * _S373 + _S366 * _S374 + _S365 * _S375 + (*dist_coeffs_9)[int(3)] * (r2_20 * _S375);
    float _S377 = v_20 * _S376;
    float _S378 = u_20 * _S376;
    float2  _S379 = make_float2 (1.0f + r2_20 * _S367) * _S369 + make_float2 (_S327 * (v_20 * _S369.y) + 2.0f * u_20 * _S372 + 2.0f * (u_20 * _S372) + _S326 * (v_20 * _S369.x) + _S378 + _S378, 2.0f * v_20 * _S371 + 2.0f * (v_20 * _S371) + _S327 * u_20 * _S369.y + _S326 * u_20 * _S369.x + _S377 + _S377);
    float2  _S380 = _S300 * _S379;
    float2  _S381 = _S363 * _S379;
    float _S382 = _S380.x + _S380.y;
    if(_S314)
    {
        float _S383 = _S382 / _S316;
        float _S384 = _S317 * - _S383;
        float _S385 = _S313 * (0.3333333432674408f * - (_S302 * _S383));
        k_0 = _S385 + _S385;
        _S315 = _S384;
        _S316 = 0.0f;
    }
    else
    {
        float _S386 = _S382 / _S315;
        float _S387 = _S313 * - _S386;
        k_0 = _S301 * _S386;
        _S315 = 0.0f;
        _S316 = _S387;
    }
    DiffPair_float_0 _S388;
    (&_S388)->primal_0 = _S348;
    (&_S388)->differential_0 = 0.0f;
    DiffPair_float_0 _S389;
    (&_S389)->primal_0 = _S350;
    (&_S389)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S388, &_S389, k_0);
    float _S390 = _S389.differential_0 + _S315;
    float _S391 = _S388.differential_0 + _S316;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S392;
    (&_S392)->primal_0 = _S355;
    (&_S392)->differential_0 = _S354;
    s_bwd_length_impl_0(&_S392, _S391);
    float2  _S393 = _S392.differential_0 + _S381;
    float3  _S394 = make_float3 (_S393.x, _S393.y, _S390);
    J_1[int(1)] = _S394;
    *cov2d_3 = mul_6(mul_5(J_1, cov3d_1), transpose_1(J_1));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_1(float3  mean3d_2, Matrix<float, 3, 3>  cov3d_2, float4  intrins_3, FixedArray<float, 10>  * dist_coeffs_10, Matrix<float, 2, 2>  * cov2d_4, float2  * mean2d_4)
{
    float2  _S395 = float2 {mean3d_2.x, mean3d_2.y};
    float r_9 = length_0(_S395);
    float _S396 = mean3d_2.z;
    float theta_5 = (F32_atan2((r_9), (_S396)));
    float k_1;
    if(theta_5 < 0.00100000004749745f)
    {
        k_1 = (1.0f - theta_5 * theta_5 / 3.0f) / _S396;
    }
    else
    {
        k_1 = theta_5 / r_9;
    }
    float2  _S397 = _S395 * make_float2 (k_1);
    float u_21 = _S397.x;
    float v_21 = _S397.y;
    float r2_21 = u_21 * u_21 + v_21 * v_21;
    float _S398 = 2.0f * (*dist_coeffs_10)[int(4)];
    float _S399 = 2.0f * (*dist_coeffs_10)[int(5)];
    float2  _S400 = _S397 * make_float2 (1.0f + r2_21 * ((*dist_coeffs_10)[int(0)] + r2_21 * ((*dist_coeffs_10)[int(1)] + r2_21 * ((*dist_coeffs_10)[int(2)] + r2_21 * (*dist_coeffs_10)[int(3)])))) + make_float2 (_S398 * u_21 * v_21 + (*dist_coeffs_10)[int(5)] * (r2_21 + 2.0f * u_21 * u_21) + (*dist_coeffs_10)[int(6)] * r2_21, _S399 * u_21 * v_21 + (*dist_coeffs_10)[int(4)] * (r2_21 + 2.0f * v_21 * v_21) + (*dist_coeffs_10)[int(7)] * r2_21);
    float2  _S401 = _S400 + make_float2 ((*dist_coeffs_10)[int(8)] * _S400.x + (*dist_coeffs_10)[int(9)] * _S400.y, 0.0f);
    float fx_4 = intrins_3.x;
    float fy_4 = intrins_3.y;
    *mean2d_4 = make_float2 (fx_4 * _S401.x + intrins_3.z, fy_4 * _S401.y + intrins_3.w);
    Matrix<float, 2, 3>  J_2;
    float2  _S402 = make_float2 (0.0f);
    float2  seed_2 = _S402;
    *&((&seed_2)->x) = 1.0f;
    float2  _S403 = seed_2;
    float _S404 = s_primal_ctx_atan2_0(r_9, _S396);
    bool _S405 = _S404 < 0.00100000004749745f;
    float _S406;
    float _S407;
    float _S408;
    if(_S405)
    {
        float _S409 = 1.0f - _S404 * _S404 / 3.0f;
        float _S410 = _S396 * _S396;
        k_1 = _S409 / _S396;
        _S406 = 0.0f;
        _S407 = _S410;
        _S408 = _S409;
    }
    else
    {
        float _S411 = r_9 * r_9;
        k_1 = _S404 / r_9;
        _S406 = _S411;
        _S407 = 0.0f;
        _S408 = 0.0f;
    }
    float2  _S412 = make_float2 (k_1);
    float2  _S413 = _S395 * make_float2 (k_1);
    float u_22 = _S413.x;
    float v_22 = _S413.y;
    float r2_22 = u_22 * u_22 + v_22 * v_22;
    float _S414 = (*dist_coeffs_10)[int(2)] + r2_22 * (*dist_coeffs_10)[int(3)];
    float _S415 = (*dist_coeffs_10)[int(1)] + r2_22 * _S414;
    float _S416 = (*dist_coeffs_10)[int(0)] + r2_22 * _S415;
    float _S417 = fx_4 * _S403.x;
    float2  _S418 = make_float2 (_S417, fy_4 * _S403.y) + make_float2 ((*dist_coeffs_10)[int(8)] * _S417, (*dist_coeffs_10)[int(9)] * _S417);
    float2  _S419 = _S413 * _S418;
    float _S420 = (*dist_coeffs_10)[int(4)] * _S418.y;
    float _S421 = (*dist_coeffs_10)[int(5)] * _S418.x;
    float _S422 = _S419.x + _S419.y;
    float _S423 = r2_22 * _S422;
    float _S424 = r2_22 * _S423;
    float _S425 = (*dist_coeffs_10)[int(7)] * _S418.y + _S420 + (*dist_coeffs_10)[int(6)] * _S418.x + _S421 + _S416 * _S422 + _S415 * _S423 + _S414 * _S424 + (*dist_coeffs_10)[int(3)] * (r2_22 * _S424);
    float _S426 = v_22 * _S425;
    float _S427 = u_22 * _S425;
    float2  _S428 = make_float2 (1.0f + r2_22 * _S416) * _S418 + make_float2 (_S399 * (v_22 * _S418.y) + 2.0f * u_22 * _S421 + 2.0f * (u_22 * _S421) + _S398 * (v_22 * _S418.x) + _S427 + _S427, 2.0f * v_22 * _S420 + 2.0f * (v_22 * _S420) + _S399 * u_22 * _S418.y + _S398 * u_22 * _S418.x + _S426 + _S426);
    float2  _S429 = _S395 * _S428;
    float2  _S430 = _S412 * _S428;
    float _S431 = _S429.x + _S429.y;
    if(_S405)
    {
        float _S432 = _S431 / _S407;
        float _S433 = _S408 * - _S432;
        float _S434 = _S404 * (0.3333333432674408f * - (_S396 * _S432));
        k_1 = _S434 + _S434;
        _S406 = _S433;
        _S407 = 0.0f;
    }
    else
    {
        float _S435 = _S431 / _S406;
        float _S436 = _S404 * - _S435;
        k_1 = r_9 * _S435;
        _S406 = 0.0f;
        _S407 = _S436;
    }
    DiffPair_float_0 _S437;
    (&_S437)->primal_0 = r_9;
    (&_S437)->differential_0 = 0.0f;
    DiffPair_float_0 _S438;
    (&_S438)->primal_0 = _S396;
    (&_S438)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S437, &_S438, k_1);
    float _S439 = _S438.differential_0 + _S406;
    float _S440 = _S437.differential_0 + _S407;
    float2  _S441 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S442;
    (&_S442)->primal_0 = _S395;
    (&_S442)->differential_0 = _S441;
    s_bwd_length_impl_0(&_S442, _S440);
    float2  _S443 = _S442.differential_0 + _S430;
    float3  _S444 = make_float3 (_S443.x, _S443.y, _S439);
    J_2[int(0)] = _S444;
    float2  seed_3 = _S402;
    *&((&seed_3)->y) = 1.0f;
    float2  _S445 = seed_3;
    if(_S405)
    {
        float _S446 = 1.0f - _S404 * _S404 / 3.0f;
        float _S447 = _S396 * _S396;
        k_1 = _S446 / _S396;
        _S406 = 0.0f;
        _S407 = _S447;
        _S408 = _S446;
    }
    else
    {
        float _S448 = r_9 * r_9;
        k_1 = _S404 / r_9;
        _S406 = _S448;
        _S407 = 0.0f;
        _S408 = 0.0f;
    }
    float2  _S449 = make_float2 (k_1);
    float2  _S450 = _S395 * make_float2 (k_1);
    float u_23 = _S450.x;
    float v_23 = _S450.y;
    float r2_23 = u_23 * u_23 + v_23 * v_23;
    float _S451 = (*dist_coeffs_10)[int(2)] + r2_23 * (*dist_coeffs_10)[int(3)];
    float _S452 = (*dist_coeffs_10)[int(1)] + r2_23 * _S451;
    float _S453 = (*dist_coeffs_10)[int(0)] + r2_23 * _S452;
    float _S454 = fx_4 * _S445.x;
    float2  _S455 = make_float2 (_S454, fy_4 * _S445.y) + make_float2 ((*dist_coeffs_10)[int(8)] * _S454, (*dist_coeffs_10)[int(9)] * _S454);
    float2  _S456 = _S450 * _S455;
    float _S457 = (*dist_coeffs_10)[int(4)] * _S455.y;
    float _S458 = (*dist_coeffs_10)[int(5)] * _S455.x;
    float _S459 = _S456.x + _S456.y;
    float _S460 = r2_23 * _S459;
    float _S461 = r2_23 * _S460;
    float _S462 = (*dist_coeffs_10)[int(7)] * _S455.y + _S457 + (*dist_coeffs_10)[int(6)] * _S455.x + _S458 + _S453 * _S459 + _S452 * _S460 + _S451 * _S461 + (*dist_coeffs_10)[int(3)] * (r2_23 * _S461);
    float _S463 = v_23 * _S462;
    float _S464 = u_23 * _S462;
    float2  _S465 = make_float2 (1.0f + r2_23 * _S453) * _S455 + make_float2 (_S399 * (v_23 * _S455.y) + 2.0f * u_23 * _S458 + 2.0f * (u_23 * _S458) + _S398 * (v_23 * _S455.x) + _S464 + _S464, 2.0f * v_23 * _S457 + 2.0f * (v_23 * _S457) + _S399 * u_23 * _S455.y + _S398 * u_23 * _S455.x + _S463 + _S463);
    float2  _S466 = _S395 * _S465;
    float2  _S467 = _S449 * _S465;
    float _S468 = _S466.x + _S466.y;
    if(_S405)
    {
        float _S469 = _S468 / _S407;
        float _S470 = _S408 * - _S469;
        float _S471 = _S404 * (0.3333333432674408f * - (_S396 * _S469));
        k_1 = _S471 + _S471;
        _S406 = _S470;
        _S407 = 0.0f;
    }
    else
    {
        float _S472 = _S468 / _S406;
        float _S473 = _S404 * - _S472;
        k_1 = r_9 * _S472;
        _S406 = 0.0f;
        _S407 = _S473;
    }
    DiffPair_float_0 _S474;
    (&_S474)->primal_0 = r_9;
    (&_S474)->differential_0 = 0.0f;
    DiffPair_float_0 _S475;
    (&_S475)->primal_0 = _S396;
    (&_S475)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S474, &_S475, k_1);
    float _S476 = _S475.differential_0 + _S406;
    float _S477 = _S474.differential_0 + _S407;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S478;
    (&_S478)->primal_0 = _S395;
    (&_S478)->differential_0 = _S441;
    s_bwd_length_impl_0(&_S478, _S477);
    float2  _S479 = _S478.differential_0 + _S467;
    float3  _S480 = make_float3 (_S479.x, _S479.y, _S476);
    J_2[int(1)] = _S480;
    *cov2d_4 = mul_6(mul_5(J_2, cov3d_2), transpose_1(J_2));
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_0(SigmaPoints_0 * sigmas_2, float4  intrins_4, FixedArray<float, 10>  * dist_coeffs_11, Matrix<float, 2, 2>  * cov2d_5, float2  * mean2d_5)
{
    float2  * _S481;
    bool _S482;
    float2  * _S483;
    bool _S484;
    float2  * _S485;
    bool _S486;
    bool _S487;
    float2  * _S488;
    bool _S489;
    float2  * _S490;
    bool _S491;
    float2  * _S492;
    bool _S493;
    bool _S494;
    float2  * _S495;
    bool _S496;
    bool _S497;
    int2  _S498 = make_int2 (int(0));
    float2  _S499 = make_float2 ((float)_S498.x, (float)_S498.y);
    *mean2d_5 = _S499;
    *cov2d_5 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_2;
    for(;;)
    {
        float k_2;
        _S481 = &proj_points_2[int(0)];
        for(;;)
        {
            float2  _S500 = float2 {sigmas_2->p_0[int(0)].x, sigmas_2->p_0[int(0)].y};
            float r_10 = length_0(_S500);
            float _S501 = sigmas_2->p_0[int(0)].z;
            float theta_6 = (F32_atan2((r_10), (_S501)));
            if(theta_6 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_6 * theta_6 / 3.0f) / _S501;
            }
            else
            {
                k_2 = theta_6 / r_10;
            }
            float2  _S502 = _S500 * make_float2 (k_2);
            proj_points_2[int(0)] = _S502;
            Matrix<float, 2, 2>  _S503 = camera_distortion_jac_0(_S502, dist_coeffs_11);
            bool _S504 = !((F32_min((determinant_0(_S503)), ((F32_min((_S503.rows[int(0)].x), (_S503.rows[int(1)].y)))))) > 0.0f);
            _S482 = _S504;
            if(_S504)
            {
                break;
            }
            float u_24 = proj_points_2[int(0)].x;
            float v_24 = proj_points_2[int(0)].y;
            float r2_24 = u_24 * u_24 + v_24 * v_24;
            float2  _S505 = proj_points_2[int(0)] * make_float2 (1.0f + r2_24 * ((*dist_coeffs_11)[int(0)] + r2_24 * ((*dist_coeffs_11)[int(1)] + r2_24 * ((*dist_coeffs_11)[int(2)] + r2_24 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_24 * v_24 + (*dist_coeffs_11)[int(5)] * (r2_24 + 2.0f * u_24 * u_24) + (*dist_coeffs_11)[int(6)] * r2_24, 2.0f * (*dist_coeffs_11)[int(5)] * u_24 * v_24 + (*dist_coeffs_11)[int(4)] * (r2_24 + 2.0f * v_24 * v_24) + (*dist_coeffs_11)[int(7)] * r2_24);
            float2  _S506 = _S505 + make_float2 ((*dist_coeffs_11)[int(8)] * _S505.x + (*dist_coeffs_11)[int(9)] * _S505.y, 0.0f);
            proj_points_2[int(0)] = make_float2 (intrins_4.x * _S506.x + intrins_4.z, intrins_4.y * _S506.y + intrins_4.w);
            break;
        }
        bool all_valid_4 = true & (!_S482);
        _S483 = &proj_points_2[int(1)];
        for(;;)
        {
            float2  _S507 = float2 {sigmas_2->p_0[int(1)].x, sigmas_2->p_0[int(1)].y};
            float r_11 = length_0(_S507);
            float _S508 = sigmas_2->p_0[int(1)].z;
            float theta_7 = (F32_atan2((r_11), (_S508)));
            if(theta_7 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_7 * theta_7 / 3.0f) / _S508;
            }
            else
            {
                k_2 = theta_7 / r_11;
            }
            float2  _S509 = _S507 * make_float2 (k_2);
            proj_points_2[int(1)] = _S509;
            Matrix<float, 2, 2>  _S510 = camera_distortion_jac_0(_S509, dist_coeffs_11);
            bool _S511 = !((F32_min((determinant_0(_S510)), ((F32_min((_S510.rows[int(0)].x), (_S510.rows[int(1)].y)))))) > 0.0f);
            _S484 = _S511;
            if(_S511)
            {
                break;
            }
            float u_25 = proj_points_2[int(1)].x;
            float v_25 = proj_points_2[int(1)].y;
            float r2_25 = u_25 * u_25 + v_25 * v_25;
            float2  _S512 = proj_points_2[int(1)] * make_float2 (1.0f + r2_25 * ((*dist_coeffs_11)[int(0)] + r2_25 * ((*dist_coeffs_11)[int(1)] + r2_25 * ((*dist_coeffs_11)[int(2)] + r2_25 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_25 * v_25 + (*dist_coeffs_11)[int(5)] * (r2_25 + 2.0f * u_25 * u_25) + (*dist_coeffs_11)[int(6)] * r2_25, 2.0f * (*dist_coeffs_11)[int(5)] * u_25 * v_25 + (*dist_coeffs_11)[int(4)] * (r2_25 + 2.0f * v_25 * v_25) + (*dist_coeffs_11)[int(7)] * r2_25);
            float2  _S513 = _S512 + make_float2 ((*dist_coeffs_11)[int(8)] * _S512.x + (*dist_coeffs_11)[int(9)] * _S512.y, 0.0f);
            proj_points_2[int(1)] = make_float2 (intrins_4.x * _S513.x + intrins_4.z, intrins_4.y * _S513.y + intrins_4.w);
            break;
        }
        bool all_valid_5 = all_valid_4 & (!_S484);
        for(;;)
        {
            _S485 = &proj_points_2[int(2)];
            for(;;)
            {
                float2  _S514 = float2 {sigmas_2->p_0[int(2)].x, sigmas_2->p_0[int(2)].y};
                float r_12 = length_0(_S514);
                float _S515 = sigmas_2->p_0[int(2)].z;
                float theta_8 = (F32_atan2((r_12), (_S515)));
                if(theta_8 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_8 * theta_8 / 3.0f) / _S515;
                }
                else
                {
                    k_2 = theta_8 / r_12;
                }
                float2  _S516 = _S514 * make_float2 (k_2);
                proj_points_2[int(2)] = _S516;
                Matrix<float, 2, 2>  _S517 = camera_distortion_jac_0(_S516, dist_coeffs_11);
                bool _S518 = !((F32_min((determinant_0(_S517)), ((F32_min((_S517.rows[int(0)].x), (_S517.rows[int(1)].y)))))) > 0.0f);
                _S486 = _S518;
                if(_S518)
                {
                    break;
                }
                float u_26 = proj_points_2[int(2)].x;
                float v_26 = proj_points_2[int(2)].y;
                float r2_26 = u_26 * u_26 + v_26 * v_26;
                float2  _S519 = proj_points_2[int(2)] * make_float2 (1.0f + r2_26 * ((*dist_coeffs_11)[int(0)] + r2_26 * ((*dist_coeffs_11)[int(1)] + r2_26 * ((*dist_coeffs_11)[int(2)] + r2_26 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_26 * v_26 + (*dist_coeffs_11)[int(5)] * (r2_26 + 2.0f * u_26 * u_26) + (*dist_coeffs_11)[int(6)] * r2_26, 2.0f * (*dist_coeffs_11)[int(5)] * u_26 * v_26 + (*dist_coeffs_11)[int(4)] * (r2_26 + 2.0f * v_26 * v_26) + (*dist_coeffs_11)[int(7)] * r2_26);
                float2  _S520 = _S519 + make_float2 ((*dist_coeffs_11)[int(8)] * _S519.x + (*dist_coeffs_11)[int(9)] * _S519.y, 0.0f);
                proj_points_2[int(2)] = make_float2 (intrins_4.x * _S520.x + intrins_4.z, intrins_4.y * _S520.y + intrins_4.w);
                break;
            }
            _S487 = all_valid_5 & (!_S486);
            break;
        }
        _S488 = &proj_points_2[int(3)];
        for(;;)
        {
            float2  _S521 = float2 {sigmas_2->p_0[int(3)].x, sigmas_2->p_0[int(3)].y};
            float r_13 = length_0(_S521);
            float _S522 = sigmas_2->p_0[int(3)].z;
            float theta_9 = (F32_atan2((r_13), (_S522)));
            if(theta_9 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_9 * theta_9 / 3.0f) / _S522;
            }
            else
            {
                k_2 = theta_9 / r_13;
            }
            float2  _S523 = _S521 * make_float2 (k_2);
            proj_points_2[int(3)] = _S523;
            Matrix<float, 2, 2>  _S524 = camera_distortion_jac_0(_S523, dist_coeffs_11);
            bool _S525 = !((F32_min((determinant_0(_S524)), ((F32_min((_S524.rows[int(0)].x), (_S524.rows[int(1)].y)))))) > 0.0f);
            _S489 = _S525;
            if(_S525)
            {
                break;
            }
            float u_27 = proj_points_2[int(3)].x;
            float v_27 = proj_points_2[int(3)].y;
            float r2_27 = u_27 * u_27 + v_27 * v_27;
            float2  _S526 = proj_points_2[int(3)] * make_float2 (1.0f + r2_27 * ((*dist_coeffs_11)[int(0)] + r2_27 * ((*dist_coeffs_11)[int(1)] + r2_27 * ((*dist_coeffs_11)[int(2)] + r2_27 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_27 * v_27 + (*dist_coeffs_11)[int(5)] * (r2_27 + 2.0f * u_27 * u_27) + (*dist_coeffs_11)[int(6)] * r2_27, 2.0f * (*dist_coeffs_11)[int(5)] * u_27 * v_27 + (*dist_coeffs_11)[int(4)] * (r2_27 + 2.0f * v_27 * v_27) + (*dist_coeffs_11)[int(7)] * r2_27);
            float2  _S527 = _S526 + make_float2 ((*dist_coeffs_11)[int(8)] * _S526.x + (*dist_coeffs_11)[int(9)] * _S526.y, 0.0f);
            proj_points_2[int(3)] = make_float2 (intrins_4.x * _S527.x + intrins_4.z, intrins_4.y * _S527.y + intrins_4.w);
            break;
        }
        bool all_valid_6 = _S487 & (!_S489);
        _S490 = &proj_points_2[int(4)];
        for(;;)
        {
            float2  _S528 = float2 {sigmas_2->p_0[int(4)].x, sigmas_2->p_0[int(4)].y};
            float r_14 = length_0(_S528);
            float _S529 = sigmas_2->p_0[int(4)].z;
            float theta_10 = (F32_atan2((r_14), (_S529)));
            if(theta_10 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_10 * theta_10 / 3.0f) / _S529;
            }
            else
            {
                k_2 = theta_10 / r_14;
            }
            float2  _S530 = _S528 * make_float2 (k_2);
            proj_points_2[int(4)] = _S530;
            Matrix<float, 2, 2>  _S531 = camera_distortion_jac_0(_S530, dist_coeffs_11);
            bool _S532 = !((F32_min((determinant_0(_S531)), ((F32_min((_S531.rows[int(0)].x), (_S531.rows[int(1)].y)))))) > 0.0f);
            _S491 = _S532;
            if(_S532)
            {
                break;
            }
            float u_28 = proj_points_2[int(4)].x;
            float v_28 = proj_points_2[int(4)].y;
            float r2_28 = u_28 * u_28 + v_28 * v_28;
            float2  _S533 = proj_points_2[int(4)] * make_float2 (1.0f + r2_28 * ((*dist_coeffs_11)[int(0)] + r2_28 * ((*dist_coeffs_11)[int(1)] + r2_28 * ((*dist_coeffs_11)[int(2)] + r2_28 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_28 * v_28 + (*dist_coeffs_11)[int(5)] * (r2_28 + 2.0f * u_28 * u_28) + (*dist_coeffs_11)[int(6)] * r2_28, 2.0f * (*dist_coeffs_11)[int(5)] * u_28 * v_28 + (*dist_coeffs_11)[int(4)] * (r2_28 + 2.0f * v_28 * v_28) + (*dist_coeffs_11)[int(7)] * r2_28);
            float2  _S534 = _S533 + make_float2 ((*dist_coeffs_11)[int(8)] * _S533.x + (*dist_coeffs_11)[int(9)] * _S533.y, 0.0f);
            proj_points_2[int(4)] = make_float2 (intrins_4.x * _S534.x + intrins_4.z, intrins_4.y * _S534.y + intrins_4.w);
            break;
        }
        bool all_valid_7 = all_valid_6 & (!_S491);
        for(;;)
        {
            _S492 = &proj_points_2[int(5)];
            for(;;)
            {
                float2  _S535 = float2 {sigmas_2->p_0[int(5)].x, sigmas_2->p_0[int(5)].y};
                float r_15 = length_0(_S535);
                float _S536 = sigmas_2->p_0[int(5)].z;
                float theta_11 = (F32_atan2((r_15), (_S536)));
                if(theta_11 < 0.00100000004749745f)
                {
                    k_2 = (1.0f - theta_11 * theta_11 / 3.0f) / _S536;
                }
                else
                {
                    k_2 = theta_11 / r_15;
                }
                float2  _S537 = _S535 * make_float2 (k_2);
                proj_points_2[int(5)] = _S537;
                Matrix<float, 2, 2>  _S538 = camera_distortion_jac_0(_S537, dist_coeffs_11);
                bool _S539 = !((F32_min((determinant_0(_S538)), ((F32_min((_S538.rows[int(0)].x), (_S538.rows[int(1)].y)))))) > 0.0f);
                _S493 = _S539;
                if(_S539)
                {
                    break;
                }
                float u_29 = proj_points_2[int(5)].x;
                float v_29 = proj_points_2[int(5)].y;
                float r2_29 = u_29 * u_29 + v_29 * v_29;
                float2  _S540 = proj_points_2[int(5)] * make_float2 (1.0f + r2_29 * ((*dist_coeffs_11)[int(0)] + r2_29 * ((*dist_coeffs_11)[int(1)] + r2_29 * ((*dist_coeffs_11)[int(2)] + r2_29 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_29 * v_29 + (*dist_coeffs_11)[int(5)] * (r2_29 + 2.0f * u_29 * u_29) + (*dist_coeffs_11)[int(6)] * r2_29, 2.0f * (*dist_coeffs_11)[int(5)] * u_29 * v_29 + (*dist_coeffs_11)[int(4)] * (r2_29 + 2.0f * v_29 * v_29) + (*dist_coeffs_11)[int(7)] * r2_29);
                float2  _S541 = _S540 + make_float2 ((*dist_coeffs_11)[int(8)] * _S540.x + (*dist_coeffs_11)[int(9)] * _S540.y, 0.0f);
                proj_points_2[int(5)] = make_float2 (intrins_4.x * _S541.x + intrins_4.z, intrins_4.y * _S541.y + intrins_4.w);
                break;
            }
            _S494 = all_valid_7 & (!_S493);
            break;
        }
        _S495 = &proj_points_2[int(6)];
        for(;;)
        {
            float2  _S542 = float2 {sigmas_2->p_0[int(6)].x, sigmas_2->p_0[int(6)].y};
            float r_16 = length_0(_S542);
            float _S543 = sigmas_2->p_0[int(6)].z;
            float theta_12 = (F32_atan2((r_16), (_S543)));
            if(theta_12 < 0.00100000004749745f)
            {
                k_2 = (1.0f - theta_12 * theta_12 / 3.0f) / _S543;
            }
            else
            {
                k_2 = theta_12 / r_16;
            }
            float2  _S544 = _S542 * make_float2 (k_2);
            proj_points_2[int(6)] = _S544;
            Matrix<float, 2, 2>  _S545 = camera_distortion_jac_0(_S544, dist_coeffs_11);
            bool _S546 = !((F32_min((determinant_0(_S545)), ((F32_min((_S545.rows[int(0)].x), (_S545.rows[int(1)].y)))))) > 0.0f);
            _S496 = _S546;
            if(_S546)
            {
                break;
            }
            float u_30 = proj_points_2[int(6)].x;
            float v_30 = proj_points_2[int(6)].y;
            float r2_30 = u_30 * u_30 + v_30 * v_30;
            float2  _S547 = proj_points_2[int(6)] * make_float2 (1.0f + r2_30 * ((*dist_coeffs_11)[int(0)] + r2_30 * ((*dist_coeffs_11)[int(1)] + r2_30 * ((*dist_coeffs_11)[int(2)] + r2_30 * (*dist_coeffs_11)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_11)[int(4)] * u_30 * v_30 + (*dist_coeffs_11)[int(5)] * (r2_30 + 2.0f * u_30 * u_30) + (*dist_coeffs_11)[int(6)] * r2_30, 2.0f * (*dist_coeffs_11)[int(5)] * u_30 * v_30 + (*dist_coeffs_11)[int(4)] * (r2_30 + 2.0f * v_30 * v_30) + (*dist_coeffs_11)[int(7)] * r2_30);
            float2  _S548 = _S547 + make_float2 ((*dist_coeffs_11)[int(8)] * _S547.x + (*dist_coeffs_11)[int(9)] * _S547.y, 0.0f);
            proj_points_2[int(6)] = make_float2 (intrins_4.x * _S548.x + intrins_4.z, intrins_4.y * _S548.y + intrins_4.w);
            break;
        }
        _S497 = _S494 & (!_S496);
        break;
    }
    if(!_S497)
    {
        return false;
    }
    float2  _S549 = *mean2d_5 + make_float2 (sigmas_2->w_mean_0[int(0)]) * *_S481;
    *mean2d_5 = _S549;
    float2  _S550 = _S549 + make_float2 (sigmas_2->w_mean_0[int(1)]) * *_S483;
    *mean2d_5 = _S550;
    float2  _S551 = _S550 + make_float2 (sigmas_2->w_mean_0[int(2)]) * *_S485;
    *mean2d_5 = _S551;
    float2  _S552 = _S551 + make_float2 (sigmas_2->w_mean_0[int(3)]) * *_S488;
    *mean2d_5 = _S552;
    float2  _S553 = _S552 + make_float2 (sigmas_2->w_mean_0[int(4)]) * *_S490;
    *mean2d_5 = _S553;
    float2  _S554 = _S553 + make_float2 (sigmas_2->w_mean_0[int(5)]) * *_S492;
    *mean2d_5 = _S554;
    float2  _S555 = _S554 + make_float2 (sigmas_2->w_mean_0[int(6)]) * *_S495;
    *mean2d_5 = _S555;
    float2  d_14 = *_S481 - _S555;
    float _S556 = d_14.x;
    float _S557 = d_14.y;
    float _S558 = _S556 * _S557;
    Matrix<float, 2, 2>  _S559 = *cov2d_5 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S556 * _S556, _S558, _S558, _S557 * _S557);
    *cov2d_5 = _S559;
    float2  d_15 = *_S483 - *mean2d_5;
    float _S560 = d_15.x;
    float _S561 = d_15.y;
    float _S562 = _S560 * _S561;
    Matrix<float, 2, 2>  _S563 = _S559 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S560 * _S560, _S562, _S562, _S561 * _S561);
    *cov2d_5 = _S563;
    float2  d_16 = *_S485 - *mean2d_5;
    float _S564 = d_16.x;
    float _S565 = d_16.y;
    float _S566 = _S564 * _S565;
    Matrix<float, 2, 2>  _S567 = _S563 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S564 * _S564, _S566, _S566, _S565 * _S565);
    *cov2d_5 = _S567;
    float2  d_17 = *_S488 - *mean2d_5;
    float _S568 = d_17.x;
    float _S569 = d_17.y;
    float _S570 = _S568 * _S569;
    Matrix<float, 2, 2>  _S571 = _S567 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S568 * _S568, _S570, _S570, _S569 * _S569);
    *cov2d_5 = _S571;
    float2  d_18 = *_S490 - *mean2d_5;
    float _S572 = d_18.x;
    float _S573 = d_18.y;
    float _S574 = _S572 * _S573;
    Matrix<float, 2, 2>  _S575 = _S571 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S572 * _S572, _S574, _S574, _S573 * _S573);
    *cov2d_5 = _S575;
    float2  d_19 = *_S492 - *mean2d_5;
    float _S576 = d_19.x;
    float _S577 = d_19.y;
    float _S578 = _S576 * _S577;
    Matrix<float, 2, 2>  _S579 = _S575 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S576 * _S576, _S578, _S578, _S577 * _S577);
    *cov2d_5 = _S579;
    float2  d_20 = *_S495 - *mean2d_5;
    float _S580 = d_20.x;
    float _S581 = d_20.y;
    float _S582 = _S580 * _S581;
    *cov2d_5 = _S579 + makeMatrix<float, 2, 2> (sigmas_2->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S580 * _S580, _S582, _S582, _S581 * _S581);
    return true;
}

inline __device__ bool fisheye_proj_3dgs_ut_1(SigmaPoints_0 * sigmas_3, float4  intrins_5, FixedArray<float, 10>  * dist_coeffs_12, Matrix<float, 2, 2>  * cov2d_6, float2  * mean2d_6)
{
    int2  _S583 = make_int2 (int(0));
    float2  _S584 = make_float2 ((float)_S583.x, (float)_S583.y);
    *mean2d_6 = _S584;
    *cov2d_6 = makeMatrix<float, 2, 2> (0.0f);
    FixedArray<float2 , 7>  proj_points_3;
    float2  _S585 = float2 {sigmas_3->p_0[int(0)].x, sigmas_3->p_0[int(0)].y};
    float r_17 = length_0(_S585);
    float _S586 = sigmas_3->p_0[int(0)].z;
    float theta_13 = (F32_atan2((r_17), (_S586)));
    float k_3;
    if(theta_13 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_13 * theta_13 / 3.0f) / _S586;
    }
    else
    {
        k_3 = theta_13 / r_17;
    }
    float2  _S587 = _S585 * make_float2 (k_3);
    float u_31 = _S587.x;
    float v_31 = _S587.y;
    float r2_31 = u_31 * u_31 + v_31 * v_31;
    float _S588 = 2.0f * (*dist_coeffs_12)[int(4)];
    float _S589 = 2.0f * (*dist_coeffs_12)[int(5)];
    float2  _S590 = _S587 * make_float2 (1.0f + r2_31 * ((*dist_coeffs_12)[int(0)] + r2_31 * ((*dist_coeffs_12)[int(1)] + r2_31 * ((*dist_coeffs_12)[int(2)] + r2_31 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S588 * u_31 * v_31 + (*dist_coeffs_12)[int(5)] * (r2_31 + 2.0f * u_31 * u_31) + (*dist_coeffs_12)[int(6)] * r2_31, _S589 * u_31 * v_31 + (*dist_coeffs_12)[int(4)] * (r2_31 + 2.0f * v_31 * v_31) + (*dist_coeffs_12)[int(7)] * r2_31);
    float2  _S591 = _S590 + make_float2 ((*dist_coeffs_12)[int(8)] * _S590.x + (*dist_coeffs_12)[int(9)] * _S590.y, 0.0f);
    float fx_5 = intrins_5.x;
    float fy_5 = intrins_5.y;
    float cx_2 = intrins_5.z;
    float cy_2 = intrins_5.w;
    proj_points_3[int(0)] = make_float2 (fx_5 * _S591.x + cx_2, fy_5 * _S591.y + cy_2);
    float2  _S592 = float2 {sigmas_3->p_0[int(1)].x, sigmas_3->p_0[int(1)].y};
    float r_18 = length_0(_S592);
    float _S593 = sigmas_3->p_0[int(1)].z;
    float theta_14 = (F32_atan2((r_18), (_S593)));
    if(theta_14 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_14 * theta_14 / 3.0f) / _S593;
    }
    else
    {
        k_3 = theta_14 / r_18;
    }
    float2  _S594 = _S592 * make_float2 (k_3);
    float u_32 = _S594.x;
    float v_32 = _S594.y;
    float r2_32 = u_32 * u_32 + v_32 * v_32;
    float2  _S595 = _S594 * make_float2 (1.0f + r2_32 * ((*dist_coeffs_12)[int(0)] + r2_32 * ((*dist_coeffs_12)[int(1)] + r2_32 * ((*dist_coeffs_12)[int(2)] + r2_32 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S588 * u_32 * v_32 + (*dist_coeffs_12)[int(5)] * (r2_32 + 2.0f * u_32 * u_32) + (*dist_coeffs_12)[int(6)] * r2_32, _S589 * u_32 * v_32 + (*dist_coeffs_12)[int(4)] * (r2_32 + 2.0f * v_32 * v_32) + (*dist_coeffs_12)[int(7)] * r2_32);
    float2  _S596 = _S595 + make_float2 ((*dist_coeffs_12)[int(8)] * _S595.x + (*dist_coeffs_12)[int(9)] * _S595.y, 0.0f);
    proj_points_3[int(1)] = make_float2 (fx_5 * _S596.x + cx_2, fy_5 * _S596.y + cy_2);
    float2  _S597 = float2 {sigmas_3->p_0[int(2)].x, sigmas_3->p_0[int(2)].y};
    float r_19 = length_0(_S597);
    float _S598 = sigmas_3->p_0[int(2)].z;
    float theta_15 = (F32_atan2((r_19), (_S598)));
    if(theta_15 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_15 * theta_15 / 3.0f) / _S598;
    }
    else
    {
        k_3 = theta_15 / r_19;
    }
    float2  _S599 = _S597 * make_float2 (k_3);
    float u_33 = _S599.x;
    float v_33 = _S599.y;
    float r2_33 = u_33 * u_33 + v_33 * v_33;
    float2  _S600 = _S599 * make_float2 (1.0f + r2_33 * ((*dist_coeffs_12)[int(0)] + r2_33 * ((*dist_coeffs_12)[int(1)] + r2_33 * ((*dist_coeffs_12)[int(2)] + r2_33 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S588 * u_33 * v_33 + (*dist_coeffs_12)[int(5)] * (r2_33 + 2.0f * u_33 * u_33) + (*dist_coeffs_12)[int(6)] * r2_33, _S589 * u_33 * v_33 + (*dist_coeffs_12)[int(4)] * (r2_33 + 2.0f * v_33 * v_33) + (*dist_coeffs_12)[int(7)] * r2_33);
    float2  _S601 = _S600 + make_float2 ((*dist_coeffs_12)[int(8)] * _S600.x + (*dist_coeffs_12)[int(9)] * _S600.y, 0.0f);
    proj_points_3[int(2)] = make_float2 (fx_5 * _S601.x + cx_2, fy_5 * _S601.y + cy_2);
    float2  _S602 = float2 {sigmas_3->p_0[int(3)].x, sigmas_3->p_0[int(3)].y};
    float r_20 = length_0(_S602);
    float _S603 = sigmas_3->p_0[int(3)].z;
    float theta_16 = (F32_atan2((r_20), (_S603)));
    if(theta_16 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_16 * theta_16 / 3.0f) / _S603;
    }
    else
    {
        k_3 = theta_16 / r_20;
    }
    float2  _S604 = _S602 * make_float2 (k_3);
    float u_34 = _S604.x;
    float v_34 = _S604.y;
    float r2_34 = u_34 * u_34 + v_34 * v_34;
    float2  _S605 = _S604 * make_float2 (1.0f + r2_34 * ((*dist_coeffs_12)[int(0)] + r2_34 * ((*dist_coeffs_12)[int(1)] + r2_34 * ((*dist_coeffs_12)[int(2)] + r2_34 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S588 * u_34 * v_34 + (*dist_coeffs_12)[int(5)] * (r2_34 + 2.0f * u_34 * u_34) + (*dist_coeffs_12)[int(6)] * r2_34, _S589 * u_34 * v_34 + (*dist_coeffs_12)[int(4)] * (r2_34 + 2.0f * v_34 * v_34) + (*dist_coeffs_12)[int(7)] * r2_34);
    float2  _S606 = _S605 + make_float2 ((*dist_coeffs_12)[int(8)] * _S605.x + (*dist_coeffs_12)[int(9)] * _S605.y, 0.0f);
    proj_points_3[int(3)] = make_float2 (fx_5 * _S606.x + cx_2, fy_5 * _S606.y + cy_2);
    float2  _S607 = float2 {sigmas_3->p_0[int(4)].x, sigmas_3->p_0[int(4)].y};
    float r_21 = length_0(_S607);
    float _S608 = sigmas_3->p_0[int(4)].z;
    float theta_17 = (F32_atan2((r_21), (_S608)));
    if(theta_17 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_17 * theta_17 / 3.0f) / _S608;
    }
    else
    {
        k_3 = theta_17 / r_21;
    }
    float2  _S609 = _S607 * make_float2 (k_3);
    float u_35 = _S609.x;
    float v_35 = _S609.y;
    float r2_35 = u_35 * u_35 + v_35 * v_35;
    float2  _S610 = _S609 * make_float2 (1.0f + r2_35 * ((*dist_coeffs_12)[int(0)] + r2_35 * ((*dist_coeffs_12)[int(1)] + r2_35 * ((*dist_coeffs_12)[int(2)] + r2_35 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S588 * u_35 * v_35 + (*dist_coeffs_12)[int(5)] * (r2_35 + 2.0f * u_35 * u_35) + (*dist_coeffs_12)[int(6)] * r2_35, _S589 * u_35 * v_35 + (*dist_coeffs_12)[int(4)] * (r2_35 + 2.0f * v_35 * v_35) + (*dist_coeffs_12)[int(7)] * r2_35);
    float2  _S611 = _S610 + make_float2 ((*dist_coeffs_12)[int(8)] * _S610.x + (*dist_coeffs_12)[int(9)] * _S610.y, 0.0f);
    proj_points_3[int(4)] = make_float2 (fx_5 * _S611.x + cx_2, fy_5 * _S611.y + cy_2);
    float2  _S612 = float2 {sigmas_3->p_0[int(5)].x, sigmas_3->p_0[int(5)].y};
    float r_22 = length_0(_S612);
    float _S613 = sigmas_3->p_0[int(5)].z;
    float theta_18 = (F32_atan2((r_22), (_S613)));
    if(theta_18 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_18 * theta_18 / 3.0f) / _S613;
    }
    else
    {
        k_3 = theta_18 / r_22;
    }
    float2  _S614 = _S612 * make_float2 (k_3);
    float u_36 = _S614.x;
    float v_36 = _S614.y;
    float r2_36 = u_36 * u_36 + v_36 * v_36;
    float2  _S615 = _S614 * make_float2 (1.0f + r2_36 * ((*dist_coeffs_12)[int(0)] + r2_36 * ((*dist_coeffs_12)[int(1)] + r2_36 * ((*dist_coeffs_12)[int(2)] + r2_36 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S588 * u_36 * v_36 + (*dist_coeffs_12)[int(5)] * (r2_36 + 2.0f * u_36 * u_36) + (*dist_coeffs_12)[int(6)] * r2_36, _S589 * u_36 * v_36 + (*dist_coeffs_12)[int(4)] * (r2_36 + 2.0f * v_36 * v_36) + (*dist_coeffs_12)[int(7)] * r2_36);
    float2  _S616 = _S615 + make_float2 ((*dist_coeffs_12)[int(8)] * _S615.x + (*dist_coeffs_12)[int(9)] * _S615.y, 0.0f);
    proj_points_3[int(5)] = make_float2 (fx_5 * _S616.x + cx_2, fy_5 * _S616.y + cy_2);
    float2  _S617 = float2 {sigmas_3->p_0[int(6)].x, sigmas_3->p_0[int(6)].y};
    float r_23 = length_0(_S617);
    float _S618 = sigmas_3->p_0[int(6)].z;
    float theta_19 = (F32_atan2((r_23), (_S618)));
    if(theta_19 < 0.00100000004749745f)
    {
        k_3 = (1.0f - theta_19 * theta_19 / 3.0f) / _S618;
    }
    else
    {
        k_3 = theta_19 / r_23;
    }
    float2  _S619 = _S617 * make_float2 (k_3);
    float u_37 = _S619.x;
    float v_37 = _S619.y;
    float r2_37 = u_37 * u_37 + v_37 * v_37;
    float2  _S620 = _S619 * make_float2 (1.0f + r2_37 * ((*dist_coeffs_12)[int(0)] + r2_37 * ((*dist_coeffs_12)[int(1)] + r2_37 * ((*dist_coeffs_12)[int(2)] + r2_37 * (*dist_coeffs_12)[int(3)])))) + make_float2 (_S588 * u_37 * v_37 + (*dist_coeffs_12)[int(5)] * (r2_37 + 2.0f * u_37 * u_37) + (*dist_coeffs_12)[int(6)] * r2_37, _S589 * u_37 * v_37 + (*dist_coeffs_12)[int(4)] * (r2_37 + 2.0f * v_37 * v_37) + (*dist_coeffs_12)[int(7)] * r2_37);
    float2  _S621 = _S620 + make_float2 ((*dist_coeffs_12)[int(8)] * _S620.x + (*dist_coeffs_12)[int(9)] * _S620.y, 0.0f);
    float2  _S622 = make_float2 (fx_5 * _S621.x + cx_2, fy_5 * _S621.y + cy_2);
    proj_points_3[int(6)] = _S622;
    float2  _S623 = *mean2d_6 + make_float2 (sigmas_3->w_mean_0[int(0)]) * proj_points_3[int(0)];
    *mean2d_6 = _S623;
    float2  _S624 = _S623 + make_float2 (sigmas_3->w_mean_0[int(1)]) * proj_points_3[int(1)];
    *mean2d_6 = _S624;
    float2  _S625 = _S624 + make_float2 (sigmas_3->w_mean_0[int(2)]) * proj_points_3[int(2)];
    *mean2d_6 = _S625;
    float2  _S626 = _S625 + make_float2 (sigmas_3->w_mean_0[int(3)]) * proj_points_3[int(3)];
    *mean2d_6 = _S626;
    float2  _S627 = _S626 + make_float2 (sigmas_3->w_mean_0[int(4)]) * proj_points_3[int(4)];
    *mean2d_6 = _S627;
    float2  _S628 = _S627 + make_float2 (sigmas_3->w_mean_0[int(5)]) * proj_points_3[int(5)];
    *mean2d_6 = _S628;
    float2  _S629 = _S628 + make_float2 (sigmas_3->w_mean_0[int(6)]) * _S622;
    *mean2d_6 = _S629;
    float2  d_21 = proj_points_3[int(0)] - _S629;
    float _S630 = d_21.x;
    float _S631 = d_21.y;
    float _S632 = _S630 * _S631;
    Matrix<float, 2, 2>  _S633 = *cov2d_6 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(0)]) * makeMatrix<float, 2, 2> (_S630 * _S630, _S632, _S632, _S631 * _S631);
    *cov2d_6 = _S633;
    float2  d_22 = proj_points_3[int(1)] - *mean2d_6;
    float _S634 = d_22.x;
    float _S635 = d_22.y;
    float _S636 = _S634 * _S635;
    Matrix<float, 2, 2>  _S637 = _S633 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(1)]) * makeMatrix<float, 2, 2> (_S634 * _S634, _S636, _S636, _S635 * _S635);
    *cov2d_6 = _S637;
    float2  d_23 = proj_points_3[int(2)] - *mean2d_6;
    float _S638 = d_23.x;
    float _S639 = d_23.y;
    float _S640 = _S638 * _S639;
    Matrix<float, 2, 2>  _S641 = _S637 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(2)]) * makeMatrix<float, 2, 2> (_S638 * _S638, _S640, _S640, _S639 * _S639);
    *cov2d_6 = _S641;
    float2  d_24 = proj_points_3[int(3)] - *mean2d_6;
    float _S642 = d_24.x;
    float _S643 = d_24.y;
    float _S644 = _S642 * _S643;
    Matrix<float, 2, 2>  _S645 = _S641 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(3)]) * makeMatrix<float, 2, 2> (_S642 * _S642, _S644, _S644, _S643 * _S643);
    *cov2d_6 = _S645;
    float2  d_25 = proj_points_3[int(4)] - *mean2d_6;
    float _S646 = d_25.x;
    float _S647 = d_25.y;
    float _S648 = _S646 * _S647;
    Matrix<float, 2, 2>  _S649 = _S645 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(4)]) * makeMatrix<float, 2, 2> (_S646 * _S646, _S648, _S648, _S647 * _S647);
    *cov2d_6 = _S649;
    float2  d_26 = proj_points_3[int(5)] - *mean2d_6;
    float _S650 = d_26.x;
    float _S651 = d_26.y;
    float _S652 = _S650 * _S651;
    Matrix<float, 2, 2>  _S653 = _S649 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(5)]) * makeMatrix<float, 2, 2> (_S650 * _S650, _S652, _S652, _S651 * _S651);
    *cov2d_6 = _S653;
    float2  d_27 = _S622 - *mean2d_6;
    float _S654 = d_27.x;
    float _S655 = d_27.y;
    float _S656 = _S654 * _S655;
    *cov2d_6 = _S653 + makeMatrix<float, 2, 2> (sigmas_3->w_cov_0[int(6)]) * makeMatrix<float, 2, 2> (_S654 * _S654, _S656, _S656, _S655 * _S655);
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
    float _S657 = *&((covar_1->rows + (int(0)))->x) + eps2d_0;
    *&((covar_1->rows + (int(0)))->x) = _S657;
    float _S658 = *&((covar_1->rows + (int(1)))->y) + eps2d_0;
    *&((covar_1->rows + (int(1)))->y) = _S658;
    float det_blur_0 = _S657 * _S658 - *&((covar_1->rows + (int(0)))->y) * *&((covar_1->rows + (int(1)))->x);
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

inline __device__ void _d_exp_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_9, float3  dOut_13)
{
    float3  _S659 = exp_0((*dpx_9).primal_0) * dOut_13;
    dpx_9->primal_0 = (*dpx_9).primal_0;
    dpx_9->differential_0 = _S659;
    return;
}

inline __device__ void _d_log_0(DiffPair_float_0 * dpx_10, float dOut_14)
{
    float _S660 = 1.0f / (*dpx_10).primal_0 * dOut_14;
    dpx_10->primal_0 = (*dpx_10).primal_0;
    dpx_10->differential_0 = _S660;
    return;
}

inline __device__ void _d_max_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_11, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpy_5, float3  dOut_15)
{
    DiffPair_float_0 left_dp_0;
    (&left_dp_0)->primal_0 = (*dpx_11).primal_0.x;
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
    (&left_dp_1)->primal_0 = (*dpx_11).primal_0.y;
    (&left_dp_1)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_1;
    (&right_dp_1)->primal_0 = (*dpy_5).primal_0.y;
    (&right_dp_1)->differential_0 = 0.0f;
    _d_max_0(&left_dp_1, &right_dp_1, dOut_15.y);
    *&((&left_d_result_5)->y) = left_dp_1.differential_0;
    *&((&right_d_result_5)->y) = right_dp_1.differential_0;
    DiffPair_float_0 left_dp_2;
    (&left_dp_2)->primal_0 = (*dpx_11).primal_0.z;
    (&left_dp_2)->differential_0 = 0.0f;
    DiffPair_float_0 right_dp_2;
    (&right_dp_2)->primal_0 = (*dpy_5).primal_0.z;
    (&right_dp_2)->differential_0 = 0.0f;
    _d_max_0(&left_dp_2, &right_dp_2, dOut_15.z);
    *&((&left_d_result_5)->z) = left_dp_2.differential_0;
    *&((&right_d_result_5)->z) = right_dp_2.differential_0;
    dpx_11->primal_0 = (*dpx_11).primal_0;
    dpx_11->differential_0 = left_d_result_5;
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

inline __device__ void projection_3dgs_persp(bool antialiased_0, float3  mean_1, float4  quat_4, float3  scale_3, float in_opacity_0, FixedArray<float3 , 16>  * sh_coeffs_0, Matrix<float, 3, 3>  R_6, float3  t_3, float fx_7, float fy_7, float cx_4, float cy_4, FixedArray<float, 10>  * dist_coeffs_13, uint image_width_0, uint image_height_0, float near_plane_0, float far_plane_0, int4  * aabb_xyxy_0, float2  * mean2d_8, float * depth_0, float3  * conic_0, float * opacity_0, float3  * rgb_0)
{
    for(;;)
    {
        float3  mean_c_0 = mul_0(R_6, mean_1) + t_3;
        float _S661 = mean_c_0.z;
        bool _S662;
        if(_S661 < near_plane_0)
        {
            _S662 = true;
        }
        else
        {
            _S662 = _S661 > far_plane_0;
        }
        if(_S662)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_0;
        float3  _S663 = exp_0(scale_3);
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
        Matrix<float, 3, 3>  M_2 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_4 + z2_4), 2.0f * (xy_4 + wz_4), 2.0f * (xz_4 - wy_4), 2.0f * (xy_4 - wz_4), 1.0f - 2.0f * (x2_4 + z2_4), 2.0f * (yz_4 + wx_4), 2.0f * (xz_4 + wy_4), 2.0f * (yz_4 - wx_4), 1.0f - 2.0f * (x2_4 + y2_4))), makeMatrix<float, 3, 3> (_S663.x, 0.0f, 0.0f, 0.0f, _S663.y, 0.0f, 0.0f, 0.0f, _S663.z));
        Matrix<float, 3, 3>  _S664 = transpose_0(R_6);
        float _S665 = float(image_width_0);
        float _S666 = float(image_height_0);
        float _S667 = 0.30000001192092896f * (0.5f * _S665 / fx_7);
        float _S668 = 0.30000001192092896f * (0.5f * _S666 / fy_7);
        float rz_1 = 1.0f / mean_c_0.z;
        float rz2_1 = rz_1 * rz_1;
        Matrix<float, 2, 3>  J_4 = makeMatrix<float, 2, 3> (fx_7 * rz_1, 0.0f, - fx_7 * (mean_c_0.z * (F32_min(((_S665 - cx_4) / fx_7 + _S667), ((F32_max((- (cx_4 / fx_7 + _S667)), (mean_c_0.x * rz_1))))))) * rz2_1, 0.0f, fy_7 * rz_1, - fy_7 * (mean_c_0.z * (F32_min(((_S666 - cy_4) / fy_7 + _S668), ((F32_max((- (cy_4 / fy_7 + _S668)), (mean_c_0.y * rz_1))))))) * rz2_1);
        covar2d_0 = mul_6(mul_5(J_4, mul_4(mul_4(R_6, mul_4(M_2, transpose_0(M_2))), _S664)), transpose_1(J_4));
        *mean2d_8 = make_float2 (fx_7 * mean_c_0.x * rz_1 + cx_4, fy_7 * mean_c_0.y * rz_1 + cy_4);
        float det_orig_1 = *&(((&covar2d_0)->rows + (int(0)))->x) * *&(((&covar2d_0)->rows + (int(1)))->y) - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float _S669 = *&(((&covar2d_0)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(0)))->x) = _S669;
        float _S670 = *&(((&covar2d_0)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_0)->rows + (int(1)))->y) = _S670;
        float det_blur_1 = _S669 * _S670 - *&(((&covar2d_0)->rows + (int(0)))->y) * *&(((&covar2d_0)->rows + (int(1)))->x);
        float compensation_1 = (F32_sqrt(((F32_max((0.0f), (det_orig_1 / det_blur_1))))));
        if(det_blur_1 <= 0.0f)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_1 = 1.0f / (covar2d_0.rows[int(0)].x * covar2d_0.rows[int(1)].y - covar2d_0.rows[int(0)].y * covar2d_0.rows[int(1)].x);
        Matrix<float, 2, 2>  _S671 = makeMatrix<float, 2, 2> (covar2d_0.rows[int(1)].y * invdet_1, - covar2d_0.rows[int(0)].y * invdet_1, - covar2d_0.rows[int(1)].x * invdet_1, covar2d_0.rows[int(0)].x * invdet_1);
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
            _S662 = true;
        }
        else
        {
            _S662 = xmin_0 >= _S665;
        }
        if(_S662)
        {
            _S662 = true;
        }
        else
        {
            _S662 = ymax_0 <= 0.0f;
        }
        if(_S662)
        {
            _S662 = true;
        }
        else
        {
            _S662 = ymin_0 >= _S666;
        }
        if(_S662)
        {
            *aabb_xyxy_0 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_0 = make_int4 (int(xmin_0), int(ymin_0), int(xmax_0), int(ymax_0));
        *depth_0 = 0.5f * (F32_log((dot_0(mean_c_0, mean_c_0) + 9.99999997475242708e-07f)));
        *conic_0 = make_float3 (_S671.rows[int(0)].x, _S671.rows[int(0)].y, _S671.rows[int(1)].y);
        float3  _S672 = mean_1 - - mul_0(_S664, t_3);
        float3  _S673 = make_float3 (0.282094806432724f) * (*sh_coeffs_0)[int(0)];
        *rgb_0 = _S673;
        float _S674 = _S672.x;
        float _S675 = _S672.y;
        float _S676 = _S672.z;
        float norm_0 = (F32_sqrt((_S674 * _S674 + _S675 * _S675 + _S676 * _S676)));
        float x_17 = _S674 / norm_0;
        float y_3 = _S675 / norm_0;
        float z_0 = _S676 / norm_0;
        float3  _S677 = _S673 + make_float3 (0.48860251903533936f) * (make_float3 (- y_3) * (*sh_coeffs_0)[int(1)] + make_float3 (z_0) * (*sh_coeffs_0)[int(2)] - make_float3 (x_17) * (*sh_coeffs_0)[int(3)]);
        *rgb_0 = _S677;
        float z2_5 = z_0 * z_0;
        float fTmp0B_0 = -1.09254848957061768f * z_0;
        float fC1_0 = x_17 * x_17 - y_3 * y_3;
        float fS1_0 = 2.0f * x_17 * y_3;
        float3  _S678 = _S677 + (make_float3 (0.54627424478530884f * fS1_0) * (*sh_coeffs_0)[int(4)] + make_float3 (fTmp0B_0 * y_3) * (*sh_coeffs_0)[int(5)] + make_float3 (0.94617468118667603f * z2_5 - 0.31539157032966614f) * (*sh_coeffs_0)[int(6)] + make_float3 (fTmp0B_0 * x_17) * (*sh_coeffs_0)[int(7)] + make_float3 (0.54627424478530884f * fC1_0) * (*sh_coeffs_0)[int(8)]);
        *rgb_0 = _S678;
        float fTmp0C_0 = -2.28522896766662598f * z2_5 + 0.4570457935333252f;
        float fTmp1B_0 = 1.44530570507049561f * z_0;
        *rgb_0 = max_0(_S678 + (make_float3 (-0.59004360437393188f * (x_17 * fS1_0 + y_3 * fC1_0)) * (*sh_coeffs_0)[int(9)] + make_float3 (fTmp1B_0 * fS1_0) * (*sh_coeffs_0)[int(10)] + make_float3 (fTmp0C_0 * y_3) * (*sh_coeffs_0)[int(11)] + make_float3 (z_0 * (1.86588168144226074f * z2_5 - 1.11952900886535645f)) * (*sh_coeffs_0)[int(12)] + make_float3 (fTmp0C_0 * x_17) * (*sh_coeffs_0)[int(13)] + make_float3 (fTmp1B_0 * fC1_0) * (*sh_coeffs_0)[int(14)] + make_float3 (-0.59004360437393188f * (x_17 * fC1_0 - y_3 * fS1_0)) * (*sh_coeffs_0)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_fisheye(bool antialiased_1, float3  mean_2, float4  quat_5, float3  scale_4, float in_opacity_1, FixedArray<float3 , 16>  * sh_coeffs_1, Matrix<float, 3, 3>  R_7, float3  t_4, float fx_8, float fy_8, float cx_5, float cy_5, FixedArray<float, 10>  * dist_coeffs_14, uint image_width_1, uint image_height_1, float near_plane_1, float far_plane_1, int4  * aabb_xyxy_1, float2  * mean2d_9, float * depth_1, float3  * conic_1, float * opacity_1, float3  * rgb_1)
{
    for(;;)
    {
        float3  mean_c_1 = mul_0(R_7, mean_2) + t_4;
        float _S679 = length_1(mean_c_1);
        bool _S680;
        if(_S679 < near_plane_1)
        {
            _S680 = true;
        }
        else
        {
            _S680 = _S679 > far_plane_1;
        }
        if(_S680)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_6 = make_float4 (fx_8, fy_8, cx_5, cy_5);
        Matrix<float, 2, 2>  covar2d_1;
        float3  _S681 = exp_0(scale_4);
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
        Matrix<float, 3, 3>  M_3 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_5 + z2_6), 2.0f * (xy_5 + wz_5), 2.0f * (xz_5 - wy_5), 2.0f * (xy_5 - wz_5), 1.0f - 2.0f * (x2_5 + z2_6), 2.0f * (yz_5 + wx_5), 2.0f * (xz_5 + wy_5), 2.0f * (yz_5 - wx_5), 1.0f - 2.0f * (x2_5 + y2_5))), makeMatrix<float, 3, 3> (_S681.x, 0.0f, 0.0f, 0.0f, _S681.y, 0.0f, 0.0f, 0.0f, _S681.z));
        Matrix<float, 3, 3>  _S682 = transpose_0(R_7);
        bool _S683 = fisheye_proj_3dgs_0(mean_c_1, mul_4(mul_4(R_7, mul_4(M_3, transpose_0(M_3))), _S682), intrins_6, dist_coeffs_14, &covar2d_1, mean2d_9);
        if(!(true & _S683))
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_2 = *&(((&covar2d_1)->rows + (int(0)))->x) * *&(((&covar2d_1)->rows + (int(1)))->y) - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float _S684 = *&(((&covar2d_1)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(0)))->x) = _S684;
        float _S685 = *&(((&covar2d_1)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_1)->rows + (int(1)))->y) = _S685;
        float det_blur_2 = _S684 * _S685 - *&(((&covar2d_1)->rows + (int(0)))->y) * *&(((&covar2d_1)->rows + (int(1)))->x);
        float compensation_2 = (F32_sqrt(((F32_max((0.0f), (det_orig_2 / det_blur_2))))));
        if(det_blur_2 <= 0.0f)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_2 = 1.0f / (covar2d_1.rows[int(0)].x * covar2d_1.rows[int(1)].y - covar2d_1.rows[int(0)].y * covar2d_1.rows[int(1)].x);
        Matrix<float, 2, 2>  _S686 = makeMatrix<float, 2, 2> (covar2d_1.rows[int(1)].y * invdet_2, - covar2d_1.rows[int(0)].y * invdet_2, - covar2d_1.rows[int(1)].x * invdet_2, covar2d_1.rows[int(0)].x * invdet_2);
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
            _S680 = true;
        }
        else
        {
            _S680 = xmin_1 >= float(image_width_1);
        }
        if(_S680)
        {
            _S680 = true;
        }
        else
        {
            _S680 = ymax_1 <= 0.0f;
        }
        if(_S680)
        {
            _S680 = true;
        }
        else
        {
            _S680 = ymin_1 >= float(image_height_1);
        }
        if(_S680)
        {
            *aabb_xyxy_1 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_1 = make_int4 (int(xmin_1), int(ymin_1), int(xmax_1), int(ymax_1));
        *depth_1 = 0.5f * (F32_log((dot_0(mean_c_1, mean_c_1) + 9.99999997475242708e-07f)));
        *conic_1 = make_float3 (_S686.rows[int(0)].x, _S686.rows[int(0)].y, _S686.rows[int(1)].y);
        float3  _S687 = mean_2 - - mul_0(_S682, t_4);
        float3  _S688 = make_float3 (0.282094806432724f) * (*sh_coeffs_1)[int(0)];
        *rgb_1 = _S688;
        float _S689 = _S687.x;
        float _S690 = _S687.y;
        float _S691 = _S687.z;
        float norm_1 = (F32_sqrt((_S689 * _S689 + _S690 * _S690 + _S691 * _S691)));
        float x_19 = _S689 / norm_1;
        float y_4 = _S690 / norm_1;
        float z_1 = _S691 / norm_1;
        float3  _S692 = _S688 + make_float3 (0.48860251903533936f) * (make_float3 (- y_4) * (*sh_coeffs_1)[int(1)] + make_float3 (z_1) * (*sh_coeffs_1)[int(2)] - make_float3 (x_19) * (*sh_coeffs_1)[int(3)]);
        *rgb_1 = _S692;
        float z2_7 = z_1 * z_1;
        float fTmp0B_1 = -1.09254848957061768f * z_1;
        float fC1_1 = x_19 * x_19 - y_4 * y_4;
        float fS1_1 = 2.0f * x_19 * y_4;
        float3  _S693 = _S692 + (make_float3 (0.54627424478530884f * fS1_1) * (*sh_coeffs_1)[int(4)] + make_float3 (fTmp0B_1 * y_4) * (*sh_coeffs_1)[int(5)] + make_float3 (0.94617468118667603f * z2_7 - 0.31539157032966614f) * (*sh_coeffs_1)[int(6)] + make_float3 (fTmp0B_1 * x_19) * (*sh_coeffs_1)[int(7)] + make_float3 (0.54627424478530884f * fC1_1) * (*sh_coeffs_1)[int(8)]);
        *rgb_1 = _S693;
        float fTmp0C_1 = -2.28522896766662598f * z2_7 + 0.4570457935333252f;
        float fTmp1B_1 = 1.44530570507049561f * z_1;
        *rgb_1 = max_0(_S693 + (make_float3 (-0.59004360437393188f * (x_19 * fS1_1 + y_4 * fC1_1)) * (*sh_coeffs_1)[int(9)] + make_float3 (fTmp1B_1 * fS1_1) * (*sh_coeffs_1)[int(10)] + make_float3 (fTmp0C_1 * y_4) * (*sh_coeffs_1)[int(11)] + make_float3 (z_1 * (1.86588168144226074f * z2_7 - 1.11952900886535645f)) * (*sh_coeffs_1)[int(12)] + make_float3 (fTmp0C_1 * x_19) * (*sh_coeffs_1)[int(13)] + make_float3 (fTmp1B_1 * fC1_1) * (*sh_coeffs_1)[int(14)] + make_float3 (-0.59004360437393188f * (x_19 * fC1_1 - y_4 * fS1_1)) * (*sh_coeffs_1)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_ortho(bool antialiased_2, float3  mean_3, float4  quat_6, float3  scale_5, float in_opacity_2, FixedArray<float3 , 16>  * sh_coeffs_2, Matrix<float, 3, 3>  R_8, float3  t_5, float fx_9, float fy_9, float cx_6, float cy_6, FixedArray<float, 10>  * dist_coeffs_15, uint image_width_2, uint image_height_2, float near_plane_2, float far_plane_2, int4  * aabb_xyxy_2, float2  * mean2d_10, float * depth_2, float3  * conic_2, float * opacity_2, float3  * rgb_2)
{
    for(;;)
    {
        float3  mean_c_2 = mul_0(R_8, mean_3) + t_5;
        float _S694 = mean_c_2.z;
        bool _S695;
        if(_S694 < near_plane_2)
        {
            _S695 = true;
        }
        else
        {
            _S695 = _S694 > far_plane_2;
        }
        if(_S695)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        Matrix<float, 2, 2>  covar2d_2;
        float3  _S696 = exp_0(scale_5);
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
        Matrix<float, 3, 3>  M_4 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_6 + z2_8), 2.0f * (xy_6 + wz_6), 2.0f * (xz_6 - wy_6), 2.0f * (xy_6 - wz_6), 1.0f - 2.0f * (x2_6 + z2_8), 2.0f * (yz_6 + wx_6), 2.0f * (xz_6 + wy_6), 2.0f * (yz_6 - wx_6), 1.0f - 2.0f * (x2_6 + y2_6))), makeMatrix<float, 3, 3> (_S696.x, 0.0f, 0.0f, 0.0f, _S696.y, 0.0f, 0.0f, 0.0f, _S696.z));
        Matrix<float, 3, 3>  _S697 = transpose_0(R_8);
        Matrix<float, 2, 3>  J_5 = makeMatrix<float, 2, 3> (fx_9, 0.0f, 0.0f, 0.0f, fy_9, 0.0f);
        covar2d_2 = mul_6(mul_5(J_5, mul_4(mul_4(R_8, mul_4(M_4, transpose_0(M_4))), _S697)), transpose_1(J_5));
        *mean2d_10 = make_float2 (fx_9 * mean_c_2.x + cx_6, fy_9 * mean_c_2.y + cy_6);
        float det_orig_3 = *&(((&covar2d_2)->rows + (int(0)))->x) * *&(((&covar2d_2)->rows + (int(1)))->y) - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float _S698 = *&(((&covar2d_2)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(0)))->x) = _S698;
        float _S699 = *&(((&covar2d_2)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_2)->rows + (int(1)))->y) = _S699;
        float det_blur_3 = _S698 * _S699 - *&(((&covar2d_2)->rows + (int(0)))->y) * *&(((&covar2d_2)->rows + (int(1)))->x);
        float compensation_3 = (F32_sqrt(((F32_max((0.0f), (det_orig_3 / det_blur_3))))));
        if(det_blur_3 <= 0.0f)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float invdet_3 = 1.0f / (covar2d_2.rows[int(0)].x * covar2d_2.rows[int(1)].y - covar2d_2.rows[int(0)].y * covar2d_2.rows[int(1)].x);
        Matrix<float, 2, 2>  _S700 = makeMatrix<float, 2, 2> (covar2d_2.rows[int(1)].y * invdet_3, - covar2d_2.rows[int(0)].y * invdet_3, - covar2d_2.rows[int(1)].x * invdet_3, covar2d_2.rows[int(0)].x * invdet_3);
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
            _S695 = true;
        }
        else
        {
            _S695 = xmin_2 >= float(image_width_2);
        }
        if(_S695)
        {
            _S695 = true;
        }
        else
        {
            _S695 = ymax_2 <= 0.0f;
        }
        if(_S695)
        {
            _S695 = true;
        }
        else
        {
            _S695 = ymin_2 >= float(image_height_2);
        }
        if(_S695)
        {
            *aabb_xyxy_2 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_2 = make_int4 (int(xmin_2), int(ymin_2), int(xmax_2), int(ymax_2));
        *depth_2 = 0.5f * (F32_log((dot_0(mean_c_2, mean_c_2) + 9.99999997475242708e-07f)));
        *conic_2 = make_float3 (_S700.rows[int(0)].x, _S700.rows[int(0)].y, _S700.rows[int(1)].y);
        float3  _S701 = mean_3 - - mul_0(_S697, t_5);
        float3  _S702 = make_float3 (0.282094806432724f) * (*sh_coeffs_2)[int(0)];
        *rgb_2 = _S702;
        float _S703 = _S701.x;
        float _S704 = _S701.y;
        float _S705 = _S701.z;
        float norm_2 = (F32_sqrt((_S703 * _S703 + _S704 * _S704 + _S705 * _S705)));
        float x_21 = _S703 / norm_2;
        float y_5 = _S704 / norm_2;
        float z_2 = _S705 / norm_2;
        float3  _S706 = _S702 + make_float3 (0.48860251903533936f) * (make_float3 (- y_5) * (*sh_coeffs_2)[int(1)] + make_float3 (z_2) * (*sh_coeffs_2)[int(2)] - make_float3 (x_21) * (*sh_coeffs_2)[int(3)]);
        *rgb_2 = _S706;
        float z2_9 = z_2 * z_2;
        float fTmp0B_2 = -1.09254848957061768f * z_2;
        float fC1_2 = x_21 * x_21 - y_5 * y_5;
        float fS1_2 = 2.0f * x_21 * y_5;
        float3  _S707 = _S706 + (make_float3 (0.54627424478530884f * fS1_2) * (*sh_coeffs_2)[int(4)] + make_float3 (fTmp0B_2 * y_5) * (*sh_coeffs_2)[int(5)] + make_float3 (0.94617468118667603f * z2_9 - 0.31539157032966614f) * (*sh_coeffs_2)[int(6)] + make_float3 (fTmp0B_2 * x_21) * (*sh_coeffs_2)[int(7)] + make_float3 (0.54627424478530884f * fC1_2) * (*sh_coeffs_2)[int(8)]);
        *rgb_2 = _S707;
        float fTmp0C_2 = -2.28522896766662598f * z2_9 + 0.4570457935333252f;
        float fTmp1B_2 = 1.44530570507049561f * z_2;
        *rgb_2 = max_0(_S707 + (make_float3 (-0.59004360437393188f * (x_21 * fS1_2 + y_5 * fC1_2)) * (*sh_coeffs_2)[int(9)] + make_float3 (fTmp1B_2 * fS1_2) * (*sh_coeffs_2)[int(10)] + make_float3 (fTmp0C_2 * y_5) * (*sh_coeffs_2)[int(11)] + make_float3 (z_2 * (1.86588168144226074f * z2_9 - 1.11952900886535645f)) * (*sh_coeffs_2)[int(12)] + make_float3 (fTmp0C_2 * x_21) * (*sh_coeffs_2)[int(13)] + make_float3 (fTmp1B_2 * fC1_2) * (*sh_coeffs_2)[int(14)] + make_float3 (-0.59004360437393188f * (x_21 * fC1_2 - y_5 * fS1_2)) * (*sh_coeffs_2)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_persp(bool antialiased_3, float3  mean_4, float4  quat_7, float3  scale_6, float in_opacity_3, FixedArray<float3 , 16>  * sh_coeffs_3, Matrix<float, 3, 3>  R_9, float3  t_6, float fx_10, float fy_10, float cx_7, float cy_7, FixedArray<float, 10>  * dist_coeffs_16, uint image_width_3, uint image_height_3, float near_plane_3, float far_plane_3, int4  * aabb_xyxy_3, float2  * mean2d_11, float * depth_3, float3  * conic_3, float * opacity_3, float3  * rgb_3)
{
    for(;;)
    {
        float3  mean_c_3 = mul_0(R_9, mean_4) + t_6;
        float _S708 = mean_c_3.z;
        bool _S709;
        if(_S708 < near_plane_3)
        {
            _S709 = true;
        }
        else
        {
            _S709 = _S708 > far_plane_3;
        }
        if(_S709)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_7 = make_float4 (fx_10, fy_10, cx_7, cy_7);
        Matrix<float, 2, 2>  covar2d_3;
        float3  _S710 = exp_0(scale_6);
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
        Matrix<float, 3, 3>  _S711 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_7 + z2_10), 2.0f * (xy_7 + wz_7), 2.0f * (xz_7 - wy_7), 2.0f * (xy_7 - wz_7), 1.0f - 2.0f * (x2_7 + z2_10), 2.0f * (yz_7 + wx_7), 2.0f * (xz_7 + wy_7), 2.0f * (yz_7 - wx_7), 1.0f - 2.0f * (x2_7 + y2_7))));
        SigmaPoints_0 ret_0;
        (&ret_0)->p_0[int(0)] = mean_4;
        (&ret_0)->w_mean_0[int(0)] = 0.0f;
        (&ret_0)->w_cov_0[int(0)] = 2.0f;
        float _S712 = (F32_sqrt((3.0f)));
        float3  delta_0 = make_float3 (_S712 * _S710.x) * _S711.rows[0U];
        float3  _S713 = mean_4 + delta_0;
        float3  _S714 = mean_4 - delta_0;
        float3  delta_1 = make_float3 (_S712 * _S710.y) * _S711.rows[1U];
        float3  _S715 = mean_4 + delta_1;
        float3  _S716 = mean_4 - delta_1;
        float3  delta_2 = make_float3 (_S712 * _S710.z) * _S711.rows[2U];
        float3  _S717 = mean_4 + delta_2;
        float3  _S718 = mean_4 - delta_2;
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
        (&ret_0)->p_0[0U] = mul_0(R_9, (&ret_0)->p_0[0U]) + t_6;
        (&ret_0)->p_0[1U] = mul_0(R_9, _S713) + t_6;
        (&ret_0)->p_0[2U] = mul_0(R_9, _S715) + t_6;
        (&ret_0)->p_0[3U] = mul_0(R_9, _S717) + t_6;
        (&ret_0)->p_0[4U] = mul_0(R_9, _S714) + t_6;
        (&ret_0)->p_0[5U] = mul_0(R_9, _S716) + t_6;
        (&ret_0)->p_0[6U] = mul_0(R_9, _S718) + t_6;
        SigmaPoints_0 _S719 = ret_0;
        bool _S720 = persp_proj_3dgs_ut_0(&_S719, intrins_7, dist_coeffs_16, image_width_3, image_height_3, &covar2d_3, mean2d_11);
        if(!(true & _S720))
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_4 = *&(((&covar2d_3)->rows + (int(0)))->x) * *&(((&covar2d_3)->rows + (int(1)))->y) - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
        float _S721 = *&(((&covar2d_3)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(0)))->x) = _S721;
        float _S722 = *&(((&covar2d_3)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_3)->rows + (int(1)))->y) = _S722;
        float det_blur_4 = _S721 * _S722 - *&(((&covar2d_3)->rows + (int(0)))->y) * *&(((&covar2d_3)->rows + (int(1)))->x);
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
            _S709 = true;
        }
        else
        {
            _S709 = xmin_3 >= float(image_width_3);
        }
        if(_S709)
        {
            _S709 = true;
        }
        else
        {
            _S709 = ymax_3 <= 0.0f;
        }
        if(_S709)
        {
            _S709 = true;
        }
        else
        {
            _S709 = ymin_3 >= float(image_height_3);
        }
        if(_S709)
        {
            *aabb_xyxy_3 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_3 = make_int4 (int(xmin_3), int(ymin_3), int(xmax_3), int(ymax_3));
        *depth_3 = 0.5f * (F32_log((dot_0(mean_c_3, mean_c_3) + 9.99999997475242708e-07f)));
        *conic_3 = exp_0(- scale_6);
        float3  _S723 = mean_4 - - mul_0(transpose_0(R_9), t_6);
        float3  _S724 = make_float3 (0.282094806432724f) * (*sh_coeffs_3)[int(0)];
        *rgb_3 = _S724;
        float _S725 = _S723.x;
        float _S726 = _S723.y;
        float _S727 = _S723.z;
        float norm_3 = (F32_sqrt((_S725 * _S725 + _S726 * _S726 + _S727 * _S727)));
        float x_23 = _S725 / norm_3;
        float y_6 = _S726 / norm_3;
        float z_3 = _S727 / norm_3;
        float3  _S728 = _S724 + make_float3 (0.48860251903533936f) * (make_float3 (- y_6) * (*sh_coeffs_3)[int(1)] + make_float3 (z_3) * (*sh_coeffs_3)[int(2)] - make_float3 (x_23) * (*sh_coeffs_3)[int(3)]);
        *rgb_3 = _S728;
        float z2_11 = z_3 * z_3;
        float fTmp0B_3 = -1.09254848957061768f * z_3;
        float fC1_3 = x_23 * x_23 - y_6 * y_6;
        float fS1_3 = 2.0f * x_23 * y_6;
        float3  _S729 = _S728 + (make_float3 (0.54627424478530884f * fS1_3) * (*sh_coeffs_3)[int(4)] + make_float3 (fTmp0B_3 * y_6) * (*sh_coeffs_3)[int(5)] + make_float3 (0.94617468118667603f * z2_11 - 0.31539157032966614f) * (*sh_coeffs_3)[int(6)] + make_float3 (fTmp0B_3 * x_23) * (*sh_coeffs_3)[int(7)] + make_float3 (0.54627424478530884f * fC1_3) * (*sh_coeffs_3)[int(8)]);
        *rgb_3 = _S729;
        float fTmp0C_3 = -2.28522896766662598f * z2_11 + 0.4570457935333252f;
        float fTmp1B_3 = 1.44530570507049561f * z_3;
        *rgb_3 = max_0(_S729 + (make_float3 (-0.59004360437393188f * (x_23 * fS1_3 + y_6 * fC1_3)) * (*sh_coeffs_3)[int(9)] + make_float3 (fTmp1B_3 * fS1_3) * (*sh_coeffs_3)[int(10)] + make_float3 (fTmp0C_3 * y_6) * (*sh_coeffs_3)[int(11)] + make_float3 (z_3 * (1.86588168144226074f * z2_11 - 1.11952900886535645f)) * (*sh_coeffs_3)[int(12)] + make_float3 (fTmp0C_3 * x_23) * (*sh_coeffs_3)[int(13)] + make_float3 (fTmp1B_3 * fC1_3) * (*sh_coeffs_3)[int(14)] + make_float3 (-0.59004360437393188f * (x_23 * fC1_3 - y_6 * fS1_3)) * (*sh_coeffs_3)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_3dgs_eval3d_fisheye(bool antialiased_4, float3  mean_5, float4  quat_8, float3  scale_7, float in_opacity_4, FixedArray<float3 , 16>  * sh_coeffs_4, Matrix<float, 3, 3>  R_10, float3  t_7, float fx_11, float fy_11, float cx_8, float cy_8, FixedArray<float, 10>  * dist_coeffs_17, uint image_width_4, uint image_height_4, float near_plane_4, float far_plane_4, int4  * aabb_xyxy_4, float2  * mean2d_12, float * depth_4, float3  * conic_4, float * opacity_4, float3  * rgb_4)
{
    for(;;)
    {
        float3  mean_c_4 = mul_0(R_10, mean_5) + t_7;
        float _S730 = length_1(mean_c_4);
        bool _S731;
        if(_S730 < near_plane_4)
        {
            _S731 = true;
        }
        else
        {
            _S731 = _S730 > far_plane_4;
        }
        if(_S731)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float4  intrins_8 = make_float4 (fx_11, fy_11, cx_8, cy_8);
        Matrix<float, 2, 2>  covar2d_4;
        float3  _S732 = exp_0(scale_7);
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
        Matrix<float, 3, 3>  _S733 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_8 + z2_12), 2.0f * (xy_8 + wz_8), 2.0f * (xz_8 - wy_8), 2.0f * (xy_8 - wz_8), 1.0f - 2.0f * (x2_8 + z2_12), 2.0f * (yz_8 + wx_8), 2.0f * (xz_8 + wy_8), 2.0f * (yz_8 - wx_8), 1.0f - 2.0f * (x2_8 + y2_8))));
        SigmaPoints_0 ret_1;
        (&ret_1)->p_0[int(0)] = mean_5;
        (&ret_1)->w_mean_0[int(0)] = 0.0f;
        (&ret_1)->w_cov_0[int(0)] = 2.0f;
        float _S734 = (F32_sqrt((3.0f)));
        float3  delta_3 = make_float3 (_S734 * _S732.x) * _S733.rows[0U];
        float3  _S735 = mean_5 + delta_3;
        float3  _S736 = mean_5 - delta_3;
        float3  delta_4 = make_float3 (_S734 * _S732.y) * _S733.rows[1U];
        float3  _S737 = mean_5 + delta_4;
        float3  _S738 = mean_5 - delta_4;
        float3  delta_5 = make_float3 (_S734 * _S732.z) * _S733.rows[2U];
        float3  _S739 = mean_5 + delta_5;
        float3  _S740 = mean_5 - delta_5;
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
        (&ret_1)->p_0[0U] = mul_0(R_10, (&ret_1)->p_0[0U]) + t_7;
        (&ret_1)->p_0[1U] = mul_0(R_10, _S735) + t_7;
        (&ret_1)->p_0[2U] = mul_0(R_10, _S737) + t_7;
        (&ret_1)->p_0[3U] = mul_0(R_10, _S739) + t_7;
        (&ret_1)->p_0[4U] = mul_0(R_10, _S736) + t_7;
        (&ret_1)->p_0[5U] = mul_0(R_10, _S738) + t_7;
        (&ret_1)->p_0[6U] = mul_0(R_10, _S740) + t_7;
        SigmaPoints_0 _S741 = ret_1;
        bool _S742 = fisheye_proj_3dgs_ut_0(&_S741, intrins_8, dist_coeffs_17, &covar2d_4, mean2d_12);
        if(!(true & _S742))
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float det_orig_5 = *&(((&covar2d_4)->rows + (int(0)))->x) * *&(((&covar2d_4)->rows + (int(1)))->y) - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
        float _S743 = *&(((&covar2d_4)->rows + (int(0)))->x) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(0)))->x) = _S743;
        float _S744 = *&(((&covar2d_4)->rows + (int(1)))->y) + 0.30000001192092896f;
        *&(((&covar2d_4)->rows + (int(1)))->y) = _S744;
        float det_blur_5 = _S743 * _S744 - *&(((&covar2d_4)->rows + (int(0)))->y) * *&(((&covar2d_4)->rows + (int(1)))->x);
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
            _S731 = true;
        }
        else
        {
            _S731 = xmin_4 >= float(image_width_4);
        }
        if(_S731)
        {
            _S731 = true;
        }
        else
        {
            _S731 = ymax_4 <= 0.0f;
        }
        if(_S731)
        {
            _S731 = true;
        }
        else
        {
            _S731 = ymin_4 >= float(image_height_4);
        }
        if(_S731)
        {
            *aabb_xyxy_4 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_4 = make_int4 (int(xmin_4), int(ymin_4), int(xmax_4), int(ymax_4));
        *depth_4 = 0.5f * (F32_log((dot_0(mean_c_4, mean_c_4) + 9.99999997475242708e-07f)));
        *conic_4 = exp_0(- scale_7);
        float3  _S745 = mean_5 - - mul_0(transpose_0(R_10), t_7);
        float3  _S746 = make_float3 (0.282094806432724f) * (*sh_coeffs_4)[int(0)];
        *rgb_4 = _S746;
        float _S747 = _S745.x;
        float _S748 = _S745.y;
        float _S749 = _S745.z;
        float norm_4 = (F32_sqrt((_S747 * _S747 + _S748 * _S748 + _S749 * _S749)));
        float x_25 = _S747 / norm_4;
        float y_7 = _S748 / norm_4;
        float z_4 = _S749 / norm_4;
        float3  _S750 = _S746 + make_float3 (0.48860251903533936f) * (make_float3 (- y_7) * (*sh_coeffs_4)[int(1)] + make_float3 (z_4) * (*sh_coeffs_4)[int(2)] - make_float3 (x_25) * (*sh_coeffs_4)[int(3)]);
        *rgb_4 = _S750;
        float z2_13 = z_4 * z_4;
        float fTmp0B_4 = -1.09254848957061768f * z_4;
        float fC1_4 = x_25 * x_25 - y_7 * y_7;
        float fS1_4 = 2.0f * x_25 * y_7;
        float3  _S751 = _S750 + (make_float3 (0.54627424478530884f * fS1_4) * (*sh_coeffs_4)[int(4)] + make_float3 (fTmp0B_4 * y_7) * (*sh_coeffs_4)[int(5)] + make_float3 (0.94617468118667603f * z2_13 - 0.31539157032966614f) * (*sh_coeffs_4)[int(6)] + make_float3 (fTmp0B_4 * x_25) * (*sh_coeffs_4)[int(7)] + make_float3 (0.54627424478530884f * fC1_4) * (*sh_coeffs_4)[int(8)]);
        *rgb_4 = _S751;
        float fTmp0C_4 = -2.28522896766662598f * z2_13 + 0.4570457935333252f;
        float fTmp1B_4 = 1.44530570507049561f * z_4;
        *rgb_4 = max_0(_S751 + (make_float3 (-0.59004360437393188f * (x_25 * fS1_4 + y_7 * fC1_4)) * (*sh_coeffs_4)[int(9)] + make_float3 (fTmp1B_4 * fS1_4) * (*sh_coeffs_4)[int(10)] + make_float3 (fTmp0C_4 * y_7) * (*sh_coeffs_4)[int(11)] + make_float3 (z_4 * (1.86588168144226074f * z2_13 - 1.11952900886535645f)) * (*sh_coeffs_4)[int(12)] + make_float3 (fTmp0C_4 * x_25) * (*sh_coeffs_4)[int(13)] + make_float3 (fTmp1B_4 * fC1_4) * (*sh_coeffs_4)[int(14)] + make_float3 (-0.59004360437393188f * (x_25 * fC1_4 - y_7 * fS1_4)) * (*sh_coeffs_4)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_3dgs_persp_differentiable(bool antialiased_5, float3  mean_6, float4  quat_9, float3  scale_8, float in_opacity_5, FixedArray<float3 , 16>  * sh_coeffs_5, Matrix<float, 3, 3>  R_11, float3  t_8, float fx_12, float fy_12, float cx_9, float cy_9, FixedArray<float, 10>  * dist_coeffs_18, uint image_width_5, uint image_height_5, float near_plane_5, float far_plane_5, int4  * aabb_xyxy_5, float2  * mean2d_13, float * depth_5, float3  * conic_5, float * opacity_5, float3  * rgb_5)
{
    float3  mean_c_5 = mul_0(R_11, mean_6) + t_8;
    float3  _S752 = exp_0(scale_8);
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
    Matrix<float, 3, 3>  M_5 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_9 + z2_14), 2.0f * (xy_9 + wz_9), 2.0f * (xz_9 - wy_9), 2.0f * (xy_9 - wz_9), 1.0f - 2.0f * (x2_9 + z2_14), 2.0f * (yz_9 + wx_9), 2.0f * (xz_9 + wy_9), 2.0f * (yz_9 - wx_9), 1.0f - 2.0f * (x2_9 + y2_9))), makeMatrix<float, 3, 3> (_S752.x, 0.0f, 0.0f, 0.0f, _S752.y, 0.0f, 0.0f, 0.0f, _S752.z));
    Matrix<float, 3, 3>  _S753 = transpose_0(R_11);
    float _S754 = float(image_width_5);
    float _S755 = float(image_height_5);
    float _S756 = 0.30000001192092896f * (0.5f * _S754 / fx_12);
    float _S757 = 0.30000001192092896f * (0.5f * _S755 / fy_12);
    float rz_2 = 1.0f / mean_c_5.z;
    float rz2_2 = rz_2 * rz_2;
    Matrix<float, 2, 3>  J_6 = makeMatrix<float, 2, 3> (fx_12 * rz_2, 0.0f, - fx_12 * (mean_c_5.z * (F32_min(((_S754 - cx_9) / fx_12 + _S756), ((F32_max((- (cx_9 / fx_12 + _S756)), (mean_c_5.x * rz_2))))))) * rz2_2, 0.0f, fy_12 * rz_2, - fy_12 * (mean_c_5.z * (F32_min(((_S755 - cy_9) / fy_12 + _S757), ((F32_max((- (cy_9 / fy_12 + _S757)), (mean_c_5.y * rz_2))))))) * rz2_2);
    Matrix<float, 2, 2>  covar2d_5 = mul_6(mul_5(J_6, mul_4(mul_4(R_11, mul_4(M_5, transpose_0(M_5))), _S753)), transpose_1(J_6));
    *mean2d_13 = make_float2 (fx_12 * mean_c_5.x * rz_2 + cx_9, fy_12 * mean_c_5.y * rz_2 + cy_9);
    float det_orig_6 = *&(((&covar2d_5)->rows + (int(0)))->x) * *&(((&covar2d_5)->rows + (int(1)))->y) - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x);
    float _S758 = *&(((&covar2d_5)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(0)))->x) = _S758;
    float _S759 = *&(((&covar2d_5)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_5)->rows + (int(1)))->y) = _S759;
    float compensation_6 = (F32_sqrt(((F32_max((0.0f), (det_orig_6 / (_S758 * _S759 - *&(((&covar2d_5)->rows + (int(0)))->y) * *&(((&covar2d_5)->rows + (int(1)))->x))))))));
    float invdet_4 = 1.0f / (covar2d_5.rows[int(0)].x * covar2d_5.rows[int(1)].y - covar2d_5.rows[int(0)].y * covar2d_5.rows[int(1)].x);
    Matrix<float, 2, 2>  _S760 = makeMatrix<float, 2, 2> (covar2d_5.rows[int(1)].y * invdet_4, - covar2d_5.rows[int(0)].y * invdet_4, - covar2d_5.rows[int(1)].x * invdet_4, covar2d_5.rows[int(0)].x * invdet_4);
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
    *conic_5 = make_float3 (_S760.rows[int(0)].x, _S760.rows[int(0)].y, _S760.rows[int(1)].y);
    float3  _S761 = mean_6 - - mul_0(_S753, t_8);
    float3  _S762 = make_float3 (0.282094806432724f) * (*sh_coeffs_5)[int(0)];
    *rgb_5 = _S762;
    float _S763 = _S761.x;
    float _S764 = _S761.y;
    float _S765 = _S761.z;
    float norm_5 = (F32_sqrt((_S763 * _S763 + _S764 * _S764 + _S765 * _S765)));
    float x_27 = _S763 / norm_5;
    float y_8 = _S764 / norm_5;
    float z_5 = _S765 / norm_5;
    float3  _S766 = _S762 + make_float3 (0.48860251903533936f) * (make_float3 (- y_8) * (*sh_coeffs_5)[int(1)] + make_float3 (z_5) * (*sh_coeffs_5)[int(2)] - make_float3 (x_27) * (*sh_coeffs_5)[int(3)]);
    *rgb_5 = _S766;
    float z2_15 = z_5 * z_5;
    float fTmp0B_5 = -1.09254848957061768f * z_5;
    float fC1_5 = x_27 * x_27 - y_8 * y_8;
    float fS1_5 = 2.0f * x_27 * y_8;
    float3  _S767 = _S766 + (make_float3 (0.54627424478530884f * fS1_5) * (*sh_coeffs_5)[int(4)] + make_float3 (fTmp0B_5 * y_8) * (*sh_coeffs_5)[int(5)] + make_float3 (0.94617468118667603f * z2_15 - 0.31539157032966614f) * (*sh_coeffs_5)[int(6)] + make_float3 (fTmp0B_5 * x_27) * (*sh_coeffs_5)[int(7)] + make_float3 (0.54627424478530884f * fC1_5) * (*sh_coeffs_5)[int(8)]);
    *rgb_5 = _S767;
    float fTmp0C_5 = -2.28522896766662598f * z2_15 + 0.4570457935333252f;
    float fTmp1B_5 = 1.44530570507049561f * z_5;
    *rgb_5 = max_0(_S767 + (make_float3 (-0.59004360437393188f * (x_27 * fS1_5 + y_8 * fC1_5)) * (*sh_coeffs_5)[int(9)] + make_float3 (fTmp1B_5 * fS1_5) * (*sh_coeffs_5)[int(10)] + make_float3 (fTmp0C_5 * y_8) * (*sh_coeffs_5)[int(11)] + make_float3 (z_5 * (1.86588168144226074f * z2_15 - 1.11952900886535645f)) * (*sh_coeffs_5)[int(12)] + make_float3 (fTmp0C_5 * x_27) * (*sh_coeffs_5)[int(13)] + make_float3 (fTmp1B_5 * fC1_5) * (*sh_coeffs_5)[int(14)] + make_float3 (-0.59004360437393188f * (x_27 * fC1_5 - y_8 * fS1_5)) * (*sh_coeffs_5)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_fisheye_differentiable(bool antialiased_6, float3  mean_7, float4  quat_10, float3  scale_9, float in_opacity_6, FixedArray<float3 , 16>  * sh_coeffs_6, Matrix<float, 3, 3>  R_12, float3  t_9, float fx_13, float fy_13, float cx_10, float cy_10, FixedArray<float, 10>  * dist_coeffs_19, uint image_width_6, uint image_height_6, float near_plane_6, float far_plane_6, int4  * aabb_xyxy_6, float2  * mean2d_14, float * depth_6, float3  * conic_6, float * opacity_6, float3  * rgb_6)
{
    float3  mean_c_6 = mul_0(R_12, mean_7) + t_9;
    float3  _S768 = exp_0(scale_9);
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
    Matrix<float, 3, 3>  M_6 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_10 + z2_16), 2.0f * (xy_10 + wz_10), 2.0f * (xz_10 - wy_10), 2.0f * (xy_10 - wz_10), 1.0f - 2.0f * (x2_10 + z2_16), 2.0f * (yz_10 + wx_10), 2.0f * (xz_10 + wy_10), 2.0f * (yz_10 - wx_10), 1.0f - 2.0f * (x2_10 + y2_10))), makeMatrix<float, 3, 3> (_S768.x, 0.0f, 0.0f, 0.0f, _S768.y, 0.0f, 0.0f, 0.0f, _S768.z));
    Matrix<float, 3, 3>  _S769 = transpose_0(R_12);
    Matrix<float, 2, 2>  covar2d_6;
    bool _S770 = fisheye_proj_3dgs_1(mean_c_6, mul_4(mul_4(R_12, mul_4(M_6, transpose_0(M_6))), _S769), make_float4 (fx_13, fy_13, cx_10, cy_10), dist_coeffs_19, &covar2d_6, mean2d_14);
    float det_orig_7 = *&(((&covar2d_6)->rows + (int(0)))->x) * *&(((&covar2d_6)->rows + (int(1)))->y) - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x);
    float _S771 = *&(((&covar2d_6)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(0)))->x) = _S771;
    float _S772 = *&(((&covar2d_6)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_6)->rows + (int(1)))->y) = _S772;
    float compensation_7 = (F32_sqrt(((F32_max((0.0f), (det_orig_7 / (_S771 * _S772 - *&(((&covar2d_6)->rows + (int(0)))->y) * *&(((&covar2d_6)->rows + (int(1)))->x))))))));
    float invdet_5 = 1.0f / (covar2d_6.rows[int(0)].x * covar2d_6.rows[int(1)].y - covar2d_6.rows[int(0)].y * covar2d_6.rows[int(1)].x);
    Matrix<float, 2, 2>  _S773 = makeMatrix<float, 2, 2> (covar2d_6.rows[int(1)].y * invdet_5, - covar2d_6.rows[int(0)].y * invdet_5, - covar2d_6.rows[int(1)].x * invdet_5, covar2d_6.rows[int(0)].x * invdet_5);
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
    *conic_6 = make_float3 (_S773.rows[int(0)].x, _S773.rows[int(0)].y, _S773.rows[int(1)].y);
    float3  _S774 = mean_7 - - mul_0(_S769, t_9);
    float3  _S775 = make_float3 (0.282094806432724f) * (*sh_coeffs_6)[int(0)];
    *rgb_6 = _S775;
    float _S776 = _S774.x;
    float _S777 = _S774.y;
    float _S778 = _S774.z;
    float norm_6 = (F32_sqrt((_S776 * _S776 + _S777 * _S777 + _S778 * _S778)));
    float x_29 = _S776 / norm_6;
    float y_9 = _S777 / norm_6;
    float z_6 = _S778 / norm_6;
    float3  _S779 = _S775 + make_float3 (0.48860251903533936f) * (make_float3 (- y_9) * (*sh_coeffs_6)[int(1)] + make_float3 (z_6) * (*sh_coeffs_6)[int(2)] - make_float3 (x_29) * (*sh_coeffs_6)[int(3)]);
    *rgb_6 = _S779;
    float z2_17 = z_6 * z_6;
    float fTmp0B_6 = -1.09254848957061768f * z_6;
    float fC1_6 = x_29 * x_29 - y_9 * y_9;
    float fS1_6 = 2.0f * x_29 * y_9;
    float3  _S780 = _S779 + (make_float3 (0.54627424478530884f * fS1_6) * (*sh_coeffs_6)[int(4)] + make_float3 (fTmp0B_6 * y_9) * (*sh_coeffs_6)[int(5)] + make_float3 (0.94617468118667603f * z2_17 - 0.31539157032966614f) * (*sh_coeffs_6)[int(6)] + make_float3 (fTmp0B_6 * x_29) * (*sh_coeffs_6)[int(7)] + make_float3 (0.54627424478530884f * fC1_6) * (*sh_coeffs_6)[int(8)]);
    *rgb_6 = _S780;
    float fTmp0C_6 = -2.28522896766662598f * z2_17 + 0.4570457935333252f;
    float fTmp1B_6 = 1.44530570507049561f * z_6;
    *rgb_6 = max_0(_S780 + (make_float3 (-0.59004360437393188f * (x_29 * fS1_6 + y_9 * fC1_6)) * (*sh_coeffs_6)[int(9)] + make_float3 (fTmp1B_6 * fS1_6) * (*sh_coeffs_6)[int(10)] + make_float3 (fTmp0C_6 * y_9) * (*sh_coeffs_6)[int(11)] + make_float3 (z_6 * (1.86588168144226074f * z2_17 - 1.11952900886535645f)) * (*sh_coeffs_6)[int(12)] + make_float3 (fTmp0C_6 * x_29) * (*sh_coeffs_6)[int(13)] + make_float3 (fTmp1B_6 * fC1_6) * (*sh_coeffs_6)[int(14)] + make_float3 (-0.59004360437393188f * (x_29 * fC1_6 - y_9 * fS1_6)) * (*sh_coeffs_6)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_ortho_differentiable(bool antialiased_7, float3  mean_8, float4  quat_11, float3  scale_10, float in_opacity_7, FixedArray<float3 , 16>  * sh_coeffs_7, Matrix<float, 3, 3>  R_13, float3  t_10, float fx_14, float fy_14, float cx_11, float cy_11, FixedArray<float, 10>  * dist_coeffs_20, uint image_width_7, uint image_height_7, float near_plane_7, float far_plane_7, int4  * aabb_xyxy_7, float2  * mean2d_15, float * depth_7, float3  * conic_7, float * opacity_7, float3  * rgb_7)
{
    float3  mean_c_7 = mul_0(R_13, mean_8) + t_10;
    float3  _S781 = exp_0(scale_10);
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
    Matrix<float, 3, 3>  M_7 = mul_4(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_11 + z2_18), 2.0f * (xy_11 + wz_11), 2.0f * (xz_11 - wy_11), 2.0f * (xy_11 - wz_11), 1.0f - 2.0f * (x2_11 + z2_18), 2.0f * (yz_11 + wx_11), 2.0f * (xz_11 + wy_11), 2.0f * (yz_11 - wx_11), 1.0f - 2.0f * (x2_11 + y2_11))), makeMatrix<float, 3, 3> (_S781.x, 0.0f, 0.0f, 0.0f, _S781.y, 0.0f, 0.0f, 0.0f, _S781.z));
    Matrix<float, 3, 3>  _S782 = transpose_0(R_13);
    Matrix<float, 2, 3>  J_7 = makeMatrix<float, 2, 3> (fx_14, 0.0f, 0.0f, 0.0f, fy_14, 0.0f);
    Matrix<float, 2, 2>  covar2d_7 = mul_6(mul_5(J_7, mul_4(mul_4(R_13, mul_4(M_7, transpose_0(M_7))), _S782)), transpose_1(J_7));
    *mean2d_15 = make_float2 (fx_14 * mean_c_7.x + cx_11, fy_14 * mean_c_7.y + cy_11);
    float det_orig_8 = *&(((&covar2d_7)->rows + (int(0)))->x) * *&(((&covar2d_7)->rows + (int(1)))->y) - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x);
    float _S783 = *&(((&covar2d_7)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(0)))->x) = _S783;
    float _S784 = *&(((&covar2d_7)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_7)->rows + (int(1)))->y) = _S784;
    float compensation_8 = (F32_sqrt(((F32_max((0.0f), (det_orig_8 / (_S783 * _S784 - *&(((&covar2d_7)->rows + (int(0)))->y) * *&(((&covar2d_7)->rows + (int(1)))->x))))))));
    float invdet_6 = 1.0f / (covar2d_7.rows[int(0)].x * covar2d_7.rows[int(1)].y - covar2d_7.rows[int(0)].y * covar2d_7.rows[int(1)].x);
    Matrix<float, 2, 2>  _S785 = makeMatrix<float, 2, 2> (covar2d_7.rows[int(1)].y * invdet_6, - covar2d_7.rows[int(0)].y * invdet_6, - covar2d_7.rows[int(1)].x * invdet_6, covar2d_7.rows[int(0)].x * invdet_6);
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
    *conic_7 = make_float3 (_S785.rows[int(0)].x, _S785.rows[int(0)].y, _S785.rows[int(1)].y);
    float3  _S786 = mean_8 - - mul_0(_S782, t_10);
    float3  _S787 = make_float3 (0.282094806432724f) * (*sh_coeffs_7)[int(0)];
    *rgb_7 = _S787;
    float _S788 = _S786.x;
    float _S789 = _S786.y;
    float _S790 = _S786.z;
    float norm_7 = (F32_sqrt((_S788 * _S788 + _S789 * _S789 + _S790 * _S790)));
    float x_31 = _S788 / norm_7;
    float y_10 = _S789 / norm_7;
    float z_7 = _S790 / norm_7;
    float3  _S791 = _S787 + make_float3 (0.48860251903533936f) * (make_float3 (- y_10) * (*sh_coeffs_7)[int(1)] + make_float3 (z_7) * (*sh_coeffs_7)[int(2)] - make_float3 (x_31) * (*sh_coeffs_7)[int(3)]);
    *rgb_7 = _S791;
    float z2_19 = z_7 * z_7;
    float fTmp0B_7 = -1.09254848957061768f * z_7;
    float fC1_7 = x_31 * x_31 - y_10 * y_10;
    float fS1_7 = 2.0f * x_31 * y_10;
    float3  _S792 = _S791 + (make_float3 (0.54627424478530884f * fS1_7) * (*sh_coeffs_7)[int(4)] + make_float3 (fTmp0B_7 * y_10) * (*sh_coeffs_7)[int(5)] + make_float3 (0.94617468118667603f * z2_19 - 0.31539157032966614f) * (*sh_coeffs_7)[int(6)] + make_float3 (fTmp0B_7 * x_31) * (*sh_coeffs_7)[int(7)] + make_float3 (0.54627424478530884f * fC1_7) * (*sh_coeffs_7)[int(8)]);
    *rgb_7 = _S792;
    float fTmp0C_7 = -2.28522896766662598f * z2_19 + 0.4570457935333252f;
    float fTmp1B_7 = 1.44530570507049561f * z_7;
    *rgb_7 = max_0(_S792 + (make_float3 (-0.59004360437393188f * (x_31 * fS1_7 + y_10 * fC1_7)) * (*sh_coeffs_7)[int(9)] + make_float3 (fTmp1B_7 * fS1_7) * (*sh_coeffs_7)[int(10)] + make_float3 (fTmp0C_7 * y_10) * (*sh_coeffs_7)[int(11)] + make_float3 (z_7 * (1.86588168144226074f * z2_19 - 1.11952900886535645f)) * (*sh_coeffs_7)[int(12)] + make_float3 (fTmp0C_7 * x_31) * (*sh_coeffs_7)[int(13)] + make_float3 (fTmp1B_7 * fC1_7) * (*sh_coeffs_7)[int(14)] + make_float3 (-0.59004360437393188f * (x_31 * fC1_7 - y_10 * fS1_7)) * (*sh_coeffs_7)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_persp_differentiable(bool antialiased_8, float3  mean_9, float4  quat_12, float3  scale_11, float in_opacity_8, FixedArray<float3 , 16>  * sh_coeffs_8, Matrix<float, 3, 3>  R_14, float3  t_11, float fx_15, float fy_15, float cx_12, float cy_12, FixedArray<float, 10>  * dist_coeffs_21, uint image_width_8, uint image_height_8, float near_plane_8, float far_plane_8, int4  * aabb_xyxy_8, float2  * mean2d_16, float * depth_8, float3  * conic_8, float * opacity_8, float3  * rgb_8)
{
    float3  mean_c_8 = mul_0(R_14, mean_9) + t_11;
    float4  intrins_9 = make_float4 (fx_15, fy_15, cx_12, cy_12);
    float3  _S793 = exp_0(scale_11);
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
    Matrix<float, 3, 3>  _S794 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_12 + z2_20), 2.0f * (xy_12 + wz_12), 2.0f * (xz_12 - wy_12), 2.0f * (xy_12 - wz_12), 1.0f - 2.0f * (x2_12 + z2_20), 2.0f * (yz_12 + wx_12), 2.0f * (xz_12 + wy_12), 2.0f * (yz_12 - wx_12), 1.0f - 2.0f * (x2_12 + y2_12))));
    SigmaPoints_0 ret_2;
    (&ret_2)->p_0[int(0)] = mean_9;
    (&ret_2)->w_mean_0[int(0)] = 0.0f;
    (&ret_2)->w_cov_0[int(0)] = 2.0f;
    float _S795 = (F32_sqrt((3.0f)));
    float3  delta_6 = make_float3 (_S795 * _S793.x) * _S794.rows[0U];
    float3  _S796 = mean_9 + delta_6;
    float3  _S797 = mean_9 - delta_6;
    float3  delta_7 = make_float3 (_S795 * _S793.y) * _S794.rows[1U];
    float3  _S798 = mean_9 + delta_7;
    float3  _S799 = mean_9 - delta_7;
    float3  delta_8 = make_float3 (_S795 * _S793.z) * _S794.rows[2U];
    float3  _S800 = mean_9 + delta_8;
    float3  _S801 = mean_9 - delta_8;
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
    (&ret_2)->p_0[0U] = mul_0(R_14, (&ret_2)->p_0[0U]) + t_11;
    (&ret_2)->p_0[1U] = mul_0(R_14, _S796) + t_11;
    (&ret_2)->p_0[2U] = mul_0(R_14, _S798) + t_11;
    (&ret_2)->p_0[3U] = mul_0(R_14, _S800) + t_11;
    (&ret_2)->p_0[4U] = mul_0(R_14, _S797) + t_11;
    (&ret_2)->p_0[5U] = mul_0(R_14, _S799) + t_11;
    (&ret_2)->p_0[6U] = mul_0(R_14, _S801) + t_11;
    SigmaPoints_0 _S802 = ret_2;
    Matrix<float, 2, 2>  covar2d_8;
    bool _S803 = persp_proj_3dgs_ut_1(&_S802, intrins_9, dist_coeffs_21, image_width_8, image_height_8, &covar2d_8, mean2d_16);
    float det_orig_9 = *&(((&covar2d_8)->rows + (int(0)))->x) * *&(((&covar2d_8)->rows + (int(1)))->y) - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x);
    float _S804 = *&(((&covar2d_8)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(0)))->x) = _S804;
    float _S805 = *&(((&covar2d_8)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_8)->rows + (int(1)))->y) = _S805;
    float compensation_9 = (F32_sqrt(((F32_max((0.0f), (det_orig_9 / (_S804 * _S805 - *&(((&covar2d_8)->rows + (int(0)))->y) * *&(((&covar2d_8)->rows + (int(1)))->x))))))));
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
    float3  _S806 = mean_9 - - mul_0(transpose_0(R_14), t_11);
    float3  _S807 = make_float3 (0.282094806432724f) * (*sh_coeffs_8)[int(0)];
    *rgb_8 = _S807;
    float _S808 = _S806.x;
    float _S809 = _S806.y;
    float _S810 = _S806.z;
    float norm_8 = (F32_sqrt((_S808 * _S808 + _S809 * _S809 + _S810 * _S810)));
    float x_33 = _S808 / norm_8;
    float y_11 = _S809 / norm_8;
    float z_8 = _S810 / norm_8;
    float3  _S811 = _S807 + make_float3 (0.48860251903533936f) * (make_float3 (- y_11) * (*sh_coeffs_8)[int(1)] + make_float3 (z_8) * (*sh_coeffs_8)[int(2)] - make_float3 (x_33) * (*sh_coeffs_8)[int(3)]);
    *rgb_8 = _S811;
    float z2_21 = z_8 * z_8;
    float fTmp0B_8 = -1.09254848957061768f * z_8;
    float fC1_8 = x_33 * x_33 - y_11 * y_11;
    float fS1_8 = 2.0f * x_33 * y_11;
    float3  _S812 = _S811 + (make_float3 (0.54627424478530884f * fS1_8) * (*sh_coeffs_8)[int(4)] + make_float3 (fTmp0B_8 * y_11) * (*sh_coeffs_8)[int(5)] + make_float3 (0.94617468118667603f * z2_21 - 0.31539157032966614f) * (*sh_coeffs_8)[int(6)] + make_float3 (fTmp0B_8 * x_33) * (*sh_coeffs_8)[int(7)] + make_float3 (0.54627424478530884f * fC1_8) * (*sh_coeffs_8)[int(8)]);
    *rgb_8 = _S812;
    float fTmp0C_8 = -2.28522896766662598f * z2_21 + 0.4570457935333252f;
    float fTmp1B_8 = 1.44530570507049561f * z_8;
    *rgb_8 = max_0(_S812 + (make_float3 (-0.59004360437393188f * (x_33 * fS1_8 + y_11 * fC1_8)) * (*sh_coeffs_8)[int(9)] + make_float3 (fTmp1B_8 * fS1_8) * (*sh_coeffs_8)[int(10)] + make_float3 (fTmp0C_8 * y_11) * (*sh_coeffs_8)[int(11)] + make_float3 (z_8 * (1.86588168144226074f * z2_21 - 1.11952900886535645f)) * (*sh_coeffs_8)[int(12)] + make_float3 (fTmp0C_8 * x_33) * (*sh_coeffs_8)[int(13)] + make_float3 (fTmp1B_8 * fC1_8) * (*sh_coeffs_8)[int(14)] + make_float3 (-0.59004360437393188f * (x_33 * fC1_8 - y_11 * fS1_8)) * (*sh_coeffs_8)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_3dgs_eval3d_fisheye_differentiable(bool antialiased_9, float3  mean_10, float4  quat_13, float3  scale_12, float in_opacity_9, FixedArray<float3 , 16>  * sh_coeffs_9, Matrix<float, 3, 3>  R_15, float3  t_12, float fx_16, float fy_16, float cx_13, float cy_13, FixedArray<float, 10>  * dist_coeffs_22, uint image_width_9, uint image_height_9, float near_plane_9, float far_plane_9, int4  * aabb_xyxy_9, float2  * mean2d_17, float * depth_9, float3  * conic_9, float * opacity_9, float3  * rgb_9)
{
    float3  mean_c_9 = mul_0(R_15, mean_10) + t_12;
    float4  intrins_10 = make_float4 (fx_16, fy_16, cx_13, cy_13);
    float3  _S813 = exp_0(scale_12);
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
    Matrix<float, 3, 3>  _S814 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_13 + z2_22), 2.0f * (xy_13 + wz_13), 2.0f * (xz_13 - wy_13), 2.0f * (xy_13 - wz_13), 1.0f - 2.0f * (x2_13 + z2_22), 2.0f * (yz_13 + wx_13), 2.0f * (xz_13 + wy_13), 2.0f * (yz_13 - wx_13), 1.0f - 2.0f * (x2_13 + y2_13))));
    SigmaPoints_0 ret_3;
    (&ret_3)->p_0[int(0)] = mean_10;
    (&ret_3)->w_mean_0[int(0)] = 0.0f;
    (&ret_3)->w_cov_0[int(0)] = 2.0f;
    float _S815 = (F32_sqrt((3.0f)));
    float3  delta_9 = make_float3 (_S815 * _S813.x) * _S814.rows[0U];
    float3  _S816 = mean_10 + delta_9;
    float3  _S817 = mean_10 - delta_9;
    float3  delta_10 = make_float3 (_S815 * _S813.y) * _S814.rows[1U];
    float3  _S818 = mean_10 + delta_10;
    float3  _S819 = mean_10 - delta_10;
    float3  delta_11 = make_float3 (_S815 * _S813.z) * _S814.rows[2U];
    float3  _S820 = mean_10 + delta_11;
    float3  _S821 = mean_10 - delta_11;
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
    (&ret_3)->p_0[0U] = mul_0(R_15, (&ret_3)->p_0[0U]) + t_12;
    (&ret_3)->p_0[1U] = mul_0(R_15, _S816) + t_12;
    (&ret_3)->p_0[2U] = mul_0(R_15, _S818) + t_12;
    (&ret_3)->p_0[3U] = mul_0(R_15, _S820) + t_12;
    (&ret_3)->p_0[4U] = mul_0(R_15, _S817) + t_12;
    (&ret_3)->p_0[5U] = mul_0(R_15, _S819) + t_12;
    (&ret_3)->p_0[6U] = mul_0(R_15, _S821) + t_12;
    SigmaPoints_0 _S822 = ret_3;
    Matrix<float, 2, 2>  covar2d_9;
    bool _S823 = fisheye_proj_3dgs_ut_1(&_S822, intrins_10, dist_coeffs_22, &covar2d_9, mean2d_17);
    float det_orig_10 = *&(((&covar2d_9)->rows + (int(0)))->x) * *&(((&covar2d_9)->rows + (int(1)))->y) - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x);
    float _S824 = *&(((&covar2d_9)->rows + (int(0)))->x) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(0)))->x) = _S824;
    float _S825 = *&(((&covar2d_9)->rows + (int(1)))->y) + 0.30000001192092896f;
    *&(((&covar2d_9)->rows + (int(1)))->y) = _S825;
    float compensation_10 = (F32_sqrt(((F32_max((0.0f), (det_orig_10 / (_S824 * _S825 - *&(((&covar2d_9)->rows + (int(0)))->y) * *&(((&covar2d_9)->rows + (int(1)))->x))))))));
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
    float3  _S826 = mean_10 - - mul_0(transpose_0(R_15), t_12);
    float3  _S827 = make_float3 (0.282094806432724f) * (*sh_coeffs_9)[int(0)];
    *rgb_9 = _S827;
    float _S828 = _S826.x;
    float _S829 = _S826.y;
    float _S830 = _S826.z;
    float norm_9 = (F32_sqrt((_S828 * _S828 + _S829 * _S829 + _S830 * _S830)));
    float x_35 = _S828 / norm_9;
    float y_12 = _S829 / norm_9;
    float z_9 = _S830 / norm_9;
    float3  _S831 = _S827 + make_float3 (0.48860251903533936f) * (make_float3 (- y_12) * (*sh_coeffs_9)[int(1)] + make_float3 (z_9) * (*sh_coeffs_9)[int(2)] - make_float3 (x_35) * (*sh_coeffs_9)[int(3)]);
    *rgb_9 = _S831;
    float z2_23 = z_9 * z_9;
    float fTmp0B_9 = -1.09254848957061768f * z_9;
    float fC1_9 = x_35 * x_35 - y_12 * y_12;
    float fS1_9 = 2.0f * x_35 * y_12;
    float3  _S832 = _S831 + (make_float3 (0.54627424478530884f * fS1_9) * (*sh_coeffs_9)[int(4)] + make_float3 (fTmp0B_9 * y_12) * (*sh_coeffs_9)[int(5)] + make_float3 (0.94617468118667603f * z2_23 - 0.31539157032966614f) * (*sh_coeffs_9)[int(6)] + make_float3 (fTmp0B_9 * x_35) * (*sh_coeffs_9)[int(7)] + make_float3 (0.54627424478530884f * fC1_9) * (*sh_coeffs_9)[int(8)]);
    *rgb_9 = _S832;
    float fTmp0C_9 = -2.28522896766662598f * z2_23 + 0.4570457935333252f;
    float fTmp1B_9 = 1.44530570507049561f * z_9;
    *rgb_9 = max_0(_S832 + (make_float3 (-0.59004360437393188f * (x_35 * fS1_9 + y_12 * fC1_9)) * (*sh_coeffs_9)[int(9)] + make_float3 (fTmp1B_9 * fS1_9) * (*sh_coeffs_9)[int(10)] + make_float3 (fTmp0C_9 * y_12) * (*sh_coeffs_9)[int(11)] + make_float3 (z_9 * (1.86588168144226074f * z2_23 - 1.11952900886535645f)) * (*sh_coeffs_9)[int(12)] + make_float3 (fTmp0C_9 * x_35) * (*sh_coeffs_9)[int(13)] + make_float3 (fTmp1B_9 * fC1_9) * (*sh_coeffs_9)[int(14)] + make_float3 (-0.59004360437393188f * (x_35 * fC1_9 - y_12 * fS1_9)) * (*sh_coeffs_9)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ float3  s_primal_ctx_mul_0(Matrix<float, 3, 3>  _S833, float3  _S834)
{
    return mul_0(_S833, _S834);
}

inline __device__ float3  s_primal_ctx_exp_0(float3  _S835)
{
    return exp_0(_S835);
}

inline __device__ Matrix<float, 3, 3>  s_primal_ctx_mul_1(Matrix<float, 3, 3>  _S836, Matrix<float, 3, 3>  _S837)
{
    return mul_4(_S836, _S837);
}

inline __device__ float s_primal_ctx_max_0(float _S838, float _S839)
{
    return (F32_max((_S838), (_S839)));
}

inline __device__ float s_primal_ctx_min_0(float _S840, float _S841)
{
    return (F32_min((_S840), (_S841)));
}

inline __device__ Matrix<float, 2, 3>  s_primal_ctx_mul_2(Matrix<float, 2, 3>  _S842, Matrix<float, 3, 3>  _S843)
{
    return mul_5(_S842, _S843);
}

inline __device__ Matrix<float, 2, 2>  s_primal_ctx_mul_3(Matrix<float, 2, 3>  _S844, Matrix<float, 3, 2>  _S845)
{
    return mul_6(_S844, _S845);
}

inline __device__ float s_primal_ctx_sqrt_0(float _S846)
{
    return (F32_sqrt((_S846)));
}

inline __device__ float s_primal_ctx_exp_1(float _S847)
{
    return (F32_exp((_S847)));
}

inline __device__ float s_primal_ctx_log_0(float _S848)
{
    return (F32_log((_S848)));
}

inline __device__ float s_primal_ctx_dot_0(float3  _S849, float3  _S850)
{
    return dot_0(_S849, _S850);
}

inline __device__ void s_bwd_prop_max_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S851, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S852, float3  _S853)
{
    _d_max_vector_0(_S851, _S852, _S853);
    return;
}

inline __device__ void s_bwd_prop_mul_1(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S854, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S855, float3  _S856)
{
    _d_mul_0(_S854, _S855, _S856);
    return;
}

inline __device__ void s_bwd_prop_log_0(DiffPair_float_0 * _S857, float _S858)
{
    _d_log_0(_S857, _S858);
    return;
}

inline __device__ void s_bwd_prop_dot_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S859, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S860, float _S861)
{
    _d_dot_0(_S859, _S860, _S861);
    return;
}

inline __device__ void s_bwd_prop_min_0(DiffPair_float_0 * _S862, DiffPair_float_0 * _S863, float _S864)
{
    _d_min_0(_S862, _S863, _S864);
    return;
}

inline __device__ void s_bwd_prop_exp_0(DiffPair_float_0 * _S865, float _S866)
{
    _d_exp_0(_S865, _S866);
    return;
}

inline __device__ void s_bwd_prop_max_1(DiffPair_float_0 * _S867, DiffPair_float_0 * _S868, float _S869)
{
    _d_max_0(_S867, _S868, _S869);
    return;
}

inline __device__ void s_bwd_prop_mul_2(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S870, DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 * _S871, Matrix<float, 2, 2>  _S872)
{
    mul_3(_S870, _S871, _S872);
    return;
}

inline __device__ void s_bwd_prop_mul_3(DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 * _S873, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S874, Matrix<float, 2, 3>  _S875)
{
    mul_2(_S873, _S874, _S875);
    return;
}

inline __device__ void s_bwd_prop_mul_4(DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S876, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S877, Matrix<float, 3, 3>  _S878)
{
    mul_1(_S876, _S877, _S878);
    return;
}

inline __device__ void s_bwd_prop_exp_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S879, float3  _S880)
{
    _d_exp_vector_0(_S879, _S880);
    return;
}

inline __device__ void projection_3dgs_persp_vjp(bool antialiased_10, float3  mean_11, float4  quat_14, float3  scale_13, float in_opacity_10, FixedArray<float3 , 16>  * sh_coeffs_10, Matrix<float, 3, 3>  R_16, float3  t_13, float fx_17, float fy_17, float cx_14, float cy_14, FixedArray<float, 10>  * dist_coeffs_23, uint image_width_10, uint image_height_10, float2  v_mean2d_0, float v_depth_0, float3  v_conic_0, float v_opacity_0, float3  v_rgb_0, float3  * v_mean_0, float4  * v_quat_0, float3  * v_scale_0, float * v_in_opacity_0, FixedArray<float3 , 16>  * v_sh_coeffs_0, Matrix<float, 3, 3>  * v_R_2, float3  * v_t_1)
{
    float3  mean_c_10 = s_primal_ctx_mul_0(R_16, mean_11) + t_13;
    float3  _S881 = s_primal_ctx_exp_0(scale_13);
    float _S882 = quat_14.y;
    float x2_14 = _S882 * _S882;
    float y2_14 = quat_14.z * quat_14.z;
    float z2_24 = quat_14.w * quat_14.w;
    float xy_14 = quat_14.y * quat_14.z;
    float xz_14 = quat_14.y * quat_14.w;
    float yz_14 = quat_14.z * quat_14.w;
    float wx_14 = quat_14.x * quat_14.y;
    float wy_14 = quat_14.x * quat_14.z;
    float wz_14 = quat_14.x * quat_14.w;
    Matrix<float, 3, 3>  _S883 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_14 + z2_24), 2.0f * (xy_14 + wz_14), 2.0f * (xz_14 - wy_14), 2.0f * (xy_14 - wz_14), 1.0f - 2.0f * (x2_14 + z2_24), 2.0f * (yz_14 + wx_14), 2.0f * (xz_14 + wy_14), 2.0f * (yz_14 - wx_14), 1.0f - 2.0f * (x2_14 + y2_14)));
    Matrix<float, 3, 3>  S_0 = makeMatrix<float, 3, 3> (_S881.x, 0.0f, 0.0f, 0.0f, _S881.y, 0.0f, 0.0f, 0.0f, _S881.z);
    Matrix<float, 3, 3>  _S884 = s_primal_ctx_mul_1(_S883, S_0);
    Matrix<float, 3, 3>  _S885 = transpose_0(_S884);
    Matrix<float, 3, 3>  _S886 = s_primal_ctx_mul_1(_S884, _S885);
    Matrix<float, 3, 3>  _S887 = s_primal_ctx_mul_1(R_16, _S886);
    Matrix<float, 3, 3>  _S888 = transpose_0(R_16);
    Matrix<float, 3, 3>  _S889 = s_primal_ctx_mul_1(_S887, _S888);
    float _S890 = float(image_width_10);
    float _S891 = float(image_height_10);
    float _S892 = 0.30000001192092896f * (0.5f * _S890 / fx_17);
    float lim_x_pos_2 = (_S890 - cx_14) / fx_17 + _S892;
    float _S893 = 0.30000001192092896f * (0.5f * _S891 / fy_17);
    float lim_y_pos_2 = (_S891 - cy_14) / fy_17 + _S893;
    float rz_3 = 1.0f / mean_c_10.z;
    float _S894 = mean_c_10.z * mean_c_10.z;
    float rz2_3 = rz_3 * rz_3;
    float _S895 = - (cx_14 / fx_17 + _S892);
    float _S896 = mean_c_10.x * rz_3;
    float _S897 = s_primal_ctx_max_0(_S895, _S896);
    float _S898 = s_primal_ctx_min_0(lim_x_pos_2, _S897);
    float _S899 = - (cy_14 / fy_17 + _S893);
    float _S900 = mean_c_10.y * rz_3;
    float _S901 = s_primal_ctx_max_0(_S899, _S900);
    float _S902 = s_primal_ctx_min_0(lim_y_pos_2, _S901);
    float _S903 = - fx_17;
    float _S904 = _S903 * (mean_c_10.z * _S898);
    float _S905 = - fy_17;
    float _S906 = _S905 * (mean_c_10.z * _S902);
    Matrix<float, 2, 3>  J_8 = makeMatrix<float, 2, 3> (fx_17 * rz_3, 0.0f, _S904 * rz2_3, 0.0f, fy_17 * rz_3, _S906 * rz2_3);
    Matrix<float, 2, 3>  _S907 = s_primal_ctx_mul_2(J_8, _S889);
    Matrix<float, 3, 2>  _S908 = transpose_1(J_8);
    Matrix<float, 2, 2>  _S909 = s_primal_ctx_mul_3(_S907, _S908);
    float _S910 = fx_17 * mean_c_10.x;
    float _S911 = fy_17 * mean_c_10.y;
    float _S912 = _S909.rows[int(0)].y * _S909.rows[int(1)].x;
    float det_orig_11 = _S909.rows[int(0)].x * _S909.rows[int(1)].y - _S912;
    float _S913 = _S909.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S914 = _S909;
    *&(((&_S914)->rows + (int(0)))->x) = _S913;
    float _S915 = _S909.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S914)->rows + (int(1)))->y) = _S915;
    Matrix<float, 2, 2>  _S916 = _S914;
    Matrix<float, 2, 2>  _S917 = _S914;
    float det_blur_6 = _S913 * _S915 - _S912;
    float _S918 = det_orig_11 / det_blur_6;
    float _S919 = det_blur_6 * det_blur_6;
    float _S920 = s_primal_ctx_max_0(0.0f, _S918);
    float _S921 = s_primal_ctx_sqrt_0(_S920);
    float invdet_7 = 1.0f / det_blur_6;
    float _S922 = - _S909.rows[int(0)].y;
    float _S923 = - _S909.rows[int(1)].x;
    float _S924 = - in_opacity_10;
    float _S925 = 1.0f + s_primal_ctx_exp_1(_S924);
    float _S926 = 1.0f / _S925;
    float _S927 = _S925 * _S925;
    float _S928;
    if(antialiased_10)
    {
        _S928 = _S926 * _S921;
    }
    else
    {
        _S928 = _S926;
    }
    float _S929 = _S928 / 0.00392156885936856f;
    float _S930 = 2.0f * s_primal_ctx_log_0(_S929);
    float _S931 = s_primal_ctx_sqrt_0(_S930);
    float _S932 = _S916.rows[int(0)].x;
    float _S933 = _S917.rows[int(1)].y;
    float _S934 = s_primal_ctx_dot_0(mean_c_10, mean_c_10) + 9.99999997475242708e-07f;
    float3  _S935 = mean_11 - - s_primal_ctx_mul_0(_S888, t_13);
    float _S936 = _S935.x;
    float _S937 = _S935.y;
    float _S938 = _S935.z;
    float _S939 = _S936 * _S936 + _S937 * _S937 + _S938 * _S938;
    float _S940 = s_primal_ctx_sqrt_0(_S939);
    float x_36 = _S936 / _S940;
    float3  _S941 = make_float3 (x_36);
    float _S942 = _S940 * _S940;
    float y_13 = _S937 / _S940;
    float z_10 = _S938 / _S940;
    float3  _S943 = make_float3 (z_10);
    float _S944 = - y_13;
    float3  _S945 = make_float3 (_S944);
    float z2_25 = z_10 * z_10;
    float fTmp0B_10 = -1.09254848957061768f * z_10;
    float fC1_10 = x_36 * x_36 - y_13 * y_13;
    float _S946 = 2.0f * x_36;
    float fS1_10 = _S946 * y_13;
    float pSH6_0 = 0.94617468118667603f * z2_25 - 0.31539157032966614f;
    float3  _S947 = make_float3 (pSH6_0);
    float pSH7_0 = fTmp0B_10 * x_36;
    float3  _S948 = make_float3 (pSH7_0);
    float pSH5_0 = fTmp0B_10 * y_13;
    float3  _S949 = make_float3 (pSH5_0);
    float pSH8_0 = 0.54627424478530884f * fC1_10;
    float3  _S950 = make_float3 (pSH8_0);
    float pSH4_0 = 0.54627424478530884f * fS1_10;
    float3  _S951 = make_float3 (pSH4_0);
    float fTmp0C_10 = -2.28522896766662598f * z2_25 + 0.4570457935333252f;
    float fTmp1B_10 = 1.44530570507049561f * z_10;
    float _S952 = 1.86588168144226074f * z2_25 - 1.11952900886535645f;
    float pSH12_0 = z_10 * _S952;
    float3  _S953 = make_float3 (pSH12_0);
    float pSH13_0 = fTmp0C_10 * x_36;
    float3  _S954 = make_float3 (pSH13_0);
    float pSH11_0 = fTmp0C_10 * y_13;
    float3  _S955 = make_float3 (pSH11_0);
    float pSH14_0 = fTmp1B_10 * fC1_10;
    float3  _S956 = make_float3 (pSH14_0);
    float pSH10_0 = fTmp1B_10 * fS1_10;
    float3  _S957 = make_float3 (pSH10_0);
    float pSH15_0 = -0.59004360437393188f * (x_36 * fC1_10 - y_13 * fS1_10);
    float3  _S958 = make_float3 (pSH15_0);
    float pSH9_0 = -0.59004360437393188f * (x_36 * fS1_10 + y_13 * fC1_10);
    float3  _S959 = make_float3 (pSH9_0);
    float3  _S960 = make_float3 (0.0f);
    float3  _S961 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S962;
    (&_S962)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_10)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S944) * (*sh_coeffs_10)[int(1)] + make_float3 (z_10) * (*sh_coeffs_10)[int(2)] - make_float3 (x_36) * (*sh_coeffs_10)[int(3)]) + (make_float3 (pSH4_0) * (*sh_coeffs_10)[int(4)] + make_float3 (pSH5_0) * (*sh_coeffs_10)[int(5)] + make_float3 (pSH6_0) * (*sh_coeffs_10)[int(6)] + make_float3 (pSH7_0) * (*sh_coeffs_10)[int(7)] + make_float3 (pSH8_0) * (*sh_coeffs_10)[int(8)]) + (make_float3 (pSH9_0) * (*sh_coeffs_10)[int(9)] + make_float3 (pSH10_0) * (*sh_coeffs_10)[int(10)] + make_float3 (pSH11_0) * (*sh_coeffs_10)[int(11)] + make_float3 (pSH12_0) * (*sh_coeffs_10)[int(12)] + make_float3 (pSH13_0) * (*sh_coeffs_10)[int(13)] + make_float3 (pSH14_0) * (*sh_coeffs_10)[int(14)] + make_float3 (pSH15_0) * (*sh_coeffs_10)[int(15)]) + make_float3 (0.5f);
    (&_S962)->differential_0 = _S961;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S963;
    (&_S963)->primal_0 = _S960;
    (&_S963)->differential_0 = _S961;
    s_bwd_prop_max_0(&_S962, &_S963, v_rgb_0);
    float3  _S964 = _S958 * _S962.differential_0;
    float3  _S965 = (*sh_coeffs_10)[int(15)] * _S962.differential_0;
    float3  _S966 = _S956 * _S962.differential_0;
    float3  _S967 = (*sh_coeffs_10)[int(14)] * _S962.differential_0;
    float3  _S968 = _S954 * _S962.differential_0;
    float3  _S969 = (*sh_coeffs_10)[int(13)] * _S962.differential_0;
    float3  _S970 = _S953 * _S962.differential_0;
    float3  _S971 = (*sh_coeffs_10)[int(12)] * _S962.differential_0;
    float3  _S972 = _S955 * _S962.differential_0;
    float3  _S973 = (*sh_coeffs_10)[int(11)] * _S962.differential_0;
    float3  _S974 = _S957 * _S962.differential_0;
    float3  _S975 = (*sh_coeffs_10)[int(10)] * _S962.differential_0;
    float3  _S976 = _S959 * _S962.differential_0;
    float3  _S977 = (*sh_coeffs_10)[int(9)] * _S962.differential_0;
    float s_diff_fS2_T_0 = -0.59004360437393188f * (_S977.x + _S977.y + _S977.z);
    float s_diff_fC2_T_0 = -0.59004360437393188f * (_S965.x + _S965.y + _S965.z);
    float _S978 = _S975.x + _S975.y + _S975.z;
    float _S979 = _S967.x + _S967.y + _S967.z;
    float _S980 = _S973.x + _S973.y + _S973.z;
    float _S981 = _S969.x + _S969.y + _S969.z;
    float _S982 = _S971.x + _S971.y + _S971.z;
    float _S983 = - s_diff_fC2_T_0;
    float3  _S984 = _S950 * _S962.differential_0;
    float3  _S985 = (*sh_coeffs_10)[int(8)] * _S962.differential_0;
    float3  _S986 = _S948 * _S962.differential_0;
    float3  _S987 = (*sh_coeffs_10)[int(7)] * _S962.differential_0;
    float3  _S988 = _S947 * _S962.differential_0;
    float3  _S989 = (*sh_coeffs_10)[int(6)] * _S962.differential_0;
    float3  _S990 = _S949 * _S962.differential_0;
    float3  _S991 = (*sh_coeffs_10)[int(5)] * _S962.differential_0;
    float3  _S992 = _S951 * _S962.differential_0;
    float3  _S993 = (*sh_coeffs_10)[int(4)] * _S962.differential_0;
    float _S994 = _S991.x + _S991.y + _S991.z;
    float _S995 = _S987.x + _S987.y + _S987.z;
    float _S996 = fTmp1B_10 * _S978 + x_36 * s_diff_fS2_T_0 + y_13 * _S983 + 0.54627424478530884f * (_S993.x + _S993.y + _S993.z);
    float _S997 = fTmp1B_10 * _S979 + y_13 * s_diff_fS2_T_0 + x_36 * s_diff_fC2_T_0 + 0.54627424478530884f * (_S985.x + _S985.y + _S985.z);
    float _S998 = y_13 * - _S997;
    float _S999 = x_36 * _S997;
    float _S1000 = z_10 * (1.86588168144226074f * (z_10 * _S982) + -2.28522896766662598f * (y_13 * _S980 + x_36 * _S981) + 0.94617468118667603f * (_S989.x + _S989.y + _S989.z));
    float3  _S1001 = make_float3 (0.48860251903533936f) * _S962.differential_0;
    float3  _S1002 = - _S1001;
    float3  _S1003 = _S941 * _S1002;
    float3  _S1004 = (*sh_coeffs_10)[int(3)] * _S1002;
    float3  _S1005 = _S943 * _S1001;
    float3  _S1006 = (*sh_coeffs_10)[int(2)] * _S1001;
    float3  _S1007 = _S945 * _S1001;
    float3  _S1008 = (*sh_coeffs_10)[int(1)] * _S1001;
    float _S1009 = (_S952 * _S982 + 1.44530570507049561f * (fS1_10 * _S978 + fC1_10 * _S979) + -1.09254848957061768f * (y_13 * _S994 + x_36 * _S995) + _S1000 + _S1000 + _S1006.x + _S1006.y + _S1006.z) / _S942;
    float _S1010 = _S940 * _S1009;
    float _S1011 = (fTmp0C_10 * _S980 + fC1_10 * s_diff_fS2_T_0 + fS1_10 * _S983 + fTmp0B_10 * _S994 + _S946 * _S996 + _S998 + _S998 + - (_S1008.x + _S1008.y + _S1008.z)) / _S942;
    float _S1012 = _S940 * _S1011;
    float _S1013 = (fTmp0C_10 * _S981 + fS1_10 * s_diff_fS2_T_0 + fC1_10 * s_diff_fC2_T_0 + fTmp0B_10 * _S995 + 2.0f * (y_13 * _S996) + _S999 + _S999 + _S1004.x + _S1004.y + _S1004.z) / _S942;
    float _S1014 = _S940 * _S1013;
    float _S1015 = _S938 * - _S1009 + _S937 * - _S1011 + _S936 * - _S1013;
    DiffPair_float_0 _S1016;
    (&_S1016)->primal_0 = _S939;
    (&_S1016)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1016, _S1015);
    float _S1017 = _S938 * _S1016.differential_0;
    float _S1018 = _S937 * _S1016.differential_0;
    float _S1019 = _S936 * _S1016.differential_0;
    float3  _S1020 = make_float3 (0.282094806432724f) * _S962.differential_0;
    float3  _S1021 = make_float3 (_S1014 + _S1019 + _S1019, _S1012 + _S1018 + _S1018, _S1010 + _S1017 + _S1017);
    float3  _S1022 = - - _S1021;
    Matrix<float, 3, 3>  _S1023 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1024;
    (&_S1024)->primal_0 = _S888;
    (&_S1024)->differential_0 = _S1023;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1025;
    (&_S1025)->primal_0 = t_13;
    (&_S1025)->differential_0 = _S961;
    s_bwd_prop_mul_1(&_S1024, &_S1025, _S1022);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1026 = _S1024;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1027 = _S1025;
    float2  _S1028 = make_float2 (0.0f);
    float2  _S1029 = _S1028;
    *&((&_S1029)->y) = v_conic_0.z;
    float2  _S1030 = _S1028;
    *&((&_S1030)->y) = v_conic_0.y;
    *&((&_S1030)->x) = v_conic_0.x;
    float _S1031 = 0.5f * v_depth_0;
    DiffPair_float_0 _S1032;
    (&_S1032)->primal_0 = _S934;
    (&_S1032)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1032, _S1031);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1033;
    (&_S1033)->primal_0 = mean_c_10;
    (&_S1033)->differential_0 = _S961;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1034;
    (&_S1034)->primal_0 = mean_c_10;
    (&_S1034)->differential_0 = _S961;
    s_bwd_prop_dot_0(&_S1033, &_S1034, _S1032.differential_0);
    DiffPair_float_0 _S1035;
    (&_S1035)->primal_0 = _S933;
    (&_S1035)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1035, 0.0f);
    DiffPair_float_0 _S1036;
    (&_S1036)->primal_0 = _S932;
    (&_S1036)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1036, 0.0f);
    DiffPair_float_0 _S1037;
    (&_S1037)->primal_0 = 3.32999992370605469f;
    (&_S1037)->differential_0 = 0.0f;
    DiffPair_float_0 _S1038;
    (&_S1038)->primal_0 = _S931;
    (&_S1038)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1037, &_S1038, 0.0f);
    DiffPair_float_0 _S1039;
    (&_S1039)->primal_0 = _S930;
    (&_S1039)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1039, _S1038.differential_0);
    float _S1040 = 2.0f * _S1039.differential_0;
    DiffPair_float_0 _S1041;
    (&_S1041)->primal_0 = _S929;
    (&_S1041)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S1041, _S1040);
    float _S1042 = v_opacity_0 + 254.9999847412109375f * _S1041.differential_0;
    Matrix<float, 2, 2>  _S1043 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S1044 = _S1043;
    _S1044[int(1)] = _S1029;
    _S1044[int(0)] = _S1030;
    Matrix<float, 2, 2>  _S1045 = _S1044;
    FixedArray<float3 , 16>  _S1046;
    _S1046[int(0)] = _S961;
    _S1046[int(1)] = _S961;
    _S1046[int(2)] = _S961;
    _S1046[int(3)] = _S961;
    _S1046[int(4)] = _S961;
    _S1046[int(5)] = _S961;
    _S1046[int(6)] = _S961;
    _S1046[int(7)] = _S961;
    _S1046[int(8)] = _S961;
    _S1046[int(9)] = _S961;
    _S1046[int(10)] = _S961;
    _S1046[int(11)] = _S961;
    _S1046[int(12)] = _S961;
    _S1046[int(13)] = _S961;
    _S1046[int(14)] = _S961;
    _S1046[int(15)] = _S961;
    _S1046[int(7)] = _S986;
    _S1046[int(0)] = _S1020;
    _S1046[int(1)] = _S1007;
    _S1046[int(2)] = _S1005;
    _S1046[int(3)] = _S1003;
    _S1046[int(4)] = _S992;
    _S1046[int(5)] = _S990;
    _S1046[int(6)] = _S988;
    _S1046[int(15)] = _S964;
    _S1046[int(8)] = _S984;
    _S1046[int(9)] = _S976;
    _S1046[int(10)] = _S974;
    _S1046[int(11)] = _S972;
    _S1046[int(12)] = _S970;
    _S1046[int(13)] = _S968;
    _S1046[int(14)] = _S966;
    float3  _S1047 = _S1046[int(0)];
    float3  _S1048 = _S1046[int(1)];
    float3  _S1049 = _S1046[int(2)];
    float3  _S1050 = _S1046[int(3)];
    float3  _S1051 = _S1046[int(4)];
    float3  _S1052 = _S1046[int(5)];
    float3  _S1053 = _S1046[int(6)];
    float3  _S1054 = _S1046[int(7)];
    float3  _S1055 = _S1046[int(8)];
    float3  _S1056 = _S1046[int(9)];
    float3  _S1057 = _S1046[int(10)];
    float3  _S1058 = _S1046[int(11)];
    float3  _S1059 = _S1046[int(12)];
    float3  _S1060 = _S1046[int(13)];
    float3  _S1061 = _S1046[int(14)];
    float3  _S1062 = _S1046[int(15)];
    float3  _S1063 = _S1034.differential_0 + _S1033.differential_0;
    float2  _S1064 = make_float2 (0.0f, _S1035.differential_0);
    float2  _S1065 = make_float2 (_S1036.differential_0, 0.0f);
    float _S1066;
    if(antialiased_10)
    {
        float _S1067 = _S926 * _S1042;
        _S928 = _S921 * _S1042;
        _S1066 = _S1067;
    }
    else
    {
        _S928 = _S1042;
        _S1066 = 0.0f;
    }
    float _S1068 = - (_S928 / _S927);
    DiffPair_float_0 _S1069;
    (&_S1069)->primal_0 = _S924;
    (&_S1069)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S1069, _S1068);
    float _S1070 = - _S1069.differential_0;
    float _S1071 = invdet_7 * _S1045.rows[int(1)].y;
    float _S1072 = - (invdet_7 * _S1045.rows[int(1)].x);
    float _S1073 = - (invdet_7 * _S1045.rows[int(0)].y);
    float _S1074 = invdet_7 * _S1045.rows[int(0)].x;
    float _S1075 = - ((_S913 * _S1045.rows[int(1)].y + _S923 * _S1045.rows[int(1)].x + _S922 * _S1045.rows[int(0)].y + _S915 * _S1045.rows[int(0)].x) / _S919);
    DiffPair_float_0 _S1076;
    (&_S1076)->primal_0 = _S920;
    (&_S1076)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1076, _S1066);
    DiffPair_float_0 _S1077;
    (&_S1077)->primal_0 = 0.0f;
    (&_S1077)->differential_0 = 0.0f;
    DiffPair_float_0 _S1078;
    (&_S1078)->primal_0 = _S918;
    (&_S1078)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1077, &_S1078, _S1076.differential_0);
    float _S1079 = _S1078.differential_0 / _S919;
    float s_diff_det_orig_T_0 = det_blur_6 * _S1079;
    float _S1080 = _S1075 + det_orig_11 * - _S1079;
    float _S1081 = - _S1080;
    float _S1082 = _S913 * _S1080;
    float _S1083 = _S915 * _S1080;
    Matrix<float, 2, 2>  _S1084 = _S1043;
    _S1084[int(1)] = _S1064;
    _S1084[int(0)] = _S1065;
    _S914 = _S1084;
    *&(((&_S914)->rows + (int(1)))->y) = 0.0f;
    float _S1085 = _S1074 + _S1082 + _S1084.rows[int(1)].y;
    *&(((&_S914)->rows + (int(0)))->x) = 0.0f;
    float _S1086 = _S1071 + _S1083 + _S1084.rows[int(0)].x;
    float _S1087 = _S1081 + - s_diff_det_orig_T_0;
    float _S1088 = _S1072 + _S909.rows[int(0)].y * _S1087;
    float _S1089 = _S1073 + _S909.rows[int(1)].x * _S1087;
    float _S1090 = _S909.rows[int(1)].y * s_diff_det_orig_T_0;
    float _S1091 = _S1085 + _S909.rows[int(0)].x * s_diff_det_orig_T_0;
    float2  _S1092 = _S1028;
    *&((&_S1092)->x) = _S1088;
    *&((&_S1092)->y) = _S1091;
    float _S1093 = _S1086 + _S1090;
    float2  _S1094 = _S1028;
    *&((&_S1094)->y) = _S1089;
    *&((&_S1094)->x) = _S1093;
    float _S1095 = _S911 * v_mean2d_0.y;
    float _S1096 = fy_17 * (rz_3 * v_mean2d_0.y);
    float _S1097 = _S910 * v_mean2d_0.x;
    float _S1098 = fx_17 * (rz_3 * v_mean2d_0.x);
    Matrix<float, 2, 2>  _S1099 = _S1043;
    _S1099[int(1)] = _S1092;
    _S1099[int(0)] = _S1094;
    Matrix<float, 2, 2>  _S1100 = _S914 + _S1099;
    Matrix<float, 2, 3>  _S1101 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1102;
    (&_S1102)->primal_0 = _S907;
    (&_S1102)->differential_0 = _S1101;
    Matrix<float, 3, 2>  _S1103 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1104;
    (&_S1104)->primal_0 = _S908;
    (&_S1104)->differential_0 = _S1103;
    s_bwd_prop_mul_2(&_S1102, &_S1104, _S1100);
    Matrix<float, 2, 3>  _S1105 = transpose_2(_S1104.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1106;
    (&_S1106)->primal_0 = J_8;
    (&_S1106)->differential_0 = _S1101;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1107;
    (&_S1107)->primal_0 = _S889;
    (&_S1107)->differential_0 = _S1023;
    s_bwd_prop_mul_3(&_S1106, &_S1107, _S1102.differential_0);
    Matrix<float, 2, 3>  _S1108 = _S1105 + _S1106.differential_0;
    float _S1109 = _S906 * _S1108.rows[int(1)].z;
    float s_diff_ty_T_0 = _S905 * (rz2_3 * _S1108.rows[int(1)].z);
    float _S1110 = fy_17 * _S1108.rows[int(1)].y;
    float _S1111 = _S904 * _S1108.rows[int(0)].z;
    float s_diff_tx_T_0 = _S903 * (rz2_3 * _S1108.rows[int(0)].z);
    float _S1112 = fx_17 * _S1108.rows[int(0)].x;
    float _S1113 = mean_c_10.z * s_diff_ty_T_0;
    float _S1114 = _S902 * s_diff_ty_T_0;
    DiffPair_float_0 _S1115;
    (&_S1115)->primal_0 = lim_y_pos_2;
    (&_S1115)->differential_0 = 0.0f;
    DiffPair_float_0 _S1116;
    (&_S1116)->primal_0 = _S901;
    (&_S1116)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1115, &_S1116, _S1113);
    DiffPair_float_0 _S1117;
    (&_S1117)->primal_0 = _S899;
    (&_S1117)->differential_0 = 0.0f;
    DiffPair_float_0 _S1118;
    (&_S1118)->primal_0 = _S900;
    (&_S1118)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1117, &_S1118, _S1116.differential_0);
    float _S1119 = mean_c_10.y * _S1118.differential_0;
    float _S1120 = rz_3 * _S1118.differential_0;
    float _S1121 = mean_c_10.z * s_diff_tx_T_0;
    float _S1122 = _S898 * s_diff_tx_T_0;
    DiffPair_float_0 _S1123;
    (&_S1123)->primal_0 = lim_x_pos_2;
    (&_S1123)->differential_0 = 0.0f;
    DiffPair_float_0 _S1124;
    (&_S1124)->primal_0 = _S897;
    (&_S1124)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S1123, &_S1124, _S1121);
    DiffPair_float_0 _S1125;
    (&_S1125)->primal_0 = _S895;
    (&_S1125)->differential_0 = 0.0f;
    DiffPair_float_0 _S1126;
    (&_S1126)->primal_0 = _S896;
    (&_S1126)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1125, &_S1126, _S1124.differential_0);
    float _S1127 = rz_3 * (_S1109 + _S1111);
    float _S1128 = _S1114 + _S1122 + - ((_S1095 + _S1097 + _S1110 + _S1112 + _S1119 + mean_c_10.x * _S1126.differential_0 + _S1127 + _S1127) / _S894);
    float _S1129 = _S1096 + _S1120;
    float _S1130 = _S1098 + rz_3 * _S1126.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1131;
    (&_S1131)->primal_0 = _S887;
    (&_S1131)->differential_0 = _S1023;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1132;
    (&_S1132)->primal_0 = _S888;
    (&_S1132)->differential_0 = _S1023;
    s_bwd_prop_mul_4(&_S1131, &_S1132, _S1107.differential_0);
    Matrix<float, 3, 3>  _S1133 = transpose_0(_S1132.differential_0 + _S1026.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1134;
    (&_S1134)->primal_0 = R_16;
    (&_S1134)->differential_0 = _S1023;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1135;
    (&_S1135)->primal_0 = _S886;
    (&_S1135)->differential_0 = _S1023;
    s_bwd_prop_mul_4(&_S1134, &_S1135, _S1131.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1136;
    (&_S1136)->primal_0 = _S884;
    (&_S1136)->differential_0 = _S1023;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1137;
    (&_S1137)->primal_0 = _S885;
    (&_S1137)->differential_0 = _S1023;
    s_bwd_prop_mul_4(&_S1136, &_S1137, _S1135.differential_0);
    Matrix<float, 3, 3>  _S1138 = _S1136.differential_0 + transpose_0(_S1137.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1139;
    (&_S1139)->primal_0 = _S883;
    (&_S1139)->differential_0 = _S1023;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1140;
    (&_S1140)->primal_0 = S_0;
    (&_S1140)->differential_0 = _S1023;
    s_bwd_prop_mul_4(&_S1139, &_S1140, _S1138);
    Matrix<float, 3, 3>  _S1141 = transpose_0(_S1139.differential_0);
    float _S1142 = 2.0f * - _S1141.rows[int(2)].z;
    float _S1143 = 2.0f * _S1141.rows[int(2)].y;
    float _S1144 = 2.0f * _S1141.rows[int(2)].x;
    float _S1145 = 2.0f * _S1141.rows[int(1)].z;
    float _S1146 = 2.0f * - _S1141.rows[int(1)].y;
    float _S1147 = 2.0f * _S1141.rows[int(1)].x;
    float _S1148 = 2.0f * _S1141.rows[int(0)].z;
    float _S1149 = 2.0f * _S1141.rows[int(0)].y;
    float _S1150 = 2.0f * - _S1141.rows[int(0)].x;
    float _S1151 = - _S1147 + _S1149;
    float _S1152 = _S1144 + - _S1148;
    float _S1153 = - _S1143 + _S1145;
    float _S1154 = _S1143 + _S1145;
    float _S1155 = _S1144 + _S1148;
    float _S1156 = _S1147 + _S1149;
    float _S1157 = quat_14.w * (_S1146 + _S1150);
    float _S1158 = quat_14.z * (_S1142 + _S1150);
    float _S1159 = quat_14.y * (_S1142 + _S1146);
    float _S1160 = quat_14.x * _S1151 + quat_14.z * _S1154 + quat_14.y * _S1155 + _S1157 + _S1157;
    float _S1161 = quat_14.x * _S1152 + quat_14.w * _S1154 + quat_14.y * _S1156 + _S1158 + _S1158;
    float _S1162 = quat_14.x * _S1153 + quat_14.w * _S1155 + quat_14.z * _S1156 + _S1159 + _S1159;
    float _S1163 = quat_14.w * _S1151 + quat_14.z * _S1152 + quat_14.y * _S1153;
    float3  _S1164 = _S961;
    *&((&_S1164)->z) = _S1140.differential_0.rows[int(2)].z;
    *&((&_S1164)->y) = _S1140.differential_0.rows[int(1)].y;
    *&((&_S1164)->x) = _S1140.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1165;
    (&_S1165)->primal_0 = scale_13;
    (&_S1165)->differential_0 = _S961;
    s_bwd_prop_exp_1(&_S1165, _S1164);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1166 = _S1165;
    float3  _S1167 = _S961;
    *&((&_S1167)->z) = _S1128;
    *&((&_S1167)->y) = _S1129;
    *&((&_S1167)->x) = _S1130;
    float3  _S1168 = _S1063 + _S1167;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1169;
    (&_S1169)->primal_0 = R_16;
    (&_S1169)->differential_0 = _S1023;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1170;
    (&_S1170)->primal_0 = mean_11;
    (&_S1170)->differential_0 = _S961;
    s_bwd_prop_mul_1(&_S1169, &_S1170, _S1168);
    float3  _S1171 = _S1168 + _S1027.differential_0;
    Matrix<float, 3, 3>  _S1172 = _S1133 + _S1134.differential_0 + _S1169.differential_0;
    float4  _S1173 = make_float4 (0.0f);
    *&((&_S1173)->w) = _S1160;
    *&((&_S1173)->z) = _S1161;
    *&((&_S1173)->y) = _S1162;
    *&((&_S1173)->x) = _S1163;
    float4  _S1174 = _S1173;
    float3  _S1175 = _S1170.differential_0 + _S1021;
    *v_mean_0 = _S1175;
    *v_quat_0 = _S1174;
    *v_scale_0 = _S1166.differential_0;
    *v_in_opacity_0 = _S1070;
    (*v_sh_coeffs_0)[int(0)] = _S1047;
    (*v_sh_coeffs_0)[int(1)] = _S1048;
    (*v_sh_coeffs_0)[int(2)] = _S1049;
    (*v_sh_coeffs_0)[int(3)] = _S1050;
    (*v_sh_coeffs_0)[int(4)] = _S1051;
    (*v_sh_coeffs_0)[int(5)] = _S1052;
    (*v_sh_coeffs_0)[int(6)] = _S1053;
    (*v_sh_coeffs_0)[int(7)] = _S1054;
    (*v_sh_coeffs_0)[int(8)] = _S1055;
    (*v_sh_coeffs_0)[int(9)] = _S1056;
    (*v_sh_coeffs_0)[int(10)] = _S1057;
    (*v_sh_coeffs_0)[int(11)] = _S1058;
    (*v_sh_coeffs_0)[int(12)] = _S1059;
    (*v_sh_coeffs_0)[int(13)] = _S1060;
    (*v_sh_coeffs_0)[int(14)] = _S1061;
    (*v_sh_coeffs_0)[int(15)] = _S1062;
    *v_R_2 = _S1172;
    *v_t_1 = _S1171;
    return;
}

struct s_bwd_prop_s_bwd_prop_atan2_Intermediates_0
{
    DiffPair_float_0 _S1176;
    DiffPair_float_0 _S1177;
};

struct s_bwd_prop_fisheye_proj_3dgs_Intermediates_0
{
    DiffPair_float_0 _S1178;
    DiffPair_float_0 _S1179;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1180;
    DiffPair_float_0 _S1181;
    DiffPair_float_0 _S1182;
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1183;
};

struct s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S1184;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1185;
};

inline __device__ float s_primal_ctx_s_primal_ctx_atan2_0(float _S1186, float _S1187)
{
    return s_primal_ctx_atan2_0(_S1186, _S1187);
}

struct s_bwd_prop_d_atan2_Intermediates_0
{
    DiffPair_float_0 _S1188;
    DiffPair_float_0 _S1189;
};

inline __device__ void s_primal_ctx_d_atan2_0(DiffPair_float_0 * dpdpy_0, DiffPair_float_0 * dpdpx_0, float dpdOut_0, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_0)
{
    DiffPair_float_0 _S1190 = { 0.0f, 0.0f };
    _s_diff_ctx_0->_S1188 = _S1190;
    _s_diff_ctx_0->_S1189 = _S1190;
    (&_s_diff_ctx_0->_S1188)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S1188)->differential_0 = 0.0f;
    (&_s_diff_ctx_0->_S1189)->primal_0 = 0.0f;
    (&_s_diff_ctx_0->_S1189)->differential_0 = 0.0f;
    DiffPair_float_0 _S1191 = *dpdpy_0;
    _s_diff_ctx_0->_S1188 = *dpdpy_0;
    DiffPair_float_0 _S1192 = *dpdpx_0;
    _s_diff_ctx_0->_S1189 = *dpdpx_0;
    float _S1193 = _S1192.primal_0 * _S1192.primal_0 + _S1191.primal_0 * _S1191.primal_0;
    float _S1194 = - _S1191.primal_0 / _S1193 * dpdOut_0;
    float _S1195 = _S1192.primal_0 / _S1193 * dpdOut_0;
    dpdpy_0->primal_0 = _S1191.primal_0;
    dpdpy_0->differential_0 = _S1195;
    dpdpx_0->primal_0 = _S1192.primal_0;
    dpdpx_0->differential_0 = _S1194;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_atan2_0(DiffPair_float_0 * _S1196, DiffPair_float_0 * _S1197, float _S1198, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_1)
{
    DiffPair_float_0 _S1199 = { 0.0f, 0.0f };
    _s_diff_ctx_1->_S1176 = _S1199;
    _s_diff_ctx_1->_S1177 = _S1199;
    (&_s_diff_ctx_1->_S1176)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S1176)->differential_0 = 0.0f;
    (&_s_diff_ctx_1->_S1177)->primal_0 = 0.0f;
    (&_s_diff_ctx_1->_S1177)->differential_0 = 0.0f;
    DiffPair_float_0 _S1200 = *_S1196;
    _s_diff_ctx_1->_S1176 = *_S1196;
    DiffPair_float_0 _S1201 = *_S1197;
    _s_diff_ctx_1->_S1177 = *_S1197;
    DiffPair_float_0 _S1202 = _S1200;
    DiffPair_float_0 _S1203 = _S1201;
    s_bwd_prop_d_atan2_Intermediates_0 _S1204;
    (&_S1204)->_S1188 = _S1199;
    (&_S1204)->_S1189 = _S1199;
    s_primal_ctx_d_atan2_0(&_S1202, &_S1203, _S1198, &_S1204);
    *_S1196 = _S1202;
    *_S1197 = _S1203;
    return;
}

struct s_bwd_prop_s_bwd_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1205;
};

struct s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1206;
};

struct s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1207;
    DiffPair_float_0 _S1208;
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1209;
};

struct s_bwd_prop_d_sqrt_Intermediates_0
{
    DiffPair_float_0 _S1210;
};

inline __device__ void s_primal_ctx_d_sqrt_0(DiffPair_float_0 * dpdpx_1, float dpdOut_1, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_2)
{
    DiffPair_float_0 _S1211 = { 0.0f, 0.0f };
    _s_diff_ctx_2->_S1210 = _S1211;
    (&_s_diff_ctx_2->_S1210)->primal_0 = 0.0f;
    (&_s_diff_ctx_2->_S1210)->differential_0 = 0.0f;
    DiffPair_float_0 _S1212 = *dpdpx_1;
    _s_diff_ctx_2->_S1210 = *dpdpx_1;
    float _S1213 = 0.5f / s_primal_ctx_sqrt_0(s_primal_ctx_max_0(1.00000001168609742e-07f, _S1212.primal_0)) * dpdOut_1;
    dpdpx_1->primal_0 = _S1212.primal_0;
    dpdpx_1->differential_0 = _S1213;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_sqrt_0(DiffPair_float_0 * _S1214, float _S1215, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_3)
{
    DiffPair_float_0 _S1216 = { 0.0f, 0.0f };
    _s_diff_ctx_3->_S1206 = _S1216;
    (&_s_diff_ctx_3->_S1206)->primal_0 = 0.0f;
    (&_s_diff_ctx_3->_S1206)->differential_0 = 0.0f;
    DiffPair_float_0 _S1217 = *_S1214;
    _s_diff_ctx_3->_S1206 = *_S1214;
    DiffPair_float_0 _S1218 = _S1217;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1219;
    (&_S1219)->_S1210 = _S1216;
    s_primal_ctx_d_sqrt_0(&_S1218, _S1215, &_S1219);
    *_S1214 = _S1218;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_prop_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_2, float dp_s_dOut_0, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_4)
{
    float2  _S1220 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1221 = { _S1220, _S1220 };
    DiffPair_float_0 _S1222 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1223 = { _S1222 };
    _s_diff_ctx_4->_S1207 = _S1221;
    _s_diff_ctx_4->_S1208 = _S1222;
    _s_diff_ctx_4->_S1209 = _S1223;
    (&_s_diff_ctx_4->_S1207)->primal_0 = _S1220;
    (&_s_diff_ctx_4->_S1207)->differential_0 = _S1220;
    (&_s_diff_ctx_4->_S1208)->primal_0 = 0.0f;
    (&_s_diff_ctx_4->_S1208)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1224 = *dpdpx_2;
    _s_diff_ctx_4->_S1207 = *dpdpx_2;
    float _S1225 = _S1224.primal_0.x;
    float _S1226 = _S1224.primal_0.y;
    DiffPair_float_0 _S1227;
    (&_S1227)->primal_0 = _S1225 * _S1225 + _S1226 * _S1226;
    (&_S1227)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_sqrt_0(&_S1227, dp_s_dOut_0, &_s_diff_ctx_4->_S1209);
    _s_diff_ctx_4->_S1208 = _S1227;
    float _S1228 = _S1224.primal_0.y * _S1227.differential_0;
    float _S1229 = _S1228 + _S1228;
    float _S1230 = _S1224.primal_0.x * _S1227.differential_0;
    float _S1231 = _S1230 + _S1230;
    float2  _S1232 = _S1220;
    *&((&_S1232)->y) = _S1229;
    *&((&_S1232)->x) = _S1231;
    dpdpx_2->primal_0 = _S1224.primal_0;
    dpdpx_2->differential_0 = _S1232;
    return;
}

inline __device__ void s_primal_ctx_s_bwd_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S1233, float _S1234, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_5)
{
    float2  _S1235 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1236 = { _S1235, _S1235 };
    _s_diff_ctx_5->_S1205 = _S1236;
    (&_s_diff_ctx_5->_S1205)->primal_0 = _S1235;
    (&_s_diff_ctx_5->_S1205)->differential_0 = _S1235;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1237 = *_S1233;
    _s_diff_ctx_5->_S1205 = *_S1233;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1238 = _S1237;
    DiffPair_float_0 _S1239 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1240 = { _S1239 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1241;
    (&_S1241)->_S1207 = _S1236;
    (&_S1241)->_S1208 = _S1239;
    (&_S1241)->_S1209 = _S1240;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1238, _S1234, &_S1241);
    *_S1233 = _S1238;
    return;
}

inline __device__ bool s_primal_ctx_fisheye_proj_3dgs_0(float3  dpmean3d_0, Matrix<float, 3, 3>  dpcov3d_0, float4  dpintrins_0, FixedArray<float, 10>  * dpdist_coeffs_0, Matrix<float, 2, 2>  * dpcov2d_0, float2  * dpmean2d_0, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_6)
{
    DiffPair_float_0 _S1242 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1243 = { _S1242, _S1242 };
    _s_diff_ctx_6->_S1178 = _S1242;
    _s_diff_ctx_6->_S1179 = _S1242;
    _s_diff_ctx_6->_S1180 = _S1243;
    _s_diff_ctx_6->_S1181 = _S1242;
    _s_diff_ctx_6->_S1182 = _S1242;
    _s_diff_ctx_6->_S1183 = _S1243;
    (&_s_diff_ctx_6->_S1178)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1178)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S1179)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1179)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S1181)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1181)->differential_0 = 0.0f;
    (&_s_diff_ctx_6->_S1182)->primal_0 = 0.0f;
    (&_s_diff_ctx_6->_S1182)->differential_0 = 0.0f;
    float2  _S1244 = make_float2 (0.0f);
    float2  _S1245 = float2 {dpmean3d_0.x, dpmean3d_0.y};
    float _S1246 = length_0(_S1245);
    float _S1247 = dpmean3d_0.z;
    float _S1248 = s_primal_ctx_atan2_0(_S1246, _S1247);
    float k_4;
    if(_S1248 < 0.00100000004749745f)
    {
        k_4 = (1.0f - _S1248 * _S1248 / 3.0f) / _S1247;
    }
    else
    {
        k_4 = _S1248 / _S1246;
    }
    float2  _S1249 = _S1245 * make_float2 (k_4);
    float u_38 = _S1249.x;
    float v_38 = _S1249.y;
    float r2_38 = u_38 * u_38 + v_38 * v_38;
    float _S1250 = 2.0f * (*dpdist_coeffs_0)[int(4)];
    float _S1251 = 2.0f * (*dpdist_coeffs_0)[int(5)];
    float2  _S1252 = _S1249 * make_float2 (1.0f + r2_38 * ((*dpdist_coeffs_0)[int(0)] + r2_38 * ((*dpdist_coeffs_0)[int(1)] + r2_38 * ((*dpdist_coeffs_0)[int(2)] + r2_38 * (*dpdist_coeffs_0)[int(3)])))) + make_float2 (_S1250 * u_38 * v_38 + (*dpdist_coeffs_0)[int(5)] * (r2_38 + 2.0f * u_38 * u_38) + (*dpdist_coeffs_0)[int(6)] * r2_38, _S1251 * u_38 * v_38 + (*dpdist_coeffs_0)[int(4)] * (r2_38 + 2.0f * v_38 * v_38) + (*dpdist_coeffs_0)[int(7)] * r2_38);
    float2  _S1253 = _S1252 + make_float2 ((*dpdist_coeffs_0)[int(8)] * _S1252.x + (*dpdist_coeffs_0)[int(9)] * _S1252.y, 0.0f);
    float fx_18 = dpintrins_0.x;
    float fy_18 = dpintrins_0.y;
    float2  _S1254 = make_float2 (fx_18 * _S1253.x + dpintrins_0.z, fy_18 * _S1253.y + dpintrins_0.w);
    Matrix<float, 2, 3>  J_9 = makeMatrix<float, 2, 3> (0.0f);
    float _S1255 = s_primal_ctx_s_primal_ctx_atan2_0(_S1246, _S1247);
    bool _S1256 = _S1255 < 0.00100000004749745f;
    float _S1257;
    float _S1258;
    float _S1259;
    if(_S1256)
    {
        float _S1260 = 1.0f - _S1255 * _S1255 / 3.0f;
        float _S1261 = _S1247 * _S1247;
        k_4 = _S1260 / _S1247;
        _S1257 = 0.0f;
        _S1258 = _S1261;
        _S1259 = _S1260;
    }
    else
    {
        float _S1262 = _S1246 * _S1246;
        k_4 = _S1255 / _S1246;
        _S1257 = _S1262;
        _S1258 = 0.0f;
        _S1259 = 0.0f;
    }
    float2  _S1263 = make_float2 (k_4);
    float2  _S1264 = _S1245 * make_float2 (k_4);
    float u_39 = _S1264.x;
    float v_39 = _S1264.y;
    float r2_39 = u_39 * u_39 + v_39 * v_39;
    float _S1265 = (*dpdist_coeffs_0)[int(2)] + r2_39 * (*dpdist_coeffs_0)[int(3)];
    float _S1266 = (*dpdist_coeffs_0)[int(1)] + r2_39 * _S1265;
    float _S1267 = (*dpdist_coeffs_0)[int(0)] + r2_39 * _S1266;
    float2  _S1268 = make_float2 (fx_18, 0.0f) + make_float2 ((*dpdist_coeffs_0)[int(8)] * fx_18, (*dpdist_coeffs_0)[int(9)] * fx_18);
    float2  _S1269 = _S1264 * _S1268;
    float _S1270 = (*dpdist_coeffs_0)[int(4)] * _S1268.y;
    float _S1271 = (*dpdist_coeffs_0)[int(5)] * _S1268.x;
    float _S1272 = _S1269.x + _S1269.y;
    float _S1273 = r2_39 * _S1272;
    float _S1274 = r2_39 * _S1273;
    float _S1275 = (*dpdist_coeffs_0)[int(7)] * _S1268.y + _S1270 + (*dpdist_coeffs_0)[int(6)] * _S1268.x + _S1271 + _S1267 * _S1272 + _S1266 * _S1273 + _S1265 * _S1274 + (*dpdist_coeffs_0)[int(3)] * (r2_39 * _S1274);
    float _S1276 = v_39 * _S1275;
    float _S1277 = u_39 * _S1275;
    float2  _S1278 = make_float2 (1.0f + r2_39 * _S1267) * _S1268 + make_float2 (_S1251 * (v_39 * _S1268.y) + 2.0f * u_39 * _S1271 + 2.0f * (u_39 * _S1271) + _S1250 * (v_39 * _S1268.x) + _S1277 + _S1277, 2.0f * v_39 * _S1270 + 2.0f * (v_39 * _S1270) + _S1251 * u_39 * _S1268.y + _S1250 * u_39 * _S1268.x + _S1276 + _S1276);
    float2  _S1279 = _S1245 * _S1278;
    float2  _S1280 = _S1263 * _S1278;
    float _S1281 = _S1279.x + _S1279.y;
    if(_S1256)
    {
        float _S1282 = _S1281 / _S1258;
        float _S1283 = _S1259 * - _S1282;
        float _S1284 = _S1255 * (0.3333333432674408f * - (_S1247 * _S1282));
        k_4 = _S1284 + _S1284;
        _S1257 = _S1283;
        _S1258 = 0.0f;
    }
    else
    {
        float _S1285 = _S1281 / _S1257;
        float _S1286 = _S1255 * - _S1285;
        k_4 = _S1246 * _S1285;
        _S1257 = 0.0f;
        _S1258 = _S1286;
    }
    DiffPair_float_0 _S1287;
    (&_S1287)->primal_0 = _S1246;
    (&_S1287)->differential_0 = 0.0f;
    DiffPair_float_0 _S1288;
    (&_S1288)->primal_0 = _S1247;
    (&_S1288)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1287, &_S1288, k_4, &_s_diff_ctx_6->_S1180);
    _s_diff_ctx_6->_S1178 = _S1287;
    _s_diff_ctx_6->_S1179 = _S1288;
    float _S1289 = _S1288.differential_0 + _S1257;
    float _S1290 = _S1287.differential_0 + _S1258;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1291;
    (&_S1291)->primal_0 = _S1245;
    (&_S1291)->differential_0 = _S1244;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1292 = { _S1244, _S1244 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1293;
    (&_S1293)->_S1205 = _S1292;
    s_primal_ctx_s_bwd_length_impl_0(&_S1291, _S1290, &_S1293);
    float2  _S1294 = _S1291.differential_0 + _S1280;
    float3  _S1295 = make_float3 (_S1294.x, _S1294.y, _S1289);
    Matrix<float, 2, 3>  _S1296 = J_9;
    _S1296[int(0)] = _S1295;
    if(_S1256)
    {
        float _S1297 = 1.0f - _S1255 * _S1255 / 3.0f;
        float _S1298 = _S1247 * _S1247;
        k_4 = _S1297 / _S1247;
        _S1257 = 0.0f;
        _S1258 = _S1298;
        _S1259 = _S1297;
    }
    else
    {
        float _S1299 = _S1246 * _S1246;
        k_4 = _S1255 / _S1246;
        _S1257 = _S1299;
        _S1258 = 0.0f;
        _S1259 = 0.0f;
    }
    float2  _S1300 = make_float2 (k_4);
    float2  _S1301 = _S1245 * make_float2 (k_4);
    float u_40 = _S1301.x;
    float v_40 = _S1301.y;
    float r2_40 = u_40 * u_40 + v_40 * v_40;
    float _S1302 = (*dpdist_coeffs_0)[int(2)] + r2_40 * (*dpdist_coeffs_0)[int(3)];
    float _S1303 = (*dpdist_coeffs_0)[int(1)] + r2_40 * _S1302;
    float _S1304 = (*dpdist_coeffs_0)[int(0)] + r2_40 * _S1303;
    float2  _S1305 = make_float2 (0.0f, fy_18);
    float2  _S1306 = _S1301 * _S1305;
    float _S1307 = (*dpdist_coeffs_0)[int(4)] * fy_18;
    float _S1308 = _S1306.x + _S1306.y;
    float _S1309 = r2_40 * _S1308;
    float _S1310 = r2_40 * _S1309;
    float _S1311 = (*dpdist_coeffs_0)[int(7)] * fy_18 + _S1307 + _S1304 * _S1308 + _S1303 * _S1309 + _S1302 * _S1310 + (*dpdist_coeffs_0)[int(3)] * (r2_40 * _S1310);
    float _S1312 = v_40 * _S1311;
    float _S1313 = u_40 * _S1311;
    float2  _S1314 = make_float2 (1.0f + r2_40 * _S1304) * _S1305 + make_float2 (_S1251 * (v_40 * fy_18) + _S1313 + _S1313, 2.0f * v_40 * _S1307 + 2.0f * (v_40 * _S1307) + _S1251 * u_40 * fy_18 + _S1312 + _S1312);
    float2  _S1315 = _S1245 * _S1314;
    float2  _S1316 = _S1300 * _S1314;
    float _S1317 = _S1315.x + _S1315.y;
    if(_S1256)
    {
        float _S1318 = _S1317 / _S1258;
        float _S1319 = _S1259 * - _S1318;
        float _S1320 = _S1255 * (0.3333333432674408f * - (_S1247 * _S1318));
        k_4 = _S1320 + _S1320;
        _S1257 = _S1319;
        _S1258 = 0.0f;
    }
    else
    {
        float _S1321 = _S1317 / _S1257;
        float _S1322 = _S1255 * - _S1321;
        k_4 = _S1246 * _S1321;
        _S1257 = 0.0f;
        _S1258 = _S1322;
    }
    DiffPair_float_0 _S1323;
    (&_S1323)->primal_0 = _S1246;
    (&_S1323)->differential_0 = 0.0f;
    DiffPair_float_0 _S1324;
    (&_S1324)->primal_0 = _S1247;
    (&_S1324)->differential_0 = 0.0f;
    s_primal_ctx_s_bwd_prop_atan2_0(&_S1323, &_S1324, k_4, &_s_diff_ctx_6->_S1183);
    _s_diff_ctx_6->_S1181 = _S1323;
    _s_diff_ctx_6->_S1182 = _S1324;
    float _S1325 = _S1324.differential_0 + _S1257;
    float _S1326 = _S1323.differential_0 + _S1258;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1327;
    (&_S1327)->primal_0 = _S1245;
    (&_S1327)->differential_0 = _S1244;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1328;
    (&_S1328)->_S1205 = _S1292;
    s_primal_ctx_s_bwd_length_impl_0(&_S1327, _S1326, &_S1328);
    float2  _S1329 = _S1327.differential_0 + _S1316;
    float3  _S1330 = make_float3 (_S1329.x, _S1329.y, _S1325);
    _S1296[int(1)] = _S1330;
    *dpcov2d_0 = s_primal_ctx_mul_3(s_primal_ctx_mul_2(_S1296, dpcov3d_0), transpose_1(_S1296));
    *dpmean2d_0 = _S1254;
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

inline __device__ void s_bwd_prop_d_sqrt_0(DiffPair_1 * dpdpx_3, DiffPair_float_0 * dpdOut_2, s_bwd_prop_d_sqrt_Intermediates_0 * _s_diff_ctx_7)
{
    DiffPair_1 _S1331 = *dpdpx_3;
    float _S1332 = s_primal_ctx_max_0(1.00000001168609742e-07f, (&_s_diff_ctx_7->_S1210)->primal_0);
    float _S1333 = s_primal_ctx_sqrt_0(_S1332);
    float _S1334 = 0.5f / _S1333 * (*dpdpx_3).differential_0.differential_0;
    float _S1335 = 0.5f * - ((*dpdOut_2).primal_0 * (*dpdpx_3).differential_0.differential_0 / (_S1333 * _S1333));
    DiffPair_float_0 _S1336;
    (&_S1336)->primal_0 = _S1332;
    (&_S1336)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1336, _S1335);
    DiffPair_float_0 _S1337;
    (&_S1337)->primal_0 = 1.00000001168609742e-07f;
    (&_S1337)->differential_0 = 0.0f;
    DiffPair_float_0 _S1338;
    (&_S1338)->primal_0 = (&_s_diff_ctx_7->_S1210)->primal_0;
    (&_S1338)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S1337, &_S1338, _S1336.differential_0);
    DiffPair_float_0 dpdpx_4 = { _S1338.differential_0 + (*dpdpx_3).differential_0.primal_0, 0.0f };
    dpdOut_2->primal_0 = (*dpdOut_2).primal_0;
    dpdOut_2->differential_0 = _S1334;
    dpdpx_3->primal_0 = _S1331.primal_0;
    dpdpx_3->differential_0 = dpdpx_4;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_sqrt_0(DiffPair_1 * _S1339, DiffPair_float_0 * _S1340, s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 * _s_diff_ctx_8)
{
    DiffPair_1 _S1341 = *_S1339;
    DiffPair_float_0 _S1342 = _s_diff_ctx_8->_S1206;
    DiffPair_float_0 _S1343 = { 0.0f, 0.0f };
    s_bwd_prop_d_sqrt_Intermediates_0 _S1344;
    (&_S1344)->_S1210 = _S1343;
    s_primal_ctx_d_sqrt_0(&_S1342, (*_S1340).primal_0, &_S1344);
    DiffPair_float_0 _S1345 = { (*_S1339).differential_0.primal_0, (*_S1339).differential_0.differential_0 };
    DiffPair_1 _S1346;
    (&_S1346)->primal_0 = _s_diff_ctx_8->_S1206;
    (&_S1346)->differential_0 = _S1345;
    DiffPair_float_0 _S1347;
    (&_S1347)->primal_0 = (*_S1340).primal_0;
    (&_S1347)->differential_0 = 0.0f;
    s_bwd_prop_d_sqrt_Intermediates_0 _S1348 = _S1344;
    s_bwd_prop_d_sqrt_0(&_S1346, &_S1347, &_S1348);
    DiffPair_float_0 _S1349 = { _S1346.differential_0.primal_0, _S1346.differential_0.differential_0 };
    _S1340->primal_0 = (*_S1340).primal_0;
    _S1340->differential_0 = _S1347.differential_0;
    _S1339->primal_0 = _S1341.primal_0;
    _S1339->differential_0 = _S1349;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_sqrt_0(DiffPair_float_0 * _S1350, float _s_dOut_3)
{
    DiffPair_float_0 _S1351;
    (&_S1351)->primal_0 = (*_S1350).primal_0;
    (&_S1351)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1351, _s_dOut_3);
    _S1350->primal_0 = (*_S1350).primal_0;
    _S1350->differential_0 = _S1351.differential_0;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_length_impl_0(DiffPair_0 * dpdpx_5, DiffPair_float_0 * dp_s_dOut_1, s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 * _s_diff_ctx_9)
{
    DiffPair_0 _S1352 = *dpdpx_5;
    float len_0 = *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->x) * *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->x) + *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->y) * *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->y);
    DiffPair_float_0 _S1353 = { len_0, 0.0f };
    float2  _S1354 = make_float2 (0.0f);
    float _S1355 = (*dpdpx_5).differential_0.differential_0.x;
    float _S1356 = _S1355 + _S1355;
    float _S1357 = (&_s_diff_ctx_9->_S1208)->differential_0 * _S1356;
    float _S1358 = (*dpdpx_5).differential_0.differential_0.y + (*dpdpx_5).differential_0.differential_0.y;
    float _S1359 = (&_s_diff_ctx_9->_S1208)->differential_0 * _S1358;
    DiffPair_float_0 _S1360 = { 0.0f, *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->x) * _S1356 + *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->y) * _S1358 };
    DiffPair_1 _S1361;
    (&_S1361)->primal_0 = _S1353;
    (&_S1361)->differential_0 = _S1360;
    DiffPair_float_0 _S1362;
    (&_S1362)->primal_0 = (*dp_s_dOut_1).primal_0;
    (&_S1362)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_sqrt_0(&_S1361, &_S1362, &_s_diff_ctx_9->_S1209);
    DiffPair_float_0 _S1363;
    (&_S1363)->primal_0 = len_0;
    (&_S1363)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1363, 0.0f);
    float _S1364 = _S1361.differential_0.primal_0 + _S1363.differential_0;
    float _S1365 = *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->y) * _S1364;
    float _S1366 = _S1359 + _S1365 + _S1365;
    float _S1367 = *&((&(&_s_diff_ctx_9->_S1207)->primal_0)->x) * _S1364;
    float _S1368 = _S1357 + _S1367 + _S1367;
    float2  _S1369 = _S1354;
    *&((&_S1369)->y) = _S1366;
    *&((&_S1369)->x) = _S1368;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dpdpx_6 = { _S1352.differential_0.primal_0 + _S1369, _S1354 };
    dp_s_dOut_1->primal_0 = (*dp_s_dOut_1).primal_0;
    dp_s_dOut_1->differential_0 = _S1362.differential_0;
    dpdpx_5->primal_0 = _S1352.primal_0;
    dpdpx_5->differential_0 = dpdpx_6;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_length_impl_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpdpx_7, float _s_dOut_4)
{
    float _S1370 = (*dpdpx_7).primal_0.x;
    float _S1371 = (*dpdpx_7).primal_0.y;
    DiffPair_float_0 _S1372;
    (&_S1372)->primal_0 = _S1370 * _S1370 + _S1371 * _S1371;
    (&_S1372)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_sqrt_0(&_S1372, _s_dOut_4);
    float _S1373 = (*dpdpx_7).primal_0.y * _S1372.differential_0;
    float _S1374 = _S1373 + _S1373;
    float _S1375 = (*dpdpx_7).primal_0.x * _S1372.differential_0;
    float _S1376 = _S1375 + _S1375;
    float2  _S1377 = make_float2 (0.0f);
    *&((&_S1377)->y) = _S1374;
    *&((&_S1377)->x) = _S1376;
    dpdpx_7->primal_0 = (*dpdpx_7).primal_0;
    dpdpx_7->differential_0 = _S1377;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_length_impl_0(DiffPair_0 * _S1378, DiffPair_float_0 * _S1379, s_bwd_prop_s_bwd_length_impl_Intermediates_0 * _s_diff_ctx_10)
{
    DiffPair_0 _S1380 = *_S1378;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1381 = _s_diff_ctx_10->_S1205;
    float2  _S1382 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1383 = { _S1382, _S1382 };
    DiffPair_float_0 _S1384 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_sqrt_Intermediates_0 _S1385 = { _S1384 };
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1386;
    (&_S1386)->_S1207 = _S1383;
    (&_S1386)->_S1208 = _S1384;
    (&_S1386)->_S1209 = _S1385;
    s_primal_ctx_s_bwd_prop_length_impl_0(&_S1381, (*_S1379).primal_0, &_S1386);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1387 = { (*_S1378).differential_0.primal_0, (*_S1378).differential_0.differential_0 };
    DiffPair_0 _S1388;
    (&_S1388)->primal_0 = _s_diff_ctx_10->_S1205;
    (&_S1388)->differential_0 = _S1387;
    DiffPair_float_0 _S1389;
    (&_S1389)->primal_0 = (*_S1379).primal_0;
    (&_S1389)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_length_impl_Intermediates_0 _S1390 = _S1386;
    s_bwd_prop_s_bwd_prop_length_impl_0(&_S1388, &_S1389, &_S1390);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1391;
    (&_S1391)->primal_0 = (&_s_diff_ctx_10->_S1205)->primal_0;
    (&_S1391)->differential_0 = _S1382;
    s_bwd_prop_s_primal_ctx_length_impl_0(&_S1391, 0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1392 = { _S1388.differential_0.primal_0 + _S1391.differential_0, _S1388.differential_0.differential_0 };
    _S1379->primal_0 = (*_S1379).primal_0;
    _S1379->differential_0 = _S1389.differential_0;
    _S1378->primal_0 = _S1380.primal_0;
    _S1378->differential_0 = _S1392;
    return;
}

inline __device__ void s_bwd_prop_d_atan2_0(DiffPair_1 * dpdpy_1, DiffPair_1 * dpdpx_8, DiffPair_float_0 * dpdOut_3, s_bwd_prop_d_atan2_Intermediates_0 * _s_diff_ctx_11)
{
    DiffPair_1 _S1393 = *dpdpy_1;
    DiffPair_1 _S1394 = *dpdpx_8;
    float _S1395 = - (&_s_diff_ctx_11->_S1188)->primal_0;
    float _S1396 = (&_s_diff_ctx_11->_S1189)->primal_0 * (&_s_diff_ctx_11->_S1189)->primal_0 + (&_s_diff_ctx_11->_S1188)->primal_0 * (&_s_diff_ctx_11->_S1188)->primal_0;
    float _S1397 = _S1396 * _S1396;
    float _S1398 = (*dpdOut_3).primal_0 * (*dpdpy_1).differential_0.differential_0 / _S1397;
    float _S1399 = (&_s_diff_ctx_11->_S1189)->primal_0 * - _S1398;
    float _S1400 = (&_s_diff_ctx_11->_S1188)->primal_0 * _S1399;
    float _S1401 = (&_s_diff_ctx_11->_S1189)->primal_0 * _S1399;
    float _S1402 = (*dpdOut_3).primal_0 * (*dpdpx_8).differential_0.differential_0 / _S1397;
    float _S1403 = _S1395 * - _S1402;
    float _S1404 = (&_s_diff_ctx_11->_S1188)->primal_0 * _S1403;
    float _S1405 = (&_s_diff_ctx_11->_S1189)->primal_0 * _S1403;
    DiffPair_float_0 dpdpx_9 = { _S1405 + _S1405 + ((*dpdpx_8).differential_0.primal_0 + (_S1401 + _S1401 + _S1396 * _S1398)), 0.0f };
    DiffPair_float_0 dpdpy_2 = { _S1400 + _S1400 + (*dpdpy_1).differential_0.primal_0 + _S1404 + _S1404 + - (_S1396 * _S1402), 0.0f };
    float _S1406 = (&_s_diff_ctx_11->_S1189)->primal_0 / _S1396 * (*dpdpy_1).differential_0.differential_0 + _S1395 / _S1396 * (*dpdpx_8).differential_0.differential_0;
    dpdOut_3->primal_0 = (*dpdOut_3).primal_0;
    dpdOut_3->differential_0 = _S1406;
    dpdpy_1->primal_0 = _S1393.primal_0;
    dpdpy_1->differential_0 = dpdpy_2;
    dpdpx_8->primal_0 = _S1394.primal_0;
    dpdpx_8->differential_0 = dpdpx_9;
    return;
}

inline __device__ void s_bwd_prop_s_bwd_prop_atan2_0(DiffPair_1 * _S1407, DiffPair_1 * _S1408, DiffPair_float_0 * _S1409, s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 * _s_diff_ctx_12)
{
    DiffPair_1 _S1410 = *_S1407;
    DiffPair_1 _S1411 = *_S1408;
    DiffPair_float_0 _S1412 = _s_diff_ctx_12->_S1176;
    DiffPair_float_0 _S1413 = _s_diff_ctx_12->_S1177;
    DiffPair_float_0 _S1414 = { 0.0f, 0.0f };
    s_bwd_prop_d_atan2_Intermediates_0 _S1415;
    (&_S1415)->_S1188 = _S1414;
    (&_S1415)->_S1189 = _S1414;
    s_primal_ctx_d_atan2_0(&_S1412, &_S1413, (*_S1409).primal_0, &_S1415);
    DiffPair_float_0 _S1416 = { (*_S1408).differential_0.primal_0, (*_S1408).differential_0.differential_0 };
    DiffPair_float_0 _S1417 = { (*_S1407).differential_0.primal_0, (*_S1407).differential_0.differential_0 };
    DiffPair_1 _S1418;
    (&_S1418)->primal_0 = _s_diff_ctx_12->_S1176;
    (&_S1418)->differential_0 = _S1417;
    DiffPair_1 _S1419;
    (&_S1419)->primal_0 = _s_diff_ctx_12->_S1177;
    (&_S1419)->differential_0 = _S1416;
    DiffPair_float_0 _S1420;
    (&_S1420)->primal_0 = (*_S1409).primal_0;
    (&_S1420)->differential_0 = 0.0f;
    s_bwd_prop_d_atan2_Intermediates_0 _S1421 = _S1415;
    s_bwd_prop_d_atan2_0(&_S1418, &_S1419, &_S1420, &_S1421);
    DiffPair_float_0 _S1422 = { _S1419.differential_0.primal_0, _S1419.differential_0.differential_0 };
    DiffPair_float_0 _S1423 = { _S1418.differential_0.primal_0, _S1418.differential_0.differential_0 };
    _S1409->primal_0 = (*_S1409).primal_0;
    _S1409->differential_0 = _S1420.differential_0;
    _S1407->primal_0 = _S1410.primal_0;
    _S1407->differential_0 = _S1423;
    _S1408->primal_0 = _S1411.primal_0;
    _S1408->differential_0 = _S1422;
    return;
}

inline __device__ void s_bwd_prop_s_primal_ctx_atan2_0(DiffPair_float_0 * _S1424, DiffPair_float_0 * _S1425, float _s_dOut_5)
{
    DiffPair_float_0 _S1426;
    (&_S1426)->primal_0 = (*_S1424).primal_0;
    (&_S1426)->differential_0 = 0.0f;
    DiffPair_float_0 _S1427;
    (&_S1427)->primal_0 = (*_S1425).primal_0;
    (&_S1427)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1426, &_S1427, _s_dOut_5);
    _S1425->primal_0 = (*_S1425).primal_0;
    _S1425->differential_0 = _S1427.differential_0;
    _S1424->primal_0 = (*_S1424).primal_0;
    _S1424->differential_0 = _S1426.differential_0;
    return;
}

inline __device__ void s_bwd_prop_fisheye_proj_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean3d_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpcov3d_1, DiffPair_vectorx3Cfloatx2C4x3E_0 * dpintrins_1, DiffPair_arrayx3Cfloatx2C10x3E_0 * dpdist_coeffs_1, Matrix<float, 2, 2>  dpcov2d_1, float2  dpmean2d_1, s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 * _s_diff_ctx_13)
{
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1428 = *dpcov3d_1;
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S1429 = *dpintrins_1;
    FixedArray<float, 10>  _S1430 = dpdist_coeffs_1->primal_0;
    float2  _S1431 = make_float2 (0.0f);
    float2  _S1432 = float2 {(*dpmean3d_1).primal_0.x, (*dpmean3d_1).primal_0.y};
    float _S1433 = length_0(_S1432);
    float _S1434 = (*dpmean3d_1).primal_0.z;
    float _S1435 = s_primal_ctx_atan2_0(_S1433, _S1434);
    bool _S1436 = _S1435 < 0.00100000004749745f;
    float k_5;
    float _S1437;
    float _S1438;
    float _S1439;
    if(_S1436)
    {
        float _S1440 = 1.0f - _S1435 * _S1435 / 3.0f;
        float _S1441 = _S1434 * _S1434;
        k_5 = _S1440 / _S1434;
        _S1437 = 0.0f;
        _S1438 = _S1441;
        _S1439 = _S1440;
    }
    else
    {
        float _S1442 = _S1433 * _S1433;
        k_5 = _S1435 / _S1433;
        _S1437 = _S1442;
        _S1438 = 0.0f;
        _S1439 = 0.0f;
    }
    float2  _S1443 = make_float2 (k_5);
    float2  _S1444 = _S1432 * make_float2 (k_5);
    float u_41 = _S1444.x;
    float v_41 = _S1444.y;
    float r2_41 = u_41 * u_41 + v_41 * v_41;
    float _S1445 = _S1430[int(2)] + r2_41 * _S1430[int(3)];
    float _S1446 = _S1430[int(1)] + r2_41 * _S1445;
    float _S1447 = _S1430[int(0)] + r2_41 * _S1446;
    float radial_0 = 1.0f + r2_41 * _S1447;
    float2  _S1448 = make_float2 (radial_0);
    float _S1449 = 2.0f * _S1430[int(4)];
    float _S1450 = _S1449 * u_41;
    float _S1451 = 2.0f * u_41;
    float _S1452 = r2_41 + _S1451 * u_41;
    float _S1453 = 2.0f * _S1430[int(5)];
    float _S1454 = _S1453 * u_41;
    float _S1455 = 2.0f * v_41;
    float _S1456 = r2_41 + _S1455 * v_41;
    float2  _S1457 = _S1444 * make_float2 (radial_0) + make_float2 (_S1450 * v_41 + _S1430[int(5)] * _S1452 + _S1430[int(6)] * r2_41, _S1454 * v_41 + _S1430[int(4)] * _S1456 + _S1430[int(7)] * r2_41);
    float _S1458 = _S1457.x;
    float _S1459 = _S1457.y;
    float2  _S1460 = _S1457 + make_float2 (_S1430[int(8)] * _S1458 + _S1430[int(9)] * _S1459, 0.0f);
    float fx_19 = _S1429.primal_0.x;
    float fy_19 = _S1429.primal_0.y;
    float _S1461 = _S1460.x;
    float _S1462 = _S1460.y;
    Matrix<float, 2, 3>  J_10 = makeMatrix<float, 2, 3> (0.0f);
    float _S1463 = s_primal_ctx_s_primal_ctx_atan2_0(_S1433, _S1434);
    bool _S1464 = _S1463 < 0.00100000004749745f;
    float _S1465;
    float _S1466;
    float _S1467;
    if(_S1464)
    {
        float _S1468 = 1.0f - _S1463 * _S1463 / 3.0f;
        float _S1469 = _S1434 * _S1434;
        k_5 = _S1468 / _S1434;
        _S1465 = 0.0f;
        _S1466 = _S1469;
        _S1467 = _S1468;
    }
    else
    {
        float _S1470 = _S1433 * _S1433;
        k_5 = _S1463 / _S1433;
        _S1465 = _S1470;
        _S1466 = 0.0f;
        _S1467 = 0.0f;
    }
    float2  _S1471 = make_float2 (k_5);
    float2  _S1472 = _S1432 * make_float2 (k_5);
    float u_42 = _S1472.x;
    float v_42 = _S1472.y;
    float r2_42 = u_42 * u_42 + v_42 * v_42;
    float _S1473 = _S1430[int(2)] + r2_42 * _S1430[int(3)];
    float _S1474 = _S1430[int(1)] + r2_42 * _S1473;
    float _S1475 = _S1430[int(0)] + r2_42 * _S1474;
    float2  _S1476 = make_float2 (1.0f + r2_42 * _S1475);
    float _S1477 = _S1449 * u_42;
    float _S1478 = 2.0f * u_42;
    float _S1479 = _S1453 * u_42;
    float _S1480 = 2.0f * v_42;
    float2  _S1481 = make_float2 (fx_19, 0.0f) + make_float2 (_S1430[int(8)] * fx_19, _S1430[int(9)] * fx_19);
    float2  _S1482 = _S1472 * _S1481;
    float _S1483 = _S1430[int(4)] * _S1481.y;
    float _S1484 = v_42 * _S1481.y;
    float _S1485 = _S1430[int(5)] * _S1481.x;
    float _S1486 = v_42 * _S1481.x;
    float _S1487 = _S1482.x + _S1482.y;
    float _S1488 = r2_42 * _S1487;
    float _S1489 = r2_42 * _S1488;
    float _S1490 = r2_42 * _S1489;
    float _S1491 = _S1430[int(7)] * _S1481.y + _S1483 + _S1430[int(6)] * _S1481.x + _S1485 + _S1475 * _S1487 + _S1474 * _S1488 + _S1473 * _S1489 + _S1430[int(3)] * _S1490;
    float _S1492 = v_42 * _S1491;
    float _S1493 = u_42 * _S1491;
    float2  _S1494 = _S1476 * _S1481 + make_float2 (_S1453 * _S1484 + _S1478 * _S1485 + 2.0f * (u_42 * _S1485) + _S1449 * _S1486 + _S1493 + _S1493, _S1480 * _S1483 + 2.0f * (v_42 * _S1483) + _S1479 * _S1481.y + _S1477 * _S1481.x + _S1492 + _S1492);
    float2  _S1495 = _S1432 * _S1494;
    float2  _S1496 = _S1471 * _S1494;
    float _S1497 = _S1495.x + _S1495.y;
    float k_6;
    float _S1498;
    float _S1499;
    float _S1500;
    float _S1501;
    float _S1502;
    float _S1503;
    float _S1504;
    float _S1505;
    if(_S1464)
    {
        float _S1506 = _S1497 / _S1466;
        float _S1507 = _S1466 * _S1466;
        float _S1508 = - _S1506;
        float _S1509 = _S1467 * _S1508;
        float _S1510 = 0.3333333432674408f * - (_S1434 * _S1506);
        float _S1511 = _S1463 * _S1510;
        k_5 = _S1511 + _S1511;
        k_6 = _S1509;
        _S1498 = 0.0f;
        _S1499 = 0.0f;
        _S1500 = 0.0f;
        _S1501 = 0.0f;
        _S1502 = _S1510;
        _S1503 = _S1506;
        _S1504 = _S1508;
        _S1505 = _S1507;
    }
    else
    {
        float _S1512 = _S1497 / _S1465;
        float _S1513 = _S1465 * _S1465;
        float _S1514 = - _S1512;
        float _S1515 = _S1463 * _S1514;
        k_5 = _S1433 * _S1512;
        k_6 = 0.0f;
        _S1498 = _S1515;
        _S1499 = _S1512;
        _S1500 = _S1514;
        _S1501 = _S1513;
        _S1502 = 0.0f;
        _S1503 = 0.0f;
        _S1504 = 0.0f;
        _S1505 = 0.0f;
    }
    DiffPair_float_0 _S1516 = { _S1433, 0.0f };
    DiffPair_float_0 _S1517 = { _S1434, 0.0f };
    float _S1518 = (&_s_diff_ctx_13->_S1179)->differential_0 + k_6;
    float _S1519 = (&_s_diff_ctx_13->_S1178)->differential_0 + _S1498;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1520 = { _S1432, _S1431 };
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1521;
    (&_S1521)->primal_0 = _S1432;
    (&_S1521)->differential_0 = _S1431;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1522 = { _S1431, _S1431 };
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1523;
    (&_S1523)->_S1205 = _S1522;
    s_primal_ctx_s_bwd_length_impl_0(&_S1521, _S1519, &_S1523);
    float2  _S1524 = _S1521.differential_0 + _S1496;
    float3  _S1525 = make_float3 (_S1524.x, _S1524.y, _S1518);
    Matrix<float, 2, 3>  _S1526 = J_10;
    _S1526[int(0)] = _S1525;
    float _S1527;
    float _S1528;
    if(_S1464)
    {
        float _S1529 = 1.0f - _S1463 * _S1463 / 3.0f;
        float _S1530 = _S1434 * _S1434;
        k_6 = _S1529 / _S1434;
        _S1498 = 0.0f;
        _S1527 = _S1530;
        _S1528 = _S1529;
    }
    else
    {
        float _S1531 = _S1433 * _S1433;
        k_6 = _S1463 / _S1433;
        _S1498 = _S1531;
        _S1527 = 0.0f;
        _S1528 = 0.0f;
    }
    float2  _S1532 = make_float2 (k_6);
    float2  _S1533 = _S1432 * make_float2 (k_6);
    float u_43 = _S1533.x;
    float v_43 = _S1533.y;
    float r2_43 = u_43 * u_43 + v_43 * v_43;
    float _S1534 = _S1430[int(2)] + r2_43 * _S1430[int(3)];
    float _S1535 = _S1430[int(1)] + r2_43 * _S1534;
    float _S1536 = _S1430[int(0)] + r2_43 * _S1535;
    float2  _S1537 = make_float2 (1.0f + r2_43 * _S1536);
    float _S1538 = _S1449 * u_43;
    float _S1539 = 2.0f * u_43;
    float _S1540 = _S1453 * u_43;
    float _S1541 = 2.0f * v_43;
    float2  _S1542 = make_float2 (0.0f, fy_19);
    float2  _S1543 = _S1533 * _S1542;
    float _S1544 = _S1430[int(4)] * fy_19;
    float _S1545 = v_43 * fy_19;
    float _S1546 = _S1543.x + _S1543.y;
    float _S1547 = r2_43 * _S1546;
    float _S1548 = r2_43 * _S1547;
    float _S1549 = r2_43 * _S1548;
    float _S1550 = _S1430[int(7)] * fy_19 + _S1544 + _S1536 * _S1546 + _S1535 * _S1547 + _S1534 * _S1548 + _S1430[int(3)] * _S1549;
    float _S1551 = v_43 * _S1550;
    float _S1552 = u_43 * _S1550;
    float2  _S1553 = _S1537 * _S1542 + make_float2 (_S1453 * _S1545 + _S1552 + _S1552, _S1541 * _S1544 + 2.0f * (v_43 * _S1544) + _S1540 * fy_19 + _S1551 + _S1551);
    float2  _S1554 = _S1432 * _S1553;
    float2  _S1555 = _S1532 * _S1553;
    float _S1556 = _S1554.x + _S1554.y;
    float _S1557;
    float _S1558;
    float _S1559;
    float _S1560;
    float _S1561;
    float _S1562;
    float _S1563;
    float _S1564;
    float _S1565;
    if(_S1464)
    {
        float _S1566 = _S1556 / _S1527;
        float _S1567 = _S1527 * _S1527;
        float _S1568 = - _S1566;
        float _S1569 = _S1528 * _S1568;
        float _S1570 = 0.3333333432674408f * - (_S1434 * _S1566);
        float _S1571 = _S1463 * _S1570;
        k_6 = _S1571 + _S1571;
        _S1557 = _S1569;
        _S1558 = 0.0f;
        _S1559 = 0.0f;
        _S1560 = 0.0f;
        _S1561 = 0.0f;
        _S1562 = _S1570;
        _S1563 = _S1566;
        _S1564 = _S1568;
        _S1565 = _S1567;
    }
    else
    {
        float _S1572 = _S1556 / _S1498;
        float _S1573 = _S1498 * _S1498;
        float _S1574 = - _S1572;
        float _S1575 = _S1463 * _S1574;
        k_6 = _S1433 * _S1572;
        _S1557 = 0.0f;
        _S1558 = _S1575;
        _S1559 = _S1572;
        _S1560 = _S1574;
        _S1561 = _S1573;
        _S1562 = 0.0f;
        _S1563 = 0.0f;
        _S1564 = 0.0f;
        _S1565 = 0.0f;
    }
    float _S1576 = (&_s_diff_ctx_13->_S1182)->differential_0 + _S1557;
    float _S1577 = (&_s_diff_ctx_13->_S1181)->differential_0 + _S1558;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1578;
    (&_S1578)->primal_0 = _S1432;
    (&_S1578)->differential_0 = _S1431;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1579;
    (&_S1579)->_S1205 = _S1522;
    s_primal_ctx_s_bwd_length_impl_0(&_S1578, _S1577, &_S1579);
    float2  _S1580 = _S1578.differential_0 + _S1555;
    float3  _S1581 = make_float3 (_S1580.x, _S1580.y, _S1576);
    _S1526[int(1)] = _S1581;
    Matrix<float, 3, 2>  _S1582 = transpose_1(_S1526);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1583;
    (&_S1583)->primal_0 = s_primal_ctx_mul_2(_S1526, _S1428.primal_0);
    (&_S1583)->differential_0 = J_10;
    Matrix<float, 3, 2>  _S1584 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S1585;
    (&_S1585)->primal_0 = _S1582;
    (&_S1585)->differential_0 = _S1584;
    s_bwd_prop_mul_2(&_S1583, &_S1585, dpcov2d_1);
    Matrix<float, 2, 3>  _S1586 = transpose_2(_S1585.differential_0);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S1587;
    (&_S1587)->primal_0 = _S1526;
    (&_S1587)->differential_0 = J_10;
    Matrix<float, 3, 3>  _S1588 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1589;
    (&_S1589)->primal_0 = _S1428.primal_0;
    (&_S1589)->differential_0 = _S1588;
    s_bwd_prop_mul_3(&_S1587, &_S1589, _S1583.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S1590 = _S1589;
    Matrix<float, 2, 3>  _S1591 = _S1586 + _S1587.differential_0;
    float2  _S1592 = _S1431;
    *&((&_S1592)->y) = _S1591.rows[int(1)].y;
    *&((&_S1592)->x) = _S1591.rows[int(1)].x;
    float2  _S1593 = _S1592;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1594 = { _S1431, _S1592 };
    DiffPair_0 _S1595;
    (&_S1595)->primal_0 = _S1520;
    (&_S1595)->differential_0 = _S1594;
    DiffPair_float_0 _S1596;
    (&_S1596)->primal_0 = _S1577;
    (&_S1596)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1597 = _S1579;
    s_bwd_prop_s_bwd_length_impl_0(&_S1595, &_S1596, &_S1597);
    DiffPair_0 _S1598 = _S1595;
    DiffPair_float_0 _S1599 = _S1596;
    DiffPair_float_0 _S1600 = { 0.0f, _S1591.rows[int(1)].z };
    DiffPair_float_0 _S1601 = { 0.0f, _S1596.differential_0 };
    DiffPair_1 _S1602;
    (&_S1602)->primal_0 = _S1516;
    (&_S1602)->differential_0 = _S1601;
    DiffPair_1 _S1603;
    (&_S1603)->primal_0 = _S1517;
    (&_S1603)->differential_0 = _S1600;
    DiffPair_float_0 _S1604;
    (&_S1604)->primal_0 = k_6;
    (&_S1604)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1602, &_S1603, &_S1604, &_s_diff_ctx_13->_S1183);
    DiffPair_1 _S1605 = _S1602;
    DiffPair_1 _S1606 = _S1603;
    DiffPair_float_0 _S1607 = _S1604;
    if(_S1464)
    {
        float _S1608 = _S1607.differential_0 + _S1607.differential_0;
        float _S1609 = _S1562 * _S1608;
        float _S1610 = - (0.3333333432674408f * (_S1463 * _S1608));
        float _S1611 = _S1564 * _S1591.rows[int(1)].z;
        float _S1612 = (_S1434 * _S1610 + - (_S1528 * _S1591.rows[int(1)].z)) / _S1565;
        float _S1613 = _S1556 * - _S1612;
        float _S1614 = _S1563 * _S1610 + _S1606.differential_0.primal_0;
        k_6 = _S1527 * _S1612;
        _S1557 = 0.0f;
        _S1558 = _S1613;
        _S1559 = _S1611;
        _S1560 = _S1605.differential_0.primal_0;
        _S1561 = _S1609;
        _S1562 = _S1614;
    }
    else
    {
        float _S1615 = _S1560 * _S1599.differential_0;
        float _S1616 = (_S1433 * _S1607.differential_0 + - (_S1463 * _S1599.differential_0)) / _S1561;
        float _S1617 = _S1556 * - _S1616;
        float _S1618 = _S1559 * _S1607.differential_0 + _S1605.differential_0.primal_0;
        k_6 = _S1498 * _S1616;
        _S1557 = _S1617;
        _S1558 = 0.0f;
        _S1559 = 0.0f;
        _S1560 = _S1618;
        _S1561 = _S1615;
        _S1562 = _S1606.differential_0.primal_0;
    }
    float2  _S1619 = _S1532 * _S1593;
    float2  _S1620 = _S1553 * _S1593;
    float2  _S1621 = _S1431;
    *&((&_S1621)->y) = k_6;
    *&((&_S1621)->x) = k_6;
    float2  _S1622 = _S1553 * _S1621;
    float2  _S1623 = _S1619 + _S1432 * _S1621;
    float _S1624 = _S1623.x;
    float _S1625 = _S1624 + _S1624;
    float _S1626 = _S1550 * _S1625;
    float _S1627 = _S1623.y + _S1623.y;
    float _S1628 = _S1550 * _S1627;
    float _S1629 = u_43 * _S1625 + v_43 * _S1627;
    float _S1630 = _S1430[int(3)] * _S1629;
    float _S1631 = _S1549 * _S1629;
    float _S1632 = _S1548 * _S1629;
    float _S1633 = _S1548 * _S1630;
    float _S1634 = _S1547 * _S1629;
    float _S1635 = _S1534 * _S1629 + r2_43 * _S1630;
    float _S1636 = _S1547 * _S1635;
    float _S1637 = _S1546 * _S1629;
    float _S1638 = _S1535 * _S1629 + r2_43 * _S1635;
    float _S1639 = _S1546 * _S1638;
    float _S1640 = _S1536 * _S1629 + r2_43 * _S1638;
    float _S1641 = v_43 * (_S1449 * _S1623.x);
    float _S1642 = _S1538 * _S1623.y;
    float _S1643 = _S1430[int(5)] * (_S1629 + u_43 * (2.0f * _S1623.x) + _S1539 * _S1623.x);
    float _S1644 = _S1430[int(6)] * _S1629;
    float _S1645 = _S1453 * _S1623.x;
    float _S1646 = _S1545 * _S1623.x;
    float _S1647 = v_43 * _S1645;
    float _S1648 = fy_19 * _S1645;
    float _S1649 = _S1540 * _S1623.y;
    float _S1650 = fy_19 * _S1623.y;
    float _S1651 = 2.0f * _S1623.y;
    float _S1652 = _S1544 * _S1651;
    float _S1653 = _S1544 * _S1623.y;
    float _S1654 = _S1629 + v_43 * _S1651 + _S1541 * _S1623.y;
    float _S1655 = _S1430[int(4)] * _S1654;
    float _S1656 = fy_19 * _S1654;
    float _S1657 = _S1430[int(7)] * _S1629;
    float _S1658 = fy_19 * _S1629;
    float2  _S1659 = _S1537 * _S1623;
    float2  _S1660 = _S1542 * _S1623;
    float2  _S1661 = _S1431;
    *&((&_S1661)->y) = _S1640;
    *&((&_S1661)->x) = _S1640;
    float2  _S1662 = _S1542 * _S1661;
    float _S1663 = _S1647 + _S1649 + _S1655 + _S1657;
    float _S1664 = _S1641 + _S1642 + _S1643 + _S1644;
    float2  _S1665 = _S1659 + _S1533 * _S1661;
    float2  _S1666 = _S1431;
    *&((&_S1666)->y) = _S1663;
    *&((&_S1666)->x) = _S1664;
    float2  _S1667 = _S1665 + _S1666;
    float _S1668 = _S1660.x + _S1660.y;
    float _S1669 = _S1637 + r2_43 * _S1668;
    float _S1670 = _S1634 + r2_43 * _S1669;
    float _S1671 = _S1632 + r2_43 * _S1670;
    float _S1672 = _S1633 + _S1636 + _S1639 + _S1536 * _S1668 + _S1535 * _S1669 + _S1534 * _S1670 + _S1430[int(3)] * _S1671;
    float _S1673 = v_43 * _S1672;
    float _S1674 = u_43 * _S1672;
    float2  _S1675 = _S1622 + _S1598.differential_0.primal_0;
    float2  _S1676 = _S1662 + make_float2 (_S1626 + _S1453 * _S1650 + _S1674 + _S1674, _S1628 + _S1648 + _S1652 + 2.0f * _S1653 + _S1673 + _S1673);
    float _S1677 = _S1646 + u_43 * _S1650;
    float _S1678 = _S1631 + r2_43 * _S1671;
    float2  _S1679 = _S1432 * _S1676;
    float _S1680 = _S1620.x + _S1620.y + _S1679.x + _S1679.y;
    float2  _S1681 = _S1532 * _S1676 + _S1675;
    if(_S1464)
    {
        float _S1682 = _S1434 * _S1558;
        float _S1683 = _S1680 / _S1527;
        float _S1684 = _S1463 * (0.3333333432674408f * - (_S1559 + _S1434 * _S1683));
        float _S1685 = _S1682 + _S1682 + _S1528 * - _S1683 + _S1562;
        k_6 = _S1684 + _S1684 + _S1561;
        _S1498 = _S1685;
        _S1527 = _S1560;
    }
    else
    {
        float _S1686 = _S1433 * _S1557;
        float _S1687 = _S1680 / _S1498;
        float _S1688 = _S1686 + _S1686 + _S1463 * - _S1687 + _S1560;
        k_6 = _S1433 * _S1687 + _S1561;
        _S1498 = _S1562;
        _S1527 = _S1688;
    }
    DiffPair_float_0 _S1689;
    (&_S1689)->primal_0 = _S1433;
    (&_S1689)->differential_0 = 0.0f;
    DiffPair_float_0 _S1690;
    (&_S1690)->primal_0 = _S1434;
    (&_S1690)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1689, &_S1690, k_6);
    float _S1691 = _S1690.differential_0 + _S1498;
    float _S1692 = _S1689.differential_0 + _S1527;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1693;
    (&_S1693)->primal_0 = _S1432;
    (&_S1693)->differential_0 = _S1431;
    s_bwd_length_impl_0(&_S1693, _S1692);
    float2  _S1694 = _S1693.differential_0 + _S1681;
    float3  _S1695 = make_float3 (_S1694.x, _S1694.y, _S1691);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1696;
    (&_S1696)->primal_0 = _S1432;
    (&_S1696)->differential_0 = _S1431;
    s_bwd_length_impl_0(&_S1696, 0.0f);
    float3  _S1697 = _S1695 + make_float3 (_S1696.differential_0.x, _S1696.differential_0.y, 0.0f);
    float2  _S1698 = _S1431;
    *&((&_S1698)->y) = _S1591.rows[int(0)].y;
    *&((&_S1698)->x) = _S1591.rows[int(0)].x;
    float2  _S1699 = _S1698;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1700 = { _S1431, _S1698 };
    DiffPair_0 _S1701;
    (&_S1701)->primal_0 = _S1520;
    (&_S1701)->differential_0 = _S1700;
    DiffPair_float_0 _S1702;
    (&_S1702)->primal_0 = _S1519;
    (&_S1702)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_length_impl_Intermediates_0 _S1703 = _S1523;
    s_bwd_prop_s_bwd_length_impl_0(&_S1701, &_S1702, &_S1703);
    DiffPair_0 _S1704 = _S1701;
    DiffPair_float_0 _S1705 = _S1702;
    DiffPair_float_0 _S1706 = { 0.0f, _S1591.rows[int(0)].z };
    DiffPair_float_0 _S1707 = { 0.0f, _S1702.differential_0 };
    DiffPair_1 _S1708;
    (&_S1708)->primal_0 = _S1516;
    (&_S1708)->differential_0 = _S1707;
    DiffPair_1 _S1709;
    (&_S1709)->primal_0 = _S1517;
    (&_S1709)->differential_0 = _S1706;
    DiffPair_float_0 _S1710;
    (&_S1710)->primal_0 = k_5;
    (&_S1710)->differential_0 = 0.0f;
    s_bwd_prop_s_bwd_prop_atan2_0(&_S1708, &_S1709, &_S1710, &_s_diff_ctx_13->_S1180);
    DiffPair_1 _S1711 = _S1708;
    DiffPair_1 _S1712 = _S1709;
    DiffPair_float_0 _S1713 = _S1710;
    if(_S1464)
    {
        float _S1714 = _S1713.differential_0 + _S1713.differential_0;
        float _S1715 = _S1502 * _S1714;
        float _S1716 = - (0.3333333432674408f * (_S1463 * _S1714));
        float _S1717 = _S1504 * _S1591.rows[int(0)].z;
        float _S1718 = (_S1434 * _S1716 + - (_S1467 * _S1591.rows[int(0)].z)) / _S1505;
        float _S1719 = _S1497 * - _S1718;
        float _S1720 = _S1503 * _S1716 + _S1712.differential_0.primal_0;
        k_5 = _S1466 * _S1718;
        k_6 = 0.0f;
        _S1498 = _S1719;
        _S1499 = _S1717;
        _S1500 = _S1711.differential_0.primal_0;
        _S1501 = _S1715;
        _S1502 = _S1720;
    }
    else
    {
        float _S1721 = _S1500 * _S1705.differential_0;
        float _S1722 = (_S1433 * _S1713.differential_0 + - (_S1463 * _S1705.differential_0)) / _S1501;
        float _S1723 = _S1497 * - _S1722;
        float _S1724 = _S1499 * _S1713.differential_0 + _S1711.differential_0.primal_0;
        k_5 = _S1465 * _S1722;
        k_6 = _S1723;
        _S1498 = 0.0f;
        _S1499 = 0.0f;
        _S1500 = _S1724;
        _S1501 = _S1721;
        _S1502 = _S1712.differential_0.primal_0;
    }
    float2  _S1725 = _S1471 * _S1699;
    float2  _S1726 = _S1494 * _S1699;
    float2  _S1727 = _S1431;
    *&((&_S1727)->y) = k_5;
    *&((&_S1727)->x) = k_5;
    float2  _S1728 = _S1494 * _S1727;
    float2  _S1729 = _S1725 + _S1432 * _S1727;
    float _S1730 = _S1729.x;
    float _S1731 = _S1730 + _S1730;
    float _S1732 = _S1491 * _S1731;
    float _S1733 = _S1729.y + _S1729.y;
    float _S1734 = _S1491 * _S1733;
    float _S1735 = u_42 * _S1731 + v_42 * _S1733;
    float _S1736 = _S1430[int(3)] * _S1735;
    float _S1737 = _S1490 * _S1735;
    float _S1738 = _S1489 * _S1735;
    float _S1739 = _S1489 * _S1736;
    float _S1740 = _S1488 * _S1735;
    float _S1741 = _S1473 * _S1735 + r2_42 * _S1736;
    float _S1742 = _S1488 * _S1741;
    float _S1743 = _S1487 * _S1735;
    float _S1744 = _S1474 * _S1735 + r2_42 * _S1741;
    float _S1745 = _S1487 * _S1744;
    float _S1746 = _S1475 * _S1735 + r2_42 * _S1744;
    float _S1747 = _S1449 * _S1729.x;
    float _S1748 = _S1486 * _S1729.x;
    float _S1749 = v_42 * _S1747;
    float _S1750 = _S1481.x * _S1747;
    float _S1751 = _S1477 * _S1729.y;
    float _S1752 = _S1481.x * _S1729.y;
    float _S1753 = 2.0f * _S1729.x;
    float _S1754 = _S1485 * _S1753;
    float _S1755 = _S1485 * _S1729.x;
    float _S1756 = _S1735 + u_42 * _S1753 + _S1478 * _S1729.x;
    float _S1757 = _S1430[int(5)] * _S1756;
    float _S1758 = _S1481.x * _S1756;
    float _S1759 = _S1430[int(6)] * _S1735;
    float _S1760 = _S1481.x * _S1735;
    float _S1761 = _S1453 * _S1729.x;
    float _S1762 = _S1484 * _S1729.x;
    float _S1763 = v_42 * _S1761;
    float _S1764 = _S1481.y * _S1761;
    float _S1765 = _S1479 * _S1729.y;
    float _S1766 = _S1481.y * _S1729.y;
    float _S1767 = 2.0f * _S1729.y;
    float _S1768 = _S1483 * _S1767;
    float _S1769 = _S1483 * _S1729.y;
    float _S1770 = _S1735 + v_42 * _S1767 + _S1480 * _S1729.y;
    float _S1771 = _S1430[int(4)] * _S1770;
    float _S1772 = _S1481.y * _S1770;
    float _S1773 = _S1430[int(7)] * _S1735;
    float _S1774 = _S1481.y * _S1735;
    float2  _S1775 = _S1476 * _S1729;
    float2  _S1776 = _S1481 * _S1729;
    float2  _S1777 = _S1431;
    *&((&_S1777)->y) = _S1746;
    *&((&_S1777)->x) = _S1746;
    float2  _S1778 = _S1481 * _S1777;
    float _S1779 = _S1763 + _S1765 + _S1771 + _S1773;
    float _S1780 = _S1749 + _S1751 + _S1757 + _S1759;
    float2  _S1781 = _S1775 + _S1472 * _S1777;
    float2  _S1782 = _S1431;
    *&((&_S1782)->y) = _S1779;
    *&((&_S1782)->x) = _S1780;
    float2  _S1783 = _S1781 + _S1782;
    float _S1784 = fx_19 * _S1783.x;
    float _S1785 = fx_19 * _S1783.y;
    float _S1786 = _S1776.x + _S1776.y;
    float _S1787 = _S1743 + r2_42 * _S1786;
    float _S1788 = _S1740 + r2_42 * _S1787;
    float _S1789 = _S1738 + r2_42 * _S1788;
    float _S1790 = _S1739 + _S1742 + _S1745 + _S1475 * _S1786 + _S1474 * _S1787 + _S1473 * _S1788 + _S1430[int(3)] * _S1789;
    float _S1791 = v_42 * _S1790;
    float _S1792 = u_42 * _S1790;
    float2  _S1793 = _S1728 + _S1704.differential_0.primal_0;
    float _S1794 = _S1774 + _S1658;
    float _S1795 = _S1748 + u_42 * _S1752;
    float _S1796 = _S1430[int(8)] * _S1783.x + _S1430[int(9)] * _S1783.y + _S1783.x;
    float _S1797 = _S1788 + _S1670;
    float _S1798 = _S1787 + _S1669;
    float2  _S1799 = _S1778 + make_float2 (_S1732 + _S1754 + _S1453 * _S1766 + 2.0f * _S1755 + _S1449 * _S1752 + _S1792 + _S1792, _S1734 + _S1750 + _S1764 + _S1768 + 2.0f * _S1769 + _S1791 + _S1791);
    float _S1800 = _S1762 + u_42 * _S1766 + _S1677;
    float _S1801 = _S1737 + r2_42 * _S1789 + _S1678;
    float _S1802 = _S1789 + _S1671;
    float _S1803 = _S1772 + _S1656;
    float2  _S1804 = _S1432 * _S1799;
    float _S1805 = _S1726.x + _S1726.y + _S1804.x + _S1804.y;
    float2  _S1806 = _S1471 * _S1799 + _S1793;
    if(_S1464)
    {
        float _S1807 = _S1434 * _S1498;
        float _S1808 = _S1805 / _S1466;
        float _S1809 = _S1463 * (0.3333333432674408f * - (_S1499 + _S1434 * _S1808));
        float _S1810 = _S1807 + _S1807 + _S1467 * - _S1808 + _S1502;
        k_5 = _S1809 + _S1809 + _S1501;
        _S1465 = _S1810;
        _S1466 = _S1500;
    }
    else
    {
        float _S1811 = _S1433 * k_6;
        float _S1812 = _S1805 / _S1465;
        float _S1813 = _S1811 + _S1811 + _S1463 * - _S1812 + _S1500;
        k_5 = _S1433 * _S1812 + _S1501;
        _S1465 = _S1502;
        _S1466 = _S1813;
    }
    DiffPair_float_0 _S1814;
    (&_S1814)->primal_0 = _S1433;
    (&_S1814)->differential_0 = 0.0f;
    DiffPair_float_0 _S1815;
    (&_S1815)->primal_0 = _S1434;
    (&_S1815)->differential_0 = 0.0f;
    s_bwd_prop_s_primal_ctx_atan2_0(&_S1814, &_S1815, k_5);
    float _S1816 = _S1815.differential_0 + _S1465;
    float _S1817 = _S1814.differential_0 + _S1466;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1818;
    (&_S1818)->primal_0 = _S1432;
    (&_S1818)->differential_0 = _S1431;
    s_bwd_length_impl_0(&_S1818, _S1817);
    float2  _S1819 = _S1818.differential_0 + _S1806;
    float3  _S1820 = make_float3 (_S1819.x, _S1819.y, _S1816);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1821;
    (&_S1821)->primal_0 = _S1432;
    (&_S1821)->differential_0 = _S1431;
    s_bwd_length_impl_0(&_S1821, 0.0f);
    float _S1822 = fx_19 * dpmean2d_1.x;
    float2  _S1823 = make_float2 (_S1822, fy_19 * dpmean2d_1.y) + make_float2 (_S1430[int(8)] * _S1822, _S1430[int(9)] * _S1822);
    float2  _S1824 = _S1444 * _S1823;
    float2  _S1825 = _S1448 * _S1823;
    float _S1826 = _S1430[int(4)] * _S1823.y;
    float _S1827 = v_41 * _S1823.y;
    float _S1828 = _S1430[int(5)] * _S1823.x;
    float _S1829 = v_41 * _S1823.x;
    float _S1830 = _S1824.x + _S1824.y;
    float _S1831 = r2_41 * _S1830;
    float _S1832 = r2_41 * _S1831;
    float _S1833 = r2_41 * _S1832;
    float _S1834 = _S1430[int(7)] * _S1823.y + _S1826 + _S1430[int(6)] * _S1823.x + _S1828 + _S1447 * _S1830 + _S1446 * _S1831 + _S1445 * _S1832 + _S1430[int(3)] * _S1833;
    float _S1835 = v_41 * _S1834;
    float _S1836 = u_41 * _S1834;
    float _S1837 = _S1455 * _S1826 + 2.0f * (v_41 * _S1826) + _S1454 * _S1823.y + _S1450 * _S1823.x + _S1835 + _S1835;
    float _S1838 = _S1453 * _S1827 + _S1451 * _S1828 + 2.0f * (u_41 * _S1828) + _S1449 * _S1829 + _S1836 + _S1836;
    float _S1839 = _S1459 * _S1822 + _S1785;
    float _S1840 = _S1458 * _S1822 + _S1784;
    float _S1841 = r2_41 * _S1823.y + _S1794;
    float _S1842 = r2_41 * _S1823.x + _S1760;
    float _S1843 = 2.0f * (u_41 * _S1827 + _S1800) + _S1452 * _S1823.x + _S1758;
    float _S1844 = _S1456 * _S1823.y + 2.0f * (u_41 * _S1829 + _S1795) + _S1803;
    float _S1845 = r2_41 * _S1833 + _S1801;
    float _S1846 = _S1833 + _S1802;
    float _S1847 = _S1832 + _S1797;
    float _S1848 = _S1831 + _S1798;
    float3  _S1849 = _S1820 + make_float3 (_S1821.differential_0.x, _S1821.differential_0.y, 0.0f) + _S1697;
    float4  _S1850 = make_float4 (_S1461 * dpmean2d_1.x + _S1796, _S1462 * dpmean2d_1.y + _S1667.y, dpmean2d_1.x, dpmean2d_1.y);
    FixedArray<float, 10>  _S1851;
    _S1851[int(0)] = 0.0f;
    _S1851[int(1)] = 0.0f;
    _S1851[int(2)] = 0.0f;
    _S1851[int(3)] = 0.0f;
    _S1851[int(4)] = 0.0f;
    _S1851[int(5)] = 0.0f;
    _S1851[int(6)] = 0.0f;
    _S1851[int(7)] = 0.0f;
    _S1851[int(8)] = 0.0f;
    _S1851[int(9)] = 0.0f;
    _S1851[int(9)] = _S1839;
    _S1851[int(8)] = _S1840;
    _S1851[int(7)] = _S1841;
    _S1851[int(6)] = _S1842;
    _S1851[int(5)] = _S1843;
    _S1851[int(4)] = _S1844;
    _S1851[int(3)] = _S1845;
    _S1851[int(2)] = _S1846;
    _S1851[int(1)] = _S1847;
    _S1851[int(0)] = _S1848;
    FixedArray<float, 10>  _S1852 = {
        _S1851[int(0)], _S1851[int(1)], _S1851[int(2)], _S1851[int(3)], _S1851[int(4)], _S1851[int(5)], _S1851[int(6)], _S1851[int(7)], _S1851[int(8)], _S1851[int(9)]
    };
    float2  _S1853 = _S1825 + make_float2 (_S1838, _S1837);
    float2  _S1854 = _S1432 * _S1853;
    float2  _S1855 = _S1443 * _S1853;
    float _S1856 = _S1854.x + _S1854.y;
    if(_S1436)
    {
        float _S1857 = _S1856 / _S1438;
        float _S1858 = _S1439 * - _S1857;
        float _S1859 = _S1435 * (0.3333333432674408f * - (_S1434 * _S1857));
        k_5 = _S1859 + _S1859;
        _S1437 = _S1858;
        _S1438 = 0.0f;
    }
    else
    {
        float _S1860 = _S1856 / _S1437;
        float _S1861 = _S1435 * - _S1860;
        k_5 = _S1433 * _S1860;
        _S1437 = 0.0f;
        _S1438 = _S1861;
    }
    DiffPair_float_0 _S1862;
    (&_S1862)->primal_0 = _S1433;
    (&_S1862)->differential_0 = 0.0f;
    DiffPair_float_0 _S1863;
    (&_S1863)->primal_0 = _S1434;
    (&_S1863)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S1862, &_S1863, k_5);
    float _S1864 = _S1863.differential_0 + _S1437;
    float _S1865 = _S1862.differential_0 + _S1438;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S1866;
    (&_S1866)->primal_0 = _S1432;
    (&_S1866)->differential_0 = _S1431;
    s_bwd_length_impl_0(&_S1866, _S1865);
    float2  _S1867 = _S1866.differential_0 + _S1855;
    dpdist_coeffs_1->primal_0 = dpdist_coeffs_1->primal_0;
    dpdist_coeffs_1->differential_0 = _S1852;
    dpintrins_1->primal_0 = (*dpintrins_1).primal_0;
    dpintrins_1->differential_0 = _S1850;
    dpcov3d_1->primal_0 = (*dpcov3d_1).primal_0;
    dpcov3d_1->differential_0 = _S1590.differential_0;
    float3  _S1868 = _S1849 + make_float3 (_S1867.x, _S1867.y, _S1864);
    dpmean3d_1->primal_0 = (*dpmean3d_1).primal_0;
    dpmean3d_1->differential_0 = _S1868;
    return;
}

inline __device__ void s_bwd_prop_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_12, float _s_dOut_6)
{
    float _S1869 = (*dpx_12).primal_0.x;
    float _S1870 = (*dpx_12).primal_0.y;
    float _S1871 = (*dpx_12).primal_0.z;
    DiffPair_float_0 _S1872;
    (&_S1872)->primal_0 = _S1869 * _S1869 + _S1870 * _S1870 + _S1871 * _S1871;
    (&_S1872)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S1872, _s_dOut_6);
    float _S1873 = (*dpx_12).primal_0.z * _S1872.differential_0;
    float _S1874 = _S1873 + _S1873;
    float _S1875 = (*dpx_12).primal_0.y * _S1872.differential_0;
    float _S1876 = _S1875 + _S1875;
    float _S1877 = (*dpx_12).primal_0.x * _S1872.differential_0;
    float _S1878 = _S1877 + _S1877;
    float3  _S1879 = make_float3 (0.0f);
    *&((&_S1879)->z) = _S1874;
    *&((&_S1879)->y) = _S1876;
    *&((&_S1879)->x) = _S1878;
    dpx_12->primal_0 = (*dpx_12).primal_0;
    dpx_12->differential_0 = _S1879;
    return;
}

inline __device__ void s_bwd_length_impl_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S1880, float _S1881)
{
    s_bwd_prop_length_impl_1(_S1880, _S1881);
    return;
}

inline __device__ void projection_3dgs_fisheye_vjp(bool antialiased_11, float3  mean_12, float4  quat_15, float3  scale_14, float in_opacity_11, FixedArray<float3 , 16>  * sh_coeffs_11, Matrix<float, 3, 3>  R_17, float3  t_14, float fx_20, float fy_20, float cx_15, float cy_15, FixedArray<float, 10>  * dist_coeffs_24, uint image_width_11, uint image_height_11, float2  v_mean2d_1, float v_depth_1, float3  v_conic_1, float v_opacity_1, float3  v_rgb_1, float3  * v_mean_1, float4  * v_quat_1, float3  * v_scale_1, float * v_in_opacity_1, FixedArray<float3 , 16>  * v_sh_coeffs_1, Matrix<float, 3, 3>  * v_R_3, float3  * v_t_2)
{
    Matrix<float, 2, 2>  _S1882 = makeMatrix<float, 2, 2> (0.0f);
    DiffPair_float_0 _S1883 = { 0.0f, 0.0f };
    s_bwd_prop_s_bwd_prop_atan2_Intermediates_0 _S1884 = { _S1883, _S1883 };
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S1885 = { _S1883, _S1883, _S1884, _S1883, _S1883, _S1884 };
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1886;
    (&_S1886)->_S1184 = _S1882;
    (&_S1886)->_S1185 = _S1885;
    float3  mean_c_11 = s_primal_ctx_mul_0(R_17, mean_12) + t_14;
    float4  intrins_11 = make_float4 (fx_20, fy_20, cx_15, cy_15);
    float3  _S1887 = s_primal_ctx_exp_0(scale_14);
    float _S1888 = quat_15.y;
    float x2_15 = _S1888 * _S1888;
    float y2_15 = quat_15.z * quat_15.z;
    float z2_26 = quat_15.w * quat_15.w;
    float xy_15 = quat_15.y * quat_15.z;
    float xz_15 = quat_15.y * quat_15.w;
    float yz_15 = quat_15.z * quat_15.w;
    float wx_15 = quat_15.x * quat_15.y;
    float wy_15 = quat_15.x * quat_15.z;
    float wz_15 = quat_15.x * quat_15.w;
    Matrix<float, 3, 3>  _S1889 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_15 + z2_26), 2.0f * (xy_15 + wz_15), 2.0f * (xz_15 - wy_15), 2.0f * (xy_15 - wz_15), 1.0f - 2.0f * (x2_15 + z2_26), 2.0f * (yz_15 + wx_15), 2.0f * (xz_15 + wy_15), 2.0f * (yz_15 - wx_15), 1.0f - 2.0f * (x2_15 + y2_15)));
    Matrix<float, 3, 3>  S_1 = makeMatrix<float, 3, 3> (_S1887.x, 0.0f, 0.0f, 0.0f, _S1887.y, 0.0f, 0.0f, 0.0f, _S1887.z);
    Matrix<float, 3, 3>  _S1890 = s_primal_ctx_mul_1(_S1889, S_1);
    Matrix<float, 3, 3>  _S1891 = transpose_0(_S1890);
    Matrix<float, 3, 3>  _S1892 = s_primal_ctx_mul_1(_S1890, _S1891);
    Matrix<float, 3, 3>  _S1893 = s_primal_ctx_mul_1(R_17, _S1892);
    Matrix<float, 3, 3>  _S1894 = transpose_0(R_17);
    Matrix<float, 3, 3>  _S1895 = s_primal_ctx_mul_1(_S1893, _S1894);
    Matrix<float, 2, 2>  _S1896 = _S1882;
    float2  _S1897 = make_float2 (0.0f);
    float2  _S1898 = _S1897;
    bool _S1899 = s_primal_ctx_fisheye_proj_3dgs_0(mean_c_11, _S1895, intrins_11, dist_coeffs_24, &_S1896, &_S1898, &(&_S1886)->_S1185);
    (&_S1886)->_S1184 = _S1896;
    s_bwd_prop_projection_3dgs_fisheye_differentiable_Intermediates_0 _S1900 = _S1886;
    float _S1901 = _S1886._S1184.rows[int(0)].y * _S1886._S1184.rows[int(1)].x;
    float det_orig_12 = _S1886._S1184.rows[int(0)].x * _S1886._S1184.rows[int(1)].y - _S1901;
    float _S1902 = _S1886._S1184.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S1903 = _S1886._S1184;
    *&(((&_S1903)->rows + (int(0)))->x) = _S1902;
    float _S1904 = _S1886._S1184.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S1903)->rows + (int(1)))->y) = _S1904;
    Matrix<float, 2, 2>  _S1905 = _S1903;
    Matrix<float, 2, 2>  _S1906 = _S1903;
    float det_blur_7 = _S1902 * _S1904 - _S1901;
    float _S1907 = det_orig_12 / det_blur_7;
    float _S1908 = det_blur_7 * det_blur_7;
    float _S1909 = s_primal_ctx_max_0(0.0f, _S1907);
    float _S1910 = s_primal_ctx_sqrt_0(_S1909);
    float invdet_8 = 1.0f / det_blur_7;
    float _S1911 = - _S1886._S1184.rows[int(0)].y;
    float _S1912 = - _S1886._S1184.rows[int(1)].x;
    float _S1913 = - in_opacity_11;
    float _S1914 = 1.0f + s_primal_ctx_exp_1(_S1913);
    float _S1915 = 1.0f / _S1914;
    float _S1916 = _S1914 * _S1914;
    float _S1917;
    if(antialiased_11)
    {
        _S1917 = _S1915 * _S1910;
    }
    else
    {
        _S1917 = _S1915;
    }
    float _S1918 = _S1917 / 0.00392156885936856f;
    float _S1919 = 2.0f * s_primal_ctx_log_0(_S1918);
    float _S1920 = s_primal_ctx_sqrt_0(_S1919);
    float _S1921 = _S1905.rows[int(0)].x;
    float _S1922 = _S1906.rows[int(1)].y;
    float _S1923 = s_primal_ctx_dot_0(mean_c_11, mean_c_11) + 9.99999997475242708e-07f;
    float3  _S1924 = mean_12 - - s_primal_ctx_mul_0(_S1894, t_14);
    float _S1925 = _S1924.x;
    float _S1926 = _S1924.y;
    float _S1927 = _S1924.z;
    float _S1928 = _S1925 * _S1925 + _S1926 * _S1926 + _S1927 * _S1927;
    float _S1929 = s_primal_ctx_sqrt_0(_S1928);
    float x_37 = _S1925 / _S1929;
    float3  _S1930 = make_float3 (x_37);
    float _S1931 = _S1929 * _S1929;
    float y_14 = _S1926 / _S1929;
    float z_11 = _S1927 / _S1929;
    float3  _S1932 = make_float3 (z_11);
    float _S1933 = - y_14;
    float3  _S1934 = make_float3 (_S1933);
    float z2_27 = z_11 * z_11;
    float fTmp0B_11 = -1.09254848957061768f * z_11;
    float fC1_11 = x_37 * x_37 - y_14 * y_14;
    float _S1935 = 2.0f * x_37;
    float fS1_11 = _S1935 * y_14;
    float pSH6_1 = 0.94617468118667603f * z2_27 - 0.31539157032966614f;
    float3  _S1936 = make_float3 (pSH6_1);
    float pSH7_1 = fTmp0B_11 * x_37;
    float3  _S1937 = make_float3 (pSH7_1);
    float pSH5_1 = fTmp0B_11 * y_14;
    float3  _S1938 = make_float3 (pSH5_1);
    float pSH8_1 = 0.54627424478530884f * fC1_11;
    float3  _S1939 = make_float3 (pSH8_1);
    float pSH4_1 = 0.54627424478530884f * fS1_11;
    float3  _S1940 = make_float3 (pSH4_1);
    float fTmp0C_11 = -2.28522896766662598f * z2_27 + 0.4570457935333252f;
    float fTmp1B_11 = 1.44530570507049561f * z_11;
    float _S1941 = 1.86588168144226074f * z2_27 - 1.11952900886535645f;
    float pSH12_1 = z_11 * _S1941;
    float3  _S1942 = make_float3 (pSH12_1);
    float pSH13_1 = fTmp0C_11 * x_37;
    float3  _S1943 = make_float3 (pSH13_1);
    float pSH11_1 = fTmp0C_11 * y_14;
    float3  _S1944 = make_float3 (pSH11_1);
    float pSH14_1 = fTmp1B_11 * fC1_11;
    float3  _S1945 = make_float3 (pSH14_1);
    float pSH10_1 = fTmp1B_11 * fS1_11;
    float3  _S1946 = make_float3 (pSH10_1);
    float pSH15_1 = -0.59004360437393188f * (x_37 * fC1_11 - y_14 * fS1_11);
    float3  _S1947 = make_float3 (pSH15_1);
    float pSH9_1 = -0.59004360437393188f * (x_37 * fS1_11 + y_14 * fC1_11);
    float3  _S1948 = make_float3 (pSH9_1);
    float3  _S1949 = make_float3 (0.0f);
    float3  _S1950 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1951;
    (&_S1951)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_11)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S1933) * (*sh_coeffs_11)[int(1)] + make_float3 (z_11) * (*sh_coeffs_11)[int(2)] - make_float3 (x_37) * (*sh_coeffs_11)[int(3)]) + (make_float3 (pSH4_1) * (*sh_coeffs_11)[int(4)] + make_float3 (pSH5_1) * (*sh_coeffs_11)[int(5)] + make_float3 (pSH6_1) * (*sh_coeffs_11)[int(6)] + make_float3 (pSH7_1) * (*sh_coeffs_11)[int(7)] + make_float3 (pSH8_1) * (*sh_coeffs_11)[int(8)]) + (make_float3 (pSH9_1) * (*sh_coeffs_11)[int(9)] + make_float3 (pSH10_1) * (*sh_coeffs_11)[int(10)] + make_float3 (pSH11_1) * (*sh_coeffs_11)[int(11)] + make_float3 (pSH12_1) * (*sh_coeffs_11)[int(12)] + make_float3 (pSH13_1) * (*sh_coeffs_11)[int(13)] + make_float3 (pSH14_1) * (*sh_coeffs_11)[int(14)] + make_float3 (pSH15_1) * (*sh_coeffs_11)[int(15)]) + make_float3 (0.5f);
    (&_S1951)->differential_0 = _S1950;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S1952;
    (&_S1952)->primal_0 = _S1949;
    (&_S1952)->differential_0 = _S1950;
    s_bwd_prop_max_0(&_S1951, &_S1952, v_rgb_1);
    float3  _S1953 = _S1947 * _S1951.differential_0;
    float3  _S1954 = (*sh_coeffs_11)[int(15)] * _S1951.differential_0;
    float3  _S1955 = _S1945 * _S1951.differential_0;
    float3  _S1956 = (*sh_coeffs_11)[int(14)] * _S1951.differential_0;
    float3  _S1957 = _S1943 * _S1951.differential_0;
    float3  _S1958 = (*sh_coeffs_11)[int(13)] * _S1951.differential_0;
    float3  _S1959 = _S1942 * _S1951.differential_0;
    float3  _S1960 = (*sh_coeffs_11)[int(12)] * _S1951.differential_0;
    float3  _S1961 = _S1944 * _S1951.differential_0;
    float3  _S1962 = (*sh_coeffs_11)[int(11)] * _S1951.differential_0;
    float3  _S1963 = _S1946 * _S1951.differential_0;
    float3  _S1964 = (*sh_coeffs_11)[int(10)] * _S1951.differential_0;
    float3  _S1965 = _S1948 * _S1951.differential_0;
    float3  _S1966 = (*sh_coeffs_11)[int(9)] * _S1951.differential_0;
    float s_diff_fS2_T_1 = -0.59004360437393188f * (_S1966.x + _S1966.y + _S1966.z);
    float s_diff_fC2_T_1 = -0.59004360437393188f * (_S1954.x + _S1954.y + _S1954.z);
    float _S1967 = _S1964.x + _S1964.y + _S1964.z;
    float _S1968 = _S1956.x + _S1956.y + _S1956.z;
    float _S1969 = _S1962.x + _S1962.y + _S1962.z;
    float _S1970 = _S1958.x + _S1958.y + _S1958.z;
    float _S1971 = _S1960.x + _S1960.y + _S1960.z;
    float _S1972 = - s_diff_fC2_T_1;
    float3  _S1973 = _S1939 * _S1951.differential_0;
    float3  _S1974 = (*sh_coeffs_11)[int(8)] * _S1951.differential_0;
    float3  _S1975 = _S1937 * _S1951.differential_0;
    float3  _S1976 = (*sh_coeffs_11)[int(7)] * _S1951.differential_0;
    float3  _S1977 = _S1936 * _S1951.differential_0;
    float3  _S1978 = (*sh_coeffs_11)[int(6)] * _S1951.differential_0;
    float3  _S1979 = _S1938 * _S1951.differential_0;
    float3  _S1980 = (*sh_coeffs_11)[int(5)] * _S1951.differential_0;
    float3  _S1981 = _S1940 * _S1951.differential_0;
    float3  _S1982 = (*sh_coeffs_11)[int(4)] * _S1951.differential_0;
    float _S1983 = _S1980.x + _S1980.y + _S1980.z;
    float _S1984 = _S1976.x + _S1976.y + _S1976.z;
    float _S1985 = fTmp1B_11 * _S1967 + x_37 * s_diff_fS2_T_1 + y_14 * _S1972 + 0.54627424478530884f * (_S1982.x + _S1982.y + _S1982.z);
    float _S1986 = fTmp1B_11 * _S1968 + y_14 * s_diff_fS2_T_1 + x_37 * s_diff_fC2_T_1 + 0.54627424478530884f * (_S1974.x + _S1974.y + _S1974.z);
    float _S1987 = y_14 * - _S1986;
    float _S1988 = x_37 * _S1986;
    float _S1989 = z_11 * (1.86588168144226074f * (z_11 * _S1971) + -2.28522896766662598f * (y_14 * _S1969 + x_37 * _S1970) + 0.94617468118667603f * (_S1978.x + _S1978.y + _S1978.z));
    float3  _S1990 = make_float3 (0.48860251903533936f) * _S1951.differential_0;
    float3  _S1991 = - _S1990;
    float3  _S1992 = _S1930 * _S1991;
    float3  _S1993 = (*sh_coeffs_11)[int(3)] * _S1991;
    float3  _S1994 = _S1932 * _S1990;
    float3  _S1995 = (*sh_coeffs_11)[int(2)] * _S1990;
    float3  _S1996 = _S1934 * _S1990;
    float3  _S1997 = (*sh_coeffs_11)[int(1)] * _S1990;
    float _S1998 = (_S1941 * _S1971 + 1.44530570507049561f * (fS1_11 * _S1967 + fC1_11 * _S1968) + -1.09254848957061768f * (y_14 * _S1983 + x_37 * _S1984) + _S1989 + _S1989 + _S1995.x + _S1995.y + _S1995.z) / _S1931;
    float _S1999 = _S1929 * _S1998;
    float _S2000 = (fTmp0C_11 * _S1969 + fC1_11 * s_diff_fS2_T_1 + fS1_11 * _S1972 + fTmp0B_11 * _S1983 + _S1935 * _S1985 + _S1987 + _S1987 + - (_S1997.x + _S1997.y + _S1997.z)) / _S1931;
    float _S2001 = _S1929 * _S2000;
    float _S2002 = (fTmp0C_11 * _S1970 + fS1_11 * s_diff_fS2_T_1 + fC1_11 * s_diff_fC2_T_1 + fTmp0B_11 * _S1984 + 2.0f * (y_14 * _S1985) + _S1988 + _S1988 + _S1993.x + _S1993.y + _S1993.z) / _S1931;
    float _S2003 = _S1929 * _S2002;
    float _S2004 = _S1927 * - _S1998 + _S1926 * - _S2000 + _S1925 * - _S2002;
    DiffPair_float_0 _S2005;
    (&_S2005)->primal_0 = _S1928;
    (&_S2005)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2005, _S2004);
    float _S2006 = _S1927 * _S2005.differential_0;
    float _S2007 = _S1926 * _S2005.differential_0;
    float _S2008 = _S1925 * _S2005.differential_0;
    float3  _S2009 = make_float3 (0.282094806432724f) * _S1951.differential_0;
    float3  _S2010 = make_float3 (_S2003 + _S2008 + _S2008, _S2001 + _S2007 + _S2007, _S1999 + _S2006 + _S2006);
    float3  _S2011 = - - _S2010;
    Matrix<float, 3, 3>  _S2012 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2013;
    (&_S2013)->primal_0 = _S1894;
    (&_S2013)->differential_0 = _S2012;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2014;
    (&_S2014)->primal_0 = t_14;
    (&_S2014)->differential_0 = _S1950;
    s_bwd_prop_mul_1(&_S2013, &_S2014, _S2011);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2015 = _S2013;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2016 = _S2014;
    float2  _S2017 = _S1897;
    *&((&_S2017)->y) = v_conic_1.z;
    float2  _S2018 = _S1897;
    *&((&_S2018)->y) = v_conic_1.y;
    *&((&_S2018)->x) = v_conic_1.x;
    float _S2019 = 0.5f * v_depth_1;
    DiffPair_float_0 _S2020;
    (&_S2020)->primal_0 = _S1923;
    (&_S2020)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2020, _S2019);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2021;
    (&_S2021)->primal_0 = mean_c_11;
    (&_S2021)->differential_0 = _S1950;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2022;
    (&_S2022)->primal_0 = mean_c_11;
    (&_S2022)->differential_0 = _S1950;
    s_bwd_prop_dot_0(&_S2021, &_S2022, _S2020.differential_0);
    DiffPair_float_0 _S2023;
    (&_S2023)->primal_0 = _S1922;
    (&_S2023)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2023, 0.0f);
    DiffPair_float_0 _S2024;
    (&_S2024)->primal_0 = _S1921;
    (&_S2024)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2024, 0.0f);
    DiffPair_float_0 _S2025;
    (&_S2025)->primal_0 = 3.32999992370605469f;
    (&_S2025)->differential_0 = 0.0f;
    DiffPair_float_0 _S2026;
    (&_S2026)->primal_0 = _S1920;
    (&_S2026)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2025, &_S2026, 0.0f);
    DiffPair_float_0 _S2027;
    (&_S2027)->primal_0 = _S1919;
    (&_S2027)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2027, _S2026.differential_0);
    float _S2028 = 2.0f * _S2027.differential_0;
    DiffPair_float_0 _S2029;
    (&_S2029)->primal_0 = _S1918;
    (&_S2029)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2029, _S2028);
    float2  _S2030 = make_float2 (_S2024.differential_0, 0.0f);
    float _S2031 = v_opacity_1 + 254.9999847412109375f * _S2029.differential_0;
    Matrix<float, 2, 2>  _S2032 = _S1882;
    _S2032[int(1)] = _S2017;
    _S2032[int(0)] = _S2018;
    Matrix<float, 2, 2>  _S2033 = _S2032;
    FixedArray<float3 , 16>  _S2034;
    _S2034[int(0)] = _S1950;
    _S2034[int(1)] = _S1950;
    _S2034[int(2)] = _S1950;
    _S2034[int(3)] = _S1950;
    _S2034[int(4)] = _S1950;
    _S2034[int(5)] = _S1950;
    _S2034[int(6)] = _S1950;
    _S2034[int(7)] = _S1950;
    _S2034[int(8)] = _S1950;
    _S2034[int(9)] = _S1950;
    _S2034[int(10)] = _S1950;
    _S2034[int(11)] = _S1950;
    _S2034[int(12)] = _S1950;
    _S2034[int(13)] = _S1950;
    _S2034[int(14)] = _S1950;
    _S2034[int(15)] = _S1950;
    _S2034[int(7)] = _S1975;
    _S2034[int(0)] = _S2009;
    _S2034[int(1)] = _S1996;
    _S2034[int(2)] = _S1994;
    _S2034[int(3)] = _S1992;
    _S2034[int(4)] = _S1981;
    _S2034[int(5)] = _S1979;
    _S2034[int(6)] = _S1977;
    _S2034[int(15)] = _S1953;
    _S2034[int(8)] = _S1973;
    _S2034[int(9)] = _S1965;
    _S2034[int(10)] = _S1963;
    _S2034[int(11)] = _S1961;
    _S2034[int(12)] = _S1959;
    _S2034[int(13)] = _S1957;
    _S2034[int(14)] = _S1955;
    float3  _S2035 = _S2034[int(0)];
    float3  _S2036 = _S2034[int(1)];
    float3  _S2037 = _S2034[int(2)];
    float3  _S2038 = _S2034[int(3)];
    float3  _S2039 = _S2034[int(4)];
    float3  _S2040 = _S2034[int(5)];
    float3  _S2041 = _S2034[int(6)];
    float3  _S2042 = _S2034[int(7)];
    float3  _S2043 = _S2034[int(8)];
    float3  _S2044 = _S2034[int(9)];
    float3  _S2045 = _S2034[int(10)];
    float3  _S2046 = _S2034[int(11)];
    float3  _S2047 = _S2034[int(12)];
    float3  _S2048 = _S2034[int(13)];
    float3  _S2049 = _S2034[int(14)];
    float3  _S2050 = _S2034[int(15)];
    float3  _S2051 = _S2022.differential_0 + _S2021.differential_0;
    float2  _S2052 = make_float2 (0.0f, _S2023.differential_0);
    float _S2053;
    if(antialiased_11)
    {
        float _S2054 = _S1915 * _S2031;
        _S1917 = _S1910 * _S2031;
        _S2053 = _S2054;
    }
    else
    {
        _S1917 = _S2031;
        _S2053 = 0.0f;
    }
    float _S2055 = - (_S1917 / _S1916);
    DiffPair_float_0 _S2056;
    (&_S2056)->primal_0 = _S1913;
    (&_S2056)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2056, _S2055);
    float _S2057 = - _S2056.differential_0;
    float _S2058 = invdet_8 * _S2033.rows[int(1)].y;
    float _S2059 = - (invdet_8 * _S2033.rows[int(1)].x);
    float _S2060 = - (invdet_8 * _S2033.rows[int(0)].y);
    float _S2061 = invdet_8 * _S2033.rows[int(0)].x;
    float _S2062 = - ((_S1902 * _S2033.rows[int(1)].y + _S1912 * _S2033.rows[int(1)].x + _S1911 * _S2033.rows[int(0)].y + _S1904 * _S2033.rows[int(0)].x) / _S1908);
    DiffPair_float_0 _S2063;
    (&_S2063)->primal_0 = _S1909;
    (&_S2063)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2063, _S2053);
    DiffPair_float_0 _S2064;
    (&_S2064)->primal_0 = 0.0f;
    (&_S2064)->differential_0 = 0.0f;
    DiffPair_float_0 _S2065;
    (&_S2065)->primal_0 = _S1907;
    (&_S2065)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2064, &_S2065, _S2063.differential_0);
    float _S2066 = _S2065.differential_0 / _S1908;
    float s_diff_det_orig_T_1 = det_blur_7 * _S2066;
    float _S2067 = _S2062 + det_orig_12 * - _S2066;
    float _S2068 = - _S2067;
    float _S2069 = _S1902 * _S2067;
    float _S2070 = _S1904 * _S2067;
    Matrix<float, 2, 2>  _S2071 = _S1882;
    _S2071[int(1)] = _S2052;
    _S2071[int(0)] = _S2030;
    _S1903 = _S2071;
    *&(((&_S1903)->rows + (int(1)))->y) = 0.0f;
    float _S2072 = _S2061 + _S2069 + _S2071.rows[int(1)].y;
    *&(((&_S1903)->rows + (int(0)))->x) = 0.0f;
    float _S2073 = _S2058 + _S2070 + _S2071.rows[int(0)].x;
    float _S2074 = _S2068 + - s_diff_det_orig_T_1;
    float _S2075 = _S2059 + _S1900._S1184.rows[int(0)].y * _S2074;
    float _S2076 = _S2060 + _S1900._S1184.rows[int(1)].x * _S2074;
    float _S2077 = _S1900._S1184.rows[int(1)].y * s_diff_det_orig_T_1;
    float _S2078 = _S2072 + _S1900._S1184.rows[int(0)].x * s_diff_det_orig_T_1;
    float2  _S2079 = _S1897;
    *&((&_S2079)->x) = _S2075;
    *&((&_S2079)->y) = _S2078;
    float _S2080 = _S2073 + _S2077;
    float2  _S2081 = _S1897;
    *&((&_S2081)->y) = _S2076;
    *&((&_S2081)->x) = _S2080;
    Matrix<float, 2, 2>  _S2082 = _S1882;
    _S2082[int(1)] = _S2079;
    _S2082[int(0)] = _S2081;
    Matrix<float, 2, 2>  _S2083 = _S1903 + _S2082;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2084;
    (&_S2084)->primal_0 = mean_c_11;
    (&_S2084)->differential_0 = _S1950;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2085;
    (&_S2085)->primal_0 = _S1895;
    (&_S2085)->differential_0 = _S2012;
    float4  _S2086 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 _S2087;
    (&_S2087)->primal_0 = intrins_11;
    (&_S2087)->differential_0 = _S2086;
    FixedArray<float, 10>  _S2088 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C10x3E_0 _S2089;
    (&_S2089)->primal_0 = *dist_coeffs_24;
    (&_S2089)->differential_0 = _S2088;
    s_bwd_prop_fisheye_proj_3dgs_Intermediates_0 _S2090 = _S1900._S1185;
    s_bwd_prop_fisheye_proj_3dgs_0(&_S2084, &_S2085, &_S2087, &_S2089, _S2083, v_mean2d_1, &_S2090);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2091;
    (&_S2091)->primal_0 = _S1893;
    (&_S2091)->differential_0 = _S2012;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2092;
    (&_S2092)->primal_0 = _S1894;
    (&_S2092)->differential_0 = _S2012;
    s_bwd_prop_mul_4(&_S2091, &_S2092, _S2085.differential_0);
    Matrix<float, 3, 3>  _S2093 = transpose_0(_S2092.differential_0 + _S2015.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2094;
    (&_S2094)->primal_0 = R_17;
    (&_S2094)->differential_0 = _S2012;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2095;
    (&_S2095)->primal_0 = _S1892;
    (&_S2095)->differential_0 = _S2012;
    s_bwd_prop_mul_4(&_S2094, &_S2095, _S2091.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2096;
    (&_S2096)->primal_0 = _S1890;
    (&_S2096)->differential_0 = _S2012;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2097;
    (&_S2097)->primal_0 = _S1891;
    (&_S2097)->differential_0 = _S2012;
    s_bwd_prop_mul_4(&_S2096, &_S2097, _S2095.differential_0);
    Matrix<float, 3, 3>  _S2098 = _S2096.differential_0 + transpose_0(_S2097.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2099;
    (&_S2099)->primal_0 = _S1889;
    (&_S2099)->differential_0 = _S2012;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2100;
    (&_S2100)->primal_0 = S_1;
    (&_S2100)->differential_0 = _S2012;
    s_bwd_prop_mul_4(&_S2099, &_S2100, _S2098);
    Matrix<float, 3, 3>  _S2101 = transpose_0(_S2099.differential_0);
    float _S2102 = 2.0f * - _S2101.rows[int(2)].z;
    float _S2103 = 2.0f * _S2101.rows[int(2)].y;
    float _S2104 = 2.0f * _S2101.rows[int(2)].x;
    float _S2105 = 2.0f * _S2101.rows[int(1)].z;
    float _S2106 = 2.0f * - _S2101.rows[int(1)].y;
    float _S2107 = 2.0f * _S2101.rows[int(1)].x;
    float _S2108 = 2.0f * _S2101.rows[int(0)].z;
    float _S2109 = 2.0f * _S2101.rows[int(0)].y;
    float _S2110 = 2.0f * - _S2101.rows[int(0)].x;
    float _S2111 = - _S2107 + _S2109;
    float _S2112 = _S2104 + - _S2108;
    float _S2113 = - _S2103 + _S2105;
    float _S2114 = _S2103 + _S2105;
    float _S2115 = _S2104 + _S2108;
    float _S2116 = _S2107 + _S2109;
    float _S2117 = quat_15.w * (_S2106 + _S2110);
    float _S2118 = quat_15.z * (_S2102 + _S2110);
    float _S2119 = quat_15.y * (_S2102 + _S2106);
    float _S2120 = quat_15.x * _S2111 + quat_15.z * _S2114 + quat_15.y * _S2115 + _S2117 + _S2117;
    float _S2121 = quat_15.x * _S2112 + quat_15.w * _S2114 + quat_15.y * _S2116 + _S2118 + _S2118;
    float _S2122 = quat_15.x * _S2113 + quat_15.w * _S2115 + quat_15.z * _S2116 + _S2119 + _S2119;
    float _S2123 = quat_15.w * _S2111 + quat_15.z * _S2112 + quat_15.y * _S2113;
    float3  _S2124 = _S1950;
    *&((&_S2124)->z) = _S2100.differential_0.rows[int(2)].z;
    *&((&_S2124)->y) = _S2100.differential_0.rows[int(1)].y;
    *&((&_S2124)->x) = _S2100.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2125;
    (&_S2125)->primal_0 = scale_14;
    (&_S2125)->differential_0 = _S1950;
    s_bwd_prop_exp_1(&_S2125, _S2124);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2126 = _S2125;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2127;
    (&_S2127)->primal_0 = mean_c_11;
    (&_S2127)->differential_0 = _S1950;
    s_bwd_length_impl_1(&_S2127, 0.0f);
    float3  _S2128 = _S2084.differential_0 + _S2127.differential_0 + _S2051;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2129;
    (&_S2129)->primal_0 = R_17;
    (&_S2129)->differential_0 = _S2012;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2130;
    (&_S2130)->primal_0 = mean_12;
    (&_S2130)->differential_0 = _S1950;
    s_bwd_prop_mul_1(&_S2129, &_S2130, _S2128);
    float3  _S2131 = _S2128 + _S2016.differential_0;
    Matrix<float, 3, 3>  _S2132 = _S2093 + _S2094.differential_0 + _S2129.differential_0;
    float4  _S2133 = _S2086;
    *&((&_S2133)->w) = _S2120;
    *&((&_S2133)->z) = _S2121;
    *&((&_S2133)->y) = _S2122;
    *&((&_S2133)->x) = _S2123;
    float4  _S2134 = _S2133;
    float3  _S2135 = _S2130.differential_0 + _S2010;
    *v_mean_1 = _S2135;
    *v_quat_1 = _S2134;
    *v_scale_1 = _S2126.differential_0;
    *v_in_opacity_1 = _S2057;
    (*v_sh_coeffs_1)[int(0)] = _S2035;
    (*v_sh_coeffs_1)[int(1)] = _S2036;
    (*v_sh_coeffs_1)[int(2)] = _S2037;
    (*v_sh_coeffs_1)[int(3)] = _S2038;
    (*v_sh_coeffs_1)[int(4)] = _S2039;
    (*v_sh_coeffs_1)[int(5)] = _S2040;
    (*v_sh_coeffs_1)[int(6)] = _S2041;
    (*v_sh_coeffs_1)[int(7)] = _S2042;
    (*v_sh_coeffs_1)[int(8)] = _S2043;
    (*v_sh_coeffs_1)[int(9)] = _S2044;
    (*v_sh_coeffs_1)[int(10)] = _S2045;
    (*v_sh_coeffs_1)[int(11)] = _S2046;
    (*v_sh_coeffs_1)[int(12)] = _S2047;
    (*v_sh_coeffs_1)[int(13)] = _S2048;
    (*v_sh_coeffs_1)[int(14)] = _S2049;
    (*v_sh_coeffs_1)[int(15)] = _S2050;
    *v_R_3 = _S2132;
    *v_t_2 = _S2131;
    return;
}

inline __device__ void projection_3dgs_ortho_vjp(bool antialiased_12, float3  mean_13, float4  quat_16, float3  scale_15, float in_opacity_12, FixedArray<float3 , 16>  * sh_coeffs_12, Matrix<float, 3, 3>  R_18, float3  t_15, float fx_21, float fy_21, float cx_16, float cy_16, FixedArray<float, 10>  * dist_coeffs_25, uint image_width_12, uint image_height_12, float2  v_mean2d_2, float v_depth_2, float3  v_conic_2, float v_opacity_2, float3  v_rgb_2, float3  * v_mean_2, float4  * v_quat_2, float3  * v_scale_2, float * v_in_opacity_2, FixedArray<float3 , 16>  * v_sh_coeffs_2, Matrix<float, 3, 3>  * v_R_4, float3  * v_t_3)
{
    float3  mean_c_12 = s_primal_ctx_mul_0(R_18, mean_13) + t_15;
    float3  _S2136 = s_primal_ctx_exp_0(scale_15);
    float _S2137 = quat_16.y;
    float x2_16 = _S2137 * _S2137;
    float y2_16 = quat_16.z * quat_16.z;
    float z2_28 = quat_16.w * quat_16.w;
    float xy_16 = quat_16.y * quat_16.z;
    float xz_16 = quat_16.y * quat_16.w;
    float yz_16 = quat_16.z * quat_16.w;
    float wx_16 = quat_16.x * quat_16.y;
    float wy_16 = quat_16.x * quat_16.z;
    float wz_16 = quat_16.x * quat_16.w;
    Matrix<float, 3, 3>  _S2138 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_16 + z2_28), 2.0f * (xy_16 + wz_16), 2.0f * (xz_16 - wy_16), 2.0f * (xy_16 - wz_16), 1.0f - 2.0f * (x2_16 + z2_28), 2.0f * (yz_16 + wx_16), 2.0f * (xz_16 + wy_16), 2.0f * (yz_16 - wx_16), 1.0f - 2.0f * (x2_16 + y2_16)));
    Matrix<float, 3, 3>  S_2 = makeMatrix<float, 3, 3> (_S2136.x, 0.0f, 0.0f, 0.0f, _S2136.y, 0.0f, 0.0f, 0.0f, _S2136.z);
    Matrix<float, 3, 3>  _S2139 = s_primal_ctx_mul_1(_S2138, S_2);
    Matrix<float, 3, 3>  _S2140 = transpose_0(_S2139);
    Matrix<float, 3, 3>  _S2141 = s_primal_ctx_mul_1(_S2139, _S2140);
    Matrix<float, 3, 3>  _S2142 = s_primal_ctx_mul_1(R_18, _S2141);
    Matrix<float, 3, 3>  _S2143 = transpose_0(R_18);
    Matrix<float, 3, 3>  _S2144 = s_primal_ctx_mul_1(_S2142, _S2143);
    Matrix<float, 2, 3>  J_11 = makeMatrix<float, 2, 3> (fx_21, 0.0f, 0.0f, 0.0f, fy_21, 0.0f);
    Matrix<float, 2, 3>  _S2145 = s_primal_ctx_mul_2(J_11, _S2144);
    Matrix<float, 3, 2>  _S2146 = transpose_1(J_11);
    Matrix<float, 2, 2>  _S2147 = s_primal_ctx_mul_3(_S2145, _S2146);
    float _S2148 = _S2147.rows[int(0)].y * _S2147.rows[int(1)].x;
    float det_orig_13 = _S2147.rows[int(0)].x * _S2147.rows[int(1)].y - _S2148;
    float _S2149 = _S2147.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2150 = _S2147;
    *&(((&_S2150)->rows + (int(0)))->x) = _S2149;
    float _S2151 = _S2147.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2150)->rows + (int(1)))->y) = _S2151;
    Matrix<float, 2, 2>  _S2152 = _S2150;
    Matrix<float, 2, 2>  _S2153 = _S2150;
    float det_blur_8 = _S2149 * _S2151 - _S2148;
    float _S2154 = det_orig_13 / det_blur_8;
    float _S2155 = det_blur_8 * det_blur_8;
    float _S2156 = s_primal_ctx_max_0(0.0f, _S2154);
    float _S2157 = s_primal_ctx_sqrt_0(_S2156);
    float invdet_9 = 1.0f / det_blur_8;
    float _S2158 = - _S2147.rows[int(0)].y;
    float _S2159 = - _S2147.rows[int(1)].x;
    float _S2160 = - in_opacity_12;
    float _S2161 = 1.0f + s_primal_ctx_exp_1(_S2160);
    float _S2162 = 1.0f / _S2161;
    float _S2163 = _S2161 * _S2161;
    float _S2164;
    if(antialiased_12)
    {
        _S2164 = _S2162 * _S2157;
    }
    else
    {
        _S2164 = _S2162;
    }
    float _S2165 = _S2164 / 0.00392156885936856f;
    float _S2166 = 2.0f * s_primal_ctx_log_0(_S2165);
    float _S2167 = s_primal_ctx_sqrt_0(_S2166);
    float _S2168 = _S2152.rows[int(0)].x;
    float _S2169 = _S2153.rows[int(1)].y;
    float _S2170 = s_primal_ctx_dot_0(mean_c_12, mean_c_12) + 9.99999997475242708e-07f;
    float3  _S2171 = mean_13 - - s_primal_ctx_mul_0(_S2143, t_15);
    float _S2172 = _S2171.x;
    float _S2173 = _S2171.y;
    float _S2174 = _S2171.z;
    float _S2175 = _S2172 * _S2172 + _S2173 * _S2173 + _S2174 * _S2174;
    float _S2176 = s_primal_ctx_sqrt_0(_S2175);
    float x_38 = _S2172 / _S2176;
    float3  _S2177 = make_float3 (x_38);
    float _S2178 = _S2176 * _S2176;
    float y_15 = _S2173 / _S2176;
    float z_12 = _S2174 / _S2176;
    float3  _S2179 = make_float3 (z_12);
    float _S2180 = - y_15;
    float3  _S2181 = make_float3 (_S2180);
    float z2_29 = z_12 * z_12;
    float fTmp0B_12 = -1.09254848957061768f * z_12;
    float fC1_12 = x_38 * x_38 - y_15 * y_15;
    float _S2182 = 2.0f * x_38;
    float fS1_12 = _S2182 * y_15;
    float pSH6_2 = 0.94617468118667603f * z2_29 - 0.31539157032966614f;
    float3  _S2183 = make_float3 (pSH6_2);
    float pSH7_2 = fTmp0B_12 * x_38;
    float3  _S2184 = make_float3 (pSH7_2);
    float pSH5_2 = fTmp0B_12 * y_15;
    float3  _S2185 = make_float3 (pSH5_2);
    float pSH8_2 = 0.54627424478530884f * fC1_12;
    float3  _S2186 = make_float3 (pSH8_2);
    float pSH4_2 = 0.54627424478530884f * fS1_12;
    float3  _S2187 = make_float3 (pSH4_2);
    float fTmp0C_12 = -2.28522896766662598f * z2_29 + 0.4570457935333252f;
    float fTmp1B_12 = 1.44530570507049561f * z_12;
    float _S2188 = 1.86588168144226074f * z2_29 - 1.11952900886535645f;
    float pSH12_2 = z_12 * _S2188;
    float3  _S2189 = make_float3 (pSH12_2);
    float pSH13_2 = fTmp0C_12 * x_38;
    float3  _S2190 = make_float3 (pSH13_2);
    float pSH11_2 = fTmp0C_12 * y_15;
    float3  _S2191 = make_float3 (pSH11_2);
    float pSH14_2 = fTmp1B_12 * fC1_12;
    float3  _S2192 = make_float3 (pSH14_2);
    float pSH10_2 = fTmp1B_12 * fS1_12;
    float3  _S2193 = make_float3 (pSH10_2);
    float pSH15_2 = -0.59004360437393188f * (x_38 * fC1_12 - y_15 * fS1_12);
    float3  _S2194 = make_float3 (pSH15_2);
    float pSH9_2 = -0.59004360437393188f * (x_38 * fS1_12 + y_15 * fC1_12);
    float3  _S2195 = make_float3 (pSH9_2);
    float3  _S2196 = make_float3 (0.0f);
    float3  _S2197 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2198;
    (&_S2198)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_12)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2180) * (*sh_coeffs_12)[int(1)] + make_float3 (z_12) * (*sh_coeffs_12)[int(2)] - make_float3 (x_38) * (*sh_coeffs_12)[int(3)]) + (make_float3 (pSH4_2) * (*sh_coeffs_12)[int(4)] + make_float3 (pSH5_2) * (*sh_coeffs_12)[int(5)] + make_float3 (pSH6_2) * (*sh_coeffs_12)[int(6)] + make_float3 (pSH7_2) * (*sh_coeffs_12)[int(7)] + make_float3 (pSH8_2) * (*sh_coeffs_12)[int(8)]) + (make_float3 (pSH9_2) * (*sh_coeffs_12)[int(9)] + make_float3 (pSH10_2) * (*sh_coeffs_12)[int(10)] + make_float3 (pSH11_2) * (*sh_coeffs_12)[int(11)] + make_float3 (pSH12_2) * (*sh_coeffs_12)[int(12)] + make_float3 (pSH13_2) * (*sh_coeffs_12)[int(13)] + make_float3 (pSH14_2) * (*sh_coeffs_12)[int(14)] + make_float3 (pSH15_2) * (*sh_coeffs_12)[int(15)]) + make_float3 (0.5f);
    (&_S2198)->differential_0 = _S2197;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2199;
    (&_S2199)->primal_0 = _S2196;
    (&_S2199)->differential_0 = _S2197;
    s_bwd_prop_max_0(&_S2198, &_S2199, v_rgb_2);
    float3  _S2200 = _S2194 * _S2198.differential_0;
    float3  _S2201 = (*sh_coeffs_12)[int(15)] * _S2198.differential_0;
    float3  _S2202 = _S2192 * _S2198.differential_0;
    float3  _S2203 = (*sh_coeffs_12)[int(14)] * _S2198.differential_0;
    float3  _S2204 = _S2190 * _S2198.differential_0;
    float3  _S2205 = (*sh_coeffs_12)[int(13)] * _S2198.differential_0;
    float3  _S2206 = _S2189 * _S2198.differential_0;
    float3  _S2207 = (*sh_coeffs_12)[int(12)] * _S2198.differential_0;
    float3  _S2208 = _S2191 * _S2198.differential_0;
    float3  _S2209 = (*sh_coeffs_12)[int(11)] * _S2198.differential_0;
    float3  _S2210 = _S2193 * _S2198.differential_0;
    float3  _S2211 = (*sh_coeffs_12)[int(10)] * _S2198.differential_0;
    float3  _S2212 = _S2195 * _S2198.differential_0;
    float3  _S2213 = (*sh_coeffs_12)[int(9)] * _S2198.differential_0;
    float s_diff_fS2_T_2 = -0.59004360437393188f * (_S2213.x + _S2213.y + _S2213.z);
    float s_diff_fC2_T_2 = -0.59004360437393188f * (_S2201.x + _S2201.y + _S2201.z);
    float _S2214 = _S2211.x + _S2211.y + _S2211.z;
    float _S2215 = _S2203.x + _S2203.y + _S2203.z;
    float _S2216 = _S2209.x + _S2209.y + _S2209.z;
    float _S2217 = _S2205.x + _S2205.y + _S2205.z;
    float _S2218 = _S2207.x + _S2207.y + _S2207.z;
    float _S2219 = - s_diff_fC2_T_2;
    float3  _S2220 = _S2186 * _S2198.differential_0;
    float3  _S2221 = (*sh_coeffs_12)[int(8)] * _S2198.differential_0;
    float3  _S2222 = _S2184 * _S2198.differential_0;
    float3  _S2223 = (*sh_coeffs_12)[int(7)] * _S2198.differential_0;
    float3  _S2224 = _S2183 * _S2198.differential_0;
    float3  _S2225 = (*sh_coeffs_12)[int(6)] * _S2198.differential_0;
    float3  _S2226 = _S2185 * _S2198.differential_0;
    float3  _S2227 = (*sh_coeffs_12)[int(5)] * _S2198.differential_0;
    float3  _S2228 = _S2187 * _S2198.differential_0;
    float3  _S2229 = (*sh_coeffs_12)[int(4)] * _S2198.differential_0;
    float _S2230 = _S2227.x + _S2227.y + _S2227.z;
    float _S2231 = _S2223.x + _S2223.y + _S2223.z;
    float _S2232 = fTmp1B_12 * _S2214 + x_38 * s_diff_fS2_T_2 + y_15 * _S2219 + 0.54627424478530884f * (_S2229.x + _S2229.y + _S2229.z);
    float _S2233 = fTmp1B_12 * _S2215 + y_15 * s_diff_fS2_T_2 + x_38 * s_diff_fC2_T_2 + 0.54627424478530884f * (_S2221.x + _S2221.y + _S2221.z);
    float _S2234 = y_15 * - _S2233;
    float _S2235 = x_38 * _S2233;
    float _S2236 = z_12 * (1.86588168144226074f * (z_12 * _S2218) + -2.28522896766662598f * (y_15 * _S2216 + x_38 * _S2217) + 0.94617468118667603f * (_S2225.x + _S2225.y + _S2225.z));
    float3  _S2237 = make_float3 (0.48860251903533936f) * _S2198.differential_0;
    float3  _S2238 = - _S2237;
    float3  _S2239 = _S2177 * _S2238;
    float3  _S2240 = (*sh_coeffs_12)[int(3)] * _S2238;
    float3  _S2241 = _S2179 * _S2237;
    float3  _S2242 = (*sh_coeffs_12)[int(2)] * _S2237;
    float3  _S2243 = _S2181 * _S2237;
    float3  _S2244 = (*sh_coeffs_12)[int(1)] * _S2237;
    float _S2245 = (_S2188 * _S2218 + 1.44530570507049561f * (fS1_12 * _S2214 + fC1_12 * _S2215) + -1.09254848957061768f * (y_15 * _S2230 + x_38 * _S2231) + _S2236 + _S2236 + _S2242.x + _S2242.y + _S2242.z) / _S2178;
    float _S2246 = _S2176 * _S2245;
    float _S2247 = (fTmp0C_12 * _S2216 + fC1_12 * s_diff_fS2_T_2 + fS1_12 * _S2219 + fTmp0B_12 * _S2230 + _S2182 * _S2232 + _S2234 + _S2234 + - (_S2244.x + _S2244.y + _S2244.z)) / _S2178;
    float _S2248 = _S2176 * _S2247;
    float _S2249 = (fTmp0C_12 * _S2217 + fS1_12 * s_diff_fS2_T_2 + fC1_12 * s_diff_fC2_T_2 + fTmp0B_12 * _S2231 + 2.0f * (y_15 * _S2232) + _S2235 + _S2235 + _S2240.x + _S2240.y + _S2240.z) / _S2178;
    float _S2250 = _S2176 * _S2249;
    float _S2251 = _S2174 * - _S2245 + _S2173 * - _S2247 + _S2172 * - _S2249;
    DiffPair_float_0 _S2252;
    (&_S2252)->primal_0 = _S2175;
    (&_S2252)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2252, _S2251);
    float _S2253 = _S2174 * _S2252.differential_0;
    float _S2254 = _S2173 * _S2252.differential_0;
    float _S2255 = _S2172 * _S2252.differential_0;
    float3  _S2256 = make_float3 (0.282094806432724f) * _S2198.differential_0;
    float3  _S2257 = make_float3 (_S2250 + _S2255 + _S2255, _S2248 + _S2254 + _S2254, _S2246 + _S2253 + _S2253);
    float3  _S2258 = - - _S2257;
    Matrix<float, 3, 3>  _S2259 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2260;
    (&_S2260)->primal_0 = _S2143;
    (&_S2260)->differential_0 = _S2259;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2261;
    (&_S2261)->primal_0 = t_15;
    (&_S2261)->differential_0 = _S2197;
    s_bwd_prop_mul_1(&_S2260, &_S2261, _S2258);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2262 = _S2260;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2263 = _S2261;
    float2  _S2264 = make_float2 (0.0f);
    float2  _S2265 = _S2264;
    *&((&_S2265)->y) = v_conic_2.z;
    float2  _S2266 = _S2264;
    *&((&_S2266)->y) = v_conic_2.y;
    *&((&_S2266)->x) = v_conic_2.x;
    float _S2267 = 0.5f * v_depth_2;
    DiffPair_float_0 _S2268;
    (&_S2268)->primal_0 = _S2170;
    (&_S2268)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2268, _S2267);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2269;
    (&_S2269)->primal_0 = mean_c_12;
    (&_S2269)->differential_0 = _S2197;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2270;
    (&_S2270)->primal_0 = mean_c_12;
    (&_S2270)->differential_0 = _S2197;
    s_bwd_prop_dot_0(&_S2269, &_S2270, _S2268.differential_0);
    DiffPair_float_0 _S2271;
    (&_S2271)->primal_0 = _S2169;
    (&_S2271)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2271, 0.0f);
    DiffPair_float_0 _S2272;
    (&_S2272)->primal_0 = _S2168;
    (&_S2272)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2272, 0.0f);
    DiffPair_float_0 _S2273;
    (&_S2273)->primal_0 = 3.32999992370605469f;
    (&_S2273)->differential_0 = 0.0f;
    DiffPair_float_0 _S2274;
    (&_S2274)->primal_0 = _S2167;
    (&_S2274)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2273, &_S2274, 0.0f);
    DiffPair_float_0 _S2275;
    (&_S2275)->primal_0 = _S2166;
    (&_S2275)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2275, _S2274.differential_0);
    float _S2276 = 2.0f * _S2275.differential_0;
    DiffPair_float_0 _S2277;
    (&_S2277)->primal_0 = _S2165;
    (&_S2277)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2277, _S2276);
    float _S2278 = v_opacity_2 + 254.9999847412109375f * _S2277.differential_0;
    Matrix<float, 2, 2>  _S2279 = makeMatrix<float, 2, 2> (0.0f);
    Matrix<float, 2, 2>  _S2280 = _S2279;
    _S2280[int(1)] = _S2265;
    _S2280[int(0)] = _S2266;
    Matrix<float, 2, 2>  _S2281 = _S2280;
    FixedArray<float3 , 16>  _S2282;
    _S2282[int(0)] = _S2197;
    _S2282[int(1)] = _S2197;
    _S2282[int(2)] = _S2197;
    _S2282[int(3)] = _S2197;
    _S2282[int(4)] = _S2197;
    _S2282[int(5)] = _S2197;
    _S2282[int(6)] = _S2197;
    _S2282[int(7)] = _S2197;
    _S2282[int(8)] = _S2197;
    _S2282[int(9)] = _S2197;
    _S2282[int(10)] = _S2197;
    _S2282[int(11)] = _S2197;
    _S2282[int(12)] = _S2197;
    _S2282[int(13)] = _S2197;
    _S2282[int(14)] = _S2197;
    _S2282[int(15)] = _S2197;
    _S2282[int(7)] = _S2222;
    _S2282[int(0)] = _S2256;
    _S2282[int(1)] = _S2243;
    _S2282[int(2)] = _S2241;
    _S2282[int(3)] = _S2239;
    _S2282[int(4)] = _S2228;
    _S2282[int(5)] = _S2226;
    _S2282[int(6)] = _S2224;
    _S2282[int(15)] = _S2200;
    _S2282[int(8)] = _S2220;
    _S2282[int(9)] = _S2212;
    _S2282[int(10)] = _S2210;
    _S2282[int(11)] = _S2208;
    _S2282[int(12)] = _S2206;
    _S2282[int(13)] = _S2204;
    _S2282[int(14)] = _S2202;
    float3  _S2283 = _S2282[int(0)];
    float3  _S2284 = _S2282[int(1)];
    float3  _S2285 = _S2282[int(2)];
    float3  _S2286 = _S2282[int(3)];
    float3  _S2287 = _S2282[int(4)];
    float3  _S2288 = _S2282[int(5)];
    float3  _S2289 = _S2282[int(6)];
    float3  _S2290 = _S2282[int(7)];
    float3  _S2291 = _S2282[int(8)];
    float3  _S2292 = _S2282[int(9)];
    float3  _S2293 = _S2282[int(10)];
    float3  _S2294 = _S2282[int(11)];
    float3  _S2295 = _S2282[int(12)];
    float3  _S2296 = _S2282[int(13)];
    float3  _S2297 = _S2282[int(14)];
    float3  _S2298 = _S2282[int(15)];
    float3  _S2299 = _S2270.differential_0 + _S2269.differential_0;
    float2  _S2300 = make_float2 (0.0f, _S2271.differential_0);
    float2  _S2301 = make_float2 (_S2272.differential_0, 0.0f);
    float _S2302;
    if(antialiased_12)
    {
        float _S2303 = _S2162 * _S2278;
        _S2164 = _S2157 * _S2278;
        _S2302 = _S2303;
    }
    else
    {
        _S2164 = _S2278;
        _S2302 = 0.0f;
    }
    float _S2304 = - (_S2164 / _S2163);
    DiffPair_float_0 _S2305;
    (&_S2305)->primal_0 = _S2160;
    (&_S2305)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2305, _S2304);
    float _S2306 = - _S2305.differential_0;
    float _S2307 = invdet_9 * _S2281.rows[int(1)].y;
    float _S2308 = - (invdet_9 * _S2281.rows[int(1)].x);
    float _S2309 = - (invdet_9 * _S2281.rows[int(0)].y);
    float _S2310 = invdet_9 * _S2281.rows[int(0)].x;
    float _S2311 = - ((_S2149 * _S2281.rows[int(1)].y + _S2159 * _S2281.rows[int(1)].x + _S2158 * _S2281.rows[int(0)].y + _S2151 * _S2281.rows[int(0)].x) / _S2155);
    DiffPair_float_0 _S2312;
    (&_S2312)->primal_0 = _S2156;
    (&_S2312)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2312, _S2302);
    DiffPair_float_0 _S2313;
    (&_S2313)->primal_0 = 0.0f;
    (&_S2313)->differential_0 = 0.0f;
    DiffPair_float_0 _S2314;
    (&_S2314)->primal_0 = _S2154;
    (&_S2314)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2313, &_S2314, _S2312.differential_0);
    float _S2315 = _S2314.differential_0 / _S2155;
    float s_diff_det_orig_T_2 = det_blur_8 * _S2315;
    float _S2316 = _S2311 + det_orig_13 * - _S2315;
    float _S2317 = - _S2316;
    float _S2318 = _S2149 * _S2316;
    float _S2319 = _S2151 * _S2316;
    Matrix<float, 2, 2>  _S2320 = _S2279;
    _S2320[int(1)] = _S2300;
    _S2320[int(0)] = _S2301;
    _S2150 = _S2320;
    *&(((&_S2150)->rows + (int(1)))->y) = 0.0f;
    float _S2321 = _S2310 + _S2318 + _S2320.rows[int(1)].y;
    *&(((&_S2150)->rows + (int(0)))->x) = 0.0f;
    float _S2322 = _S2307 + _S2319 + _S2320.rows[int(0)].x;
    float _S2323 = _S2317 + - s_diff_det_orig_T_2;
    float _S2324 = _S2308 + _S2147.rows[int(0)].y * _S2323;
    float _S2325 = _S2309 + _S2147.rows[int(1)].x * _S2323;
    float _S2326 = _S2147.rows[int(1)].y * s_diff_det_orig_T_2;
    float _S2327 = _S2321 + _S2147.rows[int(0)].x * s_diff_det_orig_T_2;
    float2  _S2328 = _S2264;
    *&((&_S2328)->x) = _S2324;
    *&((&_S2328)->y) = _S2327;
    float _S2329 = _S2322 + _S2326;
    float2  _S2330 = _S2264;
    *&((&_S2330)->y) = _S2325;
    *&((&_S2330)->x) = _S2329;
    float _S2331 = fy_21 * v_mean2d_2.y;
    float _S2332 = fx_21 * v_mean2d_2.x;
    Matrix<float, 2, 2>  _S2333 = _S2279;
    _S2333[int(1)] = _S2328;
    _S2333[int(0)] = _S2330;
    Matrix<float, 2, 2>  _S2334 = _S2150 + _S2333;
    Matrix<float, 2, 3>  _S2335 = makeMatrix<float, 2, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2336;
    (&_S2336)->primal_0 = _S2145;
    (&_S2336)->differential_0 = _S2335;
    Matrix<float, 3, 2>  _S2337 = makeMatrix<float, 3, 2> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C2x3E_0 _S2338;
    (&_S2338)->primal_0 = _S2146;
    (&_S2338)->differential_0 = _S2337;
    s_bwd_prop_mul_2(&_S2336, &_S2338, _S2334);
    DiffPair_matrixx3Cfloatx2C2x2C3x3E_0 _S2339;
    (&_S2339)->primal_0 = J_11;
    (&_S2339)->differential_0 = _S2335;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2340;
    (&_S2340)->primal_0 = _S2144;
    (&_S2340)->differential_0 = _S2259;
    s_bwd_prop_mul_3(&_S2339, &_S2340, _S2336.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2341;
    (&_S2341)->primal_0 = _S2142;
    (&_S2341)->differential_0 = _S2259;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2342;
    (&_S2342)->primal_0 = _S2143;
    (&_S2342)->differential_0 = _S2259;
    s_bwd_prop_mul_4(&_S2341, &_S2342, _S2340.differential_0);
    Matrix<float, 3, 3>  _S2343 = transpose_0(_S2342.differential_0 + _S2262.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2344;
    (&_S2344)->primal_0 = R_18;
    (&_S2344)->differential_0 = _S2259;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2345;
    (&_S2345)->primal_0 = _S2141;
    (&_S2345)->differential_0 = _S2259;
    s_bwd_prop_mul_4(&_S2344, &_S2345, _S2341.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2346;
    (&_S2346)->primal_0 = _S2139;
    (&_S2346)->differential_0 = _S2259;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2347;
    (&_S2347)->primal_0 = _S2140;
    (&_S2347)->differential_0 = _S2259;
    s_bwd_prop_mul_4(&_S2346, &_S2347, _S2345.differential_0);
    Matrix<float, 3, 3>  _S2348 = _S2346.differential_0 + transpose_0(_S2347.differential_0);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2349;
    (&_S2349)->primal_0 = _S2138;
    (&_S2349)->differential_0 = _S2259;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2350;
    (&_S2350)->primal_0 = S_2;
    (&_S2350)->differential_0 = _S2259;
    s_bwd_prop_mul_4(&_S2349, &_S2350, _S2348);
    Matrix<float, 3, 3>  _S2351 = transpose_0(_S2349.differential_0);
    float _S2352 = 2.0f * - _S2351.rows[int(2)].z;
    float _S2353 = 2.0f * _S2351.rows[int(2)].y;
    float _S2354 = 2.0f * _S2351.rows[int(2)].x;
    float _S2355 = 2.0f * _S2351.rows[int(1)].z;
    float _S2356 = 2.0f * - _S2351.rows[int(1)].y;
    float _S2357 = 2.0f * _S2351.rows[int(1)].x;
    float _S2358 = 2.0f * _S2351.rows[int(0)].z;
    float _S2359 = 2.0f * _S2351.rows[int(0)].y;
    float _S2360 = 2.0f * - _S2351.rows[int(0)].x;
    float _S2361 = - _S2357 + _S2359;
    float _S2362 = _S2354 + - _S2358;
    float _S2363 = - _S2353 + _S2355;
    float _S2364 = _S2353 + _S2355;
    float _S2365 = _S2354 + _S2358;
    float _S2366 = _S2357 + _S2359;
    float _S2367 = quat_16.w * (_S2356 + _S2360);
    float _S2368 = quat_16.z * (_S2352 + _S2360);
    float _S2369 = quat_16.y * (_S2352 + _S2356);
    float _S2370 = quat_16.x * _S2361 + quat_16.z * _S2364 + quat_16.y * _S2365 + _S2367 + _S2367;
    float _S2371 = quat_16.x * _S2362 + quat_16.w * _S2364 + quat_16.y * _S2366 + _S2368 + _S2368;
    float _S2372 = quat_16.x * _S2363 + quat_16.w * _S2365 + quat_16.z * _S2366 + _S2369 + _S2369;
    float _S2373 = quat_16.w * _S2361 + quat_16.z * _S2362 + quat_16.y * _S2363;
    float3  _S2374 = _S2197;
    *&((&_S2374)->z) = _S2350.differential_0.rows[int(2)].z;
    *&((&_S2374)->y) = _S2350.differential_0.rows[int(1)].y;
    *&((&_S2374)->x) = _S2350.differential_0.rows[int(0)].x;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2375;
    (&_S2375)->primal_0 = scale_15;
    (&_S2375)->differential_0 = _S2197;
    s_bwd_prop_exp_1(&_S2375, _S2374);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2376 = _S2375;
    float3  _S2377 = _S2197;
    *&((&_S2377)->y) = _S2331;
    *&((&_S2377)->x) = _S2332;
    float3  _S2378 = _S2299 + _S2377;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2379;
    (&_S2379)->primal_0 = R_18;
    (&_S2379)->differential_0 = _S2259;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2380;
    (&_S2380)->primal_0 = mean_13;
    (&_S2380)->differential_0 = _S2197;
    s_bwd_prop_mul_1(&_S2379, &_S2380, _S2378);
    float3  _S2381 = _S2378 + _S2263.differential_0;
    Matrix<float, 3, 3>  _S2382 = _S2343 + _S2344.differential_0 + _S2379.differential_0;
    float4  _S2383 = make_float4 (0.0f);
    *&((&_S2383)->w) = _S2370;
    *&((&_S2383)->z) = _S2371;
    *&((&_S2383)->y) = _S2372;
    *&((&_S2383)->x) = _S2373;
    float4  _S2384 = _S2383;
    float3  _S2385 = _S2380.differential_0 + _S2257;
    *v_mean_2 = _S2385;
    *v_quat_2 = _S2384;
    *v_scale_2 = _S2376.differential_0;
    *v_in_opacity_2 = _S2306;
    (*v_sh_coeffs_2)[int(0)] = _S2283;
    (*v_sh_coeffs_2)[int(1)] = _S2284;
    (*v_sh_coeffs_2)[int(2)] = _S2285;
    (*v_sh_coeffs_2)[int(3)] = _S2286;
    (*v_sh_coeffs_2)[int(4)] = _S2287;
    (*v_sh_coeffs_2)[int(5)] = _S2288;
    (*v_sh_coeffs_2)[int(6)] = _S2289;
    (*v_sh_coeffs_2)[int(7)] = _S2290;
    (*v_sh_coeffs_2)[int(8)] = _S2291;
    (*v_sh_coeffs_2)[int(9)] = _S2292;
    (*v_sh_coeffs_2)[int(10)] = _S2293;
    (*v_sh_coeffs_2)[int(11)] = _S2294;
    (*v_sh_coeffs_2)[int(12)] = _S2295;
    (*v_sh_coeffs_2)[int(13)] = _S2296;
    (*v_sh_coeffs_2)[int(14)] = _S2297;
    (*v_sh_coeffs_2)[int(15)] = _S2298;
    *v_R_4 = _S2382;
    *v_t_3 = _S2381;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2386;
};

inline __device__ void projection_3dgs_eval3d_persp_vjp(bool antialiased_13, float3  mean_14, float4  quat_17, float3  scale_16, float in_opacity_13, FixedArray<float3 , 16>  * sh_coeffs_13, Matrix<float, 3, 3>  R_19, float3  t_16, float fx_22, float fy_22, float cx_17, float cy_17, FixedArray<float, 10>  * dist_coeffs_26, uint image_width_13, uint image_height_13, float2  v_mean2d_3, float v_depth_3, float3  v_conic_3, float v_opacity_3, float3  v_rgb_3, float3  * v_mean_3, float4  * v_quat_3, float3  * v_scale_3, float * v_in_opacity_3, FixedArray<float3 , 16>  * v_sh_coeffs_3, Matrix<float, 3, 3>  * v_R_5, float3  * v_t_4)
{
    Matrix<float, 2, 2>  _S2387 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2388;
    (&_S2388)->_S2386 = _S2387;
    float2  _S2389 = make_float2 (0.0f);
    float3  _S2390 = make_float3 (0.0f);
    float4  intrins_12 = make_float4 (fx_22, fy_22, cx_17, cy_17);
    float3  _S2391 = s_primal_ctx_exp_0(scale_16);
    float _S2392 = quat_17.y;
    float x2_17 = _S2392 * _S2392;
    float y2_17 = quat_17.z * quat_17.z;
    float z2_30 = quat_17.w * quat_17.w;
    float xy_17 = quat_17.y * quat_17.z;
    float xz_17 = quat_17.y * quat_17.w;
    float yz_17 = quat_17.z * quat_17.w;
    float wx_17 = quat_17.x * quat_17.y;
    float wy_17 = quat_17.x * quat_17.z;
    float wz_17 = quat_17.x * quat_17.w;
    Matrix<float, 3, 3>  _S2393 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_17 + z2_30), 2.0f * (xy_17 + wz_17), 2.0f * (xz_17 - wy_17), 2.0f * (xy_17 - wz_17), 1.0f - 2.0f * (x2_17 + z2_30), 2.0f * (yz_17 + wx_17), 2.0f * (xz_17 + wy_17), 2.0f * (yz_17 - wx_17), 1.0f - 2.0f * (x2_17 + y2_17))));
    FixedArray<float3 , 7>  _S2394 = {
        _S2390, _S2390, _S2390, _S2390, _S2390, _S2390, _S2390
    };
    FixedArray<float, 7>  _S2395 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2396;
    (&_S2396)->p_0 = _S2394;
    (&_S2396)->w_mean_0 = _S2395;
    (&_S2396)->w_cov_0 = _S2395;
    (&_S2396)->p_0[int(0)] = mean_14;
    SigmaPoints_0 _S2397 = _S2396;
    (&_S2397)->w_mean_0[int(0)] = 0.0f;
    (&_S2397)->w_cov_0[int(0)] = 2.0f;
    float _S2398 = s_primal_ctx_sqrt_0(3.0f);
    float _S2399 = _S2398 * _S2391.x;
    float3  delta_12 = make_float3 (_S2399) * _S2393.rows[0U];
    float3  _S2400 = mean_14 + delta_12;
    float3  _S2401 = mean_14 - delta_12;
    float _S2402 = _S2398 * _S2391.y;
    float3  delta_13 = make_float3 (_S2402) * _S2393.rows[1U];
    float3  _S2403 = mean_14 + delta_13;
    float3  _S2404 = mean_14 - delta_13;
    float _S2405 = _S2398 * _S2391.z;
    float3  delta_14 = make_float3 (_S2405) * _S2393.rows[2U];
    float3  _S2406 = mean_14 + delta_14;
    float3  _S2407 = mean_14 - delta_14;
    (&_S2397)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2397)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2397)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2397)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2397)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2397)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2397)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2397)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2397)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2397)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2397)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2397)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2408 = _S2396;
    (&_S2397)->p_0[0U] = s_primal_ctx_mul_0(R_19, _S2396.p_0[0U]) + t_16;
    (&_S2397)->p_0[1U] = s_primal_ctx_mul_0(R_19, _S2400) + t_16;
    (&_S2397)->p_0[2U] = s_primal_ctx_mul_0(R_19, _S2403) + t_16;
    (&_S2397)->p_0[3U] = s_primal_ctx_mul_0(R_19, _S2406) + t_16;
    (&_S2397)->p_0[4U] = s_primal_ctx_mul_0(R_19, _S2401) + t_16;
    (&_S2397)->p_0[5U] = s_primal_ctx_mul_0(R_19, _S2404) + t_16;
    (&_S2397)->p_0[6U] = s_primal_ctx_mul_0(R_19, _S2407) + t_16;
    float2  _S2409 = _S2389;
    Matrix<float, 2, 2>  _S2410 = _S2387;
    SigmaPoints_0 _S2411 = _S2397;
    bool _S2412 = persp_proj_3dgs_ut_1(&_S2411, intrins_12, dist_coeffs_26, image_width_13, image_height_13, &_S2410, &_S2409);
    (&_S2388)->_S2386 = _S2410;
    s_bwd_prop_projection_3dgs_eval3d_persp_differentiable_Intermediates_0 _S2413 = _S2388;
    float3  mean_c_13 = s_primal_ctx_mul_0(R_19, mean_14) + t_16;
    float3  _S2414 = make_float3 (_S2399);
    float3  _S2415 = make_float3 (_S2402);
    float3  _S2416 = make_float3 (_S2405);
    float _S2417 = _S2388._S2386.rows[int(0)].y * _S2388._S2386.rows[int(1)].x;
    float det_orig_14 = _S2388._S2386.rows[int(0)].x * _S2388._S2386.rows[int(1)].y - _S2417;
    float _S2418 = _S2388._S2386.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2419 = _S2388._S2386;
    *&(((&_S2419)->rows + (int(0)))->x) = _S2418;
    float _S2420 = _S2388._S2386.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2419)->rows + (int(1)))->y) = _S2420;
    Matrix<float, 2, 2>  _S2421 = _S2419;
    Matrix<float, 2, 2>  _S2422 = _S2419;
    float det_blur_9 = _S2418 * _S2420 - _S2417;
    float _S2423 = det_orig_14 / det_blur_9;
    float _S2424 = det_blur_9 * det_blur_9;
    float _S2425 = s_primal_ctx_max_0(0.0f, _S2423);
    float _S2426 = s_primal_ctx_sqrt_0(_S2425);
    float _S2427 = - in_opacity_13;
    float _S2428 = 1.0f + s_primal_ctx_exp_1(_S2427);
    float _S2429 = 1.0f / _S2428;
    float _S2430 = _S2428 * _S2428;
    float _S2431;
    if(antialiased_13)
    {
        _S2431 = _S2429 * _S2426;
    }
    else
    {
        _S2431 = _S2429;
    }
    float _S2432 = _S2431 / 0.00392156885936856f;
    float _S2433 = 2.0f * s_primal_ctx_log_0(_S2432);
    float _S2434 = s_primal_ctx_sqrt_0(_S2433);
    float _S2435 = _S2421.rows[int(0)].x;
    float _S2436 = _S2422.rows[int(1)].y;
    float _S2437 = s_primal_ctx_dot_0(mean_c_13, mean_c_13) + 9.99999997475242708e-07f;
    float3  _S2438 = - scale_16;
    Matrix<float, 3, 3>  _S2439 = transpose_0(R_19);
    float3  _S2440 = mean_14 - - s_primal_ctx_mul_0(_S2439, t_16);
    float _S2441 = _S2440.x;
    float _S2442 = _S2440.y;
    float _S2443 = _S2440.z;
    float _S2444 = _S2441 * _S2441 + _S2442 * _S2442 + _S2443 * _S2443;
    float _S2445 = s_primal_ctx_sqrt_0(_S2444);
    float x_39 = _S2441 / _S2445;
    float3  _S2446 = make_float3 (x_39);
    float _S2447 = _S2445 * _S2445;
    float y_16 = _S2442 / _S2445;
    float z_13 = _S2443 / _S2445;
    float3  _S2448 = make_float3 (z_13);
    float _S2449 = - y_16;
    float3  _S2450 = make_float3 (_S2449);
    float z2_31 = z_13 * z_13;
    float fTmp0B_13 = -1.09254848957061768f * z_13;
    float fC1_13 = x_39 * x_39 - y_16 * y_16;
    float _S2451 = 2.0f * x_39;
    float fS1_13 = _S2451 * y_16;
    float pSH6_3 = 0.94617468118667603f * z2_31 - 0.31539157032966614f;
    float3  _S2452 = make_float3 (pSH6_3);
    float pSH7_3 = fTmp0B_13 * x_39;
    float3  _S2453 = make_float3 (pSH7_3);
    float pSH5_3 = fTmp0B_13 * y_16;
    float3  _S2454 = make_float3 (pSH5_3);
    float pSH8_3 = 0.54627424478530884f * fC1_13;
    float3  _S2455 = make_float3 (pSH8_3);
    float pSH4_3 = 0.54627424478530884f * fS1_13;
    float3  _S2456 = make_float3 (pSH4_3);
    float fTmp0C_13 = -2.28522896766662598f * z2_31 + 0.4570457935333252f;
    float fTmp1B_13 = 1.44530570507049561f * z_13;
    float _S2457 = 1.86588168144226074f * z2_31 - 1.11952900886535645f;
    float pSH12_3 = z_13 * _S2457;
    float3  _S2458 = make_float3 (pSH12_3);
    float pSH13_3 = fTmp0C_13 * x_39;
    float3  _S2459 = make_float3 (pSH13_3);
    float pSH11_3 = fTmp0C_13 * y_16;
    float3  _S2460 = make_float3 (pSH11_3);
    float pSH14_3 = fTmp1B_13 * fC1_13;
    float3  _S2461 = make_float3 (pSH14_3);
    float pSH10_3 = fTmp1B_13 * fS1_13;
    float3  _S2462 = make_float3 (pSH10_3);
    float pSH15_3 = -0.59004360437393188f * (x_39 * fC1_13 - y_16 * fS1_13);
    float3  _S2463 = make_float3 (pSH15_3);
    float pSH9_3 = -0.59004360437393188f * (x_39 * fS1_13 + y_16 * fC1_13);
    float3  _S2464 = make_float3 (pSH9_3);
    float3  _S2465 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2466;
    (&_S2466)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_13)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2449) * (*sh_coeffs_13)[int(1)] + make_float3 (z_13) * (*sh_coeffs_13)[int(2)] - make_float3 (x_39) * (*sh_coeffs_13)[int(3)]) + (make_float3 (pSH4_3) * (*sh_coeffs_13)[int(4)] + make_float3 (pSH5_3) * (*sh_coeffs_13)[int(5)] + make_float3 (pSH6_3) * (*sh_coeffs_13)[int(6)] + make_float3 (pSH7_3) * (*sh_coeffs_13)[int(7)] + make_float3 (pSH8_3) * (*sh_coeffs_13)[int(8)]) + (make_float3 (pSH9_3) * (*sh_coeffs_13)[int(9)] + make_float3 (pSH10_3) * (*sh_coeffs_13)[int(10)] + make_float3 (pSH11_3) * (*sh_coeffs_13)[int(11)] + make_float3 (pSH12_3) * (*sh_coeffs_13)[int(12)] + make_float3 (pSH13_3) * (*sh_coeffs_13)[int(13)] + make_float3 (pSH14_3) * (*sh_coeffs_13)[int(14)] + make_float3 (pSH15_3) * (*sh_coeffs_13)[int(15)]) + make_float3 (0.5f);
    (&_S2466)->differential_0 = _S2390;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2467;
    (&_S2467)->primal_0 = _S2465;
    (&_S2467)->differential_0 = _S2390;
    s_bwd_prop_max_0(&_S2466, &_S2467, v_rgb_3);
    float3  _S2468 = _S2463 * _S2466.differential_0;
    float3  _S2469 = (*sh_coeffs_13)[int(15)] * _S2466.differential_0;
    float3  _S2470 = _S2461 * _S2466.differential_0;
    float3  _S2471 = (*sh_coeffs_13)[int(14)] * _S2466.differential_0;
    float3  _S2472 = _S2459 * _S2466.differential_0;
    float3  _S2473 = (*sh_coeffs_13)[int(13)] * _S2466.differential_0;
    float3  _S2474 = _S2458 * _S2466.differential_0;
    float3  _S2475 = (*sh_coeffs_13)[int(12)] * _S2466.differential_0;
    float3  _S2476 = _S2460 * _S2466.differential_0;
    float3  _S2477 = (*sh_coeffs_13)[int(11)] * _S2466.differential_0;
    float3  _S2478 = _S2462 * _S2466.differential_0;
    float3  _S2479 = (*sh_coeffs_13)[int(10)] * _S2466.differential_0;
    float3  _S2480 = _S2464 * _S2466.differential_0;
    float3  _S2481 = (*sh_coeffs_13)[int(9)] * _S2466.differential_0;
    float s_diff_fS2_T_3 = -0.59004360437393188f * (_S2481.x + _S2481.y + _S2481.z);
    float s_diff_fC2_T_3 = -0.59004360437393188f * (_S2469.x + _S2469.y + _S2469.z);
    float _S2482 = _S2479.x + _S2479.y + _S2479.z;
    float _S2483 = _S2471.x + _S2471.y + _S2471.z;
    float _S2484 = _S2477.x + _S2477.y + _S2477.z;
    float _S2485 = _S2473.x + _S2473.y + _S2473.z;
    float _S2486 = _S2475.x + _S2475.y + _S2475.z;
    float _S2487 = - s_diff_fC2_T_3;
    float3  _S2488 = _S2455 * _S2466.differential_0;
    float3  _S2489 = (*sh_coeffs_13)[int(8)] * _S2466.differential_0;
    float3  _S2490 = _S2453 * _S2466.differential_0;
    float3  _S2491 = (*sh_coeffs_13)[int(7)] * _S2466.differential_0;
    float3  _S2492 = _S2452 * _S2466.differential_0;
    float3  _S2493 = (*sh_coeffs_13)[int(6)] * _S2466.differential_0;
    float3  _S2494 = _S2454 * _S2466.differential_0;
    float3  _S2495 = (*sh_coeffs_13)[int(5)] * _S2466.differential_0;
    float3  _S2496 = _S2456 * _S2466.differential_0;
    float3  _S2497 = (*sh_coeffs_13)[int(4)] * _S2466.differential_0;
    float _S2498 = _S2495.x + _S2495.y + _S2495.z;
    float _S2499 = _S2491.x + _S2491.y + _S2491.z;
    float _S2500 = fTmp1B_13 * _S2482 + x_39 * s_diff_fS2_T_3 + y_16 * _S2487 + 0.54627424478530884f * (_S2497.x + _S2497.y + _S2497.z);
    float _S2501 = fTmp1B_13 * _S2483 + y_16 * s_diff_fS2_T_3 + x_39 * s_diff_fC2_T_3 + 0.54627424478530884f * (_S2489.x + _S2489.y + _S2489.z);
    float _S2502 = y_16 * - _S2501;
    float _S2503 = x_39 * _S2501;
    float _S2504 = z_13 * (1.86588168144226074f * (z_13 * _S2486) + -2.28522896766662598f * (y_16 * _S2484 + x_39 * _S2485) + 0.94617468118667603f * (_S2493.x + _S2493.y + _S2493.z));
    float3  _S2505 = make_float3 (0.48860251903533936f) * _S2466.differential_0;
    float3  _S2506 = - _S2505;
    float3  _S2507 = _S2446 * _S2506;
    float3  _S2508 = (*sh_coeffs_13)[int(3)] * _S2506;
    float3  _S2509 = _S2448 * _S2505;
    float3  _S2510 = (*sh_coeffs_13)[int(2)] * _S2505;
    float3  _S2511 = _S2450 * _S2505;
    float3  _S2512 = (*sh_coeffs_13)[int(1)] * _S2505;
    float _S2513 = (_S2457 * _S2486 + 1.44530570507049561f * (fS1_13 * _S2482 + fC1_13 * _S2483) + -1.09254848957061768f * (y_16 * _S2498 + x_39 * _S2499) + _S2504 + _S2504 + _S2510.x + _S2510.y + _S2510.z) / _S2447;
    float _S2514 = _S2445 * _S2513;
    float _S2515 = (fTmp0C_13 * _S2484 + fC1_13 * s_diff_fS2_T_3 + fS1_13 * _S2487 + fTmp0B_13 * _S2498 + _S2451 * _S2500 + _S2502 + _S2502 + - (_S2512.x + _S2512.y + _S2512.z)) / _S2447;
    float _S2516 = _S2445 * _S2515;
    float _S2517 = (fTmp0C_13 * _S2485 + fS1_13 * s_diff_fS2_T_3 + fC1_13 * s_diff_fC2_T_3 + fTmp0B_13 * _S2499 + 2.0f * (y_16 * _S2500) + _S2503 + _S2503 + _S2508.x + _S2508.y + _S2508.z) / _S2447;
    float _S2518 = _S2445 * _S2517;
    float _S2519 = _S2443 * - _S2513 + _S2442 * - _S2515 + _S2441 * - _S2517;
    DiffPair_float_0 _S2520;
    (&_S2520)->primal_0 = _S2444;
    (&_S2520)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2520, _S2519);
    float _S2521 = _S2443 * _S2520.differential_0;
    float _S2522 = _S2442 * _S2520.differential_0;
    float _S2523 = _S2441 * _S2520.differential_0;
    float3  _S2524 = make_float3 (0.282094806432724f) * _S2466.differential_0;
    float3  _S2525 = make_float3 (_S2518 + _S2523 + _S2523, _S2516 + _S2522 + _S2522, _S2514 + _S2521 + _S2521);
    float3  _S2526 = - - _S2525;
    Matrix<float, 3, 3>  _S2527 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2528;
    (&_S2528)->primal_0 = _S2439;
    (&_S2528)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2529;
    (&_S2529)->primal_0 = t_16;
    (&_S2529)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2528, &_S2529, _S2526);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2530 = _S2529;
    Matrix<float, 3, 3>  _S2531 = transpose_0(_S2528.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2532;
    (&_S2532)->primal_0 = _S2438;
    (&_S2532)->differential_0 = _S2390;
    s_bwd_prop_exp_1(&_S2532, v_conic_3);
    float3  _S2533 = - _S2532.differential_0;
    float _S2534 = 0.5f * v_depth_3;
    DiffPair_float_0 _S2535;
    (&_S2535)->primal_0 = _S2437;
    (&_S2535)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2535, _S2534);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2536;
    (&_S2536)->primal_0 = mean_c_13;
    (&_S2536)->differential_0 = _S2390;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2537;
    (&_S2537)->primal_0 = mean_c_13;
    (&_S2537)->differential_0 = _S2390;
    s_bwd_prop_dot_0(&_S2536, &_S2537, _S2535.differential_0);
    DiffPair_float_0 _S2538;
    (&_S2538)->primal_0 = _S2436;
    (&_S2538)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2538, 0.0f);
    DiffPair_float_0 _S2539;
    (&_S2539)->primal_0 = _S2435;
    (&_S2539)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2539, 0.0f);
    DiffPair_float_0 _S2540;
    (&_S2540)->primal_0 = 3.32999992370605469f;
    (&_S2540)->differential_0 = 0.0f;
    DiffPair_float_0 _S2541;
    (&_S2541)->primal_0 = _S2434;
    (&_S2541)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2540, &_S2541, 0.0f);
    DiffPair_float_0 _S2542;
    (&_S2542)->primal_0 = _S2433;
    (&_S2542)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2542, _S2541.differential_0);
    float _S2543 = 2.0f * _S2542.differential_0;
    DiffPair_float_0 _S2544;
    (&_S2544)->primal_0 = _S2432;
    (&_S2544)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2544, _S2543);
    float _S2545 = v_opacity_3 + 254.9999847412109375f * _S2544.differential_0;
    FixedArray<float3 , 16>  _S2546;
    _S2546[int(0)] = _S2390;
    _S2546[int(1)] = _S2390;
    _S2546[int(2)] = _S2390;
    _S2546[int(3)] = _S2390;
    _S2546[int(4)] = _S2390;
    _S2546[int(5)] = _S2390;
    _S2546[int(6)] = _S2390;
    _S2546[int(7)] = _S2390;
    _S2546[int(8)] = _S2390;
    _S2546[int(9)] = _S2390;
    _S2546[int(10)] = _S2390;
    _S2546[int(11)] = _S2390;
    _S2546[int(12)] = _S2390;
    _S2546[int(13)] = _S2390;
    _S2546[int(14)] = _S2390;
    _S2546[int(15)] = _S2390;
    _S2546[int(7)] = _S2490;
    _S2546[int(0)] = _S2524;
    _S2546[int(1)] = _S2511;
    _S2546[int(2)] = _S2509;
    _S2546[int(3)] = _S2507;
    _S2546[int(4)] = _S2496;
    _S2546[int(5)] = _S2494;
    _S2546[int(6)] = _S2492;
    _S2546[int(15)] = _S2468;
    _S2546[int(8)] = _S2488;
    _S2546[int(9)] = _S2480;
    _S2546[int(10)] = _S2478;
    _S2546[int(11)] = _S2476;
    _S2546[int(12)] = _S2474;
    _S2546[int(13)] = _S2472;
    _S2546[int(14)] = _S2470;
    float3  _S2547 = _S2546[int(0)];
    float3  _S2548 = _S2546[int(1)];
    float3  _S2549 = _S2546[int(2)];
    float3  _S2550 = _S2546[int(3)];
    float3  _S2551 = _S2546[int(4)];
    float3  _S2552 = _S2546[int(5)];
    float3  _S2553 = _S2546[int(6)];
    float3  _S2554 = _S2546[int(7)];
    float3  _S2555 = _S2546[int(8)];
    float3  _S2556 = _S2546[int(9)];
    float3  _S2557 = _S2546[int(10)];
    float3  _S2558 = _S2546[int(11)];
    float3  _S2559 = _S2546[int(12)];
    float3  _S2560 = _S2546[int(13)];
    float3  _S2561 = _S2546[int(14)];
    float3  _S2562 = _S2546[int(15)];
    float3  _S2563 = _S2537.differential_0 + _S2536.differential_0;
    float2  _S2564 = make_float2 (0.0f, _S2538.differential_0);
    float2  _S2565 = make_float2 (_S2539.differential_0, 0.0f);
    float _S2566;
    if(antialiased_13)
    {
        float _S2567 = _S2429 * _S2545;
        _S2431 = _S2426 * _S2545;
        _S2566 = _S2567;
    }
    else
    {
        _S2431 = _S2545;
        _S2566 = 0.0f;
    }
    float _S2568 = - (_S2431 / _S2430);
    DiffPair_float_0 _S2569;
    (&_S2569)->primal_0 = _S2427;
    (&_S2569)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2569, _S2568);
    float _S2570 = - _S2569.differential_0;
    DiffPair_float_0 _S2571;
    (&_S2571)->primal_0 = _S2425;
    (&_S2571)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2571, _S2566);
    DiffPair_float_0 _S2572;
    (&_S2572)->primal_0 = 0.0f;
    (&_S2572)->differential_0 = 0.0f;
    DiffPair_float_0 _S2573;
    (&_S2573)->primal_0 = _S2423;
    (&_S2573)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2572, &_S2573, _S2571.differential_0);
    float _S2574 = _S2573.differential_0 / _S2424;
    float s_diff_det_blur_T_0 = det_orig_14 * - _S2574;
    float s_diff_det_orig_T_3 = det_blur_9 * _S2574;
    float _S2575 = - s_diff_det_blur_T_0;
    float _S2576 = _S2418 * s_diff_det_blur_T_0;
    float _S2577 = _S2420 * s_diff_det_blur_T_0;
    Matrix<float, 2, 2>  _S2578 = _S2387;
    _S2578[int(1)] = _S2564;
    _S2578[int(0)] = _S2565;
    float _S2579 = _S2577 + _S2578.rows[int(0)].x;
    float _S2580 = _S2575 + - s_diff_det_orig_T_3;
    float _S2581 = _S2413._S2386.rows[int(0)].y * _S2580;
    float _S2582 = _S2413._S2386.rows[int(1)].x * _S2580;
    float _S2583 = _S2413._S2386.rows[int(1)].y * s_diff_det_orig_T_3;
    float _S2584 = _S2576 + _S2578.rows[int(1)].y + _S2413._S2386.rows[int(0)].x * s_diff_det_orig_T_3;
    float2  _S2585 = _S2389;
    *&((&_S2585)->x) = _S2581;
    *&((&_S2585)->y) = _S2584;
    float _S2586 = _S2579 + _S2583;
    float2  _S2587 = _S2389;
    *&((&_S2587)->y) = _S2582;
    *&((&_S2587)->x) = _S2586;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2588;
    (&_S2588)->primal_0 = R_19;
    (&_S2588)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2589;
    (&_S2589)->primal_0 = _S2407;
    (&_S2589)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2588, &_S2589, _S2390);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2590;
    (&_S2590)->primal_0 = R_19;
    (&_S2590)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2591;
    (&_S2591)->primal_0 = _S2404;
    (&_S2591)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2590, &_S2591, _S2390);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2592;
    (&_S2592)->primal_0 = R_19;
    (&_S2592)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2593;
    (&_S2593)->primal_0 = _S2401;
    (&_S2593)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2592, &_S2593, _S2390);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2594;
    (&_S2594)->primal_0 = R_19;
    (&_S2594)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2595;
    (&_S2595)->primal_0 = _S2406;
    (&_S2595)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2594, &_S2595, _S2390);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2596;
    (&_S2596)->primal_0 = R_19;
    (&_S2596)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2597;
    (&_S2597)->primal_0 = _S2403;
    (&_S2597)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2596, &_S2597, _S2390);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2598;
    (&_S2598)->primal_0 = R_19;
    (&_S2598)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2599;
    (&_S2599)->primal_0 = _S2400;
    (&_S2599)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2598, &_S2599, _S2390);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2600;
    (&_S2600)->primal_0 = R_19;
    (&_S2600)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2601;
    (&_S2601)->primal_0 = _S2408.p_0[0U];
    (&_S2601)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2600, &_S2601, _S2390);
    float3  _S2602 = - _S2589.differential_0 + _S2595.differential_0;
    float3  _S2603 = _S2416 * _S2602;
    float3  _S2604 = _S2393.rows[2U] * _S2602;
    float _S2605 = _S2398 * (_S2604.x + _S2604.y + _S2604.z);
    float3  _S2606 = - _S2591.differential_0 + _S2597.differential_0;
    float3  _S2607 = _S2415 * _S2606;
    float3  _S2608 = _S2393.rows[1U] * _S2606;
    float _S2609 = _S2398 * (_S2608.x + _S2608.y + _S2608.z);
    float3  _S2610 = - _S2593.differential_0 + _S2599.differential_0;
    float3  _S2611 = _S2414 * _S2610;
    float3  _S2612 = _S2393.rows[0U] * _S2610;
    float _S2613 = _S2398 * (_S2612.x + _S2612.y + _S2612.z);
    Matrix<float, 3, 3>  _S2614 = _S2527;
    _S2614[2U] = _S2603;
    _S2614[1U] = _S2607;
    _S2614[0U] = _S2611;
    Matrix<float, 3, 3>  _S2615 = transpose_0(transpose_0(_S2614));
    float _S2616 = 2.0f * - _S2615.rows[int(2)].z;
    float _S2617 = 2.0f * _S2615.rows[int(2)].y;
    float _S2618 = 2.0f * _S2615.rows[int(2)].x;
    float _S2619 = 2.0f * _S2615.rows[int(1)].z;
    float _S2620 = 2.0f * - _S2615.rows[int(1)].y;
    float _S2621 = 2.0f * _S2615.rows[int(1)].x;
    float _S2622 = 2.0f * _S2615.rows[int(0)].z;
    float _S2623 = 2.0f * _S2615.rows[int(0)].y;
    float _S2624 = 2.0f * - _S2615.rows[int(0)].x;
    float _S2625 = - _S2621 + _S2623;
    float _S2626 = _S2618 + - _S2622;
    float _S2627 = - _S2617 + _S2619;
    float _S2628 = _S2617 + _S2619;
    float _S2629 = _S2618 + _S2622;
    float _S2630 = _S2621 + _S2623;
    float _S2631 = quat_17.w * (_S2620 + _S2624);
    float _S2632 = quat_17.z * (_S2616 + _S2624);
    float _S2633 = quat_17.y * (_S2616 + _S2620);
    float _S2634 = quat_17.x * _S2625 + quat_17.z * _S2628 + quat_17.y * _S2629 + _S2631 + _S2631;
    float _S2635 = quat_17.x * _S2626 + quat_17.w * _S2628 + quat_17.y * _S2630 + _S2632 + _S2632;
    float _S2636 = quat_17.x * _S2627 + quat_17.w * _S2629 + quat_17.z * _S2630 + _S2633 + _S2633;
    float _S2637 = quat_17.w * _S2625 + quat_17.z * _S2626 + quat_17.y * _S2627;
    float3  _S2638 = _S2390;
    *&((&_S2638)->z) = _S2605;
    *&((&_S2638)->y) = _S2609;
    *&((&_S2638)->x) = _S2613;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2639;
    (&_S2639)->primal_0 = scale_16;
    (&_S2639)->differential_0 = _S2390;
    s_bwd_prop_exp_1(&_S2639, _S2638);
    Matrix<float, 2, 2>  _S2640 = _S2387;
    _S2640[int(1)] = _S2585;
    _S2640[int(0)] = _S2587;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2641;
    (&_S2641)->primal_0 = R_19;
    (&_S2641)->differential_0 = _S2527;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2642;
    (&_S2642)->primal_0 = mean_14;
    (&_S2642)->differential_0 = _S2390;
    s_bwd_prop_mul_1(&_S2641, &_S2642, _S2563);
    float3  _S2643 = _S2563 + _S2530.differential_0;
    Matrix<float, 3, 3>  _S2644 = _S2588.differential_0 + _S2590.differential_0 + _S2592.differential_0 + _S2594.differential_0 + _S2596.differential_0 + _S2598.differential_0 + _S2600.differential_0 + _S2641.differential_0 + _S2531;
    float3  _S2645 = _S2639.differential_0 + _S2533;
    float4  _S2646 = make_float4 (0.0f);
    *&((&_S2646)->w) = _S2634;
    *&((&_S2646)->z) = _S2635;
    *&((&_S2646)->y) = _S2636;
    *&((&_S2646)->x) = _S2637;
    float4  _S2647 = _S2646;
    float3  _S2648 = _S2589.differential_0 + _S2595.differential_0 + _S2591.differential_0 + _S2597.differential_0 + _S2593.differential_0 + _S2599.differential_0 + _S2642.differential_0 + _S2525;
    *v_mean_3 = _S2648;
    *v_quat_3 = _S2647;
    *v_scale_3 = _S2645;
    *v_in_opacity_3 = _S2570;
    (*v_sh_coeffs_3)[int(0)] = _S2547;
    (*v_sh_coeffs_3)[int(1)] = _S2548;
    (*v_sh_coeffs_3)[int(2)] = _S2549;
    (*v_sh_coeffs_3)[int(3)] = _S2550;
    (*v_sh_coeffs_3)[int(4)] = _S2551;
    (*v_sh_coeffs_3)[int(5)] = _S2552;
    (*v_sh_coeffs_3)[int(6)] = _S2553;
    (*v_sh_coeffs_3)[int(7)] = _S2554;
    (*v_sh_coeffs_3)[int(8)] = _S2555;
    (*v_sh_coeffs_3)[int(9)] = _S2556;
    (*v_sh_coeffs_3)[int(10)] = _S2557;
    (*v_sh_coeffs_3)[int(11)] = _S2558;
    (*v_sh_coeffs_3)[int(12)] = _S2559;
    (*v_sh_coeffs_3)[int(13)] = _S2560;
    (*v_sh_coeffs_3)[int(14)] = _S2561;
    (*v_sh_coeffs_3)[int(15)] = _S2562;
    *v_R_5 = _S2644;
    *v_t_4 = _S2643;
    return;
}

struct s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0
{
    Matrix<float, 2, 2>  _S2649;
};

inline __device__ void projection_3dgs_eval3d_fisheye_vjp(bool antialiased_14, float3  mean_15, float4  quat_18, float3  scale_17, float in_opacity_14, FixedArray<float3 , 16>  * sh_coeffs_14, Matrix<float, 3, 3>  R_20, float3  t_17, float fx_23, float fy_23, float cx_18, float cy_18, FixedArray<float, 10>  * dist_coeffs_27, uint image_width_14, uint image_height_14, float2  v_mean2d_4, float v_depth_4, float3  v_conic_4, float v_opacity_4, float3  v_rgb_4, float3  * v_mean_4, float4  * v_quat_4, float3  * v_scale_4, float * v_in_opacity_4, FixedArray<float3 , 16>  * v_sh_coeffs_4, Matrix<float, 3, 3>  * v_R_6, float3  * v_t_5)
{
    Matrix<float, 2, 2>  _S2650 = makeMatrix<float, 2, 2> (0.0f);
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2651;
    (&_S2651)->_S2649 = _S2650;
    float2  _S2652 = make_float2 (0.0f);
    float3  _S2653 = make_float3 (0.0f);
    float4  intrins_13 = make_float4 (fx_23, fy_23, cx_18, cy_18);
    float3  _S2654 = s_primal_ctx_exp_0(scale_17);
    float _S2655 = quat_18.y;
    float x2_18 = _S2655 * _S2655;
    float y2_18 = quat_18.z * quat_18.z;
    float z2_32 = quat_18.w * quat_18.w;
    float xy_18 = quat_18.y * quat_18.z;
    float xz_18 = quat_18.y * quat_18.w;
    float yz_18 = quat_18.z * quat_18.w;
    float wx_18 = quat_18.x * quat_18.y;
    float wy_18 = quat_18.x * quat_18.z;
    float wz_18 = quat_18.x * quat_18.w;
    Matrix<float, 3, 3>  _S2656 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_18 + z2_32), 2.0f * (xy_18 + wz_18), 2.0f * (xz_18 - wy_18), 2.0f * (xy_18 - wz_18), 1.0f - 2.0f * (x2_18 + z2_32), 2.0f * (yz_18 + wx_18), 2.0f * (xz_18 + wy_18), 2.0f * (yz_18 - wx_18), 1.0f - 2.0f * (x2_18 + y2_18))));
    FixedArray<float3 , 7>  _S2657 = {
        _S2653, _S2653, _S2653, _S2653, _S2653, _S2653, _S2653
    };
    FixedArray<float, 7>  _S2658 = {
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f
    };
    SigmaPoints_0 _S2659;
    (&_S2659)->p_0 = _S2657;
    (&_S2659)->w_mean_0 = _S2658;
    (&_S2659)->w_cov_0 = _S2658;
    (&_S2659)->p_0[int(0)] = mean_15;
    SigmaPoints_0 _S2660 = _S2659;
    (&_S2660)->w_mean_0[int(0)] = 0.0f;
    (&_S2660)->w_cov_0[int(0)] = 2.0f;
    float _S2661 = s_primal_ctx_sqrt_0(3.0f);
    float _S2662 = _S2661 * _S2654.x;
    float3  delta_15 = make_float3 (_S2662) * _S2656.rows[0U];
    float3  _S2663 = mean_15 + delta_15;
    float3  _S2664 = mean_15 - delta_15;
    float _S2665 = _S2661 * _S2654.y;
    float3  delta_16 = make_float3 (_S2665) * _S2656.rows[1U];
    float3  _S2666 = mean_15 + delta_16;
    float3  _S2667 = mean_15 - delta_16;
    float _S2668 = _S2661 * _S2654.z;
    float3  delta_17 = make_float3 (_S2668) * _S2656.rows[2U];
    float3  _S2669 = mean_15 + delta_17;
    float3  _S2670 = mean_15 - delta_17;
    (&_S2660)->w_mean_0[1U] = 0.1666666716337204f;
    (&_S2660)->w_cov_0[1U] = 0.1666666716337204f;
    (&_S2660)->w_mean_0[2U] = 0.1666666716337204f;
    (&_S2660)->w_cov_0[2U] = 0.1666666716337204f;
    (&_S2660)->w_mean_0[3U] = 0.1666666716337204f;
    (&_S2660)->w_cov_0[3U] = 0.1666666716337204f;
    (&_S2660)->w_mean_0[4U] = 0.1666666716337204f;
    (&_S2660)->w_cov_0[4U] = 0.1666666716337204f;
    (&_S2660)->w_mean_0[5U] = 0.1666666716337204f;
    (&_S2660)->w_cov_0[5U] = 0.1666666716337204f;
    (&_S2660)->w_mean_0[6U] = 0.1666666716337204f;
    (&_S2660)->w_cov_0[6U] = 0.1666666716337204f;
    SigmaPoints_0 _S2671 = _S2659;
    (&_S2660)->p_0[0U] = s_primal_ctx_mul_0(R_20, _S2659.p_0[0U]) + t_17;
    (&_S2660)->p_0[1U] = s_primal_ctx_mul_0(R_20, _S2663) + t_17;
    (&_S2660)->p_0[2U] = s_primal_ctx_mul_0(R_20, _S2666) + t_17;
    (&_S2660)->p_0[3U] = s_primal_ctx_mul_0(R_20, _S2669) + t_17;
    (&_S2660)->p_0[4U] = s_primal_ctx_mul_0(R_20, _S2664) + t_17;
    (&_S2660)->p_0[5U] = s_primal_ctx_mul_0(R_20, _S2667) + t_17;
    (&_S2660)->p_0[6U] = s_primal_ctx_mul_0(R_20, _S2670) + t_17;
    float2  _S2672 = _S2652;
    Matrix<float, 2, 2>  _S2673 = _S2650;
    SigmaPoints_0 _S2674 = _S2660;
    bool _S2675 = fisheye_proj_3dgs_ut_1(&_S2674, intrins_13, dist_coeffs_27, &_S2673, &_S2672);
    (&_S2651)->_S2649 = _S2673;
    s_bwd_prop_projection_3dgs_eval3d_fisheye_differentiable_Intermediates_0 _S2676 = _S2651;
    float3  mean_c_14 = s_primal_ctx_mul_0(R_20, mean_15) + t_17;
    float3  _S2677 = make_float3 (_S2662);
    float3  _S2678 = make_float3 (_S2665);
    float3  _S2679 = make_float3 (_S2668);
    float _S2680 = _S2651._S2649.rows[int(0)].y * _S2651._S2649.rows[int(1)].x;
    float det_orig_15 = _S2651._S2649.rows[int(0)].x * _S2651._S2649.rows[int(1)].y - _S2680;
    float _S2681 = _S2651._S2649.rows[int(0)].x + 0.30000001192092896f;
    Matrix<float, 2, 2>  _S2682 = _S2651._S2649;
    *&(((&_S2682)->rows + (int(0)))->x) = _S2681;
    float _S2683 = _S2651._S2649.rows[int(1)].y + 0.30000001192092896f;
    *&(((&_S2682)->rows + (int(1)))->y) = _S2683;
    Matrix<float, 2, 2>  _S2684 = _S2682;
    Matrix<float, 2, 2>  _S2685 = _S2682;
    float det_blur_10 = _S2681 * _S2683 - _S2680;
    float _S2686 = det_orig_15 / det_blur_10;
    float _S2687 = det_blur_10 * det_blur_10;
    float _S2688 = s_primal_ctx_max_0(0.0f, _S2686);
    float _S2689 = s_primal_ctx_sqrt_0(_S2688);
    float _S2690 = - in_opacity_14;
    float _S2691 = 1.0f + s_primal_ctx_exp_1(_S2690);
    float _S2692 = 1.0f / _S2691;
    float _S2693 = _S2691 * _S2691;
    float _S2694;
    if(antialiased_14)
    {
        _S2694 = _S2692 * _S2689;
    }
    else
    {
        _S2694 = _S2692;
    }
    float _S2695 = _S2694 / 0.00392156885936856f;
    float _S2696 = 2.0f * s_primal_ctx_log_0(_S2695);
    float _S2697 = s_primal_ctx_sqrt_0(_S2696);
    float _S2698 = _S2684.rows[int(0)].x;
    float _S2699 = _S2685.rows[int(1)].y;
    float _S2700 = s_primal_ctx_dot_0(mean_c_14, mean_c_14) + 9.99999997475242708e-07f;
    float3  _S2701 = - scale_17;
    Matrix<float, 3, 3>  _S2702 = transpose_0(R_20);
    float3  _S2703 = mean_15 - - s_primal_ctx_mul_0(_S2702, t_17);
    float _S2704 = _S2703.x;
    float _S2705 = _S2703.y;
    float _S2706 = _S2703.z;
    float _S2707 = _S2704 * _S2704 + _S2705 * _S2705 + _S2706 * _S2706;
    float _S2708 = s_primal_ctx_sqrt_0(_S2707);
    float x_40 = _S2704 / _S2708;
    float3  _S2709 = make_float3 (x_40);
    float _S2710 = _S2708 * _S2708;
    float y_17 = _S2705 / _S2708;
    float z_14 = _S2706 / _S2708;
    float3  _S2711 = make_float3 (z_14);
    float _S2712 = - y_17;
    float3  _S2713 = make_float3 (_S2712);
    float z2_33 = z_14 * z_14;
    float fTmp0B_14 = -1.09254848957061768f * z_14;
    float fC1_14 = x_40 * x_40 - y_17 * y_17;
    float _S2714 = 2.0f * x_40;
    float fS1_14 = _S2714 * y_17;
    float pSH6_4 = 0.94617468118667603f * z2_33 - 0.31539157032966614f;
    float3  _S2715 = make_float3 (pSH6_4);
    float pSH7_4 = fTmp0B_14 * x_40;
    float3  _S2716 = make_float3 (pSH7_4);
    float pSH5_4 = fTmp0B_14 * y_17;
    float3  _S2717 = make_float3 (pSH5_4);
    float pSH8_4 = 0.54627424478530884f * fC1_14;
    float3  _S2718 = make_float3 (pSH8_4);
    float pSH4_4 = 0.54627424478530884f * fS1_14;
    float3  _S2719 = make_float3 (pSH4_4);
    float fTmp0C_14 = -2.28522896766662598f * z2_33 + 0.4570457935333252f;
    float fTmp1B_14 = 1.44530570507049561f * z_14;
    float _S2720 = 1.86588168144226074f * z2_33 - 1.11952900886535645f;
    float pSH12_4 = z_14 * _S2720;
    float3  _S2721 = make_float3 (pSH12_4);
    float pSH13_4 = fTmp0C_14 * x_40;
    float3  _S2722 = make_float3 (pSH13_4);
    float pSH11_4 = fTmp0C_14 * y_17;
    float3  _S2723 = make_float3 (pSH11_4);
    float pSH14_4 = fTmp1B_14 * fC1_14;
    float3  _S2724 = make_float3 (pSH14_4);
    float pSH10_4 = fTmp1B_14 * fS1_14;
    float3  _S2725 = make_float3 (pSH10_4);
    float pSH15_4 = -0.59004360437393188f * (x_40 * fC1_14 - y_17 * fS1_14);
    float3  _S2726 = make_float3 (pSH15_4);
    float pSH9_4 = -0.59004360437393188f * (x_40 * fS1_14 + y_17 * fC1_14);
    float3  _S2727 = make_float3 (pSH9_4);
    float3  _S2728 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2729;
    (&_S2729)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_14)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S2712) * (*sh_coeffs_14)[int(1)] + make_float3 (z_14) * (*sh_coeffs_14)[int(2)] - make_float3 (x_40) * (*sh_coeffs_14)[int(3)]) + (make_float3 (pSH4_4) * (*sh_coeffs_14)[int(4)] + make_float3 (pSH5_4) * (*sh_coeffs_14)[int(5)] + make_float3 (pSH6_4) * (*sh_coeffs_14)[int(6)] + make_float3 (pSH7_4) * (*sh_coeffs_14)[int(7)] + make_float3 (pSH8_4) * (*sh_coeffs_14)[int(8)]) + (make_float3 (pSH9_4) * (*sh_coeffs_14)[int(9)] + make_float3 (pSH10_4) * (*sh_coeffs_14)[int(10)] + make_float3 (pSH11_4) * (*sh_coeffs_14)[int(11)] + make_float3 (pSH12_4) * (*sh_coeffs_14)[int(12)] + make_float3 (pSH13_4) * (*sh_coeffs_14)[int(13)] + make_float3 (pSH14_4) * (*sh_coeffs_14)[int(14)] + make_float3 (pSH15_4) * (*sh_coeffs_14)[int(15)]) + make_float3 (0.5f);
    (&_S2729)->differential_0 = _S2653;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2730;
    (&_S2730)->primal_0 = _S2728;
    (&_S2730)->differential_0 = _S2653;
    s_bwd_prop_max_0(&_S2729, &_S2730, v_rgb_4);
    float3  _S2731 = _S2726 * _S2729.differential_0;
    float3  _S2732 = (*sh_coeffs_14)[int(15)] * _S2729.differential_0;
    float3  _S2733 = _S2724 * _S2729.differential_0;
    float3  _S2734 = (*sh_coeffs_14)[int(14)] * _S2729.differential_0;
    float3  _S2735 = _S2722 * _S2729.differential_0;
    float3  _S2736 = (*sh_coeffs_14)[int(13)] * _S2729.differential_0;
    float3  _S2737 = _S2721 * _S2729.differential_0;
    float3  _S2738 = (*sh_coeffs_14)[int(12)] * _S2729.differential_0;
    float3  _S2739 = _S2723 * _S2729.differential_0;
    float3  _S2740 = (*sh_coeffs_14)[int(11)] * _S2729.differential_0;
    float3  _S2741 = _S2725 * _S2729.differential_0;
    float3  _S2742 = (*sh_coeffs_14)[int(10)] * _S2729.differential_0;
    float3  _S2743 = _S2727 * _S2729.differential_0;
    float3  _S2744 = (*sh_coeffs_14)[int(9)] * _S2729.differential_0;
    float s_diff_fS2_T_4 = -0.59004360437393188f * (_S2744.x + _S2744.y + _S2744.z);
    float s_diff_fC2_T_4 = -0.59004360437393188f * (_S2732.x + _S2732.y + _S2732.z);
    float _S2745 = _S2742.x + _S2742.y + _S2742.z;
    float _S2746 = _S2734.x + _S2734.y + _S2734.z;
    float _S2747 = _S2740.x + _S2740.y + _S2740.z;
    float _S2748 = _S2736.x + _S2736.y + _S2736.z;
    float _S2749 = _S2738.x + _S2738.y + _S2738.z;
    float _S2750 = - s_diff_fC2_T_4;
    float3  _S2751 = _S2718 * _S2729.differential_0;
    float3  _S2752 = (*sh_coeffs_14)[int(8)] * _S2729.differential_0;
    float3  _S2753 = _S2716 * _S2729.differential_0;
    float3  _S2754 = (*sh_coeffs_14)[int(7)] * _S2729.differential_0;
    float3  _S2755 = _S2715 * _S2729.differential_0;
    float3  _S2756 = (*sh_coeffs_14)[int(6)] * _S2729.differential_0;
    float3  _S2757 = _S2717 * _S2729.differential_0;
    float3  _S2758 = (*sh_coeffs_14)[int(5)] * _S2729.differential_0;
    float3  _S2759 = _S2719 * _S2729.differential_0;
    float3  _S2760 = (*sh_coeffs_14)[int(4)] * _S2729.differential_0;
    float _S2761 = _S2758.x + _S2758.y + _S2758.z;
    float _S2762 = _S2754.x + _S2754.y + _S2754.z;
    float _S2763 = fTmp1B_14 * _S2745 + x_40 * s_diff_fS2_T_4 + y_17 * _S2750 + 0.54627424478530884f * (_S2760.x + _S2760.y + _S2760.z);
    float _S2764 = fTmp1B_14 * _S2746 + y_17 * s_diff_fS2_T_4 + x_40 * s_diff_fC2_T_4 + 0.54627424478530884f * (_S2752.x + _S2752.y + _S2752.z);
    float _S2765 = y_17 * - _S2764;
    float _S2766 = x_40 * _S2764;
    float _S2767 = z_14 * (1.86588168144226074f * (z_14 * _S2749) + -2.28522896766662598f * (y_17 * _S2747 + x_40 * _S2748) + 0.94617468118667603f * (_S2756.x + _S2756.y + _S2756.z));
    float3  _S2768 = make_float3 (0.48860251903533936f) * _S2729.differential_0;
    float3  _S2769 = - _S2768;
    float3  _S2770 = _S2709 * _S2769;
    float3  _S2771 = (*sh_coeffs_14)[int(3)] * _S2769;
    float3  _S2772 = _S2711 * _S2768;
    float3  _S2773 = (*sh_coeffs_14)[int(2)] * _S2768;
    float3  _S2774 = _S2713 * _S2768;
    float3  _S2775 = (*sh_coeffs_14)[int(1)] * _S2768;
    float _S2776 = (_S2720 * _S2749 + 1.44530570507049561f * (fS1_14 * _S2745 + fC1_14 * _S2746) + -1.09254848957061768f * (y_17 * _S2761 + x_40 * _S2762) + _S2767 + _S2767 + _S2773.x + _S2773.y + _S2773.z) / _S2710;
    float _S2777 = _S2708 * _S2776;
    float _S2778 = (fTmp0C_14 * _S2747 + fC1_14 * s_diff_fS2_T_4 + fS1_14 * _S2750 + fTmp0B_14 * _S2761 + _S2714 * _S2763 + _S2765 + _S2765 + - (_S2775.x + _S2775.y + _S2775.z)) / _S2710;
    float _S2779 = _S2708 * _S2778;
    float _S2780 = (fTmp0C_14 * _S2748 + fS1_14 * s_diff_fS2_T_4 + fC1_14 * s_diff_fC2_T_4 + fTmp0B_14 * _S2762 + 2.0f * (y_17 * _S2763) + _S2766 + _S2766 + _S2771.x + _S2771.y + _S2771.z) / _S2710;
    float _S2781 = _S2708 * _S2780;
    float _S2782 = _S2706 * - _S2776 + _S2705 * - _S2778 + _S2704 * - _S2780;
    DiffPair_float_0 _S2783;
    (&_S2783)->primal_0 = _S2707;
    (&_S2783)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2783, _S2782);
    float _S2784 = _S2706 * _S2783.differential_0;
    float _S2785 = _S2705 * _S2783.differential_0;
    float _S2786 = _S2704 * _S2783.differential_0;
    float3  _S2787 = make_float3 (0.282094806432724f) * _S2729.differential_0;
    float3  _S2788 = make_float3 (_S2781 + _S2786 + _S2786, _S2779 + _S2785 + _S2785, _S2777 + _S2784 + _S2784);
    float3  _S2789 = - - _S2788;
    Matrix<float, 3, 3>  _S2790 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2791;
    (&_S2791)->primal_0 = _S2702;
    (&_S2791)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2792;
    (&_S2792)->primal_0 = t_17;
    (&_S2792)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2791, &_S2792, _S2789);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2793 = _S2792;
    Matrix<float, 3, 3>  _S2794 = transpose_0(_S2791.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2795;
    (&_S2795)->primal_0 = _S2701;
    (&_S2795)->differential_0 = _S2653;
    s_bwd_prop_exp_1(&_S2795, v_conic_4);
    float3  _S2796 = - _S2795.differential_0;
    float _S2797 = 0.5f * v_depth_4;
    DiffPair_float_0 _S2798;
    (&_S2798)->primal_0 = _S2700;
    (&_S2798)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2798, _S2797);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2799;
    (&_S2799)->primal_0 = mean_c_14;
    (&_S2799)->differential_0 = _S2653;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2800;
    (&_S2800)->primal_0 = mean_c_14;
    (&_S2800)->differential_0 = _S2653;
    s_bwd_prop_dot_0(&_S2799, &_S2800, _S2798.differential_0);
    DiffPair_float_0 _S2801;
    (&_S2801)->primal_0 = _S2699;
    (&_S2801)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2801, 0.0f);
    DiffPair_float_0 _S2802;
    (&_S2802)->primal_0 = _S2698;
    (&_S2802)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2802, 0.0f);
    DiffPair_float_0 _S2803;
    (&_S2803)->primal_0 = 3.32999992370605469f;
    (&_S2803)->differential_0 = 0.0f;
    DiffPair_float_0 _S2804;
    (&_S2804)->primal_0 = _S2697;
    (&_S2804)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S2803, &_S2804, 0.0f);
    DiffPair_float_0 _S2805;
    (&_S2805)->primal_0 = _S2696;
    (&_S2805)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2805, _S2804.differential_0);
    float _S2806 = 2.0f * _S2805.differential_0;
    DiffPair_float_0 _S2807;
    (&_S2807)->primal_0 = _S2695;
    (&_S2807)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S2807, _S2806);
    float _S2808 = v_opacity_4 + 254.9999847412109375f * _S2807.differential_0;
    FixedArray<float3 , 16>  _S2809;
    _S2809[int(0)] = _S2653;
    _S2809[int(1)] = _S2653;
    _S2809[int(2)] = _S2653;
    _S2809[int(3)] = _S2653;
    _S2809[int(4)] = _S2653;
    _S2809[int(5)] = _S2653;
    _S2809[int(6)] = _S2653;
    _S2809[int(7)] = _S2653;
    _S2809[int(8)] = _S2653;
    _S2809[int(9)] = _S2653;
    _S2809[int(10)] = _S2653;
    _S2809[int(11)] = _S2653;
    _S2809[int(12)] = _S2653;
    _S2809[int(13)] = _S2653;
    _S2809[int(14)] = _S2653;
    _S2809[int(15)] = _S2653;
    _S2809[int(7)] = _S2753;
    _S2809[int(0)] = _S2787;
    _S2809[int(1)] = _S2774;
    _S2809[int(2)] = _S2772;
    _S2809[int(3)] = _S2770;
    _S2809[int(4)] = _S2759;
    _S2809[int(5)] = _S2757;
    _S2809[int(6)] = _S2755;
    _S2809[int(15)] = _S2731;
    _S2809[int(8)] = _S2751;
    _S2809[int(9)] = _S2743;
    _S2809[int(10)] = _S2741;
    _S2809[int(11)] = _S2739;
    _S2809[int(12)] = _S2737;
    _S2809[int(13)] = _S2735;
    _S2809[int(14)] = _S2733;
    float3  _S2810 = _S2809[int(0)];
    float3  _S2811 = _S2809[int(1)];
    float3  _S2812 = _S2809[int(2)];
    float3  _S2813 = _S2809[int(3)];
    float3  _S2814 = _S2809[int(4)];
    float3  _S2815 = _S2809[int(5)];
    float3  _S2816 = _S2809[int(6)];
    float3  _S2817 = _S2809[int(7)];
    float3  _S2818 = _S2809[int(8)];
    float3  _S2819 = _S2809[int(9)];
    float3  _S2820 = _S2809[int(10)];
    float3  _S2821 = _S2809[int(11)];
    float3  _S2822 = _S2809[int(12)];
    float3  _S2823 = _S2809[int(13)];
    float3  _S2824 = _S2809[int(14)];
    float3  _S2825 = _S2809[int(15)];
    float3  _S2826 = _S2800.differential_0 + _S2799.differential_0;
    float2  _S2827 = make_float2 (0.0f, _S2801.differential_0);
    float2  _S2828 = make_float2 (_S2802.differential_0, 0.0f);
    float _S2829;
    if(antialiased_14)
    {
        float _S2830 = _S2692 * _S2808;
        _S2694 = _S2689 * _S2808;
        _S2829 = _S2830;
    }
    else
    {
        _S2694 = _S2808;
        _S2829 = 0.0f;
    }
    float _S2831 = - (_S2694 / _S2693);
    DiffPair_float_0 _S2832;
    (&_S2832)->primal_0 = _S2690;
    (&_S2832)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2832, _S2831);
    float _S2833 = - _S2832.differential_0;
    DiffPair_float_0 _S2834;
    (&_S2834)->primal_0 = _S2688;
    (&_S2834)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S2834, _S2829);
    DiffPair_float_0 _S2835;
    (&_S2835)->primal_0 = 0.0f;
    (&_S2835)->differential_0 = 0.0f;
    DiffPair_float_0 _S2836;
    (&_S2836)->primal_0 = _S2686;
    (&_S2836)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S2835, &_S2836, _S2834.differential_0);
    float _S2837 = _S2836.differential_0 / _S2687;
    float s_diff_det_blur_T_1 = det_orig_15 * - _S2837;
    float s_diff_det_orig_T_4 = det_blur_10 * _S2837;
    float _S2838 = - s_diff_det_blur_T_1;
    float _S2839 = _S2681 * s_diff_det_blur_T_1;
    float _S2840 = _S2683 * s_diff_det_blur_T_1;
    Matrix<float, 2, 2>  _S2841 = _S2650;
    _S2841[int(1)] = _S2827;
    _S2841[int(0)] = _S2828;
    float _S2842 = _S2840 + _S2841.rows[int(0)].x;
    float _S2843 = _S2838 + - s_diff_det_orig_T_4;
    float _S2844 = _S2676._S2649.rows[int(0)].y * _S2843;
    float _S2845 = _S2676._S2649.rows[int(1)].x * _S2843;
    float _S2846 = _S2676._S2649.rows[int(1)].y * s_diff_det_orig_T_4;
    float _S2847 = _S2839 + _S2841.rows[int(1)].y + _S2676._S2649.rows[int(0)].x * s_diff_det_orig_T_4;
    float2  _S2848 = _S2652;
    *&((&_S2848)->x) = _S2844;
    *&((&_S2848)->y) = _S2847;
    float _S2849 = _S2842 + _S2846;
    float2  _S2850 = _S2652;
    *&((&_S2850)->y) = _S2845;
    *&((&_S2850)->x) = _S2849;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2851;
    (&_S2851)->primal_0 = R_20;
    (&_S2851)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2852;
    (&_S2852)->primal_0 = _S2670;
    (&_S2852)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2851, &_S2852, _S2653);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2853;
    (&_S2853)->primal_0 = R_20;
    (&_S2853)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2854;
    (&_S2854)->primal_0 = _S2667;
    (&_S2854)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2853, &_S2854, _S2653);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2855;
    (&_S2855)->primal_0 = R_20;
    (&_S2855)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2856;
    (&_S2856)->primal_0 = _S2664;
    (&_S2856)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2855, &_S2856, _S2653);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2857;
    (&_S2857)->primal_0 = R_20;
    (&_S2857)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2858;
    (&_S2858)->primal_0 = _S2669;
    (&_S2858)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2857, &_S2858, _S2653);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2859;
    (&_S2859)->primal_0 = R_20;
    (&_S2859)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2860;
    (&_S2860)->primal_0 = _S2666;
    (&_S2860)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2859, &_S2860, _S2653);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2861;
    (&_S2861)->primal_0 = R_20;
    (&_S2861)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2862;
    (&_S2862)->primal_0 = _S2663;
    (&_S2862)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2861, &_S2862, _S2653);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2863;
    (&_S2863)->primal_0 = R_20;
    (&_S2863)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2864;
    (&_S2864)->primal_0 = _S2671.p_0[0U];
    (&_S2864)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2863, &_S2864, _S2653);
    float3  _S2865 = - _S2852.differential_0 + _S2858.differential_0;
    float3  _S2866 = _S2679 * _S2865;
    float3  _S2867 = _S2656.rows[2U] * _S2865;
    float _S2868 = _S2661 * (_S2867.x + _S2867.y + _S2867.z);
    float3  _S2869 = - _S2854.differential_0 + _S2860.differential_0;
    float3  _S2870 = _S2678 * _S2869;
    float3  _S2871 = _S2656.rows[1U] * _S2869;
    float _S2872 = _S2661 * (_S2871.x + _S2871.y + _S2871.z);
    float3  _S2873 = - _S2856.differential_0 + _S2862.differential_0;
    float3  _S2874 = _S2677 * _S2873;
    float3  _S2875 = _S2656.rows[0U] * _S2873;
    float _S2876 = _S2661 * (_S2875.x + _S2875.y + _S2875.z);
    Matrix<float, 3, 3>  _S2877 = _S2790;
    _S2877[2U] = _S2866;
    _S2877[1U] = _S2870;
    _S2877[0U] = _S2874;
    Matrix<float, 3, 3>  _S2878 = transpose_0(transpose_0(_S2877));
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
    float _S2894 = quat_18.w * (_S2883 + _S2887);
    float _S2895 = quat_18.z * (_S2879 + _S2887);
    float _S2896 = quat_18.y * (_S2879 + _S2883);
    float _S2897 = quat_18.x * _S2888 + quat_18.z * _S2891 + quat_18.y * _S2892 + _S2894 + _S2894;
    float _S2898 = quat_18.x * _S2889 + quat_18.w * _S2891 + quat_18.y * _S2893 + _S2895 + _S2895;
    float _S2899 = quat_18.x * _S2890 + quat_18.w * _S2892 + quat_18.z * _S2893 + _S2896 + _S2896;
    float _S2900 = quat_18.w * _S2888 + quat_18.z * _S2889 + quat_18.y * _S2890;
    float3  _S2901 = _S2653;
    *&((&_S2901)->z) = _S2868;
    *&((&_S2901)->y) = _S2872;
    *&((&_S2901)->x) = _S2876;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2902;
    (&_S2902)->primal_0 = scale_17;
    (&_S2902)->differential_0 = _S2653;
    s_bwd_prop_exp_1(&_S2902, _S2901);
    Matrix<float, 2, 2>  _S2903 = _S2650;
    _S2903[int(1)] = _S2848;
    _S2903[int(0)] = _S2850;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2904;
    (&_S2904)->primal_0 = mean_c_14;
    (&_S2904)->differential_0 = _S2653;
    s_bwd_length_impl_1(&_S2904, 0.0f);
    float3  _S2905 = _S2904.differential_0 + _S2826;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2906;
    (&_S2906)->primal_0 = R_20;
    (&_S2906)->differential_0 = _S2790;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2907;
    (&_S2907)->primal_0 = mean_15;
    (&_S2907)->differential_0 = _S2653;
    s_bwd_prop_mul_1(&_S2906, &_S2907, _S2905);
    float3  _S2908 = _S2905 + _S2793.differential_0;
    Matrix<float, 3, 3>  _S2909 = _S2851.differential_0 + _S2853.differential_0 + _S2855.differential_0 + _S2857.differential_0 + _S2859.differential_0 + _S2861.differential_0 + _S2863.differential_0 + _S2906.differential_0 + _S2794;
    float3  _S2910 = _S2902.differential_0 + _S2796;
    float4  _S2911 = make_float4 (0.0f);
    *&((&_S2911)->w) = _S2897;
    *&((&_S2911)->z) = _S2898;
    *&((&_S2911)->y) = _S2899;
    *&((&_S2911)->x) = _S2900;
    float4  _S2912 = _S2911;
    float3  _S2913 = _S2852.differential_0 + _S2858.differential_0 + _S2854.differential_0 + _S2860.differential_0 + _S2856.differential_0 + _S2862.differential_0 + _S2907.differential_0 + _S2788;
    *v_mean_4 = _S2913;
    *v_quat_4 = _S2912;
    *v_scale_4 = _S2910;
    *v_in_opacity_4 = _S2833;
    (*v_sh_coeffs_4)[int(0)] = _S2810;
    (*v_sh_coeffs_4)[int(1)] = _S2811;
    (*v_sh_coeffs_4)[int(2)] = _S2812;
    (*v_sh_coeffs_4)[int(3)] = _S2813;
    (*v_sh_coeffs_4)[int(4)] = _S2814;
    (*v_sh_coeffs_4)[int(5)] = _S2815;
    (*v_sh_coeffs_4)[int(6)] = _S2816;
    (*v_sh_coeffs_4)[int(7)] = _S2817;
    (*v_sh_coeffs_4)[int(8)] = _S2818;
    (*v_sh_coeffs_4)[int(9)] = _S2819;
    (*v_sh_coeffs_4)[int(10)] = _S2820;
    (*v_sh_coeffs_4)[int(11)] = _S2821;
    (*v_sh_coeffs_4)[int(12)] = _S2822;
    (*v_sh_coeffs_4)[int(13)] = _S2823;
    (*v_sh_coeffs_4)[int(14)] = _S2824;
    (*v_sh_coeffs_4)[int(15)] = _S2825;
    *v_R_6 = _S2909;
    *v_t_5 = _S2908;
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

inline __device__ void s_bwd_prop_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * dpquat_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpscale_0, Matrix<float, 3, 3>  _s_dOut_7)
{
    float _S2914 = (*dpquat_0).primal_0.y;
    float x2_20 = _S2914 * _S2914;
    float y2_20 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.z;
    float z2_35 = (*dpquat_0).primal_0.w * (*dpquat_0).primal_0.w;
    float xy_20 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.z;
    float xz_20 = (*dpquat_0).primal_0.y * (*dpquat_0).primal_0.w;
    float yz_20 = (*dpquat_0).primal_0.z * (*dpquat_0).primal_0.w;
    float wx_20 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.y;
    float wy_20 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.z;
    float wz_20 = (*dpquat_0).primal_0.x * (*dpquat_0).primal_0.w;
    Matrix<float, 3, 3>  _S2915 = transpose_0(transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_20 + z2_35), 2.0f * (xy_20 + wz_20), 2.0f * (xz_20 - wy_20), 2.0f * (xy_20 - wz_20), 1.0f - 2.0f * (x2_20 + z2_35), 2.0f * (yz_20 + wx_20), 2.0f * (xz_20 + wy_20), 2.0f * (yz_20 - wx_20), 1.0f - 2.0f * (x2_20 + y2_20))));
    Matrix<float, 3, 3>  _S2916 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2917;
    (&_S2917)->primal_0 = makeMatrix<float, 3, 3> ((*dpscale_0).primal_0.x, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.y, 0.0f, 0.0f, 0.0f, (*dpscale_0).primal_0.z);
    (&_S2917)->differential_0 = _S2916;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2918;
    (&_S2918)->primal_0 = _S2915;
    (&_S2918)->differential_0 = _S2916;
    s_bwd_prop_mul_4(&_S2917, &_S2918, _s_dOut_7);
    Matrix<float, 3, 3>  _S2919 = transpose_0(transpose_0(_S2918.differential_0));
    float _S2920 = 2.0f * - _S2919.rows[int(2)].z;
    float _S2921 = 2.0f * _S2919.rows[int(2)].y;
    float _S2922 = 2.0f * _S2919.rows[int(2)].x;
    float _S2923 = 2.0f * _S2919.rows[int(1)].z;
    float _S2924 = 2.0f * - _S2919.rows[int(1)].y;
    float _S2925 = 2.0f * _S2919.rows[int(1)].x;
    float _S2926 = 2.0f * _S2919.rows[int(0)].z;
    float _S2927 = 2.0f * _S2919.rows[int(0)].y;
    float _S2928 = 2.0f * - _S2919.rows[int(0)].x;
    float _S2929 = - _S2925 + _S2927;
    float _S2930 = _S2922 + - _S2926;
    float _S2931 = - _S2921 + _S2923;
    float _S2932 = _S2921 + _S2923;
    float _S2933 = _S2922 + _S2926;
    float _S2934 = _S2925 + _S2927;
    float _S2935 = (*dpquat_0).primal_0.w * (_S2924 + _S2928);
    float _S2936 = (*dpquat_0).primal_0.z * (_S2920 + _S2928);
    float _S2937 = (*dpquat_0).primal_0.y * (_S2920 + _S2924);
    float _S2938 = (*dpquat_0).primal_0.x * _S2929 + (*dpquat_0).primal_0.z * _S2932 + (*dpquat_0).primal_0.y * _S2933 + _S2935 + _S2935;
    float _S2939 = (*dpquat_0).primal_0.x * _S2930 + (*dpquat_0).primal_0.w * _S2932 + (*dpquat_0).primal_0.y * _S2934 + _S2936 + _S2936;
    float _S2940 = (*dpquat_0).primal_0.x * _S2931 + (*dpquat_0).primal_0.w * _S2933 + (*dpquat_0).primal_0.z * _S2934 + _S2937 + _S2937;
    float _S2941 = (*dpquat_0).primal_0.w * _S2929 + (*dpquat_0).primal_0.z * _S2930 + (*dpquat_0).primal_0.y * _S2931;
    float3  _S2942 = make_float3 (_S2917.differential_0.rows[int(0)].x, _S2917.differential_0.rows[int(1)].y, _S2917.differential_0.rows[int(2)].z);
    dpscale_0->primal_0 = (*dpscale_0).primal_0;
    dpscale_0->differential_0 = _S2942;
    float4  _S2943 = make_float4 (0.0f);
    *&((&_S2943)->w) = _S2938;
    *&((&_S2943)->z) = _S2939;
    *&((&_S2943)->y) = _S2940;
    *&((&_S2943)->x) = _S2941;
    dpquat_0->primal_0 = (*dpquat_0).primal_0;
    dpquat_0->differential_0 = _S2943;
    return;
}

inline __device__ void s_bwd_compute_3dgut_iscl_rot_0(DiffPair_vectorx3Cfloatx2C4x3E_0 * _S2944, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2945, Matrix<float, 3, 3>  _S2946)
{
    s_bwd_prop_compute_3dgut_iscl_rot_0(_S2944, _S2945, _S2946);
    return;
}

inline __device__ void compute_3dgut_iscl_rot_vjp(float4  quat_20, float3  scale_19, Matrix<float, 3, 3>  v_iscl_rot_0, float4  * v_quat_5, float3  * v_scale_5)
{
    float4  _S2947 = make_float4 (0.0f);
    DiffPair_vectorx3Cfloatx2C4x3E_0 dp_quat_0;
    (&dp_quat_0)->primal_0 = quat_20;
    (&dp_quat_0)->differential_0 = _S2947;
    float3  _S2948 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_scale_0;
    (&dp_scale_0)->primal_0 = scale_19;
    (&dp_scale_0)->differential_0 = _S2948;
    s_bwd_compute_3dgut_iscl_rot_0(&dp_quat_0, &dp_scale_0, v_iscl_rot_0);
    *v_quat_5 = dp_quat_0.differential_0;
    *v_scale_5 = dp_scale_0.differential_0;
    return;
}

inline __device__ void _d_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * a_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * b_0, float3  dOut_16)
{
    float _S2949 = dOut_16.y;
    float _S2950 = dOut_16.z;
    float _S2951 = dOut_16.x;
    float _S2952 = (*a_0).primal_0.z * _S2949 + - (*a_0).primal_0.y * _S2950;
    float _S2953 = - (*a_0).primal_0.z * _S2951 + (*a_0).primal_0.x * _S2950;
    float _S2954 = (*a_0).primal_0.y * _S2951 + - (*a_0).primal_0.x * _S2949;
    float3  _S2955 = make_float3 (- (*b_0).primal_0.z * _S2949 + (*b_0).primal_0.y * _S2950, (*b_0).primal_0.z * _S2951 + - (*b_0).primal_0.x * _S2950, - (*b_0).primal_0.y * _S2951 + (*b_0).primal_0.x * _S2949);
    a_0->primal_0 = (*a_0).primal_0;
    a_0->differential_0 = _S2955;
    float3  _S2956 = make_float3 (_S2952, _S2953, _S2954);
    b_0->primal_0 = (*b_0).primal_0;
    b_0->differential_0 = _S2956;
    return;
}

inline __device__ float3  cross_0(float3  left_10, float3  right_10)
{
    float _S2957 = left_10.y;
    float _S2958 = right_10.z;
    float _S2959 = left_10.z;
    float _S2960 = right_10.y;
    float _S2961 = right_10.x;
    float _S2962 = left_10.x;
    return make_float3 (_S2957 * _S2958 - _S2959 * _S2960, _S2959 * _S2961 - _S2962 * _S2958, _S2962 * _S2960 - _S2957 * _S2961);
}

inline __device__ float evaluate_alpha_3dgs(float3  mean_16, Matrix<float, 3, 3>  iscl_rot_0, float opacity_10, float3  ray_o_0, float3  ray_d_0)
{
    float3  grd_0 = mul_0(iscl_rot_0, ray_d_0);
    float3  gcrod_0 = cross_0(grd_0, mul_0(iscl_rot_0, ray_o_0 - mean_16));
    return opacity_10 * (F32_exp((-0.5f * dot_0(gcrod_0, gcrod_0) / dot_0(grd_0, grd_0))));
}

inline __device__ float3  s_primal_ctx_cross_0(float3  _S2963, float3  _S2964)
{
    return cross_0(_S2963, _S2964);
}

inline __device__ void s_bwd_prop_cross_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2965, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2966, float3  _S2967)
{
    _d_cross_0(_S2965, _S2966, _S2967);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_0, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_0, DiffPair_float_0 * dpopacity_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_0, float _s_dOut_8)
{
    float3  _S2968 = (*dpray_o_0).primal_0 - (*dpmean_0).primal_0;
    float3  _S2969 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, _S2968);
    float3  _S2970 = s_primal_ctx_mul_0((*dpiscl_rot_0).primal_0, (*dpray_d_0).primal_0);
    float3  _S2971 = s_primal_ctx_cross_0(_S2970, _S2969);
    float _S2972 = -0.5f * s_primal_ctx_dot_0(_S2971, _S2971);
    float _S2973 = s_primal_ctx_dot_0(_S2970, _S2970);
    float _S2974 = _S2972 / _S2973;
    float _S2975 = _S2973 * _S2973;
    float _S2976 = (*dpopacity_0).primal_0 * _s_dOut_8;
    float _S2977 = s_primal_ctx_exp_1(_S2974) * _s_dOut_8;
    DiffPair_float_0 _S2978;
    (&_S2978)->primal_0 = _S2974;
    (&_S2978)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S2978, _S2976);
    float _S2979 = _S2978.differential_0 / _S2975;
    float _S2980 = _S2972 * - _S2979;
    float _S2981 = _S2973 * _S2979;
    float3  _S2982 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2983;
    (&_S2983)->primal_0 = _S2970;
    (&_S2983)->differential_0 = _S2982;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2984;
    (&_S2984)->primal_0 = _S2970;
    (&_S2984)->differential_0 = _S2982;
    s_bwd_prop_dot_0(&_S2983, &_S2984, _S2980);
    float _S2985 = -0.5f * _S2981;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2986;
    (&_S2986)->primal_0 = _S2971;
    (&_S2986)->differential_0 = _S2982;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2987;
    (&_S2987)->primal_0 = _S2971;
    (&_S2987)->differential_0 = _S2982;
    s_bwd_prop_dot_0(&_S2986, &_S2987, _S2985);
    float3  _S2988 = _S2987.differential_0 + _S2986.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2989;
    (&_S2989)->primal_0 = _S2970;
    (&_S2989)->differential_0 = _S2982;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2990;
    (&_S2990)->primal_0 = _S2969;
    (&_S2990)->differential_0 = _S2982;
    s_bwd_prop_cross_0(&_S2989, &_S2990, _S2988);
    float3  _S2991 = _S2984.differential_0 + _S2983.differential_0 + _S2989.differential_0;
    Matrix<float, 3, 3>  _S2992 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2993;
    (&_S2993)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2993)->differential_0 = _S2992;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2994;
    (&_S2994)->primal_0 = (*dpray_d_0).primal_0;
    (&_S2994)->differential_0 = _S2982;
    s_bwd_prop_mul_1(&_S2993, &_S2994, _S2991);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S2995;
    (&_S2995)->primal_0 = (*dpiscl_rot_0).primal_0;
    (&_S2995)->differential_0 = _S2992;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S2996;
    (&_S2996)->primal_0 = _S2968;
    (&_S2996)->differential_0 = _S2982;
    s_bwd_prop_mul_1(&_S2995, &_S2996, _S2990.differential_0);
    float3  _S2997 = - _S2996.differential_0;
    dpray_d_0->primal_0 = (*dpray_d_0).primal_0;
    dpray_d_0->differential_0 = _S2994.differential_0;
    dpray_o_0->primal_0 = (*dpray_o_0).primal_0;
    dpray_o_0->differential_0 = _S2996.differential_0;
    dpopacity_0->primal_0 = (*dpopacity_0).primal_0;
    dpopacity_0->differential_0 = _S2977;
    Matrix<float, 3, 3>  _S2998 = _S2993.differential_0 + _S2995.differential_0;
    dpiscl_rot_0->primal_0 = (*dpiscl_rot_0).primal_0;
    dpiscl_rot_0->differential_0 = _S2998;
    dpmean_0->primal_0 = (*dpmean_0).primal_0;
    dpmean_0->differential_0 = _S2997;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S2999, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3000, DiffPair_float_0 * _S3001, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3002, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3003, float _S3004)
{
    s_bwd_prop_evaluate_alpha_3dgs_0(_S2999, _S3000, _S3001, _S3002, _S3003, _S3004);
    return;
}

inline __device__ void evaluate_alpha_3dgs_vjp(float3  mean_17, Matrix<float, 3, 3>  iscl_rot_1, float opacity_11, float3  ray_o_1, float3  ray_d_1, float v_alpha_0, float3  * v_mean_5, Matrix<float, 3, 3>  * v_iscl_rot_1, float * v_opacity_5, float3  * v_ray_o_1, float3  * v_ray_d_1)
{
    float3  _S3005 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_0;
    (&dp_mean_0)->primal_0 = mean_17;
    (&dp_mean_0)->differential_0 = _S3005;
    Matrix<float, 3, 3>  _S3006 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_0;
    (&dp_iscl_rot_0)->primal_0 = iscl_rot_1;
    (&dp_iscl_rot_0)->differential_0 = _S3006;
    DiffPair_float_0 dp_opacity_0;
    (&dp_opacity_0)->primal_0 = opacity_11;
    (&dp_opacity_0)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_0;
    (&dp_ray_o_0)->primal_0 = ray_o_1;
    (&dp_ray_o_0)->differential_0 = _S3005;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_0;
    (&dp_ray_d_0)->primal_0 = ray_d_1;
    (&dp_ray_d_0)->differential_0 = _S3005;
    s_bwd_evaluate_alpha_3dgs_0(&dp_mean_0, &dp_iscl_rot_0, &dp_opacity_0, &dp_ray_o_0, &dp_ray_d_0, v_alpha_0);
    *v_mean_5 = dp_mean_0.differential_0;
    *v_iscl_rot_1 = dp_iscl_rot_0.differential_0;
    *v_opacity_5 = dp_opacity_0.differential_0;
    *v_ray_o_1 = dp_ray_o_0.differential_0;
    *v_ray_d_1 = dp_ray_d_0.differential_0;
    return;
}

inline __device__ void evaluate_color_3dgs(float3  mean_18, Matrix<float, 3, 3>  iscl_rot_2, float opacity_12, float3  rgb_10, float3  ray_o_2, float3  ray_d_2, float3  * out_rgb_0, float * depth_10)
{
    *out_rgb_0 = rgb_10;
    float3  _S3007 = mean_18 - ray_o_2;
    *depth_10 = 0.5f * (F32_log((dot_0(_S3007, _S3007) + 9.99999997475242708e-07f)));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpmean_1, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * dpiscl_rot_1, DiffPair_float_0 * dpopacity_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_1, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_1, float3  dpout_rgb_0, float dpdepth_0)
{
    float3  _S3008 = (*dpmean_1).primal_0 - (*dpray_o_1).primal_0;
    float _S3009 = 0.5f * dpdepth_0;
    DiffPair_float_0 _S3010;
    (&_S3010)->primal_0 = s_primal_ctx_dot_0(_S3008, _S3008) + 9.99999997475242708e-07f;
    (&_S3010)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S3010, _S3009);
    float3  _S3011 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3012;
    (&_S3012)->primal_0 = _S3008;
    (&_S3012)->differential_0 = _S3011;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3013;
    (&_S3013)->primal_0 = _S3008;
    (&_S3013)->differential_0 = _S3011;
    s_bwd_prop_dot_0(&_S3012, &_S3013, _S3010.differential_0);
    float3  _S3014 = _S3013.differential_0 + _S3012.differential_0;
    float3  _S3015 = - _S3014;
    dpray_d_1->primal_0 = (*dpray_d_1).primal_0;
    dpray_d_1->differential_0 = _S3011;
    dpray_o_1->primal_0 = (*dpray_o_1).primal_0;
    dpray_o_1->differential_0 = _S3015;
    dprgb_0->primal_0 = (*dprgb_0).primal_0;
    dprgb_0->differential_0 = dpout_rgb_0;
    dpopacity_1->primal_0 = (*dpopacity_1).primal_0;
    dpopacity_1->differential_0 = 0.0f;
    Matrix<float, 3, 3>  _S3016 = makeMatrix<float, 3, 3> (0.0f);
    dpiscl_rot_1->primal_0 = (*dpiscl_rot_1).primal_0;
    dpiscl_rot_1->differential_0 = _S3016;
    dpmean_1->primal_0 = (*dpmean_1).primal_0;
    dpmean_1->differential_0 = _S3014;
    return;
}

inline __device__ void s_bwd_evaluate_color_3dgs_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3017, DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 * _S3018, DiffPair_float_0 * _S3019, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3020, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3021, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3022, float3  _S3023, float _S3024)
{
    s_bwd_prop_evaluate_color_3dgs_0(_S3017, _S3018, _S3019, _S3020, _S3021, _S3022, _S3023, _S3024);
    return;
}

inline __device__ void evaluate_color_3dgs_vjp(float3  mean_19, Matrix<float, 3, 3>  iscl_rot_3, float opacity_13, float3  rgb_11, float3  ray_o_3, float3  ray_d_3, float3  v_out_rgb_0, float v_depth_5, float3  * v_mean_6, Matrix<float, 3, 3>  * v_iscl_rot_2, float * v_opacity_6, float3  * v_rgb_5, float3  * v_ray_o_2, float3  * v_ray_d_2)
{
    float3  _S3025 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_mean_1;
    (&dp_mean_1)->primal_0 = mean_19;
    (&dp_mean_1)->differential_0 = _S3025;
    Matrix<float, 3, 3>  _S3026 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 dp_iscl_rot_1;
    (&dp_iscl_rot_1)->primal_0 = iscl_rot_3;
    (&dp_iscl_rot_1)->differential_0 = _S3026;
    DiffPair_float_0 dp_opacity_1;
    (&dp_opacity_1)->primal_0 = opacity_13;
    (&dp_opacity_1)->differential_0 = 0.0f;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_0;
    (&dp_rgb_0)->primal_0 = rgb_11;
    (&dp_rgb_0)->differential_0 = _S3025;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_1;
    (&dp_ray_o_1)->primal_0 = ray_o_3;
    (&dp_ray_o_1)->differential_0 = _S3025;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_1;
    (&dp_ray_d_1)->primal_0 = ray_d_3;
    (&dp_ray_d_1)->differential_0 = _S3025;
    s_bwd_evaluate_color_3dgs_0(&dp_mean_1, &dp_iscl_rot_1, &dp_opacity_1, &dp_rgb_0, &dp_ray_o_1, &dp_ray_d_1, v_out_rgb_0, v_depth_5);
    *v_mean_6 = dp_mean_1.differential_0;
    *v_iscl_rot_2 = dp_iscl_rot_1.differential_0;
    *v_opacity_6 = dp_opacity_1.differential_0;
    *v_rgb_5 = dp_rgb_0.differential_0;
    *v_ray_o_2 = dp_ray_o_1.differential_0;
    *v_ray_d_2 = dp_ray_d_1.differential_0;
    return;
}

inline __device__ void _d_abs_0(DiffPair_float_0 * dpx_13, float dOut_17)
{
    float _S3027 = _slang_select(((*dpx_13).primal_0) > 0.0f, 1.0f,_slang_select(((*dpx_13).primal_0) == 0.0f, 0.0f,-1.0f)) * dOut_17;
    dpx_13->primal_0 = (*dpx_13).primal_0;
    dpx_13->differential_0 = _S3027;
    return;
}

inline __device__ void _d_exp2_0(DiffPair_float_0 * dpx_14, float dOut_18)
{
    float _S3028 = (F32_exp2(((*dpx_14).primal_0))) * 50.693145751953125f * dOut_18;
    dpx_14->primal_0 = (*dpx_14).primal_0;
    dpx_14->differential_0 = _S3028;
    return;
}

inline __device__ void _d_log_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_15, float3  dOut_19)
{
    float3  _S3029 = make_float3 (1.0f) / (*dpx_15).primal_0 * dOut_19;
    dpx_15->primal_0 = (*dpx_15).primal_0;
    dpx_15->differential_0 = _S3029;
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

inline __device__ void projection_opaque_triangle_persp(float3  mean_20, float4  quat_21, float3  scale_20, float2  hardness_0, FixedArray<float3 , 16>  * sh_coeffs_15, FixedArray<float3 , 2>  * ch_coeffs_0, Matrix<float, 3, 3>  R_21, float3  t_18, float fx_24, float fy_24, float cx_19, float cy_19, FixedArray<float, 10>  * dist_coeffs_28, uint image_width_15, uint image_height_15, float near_plane_10, float far_plane_10, int4  * aabb_xyxy_10, float2  * uv0_0, float2  * uv1_0, float2  * uv2_0, float3  * depth_11, float2  * out_hardness_0, FixedArray<float3 , 3>  * rgb_12, float3  * normal_0)
{
    for(;;)
    {
        float3  mean_c_15 = mul_0(R_21, mean_20) + t_18;
        float _S3030 = mean_c_15.z;
        bool _S3031;
        if(_S3030 < near_plane_10)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = _S3030 > far_plane_10;
        }
        if(_S3031)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3032 = scale_20.x;
        float sx_1 = (F32_exp((_S3032)));
        float _S3033 = scale_20.y;
        float sy_1 = (F32_exp((_S3033)));
        float sz_1 = scale_20.z - 0.5f * (_S3032 + _S3033);
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
        Matrix<float, 3, 3>  _S3034 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_21 + z2_36), 2.0f * (xy_21 + wz_21), 2.0f * (xz_21 - wy_21), 2.0f * (xy_21 - wz_21), 1.0f - 2.0f * (x2_21 + z2_36), 2.0f * (yz_21 + wx_21), 2.0f * (xz_21 + wy_21), 2.0f * (yz_21 - wx_21), 1.0f - 2.0f * (x2_21 + y2_21)));
        float3  vert0_c_0 = mul_0(R_21, mul_0(_S3034, make_float3 (sx_1, 0.0f, 0.0f)) + mean_20) + t_18;
        float3  vert1_c_0 = mul_0(R_21, mul_0(_S3034, make_float3 (sx_1 * (-0.5f + sz_1), sy_1, 0.0f)) + mean_20) + t_18;
        float3  vert2_c_0 = mul_0(R_21, mul_0(_S3034, make_float3 (sx_1 * (-0.5f - sz_1), - sy_1, 0.0f)) + mean_20) + t_18;
        float _S3035 = vert0_c_0.z;
        float _S3036 = vert1_c_0.z;
        float _S3037 = vert2_c_0.z;
        if(_S3035 < near_plane_10)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = _S3035 > far_plane_10;
        }
        if(_S3031)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = _S3036 < near_plane_10;
        }
        if(_S3031)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = _S3036 > far_plane_10;
        }
        if(_S3031)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = _S3037 < near_plane_10;
        }
        if(_S3031)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = _S3037 > far_plane_10;
        }
        if(_S3031)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        for(;;)
        {
            *uv0_0 = float2 {vert0_c_0.x, vert0_c_0.y} / make_float2 (_S3035);
            if(_S3035 < 0.0f)
            {
                _S3031 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3038 = camera_distortion_jac_0(*uv0_0, dist_coeffs_28);
                _S3031 = !((F32_min((determinant_0(_S3038)), ((F32_min((_S3038.rows[int(0)].x), (_S3038.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3031)
            {
                break;
            }
            float u_44 = (*uv0_0).x;
            float v_44 = (*uv0_0).y;
            float r2_44 = u_44 * u_44 + v_44 * v_44;
            float2  _S3039 = *uv0_0 * make_float2 (1.0f + r2_44 * ((*dist_coeffs_28)[int(0)] + r2_44 * ((*dist_coeffs_28)[int(1)] + r2_44 * ((*dist_coeffs_28)[int(2)] + r2_44 * (*dist_coeffs_28)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_28)[int(4)] * u_44 * v_44 + (*dist_coeffs_28)[int(5)] * (r2_44 + 2.0f * u_44 * u_44) + (*dist_coeffs_28)[int(6)] * r2_44, 2.0f * (*dist_coeffs_28)[int(5)] * u_44 * v_44 + (*dist_coeffs_28)[int(4)] * (r2_44 + 2.0f * v_44 * v_44) + (*dist_coeffs_28)[int(7)] * r2_44);
            float2  _S3040 = _S3039 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3039.x + (*dist_coeffs_28)[int(9)] * _S3039.y, 0.0f);
            *uv0_0 = make_float2 (fx_24 * _S3040.x + cx_19, fy_24 * _S3040.y + cy_19);
            break;
        }
        bool all_valid_8 = true & (!_S3031);
        for(;;)
        {
            *uv1_0 = float2 {vert1_c_0.x, vert1_c_0.y} / make_float2 (_S3036);
            if(_S3036 < 0.0f)
            {
                _S3031 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3041 = camera_distortion_jac_0(*uv1_0, dist_coeffs_28);
                _S3031 = !((F32_min((determinant_0(_S3041)), ((F32_min((_S3041.rows[int(0)].x), (_S3041.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3031)
            {
                break;
            }
            float u_45 = (*uv1_0).x;
            float v_45 = (*uv1_0).y;
            float r2_45 = u_45 * u_45 + v_45 * v_45;
            float2  _S3042 = *uv1_0 * make_float2 (1.0f + r2_45 * ((*dist_coeffs_28)[int(0)] + r2_45 * ((*dist_coeffs_28)[int(1)] + r2_45 * ((*dist_coeffs_28)[int(2)] + r2_45 * (*dist_coeffs_28)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_28)[int(4)] * u_45 * v_45 + (*dist_coeffs_28)[int(5)] * (r2_45 + 2.0f * u_45 * u_45) + (*dist_coeffs_28)[int(6)] * r2_45, 2.0f * (*dist_coeffs_28)[int(5)] * u_45 * v_45 + (*dist_coeffs_28)[int(4)] * (r2_45 + 2.0f * v_45 * v_45) + (*dist_coeffs_28)[int(7)] * r2_45);
            float2  _S3043 = _S3042 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3042.x + (*dist_coeffs_28)[int(9)] * _S3042.y, 0.0f);
            *uv1_0 = make_float2 (fx_24 * _S3043.x + cx_19, fy_24 * _S3043.y + cy_19);
            break;
        }
        bool all_valid_9 = all_valid_8 & (!_S3031);
        for(;;)
        {
            *uv2_0 = float2 {vert2_c_0.x, vert2_c_0.y} / make_float2 (_S3037);
            if(_S3037 < 0.0f)
            {
                _S3031 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S3044 = camera_distortion_jac_0(*uv2_0, dist_coeffs_28);
                _S3031 = !((F32_min((determinant_0(_S3044)), ((F32_min((_S3044.rows[int(0)].x), (_S3044.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S3031)
            {
                break;
            }
            float u_46 = (*uv2_0).x;
            float v_46 = (*uv2_0).y;
            float r2_46 = u_46 * u_46 + v_46 * v_46;
            float2  _S3045 = *uv2_0 * make_float2 (1.0f + r2_46 * ((*dist_coeffs_28)[int(0)] + r2_46 * ((*dist_coeffs_28)[int(1)] + r2_46 * ((*dist_coeffs_28)[int(2)] + r2_46 * (*dist_coeffs_28)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_28)[int(4)] * u_46 * v_46 + (*dist_coeffs_28)[int(5)] * (r2_46 + 2.0f * u_46 * u_46) + (*dist_coeffs_28)[int(6)] * r2_46, 2.0f * (*dist_coeffs_28)[int(5)] * u_46 * v_46 + (*dist_coeffs_28)[int(4)] * (r2_46 + 2.0f * v_46 * v_46) + (*dist_coeffs_28)[int(7)] * r2_46);
            float2  _S3046 = _S3045 + make_float2 ((*dist_coeffs_28)[int(8)] * _S3045.x + (*dist_coeffs_28)[int(9)] * _S3045.y, 0.0f);
            *uv2_0 = make_float2 (fx_24 * _S3046.x + cx_19, fy_24 * _S3046.y + cy_19);
            break;
        }
        if(!(all_valid_9 & (!_S3031)))
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
            _S3031 = true;
        }
        else
        {
            _S3031 = xmin_5 >= float(image_width_15);
        }
        if(_S3031)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = ymax_5 <= 0.0f;
        }
        if(_S3031)
        {
            _S3031 = true;
        }
        else
        {
            _S3031 = ymin_5 >= float(image_height_15);
        }
        if(_S3031)
        {
            _S3031 = true;
        }
        else
        {
            if(_S3030 <= 0.0f)
            {
                if(xmin_5 <= 0.0f)
                {
                    _S3031 = xmax_5 >= float(image_width_15);
                }
                else
                {
                    _S3031 = false;
                }
                if(_S3031)
                {
                    _S3031 = true;
                }
                else
                {
                    if(ymin_5 <= 0.0f)
                    {
                        _S3031 = ymax_5 >= float(image_width_15);
                    }
                    else
                    {
                        _S3031 = false;
                    }
                }
            }
            else
            {
                _S3031 = false;
            }
        }
        if(_S3031)
        {
            *aabb_xyxy_10 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_10 = make_int4 (int((F32_floor((xmin_5)))), int((F32_floor((ymin_5)))), int((F32_ceil((xmax_5)))), int((F32_ceil((ymax_5)))));
        *depth_11 = log_0(make_float3 (length_1(vert0_c_0), length_1(vert1_c_0), length_1(vert2_c_0)) + make_float3 (9.999999960041972e-13f));
        *out_hardness_0 = hardness_0;
        float3  _S3047 = mean_20 - - mul_0(transpose_0(R_21), t_18);
        float _S3048 = _S3047.x;
        float _S3049 = _S3047.y;
        float _S3050 = _S3047.z;
        float norm_10 = (F32_sqrt((_S3048 * _S3048 + _S3049 * _S3049 + _S3050 * _S3050)));
        float x_44 = _S3048 / norm_10;
        float y_18 = _S3049 / norm_10;
        float z_15 = _S3050 / norm_10;
        float z2_37 = z_15 * z_15;
        float fTmp0B_15 = -1.09254848957061768f * z_15;
        float fC1_15 = x_44 * x_44 - y_18 * y_18;
        float fS1_15 = 2.0f * x_44 * y_18;
        float fTmp0C_15 = -2.28522896766662598f * z2_37 + 0.4570457935333252f;
        float fTmp1B_15 = 1.44530570507049561f * z_15;
        float3  color_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_15)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_18) * (*sh_coeffs_15)[int(1)] + make_float3 (z_15) * (*sh_coeffs_15)[int(2)] - make_float3 (x_44) * (*sh_coeffs_15)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_15) * (*sh_coeffs_15)[int(4)] + make_float3 (fTmp0B_15 * y_18) * (*sh_coeffs_15)[int(5)] + make_float3 (0.94617468118667603f * z2_37 - 0.31539157032966614f) * (*sh_coeffs_15)[int(6)] + make_float3 (fTmp0B_15 * x_44) * (*sh_coeffs_15)[int(7)] + make_float3 (0.54627424478530884f * fC1_15) * (*sh_coeffs_15)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_44 * fS1_15 + y_18 * fC1_15)) * (*sh_coeffs_15)[int(9)] + make_float3 (fTmp1B_15 * fS1_15) * (*sh_coeffs_15)[int(10)] + make_float3 (fTmp0C_15 * y_18) * (*sh_coeffs_15)[int(11)] + make_float3 (z_15 * (1.86588168144226074f * z2_37 - 1.11952900886535645f)) * (*sh_coeffs_15)[int(12)] + make_float3 (fTmp0C_15 * x_44) * (*sh_coeffs_15)[int(13)] + make_float3 (fTmp1B_15 * fC1_15) * (*sh_coeffs_15)[int(14)] + make_float3 (-0.59004360437393188f * (x_44 * fC1_15 - y_18 * fS1_15)) * (*sh_coeffs_15)[int(15)]);
        float3  _S3051 = make_float3 (0.0f);
        (*rgb_12)[int(0)] = max_0(color_0 + (*ch_coeffs_0)[int(0)] + make_float3 (0.5f), _S3051);
        float3  _S3052 = color_0 - (*ch_coeffs_0)[int(0)] * make_float3 (0.5f);
        float3  _S3053 = (*ch_coeffs_0)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_12)[int(1)] = max_0(_S3052 + _S3053 + make_float3 (0.5f), _S3051);
        (*rgb_12)[int(2)] = max_0(_S3052 - _S3053 + make_float3 (0.5f), _S3051);
        float3  _S3054 = normalize_0(cross_0(vert1_c_0 - vert0_c_0, vert2_c_0 - vert0_c_0));
        *normal_0 = _S3054 * make_float3 (float(- (F32_sign((dot_0(_S3054, mean_c_15))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_fisheye(float3  mean_21, float4  quat_22, float3  scale_21, float2  hardness_1, FixedArray<float3 , 16>  * sh_coeffs_16, FixedArray<float3 , 2>  * ch_coeffs_1, Matrix<float, 3, 3>  R_22, float3  t_19, float fx_25, float fy_25, float cx_20, float cy_20, FixedArray<float, 10>  * dist_coeffs_29, uint image_width_16, uint image_height_16, float near_plane_11, float far_plane_11, int4  * aabb_xyxy_11, float2  * uv0_1, float2  * uv1_1, float2  * uv2_1, float3  * depth_12, float2  * out_hardness_1, FixedArray<float3 , 3>  * rgb_13, float3  * normal_1)
{
    bool _S3055;
    bool _S3056;
    bool _S3057;
    for(;;)
    {
        float3  mean_c_16 = mul_0(R_22, mean_21) + t_19;
        float _S3058 = length_1(mean_c_16);
        bool _S3059;
        if(_S3058 < near_plane_11)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = _S3058 > far_plane_11;
        }
        if(_S3059)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S3060 = scale_21.x;
        float sx_2 = (F32_exp((_S3060)));
        float _S3061 = scale_21.y;
        float sy_2 = (F32_exp((_S3061)));
        float sz_2 = scale_21.z - 0.5f * (_S3060 + _S3061);
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
        Matrix<float, 3, 3>  _S3062 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_22 + z2_38), 2.0f * (xy_22 + wz_22), 2.0f * (xz_22 - wy_22), 2.0f * (xy_22 - wz_22), 1.0f - 2.0f * (x2_22 + z2_38), 2.0f * (yz_22 + wx_22), 2.0f * (xz_22 + wy_22), 2.0f * (yz_22 - wx_22), 1.0f - 2.0f * (x2_22 + y2_22)));
        float3  vert0_c_1 = mul_0(R_22, mul_0(_S3062, make_float3 (sx_2, 0.0f, 0.0f)) + mean_21) + t_19;
        float3  vert1_c_1 = mul_0(R_22, mul_0(_S3062, make_float3 (sx_2 * (-0.5f + sz_2), sy_2, 0.0f)) + mean_21) + t_19;
        float3  vert2_c_1 = mul_0(R_22, mul_0(_S3062, make_float3 (sx_2 * (-0.5f - sz_2), - sy_2, 0.0f)) + mean_21) + t_19;
        float _S3063 = length_1(vert0_c_1);
        float _S3064 = length_1(vert1_c_1);
        float _S3065 = length_1(vert2_c_1);
        if(_S3063 < near_plane_11)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = _S3063 > far_plane_11;
        }
        if(_S3059)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = _S3064 < near_plane_11;
        }
        if(_S3059)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = _S3064 > far_plane_11;
        }
        if(_S3059)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = _S3065 < near_plane_11;
        }
        if(_S3059)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = _S3065 > far_plane_11;
        }
        if(_S3059)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float k_7;
        for(;;)
        {
            float2  _S3066 = float2 {vert0_c_1.x, vert0_c_1.y};
            float r_24 = length_0(_S3066);
            float _S3067 = vert0_c_1.z;
            float theta_20 = (F32_atan2((r_24), (_S3067)));
            if(theta_20 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_20 * theta_20 / 3.0f) / _S3067;
            }
            else
            {
                k_7 = theta_20 / r_24;
            }
            float2  _S3068 = _S3066 * make_float2 (k_7);
            *uv0_1 = _S3068;
            Matrix<float, 2, 2>  _S3069 = camera_distortion_jac_0(_S3068, dist_coeffs_29);
            bool _S3070 = !((F32_min((determinant_0(_S3069)), ((F32_min((_S3069.rows[int(0)].x), (_S3069.rows[int(1)].y)))))) > 0.0f);
            _S3055 = _S3070;
            if(_S3070)
            {
                break;
            }
            float u_47 = (*uv0_1).x;
            float v_47 = (*uv0_1).y;
            float r2_47 = u_47 * u_47 + v_47 * v_47;
            float2  _S3071 = *uv0_1 * make_float2 (1.0f + r2_47 * ((*dist_coeffs_29)[int(0)] + r2_47 * ((*dist_coeffs_29)[int(1)] + r2_47 * ((*dist_coeffs_29)[int(2)] + r2_47 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_47 * v_47 + (*dist_coeffs_29)[int(5)] * (r2_47 + 2.0f * u_47 * u_47) + (*dist_coeffs_29)[int(6)] * r2_47, 2.0f * (*dist_coeffs_29)[int(5)] * u_47 * v_47 + (*dist_coeffs_29)[int(4)] * (r2_47 + 2.0f * v_47 * v_47) + (*dist_coeffs_29)[int(7)] * r2_47);
            float2  _S3072 = _S3071 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3071.x + (*dist_coeffs_29)[int(9)] * _S3071.y, 0.0f);
            *uv0_1 = make_float2 (fx_25 * _S3072.x + cx_20, fy_25 * _S3072.y + cy_20);
            break;
        }
        bool all_valid_10 = true & (!_S3055);
        for(;;)
        {
            float2  _S3073 = float2 {vert1_c_1.x, vert1_c_1.y};
            float r_25 = length_0(_S3073);
            float _S3074 = vert1_c_1.z;
            float theta_21 = (F32_atan2((r_25), (_S3074)));
            if(theta_21 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_21 * theta_21 / 3.0f) / _S3074;
            }
            else
            {
                k_7 = theta_21 / r_25;
            }
            float2  _S3075 = _S3073 * make_float2 (k_7);
            *uv1_1 = _S3075;
            Matrix<float, 2, 2>  _S3076 = camera_distortion_jac_0(_S3075, dist_coeffs_29);
            bool _S3077 = !((F32_min((determinant_0(_S3076)), ((F32_min((_S3076.rows[int(0)].x), (_S3076.rows[int(1)].y)))))) > 0.0f);
            _S3056 = _S3077;
            if(_S3077)
            {
                break;
            }
            float u_48 = (*uv1_1).x;
            float v_48 = (*uv1_1).y;
            float r2_48 = u_48 * u_48 + v_48 * v_48;
            float2  _S3078 = *uv1_1 * make_float2 (1.0f + r2_48 * ((*dist_coeffs_29)[int(0)] + r2_48 * ((*dist_coeffs_29)[int(1)] + r2_48 * ((*dist_coeffs_29)[int(2)] + r2_48 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_48 * v_48 + (*dist_coeffs_29)[int(5)] * (r2_48 + 2.0f * u_48 * u_48) + (*dist_coeffs_29)[int(6)] * r2_48, 2.0f * (*dist_coeffs_29)[int(5)] * u_48 * v_48 + (*dist_coeffs_29)[int(4)] * (r2_48 + 2.0f * v_48 * v_48) + (*dist_coeffs_29)[int(7)] * r2_48);
            float2  _S3079 = _S3078 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3078.x + (*dist_coeffs_29)[int(9)] * _S3078.y, 0.0f);
            *uv1_1 = make_float2 (fx_25 * _S3079.x + cx_20, fy_25 * _S3079.y + cy_20);
            break;
        }
        bool all_valid_11 = all_valid_10 & (!_S3056);
        for(;;)
        {
            float2  _S3080 = float2 {vert2_c_1.x, vert2_c_1.y};
            float r_26 = length_0(_S3080);
            float _S3081 = vert2_c_1.z;
            float theta_22 = (F32_atan2((r_26), (_S3081)));
            if(theta_22 < 0.00100000004749745f)
            {
                k_7 = (1.0f - theta_22 * theta_22 / 3.0f) / _S3081;
            }
            else
            {
                k_7 = theta_22 / r_26;
            }
            float2  _S3082 = _S3080 * make_float2 (k_7);
            *uv2_1 = _S3082;
            Matrix<float, 2, 2>  _S3083 = camera_distortion_jac_0(_S3082, dist_coeffs_29);
            bool _S3084 = !((F32_min((determinant_0(_S3083)), ((F32_min((_S3083.rows[int(0)].x), (_S3083.rows[int(1)].y)))))) > 0.0f);
            _S3057 = _S3084;
            if(_S3084)
            {
                break;
            }
            float u_49 = (*uv2_1).x;
            float v_49 = (*uv2_1).y;
            float r2_49 = u_49 * u_49 + v_49 * v_49;
            float2  _S3085 = *uv2_1 * make_float2 (1.0f + r2_49 * ((*dist_coeffs_29)[int(0)] + r2_49 * ((*dist_coeffs_29)[int(1)] + r2_49 * ((*dist_coeffs_29)[int(2)] + r2_49 * (*dist_coeffs_29)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_29)[int(4)] * u_49 * v_49 + (*dist_coeffs_29)[int(5)] * (r2_49 + 2.0f * u_49 * u_49) + (*dist_coeffs_29)[int(6)] * r2_49, 2.0f * (*dist_coeffs_29)[int(5)] * u_49 * v_49 + (*dist_coeffs_29)[int(4)] * (r2_49 + 2.0f * v_49 * v_49) + (*dist_coeffs_29)[int(7)] * r2_49);
            float2  _S3086 = _S3085 + make_float2 ((*dist_coeffs_29)[int(8)] * _S3085.x + (*dist_coeffs_29)[int(9)] * _S3085.y, 0.0f);
            *uv2_1 = make_float2 (fx_25 * _S3086.x + cx_20, fy_25 * _S3086.y + cy_20);
            break;
        }
        if(!(all_valid_11 & (!_S3057)))
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
            _S3059 = true;
        }
        else
        {
            _S3059 = xmin_6 >= float(image_width_16);
        }
        if(_S3059)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = ymax_6 <= 0.0f;
        }
        if(_S3059)
        {
            _S3059 = true;
        }
        else
        {
            _S3059 = ymin_6 >= float(image_height_16);
        }
        if(_S3059)
        {
            _S3059 = true;
        }
        else
        {
            if((mean_c_16.z) <= 0.0f)
            {
                if(xmin_6 <= 0.0f)
                {
                    _S3059 = xmax_6 >= float(image_width_16);
                }
                else
                {
                    _S3059 = false;
                }
                if(_S3059)
                {
                    _S3059 = true;
                }
                else
                {
                    if(ymin_6 <= 0.0f)
                    {
                        _S3059 = ymax_6 >= float(image_width_16);
                    }
                    else
                    {
                        _S3059 = false;
                    }
                }
            }
            else
            {
                _S3059 = false;
            }
        }
        if(_S3059)
        {
            *aabb_xyxy_11 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_11 = make_int4 (int((F32_floor((xmin_6)))), int((F32_floor((ymin_6)))), int((F32_ceil((xmax_6)))), int((F32_ceil((ymax_6)))));
        *depth_12 = log_0(make_float3 (_S3063, _S3064, _S3065) + make_float3 (9.999999960041972e-13f));
        *out_hardness_1 = hardness_1;
        float3  _S3087 = mean_21 - - mul_0(transpose_0(R_22), t_19);
        float _S3088 = _S3087.x;
        float _S3089 = _S3087.y;
        float _S3090 = _S3087.z;
        float norm_11 = (F32_sqrt((_S3088 * _S3088 + _S3089 * _S3089 + _S3090 * _S3090)));
        float x_46 = _S3088 / norm_11;
        float y_19 = _S3089 / norm_11;
        float z_16 = _S3090 / norm_11;
        float z2_39 = z_16 * z_16;
        float fTmp0B_16 = -1.09254848957061768f * z_16;
        float fC1_16 = x_46 * x_46 - y_19 * y_19;
        float fS1_16 = 2.0f * x_46 * y_19;
        float fTmp0C_16 = -2.28522896766662598f * z2_39 + 0.4570457935333252f;
        float fTmp1B_16 = 1.44530570507049561f * z_16;
        float3  color_1 = make_float3 (0.282094806432724f) * (*sh_coeffs_16)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_19) * (*sh_coeffs_16)[int(1)] + make_float3 (z_16) * (*sh_coeffs_16)[int(2)] - make_float3 (x_46) * (*sh_coeffs_16)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_16) * (*sh_coeffs_16)[int(4)] + make_float3 (fTmp0B_16 * y_19) * (*sh_coeffs_16)[int(5)] + make_float3 (0.94617468118667603f * z2_39 - 0.31539157032966614f) * (*sh_coeffs_16)[int(6)] + make_float3 (fTmp0B_16 * x_46) * (*sh_coeffs_16)[int(7)] + make_float3 (0.54627424478530884f * fC1_16) * (*sh_coeffs_16)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_46 * fS1_16 + y_19 * fC1_16)) * (*sh_coeffs_16)[int(9)] + make_float3 (fTmp1B_16 * fS1_16) * (*sh_coeffs_16)[int(10)] + make_float3 (fTmp0C_16 * y_19) * (*sh_coeffs_16)[int(11)] + make_float3 (z_16 * (1.86588168144226074f * z2_39 - 1.11952900886535645f)) * (*sh_coeffs_16)[int(12)] + make_float3 (fTmp0C_16 * x_46) * (*sh_coeffs_16)[int(13)] + make_float3 (fTmp1B_16 * fC1_16) * (*sh_coeffs_16)[int(14)] + make_float3 (-0.59004360437393188f * (x_46 * fC1_16 - y_19 * fS1_16)) * (*sh_coeffs_16)[int(15)]);
        float3  _S3091 = make_float3 (0.0f);
        (*rgb_13)[int(0)] = max_0(color_1 + (*ch_coeffs_1)[int(0)] + make_float3 (0.5f), _S3091);
        float3  _S3092 = color_1 - (*ch_coeffs_1)[int(0)] * make_float3 (0.5f);
        float3  _S3093 = (*ch_coeffs_1)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgb_13)[int(1)] = max_0(_S3092 + _S3093 + make_float3 (0.5f), _S3091);
        (*rgb_13)[int(2)] = max_0(_S3092 - _S3093 + make_float3 (0.5f), _S3091);
        float3  _S3094 = normalize_0(cross_0(vert1_c_1 - vert0_c_1, vert2_c_1 - vert0_c_1));
        *normal_1 = _S3094 * make_float3 (float(- (F32_sign((dot_0(_S3094, mean_c_16))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_persp_differentiable(float3  mean_22, float4  quat_23, float3  scale_22, float2  hardness_2, FixedArray<float3 , 16>  * sh_coeffs_17, FixedArray<float3 , 2>  * ch_coeffs_2, Matrix<float, 3, 3>  R_23, float3  t_20, float fx_26, float fy_26, float cx_21, float cy_21, FixedArray<float, 10>  * dist_coeffs_30, uint image_width_17, uint image_height_17, float near_plane_12, float far_plane_12, int4  * aabb_xyxy_12, float2  * uv0_2, float2  * uv1_2, float2  * uv2_2, float3  * depth_13, float2  * out_hardness_2, FixedArray<float3 , 3>  * rgb_14, float3  * normal_2)
{
    float3  mean_c_17 = mul_0(R_23, mean_22) + t_20;
    float _S3095 = scale_22.x;
    float sx_3 = (F32_exp((_S3095)));
    float _S3096 = scale_22.y;
    float sy_3 = (F32_exp((_S3096)));
    float sz_3 = scale_22.z - 0.5f * (_S3095 + _S3096);
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
    Matrix<float, 3, 3>  _S3097 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_23 + z2_40), 2.0f * (xy_23 + wz_23), 2.0f * (xz_23 - wy_23), 2.0f * (xy_23 - wz_23), 1.0f - 2.0f * (x2_23 + z2_40), 2.0f * (yz_23 + wx_23), 2.0f * (xz_23 + wy_23), 2.0f * (yz_23 - wx_23), 1.0f - 2.0f * (x2_23 + y2_23)));
    float3  vert0_c_2 = mul_0(R_23, mul_0(_S3097, make_float3 (sx_3, 0.0f, 0.0f)) + mean_22) + t_20;
    float3  vert1_c_2 = mul_0(R_23, mul_0(_S3097, make_float3 (sx_3 * (-0.5f + sz_3), sy_3, 0.0f)) + mean_22) + t_20;
    float3  vert2_c_2 = mul_0(R_23, mul_0(_S3097, make_float3 (sx_3 * (-0.5f - sz_3), - sy_3, 0.0f)) + mean_22) + t_20;
    float2  _S3098 = float2 {vert0_c_2.x, vert0_c_2.y} / make_float2 (vert0_c_2.z);
    float u_50 = _S3098.x;
    float v_50 = _S3098.y;
    float r2_50 = u_50 * u_50 + v_50 * v_50;
    float _S3099 = 2.0f * (*dist_coeffs_30)[int(4)];
    float _S3100 = 2.0f * (*dist_coeffs_30)[int(5)];
    float2  _S3101 = _S3098 * make_float2 (1.0f + r2_50 * ((*dist_coeffs_30)[int(0)] + r2_50 * ((*dist_coeffs_30)[int(1)] + r2_50 * ((*dist_coeffs_30)[int(2)] + r2_50 * (*dist_coeffs_30)[int(3)])))) + make_float2 (_S3099 * u_50 * v_50 + (*dist_coeffs_30)[int(5)] * (r2_50 + 2.0f * u_50 * u_50) + (*dist_coeffs_30)[int(6)] * r2_50, _S3100 * u_50 * v_50 + (*dist_coeffs_30)[int(4)] * (r2_50 + 2.0f * v_50 * v_50) + (*dist_coeffs_30)[int(7)] * r2_50);
    float2  _S3102 = _S3101 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3101.x + (*dist_coeffs_30)[int(9)] * _S3101.y, 0.0f);
    *uv0_2 = make_float2 (fx_26 * _S3102.x + cx_21, fy_26 * _S3102.y + cy_21);
    float2  _S3103 = float2 {vert1_c_2.x, vert1_c_2.y} / make_float2 (vert1_c_2.z);
    float u_51 = _S3103.x;
    float v_51 = _S3103.y;
    float r2_51 = u_51 * u_51 + v_51 * v_51;
    float2  _S3104 = _S3103 * make_float2 (1.0f + r2_51 * ((*dist_coeffs_30)[int(0)] + r2_51 * ((*dist_coeffs_30)[int(1)] + r2_51 * ((*dist_coeffs_30)[int(2)] + r2_51 * (*dist_coeffs_30)[int(3)])))) + make_float2 (_S3099 * u_51 * v_51 + (*dist_coeffs_30)[int(5)] * (r2_51 + 2.0f * u_51 * u_51) + (*dist_coeffs_30)[int(6)] * r2_51, _S3100 * u_51 * v_51 + (*dist_coeffs_30)[int(4)] * (r2_51 + 2.0f * v_51 * v_51) + (*dist_coeffs_30)[int(7)] * r2_51);
    float2  _S3105 = _S3104 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3104.x + (*dist_coeffs_30)[int(9)] * _S3104.y, 0.0f);
    *uv1_2 = make_float2 (fx_26 * _S3105.x + cx_21, fy_26 * _S3105.y + cy_21);
    float2  _S3106 = float2 {vert2_c_2.x, vert2_c_2.y} / make_float2 (vert2_c_2.z);
    float u_52 = _S3106.x;
    float v_52 = _S3106.y;
    float r2_52 = u_52 * u_52 + v_52 * v_52;
    float2  _S3107 = _S3106 * make_float2 (1.0f + r2_52 * ((*dist_coeffs_30)[int(0)] + r2_52 * ((*dist_coeffs_30)[int(1)] + r2_52 * ((*dist_coeffs_30)[int(2)] + r2_52 * (*dist_coeffs_30)[int(3)])))) + make_float2 (_S3099 * u_52 * v_52 + (*dist_coeffs_30)[int(5)] * (r2_52 + 2.0f * u_52 * u_52) + (*dist_coeffs_30)[int(6)] * r2_52, _S3100 * u_52 * v_52 + (*dist_coeffs_30)[int(4)] * (r2_52 + 2.0f * v_52 * v_52) + (*dist_coeffs_30)[int(7)] * r2_52);
    float2  _S3108 = _S3107 + make_float2 ((*dist_coeffs_30)[int(8)] * _S3107.x + (*dist_coeffs_30)[int(9)] * _S3107.y, 0.0f);
    float _S3109 = fx_26 * _S3108.x + cx_21;
    float _S3110 = fy_26 * _S3108.y + cy_21;
    float2  _S3111 = make_float2 (_S3109, _S3110);
    *uv2_2 = _S3111;
    float2  e0_2 = *uv1_2 - *uv0_2;
    float2  e1_2 = _S3111 - *uv1_2;
    float offset_2 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_2.y))))) - 1.0f) * ((F32_abs((e0_2.x * e1_2.y - e0_2.y * e1_2.x))) / (length_0(e0_2) + length_0(e1_2) + length_0(*uv0_2 - _S3111)));
    *aabb_xyxy_12 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_2).x), ((*uv1_2).x)))), (_S3109))) - offset_2)))), int((F32_floor(((F32_min(((F32_min(((*uv0_2).y), ((*uv1_2).y)))), (_S3110))) - offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).x), ((*uv1_2).x)))), (_S3109))) + offset_2)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_2).y), ((*uv1_2).y)))), (_S3110))) + offset_2)))));
    *depth_13 = log_0(make_float3 (length_1(vert0_c_2), length_1(vert1_c_2), length_1(vert2_c_2)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_2 = hardness_2;
    float3  _S3112 = mean_22 - - mul_0(transpose_0(R_23), t_20);
    float _S3113 = _S3112.x;
    float _S3114 = _S3112.y;
    float _S3115 = _S3112.z;
    float norm_12 = (F32_sqrt((_S3113 * _S3113 + _S3114 * _S3114 + _S3115 * _S3115)));
    float x_48 = _S3113 / norm_12;
    float y_20 = _S3114 / norm_12;
    float z_17 = _S3115 / norm_12;
    float z2_41 = z_17 * z_17;
    float fTmp0B_17 = -1.09254848957061768f * z_17;
    float fC1_17 = x_48 * x_48 - y_20 * y_20;
    float fS1_17 = 2.0f * x_48 * y_20;
    float fTmp0C_17 = -2.28522896766662598f * z2_41 + 0.4570457935333252f;
    float fTmp1B_17 = 1.44530570507049561f * z_17;
    float3  color_2 = make_float3 (0.282094806432724f) * (*sh_coeffs_17)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_20) * (*sh_coeffs_17)[int(1)] + make_float3 (z_17) * (*sh_coeffs_17)[int(2)] - make_float3 (x_48) * (*sh_coeffs_17)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_17) * (*sh_coeffs_17)[int(4)] + make_float3 (fTmp0B_17 * y_20) * (*sh_coeffs_17)[int(5)] + make_float3 (0.94617468118667603f * z2_41 - 0.31539157032966614f) * (*sh_coeffs_17)[int(6)] + make_float3 (fTmp0B_17 * x_48) * (*sh_coeffs_17)[int(7)] + make_float3 (0.54627424478530884f * fC1_17) * (*sh_coeffs_17)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_48 * fS1_17 + y_20 * fC1_17)) * (*sh_coeffs_17)[int(9)] + make_float3 (fTmp1B_17 * fS1_17) * (*sh_coeffs_17)[int(10)] + make_float3 (fTmp0C_17 * y_20) * (*sh_coeffs_17)[int(11)] + make_float3 (z_17 * (1.86588168144226074f * z2_41 - 1.11952900886535645f)) * (*sh_coeffs_17)[int(12)] + make_float3 (fTmp0C_17 * x_48) * (*sh_coeffs_17)[int(13)] + make_float3 (fTmp1B_17 * fC1_17) * (*sh_coeffs_17)[int(14)] + make_float3 (-0.59004360437393188f * (x_48 * fC1_17 - y_20 * fS1_17)) * (*sh_coeffs_17)[int(15)]);
    float3  _S3116 = make_float3 (0.0f);
    (*rgb_14)[int(0)] = max_0(color_2 + (*ch_coeffs_2)[int(0)] + make_float3 (0.5f), _S3116);
    float3  _S3117 = color_2 - (*ch_coeffs_2)[int(0)] * make_float3 (0.5f);
    float3  _S3118 = (*ch_coeffs_2)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_14)[int(1)] = max_0(_S3117 + _S3118 + make_float3 (0.5f), _S3116);
    (*rgb_14)[int(2)] = max_0(_S3117 - _S3118 + make_float3 (0.5f), _S3116);
    float3  _S3119 = normalize_0(cross_0(vert1_c_2 - vert0_c_2, vert2_c_2 - vert0_c_2));
    *normal_2 = _S3119 * make_float3 (float(- (F32_sign((dot_0(_S3119, mean_c_17))))));
    return;
}

inline __device__ void _projection_opaque_triangle_fisheye_differentiable(float3  mean_23, float4  quat_24, float3  scale_23, float2  hardness_3, FixedArray<float3 , 16>  * sh_coeffs_18, FixedArray<float3 , 2>  * ch_coeffs_3, Matrix<float, 3, 3>  R_24, float3  t_21, float fx_27, float fy_27, float cx_22, float cy_22, FixedArray<float, 10>  * dist_coeffs_31, uint image_width_18, uint image_height_18, float near_plane_13, float far_plane_13, int4  * aabb_xyxy_13, float2  * uv0_3, float2  * uv1_3, float2  * uv2_3, float3  * depth_14, float2  * out_hardness_3, FixedArray<float3 , 3>  * rgb_15, float3  * normal_3)
{
    float3  mean_c_18 = mul_0(R_24, mean_23) + t_21;
    float _S3120 = scale_23.x;
    float sx_4 = (F32_exp((_S3120)));
    float _S3121 = scale_23.y;
    float sy_4 = (F32_exp((_S3121)));
    float sz_4 = scale_23.z - 0.5f * (_S3120 + _S3121);
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
    Matrix<float, 3, 3>  _S3122 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_24 + z2_42), 2.0f * (xy_24 + wz_24), 2.0f * (xz_24 - wy_24), 2.0f * (xy_24 - wz_24), 1.0f - 2.0f * (x2_24 + z2_42), 2.0f * (yz_24 + wx_24), 2.0f * (xz_24 + wy_24), 2.0f * (yz_24 - wx_24), 1.0f - 2.0f * (x2_24 + y2_24)));
    float3  vert0_c_3 = mul_0(R_24, mul_0(_S3122, make_float3 (sx_4, 0.0f, 0.0f)) + mean_23) + t_21;
    float3  vert1_c_3 = mul_0(R_24, mul_0(_S3122, make_float3 (sx_4 * (-0.5f + sz_4), sy_4, 0.0f)) + mean_23) + t_21;
    float3  vert2_c_3 = mul_0(R_24, mul_0(_S3122, make_float3 (sx_4 * (-0.5f - sz_4), - sy_4, 0.0f)) + mean_23) + t_21;
    float2  _S3123 = float2 {vert0_c_3.x, vert0_c_3.y};
    float r_27 = length_0(_S3123);
    float _S3124 = vert0_c_3.z;
    float theta_23 = (F32_atan2((r_27), (_S3124)));
    float k_8;
    if(theta_23 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_23 * theta_23 / 3.0f) / _S3124;
    }
    else
    {
        k_8 = theta_23 / r_27;
    }
    float2  _S3125 = _S3123 * make_float2 (k_8);
    float u_53 = _S3125.x;
    float v_53 = _S3125.y;
    float r2_53 = u_53 * u_53 + v_53 * v_53;
    float _S3126 = 2.0f * (*dist_coeffs_31)[int(4)];
    float _S3127 = 2.0f * (*dist_coeffs_31)[int(5)];
    float2  _S3128 = _S3125 * make_float2 (1.0f + r2_53 * ((*dist_coeffs_31)[int(0)] + r2_53 * ((*dist_coeffs_31)[int(1)] + r2_53 * ((*dist_coeffs_31)[int(2)] + r2_53 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3126 * u_53 * v_53 + (*dist_coeffs_31)[int(5)] * (r2_53 + 2.0f * u_53 * u_53) + (*dist_coeffs_31)[int(6)] * r2_53, _S3127 * u_53 * v_53 + (*dist_coeffs_31)[int(4)] * (r2_53 + 2.0f * v_53 * v_53) + (*dist_coeffs_31)[int(7)] * r2_53);
    float2  _S3129 = _S3128 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3128.x + (*dist_coeffs_31)[int(9)] * _S3128.y, 0.0f);
    *uv0_3 = make_float2 (fx_27 * _S3129.x + cx_22, fy_27 * _S3129.y + cy_22);
    float2  _S3130 = float2 {vert1_c_3.x, vert1_c_3.y};
    float r_28 = length_0(_S3130);
    float _S3131 = vert1_c_3.z;
    float theta_24 = (F32_atan2((r_28), (_S3131)));
    if(theta_24 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_24 * theta_24 / 3.0f) / _S3131;
    }
    else
    {
        k_8 = theta_24 / r_28;
    }
    float2  _S3132 = _S3130 * make_float2 (k_8);
    float u_54 = _S3132.x;
    float v_54 = _S3132.y;
    float r2_54 = u_54 * u_54 + v_54 * v_54;
    float2  _S3133 = _S3132 * make_float2 (1.0f + r2_54 * ((*dist_coeffs_31)[int(0)] + r2_54 * ((*dist_coeffs_31)[int(1)] + r2_54 * ((*dist_coeffs_31)[int(2)] + r2_54 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3126 * u_54 * v_54 + (*dist_coeffs_31)[int(5)] * (r2_54 + 2.0f * u_54 * u_54) + (*dist_coeffs_31)[int(6)] * r2_54, _S3127 * u_54 * v_54 + (*dist_coeffs_31)[int(4)] * (r2_54 + 2.0f * v_54 * v_54) + (*dist_coeffs_31)[int(7)] * r2_54);
    float2  _S3134 = _S3133 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3133.x + (*dist_coeffs_31)[int(9)] * _S3133.y, 0.0f);
    *uv1_3 = make_float2 (fx_27 * _S3134.x + cx_22, fy_27 * _S3134.y + cy_22);
    float2  _S3135 = float2 {vert2_c_3.x, vert2_c_3.y};
    float r_29 = length_0(_S3135);
    float _S3136 = vert2_c_3.z;
    float theta_25 = (F32_atan2((r_29), (_S3136)));
    if(theta_25 < 0.00100000004749745f)
    {
        k_8 = (1.0f - theta_25 * theta_25 / 3.0f) / _S3136;
    }
    else
    {
        k_8 = theta_25 / r_29;
    }
    float2  _S3137 = _S3135 * make_float2 (k_8);
    float u_55 = _S3137.x;
    float v_55 = _S3137.y;
    float r2_55 = u_55 * u_55 + v_55 * v_55;
    float2  _S3138 = _S3137 * make_float2 (1.0f + r2_55 * ((*dist_coeffs_31)[int(0)] + r2_55 * ((*dist_coeffs_31)[int(1)] + r2_55 * ((*dist_coeffs_31)[int(2)] + r2_55 * (*dist_coeffs_31)[int(3)])))) + make_float2 (_S3126 * u_55 * v_55 + (*dist_coeffs_31)[int(5)] * (r2_55 + 2.0f * u_55 * u_55) + (*dist_coeffs_31)[int(6)] * r2_55, _S3127 * u_55 * v_55 + (*dist_coeffs_31)[int(4)] * (r2_55 + 2.0f * v_55 * v_55) + (*dist_coeffs_31)[int(7)] * r2_55);
    float2  _S3139 = _S3138 + make_float2 ((*dist_coeffs_31)[int(8)] * _S3138.x + (*dist_coeffs_31)[int(9)] * _S3138.y, 0.0f);
    float _S3140 = fx_27 * _S3139.x + cx_22;
    float _S3141 = fy_27 * _S3139.y + cy_22;
    float2  _S3142 = make_float2 (_S3140, _S3141);
    *uv2_3 = _S3142;
    float2  e0_3 = *uv1_3 - *uv0_3;
    float2  e1_3 = _S3142 - *uv1_3;
    float offset_3 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_3.y))))) - 1.0f) * ((F32_abs((e0_3.x * e1_3.y - e0_3.y * e1_3.x))) / (length_0(e0_3) + length_0(e1_3) + length_0(*uv0_3 - _S3142)));
    *aabb_xyxy_13 = make_int4 (int((F32_floor(((F32_min(((F32_min(((*uv0_3).x), ((*uv1_3).x)))), (_S3140))) - offset_3)))), int((F32_floor(((F32_min(((F32_min(((*uv0_3).y), ((*uv1_3).y)))), (_S3141))) - offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).x), ((*uv1_3).x)))), (_S3140))) + offset_3)))), int((F32_ceil(((F32_max(((F32_max(((*uv0_3).y), ((*uv1_3).y)))), (_S3141))) + offset_3)))));
    *depth_14 = log_0(make_float3 (length_1(vert0_c_3), length_1(vert1_c_3), length_1(vert2_c_3)) + make_float3 (9.999999960041972e-13f));
    *out_hardness_3 = hardness_3;
    float3  _S3143 = mean_23 - - mul_0(transpose_0(R_24), t_21);
    float _S3144 = _S3143.x;
    float _S3145 = _S3143.y;
    float _S3146 = _S3143.z;
    float norm_13 = (F32_sqrt((_S3144 * _S3144 + _S3145 * _S3145 + _S3146 * _S3146)));
    float x_50 = _S3144 / norm_13;
    float y_21 = _S3145 / norm_13;
    float z_18 = _S3146 / norm_13;
    float z2_43 = z_18 * z_18;
    float fTmp0B_18 = -1.09254848957061768f * z_18;
    float fC1_18 = x_50 * x_50 - y_21 * y_21;
    float fS1_18 = 2.0f * x_50 * y_21;
    float fTmp0C_18 = -2.28522896766662598f * z2_43 + 0.4570457935333252f;
    float fTmp1B_18 = 1.44530570507049561f * z_18;
    float3  color_3 = make_float3 (0.282094806432724f) * (*sh_coeffs_18)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_21) * (*sh_coeffs_18)[int(1)] + make_float3 (z_18) * (*sh_coeffs_18)[int(2)] - make_float3 (x_50) * (*sh_coeffs_18)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_18) * (*sh_coeffs_18)[int(4)] + make_float3 (fTmp0B_18 * y_21) * (*sh_coeffs_18)[int(5)] + make_float3 (0.94617468118667603f * z2_43 - 0.31539157032966614f) * (*sh_coeffs_18)[int(6)] + make_float3 (fTmp0B_18 * x_50) * (*sh_coeffs_18)[int(7)] + make_float3 (0.54627424478530884f * fC1_18) * (*sh_coeffs_18)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_50 * fS1_18 + y_21 * fC1_18)) * (*sh_coeffs_18)[int(9)] + make_float3 (fTmp1B_18 * fS1_18) * (*sh_coeffs_18)[int(10)] + make_float3 (fTmp0C_18 * y_21) * (*sh_coeffs_18)[int(11)] + make_float3 (z_18 * (1.86588168144226074f * z2_43 - 1.11952900886535645f)) * (*sh_coeffs_18)[int(12)] + make_float3 (fTmp0C_18 * x_50) * (*sh_coeffs_18)[int(13)] + make_float3 (fTmp1B_18 * fC1_18) * (*sh_coeffs_18)[int(14)] + make_float3 (-0.59004360437393188f * (x_50 * fC1_18 - y_21 * fS1_18)) * (*sh_coeffs_18)[int(15)]);
    float3  _S3147 = make_float3 (0.0f);
    (*rgb_15)[int(0)] = max_0(color_3 + (*ch_coeffs_3)[int(0)] + make_float3 (0.5f), _S3147);
    float3  _S3148 = color_3 - (*ch_coeffs_3)[int(0)] * make_float3 (0.5f);
    float3  _S3149 = (*ch_coeffs_3)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgb_15)[int(1)] = max_0(_S3148 + _S3149 + make_float3 (0.5f), _S3147);
    (*rgb_15)[int(2)] = max_0(_S3148 - _S3149 + make_float3 (0.5f), _S3147);
    float3  _S3150 = normalize_0(cross_0(vert1_c_3 - vert0_c_3, vert2_c_3 - vert0_c_3));
    *normal_3 = _S3150 * make_float3 (float(- (F32_sign((dot_0(_S3150, mean_c_18))))));
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_16, float3  _s_dOut_9)
{
    float _S3151 = length_1((*dpx_16).primal_0);
    float3  _S3152 = (*dpx_16).primal_0 * _s_dOut_9;
    float3  _S3153 = make_float3 (1.0f / _S3151) * _s_dOut_9;
    float _S3154 = - ((_S3152.x + _S3152.y + _S3152.z) / (_S3151 * _S3151));
    float3  _S3155 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3156;
    (&_S3156)->primal_0 = (*dpx_16).primal_0;
    (&_S3156)->differential_0 = _S3155;
    s_bwd_length_impl_1(&_S3156, _S3154);
    float3  _S3157 = _S3153 + _S3156.differential_0;
    dpx_16->primal_0 = (*dpx_16).primal_0;
    dpx_16->differential_0 = _S3157;
    return;
}

inline __device__ void s_bwd_normalize_impl_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3158, float3  _S3159)
{
    s_bwd_prop_normalize_impl_0(_S3158, _S3159);
    return;
}

inline __device__ void s_bwd_prop_log_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S3160, float3  _S3161)
{
    _d_log_vector_0(_S3160, _S3161);
    return;
}

inline __device__ void s_bwd_prop_exp2_0(DiffPair_float_0 * _S3162, float _S3163)
{
    _d_exp2_0(_S3162, _S3163);
    return;
}

inline __device__ void s_bwd_prop_abs_0(DiffPair_float_0 * _S3164, float _S3165)
{
    _d_abs_0(_S3164, _S3165);
    return;
}

inline __device__ void projection_opaque_triangle_persp_vjp(float3  mean_24, float4  quat_25, float3  scale_24, float2  hardness_4, FixedArray<float3 , 16>  * sh_coeffs_19, FixedArray<float3 , 2>  * ch_coeffs_4, Matrix<float, 3, 3>  R_25, float3  t_22, float fx_28, float fy_28, float cx_23, float cy_23, FixedArray<float, 10>  * dist_coeffs_32, uint image_width_19, uint image_height_19, float2  v_uv0_0, float2  v_uv1_0, float2  v_uv2_0, float3  v_depth_6, float2  v_out_hardness_0, FixedArray<float3 , 3>  * v_rgb_6, float3  v_normal_0, float3  * v_mean_7, float4  * v_quat_6, float3  * v_scale_6, float2  * v_hardness_0, FixedArray<float3 , 16>  * v_sh_coeffs_5, FixedArray<float3 , 2>  * v_ch_coeffs_0, Matrix<float, 3, 3>  * v_R_7, float3  * v_t_6)
{
    float3  mean_c_19 = s_primal_ctx_mul_0(R_25, mean_24) + t_22;
    float _S3166 = scale_24.x;
    float _S3167 = s_primal_ctx_exp_1(_S3166);
    float _S3168 = scale_24.y;
    float _S3169 = s_primal_ctx_exp_1(_S3168);
    float sz_5 = scale_24.z - 0.5f * (_S3166 + _S3168);
    float _S3170 = quat_25.y;
    float x2_25 = _S3170 * _S3170;
    float y2_25 = quat_25.z * quat_25.z;
    float z2_44 = quat_25.w * quat_25.w;
    float xy_25 = quat_25.y * quat_25.z;
    float xz_25 = quat_25.y * quat_25.w;
    float yz_25 = quat_25.z * quat_25.w;
    float wx_25 = quat_25.x * quat_25.y;
    float wy_25 = quat_25.x * quat_25.z;
    float wz_25 = quat_25.x * quat_25.w;
    Matrix<float, 3, 3>  _S3171 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_25 + z2_44), 2.0f * (xy_25 + wz_25), 2.0f * (xz_25 - wy_25), 2.0f * (xy_25 - wz_25), 1.0f - 2.0f * (x2_25 + z2_44), 2.0f * (yz_25 + wx_25), 2.0f * (xz_25 + wy_25), 2.0f * (yz_25 - wx_25), 1.0f - 2.0f * (x2_25 + y2_25)));
    float3  _S3172 = make_float3 (_S3167, 0.0f, 0.0f);
    float3  vert0_1 = s_primal_ctx_mul_0(_S3171, _S3172) + mean_24;
    float _S3173 = -0.5f + sz_5;
    float3  _S3174 = make_float3 (_S3167 * _S3173, _S3169, 0.0f);
    float3  vert1_1 = s_primal_ctx_mul_0(_S3171, _S3174) + mean_24;
    float _S3175 = -0.5f - sz_5;
    float3  _S3176 = make_float3 (_S3167 * _S3175, - _S3169, 0.0f);
    float3  vert2_1 = s_primal_ctx_mul_0(_S3171, _S3176) + mean_24;
    float3  vert0_c_4 = s_primal_ctx_mul_0(R_25, vert0_1) + t_22;
    float3  vert1_c_4 = s_primal_ctx_mul_0(R_25, vert1_1) + t_22;
    float3  vert2_c_4 = s_primal_ctx_mul_0(R_25, vert2_1) + t_22;
    float2  _S3177 = float2 {vert0_c_4.x, vert0_c_4.y};
    float _S3178 = vert0_c_4.z;
    float2  _S3179 = make_float2 (_S3178);
    float2  _S3180 = _S3177 / make_float2 (_S3178);
    float2  _S3181 = make_float2 (_S3178 * _S3178);
    float u_56 = _S3180.x;
    float v_56 = _S3180.y;
    float r2_56 = u_56 * u_56 + v_56 * v_56;
    float _S3182 = (*dist_coeffs_32)[int(2)] + r2_56 * (*dist_coeffs_32)[int(3)];
    float _S3183 = (*dist_coeffs_32)[int(1)] + r2_56 * _S3182;
    float _S3184 = (*dist_coeffs_32)[int(0)] + r2_56 * _S3183;
    float radial_1 = 1.0f + r2_56 * _S3184;
    float2  _S3185 = make_float2 (radial_1);
    float _S3186 = 2.0f * (*dist_coeffs_32)[int(4)];
    float _S3187 = _S3186 * u_56;
    float _S3188 = 2.0f * u_56;
    float _S3189 = 2.0f * (*dist_coeffs_32)[int(5)];
    float _S3190 = _S3189 * u_56;
    float _S3191 = 2.0f * v_56;
    float2  _S3192 = _S3180 * make_float2 (radial_1) + make_float2 (_S3187 * v_56 + (*dist_coeffs_32)[int(5)] * (r2_56 + _S3188 * u_56) + (*dist_coeffs_32)[int(6)] * r2_56, _S3190 * v_56 + (*dist_coeffs_32)[int(4)] * (r2_56 + _S3191 * v_56) + (*dist_coeffs_32)[int(7)] * r2_56);
    float2  _S3193 = _S3192 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3192.x + (*dist_coeffs_32)[int(9)] * _S3192.y, 0.0f);
    float _S3194 = fx_28 * _S3193.x + cx_23;
    float _S3195 = fy_28 * _S3193.y + cy_23;
    float2  _S3196 = make_float2 (_S3194, _S3195);
    float2  _S3197 = float2 {vert1_c_4.x, vert1_c_4.y};
    float _S3198 = vert1_c_4.z;
    float2  _S3199 = make_float2 (_S3198);
    float2  _S3200 = _S3197 / make_float2 (_S3198);
    float2  _S3201 = make_float2 (_S3198 * _S3198);
    float u_57 = _S3200.x;
    float v_57 = _S3200.y;
    float r2_57 = u_57 * u_57 + v_57 * v_57;
    float _S3202 = (*dist_coeffs_32)[int(2)] + r2_57 * (*dist_coeffs_32)[int(3)];
    float _S3203 = (*dist_coeffs_32)[int(1)] + r2_57 * _S3202;
    float _S3204 = (*dist_coeffs_32)[int(0)] + r2_57 * _S3203;
    float radial_2 = 1.0f + r2_57 * _S3204;
    float2  _S3205 = make_float2 (radial_2);
    float _S3206 = _S3186 * u_57;
    float _S3207 = 2.0f * u_57;
    float _S3208 = _S3189 * u_57;
    float _S3209 = 2.0f * v_57;
    float2  _S3210 = _S3200 * make_float2 (radial_2) + make_float2 (_S3206 * v_57 + (*dist_coeffs_32)[int(5)] * (r2_57 + _S3207 * u_57) + (*dist_coeffs_32)[int(6)] * r2_57, _S3208 * v_57 + (*dist_coeffs_32)[int(4)] * (r2_57 + _S3209 * v_57) + (*dist_coeffs_32)[int(7)] * r2_57);
    float2  _S3211 = _S3210 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3210.x + (*dist_coeffs_32)[int(9)] * _S3210.y, 0.0f);
    float _S3212 = fx_28 * _S3211.x + cx_23;
    float _S3213 = fy_28 * _S3211.y + cy_23;
    float2  _S3214 = make_float2 (_S3212, _S3213);
    float2  _S3215 = float2 {vert2_c_4.x, vert2_c_4.y};
    float _S3216 = vert2_c_4.z;
    float2  _S3217 = make_float2 (_S3216);
    float2  _S3218 = _S3215 / make_float2 (_S3216);
    float2  _S3219 = make_float2 (_S3216 * _S3216);
    float u_58 = _S3218.x;
    float v_58 = _S3218.y;
    float r2_58 = u_58 * u_58 + v_58 * v_58;
    float _S3220 = (*dist_coeffs_32)[int(2)] + r2_58 * (*dist_coeffs_32)[int(3)];
    float _S3221 = (*dist_coeffs_32)[int(1)] + r2_58 * _S3220;
    float _S3222 = (*dist_coeffs_32)[int(0)] + r2_58 * _S3221;
    float radial_3 = 1.0f + r2_58 * _S3222;
    float2  _S3223 = make_float2 (radial_3);
    float _S3224 = _S3186 * u_58;
    float _S3225 = 2.0f * u_58;
    float _S3226 = _S3189 * u_58;
    float _S3227 = 2.0f * v_58;
    float2  _S3228 = _S3218 * make_float2 (radial_3) + make_float2 (_S3224 * v_58 + (*dist_coeffs_32)[int(5)] * (r2_58 + _S3225 * u_58) + (*dist_coeffs_32)[int(6)] * r2_58, _S3226 * v_58 + (*dist_coeffs_32)[int(4)] * (r2_58 + _S3227 * v_58) + (*dist_coeffs_32)[int(7)] * r2_58);
    float2  _S3229 = _S3228 + make_float2 ((*dist_coeffs_32)[int(8)] * _S3228.x + (*dist_coeffs_32)[int(9)] * _S3228.y, 0.0f);
    float _S3230 = fx_28 * _S3229.x + cx_23;
    float _S3231 = fy_28 * _S3229.y + cy_23;
    float2  _S3232 = make_float2 (_S3230, _S3231);
    float2  e0_4 = _S3214 - _S3196;
    float2  e1_4 = _S3232 - _S3214;
    float2  e2_0 = _S3196 - _S3232;
    float _S3233 = e0_4.x;
    float _S3234 = e1_4.y;
    float _S3235 = e0_4.y;
    float _S3236 = e1_4.x;
    float _S3237 = _S3233 * _S3234 - _S3235 * _S3236;
    float _S3238 = 1.0f - hardness_4.y;
    float _S3239 = -1.0f / _S3238;
    float _S3240 = _S3238 * _S3238;
    float _S3241 = s_primal_ctx_max_0(_S3194, _S3212);
    float _S3242 = s_primal_ctx_min_0(_S3194, _S3212);
    float _S3243 = s_primal_ctx_max_0(_S3195, _S3213);
    float _S3244 = s_primal_ctx_min_0(_S3195, _S3213);
    float3  _S3245 = make_float3 (length_1(vert0_c_4), length_1(vert1_c_4), length_1(vert2_c_4)) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3246 = transpose_0(R_25);
    float3  _S3247 = mean_24 - - s_primal_ctx_mul_0(_S3246, t_22);
    float _S3248 = _S3247.x;
    float _S3249 = _S3247.y;
    float _S3250 = _S3247.z;
    float _S3251 = _S3248 * _S3248 + _S3249 * _S3249 + _S3250 * _S3250;
    float _S3252 = s_primal_ctx_sqrt_0(_S3251);
    float x_51 = _S3248 / _S3252;
    float3  _S3253 = make_float3 (x_51);
    float _S3254 = _S3252 * _S3252;
    float y_22 = _S3249 / _S3252;
    float z_19 = _S3250 / _S3252;
    float3  _S3255 = make_float3 (z_19);
    float _S3256 = - y_22;
    float3  _S3257 = make_float3 (_S3256);
    float z2_45 = z_19 * z_19;
    float fTmp0B_19 = -1.09254848957061768f * z_19;
    float fC1_19 = x_51 * x_51 - y_22 * y_22;
    float _S3258 = 2.0f * x_51;
    float fS1_19 = _S3258 * y_22;
    float pSH6_5 = 0.94617468118667603f * z2_45 - 0.31539157032966614f;
    float3  _S3259 = make_float3 (pSH6_5);
    float pSH7_5 = fTmp0B_19 * x_51;
    float3  _S3260 = make_float3 (pSH7_5);
    float pSH5_5 = fTmp0B_19 * y_22;
    float3  _S3261 = make_float3 (pSH5_5);
    float pSH8_5 = 0.54627424478530884f * fC1_19;
    float3  _S3262 = make_float3 (pSH8_5);
    float pSH4_5 = 0.54627424478530884f * fS1_19;
    float3  _S3263 = make_float3 (pSH4_5);
    float fTmp0C_19 = -2.28522896766662598f * z2_45 + 0.4570457935333252f;
    float fTmp1B_19 = 1.44530570507049561f * z_19;
    float _S3264 = 1.86588168144226074f * z2_45 - 1.11952900886535645f;
    float pSH12_5 = z_19 * _S3264;
    float3  _S3265 = make_float3 (pSH12_5);
    float pSH13_5 = fTmp0C_19 * x_51;
    float3  _S3266 = make_float3 (pSH13_5);
    float pSH11_5 = fTmp0C_19 * y_22;
    float3  _S3267 = make_float3 (pSH11_5);
    float pSH14_5 = fTmp1B_19 * fC1_19;
    float3  _S3268 = make_float3 (pSH14_5);
    float pSH10_5 = fTmp1B_19 * fS1_19;
    float3  _S3269 = make_float3 (pSH10_5);
    float pSH15_5 = -0.59004360437393188f * (x_51 * fC1_19 - y_22 * fS1_19);
    float3  _S3270 = make_float3 (pSH15_5);
    float pSH9_5 = -0.59004360437393188f * (x_51 * fS1_19 + y_22 * fC1_19);
    float3  _S3271 = make_float3 (pSH9_5);
    float3  color_4 = make_float3 (0.282094806432724f) * (*sh_coeffs_19)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3256) * (*sh_coeffs_19)[int(1)] + make_float3 (z_19) * (*sh_coeffs_19)[int(2)] - make_float3 (x_51) * (*sh_coeffs_19)[int(3)]) + (make_float3 (pSH4_5) * (*sh_coeffs_19)[int(4)] + make_float3 (pSH5_5) * (*sh_coeffs_19)[int(5)] + make_float3 (pSH6_5) * (*sh_coeffs_19)[int(6)] + make_float3 (pSH7_5) * (*sh_coeffs_19)[int(7)] + make_float3 (pSH8_5) * (*sh_coeffs_19)[int(8)]) + (make_float3 (pSH9_5) * (*sh_coeffs_19)[int(9)] + make_float3 (pSH10_5) * (*sh_coeffs_19)[int(10)] + make_float3 (pSH11_5) * (*sh_coeffs_19)[int(11)] + make_float3 (pSH12_5) * (*sh_coeffs_19)[int(12)] + make_float3 (pSH13_5) * (*sh_coeffs_19)[int(13)] + make_float3 (pSH14_5) * (*sh_coeffs_19)[int(14)] + make_float3 (pSH15_5) * (*sh_coeffs_19)[int(15)]);
    float3  _S3272 = color_4 + (*ch_coeffs_4)[int(0)] + make_float3 (0.5f);
    float3  _S3273 = make_float3 (0.0f);
    float3  _S3274 = color_4 - (*ch_coeffs_4)[int(0)] * make_float3 (0.5f);
    float _S3275 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3276 = make_float3 (_S3275);
    float3  _S3277 = (*ch_coeffs_4)[int(1)] * make_float3 (_S3275);
    float3  _S3278 = _S3274 + _S3277 + make_float3 (0.5f);
    float3  _S3279 = _S3274 - _S3277 + make_float3 (0.5f);
    float3  _S3280 = vert1_c_4 - vert0_c_4;
    float3  _S3281 = vert2_c_4 - vert0_c_4;
    float3  _S3282 = s_primal_ctx_cross_0(_S3280, _S3281);
    float3  _S3283 = normalize_0(_S3282);
    float3  _S3284 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3283, mean_c_19)))))) * v_normal_0;
    float3  _S3285 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3286;
    (&_S3286)->primal_0 = _S3283;
    (&_S3286)->differential_0 = _S3285;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3287;
    (&_S3287)->primal_0 = mean_c_19;
    (&_S3287)->differential_0 = _S3285;
    s_bwd_prop_dot_0(&_S3286, &_S3287, 0.0f);
    float3  _S3288 = _S3284 + _S3286.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3289;
    (&_S3289)->primal_0 = _S3282;
    (&_S3289)->differential_0 = _S3285;
    s_bwd_normalize_impl_0(&_S3289, _S3288);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3290;
    (&_S3290)->primal_0 = _S3280;
    (&_S3290)->differential_0 = _S3285;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3291;
    (&_S3291)->primal_0 = _S3281;
    (&_S3291)->differential_0 = _S3285;
    s_bwd_prop_cross_0(&_S3290, &_S3291, _S3289.differential_0);
    float3  _S3292 = - _S3291.differential_0;
    float3  _S3293 = - _S3290.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3294;
    (&_S3294)->primal_0 = _S3279;
    (&_S3294)->differential_0 = _S3285;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3295;
    (&_S3295)->primal_0 = _S3273;
    (&_S3295)->differential_0 = _S3285;
    s_bwd_prop_max_0(&_S3294, &_S3295, (*v_rgb_6)[int(2)]);
    float3  _S3296 = - _S3294.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3297;
    (&_S3297)->primal_0 = _S3278;
    (&_S3297)->differential_0 = _S3285;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3298;
    (&_S3298)->primal_0 = _S3273;
    (&_S3298)->differential_0 = _S3285;
    s_bwd_prop_max_0(&_S3297, &_S3298, (*v_rgb_6)[int(1)]);
    float3  _S3299 = _S3276 * (_S3296 + _S3297.differential_0);
    float3  _S3300 = _S3294.differential_0 + _S3297.differential_0;
    float3  _S3301 = make_float3 (0.5f) * - _S3300;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3302;
    (&_S3302)->primal_0 = _S3272;
    (&_S3302)->differential_0 = _S3285;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3303;
    (&_S3303)->primal_0 = _S3273;
    (&_S3303)->differential_0 = _S3285;
    s_bwd_prop_max_0(&_S3302, &_S3303, (*v_rgb_6)[int(0)]);
    float3  _S3304 = _S3301 + _S3302.differential_0;
    float3  _S3305 = _S3300 + _S3302.differential_0;
    float3  _S3306 = _S3270 * _S3305;
    float3  _S3307 = (*sh_coeffs_19)[int(15)] * _S3305;
    float3  _S3308 = _S3268 * _S3305;
    float3  _S3309 = (*sh_coeffs_19)[int(14)] * _S3305;
    float3  _S3310 = _S3266 * _S3305;
    float3  _S3311 = (*sh_coeffs_19)[int(13)] * _S3305;
    float3  _S3312 = _S3265 * _S3305;
    float3  _S3313 = (*sh_coeffs_19)[int(12)] * _S3305;
    float3  _S3314 = _S3267 * _S3305;
    float3  _S3315 = (*sh_coeffs_19)[int(11)] * _S3305;
    float3  _S3316 = _S3269 * _S3305;
    float3  _S3317 = (*sh_coeffs_19)[int(10)] * _S3305;
    float3  _S3318 = _S3271 * _S3305;
    float3  _S3319 = (*sh_coeffs_19)[int(9)] * _S3305;
    float s_diff_fS2_T_5 = -0.59004360437393188f * (_S3319.x + _S3319.y + _S3319.z);
    float s_diff_fC2_T_5 = -0.59004360437393188f * (_S3307.x + _S3307.y + _S3307.z);
    float _S3320 = _S3317.x + _S3317.y + _S3317.z;
    float _S3321 = _S3309.x + _S3309.y + _S3309.z;
    float _S3322 = _S3315.x + _S3315.y + _S3315.z;
    float _S3323 = _S3311.x + _S3311.y + _S3311.z;
    float _S3324 = _S3313.x + _S3313.y + _S3313.z;
    float _S3325 = - s_diff_fC2_T_5;
    float3  _S3326 = _S3262 * _S3305;
    float3  _S3327 = (*sh_coeffs_19)[int(8)] * _S3305;
    float3  _S3328 = _S3260 * _S3305;
    float3  _S3329 = (*sh_coeffs_19)[int(7)] * _S3305;
    float3  _S3330 = _S3259 * _S3305;
    float3  _S3331 = (*sh_coeffs_19)[int(6)] * _S3305;
    float3  _S3332 = _S3261 * _S3305;
    float3  _S3333 = (*sh_coeffs_19)[int(5)] * _S3305;
    float3  _S3334 = _S3263 * _S3305;
    float3  _S3335 = (*sh_coeffs_19)[int(4)] * _S3305;
    float _S3336 = _S3333.x + _S3333.y + _S3333.z;
    float _S3337 = _S3329.x + _S3329.y + _S3329.z;
    float _S3338 = fTmp1B_19 * _S3320 + x_51 * s_diff_fS2_T_5 + y_22 * _S3325 + 0.54627424478530884f * (_S3335.x + _S3335.y + _S3335.z);
    float _S3339 = fTmp1B_19 * _S3321 + y_22 * s_diff_fS2_T_5 + x_51 * s_diff_fC2_T_5 + 0.54627424478530884f * (_S3327.x + _S3327.y + _S3327.z);
    float _S3340 = y_22 * - _S3339;
    float _S3341 = x_51 * _S3339;
    float _S3342 = z_19 * (1.86588168144226074f * (z_19 * _S3324) + -2.28522896766662598f * (y_22 * _S3322 + x_51 * _S3323) + 0.94617468118667603f * (_S3331.x + _S3331.y + _S3331.z));
    float3  _S3343 = make_float3 (0.48860251903533936f) * _S3305;
    float3  _S3344 = - _S3343;
    float3  _S3345 = _S3253 * _S3344;
    float3  _S3346 = (*sh_coeffs_19)[int(3)] * _S3344;
    float3  _S3347 = _S3255 * _S3343;
    float3  _S3348 = (*sh_coeffs_19)[int(2)] * _S3343;
    float3  _S3349 = _S3257 * _S3343;
    float3  _S3350 = (*sh_coeffs_19)[int(1)] * _S3343;
    float _S3351 = (_S3264 * _S3324 + 1.44530570507049561f * (fS1_19 * _S3320 + fC1_19 * _S3321) + -1.09254848957061768f * (y_22 * _S3336 + x_51 * _S3337) + _S3342 + _S3342 + _S3348.x + _S3348.y + _S3348.z) / _S3254;
    float _S3352 = _S3252 * _S3351;
    float _S3353 = (fTmp0C_19 * _S3322 + fC1_19 * s_diff_fS2_T_5 + fS1_19 * _S3325 + fTmp0B_19 * _S3336 + _S3258 * _S3338 + _S3340 + _S3340 + - (_S3350.x + _S3350.y + _S3350.z)) / _S3254;
    float _S3354 = _S3252 * _S3353;
    float _S3355 = (fTmp0C_19 * _S3323 + fS1_19 * s_diff_fS2_T_5 + fC1_19 * s_diff_fC2_T_5 + fTmp0B_19 * _S3337 + 2.0f * (y_22 * _S3338) + _S3341 + _S3341 + _S3346.x + _S3346.y + _S3346.z) / _S3254;
    float _S3356 = _S3252 * _S3355;
    float _S3357 = _S3250 * - _S3351 + _S3249 * - _S3353 + _S3248 * - _S3355;
    DiffPair_float_0 _S3358;
    (&_S3358)->primal_0 = _S3251;
    (&_S3358)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3358, _S3357);
    float _S3359 = _S3250 * _S3358.differential_0;
    float _S3360 = _S3249 * _S3358.differential_0;
    float _S3361 = _S3248 * _S3358.differential_0;
    float3  _S3362 = make_float3 (0.282094806432724f) * _S3305;
    float3  _S3363 = make_float3 (_S3356 + _S3361 + _S3361, _S3354 + _S3360 + _S3360, _S3352 + _S3359 + _S3359);
    float3  _S3364 = - - _S3363;
    Matrix<float, 3, 3>  _S3365 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3366;
    (&_S3366)->primal_0 = _S3246;
    (&_S3366)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3367;
    (&_S3367)->primal_0 = t_22;
    (&_S3367)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3366, &_S3367, _S3364);
    Matrix<float, 3, 3>  _S3368 = transpose_0(_S3366.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3369;
    (&_S3369)->primal_0 = _S3245;
    (&_S3369)->differential_0 = _S3285;
    s_bwd_prop_log_1(&_S3369, v_depth_6);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3370;
    (&_S3370)->primal_0 = vert2_c_4;
    (&_S3370)->differential_0 = _S3285;
    s_bwd_length_impl_1(&_S3370, _S3369.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3371;
    (&_S3371)->primal_0 = vert1_c_4;
    (&_S3371)->differential_0 = _S3285;
    s_bwd_length_impl_1(&_S3371, _S3369.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3372;
    (&_S3372)->primal_0 = vert0_c_4;
    (&_S3372)->differential_0 = _S3285;
    s_bwd_length_impl_1(&_S3372, _S3369.differential_0.x);
    DiffPair_float_0 _S3373;
    (&_S3373)->primal_0 = _S3244;
    (&_S3373)->differential_0 = 0.0f;
    DiffPair_float_0 _S3374;
    (&_S3374)->primal_0 = _S3231;
    (&_S3374)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3373, &_S3374, 0.0f);
    DiffPair_float_0 _S3375;
    (&_S3375)->primal_0 = _S3195;
    (&_S3375)->differential_0 = 0.0f;
    DiffPair_float_0 _S3376;
    (&_S3376)->primal_0 = _S3213;
    (&_S3376)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3375, &_S3376, _S3373.differential_0);
    DiffPair_float_0 _S3377;
    (&_S3377)->primal_0 = _S3243;
    (&_S3377)->differential_0 = 0.0f;
    DiffPair_float_0 _S3378;
    (&_S3378)->primal_0 = _S3231;
    (&_S3378)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3377, &_S3378, 0.0f);
    DiffPair_float_0 _S3379;
    (&_S3379)->primal_0 = _S3195;
    (&_S3379)->differential_0 = 0.0f;
    DiffPair_float_0 _S3380;
    (&_S3380)->primal_0 = _S3213;
    (&_S3380)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3379, &_S3380, _S3377.differential_0);
    DiffPair_float_0 _S3381;
    (&_S3381)->primal_0 = _S3242;
    (&_S3381)->differential_0 = 0.0f;
    DiffPair_float_0 _S3382;
    (&_S3382)->primal_0 = _S3230;
    (&_S3382)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3381, &_S3382, 0.0f);
    DiffPair_float_0 _S3383;
    (&_S3383)->primal_0 = _S3194;
    (&_S3383)->differential_0 = 0.0f;
    DiffPair_float_0 _S3384;
    (&_S3384)->primal_0 = _S3212;
    (&_S3384)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3383, &_S3384, _S3381.differential_0);
    DiffPair_float_0 _S3385;
    (&_S3385)->primal_0 = _S3241;
    (&_S3385)->differential_0 = 0.0f;
    DiffPair_float_0 _S3386;
    (&_S3386)->primal_0 = _S3230;
    (&_S3386)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3385, &_S3386, 0.0f);
    DiffPair_float_0 _S3387;
    (&_S3387)->primal_0 = _S3194;
    (&_S3387)->differential_0 = 0.0f;
    DiffPair_float_0 _S3388;
    (&_S3388)->primal_0 = _S3212;
    (&_S3388)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3387, &_S3388, _S3385.differential_0);
    DiffPair_float_0 _S3389;
    (&_S3389)->primal_0 = _S3239;
    (&_S3389)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3389, 0.0f);
    float _S3390 = - (-1.0f * - (_S3389.differential_0 / _S3240));
    float2  _S3391 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3392;
    (&_S3392)->primal_0 = e2_0;
    (&_S3392)->differential_0 = _S3391;
    s_bwd_length_impl_0(&_S3392, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3393;
    (&_S3393)->primal_0 = e1_4;
    (&_S3393)->differential_0 = _S3391;
    s_bwd_length_impl_0(&_S3393, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3394;
    (&_S3394)->primal_0 = e0_4;
    (&_S3394)->differential_0 = _S3391;
    s_bwd_length_impl_0(&_S3394, -0.0f);
    DiffPair_float_0 _S3395;
    (&_S3395)->primal_0 = _S3237;
    (&_S3395)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3395, 0.0f);
    float _S3396 = - _S3395.differential_0;
    float2  _S3397 = _S3393.differential_0 + make_float2 (_S3235 * _S3396, _S3233 * _S3395.differential_0);
    float2  _S3398 = _S3394.differential_0 + make_float2 (_S3234 * _S3395.differential_0, _S3236 * _S3396);
    float2  _S3399 = v_uv2_0 + - _S3392.differential_0 + _S3397;
    float _S3400 = fx_28 * (_S3382.differential_0 + _S3386.differential_0 + _S3399.x);
    float2  _S3401 = make_float2 (_S3400, fy_28 * (_S3374.differential_0 + _S3378.differential_0 + _S3399.y)) + make_float2 ((*dist_coeffs_32)[int(8)] * _S3400, (*dist_coeffs_32)[int(9)] * _S3400);
    float2  _S3402 = _S3218 * _S3401;
    float _S3403 = (*dist_coeffs_32)[int(4)] * _S3401.y;
    float _S3404 = (*dist_coeffs_32)[int(5)] * _S3401.x;
    float _S3405 = _S3402.x + _S3402.y;
    float _S3406 = r2_58 * _S3405;
    float _S3407 = r2_58 * _S3406;
    float _S3408 = (*dist_coeffs_32)[int(7)] * _S3401.y + _S3403 + (*dist_coeffs_32)[int(6)] * _S3401.x + _S3404 + _S3222 * _S3405 + _S3221 * _S3406 + _S3220 * _S3407 + (*dist_coeffs_32)[int(3)] * (r2_58 * _S3407);
    float _S3409 = v_58 * _S3408;
    float _S3410 = u_58 * _S3408;
    float2  _S3411 = (_S3223 * _S3401 + make_float2 (_S3189 * (v_58 * _S3401.y) + _S3225 * _S3404 + 2.0f * (u_58 * _S3404) + _S3186 * (v_58 * _S3401.x) + _S3410 + _S3410, _S3227 * _S3403 + 2.0f * (v_58 * _S3403) + _S3226 * _S3401.y + _S3224 * _S3401.x + _S3409 + _S3409)) / _S3219;
    float2  _S3412 = _S3215 * - _S3411;
    float2  _S3413 = _S3217 * _S3411;
    float2  _S3414 = v_uv1_0 + - _S3397 + _S3398;
    float _S3415 = fx_28 * (_S3384.differential_0 + _S3388.differential_0 + _S3414.x);
    float2  _S3416 = make_float2 (_S3415, fy_28 * (_S3376.differential_0 + _S3380.differential_0 + _S3414.y)) + make_float2 ((*dist_coeffs_32)[int(8)] * _S3415, (*dist_coeffs_32)[int(9)] * _S3415);
    float2  _S3417 = _S3200 * _S3416;
    float _S3418 = (*dist_coeffs_32)[int(4)] * _S3416.y;
    float _S3419 = (*dist_coeffs_32)[int(5)] * _S3416.x;
    float _S3420 = _S3417.x + _S3417.y;
    float _S3421 = r2_57 * _S3420;
    float _S3422 = r2_57 * _S3421;
    float _S3423 = (*dist_coeffs_32)[int(7)] * _S3416.y + _S3418 + (*dist_coeffs_32)[int(6)] * _S3416.x + _S3419 + _S3204 * _S3420 + _S3203 * _S3421 + _S3202 * _S3422 + (*dist_coeffs_32)[int(3)] * (r2_57 * _S3422);
    float _S3424 = v_57 * _S3423;
    float _S3425 = u_57 * _S3423;
    float2  _S3426 = (_S3205 * _S3416 + make_float2 (_S3189 * (v_57 * _S3416.y) + _S3207 * _S3419 + 2.0f * (u_57 * _S3419) + _S3186 * (v_57 * _S3416.x) + _S3425 + _S3425, _S3209 * _S3418 + 2.0f * (v_57 * _S3418) + _S3208 * _S3416.y + _S3206 * _S3416.x + _S3424 + _S3424)) / _S3201;
    float2  _S3427 = _S3197 * - _S3426;
    float2  _S3428 = _S3199 * _S3426;
    float _S3429 = _S3427.x + _S3427.y;
    float2  _S3430 = v_uv0_0 + _S3392.differential_0 + - _S3398;
    float _S3431 = fx_28 * (_S3383.differential_0 + _S3387.differential_0 + _S3430.x);
    float2  _S3432 = make_float2 (_S3431, fy_28 * (_S3375.differential_0 + _S3379.differential_0 + _S3430.y)) + make_float2 ((*dist_coeffs_32)[int(8)] * _S3431, (*dist_coeffs_32)[int(9)] * _S3431);
    float2  _S3433 = _S3180 * _S3432;
    float _S3434 = (*dist_coeffs_32)[int(4)] * _S3432.y;
    float _S3435 = (*dist_coeffs_32)[int(5)] * _S3432.x;
    float _S3436 = _S3433.x + _S3433.y;
    float _S3437 = r2_56 * _S3436;
    float _S3438 = r2_56 * _S3437;
    float _S3439 = (*dist_coeffs_32)[int(7)] * _S3432.y + _S3434 + (*dist_coeffs_32)[int(6)] * _S3432.x + _S3435 + _S3184 * _S3436 + _S3183 * _S3437 + _S3182 * _S3438 + (*dist_coeffs_32)[int(3)] * (r2_56 * _S3438);
    float _S3440 = v_56 * _S3439;
    float _S3441 = u_56 * _S3439;
    float2  _S3442 = (_S3185 * _S3432 + make_float2 (_S3189 * (v_56 * _S3432.y) + _S3188 * _S3435 + 2.0f * (u_56 * _S3435) + _S3186 * (v_56 * _S3432.x) + _S3441 + _S3441, _S3191 * _S3434 + 2.0f * (v_56 * _S3434) + _S3190 * _S3432.y + _S3187 * _S3432.x + _S3440 + _S3440)) / _S3181;
    float2  _S3443 = _S3177 * - _S3442;
    float2  _S3444 = _S3179 * _S3442;
    float _S3445 = _S3443.x + _S3443.y;
    float3  _S3446 = _S3291.differential_0 + _S3370.differential_0 + make_float3 (_S3413.x, _S3413.y, _S3412.x + _S3412.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3447;
    (&_S3447)->primal_0 = R_25;
    (&_S3447)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3448;
    (&_S3448)->primal_0 = vert2_1;
    (&_S3448)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3447, &_S3448, _S3446);
    float3  _S3449 = _S3290.differential_0 + _S3371.differential_0 + make_float3 (_S3428.x, _S3428.y, _S3429);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3450;
    (&_S3450)->primal_0 = R_25;
    (&_S3450)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3451;
    (&_S3451)->primal_0 = vert1_1;
    (&_S3451)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3450, &_S3451, _S3449);
    float3  _S3452 = _S3292 + _S3293 + _S3372.differential_0 + make_float3 (_S3444.x, _S3444.y, _S3445);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3453;
    (&_S3453)->primal_0 = R_25;
    (&_S3453)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3454;
    (&_S3454)->primal_0 = vert0_1;
    (&_S3454)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3453, &_S3454, _S3452);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3455;
    (&_S3455)->primal_0 = _S3171;
    (&_S3455)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3456;
    (&_S3456)->primal_0 = _S3176;
    (&_S3456)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3455, &_S3456, _S3448.differential_0);
    float _S3457 = - _S3456.differential_0.y;
    float _S3458 = _S3175 * _S3456.differential_0.x;
    float _S3459 = - (_S3167 * _S3456.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3460;
    (&_S3460)->primal_0 = _S3171;
    (&_S3460)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3461;
    (&_S3461)->primal_0 = _S3174;
    (&_S3461)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3460, &_S3461, _S3451.differential_0);
    float _S3462 = _S3167 * _S3461.differential_0.x;
    float _S3463 = _S3173 * _S3461.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3464;
    (&_S3464)->primal_0 = _S3171;
    (&_S3464)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3465;
    (&_S3465)->primal_0 = _S3172;
    (&_S3465)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3464, &_S3465, _S3454.differential_0);
    Matrix<float, 3, 3>  _S3466 = transpose_0(_S3455.differential_0 + _S3460.differential_0 + _S3464.differential_0);
    float _S3467 = 2.0f * - _S3466.rows[int(2)].z;
    float _S3468 = 2.0f * _S3466.rows[int(2)].y;
    float _S3469 = 2.0f * _S3466.rows[int(2)].x;
    float _S3470 = 2.0f * _S3466.rows[int(1)].z;
    float _S3471 = 2.0f * - _S3466.rows[int(1)].y;
    float _S3472 = 2.0f * _S3466.rows[int(1)].x;
    float _S3473 = 2.0f * _S3466.rows[int(0)].z;
    float _S3474 = 2.0f * _S3466.rows[int(0)].y;
    float _S3475 = 2.0f * - _S3466.rows[int(0)].x;
    float _S3476 = - _S3472 + _S3474;
    float _S3477 = _S3469 + - _S3473;
    float _S3478 = - _S3468 + _S3470;
    float _S3479 = _S3468 + _S3470;
    float _S3480 = _S3469 + _S3473;
    float _S3481 = _S3472 + _S3474;
    float _S3482 = quat_25.w * (_S3471 + _S3475);
    float _S3483 = quat_25.z * (_S3467 + _S3475);
    float _S3484 = quat_25.y * (_S3467 + _S3471);
    float _S3485 = quat_25.x * _S3476 + quat_25.z * _S3479 + quat_25.y * _S3480 + _S3482 + _S3482;
    float _S3486 = quat_25.x * _S3477 + quat_25.w * _S3479 + quat_25.y * _S3481 + _S3483 + _S3483;
    float _S3487 = quat_25.x * _S3478 + quat_25.w * _S3480 + quat_25.z * _S3481 + _S3484 + _S3484;
    float _S3488 = quat_25.w * _S3476 + quat_25.z * _S3477 + quat_25.y * _S3478;
    float _S3489 = _S3459 + _S3462;
    float _S3490 = 0.5f * - _S3489;
    float _S3491 = _S3457 + _S3461.differential_0.y;
    DiffPair_float_0 _S3492;
    (&_S3492)->primal_0 = _S3168;
    (&_S3492)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3492, _S3491);
    float _S3493 = _S3490 + _S3492.differential_0;
    float _S3494 = _S3458 + _S3463 + _S3465.differential_0.x;
    DiffPair_float_0 _S3495;
    (&_S3495)->primal_0 = _S3166;
    (&_S3495)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3495, _S3494);
    float _S3496 = _S3490 + _S3495.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3497;
    (&_S3497)->primal_0 = R_25;
    (&_S3497)->differential_0 = _S3365;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3498;
    (&_S3498)->primal_0 = mean_24;
    (&_S3498)->differential_0 = _S3285;
    s_bwd_prop_mul_1(&_S3497, &_S3498, _S3287.differential_0);
    float3  _S3499 = _S3367.differential_0 + _S3446 + _S3449 + _S3452 + _S3287.differential_0;
    Matrix<float, 3, 3>  _S3500 = _S3368 + _S3447.differential_0 + _S3450.differential_0 + _S3453.differential_0 + _S3497.differential_0;
    FixedArray<float3 , 2>  _S3501;
    _S3501[int(0)] = _S3285;
    _S3501[int(1)] = _S3285;
    _S3501[int(1)] = _S3299;
    _S3501[int(0)] = _S3304;
    FixedArray<float3 , 16>  _S3502;
    _S3502[int(0)] = _S3285;
    _S3502[int(1)] = _S3285;
    _S3502[int(2)] = _S3285;
    _S3502[int(3)] = _S3285;
    _S3502[int(4)] = _S3285;
    _S3502[int(5)] = _S3285;
    _S3502[int(6)] = _S3285;
    _S3502[int(7)] = _S3285;
    _S3502[int(8)] = _S3285;
    _S3502[int(9)] = _S3285;
    _S3502[int(10)] = _S3285;
    _S3502[int(11)] = _S3285;
    _S3502[int(12)] = _S3285;
    _S3502[int(13)] = _S3285;
    _S3502[int(14)] = _S3285;
    _S3502[int(15)] = _S3285;
    _S3502[int(15)] = _S3306;
    _S3502[int(14)] = _S3308;
    _S3502[int(13)] = _S3310;
    _S3502[int(12)] = _S3312;
    _S3502[int(11)] = _S3314;
    _S3502[int(10)] = _S3316;
    _S3502[int(9)] = _S3318;
    _S3502[int(8)] = _S3326;
    _S3502[int(7)] = _S3328;
    _S3502[int(6)] = _S3330;
    _S3502[int(5)] = _S3332;
    _S3502[int(4)] = _S3334;
    _S3502[int(3)] = _S3345;
    _S3502[int(2)] = _S3347;
    _S3502[int(1)] = _S3349;
    _S3502[int(0)] = _S3362;
    float2  _S3503 = v_out_hardness_0 + make_float2 (0.0f, _S3390);
    float3  _S3504 = make_float3 (_S3496, _S3493, _S3489);
    float4  _S3505 = make_float4 (0.0f);
    *&((&_S3505)->w) = _S3485;
    *&((&_S3505)->z) = _S3486;
    *&((&_S3505)->y) = _S3487;
    *&((&_S3505)->x) = _S3488;
    *v_mean_7 = _S3363 + _S3448.differential_0 + _S3451.differential_0 + _S3454.differential_0 + _S3498.differential_0;
    *v_quat_6 = _S3505;
    *v_scale_6 = _S3504;
    *v_hardness_0 = _S3503;
    *v_sh_coeffs_5 = _S3502;
    *v_ch_coeffs_0 = _S3501;
    *v_R_7 = _S3500;
    *v_t_6 = _S3499;
    return;
}

inline __device__ void projection_opaque_triangle_fisheye_vjp(float3  mean_25, float4  quat_26, float3  scale_25, float2  hardness_5, FixedArray<float3 , 16>  * sh_coeffs_20, FixedArray<float3 , 2>  * ch_coeffs_5, Matrix<float, 3, 3>  R_26, float3  t_23, float fx_29, float fy_29, float cx_24, float cy_24, FixedArray<float, 10>  * dist_coeffs_33, uint image_width_20, uint image_height_20, float2  v_uv0_1, float2  v_uv1_1, float2  v_uv2_1, float3  v_depth_7, float2  v_out_hardness_1, FixedArray<float3 , 3>  * v_rgb_7, float3  v_normal_1, float3  * v_mean_8, float4  * v_quat_7, float3  * v_scale_7, float2  * v_hardness_1, FixedArray<float3 , 16>  * v_sh_coeffs_6, FixedArray<float3 , 2>  * v_ch_coeffs_1, Matrix<float, 3, 3>  * v_R_8, float3  * v_t_7)
{
    float3  mean_c_20 = s_primal_ctx_mul_0(R_26, mean_25) + t_23;
    float _S3506 = scale_25.x;
    float _S3507 = s_primal_ctx_exp_1(_S3506);
    float _S3508 = scale_25.y;
    float _S3509 = s_primal_ctx_exp_1(_S3508);
    float sz_6 = scale_25.z - 0.5f * (_S3506 + _S3508);
    float _S3510 = quat_26.y;
    float x2_26 = _S3510 * _S3510;
    float y2_26 = quat_26.z * quat_26.z;
    float z2_46 = quat_26.w * quat_26.w;
    float xy_26 = quat_26.y * quat_26.z;
    float xz_26 = quat_26.y * quat_26.w;
    float yz_26 = quat_26.z * quat_26.w;
    float wx_26 = quat_26.x * quat_26.y;
    float wy_26 = quat_26.x * quat_26.z;
    float wz_26 = quat_26.x * quat_26.w;
    Matrix<float, 3, 3>  _S3511 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_26 + z2_46), 2.0f * (xy_26 + wz_26), 2.0f * (xz_26 - wy_26), 2.0f * (xy_26 - wz_26), 1.0f - 2.0f * (x2_26 + z2_46), 2.0f * (yz_26 + wx_26), 2.0f * (xz_26 + wy_26), 2.0f * (yz_26 - wx_26), 1.0f - 2.0f * (x2_26 + y2_26)));
    float3  _S3512 = make_float3 (_S3507, 0.0f, 0.0f);
    float3  vert0_2 = s_primal_ctx_mul_0(_S3511, _S3512) + mean_25;
    float _S3513 = -0.5f + sz_6;
    float3  _S3514 = make_float3 (_S3507 * _S3513, _S3509, 0.0f);
    float3  vert1_2 = s_primal_ctx_mul_0(_S3511, _S3514) + mean_25;
    float _S3515 = -0.5f - sz_6;
    float3  _S3516 = make_float3 (_S3507 * _S3515, - _S3509, 0.0f);
    float3  vert2_2 = s_primal_ctx_mul_0(_S3511, _S3516) + mean_25;
    float3  vert0_c_5 = s_primal_ctx_mul_0(R_26, vert0_2) + t_23;
    float3  vert1_c_5 = s_primal_ctx_mul_0(R_26, vert1_2) + t_23;
    float3  vert2_c_5 = s_primal_ctx_mul_0(R_26, vert2_2) + t_23;
    float _S3517 = length_1(vert0_c_5);
    float _S3518 = length_1(vert1_c_5);
    float _S3519 = length_1(vert2_c_5);
    float2  _S3520 = float2 {vert0_c_5.x, vert0_c_5.y};
    float _S3521 = length_0(_S3520);
    float _S3522 = vert0_c_5.z;
    float _S3523 = s_primal_ctx_atan2_0(_S3521, _S3522);
    bool _S3524 = _S3523 < 0.00100000004749745f;
    float k_9;
    float _S3525;
    float _S3526;
    float _S3527;
    if(_S3524)
    {
        float _S3528 = 1.0f - _S3523 * _S3523 / 3.0f;
        float _S3529 = _S3522 * _S3522;
        k_9 = _S3528 / _S3522;
        _S3525 = 0.0f;
        _S3526 = _S3529;
        _S3527 = _S3528;
    }
    else
    {
        float _S3530 = _S3521 * _S3521;
        k_9 = _S3523 / _S3521;
        _S3525 = _S3530;
        _S3526 = 0.0f;
        _S3527 = 0.0f;
    }
    float2  _S3531 = make_float2 (k_9);
    float2  _S3532 = _S3520 * make_float2 (k_9);
    float u_59 = _S3532.x;
    float v_59 = _S3532.y;
    float r2_59 = u_59 * u_59 + v_59 * v_59;
    float _S3533 = (*dist_coeffs_33)[int(2)] + r2_59 * (*dist_coeffs_33)[int(3)];
    float _S3534 = (*dist_coeffs_33)[int(1)] + r2_59 * _S3533;
    float _S3535 = (*dist_coeffs_33)[int(0)] + r2_59 * _S3534;
    float radial_4 = 1.0f + r2_59 * _S3535;
    float2  _S3536 = make_float2 (radial_4);
    float _S3537 = 2.0f * (*dist_coeffs_33)[int(4)];
    float _S3538 = _S3537 * u_59;
    float _S3539 = 2.0f * u_59;
    float _S3540 = 2.0f * (*dist_coeffs_33)[int(5)];
    float _S3541 = _S3540 * u_59;
    float _S3542 = 2.0f * v_59;
    float2  _S3543 = _S3532 * make_float2 (radial_4) + make_float2 (_S3538 * v_59 + (*dist_coeffs_33)[int(5)] * (r2_59 + _S3539 * u_59) + (*dist_coeffs_33)[int(6)] * r2_59, _S3541 * v_59 + (*dist_coeffs_33)[int(4)] * (r2_59 + _S3542 * v_59) + (*dist_coeffs_33)[int(7)] * r2_59);
    float2  _S3544 = _S3543 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3543.x + (*dist_coeffs_33)[int(9)] * _S3543.y, 0.0f);
    float _S3545 = fx_29 * _S3544.x + cx_24;
    float _S3546 = fy_29 * _S3544.y + cy_24;
    float2  _S3547 = make_float2 (_S3545, _S3546);
    float2  _S3548 = float2 {vert1_c_5.x, vert1_c_5.y};
    float _S3549 = length_0(_S3548);
    float _S3550 = vert1_c_5.z;
    float _S3551 = s_primal_ctx_atan2_0(_S3549, _S3550);
    bool _S3552 = _S3551 < 0.00100000004749745f;
    float _S3553;
    float _S3554;
    float _S3555;
    if(_S3552)
    {
        float _S3556 = 1.0f - _S3551 * _S3551 / 3.0f;
        float _S3557 = _S3550 * _S3550;
        k_9 = _S3556 / _S3550;
        _S3553 = 0.0f;
        _S3554 = _S3557;
        _S3555 = _S3556;
    }
    else
    {
        float _S3558 = _S3549 * _S3549;
        k_9 = _S3551 / _S3549;
        _S3553 = _S3558;
        _S3554 = 0.0f;
        _S3555 = 0.0f;
    }
    float2  _S3559 = make_float2 (k_9);
    float2  _S3560 = _S3548 * make_float2 (k_9);
    float u_60 = _S3560.x;
    float v_60 = _S3560.y;
    float r2_60 = u_60 * u_60 + v_60 * v_60;
    float _S3561 = (*dist_coeffs_33)[int(2)] + r2_60 * (*dist_coeffs_33)[int(3)];
    float _S3562 = (*dist_coeffs_33)[int(1)] + r2_60 * _S3561;
    float _S3563 = (*dist_coeffs_33)[int(0)] + r2_60 * _S3562;
    float radial_5 = 1.0f + r2_60 * _S3563;
    float2  _S3564 = make_float2 (radial_5);
    float _S3565 = _S3537 * u_60;
    float _S3566 = 2.0f * u_60;
    float _S3567 = _S3540 * u_60;
    float _S3568 = 2.0f * v_60;
    float2  _S3569 = _S3560 * make_float2 (radial_5) + make_float2 (_S3565 * v_60 + (*dist_coeffs_33)[int(5)] * (r2_60 + _S3566 * u_60) + (*dist_coeffs_33)[int(6)] * r2_60, _S3567 * v_60 + (*dist_coeffs_33)[int(4)] * (r2_60 + _S3568 * v_60) + (*dist_coeffs_33)[int(7)] * r2_60);
    float2  _S3570 = _S3569 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3569.x + (*dist_coeffs_33)[int(9)] * _S3569.y, 0.0f);
    float _S3571 = fx_29 * _S3570.x + cx_24;
    float _S3572 = fy_29 * _S3570.y + cy_24;
    float2  _S3573 = make_float2 (_S3571, _S3572);
    float2  _S3574 = float2 {vert2_c_5.x, vert2_c_5.y};
    float _S3575 = length_0(_S3574);
    float _S3576 = vert2_c_5.z;
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
    float u_61 = _S3586.x;
    float v_61 = _S3586.y;
    float r2_61 = u_61 * u_61 + v_61 * v_61;
    float _S3587 = (*dist_coeffs_33)[int(2)] + r2_61 * (*dist_coeffs_33)[int(3)];
    float _S3588 = (*dist_coeffs_33)[int(1)] + r2_61 * _S3587;
    float _S3589 = (*dist_coeffs_33)[int(0)] + r2_61 * _S3588;
    float radial_6 = 1.0f + r2_61 * _S3589;
    float2  _S3590 = make_float2 (radial_6);
    float _S3591 = _S3537 * u_61;
    float _S3592 = 2.0f * u_61;
    float _S3593 = _S3540 * u_61;
    float _S3594 = 2.0f * v_61;
    float2  _S3595 = _S3586 * make_float2 (radial_6) + make_float2 (_S3591 * v_61 + (*dist_coeffs_33)[int(5)] * (r2_61 + _S3592 * u_61) + (*dist_coeffs_33)[int(6)] * r2_61, _S3593 * v_61 + (*dist_coeffs_33)[int(4)] * (r2_61 + _S3594 * v_61) + (*dist_coeffs_33)[int(7)] * r2_61);
    float2  _S3596 = _S3595 + make_float2 ((*dist_coeffs_33)[int(8)] * _S3595.x + (*dist_coeffs_33)[int(9)] * _S3595.y, 0.0f);
    float _S3597 = fx_29 * _S3596.x + cx_24;
    float _S3598 = fy_29 * _S3596.y + cy_24;
    float2  _S3599 = make_float2 (_S3597, _S3598);
    float2  e0_5 = _S3573 - _S3547;
    float2  e1_5 = _S3599 - _S3573;
    float2  e2_1 = _S3547 - _S3599;
    float _S3600 = e0_5.x;
    float _S3601 = e1_5.y;
    float _S3602 = e0_5.y;
    float _S3603 = e1_5.x;
    float _S3604 = _S3600 * _S3601 - _S3602 * _S3603;
    float _S3605 = 1.0f - hardness_5.y;
    float _S3606 = -1.0f / _S3605;
    float _S3607 = _S3605 * _S3605;
    float _S3608 = s_primal_ctx_max_0(_S3545, _S3571);
    float _S3609 = s_primal_ctx_min_0(_S3545, _S3571);
    float _S3610 = s_primal_ctx_max_0(_S3546, _S3572);
    float _S3611 = s_primal_ctx_min_0(_S3546, _S3572);
    float3  _S3612 = make_float3 (_S3517, _S3518, _S3519) + make_float3 (9.999999960041972e-13f);
    Matrix<float, 3, 3>  _S3613 = transpose_0(R_26);
    float3  _S3614 = mean_25 - - s_primal_ctx_mul_0(_S3613, t_23);
    float _S3615 = _S3614.x;
    float _S3616 = _S3614.y;
    float _S3617 = _S3614.z;
    float _S3618 = _S3615 * _S3615 + _S3616 * _S3616 + _S3617 * _S3617;
    float _S3619 = s_primal_ctx_sqrt_0(_S3618);
    float x_52 = _S3615 / _S3619;
    float3  _S3620 = make_float3 (x_52);
    float _S3621 = _S3619 * _S3619;
    float y_23 = _S3616 / _S3619;
    float z_20 = _S3617 / _S3619;
    float3  _S3622 = make_float3 (z_20);
    float _S3623 = - y_23;
    float3  _S3624 = make_float3 (_S3623);
    float z2_47 = z_20 * z_20;
    float fTmp0B_20 = -1.09254848957061768f * z_20;
    float fC1_20 = x_52 * x_52 - y_23 * y_23;
    float _S3625 = 2.0f * x_52;
    float fS1_20 = _S3625 * y_23;
    float pSH6_6 = 0.94617468118667603f * z2_47 - 0.31539157032966614f;
    float3  _S3626 = make_float3 (pSH6_6);
    float pSH7_6 = fTmp0B_20 * x_52;
    float3  _S3627 = make_float3 (pSH7_6);
    float pSH5_6 = fTmp0B_20 * y_23;
    float3  _S3628 = make_float3 (pSH5_6);
    float pSH8_6 = 0.54627424478530884f * fC1_20;
    float3  _S3629 = make_float3 (pSH8_6);
    float pSH4_6 = 0.54627424478530884f * fS1_20;
    float3  _S3630 = make_float3 (pSH4_6);
    float fTmp0C_20 = -2.28522896766662598f * z2_47 + 0.4570457935333252f;
    float fTmp1B_20 = 1.44530570507049561f * z_20;
    float _S3631 = 1.86588168144226074f * z2_47 - 1.11952900886535645f;
    float pSH12_6 = z_20 * _S3631;
    float3  _S3632 = make_float3 (pSH12_6);
    float pSH13_6 = fTmp0C_20 * x_52;
    float3  _S3633 = make_float3 (pSH13_6);
    float pSH11_6 = fTmp0C_20 * y_23;
    float3  _S3634 = make_float3 (pSH11_6);
    float pSH14_6 = fTmp1B_20 * fC1_20;
    float3  _S3635 = make_float3 (pSH14_6);
    float pSH10_6 = fTmp1B_20 * fS1_20;
    float3  _S3636 = make_float3 (pSH10_6);
    float pSH15_6 = -0.59004360437393188f * (x_52 * fC1_20 - y_23 * fS1_20);
    float3  _S3637 = make_float3 (pSH15_6);
    float pSH9_6 = -0.59004360437393188f * (x_52 * fS1_20 + y_23 * fC1_20);
    float3  _S3638 = make_float3 (pSH9_6);
    float3  color_5 = make_float3 (0.282094806432724f) * (*sh_coeffs_20)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S3623) * (*sh_coeffs_20)[int(1)] + make_float3 (z_20) * (*sh_coeffs_20)[int(2)] - make_float3 (x_52) * (*sh_coeffs_20)[int(3)]) + (make_float3 (pSH4_6) * (*sh_coeffs_20)[int(4)] + make_float3 (pSH5_6) * (*sh_coeffs_20)[int(5)] + make_float3 (pSH6_6) * (*sh_coeffs_20)[int(6)] + make_float3 (pSH7_6) * (*sh_coeffs_20)[int(7)] + make_float3 (pSH8_6) * (*sh_coeffs_20)[int(8)]) + (make_float3 (pSH9_6) * (*sh_coeffs_20)[int(9)] + make_float3 (pSH10_6) * (*sh_coeffs_20)[int(10)] + make_float3 (pSH11_6) * (*sh_coeffs_20)[int(11)] + make_float3 (pSH12_6) * (*sh_coeffs_20)[int(12)] + make_float3 (pSH13_6) * (*sh_coeffs_20)[int(13)] + make_float3 (pSH14_6) * (*sh_coeffs_20)[int(14)] + make_float3 (pSH15_6) * (*sh_coeffs_20)[int(15)]);
    float3  _S3639 = color_5 + (*ch_coeffs_5)[int(0)] + make_float3 (0.5f);
    float3  _S3640 = make_float3 (0.0f);
    float3  _S3641 = color_5 - (*ch_coeffs_5)[int(0)] * make_float3 (0.5f);
    float _S3642 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S3643 = make_float3 (_S3642);
    float3  _S3644 = (*ch_coeffs_5)[int(1)] * make_float3 (_S3642);
    float3  _S3645 = _S3641 + _S3644 + make_float3 (0.5f);
    float3  _S3646 = _S3641 - _S3644 + make_float3 (0.5f);
    float3  _S3647 = vert1_c_5 - vert0_c_5;
    float3  _S3648 = vert2_c_5 - vert0_c_5;
    float3  _S3649 = s_primal_ctx_cross_0(_S3647, _S3648);
    float3  _S3650 = normalize_0(_S3649);
    float3  _S3651 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S3650, mean_c_20)))))) * v_normal_1;
    float3  _S3652 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3653;
    (&_S3653)->primal_0 = _S3650;
    (&_S3653)->differential_0 = _S3652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3654;
    (&_S3654)->primal_0 = mean_c_20;
    (&_S3654)->differential_0 = _S3652;
    s_bwd_prop_dot_0(&_S3653, &_S3654, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3655 = _S3654;
    float3  _S3656 = _S3651 + _S3653.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3657;
    (&_S3657)->primal_0 = _S3649;
    (&_S3657)->differential_0 = _S3652;
    s_bwd_normalize_impl_0(&_S3657, _S3656);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3658;
    (&_S3658)->primal_0 = _S3647;
    (&_S3658)->differential_0 = _S3652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3659;
    (&_S3659)->primal_0 = _S3648;
    (&_S3659)->differential_0 = _S3652;
    s_bwd_prop_cross_0(&_S3658, &_S3659, _S3657.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3660 = _S3658;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3661 = _S3659;
    float3  _S3662 = - _S3659.differential_0;
    float3  _S3663 = - _S3658.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3664;
    (&_S3664)->primal_0 = _S3646;
    (&_S3664)->differential_0 = _S3652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3665;
    (&_S3665)->primal_0 = _S3640;
    (&_S3665)->differential_0 = _S3652;
    s_bwd_prop_max_0(&_S3664, &_S3665, (*v_rgb_7)[int(2)]);
    float3  _S3666 = - _S3664.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3667;
    (&_S3667)->primal_0 = _S3645;
    (&_S3667)->differential_0 = _S3652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3668;
    (&_S3668)->primal_0 = _S3640;
    (&_S3668)->differential_0 = _S3652;
    s_bwd_prop_max_0(&_S3667, &_S3668, (*v_rgb_7)[int(1)]);
    float3  _S3669 = _S3643 * (_S3666 + _S3667.differential_0);
    float3  _S3670 = _S3664.differential_0 + _S3667.differential_0;
    float3  _S3671 = make_float3 (0.5f) * - _S3670;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3672;
    (&_S3672)->primal_0 = _S3639;
    (&_S3672)->differential_0 = _S3652;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3673;
    (&_S3673)->primal_0 = _S3640;
    (&_S3673)->differential_0 = _S3652;
    s_bwd_prop_max_0(&_S3672, &_S3673, (*v_rgb_7)[int(0)]);
    float3  _S3674 = _S3671 + _S3672.differential_0;
    float3  _S3675 = _S3670 + _S3672.differential_0;
    float3  _S3676 = _S3637 * _S3675;
    float3  _S3677 = (*sh_coeffs_20)[int(15)] * _S3675;
    float3  _S3678 = _S3635 * _S3675;
    float3  _S3679 = (*sh_coeffs_20)[int(14)] * _S3675;
    float3  _S3680 = _S3633 * _S3675;
    float3  _S3681 = (*sh_coeffs_20)[int(13)] * _S3675;
    float3  _S3682 = _S3632 * _S3675;
    float3  _S3683 = (*sh_coeffs_20)[int(12)] * _S3675;
    float3  _S3684 = _S3634 * _S3675;
    float3  _S3685 = (*sh_coeffs_20)[int(11)] * _S3675;
    float3  _S3686 = _S3636 * _S3675;
    float3  _S3687 = (*sh_coeffs_20)[int(10)] * _S3675;
    float3  _S3688 = _S3638 * _S3675;
    float3  _S3689 = (*sh_coeffs_20)[int(9)] * _S3675;
    float s_diff_fS2_T_6 = -0.59004360437393188f * (_S3689.x + _S3689.y + _S3689.z);
    float s_diff_fC2_T_6 = -0.59004360437393188f * (_S3677.x + _S3677.y + _S3677.z);
    float _S3690 = _S3687.x + _S3687.y + _S3687.z;
    float _S3691 = _S3679.x + _S3679.y + _S3679.z;
    float _S3692 = _S3685.x + _S3685.y + _S3685.z;
    float _S3693 = _S3681.x + _S3681.y + _S3681.z;
    float _S3694 = _S3683.x + _S3683.y + _S3683.z;
    float _S3695 = - s_diff_fC2_T_6;
    float3  _S3696 = _S3629 * _S3675;
    float3  _S3697 = (*sh_coeffs_20)[int(8)] * _S3675;
    float3  _S3698 = _S3627 * _S3675;
    float3  _S3699 = (*sh_coeffs_20)[int(7)] * _S3675;
    float3  _S3700 = _S3626 * _S3675;
    float3  _S3701 = (*sh_coeffs_20)[int(6)] * _S3675;
    float3  _S3702 = _S3628 * _S3675;
    float3  _S3703 = (*sh_coeffs_20)[int(5)] * _S3675;
    float3  _S3704 = _S3630 * _S3675;
    float3  _S3705 = (*sh_coeffs_20)[int(4)] * _S3675;
    float _S3706 = _S3703.x + _S3703.y + _S3703.z;
    float _S3707 = _S3699.x + _S3699.y + _S3699.z;
    float _S3708 = fTmp1B_20 * _S3690 + x_52 * s_diff_fS2_T_6 + y_23 * _S3695 + 0.54627424478530884f * (_S3705.x + _S3705.y + _S3705.z);
    float _S3709 = fTmp1B_20 * _S3691 + y_23 * s_diff_fS2_T_6 + x_52 * s_diff_fC2_T_6 + 0.54627424478530884f * (_S3697.x + _S3697.y + _S3697.z);
    float _S3710 = y_23 * - _S3709;
    float _S3711 = x_52 * _S3709;
    float _S3712 = z_20 * (1.86588168144226074f * (z_20 * _S3694) + -2.28522896766662598f * (y_23 * _S3692 + x_52 * _S3693) + 0.94617468118667603f * (_S3701.x + _S3701.y + _S3701.z));
    float3  _S3713 = make_float3 (0.48860251903533936f) * _S3675;
    float3  _S3714 = - _S3713;
    float3  _S3715 = _S3620 * _S3714;
    float3  _S3716 = (*sh_coeffs_20)[int(3)] * _S3714;
    float3  _S3717 = _S3622 * _S3713;
    float3  _S3718 = (*sh_coeffs_20)[int(2)] * _S3713;
    float3  _S3719 = _S3624 * _S3713;
    float3  _S3720 = (*sh_coeffs_20)[int(1)] * _S3713;
    float _S3721 = (_S3631 * _S3694 + 1.44530570507049561f * (fS1_20 * _S3690 + fC1_20 * _S3691) + -1.09254848957061768f * (y_23 * _S3706 + x_52 * _S3707) + _S3712 + _S3712 + _S3718.x + _S3718.y + _S3718.z) / _S3621;
    float _S3722 = _S3619 * _S3721;
    float _S3723 = (fTmp0C_20 * _S3692 + fC1_20 * s_diff_fS2_T_6 + fS1_20 * _S3695 + fTmp0B_20 * _S3706 + _S3625 * _S3708 + _S3710 + _S3710 + - (_S3720.x + _S3720.y + _S3720.z)) / _S3621;
    float _S3724 = _S3619 * _S3723;
    float _S3725 = (fTmp0C_20 * _S3693 + fS1_20 * s_diff_fS2_T_6 + fC1_20 * s_diff_fC2_T_6 + fTmp0B_20 * _S3707 + 2.0f * (y_23 * _S3708) + _S3711 + _S3711 + _S3716.x + _S3716.y + _S3716.z) / _S3621;
    float _S3726 = _S3619 * _S3725;
    float _S3727 = _S3617 * - _S3721 + _S3616 * - _S3723 + _S3615 * - _S3725;
    DiffPair_float_0 _S3728;
    (&_S3728)->primal_0 = _S3618;
    (&_S3728)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S3728, _S3727);
    float _S3729 = _S3617 * _S3728.differential_0;
    float _S3730 = _S3616 * _S3728.differential_0;
    float _S3731 = _S3615 * _S3728.differential_0;
    float3  _S3732 = make_float3 (0.282094806432724f) * _S3675;
    float3  _S3733 = make_float3 (_S3726 + _S3731 + _S3731, _S3724 + _S3730 + _S3730, _S3722 + _S3729 + _S3729);
    float3  _S3734 = - - _S3733;
    Matrix<float, 3, 3>  _S3735 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3736;
    (&_S3736)->primal_0 = _S3613;
    (&_S3736)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3737;
    (&_S3737)->primal_0 = t_23;
    (&_S3737)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3736, &_S3737, _S3734);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3738 = _S3737;
    Matrix<float, 3, 3>  _S3739 = transpose_0(_S3736.differential_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3740;
    (&_S3740)->primal_0 = _S3612;
    (&_S3740)->differential_0 = _S3652;
    s_bwd_prop_log_1(&_S3740, v_depth_7);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3741 = _S3740;
    DiffPair_float_0 _S3742;
    (&_S3742)->primal_0 = _S3611;
    (&_S3742)->differential_0 = 0.0f;
    DiffPair_float_0 _S3743;
    (&_S3743)->primal_0 = _S3598;
    (&_S3743)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3742, &_S3743, 0.0f);
    DiffPair_float_0 _S3744;
    (&_S3744)->primal_0 = _S3546;
    (&_S3744)->differential_0 = 0.0f;
    DiffPair_float_0 _S3745;
    (&_S3745)->primal_0 = _S3572;
    (&_S3745)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3744, &_S3745, _S3742.differential_0);
    DiffPair_float_0 _S3746;
    (&_S3746)->primal_0 = _S3610;
    (&_S3746)->differential_0 = 0.0f;
    DiffPair_float_0 _S3747;
    (&_S3747)->primal_0 = _S3598;
    (&_S3747)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3746, &_S3747, 0.0f);
    DiffPair_float_0 _S3748;
    (&_S3748)->primal_0 = _S3546;
    (&_S3748)->differential_0 = 0.0f;
    DiffPair_float_0 _S3749;
    (&_S3749)->primal_0 = _S3572;
    (&_S3749)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3748, &_S3749, _S3746.differential_0);
    DiffPair_float_0 _S3750;
    (&_S3750)->primal_0 = _S3609;
    (&_S3750)->differential_0 = 0.0f;
    DiffPair_float_0 _S3751;
    (&_S3751)->primal_0 = _S3597;
    (&_S3751)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3750, &_S3751, 0.0f);
    DiffPair_float_0 _S3752;
    (&_S3752)->primal_0 = _S3545;
    (&_S3752)->differential_0 = 0.0f;
    DiffPair_float_0 _S3753;
    (&_S3753)->primal_0 = _S3571;
    (&_S3753)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S3752, &_S3753, _S3750.differential_0);
    DiffPair_float_0 _S3754;
    (&_S3754)->primal_0 = _S3608;
    (&_S3754)->differential_0 = 0.0f;
    DiffPair_float_0 _S3755;
    (&_S3755)->primal_0 = _S3597;
    (&_S3755)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3754, &_S3755, 0.0f);
    DiffPair_float_0 _S3756;
    (&_S3756)->primal_0 = _S3545;
    (&_S3756)->differential_0 = 0.0f;
    DiffPair_float_0 _S3757;
    (&_S3757)->primal_0 = _S3571;
    (&_S3757)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S3756, &_S3757, _S3754.differential_0);
    DiffPair_float_0 _S3758;
    (&_S3758)->primal_0 = _S3606;
    (&_S3758)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S3758, 0.0f);
    float _S3759 = - (-1.0f * - (_S3758.differential_0 / _S3607));
    float2  _S3760 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3761;
    (&_S3761)->primal_0 = e2_1;
    (&_S3761)->differential_0 = _S3760;
    s_bwd_length_impl_0(&_S3761, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3762;
    (&_S3762)->primal_0 = e1_5;
    (&_S3762)->differential_0 = _S3760;
    s_bwd_length_impl_0(&_S3762, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3763;
    (&_S3763)->primal_0 = e0_5;
    (&_S3763)->differential_0 = _S3760;
    s_bwd_length_impl_0(&_S3763, -0.0f);
    DiffPair_float_0 _S3764;
    (&_S3764)->primal_0 = _S3604;
    (&_S3764)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S3764, 0.0f);
    float _S3765 = - _S3764.differential_0;
    float2  _S3766 = _S3762.differential_0 + make_float2 (_S3602 * _S3765, _S3600 * _S3764.differential_0);
    float2  _S3767 = _S3763.differential_0 + make_float2 (_S3601 * _S3764.differential_0, _S3603 * _S3765);
    float2  _S3768 = v_uv2_1 + - _S3761.differential_0 + _S3766;
    float _S3769 = fx_29 * (_S3751.differential_0 + _S3755.differential_0 + _S3768.x);
    float2  _S3770 = make_float2 (_S3769, fy_29 * (_S3743.differential_0 + _S3747.differential_0 + _S3768.y)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3769, (*dist_coeffs_33)[int(9)] * _S3769);
    float2  _S3771 = _S3586 * _S3770;
    float2  _S3772 = _S3590 * _S3770;
    float _S3773 = (*dist_coeffs_33)[int(4)] * _S3770.y;
    float _S3774 = (*dist_coeffs_33)[int(5)] * _S3770.x;
    float _S3775 = _S3771.x + _S3771.y;
    float _S3776 = r2_61 * _S3775;
    float _S3777 = r2_61 * _S3776;
    float _S3778 = (*dist_coeffs_33)[int(7)] * _S3770.y + _S3773 + (*dist_coeffs_33)[int(6)] * _S3770.x + _S3774 + _S3589 * _S3775 + _S3588 * _S3776 + _S3587 * _S3777 + (*dist_coeffs_33)[int(3)] * (r2_61 * _S3777);
    float _S3779 = v_61 * _S3778;
    float _S3780 = u_61 * _S3778;
    float _S3781 = _S3594 * _S3773 + 2.0f * (v_61 * _S3773) + _S3593 * _S3770.y + _S3591 * _S3770.x + _S3779 + _S3779;
    float _S3782 = _S3540 * (v_61 * _S3770.y) + _S3592 * _S3774 + 2.0f * (u_61 * _S3774) + _S3537 * (v_61 * _S3770.x) + _S3780 + _S3780;
    float2  _S3783 = v_uv0_1 + _S3761.differential_0 + - _S3767;
    float2  _S3784 = v_out_hardness_1 + make_float2 (0.0f, _S3759);
    float _S3785 = _S3753.differential_0 + _S3757.differential_0;
    float2  _S3786 = v_uv1_1 + - _S3766 + _S3767;
    float3  _S3787 = _S3662 + _S3663;
    FixedArray<float3 , 2>  _S3788;
    _S3788[int(0)] = _S3652;
    _S3788[int(1)] = _S3652;
    _S3788[int(1)] = _S3669;
    _S3788[int(0)] = _S3674;
    float3  _S3789 = _S3788[int(0)];
    float3  _S3790 = _S3788[int(1)];
    FixedArray<float3 , 16>  _S3791;
    _S3791[int(0)] = _S3652;
    _S3791[int(1)] = _S3652;
    _S3791[int(2)] = _S3652;
    _S3791[int(3)] = _S3652;
    _S3791[int(4)] = _S3652;
    _S3791[int(5)] = _S3652;
    _S3791[int(6)] = _S3652;
    _S3791[int(7)] = _S3652;
    _S3791[int(8)] = _S3652;
    _S3791[int(9)] = _S3652;
    _S3791[int(10)] = _S3652;
    _S3791[int(11)] = _S3652;
    _S3791[int(12)] = _S3652;
    _S3791[int(13)] = _S3652;
    _S3791[int(14)] = _S3652;
    _S3791[int(15)] = _S3652;
    _S3791[int(7)] = _S3698;
    _S3791[int(0)] = _S3732;
    _S3791[int(1)] = _S3719;
    _S3791[int(2)] = _S3717;
    _S3791[int(3)] = _S3715;
    _S3791[int(4)] = _S3704;
    _S3791[int(5)] = _S3702;
    _S3791[int(6)] = _S3700;
    _S3791[int(15)] = _S3676;
    _S3791[int(8)] = _S3696;
    _S3791[int(9)] = _S3688;
    _S3791[int(10)] = _S3686;
    _S3791[int(11)] = _S3684;
    _S3791[int(12)] = _S3682;
    _S3791[int(13)] = _S3680;
    _S3791[int(14)] = _S3678;
    float3  _S3792 = _S3791[int(0)];
    float3  _S3793 = _S3791[int(1)];
    float3  _S3794 = _S3791[int(2)];
    float3  _S3795 = _S3791[int(3)];
    float3  _S3796 = _S3791[int(4)];
    float3  _S3797 = _S3791[int(5)];
    float3  _S3798 = _S3791[int(6)];
    float3  _S3799 = _S3791[int(7)];
    float3  _S3800 = _S3791[int(8)];
    float3  _S3801 = _S3791[int(9)];
    float3  _S3802 = _S3791[int(10)];
    float3  _S3803 = _S3791[int(11)];
    float3  _S3804 = _S3791[int(12)];
    float3  _S3805 = _S3791[int(13)];
    float3  _S3806 = _S3791[int(14)];
    float3  _S3807 = _S3791[int(15)];
    float _S3808 = _S3752.differential_0 + _S3756.differential_0;
    float _S3809 = _S3744.differential_0 + _S3748.differential_0;
    float _S3810 = _S3745.differential_0 + _S3749.differential_0;
    float2  _S3811 = _S3772 + make_float2 (_S3782, _S3781);
    float2  _S3812 = _S3574 * _S3811;
    float2  _S3813 = _S3585 * _S3811;
    float _S3814 = _S3812.x + _S3812.y;
    if(_S3578)
    {
        float _S3815 = _S3814 / _S3580;
        float _S3816 = _S3581 * - _S3815;
        float _S3817 = _S3577 * (0.3333333432674408f * - (_S3576 * _S3815));
        k_9 = _S3817 + _S3817;
        _S3579 = _S3816;
        _S3580 = 0.0f;
    }
    else
    {
        float _S3818 = _S3814 / _S3579;
        float _S3819 = _S3577 * - _S3818;
        k_9 = _S3575 * _S3818;
        _S3579 = 0.0f;
        _S3580 = _S3819;
    }
    DiffPair_float_0 _S3820;
    (&_S3820)->primal_0 = _S3575;
    (&_S3820)->differential_0 = 0.0f;
    DiffPair_float_0 _S3821;
    (&_S3821)->primal_0 = _S3576;
    (&_S3821)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3820, &_S3821, k_9);
    float _S3822 = _S3821.differential_0 + _S3579;
    float _S3823 = _S3820.differential_0 + _S3580;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3824;
    (&_S3824)->primal_0 = _S3574;
    (&_S3824)->differential_0 = _S3760;
    s_bwd_length_impl_0(&_S3824, _S3823);
    float2  _S3825 = _S3824.differential_0 + _S3813;
    float _S3826 = fx_29 * (_S3786.x + _S3785);
    float2  _S3827 = make_float2 (_S3826, fy_29 * (_S3786.y + _S3810)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3826, (*dist_coeffs_33)[int(9)] * _S3826);
    float2  _S3828 = _S3560 * _S3827;
    float _S3829 = (*dist_coeffs_33)[int(4)] * _S3827.y;
    float _S3830 = (*dist_coeffs_33)[int(5)] * _S3827.x;
    float _S3831 = _S3828.x + _S3828.y;
    float _S3832 = r2_60 * _S3831;
    float _S3833 = r2_60 * _S3832;
    float _S3834 = (*dist_coeffs_33)[int(7)] * _S3827.y + _S3829 + (*dist_coeffs_33)[int(6)] * _S3827.x + _S3830 + _S3563 * _S3831 + _S3562 * _S3832 + _S3561 * _S3833 + (*dist_coeffs_33)[int(3)] * (r2_60 * _S3833);
    float _S3835 = v_60 * _S3834;
    float _S3836 = u_60 * _S3834;
    float3  _S3837 = _S3661.differential_0 + make_float3 (_S3825.x, _S3825.y, _S3822);
    float2  _S3838 = _S3564 * _S3827 + make_float2 (_S3540 * (v_60 * _S3827.y) + _S3566 * _S3830 + 2.0f * (u_60 * _S3830) + _S3537 * (v_60 * _S3827.x) + _S3836 + _S3836, _S3568 * _S3829 + 2.0f * (v_60 * _S3829) + _S3567 * _S3827.y + _S3565 * _S3827.x + _S3835 + _S3835);
    float2  _S3839 = _S3548 * _S3838;
    float2  _S3840 = _S3559 * _S3838;
    float _S3841 = _S3839.x + _S3839.y;
    if(_S3552)
    {
        float _S3842 = _S3841 / _S3554;
        float _S3843 = _S3555 * - _S3842;
        float _S3844 = _S3551 * (0.3333333432674408f * - (_S3550 * _S3842));
        k_9 = _S3844 + _S3844;
        _S3553 = _S3843;
        _S3554 = 0.0f;
    }
    else
    {
        float _S3845 = _S3841 / _S3553;
        float _S3846 = _S3551 * - _S3845;
        k_9 = _S3549 * _S3845;
        _S3553 = 0.0f;
        _S3554 = _S3846;
    }
    DiffPair_float_0 _S3847;
    (&_S3847)->primal_0 = _S3549;
    (&_S3847)->differential_0 = 0.0f;
    DiffPair_float_0 _S3848;
    (&_S3848)->primal_0 = _S3550;
    (&_S3848)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3847, &_S3848, k_9);
    float _S3849 = _S3848.differential_0 + _S3553;
    float _S3850 = _S3847.differential_0 + _S3554;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3851;
    (&_S3851)->primal_0 = _S3548;
    (&_S3851)->differential_0 = _S3760;
    s_bwd_length_impl_0(&_S3851, _S3850);
    float2  _S3852 = _S3851.differential_0 + _S3840;
    float _S3853 = fx_29 * (_S3783.x + _S3808);
    float2  _S3854 = make_float2 (_S3853, fy_29 * (_S3783.y + _S3809)) + make_float2 ((*dist_coeffs_33)[int(8)] * _S3853, (*dist_coeffs_33)[int(9)] * _S3853);
    float2  _S3855 = _S3532 * _S3854;
    float _S3856 = (*dist_coeffs_33)[int(4)] * _S3854.y;
    float _S3857 = (*dist_coeffs_33)[int(5)] * _S3854.x;
    float _S3858 = _S3855.x + _S3855.y;
    float _S3859 = r2_59 * _S3858;
    float _S3860 = r2_59 * _S3859;
    float _S3861 = (*dist_coeffs_33)[int(7)] * _S3854.y + _S3856 + (*dist_coeffs_33)[int(6)] * _S3854.x + _S3857 + _S3535 * _S3858 + _S3534 * _S3859 + _S3533 * _S3860 + (*dist_coeffs_33)[int(3)] * (r2_59 * _S3860);
    float _S3862 = v_59 * _S3861;
    float _S3863 = u_59 * _S3861;
    float3  _S3864 = _S3660.differential_0 + make_float3 (_S3852.x, _S3852.y, _S3849);
    float2  _S3865 = _S3536 * _S3854 + make_float2 (_S3540 * (v_59 * _S3854.y) + _S3539 * _S3857 + 2.0f * (u_59 * _S3857) + _S3537 * (v_59 * _S3854.x) + _S3863 + _S3863, _S3542 * _S3856 + 2.0f * (v_59 * _S3856) + _S3541 * _S3854.y + _S3538 * _S3854.x + _S3862 + _S3862);
    float2  _S3866 = _S3520 * _S3865;
    float2  _S3867 = _S3531 * _S3865;
    float _S3868 = _S3866.x + _S3866.y;
    if(_S3524)
    {
        float _S3869 = _S3868 / _S3526;
        float _S3870 = _S3527 * - _S3869;
        float _S3871 = _S3523 * (0.3333333432674408f * - (_S3522 * _S3869));
        k_9 = _S3871 + _S3871;
        _S3525 = _S3870;
        _S3526 = 0.0f;
    }
    else
    {
        float _S3872 = _S3868 / _S3525;
        float _S3873 = _S3523 * - _S3872;
        k_9 = _S3521 * _S3872;
        _S3525 = 0.0f;
        _S3526 = _S3873;
    }
    DiffPair_float_0 _S3874;
    (&_S3874)->primal_0 = _S3521;
    (&_S3874)->differential_0 = 0.0f;
    DiffPair_float_0 _S3875;
    (&_S3875)->primal_0 = _S3522;
    (&_S3875)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S3874, &_S3875, k_9);
    float _S3876 = _S3875.differential_0 + _S3525;
    float _S3877 = _S3874.differential_0 + _S3526;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3878;
    (&_S3878)->primal_0 = _S3520;
    (&_S3878)->differential_0 = _S3760;
    s_bwd_length_impl_0(&_S3878, _S3877);
    float2  _S3879 = _S3878.differential_0 + _S3867;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3880;
    (&_S3880)->primal_0 = vert2_c_5;
    (&_S3880)->differential_0 = _S3652;
    s_bwd_length_impl_1(&_S3880, _S3741.differential_0.z);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3881;
    (&_S3881)->primal_0 = vert1_c_5;
    (&_S3881)->differential_0 = _S3652;
    s_bwd_length_impl_1(&_S3881, _S3741.differential_0.y);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3882;
    (&_S3882)->primal_0 = vert0_c_5;
    (&_S3882)->differential_0 = _S3652;
    s_bwd_length_impl_1(&_S3882, _S3741.differential_0.x);
    float3  _S3883 = _S3880.differential_0 + _S3837;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3884;
    (&_S3884)->primal_0 = R_26;
    (&_S3884)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3885;
    (&_S3885)->primal_0 = vert2_2;
    (&_S3885)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3884, &_S3885, _S3883);
    float3  _S3886 = _S3881.differential_0 + _S3864;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3887;
    (&_S3887)->primal_0 = R_26;
    (&_S3887)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3888;
    (&_S3888)->primal_0 = vert1_2;
    (&_S3888)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3887, &_S3888, _S3886);
    float3  _S3889 = _S3882.differential_0 + _S3787 + make_float3 (_S3879.x, _S3879.y, _S3876);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3890;
    (&_S3890)->primal_0 = R_26;
    (&_S3890)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3891;
    (&_S3891)->primal_0 = vert0_2;
    (&_S3891)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3890, &_S3891, _S3889);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3892;
    (&_S3892)->primal_0 = _S3511;
    (&_S3892)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3893;
    (&_S3893)->primal_0 = _S3516;
    (&_S3893)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3892, &_S3893, _S3885.differential_0);
    float _S3894 = - _S3893.differential_0.y;
    float _S3895 = _S3515 * _S3893.differential_0.x;
    float _S3896 = - (_S3507 * _S3893.differential_0.x);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3897;
    (&_S3897)->primal_0 = _S3511;
    (&_S3897)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3898;
    (&_S3898)->primal_0 = _S3514;
    (&_S3898)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3897, &_S3898, _S3888.differential_0);
    float _S3899 = _S3507 * _S3898.differential_0.x;
    float _S3900 = _S3513 * _S3898.differential_0.x;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3901;
    (&_S3901)->primal_0 = _S3511;
    (&_S3901)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3902;
    (&_S3902)->primal_0 = _S3512;
    (&_S3902)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3901, &_S3902, _S3891.differential_0);
    Matrix<float, 3, 3>  _S3903 = transpose_0(_S3892.differential_0 + _S3897.differential_0 + _S3901.differential_0);
    float _S3904 = 2.0f * - _S3903.rows[int(2)].z;
    float _S3905 = 2.0f * _S3903.rows[int(2)].y;
    float _S3906 = 2.0f * _S3903.rows[int(2)].x;
    float _S3907 = 2.0f * _S3903.rows[int(1)].z;
    float _S3908 = 2.0f * - _S3903.rows[int(1)].y;
    float _S3909 = 2.0f * _S3903.rows[int(1)].x;
    float _S3910 = 2.0f * _S3903.rows[int(0)].z;
    float _S3911 = 2.0f * _S3903.rows[int(0)].y;
    float _S3912 = 2.0f * - _S3903.rows[int(0)].x;
    float _S3913 = - _S3909 + _S3911;
    float _S3914 = _S3906 + - _S3910;
    float _S3915 = - _S3905 + _S3907;
    float _S3916 = _S3905 + _S3907;
    float _S3917 = _S3906 + _S3910;
    float _S3918 = _S3909 + _S3911;
    float _S3919 = quat_26.w * (_S3908 + _S3912);
    float _S3920 = quat_26.z * (_S3904 + _S3912);
    float _S3921 = quat_26.y * (_S3904 + _S3908);
    float _S3922 = quat_26.x * _S3913 + quat_26.z * _S3916 + quat_26.y * _S3917 + _S3919 + _S3919;
    float _S3923 = quat_26.x * _S3914 + quat_26.w * _S3916 + quat_26.y * _S3918 + _S3920 + _S3920;
    float _S3924 = quat_26.x * _S3915 + quat_26.w * _S3917 + quat_26.z * _S3918 + _S3921 + _S3921;
    float _S3925 = quat_26.w * _S3913 + quat_26.z * _S3914 + quat_26.y * _S3915;
    float _S3926 = _S3896 + _S3899;
    float _S3927 = 0.5f * - _S3926;
    float _S3928 = _S3894 + _S3898.differential_0.y;
    DiffPair_float_0 _S3929;
    (&_S3929)->primal_0 = _S3508;
    (&_S3929)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3929, _S3928);
    float _S3930 = _S3927 + _S3929.differential_0;
    float _S3931 = _S3895 + _S3900 + _S3902.differential_0.x;
    DiffPair_float_0 _S3932;
    (&_S3932)->primal_0 = _S3506;
    (&_S3932)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S3932, _S3931);
    float _S3933 = _S3927 + _S3932.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3934;
    (&_S3934)->primal_0 = mean_c_20;
    (&_S3934)->differential_0 = _S3652;
    s_bwd_length_impl_1(&_S3934, 0.0f);
    float3  _S3935 = _S3934.differential_0 + _S3655.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S3936;
    (&_S3936)->primal_0 = R_26;
    (&_S3936)->differential_0 = _S3735;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S3937;
    (&_S3937)->primal_0 = mean_25;
    (&_S3937)->differential_0 = _S3652;
    s_bwd_prop_mul_1(&_S3936, &_S3937, _S3935);
    float3  _S3938 = _S3883 + _S3886 + _S3889 + _S3935 + _S3738.differential_0;
    Matrix<float, 3, 3>  _S3939 = _S3884.differential_0 + _S3887.differential_0 + _S3890.differential_0 + _S3936.differential_0 + _S3739;
    float3  _S3940 = make_float3 (_S3933, _S3930, _S3926);
    float4  _S3941 = make_float4 (0.0f);
    *&((&_S3941)->w) = _S3922;
    *&((&_S3941)->z) = _S3923;
    *&((&_S3941)->y) = _S3924;
    *&((&_S3941)->x) = _S3925;
    float4  _S3942 = _S3941;
    float3  _S3943 = _S3885.differential_0 + _S3888.differential_0 + _S3891.differential_0 + _S3937.differential_0 + _S3733;
    *v_mean_8 = _S3943;
    *v_quat_7 = _S3942;
    *v_scale_7 = _S3940;
    *v_hardness_1 = _S3784;
    (*v_sh_coeffs_6)[int(0)] = _S3792;
    (*v_sh_coeffs_6)[int(1)] = _S3793;
    (*v_sh_coeffs_6)[int(2)] = _S3794;
    (*v_sh_coeffs_6)[int(3)] = _S3795;
    (*v_sh_coeffs_6)[int(4)] = _S3796;
    (*v_sh_coeffs_6)[int(5)] = _S3797;
    (*v_sh_coeffs_6)[int(6)] = _S3798;
    (*v_sh_coeffs_6)[int(7)] = _S3799;
    (*v_sh_coeffs_6)[int(8)] = _S3800;
    (*v_sh_coeffs_6)[int(9)] = _S3801;
    (*v_sh_coeffs_6)[int(10)] = _S3802;
    (*v_sh_coeffs_6)[int(11)] = _S3803;
    (*v_sh_coeffs_6)[int(12)] = _S3804;
    (*v_sh_coeffs_6)[int(13)] = _S3805;
    (*v_sh_coeffs_6)[int(14)] = _S3806;
    (*v_sh_coeffs_6)[int(15)] = _S3807;
    (*v_ch_coeffs_1)[int(0)] = _S3789;
    (*v_ch_coeffs_1)[int(1)] = _S3790;
    *v_R_8 = _S3939;
    *v_t_7 = _S3938;
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
        DiffPair_float_0 _S3944 = *dpx_17;
        float _S3945 = val_0 * (*dpy_6).primal_0 / (*dpx_17).primal_0 * dOut_20;
        dpx_17->primal_0 = (*dpx_17).primal_0;
        dpx_17->differential_0 = _S3945;
        float _S3946 = val_0 * (F32_log((_S3944.primal_0))) * dOut_20;
        dpy_6->primal_0 = (*dpy_6).primal_0;
        dpy_6->differential_0 = _S3946;
    }
    return;
}

inline __device__ float evaluate_alpha_opaque_triangle_fast(float2  v0_0, float2  v1_0, float2  v2_0, float2  hardness_6, float2  p_1)
{
    float2  e0_6 = v1_0 - v0_0;
    float2  e1_6 = v2_0 - v1_0;
    float2  e2_2 = v0_0 - v2_0;
    float _S3947 = e0_6.x * e1_6.y - e0_6.y * e1_6.x;
    float se_0 = float((F32_sign((_S3947))));
    float2  _S3948 = p_1 - v0_0;
    float2  _S3949 = normalize_1(e0_6);
    float2  _S3950 = p_1 - v1_0;
    float2  _S3951 = normalize_1(e1_6);
    float2  _S3952 = p_1 - v2_0;
    float2  _S3953 = normalize_1(e2_2);
    float _S3954 = hardness_6.x;
    float _S3955 = 1.0f - clamp_0(hardness_6.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_1 = 1.0f - (1.0f + (F32_max(((F32_max((se_0 * (_S3948.x * _S3949.y - _S3948.y * _S3949.x)), (se_0 * (_S3950.x * _S3951.y - _S3950.y * _S3951.x))))), (se_0 * (_S3952.x * _S3953.y - _S3952.y * _S3953.x)))) / ((F32_abs((_S3947))) / (length_0(e0_6) + length_0(e1_6) + length_0(e2_2)))) * (1.0f - (F32_exp2((-1.0f / _S3955))));
    float _S3956;
    if(a_1 <= 0.0f)
    {
        _S3956 = 0.0f;
    }
    else
    {
        _S3956 = (F32_min(((F32_pow((a_1), (_S3955)))), (0.99900001287460327f)));
    }
    return _S3954 * _S3956;
}

inline __device__ float s_primal_ctx_abs_0(float _S3957)
{
    return (F32_abs((_S3957)));
}

inline __device__ float s_primal_ctx_clamp_0(float _S3958, float _S3959, float _S3960)
{
    return clamp_0(_S3958, _S3959, _S3960);
}

inline __device__ float s_primal_ctx_exp2_0(float _S3961)
{
    return (F32_exp2((_S3961)));
}

inline __device__ float s_primal_ctx_pow_0(float _S3962, float _S3963)
{
    return (F32_pow((_S3962), (_S3963)));
}

inline __device__ void s_bwd_prop_pow_0(DiffPair_float_0 * _S3964, DiffPair_float_0 * _S3965, float _S3966)
{
    _d_pow_0(_S3964, _S3965, _S3966);
    return;
}

inline __device__ void s_bwd_prop_clamp_0(DiffPair_float_0 * _S3967, DiffPair_float_0 * _S3968, DiffPair_float_0 * _S3969, float _S3970)
{
    _d_clamp_0(_S3967, _S3968, _S3969, _S3970);
    return;
}

inline __device__ void s_bwd_prop_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpx_18, float2  _s_dOut_10)
{
    float _S3971 = length_0((*dpx_18).primal_0);
    float2  _S3972 = (*dpx_18).primal_0 * _s_dOut_10;
    float2  _S3973 = make_float2 (1.0f / _S3971) * _s_dOut_10;
    float _S3974 = - ((_S3972.x + _S3972.y) / (_S3971 * _S3971));
    float2  _S3975 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S3976;
    (&_S3976)->primal_0 = (*dpx_18).primal_0;
    (&_S3976)->differential_0 = _S3975;
    s_bwd_length_impl_0(&_S3976, _S3974);
    float2  _S3977 = _S3973 + _S3976.differential_0;
    dpx_18->primal_0 = (*dpx_18).primal_0;
    dpx_18->differential_0 = _S3977;
    return;
}

inline __device__ void s_bwd_normalize_impl_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S3978, float2  _S3979)
{
    s_bwd_prop_normalize_impl_1(_S3978, _S3979);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_0, float2  p_2, float _s_dOut_11)
{
    float2  e0_7 = (*dpv1_0).primal_0 - (*dpv0_0).primal_0;
    float2  e1_7 = (*dpv2_0).primal_0 - (*dpv1_0).primal_0;
    float2  e2_3 = (*dpv0_0).primal_0 - (*dpv2_0).primal_0;
    float _S3980 = e0_7.x;
    float _S3981 = e1_7.y;
    float _S3982 = e0_7.y;
    float _S3983 = e1_7.x;
    float _S3984 = _S3980 * _S3981 - _S3982 * _S3983;
    float se_1 = float((F32_sign((_S3984))));
    float2  _S3985 = p_2 - (*dpv0_0).primal_0;
    float2  _S3986 = normalize_1(e0_7);
    float _S3987 = _S3985.x;
    float _S3988 = _S3986.y;
    float _S3989 = _S3985.y;
    float _S3990 = _S3986.x;
    float de0_0 = se_1 * (_S3987 * _S3988 - _S3989 * _S3990);
    float2  _S3991 = p_2 - (*dpv1_0).primal_0;
    float2  _S3992 = normalize_1(e1_7);
    float _S3993 = _S3991.x;
    float _S3994 = _S3992.y;
    float _S3995 = _S3991.y;
    float _S3996 = _S3992.x;
    float de1_0 = se_1 * (_S3993 * _S3994 - _S3995 * _S3996);
    float2  _S3997 = p_2 - (*dpv2_0).primal_0;
    float2  _S3998 = normalize_1(e2_3);
    float _S3999 = _S3997.x;
    float _S4000 = _S3998.y;
    float _S4001 = _S3997.y;
    float _S4002 = _S3998.x;
    float de2_0 = se_1 * (_S3999 * _S4000 - _S4001 * _S4002);
    float _S4003 = s_primal_ctx_max_0(de0_0, de1_0);
    float _S4004 = s_primal_ctx_max_0(_S4003, de2_0);
    float _S4005 = s_primal_ctx_abs_0(_S3984);
    float _S4006 = length_0(e0_7) + length_0(e1_7) + length_0(e2_3);
    float dmax_0 = _S4005 / _S4006;
    float _S4007 = _S4006 * _S4006;
    float _S4008 = (*dphardness_0).primal_0.x;
    float _S4009 = (*dphardness_0).primal_0.y;
    float _S4010 = dmax_0 * dmax_0;
    float _S4011 = 1.0f + _S4004 / dmax_0;
    float _S4012 = 1.0f - s_primal_ctx_clamp_0(_S4009, 0.00499999988824129f, 0.98000001907348633f);
    float _S4013 = -1.0f / _S4012;
    float _S4014 = _S4012 * _S4012;
    float _S4015 = 1.0f - s_primal_ctx_exp2_0(_S4013);
    float a_2 = 1.0f - _S4011 * _S4015;
    bool _S4016 = a_2 <= 0.0f;
    float _S4017;
    float _S4018;
    if(_S4016)
    {
        _S4017 = 0.0f;
        _S4018 = 0.0f;
    }
    else
    {
        float _S4019 = s_primal_ctx_pow_0(a_2, _S4012);
        _S4017 = s_primal_ctx_min_0(_S4019, 0.99900001287460327f);
        _S4018 = _S4019;
    }
    float _S4020 = _S4008 * _s_dOut_11;
    float _S4021 = _S4017 * _s_dOut_11;
    if(_S4016)
    {
        _S4017 = 0.0f;
        _S4018 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4022;
        (&_S4022)->primal_0 = _S4018;
        (&_S4022)->differential_0 = 0.0f;
        DiffPair_float_0 _S4023;
        (&_S4023)->primal_0 = 0.99900001287460327f;
        (&_S4023)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4022, &_S4023, _S4020);
        DiffPair_float_0 _S4024;
        (&_S4024)->primal_0 = a_2;
        (&_S4024)->differential_0 = 0.0f;
        DiffPair_float_0 _S4025;
        (&_S4025)->primal_0 = _S4012;
        (&_S4025)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4024, &_S4025, _S4022.differential_0);
        _S4017 = _S4024.differential_0;
        _S4018 = _S4025.differential_0;
    }
    float _S4026 = - _S4017;
    float _S4027 = _S4015 * _S4026;
    float _S4028 = - (_S4011 * _S4026);
    DiffPair_float_0 _S4029;
    (&_S4029)->primal_0 = _S4013;
    (&_S4029)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4029, _S4028);
    float _S4030 = - (-1.0f * - (_S4029.differential_0 / _S4014) + _S4018);
    float _S4031 = _S4027 / _S4010;
    float s_diff_dmax_T_0 = _S4004 * - _S4031;
    float _S4032 = dmax_0 * _S4031;
    DiffPair_float_0 _S4033;
    (&_S4033)->primal_0 = _S4009;
    (&_S4033)->differential_0 = 0.0f;
    DiffPair_float_0 _S4034;
    (&_S4034)->primal_0 = 0.00499999988824129f;
    (&_S4034)->differential_0 = 0.0f;
    DiffPair_float_0 _S4035;
    (&_S4035)->primal_0 = 0.98000001907348633f;
    (&_S4035)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4033, &_S4034, &_S4035, _S4030);
    float _S4036 = s_diff_dmax_T_0 / _S4007;
    float _S4037 = _S4005 * - _S4036;
    float _S4038 = _S4006 * _S4036;
    float2  _S4039 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4040;
    (&_S4040)->primal_0 = e2_3;
    (&_S4040)->differential_0 = _S4039;
    s_bwd_length_impl_0(&_S4040, _S4037);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4041;
    (&_S4041)->primal_0 = e1_7;
    (&_S4041)->differential_0 = _S4039;
    s_bwd_length_impl_0(&_S4041, _S4037);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4042;
    (&_S4042)->primal_0 = e0_7;
    (&_S4042)->differential_0 = _S4039;
    s_bwd_length_impl_0(&_S4042, _S4037);
    DiffPair_float_0 _S4043;
    (&_S4043)->primal_0 = _S3984;
    (&_S4043)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4043, _S4038);
    DiffPair_float_0 _S4044;
    (&_S4044)->primal_0 = _S4003;
    (&_S4044)->differential_0 = 0.0f;
    DiffPair_float_0 _S4045;
    (&_S4045)->primal_0 = de2_0;
    (&_S4045)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4044, &_S4045, _S4032);
    DiffPair_float_0 _S4046;
    (&_S4046)->primal_0 = de0_0;
    (&_S4046)->differential_0 = 0.0f;
    DiffPair_float_0 _S4047;
    (&_S4047)->primal_0 = de1_0;
    (&_S4047)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4046, &_S4047, _S4044.differential_0);
    float _S4048 = se_1 * _S4045.differential_0;
    float _S4049 = - _S4048;
    float _S4050 = _S4002 * _S4049;
    float _S4051 = _S4000 * _S4048;
    float2  _S4052 = make_float2 (_S4001 * _S4049, _S3999 * _S4048);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4053;
    (&_S4053)->primal_0 = e2_3;
    (&_S4053)->differential_0 = _S4039;
    s_bwd_normalize_impl_1(&_S4053, _S4052);
    float2  _S4054 = - make_float2 (_S4051, _S4050);
    float _S4055 = se_1 * _S4047.differential_0;
    float _S4056 = - _S4055;
    float _S4057 = _S3996 * _S4056;
    float _S4058 = _S3994 * _S4055;
    float2  _S4059 = make_float2 (_S3995 * _S4056, _S3993 * _S4055);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4060;
    (&_S4060)->primal_0 = e1_7;
    (&_S4060)->differential_0 = _S4039;
    s_bwd_normalize_impl_1(&_S4060, _S4059);
    float2  _S4061 = - make_float2 (_S4058, _S4057);
    float _S4062 = se_1 * _S4046.differential_0;
    float _S4063 = - _S4062;
    float _S4064 = _S3990 * _S4063;
    float _S4065 = _S3988 * _S4062;
    float2  _S4066 = make_float2 (_S3989 * _S4063, _S3987 * _S4062);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4067;
    (&_S4067)->primal_0 = e0_7;
    (&_S4067)->differential_0 = _S4039;
    s_bwd_normalize_impl_1(&_S4067, _S4066);
    float2  _S4068 = - make_float2 (_S4065, _S4064);
    float _S4069 = - _S4043.differential_0;
    float2  _S4070 = _S4040.differential_0 + _S4053.differential_0;
    float2  _S4071 = - _S4070;
    float2  _S4072 = _S4041.differential_0 + _S4060.differential_0 + make_float2 (_S3982 * _S4069, _S3980 * _S4043.differential_0);
    float2  _S4073 = - _S4072;
    float2  _S4074 = _S4042.differential_0 + _S4067.differential_0 + make_float2 (_S3981 * _S4043.differential_0, _S3983 * _S4069);
    float2  _S4075 = - _S4074;
    float2  _S4076 = make_float2 (_S4021, _S4033.differential_0);
    dphardness_0->primal_0 = (*dphardness_0).primal_0;
    dphardness_0->differential_0 = _S4076;
    float2  _S4077 = _S4054 + _S4071 + _S4072;
    dpv2_0->primal_0 = (*dpv2_0).primal_0;
    dpv2_0->differential_0 = _S4077;
    float2  _S4078 = _S4061 + _S4073 + _S4074;
    dpv1_0->primal_0 = (*dpv1_0).primal_0;
    dpv1_0->differential_0 = _S4078;
    float2  _S4079 = _S4068 + _S4070 + _S4075;
    dpv0_0->primal_0 = (*dpv0_0).primal_0;
    dpv0_0->differential_0 = _S4079;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_fast_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4080, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4081, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4082, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4083, float2  _S4084, float _S4085)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_fast_0(_S4080, _S4081, _S4082, _S4083, _S4084, _S4085);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_fast_vjp(float2  v0_1, float2  v1_1, float2  v2_1, float2  hardness_7, float2  p_3, float v_alpha_1, float2  * v_v0_0, float2  * v_v1_0, float2  * v_v2_0, float2  * v_hardness_2)
{
    float2  _S4086 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_0;
    (&dp_v0_0)->primal_0 = v0_1;
    (&dp_v0_0)->differential_0 = _S4086;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_0;
    (&dp_v1_0)->primal_0 = v1_1;
    (&dp_v1_0)->differential_0 = _S4086;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_0;
    (&dp_v2_0)->primal_0 = v2_1;
    (&dp_v2_0)->differential_0 = _S4086;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_0;
    (&dp_hardness_0)->primal_0 = hardness_7;
    (&dp_hardness_0)->differential_0 = _S4086;
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
    float2  _S4087 = p_4 - v0_2;
    float2  _S4088 = p_4 - v1_2;
    float2  _S4089 = p_4 - v2_2;
    float _S4090 = e0_8.x;
    float _S4091 = e1_8.y;
    float _S4092 = e0_8.y;
    float _S4093 = e1_8.x;
    float _S4094 = _S4090 * _S4091 - _S4092 * _S4093;
    float se_2 = float((F32_sign((_S4094))));
    float _S4095 = hardness_8.x;
    float _S4096 = 1.0f - clamp_0(hardness_8.y, 0.00499999988824129f, 0.98000001907348633f);
    float a_3 = 1.0f - (1.0f + float((F32_sign(((F32_max(((F32_max((se_2 * (_S4087.x * _S4092 - _S4087.y * _S4090)), (se_2 * (_S4088.x * _S4091 - _S4088.y * _S4093))))), (se_2 * (_S4089.x * e2_4.y - _S4089.y * e2_4.x)))))))) * (F32_min(((F32_min((length_0(_S4087 - e0_8 * make_float2 (clamp_0(dot_1(_S4087, e0_8) / dot_1(e0_8, e0_8), 0.0f, 1.0f)))), (length_0(_S4088 - e1_8 * make_float2 (clamp_0(dot_1(_S4088, e1_8) / dot_1(e1_8, e1_8), 0.0f, 1.0f))))))), (length_0(_S4089 - e2_4 * make_float2 (clamp_0(dot_1(_S4089, e2_4) / dot_1(e2_4, e2_4), 0.0f, 1.0f)))))) / ((F32_abs((_S4094))) / (length_0(e0_8) + length_0(e1_8) + length_0(e2_4)))) * (1.0f - (F32_exp2((-1.0f / _S4096))));
    float _S4097;
    if(a_3 <= 0.0f)
    {
        _S4097 = 0.0f;
    }
    else
    {
        _S4097 = (F32_min(((F32_pow((a_3), (_S4096)))), (0.99900001287460327f)));
    }
    return _S4095 * _S4097;
}

inline __device__ float s_primal_ctx_dot_1(float2  _S4098, float2  _S4099)
{
    return dot_1(_S4098, _S4099);
}

inline __device__ void s_bwd_prop_dot_1(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4100, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4101, float _S4102)
{
    _d_dot_1(_S4100, _S4101, _S4102);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv0_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv1_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dpv2_1, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_1, float2  p_5, float _s_dOut_12)
{
    float2  e0_9 = (*dpv1_1).primal_0 - (*dpv0_1).primal_0;
    float2  e1_9 = (*dpv2_1).primal_0 - (*dpv1_1).primal_0;
    float2  e2_5 = (*dpv0_1).primal_0 - (*dpv2_1).primal_0;
    float2  _S4103 = p_5 - (*dpv0_1).primal_0;
    float _S4104 = s_primal_ctx_dot_1(_S4103, e0_9);
    float _S4105 = s_primal_ctx_dot_1(e0_9, e0_9);
    float _S4106 = _S4104 / _S4105;
    float _S4107 = _S4105 * _S4105;
    float _S4108 = s_primal_ctx_clamp_0(_S4106, 0.0f, 1.0f);
    float2  _S4109 = make_float2 (_S4108);
    float2  _S4110 = _S4103 - e0_9 * make_float2 (_S4108);
    float _S4111 = length_0(_S4110);
    float2  _S4112 = p_5 - (*dpv1_1).primal_0;
    float _S4113 = s_primal_ctx_dot_1(_S4112, e1_9);
    float _S4114 = s_primal_ctx_dot_1(e1_9, e1_9);
    float _S4115 = _S4113 / _S4114;
    float _S4116 = _S4114 * _S4114;
    float _S4117 = s_primal_ctx_clamp_0(_S4115, 0.0f, 1.0f);
    float2  _S4118 = make_float2 (_S4117);
    float2  _S4119 = _S4112 - e1_9 * make_float2 (_S4117);
    float _S4120 = length_0(_S4119);
    float2  _S4121 = p_5 - (*dpv2_1).primal_0;
    float _S4122 = s_primal_ctx_dot_1(_S4121, e2_5);
    float _S4123 = s_primal_ctx_dot_1(e2_5, e2_5);
    float _S4124 = _S4122 / _S4123;
    float _S4125 = _S4123 * _S4123;
    float _S4126 = s_primal_ctx_clamp_0(_S4124, 0.0f, 1.0f);
    float2  _S4127 = make_float2 (_S4126);
    float2  _S4128 = _S4121 - e2_5 * make_float2 (_S4126);
    float _S4129 = length_0(_S4128);
    float _S4130 = e0_9.x;
    float _S4131 = e1_9.y;
    float _S4132 = e0_9.y;
    float _S4133 = e1_9.x;
    float _S4134 = _S4130 * _S4131 - _S4132 * _S4133;
    float se_3 = float((F32_sign((_S4134))));
    float _S4135 = _S4103.x;
    float _S4136 = _S4103.y;
    float s0_0 = se_3 * (_S4135 * _S4132 - _S4136 * _S4130);
    float _S4137 = _S4112.x;
    float _S4138 = _S4112.y;
    float s1_0 = se_3 * (_S4137 * _S4131 - _S4138 * _S4133);
    float _S4139 = _S4121.x;
    float _S4140 = e2_5.y;
    float _S4141 = _S4121.y;
    float _S4142 = e2_5.x;
    float s2_0 = se_3 * (_S4139 * _S4140 - _S4141 * _S4142);
    float _S4143 = s_primal_ctx_max_0(s0_0, s1_0);
    float sv_0 = float((F32_sign((s_primal_ctx_max_0(_S4143, s2_0)))));
    float _S4144 = s_primal_ctx_min_0(_S4111, _S4120);
    float dv_0 = sv_0 * s_primal_ctx_min_0(_S4144, _S4129);
    float _S4145 = s_primal_ctx_abs_0(_S4134);
    float _S4146 = length_0(e0_9) + length_0(e1_9) + length_0(e2_5);
    float dmax_1 = _S4145 / _S4146;
    float _S4147 = _S4146 * _S4146;
    float _S4148 = (*dphardness_1).primal_0.x;
    float _S4149 = (*dphardness_1).primal_0.y;
    float _S4150 = dmax_1 * dmax_1;
    float _S4151 = 1.0f + dv_0 / dmax_1;
    float _S4152 = 1.0f - s_primal_ctx_clamp_0(_S4149, 0.00499999988824129f, 0.98000001907348633f);
    float _S4153 = -1.0f / _S4152;
    float _S4154 = _S4152 * _S4152;
    float _S4155 = 1.0f - s_primal_ctx_exp2_0(_S4153);
    float a_4 = 1.0f - _S4151 * _S4155;
    bool _S4156 = a_4 <= 0.0f;
    float _S4157;
    float _S4158;
    if(_S4156)
    {
        _S4157 = 0.0f;
        _S4158 = 0.0f;
    }
    else
    {
        float _S4159 = s_primal_ctx_pow_0(a_4, _S4152);
        _S4157 = s_primal_ctx_min_0(_S4159, 0.99900001287460327f);
        _S4158 = _S4159;
    }
    float _S4160 = _S4148 * _s_dOut_12;
    float _S4161 = _S4157 * _s_dOut_12;
    if(_S4156)
    {
        _S4157 = 0.0f;
        _S4158 = 0.0f;
    }
    else
    {
        DiffPair_float_0 _S4162;
        (&_S4162)->primal_0 = _S4158;
        (&_S4162)->differential_0 = 0.0f;
        DiffPair_float_0 _S4163;
        (&_S4163)->primal_0 = 0.99900001287460327f;
        (&_S4163)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S4162, &_S4163, _S4160);
        DiffPair_float_0 _S4164;
        (&_S4164)->primal_0 = a_4;
        (&_S4164)->differential_0 = 0.0f;
        DiffPair_float_0 _S4165;
        (&_S4165)->primal_0 = _S4152;
        (&_S4165)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S4164, &_S4165, _S4162.differential_0);
        _S4157 = _S4164.differential_0;
        _S4158 = _S4165.differential_0;
    }
    float _S4166 = - _S4157;
    float _S4167 = _S4155 * _S4166;
    float _S4168 = - (_S4151 * _S4166);
    DiffPair_float_0 _S4169;
    (&_S4169)->primal_0 = _S4153;
    (&_S4169)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4169, _S4168);
    float _S4170 = - (-1.0f * - (_S4169.differential_0 / _S4154) + _S4158);
    float _S4171 = _S4167 / _S4150;
    float s_diff_dmax_T_1 = dv_0 * - _S4171;
    float s_diff_dv_T_0 = dmax_1 * _S4171;
    DiffPair_float_0 _S4172;
    (&_S4172)->primal_0 = _S4149;
    (&_S4172)->differential_0 = 0.0f;
    DiffPair_float_0 _S4173;
    (&_S4173)->primal_0 = 0.00499999988824129f;
    (&_S4173)->differential_0 = 0.0f;
    DiffPair_float_0 _S4174;
    (&_S4174)->primal_0 = 0.98000001907348633f;
    (&_S4174)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4172, &_S4173, &_S4174, _S4170);
    float _S4175 = s_diff_dmax_T_1 / _S4147;
    float _S4176 = _S4145 * - _S4175;
    float _S4177 = _S4146 * _S4175;
    float2  _S4178 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4179;
    (&_S4179)->primal_0 = e2_5;
    (&_S4179)->differential_0 = _S4178;
    s_bwd_length_impl_0(&_S4179, _S4176);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4180;
    (&_S4180)->primal_0 = e1_9;
    (&_S4180)->differential_0 = _S4178;
    s_bwd_length_impl_0(&_S4180, _S4176);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4181;
    (&_S4181)->primal_0 = e0_9;
    (&_S4181)->differential_0 = _S4178;
    s_bwd_length_impl_0(&_S4181, _S4176);
    DiffPair_float_0 _S4182;
    (&_S4182)->primal_0 = _S4134;
    (&_S4182)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4182, _S4177);
    float _S4183 = sv_0 * s_diff_dv_T_0;
    DiffPair_float_0 _S4184;
    (&_S4184)->primal_0 = _S4144;
    (&_S4184)->differential_0 = 0.0f;
    DiffPair_float_0 _S4185;
    (&_S4185)->primal_0 = _S4129;
    (&_S4185)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4184, &_S4185, _S4183);
    DiffPair_float_0 _S4186;
    (&_S4186)->primal_0 = _S4111;
    (&_S4186)->differential_0 = 0.0f;
    DiffPair_float_0 _S4187;
    (&_S4187)->primal_0 = _S4120;
    (&_S4187)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4186, &_S4187, _S4184.differential_0);
    DiffPair_float_0 _S4188;
    (&_S4188)->primal_0 = _S4143;
    (&_S4188)->differential_0 = 0.0f;
    DiffPair_float_0 _S4189;
    (&_S4189)->primal_0 = s2_0;
    (&_S4189)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4188, &_S4189, 0.0f);
    DiffPair_float_0 _S4190;
    (&_S4190)->primal_0 = s0_0;
    (&_S4190)->differential_0 = 0.0f;
    DiffPair_float_0 _S4191;
    (&_S4191)->primal_0 = s1_0;
    (&_S4191)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4190, &_S4191, _S4188.differential_0);
    float _S4192 = se_3 * _S4189.differential_0;
    float _S4193 = - _S4192;
    float _S4194 = _S4141 * _S4193;
    float _S4195 = _S4142 * _S4193;
    float _S4196 = _S4139 * _S4192;
    float _S4197 = _S4140 * _S4192;
    float _S4198 = se_3 * _S4191.differential_0;
    float _S4199 = - _S4198;
    float _S4200 = _S4133 * _S4199;
    float _S4201 = _S4131 * _S4198;
    float _S4202 = se_3 * _S4190.differential_0;
    float _S4203 = - _S4202;
    float _S4204 = _S4130 * _S4203;
    float _S4205 = _S4132 * _S4202;
    float _S4206 = - _S4182.differential_0;
    float _S4207 = _S4138 * _S4199 + _S4132 * _S4206;
    float _S4208 = _S4135 * _S4202 + _S4133 * _S4206;
    float _S4209 = _S4137 * _S4198 + _S4130 * _S4182.differential_0;
    float _S4210 = _S4136 * _S4203 + _S4131 * _S4182.differential_0;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4211;
    (&_S4211)->primal_0 = _S4128;
    (&_S4211)->differential_0 = _S4178;
    s_bwd_length_impl_0(&_S4211, _S4185.differential_0);
    float2  _S4212 = - _S4211.differential_0;
    float2  _S4213 = e2_5 * _S4212;
    float2  _S4214 = _S4127 * _S4212;
    float _S4215 = _S4213.x + _S4213.y;
    DiffPair_float_0 _S4216;
    (&_S4216)->primal_0 = _S4124;
    (&_S4216)->differential_0 = 0.0f;
    DiffPair_float_0 _S4217;
    (&_S4217)->primal_0 = 0.0f;
    (&_S4217)->differential_0 = 0.0f;
    DiffPair_float_0 _S4218;
    (&_S4218)->primal_0 = 1.0f;
    (&_S4218)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4216, &_S4217, &_S4218, _S4215);
    float _S4219 = _S4216.differential_0 / _S4125;
    float _S4220 = _S4122 * - _S4219;
    float _S4221 = _S4123 * _S4219;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4222;
    (&_S4222)->primal_0 = e2_5;
    (&_S4222)->differential_0 = _S4178;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4223;
    (&_S4223)->primal_0 = e2_5;
    (&_S4223)->differential_0 = _S4178;
    s_bwd_prop_dot_1(&_S4222, &_S4223, _S4220);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4224;
    (&_S4224)->primal_0 = _S4121;
    (&_S4224)->differential_0 = _S4178;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4225;
    (&_S4225)->primal_0 = e2_5;
    (&_S4225)->differential_0 = _S4178;
    s_bwd_prop_dot_1(&_S4224, &_S4225, _S4221);
    float2  _S4226 = - (_S4211.differential_0 + _S4224.differential_0 + make_float2 (_S4197, _S4195));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4227;
    (&_S4227)->primal_0 = _S4119;
    (&_S4227)->differential_0 = _S4178;
    s_bwd_length_impl_0(&_S4227, _S4187.differential_0);
    float2  _S4228 = - _S4227.differential_0;
    float2  _S4229 = e1_9 * _S4228;
    float2  _S4230 = _S4118 * _S4228;
    float _S4231 = _S4229.x + _S4229.y;
    DiffPair_float_0 _S4232;
    (&_S4232)->primal_0 = _S4115;
    (&_S4232)->differential_0 = 0.0f;
    DiffPair_float_0 _S4233;
    (&_S4233)->primal_0 = 0.0f;
    (&_S4233)->differential_0 = 0.0f;
    DiffPair_float_0 _S4234;
    (&_S4234)->primal_0 = 1.0f;
    (&_S4234)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4232, &_S4233, &_S4234, _S4231);
    float _S4235 = _S4232.differential_0 / _S4116;
    float _S4236 = _S4113 * - _S4235;
    float _S4237 = _S4114 * _S4235;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4238;
    (&_S4238)->primal_0 = e1_9;
    (&_S4238)->differential_0 = _S4178;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4239;
    (&_S4239)->primal_0 = e1_9;
    (&_S4239)->differential_0 = _S4178;
    s_bwd_prop_dot_1(&_S4238, &_S4239, _S4236);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4240;
    (&_S4240)->primal_0 = _S4112;
    (&_S4240)->differential_0 = _S4178;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4241;
    (&_S4241)->primal_0 = e1_9;
    (&_S4241)->differential_0 = _S4178;
    s_bwd_prop_dot_1(&_S4240, &_S4241, _S4237);
    float2  _S4242 = - (_S4227.differential_0 + _S4240.differential_0 + make_float2 (_S4201, _S4200));
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4243;
    (&_S4243)->primal_0 = _S4110;
    (&_S4243)->differential_0 = _S4178;
    s_bwd_length_impl_0(&_S4243, _S4186.differential_0);
    float2  _S4244 = - _S4243.differential_0;
    float2  _S4245 = e0_9 * _S4244;
    float2  _S4246 = _S4109 * _S4244;
    float _S4247 = _S4245.x + _S4245.y;
    DiffPair_float_0 _S4248;
    (&_S4248)->primal_0 = _S4106;
    (&_S4248)->differential_0 = 0.0f;
    DiffPair_float_0 _S4249;
    (&_S4249)->primal_0 = 0.0f;
    (&_S4249)->differential_0 = 0.0f;
    DiffPair_float_0 _S4250;
    (&_S4250)->primal_0 = 1.0f;
    (&_S4250)->differential_0 = 0.0f;
    s_bwd_prop_clamp_0(&_S4248, &_S4249, &_S4250, _S4247);
    float _S4251 = _S4248.differential_0 / _S4107;
    float _S4252 = _S4104 * - _S4251;
    float _S4253 = _S4105 * _S4251;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4254;
    (&_S4254)->primal_0 = e0_9;
    (&_S4254)->differential_0 = _S4178;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4255;
    (&_S4255)->primal_0 = e0_9;
    (&_S4255)->differential_0 = _S4178;
    s_bwd_prop_dot_1(&_S4254, &_S4255, _S4252);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4256;
    (&_S4256)->primal_0 = _S4103;
    (&_S4256)->differential_0 = _S4178;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4257;
    (&_S4257)->primal_0 = e0_9;
    (&_S4257)->differential_0 = _S4178;
    s_bwd_prop_dot_1(&_S4256, &_S4257, _S4253);
    float2  _S4258 = - (_S4243.differential_0 + _S4256.differential_0 + make_float2 (_S4205, _S4204));
    float2  _S4259 = _S4179.differential_0 + _S4214 + _S4223.differential_0 + _S4222.differential_0 + _S4225.differential_0 + make_float2 (_S4194, _S4196);
    float2  _S4260 = - _S4259;
    float2  _S4261 = _S4180.differential_0 + _S4230 + _S4239.differential_0 + _S4238.differential_0 + _S4241.differential_0 + make_float2 (_S4207, _S4209);
    float2  _S4262 = - _S4261;
    float2  _S4263 = _S4181.differential_0 + _S4246 + _S4255.differential_0 + _S4254.differential_0 + _S4257.differential_0 + make_float2 (_S4210, _S4208);
    float2  _S4264 = - _S4263;
    float2  _S4265 = make_float2 (_S4161, _S4172.differential_0);
    dphardness_1->primal_0 = (*dphardness_1).primal_0;
    dphardness_1->differential_0 = _S4265;
    float2  _S4266 = _S4226 + _S4260 + _S4261;
    dpv2_1->primal_0 = (*dpv2_1).primal_0;
    dpv2_1->differential_0 = _S4266;
    float2  _S4267 = _S4242 + _S4262 + _S4263;
    dpv1_1->primal_0 = (*dpv1_1).primal_0;
    dpv1_1->differential_0 = _S4267;
    float2  _S4268 = _S4258 + _S4259 + _S4264;
    dpv0_1->primal_0 = (*dpv0_1).primal_0;
    dpv0_1->differential_0 = _S4268;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_precise_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4269, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4270, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4271, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4272, float2  _S4273, float _S4274)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_precise_0(_S4269, _S4270, _S4271, _S4272, _S4273, _S4274);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_precise_vjp(float2  v0_3, float2  v1_3, float2  v2_3, float2  hardness_9, float2  p_6, float v_alpha_2, float2  * v_v0_1, float2  * v_v1_1, float2  * v_v2_1, float2  * v_hardness_3)
{
    float2  _S4275 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_1;
    (&dp_v0_1)->primal_0 = v0_3;
    (&dp_v0_1)->differential_0 = _S4275;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_1;
    (&dp_v1_1)->primal_0 = v1_3;
    (&dp_v1_1)->differential_0 = _S4275;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_1;
    (&dp_v2_1)->primal_0 = v2_3;
    (&dp_v2_1)->differential_0 = _S4275;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_1;
    (&dp_hardness_1)->primal_0 = hardness_9;
    (&dp_hardness_1)->differential_0 = _S4275;
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
    float _S4276 = 0.3333333432674408f * dpdepth_1;
    float3  _S4277 = make_float3 (0.3333333432674408f) * dpcolor_0;
    float3  _S4278 = make_float3 (0.0f);
    float3  _S4279 = _S4278;
    *&((&_S4279)->z) = _S4276;
    *&((&_S4279)->y) = _S4276;
    *&((&_S4279)->x) = _S4276;
    dpdepths_0->primal_0 = (*dpdepths_0).primal_0;
    dpdepths_0->differential_0 = _S4279;
    FixedArray<float3 , 3>  _S4280;
    _S4280[int(0)] = _S4278;
    _S4280[int(1)] = _S4278;
    _S4280[int(2)] = _S4278;
    _S4280[int(2)] = _S4277;
    _S4280[int(1)] = _S4277;
    _S4280[int(0)] = _S4277;
    dpcolors_0->primal_0 = dpcolors_0->primal_0;
    dpcolors_0->differential_0 = _S4280;
    float2  _S4281 = make_float2 (0.0f);
    dpv2_2->primal_0 = (*dpv2_2).primal_0;
    dpv2_2->differential_0 = _S4281;
    dpv1_2->primal_0 = (*dpv1_2).primal_0;
    dpv1_2->differential_0 = _S4281;
    dpv0_2->primal_0 = (*dpv0_2).primal_0;
    dpv0_2->differential_0 = _S4281;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_0(DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4282, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4283, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S4284, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S4285, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S4286, float2  _S4287, float3  _S4288, float _S4289)
{
    s_bwd_prop_evaluate_color_opaque_triangle_0(_S4282, _S4283, _S4284, _S4285, _S4286, _S4287, _S4288, _S4289);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(float2  v0_5, float2  v1_5, float2  v2_5, FixedArray<float3 , 3>  * colors_1, float3  depths_1, float2  p_9, float3  v_color_0, float v_depth_8, float2  * v_v0_2, float2  * v_v1_2, float2  * v_v2_2, FixedArray<float3 , 3>  * v_colors_0, float3  * v_depths_0)
{
    float2  _S4290 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v0_2;
    (&dp_v0_2)->primal_0 = v0_5;
    (&dp_v0_2)->differential_0 = _S4290;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v1_2;
    (&dp_v1_2)->primal_0 = v1_5;
    (&dp_v1_2)->differential_0 = _S4290;
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_v2_2;
    (&dp_v2_2)->primal_0 = v2_5;
    (&dp_v2_2)->differential_0 = _S4290;
    float3  _S4291 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S4292 = { _S4291, _S4291, _S4291 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_colors_0;
    (&dp_colors_0)->primal_0 = *colors_1;
    (&dp_colors_0)->differential_0 = _S4292;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_depths_0;
    (&dp_depths_0)->primal_0 = depths_1;
    (&dp_depths_0)->differential_0 = _S4291;
    s_bwd_evaluate_color_opaque_triangle_0(&dp_v0_2, &dp_v1_2, &dp_v2_2, &dp_colors_0, &dp_depths_0, p_9, v_color_0, v_depth_8);
    *v_v0_2 = dp_v0_2.differential_0;
    *v_v1_2 = dp_v2_2.differential_0;
    *v_v2_2 = dp_v1_2.differential_0;
    *v_colors_0 = (&dp_colors_0)->differential_0;
    *v_depths_0 = dp_depths_0.differential_0;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp(float3  mean_26, float4  quat_27, float3  scale_26, float2  hardness_10, FixedArray<float3 , 16>  * sh_coeffs_21, FixedArray<float3 , 2>  * ch_coeffs_6, Matrix<float, 3, 3>  R_27, float3  t_24, float fx_30, float fy_30, float cx_25, float cy_25, FixedArray<float, 10>  * dist_coeffs_34, uint image_width_21, uint image_height_21, float near_plane_14, float far_plane_14, int4  * aabb_xyxy_14, float * depth_16, FixedArray<float3 , 3>  * verts_0, FixedArray<float3 , 3>  * rgbs_0, float3  * normal_4)
{
    for(;;)
    {
        float3  mean_c_21 = mul_0(R_27, mean_26) + t_24;
        float _S4293 = mean_c_21.z;
        bool _S4294;
        if(_S4293 < near_plane_14)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = _S4293 > far_plane_14;
        }
        if(_S4294)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4295 = scale_26.x;
        float sx_5 = (F32_exp((_S4295)));
        float _S4296 = scale_26.y;
        float sy_5 = (F32_exp((_S4296)));
        float sz_7 = scale_26.z - 0.5f * (_S4295 + _S4296);
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
        Matrix<float, 3, 3>  _S4297 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_27 + z2_48), 2.0f * (xy_27 + wz_27), 2.0f * (xz_27 - wy_27), 2.0f * (xy_27 - wz_27), 1.0f - 2.0f * (x2_27 + z2_48), 2.0f * (yz_27 + wx_27), 2.0f * (xz_27 + wy_27), 2.0f * (yz_27 - wx_27), 1.0f - 2.0f * (x2_27 + y2_27)));
        float3  vert0_3 = mul_0(_S4297, make_float3 (sx_5, 0.0f, 0.0f)) + mean_26;
        float3  vert1_3 = mul_0(_S4297, make_float3 (sx_5 * (-0.5f + sz_7), sy_5, 0.0f)) + mean_26;
        float3  vert2_3 = mul_0(_S4297, make_float3 (sx_5 * (-0.5f - sz_7), - sy_5, 0.0f)) + mean_26;
        float3  vert0_c_6 = mul_0(R_27, vert0_3) + t_24;
        float3  vert1_c_6 = mul_0(R_27, vert1_3) + t_24;
        float3  vert2_c_6 = mul_0(R_27, vert2_3) + t_24;
        float _S4298 = vert0_c_6.z;
        float _S4299 = vert1_c_6.z;
        float _S4300 = vert2_c_6.z;
        if(_S4298 < near_plane_14)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = _S4298 > far_plane_14;
        }
        if(_S4294)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = _S4299 < near_plane_14;
        }
        if(_S4294)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = _S4299 > far_plane_14;
        }
        if(_S4294)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = _S4300 < near_plane_14;
        }
        if(_S4294)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = _S4300 > far_plane_14;
        }
        if(_S4294)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_4;
        for(;;)
        {
            float2  uv0_5 = float2 {vert0_c_6.x, vert0_c_6.y} / make_float2 (_S4298);
            if(_S4298 < 0.0f)
            {
                _S4294 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4301 = camera_distortion_jac_0(uv0_5, dist_coeffs_34);
                _S4294 = !((F32_min((determinant_0(_S4301)), ((F32_min((_S4301.rows[int(0)].x), (_S4301.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4294)
            {
                uv0_4 = uv0_5;
                break;
            }
            float u_62 = uv0_5.x;
            float v_62 = uv0_5.y;
            float r2_62 = u_62 * u_62 + v_62 * v_62;
            float2  _S4302 = uv0_5 * make_float2 (1.0f + r2_62 * ((*dist_coeffs_34)[int(0)] + r2_62 * ((*dist_coeffs_34)[int(1)] + r2_62 * ((*dist_coeffs_34)[int(2)] + r2_62 * (*dist_coeffs_34)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_34)[int(4)] * u_62 * v_62 + (*dist_coeffs_34)[int(5)] * (r2_62 + 2.0f * u_62 * u_62) + (*dist_coeffs_34)[int(6)] * r2_62, 2.0f * (*dist_coeffs_34)[int(5)] * u_62 * v_62 + (*dist_coeffs_34)[int(4)] * (r2_62 + 2.0f * v_62 * v_62) + (*dist_coeffs_34)[int(7)] * r2_62);
            float2  _S4303 = _S4302 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4302.x + (*dist_coeffs_34)[int(9)] * _S4302.y, 0.0f);
            uv0_4 = make_float2 (fx_30 * _S4303.x + cx_25, fy_30 * _S4303.y + cy_25);
            break;
        }
        float2  uv1_4;
        bool all_valid_12 = true & (!_S4294);
        for(;;)
        {
            float2  uv1_5 = float2 {vert1_c_6.x, vert1_c_6.y} / make_float2 (_S4299);
            if(_S4299 < 0.0f)
            {
                _S4294 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4304 = camera_distortion_jac_0(uv1_5, dist_coeffs_34);
                _S4294 = !((F32_min((determinant_0(_S4304)), ((F32_min((_S4304.rows[int(0)].x), (_S4304.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4294)
            {
                uv1_4 = uv1_5;
                break;
            }
            float u_63 = uv1_5.x;
            float v_63 = uv1_5.y;
            float r2_63 = u_63 * u_63 + v_63 * v_63;
            float2  _S4305 = uv1_5 * make_float2 (1.0f + r2_63 * ((*dist_coeffs_34)[int(0)] + r2_63 * ((*dist_coeffs_34)[int(1)] + r2_63 * ((*dist_coeffs_34)[int(2)] + r2_63 * (*dist_coeffs_34)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_34)[int(4)] * u_63 * v_63 + (*dist_coeffs_34)[int(5)] * (r2_63 + 2.0f * u_63 * u_63) + (*dist_coeffs_34)[int(6)] * r2_63, 2.0f * (*dist_coeffs_34)[int(5)] * u_63 * v_63 + (*dist_coeffs_34)[int(4)] * (r2_63 + 2.0f * v_63 * v_63) + (*dist_coeffs_34)[int(7)] * r2_63);
            float2  _S4306 = _S4305 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4305.x + (*dist_coeffs_34)[int(9)] * _S4305.y, 0.0f);
            uv1_4 = make_float2 (fx_30 * _S4306.x + cx_25, fy_30 * _S4306.y + cy_25);
            break;
        }
        float2  uv2_4;
        bool all_valid_13 = all_valid_12 & (!_S4294);
        for(;;)
        {
            float2  uv2_5 = float2 {vert2_c_6.x, vert2_c_6.y} / make_float2 (_S4300);
            if(_S4300 < 0.0f)
            {
                _S4294 = true;
            }
            else
            {
                Matrix<float, 2, 2>  _S4307 = camera_distortion_jac_0(uv2_5, dist_coeffs_34);
                _S4294 = !((F32_min((determinant_0(_S4307)), ((F32_min((_S4307.rows[int(0)].x), (_S4307.rows[int(1)].y)))))) > 0.0f);
            }
            if(_S4294)
            {
                uv2_4 = uv2_5;
                break;
            }
            float u_64 = uv2_5.x;
            float v_64 = uv2_5.y;
            float r2_64 = u_64 * u_64 + v_64 * v_64;
            float2  _S4308 = uv2_5 * make_float2 (1.0f + r2_64 * ((*dist_coeffs_34)[int(0)] + r2_64 * ((*dist_coeffs_34)[int(1)] + r2_64 * ((*dist_coeffs_34)[int(2)] + r2_64 * (*dist_coeffs_34)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_34)[int(4)] * u_64 * v_64 + (*dist_coeffs_34)[int(5)] * (r2_64 + 2.0f * u_64 * u_64) + (*dist_coeffs_34)[int(6)] * r2_64, 2.0f * (*dist_coeffs_34)[int(5)] * u_64 * v_64 + (*dist_coeffs_34)[int(4)] * (r2_64 + 2.0f * v_64 * v_64) + (*dist_coeffs_34)[int(7)] * r2_64);
            float2  _S4309 = _S4308 + make_float2 ((*dist_coeffs_34)[int(8)] * _S4308.x + (*dist_coeffs_34)[int(9)] * _S4308.y, 0.0f);
            uv2_4 = make_float2 (fx_30 * _S4309.x + cx_25, fy_30 * _S4309.y + cy_25);
            break;
        }
        if(!(all_valid_13 & (!_S4294)))
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_10 = uv1_4 - uv0_4;
        float2  e1_10 = uv2_4 - uv1_4;
        float offset_4 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_10.y))))) - 1.0f) * ((F32_abs((e0_10.x * e1_10.y - e0_10.y * e1_10.x))) / (length_0(e0_10) + length_0(e1_10) + length_0(uv0_4 - uv2_4)));
        float _S4310 = uv0_4.x;
        float _S4311 = uv1_4.x;
        float _S4312 = uv2_4.x;
        float xmax_7 = (F32_max(((F32_max((_S4310), (_S4311)))), (_S4312))) + offset_4;
        float xmin_7 = (F32_min(((F32_min((_S4310), (_S4311)))), (_S4312))) - offset_4;
        float _S4313 = uv0_4.y;
        float _S4314 = uv1_4.y;
        float _S4315 = uv2_4.y;
        float ymax_7 = (F32_max(((F32_max((_S4313), (_S4314)))), (_S4315))) + offset_4;
        float ymin_7 = (F32_min(((F32_min((_S4313), (_S4314)))), (_S4315))) - offset_4;
        if(xmax_7 <= 0.0f)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = xmin_7 >= float(image_width_21);
        }
        if(_S4294)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = ymax_7 <= 0.0f;
        }
        if(_S4294)
        {
            _S4294 = true;
        }
        else
        {
            _S4294 = ymin_7 >= float(image_height_21);
        }
        if(_S4294)
        {
            _S4294 = true;
        }
        else
        {
            if(_S4293 <= 0.0f)
            {
                if(xmin_7 <= 0.0f)
                {
                    _S4294 = xmax_7 >= float(image_width_21);
                }
                else
                {
                    _S4294 = false;
                }
                if(_S4294)
                {
                    _S4294 = true;
                }
                else
                {
                    if(ymin_7 <= 0.0f)
                    {
                        _S4294 = ymax_7 >= float(image_width_21);
                    }
                    else
                    {
                        _S4294 = false;
                    }
                }
            }
            else
            {
                _S4294 = false;
            }
        }
        if(_S4294)
        {
            *aabb_xyxy_14 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_14 = make_int4 (int((F32_floor((xmin_7)))), int((F32_floor((ymin_7)))), int((F32_ceil((xmax_7)))), int((F32_ceil((ymax_7)))));
        *depth_16 = (F32_log((length_1(vert0_c_6 + vert1_c_6 + vert2_c_6) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4316 = mean_26 - - mul_0(transpose_0(R_27), t_24);
        float _S4317 = _S4316.x;
        float _S4318 = _S4316.y;
        float _S4319 = _S4316.z;
        float norm_14 = (F32_sqrt((_S4317 * _S4317 + _S4318 * _S4318 + _S4319 * _S4319)));
        float x_54 = _S4317 / norm_14;
        float y_24 = _S4318 / norm_14;
        float z_21 = _S4319 / norm_14;
        float z2_49 = z_21 * z_21;
        float fTmp0B_21 = -1.09254848957061768f * z_21;
        float fC1_21 = x_54 * x_54 - y_24 * y_24;
        float fS1_21 = 2.0f * x_54 * y_24;
        float fTmp0C_21 = -2.28522896766662598f * z2_49 + 0.4570457935333252f;
        float fTmp1B_21 = 1.44530570507049561f * z_21;
        float3  color_7 = make_float3 (0.282094806432724f) * (*sh_coeffs_21)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_24) * (*sh_coeffs_21)[int(1)] + make_float3 (z_21) * (*sh_coeffs_21)[int(2)] - make_float3 (x_54) * (*sh_coeffs_21)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_21) * (*sh_coeffs_21)[int(4)] + make_float3 (fTmp0B_21 * y_24) * (*sh_coeffs_21)[int(5)] + make_float3 (0.94617468118667603f * z2_49 - 0.31539157032966614f) * (*sh_coeffs_21)[int(6)] + make_float3 (fTmp0B_21 * x_54) * (*sh_coeffs_21)[int(7)] + make_float3 (0.54627424478530884f * fC1_21) * (*sh_coeffs_21)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_54 * fS1_21 + y_24 * fC1_21)) * (*sh_coeffs_21)[int(9)] + make_float3 (fTmp1B_21 * fS1_21) * (*sh_coeffs_21)[int(10)] + make_float3 (fTmp0C_21 * y_24) * (*sh_coeffs_21)[int(11)] + make_float3 (z_21 * (1.86588168144226074f * z2_49 - 1.11952900886535645f)) * (*sh_coeffs_21)[int(12)] + make_float3 (fTmp0C_21 * x_54) * (*sh_coeffs_21)[int(13)] + make_float3 (fTmp1B_21 * fC1_21) * (*sh_coeffs_21)[int(14)] + make_float3 (-0.59004360437393188f * (x_54 * fC1_21 - y_24 * fS1_21)) * (*sh_coeffs_21)[int(15)]);
        float3  _S4320 = make_float3 (0.0f);
        (*rgbs_0)[int(0)] = max_0(color_7 + (*ch_coeffs_6)[int(0)] + make_float3 (0.5f), _S4320);
        float3  _S4321 = color_7 - (*ch_coeffs_6)[int(0)] * make_float3 (0.5f);
        float3  _S4322 = (*ch_coeffs_6)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_0)[int(1)] = max_0(_S4321 + _S4322 + make_float3 (0.5f), _S4320);
        (*rgbs_0)[int(2)] = max_0(_S4321 - _S4322 + make_float3 (0.5f), _S4320);
        (*verts_0)[int(0)] = vert0_3;
        (*verts_0)[int(1)] = vert1_3;
        (*verts_0)[int(2)] = vert2_3;
        float3  _S4323 = normalize_0(cross_0(vert1_c_6 - vert0_c_6, vert2_c_6 - vert0_c_6));
        *normal_4 = _S4323 * make_float3 (float(- (F32_sign((dot_0(_S4323, mean_c_21))))));
        break;
    }
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye(float3  mean_27, float4  quat_28, float3  scale_27, float2  hardness_11, FixedArray<float3 , 16>  * sh_coeffs_22, FixedArray<float3 , 2>  * ch_coeffs_7, Matrix<float, 3, 3>  R_28, float3  t_25, float fx_31, float fy_31, float cx_26, float cy_26, FixedArray<float, 10>  * dist_coeffs_35, uint image_width_22, uint image_height_22, float near_plane_15, float far_plane_15, int4  * aabb_xyxy_15, float * depth_17, FixedArray<float3 , 3>  * verts_1, FixedArray<float3 , 3>  * rgbs_1, float3  * normal_5)
{
    bool _S4324;
    bool _S4325;
    bool _S4326;
    for(;;)
    {
        float3  mean_c_22 = mul_0(R_28, mean_27) + t_25;
        float _S4327 = length_1(mean_c_22);
        bool _S4328;
        if(_S4327 < near_plane_15)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = _S4327 > far_plane_15;
        }
        if(_S4328)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S4329 = scale_27.x;
        float sx_6 = (F32_exp((_S4329)));
        float _S4330 = scale_27.y;
        float sy_6 = (F32_exp((_S4330)));
        float sz_8 = scale_27.z - 0.5f * (_S4329 + _S4330);
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
        Matrix<float, 3, 3>  _S4331 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_28 + z2_50), 2.0f * (xy_28 + wz_28), 2.0f * (xz_28 - wy_28), 2.0f * (xy_28 - wz_28), 1.0f - 2.0f * (x2_28 + z2_50), 2.0f * (yz_28 + wx_28), 2.0f * (xz_28 + wy_28), 2.0f * (yz_28 - wx_28), 1.0f - 2.0f * (x2_28 + y2_28)));
        float3  vert0_4 = mul_0(_S4331, make_float3 (sx_6, 0.0f, 0.0f)) + mean_27;
        float3  vert1_4 = mul_0(_S4331, make_float3 (sx_6 * (-0.5f + sz_8), sy_6, 0.0f)) + mean_27;
        float3  vert2_4 = mul_0(_S4331, make_float3 (sx_6 * (-0.5f - sz_8), - sy_6, 0.0f)) + mean_27;
        float3  vert0_c_7 = mul_0(R_28, vert0_4) + t_25;
        float3  vert1_c_7 = mul_0(R_28, vert1_4) + t_25;
        float3  vert2_c_7 = mul_0(R_28, vert2_4) + t_25;
        float _S4332 = length_1(vert0_c_7);
        float _S4333 = length_1(vert1_c_7);
        float _S4334 = length_1(vert2_c_7);
        if(_S4332 < near_plane_15)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = _S4332 > far_plane_15;
        }
        if(_S4328)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = _S4333 < near_plane_15;
        }
        if(_S4328)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = _S4333 > far_plane_15;
        }
        if(_S4328)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = _S4334 < near_plane_15;
        }
        if(_S4328)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = _S4334 > far_plane_15;
        }
        if(_S4328)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  uv0_6;
        float k_10;
        for(;;)
        {
            float2  _S4335 = float2 {vert0_c_7.x, vert0_c_7.y};
            float r_30 = length_0(_S4335);
            float _S4336 = vert0_c_7.z;
            float theta_26 = (F32_atan2((r_30), (_S4336)));
            if(theta_26 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_26 * theta_26 / 3.0f) / _S4336;
            }
            else
            {
                k_10 = theta_26 / r_30;
            }
            float2  uv0_7 = _S4335 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4337 = camera_distortion_jac_0(uv0_7, dist_coeffs_35);
            bool _S4338 = !((F32_min((determinant_0(_S4337)), ((F32_min((_S4337.rows[int(0)].x), (_S4337.rows[int(1)].y)))))) > 0.0f);
            _S4324 = _S4338;
            if(_S4338)
            {
                uv0_6 = uv0_7;
                break;
            }
            float u_65 = uv0_7.x;
            float v_65 = uv0_7.y;
            float r2_65 = u_65 * u_65 + v_65 * v_65;
            float2  _S4339 = uv0_7 * make_float2 (1.0f + r2_65 * ((*dist_coeffs_35)[int(0)] + r2_65 * ((*dist_coeffs_35)[int(1)] + r2_65 * ((*dist_coeffs_35)[int(2)] + r2_65 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_65 * v_65 + (*dist_coeffs_35)[int(5)] * (r2_65 + 2.0f * u_65 * u_65) + (*dist_coeffs_35)[int(6)] * r2_65, 2.0f * (*dist_coeffs_35)[int(5)] * u_65 * v_65 + (*dist_coeffs_35)[int(4)] * (r2_65 + 2.0f * v_65 * v_65) + (*dist_coeffs_35)[int(7)] * r2_65);
            float2  _S4340 = _S4339 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4339.x + (*dist_coeffs_35)[int(9)] * _S4339.y, 0.0f);
            uv0_6 = make_float2 (fx_31 * _S4340.x + cx_26, fy_31 * _S4340.y + cy_26);
            break;
        }
        float2  uv1_6;
        bool all_valid_14 = true & (!_S4324);
        for(;;)
        {
            float2  _S4341 = float2 {vert1_c_7.x, vert1_c_7.y};
            float r_31 = length_0(_S4341);
            float _S4342 = vert1_c_7.z;
            float theta_27 = (F32_atan2((r_31), (_S4342)));
            if(theta_27 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_27 * theta_27 / 3.0f) / _S4342;
            }
            else
            {
                k_10 = theta_27 / r_31;
            }
            float2  uv1_7 = _S4341 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4343 = camera_distortion_jac_0(uv1_7, dist_coeffs_35);
            bool _S4344 = !((F32_min((determinant_0(_S4343)), ((F32_min((_S4343.rows[int(0)].x), (_S4343.rows[int(1)].y)))))) > 0.0f);
            _S4325 = _S4344;
            if(_S4344)
            {
                uv1_6 = uv1_7;
                break;
            }
            float u_66 = uv1_7.x;
            float v_66 = uv1_7.y;
            float r2_66 = u_66 * u_66 + v_66 * v_66;
            float2  _S4345 = uv1_7 * make_float2 (1.0f + r2_66 * ((*dist_coeffs_35)[int(0)] + r2_66 * ((*dist_coeffs_35)[int(1)] + r2_66 * ((*dist_coeffs_35)[int(2)] + r2_66 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_66 * v_66 + (*dist_coeffs_35)[int(5)] * (r2_66 + 2.0f * u_66 * u_66) + (*dist_coeffs_35)[int(6)] * r2_66, 2.0f * (*dist_coeffs_35)[int(5)] * u_66 * v_66 + (*dist_coeffs_35)[int(4)] * (r2_66 + 2.0f * v_66 * v_66) + (*dist_coeffs_35)[int(7)] * r2_66);
            float2  _S4346 = _S4345 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4345.x + (*dist_coeffs_35)[int(9)] * _S4345.y, 0.0f);
            uv1_6 = make_float2 (fx_31 * _S4346.x + cx_26, fy_31 * _S4346.y + cy_26);
            break;
        }
        float2  uv2_6;
        bool all_valid_15 = all_valid_14 & (!_S4325);
        for(;;)
        {
            float2  _S4347 = float2 {vert2_c_7.x, vert2_c_7.y};
            float r_32 = length_0(_S4347);
            float _S4348 = vert2_c_7.z;
            float theta_28 = (F32_atan2((r_32), (_S4348)));
            if(theta_28 < 0.00100000004749745f)
            {
                k_10 = (1.0f - theta_28 * theta_28 / 3.0f) / _S4348;
            }
            else
            {
                k_10 = theta_28 / r_32;
            }
            float2  uv2_7 = _S4347 * make_float2 (k_10);
            Matrix<float, 2, 2>  _S4349 = camera_distortion_jac_0(uv2_7, dist_coeffs_35);
            bool _S4350 = !((F32_min((determinant_0(_S4349)), ((F32_min((_S4349.rows[int(0)].x), (_S4349.rows[int(1)].y)))))) > 0.0f);
            _S4326 = _S4350;
            if(_S4350)
            {
                uv2_6 = uv2_7;
                break;
            }
            float u_67 = uv2_7.x;
            float v_67 = uv2_7.y;
            float r2_67 = u_67 * u_67 + v_67 * v_67;
            float2  _S4351 = uv2_7 * make_float2 (1.0f + r2_67 * ((*dist_coeffs_35)[int(0)] + r2_67 * ((*dist_coeffs_35)[int(1)] + r2_67 * ((*dist_coeffs_35)[int(2)] + r2_67 * (*dist_coeffs_35)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_35)[int(4)] * u_67 * v_67 + (*dist_coeffs_35)[int(5)] * (r2_67 + 2.0f * u_67 * u_67) + (*dist_coeffs_35)[int(6)] * r2_67, 2.0f * (*dist_coeffs_35)[int(5)] * u_67 * v_67 + (*dist_coeffs_35)[int(4)] * (r2_67 + 2.0f * v_67 * v_67) + (*dist_coeffs_35)[int(7)] * r2_67);
            float2  _S4352 = _S4351 + make_float2 ((*dist_coeffs_35)[int(8)] * _S4351.x + (*dist_coeffs_35)[int(9)] * _S4351.y, 0.0f);
            uv2_6 = make_float2 (fx_31 * _S4352.x + cx_26, fy_31 * _S4352.y + cy_26);
            break;
        }
        if(!(all_valid_15 & (!_S4326)))
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float2  e0_11 = uv1_6 - uv0_6;
        float2  e1_11 = uv2_6 - uv1_6;
        float offset_5 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_11.y))))) - 1.0f) * ((F32_abs((e0_11.x * e1_11.y - e0_11.y * e1_11.x))) / (length_0(e0_11) + length_0(e1_11) + length_0(uv0_6 - uv2_6)));
        float _S4353 = uv0_6.x;
        float _S4354 = uv1_6.x;
        float _S4355 = uv2_6.x;
        float xmax_8 = (F32_max(((F32_max((_S4353), (_S4354)))), (_S4355))) + offset_5;
        float xmin_8 = (F32_min(((F32_min((_S4353), (_S4354)))), (_S4355))) - offset_5;
        float _S4356 = uv0_6.y;
        float _S4357 = uv1_6.y;
        float _S4358 = uv2_6.y;
        float ymax_8 = (F32_max(((F32_max((_S4356), (_S4357)))), (_S4358))) + offset_5;
        float ymin_8 = (F32_min(((F32_min((_S4356), (_S4357)))), (_S4358))) - offset_5;
        if(xmax_8 <= 0.0f)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = xmin_8 >= float(image_width_22);
        }
        if(_S4328)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = ymax_8 <= 0.0f;
        }
        if(_S4328)
        {
            _S4328 = true;
        }
        else
        {
            _S4328 = ymin_8 >= float(image_height_22);
        }
        if(_S4328)
        {
            _S4328 = true;
        }
        else
        {
            if((mean_c_22.z) <= 0.0f)
            {
                if(xmin_8 <= 0.0f)
                {
                    _S4328 = xmax_8 >= float(image_width_22);
                }
                else
                {
                    _S4328 = false;
                }
                if(_S4328)
                {
                    _S4328 = true;
                }
                else
                {
                    if(ymin_8 <= 0.0f)
                    {
                        _S4328 = ymax_8 >= float(image_width_22);
                    }
                    else
                    {
                        _S4328 = false;
                    }
                }
            }
            else
            {
                _S4328 = false;
            }
        }
        if(_S4328)
        {
            *aabb_xyxy_15 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_15 = make_int4 (int((F32_floor((xmin_8)))), int((F32_floor((ymin_8)))), int((F32_ceil((xmax_8)))), int((F32_ceil((ymax_8)))));
        *depth_17 = (F32_log((length_1(vert0_c_7 + vert1_c_7 + vert2_c_7) / 3.0f + 9.999999960041972e-13f)));
        float3  _S4359 = mean_27 - - mul_0(transpose_0(R_28), t_25);
        float _S4360 = _S4359.x;
        float _S4361 = _S4359.y;
        float _S4362 = _S4359.z;
        float norm_15 = (F32_sqrt((_S4360 * _S4360 + _S4361 * _S4361 + _S4362 * _S4362)));
        float x_56 = _S4360 / norm_15;
        float y_25 = _S4361 / norm_15;
        float z_22 = _S4362 / norm_15;
        float z2_51 = z_22 * z_22;
        float fTmp0B_22 = -1.09254848957061768f * z_22;
        float fC1_22 = x_56 * x_56 - y_25 * y_25;
        float fS1_22 = 2.0f * x_56 * y_25;
        float fTmp0C_22 = -2.28522896766662598f * z2_51 + 0.4570457935333252f;
        float fTmp1B_22 = 1.44530570507049561f * z_22;
        float3  color_8 = make_float3 (0.282094806432724f) * (*sh_coeffs_22)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_25) * (*sh_coeffs_22)[int(1)] + make_float3 (z_22) * (*sh_coeffs_22)[int(2)] - make_float3 (x_56) * (*sh_coeffs_22)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_22) * (*sh_coeffs_22)[int(4)] + make_float3 (fTmp0B_22 * y_25) * (*sh_coeffs_22)[int(5)] + make_float3 (0.94617468118667603f * z2_51 - 0.31539157032966614f) * (*sh_coeffs_22)[int(6)] + make_float3 (fTmp0B_22 * x_56) * (*sh_coeffs_22)[int(7)] + make_float3 (0.54627424478530884f * fC1_22) * (*sh_coeffs_22)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_56 * fS1_22 + y_25 * fC1_22)) * (*sh_coeffs_22)[int(9)] + make_float3 (fTmp1B_22 * fS1_22) * (*sh_coeffs_22)[int(10)] + make_float3 (fTmp0C_22 * y_25) * (*sh_coeffs_22)[int(11)] + make_float3 (z_22 * (1.86588168144226074f * z2_51 - 1.11952900886535645f)) * (*sh_coeffs_22)[int(12)] + make_float3 (fTmp0C_22 * x_56) * (*sh_coeffs_22)[int(13)] + make_float3 (fTmp1B_22 * fC1_22) * (*sh_coeffs_22)[int(14)] + make_float3 (-0.59004360437393188f * (x_56 * fC1_22 - y_25 * fS1_22)) * (*sh_coeffs_22)[int(15)]);
        float3  _S4363 = make_float3 (0.0f);
        (*rgbs_1)[int(0)] = max_0(color_8 + (*ch_coeffs_7)[int(0)] + make_float3 (0.5f), _S4363);
        float3  _S4364 = color_8 - (*ch_coeffs_7)[int(0)] * make_float3 (0.5f);
        float3  _S4365 = (*ch_coeffs_7)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
        (*rgbs_1)[int(1)] = max_0(_S4364 + _S4365 + make_float3 (0.5f), _S4363);
        (*rgbs_1)[int(2)] = max_0(_S4364 - _S4365 + make_float3 (0.5f), _S4363);
        (*verts_1)[int(0)] = vert0_4;
        (*verts_1)[int(1)] = vert1_4;
        (*verts_1)[int(2)] = vert2_4;
        float3  _S4366 = normalize_0(cross_0(vert1_c_7 - vert0_c_7, vert2_c_7 - vert0_c_7));
        *normal_5 = _S4366 * make_float3 (float(- (F32_sign((dot_0(_S4366, mean_c_22))))));
        break;
    }
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_persp_differentiable(float3  mean_28, float4  quat_29, float3  scale_28, float2  hardness_12, FixedArray<float3 , 16>  * sh_coeffs_23, FixedArray<float3 , 2>  * ch_coeffs_8, Matrix<float, 3, 3>  R_29, float3  t_26, float fx_32, float fy_32, float cx_27, float cy_27, FixedArray<float, 10>  * dist_coeffs_36, uint image_width_23, uint image_height_23, float near_plane_16, float far_plane_16, int4  * aabb_xyxy_16, float * depth_18, FixedArray<float3 , 3>  * verts_2, FixedArray<float3 , 3>  * rgbs_2, float3  * normal_6)
{
    float3  mean_c_23 = mul_0(R_29, mean_28) + t_26;
    float _S4367 = scale_28.x;
    float sx_7 = (F32_exp((_S4367)));
    float _S4368 = scale_28.y;
    float sy_7 = (F32_exp((_S4368)));
    float sz_9 = scale_28.z - 0.5f * (_S4367 + _S4368);
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
    Matrix<float, 3, 3>  _S4369 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_29 + z2_52), 2.0f * (xy_29 + wz_29), 2.0f * (xz_29 - wy_29), 2.0f * (xy_29 - wz_29), 1.0f - 2.0f * (x2_29 + z2_52), 2.0f * (yz_29 + wx_29), 2.0f * (xz_29 + wy_29), 2.0f * (yz_29 - wx_29), 1.0f - 2.0f * (x2_29 + y2_29)));
    float3  vert0_5 = mul_0(_S4369, make_float3 (sx_7, 0.0f, 0.0f)) + mean_28;
    float3  vert1_5 = mul_0(_S4369, make_float3 (sx_7 * (-0.5f + sz_9), sy_7, 0.0f)) + mean_28;
    float3  vert2_5 = mul_0(_S4369, make_float3 (sx_7 * (-0.5f - sz_9), - sy_7, 0.0f)) + mean_28;
    float3  vert0_c_8 = mul_0(R_29, vert0_5) + t_26;
    float3  vert1_c_8 = mul_0(R_29, vert1_5) + t_26;
    float3  vert2_c_8 = mul_0(R_29, vert2_5) + t_26;
    float2  _S4370 = float2 {vert0_c_8.x, vert0_c_8.y} / make_float2 (vert0_c_8.z);
    float u_68 = _S4370.x;
    float v_68 = _S4370.y;
    float r2_68 = u_68 * u_68 + v_68 * v_68;
    float _S4371 = 2.0f * (*dist_coeffs_36)[int(4)];
    float _S4372 = 2.0f * (*dist_coeffs_36)[int(5)];
    float2  _S4373 = _S4370 * make_float2 (1.0f + r2_68 * ((*dist_coeffs_36)[int(0)] + r2_68 * ((*dist_coeffs_36)[int(1)] + r2_68 * ((*dist_coeffs_36)[int(2)] + r2_68 * (*dist_coeffs_36)[int(3)])))) + make_float2 (_S4371 * u_68 * v_68 + (*dist_coeffs_36)[int(5)] * (r2_68 + 2.0f * u_68 * u_68) + (*dist_coeffs_36)[int(6)] * r2_68, _S4372 * u_68 * v_68 + (*dist_coeffs_36)[int(4)] * (r2_68 + 2.0f * v_68 * v_68) + (*dist_coeffs_36)[int(7)] * r2_68);
    float2  _S4374 = _S4373 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4373.x + (*dist_coeffs_36)[int(9)] * _S4373.y, 0.0f);
    float _S4375 = fx_32 * _S4374.x + cx_27;
    float _S4376 = fy_32 * _S4374.y + cy_27;
    float2  uv0_8 = make_float2 (_S4375, _S4376);
    float2  _S4377 = float2 {vert1_c_8.x, vert1_c_8.y} / make_float2 (vert1_c_8.z);
    float u_69 = _S4377.x;
    float v_69 = _S4377.y;
    float r2_69 = u_69 * u_69 + v_69 * v_69;
    float2  _S4378 = _S4377 * make_float2 (1.0f + r2_69 * ((*dist_coeffs_36)[int(0)] + r2_69 * ((*dist_coeffs_36)[int(1)] + r2_69 * ((*dist_coeffs_36)[int(2)] + r2_69 * (*dist_coeffs_36)[int(3)])))) + make_float2 (_S4371 * u_69 * v_69 + (*dist_coeffs_36)[int(5)] * (r2_69 + 2.0f * u_69 * u_69) + (*dist_coeffs_36)[int(6)] * r2_69, _S4372 * u_69 * v_69 + (*dist_coeffs_36)[int(4)] * (r2_69 + 2.0f * v_69 * v_69) + (*dist_coeffs_36)[int(7)] * r2_69);
    float2  _S4379 = _S4378 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4378.x + (*dist_coeffs_36)[int(9)] * _S4378.y, 0.0f);
    float _S4380 = fx_32 * _S4379.x + cx_27;
    float _S4381 = fy_32 * _S4379.y + cy_27;
    float2  uv1_8 = make_float2 (_S4380, _S4381);
    float2  _S4382 = float2 {vert2_c_8.x, vert2_c_8.y} / make_float2 (vert2_c_8.z);
    float u_70 = _S4382.x;
    float v_70 = _S4382.y;
    float r2_70 = u_70 * u_70 + v_70 * v_70;
    float2  _S4383 = _S4382 * make_float2 (1.0f + r2_70 * ((*dist_coeffs_36)[int(0)] + r2_70 * ((*dist_coeffs_36)[int(1)] + r2_70 * ((*dist_coeffs_36)[int(2)] + r2_70 * (*dist_coeffs_36)[int(3)])))) + make_float2 (_S4371 * u_70 * v_70 + (*dist_coeffs_36)[int(5)] * (r2_70 + 2.0f * u_70 * u_70) + (*dist_coeffs_36)[int(6)] * r2_70, _S4372 * u_70 * v_70 + (*dist_coeffs_36)[int(4)] * (r2_70 + 2.0f * v_70 * v_70) + (*dist_coeffs_36)[int(7)] * r2_70);
    float2  _S4384 = _S4383 + make_float2 ((*dist_coeffs_36)[int(8)] * _S4383.x + (*dist_coeffs_36)[int(9)] * _S4383.y, 0.0f);
    float _S4385 = fx_32 * _S4384.x + cx_27;
    float _S4386 = fy_32 * _S4384.y + cy_27;
    float2  uv2_8 = make_float2 (_S4385, _S4386);
    float2  e0_12 = uv1_8 - uv0_8;
    float2  e1_12 = uv2_8 - uv1_8;
    float offset_6 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_12.y))))) - 1.0f) * ((F32_abs((e0_12.x * e1_12.y - e0_12.y * e1_12.x))) / (length_0(e0_12) + length_0(e1_12) + length_0(uv0_8 - uv2_8)));
    *aabb_xyxy_16 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4375), (_S4380)))), (_S4385))) - offset_6)))), int((F32_floor(((F32_min(((F32_min((_S4376), (_S4381)))), (_S4386))) - offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4375), (_S4380)))), (_S4385))) + offset_6)))), int((F32_ceil(((F32_max(((F32_max((_S4376), (_S4381)))), (_S4386))) + offset_6)))));
    *depth_18 = (F32_log((length_1(vert0_c_8 + vert1_c_8 + vert2_c_8) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4387 = mean_28 - - mul_0(transpose_0(R_29), t_26);
    float _S4388 = _S4387.x;
    float _S4389 = _S4387.y;
    float _S4390 = _S4387.z;
    float norm_16 = (F32_sqrt((_S4388 * _S4388 + _S4389 * _S4389 + _S4390 * _S4390)));
    float x_58 = _S4388 / norm_16;
    float y_26 = _S4389 / norm_16;
    float z_23 = _S4390 / norm_16;
    float z2_53 = z_23 * z_23;
    float fTmp0B_23 = -1.09254848957061768f * z_23;
    float fC1_23 = x_58 * x_58 - y_26 * y_26;
    float fS1_23 = 2.0f * x_58 * y_26;
    float fTmp0C_23 = -2.28522896766662598f * z2_53 + 0.4570457935333252f;
    float fTmp1B_23 = 1.44530570507049561f * z_23;
    float3  color_9 = make_float3 (0.282094806432724f) * (*sh_coeffs_23)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_26) * (*sh_coeffs_23)[int(1)] + make_float3 (z_23) * (*sh_coeffs_23)[int(2)] - make_float3 (x_58) * (*sh_coeffs_23)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_23) * (*sh_coeffs_23)[int(4)] + make_float3 (fTmp0B_23 * y_26) * (*sh_coeffs_23)[int(5)] + make_float3 (0.94617468118667603f * z2_53 - 0.31539157032966614f) * (*sh_coeffs_23)[int(6)] + make_float3 (fTmp0B_23 * x_58) * (*sh_coeffs_23)[int(7)] + make_float3 (0.54627424478530884f * fC1_23) * (*sh_coeffs_23)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_58 * fS1_23 + y_26 * fC1_23)) * (*sh_coeffs_23)[int(9)] + make_float3 (fTmp1B_23 * fS1_23) * (*sh_coeffs_23)[int(10)] + make_float3 (fTmp0C_23 * y_26) * (*sh_coeffs_23)[int(11)] + make_float3 (z_23 * (1.86588168144226074f * z2_53 - 1.11952900886535645f)) * (*sh_coeffs_23)[int(12)] + make_float3 (fTmp0C_23 * x_58) * (*sh_coeffs_23)[int(13)] + make_float3 (fTmp1B_23 * fC1_23) * (*sh_coeffs_23)[int(14)] + make_float3 (-0.59004360437393188f * (x_58 * fC1_23 - y_26 * fS1_23)) * (*sh_coeffs_23)[int(15)]);
    float3  _S4391 = make_float3 (0.0f);
    (*rgbs_2)[int(0)] = max_0(color_9 + (*ch_coeffs_8)[int(0)] + make_float3 (0.5f), _S4391);
    float3  _S4392 = color_9 - (*ch_coeffs_8)[int(0)] * make_float3 (0.5f);
    float3  _S4393 = (*ch_coeffs_8)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_2)[int(1)] = max_0(_S4392 + _S4393 + make_float3 (0.5f), _S4391);
    (*rgbs_2)[int(2)] = max_0(_S4392 - _S4393 + make_float3 (0.5f), _S4391);
    (*verts_2)[int(0)] = vert0_5;
    (*verts_2)[int(1)] = vert1_5;
    (*verts_2)[int(2)] = vert2_5;
    float3  _S4394 = normalize_0(cross_0(vert1_c_8 - vert0_c_8, vert2_c_8 - vert0_c_8));
    *normal_6 = _S4394 * make_float3 (float(- (F32_sign((dot_0(_S4394, mean_c_23))))));
    return;
}

inline __device__ void _projection_opaque_triangle_eval3d_fisheye_differentiable(float3  mean_29, float4  quat_30, float3  scale_29, float2  hardness_13, FixedArray<float3 , 16>  * sh_coeffs_24, FixedArray<float3 , 2>  * ch_coeffs_9, Matrix<float, 3, 3>  R_30, float3  t_27, float fx_33, float fy_33, float cx_28, float cy_28, FixedArray<float, 10>  * dist_coeffs_37, uint image_width_24, uint image_height_24, float near_plane_17, float far_plane_17, int4  * aabb_xyxy_17, float * depth_19, FixedArray<float3 , 3>  * verts_3, FixedArray<float3 , 3>  * rgbs_3, float3  * normal_7)
{
    float3  mean_c_24 = mul_0(R_30, mean_29) + t_27;
    float _S4395 = scale_29.x;
    float sx_8 = (F32_exp((_S4395)));
    float _S4396 = scale_29.y;
    float sy_8 = (F32_exp((_S4396)));
    float sz_10 = scale_29.z - 0.5f * (_S4395 + _S4396);
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
    Matrix<float, 3, 3>  _S4397 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_30 + z2_54), 2.0f * (xy_30 + wz_30), 2.0f * (xz_30 - wy_30), 2.0f * (xy_30 - wz_30), 1.0f - 2.0f * (x2_30 + z2_54), 2.0f * (yz_30 + wx_30), 2.0f * (xz_30 + wy_30), 2.0f * (yz_30 - wx_30), 1.0f - 2.0f * (x2_30 + y2_30)));
    float3  vert0_6 = mul_0(_S4397, make_float3 (sx_8, 0.0f, 0.0f)) + mean_29;
    float3  vert1_6 = mul_0(_S4397, make_float3 (sx_8 * (-0.5f + sz_10), sy_8, 0.0f)) + mean_29;
    float3  vert2_6 = mul_0(_S4397, make_float3 (sx_8 * (-0.5f - sz_10), - sy_8, 0.0f)) + mean_29;
    float3  vert0_c_9 = mul_0(R_30, vert0_6) + t_27;
    float3  vert1_c_9 = mul_0(R_30, vert1_6) + t_27;
    float3  vert2_c_9 = mul_0(R_30, vert2_6) + t_27;
    float2  _S4398 = float2 {vert0_c_9.x, vert0_c_9.y};
    float r_33 = length_0(_S4398);
    float _S4399 = vert0_c_9.z;
    float theta_29 = (F32_atan2((r_33), (_S4399)));
    float k_11;
    if(theta_29 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_29 * theta_29 / 3.0f) / _S4399;
    }
    else
    {
        k_11 = theta_29 / r_33;
    }
    float2  _S4400 = _S4398 * make_float2 (k_11);
    float u_71 = _S4400.x;
    float v_71 = _S4400.y;
    float r2_71 = u_71 * u_71 + v_71 * v_71;
    float _S4401 = 2.0f * (*dist_coeffs_37)[int(4)];
    float _S4402 = 2.0f * (*dist_coeffs_37)[int(5)];
    float2  _S4403 = _S4400 * make_float2 (1.0f + r2_71 * ((*dist_coeffs_37)[int(0)] + r2_71 * ((*dist_coeffs_37)[int(1)] + r2_71 * ((*dist_coeffs_37)[int(2)] + r2_71 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4401 * u_71 * v_71 + (*dist_coeffs_37)[int(5)] * (r2_71 + 2.0f * u_71 * u_71) + (*dist_coeffs_37)[int(6)] * r2_71, _S4402 * u_71 * v_71 + (*dist_coeffs_37)[int(4)] * (r2_71 + 2.0f * v_71 * v_71) + (*dist_coeffs_37)[int(7)] * r2_71);
    float2  _S4404 = _S4403 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4403.x + (*dist_coeffs_37)[int(9)] * _S4403.y, 0.0f);
    float _S4405 = fx_33 * _S4404.x + cx_28;
    float _S4406 = fy_33 * _S4404.y + cy_28;
    float2  uv0_9 = make_float2 (_S4405, _S4406);
    float2  _S4407 = float2 {vert1_c_9.x, vert1_c_9.y};
    float r_34 = length_0(_S4407);
    float _S4408 = vert1_c_9.z;
    float theta_30 = (F32_atan2((r_34), (_S4408)));
    if(theta_30 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_30 * theta_30 / 3.0f) / _S4408;
    }
    else
    {
        k_11 = theta_30 / r_34;
    }
    float2  _S4409 = _S4407 * make_float2 (k_11);
    float u_72 = _S4409.x;
    float v_72 = _S4409.y;
    float r2_72 = u_72 * u_72 + v_72 * v_72;
    float2  _S4410 = _S4409 * make_float2 (1.0f + r2_72 * ((*dist_coeffs_37)[int(0)] + r2_72 * ((*dist_coeffs_37)[int(1)] + r2_72 * ((*dist_coeffs_37)[int(2)] + r2_72 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4401 * u_72 * v_72 + (*dist_coeffs_37)[int(5)] * (r2_72 + 2.0f * u_72 * u_72) + (*dist_coeffs_37)[int(6)] * r2_72, _S4402 * u_72 * v_72 + (*dist_coeffs_37)[int(4)] * (r2_72 + 2.0f * v_72 * v_72) + (*dist_coeffs_37)[int(7)] * r2_72);
    float2  _S4411 = _S4410 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4410.x + (*dist_coeffs_37)[int(9)] * _S4410.y, 0.0f);
    float _S4412 = fx_33 * _S4411.x + cx_28;
    float _S4413 = fy_33 * _S4411.y + cy_28;
    float2  uv1_9 = make_float2 (_S4412, _S4413);
    float2  _S4414 = float2 {vert2_c_9.x, vert2_c_9.y};
    float r_35 = length_0(_S4414);
    float _S4415 = vert2_c_9.z;
    float theta_31 = (F32_atan2((r_35), (_S4415)));
    if(theta_31 < 0.00100000004749745f)
    {
        k_11 = (1.0f - theta_31 * theta_31 / 3.0f) / _S4415;
    }
    else
    {
        k_11 = theta_31 / r_35;
    }
    float2  _S4416 = _S4414 * make_float2 (k_11);
    float u_73 = _S4416.x;
    float v_73 = _S4416.y;
    float r2_73 = u_73 * u_73 + v_73 * v_73;
    float2  _S4417 = _S4416 * make_float2 (1.0f + r2_73 * ((*dist_coeffs_37)[int(0)] + r2_73 * ((*dist_coeffs_37)[int(1)] + r2_73 * ((*dist_coeffs_37)[int(2)] + r2_73 * (*dist_coeffs_37)[int(3)])))) + make_float2 (_S4401 * u_73 * v_73 + (*dist_coeffs_37)[int(5)] * (r2_73 + 2.0f * u_73 * u_73) + (*dist_coeffs_37)[int(6)] * r2_73, _S4402 * u_73 * v_73 + (*dist_coeffs_37)[int(4)] * (r2_73 + 2.0f * v_73 * v_73) + (*dist_coeffs_37)[int(7)] * r2_73);
    float2  _S4418 = _S4417 + make_float2 ((*dist_coeffs_37)[int(8)] * _S4417.x + (*dist_coeffs_37)[int(9)] * _S4417.y, 0.0f);
    float _S4419 = fx_33 * _S4418.x + cx_28;
    float _S4420 = fy_33 * _S4418.y + cy_28;
    float2  uv2_9 = make_float2 (_S4419, _S4420);
    float2  e0_13 = uv1_9 - uv0_9;
    float2  e1_13 = uv2_9 - uv1_9;
    float offset_7 = (1.0f / (1.0f - (F32_exp2((-1.0f / (1.0f - hardness_13.y))))) - 1.0f) * ((F32_abs((e0_13.x * e1_13.y - e0_13.y * e1_13.x))) / (length_0(e0_13) + length_0(e1_13) + length_0(uv0_9 - uv2_9)));
    *aabb_xyxy_17 = make_int4 (int((F32_floor(((F32_min(((F32_min((_S4405), (_S4412)))), (_S4419))) - offset_7)))), int((F32_floor(((F32_min(((F32_min((_S4406), (_S4413)))), (_S4420))) - offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4405), (_S4412)))), (_S4419))) + offset_7)))), int((F32_ceil(((F32_max(((F32_max((_S4406), (_S4413)))), (_S4420))) + offset_7)))));
    *depth_19 = (F32_log((length_1(vert0_c_9 + vert1_c_9 + vert2_c_9) / 3.0f + 9.999999960041972e-13f)));
    float3  _S4421 = mean_29 - - mul_0(transpose_0(R_30), t_27);
    float _S4422 = _S4421.x;
    float _S4423 = _S4421.y;
    float _S4424 = _S4421.z;
    float norm_17 = (F32_sqrt((_S4422 * _S4422 + _S4423 * _S4423 + _S4424 * _S4424)));
    float x_60 = _S4422 / norm_17;
    float y_27 = _S4423 / norm_17;
    float z_24 = _S4424 / norm_17;
    float z2_55 = z_24 * z_24;
    float fTmp0B_24 = -1.09254848957061768f * z_24;
    float fC1_24 = x_60 * x_60 - y_27 * y_27;
    float fS1_24 = 2.0f * x_60 * y_27;
    float fTmp0C_24 = -2.28522896766662598f * z2_55 + 0.4570457935333252f;
    float fTmp1B_24 = 1.44530570507049561f * z_24;
    float3  color_10 = make_float3 (0.282094806432724f) * (*sh_coeffs_24)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (- y_27) * (*sh_coeffs_24)[int(1)] + make_float3 (z_24) * (*sh_coeffs_24)[int(2)] - make_float3 (x_60) * (*sh_coeffs_24)[int(3)]) + (make_float3 (0.54627424478530884f * fS1_24) * (*sh_coeffs_24)[int(4)] + make_float3 (fTmp0B_24 * y_27) * (*sh_coeffs_24)[int(5)] + make_float3 (0.94617468118667603f * z2_55 - 0.31539157032966614f) * (*sh_coeffs_24)[int(6)] + make_float3 (fTmp0B_24 * x_60) * (*sh_coeffs_24)[int(7)] + make_float3 (0.54627424478530884f * fC1_24) * (*sh_coeffs_24)[int(8)]) + (make_float3 (-0.59004360437393188f * (x_60 * fS1_24 + y_27 * fC1_24)) * (*sh_coeffs_24)[int(9)] + make_float3 (fTmp1B_24 * fS1_24) * (*sh_coeffs_24)[int(10)] + make_float3 (fTmp0C_24 * y_27) * (*sh_coeffs_24)[int(11)] + make_float3 (z_24 * (1.86588168144226074f * z2_55 - 1.11952900886535645f)) * (*sh_coeffs_24)[int(12)] + make_float3 (fTmp0C_24 * x_60) * (*sh_coeffs_24)[int(13)] + make_float3 (fTmp1B_24 * fC1_24) * (*sh_coeffs_24)[int(14)] + make_float3 (-0.59004360437393188f * (x_60 * fC1_24 - y_27 * fS1_24)) * (*sh_coeffs_24)[int(15)]);
    float3  _S4425 = make_float3 (0.0f);
    (*rgbs_3)[int(0)] = max_0(color_10 + (*ch_coeffs_9)[int(0)] + make_float3 (0.5f), _S4425);
    float3  _S4426 = color_10 - (*ch_coeffs_9)[int(0)] * make_float3 (0.5f);
    float3  _S4427 = (*ch_coeffs_9)[int(1)] * make_float3 ((F32_sqrt((0.75f))));
    (*rgbs_3)[int(1)] = max_0(_S4426 + _S4427 + make_float3 (0.5f), _S4425);
    (*rgbs_3)[int(2)] = max_0(_S4426 - _S4427 + make_float3 (0.5f), _S4425);
    (*verts_3)[int(0)] = vert0_6;
    (*verts_3)[int(1)] = vert1_6;
    (*verts_3)[int(2)] = vert2_6;
    float3  _S4428 = normalize_0(cross_0(vert1_c_9 - vert0_c_9, vert2_c_9 - vert0_c_9));
    *normal_7 = _S4428 * make_float3 (float(- (F32_sign((dot_0(_S4428, mean_c_24))))));
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_persp_vjp(float3  mean_30, float4  quat_31, float3  scale_30, float2  hardness_14, FixedArray<float3 , 16>  * sh_coeffs_25, FixedArray<float3 , 2>  * ch_coeffs_10, Matrix<float, 3, 3>  R_31, float3  t_28, float fx_34, float fy_34, float cx_29, float cy_29, FixedArray<float, 10>  * dist_coeffs_38, uint image_width_25, uint image_height_25, float v_depth_9, FixedArray<float3 , 3>  * v_verts_0, FixedArray<float3 , 3>  * v_rgbs_0, float3  v_normal_2, float3  * v_mean_9, float4  * v_quat_8, float3  * v_scale_8, float2  * v_hardness_4, FixedArray<float3 , 16>  * v_sh_coeffs_7, FixedArray<float3 , 2>  * v_ch_coeffs_2, Matrix<float, 3, 3>  * v_R_9, float3  * v_t_8)
{
    float3  mean_c_25 = s_primal_ctx_mul_0(R_31, mean_30) + t_28;
    float _S4429 = scale_30.x;
    float _S4430 = s_primal_ctx_exp_1(_S4429);
    float _S4431 = scale_30.y;
    float _S4432 = s_primal_ctx_exp_1(_S4431);
    float sz_11 = scale_30.z - 0.5f * (_S4429 + _S4431);
    float _S4433 = quat_31.y;
    float x2_31 = _S4433 * _S4433;
    float y2_31 = quat_31.z * quat_31.z;
    float z2_56 = quat_31.w * quat_31.w;
    float xy_31 = quat_31.y * quat_31.z;
    float xz_31 = quat_31.y * quat_31.w;
    float yz_31 = quat_31.z * quat_31.w;
    float wx_31 = quat_31.x * quat_31.y;
    float wy_31 = quat_31.x * quat_31.z;
    float wz_31 = quat_31.x * quat_31.w;
    Matrix<float, 3, 3>  _S4434 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_31 + z2_56), 2.0f * (xy_31 + wz_31), 2.0f * (xz_31 - wy_31), 2.0f * (xy_31 - wz_31), 1.0f - 2.0f * (x2_31 + z2_56), 2.0f * (yz_31 + wx_31), 2.0f * (xz_31 + wy_31), 2.0f * (yz_31 - wx_31), 1.0f - 2.0f * (x2_31 + y2_31)));
    float3  _S4435 = make_float3 (_S4430, 0.0f, 0.0f);
    float3  vert0_7 = s_primal_ctx_mul_0(_S4434, _S4435) + mean_30;
    float _S4436 = -0.5f + sz_11;
    float3  _S4437 = make_float3 (_S4430 * _S4436, _S4432, 0.0f);
    float3  vert1_7 = s_primal_ctx_mul_0(_S4434, _S4437) + mean_30;
    float _S4438 = -0.5f - sz_11;
    float3  _S4439 = make_float3 (_S4430 * _S4438, - _S4432, 0.0f);
    float3  vert2_7 = s_primal_ctx_mul_0(_S4434, _S4439) + mean_30;
    float3  vert0_c_10 = s_primal_ctx_mul_0(R_31, vert0_7) + t_28;
    float3  vert1_c_10 = s_primal_ctx_mul_0(R_31, vert1_7) + t_28;
    float3  vert2_c_10 = s_primal_ctx_mul_0(R_31, vert2_7) + t_28;
    float2  _S4440 = float2 {vert0_c_10.x, vert0_c_10.y};
    float _S4441 = vert0_c_10.z;
    float2  _S4442 = make_float2 (_S4441);
    float2  _S4443 = _S4440 / make_float2 (_S4441);
    float2  _S4444 = make_float2 (_S4441 * _S4441);
    float u_74 = _S4443.x;
    float v_74 = _S4443.y;
    float r2_74 = u_74 * u_74 + v_74 * v_74;
    float _S4445 = (*dist_coeffs_38)[int(2)] + r2_74 * (*dist_coeffs_38)[int(3)];
    float _S4446 = (*dist_coeffs_38)[int(1)] + r2_74 * _S4445;
    float _S4447 = (*dist_coeffs_38)[int(0)] + r2_74 * _S4446;
    float radial_7 = 1.0f + r2_74 * _S4447;
    float2  _S4448 = make_float2 (radial_7);
    float _S4449 = 2.0f * (*dist_coeffs_38)[int(4)];
    float _S4450 = _S4449 * u_74;
    float _S4451 = 2.0f * u_74;
    float _S4452 = 2.0f * (*dist_coeffs_38)[int(5)];
    float _S4453 = _S4452 * u_74;
    float _S4454 = 2.0f * v_74;
    float2  _S4455 = _S4443 * make_float2 (radial_7) + make_float2 (_S4450 * v_74 + (*dist_coeffs_38)[int(5)] * (r2_74 + _S4451 * u_74) + (*dist_coeffs_38)[int(6)] * r2_74, _S4453 * v_74 + (*dist_coeffs_38)[int(4)] * (r2_74 + _S4454 * v_74) + (*dist_coeffs_38)[int(7)] * r2_74);
    float2  _S4456 = _S4455 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4455.x + (*dist_coeffs_38)[int(9)] * _S4455.y, 0.0f);
    float _S4457 = fx_34 * _S4456.x + cx_29;
    float _S4458 = fy_34 * _S4456.y + cy_29;
    float2  uv0_10 = make_float2 (_S4457, _S4458);
    float2  _S4459 = float2 {vert1_c_10.x, vert1_c_10.y};
    float _S4460 = vert1_c_10.z;
    float2  _S4461 = make_float2 (_S4460);
    float2  _S4462 = _S4459 / make_float2 (_S4460);
    float2  _S4463 = make_float2 (_S4460 * _S4460);
    float u_75 = _S4462.x;
    float v_75 = _S4462.y;
    float r2_75 = u_75 * u_75 + v_75 * v_75;
    float _S4464 = (*dist_coeffs_38)[int(2)] + r2_75 * (*dist_coeffs_38)[int(3)];
    float _S4465 = (*dist_coeffs_38)[int(1)] + r2_75 * _S4464;
    float _S4466 = (*dist_coeffs_38)[int(0)] + r2_75 * _S4465;
    float radial_8 = 1.0f + r2_75 * _S4466;
    float2  _S4467 = make_float2 (radial_8);
    float _S4468 = _S4449 * u_75;
    float _S4469 = 2.0f * u_75;
    float _S4470 = _S4452 * u_75;
    float _S4471 = 2.0f * v_75;
    float2  _S4472 = _S4462 * make_float2 (radial_8) + make_float2 (_S4468 * v_75 + (*dist_coeffs_38)[int(5)] * (r2_75 + _S4469 * u_75) + (*dist_coeffs_38)[int(6)] * r2_75, _S4470 * v_75 + (*dist_coeffs_38)[int(4)] * (r2_75 + _S4471 * v_75) + (*dist_coeffs_38)[int(7)] * r2_75);
    float2  _S4473 = _S4472 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4472.x + (*dist_coeffs_38)[int(9)] * _S4472.y, 0.0f);
    float _S4474 = fx_34 * _S4473.x + cx_29;
    float _S4475 = fy_34 * _S4473.y + cy_29;
    float2  uv1_10 = make_float2 (_S4474, _S4475);
    float2  _S4476 = float2 {vert2_c_10.x, vert2_c_10.y};
    float _S4477 = vert2_c_10.z;
    float2  _S4478 = make_float2 (_S4477);
    float2  _S4479 = _S4476 / make_float2 (_S4477);
    float2  _S4480 = make_float2 (_S4477 * _S4477);
    float u_76 = _S4479.x;
    float v_76 = _S4479.y;
    float r2_76 = u_76 * u_76 + v_76 * v_76;
    float _S4481 = (*dist_coeffs_38)[int(2)] + r2_76 * (*dist_coeffs_38)[int(3)];
    float _S4482 = (*dist_coeffs_38)[int(1)] + r2_76 * _S4481;
    float _S4483 = (*dist_coeffs_38)[int(0)] + r2_76 * _S4482;
    float radial_9 = 1.0f + r2_76 * _S4483;
    float2  _S4484 = make_float2 (radial_9);
    float _S4485 = _S4449 * u_76;
    float _S4486 = 2.0f * u_76;
    float _S4487 = _S4452 * u_76;
    float _S4488 = 2.0f * v_76;
    float2  _S4489 = _S4479 * make_float2 (radial_9) + make_float2 (_S4485 * v_76 + (*dist_coeffs_38)[int(5)] * (r2_76 + _S4486 * u_76) + (*dist_coeffs_38)[int(6)] * r2_76, _S4487 * v_76 + (*dist_coeffs_38)[int(4)] * (r2_76 + _S4488 * v_76) + (*dist_coeffs_38)[int(7)] * r2_76);
    float2  _S4490 = _S4489 + make_float2 ((*dist_coeffs_38)[int(8)] * _S4489.x + (*dist_coeffs_38)[int(9)] * _S4489.y, 0.0f);
    float _S4491 = fx_34 * _S4490.x + cx_29;
    float _S4492 = fy_34 * _S4490.y + cy_29;
    float2  uv2_10 = make_float2 (_S4491, _S4492);
    float2  e0_14 = uv1_10 - uv0_10;
    float2  e1_14 = uv2_10 - uv1_10;
    float2  e2_6 = uv0_10 - uv2_10;
    float _S4493 = e0_14.x;
    float _S4494 = e1_14.y;
    float _S4495 = e0_14.y;
    float _S4496 = e1_14.x;
    float _S4497 = _S4493 * _S4494 - _S4495 * _S4496;
    float _S4498 = 1.0f - hardness_14.y;
    float _S4499 = -1.0f / _S4498;
    float _S4500 = _S4498 * _S4498;
    float _S4501 = s_primal_ctx_max_0(_S4457, _S4474);
    float _S4502 = s_primal_ctx_min_0(_S4457, _S4474);
    float _S4503 = s_primal_ctx_max_0(_S4458, _S4475);
    float _S4504 = s_primal_ctx_min_0(_S4458, _S4475);
    float3  _S4505 = vert0_c_10 + vert1_c_10 + vert2_c_10;
    float _S4506 = length_1(_S4505) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4507 = transpose_0(R_31);
    float3  _S4508 = mean_30 - - s_primal_ctx_mul_0(_S4507, t_28);
    float _S4509 = _S4508.x;
    float _S4510 = _S4508.y;
    float _S4511 = _S4508.z;
    float _S4512 = _S4509 * _S4509 + _S4510 * _S4510 + _S4511 * _S4511;
    float _S4513 = s_primal_ctx_sqrt_0(_S4512);
    float x_61 = _S4509 / _S4513;
    float3  _S4514 = make_float3 (x_61);
    float _S4515 = _S4513 * _S4513;
    float y_28 = _S4510 / _S4513;
    float z_25 = _S4511 / _S4513;
    float3  _S4516 = make_float3 (z_25);
    float _S4517 = - y_28;
    float3  _S4518 = make_float3 (_S4517);
    float z2_57 = z_25 * z_25;
    float fTmp0B_25 = -1.09254848957061768f * z_25;
    float fC1_25 = x_61 * x_61 - y_28 * y_28;
    float _S4519 = 2.0f * x_61;
    float fS1_25 = _S4519 * y_28;
    float pSH6_7 = 0.94617468118667603f * z2_57 - 0.31539157032966614f;
    float3  _S4520 = make_float3 (pSH6_7);
    float pSH7_7 = fTmp0B_25 * x_61;
    float3  _S4521 = make_float3 (pSH7_7);
    float pSH5_7 = fTmp0B_25 * y_28;
    float3  _S4522 = make_float3 (pSH5_7);
    float pSH8_7 = 0.54627424478530884f * fC1_25;
    float3  _S4523 = make_float3 (pSH8_7);
    float pSH4_7 = 0.54627424478530884f * fS1_25;
    float3  _S4524 = make_float3 (pSH4_7);
    float fTmp0C_25 = -2.28522896766662598f * z2_57 + 0.4570457935333252f;
    float fTmp1B_25 = 1.44530570507049561f * z_25;
    float _S4525 = 1.86588168144226074f * z2_57 - 1.11952900886535645f;
    float pSH12_7 = z_25 * _S4525;
    float3  _S4526 = make_float3 (pSH12_7);
    float pSH13_7 = fTmp0C_25 * x_61;
    float3  _S4527 = make_float3 (pSH13_7);
    float pSH11_7 = fTmp0C_25 * y_28;
    float3  _S4528 = make_float3 (pSH11_7);
    float pSH14_7 = fTmp1B_25 * fC1_25;
    float3  _S4529 = make_float3 (pSH14_7);
    float pSH10_7 = fTmp1B_25 * fS1_25;
    float3  _S4530 = make_float3 (pSH10_7);
    float pSH15_7 = -0.59004360437393188f * (x_61 * fC1_25 - y_28 * fS1_25);
    float3  _S4531 = make_float3 (pSH15_7);
    float pSH9_7 = -0.59004360437393188f * (x_61 * fS1_25 + y_28 * fC1_25);
    float3  _S4532 = make_float3 (pSH9_7);
    float3  color_11 = make_float3 (0.282094806432724f) * (*sh_coeffs_25)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4517) * (*sh_coeffs_25)[int(1)] + make_float3 (z_25) * (*sh_coeffs_25)[int(2)] - make_float3 (x_61) * (*sh_coeffs_25)[int(3)]) + (make_float3 (pSH4_7) * (*sh_coeffs_25)[int(4)] + make_float3 (pSH5_7) * (*sh_coeffs_25)[int(5)] + make_float3 (pSH6_7) * (*sh_coeffs_25)[int(6)] + make_float3 (pSH7_7) * (*sh_coeffs_25)[int(7)] + make_float3 (pSH8_7) * (*sh_coeffs_25)[int(8)]) + (make_float3 (pSH9_7) * (*sh_coeffs_25)[int(9)] + make_float3 (pSH10_7) * (*sh_coeffs_25)[int(10)] + make_float3 (pSH11_7) * (*sh_coeffs_25)[int(11)] + make_float3 (pSH12_7) * (*sh_coeffs_25)[int(12)] + make_float3 (pSH13_7) * (*sh_coeffs_25)[int(13)] + make_float3 (pSH14_7) * (*sh_coeffs_25)[int(14)] + make_float3 (pSH15_7) * (*sh_coeffs_25)[int(15)]);
    float3  _S4533 = color_11 + (*ch_coeffs_10)[int(0)] + make_float3 (0.5f);
    float3  _S4534 = make_float3 (0.0f);
    float3  _S4535 = color_11 - (*ch_coeffs_10)[int(0)] * make_float3 (0.5f);
    float _S4536 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4537 = make_float3 (_S4536);
    float3  _S4538 = (*ch_coeffs_10)[int(1)] * make_float3 (_S4536);
    float3  _S4539 = _S4535 + _S4538 + make_float3 (0.5f);
    float3  _S4540 = _S4535 - _S4538 + make_float3 (0.5f);
    float3  _S4541 = vert1_c_10 - vert0_c_10;
    float3  _S4542 = vert2_c_10 - vert0_c_10;
    float3  _S4543 = s_primal_ctx_cross_0(_S4541, _S4542);
    float3  _S4544 = normalize_0(_S4543);
    float3  _S4545 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4544, mean_c_25)))))) * v_normal_2;
    float3  _S4546 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4547;
    (&_S4547)->primal_0 = _S4544;
    (&_S4547)->differential_0 = _S4546;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4548;
    (&_S4548)->primal_0 = mean_c_25;
    (&_S4548)->differential_0 = _S4546;
    s_bwd_prop_dot_0(&_S4547, &_S4548, 0.0f);
    float3  _S4549 = _S4545 + _S4547.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4550;
    (&_S4550)->primal_0 = _S4543;
    (&_S4550)->differential_0 = _S4546;
    s_bwd_normalize_impl_0(&_S4550, _S4549);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4551;
    (&_S4551)->primal_0 = _S4541;
    (&_S4551)->differential_0 = _S4546;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4552;
    (&_S4552)->primal_0 = _S4542;
    (&_S4552)->differential_0 = _S4546;
    s_bwd_prop_cross_0(&_S4551, &_S4552, _S4550.differential_0);
    float3  _S4553 = - _S4552.differential_0;
    float3  _S4554 = - _S4551.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4555;
    (&_S4555)->primal_0 = _S4540;
    (&_S4555)->differential_0 = _S4546;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4556;
    (&_S4556)->primal_0 = _S4534;
    (&_S4556)->differential_0 = _S4546;
    s_bwd_prop_max_0(&_S4555, &_S4556, (*v_rgbs_0)[int(2)]);
    float3  _S4557 = - _S4555.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4558;
    (&_S4558)->primal_0 = _S4539;
    (&_S4558)->differential_0 = _S4546;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4559;
    (&_S4559)->primal_0 = _S4534;
    (&_S4559)->differential_0 = _S4546;
    s_bwd_prop_max_0(&_S4558, &_S4559, (*v_rgbs_0)[int(1)]);
    float3  _S4560 = _S4537 * (_S4557 + _S4558.differential_0);
    float3  _S4561 = _S4555.differential_0 + _S4558.differential_0;
    float3  _S4562 = make_float3 (0.5f) * - _S4561;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4563;
    (&_S4563)->primal_0 = _S4533;
    (&_S4563)->differential_0 = _S4546;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4564;
    (&_S4564)->primal_0 = _S4534;
    (&_S4564)->differential_0 = _S4546;
    s_bwd_prop_max_0(&_S4563, &_S4564, (*v_rgbs_0)[int(0)]);
    float3  _S4565 = _S4562 + _S4563.differential_0;
    float3  _S4566 = _S4561 + _S4563.differential_0;
    float3  _S4567 = _S4531 * _S4566;
    float3  _S4568 = (*sh_coeffs_25)[int(15)] * _S4566;
    float3  _S4569 = _S4529 * _S4566;
    float3  _S4570 = (*sh_coeffs_25)[int(14)] * _S4566;
    float3  _S4571 = _S4527 * _S4566;
    float3  _S4572 = (*sh_coeffs_25)[int(13)] * _S4566;
    float3  _S4573 = _S4526 * _S4566;
    float3  _S4574 = (*sh_coeffs_25)[int(12)] * _S4566;
    float3  _S4575 = _S4528 * _S4566;
    float3  _S4576 = (*sh_coeffs_25)[int(11)] * _S4566;
    float3  _S4577 = _S4530 * _S4566;
    float3  _S4578 = (*sh_coeffs_25)[int(10)] * _S4566;
    float3  _S4579 = _S4532 * _S4566;
    float3  _S4580 = (*sh_coeffs_25)[int(9)] * _S4566;
    float s_diff_fS2_T_7 = -0.59004360437393188f * (_S4580.x + _S4580.y + _S4580.z);
    float s_diff_fC2_T_7 = -0.59004360437393188f * (_S4568.x + _S4568.y + _S4568.z);
    float _S4581 = _S4578.x + _S4578.y + _S4578.z;
    float _S4582 = _S4570.x + _S4570.y + _S4570.z;
    float _S4583 = _S4576.x + _S4576.y + _S4576.z;
    float _S4584 = _S4572.x + _S4572.y + _S4572.z;
    float _S4585 = _S4574.x + _S4574.y + _S4574.z;
    float _S4586 = - s_diff_fC2_T_7;
    float3  _S4587 = _S4523 * _S4566;
    float3  _S4588 = (*sh_coeffs_25)[int(8)] * _S4566;
    float3  _S4589 = _S4521 * _S4566;
    float3  _S4590 = (*sh_coeffs_25)[int(7)] * _S4566;
    float3  _S4591 = _S4520 * _S4566;
    float3  _S4592 = (*sh_coeffs_25)[int(6)] * _S4566;
    float3  _S4593 = _S4522 * _S4566;
    float3  _S4594 = (*sh_coeffs_25)[int(5)] * _S4566;
    float3  _S4595 = _S4524 * _S4566;
    float3  _S4596 = (*sh_coeffs_25)[int(4)] * _S4566;
    float _S4597 = _S4594.x + _S4594.y + _S4594.z;
    float _S4598 = _S4590.x + _S4590.y + _S4590.z;
    float _S4599 = fTmp1B_25 * _S4581 + x_61 * s_diff_fS2_T_7 + y_28 * _S4586 + 0.54627424478530884f * (_S4596.x + _S4596.y + _S4596.z);
    float _S4600 = fTmp1B_25 * _S4582 + y_28 * s_diff_fS2_T_7 + x_61 * s_diff_fC2_T_7 + 0.54627424478530884f * (_S4588.x + _S4588.y + _S4588.z);
    float _S4601 = y_28 * - _S4600;
    float _S4602 = x_61 * _S4600;
    float _S4603 = z_25 * (1.86588168144226074f * (z_25 * _S4585) + -2.28522896766662598f * (y_28 * _S4583 + x_61 * _S4584) + 0.94617468118667603f * (_S4592.x + _S4592.y + _S4592.z));
    float3  _S4604 = make_float3 (0.48860251903533936f) * _S4566;
    float3  _S4605 = - _S4604;
    float3  _S4606 = _S4514 * _S4605;
    float3  _S4607 = (*sh_coeffs_25)[int(3)] * _S4605;
    float3  _S4608 = _S4516 * _S4604;
    float3  _S4609 = (*sh_coeffs_25)[int(2)] * _S4604;
    float3  _S4610 = _S4518 * _S4604;
    float3  _S4611 = (*sh_coeffs_25)[int(1)] * _S4604;
    float _S4612 = (_S4525 * _S4585 + 1.44530570507049561f * (fS1_25 * _S4581 + fC1_25 * _S4582) + -1.09254848957061768f * (y_28 * _S4597 + x_61 * _S4598) + _S4603 + _S4603 + _S4609.x + _S4609.y + _S4609.z) / _S4515;
    float _S4613 = _S4513 * _S4612;
    float _S4614 = (fTmp0C_25 * _S4583 + fC1_25 * s_diff_fS2_T_7 + fS1_25 * _S4586 + fTmp0B_25 * _S4597 + _S4519 * _S4599 + _S4601 + _S4601 + - (_S4611.x + _S4611.y + _S4611.z)) / _S4515;
    float _S4615 = _S4513 * _S4614;
    float _S4616 = (fTmp0C_25 * _S4584 + fS1_25 * s_diff_fS2_T_7 + fC1_25 * s_diff_fC2_T_7 + fTmp0B_25 * _S4598 + 2.0f * (y_28 * _S4599) + _S4602 + _S4602 + _S4607.x + _S4607.y + _S4607.z) / _S4515;
    float _S4617 = _S4513 * _S4616;
    float _S4618 = _S4511 * - _S4612 + _S4510 * - _S4614 + _S4509 * - _S4616;
    DiffPair_float_0 _S4619;
    (&_S4619)->primal_0 = _S4512;
    (&_S4619)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4619, _S4618);
    float _S4620 = _S4511 * _S4619.differential_0;
    float _S4621 = _S4510 * _S4619.differential_0;
    float _S4622 = _S4509 * _S4619.differential_0;
    float3  _S4623 = make_float3 (0.282094806432724f) * _S4566;
    float3  _S4624 = make_float3 (_S4617 + _S4622 + _S4622, _S4615 + _S4621 + _S4621, _S4613 + _S4620 + _S4620);
    float3  _S4625 = - - _S4624;
    Matrix<float, 3, 3>  _S4626 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4627;
    (&_S4627)->primal_0 = _S4507;
    (&_S4627)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4628;
    (&_S4628)->primal_0 = t_28;
    (&_S4628)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4627, &_S4628, _S4625);
    Matrix<float, 3, 3>  _S4629 = transpose_0(_S4627.differential_0);
    DiffPair_float_0 _S4630;
    (&_S4630)->primal_0 = _S4506;
    (&_S4630)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4630, v_depth_9);
    float _S4631 = 0.3333333432674408f * _S4630.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4632;
    (&_S4632)->primal_0 = _S4505;
    (&_S4632)->differential_0 = _S4546;
    s_bwd_length_impl_1(&_S4632, _S4631);
    DiffPair_float_0 _S4633;
    (&_S4633)->primal_0 = _S4504;
    (&_S4633)->differential_0 = 0.0f;
    DiffPair_float_0 _S4634;
    (&_S4634)->primal_0 = _S4492;
    (&_S4634)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4633, &_S4634, 0.0f);
    DiffPair_float_0 _S4635;
    (&_S4635)->primal_0 = _S4458;
    (&_S4635)->differential_0 = 0.0f;
    DiffPair_float_0 _S4636;
    (&_S4636)->primal_0 = _S4475;
    (&_S4636)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4635, &_S4636, _S4633.differential_0);
    DiffPair_float_0 _S4637;
    (&_S4637)->primal_0 = _S4503;
    (&_S4637)->differential_0 = 0.0f;
    DiffPair_float_0 _S4638;
    (&_S4638)->primal_0 = _S4492;
    (&_S4638)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4637, &_S4638, 0.0f);
    DiffPair_float_0 _S4639;
    (&_S4639)->primal_0 = _S4458;
    (&_S4639)->differential_0 = 0.0f;
    DiffPair_float_0 _S4640;
    (&_S4640)->primal_0 = _S4475;
    (&_S4640)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4639, &_S4640, _S4637.differential_0);
    DiffPair_float_0 _S4641;
    (&_S4641)->primal_0 = _S4502;
    (&_S4641)->differential_0 = 0.0f;
    DiffPair_float_0 _S4642;
    (&_S4642)->primal_0 = _S4491;
    (&_S4642)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4641, &_S4642, 0.0f);
    DiffPair_float_0 _S4643;
    (&_S4643)->primal_0 = _S4457;
    (&_S4643)->differential_0 = 0.0f;
    DiffPair_float_0 _S4644;
    (&_S4644)->primal_0 = _S4474;
    (&_S4644)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4643, &_S4644, _S4641.differential_0);
    DiffPair_float_0 _S4645;
    (&_S4645)->primal_0 = _S4501;
    (&_S4645)->differential_0 = 0.0f;
    DiffPair_float_0 _S4646;
    (&_S4646)->primal_0 = _S4491;
    (&_S4646)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4645, &_S4646, 0.0f);
    DiffPair_float_0 _S4647;
    (&_S4647)->primal_0 = _S4457;
    (&_S4647)->differential_0 = 0.0f;
    DiffPair_float_0 _S4648;
    (&_S4648)->primal_0 = _S4474;
    (&_S4648)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S4647, &_S4648, _S4645.differential_0);
    DiffPair_float_0 _S4649;
    (&_S4649)->primal_0 = _S4499;
    (&_S4649)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S4649, 0.0f);
    float _S4650 = - (-1.0f * - (_S4649.differential_0 / _S4500));
    float2  _S4651 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4652;
    (&_S4652)->primal_0 = e2_6;
    (&_S4652)->differential_0 = _S4651;
    s_bwd_length_impl_0(&_S4652, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4653;
    (&_S4653)->primal_0 = e1_14;
    (&_S4653)->differential_0 = _S4651;
    s_bwd_length_impl_0(&_S4653, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S4654;
    (&_S4654)->primal_0 = e0_14;
    (&_S4654)->differential_0 = _S4651;
    s_bwd_length_impl_0(&_S4654, -0.0f);
    DiffPair_float_0 _S4655;
    (&_S4655)->primal_0 = _S4497;
    (&_S4655)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S4655, 0.0f);
    float _S4656 = - _S4655.differential_0;
    float2  _S4657 = _S4653.differential_0 + make_float2 (_S4495 * _S4656, _S4493 * _S4655.differential_0);
    float2  _S4658 = _S4654.differential_0 + make_float2 (_S4494 * _S4655.differential_0, _S4496 * _S4656);
    float2  _S4659 = - _S4652.differential_0 + _S4657;
    float _S4660 = fx_34 * (_S4642.differential_0 + _S4646.differential_0 + _S4659.x);
    float2  _S4661 = make_float2 (_S4660, fy_34 * (_S4634.differential_0 + _S4638.differential_0 + _S4659.y)) + make_float2 ((*dist_coeffs_38)[int(8)] * _S4660, (*dist_coeffs_38)[int(9)] * _S4660);
    float2  _S4662 = _S4479 * _S4661;
    float _S4663 = (*dist_coeffs_38)[int(4)] * _S4661.y;
    float _S4664 = (*dist_coeffs_38)[int(5)] * _S4661.x;
    float _S4665 = _S4662.x + _S4662.y;
    float _S4666 = r2_76 * _S4665;
    float _S4667 = r2_76 * _S4666;
    float _S4668 = (*dist_coeffs_38)[int(7)] * _S4661.y + _S4663 + (*dist_coeffs_38)[int(6)] * _S4661.x + _S4664 + _S4483 * _S4665 + _S4482 * _S4666 + _S4481 * _S4667 + (*dist_coeffs_38)[int(3)] * (r2_76 * _S4667);
    float _S4669 = v_76 * _S4668;
    float _S4670 = u_76 * _S4668;
    float2  _S4671 = (_S4484 * _S4661 + make_float2 (_S4452 * (v_76 * _S4661.y) + _S4486 * _S4664 + 2.0f * (u_76 * _S4664) + _S4449 * (v_76 * _S4661.x) + _S4670 + _S4670, _S4488 * _S4663 + 2.0f * (v_76 * _S4663) + _S4487 * _S4661.y + _S4485 * _S4661.x + _S4669 + _S4669)) / _S4480;
    float2  _S4672 = _S4476 * - _S4671;
    float2  _S4673 = _S4478 * _S4671;
    float2  _S4674 = - _S4657 + _S4658;
    float _S4675 = fx_34 * (_S4644.differential_0 + _S4648.differential_0 + _S4674.x);
    float2  _S4676 = make_float2 (_S4675, fy_34 * (_S4636.differential_0 + _S4640.differential_0 + _S4674.y)) + make_float2 ((*dist_coeffs_38)[int(8)] * _S4675, (*dist_coeffs_38)[int(9)] * _S4675);
    float2  _S4677 = _S4462 * _S4676;
    float _S4678 = (*dist_coeffs_38)[int(4)] * _S4676.y;
    float _S4679 = (*dist_coeffs_38)[int(5)] * _S4676.x;
    float _S4680 = _S4677.x + _S4677.y;
    float _S4681 = r2_75 * _S4680;
    float _S4682 = r2_75 * _S4681;
    float _S4683 = (*dist_coeffs_38)[int(7)] * _S4676.y + _S4678 + (*dist_coeffs_38)[int(6)] * _S4676.x + _S4679 + _S4466 * _S4680 + _S4465 * _S4681 + _S4464 * _S4682 + (*dist_coeffs_38)[int(3)] * (r2_75 * _S4682);
    float _S4684 = v_75 * _S4683;
    float _S4685 = u_75 * _S4683;
    float2  _S4686 = (_S4467 * _S4676 + make_float2 (_S4452 * (v_75 * _S4676.y) + _S4469 * _S4679 + 2.0f * (u_75 * _S4679) + _S4449 * (v_75 * _S4676.x) + _S4685 + _S4685, _S4471 * _S4678 + 2.0f * (v_75 * _S4678) + _S4470 * _S4676.y + _S4468 * _S4676.x + _S4684 + _S4684)) / _S4463;
    float2  _S4687 = _S4459 * - _S4686;
    float2  _S4688 = _S4461 * _S4686;
    float _S4689 = _S4687.x + _S4687.y;
    float2  _S4690 = _S4652.differential_0 + - _S4658;
    float _S4691 = fx_34 * (_S4643.differential_0 + _S4647.differential_0 + _S4690.x);
    float2  _S4692 = make_float2 (_S4691, fy_34 * (_S4635.differential_0 + _S4639.differential_0 + _S4690.y)) + make_float2 ((*dist_coeffs_38)[int(8)] * _S4691, (*dist_coeffs_38)[int(9)] * _S4691);
    float2  _S4693 = _S4443 * _S4692;
    float _S4694 = (*dist_coeffs_38)[int(4)] * _S4692.y;
    float _S4695 = (*dist_coeffs_38)[int(5)] * _S4692.x;
    float _S4696 = _S4693.x + _S4693.y;
    float _S4697 = r2_74 * _S4696;
    float _S4698 = r2_74 * _S4697;
    float _S4699 = (*dist_coeffs_38)[int(7)] * _S4692.y + _S4694 + (*dist_coeffs_38)[int(6)] * _S4692.x + _S4695 + _S4447 * _S4696 + _S4446 * _S4697 + _S4445 * _S4698 + (*dist_coeffs_38)[int(3)] * (r2_74 * _S4698);
    float _S4700 = v_74 * _S4699;
    float _S4701 = u_74 * _S4699;
    float2  _S4702 = (_S4448 * _S4692 + make_float2 (_S4452 * (v_74 * _S4692.y) + _S4451 * _S4695 + 2.0f * (u_74 * _S4695) + _S4449 * (v_74 * _S4692.x) + _S4701 + _S4701, _S4454 * _S4694 + 2.0f * (v_74 * _S4694) + _S4453 * _S4692.y + _S4450 * _S4692.x + _S4700 + _S4700)) / _S4444;
    float2  _S4703 = _S4440 * - _S4702;
    float2  _S4704 = _S4442 * _S4702;
    float _S4705 = _S4703.x + _S4703.y;
    float3  _S4706 = _S4552.differential_0 + _S4632.differential_0 + make_float3 (_S4673.x, _S4673.y, _S4672.x + _S4672.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4707;
    (&_S4707)->primal_0 = R_31;
    (&_S4707)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4708;
    (&_S4708)->primal_0 = vert2_7;
    (&_S4708)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4707, &_S4708, _S4706);
    float3  _S4709 = _S4551.differential_0 + _S4632.differential_0 + make_float3 (_S4688.x, _S4688.y, _S4689);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4710;
    (&_S4710)->primal_0 = R_31;
    (&_S4710)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4711;
    (&_S4711)->primal_0 = vert1_7;
    (&_S4711)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4710, &_S4711, _S4709);
    float3  _S4712 = _S4553 + _S4554 + _S4632.differential_0 + make_float3 (_S4704.x, _S4704.y, _S4705);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4713;
    (&_S4713)->primal_0 = R_31;
    (&_S4713)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4714;
    (&_S4714)->primal_0 = vert0_7;
    (&_S4714)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4713, &_S4714, _S4712);
    float3  _S4715 = (*v_verts_0)[int(2)] + _S4708.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4716;
    (&_S4716)->primal_0 = _S4434;
    (&_S4716)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4717;
    (&_S4717)->primal_0 = _S4439;
    (&_S4717)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4716, &_S4717, _S4715);
    float _S4718 = - _S4717.differential_0.y;
    float _S4719 = _S4438 * _S4717.differential_0.x;
    float _S4720 = - (_S4430 * _S4717.differential_0.x);
    float3  _S4721 = (*v_verts_0)[int(1)] + _S4711.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4722;
    (&_S4722)->primal_0 = _S4434;
    (&_S4722)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4723;
    (&_S4723)->primal_0 = _S4437;
    (&_S4723)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4722, &_S4723, _S4721);
    float _S4724 = _S4430 * _S4723.differential_0.x;
    float _S4725 = _S4436 * _S4723.differential_0.x;
    float3  _S4726 = (*v_verts_0)[int(0)] + _S4714.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4727;
    (&_S4727)->primal_0 = _S4434;
    (&_S4727)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4728;
    (&_S4728)->primal_0 = _S4435;
    (&_S4728)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4727, &_S4728, _S4726);
    Matrix<float, 3, 3>  _S4729 = transpose_0(_S4716.differential_0 + _S4722.differential_0 + _S4727.differential_0);
    float _S4730 = 2.0f * - _S4729.rows[int(2)].z;
    float _S4731 = 2.0f * _S4729.rows[int(2)].y;
    float _S4732 = 2.0f * _S4729.rows[int(2)].x;
    float _S4733 = 2.0f * _S4729.rows[int(1)].z;
    float _S4734 = 2.0f * - _S4729.rows[int(1)].y;
    float _S4735 = 2.0f * _S4729.rows[int(1)].x;
    float _S4736 = 2.0f * _S4729.rows[int(0)].z;
    float _S4737 = 2.0f * _S4729.rows[int(0)].y;
    float _S4738 = 2.0f * - _S4729.rows[int(0)].x;
    float _S4739 = - _S4735 + _S4737;
    float _S4740 = _S4732 + - _S4736;
    float _S4741 = - _S4731 + _S4733;
    float _S4742 = _S4731 + _S4733;
    float _S4743 = _S4732 + _S4736;
    float _S4744 = _S4735 + _S4737;
    float _S4745 = quat_31.w * (_S4734 + _S4738);
    float _S4746 = quat_31.z * (_S4730 + _S4738);
    float _S4747 = quat_31.y * (_S4730 + _S4734);
    float _S4748 = quat_31.x * _S4739 + quat_31.z * _S4742 + quat_31.y * _S4743 + _S4745 + _S4745;
    float _S4749 = quat_31.x * _S4740 + quat_31.w * _S4742 + quat_31.y * _S4744 + _S4746 + _S4746;
    float _S4750 = quat_31.x * _S4741 + quat_31.w * _S4743 + quat_31.z * _S4744 + _S4747 + _S4747;
    float _S4751 = quat_31.w * _S4739 + quat_31.z * _S4740 + quat_31.y * _S4741;
    float _S4752 = _S4720 + _S4724;
    float _S4753 = 0.5f * - _S4752;
    float _S4754 = _S4718 + _S4723.differential_0.y;
    DiffPair_float_0 _S4755;
    (&_S4755)->primal_0 = _S4431;
    (&_S4755)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4755, _S4754);
    float _S4756 = _S4753 + _S4755.differential_0;
    float _S4757 = _S4719 + _S4725 + _S4728.differential_0.x;
    DiffPair_float_0 _S4758;
    (&_S4758)->primal_0 = _S4429;
    (&_S4758)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S4758, _S4757);
    float _S4759 = _S4753 + _S4758.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4760;
    (&_S4760)->primal_0 = R_31;
    (&_S4760)->differential_0 = _S4626;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4761;
    (&_S4761)->primal_0 = mean_30;
    (&_S4761)->differential_0 = _S4546;
    s_bwd_prop_mul_1(&_S4760, &_S4761, _S4548.differential_0);
    float3  _S4762 = _S4628.differential_0 + _S4706 + _S4709 + _S4712 + _S4548.differential_0;
    Matrix<float, 3, 3>  _S4763 = _S4629 + _S4707.differential_0 + _S4710.differential_0 + _S4713.differential_0 + _S4760.differential_0;
    FixedArray<float3 , 2>  _S4764;
    _S4764[int(0)] = _S4546;
    _S4764[int(1)] = _S4546;
    _S4764[int(1)] = _S4560;
    _S4764[int(0)] = _S4565;
    FixedArray<float3 , 16>  _S4765;
    _S4765[int(0)] = _S4546;
    _S4765[int(1)] = _S4546;
    _S4765[int(2)] = _S4546;
    _S4765[int(3)] = _S4546;
    _S4765[int(4)] = _S4546;
    _S4765[int(5)] = _S4546;
    _S4765[int(6)] = _S4546;
    _S4765[int(7)] = _S4546;
    _S4765[int(8)] = _S4546;
    _S4765[int(9)] = _S4546;
    _S4765[int(10)] = _S4546;
    _S4765[int(11)] = _S4546;
    _S4765[int(12)] = _S4546;
    _S4765[int(13)] = _S4546;
    _S4765[int(14)] = _S4546;
    _S4765[int(15)] = _S4546;
    _S4765[int(15)] = _S4567;
    _S4765[int(14)] = _S4569;
    _S4765[int(13)] = _S4571;
    _S4765[int(12)] = _S4573;
    _S4765[int(11)] = _S4575;
    _S4765[int(10)] = _S4577;
    _S4765[int(9)] = _S4579;
    _S4765[int(8)] = _S4587;
    _S4765[int(7)] = _S4589;
    _S4765[int(6)] = _S4591;
    _S4765[int(5)] = _S4593;
    _S4765[int(4)] = _S4595;
    _S4765[int(3)] = _S4606;
    _S4765[int(2)] = _S4608;
    _S4765[int(1)] = _S4610;
    _S4765[int(0)] = _S4623;
    float2  _S4766 = make_float2 (0.0f, _S4650);
    float3  _S4767 = make_float3 (_S4759, _S4756, _S4752);
    float4  _S4768 = make_float4 (0.0f);
    *&((&_S4768)->w) = _S4748;
    *&((&_S4768)->z) = _S4749;
    *&((&_S4768)->y) = _S4750;
    *&((&_S4768)->x) = _S4751;
    *v_mean_9 = _S4624 + _S4715 + _S4721 + _S4726 + _S4761.differential_0;
    *v_quat_8 = _S4768;
    *v_scale_8 = _S4767;
    *v_hardness_4 = _S4766;
    *v_sh_coeffs_7 = _S4765;
    *v_ch_coeffs_2 = _S4764;
    *v_R_9 = _S4763;
    *v_t_8 = _S4762;
    return;
}

inline __device__ void projection_opaque_triangle_eval3d_fisheye_vjp(float3  mean_31, float4  quat_32, float3  scale_31, float2  hardness_15, FixedArray<float3 , 16>  * sh_coeffs_26, FixedArray<float3 , 2>  * ch_coeffs_11, Matrix<float, 3, 3>  R_32, float3  t_29, float fx_35, float fy_35, float cx_30, float cy_30, FixedArray<float, 10>  * dist_coeffs_39, uint image_width_26, uint image_height_26, float v_depth_10, FixedArray<float3 , 3>  * v_verts_1, FixedArray<float3 , 3>  * v_rgbs_1, float3  v_normal_3, float3  * v_mean_10, float4  * v_quat_9, float3  * v_scale_9, float2  * v_hardness_5, FixedArray<float3 , 16>  * v_sh_coeffs_8, FixedArray<float3 , 2>  * v_ch_coeffs_3, Matrix<float, 3, 3>  * v_R_10, float3  * v_t_9)
{
    float3  mean_c_26 = s_primal_ctx_mul_0(R_32, mean_31) + t_29;
    float _S4769 = scale_31.x;
    float _S4770 = s_primal_ctx_exp_1(_S4769);
    float _S4771 = scale_31.y;
    float _S4772 = s_primal_ctx_exp_1(_S4771);
    float sz_12 = scale_31.z - 0.5f * (_S4769 + _S4771);
    float _S4773 = quat_32.y;
    float x2_32 = _S4773 * _S4773;
    float y2_32 = quat_32.z * quat_32.z;
    float z2_58 = quat_32.w * quat_32.w;
    float xy_32 = quat_32.y * quat_32.z;
    float xz_32 = quat_32.y * quat_32.w;
    float yz_32 = quat_32.z * quat_32.w;
    float wx_32 = quat_32.x * quat_32.y;
    float wy_32 = quat_32.x * quat_32.z;
    float wz_32 = quat_32.x * quat_32.w;
    Matrix<float, 3, 3>  _S4774 = transpose_0(makeMatrix<float, 3, 3> (1.0f - 2.0f * (y2_32 + z2_58), 2.0f * (xy_32 + wz_32), 2.0f * (xz_32 - wy_32), 2.0f * (xy_32 - wz_32), 1.0f - 2.0f * (x2_32 + z2_58), 2.0f * (yz_32 + wx_32), 2.0f * (xz_32 + wy_32), 2.0f * (yz_32 - wx_32), 1.0f - 2.0f * (x2_32 + y2_32)));
    float3  _S4775 = make_float3 (_S4770, 0.0f, 0.0f);
    float3  vert0_8 = s_primal_ctx_mul_0(_S4774, _S4775) + mean_31;
    float _S4776 = -0.5f + sz_12;
    float3  _S4777 = make_float3 (_S4770 * _S4776, _S4772, 0.0f);
    float3  vert1_8 = s_primal_ctx_mul_0(_S4774, _S4777) + mean_31;
    float _S4778 = -0.5f - sz_12;
    float3  _S4779 = make_float3 (_S4770 * _S4778, - _S4772, 0.0f);
    float3  vert2_8 = s_primal_ctx_mul_0(_S4774, _S4779) + mean_31;
    float3  vert0_c_11 = s_primal_ctx_mul_0(R_32, vert0_8) + t_29;
    float3  vert1_c_11 = s_primal_ctx_mul_0(R_32, vert1_8) + t_29;
    float3  vert2_c_11 = s_primal_ctx_mul_0(R_32, vert2_8) + t_29;
    float2  _S4780 = float2 {vert0_c_11.x, vert0_c_11.y};
    float _S4781 = length_0(_S4780);
    float _S4782 = vert0_c_11.z;
    float _S4783 = s_primal_ctx_atan2_0(_S4781, _S4782);
    bool _S4784 = _S4783 < 0.00100000004749745f;
    float k_12;
    float _S4785;
    float _S4786;
    float _S4787;
    if(_S4784)
    {
        float _S4788 = 1.0f - _S4783 * _S4783 / 3.0f;
        float _S4789 = _S4782 * _S4782;
        k_12 = _S4788 / _S4782;
        _S4785 = 0.0f;
        _S4786 = _S4789;
        _S4787 = _S4788;
    }
    else
    {
        float _S4790 = _S4781 * _S4781;
        k_12 = _S4783 / _S4781;
        _S4785 = _S4790;
        _S4786 = 0.0f;
        _S4787 = 0.0f;
    }
    float2  _S4791 = make_float2 (k_12);
    float2  _S4792 = _S4780 * make_float2 (k_12);
    float u_77 = _S4792.x;
    float v_77 = _S4792.y;
    float r2_77 = u_77 * u_77 + v_77 * v_77;
    float _S4793 = (*dist_coeffs_39)[int(2)] + r2_77 * (*dist_coeffs_39)[int(3)];
    float _S4794 = (*dist_coeffs_39)[int(1)] + r2_77 * _S4793;
    float _S4795 = (*dist_coeffs_39)[int(0)] + r2_77 * _S4794;
    float radial_10 = 1.0f + r2_77 * _S4795;
    float2  _S4796 = make_float2 (radial_10);
    float _S4797 = 2.0f * (*dist_coeffs_39)[int(4)];
    float _S4798 = _S4797 * u_77;
    float _S4799 = 2.0f * u_77;
    float _S4800 = 2.0f * (*dist_coeffs_39)[int(5)];
    float _S4801 = _S4800 * u_77;
    float _S4802 = 2.0f * v_77;
    float2  _S4803 = _S4792 * make_float2 (radial_10) + make_float2 (_S4798 * v_77 + (*dist_coeffs_39)[int(5)] * (r2_77 + _S4799 * u_77) + (*dist_coeffs_39)[int(6)] * r2_77, _S4801 * v_77 + (*dist_coeffs_39)[int(4)] * (r2_77 + _S4802 * v_77) + (*dist_coeffs_39)[int(7)] * r2_77);
    float2  _S4804 = _S4803 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4803.x + (*dist_coeffs_39)[int(9)] * _S4803.y, 0.0f);
    float _S4805 = fx_35 * _S4804.x + cx_30;
    float _S4806 = fy_35 * _S4804.y + cy_30;
    float2  uv0_11 = make_float2 (_S4805, _S4806);
    float2  _S4807 = float2 {vert1_c_11.x, vert1_c_11.y};
    float _S4808 = length_0(_S4807);
    float _S4809 = vert1_c_11.z;
    float _S4810 = s_primal_ctx_atan2_0(_S4808, _S4809);
    bool _S4811 = _S4810 < 0.00100000004749745f;
    float _S4812;
    float _S4813;
    float _S4814;
    if(_S4811)
    {
        float _S4815 = 1.0f - _S4810 * _S4810 / 3.0f;
        float _S4816 = _S4809 * _S4809;
        k_12 = _S4815 / _S4809;
        _S4812 = 0.0f;
        _S4813 = _S4816;
        _S4814 = _S4815;
    }
    else
    {
        float _S4817 = _S4808 * _S4808;
        k_12 = _S4810 / _S4808;
        _S4812 = _S4817;
        _S4813 = 0.0f;
        _S4814 = 0.0f;
    }
    float2  _S4818 = make_float2 (k_12);
    float2  _S4819 = _S4807 * make_float2 (k_12);
    float u_78 = _S4819.x;
    float v_78 = _S4819.y;
    float r2_78 = u_78 * u_78 + v_78 * v_78;
    float _S4820 = (*dist_coeffs_39)[int(2)] + r2_78 * (*dist_coeffs_39)[int(3)];
    float _S4821 = (*dist_coeffs_39)[int(1)] + r2_78 * _S4820;
    float _S4822 = (*dist_coeffs_39)[int(0)] + r2_78 * _S4821;
    float radial_11 = 1.0f + r2_78 * _S4822;
    float2  _S4823 = make_float2 (radial_11);
    float _S4824 = _S4797 * u_78;
    float _S4825 = 2.0f * u_78;
    float _S4826 = _S4800 * u_78;
    float _S4827 = 2.0f * v_78;
    float2  _S4828 = _S4819 * make_float2 (radial_11) + make_float2 (_S4824 * v_78 + (*dist_coeffs_39)[int(5)] * (r2_78 + _S4825 * u_78) + (*dist_coeffs_39)[int(6)] * r2_78, _S4826 * v_78 + (*dist_coeffs_39)[int(4)] * (r2_78 + _S4827 * v_78) + (*dist_coeffs_39)[int(7)] * r2_78);
    float2  _S4829 = _S4828 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4828.x + (*dist_coeffs_39)[int(9)] * _S4828.y, 0.0f);
    float _S4830 = fx_35 * _S4829.x + cx_30;
    float _S4831 = fy_35 * _S4829.y + cy_30;
    float2  uv1_11 = make_float2 (_S4830, _S4831);
    float2  _S4832 = float2 {vert2_c_11.x, vert2_c_11.y};
    float _S4833 = length_0(_S4832);
    float _S4834 = vert2_c_11.z;
    float _S4835 = s_primal_ctx_atan2_0(_S4833, _S4834);
    bool _S4836 = _S4835 < 0.00100000004749745f;
    float _S4837;
    float _S4838;
    float _S4839;
    if(_S4836)
    {
        float _S4840 = 1.0f - _S4835 * _S4835 / 3.0f;
        float _S4841 = _S4834 * _S4834;
        k_12 = _S4840 / _S4834;
        _S4837 = 0.0f;
        _S4838 = _S4841;
        _S4839 = _S4840;
    }
    else
    {
        float _S4842 = _S4833 * _S4833;
        k_12 = _S4835 / _S4833;
        _S4837 = _S4842;
        _S4838 = 0.0f;
        _S4839 = 0.0f;
    }
    float2  _S4843 = make_float2 (k_12);
    float2  _S4844 = _S4832 * make_float2 (k_12);
    float u_79 = _S4844.x;
    float v_79 = _S4844.y;
    float r2_79 = u_79 * u_79 + v_79 * v_79;
    float _S4845 = (*dist_coeffs_39)[int(2)] + r2_79 * (*dist_coeffs_39)[int(3)];
    float _S4846 = (*dist_coeffs_39)[int(1)] + r2_79 * _S4845;
    float _S4847 = (*dist_coeffs_39)[int(0)] + r2_79 * _S4846;
    float radial_12 = 1.0f + r2_79 * _S4847;
    float2  _S4848 = make_float2 (radial_12);
    float _S4849 = _S4797 * u_79;
    float _S4850 = 2.0f * u_79;
    float _S4851 = _S4800 * u_79;
    float _S4852 = 2.0f * v_79;
    float2  _S4853 = _S4844 * make_float2 (radial_12) + make_float2 (_S4849 * v_79 + (*dist_coeffs_39)[int(5)] * (r2_79 + _S4850 * u_79) + (*dist_coeffs_39)[int(6)] * r2_79, _S4851 * v_79 + (*dist_coeffs_39)[int(4)] * (r2_79 + _S4852 * v_79) + (*dist_coeffs_39)[int(7)] * r2_79);
    float2  _S4854 = _S4853 + make_float2 ((*dist_coeffs_39)[int(8)] * _S4853.x + (*dist_coeffs_39)[int(9)] * _S4853.y, 0.0f);
    float _S4855 = fx_35 * _S4854.x + cx_30;
    float _S4856 = fy_35 * _S4854.y + cy_30;
    float2  uv2_11 = make_float2 (_S4855, _S4856);
    float2  e0_15 = uv1_11 - uv0_11;
    float2  e1_15 = uv2_11 - uv1_11;
    float2  e2_7 = uv0_11 - uv2_11;
    float _S4857 = e0_15.x;
    float _S4858 = e1_15.y;
    float _S4859 = e0_15.y;
    float _S4860 = e1_15.x;
    float _S4861 = _S4857 * _S4858 - _S4859 * _S4860;
    float _S4862 = 1.0f - hardness_15.y;
    float _S4863 = -1.0f / _S4862;
    float _S4864 = _S4862 * _S4862;
    float _S4865 = s_primal_ctx_max_0(_S4805, _S4830);
    float _S4866 = s_primal_ctx_min_0(_S4805, _S4830);
    float _S4867 = s_primal_ctx_max_0(_S4806, _S4831);
    float _S4868 = s_primal_ctx_min_0(_S4806, _S4831);
    float3  _S4869 = vert0_c_11 + vert1_c_11 + vert2_c_11;
    float _S4870 = length_1(_S4869) / 3.0f + 9.999999960041972e-13f;
    Matrix<float, 3, 3>  _S4871 = transpose_0(R_32);
    float3  _S4872 = mean_31 - - s_primal_ctx_mul_0(_S4871, t_29);
    float _S4873 = _S4872.x;
    float _S4874 = _S4872.y;
    float _S4875 = _S4872.z;
    float _S4876 = _S4873 * _S4873 + _S4874 * _S4874 + _S4875 * _S4875;
    float _S4877 = s_primal_ctx_sqrt_0(_S4876);
    float x_62 = _S4873 / _S4877;
    float3  _S4878 = make_float3 (x_62);
    float _S4879 = _S4877 * _S4877;
    float y_29 = _S4874 / _S4877;
    float z_26 = _S4875 / _S4877;
    float3  _S4880 = make_float3 (z_26);
    float _S4881 = - y_29;
    float3  _S4882 = make_float3 (_S4881);
    float z2_59 = z_26 * z_26;
    float fTmp0B_26 = -1.09254848957061768f * z_26;
    float fC1_26 = x_62 * x_62 - y_29 * y_29;
    float _S4883 = 2.0f * x_62;
    float fS1_26 = _S4883 * y_29;
    float pSH6_8 = 0.94617468118667603f * z2_59 - 0.31539157032966614f;
    float3  _S4884 = make_float3 (pSH6_8);
    float pSH7_8 = fTmp0B_26 * x_62;
    float3  _S4885 = make_float3 (pSH7_8);
    float pSH5_8 = fTmp0B_26 * y_29;
    float3  _S4886 = make_float3 (pSH5_8);
    float pSH8_8 = 0.54627424478530884f * fC1_26;
    float3  _S4887 = make_float3 (pSH8_8);
    float pSH4_8 = 0.54627424478530884f * fS1_26;
    float3  _S4888 = make_float3 (pSH4_8);
    float fTmp0C_26 = -2.28522896766662598f * z2_59 + 0.4570457935333252f;
    float fTmp1B_26 = 1.44530570507049561f * z_26;
    float _S4889 = 1.86588168144226074f * z2_59 - 1.11952900886535645f;
    float pSH12_8 = z_26 * _S4889;
    float3  _S4890 = make_float3 (pSH12_8);
    float pSH13_8 = fTmp0C_26 * x_62;
    float3  _S4891 = make_float3 (pSH13_8);
    float pSH11_8 = fTmp0C_26 * y_29;
    float3  _S4892 = make_float3 (pSH11_8);
    float pSH14_8 = fTmp1B_26 * fC1_26;
    float3  _S4893 = make_float3 (pSH14_8);
    float pSH10_8 = fTmp1B_26 * fS1_26;
    float3  _S4894 = make_float3 (pSH10_8);
    float pSH15_8 = -0.59004360437393188f * (x_62 * fC1_26 - y_29 * fS1_26);
    float3  _S4895 = make_float3 (pSH15_8);
    float pSH9_8 = -0.59004360437393188f * (x_62 * fS1_26 + y_29 * fC1_26);
    float3  _S4896 = make_float3 (pSH9_8);
    float3  color_12 = make_float3 (0.282094806432724f) * (*sh_coeffs_26)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S4881) * (*sh_coeffs_26)[int(1)] + make_float3 (z_26) * (*sh_coeffs_26)[int(2)] - make_float3 (x_62) * (*sh_coeffs_26)[int(3)]) + (make_float3 (pSH4_8) * (*sh_coeffs_26)[int(4)] + make_float3 (pSH5_8) * (*sh_coeffs_26)[int(5)] + make_float3 (pSH6_8) * (*sh_coeffs_26)[int(6)] + make_float3 (pSH7_8) * (*sh_coeffs_26)[int(7)] + make_float3 (pSH8_8) * (*sh_coeffs_26)[int(8)]) + (make_float3 (pSH9_8) * (*sh_coeffs_26)[int(9)] + make_float3 (pSH10_8) * (*sh_coeffs_26)[int(10)] + make_float3 (pSH11_8) * (*sh_coeffs_26)[int(11)] + make_float3 (pSH12_8) * (*sh_coeffs_26)[int(12)] + make_float3 (pSH13_8) * (*sh_coeffs_26)[int(13)] + make_float3 (pSH14_8) * (*sh_coeffs_26)[int(14)] + make_float3 (pSH15_8) * (*sh_coeffs_26)[int(15)]);
    float3  _S4897 = color_12 + (*ch_coeffs_11)[int(0)] + make_float3 (0.5f);
    float3  _S4898 = make_float3 (0.0f);
    float3  _S4899 = color_12 - (*ch_coeffs_11)[int(0)] * make_float3 (0.5f);
    float _S4900 = s_primal_ctx_sqrt_0(0.75f);
    float3  _S4901 = make_float3 (_S4900);
    float3  _S4902 = (*ch_coeffs_11)[int(1)] * make_float3 (_S4900);
    float3  _S4903 = _S4899 + _S4902 + make_float3 (0.5f);
    float3  _S4904 = _S4899 - _S4902 + make_float3 (0.5f);
    float3  _S4905 = vert1_c_11 - vert0_c_11;
    float3  _S4906 = vert2_c_11 - vert0_c_11;
    float3  _S4907 = s_primal_ctx_cross_0(_S4905, _S4906);
    float3  _S4908 = normalize_0(_S4907);
    float3  _S4909 = make_float3 (float(- (F32_sign((s_primal_ctx_dot_0(_S4908, mean_c_26)))))) * v_normal_3;
    float3  _S4910 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4911;
    (&_S4911)->primal_0 = _S4908;
    (&_S4911)->differential_0 = _S4910;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4912;
    (&_S4912)->primal_0 = mean_c_26;
    (&_S4912)->differential_0 = _S4910;
    s_bwd_prop_dot_0(&_S4911, &_S4912, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4913 = _S4912;
    float3  _S4914 = _S4909 + _S4911.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4915;
    (&_S4915)->primal_0 = _S4907;
    (&_S4915)->differential_0 = _S4910;
    s_bwd_normalize_impl_0(&_S4915, _S4914);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4916;
    (&_S4916)->primal_0 = _S4905;
    (&_S4916)->differential_0 = _S4910;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4917;
    (&_S4917)->primal_0 = _S4906;
    (&_S4917)->differential_0 = _S4910;
    s_bwd_prop_cross_0(&_S4916, &_S4917, _S4915.differential_0);
    float3  _S4918 = - _S4917.differential_0;
    float3  _S4919 = - _S4916.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4920;
    (&_S4920)->primal_0 = _S4904;
    (&_S4920)->differential_0 = _S4910;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4921;
    (&_S4921)->primal_0 = _S4898;
    (&_S4921)->differential_0 = _S4910;
    s_bwd_prop_max_0(&_S4920, &_S4921, (*v_rgbs_1)[int(2)]);
    float3  _S4922 = - _S4920.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4923;
    (&_S4923)->primal_0 = _S4903;
    (&_S4923)->differential_0 = _S4910;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4924;
    (&_S4924)->primal_0 = _S4898;
    (&_S4924)->differential_0 = _S4910;
    s_bwd_prop_max_0(&_S4923, &_S4924, (*v_rgbs_1)[int(1)]);
    float3  _S4925 = _S4901 * (_S4922 + _S4923.differential_0);
    float3  _S4926 = _S4920.differential_0 + _S4923.differential_0;
    float3  _S4927 = make_float3 (0.5f) * - _S4926;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4928;
    (&_S4928)->primal_0 = _S4897;
    (&_S4928)->differential_0 = _S4910;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4929;
    (&_S4929)->primal_0 = _S4898;
    (&_S4929)->differential_0 = _S4910;
    s_bwd_prop_max_0(&_S4928, &_S4929, (*v_rgbs_1)[int(0)]);
    float3  _S4930 = _S4927 + _S4928.differential_0;
    float3  _S4931 = _S4926 + _S4928.differential_0;
    float3  _S4932 = _S4895 * _S4931;
    float3  _S4933 = (*sh_coeffs_26)[int(15)] * _S4931;
    float3  _S4934 = _S4893 * _S4931;
    float3  _S4935 = (*sh_coeffs_26)[int(14)] * _S4931;
    float3  _S4936 = _S4891 * _S4931;
    float3  _S4937 = (*sh_coeffs_26)[int(13)] * _S4931;
    float3  _S4938 = _S4890 * _S4931;
    float3  _S4939 = (*sh_coeffs_26)[int(12)] * _S4931;
    float3  _S4940 = _S4892 * _S4931;
    float3  _S4941 = (*sh_coeffs_26)[int(11)] * _S4931;
    float3  _S4942 = _S4894 * _S4931;
    float3  _S4943 = (*sh_coeffs_26)[int(10)] * _S4931;
    float3  _S4944 = _S4896 * _S4931;
    float3  _S4945 = (*sh_coeffs_26)[int(9)] * _S4931;
    float s_diff_fS2_T_8 = -0.59004360437393188f * (_S4945.x + _S4945.y + _S4945.z);
    float s_diff_fC2_T_8 = -0.59004360437393188f * (_S4933.x + _S4933.y + _S4933.z);
    float _S4946 = _S4943.x + _S4943.y + _S4943.z;
    float _S4947 = _S4935.x + _S4935.y + _S4935.z;
    float _S4948 = _S4941.x + _S4941.y + _S4941.z;
    float _S4949 = _S4937.x + _S4937.y + _S4937.z;
    float _S4950 = _S4939.x + _S4939.y + _S4939.z;
    float _S4951 = - s_diff_fC2_T_8;
    float3  _S4952 = _S4887 * _S4931;
    float3  _S4953 = (*sh_coeffs_26)[int(8)] * _S4931;
    float3  _S4954 = _S4885 * _S4931;
    float3  _S4955 = (*sh_coeffs_26)[int(7)] * _S4931;
    float3  _S4956 = _S4884 * _S4931;
    float3  _S4957 = (*sh_coeffs_26)[int(6)] * _S4931;
    float3  _S4958 = _S4886 * _S4931;
    float3  _S4959 = (*sh_coeffs_26)[int(5)] * _S4931;
    float3  _S4960 = _S4888 * _S4931;
    float3  _S4961 = (*sh_coeffs_26)[int(4)] * _S4931;
    float _S4962 = _S4959.x + _S4959.y + _S4959.z;
    float _S4963 = _S4955.x + _S4955.y + _S4955.z;
    float _S4964 = fTmp1B_26 * _S4946 + x_62 * s_diff_fS2_T_8 + y_29 * _S4951 + 0.54627424478530884f * (_S4961.x + _S4961.y + _S4961.z);
    float _S4965 = fTmp1B_26 * _S4947 + y_29 * s_diff_fS2_T_8 + x_62 * s_diff_fC2_T_8 + 0.54627424478530884f * (_S4953.x + _S4953.y + _S4953.z);
    float _S4966 = y_29 * - _S4965;
    float _S4967 = x_62 * _S4965;
    float _S4968 = z_26 * (1.86588168144226074f * (z_26 * _S4950) + -2.28522896766662598f * (y_29 * _S4948 + x_62 * _S4949) + 0.94617468118667603f * (_S4957.x + _S4957.y + _S4957.z));
    float3  _S4969 = make_float3 (0.48860251903533936f) * _S4931;
    float3  _S4970 = - _S4969;
    float3  _S4971 = _S4878 * _S4970;
    float3  _S4972 = (*sh_coeffs_26)[int(3)] * _S4970;
    float3  _S4973 = _S4880 * _S4969;
    float3  _S4974 = (*sh_coeffs_26)[int(2)] * _S4969;
    float3  _S4975 = _S4882 * _S4969;
    float3  _S4976 = (*sh_coeffs_26)[int(1)] * _S4969;
    float _S4977 = (_S4889 * _S4950 + 1.44530570507049561f * (fS1_26 * _S4946 + fC1_26 * _S4947) + -1.09254848957061768f * (y_29 * _S4962 + x_62 * _S4963) + _S4968 + _S4968 + _S4974.x + _S4974.y + _S4974.z) / _S4879;
    float _S4978 = _S4877 * _S4977;
    float _S4979 = (fTmp0C_26 * _S4948 + fC1_26 * s_diff_fS2_T_8 + fS1_26 * _S4951 + fTmp0B_26 * _S4962 + _S4883 * _S4964 + _S4966 + _S4966 + - (_S4976.x + _S4976.y + _S4976.z)) / _S4879;
    float _S4980 = _S4877 * _S4979;
    float _S4981 = (fTmp0C_26 * _S4949 + fS1_26 * s_diff_fS2_T_8 + fC1_26 * s_diff_fC2_T_8 + fTmp0B_26 * _S4963 + 2.0f * (y_29 * _S4964) + _S4967 + _S4967 + _S4972.x + _S4972.y + _S4972.z) / _S4879;
    float _S4982 = _S4877 * _S4981;
    float _S4983 = _S4875 * - _S4977 + _S4874 * - _S4979 + _S4873 * - _S4981;
    DiffPair_float_0 _S4984;
    (&_S4984)->primal_0 = _S4876;
    (&_S4984)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S4984, _S4983);
    float _S4985 = _S4875 * _S4984.differential_0;
    float _S4986 = _S4874 * _S4984.differential_0;
    float _S4987 = _S4873 * _S4984.differential_0;
    float3  _S4988 = make_float3 (0.282094806432724f) * _S4931;
    float3  _S4989 = make_float3 (_S4982 + _S4987 + _S4987, _S4980 + _S4986 + _S4986, _S4978 + _S4985 + _S4985);
    float3  _S4990 = - - _S4989;
    Matrix<float, 3, 3>  _S4991 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S4992;
    (&_S4992)->primal_0 = _S4871;
    (&_S4992)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4993;
    (&_S4993)->primal_0 = t_29;
    (&_S4993)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S4992, &_S4993, _S4990);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4994 = _S4993;
    Matrix<float, 3, 3>  _S4995 = transpose_0(_S4992.differential_0);
    DiffPair_float_0 _S4996;
    (&_S4996)->primal_0 = _S4870;
    (&_S4996)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S4996, v_depth_10);
    float _S4997 = 0.3333333432674408f * _S4996.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S4998;
    (&_S4998)->primal_0 = _S4869;
    (&_S4998)->differential_0 = _S4910;
    s_bwd_length_impl_1(&_S4998, _S4997);
    DiffPair_float_0 _S4999;
    (&_S4999)->primal_0 = _S4868;
    (&_S4999)->differential_0 = 0.0f;
    DiffPair_float_0 _S5000;
    (&_S5000)->primal_0 = _S4856;
    (&_S5000)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S4999, &_S5000, 0.0f);
    DiffPair_float_0 _S5001;
    (&_S5001)->primal_0 = _S4806;
    (&_S5001)->differential_0 = 0.0f;
    DiffPair_float_0 _S5002;
    (&_S5002)->primal_0 = _S4831;
    (&_S5002)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5001, &_S5002, _S4999.differential_0);
    DiffPair_float_0 _S5003;
    (&_S5003)->primal_0 = _S4867;
    (&_S5003)->differential_0 = 0.0f;
    DiffPair_float_0 _S5004;
    (&_S5004)->primal_0 = _S4856;
    (&_S5004)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5003, &_S5004, 0.0f);
    DiffPair_float_0 _S5005;
    (&_S5005)->primal_0 = _S4806;
    (&_S5005)->differential_0 = 0.0f;
    DiffPair_float_0 _S5006;
    (&_S5006)->primal_0 = _S4831;
    (&_S5006)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5005, &_S5006, _S5003.differential_0);
    DiffPair_float_0 _S5007;
    (&_S5007)->primal_0 = _S4866;
    (&_S5007)->differential_0 = 0.0f;
    DiffPair_float_0 _S5008;
    (&_S5008)->primal_0 = _S4855;
    (&_S5008)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5007, &_S5008, 0.0f);
    DiffPair_float_0 _S5009;
    (&_S5009)->primal_0 = _S4805;
    (&_S5009)->differential_0 = 0.0f;
    DiffPair_float_0 _S5010;
    (&_S5010)->primal_0 = _S4830;
    (&_S5010)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5009, &_S5010, _S5007.differential_0);
    DiffPair_float_0 _S5011;
    (&_S5011)->primal_0 = _S4865;
    (&_S5011)->differential_0 = 0.0f;
    DiffPair_float_0 _S5012;
    (&_S5012)->primal_0 = _S4855;
    (&_S5012)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5011, &_S5012, 0.0f);
    DiffPair_float_0 _S5013;
    (&_S5013)->primal_0 = _S4805;
    (&_S5013)->differential_0 = 0.0f;
    DiffPair_float_0 _S5014;
    (&_S5014)->primal_0 = _S4830;
    (&_S5014)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5013, &_S5014, _S5011.differential_0);
    DiffPair_float_0 _S5015;
    (&_S5015)->primal_0 = _S4863;
    (&_S5015)->differential_0 = 0.0f;
    s_bwd_prop_exp2_0(&_S5015, 0.0f);
    float _S5016 = - (-1.0f * - (_S5015.differential_0 / _S4864));
    float2  _S5017 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5018;
    (&_S5018)->primal_0 = e2_7;
    (&_S5018)->differential_0 = _S5017;
    s_bwd_length_impl_0(&_S5018, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5019;
    (&_S5019)->primal_0 = e1_15;
    (&_S5019)->differential_0 = _S5017;
    s_bwd_length_impl_0(&_S5019, -0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5020;
    (&_S5020)->primal_0 = e0_15;
    (&_S5020)->differential_0 = _S5017;
    s_bwd_length_impl_0(&_S5020, -0.0f);
    DiffPair_float_0 _S5021;
    (&_S5021)->primal_0 = _S4861;
    (&_S5021)->differential_0 = 0.0f;
    s_bwd_prop_abs_0(&_S5021, 0.0f);
    float _S5022 = - _S5021.differential_0;
    float2  _S5023 = _S5019.differential_0 + make_float2 (_S4859 * _S5022, _S4857 * _S5021.differential_0);
    float2  _S5024 = - _S5023;
    float2  _S5025 = _S5020.differential_0 + make_float2 (_S4858 * _S5021.differential_0, _S4860 * _S5022);
    float2  _S5026 = - _S5025;
    float2  _S5027 = - _S5018.differential_0 + _S5023;
    float _S5028 = fx_35 * (_S5008.differential_0 + _S5012.differential_0 + _S5027.x);
    float2  _S5029 = make_float2 (_S5028, fy_35 * (_S5000.differential_0 + _S5004.differential_0 + _S5027.y)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S5028, (*dist_coeffs_39)[int(9)] * _S5028);
    float2  _S5030 = _S4844 * _S5029;
    float2  _S5031 = _S4848 * _S5029;
    float _S5032 = (*dist_coeffs_39)[int(4)] * _S5029.y;
    float _S5033 = (*dist_coeffs_39)[int(5)] * _S5029.x;
    float _S5034 = _S5030.x + _S5030.y;
    float _S5035 = r2_79 * _S5034;
    float _S5036 = r2_79 * _S5035;
    float _S5037 = (*dist_coeffs_39)[int(7)] * _S5029.y + _S5032 + (*dist_coeffs_39)[int(6)] * _S5029.x + _S5033 + _S4847 * _S5034 + _S4846 * _S5035 + _S4845 * _S5036 + (*dist_coeffs_39)[int(3)] * (r2_79 * _S5036);
    float _S5038 = v_79 * _S5037;
    float _S5039 = u_79 * _S5037;
    float _S5040 = _S4852 * _S5032 + 2.0f * (v_79 * _S5032) + _S4851 * _S5029.y + _S4849 * _S5029.x + _S5038 + _S5038;
    float _S5041 = _S4800 * (v_79 * _S5029.y) + _S4850 * _S5033 + 2.0f * (u_79 * _S5033) + _S4797 * (v_79 * _S5029.x) + _S5039 + _S5039;
    float3  _S5042 = _S4916.differential_0 + _S4998.differential_0;
    float3  _S5043 = _S4918 + _S4919 + _S4998.differential_0;
    float3  _S5044 = _S4917.differential_0 + _S4998.differential_0;
    FixedArray<float3 , 2>  _S5045;
    _S5045[int(0)] = _S4910;
    _S5045[int(1)] = _S4910;
    _S5045[int(1)] = _S4925;
    _S5045[int(0)] = _S4930;
    float3  _S5046 = _S5045[int(0)];
    float3  _S5047 = _S5045[int(1)];
    FixedArray<float3 , 16>  _S5048;
    _S5048[int(0)] = _S4910;
    _S5048[int(1)] = _S4910;
    _S5048[int(2)] = _S4910;
    _S5048[int(3)] = _S4910;
    _S5048[int(4)] = _S4910;
    _S5048[int(5)] = _S4910;
    _S5048[int(6)] = _S4910;
    _S5048[int(7)] = _S4910;
    _S5048[int(8)] = _S4910;
    _S5048[int(9)] = _S4910;
    _S5048[int(10)] = _S4910;
    _S5048[int(11)] = _S4910;
    _S5048[int(12)] = _S4910;
    _S5048[int(13)] = _S4910;
    _S5048[int(14)] = _S4910;
    _S5048[int(15)] = _S4910;
    _S5048[int(7)] = _S4954;
    _S5048[int(0)] = _S4988;
    _S5048[int(1)] = _S4975;
    _S5048[int(2)] = _S4973;
    _S5048[int(3)] = _S4971;
    _S5048[int(4)] = _S4960;
    _S5048[int(5)] = _S4958;
    _S5048[int(6)] = _S4956;
    _S5048[int(15)] = _S4932;
    _S5048[int(8)] = _S4952;
    _S5048[int(9)] = _S4944;
    _S5048[int(10)] = _S4942;
    _S5048[int(11)] = _S4940;
    _S5048[int(12)] = _S4938;
    _S5048[int(13)] = _S4936;
    _S5048[int(14)] = _S4934;
    float3  _S5049 = _S5048[int(0)];
    float3  _S5050 = _S5048[int(1)];
    float3  _S5051 = _S5048[int(2)];
    float3  _S5052 = _S5048[int(3)];
    float3  _S5053 = _S5048[int(4)];
    float3  _S5054 = _S5048[int(5)];
    float3  _S5055 = _S5048[int(6)];
    float3  _S5056 = _S5048[int(7)];
    float3  _S5057 = _S5048[int(8)];
    float3  _S5058 = _S5048[int(9)];
    float3  _S5059 = _S5048[int(10)];
    float3  _S5060 = _S5048[int(11)];
    float3  _S5061 = _S5048[int(12)];
    float3  _S5062 = _S5048[int(13)];
    float3  _S5063 = _S5048[int(14)];
    float3  _S5064 = _S5048[int(15)];
    float _S5065 = _S5009.differential_0 + _S5013.differential_0;
    float2  _S5066 = _S5018.differential_0 + _S5026;
    float _S5067 = _S5002.differential_0 + _S5006.differential_0;
    float _S5068 = _S5001.differential_0 + _S5005.differential_0;
    float2  _S5069 = _S5024 + _S5025;
    float _S5070 = _S5010.differential_0 + _S5014.differential_0;
    float2  _S5071 = make_float2 (0.0f, _S5016);
    float2  _S5072 = _S5031 + make_float2 (_S5041, _S5040);
    float2  _S5073 = _S4832 * _S5072;
    float2  _S5074 = _S4843 * _S5072;
    float _S5075 = _S5073.x + _S5073.y;
    if(_S4836)
    {
        float _S5076 = _S5075 / _S4838;
        float _S5077 = _S4839 * - _S5076;
        float _S5078 = _S4835 * (0.3333333432674408f * - (_S4834 * _S5076));
        k_12 = _S5078 + _S5078;
        _S4837 = _S5077;
        _S4838 = 0.0f;
    }
    else
    {
        float _S5079 = _S5075 / _S4837;
        float _S5080 = _S4835 * - _S5079;
        k_12 = _S4833 * _S5079;
        _S4837 = 0.0f;
        _S4838 = _S5080;
    }
    DiffPair_float_0 _S5081;
    (&_S5081)->primal_0 = _S4833;
    (&_S5081)->differential_0 = 0.0f;
    DiffPair_float_0 _S5082;
    (&_S5082)->primal_0 = _S4834;
    (&_S5082)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5081, &_S5082, k_12);
    float _S5083 = _S5082.differential_0 + _S4837;
    float _S5084 = _S5081.differential_0 + _S4838;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5085;
    (&_S5085)->primal_0 = _S4832;
    (&_S5085)->differential_0 = _S5017;
    s_bwd_length_impl_0(&_S5085, _S5084);
    float2  _S5086 = _S5085.differential_0 + _S5074;
    float _S5087 = fx_35 * (_S5069.x + _S5070);
    float2  _S5088 = make_float2 (_S5087, fy_35 * (_S5069.y + _S5067)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S5087, (*dist_coeffs_39)[int(9)] * _S5087);
    float2  _S5089 = _S4819 * _S5088;
    float _S5090 = (*dist_coeffs_39)[int(4)] * _S5088.y;
    float _S5091 = (*dist_coeffs_39)[int(5)] * _S5088.x;
    float _S5092 = _S5089.x + _S5089.y;
    float _S5093 = r2_78 * _S5092;
    float _S5094 = r2_78 * _S5093;
    float _S5095 = (*dist_coeffs_39)[int(7)] * _S5088.y + _S5090 + (*dist_coeffs_39)[int(6)] * _S5088.x + _S5091 + _S4822 * _S5092 + _S4821 * _S5093 + _S4820 * _S5094 + (*dist_coeffs_39)[int(3)] * (r2_78 * _S5094);
    float _S5096 = v_78 * _S5095;
    float _S5097 = u_78 * _S5095;
    float3  _S5098 = _S5044 + make_float3 (_S5086.x, _S5086.y, _S5083);
    float2  _S5099 = _S4823 * _S5088 + make_float2 (_S4800 * (v_78 * _S5088.y) + _S4825 * _S5091 + 2.0f * (u_78 * _S5091) + _S4797 * (v_78 * _S5088.x) + _S5097 + _S5097, _S4827 * _S5090 + 2.0f * (v_78 * _S5090) + _S4826 * _S5088.y + _S4824 * _S5088.x + _S5096 + _S5096);
    float2  _S5100 = _S4807 * _S5099;
    float2  _S5101 = _S4818 * _S5099;
    float _S5102 = _S5100.x + _S5100.y;
    if(_S4811)
    {
        float _S5103 = _S5102 / _S4813;
        float _S5104 = _S4814 * - _S5103;
        float _S5105 = _S4810 * (0.3333333432674408f * - (_S4809 * _S5103));
        k_12 = _S5105 + _S5105;
        _S4812 = _S5104;
        _S4813 = 0.0f;
    }
    else
    {
        float _S5106 = _S5102 / _S4812;
        float _S5107 = _S4810 * - _S5106;
        k_12 = _S4808 * _S5106;
        _S4812 = 0.0f;
        _S4813 = _S5107;
    }
    DiffPair_float_0 _S5108;
    (&_S5108)->primal_0 = _S4808;
    (&_S5108)->differential_0 = 0.0f;
    DiffPair_float_0 _S5109;
    (&_S5109)->primal_0 = _S4809;
    (&_S5109)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5108, &_S5109, k_12);
    float _S5110 = _S5109.differential_0 + _S4812;
    float _S5111 = _S5108.differential_0 + _S4813;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5112;
    (&_S5112)->primal_0 = _S4807;
    (&_S5112)->differential_0 = _S5017;
    s_bwd_length_impl_0(&_S5112, _S5111);
    float2  _S5113 = _S5112.differential_0 + _S5101;
    float _S5114 = fx_35 * (_S5066.x + _S5065);
    float2  _S5115 = make_float2 (_S5114, fy_35 * (_S5066.y + _S5068)) + make_float2 ((*dist_coeffs_39)[int(8)] * _S5114, (*dist_coeffs_39)[int(9)] * _S5114);
    float2  _S5116 = _S4792 * _S5115;
    float _S5117 = (*dist_coeffs_39)[int(4)] * _S5115.y;
    float _S5118 = (*dist_coeffs_39)[int(5)] * _S5115.x;
    float _S5119 = _S5116.x + _S5116.y;
    float _S5120 = r2_77 * _S5119;
    float _S5121 = r2_77 * _S5120;
    float _S5122 = (*dist_coeffs_39)[int(7)] * _S5115.y + _S5117 + (*dist_coeffs_39)[int(6)] * _S5115.x + _S5118 + _S4795 * _S5119 + _S4794 * _S5120 + _S4793 * _S5121 + (*dist_coeffs_39)[int(3)] * (r2_77 * _S5121);
    float _S5123 = v_77 * _S5122;
    float _S5124 = u_77 * _S5122;
    float3  _S5125 = _S5042 + make_float3 (_S5113.x, _S5113.y, _S5110);
    float2  _S5126 = _S4796 * _S5115 + make_float2 (_S4800 * (v_77 * _S5115.y) + _S4799 * _S5118 + 2.0f * (u_77 * _S5118) + _S4797 * (v_77 * _S5115.x) + _S5124 + _S5124, _S4802 * _S5117 + 2.0f * (v_77 * _S5117) + _S4801 * _S5115.y + _S4798 * _S5115.x + _S5123 + _S5123);
    float2  _S5127 = _S4780 * _S5126;
    float2  _S5128 = _S4791 * _S5126;
    float _S5129 = _S5127.x + _S5127.y;
    if(_S4784)
    {
        float _S5130 = _S5129 / _S4786;
        float _S5131 = _S4787 * - _S5130;
        float _S5132 = _S4783 * (0.3333333432674408f * - (_S4782 * _S5130));
        k_12 = _S5132 + _S5132;
        _S4785 = _S5131;
        _S4786 = 0.0f;
    }
    else
    {
        float _S5133 = _S5129 / _S4785;
        float _S5134 = _S4783 * - _S5133;
        k_12 = _S4781 * _S5133;
        _S4785 = 0.0f;
        _S4786 = _S5134;
    }
    DiffPair_float_0 _S5135;
    (&_S5135)->primal_0 = _S4781;
    (&_S5135)->differential_0 = 0.0f;
    DiffPair_float_0 _S5136;
    (&_S5136)->primal_0 = _S4782;
    (&_S5136)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S5135, &_S5136, k_12);
    float _S5137 = _S5136.differential_0 + _S4785;
    float _S5138 = _S5135.differential_0 + _S4786;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5139;
    (&_S5139)->primal_0 = _S4780;
    (&_S5139)->differential_0 = _S5017;
    s_bwd_length_impl_0(&_S5139, _S5138);
    float2  _S5140 = _S5139.differential_0 + _S5128;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5141;
    (&_S5141)->primal_0 = vert2_c_11;
    (&_S5141)->differential_0 = _S4910;
    s_bwd_length_impl_1(&_S5141, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5142;
    (&_S5142)->primal_0 = vert1_c_11;
    (&_S5142)->differential_0 = _S4910;
    s_bwd_length_impl_1(&_S5142, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5143;
    (&_S5143)->primal_0 = vert0_c_11;
    (&_S5143)->differential_0 = _S4910;
    s_bwd_length_impl_1(&_S5143, 0.0f);
    float3  _S5144 = _S5141.differential_0 + _S5098;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5145;
    (&_S5145)->primal_0 = R_32;
    (&_S5145)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5146;
    (&_S5146)->primal_0 = vert2_8;
    (&_S5146)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S5145, &_S5146, _S5144);
    float3  _S5147 = _S5142.differential_0 + _S5125;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5148;
    (&_S5148)->primal_0 = R_32;
    (&_S5148)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5149;
    (&_S5149)->primal_0 = vert1_8;
    (&_S5149)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S5148, &_S5149, _S5147);
    float3  _S5150 = _S5143.differential_0 + _S5043 + make_float3 (_S5140.x, _S5140.y, _S5137);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5151;
    (&_S5151)->primal_0 = R_32;
    (&_S5151)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5152;
    (&_S5152)->primal_0 = vert0_8;
    (&_S5152)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S5151, &_S5152, _S5150);
    float3  _S5153 = _S5146.differential_0 + (*v_verts_1)[int(2)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5154;
    (&_S5154)->primal_0 = _S4774;
    (&_S5154)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5155;
    (&_S5155)->primal_0 = _S4779;
    (&_S5155)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S5154, &_S5155, _S5153);
    float _S5156 = - _S5155.differential_0.y;
    float _S5157 = _S4778 * _S5155.differential_0.x;
    float _S5158 = - (_S4770 * _S5155.differential_0.x);
    float3  _S5159 = _S5149.differential_0 + (*v_verts_1)[int(1)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5160;
    (&_S5160)->primal_0 = _S4774;
    (&_S5160)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5161;
    (&_S5161)->primal_0 = _S4777;
    (&_S5161)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S5160, &_S5161, _S5159);
    float _S5162 = _S4770 * _S5161.differential_0.x;
    float _S5163 = _S4776 * _S5161.differential_0.x;
    float3  _S5164 = _S5152.differential_0 + (*v_verts_1)[int(0)];
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5165;
    (&_S5165)->primal_0 = _S4774;
    (&_S5165)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5166;
    (&_S5166)->primal_0 = _S4775;
    (&_S5166)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S5165, &_S5166, _S5164);
    Matrix<float, 3, 3>  _S5167 = transpose_0(_S5154.differential_0 + _S5160.differential_0 + _S5165.differential_0);
    float _S5168 = 2.0f * - _S5167.rows[int(2)].z;
    float _S5169 = 2.0f * _S5167.rows[int(2)].y;
    float _S5170 = 2.0f * _S5167.rows[int(2)].x;
    float _S5171 = 2.0f * _S5167.rows[int(1)].z;
    float _S5172 = 2.0f * - _S5167.rows[int(1)].y;
    float _S5173 = 2.0f * _S5167.rows[int(1)].x;
    float _S5174 = 2.0f * _S5167.rows[int(0)].z;
    float _S5175 = 2.0f * _S5167.rows[int(0)].y;
    float _S5176 = 2.0f * - _S5167.rows[int(0)].x;
    float _S5177 = - _S5173 + _S5175;
    float _S5178 = _S5170 + - _S5174;
    float _S5179 = - _S5169 + _S5171;
    float _S5180 = _S5169 + _S5171;
    float _S5181 = _S5170 + _S5174;
    float _S5182 = _S5173 + _S5175;
    float _S5183 = quat_32.w * (_S5172 + _S5176);
    float _S5184 = quat_32.z * (_S5168 + _S5176);
    float _S5185 = quat_32.y * (_S5168 + _S5172);
    float _S5186 = quat_32.x * _S5177 + quat_32.z * _S5180 + quat_32.y * _S5181 + _S5183 + _S5183;
    float _S5187 = quat_32.x * _S5178 + quat_32.w * _S5180 + quat_32.y * _S5182 + _S5184 + _S5184;
    float _S5188 = quat_32.x * _S5179 + quat_32.w * _S5181 + quat_32.z * _S5182 + _S5185 + _S5185;
    float _S5189 = quat_32.w * _S5177 + quat_32.z * _S5178 + quat_32.y * _S5179;
    float _S5190 = _S5158 + _S5162;
    float _S5191 = 0.5f * - _S5190;
    float _S5192 = _S5156 + _S5161.differential_0.y;
    DiffPair_float_0 _S5193;
    (&_S5193)->primal_0 = _S4771;
    (&_S5193)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5193, _S5192);
    float _S5194 = _S5191 + _S5193.differential_0;
    float _S5195 = _S5157 + _S5163 + _S5166.differential_0.x;
    DiffPair_float_0 _S5196;
    (&_S5196)->primal_0 = _S4769;
    (&_S5196)->differential_0 = 0.0f;
    s_bwd_prop_exp_0(&_S5196, _S5195);
    float _S5197 = _S5191 + _S5196.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5198;
    (&_S5198)->primal_0 = mean_c_26;
    (&_S5198)->differential_0 = _S4910;
    s_bwd_length_impl_1(&_S5198, 0.0f);
    float3  _S5199 = _S5198.differential_0 + _S4913.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5200;
    (&_S5200)->primal_0 = R_32;
    (&_S5200)->differential_0 = _S4991;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5201;
    (&_S5201)->primal_0 = mean_31;
    (&_S5201)->differential_0 = _S4910;
    s_bwd_prop_mul_1(&_S5200, &_S5201, _S5199);
    float3  _S5202 = _S5144 + _S5147 + _S5150 + _S5199 + _S4994.differential_0;
    Matrix<float, 3, 3>  _S5203 = _S5145.differential_0 + _S5148.differential_0 + _S5151.differential_0 + _S5200.differential_0 + _S4995;
    float3  _S5204 = make_float3 (_S5197, _S5194, _S5190);
    float4  _S5205 = make_float4 (0.0f);
    *&((&_S5205)->w) = _S5186;
    *&((&_S5205)->z) = _S5187;
    *&((&_S5205)->y) = _S5188;
    *&((&_S5205)->x) = _S5189;
    float4  _S5206 = _S5205;
    float3  _S5207 = _S5153 + _S5159 + _S5164 + _S5201.differential_0 + _S4989;
    *v_mean_10 = _S5207;
    *v_quat_9 = _S5206;
    *v_scale_9 = _S5204;
    *v_hardness_5 = _S5071;
    (*v_sh_coeffs_8)[int(0)] = _S5049;
    (*v_sh_coeffs_8)[int(1)] = _S5050;
    (*v_sh_coeffs_8)[int(2)] = _S5051;
    (*v_sh_coeffs_8)[int(3)] = _S5052;
    (*v_sh_coeffs_8)[int(4)] = _S5053;
    (*v_sh_coeffs_8)[int(5)] = _S5054;
    (*v_sh_coeffs_8)[int(6)] = _S5055;
    (*v_sh_coeffs_8)[int(7)] = _S5056;
    (*v_sh_coeffs_8)[int(8)] = _S5057;
    (*v_sh_coeffs_8)[int(9)] = _S5058;
    (*v_sh_coeffs_8)[int(10)] = _S5059;
    (*v_sh_coeffs_8)[int(11)] = _S5060;
    (*v_sh_coeffs_8)[int(12)] = _S5061;
    (*v_sh_coeffs_8)[int(13)] = _S5062;
    (*v_sh_coeffs_8)[int(14)] = _S5063;
    (*v_sh_coeffs_8)[int(15)] = _S5064;
    (*v_ch_coeffs_3)[int(0)] = _S5046;
    (*v_ch_coeffs_3)[int(1)] = _S5047;
    *v_R_10 = _S5203;
    *v_t_9 = _S5202;
    return;
}

inline __device__ bool ray_triangle_intersection_uvt(float3  ray_o_4, float3  ray_d_4, FixedArray<float3 , 3>  * verts_4, float * u_80, float * v_80, float * t_30)
{
    float3  v1v0_0 = (*verts_4)[int(1)] - (*verts_4)[int(0)];
    float3  v2v0_0 = (*verts_4)[int(2)] - (*verts_4)[int(0)];
    float3  rov0_0 = ray_o_4 - (*verts_4)[int(0)];
    float3  n_0 = cross_0(v1v0_0, v2v0_0);
    float3  q_2 = cross_0(rov0_0, ray_d_4);
    float d_28 = 1.0f / dot_0(ray_d_4, n_0);
    *u_80 = d_28 * dot_0(- q_2, v2v0_0);
    *v_80 = d_28 * dot_0(q_2, v1v0_0);
    *t_30 = d_28 * dot_0(- n_0, rov0_0);
    bool _S5208;
    if((*u_80) >= 0.0f)
    {
        _S5208 = (*v_80) >= 0.0f;
    }
    else
    {
        _S5208 = false;
    }
    if(_S5208)
    {
        _S5208 = (*u_80 + *v_80) <= 1.0f;
    }
    else
    {
        _S5208 = false;
    }
    if(_S5208)
    {
        _S5208 = (*t_30) >= 0.0f;
    }
    else
    {
        _S5208 = false;
    }
    return _S5208;
}

inline __device__ float evaluate_alpha_opaque_triangle(FixedArray<float3 , 3>  * verts_5, float2  hardness_16, float3  ray_o_5, float3  ray_d_5)
{
    float3  v1v0_1 = (*verts_5)[int(1)] - (*verts_5)[int(0)];
    float3  v2v0_1 = (*verts_5)[int(2)] - (*verts_5)[int(0)];
    float3  rov0_1 = ray_o_5 - (*verts_5)[int(0)];
    float3  n_1 = cross_0(v1v0_1, v2v0_1);
    float3  q_3 = cross_0(rov0_1, ray_d_5);
    float d_29 = 1.0f / dot_0(ray_d_5, n_1);
    float u_81 = d_29 * dot_0(- q_3, v2v0_1);
    float v_81 = d_29 * dot_0(q_3, v1v0_1);
    float t_31 = d_29 * dot_0(- n_1, rov0_1);
    bool _S5209;
    if(u_81 >= 0.0f)
    {
        _S5209 = v_81 >= 0.0f;
    }
    else
    {
        _S5209 = false;
    }
    if(_S5209)
    {
        _S5209 = (u_81 + v_81) <= 1.0f;
    }
    else
    {
        _S5209 = false;
    }
    if(_S5209)
    {
        _S5209 = t_31 >= 0.0f;
    }
    else
    {
        _S5209 = false;
    }
    if(!_S5209)
    {
        return 0.0f;
    }
    float opac_0 = (F32_min(((F32_min((u_81), (v_81)))), ((F32_sqrt((0.5f))) * (1.0f - u_81 - v_81)))) * (2.0f + (F32_sqrt((2.0f))));
    float w_0 = 1.0f - (F32_pow((1.0f - opac_0), (1.0f / (1.0f - clamp_0(hardness_16.y, 0.0f, 0.99989998340606689f)))));
    float o_0 = hardness_16.x;
    float _S5210;
    if(opac_0 < 0.0f)
    {
        _S5210 = 0.0f;
    }
    else
    {
        _S5210 = (F32_min((o_0 * w_0), (0.99500000476837158f)));
    }
    return _S5210;
}

inline __device__ void s_bwd_prop_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_0, DiffPair_vectorx3Cfloatx2C2x3E_0 * dphardness_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_2, float _s_dOut_13)
{
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S5211 = *dphardness_2;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5212 = *dpray_d_2;
    float3  v1v0_2 = dpverts_0->primal_0[int(1)] - dpverts_0->primal_0[int(0)];
    float3  v2v0_2 = dpverts_0->primal_0[int(2)] - dpverts_0->primal_0[int(0)];
    float3  rov0_2 = (*dpray_o_2).primal_0 - dpverts_0->primal_0[int(0)];
    float3  _S5213 = s_primal_ctx_cross_0(v1v0_2, v2v0_2);
    float3  _S5214 = s_primal_ctx_cross_0(rov0_2, (*dpray_d_2).primal_0);
    float _S5215 = s_primal_ctx_dot_0((*dpray_d_2).primal_0, _S5213);
    float d_30 = 1.0f / _S5215;
    float _S5216 = _S5215 * _S5215;
    float3  _S5217 = - _S5214;
    float _S5218 = s_primal_ctx_dot_0(_S5217, v2v0_2);
    float u_82 = d_30 * _S5218;
    float _S5219 = s_primal_ctx_dot_0(_S5214, v1v0_2);
    float v_82 = d_30 * _S5219;
    float3  _S5220 = - _S5213;
    float t_32 = d_30 * s_primal_ctx_dot_0(_S5220, rov0_2);
    bool _S5221;
    if(u_82 >= 0.0f)
    {
        _S5221 = v_82 >= 0.0f;
    }
    else
    {
        _S5221 = false;
    }
    if(_S5221)
    {
        _S5221 = (u_82 + v_82) <= 1.0f;
    }
    else
    {
        _S5221 = false;
    }
    if(_S5221)
    {
        _S5221 = t_32 >= 0.0f;
    }
    else
    {
        _S5221 = false;
    }
    bool _S5222 = !!_S5221;
    float _S5223;
    float _S5224;
    float _S5225;
    float _S5226;
    float _S5227;
    float _S5228;
    float _S5229;
    float _S5230;
    float _S5231;
    float _S5232;
    float _S5233;
    if(_S5222)
    {
        float _S5234 = s_primal_ctx_min_0(u_82, v_82);
        float _S5235 = s_primal_ctx_sqrt_0(0.5f);
        float _S5236 = _S5235 * (1.0f - u_82 - v_82);
        float _S5237 = 2.0f + s_primal_ctx_sqrt_0(2.0f);
        float opac_1 = s_primal_ctx_min_0(_S5234, _S5236) * _S5237;
        float _S5238 = _S5211.primal_0.y;
        float _S5239 = 1.0f - opac_1;
        float _S5240 = 1.0f - s_primal_ctx_clamp_0(_S5238, 0.0f, 0.99989998340606689f);
        float _S5241 = 1.0f / _S5240;
        float _S5242 = _S5240 * _S5240;
        float w_1 = 1.0f - s_primal_ctx_pow_0(_S5239, _S5241);
        float o_1 = _S5211.primal_0.x;
        bool _S5243 = opac_1 < 0.0f;
        if(_S5243)
        {
            _S5223 = 0.0f;
        }
        else
        {
            _S5223 = o_1 * w_1;
        }
        _S5221 = _S5243;
        _S5224 = o_1;
        _S5225 = w_1;
        _S5226 = _S5239;
        _S5227 = _S5241;
        _S5228 = _S5242;
        _S5229 = _S5238;
        _S5230 = _S5237;
        _S5231 = _S5234;
        _S5232 = _S5236;
        _S5233 = _S5235;
    }
    else
    {
        _S5221 = false;
        _S5223 = 0.0f;
        _S5224 = 0.0f;
        _S5225 = 0.0f;
        _S5226 = 0.0f;
        _S5227 = 0.0f;
        _S5228 = 0.0f;
        _S5229 = 0.0f;
        _S5230 = 0.0f;
        _S5231 = 0.0f;
        _S5232 = 0.0f;
        _S5233 = 0.0f;
    }
    float2  _S5244 = make_float2 (0.0f);
    float2  _S5245;
    if(_S5222)
    {
        if(_S5221)
        {
            _S5223 = 0.0f;
            _S5224 = 0.0f;
        }
        else
        {
            DiffPair_float_0 _S5246;
            (&_S5246)->primal_0 = _S5223;
            (&_S5246)->differential_0 = 0.0f;
            DiffPair_float_0 _S5247;
            (&_S5247)->primal_0 = 0.99500000476837158f;
            (&_S5247)->differential_0 = 0.0f;
            s_bwd_prop_min_0(&_S5246, &_S5247, _s_dOut_13);
            float _S5248 = _S5224 * _S5246.differential_0;
            _S5223 = _S5225 * _S5246.differential_0;
            _S5224 = _S5248;
        }
        float _S5249 = - _S5224;
        DiffPair_float_0 _S5250;
        (&_S5250)->primal_0 = _S5226;
        (&_S5250)->differential_0 = 0.0f;
        DiffPair_float_0 _S5251;
        (&_S5251)->primal_0 = _S5227;
        (&_S5251)->differential_0 = 0.0f;
        s_bwd_prop_pow_0(&_S5250, &_S5251, _S5249);
        float _S5252 = - - (_S5251.differential_0 / _S5228);
        float s_diff_opac_T_0 = - _S5250.differential_0;
        DiffPair_float_0 _S5253;
        (&_S5253)->primal_0 = _S5229;
        (&_S5253)->differential_0 = 0.0f;
        DiffPair_float_0 _S5254;
        (&_S5254)->primal_0 = 0.0f;
        (&_S5254)->differential_0 = 0.0f;
        DiffPair_float_0 _S5255;
        (&_S5255)->primal_0 = 0.99989998340606689f;
        (&_S5255)->differential_0 = 0.0f;
        s_bwd_prop_clamp_0(&_S5253, &_S5254, &_S5255, _S5252);
        float _S5256 = _S5230 * s_diff_opac_T_0;
        DiffPair_float_0 _S5257;
        (&_S5257)->primal_0 = _S5231;
        (&_S5257)->differential_0 = 0.0f;
        DiffPair_float_0 _S5258;
        (&_S5258)->primal_0 = _S5232;
        (&_S5258)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5257, &_S5258, _S5256);
        float _S5259 = - (_S5233 * _S5258.differential_0);
        DiffPair_float_0 _S5260;
        (&_S5260)->primal_0 = u_82;
        (&_S5260)->differential_0 = 0.0f;
        DiffPair_float_0 _S5261;
        (&_S5261)->primal_0 = v_82;
        (&_S5261)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S5260, &_S5261, _S5257.differential_0);
        float2  _S5262 = make_float2 (_S5223, _S5253.differential_0);
        float _S5263 = _S5259 + _S5261.differential_0;
        _S5223 = _S5259 + _S5260.differential_0;
        _S5224 = _S5263;
        _S5245 = _S5262;
    }
    else
    {
        _S5223 = 0.0f;
        _S5224 = 0.0f;
        _S5245 = _S5244;
    }
    float3  _S5264 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5265;
    (&_S5265)->primal_0 = _S5220;
    (&_S5265)->differential_0 = _S5264;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5266;
    (&_S5266)->primal_0 = rov0_2;
    (&_S5266)->differential_0 = _S5264;
    s_bwd_prop_dot_0(&_S5265, &_S5266, 0.0f);
    float3  _S5267 = - _S5265.differential_0;
    float _S5268 = d_30 * _S5224;
    float _S5269 = _S5219 * _S5224;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5270;
    (&_S5270)->primal_0 = _S5214;
    (&_S5270)->differential_0 = _S5264;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5271;
    (&_S5271)->primal_0 = v1v0_2;
    (&_S5271)->differential_0 = _S5264;
    s_bwd_prop_dot_0(&_S5270, &_S5271, _S5268);
    float _S5272 = d_30 * _S5223;
    float _S5273 = _S5218 * _S5223;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5274;
    (&_S5274)->primal_0 = _S5217;
    (&_S5274)->differential_0 = _S5264;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5275;
    (&_S5275)->primal_0 = v2v0_2;
    (&_S5275)->differential_0 = _S5264;
    s_bwd_prop_dot_0(&_S5274, &_S5275, _S5272);
    float3  _S5276 = - _S5274.differential_0;
    float _S5277 = - ((_S5269 + _S5273) / _S5216);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5278;
    (&_S5278)->primal_0 = _S5212.primal_0;
    (&_S5278)->differential_0 = _S5264;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5279;
    (&_S5279)->primal_0 = _S5213;
    (&_S5279)->differential_0 = _S5264;
    s_bwd_prop_dot_0(&_S5278, &_S5279, _S5277);
    float3  _S5280 = _S5270.differential_0 + _S5276;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5281;
    (&_S5281)->primal_0 = rov0_2;
    (&_S5281)->differential_0 = _S5264;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5282;
    (&_S5282)->primal_0 = _S5212.primal_0;
    (&_S5282)->differential_0 = _S5264;
    s_bwd_prop_cross_0(&_S5281, &_S5282, _S5280);
    float3  _S5283 = _S5267 + _S5279.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5284;
    (&_S5284)->primal_0 = v1v0_2;
    (&_S5284)->differential_0 = _S5264;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5285;
    (&_S5285)->primal_0 = v2v0_2;
    (&_S5285)->differential_0 = _S5264;
    s_bwd_prop_cross_0(&_S5284, &_S5285, _S5283);
    float3  _S5286 = _S5266.differential_0 + _S5281.differential_0;
    float3  _S5287 = _S5275.differential_0 + _S5285.differential_0;
    float3  _S5288 = _S5271.differential_0 + _S5284.differential_0;
    float3  _S5289 = - _S5286 + - _S5287 + - _S5288;
    float3  _S5290 = _S5278.differential_0 + _S5282.differential_0;
    dpray_d_2->primal_0 = (*dpray_d_2).primal_0;
    dpray_d_2->differential_0 = _S5290;
    dpray_o_2->primal_0 = (*dpray_o_2).primal_0;
    dpray_o_2->differential_0 = _S5286;
    dphardness_2->primal_0 = (*dphardness_2).primal_0;
    dphardness_2->differential_0 = _S5245;
    FixedArray<float3 , 3>  _S5291;
    _S5291[int(0)] = _S5264;
    _S5291[int(1)] = _S5264;
    _S5291[int(2)] = _S5264;
    _S5291[int(2)] = _S5287;
    _S5291[int(0)] = _S5289;
    _S5291[int(1)] = _S5288;
    dpverts_0->primal_0 = dpverts_0->primal_0;
    dpverts_0->differential_0 = _S5291;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_opaque_triangle_0(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5292, DiffPair_vectorx3Cfloatx2C2x3E_0 * _S5293, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5294, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5295, float _S5296)
{
    s_bwd_prop_evaluate_alpha_opaque_triangle_0(_S5292, _S5293, _S5294, _S5295, _S5296);
    return;
}

inline __device__ void evaluate_alpha_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_6, float2  hardness_17, float3  ray_o_6, float3  ray_d_6, float v_alpha_3, FixedArray<float3 , 3>  * v_verts_2, float2  * v_hardness_6, float3  * v_ray_o_3, float3  * v_ray_d_3)
{
    float3  _S5297 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5298 = { _S5297, _S5297, _S5297 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_0;
    (&dp_verts_0)->primal_0 = *verts_6;
    (&dp_verts_0)->differential_0 = _S5298;
    float2  _S5299 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 dp_hardness_2;
    (&dp_hardness_2)->primal_0 = hardness_17;
    (&dp_hardness_2)->differential_0 = _S5299;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_2;
    (&dp_ray_o_2)->primal_0 = ray_o_6;
    (&dp_ray_o_2)->differential_0 = _S5297;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_2;
    (&dp_ray_d_2)->primal_0 = ray_d_6;
    (&dp_ray_d_2)->differential_0 = _S5297;
    s_bwd_evaluate_alpha_opaque_triangle_0(&dp_verts_0, &dp_hardness_2, &dp_ray_o_2, &dp_ray_d_2, v_alpha_3);
    *v_verts_2 = (&dp_verts_0)->differential_0;
    *v_hardness_6 = dp_hardness_2.differential_0;
    *v_ray_o_3 = dp_ray_o_2.differential_0;
    *v_ray_d_3 = dp_ray_d_2.differential_0;
    return;
}

inline __device__ float evaluate_sorting_depth_opaque_triangle(FixedArray<float3 , 3>  * verts_7, FixedArray<float3 , 3>  * rgbs_4, float3  ray_o_7, float3  ray_d_7)
{
    float3  n_2 = cross_0((*verts_7)[int(1)] - (*verts_7)[int(0)], (*verts_7)[int(2)] - (*verts_7)[int(0)]);
    return 1.0f / dot_0(ray_d_7, n_2) * dot_0(- n_2, ray_o_7 - (*verts_7)[int(0)]);
}

inline __device__ void evaluate_color_opaque_triangle(FixedArray<float3 , 3>  * verts_8, FixedArray<float3 , 3>  * rgbs_5, float3  ray_o_8, float3  ray_d_8, float3  * color_13, float * depth_20)
{
    float3  v1v0_3 = (*verts_8)[int(1)] - (*verts_8)[int(0)];
    float3  v2v0_3 = (*verts_8)[int(2)] - (*verts_8)[int(0)];
    float3  rov0_3 = ray_o_8 - (*verts_8)[int(0)];
    float3  n_3 = cross_0(v1v0_3, v2v0_3);
    float3  q_4 = cross_0(rov0_3, ray_d_8);
    float d_31 = 1.0f / dot_0(ray_d_8, n_3);
    float u_83 = d_31 * dot_0(- q_4, v2v0_3);
    float v_83 = d_31 * dot_0(q_4, v1v0_3);
    *depth_20 = d_31 * dot_0(- n_3, rov0_3);
    *color_13 = (*rgbs_5)[int(0)] * make_float3 (1.0f - u_83 - v_83) + (*rgbs_5)[int(1)] * make_float3 (u_83) + (*rgbs_5)[int(2)] * make_float3 (v_83);
    *depth_20 = (F32_log(((F32_max((*depth_20), (9.999999960041972e-13f))))));
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dpverts_1, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * dprgbs_0, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_3, float3  dpcolor_1, float dpdepth_2)
{
    float3  v1v0_4 = dpverts_1->primal_0[int(1)] - dpverts_1->primal_0[int(0)];
    float3  v2v0_4 = dpverts_1->primal_0[int(2)] - dpverts_1->primal_0[int(0)];
    float3  rov0_4 = (*dpray_o_3).primal_0 - dpverts_1->primal_0[int(0)];
    float3  _S5300 = s_primal_ctx_cross_0(v1v0_4, v2v0_4);
    float3  _S5301 = s_primal_ctx_cross_0(rov0_4, (*dpray_d_3).primal_0);
    float _S5302 = s_primal_ctx_dot_0((*dpray_d_3).primal_0, _S5300);
    float d_32 = 1.0f / _S5302;
    float _S5303 = _S5302 * _S5302;
    float3  _S5304 = - _S5301;
    float _S5305 = s_primal_ctx_dot_0(_S5304, v2v0_4);
    float u_84 = d_32 * _S5305;
    float3  _S5306 = make_float3 (u_84);
    float _S5307 = s_primal_ctx_dot_0(_S5301, v1v0_4);
    float v_84 = d_32 * _S5307;
    float3  _S5308 = make_float3 (v_84);
    float3  _S5309 = - _S5300;
    float _S5310 = s_primal_ctx_dot_0(_S5309, rov0_4);
    float _S5311 = d_32 * _S5310;
    float3  _S5312 = make_float3 (1.0f - u_84 - v_84);
    DiffPair_float_0 _S5313;
    (&_S5313)->primal_0 = s_primal_ctx_max_0(_S5311, 9.999999960041972e-13f);
    (&_S5313)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5313, dpdepth_2);
    DiffPair_float_0 _S5314;
    (&_S5314)->primal_0 = _S5311;
    (&_S5314)->differential_0 = 0.0f;
    DiffPair_float_0 _S5315;
    (&_S5315)->primal_0 = 9.999999960041972e-13f;
    (&_S5315)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5314, &_S5315, _S5313.differential_0);
    float3  _S5316 = dprgbs_0->primal_0[int(2)] * dpcolor_1;
    float3  _S5317 = _S5308 * dpcolor_1;
    float3  _S5318 = dprgbs_0->primal_0[int(1)] * dpcolor_1;
    float3  _S5319 = _S5306 * dpcolor_1;
    float3  _S5320 = dprgbs_0->primal_0[int(0)] * dpcolor_1;
    float3  _S5321 = _S5312 * dpcolor_1;
    float _S5322 = - (_S5320.x + _S5320.y + _S5320.z);
    float _S5323 = d_32 * _S5314.differential_0;
    float _S5324 = _S5310 * _S5314.differential_0;
    float3  _S5325 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5326;
    (&_S5326)->primal_0 = _S5309;
    (&_S5326)->differential_0 = _S5325;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5327;
    (&_S5327)->primal_0 = rov0_4;
    (&_S5327)->differential_0 = _S5325;
    s_bwd_prop_dot_0(&_S5326, &_S5327, _S5323);
    float3  _S5328 = - _S5326.differential_0;
    float _S5329 = _S5322 + _S5316.x + _S5316.y + _S5316.z;
    float _S5330 = d_32 * _S5329;
    float _S5331 = _S5307 * _S5329;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5332;
    (&_S5332)->primal_0 = _S5301;
    (&_S5332)->differential_0 = _S5325;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5333;
    (&_S5333)->primal_0 = v1v0_4;
    (&_S5333)->differential_0 = _S5325;
    s_bwd_prop_dot_0(&_S5332, &_S5333, _S5330);
    float _S5334 = _S5322 + _S5318.x + _S5318.y + _S5318.z;
    float _S5335 = d_32 * _S5334;
    float _S5336 = _S5305 * _S5334;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5337;
    (&_S5337)->primal_0 = _S5304;
    (&_S5337)->differential_0 = _S5325;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5338;
    (&_S5338)->primal_0 = v2v0_4;
    (&_S5338)->differential_0 = _S5325;
    s_bwd_prop_dot_0(&_S5337, &_S5338, _S5335);
    float3  _S5339 = - _S5337.differential_0;
    float _S5340 = - ((_S5324 + _S5331 + _S5336) / _S5303);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5341;
    (&_S5341)->primal_0 = (*dpray_d_3).primal_0;
    (&_S5341)->differential_0 = _S5325;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5342;
    (&_S5342)->primal_0 = _S5300;
    (&_S5342)->differential_0 = _S5325;
    s_bwd_prop_dot_0(&_S5341, &_S5342, _S5340);
    float3  _S5343 = _S5332.differential_0 + _S5339;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5344;
    (&_S5344)->primal_0 = rov0_4;
    (&_S5344)->differential_0 = _S5325;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5345;
    (&_S5345)->primal_0 = (*dpray_d_3).primal_0;
    (&_S5345)->differential_0 = _S5325;
    s_bwd_prop_cross_0(&_S5344, &_S5345, _S5343);
    float3  _S5346 = _S5328 + _S5342.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5347;
    (&_S5347)->primal_0 = v1v0_4;
    (&_S5347)->differential_0 = _S5325;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5348;
    (&_S5348)->primal_0 = v2v0_4;
    (&_S5348)->differential_0 = _S5325;
    s_bwd_prop_cross_0(&_S5347, &_S5348, _S5346);
    float3  _S5349 = _S5327.differential_0 + _S5344.differential_0;
    float3  _S5350 = _S5338.differential_0 + _S5348.differential_0;
    float3  _S5351 = _S5333.differential_0 + _S5347.differential_0;
    float3  _S5352 = - _S5349 + - _S5350 + - _S5351;
    float3  _S5353 = _S5341.differential_0 + _S5345.differential_0;
    dpray_d_3->primal_0 = (*dpray_d_3).primal_0;
    dpray_d_3->differential_0 = _S5353;
    dpray_o_3->primal_0 = (*dpray_o_3).primal_0;
    dpray_o_3->differential_0 = _S5349;
    FixedArray<float3 , 3>  _S5354;
    _S5354[int(0)] = _S5325;
    _S5354[int(1)] = _S5325;
    _S5354[int(2)] = _S5325;
    _S5354[int(2)] = _S5317;
    _S5354[int(1)] = _S5319;
    _S5354[int(0)] = _S5321;
    dprgbs_0->primal_0 = dprgbs_0->primal_0;
    dprgbs_0->differential_0 = _S5354;
    FixedArray<float3 , 3>  _S5355;
    _S5355[int(0)] = _S5325;
    _S5355[int(1)] = _S5325;
    _S5355[int(2)] = _S5325;
    _S5355[int(2)] = _S5350;
    _S5355[int(0)] = _S5352;
    _S5355[int(1)] = _S5351;
    dpverts_1->primal_0 = dpverts_1->primal_0;
    dpverts_1->differential_0 = _S5355;
    return;
}

inline __device__ void s_bwd_evaluate_color_opaque_triangle_1(DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5356, DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 * _S5357, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5358, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S5359, float3  _S5360, float _S5361)
{
    s_bwd_prop_evaluate_color_opaque_triangle_1(_S5356, _S5357, _S5358, _S5359, _S5360, _S5361);
    return;
}

inline __device__ void evaluate_color_opaque_triangle_vjp(FixedArray<float3 , 3>  * verts_9, FixedArray<float3 , 3>  * rgbs_6, float3  ray_o_9, float3  ray_d_9, float3  v_color_1, float v_depth_11, FixedArray<float3 , 3>  * v_verts_3, FixedArray<float3 , 3>  * v_rgbs_2, float3  * v_ray_o_4, float3  * v_ray_d_4)
{
    float3  _S5362 = make_float3 (0.0f);
    FixedArray<float3 , 3>  _S5363 = { _S5362, _S5362, _S5362 };
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_verts_1;
    (&dp_verts_1)->primal_0 = *verts_9;
    (&dp_verts_1)->differential_0 = _S5363;
    DiffPair_arrayx3Cvectorx3Cfloatx2C3x3Ex2C3x3E_0 dp_rgbs_0;
    (&dp_rgbs_0)->primal_0 = *rgbs_6;
    (&dp_rgbs_0)->differential_0 = _S5363;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_3;
    (&dp_ray_o_3)->primal_0 = ray_o_9;
    (&dp_ray_o_3)->differential_0 = _S5362;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_3;
    (&dp_ray_d_3)->primal_0 = ray_d_9;
    (&dp_ray_d_3)->differential_0 = _S5362;
    s_bwd_evaluate_color_opaque_triangle_1(&dp_verts_1, &dp_rgbs_0, &dp_ray_o_3, &dp_ray_d_3, v_color_1, v_depth_11);
    *v_verts_3 = (&dp_verts_1)->differential_0;
    *v_rgbs_2 = (&dp_rgbs_0)->differential_0;
    *v_ray_o_4 = dp_ray_o_3.differential_0;
    *v_ray_d_4 = dp_ray_d_3.differential_0;
    return;
}

inline __device__ void projection_voxel_eval3d_persp(float3  pos_0, float size_0, FixedArray<float, 8>  * densities_0, FixedArray<float3 , 16>  * sh_coeffs_27, Matrix<float, 3, 3>  R_33, float3  t_33, float fx_36, float fy_36, float cx_31, float cy_31, FixedArray<float, 10>  * dist_coeffs_40, uint image_width_27, uint image_height_27, float near_plane_18, float far_plane_18, int4  * aabb_xyxy_18, float * depth_21, float3  * rgbs_7)
{
    float2  * _S5364;
    float2  * _S5365;
    float2  * _S5366;
    float2  * _S5367;
    float2  * _S5368;
    float2  * _S5369;
    float2  * _S5370;
    float2  * _S5371;
    bool _S5372;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_0;
        float3  _S5373 = mul_0(R_33, pos_0) + t_33;
        pos_c_0[int(0)] = _S5373;
        float _S5374 = _S5373.z;
        float _S5375 = (F32_min((far_plane_18), (_S5374)));
        float _S5376 = (F32_max((near_plane_18), (_S5374)));
        float3  _S5377 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 0.0f)) + t_33;
        pos_c_0[int(1)] = _S5377;
        float _S5378 = _S5377.z;
        float _S5379 = (F32_min((_S5375), (_S5378)));
        float _S5380 = (F32_max((_S5376), (_S5378)));
        float3  _S5381 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 0.0f)) + t_33;
        pos_c_0[int(2)] = _S5381;
        float _S5382 = _S5381.z;
        float _S5383 = (F32_min((_S5379), (_S5382)));
        float _S5384 = (F32_max((_S5380), (_S5382)));
        float3  _S5385 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 1.0f, 0.0f)) + t_33;
        pos_c_0[int(3)] = _S5385;
        float _S5386 = _S5385.z;
        float _S5387 = (F32_min((_S5383), (_S5386)));
        float _S5388 = (F32_max((_S5384), (_S5386)));
        float3  _S5389 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 0.0f, 1.0f)) + t_33;
        pos_c_0[int(4)] = _S5389;
        float _S5390 = _S5389.z;
        float _S5391 = (F32_min((_S5387), (_S5390)));
        float _S5392 = (F32_max((_S5388), (_S5390)));
        float3  _S5393 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (1.0f, 0.0f, 1.0f)) + t_33;
        pos_c_0[int(5)] = _S5393;
        float _S5394 = _S5393.z;
        float _S5395 = (F32_min((_S5391), (_S5394)));
        float _S5396 = (F32_max((_S5392), (_S5394)));
        float3  _S5397 = mul_0(R_33, pos_0 + make_float3 (size_0) * make_float3 (0.0f, 1.0f, 1.0f)) + t_33;
        pos_c_0[int(6)] = _S5397;
        float _S5398 = _S5397.z;
        float _S5399 = (F32_min((_S5395), (_S5398)));
        float _S5400 = (F32_max((_S5396), (_S5398)));
        float3  _S5401 = mul_0(R_33, pos_0 + make_float3 (size_0)) + t_33;
        pos_c_0[int(7)] = _S5401;
        float _S5402 = _S5401.z;
        float _S5403 = (F32_min((_S5399), (_S5402)));
        float _S5404 = (F32_max((_S5400), (_S5402)));
        bool _S5405;
        if(_S5403 < near_plane_18)
        {
            _S5405 = true;
        }
        else
        {
            _S5405 = _S5404 > far_plane_18;
        }
        if(_S5405)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  mean_c_27 = mul_0(R_33, pos_0 + make_float3 (0.5f * size_0)) + t_33;
        FixedArray<float2 , 8>  uv_7;
        for(;;)
        {
            float3  _S5406 = pos_c_0[int(0)];
            _S5364 = &uv_7[int(0)];
            for(;;)
            {
                float _S5407 = _S5406.z;
                uv_7[int(0)] = float2 {_S5406.x, _S5406.y} / make_float2 (_S5407);
                if(_S5407 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5408 = camera_distortion_jac_0(uv_7[int(0)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5408)), ((F32_min((_S5408.rows[int(0)].x), (_S5408.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_85 = uv_7[int(0)].x;
                float v_85 = uv_7[int(0)].y;
                float r2_80 = u_85 * u_85 + v_85 * v_85;
                float2  _S5409 = uv_7[int(0)] * make_float2 (1.0f + r2_80 * ((*dist_coeffs_40)[int(0)] + r2_80 * ((*dist_coeffs_40)[int(1)] + r2_80 * ((*dist_coeffs_40)[int(2)] + r2_80 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_85 * v_85 + (*dist_coeffs_40)[int(5)] * (r2_80 + 2.0f * u_85 * u_85) + (*dist_coeffs_40)[int(6)] * r2_80, 2.0f * (*dist_coeffs_40)[int(5)] * u_85 * v_85 + (*dist_coeffs_40)[int(4)] * (r2_80 + 2.0f * v_85 * v_85) + (*dist_coeffs_40)[int(7)] * r2_80);
                float2  _S5410 = _S5409 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5409.x + (*dist_coeffs_40)[int(9)] * _S5409.y, 0.0f);
                uv_7[int(0)] = make_float2 (fx_36 * _S5410.x + cx_31, fy_36 * _S5410.y + cy_31);
                break;
            }
            bool all_valid_16 = true & (!_S5405);
            float3  _S5411 = pos_c_0[int(1)];
            _S5365 = &uv_7[int(1)];
            for(;;)
            {
                float _S5412 = _S5411.z;
                uv_7[int(1)] = float2 {_S5411.x, _S5411.y} / make_float2 (_S5412);
                if(_S5412 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5413 = camera_distortion_jac_0(uv_7[int(1)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5413)), ((F32_min((_S5413.rows[int(0)].x), (_S5413.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_86 = uv_7[int(1)].x;
                float v_86 = uv_7[int(1)].y;
                float r2_81 = u_86 * u_86 + v_86 * v_86;
                float2  _S5414 = uv_7[int(1)] * make_float2 (1.0f + r2_81 * ((*dist_coeffs_40)[int(0)] + r2_81 * ((*dist_coeffs_40)[int(1)] + r2_81 * ((*dist_coeffs_40)[int(2)] + r2_81 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_86 * v_86 + (*dist_coeffs_40)[int(5)] * (r2_81 + 2.0f * u_86 * u_86) + (*dist_coeffs_40)[int(6)] * r2_81, 2.0f * (*dist_coeffs_40)[int(5)] * u_86 * v_86 + (*dist_coeffs_40)[int(4)] * (r2_81 + 2.0f * v_86 * v_86) + (*dist_coeffs_40)[int(7)] * r2_81);
                float2  _S5415 = _S5414 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5414.x + (*dist_coeffs_40)[int(9)] * _S5414.y, 0.0f);
                uv_7[int(1)] = make_float2 (fx_36 * _S5415.x + cx_31, fy_36 * _S5415.y + cy_31);
                break;
            }
            bool all_valid_17 = all_valid_16 & (!_S5405);
            float3  _S5416 = pos_c_0[int(2)];
            _S5366 = &uv_7[int(2)];
            for(;;)
            {
                float _S5417 = _S5416.z;
                uv_7[int(2)] = float2 {_S5416.x, _S5416.y} / make_float2 (_S5417);
                if(_S5417 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5418 = camera_distortion_jac_0(uv_7[int(2)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5418)), ((F32_min((_S5418.rows[int(0)].x), (_S5418.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_87 = uv_7[int(2)].x;
                float v_87 = uv_7[int(2)].y;
                float r2_82 = u_87 * u_87 + v_87 * v_87;
                float2  _S5419 = uv_7[int(2)] * make_float2 (1.0f + r2_82 * ((*dist_coeffs_40)[int(0)] + r2_82 * ((*dist_coeffs_40)[int(1)] + r2_82 * ((*dist_coeffs_40)[int(2)] + r2_82 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_87 * v_87 + (*dist_coeffs_40)[int(5)] * (r2_82 + 2.0f * u_87 * u_87) + (*dist_coeffs_40)[int(6)] * r2_82, 2.0f * (*dist_coeffs_40)[int(5)] * u_87 * v_87 + (*dist_coeffs_40)[int(4)] * (r2_82 + 2.0f * v_87 * v_87) + (*dist_coeffs_40)[int(7)] * r2_82);
                float2  _S5420 = _S5419 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5419.x + (*dist_coeffs_40)[int(9)] * _S5419.y, 0.0f);
                uv_7[int(2)] = make_float2 (fx_36 * _S5420.x + cx_31, fy_36 * _S5420.y + cy_31);
                break;
            }
            bool all_valid_18 = all_valid_17 & (!_S5405);
            float3  _S5421 = pos_c_0[int(3)];
            _S5367 = &uv_7[int(3)];
            for(;;)
            {
                float _S5422 = _S5421.z;
                uv_7[int(3)] = float2 {_S5421.x, _S5421.y} / make_float2 (_S5422);
                if(_S5422 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5423 = camera_distortion_jac_0(uv_7[int(3)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5423)), ((F32_min((_S5423.rows[int(0)].x), (_S5423.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_88 = uv_7[int(3)].x;
                float v_88 = uv_7[int(3)].y;
                float r2_83 = u_88 * u_88 + v_88 * v_88;
                float2  _S5424 = uv_7[int(3)] * make_float2 (1.0f + r2_83 * ((*dist_coeffs_40)[int(0)] + r2_83 * ((*dist_coeffs_40)[int(1)] + r2_83 * ((*dist_coeffs_40)[int(2)] + r2_83 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_88 * v_88 + (*dist_coeffs_40)[int(5)] * (r2_83 + 2.0f * u_88 * u_88) + (*dist_coeffs_40)[int(6)] * r2_83, 2.0f * (*dist_coeffs_40)[int(5)] * u_88 * v_88 + (*dist_coeffs_40)[int(4)] * (r2_83 + 2.0f * v_88 * v_88) + (*dist_coeffs_40)[int(7)] * r2_83);
                float2  _S5425 = _S5424 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5424.x + (*dist_coeffs_40)[int(9)] * _S5424.y, 0.0f);
                uv_7[int(3)] = make_float2 (fx_36 * _S5425.x + cx_31, fy_36 * _S5425.y + cy_31);
                break;
            }
            bool all_valid_19 = all_valid_18 & (!_S5405);
            float3  _S5426 = pos_c_0[int(4)];
            _S5368 = &uv_7[int(4)];
            for(;;)
            {
                float _S5427 = _S5426.z;
                uv_7[int(4)] = float2 {_S5426.x, _S5426.y} / make_float2 (_S5427);
                if(_S5427 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5428 = camera_distortion_jac_0(uv_7[int(4)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5428)), ((F32_min((_S5428.rows[int(0)].x), (_S5428.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_89 = uv_7[int(4)].x;
                float v_89 = uv_7[int(4)].y;
                float r2_84 = u_89 * u_89 + v_89 * v_89;
                float2  _S5429 = uv_7[int(4)] * make_float2 (1.0f + r2_84 * ((*dist_coeffs_40)[int(0)] + r2_84 * ((*dist_coeffs_40)[int(1)] + r2_84 * ((*dist_coeffs_40)[int(2)] + r2_84 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_89 * v_89 + (*dist_coeffs_40)[int(5)] * (r2_84 + 2.0f * u_89 * u_89) + (*dist_coeffs_40)[int(6)] * r2_84, 2.0f * (*dist_coeffs_40)[int(5)] * u_89 * v_89 + (*dist_coeffs_40)[int(4)] * (r2_84 + 2.0f * v_89 * v_89) + (*dist_coeffs_40)[int(7)] * r2_84);
                float2  _S5430 = _S5429 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5429.x + (*dist_coeffs_40)[int(9)] * _S5429.y, 0.0f);
                uv_7[int(4)] = make_float2 (fx_36 * _S5430.x + cx_31, fy_36 * _S5430.y + cy_31);
                break;
            }
            bool all_valid_20 = all_valid_19 & (!_S5405);
            float3  _S5431 = pos_c_0[int(5)];
            _S5369 = &uv_7[int(5)];
            for(;;)
            {
                float _S5432 = _S5431.z;
                uv_7[int(5)] = float2 {_S5431.x, _S5431.y} / make_float2 (_S5432);
                if(_S5432 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5433 = camera_distortion_jac_0(uv_7[int(5)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5433)), ((F32_min((_S5433.rows[int(0)].x), (_S5433.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_90 = uv_7[int(5)].x;
                float v_90 = uv_7[int(5)].y;
                float r2_85 = u_90 * u_90 + v_90 * v_90;
                float2  _S5434 = uv_7[int(5)] * make_float2 (1.0f + r2_85 * ((*dist_coeffs_40)[int(0)] + r2_85 * ((*dist_coeffs_40)[int(1)] + r2_85 * ((*dist_coeffs_40)[int(2)] + r2_85 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_90 * v_90 + (*dist_coeffs_40)[int(5)] * (r2_85 + 2.0f * u_90 * u_90) + (*dist_coeffs_40)[int(6)] * r2_85, 2.0f * (*dist_coeffs_40)[int(5)] * u_90 * v_90 + (*dist_coeffs_40)[int(4)] * (r2_85 + 2.0f * v_90 * v_90) + (*dist_coeffs_40)[int(7)] * r2_85);
                float2  _S5435 = _S5434 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5434.x + (*dist_coeffs_40)[int(9)] * _S5434.y, 0.0f);
                uv_7[int(5)] = make_float2 (fx_36 * _S5435.x + cx_31, fy_36 * _S5435.y + cy_31);
                break;
            }
            bool all_valid_21 = all_valid_20 & (!_S5405);
            float3  _S5436 = pos_c_0[int(6)];
            _S5370 = &uv_7[int(6)];
            for(;;)
            {
                float _S5437 = _S5436.z;
                uv_7[int(6)] = float2 {_S5436.x, _S5436.y} / make_float2 (_S5437);
                if(_S5437 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5438 = camera_distortion_jac_0(uv_7[int(6)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5438)), ((F32_min((_S5438.rows[int(0)].x), (_S5438.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_91 = uv_7[int(6)].x;
                float v_91 = uv_7[int(6)].y;
                float r2_86 = u_91 * u_91 + v_91 * v_91;
                float2  _S5439 = uv_7[int(6)] * make_float2 (1.0f + r2_86 * ((*dist_coeffs_40)[int(0)] + r2_86 * ((*dist_coeffs_40)[int(1)] + r2_86 * ((*dist_coeffs_40)[int(2)] + r2_86 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_91 * v_91 + (*dist_coeffs_40)[int(5)] * (r2_86 + 2.0f * u_91 * u_91) + (*dist_coeffs_40)[int(6)] * r2_86, 2.0f * (*dist_coeffs_40)[int(5)] * u_91 * v_91 + (*dist_coeffs_40)[int(4)] * (r2_86 + 2.0f * v_91 * v_91) + (*dist_coeffs_40)[int(7)] * r2_86);
                float2  _S5440 = _S5439 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5439.x + (*dist_coeffs_40)[int(9)] * _S5439.y, 0.0f);
                uv_7[int(6)] = make_float2 (fx_36 * _S5440.x + cx_31, fy_36 * _S5440.y + cy_31);
                break;
            }
            bool all_valid_22 = all_valid_21 & (!_S5405);
            float3  _S5441 = pos_c_0[int(7)];
            _S5371 = &uv_7[int(7)];
            for(;;)
            {
                float _S5442 = _S5441.z;
                uv_7[int(7)] = float2 {_S5441.x, _S5441.y} / make_float2 (_S5442);
                if(_S5442 < 0.0f)
                {
                    _S5405 = true;
                }
                else
                {
                    Matrix<float, 2, 2>  _S5443 = camera_distortion_jac_0(uv_7[int(7)], dist_coeffs_40);
                    _S5405 = !((F32_min((determinant_0(_S5443)), ((F32_min((_S5443.rows[int(0)].x), (_S5443.rows[int(1)].y)))))) > 0.0f);
                }
                if(_S5405)
                {
                    break;
                }
                float u_92 = uv_7[int(7)].x;
                float v_92 = uv_7[int(7)].y;
                float r2_87 = u_92 * u_92 + v_92 * v_92;
                float2  _S5444 = uv_7[int(7)] * make_float2 (1.0f + r2_87 * ((*dist_coeffs_40)[int(0)] + r2_87 * ((*dist_coeffs_40)[int(1)] + r2_87 * ((*dist_coeffs_40)[int(2)] + r2_87 * (*dist_coeffs_40)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_40)[int(4)] * u_92 * v_92 + (*dist_coeffs_40)[int(5)] * (r2_87 + 2.0f * u_92 * u_92) + (*dist_coeffs_40)[int(6)] * r2_87, 2.0f * (*dist_coeffs_40)[int(5)] * u_92 * v_92 + (*dist_coeffs_40)[int(4)] * (r2_87 + 2.0f * v_92 * v_92) + (*dist_coeffs_40)[int(7)] * r2_87);
                float2  _S5445 = _S5444 + make_float2 ((*dist_coeffs_40)[int(8)] * _S5444.x + (*dist_coeffs_40)[int(9)] * _S5444.y, 0.0f);
                uv_7[int(7)] = make_float2 (fx_36 * _S5445.x + cx_31, fy_36 * _S5445.y + cy_31);
                break;
            }
            _S5372 = all_valid_22 & (!_S5405);
            break;
        }
        if(!_S5372)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*densities_0)[int(0)]), ((*densities_0)[int(1)])))), ((*densities_0)[int(2)])))), ((*densities_0)[int(3)])))), ((*densities_0)[int(4)])))), ((*densities_0)[int(5)])))), ((*densities_0)[int(6)])))), ((*densities_0)[int(7)]))) * size_0 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S5446 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5364).x), ((*_S5365).x)))), ((*_S5366).x)))), ((*_S5367).x)))), ((*_S5368).x)))), ((*_S5369).x)))), ((*_S5370).x)))), ((*_S5371).x)));
        float _S5447 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5364).x), ((*_S5365).x)))), ((*_S5366).x)))), ((*_S5367).x)))), ((*_S5368).x)))), ((*_S5369).x)))), ((*_S5370).x)))), ((*_S5371).x)));
        float _S5448 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5364).y), ((*_S5365).y)))), ((*_S5366).y)))), ((*_S5367).y)))), ((*_S5368).y)))), ((*_S5369).y)))), ((*_S5370).y)))), ((*_S5371).y)));
        float _S5449 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5364).y), ((*_S5365).y)))), ((*_S5366).y)))), ((*_S5367).y)))), ((*_S5368).y)))), ((*_S5369).y)))), ((*_S5370).y)))), ((*_S5371).y)));
        if(_S5446 <= 0.0f)
        {
            _S5405 = true;
        }
        else
        {
            _S5405 = _S5447 >= float(image_width_27);
        }
        if(_S5405)
        {
            _S5405 = true;
        }
        else
        {
            _S5405 = _S5448 <= 0.0f;
        }
        if(_S5405)
        {
            _S5405 = true;
        }
        else
        {
            _S5405 = _S5449 >= float(image_height_27);
        }
        if(_S5405)
        {
            _S5405 = true;
        }
        else
        {
            if(_S5403 <= 0.0f)
            {
                if(_S5447 <= 0.0f)
                {
                    _S5405 = _S5446 >= float(image_width_27);
                }
                else
                {
                    _S5405 = false;
                }
                if(_S5405)
                {
                    _S5405 = true;
                }
                else
                {
                    if(_S5449 <= 0.0f)
                    {
                        _S5405 = _S5448 >= float(image_width_27);
                    }
                    else
                    {
                        _S5405 = false;
                    }
                }
            }
            else
            {
                _S5405 = false;
            }
        }
        if(_S5405)
        {
            *aabb_xyxy_18 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_18 = make_int4 (int((F32_floor((_S5447)))), int((F32_floor((_S5449)))), int((F32_ceil((_S5446)))), int((F32_ceil((_S5448)))));
        *depth_21 = 0.5f * (F32_log((dot_0(mean_c_27, mean_c_27) + 9.99999997475242708e-07f)));
        float3  _S5450 = mean_c_27 - - mul_0(transpose_0(R_33), t_33);
        float3  _S5451 = make_float3 (0.282094806432724f) * (*sh_coeffs_27)[int(0)];
        *rgbs_7 = _S5451;
        float _S5452 = _S5450.x;
        float _S5453 = _S5450.y;
        float _S5454 = _S5450.z;
        float norm_18 = (F32_sqrt((_S5452 * _S5452 + _S5453 * _S5453 + _S5454 * _S5454)));
        float x_63 = _S5452 / norm_18;
        float y_30 = _S5453 / norm_18;
        float z_27 = _S5454 / norm_18;
        float3  _S5455 = _S5451 + make_float3 (0.48860251903533936f) * (make_float3 (- y_30) * (*sh_coeffs_27)[int(1)] + make_float3 (z_27) * (*sh_coeffs_27)[int(2)] - make_float3 (x_63) * (*sh_coeffs_27)[int(3)]);
        *rgbs_7 = _S5455;
        float z2_60 = z_27 * z_27;
        float fTmp0B_27 = -1.09254848957061768f * z_27;
        float fC1_27 = x_63 * x_63 - y_30 * y_30;
        float fS1_27 = 2.0f * x_63 * y_30;
        float3  _S5456 = _S5455 + (make_float3 (0.54627424478530884f * fS1_27) * (*sh_coeffs_27)[int(4)] + make_float3 (fTmp0B_27 * y_30) * (*sh_coeffs_27)[int(5)] + make_float3 (0.94617468118667603f * z2_60 - 0.31539157032966614f) * (*sh_coeffs_27)[int(6)] + make_float3 (fTmp0B_27 * x_63) * (*sh_coeffs_27)[int(7)] + make_float3 (0.54627424478530884f * fC1_27) * (*sh_coeffs_27)[int(8)]);
        *rgbs_7 = _S5456;
        float fTmp0C_27 = -2.28522896766662598f * z2_60 + 0.4570457935333252f;
        float fTmp1B_27 = 1.44530570507049561f * z_27;
        *rgbs_7 = max_0(_S5456 + (make_float3 (-0.59004360437393188f * (x_63 * fS1_27 + y_30 * fC1_27)) * (*sh_coeffs_27)[int(9)] + make_float3 (fTmp1B_27 * fS1_27) * (*sh_coeffs_27)[int(10)] + make_float3 (fTmp0C_27 * y_30) * (*sh_coeffs_27)[int(11)] + make_float3 (z_27 * (1.86588168144226074f * z2_60 - 1.11952900886535645f)) * (*sh_coeffs_27)[int(12)] + make_float3 (fTmp0C_27 * x_63) * (*sh_coeffs_27)[int(13)] + make_float3 (fTmp1B_27 * fC1_27) * (*sh_coeffs_27)[int(14)] + make_float3 (-0.59004360437393188f * (x_63 * fC1_27 - y_30 * fS1_27)) * (*sh_coeffs_27)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye(float3  pos_1, float size_1, FixedArray<float, 8>  * densities_1, FixedArray<float3 , 16>  * sh_coeffs_28, Matrix<float, 3, 3>  R_34, float3  t_34, float fx_37, float fy_37, float cx_32, float cy_32, FixedArray<float, 10>  * dist_coeffs_41, uint image_width_28, uint image_height_28, float near_plane_19, float far_plane_19, int4  * aabb_xyxy_19, float * depth_22, float3  * rgbs_8)
{
    float2  * _S5457;
    bool _S5458;
    float2  * _S5459;
    bool _S5460;
    float2  * _S5461;
    bool _S5462;
    float2  * _S5463;
    bool _S5464;
    float2  * _S5465;
    bool _S5466;
    float2  * _S5467;
    bool _S5468;
    float2  * _S5469;
    bool _S5470;
    float2  * _S5471;
    bool _S5472;
    bool _S5473;
    for(;;)
    {
        FixedArray<float3 , 8>  pos_c_1;
        float3  _S5474 = mul_0(R_34, pos_1) + t_34;
        pos_c_1[int(0)] = _S5474;
        float _S5475 = length_1(_S5474);
        float _S5476 = (F32_min((far_plane_19), (_S5475)));
        float _S5477 = (F32_max((near_plane_19), (_S5475)));
        float3  _S5478 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 0.0f)) + t_34;
        pos_c_1[int(1)] = _S5478;
        float _S5479 = length_1(_S5478);
        float _S5480 = (F32_min((_S5476), (_S5479)));
        float _S5481 = (F32_max((_S5477), (_S5479)));
        float3  _S5482 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 0.0f)) + t_34;
        pos_c_1[int(2)] = _S5482;
        float _S5483 = length_1(_S5482);
        float _S5484 = (F32_min((_S5480), (_S5483)));
        float _S5485 = (F32_max((_S5481), (_S5483)));
        float3  _S5486 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 1.0f, 0.0f)) + t_34;
        pos_c_1[int(3)] = _S5486;
        float _S5487 = length_1(_S5486);
        float _S5488 = (F32_min((_S5484), (_S5487)));
        float _S5489 = (F32_max((_S5485), (_S5487)));
        float3  _S5490 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 0.0f, 1.0f)) + t_34;
        pos_c_1[int(4)] = _S5490;
        float _S5491 = length_1(_S5490);
        float _S5492 = (F32_min((_S5488), (_S5491)));
        float _S5493 = (F32_max((_S5489), (_S5491)));
        float3  _S5494 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (1.0f, 0.0f, 1.0f)) + t_34;
        pos_c_1[int(5)] = _S5494;
        float _S5495 = length_1(_S5494);
        float _S5496 = (F32_min((_S5492), (_S5495)));
        float _S5497 = (F32_max((_S5493), (_S5495)));
        float3  _S5498 = mul_0(R_34, pos_1 + make_float3 (size_1) * make_float3 (0.0f, 1.0f, 1.0f)) + t_34;
        pos_c_1[int(6)] = _S5498;
        float _S5499 = length_1(_S5498);
        float _S5500 = (F32_min((_S5496), (_S5499)));
        float _S5501 = (F32_max((_S5497), (_S5499)));
        float3  _S5502 = mul_0(R_34, pos_1 + make_float3 (size_1)) + t_34;
        pos_c_1[int(7)] = _S5502;
        float _S5503 = length_1(_S5502);
        float _S5504 = (F32_min((_S5500), (_S5503)));
        float _S5505 = (F32_max((_S5501), (_S5503)));
        bool _S5506;
        if(_S5504 < near_plane_19)
        {
            _S5506 = true;
        }
        else
        {
            _S5506 = _S5505 > far_plane_19;
        }
        if(_S5506)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float3  mean_c_28 = mul_0(R_34, pos_1 + make_float3 (0.5f * size_1)) + t_34;
        FixedArray<float2 , 8>  uv_8;
        for(;;)
        {
            float k_13;
            float3  _S5507 = pos_c_1[int(0)];
            _S5457 = &uv_8[int(0)];
            for(;;)
            {
                float2  _S5508 = float2 {_S5507.x, _S5507.y};
                float r_36 = length_0(_S5508);
                float _S5509 = _S5507.z;
                float theta_32 = (F32_atan2((r_36), (_S5509)));
                if(theta_32 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_32 * theta_32 / 3.0f) / _S5509;
                }
                else
                {
                    k_13 = theta_32 / r_36;
                }
                float2  _S5510 = _S5508 * make_float2 (k_13);
                uv_8[int(0)] = _S5510;
                Matrix<float, 2, 2>  _S5511 = camera_distortion_jac_0(_S5510, dist_coeffs_41);
                bool _S5512 = !((F32_min((determinant_0(_S5511)), ((F32_min((_S5511.rows[int(0)].x), (_S5511.rows[int(1)].y)))))) > 0.0f);
                _S5458 = _S5512;
                if(_S5512)
                {
                    break;
                }
                float u_93 = uv_8[int(0)].x;
                float v_93 = uv_8[int(0)].y;
                float r2_88 = u_93 * u_93 + v_93 * v_93;
                float2  _S5513 = uv_8[int(0)] * make_float2 (1.0f + r2_88 * ((*dist_coeffs_41)[int(0)] + r2_88 * ((*dist_coeffs_41)[int(1)] + r2_88 * ((*dist_coeffs_41)[int(2)] + r2_88 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_93 * v_93 + (*dist_coeffs_41)[int(5)] * (r2_88 + 2.0f * u_93 * u_93) + (*dist_coeffs_41)[int(6)] * r2_88, 2.0f * (*dist_coeffs_41)[int(5)] * u_93 * v_93 + (*dist_coeffs_41)[int(4)] * (r2_88 + 2.0f * v_93 * v_93) + (*dist_coeffs_41)[int(7)] * r2_88);
                float2  _S5514 = _S5513 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5513.x + (*dist_coeffs_41)[int(9)] * _S5513.y, 0.0f);
                uv_8[int(0)] = make_float2 (fx_37 * _S5514.x + cx_32, fy_37 * _S5514.y + cy_32);
                break;
            }
            bool all_valid_23 = true & (!_S5458);
            float3  _S5515 = pos_c_1[int(1)];
            _S5459 = &uv_8[int(1)];
            for(;;)
            {
                float2  _S5516 = float2 {_S5515.x, _S5515.y};
                float r_37 = length_0(_S5516);
                float _S5517 = _S5515.z;
                float theta_33 = (F32_atan2((r_37), (_S5517)));
                if(theta_33 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_33 * theta_33 / 3.0f) / _S5517;
                }
                else
                {
                    k_13 = theta_33 / r_37;
                }
                float2  _S5518 = _S5516 * make_float2 (k_13);
                uv_8[int(1)] = _S5518;
                Matrix<float, 2, 2>  _S5519 = camera_distortion_jac_0(_S5518, dist_coeffs_41);
                bool _S5520 = !((F32_min((determinant_0(_S5519)), ((F32_min((_S5519.rows[int(0)].x), (_S5519.rows[int(1)].y)))))) > 0.0f);
                _S5460 = _S5520;
                if(_S5520)
                {
                    break;
                }
                float u_94 = uv_8[int(1)].x;
                float v_94 = uv_8[int(1)].y;
                float r2_89 = u_94 * u_94 + v_94 * v_94;
                float2  _S5521 = uv_8[int(1)] * make_float2 (1.0f + r2_89 * ((*dist_coeffs_41)[int(0)] + r2_89 * ((*dist_coeffs_41)[int(1)] + r2_89 * ((*dist_coeffs_41)[int(2)] + r2_89 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_94 * v_94 + (*dist_coeffs_41)[int(5)] * (r2_89 + 2.0f * u_94 * u_94) + (*dist_coeffs_41)[int(6)] * r2_89, 2.0f * (*dist_coeffs_41)[int(5)] * u_94 * v_94 + (*dist_coeffs_41)[int(4)] * (r2_89 + 2.0f * v_94 * v_94) + (*dist_coeffs_41)[int(7)] * r2_89);
                float2  _S5522 = _S5521 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5521.x + (*dist_coeffs_41)[int(9)] * _S5521.y, 0.0f);
                uv_8[int(1)] = make_float2 (fx_37 * _S5522.x + cx_32, fy_37 * _S5522.y + cy_32);
                break;
            }
            bool all_valid_24 = all_valid_23 & (!_S5460);
            float3  _S5523 = pos_c_1[int(2)];
            _S5461 = &uv_8[int(2)];
            for(;;)
            {
                float2  _S5524 = float2 {_S5523.x, _S5523.y};
                float r_38 = length_0(_S5524);
                float _S5525 = _S5523.z;
                float theta_34 = (F32_atan2((r_38), (_S5525)));
                if(theta_34 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_34 * theta_34 / 3.0f) / _S5525;
                }
                else
                {
                    k_13 = theta_34 / r_38;
                }
                float2  _S5526 = _S5524 * make_float2 (k_13);
                uv_8[int(2)] = _S5526;
                Matrix<float, 2, 2>  _S5527 = camera_distortion_jac_0(_S5526, dist_coeffs_41);
                bool _S5528 = !((F32_min((determinant_0(_S5527)), ((F32_min((_S5527.rows[int(0)].x), (_S5527.rows[int(1)].y)))))) > 0.0f);
                _S5462 = _S5528;
                if(_S5528)
                {
                    break;
                }
                float u_95 = uv_8[int(2)].x;
                float v_95 = uv_8[int(2)].y;
                float r2_90 = u_95 * u_95 + v_95 * v_95;
                float2  _S5529 = uv_8[int(2)] * make_float2 (1.0f + r2_90 * ((*dist_coeffs_41)[int(0)] + r2_90 * ((*dist_coeffs_41)[int(1)] + r2_90 * ((*dist_coeffs_41)[int(2)] + r2_90 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_95 * v_95 + (*dist_coeffs_41)[int(5)] * (r2_90 + 2.0f * u_95 * u_95) + (*dist_coeffs_41)[int(6)] * r2_90, 2.0f * (*dist_coeffs_41)[int(5)] * u_95 * v_95 + (*dist_coeffs_41)[int(4)] * (r2_90 + 2.0f * v_95 * v_95) + (*dist_coeffs_41)[int(7)] * r2_90);
                float2  _S5530 = _S5529 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5529.x + (*dist_coeffs_41)[int(9)] * _S5529.y, 0.0f);
                uv_8[int(2)] = make_float2 (fx_37 * _S5530.x + cx_32, fy_37 * _S5530.y + cy_32);
                break;
            }
            bool all_valid_25 = all_valid_24 & (!_S5462);
            float3  _S5531 = pos_c_1[int(3)];
            _S5463 = &uv_8[int(3)];
            for(;;)
            {
                float2  _S5532 = float2 {_S5531.x, _S5531.y};
                float r_39 = length_0(_S5532);
                float _S5533 = _S5531.z;
                float theta_35 = (F32_atan2((r_39), (_S5533)));
                if(theta_35 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_35 * theta_35 / 3.0f) / _S5533;
                }
                else
                {
                    k_13 = theta_35 / r_39;
                }
                float2  _S5534 = _S5532 * make_float2 (k_13);
                uv_8[int(3)] = _S5534;
                Matrix<float, 2, 2>  _S5535 = camera_distortion_jac_0(_S5534, dist_coeffs_41);
                bool _S5536 = !((F32_min((determinant_0(_S5535)), ((F32_min((_S5535.rows[int(0)].x), (_S5535.rows[int(1)].y)))))) > 0.0f);
                _S5464 = _S5536;
                if(_S5536)
                {
                    break;
                }
                float u_96 = uv_8[int(3)].x;
                float v_96 = uv_8[int(3)].y;
                float r2_91 = u_96 * u_96 + v_96 * v_96;
                float2  _S5537 = uv_8[int(3)] * make_float2 (1.0f + r2_91 * ((*dist_coeffs_41)[int(0)] + r2_91 * ((*dist_coeffs_41)[int(1)] + r2_91 * ((*dist_coeffs_41)[int(2)] + r2_91 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_96 * v_96 + (*dist_coeffs_41)[int(5)] * (r2_91 + 2.0f * u_96 * u_96) + (*dist_coeffs_41)[int(6)] * r2_91, 2.0f * (*dist_coeffs_41)[int(5)] * u_96 * v_96 + (*dist_coeffs_41)[int(4)] * (r2_91 + 2.0f * v_96 * v_96) + (*dist_coeffs_41)[int(7)] * r2_91);
                float2  _S5538 = _S5537 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5537.x + (*dist_coeffs_41)[int(9)] * _S5537.y, 0.0f);
                uv_8[int(3)] = make_float2 (fx_37 * _S5538.x + cx_32, fy_37 * _S5538.y + cy_32);
                break;
            }
            bool all_valid_26 = all_valid_25 & (!_S5464);
            float3  _S5539 = pos_c_1[int(4)];
            _S5465 = &uv_8[int(4)];
            for(;;)
            {
                float2  _S5540 = float2 {_S5539.x, _S5539.y};
                float r_40 = length_0(_S5540);
                float _S5541 = _S5539.z;
                float theta_36 = (F32_atan2((r_40), (_S5541)));
                if(theta_36 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_36 * theta_36 / 3.0f) / _S5541;
                }
                else
                {
                    k_13 = theta_36 / r_40;
                }
                float2  _S5542 = _S5540 * make_float2 (k_13);
                uv_8[int(4)] = _S5542;
                Matrix<float, 2, 2>  _S5543 = camera_distortion_jac_0(_S5542, dist_coeffs_41);
                bool _S5544 = !((F32_min((determinant_0(_S5543)), ((F32_min((_S5543.rows[int(0)].x), (_S5543.rows[int(1)].y)))))) > 0.0f);
                _S5466 = _S5544;
                if(_S5544)
                {
                    break;
                }
                float u_97 = uv_8[int(4)].x;
                float v_97 = uv_8[int(4)].y;
                float r2_92 = u_97 * u_97 + v_97 * v_97;
                float2  _S5545 = uv_8[int(4)] * make_float2 (1.0f + r2_92 * ((*dist_coeffs_41)[int(0)] + r2_92 * ((*dist_coeffs_41)[int(1)] + r2_92 * ((*dist_coeffs_41)[int(2)] + r2_92 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_97 * v_97 + (*dist_coeffs_41)[int(5)] * (r2_92 + 2.0f * u_97 * u_97) + (*dist_coeffs_41)[int(6)] * r2_92, 2.0f * (*dist_coeffs_41)[int(5)] * u_97 * v_97 + (*dist_coeffs_41)[int(4)] * (r2_92 + 2.0f * v_97 * v_97) + (*dist_coeffs_41)[int(7)] * r2_92);
                float2  _S5546 = _S5545 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5545.x + (*dist_coeffs_41)[int(9)] * _S5545.y, 0.0f);
                uv_8[int(4)] = make_float2 (fx_37 * _S5546.x + cx_32, fy_37 * _S5546.y + cy_32);
                break;
            }
            bool all_valid_27 = all_valid_26 & (!_S5466);
            float3  _S5547 = pos_c_1[int(5)];
            _S5467 = &uv_8[int(5)];
            for(;;)
            {
                float2  _S5548 = float2 {_S5547.x, _S5547.y};
                float r_41 = length_0(_S5548);
                float _S5549 = _S5547.z;
                float theta_37 = (F32_atan2((r_41), (_S5549)));
                if(theta_37 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_37 * theta_37 / 3.0f) / _S5549;
                }
                else
                {
                    k_13 = theta_37 / r_41;
                }
                float2  _S5550 = _S5548 * make_float2 (k_13);
                uv_8[int(5)] = _S5550;
                Matrix<float, 2, 2>  _S5551 = camera_distortion_jac_0(_S5550, dist_coeffs_41);
                bool _S5552 = !((F32_min((determinant_0(_S5551)), ((F32_min((_S5551.rows[int(0)].x), (_S5551.rows[int(1)].y)))))) > 0.0f);
                _S5468 = _S5552;
                if(_S5552)
                {
                    break;
                }
                float u_98 = uv_8[int(5)].x;
                float v_98 = uv_8[int(5)].y;
                float r2_93 = u_98 * u_98 + v_98 * v_98;
                float2  _S5553 = uv_8[int(5)] * make_float2 (1.0f + r2_93 * ((*dist_coeffs_41)[int(0)] + r2_93 * ((*dist_coeffs_41)[int(1)] + r2_93 * ((*dist_coeffs_41)[int(2)] + r2_93 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_98 * v_98 + (*dist_coeffs_41)[int(5)] * (r2_93 + 2.0f * u_98 * u_98) + (*dist_coeffs_41)[int(6)] * r2_93, 2.0f * (*dist_coeffs_41)[int(5)] * u_98 * v_98 + (*dist_coeffs_41)[int(4)] * (r2_93 + 2.0f * v_98 * v_98) + (*dist_coeffs_41)[int(7)] * r2_93);
                float2  _S5554 = _S5553 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5553.x + (*dist_coeffs_41)[int(9)] * _S5553.y, 0.0f);
                uv_8[int(5)] = make_float2 (fx_37 * _S5554.x + cx_32, fy_37 * _S5554.y + cy_32);
                break;
            }
            bool all_valid_28 = all_valid_27 & (!_S5468);
            float3  _S5555 = pos_c_1[int(6)];
            _S5469 = &uv_8[int(6)];
            for(;;)
            {
                float2  _S5556 = float2 {_S5555.x, _S5555.y};
                float r_42 = length_0(_S5556);
                float _S5557 = _S5555.z;
                float theta_38 = (F32_atan2((r_42), (_S5557)));
                if(theta_38 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_38 * theta_38 / 3.0f) / _S5557;
                }
                else
                {
                    k_13 = theta_38 / r_42;
                }
                float2  _S5558 = _S5556 * make_float2 (k_13);
                uv_8[int(6)] = _S5558;
                Matrix<float, 2, 2>  _S5559 = camera_distortion_jac_0(_S5558, dist_coeffs_41);
                bool _S5560 = !((F32_min((determinant_0(_S5559)), ((F32_min((_S5559.rows[int(0)].x), (_S5559.rows[int(1)].y)))))) > 0.0f);
                _S5470 = _S5560;
                if(_S5560)
                {
                    break;
                }
                float u_99 = uv_8[int(6)].x;
                float v_99 = uv_8[int(6)].y;
                float r2_94 = u_99 * u_99 + v_99 * v_99;
                float2  _S5561 = uv_8[int(6)] * make_float2 (1.0f + r2_94 * ((*dist_coeffs_41)[int(0)] + r2_94 * ((*dist_coeffs_41)[int(1)] + r2_94 * ((*dist_coeffs_41)[int(2)] + r2_94 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_99 * v_99 + (*dist_coeffs_41)[int(5)] * (r2_94 + 2.0f * u_99 * u_99) + (*dist_coeffs_41)[int(6)] * r2_94, 2.0f * (*dist_coeffs_41)[int(5)] * u_99 * v_99 + (*dist_coeffs_41)[int(4)] * (r2_94 + 2.0f * v_99 * v_99) + (*dist_coeffs_41)[int(7)] * r2_94);
                float2  _S5562 = _S5561 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5561.x + (*dist_coeffs_41)[int(9)] * _S5561.y, 0.0f);
                uv_8[int(6)] = make_float2 (fx_37 * _S5562.x + cx_32, fy_37 * _S5562.y + cy_32);
                break;
            }
            bool all_valid_29 = all_valid_28 & (!_S5470);
            float3  _S5563 = pos_c_1[int(7)];
            _S5471 = &uv_8[int(7)];
            for(;;)
            {
                float2  _S5564 = float2 {_S5563.x, _S5563.y};
                float r_43 = length_0(_S5564);
                float _S5565 = _S5563.z;
                float theta_39 = (F32_atan2((r_43), (_S5565)));
                if(theta_39 < 0.00100000004749745f)
                {
                    k_13 = (1.0f - theta_39 * theta_39 / 3.0f) / _S5565;
                }
                else
                {
                    k_13 = theta_39 / r_43;
                }
                float2  _S5566 = _S5564 * make_float2 (k_13);
                uv_8[int(7)] = _S5566;
                Matrix<float, 2, 2>  _S5567 = camera_distortion_jac_0(_S5566, dist_coeffs_41);
                bool _S5568 = !((F32_min((determinant_0(_S5567)), ((F32_min((_S5567.rows[int(0)].x), (_S5567.rows[int(1)].y)))))) > 0.0f);
                _S5472 = _S5568;
                if(_S5568)
                {
                    break;
                }
                float u_100 = uv_8[int(7)].x;
                float v_100 = uv_8[int(7)].y;
                float r2_95 = u_100 * u_100 + v_100 * v_100;
                float2  _S5569 = uv_8[int(7)] * make_float2 (1.0f + r2_95 * ((*dist_coeffs_41)[int(0)] + r2_95 * ((*dist_coeffs_41)[int(1)] + r2_95 * ((*dist_coeffs_41)[int(2)] + r2_95 * (*dist_coeffs_41)[int(3)])))) + make_float2 (2.0f * (*dist_coeffs_41)[int(4)] * u_100 * v_100 + (*dist_coeffs_41)[int(5)] * (r2_95 + 2.0f * u_100 * u_100) + (*dist_coeffs_41)[int(6)] * r2_95, 2.0f * (*dist_coeffs_41)[int(5)] * u_100 * v_100 + (*dist_coeffs_41)[int(4)] * (r2_95 + 2.0f * v_100 * v_100) + (*dist_coeffs_41)[int(7)] * r2_95);
                float2  _S5570 = _S5569 + make_float2 ((*dist_coeffs_41)[int(8)] * _S5569.x + (*dist_coeffs_41)[int(9)] * _S5569.y, 0.0f);
                uv_8[int(7)] = make_float2 (fx_37 * _S5570.x + cx_32, fy_37 * _S5570.y + cy_32);
                break;
            }
            _S5473 = all_valid_29 & (!_S5472);
            break;
        }
        if(!_S5473)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        if((1.0f - (F32_exp((- (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*densities_1)[int(0)]), ((*densities_1)[int(1)])))), ((*densities_1)[int(2)])))), ((*densities_1)[int(3)])))), ((*densities_1)[int(4)])))), ((*densities_1)[int(5)])))), ((*densities_1)[int(6)])))), ((*densities_1)[int(7)]))) * size_1 * (F32_sqrt((3.0f))))))) <= 0.00392156885936856f)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        float _S5571 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5457).x), ((*_S5459).x)))), ((*_S5461).x)))), ((*_S5463).x)))), ((*_S5465).x)))), ((*_S5467).x)))), ((*_S5469).x)))), ((*_S5471).x)));
        float _S5572 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5457).x), ((*_S5459).x)))), ((*_S5461).x)))), ((*_S5463).x)))), ((*_S5465).x)))), ((*_S5467).x)))), ((*_S5469).x)))), ((*_S5471).x)));
        float _S5573 = (F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((*_S5457).y), ((*_S5459).y)))), ((*_S5461).y)))), ((*_S5463).y)))), ((*_S5465).y)))), ((*_S5467).y)))), ((*_S5469).y)))), ((*_S5471).y)));
        float _S5574 = (F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((*_S5457).y), ((*_S5459).y)))), ((*_S5461).y)))), ((*_S5463).y)))), ((*_S5465).y)))), ((*_S5467).y)))), ((*_S5469).y)))), ((*_S5471).y)));
        if(_S5571 <= 0.0f)
        {
            _S5506 = true;
        }
        else
        {
            _S5506 = _S5572 >= float(image_width_28);
        }
        if(_S5506)
        {
            _S5506 = true;
        }
        else
        {
            _S5506 = _S5573 <= 0.0f;
        }
        if(_S5506)
        {
            _S5506 = true;
        }
        else
        {
            _S5506 = _S5574 >= float(image_height_28);
        }
        if(_S5506)
        {
            _S5506 = true;
        }
        else
        {
            if(_S5504 <= 0.0f)
            {
                if(_S5572 <= 0.0f)
                {
                    _S5506 = _S5571 >= float(image_width_28);
                }
                else
                {
                    _S5506 = false;
                }
                if(_S5506)
                {
                    _S5506 = true;
                }
                else
                {
                    if(_S5574 <= 0.0f)
                    {
                        _S5506 = _S5573 >= float(image_width_28);
                    }
                    else
                    {
                        _S5506 = false;
                    }
                }
            }
            else
            {
                _S5506 = false;
            }
        }
        if(_S5506)
        {
            *aabb_xyxy_19 = make_int4 (int(0), int(0), int(0), int(0));
            break;
        }
        *aabb_xyxy_19 = make_int4 (int((F32_floor((_S5572)))), int((F32_floor((_S5574)))), int((F32_ceil((_S5571)))), int((F32_ceil((_S5573)))));
        *depth_22 = 0.5f * (F32_log((dot_0(mean_c_28, mean_c_28) + 9.99999997475242708e-07f)));
        float3  _S5575 = mean_c_28 - - mul_0(transpose_0(R_34), t_34);
        float3  _S5576 = make_float3 (0.282094806432724f) * (*sh_coeffs_28)[int(0)];
        *rgbs_8 = _S5576;
        float _S5577 = _S5575.x;
        float _S5578 = _S5575.y;
        float _S5579 = _S5575.z;
        float norm_19 = (F32_sqrt((_S5577 * _S5577 + _S5578 * _S5578 + _S5579 * _S5579)));
        float x_64 = _S5577 / norm_19;
        float y_31 = _S5578 / norm_19;
        float z_28 = _S5579 / norm_19;
        float3  _S5580 = _S5576 + make_float3 (0.48860251903533936f) * (make_float3 (- y_31) * (*sh_coeffs_28)[int(1)] + make_float3 (z_28) * (*sh_coeffs_28)[int(2)] - make_float3 (x_64) * (*sh_coeffs_28)[int(3)]);
        *rgbs_8 = _S5580;
        float z2_61 = z_28 * z_28;
        float fTmp0B_28 = -1.09254848957061768f * z_28;
        float fC1_28 = x_64 * x_64 - y_31 * y_31;
        float fS1_28 = 2.0f * x_64 * y_31;
        float3  _S5581 = _S5580 + (make_float3 (0.54627424478530884f * fS1_28) * (*sh_coeffs_28)[int(4)] + make_float3 (fTmp0B_28 * y_31) * (*sh_coeffs_28)[int(5)] + make_float3 (0.94617468118667603f * z2_61 - 0.31539157032966614f) * (*sh_coeffs_28)[int(6)] + make_float3 (fTmp0B_28 * x_64) * (*sh_coeffs_28)[int(7)] + make_float3 (0.54627424478530884f * fC1_28) * (*sh_coeffs_28)[int(8)]);
        *rgbs_8 = _S5581;
        float fTmp0C_28 = -2.28522896766662598f * z2_61 + 0.4570457935333252f;
        float fTmp1B_28 = 1.44530570507049561f * z_28;
        *rgbs_8 = max_0(_S5581 + (make_float3 (-0.59004360437393188f * (x_64 * fS1_28 + y_31 * fC1_28)) * (*sh_coeffs_28)[int(9)] + make_float3 (fTmp1B_28 * fS1_28) * (*sh_coeffs_28)[int(10)] + make_float3 (fTmp0C_28 * y_31) * (*sh_coeffs_28)[int(11)] + make_float3 (z_28 * (1.86588168144226074f * z2_61 - 1.11952900886535645f)) * (*sh_coeffs_28)[int(12)] + make_float3 (fTmp0C_28 * x_64) * (*sh_coeffs_28)[int(13)] + make_float3 (fTmp1B_28 * fC1_28) * (*sh_coeffs_28)[int(14)] + make_float3 (-0.59004360437393188f * (x_64 * fC1_28 - y_31 * fS1_28)) * (*sh_coeffs_28)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
        break;
    }
    return;
}

inline __device__ void _projection_voxel_eval3d_persp_differentiable(float3  pos_2, float size_2, FixedArray<float, 8>  * densities_2, FixedArray<float3 , 16>  * sh_coeffs_29, Matrix<float, 3, 3>  R_35, float3  t_35, float fx_38, float fy_38, float cx_33, float cy_33, FixedArray<float, 10>  * dist_coeffs_42, uint image_width_29, uint image_height_29, float near_plane_20, float far_plane_20, int4  * aabb_xyxy_20, float * depth_23, float3  * rgbs_9)
{
    FixedArray<float3 , 8>  pos_c_2;
    float3  _S5582 = mul_0(R_35, pos_2) + t_35;
    pos_c_2[int(0)] = _S5582;
    float3  _S5583 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 0.0f)) + t_35;
    pos_c_2[int(1)] = _S5583;
    float3  _S5584 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 0.0f)) + t_35;
    pos_c_2[int(2)] = _S5584;
    float3  _S5585 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 1.0f, 0.0f)) + t_35;
    pos_c_2[int(3)] = _S5585;
    float3  _S5586 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 0.0f, 1.0f)) + t_35;
    pos_c_2[int(4)] = _S5586;
    float3  _S5587 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (1.0f, 0.0f, 1.0f)) + t_35;
    pos_c_2[int(5)] = _S5587;
    float3  _S5588 = mul_0(R_35, pos_2 + make_float3 (size_2) * make_float3 (0.0f, 1.0f, 1.0f)) + t_35;
    pos_c_2[int(6)] = _S5588;
    float3  _S5589 = mul_0(R_35, pos_2 + make_float3 (size_2)) + t_35;
    pos_c_2[int(7)] = _S5589;
    float3  mean_c_29 = mul_0(R_35, pos_2 + make_float3 (0.5f * size_2)) + t_35;
    FixedArray<float2 , 8>  uv_9;
    float2  _S5590 = float2 {_S5582.x, _S5582.y} / make_float2 (_S5582.z);
    float u_101 = _S5590.x;
    float v_101 = _S5590.y;
    float r2_96 = u_101 * u_101 + v_101 * v_101;
    float _S5591 = 2.0f * (*dist_coeffs_42)[int(4)];
    float _S5592 = 2.0f * (*dist_coeffs_42)[int(5)];
    float2  _S5593 = _S5590 * make_float2 (1.0f + r2_96 * ((*dist_coeffs_42)[int(0)] + r2_96 * ((*dist_coeffs_42)[int(1)] + r2_96 * ((*dist_coeffs_42)[int(2)] + r2_96 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_101 * v_101 + (*dist_coeffs_42)[int(5)] * (r2_96 + 2.0f * u_101 * u_101) + (*dist_coeffs_42)[int(6)] * r2_96, _S5592 * u_101 * v_101 + (*dist_coeffs_42)[int(4)] * (r2_96 + 2.0f * v_101 * v_101) + (*dist_coeffs_42)[int(7)] * r2_96);
    float2  _S5594 = _S5593 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5593.x + (*dist_coeffs_42)[int(9)] * _S5593.y, 0.0f);
    float _S5595 = fx_38 * _S5594.x + cx_33;
    float _S5596 = fy_38 * _S5594.y + cy_33;
    uv_9[int(0)] = make_float2 (_S5595, _S5596);
    float2  _S5597 = float2 {_S5583.x, _S5583.y} / make_float2 (_S5583.z);
    float u_102 = _S5597.x;
    float v_102 = _S5597.y;
    float r2_97 = u_102 * u_102 + v_102 * v_102;
    float2  _S5598 = _S5597 * make_float2 (1.0f + r2_97 * ((*dist_coeffs_42)[int(0)] + r2_97 * ((*dist_coeffs_42)[int(1)] + r2_97 * ((*dist_coeffs_42)[int(2)] + r2_97 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_102 * v_102 + (*dist_coeffs_42)[int(5)] * (r2_97 + 2.0f * u_102 * u_102) + (*dist_coeffs_42)[int(6)] * r2_97, _S5592 * u_102 * v_102 + (*dist_coeffs_42)[int(4)] * (r2_97 + 2.0f * v_102 * v_102) + (*dist_coeffs_42)[int(7)] * r2_97);
    float2  _S5599 = _S5598 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5598.x + (*dist_coeffs_42)[int(9)] * _S5598.y, 0.0f);
    float _S5600 = fx_38 * _S5599.x + cx_33;
    float _S5601 = fy_38 * _S5599.y + cy_33;
    uv_9[int(1)] = make_float2 (_S5600, _S5601);
    float2  _S5602 = float2 {_S5584.x, _S5584.y} / make_float2 (_S5584.z);
    float u_103 = _S5602.x;
    float v_103 = _S5602.y;
    float r2_98 = u_103 * u_103 + v_103 * v_103;
    float2  _S5603 = _S5602 * make_float2 (1.0f + r2_98 * ((*dist_coeffs_42)[int(0)] + r2_98 * ((*dist_coeffs_42)[int(1)] + r2_98 * ((*dist_coeffs_42)[int(2)] + r2_98 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_103 * v_103 + (*dist_coeffs_42)[int(5)] * (r2_98 + 2.0f * u_103 * u_103) + (*dist_coeffs_42)[int(6)] * r2_98, _S5592 * u_103 * v_103 + (*dist_coeffs_42)[int(4)] * (r2_98 + 2.0f * v_103 * v_103) + (*dist_coeffs_42)[int(7)] * r2_98);
    float2  _S5604 = _S5603 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5603.x + (*dist_coeffs_42)[int(9)] * _S5603.y, 0.0f);
    float _S5605 = fx_38 * _S5604.x + cx_33;
    float _S5606 = fy_38 * _S5604.y + cy_33;
    uv_9[int(2)] = make_float2 (_S5605, _S5606);
    float2  _S5607 = float2 {_S5585.x, _S5585.y} / make_float2 (_S5585.z);
    float u_104 = _S5607.x;
    float v_104 = _S5607.y;
    float r2_99 = u_104 * u_104 + v_104 * v_104;
    float2  _S5608 = _S5607 * make_float2 (1.0f + r2_99 * ((*dist_coeffs_42)[int(0)] + r2_99 * ((*dist_coeffs_42)[int(1)] + r2_99 * ((*dist_coeffs_42)[int(2)] + r2_99 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_104 * v_104 + (*dist_coeffs_42)[int(5)] * (r2_99 + 2.0f * u_104 * u_104) + (*dist_coeffs_42)[int(6)] * r2_99, _S5592 * u_104 * v_104 + (*dist_coeffs_42)[int(4)] * (r2_99 + 2.0f * v_104 * v_104) + (*dist_coeffs_42)[int(7)] * r2_99);
    float2  _S5609 = _S5608 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5608.x + (*dist_coeffs_42)[int(9)] * _S5608.y, 0.0f);
    float _S5610 = fx_38 * _S5609.x + cx_33;
    float _S5611 = fy_38 * _S5609.y + cy_33;
    uv_9[int(3)] = make_float2 (_S5610, _S5611);
    float2  _S5612 = float2 {_S5586.x, _S5586.y} / make_float2 (_S5586.z);
    float u_105 = _S5612.x;
    float v_105 = _S5612.y;
    float r2_100 = u_105 * u_105 + v_105 * v_105;
    float2  _S5613 = _S5612 * make_float2 (1.0f + r2_100 * ((*dist_coeffs_42)[int(0)] + r2_100 * ((*dist_coeffs_42)[int(1)] + r2_100 * ((*dist_coeffs_42)[int(2)] + r2_100 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_105 * v_105 + (*dist_coeffs_42)[int(5)] * (r2_100 + 2.0f * u_105 * u_105) + (*dist_coeffs_42)[int(6)] * r2_100, _S5592 * u_105 * v_105 + (*dist_coeffs_42)[int(4)] * (r2_100 + 2.0f * v_105 * v_105) + (*dist_coeffs_42)[int(7)] * r2_100);
    float2  _S5614 = _S5613 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5613.x + (*dist_coeffs_42)[int(9)] * _S5613.y, 0.0f);
    float _S5615 = fx_38 * _S5614.x + cx_33;
    float _S5616 = fy_38 * _S5614.y + cy_33;
    uv_9[int(4)] = make_float2 (_S5615, _S5616);
    float2  _S5617 = float2 {_S5587.x, _S5587.y} / make_float2 (_S5587.z);
    float u_106 = _S5617.x;
    float v_106 = _S5617.y;
    float r2_101 = u_106 * u_106 + v_106 * v_106;
    float2  _S5618 = _S5617 * make_float2 (1.0f + r2_101 * ((*dist_coeffs_42)[int(0)] + r2_101 * ((*dist_coeffs_42)[int(1)] + r2_101 * ((*dist_coeffs_42)[int(2)] + r2_101 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_106 * v_106 + (*dist_coeffs_42)[int(5)] * (r2_101 + 2.0f * u_106 * u_106) + (*dist_coeffs_42)[int(6)] * r2_101, _S5592 * u_106 * v_106 + (*dist_coeffs_42)[int(4)] * (r2_101 + 2.0f * v_106 * v_106) + (*dist_coeffs_42)[int(7)] * r2_101);
    float2  _S5619 = _S5618 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5618.x + (*dist_coeffs_42)[int(9)] * _S5618.y, 0.0f);
    float _S5620 = fx_38 * _S5619.x + cx_33;
    float _S5621 = fy_38 * _S5619.y + cy_33;
    uv_9[int(5)] = make_float2 (_S5620, _S5621);
    float2  _S5622 = float2 {_S5588.x, _S5588.y} / make_float2 (_S5588.z);
    float u_107 = _S5622.x;
    float v_107 = _S5622.y;
    float r2_102 = u_107 * u_107 + v_107 * v_107;
    float2  _S5623 = _S5622 * make_float2 (1.0f + r2_102 * ((*dist_coeffs_42)[int(0)] + r2_102 * ((*dist_coeffs_42)[int(1)] + r2_102 * ((*dist_coeffs_42)[int(2)] + r2_102 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_107 * v_107 + (*dist_coeffs_42)[int(5)] * (r2_102 + 2.0f * u_107 * u_107) + (*dist_coeffs_42)[int(6)] * r2_102, _S5592 * u_107 * v_107 + (*dist_coeffs_42)[int(4)] * (r2_102 + 2.0f * v_107 * v_107) + (*dist_coeffs_42)[int(7)] * r2_102);
    float2  _S5624 = _S5623 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5623.x + (*dist_coeffs_42)[int(9)] * _S5623.y, 0.0f);
    float _S5625 = fx_38 * _S5624.x + cx_33;
    float _S5626 = fy_38 * _S5624.y + cy_33;
    uv_9[int(6)] = make_float2 (_S5625, _S5626);
    float2  _S5627 = float2 {_S5589.x, _S5589.y} / make_float2 (_S5589.z);
    float u_108 = _S5627.x;
    float v_108 = _S5627.y;
    float r2_103 = u_108 * u_108 + v_108 * v_108;
    float2  _S5628 = _S5627 * make_float2 (1.0f + r2_103 * ((*dist_coeffs_42)[int(0)] + r2_103 * ((*dist_coeffs_42)[int(1)] + r2_103 * ((*dist_coeffs_42)[int(2)] + r2_103 * (*dist_coeffs_42)[int(3)])))) + make_float2 (_S5591 * u_108 * v_108 + (*dist_coeffs_42)[int(5)] * (r2_103 + 2.0f * u_108 * u_108) + (*dist_coeffs_42)[int(6)] * r2_103, _S5592 * u_108 * v_108 + (*dist_coeffs_42)[int(4)] * (r2_103 + 2.0f * v_108 * v_108) + (*dist_coeffs_42)[int(7)] * r2_103);
    float2  _S5629 = _S5628 + make_float2 ((*dist_coeffs_42)[int(8)] * _S5628.x + (*dist_coeffs_42)[int(9)] * _S5628.y, 0.0f);
    float _S5630 = fx_38 * _S5629.x + cx_33;
    float _S5631 = fy_38 * _S5629.y + cy_33;
    uv_9[int(7)] = make_float2 (_S5630, _S5631);
    *aabb_xyxy_20 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5595), (_S5600)))), (_S5605)))), (_S5610)))), (_S5615)))), (_S5620)))), (_S5625)))), (_S5630))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((_S5596), (_S5601)))), (_S5606)))), (_S5611)))), (_S5616)))), (_S5621)))), (_S5626)))), (_S5631))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5595), (_S5600)))), (_S5605)))), (_S5610)))), (_S5615)))), (_S5620)))), (_S5625)))), (_S5630))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((_S5596), (_S5601)))), (_S5606)))), (_S5611)))), (_S5616)))), (_S5621)))), (_S5626)))), (_S5631))))))));
    *depth_23 = 0.5f * (F32_log((dot_0(mean_c_29, mean_c_29) + 9.99999997475242708e-07f)));
    float3  _S5632 = mean_c_29 - - mul_0(transpose_0(R_35), t_35);
    float3  _S5633 = make_float3 (0.282094806432724f) * (*sh_coeffs_29)[int(0)];
    *rgbs_9 = _S5633;
    float _S5634 = _S5632.x;
    float _S5635 = _S5632.y;
    float _S5636 = _S5632.z;
    float norm_20 = (F32_sqrt((_S5634 * _S5634 + _S5635 * _S5635 + _S5636 * _S5636)));
    float x_65 = _S5634 / norm_20;
    float y_32 = _S5635 / norm_20;
    float z_29 = _S5636 / norm_20;
    float3  _S5637 = _S5633 + make_float3 (0.48860251903533936f) * (make_float3 (- y_32) * (*sh_coeffs_29)[int(1)] + make_float3 (z_29) * (*sh_coeffs_29)[int(2)] - make_float3 (x_65) * (*sh_coeffs_29)[int(3)]);
    *rgbs_9 = _S5637;
    float z2_62 = z_29 * z_29;
    float fTmp0B_29 = -1.09254848957061768f * z_29;
    float fC1_29 = x_65 * x_65 - y_32 * y_32;
    float fS1_29 = 2.0f * x_65 * y_32;
    float3  _S5638 = _S5637 + (make_float3 (0.54627424478530884f * fS1_29) * (*sh_coeffs_29)[int(4)] + make_float3 (fTmp0B_29 * y_32) * (*sh_coeffs_29)[int(5)] + make_float3 (0.94617468118667603f * z2_62 - 0.31539157032966614f) * (*sh_coeffs_29)[int(6)] + make_float3 (fTmp0B_29 * x_65) * (*sh_coeffs_29)[int(7)] + make_float3 (0.54627424478530884f * fC1_29) * (*sh_coeffs_29)[int(8)]);
    *rgbs_9 = _S5638;
    float fTmp0C_29 = -2.28522896766662598f * z2_62 + 0.4570457935333252f;
    float fTmp1B_29 = 1.44530570507049561f * z_29;
    *rgbs_9 = max_0(_S5638 + (make_float3 (-0.59004360437393188f * (x_65 * fS1_29 + y_32 * fC1_29)) * (*sh_coeffs_29)[int(9)] + make_float3 (fTmp1B_29 * fS1_29) * (*sh_coeffs_29)[int(10)] + make_float3 (fTmp0C_29 * y_32) * (*sh_coeffs_29)[int(11)] + make_float3 (z_29 * (1.86588168144226074f * z2_62 - 1.11952900886535645f)) * (*sh_coeffs_29)[int(12)] + make_float3 (fTmp0C_29 * x_65) * (*sh_coeffs_29)[int(13)] + make_float3 (fTmp1B_29 * fC1_29) * (*sh_coeffs_29)[int(14)] + make_float3 (-0.59004360437393188f * (x_65 * fC1_29 - y_32 * fS1_29)) * (*sh_coeffs_29)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void _projection_voxel_eval3d_fisheye_differentiable(float3  pos_3, float size_3, FixedArray<float, 8>  * densities_3, FixedArray<float3 , 16>  * sh_coeffs_30, Matrix<float, 3, 3>  R_36, float3  t_36, float fx_39, float fy_39, float cx_34, float cy_34, FixedArray<float, 10>  * dist_coeffs_43, uint image_width_30, uint image_height_30, float near_plane_21, float far_plane_21, int4  * aabb_xyxy_21, float * depth_24, float3  * rgbs_10)
{
    FixedArray<float3 , 8>  pos_c_3;
    float3  _S5639 = mul_0(R_36, pos_3) + t_36;
    pos_c_3[int(0)] = _S5639;
    pos_c_3[int(1)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 0.0f)) + t_36;
    pos_c_3[int(2)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 0.0f)) + t_36;
    pos_c_3[int(3)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 1.0f, 0.0f)) + t_36;
    pos_c_3[int(4)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 0.0f, 1.0f)) + t_36;
    pos_c_3[int(5)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (1.0f, 0.0f, 1.0f)) + t_36;
    pos_c_3[int(6)] = mul_0(R_36, pos_3 + make_float3 (size_3) * make_float3 (0.0f, 1.0f, 1.0f)) + t_36;
    pos_c_3[int(7)] = mul_0(R_36, pos_3 + make_float3 (size_3)) + t_36;
    float3  mean_c_30 = mul_0(R_36, pos_3 + make_float3 (0.5f * size_3)) + t_36;
    FixedArray<float2 , 8>  uv_10;
    float2  _S5640 = float2 {_S5639.x, _S5639.y};
    float r_44 = length_0(_S5640);
    float _S5641 = _S5639.z;
    float theta_40 = (F32_atan2((r_44), (_S5641)));
    float k_14;
    if(theta_40 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_40 * theta_40 / 3.0f) / _S5641;
    }
    else
    {
        k_14 = theta_40 / r_44;
    }
    float2  _S5642 = _S5640 * make_float2 (k_14);
    float u_109 = _S5642.x;
    float v_109 = _S5642.y;
    float r2_104 = u_109 * u_109 + v_109 * v_109;
    float _S5643 = 2.0f * (*dist_coeffs_43)[int(4)];
    float _S5644 = 2.0f * (*dist_coeffs_43)[int(5)];
    float2  _S5645 = _S5642 * make_float2 (1.0f + r2_104 * ((*dist_coeffs_43)[int(0)] + r2_104 * ((*dist_coeffs_43)[int(1)] + r2_104 * ((*dist_coeffs_43)[int(2)] + r2_104 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_109 * v_109 + (*dist_coeffs_43)[int(5)] * (r2_104 + 2.0f * u_109 * u_109) + (*dist_coeffs_43)[int(6)] * r2_104, _S5644 * u_109 * v_109 + (*dist_coeffs_43)[int(4)] * (r2_104 + 2.0f * v_109 * v_109) + (*dist_coeffs_43)[int(7)] * r2_104);
    float2  _S5646 = _S5645 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5645.x + (*dist_coeffs_43)[int(9)] * _S5645.y, 0.0f);
    uv_10[int(0)] = make_float2 (fx_39 * _S5646.x + cx_34, fy_39 * _S5646.y + cy_34);
    float2  _S5647 = float2 {pos_c_3[int(1)].x, pos_c_3[int(1)].y};
    float r_45 = length_0(_S5647);
    float _S5648 = pos_c_3[int(1)].z;
    float theta_41 = (F32_atan2((r_45), (_S5648)));
    if(theta_41 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_41 * theta_41 / 3.0f) / _S5648;
    }
    else
    {
        k_14 = theta_41 / r_45;
    }
    float2  _S5649 = _S5647 * make_float2 (k_14);
    float u_110 = _S5649.x;
    float v_110 = _S5649.y;
    float r2_105 = u_110 * u_110 + v_110 * v_110;
    float2  _S5650 = _S5649 * make_float2 (1.0f + r2_105 * ((*dist_coeffs_43)[int(0)] + r2_105 * ((*dist_coeffs_43)[int(1)] + r2_105 * ((*dist_coeffs_43)[int(2)] + r2_105 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_110 * v_110 + (*dist_coeffs_43)[int(5)] * (r2_105 + 2.0f * u_110 * u_110) + (*dist_coeffs_43)[int(6)] * r2_105, _S5644 * u_110 * v_110 + (*dist_coeffs_43)[int(4)] * (r2_105 + 2.0f * v_110 * v_110) + (*dist_coeffs_43)[int(7)] * r2_105);
    float2  _S5651 = _S5650 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5650.x + (*dist_coeffs_43)[int(9)] * _S5650.y, 0.0f);
    uv_10[int(1)] = make_float2 (fx_39 * _S5651.x + cx_34, fy_39 * _S5651.y + cy_34);
    float2  _S5652 = float2 {pos_c_3[int(2)].x, pos_c_3[int(2)].y};
    float r_46 = length_0(_S5652);
    float _S5653 = pos_c_3[int(2)].z;
    float theta_42 = (F32_atan2((r_46), (_S5653)));
    if(theta_42 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_42 * theta_42 / 3.0f) / _S5653;
    }
    else
    {
        k_14 = theta_42 / r_46;
    }
    float2  _S5654 = _S5652 * make_float2 (k_14);
    float u_111 = _S5654.x;
    float v_111 = _S5654.y;
    float r2_106 = u_111 * u_111 + v_111 * v_111;
    float2  _S5655 = _S5654 * make_float2 (1.0f + r2_106 * ((*dist_coeffs_43)[int(0)] + r2_106 * ((*dist_coeffs_43)[int(1)] + r2_106 * ((*dist_coeffs_43)[int(2)] + r2_106 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_111 * v_111 + (*dist_coeffs_43)[int(5)] * (r2_106 + 2.0f * u_111 * u_111) + (*dist_coeffs_43)[int(6)] * r2_106, _S5644 * u_111 * v_111 + (*dist_coeffs_43)[int(4)] * (r2_106 + 2.0f * v_111 * v_111) + (*dist_coeffs_43)[int(7)] * r2_106);
    float2  _S5656 = _S5655 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5655.x + (*dist_coeffs_43)[int(9)] * _S5655.y, 0.0f);
    uv_10[int(2)] = make_float2 (fx_39 * _S5656.x + cx_34, fy_39 * _S5656.y + cy_34);
    float2  _S5657 = float2 {pos_c_3[int(3)].x, pos_c_3[int(3)].y};
    float r_47 = length_0(_S5657);
    float _S5658 = pos_c_3[int(3)].z;
    float theta_43 = (F32_atan2((r_47), (_S5658)));
    if(theta_43 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_43 * theta_43 / 3.0f) / _S5658;
    }
    else
    {
        k_14 = theta_43 / r_47;
    }
    float2  _S5659 = _S5657 * make_float2 (k_14);
    float u_112 = _S5659.x;
    float v_112 = _S5659.y;
    float r2_107 = u_112 * u_112 + v_112 * v_112;
    float2  _S5660 = _S5659 * make_float2 (1.0f + r2_107 * ((*dist_coeffs_43)[int(0)] + r2_107 * ((*dist_coeffs_43)[int(1)] + r2_107 * ((*dist_coeffs_43)[int(2)] + r2_107 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_112 * v_112 + (*dist_coeffs_43)[int(5)] * (r2_107 + 2.0f * u_112 * u_112) + (*dist_coeffs_43)[int(6)] * r2_107, _S5644 * u_112 * v_112 + (*dist_coeffs_43)[int(4)] * (r2_107 + 2.0f * v_112 * v_112) + (*dist_coeffs_43)[int(7)] * r2_107);
    float2  _S5661 = _S5660 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5660.x + (*dist_coeffs_43)[int(9)] * _S5660.y, 0.0f);
    uv_10[int(3)] = make_float2 (fx_39 * _S5661.x + cx_34, fy_39 * _S5661.y + cy_34);
    float2  _S5662 = float2 {pos_c_3[int(4)].x, pos_c_3[int(4)].y};
    float r_48 = length_0(_S5662);
    float _S5663 = pos_c_3[int(4)].z;
    float theta_44 = (F32_atan2((r_48), (_S5663)));
    if(theta_44 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_44 * theta_44 / 3.0f) / _S5663;
    }
    else
    {
        k_14 = theta_44 / r_48;
    }
    float2  _S5664 = _S5662 * make_float2 (k_14);
    float u_113 = _S5664.x;
    float v_113 = _S5664.y;
    float r2_108 = u_113 * u_113 + v_113 * v_113;
    float2  _S5665 = _S5664 * make_float2 (1.0f + r2_108 * ((*dist_coeffs_43)[int(0)] + r2_108 * ((*dist_coeffs_43)[int(1)] + r2_108 * ((*dist_coeffs_43)[int(2)] + r2_108 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_113 * v_113 + (*dist_coeffs_43)[int(5)] * (r2_108 + 2.0f * u_113 * u_113) + (*dist_coeffs_43)[int(6)] * r2_108, _S5644 * u_113 * v_113 + (*dist_coeffs_43)[int(4)] * (r2_108 + 2.0f * v_113 * v_113) + (*dist_coeffs_43)[int(7)] * r2_108);
    float2  _S5666 = _S5665 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5665.x + (*dist_coeffs_43)[int(9)] * _S5665.y, 0.0f);
    uv_10[int(4)] = make_float2 (fx_39 * _S5666.x + cx_34, fy_39 * _S5666.y + cy_34);
    float2  _S5667 = float2 {pos_c_3[int(5)].x, pos_c_3[int(5)].y};
    float r_49 = length_0(_S5667);
    float _S5668 = pos_c_3[int(5)].z;
    float theta_45 = (F32_atan2((r_49), (_S5668)));
    if(theta_45 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_45 * theta_45 / 3.0f) / _S5668;
    }
    else
    {
        k_14 = theta_45 / r_49;
    }
    float2  _S5669 = _S5667 * make_float2 (k_14);
    float u_114 = _S5669.x;
    float v_114 = _S5669.y;
    float r2_109 = u_114 * u_114 + v_114 * v_114;
    float2  _S5670 = _S5669 * make_float2 (1.0f + r2_109 * ((*dist_coeffs_43)[int(0)] + r2_109 * ((*dist_coeffs_43)[int(1)] + r2_109 * ((*dist_coeffs_43)[int(2)] + r2_109 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_114 * v_114 + (*dist_coeffs_43)[int(5)] * (r2_109 + 2.0f * u_114 * u_114) + (*dist_coeffs_43)[int(6)] * r2_109, _S5644 * u_114 * v_114 + (*dist_coeffs_43)[int(4)] * (r2_109 + 2.0f * v_114 * v_114) + (*dist_coeffs_43)[int(7)] * r2_109);
    float2  _S5671 = _S5670 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5670.x + (*dist_coeffs_43)[int(9)] * _S5670.y, 0.0f);
    uv_10[int(5)] = make_float2 (fx_39 * _S5671.x + cx_34, fy_39 * _S5671.y + cy_34);
    float2  _S5672 = float2 {pos_c_3[int(6)].x, pos_c_3[int(6)].y};
    float r_50 = length_0(_S5672);
    float _S5673 = pos_c_3[int(6)].z;
    float theta_46 = (F32_atan2((r_50), (_S5673)));
    if(theta_46 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_46 * theta_46 / 3.0f) / _S5673;
    }
    else
    {
        k_14 = theta_46 / r_50;
    }
    float2  _S5674 = _S5672 * make_float2 (k_14);
    float u_115 = _S5674.x;
    float v_115 = _S5674.y;
    float r2_110 = u_115 * u_115 + v_115 * v_115;
    float2  _S5675 = _S5674 * make_float2 (1.0f + r2_110 * ((*dist_coeffs_43)[int(0)] + r2_110 * ((*dist_coeffs_43)[int(1)] + r2_110 * ((*dist_coeffs_43)[int(2)] + r2_110 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_115 * v_115 + (*dist_coeffs_43)[int(5)] * (r2_110 + 2.0f * u_115 * u_115) + (*dist_coeffs_43)[int(6)] * r2_110, _S5644 * u_115 * v_115 + (*dist_coeffs_43)[int(4)] * (r2_110 + 2.0f * v_115 * v_115) + (*dist_coeffs_43)[int(7)] * r2_110);
    float2  _S5676 = _S5675 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5675.x + (*dist_coeffs_43)[int(9)] * _S5675.y, 0.0f);
    uv_10[int(6)] = make_float2 (fx_39 * _S5676.x + cx_34, fy_39 * _S5676.y + cy_34);
    float2  _S5677 = float2 {pos_c_3[int(7)].x, pos_c_3[int(7)].y};
    float r_51 = length_0(_S5677);
    float _S5678 = pos_c_3[int(7)].z;
    float theta_47 = (F32_atan2((r_51), (_S5678)));
    if(theta_47 < 0.00100000004749745f)
    {
        k_14 = (1.0f - theta_47 * theta_47 / 3.0f) / _S5678;
    }
    else
    {
        k_14 = theta_47 / r_51;
    }
    float2  _S5679 = _S5677 * make_float2 (k_14);
    float u_116 = _S5679.x;
    float v_116 = _S5679.y;
    float r2_111 = u_116 * u_116 + v_116 * v_116;
    float2  _S5680 = _S5679 * make_float2 (1.0f + r2_111 * ((*dist_coeffs_43)[int(0)] + r2_111 * ((*dist_coeffs_43)[int(1)] + r2_111 * ((*dist_coeffs_43)[int(2)] + r2_111 * (*dist_coeffs_43)[int(3)])))) + make_float2 (_S5643 * u_116 * v_116 + (*dist_coeffs_43)[int(5)] * (r2_111 + 2.0f * u_116 * u_116) + (*dist_coeffs_43)[int(6)] * r2_111, _S5644 * u_116 * v_116 + (*dist_coeffs_43)[int(4)] * (r2_111 + 2.0f * v_116 * v_116) + (*dist_coeffs_43)[int(7)] * r2_111);
    float2  _S5681 = _S5680 + make_float2 ((*dist_coeffs_43)[int(8)] * _S5680.x + (*dist_coeffs_43)[int(9)] * _S5680.y, 0.0f);
    float _S5682 = fx_39 * _S5681.x + cx_34;
    float _S5683 = fy_39 * _S5681.y + cy_34;
    uv_10[int(7)] = make_float2 (_S5682, _S5683);
    *aabb_xyxy_21 = make_int4 (int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_10[int(0)].x), (uv_10[int(1)].x)))), (uv_10[int(2)].x)))), (uv_10[int(3)].x)))), (uv_10[int(4)].x)))), (uv_10[int(5)].x)))), (uv_10[int(6)].x)))), (_S5682))))))), int((F32_floor(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min(((F32_min((uv_10[int(0)].y), (uv_10[int(1)].y)))), (uv_10[int(2)].y)))), (uv_10[int(3)].y)))), (uv_10[int(4)].y)))), (uv_10[int(5)].y)))), (uv_10[int(6)].y)))), (_S5683))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_10[int(0)].x), (uv_10[int(1)].x)))), (uv_10[int(2)].x)))), (uv_10[int(3)].x)))), (uv_10[int(4)].x)))), (uv_10[int(5)].x)))), (uv_10[int(6)].x)))), (_S5682))))))), int((F32_ceil(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max(((F32_max((uv_10[int(0)].y), (uv_10[int(1)].y)))), (uv_10[int(2)].y)))), (uv_10[int(3)].y)))), (uv_10[int(4)].y)))), (uv_10[int(5)].y)))), (uv_10[int(6)].y)))), (_S5683))))))));
    *depth_24 = 0.5f * (F32_log((dot_0(mean_c_30, mean_c_30) + 9.99999997475242708e-07f)));
    float3  _S5684 = mean_c_30 - - mul_0(transpose_0(R_36), t_36);
    float3  _S5685 = make_float3 (0.282094806432724f) * (*sh_coeffs_30)[int(0)];
    *rgbs_10 = _S5685;
    float _S5686 = _S5684.x;
    float _S5687 = _S5684.y;
    float _S5688 = _S5684.z;
    float norm_21 = (F32_sqrt((_S5686 * _S5686 + _S5687 * _S5687 + _S5688 * _S5688)));
    float x_66 = _S5686 / norm_21;
    float y_33 = _S5687 / norm_21;
    float z_30 = _S5688 / norm_21;
    float3  _S5689 = _S5685 + make_float3 (0.48860251903533936f) * (make_float3 (- y_33) * (*sh_coeffs_30)[int(1)] + make_float3 (z_30) * (*sh_coeffs_30)[int(2)] - make_float3 (x_66) * (*sh_coeffs_30)[int(3)]);
    *rgbs_10 = _S5689;
    float z2_63 = z_30 * z_30;
    float fTmp0B_30 = -1.09254848957061768f * z_30;
    float fC1_30 = x_66 * x_66 - y_33 * y_33;
    float fS1_30 = 2.0f * x_66 * y_33;
    float3  _S5690 = _S5689 + (make_float3 (0.54627424478530884f * fS1_30) * (*sh_coeffs_30)[int(4)] + make_float3 (fTmp0B_30 * y_33) * (*sh_coeffs_30)[int(5)] + make_float3 (0.94617468118667603f * z2_63 - 0.31539157032966614f) * (*sh_coeffs_30)[int(6)] + make_float3 (fTmp0B_30 * x_66) * (*sh_coeffs_30)[int(7)] + make_float3 (0.54627424478530884f * fC1_30) * (*sh_coeffs_30)[int(8)]);
    *rgbs_10 = _S5690;
    float fTmp0C_30 = -2.28522896766662598f * z2_63 + 0.4570457935333252f;
    float fTmp1B_30 = 1.44530570507049561f * z_30;
    *rgbs_10 = max_0(_S5690 + (make_float3 (-0.59004360437393188f * (x_66 * fS1_30 + y_33 * fC1_30)) * (*sh_coeffs_30)[int(9)] + make_float3 (fTmp1B_30 * fS1_30) * (*sh_coeffs_30)[int(10)] + make_float3 (fTmp0C_30 * y_33) * (*sh_coeffs_30)[int(11)] + make_float3 (z_30 * (1.86588168144226074f * z2_63 - 1.11952900886535645f)) * (*sh_coeffs_30)[int(12)] + make_float3 (fTmp0C_30 * x_66) * (*sh_coeffs_30)[int(13)] + make_float3 (fTmp1B_30 * fC1_30) * (*sh_coeffs_30)[int(14)] + make_float3 (-0.59004360437393188f * (x_66 * fC1_30 - y_33 * fS1_30)) * (*sh_coeffs_30)[int(15)]) + make_float3 (0.5f), make_float3 (0.0f));
    return;
}

inline __device__ void projection_voxel_eval3d_persp_vjp(float3  pos_4, float size_4, FixedArray<float, 8>  * densities_4, FixedArray<float3 , 16>  * sh_coeffs_31, Matrix<float, 3, 3>  R_37, float3  t_37, float fx_40, float fy_40, float cx_35, float cy_35, FixedArray<float, 10>  * dist_coeffs_44, uint image_width_31, uint image_height_31, float3  v_rgb_8, FixedArray<float, 8>  * v_densities_0, FixedArray<float3 , 16>  * v_sh_coeffs_9, Matrix<float, 3, 3>  * v_R_11, float3  * v_t_10)
{
    float3  _S5691 = s_primal_ctx_mul_0(R_37, pos_4) + t_37;
    float _S5692 = _S5691.z;
    float2  _S5693 = make_float2 (_S5692);
    float _S5694 = s_primal_ctx_min_0(1.00000001504746622e+30f, _S5692);
    float _S5695 = s_primal_ctx_max_0(0.0f, _S5692);
    float3  pos_i_0 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S5696 = s_primal_ctx_mul_0(R_37, pos_i_0) + t_37;
    float _S5697 = _S5696.z;
    float2  _S5698 = make_float2 (_S5697);
    float _S5699 = s_primal_ctx_min_0(_S5694, _S5697);
    float _S5700 = s_primal_ctx_max_0(_S5695, _S5697);
    float3  pos_i_1 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S5701 = s_primal_ctx_mul_0(R_37, pos_i_1) + t_37;
    float _S5702 = _S5701.z;
    float2  _S5703 = make_float2 (_S5702);
    float _S5704 = s_primal_ctx_min_0(_S5699, _S5702);
    float _S5705 = s_primal_ctx_max_0(_S5700, _S5702);
    float3  pos_i_2 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S5706 = s_primal_ctx_mul_0(R_37, pos_i_2) + t_37;
    float _S5707 = _S5706.z;
    float2  _S5708 = make_float2 (_S5707);
    float _S5709 = s_primal_ctx_min_0(_S5704, _S5707);
    float _S5710 = s_primal_ctx_max_0(_S5705, _S5707);
    float3  pos_i_3 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S5711 = s_primal_ctx_mul_0(R_37, pos_i_3) + t_37;
    float _S5712 = _S5711.z;
    float2  _S5713 = make_float2 (_S5712);
    float _S5714 = s_primal_ctx_min_0(_S5709, _S5712);
    float _S5715 = s_primal_ctx_max_0(_S5710, _S5712);
    float3  pos_i_4 = pos_4 + make_float3 (size_4) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S5716 = s_primal_ctx_mul_0(R_37, pos_i_4) + t_37;
    float _S5717 = _S5716.z;
    float2  _S5718 = make_float2 (_S5717);
    float _S5719 = s_primal_ctx_min_0(_S5714, _S5717);
    float _S5720 = s_primal_ctx_max_0(_S5715, _S5717);
    float3  pos_i_5 = pos_4 + make_float3 (size_4) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S5721 = s_primal_ctx_mul_0(R_37, pos_i_5) + t_37;
    float _S5722 = _S5721.z;
    float2  _S5723 = make_float2 (_S5722);
    float _S5724 = s_primal_ctx_min_0(_S5719, _S5722);
    float _S5725 = s_primal_ctx_max_0(_S5720, _S5722);
    float3  pos_i_6 = pos_4 + make_float3 (size_4);
    float3  _S5726 = s_primal_ctx_mul_0(R_37, pos_i_6) + t_37;
    float _S5727 = _S5726.z;
    float2  _S5728 = make_float2 (_S5727);
    float3  _S5729 = pos_4 + make_float3 (0.5f * size_4);
    float3  mean_c_31 = s_primal_ctx_mul_0(R_37, _S5729) + t_37;
    float2  _S5730 = float2 {_S5691.x, _S5691.y};
    float2  _S5731 = _S5730 / make_float2 (_S5692);
    float2  _S5732 = make_float2 (_S5692 * _S5692);
    float u_117 = _S5731.x;
    float v_117 = _S5731.y;
    float r2_112 = u_117 * u_117 + v_117 * v_117;
    float _S5733 = (*dist_coeffs_44)[int(2)] + r2_112 * (*dist_coeffs_44)[int(3)];
    float _S5734 = (*dist_coeffs_44)[int(1)] + r2_112 * _S5733;
    float _S5735 = (*dist_coeffs_44)[int(0)] + r2_112 * _S5734;
    float radial_13 = 1.0f + r2_112 * _S5735;
    float2  _S5736 = make_float2 (radial_13);
    float _S5737 = 2.0f * (*dist_coeffs_44)[int(4)];
    float _S5738 = _S5737 * u_117;
    float _S5739 = 2.0f * u_117;
    float _S5740 = 2.0f * (*dist_coeffs_44)[int(5)];
    float _S5741 = _S5740 * u_117;
    float _S5742 = 2.0f * v_117;
    float2  _S5743 = _S5731 * make_float2 (radial_13) + make_float2 (_S5738 * v_117 + (*dist_coeffs_44)[int(5)] * (r2_112 + _S5739 * u_117) + (*dist_coeffs_44)[int(6)] * r2_112, _S5741 * v_117 + (*dist_coeffs_44)[int(4)] * (r2_112 + _S5742 * v_117) + (*dist_coeffs_44)[int(7)] * r2_112);
    float2  _S5744 = _S5743 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5743.x + (*dist_coeffs_44)[int(9)] * _S5743.y, 0.0f);
    float _S5745 = fx_40 * _S5744.x + cx_35;
    float _S5746 = fy_40 * _S5744.y + cy_35;
    float2  _S5747 = float2 {_S5696.x, _S5696.y};
    float2  _S5748 = _S5747 / make_float2 (_S5697);
    float2  _S5749 = make_float2 (_S5697 * _S5697);
    float u_118 = _S5748.x;
    float v_118 = _S5748.y;
    float r2_113 = u_118 * u_118 + v_118 * v_118;
    float _S5750 = (*dist_coeffs_44)[int(2)] + r2_113 * (*dist_coeffs_44)[int(3)];
    float _S5751 = (*dist_coeffs_44)[int(1)] + r2_113 * _S5750;
    float _S5752 = (*dist_coeffs_44)[int(0)] + r2_113 * _S5751;
    float radial_14 = 1.0f + r2_113 * _S5752;
    float2  _S5753 = make_float2 (radial_14);
    float _S5754 = _S5737 * u_118;
    float _S5755 = 2.0f * u_118;
    float _S5756 = _S5740 * u_118;
    float _S5757 = 2.0f * v_118;
    float2  _S5758 = _S5748 * make_float2 (radial_14) + make_float2 (_S5754 * v_118 + (*dist_coeffs_44)[int(5)] * (r2_113 + _S5755 * u_118) + (*dist_coeffs_44)[int(6)] * r2_113, _S5756 * v_118 + (*dist_coeffs_44)[int(4)] * (r2_113 + _S5757 * v_118) + (*dist_coeffs_44)[int(7)] * r2_113);
    float2  _S5759 = _S5758 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5758.x + (*dist_coeffs_44)[int(9)] * _S5758.y, 0.0f);
    float _S5760 = fx_40 * _S5759.x + cx_35;
    float _S5761 = fy_40 * _S5759.y + cy_35;
    float2  _S5762 = float2 {_S5701.x, _S5701.y};
    float2  _S5763 = _S5762 / make_float2 (_S5702);
    float2  _S5764 = make_float2 (_S5702 * _S5702);
    float u_119 = _S5763.x;
    float v_119 = _S5763.y;
    float r2_114 = u_119 * u_119 + v_119 * v_119;
    float _S5765 = (*dist_coeffs_44)[int(2)] + r2_114 * (*dist_coeffs_44)[int(3)];
    float _S5766 = (*dist_coeffs_44)[int(1)] + r2_114 * _S5765;
    float _S5767 = (*dist_coeffs_44)[int(0)] + r2_114 * _S5766;
    float radial_15 = 1.0f + r2_114 * _S5767;
    float2  _S5768 = make_float2 (radial_15);
    float _S5769 = _S5737 * u_119;
    float _S5770 = 2.0f * u_119;
    float _S5771 = _S5740 * u_119;
    float _S5772 = 2.0f * v_119;
    float2  _S5773 = _S5763 * make_float2 (radial_15) + make_float2 (_S5769 * v_119 + (*dist_coeffs_44)[int(5)] * (r2_114 + _S5770 * u_119) + (*dist_coeffs_44)[int(6)] * r2_114, _S5771 * v_119 + (*dist_coeffs_44)[int(4)] * (r2_114 + _S5772 * v_119) + (*dist_coeffs_44)[int(7)] * r2_114);
    float2  _S5774 = _S5773 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5773.x + (*dist_coeffs_44)[int(9)] * _S5773.y, 0.0f);
    float _S5775 = fx_40 * _S5774.x + cx_35;
    float _S5776 = fy_40 * _S5774.y + cy_35;
    float2  _S5777 = float2 {_S5706.x, _S5706.y};
    float2  _S5778 = _S5777 / make_float2 (_S5707);
    float2  _S5779 = make_float2 (_S5707 * _S5707);
    float u_120 = _S5778.x;
    float v_120 = _S5778.y;
    float r2_115 = u_120 * u_120 + v_120 * v_120;
    float _S5780 = (*dist_coeffs_44)[int(2)] + r2_115 * (*dist_coeffs_44)[int(3)];
    float _S5781 = (*dist_coeffs_44)[int(1)] + r2_115 * _S5780;
    float _S5782 = (*dist_coeffs_44)[int(0)] + r2_115 * _S5781;
    float radial_16 = 1.0f + r2_115 * _S5782;
    float2  _S5783 = make_float2 (radial_16);
    float _S5784 = _S5737 * u_120;
    float _S5785 = 2.0f * u_120;
    float _S5786 = _S5740 * u_120;
    float _S5787 = 2.0f * v_120;
    float2  _S5788 = _S5778 * make_float2 (radial_16) + make_float2 (_S5784 * v_120 + (*dist_coeffs_44)[int(5)] * (r2_115 + _S5785 * u_120) + (*dist_coeffs_44)[int(6)] * r2_115, _S5786 * v_120 + (*dist_coeffs_44)[int(4)] * (r2_115 + _S5787 * v_120) + (*dist_coeffs_44)[int(7)] * r2_115);
    float2  _S5789 = _S5788 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5788.x + (*dist_coeffs_44)[int(9)] * _S5788.y, 0.0f);
    float _S5790 = fx_40 * _S5789.x + cx_35;
    float _S5791 = fy_40 * _S5789.y + cy_35;
    float2  _S5792 = float2 {_S5711.x, _S5711.y};
    float2  _S5793 = _S5792 / make_float2 (_S5712);
    float2  _S5794 = make_float2 (_S5712 * _S5712);
    float u_121 = _S5793.x;
    float v_121 = _S5793.y;
    float r2_116 = u_121 * u_121 + v_121 * v_121;
    float _S5795 = (*dist_coeffs_44)[int(2)] + r2_116 * (*dist_coeffs_44)[int(3)];
    float _S5796 = (*dist_coeffs_44)[int(1)] + r2_116 * _S5795;
    float _S5797 = (*dist_coeffs_44)[int(0)] + r2_116 * _S5796;
    float radial_17 = 1.0f + r2_116 * _S5797;
    float2  _S5798 = make_float2 (radial_17);
    float _S5799 = _S5737 * u_121;
    float _S5800 = 2.0f * u_121;
    float _S5801 = _S5740 * u_121;
    float _S5802 = 2.0f * v_121;
    float2  _S5803 = _S5793 * make_float2 (radial_17) + make_float2 (_S5799 * v_121 + (*dist_coeffs_44)[int(5)] * (r2_116 + _S5800 * u_121) + (*dist_coeffs_44)[int(6)] * r2_116, _S5801 * v_121 + (*dist_coeffs_44)[int(4)] * (r2_116 + _S5802 * v_121) + (*dist_coeffs_44)[int(7)] * r2_116);
    float2  _S5804 = _S5803 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5803.x + (*dist_coeffs_44)[int(9)] * _S5803.y, 0.0f);
    float _S5805 = fx_40 * _S5804.x + cx_35;
    float _S5806 = fy_40 * _S5804.y + cy_35;
    float2  _S5807 = float2 {_S5716.x, _S5716.y};
    float2  _S5808 = _S5807 / make_float2 (_S5717);
    float2  _S5809 = make_float2 (_S5717 * _S5717);
    float u_122 = _S5808.x;
    float v_122 = _S5808.y;
    float r2_117 = u_122 * u_122 + v_122 * v_122;
    float _S5810 = (*dist_coeffs_44)[int(2)] + r2_117 * (*dist_coeffs_44)[int(3)];
    float _S5811 = (*dist_coeffs_44)[int(1)] + r2_117 * _S5810;
    float _S5812 = (*dist_coeffs_44)[int(0)] + r2_117 * _S5811;
    float radial_18 = 1.0f + r2_117 * _S5812;
    float2  _S5813 = make_float2 (radial_18);
    float _S5814 = _S5737 * u_122;
    float _S5815 = 2.0f * u_122;
    float _S5816 = _S5740 * u_122;
    float _S5817 = 2.0f * v_122;
    float2  _S5818 = _S5808 * make_float2 (radial_18) + make_float2 (_S5814 * v_122 + (*dist_coeffs_44)[int(5)] * (r2_117 + _S5815 * u_122) + (*dist_coeffs_44)[int(6)] * r2_117, _S5816 * v_122 + (*dist_coeffs_44)[int(4)] * (r2_117 + _S5817 * v_122) + (*dist_coeffs_44)[int(7)] * r2_117);
    float2  _S5819 = _S5818 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5818.x + (*dist_coeffs_44)[int(9)] * _S5818.y, 0.0f);
    float _S5820 = fx_40 * _S5819.x + cx_35;
    float _S5821 = fy_40 * _S5819.y + cy_35;
    float2  _S5822 = float2 {_S5721.x, _S5721.y};
    float2  _S5823 = _S5822 / make_float2 (_S5722);
    float2  _S5824 = make_float2 (_S5722 * _S5722);
    float u_123 = _S5823.x;
    float v_123 = _S5823.y;
    float r2_118 = u_123 * u_123 + v_123 * v_123;
    float _S5825 = (*dist_coeffs_44)[int(2)] + r2_118 * (*dist_coeffs_44)[int(3)];
    float _S5826 = (*dist_coeffs_44)[int(1)] + r2_118 * _S5825;
    float _S5827 = (*dist_coeffs_44)[int(0)] + r2_118 * _S5826;
    float radial_19 = 1.0f + r2_118 * _S5827;
    float2  _S5828 = make_float2 (radial_19);
    float _S5829 = _S5737 * u_123;
    float _S5830 = 2.0f * u_123;
    float _S5831 = _S5740 * u_123;
    float _S5832 = 2.0f * v_123;
    float2  _S5833 = _S5823 * make_float2 (radial_19) + make_float2 (_S5829 * v_123 + (*dist_coeffs_44)[int(5)] * (r2_118 + _S5830 * u_123) + (*dist_coeffs_44)[int(6)] * r2_118, _S5831 * v_123 + (*dist_coeffs_44)[int(4)] * (r2_118 + _S5832 * v_123) + (*dist_coeffs_44)[int(7)] * r2_118);
    float2  _S5834 = _S5833 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5833.x + (*dist_coeffs_44)[int(9)] * _S5833.y, 0.0f);
    float _S5835 = fx_40 * _S5834.x + cx_35;
    float _S5836 = fy_40 * _S5834.y + cy_35;
    float2  _S5837 = float2 {_S5726.x, _S5726.y};
    float2  _S5838 = _S5837 / make_float2 (_S5727);
    float2  _S5839 = make_float2 (_S5727 * _S5727);
    float u_124 = _S5838.x;
    float v_124 = _S5838.y;
    float r2_119 = u_124 * u_124 + v_124 * v_124;
    float _S5840 = (*dist_coeffs_44)[int(2)] + r2_119 * (*dist_coeffs_44)[int(3)];
    float _S5841 = (*dist_coeffs_44)[int(1)] + r2_119 * _S5840;
    float _S5842 = (*dist_coeffs_44)[int(0)] + r2_119 * _S5841;
    float radial_20 = 1.0f + r2_119 * _S5842;
    float2  _S5843 = make_float2 (radial_20);
    float _S5844 = _S5737 * u_124;
    float _S5845 = 2.0f * u_124;
    float _S5846 = _S5740 * u_124;
    float _S5847 = 2.0f * v_124;
    float2  _S5848 = _S5838 * make_float2 (radial_20) + make_float2 (_S5844 * v_124 + (*dist_coeffs_44)[int(5)] * (r2_119 + _S5845 * u_124) + (*dist_coeffs_44)[int(6)] * r2_119, _S5846 * v_124 + (*dist_coeffs_44)[int(4)] * (r2_119 + _S5847 * v_124) + (*dist_coeffs_44)[int(7)] * r2_119);
    float2  _S5849 = _S5848 + make_float2 ((*dist_coeffs_44)[int(8)] * _S5848.x + (*dist_coeffs_44)[int(9)] * _S5848.y, 0.0f);
    float _S5850 = fx_40 * _S5849.x + cx_35;
    float _S5851 = fy_40 * _S5849.y + cy_35;
    float _S5852 = s_primal_ctx_max_0(_S5745, _S5760);
    float _S5853 = s_primal_ctx_min_0(_S5745, _S5760);
    float _S5854 = s_primal_ctx_max_0(_S5746, _S5761);
    float _S5855 = s_primal_ctx_min_0(_S5746, _S5761);
    float _S5856 = s_primal_ctx_max_0(_S5852, _S5775);
    float _S5857 = s_primal_ctx_min_0(_S5853, _S5775);
    float _S5858 = s_primal_ctx_max_0(_S5854, _S5776);
    float _S5859 = s_primal_ctx_min_0(_S5855, _S5776);
    float _S5860 = s_primal_ctx_max_0(_S5856, _S5790);
    float _S5861 = s_primal_ctx_min_0(_S5857, _S5790);
    float _S5862 = s_primal_ctx_max_0(_S5858, _S5791);
    float _S5863 = s_primal_ctx_min_0(_S5859, _S5791);
    float _S5864 = s_primal_ctx_max_0(_S5860, _S5805);
    float _S5865 = s_primal_ctx_min_0(_S5861, _S5805);
    float _S5866 = s_primal_ctx_max_0(_S5862, _S5806);
    float _S5867 = s_primal_ctx_min_0(_S5863, _S5806);
    float _S5868 = s_primal_ctx_max_0(_S5864, _S5820);
    float _S5869 = s_primal_ctx_min_0(_S5865, _S5820);
    float _S5870 = s_primal_ctx_max_0(_S5866, _S5821);
    float _S5871 = s_primal_ctx_min_0(_S5867, _S5821);
    float _S5872 = s_primal_ctx_max_0(_S5868, _S5835);
    float _S5873 = s_primal_ctx_min_0(_S5869, _S5835);
    float _S5874 = s_primal_ctx_max_0(_S5870, _S5836);
    float _S5875 = s_primal_ctx_min_0(_S5871, _S5836);
    float _S5876 = s_primal_ctx_dot_0(mean_c_31, mean_c_31) + 9.99999997475242708e-07f;
    Matrix<float, 3, 3>  _S5877 = transpose_0(R_37);
    float3  _S5878 = mean_c_31 - - s_primal_ctx_mul_0(_S5877, t_37);
    float _S5879 = _S5878.x;
    float _S5880 = _S5878.y;
    float _S5881 = _S5878.z;
    float _S5882 = _S5879 * _S5879 + _S5880 * _S5880 + _S5881 * _S5881;
    float _S5883 = s_primal_ctx_sqrt_0(_S5882);
    float x_67 = _S5879 / _S5883;
    float3  _S5884 = make_float3 (x_67);
    float _S5885 = _S5883 * _S5883;
    float y_34 = _S5880 / _S5883;
    float z_31 = _S5881 / _S5883;
    float3  _S5886 = make_float3 (z_31);
    float _S5887 = - y_34;
    float3  _S5888 = make_float3 (_S5887);
    float z2_64 = z_31 * z_31;
    float fTmp0B_31 = -1.09254848957061768f * z_31;
    float fC1_31 = x_67 * x_67 - y_34 * y_34;
    float _S5889 = 2.0f * x_67;
    float fS1_31 = _S5889 * y_34;
    float pSH6_9 = 0.94617468118667603f * z2_64 - 0.31539157032966614f;
    float3  _S5890 = make_float3 (pSH6_9);
    float pSH7_9 = fTmp0B_31 * x_67;
    float3  _S5891 = make_float3 (pSH7_9);
    float pSH5_9 = fTmp0B_31 * y_34;
    float3  _S5892 = make_float3 (pSH5_9);
    float pSH8_9 = 0.54627424478530884f * fC1_31;
    float3  _S5893 = make_float3 (pSH8_9);
    float pSH4_9 = 0.54627424478530884f * fS1_31;
    float3  _S5894 = make_float3 (pSH4_9);
    float fTmp0C_31 = -2.28522896766662598f * z2_64 + 0.4570457935333252f;
    float fTmp1B_31 = 1.44530570507049561f * z_31;
    float _S5895 = 1.86588168144226074f * z2_64 - 1.11952900886535645f;
    float pSH12_9 = z_31 * _S5895;
    float3  _S5896 = make_float3 (pSH12_9);
    float pSH13_9 = fTmp0C_31 * x_67;
    float3  _S5897 = make_float3 (pSH13_9);
    float pSH11_9 = fTmp0C_31 * y_34;
    float3  _S5898 = make_float3 (pSH11_9);
    float pSH14_9 = fTmp1B_31 * fC1_31;
    float3  _S5899 = make_float3 (pSH14_9);
    float pSH10_9 = fTmp1B_31 * fS1_31;
    float3  _S5900 = make_float3 (pSH10_9);
    float pSH15_9 = -0.59004360437393188f * (x_67 * fC1_31 - y_34 * fS1_31);
    float3  _S5901 = make_float3 (pSH15_9);
    float pSH9_9 = -0.59004360437393188f * (x_67 * fS1_31 + y_34 * fC1_31);
    float3  _S5902 = make_float3 (pSH9_9);
    float3  _S5903 = make_float3 (0.0f);
    float3  _S5904 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5905;
    (&_S5905)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_31)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S5887) * (*sh_coeffs_31)[int(1)] + make_float3 (z_31) * (*sh_coeffs_31)[int(2)] - make_float3 (x_67) * (*sh_coeffs_31)[int(3)]) + (make_float3 (pSH4_9) * (*sh_coeffs_31)[int(4)] + make_float3 (pSH5_9) * (*sh_coeffs_31)[int(5)] + make_float3 (pSH6_9) * (*sh_coeffs_31)[int(6)] + make_float3 (pSH7_9) * (*sh_coeffs_31)[int(7)] + make_float3 (pSH8_9) * (*sh_coeffs_31)[int(8)]) + (make_float3 (pSH9_9) * (*sh_coeffs_31)[int(9)] + make_float3 (pSH10_9) * (*sh_coeffs_31)[int(10)] + make_float3 (pSH11_9) * (*sh_coeffs_31)[int(11)] + make_float3 (pSH12_9) * (*sh_coeffs_31)[int(12)] + make_float3 (pSH13_9) * (*sh_coeffs_31)[int(13)] + make_float3 (pSH14_9) * (*sh_coeffs_31)[int(14)] + make_float3 (pSH15_9) * (*sh_coeffs_31)[int(15)]) + make_float3 (0.5f);
    (&_S5905)->differential_0 = _S5904;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5906;
    (&_S5906)->primal_0 = _S5903;
    (&_S5906)->differential_0 = _S5904;
    s_bwd_prop_max_0(&_S5905, &_S5906, v_rgb_8);
    float3  _S5907 = _S5901 * _S5905.differential_0;
    float3  _S5908 = (*sh_coeffs_31)[int(15)] * _S5905.differential_0;
    float3  _S5909 = _S5899 * _S5905.differential_0;
    float3  _S5910 = (*sh_coeffs_31)[int(14)] * _S5905.differential_0;
    float3  _S5911 = _S5897 * _S5905.differential_0;
    float3  _S5912 = (*sh_coeffs_31)[int(13)] * _S5905.differential_0;
    float3  _S5913 = _S5896 * _S5905.differential_0;
    float3  _S5914 = (*sh_coeffs_31)[int(12)] * _S5905.differential_0;
    float3  _S5915 = _S5898 * _S5905.differential_0;
    float3  _S5916 = (*sh_coeffs_31)[int(11)] * _S5905.differential_0;
    float3  _S5917 = _S5900 * _S5905.differential_0;
    float3  _S5918 = (*sh_coeffs_31)[int(10)] * _S5905.differential_0;
    float3  _S5919 = _S5902 * _S5905.differential_0;
    float3  _S5920 = (*sh_coeffs_31)[int(9)] * _S5905.differential_0;
    float s_diff_fS2_T_9 = -0.59004360437393188f * (_S5920.x + _S5920.y + _S5920.z);
    float s_diff_fC2_T_9 = -0.59004360437393188f * (_S5908.x + _S5908.y + _S5908.z);
    float _S5921 = _S5918.x + _S5918.y + _S5918.z;
    float _S5922 = _S5910.x + _S5910.y + _S5910.z;
    float _S5923 = _S5916.x + _S5916.y + _S5916.z;
    float _S5924 = _S5912.x + _S5912.y + _S5912.z;
    float _S5925 = _S5914.x + _S5914.y + _S5914.z;
    float _S5926 = - s_diff_fC2_T_9;
    float3  _S5927 = _S5893 * _S5905.differential_0;
    float3  _S5928 = (*sh_coeffs_31)[int(8)] * _S5905.differential_0;
    float3  _S5929 = _S5891 * _S5905.differential_0;
    float3  _S5930 = (*sh_coeffs_31)[int(7)] * _S5905.differential_0;
    float3  _S5931 = _S5890 * _S5905.differential_0;
    float3  _S5932 = (*sh_coeffs_31)[int(6)] * _S5905.differential_0;
    float3  _S5933 = _S5892 * _S5905.differential_0;
    float3  _S5934 = (*sh_coeffs_31)[int(5)] * _S5905.differential_0;
    float3  _S5935 = _S5894 * _S5905.differential_0;
    float3  _S5936 = (*sh_coeffs_31)[int(4)] * _S5905.differential_0;
    float _S5937 = _S5934.x + _S5934.y + _S5934.z;
    float _S5938 = _S5930.x + _S5930.y + _S5930.z;
    float _S5939 = fTmp1B_31 * _S5921 + x_67 * s_diff_fS2_T_9 + y_34 * _S5926 + 0.54627424478530884f * (_S5936.x + _S5936.y + _S5936.z);
    float _S5940 = fTmp1B_31 * _S5922 + y_34 * s_diff_fS2_T_9 + x_67 * s_diff_fC2_T_9 + 0.54627424478530884f * (_S5928.x + _S5928.y + _S5928.z);
    float _S5941 = y_34 * - _S5940;
    float _S5942 = x_67 * _S5940;
    float _S5943 = z_31 * (1.86588168144226074f * (z_31 * _S5925) + -2.28522896766662598f * (y_34 * _S5923 + x_67 * _S5924) + 0.94617468118667603f * (_S5932.x + _S5932.y + _S5932.z));
    float3  _S5944 = make_float3 (0.48860251903533936f) * _S5905.differential_0;
    float3  _S5945 = - _S5944;
    float3  _S5946 = _S5884 * _S5945;
    float3  _S5947 = (*sh_coeffs_31)[int(3)] * _S5945;
    float3  _S5948 = _S5886 * _S5944;
    float3  _S5949 = (*sh_coeffs_31)[int(2)] * _S5944;
    float3  _S5950 = _S5888 * _S5944;
    float3  _S5951 = (*sh_coeffs_31)[int(1)] * _S5944;
    float _S5952 = (_S5895 * _S5925 + 1.44530570507049561f * (fS1_31 * _S5921 + fC1_31 * _S5922) + -1.09254848957061768f * (y_34 * _S5937 + x_67 * _S5938) + _S5943 + _S5943 + _S5949.x + _S5949.y + _S5949.z) / _S5885;
    float _S5953 = _S5883 * _S5952;
    float _S5954 = (fTmp0C_31 * _S5923 + fC1_31 * s_diff_fS2_T_9 + fS1_31 * _S5926 + fTmp0B_31 * _S5937 + _S5889 * _S5939 + _S5941 + _S5941 + - (_S5951.x + _S5951.y + _S5951.z)) / _S5885;
    float _S5955 = _S5883 * _S5954;
    float _S5956 = (fTmp0C_31 * _S5924 + fS1_31 * s_diff_fS2_T_9 + fC1_31 * s_diff_fC2_T_9 + fTmp0B_31 * _S5938 + 2.0f * (y_34 * _S5939) + _S5942 + _S5942 + _S5947.x + _S5947.y + _S5947.z) / _S5885;
    float _S5957 = _S5883 * _S5956;
    float _S5958 = _S5881 * - _S5952 + _S5880 * - _S5954 + _S5879 * - _S5956;
    DiffPair_float_0 _S5959;
    (&_S5959)->primal_0 = _S5882;
    (&_S5959)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S5959, _S5958);
    float _S5960 = _S5881 * _S5959.differential_0;
    float _S5961 = _S5880 * _S5959.differential_0;
    float _S5962 = _S5879 * _S5959.differential_0;
    float3  _S5963 = make_float3 (0.282094806432724f) * _S5905.differential_0;
    float3  _S5964 = make_float3 (_S5957 + _S5962 + _S5962, _S5955 + _S5961 + _S5961, _S5953 + _S5960 + _S5960);
    float3  _S5965 = - - _S5964;
    Matrix<float, 3, 3>  _S5966 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S5967;
    (&_S5967)->primal_0 = _S5877;
    (&_S5967)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5968;
    (&_S5968)->primal_0 = t_37;
    (&_S5968)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S5967, &_S5968, _S5965);
    Matrix<float, 3, 3>  _S5969 = transpose_0(_S5967.differential_0);
    DiffPair_float_0 _S5970;
    (&_S5970)->primal_0 = _S5876;
    (&_S5970)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S5970, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5971;
    (&_S5971)->primal_0 = mean_c_31;
    (&_S5971)->differential_0 = _S5904;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S5972;
    (&_S5972)->primal_0 = mean_c_31;
    (&_S5972)->differential_0 = _S5904;
    s_bwd_prop_dot_0(&_S5971, &_S5972, _S5970.differential_0);
    DiffPair_float_0 _S5973;
    (&_S5973)->primal_0 = _S5875;
    (&_S5973)->differential_0 = 0.0f;
    DiffPair_float_0 _S5974;
    (&_S5974)->primal_0 = _S5851;
    (&_S5974)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5973, &_S5974, 0.0f);
    DiffPair_float_0 _S5975;
    (&_S5975)->primal_0 = _S5874;
    (&_S5975)->differential_0 = 0.0f;
    DiffPair_float_0 _S5976;
    (&_S5976)->primal_0 = _S5851;
    (&_S5976)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5975, &_S5976, 0.0f);
    DiffPair_float_0 _S5977;
    (&_S5977)->primal_0 = _S5873;
    (&_S5977)->differential_0 = 0.0f;
    DiffPair_float_0 _S5978;
    (&_S5978)->primal_0 = _S5850;
    (&_S5978)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5977, &_S5978, 0.0f);
    DiffPair_float_0 _S5979;
    (&_S5979)->primal_0 = _S5872;
    (&_S5979)->differential_0 = 0.0f;
    DiffPair_float_0 _S5980;
    (&_S5980)->primal_0 = _S5850;
    (&_S5980)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5979, &_S5980, 0.0f);
    DiffPair_float_0 _S5981;
    (&_S5981)->primal_0 = _S5871;
    (&_S5981)->differential_0 = 0.0f;
    DiffPair_float_0 _S5982;
    (&_S5982)->primal_0 = _S5836;
    (&_S5982)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5981, &_S5982, _S5973.differential_0);
    DiffPair_float_0 _S5983;
    (&_S5983)->primal_0 = _S5870;
    (&_S5983)->differential_0 = 0.0f;
    DiffPair_float_0 _S5984;
    (&_S5984)->primal_0 = _S5836;
    (&_S5984)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5983, &_S5984, _S5975.differential_0);
    DiffPair_float_0 _S5985;
    (&_S5985)->primal_0 = _S5869;
    (&_S5985)->differential_0 = 0.0f;
    DiffPair_float_0 _S5986;
    (&_S5986)->primal_0 = _S5835;
    (&_S5986)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5985, &_S5986, _S5977.differential_0);
    DiffPair_float_0 _S5987;
    (&_S5987)->primal_0 = _S5868;
    (&_S5987)->differential_0 = 0.0f;
    DiffPair_float_0 _S5988;
    (&_S5988)->primal_0 = _S5835;
    (&_S5988)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5987, &_S5988, _S5979.differential_0);
    DiffPair_float_0 _S5989;
    (&_S5989)->primal_0 = _S5867;
    (&_S5989)->differential_0 = 0.0f;
    DiffPair_float_0 _S5990;
    (&_S5990)->primal_0 = _S5821;
    (&_S5990)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5989, &_S5990, _S5981.differential_0);
    DiffPair_float_0 _S5991;
    (&_S5991)->primal_0 = _S5866;
    (&_S5991)->differential_0 = 0.0f;
    DiffPair_float_0 _S5992;
    (&_S5992)->primal_0 = _S5821;
    (&_S5992)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5991, &_S5992, _S5983.differential_0);
    DiffPair_float_0 _S5993;
    (&_S5993)->primal_0 = _S5865;
    (&_S5993)->differential_0 = 0.0f;
    DiffPair_float_0 _S5994;
    (&_S5994)->primal_0 = _S5820;
    (&_S5994)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5993, &_S5994, _S5985.differential_0);
    DiffPair_float_0 _S5995;
    (&_S5995)->primal_0 = _S5864;
    (&_S5995)->differential_0 = 0.0f;
    DiffPair_float_0 _S5996;
    (&_S5996)->primal_0 = _S5820;
    (&_S5996)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5995, &_S5996, _S5987.differential_0);
    DiffPair_float_0 _S5997;
    (&_S5997)->primal_0 = _S5863;
    (&_S5997)->differential_0 = 0.0f;
    DiffPair_float_0 _S5998;
    (&_S5998)->primal_0 = _S5806;
    (&_S5998)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S5997, &_S5998, _S5989.differential_0);
    DiffPair_float_0 _S5999;
    (&_S5999)->primal_0 = _S5862;
    (&_S5999)->differential_0 = 0.0f;
    DiffPair_float_0 _S6000;
    (&_S6000)->primal_0 = _S5806;
    (&_S6000)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S5999, &_S6000, _S5991.differential_0);
    DiffPair_float_0 _S6001;
    (&_S6001)->primal_0 = _S5861;
    (&_S6001)->differential_0 = 0.0f;
    DiffPair_float_0 _S6002;
    (&_S6002)->primal_0 = _S5805;
    (&_S6002)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6001, &_S6002, _S5993.differential_0);
    DiffPair_float_0 _S6003;
    (&_S6003)->primal_0 = _S5860;
    (&_S6003)->differential_0 = 0.0f;
    DiffPair_float_0 _S6004;
    (&_S6004)->primal_0 = _S5805;
    (&_S6004)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6003, &_S6004, _S5995.differential_0);
    DiffPair_float_0 _S6005;
    (&_S6005)->primal_0 = _S5859;
    (&_S6005)->differential_0 = 0.0f;
    DiffPair_float_0 _S6006;
    (&_S6006)->primal_0 = _S5791;
    (&_S6006)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6005, &_S6006, _S5997.differential_0);
    DiffPair_float_0 _S6007;
    (&_S6007)->primal_0 = _S5858;
    (&_S6007)->differential_0 = 0.0f;
    DiffPair_float_0 _S6008;
    (&_S6008)->primal_0 = _S5791;
    (&_S6008)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6007, &_S6008, _S5999.differential_0);
    DiffPair_float_0 _S6009;
    (&_S6009)->primal_0 = _S5857;
    (&_S6009)->differential_0 = 0.0f;
    DiffPair_float_0 _S6010;
    (&_S6010)->primal_0 = _S5790;
    (&_S6010)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6009, &_S6010, _S6001.differential_0);
    DiffPair_float_0 _S6011;
    (&_S6011)->primal_0 = _S5856;
    (&_S6011)->differential_0 = 0.0f;
    DiffPair_float_0 _S6012;
    (&_S6012)->primal_0 = _S5790;
    (&_S6012)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6011, &_S6012, _S6003.differential_0);
    DiffPair_float_0 _S6013;
    (&_S6013)->primal_0 = _S5855;
    (&_S6013)->differential_0 = 0.0f;
    DiffPair_float_0 _S6014;
    (&_S6014)->primal_0 = _S5776;
    (&_S6014)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6013, &_S6014, _S6005.differential_0);
    DiffPair_float_0 _S6015;
    (&_S6015)->primal_0 = _S5854;
    (&_S6015)->differential_0 = 0.0f;
    DiffPair_float_0 _S6016;
    (&_S6016)->primal_0 = _S5776;
    (&_S6016)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6015, &_S6016, _S6007.differential_0);
    DiffPair_float_0 _S6017;
    (&_S6017)->primal_0 = _S5853;
    (&_S6017)->differential_0 = 0.0f;
    DiffPair_float_0 _S6018;
    (&_S6018)->primal_0 = _S5775;
    (&_S6018)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6017, &_S6018, _S6009.differential_0);
    DiffPair_float_0 _S6019;
    (&_S6019)->primal_0 = _S5852;
    (&_S6019)->differential_0 = 0.0f;
    DiffPair_float_0 _S6020;
    (&_S6020)->primal_0 = _S5775;
    (&_S6020)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6019, &_S6020, _S6011.differential_0);
    DiffPair_float_0 _S6021;
    (&_S6021)->primal_0 = _S5746;
    (&_S6021)->differential_0 = 0.0f;
    DiffPair_float_0 _S6022;
    (&_S6022)->primal_0 = _S5761;
    (&_S6022)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6021, &_S6022, _S6013.differential_0);
    DiffPair_float_0 _S6023;
    (&_S6023)->primal_0 = _S5746;
    (&_S6023)->differential_0 = 0.0f;
    DiffPair_float_0 _S6024;
    (&_S6024)->primal_0 = _S5761;
    (&_S6024)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6023, &_S6024, _S6015.differential_0);
    DiffPair_float_0 _S6025;
    (&_S6025)->primal_0 = _S5745;
    (&_S6025)->differential_0 = 0.0f;
    DiffPair_float_0 _S6026;
    (&_S6026)->primal_0 = _S5760;
    (&_S6026)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6025, &_S6026, _S6017.differential_0);
    DiffPair_float_0 _S6027;
    (&_S6027)->primal_0 = _S5745;
    (&_S6027)->differential_0 = 0.0f;
    DiffPair_float_0 _S6028;
    (&_S6028)->primal_0 = _S5760;
    (&_S6028)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6027, &_S6028, _S6019.differential_0);
    float _S6029 = fx_40 * (_S5978.differential_0 + _S5980.differential_0);
    float2  _S6030 = make_float2 (_S6029, fy_40 * (_S5974.differential_0 + _S5976.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6029, (*dist_coeffs_44)[int(9)] * _S6029);
    float2  _S6031 = _S5838 * _S6030;
    float _S6032 = (*dist_coeffs_44)[int(4)] * _S6030.y;
    float _S6033 = (*dist_coeffs_44)[int(5)] * _S6030.x;
    float _S6034 = _S6031.x + _S6031.y;
    float _S6035 = r2_119 * _S6034;
    float _S6036 = r2_119 * _S6035;
    float _S6037 = (*dist_coeffs_44)[int(7)] * _S6030.y + _S6032 + (*dist_coeffs_44)[int(6)] * _S6030.x + _S6033 + _S5842 * _S6034 + _S5841 * _S6035 + _S5840 * _S6036 + (*dist_coeffs_44)[int(3)] * (r2_119 * _S6036);
    float _S6038 = v_124 * _S6037;
    float _S6039 = u_124 * _S6037;
    float2  _S6040 = (_S5843 * _S6030 + make_float2 (_S5740 * (v_124 * _S6030.y) + _S5845 * _S6033 + 2.0f * (u_124 * _S6033) + _S5737 * (v_124 * _S6030.x) + _S6039 + _S6039, _S5847 * _S6032 + 2.0f * (v_124 * _S6032) + _S5846 * _S6030.y + _S5844 * _S6030.x + _S6038 + _S6038)) / _S5839;
    float2  _S6041 = _S5837 * - _S6040;
    float2  _S6042 = _S5728 * _S6040;
    float _S6043 = fx_40 * (_S5986.differential_0 + _S5988.differential_0);
    float2  _S6044 = make_float2 (_S6043, fy_40 * (_S5982.differential_0 + _S5984.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6043, (*dist_coeffs_44)[int(9)] * _S6043);
    float2  _S6045 = _S5823 * _S6044;
    float _S6046 = (*dist_coeffs_44)[int(4)] * _S6044.y;
    float _S6047 = (*dist_coeffs_44)[int(5)] * _S6044.x;
    float _S6048 = _S6045.x + _S6045.y;
    float _S6049 = r2_118 * _S6048;
    float _S6050 = r2_118 * _S6049;
    float _S6051 = (*dist_coeffs_44)[int(7)] * _S6044.y + _S6046 + (*dist_coeffs_44)[int(6)] * _S6044.x + _S6047 + _S5827 * _S6048 + _S5826 * _S6049 + _S5825 * _S6050 + (*dist_coeffs_44)[int(3)] * (r2_118 * _S6050);
    float _S6052 = v_123 * _S6051;
    float _S6053 = u_123 * _S6051;
    float2  _S6054 = (_S5828 * _S6044 + make_float2 (_S5740 * (v_123 * _S6044.y) + _S5830 * _S6047 + 2.0f * (u_123 * _S6047) + _S5737 * (v_123 * _S6044.x) + _S6053 + _S6053, _S5832 * _S6046 + 2.0f * (v_123 * _S6046) + _S5831 * _S6044.y + _S5829 * _S6044.x + _S6052 + _S6052)) / _S5824;
    float2  _S6055 = _S5822 * - _S6054;
    float2  _S6056 = _S5723 * _S6054;
    float _S6057 = fx_40 * (_S5994.differential_0 + _S5996.differential_0);
    float2  _S6058 = make_float2 (_S6057, fy_40 * (_S5990.differential_0 + _S5992.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6057, (*dist_coeffs_44)[int(9)] * _S6057);
    float2  _S6059 = _S5808 * _S6058;
    float _S6060 = (*dist_coeffs_44)[int(4)] * _S6058.y;
    float _S6061 = (*dist_coeffs_44)[int(5)] * _S6058.x;
    float _S6062 = _S6059.x + _S6059.y;
    float _S6063 = r2_117 * _S6062;
    float _S6064 = r2_117 * _S6063;
    float _S6065 = (*dist_coeffs_44)[int(7)] * _S6058.y + _S6060 + (*dist_coeffs_44)[int(6)] * _S6058.x + _S6061 + _S5812 * _S6062 + _S5811 * _S6063 + _S5810 * _S6064 + (*dist_coeffs_44)[int(3)] * (r2_117 * _S6064);
    float _S6066 = v_122 * _S6065;
    float _S6067 = u_122 * _S6065;
    float2  _S6068 = (_S5813 * _S6058 + make_float2 (_S5740 * (v_122 * _S6058.y) + _S5815 * _S6061 + 2.0f * (u_122 * _S6061) + _S5737 * (v_122 * _S6058.x) + _S6067 + _S6067, _S5817 * _S6060 + 2.0f * (v_122 * _S6060) + _S5816 * _S6058.y + _S5814 * _S6058.x + _S6066 + _S6066)) / _S5809;
    float2  _S6069 = _S5807 * - _S6068;
    float2  _S6070 = _S5718 * _S6068;
    float _S6071 = fx_40 * (_S6002.differential_0 + _S6004.differential_0);
    float2  _S6072 = make_float2 (_S6071, fy_40 * (_S5998.differential_0 + _S6000.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6071, (*dist_coeffs_44)[int(9)] * _S6071);
    float2  _S6073 = _S5793 * _S6072;
    float _S6074 = (*dist_coeffs_44)[int(4)] * _S6072.y;
    float _S6075 = (*dist_coeffs_44)[int(5)] * _S6072.x;
    float _S6076 = _S6073.x + _S6073.y;
    float _S6077 = r2_116 * _S6076;
    float _S6078 = r2_116 * _S6077;
    float _S6079 = (*dist_coeffs_44)[int(7)] * _S6072.y + _S6074 + (*dist_coeffs_44)[int(6)] * _S6072.x + _S6075 + _S5797 * _S6076 + _S5796 * _S6077 + _S5795 * _S6078 + (*dist_coeffs_44)[int(3)] * (r2_116 * _S6078);
    float _S6080 = v_121 * _S6079;
    float _S6081 = u_121 * _S6079;
    float2  _S6082 = (_S5798 * _S6072 + make_float2 (_S5740 * (v_121 * _S6072.y) + _S5800 * _S6075 + 2.0f * (u_121 * _S6075) + _S5737 * (v_121 * _S6072.x) + _S6081 + _S6081, _S5802 * _S6074 + 2.0f * (v_121 * _S6074) + _S5801 * _S6072.y + _S5799 * _S6072.x + _S6080 + _S6080)) / _S5794;
    float2  _S6083 = _S5792 * - _S6082;
    float2  _S6084 = _S5713 * _S6082;
    float _S6085 = fx_40 * (_S6010.differential_0 + _S6012.differential_0);
    float2  _S6086 = make_float2 (_S6085, fy_40 * (_S6006.differential_0 + _S6008.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6085, (*dist_coeffs_44)[int(9)] * _S6085);
    float2  _S6087 = _S5778 * _S6086;
    float _S6088 = (*dist_coeffs_44)[int(4)] * _S6086.y;
    float _S6089 = (*dist_coeffs_44)[int(5)] * _S6086.x;
    float _S6090 = _S6087.x + _S6087.y;
    float _S6091 = r2_115 * _S6090;
    float _S6092 = r2_115 * _S6091;
    float _S6093 = (*dist_coeffs_44)[int(7)] * _S6086.y + _S6088 + (*dist_coeffs_44)[int(6)] * _S6086.x + _S6089 + _S5782 * _S6090 + _S5781 * _S6091 + _S5780 * _S6092 + (*dist_coeffs_44)[int(3)] * (r2_115 * _S6092);
    float _S6094 = v_120 * _S6093;
    float _S6095 = u_120 * _S6093;
    float2  _S6096 = (_S5783 * _S6086 + make_float2 (_S5740 * (v_120 * _S6086.y) + _S5785 * _S6089 + 2.0f * (u_120 * _S6089) + _S5737 * (v_120 * _S6086.x) + _S6095 + _S6095, _S5787 * _S6088 + 2.0f * (v_120 * _S6088) + _S5786 * _S6086.y + _S5784 * _S6086.x + _S6094 + _S6094)) / _S5779;
    float2  _S6097 = _S5777 * - _S6096;
    float2  _S6098 = _S5708 * _S6096;
    float _S6099 = fx_40 * (_S6018.differential_0 + _S6020.differential_0);
    float2  _S6100 = make_float2 (_S6099, fy_40 * (_S6014.differential_0 + _S6016.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6099, (*dist_coeffs_44)[int(9)] * _S6099);
    float2  _S6101 = _S5763 * _S6100;
    float _S6102 = (*dist_coeffs_44)[int(4)] * _S6100.y;
    float _S6103 = (*dist_coeffs_44)[int(5)] * _S6100.x;
    float _S6104 = _S6101.x + _S6101.y;
    float _S6105 = r2_114 * _S6104;
    float _S6106 = r2_114 * _S6105;
    float _S6107 = (*dist_coeffs_44)[int(7)] * _S6100.y + _S6102 + (*dist_coeffs_44)[int(6)] * _S6100.x + _S6103 + _S5767 * _S6104 + _S5766 * _S6105 + _S5765 * _S6106 + (*dist_coeffs_44)[int(3)] * (r2_114 * _S6106);
    float _S6108 = v_119 * _S6107;
    float _S6109 = u_119 * _S6107;
    float2  _S6110 = (_S5768 * _S6100 + make_float2 (_S5740 * (v_119 * _S6100.y) + _S5770 * _S6103 + 2.0f * (u_119 * _S6103) + _S5737 * (v_119 * _S6100.x) + _S6109 + _S6109, _S5772 * _S6102 + 2.0f * (v_119 * _S6102) + _S5771 * _S6100.y + _S5769 * _S6100.x + _S6108 + _S6108)) / _S5764;
    float2  _S6111 = _S5762 * - _S6110;
    float2  _S6112 = _S5703 * _S6110;
    float _S6113 = fx_40 * (_S6026.differential_0 + _S6028.differential_0);
    float2  _S6114 = make_float2 (_S6113, fy_40 * (_S6022.differential_0 + _S6024.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6113, (*dist_coeffs_44)[int(9)] * _S6113);
    float2  _S6115 = _S5748 * _S6114;
    float _S6116 = (*dist_coeffs_44)[int(4)] * _S6114.y;
    float _S6117 = (*dist_coeffs_44)[int(5)] * _S6114.x;
    float _S6118 = _S6115.x + _S6115.y;
    float _S6119 = r2_113 * _S6118;
    float _S6120 = r2_113 * _S6119;
    float _S6121 = (*dist_coeffs_44)[int(7)] * _S6114.y + _S6116 + (*dist_coeffs_44)[int(6)] * _S6114.x + _S6117 + _S5752 * _S6118 + _S5751 * _S6119 + _S5750 * _S6120 + (*dist_coeffs_44)[int(3)] * (r2_113 * _S6120);
    float _S6122 = v_118 * _S6121;
    float _S6123 = u_118 * _S6121;
    float2  _S6124 = (_S5753 * _S6114 + make_float2 (_S5740 * (v_118 * _S6114.y) + _S5755 * _S6117 + 2.0f * (u_118 * _S6117) + _S5737 * (v_118 * _S6114.x) + _S6123 + _S6123, _S5757 * _S6116 + 2.0f * (v_118 * _S6116) + _S5756 * _S6114.y + _S5754 * _S6114.x + _S6122 + _S6122)) / _S5749;
    float2  _S6125 = _S5747 * - _S6124;
    float2  _S6126 = _S5698 * _S6124;
    float _S6127 = fx_40 * (_S6025.differential_0 + _S6027.differential_0);
    float2  _S6128 = make_float2 (_S6127, fy_40 * (_S6021.differential_0 + _S6023.differential_0)) + make_float2 ((*dist_coeffs_44)[int(8)] * _S6127, (*dist_coeffs_44)[int(9)] * _S6127);
    float2  _S6129 = _S5731 * _S6128;
    float _S6130 = (*dist_coeffs_44)[int(4)] * _S6128.y;
    float _S6131 = (*dist_coeffs_44)[int(5)] * _S6128.x;
    float _S6132 = _S6129.x + _S6129.y;
    float _S6133 = r2_112 * _S6132;
    float _S6134 = r2_112 * _S6133;
    float _S6135 = (*dist_coeffs_44)[int(7)] * _S6128.y + _S6130 + (*dist_coeffs_44)[int(6)] * _S6128.x + _S6131 + _S5735 * _S6132 + _S5734 * _S6133 + _S5733 * _S6134 + (*dist_coeffs_44)[int(3)] * (r2_112 * _S6134);
    float _S6136 = v_117 * _S6135;
    float _S6137 = u_117 * _S6135;
    float2  _S6138 = (_S5736 * _S6128 + make_float2 (_S5740 * (v_117 * _S6128.y) + _S5739 * _S6131 + 2.0f * (u_117 * _S6131) + _S5737 * (v_117 * _S6128.x) + _S6137 + _S6137, _S5742 * _S6130 + 2.0f * (v_117 * _S6130) + _S5741 * _S6128.y + _S5738 * _S6128.x + _S6136 + _S6136)) / _S5732;
    float2  _S6139 = _S5730 * - _S6138;
    float2  _S6140 = _S5693 * _S6138;
    float3  _S6141 = _S5964 + _S5972.differential_0 + _S5971.differential_0;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6142;
    (&_S6142)->primal_0 = R_37;
    (&_S6142)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6143;
    (&_S6143)->primal_0 = _S5729;
    (&_S6143)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6142, &_S6143, _S6141);
    DiffPair_float_0 _S6144;
    (&_S6144)->primal_0 = _S5725;
    (&_S6144)->differential_0 = 0.0f;
    DiffPair_float_0 _S6145;
    (&_S6145)->primal_0 = _S5727;
    (&_S6145)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6144, &_S6145, 0.0f);
    DiffPair_float_0 _S6146;
    (&_S6146)->primal_0 = _S5724;
    (&_S6146)->differential_0 = 0.0f;
    DiffPair_float_0 _S6147;
    (&_S6147)->primal_0 = _S5727;
    (&_S6147)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6146, &_S6147, 0.0f);
    float3  _S6148 = make_float3 (_S6042.x, _S6042.y, _S6145.differential_0 + _S6147.differential_0 + _S6041.x + _S6041.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6149;
    (&_S6149)->primal_0 = R_37;
    (&_S6149)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6150;
    (&_S6150)->primal_0 = pos_i_6;
    (&_S6150)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6149, &_S6150, _S6148);
    DiffPair_float_0 _S6151;
    (&_S6151)->primal_0 = _S5720;
    (&_S6151)->differential_0 = 0.0f;
    DiffPair_float_0 _S6152;
    (&_S6152)->primal_0 = _S5722;
    (&_S6152)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6151, &_S6152, _S6144.differential_0);
    DiffPair_float_0 _S6153;
    (&_S6153)->primal_0 = _S5719;
    (&_S6153)->differential_0 = 0.0f;
    DiffPair_float_0 _S6154;
    (&_S6154)->primal_0 = _S5722;
    (&_S6154)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6153, &_S6154, _S6146.differential_0);
    float3  _S6155 = make_float3 (_S6056.x, _S6056.y, _S6152.differential_0 + _S6154.differential_0 + _S6055.x + _S6055.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6156;
    (&_S6156)->primal_0 = R_37;
    (&_S6156)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6157;
    (&_S6157)->primal_0 = pos_i_5;
    (&_S6157)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6156, &_S6157, _S6155);
    DiffPair_float_0 _S6158;
    (&_S6158)->primal_0 = _S5715;
    (&_S6158)->differential_0 = 0.0f;
    DiffPair_float_0 _S6159;
    (&_S6159)->primal_0 = _S5717;
    (&_S6159)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6158, &_S6159, _S6151.differential_0);
    DiffPair_float_0 _S6160;
    (&_S6160)->primal_0 = _S5714;
    (&_S6160)->differential_0 = 0.0f;
    DiffPair_float_0 _S6161;
    (&_S6161)->primal_0 = _S5717;
    (&_S6161)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6160, &_S6161, _S6153.differential_0);
    float3  _S6162 = make_float3 (_S6070.x, _S6070.y, _S6159.differential_0 + _S6161.differential_0 + _S6069.x + _S6069.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6163;
    (&_S6163)->primal_0 = R_37;
    (&_S6163)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6164;
    (&_S6164)->primal_0 = pos_i_4;
    (&_S6164)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6163, &_S6164, _S6162);
    DiffPair_float_0 _S6165;
    (&_S6165)->primal_0 = _S5710;
    (&_S6165)->differential_0 = 0.0f;
    DiffPair_float_0 _S6166;
    (&_S6166)->primal_0 = _S5712;
    (&_S6166)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6165, &_S6166, _S6158.differential_0);
    DiffPair_float_0 _S6167;
    (&_S6167)->primal_0 = _S5709;
    (&_S6167)->differential_0 = 0.0f;
    DiffPair_float_0 _S6168;
    (&_S6168)->primal_0 = _S5712;
    (&_S6168)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6167, &_S6168, _S6160.differential_0);
    float3  _S6169 = make_float3 (_S6084.x, _S6084.y, _S6166.differential_0 + _S6168.differential_0 + _S6083.x + _S6083.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6170;
    (&_S6170)->primal_0 = R_37;
    (&_S6170)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6171;
    (&_S6171)->primal_0 = pos_i_3;
    (&_S6171)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6170, &_S6171, _S6169);
    DiffPair_float_0 _S6172;
    (&_S6172)->primal_0 = _S5705;
    (&_S6172)->differential_0 = 0.0f;
    DiffPair_float_0 _S6173;
    (&_S6173)->primal_0 = _S5707;
    (&_S6173)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6172, &_S6173, _S6165.differential_0);
    DiffPair_float_0 _S6174;
    (&_S6174)->primal_0 = _S5704;
    (&_S6174)->differential_0 = 0.0f;
    DiffPair_float_0 _S6175;
    (&_S6175)->primal_0 = _S5707;
    (&_S6175)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6174, &_S6175, _S6167.differential_0);
    float3  _S6176 = make_float3 (_S6098.x, _S6098.y, _S6173.differential_0 + _S6175.differential_0 + _S6097.x + _S6097.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6177;
    (&_S6177)->primal_0 = R_37;
    (&_S6177)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6178;
    (&_S6178)->primal_0 = pos_i_2;
    (&_S6178)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6177, &_S6178, _S6176);
    DiffPair_float_0 _S6179;
    (&_S6179)->primal_0 = _S5700;
    (&_S6179)->differential_0 = 0.0f;
    DiffPair_float_0 _S6180;
    (&_S6180)->primal_0 = _S5702;
    (&_S6180)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6179, &_S6180, _S6172.differential_0);
    DiffPair_float_0 _S6181;
    (&_S6181)->primal_0 = _S5699;
    (&_S6181)->differential_0 = 0.0f;
    DiffPair_float_0 _S6182;
    (&_S6182)->primal_0 = _S5702;
    (&_S6182)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6181, &_S6182, _S6174.differential_0);
    float3  _S6183 = make_float3 (_S6112.x, _S6112.y, _S6180.differential_0 + _S6182.differential_0 + _S6111.x + _S6111.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6184;
    (&_S6184)->primal_0 = R_37;
    (&_S6184)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6185;
    (&_S6185)->primal_0 = pos_i_1;
    (&_S6185)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6184, &_S6185, _S6183);
    DiffPair_float_0 _S6186;
    (&_S6186)->primal_0 = _S5695;
    (&_S6186)->differential_0 = 0.0f;
    DiffPair_float_0 _S6187;
    (&_S6187)->primal_0 = _S5697;
    (&_S6187)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6186, &_S6187, _S6179.differential_0);
    DiffPair_float_0 _S6188;
    (&_S6188)->primal_0 = _S5694;
    (&_S6188)->differential_0 = 0.0f;
    DiffPair_float_0 _S6189;
    (&_S6189)->primal_0 = _S5697;
    (&_S6189)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6188, &_S6189, _S6181.differential_0);
    float3  _S6190 = make_float3 (_S6126.x, _S6126.y, _S6187.differential_0 + _S6189.differential_0 + _S6125.x + _S6125.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6191;
    (&_S6191)->primal_0 = R_37;
    (&_S6191)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6192;
    (&_S6192)->primal_0 = pos_i_0;
    (&_S6192)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6191, &_S6192, _S6190);
    DiffPair_float_0 _S6193;
    (&_S6193)->primal_0 = 0.0f;
    (&_S6193)->differential_0 = 0.0f;
    DiffPair_float_0 _S6194;
    (&_S6194)->primal_0 = _S5692;
    (&_S6194)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6193, &_S6194, _S6186.differential_0);
    DiffPair_float_0 _S6195;
    (&_S6195)->primal_0 = 1.00000001504746622e+30f;
    (&_S6195)->differential_0 = 0.0f;
    DiffPair_float_0 _S6196;
    (&_S6196)->primal_0 = _S5692;
    (&_S6196)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6195, &_S6196, _S6188.differential_0);
    float3  _S6197 = make_float3 (_S6140.x, _S6140.y, _S6194.differential_0 + _S6196.differential_0 + _S6139.x + _S6139.y);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6198;
    (&_S6198)->primal_0 = R_37;
    (&_S6198)->differential_0 = _S5966;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6199;
    (&_S6199)->primal_0 = pos_4;
    (&_S6199)->differential_0 = _S5904;
    s_bwd_prop_mul_1(&_S6198, &_S6199, _S6197);
    float3  _S6200 = _S5968.differential_0 + _S6141 + _S6148 + _S6155 + _S6162 + _S6169 + _S6176 + _S6183 + _S6190 + _S6197;
    Matrix<float, 3, 3>  _S6201 = _S5969 + _S6142.differential_0 + _S6149.differential_0 + _S6156.differential_0 + _S6163.differential_0 + _S6170.differential_0 + _S6177.differential_0 + _S6184.differential_0 + _S6191.differential_0 + _S6198.differential_0;
    FixedArray<float3 , 16>  _S6202;
    _S6202[int(0)] = _S5904;
    _S6202[int(1)] = _S5904;
    _S6202[int(2)] = _S5904;
    _S6202[int(3)] = _S5904;
    _S6202[int(4)] = _S5904;
    _S6202[int(5)] = _S5904;
    _S6202[int(6)] = _S5904;
    _S6202[int(7)] = _S5904;
    _S6202[int(8)] = _S5904;
    _S6202[int(9)] = _S5904;
    _S6202[int(10)] = _S5904;
    _S6202[int(11)] = _S5904;
    _S6202[int(12)] = _S5904;
    _S6202[int(13)] = _S5904;
    _S6202[int(14)] = _S5904;
    _S6202[int(15)] = _S5904;
    _S6202[int(15)] = _S5907;
    _S6202[int(14)] = _S5909;
    _S6202[int(13)] = _S5911;
    _S6202[int(12)] = _S5913;
    _S6202[int(11)] = _S5915;
    _S6202[int(10)] = _S5917;
    _S6202[int(9)] = _S5919;
    _S6202[int(8)] = _S5927;
    _S6202[int(7)] = _S5929;
    _S6202[int(6)] = _S5931;
    _S6202[int(5)] = _S5933;
    _S6202[int(4)] = _S5935;
    _S6202[int(3)] = _S5946;
    _S6202[int(2)] = _S5948;
    _S6202[int(1)] = _S5950;
    _S6202[int(0)] = _S5963;
    (*v_densities_0)[int(0)] = 0.0f;
    (*v_densities_0)[int(1)] = 0.0f;
    (*v_densities_0)[int(2)] = 0.0f;
    (*v_densities_0)[int(3)] = 0.0f;
    (*v_densities_0)[int(4)] = 0.0f;
    (*v_densities_0)[int(5)] = 0.0f;
    (*v_densities_0)[int(6)] = 0.0f;
    (*v_densities_0)[int(7)] = 0.0f;
    *v_sh_coeffs_9 = _S6202;
    *v_R_11 = _S6201;
    *v_t_10 = _S6200;
    return;
}

inline __device__ void projection_voxel_eval3d_fisheye_vjp(float3  pos_5, float size_5, FixedArray<float, 8>  * densities_5, FixedArray<float3 , 16>  * sh_coeffs_32, Matrix<float, 3, 3>  R_38, float3  t_38, float fx_41, float fy_41, float cx_36, float cy_36, FixedArray<float, 10>  * dist_coeffs_45, uint image_width_32, uint image_height_32, float3  v_rgb_9, FixedArray<float, 8>  * v_densities_1, FixedArray<float3 , 16>  * v_sh_coeffs_10, Matrix<float, 3, 3>  * v_R_12, float3  * v_t_11)
{
    float3  _S6203 = s_primal_ctx_mul_0(R_38, pos_5) + t_38;
    float _S6204 = length_1(_S6203);
    float _S6205 = s_primal_ctx_min_0(1.00000001504746622e+30f, _S6204);
    float _S6206 = s_primal_ctx_max_0(0.0f, _S6204);
    float3  pos_i_7 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 0.0f);
    float3  _S6207 = s_primal_ctx_mul_0(R_38, pos_i_7) + t_38;
    float _S6208 = length_1(_S6207);
    float _S6209 = s_primal_ctx_min_0(_S6205, _S6208);
    float _S6210 = s_primal_ctx_max_0(_S6206, _S6208);
    float3  pos_i_8 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 0.0f);
    float3  _S6211 = s_primal_ctx_mul_0(R_38, pos_i_8) + t_38;
    float _S6212 = length_1(_S6211);
    float _S6213 = s_primal_ctx_min_0(_S6209, _S6212);
    float _S6214 = s_primal_ctx_max_0(_S6210, _S6212);
    float3  pos_i_9 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 1.0f, 0.0f);
    float3  _S6215 = s_primal_ctx_mul_0(R_38, pos_i_9) + t_38;
    float _S6216 = length_1(_S6215);
    float _S6217 = s_primal_ctx_min_0(_S6213, _S6216);
    float _S6218 = s_primal_ctx_max_0(_S6214, _S6216);
    float3  pos_i_10 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 0.0f, 1.0f);
    float3  _S6219 = s_primal_ctx_mul_0(R_38, pos_i_10) + t_38;
    float _S6220 = length_1(_S6219);
    float _S6221 = s_primal_ctx_min_0(_S6217, _S6220);
    float _S6222 = s_primal_ctx_max_0(_S6218, _S6220);
    float3  pos_i_11 = pos_5 + make_float3 (size_5) * make_float3 (1.0f, 0.0f, 1.0f);
    float3  _S6223 = s_primal_ctx_mul_0(R_38, pos_i_11) + t_38;
    float _S6224 = length_1(_S6223);
    float _S6225 = s_primal_ctx_min_0(_S6221, _S6224);
    float _S6226 = s_primal_ctx_max_0(_S6222, _S6224);
    float3  pos_i_12 = pos_5 + make_float3 (size_5) * make_float3 (0.0f, 1.0f, 1.0f);
    float3  _S6227 = s_primal_ctx_mul_0(R_38, pos_i_12) + t_38;
    float _S6228 = length_1(_S6227);
    float _S6229 = s_primal_ctx_min_0(_S6225, _S6228);
    float _S6230 = s_primal_ctx_max_0(_S6226, _S6228);
    float3  pos_i_13 = pos_5 + make_float3 (size_5);
    float3  _S6231 = s_primal_ctx_mul_0(R_38, pos_i_13) + t_38;
    float _S6232 = length_1(_S6231);
    float3  _S6233 = pos_5 + make_float3 (0.5f * size_5);
    float3  mean_c_32 = s_primal_ctx_mul_0(R_38, _S6233) + t_38;
    float2  _S6234 = float2 {_S6203.x, _S6203.y};
    float _S6235 = length_0(_S6234);
    float _S6236 = _S6203.z;
    float _S6237 = s_primal_ctx_atan2_0(_S6235, _S6236);
    bool _S6238 = _S6237 < 0.00100000004749745f;
    float k_15;
    float _S6239;
    float _S6240;
    float _S6241;
    if(_S6238)
    {
        float _S6242 = 1.0f - _S6237 * _S6237 / 3.0f;
        float _S6243 = _S6236 * _S6236;
        k_15 = _S6242 / _S6236;
        _S6239 = 0.0f;
        _S6240 = _S6243;
        _S6241 = _S6242;
    }
    else
    {
        float _S6244 = _S6235 * _S6235;
        k_15 = _S6237 / _S6235;
        _S6239 = _S6244;
        _S6240 = 0.0f;
        _S6241 = 0.0f;
    }
    float2  _S6245 = make_float2 (k_15);
    float2  _S6246 = _S6234 * make_float2 (k_15);
    float u_125 = _S6246.x;
    float v_125 = _S6246.y;
    float r2_120 = u_125 * u_125 + v_125 * v_125;
    float _S6247 = (*dist_coeffs_45)[int(2)] + r2_120 * (*dist_coeffs_45)[int(3)];
    float _S6248 = (*dist_coeffs_45)[int(1)] + r2_120 * _S6247;
    float _S6249 = (*dist_coeffs_45)[int(0)] + r2_120 * _S6248;
    float radial_21 = 1.0f + r2_120 * _S6249;
    float2  _S6250 = make_float2 (radial_21);
    float _S6251 = 2.0f * (*dist_coeffs_45)[int(4)];
    float _S6252 = _S6251 * u_125;
    float _S6253 = 2.0f * u_125;
    float _S6254 = 2.0f * (*dist_coeffs_45)[int(5)];
    float _S6255 = _S6254 * u_125;
    float _S6256 = 2.0f * v_125;
    float2  _S6257 = _S6246 * make_float2 (radial_21) + make_float2 (_S6252 * v_125 + (*dist_coeffs_45)[int(5)] * (r2_120 + _S6253 * u_125) + (*dist_coeffs_45)[int(6)] * r2_120, _S6255 * v_125 + (*dist_coeffs_45)[int(4)] * (r2_120 + _S6256 * v_125) + (*dist_coeffs_45)[int(7)] * r2_120);
    float2  _S6258 = _S6257 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6257.x + (*dist_coeffs_45)[int(9)] * _S6257.y, 0.0f);
    float _S6259 = fx_41 * _S6258.x + cx_36;
    float _S6260 = fy_41 * _S6258.y + cy_36;
    float2  _S6261 = float2 {_S6207.x, _S6207.y};
    float _S6262 = length_0(_S6261);
    float _S6263 = _S6207.z;
    float _S6264 = s_primal_ctx_atan2_0(_S6262, _S6263);
    bool _S6265 = _S6264 < 0.00100000004749745f;
    float _S6266;
    float _S6267;
    float _S6268;
    if(_S6265)
    {
        float _S6269 = 1.0f - _S6264 * _S6264 / 3.0f;
        float _S6270 = _S6263 * _S6263;
        k_15 = _S6269 / _S6263;
        _S6266 = 0.0f;
        _S6267 = _S6270;
        _S6268 = _S6269;
    }
    else
    {
        float _S6271 = _S6262 * _S6262;
        k_15 = _S6264 / _S6262;
        _S6266 = _S6271;
        _S6267 = 0.0f;
        _S6268 = 0.0f;
    }
    float2  _S6272 = make_float2 (k_15);
    float2  _S6273 = _S6261 * make_float2 (k_15);
    float u_126 = _S6273.x;
    float v_126 = _S6273.y;
    float r2_121 = u_126 * u_126 + v_126 * v_126;
    float _S6274 = (*dist_coeffs_45)[int(2)] + r2_121 * (*dist_coeffs_45)[int(3)];
    float _S6275 = (*dist_coeffs_45)[int(1)] + r2_121 * _S6274;
    float _S6276 = (*dist_coeffs_45)[int(0)] + r2_121 * _S6275;
    float radial_22 = 1.0f + r2_121 * _S6276;
    float2  _S6277 = make_float2 (radial_22);
    float _S6278 = _S6251 * u_126;
    float _S6279 = 2.0f * u_126;
    float _S6280 = _S6254 * u_126;
    float _S6281 = 2.0f * v_126;
    float2  _S6282 = _S6273 * make_float2 (radial_22) + make_float2 (_S6278 * v_126 + (*dist_coeffs_45)[int(5)] * (r2_121 + _S6279 * u_126) + (*dist_coeffs_45)[int(6)] * r2_121, _S6280 * v_126 + (*dist_coeffs_45)[int(4)] * (r2_121 + _S6281 * v_126) + (*dist_coeffs_45)[int(7)] * r2_121);
    float2  _S6283 = _S6282 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6282.x + (*dist_coeffs_45)[int(9)] * _S6282.y, 0.0f);
    float _S6284 = fx_41 * _S6283.x + cx_36;
    float _S6285 = fy_41 * _S6283.y + cy_36;
    float2  _S6286 = float2 {_S6211.x, _S6211.y};
    float _S6287 = length_0(_S6286);
    float _S6288 = _S6211.z;
    float _S6289 = s_primal_ctx_atan2_0(_S6287, _S6288);
    bool _S6290 = _S6289 < 0.00100000004749745f;
    float _S6291;
    float _S6292;
    float _S6293;
    if(_S6290)
    {
        float _S6294 = 1.0f - _S6289 * _S6289 / 3.0f;
        float _S6295 = _S6288 * _S6288;
        k_15 = _S6294 / _S6288;
        _S6291 = 0.0f;
        _S6292 = _S6295;
        _S6293 = _S6294;
    }
    else
    {
        float _S6296 = _S6287 * _S6287;
        k_15 = _S6289 / _S6287;
        _S6291 = _S6296;
        _S6292 = 0.0f;
        _S6293 = 0.0f;
    }
    float2  _S6297 = make_float2 (k_15);
    float2  _S6298 = _S6286 * make_float2 (k_15);
    float u_127 = _S6298.x;
    float v_127 = _S6298.y;
    float r2_122 = u_127 * u_127 + v_127 * v_127;
    float _S6299 = (*dist_coeffs_45)[int(2)] + r2_122 * (*dist_coeffs_45)[int(3)];
    float _S6300 = (*dist_coeffs_45)[int(1)] + r2_122 * _S6299;
    float _S6301 = (*dist_coeffs_45)[int(0)] + r2_122 * _S6300;
    float radial_23 = 1.0f + r2_122 * _S6301;
    float2  _S6302 = make_float2 (radial_23);
    float _S6303 = _S6251 * u_127;
    float _S6304 = 2.0f * u_127;
    float _S6305 = _S6254 * u_127;
    float _S6306 = 2.0f * v_127;
    float2  _S6307 = _S6298 * make_float2 (radial_23) + make_float2 (_S6303 * v_127 + (*dist_coeffs_45)[int(5)] * (r2_122 + _S6304 * u_127) + (*dist_coeffs_45)[int(6)] * r2_122, _S6305 * v_127 + (*dist_coeffs_45)[int(4)] * (r2_122 + _S6306 * v_127) + (*dist_coeffs_45)[int(7)] * r2_122);
    float2  _S6308 = _S6307 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6307.x + (*dist_coeffs_45)[int(9)] * _S6307.y, 0.0f);
    float _S6309 = fx_41 * _S6308.x + cx_36;
    float _S6310 = fy_41 * _S6308.y + cy_36;
    float2  _S6311 = float2 {_S6215.x, _S6215.y};
    float _S6312 = length_0(_S6311);
    float _S6313 = _S6215.z;
    float _S6314 = s_primal_ctx_atan2_0(_S6312, _S6313);
    bool _S6315 = _S6314 < 0.00100000004749745f;
    float _S6316;
    float _S6317;
    float _S6318;
    if(_S6315)
    {
        float _S6319 = 1.0f - _S6314 * _S6314 / 3.0f;
        float _S6320 = _S6313 * _S6313;
        k_15 = _S6319 / _S6313;
        _S6316 = 0.0f;
        _S6317 = _S6320;
        _S6318 = _S6319;
    }
    else
    {
        float _S6321 = _S6312 * _S6312;
        k_15 = _S6314 / _S6312;
        _S6316 = _S6321;
        _S6317 = 0.0f;
        _S6318 = 0.0f;
    }
    float2  _S6322 = make_float2 (k_15);
    float2  _S6323 = _S6311 * make_float2 (k_15);
    float u_128 = _S6323.x;
    float v_128 = _S6323.y;
    float r2_123 = u_128 * u_128 + v_128 * v_128;
    float _S6324 = (*dist_coeffs_45)[int(2)] + r2_123 * (*dist_coeffs_45)[int(3)];
    float _S6325 = (*dist_coeffs_45)[int(1)] + r2_123 * _S6324;
    float _S6326 = (*dist_coeffs_45)[int(0)] + r2_123 * _S6325;
    float radial_24 = 1.0f + r2_123 * _S6326;
    float2  _S6327 = make_float2 (radial_24);
    float _S6328 = _S6251 * u_128;
    float _S6329 = 2.0f * u_128;
    float _S6330 = _S6254 * u_128;
    float _S6331 = 2.0f * v_128;
    float2  _S6332 = _S6323 * make_float2 (radial_24) + make_float2 (_S6328 * v_128 + (*dist_coeffs_45)[int(5)] * (r2_123 + _S6329 * u_128) + (*dist_coeffs_45)[int(6)] * r2_123, _S6330 * v_128 + (*dist_coeffs_45)[int(4)] * (r2_123 + _S6331 * v_128) + (*dist_coeffs_45)[int(7)] * r2_123);
    float2  _S6333 = _S6332 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6332.x + (*dist_coeffs_45)[int(9)] * _S6332.y, 0.0f);
    float _S6334 = fx_41 * _S6333.x + cx_36;
    float _S6335 = fy_41 * _S6333.y + cy_36;
    float2  _S6336 = float2 {_S6219.x, _S6219.y};
    float _S6337 = length_0(_S6336);
    float _S6338 = _S6219.z;
    float _S6339 = s_primal_ctx_atan2_0(_S6337, _S6338);
    bool _S6340 = _S6339 < 0.00100000004749745f;
    float _S6341;
    float _S6342;
    float _S6343;
    if(_S6340)
    {
        float _S6344 = 1.0f - _S6339 * _S6339 / 3.0f;
        float _S6345 = _S6338 * _S6338;
        k_15 = _S6344 / _S6338;
        _S6341 = 0.0f;
        _S6342 = _S6345;
        _S6343 = _S6344;
    }
    else
    {
        float _S6346 = _S6337 * _S6337;
        k_15 = _S6339 / _S6337;
        _S6341 = _S6346;
        _S6342 = 0.0f;
        _S6343 = 0.0f;
    }
    float2  _S6347 = make_float2 (k_15);
    float2  _S6348 = _S6336 * make_float2 (k_15);
    float u_129 = _S6348.x;
    float v_129 = _S6348.y;
    float r2_124 = u_129 * u_129 + v_129 * v_129;
    float _S6349 = (*dist_coeffs_45)[int(2)] + r2_124 * (*dist_coeffs_45)[int(3)];
    float _S6350 = (*dist_coeffs_45)[int(1)] + r2_124 * _S6349;
    float _S6351 = (*dist_coeffs_45)[int(0)] + r2_124 * _S6350;
    float radial_25 = 1.0f + r2_124 * _S6351;
    float2  _S6352 = make_float2 (radial_25);
    float _S6353 = _S6251 * u_129;
    float _S6354 = 2.0f * u_129;
    float _S6355 = _S6254 * u_129;
    float _S6356 = 2.0f * v_129;
    float2  _S6357 = _S6348 * make_float2 (radial_25) + make_float2 (_S6353 * v_129 + (*dist_coeffs_45)[int(5)] * (r2_124 + _S6354 * u_129) + (*dist_coeffs_45)[int(6)] * r2_124, _S6355 * v_129 + (*dist_coeffs_45)[int(4)] * (r2_124 + _S6356 * v_129) + (*dist_coeffs_45)[int(7)] * r2_124);
    float2  _S6358 = _S6357 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6357.x + (*dist_coeffs_45)[int(9)] * _S6357.y, 0.0f);
    float _S6359 = fx_41 * _S6358.x + cx_36;
    float _S6360 = fy_41 * _S6358.y + cy_36;
    float2  _S6361 = float2 {_S6223.x, _S6223.y};
    float _S6362 = length_0(_S6361);
    float _S6363 = _S6223.z;
    float _S6364 = s_primal_ctx_atan2_0(_S6362, _S6363);
    bool _S6365 = _S6364 < 0.00100000004749745f;
    float _S6366;
    float _S6367;
    float _S6368;
    if(_S6365)
    {
        float _S6369 = 1.0f - _S6364 * _S6364 / 3.0f;
        float _S6370 = _S6363 * _S6363;
        k_15 = _S6369 / _S6363;
        _S6366 = 0.0f;
        _S6367 = _S6370;
        _S6368 = _S6369;
    }
    else
    {
        float _S6371 = _S6362 * _S6362;
        k_15 = _S6364 / _S6362;
        _S6366 = _S6371;
        _S6367 = 0.0f;
        _S6368 = 0.0f;
    }
    float2  _S6372 = make_float2 (k_15);
    float2  _S6373 = _S6361 * make_float2 (k_15);
    float u_130 = _S6373.x;
    float v_130 = _S6373.y;
    float r2_125 = u_130 * u_130 + v_130 * v_130;
    float _S6374 = (*dist_coeffs_45)[int(2)] + r2_125 * (*dist_coeffs_45)[int(3)];
    float _S6375 = (*dist_coeffs_45)[int(1)] + r2_125 * _S6374;
    float _S6376 = (*dist_coeffs_45)[int(0)] + r2_125 * _S6375;
    float radial_26 = 1.0f + r2_125 * _S6376;
    float2  _S6377 = make_float2 (radial_26);
    float _S6378 = _S6251 * u_130;
    float _S6379 = 2.0f * u_130;
    float _S6380 = _S6254 * u_130;
    float _S6381 = 2.0f * v_130;
    float2  _S6382 = _S6373 * make_float2 (radial_26) + make_float2 (_S6378 * v_130 + (*dist_coeffs_45)[int(5)] * (r2_125 + _S6379 * u_130) + (*dist_coeffs_45)[int(6)] * r2_125, _S6380 * v_130 + (*dist_coeffs_45)[int(4)] * (r2_125 + _S6381 * v_130) + (*dist_coeffs_45)[int(7)] * r2_125);
    float2  _S6383 = _S6382 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6382.x + (*dist_coeffs_45)[int(9)] * _S6382.y, 0.0f);
    float _S6384 = fx_41 * _S6383.x + cx_36;
    float _S6385 = fy_41 * _S6383.y + cy_36;
    float2  _S6386 = float2 {_S6227.x, _S6227.y};
    float _S6387 = length_0(_S6386);
    float _S6388 = _S6227.z;
    float _S6389 = s_primal_ctx_atan2_0(_S6387, _S6388);
    bool _S6390 = _S6389 < 0.00100000004749745f;
    float _S6391;
    float _S6392;
    float _S6393;
    if(_S6390)
    {
        float _S6394 = 1.0f - _S6389 * _S6389 / 3.0f;
        float _S6395 = _S6388 * _S6388;
        k_15 = _S6394 / _S6388;
        _S6391 = 0.0f;
        _S6392 = _S6395;
        _S6393 = _S6394;
    }
    else
    {
        float _S6396 = _S6387 * _S6387;
        k_15 = _S6389 / _S6387;
        _S6391 = _S6396;
        _S6392 = 0.0f;
        _S6393 = 0.0f;
    }
    float2  _S6397 = make_float2 (k_15);
    float2  _S6398 = _S6386 * make_float2 (k_15);
    float u_131 = _S6398.x;
    float v_131 = _S6398.y;
    float r2_126 = u_131 * u_131 + v_131 * v_131;
    float _S6399 = (*dist_coeffs_45)[int(2)] + r2_126 * (*dist_coeffs_45)[int(3)];
    float _S6400 = (*dist_coeffs_45)[int(1)] + r2_126 * _S6399;
    float _S6401 = (*dist_coeffs_45)[int(0)] + r2_126 * _S6400;
    float radial_27 = 1.0f + r2_126 * _S6401;
    float2  _S6402 = make_float2 (radial_27);
    float _S6403 = _S6251 * u_131;
    float _S6404 = 2.0f * u_131;
    float _S6405 = _S6254 * u_131;
    float _S6406 = 2.0f * v_131;
    float2  _S6407 = _S6398 * make_float2 (radial_27) + make_float2 (_S6403 * v_131 + (*dist_coeffs_45)[int(5)] * (r2_126 + _S6404 * u_131) + (*dist_coeffs_45)[int(6)] * r2_126, _S6405 * v_131 + (*dist_coeffs_45)[int(4)] * (r2_126 + _S6406 * v_131) + (*dist_coeffs_45)[int(7)] * r2_126);
    float2  _S6408 = _S6407 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6407.x + (*dist_coeffs_45)[int(9)] * _S6407.y, 0.0f);
    float _S6409 = fx_41 * _S6408.x + cx_36;
    float _S6410 = fy_41 * _S6408.y + cy_36;
    float2  _S6411 = float2 {_S6231.x, _S6231.y};
    float _S6412 = length_0(_S6411);
    float _S6413 = _S6231.z;
    float _S6414 = s_primal_ctx_atan2_0(_S6412, _S6413);
    bool _S6415 = _S6414 < 0.00100000004749745f;
    float _S6416;
    float _S6417;
    float _S6418;
    if(_S6415)
    {
        float _S6419 = 1.0f - _S6414 * _S6414 / 3.0f;
        float _S6420 = _S6413 * _S6413;
        k_15 = _S6419 / _S6413;
        _S6416 = 0.0f;
        _S6417 = _S6420;
        _S6418 = _S6419;
    }
    else
    {
        float _S6421 = _S6412 * _S6412;
        k_15 = _S6414 / _S6412;
        _S6416 = _S6421;
        _S6417 = 0.0f;
        _S6418 = 0.0f;
    }
    float2  _S6422 = make_float2 (k_15);
    float2  _S6423 = _S6411 * make_float2 (k_15);
    float u_132 = _S6423.x;
    float v_132 = _S6423.y;
    float r2_127 = u_132 * u_132 + v_132 * v_132;
    float _S6424 = (*dist_coeffs_45)[int(2)] + r2_127 * (*dist_coeffs_45)[int(3)];
    float _S6425 = (*dist_coeffs_45)[int(1)] + r2_127 * _S6424;
    float _S6426 = (*dist_coeffs_45)[int(0)] + r2_127 * _S6425;
    float radial_28 = 1.0f + r2_127 * _S6426;
    float2  _S6427 = make_float2 (radial_28);
    float _S6428 = _S6251 * u_132;
    float _S6429 = 2.0f * u_132;
    float _S6430 = _S6254 * u_132;
    float _S6431 = 2.0f * v_132;
    float2  _S6432 = _S6423 * make_float2 (radial_28) + make_float2 (_S6428 * v_132 + (*dist_coeffs_45)[int(5)] * (r2_127 + _S6429 * u_132) + (*dist_coeffs_45)[int(6)] * r2_127, _S6430 * v_132 + (*dist_coeffs_45)[int(4)] * (r2_127 + _S6431 * v_132) + (*dist_coeffs_45)[int(7)] * r2_127);
    float2  _S6433 = _S6432 + make_float2 ((*dist_coeffs_45)[int(8)] * _S6432.x + (*dist_coeffs_45)[int(9)] * _S6432.y, 0.0f);
    float _S6434 = fx_41 * _S6433.x + cx_36;
    float _S6435 = fy_41 * _S6433.y + cy_36;
    float _S6436 = s_primal_ctx_max_0(_S6259, _S6284);
    float _S6437 = s_primal_ctx_min_0(_S6259, _S6284);
    float _S6438 = s_primal_ctx_max_0(_S6260, _S6285);
    float _S6439 = s_primal_ctx_min_0(_S6260, _S6285);
    float _S6440 = s_primal_ctx_max_0(_S6436, _S6309);
    float _S6441 = s_primal_ctx_min_0(_S6437, _S6309);
    float _S6442 = s_primal_ctx_max_0(_S6438, _S6310);
    float _S6443 = s_primal_ctx_min_0(_S6439, _S6310);
    float _S6444 = s_primal_ctx_max_0(_S6440, _S6334);
    float _S6445 = s_primal_ctx_min_0(_S6441, _S6334);
    float _S6446 = s_primal_ctx_max_0(_S6442, _S6335);
    float _S6447 = s_primal_ctx_min_0(_S6443, _S6335);
    float _S6448 = s_primal_ctx_max_0(_S6444, _S6359);
    float _S6449 = s_primal_ctx_min_0(_S6445, _S6359);
    float _S6450 = s_primal_ctx_max_0(_S6446, _S6360);
    float _S6451 = s_primal_ctx_min_0(_S6447, _S6360);
    float _S6452 = s_primal_ctx_max_0(_S6448, _S6384);
    float _S6453 = s_primal_ctx_min_0(_S6449, _S6384);
    float _S6454 = s_primal_ctx_max_0(_S6450, _S6385);
    float _S6455 = s_primal_ctx_min_0(_S6451, _S6385);
    float _S6456 = s_primal_ctx_max_0(_S6452, _S6409);
    float _S6457 = s_primal_ctx_min_0(_S6453, _S6409);
    float _S6458 = s_primal_ctx_max_0(_S6454, _S6410);
    float _S6459 = s_primal_ctx_min_0(_S6455, _S6410);
    float _S6460 = s_primal_ctx_dot_0(mean_c_32, mean_c_32) + 9.99999997475242708e-07f;
    Matrix<float, 3, 3>  _S6461 = transpose_0(R_38);
    float3  _S6462 = mean_c_32 - - s_primal_ctx_mul_0(_S6461, t_38);
    float _S6463 = _S6462.x;
    float _S6464 = _S6462.y;
    float _S6465 = _S6462.z;
    float _S6466 = _S6463 * _S6463 + _S6464 * _S6464 + _S6465 * _S6465;
    float _S6467 = s_primal_ctx_sqrt_0(_S6466);
    float x_68 = _S6463 / _S6467;
    float3  _S6468 = make_float3 (x_68);
    float _S6469 = _S6467 * _S6467;
    float y_35 = _S6464 / _S6467;
    float z_32 = _S6465 / _S6467;
    float3  _S6470 = make_float3 (z_32);
    float _S6471 = - y_35;
    float3  _S6472 = make_float3 (_S6471);
    float z2_65 = z_32 * z_32;
    float fTmp0B_32 = -1.09254848957061768f * z_32;
    float fC1_32 = x_68 * x_68 - y_35 * y_35;
    float _S6473 = 2.0f * x_68;
    float fS1_32 = _S6473 * y_35;
    float pSH6_10 = 0.94617468118667603f * z2_65 - 0.31539157032966614f;
    float3  _S6474 = make_float3 (pSH6_10);
    float pSH7_10 = fTmp0B_32 * x_68;
    float3  _S6475 = make_float3 (pSH7_10);
    float pSH5_10 = fTmp0B_32 * y_35;
    float3  _S6476 = make_float3 (pSH5_10);
    float pSH8_10 = 0.54627424478530884f * fC1_32;
    float3  _S6477 = make_float3 (pSH8_10);
    float pSH4_10 = 0.54627424478530884f * fS1_32;
    float3  _S6478 = make_float3 (pSH4_10);
    float fTmp0C_32 = -2.28522896766662598f * z2_65 + 0.4570457935333252f;
    float fTmp1B_32 = 1.44530570507049561f * z_32;
    float _S6479 = 1.86588168144226074f * z2_65 - 1.11952900886535645f;
    float pSH12_10 = z_32 * _S6479;
    float3  _S6480 = make_float3 (pSH12_10);
    float pSH13_10 = fTmp0C_32 * x_68;
    float3  _S6481 = make_float3 (pSH13_10);
    float pSH11_10 = fTmp0C_32 * y_35;
    float3  _S6482 = make_float3 (pSH11_10);
    float pSH14_10 = fTmp1B_32 * fC1_32;
    float3  _S6483 = make_float3 (pSH14_10);
    float pSH10_10 = fTmp1B_32 * fS1_32;
    float3  _S6484 = make_float3 (pSH10_10);
    float pSH15_10 = -0.59004360437393188f * (x_68 * fC1_32 - y_35 * fS1_32);
    float3  _S6485 = make_float3 (pSH15_10);
    float pSH9_10 = -0.59004360437393188f * (x_68 * fS1_32 + y_35 * fC1_32);
    float3  _S6486 = make_float3 (pSH9_10);
    float3  _S6487 = make_float3 (0.0f);
    float3  _S6488 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6489;
    (&_S6489)->primal_0 = make_float3 (0.282094806432724f) * (*sh_coeffs_32)[int(0)] + make_float3 (0.48860251903533936f) * (make_float3 (_S6471) * (*sh_coeffs_32)[int(1)] + make_float3 (z_32) * (*sh_coeffs_32)[int(2)] - make_float3 (x_68) * (*sh_coeffs_32)[int(3)]) + (make_float3 (pSH4_10) * (*sh_coeffs_32)[int(4)] + make_float3 (pSH5_10) * (*sh_coeffs_32)[int(5)] + make_float3 (pSH6_10) * (*sh_coeffs_32)[int(6)] + make_float3 (pSH7_10) * (*sh_coeffs_32)[int(7)] + make_float3 (pSH8_10) * (*sh_coeffs_32)[int(8)]) + (make_float3 (pSH9_10) * (*sh_coeffs_32)[int(9)] + make_float3 (pSH10_10) * (*sh_coeffs_32)[int(10)] + make_float3 (pSH11_10) * (*sh_coeffs_32)[int(11)] + make_float3 (pSH12_10) * (*sh_coeffs_32)[int(12)] + make_float3 (pSH13_10) * (*sh_coeffs_32)[int(13)] + make_float3 (pSH14_10) * (*sh_coeffs_32)[int(14)] + make_float3 (pSH15_10) * (*sh_coeffs_32)[int(15)]) + make_float3 (0.5f);
    (&_S6489)->differential_0 = _S6488;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6490;
    (&_S6490)->primal_0 = _S6487;
    (&_S6490)->differential_0 = _S6488;
    s_bwd_prop_max_0(&_S6489, &_S6490, v_rgb_9);
    float3  _S6491 = _S6485 * _S6489.differential_0;
    float3  _S6492 = (*sh_coeffs_32)[int(15)] * _S6489.differential_0;
    float3  _S6493 = _S6483 * _S6489.differential_0;
    float3  _S6494 = (*sh_coeffs_32)[int(14)] * _S6489.differential_0;
    float3  _S6495 = _S6481 * _S6489.differential_0;
    float3  _S6496 = (*sh_coeffs_32)[int(13)] * _S6489.differential_0;
    float3  _S6497 = _S6480 * _S6489.differential_0;
    float3  _S6498 = (*sh_coeffs_32)[int(12)] * _S6489.differential_0;
    float3  _S6499 = _S6482 * _S6489.differential_0;
    float3  _S6500 = (*sh_coeffs_32)[int(11)] * _S6489.differential_0;
    float3  _S6501 = _S6484 * _S6489.differential_0;
    float3  _S6502 = (*sh_coeffs_32)[int(10)] * _S6489.differential_0;
    float3  _S6503 = _S6486 * _S6489.differential_0;
    float3  _S6504 = (*sh_coeffs_32)[int(9)] * _S6489.differential_0;
    float s_diff_fS2_T_10 = -0.59004360437393188f * (_S6504.x + _S6504.y + _S6504.z);
    float s_diff_fC2_T_10 = -0.59004360437393188f * (_S6492.x + _S6492.y + _S6492.z);
    float _S6505 = _S6502.x + _S6502.y + _S6502.z;
    float _S6506 = _S6494.x + _S6494.y + _S6494.z;
    float _S6507 = _S6500.x + _S6500.y + _S6500.z;
    float _S6508 = _S6496.x + _S6496.y + _S6496.z;
    float _S6509 = _S6498.x + _S6498.y + _S6498.z;
    float _S6510 = - s_diff_fC2_T_10;
    float3  _S6511 = _S6477 * _S6489.differential_0;
    float3  _S6512 = (*sh_coeffs_32)[int(8)] * _S6489.differential_0;
    float3  _S6513 = _S6475 * _S6489.differential_0;
    float3  _S6514 = (*sh_coeffs_32)[int(7)] * _S6489.differential_0;
    float3  _S6515 = _S6474 * _S6489.differential_0;
    float3  _S6516 = (*sh_coeffs_32)[int(6)] * _S6489.differential_0;
    float3  _S6517 = _S6476 * _S6489.differential_0;
    float3  _S6518 = (*sh_coeffs_32)[int(5)] * _S6489.differential_0;
    float3  _S6519 = _S6478 * _S6489.differential_0;
    float3  _S6520 = (*sh_coeffs_32)[int(4)] * _S6489.differential_0;
    float _S6521 = _S6518.x + _S6518.y + _S6518.z;
    float _S6522 = _S6514.x + _S6514.y + _S6514.z;
    float _S6523 = fTmp1B_32 * _S6505 + x_68 * s_diff_fS2_T_10 + y_35 * _S6510 + 0.54627424478530884f * (_S6520.x + _S6520.y + _S6520.z);
    float _S6524 = fTmp1B_32 * _S6506 + y_35 * s_diff_fS2_T_10 + x_68 * s_diff_fC2_T_10 + 0.54627424478530884f * (_S6512.x + _S6512.y + _S6512.z);
    float _S6525 = y_35 * - _S6524;
    float _S6526 = x_68 * _S6524;
    float _S6527 = z_32 * (1.86588168144226074f * (z_32 * _S6509) + -2.28522896766662598f * (y_35 * _S6507 + x_68 * _S6508) + 0.94617468118667603f * (_S6516.x + _S6516.y + _S6516.z));
    float3  _S6528 = make_float3 (0.48860251903533936f) * _S6489.differential_0;
    float3  _S6529 = - _S6528;
    float3  _S6530 = _S6468 * _S6529;
    float3  _S6531 = (*sh_coeffs_32)[int(3)] * _S6529;
    float3  _S6532 = _S6470 * _S6528;
    float3  _S6533 = (*sh_coeffs_32)[int(2)] * _S6528;
    float3  _S6534 = _S6472 * _S6528;
    float3  _S6535 = (*sh_coeffs_32)[int(1)] * _S6528;
    float _S6536 = (_S6479 * _S6509 + 1.44530570507049561f * (fS1_32 * _S6505 + fC1_32 * _S6506) + -1.09254848957061768f * (y_35 * _S6521 + x_68 * _S6522) + _S6527 + _S6527 + _S6533.x + _S6533.y + _S6533.z) / _S6469;
    float _S6537 = _S6467 * _S6536;
    float _S6538 = (fTmp0C_32 * _S6507 + fC1_32 * s_diff_fS2_T_10 + fS1_32 * _S6510 + fTmp0B_32 * _S6521 + _S6473 * _S6523 + _S6525 + _S6525 + - (_S6535.x + _S6535.y + _S6535.z)) / _S6469;
    float _S6539 = _S6467 * _S6538;
    float _S6540 = (fTmp0C_32 * _S6508 + fS1_32 * s_diff_fS2_T_10 + fC1_32 * s_diff_fC2_T_10 + fTmp0B_32 * _S6522 + 2.0f * (y_35 * _S6523) + _S6526 + _S6526 + _S6531.x + _S6531.y + _S6531.z) / _S6469;
    float _S6541 = _S6467 * _S6540;
    float _S6542 = _S6465 * - _S6536 + _S6464 * - _S6538 + _S6463 * - _S6540;
    DiffPair_float_0 _S6543;
    (&_S6543)->primal_0 = _S6466;
    (&_S6543)->differential_0 = 0.0f;
    s_bwd_prop_sqrt_0(&_S6543, _S6542);
    float _S6544 = _S6465 * _S6543.differential_0;
    float _S6545 = _S6464 * _S6543.differential_0;
    float _S6546 = _S6463 * _S6543.differential_0;
    float3  _S6547 = make_float3 (0.282094806432724f) * _S6489.differential_0;
    float3  _S6548 = make_float3 (_S6541 + _S6546 + _S6546, _S6539 + _S6545 + _S6545, _S6537 + _S6544 + _S6544);
    float3  _S6549 = - - _S6548;
    Matrix<float, 3, 3>  _S6550 = makeMatrix<float, 3, 3> (0.0f);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6551;
    (&_S6551)->primal_0 = _S6461;
    (&_S6551)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6552;
    (&_S6552)->primal_0 = t_38;
    (&_S6552)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6551, &_S6552, _S6549);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6553 = _S6552;
    Matrix<float, 3, 3>  _S6554 = transpose_0(_S6551.differential_0);
    DiffPair_float_0 _S6555;
    (&_S6555)->primal_0 = _S6460;
    (&_S6555)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S6555, 0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6556;
    (&_S6556)->primal_0 = mean_c_32;
    (&_S6556)->differential_0 = _S6488;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6557;
    (&_S6557)->primal_0 = mean_c_32;
    (&_S6557)->differential_0 = _S6488;
    s_bwd_prop_dot_0(&_S6556, &_S6557, _S6555.differential_0);
    DiffPair_float_0 _S6558;
    (&_S6558)->primal_0 = _S6459;
    (&_S6558)->differential_0 = 0.0f;
    DiffPair_float_0 _S6559;
    (&_S6559)->primal_0 = _S6435;
    (&_S6559)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6558, &_S6559, 0.0f);
    DiffPair_float_0 _S6560;
    (&_S6560)->primal_0 = _S6458;
    (&_S6560)->differential_0 = 0.0f;
    DiffPair_float_0 _S6561;
    (&_S6561)->primal_0 = _S6435;
    (&_S6561)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6560, &_S6561, 0.0f);
    DiffPair_float_0 _S6562;
    (&_S6562)->primal_0 = _S6457;
    (&_S6562)->differential_0 = 0.0f;
    DiffPair_float_0 _S6563;
    (&_S6563)->primal_0 = _S6434;
    (&_S6563)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6562, &_S6563, 0.0f);
    DiffPair_float_0 _S6564;
    (&_S6564)->primal_0 = _S6456;
    (&_S6564)->differential_0 = 0.0f;
    DiffPair_float_0 _S6565;
    (&_S6565)->primal_0 = _S6434;
    (&_S6565)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6564, &_S6565, 0.0f);
    DiffPair_float_0 _S6566;
    (&_S6566)->primal_0 = _S6455;
    (&_S6566)->differential_0 = 0.0f;
    DiffPair_float_0 _S6567;
    (&_S6567)->primal_0 = _S6410;
    (&_S6567)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6566, &_S6567, _S6558.differential_0);
    DiffPair_float_0 _S6568;
    (&_S6568)->primal_0 = _S6454;
    (&_S6568)->differential_0 = 0.0f;
    DiffPair_float_0 _S6569;
    (&_S6569)->primal_0 = _S6410;
    (&_S6569)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6568, &_S6569, _S6560.differential_0);
    DiffPair_float_0 _S6570;
    (&_S6570)->primal_0 = _S6453;
    (&_S6570)->differential_0 = 0.0f;
    DiffPair_float_0 _S6571;
    (&_S6571)->primal_0 = _S6409;
    (&_S6571)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6570, &_S6571, _S6562.differential_0);
    DiffPair_float_0 _S6572;
    (&_S6572)->primal_0 = _S6452;
    (&_S6572)->differential_0 = 0.0f;
    DiffPair_float_0 _S6573;
    (&_S6573)->primal_0 = _S6409;
    (&_S6573)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6572, &_S6573, _S6564.differential_0);
    DiffPair_float_0 _S6574;
    (&_S6574)->primal_0 = _S6451;
    (&_S6574)->differential_0 = 0.0f;
    DiffPair_float_0 _S6575;
    (&_S6575)->primal_0 = _S6385;
    (&_S6575)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6574, &_S6575, _S6566.differential_0);
    DiffPair_float_0 _S6576;
    (&_S6576)->primal_0 = _S6450;
    (&_S6576)->differential_0 = 0.0f;
    DiffPair_float_0 _S6577;
    (&_S6577)->primal_0 = _S6385;
    (&_S6577)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6576, &_S6577, _S6568.differential_0);
    DiffPair_float_0 _S6578;
    (&_S6578)->primal_0 = _S6449;
    (&_S6578)->differential_0 = 0.0f;
    DiffPair_float_0 _S6579;
    (&_S6579)->primal_0 = _S6384;
    (&_S6579)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6578, &_S6579, _S6570.differential_0);
    DiffPair_float_0 _S6580;
    (&_S6580)->primal_0 = _S6448;
    (&_S6580)->differential_0 = 0.0f;
    DiffPair_float_0 _S6581;
    (&_S6581)->primal_0 = _S6384;
    (&_S6581)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6580, &_S6581, _S6572.differential_0);
    DiffPair_float_0 _S6582;
    (&_S6582)->primal_0 = _S6447;
    (&_S6582)->differential_0 = 0.0f;
    DiffPair_float_0 _S6583;
    (&_S6583)->primal_0 = _S6360;
    (&_S6583)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6582, &_S6583, _S6574.differential_0);
    DiffPair_float_0 _S6584;
    (&_S6584)->primal_0 = _S6446;
    (&_S6584)->differential_0 = 0.0f;
    DiffPair_float_0 _S6585;
    (&_S6585)->primal_0 = _S6360;
    (&_S6585)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6584, &_S6585, _S6576.differential_0);
    DiffPair_float_0 _S6586;
    (&_S6586)->primal_0 = _S6445;
    (&_S6586)->differential_0 = 0.0f;
    DiffPair_float_0 _S6587;
    (&_S6587)->primal_0 = _S6359;
    (&_S6587)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6586, &_S6587, _S6578.differential_0);
    DiffPair_float_0 _S6588;
    (&_S6588)->primal_0 = _S6444;
    (&_S6588)->differential_0 = 0.0f;
    DiffPair_float_0 _S6589;
    (&_S6589)->primal_0 = _S6359;
    (&_S6589)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6588, &_S6589, _S6580.differential_0);
    DiffPair_float_0 _S6590;
    (&_S6590)->primal_0 = _S6443;
    (&_S6590)->differential_0 = 0.0f;
    DiffPair_float_0 _S6591;
    (&_S6591)->primal_0 = _S6335;
    (&_S6591)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6590, &_S6591, _S6582.differential_0);
    DiffPair_float_0 _S6592;
    (&_S6592)->primal_0 = _S6442;
    (&_S6592)->differential_0 = 0.0f;
    DiffPair_float_0 _S6593;
    (&_S6593)->primal_0 = _S6335;
    (&_S6593)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6592, &_S6593, _S6584.differential_0);
    DiffPair_float_0 _S6594;
    (&_S6594)->primal_0 = _S6441;
    (&_S6594)->differential_0 = 0.0f;
    DiffPair_float_0 _S6595;
    (&_S6595)->primal_0 = _S6334;
    (&_S6595)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6594, &_S6595, _S6586.differential_0);
    DiffPair_float_0 _S6596;
    (&_S6596)->primal_0 = _S6440;
    (&_S6596)->differential_0 = 0.0f;
    DiffPair_float_0 _S6597;
    (&_S6597)->primal_0 = _S6334;
    (&_S6597)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6596, &_S6597, _S6588.differential_0);
    DiffPair_float_0 _S6598;
    (&_S6598)->primal_0 = _S6439;
    (&_S6598)->differential_0 = 0.0f;
    DiffPair_float_0 _S6599;
    (&_S6599)->primal_0 = _S6310;
    (&_S6599)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6598, &_S6599, _S6590.differential_0);
    DiffPair_float_0 _S6600;
    (&_S6600)->primal_0 = _S6438;
    (&_S6600)->differential_0 = 0.0f;
    DiffPair_float_0 _S6601;
    (&_S6601)->primal_0 = _S6310;
    (&_S6601)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6600, &_S6601, _S6592.differential_0);
    DiffPair_float_0 _S6602;
    (&_S6602)->primal_0 = _S6437;
    (&_S6602)->differential_0 = 0.0f;
    DiffPair_float_0 _S6603;
    (&_S6603)->primal_0 = _S6309;
    (&_S6603)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6602, &_S6603, _S6594.differential_0);
    DiffPair_float_0 _S6604;
    (&_S6604)->primal_0 = _S6436;
    (&_S6604)->differential_0 = 0.0f;
    DiffPair_float_0 _S6605;
    (&_S6605)->primal_0 = _S6309;
    (&_S6605)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6604, &_S6605, _S6596.differential_0);
    DiffPair_float_0 _S6606;
    (&_S6606)->primal_0 = _S6260;
    (&_S6606)->differential_0 = 0.0f;
    DiffPair_float_0 _S6607;
    (&_S6607)->primal_0 = _S6285;
    (&_S6607)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6606, &_S6607, _S6598.differential_0);
    DiffPair_float_0 _S6608;
    (&_S6608)->primal_0 = _S6260;
    (&_S6608)->differential_0 = 0.0f;
    DiffPair_float_0 _S6609;
    (&_S6609)->primal_0 = _S6285;
    (&_S6609)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6608, &_S6609, _S6600.differential_0);
    DiffPair_float_0 _S6610;
    (&_S6610)->primal_0 = _S6259;
    (&_S6610)->differential_0 = 0.0f;
    DiffPair_float_0 _S6611;
    (&_S6611)->primal_0 = _S6284;
    (&_S6611)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6610, &_S6611, _S6602.differential_0);
    DiffPair_float_0 _S6612;
    (&_S6612)->primal_0 = _S6259;
    (&_S6612)->differential_0 = 0.0f;
    DiffPair_float_0 _S6613;
    (&_S6613)->primal_0 = _S6284;
    (&_S6613)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6612, &_S6613, _S6604.differential_0);
    float _S6614 = fx_41 * (_S6563.differential_0 + _S6565.differential_0);
    float2  _S6615 = make_float2 (_S6614, fy_41 * (_S6559.differential_0 + _S6561.differential_0)) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6614, (*dist_coeffs_45)[int(9)] * _S6614);
    float2  _S6616 = _S6423 * _S6615;
    float2  _S6617 = _S6427 * _S6615;
    float _S6618 = (*dist_coeffs_45)[int(4)] * _S6615.y;
    float _S6619 = (*dist_coeffs_45)[int(5)] * _S6615.x;
    float _S6620 = _S6616.x + _S6616.y;
    float _S6621 = r2_127 * _S6620;
    float _S6622 = r2_127 * _S6621;
    float _S6623 = (*dist_coeffs_45)[int(7)] * _S6615.y + _S6618 + (*dist_coeffs_45)[int(6)] * _S6615.x + _S6619 + _S6426 * _S6620 + _S6425 * _S6621 + _S6424 * _S6622 + (*dist_coeffs_45)[int(3)] * (r2_127 * _S6622);
    float _S6624 = v_132 * _S6623;
    float _S6625 = u_132 * _S6623;
    float _S6626 = _S6431 * _S6618 + 2.0f * (v_132 * _S6618) + _S6430 * _S6615.y + _S6428 * _S6615.x + _S6624 + _S6624;
    float _S6627 = _S6254 * (v_132 * _S6615.y) + _S6429 * _S6619 + 2.0f * (u_132 * _S6619) + _S6251 * (v_132 * _S6615.x) + _S6625 + _S6625;
    FixedArray<float3 , 16>  _S6628;
    _S6628[int(0)] = _S6488;
    _S6628[int(1)] = _S6488;
    _S6628[int(2)] = _S6488;
    _S6628[int(3)] = _S6488;
    _S6628[int(4)] = _S6488;
    _S6628[int(5)] = _S6488;
    _S6628[int(6)] = _S6488;
    _S6628[int(7)] = _S6488;
    _S6628[int(8)] = _S6488;
    _S6628[int(9)] = _S6488;
    _S6628[int(10)] = _S6488;
    _S6628[int(11)] = _S6488;
    _S6628[int(12)] = _S6488;
    _S6628[int(13)] = _S6488;
    _S6628[int(14)] = _S6488;
    _S6628[int(15)] = _S6488;
    _S6628[int(7)] = _S6513;
    _S6628[int(0)] = _S6547;
    _S6628[int(1)] = _S6534;
    _S6628[int(2)] = _S6532;
    _S6628[int(3)] = _S6530;
    _S6628[int(4)] = _S6519;
    _S6628[int(5)] = _S6517;
    _S6628[int(6)] = _S6515;
    _S6628[int(15)] = _S6491;
    _S6628[int(8)] = _S6511;
    _S6628[int(9)] = _S6503;
    _S6628[int(10)] = _S6501;
    _S6628[int(11)] = _S6499;
    _S6628[int(12)] = _S6497;
    _S6628[int(13)] = _S6495;
    _S6628[int(14)] = _S6493;
    float3  _S6629 = _S6628[int(0)];
    float3  _S6630 = _S6628[int(1)];
    float3  _S6631 = _S6628[int(2)];
    float3  _S6632 = _S6628[int(3)];
    float3  _S6633 = _S6628[int(4)];
    float3  _S6634 = _S6628[int(5)];
    float3  _S6635 = _S6628[int(6)];
    float3  _S6636 = _S6628[int(7)];
    float3  _S6637 = _S6628[int(8)];
    float3  _S6638 = _S6628[int(9)];
    float3  _S6639 = _S6628[int(10)];
    float3  _S6640 = _S6628[int(11)];
    float3  _S6641 = _S6628[int(12)];
    float3  _S6642 = _S6628[int(13)];
    float3  _S6643 = _S6628[int(14)];
    float3  _S6644 = _S6628[int(15)];
    float3  _S6645 = _S6548 + _S6557.differential_0 + _S6556.differential_0;
    float _S6646 = _S6611.differential_0 + _S6613.differential_0;
    float _S6647 = _S6606.differential_0 + _S6608.differential_0;
    float _S6648 = _S6607.differential_0 + _S6609.differential_0;
    float _S6649 = _S6603.differential_0 + _S6605.differential_0;
    float _S6650 = _S6610.differential_0 + _S6612.differential_0;
    float _S6651 = _S6567.differential_0 + _S6569.differential_0;
    float _S6652 = _S6571.differential_0 + _S6573.differential_0;
    float _S6653 = _S6575.differential_0 + _S6577.differential_0;
    float _S6654 = _S6579.differential_0 + _S6581.differential_0;
    float _S6655 = _S6583.differential_0 + _S6585.differential_0;
    float _S6656 = _S6587.differential_0 + _S6589.differential_0;
    float _S6657 = _S6591.differential_0 + _S6593.differential_0;
    float _S6658 = _S6595.differential_0 + _S6597.differential_0;
    float _S6659 = _S6599.differential_0 + _S6601.differential_0;
    float2  _S6660 = _S6617 + make_float2 (_S6627, _S6626);
    float2  _S6661 = _S6411 * _S6660;
    float2  _S6662 = _S6422 * _S6660;
    float _S6663 = _S6661.x + _S6661.y;
    if(_S6415)
    {
        float _S6664 = _S6663 / _S6417;
        float _S6665 = _S6418 * - _S6664;
        float _S6666 = _S6414 * (0.3333333432674408f * - (_S6413 * _S6664));
        k_15 = _S6666 + _S6666;
        _S6416 = _S6665;
        _S6417 = 0.0f;
    }
    else
    {
        float _S6667 = _S6663 / _S6416;
        float _S6668 = _S6414 * - _S6667;
        k_15 = _S6412 * _S6667;
        _S6416 = 0.0f;
        _S6417 = _S6668;
    }
    DiffPair_float_0 _S6669;
    (&_S6669)->primal_0 = _S6412;
    (&_S6669)->differential_0 = 0.0f;
    DiffPair_float_0 _S6670;
    (&_S6670)->primal_0 = _S6413;
    (&_S6670)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6669, &_S6670, k_15);
    float _S6671 = _S6670.differential_0 + _S6416;
    float _S6672 = _S6669.differential_0 + _S6417;
    float2  _S6673 = make_float2 (0.0f);
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6674;
    (&_S6674)->primal_0 = _S6411;
    (&_S6674)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6674, _S6672);
    float2  _S6675 = _S6674.differential_0 + _S6662;
    float _S6676 = fx_41 * _S6652;
    float2  _S6677 = make_float2 (_S6676, fy_41 * _S6651) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6676, (*dist_coeffs_45)[int(9)] * _S6676);
    float2  _S6678 = _S6398 * _S6677;
    float _S6679 = (*dist_coeffs_45)[int(4)] * _S6677.y;
    float _S6680 = (*dist_coeffs_45)[int(5)] * _S6677.x;
    float _S6681 = _S6678.x + _S6678.y;
    float _S6682 = r2_126 * _S6681;
    float _S6683 = r2_126 * _S6682;
    float _S6684 = (*dist_coeffs_45)[int(7)] * _S6677.y + _S6679 + (*dist_coeffs_45)[int(6)] * _S6677.x + _S6680 + _S6401 * _S6681 + _S6400 * _S6682 + _S6399 * _S6683 + (*dist_coeffs_45)[int(3)] * (r2_126 * _S6683);
    float _S6685 = v_131 * _S6684;
    float _S6686 = u_131 * _S6684;
    float3  _S6687 = make_float3 (_S6675.x, _S6675.y, _S6671);
    float2  _S6688 = _S6402 * _S6677 + make_float2 (_S6254 * (v_131 * _S6677.y) + _S6404 * _S6680 + 2.0f * (u_131 * _S6680) + _S6251 * (v_131 * _S6677.x) + _S6686 + _S6686, _S6406 * _S6679 + 2.0f * (v_131 * _S6679) + _S6405 * _S6677.y + _S6403 * _S6677.x + _S6685 + _S6685);
    float2  _S6689 = _S6386 * _S6688;
    float2  _S6690 = _S6397 * _S6688;
    float _S6691 = _S6689.x + _S6689.y;
    if(_S6390)
    {
        float _S6692 = _S6691 / _S6392;
        float _S6693 = _S6393 * - _S6692;
        float _S6694 = _S6389 * (0.3333333432674408f * - (_S6388 * _S6692));
        k_15 = _S6694 + _S6694;
        _S6391 = _S6693;
        _S6392 = 0.0f;
    }
    else
    {
        float _S6695 = _S6691 / _S6391;
        float _S6696 = _S6389 * - _S6695;
        k_15 = _S6387 * _S6695;
        _S6391 = 0.0f;
        _S6392 = _S6696;
    }
    DiffPair_float_0 _S6697;
    (&_S6697)->primal_0 = _S6387;
    (&_S6697)->differential_0 = 0.0f;
    DiffPair_float_0 _S6698;
    (&_S6698)->primal_0 = _S6388;
    (&_S6698)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6697, &_S6698, k_15);
    float _S6699 = _S6698.differential_0 + _S6391;
    float _S6700 = _S6697.differential_0 + _S6392;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6701;
    (&_S6701)->primal_0 = _S6386;
    (&_S6701)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6701, _S6700);
    float2  _S6702 = _S6701.differential_0 + _S6690;
    float _S6703 = fx_41 * _S6654;
    float2  _S6704 = make_float2 (_S6703, fy_41 * _S6653) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6703, (*dist_coeffs_45)[int(9)] * _S6703);
    float2  _S6705 = _S6373 * _S6704;
    float _S6706 = (*dist_coeffs_45)[int(4)] * _S6704.y;
    float _S6707 = (*dist_coeffs_45)[int(5)] * _S6704.x;
    float _S6708 = _S6705.x + _S6705.y;
    float _S6709 = r2_125 * _S6708;
    float _S6710 = r2_125 * _S6709;
    float _S6711 = (*dist_coeffs_45)[int(7)] * _S6704.y + _S6706 + (*dist_coeffs_45)[int(6)] * _S6704.x + _S6707 + _S6376 * _S6708 + _S6375 * _S6709 + _S6374 * _S6710 + (*dist_coeffs_45)[int(3)] * (r2_125 * _S6710);
    float _S6712 = v_130 * _S6711;
    float _S6713 = u_130 * _S6711;
    float3  _S6714 = make_float3 (_S6702.x, _S6702.y, _S6699);
    float2  _S6715 = _S6377 * _S6704 + make_float2 (_S6254 * (v_130 * _S6704.y) + _S6379 * _S6707 + 2.0f * (u_130 * _S6707) + _S6251 * (v_130 * _S6704.x) + _S6713 + _S6713, _S6381 * _S6706 + 2.0f * (v_130 * _S6706) + _S6380 * _S6704.y + _S6378 * _S6704.x + _S6712 + _S6712);
    float2  _S6716 = _S6361 * _S6715;
    float2  _S6717 = _S6372 * _S6715;
    float _S6718 = _S6716.x + _S6716.y;
    if(_S6365)
    {
        float _S6719 = _S6718 / _S6367;
        float _S6720 = _S6368 * - _S6719;
        float _S6721 = _S6364 * (0.3333333432674408f * - (_S6363 * _S6719));
        k_15 = _S6721 + _S6721;
        _S6366 = _S6720;
        _S6367 = 0.0f;
    }
    else
    {
        float _S6722 = _S6718 / _S6366;
        float _S6723 = _S6364 * - _S6722;
        k_15 = _S6362 * _S6722;
        _S6366 = 0.0f;
        _S6367 = _S6723;
    }
    DiffPair_float_0 _S6724;
    (&_S6724)->primal_0 = _S6362;
    (&_S6724)->differential_0 = 0.0f;
    DiffPair_float_0 _S6725;
    (&_S6725)->primal_0 = _S6363;
    (&_S6725)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6724, &_S6725, k_15);
    float _S6726 = _S6725.differential_0 + _S6366;
    float _S6727 = _S6724.differential_0 + _S6367;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6728;
    (&_S6728)->primal_0 = _S6361;
    (&_S6728)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6728, _S6727);
    float2  _S6729 = _S6728.differential_0 + _S6717;
    float _S6730 = fx_41 * _S6656;
    float2  _S6731 = make_float2 (_S6730, fy_41 * _S6655) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6730, (*dist_coeffs_45)[int(9)] * _S6730);
    float2  _S6732 = _S6348 * _S6731;
    float _S6733 = (*dist_coeffs_45)[int(4)] * _S6731.y;
    float _S6734 = (*dist_coeffs_45)[int(5)] * _S6731.x;
    float _S6735 = _S6732.x + _S6732.y;
    float _S6736 = r2_124 * _S6735;
    float _S6737 = r2_124 * _S6736;
    float _S6738 = (*dist_coeffs_45)[int(7)] * _S6731.y + _S6733 + (*dist_coeffs_45)[int(6)] * _S6731.x + _S6734 + _S6351 * _S6735 + _S6350 * _S6736 + _S6349 * _S6737 + (*dist_coeffs_45)[int(3)] * (r2_124 * _S6737);
    float _S6739 = v_129 * _S6738;
    float _S6740 = u_129 * _S6738;
    float3  _S6741 = make_float3 (_S6729.x, _S6729.y, _S6726);
    float2  _S6742 = _S6352 * _S6731 + make_float2 (_S6254 * (v_129 * _S6731.y) + _S6354 * _S6734 + 2.0f * (u_129 * _S6734) + _S6251 * (v_129 * _S6731.x) + _S6740 + _S6740, _S6356 * _S6733 + 2.0f * (v_129 * _S6733) + _S6355 * _S6731.y + _S6353 * _S6731.x + _S6739 + _S6739);
    float2  _S6743 = _S6336 * _S6742;
    float2  _S6744 = _S6347 * _S6742;
    float _S6745 = _S6743.x + _S6743.y;
    if(_S6340)
    {
        float _S6746 = _S6745 / _S6342;
        float _S6747 = _S6343 * - _S6746;
        float _S6748 = _S6339 * (0.3333333432674408f * - (_S6338 * _S6746));
        k_15 = _S6748 + _S6748;
        _S6341 = _S6747;
        _S6342 = 0.0f;
    }
    else
    {
        float _S6749 = _S6745 / _S6341;
        float _S6750 = _S6339 * - _S6749;
        k_15 = _S6337 * _S6749;
        _S6341 = 0.0f;
        _S6342 = _S6750;
    }
    DiffPair_float_0 _S6751;
    (&_S6751)->primal_0 = _S6337;
    (&_S6751)->differential_0 = 0.0f;
    DiffPair_float_0 _S6752;
    (&_S6752)->primal_0 = _S6338;
    (&_S6752)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6751, &_S6752, k_15);
    float _S6753 = _S6752.differential_0 + _S6341;
    float _S6754 = _S6751.differential_0 + _S6342;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6755;
    (&_S6755)->primal_0 = _S6336;
    (&_S6755)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6755, _S6754);
    float2  _S6756 = _S6755.differential_0 + _S6744;
    float _S6757 = fx_41 * _S6658;
    float2  _S6758 = make_float2 (_S6757, fy_41 * _S6657) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6757, (*dist_coeffs_45)[int(9)] * _S6757);
    float2  _S6759 = _S6323 * _S6758;
    float _S6760 = (*dist_coeffs_45)[int(4)] * _S6758.y;
    float _S6761 = (*dist_coeffs_45)[int(5)] * _S6758.x;
    float _S6762 = _S6759.x + _S6759.y;
    float _S6763 = r2_123 * _S6762;
    float _S6764 = r2_123 * _S6763;
    float _S6765 = (*dist_coeffs_45)[int(7)] * _S6758.y + _S6760 + (*dist_coeffs_45)[int(6)] * _S6758.x + _S6761 + _S6326 * _S6762 + _S6325 * _S6763 + _S6324 * _S6764 + (*dist_coeffs_45)[int(3)] * (r2_123 * _S6764);
    float _S6766 = v_128 * _S6765;
    float _S6767 = u_128 * _S6765;
    float3  _S6768 = make_float3 (_S6756.x, _S6756.y, _S6753);
    float2  _S6769 = _S6327 * _S6758 + make_float2 (_S6254 * (v_128 * _S6758.y) + _S6329 * _S6761 + 2.0f * (u_128 * _S6761) + _S6251 * (v_128 * _S6758.x) + _S6767 + _S6767, _S6331 * _S6760 + 2.0f * (v_128 * _S6760) + _S6330 * _S6758.y + _S6328 * _S6758.x + _S6766 + _S6766);
    float2  _S6770 = _S6311 * _S6769;
    float2  _S6771 = _S6322 * _S6769;
    float _S6772 = _S6770.x + _S6770.y;
    if(_S6315)
    {
        float _S6773 = _S6772 / _S6317;
        float _S6774 = _S6318 * - _S6773;
        float _S6775 = _S6314 * (0.3333333432674408f * - (_S6313 * _S6773));
        k_15 = _S6775 + _S6775;
        _S6316 = _S6774;
        _S6317 = 0.0f;
    }
    else
    {
        float _S6776 = _S6772 / _S6316;
        float _S6777 = _S6314 * - _S6776;
        k_15 = _S6312 * _S6776;
        _S6316 = 0.0f;
        _S6317 = _S6777;
    }
    DiffPair_float_0 _S6778;
    (&_S6778)->primal_0 = _S6312;
    (&_S6778)->differential_0 = 0.0f;
    DiffPair_float_0 _S6779;
    (&_S6779)->primal_0 = _S6313;
    (&_S6779)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6778, &_S6779, k_15);
    float _S6780 = _S6779.differential_0 + _S6316;
    float _S6781 = _S6778.differential_0 + _S6317;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6782;
    (&_S6782)->primal_0 = _S6311;
    (&_S6782)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6782, _S6781);
    float2  _S6783 = _S6782.differential_0 + _S6771;
    float _S6784 = fx_41 * _S6649;
    float2  _S6785 = make_float2 (_S6784, fy_41 * _S6659) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6784, (*dist_coeffs_45)[int(9)] * _S6784);
    float2  _S6786 = _S6298 * _S6785;
    float _S6787 = (*dist_coeffs_45)[int(4)] * _S6785.y;
    float _S6788 = (*dist_coeffs_45)[int(5)] * _S6785.x;
    float _S6789 = _S6786.x + _S6786.y;
    float _S6790 = r2_122 * _S6789;
    float _S6791 = r2_122 * _S6790;
    float _S6792 = (*dist_coeffs_45)[int(7)] * _S6785.y + _S6787 + (*dist_coeffs_45)[int(6)] * _S6785.x + _S6788 + _S6301 * _S6789 + _S6300 * _S6790 + _S6299 * _S6791 + (*dist_coeffs_45)[int(3)] * (r2_122 * _S6791);
    float _S6793 = v_127 * _S6792;
    float _S6794 = u_127 * _S6792;
    float3  _S6795 = make_float3 (_S6783.x, _S6783.y, _S6780);
    float2  _S6796 = _S6302 * _S6785 + make_float2 (_S6254 * (v_127 * _S6785.y) + _S6304 * _S6788 + 2.0f * (u_127 * _S6788) + _S6251 * (v_127 * _S6785.x) + _S6794 + _S6794, _S6306 * _S6787 + 2.0f * (v_127 * _S6787) + _S6305 * _S6785.y + _S6303 * _S6785.x + _S6793 + _S6793);
    float2  _S6797 = _S6286 * _S6796;
    float2  _S6798 = _S6297 * _S6796;
    float _S6799 = _S6797.x + _S6797.y;
    if(_S6290)
    {
        float _S6800 = _S6799 / _S6292;
        float _S6801 = _S6293 * - _S6800;
        float _S6802 = _S6289 * (0.3333333432674408f * - (_S6288 * _S6800));
        k_15 = _S6802 + _S6802;
        _S6291 = _S6801;
        _S6292 = 0.0f;
    }
    else
    {
        float _S6803 = _S6799 / _S6291;
        float _S6804 = _S6289 * - _S6803;
        k_15 = _S6287 * _S6803;
        _S6291 = 0.0f;
        _S6292 = _S6804;
    }
    DiffPair_float_0 _S6805;
    (&_S6805)->primal_0 = _S6287;
    (&_S6805)->differential_0 = 0.0f;
    DiffPair_float_0 _S6806;
    (&_S6806)->primal_0 = _S6288;
    (&_S6806)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6805, &_S6806, k_15);
    float _S6807 = _S6806.differential_0 + _S6291;
    float _S6808 = _S6805.differential_0 + _S6292;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6809;
    (&_S6809)->primal_0 = _S6286;
    (&_S6809)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6809, _S6808);
    float2  _S6810 = _S6809.differential_0 + _S6798;
    float _S6811 = fx_41 * _S6646;
    float2  _S6812 = make_float2 (_S6811, fy_41 * _S6648) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6811, (*dist_coeffs_45)[int(9)] * _S6811);
    float2  _S6813 = _S6273 * _S6812;
    float _S6814 = (*dist_coeffs_45)[int(4)] * _S6812.y;
    float _S6815 = (*dist_coeffs_45)[int(5)] * _S6812.x;
    float _S6816 = _S6813.x + _S6813.y;
    float _S6817 = r2_121 * _S6816;
    float _S6818 = r2_121 * _S6817;
    float _S6819 = (*dist_coeffs_45)[int(7)] * _S6812.y + _S6814 + (*dist_coeffs_45)[int(6)] * _S6812.x + _S6815 + _S6276 * _S6816 + _S6275 * _S6817 + _S6274 * _S6818 + (*dist_coeffs_45)[int(3)] * (r2_121 * _S6818);
    float _S6820 = v_126 * _S6819;
    float _S6821 = u_126 * _S6819;
    float3  _S6822 = make_float3 (_S6810.x, _S6810.y, _S6807);
    float2  _S6823 = _S6277 * _S6812 + make_float2 (_S6254 * (v_126 * _S6812.y) + _S6279 * _S6815 + 2.0f * (u_126 * _S6815) + _S6251 * (v_126 * _S6812.x) + _S6821 + _S6821, _S6281 * _S6814 + 2.0f * (v_126 * _S6814) + _S6280 * _S6812.y + _S6278 * _S6812.x + _S6820 + _S6820);
    float2  _S6824 = _S6261 * _S6823;
    float2  _S6825 = _S6272 * _S6823;
    float _S6826 = _S6824.x + _S6824.y;
    if(_S6265)
    {
        float _S6827 = _S6826 / _S6267;
        float _S6828 = _S6268 * - _S6827;
        float _S6829 = _S6264 * (0.3333333432674408f * - (_S6263 * _S6827));
        k_15 = _S6829 + _S6829;
        _S6266 = _S6828;
        _S6267 = 0.0f;
    }
    else
    {
        float _S6830 = _S6826 / _S6266;
        float _S6831 = _S6264 * - _S6830;
        k_15 = _S6262 * _S6830;
        _S6266 = 0.0f;
        _S6267 = _S6831;
    }
    DiffPair_float_0 _S6832;
    (&_S6832)->primal_0 = _S6262;
    (&_S6832)->differential_0 = 0.0f;
    DiffPair_float_0 _S6833;
    (&_S6833)->primal_0 = _S6263;
    (&_S6833)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6832, &_S6833, k_15);
    float _S6834 = _S6833.differential_0 + _S6266;
    float _S6835 = _S6832.differential_0 + _S6267;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6836;
    (&_S6836)->primal_0 = _S6261;
    (&_S6836)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6836, _S6835);
    float2  _S6837 = _S6836.differential_0 + _S6825;
    float _S6838 = fx_41 * _S6650;
    float2  _S6839 = make_float2 (_S6838, fy_41 * _S6647) + make_float2 ((*dist_coeffs_45)[int(8)] * _S6838, (*dist_coeffs_45)[int(9)] * _S6838);
    float2  _S6840 = _S6246 * _S6839;
    float _S6841 = (*dist_coeffs_45)[int(4)] * _S6839.y;
    float _S6842 = (*dist_coeffs_45)[int(5)] * _S6839.x;
    float _S6843 = _S6840.x + _S6840.y;
    float _S6844 = r2_120 * _S6843;
    float _S6845 = r2_120 * _S6844;
    float _S6846 = (*dist_coeffs_45)[int(7)] * _S6839.y + _S6841 + (*dist_coeffs_45)[int(6)] * _S6839.x + _S6842 + _S6249 * _S6843 + _S6248 * _S6844 + _S6247 * _S6845 + (*dist_coeffs_45)[int(3)] * (r2_120 * _S6845);
    float _S6847 = v_125 * _S6846;
    float _S6848 = u_125 * _S6846;
    float3  _S6849 = make_float3 (_S6837.x, _S6837.y, _S6834);
    float2  _S6850 = _S6250 * _S6839 + make_float2 (_S6254 * (v_125 * _S6839.y) + _S6253 * _S6842 + 2.0f * (u_125 * _S6842) + _S6251 * (v_125 * _S6839.x) + _S6848 + _S6848, _S6256 * _S6841 + 2.0f * (v_125 * _S6841) + _S6255 * _S6839.y + _S6252 * _S6839.x + _S6847 + _S6847);
    float2  _S6851 = _S6234 * _S6850;
    float2  _S6852 = _S6245 * _S6850;
    float _S6853 = _S6851.x + _S6851.y;
    if(_S6238)
    {
        float _S6854 = _S6853 / _S6240;
        float _S6855 = _S6241 * - _S6854;
        float _S6856 = _S6237 * (0.3333333432674408f * - (_S6236 * _S6854));
        k_15 = _S6856 + _S6856;
        _S6239 = _S6855;
        _S6240 = 0.0f;
    }
    else
    {
        float _S6857 = _S6853 / _S6239;
        float _S6858 = _S6237 * - _S6857;
        k_15 = _S6235 * _S6857;
        _S6239 = 0.0f;
        _S6240 = _S6858;
    }
    DiffPair_float_0 _S6859;
    (&_S6859)->primal_0 = _S6235;
    (&_S6859)->differential_0 = 0.0f;
    DiffPair_float_0 _S6860;
    (&_S6860)->primal_0 = _S6236;
    (&_S6860)->differential_0 = 0.0f;
    s_bwd_prop_atan2_0(&_S6859, &_S6860, k_15);
    float _S6861 = _S6860.differential_0 + _S6239;
    float _S6862 = _S6859.differential_0 + _S6240;
    DiffPair_vectorx3Cfloatx2C2x3E_0 _S6863;
    (&_S6863)->primal_0 = _S6234;
    (&_S6863)->differential_0 = _S6673;
    s_bwd_length_impl_0(&_S6863, _S6862);
    float2  _S6864 = _S6863.differential_0 + _S6852;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6865;
    (&_S6865)->primal_0 = R_38;
    (&_S6865)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6866;
    (&_S6866)->primal_0 = _S6233;
    (&_S6866)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6865, &_S6866, _S6645);
    DiffPair_float_0 _S6867;
    (&_S6867)->primal_0 = _S6230;
    (&_S6867)->differential_0 = 0.0f;
    DiffPair_float_0 _S6868;
    (&_S6868)->primal_0 = _S6232;
    (&_S6868)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6867, &_S6868, 0.0f);
    DiffPair_float_0 _S6869;
    (&_S6869)->primal_0 = _S6229;
    (&_S6869)->differential_0 = 0.0f;
    DiffPair_float_0 _S6870;
    (&_S6870)->primal_0 = _S6232;
    (&_S6870)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6869, &_S6870, 0.0f);
    float _S6871 = _S6868.differential_0 + _S6870.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6872;
    (&_S6872)->primal_0 = _S6231;
    (&_S6872)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6872, _S6871);
    float3  _S6873 = _S6872.differential_0 + _S6687;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6874;
    (&_S6874)->primal_0 = R_38;
    (&_S6874)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6875;
    (&_S6875)->primal_0 = pos_i_13;
    (&_S6875)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6874, &_S6875, _S6873);
    DiffPair_float_0 _S6876;
    (&_S6876)->primal_0 = _S6226;
    (&_S6876)->differential_0 = 0.0f;
    DiffPair_float_0 _S6877;
    (&_S6877)->primal_0 = _S6228;
    (&_S6877)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6876, &_S6877, _S6867.differential_0);
    DiffPair_float_0 _S6878;
    (&_S6878)->primal_0 = _S6225;
    (&_S6878)->differential_0 = 0.0f;
    DiffPair_float_0 _S6879;
    (&_S6879)->primal_0 = _S6228;
    (&_S6879)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6878, &_S6879, _S6869.differential_0);
    float _S6880 = _S6877.differential_0 + _S6879.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6881;
    (&_S6881)->primal_0 = _S6227;
    (&_S6881)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6881, _S6880);
    float3  _S6882 = _S6881.differential_0 + _S6714;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6883;
    (&_S6883)->primal_0 = R_38;
    (&_S6883)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6884;
    (&_S6884)->primal_0 = pos_i_12;
    (&_S6884)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6883, &_S6884, _S6882);
    DiffPair_float_0 _S6885;
    (&_S6885)->primal_0 = _S6222;
    (&_S6885)->differential_0 = 0.0f;
    DiffPair_float_0 _S6886;
    (&_S6886)->primal_0 = _S6224;
    (&_S6886)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6885, &_S6886, _S6876.differential_0);
    DiffPair_float_0 _S6887;
    (&_S6887)->primal_0 = _S6221;
    (&_S6887)->differential_0 = 0.0f;
    DiffPair_float_0 _S6888;
    (&_S6888)->primal_0 = _S6224;
    (&_S6888)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6887, &_S6888, _S6878.differential_0);
    float _S6889 = _S6886.differential_0 + _S6888.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6890;
    (&_S6890)->primal_0 = _S6223;
    (&_S6890)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6890, _S6889);
    float3  _S6891 = _S6890.differential_0 + _S6741;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6892;
    (&_S6892)->primal_0 = R_38;
    (&_S6892)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6893;
    (&_S6893)->primal_0 = pos_i_11;
    (&_S6893)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6892, &_S6893, _S6891);
    DiffPair_float_0 _S6894;
    (&_S6894)->primal_0 = _S6218;
    (&_S6894)->differential_0 = 0.0f;
    DiffPair_float_0 _S6895;
    (&_S6895)->primal_0 = _S6220;
    (&_S6895)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6894, &_S6895, _S6885.differential_0);
    DiffPair_float_0 _S6896;
    (&_S6896)->primal_0 = _S6217;
    (&_S6896)->differential_0 = 0.0f;
    DiffPair_float_0 _S6897;
    (&_S6897)->primal_0 = _S6220;
    (&_S6897)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6896, &_S6897, _S6887.differential_0);
    float _S6898 = _S6895.differential_0 + _S6897.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6899;
    (&_S6899)->primal_0 = _S6219;
    (&_S6899)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6899, _S6898);
    float3  _S6900 = _S6899.differential_0 + _S6768;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6901;
    (&_S6901)->primal_0 = R_38;
    (&_S6901)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6902;
    (&_S6902)->primal_0 = pos_i_10;
    (&_S6902)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6901, &_S6902, _S6900);
    DiffPair_float_0 _S6903;
    (&_S6903)->primal_0 = _S6214;
    (&_S6903)->differential_0 = 0.0f;
    DiffPair_float_0 _S6904;
    (&_S6904)->primal_0 = _S6216;
    (&_S6904)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6903, &_S6904, _S6894.differential_0);
    DiffPair_float_0 _S6905;
    (&_S6905)->primal_0 = _S6213;
    (&_S6905)->differential_0 = 0.0f;
    DiffPair_float_0 _S6906;
    (&_S6906)->primal_0 = _S6216;
    (&_S6906)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6905, &_S6906, _S6896.differential_0);
    float _S6907 = _S6904.differential_0 + _S6906.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6908;
    (&_S6908)->primal_0 = _S6215;
    (&_S6908)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6908, _S6907);
    float3  _S6909 = _S6908.differential_0 + _S6795;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6910;
    (&_S6910)->primal_0 = R_38;
    (&_S6910)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6911;
    (&_S6911)->primal_0 = pos_i_9;
    (&_S6911)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6910, &_S6911, _S6909);
    DiffPair_float_0 _S6912;
    (&_S6912)->primal_0 = _S6210;
    (&_S6912)->differential_0 = 0.0f;
    DiffPair_float_0 _S6913;
    (&_S6913)->primal_0 = _S6212;
    (&_S6913)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6912, &_S6913, _S6903.differential_0);
    DiffPair_float_0 _S6914;
    (&_S6914)->primal_0 = _S6209;
    (&_S6914)->differential_0 = 0.0f;
    DiffPair_float_0 _S6915;
    (&_S6915)->primal_0 = _S6212;
    (&_S6915)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6914, &_S6915, _S6905.differential_0);
    float _S6916 = _S6913.differential_0 + _S6915.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6917;
    (&_S6917)->primal_0 = _S6211;
    (&_S6917)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6917, _S6916);
    float3  _S6918 = _S6917.differential_0 + _S6822;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6919;
    (&_S6919)->primal_0 = R_38;
    (&_S6919)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6920;
    (&_S6920)->primal_0 = pos_i_8;
    (&_S6920)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6919, &_S6920, _S6918);
    DiffPair_float_0 _S6921;
    (&_S6921)->primal_0 = _S6206;
    (&_S6921)->differential_0 = 0.0f;
    DiffPair_float_0 _S6922;
    (&_S6922)->primal_0 = _S6208;
    (&_S6922)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6921, &_S6922, _S6912.differential_0);
    DiffPair_float_0 _S6923;
    (&_S6923)->primal_0 = _S6205;
    (&_S6923)->differential_0 = 0.0f;
    DiffPair_float_0 _S6924;
    (&_S6924)->primal_0 = _S6208;
    (&_S6924)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6923, &_S6924, _S6914.differential_0);
    float _S6925 = _S6922.differential_0 + _S6924.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6926;
    (&_S6926)->primal_0 = _S6207;
    (&_S6926)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6926, _S6925);
    float3  _S6927 = _S6926.differential_0 + _S6849;
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6928;
    (&_S6928)->primal_0 = R_38;
    (&_S6928)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6929;
    (&_S6929)->primal_0 = pos_i_7;
    (&_S6929)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6928, &_S6929, _S6927);
    DiffPair_float_0 _S6930;
    (&_S6930)->primal_0 = 0.0f;
    (&_S6930)->differential_0 = 0.0f;
    DiffPair_float_0 _S6931;
    (&_S6931)->primal_0 = _S6204;
    (&_S6931)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S6930, &_S6931, _S6921.differential_0);
    DiffPair_float_0 _S6932;
    (&_S6932)->primal_0 = 1.00000001504746622e+30f;
    (&_S6932)->differential_0 = 0.0f;
    DiffPair_float_0 _S6933;
    (&_S6933)->primal_0 = _S6204;
    (&_S6933)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S6932, &_S6933, _S6923.differential_0);
    float _S6934 = _S6931.differential_0 + _S6933.differential_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6935;
    (&_S6935)->primal_0 = _S6203;
    (&_S6935)->differential_0 = _S6488;
    s_bwd_length_impl_1(&_S6935, _S6934);
    float3  _S6936 = _S6935.differential_0 + make_float3 (_S6864.x, _S6864.y, _S6861);
    DiffPair_matrixx3Cfloatx2C3x2C3x3E_0 _S6937;
    (&_S6937)->primal_0 = R_38;
    (&_S6937)->differential_0 = _S6550;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S6938;
    (&_S6938)->primal_0 = pos_5;
    (&_S6938)->differential_0 = _S6488;
    s_bwd_prop_mul_1(&_S6937, &_S6938, _S6936);
    float3  _S6939 = _S6645 + _S6873 + _S6882 + _S6891 + _S6900 + _S6909 + _S6918 + _S6927 + _S6936 + _S6553.differential_0;
    Matrix<float, 3, 3>  _S6940 = _S6865.differential_0 + _S6874.differential_0 + _S6883.differential_0 + _S6892.differential_0 + _S6901.differential_0 + _S6910.differential_0 + _S6919.differential_0 + _S6928.differential_0 + _S6937.differential_0 + _S6554;
    (*v_densities_1)[int(0)] = 0.0f;
    (*v_densities_1)[int(1)] = 0.0f;
    (*v_densities_1)[int(2)] = 0.0f;
    (*v_densities_1)[int(3)] = 0.0f;
    (*v_densities_1)[int(4)] = 0.0f;
    (*v_densities_1)[int(5)] = 0.0f;
    (*v_densities_1)[int(6)] = 0.0f;
    (*v_densities_1)[int(7)] = 0.0f;
    (*v_sh_coeffs_10)[int(0)] = _S6629;
    (*v_sh_coeffs_10)[int(1)] = _S6630;
    (*v_sh_coeffs_10)[int(2)] = _S6631;
    (*v_sh_coeffs_10)[int(3)] = _S6632;
    (*v_sh_coeffs_10)[int(4)] = _S6633;
    (*v_sh_coeffs_10)[int(5)] = _S6634;
    (*v_sh_coeffs_10)[int(6)] = _S6635;
    (*v_sh_coeffs_10)[int(7)] = _S6636;
    (*v_sh_coeffs_10)[int(8)] = _S6637;
    (*v_sh_coeffs_10)[int(9)] = _S6638;
    (*v_sh_coeffs_10)[int(10)] = _S6639;
    (*v_sh_coeffs_10)[int(11)] = _S6640;
    (*v_sh_coeffs_10)[int(12)] = _S6641;
    (*v_sh_coeffs_10)[int(13)] = _S6642;
    (*v_sh_coeffs_10)[int(14)] = _S6643;
    (*v_sh_coeffs_10)[int(15)] = _S6644;
    *v_R_12 = _S6940;
    *v_t_11 = _S6939;
    return;
}

inline __device__ void _d_abs_vector_0(DiffPair_vectorx3Cfloatx2C3x3E_0 * dpx_19, float3  dOut_21)
{
    float3  _S6941 = _slang_select(((*dpx_19).primal_0) > make_float3 (0.0f), make_float3 (1.0f),_slang_select(((*dpx_19).primal_0) == make_float3 (0.0f), make_float3 (0.0f),make_float3 (-1.0f))) * dOut_21;
    dpx_19->primal_0 = (*dpx_19).primal_0;
    dpx_19->differential_0 = _S6941;
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

inline __device__ bool ray_aabb_intersection(float3  ray_o_10, float3  ray_d_10, float3  center_0, float radius_0, float * t0_0, float * t1_0)
{
    float3  m_2 = make_float3 (1.0f) / ray_d_10;
    float3  k_16 = abs_0(m_2) * make_float3 (radius_0);
    float3  _S6942 = - (m_2 * (ray_o_10 - center_0));
    float3  ta_0 = _S6942 - k_16;
    float3  tb_0 = _S6942 + k_16;
    *t0_0 = (F32_max(((F32_max((ta_0.x), (ta_0.y)))), ((F32_max((ta_0.z), (0.0f))))));
    float _S6943 = (F32_min(((F32_min((tb_0.x), (tb_0.y)))), (tb_0.z)));
    *t1_0 = _S6943;
    return (*t0_0) < _S6943;
}

inline __device__ void _d_lerp_0(DiffPair_float_0 * dpx_20, DiffPair_float_0 * dpy_7, DiffPair_float_0 * dps_0, float dOut_22)
{
    float _S6944 = (1.0f - (*dps_0).primal_0) * dOut_22;
    dpx_20->primal_0 = (*dpx_20).primal_0;
    dpx_20->differential_0 = _S6944;
    DiffPair_float_0 _S6945 = *dpy_7;
    float _S6946 = (*dps_0).primal_0 * dOut_22;
    dpy_7->primal_0 = (*dpy_7).primal_0;
    dpy_7->differential_0 = _S6946;
    float _S6947 = (_S6945.primal_0 - (*dpx_20).primal_0) * dOut_22;
    dps_0->primal_0 = _S6945.primal_0;
    dps_0->differential_0 = _S6947;
    return;
}

inline __device__ float lerp_0(float x_70, float y_36, float s_0)
{
    return x_70 + (y_36 - x_70) * s_0;
}

inline __device__ float interp_0(FixedArray<float, 8>  * densities_6, float3  w_2)
{
    float _S6948 = w_2.z;
    float _S6949 = 1.0f - _S6948;
    float _S6950 = w_2.y;
    float _S6951 = 1.0f - _S6950;
    float _S6952 = _S6949 * _S6951;
    float _S6953 = w_2.x;
    float _S6954 = 1.0f - _S6953;
    float _S6955 = _S6949 * _S6950;
    float _S6956 = _S6948 * _S6951;
    float _S6957 = _S6948 * _S6950;
    return _S6952 * _S6954 * (*densities_6)[int(0)] + _S6952 * _S6953 * (*densities_6)[int(1)] + _S6955 * _S6954 * (*densities_6)[int(2)] + _S6955 * _S6953 * (*densities_6)[int(3)] + _S6956 * _S6954 * (*densities_6)[int(4)] + _S6956 * _S6953 * (*densities_6)[int(5)] + _S6957 * _S6954 * (*densities_6)[int(6)] + _S6957 * _S6953 * (*densities_6)[int(7)];
}

inline __device__ float evaluate_alpha_voxel(float3  pos_6, float size_6, FixedArray<float, 8>  * densities_7, float3  ray_o_11, float3  ray_d_11)
{
    float _S6958 = 0.5f * size_6;
    float3  m_3 = make_float3 (1.0f) / ray_d_11;
    float3  k_17 = abs_0(m_3) * make_float3 (_S6958);
    float3  _S6959 = - (m_3 * (ray_o_11 - (pos_6 + make_float3 (_S6958))));
    float3  ta_1 = _S6959 - k_17;
    float3  tb_1 = _S6959 + k_17;
    float t0_1 = (F32_max(((F32_max((ta_1.x), (ta_1.y)))), ((F32_max((ta_1.z), (0.0f))))));
    float t1_1 = (F32_min(((F32_min((tb_1.x), (tb_1.y)))), (tb_1.z)));
    if(!(t0_1 < t1_1))
    {
        return 0.0f;
    }
    int i_12 = int(0);
    float accum_0 = 0.0f;
    for(;;)
    {
        if(i_12 < int(8))
        {
        }
        else
        {
            break;
        }
        float _S6960 = interp_0(densities_7, (ray_o_11 + ray_d_11 * make_float3 (lerp_0(t0_1, t1_1, (float(i_12) + 0.5f) / 8.0f)) - pos_6) / make_float3 (size_6));
        float _S6961;
        if(_S6960 > 1.10000002384185791f)
        {
            _S6961 = _S6960;
        }
        else
        {
            _S6961 = (F32_exp((0.90909093618392944f * _S6960 - 0.90468984842300415f)));
        }
        float accum_1 = accum_0 + _S6961;
        i_12 = i_12 + int(1);
        accum_0 = accum_1;
    }
    return (F32_min((1.0f - (F32_exp((- (t1_1 - t0_1) / 8.0f * accum_0)))), (0.99900001287460327f)));
}

struct DiffPair_arrayx3Cfloatx2C8x3E_0
{
    FixedArray<float, 8>  primal_0;
    FixedArray<float, 8>  differential_0;
};

struct s_bwd_prop_evaluate_alpha_voxel_Intermediates_0
{
    float _S6962;
};

inline __device__ float3  s_primal_ctx_abs_1(float3  _S6963)
{
    return abs_0(_S6963);
}

inline __device__ float s_primal_ctx_lerp_0(float _S6964, float _S6965, float _S6966)
{
    return lerp_0(_S6964, _S6965, _S6966);
}

inline __device__ float s_primal_ctx_interp_0(FixedArray<float, 8>  * dpdensities_0, float3  dpw_0)
{
    float _S6967 = dpw_0.z;
    float _S6968 = 1.0f - _S6967;
    float _S6969 = dpw_0.y;
    float _S6970 = 1.0f - _S6969;
    float _S6971 = _S6968 * _S6970;
    float _S6972 = dpw_0.x;
    float _S6973 = 1.0f - _S6972;
    float _S6974 = _S6968 * _S6969;
    float _S6975 = _S6967 * _S6970;
    float _S6976 = _S6967 * _S6969;
    return _S6971 * _S6973 * (*dpdensities_0)[int(0)] + _S6971 * _S6972 * (*dpdensities_0)[int(1)] + _S6974 * _S6973 * (*dpdensities_0)[int(2)] + _S6974 * _S6972 * (*dpdensities_0)[int(3)] + _S6975 * _S6973 * (*dpdensities_0)[int(4)] + _S6975 * _S6972 * (*dpdensities_0)[int(5)] + _S6976 * _S6973 * (*dpdensities_0)[int(6)] + _S6976 * _S6972 * (*dpdensities_0)[int(7)];
}

inline __device__ float s_primal_ctx_evaluate_alpha_voxel_0(float3  pos_7, float size_7, FixedArray<float, 8>  * dpdensities_1, float3  dpray_o_4, float3  dpray_d_4, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_14)
{
    _s_diff_ctx_14->_S6962 = 0.0f;
    _s_diff_ctx_14->_S6962 = 0.0f;
    float _S6977 = 0.5f * size_7;
    float3  m_4 = make_float3 (1.0f) / dpray_d_4;
    float3  k_18 = s_primal_ctx_abs_1(m_4) * make_float3 (_S6977);
    float3  _S6978 = - (m_4 * (dpray_o_4 - (pos_7 + make_float3 (_S6977))));
    float3  ta_2 = _S6978 - k_18;
    float3  tb_2 = _S6978 + k_18;
    float _S6979 = s_primal_ctx_max_0(s_primal_ctx_max_0(ta_2.x, ta_2.y), s_primal_ctx_max_0(ta_2.z, 0.0f));
    float _S6980 = s_primal_ctx_min_0(s_primal_ctx_min_0(tb_2.x, tb_2.y), tb_2.z);
    float accum_2;
    if(!!(_S6979 < _S6980))
    {
        float _S6981 = - (_S6980 - _S6979);
        bool _runFlag_0 = true;
        accum_2 = 0.0f;
        int i_13 = int(0);
        int _pc_0 = int(0);
        for(;;)
        {
            _s_diff_ctx_14->_S6962 = accum_2;
            if(_runFlag_0)
            {
            }
            else
            {
                break;
            }
            float _S6982;
            int _S6983;
            if(i_13 < int(8))
            {
                float _S6984 = s_primal_ctx_interp_0(dpdensities_1, (dpray_o_4 + dpray_d_4 * make_float3 (s_primal_ctx_lerp_0(_S6979, _S6980, (float(i_13) + 0.5f) / 8.0f)) - pos_7) / make_float3 (size_7));
                if(_S6984 > 1.10000002384185791f)
                {
                    _S6982 = _S6984;
                }
                else
                {
                    _S6982 = s_primal_ctx_exp_1(0.90909093618392944f * _S6984 - 0.90468984842300415f);
                }
                float accum_3 = accum_2 + _S6982;
                _S6983 = int(2);
                _S6982 = accum_3;
            }
            else
            {
                _S6983 = int(1);
                _S6982 = 0.0f;
            }
            if(_S6983 != int(2))
            {
                _runFlag_0 = false;
            }
            if(_runFlag_0)
            {
                int _S6985 = i_13 + int(1);
                accum_2 = _S6982;
                i_13 = _S6985;
            }
            _pc_0 = _pc_0 + int(1);
        }
        accum_2 = s_primal_ctx_min_0(1.0f - s_primal_ctx_exp_1(_S6981 / 8.0f * accum_2), 0.99900001287460327f);
    }
    else
    {
        accum_2 = 0.0f;
    }
    return accum_2;
}

inline __device__ void s_bwd_prop_abs_1(DiffPair_vectorx3Cfloatx2C3x3E_0 * _S6986, float3  _S6987)
{
    _d_abs_vector_0(_S6986, _S6987);
    return;
}

inline __device__ void s_bwd_prop_interp_0(DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpw_1, float _s_dOut_14)
{
    float _S6988 = (*dpw_1).primal_0.z;
    float _S6989 = 1.0f - _S6988;
    float _S6990 = (*dpw_1).primal_0.y;
    float _S6991 = 1.0f - _S6990;
    float _S6992 = _S6989 * _S6991;
    float _S6993 = (*dpw_1).primal_0.x;
    float _S6994 = 1.0f - _S6993;
    float _S6995 = _S6989 * _S6990;
    float _S6996 = _S6988 * _S6991;
    float _S6997 = _S6988 * _S6990;
    float _S6998 = _S6997 * _S6993 * _s_dOut_14;
    float s_diff_w7_T_0 = dpdensities_2->primal_0[int(7)] * _s_dOut_14;
    float _S6999 = _S6997 * _S6994 * _s_dOut_14;
    float s_diff_w6_T_0 = dpdensities_2->primal_0[int(6)] * _s_dOut_14;
    float _S7000 = _S6996 * _S6993 * _s_dOut_14;
    float s_diff_w5_T_0 = dpdensities_2->primal_0[int(5)] * _s_dOut_14;
    float _S7001 = _S6996 * _S6994 * _s_dOut_14;
    float s_diff_w4_T_0 = dpdensities_2->primal_0[int(4)] * _s_dOut_14;
    float _S7002 = _S6995 * _S6993 * _s_dOut_14;
    float s_diff_w3_T_0 = dpdensities_2->primal_0[int(3)] * _s_dOut_14;
    float _S7003 = _S6995 * _S6994 * _s_dOut_14;
    float s_diff_w2_T_0 = dpdensities_2->primal_0[int(2)] * _s_dOut_14;
    float _S7004 = _S6992 * _S6993 * _s_dOut_14;
    float s_diff_w1_T_0 = dpdensities_2->primal_0[int(1)] * _s_dOut_14;
    float _S7005 = _S6992 * _S6994 * _s_dOut_14;
    float s_diff_w0_T_0 = dpdensities_2->primal_0[int(0)] * _s_dOut_14;
    float _S7006 = _S6993 * s_diff_w7_T_0 + _S6994 * s_diff_w6_T_0;
    float _S7007 = _S6993 * s_diff_w5_T_0 + _S6994 * s_diff_w4_T_0;
    float _S7008 = _S6993 * s_diff_w3_T_0 + _S6994 * s_diff_w2_T_0;
    float _S7009 = _S6993 * s_diff_w1_T_0 + _S6994 * s_diff_w0_T_0;
    float3  _S7010 = make_float3 (_S6997 * s_diff_w7_T_0 + _S6996 * s_diff_w5_T_0 + _S6995 * s_diff_w3_T_0 + _S6992 * s_diff_w1_T_0 + - (_S6997 * s_diff_w6_T_0 + _S6996 * s_diff_w4_T_0 + _S6995 * s_diff_w2_T_0 + _S6992 * s_diff_w0_T_0), _S6988 * _S7006 + _S6989 * _S7008 + - (_S6988 * _S7007 + _S6989 * _S7009), _S6990 * _S7006 + _S6991 * _S7007 + - (_S6990 * _S7008 + _S6991 * _S7009));
    dpw_1->primal_0 = (*dpw_1).primal_0;
    dpw_1->differential_0 = _S7010;
    FixedArray<float, 8>  _S7011;
    _S7011[int(0)] = 0.0f;
    _S7011[int(1)] = 0.0f;
    _S7011[int(2)] = 0.0f;
    _S7011[int(3)] = 0.0f;
    _S7011[int(4)] = 0.0f;
    _S7011[int(5)] = 0.0f;
    _S7011[int(6)] = 0.0f;
    _S7011[int(7)] = 0.0f;
    _S7011[int(7)] = _S6998;
    _S7011[int(6)] = _S6999;
    _S7011[int(5)] = _S7000;
    _S7011[int(4)] = _S7001;
    _S7011[int(3)] = _S7002;
    _S7011[int(2)] = _S7003;
    _S7011[int(1)] = _S7004;
    _S7011[int(0)] = _S7005;
    dpdensities_2->primal_0 = dpdensities_2->primal_0;
    dpdensities_2->differential_0 = _S7011;
    return;
}

inline __device__ void s_bwd_prop_lerp_0(DiffPair_float_0 * _S7012, DiffPair_float_0 * _S7013, DiffPair_float_0 * _S7014, float _S7015)
{
    _d_lerp_0(_S7012, _S7013, _S7014, _S7015);
    return;
}

inline __device__ void s_bwd_prop_evaluate_alpha_voxel_0(float3  pos_8, float size_8, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_3, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_5, float _s_dOut_15, s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 * _s_diff_ctx_15)
{
    FixedArray<float, 8>  _S7016 = dpdensities_3->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7017 = *dpray_o_5;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7018 = *dpray_d_5;
    float _S7019 = 0.5f * size_8;
    float3  _S7020 = make_float3 (_S7019);
    float3  _S7021 = make_float3 (size_8);
    float3  m_5 = make_float3 (1.0f) / (*dpray_d_5).primal_0;
    float3  _S7022 = (*dpray_d_5).primal_0 * (*dpray_d_5).primal_0;
    float3  _S7023 = (*dpray_o_5).primal_0 - (pos_8 + make_float3 (_S7019));
    float3  k_19 = s_primal_ctx_abs_1(m_5) * make_float3 (_S7019);
    float3  _S7024 = - (m_5 * _S7023);
    float3  ta_3 = _S7024 - k_19;
    float3  tb_3 = _S7024 + k_19;
    float _S7025 = ta_3.x;
    float _S7026 = ta_3.y;
    float _S7027 = s_primal_ctx_max_0(_S7025, _S7026);
    float _S7028 = ta_3.z;
    float _S7029 = s_primal_ctx_max_0(_S7028, 0.0f);
    float _S7030 = s_primal_ctx_max_0(_S7027, _S7029);
    float _S7031 = tb_3.x;
    float _S7032 = tb_3.y;
    float _S7033 = s_primal_ctx_min_0(_S7031, _S7032);
    float _S7034 = tb_3.z;
    float _S7035 = s_primal_ctx_min_0(_S7033, _S7034);
    bool _S7036 = !!(_S7030 < _S7035);
    float _S7037;
    float _S7038;
    float _S7039;
    if(_S7036)
    {
        float _S7040 = - (_S7035 - _S7030) / 8.0f;
        float _S7041 = _S7040 * _s_diff_ctx_15->_S6962;
        _S7037 = 1.0f - s_primal_ctx_exp_1(_S7041);
        _S7038 = _S7041;
        _S7039 = _S7040;
    }
    else
    {
        _S7037 = 0.0f;
        _S7038 = 0.0f;
        _S7039 = 0.0f;
    }
    float3  _S7042 = make_float3 (0.0f);
    FixedArray<float, 8>  _S7043 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    float3  _S7044;
    float3  _S7045;
    FixedArray<float, 8>  _S7046;
    if(_S7036)
    {
        DiffPair_float_0 _S7047;
        (&_S7047)->primal_0 = _S7037;
        (&_S7047)->differential_0 = 0.0f;
        DiffPair_float_0 _S7048;
        (&_S7048)->primal_0 = 0.99900001287460327f;
        (&_S7048)->differential_0 = 0.0f;
        s_bwd_prop_min_0(&_S7047, &_S7048, _s_dOut_15);
        float _S7049 = - _S7047.differential_0;
        DiffPair_float_0 _S7050;
        (&_S7050)->primal_0 = _S7038;
        (&_S7050)->differential_0 = 0.0f;
        s_bwd_prop_exp_0(&_S7050, _S7049);
        float _S7051 = _S7039 * _S7050.differential_0;
        float _S7052 = 0.125f * (_s_diff_ctx_15->_S6962 * _S7050.differential_0);
        int _dc_0 = int(8);
        _S7037 = _S7051;
        _S7038 = 0.0f;
        _S7039 = 0.0f;
        _S7044 = _S7042;
        _S7045 = _S7042;
        _S7046[int(0)] = 0.0f;
        _S7046[int(1)] = 0.0f;
        _S7046[int(2)] = 0.0f;
        _S7046[int(3)] = 0.0f;
        _S7046[int(4)] = 0.0f;
        _S7046[int(5)] = 0.0f;
        _S7046[int(6)] = 0.0f;
        _S7046[int(7)] = 0.0f;
        for(;;)
        {
            if(_dc_0 >= int(0))
            {
            }
            else
            {
                break;
            }
            bool _S7053 = _dc_0 < int(8);
            float _S7054;
            float _S7055;
            int _S7056;
            float3  _S7057;
            float3  _S7058;
            bool _S7059;
            if(_S7053)
            {
                float _S7060 = (float(_dc_0) + 0.5f) / 8.0f;
                float _S7061 = s_primal_ctx_lerp_0(_S7030, _S7035, _S7060);
                float3  _S7062 = make_float3 (_S7061);
                float3  _S7063 = (_S7017.primal_0 + _S7018.primal_0 * make_float3 (_S7061) - pos_8) / make_float3 (size_8);
                FixedArray<float, 8>  _S7064 = _S7016;
                float _S7065 = s_primal_ctx_interp_0(&_S7064, _S7063);
                bool _S7066 = _S7065 > 1.10000002384185791f;
                if(_S7066)
                {
                    _S7054 = 0.0f;
                }
                else
                {
                    _S7054 = 0.90909093618392944f * _S7065 - 0.90468984842300415f;
                }
                _S7056 = int(2);
                _S7059 = _S7066;
                _S7057 = _S7063;
                _S7058 = _S7062;
                _S7055 = _S7060;
            }
            else
            {
                _S7056 = int(1);
                _S7059 = false;
                _S7054 = 0.0f;
                _S7057 = _S7042;
                _S7058 = _S7042;
                _S7055 = 0.0f;
            }
            float _S7067;
            float _S7068;
            if(!(_S7056 != int(2)))
            {
                _S7067 = _S7037;
                _S7068 = 0.0f;
            }
            else
            {
                _S7067 = 0.0f;
                _S7068 = _S7037;
            }
            if(_S7053)
            {
                float _S7069 = _S7067 + _S7068;
                float _S7070;
                if(_S7059)
                {
                    _S7070 = _S7067;
                }
                else
                {
                    DiffPair_float_0 _S7071;
                    (&_S7071)->primal_0 = _S7054;
                    (&_S7071)->differential_0 = 0.0f;
                    s_bwd_prop_exp_0(&_S7071, _S7067);
                    _S7070 = 0.90909093618392944f * _S7071.differential_0;
                }
                DiffPair_arrayx3Cfloatx2C8x3E_0 _S7072;
                (&_S7072)->primal_0 = _S7016;
                (&_S7072)->differential_0 = _S7043;
                DiffPair_vectorx3Cfloatx2C3x3E_0 _S7073;
                (&_S7073)->primal_0 = _S7057;
                (&_S7073)->differential_0 = _S7042;
                s_bwd_prop_interp_0(&_S7072, &_S7073, _S7070);
                float3  _S7074 = _S7073.differential_0 / _S7021;
                float3  _S7075 = _S7018.primal_0 * _S7074;
                float3  _S7076 = _S7058 * _S7074;
                float _S7077 = _S7075.x + _S7075.y + _S7075.z;
                DiffPair_float_0 _S7078;
                (&_S7078)->primal_0 = _S7030;
                (&_S7078)->differential_0 = 0.0f;
                DiffPair_float_0 _S7079;
                (&_S7079)->primal_0 = _S7035;
                (&_S7079)->differential_0 = 0.0f;
                DiffPair_float_0 _S7080;
                (&_S7080)->primal_0 = _S7055;
                (&_S7080)->differential_0 = 0.0f;
                s_bwd_prop_lerp_0(&_S7078, &_S7079, &_S7080, _S7077);
                float _S7081 = (&_S7072)->differential_0[int(0)] + _S7046[int(0)];
                float _S7082 = (&_S7072)->differential_0[int(1)] + _S7046[int(1)];
                float _S7083 = (&_S7072)->differential_0[int(2)] + _S7046[int(2)];
                float _S7084 = (&_S7072)->differential_0[int(3)] + _S7046[int(3)];
                float _S7085 = (&_S7072)->differential_0[int(4)] + _S7046[int(4)];
                float _S7086 = (&_S7072)->differential_0[int(5)] + _S7046[int(5)];
                float _S7087 = (&_S7072)->differential_0[int(6)] + _S7046[int(6)];
                float _S7088 = (&_S7072)->differential_0[int(7)] + _S7046[int(7)];
                float3  _S7089 = _S7074 + _S7045;
                float3  _S7090 = _S7076 + _S7044;
                float _S7091 = _S7079.differential_0 + _S7038;
                float _S7092 = _S7078.differential_0 + _S7039;
                _S7037 = _S7069;
                _S7038 = _S7091;
                _S7039 = _S7092;
                _S7044 = _S7090;
                _S7045 = _S7089;
                _S7046[int(0)] = _S7081;
                _S7046[int(1)] = _S7082;
                _S7046[int(2)] = _S7083;
                _S7046[int(3)] = _S7084;
                _S7046[int(4)] = _S7085;
                _S7046[int(5)] = _S7086;
                _S7046[int(6)] = _S7087;
                _S7046[int(7)] = _S7088;
            }
            else
            {
                _S7037 = _S7068;
            }
            _dc_0 = _dc_0 - int(1);
        }
        float _S7093 = - _S7052;
        float _S7094 = - _S7093 + _S7039;
        _S7037 = _S7093 + _S7038;
        _S7038 = _S7094;
    }
    else
    {
        _S7037 = 0.0f;
        _S7038 = 0.0f;
        _S7044 = _S7042;
        _S7045 = _S7042;
        _S7046[int(0)] = 0.0f;
        _S7046[int(1)] = 0.0f;
        _S7046[int(2)] = 0.0f;
        _S7046[int(3)] = 0.0f;
        _S7046[int(4)] = 0.0f;
        _S7046[int(5)] = 0.0f;
        _S7046[int(6)] = 0.0f;
        _S7046[int(7)] = 0.0f;
    }
    DiffPair_float_0 _S7095;
    (&_S7095)->primal_0 = _S7033;
    (&_S7095)->differential_0 = 0.0f;
    DiffPair_float_0 _S7096;
    (&_S7096)->primal_0 = _S7034;
    (&_S7096)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7095, &_S7096, _S7037);
    DiffPair_float_0 _S7097;
    (&_S7097)->primal_0 = _S7031;
    (&_S7097)->differential_0 = 0.0f;
    DiffPair_float_0 _S7098;
    (&_S7098)->primal_0 = _S7032;
    (&_S7098)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7097, &_S7098, _S7095.differential_0);
    DiffPair_float_0 _S7099;
    (&_S7099)->primal_0 = _S7027;
    (&_S7099)->differential_0 = 0.0f;
    DiffPair_float_0 _S7100;
    (&_S7100)->primal_0 = _S7029;
    (&_S7100)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7099, &_S7100, _S7038);
    DiffPair_float_0 _S7101;
    (&_S7101)->primal_0 = _S7028;
    (&_S7101)->differential_0 = 0.0f;
    DiffPair_float_0 _S7102;
    (&_S7102)->primal_0 = 0.0f;
    (&_S7102)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7101, &_S7102, _S7100.differential_0);
    DiffPair_float_0 _S7103;
    (&_S7103)->primal_0 = _S7025;
    (&_S7103)->differential_0 = 0.0f;
    DiffPair_float_0 _S7104;
    (&_S7104)->primal_0 = _S7026;
    (&_S7104)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7103, &_S7104, _S7099.differential_0);
    float3  s_diff_tb_T_0 = make_float3 (_S7097.differential_0, _S7098.differential_0, _S7096.differential_0);
    float3  s_diff_ta_T_0 = make_float3 (_S7103.differential_0, _S7104.differential_0, _S7101.differential_0);
    float3  s_diff_n_T_0 = - (s_diff_tb_T_0 + s_diff_ta_T_0);
    float3  _S7105 = _S7020 * (s_diff_tb_T_0 + - s_diff_ta_T_0);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7106;
    (&_S7106)->primal_0 = m_5;
    (&_S7106)->differential_0 = _S7042;
    s_bwd_prop_abs_1(&_S7106, _S7105);
    float3  _S7107 = m_5 * s_diff_n_T_0;
    float3  _S7108 = - ((_S7106.differential_0 + _S7023 * s_diff_n_T_0) / _S7022) + _S7044;
    dpray_d_5->primal_0 = (*dpray_d_5).primal_0;
    dpray_d_5->differential_0 = _S7108;
    float3  _S7109 = _S7107 + _S7045;
    dpray_o_5->primal_0 = (*dpray_o_5).primal_0;
    dpray_o_5->differential_0 = _S7109;
    dpdensities_3->primal_0 = dpdensities_3->primal_0;
    dpdensities_3->differential_0 = _S7046;
    return;
}

inline __device__ void s_bwd_evaluate_alpha_voxel_0(float3  _S7110, float _S7111, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S7112, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7113, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7114, float _S7115)
{
    FixedArray<float, 8>  _S7116 = _S7112->primal_0;
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S7117;
    float _S7118 = s_primal_ctx_evaluate_alpha_voxel_0(_S7110, _S7111, &_S7116, (*_S7113).primal_0, (*_S7114).primal_0, &_S7117);
    s_bwd_prop_evaluate_alpha_voxel_Intermediates_0 _S7119 = _S7117;
    s_bwd_prop_evaluate_alpha_voxel_0(_S7110, _S7111, _S7112, _S7113, _S7114, _S7115, &_S7119);
    return;
}

inline __device__ void evaluate_alpha_voxel_vjp(float3  pos_9, float size_9, FixedArray<float, 8>  * densities_8, float3  ray_o_12, float3  ray_d_12, float v_alpha_4, FixedArray<float, 8>  * v_densities_2, float3  * v_ray_o_5, float3  * v_ray_d_5)
{
    FixedArray<float, 8>  _S7120 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_0;
    (&dp_densities_0)->primal_0 = *densities_8;
    (&dp_densities_0)->differential_0 = _S7120;
    float3  _S7121 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_4;
    (&dp_ray_o_4)->primal_0 = ray_o_12;
    (&dp_ray_o_4)->differential_0 = _S7121;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_4;
    (&dp_ray_d_4)->primal_0 = ray_d_12;
    (&dp_ray_d_4)->differential_0 = _S7121;
    s_bwd_evaluate_alpha_voxel_0(pos_9, size_9, &dp_densities_0, &dp_ray_o_4, &dp_ray_d_4, v_alpha_4);
    *v_densities_2 = (&dp_densities_0)->differential_0;
    *v_ray_o_5 = dp_ray_o_4.differential_0;
    *v_ray_d_5 = dp_ray_d_4.differential_0;
    return;
}

inline __device__ void evaluate_color_voxel(float3  pos_10, float size_10, FixedArray<float, 8>  * densities_9, float3  rgb_16, float3  ray_o_13, float3  ray_d_13, float3  * out_rgb_1, float * depth_25)
{
    *out_rgb_1 = rgb_16;
    float _S7122 = 0.5f * size_10;
    float3  m_6 = make_float3 (1.0f) / ray_d_13;
    float3  k_20 = abs_0(m_6) * make_float3 (_S7122);
    float3  _S7123 = - (m_6 * (ray_o_13 - (pos_10 + make_float3 (_S7122))));
    float3  ta_4 = _S7123 - k_20;
    float3  tb_4 = _S7123 + k_20;
    float _S7124 = (F32_max(((F32_max((ta_4.x), (ta_4.y)))), ((F32_max((ta_4.z), (0.0f))))));
    float _S7125 = (F32_min(((F32_min((tb_4.x), (tb_4.y)))), (tb_4.z)));
    int i_14 = int(0);
    float accum_4 = 0.0f;
    float depth_accum_0 = 0.0f;
    for(;;)
    {
        if(i_14 < int(8))
        {
        }
        else
        {
            break;
        }
        float t_39 = lerp_0(_S7124, _S7125, (float(i_14) + 0.5f) / 8.0f);
        float _S7126 = interp_0(densities_9, (ray_o_13 + ray_d_13 * make_float3 (t_39) - pos_10) / make_float3 (size_10));
        float _S7127;
        if(_S7126 > 1.10000002384185791f)
        {
            _S7127 = _S7126;
        }
        else
        {
            _S7127 = (F32_exp((0.90909093618392944f * _S7126 - 0.90468984842300415f)));
        }
        float accum_5 = accum_4 + _S7127;
        float depth_accum_1 = depth_accum_0 + t_39 * _S7127;
        i_14 = i_14 + int(1);
        accum_4 = accum_5;
        depth_accum_0 = depth_accum_1;
    }
    *depth_25 = (F32_log(((F32_max((depth_accum_0 / accum_4), (0.0f))) + 9.99999997475242708e-07f)));
    return;
}

struct s_bwd_prop_evaluate_color_voxel_Intermediates_0
{
    float _S7128;
    float _S7129;
};

inline __device__ void s_primal_ctx_evaluate_color_voxel_0(float3  pos_11, float size_11, FixedArray<float, 8>  * dpdensities_4, float3  dprgb_1, float3  dpray_o_6, float3  dpray_d_6, float3  * dpout_rgb_1, float * dpdepth_3, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_16)
{
    _s_diff_ctx_16->_S7128 = 0.0f;
    _s_diff_ctx_16->_S7129 = 0.0f;
    float _S7130 = 0.5f * size_11;
    float3  m_7 = make_float3 (1.0f) / dpray_d_6;
    float3  k_21 = s_primal_ctx_abs_1(m_7) * make_float3 (_S7130);
    float3  _S7131 = - (m_7 * (dpray_o_6 - (pos_11 + make_float3 (_S7130))));
    float3  ta_5 = _S7131 - k_21;
    float3  tb_5 = _S7131 + k_21;
    float _S7132 = s_primal_ctx_max_0(s_primal_ctx_max_0(ta_5.x, ta_5.y), s_primal_ctx_max_0(ta_5.z, 0.0f));
    float _S7133 = s_primal_ctx_min_0(s_primal_ctx_min_0(tb_5.x, tb_5.y), tb_5.z);
    bool _runFlag_1 = true;
    float accum_6 = 0.0f;
    float depth_accum_2 = 0.0f;
    int i_15 = int(0);
    int _pc_1 = int(0);
    for(;;)
    {
        _s_diff_ctx_16->_S7128 = depth_accum_2;
        _s_diff_ctx_16->_S7129 = accum_6;
        if(_runFlag_1)
        {
        }
        else
        {
            break;
        }
        float _S7134;
        float _S7135;
        int _S7136;
        if(i_15 < int(8))
        {
            float _S7137 = s_primal_ctx_lerp_0(_S7132, _S7133, (float(i_15) + 0.5f) / 8.0f);
            float _S7138 = s_primal_ctx_interp_0(dpdensities_4, (dpray_o_6 + dpray_d_6 * make_float3 (_S7137) - pos_11) / make_float3 (size_11));
            if(_S7138 > 1.10000002384185791f)
            {
                _S7134 = _S7138;
            }
            else
            {
                _S7134 = s_primal_ctx_exp_1(0.90909093618392944f * _S7138 - 0.90468984842300415f);
            }
            float accum_7 = accum_6 + _S7134;
            float depth_accum_3 = depth_accum_2 + _S7137 * _S7134;
            _S7136 = int(1);
            _S7134 = accum_7;
            _S7135 = depth_accum_3;
        }
        else
        {
            _S7136 = int(0);
            _S7134 = 0.0f;
            _S7135 = 0.0f;
        }
        if(_S7136 != int(1))
        {
            _runFlag_1 = false;
        }
        if(_runFlag_1)
        {
            int _S7139 = i_15 + int(1);
            accum_6 = _S7134;
            depth_accum_2 = _S7135;
            i_15 = _S7139;
        }
        _pc_1 = _pc_1 + int(1);
    }
    float _S7140 = s_primal_ctx_log_0(s_primal_ctx_max_0(depth_accum_2 / accum_6, 0.0f) + 9.99999997475242708e-07f);
    *dpout_rgb_1 = dprgb_1;
    *dpdepth_3 = _S7140;
    return;
}

inline __device__ void s_bwd_prop_evaluate_color_voxel_0(float3  pos_12, float size_12, DiffPair_arrayx3Cfloatx2C8x3E_0 * dpdensities_5, DiffPair_vectorx3Cfloatx2C3x3E_0 * dprgb_2, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_o_7, DiffPair_vectorx3Cfloatx2C3x3E_0 * dpray_d_7, float3  dpout_rgb_2, float dpdepth_4, s_bwd_prop_evaluate_color_voxel_Intermediates_0 * _s_diff_ctx_17)
{
    FixedArray<float, 8>  _S7141 = dpdensities_5->primal_0;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7142 = *dpray_o_7;
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7143 = *dpray_d_7;
    float3  _S7144 = make_float3 (size_12);
    float _S7145 = 0.5f * size_12;
    float3  _S7146 = make_float3 (_S7145);
    float3  m_8 = make_float3 (1.0f) / (*dpray_d_7).primal_0;
    float3  _S7147 = (*dpray_d_7).primal_0 * (*dpray_d_7).primal_0;
    float3  _S7148 = (*dpray_o_7).primal_0 - (pos_12 + make_float3 (_S7145));
    float3  k_22 = s_primal_ctx_abs_1(m_8) * make_float3 (_S7145);
    float3  _S7149 = - (m_8 * _S7148);
    float3  ta_6 = _S7149 - k_22;
    float3  tb_6 = _S7149 + k_22;
    float _S7150 = ta_6.x;
    float _S7151 = ta_6.y;
    float _S7152 = s_primal_ctx_max_0(_S7150, _S7151);
    float _S7153 = ta_6.z;
    float _S7154 = s_primal_ctx_max_0(_S7153, 0.0f);
    float _S7155 = s_primal_ctx_max_0(_S7152, _S7154);
    float _S7156 = tb_6.x;
    float _S7157 = tb_6.y;
    float _S7158 = s_primal_ctx_min_0(_S7156, _S7157);
    float _S7159 = tb_6.z;
    float _S7160 = s_primal_ctx_min_0(_S7158, _S7159);
    float _S7161 = _s_diff_ctx_17->_S7128 / _s_diff_ctx_17->_S7129;
    float _S7162 = _s_diff_ctx_17->_S7129 * _s_diff_ctx_17->_S7129;
    float3  _S7163 = make_float3 (0.0f);
    FixedArray<float, 8>  _S7164 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_float_0 _S7165;
    (&_S7165)->primal_0 = s_primal_ctx_max_0(_S7161, 0.0f) + 9.99999997475242708e-07f;
    (&_S7165)->differential_0 = 0.0f;
    s_bwd_prop_log_0(&_S7165, dpdepth_4);
    DiffPair_float_0 _S7166;
    (&_S7166)->primal_0 = _S7161;
    (&_S7166)->differential_0 = 0.0f;
    DiffPair_float_0 _S7167;
    (&_S7167)->primal_0 = 0.0f;
    (&_S7167)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7166, &_S7167, _S7165.differential_0);
    float _S7168 = _S7166.differential_0 / _S7162;
    float _S7169 = _s_diff_ctx_17->_S7128 * - _S7168;
    float _S7170 = _s_diff_ctx_17->_S7129 * _S7168;
    int _dc_1 = int(8);
    float _S7171 = _S7169;
    float _S7172 = _S7170;
    float _S7173 = 0.0f;
    float _S7174 = 0.0f;
    float3  _S7175 = _S7163;
    float3  _S7176 = _S7163;
    FixedArray<float, 8>  _S7177;
    _S7177[int(0)] = 0.0f;
    _S7177[int(1)] = 0.0f;
    _S7177[int(2)] = 0.0f;
    _S7177[int(3)] = 0.0f;
    _S7177[int(4)] = 0.0f;
    _S7177[int(5)] = 0.0f;
    _S7177[int(6)] = 0.0f;
    _S7177[int(7)] = 0.0f;
    for(;;)
    {
        if(_dc_1 >= int(0))
        {
        }
        else
        {
            break;
        }
        bool _S7178 = _dc_1 < int(8);
        int _S7179;
        float _S7180;
        float _S7181;
        float _S7182;
        float _S7183;
        float3  _S7184;
        float3  _S7185;
        bool _S7186;
        if(_S7178)
        {
            float _S7187 = (float(_dc_1) + 0.5f) / 8.0f;
            float _S7188 = s_primal_ctx_lerp_0(_S7155, _S7160, _S7187);
            float3  _S7189 = make_float3 (_S7188);
            float3  _S7190 = (_S7142.primal_0 + _S7143.primal_0 * make_float3 (_S7188) - pos_12) / make_float3 (size_12);
            FixedArray<float, 8>  _S7191 = _S7141;
            float _S7192 = s_primal_ctx_interp_0(&_S7191, _S7190);
            bool _S7193 = _S7192 > 1.10000002384185791f;
            if(_S7193)
            {
                _S7180 = _S7192;
                _S7181 = 0.0f;
            }
            else
            {
                float _S7194 = 0.90909093618392944f * _S7192 - 0.90468984842300415f;
                _S7180 = s_primal_ctx_exp_1(_S7194);
                _S7181 = _S7194;
            }
            float _S7195 = _S7180;
            float _S7196 = _S7181;
            _S7179 = int(1);
            _S7180 = _S7188;
            _S7181 = _S7195;
            _S7186 = _S7193;
            _S7182 = _S7196;
            _S7184 = _S7190;
            _S7185 = _S7189;
            _S7183 = _S7187;
        }
        else
        {
            _S7179 = int(0);
            _S7180 = 0.0f;
            _S7181 = 0.0f;
            _S7186 = false;
            _S7182 = 0.0f;
            _S7184 = _S7163;
            _S7185 = _S7163;
            _S7183 = 0.0f;
        }
        float _S7197;
        float _S7198;
        float _S7199;
        float _S7200;
        if(!(_S7179 != int(1)))
        {
            _S7197 = _S7171;
            _S7198 = _S7172;
            _S7199 = 0.0f;
            _S7200 = 0.0f;
        }
        else
        {
            _S7197 = 0.0f;
            _S7198 = 0.0f;
            _S7199 = _S7172;
            _S7200 = _S7171;
        }
        if(_S7178)
        {
            float _S7201 = _S7181 * _S7198;
            float _S7202 = _S7198 + _S7199;
            float _S7203 = _S7180 * _S7198 + _S7197;
            float _S7204 = _S7197 + _S7200;
            float _S7205;
            if(_S7186)
            {
                _S7205 = _S7203;
            }
            else
            {
                DiffPair_float_0 _S7206;
                (&_S7206)->primal_0 = _S7182;
                (&_S7206)->differential_0 = 0.0f;
                s_bwd_prop_exp_0(&_S7206, _S7203);
                _S7205 = 0.90909093618392944f * _S7206.differential_0;
            }
            DiffPair_arrayx3Cfloatx2C8x3E_0 _S7207;
            (&_S7207)->primal_0 = _S7141;
            (&_S7207)->differential_0 = _S7164;
            DiffPair_vectorx3Cfloatx2C3x3E_0 _S7208;
            (&_S7208)->primal_0 = _S7184;
            (&_S7208)->differential_0 = _S7163;
            s_bwd_prop_interp_0(&_S7207, &_S7208, _S7205);
            float3  _S7209 = _S7208.differential_0 / _S7144;
            float3  _S7210 = _S7143.primal_0 * _S7209;
            float3  _S7211 = _S7185 * _S7209;
            float _S7212 = _S7210.x + _S7210.y + _S7210.z + _S7201;
            DiffPair_float_0 _S7213;
            (&_S7213)->primal_0 = _S7155;
            (&_S7213)->differential_0 = 0.0f;
            DiffPair_float_0 _S7214;
            (&_S7214)->primal_0 = _S7160;
            (&_S7214)->differential_0 = 0.0f;
            DiffPair_float_0 _S7215;
            (&_S7215)->primal_0 = _S7183;
            (&_S7215)->differential_0 = 0.0f;
            s_bwd_prop_lerp_0(&_S7213, &_S7214, &_S7215, _S7212);
            float _S7216 = (&_S7207)->differential_0[int(0)] + _S7177[int(0)];
            float _S7217 = (&_S7207)->differential_0[int(1)] + _S7177[int(1)];
            float _S7218 = (&_S7207)->differential_0[int(2)] + _S7177[int(2)];
            float _S7219 = (&_S7207)->differential_0[int(3)] + _S7177[int(3)];
            float _S7220 = (&_S7207)->differential_0[int(4)] + _S7177[int(4)];
            float _S7221 = (&_S7207)->differential_0[int(5)] + _S7177[int(5)];
            float _S7222 = (&_S7207)->differential_0[int(6)] + _S7177[int(6)];
            float _S7223 = (&_S7207)->differential_0[int(7)] + _S7177[int(7)];
            float3  _S7224 = _S7209 + _S7176;
            float3  _S7225 = _S7211 + _S7175;
            float _S7226 = _S7214.differential_0 + _S7173;
            float _S7227 = _S7213.differential_0 + _S7174;
            _S7171 = _S7204;
            _S7172 = _S7202;
            _S7173 = _S7226;
            _S7174 = _S7227;
            _S7175 = _S7225;
            _S7176 = _S7224;
            _S7177[int(0)] = _S7216;
            _S7177[int(1)] = _S7217;
            _S7177[int(2)] = _S7218;
            _S7177[int(3)] = _S7219;
            _S7177[int(4)] = _S7220;
            _S7177[int(5)] = _S7221;
            _S7177[int(6)] = _S7222;
            _S7177[int(7)] = _S7223;
        }
        else
        {
            _S7171 = _S7200;
            _S7172 = _S7199;
        }
        _dc_1 = _dc_1 - int(1);
    }
    DiffPair_float_0 _S7228;
    (&_S7228)->primal_0 = _S7158;
    (&_S7228)->differential_0 = 0.0f;
    DiffPair_float_0 _S7229;
    (&_S7229)->primal_0 = _S7159;
    (&_S7229)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7228, &_S7229, _S7173);
    DiffPair_float_0 _S7230;
    (&_S7230)->primal_0 = _S7156;
    (&_S7230)->differential_0 = 0.0f;
    DiffPair_float_0 _S7231;
    (&_S7231)->primal_0 = _S7157;
    (&_S7231)->differential_0 = 0.0f;
    s_bwd_prop_min_0(&_S7230, &_S7231, _S7228.differential_0);
    DiffPair_float_0 _S7232;
    (&_S7232)->primal_0 = _S7152;
    (&_S7232)->differential_0 = 0.0f;
    DiffPair_float_0 _S7233;
    (&_S7233)->primal_0 = _S7154;
    (&_S7233)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7232, &_S7233, _S7174);
    DiffPair_float_0 _S7234;
    (&_S7234)->primal_0 = _S7153;
    (&_S7234)->differential_0 = 0.0f;
    DiffPair_float_0 _S7235;
    (&_S7235)->primal_0 = 0.0f;
    (&_S7235)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7234, &_S7235, _S7233.differential_0);
    DiffPair_float_0 _S7236;
    (&_S7236)->primal_0 = _S7150;
    (&_S7236)->differential_0 = 0.0f;
    DiffPair_float_0 _S7237;
    (&_S7237)->primal_0 = _S7151;
    (&_S7237)->differential_0 = 0.0f;
    s_bwd_prop_max_1(&_S7236, &_S7237, _S7232.differential_0);
    float3  s_diff_tb_T_1 = make_float3 (_S7230.differential_0, _S7231.differential_0, _S7229.differential_0);
    float3  s_diff_ta_T_1 = make_float3 (_S7236.differential_0, _S7237.differential_0, _S7234.differential_0);
    float3  s_diff_n_T_1 = - (s_diff_tb_T_1 + s_diff_ta_T_1);
    float3  _S7238 = _S7146 * (s_diff_tb_T_1 + - s_diff_ta_T_1);
    DiffPair_vectorx3Cfloatx2C3x3E_0 _S7239;
    (&_S7239)->primal_0 = m_8;
    (&_S7239)->differential_0 = _S7163;
    s_bwd_prop_abs_1(&_S7239, _S7238);
    float3  _S7240 = m_8 * s_diff_n_T_1;
    float3  _S7241 = - ((_S7239.differential_0 + _S7148 * s_diff_n_T_1) / _S7147) + _S7175;
    dpray_d_7->primal_0 = (*dpray_d_7).primal_0;
    dpray_d_7->differential_0 = _S7241;
    float3  _S7242 = _S7240 + _S7176;
    dpray_o_7->primal_0 = (*dpray_o_7).primal_0;
    dpray_o_7->differential_0 = _S7242;
    dprgb_2->primal_0 = (*dprgb_2).primal_0;
    dprgb_2->differential_0 = dpout_rgb_2;
    dpdensities_5->primal_0 = dpdensities_5->primal_0;
    dpdensities_5->differential_0 = _S7177;
    return;
}

inline __device__ void s_bwd_evaluate_color_voxel_0(float3  _S7243, float _S7244, DiffPair_arrayx3Cfloatx2C8x3E_0 * _S7245, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7246, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7247, DiffPair_vectorx3Cfloatx2C3x3E_0 * _S7248, float3  _S7249, float _S7250)
{
    FixedArray<float, 8>  _S7251 = _S7245->primal_0;
    float3  _S7252;
    float _S7253;
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S7254;
    s_primal_ctx_evaluate_color_voxel_0(_S7243, _S7244, &_S7251, (*_S7246).primal_0, (*_S7247).primal_0, (*_S7248).primal_0, &_S7252, &_S7253, &_S7254);
    s_bwd_prop_evaluate_color_voxel_Intermediates_0 _S7255 = _S7254;
    s_bwd_prop_evaluate_color_voxel_0(_S7243, _S7244, _S7245, _S7246, _S7247, _S7248, _S7249, _S7250, &_S7255);
    return;
}

inline __device__ void evaluate_color_voxel_vjp(float3  pos_13, float size_13, FixedArray<float, 8>  * densities_10, float3  rgb_17, float3  ray_o_14, float3  ray_d_14, float3  v_out_rgb_1, float v_depth_12, FixedArray<float, 8>  * v_densities_3, float3  * v_rgb_10, float3  * v_ray_o_6, float3  * v_ray_d_6)
{
    FixedArray<float, 8>  _S7256 = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
    DiffPair_arrayx3Cfloatx2C8x3E_0 dp_densities_1;
    (&dp_densities_1)->primal_0 = *densities_10;
    (&dp_densities_1)->differential_0 = _S7256;
    float3  _S7257 = make_float3 (0.0f);
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_rgb_1;
    (&dp_rgb_1)->primal_0 = rgb_17;
    (&dp_rgb_1)->differential_0 = _S7257;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_o_5;
    (&dp_ray_o_5)->primal_0 = ray_o_14;
    (&dp_ray_o_5)->differential_0 = _S7257;
    DiffPair_vectorx3Cfloatx2C3x3E_0 dp_ray_d_5;
    (&dp_ray_d_5)->primal_0 = ray_d_14;
    (&dp_ray_d_5)->differential_0 = _S7257;
    s_bwd_evaluate_color_voxel_0(pos_13, size_13, &dp_densities_1, &dp_rgb_1, &dp_ray_o_5, &dp_ray_d_5, v_out_rgb_1, v_depth_12);
    *v_densities_3 = (&dp_densities_1)->differential_0;
    *v_rgb_10 = dp_rgb_1.differential_0;
    *v_ray_o_6 = dp_ray_o_5.differential_0;
    *v_ray_d_6 = dp_ray_d_5.differential_0;
    return;
}

